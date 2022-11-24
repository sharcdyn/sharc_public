program spectra
implicit none
integer nat,npot,ntraj
real*8, allocatable :: geom(:,:), vel(:,:),mass(:)
character*2, allocatable :: sat(:)
real*8 Ekin,Eref,Ehar, dip(3)
real*8, allocatable :: Epot(:),prob(:)
complex*16, allocatable :: dipole(:,:,:),R(:,:)
integer pot,n,nseed,ntrajf
integer, allocatable :: seed(:)
integer printlevel

logical velcorrection,acceptvel,acceptsp
integer nomodel

real*8 Ei,Ef,rando(1),ren
real*8 E0,Et,gwidth
integer nEn
real*8, allocatable :: spec(:,:,:)
real*8 maxprob,maxprobt

character*2 sat2, sat3
real*8 kkat
integer, allocatable ::  noat(:)
integer noat2,noat3
real*8 velcorr


logical ex,res
integer i,j,k

character*1 kkchar
real*8 kk
real*8, allocatable :: kk2(:)


write(6,*) "Number of atoms and number of potentials"
read(5,*) nat,npot
write(6,*) "Potential where the initial conditions where created (usually 1)",&
" and printlevel (0 minimum, 1 (+energy.out), 2 (+individual.out) or 3 (+dipole*.out)"
read(5,*) pot,printlevel
write(6,*) "Last trajectory employed"
read(5,*) ntrajf

allocate (geom(nat,3),vel(nat,3),sat(nat),noat(nat),mass(nat) )
allocate (Epot(npot),dipole(3,npot,npot),prob(npot) )
allocate (R(npot,npot))

allocate (kk2(npot) )

write(6,*) "Do you want kinetic energy correction (y/n)?"
read(5,*) kkchar
velcorrection=.true.
if (kkchar.eq."n") then
 velcorrection=.false.
endif

write(6,*) "Range of energies and number of points for the spectrum (eV)"
read(5,*) E0,Et,nEn
E0=E0/27.211
Et=Et/27.211
write(6,*) E0,Et,nEn

write(6,*) "Dipoles to scan X,Y,Z, (scale factor)"
read(5,*) dip(1),dip(2),dip(3)
write(6,*) "X  ", dip(1),"Y  ", dip(2),"Z  ",dip(3)

allocate (spec(3,npot,nEn) )

write(6,*) "Gaussian width of the individual trajectories (eV)"
read(5,*) gwidth
gwidth=gwidth/27.211
write(6,*) gwidth

write(6,*) "Range where the trajectories are gonna be allowed (eV) and zero if you do not want to restrict (1 otherwise)"
read(5,*) Ei,Ef,i
Ei=Ei/27.211
Ef=Ef/27.211
write(6,*) Ei,Ef,i

if (i.eq.0) then
 nseed=1
 allocate(seed(nseed))
 seed(1)=0
else
 inquire(file="seeds.dat",exist=ex)
 if (ex .eqv. .false.) then
  open(1,file="seeds.dat")
  write(6,*) "Generating seeds and storage them in seeds.dat"
  call random_seed(size=nseed)
  write(1,*) nseed
  allocate (seed(nseed))
  inquire(file="/dev/urandom",exist=res)
  if (res) then
   open(2, file='/dev/urandom', access='stream', form='UNFORMATTED')
   read(2) seed
   close(2)
  else
   call system_clock(count=j)
   seed = j + 37 * (/ (i - 1, i = 1, nseed) /)
  endif
  write(1,*) (seed(i),i=1,nseed)
  close(1)
  deallocate(seed)
 endif
 open(1,file="seeds.dat")
 read(1,*) nseed
 allocate(seed(nseed))
 read(1,*) (seed(i),i=1,nseed)
 close(1)
endif


ntraj=0
inquire(file="geom",exist=ex)
if (ex) then
 open(1,file="geom")
 do i=1,nat
  read(1,*) sat(i),kkat,(geom(i,j),j=1,3),mass(i)
  noat(i)=int(kkat+0.5)
  do j=1,3
   vel(i,j)=0.d0
  enddo
 enddo
 close(1)

else
 write(6,*) "File geom with the equilibrium geometry must exist"
 stop
endif

open(11,file="geoms.init")
open(13,file="geoms")
open(15,file="selection.xyz")
open(16,file="problems.out")
open(17,file="selection.diab")
if (printlevel.ge.1) then
 open(20,file="energy.out")
endif
if (printlevel.ge.2) then
 open(18,file="individual.out")
endif
if (printlevel.ge.3) then
 open(21,file="dipole1.out")
 open(22,file="dipole2.out")
 open(23,file="dipole3.out")
endif
write(16,*) ""

inquire(file="model",exist=ex)
nomodel=0
if (ex) then
 write(6,*) "Using model potentials"
 nomodel=1
 open(1,file="model")
 kkchar=""
 read(1,*,iostat=i) kkchar
 if (i.eq.0) nomodel=2
endif

call getdata(npot,Epot,dipole,0,nat,geom,sat,R,nomodel,0.d0)

Eref=Epot(pot)
write(6,*) "Reference energy ",pot,Eref

!!! Writting dipole of equilibrium geometry
open(1,file="eq_dipole1.out")
open(2,file="eq_dipole2.out")
open(3,file="eq_dipole3.out")
do i=1,3
 do j=1,npot
  do k=1,npot
   write(i,"(2(x,I6.6),5(x,F10.5))") j,k,(Epot(j)-Eref)*27.211,(Epot(k)-Eref)*27.211,&
                                       cdabs(dipole(i,j,k)),dipole(i,j,k)
  enddo
 enddo
 close(i)
enddo

do i=1,3
 do j=1,npot
  do k=1,npot
   dipole(i,j,k)=dipole(i,j,k)*dip(i)
  enddo
 enddo
enddo

! First point of the spectrum created with the equilibrium geom
do i=1,npot
 do j=1,nEn
  do k=1,3
   spec(k,i,j)=0.d0
  enddo
 enddo
enddo
i=0
if (printlevel.ge.2) i=18
call gauss(npot,pot,dipole,Epot,nEn,E0,Et,spec,gwidth,prob,i)
maxprob=0.
maxprobt=0.
write(6,*) "Probability to be excited at equilibrium"
do i=1,npot
 write(6,"(A,I3,F7.2,A,5(x,E20.10e3))") "Pot ",i,(Epot(i)-Eref)*27.211," eV ",&
  Epot(i)-Eref,prob(i),&
  (2./3.*(Epot(i)-Epot(pot))*cdabs(dipole(j,pot,i))**2,j=1,3)
 if (prob(i).ge.maxprobt) maxprobt=prob(i)
 if (Epot(i)-Eref.ge.Ei .and. Epot(i)-Eref.le.Ef .and. prob(i).ge.maxprob) maxprob=prob(i)
enddo
write(6,*) "Write renormalization factor"
read(5,*) ren
write(6,*) "Renormalized probability ",ren
do i=1,npot
 write(6,*) "Pot ",i,prob(i)/ren
enddo
!if (seed(1).ne.0) call srand(seed(1))
if (seed(1).ne.0) call random_seed(put=seed)



read(11,*) n
if (n.ne.nat) then
 write(6,*) "The specified number of atoms and the read from geoms.init are not the same"
 write(6,*) nat,n
 stop
endif

do i=1,nat
 read(11,*) sat2,kkat
 noat2=int(kkat+0.5)
 if (noat2.ne.noat(i))  then
  write(6,*) "Atomic numbers of geom, geoms.init and veloc.init are not coincident"
  write(6,*) "geom ",i,noat(i)
  write(6,*) "geoms.init ",i,noat2
  stop
 endif
enddo

n=0
Ekin=0.
if (printlevel.ge.1) write(20,905) n,Ekin,(Epot(i)-Eref,i=1,npot)
if (printlevel.ge.3) then
 do i=1,3
  write(20+i,905) n,((dipole(i,j,k),k=1,npot),j=1,npot)
 enddo
endif


 read(11,*,END=10) i,kkchar,Ehar,Ekin
 if (i.eq.0) then
  do i=1,nat
   read(11,*) (geom(i,j),j=1,3)
  enddo
  do i=1,nat
   read(11,*) (vel(i,j),j=1,3)
  enddo
 else
  rewind(11)
 endif

n=0
do while(n.lt.ntrajf)
 n=n+1
 read(11,*,END=10) i,kkchar,Ehar,Ekin
 if (i.ne.n) then
  write(16,*) "Trajectory number unexpected"
  write(16,*) "Expected ",n
  write(16,*) "Reading ",i
  stop "Trajectory number unexpected"
 endif
 do i=1,nat
  read(11,*) (geom(i,j),j=1,3)
 enddo
 do i=1,nat
  read(11,*) (vel(i,j),j=1,3)
 enddo

 call getdata(npot,Epot,dipole,n,nat,geom,sat,R,nomodel,Ehar)

 do i=1,3
  do j=1,npot
   do k=1,npot
    dipole(i,j,k)=dipole(i,j,k)*dip(i)
   enddo
  enddo
 enddo


 do i=1,npot
  Epot(i)=Epot(i)-Eref
 enddo


 kk=Ekin+Ehar-Epot(pot)
 if (kk.ge.0) then
  velcorr=dsqrt(kk/Ekin)
 else
  velcorr=-dsqrt(-kk/Ekin)
 endif
 write(6,"(4(A,F12.8))") "Eharm ",Ehar," Eharm+Ekin",Ekin+Ehar," Epot ",Epot(pot)," Vel corr ",velcorr
 acceptvel=.true.
 if (velcorrection) then
  if (velcorr.lt.0) then
   write(6,*) "It is not possible to adjust velocity in traj (skipping)",n
   write(16,*) "It is not possible to adjust velocity in traj (skipping)",n
   acceptvel=.false.
  else
   write(6,*) "Adjusting velocity traj ",n,velcorr
   Ekin=0.d0
   do i=1,nat
    do j=1,3
     vel(i,j)=velcorr*vel(i,j)
     Ekin=Ekin+0.5*mass(i)*1822.888*vel(i,j)*vel(i,j)
    enddo
   enddo
  endif
 endif

 if (printlevel.ge.1) write(20,905) n,Ekin,(Epot(i),i=1,npot)
 if (printlevel.ge.3) then
  do i=1,3
   write(20+i,905) n,((dipole(i,j,k),k=1,npot),j=1,npot)
  enddo
 endif

 if (acceptvel) then
  i=0
  if (printlevel.ge.2) i=18
  call gauss(npot,pot,dipole,Epot,nEn,E0,Et,spec,gwidth,prob,i)
 endif

 do k=1,npot
  acceptsp=.true.
  ex=.true.
  ! Checking the excitation interval
  if (prob(k).ge.maxprobt) maxprobt=prob(k)
  if (Epot(k)-Epot(pot).le.Ei .or. Epot(k)-Epot(pot).gt.Ef) then
   ex=.false.
  else
   if (prob(k).ge.maxprob) maxprob=prob(k)
  endif

  if (seed(1).ne.0.and.ex) then
   if (prob(k)/ren.gt.1) then
    write(6,"(A,I10,E20.10e3)") "Warning probability more than 1... selection of trajectories is meanless",n,prob(k)
    write(16,"(A,I10,E20.10e3)") "Warning probability more than 1... selection of trajectories is meanless",n,prob(k)
   endif
  !  rando(1)=rand()
   call random_number(rando)
   if (rando(1).gt.prob(k)/ren) acceptsp=.false.
  endif

  if (acceptsp.and.acceptvel.and.ex) then
   ntraj=ntraj+1
   if (ntraj.eq.1) then
    write(13,*) nat
    do i=1,nat
     write(13,"(A2,x,I3,x,E20.10e3)") sat(i),noat(i),mass(i)
    enddo
    write(13,*) "newtraj Traj Old oldtraj, Initial pot,Exciation pot, Initial pot energy, Final pot energy, probability"
   endif
   write(13,901) ntraj," Traj Old ",n,pot,k,Epot(pot),Epot(k),Ekin,prob(k),velcorr
   write(15,*) nat
   write(15,902) k,Epot(pot)+Ekin,Ekin,Epot(pot),Epot(k), prob(k),velcorr
   do i=1,nat
    write(13,903) (geom(i,j),j=1,3)
    write(15,904) sat(i),(geom(i,j)*.5292,j=1,3),(vel(i,j)*1000.,j=1,3)
   enddo
   do i=1,nat
    write(13,903) (vel(i,j),j=1,3)
   enddo
   write(17,"(A,I0,A,2(I0,x),E20.10e3,A,E20.10e3,A,(10000(x,I5.5,x,E30.20e3,x,E30.20e3)))") &
              "Traj ",ntraj," Old ",n,k,Epot(k)-Epot(pot)," prob ",prob(k),&
              " diabatic ",(i,R(i,k),i=1,npot)
  endif

  write(6,*) "Done traj ",n,ntrajf,ntraj,k,acceptvel,ex,acceptsp

 enddo
enddo

10 continue

write(6,"(2(A,F10.6))") "Max probability found in the energy interval ",maxprob, " and total ",maxprobt

open(20,file="spectrum.dat")
open(21,file="spectrum1.dat")
open(22,file="spectrum2.dat")
open(23,file="spectrum3.dat")
do i=1,nEn
 Eref=E0+(i-1)*(Et-E0)/float(nEn-1)
 do k=1,3
  kk=0.d0
  do j=1,npot
   spec(k,j,i)=spec(k,j,i)/float(n)
   kk=kk+spec(k,j,i)
  enddo
  write(20+k,900) Eref,kk,(spec(k,j,i),j=1,npot)
 enddo
 kk=0.d0
 do j=1,npot
  kk2(j)=0.d0
  do k=1,3
   kk2(j)=kk2(j)+spec(k,j,i)
   kk=kk+spec(k,j,i)
  enddo
 enddo
 write(20,900) Eref,kk,(kk2(j),j=1,npot)
enddo

900 format(1000(x,E20.10e3))
901 format(I0,A,I0,2(x,I4.4),5(x,E20.10e3) )
902 format(I0,20(x,E20.10e3))
903 format(3(x,E20.10e3))
904 format(A2,6(x,E20.10e3))
905 format(I10.10,100000(x,E20.10e3))
end


