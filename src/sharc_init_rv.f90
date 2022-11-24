program init_rv
implicit none
integer nat
real*8, allocatable :: geom(:,:),vel(:,:),noat(:),mass(:)
character*2, allocatable :: sat(:)
real*8 Esamp
integer traj0,trajf
integer nseed
integer, allocatable :: seed(:)
logical res,ex
logical, allocatable :: atmask(:)

real*8, allocatable :: prob(:)
real*8 sumprob
real*8 Ek,mom(3),cm(3),masscm


integer i,j,k
integer itraj

write(6,*) "Number of atoms, and kinetic energy to be sampled (in eV), if negative Temperature is sampled (in K)"
read(5,*) nat,Esamp
if (Esamp.lt.0) then
 write(6,*) "Sampling temperature ",-Esamp
 !!! Ekin=3/2 N kB T
 Esamp=-Esamp*3.d0/2.d0          !!! 3/2 T
 Esamp=Esamp*8.617333262145d-5   !!! kB (eV/K/º)
 Esamp=Esamp*float(nat)
endif
Esamp=Esamp/27.211
write(6,*) "Sampling kinetic energy (H,eV)",Esamp,Esamp*27.211
open(1,file="geom",status="old",iostat=i)
if (i.ne.0) stop "File geom does not exist"
allocate(sat(nat),noat(nat),mass(nat),geom(nat,3),vel(nat,3))
do i=1,nat
 read(1,*) sat(i),noat(i),(geom(i,j),j=1,3),mass(i)
enddo
close(1)
write(6,*) "Equilibrium geometry ready"

write(6,*) "First trajectory and last"
read(5,*) traj0,trajf

allocate(atmask(nat))
atmask=.true.
write(6,*) "Number of atoms to consider (0 all, if negative the first ones are selected)"
read(5,*) i
if (Esamp.lt.0) then
 Esamp=Esamp/float(nat)*float(i)
endif
if (i.ne.0) atmask=.false.
if (i.lt.0) then
 do j=1,abs(i)
  atmask(j)=.true.
 enddo
endif
if (i.gt.0) then
 allocate(seed(i))
 write(6,*) "List of atoms to select"
 read(5,*) seed(:)
 do j=1,i
  atmask(seed(j))=.true.
 enddo
 deallocate(seed)
endif

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
call random_seed(put=seed)
write(6,*) "Number of seeds ",nseed
write(6,*) "Seeds ",(seed(i),i=1,nseed)

allocate(prob(3*nat))
open(1,file="geoms.xyz")
open(2,file="geoms.init")
vel=0.
sumprob=0.
i=0
write(1,*) nat
write(1,903) sumprob,sumprob,sumprob,Esamp," Traj ",i
do i=1,nat
 write(1,901) sat(i),(geom(i,j)*.5292,j=1,3),(vel(i,j),j=1,3)
enddo
write(2,*) nat
do i=1,nat
 write(2,901) sat(i),noat(i),mass(i)
 mass(i)=mass(i)*1822.888
enddo

!!! Center of masses
masscm=0.
do i=1,nat
 masscm=masscm+mass(i)
enddo
do j=1,3
 cm(j)=0.
 do i=1,nat
  cm(j)=cm(j)+mass(i)*geom(i,j)
 enddo
 cm(j)=cm(j)/masscm
enddo

 

call random_seed(put=seed)
do itraj=traj0,trajf

!! Creating random velocities int he selected atoms
 call random_number(prob)
 sumprob=0.
 do i=1,3*nat
  prob(i)=prob(i)-.5
  if (atmask(i)) sumprob=sumprob+abs(prob(i))
 enddo
 do i=1,3*nat
  prob(i)=prob(i)/sumprob
 enddo
 sumprob=0.
 do i=1,3*nat
  if (atmask(i)) sumprob=sumprob+abs(prob(i))
 enddo
 k=0
 vel=0.d0
 do i=1,nat
  do j=1,3
   k=k+1
   if (atmask(i))  then
    vel(i,j)=abs(prob(k))  !!! Energy for the coordinate with no units
    vel(i,j)=sqrt(2.*vel(i,j)/mass(i))
    if (prob(k).lt.0) vel(i,j)=-vel(i,j)
   endif
  enddo
 enddo

!!! Removing velocity of the center of masses using the total momentum
 do j=1,3
  mom(j)=0.
  do i=1,nat
   mom(j)=mom(j)+mass(i)*vel(i,j)
  enddo
 enddo
 do i=1,nat
  do j=1,3
   vel(i,j)=vel(i,j)-mom(j)/masscm
  enddo
 enddo

!!! Kinetic energy
 Ek=0.
 do i=1,nat
  do j=1,3
   Ek=Ek+mass(i)*vel(i,j)*vel(i,j)*.5
  enddo
 enddo
!!! Rescaling the velocity
 do i=1,nat
  do j=1,3
   vel(i,j)=sqrt(Esamp/Ek)*vel(i,j)
  enddo
 enddo
!!! Kinetic energy
 Ek=0.
 do i=1,nat
  do j=1,3
   Ek=Ek+mass(i)*vel(i,j)*vel(i,j)*.5
  enddo
 enddo

 sumprob=0.
 write(1,*) nat
 write(1,903) Esamp,sumprob,Esamp,Esamp," Traj ",itraj
 write(2,"(I6.6,A,3(x,E20.10e3))") itraj," Traj ",sumprob,Esamp,Ek
 do i=1,nat
  write(1,901) sat(i),(geom(i,j)*.5292,j=1,3),(vel(i,j)*1000.,j=1,3)
  write(2,900) (geom(i,j),j=1,3)
 enddo
 do i=1,nat
  write(2,900) (vel(i,j),j=1,3)
 enddo
 
enddo

900 format (1000(x,E20.10e3))
901 format (A2,1000(x,E20.10e3))
903 format (4(x,F20.10),A,I6.6)
end
   
