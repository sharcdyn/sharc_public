! Program to generate initial conditions in files geoms.init and
! xyz file (geoms.xyz)


program initial_conditions
implicit none
real*8, allocatable :: x(:), p(:), mass(:), coord0(:,:),coord(:,:),dx(:),dp(:)
real*8, allocatable :: massred(:), frec (:), modno(:,:,:), veloc(:,:)
real*8, allocatable :: probmod(:), probmaxmod(:)
real*8 rand1, rand2, rand3
integer traj, ntraj, i,j,k,l, nmod, nat, np,init, inicio, nseed
integer, allocatable :: seed(:)
integer, allocatable :: noat(:), mode(:), quantum(:)
character*2, allocatable :: sat(:)
real*8 prob, kk, delta1, delta2, xmax
logical writ, nm, hessian, accepted, ex, res
complex*16, allocatable :: wf(:,:)
real*8, allocatable :: r(:,:)
real*8 ZPE, total, kintot, pottot, kintot2

real*8 Ekin,Epot, Ekin2, kkat
character*100 fmtd

write(6,*) "Number of atoms, number of modes (total number of them) and number of points"
read(5,*) nat,nmod, np
write(6,*) nat,nmod,np
allocate (mode(nmod), quantum(nmod))
if (nmod.ne.0) then
 write(6,*) "Deltas required (grid and sampling), usual good numbers: 70. 10."
 read(5,*) delta1,delta2
 write(6,*) delta1,delta2
! Delta1 -> grid
! Delta2 -> sampling

 write(6,*) "Initial trajectory and final trajectory"
 read(5,*) init,ntraj
 write(6,*) init,ntraj
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

 write(6,*) "Normal modes to scan"
 read(5,*) (mode(i),i=1,nmod)
 write(6,*) "Quantum for the normal modes"
 read(5,*) (quantum(i),i=1,nmod)
 do i=1,nmod
  write(6,*) mode(i),quantum(i)
 enddo
endif

allocate ( massred (nmod), frec(nmod), modno(nmod,nat,3) )
allocate ( x(nmod), p(nmod), mass(nat), coord(nat,3), noat(nat) )
allocate ( coord0(nat,3),veloc(nat,3), sat(nat) )
allocate ( dx(nmod),dp(nmod), probmod(nmod), probmaxmod(nmod) )
allocate ( wf(nmod,np),r(nmod,np) )


!   !!! Reading equilibrium coordinates

open(1,file="geom")
do i=1,nat
 read(1,*) sat(i), kkat, coord0(i,1), coord0(i,2), coord0(i,3), mass(i)
 noat(i)=int(kkat+0.5)
enddo
close(1)


inquire (file="nm",exist=nm)
inquire (file="hessian.dat",exist=hessian)
if(nm) then
 open(1,file="nm")
 do i=1,nmod
  kk=0.d0
  read(1,*) frec(i)
   frec(i)=frec(i)/2.194746E5
  do j=1,nat
   read(1,*) modno(i,j,1), modno(i,j,2) ,modno(i,j,3)
  enddo
 enddo
 close (1)
elseif(hessian) then
 call hessiano (nat,sat,coord0,nmod,mode, modno,frec, mass,noat)
else
 stop 'nm or hessian.dat are required'
endif

if (nmod.eq.0) stop "Finishing here, no sampling required"

open(11,file="freqs.xyz")

!Writting the eigenvectors for each normal mode on the file freqs.xyz
do i=1,nmod
 write(11,*) nat
 write(11,*) frec(i)
 do j=1,nat
  write(11,801) sat(j), coord0(j,1)/1.889726,coord0(j,2)/1.889726,coord0(j,3)/1.889726,&
   modno(i,j,1)/1.889726, modno(i,j,2)/1.889726, modno(i,j,3)/1.889726
 enddo
enddo
close(11)

801 format(A2,x,6(E20.8,x))



ZPE=0.d0
Epot=0.d0
do i=1,nmod
 if (frec(i).le.0.d0) then
  frec(i)=-frec(i)
  write(6,*) "WARNING, frec(",i,") is imag"
 else
  ZPE=ZPE+ (0.5)*frec(i)
  Epot=Epot+ (float(quantum(i))+0.5)*frec(i)
 endif
 massred(i)=0.d0
 do j=1, nat
  do k=1, 3
   massred(i)=massred(i)+( modno(i,j,k)**2 /mass(j)/1822.888 )
  enddo
 enddo
 massred(i)= 1.d0/ massred(i)
enddo
write(6,"(A,3(X,F20.10))") "ZPE (H,cm-1,eV)",ZPE,ZPE*2.194746E5,ZPE*27.211
write(6,"(A,3(X,F20.10))") "Harm energy (H,cm-1,eV)",Epot,Epot*2.194746E5,Epot*27.211


!   !!! Recalculating normal modes
do i=1,nmod
 do j=1,nat
  do k=1,3
   modno(i,j,k)=modno(i,j,k)*dsqrt( massred(i)/( mass(j)*1822.888 ) )
  enddo
 enddo
 dx(i)= ( massred(i)*frec(i) )**(-0.5)
 dp(i)=1.d0/(dx(i))
 write(6,"(a,I6.6,3(x,F30.10))") "Mode, FC, freq and mass: ",&
  mode(i), frec(i)*frec(i)*massred(i), frec(i)*219500, massred(i)/1822.888
 write(6,*) "Widths ",dx(i),dp(i)
 write(6,*) "Scanning x ",-delta2*dx(i)/2.,delta2*dx(i)/2.
 write(6,*) "Scanning p ",-delta2*dp(i)/2.,delta2*dp(i)/2.
enddo

!call srand(seed(1))
call random_seed (put=seed)

!   !!! Begin!
open(11,file="geoms.init")
open(21,file="geoms.nm")
open(22,file="moms.nm")
open(13,file="geoms.xyz")
open(14,file="harm_energy.out")

x(:)=0.d0
p(:)=0.d0

do i=1,nmod
 probmaxmod(i)=0.d0
 write(6,"(A,I6.6,A,I3)") "Checking maximum prob mode ",i," v=",quantum(i)
 do j=1,np
!   !!! We scan only in x, because probability is symmetric for x and p in a Wigner distribution,
!   !!! it doesn't matter the quantum number
  x(i)=-delta2*dx(i)/2.+(j-1)*(delta2*dx(i))/float(np-1)
  call probab (x(i),0.d0,np,frec(i),massred(i), delta1, &
  .false.,mode(i),prob,wf(i,:),r(i,:),.true., quantum(i))
!  write(98,*) x(i),prob
  if (prob.ge.probmaxmod(i)) then
   probmaxmod(i)=prob
   xmax=x(i)
  endif
 enddo
 x(i)=xmax
! stop
 probmaxmod(i)=1.1*probmaxmod(i)
 call probab (x(i),p(i),np,frec(i),massred(i), delta1, &
 .true.,mode(i),prob,wf(i,:),r(i,:),.true., quantum(i))
 write(6,*) "Xmax and Maximum prob ",xmax,probmaxmod(i)
enddo



! Number of atoms at the beginning to check
write(11,*) nat
do i=1,nat
 write(11,901) sat(i), float(noat(i)),mass(i)
enddo


write(13,*) nat
kk=0.
i=0
write(13,903) kk, kk,kk,ZPE, " Traj ",i
do j=1,nat
 write(13,901) sat(j),(coord0(j,k)*0.5292,k=1,3),kk,kk,kk
enddo

traj=init-1
total=0.d0
pottot=0.d0
kintot=0.d0
kintot2=0.d0
do while(traj.lt.ntraj)
 write(6,*) "Doing Traj", traj+1
 do i=1,nmod
  accepted=.true.
  do while (accepted)
   !rand1=rand()
   call random_number (rand1)
   x(i)=(rand1-0.5)*dx(i)*delta2
   !rand2=rand()
   call random_number (rand2)
   p(i)=(rand2-0.5)*dp(i)*delta2
  !enddo

   call probab(x(i),p(i),np,frec(i),massred(i),delta1,.false.,mode(i) &
   ,prob,wf(i,:),r(i,:),.false., quantum(i))
   probmod(i)=prob/probmaxmod(i)

   if (probmod(i).gt.1.d0) then
    write(6,*) "WARNING Mode prob. and maximum",probmod(i), probmaxmod(i)
    stop
   endif

   !rand3=rand()
   call random_number (rand3)
   if (rand3.le.probmod(i)) accepted=.false.
  enddo
 enddo
 traj=traj+1
 do j=1,nat
  do k=1,3
   coord(j,k)=coord0(j,k)
   veloc(j,k)=0.d0
   do l=1,nmod
    coord(j,k)=coord(j,k)+x(l)*modno(l,j,k)
    veloc(j,k)=veloc(j,k)+p(l)*modno(l,j,k)/massred(l)
   enddo
  enddo
 enddo
 Ekin=0.d0
 do j=1,nat
  do k=1,3
   Ekin=Ekin+mass(j)*1822.888*veloc(j,k)**2/2.
  enddo
 enddo
 Ekin2=0.d0
 Epot=0.d0

 do l=1,nmod
  Epot=Epot+x(l)**2*frec(l)**2*massred(l)/2.
  Ekin2=Ekin2+p(l)**2/2.d0/massred(l)
 enddo

 total=total+Epot+Ekin
 pottot=pottot+Epot
 kintot=kintot+Ekin
 kintot2=kintot2+Ekin2

 !Geoms.init
 write(11,"(I6.6,A,3(x,E20.10e3))") traj," Traj ",Epot,Ekin,Ekin2
 !Geoms.xyz
 write(13,*) nat
 write(13,903) Epot+Ekin, Epot,Ekin,ZPE," Traj ",traj
 write(14,900) Epot,Ekin,ZPE


 do j=1,nat
  write(11,900) coord(j,1),coord(j,2),coord(j,3)
  write(13,901) sat(j),(coord(j,k)*0.5292,k=1,3),(veloc(j,k)*1000.,k=1,3)
 enddo
 do j=1,nat
  write(11,900) veloc(j,1),veloc(j,2),veloc(j,3)
 enddo

 write(21,900) (x(j),j=1,nmod)
 write(22,900) (p(j),j=1,nmod)
 call now(6)
enddo

write(6,*) "Total Energ.", total/float(ntraj)*219500, "ZPE:",ZPE*219500
write(6,*) "Kinetic=", kintot/float(ntraj)*219500, "Weighted=", kintot2/float(ntraj)*219500
write(6,*) "Potential=", pottot/float(ntraj)*219500

close(11)
close(13)

900 format (1000(x,E20.10e3))
901 format (A2,1000(x,E20.10e3))
903 format (4(x,F20.10),A,I6.6)

end program initial_conditions

!!! Subroutine to read a triangular hessian matrix (molpro output type) and to obtain Normal modes
!!! and Frequencies (eigenvalues and eigenvectors)
!!! Needed things:
!!! -diag.f

subroutine hessiano (n,sat,coord0,nmod,mode,modno1, frec1, mass,noat)
implicit none
integer, intent(in) :: n,nmod
character*2, intent(in) :: sat(n)
real*8, intent(in) :: mass(n),coord0(n,3)
integer, intent(in) :: mode(nmod),noat(n)

real*8, intent(out) :: modno1(nmod,n,3),frec1(nmod)

real*8 Hess(3*n,3*n),Hessd(3*n),x(n,3),kk(3*n), frec(3*n), evect(3*n,3*n)
real*8 Kfor(3*n), modno(3*n,n,3)
integer n3, inicio
real*8 popo,massred

integer i,j,k,l,m

n3=3*n

!!! Creating virtual masses for 3N-atoms
k=0
do j=1,n
 do i=1,3
  k=k+1
  kk(k)=mass(j)
 enddo
enddo

!!! Reading the Hessian triangular Matrix
open(1,file="hessian.dat")
read (1,*) ( ( Hess(i,j) ,j=1,i ),i=1,n3 )
close(1)

!!! Mirroring the Hessian upper part
do i=1,n3
 do j=1,i
  Hess(j,i)=Hess(i,j)
 enddo
enddo

!!! Dividing over rootsquare of mass
do i=1,n3
 do j=1,n3
  Hess(i,j)=Hess(i,j) / dsqrt(kk(i)) / dsqrt(kk(j))
  Hess(i,j)=Hess(i,j) / 1822.888
 enddo
enddo

!!! Diagonalizing the Hessian to obtain eigenvalues and eigenvectors
call diag (Hess,n3,Hessd)

!!! Doing the rootsquare of the diagonal terms, bewaring of imaginary frecuencies
do i=1,n3
 if (Hessd(i).ge.0) then
  Hessd(i)=dsqrt(Hessd(i))
 else
  Hessd(i)=-dsqrt(-Hessd(i))
 endif
enddo

do i=1,n3
 frec(i)=Hessd(i)
 l=0
 do j=1,n
  do k=1,3
   l=l+1
   modno(i,j,k)=hess(l,i)
  enddo
 enddo
enddo

!!! Reordering frequencies and eigenvectors

do i=2,n3
 do j=1,i
  if (frec(i).lt.frec(j)) then
   popo=frec(j)
   frec(j)=frec(i)
   frec(i)=popo
   do k=1,n
    do l=1,3
     popo=modno(j,k,l)
     modno(j,k,l)=modno(i,k,l)
     modno(i,k,l)=popo
    enddo
   enddo
  endif
 enddo
enddo
write(6,*) "Frequencies obtained in the diagonalization "
do i=1,n3
 write(6,*) i,frec(i)*219500
enddo

open(1,file="freqs_all.xyz")
do i=1,n3
 massred=0.
 do j=1,n
  do k=1,3
   massred=massred+( modno(i,j,k)**2 /mass(j)/1822.888 )
  enddo
 enddo
 massred=1./massred
 popo=frec(i)**2*massred
 if (frec(i).le.0) popo=-popo
 write(1,*) n
 write(1,*) frec(i)*219475,massred/1822.888,popo
 do j=1,n
  write(1,"(A2,6(x,E25.15e3))") &
   sat(j),(coord0(j,k)*.5292,k=1,3),(modno(i,j,k)*.5292,k=1,3)
 enddo
enddo
close(1)
open(2,file="freqs_all.molden")
call molden (n,sat,noat,coord0,n3,frec,modno,2)
close(2)

!! Taking only the ones selected by the user
if (nmod.ne.0) then
 do i=1,nmod
  frec1(i)=frec(mode(i))
  do j=1,n
   do k=1,3
    modno1(i,j,k)=modno(mode(i),j,k)
   enddo
  enddo
 enddo
 open(2,file="freqs.molden")
 call molden (n,sat,noat,coord0,nmod,frec1,modno1,2)
 close(2)
endif

end subroutine


subroutine probab(x,p,np,frec,massred,delta,writ,mode,prob,wf,r,create, quantum)
implicit none
integer np, mode, quantum
real*8 x,p, prob, frec, massred
real*8 delta
logical writ,create

real*8 dr,dk,dx
real*8 r(np),k(np)
complex*16 wf(np),wfp(np)

integer i,j
character*20 fich
real*8 kk, hermit

if (create) then
 dx=(massred*frec)**(-.5d0)
 dr=delta*dx/(np-1)
 dk=2*dacos(-1.d0)/np/dr
 do i=1, np
  r(i)=-delta/2.d0*dx+(dr*(i-1))
  call hermite (quantum,r(i)/dx,hermit)
  wf(i)= dacos(-1.d0)**(0.25)* dexp( -massred*frec*r(i)**2/2 )*hermit
 enddo
endif
!write(99,909) x,p, prob
call wigner(x,p,np,r,wf,prob)

if (writ) then
! Evaluating the FT
 do i=1,np
  k(i)=-dacos(-1.d0)/dr+dk
  k(i)=k(i)+(i-1)*dk
  wfp(i)=dcmplx(0.d0,0.d0)
  do j=1,np
   kk=k(i)*r(j)
   wfp(i)=wfp(i)+wf(j)*cdexp(dcmplx(0.d0,-kk))
  enddo
 enddo
 write(6,*) "NM ",mode,frec*219474.6,massred/1822.888
 write(fich,'("nm",I6.6)') mode
 open(1,file=fich)
 do i=1,np
  write(1,900) r(i),dreal(wf(i)),dimag(wf(i)),k(i),dreal(wfp(i)),dimag(wfp(i))
 enddo
 close(1)
endif


900 format(1000(x,E20.10e3))
909 format(6(F10.6))
return
end



subroutine wigner(x, p, np, r, wf, prob)
  implicit none
  integer np
  real*8 r(np),x,p
  complex*16 wf(np)
  real*8 prob
  complex*16 probability

  integer i, ii
  complex*16 kk2, val1, val2
  real*8 lr, lr1, lr2, dr, valdata, rmax, rmin

  rmax=r(np)
  rmin=r(1)
  dr=(rmax-rmin)/(np-1)
  probability=dcmplx(0.d0,0.d0)
  do ii=-2*np,2*np
   lr=ii*dr
   !sum, value of the wavefunction in lr1
   lr1= x+lr
   !difference
   lr2= x-lr
   call valr(np,r,wf,lr1,val1)
   call valr(np,r,wf,lr2,val2)
   kk2=cdexp(dcmplx(0.d0,2.d0*p*lr))
   kk2=dconjg(val1)*val2*kk2
   probability=probability+kk2
  enddo
  probability=1.d0/dacos(-1.d0)*probability*dr

  prob=real(probability)
return
end subroutine


subroutine valr(np,r,wf,lr,valdata)
  !Gives the value of the wavefunction at lr
  implicit none
  integer np
  complex*16 wf(np),valdata
  real*8 r(np),lr

  integer pointbefore, pointafter
  integer i
  complex*16 slope

  real*8 minimum
  integer imin

  do i=1,np
   if (r(i).le.lr) then
    pointbefore=i
    pointafter=i+1
   endif
  enddo
  if (lr.le.r(1)) then
   pointbefore=1
   pointafter=2
  endif
  if (lr.ge.r(np)) then
   pointbefore=np-1
   pointafter=np
  endif
  slope=(wf(pointafter)-wf(pointbefore))/(r(pointafter)-r(pointbefore))
  valdata=wf(pointbefore)+slope*(lr-r(pointbefore))

!minimum=2*(r(np)-r(1))/(np-1)
!do i=1,npdacos
! if (dabs(r(i)-lr).le.minimum) then
!  minimum=dabs(r(i)-lr)
!  imin=i
! endif
!enddo
!
!if (lr.lt.r(1)) imin=1
!if (lr.gt.r(np)) imin=np
!valdata=wf(imin)

return
end subroutine


subroutine hermite (n,x, hermit)
implicit none
integer n, i
real*8 x, her(0:n), hermit

if (n.lt.0) stop "Not possible to create a negative Hermite polynomial"
her(0)=1.d0
if (n.ge.1) her(1)=2*x
if (n.gt.1) then
 do i=2,n
  her(i)= (x*2.d0*her(i-1)) - (2.d0*float(i-1)*her(i-2))
 enddo
endif

hermit=her(n)

return
end subroutine

subroutine now(unit)
implicit none
character*8 date
character*10 time
character*5 zone
integer values(8)

integer unit

call DATE_AND_TIME(date, time, zone, values)
write(6,*) date," ",time

return
end subroutine

subroutine molden (nat,sat,noat,coord0,nmod,frec,modno,un)
 integer, intent(in) :: nat,nmod,un
 integer, intent(in) :: noat(nat)
 character*2,intent(in) :: sat(nat)
 real*8, intent(in) :: coord0(nat,3),modno(nmod,nat,3),frec(nmod)

 integer i,j
 
 write(un,"(A)") "[Molden Format]" 
 
 write(un,"(A)") "[Atoms] AU"
 do i=1,nat
  write(un,"(A2,2(x,I5),3(x,E20.10e3))") &
   sat(i), i, noat(i), coord0(i,1), coord0(i,2), coord0(i,3)
 enddo
 
 write(un,"(A)") " [FREQ]"
 do i=1,nmod
  write(un,*) frec(i)*2.194746E5
 enddo
 
 write(un,"(A)") " [FR-COORD]"
 do i=1,nat
  write(un,"(A2,3(x,E20.10e3))") &
   sat(i), coord0(i,1), coord0(i,2), coord0(i,3)
 enddo

 write(un,"(A)") " [FR-NORM-COORD]"
 do i=1,nmod
  write(un,"(A,I10)") "Vibration                     ", i 
  do j=1,nat
   write(un,*) modno(i,j,1), modno(i,j,2), modno(i,j,3)
  enddo
 enddo
end subroutine
