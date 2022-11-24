program initial_conditions_1
!!! Program to calculate a set of conditions following a 1D Wigner distribution
!!! for an arbitrary potential described in a grid
!!! The wavefunction is obtained with the FGH method
implicit none
integer np,imod,ieig
real*8 probmin,mass,probmax,prob,frec
real*8, allocatable :: wf(:,:),PES(:),x(:),eig(:)

character*100 fich
integer i,j,ntraj
real*8 pos,mom,dx,dp
real*8 Epot,Ekin,Epav,Ekav
real*8 x0,xf,p0,pf

integer nseed,ntrials
integer, allocatable :: seed(:)
logical res,ex
real*8 rando

write(6,*) "Normal mode that you want to sample (info is gonna be taken from the nm files) and eigenfunction (ground=0)"
read(5,*) imod,ieig
ieig=ieig+1
write(6,*) "Minimum probability to consider in the sampling and number of trajectories"
read(5,*) probmin,ntraj

write(fich,"(A,I0,A)") "nm_",imod,"-info.xyz"
write(6,*) "Reading the mass from ",trim(fich)
open(1,file=fich,status="old",iostat=i)
if (i.ne.0) stop "Info file does not exist"
read(1,*) i
read(1,*) np,frec,mass
close(1)
write(6,*) "Reduced mass of the normal mode ",imod,mass/1822.888
write(6,*) "Number of grid points ",np


write(fich,"(A,I0,A)") "nm_",imod,"-harmPES.dat"
open(1,file=fich,status="old",iostat=i)
if (i.ne.0) stop "harmPES.dat file does not exist"
allocate(x(np),PES(np),wf(np,np),eig(np))
do i=1,np
 read(1,*) x(i),PES(i)
enddo

write(fich,"(A,I0,A)") "nm_",imod,"-PES.dat"
open(1,file=fich,status="old",iostat=i)
if (i.eq.0) then 
 write(6,*) "File with a custom potential does exist ",trim(fich)
 write(6,*) "Replacing the harmonic potential"
 call read_spl_PES(np,x,1,PES)
 write(fich,"(A,I0,A)") "nm_",imod,"-newPES.dat"
 write(6,*) "Copying the new potential to ",trim(fich)
 open(1,file=fich)
 do i=1,np
  write(1,"(2(X,E30.20))") x(i),PES(i)
 enddo
 close(1)
endif

write(6,*) "Calculating the eigenfunctions of the potential"
write(6,*) "Using the Fourier Grid Hamiltonian method:"
write(6,*) "J. Chem. Phys. 91, 3571 (1989)"
call FGH1D(np,mass,x,PES,wf,eig)
write(6,*) "Eigenfunctions up to the selected one ",ieig-1
do i=1,ieig
 write(6,*) "Eigenvalue ",i-1,eig(i),eig(i)*2.194746E5
enddo
write(fich,"(A,I0,A)") "nm_",imod,"-wf.out"
open(1,file=fich)
call writewf(np,x,wf(:,ieig),1)
close(1)

write(fich,"(A,I0,A)") "nm_",imod,"-eig.dat"
write(6,*) "Writing all the eigenvalues from the FGH in ",trim(fich)
open(1,file=fich)
do i=1,np
 write(1,*) eig(i)
enddo
close(1)

call check_limits(np,x,wf(:,ieig),probmin,probmax,x0,xf,p0,pf)

!!! To avoid interface, the seed generation is the main program
write(fich,"(A,I0,A)") "nm_",imod,"-seeds.dat"
open(1,file=fich,status="old",iostat=i)
if (i.ne.0) then
 write(6,*) "Creation of random seeds in file ",trim(fich)
 open(1,file=fich)
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
 open(1,file=fich,status="old")
endif
read(1,*) nseed
allocate(seed(nseed))
read(1,*) (seed(i),i=1,nseed)
close(1)
call random_seed(put=seed)

write(fich,"(A,I0,A)") "nm_",imod,"-trajs.out"
write(6,*) "Starting the mining in file ",trim(fich)
open(1,file=fich)
write(1,"(A,1(I0,x),E30.20e3)") "# ",ntraj,eig(ieig)
Epav=0.
Ekav=0.
dx=( mass*frec )**(-0.5)
dp=1./dx
do i=1,ntraj
 res=.true.
 ntrials=0
 do while (res)
  ntrials=ntrials+1
  call random_number(rando)
  pos=x0+(xf-x0)*rando
  call random_number(rando)
  mom=p0+(pf-p0)*rando
  call probab (pos,mom,np, wf(:,ieig),x,prob)
  call random_number(rando)
  if (prob/probmax*.9.ge.1) then
   write(6,*) "Obtained a probability bigger than max ",prob,probmax
   stop "Problems with max prob"
  endif
  if (rando.le.prob/probmax*.9) res=.false.
 enddo
 call dvalr(np,x,PES,pos,Epot)
 Ekin=mom*mom/2/mass
 write(1,"(4(X,E30.20e3))") pos,mom,Epot,Ekin
 write(6,"(A,I0,2(x,F10.2),X,I0)") "Accepted traj ",i,(Epot+Ekin)*2.194746E5,eig(ieig)*2.194746E5,ntrials
 Epav=Epav+Epot
 Ekav=Ekav+Ekin
enddo
Ekav=Ekav/float(ntraj)
Epav=Epav/float(ntraj)
write(6,*) "Done, average energy :",(Epav+Ekav)*2.194746E5
write(6,*) "Expected energy and half:", eig(ieig)*2.194746E5,eig(ieig)*2.194746E5/2.
write(6,*) "Potential and kinetic average :", Epav*2.194746E5,Ekav*2.194746E5
write(6,*) "Virial (2*kin/Etot): ",2*Ekav/(Ekav+Epav)

end

subroutine FGH1D(np,mass,x,PES,wf,eig)
implicit none
integer, intent(in) :: np
real*8, intent(in) :: mass,PES(np),x(np)
real*8, intent(out) :: wf(np,np),eig(np)

real*8 dx,tl,dk
integer i,j,l

dx=(x(np)-x(1))/float(np-1)
!! Evaluating the kinetic part
wf=0.
do l=1,(np-1)/2
!write(6,*) "Doing ",l,(np-1)/2
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(tl)
 do i=1,np
  do j=1,np
   tl=float(l)*2.*dacos(-1.d0)*float(i-j)/float(np)
   tl=cos(tl)*(float(l)*dacos(-1.d0)/float(np)/dx)**2
   wf(j,i)=wf(j,i)+tl
  enddo
 enddo
!$OMP END PARALLEL DO
enddo

!$OMP PARALLEL DO COLLAPSE(2)
do i=1,np
 do j=1,np
  wf(j,i)=wf(j,i)*2/float(np) *2./mass
 enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,np
 wf(i,i)=wf(i,i)+PES(i)
enddo
!$OMP END PARALLEL DO

!write(6,*) "Diagonalization"
call diag(wf,np,eig)

end subroutine

subroutine probab(x,p,np,wf,r,prob)
implicit none
integer, intent(in) :: np
real*8, intent(in) :: x,p,wf(np),r(np)
real*8, intent(out) :: prob

complex*16 wfr(np),wfp(np)

wfr=wf
call wigner(x,p,np,r,wfr,prob)

end subroutine probab

subroutine writewf(np,r,wf,un)
implicit none
integer, intent(in) :: np,un
real*8, intent(in) :: r(np),wf(np)

complex*16 wfr(np),wfp(np)
real*8 k(np),dr,dk,kk
integer i,j

wfr=wf
dr=(r(np)-r(1))/float(np-1)
dk=2*dacos(-1.d0)/float(np)/dr
! Evaluating the FT
do i=1,np
 k(i)=-dacos(-1.d0)/dr+dk
 k(i)=k(i)+(i-1)*dk
 wfp(i)=dcmplx(0.d0,0.d0)
 do j=1,np
  kk=k(i)*r(j)
  wfp(i)=wfp(i)+wfr(j)*cdexp(dcmplx(0.d0,-kk))
 enddo
enddo
do i=1,np
 write(un,900) r(i),dreal(wfr(i)),dimag(wfr(i)),k(i),dreal(wfp(i)),dimag(wfp(i))
enddo


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
  complex*16 probability(-2*np:2*np)

  integer i, ii
  complex*16 kk2, val1, val2
  real*8 lr, lr1, lr2, dr, valdata, rmax, rmin

  rmax=r(np)
  rmin=r(1)
  dr=(rmax-rmin)/(np-1)

!$OMP PARALLEL DO PRIVATE(lr,lr1,lr2,val1,val2,kk2)
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
   probability(ii)=kk2
  enddo
!$OMP END PARALLEL DO

  prob=0.
  do ii=-2*np,2*np
   prob=prob+real(probability(ii))
  enddo
  prob=1.d0/dacos(-1.d0)*prob*dr

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

subroutine dvalr(np,r,PES,lr,valdata)
  !Gives the value of the potential at lr
  implicit none
  integer, intent(in) :: np
  real*8, intent(in) :: PES(np),r(np),lr
  real*8, intent(out) :: valdata

  integer pointbefore, pointafter
  integer i
  real*8 slope

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
  slope=(PES(pointafter)-PES(pointbefore))/(r(pointafter)-r(pointbefore))
  valdata=PES(pointbefore)+slope*(lr-r(pointbefore))

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

subroutine read_spl_PES(np,x,un,PES)
implicit none
integer, intent(in) :: un,np
real*8, intent(in) :: x(np)
real*8, intent(inout) :: PES(np)

integer nred,i,j
real*8, allocatable :: redPES(:),redx(:),a(:),b(:),c(:),d(:)
real*8 kk

nred=0
do 
 read(1,*,end=1) kk
 nred=nred+1
enddo
1 continue
rewind(1)
write(6,*) "Custom potential points: ",nred
allocate(redPES(nred),redx(nred))
do i=1,nred
 read(1,*) redx(i),redPES(i)
enddo
close(1)
allocate(a(nred),b(nred),c(nred),d(nred))
write(6,*) "Splinning to ",np
call cspline(nred,nred,redx,redPES,a,b,c,d)
do i=1,np
 if (x(i).le.redx(1)) then
  PES(i)=redPES(1)+(redPES(2)-redPES(1))/(redx(2)-redx(1))*(x(i)-redx(1))
 else if (x(i).ge.redx(nred)) then
  PES(i)=redPES(nred)+(redPES(nred)-redPES(nred-1))/(redx(nred)-redx(nred-1))*(x(i)-redx(nred))
 else
  do j=1,nred-1
   if (x(i).ge.redx(j).and.x(i).lt.redx(j+1)) then
    kk=x(i)-redx(j)
    PES(i)=d(j)+c(j)*kk+b(j)*kk*kk+a(j)*kk*kk*kk
   endif
  enddo
 endif
enddo

kk=PES((np+1)/2)
PES=PES-kk

end subroutine read_spl_PES
 

subroutine check_limits(np,x,wf,probmin,probmax,x0,xf,p0,pf)
implicit none
integer, intent(in) :: np
real*8, intent(in) :: x(np),wf(np),probmin
real*8, intent(out) :: probmax,x0,xf,p0,pf

integer i
real*8 prob,xmax
real*8 dx,dp
real*8 prob0,probf
real*8 :: zero=0.

!!! Looking for the maximum
probmax=0.
do i=1,np
 call probab (x(i),zero,np, wf,x,prob)
 if (prob.ge.probmax) then
  xmax=x(i)
  probmax=prob
 endif
enddo
write(6,*) "Maximum probability at ",xmax,probmax

dx=(x(np)-x(1))/float(np-1)
x0=x(1)-dx
prob0=0.
do while (prob0.le.probmin)
 x0=x0+dx
 call probab (x0,zero,np, wf,x,prob)
 prob0=prob/probmax
enddo

prob=0.
xf=x(np)+dx
probf=0.
do while (probf.le.probmin)
 xf=xf-dx
 call probab (xf,zero,np, wf,x,prob)
 probf=prob/probmax
enddo
write(6,*) "Sampling of x between: ",x0,xf,prob0,probf

dp=2*acos(-1.)/(dx*np)
! pmin=-(np-1)/2.*dp
p0=-dp*float(np)/2.
!! Let's remove for the sampling 1/10 of the grid in momentum
p0=p0+float(np)/10.*dp
prob0=0.
do while (prob0.le.probmin)
 p0=p0+dp
 call probab (zero,p0,np, wf,x,prob)
 prob0=prob/probmax
enddo
! pmin=-(np-1)/2.*dp
pf=+dp*float(np)/2.
!! Let's remove for the sampling 1/10 of the grid in momentum
pf=pf-float(np)/10.*dp
probf=0.
do while (probf.le.probmin)
 pf=pf-dp
 call probab (zero,pf,np, wf,x,prob)
 probf=prob/probmax
enddo
write(6,*) "Sampling of p between: ",p0,pf,prob0,probf

end subroutine check_limits
