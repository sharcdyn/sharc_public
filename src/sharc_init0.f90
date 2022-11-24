! Program to generate initial conditions in files geoms.init and
! xyz file (geoms.xyz)


program initial_conditions_0
implicit none
real*8, allocatable :: mass(:), coord0(:,:),x(:),coord(:,:),vel(:,:)
real*8, allocatable :: massred(:), frec (:), modno(:,:,:)
integer, allocatable :: noat(:),mode(:)
character*2, allocatable :: sat(:)
integer np,nmod,nat
real*8 delta1, delta2
logical hessian
real*8 ZPE, total, kintot, pottot, kintot2

real*8 Ekin,Epot, Ekin2, kkat
character*100 fmtd
character*100000 fich

real*8 kk
integer k,j,i,l

write(6,*) "Number of atoms and number of points (odd number, 0 if only freqs)"
read(5,*) nat, np
write(6,*) nat,np
if (np.gt.0) then
 if (int(np/2).eq.int(float(np)/2.+.6) ) stop "Not an odd number"
endif
write(6,*) "File with the force constant matrix"
read(5,*) fich
open(2,file=fich,status="old",iostat=i)
if (i.ne.0) stop "File with force constant matrix does not exist"
if (np.ne.0) then
 write(6,*) "Delta required (grid), usual good number is 70."
 read(5,*) delta1
 write(6,*) delta1
! Delta1 -> grid
endif

!! To recycle part of the code let's create a mode vector 1,2,3,4...
nmod=3*nat
allocate(mode(nmod))
do i=1,nmod
 mode(i)=i
enddo

allocate ( massred (nmod), frec(nmod), modno(nmod,nat,3), x(nmod) )
allocate ( mass(nat), noat(nat) )
allocate ( coord0(nat,3), sat(nat) )


!   !!! Reading equilibrium coordinates

open(1,file="geom",status="old",iostat=i)
if (i.ne.0) stop "geom file must exist"
do i=1,nat
 read(1,*) sat(i), kkat, coord0(i,1), coord0(i,2), coord0(i,3), mass(i)
 noat(i)=int(kkat+0.5)
enddo
close(1)


call hessiano (nat,sat,coord0,nmod,mode, modno,frec, mass,noat,2)



ZPE=0.d0
do i=1,nmod
 if (frec(i).le.0.d0) then
  write(6,*) "WARNING, frec(",i,") is imag, not included in ZPE"
 else
  ZPE=ZPE+ (+0.5)*frec(i)
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

open(11,file="freqs.xyz")

!Writting the eigenvectors for each normal mode on the file freqs.xyz
do k=1,nmod
 i=mode(k)
 write(11,*) nat
 write(11,*) frec(i)*2.194746E5,massred(i)/1822.888,frec(i)
 do j=1,nat
  write(11,801) sat(j), coord0(j,1)/1.889726,coord0(j,2)/1.889726,coord0(j,3)/1.889726,&
   modno(i,j,1)/1.889726, modno(i,j,2)/1.889726, modno(i,j,3)/1.889726
 enddo
enddo
close(11)

801 format(A2,x,6(E20.8,x))

if (np.gt.0) then
 allocate(coord(nat,3),vel(nat,3))
 do i=1,3*nat
  write(6,*) "Writing stuff for nm ",i
!!!!!! Recalculating normal modes
  do j=1,nat
   do k=1,3
    modno(i,j,k)=modno(i,j,k)*dsqrt( massred(i)/( mass(j)*1822.888 ) )
   enddo
  enddo
!!!!!! File with the nm information
  write(fich,"(A,I0,A)") "nm_",i,"-info.xyz" 
  open(1,file=fich)
  write(1,*) nat
  write(1,"(I0,2(x,E30.20e3),A)") np,frec(i),massred(i)," Number of points, freq and reduced mass"
  do j=1,nat
   write(1,801) sat(j), coord0(j,1)/1.889726,coord0(j,2)/1.889726,coord0(j,3)/1.889726,&
   modno(i,j,1)/1.889726, modno(i,j,2)/1.889726, modno(i,j,3)/1.889726
  enddo
  close(1)
!!!!!! File with the harmonic grid
  x=0.
  kk=( massred(i)*abs(frec(i)) )**(-0.5)
  write(fich,"(A,I0,A)") "nm_",i,"-harmPES.dat"
  open(1,file=fich)
! write(fich,"(A,I0,A)") "nm_",i,"-grid.xyz"
! open(2,file=fich)
  do j=1,np
   x(i)=-delta1*kk/2.+(j-1)*(delta1*kk)/float(np-1)
   Epot=massred(i)*frec(i)**2/2.*x(i)*x(i)
   write(1,"(2(x,E30.20e3))") x(i),Epot
!  call creategeom(nat,nmod,coord0,modno,x,coord,vel)
!  write(fich,*) Epot,(x(k),k=1,nmod)
!  call writegeom(nat,sat,coord,2,fich)
  enddo
  close(1)
! close(2)
 enddo
endif

end program initial_conditions_0

!!! Subroutine to read a triangular hessian matrix (molpro output type) and to obtain Normal modes
!!! and Frequencies (eigenvalues and eigenvectors)
!!! Needed things:
!!! -diag.f

subroutine hessiano (n,sat,coord0,nmod,mode,modno1, frec1, mass,noat,un)
implicit none
integer, intent(in) :: n,nmod,un
character*2, intent(in) :: sat(n)
real*8, intent(in) :: mass(n),coord0(n,3)
integer, intent(in) :: mode(nmod),noat(n)

real*8, intent(out) :: modno1(nmod,n,3),frec1(nmod)

real*8 Hess(3*n,3*n),Hessd(3*n),x(n,3),kk(3*n), frec(3*n)
real*8 modno(3*n,n,3)
integer n3, inicio
real*8 popo,massred
real*8 coord(n,3),cm(3),totmass
real*8 proj(6,3*n),proj2(3*n,3*n)
real*8 Hess2(3*n,3*n)

integer i,j,k,l,m,ii,ij

n3=3*n

!!! Creating virtual masses for 3N-atoms and cm
k=0
cm=0.d0
totmass=0.d0
do j=1,n
 totmass=totmass+mass(j)
 do i=1,3
  cm(i)=cm(i)+coord0(j,i)*mass(j)
  k=k+1
  kk(k)=mass(j)
 enddo
enddo
!!! New coordinates with the cm at 0
cm=cm/totmass
do i=1,n
 do j=1,3
  coord(i,j)=coord0(i,j)-cm(j)
 enddo
enddo

!!! Reading the Hessian triangular Matrix
read (un,*) ( ( Hess(i,j) ,j=1,i ),i=1,n3 )
close(un)

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

Hess2=Hess
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
 write(6,*) i,frec(i)*2.194746E5
enddo

write(6,*) "Projecting out translation and rotation"
Hess=Hess2
!! Creating a projector with the translation and rotation
!! proj(6,3*n)
proj=0.d0
ii=0
do i=1,n
 popo=sqrt(mass(i))
 proj(1,ii+1)=popo
 proj(2,ii+2)=popo
 proj(3,ii+3)=popo

 proj(4,ii+2)=-popo*coord(i,3)
 proj(4,ii+3)= popo*coord(i,2)

 proj(5,ii+1)= popo*coord(i,3)
 proj(5,ii+3)=-popo*coord(i,1)

 proj(6,ii+1)=-popo*coord(i,2)
 proj(6,ii+2)= popo*coord(i,1)

 ii=ii+3
enddo

! Normalizing the projector
do i=1,6
 popo=0.d0
 do j=1,n3
  popo=popo+proj(i,j)**2
 enddo
! write(6,*) "Norm of projector ",i,popo
 do j=1,n3
  if (popo.ge.1e-10) proj(i,j)=proj(i,j)/sqrt(popo)
 enddo
enddo

call dmatmat(n3,6,n3,transpose(proj),proj,proj2)
proj2=-proj2
do i=1,n3
 proj2(i,i)=1.+proj2(i,i)
enddo

call dmatmat(n3,n3,n3,transpose(proj2),Hess,Hess2)
call dmatmat(n3,n3,n3,Hess2,proj2,Hess)


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
 write(6,*) i,frec(i)*2.194746E5
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
 write(1,*) frec(i)*2.194746E5,massred/1822.888,popo
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
