subroutine evaluate(printlevel,sat,nat,npot,x,V,U,W,deriv1,deriv2,dipole,writ,las,typdyn,typmodel,&
vel, T, gradV, egradV, pot)

use molecule_

implicit none
integer, intent(in) :: pot
integer, intent(in) :: nat,npot,typdyn,typmodel
complex*16 V(npot,npot),dipole(3,npot,npot)
real*8 deriv1(npot,nat,3)
complex*16 deriv2(npot,npot),U(npot,npot),T(npot,npot)
real*8 x(nat,3),vel(nat,3),W(npot)
logical writ
real*8 las(3)
character*2 sat(nat)
integer printlevel
logical egradV(npot,npot)


real*8 :: dh=0.01
complex*16 dipoletmp(3,npot,npot)
complex*16 Vm(npot,npot),Vp(npot,npot)
complex*16 Um(npot,npot),Up(npot,npot)
complex*16 gradV(npot,npot,nat,3)
complex*16 graddipole(3,npot,npot,nat,3)
complex*16 Wp(npot,npot),Wm(npot,npot)

integer i,j,k,l,q
integer ii,ij,ik



select case(typmodel)
 case(0)
  call molecule (nat, npot, x, sat, V, dipole, gradV, egradV, graddipole, dh, vel,T)
 case(1)
  call analytical (nat, npot, x, sat, V, dipole, gradV, egradV, graddipole, T, vel)
 case(2)
  T=dcmplx(0.d0,0.d0)
  call numerical (nat, npot, x, sat, V, dipole, gradV, egradV, graddipole, T, dh, vel)
 case default
  stop 'typmodel not valid in evaluation_xyz'
end select


do i=1,npot
 do j=1,npot
  Wp(i,j)=V(i,j)
 enddo
enddo

if (printlevel.ge.3) write(68,*) "Potential"
!!! Gradient is evaluated in adiabatic... so it must be reorder the potential (same than U0 in singlepoint)
!!! Here to obtain a good U matrix is neccesary, so the diagonalization is done using the previous one
call diago(npot,Wp,dipole,las,U,typdyn,printlevel,.true.,1,68,.false.)
call adiabatize(npot,U,Wp)

do i=1,npot
 W(i)=dreal(Wp(i,i))
enddo


Up=U
Um=U

if (printlevel.ge.3) then
 j=0
 k=0
 write(200,*) "Evaluating gradient"
 do ii=1,npot
  write(200,900) float(j),float(j),float(k),(Wp(ii,ij),ij=1,npot)
 enddo
 write(202,*) "Evaluating gradient"
endif
 

do q=1,nat
 do i=1,3

  do ii=1,npot
   do ij=1,npot
    Wm(ii,ij)=V(ii,ij)-dh*gradV(ii,ij,q,i)
    do ik=1,3
     dipoletmp(ik,ii,ij)=dipole(ik,ii,ij)-dh*graddipole(ik,ii,ij,q,i)
    enddo
   enddo
  enddo
  if (printlevel.ge.3) write(68,*) "Gradients", q,i
  Um=U
  call diago(npot,Wm,dipoletmp,las,Um,typdyn,printlevel,.false.,1,68,.false.)
  call adiabatize(npot,Um,Wm)



  do ii=1,npot
   do ij=1,npot
    Wp(ii,ij)=V(ii,ij)+dh*gradV(ii,ij,q,i)
    do ik=1,3
     dipoletmp(ik,ii,ij)=dipole(ik,ii,ij)+dh*graddipole(ik,ii,ij,q,i)
    enddo
   enddo
  enddo
  Up=U
  call diago(npot,Wp,dipoletmp,las,Up,typdyn,printlevel,.false.,1,68,.false.)
  call adiabatize(npot,Up,Wp)

  do j=1,npot
   deriv1(j,q,i)=( dreal(Wp(j,j))-dreal(Wm(j,j)) )/2./dh

   do k=1,npot
    if (iabs(typdyn).eq.23) deriv2(j,k)=(Up(j,k)-Um(j,k))/2./dh*vel(q,i)
   enddo
  enddo

  if (printlevel.ge.3) then
   write(202,900) (deriv1(ii,q,i),ii=1,npot)
  endif

  if (printlevel.ge.3) then
   do ii=1,npot
    do ij=1,npot
     write(200,900) float(q),float(i),float(ii),float(ij),Wp(ii,ij),Wm(ii,ij),(Wp(ii,ij)-Wm(ii,ij))/2./dh
    enddo
   enddo
  endif

 enddo
enddo




!do i=1,npot
! write(6,*) (Wp(i,j),j=1,npot)
!enddo
!stop


900 format (100(x,E20.8))

return
end subroutine



