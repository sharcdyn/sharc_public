subroutine rk4(np,c,Ham1,Ham2,Ham3,h)
implicit none
integer np
! dt
real*8, intent(in) :: h
!-i/hbarH Hamiltonians: t  / t+dt/2   /    t+dt
complex*16, intent(in) :: Ham1(np,np), Ham2(np,np), Ham3(np,np)
complex*16, intent(inout) :: c(np)
complex*16 :: K1(np),K2(np),K3(np),K4(np)
complex*16 kk(np,np)
integer i
real*8 Eref,Emin,Emax

Emin=huge(0.d0)
Emax=-Emin
do i=1,np
 Eref=-dimag(Ham1(i,i))
 if (Eref.le.Emin) Emin=Eref
 if (Eref.ge.Emax) Emax=Eref
enddo

!!! Creating a reference energy based in the Ham1
Eref=(Emin+Emax)/2.

!write(6,*) h,Eref,Emin,Emax

kk=Ham1
!$OMP PARALLEL DO
do i=1,np
 kk(i,i)=kk(i,i)+dcmplx(0.0d0,Eref)
enddo
!$OMP END PARALLEL DO
K1=c*dcmplx(h,0.d0)
call cmatvec(np,np,kk,K1)
!do i=1,np
! K1(i)=dcmplx(0.d0,0.d0)
! do j=1,np
! K1(i)=K1(i) + dcmplx(h,0.d0)*Ham1(i,j)*c(j)
! enddo
!enddo

kk=Ham2
!$OMP PARALLEL DO
do i=1,np
 kk(i,i)=kk(i,i)+dcmplx(0.d0,Eref)
enddo
!$OMP END PARALLEL DO

K2=dcmplx(h,0.d0)*(c+K1/2.)
call cmatvec(np,np,kk,K2)
!do i=1,np
! K2(i)=dcmplx(0.d0,0.d0)
! do j=1,np
!  K2(i)=K2(i) + dcmplx(h,0.d0)*Ham2(i,j)*(c(j)+K1(j)/2.)
! enddo
!enddo

K3=dcmplx(h,0.d0)*(c+K2/2.)
call cmatvec(np,np,kk,K3)
!do i=1,np
! K3(i)=dcmplx(0.d0,0.d0)
! do j=1,np
!  K3(i)=K3(i) + dcmplx(h,0.d0)*Ham2(i,j)*(c(j)+K2(j)/2.)
! enddo
!enddo

kk=Ham3
!$OMP PARALLEL DO
do i=1,np
 kk(i,i)=kk(i,i)+dcmplx(0.0d0,Eref)
enddo
!$OMP END PARALLEL DO
K4=dcmplx(h,0.d0)*(c+K3)
call cmatvec(np,np,kk,K4)
!do i=1,np
! K4(i)=dcmplx(0.d0,0.d0)
! do j=1,np
!  K4(i)=K4(i) + dcmplx(h,0.d0)*Ham3(i,j)*(c(j)+K3(j))
! enddo
!enddo

!! New coefficients including the extra phase
c=c+(K1+2.*K2+2.*K3+K4)/6.
c=c*cdexp(dcmplx(0.d0,-Eref*h))
!do i=1,np
! c(i)=c(i)+(1./6.)*(K1(i)+2.*K2(i)+2.*K3(i)+K4(i))
!enddo

end
		
		
