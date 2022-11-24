subroutine geomwrit(sat,nat,time,x,vel,mass,npot,E,pot)
implicit none
integer, intent(in) :: npot,pot,nat
real*8, intent(in) :: E(npot),time
real*8, intent(in) :: x(nat,3),vel(nat,3),mass(nat)
character*2, intent(in) :: sat(nat)

integer i,j

real*8 Ekt

Ekt=0.
do i=1,nat
 do j=1,3
  Ekt=Ekt+0.5*mass(i)*vel(i,j)**2
 enddo
enddo

write(11,*) nat
write(11,900) Ekt," TIME ",time/41.32,pot,(E(i),i=1,npot)
do j=1,nat
 write(11,901) sat(j),(x(j,i)*.5292,i=1,3),(vel(j,i)*1000.,i=1,3)
enddo

900 format(E20.10E3,A6,E20.10E3,x,I3.3,1000(x,E20.10e3))
901 format(A2,8(x,E20.10e3))
return
end subroutine

subroutine coeffwrit(time,npot,U,coef,pot,Uprop,deco)
implicit none
integer, intent(in) :: npot,pot
complex*16, intent(in) :: coef(npot),U(npot,npot),Uprop(npot,npot)
real*8, intent(in) :: time,deco

real*8 suma

integer i
complex*16 coefd(npot)

coefd=coef
call diabat_coef (npot,U,coefd)

suma=0.d0
do i=1,npot
 suma=suma+cdabs(coef(i))**2
enddo

! open(12,file="coef.out",   access="append")
! open(14,file="pop_a.out",  access="append")


write(12,900) time/41.32,suma,(dreal(coef(i)),dimag(coef(i)),i=1,npot)
if (isnan(suma)) stop "The sum of the coefficients is Nan"

if (deco.lt.0) suma=1.0
write(14,900) time/41.32,(cdabs(coef(i))**2/suma,i=1,npot)


! open(67,file="coef_d.out", access="append")
! open(13,file="pop_d.out",  access="append")

suma=0.d0
do i=1,npot
 suma=suma+cdabs(coefd(i))**2
enddo
write(67,900) time/41.32,suma,(dreal(coefd(i)),dimag(coefd(i)),i=1,npot)
if (deco.lt.0) suma=1.0
write(13,900) time/41.32,(cdabs(coefd(i))**2/suma,i=1,npot)

! open(65,file="coef_raw.out",access="append")
! open(64,file="pop_raw.out", access="append")

!! This rotates from diabatic representation to raw data
call diabat_coef (npot,Uprop,coefd)
suma=0.d0
do i=1,npot
 suma=suma+cdabs(coefd(i))**2
enddo
write(65,900) time/41.32,suma,(dreal(coefd(i)),dimag(coefd(i)),i=1,npot)
if (deco.lt.0) suma=1.0
write(64,900) time/41.32,(cdabs(coefd(i))**2/suma,i=1,npot)


! open(16,file="pop_d2.out", access="append")

do i=1,npot
 coefd(i)=dcmplx(0.d0,0.d0)
enddo
coefd(pot)=dcmplx(1.d0,0.d0)
call diabat_coef(npot,U,coefd)
suma=0.d0
do i=1,npot
 suma=suma+cdabs(coefd(i))**2
enddo
if (deco.lt.0) suma=1.0
write(16,900) time/41.32,(cdabs(coefd(i))**2/suma,i=1,npot)



900 format(2000(x,E20.10e3))
return
end





subroutine properties(nat,time,npot,pot,Ek,Ep,dipole, Epot, deriv2,U)
implicit none
integer, intent(in) :: npot,pot, nat
real*8, intent(in) :: time
real*8, intent(in) :: Ek(nat,3),Ep(npot), Epot
complex*16, intent(in) :: dipole(3,npot,npot),deriv2(npot,npot)
complex*16, intent(in) :: U(npot,npot)

integer i,j,k
real*8 Ekt
complex*16 kk(npot,npot)

Ekt=0.
do i=1,nat
 do j=1,3
  Ekt=Ekt+Ek(i,j)
 enddo
enddo

write(15,900) time/41.32,Ekt,Epot,(Ep(i),i=1,npot)
do i=1,3
 kk=dipole(i,:,:)
 call adiabatize(npot,U,kk)
 write(30+i,900)  time/41.32,((kk(j,k),k=1,npot),j=1,npot)
enddo
write(40,900) time/41.32,((deriv2(j,k),k=1,npot),j=1,npot)

900 format(200000(x,E20.10e3))
end

subroutine stepwrit(npot,pot,npoint,U)
implicit none
integer, intent(in) :: pot,npot,npoint
complex*16, intent(in) :: U(npot,npot)

real*8 popd(npot)
integer st(npot)

integer i,j,k
real*8 kk

!!! Getting diabatic (input) population
do i=1,npot
 popd(i)=cdabs(U(i,pot))**2
 st(i)=i
enddo

!!! Reordering the population to know what is the most important diabatic contribution
do i=1,npot-1
 do j=i+1,npot
  if (popd(j).ge.popd(i)) then
   kk=popd(i)
   popd(i)=popd(j)
   popd(j)=kk
   k=st(i)
   st(i)=st(j)
   st(j)=k
  endif
 enddo
enddo

open(1,file="npoint.dat")
write(1,"(I9.9,x,I9.9)") npoint,pot
do i=1,npot
 write(1,"(I9.9,x,E30.20e3)") st(i),popd(i)
enddo
close(1)

end


subroutine gradwrit(nat,deriv,un)
implicit none
integer, intent(in) :: nat,un
real*8, intent(in) :: deriv(nat,3)

integer i

do i=1,nat
 write(un,"(3(x,E20.10e3))") deriv(i,:)
enddo

return
end

subroutine gradwritc(nat,deriv,un)
implicit none
integer, intent(in) :: nat,un
complex*16, intent(in) :: deriv(nat,3)

integer i,j

do i=1,nat
 write(un,"(6(x,E20.10e3))") (deriv(i,j),j=1,3)
enddo

return
end
