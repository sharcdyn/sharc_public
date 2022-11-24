!Minimizing UºTU+UºdU/dt by changing signs and permutations of U matrix
subroutine permuting (printlevel,npot,deriv2,U,U_prev,T,dt, pot, coef, time, change)
implicit none
integer npot, pot
complex*16 U(npot,npot),U_prev(npot,npot)
complex*16 coef(npot), deriv2(npot,npot)
complex*16 T(npot, npot)
real*8 dt, time, kk2

complex*16 kkm(npot,npot),kkp(npot,npot),kk(npot,npot),kkperp(npot,npot),kkperm(npot,npot)
integer i,j,k,l, a, b, c, n1, kkpot, kk2pot
complex*16 Um(npot,npot),Up(npot,npot), Uper(npot,npot), Uchachi(npot,npot)
complex*16 Uperp(npot,npot), Uperm(npot,npot), Tp(npot,npot),Tm(npot,npot),Tperp(npot,npot)
complex*16 Tperm(npot,npot)
real*8 sqm,sqp,sqperp,sqperm
complex*16 kkkk
logical permut, change 
integer printlevel

if (printlevel.ge.3) write(99,*) "Potential:", pot
kk2=0.d0
kk2pot=pot
Uchachi=U_prev

call adiabatize(npot,U,Tp)
do b=1,npot
 do c=b+1,npot
  kkpot=pot
  Uper=Uchachi
  do a=1,npot
   Uper(b,a)=Uchachi(c,a)
   Uper(c,a)=Uchachi(b,a)
  enddo
  Tp=T
!  Tm=T
!  Tperp=T
!  Tperm=T
!   !!! Minimizing the derivative by changing the sign at every state
!   !!! Here, we change the sign of "l" column for both: normal U and permuted U
  Up=Uchachi
!  call adiabatize(npot,Up,Tp)
  Uperp=Uper
!  call adiabatize(npot,Uperp,Tperp)
  do i=1,npot
   do j=1,npot
    kkp(i,j)=dcmplx(0.d0,0.d0) + Tp(i,j)
    kkperp(i,j)=dcmplx(0.d0,0.d0) + Tp(i,j)
    do k=1,npot
     kkp(i,j)=kkp(i,j)+dconjg(U(k,i))*(U(k,j)-Up(k,j))/dt
     kkperp(i,j)=kkperp(i,j)+dconjg(U(k,i))*(U(k,j)-Uperp(k,j))/dt
    enddo
   enddo
  enddo
  sqp=0.
  sqperp=0.
  do i=1,npot
   do j=1,npot
    sqp=sqp+cdabs(kkp(i,j))**2
    sqperp=sqperp+cdabs(kkperp(i,j))**2
   enddo
  enddo

  do l=1,npot
   Um=Uchachi
  ! call adiabatize (npot,Um, Tm)
   Uperm=Uper
  ! call adiabatize (npot,Uperm, Tperm)
   do j=1,npot
    Um(j,l)=dcmplx(-1.d0,0.d0)*Uchachi(j,l)
    Uperm(j,l)=dcmplx(-1.d0,0.d0)*Uper(j,l)
   enddo
!   !!! Here, once we have the 4 matrices, we test what is the | |^2 value of Uº(dU/dt)
   do i=1,npot
    do j=1,npot
     kkm(i,j)=dcmplx(0.d0,0.d0) + Tp(i,j)
     kkperm(i,j)=dcmplx(0.d0,0.d0) + Tp(i,j)
     do k=1,npot
      kkm(i,j)=kkm(i,j)+dconjg(U(k,i))*(U(k,j)-Um(k,j))/dt
      kkperm(i,j)=kkperm(i,j)+dconjg(U(k,i))*(U(k,j)-Uperm(k,j))/dt
     enddo
    enddo
   enddo

!  enddo
   sqm=0.
   sqperm=0.
   do i=1,npot
    do j=1,npot
     sqm=sqm+cdabs(kkm(i,j))**2
     sqperm=sqperm+cdabs(kkperm(i,j))**2
    enddo
   enddo
   if (printlevel.ge.3) then
    write(99,"(F15.10,x,A21,2(x,I2.2),4F15.4)") time/41.32 ,&
    "Checking permutation " ,b,c, sqp, sqm, sqperp, sqperm
   endif
!   !!! We test the minimum U matrix of our 4
   if(sqp.lt.sqm.and.sqp.lt.sqperp.and.sqp.lt.sqperm) then
    Uchachi=Up
    permut=.false.
   endif
   if(sqm.lt.sqp.and.sqm.lt.sqperp.and.sqm.lt.sqperm) then
    Uchachi=Um
    permut=.false.
   endif
   if(sqperp.lt.sqp.and.sqperp.lt.sqm.and.sqperp.lt.sqperm) then
    Uchachi=Uperp
    permut=.true.
   endif
   if(sqperm.lt.sqperp.and.sqperm.lt.sqp.and.sqperm.lt.sqm) then
    Uchachi=Uperm
    permut=.true.
   endif
   if(permut) then
    if (printlevel.ge.3) write(99,*) time/41.32 ,"Ha permutao",b, c, "Signo", l
    U_prev=Uchachi
    if (change) then
     if(printlevel.ge.3) write(99,*) time/41.32 ,"Coeficientes",b, c, "Signo", l
     kkkk=dcmplx(1.d0,0.d0)*coef(b)
     coef(b)=coef(c)
     coef(c)=kkkk
     kkpot=pot
     if (b.eq.pot) kkpot=c
     if (c.eq.pot) kkpot=b
     pot=kkpot
    endif
   else
    if (printlevel.ge.3) write(99,*) time/41.32, "No permutao", b, c, "Signo", l
   endif
  enddo
  if ((permut).and.change.eqv..false.) write(99,*) time/41.32, "Permutao", b,c
  if ((permut).and.change.eqv..true.) write(99,*) time/41.32, "Coefis", b,c
 enddo
enddo

!   !!! Testing for a diabatic jump
if (kk2pot.ne.pot) then
 write(66,901) time/41.32, kk2pot,pot, kk2,kk2, kk2, kk2,kk2, kk2,"Diabatic","Norma=", kk2
endif

if (printlevel.ge.3) write(99,*)"Salida:" ,pot

900 format(100(x,E20.10e3))
901 format(F20.10,x,I3.3,x,I3.3,x,6(E20.10e3,x),2(A10,x),E20.10E3)


return
end subroutine



!Minimizing UºdU/dt by changing signs of U matrix



subroutine minimizeT(npot,T, T_prev)
implicit none
integer npot
complex*16 U(npot,npot),T(npot,npot), T_prev(npot,npot)
complex*16 Um(npot,npot),Tm(npot,npot)
complex*16 Up(npot,npot),Tp(npot,npot)
real*8 kkm,kkp,dt

integer i,j,k,l

do i=1,npot
! Changing sign in state i
 do j=1,npot
  do k=1,npot
   Tp(j,k)=T(j,k)
   Tm(j,k)=T(j,k)
   Up(j,k)=dcmplx(0.d0,0.d0)
   Um(j,k)=dcmplx(0.d0,0.d0)
  enddo
  Up(j,j)=dcmplx(1.d0,0.d0)
  Um(j,j)=dcmplx(1.d0,0.d0)
 enddo
 Um(i,i)=dcmplx(-1.d0,0.d0)
 call adiabatize(npot,Up,Tp)
 call adiabatize(npot,Um,Tm)
 kkp=0.d0
 kkm=0.d0
 do j=1,npot
  do k=1,npot
   kkp=kkp+cdabs(Tp(j,k)-T_prev(j,k))**2
   kkm=kkm+cdabs(Tm(j,k)-T_prev(j,k))**2
  enddo
 enddo

 if (kkm.le.kkp) then
  T=Tm
 endif
enddo

return
end subroutine

subroutine adiabatize(npot,U,A)
implicit none
integer, intent(in) :: npot
complex*16,intent(in) :: U(npot,npot)
complex*16, intent(inout) :: A(npot,npot)
complex*16 kk(npot,npot),U2(npot,npot)

U2=dconjg(transpose(U))
call cmatmat(npot,npot,npot,A,U,kk)
call cmatmat(npot,npot,npot,U2,kk,A)

return
end subroutine

subroutine diabatize(npot,U,A)
implicit none
integer, intent(in) :: npot
complex*16, intent(in) :: U(npot,npot)
complex*16, intent(inout) :: A(npot,npot)
complex*16 kk(npot,npot),U2(npot,npot)

U2=dconjg(transpose(U))
call cmatmat(npot,npot,npot,A,U2,kk)
call cmatmat(npot,npot,npot,U,kk,A)

return
end subroutine

subroutine adiabat_coef (npot,U,coef)
implicit none
integer, intent(in) :: npot
complex*16, intent(in) :: U(npot,npot)
complex*16, intent(inout) :: coef(npot)
complex*16 A(npot)
integer i,j

do i=1,npot
 A(i)=dcmplx(0.d0,0.d0)
 do j=1,npot
  A(i)=A(i)+coef(j)*U(j,i)
 enddo
enddo

do i=1,npot
 coef(i)=A(i)
enddo

return
end subroutine


subroutine diabat_coef (npot,U,coef)
implicit none
integer, intent(in) :: npot
complex*16, intent(in) :: U(npot,npot)
complex*16, intent(inout) :: coef(npot)
complex*16 A(npot)
integer i,j

do i=1,npot
 A(i)=dcmplx(0.d0,0.d0)
 do j=1,npot
  A(i)=A(i)+coef(j)*dconjg(U(i,j))
 enddo
enddo

do i=1,npot
 coef(i)=A(i)
enddo

end subroutine

subroutine minimizeVdipT(npot,V,dip,T, Umini, time, npol,dt,Uprop)
implicit none
integer npot, npol
complex*16 V(npol,npot,npot), dip(npol,3,npot,npot)
complex*16 Vp(npot,npot),Vm(npot,npot),T(npol,npot,npot)
complex*16 dipp(3,npot,npot), dipm(3,npot,npot),Tm(npot,npot), Tp(npot,npot)
real*8 Umini(npot,npot), Utmp(npot,npot), U2(npot,npot)

integer i,j,k,l,n, a,b
complex*16 Up(npot,npot),Um(npot,npot),Uprop(npot,npot)
real*8 kkp,kkm, time, over(npot,npot),dt
real*8 kk, maxim


do i=1,npot
 do j=1,npot
  Umini(i,j)=0.0d0
 enddo
 Umini(i,i)=1.d0
enddo




do n=1,npot
 do i=1,npot
! Changing sign in state i
  do j=1,npot
   do k=1,npot
    do l=1,3
     dipp(l,j,k)=dip(1,l,j,k)
     dipm(l,j,k)=dip(1,l,j,k)
    enddo
    Vp(j,k)=V(1,j,k)
    Vm(j,k)=V(1,j,k)
    Tp(j,k)=T(1,j,k)
    Tm(j,k)=T(1,j,k)
    Up(j,k)=dcmplx(Umini(i,j),0.d0)
    Um(j,k)=dcmplx(Umini(i,j),0.d0)
   enddo
  enddo
  Up(i,i)=dcmplx(+Umini(i,i),0.d0)
  Um(i,i)=dcmplx(-Umini(i,i),0.d0)
  call adiabatize(npot,Up,Tp)
  call adiabatize(npot,Um,Tm)
  call adiabatize(npot,Up,Vp)
  call adiabatize(npot,Um,Vm)
  do l=1,3
   call adiabatize(npot,Up,dipp(l,:,:))
   call adiabatize(npot,Um,dipm(l,:,:))
  enddo
  kkp=0.d0
  kkm=0.d0
  do j=1,npot
   do k=1,npot
    kkp=kkp+cdabs(Vp(j,k)-V(2,j,k))**2+0.5*cdabs(Tp(j,k)-T(2,j,k))**2
    kkm=kkm+cdabs(Vm(j,k)-V(2,j,k))**2+0.5*cdabs(Tm(j,k)-T(2,j,k))**2
    do l=1,3
     kkp=kkp+0.1*cdabs(dipp(l,j,k)-dip(2,l,j,k))**2
     kkm=kkm+0.1*cdabs(dipm(l,j,k)-dip(2,l,j,k))**2
    enddo
   enddo
  enddo

  if (kkm.lt.kkp) then
   write(68,"(a36,2(x,I3.3),3(x,F20.10),x(F10.2))") "Minimize has changed the sign of ", i,n,kkp,kkm,kkp-kkm,time/41.32
   Umini(i,i)=-Umini(i,i)
  endif
 enddo
enddo

do i=1,npot
 do j=1,npot
  Up(i,j)=dcmplx(Umini(i,j),0.d0)
 enddo
enddo

call adiabatize(npot,Up,T(1,:,:))
call adiabatize(npot,Up,V(1,:,:))
do i=1,3
 call adiabatize(npot,Up,dip(1,i,:,:))
enddo
Uprop=Up

return
end subroutine

subroutine createT(printlevel,npot, V,dip, T, Umini,time,dt,Uprop)
implicit none
integer, intent(in) :: printlevel
integer, intent(in) :: npot
real*8, intent(in) :: dt,time

real*8, intent(inout) :: Umini(npot,npot)
complex*16, intent(out) :: T(npot,npot)
complex*16, intent(inout) :: V(npot,npot),dip(3,npot,npot)
complex*16, intent(out) :: Uprop(npot,npot)

real*8 Utmp(npot,npot), Uprev(npot,npot)

integer i,j,k,l,n, a,b
complex*16 Up(npot,npot),Um(npot,npot),over2(npot,npot)
real*8 kkp,kkm, over(npot,npot)
real*8 kk, maxim
complex*16 tmp(npot,npot)

logical get
real*8 over_tmp(npot,npot)
logical checked(npot,npot)
integer maxv(2)

logical ex

Uprev=Umini
Utmp=0.d0

open(1,file="overlap.dat")
do i=1,npot
 ! Reading overlap < act | prev >
 read (1,*) (over(i,j),j=1,npot)
 !!! if you get < prev| act >  MOLCAS?
! read (1,*) (over(j,i),j=1,npot)
enddo
close(1)

over_tmp=transpose(over)
! Defining the change of Umini with respect to the previous one (Utmp) by getting maxima
checked=.true.
do k=1,npot
 maxv=maxloc(over_tmp**2,checked)
 i=maxv(1)
 j=maxv(2)
 checked(i,:)=.false.
 checked(:,j)=.false.
 Utmp(i,j)=1.d0
 if (maxim.lt.0.5) write(68,"(a,2(I5.4,x),2(x,F20.10))") "createT warning small biggest overlap ",j,i,over_tmp(i,j),time/41.32
 if (over_tmp(i,j).lt.0) then
  Utmp(i,j)=-Utmp(i,j)
  write(68,"(A,2(x,I2),x,F7.4,x,F20.10,x,F10.2)") "createT has changed the sign of ", j,i, Umini(i,j), over_tmp(i,j),time/41.32
 endif
 if (i.ne.j) write(68,"(A,2(x,I5),2(x,F20.10))") "createT has interchanged the states:",j,i, over_tmp(i,j), time/41.32
enddo

Utmp=transpose(Utmp)
call dmatmat(npot,npot,npot,Utmp,Uprev,Umini)
Uprop=dcmplx(Umini,0.d0)
Um=dcmplx(Uprev,0.d0)
Up=dcmplx(Utmp,0.d0)


! Rotating V and dip...
call adiabatize(npot,Uprop,V)
do l=1,3
 call adiabatize(npot,Uprop,dip(l,:,:))
enddo

!Calculating T
over2=dcmplx(over,0.d0)
call cmatmat(npot,npot,npot,transpose(dconjg(Up)),over2,tmp)
over2=tmp
call adiabatize(npot,Um,over2)
do i=1,npot
 do j=1,npot
  if (i.eq.j) then
   T(i,j)=dcmplx(1.d0,0.d0)
  else
   T(i,j)=dcmplx(0.d0,0.d0)
  endif
  T(i,j)=(T(i,j)-over2(i,j))/dt
 enddo
enddo

open(1,file="overlap.new")
write(1,*) "Uprev"
do i=1,npot
 write(1,901) (Uprev(i,j),j=1,npot)
enddo
write(1,"(A,10000(x,F6.2))") "overlap ",(over(i,i),i=1,npot)
do i=1,npot
 write(1,900) ((over(i,j)),j=1,npot)
enddo
write(1,*) "Uchange"
do i=1,npot
 write(1,901) (Utmp(i,j),j=1,npot)
enddo
write(1,*) "Unew"
do i=1,npot
 write(1,901) (Umini(i,j),j=1,npot)
enddo
write(1,"(A,10000(x,F6.2))") "overlap with previous ",(dreal(tmp(i,i)),i=1,npot)
do i=1,npot
 write(1,900) (dreal(over2(i,j)),j=1,npot)
enddo
write(1,"(A,10000(x,F6.2))") "New overlap with start ",(dreal(over2(i,i)),i=1,npot)
write(1,*) "created T"
do i=1,npot
 write(1,900) (dreal(T(i,j)),j=1,npot)
enddo
write(1,*) "created T imag (should be zero)"
do i=1,npot
 write(1,900) (dimag(T(i,j)),j=1,npot)
enddo
close(1)

if (printlevel.ge.3) then
 write(104,900) time/41.32,((cdabs(T(i,j)),j=1,npot),i=1,npot)
 write(105,900) time/41.32,((Umini(i,j),j=1,npot),i=1,npot)
endif

!!! TESTING Utmp is unitary
call dmatmat(npot,npot,npot,transpose(Utmp),Utmp,Uprev)
!Uprev=matmul(transpose(Utmp),Utmp)
do j=1,npot
 do i=1,npot
!  Uprev(i,j)=0
!  do k=1,npot
!   Uprev(i,j)=Uprev(i,j)+Utmp(k,i)*Utmp(k,j)
!  enddo
  if (i.eq.j.and.Uprev(i,j).ne.1) then
   write(6,*) "Problem with Utmp ",i,Uprev(i,j)
   do k=1,npot
    write(6,*) Utmp(k,:)
   enddo
   stop
  endif
  if (i.ne.j.and.Uprev(i,j).ne.0) then
   write(6,*) "Problem with Utmp ",i,j,Uprev(i,j)
   do k=1,npot
    write(6,*) Utmp(k,:)
   enddo
   stop
  endif
 enddo
enddo



900 format(1000(x,e20.10e3))
901 format(1000(x,f4.1))

end subroutine

subroutine createT_diab(printlevel,npot, V,dip, T, Umini,time,dt,Uprop,Vprev)
use diab_
implicit none
integer, intent(in) :: printlevel
integer, intent(in) :: npot
real*8, intent(in) :: dt,time

complex*16, intent(in) :: Vprev(npot,npot)
real*8, intent(inout) :: Umini(npot,npot)
complex*16, intent(out) :: T(npot,npot)
complex*16, intent(inout) :: V(npot,npot),dip(3,npot,npot)
complex*16, intent(out) :: Uprop(npot,npot)

real*8 Utmp(npot,npot), Uprev(npot,npot)

real*8 over(npot,npot),over_tmp(npot,npot)
complex*16 Up(npot,npot),Um(npot,npot),Uk(npot,npot)
complex*16 over2(npot,npot),tmp(npot,npot),over3(npot,npot)
integer diab_order(npot)

integer i,j,k,l

Uprev=Umini
Utmp=0.d0
do i=1,npot
 Utmp(i,i)=1.0d0
enddo


open(1,file="overlap.dat")
do i=1,npot
 ! Reading overlap < act | prev >
 read (1,*) (over(i,j),j=1,npot)
 !!! if you get < prev| act >  MOLCAS?
! read (1,*) (over(j,i),j=1,npot)
enddo
close(1)

do i=1,npot
 diab_order(i)=i
enddo
over_tmp=transpose(over)
!! Diabatizing over for each potential with respect to the previous one
call diab_projord(npot,npot,over_tmp,diab_order,Utmp)
call dmatmat(npot,npot,npot,Utmp,Uprev,Umini)

Up=dcmplx(Umini,0.d0)
Um=dcmplx(Uprev,0.d0)
Uk=dcmplx(Utmp,0.d0)


! Rotating V and dip...
call adiabatize(npot,Up,V)
do l=1,3
 call adiabatize(npot,Up,dip(l,:,:))
enddo
!! When diabatizing with Up the raw representation is obtained
!Uprop=dconjg(transpose(Up))
Uprop=Up


!Calculating T
over2=dcmplx(over,0.d0)
call cmatmat(npot,npot,npot,transpose(dconjg(Uk)),over2,tmp)
over3=tmp
over2=tmp
call adiabatize(npot,Um,over2)


do i=1,npot
 do j=1,npot
  if (i.eq.j) then
   T(i,j)=dcmplx(1.d0,0.d0)
  else
   T(i,j)=dcmplx(0.d0,0.d0)
  endif
  T(i,j)=(T(i,j)-over2(i,j))/dt
 enddo
enddo

open(1,file="overlap.new")
write(1,*) "Uprev"
do i=1,npot
 write(1,900) (Uprev(i,j),j=1,npot)
enddo
write(1,"(A,1000(x,f5.2))") "overlap ",(over(i,i),i=1,npot)
do i=1,npot
 write(1,900) ((over(i,j)),j=1,npot)
enddo
write(1,*) "U change"
do i=1,npot
 write(1,900) (Utmp(i,j),j=1,npot)
enddo
write(1,"(A,1000(x,f5.2))") "New overlap with previous",(dreal(over3(i,i)),i=1,npot)
do i=1,npot
 write(1,900) (dreal(over3(i,j)),j=1,npot)
enddo
write(1,*) "U new"
do i=1,npot
 write(1,900) (Umini(i,j),j=1,npot)
enddo
write(1,"(A,1000(x,f5.2))") "New overlap with start",(dreal(over2(i,i)),i=1,npot)
do i=1,npot
 write(1,900) (dreal(over2(i,j)),j=1,npot)
enddo
write(1,*) "created T"
do i=1,npot
 write(1,900) (dreal(T(i,j)),j=1,npot)
enddo
write(1,*) "created T imag (should be zero)"
do i=1,npot
 write(1,900) (dimag(T(i,j)),j=1,npot)
enddo
close(1)

do i=1,npot
 if (dreal(over2(i,i)).lt.0.9) then
  write(68,"(A,I0,3(X,F5.2))") "Diabatization of the overlap matrix can only be done up to ",i,over2(i,i),time/41.32
 endif
enddo

!call use_old_phase(npot,Vprev,V,dip,T,Um,Up)

if (printlevel.ge.3) then
 write(104,900) time/41.32,((cdabs(T(i,j)),j=1,npot),i=1,npot)
 write(105,900) time/41.32,((Umini(i,j),j=1,npot),i=1,npot)
endif



900 format(1000(x,e20.10e3))
901 format(1000(x,i2))

end subroutine

!!! Diabatizationg based in the SO-overlap
subroutine createT_diab2(printlevel,npot,npol, Vt,dipt, Tt, Umini,time,dt,Uprop)
use diab_
implicit none
integer, intent(in) :: printlevel
integer, intent(in) :: npot,npol
real*8, intent(in) :: dt,time

real*8, intent(inout) :: Umini(npot,npot)
complex*16, intent(inout) :: Tt(npol,npot,npot)
complex*16, intent(inout) :: Vt(npol,npot,npot),dipt(npol,3,npot,npot)
complex*16, intent(inout) :: Uprop(npot,npot)

real*8 Utmp(npot,npot), Uprev(npot,npot)

real*8 over(npot,npot),over_tmp(npot,npot)
complex*16 Up(npot,npot),Um(npot,npot),Uk(npot,npot)
complex*16 over2(npot,npot),tmp(npot,npot),over3(npot,npot)
integer diab_order(npot)
complex*16 V(npot,npot),dip(3,npot,npot),T(npot,npot)
real*8 W(npot)

integer i,j,k,l

V=Vt(1,:,:)
dip=dipt(1,:,:,:)
T=Tt(1,:,:)

! Previous rotation matrix
Um=Uprop

open(1,file="overlap.dat")
do i=1,npot
 ! Reading overlap < act | prev >
 read (1,*) (over(i,j),j=1,npot)
 !!! if you get < prev| act >  MOLCAS?
! read (1,*) (over(j,i),j=1,npot)
enddo
close(1)

do i=1,npot
 diab_order(i)=i
enddo

over2=dcmplx(transpose(over),0.d0)
! SO-Overlap with respect to the previous function
call cmatmat(npot,npot,npot,over2,Uk,over3)


!! Diabatizing over for each potential with respect to the previous one
tmp=Um
call cdiab_projord(npot,npot,over3,diab_order,tmp)
call cmatmat(npot,npot,npot,tmp,Um,Up)
Uprop=Up

! Rotating V and dip...
call adiabatize(npot,Up,V)
do l=1,3
 call adiabatize(npot,Up,dip(l,:,:))
enddo
!! When diabatizing with Up the raw representation is obtained
! Uprop=dconjg(transpose(Up))


!Calculating T
over2=dcmplx(over,0.d0)
call cmatmat(npot,npot,npot,dconjg(transpose(Up)),over2,over3)
call cmatmat(npot,npot,npot,over3,Um,over2)

do i=1,npot
 do j=1,npot
  if (i.eq.j) then
   T(i,j)=dcmplx(1.d0,0.d0)
  else
   T(i,j)=dcmplx(0.d0,0.d0)
  endif
  T(i,j)=(T(i,j)-over2(i,j))/dt
 enddo
enddo

open(1,file="overlap.new")
write(1,*) "Uprev real"
do i=1,npot
 write(1,900) (dreal(Um(i,j)),j=1,npot)
enddo
write(1,*) "Uprev imag"
do i=1,npot
 write(1,900) (dimag(Um(i,j)),j=1,npot)
enddo
write(1,"(A,1000(x,f5.2))") "overlap ",(over(i,i),i=1,npot)
do i=1,npot
 write(1,900) (over(i,j),j=1,npot)
enddo
write(1,*) "U new real"
do i=1,npot
 write(1,900) (dreal(Up(i,j)),j=1,npot)
enddo
write(1,*) "U new imag"
do i=1,npot
 write(1,900) (dimag(Up(i,j)),j=1,npot)
enddo
write(1,"(A,1000(x,f5.2))") "New overlap with start",(dreal(over2(i,i)),i=1,npot)
write(1,*) "real"
do i=1,npot
 write(1,900) (dreal(over2(i,j)),j=1,npot)
enddo
write(1,*) "imag"
do i=1,npot
 write(1,900) (dimag(over2(i,j)),j=1,npot)
enddo
write(1,*) "created T"
do i=1,npot
 write(1,900) (dreal(T(i,j)),j=1,npot)
enddo
write(1,*) "created T imag (should be zero)"
do i=1,npot
 write(1,900) (dimag(T(i,j)),j=1,npot)
enddo
close(1)

do i=1,npot
 if (cdabs(over2(i,i))**2.lt.0.9) then
  write(68,"(A,I0,3(X,F5.2))") "Diabatization of the overlap matrix can only be done up to ",i,over2(i,i),time/41.32
 endif
enddo

if (printlevel.ge.3) then
 write(104,900) time/41.32,((cdabs(T(i,j)),j=1,npot),i=1,npot)
 write(105,900) time/41.32,((Umini(i,j),j=1,npot),i=1,npot)
endif

Vt(1,:,:)=V
Tt(1,:,:)=T
dipt(1,:,:,:)=dip


900 format(1000(x,e20.10e3))
901 format(1000(x,i2))

end subroutine


subroutine use_old_phase(npot,Vprev,V,dipole,T,Um,Up)

implicit none
integer, intent(in) :: npot
complex*16, intent(in) :: Vprev(npot,npot)
complex*16, intent(in) :: Um(npot,npot),Up(npot,npot)
complex*16, intent(inout) :: V(npot,npot),dipole(3,npot,npot),T(npot,npot)

complex*16 SOCp(npot,npot)
complex*16 SOC(npot,npot)
complex*16 phase

integer i,j,k


!! Going to the input representation
SOCp=Vprev
SOC=V
call diabatize(npot,Up,SOC)
call diabatize(npot,Um,SOCp)

!!! Removing diagonal terms to evaluate only the phase of the SOC
!$OMP PARALLEL DO
do i=1,npot
 SOC(i,i)=dcmplx(0.d0,0.d0)
 SOCp(i,i)=dcmplx(0.d0,0.d0)
enddo
!$OMP END PARALLEL DO

!!! Going back to the diabatic representation
call adiabatize(npot,Up,SOC)
call adiabatize(npot,Um,SOCp)

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(phase,k)
do j=1,npot
 do i=1,npot
  phase=dcmplx(1.d0,0.d0)
  if (cdabs(SOC(i,j)).ge.1e-8 ) then
   phase=phase*cdexp(dcmplx(0.d0, atan2(dimag(SOCp(i,j)),dreal(SOCp(i,j)) ) ))
   phase=phase*cdexp(dcmplx(0.d0,-atan2(dimag(SOC(i,j)), dreal(SOC(i,j))  ) ))
  endif
  do k=1,5
   if (k.eq.4) then 
    V(i,j)=V(i,j)*phase
   else if (k.eq.5) then
    T(i,j)=T(i,j)*phase
   else
    dipole(k,i,j)=dipole(k,i,j)*phase
   endif
  enddo

 enddo
enddo
!$OMP END PARALLEL DO

return
end subroutine

subroutine diab_all(fich,npot,V,dip,T)
!!! Subroutine to diabatize the potential and dipole and reset the T to zero
!!! It reads the <act|ref> overlap, creates a rotation matrix to minimize the off-diagonal elements
!!! and rotates everything with it
use diab_
implicit none
integer, intent(in) :: npot
character*100, intent(in) :: fich
complex*16, intent(inout) :: V(npot,npot),dip(3,npot,npot),T(npot,npot)

real*8 overlap(npot,npot)
real*8 U(npot,npot)
complex*16 U2(npot,npot)
integer diab_order(npot)
integer i,j

U=0.d0
do i=1,npot
 diab_order(i)=i
 U(i,i)=1.d0
enddo

open(1,file=fich)
do i=1,npot
!! Reading as transpose <act|ref>
 read(1,*) (overlap(j,i),j=1,npot)
enddo
close(1)

call diab_projord(npot,npot,overlap,diab_order,U)
U2=dcmplx(U,0.d0)

! Rotating
call adiabatize(npot,U2,V)
do i=1,3
 call adiabatize(npot,U2,dip(i,:,:))
enddo
T=dcmplx(0.d0,0.d0)

end subroutine diab_all
