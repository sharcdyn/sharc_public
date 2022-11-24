subroutine singlepoint(printlevel,nat,npot,x,sat,V,dipole,Tpart,vel,U,las,typdyn)
implicit none
integer, intent(in) :: npot,nat,printlevel,typdyn
real*8, intent(in) :: x(nat,3),vel(nat,3),las(3)
character*2, intent(in) :: sat(nat)

complex*16, intent(out) :: V(npot,npot),dipole(3,npot,npot),Tpart(npot,npot)
complex*16, intent(out) :: U(npot,npot)

integer i
character*100 char_tmp
logical ex
complex*16 Vtmp(npot,npot)

V=dcmplx(0.d0,0.d0)
dipole=dcmplx(0.d0,0.d0)
Tpart=dcmplx(0.d0,0.d0)

call senddata(x,nat,npot,sat,vel,0)
i=0
inquire (file="run.sh",exist=ex)
if (ex) then
 write(char_tmp,"(A,I1)") "bash run.sh ",i
 call system(char_tmp)
else
 stop 'File run.sh does not exist'
endif
inquire (file="take.sh",exist=ex)
if (ex) then
 call takedata(npot,V,dipole,Tpart,0)
else
 stop 'File take.sh does not exist'
endif

Vtmp=V
!!! Getting U-matrix to describe the adiabatic states in terms of input ones (useful to know the required gradients)
U=dcmplx(0.d0,0.d0)
do i=1,npot
 U(i,i)=dcmplx(1.d0,0.d0)
enddo
call diago(npot,Vtmp,dipole,las,U,typdyn,printlevel,.true.,2,68,.false.)

return
end subroutine

subroutine analytical (nat, npot, x, sat, V, dipole, gradV, egradV, graddipole, Tpart, vel)
implicit none
integer, intent(in) :: npot,nat
real*8, intent(in) :: x(nat,3),vel(nat,3)
character*2, intent(in) :: sat(nat)
complex*16, intent(in) :: V(npot,npot)
complex*16, intent(in) :: dipole(3,npot,npot),Tpart(npot,npot)
complex*16, intent(inout) :: gradV(npot,npot,nat,3), graddipole(3,npot,npot,nat,3)
logical, intent(inout) :: egradV(npot,npot)

integer i,j,k
integer ii,ij
logical ex
character*20 fich
real*8 kk(nat,3)

complex*16 oldgradV(npot,npot,nat,3)
character*100 char_tmp

oldgradV=gradV

graddipole=dcmplx(0.d0,0.d0)

!!! Potential was evaluated by singlepoint
!call takedata(x,nat,npot,sat,V,dipole,Tpart,vel,0)
!call senddata(x,nat,npot,sat,vel,0)
!i=0
!write(char_tmp,"(A,I1)") "bash run.sh ",i
!call system(char_tmp)
!call takedata(npot,V,dipole,Tpart,0)

!Potential gradients
inquire(file="run_grad.sh",exist=ex)
if (ex) then
 call system("bash run_grad.sh")
else
 stop 'File run_grad.sh does not exist'
endif

open(1,file="grad.dat")
do 
 read(1,*,end=1) i,j

! write(6,*) "Reading grads ",i,j
 do ii=1,nat
  read(1,*) (kk(ii,ij),ij=1,3)
 enddo
!write(6,*) "Diagonal gradients read"
 if (i.gt.npot .or. j.gt.npot .or. i.le.0 .or. j.le.0) then
  write(6,*) "Requiring gradient for matrix element ",i,j,npot
!  stop 'Gradient for potential not defined'
 else
!!! Initialize it
  egradV(i,j)=.true.
  do ii=1,nat
   do ij=1,3
    gradV(i,j,ii,ij)=dcmplx(0.0d0,0.0d0)
    gradV(i,j,ii,ij)=gradV(i,j,ii,ij)+dcmplx(kk(ii,ij),0.0d0)
   enddo
  enddo
 endif

enddo

1 continue
close(1)

inquire(file="grad.cmp",exist=ex)
if (ex) then
 open(1,file="grad.cmp")
 do
  read(1,*,END=2) i,j

  do ii=1,nat
   read(1,*) (kk(ii,ij),ij=1,3)
  enddo

  if (i.gt.npot .or. j.gt.npot .or. i.le.0 .or. j.le.0) then
   write(6,*) "Requiring gradient for matrix element ",i,j,npot
!  stop 'Gradient for potential not defined'
  else
!!! Initialize it
   egradV(i,j)=.true.
   do ii=1,nat
    do ij=1,3
     gradV(i,j,ii,ij)=gradV(i,j,ii,ij)-dcmplx(0.d0,dimag(gradV(i,j,ii,ij)))
     gradV(i,j,ii,ij)=gradV(i,j,ii,ij)+dcmplx(0.d0,kk(ii,ij))
    enddo
   enddo
  endif


 enddo
endif
2 continue
if (ex) close(1)


901 format(A6,I2.2,A1,I2.2,A4)

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine numerical (nat, npot, x, sat, V, dipole, gradV, egradV, graddipole, Tpart, dh,vel)
implicit none
integer, intent(in) :: npot,nat
real*8, intent(in) :: x(nat,3),vel(nat,3),dh
character*2, intent(in) :: sat(nat)
complex*16, intent(out) :: V(npot,npot),dipole(3,npot,npot),Tpart(npot,npot)
complex*16, intent(out) :: gradV(npot,npot,nat,3)
complex*16, intent(out) :: graddipole(3,npot,npot,nat,3)
logical, intent(out) :: egradV(npot,npot)

complex*16 Vm(npot,npot),Vp(npot,npot)
complex*16 dipolem(3,npot,npot),dipolep(3,npot,npot)
complex*16 Tpartm(npot,npot),Tpartp(npot,npot)

integer i,j,k
integer ii,ij,ntime
real*8 x2(nat,3)

character*100 char_tmp
logical ex

egradV=.true.

ntime=0
do i=1,nat
 do j=1,3
  x2(i,j)=x(i,j)
 enddo
enddo
call senddata(x2,nat,npot,sat,vel,ntime)

do ii=1,nat
 do ij=1,3

  do i=1,nat
   do j=1,3
    x2(i,j)=x(i,j)
   enddo
  enddo
  x2(ii,ij)=x(ii,ij)-dh
  ntime=ntime+1
  call senddata(x2,nat,npot,sat,vel,ntime)

  do i=1,nat
   do j=1,3
    x2(i,j)=x(i,j)
   enddo
  enddo
  x2(ii,ij)=x(ii,ij)+dh
  ntime=ntime+1
  call senddata(x2,nat,npot,sat,vel,ntime)

 enddo
enddo
inquire(file="run.sh",exist=ex)
if (ex) then
 write(char_tmp,"(A,I5)") "bash run.sh ",ntime
 call system(char_tmp)
else
 stop 'File run.sh does not exist'
endif

ntime=0
inquire(file="take.sh",exist=ex)
if (ex) then
 call takedata(npot,V,dipole,Tpart,ntime)
else
 stop 'File take.sh does not exist'
endif

do ii=1,nat
 do ij=1,3
  ntime=ntime+1
  call takedata(npot,Vm,dipolem,Tpart,ntime)
  write(char_tmp,"(A,I0,A)") "overlap",ntime,".dat"
  inquire(file=char_tmp,exist=ex)
  if (ex) then
   call diab_all(char_tmp,npot,Vm,dipolem,Tpart)
  endif
  ntime=ntime+1
  call takedata(npot,Vp,dipolep,Tpart,ntime)
  write(char_tmp,"(A,I0,A)") "overlap",ntime,".dat"
  inquire(file=char_tmp,exist=ex)
  if (ex) then
   call diab_all(char_tmp,npot,Vm,dipolem,Tpart)
  endif

  do i=1,npot
   do j=1,npot
    gradV(i,j,ii,ij)=(Vp(i,j)-Vm(i,j))/2./dh
    do k=1,3
     graddipole(k,i,j,ii,ij)=(dipolep(k,i,j)-dipolem(k,i,j))/2./dh
    enddo
   enddo
  enddo
 enddo
enddo

open(1,file="numgrad.out")
do i=1,npot
 do j=1,npot
  write(1,*) i,j
  do ii=1,nat
   write(1,"(3(x,F10.5))") (gradV(i,j,ii,ij),ij=1,3)
  enddo
 enddo
enddo
close(1)


901 format(A6,I2.2,A1,I2.2,A4)

return
end subroutine

subroutine senddata(x2,nat,npot,sat,vel,ntime)
implicit none
integer nat,npot
real*8 x2(nat,3),vel(nat,3)
character*2 sat(nat)
integer ntime

integer i,j
character*100 char_tmp

if (ntime.le.9) then
 write(char_tmp,"(A,I1.1,A)") "geom",ntime,".xyz"
else if (ntime.le.99) then
 write(char_tmp,"(A,I2.2,A)") "geom",ntime,".xyz"
else
 write(char_tmp,"(A,I3.3,A)") "geom",ntime,".xyz"
endif
open(1,file=char_tmp)

if (ntime.le.9) then
 write(char_tmp,"(A,I1.1,A)") "veloc",ntime,".xyz"
else if (ntime.le.99) then
 write(char_tmp,"(A,I2.2,A)") "veloc",ntime,".xyz"
else
 write(char_tmp,"(A,I3.3,A)") "veloc",ntime,".xyz"
endif
open(2,file=char_tmp)

write(1,*) nat
write(1,*) ""
write(2,*) nat
write(2,*) ""
do i=1,nat
 do j=1,3
  if (dabs(x2(i,j)).le.1e-10) x2(i,j)=0.d0
 enddo
 write(1,900) sat(i),(x2(i,j)*.5292,j=1,3)
 write(2,900) sat(i),(vel(i,j),j=1,3)
enddo
close(1)
close(2)
900 format(A2,3(x,E20.10e3))

return
end

subroutine takedata(npot,V,dipole,Tpart,ntime)
implicit none
integer, intent(in) :: npot
integer, intent(in) :: ntime
complex*16, intent(out) :: V(npot,npot),dipole(3,npot,npot),Tpart(npot,npot)

real*8 kkr(npot),kki(npot)
logical ex1(3),ex2(3),diab
integer i,j,k

! Used to readiabatize the dipoles
real*8 Ar(npot,npot),Ai(npot,npot)
real*8 Ur(npot,npot),Ui(npot,npot)
complex*16 U(npot,npot)
real*8 fv1(npot),fv2(npot),fm1(2,npot)

character*100 char_tmp

!write(6,*) "Running external program"
write(char_tmp,"(A,I5)") "bash take.sh ",ntime
call system(char_tmp)
!write(6,*) "Reading data"

V=dcmplx(0.d0,0.d0)
dipole=dcmplx(0.d0,0.d0)
Tpart=dcmplx(0.d0,0.d0)

!Potential H0
open(1,file="H0.dat")
inquire(file="H0.cmp",exist=ex2(1))
if (ex2(1)) open(2,file="H0.cmp")
do i=1,npot
 do j=1,npot
  kkr(j)=0.d0
  kki(j)=0.d0
 enddo
 read(1,*) (kkr(j),j=1,npot)
 if (ex2(1)) read(2,*) (kki(j),j=1,npot)
 do j=1,npot
  V(i,j)=dcmplx(kkr(j),kki(j))
 enddo
enddo
close(1)
if (ex2(1)) close(2)
!write(6,*) "Potential read"

!NACME
inquire(file="T.dat",exist=ex1(1))
inquire(file="T.cmp",exist=ex2(1))
if (ex1(1)) open(1,file="T.dat")
if (ex2(1)) open(2,file="T.cmp")
do i=1,npot
 do j=1,npot
  kkr(j)=0.d0
  kki(j)=0.d0
 enddo
 if (ex1(1)) read(1,*) (kkr(j),j=1,npot)
 if (ex2(1)) read(2,*) (kki(j),j=1,npot)
 do j=1,npot
  Tpart(i,j)=dcmplx(kkr(j),kki(j))
 enddo
enddo
if (ex1(1)) close(1)
if (ex2(1)) close(2)
!write(6,*) "NACME read"

do k=1,3
 !dipole moments in diabatic
 diab=.false.
 if (k.eq.1) then
  inquire(file="dipole1.dat",exist=ex1(1))
  if (ex1(1)) open(1,file="dipole1.dat")
  inquire(file="dipole1.cmp",exist=ex2(1))
  if (ex2(1)) open(1,file="dipole1.cmp")
 endif
 if (k.eq.2) then
  inquire(file="dipole2.dat",exist=ex1(2))
  if (ex1(2)) open(1,file="dipole2.dat")
  inquire(file="dipole2.cmp",exist=ex2(2))
  if (ex2(2)) open(1,file="dipole2.cmp")
 endif
 if (k.eq.3) then
  inquire(file="dipole3.dat",exist=ex1(3))
  if (ex1(3)) open(1,file="dipole3.dat")
  inquire(file="dipole3.cmp",exist=ex2(3))
  if (ex2(3)) open(1,file="dipole3.cmp")
 endif
 if (ex1(k)) diab=.true.
 if (ex2(k)) diab=.true.
 if (diab) then
  do i=1,npot
   do j=1,npot
    kkr(j)=0.d0
    kki(j)=0.d0
   enddo
   if (ex1(k)) read(1,*) (kkr(j),j=1,npot)
   if (ex2(k)) read(2,*) (kki(j),j=1,npot)
   do j=1,npot
    dipole(k,i,j)=dcmplx(kkr(j),kki(j))
   enddo
  enddo
 else
 ! dipole moments in adiabatic
  if (k.eq.1) then
   inquire(file="dip1r.dat",exist=ex1(1))
   inquire(file="dip1i.dat",exist=ex2(1))
   if (ex1(1)) open(1,file="dip1r.dat")
   if (ex2(1)) open(2,file="dip1i.dat")
  endif
  if (k.eq.2) then
   inquire(file="dip2r.dat",exist=ex1(2))
   inquire(file="dip2i.dat",exist=ex2(2))
   if (ex1(2)) open(1,file="dip2r.dat")
   if (ex2(2)) open(2,file="dip2i.dat")
  endif
  if (k.eq.3) then
   inquire(file="dip3r.dat",exist=ex1(3))
   inquire(file="dip3i.dat",exist=ex2(3))
   if (ex1(3)) open(1,file="dip3r.dat")
   if (ex2(3)) open(2,file="dip3i.dat")
  endif
  do i=1,npot
   do j=1,npot
    kkr(j)=0.d0
    kki(j)=0.d0
   enddo
   if (ex1(k)) read(1,*) (kkr(j),j=1,npot)
   if (ex2(k)) read(2,*) (kki(j),j=1,npot)
   do j=1,npot
    dipole(k,i,j)=dcmplx(kkr(j),kki(j))
   enddo
  enddo
  if (k.eq.1) then
   do i=1,npot
    do j=1,npot
     Ar(i,j)=dreal(V(i,j))
     Ai(i,j)=dimag(V(i,j))
    enddo
   enddo
   !call ch(npot,npot,Ar,Ai,kkr,1,Ur,Ui,fv1,fv2,fm1,i)
   U=V
   call cdiag(U,npot,kkr)
!   do i=1,npot
!    do j=1,npot
!     U(i,j)=dcmplx(Ur(i,j),Ui(i,j))
!    enddo
!   enddo
  endif
  call diabatize(npot,U,dipole(k,:,:))
 endif
enddo
!write(6,*) "Dipole read"

!Potential H1
inquire(file="H1.dat",exist=ex1(1))
inquire(file="H1.cmp",exist=ex2(1))
if (ex1(1)) open(1,file="H1.dat")
if (ex2(1)) open(2,file="H1.cmp")
do i=1,npot
 do j=1,npot
  kkr(j)=0.d0
  kki(j)=0.d0
 enddo
 kkr(i)=dreal(V(i,i))
 kki(i)=dimag(V(i,i))
 if (ex1(1)) read(1,*) (kkr(j),j=1,npot)
 if (ex2(1)) read(2,*) (kki(j),j=1,npot)
 do j=1,npot
  if (i.eq.j) then
   V(i,j)=dcmplx(kkr(j),kki(j))
  else
   V(i,j)=V(i,j)+dcmplx(kkr(j),kki(j))
  endif
 enddo
enddo
if (ex1(1)) close(1)
if (ex2(1)) close(2)

return
end subroutine

