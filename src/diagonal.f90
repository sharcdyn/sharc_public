subroutine diago(npot,V,dipole,las,U,typdyn,printlevel,reorder,idiag_type,un,partdiag)
implicit none
integer, intent(in) ::npot,typdyn, un
complex*16, intent(in) :: dipole(3,npot,npot)
real*8, intent (in) :: las(3)
integer, intent(in) :: idiag_type
logical, intent(in) :: reorder, partdiag
integer, intent(in) :: printlevel
complex*16, intent(inout) :: V(npot,npot)
complex*16, intent(inout) :: U(npot,npot)

logical writtmp

real*8 ignacios, ignacios2
real*8 VW1r(npot,npot),VW1i(npot,npot)
real*8 VW2r(npot,npot),VW2i(npot,npot)
real*8 VW(npot)
real*8 fv1(npot),fv2(npot),fm1(2,npot)
integer ierr,diag_type
complex*16 U1(npot,npot),U2(npot,npot),W(npot,npot)

complex*16 kki
real*8 kkr



integer i,j,k

diag_type=idiag_type

!ch
!diag_type=1
!jacobi
!diag_type=2

ignacios=0.d0


writtmp=.false.
if (printlevel.ge.3) writtmp=.true.

do i=1,3
 do j=1,npot
  do k=1,npot
   V(j,k)=V(j,k)-dipole(i,j,k)*dcmplx(las(i),0.d0)
  enddo
 enddo
enddo

!   !!! For FISH dynamics
if (iabs(typdyn).lt.20.and.iabs(typdyn).ge.10) then
 do i=1,npot
  do j=1,npot
   U(i,j)=dcmplx(0.d0,0.d0)
  enddo
  U(i,i)=dcmplx(1.d0,0.d0)
 enddo
else

!   !!! Diagonalizing SHARC method and Born-Oppenheimer Dynamics
 U1=U
 W=V
 call adiabatize(npot,U1,W)

 U2=W
 do i=1,npot
  do j=1,npot
!!! Ignacio's suggestion, let's remove from the diagonalization small couplings (1e-6)
   if (partdiag.and.i.ne.j) then
    if (cdabs(U2(i,j)).lt.ignacios) then
     if (cdabs(U2(i,i)-U2(j,j)).lt.ignacios) then
      if (printlevel.ge.3) write(999,901) "Ignacio aplied", i, j, U2(i,i), U2(j,j), U2(i,j)
      U2(i,j)=dcmplx(0.d0,0.d0)
     endif
    endif
   endif
  enddo
 enddo

 if (diag_type.eq.2) then
!!! U2 potential adiabatized with previous U
!!! VW the eigenvalues
!!! U is the rotation matrix
!!! Jacobi diagonalization without reorder the energies
  call Heigensystem(npot,U2,npot,VW,U,npot,0)
!!! The eigenvalues are on rows... let's change to columns to be
!!! consistent with the rest of the program
  do i=1,npot
   do j=1,npot
    U2(i,j)=dconjg(U(j,i))
   enddo
  enddo
 endif

 if (diag_type.eq.1) then

  !call ch(npot,npot,VW1r,VW1i,VW,1,VW2r,VW2i,fv1,fv2,fm1,ierr)
  call cdiag(U2,npot,VW)

 endif

 if (writtmp) then
  write(667,"(A,1000(x,F4.2))") "U2 bef ",(cdabs(U2(i,i))**2,i=1,npot)
  do i=1,npot
   write(667,900) U2(i,:)
  enddo

  W=V
  call cmatmat(npot,npot,npot,U1,U2,U)
  call adiabatize (npot,U,W)
  write(666,*) "Adiabatic before order U2"
  do i=1,npot
   write(666,900) W(i,:)
  enddo
 endif


! call degener(npot,U2,VW)
! call orderU2(npot,U2,un)
! call cmatmat(npot,npot,npot,U1,U2,U)
 call orderU2_old(npot,U2,un)
 call signo(npot,U1,U2,U,writtmp)
 if (printlevel.ge.3) then
  W=V
  call adiabatize (npot,U,W)
  write(666,*) "Adiabatic after order and sign U2"
  do i=1,npot
   write(666,900) (dreal(W(i,j)),dimag(W(i,j)),j=1,npot)
  enddo
 endif
 if (reorder) then
  W=V
  call adiabatize(npot,U,W)
  do i=1,npot
   VW(i)=dreal(W(i,i))
  enddo

! Reordering
  do i=1,npot
   do j=i,npot
    if ( (VW(j).le.VW(i)) .and. (i.ne.j) ) then
     kkr=VW(i)
     VW(i)=VW(j)
     VW(j)=kkr
     do k=1,npot
      kki=U(k,i)
      U(k,i)=U(k,j)
      U(k,j)=kki
     enddo
    endif

   enddo
  enddo
! else
!  call diab_potential(npot,U2)
 endif

 if (writtmp) then
  write(667,"(A,1000(x,F4.2))") "U2 after ",(cdabs(U2(i,i))**2,i=1,npot)
  do i=1,npot
   write(667,900) (dreal(U2(i,j)),dimag(U2(i,j)),j=1,npot)
  enddo
 endif



endif

if (printlevel.ge.3) then
 W=V
 call adiabatize(npot,U,W)
 write(666,*) "Hamiltonian after diago"
 do i=1,npot
  write(666,*) (W(i,j),j=1,npot)
 enddo
endif

900 format(100000(E20.10e3))
901 format(A35,x,I2,x,I2,100000(E20.10e3))

return
end subroutine

subroutine orderU2_old(npot,U2,un)
implicit none
integer, intent(in) :: npot,un
complex*16, intent(inout) :: U2(npot,npot)

complex*16 U(npot,npot)
complex*16 kk
integer i,j,k
real*8 maxim
integer ibest


do i=1,npot
!!! Taking maximum overlap with old state i (bra)
 maxim=cdabs(U2(i,i))**2
 do j=1,npot
  if (cdabs(U2(i,j))**2.ge.maxim) then
   maxim=cdabs(U2(i,j))**2
   ibest=j
  endif
 enddo
! Check if it is at least 50% better than original
 if ( ibest.ne.i.and. cdabs(U2(i,ibest))**2.gt.1.5*cdabs(U2(i,i))**2 ) then
  U=U2
 ! Exchange i and ibest
  do j=1,npot
   kk=U2(j,i)
   U2(j,i)=U2(j,ibest)
   U2(j,ibest)=kk
  enddo
  write(un,"(A9,2(x,I4.4),4(x,E20.10e3))") "Changing ",i,ibest,U(i,i),U2(i,i)
 endif
enddo
do i=1,npot
 if (cdabs(U2(i,i))**2.le. 0.8) write(un,*) "Small overlap between consecutive steps ",i,cdabs(U2(i,i))**2,U2(i,i)
enddo


900 format(1000(x,F10.4))
end subroutine orderU2_old


subroutine orderU2(npot,U2,un)
implicit none
integer, intent(in) :: npot,un
complex*16, intent(inout) :: U2(npot,npot)

complex*16 tmp1(npot,npot),tmp2(npot,npot)
integer i,j,k
real*8 maxim
integer irow,icol

logical checked(npot,npot)
integer maxv(2)
real*8 kk(npot,npot)

kk=cdabs(U2)**2
tmp2=dcmplx(0.d0,0.d0)
checked=.true.
do k=1,npot
!! Taking maximum value
 maxv=maxloc(kk,checked)
 i=maxv(1)
 j=maxv(2)
 checked(i,:)=.false.
 checked(:,j)=.false.
 tmp2(i,j)=dcmplx(1.d0,0.d0)
 if (i.ne.j) then
  write(un,"(A9,2(x,I4.4),4(x,E20.10e3))") "Changing ",i,j,kk(i,j)
 endif
enddo

!! Lets put in tmp1 the U2 permutated
call cmatmat(npot,npot,npot,U2,transpose(dconjg(tmp2)),tmp1)

!Removing phase
do i=1,npot
 maxim=datan2( dimag(tmp1(i,i)),dreal(tmp1(i,i)) )
 U2(:,i)=tmp1(:,i)*cdexp(dcmplx(0.d0,-maxim))
enddo

do i=1,npot
 if (cdabs(U2(i,i))**2.le. 0.8) write(un,*) "Small overlap between consecutive steps ",i,cdabs(U2(i,i))**2,U2(i,i)
enddo


end subroutine orderU2
   

subroutine signo(npot,U1,U2,U,writU)
implicit none
integer npot
complex*16 U1(npot,npot),U2(npot,npot),U(npot,npot)

real*8 best,phase(npot)
real*8 best2,phase2(npot)
integer ibest,ibest2
real*8 phase3,phase4,r

integer i,j,k
logical writU

do i=1,npot
 best=0
 best2=0
 do j=1,npot
  if (cdabs(U1(j,i))**2.ge.best) then
   best=cdabs(U1(j,i))**2
   ibest=j
  endif
  if (cdabs(U2(j,i))**2.ge.best2) then
   best2=cdabs(U2(j,i))**2
   ibest2=j
  endif
 enddo

 ! Phase of U1
! r=cdabs(U1(ibest,i))
 phase(i)=datan2(dimag(U1(ibest,i)),dreal(U1(ibest,i)))
! if (dimag(U1(ibest,i)).lt.0) phase(i)=-phase(i)
 ! Pĥase of U2
! r=cdabs(U2(ibest2,i))
 phase2(i)=datan2(dimag(U2(ibest2,i)),dreal(U2(ibest2,i)))
! if (dimag(U2(ibest2,i)).lt.0) phase2(i)=-phase2(i)

!  Checking sign by comparing phase 0 with phase pi/2
! if ( phase2(i).gt.dacos(-1.d0) ) phase2(i)=phase2(i)-2.d0*dacos(-1.d0)
! if ( phase2(i)**2.gt.(phase2(i)-dacos(-1.d0))**2 )  then
  do j=1,npot
!   U2(j,i)=U2(j,i)*cdexp(dcmplx(0.d0,-dacos(-1.d0)))
!   U2(j,i)=U2(j,i)
   U2(j,i)=U2(j,i)*cdexp(dcmplx(0.d0,-phase2(i)))
  enddo
!  write(69,*) "Changed sign ",U2(ibest2,i)
! endif
enddo
do i=1,npot
 do j=1,npot
  U(i,j)=dcmplx(0.d0,0.d0)
  do k=1,npot
   U(i,j)=U(i,j)+U1(i,k)*U2(k,j)
  enddo
 enddo
enddo

 if (writU) write(666,*) "pot prev bef after total"
do i=1,npot
 best=0.
 best2=0.
 do j=1,npot
  if (cdabs(U(j,i))**2.ge.best) then
   best=cdabs(U(j,i))**2
   ibest=j
  endif
  if (cdabs(U2(j,i))**2.ge.best2) then
   best2=cdabs(U2(j,i))**2
   ibest2=j
  endif
 enddo
 phase3=datan2(dimag(U2(ibest2,i)),dreal(U2(ibest2,i)))

 phase4=datan2(dimag(U(ibest,i)),dreal(U(ibest,i)))

 if (writU) write(666,900) i,phase(i)/dacos(-1.d0),phase2(i)/dacos(-1.d0),phase3/dacos(-1.d0),phase4/dacos(-1.d0)
enddo

900 format(I3.3,4(x,F5.2))
return
end subroutine

subroutine timederiv(npot,deriv2,U,U_prev,dt)
implicit none
integer, intent(in) :: npot
complex*16, intent(in) :: U(npot,npot),U_prev(npot,npot)
complex*16, intent(out) :: deriv2(npot,npot)
real*8, intent(in) :: dt
complex*16 kk(npot,npot)

kk=U-U_prev
call cmatmat(npot,npot,npot,transpose(dconjg(U)),kk,deriv2)
deriv2=deriv2/dt


900 format(100(x,E20.10e3))

return
end subroutine

subroutine reordering(npot,V,coef,pot,U_prev,U,time)
implicit none
integer, intent(in) :: npot
integer, intent(inout) :: pot
real*8, intent(in) :: time
complex*16, intent(inout) :: V(npot,npot),coef(npot)
complex*16, intent(inout) :: U_prev(npot,npot),U(npot,npot)

complex*16 kkmat(npot,npot)
real*8 W(npot)

integer i,j,k,l
complex*16 kk
real*8 mini
integer ibest,ipot


!write(100,*) "Diabatic pots"
!do i=1,npot
! write(100,"(100(E20.10e3))") V(i,:)
!enddo

kkmat=V
! Getting adiabatic potentials
call adiabatize(npot,U,kkmat)
mini=0.d0
do i=1,npot
 W(i)=dreal(kkmat(i,i))
 kkmat(i,i)=kkmat(i,i)-dcmplx(W(i),0.d0)
 do j=1,npot
  mini=mini+cdabs(kkmat(i,j))**2
 enddo
enddo

if (mini.ge.1e-8) then
 write(69,"(A,E20.10)") &
  "Warning the reordering is not happening in an adiabatic potential, off diag square ",mini
endif

!write(100,*) "Adiabatic with no reordering pot"
!do i=1,npot
! write(100,"(100(E20.10e3))") W(i,:)
!enddo

!mini=0.d0
!do i=1,npot
! do j=1,npot
!  if (i.ne.j) mini=mini+cdabs(W(i,j))
! enddo
!enddo
!write(668,*) "Off diagonal elements", mini

ipot=pot
do i=1,npot
 ibest=i
 mini=W(ibest)
 do j=i+1,npot
  if (W(j).lt.mini) then
   mini=W(j)
   ibest=j
  endif
 enddo
 if (i.ne.ibest) then
 ! Exchanging i,ibest
  j=ipot
  k=ipot
  if (ipot.eq.i) k=ibest
  if (ipot.eq.ibest) k=i
  ipot=k
  write(69,"(A11,F20.10,4(x,I4.4))") "Reordering ",time/41.32,i,ibest,j,ipot
! Changing adiabatic properties
  mini=W(i)
  W(i)=W(ibest)
  W(ibest)=mini
  kk=coef(i)
  coef(i)=coef(ibest)
  coef(ibest)=kk
  do k=1,npot
   kk=U(k,i)
   U(k,i)=U(k,ibest)
   U(k,ibest)=kk
   kk=U_prev(k,i)
   U_prev(k,i)=U_prev(k,ibest)
   U_prev(k,ibest)=kk
  enddo
  write(69,"(4(x,E20.10e3))") W(i),W(ibest)
 endif
enddo
pot=ipot

return
end

subroutine degener(npot,U2,W)
use diab_
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: W(npot)
complex*16, intent(inout) :: U2(npot,npot)

complex*16 V(npot,npot)
complex*16 rot(npot,npot)
integer degen(npot),ndeg
real*8 Edeg(npot)
integer deg_pot(npot)
integer i,j,k
logical nodeg

integer, allocatable :: diaborder(:),pots(:)
complex*16, allocatable :: over(:,:),minirot(:,:)
integer minipot

! do i=1,npot
!  write(6,*) "V ",i,W(i)
! enddo
 
degen(1)=1
ndeg=1
Edeg(1)=W(1)
deg_pot(1)=1
do i=2,npot
 nodeg=.true.
 do j=1,ndeg
  if (dabs(Edeg(j)-W(i)) .le.1e-4 ) then
   nodeg=.false.
   degen(i)=j
   deg_pot(j)=deg_pot(j)+1
  endif
 enddo
 if (nodeg) then
  ndeg=ndeg+1
  Edeg(ndeg)=W(i)
  deg_pot(ndeg)=1
  degen(i)=ndeg
 endif
enddo

! write(6,*) "Degeneration ",ndeg
! do i=1,ndeg
!  write(6,*) i,Edeg(i),deg_pot(i)
! enddo
! do i=1,npot
!  write(6,*) "State ",i,degen(i),W(i),"diag",U2(i,i),"rest",U2(:,i)
! enddo
 
!! In general nothing is gonna be rotated, so an unit matrix is created (rot)
rot=dcmplx(0.d0,0.d0)
do j=1,npot
 rot(j,j)=dcmplx(1.d0,0.d0)
enddo

do i=1,ndeg
!! Lets do a reduced overlap to diabatize
 if (deg_pot(i).ne.1) then
!! A global rotation matrix

  j=deg_pot(i)
  allocate(over(j,j),minirot(j,j))
  allocate(diaborder(j),pots(j))
  minipot=j
  minirot=dcmplx(0.d0,0.d0)
!! A reduced rotation matrix
  do j=1,minipot
   diaborder(j)=j
   minirot(j,j)=dcmplx(1.d0,0.d0)
  enddo
  k=0
!! Potentials in the subset
  do j=1,npot
   if (degen(j).eq.i) then
    k=k+1
    pots(k)=j
   endif
  enddo
!! Overlap in the subset
  do j=1,minipot
   do k=1,minipot
    over(j,k)=U2(pots(j),pots(k))
   enddo
  enddo

!  write(6,*) "Subset ",i,deg_pot(i),pots(:)
!   write(6,*) "Over"
!   do j=1,minipot
!    write(6,*) (over(j,:))
!   enddo

  call cdiab_projord(minipot,minipot,over,diaborder,minirot)
!   write(6,*) "Rot"
!   do j=1,minipot
!    write(6,*) (minirot(j,:))
!   enddo

  do j=1,minipot
   do k=1,minipot
    rot(pots(j),pots(k))=minirot(j,k)
   enddo
  enddo

  deallocate(over,minirot,diaborder,pots)
 endif

enddo

call cmatmat(npot,npot,npot,U2,rot,V)

! write(6,*) "new U2"
! do i=1,npot
!  write(6,*) "U2",i,degen(i),"diag",V(i,i),"rest",V(i,:)
! enddo
U2=V

end subroutine degener

subroutine diab_potential(npot,U2)
use diab_
implicit none
integer, intent(in) :: npot
complex*16, intent(inout) :: U2(npot,npot)

integer i
integer ord(npot) 
complex*16 rot(npot,npot)

rot=dcmplx(0.d0,0.d0)
do i=1,npot
 ord(i)=i
 rot(i,i)=dcmplx(1.d0,0.d0)
enddo

call cdiab_projord(npot,npot,U2,ord,rot)
U2=rot

end subroutine diab_potential
