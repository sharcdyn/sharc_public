module diab_
implicit none
contains

subroutine diab_projord(norb,norb2,overlap,ord,rot)
implicit none
integer, intent(in) :: norb,norb2
integer, intent(in) :: ord(norb2)
real*8, intent(in) :: overlap(norb2,norb)
real*8, intent(inout) :: rot(norb,norb)

integer i,j,k,restorb
real*8 ovrot(norb2,norb),rot2(norb,norb)
real*8, allocatable :: proj(:,:),dproj(:)
real*8 tmp_rot(norb,norb)

do k=1,norb2
!!! Creating the new rotated overlap matrix ovrot
! ovrot=matmul(overlap,rot)
 call dmatmat(norb2,norb,norb,overlap,rot,ovrot)
 restorb=norb-k+1
!!! Projecting over the k value, except for the already done
 allocate(proj(restorb,restorb),dproj(restorb))
!$OMP PARALLEL DO COLLAPSE(2)
 do i=1,restorb
  do j=1,restorb
   proj(j,i)=ovrot(ord(k),j)*ovrot(ord(k),i)
  enddo
 enddo
!$OMP END PARALLEL DO
 call diag(proj,restorb,dproj)
 rot2=0.
!!! Creating a rotation matrix wiht 1 in the diagonal and putting the new reduced rotation matrix
!$OMP PARALLEL DO
 do i=1,norb
  rot2(i,i)=1.
 enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO COLLAPSE(2)
 do i=1,restorb
  do j=1,restorb
   rot2(j,i)=proj(j,i)
  enddo
 enddo
!$OMP END PARALLEL DO
! write(6,"(A,2(I0,x),A,F20.8)") "Taking maximum projection of ",k,ord(k)," max proj ",dproj(restorb)
 deallocate(proj,dproj)
! rot=matmul(rot,rot2)
 call dmatmat(norb,norb,norb,rot,rot2,tmp_rot)
 rot=tmp_rot
enddo

if (allocated(proj))  deallocate(proj)
if (allocated(dproj)) deallocate(dproj)


!! Reordering rot --> rot2
!$OMP PARALLEL DO COLLAPSE(2)
do i=1,norb2
! write(6,*) "Possition ",norb-i+1," to ",ord(i)
 do j=1,norb
  rot2(j,ord(i))=rot(j,norb-i+1)
 enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
do i=norb2+1,norb
! write(6,*) "Possition ",i-norb2," to ",i
 do j=1,norb
  rot2(j,i)=rot(j,i-norb2)
 enddo
enddo
!$OMP END PARALLEL DO
!ovrot=matmul(overlap,rot2)
call dmatmat(norb2,norb,norb,overlap,rot2,ovrot)
!$OMP PARALLEL DO PRIVATE(j)
do i=1,norb2
 if (ovrot(i,i).le.0) then
  do j=1,norb
   rot2(j,i)=-rot2(j,i)
  enddo
 endif
enddo
!$OMP END PARALLEL DO
!ovrot=matmul(overlap,rot2)
call dmatmat(norb2,norb,norb,overlap,rot2,ovrot)
!do i=1,norb2
!  write(6,"(A,2(x,I0),2(x,F20.8))") "Final ",i,i,ovrot(i,i),ovrot(i,i)**2
!enddo
rot=rot2

end

subroutine cdiab_projord(norb,norb2,overlap,ord,rot)
implicit none
integer, intent(in) :: norb,norb2
integer, intent(in) :: ord(norb2)
complex*16, intent(in) :: overlap(norb2,norb)
complex*16, intent(inout) :: rot(norb,norb)

integer i,j,k,restorb
complex*16 ovrot(norb2,norb),rot2(norb,norb)
real*8, allocatable :: dproj(:)
complex*16, allocatable :: proj(:,:)
complex*16 tmp_rot(norb,norb)
real*8 phase

do k=1,norb2
!!! Creating the new rotated overlap matrix ovrot
! ovrot=matmul(overlap,rot)
 call cmatmat(norb2,norb,norb,overlap,rot,ovrot)
 restorb=norb-k+1
!!! Projecting over the k value, except for the already done
 allocate(proj(restorb,restorb),dproj(restorb))
!$OMP PARALLEL DO COLLAPSE(2)
 do i=1,restorb
  do j=1,restorb
   proj(j,i)=dconjg(ovrot(ord(k),j))*ovrot(ord(k),i)
  enddo
 enddo
!$OMP END PARALLEL DO
 call cdiag(proj,restorb,dproj)
 rot2=dcmplx(0.d0,0.d0)
!!! Creating a rotation matrix wiht 1 in the diagonal and putting the new reduced rotation matrix
!$OMP PARALLEL DO
 do i=1,norb
  rot2(i,i)=dcmplx(1.d0,0.d0)
 enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO COLLAPSE(2)
 do i=1,restorb
  do j=1,restorb
   rot2(j,i)=proj(j,i)
  enddo
 enddo
!$OMP END PARALLEL DO
! write(6,"(A,2(I0,x),A,F20.8)") "Taking maximum projection of ",k,ord(k)," max proj ",dproj(restorb)
 deallocate(proj,dproj)
! rot=matmul(rot,rot2)
 call cmatmat(norb,norb,norb,rot,rot2,tmp_rot)
 rot=tmp_rot
enddo

!! Reordering rot --> rot2
!$OMP PARALLEL DO COLLAPSE(2)
do i=1,norb2
! write(6,*) "Possition ",norb-i+1," to ",ord(i)
 do j=1,norb
  rot2(j,ord(i))=rot(j,norb-i+1)
 enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
do i=norb2+1,norb
! write(6,*) "Possition ",i-norb2," to ",i
 do j=1,norb
  rot2(j,i)=rot(j,i-norb2)
 enddo
enddo
!$OMP END PARALLEL DO
!ovrot=matmul(overlap,rot2)
call cmatmat(norb2,norb,norb,overlap,rot2,ovrot)
!$OMP PARALLEL DO PRIVATE(j,phase)
do i=1,norb2
 phase=datan2(dimag(ovrot(i,i)),dreal(ovrot(i,i)))
! write(6,*) "Overlap rotated ",i,ovrot(i,i),phase
 do j=1,norb
  rot2(j,i)=rot2(j,i)*cdexp(dcmplx(0.d0,-phase))
 enddo
enddo
!$OMP END PARALLEL DO
!ovrot=matmul(overlap,rot2)
call cmatmat(norb2,norb,norb,overlap,rot2,ovrot)
!do i=1,norb2
!  write(6,"(A,2(x,I0),2(x,F20.8))") "Final ",i,i,ovrot(i,i),ovrot(i,i)**2
!enddo
rot=rot2

end

end module diab_
