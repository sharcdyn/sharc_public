subroutine diag(A,npot,eig)
implicit none
integer, intent(in) :: npot
real*8,intent(inout) :: A(npot,npot)
real*8, intent(out) :: eig(npot)

integer info,lwork
real*8 work(1)
real*8, allocatable :: work2(:)


lwork=-1
call dsyev("V", "U", npot, A, npot, eig, WORK,LWORK,info)
lwork=int(work(1))
allocate (work2(lwork))
call dsyev("V", "U", npot, A, npot, eig, WORK2,LWORK,info)
deallocate (work2)

return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diag_nfunc(A,n,eig,nf)
implicit none
integer, intent(in) :: n,nf
real*8, intent(inout) :: A(n,n)
real*8, intent(out) :: eig(n)

integer info,lwork
real*8, allocatable :: work(:)
real*8 eigv(n,nf)
integer iwork(5*n),ifail(n)
integer i,j

lwork=-1
allocate(work(1))
call dsyevx ('V', 'I', 'U', n, A, n, 0., 0., 1, nf, -1., i, eig, eigv, n, WORK, LWORK, IWORK, IFAIL, INFO)
lwork=int(work(1))
deallocate(work)
allocate(work(lwork))
call dsyevx ('V', 'I', 'U', n, A, n, 0., 0., 1, nf, -1., i, eig, eigv, n, WORK, LWORK, IWORK, IFAIL, INFO)
deallocate(work)

A=0.
do i=1,nf
 do j=1,n
  A(j,i)=eigv(j,i)
 enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine order(A,n,v,desc)
implicit none
integer, intent(in) :: n
real*8, intent(inout) :: A(n,n)
real*8, intent(inout) ::  v(n)
logical, intent(in) :: desc

real*8 kk
integer n2

integer i,j,k

do i=1,n-1
 do j=i+1,n
  if (v(j).le.v(i)) then
   kk=v(j)
   v(j)=v(i)
   v(i)=kk
   do k=1,n
    kk=A(k,j)
    A(k,j)=A(k,i)
    A(k,i)=kk
   enddo
  endif
 enddo
enddo

if (desc) then
 n2=n
 if (mod(n,2).eq.1) n2=n-1
 do i=1,n2/2
  kk=v(i)
  v(i)=v(n-i+1)
  v(n-i+1)=kk
  do k=1,n
   kk=A(k,i)
   A(k,i)=A(k,n-i+1)
   A(k,n-i+1)=kk
  enddo
 enddo
endif


return
end
!!!!!!!!!! COMPLEX
subroutine cdiag(A,npot,eig)
implicit none
integer, intent(in) :: npot
complex*16,intent(inout) :: A(npot,npot)
real*8, intent(out) :: eig(npot)

integer info,lwork
complex*16 work(1)
complex*16, allocatable :: work2(:)
real*8 rwork(3*npot-2)



lwork=-1
call zheev("V", "U", npot, A, npot, eig, work,lwork,rwork,info)
lwork=int(real(work(1)))
allocate (work2(lwork))
call zheev("V", "U", npot, A, npot, eig, work2,lwork,rwork,info)
deallocate(work2)

return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine corder(A,n,v,desc)
implicit none
integer, intent(in) :: n
complex*16, intent(inout) :: A(n,n)
real*8, intent(inout) ::  v(n)
logical, intent(in) :: desc

real*8 kk
complex*16 ckk
integer n2

integer i,j,k

do i=1,n-1
 do j=i+1,n
  if (v(j).le.v(i)) then
   kk=v(j)
   v(j)=v(i)
   v(i)=kk
   do k=1,n
    ckk=A(k,j)
    A(k,j)=A(k,i)
    A(k,i)=kk
   enddo
  endif
 enddo
enddo

if (desc) then
 n2=n
 if (mod(n,2).eq.1) n2=n-1
 do i=1,n2/2
  kk=v(i)
  v(i)=v(n-i+1)
  v(n-i+1)=kk
  do k=1,n
   ckk=A(k,i)
   A(k,i)=A(k,n-i+1)
   A(k,n-i+1)=kk
  enddo
 enddo
endif


return
end


!!!!!!!!!!!
! Calculating eigenvalues from a non-orthonormal basis using the overlap
!! n is the size
!! ov enters as the overlap and leaves as the transformation matrix to a orthonormal basis
!! A enters as the Hamiltonian and leaves as the rotation matrix
!! w contains the eigenvalues
subroutine eig(n,A,ov,w)
implicit none
integer, intent(in) :: n
real*8, intent(inout) :: A(n,n),ov(n,n)
real*8, intent(out) :: w(n)

real*8 :: thr=1e-5
real*8 ov0(n,n)

integer i,j,k,l
integer nli
real*8, allocatable :: tmpmat(:,:),tmpvec(:)
real*8, allocatable :: tmpmat2(:,:),tmpmat3(:,:)
real*8 kk

ov0=ov
call diag(ov,n,w)
call order(ov,n,w,.true.)
nli=0
do i=1,n
 if (w(i).ge.thr) then
  nli=nli+1
  do j=1,n
   ov(j,i)=ov(j,i)/sqrt(w(i))
  enddo
 else
  if (.false.) write(6,*) "Remove linear dependencies ",w(i)
  do j=1,n
   ov(j,i)=0.
  enddo
 endif
enddo

allocate(tmpmat(nli,nli),tmpvec(nli) )
allocate(tmpmat2(n,nli),tmpmat3(n,nli) )


if (.false.) then
!!!Checking othonormality

 write(6,*) "Linear independent ",nli," from ",n
 do i=1,nli
  do j=1,nli
   tmpmat(i,j)=0.
   do k=1,n
    do l=1,n
     tmpmat(i,j)=tmpmat(i,j)+ov(k,i)*ov(l,j)*ov0(k,l)
    enddo
   enddo
   write(6,*) "Done ",i,j,tmpmat(i,j)
  enddo
 enddo
 kk=0.
 do i=1,nli
  do j=1,nli
   if (i.eq.j) then
    kk=kk+(1-tmpmat(i,j))**2
   else
    kk=kk+tmpmat(i,j)**2
   endif
  enddo
 enddo
 write(6,*) "New basis functions ",kk

endif

do i=1,nli
 do j=1,n
  tmpmat2(j,i)=ov(j,i)
 enddo
enddo

call dmatmat(n,n,nli,A,tmpmat2,tmpmat3)
call dmatmat(nli,n,nli, transpose(tmpmat2),tmpmat3,tmpmat)
tmpmat=matmul(transpose(tmpmat2),tmpmat3)


!do i=1,nli
! do j=1,nli
!  tmpmat(i,j)=0.
!  do k=1,n
!   do l=1,n
!    tmpmat(i,j)=tmpmat(i,j)+ov(k,i)*ov(l,j)*A(k,l)
!   enddo
!  enddo
! enddo
!enddo

call diag(tmpmat,nli,tmpvec)
A=0.
w=0.
do i=1,nli
 w(i)=tmpvec(i)
 do j=1,nli
  A(j,i)=tmpmat(j,i)
 enddo
enddo


end


subroutine cmatmat(n1,n2,n3,mat1,mat2,matr)
implicit none
integer, intent(in) :: n1,n2,n3
complex*16, intent(in) :: mat1(n1,n2),mat2(n2,n3)
complex*16, intent(out) :: matr(n1,n3)

complex*16 zone,zzero
zone=dcmplx(1.d0,0.d0)
zzero=dcmplx(0.d0,0.d0)

!matr=matmul(mat1,mat2)
call zgemm("N","N",n1,n3,n2,zone,mat1,n1,mat2,n2,zzero,matr,n1)
!call zhemm('L','U',n1,n3,zone,mat1,n1,mat2,n2,zzero,matr,n1)

end subroutine cmatmat


subroutine dmatmat(n1,n2,n3,mat1,mat2,matr)
implicit none
integer, intent(in) :: n1,n2,n3
real*8, intent(in) :: mat1(n1,n2),mat2(n2,n3)
real*8, intent(out) :: matr(n1,n3)

real*8 one,zero
one=1.d0
zero=0.d0

!matr=matmul(mat1,mat2)
call dgemm("N","N",n1,n3,n2,one,mat1,n1,mat2,n2,zero,matr,n1)
!call zhemm('L','U',n1,n3,zone,mat1,n1,mat2,n2,zzero,matr,n1)

end subroutine dmatmat

subroutine cmatvec(n1,n2,mat,vec)
implicit none
integer, intent(in) :: n1,n2
complex*16, intent(in) :: mat(n1,n2)
complex*16, intent(inout) :: vec(n2)

complex*16 vecr(n2),zone,zzero

zone=dcmplx(1.d0,0.0d0)
zzero=dcmplx(0.d0,0.d0)

call zgemv("N",n1,n2,zone,mat,n1,vec,1,zzero,vecr,1)

vec=vecr

end subroutine cmatvec

