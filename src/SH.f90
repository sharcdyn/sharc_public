subroutine surfacehopp (npot, pot, coef, delta, inter, prob, step, typsh,time)
implicit none
integer, intent(in) :: npot, step, typsh, pot
real*8, intent(in) :: delta, time
complex*16, intent(in) :: coef(npot), inter(npot, npot)
real*8, intent(inout) ::prob(npot)  !!! In the first step is output only, in the second input (output of the first step) and output (probability)

complex*16 kk
real*8 ropunt(npot), suma
real*8 prob1(npot), probt
integer i,j, teta(npot)


if (step.eq.0) then
 do i=1,npot
  prob(i)=0.d0
 enddo
endif

! advance Tully SH
if (iabs(typsh).eq.1) then
 if (step.eq.1) then
  do i=1,npot
   prob(i)=0.d0
   kk=dconjg(coef(pot))*coef(i)*(inter(pot,i))
   prob(i)=-2*dreal(kk)
   prob(i)=prob(i)/dreal(dconjg(coef(pot))*coef(pot))*delta
  enddo
 else
  probt=0.d0
  do i=1,npot
   prob1(i)=prob(i)
   kk=dconjg(coef(pot))*coef(i)*(inter(pot,i))
   prob(i)=-2*dreal(kk)
   prob(i)=prob(i)/dreal(dconjg(coef(pot))*coef(pot))*delta

!   !!! Integral before and after coefficient propagation
   prob(i)=(prob1(i)+prob(i))/2.

   if (prob(i).lt.0) prob(i)=0.d0
   if (i.eq.pot) prob(i)=0.d0
  enddo
 endif
endif


! Regular Tully
if (iabs(typsh).eq.2.and.step.eq.2) then
 probt=0.d0
 do i=1,npot
  prob(i)=0.d0
  kk=dconjg(coef(pot))*coef(i)*(inter(pot,i))
  prob(i)=-2*dreal(kk)
  prob(i)=prob(i)/dreal(dconjg(coef(pot))*coef(pot))*delta
!   !!! Integral after coefficient propagation
  if (prob(i).lt.0) prob(i)=0.d0
  if (i.eq.pot) prob(i)=0.d0
 enddo
endif

! Mitric SH

if(iabs(typsh).eq.3.and.step.eq.2) then
 do i=1,npot
  kk=dcmplx(0.d0,0.d0)
  do j=1,npot
   kk=kk+inter(i,j)*coef(j)
  enddo
  ropunt(i)=2*dreal(dconjg(coef(i))*kk)
 enddo
endif

if(iabs(typsh).eq.4) then
! Storing ro(ii(t)) as prob
 if(step.eq.1) then
  do i=1, npot
   prob(i)=dreal(dconjg(coef(i))*coef(i))
  enddo
 else
  do i=1,npot
   ropunt(i)=dreal(dconjg(coef(i))*coef(i)) -  prob(i)
   ropunt(i)=ropunt(i)/delta
  enddo
 endif
endif

if (iabs(typsh).eq.3.or.iabs(typsh).eq.4) then
 if (step.eq.2) then
  suma=0.d0
  do i=1,npot
   if(ropunt(i).gt.0.d0) then
    teta(i)=1
    suma=suma+ropunt(i)
   else
    teta(i)=0
   endif
  enddo
  do i=1,npot
   if(-ropunt(pot).gt.0.d0.and.teta(i).eq.1) then
    kk=( -ropunt(pot)/( dconjg(coef(pot))*coef(pot) ) )
    kk=kk*ropunt(i)/suma*delta
    prob(i)=dreal(kk)
   else
    prob(i)=0.d0
   endif
  enddo
 endif
endif

!   !!! Renormalization removing impossible things
if (step.eq.2) then
 probt=0.d0
 do i=1,npot
  probt=probt+prob(i)
 enddo
 if (probt.gt.1) then
  write(69,*) "Warning total probability more than 1 ",probt,time/41.32
  do i=1,npot
   prob(i)=prob(i)/probt
  enddo
 endif
endif


return
end subroutine
