subroutine splining(time,dt,V,T,dipole,npot,npol,Vsp,Tsp,dipsp,printlevel,nsp)
implicit none
integer,intent(in) :: npot,npol,nsp,printlevel
real*8, intent(in) :: time,dt
complex*16, intent(in) :: V(npol,npot,npot),T(npol,npot,npot),dipole(npol,3,npot,npot)
complex*16,intent(out) :: Vsp(4,npot,npot),Tsp(4,npot,npot),dipsp(3,4,npot,npot)

real*8 timeor(npol),kkor(npol)


integer i,j,k,l
real*8 asp(npol),bsp(npol),csp(npol),dsp(npol),asp2(npol),bsp2(npol),csp2(npol),dsp2(npol)

real*8 phase(npol),ampl(npol),aa,ap

complex*16 mattmp(npot,npot)
real*8 timesp


! Splines
do i=1,npol
 timeor(npol-i+1)=time-(i-1)*dt
enddo

Vsp=dcmplx(0.d0,0.d0)
Tsp=dcmplx(0.d0,0.d0)
dipsp=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO COLLAPSE(3) 
do i=1,npot
 do j=1,npot
  do k=-1,3
   if (k.eq.-1) Vsp(1,i,j)=V(1,i,j)
   if (k.eq.0)  Tsp(1,i,j)=T(1,i,j)
   if (k.gt.0)  dipsp(k,1,i,j)=dipole(1,k,i,j)
  enddo
 enddo
enddo
!$OMP END PARALLEL DO

! write(100,900) time,dt,timeor(:)
! write(101,900) time,delta/2.,timesp(:)

  ! y =ax+b
if (npol.ge.2) then
 !$OMP PARALLEL DO COLLAPSE(3) 
 do i=1,npot
  do j=1,npot
   do k=-1,3
    if (k.eq.-1) Vsp(1,i,j)=V(2,i,j)
    if (k.eq.0)  Tsp(1,i,j)=T(2,i,j)
    if (k.gt.0)  dipsp(k,1,i,j)=dipole(2,k,i,j)
   enddo
  enddo
 enddo
 !$OMP END PARALLEL DO
 do i=1,npot
  do j=1,npot

  ! Potential
!  call ampphase(npol,V(:,i,j),ampl(:),phase(:) )
   if (npol.eq.2) then
!   aa=(ampl(1)-ampl(2))/dt
!   ap=(phase(1)-phase(2))/dt
!   a=dcmplx(aa*dcos(ap),aa*dsin(ap))
    
    Vsp(2,i,j)=(V(1,i,j)-V(2,i,j))/dt
   else
  ! Potential
    do k=1,npol
     kkor(npol-k+1)=dreal(V(k,i,j))
!   kkor(npol-k+1)=ampl(k)
    enddo
    call cspline(npol,npol,timeor,kkor,asp,bsp,csp,dsp)
    do k=1,npol
     kkor(npol-k+1)=dimag(V(k,i,j))
!    kkor(npol-k+1)=phase(k)
    enddo
    call cspline(npol,npol,timeor,kkor,asp2,bsp2,csp2,dsp2)
    Vsp(4,i,j)=dcmplx( asp(npol-1),asp2(npol-1) )
    Vsp(3,i,j)=dcmplx( bsp(npol-1),bsp2(npol-1) )
    Vsp(2,i,j)=dcmplx( csp(npol-1),csp2(npol-1) )
    Vsp(1,i,j)=dcmplx( dsp(npol-1),dsp2(npol-1) )
   endif

   ! Dipoles
   do l=1,3
!   call ampphase(npol,dipole(:,l,i,j),ampl(:),phase(:) )
    if (npol.eq.2) then
!    aa=(ampl(1)-ampl(2))/dt
!    ap=(phase(1)-phase(2))/dt
!    a=dcmplx(aa*dcos(ap),aa*dsin(ap))
     dipsp(l,2,i,j)=(dipole(1,l,i,j)-dipole(2,l,i,j))/dt
     dipsp(l,1,i,j)=dipole(2,l,i,j)
    else
     do k=1,npol
!     kkor(npol-k+1)=ampl(k)
      kkor(npol-k+1)=dreal(dipole(k,l,i,j))
     enddo
     call cspline(npol,npol,timeor,kkor,asp,bsp,csp,dsp)
     do k=1,npol
 !     kkor(npol-k+1)=phase(k)
      kkor(npol-k+1)=dimag(dipole(k,l,i,j))
     enddo
     call cspline(npol,npol,timeor,kkor,asp2,bsp2,csp2,dsp2)
     dipsp(l,4,i,j)=dcmplx(asp(npol-1),asp2(npol-1))
     dipsp(l,3,i,j)=dcmplx(bsp(npol-1),bsp2(npol-1))
     dipsp(l,2,i,j)=dcmplx(csp(npol-1),csp2(npol-1))
     dipsp(l,1,i,j)=dcmplx(dsp(npol-1),dsp2(npol-1))
    endif
   enddo

   ! NACMEs
!  call ampphase(npol,T(:,i,j),ampl(:),phase(:) )
   if (npol.eq.2) then
!   aa=(ampl(1)-ampl(2))/dt
!   ap=(phase(1)-phase(2))/dt
!   a=dcmplx(aa*dcos(ap),aa*dsin(ap))
    Tsp(3,i,j)=(T(1,i,j)-T(2,i,j))/dt
    Tsp(4,i,j)=T(2,i,j)
   else
    do k=1,npol
 !    kkor(npol-k+1)=ampl(k)
     kkor(npol-k+1)=real(T(k,i,j))
    enddo
    call cspline(npol,npol,timeor,kkor,asp,bsp,csp,dsp)
    do k=1,npol
 !    kkor(npol-k+1)=phase(k)
     kkor(npol-k+1)=dimag(T(k,i,j))
    enddo
    call cspline(npol,npol,timeor,kkor,asp2,bsp2,csp2,dsp2)
    Tsp(4,i,j)=dcmplx(asp(npol-1),asp2(npol-1))
    Tsp(3,i,j)=dcmplx(bsp(npol-1),bsp2(npol-1))
    Tsp(2,i,j)=dcmplx(csp(npol-1),csp2(npol-1))
    Tsp(1,i,j)=dcmplx(dsp(npol-1),dsp2(npol-1))
   endif
  enddo
 enddo
endif



!!!!!!!!!!!!!!!!!!!!!!
if (.false.) then
do k=1,4
 ! Creating an hermitian matrix for V and dip
 do i=1,npot
  do j=1,npot
   mattmp(i,j)=Vsp(k,i,j)
  enddo
 enddo
 do i=1,npot
  do j=1,i
   Vsp(k,i,j)=(mattmp(i,j)+dconjg(mattmp(j,i)))/2.
   Vsp(k,j,i)=dconjg(Vsp(k,i,j))
  enddo
 enddo
 do l=1,3
  do i=1,npot
   do j=1,npot
    mattmp(i,j)=dipsp(l,k,i,j)
   enddo
  enddo
  do i=1,npot
   do j=1,i
    dipsp(l,k,i,j)=(mattmp(i,j)+dconjg(mattmp(j,i)))/2.
    dipsp(l,k,j,i)=dconjg(dipsp(l,k,i,j))
   enddo
  enddo
 enddo
 ! Creating antihermitian for the coupling
 do i=1,npot
  do j=1,npot
   mattmp(i,j)=Tsp(k,i,j)
  enddo
 enddo
 do i=1,npot
  do j=1,i
   Tsp(k,i,j)=(mattmp(i,j)-dconjg(mattmp(j,i)))/2.
   Tsp(k,j,i)=-dconjg(Tsp(k,i,j))
  enddo
 enddo

enddo
endif

if (printlevel.ge.2) then
 write(100,900) time-dt,(dreal(V(2,i,i)),i=1,npot)
 write(101,900) time-dt,((V(2,i,j),j=1,npot),i=1,npot)
 do i=0,nsp
  timesp=time-dt+i*dt/float(nsp)
  call evaluatespline(npot,Vsp,timesp,time-dt,mattmp)
  write(102,900) timesp,(dreal(mattmp(j,j)),j=1,npot)
  write(103,900) timesp,((mattmp(j,k),k=1,npot),j=1,npot)
 enddo
 write(100,900) time,(dreal(V(1,i,i)),i=1,npot)
 write(101,900) time,((V(1,i,j),j=1,npot),i=1,npot)
endif

 900 format(1000000(E20.10e3))

return
end subroutine splining

subroutine ampphase(npol,vec,ampl,phase )
implicit none
integer, intent(in) :: npol
complex*16, intent(in) :: vec(npol)
real*8, intent(out) :: ampl(npol),phase(npol)

integer k
real*8 dif,dtest
real*8 pi

pi=dacos(-1.d0)

!$OMP PARALLEL DO
do k=1,npol
 ampl(k)=cdabs(vec(k))
 phase(k)=datan2( dimag(vec(k)),dreal(vec(k)) )
enddo
!$OMP END PARALLEL DO

dif=dabs(phase(1)-phase(2))
dtest=dabs(phase(1)-2*pi-phase(2))
if (dtest.le.dif) phase(1)=phase(1)-2*pi
dtest=dabs(phase(1)+2*pi-phase(2))
if (dtest.le.dif) phase(1)=phase(1)+2*pi

!do k=npol,2,-1

! dif=phase(k-1)-phase(k)
! if (dif.gt.+2.*dacos(-1.d0) ) phase(k)=phase(k)-2.*dacos(-1.d0)
! if (dif.lt.-2.*dacos(-1.d0) ) phase(k)=phase(k)+2.*dacos(-1.d0)
!
! !! Checking the sign
! if ((dif-dacos(-1.d0))**2.lt.dif**2) then
!  phase(k-1)=phase(k-1)-dacos(-1.d0)
!  ampl(k-1)=-ampl(k-1)
! endif
! ! Changing phase +pi
! if ((dif+dacos(-1.d0))**2.lt.dif**2) then
!  phase(k-1)=phase(k-1)+dacos(-1.d0)
!  ampl(k-1)=-ampl(k-1)
! endif
!
!enddo


return
end subroutine ampphase

subroutine evaluatespline(npot,datsp,time,time0,mat)
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: time,time0
complex*16, intent(in) :: datsp(4,npot,npot)
complex*16, intent(out) :: mat(npot,npot)

integer i,j,l

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l)
do i=1,npot
 do j=1,npot
  mat(i,j)=dcmplx(0.d0,0.d0)
  do l=1,4
   mat(i,j)=mat(i,j)+datsp(l,i,j)*(time-time0)**float(l-1)
  enddo
 enddo
enddo
!$OMP END PARALLEL DO

end subroutine evaluatespline

subroutine evaluatespline_ampphase(npot,datsp,time,time0,mat)
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: time,time0
complex*16, intent(in) :: datsp(4,npot,npot)
complex*16, intent(out) :: mat(npot,npot)

real*8 amp,phase
integer i,j,l

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l,amp,phase)
do i=1,npot
 do j=1,npot
  amp=0.0d0
  phase=0.d0
  do l=1,4
   amp=amp+dreal(datsp(l,i,j))*(time-time0)**float(l-1)
   phase=phase+dimag(datsp(l,i,j))*(time-time0)**float(l-1)
  enddo
  mat(i,j)=dcmplx(amp,0.d0)*cdexp( dcmplx(0.d0,phase) )
 enddo
enddo
!$OMP END PARALLEL DO

end subroutine evaluatespline_ampphase
