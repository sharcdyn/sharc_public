subroutine laser(nlas,tlas,npar,dlas,npart,laspar,las,time)
implicit none
integer nlas,npart
integer tlas(nlas),npar(nlas),dlas(nlas)
real*8 laspar(npart),las(3),time

integer i,j
real*8 kk

las=0.d0
j=0
do i=1,nlas

 if (tlas(i).eq.0) then
 ! cosine laser
 ! par 1 frequency
 ! par 2 amplitude
  kk=laspar(j+2)*dcos(laspar(j+1)*time)
 endif

 if (tlas(i).eq.1) then
! sin**2 cos**2 laser field
! par 1 frequency
! par 2 max amplitude
! par 3 start laser
! par 4 start of platoo
! par 5 end of platoo
! par 6 end of laser
  kk=0.
  if ( time.ge.laspar(j+3)*41.32.and.time.le.laspar(j+4)*41.32 ) then
   kk=time-laspar(j+3)*41.32
   kk=kk/(laspar(j+4)*41.32-laspar(j+3)*41.32)
   kk=kk*dacos(-1.d0)/2.
   kk=dsin(kk)**2
  endif
  if ( time.ge.laspar(j+4)*41.32.and.time.le.laspar(j+5)*41.32 ) kk=1.d0
  if ( time.ge.laspar(j+5)*41.32.and.time.le.laspar(j+6)*41.32 ) then
   kk=time-laspar(j+5)*41.32
   kk=kk/(laspar(j+6)*41.32-laspar(j+5)*41.32)
   kk=kk*dacos(-1.d0)/2.
   kk=dcos(kk)**2
  endif
  kk=kk*laspar(j+2)*dcos(laspar(j+1)*time)
 endif
 if (tlas(i).eq.2) then
 ! cosh**-2
 ! par 1 frequency
 ! par 2 max amplitude
 ! par 3 time of maximum
 ! par 4 width
  kk=(time-laspar(j+3)*41.32)/(laspar(j+4)*41.32)
  kk=1.d0/dcosh(kk)/dcosh(kk)
  kk=kk*laspar(j+2)*dcos(laspar(j+1)*time)
 endif
 las(dlas(i))=las(dlas(i))+kk
 j=j+npar(i)
enddo
end subroutine laser
