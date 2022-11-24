! program to create displaced geometries using the normal mode info
program sharc_init_prepgeom
implicit none
integer nat,ngeo
real*8, allocatable :: geom0(:,:)
character*2, allocatable :: sat(:)
character*100000 fich

real*8, allocatable :: geom(:,:),modno(:,:,:),freq(:),massred(:)
real*8, allocatable :: displ(:),x(:),dispfactor(:),fc(:)

integer nmod
integer, allocatable :: lmode(:)
real*8 Eharm

real*8 kk
integer i,j,k


write(6,*) "Number of modes to be used and xyz file to storage the geometries"
read(5,*) nmod,fich
allocate(lmode(nmod),x(nmod))
open(1,file=fich)
write(6,*) "List of normal modes to be used (if negative the movement will be relative to the harmonic width)"
read(5,*) (lmode(i),i=1,nmod)
do i=1,nmod
 write(6,*) "Reading data for normal mode ",lmode(i)
 
 write(fich,"(A,I0,A)") "nm_",iabs(lmode(i)),"-info.xyz"
 open(2,file=fich,status="old")
 read(2,*) j
 if (i.eq.1) then
  nat=j
  allocate(geom(nat,3),geom0(nat,3),sat(nat))
  allocate(modno(3*nat,nat,3),freq(3*nat),displ(3*nat),massred(3*nat),dispfactor(3*nat),fc(3*nat))
  modno=0.
  massred=0.
  freq=0.
  displ=0.
  fc=0.
  dispfactor=0.
 else
  if (j.ne.nat) stop "All normal modes must have the same number of atoms"
 endif
 read(2,*) j,freq(iabs(lmode(i))),massred(iabs(lmode(i)))
 do j=1,nat
  read(2,*) fich,(geom(j,k),k=1,3),(modno(iabs(lmode(i)),j,k),k=1,3)
  if (i.eq.1) then
   sat(j)=trim(fich)
  else
   if (trim(fich)==sat(j)) stop "Atoms must be the same in all the normal modes"
  endif
  do k=1,3
   geom(j,k)=geom(j,k)*1.889726
   modno(iabs(lmode(i)),j,k)=modno(iabs(lmode(i)),j,k)*1.889726
  enddo
 enddo
 if (i.eq.1) then
  geom0=geom
 endif
 close(2)
 fc(iabs(lmode(i)))=massred(iabs(lmode(i)))*freq(iabs(lmode(i)))**2
 dispfactor(iabs(lmode(i)))=1.d0
 if (lmode(i).le.0) then
  kk=massred(iabs(lmode(i)))*abs(freq(iabs(lmode(i))))
  dispfactor(iabs(lmode(i)))=1.d0/sqrt(kk)
  write(6,"(A,I0,2(X,F10.6))") "Applying conversion factor ",lmode(i),dispfactor(iabs(lmode(i))),fc(iabs(lmode(i)))
 endif
enddo
 
write(6,*) "Number of geometries (if negative an harmonic energy criteria will be used)"
read(5,*) ngeo
do i=1,iabs(ngeo)
 read(5,*) (x(j),j=1,nmod)
 if (ngeo.le.0) then
  do j=1,nmod
   if (lmode(j).le.0) x(j)=x(j)*freq(iabs(lmode(j)))
   kk=dsqrt( dabs(x(j))*2.d0/fc(iabs(lmode(j))) )
   if (x(j).le.0) kk=-kk
   x(j)=kk
  enddo
 endif
 kk=0.
 Eharm=0.
 do j=1,nmod
  displ(iabs(lmode(j)))=x(j)*dispfactor(iabs(lmode(j)))
  Eharm=Eharm+fc(iabs(lmode(j)))/2.d0*displ(iabs(lmode(j)))*displ(iabs(lmode(j)))
  kk=kk+displ(iabs(lmode(j)))*displ(iabs(lmode(j)))
 enddo
 kk=sqrt(kk)
 if (nmod.eq.1 .and. x(1).lt.0.d0) kk=-kk
 call creategeom(nat,3*nat,geom0,modno,displ,geom)
 write(fich,"(2(F10.6,X),A,I0,100000(X,I0,X,F10.6))") Eharm,kk,"Geometry ",i,(lmode(j),x(j),j=1,nmod)
 call writegeom(nat,sat,geom,1,fich)
enddo

end
 
