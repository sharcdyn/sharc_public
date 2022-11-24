program internal
implicit none
integer nat
character*2, allocatable :: at(:)
real*8, allocatable :: geom(:,:), mass(:),vel(:,:)
character*10000 kkchar

!!!! Center of mass
integer ncm
logical, allocatable :: atcm(:,:)
type cmat
 integer natcm
 integer, allocatable :: atcm(:)
 real*8 mass
end type
type(cmat), allocatable :: cm(:)

!!!! Distances
integer ndist
real*8, allocatable :: dist0(:),dist(:)
integer, allocatable :: idist(:,:)
!!!! Angles
integer nang
real*8, allocatable :: ang0(:),ang(:)
integer, allocatable :: iang(:,:)
!!!! Dihedral
integer ndie
real*8, allocatable :: die0(:),die(:)
integer, allocatable :: idie(:,:)

logical ex

character*20 fich
integer ngeom
integer i,j,k

inquire (file="traj.xyz",exist=ex)
if (ex) then
 write(6,*) "Reading file traj.xyz"
 open(1,file="traj.xyz")
 read(1,*) nat
 rewind(1)
else
 write(6,*) "traj.xyz does not exist. Stopping"
 stop "traj.xyz does not exist. Stopping"
endif

write(6,*) "Number of center of masses"
read(5,*) ncm
if (ncm.ne.0) then
 allocate(cm(ncm),atcm(nat,ncm) )

 do i=1,ncm
  write(fich,"(A,I3.3,A)") "fragment",i,".out"
  open(100+i,file=fich)
  write(100+i,"(A)") "# Ekin, Et, Erot, dist (A), Etr, ..."
 enddo

 inquire(file="geom",exist=ex)
 if (ex) then
  write(6,*) "Masses will be read from the geom file"
 else
  stop "File geom does not exist and it is required to read masses"
 endif


 do i=1,ncm
  write(6,*) "Number of atoms in the center of mass ",i
  read(5,*) cm(i)%natcm
  allocate( cm(i)%atcm(cm(i)%natcm) )
  write(6,*) "Which atoms belongs to that center"
  read(5,*) (cm(i)%atcm(j),j=1,cm(i)%natcm)
  do j=1,nat
   atcm(j,i)=.false.
  enddo
  do j=1,cm(i)%natcm
   atcm(cm(i)%atcm(j),i)=.true.
  enddo
  write(6,*) "This center of mass will correspond with atom ",nat+i
 enddo
endif

write(6,*) "Number of distances"
read(5,*) ndist
allocate( dist(ndist),idist(ndist,2),dist0(ndist) )
do i=1,ndist
 write(6,*) "Distance ",i
 write(6,*) "Atoms in the distance "
 read(5,*) idist(i,:)
enddo

write(6,*) "Number of angles"
read(5,*) nang
allocate( ang(nang),iang(nang,3),ang0(nang) )
do i=1,nang
 write(6,*) "Angle ",i
 write(6,*) "Atoms in the angle "
 read(5,*) iang(i,:)
enddo

write(6,*) "Number of dihedral angles"
read(5,*) ndie
allocate( die(ndie),idie(ndie,4),die0(ndie) )
do i=1,ndie
 write(6,*) "Dihedral angle ",i
 write(6,*) "Atoms in the dihedral "
 read(5,*) idie(i,:)
enddo

open(10,file="distance.time")
open(11,file="angle.time")
open(12,file="dihedral.time")
dist=0.d0
ngeom=0
do while (.true.)
 read(1,*,END=10) nat
 ngeom=ngeom+1
 if (ngeom.eq.1) then
  if (ncm.ne.0) then
   allocate(mass(nat) )
   write(6,*) "Reading masses from geom"
   inquire (file="geom",exist=ex)
   if (ex) then
    open(3,file="geom")
    do i=1,nat
     read(3,*) (fich,j=1,5),mass(i)
    enddo
    close(3)
   else
    write(6,*) "geom file does not exist, it is required to create a center of mass"
    stop  "geom file does not exist, it is required to create a center of mass"
   endif

   write(6,*) "xyz out file, containing the center of masses"
   read(5,*) fich
   open(2,file=fich)

   do i=1,ncm
    write(6,*) "Defining atom ",i+nat
    write(6,*) (cm(i)%atcm(:))
    cm(i)%mass=0.d0
    do j=1,cm(i)%natcm
     cm(i)%mass=cm(i)%mass+mass(cm(i)%atcm(j))
    enddo
    write(6,*) "Mass ",cm(i)%mass
   enddo
  endif
  allocate( geom(nat+ncm,3),at(nat),vel(nat+ncm,3) )
 endif

 read(1,*) fich
 do i=1,nat
  read(1,"(A)") kkchar
  kkchar=trim(kkchar)//" 0. 0. 0."
  read(kkchar,*) at(i),geom(i,:),vel(i,:)
 enddo
 do i=1,ncm
  do j=1,3
   geom(i+nat,j)=0.d0
   vel(i+nat,j)=0.d0
   do k=1,cm(i)%natcm

!    write(6,*) "Accumulating ",i,j,k,cm(i)%atcm(k)
    geom(i+nat,j)=geom(i+nat,j)+mass(cm(i)%atcm(k))*geom(cm(i)%atcm(k),j)

!!! Total momenta
    vel(i+nat,j)=vel(i+nat,j)+mass(cm(i)%atcm(k))*vel(cm(i)%atcm(k),j)
!    write(6,*) geom(i+nat,j),vel(i+nat,j)
   enddo
   geom(i+nat,j)=geom(i+nat,j)/cm(i)%mass
   vel(i+nat,j)=vel(i+nat,j)/cm(i)%mass
!   write(6,*) "Final ",geom(i+nat,j),vel(i+nat,j)
  enddo
 enddo

 if (ncm.ne.0) then
  write(2,*) nat+ncm
  write(2,"(A)") trim(fich)
  do i=1,nat
   write(2,"(A2,6(x,E20.10e3))") at(i),geom(i,:),vel(i,:)
  enddo
  do i=1,ncm
   write(2,"(A2,6(x,E20.10e3))") "XX",geom(i+nat,:),vel(i+nat,:)
  enddo
 endif

 do i=1,ndist
  call distance(idist(i,:),nat+ncm,geom,dist0(i))
  dist(i)=dist(i)+dist0(i)
 enddo
 write(10,900) dist0(:)

 do i=1,nang
  call angle(iang(i,:),nat+ncm,geom,ang0(i))
  ang(i)=ang(i)+ang0(i)
 enddo
 write(11,900) ang0(:)

 do i=1,ndie
  call dihedral(idie(i,:),nat+ncm,geom,die0(i))
  die(i)=die(i)+die0(i)
 enddo
 write(12,900) die0(:)

 if (ncm.ne.0) call fragment(nat,ncm,geom,vel,atcm,mass)


900 format(100000(x,F20.10))
enddo

10 continue

write(6,*) "Average distances"
do i=1,ndist
 write(6,901) "Distance ",i,idist(i,:),dist(i)/float(ngeom)
901 format(a9,3(I4.4,x),F20.10)
enddo

write(6,*) "Average angles"
do i=1,nang
 write(6,902) "Angle ",i,iang(i,:),ang(i)/float(ngeom)*180.d0/dacos(-1.d0)
902 format(a9,4(I4.4,x),F20.10)
enddo

write(6,*) "Average dihedral angles"
do i=1,ndie
 write(6,903) "Dihedral ",i,idie(i,:),die(i)/float(ngeom)*180.d0/dacos(-1.d0)
903 format(a9,5(I4.4,x),F20.10)
enddo

end

subroutine distance(atoms,nat,geom,r)
implicit none
integer, intent (in) :: atoms(2),nat
real*8, intent(in) :: geom(nat,3)
real*8, intent(out) :: r

integer i,j

r=0.d0
do i=1,3
 r=r+(geom(atoms(1),i)-geom(atoms(2),i))**2
enddo
r=dsqrt(r)

return
end


subroutine angle(atoms,nat,geom,r)
implicit none
integer, intent (in) :: atoms(3),nat
real*8, intent(in) :: geom(nat,3)
real*8, intent(out) :: r

integer i,j
real*8 vec1(3),vec2(3),r1,r2

r1=0.d0
r2=0.d0
do i=1,3
 vec1(i)=geom(atoms(1),i)-geom(atoms(2),i)
 vec2(i)=geom(atoms(3),i)-geom(atoms(2),i)
 r1=r1+vec1(i)**2
 r2=r2+vec2(i)**2
enddo
r1=dsqrt(r1)
r2=dsqrt(r2)

r=0.d0
if (r1.gt.0.and.r2.gt.0) then
 do i=1,3
  r=r+vec1(i)*vec2(i)
 enddo
 r=r/r1/r2
 r=dacos(r)
endif

return
end


subroutine dihedral(atoms,nat,geom,r)
implicit none
integer, intent (in) :: atoms(4),nat
real*8, intent(in) :: geom(nat,3)
real*8, intent(out) :: r

integer i,j
real*8 b(3,3),r0(3)
real*8 b12(3),b23(3),r12,r23
real*8 r1,r2

r0=0.d0
do i=1,3
 b(1,i)=geom(atoms(2),i)-geom(atoms(1),i)
 b(2,i)=geom(atoms(3),i)-geom(atoms(2),i)
 b(3,i)=geom(atoms(4),i)-geom(atoms(3),i)
 do j=1,3
  r0(j)=r0(j)+b(j,i)**2
 enddo
enddo
do j=1,3
 r0(j)=dsqrt(r0(j))
enddo

b12(1)=+(b(1,2)*b(2,3)-b(1,3)*b(2,2))
b12(2)=-(b(1,1)*b(2,3)-b(1,3)*b(2,1))
b12(3)=+(b(1,1)*b(2,2)-b(1,2)*b(2,1))
r12=dsqrt(b12(1)**2+b12(2)**2+b12(3)**2)

b23(1)=+(b(2,2)*b(3,3)-b(2,3)*b(3,2))
b23(2)=-(b(2,1)*b(3,3)-b(2,3)*b(3,1))
b23(3)=+(b(2,1)*b(3,2)-b(2,2)*b(3,1))
r23=dsqrt(b23(1)**2+b23(2)**2+b23(3)**2)

r1=0.d0
r2=0.d0
do i=1,3
 r1=r1+b(1,i)*b23(i)
 r2=r2+b12(i)*b23(i)
enddo
r1=r0(2)*r1
r=datan2(r1,r2)
!write(6,*) r1,r2,r

return
end

subroutine fragment(nat,ncm,geom,vel,atcm,mass)
implicit none
integer, intent(in) :: nat,ncm
logical, intent(in) :: atcm(nat,ncm)
real*8, intent(in) :: mass(nat)
real*8, intent(in) :: geom(nat+ncm,3),vel(nat+ncm,3)

real*8 x(nat,3),v(nat,3),xrel(3),Lfrag(3)

real*8 mom(ncm,3),cm(ncm,3)
real*8 Ekin(ncm),masscm(ncm)

real*8 inert(3,3),dinert(3),mat(3,3)

real*8 Et,Erot
real*8 dist(ncm),Etr(ncm)
real*8 vec(3)
integer i,j,k,l
real*8 pmom


do i=1,nat
 do j=1,3
  x(i,j)=geom(i,j)/.5292
  v(i,j)=vel(i,j)/1000.
 enddo
enddo

!! Properties of the center of masses
do i=1,ncm
 masscm(i)=0.d0
 Ekin(i)=0.d0
 do k=1,3
  mom(i,k)=0.d0
  cm(i,k)=0.d0
 enddo
 do j=1,nat
  if (atcm(j,i)) then
   masscm(i)=masscm(i)+mass(j)
   do k=1,3
    cm(i,k)=cm(i,k)+mass(j)*x(j,k)
    mom(i,k)=mom(i,k)+mass(j)*v(j,k)
    Ekin(i)=Ekin(i)+mass(j)*v(j,k)*v(j,k)
   enddo
  endif
 enddo
 masscm(i)=masscm(i)*1822.888
 Ekin(i)=Ekin(i)*1822.888*0.5
 do k=1,3
  cm(i,k)=cm(i,k)*1822.888/masscm(i)
  mom(i,k)=mom(i,k)*1822.888
 enddo
enddo

do i=1,ncm

!!! Translational energy
 Et=0.d0
 do k=1,3
  Et=Et+mom(i,k)*mom(i,k)
 enddo
 Et=Et*0.5/masscm(i)

!!! Rotational energy!!!!
 do k=1,3
  Lfrag(k)=0.d0
  do l=1,3
   inert(k,l)=0.d0
  enddo
 enddo
 do j=1,nat
  if (atcm(j,i)) then
!!! Relative momentum
   do k=1,3
    vec(k)=mass(j)*1822.888*( v(j,k) -mom(i,k)/masscm(i) )
   enddo
!!! Relative position
   do k=1,3
    xrel(k)=x(j,k)-cm(i,k)
   enddo

!!! Angular momentum
   Lfrag(1)=Lfrag(1)+  xrel(2)*vec(3)-xrel(3)*vec(2)
   Lfrag(2)=Lfrag(2)+  xrel(3)*vec(1)-xrel(1)*vec(3)
   Lfrag(3)=Lfrag(3)+  xrel(1)*vec(2)-xrel(2)*vec(1)

!!! Inertia momentum
   inert(1,1)=inert(1,1) + mass(j)*1822.888*(xrel(2)*xrel(2)+xrel(3)*xrel(3))
   inert(1,2)=inert(1,2) - mass(j)*1822.888* xrel(1)*xrel(2)
   inert(1,3)=inert(1,3) - mass(j)*1822.888* xrel(1)*xrel(3)
   inert(2,1)=inert(2,1) - mass(j)*1822.888* xrel(2)*xrel(1)
   inert(2,2)=inert(2,2) + mass(j)*1822.888*(xrel(1)*xrel(1)+xrel(3)*xrel(3))
   inert(2,3)=inert(2,3) - mass(j)*1822.888* xrel(2)*xrel(3)
   inert(3,1)=inert(3,1) - mass(j)*1822.888* xrel(3)*xrel(1)
   inert(3,2)=inert(3,2) - mass(j)*1822.888* xrel(3)*xrel(2)
   inert(3,3)=inert(3,3) + mass(j)*1822.888*(xrel(1)*xrel(1)+xrel(2)*xrel(2))
  endif
 enddo

 call diag(inert,3,dinert)
! call diag_order(inert,3,dinert)
!!! Calculating the inverse I'
 mat=0.
 do j=1,3
  if (dinert(j).le.1e-8) then
   mat(j,j)=0.
   do k=1,3
    inert(k,j)=0.
   enddo
  else
   mat(j,j)=1./dinert(j)
  endif
 enddo
!!!! Calculating I'=U Id-1 Udag
 mat=matmul(mat,transpose(inert))
 mat=matmul(inert,mat)

!!!! Calculating L\dag I' L
 vec=matmul(mat,Lfrag)
 Erot=vec(1)*Lfrag(1)+vec(2)*Lfrag(2)+vec(3)*Lfrag(3)
 Erot=Erot*.5

!!! Properties with respect to the other centers of mass
 do j=1,ncm
  dist(j)=0.d0
  do k=1,3
   vec(k)=cm(j,k)-cm(i,k)
   dist(j)=dist(j)+vec(k)**2
  enddo
  dist(j)=sqrt(dist(j))
  if (i.ne.j) then
   do k=1,3
    vec(k)=vec(k)/dist(j)
   enddo
  endif
!!! Relative translational energy
  pmom=0
  do k=1,3
!!! Momentum projected between the two center of masses
   pmom=pmom+mom(i,k)*vec(k)
  enddo
  Etr(j)=pmom**2/masscm(i)*0.5
  if (pmom.lt.0) Etr(j)=-Etr(j)
 enddo
 write(100+i,900) Ekin(i),Et,Erot,(dist(j),Etr(j),j=1,ncm)
enddo


900 format(10000000(x,E20.10e3))
end subroutine fragment
