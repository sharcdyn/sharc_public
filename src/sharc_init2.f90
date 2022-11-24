program initial_conditions_2
!!! This program will get the info from the 0,
!!! and the initial conditions from the 1, divided by normal mode, to
!!! create the geoms.init and geoms.xyz files
implicit none
integer nmod,nat,ntrajs
character*2, allocatable :: sat(:)
real*8, allocatable :: coord0(:,:), coord(:,:),vel(:,:),noat(:),mass(:)
real*8, allocatable :: modno(:,:,:),massred(:)
real*8, allocatable :: x(:),p(:)
integer, allocatable :: mode(:)
type trajs
 integer ntrajs
 real*8 eig
 real*8, allocatable :: x(:),p(:),Ep(:)
end type trajs
type(trajs), allocatable :: traj(:)

real*8 Epot,Ekin,Ekin2,Eq

integer i,j,k
character*100000 fich
real*8 kk
logical ex

write(6,*) "Number of atoms and number of normal modes sampled"
read(5,*) nat,nmod
if (nmod.eq.0) write(6,*) "All normal modes with the trajs file will be used"
write(6,*) "Number of common trajectories (0 to check it)"
read(5,*) ntrajs
if (nmod.eq.0) then
 nmod=0
 do i=1,3*nat
  write(fich,"(A,I0,A)") "nm_",i,"-trajs.out"
  inquire(file=fich,exist=ex) 
  if (ex) nmod=nmod+1
 enddo
 allocate(mode(nmod)) 
 nmod=0
 do i=1,3*nat
  write(fich,"(A,I0,A)") "nm_",i,"-trajs.out"
  inquire(file=fich,exist=ex) 
  if (ex) then
   nmod=nmod+1
   mode(nmod)=i
  endif
 enddo
 write(6,"(A,1000(x,I0))") "Normal modes to use ",mode
else if (nmod.lt.0) then
 nmod=-nmod
 write(6,*) "The last ",nmod," normal modes will consider"
 allocate(mode(nmod))
 j=3*nat
 do i=1,nmod
  mode(i)=j
  j=j-1
 enddo
else
 allocate(mode(nmod))
 write(6,*) "List of normal modes"
 read(5,*) (mode(i),i=1,nmod)
endif

allocate(sat(nat),coord0(nat,3),vel(nat,3),coord(nat,3),mass(nat),noat(nat))
open(1,file="geom",status="old",iostat=i)
if (i.ne.0) stop "File geom is compulsory"
do i=1,nat
 read(1,*) sat(i),noat(i),(coord0(i,j),j=1,3),mass(i)
enddo
close(1)

allocate(modno(nmod,nat,3),massred(nmod))
allocate(traj(nmod))
Eq=0.
do i=1,nmod
 write(fich,"(A,I0,A)") "nm_",mode(i),"-info.xyz"
 open(1,file=fich,status="old",iostat=j)
 if (j.ne.0) then
  write(6,*) "Problems with file with info for nm ",trim(fich)
  stop "Info file does not exist"
 endif
 read(1,*) j
 if (nat.ne.j) then
  write(6,*) "Number of atoms is not correct in file ",trim(fich)
  stop "Number of atoms incorrect in info file"
 endif
 read(1,*) kk,kk,massred(i)
 do j=1,nat
  read(1,*) fich,(coord(j,k),k=1,3),(modno(i,j,k),k=1,3)
  do k=1,3
   coord(j,k)=coord(j,k)*1.889726
   modno(i,j,k)=modno(i,j,k)*1.889726
   if (abs(coord(j,k)-coord0(j,k)).ge.1e-5) stop "Different equilibrium distances"
  enddo
 enddo
 close(1)
 write(fich,"(A,I0,A)") "nm_",mode(i),"-trajs.out"
 open(1,file=fich,status="old",iostat=j)
 if (j.ne.0) then
  write(6,*) "Problems with file with trajs for nm ",trim(fich)
  stop "Traj file does not exist"
 endif
 read(1,*) fich,j,kk
 traj(i)%ntrajs=ntrajs
 if (ntrajs.eq.0) traj(i)%ntrajs=j
 traj(i)%eig=kk
 Eq=Eq+traj(i)%eig
 allocate(traj(i)%x(j),traj(i)%p(j),traj(i)%Ep(j))
 write(6,*) "Reading ",traj(i)%ntrajs
 do j=1,traj(i)%ntrajs
  read(1,*) traj(i)%x(j),traj(i)%p(j),traj(i)%Ep(j)
 enddo
 close(1)
enddo

do i=1,nmod
 if (i.eq.1) ntrajs=traj(i)%ntrajs
 if (traj(i)%ntrajs.le.ntrajs) ntrajs=traj(i)%ntrajs
enddo
write(6,*) "Calculating the common set of trajectories available ",ntrajs

open(1,file="geoms.init")
open(2,file="geoms.xyz")
j=0
Epot=0.
Ekin=0.
Ekin2=0.
vel=0.


write(1,*) nat
do i=1,nat
 write(1,900) sat(i),noat(i),mass(i)
enddo
write(fich,901) Eq," Traj ",j,Epot,Ekin,Ekin2
call writegeomvel(nat,sat,coord0,2,fich,vel)
write(fich,902) j," Traj ",Epot+Ekin,Ekin,Ekin2,Eq
call writegeomvel2(nat,sat,coord0,1,fich,vel)

allocate(x(nmod),p(nmod))
do i=1,ntrajs
 Epot=0.
 Ekin=0.
 Ekin2=0.
 do j=1,nmod
  x(j)=traj(j)%x(i)
  p(j)=traj(j)%p(i)
  Ekin=Ekin+0.5*p(j)**2/massred(j)
  Epot=Epot+traj(j)%Ep(i)
 enddo
 call creategeomvel(nat,nmod,coord0,modno,x,coord,vel,p,massred)
 do j=1,nat
  do k=1,3
   Ekin2=Ekin2+0.5*mass(j)*1822.888*vel(j,k)**2
  enddo
 enddo
 write(fich,901) Epot+Ekin," Traj ",i,Epot,Ekin,Ekin2
 call writegeomvel(nat,sat,coord,2,fich,vel)
 write(fich,902) i," Traj ",Epot+Ekin,Ekin,Ekin2
 call writegeomvel2(nat,sat,coord,1,fich,vel)
enddo
close(1)
close(2)
 
 
900 format(A2,6(x,E20.10e3))
901 format((F20.10,A,I0,10(x,F20.10)))
902 format(I0,A,(1000(x,F20.10)))

end program
