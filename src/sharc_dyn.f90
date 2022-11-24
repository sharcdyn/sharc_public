program QCTnD_xyz

use molecule_

implicit none
integer npot,nat
real*8 dt,time,tf
real*8 popt, Ekt, Ek0
character*2, allocatable :: sat(:)
real*8, allocatable :: x(:,:), vel(:,:), mass(:), noat(:), Ek(:,:)
real*8, allocatable :: deriv1(:,:,:),W(:)
real*8, allocatable :: Umini(:,:)

complex*16, allocatable :: coef(:),U(:,:),deriv2(:,:)
complex*16, allocatable :: U_prev(:,:), U_2prev(:,:)
complex*16, allocatable :: Uprop(:,:)

complex*16, allocatable :: dipole(:,:,:,:),Tpart(:,:,:),V(:,:,:)
complex*16, allocatable :: gradV(:,:,:,:)
logical, allocatable :: egradV(:,:)


integer pot
integer printstep,nstep
real*8 deco,rando(1), Epot
!real*8 rando
integer, allocatable :: seed(:)

integer nlas,nseed
integer npart
integer, allocatable :: tlas(:),npar(:),dlas(:)
real*8, allocatable :: laspar(:)
real*8 las(3), Emin

integer nsteps,npoint, npol
integer typdyn,typmodel,typsh

integer i,j
logical ex,res
integer printlevel

npoint=0
write(6,*) "Number of Atoms, How many potentials and the initial one"
read(5,*) nat, npot, pot
write(6,*) nat, npot,pot
if (pot.le.0) then
 write(6,*) "Potential impossible, let's try reading it in geometry"
endif


call input(typdyn,typmodel,typsh,npol,deco,printlevel,Emin,nsteps)


allocate ( mass(nat),x(nat,3),deriv1(npot,nat,3),vel(nat,3), sat(nat))
allocate ( noat(nat),Ek(nat,3),W(npot) )
allocate ( coef(npot),U(npot,npot),deriv2(npot,npot) )
allocate ( U_prev(npot,npot), U_2prev(npot,npot), Umini(npot,npot),Uprop(npot,npot) )
allocate ( dipole(npol,3,npot,npot),V(npol,npot,npot) ,Tpart(npol,npot,npot))
allocate ( gradV(npot,npot,nat,3),egradV(npot,npot) )


write(6,*) "Times (initial, step, end)"
read(5,*) time,dt,tf
write(6,"(3(x,F10.2),x,I5)") time,dt,tf
time=time*41.32
dt=dt*41.32
tf=tf*41.32
write(6,"(3(x,F10.2),x,I5)") time,dt,tf,nsteps
write(6,*) "Print step"
read(5,*) printstep

inquire(file="seeds.dat",exist=ex)
if (ex .eqv. .false.) then
 open(1,file="seeds.dat")
 write(6,*) "Generating seeds and storage them in seeds.dat"
 call random_seed(size=nseed)
 write(1,*) nseed
 allocate (seed(nseed))
 inquire(file="/dev/urandom",exist=res)
 if (res) then
  open(2, file='/dev/urandom', access='stream', form='UNFORMATTED')
  read(2) seed
  close(2)
 else
  call system_clock(count=j)
  seed = j + 37 * (/ (i - 1, i = 1, nseed) /)
 endif
 write(1,*) (seed(i),i=1,nseed)
 close(1)
 deallocate(seed)
endif
open(1,file="seeds.dat")
read(1,*) nseed
allocate(seed(nseed))
read(1,*) (seed(i),i=1,nseed)
close(1)
call random_seed(put=seed)
write(6,*) "Number of seeds ",nseed
write(6,*) "Seeds ",(seed(i),i=1,nseed)

inquire (file="restart.res",exist=res)

!   !!! dt.dat contains the time step
open(1,file="dt.dat")
write(1,"(A5,E20.10e3)") "dt = ",dt
close(1)

nstep=0
nlas=0
npart=0
las=0.d0
inquire(file="laser.dat",exist=ex)
if (ex) then
 open(1,file="laser.dat")
 write(6,*) "Reading laser parameters from laser.dat"
 read(1,*) nlas
 allocate(tlas(nlas),npar(nlas),dlas(nlas))
 do j=1,nlas
  read(1,*) tlas(j),dlas(j)
  select case (tlas(j))
   case(0)
    npar(j)=2
   case(1)
    npar(j)=6
   case(2)
    npar(j)=4
   case default
    write(6,*) "Laser type for the number ",j," does not exist"
    stop "This laser type does not exist"
  end select

  npart=npart+npar(j)
 enddo
 allocate( laspar(npart) )
 read(1,*) (laspar(i),i=1,npart)
 close(1)
else
 nlas=1
 allocate(tlas(1),npar(1),dlas(1),laspar(2))
 tlas(1)=0
 dlas(1)=1
 npar(1)=2
 laspar(:)=0.
 npart=npar(1)
endif
call laser(nlas,tlas,npar,dlas,npart,laspar,las,time)
write(6,*) "Laser ",las

 !call srand(seed(1))

if (res) then
 call restart(.false.,sat,nat,x,deriv1,deriv2,vel,npot,pot,coef,V,U,W,U_prev,dipole,las,Tpart, &
 Epot,time,npoint,mass,noat, Umini, npol)
! Getting previous random numbers
 do i=1,npoint*printstep
  do j=1,nsteps
   call random_number(rando)
   !rando(1)=rand(0)
  enddo
 enddo
else
 call readgeom(nat,npot,sat,noat,x,mass,vel,pot,coef)
endif
Ekt=0.d0
do i=1,nat
 do j=1,3
  Ek(i,j)=1./2.*mass(i)*vel(i,j)*vel(i,j)
  Ekt=Ekt+Ek(i,j)
 enddo
enddo
Ek0=Ekt


!   !!! End of initial conditions
write(6,*) "Initial condition read"
if (typmodel.eq.0) call modeldescription(nat,npot,sat)


!   !!! Output files
if (res) then
 open(11,file="traj.xyz",   access="append")
 open(12,file="coef.out",   access="append")
 open(67,file="coef_d.out", access="append")
 open(65,file="coef_raw.out",access="append")
 open(64,file="pop_raw.out", access="append")
 open(13,file="pop_d.out",  access="append")
 open(14,file="pop_a.out",  access="append")
 open(16,file="pop_d2.out", access="append")
 open(15,file="energy.out", access="append")
 open(31,file="dipole1.out",access="append")
 open(32,file="dipole2.out",access="append")
 open(33,file="dipole3.out",access="append")
 open(10,file="prob.out"   ,access="append")
 open(20,file="laser.out"  ,access="append")
 open(66,file="jumps.out"  ,access="append")
 open(40,file="nacme.out"  ,access="append")
 open(68,file="problems1.out", access="append")
 open(69,file="problems2.out", access="append")
 open(17,file="grads.out",access="append")
 open(18,file="grads0.out",access="append")

 write(6,*) "Restating"
else
 open(11,file="traj.xyz")
 open(12,file="coef.out")
 open(67,file="coef_d.out")
 open(65,file="coef_raw.out")
 open(64,file="pop_raw.out")
 open(13,file="pop_d.out")
 open(14,file="pop_a.out")
 open(16,file="pop_d2.out")
 open(15,file="energy.out")
 open(31,file="dipole1.out")
 open(32,file="dipole2.out")
 open(33,file="dipole3.out")
 open(10,file="prob.out")
 open(20,file="laser.out")
 open(66,file="jumps.out")
 open(40,file="nacme.out")
 open(68,file="problems1.out")
 open(69,file="problems2.out")
 open(17,file="grads.out")
 open(18,file="grads0.out")
 write(66,"(A)") "# TIME prevpot jumppot random accum.prob pot.prob total.prob Kinetic Potential.difference Pot.dif.end.of step."
 write(68,*)
 write(69,*)

 if (printlevel.ge.2) then
  
  open(100,file="ener_input.out")
  open(101,file="totham_input.out")
  open(102,file="ener_spl.out")
  open(103,file="totham_spl.out")
  open(106,file="ener_raw.out")
  open(107,file="totham_raw.out")
  open(108,file="ener_adia.out")
  open(109,file="totham_adia.out")

 endif
 if (printlevel.ge.3) then

  open(98,file="deb_inte")
  open(97,file="deb_Ttmp")
  open(96,file="deb_T")
  open(95,file="deb_Vtmp")
  open(94,file="deb_W")
  open(93,file="deb_Vad")
  open(92,file="deb_Vd")
  open(91,file="deb_kinetic")

  open(200,file="deb_eval-gradients")
  open(201,file="deb_Us")
  open(202,file="deb_eval-deriv1")

  open(666,file="deb_diag-phases")
  open(667,file="deb_diag-U2")
  open(999,file="deb_diag-Ignacio")

  open(99,file="deb_variado")

  open(104,file="deb_T2")
  open(105,file="deb_Umini")
 endif

 write(6,*) "before propag0",pot

!!!!!!!!!!!!!!!!
!!! Begining !!!
!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!
!!! First Step !!!
!!!!!!!!!!!!!!!!!!

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Before propag0"
close(1)
endif

 call propag0(printlevel,sat,nat,x,deriv1,deriv2,vel,npot,pot,coef, &
 V,U,W,U_prev,Uprop, U_2prev,dipole,las,typmodel,typdyn,typsh,Tpart, Epot, &
 Umini,npol, gradV,egradV,npoint,time)
 write(6,*) "after propag0",pot

 !write(96,900) time/41.32, ((dreal(Tpart(1,i,j)), dimag(Tpart(1,i,j)),j=i,npot),i=1,npot)
 !write(96,*) ""

!   !!! Calculating the norm
 popt=0.d0
 do i=1,npot
  popt=popt+cdabs(coef(i))**2
 enddo

 call geomwrit(sat,nat,time,x,vel,mass,npot,W,pot)
 if (printlevel .ge. 1) then
  do i=1,npot
   write(17,*) i,time/41.32
   call gradwrit(nat,deriv1(i,:,:),17)
  enddo
 endif
 call coeffwrit(time,npot,U,coef,pot,Uprop,deco)
 call properties(nat,time,npot,pot,Ek,W,dipole(1,:,:,:), Epot, deriv2,U)

 write(6,*) "Beginning"

 write(6,*) "Done ",pot,time/41.32,popt
endif

!!!!!!!!!!!!!!!!!!!
!!! Other Steps !!!
!!!!!!!!!!!!!!!!!!!

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Before next steps"
close(1)
endif

do while (time.lt.tf)
 nstep=nstep+1
 time=time+dt
 popt=0.d0
!   !!! Re-calculating the norm before the new step
 do i=1,npot
  popt=popt+cdabs(coef(i))**2
 enddo

 if (nstep.eq.printstep) then
  write(6,*) "Doing ",pot,time/41.32,popt
 endif
 
 npoint=npoint+1

! Writing data for restart
 call restart(.true.,sat,nat,x,deriv1,deriv2,vel,npot,pot,coef,V,U,W,U_prev,dipole,las,Tpart, &
 Epot,time,npoint,mass,noat, Umini, npol)

if (printlevel.ge.3) then
 open(1,file="debug.dat")
 write(1,*) "Before propag"
 write(1,*) "Printlevel ",printlevel
 write(1,*) "Sat ",sat
 write(1,*) "Nat ",nat
 write(1,*) "Time ",time
 write(1,*) "x ",x
 write(1,*) "vel ",vel
 write(1,*) "deriv1 ",deriv1
 write(1,*) "deriv2 ",deriv2
 write(1,*) "mass ",mass
 write(1,*) "dt ",dt
 write(1,*) "npot pot ",npot,pot
 write(1,*) "coef ",coef
 write(1,*) "V ",V
 write(1,*) "U ",U
 write(1,*) "W ",W
 write(1,*) "U prev ",U_prev
 write(1,*) "U prev2 ",U_2prev
 write(1,*) "dipole ",dipole
 write(1,*) "typdyn ",typdyn
 write(1,*) "typmodel ",typmodel
 write(1,*) "typsh ",typsh
 write(1,*) "deco ",deco
 write(1,*) "nsteps ",nsteps
 write(1,*) "Laser ",nlas,npart,npar,las,tlas,dlas,&
 laspar
 write(1,*) "T ",Tpart
 write(1,*) "Rando",rando
 write(1,*) "Epot Emin npol", Epot, Emin, npol
 write(1,*) "Noat ", noat
 write(1,*) "pritnstep ",printstep
 write(1,*) "Umini ", Umini
close(1)
endif


 call propag(printlevel,sat,nat,time,x,vel,deriv1,deriv2,mass,dt,npot,pot,coef,V,U,W,U_prev,Uprop,&
 U_2prev,dipole,typmodel,typdyn,typsh,deco,nsteps,nlas,npart,npar,las,tlas,dlas,&
 laspar,Tpart,rando, Epot, noat,printstep,nstep, Emin, Umini,npol, gradV,egradV,npoint)

 if (printlevel.ge.3) write(96,900) time/41.32, ((dreal(Tpart(1,i,j)), dimag(Tpart(1,i,j)),j=i,npot),i=1,npot)

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Propag finished"
close(1)
endif

 popt=0.d0
 do i=1,npot
  popt=popt+cdabs(coef(i))**2
 enddo

 if (nstep.eq.printstep) then
  write(6,*) "Done ",pot,time/41.32,popt
 endif
 !write(6,*) Ek(1,1), Ek(1,2), Ek(1,3), Ek(2,1), Ek(2,2), Ek(2,3)

 if (nstep.eq.printstep) then
  nstep=0

  Ekt=0.d0
  do i=1,3
   do j=1,nat
    Ek(j,i)=1./2.*mass(j)*vel(j,i)*vel(j,i)
    Ekt=Ekt+Ek(j,i)
   enddo
   !write(79,*) (vel(j,i),j=1,nat)
  enddo
  !write(78,*) (mass(j),j=1,nat)
  call geomwrit(sat,nat,time,x,vel,mass,npot,W,pot)

  if (printlevel.ge.1) then
   do i=1,npot
    write(17,*) i,time/41.32
    call gradwrit(nat,deriv1(i,:,:),17)
   enddo
  endif

  if (printlevel.eq.0) call coeffwrit(time,npot,U,coef,pot,Uprop,deco)
  call properties(nat,time,npot,pot,Ek,W,dipole(1,:,:,:), Epot,deriv2,U)

  write(20,900) time/41.32,las

!  write(77,900) time/41.32, ((dreal(U(i,j)),j=1,npot),i=1,npot)
!  write(78,900) time/41.32, ((dimag(U(i,j)),j=1,npot),i=1,npot)

 endif

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Before laser"
close(1)
endif

 call laser(nlas,tlas,npar,dlas,npart,laspar,las,time)

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Laser finished"
close(1)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Damped Dynamics Mode !!!
!!! Sort of ...          !!!
!!! Scaling velocities to mantain the Ek0 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (typdyn.lt.0) then
 if (printlevel.ge.3) write(91,*) "Damping"
 Ekt=0.
 do i=1,nat
  do j=1,3
   Ekt=Ekt+mass(i)*vel(i,j)**2
  enddo
 enddo
 do i=1,nat
  do j=1,3
   vel(i,j)=vel(i,j)*sqrt(Ek0/Ekt)
  enddo
 enddo
endif

enddo

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Program finished"
close(1)
endif
 
900 format (10000(x,E20.10e3))
end program




subroutine propag0(printlevel,sat,nat,x,deriv1,deriv2,vel,npot,pot,coef,V,U,&
W,U_prev,Uprop, U_2prev,dipole,las,typmodel,typdyn,typsh,Tpart, Epot, &
Umini, npol, gradV,egradV,npoint,time)
implicit none
integer nat,npot,pot,typdyn,typsh,typmodel, npol,npoint
real*8 x(nat,3),vel(nat,3), Epot
complex*16 V(npol,npot,npot),U(npot,npot),dipole(npol,3,npot,npot),Tpart(npol,npot,npot)
complex*16 U_prev(npot,npot), U_2prev(npot,npot),Uprop(npot,npot)
complex*16 coef(npot),coef0(npot)
real*8 las(3),W(npot),time
real*8 Umini(npot,npot)
character*2 sat(nat)

complex*16 gradV(npot,npot,nat,3)
logical egradV(npot,npot)

!   !!! Two derivatives:
!   !!! One of the diagonal elements
real*8 deriv1(npot,nat,3)
!   !!! Another for the diagonalization matrix (used in adiabatic)
complex*16 deriv2(npot,npot)
integer i,j,k,l

real*8 rnd,prob
integer printlevel
logical ex
complex*16 Vtmp(npot,npot),diptmp(3,npot,npot)

!interface
!
!subroutine diag(npot,V,dipole,las,U,typdyn,debug,reorder,diag_type,un,ignacios)
!implicit none
!integer, intent(in) ::npot,typdyn
!integer, intent(in) :: un
!complex*16, intent(in) :: V(npot,npot),dipole(3,npot,npot)
!complex*16 U(npot,npot)
!real*8, intent (in) :: las(3)
!real(kind(1d0)), intent(in) :: ignacios
!integer, intent(in) :: diag_type
!logical, intent(in) :: debug
!logical, intent(in) :: reorder
!end subroutine diag
!end interface


deriv2=dcmplx(0.d0,0.d0)
gradV=dcmplx(0.d0,0.d0)
egradV=.false.

!   !!! Creating an Identity "U" matrix (rotation matrix)
do i=1,npot
 do j=1,npot
  U(i,j)=dcmplx(0.d0,0.d0)
 enddo
 U(i,i)=dcmplx(1.d0,0.d0)
enddo
U_prev=U
Uprop=U

write(6,*) "calling evaluation"

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Before propag0/evaluate"
close(1)
endif

call stepwrit(npot,pot,npoint,U)
!!!! Just the potential part if analyitical gradients are available
if (typmodel.eq.1) then
 call stepwrit(npot,pot,npoint,U)
 call singlepoint(printlevel,nat,npot,x,sat,V(1,:,:),dipole(1,:,:,:),Tpart(1,:,:),vel,U,las,typdyn)
 egradV=.false.
 call stepwrit(npot,pot,npoint,U)
endif

call evaluate (printlevel, sat,nat,npot,x,V(1,:,:),U,W,deriv1,deriv2, dipole(1,:,:,:), &
.true.,las,typdyn,typmodel,vel,Tpart(1,:,:), gradV,egradV,pot )
!!! Writting gradients
if (printlevel.ge.1) then
 do i=1,npot
  do j=1,npot
   write(18,*) "POT ",i,j,egradV(i,j),V(1,i,j)
   if (egradV(i,j)) call gradwritc(nat,gradV(i,j,:,:),18)
  enddo
 enddo
endif

if (printlevel.ge.2) then
 write(106,"(1000000(E20.10e3))") time,(dreal(V(1,i,i)),i=1,npot)
 write(107,"(1000000(E20.10e3))") time,((V(1,i,j),j=1,npot),i=1,npot)
endif


  
do i=1,npot
 do j=1,npot
  Vtmp(i,j)=V(1,i,j)
  do k=1,3
   diptmp(k,i,j)=dipole(1,k,i,j)
  enddo
!   !!! Setting V, dipole, T for all polynomials equal to the first value for each one
  do l=1,npol
   V(l,i,j)=Vtmp(i,j)
   Tpart(l,i,j)=dcmplx(0.d0,0.d0)
   do k=1,3
    dipole(l,k,i,j)=diptmp(k,i,j)
   enddo
  enddo
 enddo
enddo

!do i=1,npot
! do j=1,nat
!  write(6,*) deriv1(i,j,:)
! enddo
!enddo



if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "After propag0/evaluate"
close(1)
endif

call diago(npot,Vtmp,diptmp,las,U,typdyn,printlevel,.false.,1,69,.false.)

write(6,*) "Back evaluate"

U_prev=U
U_2prev=U

!do i=1,npot
! coef(i)=dcmplx(0.d0,0.d0)
! do j=1,npot
!  coef(i)=coef(i)+dconjg(U(j,i))*coef0(j)
! enddo
!enddo

coef0=coef



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! "Diabatic" potential beginning !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   !!! Redoing the coefficient by collapsing to the most feasible
!   !!! in general at the begining only one potential should be important
!   !!! This is done if the file diabatic exists

inquire(file="diabatic",exist=ex)
if (ex) then
 write(6,*) "Considering a diabatic initial state and recalculating the adiabatic one"
 call adiabat_coef(npot,U,coef0)


 coef=coef0
!!! Redefining pot by collapsing the wavefunction with probabilities and random numbers
 call random_number(rnd)
 prob=0.
 do i=1,npot
  prob=prob+cdabs(coef(i))**2
 enddo
 write(6,*) "Total probability in the adiabatic coefficients ",prob
 do i=1,npot
  coef(i)=coef(i)/sqrt(prob)
 enddo
 prob=0.
 pot=0
 do while (prob.le.rnd)
  pot=pot+1
  prob=prob+cdabs(coef(pot))**2
 enddo
 write(6,*) "Accumulated probability ",prob
 write(6,*) "State probability ",cdabs(coef(pot))**2,coef(pot)
 write(6,*) "Random number ",rnd
endif
write(6,*) "Propagating in potential ",pot


write(6,*) "Adiabatic"
do i=1,npot
 write(6,"(I6,X,F6.4,2(X,E20.12e3))") i, cdabs(coef(i))**2,dreal(coef(i)), dimag(coef(i))
enddo

write(6,*) "Original"
call diabat_coef(npot,U,coef0)
do i=1,npot
 write(6,"(I6,X,F6.4,2(X,E20.12e3))") i,cdabs(coef0(i))**2,dreal(coef0(i)), dimag(coef0(i))
enddo

!   !!! Using the gradients of Potentials and Dipoles, only for Ehrenfest
if (typsh.eq.0) then
 call ehrenfest (coef, nat,npot, deriv1, V, Epot,deriv1(pot,:,:),npol )
else
 Epot=W(pot)
endif

Umini=0.d0
do i=1,npot
 Umini(i,i)=+1.0
enddo


return
end subroutine propag0


subroutine ehrenfest (coef, nat,npot, deriv1, V, Ehrpot, Ehrgrad, npol)
implicit none
integer npot,nat,npol
real*8 deriv1(npot,nat,3), Ehrgrad(nat,3), Ehrpot
complex*16 coef(npot),V(npol,npot,npot), kk
integer i,j,k
real*8 norm


Ehrgrad=0.d0
Ehrpot=0.d0
norm=0.d0
do i=1,npot
 Ehrpot=Ehrpot+ dreal( cdabs(coef(i)**2)*V(1,i,i))
 norm=norm+cdabs(coef(i))**2
enddo
!write(6,*) norm

!   !!! WARNING Ehrenfest without NAC vectors!!!!!
do j=1,nat
 do k=1,3
  kk=dcmplx(0.d0,0.d0)
  do i=1,npot
   kk= kk+ (deriv1(i,j,k)*(cdabs(coef(i))**2))
  enddo
  Ehrgrad(j,k)=dreal(kk)/norm
 enddo
enddo

!write(6,*) Ehrpot
return
end subroutine


subroutine restart(writ,sat,nat,geom,deriv1,deriv2,vel,npot,pot,coef,V,U,W,U_prev,dipole,las,Tpart, &
Epot,time,npoint,mass,noat,Umini,npol)
implicit none
logical writ
integer nat,npot,pot,npoint, npol
character*2 sat(nat)
real*8 mass(nat),noat(nat),geom(nat,3),vel(nat,3),time
complex*16 coef(npot),V(npol,npot,npot),deriv2(npot,npot)
complex*16 U(npot,npot),U_prev(npot,npot),dipole(npol,3,npot,npot),Tpart(npol,npot,npot)
real*8 deriv1(npot,nat,3)
real*8 W(npot),las(3),Epot


integer i,j,k,ii
real*8 Umini(npot,npot)
logical reading

reading=.true.
if (writ) reading=.false.

if (reading) open(1,file="restart.res",form="unformatted")
if (writ)    open(1,file="restart",form="unformatted")
 ! Parameters
if (writ) write(1) nat,npot,pot,npoint,time,Epot,npol
if (reading) read(1) nat,npot,pot,npoint,time,Epot,npol

 ! Geoemtry and velocities
do i=1,nat
 if (writ) write(1) sat(i),noat(i),(geom(i,j),j=1,3), mass(i), (vel(i,j),j=1,3)
 if (reading) read(1) sat(i),noat(i),(geom(i,j),j=1,3), mass(i), (vel(i,j),j=1,3)
! write(6,*) sat(i),noat(i),(geom(i,j),j=1,3)
! write(6,*) mass(i), (vel(i,j),j=1,3)
enddo

! Coeficients
if (writ) write(1) (coef(i),i=1,npot)
if (reading) read(1) (coef(i),i=1,npot)

!Derivatives
do i=1,npot
 do j=1,nat
  if (writ) write(1) (deriv1(i,j,k),k=1,3)
  if (reading) read(1) (deriv1(i,j,k),k=1,3)
 enddo
enddo

!Hamiltonian matrix
do ii=1,npol
 do i=1,npot
  if (writ) write(1) (V(ii,i,j),j=1,npot)
  if (reading) read(1) (V(ii,i,j),j=1,npot)
 enddo
! Kinetic coupling
 do i=1,npot
  if (writ) write(1) (Tpart(ii,i,j),j=1,npot)
  if (reading) read(1) (Tpart(ii,i,j),j=1,npot)
 enddo
!Dipoles
 do i=1,3
  do j=1,npot
   if (writ) write(1) (dipole(ii,i,j,k),k=1,npot)
   if (reading) read(1) (dipole(ii,i,j,k),k=1,npot)
  enddo
 enddo
enddo

!Potentials
if (writ) write(1) (W(i),i=1,npot)
if (reading) read(1) (W(i),i=1,npot)

!Laser
if (writ) write(1) (las(i),i=1,3)
if (reading) read(1) (las(i),i=1,3)

! U matrix
do i=1,npot
 if (writ)   write(1) (U(i,j),j=1,npot)
 if (reading) read(1) (U(i,j),j=1,npot)
enddo
! U_prev
do i=1,npot
 if (writ)   write(1) (U_prev(i,j),j=1,npot)
 if (reading) read(1) (U_prev(i,j),j=1,npot)
enddo

! Umini
do i=1,npot
 if (writ)   write(1) (Umini(i,j),j=1,npot)
 if (reading) read(1) (Umini(i,j),j=1,npot)
enddo

end subroutine


subroutine velver1(nat,x,vel,deriv,mass,dt,printlevel,noat,time)
implicit none
integer,intent(in) :: printlevel
integer, intent(in) :: nat
real*8, intent(in) :: mass(nat),dt,noat(nat),time
real*8,intent(in) :: deriv(nat,3)
real*8,intent(inout) :: x(nat,3),vel(nat,3)

real*8 acel(nat,3),Ekt
integer i,j,l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Velocity Verlet algorithm !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   !!! First part of the algorithm
if (printlevel.ge.3) then
write(91,*) "Time ",time
write(91,*) "before verlet "
write(91,*)"X"
write(91,*) (x(:,j),j=1,3)
write(91,*)"vel"
write(91,*) (vel(:,j),j=1,3)
Ekt=0.d0
do j=1,3
 do l=1,nat
  Ekt=Ekt+mass(l)*vel(l,j)**2/2.
 enddo
enddo
write(91,*) "Kinetic ",Ekt
endif


do j=1,nat
 do i=1,3
  !   !!! Acceleration at t
  if (int(noat(j)).ne.0.and.mass(j).gt.0) acel(j,i)=-1./mass(j)*deriv(j,i)
  !   !!! Position at t+dt
  x(j,i)=x(j,i)+vel(j,i)*dt+1./2.*acel(j,i)*dt*dt
  !   !!! Velocity at t+(dt/2)
  vel(j,i)=vel(j,i)+acel(j,i)*dt/2.
 enddo
enddo

if (printlevel.ge.3) then
write(91,*) "acel"
write(91,*) (acel(:,j),j=1,3)
write(91,*) "after first verlet "
write(91,*)"X"
write(91,*) (x(:,j),j=1,3)
write(91,*)"vel"
write(91,*) (vel(:,j),j=1,3)
Ekt=0.d0
do j=1,3
 do l=1,nat
  Ekt=Ekt+mass(l)*vel(l,j)**2/2.
 enddo
enddo
write(91,*) "Kinetic ",Ekt
endif

end

subroutine velver2(nat,x,vel,deriv,mass,dt,printlevel,noat,time)
implicit none
integer,intent(in) :: printlevel
integer, intent(in) :: nat
real*8, intent(in) :: mass(nat),dt,noat(nat),time,deriv(nat,3)
real*8 x(nat,3),vel(nat,3)

real*8 acel(nat,3),Ekt
integer i,j,l


if (printlevel.ge.3) then
write(91,*) "before second verlet "
write(91,*)"X"
write(91,*) (x(:,j),j=1,3)
write(91,*)"vel"
write(91,*) (vel(:,j),j=1,3)
Ekt=0.d0
do j=1,3
 do l=1,nat
  Ekt=Ekt+mass(l)*vel(l,j)**2/2.
 enddo
enddo
write(91,*) "Kinetic ",Ekt
endif

!   !!! Second part of the VelVerlet algorithm
do j=1,nat
 do i=1,3
!   !!! acceleration at t+dt
  if (int(noat(j)).ne.0.and.mass(j).ge.1) acel(j,i)=-1./mass(j)*deriv(j,i)
!   !!! Velocity at t+dt
  vel(j,i)=vel(j,i)+acel(j,i)*dt/2.
! write(6,*) acel(j,i),1./2.*acel(j,i)*dt*dt
 enddo
enddo

if (printlevel.ge.3) then
write(91,*) "acel"
write(91,*) (acel(:,j),j=1,3)
write(91,*) "after second verlet "
write(91,*) "Time ",time
write(91,*)"X"
write(91,*) (x(:,j),j=1,3)
write(91,*)"vel"
write(91,*) (vel(:,j),j=1,3)
Ekt=0.d0
do j=1,3
 do l=1,nat
  Ekt=Ekt+mass(l)*vel(l,j)**2/2.
 enddo
enddo
write(91,*) "Kinetic ",Ekt
endif

end


