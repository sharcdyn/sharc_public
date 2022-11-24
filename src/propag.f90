
subroutine propag(printlevel,sat,nat,time,x,vel,deriv1,deriv2,mass,dt,npot,pot,coef,V,U,W,U_prev,Uprop,&
U_2prev,dipole,typmodel,typdyn,typsh, deco,nsteps,nlas,npart,npar,las,tlas,dlas,&
laspar, T,rando, Epot, noat,printstep,nstep, Emin , Umini, npol, gradV,egradV, npoint)
implicit none
integer npot,pot,typdyn, nat, typmodel,typsh,npoint
real*8 x(nat,3),vel(nat,3),mass(nat), W(npot),noat(nat)
integer npol
real*8 dt,time
complex*16 V(npol,npot,npot)
complex*16 U(npot,npot),U_prev(npot,npot),U_2prev(npot,npot)
complex*16 T(npol,npot,npot)
complex*16 dipole(npol,3,npot,npot)
complex*16 coef(npot),coef_prev(npot), coefa(npot)
real*8 popt, Ek(nat,3),tau, suma,Etot
character*2 sat(nat)
real*8 deco, Epot, Emin
complex*16 gradV(npot,npot,nat,3)
logical egradV(npot,npot)

integer kkpot, potinic

!Raw data transformation
complex*16 Uprop(npot,npot)
real*8 Umini(npot,npot)

!Runge Kutta parameters
complex*16 inte1(npot,npot), inte2(npot,npot), inte3(npot, npot), intea(npot,npot)

!Laser parameters
real*8 las(3)
integer nlas
integer npart
integer tlas(nlas),npar(nlas),dlas(nlas)
real*8 laspar(npart)

!   !!! Two derivatives:
!   !!! One of the diagonal elements

real*8 deriv1(npot,nat,3),deriv1_prev(npot,nat,3)

!   !!! Another for the diagonalization matrix (used in adiabatic)
!   !!! equivalent to <phi_i|\partial\partial t|phi_j>

complex*16 deriv2(npot,npot), deriv2_1(npot,npot)

real*8 prob1(npot),prob2(npot)
real*8 prob(npot),probt,proba(npot)
logical nojump
integer printlevel

!!!!!!!!!!!!!!!!!!!!!!!
!!! Substeps method !!!
!!!!!!!!!!!!!!!!!!!!!!!

integer nsteps,printstep,nstep
!complex*16 V1(npot,npot),dipole1(3,npot,npot), T1(npot,npot)
!complex*16 V2(npot,npot),dipole2(3,npot,npot), T2(npot,npot)
!complex*16 V3(npot,npot),dipole3(3,npot,npot), T3(npot,npot)
real*8 time1,time2,time3,time0
complex*16 U1(npot,npot),U2(npot,npot),U3(npot,npot)
complex*16 Vtmp(npot,npot), Ttmp(npot,npot),diptmp(3,npot,npot)
complex*16 V3(npot,npot)
complex*16 Vtmp0(npot,npot)

!!! Matrices related with the input states, useful when not all gradients are available
complex*16 U0(npot,npot),V0(npot,npot),dipole0(3,npot,npot),T0(npot,npot)

!Splining
!complex*16 Vsp(2*nsteps-1,npot,npot),dipsp(2*nsteps-1,3,npot,npot),Tsp(2*nsteps-1,npot,npot)
complex*16 Vsp(4,npot,npot),dipsp(3,4,npot,npot),Tsp(4,npot,npot)

real*8 rando(1)
!real*8 rando

integer i,j,k,l,n, probjump
complex*16 kk
integer pprev
real*8 Ekt, dEpot, dEpot2
real*8 delta

if (printlevel.ge.3) then
 open(1,file="debug.dat")
 write(1,*) "Entering propag"
 close(1)
endif

!!! Total energy
Etot=Epot
do i=1,nat
 do j=1,3
  Etot=Etot+mass(i)*vel(i,j)**2/2.
 enddo
enddo

!   !!! Saving previous energies and derivatives
do i=1,npot
 do j=1,npot
  do l=npol,2,-1
   V(l,i,j)=V(l-1,i,j)
   T(l,i,j)=T(l-1,i,j)
   do k=1,3
    dipole(l,k,i,j)=dipole(l-1,k,i,j)
   enddo
  enddo
 enddo
enddo

!   !!! Reseting deriv2 and deriv1
deriv1_prev=deriv1
deriv1=0.d0
deriv2=dcmplx(0.d0,0.0d0)
!write (6,*) coef

call velver1(nat,x,vel,deriv1_prev(pot,:,:),mass,dt,printlevel,noat,time-dt)

if (printlevel.ge.3) then
open(1,file="debug.dat")
write(1,*) "Before propag/evaluate"
close(1)
endif

!write(6,*) "Entering propag"
!do i=1,npot
! write(6,*) (coef(i))
!enddo
if (typsh.lt.0 ) call diabat_coef (npot,U,coef)
!write(6,*) "Entering propag"
!do i=1,npot
! write(6,*) (coef(i))
!enddo

write(68,*) "Time  " , time/41.32

U0=dcmplx(0.d0,0.d0)
U=U_prev
call stepwrit(npot,pot,npoint,U)

!!! Single point calculation for analytical gradients
if (typmodel.eq.1) then
!!! Diagonalize with respect to spin-free states to evaluate check the
!!! required gradients
 call singlepoint(printlevel,nat,npot,x,sat,V0,dipole0,T0,vel,U0,las,typdyn)
 V(1,:,:)=V0
 dipole(1,:,:,:)=dipole0
 T(1,:,:)=T0
 egradV=.false.
 gradV=(0.,0.)
 call stepwrit(npot,pot,npoint,U0)
 Vtmp=V0
 call adiabatize(npot,U0,Vtmp)
 do i=1,npot
  W(i)=dreal(Vtmp(i,i))
 enddo
!!! Gradients are going to be calculated at the end checking the last potential
else
 U0=U
 call evaluate(printlevel, sat,nat,npot,x,V0,U0,W,deriv1,deriv2,dipole0, &
 .false.,las,typdyn,typmodel,vel,T0, gradV,egradV, pot )
 V(1,:,:)=V0
 dipole(1,:,:,:)=dipole0
 T(1,:,:)=T0
 call stepwrit(npot,pot,npoint,U0)
endif
!   !!! Using the gradients of Potentials and Dipoles, only for Ehrenfest
if (typsh.eq.0) then
!call ehrenfest (coef, nat,npot, deriv1, V, Epot,deriv1(pot,:,:),npol )
else
 Epot=W(pot)
endif


if (printlevel.ge.3)  then
 open(1,file="debug.dat")
 write(1,*)"Before propag/substeps"
 close(1)
endif
if (printlevel.ge.2) then
 write(106,"(1000000(E20.10e3))") time,(dreal(V(1,i,i)),i=1,npot)
 write(107,"(1000000(E20.10e3))") time,((V(1,i,j),j=1,npot),i=1,npot)
endif

T(1,:,:)=dcmplx(0.d0,0.d0)
if (typmodel.ne.0) then
 inquire(file="overlap.dat",exist=nojump)
 if (nojump) then
  T(1,:,:)=dcmplx(0.d0,0.d0)
!  if (iabs(typsh).eq.4) then
  call createT_diab(printlevel,npot,V(1,:,:),dipole(1,:,:,:),T(1,:,:),Umini,time,dt,Uprop,V(2,:,:))
!   call createT_diab2(printlevel,npot,npol,V,dipole,T,Umini,time,dt,Uprop)
!  else
!!! createT can rotate the V and dipole matrices to take into account diabatic jumps
!   call createT(printlevel,npot,V(1,:,:),dipole(1,:,:,:),T(1,:,:),Umini,time,dt,Uprop)
!  endif
 else
  call minimizeVdipT(npot,V, dipole, T,Umini, time, npol,dt,Uprop)
 endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Beginning of substeps !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time of substeps is delta
delta=dt/float(nsteps)

call splining(time,dt,V,T,dipole,npot,npol,Vsp,Tsp,dipsp,printlevel,nsteps)

potinic=pot
time0=time-dt
do n=1,nsteps

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*) "Before propag/lasers and diag subroutines"
  close(1)
 endif

 time1=time-dt+(n-1)*delta
 time2=time1+delta/2.
 time3=time1+delta

 if (printlevel.ge.3) write(6,*) "Substep ",n,time0,time1,time2,time3


 call evaluatespline(npot,Vsp,time1,time0,Vtmp)
 call evaluatespline(npot,Tsp,time1,time0,Ttmp)
 do k=1,3
  call evaluatespline(npot,dipsp(k,:,:,:),time1,time0,diptmp(k,:,:))
 enddo
 call laser(nlas,tlas,npar,dlas,npart,laspar,las,time1)
 U1=U_prev

 if (printlevel.ge.3) then
  write(201,*) "Time ",time1/41.32
  write(201,*) "Input U1 to diag"
  do i=1,npot
   write(201,900) (U1(i,:))
  enddo
  write(666,*) time1/41.32
  write(667,*) time1/41.32
  write(999,*) time1/41.32
 endif

!!! Preparing the matrices to propagate... inte1,inte2,inte3

 if (printlevel.ge.3) write(69,*) "Time " , time1/41.32
 call diago(npot,Vtmp,diptmp,las,U1,typdyn,printlevel,.false.,2,69,.true.)
 Vtmp0=Vtmp
 call adiabatize(npot,U1,Vtmp0)

 if (printlevel.ge.3) then
  write(201,*) "Output U1 from diag"
  do i=1,npot
   write(201,900) (U1(i,:))
  enddo
 endif

 if (iabs(typsh).eq.4) then
!!! U constant during the substeps
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte1,deriv2, U1,U1, typsh)
 else
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte1,deriv2, U1,U_prev, typsh)
 endif

 if (printlevel.ge.2) then
  call adiabatize(npot,U1,Vtmp)
  write(108,"(1000000(E20.10e3))") time1,(dreal(Vtmp(i,i)),i=1,npot)
  write(109,"(1000000(E20.10e3))") time1,((Vtmp(i,j),j=1,npot),i=1,npot)
 endif

 deriv2_1=deriv2

! write(199,900) time1/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)
! write(197,900) time1/41.32,((dreal(U1(i,j)),j=1,npot),i=1,npot)
! write(196,900) time1/41.32,((dimag(U1(i,j)),j=1,npot),i=1,npot)
! call adiabatize(npot,U1,Vtmp)
! write(198,900) time1/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)

 call evaluatespline(npot,Vsp,time2,time0,Vtmp)
 call evaluatespline(npot,Tsp,time2,time0,Ttmp)
 do k=1,3
  call evaluatespline(npot,dipsp(k,:,:,:),time2,time0,diptmp(k,:,:))
 enddo
 call laser(nlas,tlas,npar,dlas,npart,laspar,las,time2)
 U2=U1

 if (printlevel.ge.3) then
  write(201,*) "Input U2 to diag"
  do i=1,npot
   write(201,900) (U2(i,:))
  enddo
  write(667,*) time2/41.32
  write(666,*) time2/41.32
  write(999,*) time2/41.32
 endif

 if (printlevel.ge.3) write(69,*) "Time " , time2/41.32
 call diago(npot,Vtmp,diptmp,las,U2,typdyn,printlevel,.false.,2,69,.true.)

 if (printlevel.ge.3.and.iabs(typsh).ne.4) then
  write(201,*) "Output U2 from diag"
  do i=1,npot
   write(201,900) (U2(i,:))
  enddo
 endif

 if (iabs(typsh).eq.4) then
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte2,deriv2, U1,U1, typsh)
 else
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte2,deriv2, U2,U1, typsh)
 endif
 if (printlevel.ge.2) then
  call adiabatize(npot,U2,Vtmp)
  write(108,"(1000000(E20.10e3))") time2,(dreal(Vtmp(i,i)),i=1,npot)
  write(109,"(1000000(E20.10e3))") time2,((Vtmp(i,j),j=1,npot),i=1,npot)
 endif
! write(199,900) time2/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)
! write(197,900) time2/41.32,((dreal(U2(i,j)),j=1,npot),i=1,npot)
! write(196,900) time2/41.32,((dimag(U2(i,j)),j=1,npot),i=1,npot)
! call adiabatize(npot,U2,Vtmp)
! write(198,900) time2/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)

 call evaluatespline(npot,Vsp,time3,time0,Vtmp)
 call evaluatespline(npot,Tsp,time3,time0,Ttmp)
 do k=1,3
  call evaluatespline(npot,dipsp(k,:,:,:),time3,time0,diptmp(k,:,:))
 enddo
 call laser(nlas,tlas,npar,dlas,npart,laspar,las,time3)
 U3=U2

 if (printlevel.ge.3) then
  write(201,*) "Input U3 to diag"
  do i=1,npot
   write(201,900) (U3(i,:))
  enddo

!Theese debug files are inside diagonal.f90 without time:
  write(666,*) time3/41.32
  write(667,*) time3/41.32
  write(999,*) time3/41.32
 endif

 if (printlevel.ge.3) write(69,*) "Time " , time3/41.32
 call diago(npot,Vtmp,diptmp,las,U3,typdyn,printlevel,.false.,2,69,.true.)

 if (printlevel.ge.3) then
  write(201,*) "Output U3 from diag"
  do i=1,npot
   write(201,900) (U3(i,:))
  enddo
 endif

 V3=Vtmp
! write(990,*) "TIME ",time3/41.32
 if (iabs(typsh).eq.4) then
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte3,deriv2, U1,U1, typsh)
 else
  call createinter (printlevel,npot, Vtmp, Ttmp,typdyn,delta/2., inte3,deriv2, U3,U2, typsh)
 endif
 if (printlevel.ge.2) then
  call adiabatize(npot,U3,Vtmp)
  write(108,"(1000000(E20.10e3))") time3,(dreal(Vtmp(i,i)),i=1,npot)
  write(109,"(1000000(E20.10e3))") time3,((Vtmp(i,j),j=1,npot),i=1,npot)
 endif

! write(199,900) time3/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)
! write(197,900) time3/41.32,((dreal(U3(i,j)),j=1,npot),i=1,npot)
! write(196,900) time3/41.32,((dimag(U3(i,j)),j=1,npot),i=1,npot)

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*)"Before propag/createinters"
  close(1)
 endif

!do i=1,npot
! write (6,910) (inte1(i,j),j=1,npot)
!enddo
!do i=1,npot
! write (6,910) (inte2(i,j),j=1,npot)
!enddo
!do i=1,npot
! write (6,910) (inte3(i,j),j=1,npot)
!enddo
!910 format (10(E20.8))
!   !!! Saving U2 as U_prev for the next step
 U_prev=U2
 U=U3
 Vtmp=V3
 call adiabatize(npot,U,Vtmp)
! write(198,900) time3/41.32,((dreal(Vtmp(i,j)),j=1,npot),i=1,npot)

!!! Surface Hopping before coefficient propagation
 prob=0.0d0

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*) "Before propag/first_surfacehopp"
  close(1)
 endif

 intea=inte1
 coefa=coef

 if (typsh.lt.0) then
  call adiabat_coef (npot,U1,coefa)
  call adiabatize (npot,U1,intea)
  do i=1,npot
   do j=1,npot
    intea(i,j)=intea(i,j)+deriv2_1(i,j)
   enddo
  enddo
 endif

!! Reset prob
 if (typsh.ne.0) call surfacehopp (npot, pot, coefa, delta, intea, prob, 0, typsh,time1)

!! For the ones that used values before and afer... calculates prob before propagate coefficients
 if (typsh.ne.0) call surfacehopp (npot, pot, coefa, delta, intea, prob, 1, typsh,time1)

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*) "Before propag/rk4"
  close(1)
 endif

 call rk4(npot,coef,inte1,inte2,inte3,delta)
 if (iabs(typsh).eq.4) then
!!! Coefficients transformed to the totally adiabatic representation
!!! As U1 was constant during propagataion, all propagation took place in the U1 representation
!!! it is neccesary to go back to the input representation (diabatize U1) and then adiabatize U3
!  call reorder_energy(npot,V3,U3)
  call diabat_coef(npot,U1,coef)
  call adiabat_coef(npot,U3,coef)
 endif

 if (printlevel.ge.1) call coeffwrit(time3,npot,U3,coef,pot,Uprop,deco)
 if (printlevel.ge.3) then

  popt=0.0d0
  do i=1,npot
   popt=popt+cdabs(coef(i))**2
  enddo
  write(6,*) "End of Subspet, norm: ", popt

  write(98,*) "inte1 ",time1/41.32
  do i=1,npot
   write(98,900) (dreal(inte1(i,j)),dimag(inte1(i,j)),j=1,npot)
  enddo
  write(98,*) "inte2 ",time2/41.32
  do i=1,npot
   write(98,900) (dreal(inte2(i,j)),dimag(inte2(i,j)),j=1,npot)
  enddo
  write(98,*) "inte3 ",time3/41.32
  do i=1,npot
   write(98,900) (dreal(inte3(i,j)),dimag(inte3(i,j)),j=1,npot)
  enddo
  write(98,*) ""

!  write(97,*) time3/41.32
  write(97,900) time3/41.32, ((dreal(Ttmp(i,j)), dimag(Ttmp(i,j)),j=i,npot),i=1,npot)
  write(97,*)""


  write(95,900) time3/41.32, (dreal(Vtmp(j,j)), dimag(Vtmp(j,j)),j=1,npot)
  write(94,900) time3/41.32, (W(j),j=1,npot)

  write(93,*) time3/41.32
  do i=1,npot
   write(93,900) time3/41.32, (dreal(Vtmp(i,j)),dimag(Vtmp(i,j)),j=1,npot)
  enddo
  write(93,*) ""

  write(92,*) time3/41.32
  do i=1,npot
   write(92,900) (dreal(V3(i,j)),dimag(V3(i,j)),j=1,npot)
  enddo
  write(92,*) ""
 endif


909 format (2(F20.4),x,2(F20.4),x,2(F20.4))

!!!!!!!!!!!!!!!!!!!!!!!
!!! Surface Hopping !!!
!!!!!!!!!!!!!!!!!!!!!!!

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*)"Before propag/second_surfacehopp"
  close(1)
 endif

 intea=inte3
 coefa=coef
 Vtmp=V3

 if (typsh.lt.0) then

  call diabatize (npot,U,intea)

  do i=1,npot
   do j=1,npot
    intea(i,j)=intea(i,j)+deriv2(i,j)
   enddo
  enddo

  call diabat_coef (npot,U,coefa)
 endif



!!! Surface hopping after propagation
 if (typsh.ne.0) call surfacehopp (npot, pot, coefa, delta, intea, prob, 2, typsh,time3)

 if (printlevel.ge.3) then
  open(1,file="debug.dat")
  write(1,*)"After propag/second_surfacehopp"
  close(1)
 endif

 i=iabs(typdyn)
 if (i.eq.1.or.typsh.eq.0) prob=0

!   !!! Evaluation of jump
!   !!!Reseting the random number just in case
 rando(1)=1.0
 call random_number(rando)
!write(79,*) rando
!rando(1)=rand(0)
 proba=0.d0
 popt=0.d0
 nojump=.true.

 pprev=pot

 probjump=0
 probt=0.d0
 do i=1,npot
  popt=popt+cdabs(coefa(i))**2
  probt=probt+prob(i)
  proba(i)=probt
 enddo


 do i=1,npot
!   !!! Trying to jump to potential i
  if (nojump) then
   if (rando(1).lt.proba(i).and.pprev.ne.i) then
    probjump=i
    nojump=.false.
   endif
  endif
 enddo

 if (probjump.ne.0) then
!!!! CHECKING with the adiabatic energy after the substeps
! potinic (for dEpot) is saved before the beginning of substeps
! pprev (for dEpot2) is saved just before trying the jumps
  Vtmp=V3
  call adiabatize(npot,U3,Vtmp)
  dEpot=dreal( Vtmp(probjump,probjump)-Vtmp(potinic,potinic) )
  dEpot2=W(probjump)-W(pprev)
 endif

 if (probjump.ne.0) then
!! Is it allowed by energy difference?
  if (dabs(dEpot2).gt.Emin) then
   write(66,901) time3/41.32, pprev,probjump,rando(1),proba(probjump),&
   prob(probjump),probt,Ekt,dEpot,dEpot2,"Forbidden jump","Norma",popt
   probjump=0
  endif
 endif

 if (probjump.ne.0) then
!   !!! Without energy correction
  if (iabs(typdyn).ne.11.and.iabs(typdyn).ne.21) pot=probjump

!   To conserve the energy, one should test if there is available energy
  if (iabs(typdyn).eq.11.or.iabs(typdyn).eq.21) then
!   !!! Energy correction for typdyns= *1
   Ekt=0.d0
   do j=1,3
    do l=1,nat
     Ekt=Ekt+mass(l)*(vel(l,j)**2)/2.
    enddo
   enddo
   if (Ekt.ge.dEpot.and.Ekt.ge.1d-8) then
    pot=probjump
    Epot=dreal(Vtmp(pot,pot))
   else
    write(66,901) time3/41.32, pprev,probjump, rando(1),proba(probjump), &
    prob(probjump), probt, Ekt, dEpot,dEpot2," Frustated","Norma=",popt
    probjump=0
   endif
  endif

 endif

!   !!! We have a jump
 if (probjump.ne.0) write(66,901) time3/41.32, pprev,probjump, rando(1),&
 proba(probjump), prob(probjump), probt, Ekt, dEpot,dEpot2," Success  ","Norma=",popt
!  write(99,*) T
!  write(99,*) deriv2
!  write(99,*)

 if (typsh.gt.0) then
  i=pot
  call reordering(npot,V3,coef,pot,U_prev,U,time3)
  if (printlevel.ge.3) then
   write(201,*) "U3 after energy reordering"
   do j=1,npot
    write(201,900) (U(j,:))
   enddo
  endif
  if (i.ne.pot) then
   write(66,901) time3/41.32,i,pot,0.0d0,1.d0,1.d0,1.d0,0.d0,0.d0,0.d0," Diabatic  ","Norma=",popt
  endif
  Vtmp=V3
  call adiabatize(npot,U,Vtmp)
  Epot=dreal(Vtmp(pot,pot))
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Decoherence factor !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (deco.ge.1e-5) then
  Ekt=0.d0
  do j=1,nat
   do l=1,3
    Ekt=Ekt+mass(j)*vel(j,l)**2/2.
   enddo
  enddo
  if (typsh.lt.0) then
   Vtmp=V3
   call adiabatize (npot,U,Vtmp)
   call adiabat_coef (npot,U,coef)
  endif
  tau=0.d0
  suma=0.d0
  do j=1,npot
   if (j.ne.pot) then
    if (Ekt.ge.1d-8) then
     tau=(1/cdabs(Vtmp(j,j)-Vtmp(pot,pot) )*(1+ abs(deco)/Ekt) )
     coef(j)=coef(j)*dexp( -dt/nsteps/tau )
    endif
    suma=suma+cdabs(coef(j))**2
   endif
  enddo
  coef(pot)=coef(pot)*dsqrt( (popt-suma)/cdabs(coef(pot))**2  )
  if (typsh.lt.0) then
   call diabat_coef(npot,U,coef)
  endif
 endif
 if (deco.le.-1e-5) then
  if (typsh.lt.0) then
   Vtmp=V3
   call adiabatize (npot,U,Vtmp)
   call adiabat_coef (npot,U,coef)
  endif
  suma=0.d0
  do j=1,npot
   if (j.ne.pot) then
!!! Using energy diferences between both paths to calculate decoherence
    tau=cdabs( (Vtmp(j,j)-Vtmp0(j,j))-(Vtmp(pot,pot)-Vtmp0(pot,pot)) )
    coef(j)=coef(j)*dexp( -abs(deco)*tau *dt/nsteps )
    suma=suma+cdabs(coef(j))**2
   endif
  enddo
  coef(pot)=coef(pot)*dsqrt( (popt-suma)/cdabs(coef(pot))**2  )
  if (typsh.lt.0) then
   call diabat_coef(npot,U,coef)
  endif
 endif

 if (printlevel.ge.1) then
  write(10,902) time3/41.32,pot, probt,rando(1),(prob(i),i=1,npot)
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finishing substeps !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo

!!! Born-Oppenheimer
if (iabs(typdyn).eq.1) pot=potinic


!write(6,*) "Exiting propag"
!do i=1,npot
! write(6,*) (coef(i))
!enddo
if (typsh.lt.0) call adiabat_coef (npot,U,coef)
!write(6,*) "Exiting propag"
!do i=1,npot
! write(6,*) (coef(i))
!enddo

!!! Calculating the gradients in analytical model
if (typmodel.eq.1) then
 call stepwrit(npot,pot,npoint,U0)
!! Do not touch U,V,dipole and T
 call evaluate(printlevel, sat,nat,npot,x,V0,U0,W,deriv1,deriv2,dipole0, &
 .false.,las,typdyn,typmodel,vel,T0, gradV,egradV, pot )
endif
call velver2(nat,x,vel,deriv1(pot,:,:),mass,dt,printlevel,noat,time3)
!!! Writting gradients
if (printlevel.ge.1) then
 do i=1,npot
  do j=1,npot
   write(18,*) "POT ",i,j,egradV(i,j),V(1,i,j)
   if (egradV(i,j)) call gradwritc(nat,gradV(i,j,:,:),18)
  enddo
 enddo
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Energy correction out of substeps !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (pot.ne.potinic) then
 if (iabs(typdyn).eq.11.or.iabs(typdyn).eq.21) then
  Ekt=0.d0
  do j=1,3
   do l=1,nat
    Ekt=Ekt+.5*mass(l)*(vel(l,j)**2)
   enddo
  enddo

! dEpot=dreal( Vtmp(pot,pot)-Vtmp(potinic,potinic) )
! dEpot= W(pot)-W(potinic)
  dEpot=W(pot)-(Etot-Ekt)    ! Scaling to maintain the total energy...
  if (printlevel.ge.3) then
  write(99,*) "Correction of energy ",dEpot
  write(99,*) potinic,W(potinic),pot,W(pot)
  endif
  if (dEpot.lt.Ekt) then
   write(68,"(A,F20.2,A,2(A,F12.6))") &
   "JUMP ",time3/41.32," correcting the velocity",&
   " Etot(t-dt)-Etot(t)", dEpot, " dEpot ",W(pot)-W(potinic)
   do j=1,3
    do l=1,nat
     vel(l,j)=dsqrt( 1.-dEpot/Ekt ) * vel(l,j)
    enddo
   enddo
  else
   write(68,"(A,F20.2,A,5(A,F12.6))") &
"JUMP ",time3/41.32," Kinetic energy lower than potential difference, vel to zero, ", &
" Ekt ",Ekt," Etot(t-dt)-Etot(t) ",dEpot," dEpot ",W(pot)-W(potinic)
   do j=1,3
    do l=1,nat
     vel(l,j)=0.d0
    enddo
   enddo
  endif

!!! New kinetic energy
  Ekt=0.d0
  do j=1,3
   do l=1,nat
    Ekt=Ekt+mass(l)*vel(l,j)**2/2.
   enddo
  enddo

 endif
endif

if (typsh.ne.0) Epot=W(pot)

if (printlevel.eq.0) then
 if (printstep.eq.nstep) then
  probt=probt*nsteps*printstep
  do i=1,npot
   prob(i)=prob(i)*nsteps*printstep
  enddo
 endif
 write(10,902) time3/41.32,pot, probt,rando(1),(prob(i),i=1,npot)
endif


!!!! Recalculating NACME with no laser to print them
!!!! Storage them in deriv2 (not used at the beginnig of propag)
las=0.d0
call diago(npot,V(1,:,:),dipole,las,U3,typdyn,printlevel,.true.,2,69,.true.)
deriv2=T(1,:,:)
call adiabatize(npot,U3,deriv2)
call laser(nlas,tlas,npar,dlas,npart,laspar,las,time3)




902 format(F20.6,x,I3.3,x,E20.10e3,1000(x,F10.8))
901 format(F20.10,x,I3.3,x,I3.3,x,7(E20.10e3,x),2(A10,x),10000(X,E20.10E3))
900 format(100000(x,E20.10e3))
903 format(F20.10,x,2(E20.10e3,x),3(I3.3,x))
910 format(F20.10,x,A30,x,I2,x,I2,x,2(F20.10),A5, x,F20.10, A5,x,F20.10)

if (printlevel.ge.3) then
 open(1,file="debug.dat")
 write(1,*) "Finishing propag"
 close(1)
endif

return
end subroutine propag




subroutine createinter (printlevel,npot, V, T,typdyn,delta, inte,deriv2, U,U_prev, typsh)
!!           createinter (debug,npot, Vtmp, Ttmp,typdyn,delta/2., inte3,deriv2, U3,U2, typsh)

implicit none
integer, intent(in) :: npot, typdyn, typsh
real*8, intent(in) :: delta
complex*16, intent(in) :: V(npot,npot),T(npot,npot)
complex*16, intent(inout) :: deriv2(npot, npot)
complex*16, intent(out) :: inte(npot,npot)
complex*16, intent(in) :: U(npot,npot), U_prev(npot, npot)

complex*16 kk(npot,npot)
integer i,j,k


complex*16  W(npot,npot),Tw(npot,npot)
integer printlevel
real*8 kkr,kki

W=V
Tw=T

if (printlevel.ge.3) then
 open(1,file="debug.dat")
 write(1,*)"Before createinter/timederiv"
 close(1)
endif

if (iabs(typdyn).ne.23) then
 call timederiv (npot,deriv2,U,U_prev,delta)
 if (printlevel.ge.3) then
  write(98,*) "Timederiv"
  do i=1,npot
   do j=1,npot
    write(98,*) i,j,U_prev(i,j),U(i,j),deriv2(i,j)
   enddo
  enddo
 endif
endif

!if (debug) then
!open(1,file="debug.dat")
!write(1,*)"Before createinter/adiabatize"
!close(1)
!endif
!   !!! Adiabatizing potential and T matrix with the new sign (if changed in timederiv)


if (typsh.gt.0) then
 call adiabatize(npot,U,W)
 call adiabatize(npot,U,Tw)
 do i=1,npot
  do j=1,npot
   inte(i,j)=deriv2(i,j)+Tw(i,j)
  enddo
 enddo
endif

if (typsh.lt.0) then
 do i=1,npot
  do j=1,npot
   inte(i,j)=Tw(i,j)
  enddo
 enddo
endif

if (printlevel.ge.3) then
 write(98,*) "T part transformed U T U"
 do i=1,npot
  write(98,"(10000(x,F7.3))") (dreal(Tw(i,j)),dimag(Tw(i,j)),j=1,npot)
 enddo
 write(98,*) "V part transformed U V U"
 do i=1,npot
  write(98,"(10000(x,F6.3))") (dreal(W(i,j)),dimag(W(i,j)),j=1,npot)
 enddo
endif

!  Antihermitian deriv2+T
do i=1,npot
 do j=1,npot
  kkr=dreal(inte(i,j))-dreal(inte(j,i))
  kki=dimag(inte(i,j))+dimag(inte(j,i))
  kk(i,j)=dcmplx(kkr/2.,kki/2.)
 enddo
 kk(i,i)=dcmplx(0.d0,0.d0)
enddo

!   !!! Interaction matrix at time t
!$OMP PARALLEL DO COLLAPSE(2)
do i=1,npot
 do j=1,npot
  inte(i,j)=dcmplx(0.d0,0.d0)
  inte(i,j)=inte(i,j)+dcmplx(0.d0,1.d0)*W(i,j)
  inte(i,j)=inte(i,j)+kk(i,j)
  inte(i,j)=-inte(i,j)
 enddo
enddo
!$OMP END PARALLEL DO



901 format(F20.10,x,I3.3,x,I3.3,x,6(E20.10e3,x),2(A10,x),E20.10E3)

return
end subroutine
