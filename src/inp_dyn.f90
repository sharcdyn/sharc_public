subroutine input(typdyn,typmodel,typsh,npol,deco,printlevel,SHminE,nsteps)

use molecule_

implicit none
integer, intent(out) :: typdyn,typmodel,typsh,npol
real*8, intent(out) :: deco,SHminE
integer, intent(out) :: printlevel,nsteps

character*100 model,dyntype,SHcorr,TSH,damp,coef_adiab,nacme

namelist /parameters/ model,dyntype,SHcorr,TSH,damp,coef_adiab,&
                      nacme,deco,npol,SHminE,printlevel,&
                      nsteps

logical ex

!!! DEFAULT VALUES
nsteps=51
npol=5
deco=0.1d0
model="analytical"
dyntype="SHARC"
SHcorr="yes"
TSH="dens_bef-after"
damp="no"
coef_adiab="yes"
nacme="overlap"
printlevel=0
!!! No limitation in the jump energy
SHminE=100000

write(6,*) "SHARC: Richter et al. J. Chem. Theory Comput.7, 1253 (2011)"

inquire (file="parameters.dat",exist=ex)
if (ex) then
 open(1,file="parameters.dat")
 read(1,NML=parameters)
 close(1)
endif

write(6,*) "Decoherence parameter ",abs(deco)
if (deco.lt.0.) write(6,*) "The wavepacket will not be renormalized, loose of norm is expected"

select case (model)
 case("model")
  typmodel=0
 case("analytical")
  typmodel=1
 case("numerical")
  typmodel=2
 case default
  write(6,*) "Possible values of model:"
  write(6,"(A)") "model: Electronic properties are evaluated using a FORTRAN subroutine, see nomodel.f90 as template."
  write(6,*) "analytical: Electronic properties are get from an ab initio code using analytical derivates"
  write(6,*) "numerical: Electronic properties are evaluated numerical doing several single point calculations. VERY EXPENSIVE"
  stop "model is not valid"
end select

select case (dyntype)
 case ("SHARC")
  if (SHcorr=="yes") then
   if (nacme=="overlap") then
    typdyn=21
   else if (nacme=="NACME") then
    typdyn=23
   else
    stop "nacme can only be overlap or NACME"
   endif
  else if (SHcorr=="no") then
   if (nacme=="overlap") then
    typdyn=22
   else if (nacme=="NACME") then
    stop "Not implemented"
   else
    stop "nacme can only be overlap or NACME"
   endif
  else
   write(6,*) "SHcorr only yes or no"
   stop
  endif
 case("FISH")
  if (SHcorr=="yes") then
   if (nacme=="overlap") then
    typdyn=11
   else if (nacme=="NACME") then
    stop "Not implemented"
   else
    stop "nacme can only be overlap or NACME"
   endif
  else if (SHcorr=="no") then
   if (nacme=="overlap") then
    typdyn=12
   else if (nacme=="NACME") then
    stop "Not implemented"
   else
    stop "nacme can only be overlap or NACME"
   endif
  else
   write(6,*) "SHcorr only yes or no"
   stop
  endif
 case("BO")
   typdyn=1
 case default
  write(6,*) "Possible values for the dyntype variable:"
  write(6,*) "BO: Born-Oppenheimer"
  write(6,*) "SHARC"
  write(6,*) "FISH"
  stop "dyntype no valid"
end select

select case (TSH)
! case ("Erhenfest")
!  typsh=0
 case ("Tully")
  typsh=2
 case ("advTully")
  typsh=1
 case ("dens_after")
  typsh=3
 case ("dens_bef-after")
  typsh=4
 case default
  write(6,*) "Possible values of TSH:"
  write(6,*) "Erhenfest, not ready yet"
  write(6,*) "Tully:"
  write(6,*) "afvTully:"
  write(6,*) "dens_after"
  write(6,*) "dens_bef-after"
  stop "TSH value not valid"
end select

select case (damp)
 case("no")
  typdyn=typdyn
 case("yes")
  typdyn=-typdyn
 case default
  stop "damp value not valid, only yes or no"
end select
select case (coef_adiab)
 case("yes")
  typsh=typsh
 case("no")
  typsh=-typsh
 case default
  stop "Propagate in adiabatic, yes or no"
end select

select case(printlevel)
 case (0)
  write(6,*) "Regular printing"
 case (1)
  write(6,*) "Printing data every substep"
 case(2)
  write(6,*) "Semi debug mode"
 case(3)
  write(6,*) "Denug mode"
 case default
  stop "printlevel can be 0, 1, 2, 3"
end select


if (typdyn.eq.0) write(6,*) "Using model potentials"
if (typdyn.eq.1) write(6,*) "Using direct dynamics with analytical gradients"
if (typdyn.eq.2) write(6,*) "Using direct dynamics with numerical gradients"


!   !!! Non Surface-Hopping methods:
if (typdyn.eq.1.or.typdyn.eq.-1) write(6,*) "Using Born-Oppenheimer Dynamics (no SH)"

!   !!! Ehrenfest
if (iabs(typdyn).ge.10.and.iabs(typdyn).lt.20.and.typsh.eq.0) write(6,*) "Using FISH with Ehrenfest Dynamics (Mean-field)"
if (iabs(typdyn).ge.20.and.typsh.eq.0) write(6,*) "Using SHARC with Ehrenfest Dynamics (Mean-field)"

!   !!! Surface Hopping type Dynamics:
if (typdyn.eq.11.or.typdyn.eq.-11) write(6,*) "Using FISH (with energy correction)"
if (typdyn.eq.12.or.typdyn.eq.-12) write(6,*) "Using FISH (without correction)"

if (typdyn.eq.21.or.typdyn.eq.-21) write(6,*) "Using SHARC (with energy correction)"
if (typdyn.eq.22.or.typdyn.eq.-22) write(6,*) "Using SHARC (without correction)"
if (typdyn.eq.23.or.typdyn.eq.-23) write(6,*) "Using SHARC (with coordinate dependency)"

if (iabs(typsh).eq.1) write(6,*) "Using Tully SH before and after coefficients propagation"
if (iabs(typsh).eq.2) write(6,*) "Using Tully SH only after coefficients propagation"
if (iabs(typsh).eq.3) then
 write(6,*) "Using density SH after coefficients propagation"
 write(6,*) "Mitric PCCP 14, 8299 (2012)"
endif
if (iabs(typsh).eq.4) then
 write(6,*) "Using density SH before and after coefficients propagation"
 write(6,*) "Mitric PCCP 14, 8299 (2012)"
 write(6,*) "In this case, the propagation is done 'diabatically' in each electronic substep"
endif
if (typsh.lt.0) write(6,*) "PROPAGATING THE COEFICIENTS WITHOUT ADIABATIZATION"

if (SHminE.ge.1000) then
 write(6,*) "No limitation in the energy difference on the Surface Hopping"
else if (SHminE.ge.0) then
 write(6,*) "Minimum energy to use SH algorithm SHminE (eV)",SHminE
else
 stop 'SHminE cannot be negative'
endif

SHminE=SHminE/27.211

return

end

subroutine readgeom(nat,npot,sat,noat,x,mass,vel,pot,coef)
implicit none
integer, intent(in) :: nat,npot
integer, intent(inout) :: pot
character*2, intent(out) :: sat(nat)
real*8, intent(out) :: x(nat,3),vel(nat,3),mass(nat),noat(nat)
complex*16, intent(out) :: coef(npot)

logical ex1,ex2
real*8 kk(npot),norm,rnd
integer i,j
integer pot_trial
real*8 prob

coef=dcmplx(0.d0,0.d0)

inquire(file="geom.sharc",exist=ex1)
inquire(file="geom",exist=ex2)
if (ex1.or.ex2) then
 if (ex1) then
  open(51,file="geom.sharc")
 else
  open(51,file="geom")
 endif
else
 stop "A geometry file is compulsary (geom or geom.sharc)"
endif

pot_trial=-1

!   !!! Read intial conditions
do i=1, nat
 read(51,*) sat(i), noat(i), x(i,1), x(i,2), x(i,3), mass(i)
 mass(i)=mass(i)*1822.888
 do j=1,3
  vel(i,j)=0.d0
 enddo
enddo
ex1=.true.
do i=1,nat
 read(51,*,END=10)  vel(i,1), vel(i,2), vel(i,3)
enddo
ex1=.false.


read(51,*,END=10) pot_trial

read(51,*,END=10) (kk(i),i=1,npot)
do i=1,npot
 coef(i)=dcmplx(kk(i),0.d0)
enddo
read(51,*,END=10) (kk(i),i=1,npot)
do i=1,npot
 coef(i)=coef(i)+dcmplx(0.d0,kk(i))
enddo

10 continue

norm=0.
do i=1,npot
 norm=norm+abs(coef(i))**2
enddo

if (pot.le.0) then
 if (pot_trial.gt.0) then
  write(6,*) "Trajectory starting at potential ",pot_trial
  pot=pot_trial
  if (norm.le.1e-3) then
   write(6,*) "No coefficients, so assuming one in the initial potential"
   coef=dcmplx(0.d0,0.d0)
   coef(pot)=dcmplx(1.,0.)
  endif
 else if (pot_trial.eq.0) then
  write(6,*) "Starting potential is 0, collapsing"
  if (norm.le.1e-3) then
   write(6,*) "No coefficients, so nothing to collapse, aborting"
   stop 'Collapsing is not possible without coefficients'
  endif
  do i=1,npot
   coef(i)=coef(i)/sqrt(norm)
  enddo
  call random_number(rnd)
  prob=0.
  do i=1,npot
   prob=prob+abs(coef(i))**2
   if (rnd.le.prob.and.rnd.gt.prob-abs(coef(i))**2) then
    pot=i
    write(6,*) "Collapsing wavefunction to potential ",pot
    write(6,*) "Prob ",abs(coef(i))**2
    write(6,*) "Accumulated prob ",prob
    write(6,*) "Random number ",rnd
   endif
  enddo
  do i=1,npot
   coef(i)=coef(i)*sqrt(norm)
  enddo
 else
  write(6,*) "Geom file has not potential... aborting"
  stop 'No potential to start'
 endif
else
 if (norm.le.1e-3) then
  write(6,*) "No coefficients, so assuming on in the initial potential ",pot
  coef=dcmplx(0.d0,0.d0)
  coef(pot)=dcmplx(1.d0,0.d0)
 endif
endif
 

if (pot_trial.eq.0) then
endif

if (ex1) write(6,*) "Warning vel are not present, considering 0"

end


subroutine inputhelp()
implicit none

write(6,*) "npol: Variable to select the number of steps used in the splinning procedure"
write(6,*) "model: Defines the type of dynamics employed (default: analytical; other possibilities: numerical model)"

end
