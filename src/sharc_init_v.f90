! Program to generate the vibrational quantum number on each vibrational mode

program sharc_init_v
implicit none
integer nvib,ntraj,nseed,un
integer, allocatable :: vib(:),seed(:),unvib(:)
real*8, allocatable :: frec(:),massred(:)
character*100 fich
integer i,j

write(6,*) "Number of vibrational modes to be sampled and number of trajectories"
read(5,*) nvib,ntraj
allocate(vib(nvib))
allocate(frec(nvib),massred(nvib),unvib(nvib))
write(6,*) "Quantum normal modes to be sampled"
write(6,*) "If you use a negative number:"
write(6,*) "If it is a imaginary mode, the sign of momentum will be reversed"
write(6,*) "If it is a bound mode, only momentum will be sampled"
read(5,*) vib(:)

write(6,*) "Reading information for nm files"
write(6,*) "Creating directly the nm -trajs files for a classical sampling"
do i=1,nvib
 write(fich,"(A,I0,A)") "nm_",iabs(vib(i)),"-info.xyz"
 open(newunit=nseed,file=fich,status="old",iostat=j)
 if (j.ne.0) then
  write(6,*) "File with normal mode does not exist ",trim(fich)
  stop
 endif
 read(nseed,*) j
 read(nseed,*) j,frec(i),massred(i)
 close(nseed)
 write(fich,"(A,I0,A)") "nm_",iabs(vib(i)),"-trajs.out"
 open(newunit=unvib(i),file=fich,status="replace")
 write(unvib(i),"(A,1(I0,x),E30.20e3)") "# ",ntraj,frec(i)
enddo

open(newunit=un,file="v-seeds.dat",status="old",iostat=i)
if (i.ne.0) then
 open(newunit=un,file="v-seeds.dat",status="new")
 write(6,*) "Creating seeds in file v-seeds.dat"
!!! To avoid interface, the seed generation is the main program
 call random_seed(size=nseed)
 write(un,*) nseed
 allocate (seed(nseed))
 open(newunit=j,file="/dev/urandom",status="old",iostat=i,access='stream', form='UNFORMATTED')
 if (i.eq.0) then
  read(j) seed
  close(j)
 else
  call system_clock(count=j)
  seed = j + 37 * (/ (i - 1, i = 1, nseed) /)
 endif
 write(un,*) (seed(i),i=1,nseed)
 close(un)
 deallocate(seed)
endif
close(un)

 

write(6,*) "Output file"
read(5,*) fich
open(newunit=un,file=fich,status="replace")
write(un,*) nvib,ntraj

write(6,*) "Type of ensemble : microcanical"
read(5,*) fich
select case(fich)
 case("microcanonical")
  call microcanonical(nvib,frec,massred,ntraj,un,unvib,vib)
 case default
  stop "Ensemble not prepared"
end select

end program

subroutine microcanonical(nvib,frec,massred,ntraj,un,unvib,vib)
implicit none
integer, intent(in) :: nvib,ntraj,un
real*8, intent(in) :: frec(nvib),massred(nvib)
integer, intent(in) :: unvib(nvib),vib(nvib)

integer nseed
integer, allocatable :: seed(:)
integer qvib(nvib)
real*8 E(nvib),Exc,Ev(nvib),Ebound

integer itraj,ivib
real*8 rando(nvib),ptot,ptot_cont
character*100 fmat
real*8 q,p,amp,ZPE

write(6,*) "Energy beyond the ZPE to be sampled (in Hartree)"
read(5,*) Exc

open(newunit=ivib,file="v-seeds.dat",status="old",iostat=itraj)
if (itraj.ne.0) stop "Bug in the code, I cannot open v-seeds.dat"
read(ivib,*) nseed
allocate(seed(nseed))
read(ivib,*) seed(:)
close(ivib)
write(fmat,"(A,I0,A)") "(",nvib,"(x,I0,2(x,E30.20e3)))"
call random_seed(put=seed)

itraj=0
do while(itraj.le.ntraj)
 call random_number(rando) 
 ptot=0.
 do ivib=1,nvib
  ptot=ptot+rando(ivib)
 enddo
 do ivib=1,nvib
  E(ivib)=Exc*rando(ivib)/ptot
 enddo
!!! First let's put vibrational energy in bound states
 Ebound=0.
 ZPE=0.
 do ivib=1,nvib
  if (frec(ivib).ge.0) then
   qvib(ivib)=int(E(ivib)/frec(ivib)+0.5)
   Ev(ivib)=(qvib(ivib)+.5)*frec(ivib)
!!! Adding the ZPE to the Energy
   E(ivib)=E(ivib)+.5*frec(ivib)
   Ebound=Ebound+E(ivib)
   ZPE=ZPE+.5*frec(ivib)
  endif
 enddo
!! Only continue if there is enough energy
 if (Ebound.le.Exc+ZPE) then
  itraj=itraj+1
!!! Energy to share between the non bound states Exc-Ebound
  ptot_cont=0.
  do ivib=1,nvib
   if (frec(ivib).le.0) then
    ptot_cont=ptot_cont+rando(ivib)
   endif
  enddo
  do ivib=1,nvib
   if (frec(ivib).le.0) then
    Ev(ivib)=(Exc+ZPE-Ebound)*rando(ivib)/ptot_cont
    qvib(ivib)=-1
   endif
  enddo
  write(un,fmat) (qvib(ivib),Ev(ivib),E(ivib),ivib=1,nvib)
!!! Classical sampling !!!
  call random_number(rando) 
  do ivib=1,nvib
!!! Sampling only momentum
   q=0.
   p=sqrt(2.*E(ivib)*massred(ivib))
   if (frec(ivib).le.0) then
    if (vib(ivib).le.0) p=-p
   else if (vib(ivib).le.0) then
    if (rando(ivib).le.0.5) p=-p
   else
!!! Classical sampling 
    amp=sqrt(2*E(ivib)/massred(ivib))/frec(ivib)
    write(6,*) "Energy amp ",E(ivib),amp
    q=amp*cos(2*acos(-1.)*rando(ivib))
    p=-frec(ivib)*amp*sin(2*acos(-1.)*rando(ivib))
    p=p*massred(ivib)
   endif
!!!!!!!!!!!!!!!!!!!!!!!
   write(unvib(ivib),"(4(x,E30.20e3))") q,p,frec(ivib)**2*massred(ivib)/2.*q*q,p*p/2./massred(ivib)
  enddo
 endif
enddo

end subroutine microcanonical
