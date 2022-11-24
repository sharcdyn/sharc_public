subroutine gauss(npot,pot,dipole,Epot,nEn,E0,Et,spec,gwidth,prob,un)
implicit none
integer npot,pot,nEn,un
complex*16 dipole(3,npot,npot)
real*8 spec(3,npot,nEn),Epot(npot)
real*8 E0,Et,gwidth,prob(npot)


integer i,j,k
real*8 dE,osc(3),Eactual,fosc


do i=1,npot
 dE=Epot(i)-Epot(pot)
 prob(i)=0.d0
 do j=1,3
  osc(j)=cdabs(dipole(j,pot,i))**2
  osc(j)=2./3.*dE*osc(j)
  prob(i)=prob(i)+osc(j)
 enddo
 do j=1,nEn
  Eactual=E0+(j-1)*(Et-E0)/float(nEn-1)
  fosc=gwidth/2./dsqrt(2.*dlog(2.d0))
  fosc=(dE-Eactual)/fosc
  fosc=dexp(-fosc*fosc/2.)
  do k=1,3
   spec(k,i,j)=spec(k,i,j)+osc(k)*fosc
  enddo
 enddo
!!! Writting probabilities by traj
 if (un.ne.0) write(un,"(2(x,E30.20e3))") dE,prob(i)
enddo
if (un.ne.0) write(un,*) ""

end
