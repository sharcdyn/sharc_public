subroutine creategeom(nat,nmod,coord0,modno,disp,coord)
implicit none
integer, intent(in) :: nat,nmod
real*8, intent(in) :: coord0(nat,3),modno(nmod,nat,3),disp(nmod)
real*8, intent(out) :: coord(nat,3)
integer i,j,k

coord=coord0
do i=1,nat
 do j=1,3
  do k=1,nmod
   coord(i,j)=coord(i,j)+modno(k,i,j)*disp(k)
  enddo
 enddo
enddo

end subroutine creategeom

subroutine creategeomvel(nat,nmod,coord0,modno,disp,coord,vel,p,massred)
implicit none
integer, intent(in) :: nat,nmod
real*8, intent(in) :: coord0(nat,3),modno(nmod,nat,3),disp(nmod)
real*8, intent(in) :: p(nmod),massred(nmod)
real*8, intent(out) :: coord(nat,3),vel(nat,3)
integer i,j,k
real*8 vdisp(nmod)

vel=0.
vdisp=0.
do i=1,nmod
 vdisp(i)=p(i)/massred(i)
enddo
 
coord=coord0
do i=1,nat
 do j=1,3
  do k=1,nmod
   vel(i,j)=vel(i,j)+modno(k,i,j)*vdisp(k)
   coord(i,j)=coord(i,j)+modno(k,i,j)*disp(k)
  enddo
 enddo
enddo

end subroutine creategeomvel

subroutine writegeom(nat,sat,coord,un,title)
implicit none
integer, intent(in) :: nat,un
real*8, intent(in) :: coord(nat,3)
character*2, intent(in) :: sat(nat)
character*100000, intent(in) :: title

integer i,j
character*100 fmt

fmt="(A2,3(x,E30.20e3))"

write(un,*) nat
write(un,"(A)") trim(title)
do i=1,nat
 write(un,fmt) sat(i),(coord(i,j)/1.889726,j=1,3)
enddo

end subroutine writegeom

subroutine writegeomvel(nat,sat,coord,un,title,vel)
implicit none
integer, intent(in) :: nat,un
real*8, intent(in) :: coord(nat,3)
character*2, intent(in) :: sat(nat)
character*100000, intent(in) :: title
real*8, intent(in) :: vel(nat,3)

integer i,j
character*100 fmt

fmt="(A2,6(x,E30.20e3))"

write(un,*) nat
write(un,"(A)") trim(title)
do i=1,nat
 write(un,fmt) sat(i),(coord(i,j)/1.889726,j=1,3),(vel(i,j)*1000,j=1,3)
enddo

end subroutine writegeomvel

subroutine writegeomvel2(nat,sat,coord,un,title,vel)
implicit none
integer, intent(in) :: nat,un
real*8, intent(in) :: coord(nat,3)
character*2, intent(in) :: sat(nat)
character*100000, intent(in) :: title
real*8, intent(in) :: vel(nat,3)

integer i,j
character*100 fmt

fmt="(3(x,E30.20e3))"

write(un,"(A)") trim(title)
do i=1,nat
 write(un,fmt) (coord(i,j),j=1,3)
enddo
do i=1,nat
 write(un,fmt) (vel(i,j),j=1,3)
enddo

end subroutine writegeomvel2
