      subroutine cspline(nor,np,ro,dat_or,a,b,c,d)
      implicit none
      integer nor,np
      real*8 ro(nor)
      real*8 dat_or(nor)

      real*8 h(np-1)
      real*8 matriz(np-2,np-2),y(np-2)
      real*8 S(np),Sb(np-2)
      real*8 a(nor),b(nor),c(nor),d(nor)

      real*8 ref

      integer i,j,k

c     Para facilitar las cosas, vamos a poner como 0 el penultimo valor y anadirlo a d al final
      ref=dat_or(np-1)
      do i=1,np
       dat_or(i)=dat_or(i)-ref
      enddo
      
c     Defino los incrementos
      do i=1,np-1
       h(i)=ro(i+1)-ro(i)
       if (h(i).le.1e-30) then
        write(6,*) "Ojo con el dato ",i, " delta ",h(i)
        stop
       endif
      enddo
c     Defino la matriz de soluciones y
      do i=1,np-2
c      if (i.eq.1) then
c       y(1)=(dat_or(2)-dat_or(1))/h(1)
c     else if (i.eq.np) then
c       y(np)=(dat_or(np)-dat_or(np-1))/h(np-1)
c      else
      y(i)=(dat_or(i+2)-dat_or(i+1))/h(i+1)-(dat_or(i+1)-dat_or(i))/h(i)
c      endif
      y(i)=6*y(i)
      enddo
c     write(6,*) "Matriz Y"
c     write(6,*) (y(i),i=1,np-2)

c     Defino la matriz
      do i=1,np-2
       do j=1,np-2
        matriz(i,j)=0.0d0
       enddo
      enddo
      do i=1,np-2
c      if (i.eq.1) then
c       matriz(1,1)=
c    =  ((h(1)+h(2))*(h(1)+2*h(2))/h(2) * (h(2)*h(2)-h(1)*h(1))/h(2)
c       matriz(1,1)=2*(0+h(1))
c       matriz(i,i+1)=h(i)
c       matriz(i+1,i)=h(i)
c      else if (i.eq.np) then
c       matriz(np,np)=
c    =  ((h(np-2)+h(np-1)*(h(np-1)+2*h(np-2))/h(np-2)
c    *  (h(np-2)*h(np-2)-h(np-1)*h(np-1))/h(np-2)
c       matriz(np,np)=2*(h(np-1)+0)
c      else  
       matriz(i,i)=2*(h(i)+h(i+1))
       if (i.ne.1) then
        matriz(i,i-1)=h(i)
        matriz(i-1,i)=h(i)
       endif
      enddo
c     write(6,*) "Matriz de intervalos"
c     do i=1,np-2
c      write(6,*) (matriz(i,j),j=1,np-2)
c     enddo
c      call matipd (matriz,np-2,np-2,i)
      call matinver(matriz,np-2,i)
     
      if (i.ne.0) then
       write(6,*) "Inverse calculation fails ",i
       stop 'Inverse calculation fails'
      endif
c     write(6,*) k
c     do i=1,np-2
c      write(6,*) (matriz(i,j),j=1,np-2)
c     enddo
      do i=1,np-2
       Sb(i)=0.d0
       do j=1,np-2
        Sb(i)=Sb(i)+matriz(i,j)*y(j)
       enddo
      enddo
      S(1)=0.d0
      S(np)=0.d0
      do i=2,np-1
       S(i)=Sb(i-1)
      enddo
c     write(6,*) "Matriz S"
c     write(6,*) (S(i),i=1,np)

c    Poner otra vez la referencia
      do i=1,np
       dat_or(i)=dat_or(i)+ref
      enddo
c     Datos del spline y las derivadas
      do k=1,np-1
       a(k)=(S(k+1)-S(k))/(6*h(k))
       b(k)=S(k)/2.
       c(k)=(dat_or(k+1)-dat_or(k))/h(k)-(2*h(k)*S(k)+h(k)*S(k+1))/6.
       d(k)=dat_or(k)
      enddo

      return
      end

c------------------------------------------------
      subroutine matinver(matriz,n,istat)
      integer n,istat
      real*8 matriz(n,n)
      real*8 eigm(n),mat2(n,n)
      integer i,j

      real*8 mat0(n,n)
c     mat0=matriz
      
      call diag(matriz,n,eigm)
       istat=0
      do i=1,n
       do j=1,n
        mat2(i,j)=0
       enddo
       if (abs(eigm(i)).le.1e-12) then
        istat=1
       endif
       mat2(i,i)=1./eigm(i)
      enddo
      call dmatmat(n,n,n,matriz,mat2,mat0)
      matriz=transpose(matriz)
      call dmatmat(n,n,n,mat0,matriz,mat2)
ccccc
c     matriz=matmul(mat0,mat2)
c     do i=1,n
c      do j=1,n
c       write(6,*) i,j,matriz(i,j)
c      enddo
c     enddo
ccccc
      matriz=mat2
      return
      end subroutine matinver
      

