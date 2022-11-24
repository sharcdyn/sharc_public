      subroutine diag(VWt,npot,VWdt)      
c     npot--->Dimensiones (cuadrada)
c     VWt --->Matriz a diagonalizar, sale matriz de diagonalizacion
c     VWdt--->Vector con la diagonal
      implicit none
      integer npot
      real*8 VWt(npot,npot),VWdt(npot)

      real*8 E(npot),kk(npot,npot)

      call TRED2(VWt,npot,npot,VWdt,E)

      call TQLI(VWdt,E,npot,npot,VWt)

      return
      end
c----------------------------------------------------------------
      subroutine diag_order(VWt,npot,VWdt)
      implicit none
      integer npot
      real*8 VWt(npot,npot),VWdt(npot)
      
      real*8 kk
      integer i,j,k

      do i=1,npot
       do j=i+1,npot
        if (VWdt(i).ge.VWdt(j)) then
         kk=VWdt(i)
         VWdt(i)=VWdt(j)
         VWdt(j)=kk
         do k=1,npot
          kk=VWt(k,i)
          VWt(k,i)=VWt(k,j)
          VWt(k,j)=kk
         enddo
        endif
       enddo
      enddo

      end subroutine diag_order

c-----------------------------------------------------------------
      SUBROUTINE TRED2(A,N,NP,D,E)
      implicit real*8 (A-H,O-Z)
      DIMENSION A(NP,NP),D(NP),E(NP)
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1  
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+DABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-DSIGN(DSQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
c-------------------------------------------------------------------
      SUBROUTINE TQLI(D,E,N,NP,Z)
      implicit real*8 (A-H,O-Z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=DABS(D(M))+DABS(D(M+1))
            IF (DABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G))THEN
                C=G/F
                R=DSQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=DSQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
