C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SCALARD(SCALEJ,A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     + A3X,A3Y,A3Z,U,V,W,C,UROT,VROT,WROT,IMAX,JMAX,KMAX,IN,JN,KN,
     + ARTSSP,IPRESC)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,KSTRID,IL,I,J,K,KK,KA,JJ,
     +     IN,JN,KN,IPRESC

      REAL :: A1K,A2K,A3K,A1XK,A1YK,A1ZK,A2XK,A2YK,A2ZK,A3XK,A3YK,A3ZK,
     +     CCII,UCII,VCII,WCII,UCROT,VCROT,WCROT,CHAT,ARTSSP

      REAL :: A1(*),A1X(*),A1Y(*),A2(*),A2X(*),A2Y(*),U(*),V(*),W(*),
     +        C(*),A1Z(*),A2Z(*),A3(*),A3X(*),A3Y(*),A3Z(*),
     +        UROT(*),VROT(*),WROT(*),SCALEJ(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID

      DO  K = 0,KMAX+1
      KA      = (KN+K-1)*IL
      DO  J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO  I = 0,IMAX+1
      KK      = JJ + I


C ... CONTRAVARIANT VELOCITY BASED SCALING

      IF(IPRESC == 1) THEN
          CHAT = 340. !ARTSSP
      ELSE IF(IPRESC == 2) THEN
          CHAT =  C(KK)
      ELSE IF(IPRESC == 0) THEN
          CHAT =  C(KK)
      ELSE
          WRITE(*,*) 'Wrong IPRESC option. Exiting..'
          STOP
      ENDIF

      A1K     = .5*(A1(KK) + A1(KK+1))
      A2K     = .5*(A2(KK) + A2(KK+ISTRID))
      A3K     = .5*(A3(KK) + A3(KK+IL))
      A1XK    = .5*(A1X(KK)*A1(KK) + A1X(KK+1)     *A1(KK+1))
      A1YK    = .5*(A1Y(KK)*A1(KK) + A1Y(KK+1)     *A1(KK+1))
      A1ZK    = .5*(A1Z(KK)*A1(KK) + A1Z(KK+1)     *A1(KK+1))
      A2XK    = .5*(A2X(KK)*A2(KK) + A2X(KK+ISTRID)*A2(KK+ISTRID))
      A2YK    = .5*(A2Y(KK)*A2(KK) + A2Y(KK+ISTRID)*A2(KK+ISTRID))
      A2ZK    = .5*(A2Z(KK)*A2(KK) + A2Z(KK+ISTRID)*A2(KK+ISTRID))
      A3XK    = .5*(A3X(KK)*A3(KK) + A3X(KK+IL)    *A3(KK+IL))
      A3YK    = .5*(A3Y(KK)*A3(KK) + A3Y(KK+IL)    *A3(KK+IL))
      A3ZK    = .5*(A3Z(KK)*A3(KK) + A3Z(KK+IL)    *A3(KK+IL))
      UCROT   = .5*(A1(KK)*UROT(KK) + A1(KK+1)     *UROT(KK+1))
      VCROT   = .5*(A2(KK)*VROT(KK) + A2(KK+ISTRID)*VROT(KK+ISTRID))
      WCROT   = .5*(A3(KK)*WROT(KK) + A3(KK+IL)    *WROT(KK+IL))
      UCII    = A1XK*U(KK) + A1YK*V(KK) + A1ZK*W(KK) - UCROT
      VCII    = A2XK*U(KK) + A2YK*V(KK) + A2ZK*W(KK) - VCROT
      WCII    = A3XK*U(KK) + A3YK*V(KK) + A3ZK*W(KK) - WCROT
      CCII    = (A1K + A2K + A3K)*CHAT
c      SCALEJ(KK)  = SQRT(UCII**2+VCII**2+WCII**2+CCII**2)
      SCALEJ(KK)  = SQRT(UCII**2+VCII**2+WCII**2)+CCII

      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE SCALARD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TIME3(A1,A2,A3,D1,D2,D3,VOL,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     + A3X,A3Y,A3Z,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDI1,IDI2,IDI3,M,TIMEL,
     + PVISC,ZZZ,IVISC,FRSPRE,CDIFF)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A1(*),A1X(*),A1Y(*),A2(*),A2X(*),A2Y(*),DTL(*),P(*),
     + VOL(*),RO(*),U(*),V(*),W(*),C(*),VIS(*),VIST(*),
     + A1Z(*),A2Z(*),A3(*),A3X(*),A3Y(*),A3Z(*),D1(*),D2(*),D3(*),
     + UROT(*),VROT(*),WROT(*),ZZZ(*)
      LOGICAL :: TIMEL

      EPS     = 1.E-10
C
C ... TIME-STEP OF THE IMPLICIT BIDIAGONAL METHOD (3-D)
C
C     CDIFF   = GAMMA/(PR+EPS)*1.
C     CDIFF   = 2.
      PCFL    = 1./CFL
C     PVISC   = 1.8            !*8.**(M-1)
      T       = 0.
      PRS     = PR/PRT

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      NTOT    = ISTRID*JSTRID*KSTRID

      CALL ZEROZZ(ZZZ,NTOT)

C ... FIND MAXIMUN SIZE OF THE CELLS
      IF(IDI1 /= 0) THEN
      DO 1000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 1000 I = 1,MIN(IDI1,IMAX)
         KK      = JJ + I
         RNUM    = .5*(A1(KK) + A1(KK+1))
         ZZZ(KK) = RNUM
 1000 CONTINUE
      ENDIF

      IF(IDI2 /= 0) THEN
      DO 2000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 2000 J = 1,MIN(IDI2,JMAX)
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 2000 I = 1,IMAX
         KK      = JJ + I
         RNUM    = .5*(A2(KK) + A2(KK+ISTRID))
         ZZZ(KK) = MAX(RNUM,ZZZ(KK))
 2000 CONTINUE
      ENDIF

      IF(IDI3 /= 0) THEN
      DO 3000 K = 1,MIN(IDI3,KMAX)
      KA      = (KN+K-1)*IL
      DO 3000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 3000 I = 1,IMAX
         KK      = JJ + I
         RNUM    = .5*(A3(KK) + A3(KK+IL))
         ZZZ(KK) = MAX(RNUM,ZZZ(KK))
 3000 CONTINUE
      ENDIF

      DO 4000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 4000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 4000 I = 1,IMAX
      KK      = JJ + I

      PVOL    = 1./VOL(KK)
      YPRO    = 1./(RO(KK)+ EPS)

C ... CONVECTIVE TIME STEP
      A1K     = .5*(A1(KK) + A1(KK+1))
      A2K     = .5*(A2(KK) + A2(KK+ISTRID))
      A3K     = .5*(A3(KK) + A3(KK+IL))
      A1XK    = .5*(A1X(KK)*A1(KK) + A1X(KK+1)     *A1(KK+1))
      A1YK    = .5*(A1Y(KK)*A1(KK) + A1Y(KK+1)     *A1(KK+1))
      A1ZK    = .5*(A1Z(KK)*A1(KK) + A1Z(KK+1)     *A1(KK+1))
      A2XK    = .5*(A2X(KK)*A2(KK) + A2X(KK+ISTRID)*A2(KK+ISTRID))
      A2YK    = .5*(A2Y(KK)*A2(KK) + A2Y(KK+ISTRID)*A2(KK+ISTRID))
      A2ZK    = .5*(A2Z(KK)*A2(KK) + A2Z(KK+ISTRID)*A2(KK+ISTRID))
      A3XK    = .5*(A3X(KK)*A3(KK) + A3X(KK+IL)    *A3(KK+IL))
      A3YK    = .5*(A3Y(KK)*A3(KK) + A3Y(KK+IL)    *A3(KK+IL))
      A3ZK    = .5*(A3Z(KK)*A3(KK) + A3Z(KK+IL)    *A3(KK+IL))
      UCROT   = .5*(A1(KK)*UROT(KK) + A1(KK+1)     *UROT(KK+1))
      VCROT   = .5*(A2(KK)*VROT(KK) + A2(KK+ISTRID)*VROT(KK+ISTRID))
      WCROT   = .5*(A3(KK)*WROT(KK) + A3(KK+IL)    *WROT(KK+IL))
      UCII    = A1XK*U(KK) + A1YK*V(KK) + A1ZK*W(KK) - UCROT
      VCII    = A2XK*U(KK) + A2YK*V(KK) + A2ZK*W(KK) - VCROT
      WCII    = A3XK*U(KK) + A3YK*V(KK) + A3ZK*W(KK) - WCROT
      IF(IVISC <= 2) THEN
      B1      = PVOL*(ABS(UCII) + ABS(VCII) + ABS(WCII)
c     +        + 25.*SQRT(A1K**2 + A2K**2 + A3K**2))*PCFL
     +        + C(KK)*SQRT(A1K**2 + A2K**2 + A3K**2))*PCFL

      ELSE ! only convection
      B1      = PVOL*(ABS(UCII) + ABS(VCII) + ABS(WCII))*PCFL
      ENDIF

C ... VISCOUS TIME STEP
      RNUM    = ZZZ(KK)**2*PVOL**2*YPRO   ! 1/d**2
      T       = RNUM*(VIS(KK)+ PRS*VIST(KK))/CDIFF ! Changed
      B2      = 2.*T*PVISC
      IF(IVISC == 2 .OR. IVISC == 4) B2 = 0. ! No viscous time-step

C ... COMBINATIONS
CC    B11     = B1 + B2
      B11     = MAX(B1,B2)
      DENLOW  = (1. + MAX(0.,(FRSDEN-2.*RO(KK))*YPRO))
c ... test
      DENLOW  = (1. + MAX(0.,(FRSPRE-2.*P(KK))/P(KK)))

      DENMAX  = (1. + MAX(0.,(RO(KK)-2.*FRSDEN)*YPRO)*5.)
      DENLOW  = MAX(DENLOW,DENMAX)
c        denlow = 1.
      DTL(KK) = 1./(B11 * DENLOW)

 4000 CONTINUE

C ... FOR TIME-ACCURATE SOLUTIONS

      IF(TIMEL) THEN
         DO 5000 K = 1,KMAX
         KA      = (KN+K-1)*IL
         DO 5000 J = 1,JMAX
         JJ      = (JN+J-1)*ISTRID + IN + KA
         DO 5000 I = 1,IMAX
            KK      = JJ + I
            DTL(KK) = DTL(KK)*DT/(1.5*DTL(KK)+DT)
 5000    CONTINUE
      ENDIF ! TIMEL

      RETURN
      END SUBROUTINE TIME3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     + A3X,A3Y,A3Z,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDI1,IDI2,IDI3,M,
     + TIMEL,PVISC,ZZZ,IVISC,CDIFF)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A1(*),A1X(*),A1Y(*),A2(*),A2X(*),A2Y(*),DTL(*),P(*),
     + VOL(*),RO(*),U(*),V(*),W(*),C(*),VIS(*),VIST(*),
     + A1Z(*),A2Z(*),A3(*),A3X(*),A3Y(*),A3Z(*),D1(*),D2(*),D3(*),
     + UROT(*),VROT(*),WROT(*),ZZZ(*)
      LOGICAL :: TIMEL

      EPS     = 1.E-10
C
C ... TIME-STEP OF THE IMPLICIT BIDIAGONAL METHOD (3-D)
C
C     CDIFF   = GAMMA/(PR+EPS)*1.
C      CDIFF   = 200.
      PCFL    = 1./CFL
C     PVISC   = 1.8            !*8.**(M-1)
      T       = 0.
      PRS     = PR/PRT

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      NTOT    = ISTRID*JSTRID*KSTRID

      CALL ZEROZZ(ZZZ,NTOT)

C ... FIND MAXIMUN SIZE OF THE CELLS
      IF(IDI1 /= 0) THEN
      DO 1000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 1000 I = 1,MIN(IDI1,IMAX)
         KK      = JJ + I
         RNUM    = .5*(A1(KK) + A1(KK+1))
         ZZZ(KK) = RNUM
 1000 CONTINUE
      ENDIF

      IF(IDI2 /= 0) THEN
      DO 2000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 2000 J = 1,MIN(IDI2,JMAX)
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 2000 I = 1,IMAX
         KK      = JJ + I
         RNUM    = .5*(A2(KK) + A2(KK+ISTRID))
         ZZZ(KK) = MAX(RNUM,ZZZ(KK))
 2000 CONTINUE
      ENDIF

      IF(IDI3 /= 0) THEN
      DO 3000 K = 1,MIN(IDI3,KMAX)
      KA      = (KN+K-1)*IL
      DO 3000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 3000 I = 1,IMAX
         KK      = JJ + I
         RNUM    = .5*(A3(KK) + A3(KK+IL))
         ZZZ(KK) = MAX(RNUM,ZZZ(KK))
 3000 CONTINUE
      ENDIF

      DO 4000 K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO 4000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 4000 I = 1,IMAX
      KK      = JJ + I

      PVOL    = 1./VOL(KK)
      YPRO    = 1./(RO(KK)+ EPS)
C ... CONVECTIVE TIME STEP
      A1K     = .5*(A1(KK) + A1(KK+1))
      A2K     = .5*(A2(KK) + A2(KK+ISTRID))
      A3K     = .5*(A3(KK) + A3(KK+IL))
      A1XK    = .5*(A1X(KK)*A1(KK) + A1X(KK+1)     *A1(KK+1))
      A1YK    = .5*(A1Y(KK)*A1(KK) + A1Y(KK+1)     *A1(KK+1))
      A1ZK    = .5*(A1Z(KK)*A1(KK) + A1Z(KK+1)     *A1(KK+1))
      A2XK    = .5*(A2X(KK)*A2(KK) + A2X(KK+ISTRID)*A2(KK+ISTRID))
      A2YK    = .5*(A2Y(KK)*A2(KK) + A2Y(KK+ISTRID)*A2(KK+ISTRID))
      A2ZK    = .5*(A2Z(KK)*A2(KK) + A2Z(KK+ISTRID)*A2(KK+ISTRID))
      A3XK    = .5*(A3X(KK)*A3(KK) + A3X(KK+IL)    *A3(KK+IL))
      A3YK    = .5*(A3Y(KK)*A3(KK) + A3Y(KK+IL)    *A3(KK+IL))
      A3ZK    = .5*(A3Z(KK)*A3(KK) + A3Z(KK+IL)    *A3(KK+IL))
      UCROT   = .5*(A1(KK)*UROT(KK) + A1(KK+1)     *UROT(KK+1))
      VCROT   = .5*(A2(KK)*VROT(KK) + A2(KK+ISTRID)*VROT(KK+ISTRID))
      WCROT   = .5*(A3(KK)*WROT(KK) + A3(KK+IL)    *WROT(KK+IL))
      UCII    = A1XK*U(KK) + A1YK*V(KK) + A1ZK*W(KK) - UCROT
      VCII    = A2XK*U(KK) + A2YK*V(KK) + A2ZK*W(KK) - VCROT
      WCII    = A3XK*U(KK) + A3YK*V(KK) + A3ZK*W(KK) - WCROT
      IF(IVISC <= 2) THEN
      B1      = PVOL*(ABS(UCII) + ABS(VCII) + ABS(WCII)
     +        + ARTSSP*SQRT(A1K**2 + A2K**2 + A3K**2))*PCFL
      ELSE ! only convection
      B1      = PVOL*(ABS(UCII) + ABS(VCII) + ABS(WCII))*PCFL
      ENDIF

C ... VISCOUS TIME STEP
      RNUM    = ZZZ(KK)**2*PVOL**2*YPRO   ! 1/d**2
      T       = RNUM*(VIS(KK)+ PRS*VIST(KK)) / CDIFF !*.1, Changed
      B2      = 2.*T*PVISC
      IF(IVISC == 2 .OR. IVISC == 4) B2 = 0. ! No viscous time-step

C ... COMBINATIONS
CC    B11     = B1 + B2
      B11     = MAX(B1,B2)
      DENLOW  = (1. + MAX(0.,(FRSDEN-2.*RO(KK))*YPRO))
      DENMAX  = (1. + MAX(0.,(RO(KK)-2.*FRSDEN)*YPRO)*5.)
      DENLOW  = MAX(DENLOW,DENMAX)
c      DENLOW  = MIN(10.,DENLOW)
      DTL(KK) = 1./(B11 * DENLOW)                   *10. ! What is this?
c      IF(A1(KK) < 1.E-20 .OR. A1(KK+1) < 1.E-20      .OR.
c     +   A2(KK) < 1.E-20 .OR. A2(KK+ISTRID) < 1.E-20 .OR.
c     +   A3(KK) < 1.E-20 .OR. A3(KK+IL) < 1.E-20) DTL(KK) = DTL(KK)*.2

 4000 CONTINUE

C ... FOR TIME-ACCURATE SOLUTIONS REMOVED

c      IF(TIMEL) THEN
c         DO 5000 K = 1,KMAX
c         KA      = (KN+K-1)*IL
c         DO 5000 J = 1,JMAX
c         JJ      = (JN+J-1)*ISTRID + IN + KA
c         DO 5000 I = 1,IMAX
c            KK      = JJ + I
c            DTL(KK) = DTL(KK)/(1.5*DTL(KK)/DT+1.)
c 5000    CONTINUE
c      ENDIF ! TIMEL

      RETURN
      END SUBROUTINE TIMEIN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TIMEHS(DTL,CH,RO,CP,D1,D2,D3,SMAX,IMAX,JMAX,KMAX,
     +  IN,JN,KN)

      IMPLICIT NONE

      REAL    :: DTL(*), CH(*), RO(*), CP(*), D1(*), D2(*), D3(*), SMAX
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN
      INTEGER :: I, J, K, KA, KK, JJ, ISTRID, JSTRID, IL
      REAL    :: DIFF, DT1, DT2, DT3

C ... Loca time step for the heat structures

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
c     SMAX    = 20. ! Large enough?

      DO 1000 K = 0,KMAX+1
      KA      = (KN+K-1)*IL
      DO 1000 J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 1000 I = 0,IMAX+1
         KK     = JJ + I
         DIFF   = CH(KK)/(RO(KK)*CP(KK))
         DT1    = SMAX*D1(KK)**2/DIFF
         DT2    = SMAX*D2(KK)**2/DIFF
         DT3    = SMAX*D3(KK)**2/DIFF
         DTL(KK)= MIN(DT1,DT2,DT3)
 1000 CONTINUE
      
      RETURN
      END SUBROUTINE TIMEHS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PREUPM(RO,PRO,DRDP,DRDH,U,V,W,FRSPRE,FRSDEN,
     +           ARTSSP,FRSSSP,NTOT,NGL,ISTATE,IPRESC,PSEUCO,
     +           IMAX,JMAX,KMAX,IN,JN,KN,IPHASE,MULPHL,INN)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES) PRO(*)

      REAL    :: RO(*),U(*),V(*),W(*),DRDP(*),DRDH(*)
      REAL    :: FRSSSP,FRSSS2,BETA,PSEUCO,GAMMA,RGAS,FRSPRE,FRSDEN,
     +           ARTSSP,DROPE,ULOC,RMA2,UMAX
      INTEGER :: NTOT,NGL,ISTATE,IPRESC,I,IMAX,JMAX,KMAX,IN,JN,KN,
     +     ISTRID,JSTRID,IL,J,K,JJ,KA,KK,IPHASE,INN
      LOGICAL :: MULPHL
        
      FRSSS2  = FRSSSP**2
      BETA    = PSEUCO**2
      UMAX    = 0.
       
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      DO K =  1-INN,KMAX+INN
      KA      = (KN+K-1)*IL
      DO J =  1-INN,JMAX+INN
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO I =  1-INN,IMAX+INN
         KK   = JJ + I
         UMAX = MAX(UMAX,U(KK)**2 + V(KK)**2 + W(KK)**2) !  Find Umax
      ENDDO
      ENDDO
      ENDDO

      UMAX = 0. !25.*UMAX/PSEUCO**2  ! Limit the maximum to 0.2

      IF(MULPHL) THEN
         DO K =  1-INN,KMAX+INN
         KA      = (KN+K-1)*IL
         DO J =  1-INN,JMAX+INN
         JJ      = (JN+J-1)*ISTRID + IN + KA
         DO I =  1-INN,IMAX+INN
         KK   = JJ + I
            ULOC = .5*(U(KK)**2 + V(KK)**2 + W(KK)**2 )
            PRO(KK)%DRODP(IPHASE) =!PRO(KK)%DRODP(IPHASE) +! In MF does not work
     +      PRO(KK)%RO(IPHASE)/RO(KK)* 
     +      1./(BETA*MAX(ARTSSP**2/BETA,ULOC,UMAX)) - 
     +      PRO(KK)%DRODH(IPHASE)/PRO(KK)%RO(IPHASE)
         ENDDO
         ENDDO
         ENDDO
      ELSE
         DO K =  1-INN,KMAX+INN
         KA      = (KN+K-1)*IL
         DO J =  1-INN,JMAX+INN
         JJ      = (JN+J-1)*ISTRID + IN + KA
         DO I =  1-INN,IMAX+INN
         KK   = JJ + I

            ULOC = U(KK)**2 + V(KK)**2 + W(KK)**2 !  Velocities
c           RMA2 = ULOC/FRSSS2
c           IF(RMA2 <= 1./BETA) THEN
            DRDP(KK) = DRDP(KK) + !MAX(DRDP(KK),
     +      1./(BETA*MAX(ARTSSP**2/BETA,ULOC,UMAX)) - DRDH(KK)/RO(KK)!)
c     +     1./(BETA*MAX(ARTSSP**2/BETA,ULOC)) - DRDH(KK)/RO(KK)!)
c          if(ngl == 1) write(5555,*)KK,i,j,k,DRDP(KK),sqrt(uloc),beta
C laivakokeilu, lirutetaan alespäin.
c           DRDP(KK) = 1./14.**2
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE PREUPM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONSVA(RO,RM,RN,RW,E,C,U,V,W,RK,P,DRDP,DRDH,XC,YC,ZC,
     + PRO,VAR,FRSPRE,FRSVEL,NTOT,ITURB,JRIMP,FI,MAXSB,NSCAL,IPRESC,
     + PSEUCO,MULPHC,NPHASE,NGL,RNUT,REPS)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND

      IMPLICIT NONE

      INTEGER :: NTOT,ITURB,JRIMP,MAXSB,NSCAL,IPRESC,NPHASE,I,IPHASE,ngl
      REAL    :: RO(*),RM(*),RN(*),RW(*),E(*),
     +           U(*),C(*),V(*),W(*),RK(*),P(*),RNUT(*),REPS(*),
     +           DRDP(*),DRDH(*),FI(MAXSB,MAX(1,NSCAL))
      REAL    :: PSEUCO,CMIN,EPS,R23,RKIN,EPOT,YPRO,
     +           XTURB,FRSPRE,FRSVEL,ETOT
      REAL  ::   XC(*), YC(*), ZC(*)

      CHARACTER(LEN=10) ::  MULPHC

      TYPE (PROPERTIES)       :: PRO(*)
      TYPE (MPHASE_VARIABLES) :: VAR(*)

      CMIN = 1.
      EPS  = 1.E-10
      R23 = 2./3.
C
C ... SOUND SPEED BASED ON THE DERIVATIVES OF DENSITY AND TOTAL
C ... ENERGY FROM THE SPECIFIC INTERNAL ENERGY
C
      
      IF(MULPHC /= 'MULTI' .OR.
     +   MULPHC /= 'CAVIT') THEN ! A single-phase flow

      IF(ITURB < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL

      DO 1200 I = 1,NTOT
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*.5
         C(I)   = SQRT(MAX(1./(DRDH(I)/RO(I)+DRDP(I)),CMIN))
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +             +GZ*(ZC(I)-GROUND))
         EPOT   = 0.
         E(I)   = (E(I) + RKIN + EPOT)*RO(I)
         RM(I)  = RO(I)*U(I)
         RN(I)  = RO(I)*V(I)
         RW(I)  = RO(I)*W(I)
1200  CONTINUE

      ELSEIF(ITURB >= 3 .AND. JRIMP == 0 .AND. ITURB /= 8) THEN  ! TWO-EQUATION
 
      DO 1000 I = 1,NTOT
         YPRO   = 1./RO(I)
         XTURB  = 1. + R23*RK(I)*YPRO*DRDP(I)
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*.5
         C(I)   = SQRT(MAX(XTURB/(DRDH(I)*YPRO + DRDP(I)),
     +            CMIN))
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +             +GZ*(ZC(I)-GROUND))
         EPOT   = 0.
         
         E(I)   = (E(I) + RKIN + EPOT)*RO(I) + RK(I)
         RM(I)  = RO(I)*U(I)
         RN(I)  = RO(I)*V(I)
         RW(I)  = RO(I)*W(I)
         IF(ITURB == 9) RNUT(I) = REPS(I)/RO(I)
1000  CONTINUE

      ELSEIF(ITURB >= 24 .AND. JRIMP == 1) THEN ! REYNOLDS STRESS

      DO 1100 I = 1,NTOT
         YPRO   = 1./RO(I)
         XTURB  = 1. + R23*RK(I)*YPRO*DRDP(I)
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*.5
         C(I)   = SQRT(MAX(XTURB/(DRDH(I)*YPRO + DRDP(I)),
     +            CMIN))
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +             +GZ*(ZC(I)-GROUND))
         EPOT   = 0.
         E(I)   = (E(I)+RKIN+EPOT)*RO(I)+.5*(FI(I,1)+FI(I,4)+FI(I,6))
         RM(I)  = RO(I)*U(I)
         RN(I)  = RO(I)*V(I)
         RW(I)  = RO(I)*W(I)
1100  CONTINUE
      ELSE
         WRITE(*,*) 'impossible jrimp and iturb values'
         WRITE(*,*) 'Exiting...'
         STOP
      ENDIF
      ENDIF ! Single-phase flow

      IF(MULPHC == 'MULTI'  .OR.
     +   MULPHC == 'CAVIT') THEN ! A multi-phase flow
         DO I = 1,NTOT
         E(I)   = 0.  !RK(I)
         C(I)   = 0.
         RM(I)  = 0. ; RN(I) = 0. ; RW(I) = 0.
         DO IPHASE = 1,NPHASE
         YPRO   = 1./PRO(I)%RO(IPHASE)
         XTURB  = 1. + R23*RK(I)*YPRO*PRO(I)%DRODP(IPHASE)
         xturb  = 1.

         PRO(I)%C(IPHASE)   = SQRT(MAX(XTURB/(PRO(I)%DRODH(IPHASE)*YPRO 
     +           + PRO(I)%DRODP(IPHASE)),CMIN))
c         VAR(I)%ARO(IPHASE) = VAR(I)%ALFA(IPHASE)*PRO(I)%RO(IPHASE)
         IF(MULPHC == 'MULTI') THEN
            VAR(I)%ARM(IPHASE) = VAR(I)%ARO(IPHASE)*VAR(I)%U(IPHASE)
            VAR(I)%ARN(IPHASE) = VAR(I)%ARO(IPHASE)*VAR(I)%V(IPHASE)
            VAR(I)%ARW(IPHASE) = VAR(I)%ARO(IPHASE)*VAR(I)%W(IPHASE)
            RM(I)  = RM(I) + VAR(I)%ARM(IPHASE)
            RN(I)  = RN(I) + VAR(I)%ARN(IPHASE)
            RW(I)  = RW(I) + VAR(I)%ARW(IPHASE)
            RKIN   = (VAR(I)%U(IPHASE)**2 + VAR(I)%V(IPHASE)**2 + 
     +                VAR(I)%W(IPHASE)**2)*.5  + RK(I)/RO(I)

            ETOT   = (PRO(I)%E(IPHASE) + RKIN) * VAR(I)%ARO(IPHASE)
            VAR(I)%ARE(IPHASE)  = ETOT
            VAR(I)%ETOT(IPHASE) = PRO(I)%E(IPHASE) + RKIN
         ELSEIF(MULPHC == 'CAVIT') THEN
            RM(I)  = RO(I)*U(I)
            RN(I)  = RO(I)*V(I)
            RW(I)  = RO(I)*W(I)
            RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*.5 + RK(I)/RO(I)
            ETOT   = (PRO(I)%E(IPHASE) + RKIN) * VAR(I)%ARO(IPHASE)
         VAR(I)%ARE(IPHASE)  = ETOT
         VAR(I)%ETOT(IPHASE) = PRO(I)%E(IPHASE) + RKIN
         ENDIF
         E(I)   = E(I) + ETOT
         C(I)   = C(I) + VAR(I)%ALFA(IPHASE)/(PRO(I)%RO(IPHASE) *
     +            PRO(I)%C(IPHASE)**2)
         ENDDO
         C(I)   = 1./SQRT(RO(I)*C(I))
c         C(I)   = 25.
      ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CONSVA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRIMVE(U,V,W,RK,RO,E,XC,YC,ZC,NTOT,FRSDEN,ITURB)

      USE NS3CO, ONLY : GX, GY, GZ,GROUND

      REAL :: RO(*), E(*), U(*), V(*), W(*), RK(*)
      REAL :: XC(*), YC(*), ZC(*)
C
C ... INTERNAL ENERGY FROM TOTAL ENERGY
C
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! TWO-EQUATION MODEL
      DO 1000 I = 1,NTOT
         RO(I)  = MAX(1.E-4*FRSDEN,RO(I))
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*RO(I)*.5
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +            +GZ*(ZC(I)-GROUND))*RO(I)
         EPOT   = 0.
         E(I)   = (E(I) - RKIN - RK(I) - EPOT)/RO(I)
1000  CONTINUE
      ENDIF

      IF(ITURB < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL 
      DO 1200 I = 1,NTOT
         RO(I)  = MAX(1.E-4*FRSDEN,RO(I))
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*RO(I)*.5
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +            +GZ*(ZC(I)-GROUND))*RO(I)
         EPOT   = 0.
         E(I)   = (E(I) - RKIN - EPOT)/RO(I)
1200  CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE PRIMVE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONSVE(E,U,V,W,RK,RO,XC,YC,ZC,
     &                  NTOT,ITURB,FI,MAXSB,NSCAL)

      USE NS3CO, ONLY : GX, GY, GZ, GROUND

      REAL :: RO(*),E(*),U(*),V(*),W(*),RK(*),FI(MAXSB,MAX(1,NSCAL))
      REAL :: XC(*), YC(*), ZC(*)

      EPS     = 1.E-20
C
C ... TOTAL ENERGY FROM SPECIFIC INTERNAL ENERGY
C
      IF(ITURB < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL
      DO 1200 I = 1,NTOT
         RKIN   = (U(I)**2 + V(I)**2 + W(I)**2)*.5
c         EPOT   =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +            +GZ*(ZC(I)-GROUND))
         EPOT   = 0.
         E(I)   = (E(I) + RKIN + EPOT)*RO(I)
1200  CONTINUE

      ELSEIF(ITURB >= 3 .AND. ITURB /= 8) THEN ! TWO-EQUATION MODEL
      DO 1000 I = 1,NTOT
      RKIN    = (U(I)**2 + V(I)**2 + W(I)**2)*.5
c      EPOT    =-(GX*(XC(I)-GROUND) + GY*(YC(I)-GROUND)
c     +          +GZ*(ZC(I)-GROUND))
      EPOT     = 0.
      E(I)    = (E(I) + RKIN + EPOT)*RO(I) + RK(I)
1000  CONTINUE

      ENDIF

      RETURN
      END SUBROUTINE CONSVE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DERIV(RO,P,DEDPH,DEDRH,DRDP,DRDH,EIP,EIM,ELR,ERL,
     +     PPI,PMI,ROIP,ROIM,NTOT,IL)

      REAL :: RO(*),P(*),DEDPH(*),DEDRH(*),DRDP(*),DRDH(*),EIP(*),
     +     EIM(*),ELR(*),ERL(*),PPI(*),PMI(*),ROIP(*),ROIM(*)
C
C ... HAT DERIVATIVES OF INTERNAL ENERGY ACCORDING TO GLAISTER
C
      DO 10 I = 1,NTOT
         I2   = I - IL
         IF(I2 <= 0) I2 = I ! On the lower boundary
         YPRO1 = 1./RO(I)
         YPRO2 = 1./RO(I2)
         DEDP1   = -DRDP(I )/DRDH(I ) - YPRO1
         DEDP2   = -DRDP(I2)/DRDH(I2) - YPRO2
         DEDPH(I) = 0.5*(DEDP1 + DEDP2)
         
         DEDR1   = 1./DRDH(I ) + P(I )*YPRO1**2
         DEDR2   = 1./DRDH(I2) + P(I2)*YPRO2**2
         DEDRH(I) = 0.5*(DEDR1+DEDR2)
 10   CONTINUE
      DO 20 I = 1,NTOT
         DP     = 2.*(PPI(I) - PMI(I))
         DRO    = 2.*(ROIP(I)- ROIM(I))

         IF(ABS(DP) >= 1.E-2*(PPI(I) + PMI(I))) THEN
            DEDPH(I) = (EIP(I)+ELR(I) - EIM(I)-ERL(I))/DP
         ENDIF
      
         IF(ABS(DRO) >= 1.E-2*(ROIP(I)+ROIM(I))) THEN
            DEDRH(I) = (EIP(I)+ERL(I) - EIM(I)-ELR(I))/DRO
         ENDIF
 20   CONTINUE

      RETURN
      END SUBROUTINE DERIV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ******************* LIMITERS AND FLUXES ******************************
      SUBROUTINE XXTRAP(RO,ROIP,ROIM,NTOT,IL,KL,IA,INTEM,RKSI,INCHIML)

      REAL    :: RO(*), ROIP(*), ROIM(*), RKSI(*)
      LOGICAL :: INCHIML
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C
C ...    RIGHT SIDE:  ROIP
C ...    LEFT SIDE :  ROIM

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-10
      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION

      DO  400 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL

         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(I)  
            ROIM(IG)= RO(IT1)
         ELSE
            DIFF0   = RO(IT1) - RO(IT0)
            DIFF1   = RO(I)   - RO(IT1)
            DIFF2   = RO(IT2) - RO(I)
            ROIP(IG)= RO(I)   - (DIFF1*R1 + DIFF2*R2)
            ROIM(IG)= RO(IT1) + (DIFF1*R1 + DIFF0*R2)
       ENDIF

 400  CONTINUE

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO  500 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(I)  
            ROIM(IG)= RO(IT1)
c            ROIP(IG)= .5*(1.95*RO(I)+.05*RO(IT1))   ! RO(I)  
c            ROIM(IG)= .5*(.05*RO(I)+1.95*RO(IT1))   ! RO(IT1)

         ELSE
            DIFF0   = RO(IT1) - RO(IT0)
            DIFF1   = RO(I)   - RO(IT1)
            DIFF2   = RO(IT2) - RO(I)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            ROIP(IG)= RO(I)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            ROIM(IG)= RO(IT1) + PHIBI*(DIFF1*R1 + DIFF0*R2)

         ENDIF

 500  CONTINUE

      ELSE ! FIRST ORDER UPWIND

         DO IG = 1,NTOT
         I       = IA + IG
         ROIP(IG)= RO(I)
         ROIM(IG)= RO(I-IL)
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE XXTRAP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE XXTRAP_SURF(RO,ROIP,ROIM,KX1,KX2,KY1,KY2,ISTR,JSTR,
     +      IN,JN,IL,KA,INTEM,RKSI,INCHIML,NFL)

      IMPLICIT NONE

      INTEGER :: INTER,INTEM,INTE2,KX1,KX2,KY1,KY2,NN,I,J,IJ,IN,JN,
     +           KA,JSTR,ISTR,IT0,IT1,IT2,IG,IL,II,NFL
      REAL    :: RK1,RK2,R1,R2,EPS,DIFF0,DIFF1,DIFF2,DIFF12,PHIBI,PHIFI

      REAL    :: RO(*), ROIP(*), ROIM(*), RKSI(*)
      LOGICAL :: INCHIML
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C
C ...    RIGHT SIDE:  ROIP
C ...    LEFT SIDE :  ROIM

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL :: RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-10
      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL                 ! 2. comp cell
         IT1     = II - IL                 ! First ghost cell
         IT0     = II - 2*IL               ! Second ghost cell
         IG      = I  + NN                ! Patch interior index
*

         IF((RKSI(II) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(II) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(II)  
            ROIM(IG)= RO(IT1)
         ELSE
            DIFF0   = RO(IT1) - RO(IT0)
            DIFF1   = RO(II)   - RO(IT1)
            DIFF2   = RO(IT2) - RO(II)
            ROIP(IG)= RO(II)   - (DIFF1*R1 + DIFF2*R2)
            ROIM(IG)= RO(IT1) + (DIFF1*R1 + DIFF0*R2)
       ENDIF

      ENDDO ; ENDDO

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL            ! 2. comp cell
         IT1     = II - IL            ! First ghost cell
         IT0     = II - 2*IL          ! Second ghost cell
         IG      = I  + NN               ! Patch interior index

         IF((RKSI(II) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(II) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(II)  
            ROIM(IG)= RO(IT1)
         ELSE
            DIFF0   = RO(IT1) - RO(IT0)
            DIFF1   = RO(II)   - RO(IT1)
            DIFF2   = RO(IT2) - RO(II)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            ROIP(IG)= RO(II)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            ROIM(IG)= RO(IT1) + PHIBI*(DIFF1*R1 + DIFF0*R2)

         ENDIF

      ENDDO ; ENDDO

      ELSE ! FIRST ORDER UPWIND

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL            ! 2. comp cell
         IT1     = II - IL            ! First ghost cell
         IT0     = II - 2*IL          ! Second ghost cell
         IG      = I  + NN               ! Patch interior index
         ROIP(IG)= RO(II)
         ROIM(IG)= RO(II-IL)
      ENDDO ; ENDDO

      ENDIF

      RETURN
      END SUBROUTINE XXTRAP_SURF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TXTRAP(TRM,G1R,G1L,RET1R,RET1L,NTOT,IL,KL,IA,INTEM,
     + RKSI,INCHIML)

      USE TYPE_ARRAYS

      INTEGER :: NTOT, IL, KL, IA, INTEM
      REAL    :: G1R(*), G1L(*), RET1R(*), RET1L(*), RKSI(*)
      LOGICAL :: INCHIML

      TYPE(INTERMITTENCY) :: TRM(*)

C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C ... TRANSITION MODEL VARIABLES
C
C ...    RIGHT SIDE:  G1R, RET1R
C ...    LEFT SIDE :  G1L, RET1L

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-10
      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION

      DO  400 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
C ... Does this make sense on the Chimera boundary?? (not tested now)

            G1R(IG)   = TRM(I)%G
            G1L(IG)   = TRM(IT1)%G
            RET1R(IG) = TRM(I)%RET
            RET1L(IG) = TRM(IT1)%RET

         ELSE

            DIFF0     = TRM(IT1)%G - TRM(IT0)%G
            DIFF1     = TRM(I)%G   - TRM(IT1)%G
            DIFF2     = TRM(IT2)%G - TRM(I)%G
            G1R(IG)   = TRM(I)%G   - (DIFF1*R1 + DIFF2*R2)
            G1L(IG)   = TRM(IT1)%G + (DIFF1*R1 + DIFF0*R2)

            DIFF0     = TRM(IT1)%RET - TRM(IT0)%RET
            DIFF1     = TRM(I)%RET   - TRM(IT1)%RET
            DIFF2     = TRM(IT2)%RET - TRM(I)%G
            RET1R(IG) = TRM(I)%RET   - (DIFF1*R1 + DIFF2*R2)
            RET1L(IG) = TRM(IT1)%RET + (DIFF1*R1 + DIFF0*R2)

         ENDIF
 400  CONTINUE

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO  500 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
C ... Does this make sense on the Chimera boundary?? (not tested now)

            G1R(IG)   = TRM(I)%G
            G1L(IG)   = TRM(IT1)%G
            RET1R(IG) = TRM(I)%RET
            RET1L(IG) = TRM(IT1)%RET

         ELSE

            DIFF0   = TRM(IT1)%G - TRM(IT0)%G
            DIFF1   = TRM(I)%G   - TRM(IT1)%G
            DIFF2   = TRM(IT2)%G - TRM(I)%G
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            G1R(IG)= TRM(I)%G   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            G1L(IG)= TRM(IT1)%G + PHIBI*(DIFF1*R1 + DIFF0*R2)

            DIFF0   = TRM(IT1)%RET - TRM(IT0)%RET
            DIFF1   = TRM(I)%RET   - TRM(IT1)%RET
            DIFF2   = TRM(IT2)%RET - TRM(I)%RET
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            RET1R(IG)= TRM(I)%RET   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            RET1L(IG)= TRM(IT1)%RET + PHIBI*(DIFF1*R1 + DIFF0*R2)
      ENDIF
 500  CONTINUE

      ELSE ! FIRST ORDER UPWIND
      DO  600 IG = 1,NTOT
         I       = IA + IG
         G1R(IG) = TRM(I)%G
         G1L(IG) = TRM(I-IL)%G

         RET1R(IG) = TRM(I)%RET
         RET1L(IG) = TRM(I-IL)%RET
 600  CONTINUE


      ENDIF

      RETURN
      END SUBROUTINE TXTRAP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DXTRAP(PRO,VAR,T1R,T1L,A1R,A1L,HSATR,HSATL,CPR,CPL,
     + TSATR,TSATL,NTOT,IL,KL,IA,INTEM,IPHASE,RKSI,INCHIML,MULPHC,INTEA)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL    :: T1R(*), T1L(*), A1R(*), A1L(*), RKSI(*),
     +           HSATR(*),HSATL(*),CPR(*),CPL(*),TSATR(*),TSATL(*)
      INTEGER :: IPHASE,INTEA,NTOT,IL,KL,IA,INTEM,INTE2,INTER,IG,I,IT0,
     +           IT1,IT2
      LOGICAL :: INCHIML
      CHARACTER(LEN=*) :: MULPHC
      REAL :: RK1,RK2,R1,R2,EPS,DIFF0,DIFF1,DIFF2,DIFF12,PHIBI,PHIFI,
     +        RIII,RII2,BETA

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)

C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL (MULTIPHASE)
C
C ...    RIGHT SIDE:  T1R, ROIP
C ...    LEFT SIDE :  T1L, ROIM

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

C ... vmv: pultattiin tanne ennen ykkosena
      INTER   = IABS(INTEM)

C ... vmv: 
ccc      INTER   = 3

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-10
      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION

      DO  400 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            T1R(IG) = PRO(I)%DTEMP(IPHASE)
            T1L(IG) = PRO(I-IL)%DTEMP(IPHASE)
            HSATR(IG) = PRO(I)%HSAT(IPHASE)
            HSATL(IG) = PRO(I-IL)%HSAT(IPHASE)
            CPR(IG)   = PRO(I)%CP(IPHASE)
            CPL(IG)   = PRO(I-IL)%CP(IPHASE)
            TSATR(IG) = PRO(I)%TSAT
            TSATL(IG) = PRO(I-IL)%TSAT
         ELSE
            DIFF0   = PRO(IT1)%DTEMP(IPHASE) - PRO(IT0)%DTEMP(IPHASE)
            DIFF1   = PRO(I)%DTEMP(IPHASE)   - PRO(IT1)%DTEMP(IPHASE)
            DIFF2   = PRO(IT2)%DTEMP(IPHASE) - PRO(I)%DTEMP(IPHASE)
            T1R(IG) = PRO(I)%DTEMP(IPHASE)   - (DIFF1*R1 + DIFF2*R2)
            T1L(IG) = PRO(IT1)%DTEMP(IPHASE) + (DIFF1*R1 + DIFF0*R2)

            DIFF0   = PRO(IT1)%HSAT(IPHASE) - PRO(IT0)%HSAT(IPHASE)
            DIFF1   = PRO(I)%HSAT(IPHASE)   - PRO(IT1)%HSAT(IPHASE)
            DIFF2   = PRO(IT2)%HSAT(IPHASE) - PRO(I)%HSAT(IPHASE)
            HSATR(IG) = PRO(I)%HSAT(IPHASE)   - (DIFF1*R1 + DIFF2*R2)
            HSATL(IG) = PRO(IT1)%HSAT(IPHASE) + (DIFF1*R1 + DIFF0*R2)

            DIFF0   = PRO(IT1)%CP(IPHASE) - PRO(IT0)%CP(IPHASE)
            DIFF1   = PRO(I)%CP(IPHASE)   - PRO(IT1)%CP(IPHASE)
            DIFF2   = PRO(IT2)%CP(IPHASE) - PRO(I)%CP(IPHASE)
            CPR(IG) = PRO(I)%CP(IPHASE)   - (DIFF1*R1 + DIFF2*R2)
            CPL(IG) = PRO(IT1)%CP(IPHASE) + (DIFF1*R1 + DIFF0*R2)

            DIFF0   = PRO(IT1)%TSAT - PRO(IT0)%TSAT
            DIFF1   = PRO(I)%TSAT   - PRO(IT1)%TSAT
            DIFF2   = PRO(IT2)%TSAT - PRO(I)%TSAT
            TSATR(IG) = PRO(I)%TSAT   - (DIFF1*R1 + DIFF2*R2)
            TSATL(IG) = PRO(IT1)%TSAT + (DIFF1*R1 + DIFF0*R2)
         ENDIF
 400  CONTINUE

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO  500 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
           T1R(IG)   = PRO(I)%DTEMP(IPHASE)
           T1L(IG)   = PRO(I-IL)%DTEMP(IPHASE)
           HSATR(IG) = PRO(I)%HSAT(IPHASE)
           HSATL(IG) = PRO(I-IL)%HSAT(IPHASE)
           CPR(IG)   = PRO(I)%CP(IPHASE)
           CPL(IG)   = PRO(I-IL)%CP(IPHASE)
           TSATR(IG) = PRO(I)%TSAT
           TSATL(IG) = PRO(I-IL)%TSAT
         ELSE
            DIFF0   = PRO(IT1)%DTEMP(IPHASE) - PRO(IT0)%DTEMP(IPHASE)
            DIFF1   = PRO(I)%DTEMP(IPHASE)   - PRO(IT1)%DTEMP(IPHASE)
            DIFF2   = PRO(IT2)%DTEMP(IPHASE) - PRO(I)%DTEMP(IPHASE)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (ONLY VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            T1R(IG)= PRO(I)%DTEMP(IPHASE)   - PHIFI*(DIFF1*R1+ DIFF2*R2)
            T1L(IG)= PRO(IT1)%DTEMP(IPHASE) + PHIBI*(DIFF1*R1+ DIFF0*R2)

            DIFF0   = PRO(IT1)%HSAT(IPHASE) - PRO(IT0)%HSAT(IPHASE)
            DIFF1   = PRO(I)%HSAT(IPHASE)   - PRO(IT1)%HSAT(IPHASE)
            DIFF2   = PRO(IT2)%HSAT(IPHASE) - PRO(I)%HSAT(IPHASE)
            DIFF12  = DIFF1**2 + EPS
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
            HSATR(IG)= PRO(I)%HSAT(IPHASE)   - PHIFI*(DIFF1*R1+DIFF2*R2)
            HSATL(IG)= PRO(IT1)%HSAT(IPHASE) + PHIBI*(DIFF1*R1+DIFF0*R2)

            DIFF0    = PRO(IT1)%CP(IPHASE) - PRO(IT0)%CP(IPHASE)
            DIFF1    = PRO(I)%CP(IPHASE)   - PRO(IT1)%CP(IPHASE)
            DIFF2    = PRO(IT2)%CP(IPHASE) - PRO(I)%CP(IPHASE)
            DIFF12   = DIFF1**2 + EPS
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
            CPR(IG)  = PRO(I)%CP(IPHASE)   - PHIFI*(DIFF1*R1+DIFF2*R2)
            CPL(IG)  = PRO(IT1)%CP(IPHASE) + PHIBI*(DIFF1*R1+DIFF0*R2)

            DIFF0    = PRO(IT1)%TSAT - PRO(IT0)%TSAT
            DIFF1    = PRO(I)%TSAT   - PRO(IT1)%TSAT
            DIFF2    = PRO(IT2)%TSAT - PRO(I)%TSAT
            DIFF12   = DIFF1**2 + EPS
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
            TSATR(IG)= PRO(I)%TSAT   - PHIFI*(DIFF1*R1+DIFF2*R2)
            TSATL(IG)= PRO(IT1)%TSAT + PHIBI*(DIFF1*R1+DIFF0*R2)
         ENDIF
 500  CONTINUE

      ELSE ! FIRST ORDER UPWIND
      DO  600 IG = 1,NTOT
         I       = IA + IG
         T1R(IG)   = PRO(I)%DTEMP(IPHASE)
         T1L(IG)   = PRO(I-IL)%DTEMP(IPHASE)
         HSATR(IG) = PRO(I)%HSAT(IPHASE)
         HSATL(IG) = PRO(I-IL)%HSAT(IPHASE)
         CPR(IG)   = PRO(I)%CP(IPHASE)
         CPL(IG)   = PRO(I-IL)%CP(IPHASE)
         TSATR(IG) = PRO(I)%TSAT
         TSATL(IG) = PRO(I-IL)%TSAT

 600  CONTINUE

      ENDIF

C ... vmv: pultattiin tanne ennen ykkosena, nyt inputista

      INTER   = IABS(INTEA)

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-4 ! More stable than 1.E-6
      IF (INTEA <= 0 .AND. INTEA >= -5) THEN ! WITHOUT LIMITATION

      DO  1400 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN

            IF(IPHASE == 1) THEN
            A1R(IG) = VAR(I)%ALFA(IPHASE) 
            A1L(IG) = VAR(IT1)%ALFA(IPHASE) 
            ELSE IF(IPHASE == 2) THEN
            A1R(IG) = 1. - A1R(IG)
            A1L(IG) = 1. - A1L(IG)
            ENDIF
         ELSE
            DIFF0   = VAR(IT1)%ALFA(IPHASE) - VAR(IT0)%ALFA(IPHASE)
            DIFF1   = VAR(I)%ALFA(IPHASE)   - VAR(IT1)%ALFA(IPHASE)
            DIFF2   = VAR(IT2)%ALFA(IPHASE) - VAR(I)%ALFA(IPHASE)
            IF(IPHASE == 1) THEN
            A1R(IG) = VAR(I)%ALFA(IPHASE)   - (DIFF1*R1 + DIFF2*R2)
            A1L(IG) = VAR(IT1)%ALFA(IPHASE) + (DIFF1*R1 + DIFF0*R2)
            ELSE IF(IPHASE == 2) THEN
            A1R(IG) = 1. - A1R(IG)
            A1L(IG) = 1. - A1L(IG)
         ENDIF
      ENDIF
 1400  CONTINUE

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO  1500 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN

            IF(IPHASE == 1) THEN
            A1R(IG) = VAR(I)%ALFA(IPHASE) 
            A1L(IG) = VAR(IT1)%ALFA(IPHASE) 
            ELSE IF(IPHASE == 2) THEN
            A1R(IG) = 1. - A1R(IG)
            A1L(IG) = 1. - A1L(IG)
            ENDIF

         ELSE

            DIFF0   = VAR(IT1)%ALFA(IPHASE) - VAR(IT0)%ALFA(IPHASE)
            DIFF1   = VAR(I)%ALFA(IPHASE)   - VAR(IT1)%ALFA(IPHASE)
            DIFF2   = VAR(IT2)%ALFA(IPHASE) - VAR(I)%ALFA(IPHASE)
            DIFF12  = DIFF1**2 + EPS ! Always positive

            SELECT CASE(INTER)
C ... LIMITERS (VAN ALBADA)
            CASE(1)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
            PHIBI    = MAX(0.,PHIBI) ; PHIFI = MAX(0.,PHIFI)

C ... LIMITERS (VAN ALBADA 2, not TVD)
            CASE(2)
            PHIBI    = (2.*DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (2.*DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
            R1       = 0.; R2 = 0.5            
            
C ... vmv: limiter Roe SUPERBEE if BETA = 2 acting on a central difference
            CASE(3)
            BETA    = 2.0          ! Superbee. If 1 < BETA < 2 Sweby
            DIFF1   = SIGN(MAX(ABS(DIFF1),EPS),DIFF1) !DIFF1+EPS
            RIII    = DIFF0/DIFF1  ! r-
            RII2    = DIFF2/DIFF1  ! r+
            PHIBI   = MAX(0.0,MIN(BETA*RIII,1.0),MIN(RIII,BETA))
            PHIFI   = MAX(0.0,MIN(BETA*RII2,1.0),MIN(RII2,BETA))
c            R1      = 0.5 ; R2 = 0. ! Comment this to retain the original
             
C ... Limiter Sweby or the Roe SUPERBEE acting on a second-order upwind
            CASE(4)
            BETA    = 1.5          ! If 1 < BETA < 2 Sweby limiter
c            DIFF1   = DIFF1+EPS
c            RIII    = DIFF0/DIFF1  ! r- Central
c            RII2    = DIFF2/DIFF1  ! r+
c            RIII    = 1./(RIII + EPS) ! Upwind, choice 1
c            RII2    = 1./(RII2 + EPS)
            DIFF0   = SIGN(MAX(ABS(DIFF0),EPS),DIFF0) !DIFF0+EPS
            DIFF2   = SIGN(MAX(ABS(DIFF2),EPS),DIFF2) !DIFF2+EPS
            RIII    = DIFF1/DIFF0  ! r- Upwind, choice 2 (mersu)
            RII2    = DIFF1/DIFF2  ! r+

            PHIBI   = MAX(0.0,MIN(BETA*RIII,1.0),MIN(RIII,BETA))
            PHIFI   = MAX(0.0,MIN(BETA*RII2,1.0),MIN(RII2,BETA))
            R1      = 0.; R2 = 0.5
            END SELECT
            
C ... EXTRAPOLATING
            IF(IPHASE == 1) THEN
            A1R(IG)= VAR(I)%ALFA(IPHASE)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            A1L(IG)= VAR(IT1)%ALFA(IPHASE) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            ELSE IF(IPHASE == 2) THEN
            A1R(IG) = 1. - A1R(IG)
            A1L(IG) = 1. - A1L(IG)
            ENDIF
c           END SELECT

      ENDIF
 1500  CONTINUE

      ELSE ! FIRST ORDER UPWIND
      DO  1600 IG = 1,NTOT
         I       = IA + IG
         A1R(IG) = VAR(I)%ALFA(IPHASE)
         A1L(IG) = VAR(I-IL)%ALFA(IPHASE)
c         IF(IPHASE == 2) THEN
c         A1R(IG) = 1. - A1R(IG)
c         A1L(IG) = 1. - A1L(IG)
c         ENDIF
 1600  CONTINUE

      ENDIF

      RETURN
      END SUBROUTINE DXTRAP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE XXTRAA(RO,ROIP,ROIM,NTOT,IL,KL,IA,INTEM,A,VOL,D2,RKSI,
     + INCHIML)

      REAL    :: RO(*), ROIP(*), ROIM(*), A(*), VOL(*), D2(*), RKSI(*)
      LOGICAL :: INCHIML
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C ... size of the zells is included.
C ...    RIGHT SIDE:  ROIP
C ...    LEFT SIDE :  ROIM
C ... IT0 IT1   I  IT2
C ... -1   I   +1   +2
      REAL RKK1(5)

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      EPS     = 1.E-10
      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      DATA RKK1/0.,1.,1.3333333,2.,1.5/
      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1


      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION
      DO  400 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN

            ROIP(IG)= RO(I)  
            ROIM(IG)= RO(IT1)

         ELSE

c           DXL     = VOL(IT1)/(EPS+A(IT1)+A(I))
c           DXR     = VOL(I  )/(EPS+A(IT1)+A(IT2))
            DXL     = .5*D2(IT1)
            DXR     = .5*D2(I)

c           DX0     = A(IT1)/(VOL(IT0) + VOL(IT1))
c           DX1     = A(I)  /(VOL(IT1) + VOL(I  ))
c           DX2     = A(IT2)/(VOL(I  ) + VOL(IT2))
     
            DX0     = 1./(D2(IT0) + D2(IT1) + EPS)
            DX1     = 1./(D2(IT1) + D2(I)   + EPS)
            DX2     = 1./(D2(I)   + D2(IT2) + EPS)

            ROIP(IG)= RO(I)   - DXR*(DX1*(RO(I  )-RO(IT1))*RK1 +
     2                               DX2*(RO(IT2)-RO(I  ))*RK2)
            ROIM(IG)= RO(IT1) + DXL*(DX1*(RO(I) - RO(IT1))*RK1 +
     2                               DX0*(RO(IT1)-RO(IT0))*RK2)
      ENDIF
 400  CONTINUE

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO  500 IG = 1,NTOT
         I       = IA + IG
         IT0     = I - KL
         IT1     = I - IL
         IT2     = I + IL
         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN

            ROIP(IG)= RO(I)  
            ROIM(IG)= RO(IT1)

         ELSE


c           DXL     = VOL(IT1)/(EPS+A(IT1)+A(I))
c           DXR     = VOL(I  )/(EPS+A(IT1)+A(IT2))
            DXL     = .5*D2(IT1)
            DXR     = .5*D2(I)

c           DX0     = A(IT1)/(VOL(IT0) + VOL(IT1))
c           DX1     = A(I)  /(VOL(IT1) + VOL(I  ))
c           DX2     = A(IT2)/(VOL(I  ) + VOL(IT2))

            DX0     = 1./(D2(IT0) + D2(IT1) + EPS)
            DX1     = 1./(D2(IT1) + D2(I)   + EPS)
            DX2     = 1./(D2(I)   + D2(IT2) + EPS)

            DIFF0   = DX0*(RO(IT1) - RO(IT0))
            DIFF1   = DX1*(RO(I)   - RO(IT1))
            DIFF2   = DX2*(RO(IT2) - RO(I))
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)

C ... EXTRAPOLATING
            ROIP(IG)= RO(I)  -DXR*PHIFI*(DIFF1*RK1+DIFF2*RK2)
            ROIM(IG)= RO(IT1)+DXL*PHIBI*(DIFF1*RK1+DIFF0*RK2)
      ENDIF
 500  CONTINUE

      ELSE ! FIRST ORDER UPWIND

      DO  600 IG = 1,NTOT
         I       = IA + IG
         ROIP(IG)= RO(I)
         ROIM(IG)= RO(I-IL)
 600  CONTINUE

      ENDIF

      RETURN
      END SUBROUTINE XXTRAA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE XXTRAA_SURF(RO,ROIP,ROIM,KX1,KX2,KY1,KY2,ISTR,JSTR,
     +      IN,JN,IL,KA,INTEM,A,VOL,D2,RKSI,INCHIML,NFL)

      IMPLICIT NONE

      INTEGER :: INTER,INTEM,INTE2,KX1,KX2,KY1,KY2,NN,I,J,IJ,IN,JN,
     +           KA,JSTR,ISTR,IT0,IT1,IT2,IG,IL,II,NFL
      REAL    :: RK1,RK2,R1,R2,EPS,DIFF0,DIFF1,DIFF2,DIFF12,PHIBI,PHIFI,
     +     DXL,DXR,DX0,DX1,DX2
      REAL    :: RO(*),ROIP(*),ROIM(*),RKSI(*),A(*),VOL(*),D2(*)
      LOGICAL :: INCHIML
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C
C ...    RIGHT SIDE:  ROIP
C ...    LEFT SIDE :  ROIM

C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2
      EPS     = 1.E-10
      IF (INTEM <= 0 .AND. INTEM >= -5) THEN ! WITHOUT LIMITATION

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL                 ! 2. comp cell
         IT1     = II - IL                 ! First ghost cell
         IT0     = II - 2*IL               ! Second ghost cell
         IG      = I  + NN                ! Patch interior index
*
         IF((RKSI(II) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(II) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(II)  
            ROIM(IG)= RO(IT1)
         ELSE
c           DXL     = VOL(IT1)/(EPS+A(IT1)+A(I))
c           DXR     = VOL(I  )/(EPS+A(IT1)+A(IT2))
            DXL     = .5*D2(IT1)
            DXR     = .5*D2(I)

c           DX0     = A(IT1)/(VOL(IT0) + VOL(IT1))
c           DX1     = A(I)  /(VOL(IT1) + VOL(I  ))
c           DX2     = A(IT2)/(VOL(I  ) + VOL(IT2))
     
            DX0     = 1./(D2(IT0) + D2(IT1) + EPS)
            DX1     = 1./(D2(IT1) + D2(I)   + EPS)
            DX2     = 1./(D2(I)   + D2(IT2) + EPS)

            ROIP(IG)= RO(I)   - DXR*(DX1*(RO(I  )-RO(IT1))*RK1 +
     2                               DX2*(RO(IT2)-RO(I  ))*RK2)
            ROIM(IG)= RO(IT1) + DXL*(DX1*(RO(I) - RO(IT1))*RK1 +
     2                               DX0*(RO(IT1)-RO(IT0))*RK2)
       ENDIF

      ENDDO ; ENDDO

      ELSE IF(INTER <= 4) THEN ! FLUX LIMITATION

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL            ! 2. comp cell
         IT1     = II - IL            ! First ghost cell
         IT0     = II - 2*IL          ! Second ghost cell
         IG      = I  + NN               ! Patch interior index

         IF((RKSI(II) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
     +       RKSI(II) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
            ROIP(IG)= RO(II)  
            ROIM(IG)= RO(IT1)
         ELSE
c           DXL     = VOL(IT1)/(EPS+A(IT1)+A(I))
c           DXR     = VOL(I  )/(EPS+A(IT1)+A(IT2))
            DXL     = .5*D2(IT1)
            DXR     = .5*D2(I)

c           DX0     = A(IT1)/(VOL(IT0) + VOL(IT1))
c           DX1     = A(I)  /(VOL(IT1) + VOL(I  ))
c           DX2     = A(IT2)/(VOL(I  ) + VOL(IT2))

            DX0     = 1./(D2(IT0) + D2(IT1) + EPS)
            DX1     = 1./(D2(IT1) + D2(I)   + EPS)
            DX2     = 1./(D2(I)   + D2(IT2) + EPS)

            DIFF0   = DX0*(RO(IT1) - RO(IT0))
            DIFF1   = DX1*(RO(I)   - RO(IT1))
            DIFF2   = DX2*(RO(IT2) - RO(I))
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
c           RIII    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF0+EPS)))
c           PHIBI   = (RIII**2+RIII)/(RIII**2+1.)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
C           RII2    = MAX(-1.E+07,MIN( 1.E+07,DIFF1/(DIFF2+EPS)))
C           PHIFI   = (RII2**2+RII2)/(RII2**2+1.)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)

C ... EXTRAPOLATING
            ROIP(IG)= RO(I)  -DXR*PHIFI*(DIFF1*RK1+DIFF2*RK2)
            ROIM(IG)= RO(IT1)+DXL*PHIBI*(DIFF1*RK1+DIFF0*RK2)
         ENDIF

      ENDDO ; ENDDO

      ELSE ! FIRST ORDER UPWIND

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA + NFL
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1

      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IT2     = II + IL            ! 2. comp cell
         IT1     = II - IL            ! First ghost cell
         IT0     = II - 2*IL          ! Second ghost cell
         IG      = I  + NN               ! Patch interior index
         ROIP(IG)= RO(II)
         ROIM(IG)= RO(II-IL)
      ENDDO ; ENDDO

      ENDIF

      RETURN
      END SUBROUTINE XXTRAA_SURF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXVL(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,RKSI,INCHIML)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),RKSI(*)


      REAL, POINTER :: 
     1     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),
     2     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),
     3     F1P(:),F2P(:),F3P(:),F4P(:),F5P(:),
     4     F1M(:),F2M(:),F3M(:),F4M(:),F5M(:),
     5     F2V(:),F3V(:),F4V(:),F5V(:),PMI(:),PPI(:)
      REAL, TARGET ::  ZZZ(MAXW)
      
      LOGICAL INCHIML

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      F1P  => ZZZ(12*MAW+1:13*MAW);F2P  => ZZZ(13*MAW+1:14*MAW)
      F3P  => ZZZ(14*MAW+1:15*MAW);F4P  => ZZZ(15*MAW+1:16*MAW)
      F1M  => ZZZ(16*MAW+1:17*MAW);F2M  => ZZZ(17*MAW+1:18*MAW)
      F3M  => ZZZ(18*MAW+1:19*MAW);F4M  => ZZZ(19*MAW+1:20*MAW)
      F2V  => ZZZ(20*MAW+1:21*MAW);F3V  => ZZZ(21*MAW+1:22*MAW)
      F4V  => ZZZ(22*MAW+1:23*MAW);F5V  => ZZZ(23*MAW+1:24*MAW)
      F5M  => ZZZ(24*MAW+1:25*MAW);F5P  => ZZZ(25*MAW+1:26*MAW)


      YGM     = 1./GAMMA
      GM1     = GAMMA - 1.
      GM1PR   = 1./(GM1*PRLAM)
      G2M1    = 2.*(GAMMA**2 - 1.)
      G2M1    = 1./G2M1


C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW. UPDATED 16.6.1995
C
      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES
C
      CALL XXTRAP(RO,ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RM,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RN,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RW,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( E, EIP, EIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)

      DO  500 IG= 1,IJTRID
         I       = IA + IG

C ... PRESSURES AND LIMITATIONS. IDEAL GAS IS ASSUMED

      PPI(IG) = (GAMMA-1.)*(EIP(IG) -.5/ROIP(IG)*(RMIP(IG)**2
     2                             +RNIP(IG)**2+RWIP(IG)**2))
      PMI(IG) = (GAMMA-1.)*(EIM(IG) -.5/ROIM(IG)*(RMIM(IG)**2
     2                             +RNIM(IG)**2+RWIM(IG)**2))

         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
 500  CONTINUE


      DO  900 IG= 1,IJTRID
      I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

*         A21      = A2ZA(I) - A2YA(I)
*         A22      = A2XA(I) - A2ZA(I)
*         A23      = A2YA(I) - A2XA(I)
*         PALPO    = A21**2 + A22**2 + A23**2
*         IF(PALPO < 1.E-3) THEN
*            A21   = A2ZA(I) - A2YA(I)
*            A22   = A2XA(I) + A2ZA(I)
*            A23   =-A2YA(I) - A2XA(I)
*            PALPO = A21**2 + A22**2 + A23**2
*         ENDIF
*         PALPO    = 1./SQRT(PALPO)
*         A21      = A21*PALPO
*         A22      = A22*PALPO
*         A23      = A23*PALPO
*         A31      = A2YA(I)*A23 - A2ZA(I)*A22
*         A32      = A2ZA(I)*A21 - A2XA(I)*A23
*         A33      = A2XA(I)*A22 - A2YA(I)*A21

         RMIP(IG) = A11*RMIPI + A12*RNIPI + A13*RWIPI
         RNIP(IG) = A21*RMIPI + A22*RNIPI + A23*RWIPI
         RWIP(IG) = A31*RMIPI + A32*RNIPI + A33*RWIPI
         RMIM(IG) = A11*RMIMI + A12*RNIMI + A13*RWIMI
         RNIM(IG) = A21*RMIMI + A22*RNIMI + A23*RWIMI
         RWIM(IG) = A31*RMIMI + A32*RNIMI + A33*RWIMI

600   CONTINUE

C*******************************************************************
C ... VAN LEER'S FLUX-VECTOR SPLITTING                             *
C*******************************************************************

      YPRO    = 1./ROIP(IG)
      UIP     = YPRO*RMIP(IG)
      VIP     = YPRO*RNIP(IG)
      WIP     = YPRO*RWIP(IG)
      PIP     = MAX(0.01,PPI(IG))
      CIP     = SQRT(GAMMA*PIP*YPRO    )
      RMAIP   = UIP/CIP
      IF(RMAIP < -1.) THEN
              F1M(IG)  = A(I)*RMIP(IG)
              F2MIG    = A(I)*(RMIP(IG)*UIP+PIP)
              F3MIG    = A(I)*RMIP(IG)*VIP
              F4MIG    = A(I)*RMIP(IG)*WIP
              F5M(IG)  = A(I)*UIP*(EIP(IG)+PIP)
      ELSE IF(RMAIP < 1.) THEN
              F1MIG    = -A(I)*ROIP(IG)*CIP*(.5*(RMAIP-1.))**2
              F2MI     = GM1*UIP - 2.*CIP
              F2MIG    = F1MIG*F2MI*YGM
              F3MIG    = F1MIG*VIP
              F4MIG    = F1MIG*WIP
              F5M(IG)  = F1MIG*(EIP(IG) + PIP)*YPRO
              F1M(IG)  = F1MIG
      ELSE
              F1M(IG)  = 0.
              F2MIG    = 0.
              F3MIG    = 0.
              F4MIG    = 0.
              F5M(IG)  = 0.
      ENDIF
      F2M(IG) = F2MIG
      F3M(IG) = F3MIG
      F4M(IG) = F4MIG

601   CONTINUE

      YPRO    = 1./ROIM(IG)
      UIM     = YPRO*RMIM(IG)
      VIM     = YPRO*RNIM(IG)
      WIM     = YPRO*RWIM(IG)
      PIM     = MAX(0.01,PMI(IG))
      CIM     = SQRT(GAMMA*PIM*YPRO    )
      RMAIM   = UIM/CIM
      IF(RMAIM > 1.) THEN
              F1P(IG)  = A(I)*RMIM(IG)
              F2PIG    = A(I)*(RMIM(IG)*UIM+PIM)
              F3PIG    = A(I)*RMIM(IG)*VIM
              F4PIG    = A(I)*RMIM(IG)*WIM
              F5P(IG)  = A(I)*UIM*(EIM(IG)+PIM)
      ELSE IF(RMAIM > -1.) THEN
              F1PIG    = +A(I)*ROIM(IG)*CIM*(.5*(RMAIM+1.))**2
              F2PI     = GM1*UIM + 2.*CIM
              F2PIG    = F1PIG*F2PI*YGM
              F3PIG    = F1PIG*VIM
              F4PIG    = F1PIG*WIM
              F5P(IG)  = F1PIG*(EIM(IG) + PIM)*YPRO
              F1P(IG)  = F1PIG
      ELSE
              F1P(IG)  = 0.
              F2PIG    = 0.
              F3PIG    = 0.
              F4PIG    = 0.
              F5P(IG)  = 0.
      ENDIF
      F2P(IG) = F2PIG
      F3P(IG) = F3PIG
      F4P(IG) = F4PIG
800   CONTINUE

C
C ... FLUXES
C
      F2PIG   = F2P(IG) + F2M(IG)
      F3PIG   = F3P(IG) + F3M(IG)
      F4PIG   = F4P(IG) + F4M(IG)

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG

      FRO(I)  = F1M(IG) + F1P(IG) + FRO(I)
      FRM(I)  = F2PM              + FRM(I)
      FRN(I)  = F3PM              + FRN(I)
      FRW(I)  = F4PM              + FRW(I)
      FE(I)   = F5M(IG) + F5P(IG) +  FE(I)
900   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE FLUXVL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXWP(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,RKSI,INCHIML)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),RKSI(*)


      REAL, POINTER ::
     1     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),
     2     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),
     3      F1P(:), F2P(:), F3P(:),F4P(:), F5P(:),
     4      F1M(:), F2M(:), F3M(:),F4M(:), F5M(:),
     5      F2V(:), F3V(:), F4V(:),F5V(:), PMI(:),PPI(:)
      REAL, TARGET ::  ZZZ(MAXW)

      LOGICAL :: INCHIML

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      F1P  => ZZZ(12*MAW+1:13*MAW);F2P  => ZZZ(13*MAW+1:14*MAW)
      F3P  => ZZZ(14*MAW+1:15*MAW);F4P  => ZZZ(15*MAW+1:16*MAW)
      F1M  => ZZZ(16*MAW+1:17*MAW);F2M  => ZZZ(17*MAW+1:18*MAW)
      F3M  => ZZZ(18*MAW+1:19*MAW);F4M  => ZZZ(19*MAW+1:20*MAW)
      F2V  => ZZZ(20*MAW+1:21*MAW);F3V  => ZZZ(21*MAW+1:22*MAW)
      F4V  => ZZZ(22*MAW+1:23*MAW);F5V  => ZZZ(23*MAW+1:24*MAW)
      F5M  => ZZZ(24*MAW+1:25*MAW);F5P  => ZZZ(25*MAW+1:26*MAW)

      YGM     = 1./GAMMA
      GM1     = GAMMA - 1.
      GM1PR   = 1./(GM1*PRLAM)
      G2M1    = 2.*(GAMMA**2 - 1.)
      G2M1    = 1./G2M1

C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW. UPDATED 16.6.1995
C
      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES
C
      CALL XXTRAP(RO,ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RM,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RN,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(RW,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( E, EIP, EIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)

      DO  500 IG= 1,IJTRID
         I       = IA + IG

C ... PRESSURES AND LIMITATIONS. IDEAL GAS IS ASSUMED

      PPI(IG) = (GAMMA-1.)*(EIP(IG) -.5/ROIP(IG)*(RMIP(IG)**2
     2                             +RNIP(IG)**2+RWIP(IG)**2))
      PMI(IG) = (GAMMA-1.)*(EIM(IG) -.5/ROIM(IG)*(RMIM(IG)**2
     2                             +RNIM(IG)**2+RWIM(IG)**2))

         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
 500  CONTINUE


      DO  900 IG= 1,IJTRID
      I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

*         A21      = A2ZA(I) - A2YA(I)
*         A22      = A2XA(I) - A2ZA(I)
*         A23      = A2YA(I) - A2XA(I)
*         PALPO    = A21**2 + A22**2 + A23**2
*         IF(PALPO < 1.E-3) THEN
*            A21   = A2ZA(I) - A2YA(I)
*            A22   = A2XA(I) + A2ZA(I)
*            A23   =-A2YA(I) - A2XA(I)
*            PALPO = A21**2 + A22**2 + A23**2
*         ENDIF
*         PALPO    = 1./SQRT(PALPO)
*         A21      = A21*PALPO
*         A22      = A22*PALPO
*         A23      = A23*PALPO
*         A31      = A2YA(I)*A23 - A2ZA(I)*A22
*         A32      = A2ZA(I)*A21 - A2XA(I)*A23
*         A33      = A2XA(I)*A22 - A2YA(I)*A21

         RMIP(IG) = A11*RMIPI + A12*RNIPI + A13*RWIPI
         RNIP(IG) = A21*RMIPI + A22*RNIPI + A23*RWIPI
         RWIP(IG) = A31*RMIPI + A32*RNIPI + A33*RWIPI
         RMIM(IG) = A11*RMIMI + A12*RNIMI + A13*RWIMI
         RNIM(IG) = A21*RMIMI + A22*RNIMI + A23*RWIMI
         RWIM(IG) = A31*RMIMI + A32*RNIMI + A33*RWIMI

600   CONTINUE

C*******************************************************************
C ... WAVE/PARTICLE SPLITTING                                      *
C*******************************************************************

      YPRO    = 1./ROIP(IG)
      UIP     = YPRO*RMIP(IG)
      VIP     = YPRO*RNIP(IG)
      WIP     = YPRO*RWIP(IG)
      PIP     = MAX(0.01,PPI(IG))
      CIP     = SQRT(GAMMA*PIP*YPRO    )
      RMAIP   = UIP/CIP
      IF(RMAIP < -1.) THEN
              F1M(IG)  = A(I)*RMIP(IG)
              F2MIG    = A(I)*(RMIP(IG)*UIP+PIP)
              F3MIG    = A(I)*RMIP(IG)*VIP
              F4MIG    = A(I)*RMIP(IG)*WIP
              F5M(IG)  = A(I)*UIP*(EIP(IG)+PIP)
      ELSE IF(RMAIP < 1.) THEN
              F1MIG    = -A(I)*ROIP(IG)*CIP*(.5*(RMAIP-1.))**2
              RMA2     = .25*(RMAIP - 1.)**2*(2.+ RMAIP)
c              RMA2     = .25*(RMAIP + 1.)**2*(2.-RMAIP)
              F2MIG    = F1MIG*UIP          + A(I)*PIP*RMA2
              F3MIG    = F1MIG*VIP
              F4MIG    = F1MIG*WIP
              F5M(IG)  = F1MIG*EIP(IG)*YPRO + A(I)*UIP*PIP*RMA2
              F1M(IG)  = F1MIG
      ELSE
              F1M(IG)  = 0.
              F2MIG    = 0.
              F3MIG    = 0.
              F4MIG    = 0.
              F5M(IG)  = 0.
      ENDIF
      F2M(IG) = F2MIG
      F3M(IG) = F3MIG
      F4M(IG) = F4MIG

601   CONTINUE

      YPRO    = 1./ROIM(IG)
      UIM     = YPRO*RMIM(IG)
      VIM     = YPRO*RNIM(IG)
      WIM     = YPRO*RWIM(IG)
      PIM     = MAX(0.01,PMI(IG))
      CIM     = SQRT(GAMMA*PIM*YPRO    )
      RMAIM   = UIM/CIM
      IF(RMAIM > 1.) THEN
              F1P(IG)  = A(I)*RMIM(IG)
              F2PIG    = A(I)*(RMIM(IG)*UIM+PIM)
              F3PIG    = A(I)*RMIM(IG)*VIM
              F4PIG    = A(I)*RMIM(IG)*WIM
              F5P(IG)  = A(I)*UIM*(EIM(IG)+PIM)
      ELSE IF(RMAIM > -1.) THEN
              F1PIG    = +A(I)*ROIM(IG)*CIM*(.5*(RMAIM+1.))**2
              RMA2     = .25*(RMAIM + 1.)**2*(2.-RMAIM)
c              RMA2     = .25*(RMAIM - 1.)**2*(2.+ RMAIM)
              F2PIG    = F1PIG*UIM          + A(I)*PIM*RMA2
              F3PIG    = F1PIG*VIM
              F4PIG    = F1PIG*WIM
              F5P(IG)  = F1PIG*EIM(IG)*YPRO + A(I)*UIM*PIM*RMA2
              F1P(IG)  = F1PIG
      ELSE
              F1P(IG)  = 0.
              F2PIG    = 0.
              F3PIG    = 0.
              F4PIG    = 0.
              F5P(IG)  = 0.
      ENDIF
      F2P(IG) = F2PIG
      F3P(IG) = F3PIG
      F4P(IG) = F4PIG
800   CONTINUE

C
C ... FLUXES
C
      F2PIG   = F2P(IG) + F2M(IG)
      F3PIG   = F3P(IG) + F3M(IG)
      F4PIG   = F4P(IG) + F4M(IG)

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG

      FRO(I)  = F1M(IG) + F1P(IG) + FRO(I)
      FRM(I)  = F2PM              + FRM(I)
      FRN(I)  = F3PM              + FRN(I)
      FRW(I)  = F4PM              + FRW(I)
      FE(I)   = F5M(IG) + F5P(IG) +  FE(I)
900   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE FLUXWP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXAU(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,RKSI,INCHIML)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),RKSI(*)
      REAL :: FRSDEN, FRSPRE, CHAT2

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:), ERL(:),ELR(:),DEDPH(:),DEDRH(:),
     +     REIM(:),REIP(:),PDMI(:),PDPI(:)
      REAL, TARGET ::  ZZZ(MAXW)
      LOGICAL :: INCHIML

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      ERL  => ZZZ(12*MAW+1:13*MAW); ELR => ZZZ(13*MAW+1:14*MAW)
      PDMI => ZZZ(14*MAW+1:15*MAW);PDMI => ZZZ(15*MAW+1:16*MAW)
      RKIM => ZZZ(16*MAW+1:17*MAW);RKIP => ZZZ(17*MAW+1:18*MAW) 
      REIM => ZZZ(18*MAW+1:19*MAW);REIP => ZZZ(19*MAW+1:20*MAW)
      DEDRH=> ZZZ(20*MAW+1:21*MAW);DEDPH=> ZZZ(21*MAW+1:22*MAW)


C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW USING AUSM
C

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE PRIMITIVE VARIABLES
C
      CALL XXTRAP(RO ,ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
cc      CALL XXTRAP(   E, EIP, EIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   U,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   V,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   W,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   P, PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( PDI,PDPI,PDMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP(  RK, RKIP,RKIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( REPS,REIP,REIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      ENDIF

      DO  500 IG= 1,IJTRID
         I       = IA + IG
C ... CONTRAVARIANT MASS FLUXES (RMIP,UIP) AND PRESSURES
         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
         PDPI(IG)    = MAX(-0.999*FRSPRE,PDPI(IG))
         PDMI(IG)    = MAX(-0.999*FRSPRE,PDMI(IG))
 500  CONTINUE


C ... SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE AND DENSITY
C     ISTATE Determines the Equation of State to be  Used.
C     1 = Perfect Gas
C     2 = Chemically Reacting Equilibrium Gas (Air)

      CALL EFPRO(EIP,PPI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(EIM,PMI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ELR,PPI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ERL,PMI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)

C ... DERIVATIVES ACCORDING TO GLAISTER

      I1     = IA + 1

cc      CALL DERIV(RO(I1),P(I1),DEDPH,DEDRH,DRDP(I1),DRDH(I1),
cc     + EIP,EIM,ELR,ERL,PPI,PMI,ROIP,ROIM,IJTRID,IL)

      DO  600 IG= 1,IJTRID
         I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

*         A21      = A2ZA(I) - A2YA(I)
*         A22      = A2XA(I) - A2ZA(I)
*         A23      = A2YA(I) - A2XA(I)
*         PALPO    = A21**2 + A22**2 + A23**2
*         IF(PALPO < 1.E-3) THEN
*            A21   = A2ZA(I) - A2YA(I)
*            A22   = A2XA(I) + A2ZA(I)
*            A23   =-A2YA(I) - A2XA(I)
*            PALPO = A21**2 + A22**2 + A23**2
*         ENDIF
*         PALPO    = 1./SQRT(PALPO)
*         A21      = A21*PALPO
*         A22      = A22*PALPO
*         A23      = A23*PALPO
*         A31      = A2YA(I)*A23 - A2ZA(I)*A22
*         A32      = A2ZA(I)*A21 - A2XA(I)*A23
*         A33      = A2XA(I)*A22 - A2YA(I)*A21

         RMIP(IG) = A11*RMIPI + A12*RNIPI + A13*RWIPI
         RNIP(IG) = A21*RMIPI + A22*RNIPI + A23*RWIPI
         RWIP(IG) = A31*RMIPI + A32*RNIPI + A33*RWIPI
         RMIM(IG) = A11*RMIMI + A12*RNIMI + A13*RWIMI
         RNIM(IG) = A21*RMIMI + A22*RNIMI + A23*RWIMI
         RWIM(IG) = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... WAVE/PARTICLE SPLITTING USING LIOU'S METHOD (AUSM)           *
C ... IN THE CALCULATION OF THE SOUND SPEED DERIVATIVES            *
C ... ARE APPROXIMATED                                             *
C*******************************************************************
      UC      = UROT(I)

      YPROR   = 1./ROIP(IG)
      UR      = RMIP(IG)
      VR      = RNIP(IG)
      WR      = RWIP(IG)

      YPROL   = 1./ROIM(IG)
      UL      = RMIM(IG)
      VL      = RNIM(IG)
      WL      = RWIM(IG)

c ... ENTROPY FIX
c      UAVE    = .5*(UL+UR)
c      DISS    = 0.
c      DELTA   = .2*18.
c      IF(ABS(UAVE)  < DELTA)DISS=.5*(UAVE**2/DELTA + DELTA)

C ... TOTAL INTERNAL ENERGY

      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2)
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2)

      PR      = PPI(IG)
      PDR     = PDPI(IG) 
      HR      = EIR + PR*YPROR
      ER      = EIP(IG)
c      CR      = (PR/ROIP(IG)-ROIP(IG)*DEDRH(IG))/(ROIP(IG)*DEDPH(IG))
c      CR      = MAX(CR,0.00001*GAMMA*FRSPRE/FRSDEN)
c      CR      = SQRT(CR) 
c     CR      = sqrt(gamma*pr/roip(ig))
      CR      = C(I)
      RMAIP   = (UR-UC)/CR
      IF(ABS(RMAIP) <= 1.) THEN
         RMAR = -(.5*(RMAIP-1.))**2
         PDR  = PDR*.50*(RMAIP - 1.)**2*(2.+ RMAIP)
         PR   = PR*.50*(RMAIP - 1.)**2*(2.+ RMAIP)
      ELSE       
         RMAR = .5*(RMAIP-ABS(RMAIP))
         PDR  = PDR*(RMAIP-ABS(RMAIP))/RMAIP
         PR   = PR*(RMAIP-ABS(RMAIP))/RMAIP
      ENDIF      
                 
      PL      = PMI(IG)
      PDL     = PDMI(IG) 
      HL      = EIL + PL*YPROL
      EL      = EIM(IG)
c      CL      = (PL/ROIM(IG)-ROIM(IG)*DEDRH(IG))/(ROIM(IG)*DEDPH(IG))
c      CL      = MAX(CL,0.00001*GAMMA*FRSPRE/FRSDEN)
c      CL      = SQRT(CL) 
c      CL      = sqrt(gamma*pl/roim(ig))
      CL      = C(I-IL)
      RMAIM   = (UL-UC)/CL
      IF(ABS(RMAIM) <= 1.) THEN
         RMAL = (.5*(RMAIM+1.))**2
         PDL  = PDL*.50*(RMAIM + 1.)**2*(2.- RMAIM)
         PL   = PL*.50*(RMAIM + 1.)**2*(2.- RMAIM)
      ELSE       
         RMAL = .5*(RMAIM+ABS(RMAIM))
         PDL  = PDL*(RMAIM+ABS(RMAIM))/RMAIM
         PL   = PL*(RMAIM+ABS(RMAIM))/RMAIM
      ENDIF      
C *********************************************************
C     THIS CAN BE ACTIVATED FOR PERFECT GAS
C     GAMMA   = 1.4
C     GM1     = GAMMA - 1.
C     CHAT    = SQRT(GAMMA*PHAT/ROHAT)
C     PDPDE   = 1./(ROHAT*GM1)
C *********************************************************

      RMA     = RMAR + RMAL
      ABSRMA  = ABS(RMA) !+.5*DISS/(C(I)+C(I-IL))



C ... FLUX DIFFERENCES

      ROUCR   = ROIP(IG)*CR
      ROUCL   = ROIM(IG)*CL
      ROUUR   = ROUCR*UR
      ROUUL   = ROUCL*UL
      ROUVR   = ROUCR*VR
      ROUVL   = ROUCL*VL
      ROUWR   = ROUCR*WR
      ROUWL   = ROUCL*WL
      ROUER   = ROUCR*EIR
      ROUEL   = ROUCL*EIL

C ... ABSOLUTE VALUES OF THE CHARACTERISTIC SPEEDS


      F1PIG   = .5*A(I)*(RMA*(ROUCR + ROUCL) - ABSRMA*(ROUCR - ROUCL))
      F2PIG   = .5*A(I)*(RMA*(ROUUR + ROUUL) - ABSRMA*(ROUUR - ROUUL)+
     +           PDR + PDL)
      F3PIG   = .5*A(I)*(RMA*(ROUVR + ROUVL) - ABSRMA*(ROUVR - ROUVL))
      F4PIG   = .5*A(I)*(RMA*(ROUWR + ROUWL) - ABSRMA*(ROUWR - ROUWL))
      F5PIG   = .5*A(I)*(RMA*(ROUER + ROUEL) - ABSRMA*(ROUER - ROUEL)+
     +           PDPI(IG)*UR + PDMI(IG)*UL)
c     +           PDR*UR + PDL*UL)

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG +  FE(I)
600   CONTINUE   
                 
1000  CONTINUE   
                 
      RETURN     
      END SUBROUTINE FLUXAU        
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
                 
C	************ MODIFIED 7.6.2012 by Juho N *****************
      SUBROUTINE FLUINC(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,JSTATE,IPRESC,ARTSSP,XXTRAL,FRESUL,M,F1H,WAVE,
     7 FRSSIE,PRC,PSEUCO,RKSI,INCHIML,NGL,TRANSL,TRM,FRSTEM,IUPPTL)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IYCLE,INTEM,INTET,IDI,ISTATE,
     +     KMAXP1,KSTR,IDIR,ITURB,MAW,MAXW,ISTRID,JSTRID,KSTRID,NTOT,
     +     IJTRID,IL,KL,KG,IA,IG,I,ICYCLE,I8,IPRESC,M,INTEMM,NGL,
     +     IUPPTL,mmm,nnn,lll

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),D2(*),F1H(*),WAVE(*),RKSI(*)

      REAL :: FRSDEN,FRSPRE,PRLAM,PRT,RGAS,GAMMA,E0REF,T0REF,
     +   RKLIM,EPSLIM,
     +   RMIPI,RNIPI,RWIPI,RMIMI,RNIMI,RWIMI,A11,A12,A13,A21,A22,A23,
     +   A31,A32,A33,PALPO,UR,UL,VR,VL,WR,WL,YPROR,YPROL,RKR,RKL,RER,
     +   REL,EIR,EIL,PDR,PDL,PKR,PKL,ROHAT,CHAT,UC,RMIPP,RMIMM,UAVE1,
     +   UAVE,UDAMP,PAVE1,PAVE,F1PI,F1PIG,F2PIG,F3PIG,F4PIG,F5PIG,F6PIG,
     +   F7PIG,F2PM,F3PM,F4PM,ARTSSP,F8PIG,FRSSIE,APU,PSEUCO,F9PIG,
     +   F10PIG,F1PIF,UCONV,RMIPP1,RMIMM1,FRSTEM

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:), ERL(:),ELR(:),DEDPH(:),DEDRH(:),
     +     REIM(:),REIP(:),PDMI(:),PDPI(:),TEIM(:),TEIP(:),
     +     WHIM(:),WHIP(:),G1R(:),G1L(:),RET1R(:),RET1L(:)
      REAL, TARGET ::  ZZZ(MAXW)

      INTEGER :: JSTATE(*)

      LOGICAL :: XXTRAL, FRESUL, INCHIML, INCHIM, TRANSL

      TYPE(PRE_COR)       :: PRC(*)
      TYPE(INTERMITTENCY) :: TRM(*)

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      PDMI => ZZZ(14*MAW+1:15*MAW);PDPI => ZZZ(15*MAW+1:16*MAW)
      RKIM => ZZZ(16*MAW+1:17*MAW);RKIP => ZZZ(17*MAW+1:18*MAW) 
      REIM => ZZZ(18*MAW+1:19*MAW);REIP => ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(22*MAW+1:23*MAW);TEIP => ZZZ(23*MAW+1:24*MAW) 
      WHIM => ZZZ(24*MAW+1:25*MAW);WHIP => ZZZ(25*MAW+1:26*MAW) 
      G1R  => ZZZ(26*MAW+1:27*MAW);G1L  => ZZZ(27*MAW+1:28*MAW) 
      RET1R=> ZZZ(28*MAW+1:29*MAW);RET1L=> ZZZ(29*MAW+1:30*MAW) 


C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW. SIMPLIFIED SPLITTING
C ... THIS SUBROUTINE UTILIZES SPECIFIC INTERNAL ENERGY
C
      KSTRID  = KMAX + 2*KN	!KN = ghost cells = 2 per side ??
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID

*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT (PRIMITIVE) VARIABLES

      IF(.NOT. XXTRAL) THEN

      INCHIM = INCHIML ! Mersu

      CALL XXTRAP(TEMP,TEIP,TEIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)
      CALL XXTRAP(   U,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)
      CALL XXTRAP(   V,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)
      CALL XXTRAP(   W,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)
      CALL XXTRAP(   P, PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)
      CALL XXTRAP( PDI,PDPI,PDMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIM)

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP(  RK, RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIM)
      CALL XXTRAP( REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIM)
      ENDIF

	IF(FRESUL .AND. M == 1) THEN
      CALL XXTRAP(WAVE,WHIP,WHIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      ENDIF

	IF(TRANSL) THEN
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
      ENDIF

      ELSE	!IF(.NOT. XXTRAL) THEN

      I8 =IJTRID
      CALL XXTRAA(TEMP,TEIP,TEIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   U,RMIP,RMIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   V,RNIP,RNIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   W,RWIP,RWIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   P, PPI, PMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA( PDI,PDPI,PDMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAA( RK, RKIP,RKIM,I8,IL,KL,IA,INTET,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(REPS,REIP,REIM,I8,IL,KL,IA,INTET,A,VOL,D2,RKSI,INCHIM)
      ENDIF

      IF(FRESUL .AND. M == 1)  THEN
      CALL XXTRAA(WAVE,WHIP,WHIM,I8,IL,KL,IA,INTET,A,VOL,D2,RKSI,INCHIM)
      ENDIF

      IF(TRANSL) THEN ! This is not distance-based. Who will write TXTRAA?
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
      ENDIF

      ENDIF ! .NOT.XXTRAL

      DO  500 IG  = 1,IJTRID
         I        = IA + IG
         IF(ISTATE /= 10) THEN ! Limit the pressure
	   ! P and PD = pressure and pres diff 
         PPI(IG)  = MAX(0.001*FRSPRE,PPI(IG)) 
         PMI(IG)  = MAX(0.001*FRSPRE,PMI(IG))
         PDPI(IG) = MAX(-0.999*FRSPRE,PDPI(IG))
         PDMI(IG) = MAX(-0.999*FRSPRE,PDMI(IG))
         ENDIF
         RKIP(IG) = MAX(0.,RKIP(IG))
         RKIM(IG) = MAX(0.,RKIM(IG))
         REIP(IG) = MAX(0.,REIP(IG)) !RE = epsilon
         REIM(IG) = MAX(0.,REIM(IG))
 500  CONTINUE

C ... DENSITY AND SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE
C     AND TEMPERATURE. ISTATE Determines the Equation of State.

      CALL ROFPT(TEIP,PPI,ROIP,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL ROFPT(TEIM,PMI,ROIM,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL EFPT(EIP,PPI,TEIP,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)
      CALL EFPT(EIM,PMI,TEIM,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)
                        
      DO  600 IG  = 1,IJTRID  !Momentum = U*RO
         I        = IA + IG

         RMIPI    = RMIP(IG)	!X plus = X right
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)	!Y minus = Y left
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)	! unit vectors of surface 2?? of the cells
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         UR       = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VR       = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WR       = A31*RMIPI + A32*RNIPI + A33*RWIPI
         UL       = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VL       = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WL       = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... SIMPLIFIED ROE BASED FLUX SPLITTING                          *
C*******************************************************************

c      YPROR   = 1./ROIP(IG)
c      YPROL   = 1./ROIM(IG)
      RKR     = RKIP(IG) !*YPROR
      RER     = REIP(IG) !*YPROR
      RKL     = RKIM(IG) !*YPROL
      REL     = REIM(IG) !*YPROL

C ... TOTAL INTERNAL ENERGY

c      PDR     = PDI(I) !  PDPI(IG)
      PDR     = PDPI(IG)		    !PD = pressure difference ?
      PKR     = .6667*ROIP(IG)*RKIP(IG) !PKR = total internal energy ?

c      PDL     = PDI(I-IL) ! PDMI(IG)
      PDL     = PDMI(IG)
      PKL     = .6667*ROIM(IG)*RKIM(IG)
      ROHAT   =.5*(ROIP(IG)+ROIM(IG))

      IF(IPRESC == 1) THEN
          CHAT = ARTSSP
      ELSE IF(IPRESC == 2) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE IF(IPRESC == 0) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE
          WRITE(*,*) 'Wrong IPRESC option. Exiting..'
          STOP
      ENDIF
      UC      = UROT(I)	!rotational velocity or static moment of the cell face ??

* Equations replaced by ESa 3.9.2001 (equations from 'propeller-finflo' 

c      RMIPP   = A11*RM(I)    + A12*RN(I)    + A13*RW(I)
c      RMIMM   = A11*RM(I-IL) + A12*RN(I-IL) + A13*RW(I-IL)
c      UAVE1   = (RMIPP/RO(I) + RMIMM/RO(I-IL))*.5 - UC
        call ijkpai(I,imax,jmax,kmax,mmm,nnn,lll)

      RMIPP   = A11*U(I)    + A12*V(I)    + A13*W(I)
      RMIMM   = A11*U(I-IL) + A12*V(I-IL) + A13*W(I-IL)

c      RMIPP = UR ; RMIMM = UL   ! Testaa joskus    

      UAVE1   = (RMIPP + RMIMM)*.5 - UC
c      APU     = 2./(1./PRC(I)%AP + 1./PRC(I-IL)%AP)
c      APU     = APU/(1.e-10+ A(I)) ! Contain AP of Rhie and Chow (UAVE1?) . A(I)*ROHAT?

      UDAMP   = MAX(PSEUCO*ABS(UAVE1),ARTSSP)     ! ,APU) !CHAT  ! Alternatives   
C laivakokeilu
c      UDAMP = 14. ! Lirutetaan alespäin
c        udamp = MIN(UDAMP,CHAT) ! For compressible flows ?
      UAVE    = UAVE1 - (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP) ! K added 14.10.08
c      IF(M >= 2) UAVE = UAVE1 ! Ei vaikuta
      PAVE1   = .5*(PDL+PKL + PDR+PKR)
c      PAVE1   = .5*(PDI(I)+PKL + PDI(I-IL)+PKR)
      PAVE    = PAVE1 !-.5*ROHAT*UDAMP*(UR-UL) ! Dual dissipation scheme


      IF(UAVE >= 0.) THEN     ! if U_average positive -> use minus/left values
c      F1PI    = UAVE*ROHAT !ROIM(IG)
      F1PI    = UAVE*ROIM(IG)
      F1PIG   = A(I)*F1PI
      F2PIG   = F1PIG*UL + A(I)*PAVE

      F3PIG   = F1PIG*VL
      F4PIG   = F1PIG*WL
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL
      F5PIG   = F1PIG*EIL +.5*A(I)*(UAVE+UC)*(PDR+PKR+PDL+PKL+2.*FRSPRE) 
c      F5PIG   = F1PIG*EIL +.5*A(I)*(UAVE+UC)*(PKR+PKL + PPI(IG)+PMI(IG)) ! Hyvät hyssykät
      F6PIG   = F1PIG*RKL
      F7PIG   = F1PIG*REL
      IF(FRESUL .AND. M == 1)  F8PIG = F1PIG*WHIM(IG)

C	For this index I use left(L) values	
      IF(TRANSL) THEN
         F9PIG  = F1PIG * G1L(IG)		
         F10PIG = F1PIG * RET1L(IG)		  
      ENDIF



      ELSE	! if U_average negative -> use plus/right values

c      F1PI    = UAVE*ROHAT! ROIP(IG)
      F1PI    = UAVE*ROIP(IG)
      F1PIG   = A(I)*F1PI
      F2PIG   = F1PIG*UR + A(I)*PAVE
      F3PIG   = F1PIG*VR
      F4PIG   = F1PIG*WR
      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
      F5PIG   = F1PIG*EIR +.5*A(I)*(UAVE+UC)*(PDR+PKR+PDL+PKL+2*FRSPRE) 
      F6PIG   = F1PIG*RKR
      F7PIG   = F1PIG*RER

      IF(FRESUL .AND. M == 1)  F8PIG = F1PIG*WHIP(IG)


C	For this index I use right(R) values
      IF(TRANSL) THEN
         F9PIG  = F1PIG * G1R(IG)		
         F10PIG = F1PIG * RET1R(IG)		    
      ENDIF


      ENDIF


c      if(ngl == 1 .and. idir == 2) then
c        write(6667,*) mmm,nnn,f1pig,f2pig,f3pig,f4pig,f5pig
c      endif

C*    HATS SHOULD BE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.94

      HAT1(I) = 0.  
      HAT2(I) = F1PI
      HAT3(I) = HAT2(I)
      HAT4(I) = SIGN(1.,F1PIG)*HAT2(I)/ROHAT

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
c      if(ngl == 1.and.idir == 2.and.kg == 1)write(6667,*) i,fro(i),f1pig
c      if(ngl == 1.and.i <= 808.and.idir == 2) write(6701,*)i,ig,
c     +roim(ig),uave,
c     + roip(ig),pdl,pkl,rmim(ig),rmip(ig),ul,ur

      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG +  FE(I)
      FRK(I)  = F6PIG + FRK(I)
      FEPS(I) = F7PIG + FEPS(I)
      IF(FRESUL .AND. M == 1)  F1H(I) = F8PIG + F1H(I)

C     Update transition model fluxes with the massflow part
      IF(TRANSL) THEN
         TRM(I)%FG   = TRM(I)%FG   + F9PIG
         TRM(I)%FRET = TRM(I)%FRET + F10PIG
      ENDIF

600   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE FLUINC
C	************** MOD ENDS 7.6.2012 by Juho N ***************
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE JAMESON(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,JSTATE,IPRESC,ARTSSP,XXTRAL,FRESUL,M,F1H,WAVE,
     7 FRSSIE,PRC,PSEUCO,RJK2,RJK4,SCALEJ)

      USE NS3CO, ONLY : IN, JN, KN
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IYCLE,INTEM,INTET,IDI,ISTATE,
     +     KMAXP1,KSTR,IDIR,ITURB,MAW,MAXW,ISTRID,JSTRID,KSTRID,NTOT,
     +     IJTRID,IL,KL,KG,IA,IG,I,ICYCLE,I8,IPRESC,M,INTEMM

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),D2(*),F1H(*),WAVE(*),SCALEJ(*)

      REAL :: FRSDEN,FRSPRE,PRLAM,PRT,RGAS,GAMMA,E0REF,T0REF,
     +   RKLIM,EPSLIM,
     +   RMIPI,RNIPI,RWIPI,RMIMI,RNIMI,RWIMI,A11,A12,A13,A21,A22,A23,
     +   A31,A32,A33,PALPO,UR,UL,VR,VL,WR,WL,YPROR,YPROL,RKR,RKL,RER,
     +   REL,EIR,EIL,PDR,PDL,PKR,PKL,ROHAT,CHAT,UC,RMIPP,RMIMM,UAVE1,
     +   UAVE,UDAMP,PAVE1,PAVE,F1PI,F1PIG,F2PIG,F3PIG,F4PIG,F5PIG,F6PIG,
     +   F7PIG,F2PM,F3PM,F4PM,ARTSSP,F8PIG,FRSSIE,APU,PSEUCO,RJK2,RJK4,
     +   EPS2I,EPS4I,RLAMI

      REAL, POINTER ::   
     +     PRES1(:),PRES2(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:), ERL(:),ELR(:),DEDPH(:),DEDRH(:),
     +     REIM(:),REIP(:),PDMI(:),PDPI(:),TEIM(:),TEIP(:),
     +     WHIM(:),WHIP(:)
      REAL, TARGET ::  ZZZ(MAXW)

      INTEGER :: JSTATE(*)

      LOGICAL :: XXTRAL, FRESUL

      TYPE(PRE_COR) :: PRC(*)

      PRES1 => ZZZ( 0*MAW+1: 1*MAW);PRES2 => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      PDMI => ZZZ(14*MAW+1:15*MAW);PDPI => ZZZ(15*MAW+1:16*MAW)
      RKIM => ZZZ(16*MAW+1:17*MAW);RKIP => ZZZ(17*MAW+1:18*MAW) 
      REIM => ZZZ(18*MAW+1:19*MAW);REIP => ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(22*MAW+1:23*MAW);TEIP => ZZZ(23*MAW+1:24*MAW) 
      WHIM => ZZZ(24*MAW+1:25*MAW);WHIP => ZZZ(25*MAW+1:26*MAW) 


CSCALEJ(I),SCALEJ(I-IL)
C ... FLUXES FOR THREE-DIMENSIONAL FLOW. SIMPLIFIED SPLITTING
C ... THIS SUBROUTINE UTILIZES SPECIFIC INTERNAL ENERGY
C
      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID

*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

C ... First pressure sensor

      KG      = 0
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID

      DO  IG  = 1,IJTRID
c         PRES1(IG) = 0.
         I        = IA + IG
         PRES1(IG) = ABS(PDI(I+IL)-2.*PDI(I)+PDI(I-IL)) /
     +                  (PDI(I+IL)+2.*PDI(I)+PDI(I-IL) + 4.*FRSPRE)
      ENDDO
       
      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... Pressure sensor


      DO  IG  = 1,IJTRID
         I        = IA + IG
         PRES2(IG) = ABS(PDI(I+IL)-2.*PDI(I)+PDI(I-IL)) /
     +                  (PDI(I+IL)+2.*PDI(I)+PDI(I-IL) + 4.*FRSPRE)

      ENDDO

C ... Calculate Jameson type dissipation

      DO  600 IG  = 1,IJTRID
         I        = IA + IG
         EPS2I    = RJK2*MAX(PRES1(IG),PRES2(IG))
cc        write(8004,*) i,ig,eps2I
cc         eps2i=0.
         EPS4I    = MAX(0.,(RJK4-EPS2I))
c          eps4i = 1.E-2 ; eps2i = 0. ! trials
C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES
cc        write(8001,*) ig,pres1(ig),pres2(ig)

c         UC      = UROT(I)
c         RMIPP   = A11*U(I)    + A12*V(I)    + A13*W(I)
c         RMIMM   = A11*U(I-IL) + A12*V(I-IL) + A13*W(I-IL)
c         UAVE1   = (RMIPP + RMIMM)*.5 - UC

      IF(IPRESC == 1) THEN
          CHAT = ARTSSP
      ELSE IF(IPRESC == 2) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE IF(IPRESC == 0) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE
          WRITE(*,*) 'Wrong IPRESC option. Exiting..'
          STOP
      ENDIF

c      RLAMI   = (ABS(UAVE1) + CHAT)
c      write(7999,*) IG,I,RLAMI,EPS2I,EPS4I,SCALEJ(I),SCALEJ(I-IL)
c      RLAMI   = 1000.*RLAMI

      RLAMI = .5*(SCALEJ(I)+SCALEJ(I-IL))
c         rlami = a(i)*15. ! A trial 
      F1PIG   = RLAMI*(EPS2I*(RO(I) - RO(I-IL))   -
     +  EPS4I*(RO(I+IL) - 3.*RO(I) + 3.*RO(I-IL) - RO(I-2*IL)))

      IF(IPRESC /= 0) THEN
      F1PIG   = RLAMI*(EPS2I*(PDI(I) - PDI(I-IL))   -
     +  EPS4I*(PDI(I+IL) - 3.*PDI(I) + 3.*PDI(I-IL) - PDI(I-2*IL))) /
     +  CHAT**2
      ENDIF

      F2PIG   = RLAMI*(EPS2I*(RM(I) - RM(I-IL))   -
     +  EPS4I*(RM(I+IL) - 3.*RM(I) + 3.*RM(I-IL) - RM(I-2*IL)))
      F3PIG   = RLAMI*(EPS2I*(RN(I) - RN(I-IL))   -
     +  EPS4I*(RN(I+IL) - 3.*RN(I) + 3.*RN(I-IL) - RN(I-2*IL)))
      F4PIG   = RLAMI*(EPS2I*(RW(I) - RW(I-IL))   -
     +  EPS4I*(RW(I+IL) - 3.*RW(I) + 3.*RW(I-IL) - RW(I-2*IL)))
c      F5PIG   = RLAMI*(EPS2I*(RO(I)*E(I) - RO(I-IL)*E(I-IL))  -
c     +  EPS4I*(RO(I+IL)*E(I+IL) - 3.*RO(I)*E(I) + 3.*RO(I-IL)*E(I-IL)
c     +       - RO(I-2*IL)*E(I-2*IL)))
      F5PIG   = RLAMI*(EPS2I*(E(I) - E(I-IL))  -
     +  EPS4I*(E(I+IL) - 3.*E(I) + 3.*E(I-IL)
     +       - E(I-2*IL)))

c      F6PIG   = RLAMI*(EPS2I*(RK(I) - RK(I-IL))   -
c     +  EPS4I*(RK(I+IL) - 3.*RK(I) + 3.*RK(I-IL) - RK(I-2*IL)))
c      F7PIG   = RLAMI*(EPS2I*(REPS(I) - REPS(I-IL))   -
c     +  EPS4I*(REPS(I+IL) - 3.*REPS(I) + 3.*REPS(I-IL) - REPS(I-2*IL)))

C
C ... FLUXES
C
c      write(8000,*) IG,I,F2pig,FRM(I),
c     + RM(I+IL),RM(I),RM(I-IL),RM(I-2*IL)
      FRO(I)  = -F1PIG  + FRO(I)
      FRM(I)  = -F2PIG  + FRM(I)
      FRN(I)  = -F3PIG  + FRN(I)
      FRW(I)  = -F4PIG  + FRW(I)
      FE(I)   = -F5PIG  +  FE(I) ! Does not work for water (Comment)
C ... For turbulence a limited second-order scheme is currently applied.
c      FRK(I)  = -F6PIG + FRK(I)
c      FEPS(I) = -F7PIG + FEPS(I)

      PRES1(IG)= PRES2(IG)

      IF(FRESUL .AND. M == 1)  THEN
         WRITE(*,*)'DO NOT USE Jameson dissipation with a free surface.'
         WRITE(13,*)'DO NOT USE Jameson dissipation with a free surface'
         STOP 'Jameson'
      ENDIF
600   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE JAMESON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUMPH(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,UG,VG,WG,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDI,PRO,VAR,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,JSTATE,IPRESC,ARTSSP,PSEUCO,RKSI,INCHIML,NGL,
     7 XC,YC,ZC,GRAVIL,MULPHC,FRSTEM,IUPPTL,INTEA,PRC,FRADEN)

      USE NS3CO, ONLY : IN, JN, KN, GX, GY, GZ, GROUND
      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : EPS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IYCLE,INTEM,INTET,IDI,ISTATE,
     +     KMAXP1,KSTR,IDIR,ITURB,MAW,MAXW,ISTRID,JSTRID,KSTRID,NTOT,
     +     IJTRID,IL,KL,KG,IA,IG,I,ICYCLE,IPHASE,II,IPRESC,NGL,IUPPTL,
     +     MMM,NNN,LLL,iii,INTEA

      REAL :: A(*),RO(*),RN(*),RW(*),UG(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),VG(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),WG(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDI(*),D2(*),RKSI(*),FRADEN(*)

      REAL :: FRSDEN,FRSPRE,PRLAM,PRT,RGAS,GAMMA,E0REF,T0REF,
     +   RKLIM,EPSLIM,
     +   RMIPI,RNIPI,RWIPI,RMIMI,RNIMI,RWIMI,A11,A12,A13,A21,A22,A23,
     +   A31,A32,A33,PALPO,UR,UL,VR,VL,WR,WL,YPROR,YPROL,RKR,RKL,RER,
     +   REL,EIR,EIL,PDR,PDL,PKR,PKL,ROHAT,CHAT,UC,RMIPP,RMIMM,UAVE1,
     +   UAVE,UDAMP,PAVE1,PAVE,F1PI,F1PIG,F2PIG,F3PIG,F4PIG,F5PIG,F6PIG,
     +   F7PIG,F2PM,F3PM,F4PM,AAVE,F1PIA,ROHATC,ARTSSP,UMAX,UMIN,UAV2,
     +   PSEUCO,XRG,YRG,ZRG,FX,FY,FZ,FRSTEM,XAVE,YMXA,F1PIX,AAVE1,
     +   VISJ,DISTJ,FREDIF,UDAMC

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:), T1L(:),T1R(:),A1L(:),  A1R(:),
     +     REIM(:),REIP(:),PDMI(:),PDPI(:),TEIM(:),TEIP(:),
     +      T2L(:), T2R(:), F1R(:),F1RM(:),F1RN(:),F1RW(:),
     +      F1E(:),F1RK(:),F1EPS(:),UAVF1(:),PAVF1(:),AAVF1(:),
     +     ULF1(:),URF1(:),VLF1(:),VRF1(:),WlF1(:),WRF1(:),
     +    HSATR(:),HSATL(:),CPR(:),CPL(:),TSATR(:),TSATL(:)

      REAL, TARGET ::  ZZZ(MAXW)

      REAL :: XC(*),YC(*),ZC(*)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(PRE_COR) :: PRC(*)

      INTEGER :: JSTATE(*)

      LOGICAL :: INCHIML, GRAVIL

      CHARACTER(LEN = *) :: MULPHC

      REAL, ALLOCATABLE, DIMENSION(:) :: U,V,W

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      A1L  => ZZZ(12*MAW+1:13*MAW); A1R => ZZZ(13*MAW+1:14*MAW)
      PDMI => ZZZ(14*MAW+1:15*MAW);PDPI => ZZZ(15*MAW+1:16*MAW)
      RKIM => ZZZ(16*MAW+1:17*MAW);RKIP => ZZZ(17*MAW+1:18*MAW) 
      REIM => ZZZ(18*MAW+1:19*MAW);REIP => ZZZ(19*MAW+1:20*MAW)
      T1L  => ZZZ(20*MAW+1:21*MAW);T1R  => ZZZ(21*MAW+1:22*MAW)
      TEIM => ZZZ(22*MAW+1:23*MAW);TEIP => ZZZ(23*MAW+1:24*MAW) 
      F1R  => ZZZ(24*MAW+1:25*MAW);F1RM => ZZZ(25*MAW+1:26*MAW) 
      F1RN => ZZZ(26*MAW+1:27*MAW);F1RW => ZZZ(27*MAW+1:28*MAW) 
      F1E  => ZZZ(28*MAW+1:29*MAW);F1RK => ZZZ(29*MAW+1:30*MAW) 
      F1EPS=> ZZZ(30*MAW+1:31*MAW);UAVF1=> ZZZ(31*MAW+1:32*MAW)
      PAVF1=> ZZZ(32*MAW+1:33*MAW);AAVF1=> ZZZ(33*MAW+1:34*MAW)
      URF1 => ZZZ(34*MAW+1:35*MAW);ULF1 => ZZZ(35*MAW+1:36*MAW)
      VRF1 => ZZZ(36*MAW+1:37*MAW);VLF1 => ZZZ(37*MAW+1:38*MAW)
      WRF1 => ZZZ(38*MAW+1:39*MAW);WLF1 => ZZZ(39*MAW+1:40*MAW)
      HSATR=> ZZZ(40*MAW+1:41*MAW);HSATL=> ZZZ(41*MAW+1:42*MAW)
      CPR  => ZZZ(42*MAW+1:43*MAW);CPL  => ZZZ(43*MAW+1:44*MAW)
      TSATR=> ZZZ(44*MAW+1:45*MAW);TSATL=> ZZZ(45*MAW+1:46*MAW)

C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW. SIMPLIFIED SPLITTING
C ... THIS SUBROUTINE UTILIZES SPECIFIC INTERNAL ENERGY
C

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

	ALLOCATE(U(NTOT),V(NTOT),W(NTOT))	

      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF
        
      DO 2000 KG = 1,KMAXP1
*      IA         = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA         = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT (PRIMITIVE) VARIABLES

      CALL XXTRAP(TEMP,TEIP,TEIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   P, PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( PDI,PDPI,PDMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP(  RK,RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP(REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      ENDIF
         
      DO  1500 IG = 1,IJTRID
         I        = IA + IG
         IF(ISTATE /= 10) THEN ! Limit the pressure
         PPI(IG)  = MAX(0.001*FRSPRE,PPI(IG))
         PMI(IG)  = MAX(0.001*FRSPRE,PMI(IG))
         PDPI(IG) = MAX(-0.999*FRSPRE,PDPI(IG))
         PDMI(IG) = MAX(-0.999*FRSPRE,PDMI(IG))
         ENDIF
         RKIP(IG) = MAX(0.,RKIP(IG))
         RKIM(IG) = MAX(0.,RKIM(IG))
         REIP(IG) = MAX(0.,REIP(IG))
         REIM(IG) = MAX(0.,REIM(IG))
1500  CONTINUE
                        

      DO 1900 IPHASE = 1,NPHASES

       IF(MULPHC == 'CAVIT' .AND. IPHASE == 1) THEN ! Homogeneus model

        U(1:NTOT) = UG(1:NTOT)
        V(1:NTOT) = VG(1:NTOT)
        W(1:NTOT) = WG(1:NTOT)

        CALL XXTRAP(   U,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
        CALL XXTRAP(   V,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
        CALL XXTRAP(   W,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)

       ELSE IF(MULPHC == 'MULTI') THEN ! Two-fluid model

        U(1:NTOT) = VAR(1:NTOT)%U(IPHASE)
        V(1:NTOT) = VAR(1:NTOT)%V(IPHASE)
        W(1:NTOT) = VAR(1:NTOT)%W(IPHASE)

        CALL XXTRAP(   U,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
        CALL XXTRAP(   V,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
        CALL XXTRAP(   W,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
       ENDIF

      CALL DXTRAP(PRO,VAR,T1R,T1L,A1R,A1L,HSATR,HSATL,CPR,CPL,TSATR,
     + TSATL,IJTRID,IL,KL,IA,1,IPHASE,RKSI,INCHIML,MULPHC,INTEA) ! Limiter is default

C ... DENSITY AND SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE
C     AND TEMPERATURE. ISTATE Determines the Equation of State.

      ISTATE = JSTATE(IPHASE)

      CALL ROFPT(T1R,PPI,ROIP,IJTRID,ISTATE,RGAS,GAMMA,FRADEN(IPHASE),
     +           FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
      CALL ROFPT(T1L,PMI,ROIM,IJTRID,ISTATE,RGAS,GAMMA,FRADEN(IPHASE),
     +           FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
cc      CALL EFPT (EIP,PPI,TEIP,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,E0REF,T0REF)
cc      CALL EFPT (EIM,PMI,TEIM,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,E0REF,T0REF)
      CALL EFPTG(EIP,PPI,T1R,HSATR,CPR,TSATR,ROIP,PRO,IJTRID,ISTATE,
     +  RGAS,GAMMA,0.,0.,E0REF,T0REF,IPHASE)
      CALL EFPTG(EIM,PMI,T1L,HSATL,CPL,TSATL,ROIM,PRO,IJTRID,ISTATE,
     + RGAS,GAMMA,0.,0.,E0REF,T0REF,IPHASE)

      DO  1800 IG = 1,IJTRID
         I        = IA + IG
c        call ijkpai(I,imax,jmax,kmax,mmm,nnn,lll)

      IF(MULPHC == 'MULTI' .OR. MULPHC == 'CAVIT')THEN

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         URF1(IG) = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VRF1(IG) = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WRF1(IG) = A31*RMIPI + A32*RNIPI + A33*RWIPI
         ULF1(IG) = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VLF1(IG) = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WLF1(IG) = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... Average face velocities for unequa velocities                *
C*******************************************************************

         UC        = UROT(I)
         RMIPP     = A11*U(I)    + A12*V(I)    + A13*W(I)
         RMIMM     = A11*U(I-IL) + A12*V(I-IL) + A13*W(I-IL)
         UAVE1     = (RMIPP + RMIMM)*.5 - UC
         UAVF1(IG) = UAVE1

      ENDIF ! MULPHC == 'MULTI' or 'CAVIT'

         A1R(IG)  = MIN(1.,A1R(IG))
         A1L(IG)  = MIN(1.,A1L(IG))
         A1R(IG)  = MAX(0.,A1R(IG))
         A1L(IG)  = MAX(0.,A1L(IG))
         IF(A1R(IG) < 2.5E-7) A1R(IG) = 0. ! Intel epsilon 6.E-8
         IF(A1L(IG) < 2.5E-7) A1L(IG) = 0.
         
         F1PIG    = 0.
         F2PIG    = 0.
         F3PIG    = 0.
         F4PIG    = 0.
         F5PIG    = 0.
         F6PIG    = 0.
         F7PIG    = 0.

C*******************************************************************
C ... SIMPLE UPWIND FLUX                        *
C*******************************************************************

      RKR     = RKIP(IG)
      RER     = REIP(IG)
      RKL     = RKIM(IG)
      REL     = REIM(IG)

      UR      = URF1(IG)
      UL      = ULF1(IG)
      VR      = VRF1(IG)
      VL      = VLF1(IG)
      WR      = WRF1(IG)
      WL      = WLF1(IG)

      PDR     = PDPI(IG)
      PKR     = .6667*RKIP(IG)*RO(I) ! ROIP(IG)

      PDL     = PDMI(IG)
      PKL     = .6667*RKIM(IG)*RO(I-IL) !ROIM(IG)

c      ROHAT   = SQRT(ROIP(IG)*ROIM(IG))
      ROHAT    =.5*(ROIP(IG)+ROIM(IG)) ! Density of IPHASE 
c      ROHATC  = 1./(.5*(1./RO(I) + 1./RO(I-IL)))
      ROHATC   = .5*(RO(I) + RO(I-IL))             ! Total density and
      AAVE    = .5*(VAR(I)%ALFA(IPHASE) + VAR(I-IL)%ALFA(IPHASE))
      XAVE     = .5*(VAR(I)%X(IPHASE) + VAR(I-IL)%X(IPHASE)) ! fraction
      YMXA     = 1. - XAVE
      YMXA     = MIN(1.,YMXA)
      YMXA     = MAX(0.,YMXA)

      IF(IPRESC == 1) THEN
          CHAT = ARTSSP
      ELSE IF(IPRESC == 2) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE IF(IPRESC == 0) THEN
          CHAT =  MAX(C(I),C(I-IL))
      ELSE
          WRITE(*,*) 'Wrong IPRESC option. Exiting..'
          STOP
      ENDIF

      UC      = UROT(I)
      UAVE1   = UAVF1(IG)
      UDAMC   = MAX(PSEUCO*ABS(UAVF1(IG)),ARTSSP)   !CHAT  ! 
      PAVE    = .5*(PDL + PDR)
c      PAVE    = .5*(PDI(I) + PDI(I-IL))
      PAVE1   = .5*(PKL + PKR)
c      PAVE    = PAVE1 !-.5*ROHAT*CHAT*(UR - UL)! Dual dissipation scheme
      AAVE1   = 2.*(VAR(I)%ALFA(IPHASE) * VAR(I-IL)%ALFA(IPHASE))/
     +             (VAR(I)%ALFA(IPHASE) + VAR(I-IL)%ALFA(IPHASE)+EPS)

      UDAMP   = (PDR+PKR - PDL-PKL)/(2.*ROHATC*UDAMC) ! K added 14.10.08
c      UDAMP   = (PDR - PDL)/(2.*ROHATC*UDAMP)
      UAVE    = UAVE1 - UDAMP 
      UAV2    = UAVE + UC
       
c     HAT4(I) = 1./ROHAT ! Scalars do not work
    
C ... Gravitational source (not activated)

      XRG     = GX*(.5*(XC(I)+XC(I-IL))-GROUND)
      YRG     = GY*(.5*(YC(I)+YC(I-IL))-GROUND)
      ZRG     = GZ*(.5*(ZC(I)+ZC(I-IL))-GROUND)
    
C ... Equal or unequal velocities, unequal temperatures. Upwind fluxes (UAV2)

      IF(UAVE >= 0.) THEN     
      F1PIG   = A(I)*UAVE*A1L(IG)*ROHAT !ROIM(IG)
      F1PIA   = A(I)*A1L(IG)*ROHAT !ROIM(IG)
C      F1PIA   = A(I)*AAVE!*ROHAT !ROIM(IG)
      F1PIX   = A(I)*UAVE*YMXA*ROHATC
      F2PIG   = F1PIG*UL + A(I)*AAVE1*PAVE1
      F3PIG   = F1PIG*VL
      F4PIG   = F1PIG*WL
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL
c1, aave c2     F5PIG   = F1PIG*EIL +A(I)*(UAVE1+UC)*AAVE1*(PAVE + FRSPRE + PAVE1)
      F5PIG   = F1PIG*EIL +A(I)*UAV2*AAVE*(PAVE + FRSPRE + PAVE1)
      F6PIG   = F1PIG*RKL
      F7PIG   = F1PIG*REL
      ELSE
      F1PIG   = A(I)*UAVE*A1R(IG)*ROHAT !ROIP(IG)
      F1PIA   = A(I)*A1R(IG)*ROHAT !ROIP(IG)
C      F1PIA   = A(I)*AAVE!*ROHAT !ROIM(IG)
      F1PIX   = A(I)*UAVE*YMXA*ROHATC
      F2PIG   = F1PIG*UR + A(I)*AAVE1*PAVE1
      F3PIG   = F1PIG*VR
      F4PIG   = F1PIG*WR
      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
c1      F5PIG   = F1PIG*EIR +A(I)*(UAVE1+UC)*AAVE1*(PAVE + FRSPRE + PAVE1)
      F5PIG   = F1PIG*EIR +A(I)*UAV2*AAVE*(PAVE + FRSPRE + PAVE1)
      F6PIG   = F1PIG*RKR
      F7PIG   = F1PIG*RER
      ENDIF

      IF(MULPHC == 'CAVIT') F2PIG   = F2PIG + A(I)*AAVE*PAVE
C      F2PIG   = F2PIG + A(I)*AAVE*PAVE

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      A11      = A2XA(I)
      A12      = A2YA(I)
      A13      = A2ZA(I)

      CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
 
c      IF(GRAVIL) THEN ! Test this someday
c         FX    = A(I)*A2XA(I)*XRG*(AAVE*ROHAT-FRSDEN)         
c         FY    = A(I)*A2YA(I)*YRG*(AAVE*ROHAT-FRSDEN)         
c         FZ    = A(I)*A2ZA(I)*ZRG*(AAVE*ROHAT-FRSDEN)
c         F2PM  = F2PM + FX     
c         F3PM  = F3PM + FY       
c         F4PM  = F4PM + FZ
c         FX    = A(I)*XRG*AAVE*ROHAT      
c         FY    = A(I)*YRG*AAVE*ROHAT
c         FZ    = A(I)*ZRG*AAVE*ROHAT
c         F5PIG = F5PIG +.5*(FX*(U(I)+U(I-IL)) + FY*(V(I)+V(I-IL)) +
c     +                      FZ*(W(I)+W(I-IL)))
c      ENDIF
C
C ... FLUXES
      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG +  FE(I)
      FRK(I)  = F6PIG + FRK(I)
      FEPS(I) = F7PIG + FEPS(I)
      VAR(I)%FRO(IPHASE) = F1PIG
      VAR(I)%FA(IPHASE)  = F1PIA
      VAR(I)%FX(IPHASE)  = F1PIX
      VAR(I)%FU(IPHASE)   = UAVE
      VAR(I)%FRM(IPHASE) = F2PM + VAR(I)%FRM(IPHASE)
      VAR(I)%FRN(IPHASE) = F3PM + VAR(I)%FRN(IPHASE)
      VAR(I)%FRW(IPHASE) = F4PM + VAR(I)%FRW(IPHASE)
      VAR(I)%FE(IPHASE)  = F5PIG + VAR(I)%FE(IPHASE)
      PRC(I)%FPM = A(I)*A2XA(I)*PAVE
      PRC(I)%FPN = A(I)*A2YA(I)*PAVE
      PRC(I)%FPW = A(I)*A2ZA(I)*PAVE

c      if(idir == 1 .and. mmm == imax+1) then
c      write(10000+idir,*) mmm,nnn,fro(i),U(I+IL),U(I),U(I-IL),U(I-2*IL)
c      endif
1800  CONTINUE
1900  CONTINUE

C*    HATS SHOULD BE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.94

      DO IG   = 1,IJTRID
      I       = IA + IG
      HAT1(I) = 0.  
      HAT2(I) = FRO(I)/(A(I)+1.E-20)
      HAT3(I) = HAT2(I)
      HAT4(I) = SIGN(1.,FRO(I))*HAT2(I)*HAT4(I)
      ENDDO
2000  CONTINUE

	DEALLOCATE(U,V,W)	

      RETURN
      END SUBROUTINE FLUMPH
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXTF(FRO,FRM,FRN,FRW,FE,VOL,A,RO,RM,RN,RW,P,UG,VG,WG,
     2 E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,EPS2,VIST,A2XA,A2YA,A2ZA,
     3 CP,PRLAM,PRT,UROT,PRO,VAR,TEMP,CH,DRDH,DRDP,ISTATE,
     4 GAMMA,FRSDEN,FRSPRE,HAT1,HAT2,HAT3,HAT4,KSTR,IDIR,MULPHL,NPHASE,
     5 FRESUL,FREDIF,F1H,WAVE,M,IDIFF,NGL)

      USE NS3CO, ONLY : IN, JN, KN
      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : EPS10,EPS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,KSTRID,ICYCLE,
     2   INTEM,IDI,ISTATE,KSTR,IDIR,IL,ILL,IA,IMINID,IMAXID,JMINID,
     3   JMAXID,KMINID,KMAXID,IG,JG,KG,IA2,I,NPHASE,IPHASE,M,IDIFF,
     4   NTOT,NGL

      REAL :: A(*),RO(*),RN(*),RW(*),UG(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),VG(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),WG(*),P(*),UROT(*),TEMP(*),CH(*),CP(*)
     4  ,DRDP(*),DRDH(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*),F1H(*),WAVE(*)

      REAL :: PRLAM,PRS,PRT,GAMMA,FRSDEN,FRSPRE,YPD2,SURG,SURF,SURFE,
     2        DUC,F2VIG,F3VIG,F4VIG,F5VIG,AAVE,FREDIF,SURFE1,SURF4,
     3        AMUP,AMUM,TIF,RLAMP,RLAMM,SURFA
      LOGICAL :: INVIS, MULPHL, FRESUL

      REAL, ALLOCATABLE, DIMENSION(:) :: U,V,W

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      KSTRID  = KMAX + 2*KN	
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN
      NTOT    = ISTRID*JSTRID*KSTRID

	ALLOCATE(U(NTOT),V(NTOT),W(NTOT))	

      PRS     = PRLAM/PRT !INPUT file laminar and turbulent 
      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)


C ... THIN-LAYER VISCOUS FLUXES ARE CALCULATED (LAMINAR AND TURBULENT)
C ... FOR A TWO_FLUID MODEL
        
      IF(IDI /= 0) THEN

            FRM(1:NTOT) = 0.
            FRN(1:NTOT) = 0.
            FRW(1:NTOT) = 0.

      DO IPHASE = 1,NPHASES

        U(1:NTOT) = VAR(1:NTOT)%U(IPHASE)
        V(1:NTOT) = VAR(1:NTOT)%V(IPHASE)
        W(1:NTOT) = VAR(1:NTOT)%W(IPHASE)

      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2
            YPD2   = A(I)/(VOL(I) + VOL(I-IL))	!1/distance = 2A/(V(I)+V(I-IL))
            SURG   = -A(I)*YPD2			!no 2, below calculation of average
C ... EPS2 turbulent over molecylar viscosity plus 1
c           SURF   = SURG*(EPS2(I)*VIS(I)+EPS2(I-IL)*VIS(I-IL)) !=SURG*turbulent_visc?
c            SURF   = SURG*(PRO(I)%VIS(IPHASE)   + VIST(I) + 
c     +                     PRO(I-IL)%VIS(IPHASE)+ VIST(I-IL)) !=SURG*turbulent_visc?

            IF(IDIFF <= 3) THEN ! Aritmethic average
            SURF   = SURG*(PRO(I)%VIS(IPHASE)   + VIST(I)
c     +              *VAR(I)%ALFA(IPHASE)  
     +             +       PRO(I-IL)%VIS(IPHASE)+ VIST(I-IL)) !=SURG*turbulent_visc?
c     +              *VAR(I-IL)%ALFA(IPHASE) )

            ELSE ! Distance-based solution from the geometrical average
            AMUP  = VAR(I)%ALFA(IPHASE)*(PRO(I)%VIS(IPHASE) + VIST(I))
            AMUM  = VAR(I-IL)%ALFA(IPHASE)*(PRO(I-IL)%VIS(IPHASE)+
     +              VIST(I-IL))
            YPD2  =-.5*(VOL(I)*(AMUM+EPS) + VOL(I-IL)*(AMUP+EPS))
            SURF  = A(I)**2*AMUP*AMUM/YPD2 ! Ei pelaa
            ENDIF

            DUC    = .333333*(A2XA(I)*(U(I) - U(I-IL)) +
     +           A2YA(I)*(V(I) - V(I-IL)) + A2ZA(I)*(W(I)-W(I-IL)))
            IF(IDIFF == 2 .OR. IDIFF == 5) DUC = 0. !             

            F2VIG  = SURF*(U(I)-U(I-IL) + A2XA(I)*DUC)
            F3VIG  = SURF*(V(I)-V(I-IL) + A2YA(I)*DUC)
            F4VIG  = SURF*(W(I)-W(I-IL) + A2ZA(I)*DUC)
            F5VIG  = .5*((U(I)+U(I-IL))*F2VIG+(V(I)+V(I-IL))*F3VIG+
     +                   (W(I)+W(I-IL))*F4VIG) 

            IF(IDIFF <= 3) THEN
            AAVE  = .5*(VAR(I)%ALFA(IPHASE) + VAR(I-IL)%ALFA(IPHASE))
c4 nan            AAVE  = 2.*VAR(I)%ALFA(IPHASE)*VAR(I-IL)%ALFA(IPHASE) /
c     +                (VAR(I)%ALFA(IPHASE) + VAR(I-IL)%ALFA(IPHASE)+EPS)
c3            IF(VAR(I)%ALFA(IPHASE) <= 1.E-6 .OR. VAR(I-IL)%ALFA(IPHASE)
c3     +       <= 1.E-6) AAVE = 0.
            VAR(I)%FRM(IPHASE) = AAVE*F2VIG
            VAR(I)%FRN(IPHASE) = AAVE*F3VIG
            VAR(I)%FRW(IPHASE) = AAVE*F4VIG
            F5VIG              = AAVE*F5VIG

C ... Sum of the individual fluxes for checking

            FRM(I) = FRM(I) + AAVE*F2VIG
            FRN(I) = FRN(I) + AAVE*F3VIG
            FRW(I) = FRW(I) + AAVE*F4VIG

            ELSEIF(IDIFF >= 4) THEN

            VAR(I)%FRM(IPHASE) = F2VIG
            VAR(I)%FRN(IPHASE) = F3VIG
            VAR(I)%FRW(IPHASE) = F4VIG

            FRM(I) = FRM(I) + F2VIG
            FRN(I) = FRN(I) + F3VIG
            FRW(I) = FRW(I) + F4VIG

            ENDIF ! IDIFF <= 3


c            FE(I)  = FE(I) + F5VIG + SURFE*(TEMP(I) - TEMP(I-IL))

C                    =====
C ... POSSIBLE DIFFUSION OF TURBULENT KINECTIC ENERGY PPR 10.8.1995
C

               IF(IDIFF <= 6) THEN ! Aritmethic average
               SURFE = SURG*(PRO(I)%CP(IPHASE)*(EPS2(I)-1.)*
     +         PRO(I)%VIS(IPHASE)/PRT + 
     +         PRO(I)%CH(IPHASE) + PRO(I-IL)%CP(IPHASE)*
     +         (EPS2(I-IL)-1.)*PRO(I-IL)%VIS(IPHASE)/PRT + 
     +         PRO(I-IL)%CH(IPHASE))
               TIF    = .5*(PRO(I)%TSAT+PRO(I-IL)%TSAT)
               VAR(I)%FE(IPHASE) = F5VIG +  
     +         SURFE*(VAR(I)%ALFA(IPHASE)*PRO(I)%DTEMP(IPHASE) -
     +         VAR(I-IL)%ALFA(IPHASE)*PRO(I-IL)%DTEMP(IPHASE)  -
     +         TIF*(VAR(I)%ALFA(IPHASE) - VAR(I-IL)%ALFA(IPHASE)))

               ELSE ! Equally Distance-based solution testaamati, works inte
               RLAMP = PRO(I)%CH(IPHASE) + PRO(I)%CP(IPHASE)*VIST(I)/PRT
               AMUP  = VAR(I)%ALFA(IPHASE)*RLAMP
               RLAMM = PRO(I-IL)%CH(IPHASE) + 
     +                 PRO(I-IL)%CP(IPHASE)*VIST(I-IL)/PRT
               AMUM  = VAR(I-IL)%ALFA(IPHASE)*RLAMM
               YPD2  =-.5*(VOL(I)*(AMUM+EPS) + VOL(I-IL)*(AMUP+EPS))
               SURFE = A(I)**2*AMUP*AMUM/YPD2
               AMUP  = (PRO(I)%TSAT - PRO(I)%DTEMP(IPHASE))*RLAMP
               AMUM  = (PRO(I-IL)%TSAT -PRO(I-IL)%DTEMP(IPHASE))*RLAMM
               YPD2  =-.5*(VOL(I)*(AMUM+EPS) + VOL(I-IL)*(AMUP+EPS))
               SURFA = A(I)**2*AMUP*AMUM/YPD2

c               TIF    = .5*(PRO(I)%TSAT+PRO(I-IL)%TSAT)
c               TIF    = .25*(PRO(I)%DTEMP(1)+PRO(I-IL)%DTEMP(1)
c     2                     + PRO(I)%DTEMP(2)+PRO(I-IL)%DTEMP(2))

               VAR(I)%FE(IPHASE) = F5VIG +  
     +         SURFE*(PRO(I)%DTEMP(IPHASE) -PRO(I-IL)%DTEMP(IPHASE))  -
     +         SURFA*(VAR(I)%ALFA(IPHASE) - VAR(I-IL)%ALFA(IPHASE))
               ENDIF

c               VAR(I)%FRO(IPHASE) = SURF*(VAR(I)%ALFA(IPHASE) - ! Artificial
c     +           VAR(I)%ALFA(IPHASE))                           ! diffusion

c     +         SURFE * (PRO(I)%TEMP(IPHASE)-PRO(I-IL)%TEMP(IPHASE))) !+
c     +         SURFE * (VAR(I)%ALFA(IPHASE)-VAR(I-IL)%ALFA(IPHASE)*
c     +         (PRO(I)%TSAT+PRO(I-IL)%TSAT - PRO(I)%TEMP(IPHASE) -
c     +         PRO(I-IL)%TEMP(IPHASE)) *.5)
C ... above does not work yet

 850     CONTINUE
 900  CONTINUE
      ENDDO

C ... Correct tiny energy fluxes

      DO KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO IG = 1,IMAXID
            I      = IG + IA2

               IF(VAR(I)%FE(1) > 0.) THEN
                  IF(VAR(I)%ALFA(1) < EPS10) THEN
                  VAR(I)%FE(2) = VAR(I)%FE(2) + VAR(I)%FE(1)
                  VAR(I)%FE(1) = 0.
                  ENDIF
               ELSE IF(VAR(I)%FE(1) < 0.) THEN
                  IF(VAR(I-IL)%ALFA(1) < EPS10) THEN
                  VAR(I)%FE(2) = VAR(I)%FE(2) + VAR(I)%FE(1)
                  VAR(I)%FE(1) = 0.
                  ENDIF
               ENDIF

               IF(VAR(I)%FE(2) > 0.) THEN
                  IF(VAR(I)%ALFA(2) < EPS10) THEN
                  VAR(I)%FE(1) = VAR(I)%FE(1) + VAR(I)%FE(2)
                  VAR(I)%FE(2) = 0.
                  ENDIF
               ELSE IF(VAR(I)%FE(2) < 0.) THEN
                  IF(VAR(I-IL)%ALFA(2) < EPS10) THEN
                  VAR(I)%FE(1) = VAR(I)%FE(1) + VAR(I)%FE(2)
                  VAR(I)%FE(2) = 0.
                  ENDIF
               ENDIF

               FE(I) = VAR(I)%FE(1) + VAR(I)%FE(2)

      ENDDO; ENDDO; ENDDO

      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION COULD BE INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
c         FRO(I)  = 0.
         FRM(I)  = 0.
         FRN(I)  = 0.
         FRW(I)  = 0.
         FE(I)   = 0.
         VAR(I)%FE(1) = 0.
         VAR(I)%FE(2) = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

	DEALLOCATE(U,V,W)	

      RETURN
      END SUBROUTINE FLUXTF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C	*********** MODIFIED 7.6.2012 by Juho N ********************
      SUBROUTINE FLUXDO(FRO,FRM,FRN,FRW,FE,VOL,A,RO,RM,RN,RW,P,U,V,W,E,
     2 C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,EPS2,VIST,A2XA,A2YA,A2ZA,
     3 CP,PRLAM,PRT,UROT,PRO,VAR,TEMP,CH,DRDH,DRDP,ISTATE,
     4 GAMMA,FRSDEN,FRSPRE,HAT1,HAT2,HAT3,HAT4,KSTR,IDIR,MULPHL,NPHASE,
     5 FRESUL,FREDIF,F1H,WAVE,M,IDIFF)

      USE NS3CO, ONLY : IN, JN, KN
      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : EPS10

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,KSTRID,ICYCLE,
     2   INTEM,IDI,ISTATE,KSTR,IDIR,IL,ILL,IA,IMINID,IMAXID,JMINID,
     3   JMAXID,KMINID,KMAXID,IG,JG,KG,IA2,I,NPHASE,IPHASE,M,IDIFF

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),UROT(*),TEMP(*),CH(*),CP(*)
     4  ,DRDP(*),DRDH(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*),F1H(*),WAVE(*)

      REAL :: PRLAM,PRS,PRT,GAMMA,FRSDEN,FRSPRE,YPD2,SURG,SURF,SURFE,
     2        DUC,F2VIG,F3VIG,F4VIG,F5VIG,AAVE,FREDIF,SURFE1,SURF4
      LOGICAL :: INVIS, MULPHL, FRESUL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      KSTRID  = KMAX + 2*KN	
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN		

      PRS     = PRLAM/PRT !INPUT file laminar and turbulent 
      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)


C ... THIN-LAYER VISCOUS FLUXES ARE CALCULATED (LAMINAR AND TURBULENT)
        
      IF(IDI /= 0) THEN

      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2
            YPD2   = A(I)/(VOL(I) + VOL(I-IL))	!1/distance = 2A/(V(I)+V(I-IL))
            SURG   = -A(I)*YPD2			!no 2 because below calculation of average
! EPS2 turbulent over molecylar viscosity plus 1
            SURF   = SURG*(EPS2(I)*VIS(I)+EPS2(I-IL)*VIS(I-IL)) !=SURG*turbulent_visc?
            SURFE  = SURG*(CP(I)*VIST(I)/PRT+CH(I) + 
     +                     CP(I-IL)*VIST(I-IL)/PRT+CH(I-IL))
            IF(IDIFF <= 3) THEN ! Aritmethic average
            SURF   = SURG*(EPS2(I)*VIS(I)+EPS2(I-IL)*VIS(I-IL))	!+ because average
				  !CP = specific heat
            SURFE  = SURG*(CP(I)*VIST(I)/PRT+CH(I) + 
     &                      CP(I-IL)*VIST(I-IL)/PRT+CH(I-IL))
            ELSE ! Equally Distance-based solution
            SURF4  = .5/(VIS(I)+VIST(I)) + .5/(VIST(I-IL)+VIS(I-IL))
            SURF   = 2.*SURG/SURF4
            SURFE1 = .5/(CP(I)*VIST(I)/PRT+CH(I)) +
     &               .5/(CP(I-IL)*VIST(I-IL)/PRT+CH(I-IL))
            SURFE  = 2.*SURG/SURFE1
            ENDIF

            DUC    = .333333*(A2XA(I)*(U(I) - U(I-IL)) +
     +           A2YA(I)*(V(I) - V(I-IL)) + A2ZA(I)*(W(I)-W(I-IL)))
            IF(IDIFF == 2 .OR. IDIFF == 5) DUC = 0. !             

            F2VIG  = SURF*(U(I)-U(I-IL) + A2XA(I)*DUC)
            F3VIG  = SURF*(V(I)-V(I-IL) + A2YA(I)*DUC)
            F4VIG  = SURF*(W(I)-W(I-IL) + A2ZA(I)*DUC)
            F5VIG  = .5*((U(I)+U(I-IL))*F2VIG+(V(I)+V(I-IL))*F3VIG+
     +       (W(I)+W(I-IL))*F4VIG) 
            
            FRM(I) = F2VIG
            FRN(I) = F3VIG
            FRW(I) = F4VIG
            FE(I)  = FE(I) + F5VIG + SURFE*(TEMP(I) - TEMP(I-IL))
C                    =====
C ... POSSIBLE DIFFUSION OF TURBULENT KINECTIC ENERGY PPR 10.8.1995
C
            IF(MULPHL) THEN ! Currently only approximative

               DO IPHASE = 1,NPHASES
               AAVE  = .5*(VAR(I)%ALFA(IPHASE) + VAR(I-IL)%ALFA(IPHASE))
               IF(IDIFF <= 3) THEN ! Aritmethic average
               SURFE = SURG*(PRO(I)%CP(IPHASE)*(EPS2(I)-1.)*
     +         PRO(I)%VIS(IPHASE)/PRT + 
     +         PRO(I)%CH(IPHASE) + PRO(I-IL)%CP(IPHASE)*
     +         (EPS2(I-IL)-1.)*PRO(I-IL)%VIS(IPHASE)/PRT + 
     +         PRO(I-IL)%CH(IPHASE))
               ELSE ! Equally Distance-based solution testaamati
               SURFE1 = .5/(PRO(I)%CP(IPHASE)*VIST(I)/PRT +
     +                  .5/(PRO(I-IL)%CP(IPHASE)*VIST(I-IL)/PRT +
     +                  .5/PRO(I)%CH(IPHASE)) + .5/PRO(I-IL)%CH(IPHASE))
     +                    
               SURFE  = 2.*SURG/SURFE1
               ENDIF

               VAR(I)%FE(IPHASE) = AAVE * F5VIG +  
     +         SURFE*(VAR(I)%ALFA(IPHASE)*PRO(I)%DTEMP(IPHASE) -
     +         VAR(I-IL)%ALFA(IPHASE)*PRO(I-IL)%DTEMP(IPHASE)  -
     +         .5*(PRO(I)%TSAT+PRO(I-IL)%TSAT)*
     +         (VAR(I)%ALFA(IPHASE) - VAR(I-IL)%ALFA(IPHASE)))

c               VAR(I)%FRO(IPHASE) = SURF*(VAR(I)%ALFA(IPHASE) - ! Artificial
c     +           VAR(I)%ALFA(IPHASE))                           ! diffusion

c     +         SURFE * (PRO(I)%TEMP(IPHASE)-PRO(I-IL)%TEMP(IPHASE))) !+
c     +         SURFE * (VAR(I)%ALFA(IPHASE)-VAR(I-IL)%ALFA(IPHASE)*
c     +         (PRO(I)%TSAT+PRO(I-IL)%TSAT - PRO(I)%TEMP(IPHASE) -
c     +         PRO(I-IL)%TEMP(IPHASE)) *.5)
C ... above does not work yet
               ENDDO

               IF(VAR(I)%FE(1) > 0.) THEN
                  IF(VAR(I)%ALFA(1) < EPS10) THEN
                  VAR(I)%FE(2) = VAR(I)%FE(2) + VAR(I)%FE(1)
                  VAR(I)%FE(1) = 0.
                  ENDIF
               ELSE IF(VAR(I)%FE(1) < 0.) THEN
                  IF(VAR(I-IL)%ALFA(1) < EPS10) THEN
                  VAR(I)%FE(2) = VAR(I)%FE(2) + VAR(I)%FE(1)
                  VAR(I)%FE(1) = 0.
                  ENDIF
               ENDIF

               IF(VAR(I)%FE(2) > 0.) THEN
                  IF(VAR(I)%ALFA(2) < EPS10) THEN
                  VAR(I)%FE(1) = VAR(I)%FE(1) + VAR(I)%FE(2)
                  VAR(I)%FE(2) = 0.
                  ENDIF
               ELSE IF(VAR(I)%FE(2) < 0.) THEN
                  IF(VAR(I-IL)%ALFA(2) < EPS10) THEN
                  VAR(I)%FE(1) = VAR(I)%FE(1) + VAR(I)%FE(2)
                  VAR(I)%FE(2) = 0.
                  ENDIF
               ENDIF

                FE(I) = VAR(I)%FE(1) + VAR(I)%FE(2)
cc               VAR(I)%FE(1) = AAVE*FE(I)
cc               VAR(I)%FE(2) = (1- AAVE)*FE(I)
c              IF(AAVE < 0.9999999) THEN
c               VAR(I)%FE(1) = AAVE*FE(I)
c               VAR(I)%FE(2) = (1- AAVE)*FE(I)
c               ELSE
c               VAR(I)%FE(1) = FE(I)
c               VAR(I)%FE(2) = 0.
c               ENDIF
            ENDIF ! MULPHL

            IF(FRESUL .AND. M == 1) THEN ! Mahdollisesti voi pudottaa pois?
                F1H(I) = 0.25*SURG*(RO(I)+RO(I-IL))*FREDIF*!kaksi nollaa
     +                   (WAVE(I) -WAVE(I-IL)) ! p.o. 0.25, WIGLEY kaatui
c     ++0.005*(WAVE(I+IL)-3.*WAVE(I)+3.*WAVE(I-IL)-WAVE(I-2*IL))!Merkki testattu
            ENDIF


 850     CONTINUE
 900  CONTINUE
      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
c         FRO(I)  = 0.
         FRM(I)  = 0.
         FRN(I)  = 0.
         FRW(I)  = 0.
         FE(I)   = 0.
         IF(MULPHL) THEN
            VAR(I)%FE(1) = 0.
            VAR(I)%FE(2) = 0.
         ENDIF
         IF(FRESUL) THEN
            F1H(I) = 0.
         ENDIF
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLUXDO 
C	*********** MOD ENDS 7.6.2012 by Juho N ********************
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

      LOGICAL :: INVIS

C ... SET THE CONTROL PARAMETER FOR A VISCOUS SWEEP

      IMINID  = 1
      IMAXID  = IMAX
      JMINID  = 1
      JMAXID  = JMAX
      KMINID  = 1
      KMAXID  = KMAX
      INVIS   = .TRUE.

      IF(IDIR == 1) THEN
         IMAXID = MIN(IDI,IMAX+1)
         IMINID = IMAXID + 1
         IF(IMAXID == (IMAX+1)) INVIS = .FALSE.
      ENDIF
      IF(IDIR == 2) THEN
         JMAXID = MIN(IDI,JMAX+1)
         JMINID = JMAXID + 1
         IF(JMAXID == (JMAX+1)) INVIS = .FALSE.
      ENDIF
      IF(IDIR == 3) THEN
         KMAXID = MIN(IDI,KMAX+1)
         KMINID = KMAXID + 1
         IF(KMAXID == (KMAX+1)) INVIS = .FALSE.
      ENDIF

      RETURN
      END SUBROUTINE SINVIS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C


C ******************************
C *    EXACT DIFFUSION TERMS   *
C ******************************

      SUBROUTINE FLUXNS(FRO,FRM,FRN,FRW,FE,VOL,A,RO,RM,RN,RW,P,U,V,W,E,
     2 C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,EPS2,VIST,
     3 W12,W13,W23,S11,S12,S13,S22,S23,S33,DTDX,DTDY,DTDZ, 
     4 A2XA,A2YA,A2ZA,PRLAM,PRT,UROT,TEMP,CH,DRDH,DRDP,
     5 ISTATE,GAMMA,FRSDEN,FRSPRE,HAT1,HAT2,HAT3,HAT4,KSTR,IDIR,
     6 ITURB) 

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*),
     2  FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*),
     3  VOL(*),A2ZA(*),FRW(*),W(*),P(*),UROT(*),TEMP(*),CH(*),
     4  DRDP(*),DRDH(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*),
     5  W12(*),W13(*),W23(*),S11(*),S12(*),S13(*),S22(*),S23(*),S33(*), 
     6  DTDX(*),DTDY(*),DTDZ(*) 
      LOGICAL INVIS

C     HUOM ! LAMPOTILA-DERIVAATAT DTDX, DTDY JA DTDZ TUODAAN VOISSA
C     W12,W13,W23 
C     Old unused version

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN

      R23     = .6666667
      PRS     = PRLAM/PRT
      IL      = KSTR


      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)


c      open(669,file='incns')
C ... COMPLETE VISCOUS FLUXES ARE CALCULATED (LAMINAR AND TURBULENT)

      IF(IDI /= 0) THEN
      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2
            DUDX   = 0.5*(S11(I) + S11(I-IL))
            DVDY   = 0.5*(S22(I) + S22(I-IL))
            DWDZ   = 0.5*(S33(I) + S33(I-IL))
c         call ijkpai(I,imax,jmax,kmax,mmm,nnn,lll)
c         write(669,*) mmm,nnn,R23*(DUDX+DVDY+DWDZ)*A2XA(I),
c     &   R23*(DUDX+DVDY+DWDZ)*A2YA(I)     
c            DUDY   = 0.5*(S12(I) + W12(I) + S12(I-IL) + W12(I-IL))
c            DUDZ   = 0.5*(S13(I) + W13(I) + S13(I-IL) + W13(I-IL))
c            DVDX   = 0.5*(S12(I) - W12(I) + S12(I-IL) - W12(I-IL))
c            DVDZ   = 0.5*(S23(I) + W23(I) + S23(I-IL) + W23(I-IL))
c            DWDX   = 0.5*(S13(I) - W13(I) + S13(I-IL) - W13(I-IL))
c            DWDY   = 0.5*(S23(I) - W23(I) + S23(I-IL) - W23(I-IL)) 
            DTDXI   = 0.5*(DTDX(I) + DTDX(I-IL))
            DTDYI   = 0.5*(DTDY(I) + DTDY(I-IL))
            DTDZI   = 0.5*(DTDZ(I) + DTDZ(I-IL))
cc       if(kg == kmaxid .or. jg == jmaxid .or. ig == imaxid)
cc     & write(77,*)ig,jg,kg,dtdx(i),dtdx(i-il),dtdy(i),dtdy(i-il)      
C      write(67,*) W12(I),W12(I-IL),W13(I),W13(I-IL),W23(I),W23(I-IL)

            
            SURF   = 0.5*A(I)*(EPS2(I)*VIS(I)+EPS2(I-IL)*VIS(I-IL))
            SURFE  = 0.5*A(I)*((1.+PRS*VIST(I)/VIS(I))*CH(I) + 
     &           (1.+PRS*VIST(I-IL)/VIS(I-IL))*CH(I-IL))

            F2VIG  = SURF*(R23*(2.*DUDX - DVDY - DWDZ)*A2XA(I) +
     &      (S12(I)+S12(I-IL))*A2YA(I) + (S13(I)+S13(I-IL))*A2ZA(I))

            F3VIG  = SURF*((S12(I)+S12(I-IL))*A2XA(I) + 
     &      R23*(2.*DVDY - DUDX - DWDZ)*A2YA(I) +
     &      (S23(I)+S23(I-IL))*A2ZA(I))

            F4VIG  = SURF*((S13(I)+S13(I-IL))*A2XA(I) +
     &      (S23(I)+S23(I-IL))*A2YA(I) + 
     &      R23*(2.*DWDZ - DUDX - DVDY)*A2ZA(I))

            F5VIG  = .5*((U(I) + U(I-IL))*F2VIG +
     &      (V(I) + V(I-IL))*F3VIG + (W(I) + W(I-IL))*F4VIG) +
     &      SURFE*(DTDXI*A2XA(I) + DTDYI*A2YA(I) + DTDZI*A2ZA(I))      

c      if(idir == 1) then      
c      write(67,*) -SURFE*(DTDXI*A2XA(I)+DTDYI*A2YA(I)+DTDZI*A2ZA(I))
c      else
c      write(68,*) -SURFE*(DTDXI*A2XA(I)+DTDYI*A2YA(I)+DTDZI*A2ZA(I))
c      endif
c            F2VIG  = SURF*(R23*(2.*DUDX - DVDY - DWDZ)*A2XA(I) +
c     &      (DUDY + DVDX)*A2YA(I) + (DWDX + DUDZ)*A2ZA(I))
c
c            F3VIG  = SURF*((DUDY + DVDX)*A2XA(I) + 
c     &      R23*(2.*DVDY - DUDX - DWDZ)*A2YA(I) +
c     &      (DVDZ + DWDY)*A2ZA(I))
c
c            F4VIG  = SURF*((DWDX + DUDZ)*A2XA(I) +
c     &      (DVDZ + DWDY)*A2YA(I) + 
c     &      R23*(2.*DWDZ - DUDX - DVDY)*A2ZA(I))
c
c            F5VIG  = .5*((U(I) + U(I-IL))*F2VIG +
c     &      (V(I) + V(I-IL))*F3VIG + (W(I) + W(I-IL))*F4VIG) +
c     &      (2.*A(I)/(VOL(I)+VOL(I-IL)))*SURFE*(TEMP(I) - TEMP(I-IL))

            FRO(I) = 0.
            FRM(I) = -F2VIG
            FRN(I) = -F3VIG
            FRW(I) = -F4VIG
            FE(I)  = -F5VIG
            
c      call ijkpai(I,imax,jmax,kmax,mm,nn,ll)
c      if(mm  == 17)then
c      write(67,*) mm,nn,ll
c      write(66,*)  FRM(I),FRN(I),FE(I)
c      endif
C      write(67,*) A(I)*(DTDX*A2XA(I)+DTDY*A2YA(I)+DTDZ*A2ZA(I)),
C     & FRM(I),FRN(I),FE(I)
C      write(67,*) DTDX,DTDY,DTDZ
  
C            FE(I)  = FE(I) - F5VIG
C                    =====
C ... POSSIBLE DIFFUSION OF TURBULENT KINECTIC ENERGY PPR 10.8.1995
C
 850     CONTINUE
 900  CONTINUE
      ENDIF

c tulostus
c      DO 500 K = 1,KMAX
c      IA  = (KN+K-1)*ILL + JN*ISTRID
c
c      DO 700 JJ = 1,JMAX*ISTRID
c         L  = IA + JJ
c         write(67,*) W12(L),W13(L),W23(L)
c         
c 700     continue
c 500     continue
c
c tulostus
c      close(669)
      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
         FRO(I)  = 0.
         FRM(I)  = 0.
         FRN(I)  = 0.
         FRW(I)  = 0.
         FE(I)   = 0.
 950  CONTINUE
 1000 CONTINUE

      ENDIF

      RETURN
      END SUBROUTINE FLUXNS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXBL(FRO,FRM,FRN,FRW,FE,VOL,A,XC,YC,ZC,RO,RM,RN,RW,
     2 P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,EPS2,VIST,
     3 W12,W13,W23,S11,S12,S13,S22,S23,S33,DTDX,DTDY,DTDZ,D2,A2XA,A2YA, 
     4 A2ZA,CP,PRLAM,PRT,UROT,PRO,VAR,TEMP,CH,DRDH,DRDP,ISTATE,
     5 GAMMA,FRSDEN,FRSPRE,HAT1,HAT2,HAT3,HAT4,KSTR,IDIR,ITURB,MULPHL,
     6 IDIFF)
cc     7 dudxi,dudyi,dudzi,dvdxi,dvdyi,dvdzi,dwdxi,dwdyi,dwdzi) 

      USE NS3CO, ONLY : IN, JN, KN     
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,ISTATE,IDIR,
     2  ITURB,ISTRID,JSTRID,KSTRID,IA,ILL,IL,IMINID,IMAXID,JMINID,
     3  JMAXID,KMINID,KMAXID,IG,JG,KG,IA2,I,KSTR,IDIFF

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*),
     2  FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*),
     3  VOL(*),A2ZA(*),FRW(*),W(*),P(*),UROT(*),TEMP(*),CH(*),CP(*),
     4  DRDP(*),DRDH(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*),
     5  W12(*),W13(*),W23(*),S11(*),S12(*),S13(*),S22(*),S23(*),S33(*), 
     6  DTDX(*),DTDY(*),DTDZ(*),D2(*)
c      real dudxi(*),dudyi(*),dudzi(*),dvdxi(*),dvdyi(*),dvdzi(*),
c     2 dwdxi(*),dwdyi(*),dwdzi(*)

      REAL :: XC(*), YC(*), ZC(*) 

      REAL :: PRLAM,PRT,GAMMA,FRSDEN,FRSPRE,R13,R23,PRS,DUDX,DUDY,DUDZ,
     2  DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,S11I,S12I,S13I,S22I,S23I,S33I,
     3  DTDXI,DTDYI,DTDZI,SXI,SYI,SZI,SLEN,ALFA,SURF,SURF2,SURFE,
     4  DUDXSJ,DVDXSJ,DWDXSJ,SDOTDU,SDOTDV,SDOTDW,ADOTU,F2VIB,F3VIB,
     5  F4VIB,F5VIB,DUVWSJ,W1,W2,AAVE,SURFE1,SURF3,SURF4,W1P,W2P

      LOGICAL :: INVIS,MULPHL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

C     Viscous flux calculation using Blazek formulae

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN
     
      R13     = .3333333
      R23     = .6666667
      PRS     = PRLAM/PRT
      IL      = KSTR
       
      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

C ... COMPLETE VISCOUS FLUXES ARE CALCULATED (LAMINAR AND TURBULENT)
      
      IF(IDI /= 0) THEN

      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2

C ... Distance weighting, although it does not has any effect

            W1  = 2.0*( D2(I) / ( D2(I) + D2(I-IL) ) ) 
            W2  = 2.0 - W1
            W1P = .5*W1
            W2P = .5*W2

            DUDX   = W2P*S11(I) + W1P*S11(I-IL)
            DVDY   = W2P*S22(I) + W1P*S22(I-IL)
            DWDZ   = W2P*S33(I) + W1P*S33(I-IL)

            DUDY   = W2P*(S12(I) + W12(I)) + 
     &               W1P*(S12(I-IL) + W12(I-IL))
            DUDZ   = W2P*(S13(I) + W13(I)) + 
     &               W1P*(S13(I-IL) + W13(I-IL))

            DVDX   = W2P*(S12(I) - W12(I)) + 
     &               W1P*(S12(I-IL) - W12(I-IL))
            DVDZ   = W2P*(S23(I) + W23(I)) + 
     &               W1P*(S23(I-IL) + W23(I-IL))

            DWDX   = W2P*(S13(I) - W13(I)) + 
     &               W1P*(S13(I-IL) - W13(I-IL))
            DWDY   = W2P*(S23(I) - W23(I)) + 
     &               W1P*(S23(I-IL) - W23(I-IL))

            S11I   = W2*S11(I) + W1*S11(I-IL)
            S12I   = W2*S12(I) + W1*S12(I-IL)
            S13I   = W2*S13(I) + W1*S13(I-IL)
            S22I   = W2*S22(I) + W1*S22(I-IL)
            S23I   = W2*S23(I) + W1*S23(I-IL)
            S33I   = W2*S33(I) + W1*S33(I-IL)

            DTDXI  = W2P*DTDX(I) + W1P*DTDX(I-IL)
            DTDYI  = W2P*DTDY(I) + W1P*DTDY(I-IL)
            DTDZI  = W2P*DTDZ(I) + W1P*DTDZ(I-IL)

C ... test weather this is the same as XFC, YFC and ZFC

            SXI    = XC(I) - XC(I-IL)
            SYI    = YC(I) - YC(I-IL)
            SZI    = ZC(I) - ZC(I-IL)
            SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
            SXI    = SXI/SLEN
            SYI    = SYI/SLEN
            SZI    = SZI/SLEN
            ALFA   = SXI*A2XA(I) + SYI*A2YA(I) + SZI*A2ZA(I)
c            ALFA   = MAX(ALFA,0.)
            
            IF(IDIFF <= 3) THEN ! Aritmethic average
            SURF   = 0.5*A(I)*(EPS2(I)*VIS(I)+EPS2(I-IL)*VIS(I-IL))
            SURF2  = 0.5*A(I)*(VIS(I)+VIST(I)+VIST(I-IL)+VIS(I-IL))
            SURFE  = 0.5*A(I)*(CP(I)*VIST(I)/PRT+CH(I) + 
     &                      CP(I-IL)*VIST(I-IL)/PRT+CH(I-IL))
c            SURF3  = A(I)*(W2P*(VIS(I)+VIST(I)) +
c     &                     W1P*(VIST(I-IL)+VIS(I-IL)))
            ELSE ! Distance based solution
            SURF4  = W1P/(VIS(I)+VIST(I)) + W2P/(VIST(I-IL)+VIS(I-IL))
            SURF   = A(I)/SURF4
            SURFE1 = W1P/(CP(I)*VIST(I)/PRT+CH(I)) +
     &               W2P/(CP(I-IL)*VIST(I-IL)/PRT+CH(I-IL))
            SURFE  = A(I)/SURFE1
            ENDIF

            DUDXSJ = DUDX*SXI + DUDY*SYI + DUDZ*SZI     
            DVDXSJ = DVDX*SXI + DVDY*SYI + DVDZ*SZI     
            DWDXSJ = DWDX*SXI + DWDY*SYI + DWDZ*SZI 
            DUVWSJ = A2XA(I)*DUDXSJ+A2YA(I)*DVDXSJ+A2ZA(I)*DWDXSJ

            SDOTDU = SXI*DUDX + SYI*DUDY + SZI*DUDZ
            SDOTDV = SXI*DVDX + SYI*DVDY + SZI*DVDZ
            SDOTDW = SXI*DWDX + SYI*DWDY + SZI*DWDZ
            ADOTU  = (A2XA(I)*(U(I)-U(I-IL))+ A2YA(I)*(V(I)-V(I-IL))    
     &      +         A2ZA(I)*(W(I)-W(I-IL)))

            F2VIB  = A2XA(I)*S11I + A2YA(I)*S12I + A2ZA(I)*S13I
            F3VIB  = A2XA(I)*S12I + A2YA(I)*S22I + A2ZA(I)*S23I
            F4VIB  = A2XA(I)*S13I + A2YA(I)*S23I + A2ZA(I)*S33I

            IF(ALFA > 0. .AND. IDIFF /= 3  .OR. ! Correct the flux
     &         ALFA > 0. .AND. IDIFF /= 6) THEN ! (Blazek)

            F2VIB  = ALFA * ((U(I)-U(I-IL))/SLEN - SDOTDU)
     &      + SXI*(ADOTU/SLEN - DUVWSJ) + F2VIB  
            F3VIB  = ALFA * ((V(I)-V(I-IL))/SLEN - SDOTDV)
     &      + SYI*(ADOTU/SLEN - DUVWSJ) + F3VIB 
            F4VIB  = ALFA * ((W(I)-W(I-IL))/SLEN - SDOTDW) 
     &      + SZI*(ADOTU/SLEN - DUVWSJ) + F4VIB

            ENDIF

            IF(IDIFF /= 2 .AND. IDIFF /= 5) THEN ! Add lambda viscosity
            F2VIB  = F2VIB - R23*A2XA(I)*(DUDX + DVDY + DWDZ)
            F3VIB  = F3VIB - R23*A2YA(I)*(DUDX + DVDY + DWDZ)
            F4VIB  = F4VIB - R23*A2ZA(I)*(DUDX + DVDY + DWDZ)
            ENDIF

            F2VIB  = SURF * F2VIB
            F3VIB  = SURF * F3VIB
            F4VIB  = SURF * F4VIB

            F5VIB  = .5*((U(I) + U(I-IL))*F2VIB +
     &      (V(I) + V(I-IL))*F3VIB + (W(I) + W(I-IL))*F4VIB) +
     &      SURFE*(DTDXI*A2XA(I) + DTDYI*A2YA(I) + DTDZI*A2ZA(I))

            IF(ALFA > 0. .AND. IDIFF /= 3  .OR. ! Correct the flux
     &         ALFA > 0. .AND. IDIFF /= 6) THEN ! (Blazek)

            F5VIB  = SURFE * ALFA * ((TEMP(I)-TEMP(I-IL))/SLEN - 
     &               SXI*DTDXI - SYI*DTDYI - SZI*DTDZI) + F5VIB
            ENDIF

c           FRO(I) = 0.
            FRM(I) = -F2VIB
            FRN(I) = -F3VIB
            FRW(I) = -F4VIB 
c           FE(I)  = -F5VIB
              
            FE(I)  = FE(I) - F5VIB 
C                    =====
C ... POSSIBLE DIFFUSION OF TURBULENT KINECTIC ENERGY PPR 10.8.1995
C
            IF(MULPHL) THEN ! Currently approximative (skeidaa, do not use)
               AAVE = .5*(VAR(I)%ALFA(1) + VAR(I-IL)%ALFA(1))
               VAR(I)%FE(1) = AAVE*FE(I)
               VAR(I)%FE(2) = (1.- AAVE)*FE(I)
            ENDIF
 850     CONTINUE
 900  CONTINUE
      ENDIF

c tulostus
c      DO 500 K = 1,KMAX
c      IA  = (KN+K-1)*ILL + JN*ISTRID
c
c      DO 700 JJ = 1,JMAX*ISTRID
c         L  = IA + JJ
c         write(67,*) W12(L),W13(L),W23(L)
c         
c 700     continue
c 500     continue
c
c tulostus
c      close(69)

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED

      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
c         FRO(I)  = 0.
         FRM(I)  = 0.
         FRN(I)  = 0.
         FRW(I)  = 0.
         FE(I)   = 0.
         IF(MULPHL) THEN ! Currently approximative
            VAR(I)%FE(1) = 0.
            VAR(I)%FE(2) = 0.
         ENDIF
 950  CONTINUE
 1000 CONTINUE

      ENDIF

      RETURN
      END SUBROUTINE FLUXBL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXHS(FE,A2,XC,YC,ZC,E,IMAX,JMAX,KMAX,
     2 ICYCLE,INTEM,IDI,DTDX,DTDY,DTDZ,D2,A2XA,A2YA,A2ZA, 
     4 TEMP,CH,KSTR,IDIR)

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      REAL :: A2(*),E(*),FE(*),A2XA(*),A2YA(*),A2ZA(*),
     1  TEMP(*),CH(*),DTDX(*),DTDY(*),DTDZ(*),D2(*)

      REAL :: XC(*), YC(*), ZC(*) 

      REAL :: DTDXI,DTDYI,DTDZI,SXI,SYI,SZI,SLEN,ALFA,SURFE,W1,W2,alpo

      INTEGER :: IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,ISTATE,IDIR,
     2  ITURB,ISTRID,JSTRID,KSTRID,IA,ILL,IL,IMINID,IMAXID,JMINID,
     3  JMAXID,KMINID,KMAXID,IG,JG,KG,IA2,I,KSTR

      LOGICAL :: INVIS

C     Heat flux calculation in solids using Blazek formulae

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN
      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

      IF(IDI /= 0) THEN

      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2

C ... Distance weighting, although it does not has any effect

            W1  = 2.0*( D2(I) / ( D2(I) + D2(I-IL) ) ) 
            W2  = 2.0 - W1

            DTDXI  = 0.5*(W2*DTDX(I) + W1*DTDX(I-IL))
            DTDYI  = 0.5*(W2*DTDY(I) + W1*DTDY(I-IL))
            DTDZI  = 0.5*(W2*DTDZ(I) + W1*DTDZ(I-IL))

            SXI    = XC(I) - XC(I-IL)
            SYI    = YC(I) - YC(I-IL)
            SZI    = ZC(I) - ZC(I-IL)
            SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
            SXI    = SXI/SLEN
            SYI    = SYI/SLEN
            SZI    = SZI/SLEN
            ALFA   = SXI*A2XA(I) + SYI*A2YA(I) + SZI*A2ZA(I)
            ALFA   = MAX(ALFA,0.)
              
            SURFE  = 0.5*(CH(I) + CH(I-IL))
            FE(I)  = -SURFE*A2(I)*(ALFA*(TEMP(I)-TEMP(I-IL))/SLEN
     &             + DTDXI*(A2XA(I)-ALFA*SXI) + DTDYI*(A2YA(I)-ALFA*SYI)
     &             + DTDZI*(A2ZA(I)-ALFA*SZI))
c ... Thin-layer approximation can be used if the next line is activated
c            FE(I)  = -SURFE*A2(I)*(TEMP(I)-TEMP(I-IL))/SLEN
 850     CONTINUE
 900  CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLUXHS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXKT(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDIFF,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,XXTRAL,RKSI,INCHIML,ENTROPY_FIX,FRSSIE,PDI,
     7 PSEUCO,ARTSSP,NGL,FRSTEM,IUPPTL)

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE
      
      INTEGER :: IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,KMAXP1,KSTR,IDIR,
     2   ITURB,MAW,MAXW,ISTRID,JSTRID,IJTRID,IL,KL,KG,IA,IG,
     3   I,I1,ISTATE,mm,nn,ll,I8,NGL,IUPPTL

      REAL :: FRSDEN,FRSPRE,GAMMA,E0REF,T0REF,RMIPI,RNIPI,RWIPI,RMIMI,
     2   RNIMI,RWIMI,RKIPI,REIPI,RKIMI,REIMI,RKLIM,A11,A12,A13,A21,
     3   A22,A23,A31,A32,A33,PALPO,UR,VR,WR,UL,VL,WL,YPROR,YPROL,RKR,
     4   RKL,EIR,EIL,EPSLIM,PR,PKR,EPSR,HR,ER,ROHAT,YPROH,PL,PKL,EPSL,
     5   HL,EL,SQRL,CQL,CQR,UHAT,VHAT,WHAT,HHAT,EHAT,RKHAT,EPSHAT,PHAT,
     6   PMIN,PMAX,DRDHH,DRDPH,XTURB,CHAT2,CMIN,CMAX,PDPDE,DRO,DU,DKDP,
     7   PCHAT2,ALFA1,ALFA2,ALFA3,ALFA4,ALFA5,ALFA6,ALFA7,UC,RL1,RL2,
     8   RL3,RL4,RL5,RL6,RL7,CHAT,DK,DP,ROUCR,ROUCL,DAMP,F1PIG,F2PIG,
     9   F3PIG,F4PIG,F5PIG,F6PIG,F7PIG,F2PM,F3PM,F4PM,PRLAM,PRT,RGAS,
     1   RL1H,RL1R,RL1L,RL2H,RL2R,RL2L,RL3H,RL3R,RL3L,EPS22,SR,FRSSIE,
     2   PDR,PDL,PSEUCO,ARTSSP,PAVE,FRSTEM

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDIFF(*),D2(*),RKSI(*),PDI(*)

      LOGICAL :: XXTRAL, INCHIML, INCHIM, ENTROPY_FIX
					
      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:),REIM(:),REIP(:),ERL(:), ELR(:),DEDPH(:),
     +     DEDRH(:),TEIM(:),TEIP(:),PDMI(:),PDPI(:)
      REAL, TARGET ::  ZZZ(MAXW)

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      RKIM => ZZZ(12*MAW+1:13*MAW);RKIP => ZZZ(13*MAW+1:14*MAW) 
      REIM => ZZZ(14*MAW+1:15*MAW);REIP => ZZZ(15*MAW+1:16*MAW)
      ERL  => ZZZ(16*MAW+1:17*MAW); ELR => ZZZ(17*MAW+1:18*MAW)
      DEDRH=> ZZZ(18*MAW+1:19*MAW);DEDPH=> ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(20*MAW+1:21*MAW);TEIP => ZZZ(21*MAW+1:22*MAW)
      PDMI => ZZZ(22*MAW+1:23*MAW);PDPI => ZZZ(23*MAW+1:24*MAW)


C ... FLUXES FOR THREE-DIMENSIONAL FLOW

      INCHIM  = INCHIML 

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL
       
      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES

      IF(.NOT. XXTRAL) THEN
       
      CALL XXTRAP(TEMP, TEIP,TEIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   U, RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   V, RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   W, RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   P,  PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP( PDI, PDPI,PDMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP( RK, RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP(REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      ENDIF

      ELSE

      I8 =IJTRID
      CALL XXTRAA(TEMP,TEIP,TEIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   U,RMIP,RMIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   V,RNIP,RNIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   W,RWIP,RWIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   P, PPI, PMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA( PDI,PDPI,PDMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAA(  RK,RKIP,RKIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(REPS,REIP,REIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      ENDIF

      ENDIF ! .NOT.XXTRAL

      DO  500 IG  = 1,IJTRID
         I        = IA + IG
         IF(ISTATE /= 10) THEN ! Limit the pressure
         PPI(IG)  = MAX(0.001*FRSPRE,PPI(IG)) 
         PMI(IG)  = MAX(0.001*FRSPRE,PMI(IG))
         PDPI(IG) = MAX(-0.999*FRSPRE,PDPI(IG))
         PDMI(IG) = MAX(-0.999*FRSPRE,PDMI(IG))
      ENDIF
         RKIP(IG) = MAX(0.,RKIP(IG))
         RKIM(IG) = MAX(0.,RKIM(IG))
         REIP(IG) = MAX(0.,REIP(IG)) !RE = epsilon
         REIM(IG) = MAX(0.,REIM(IG))
 500  CONTINUE

C ... DENSITY AND SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE
C     AND TEMPERATURE. ISTATE Determines the Equation of State.

      CALL ROFPT(TEIP,PPI,ROIP,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL ROFPT(TEIM,PMI,ROIM,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL EFPT(EIP,PPI,TEIP,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)
      CALL EFPT(EIM,PMI,TEIM,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)

C ... DERIVATIVES ACCORDING TO GLAISTER

      I1     = IA + 1
      CALL DERIV(RO(I1),P(I1),DEDPH,DEDRH,DRDP(I1),DRDH(I1),
     + EIP,EIM,ELR,ERL,PPI,PMI,ROIP,ROIM,IJTRID,IL)


      DO  600 IG= 1,IJTRID
         I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)
						
C ... LIMITATIONS OF TURBULENCE QUANTITIES

         RKIPI    = MAX(0.001*RKLIM,RKIP(IG))
         REIPI    = MAX(0.001*EPSLIM,REIP(IG))
         RKIMI    = MAX(0.001*RKLIM,RKIM(IG))
	 REIMI    = MAX(0.001*EPSLIM,REIM(IG))

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         UR       = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VR       = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WR       = A31*RMIPI + A32*RNIPI + A33*RWIPI
         UL       = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VL       = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WL       = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... ROE'S FLUX SPLITTING                                         *
C*******************************************************************

      YPROR   = 1./ROIP(IG)
      RKR     = RKIPI*YPROR

      YPROL   = 1./ROIM(IG)
      RKL     = RKIMI*YPROL

C ... TOTAL INTERNAL ENERGY

      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL

      PDR     = PDPI(IG)		    !PD = pressure difference ?
      PR      = PPI(IG)
      PKR     = .6667*RKIPI
      EPSR    = REIPI*YPROR
      HR      = EIR + PR*YPROR
      ER      = EIP(IG)
      ROHAT   = SQRT(ROIP(IG)*ROIM(IG))
      YPROH   = 1./ROHAT

      PDL     = PDMI(IG)		    !PD = pressure difference ?
      PL      = PMI(IG)
      PKL     = .6667*RKIMI
      EPSL    = REIMI*YPROL
      HL      = EIL + YPROL*PL
      EL      = EIM(IG)
      SQRL    = 1./(ROIM(IG) + ROHAT)
      CQL     = ROIM(IG)*SQRL
      CQR     = 1.- CQL

      UHAT    = CQR*UR  + CQL*UL
      VHAT    = CQR*VR  + CQL*VL
      WHAT    = CQR*WR  + CQL*WL
      HHAT    = CQR*HR  + CQL*HL
      EHAT    = CQR*ER  + CQL*EL
      RKHAT   = CQR*RKR + CQL*RKL
      PHAT    = ROHAT*(HHAT - EHAT-.5*(UHAT**2+VHAT**2+WHAT**2)-RKHAT)

      PMIN    = MIN(PR,PL)
      PMAX    = MAX(PR,PL)
      PHAT    = MIN(PMAX,PHAT)
      PHAT    = MAX(PMIN,PHAT)

      DRDHH   = 1./(DEDRH(IG)-PHAT*YPROH**2)
      DRDPH   = -DRDHH*(DEDPH(IG) + YPROH)
C *********************************************************

      XTURB   = 1. + .6667*RKHAT*YPROH*DRDPH
      CHAT2   = XTURB/(DRDHH*YPROH + DRDPH)

      CMIN    = MIN(C(I),C(I-IL))
      CMAX    = MAX(C(I),C(I-IL))
      CHAT2   = MIN(CMAX**2,CHAT2)
      CHAT2   = MAX(CMIN**2,CHAT2)    ! EVALUOITAVA
      CHAT    = SQRT(CHAT2)

C ... CHARACTERISTIC SPEEDS

       UC      = UROT(I)

       RL1H    = UHAT-UC + CHAT
       RL2H    = UHAT-UC - CHAT
       RL3H    = UHAT-UC

       RL1R    = UR-UC + C(I)
       RL2R    = UR-UC - C(I)
       RL3R    = UR-UC

       RL1L    = UL-UC + C(I-IL)
       RL2L    = UL-UC - C(I-IL)
       RL3L    = UL-UC

C ... Try an entropy fix

      RL1     = ABS(RL1H)
      RL2     = ABS(RL2H)
      RL3     = ABS(RL3H)
      SR      = MAX(RL1,RL2,RL3)
c      SR      = MAX(ABS(RL1R),ABS(RL2R),ABS(RL1L),ABS(RL2L))

      IF(ENTROPY_FIX) THEN

c         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1H-RL1R))
         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1R-RL1H))
         IF(RL1 <= EPS22) RL1 = .5*(RL1H**2+EPS22**2)/EPS22

c         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2H-RL2R))
         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2R-RL2H))
         IF(RL2 <= EPS22) RL2 = .5*(RL2H**2+EPS22**2)/EPS22

c         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3H-RL3R))
         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3R-RL3H))
         IF(RL3 <= EPS22) RL3 = .5*(RL3H**2+EPS22**2)/EPS22

      ENDIF ! ENTROPY_FIX

      RL4     = RL3
      RL5     = RL3
      RL6     = RL3
      RL7     = RL3

C ... CONVECTIVE FLUXES

c      CHAT    = MAX(PSEUCO*ABS(UHAT),ARTSSP)

c      chat = artssp   ! Incompressible, watch your step!! Just joo
c      SR      = CHAT + ABS(UHAT) ! Incompressible ?

      SR = MAX(ABS(UR)+C(I),ABS(UL)+C(I-IL))

      ROUCR   = ROIP(IG)*(UR - UC)
      ROUCL   = ROIM(IG)*(UL - UC)
c      DAMP    = SR*(ROIP(IG) - ROIM(IG)) ! Alternative for a compressible flow
      UHAT = .5*(UR+UL)
      DAMP    = SR*(PDR - PDL)/CHAT**2
      F1PIG   = .5*A(I)*(ROUCR + ROUCL - DAMP)
      F2PIG   = .5*A(I)*(ROUCR*UR+PDR+PKR  +  ROUCL*UL+PDL+PKL  -
c      F2PIG   = .5*A(I)*(ROUCR*UR+PDR +  ROUCL*UL+PDL  - ! Auttaa massataseessa
     +          SR*(ROIP(IG)*UR - ROIM(IG)*UL))
      F3PIG   = .5*A(I)*(ROUCR*VR  +  ROUCL*VL - 
     +             SR*(ROIP(IG)*VR - ROIM(IG)*VL))
c     +             ABS(UHAT)*(ROIP(IG)*VR - ROIM(IG)*VL))
      F4PIG   = .5*A(I)*(ROUCR*WR  +  ROUCL*WL - 
     +            SR*(ROIP(IG)*WR - ROIM(IG)*WL))
c     +            ABS(UHAT)*(ROIP(IG)*WR - ROIM(IG)*WL))
      F5PIG   = .5*A(I)*(ROUCR*(HR+PKR*YPROR)  +  ROUCL*(HL+PKL*YPROL)
     +          + UC*(PR+PKR + PL+PKL) -
     +            SR*(ROIP(IG)*EIR  - ROIM(IG)*EIL))
      F6PIG   = .5*A(I)*(ROUCR*RKR  + ROUCL*RKL - 
     +           ABS(UHAT)*(ROIP(IG)*RKR  - ROIM(IG)*RKL))
      F7PIG   = .5*A(I)*(ROUCR*EPSR + ROUCL*EPSL - 
     +           ABS(UHAT)*(ROIP(IG)*EPSR - ROIM(IG)*EPSL))  


C*    HATS ARE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.
      HAT1(I) = DAMP
      HAT2(I) = ROUCR
      HAT3(I) = ROUCL
      HAT4(I) = RL3

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES


      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG + FE(I)
      FRK(I)  = F6PIG + FRK(I)
      FEPS(I) = F7PIG + FEPS(I)

600   CONTINUE

 1000 CONTINUE

      RETURN
      END SUBROUTINE FLUXKT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C *************************** K-EPSILON ********************************
      SUBROUTINE FLUXKE(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDIFF,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,XXTRAL,RKSI,INCHIML,ENTROPY_FIX,TRANSL,TRM)

      USE NS3CO, ONLY : IN, JN, KN
      USE TYPE_ARRAYS
         
      IMPLICIT NONE
         
      INTEGER :: IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,KMAXP1,KSTR,IDIR,
     2   ITURB,MAW,MAXW,ISTRID,JSTRID,IJTRID,IL,KL,KG,IA,IG,
     3   I,I1,ISTATE,mm,nn,ll,I8

      REAL :: FRSDEN,FRSPRE,GAMMA,E0REF,T0REF,RMIPI,RNIPI,RWIPI,RMIMI,
     2   RNIMI,RWIMI,RKIPI,REIPI,RKIMI,REIMI,RKLIM,A11,A12,A13,A21,
     3   A22,A23,A31,A32,A33,PALPO,UR,VR,WR,UL,VL,WL,YPROR,YPROL,RKR,
     4   RKL,EIR,EIL,EPSLIM,PR,PKR,EPSR,HR,ER,ROHAT,YPROH,PL,PKL,EPSL,
     5   HL,EL,SQRL,CQL,CQR,UHAT,VHAT,WHAT,HHAT,EHAT,RKHAT,EPSHAT,PHAT,
     6   PMIN,PMAX,DRDHH,DRDPH,XTURB,CHAT2,CMIN,CMAX,PDPDE,DRO,DU,DKDP,
     7   PCHAT2,ALFA1,ALFA2,ALFA3,ALFA4,ALFA5,ALFA6,ALFA7,UC,RL1,RL2,
     8   RL3,RL4,RL5,RL6,RL7,CHAT,DK,DP,ROUCR,ROUCL,DAMP,F1PIG,F2PIG,
     9   F3PIG,F4PIG,F5PIG,F6PIG,F7PIG,F2PM,F3PM,F4PM,PRLAM,PRT,RGAS,
     1   RL1H,RL1R,RL1L,RL2H,RL2R,RL2L,RL3H,RL3R,RL3L,EPS22,F9PIG,
     2   F10PIG,GR,GL,GHAT,RER,REL,REHAT,ALFA8,ALFA9,RL8,RL9

      REAL :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDIFF(*),D2(*),RKSI(*)
         
      LOGICAL :: XXTRAL, INCHIML, INCHIM, ENTROPY_FIX, TRANSL
         
      TYPE(INTERMITTENCY) :: TRM(*)

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:),REIM(:),REIP(:),ERL(:), ELR(:),DEDPH(:),
     +     DEDRH(:),TEIM(:),TEIP(:),G1R(:),G1L(:),RET1R(:),RET1L(:)
      REAL, TARGET ::  ZZZ(MAXW)

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      RKIM => ZZZ(12*MAW+1:13*MAW);RKIP => ZZZ(13*MAW+1:14*MAW) 
      REIM => ZZZ(14*MAW+1:15*MAW);REIP => ZZZ(15*MAW+1:16*MAW)
      ERL  => ZZZ(16*MAW+1:17*MAW); ELR => ZZZ(17*MAW+1:18*MAW)
      DEDRH=> ZZZ(18*MAW+1:19*MAW);DEDPH=> ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(20*MAW+1:21*MAW);TEIP => ZZZ(21*MAW+1:22*MAW) 
      G1R  => ZZZ(22*MAW+1:23*MAW);G1L  => ZZZ(23*MAW+1:24*MAW) 
      RET1R=> ZZZ(24*MAW+1:25*MAW);RET1L=> ZZZ(25*MAW+1:26*MAW) 

C ... FLUXES FOR THREE-DIMENSIONAL FLOW

      INCHIM  = INCHIML 

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES

      IF(.NOT. XXTRAL) THEN

      CALL XXTRAP( RO, ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  U, RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  V, RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  W, RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  P,  PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP( RK, RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP(REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      ENDIF
      IF(TRANSL) THEN
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
         ENDIF

      ELSE

      I8 =IJTRID
      CALL XXTRAA(  RO,ROIP,ROIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   U,RMIP,RMIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   V,RNIP,RNIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   W,RWIP,RWIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   P, PPI, PMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAA(  RK,RKIP,RKIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(REPS,REIP,REIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      ENDIF
      IF(TRANSL) THEN ! This is not distance-based. Who will write TXTRAA?
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
         ENDIF

      ENDIF ! .NOT.XXTRAL

      DO  500 IG = 1,IJTRID
         I       = IA + IG
         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
500   CONTINUE

C ... SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE AND DENSITY
C     ISTATE Deterines the Equation of State to be  Used.

      CALL EFPRO(EIP,PPI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(EIM,PMI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ELR,PPI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ERL,PMI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
							
C ... DERIVATIVES ACCORDING TO GLAISTER

      I1     = IA + 1
      CALL DERIV(RO(I1),P(I1),DEDPH,DEDRH,DRDP(I1),DRDH(I1),
     + EIP,EIM,ELR,ERL,PPI,PMI,ROIP,ROIM,IJTRID,IL)


      DO  600 IG= 1,IJTRID
         I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... LIMITATIONS OF TURBULENCE QUANTITIES

         RKIPI    = MAX(0.001*RKLIM,RKIP(IG))
         REIPI    = MAX(0.001*EPSLIM,REIP(IG))
         RKIMI    = MAX(0.001*RKLIM,RKIM(IG))
	 REIMI    = MAX(0.001*EPSLIM,REIM(IG))

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         UR       = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VR       = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WR       = A31*RMIPI + A32*RNIPI + A33*RWIPI
         UL       = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VL       = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WL       = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... ROE'S FLUX SPLITTING                                         *
C*******************************************************************

      YPROR   = 1./ROIP(IG)
      RKR     = RKIPI*YPROR

      YPROL   = 1./ROIM(IG)
      RKL     = RKIMI*YPROL

C ... TOTAL INTERNAL ENERGY

      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL

      PR      = PPI(IG)
      PKR     = .6667*RKIPI
      EPSR    = REIPI*YPROR
      HR      = EIR + PR*YPROR
      ER      = EIP(IG)
      ROHAT   = SQRT(ROIP(IG)*ROIM(IG))
      YPROH   = 1./ROHAT

      PL      = PMI(IG)
      PKL     = .6667*RKIMI
      EPSL    = REIMI*YPROL
      HL      = EIL + YPROL*PL
      EL      = EIM(IG)
      SQRL    = 1./(ROIM(IG) + ROHAT)
      CQL     = ROIM(IG)*SQRL
      CQR     = 1.- CQL

      UHAT    = CQR*UR  + CQL*UL
      VHAT    = CQR*VR  + CQL*VL
      WHAT    = CQR*WR  + CQL*WL
      HHAT    = CQR*HR  + CQL*HL
      EHAT    = CQR*ER  + CQL*EL
      RKHAT   = CQR*RKR + CQL*RKL
      EPSHAT  = CQR*EPSR+ CQL*EPSL
      PHAT    = ROHAT*(HHAT - EHAT-.5*(UHAT**2+VHAT**2+WHAT**2)-RKHAT)

      PMIN    = MIN(PR,PL)
      PMAX    = MAX(PR,PL)
      PHAT    = MIN(PMAX,PHAT)
      PHAT    = MAX(PMIN,PHAT)

      DRDHH   = 1./(DEDRH(IG)-PHAT*YPROH**2)
      DRDPH   = -DRDHH*(DEDPH(IG) + YPROH)

      IF(TRANSL) THEN
         GR   = G1R(IG);   GL  = G1L(IG);   GHAT  = CQR*GR  + CQL*GL
         RER  = RET1R(IG); REL = RET1L(IG); REHAT = CQR*RER + CQL*REL
      ENDIF

C *********************************************************

      XTURB   = 1. + .6667*RKHAT*YPROH*DRDPH
      CHAT2   = XTURB/(DRDHH*YPROH + DRDPH)

      CMIN    = MIN(C(I),C(I-IL))
      CMAX    = MAX(C(I),C(I-IL))
      CHAT2   = MIN(CMAX**2,CHAT2)
      CHAT2   = MAX(CMIN**2,CHAT2)    ! EVALUOITAVA
      CHAT    = SQRT(CHAT2)
      PDPDE   = DEDPH(IG)

C ... MULTIPLIERS (CHARACTERISTIC VARIABLES)

      DRO     = ROIP(IG)-ROIM(IG)
      DU      = UR - UL
      DK      = RKR - RKL
      DP      = PR - PL + .6667*(RKHAT*DRO + ROHAT*DK)
      PCHAT2  = 1./CHAT**2
      ALFA1   = .5*(DP + ROHAT*CHAT*DU)*PCHAT2
      ALFA2   = .5*(DP - ROHAT*CHAT*DU)*PCHAT2
c      ALFA11   = .5*( DP)*PCHAT2
c      ALFA22   = .5*( DP)*PCHAT2
      ALFA3   = DRO - DP*PCHAT2
      ALFA4   = ROHAT*(VR - VL)
      ALFA5   = ROHAT*(WR - WL)
      ALFA6   = ROHAT*DK
      ALFA7   = ROHAT*(EPSR - EPSL)
      IF(TRANSL) THEN
      ALFA8   = ROHAT*(GR - GL)
      ALFA9   = ROHAT*(RER - REL)
      ENDIF

C ... ABSOLUTE VALUES OF THE CHARACTERISTIC SPEEDS

c      UC      = UROT(I)
c      RL1     = ABS(UHAT-UC + CHAT)
c      RL2     = ABS(UHAT-UC - CHAT)
c      RL3     = ABS(UHAT-UC)

       UC      = UROT(I)

       RL1H    = UHAT-UC + CHAT
       RL2H    = UHAT-UC - CHAT
       RL3H    = UHAT-UC

       RL1R    = UR-UC + C(I)
       RL2R    = UR-UC - C(I)
       RL3R    = UR-UC

       RL1L    = UL-UC + C(I-IL)
       RL2L    = UL-UC - C(I-IL)
       RL3L    = UL-UC

C ... Try an entropy fix

      RL1     = ABS(RL1H)
      RL2     = ABS(RL2H)
      RL3     = ABS(RL3H)

      IF(ENTROPY_FIX) THEN

c         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1H-RL1R))
         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1R-RL1H))
         IF(RL1 <= EPS22) RL1 = .5*(RL1H**2+EPS22**2)/EPS22
        
c         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2H-RL2R))
         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2R-RL2H))
         IF(RL2 <= EPS22) RL2 = .5*(RL2H**2+EPS22**2)/EPS22

c         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3H-RL3R))
         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3R-RL3H))
         IF(RL3 <= EPS22) RL3 = .5*(RL3H**2+EPS22**2)/EPS22

      ENDIF ! ENTROPY_FIX

      RL4     = RL3
      RL5     = RL3
      RL6     = RL3
      RL7     = RL3
      RL8     = RL3
      RL9     = RL3

C ... CONVECTIVE FLUXES

      ROUCR   = ROIP(IG)*(UR - UC)
      ROUCL   = ROIM(IG)*(UL - UC)
      DAMP    = RL1*ALFA1 + RL2*ALFA2 + RL3*ALFA3
      F1PIG   = .5*A(I)*(ROUCR + ROUCL - DAMP)
      F2PIG   = .5*A(I)*(ROUCR*UR+PR+PKR  +  ROUCL*UL+PL+PKL          -
     +          (DAMP*UHAT + (RL1*ALFA1-RL2*ALFA2) * CHAT))
c       if(icycle == 4) write(789,*)'f2pig,roucr,roucl,uhat,rl1,rl2',
c     +f2pig,roucr,roucl,uhat,rl1,rl2,alfa1,alfa2,chat,RL1*ALFA1-RL2*ALFA2
c     +   (DAMP*UHAT+(RL1*ALFA11-RL2*ALFA22) * CHAT+rohat*abs(uhat)*du))
      F3PIG   = .5*A(I)*(ROUCR*VR  +  ROUCL*VL - DAMP*VHAT - RL4*ALFA4)
      F4PIG   = .5*A(I)*(ROUCR*WR  +  ROUCL*WL - DAMP*WHAT - RL5*ALFA5)
      F5PIG   = .5*A(I)*(ROUCR*(HR+PKR*YPROR)  +  ROUCL*(HL+PKL*YPROL)
     +          + UC*(PR+PKR + PL+PKL)
     +          - DAMP*HHAT - ((RL1*ALFA1 - RL2*ALFA2)*UHAT*CHAT
     +          - RL3*ALFA3*ROHAT*CHAT**2*PDPDE + RL4*ALFA4*VHAT
     +          + RL5*ALFA5*WHAT) + RL6*ALFA6*(1.-.6667*ROHAT*PDPDE))
      F6PIG   = .5*A(I)*(ROUCR*RKR + ROUCL*RKL - DAMP*RKHAT - RL6*ALFA6)
      F7PIG   = .5*A(I)*(ROUCR*EPSR+ROUCL*EPSL - DAMP*EPSHAT- RL7*ALFA7)

C	For this index I use left(L) values	
      IF(TRANSL) THEN
         F9PIG  = .5*A(I)*(ROUCR*GR + ROUCL*GL - DAMP*GHAT - RL8*ALFA8)
         F10PIG = .5*A(I)*(ROUCR*RER+ ROUCL*REL- DAMP*REHAT- RL9*ALFA9) 
      ENDIF

C*    HATS ARE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.
      HAT1(I) = DAMP
      HAT2(I) = ROUCR
      HAT3(I) = ROUCL
      HAT4(I) = RL3

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES


      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG + FE(I)
      FRK(I)  = F6PIG + FRK(I)
      FEPS(I) = F7PIG + FEPS(I)

C     Update transition model fluxes with the massflow part
      IF(TRANSL) THEN
         TRM(I)%FG   = TRM(I)%FG   + F9PIG
         TRM(I)%FRET = TRM(I)%FRET + F10PIG
      ENDIF

600   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE FLUXKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C *************************** CUSP ********************************
      SUBROUTINE FLUCUSP(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDIFF,A2XA,A2YA,A2ZA,D2,PRLAM,PRT,
     4 IN,JN,KN,UROT,TEMP,CH,DRDH,DRDP,ISTATE,RGAS,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,RKLIM,EPSLIM,HAT1,HAT2,HAT3,HAT4,KMAXP1,KSTR,IDIR,
     6 ITURB,ZZZ,MAW,MAXW,XXTRAL,RKSI,INCHIML,ENTROPY_FIX,TRANSL,TRM)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,KMAXP1,KSTR,IDIR,
     2   ITURB,MAW,MAXW,ISTRID,JSTRID,IN,JN,KN,IJTRID,IL,KL,KG,IA,IG,
     3   I,I1,ISTATE,mm,nn,ll,I8

      REAL FRSDEN,FRSPRE,GAMMA,E0REF,T0REF,RMIPI,RNIPI,RWIPI,RMIMI,
     2   RNIMI,RWIMI,RKIPI,REIPI,RKIMI,REIMI,RKLIM,A11,A12,A13,A21,
     3   A22,A23,A31,A32,A33,PALPO,UR,VR,WR,UL,VL,WL,YPROR,YPROL,RKR,
     4   RKL,EIR,EIL,EPSLIM,PR,PKR,EPSR,HR,ER,ROHAT,YPROH,PL,PKL,EPSL,
     5   HL,EL,SQRL,CQL,CQR,UHAT,VHAT,WHAT,HHAT,EHAT,RKHAT,EPSHAT,PHAT,
     6   PMIN,PMAX,DRDHH,DRDPH,XTURB,CHAT2,CMIN,CMAX,PDPDE,DRO,DU,DKDP,
     7   PCHAT2,ALFA1,ALFA2,ALFA3,ALFA4,ALFA5,ALFA6,ALFA7,UC,RL1,RL2,
     8   RL3,RL4,RL5,RL6,RL7,CHAT,DK,DP,ROUCR,ROUCL,DAMP,F1PIG,F2PIG,
     9   F3PIG,F4PIG,F5PIG,F6PIG,F7PIG,F2PM,F3PM,F4PM,PRLAM,PRT,RGAS,
     1   RL1H,RL1R,RL1L,RL2H,RL2R,RL2L,RL3H,RL3R,RL3L,EPS22,RMA,ALFA,
     2   BETA,F8PIG,F9PIG,GR,GL,GHAT,RER,REL,REHAT,ALFA8,ALFA9,RL8,RL9

      REAL A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*),
     5   HAT2(*),HAT3(*),HAT4(*),PDIFF(*),D2(*),RKSI(*)

      LOGICAL XXTRAL,INCHIML,INCHIM,ENTROPY_FIX

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:),REIM(:),REIP(:),ERL(:), ELR(:),DEDPH(:),
     +     DEDRH(:),TEIM(:),TEIP(:),G1R(:),G1L(:),RET1R(:),RET1L(:)
      REAL, TARGET ::  ZZZ(MAXW)

      TYPE(INTERMITTENCY) :: TRM(*)

      LOGICAL TRANSL

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      RKIM => ZZZ(12*MAW+1:13*MAW);RKIP => ZZZ(13*MAW+1:14*MAW) 
      REIM => ZZZ(14*MAW+1:15*MAW);REIP => ZZZ(15*MAW+1:16*MAW)
      ERL  => ZZZ(16*MAW+1:17*MAW); ELR => ZZZ(17*MAW+1:18*MAW)
      DEDRH=> ZZZ(18*MAW+1:19*MAW);DEDPH=> ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(20*MAW+1:21*MAW);TEIP => ZZZ(21*MAW+1:22*MAW) 
      G1R  => ZZZ(22*MAW+1:23*MAW);G1L  => ZZZ(23*MAW+1:24*MAW) 
      RET1R=> ZZZ(24*MAW+1:25*MAW);RET1L=> ZZZ(25*MAW+1:26*MAW) 


C ... FLUXES FOR THREE-DIMENSIONAL FLOW USING A H-CUSP SCHEME

      INCHIM  = INCHIML 

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL
       
      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES

      IF(.NOT. XXTRAL) THEN

      CALL XXTRAP( RO, ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  U, RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  V, RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  W, RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  P,  PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP( RK, RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP(REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      IF(TRANSL) THEN
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
      ENDIF

      ENDIF

      ELSE

      I8 =IJTRID
      CALL XXTRAA(  RO,ROIP,ROIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   U,RMIP,RMIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   V,RNIP,RNIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   W,RWIP,RWIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(   P, PPI, PMI,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAA(  RK,RKIP,RKIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      CALL XXTRAA(REPS,REIP,REIM,I8,IL,KL,IA,INTEM,A,VOL,D2,RKSI,INCHIM)
      ENDIF
      IF(TRANSL) THEN ! This is not distance-based. Who will write TXTRAA?
       CALL TXTRAP(TRM,G1R,G1L,RET1R,RET1L,IJTRID,IL,KL,IA,INTET,RKSI,
     +  INCHIML)
      ENDIF

      ENDIF ! .NOT.XXTRAL

      DO  500 IG = 1,IJTRID
         I       = IA + IG
         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
500   CONTINUE

C ... SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE AND DENSITY
C     ISTATE Deterines the Equation of State to be  Used.

      CALL EFPRO(EIP,PPI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(EIM,PMI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ELR,PPI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ERL,PMI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)

C ... DERIVATIVES ACCORDING TO GLAISTER

      I1     = IA + 1
      CALL DERIV(RO(I1),P(I1),DEDPH,DEDRH,DRDP(I1),DRDH(I1),
     + EIP,EIM,ELR,ERL,PPI,PMI,ROIP,ROIM,IJTRID,IL)


      DO  600 IG= 1,IJTRID
         I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... LIMITATIONS OF TURBULENCE QUANTITIES

         RKIPI    = MAX(0.001*RKLIM,RKIP(IG))
         REIPI    = MAX(0.001*EPSLIM,REIP(IG))
         RKIMI    = MAX(0.001*RKLIM,RKIM(IG))
	 REIMI    = MAX(0.001*EPSLIM,REIM(IG))

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         UR       = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VR       = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WR       = A31*RMIPI + A32*RNIPI + A33*RWIPI
         UL       = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VL       = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WL       = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... ROE'S FLUX SPLITTING                                         *
C*******************************************************************

      YPROR   = 1./ROIP(IG)
      RKR     = RKIPI*YPROR

      YPROL   = 1./ROIM(IG)
      RKL     = RKIMI*YPROL

C ... TOTAL INTERNAL ENERGY

      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL

      PR      = PPI(IG)
      PKR     = .6667*RKIPI
      EPSR    = REIPI*YPROR
      HR      = EIR + PR*YPROR
      ER      = EIP(IG)
      ROHAT   = SQRT(ROIP(IG)*ROIM(IG))
      YPROH   = 1./ROHAT

      PL      = PMI(IG)
      PKL     = .6667*RKIMI
      EPSL    = REIMI*YPROL
      HL      = EIL + YPROL*PL
      EL      = EIM(IG)
      SQRL    = 1./(ROIM(IG) + ROHAT)
      CQL     = ROIM(IG)*SQRL
      CQR     = 1.- CQL

      UHAT    = CQR*UR  + CQL*UL
      VHAT    = CQR*VR  + CQL*VL
      WHAT    = CQR*WR  + CQL*WL
      HHAT    = CQR*HR  + CQL*HL
      EHAT    = CQR*ER  + CQL*EL
      RKHAT   = CQR*RKR + CQL*RKL
      EPSHAT  = CQR*EPSR+ CQL*EPSL
      PHAT    = ROHAT*(HHAT - EHAT-.5*(UHAT**2+VHAT**2+WHAT**2)-RKHAT)

      PMIN    = MIN(PR,PL)
      PMAX    = MAX(PR,PL)
      PHAT    = MIN(PMAX,PHAT)
      PHAT    = MAX(PMIN,PHAT)

      DRDHH   = 1./(DEDRH(IG)-PHAT*YPROH**2)
      DRDPH   = -DRDHH*(DEDPH(IG) + YPROH)

      IF(TRANSL) THEN
         GR   = G1R(IG);   GL  = G1L(IG);   GHAT  = CQR*GR  + CQL*GL
         RER  = RET1R(IG); REL = RET1L(IG); REHAT = CQR*RER + CQL*REL
      ENDIF

C *********************************************************

      XTURB   = 1. + .6667*RKHAT*YPROH*DRDPH
      CHAT2   = XTURB/(DRDHH*YPROH + DRDPH)

      CMIN    = MIN(C(I),C(I-IL))
      CMAX    = MAX(C(I),C(I-IL))
      CHAT2   = MIN(CMAX**2,CHAT2)
      CHAT2   = MAX(CMIN**2,CHAT2)    ! EVALUOITAVA
      CHAT    = SQRT(CHAT2)
      RMA     = UHAT/CHAT
      PDPDE   = DEDPH(IG)

C ... MULTIPLIERS (CHARACTERISTIC VARIABLES)

      DRO     = ROIP(IG)-ROIM(IG)
      DU      = UR - UL
      DK      = RKR - RKL
      DP      = PR - PL + .6667*(RKHAT*DRO + ROHAT*DK)
      PCHAT2  = 1./CHAT**2
      ALFA1   = .5*(DP + ROHAT*CHAT*DU)*PCHAT2
      ALFA2   = .5*(DP - ROHAT*CHAT*DU)*PCHAT2
c      ALFA11   = .5*( DP)*PCHAT2
c      ALFA22   = .5*( DP)*PCHAT2
      ALFA3   = DRO - DP*PCHAT2
      ALFA4   = ROHAT*(VR - VL)
      ALFA5   = ROHAT*(WR - WL)
      ALFA6   = ROHAT*DK
      ALFA7   = ROHAT*(EPSR - EPSL)
      IF(TRANSL) THEN
      ALFA8   = ROHAT*(GR - GL)
      ALFA9   = ROHAT*(RER - REL)
      ENDIF

C ... ABSOLUTE VALUES OF THE CHARACTERISTIC SPEEDS

c      UC      = UROT(I)
c      RL1     = ABS(UHAT-UC + CHAT)
c      RL2     = ABS(UHAT-UC - CHAT)
c      RL3     = ABS(UHAT-UC)

       UC      = UROT(I)

       RL1H    = UHAT-UC + CHAT
       RL2H    = UHAT-UC - CHAT
       RL3H    = UHAT-UC

       RL1R    = UR-UC + C(I)
       RL2R    = UR-UC - C(I)
       RL3R    = UR-UC

       RL1L    = UL-UC + C(I-IL)
       RL2L    = UL-UC - C(I-IL)
       RL3L    = UL-UC

C ... Try an entropy fix

      RL1     = ABS(RL1H)
      RL2     = ABS(RL2H)
      RL3     = ABS(RL3H)

      IF(ENTROPY_FIX) THEN

c         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1H-RL1R))
         EPS22   = 4.*MAX(1.E-10,(RL1H-RL1L),(RL1R-RL1H))
         IF(RL1 <= EPS22) RL1 = .5*(RL1H**2+EPS22**2)/EPS22

c         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2H-RL2R))
         EPS22   = 4.*MAX(1.E-10,(RL2H-RL2L),(RL2R-RL2H))
         IF(RL2 <= EPS22) RL2 = .5*(RL2H**2+EPS22**2)/EPS22

c         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3H-RL3R))
         EPS22   = 4.*MAX(1.E-10,(RL3H-RL3L),(RL3R-RL3H))
         IF(RL3 <= EPS22) RL3 = .5*(RL3H**2+EPS22**2)/EPS22

      ENDIF ! ENTROPY_FIX

      RL4     = RL3
      RL5     = RL3
      RL6     = RL3
      RL7     = RL3

C ... CONVECTIVE FLUXES

      ROUCR   = ROIP(IG)*(UR - UC)
      ROUCL   = ROIM(IG)*(UL - UC)


      IF(ABS(RMA) > .2) THEN
         ALFA    = ABS(RMA)
      ELSE
         ALFA    = (RMA**2 + .04)/.4
      ENDIF

      IF(RMA >= 0. .AND. RMA <= 1.) THEN
         BETA = MAX(0.,2.*RMA-1.)
      ELSEIF(RMA > -1. .AND. RMA < 0.) THEN
         BETA = MIN(0.,2*RMA+1.)
      ELSE
         BETA = SIGN(1.,RMA)
      ENDIF

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)-ROIM(IG)) +
     +           BETA*(ROUCR - ROUCL)
c      write(7999,*)I,IG,IDIR,RMA,ALFA,BETA,ALFA*CHAT-BETA*UHAT,UHAT,CHAT
      F1PIG   = .5*A(I)*(ROUCR + ROUCL - DAMP)

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*UR - ROIM(IG)*UL) +
     +           BETA*(ROUCR*UR+PR+PKR  -  ROUCL*UL-PL-PKL)
      F2PIG   = .5*A(I)*(ROUCR*UR+PR+PKR  +  ROUCL*UL+PL+PKL - DAMP)
c       if(icycle == 4) write(789,*)'f2pig,roucr,roucl,uhat,rl1,rl2',
c     +f2pig,roucr,roucl,uhat,rl1,rl2,alfa1,alfa2,chat,RL1*ALFA1-RL2*ALFA2
c     +   (DAMP*UHAT+(RL1*ALFA11-RL2*ALFA22) * CHAT+rohat*abs(uhat)*du))

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*VR - ROIM(IG)*VL) +
     +           BETA*(ROUCR*VR  -  ROUCL*VL)
      F3PIG   = .5*A(I)*(ROUCR*VR  +  ROUCL*VL - DAMP)

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*WR - ROIM(IG)*WL) +
     +           BETA*(ROUCR*WR - ROUCL*WL)
      F4PIG   = .5*A(I)*(ROUCR*WR  +  ROUCL*WL - DAMP)

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*(HR+PKR*YPROR) -
     +           ROIM(IG)*(HL+PKL*YPROL)) +
     +           BETA*(ROUCR*(HR+PKR*YPROR) - ROUCL*(HL+PKL*YPROL))
      F5PIG   = .5*A(I)*(ROUCR*(HR+PKR*YPROR)  +  ROUCL*(HL+PKL*YPROL)
     +          + UC*(PR+PKR + PL+PKL)         - DAMP)

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*RKR - ROIM(IG)*RKL) +
     +           BETA*(ROUCR*RKR - ROUCL*RKL)
      F6PIG   = .5*A(I)*(ROUCR*RKR + ROUCL*RKL - DAMP)

      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*EPSR - ROIM(IG)*EPSL)+
     +           BETA*(ROUCR*EPSR - ROUCL*EPSL)
      F7PIG   = .5*A(I)*(ROUCR*EPSR + ROUCL*EPSL - DAMP)
      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*GR - ROIM(IG)*GL)+
     +           BETA*(ROUCR*GR - ROUCL*GL)
      F8PIG   = .5*A(I)*(ROUCR*GR + ROUCL*GL - DAMP)
      DAMP    = (ALFA*CHAT - BETA*UHAT)*(ROIP(IG)*RER - ROIM(IG)*REL)+
     +           BETA*(ROUCR*RER - ROUCL*REL)
      F9PIG   = .5*A(I)*(ROUCR*RER + ROUCL*REL - DAMP)
c       write(7998,*) idir,ig,f1pig,f2pig,f3pig,f5pig,f6pig
C*    HATS ARE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.
      HAT1(I) = DAMP
      HAT2(I) = ROUCR
      HAT3(I) = ROUCL
      HAT4(I) = RL3

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES


      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
      FRO(I)  = F1PIG + FRO(I)
      FRM(I)  = F2PM  + FRM(I)
      FRN(I)  = F3PM  + FRN(I)
      FRW(I)  = F4PM  + FRW(I)
      FE(I)   = F5PIG + FE(I)
      FRK(I)  = F6PIG + FRK(I)
      FEPS(I) = F7PIG + FEPS(I)

C     Update transition model fluxes with the massflow part
      IF(TRANSL) THEN
         TRM(I)%FG   = TRM(I)%FG   + F8PIG
         TRM(I)%FRET = TRM(I)%FRET + F9PIG
      ENDIF


600   CONTINUE

1000  CONTINUE
      RETURN
c       write(7999,*) '*********************************************'
      END SUBROUTINE FLUCUSP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUDKE(FE,FRK,FEPS,VOL,A,RO,RK,REPS,IMAX,JMAX,KMAX,
     2     IDI,VIS,EPS2,VIST,FUN1,psigk1,psige1,psigk2,psige2,
     3     KSTR,IDIR,IDIFF,TRANSL,TRM)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,IDI,KSTR,IDIR,IDIFF,ISTRID,JSTRID,KSTRID,
     2     ILL,IA,IL,IG,JG,KG,I,IMINID,IMAXID,JMINID,JMAXID,KMINID,
     3     KMAXID,IA2
      REAL MU,MUT,GPM,RETPM,CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1,EPS,
     2     ONEMF1,PSIGK,PSIGK1,PSIGK2,PSIGE,PSIGE1,PSIGE2,YPD2,SURG,
     3     RMUT2,SURFK,SURFEP,YPRO1,YPRO2,F6VIG,F7VIG
      REAL :: A(*),RO(*),FE(*),VIS(*),EPS2(*),VIST(*),VOL(*),
     3     RK(*),REPS(*),FRK(*),FEPS(*),FUN1(*)
      LOGICAL :: INVIS, TRANSL
      TYPE(INTERMITTENCY)    :: TRM(*)

C
C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (TWO-EQUATION TURB.)
C
      EPS     = 1.E-10

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN
      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +            IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

      IF(TRANSL) THEN
         CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1) 
      ENDIF

      IF(IDI /= 0) THEN
C
C ... THIN-LAYER VISCOUS FLUXES ARE CALCULATED (LAMINAR AND TURBULENT)
C
      DO 900 KG = 1,KMAXID
      IA      = (KG+KN-1)*ILL
      DO 850 JG = 1,JMAXID
      IA2       = (JG+JN-1)*ISTRID + IA + IN
      DO 850 IG = 1,IMAXID
            I      = IG + IA2

C     Menter's blending of coefficients. In case of (ITURB /= 6)
C     FUN1 = 0. always and PSIG?1 = 0. and PSIG?2 = PSIG? 
            ONEMF1 = MAX(0.0,(1.0 - FUN1(I)))
            PSIGK  = FUN1(I)*PSIGK1 + ONEMF1*PSIGK2
            PSIGE  = FUN1(I)*PSIGE1 + ONEMF1*PSIGE2

            YPD2   = A(I)/(VOL(I) + VOL(I-IL))
            SURG   = -A(I)*YPD2
            RMUT2  = (EPS2(I)-1.)*VIS(I) + (EPS2(I-IL)-1.)*VIS(I-IL)
c           RMUT2  = VIST(I) + VIST(I-IL) ! vois vaihtaa tahan joskus
            SURFK  = SURG*(VIS(I)+VIS(I-IL)  + RMUT2*PSIGK)
            SURFEP = SURG*(VIS(I)+VIS(I-IL)  + RMUT2*PSIGE)
            
            YPRO1  = 1./RO(I)
            YPRO2  = 1./RO(I-IL)
            
            F6VIG  = SURFK *(  RK(I)*YPRO1  - RK(I-IL)*YPRO2)
            F7VIG  = SURFEP*(REPS(I)*YPRO1- REPS(I-IL)*YPRO2)

            FRK(I) = F6VIG
            FEPS(I)= F7VIG
            FE(I)  = F6VIG ! Added 30.12.2003

C          Intermittency variables. (Moved from FLUXDO 25.6.2015)

           IF(TRANSL) THEN 
	       MU	   = VIS(I) + VIS(I-IL)
             MUT     = (EPS2(I)-1.)*VIS(I) + 
     +                 (EPS2(I-IL)-1.)*VIS(I-IL)
	       GPM	   = TRM(I)%G   - TRM(I-IL)%G
	       RETPM  	   = TRM(I)%RET - TRM(I-IL)%RET
             TRM(I)%FG   = SURG*(MU + MUT/SIGF)*GPM
             TRM(I)%FRET = SURG*(MU + MUT)*SIGOT*RETPM
           ENDIF

 850     CONTINUE
 900  CONTINUE
      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED

      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I      = IA2 + IG
         FRK(I) = 0.
         FEPS(I)= 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLUDKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLSCSA(RO,F1RNU,REPS,VOL,A,
     2 PSIGSC,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 NGL,FRSDEN,ITURB,KSTR,IDIR)

C ... Diffusion thin-layer flux for RNUT in Spalart-Allmaras

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,NGL,ITURB,
     2   KSTR,IDIR,ISTRID,JSTRID,KSTRID,ILL,IL,IA,IA2,IG,JG,KG,I,IMINID,
     3   IMAXID,JMINID,JMAXID,KMINID,KMAXID

      REAL :: PSIGSC, FRSDEN, YPD2, SURG, SURFI1, PSIGE1, F1VIG
      REAL :: RO(*), F1RNU(*), REPS(*), VOL(*), A(*), VIS(*)

      LOGICAL :: INVIS

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID

      IL      = KSTR
      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

c      write(*,'(8I4)') nbl,IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID

      PSIGE1  = 1./PSIGSC

C ... SCHMIDT'S NUMBER  = PSIGE = 1/PSIGE1
C
C ... VISCOUS FLUXES ARE ADDED

      DO 900 KG = 1,KMAXID
         IA      = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = IA + (JG+JN-1)*ISTRID + IN
         DO 850 IG = 1,IMAXID
            I      = IA2 + IG
            YPD2    = A(I)/(VOL(I) + VOL(I-IL))
            SURG    = -A(I)*YPD2
            SURFI1  = SURG*PSIGE1*
     +          (VIS(I)+REPS(I) + VIS(I-IL)+REPS(I))
            F1VIG   = SURFI1*(REPS(I)/RO(I) - REPS(I-IL)/RO(I-IL))
            F1RNU(I)= F1VIG
 850     CONTINUE
 900  CONTINUE
      
      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I      = IA2 + IG
         F1RNU(I)  = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLSCSA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLKEBL(FE,FRK,FEPS,VOL,A,A2XA,A2YA,A2ZA,D2,
     1     RO,RK,REPS,IMAX,JMAX,KMAX,IDI,VIS,EPS2,VIST,FUN1,
     2     DKDX,DKDY,DKDZ,DEPSDX,DEPSDY,DEPSDZ,XC,YC,ZC,
     3     PSIGK1,PSIGE1,PSIGK2,PSIGE2,KSTR,IDIR,IDIFF)

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IDI,KSTR,IDIR,ISTRID,JSTRID,
     + KSTRID,ILL,IA,IL,IMINID,JMINID,KMINID,IMAXID,JMAXID,KMAXID,
     + IG,JG,KG,IA2,I,IDIFF

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: A(*),RO(*),FE(*),VIS(*),EPS2(*),VIST(*),VOL(*),D2(*),
     2 RK(*),REPS(*),FRK(*),FEPS(*),DKDX(*),DKDY(*),DKDZ(*),FUN1(*),
     3 A2XA(*),A2YA(*),A2ZA(*),DEPSDX(*),DEPSDY(*),DEPSDZ(*)

      REAL :: PSIGK1,PSIGE1,PSIGK2,PSIGE2,ALFA,SLEN,SXI,SYI,SZI,EPS,
     2 ONEMF1,PSIGK,PSIGE,W1,W2,DKDXI,DKDYI,DKDZI,DEDXI,DEDYI,DEDZI,
     3 RMUT2,SURFK,SURFEP,YPRO1,YPRO2,frki,frkj,suhde

      LOGICAL INVIS

      EPS     = 1.E-10

C     Viscous flux calculation using Blazek formulae for k-e/w variables

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN

      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

      IF(IDI /= 0) THEN
      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2

C     Menter's blending of coefficients. In case of (ITURB /= 6)
C     FUN1 = 0. always and PSIG?1 = 0. and PSIG?2 = PSIG? 
            
            ONEMF1 = MAX(0.0,(1.0 - FUN1(I)))
            PSIGK  = FUN1(I)*PSIGK1 + ONEMF1*PSIGK2
            PSIGE  = FUN1(I)*PSIGE1 + ONEMF1*PSIGE2

C ... Distance weighting, although it does not has any effect

            W1  = 2.0*( D2(I) / ( D2(I) + D2(I-IL) ) ) 
            W2  = 2.0 - W1

            DKDXI   = 0.5*(W2*DKDX(I)   + W1*DKDX(I-IL))
            DKDYI   = 0.5*(W2*DKDY(I)   + W1*DKDY(I-IL))
            DKDZI   = 0.5*(W2*DKDZ(I)   + W1*DKDZ(I-IL))
            DEDXI   = 0.5*(W2*DEPSDX(I) + W1*DEPSDX(I-IL))
            DEDYI   = 0.5*(W2*DEPSDY(I) + W1*DEPSDY(I-IL))
            DEDZI   = 0.5*(W2*DEPSDZ(I) + W1*DEPSDZ(I-IL))

c            DKDXI   = 0.5*(DKDX(I) + DKDX(I-IL))
c            DKDYI   = 0.5*(DKDY(I) + DKDY(I-IL))
c            DKDZI   = 0.5*(DKDZ(I) + DKDZ(I-IL))
c            DEDXI   = 0.5*(DEPSDX(I) + DEPSDX(I-IL))
c            DEDYI   = 0.5*(DEPSDY(I) + DEPSDY(I-IL))
c            DEDZI   = 0.5*(DEPSDZ(I) + DEPSDZ(I-IL))

            SXI    = XC(I) - XC(I-IL)
            SYI    = YC(I) - YC(I-IL)
            SZI    = ZC(I) - ZC(I-IL)
            SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
            SXI    = SXI/SLEN
            SYI    = SYI/SLEN
            SZI    = SZI/SLEN
            ALFA   = SXI*A2XA(I) + SYI*A2YA(I) + SZI*A2ZA(I)
            ALFA   = MAX(ALFA,0.)

            RMUT2  = VIST(I) + VIST(I-IL)

            SURFK  = 0.5*A(I)*(VIS(I)+VIS(I-IL) + RMUT2*PSIGK)
            SURFEP = 0.5*A(I)*(VIS(I)+VIS(I-IL) + RMUT2*PSIGE) 
            YPRO1  = 1./RO(I)
            YPRO2  = 1./RO(I-IL)  

            FRK(I) = 
     &     -SURFK*(ALFA* (RK(I)*YPRO1 - RK(I-IL)*YPRO2)  / SLEN +
     &      DKDXI*(A2XA(I)-ALFA*SXI) + DKDYI*(A2YA(I)-ALFA*SYI) +
     &      DKDZI*(A2ZA(I)-ALFA*SZI))

            FEPS(I)= 
     &     -SURFEP*(ALFA*(REPS(I)*YPRO1 - REPS(I-IL)*YPRO2)/SLEN+
     &      DEDXI*(A2XA(I)-ALFA*SXI) + DEDYI*(A2YA(I)-ALFA*SYI) +
     &      DEDZI*(A2ZA(I)-ALFA*SZI))   
            FE(I)  = FRK(I) ! Added 30.12.2003

 850     CONTINUE
 900  CONTINUE
      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
         FRK(I)  = 0.
         FEPS(I) = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLKEBL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLKENS(FE,FRK,FEPS,VOL,A,A2XA,A2YA,A2ZA,
     1     RO,RK,REPS,IMAX,JMAX,KMAX,IDI,VIS,EPS2,VIST,FUN1,
     2     DKDX,DKDY,DKDZ,DEPSDX,DEPSDY,DEPSDZ,
     3     psigk1,psige1,psigk2,psige2,KSTR,IDIR)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),FE(*),VIS(*),EPS2(*),VIST(*),VOL(*),
     2 RK(*),REPS(*),FRK(*),FEPS(*),DKDX(*),DKDY(*),DKDZ(*),FUN1(*),
     3 A2XA(*),A2YA(*),A2ZA(*),DEPSDX(*),DEPSDY(*),DEPSDZ(*)
      LOGICAL :: INVIS

      EPS     = 1.E-10
C
C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (TWO-EQUATION TURB.)
C ... Old unused version
C
      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN

      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

      IF(IDI /= 0) THEN
      DO 900 KG = 1,KMAXID
         IA        = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = (JG+JN-1)*ISTRID + IA + IN
         DO 850 IG = 1,IMAXID
            I      = IG + IA2

C     Menter's blending of coefficients. In case of (ITURB /= 6)
C     FUN1 = 0. always and PSIG?1 = 0. and PSIG?2 = PSIG? 
            
            ONEMF1 = MAX(0.0,(1.0 - FUN1(I)))
            PSIGK  = FUN1(I)*PSIGK1 + ONEMF1*PSIGK2
            PSIGE  = FUN1(I)*PSIGE1 + ONEMF1*PSIGE2

            DKDXI   = 0.5*(DKDX(I) + DKDX(I-IL))
            DKDYI   = 0.5*(DKDY(I) + DKDY(I-IL))
            DKDZI   = 0.5*(DKDZ(I) + DKDZ(I-IL))
            DEDXI   = 0.5*(DEPSDX(I) + DEPSDX(I-IL))
            DEDYI   = 0.5*(DEPSDY(I) + DEPSDY(I-IL))
            DEDZI   = 0.5*(DEPSDZ(I) + DEPSDZ(I-IL))

            RMUT2   = VIST(I) + VIST(I-IL)

            SURFK   = 0.5*A(I)*(VIS(I)+VIS(I-IL) + RMUT2*PSIGK)
            SURFEP  = 0.5*A(I)*(VIS(I)+VIS(I-IL) + RMUT2*PSIGE) 

            FRK(I)  =             
     &      -SURFK*(DKDXI*A2XA(I) + DKDYI*A2YA(I) + DKDZI*A2ZA(I))      
            FEPS(I) = 
     &      -SURFEP*(DEDXI*A2XA(I) + DEDYI*A2YA(I) + DEDZI*A2ZA(I))      
            FE(I)  = FRK(I) ! Added 30.12.2003

 850     CONTINUE
 900  CONTINUE
      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
         FRK(I)  = 0.
         FEPS(I) = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLKENS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLFUN(UTAUM,F2VWF,F3VWF,F4VWF,A3,A3X,A3Y,A3Z,VIS,
     2  VIST,EPS2,UWALL,VWALL,WWALL,TWALL,HFLUX,ISTRID,JSTRID,KSTRID,
     3  IDIM,IDIR,KBOT,KX1,KX2,KY1,KY2,ITURB,ISTR,JSTR,KSTR,IN,JN,KN,
     4  VOL,RO,RELAX,U,V,W,PTUR,REPS,EPSOLD,STRAIN,RBK,FUN1)

      USE NS3CO, ONLY : LN
      USE CONSTANTS, ONLY : EPS12, EPS

      IMPLICIT NONE

      INTEGER I,J,K,ISTRID,JSTRID,KSTRID,IDIM,IDIR,KBOT,KX1,KX2,KY1,KY2,
     2  IN,JN,KN,ISTR,JSTR,KSTR,KA,NFL,LSTRID,IJ,NN,NR,N1,II,II2,IM1,
     3  IM2,IF1,IF2,NI,NH,NP,ITER,ITERM,KOKO,ITURB

      REAL UTAUM(*),F2VWF(*),F3VWF(*),F4VWF(*),A3(*),A3X(*),A3Y(*),
     2  A3Z(*),VIS(*),VIST(*),EPS2(*),UWALL(*),VWALL(*),WWALL(*),
     3  TWALL(*),HFLUX(*),VOL(*),RO(*),U(*),V(*),W(*),PTUR(*),
     4  REPS(*),EPSOLD(*),STRAIN(*),RBK(*),FUN1(*)

      REAL YPLUS1,YPLUSR,YPLUSM,DF1DUT,DFRDUT,F1LOG,FRLOG,DUTAU1,DUTAUR,
     2  UPLUS1,UPLUSR,UPLUSM,CF1,CFR,PHIB1,DIST,CKAPPA,UHAT11,UNOR,UTX,
     3  UTY,UTZ,UTABS,RELAX,RELAR,ARG,BETA1,BSTAR,OMVIS,OMLOG,OMB1,OMB2,
     4  PHIOM,OMBC,PLUSFC,PSKP,SKPLUS,SKPMIN,SR1,SR2,SRCOEF,RNUB,UTAU2,
     5  YPLUSS,BETA2,BETA,RK,KPLUSM,DB,RL
     
      REAL, ALLOCATABLE :: UTAU1(:), UTAUR(:), U1(:), UHATX(:),
     2                     UHATY(:), UHATZ(:)

      ITERM   = 10 ! Inputtiin?
      CKAPPA  = 0.41; BETA1 = 0.075; BETA2 = 0.0828; BSTAR = 0.09
      RELAR   = RELAX !?
      
      KOKO    = (KY2-KY1+1)*(KX2-KX1+1)
      
      ALLOCATE(UTAU1(KOKO),UTAUR(KOKO),U1(KOKO),UHATX(KOKO),
     2         UHATY(KOKO),UHATZ(KOKO)) 

      LSTRID  = KSTR*IDIR          ! IDIR = 1 or -1 (ILMO in CALL WALLFUN)
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

      KA      = (KN+KBOT-1)*KSTR

      IF(IDIM > KBOT) THEN ! What is the meaning of this here and in BOTSKE?

      DO J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NR   = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      DO I = KX1,KX2
         K        = NN + I                 ! Local p-index withoud ghost cells
         II       = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         NP       = I  + NR                ! Patch index with LN ghost cells
         IF1      = II  + NFL              ! FLUX INDEX
         UTAUR(K) = UTAUM(NP)              ! Start from the existing values
         UTAU1(K) = UTAUM(NP) 	       ! -"-
         RK       = RBK(NP)	             ! Surface rougness

         UHAT11   = U(II)**2 + V(II)**2 + W(II)**2 + EPS12
         U1(K)    = SQRT(UHAT11)
         UNOR     = A3X(IF1)*U(II) + A3Y(IF1)*V(II) + A3Z(IF1)*W(II)
         UTX      = U(II) - UNOR*A3X(IF1)
         UTY      = V(II) - UNOR*A3Y(IF1)
         UTZ      = W(II) - UNOR*A3Z(IF1)
         UTABS    = SQRT(UTX**2 + UTY**2 + UTZ**2 + EPS12)
         UHATX(K) = UTX/UTABS; UHATY(K) = UTY/UTABS; UHATZ(K) =UTZ/UTABS

      ENDDO; ENDDO

      DO ITER = 1,ITERM

      DO J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO I = KX1,KX2
         K       = NN + I                 ! Local p-index without ghost cells
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! 1. ghost cell
         IM2     = IM1- LSTRID            ! 2. ghost cell
         IF1     = II  + NFL              ! FLUX INDEX
         IF2     = IF1 + LSTRID           ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP      = I  + NR                ! Patch index with LN ghost cells

         DIST    = VOL(II)/(A3(IF1)+A3(IF2))
         YPLUS1  = RO(II)*UTAU1(K)*DIST/VIS(II)
         YPLUSR  = RO(II)*UTAUR(K)*DIST/VIS(II)
         YPLUSM  = RO(II)*UTAUM(NP)*DIST/VIS(II)
         KPLUSM  = RO(II)*UTAUM(NP)*RK/VIS(II)

c         KPLUSM  = 70. ! Test
         
         DB	 = 1./CKAPPA*LOG(1. + 0.3*KPLUSM)
         RL	 = EXP(5.1-0.25*DB)/EXP(5.1)

         DF1DUT  = U1(K)/UTAU1(K)**2 + 1./(CKAPPA*UTAU1(K))
         DFRDUT  = U1(K)/UTAUR(K)**2 + YPLUSR/UTAUR(K) * (
     &   0.4/(CKAPPA*(1.+0.4*YPLUSR)) +7.8*(1./11.*EXP(-YPLUSR/11.) - 
     &   (1./11. + YPLUSR/33.)*EXP(-MIN(YPLUSR/3.,500.))) )

         F1LOG   = -U1(K)/UTAU1(K) + 1./CKAPPA*LOG(YPLUS1) + 5.1 -DB
         
c         FRLOG   = -U1(K)/UTAUR(K) + 1./CKAPPA*LOG(1+0.4*YPLUSR)+7.8*
c     &   (1.- EXP(-YPLUSR/11.) - YPLUSR/11.*EXP(-MIN(YPLUSR/3.,500.)))
         FRLOG   = -U1(K)/UTAUR(K) + 1./CKAPPA*LOG(1+0.4*YPLUSR)+7.8*
     &   (1.- EXP(-YPLUSR*RL/11.) - 
     &    YPLUSR*RL/11.*EXP(-MIN(YPLUSR*RL/3.,500.)))
     &   - DB
c         FRLOG	 = MAX(0.001,FRLOG)

         DUTAU1  = -F1LOG/DF1DUT
         DUTAUR  = -FRLOG/DFRDUT

         DUTAU1  = DUTAU1/(1.+ABS(DUTAU1)/(RELAX*UTAU1(K)))
         DUTAUR  = DUTAUR/(1.+ABS(DUTAUR)/(RELAR*UTAUR(K)))
         
         UTAU1(K)= UTAU1(K) + DUTAU1
         UTAUR(K)= UTAUR(K) + DUTAUR

         ARG     = YPLUSM/27. ; PHIB1 = TANH(ARG**4)
         UTAUM(NP) = (1.-PHIB1)*UTAUR(K) + PHIB1*UTAU1(K)

      ENDDO; ENDDO

      ENDDO  ! ITER
      
      DO J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO I = KX1,KX2
         K       = NN + I                 ! Local p-index without ghost cells
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! 1. ghost cell
         IM2     = IM1- LSTRID            ! 2. ghost cell
         IF1     = II  + NFL              ! FLUX INDEX
         IF2     = IF1 + LSTRID           ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP      = I  + NR                ! Patch index with LN ghost cells
         
         F2VWF(K) = RO(II)*UTAUM(NP)**2*UHATX(K)
         F3VWF(K) = RO(II)*UTAUM(NP)**2*UHATY(K)
         F4VWF(K) = RO(II)*UTAUM(NP)**2*UHATZ(K)
         
         DIST   = VOL(II)/(A3(IF1)+A3(IF2))
         BETA   = FUN1(II)*BETA1 + (1.- FUN1(II))*BETA2
         
C     Surface value of omega:

         RNUB    = VIS(II)/RO(II)
         UTAU2   = RNUB*STRAIN(II)
         PLUSFC  = SQRT(UTAU2)/RNUB
	    
C     This is the new varibale k+ smooth wall boundary condition: 

         YPLUS1  = 2.0*PLUSFC*DIST
         SKPMIN  = MIN(2.4*YPLUS1**0.85 , 8.0) 
         SKPMIN  = MAX(SKPMIN, EPS)
         SKPLUS  = MAX(RBK(NP)*PLUSFC , SKPMIN) ! pinnan karheus.
c         SKPLUS  = SKPMIN

C     The variable k+-method ends here.

         PSKP    = 1.0/SKPLUS
         SR1     = 2500.0*PSKP**2
         SR2     = 100.0 *PSKP
         SRCOEF  = MAX(SR1,SR2)
         OMVIS   = STRAIN(II)*SRCOEF
         OMVIS  = MIN(OMVIS,6.*VIS(II)/(BETA*RO(II)*DIST**2))  ! 6.
c         OMVIS  = 6.*VIS(II)/(BETA*RO(II)*DIST**2)  ! 6.
c                  write(668,*) OMVIS,STRAIN(II),U1(K)/DIST 

         OMLOG  = UTAUM(NP)/(SQRT(BSTAR)*CKAPPA*DIST)
         OMB1   = OMVIS + OMLOG
         OMB2   = (OMVIS**1.2 + OMLOG**1.2)**(1./1.2)
         PHIOM  = TANH((0.1*RO(II)*UTAUM(NP)*DIST/VIS(II))**4)
         OMBC   = PHIOM*OMB1 + (1.- PHIOM)*OMB2
         
         REPS(II) = RO(II)*OMBC; EPSOLD(II) = REPS(II)    ! B.C for omega
         PTUR(II) = RO(II)*UTAUM(NP)**2*U1(K)/DIST        ! Pelaako?
         
         yplus1   = RO(II)*UTAUM(NP)*DIST/VIS(II)
         uplus1   = U1(K)/UTAUM(NP)
         ypluss   = UPLUS1 + EXP(-5.22*CKAPPA) * (EXP(CKAPPA*UPLUS1)-1.
     &   -CKAPPA*UPLUS1-.5*(CKAPPA*UPLUS1)**2-.1667*(CKAPPA*UPLUS1)**3)

c         WRITE(668,*) j,NP,F2VWF(K),OMVIS,OMLOG,PHIOM   
c         WRITE(668,*) j,NP,F2VWF(K),OMB1,OMB2,OMBC*RO(II),REPS(II) 
c         IF(K == 30) write(667,*) J,YPLUS1,YPLUSS,UPLUS1
         
      ENDDO; ENDDO

      ENDIF  ! IF(IDIR > KBOT)
      
      DEALLOCATE(UTAU1,UTAUR,U1,UHATX,UHATY,UHATZ)

      END SUBROUTINE WALLFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOTSKE(F3R,F3RM,F3RN,F3RW,F3E,P,U,V,W,C,E,VOL,A3,A3X,
     2  A3Y,A3Z,OHMI,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     3  PDI,PRO,VAR,KX1,KX2,KY1,KY2,CP,PR,PRT,RO,TOM,IFLUX,
     4  CX,CY,CZ,CMX,CMY,CMZ,M,IN,JN,KN,XC,YC,ZC,RKLIM,XMOM,YMOM,ZMOM,
     5  IDIR,DX,DY,DZ,F3RK,F3EPS,RK,REPS,PTUR,UROT,TEMP,CH,QX,QY,QZ,
     6  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,
     7  ISTR,JSTR,KSTR,IHEAT,RCON,IDIS,PSIGK1,PSIGE1,PSIGK2,PSIGE2,
     9  UWALL,VWALL,WWALL,CPWALL,TWALL,QWALL,QWFRIC,HFLUX,TAUW1,TAUW2,
     +  TAUWX,TAUWY,TAUWZ,SURLE,DSURLE,WMFLU2,PORO2,WHSTA2,WTEM2,RBK,
     1  RSDIRX,RSDIRY,RSDIRZ,SURFX,SURFY,SURFZ,XCP,YCP,ZCP,
     2  MAX11,CMU,XCO,YCO,ZCO,INTEM,INTET,MULPHL,MULPHC,N,
     3  TURCOR,ICYCLE,FRESUL,F1H,IWAVEB,BUX,BUY,BUZ,VMX,VMY,VMZ,STRAIN,
     4  TRANSL,TRM,WAVE,SURFPX,SURFPY,SURFPZ,PRC,UTAUM,
     5  F2VWF,F3VWF,F4VWF,TURLIM)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY: PII,EPS10
      USE NS3CO, ONLY : LN, GX, GY, GZ, GROUND, ALTITUDE, G0,
     2  WALLFUNL

      IMPLICIT NONE

      INTEGER :: IDIR,IDIS,ITURB,MAX11,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     2  KX1,KX2,KY1,KY2,IFLUX,M,IN,JN,KN,ISTR,JSTR,KSTR,NFL,IHFL,IHT,
     3  KA,I,J,NN,NR,NI,N1,II,IJ,II2,IF1,IF2,NH,NP,IZZ,IF3,II3,IA,N,
     4  IHEAT,LSTRID,ISTATE,INTEM,INTET,IM1,IM2,ICYCLE,IP,K

      INTEGER :: IWAVEB(*)

      REAL :: F3R(*),F3RW(*),F3E(*),P(*),V(*),U(*),C(*),RO(*),A3(*),
     2  A3Y(*),A3Z(*),A3X(*),OHMI(*),VIS(*),EPS2(*),VIST(*),PDI(*),
     3  F3RN(*),W(*),VOL(*),F3RM(*),CP(*),
     4  RK(*),REPS(*),F3RK(*),F3EPS(*),UROT(*),
     5  TEMP(*),CH(*),RCON(*),
     6  PTUR(*),E(*),
     7  UWALL(*),VWALL(*),WWALL(*),CPWALL(*),TWALL(*),QWALL(*),
     8  QWFRIC(*),HFLUX(*),TAUW1(*),TAUW2(*),SURFX(*),SURFY(*),
     9  SURFZ(*),TAUWX(*),TAUWY(*),TAUWZ(*),
c    9  SURFZ(*),XCP(*),YCP(*),ZCP(*),TAUWX(*),TAUWY(*),TAUWZ(*),
     +  SURLE(*),DSURLE(*),WMFLU2(*),PORO2(*),WHSTA2(*),WTEM2(*),
     1  RBK(*),RSDIRX(*),RSDIRY(*),RSDIRZ(*),F1H(*),STRAIN(*),WAVE(*),
     2  SURFPX(*),SURFPY(*),SURFPZ(*),UTAUM(*),F2VWF(*),F3VWF(*),
     3  F4VWF(*)

      REAL :: PR,PRT,TOM,CX,CY,CZ,CMX,CMY,CMZ,RKLIM, ! XMOM,YMOM,ZMOM,
     2  DX,DY,DZ,QX,QY,QZ,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     3  PSIGK1,PSIGE1,PSIGK2,PSIGE2,CMU,WTEMP,WMFLUX,POROS,EPS3,
     4  RC1,RC2,RM1,RM2,PSIGK,PSIGE,ROMFAC,REPFAC,PRS,YMR,YPG,YP1,
     5  YP4,RSUM,EPS,TSURF,VISLIM,CHLIM,ROLIM,XAUX1,XAUX2,XAUX3,
     6  STARI,STAR2,DIST,SURG,EPS2B,VISTB,VISB,CHSB,STRAIB,OHMIB,
     7  ROB,EPSUR2,SURF,SURFE,YPRO1,YPRO2,RKSURF,SQRK,EPSUR3,EPSURF,
     8  RNUB,UTAU2,PLUSFC,YPLUS1,SKPMIN,SKPLUS,PSKP,SR1,SR2,SRCOEF,
     9  OMSURF,USURF,VSURF,WSURF,TEMP1,TEMP2,CVI2,UVI2,VVI2,WVI2,
     +  UVEL,YP,PARA,ALPHA,CVII,UVII,WVII,DUC,F2V,F3V,F4V,F5VF,F5VQ,
     1  F5V,F6V,F7V,EVII,UTAU,EVIK,RKII,ANXI,ANYI,ANZI,XNORI,XNORJ,
     2  ANXJ,ANYJ,ANZJ,STRAI,STRA2,CVI1,UVI1,VVI1,WVI1,VVII,XKAPPA,
     3  PB1,PB,PB4,YPD2,F2VIG,F3VIG,F4VIG,F5VIG,XD,YD,ZD,BX,BY,BZ,
     4  DELTAP,X,Y,Z,F2VS,F3VS,F4VS,F5VS,F6VS,F7VS,
     5  EB,TEMPB,EINSUC,USUCT,VSUCT,WSUCT,TSUCT,WHSTAG,SUCKIN,HB,
     6  CVS2,UVS2,VVS2,WVS2,F5VFS,F5VQS,WALLMF,OMPOR,RAPUB,SURFRK,
     7  SURFEP,DISTG,DISTF,SURFE1,SURGE,CHS,BUX,BUY,BUZ,GTOT,BXDP,
     8  BYDP,BZDP,VMX,VMY,VMZ,APUL,APUG,DREF,DMAX,DTREF,DTSUB,DW,
     9  DTSUP,RN,RNREF,F,GW,GWA,HFG,A2E,WAVEB,VIS1B,VIS2B,CHS1B,
     +  CHS2B,CV1(2),CV2(2),UV1(2),UV2(2),VV1(2),VV2(2),WV1(2),WV2(2),
     +  DXP,DYP,DZP,EPS2X,EPS2Y,EPS2Z,TURLIM

      REAL :: XMULTI(6)

      REAL :: F2VMF(2),F3VMF(2),F4VMF(2),F5VMF(2),F5VFMF(2),F5VQMF(2)

      REAL :: XCO(*), YCO(*), ZCO(*), XC(*), YC(*), ZC(*),
     2        XCP(*), YCP(*), ZCP(*)

      REAL ::   XMOM,YMOM,ZMOM

      LOGICAL :: MULPHL, TURCOR, FRESUL, TRANSL

      CHARACTER(LEN = *) :: MULPHC

      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(INTERMITTENCY) TRM(*)
      TYPE(PRE_COR) PRC(*)
C
C ... LU FACTORED METHOD.BOUNDARY CONDITION FOR THE SOLID BOUNDARY
C ... BOUNDARY FLUXES ARE CALCULATED EXPLICITLY. POINTS ARE NUMBERED
C ... XI-WISE
C
C     Modified in order to work also for the k-omega model. Also 
C     permeable-surface fluxes included.
C                                            AH 11.2.1997 & PR 9.3.99
      DATA XMULTI /1.,2.,2.,1.,2.,1./

      EPS3    = 1.E-20

C ... Injection or suction

       POROS = 0.0
       APUL  = 0.0 ; APUG = 0.0

C ... ONE MINUS POROSITY

       OMPOR = 1.- POROS 

C ... Testing Testing
c      OMPOR = 0.1 
c      OMPOR = 10.

      XKAPPA  = 0.4
      RC1     = 1.500
      RC2     = RC1 - 1.
      IF(INTEM == 6) THEN
         RC1  = 1.
         RC2  = RC1 - 1.
      ENDIF
C ... FOR TURBULENCE
      RM1     = 1.000
      RM2     = RM1 - 1.
C     k-epsilon or SST

      IF (IDIS == 2) THEN        ! k-omega model
         RM1    = RC1
         RM2    = RC2
         PSIGK  = PSIGK1
         PSIGE  = PSIGE1
         ROMFAC = 1.0
         REPFAC = 1.0 - ROMFAC
      ELSE                       ! k-epsilon model
C        RM1    = 1.500          ! Activate this one day
C        RM2    = RM1 - 1.0      ! Activate this one day
         RM1    = 1.000
         RM2    = RM1 - 1.0
         PSIGK  = PSIGK2
         PSIGE  = PSIGE2
         ROMFAC = 0.0
         REPFAC = 1.0 - ROMFAC
      END IF

C ... TURBULENT VISCOCITY IS ZERO FOR RSM.
      IF(ITURB > 23) THEN
         RM1 = 0.
         RM2 = 0.
      ENDIF
      IF(ITURB >= 10 .AND. ITURB <= 19) THEN
         RM1 = 0.
         RM2 = 0.
      ENDIF
C ... this is not exatly correct, but may work better in lower levels.

      PRS     = PR/PRT

      LSTRID  = KSTR*IDIR          ! IDIR = 1 or -1 (ILMO in CALL BOTSKE)
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

      YMR     = .166666667
      YPG     = DIFPRE
      YP1     = 1.
      YP4     = 0.
      RSUM    = 0.
      IF(M == 1) RSUM = 1.

C ... INCO, TWO-PHASE and Kurganov-Tadmor fluxes

      IF(IFLUX == 4 .OR. IFLUX == 6 .OR. IFLUX == 7) THEN ! PRESSURE DIFFERENCE
         YPG  = DIFPRE - FRSPRE
         YP1  = 0.
         YP4  = 1.
      ENDIF

      EPS     = 1.E-3                            !  20.12.1994 : 10
      TSURF   = T0
      IHFL    = 0
      IHT     = 1

      IF(IHEAT == 2 .OR. IHEAT >= 12) IHFL  = 1  ! Heat flux is given
      IF(IHEAT == 0)  IHT  = 0
      IF(IHEAT == 15) IHT  = 0
      IF(IHEAT == 15) IHFL = 0
       
C ... IHEAT = 0 ADIABATIC
C ...         1 WALL TEMPERATURE IS GIVEN IN TWALL
C ...         2 HEAT FLUX IS GIVEN        IN HFLUX
C ...         3 WALL TEMPERATURE IS PUT TO FREE_STREAM STAGNATION TEMEPERATURE

      VISLIM  = 0.01*FRSVIS
      CHLIM  = 100.0*FRSVIS
      ROLIM   = .001*FRSDEN

C ... EPSILON TILDE IS ZERO AT SURFACE IN CHIENS EPSILON EQ. PPR 12.94

      XAUX1 = 0.0  !diffusion of the turbulence variables
      XAUX2 = 0.0  !interpolate dissipation value at the wall
      XAUX3 = 0.0  !Speziale's dissipation value at the wall
      IF(ITURB ==  3 .OR. ITURB ==  4 .OR. ITURB == 7 .OR. 
     +   ITURB == 11 .OR. ITURB >= 21) THEN
         XAUX3 = 0.
      ELSEIF(ITURB == 10) THEN
         XAUX1 = 1.             !diffusion of the turbulence variables
         XAUX3 = 1.
      ENDIF
      IF(IDIS == 2) XAUX1 = 1.

C ... Start with viscous surface fluxes and then add the inviscid part

      KA      = (KN+KBOT-1)*KSTR

      IF(IDIM > KBOT) THEN

      DO 1000 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1000 I = KX1,KX2
         K       = NN + I                 ! Local p-index without ghost cells
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! 1. ghost cell
         IM2     = IM1- LSTRID            ! 2. ghost cell
         IF1 = II  + NFL                  ! FLUX INDEX
         IF2 = IF1 + LSTRID               ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP  = I   + NR                   ! Patch index with LN ghost cells

         STRA2   = STRAIN(II2)
         STRAI   = STRAIN(II)
         DIST   = (A3(IF1)+A3(IF1+LSTRID))/VOL(II)
         DISTG  = (A3(IF1)+A3(IF1-LSTRID))/VOL(IM1)
         SURG    =-DIST*YMR
C        EPS2B   = MAX(1.    ,RM1*EPS2(II) - RM2*EPS2(II2))
         EPS2B   = 1.           ! work better
C        VISTB   = MAX(0.    ,RM1*VIST(II) - RM2*VIST(II2))
         VISTB   = 0.           ! work better
         VISB    = MAX(VISLIM,RC1* VIS(II) - RC2* VIS(II2))
         CHSB    = MAX(EPS   ,RC1*  CH(II) - RC2*  CH(II2))
         STRAIB  = MAX(EPS   ,RC1*STRAI    - RC2*STRA2)
         OHMIB   = MAX(EPS   ,RC1*OHMI(II) - RC2*OHMI(II2))
         PB1     = MAX(0.01*FRSPRE,RC1*P(II)-RC2*P(II2))
C        ROB     = MAX(EPS   ,RC1*  RO(II) - RC2*  RO(II2))
         ROB     = MAX(ROLIM ,RC1*  RO(II) - RC2*  RO(II2))
         EPSUR2  = MAX(EPS   ,RC1*REPS(II) - RC2*REPS(II2))
c        TEMPB   = MAX(EPS   ,RC1*TEMP(II) - RC2*TEMP(II2))
         TEMPB   = MAX(0.01*T0, RC1*TEMP(II) - RC2*TEMP(II2))
         SURF    = SURG*EPS2B*VISB
         SURFE   = SURG*(PRS*VISTB/VISB+1.)*CHSB
         RAPUB   = MAX(1.0E-20,EPS2B-1.0)        ! TTB 27.3.2000
         SURFRK  = SURG*(1.0 + PSIGK*RAPUB)*VISB ! TTB 15.3.2000
         SURFEP  = SURG*(1.0 + PSIGE*RAPUB)*VISB ! TTB 15.3.2000
         
         IF(MULPHL .AND. MULPHC == 'MULTI') THEN
            VIS1B  = MAX(0.,RC1*PRO(II)%VIS(1) - RC2*PRO(II2)%VIS(1))
            VIS2B  = MAX(0.,RC1*PRO(II)%VIS(2) - RC2*PRO(II2)%VIS(2))
            CHS1B  = MAX(0.,RC1*PRO(II)%CH(1)  - RC2*PRO(II2)%CH(1))
            CHS2B  = MAX(0.,RC1*PRO(II)%CH(2)  - RC2*PRO(II2)%CH(2))
         ENDIF
         
         SURFE1  = .5*DIST*(PRS*VISTB/VISB+1.)*CHSB
         YPRO1   = 1./RO(II)
         YPRO2   = 1./RO(II2)

         CHS     = MAX(EPS   ,RC1* CH(IM1) - RC2*  CH(IM2))
         SURGE  = 0.5*DISTG*CHS

         RKSURF = 0.0
         TSUCT  = 0.0
C ... in Speziales model:
c         SQRK    = IDIR*(9.*SQRT(RK(II)*YPRO1)-1.*
c     +        SQRT(RK(II2)*YPRO2)-8.*SQRT(RKSURF/ROB))*DIST/6.
C ... in Speziales model the first-order approximation here is
C ... somewhat more stable. Both ways produce a kink near the wall
         SQRK    = (SQRT(RK(II)*YPRO1)-SQRT(RKSURF/ROB))*DIST
         EPSUR3  = 2.*VISB*SQRK**2
         EPSURF  = XAUX2*EPSUR2/ROB + XAUX3*EPSUR3/ROB
         
C     Surface value of omega:

         RNUB    = VISB/ROB
         UTAU2   = RNUB*STRAIB
         PLUSFC  = SQRT(UTAU2)/RNUB
	    
C     This is the new varibale k+ smooth wall boundary condition: 

         YPLUS1  = 2.0*PLUSFC/DIST
         SKPMIN  = MIN(2.4*YPLUS1**0.85 , 8.0) 
         SKPMIN  = MAX(SKPMIN, EPS)
         SKPLUS  = MAX(RBK(NP)*PLUSFC , SKPMIN) ! pinnan karheus.
c         SKPLUS  = SKPMIN

C     The variable k+-method ends here.

         PSKP    = 1.0/SKPLUS
         SR1     = 2500.0*PSKP**2
         SR2     = 100.0 *PSKP
         SRCOEF  = MAX(SR1,SR2)
         OMSURF  = STRAIB*SRCOEF

C     Surface value of k and omega/epsilon

         RKSURF  = 0.0

         EPSURF  = ROMFAC*OMSURF + REPFAC*EPSURF

         USURF   = UWALL(NP)
         VSURF   = VWALL(NP)
         WSURF   = WWALL(NP)

         TEMP1   = TEMP(II)
         TEMP2   = TEMP(II+LSTRID)

         IF(IHEAT == 1) THEN
            TSURF     = TWALL(NP)
         ELSE IF(IHEAT /= 3 .AND. IHEAT /= 4) THEN
            TSURF     = T0
            TWALL(NP) = TEMPB
         ELSE IF(IHEAT == 3) THEN
            TSURF     = TEMPB
            TWALL(NP) = TEMPB
         ELSE IF(IHEAT == 4) THEN
            TSURF     = (SURFE1*TEMP(II)+SURGE*TEMP(IM1))/(SURFE1+SURGE)
            TWALL(NP) = TSURF
         ENDIF

C ... Second- and first-order fluxes for the wall
         
         CVI2    = IHT*IDIR*(9.000*TEMP1 - 1.000*TEMP2  - 8.*TSURF)
         UVI2    = IDIR*(9.000*U(II) - 1.000*U(II2) - 8.*USURF)
         VVI2    = IDIR*(9.000*V(II) - 1.000*V(II2) - 8.*VSURF)
         WVI2    = IDIR*(9.000*W(II) - 1.000*W(II2) - 8.*WSURF)
         CVI1    = IHT*IDIR*(TEMP1 - TSURF)
         UVI1    = IDIR*(U(II) - USURF)*6.
         VVI1    = IDIR*(V(II) - VSURF)*6.
         WVI1    = IDIR*(W(II) - WSURF)*6.

C ... INVISCID FLUXES:  EXTRAPOLATION OF THE PRESSURE

C        PB4      = RC1*PDI(II) - RC2*PDI(II2)
C ...  kaatuu. Eika kaadu
         PB4      = MAX(-.99*FRSPRE, RC1*PDI(II) - RC2*PDI(II2))
         PB       = YP1*PB1 + YP4*PB4
c         pb = .5*(pb+cpwall(np)) ! relaksaatiokokeilu
         CPWALL(NP)= PB4

         IF(MULPHL .AND. MULPHC == 'MULTI') THEN

         PRC(IF1)%FPM = A3(IF1)*A3X(IF1)*PB4
         PRC(IF1)%FPN = A3(IF1)*A3Y(IF1)*PB4
         PRC(IF1)%FPW = A3(IF1)*A3Z(IF1)*PB4

         DO IP = 1,NPHASES
         TEMP1  = PRO(II)%TEMP(IP)
         TEMP2  = PRO(II+LSTRID)%TEMP(IP)
         CV2(IP)= IHT*IDIR*(9.000*TEMP1 - 1.000*TEMP2  - 8.*TSURF)
         UV2(IP)= IDIR*(9.000*VAR(II)%U(IP)-1.0*VAR(II2)%U(IP)-8.*USURF)
         VV2(IP)= IDIR*(9.000*VAR(II)%V(IP)-1.0*VAR(II2)%V(IP)-8.*VSURF)
         WV2(IP)= IDIR*(9.000*VAR(II)%W(IP)-1.0*VAR(II2)%W(IP)-8.*WSURF)
         CV1(IP)= IHT*IDIR*(TEMP1 - TSURF)
         UV1(IP)= IDIR*(VAR(II)%U(IP) - USURF)*6.
         VV1(IP)= IDIR*(VAR(II)%V(IP) - VSURF)*6.
         WV1(IP)= IDIR*(VAR(II)%W(IP) - WSURF)*6.
         ENDDO
         ENDIF

         IF(IHEAT == 4) THEN    ! Coupling with a heat surface
            CVII = CVI1 
            SURFE= -SURFE1      ! Mersu
         ENDIF

C ... Blend the first- and second-order discretizations

         ALPHA   = 1.
         YP      = 1./DIST

         IF(TURCOR .AND. MULPHC /= 'MULTI') THEN
         UVEL    = SQRT((U(II)-USURF)**2 + (V(II)-VSURF)**2 + 
     +                  (W(II)-WSURF)**2)
         PARA    = YP*UVEL*RO(II)/VIS(II)	!   
c        F1      = (TANH((PARA-70.)/1.) + 1.)/2.
         ALPHA   = MIN((PARA)/140.,1.)
c        ALPHA   = 1. ! second order discretization can be forced
         ENDIF
         
         CVII    = (1.-ALPHA)*CVI1 + ALPHA*CVI2
         UVII    = (1.-ALPHA)*UVI1 + ALPHA*UVI2
         VVII    = (1.-ALPHA)*VVI1 + ALPHA*VVI2
         WVII    = (1.-ALPHA)*WVI1 + ALPHA*WVI2

         IF(.NOT.WALLFUNL) THEN
         DUC     = .33333*(A3Y(IF1)*VVII + A3Z(IF1)*WVII +A3X(IF1)*UVII)
         F2V     = SURF*(UVII + A3X(IF1)*DUC)
         F3V     = SURF*(VVII + A3Y(IF1)*DUC)
         F4V     = SURF*(WVII + A3Z(IF1)*DUC)
         
         ELSE
         
         F2V     =-F2VWF(K)
         F3V     =-F3VWF(K)
         F4V     =-F4VWF(K)
         
         EPS2X   = F2V/(MAX(EPS10,SURG*VISB*(UVII + A3X(IF1)*DUC)))
         EPS2Y   = F3V/(MAX(EPS10,SURG*VISB*(VVII + A3Y(IF1)*DUC)))
         EPS2Z   = F4V/MAX(EPS10,(SURG*VISB*(WVII + A3Z(IF1)*DUC)))
c         write(666,*)i,j,eps2x,eps2y,eps2z
c         EPS2(IM1) = 2.*MIN(TURLIM,MAX(1.,EPS2X,EPS2Y,EPS2Z))

         ENDIF
         
         F5VF    = USURF*F2V  + VSURF*F3V + WSURF*F4V
         F5VQ    = (1-IHFL)*SURFE*CVII + IHFL*HFLUX(NP)
         F5V     = F5VF + F5VQ 

         IF(MULPHL .AND. MULPHC == 'MULTI') THEN

         DO IP = 1,NPHASES

         IF(TURCOR) THEN
         UVEL    =SQRT((VAR(II)%U(IP)-USURF)**2+(VAR(II)%V(IP)-VSURF)**2
     +           +     (VAR(II)%W(IP)-WSURF)**2)
         PARA    = YP*UVEL*PRO(II)%RO(IP)/PRO(II)%VIS(IP)
c        F1      = (TANH((PARA-70.)/1.) + 1.)/2.
         ALPHA   = MIN((PARA)/140.,1.)
c        ALPHA   = 1. ! second order discretization can be forced
         ENDIF

         CVII    = (1.-ALPHA)*CV1(IP) + ALPHA*CV2(IP)
         UVII    = (1.-ALPHA)*UV1(IP) + ALPHA*UV2(IP)
         VVII    = (1.-ALPHA)*VV1(IP) + ALPHA*VV2(IP)
         WVII    = (1.-ALPHA)*WV1(IP) + ALPHA*WV2(IP)

         DUC     = .33333*(A3Y(IF1)*VVII + A3Z(IF1)*WVII +A3X(IF1)*UVII)
         F2V     = SURF*(UVII + A3X(IF1)*DUC)
         F3V     = SURF*(VVII + A3Y(IF1)*DUC)
         F4V     = SURF*(WVII + A3Z(IF1)*DUC)
         F5VF    = USURF*F2V  + VSURF*F3V + WSURF*F4V
         F5VQ    = (1-IHFL)*SURFE*CVII + IHFL*HFLUX(NP)
         F5V     = F5VF + F5VQ

c         VAR(IF1)%FRM(IP) = A3(IF1)*(A3X(IF1)*PB + F2V)
c         VAR(IF1)%FRN(IP) = A3(IF1)*(A3Y(IF1)*PB + F3V)
c         VAR(IF1)%FRW(IP) = A3(IF1)*(A3Z(IF1)*PB + F4V)

C ... Pressure gradient was separated 4.6.2018

         VAR(IF1)%FRM(IP) = A3(IF1)*F2V
         VAR(IF1)%FRN(IP) = A3(IF1)*F3V
         VAR(IF1)%FRW(IP) = A3(IF1)*F4V

         PRC(IF1)%FPM = A3(IF1)*A3X(IF1)*PB
         PRC(IF1)%FPN = A3(IF1)*A3Y(IF1)*PB
         PRC(IF1)%FPW = A3(IF1)*A3Z(IF1)*PB

         VAR(IF1)%FRM(IP) = VAR(IF1)%FRM(IP)*VAR(II)%ALFA(IP)
         VAR(IF1)%FRN(IP) = VAR(IF1)%FRN(IP)*VAR(II)%ALFA(IP)
         VAR(IF1)%FRW(IP) = VAR(IF1)%FRW(IP)*VAR(II)%ALFA(IP)

         F2VMF(IP)  = F2V;  F3VMF(IP)  = F3V;  F4VMF(IP) = F4V
         F5VFMF(IP) = F5VF; F5VQMF(IP) = F5VQ; F5VMF(IP) = F5V

         VAR(IF1)%FA(IP) = A3(IF1)*VAR(II)%ALFA(IP)
         VAR(IF1)%FX(IP) = 0.

         ENDDO

C ... These are weighted by the void fraction (area)

         F2V = 0.; F3V = 0.; F4V = 0.; F5VF = 0.; F5VQ = 0.; F5V = 0.

         DO IP = 1,NPHASES
            F2V  = F2V  + F2VMF(IP)* VAR(II)%ALFA(IP)
            F3V  = F3V  + F3VMF(IP)* VAR(II)%ALFA(IP)
            F4V  = F4V  + F4VMF(IP)* VAR(II)%ALFA(IP)
            F5VF = F5VF + F5VFMF(IP)*VAR(II)%ALFA(IP)
            F5VQ = F5VQ + F5VQMF(IP)*VAR(II)%ALFA(IP)
            F5V  = F5V  + F5VMF(IP)* VAR(II)%ALFA(IP)
         ENDDO
         ENDIF ! MULPHC == 'MULTI'

C .. Store the heat fluxes

         QWFRIC(NP)= F5VF  ! heat produced by friction
         QWALL(NP) = F5VQ  ! heat flux

         QY = QY + A3(IF1)*RSUM*F5V
         QZ = QZ + A3(IF1)*RSUM*F5VQ

C ... Diffusion of turbulence quantities

            EVII    = IDIR*(9.0*REPS(II)*YPRO1 - 1.0*REPS(II2)*YPRO2
     +              - 8.*EPSURF)
            RKII    = IDIR*(9.000*RK(II)*YPRO1 - 1.000*RK(II2)*YPRO2
     +              - 8.*RKSURF)
            RKII    = IDIR*MAX(IDIR*RKII,0.0) ! The wall cannot create turbulence


         IF(.NOT.TURCOR) THEN
            UTAU    = SQRT(SQRT(F2V**2 + F3V**2 + F4V**2)/ROB)

            F6V     = SURFRK*RKII*max(ROMFAC,XAUX1) ! TTB 15.3.2000
            F7V     = SURFEP*EVII*max(ROMFAC,XAUX1) ! TTB 15.3.2000

         ELSE

            UTAU    = SQRT(SQRT(F2V**2 + F3V**2 + F4V**2)/ROB)
            EVIK    = IDIR*6.*YP*MIN(0.,-ROB*UTAU**3/XKAPPA/YP**2 + 
     +                4.*VISB*UTAU**2/SQRT(CMU)/YP**3)*
     +               (1.+VIST(II)/VIS(II)*PSIGE)

            F6V     = SURF*RKII*.5*(1. + EPS2(II))
            F7V     = SURF*(ROMFAC*EVII + REPFAC*EVIK)

C ... Modify also the production of turbulence

         IF(.NOT.WALLFUNL)
     +   PTUR(II) = VIST(II)*(RO(II)*UTAU**2)**2/(VIS(II)+VIST(II))**2
           
         ENDIF

C ... SHEAR STRESSES ALONG THE GRID LINES

         ANXI     = XC(II+ISTR) - XC(II-ISTR) 
         ANYI     = YC(II+ISTR) - YC(II-ISTR) 
         ANZI     = ZC(II+ISTR) - ZC(II-ISTR)
         XNORI    = 1./SQRT(ANXI**2 + ANYI**2 + ANZI**2+EPS3)

         ANXJ     = XC(II+JSTR) - XC(II-JSTR) 
         ANYJ     = YC(II+JSTR) - YC(II-JSTR) 
         ANZJ     = ZC(II+JSTR) - ZC(II-JSTR)
         XNORJ    = 1./SQRT(ANXJ**2 + ANYJ**2 + ANZJ**2+EPS3)

         TAUW1(NP)= -(ANXI*F2V + ANYI*F3V + ANZI*F4V)*XNORI
         TAUW2(NP)= -(ANXJ*F2V + ANYJ*F3V + ANZJ*F4V)*XNORJ

         TAUWX(NP)= -F2V
         TAUWY(NP)= -F3V
         TAUWZ(NP)= -F4V

C ... INVISCID FLUXES:  EXTRAPOLATION OF THE PRESSURE WAS MOVED UP

         F3R(IF1)  = 0.
         IF(MULPHL) THEN
           VAR(IF1)%FRO(1) = 0.
           VAR(IF1)%FRO(2) = 0.
         ENDIF

C ... Intermittency variables
         IF(TRANSL) THEN
            TRM(IF1)%FG   = 0.
            TRM(IF1)%FRET = 0.
         ENDIF

C ***** WARNING!! FOR TITAN PROBE **************************************
C     IF(KBOT > 61) A3(IF1) = 0.
C***********************************************************************

C ... Store the fluxes

         F3RM(IF1) = A3(IF1)*(A3X(IF1)*PB + F2V)
         F3RN(IF1) = A3(IF1)*(A3Y(IF1)*PB + F3V)
         F3RW(IF1) = A3(IF1)*(A3Z(IF1)*PB + F4V)
         F3E(IF1)  = A3(IF1)*(UROT(IF1)*PB1 + F5V + F6V)
         QWFRIC(NP)= QWFRIC(NP) + UROT(IF1)*PB
         QY        = QY + RSUM*UROT(IF1)*PB

         IF(FRESUL) F1H(IF1) = 0.

         IF(MULPHL) THEN
C ... Wall boiling
           DREF  = 0.6E-3; DMAX = 1.4E-3; DTREF = 45.
           RNREF = 0.8*9.922E5 
           DTSUB = MAX(0.,PRO(II)%TSAT-T0)
           DTSUP = MAX(0.,PRO(II)%DTEMP(1)-PRO(II)%TSAT)
           DW    = MIN(DREF*EXP(-DTSUB/DTREF),DMAX) ! (9)
           RN    = RNREF*(DTSUP/10.)**1.805 ! (11) muutoksia
           F     = SQRT(1.333*G0*(PRO(II)%RO(1)-PRO(II)%RO(2))/
     +             (DW*PRO(II)%RO(2)))
           A2E   = MIN(5.,PII*DW**2*RN)
           GW    = A2E*PII/6.*DW*PRO(II)%RO(2)*F*RN
           HFG   = PRO(II)%HSAT(2)-PRO(II)%HSAT(1)
           GW    = MIN(IDIR*VAR(IF1)%X(1)*(F5V + F6V),GW*HFG)/HFG
           GWA   = A3(IF1)/VOL(II)*GW
           VAR(II)%EVAPR(1) = -GWA ! RPI model for evaporation
            VAR(II)%EVAPR(1) = 0.  ! RPI model neglected

c           VAR(II)%EVAPR(1) = A3(IF1)/VOL(II)*MIN(1.E-2*VAR(II)%X(1),
c     +     F5V/(PRO(II)%HSAT(2)-PRO(II)%HSAT(1)))
c           F5V = F5V - VOL(II)/A3(IF1)*
c     +     VAR(II)%EVAPR(1)/(PRO(II)%HSAT(2)-PRO(II)%HSAT(1))
c           F3E(IF1) = F3E(IF1) - VOL(II)/A3(IF1)*
c     +     VAR(II)%EVAPR(1)/(PRO(II)%HSAT(2)-PRO(II)%HSAT(1))

           VAR(IF1)%FE(1) = A3(IF1)*(VAR(II)%ALFA(1)*UROT(IF1)*PB + 
c     +                   (F5V + F6V))
     +                   VAR(II)%X(1)*(F5V + F6V)) ! Heat to masses
           VAR(IF1)%FE(2) = A3(IF1)*(VAR(II)%ALFA(2)*UROT(IF1)*PB + 
c     +                   (F5V + F6V))
     +                   VAR(II)%X(2)*(F5V + F6V))
c           VAR(IF1)%FE(2) = F3E(IF1) - VAR(IF1)%FE(1) ! Uusi heh heh
c           VAR(IF1)%FE(1) = VAR(IF1)%FE(1) - A3(IF1)*IDIR*GW*HFG ! heh heh

c           VAR(IF1)%FE(2) = F3E(IF1) - VAR(IF1)%FE(1)
           VAR(II)%EVAPR(2) = -VAR(II)%EVAPR(1)

         ENDIF


         F3RK(IF1)  = A3(IF1)*F6V
         F3EPS(IF1) = A3(IF1)*F7V

         QWFRIC(NP)= IDIR*QWFRIC(NP)  ! The positive sign is always outwards
         QWALL(NP) = IDIR*QWALL(NP)

1000  CONTINUE

C ... Then a flux correction next to the solid surfaces

      IF(TURCOR .AND. ITURB < 9) THEN ! Does not work for the RSM-based model

      DO J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR  + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID             ! 2. comp cell
         II3     = II + 2*LSTRID           ! 3. comp cell
         IF1 = II  + NFL                ! FLUX INDEX
         IF2 = IF1 + LSTRID             ! 2. face
         IF3 = IF1 + 2*LSTRID           ! 3. face

C ... Correct the fluxes next to the surfaces 2. cell face

         YPD2   = A3(IF2)/(VOL(II2) + VOL(II))
         SURG   = -A3(IF2)*YPD2

         SURF   = SURG*(-.33333*VIST(II2)+VIST(II))
         SURFE  = SURG*((-.33333*PRS*VIST(II2)/VIS(II2))*CH(II2) + 
     +        (PRS*VIST(II)/VIS(II))* CH(II))
         DUC    = .333333*(A3X(IF2)*(U(II2) - U(II)) +
     +        A3Y(IF2)*(V(II2) - V(II)) + A3Z(IF2)*(W(II2)-W(II)))
         
         F2VIG  = SURF*(U(II2)-U(II) + A3X(IF2)*DUC)
         F3VIG  = SURF*(V(II2)-V(II) + A3Y(IF2)*DUC)
         F4VIG  = SURF*(W(II2)-W(II) + A3Z(IF2)*DUC)
         F5VIG  = .5*((U(II2)+U(II))*F2VIG+(V(II2)+V(II))*F3VIG+
     +        (W(II2)+W(II))*F4VIG) + SURFE*(TEMP(II2) - TEMP(II))
         
         F3RM(IF2) = F3RM(IF2) + F2VIG
         F3RN(IF2) = F3RN(IF2) + F3VIG
         F3RW(IF2) = F3RW(IF2) + F4VIG
         F3E(IF2)  = F3E(IF2)  + F5VIG

         IF(MULPHL .AND. MULPHC == 'MULTI') THEN

         DO IP = 1,NPHASES

         SURF   = SURG*(-.33333*VIST(II2)+VIST(II))
         SURFE  = SURG*((-.33333*PRS*VIST(II2)/VIS(II2))*CH(II2) + 
     +        (PRS*VIST(II)/VIS(II))* CH(II))
         DUC    = .333333*(A3X(IF2)*(VAR(II2)%U(IP) - VAR(II)%U(IP)) +
     +                     A3Y(IF2)*(VAR(II2)%V(IP) - VAR(II)%V(IP)) +
     +                     A3Z(IF2)*(VAR(II2)%W(IP) - VAR(II)%W(IP)))
         F2VIG  = SURF*(VAR(II2)%U(IP)-VAR(II)%U(IP) + A3X(IF2)*DUC)
         F3VIG  = SURF*(VAR(II2)%V(IP)-VAR(II)%V(IP) + A3Y(IF2)*DUC)
         F4VIG  = SURF*(VAR(II2)%W(IP)-VAR(II)%W(IP) + A3Z(IF2)*DUC)
         F5VIG  = .5*((VAR(II2)%U(IP)+VAR(II)%U(IP))*F2VIG  +
     +                (VAR(II2)%V(IP)+VAR(II)%V(IP))*F3VIG  +
     +                (VAR(II2)%W(IP)+VAR(II)%W(IP))*F4VIG) + 
     +          SURFE*(PRO(II2)%TEMP(IP) - PRO(II)%TEMP(IP))

         VAR(IF2)%FRM(IP) = VAR(IF2)%FRM(IP) + VAR(IF2)%ALFA(IP)*F2VIG
         VAR(IF2)%FRN(IP) = VAR(IF2)%FRN(IP) + VAR(IF2)%ALFA(IP)*F3VIG
         VAR(IF2)%FRW(IP) = VAR(IF2)%FRW(IP) + VAR(IF2)%ALFA(IP)*F4VIG
         VAR(IF2)%FE(IP)  = VAR(IF2)%FE(IP)  + VAR(IF2)%ALFA(IP)*F5VIG

         ENDDO

         ENDIF ! MULTI

C ... Correct the fluxes next to the surfaces 3. cell face

         YPD2   = A3(IF3)/(VOL(II3) + VOL(II2))
         SURG   = -A3(IF3)*YPD2

         SURF   = SURG*(-.2*VIST(II3)+.3333*VIST(II2))
         SURFE  = SURG*((-.2*PRS*VIST(II3)/VIS(II3))*CH(II3) + 
     +        (+.3333*PRS*VIST(II2)/VIS(II2))* CH(II2))
         DUC    = .2*(A3X(IF3)*(U(II3) - U(II2)) +
     +        A3Y(IF3)*(V(II3) - V(II2)) + A3Z(IF3)*(W(II3)-W(II2)))
         
         F2VIG  = SURF*(U(II3)-U(II2) + A3X(IF3)*DUC)
         F3VIG  = SURF*(V(II3)-V(II2) + A3Y(IF3)*DUC)
         F4VIG  = SURF*(W(II3)-W(II2) + A3Z(IF3)*DUC)
         F5VIG  = .5*((U(II3)+U(II2))*F2VIG+(V(II3)+V(II2))*F3VIG+
     +        (W(II3)+W(II2))*F4VIG) + SURFE*(TEMP(II3) - TEMP(II2))
         
         F3RM(IF3) = F3RM(IF3) + F2VIG
         F3RN(IF3) = F3RN(IF3) + F3VIG
         F3RW(IF3) = F3RW(IF3) + F4VIG
         F3E(IF3)  = F3E(IF3)  + F5VIG

         IF(MULPHL .AND. MULPHC == 'MULTI') THEN

         DO IP = 1,NPHASES

         SURF   = SURG*(-.2*VIST(II3)+.3333*VIST(II2))
         SURFE  = SURG*((-.2*PRS*VIST(II3)/VIS(II3))*CH(II3) + 
     +        (+.3333*PRS*VIST(II2)/VIS(II2))* CH(II2))

         DUC    = .2*(A3X(IF3)*(VAR(II3)%U(IP) - VAR(II2)%U(IP)) +
     +                A3Y(IF3)*(VAR(II3)%V(IP) - VAR(II2)%V(IP)) +
     +                A3Z(IF3)*(VAR(II3)%W(IP) - VAR(II2)%W(IP)))
         F2VIG  = SURF*(VAR(II3)%U(IP)-VAR(II2)%U(IP) + A3X(IF3)*DUC)
         F3VIG  = SURF*(VAR(II3)%V(IP)-VAR(II2)%V(IP) + A3Y(IF3)*DUC)
         F4VIG  = SURF*(VAR(II3)%W(IP)-VAR(II2)%W(IP) + A3Z(IF3)*DUC)
         F5VIG  = .5*((VAR(II3)%U(IP)+VAR(II2)%U(IP))*F2VIG  +
     +                (VAR(II3)%V(IP)+VAR(II2)%V(IP))*F3VIG  +
     +                (VAR(II3)%W(IP)+VAR(II2)%W(IP))*F4VIG) + 
     +          SURFE*(PRO(II3)%TEMP(IP) - PRO(II2)%TEMP(IP))

         VAR(IF3)%FRM(IP) = VAR(IF3)%FRM(IP) + VAR(IF3)%ALFA(IP)*F2VIG
         VAR(IF3)%FRN(IP) = VAR(IF3)%FRN(IP) + VAR(IF3)%ALFA(IP)*F3VIG
         VAR(IF3)%FRW(IP) = VAR(IF3)%FRW(IP) + VAR(IF3)%ALFA(IP)*F4VIG
         VAR(IF3)%FE(IP)  = VAR(IF3)%FE(IP)  + VAR(IF3)%ALFA(IP)*F5VIG

         ENDDO

         ENDIF ! MULTI

      ENDDO
      ENDDO
      ENDIF ! TURB <= 9

      ELSE ! ONLY INVISCID SURFACE FLUXES

      DO 1100 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
c     NR      = (JN+J-KY1)*(KX2-KX1+1+2*IN) + IN
*      NR = (JN+J-KY1)*(KX2-KX1+1+2*IN) - KX1 + IN + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1100 I = KX1,KX2
         II      = (IN+I-1)*ISTR  + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID             ! CELL NEXT TO SURFACE
         IF1 = II + NFL                    ! FLUX INDEX
         NI      = I  + N1                 ! New patch index
         NP  = I  + NR                     ! Patch index with LN ghost cells

C ... ONLY INVISCID FLUXES:  EXTRAPOLATION OF THE PRESSURE

         PB1     = MAX(0.01*FRSPRE, RC1*P(II) - RC2*P(II2))
         PB4     = MAX(-.99*FRSPRE, RC1*PDI(II) - RC2*PDI(II2))
         PB      = YP1*PB1 + YP4*PB4
         ROB     = MAX(EPS         ,RC1*RO(II)- RC2*RO(II2))
         CPWALL(NP) = PB4

         F3R(IF1)  = 0.
         F3RM(IF1) = A3(IF1)*A3X(IF1) *PB
         F3RN(IF1) = A3(IF1)*A3Y(IF1) *PB
         F3RW(IF1) = A3(IF1)*A3Z(IF1) *PB
         F3E(IF1)  = A3(IF1)*UROT(IF1)*PB1
         F3RK(IF1) = 0.
         F3EPS(IF1)= 0.
         IF(FRESUL) F1H(IF1) = 0.

       
         TAUW1(NP)  = 0.       ! Inviscid flow
         TAUW2(NP)  = 0.
         QWFRIC(NP) =	0.
         QWALL(NP)  =	0.
1100  CONTINUE

      ENDIF  ! IF(IDIR > KBOT)
      

C ... *******   START OF POROSITY LOOP   *******   *******

C ... Before computing the forces acting on the surface,
C ... scale the flux values with possible porosity effect

      IF(OMPOR /= 1.) THEN

      DO  J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
c     NR      = (JN+J-KY1)*(KX2-KX1+1+2*IN) + IN
*      NR = (JN+J-KY1)*(KX2-KX1+1+2*IN) - KX1 + IN + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO  I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IF1 = II  + NFL                  ! FLUX INDEX
         IF2 = IF1 + LSTRID               ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP  = I   + NR                   ! Patch index with LN ghost cells

C ... OMPOR = (1.-POROS)

C .. The heat flux components
         QWFRIC(NP)= QWFRIC(NP)*OMPOR  ! heat produced by friction
         QWALL(NP) = QWALL(NP)*OMPOR ! heat flux

C ... Shear stresses along the grid lines
         TAUW1(NP)= TAUW1(NP)*OMPOR
         TAUW2(NP)= TAUW2(NP)*OMPOR

         TAUWX(NP)= TAUWX(NP)*OMPOR
         TAUWY(NP)= TAUWY(NP)*OMPOR
         TAUWZ(NP)= TAUWZ(NP)*OMPOR

C ... The fluxes
         F3RM(IF1) = F3RM(IF1)*OMPOR
         F3RN(IF1) = F3RN(IF1)*OMPOR
         F3RW(IF1) = F3RW(IF1)*OMPOR
         F3E(IF1)  = F3E(IF1)*OMPOR
         F3RK(IF1) = F3RK(IF1)*OMPOR
         F3EPS(IF1)= F3EPS(IF1)*OMPOR

      ENDDO
      ENDDO 

C ... Total heat produced by friction and heat flux.
C ... See FORCES file for Qwall

         QY        = QY*OMPOR
         QZ        = QZ*OMPOR
       
      ENDIF  !(POROS /= 0.)

C ... *******   END OF POROSITY LOOP   *******     *******    

C ... COMPUTE FORCES ACTING ON THE SURFACE

      IF(M == 1) THEN

      DO 2000 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
c     NR      = (JN+J-KY1)*(KX2-KX1+1+2*IN) + IN
*      NR = (JN+J-KY1)*(KX2-KX1+1+2*IN) - KX1 + IN + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 2000 I = KX1,KX2
         II      = (IN+I-1)*ISTR  + IJ + 1  ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID              ! CELL NEXT TO SURFACE
         IF1 = II + NFL                     ! FLUX INDEX
         NI      = I  + N1                  ! PATCH VELOCITY INDEX
         NH      = I  + NN                  ! New patch index
         NP  = I  + NR                      ! Patch index with LN ghost cells


C ... Centerpoint of the cell surface

c         X=.25*REAL(XCO(IF1)+XCO(IF+ISTR)+XCO(IF+JSTR)+XCO(IF+ISTR+JSTR))
c         Y=.25*REAL(YCO(IF1)+YCO(IF+ISTR)+YCO(IF+JSTR)+YCO(IF+ISTR+JSTR))
c         Z=.25*REAL(ZCO(IF1)+ZCO(IF+ISTR)+ZCO(IF+JSTR)+ZCO(IF+ISTR+JSTR))

         X = XCP(NP)
         Y = YCP(NP)
         Z = ZCP(NP)

C ... VECTOR FOR THE MOMENTUM. CENTERPOINT COORDINATES ARE UTILIZED

         XD = XCP(NP) - XMOM
         YD = YCP(NP) - YMOM
         ZD = ZCP(NP) - ZMOM

         DELTAP = FRSDEN*(GX*(XCP(NP) - GROUND) +
     &                    GY*(YCP(NP) - GROUND) +
     &                    GZ*(ZCP(NP) - GROUND))

         IF(ALTITUDE >= GROUND) DELTAP = 0.0

C ... INCLUDE POROSITY 
C ... (Question about the role of pressure across porous wall!)   
         BX    = -IDIR*(F3RM(IF1)-A3(IF1)*OMPOR*A3X(IF1)*YPG) ! Oli YPG
         BY    = -IDIR*(F3RN(IF1)-A3(IF1)*OMPOR*A3Y(IF1)*YPG)
         BZ    = -IDIR*(F3RW(IF1)-A3(IF1)*OMPOR*A3Z(IF1)*YPG)
         BXDP  = -IDIR*(-A3(IF1)*OMPOR*A3X(IF1)*( -DELTAP))
         BYDP  = -IDIR*(-A3(IF1)*OMPOR*A3Y(IF1)*( -DELTAP))
         BZDP  = -IDIR*(-A3(IF1)*OMPOR*A3Z(IF1)*( -DELTAP))

         BX    = BX + BXDP
         BY    = BY + BYDP
         BZ    = BZ + BZDP

C ... Store surface forces

         SURFX(NP) = BX/(A3(IF1) + 1.E-20)
         SURFY(NP) = BY/(A3(IF1) + 1.E-20)
         SURFZ(NP) = BZ/(A3(IF1) + 1.E-20)

         CMX     = CMX + BZ*YD - BY*ZD
         CMY     = CMY + BX*ZD - BZ*XD
         CMZ     = CMZ + BY*XD - BX*YD

         CX      = CX  + BX
         CY      = CY  + BY
         CZ      = CZ  + BZ
         TOM     = TOM - IDIR*(F3E(IF1) - F3RK(IF1)) ! EFFECT OF TURBULENCE
                          
C ... AUXILIARY FORCE COMPONENTS

         PB1     = MAX(0.01*FRSPRE, RC1*P(II) - RC2*P(II2))
         PB4     = RC1*PDI(II) - RC2*PDI(II2)
         PB      = YP1*PB1 + YP4*PB4
         DXP      = - IDIR*(A3(IF1)*OMPOR*A3X(IF1)*(PB - YPG))
         DYP      = - IDIR*(A3(IF1)*OMPOR*A3Y(IF1)*(PB - YPG))
         DZP      = - IDIR*(A3(IF1)*OMPOR*A3Z(IF1)*(PB - YPG))
         DX      = DX  - IDIR*(A3(IF1)*OMPOR*A3X(IF1)*(PB - YPG))
         DY      = DY  - IDIR*(A3(IF1)*OMPOR*A3Y(IF1)*(PB - YPG))
         DZ      = DZ  - IDIR*(A3(IF1)*OMPOR*A3Z(IF1)*(PB - YPG))

C ... Pressure component of the force

         SURFPX(NP) = (DX+DXP)/(A3(IF1) + 1.E-20)
         SURFPY(NP) = (DY+DYP)/(A3(IF1) + 1.E-20)
         SURFPZ(NP) = (DZ+DZP)/(A3(IF1) + 1.E-20)

         BX      = BX + IDIR*(A3(IF1)*OMPOR*A3X(IF1)*(PB - YPG))
         BY      = BY + IDIR*(A3(IF1)*OMPOR*A3Y(IF1)*(PB - YPG))
         BZ      = BZ + IDIR*(A3(IF1)*OMPOR*A3Z(IF1)*(PB - YPG))

         VMX     = VMX + BZ*YD - BY*ZD
         VMY     = VMY + BX*ZD - BZ*XD
         VMZ     = VMZ + BY*XD - BX*YD

         BUX     = BUX + BXDP
         BUY     = BUY + BYDP
         BUZ     = BUZ + BZDP
         QX      = QX  - IDIR*F3E(IF1)

      IF(FRESUL) IWAVEB(II-LSTRID) = 8
        
2000  CONTINUE
 
      TOM     = TOM + IDIR*QZ ! EFFECT OF HEAT FLUX IS SUBSTRACTED
      QZ      = IDIR*QZ
      QY      = IDIR*QY
      DX      = DX + BUX ! Add hydrostatic component
      DY      = DY + BUY
      DZ      = DZ + BUZ
      ENDIF

      RETURN
      END SUBROUTINE BOTSKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOTSOL(F3R,F3RM,F3RN,F3RW,F3E,P,U,V,W,C,E,VOL,A3,A3X,
     2  A3Y,A3Z,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,
     3  KY1,KY2,PR,PRT,RO,IFLUX,M,IN,JN,KN,RKLIM,IDIR,F3RK,F3EPS,UROT,
     4  TEMP,CH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,ITURB,ISTR,JSTR,KSTR,
     5  IHEAT,RCON,IDIS,TWALL,QWALL,HFLUX,XCP,YCP,ZCP,INTEM,INTET,n,
     6  ihff,HAT1)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: IDIR,IDIS,ITURB,MAX11,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     2  KX1,KX2,KY1,KY2,IFLUX,M,IN,JN,KN,ISTR,JSTR,KSTR,NFL,IHFL,IHT,
     3  KA,I,J,NN,NR,NI,N1,II,IJ,II2,IF1,IF2,NH,NP,IZZ,IF3,II3,IA,N,
     4  IHEAT,LSTRID,ISTATE,INTEM,INTET,IM1,IM2,ihff

      REAL :: F3R(*),F3RW(*),F3E(*),P(*),V(*),U(*),C(*),RO(*),A3(*),
     2  A3Y(*),A3Z(*),A3X(*),VIS(*),EPS2(*),VIST(*),F3RN(*),W(*),
     3  VOL(*),F3RM(*),F3RK(*),F3EPS(*),UROT(*),TEMP(*),CH(*),RCON(*),
     4  E(*),QWALL(*),TWALL(*),HFLUX(*),HAT1(*)

      REAL :: PR,PRT,RKLIM,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     2  RC1,RC2,RM1,RM2,PRS,EPS,TSURF,VISLIM,CHLIM,ROLIM,DIST,
     3  SURG,EPS2B,VISTB,VISB,CHSB,ROB,SURF,SURFE,USURF,VSURF,
     4  WSURF,TEMP1,TEMP2,CVI2,UVI2,VVI2,WVI2,ALPHA,CVII,
     5  UVII,WVII,F5VF,F5V,CVI1,UVI1,VVI1,WVI1,VVII,PB1,PB,
     6  EB,TEMPB,HB,DISTF,CHS,SURGE

      REAL :: XMULTI(6)

      REAL :: XCP(*), YCP(*), ZCP(*)

C
C ... LU FACTORED METHOD. BOUNDARY CONDITION FOR THE SOLID BOUNDARY
C ... BOUNDARY FLUXES ARE CALCULATED EXPLICITLY. POINTS ARE NUMBERED
C ... XI-WISE. THIS IS FOR THE SOLID SIDE (TSii 8.11.2004)
C

C ... this is not exatly correct, but may work better in lower levels.

      PRS     = PR/PRT

      RC1     = 1.0  ! 1.500
      RC2     = RC1 - 1.
      IF(INTEM == 6) THEN
         RC1  = 1.
         RC2  = RC1 - 1.
      ENDIF

      IF (IDIS == 2) THEN        ! k-omega model
         RM1    = RC1
         RM2    = RC2
      ELSE                       ! k-epsilon model
C        RM1    = 1.500          ! Activate this one day
C        RM2    = RM1 - 1.0      ! Activate this one day
         RM1    = 1.000
         RM2    = RM1 - 1.0
      END IF

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

      EPS     = 1.E-3                            !  20.12.1994 : 10
      IHFL    = 0
      IHT     = 1

      IF(IHEAT == 2 .OR. IHEAT >= 12) IHFL  = 1  ! Heat flux is given
      IF(IHEAT == 0)  IHT  = 0
      IF(IHEAT == 4)  IHFL = 0
      IF(IHEAT == 15) IHT  = 0
      IF(IHEAT == 15) IHFL = 0

C ... IHEAT = 0 ADIABATIC
C ...         1 WALL TEMPERATURE IS GIVEN IN TWALL
C ...         2 HEAT FLUX IS GIVEN        IN HFLUX
C ...         3 WALL TEMPERATURE IS PUT TO FREE_STREAM STAGNATION TEMEPERATURE
C ...         4 COUPLING WITH A FLUID BLOCK

      VISLIM  = 1.E-8
      CHLIM   = 100.*FRSVIS
      ROLIM   = .001*FRSDEN

C ... Start with fluid-side contribution and then add the solid part

      KA      = (KN+KBOT-1)*KSTR
         
      IF(IDIM > KBOT) THEN

      DO 1000 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1000 I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! 1. ghost cell
         IM2     = IM1- LSTRID            ! 2. ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IF2 = IF1 + LSTRID                  ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP  = I   + NR                      ! Patch index with LN ghost cells

         DISTF   = (A3(IF1)+A3(IF1-LSTRID))/VOL(IM1)
         EPS2B   = MAX(1.    ,  RM1*EPS2(IM1) - RM2*EPS2(IM2))
         EPS2B   = 1.           ! works better
         VISTB   = MAX(0.    ,  RM1*VIST(IM1) - RM2*VIST(IM2))
         VISTB   = 0.           ! works better
         VISB    = MAX(VISLIM,  RC1* VIS(IM1) - RC2* VIS(IM2))
         CHSB    = MAX(EPS   ,  RC1*  CH(IM1) - RC2*  CH(IM2))
         TEMPB   = MAX(0.01*T0, RC1*TEMP(IM1) - RC2*TEMP(IM2))
         SURFE   = .5*DISTF*(PRS*VISTB/VISB+1.)*CHSB
         TEMP1   = TEMP(IM1)
         TEMP2   = TEMP(IM2)

         DIST    = (A3(IF1)+A3(IF1+LSTRID))/VOL(II)
         CHS     = MAX(EPS   ,  RC1*  CH(II)  - RC2*  CH(II2))
         SURGE   = .5*DIST*CHS

         IF(IHEAT == 1) THEN
            TSURF     = TWALL(NP)
         ELSE IF(IHEAT /= 3 .AND. IHEAT /= 4) THEN
            TSURF     = T0
            TWALL(NP) = TEMPB
         ELSE IF(IHEAT == 3) THEN
            TSURF     = TEMPB
            TWALL(NP) = TEMPB
         ELSE IF(IHEAT == 4) THEN
            TSURF     = (SURFE*TEMP(IM1) + SURGE*TEMP(II))/(SURFE+SURGE)
            TWALL(NP) = TSURF
         ENDIF
cc           write(997,*) J,KY1,Ky2,KX1,KX2
cc           write(997,*) TEMP(IM1),surfe,surge,Twall(np)
C ... First-order fluxes for this wall
         
         DIST    = (A3(IF1)+A3(IF1+LSTRID))/VOL(II)
         CHS     = MAX(EPS   ,RC1*  CH(II) - RC2*  CH(II2))
         TEMP1   = TEMP(II)
         TEMP2   = TEMP(II+LSTRID)
         SURGE   = -.5*DIST*CHS

C ... Heat flux on the surface

         CVI1     = IHT*IDIR*(TEMP1 - TSURF)
         F5V      = (1-IHFL)*SURGE*CVI1 + IHFL*HFLUX(NP)
          
C ... Store the fluxes

         F3R(IF1)  = 0.    ; HAT1(IF1) = F3R(IF1)   ! No suction or blowing ?
         F3RM(IF1) = 0.
         F3RN(IF1) = 0.
         F3RW(IF1) = 0.
         F3E(IF1)  = A3(IF1)*F5V
         F3RK(IF1) = 0.
         F3EPS(IF1)= 0.

1000  CONTINUE
      ENDIF

C ... STORE THE HEAT FLUXES ACTING ON THE SURFACE ? (In BOTSKE)

      RETURN
      END SUBROUTINE BOTSOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOTSFR(F3R,F3RM,F3RN,F3RW,F3E,PDIFF,P,U,V,W,C,E,VOL,
     2  A3,A3X,A3Y,A3Z,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     3  KX1,KX2,KY1,KY2,PR,PRT,RO,IFLUX,M,IN,JN,KN,RKLIM,IDIR,F3RK,
     4  F3EPS,UROT,TEMP,CH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,ITURB,ISTR,
     5  JSTR,KSTR,IHEAT,RCON,IDIS,TWALL,QWALL,HFLUX,XCP,YCP,ZCP,WMFLUX,
     6  CPWALL,SURLE,DSURLE,SURFX,SURFY,SURFZ,INTEM,INTET,IFSBC,n,HAT1,
     7  WAVEH,SURFPX,SURFPY,SURFPZ)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: IDIR,IDIS,ITURB,MAX11,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     2  KX1,KX2,KY1,KY2,IFLUX,M,IN,JN,KN,ISTR,JSTR,KSTR,NFL,IHFL,IHT,
     3  KA,I,J,NN,NR,NI,N1,II,IJ,II2,IF1,IF2,NH,NP,IZZ,IF3,II3,IA,N,
     4  IHEAT,LSTRID,ISTATE,INTEM,INTET,IM1,IM2,IFSBC

      REAL :: F3R(*),F3RW(*),F3E(*),P(*),V(*),U(*),C(*),RO(*),A3(*),
     2  A3Y(*),A3Z(*),A3X(*),VIS(*),EPS2(*),VIST(*),F3RN(*),W(*),
     3  VOL(*),F3RM(*),F3RK(*),F3EPS(*),UROT(*),TEMP(*),CH(*),RCON(*),
     4  E(*),QWALL(*),TWALL(*),HFLUX(*),WMFLUX(*),
     5  CPWALL(*),PDIFF(*),SURLE(*),DSURLE(*),HAT1(*),WAVEH(*),
     6  SURFX(*),SURFY(*),SURFZ(*),SURFPX(*),SURFPY(*),SURFPZ(*)

      REAL :: PR,PRT,RKLIM,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     2  RC1,RC2,RM1,RM2,PRS,EPS,TSURF,VISLIM,CHLIM,ROLIM,DIST,
     3  SURG,EPS2B,VISTB,VISB,CHSB,ROB,SURF,SURFE,USURF,VSURF,
     4  WSURF,TEMP1,TEMP2,CVI2,UVI2,VVI2,WVI2,ALPHA,CVII,
     5  UVII,WVII,F5VF,F5V,CVI1,UVI1,VVI1,WVI1,VVII,PB1,PB,
     6  EB,TEMPB,HB,DISTF,CHS,SURGE,PDM,PB4,DISTG,YMR,YPRO1,SURFE1,
     7  YP1,YP4,DELTAP

      REAL :: XMULTI(6)

      REAL :: XCP(*), YCP(*), ZCP(*)

C
C ... LU FACTORED METHOD. BOUNDARY CONDITION FOR THE FREE-SURFACE
C ... BOUNDARY FLUXES ARE CALCULATED EXPLICITLY. POINTS ARE NUMBERED
C ... XI-WISE. THIS IS FOR THE SOLID SIDE (TSii 2.10.2006)
C

C ... this is not exatly correct, but may work better in lower levels.

      PRS     = PR/PRT

      RC1     = 1!.500
      RC2     = RC1 - 1.
      IF(INTEM == 6) THEN
         RC1  = 1.
         RC2  = RC1 - 1.
      ENDIF

      YMR     = .166666667

      IF (IDIS == 2) THEN     ! k-omega model
         RM1    = RC1
         RM2    = RC2
      ELSE                       ! k-epsilon model
C        RM1    = 1.500          ! Activate this one day
C        RM2    = RM1 - 1.0      ! Activate this one day
         RM1    = 1.000
         RM2    = RM1 - 1.0
      END IF

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

      EPS     = 1.E-3                            !  20.12.1994 : 10
      IHFL    = 0
      IHT     = 1

      IF(IHEAT == 2 .OR. IHEAT >= 12) IHFL = 1  ! Heat flux is given
      IF(IHEAT == 0)  IHT  = 0
      IF(IHEAT == 15) IHT  = 0
      IF(IHEAT == 15) IHFL = 0

C ... IHEAT = 0 ADIABATIC
C ...         1 SURFACE TEMPERATURE IS GIVEN IN TWALL
C ...         2 SURFACE HEAT FLUX IS GIVEN   IN HFLUX
C ...         3 WALL TEMPERATURE IS PUT TO FREE_STREAM STAGNATION TEMEPERATURE

      VISLIM  = 0.01*FRSVIS
      CHLIM   = 100.*FRSVIS
      ROLIM   = .001*FRSDEN
      
C ... Start with viscous surface fluxes and then add the inviscid part

      KA      = (KN+KBOT-1)*KSTR
         
      DO 1000 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1000 I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! 1. ghost cell
         IM2     = IM1- LSTRID            ! 2. ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IF2 = IF1 + LSTRID                  ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NH      = I  + NN                ! New patch index
         NP  = I   + NR                      ! Patch index with LN ghost cells

         DIST    = (A3(IF1)+A3(IF1+LSTRID))/VOL(II)
         DISTG   = (A3(IF1)+A3(IF1-LSTRID))/VOL(IM1)
         SURG    =-DIST*YMR

         PB1     = MAX(0.01*FRSPRE,.5*(P(II)+P(IM1)))
         PB4     = MAX(-.99*FRSPRE, .5*(PDIFF(II) + PDIFF(IM1)))
         PDM     = .5*(PDIFF(II) + PDIFF(IM1))
         ROB     = MAX(ROLIM ,.5*  RO(II) - .5*  RO(IM1))
         TEMPB   = MAX(0.01*T0, .5*TEMP(II) - .5*TEMP(IM1))
         YPRO1   = 1./RO(II)
         TEMP1   = TEMP(II)
         TEMP2   = TEMP(II+LSTRID)

         IF(IHEAT == 1) THEN
            TSURF     = TWALL(NP)
         ELSE IF(IHEAT /= 3 .AND. IHEAT /= 4) THEN
            TSURF     = T0
            TWALL(NP) = TEMPB
         ELSE IF(IHEAT == 3) THEN
            TSURF     = TEMPB
            TWALL(NP) = TEMPB
         ENDIF

c      IF(IDIM > KBOT) THEN ! Viscous fluxes currently zeros
c      ENDIF


C ... INVISCID FLUXES:  EXTRAPOLATION OF THE PRESSURE

      IF(IFSBC == 4) THEN ! Evaluate the free-surface fluxes?

c         F3R(IF1)  = 0.
c         F3RM(IF1) = PB4
c         F3RN(IF1) = A3(IF1)*A3Y(IF1) *PB4
c         F3RW(IF1) = A3(IF1)*A3Z(IF1) *PB4
c         F3E(IF1)  = A3(IF1)*UROT(IF1)*PB1
c
c         DELTAP   = FRSDEN*G0*WAVEH(NP)
c         F3RM(IF1) = F3RM(IF1) + A3(IF1)*A3X(IF1)*DELTAP - HAT1(IF1)
c         F3RN(IF1) = F3RN(IF1) + A3(IF1)*A3Y(IF1)*DELTAP - HAT1(IF1)
c         F3RW(IF1) = F3RW(IF1) + A3(IF1)*A3Z(IF1)*DELTAP - HAT1(IF1)
c         F3RK(IF1) = 0.
c         F3EPS(IF1)= 0.
c         F3RM(IF1) = F3RM(IF1) - 0.05*A3(IF1)*
c     +   (-PDIFF(IM2)+3.*PDIFF(IM1)-3.*PDIFF(II)+PDIFF(II2))
c         F3RN(IF1) = F3RN(IF1) - 0.05*A3(IF1)*
c     +   (-PDIFF(IM2)+3.*PDIFF(IM1)-3.*PDIFF(II)+PDIFF(II2))
c         F3RW(IF1) = F3RW(IF1) - 0.05*A3(IF1)*
c     +   (-PDIFF(IM2)+3.*PDIFF(IM1)-3.*PDIFF(II)+PDIFF(II2))

c         F1H(IF1)  = 0.

      ENDIF  ! IFSBC == 4
       
C ... Store the fluxes
       
          WMFLUX(NP) = F3R(IF1) /(A3(IF1) + 1.E-20)
          CPWALL(NP) = PDM
c          DSURLE(NP) = P(II)
          TWALL(NP)  = TEMPB
          SURFX(NP)  = F3RM(IF1)
          SURFY(NP)  = F3RN(IF1)
          SURFZ(NP)  = F3RW(IF1)
1000  CONTINUE

C ... STORE THE HEAT FLUXES ACTING ON THE SURFACE ? (In BOTSKE)

      RETURN
      END SUBROUTINE BOTSFR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ACTUATOR_FORCES(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,UROT,POROS,
     9  RMLOSS,IF2,WMFLUX,SURLE,DSURLE,XCP,YCP,ZCP,SURFX,SURFY,SURFZ,
     1  IACTU,TOMEGA,CX,CY,CZ,CMX,CMY,CMZ,SURFA,SURFT,SURF2X,SURF2Y,
     2  SURF2Z,SURFMX,SURFMY,SURFMZ)

      USE CONSTANTS, ONLY : PII,EPS6

      USE INTEGERS, ONLY  : IREPEA

      USE NS3CO, ONLY     : ALPHA, BETA, ICYOLD, LN

      USE FLIGHT, ONLY : XCG,YCG,ZCG,RIN,ROUT,THRUST,TORQUE,ADV,
     &    ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI,ACTUA,RIN,ROUT, ! Actuator disc
     &   IFA,IFT,IFR,NSERIES,VTIP,CDBLADE,CFBLADE,SIGMA,CFTI,SHIP,
     &   OSKU,NBLADE,QFACT

      USE BLADE_VARIABLES, ONLY : CBLADE

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IF2,INTER,INTE2,IFN,IACTU,NSERIE,L,NP2,NM2P2,NINTERP,TQINT,
     6  KERTA

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  UROT(*),POROS(*),RMLOSS(*),WMFLUX(*),SURLE(*),DSURLE(*),
     1  SURFX(*),SURFY(*),SURFZ(*),SURFA(*),
     2  SURFT(*),SURF2X(*),SURF2Y(*),SURF2Z(*),SURFMX(*),SURFMY(*),
     3  SURFMZ(*)

      REAL :: UM,VM,WM,ROM,RELAX,RK1,RK2,RC1,RC2,RM1,RM2,YP1,YP4,
     1  RGAS,GAMMA,E0REF,T0REF,TURLIM,FRSSIE,PRS,PR,PRT,UC,
     2  PLOSS,UH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,PLOSS2,R1,R2,
     3  CX,CY,CZ,CMX,CMY,CMZ,TOMEGA,ALA,X,Y,Z,BX,BY,BZ,
     4  XD,YD,ZD,RX,RY,RZ,RX1,RY1,RZ1,PLOSS2M,A0,A1,B1,PSI,PLEN,RLEN,
     5  SACB,CACB,SB,RX2,RY2,RZ2,PLOSS2R,ALPHA_DISC,THETA,FORDIS,
     6  APU1,APU2,APU3,SQRTHE,CORR,APUX,APUY,APUZ,CMX1,CMY1,CMZ1,PSIL,
     7  DFA,FMU,XXX,FTPA,FTIA,QCI,QI,QP,RLAM,VIO,APU5,PHI,rR,rR1,rR2,
     8  UT,APU6,U10,V10,W10,VA,VT,VX,VY,VZ,NX,NY,NZ,NLEN,TX,TY,TZ,VLEN,
     9  AAD,BAD,CAD,DAD,EAD,FAD,GAD,DIST,DISTT1,DISTT2,DISTQ1,DISTQ2,
     1  THEANG1,Cldist,Xloc,Yloc,Zloc,Xac,Yac,Zac,alphai,VE,COM,CL,CD,
     2  NXBLE,NYBLE,NZBLE,PXBLE,PYBLE,PZBLE,RVECLEN,CBLA,VRAD

      REAL, ALLOCATABLE :: CC(:)

      REAL ::  PLOSS2X,PLOSS2Y,PLOSS2Z,
     1         PLOSS2MX,PLOSS2MY,PLOSS2MZ,PLOSS2RX,PLOSS2RY,PLOSS2RZ,
     2         TVECX,TVECY,TVECZ,Xfin,Yfin,Zfin

      REAL :: XCP(*), YCP(*), ZCP(*)
      REAL :: XMOM,YMOM,ZMOM

      CHARACTER(LEN=15) :: HEADER

      LOGICAL :: THERE

C ... Calculate actuator forces based on FLIGHT(IACTU) variables. The
C ... forces are calculated patchwise and stored as F/A. The same force
c ... is calculated for the lower and upper surfaces of the disc. These
C ... are applied in ACTUATOR_DISC

C **********************************************************************

      RELAX  = 1. ! Useless?
      kerta  = irepea(20)
c      write(6,*) kerta
C ... Actuator disc parameters

C **********************************************************************

      A0     = CONEA(IACTU)
      A1     = ROTA1(IACTU)

C ... Check maximum and minimum angles

      IF(IFA(IACTU) == 3 .AND. A1 > 1.5707 .AND. A1 < 1.5709 .AND. 
     &     .NOT.SHIP) THEN
         A1  = 1.5707
      ELSEIF(IFA(IACTU) == 3 .AND. A1 > 1.5709 .AND. .NOT.SHIP) THEN
         WRITE(*,*)'MAXIMUM VALUE FOR ROTA1 IS 90deg, CURRENT ROTA1',
     &        A1*180./PII
         WRITE(13,*)' MAXIMUM VALUE FOR ROTA1 IS 90deg, CURRENT ROTA1',
     &        A1*180./PII
         STOP
      ENDIF
      IF(IFA(IACTU) == 3 .AND. A1 < -1.5707 .AND. A1 > -1.5709 .AND.
     &     .NOT. SHIP) THEN
         A1  = -1.5707
      ELSEIF(IFA(IACTU) == 3 .AND. A1 < -1.5709 .AND. .NOT.SHIP) THEN
         WRITE(*,*)'MINIMUM VALUE FOR ROTA1 IS -90deg, CURRENT ROTA1',
     &        A1*180./PII
         WRITE(13,*)' MINIMUM VALUE FOR ROTA1 IS -90deg, CURRENT ROTA1',
     &        A1*180./PII
         STOP
      ENDIF

      B1 = ROTB1(IACTU)   
      R2 = ROUT(IACTU)
      R1 = RIN(IACTU)

C ... Free-stream velocity vector

      U10 = FRSVEL*COS(ALPHA)*COS(BETA)
      V10 = FRSVEL*SIN(ALPHA)
      W10 = FRSVEL*COS(ALPHA)*SIN(BETA)

      NSERIE = NSERIES(IACTU) ! Tälle arametri FLIGHTistä inputtiin !!

      ALLOCATE(CC(0:NSERIE)) 

      LSTRID = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      CMX1=0. ; CMY1=0 ; CMZ1=0.


*      DO J = KY1+1,KY2-1  ! These should not be extended, obs. below

      DO J = KY1,KY2  ! Nowadays these are not extended

      IJ   = (JN+J-1)*JSTR + KA
*      NR = (JN+J-(KY1+1))*((KX2-1)-(KX1+1)+1+2*IN) - (KX1+1) + IN +1 !??
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
*      DO I = KX1+1,KX2-1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IFN = IF1 + LSTRID                  ! 2. face FLUX
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

****************************************************************************

C ... Subsonic flow is assumed. Do we ever need these? (Local conditions)

          UM = .5*(U(II) + U(IM1)) - UROT(IF1)*A3X(IF1)
          VM = .5*(V(II) + V(IM1)) - UROT(IF1)*A3Y(IF1)
          WM = .5*(W(II) + W(IM1)) - UROT(IF1)*A3Z(IF1)
          ROM= .5*(RO(II)+RO(IM1))
          COM= .5*(C(II)+C(IM1))


          UC = A3X(IF1)*UM + A3Y(IF1)*VM + A3Z(IF1)*WM - UROT(IF1)
          UH = UC/(POROS(NP)+1.E-6)

c ... Thrust vector (t) in the FINFLO coord.

          CALL EULFIN(0.0,B1,A1,0.0,0.0,-1.0,TVECX,TVECY,TVECZ,2)

          IF (SHIP) THEN 
             CALL EULFIN(0.0,A1,B1,0.0,0.0,-1.0,TVECX,TVECY,TVECZ,1)
          ENDIF

          TX = TVECX
          TY = TVECY
          TZ = TVECZ

C ... Free-stream velocity components for the disc

          VA = -(TX*U10 + TY*V10 + TZ*W10)
          VT = SQRT(FRSVEL**2 - VA**2)

C ... Extra velocity

          UT        = SQRT(.5*THRUST(IACTU)/(FRSDEN*ACTUA(IACTU)))
          APU6      = -0.5*FRSVEL**2+.5*SQRT(FRSVEL**4+4*UT**4)
          VIO       = SQRT(APU6)

C ... Angle of attack for the disk

          ALPHA_DISC= ATAN((VA + VIO)/(VT+EPS6))

C ... Bramwell is inconsistent, remove some day (sob sob)
C          ALPHA_DISC= (-ALPHA-A1)*COS(B1) + BETA*SIN(B1)

C ... A vector (P) between the hub and point CP in the FINFLO coord. 
C     ??? disc coordinates ???

          RX        = XCP(NP)-XCG(IACTU)
          RY        = YCP(NP)-YCG(IACTU)
          RZ        = ZCP(NP)-ZCG(IACTU)
          RVECLEN   = SQRT(RX**2+RY**2+RZ**2)
          PLEN      = 1./RVECLEN

C ... A dimensionless distance

          RX1       = RX*PLEN
          RY1       = RY*PLEN
          RZ1       = RZ*PLEN

C ... Vector product r=txP Distance from the hub (?)

          RX2       = RZ*TY - RY*TZ
          RY2       = RX*TZ - RZ*TX
          RZ2       = RY*TX - RX*TY
          RLEN      = SQRT(RX2**2 + RY2**2 + RZ2**2)

C ... A flow field vector (V) in the FINFLO coord.

          VX        = FRSVEL*COS(BETA)*COS(ALPHA)
          VY        = FRSVEL*COS(BETA)*SIN(ALPHA)
          VZ        = FRSVEL*SIN(BETA)
          IF (IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19  .OR.
     +        IFT(IACTU) >= 11 .AND. IFT(IACTU) <= 19) THEN
             VX = 1.0
             VY = 0.0
             VZ = 0.0
          ENDIF !!! JOS TAMAN KOMMENTOI, NIIN AO. !!! PITAA OTTA KAYTTOON 

C ... Vector product n=txV

          NX        = VZ*TY - VY*TZ
          NY        = VX*TZ - VZ*TX
          NZ        = VY*TX - VX*TY
          NLEN      = SQRT(NX**2 + NY**2 + NZ**2)

C ... PSI angle

C          IF (IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19  .OR. !!!
C     +         IFT(IACTU) >= 11 .AND. IFT(IACTU) <= 19) THEN
C             IF(RZ <= 0.) THEN
C                PSIL = ACOS(RX/SQRT(RX**2+RZ**2))             
C             ELSEIF(RZ > 0.) THEN 
C                PSIL = PII+ACOS(-RX/SQRT(RX**2+RZ**2))                
C             ENDIF
C          ELSE                                               !!!
          PSI     =  ACOS((RX2*NX+RY2*NY+RZ2*NZ)/(RLEN*NLEN))
          IF ((RX2*VX+RY2*VY+RZ2*VZ) < 0) THEN
             PSIL   = PSI
          ELSEIF ((RX2*VX+RY2*VY+RZ2*VZ) >= 0) THEN
             PSIL   = 2*PII-PSI
          ENDIF
CCCC          ENDIF

C ... Coeffcients for the blade element

          IF (IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19  .OR.
     +        IFT(IACTU) >= 11 .AND. IFT(IACTU) <= 19) THEN
             CALL BLELEMENT(PSIL,RVECLEN,R1,R2,UM,VM,WM,COM,CL,CD,
     +            ALPHAI,VE,IACTU,NXBLE,NYBLE,NZBLE,PXBLE,PYBLE,PZBLE,
     +            RX,RY,RZ)
          ENDIF

C **********************************************************************

C ... Axial force

C **********************************************************************

      IF(IFA(IACTU) == 0) THEN ! No axial force (for testing)

          PLOSS2 = 0. 

      ELSEIF(IFA(IACTU) == 1) THEN ! A constant distribution

          PLOSS2    = THRUST(IACTU)/ACTUA(IACTU)  ! Force/area
                                                            !!! Yhtälö (13)
      ELSEIF(IFA(IACTU) == 2 .OR. IFA(IACTU) == 3) THEN ! Bramwell's distribution

          THETA     = SQRT(1.-(RLEN/R2)**2) ! Onko OK?
          DFA       = 0. ! Correction

C.... Momentum source

          CC(0:NSERIE) = 0.

C ... Use *SQRT(LAUSEKE) - APU to simplify 

          APU1      = 1.-THETA**2
          APU2      = (1-SIN(ALPHA_DISC))/(1+SIN(ALPHA_DISC))
          APU3      = SQRT(APU2)
          SQRTHE    = SQRT(1.-THETA**2)

          CC(0)     = 15*THETA*APU1*.125 ! Voi mersu tota nollaa
          CC(1)     =-15*PII/256*(5.-9.*THETA**2)*SQRTHE*APU3
          CC(3)     =-45*PII/256*APU1*SQRTHE*APU2*APU3

          DO L      = 2,NSERIE,2 ! Heh heh
             NP2    = L/2 ; NM2P2 = (L-2)/2
             C(L)   = (-1)**NM2P2*15./8.*((THETA+L)/(THETA**2-1)*
     +       (9.*THETA**2+L**2-6.)/(L**2-9.)+3.*THETA/(L**2-9))*
     +       ((1.-THETA)/(1.+THETA))**NP2*APU2**NP2
          ENDDO

          FORDIS    = .5*CC(0)

          DO L = 1,NSERIE
             FORDIS = FORDIS - CC(L)*COS(REAL(L)*PSIL) ! PSI vai PSIL ??
          ENDDO

          PLOSS2    = 4.*FORDIS*THRUST(IACTU)/(PII*R2**2)

          IF(IFA(IACTU) == 3) THEN
          
             XXX = (1./PLEN)/R2
             IF(XXX >= 0.22 .AND. XXX <= 0.8) THEN
                FMU = 0.9781-11.147*XXX+37.338*XXX**2-31.166*XXX**3
             ELSEIF(XXX > 0.8) THEN
                FMU = -266.28+988.25*XXX-1208.33*XXX**2+486.37*XXX**3
             ENDIF

             DFA = 0.3*ADV(IACTU)*COS(ALPHA_DISC)*SIN(PSIL)*FMU* ! PSI or PSIL?
     +                 THRUST(IACTU)/ACTUA(IACTU)

             PLOSS2    = PLOSS2 + DFA

          ENDIF ! IFA(IACTU) == 3


       ELSEIF(IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19) THEN ! Blade element
          CBLA = CBLADE(IACTU)
          IF (RVECLEN/R2 >= 0.93) THEN
             CBLA = (1.0-(RVECLEN/R2-0.93)*5.14)*CBLADE(IACTU)
          ENDIF
          PLOSS2 = NBLADE(IACTU)/(2*PII*RVECLEN)*
     +         (CL*COS(alphai)-CD*sin(alphai))*
     +         0.5*ROM*VE**2*CBLA ! Force/area  


       ELSEIF(IFA(IACTU) == 6) THEN ! Miklos diploma thesis
          PLOSS2 = SURFA(NP)


       ELSEIF(IFA(IACTU) == 7 .OR. IFA(IACTU) == 8) THEN ! Uneven distribution
          PLOSS2 = SURFA(NP)*OSKU(IACTU)%DAMPT
       ENDIF                    ! IFA(IACTU) == 0

C ... Transform PLOSS2 to FINFLO coordinates

      IF (SHIP) THEN
         CALL EULFIN(0.0,A1,B1,0.0,0.0,PLOSS2,PLOSS2X,PLOSS2Y,PLOSS2Z,1)
      ELSEIF(IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19) THEN
         PLOSS2X = PLOSS2*NXBLE
         PLOSS2Y = PLOSS2*NYBLE
         PLOSS2Z = PLOSS2*NZBLE
      ELSE
         CALL EULFIN(0.0,B1,A1,0.0,0.0,-PLOSS2,
     +               PLOSS2X,PLOSS2Y,PLOSS2Z,2)
      ENDIF

          SURFX(NP) = PLOSS2X
          SURFY(NP) = PLOSS2Y
          SURFZ(NP) = PLOSS2Z

      SURF2X(NP) = PLOSS2X
      SURF2Y(NP) = PLOSS2Y
      SURF2Z(NP) = PLOSS2Z

c       write(9000+kerta,*)i,j,np,surfx(np)*a3(if1),surfy(np)*a3(if1),
c     2               surfz(np)*a3(if1)

C **********************************************************************

C ... Tangential force

C **********************************************************************

      IF(IFT(IACTU) == 0) THEN ! No tangential force

          PLOSS2M = 0.

      ELSEIF(IFT(IACTU) == 1) THEN ! An even distribution

         PLOSS2M = 1.5*TORQUE(IACTU)/(PII*(R2**3-R1**3)) !momentti

      ELSEIF(IFT(IACTU) == 2) THEN

         XXX   = (1./PLEN)/R2
         APU5  = .5*FRSDEN*SIGMA(IACTU)*CDBLADE(IACTU)*VTIP(IACTU)**2

         FTPA  = APU5*(XXX+SIN(PSIL)*COS(ALPHA_DISC)*ADV(IACTU))**2/XXX
         QP    = .5*PII*APU5*R2**3*(1.+COS(ALPHA_DISC*ADV(IACTU)**2)**2)
         QCI   = ABS(TORQUE(IACTU)) - QP
         PHI   = (VA + VIO) /
     +           (FRSVEL*COS(ALPHA_DISC)*SIN(PSIL) + VTIP(IACTU)*XXX)
 
         IF(ABS(PHI) > 0.5) PHI = 0.5 ! Limiter for PSI

         FTIA  = PHI*PLOSS2*CFTI(IACTU) ! Tarkistettava kerroin, oletus = 1.
         RLAM  = (VA + VIO)/VTIP(IACTU)
         QI    = RLAM*TORQUE(IACTU)*R2

         PLOSS2M = ABS(TORQUE(IACTU))/TORQUE(IACTU)*(FTPA + FTIA*QCI/QI)


      ELSEIF(IFT(IACTU) >= 11 .AND. IFT(IACTU) <= 19) THEN ! Blade element
         CBLA = CBLADE(IACTU)
         IF (RVECLEN/R2 >= 0.93) THEN
            CBLA = (1.0-(RVECLEN/R2-0.93)*5.14)*CBLADE(IACTU)
         ENDIF
         PLOSS2M = NBLADE(IACTU)/(2*PII*RVECLEN)*
     +        (CD*COS(alphai)+CL*sin(alphai))*
     +        0.5*ROM*VE**2*CBLA*QFACT(IACTU) ! Force/area


      ELSEIF(IFT(IACTU) == 6) THEN
         PLOSS2M = SURFT(NP)

      ELSEIF(IFT(IACTU) == 7 .OR. IFT(IACTU) == 8) THEN
         PLOSS2M = SURFT(NP)*OSKU(IACTU)%DAMPN
         IF(IFA(IACTU) /= 7 .AND. IFA(IACTU) /= 8) THEN
            WRITE(*,*)'IFA must be 7 or 8. Change IFA=',IFA(IACTU)
            STOP
         ENDIF
c        write(7850+kerta,*) i,j,np,ploss2,ploss2m
c         write(7860+kerta,*) i,j,np,surfa(np),surft(np)

      ENDIF ! IFT(IACTU) == 0


C ... Transform PLOSS2M to FINFLO coordinates

      IF (SHIP) THEN
         PLOSS2MX = PLOSS2M*RX2/RLEN
         PLOSS2MY = PLOSS2M*RY2/RLEN
         PLOSS2MZ = PLOSS2M*RZ2/RLEN
      ELSEIF(IFA(IACTU) >= 11 .AND. IFA(IACTU) <= 19) THEN
         PLOSS2MX = PLOSS2M*PXBLE
         PLOSS2MY = PLOSS2M*PYBLE
         PLOSS2MZ = PLOSS2M*PZBLE
      ELSE
         CALL EULFIN(0.0,B1,A1,
     +   -PLOSS2M*SIN(PSIL-BETA),-PLOSS2M*COS(PSIL-BETA),
     +    0.0,PLOSS2MX,PLOSS2MY,PLOSS2MZ,2)
      ENDIF

          SURFX(NP) = SURFX(NP) + PLOSS2MX
          SURFY(NP) = SURFY(NP) + PLOSS2MY
          SURFZ(NP) = SURFZ(NP) + PLOSS2MZ

      SURFMX(NP) = PLOSS2MX
      SURFMY(NP) = PLOSS2MY
      SURFMZ(NP) = PLOSS2MZ

C **********************************************************************

C ... Radial force

C **********************************************************************

      IF(IFR(IACTU) == 0) THEN ! No radial force

          PLOSS2R = 0.

      ELSEIF(IFR(IACTU) == 1) THEN ! An even distribution

         PLOSS2R = 10. ! For testing only

      ELSEIF(IFR(IACTU) == 2) THEN

         XXX     = (1./PLEN)/R2
         PLOSS2R = .5*FRSDEN*SIGMA(IACTU)*CFBLADE(IACTU)*
     2        (FRSVEL*COS(PSIL)*COS(ALPHA_DISC))**2*XXX

      ELSEIF(IFR(IACTU) >= 11 .OR. IFR(IACTU) <= 19) THEN
         CBLA = CBLADE(IACTU)
         IF (RVECLEN/R2 >= 0.93) THEN
            CBLA = (1.0-(RVECLEN/R2-0.93)*5.14)*CBLADE(IACTU)
      ENDIF
         VRAD = UM*RX+VM*RY+WM*RZ
         PLOSS2M =  NBLADE(IACTU)/(2*PII*RVECLEN)*0.5*ROM*VRAD**2*
     +        2*CBLA*CFBLADE(IACTU)

      ENDIF

C ... To FINFLO coordinates
      
      IF (SHIP) THEN
         PLOSS2RX = PLOSS2R*RX2/RLEN
         PLOSS2RY = PLOSS2R*RY2/RLEN
         PLOSS2RZ = PLOSS2R*RZ2/RLEN
      ELSEIF(IFR(IACTU) >= 11 .AND. IFR(IACTU) <= 19) THEN
         PLOSS2RX = PLOSS2M*RX1
         PLOSS2RY = -PLOSS2M*RY1
         PLOSS2RZ = -PLOSS2M*RZ1
      ELSE
         CALL EULFIN(0.0,B1,A1,
     +        PLOSS2R*COS(PSIL-BETA),PLOSS2R*SIN(PSIL-BETA),
     +        0.0,PLOSS2RX,PLOSS2RY,PLOSS2RZ,2)
      ENDIF

          SURFX(NP) = SURFX(NP) + PLOSS2RX
          SURFY(NP) = SURFY(NP) - PLOSS2RY
          SURFZ(NP) = SURFZ(NP) - PLOSS2RZ
        
      ENDDO  ; ENDDO

      IREPEA(20) = IREPEA(20) + 1

      RETURN
      END SUBROUTINE ACTUATOR_FORCES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ACTUATOR_DISC(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,UROT,POROS,
     9  RMLOSS,IF2,WMFLUX,SURLE,DSURLE,XCP,YCP,ZCP,SURFX,SURFY,SURFZ,
     1  DX,DY,DZ,CXO,CYO,CZO,CMXO,CMYO,CMZO,IACTU,CPWALL,TOM,ACNOR,
     2  SURF2X,SURF2Y,SURF2Z,SURFMX,SURFMY,SURFMZ)

      USE NS3CO, ONLY : LN
      USE CONSTANTS, ONLY : PII,EPS10

      USE FLIGHT, ONLY : XCG,YCG,ZCG,OSKU,RIN,ROUT,THRUST,TORQUE,ADV,
     &     ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI,ACTUA,IFA,IFT,UTSP

      USE BLADE_VARIABLES, ONLY    : QZE_D

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IF2,INTER,INTE2,IFN,IACTU,ACNOR

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  UROT(*),POROS(*),RMLOSS(*),WMFLUX(*),SURLE(*),DSURLE(*),
     1  SURFX(*),SURFY(*),SURFZ(*),CPWALL(*),
     2  SURF2X(*),SURF2Y(*),SURF2Z(*),SURFMX(*),SURFMY(*),SURFMZ(*)

      REAL :: UM,VM,WM,ROM,RELAX,RK1,RK2,RC1,RC2,RM1,RM2,YP1,YP4,
     1  RGAS,GAMMA,E0REF,T0REF,TURLIM,FRSSIE,PRS,PR,PRT,UC,
     2  PLOSS,UH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,PLOSS2,R1,R2,
     3  DSURFX,DSURFY,DSURFZ,DSURFE,CX,CY,CZ,BX,BY,BZ,DX,DY,DZ,DELTAP,
     4  X,Y,Z,XD,YD,ZD,CMX,CMY,CMZ,CXO,CYO,CZO,CMXO,
     5  CMYO,CMZO,TOM,APU,APU2,B2X,B2Y,B2Z,BMX,BMY,BMZ,FAX,FAY,FAZ,
     6  FMX,FMY,FMZ,TAX,TAY,TAZ,TMX,TMY,TMZ,ACTEST,FA,TQ,FAT,TQT,
     7  FAQ,TQA,UAVE

C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL :: RKK1(5)

      REAL :: FAX2,FAY2,FAZ2,TMX2,TMY2,TMZ2,FMX2,FMY2,FMZ2,
     2        TAX2,TAY2,TAZ2

      REAL :: XCP(*), YCP(*), ZCP(*)
      REAL :: XMOM,YMOM,ZMOM

      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      APU2 = 0
      INTER   = IABS(INTEM)
      UAVE  = 0

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2

      RELAX   = 1. ! Useless?

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      IF(M == 1) THEN 

      DX   = 0.; DY   = 0.; DZ   = 0.
      CXO  = 0.; CYO  = 0.; CZO  = 0.
      CMXO = 0.; CMYO = 0.; CMZO = 0.

         FAX  = 0.; FAY  = 0.; FAZ  = 0.
         FMX  = 0.; FMY  = 0.; FMZ  = 0.
         TAX  = 0.; TAY  = 0.; TAZ  = 0.
         TMX  = 0.; TMY  = 0.; TMZ  = 0.

      ENDIF

*      DO J = KY1+1,KY2-1 ! These should not be extended, obs. below

      DO J = KY1,KY2  ! Nowadays these are not extended

      IJ   = (JN+J-1)*JSTR + KA
*      NR = (JN+J-(KY1+1))*((KX2-1)-(KX1+1)+1+2*IN) - (KX1+1) + IN +1 !??
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
*      DO I = KX1+1,KX2-1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IFN = IF1 + LSTRID                  ! 2. face FLUX
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

c         IF(IBTYPE == 8) THEN

C ... Subsonic flow is assumed

         UM = .5*(U(II) + U(IM1)) - UROT(IF1)*A3X(IF1)
         VM = .5*(V(II) + V(IM1)) - UROT(IF1)*A3Y(IF1)
         WM = .5*(W(II) + W(IM1)) - UROT(IF1)*A3Z(IF1)
         ROM= .5*(RO(II)+RO(IM1))

         UC = A3X(IF1)*UM + A3Y(IF1)*VM + A3Z(IF1)*WM - UROT(IF1)
         UH = UC/(POROS(NP)+1.E-6)
         UAVE = UAVE + UH*A3(IF1)

C ... Boundary force sources

         BX      = -SURFX(NP)*IDIR *.5
         BY      = -SURFY(NP)*IDIR *.5
         BZ      = -SURFZ(NP)*IDIR *.5

C ... Axial and tangential boundary force sources

         B2X     = -SURF2X(NP)*IDIR *.5
         B2Y     = -SURF2Y(NP)*IDIR *.5
         B2Z     = -SURF2Z(NP)*IDIR *.5

         BMX     = -SURFMX(NP)*IDIR *.5
         BMY     = -SURFMY(NP)*IDIR *.5
         BMZ     = -SURFMZ(NP)*IDIR *.5

C ... Energy

c         DSURFE   =  UM*BX + VM*BY + WM*BZ
         DSURFE   =  U(II)*BX + V(II)*BY + W(II)*BZ

         F3RM(IF1) = F3RM(IF1) + A3(IF1)*BX     ! A3(IF1)*A3X(IF1)*PLOSS
         F3RN(IF1) = F3RN(IF1) + A3(IF1)*BY     ! A3(IF1)*A3Y(IF1)*PLOSS
         F3RW(IF1) = F3RW(IF1) + A3(IF1)*BZ     ! A3(IF1)*A3Z(IF1)*PLOSS
         F3E(IF1)  = F3E(IF1)  + A3(IF1)*DSURFE

C ... Next surface (kommentoitu 1. kertaluvun menetelmää varten

c         F3RM(IFN)= F3RM(IFN) + R2*A3(IFN)*BX*2.
c         F3RN(IFN)= F3RN(IFN) + R2*A3(IFN)*BY*2.
c         F3RW(IFN)= F3RW(IFN) + R2*A3(IFN)*BZ*2.
c         F3E(IFN) = F3E(IFN)  + R2*A3(IFN)*DSURFE*2.

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Add turbulence source ?
c         F3RK(IF1)  = F3R(IF1) + ??? 
c         F3EPS(IF1) = F3R(IF1) + ???
         ENDIF

      IF(M == 1) THEN ! Store forces and moments

C ... T x OMEGA from the energy flux

         TOM     = TOM - IDIR*(F3E(IF1) - F3RK(IF1)) ! EFFECT OF TURBULENCE

C ... Pressure forces

         DELTAP  = 0.
         DELTAP  = DELTAP + PDIFF(II) - PDIFF(IM1)
         DX      = DX  - IDIR*(A3(IF1)*A3X(IF1) * DELTAP) ! Oli -IDIR*
         DY      = DY  - IDIR*(A3(IF1)*A3Y(IF1) * DELTAP)
         DZ      = DZ  - IDIR*(A3(IF1)*A3Z(IF1) * DELTAP)

C ... Vector for the momentum

         X       = XCP(NP)
         Y       = YCP(NP)
         Z       = ZCP(NP)
         XD      = X - XCG(IACTU)
         YD      = Y - YCG(IACTU)
         ZD      = Z - ZCG(IACTU)

C ... Total force from the momentum flux

         CX      = CX - F3RM(IF1)*IDIR 
!!! MCHECKiin menee pinnan yli integroidut vuot. Niiden erotus
!!! pitäisi olla kokonaisvoima. Siis haetaan sillä indeksisuunnalla
!!! kyseisten pintojen vuo
!!! Korjaus: siellä on kokonaisliikemäärä. Mukana oleva rotdia, muuttaa
!!! Siihen y-liikemäärän.
         CY      = CY - F3RN(IF1)*IDIR
         CZ      = CZ - F3RW(IF1)*IDIR
         
C ... Total torque from the momentum flux

         CMX     = CMX + F3RW(IF1)*IDIR*YD - F3RN(IF1)*IDIR*ZD
         CMY     = CMY + F3RM(IF1)*IDIR*ZD - F3RW(IF1)*IDIR*XD
         CMZ     = CMZ + F3RN(IF1)*IDIR*XD - F3RM(IF1)*IDIR*YD

c         CMX     = CMX - F3RW(IF1)*IDIR*YD + F3RN(IF1)*IDIR*ZD ! Antakee mersu
c         CMY     = CMY - F3RM(IF1)*IDIR*ZD + F3RW(IF1)*IDIR*XD
c         CMZ     = CMZ - F3RN(IF1)*IDIR*XD + F3RM(IF1)*IDIR*YD

C ... Added actuator forces

         BX      = A3(IF1)*BX
         BY      = A3(IF1)*BY
         BZ      = A3(IF1)*BZ
         CXO     = CXO  - BX
         CYO     = CYO  - BY
         CZO     = CZO  - BZ
         
C ... Added axial and tangential actuator forces

         B2X     = A3(IF1)*B2X
         B2Y     = A3(IF1)*B2Y
         B2Z     = A3(IF1)*B2Z
         FAX     = FAX  - B2X
         FAY     = FAY  - B2Y
         FAZ     = FAZ  - B2Z

         BMX     = A3(IF1)*BMX
         BMY     = A3(IF1)*BMY
         BMZ     = A3(IF1)*BMZ
         FMX     = FMX  - BMX
         FMY     = FMY  - BMY
         FMZ     = FMZ  - BMZ

C ... Added actuator torque

         CMXO    = CMXO + BZ*YD - BY*ZD
         CMYO    = CMYO + BX*ZD - BZ*XD
         CMZO    = CMZO + BY*XD - BX*YD

C ... Added axial and tangential actuator torques

         TAX     = TAX + B2Z*YD - B2Y*ZD
         TAY     = TAY + B2X*ZD - B2Z*XD
         TAZ     = TAZ + B2Y*XD - B2X*YD

         TMX     = TMX + BMZ*YD - BMY*ZD
         TMY     = TMY + BMX*ZD - BMZ*XD
         TMZ     = TMZ + BMY*XD - BMX*YD

C ... Store data for visualization

         SURLE(NP)  = DSURFE
C         DSURLE(NP) = BX
         WMFLUX(NP) = ROM*UC

         DSURLE(NP) = UH
         CPWALL(NP) = MAX(-.99*FRSPRE, PDIFF(II))
          
      ENDIF ! M == 1
          
      ENDDO
      ENDDO

      UAVE = ABS(UAVE)/ACTUA(IACTU)
      IF(IFA(IACTU) == 8) UTSP(IACTU) = UAVE
      IF(IFA(IACTU) == 8) QZE_D(IACTU) = UTSP(IACTU)     

      RETURN
      END SUBROUTINE ACTUATOR_DISC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ACTUATOR_SCALE(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,UROT,POROS,
     9  RMLOSS,IF2,WMFLUX,SURLE,DSURLE,XCP,YCP,ZCP,SURFX,SURFY,SURFZ,
     1  DX,DY,DZ,CXO,CYO,CZO,CMXO,CMYO,CMZO,IACTU,CPWALL,TOM,ACNOR,
     2  SURF2X,SURF2Y,SURF2Z,SURFMX,SURFMY,SURFMZ)

      USE NS3CO, ONLY : LN
      USE CONSTANTS, ONLY : PII,EPS10

      USE FLIGHT, ONLY    : XCG,YCG,ZCG,OSKU,RIN,ROUT,THRUST,TORQUE,ADV,
     &     ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI,ACTUA,IFA,IFT,UTSP
      USE BLADE_VARIABLES, ONLY : QZE_S1,QBE_S1

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IF2,INTER,INTE2,IFN,IACTU,ACNOR

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  UROT(*),POROS(*),RMLOSS(*),WMFLUX(*),SURLE(*),DSURLE(*),
     1  SURFX(*),SURFY(*),SURFZ(*),CPWALL(*),
     2  SURF2X(*),SURF2Y(*),SURF2Z(*),SURFMX(*),SURFMY(*),SURFMZ(*)

      REAL :: UM,VM,WM,ROM,RELAX,RK1,RK2,RC1,RC2,RM1,RM2,YP1,YP4,
     1  RGAS,GAMMA,E0REF,T0REF,TURLIM,FRSSIE,PRS,PR,PRT,UC,
     2  PLOSS,UH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,PLOSS2,R1,R2,
     3  DSURFX,DSURFY,DSURFZ,DSURFE,CX,CY,CZ,BX,BY,BZ,DX,DY,DZ,DELTAP,
     4  X,Y,Z,XD,YD,ZD,CMX,CMY,CMZ,CXO,CYO,CZO,CMXO,
     5  CMYO,CMZO,TOM,APU,APU2,B2X,B2Y,B2Z,BMX,BMY,BMZ,FAX,FAY,FAZ,
     6  FMX,FMY,FMZ,TAX,TAY,TAZ,TMX,TMY,TMZ,ACTEST,FA,TQ,FAT,TQT,
     7  FAQ,TQA,UAVE,F3RMIF,F3RNIF,F3RWIF,F3EIF,DAMPTO,DAMPNO

C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL :: RKK1(5)

      REAL ::  FAX2,FAY2,FAZ2,TMX2,TMY2,TMZ2,FMX2,FMY2,FMZ2,
     2         TAX2,TAY2,TAZ2

      REAL :: XCP(*), YCP(*), ZCP(*)
      REAL :: XMOM,YMOM,ZMOM

      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      APU2 = 0
      INTER   = IABS(INTEM)
      UAVE  = 0

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1 = RKK1(INTE2)
      RK2 = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2

      RELAX   = 1. ! Useless?

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      DX   = 0.; DY   = 0.; DZ   = 0.
      CXO  = 0.; CYO  = 0.; CZO  = 0.
      CMXO = 0.; CMYO = 0.; CMZO = 0.

      FAX  = 0.; FAY  = 0.; FAZ  = 0.
      FMX  = 0.; FMY  = 0.; FMZ  = 0.
      TAX  = 0.; TAY  = 0.; TAZ  = 0.
      TMX  = 0.; TMY  = 0.; TMZ  = 0.

*      DO J = KY1+1,KY2-1 ! These should not be extended, obs. below

      DO J = KY1,KY2  ! Nowafays these are not extended

      IJ   = (JN+J-1)*JSTR + KA
*      NR = (JN+J-(KY1+1))*((KX2-1)-(KX1+1)+1+2*IN) - (KX1+1) + IN +1 !??
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
*      DO I = KX1+1,KX2-1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IFN = IF1 + LSTRID                  ! 2. face FLUX
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

c         IF(IBTYPE == 8) THEN

C ... Subsonic flow is assumed

         UM = .5*(U(II) + U(IM1)) - UROT(IF1)*A3X(IF1)
         VM = .5*(V(II) + V(IM1)) - UROT(IF1)*A3Y(IF1)
         WM = .5*(W(II) + W(IM1)) - UROT(IF1)*A3Z(IF1)
         ROM= .5*(RO(II)+RO(IM1))

         UC = A3X(IF1)*UM + A3Y(IF1)*VM + A3Z(IF1)*WM - UROT(IF1)
         UH = UC/(POROS(NP)+1.E-6)
         UAVE = UAVE + UH*A3(IF1)

C ... Boundary force sources

         BX      = -SURFX(NP)*IDIR *.5
         BY      = -SURFY(NP)*IDIR *.5
         BZ      = -SURFZ(NP)*IDIR *.5

C ... Axial and tangential boundary force sources

         B2X     = -SURF2X(NP)*IDIR *.5
         B2Y     = -SURF2Y(NP)*IDIR *.5
         B2Z     = -SURF2Z(NP)*IDIR *.5

         BMX     = -SURFMX(NP)*IDIR *.5
         BMY     = -SURFMY(NP)*IDIR *.5
         BMZ     = -SURFMZ(NP)*IDIR *.5

C ... Vector for the momentum

         X       = XCP(NP)
         Y       = YCP(NP)
         Z       = ZCP(NP)
         XD      = X - XCG(IACTU)
         YD      = Y - YCG(IACTU)
         ZD      = Z - ZCG(IACTU)

C ... Added actuator forces

         BX      = A3(IF1)*BX
         BY      = A3(IF1)*BY
         BZ      = A3(IF1)*BZ

C ... Added axial and tangential actuator forces

         B2X     = A3(IF1)*B2X
         B2Y     = A3(IF1)*B2Y
         B2Z     = A3(IF1)*B2Z
         FAX     = FAX  - B2X
         FAY     = FAY  - B2Y
         FAZ     = FAZ  - B2Z

         BMX     = A3(IF1)*BMX
         BMY     = A3(IF1)*BMY
         BMZ     = A3(IF1)*BMZ
         FMX     = FMX  - BMX
         FMY     = FMY  - BMY
         FMZ     = FMZ  - BMZ

C ... Added axial and tangential actuator torques

         TAX     = TAX + B2Z*YD - B2Y*ZD
         TAY     = TAY + B2X*ZD - B2Z*XD
         TAZ     = TAZ + B2Y*XD - B2X*YD
c       write(9200,*) i,j,np,b2x,b2y,b2z,TAX

         TMX     = TMX + BMZ*YD - BMY*ZD
         TMY     = TMY + BMX*ZD - BMZ*XD
         TMZ     = TMZ + BMY*XD - BMX*YD

c      ENDIF ! M == 1 ei voi rajoittua tahan
          
      ENDDO
      ENDDO

      UAVE = ABS(UAVE)/ACTUA(IACTU)
      IF(IFA(IACTU) == 8) UTSP(IACTU) = UAVE

1000  CONTINUE

C ... Normalization coefficient (kakkonen pitais korvata jotenkin...)     

      ACTEST = ABS(OSKU(IACTU)%DAMPT-1.)

      IF(ACNOR == 1 .AND. ACTEST < 1.E-5  .OR.
     +   ACNOR == 1 .AND. IFA(IACTU) == 8) THEN

c      write(6100,'(7A15)') 'IACTU','FAT','TQT','FA','TQ','FAQ','TQA'

         FAT = ABS(THRUST(IACTU))
         TQT = ABS(TORQUE(IACTU))

         CALL FINEUL(0.0,ROTA1(IACTU),ROTB1(IACTU),
     +        FAX,FAY,FAZ,FAZ2,FAY2,FAX2,1)

         CALL FINEUL(0.0,ROTA1(IACTU),ROTB1(IACTU),
     +        TMX,TMY,TMZ,TMZ2,TMY2,TMX2,1)

         CALL FINEUL(0.0,ROTA1(IACTU),ROTB1(IACTU),
     +        FMX,FMY,FMZ,FMZ2,FMY2,FMX2,1)

         CALL FINEUL(0.0,ROTA1(IACTU),ROTB1(IACTU),
     +        TAX,TAY,TAZ,TAZ2,TAY2,TAX2,1)

         FA  = FAX2
         TQ  = TMX2  ! Torque
         FAQ = FMX2  ! Axial force caused byt the torque
         TQA = TAX2  ! Torque caused byt the axial force

         DAMPTO = OSKU(IACTU)%DAMPT
         DAMPNO = OSKU(IACTU)%DAMPN
         OSKU(IACTU)%DAMPT=(FAT-FAQ*TQT/(TQ+EPS10))/
     &                     (FA-FAQ*TQA/(TQ+EPS10))*.5*DAMPTO
         OSKU(IACTU)%DAMPN=(TQT-FAT*TQA/(FA+EPS10))/
     &                     (TQ-FAQ*TQA/(FA+EPS10))*.5*DAMPNO

         IF (IFA(IACTU) == 8 .OR. IFA(IACTU) == 10) THEN
            QZE_S1(IACTU) = OSKU(IACTU)%DAMPT
            QBE_S1(IACTU) = OSKU(IACTU)%DAMPN
         ENDIF
c         write(*,*)OSKU(IACTU)%DAMPT,OSKU(IACTU)%DAMPN,IACTU
      ENDIF
c      write(6100,'(I5,10X,6F15.7)')IACTU,FAT,TQT,FA,TQ,FAQ,TQA
c      ENDIF ! M == 1

      RETURN
      END SUBROUTINE ACTUATOR_SCALE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MINOR_LOSS(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,UROT,POROS,
     9  RMLOSS,IF2,WMFLUX,SURLE,DSURLE,XCP,YCP,ZCP,SURFX,SURFY,SURFZ,
     1  CX,CY,CZ,CPWALL,SURFPX,SURFPY,SURFPZ)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IF2,INTER,INTE2,IFN

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  UROT(*),POROS(*),RMLOSS(*),WMFLUX(*),SURLE(*),DSURLE(*),
     1  SURFX(*),SURFY(*),SURFZ(*),CPWALL(*),
     2  SURFPX(*),SURFPY(*),SURFPZ(*)

      REAL :: UM,VM,WM,ROM,RELAX,RK1,RK2,RC1,RC2,RM1,RM2,YP1,YP4,
     1  RGAS,GAMMA,E0REF,T0REF,TURLIM,FRSSIE,PRS,PR,PRT,UC,
     2  PLOSS,UH,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,PLOSS2,R1,R2,F,
     3  CX,CY,CZ

C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES USING MUSCL
C ... RK1=1  RK2=1 (H{NEL),    RK1=0   RK2=2   UPWIND,
C ... RK1=2  RK2=0 CENTRAL,    RK1=4/3 RK2=2/3 THIRD ORDER UPWIND BIASED
C ... RK1=0  RK2=0 FIRST ORDER RK1=3/2 RK2=1/2 QUICK

      REAL :: RKK1(5)

      REAL :: XCP(*), YCP(*), ZCP(*)

      DATA RKK1/0.,1.,1.3333333,2.,1.5/

      INTER   = IABS(INTEM)

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      INTE2   = MIN(5,MAX(1,INTER))
      RK1     = RKK1(INTE2)
      RK2     = 2.-RK1

      R1      = .25*RK1
      R2      = .25*RK2

      RELAX   = 1. ! Useless?
      F       = 0.
c      CX = 0.; CY = 0.; CZ = 0. Nollataan ns3c.f

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR
c      write(6,*) 'minor',kx1,kx2,ky1,ky2

*      DO J = KY1+1,KY2-1 ! These should not be extended, obs. below

      DO J = KY1,KY2  ! Nowadays these are not extended

      IJ   = (JN+J-1)*JSTR + KA
*      NR = (JN+J-(KY1+1))*((KX2-1)-(KX1+1)+1+2*IN) - (KX1+1) + IN +1 !??
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
*      DO I = KX1+1,KX2-1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IFN = IF1 + LSTRID                  ! 2. face FLUX
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

c         IF(IBTYPE == 7) THEN

C ... Subsonic flow is assumed

          UM = .5*(U(II) + U(IM1)) - UROT(IF1)*A3X(IF1)
          VM = .5*(V(II) + V(IM1)) - UROT(IF1)*A3Y(IF1)
          WM = .5*(W(II) + W(IM1)) - UROT(IF1)*A3Z(IF1)
          ROM= .5*(RO(II)+RO(IM1))

          UC = A3X(IF1)*UM + A3Y(IF1)*VM + A3Z(IF1)*WM - UROT(IF1)
          UH = UC/(POROS(NP)+1.E-6)

C.... Momentum source

          PLOSS2 = -IDIR * .5*RMLOSS(NP)*ROM*UH*ABS(UH)
          PLOSS  = .5*PLOSS2
          PLOSS2 = R2*PLOSS2   *.5

C ... Boundary surface

         F3RM(IF1) = F3RM(IF1) + A3(IF1)*A3X(IF1)*PLOSS 
c     &                       + A3(IFN)*A3X(IFN)*PLOSS2
         F3RN(IF1) = F3RN(IF1) + A3(IF1)*A3Y(IF1)*PLOSS
c     &                       + A3(IFN)*A3Y(IFN)*PLOSS2

         F3RW(IF1) = F3RW(IF1) + A3(IF1)*A3Z(IF1)*PLOSS
c     &                       + A3(IFN)*A3Z(IFN)*PLOSS2

         F3E(IF1)  = F3E(IF1)  + A3(IF1)*UC*PLOSS

C ... Next surface (mersu)

c         F3RM(IFN)= F3RM(IFN) + A3(IFN)*A3X(IFN)*PLOSS2
c         F3RN(IFN)= F3RN(IFN) + A3(IFN)*A3Y(IFN)*PLOSS2
c         F3RW(IFN)= F3RW(IFN) + A3(IFN)*A3Z(IFN)*PLOSS2
c         F3E(IFN) = F3E(IFN)  + A3(IFN)*UC*PLOSS2

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Add turbulence source ?
c         F3RK(IF1)  = F3R(IF1) + ??? 
c         F3EPS(IF1) = F3R(IF1) + ???
         ENDIF

C ... Store data for visualization

        IF(M == 1) THEN

          SURLE(NP)  = PLOSS
          WMFLUX(NP) = ROM*UC
          DSURLE(NP) = UH
          CPWALL(NP) = MAX(-.99*FRSPRE, PDIFF(II))

          SURFX(NP) = -A3X(IF1)*PLOSS * IDIR
          SURFY(NP) = -A3Y(IF1)*PLOSS * IDIR
          SURFZ(NP) = -A3Z(IF1)*PLOSS * IDIR

          CX = CX + A3(IF1)*SURFX(NP)!*IDIR
          CY = CY + A3(IF1)*SURFY(NP)!*IDIR
          CZ = CZ + A3(IF1)*SURFZ(NP)!*IDIR
c            write(701,*) cx,cy,cz
        ENDIF

c         write(674,*) SURFX(NP)
c         write(674,*) SURFY(NP)
c         write(674,*) SURFZ(NP)

      ENDDO
      ENDDO

1000  CONTINUE

      RETURN
      END SUBROUTINE MINOR_LOSS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INFLO(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,HAT1)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NN,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  HAT1(*)

      REAL :: PR,PRT,PRS,FRSPRE,FRSDEN,FRSVEL,FRSVIS,T0,DIFPRE,
     2  VELOC,DPKOE,RELAX,DPDIF,ROAVE,RGAS,GAMMA,E0REF,T0REF,
     3  RK1,RK2,TURLIM,RC1,RC2,RM1,RM2,F2V,F3V,F4V,F5V,F6V,YP1,YP4,
     4  PB1,PB4,PB,YPD2,SURG,SURF,SURFE,DUC,FRSSIE

      RELAX   = 1.

      RK1     = 1.
      RK2     = 1. - RK1
      RC1     = 1.     !500
      RC2     = RC1 - 1.
      IF(INTEM == 6) THEN
         RC1  = 1.
         RC2  = RC1 - 1.
      ENDIF
C ... For turbulence variables
      RM1     = 1.000
      RM2     = RM1 - 1.
C     k-epsilon or SST
      
      YP1     = 1.
      YP4     = 0.
                   
      IF(IFLUX == 4 .OR. IFLUX == 6 .OR. IFLUX == 7) THEN ! Pressure difference
         YP1  = 0.
         YP4  = 1.
      ENDIF
      PRS     = PR/PRT

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

         IF(IBTYPE == 7) THEN

         IF(M2 == 1) THEN ! Update ghost cells
C ... Subsonic flow is assumed, extrapolation of pressure
         DPKOE       = RELAX*(RK1*P(II)     + RK2*P(II2)) +
     +                (1.- RELAX)*P(IM1) 
           
         DPDIF       = RELAX*(RK1*PDIFF(II) + RK2*PDIFF(II2)) +
     +                (1.- RELAX)*PDIFF(IM1) 
         EPS2(IM1)   = RK1*EPS2(II) + RK2*EPS2(II2)
         EPS2(IM2)   = EPS2(IM1)
         VIST(IM1)   = RK1*VIST(II) + RK2*VIST(II2)
         VIST(IM2)   = VIST(IM1)
         DPKOE       = MAX(DPKOE,.01*FRSPRE) 

         RM(IM1)     = RO(NP)*BOUNU(NP)
         RN(IM1)     = RO(NP)*BOUNV(NP)
         RW(IM1)     = RO(NP)*BOUNW(NP)

         TEMP(IM1)   = BOUNT(NP)
         TEMP(IM2)   = TEMP(IM1)

         RM(IM2)     = RM(IM1)
         RN(IM2)     = RN(IM1)
         RW(IM2)     = RW(IM1)

         P(IM1)      = DPKOE
         P(IM2)      = P(IM1)
         PDIFF(IM1)  = DPDIF
         PDIFF(IM2)  = PDIFF(IM1)
         
         CALL EFPT(E(IM1),P(IM1),TEMP(IM1),1,ISTATE,RGAS,GAMMA,0.,
     +   FRSSIE,E0REF,T0REF)
         
         RO(IM2) = RO(IM1)
         E(IM2)  = E(IM1)
         U(IM1)  = RM(IM1)/RO(IM1)
         U(IM2)  = RM(IM2)/RO(IM2)
         V(IM1)  = RN(IM1)/RO(IM1)
         V(IM2)  = RN(IM2)/RO(IM2)
         W(IM1)  = RW(IM1)/RO(IM1)
         W(IM2)  = RW(IM2)/RO(IM2)
         E(IM1)  = RO(IM1)*(E(IM1)+.5*(U(IM1)**2+V(IM1)**2+W(IM1)**2))
         E(IM2)  = RO(IM2)*(E(IM2)+.5*(U(IM2)**2+V(IM2)**2+W(IM2)**2))

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         RK(IM1)     = BOUNRK(NP)
         RK(IM2)     = RK(IM1)
         REPS(IM1)   = BOUNEP(NP)
         REPS(IM2)   = REPS(IM1)
         E(IM1)      = E(IM1) + RK(IM1)
         E(IM2)      = E(IM2) + RK(IM2)
         IF(ITURB == 9) THEN
           RNUT(IM1) = BOUNEP(NP)/RO(IM1) ; RNUT(IM2) = RNUT(IM1)
         ENDIF
      ENDIF

      IF (NSCAL >= 1) THEN
      DO  560 NS = 1,NSCAL
         FI(IM1,NS)  = BOUNFI(NP,NS)
         FI(IM2,NS)  = FI(IM1,NS)
560   CONTINUE
      ENDIF

      IF (ISTRES >= 1) THEN
      DO  570 NS = 1,5
         BIJ(IM1,NS) = BOUNBI(NP,NS)
         BIJ(IM2,NS) = BIJ(IM1,NS) 
570   CONTINUE
      ENDIF

      ENDIF !(M2 == 1)  Ghost cells are updated

         YPD2    = A3(IF1)/(VOL(II) + VOL(IM1))
         SURG    = -A3(IF1)*YPD2
         SURF    = SURG*(VIS(II) + VIS(IM1))             ! LAMINAR
         SURFE    = SURG*((1.+PRS*VIST(II)/VIS(II))*CH(II) + 
     +   (1.+PRS*VIST(IM1)/VIS(IM1))* CH(IM1))        ! BOUSSINESQ
         DUC     = .333333*(A3X(IF1)*(U(II) - U(IM1))   +
     +           A3Y(IF1)*(V(II) - V(IM1)) + A3Z(IF1)*(W(II) - W(IM1)))

         F2V     = SURF*(U(II)-U(IM1) + A3X(IF1)*DUC)
         F3V     = SURF*(V(II)-V(IM1) + A3Y(IF1)*DUC)
         F4V     = SURF*(W(II)-W(IM1) + A3Z(IF1)*DUC)
         F5V     = .5*((U(II)+U(IM1))*F2V+(V(II)+V(IM1))*F3V+
     +           (W(II)+W(IM1))*F4V) + SURFE*(TEMP(II) - TEMP(IM1))

         PB1     = MAX(0.01*FRSPRE, RC1*P(II) - RC2*P(II2))
         PB4     = RC1*PDIFF(II) - RC2*PDIFF(II2)
         PB      = YP1*PB1 + YP4*PB4

         F3R(IF1) = A3(IF1)*BOUNMF(NP)
         HAT1(IF1)= F3R(IF1)
         F3RM(IF1)= F3R(IF1)*U(IM1) + A3(IF1)*A3X(IF1)*PB + F2V
         F3RN(IF1)= F3R(IF1)*V(IM1) + A3(IF1)*A3Y(IF1)*PB + F3V
         F3RW(IF1)= F3R(IF1)*W(IM1) + A3(IF1)*A3Z(IF1)*PB + F4V
         F3E(IF1) = F3R(IF1)*(E(IM1)+ PB + .667*RK(IM1))/RO(IM1) +
     +             A3(IF1)*WROT(IF1)*(PB + .667*RK(IM1) /RO(IM1))+ F5V
         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         F3RK(IF1)  = F3R(IF1)*RK(IM1)  /RO(IM1)
         F3EPS(IF1) = F3R(IF1)*REPS(IM1)/RO(IM1)
         ENDIF
      ENDIF ! IBTYPE == 7
      ENDDO
      ENDDO

1000  CONTINUE

      RETURN
      END SUBROUTINE INFLO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTFLO
      END SUBROUTINE OUTFLO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INFPER(F3R,F3RM,F3RN,F3RW,F3E,F3RK,F3EPS,
     1  P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,OHMI,S11,
     2  VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,KY1,
     2  KY2,
     3  RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,PR,PRT,IFLUX,M,M2,IN,JN,KN,
     4  FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,
     5  ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,IBTYPE,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,BOUNW,
     7  BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNMF,
     8  MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,FRSSIE,RNUT,ZZZ,RKSI,
     9  D2,UROT,HAT1,HAT2,HAT3,HAT4,XXTRAL,INCHIML,
     9  MAW,MAXW,IMAX,JMAX,KMAX,NGL,ITYP,INLRC,INLOUT,ARTSSP,PSEUCO,
     1  FRSTEM,IUPPTL)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,IBTYPE,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NN,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  MAW,MAXW,IMAX,JMAX,KMAX,IJTRID,IL,KL,IA,IG,I8,KX11,KX22,KY11,
     6  KY22,NGL,ITYP,INLRC,INLOUT,IUPPTL

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),F3R(*),F3RM(*),
     8  F3RN(*),F3RW(*),F3E(*),F3RK(*),F3EPS(*),BOUNMF(*),RNUT(*),
     9  RKSI(*),D2(*),UROT(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*)

      REAL :: PR,PRT,PRS,FRSPRE,FRSDEN,FRSVEL,FRSVIS,T0,DIFPRE,
     2  VELOC,DPKOE,RELAX,DPDIF,ROAVE,RGAS,GAMMA,E0REF,T0REF,
     3  RK1,RK2,TURLIM,RC1,RC2,RM1,RM2,F2V,F3V,F4V,F5V,F6V,YP1,YP4,
     4  PB1,PB4,PB,YPD2,SURG,SURF,SURFE,DUC,FRSSIE,RMIPI,RNIPI,RWIPI,
     5  RMIMI,RNIMI,RWIMI,A11,A12,A13,A21,A22,A23,A31,A32,A33,UR,UL,
     6  VR,VL,WR,WL,RKR,RKL,RER,REL,EIR,EIL,PDR,PDL,UC,F1PIG,PKR,PKL,
     7  ROHAT,PAVE1,UAVE,UAVE1,F1PI,F2PIG,F3PIG,F4PIG,F5PIG,F6PIG,
     8  F7PIG,F2PM,F3PM,F4PM,RMIPP,RMIMM,PAVE,AREA,UDAMP,ARTSSP,PSEUCO,
     9  FRSTEM

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:), ERL(:),ELR(:),DEDPH(:),DEDRH(:),
     +     REIM(:),REIP(:),PDMI(:),PDPI(:),TEIM(:),TEIP(:)
      REAL, TARGET ::  ZZZ(MAXW)

      LOGICAL :: XXTRAL, INCHIML, INCH

c      TYPE(PRE_COR) PRC(*)

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      PDMI => ZZZ(14*MAW+1:15*MAW);PDPI => ZZZ(15*MAW+1:16*MAW)
      RKIM => ZZZ(16*MAW+1:17*MAW);RKIP => ZZZ(17*MAW+1:18*MAW) 
      REIM => ZZZ(18*MAW+1:19*MAW);REIP => ZZZ(19*MAW+1:20*MAW)
      TEIM => ZZZ(22*MAW+1:23*MAW);TEIP => ZZZ(23*MAW+1:24*MAW) 
C
C ... FLUXES FOR PERIODIC PATCHES. SIMPLIFIED SPLITTING
C ... THIS SUBROUTINE UTILIZES SPECIFIC INTERNAL ENERGY
C

      AREA    = 0.

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN

      KX11    = KX1 !+ 1
      KX22    = KX2 !- 1
      KY11    = KY1 !+ 1
      KY22    = KY2 !- 1

      IJTRID  = (KX22-KX11+1)*(KY22-KY11+1)
      IL      = KSTR
      KL      = 2*IL

      LSTRID  = KSTR    ! KSTR*IDIR ! No BOT/TOP treatment
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = LSTRID
      KA      = (KN+KBOT-1)*KSTR !+ NFL

c      write(6666,*) ityp,inlrc,ijtrid,lstrid,il,kl,idir,ka,xxtral
      IF (ITURB  < 3 .OR. ITURB == 8) THEN ! ALGEBRAIC MODEL, ZERO K AND EPSILON
         CALL ZEROZZ(RKIP,ISTRID*JSTRID)
         CALL ZEROZZ(RKIM,ISTRID*JSTRID)
         CALL ZEROZZ(REIP,ISTRID*JSTRID)
         CALL ZEROZZ(REIM,ISTRID*JSTRID)
      ENDIF

C ... EXTRAPOLATION OF THE DEPENDENT (PRIMITIVE) VARIABLES

      IF(.NOT. XXTRAL) THEN

      INCH = INCHIML ! Mersu

      CALL XXTRAP_SURF(TEMP,TEIP,TEIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)
      CALL XXTRAP_SURF(   U,RMIP,RMIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)
c      CALL XXTRAP_SURF2(   U,RMIP,RMIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
c     +      LSTRID,KA,INTEM,RKSI,INCH,NFL,ngl)
      CALL XXTRAP_SURF(   V,RNIP,RNIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)
      CALL XXTRAP_SURF(   W,RWIP,RWIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)
      CALL XXTRAP_SURF(  P,PPI,PMI,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)
      CALL XXTRAP_SURF(PDIFF,PDPI,PDMI,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,RKSI,INCH,NFL)

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAP_SURF(  RK,RKIP,RKIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTET,RKSI,INCH,NFL)
      CALL XXTRAP_SURF(REPS,REIP,REIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTET,RKSI,INCH,NFL)
      ENDIF

      ELSE

      CALL XXTRAA_SURF(TEMP,TEIP,TEIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(   U,RMIP,RMIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(   V,RNIP,RNIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(   W,RWIP,RWIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(  P,PPI,PMI,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(PDIFF,PDPI,PDMI,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTEM,A3,VOL,D2,RKSI,INCH,NFL)

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL XXTRAA_SURF(  RK,RKIP,RKIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTET,A3,VOL,D2,RKSI,INCH,NFL)
      CALL XXTRAA_SURF(REPS,REIP,REIM,KX1,KX2,KY1,KY2,ISTR,JSTR,IN,JN,
     +      LSTRID,KA,INTET,A3,VOL,D2,RKSI,INCH,NFL)
      ENDIF

      ENDIF ! .NOT.XXTRAL
 
      DO IG  = 1,IJTRID
         IF(ISTATE /= 10) THEN ! Limit the pressure
         PPI(IG)  = MAX(0.001*FRSPRE,PPI(IG))
         PMI(IG)  = MAX(0.001*FRSPRE,PMI(IG))
         PDPI(IG) = MAX(-0.999*FRSPRE,PDPI(IG))
         PDMI(IG) = MAX(-0.999*FRSPRE,PDMI(IG))
         ENDIF
         RKIP(IG) = MAX(0.,RKIP(IG))
         RKIM(IG) = MAX(0.,RKIM(IG))
         REIP(IG) = MAX(0.,REIP(IG))
         REIM(IG) = MAX(0.,REIM(IG))
      ENDDO

C ... DENSITY AND SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE
C     AND TEMPERATURE. ISTATE Determines the Equation of State.

      CALL ROFPT(TEIP,PPI,ROIP,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL ROFPT(TEIM,PMI,ROIM,IJTRID,ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,
     +           E0REF,T0REF,FRSTEM,IUPPTL)
      CALL EFPT(EIP,PPI,TEIP,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)
      CALL EFPT(EIM,PMI,TEIM,IJTRID,ISTATE,RGAS,GAMMA,0.,FRSSIE,
     +           E0REF,T0REF)
       
C ... Indeces from patch data

      LSTRID  = KSTR   !*IDIR ! No BOT/TOP treatment
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = LSTRID
      KA      = (KN+KBOT-1)*KSTR

      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1   = (JN+J-1)*ISTRID + IN
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II      = II + NFL
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IG      = I  + NN                ! Patch interior index
         NP  = I   + NR                      ! Patch index with LN ghost cells

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A3X(II)
         A12      = A3Y(II)
         A13      = A3Z(II)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         UR       = A11*RMIPI + A12*RNIPI + A13*RWIPI
         VR       = A21*RMIPI + A22*RNIPI + A23*RWIPI
         WR       = A31*RMIPI + A32*RNIPI + A33*RWIPI
         UL       = A11*RMIMI + A12*RNIMI + A13*RWIMI
         VL       = A21*RMIMI + A22*RNIMI + A23*RWIMI
         WL       = A31*RMIMI + A32*RNIMI + A33*RWIMI

C*******************************************************************
C ... SIMPLIFIED ROE BASED FLUX SPLITTING                          *
C*******************************************************************

c      YPROR   = 1./ROIP(IG)
c      YPROL   = 1./ROIM(IG)
      RKR     = RKIP(IG) !*YPROR
      RER     = REIP(IG) !*YPROR
      RKL     = RKIM(IG) !*YPROL
      REL     = REIM(IG) !*YPROL

C ... TOTAL INTERNAL ENERGY

c      PDR     = PDI(II) !  PDPI(IG)
      PDR     = PDPI(IG)
      PKR     = .6667*ROIP(IG)*RKIP(IG)

c      PDL     = PDI(I-IL) ! PDMI(IG)
      PDL     = PDMI(IG)
      PKL     = .6667*ROIM(IG)*RKIM(IG)
      ROHAT   =.5*(ROIP(IG)+ROIM(IG))

      UC      = UROT(II)
 
* Equations replaced by ESa 3.9.2001 (equations from 'propeller-finflo' 

      RMIPP   = A11*U(II)    + A12*V(II)    + A13*W(II)

      RMIMM   = A11*U(II-IL) + A12*V(II-IL) + A13*W(II-IL)
      UAVE1   = (RMIPP + RMIMM)*.5 - UC
c     No damping term in the velocity, for testing only
      UDAMP   = MAX(PSEUCO*ABS(UAVE1),ARTSSP)     ! ,APU) !CHAT  ! Alternatives   

      SELECT CASE(INLRC)

      CASE(0,1) ! Central
      UAVE    = UAVE1 !- (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP)
      CASE(2) ! One-sided (not for periodic)
         IF(KBOT == 1 .AND. INLOUT == 1) THEN
         UAVE = UL - UC
         ELSEIF(KBOT == 1 .AND. INLOUT == 2) THEN
         UAVE = UR - UC
         ELSEIF(KBOT /= 1 .AND. INLOUT == 1) THEN
         UAVE = UR - UC
         ELSEIF(KBOT /= 1 .AND. INLOUT == 2) THEN
         UAVE = UL - UC
         ENDIF
      CASE(3) ! 1st order
         IF(KBOT == 1 .AND. INLOUT == 1) THEN
         UAVE = RMIMM - UC!- (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP)
         ELSEIF(KBOT == 1 .AND. INLOUT == 2) THEN
         UAVE = RMIPP - UC!- (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP)
         ELSEIF(KBOT /= 1 .AND. INLOUT == 1) THEN
         UAVE = RMIPP - UC!- (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP)
         ELSEIF(KBOT /= 1 .AND. INLOUT == 2) THEN
         UAVE = RMIMM - UC!- (PDR+PKR - PDL-PKL)/(2.*ROHAT*UDAMP)
         ENDIF
      END SELECT

      PAVE1   = .5*(PDL+PKL + PDR+PKR)
      PAVE    = PAVE1

      IF(UAVE >= 0.) THEN
      F1PI    = UAVE*ROIM(IG) ! ROHAT
      F1PIG   = A3(II)*F1PI
      F2PIG   = F1PIG*UL + A3(II)*PAVE

      F3PIG   = F1PIG*VL
      F4PIG   = F1PIG*WL
      EIL     = EIM(IG) +.5*(UL**2 + VL**2 + WL**2) + RKL
      F5PIG   = F1PIG*EIL+
     +          .5*A3(II)*(UAVE+UC)*(PDR+PKR+PDL+PKL+2.*FRSPRE) 
      F6PIG   = F1PIG*RKL
      F7PIG   = F1PIG*REL

      ELSE

      F1PI    = UAVE*ROIP(IG) ! ROHAT
      F1PIG   = A3(II)*F1PI
      F2PIG   = F1PIG*UR + A3(II)*PAVE
      F3PIG   = F1PIG*VR
      F4PIG   = F1PIG*WR
      EIR     = EIP(IG) +.5*(UR**2 + VR**2 + WR**2) + RKR
      F5PIG   = F1PIG*EIR+
     +         .5*A3(II)*(UAVE+UC)*(PDR+PKR+PDL+PKL+2.*FRSPRE) 
      F6PIG   = F1PIG*RKR
      F7PIG   = F1PIG*RER
      ENDIF

C*    HATS SHOULD BE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.94

      HAT1(II)= 0.  
      HAT2(II)= F1PI
      HAT3(II)= HAT2(II)
      HAT4(II)= SIGN(1.,F1PIG)*HAT2(II)/ROHAT

C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG
C
C ... FLUXES
C
c      if(ngl == 1) write(6666,*)i,j,ii,ig,ijtrid,f3r(ii),f1pig,roip(ig)
c      if(ngl == 1.and.ii <= 808) write(6702,*)ii,ig,roim(ig),uave,
c     + roip(ig),pdl,pkl,rmim(ig),rmip(ig),ul,ur
      F3R(II)   = F1PIG    ; HAT1(IF1) = F3R(IF1)   ! Not tested
      F3RM(II)  = F2PM 
      F3RN(II)  = F3PM 
      F3RW(II)  = F4PM 
      F3E(II)   = F5PIG
      F3RK(II)  = F6PIG 
      F3EPS(II) = F7PIG

      ENDDO ; ENDDO

      RETURN
      END SUBROUTINE INFPER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INLET(P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,XC,YC,ZC,
     2  OHMI,S11,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,
     3  KY1,KY2,RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,CP,PRT,IFLUX,M,M2,IN,
     4  JN,KN,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,
     5  KSTR,ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,INLOUT,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,
     7  BOUNW,BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,
     8  BOUNA1,BOUNA2,MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,IPRINT,JSTATE,
     9  PRO,VAR,NPHASE,IWAVEB,FRESUL,FRSSIE,RNUT,
     1  TRANSL,TRM,BOUNG,BOUNRET,PLE2,BOUNU1,BOUNU2,BOUNV1,BOUNV2,
     2  BOUNW1,BOUNW2,MULPHC,PRC)

      USE NS3CO, ONLY : LN, IPRESC
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,INLOUT,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NN,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IPRINT,NPHASE

      INTEGER :: JSTATE(*), IWAVEB(*)

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),CP(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),BOUNA1(*),BOUNA2(*),
     8  RNUT(*),BOUNG(*),BOUNRET(*),PLE2(*),
     9  BOUNU1(*),BOUNU2(*),BOUNV1(*),BOUNV2(*),BOUNW1(*),BOUNW2(*)
c     7  TEMP1(*),TEMP2(*),ALFA1(*),ALFA2(*),X1(*),X2(*),RO1(*),RO2(*)

      REAL :: PR,PRT,PRS,FRSPRE,FRSDEN,FRSVEL,FRSVIS,T0,DIFPRE,
     2  VELOC,DPKOE,RELAX,DPDIF,ROAVE,RGAS,GAMMA,E0REF,T0REF,
     3  RK1,RK2,TURLIM,FRSSIE
      
      TYPE(PROPERTIES),       TARGET :: PRO(*)
      TYPE(MPHASE_VARIABLES), TARGET :: VAR(*)
      TYPE(INTERMITTENCY),    TARGET :: TRM(*)
      TYPE(PRE_COR),          TARGET :: PRC(*)

      REAL :: XC(*), YC(*), ZC(*)

      LOGICAL :: FRESUL, TRANSL

      CHARACTER(*) :: MULPHC
!
      RELAX   = 1.0
      RK1     = 1.
      RK2     = 0.

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

C ... Update the ghost-cell values at the inlets and outlets

      KA      = (KN+KBOT-1)*KSTR
           
      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells
         IF(FRESUL) IWAVEB(IM1) = 3

         IF(INLOUT == 4) THEN
         
         VELOC   = A3X(IF1)*U(IM1) + A3Y(IF1)*V(IM1) + A3Z(IF1)*W(IM1)

C ... Subsonic flow is assumed, extrapolation of pressure

         DPKOE       = RELAX*(RK1*P(II)     + RK2*P(II2)) +
     +                (1.- RELAX)*P(IM1) 
         DPDIF       = RELAX*(RK1*PDIFF(II) + RK2*PDIFF(II2)) +
c         DPDIF       = RELAX*(RK1*PLE2(II) + RK2*PLE2(II2)) +
     +                (1.- RELAX)*PDIFF(IM1) 
         EPS2(IM1)   = RK1*EPS2(II) + RK2*EPS2(II2)
         EPS2(IM2)   = EPS2(IM1)
         VIST(IM1)   = RK1*VIST(II) + RK2*VIST(II2)
         VIST(IM2)   = VIST(IM1)

         IF (IPRESC == 1) THEN ! Pressure gradient
            PRC(IM1)%DPDX  = RK1*PRC(II)%DPDX + RK2*PRC(II2)%DPDX
            PRC(IM2)%DPDX  = PRC(IM1)%DPDX
            PRC(IM1)%DPDY  = RK1*PRC(II)%DPDY + RK2*PRC(II2)%DPDY
            PRC(IM2)%DPDY  = PRC(IM1)%DPDY
            PRC(IM1)%DPDZ  = RK1*PRC(II)%DPDZ + RK2*PRC(II2)%DPDZ
            PRC(IM2)%DPDZ  = PRC(IM1)%DPDZ
         ENDIF ! Pressure gradient

c        DPKOE       = MAX(DPKOE,.01*FRSPRE) 

c         IF(-VELOC*IDIR  > 0.) THEN ! Inflow

c         RM(IM1)     = 2.*BOUNU(NP)*BOUNR(NP) - RM(II)
c         RN(IM1)     = 2.*BOUNV(NP)*BOUNR(NP) - RN(II)
c         RW(IM1)     = 2.*BOUNW(NP)*BOUNR(NP) - RW(II)
         RM(IM1)     = BOUNU(NP)*BOUNR(NP)
         RN(IM1)     = BOUNV(NP)*BOUNR(NP)
         RW(IM1)     = BOUNW(NP)*BOUNR(NP)
c          write(677+icycle,*) bounr(np),bounu(np),bounv(np)
         IF(ICYCLE == 1) RO(IM1)     = BOUNR(NP)



         TEMP(IM1)   = BOUNT(NP)
         TEMP(IM2)   = TEMP(IM1)
           
         IF(NPHASE > 1) THEN
         PRO(IM1)%TEMP(1)  = BOUNT(NP)
         PRO(IM1)%TEMP(2)  = BOUNT(NP)
         PRO(IM2)%TEMP(1)  = PRO(IM1)%TEMP(1)
         PRO(IM2)%TEMP(2)  = PRO(IM1)%TEMP(2)
         PRO(IM1)%DTEMP(1) = BOUNT(NP)
         PRO(IM1)%DTEMP(2) = BOUNT(NP)
         PRO(IM2)%DTEMP(1) = PRO(IM1)%DTEMP(1)
         PRO(IM2)%DTEMP(2) = PRO(IM1)%DTEMP(2)
         VAR(IM1)%ALFA(1)  = BOUNA1(NP)
         VAR(IM1)%ALFA(2)  = BOUNA2(NP)

         IF(MULPHC == 'MULTI') THEN
            VAR(IM1)%U(1) = BOUNU1(NP);    VAR(IM1)%U(2) = BOUNU2(NP)
            VAR(IM1)%V(1) = BOUNV1(NP);    VAR(IM1)%V(2) = BOUNV2(NP)
            VAR(IM1)%W(1) = BOUNW1(NP);    VAR(IM1)%W(2) = BOUNW2(NP)
            VAR(IM2)%U(1) = VAR(IM1)%U(1); VAR(IM2)%U(2) = VAR(IM1)%U(2)
            VAR(IM2)%V(1) = VAR(IM1)%V(1); VAR(IM2)%V(2) = VAR(IM1)%V(2)
            VAR(IM2)%W(1) = VAR(IM1)%W(1); VAR(IM2)%W(2) = VAR(IM1)%W(2)
         ENDIF ! MULTI

         RO(IM1) = BOUNA1(NP)*PRO(IM1)%RO(1)+BOUNA2(NP)*PRO(IM1)%RO(2)
         VAR(IM2)%ALFA(1)  = VAR(IM1)%ALFA(1)
         VAR(IM2)%ALFA(2)  = VAR(IM1)%ALFA(2)
         RO(IM2) = BOUNA1(NP)*PRO(IM2)%RO(1)+BOUNA2(NP)*PRO(IM2)%RO(2)

         PRO(IM1)%QIF(1)   = RK1*PRO(II)%QIF(1) + RK2*PRO(II2)%QIF(1)
         PRO(IM1)%QIF(2)   = RK1*PRO(II)%QIF(2) + RK2*PRO(II2)%QIF(2)
         VAR(IM1)%EVAPR(1) = RK1*VAR(II)%EVAPR(1) +RK2*VAR(II2)%EVAPR(1)
         VAR(IM1)%EVAPR(2) = RK1*VAR(II)%EVAPR(2) +RK2*VAR(II2)%EVAPR(2)

         VAR(IM1)%X(1)     = VAR(IM1)%ALFA(1)*PRO(IM1)%RO(1)/RO(IM1) !BOUNR(NP)

         IF(VAR(IM1)%X(1) > .999999) VAR(IM1)%X(1) = 1. ! Prevents error
         VAR(IM1)%X(2)     = 1.- VAR(IM1)%X(1)
         
         VAR(IM2)%X(1)     = VAR(IM2)%ALFA(1)*PRO(IM2)%RO(1)/RO(IM2) !BOUNR(NP)
         IF(VAR(IM2)%X(1) > .999999) VAR(IM2)%X(1) = 1.
         VAR(IM2)%X(2)     = 1.- VAR(IM2)%X(1)
         ENDIF
c          write(677,*) icycle,bouna1(np),bouna2(np),PRO(IM1)%RO(1),
c     2  ro(im1),var(im1)%alfa(1),VAR(IM1)%X(1)
         IF(ICYCLE >= 1) THEN
c ... More stable without these, but this is needed for a nozzle??
c         RM(IM1)     = 2.*RM(IM1) - 1.5*RM(II) +.5*RM(II2)
c         RN(IM1)     = 2.*RN(IM1) - 1.5*RN(II) +.5*RN(II2)
c         RW(IM1)     = 2.*RW(IM1) - 1.5*RW(II) +.5*RW(II2)

         ROAVE       = RO(IM1)

         RM(IM1)     = RM(IM1) + A3X(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         RN(IM1)     = RN(IM1) + A3Y(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         RW(IM1)     = RW(IM1) + A3Z(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         
         ENDIF

         RM(IM2)     = RM(IM1)
         RN(IM2)     = RN(IM1)
         RW(IM2)     = RW(IM1)

         P(IM1)      = DPKOE
         P(IM2)      = P(IM1)
         PDIFF(IM1)  = DPDIF
         PDIFF(IM2)  = PDIFF(IM1)

         BOUNP(NP)   = DPKOE
         BOUNPD(NP)  = DPDIF

C ... to get correct internal energy in calculated bc
         
         CALL EFPT(E(IM1),P(IM1),TEMP(IM1),1,ISTATE,RGAS,GAMMA,0.,
     +   FRSSIE,E0REF,T0REF)
         
         E(IM2)  = E(IM1)
         U(IM1)  = RM(IM1)/RO(IM1)
         U(IM2)  = RM(IM2)/RO(IM2)
         V(IM1)  = RN(IM1)/RO(IM1)
         V(IM2)  = RN(IM2)/RO(IM2)
         W(IM1)  = RW(IM1)/RO(IM1)
         W(IM2)  = RW(IM2)/RO(IM2)
         E(IM1)  = RO(IM1)*(E(IM1)+.5*(U(IM1)**2+V(IM1)**2+W(IM1)**2))
         E(IM2)  = RO(IM2)*(E(IM2)+.5*(U(IM2)**2+V(IM2)**2+W(IM2)**2))


         ELSE IF(INLOUT == 5) THEN
         
         VELOC   = A3X(IF1)*U(IM1) + A3Y(IF1)*V(IM1) + A3Z(IF1)*W(IM1)

C ... Supersonic inflow is assumed, all the flow variables are set

         EPS2(IM1)   = RK1*EPS2(II) + RK2*EPS2(II2)
         EPS2(IM2)   = EPS2(IM1)
         VIST(IM1)   = RK1*VIST(II) + RK2*VIST(II2)
         VIST(IM2)   = VIST(IM1)

         P(IM1)      = BOUNP(NP)
         PDIFF(IM1)  = BOUNPD(NP)
         P(IM2)      = P(IM1)
         PDIFF(IM2)  = PDIFF(IM1)

         RM(IM1)     = BOUNU(NP)*BOUNR(NP)
         RN(IM1)     = BOUNV(NP)*BOUNR(NP)
         RW(IM1)     = BOUNW(NP)*BOUNR(NP)
c          write(677+icycle,*) bounr(np),bounu(np),bounv(np)

         IF(ICYCLE == 1) RO(IM1)     = BOUNR(NP)

         TEMP(IM1)   = BOUNT(NP)
         TEMP(IM2)   = TEMP(IM1)
           
         IF(NPHASE > 1) THEN
         PRO(IM1)%TEMP(1)  = BOUNT(NP)
         PRO(IM1)%TEMP(2)  = BOUNT(NP)
         PRO(IM2)%TEMP(1)  = PRO(IM1)%TEMP(1)
         PRO(IM2)%TEMP(2)  = PRO(IM1)%TEMP(2)
         PRO(IM1)%DTEMP(1) = BOUNT(NP)
         PRO(IM1)%DTEMP(2) = BOUNT(NP)
         PRO(IM2)%DTEMP(1) = PRO(IM1)%DTEMP(1)
         PRO(IM2)%DTEMP(2) = PRO(IM1)%DTEMP(2)

         IF(MULPHC == 'MULTI') THEN
            VAR(IM1)%U(1) = BOUNU1(NP);    VAR(IM1)%U(2) = BOUNU2(NP)
            VAR(IM1)%V(1) = BOUNV1(NP);    VAR(IM1)%V(2) = BOUNV2(NP)
            VAR(IM1)%W(1) = BOUNW1(NP);    VAR(IM1)%W(2) = BOUNW2(NP)
            VAR(IM2)%U(1) = VAR(IM1)%U(1); VAR(IM2)%U(2) = VAR(IM1)%U(2)
            VAR(IM2)%V(1) = VAR(IM1)%V(1); VAR(IM2)%V(2) = VAR(IM1)%V(2)
            VAR(IM2)%W(1) = VAR(IM1)%W(1); VAR(IM2)%W(2) = VAR(IM1)%W(2)
         ENDIF ! MULTI

         VAR(IM1)%ALFA(1)  = BOUNA1(NP)
         VAR(IM1)%ALFA(2)  = BOUNA2(NP)
         RO(IM1) = BOUNA1(NP)*PRO(IM1)%RO(1)+BOUNA2(NP)*PRO(IM1)%RO(2)
         VAR(IM2)%ALFA(1)  = VAR(IM1)%ALFA(1)
         VAR(IM2)%ALFA(2)  = VAR(IM1)%ALFA(2)
         RO(IM2) = BOUNA1(NP)*PRO(IM2)%RO(1)+BOUNA2(NP)*PRO(IM2)%RO(2)

         PRO(IM1)%QIF(1)   = RK1*PRO(II)%QIF(1) + RK2*PRO(II2)%QIF(1)
         PRO(IM1)%QIF(2)   = RK1*PRO(II)%QIF(2) + RK2*PRO(II2)%QIF(2)
         VAR(IM1)%EVAPR(1) = RK1*VAR(II)%EVAPR(1) +RK2*VAR(II2)%EVAPR(1)
         VAR(IM1)%EVAPR(2) = RK1*VAR(II)%EVAPR(2) +RK2*VAR(II2)%EVAPR(2)

         VAR(IM1)%X(1)     = VAR(IM1)%ALFA(1)*PRO(IM1)%RO(1)/RO(IM1) !BOUNR(NP)

         IF(VAR(IM1)%X(1) > .999999) VAR(IM1)%X(1) = 1. ! Prevents error
         VAR(IM1)%X(2)     = 1.- VAR(IM1)%X(1)
         
         VAR(IM2)%X(1)     = VAR(IM2)%ALFA(1)*PRO(IM2)%RO(1)/RO(IM2) !BOUNR(NP)
         IF(VAR(IM2)%X(1) > .999999) VAR(IM2)%X(1) = 1.
         VAR(IM2)%X(2)     = 1.- VAR(IM2)%X(1)
         ENDIF

         IF(ICYCLE >= 1) THEN
c ... More stable without these, but this is needed for a nozzle??
c         RM(IM1)     = 2.*RM(IM1) - 1.5*RM(II) +.5*RM(II2)
c         RN(IM1)     = 2.*RN(IM1) - 1.5*RN(II) +.5*RN(II2)
c         RW(IM1)     = 2.*RW(IM1) - 1.5*RW(II) +.5*RW(II2)

         ROAVE       = RO(IM1)

         RM(IM1)     = RM(IM1) + A3X(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         RN(IM1)     = RN(IM1) + A3Y(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         RW(IM1)     = RW(IM1) + A3Z(IF1)*WROT(IF1)*(ROAVE-BOUNR(NP))
         
         ENDIF

         RM(IM2)     = RM(IM1)
         RN(IM2)     = RN(IM1)
         RW(IM2)     = RW(IM1)

C ... to get correct internal energy in calculated bc
         
         CALL EFPT(E(IM1),P(IM1),TEMP(IM1),1,ISTATE,RGAS,GAMMA,0.,
     +   FRSSIE,E0REF,T0REF)
         
         E(IM2)  = E(IM1)
         U(IM1)  = RM(IM1)/BOUNR(NP)
         U(IM2)  = RM(IM2)/BOUNR(NP)
         V(IM1)  = RN(IM1)/BOUNR(NP)
         V(IM2)  = RN(IM2)/BOUNR(NP)
         W(IM1)  = RW(IM1)/BOUNR(NP)
         W(IM2)  = RW(IM2)/BOUNR(NP)
         E(IM1)  = RO(IM1)*(E(IM1)+.5*(U(IM1)**2+V(IM1)**2+W(IM1)**2))
         E(IM2)  = RO(IM2)*(E(IM2)+.5*(U(IM2)**2+V(IM2)**2+W(IM2)**2))

      ELSE IF(INLOUT <= 6) THEN
         WRITE(*,'(A,A,I2)') ' No such inlet type yet, contact the',
     +   ' system administrator. Type =',INLOUT
         WRITE(*,*) ' Exiting in INLET...'
         STOP
      ENDIF ! INLOUT == 4

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         RK(IM1)     = BOUNRK(NP)
         RK(IM2)     = RK(IM1)
         REPS(IM1)   = BOUNEP(NP)
         REPS(IM2)   = REPS(IM1)
         E(IM1)      = E(IM1) + RK(IM1)
         E(IM2)      = E(IM2) + RK(IM2)
         IF(ITURB == 9) THEN
          RNUT(IM1) = BOUNEP(NP)/RO(IM1) ; RNUT(IM2) = RNUT(IM1)
         ENDIF
      ENDIF

      IF(TRANSL) THEN ! Intermittency variables
         TRM(IM1)%G   = BOUNG(NP)
         TRM(IM2)%G   = TRM(IM1)%G
         TRM(IM1)%RET = BOUNRET(NP)
         TRM(IM2)%RET = TRM(IM1)%RET
      ENDIF

      IF (NSCAL >= 1) THEN
      DO  560 NS = 1,NSCAL
         FI(IM1,NS)  = BOUNFI(NP,NS)
         FI(IM2,NS)  = FI(IM1,NS)
560   CONTINUE
      ENDIF

      IF (ISTRES >= 1) THEN
      DO  570 NS = 1,5
         BIJ(IM1,NS)  = BOUNBI(NP,NS)
         BIJ(IM2,NS)  = BIJ(IM1,NS) 
570   CONTINUE
      ENDIF

      IF (ITURB >= 10.AND.ITURB <= 19 .OR. ISTRES >= 1) THEN
C ... REYNOLDS STRESSES ARE TAKEN FROM COMPUTATIONAL DOMAIN IF GATSKI-TYPE ASM IS USED
         S11(IM1,1)  = S11(II,1)
         S11(IM2,1)  = S11(IM1,1)
         S11(IM1,2)  = S11(II,2)
         S11(IM2,2)  = S11(IM1,2)
         S11(IM1,3)  = S11(II,3)
         S11(IM2,3)  = S11(IM1,3)
         S11(IM1,4)  = S11(II,4)
         S11(IM2,4)  = S11(IM1,4)
         S11(IM1,5)  = S11(II,5)
         S11(IM2,5)  = S11(IM1,5)
         S11(IM1,6)  = S11(II,6)
         S11(IM2,6)  = S11(IM1,6)
      ENDIF
********************************************************************************        
      ENDDO ! Antakee mersu
      ENDDO

1000  CONTINUE

      RETURN
      END SUBROUTINE INLET
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTLET(P,PDIFF,U,V,W,C,E,VOL,A3,A3X,A3Y,A3Z,XC,YC,ZC,
     2  OHMI,S11,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,KX1,KX2,
     3  KY1,KY2,RO,RM,RN,RW,RK,REPS,WROT,TEMP,CH,CP,PRT,IFLUX,M,M2,IN,
     4  JN,KN,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,
     5  KSTR,ISTRES,ICYCLE,ISTATE,IDIR,MAX11,INTEM,INTET,N,INLOUT,
     6  RGAS,GAMMA,E0REF,T0REF,FI,BIJ,BOUNR,BOUNU,BOUNV,
     7  BOUNW,BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,
     8  BOUNA1,BOUNA2,MAXSB,MAXSS,MAXEB,IBF,NSCAL,TURLIM,IPRINT,JSTATE,
     9  PRO,VAR,NPHASE,IWAVEB,FRESUL,RNUT,
     1  TRANSL,TRM,BOUNG,BOUNRET,PLE2,BOUNU1,BOUNU2,BOUNV1,BOUNV2,
     2  BOUNW1,BOUNW2,MULPHC,if2,PRC)

      USE NS3CO, ONLY : LN, GX, GY, GZ, GROUND, IPRESC
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KSTRID,IDIM,KBOT,IFLUX,M,M2,IN,JN,KN,
     2  ITURB,ISTR,JSTR,KSTR,IDIS,MAX11,INTEM,INTET,N,INLOUT,IDIR,
     3  NFL,KA,I,J,II,IJ,II2,IM1,IM2,IF1,NI,NN,NP,NR,KX1,KX2,KY1,KY2,
     4  N1,LSTRID,NS,IBF,NSCAL,ISTRES,MAXSB,MAXEB,ISTATE,ICYCLE,MAXSS,
     5  IPRINT,IFILE,NPHASE,IPHASE,if2

      INTEGER :: JSTATE(*), IWAVEB(*)

      REAL :: P(*),V(*),U(*),C(*),RO(*),A3(*),A3Y(*),A3Z(*),A3X(*),
     2  OHMI(*),VIS(*),EPS2(*),VIST(*),CH(*),CP(*),RM(*),RN(*),RW(*),
     3  W(*),VOL(*),E(*),S11(MAXSS,*),RK(*),REPS(*),WROT(*),TEMP(*),
     4  PDIFF(*),FI(MAXSB,MAX(1,NSCAL)),BOUNR(*),BOUNU(*),BOUNV(*),
     5  BOUNW(*),BOUNE(*),BOUNRK(*),BOUNEP(*),
     6  BOUNFI(IBF,MAX(1,NSCAL)),BOUNBI(IBF,*),
     7  BIJ(MAXEB,*),BOUNT(*),BOUNP(*),BOUNPD(*),BOUNA1(*),BOUNA2(*),
     8  RNUT(*),BOUNG(*),BOUNRET(*),PLE2(*),
     9  BOUNU1(*),BOUNU2(*),BOUNV1(*),BOUNV2(*),BOUNW1(*),BOUNW2(*)
c     7  TEMP1(*),TEMP2(*),ALFA1(*),ALFA2(*),X1(*),X2(*),RO1(*),RO2(*)

      REAL :: PR,PRT,PRS,FRSPRE,FRSDEN,FRSVEL,FRSVIS,T0,DIFPRE,
     2  VELOC,DPKOE,RELAX,DPDIF,ROAVE,RGAS,GAMMA,E0REF,T0REF,
     3  RK1,RK2,A11,A12,A13,A21,A22,A23,A31,A32,A33,PALPO,UIM1,VIM1,
     4  WIM1,PLOSS3,PLOSS4,TURLIM,UKM1,VKM1,WKM1,VELC1,VELC2,VELC3,
     5  UKP1,VKP1,WKP1,UU,VV,WW,PEFEC,E1(1),RMASS,DALFA,DXX,
     6  RK3,RK4,PAAVO

      LOGICAL :: FRESUL, TRANSL

      
      TYPE(PROPERTIES),       TARGET :: PRO(*)
      TYPE(MPHASE_VARIABLES), TARGET :: VAR(*)
      TYPE(INTERMITTENCY),    TARGET :: TRM(*)
      TYPE(PRE_COR),          TARGET :: PRC(*)

      REAL :: XC(*), YC(*), ZC(*)

      CHARACTER(LEN=3) ::  NBLC
      CHARACTER(*)     ::  MULPHC
      
      RELAX   = 1.0
      RK1     = 1.
      RK2     = 0.
      IF(INLOUT == 2) THEN
         RK3  = 0.7
         RK4  = 0.3
      ELSE IF(INLOUT == 3) THEN
         RK3  = 1.
         RK4  = 0.
      ENDIF

          
      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

C ... Update the ghost-cell values at the inlets and outlets

      KA   = (KN+KBOT-1)*KSTR
          
      DO J = KY1,KY2 ! These should not be extended

      IJ   = (JN+J-1)*JSTR + KA
      NN   = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IM1     = II - LSTRID            ! First ghost cell
         IM2     = IM1- LSTRID            ! Second ghost cell
         IF1 = II  + NFL                     ! FLUX INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

         A11 = A3X(IF1)
         A12 = A3Y(IF1)
         A13 = A3Z(IF1)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

         IF(FRESUL) IWAVEB(IM1) = 5

********************************************************************************

         IF(INLOUT == 2 .OR. INLOUT == 3) THEN

C ... Pressure or total pressure is given, velocity is extrapolated
           
c         VELOC    = A11*U(KM1) + A12*V(KM1) + A13*W(KM1)
         IF(MULPHC /= 'MULTI') THEN
         UIM1     = BOUNU(NP)!/BOUNR(NP)
         VIM1     = BOUNV(NP)!/BOUNR(NP)
         WIM1     = BOUNW(NP)!/BOUNR(NP)
         ELSE ! If the phases flow to different directions, troubles!!
         UIM1     = VAR(II)%X(1)*BOUNU1(NP)+VAR(II)%X(2)*BOUNU2(NP)
         VIM1     = VAR(II)%X(1)*BOUNV1(NP)+VAR(II)%X(2)*BOUNV2(NP)
         WIM1     = VAR(II)%X(1)*BOUNW1(NP)+VAR(II)%X(2)*BOUNW2(NP)
         ENDIF

         VELOC    = .5*(A11*(U(II)+UIM1) + A12*(V(II)+VIM1) 
     +               + A13*(W(II)+WIM1))

C         RMALOC  = VELOC/C(II)

C ... Total pressure loss from the domain (no-loss can be activated)

         PLOSS3   = 0.
         PLOSS4   = 0.
c         PLOSS3   = .5*RO(II)*(U(II)**2+V(II)**2+W(II)**2)
c         PLOSS4   = .5*RO(II)*(U(II)**2+V(II)**2+W(II)**2)

C ... Hydrostatic pressure does not work here

c         PLOSS3 =  -GTOT + FRSDEN*(GX*XC(II)+GY*YC(II)+GZ*ZC(II))
c         PLOSS4 =  -GTOT + FRSDEN*(GX*XC(II)+GY*YC(II)+GZ*ZC(II))
c         PLOSS3 =
c     &               (FRSDEN-RO(IM1))*(GX*(XC(IM1) - GROUND) ! Mersu
c     &              +                  GY*(YC(IM1) - GROUND)
c     &              +                  GZ*(ZC(IM1) - GROUND))    


C ... Static or total pressure can be activated

         P(IM1)    = RK1*(BOUNP(NP) -PLOSS3) + RK2*(BOUNP(NP) -PLOSS4)
         P(IM2)    = P(IM1)

         PDIFF(IM1)= RK1*(BOUNPD(NP)-PLOSS3) + RK2*(BOUNPD(NP)-PLOSS4)
         PDIFF(IM2)= PDIFF(IM1)

         IF(-VELOC*IDIR  > 0.) THEN ! Outflow

         U(IM1)   = (RK1*RO(II)*U(II) + RK2*RO(II2)*U(II2))/RO(IM1)
         U(IM2)   = U(IM1)
         V(IM1)   = RK1*V(II) + RK2*V(II2)
         V(IM2)   = V(IM1)
         W(IM1)   = RK1*W(II) + RK2*W(II2)
         W(IM2)   = W(IM1)

         EPS2(IM2)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM2)= RK1*VIST(II) + RK2*VIST(II2)
         EPS2(IM2)= MIN(EPS2(IM2),TURLIM)
         VIST(IM2)= MIN(VIST(IM2),TURLIM*VIS(IM2))
         EPS2(IM1)= EPS2(IM2)
         VIST(IM1)= VIST(IM2)
         TEMP(IM1)= RK1*TEMP(II) + RK2*TEMP(II2)
         TEMP(IM2)= TEMP(IM1)
c         RO(IM1)  = RK1*RO(II) + RK2*RO(II2)
c         RO(IM2)  = RO(IM1)

         RM(IM1)  = RO(IM1)*U(IM1)
         RM(IM2)  = RO(IM2)*U(IM2)
         RN(IM1)  = RO(IM1)*V(IM1)
         RN(IM2)  = RO(IM2)*V(IM2)
         RW(IM1)  = RO(IM1)*W(IM1)
         RW(IM2)  = RO(IM2)*W(IM2)

         IF(NPHASE > 1) THEN

         IF(VAR(II)%ALFA(1) <= 0.5) THEN
         VAR(IM2)%ALFA(1) = RK1*VAR(II)%ALFA(1) + RK2*VAR(II2)%ALFA(1)
         VAR(IM1)%ALFA(1) = VAR(IM2)%ALFA(1)
         VAR(IM2)%X(2)    = RK1*VAR(II)%X(2) + RK2*VAR(II2)%X(2)
         VAR(IM1)%X(2)    = VAR(IM2)%X(2)
         VAR(IM2)%ALFA(2) = 1.- VAR(IM2)%ALFA(1)
         VAR(IM1)%ALFA(2) = 1.- VAR(IM1)%ALFA(1)
c         VAR(IM2)%X(2)    = PRO(IM2)%RO(2)*VAR(IM2)%ALFA(2)/RO(IM2)
c         VAR(IM1)%X(2)    = PRO(IM1)%RO(2)*VAR(IM1)%ALFA(2)/RO(IM1)
         VAR(IM2)%X(1)    = 1.- VAR(IM2)%X(2)
         VAR(IM1)%X(1)    = 1.- VAR(IM1)%X(2)
         ELSE
         VAR(IM2)%ALFA(2) = RK1*VAR(II)%ALFA(2) + RK2*VAR(II2)%ALFA(2)
         VAR(IM1)%ALFA(2) = VAR(IM2)%ALFA(2)
         VAR(IM2)%X(2)    = RK1*VAR(II)%X(2) + RK2*VAR(II2)%X(2)
         VAR(IM1)%X(2)    = VAR(IM2)%X(2)
         VAR(IM2)%ALFA(1) = 1.- VAR(IM2)%ALFA(2)
         VAR(IM1)%ALFA(1) = 1.- VAR(IM1)%ALFA(2)
c         VAR(IM2)%X(2)    = PRO(IM2)%RO(2)*VAR(IM2)%ALFA(2)/RO(IM2)
c         VAR(IM1)%X(2)    = PRO(IM1)%RO(2)*VAR(IM1)%ALFA(2)/RO(IM1)
         VAR(IM2)%X(1)    = 1.- VAR(IM2)%X(2)
         VAR(IM1)%X(1)    = 1.- VAR(IM1)%X(2)
         ENDIF

         IF(MULPHC == 'MULTI') THEN
         DO IPHASE = 1,NPHASE

c         VAR(IM1)%U(IPHASE)=(RK1*PRO(II)%RO(IPHASE)*VAR(II)%U(IPHASE)+
c     +   RK2*PRO(II2)%RO(IPHASE)*VAR(II2)%U(IPHASE))/PRO(II)%RO(IPHASE) ! IM1
         VAR(IM1)%U(IPHASE)=RK1*VAR(II)%U(IPHASE)+RK2*VAR(II2)%U(IPHASE)
         VAR(IM2)%U(IPHASE)= VAR(IM1)%U(IPHASE)
         VAR(IM1)%V(IPHASE)=RK1*VAR(II)%V(IPHASE)+RK2*VAR(II2)%V(IPHASE)
         VAR(IM2)%V(IPHASE)= VAR(IM1)%V(IPHASE)
         VAR(IM1)%W(IPHASE)=RK1*VAR(II)%W(IPHASE)+RK2*VAR(II2)%W(IPHASE)
         VAR(IM2)%W(IPHASE)= VAR(IM1)%W(IPHASE)
         ENDDO
         ENDIF

         PRO(IM2)%TEMP(1) = RK1*PRO(II)%TEMP(1) + RK2*PRO(II2)%TEMP(1)
         PRO(IM1)%TEMP(1) = PRO(IM2)%TEMP(1)
         PRO(IM2)%TEMP(2) = RK1*PRO(II)%TEMP(2) + RK2*PRO(II2)%TEMP(2)
         PRO(IM1)%TEMP(2) = PRO(IM2)%TEMP(2)
         PRO(IM2)%DTEMP(1)= RK1*PRO(II)%DTEMP(1) + RK2*PRO(II2)%DTEMP(1)
         PRO(IM1)%DTEMP(1)= PRO(IM2)%DTEMP(1)
         PRO(IM2)%DTEMP(2)= RK1*PRO(II)%DTEMP(2) + RK2*PRO(II2)%DTEMP(2)
         PRO(IM1)%DTEMP(2)= PRO(IM2)%DTEMP(2)
         PRO(IM1)%QIF(1)  = RK1*PRO(II)%QIF(1) + RK2*PRO(II2)%QIF(1)
         PRO(IM1)%QIF(2)  = RK1*PRO(II)%QIF(2) + RK2*PRO(II2)%QIF(2)
         VAR(IM1)%EVAPR(1)= RK1*VAR(II)%EVAPR(1) +RK2*VAR(II2)%EVAPR(1)
         VAR(IM1)%EVAPR(2)= RK1*VAR(II)%EVAPR(2) +RK2*VAR(II2)%EVAPR(2)

         ENDIF ! NPHASE > 1

         ELSE ! Inflow
           
         UKM1     = RK1*U(II) + RK2*U(II2) ! From the domain
         VKM1     = RK1*V(II) + RK2*V(II2) ! From the domain
         WKM1     = RK1*W(II) + RK2*W(II2) ! From the domain

         VELC1    = A11*UKM1 + A12*VKM1 + A13*WKM1

         UKP1     = BOUNU(NP)!/BOUNR(NP)
         VKP1     = BOUNV(NP)!/BOUNR(NP)
         WKP1     = BOUNW(NP)!/BOUNR(NP)

         VELC2    = A21*UKP1 + A22*VKP1 + A23*WKP1
         VELC3    = A31*UKP1 + A32*VKP1 + A33*WKP1
cc            velc2= 0.
cc            velc3=0.
         UU       = A11*VELC1 + A21*VELC2 + A31*VELC3
         VV       = A12*VELC1 + A22*VELC2 + A32*VELC3
         WW       = A13*VELC1 + A23*VELC2 + A33*VELC3

         U(IM1)   = RK3*U(IM1) + RK4*UU
         U(IM2)   = U(IM1)
         V(IM1)   = RK3*V(IM1) + RK4*VV
         V(IM2)   = V(IM1)
         W(IM1)   = RK3*W(IM1) + RK4*WW
         W(IM2)   = W(IM1)

         EPS2(IM2)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM2)= RK1*VIST(II) + RK2*VIST(II2)
         EPS2(IM2)= MIN(EPS2(IM2),TURLIM)
         VIST(IM2)= MIN(VIST(IM2),TURLIM*VIS(IM2))
         EPS2(IM1)= EPS2(IM2)
         VIST(IM1)= VIST(IM2)
         TEMP(IM1)= BOUNT(NP)
         TEMP(IM2)= TEMP(IM1)
         RO(IM1)  = BOUNR(NP)
         RO(IM2)  = RO(IM1)
         IF(NPHASE > 1) THEN
         PRO(IM1)%TEMP(1)  = RK3*PRO(IM1)%TEMP(1) +RK4*BOUNT(NP)
         PRO(IM1)%TEMP(2)  = BOUNT(NP)
         PRO(IM2)%TEMP(1)  = PRO(IM1)%TEMP(1)
         PRO(IM2)%TEMP(2)  = PRO(IM1)%TEMP(2)
c         write(667,*) bouna1(np),bouna2(np),VAR(IM1)%ALFA(1)
         DALFA             = -RK4*VAR(IM1)%ALFA(1) +RK4*BOUNA1(NP)
         DALFA = 0. ! This does not change the void fraction
         VAR(IM1)%ALFA(1)  = VAR(IM1)%ALFA(1) + DALFA
         VAR(IM1)%ALFA(2)  = 1. - VAR(IM1)%ALFA(1)
c         VAR(IM1)%ALFA(2)  = RK3*VAR(IM1)%ALFA(2) +RK4*BOUNA2(NP)
c         write(667,*) 'alfa',VAR(IM1)%ALFA(1),VAR(IM1)%ALFA(2)
         DXX               = DALFA*PRO(IM1)%RO(1)/RO(IM1)
         VAR(IM2)%ALFA(1)  = VAR(IM1)%ALFA(1)
         VAR(IM2)%ALFA(2)  = VAR(IM1)%ALFA(2)
         VAR(IM1)%X(1)     = VAR(IM1)%X(1) + DXX
         VAR(IM1)%X(2)     = 1.- VAR(IM1)%X(1)
c         VAR(IM2)%X(1)     = VAR(IM2)%ALFA(1)*PRO(IM2)%RO(1)/RO(IM2)
         VAR(IM2)%X(1)     = VAR(IM2)%X(1) + DXX
         VAR(IM2)%X(2)     = 1.- VAR(IM2)%X(1)

         PRO(IM1)%QIF(1)   = RK1*PRO(II)%QIF(1) + RK2*PRO(II2)%QIF(1)
         PRO(IM1)%QIF(2)   = RK1*PRO(II)%QIF(2) + RK2*PRO(II2)%QIF(2)
         VAR(IM1)%EVAPR(1) = RK1*VAR(II)%EVAPR(1) +RK2*VAR(II2)%EVAPR(1)
         VAR(IM1)%EVAPR(2) = RK1*VAR(II)%EVAPR(2) +RK2*VAR(II2)%EVAPR(2)
         ENDIF
         
         RM(IM1)  = RO(IM1)*U(IM1)
         RM(IM2)  = RO(IM2)*U(IM2)
         RN(IM1)  = RO(IM1)*V(IM1)
         RN(IM2)  = RO(IM2)*V(IM2)
         RW(IM1)  = RO(IM1)*W(IM1)
         RW(IM2)  = RO(IM2)*W(IM2)

         ENDIF

      ELSE IF(INLOUT == 4) THEN

C ... OUTLET BOUNDARY IS UPDATED. None of the variables is set!!
C ... All are assumed to be zero gradient. This is use with jet outflow
C ... 13.8.96 PPR

         VELOC    = A11*U(II) + A12*V(II) + A13*W(II)
         PEFEC    = P(II)

         P(IM1)   = RK1*P(II) + RK2*P(II2)
         P(IM2)   = RK1*P(IM1)+ RK2*P(II)
         PDIFF(IM1)= RK1*PDIFF(II) +  RK2*PDIFF(II2)
         PDIFF(IM2)= RK1*PDIFF(IM1)  + RK2*PDIFF(II)
         TEMP(IM1)= RK1*TEMP(II) + RK2*TEMP(II2)
         TEMP(IM2)= RK1*TEMP(IM1)+ RK2*TEMP(II)
         EPS2(IM1)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM1)= RK1*VIST(II) + RK2*VIST(II2)
         EPS2(IM1)= MIN(EPS2(IM1),TURLIM)
         VIST(IM1)= MIN(VIST(IM1),TURLIM*VIS(IM1))
         EPS2(IM2)= EPS2(IM1)
         VIST(IM2)= VIST(IM1)
         IF(NPHASE > 1) THEN
         PRO(IM2)%TEMP(1) = RK1*PRO(II)%TEMP(1) + RK2*PRO(II2)%TEMP(1)
         PRO(IM1)%TEMP(1) = PRO(IM2)%TEMP(1)
         PRO(IM2)%TEMP(2) = RK1*PRO(II)%TEMP(2) + RK2*PRO(II2)%TEMP(2)
         PRO(IM1)%TEMP(2) = PRO(IM2)%TEMP(2)
         VAR(IM2)%ALFA(1) = RK1*VAR(II)%ALFA(1) + RK2*VAR(II2)%ALFA(1)
         VAR(IM1)%ALFA(1) = VAR(IM2)%ALFA(1)
         VAR(IM2)%X(1)    = RK1*VAR(II)%X(1) + RK2*VAR(II2)%X(1)
         VAR(IM1)%X(1)    = VAR(IM2)%X(1)
         VAR(IM2)%ALFA(2) = 1.- VAR(IM2)%ALFA(1)
         VAR(IM1)%ALFA(2) = 1.- VAR(IM1)%ALFA(1)
         VAR(IM2)%X(2)    = 1.- VAR(IM2)%X(1)
         VAR(IM1)%X(2)    = 1.- VAR(IM1)%X(1)

         IF(MULPHC == 'MULTI') THEN
            VAR(IM2)%U(1) = RK1*VAR(II)%U(1) + RK2*VAR(II2)%U(1)
            VAR(IM1)%U(1) = VAR(IM2)%U(1)
            VAR(IM2)%V(1) = RK1*VAR(II)%V(1) + RK2*VAR(II2)%V(1)
            VAR(IM1)%V(1) = VAR(IM2)%V(1)
            VAR(IM2)%W(1) = RK1*VAR(II)%W(1) + RK2*VAR(II2)%W(1)
            VAR(IM1)%W(1) = VAR(IM2)%W(1)
         ENDIF

         ENDIF

         RO(IM1)  = RK1*RO(II) + RK2*RO(II2)
         RO(IM2)  = RK1*RO(IM1)+ RK2*RO(II)
         RM(IM1)  = RK1*RM(II) + RK2*RM(II2)
         RM(IM2)  = RK1*RM(IM1)+ RK2* RM(II)
         RN(IM1)  = RK1*RN(II) + RK2*RN(II2)
         RN(IM2)  = RK1*RN(IM1)+ RK2* RN(II)
         RW(IM1)  = RK1*RW(II) + RK2*RW(II2)
         RW(IM2)  = RK1*RW(IM1)+ RK2* RW(II)

         U(IM1)   = RM(IM1)/RO(IM1)
         U(IM2)   = RM(IM2)/RO(IM2)
         V(IM1)   = RN(IM1)/RO(IM1)
         V(IM2)   = RN(IM2)/RO(IM2)
         W(IM1)   = RW(IM1)/RO(IM1)
         W(IM2)   = RW(IM2)/RO(IM2)

      IF(ISTATE >= 8) THEN
         WRITE(*,*) 'Outlet type 4 requires EFPRO'
         STOP
      ENDIF

C ... Internal and total energy

      CALL EFPRO(E1(1),P(IM1),RO(IM1),1,ISTATE,GAMMA,0.,E0REF,T0REF)
      E(IM1) = RO(IM1)*(E1(1) + .5*(U(IM1)**2 + V(IM1)**2
     +         +W(IM1)**2))
      CALL EFPRO(E1(1),P(IM2),RO(IM2),1,ISTATE,GAMMA,0.,E0REF,T0REF)
      E(IM2) = RO(IM2)*(E1(1) + .5*(U(IM2)**2 + V(IM2)**2
     +         +W(IM2)**2))

      ELSE IF(INLOUT == 5) THEN

C ... u-velocity is given. Only for special purposes as engine inlet or 
C ... boundary-layer upper surface

         VELOC    = A11*U(IM2) + A12*V(IM2) + A13*W(IM2)
c         RMALOC  = VELOC/C(KP2)

         PEFEC    = BOUNP(NP)
         P(IM2)   = 2.*PEFEC - P(II)
         P(IM1)   = P(IM2)
         PDIFF(IM2)   = 2.*BOUNPD(NP) - PDIFF(II)
         PDIFF(IM1)   = PDIFF(IM2)

         EPS2(IM2)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM2)= RK1*VIST(II) + RK2*VIST(II2)
         EPS2(IM2)= MIN(EPS2(IM2),TURLIM)
         VIST(IM2)= MIN(VIST(IM2),TURLIM*VIS(II))
         EPS2(IM1)= EPS2(IM2)
         VIST(IM1)= VIST(IM2)

         P(IM1)   = RK1*P(II) + RK2*P(II2)
         P(IM2)   = P(IM1)
         PDIFF(IM1)  = RK1*PDIFF(II) + RK2*PDIFF(II2)
         PDIFF(IM2)  = PDIFF(IM1)

         RO(IM1)  = RK1*RO(II) + RK2*RO(II2)
         RO(IM2)  = RO(IM1)
         RM(IM1)  = 2.*BOUNR(NP)*BOUNU(NP) - RM(II)
         RM(IM2)  = 2.*RM(IM1)   - RM(II)
         U(IM1)   = RM(IM1)/RO(IM1)
         U(IM2)   = RM(IM2)/RO(IM2)
         RN(IM1)  = RK1*RN(II)   + RK2*RN(II2)
         RN(IM2)  = RN(IM1)
         V(IM1)   = RN(IM1)/RO(IM1)
         V(IM2)   = RN(IM2)/RO(IM2)
         RW(IM1)  = RK1*RW(II)   + RK2*RW(II2)
         RW(IM2)  = RW(IM1)
         W(IM1)   = RW(IM1)/RO(IM1)
         W(IM2)   = RW(IM2)/RO(IM2)
      IF(ISTATE >= 8) THEN
         WRITE(*,*) 'Outlet type 5 requires EFPRO'
         STOP
      ENDIF
         CALL EFPRO(E(IM1),P(IM1),RO(IM1),1,ISTATE,GAMMA,0.,E0REF,T0REF)
         CALL CONSVE(E(IM1),U(IM1),V(IM1),W(IM1),RK(IM1),RO(IM1),XC(IM1)
     +   ,YC(IM1),ZC(IM1),1,ITURB,FI,MAXSB,NSCAL)
         E(IM2)   = E(IM1)

      ELSE IF(INLOUT == 6) THEN

C ... u-velocity is given. Only for special purposes as engine inlet or 
C ... boundary-layer upper surface. Differs from no. 5 only concerning
C ... temperature extrapolation

         VELOC    = A11*U(IM2) + A12*V(IM2) + A13*W(IM2)
c         RMALOC  = VELOC/C(KP2)
         PEFEC    = BOUNP(NP)

         EPS2(IM2)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM2)= RK1*VIST(II) + RK2*VIST(II2)
         TEMP(IM2)= RK1*TEMP(II) + RK2*TEMP(II2)
         EPS2(IM2)= MIN(EPS2(IM2),TURLIM)
         VIST(IM2)= MIN(VIST(IM2),TURLIM*VIS(II))
         EPS2(IM1)= EPS2(IM2)
         VIST(IM1)= VIST(IM2)
         TEMP(IM1)= TEMP(IM2)

         P(IM1)   = RK1*P(II) + RK2*P(II2)
         P(IM2)   = P(IM1)
         PDIFF(IM1) = RK1*PDIFF(II) + RK2*PDIFF(II2)
         PDIFF(IM2) = PDIFF(IM1)

         RO(IM1)  = RK1*RO(II) + RK2*RO(II2)
         RO(IM2)  = RO(IM1)
         RM(IM1)  = 2.*BOUNR(NP)*BOUNU(NP) - RM(II)
         RM(IM2)  = 2.*RM(IM1)   - RM(II)
         U(IM1)   = RM(IM1)/RO(IM1)
         U(IM2)   = RM(IM2)/RO(IM2)
         RN(IM1)  = RK1*RN(II)   + RK2*RN(II2)
         RN(IM2)  = RN(IM1)
         V(IM1)   = RN(IM1)/RO(IM1)
         V(IM2)   = RN(IM2)/RO(IM2)
         RW(IM1)  = RK1*RW(II)   + RK2*RW(II2)
         RW(IM2)  = RW(IM1)
         W(IM1)   = RW(IM1)/RO(IM1)
         W(IM2)   = RW(IM2)/RO(IM2)
         CALL EFPRO(E(IM1),P(IM1),RO(IM1),1,ISTATE,GAMMA,0.,E0REF,T0REF)
         CALL CONSVE(E(IM1),U(IM1),V(IM1),W(IM1),RK(IM1),RO(IM1),XC(IM1)
     +   ,YC(IM1),ZC(IM1),1,ITURB,FI,MAXSB,NSCAL)
         E(IM2)   = E(IM1)

      ELSE IF(INLOUT == 7) THEN

C ... normal velocity is given. Density is extrapolated

         VELOC    = A11*U(IM2) + A12*V(IM2) + A13*W(IM2)
c         RMALOC  = VELOC/C(KP2)

         EPS2(IM2)= RK1*EPS2(II) + RK2*EPS2(II2)
         VIST(IM2)= RK1*VIST(II) + RK2*VIST(II2)
         EPS2(IM2)= MIN(EPS2(IM2),TURLIM)
         VIST(IM2)= MIN(VIST(IM2),TURLIM*VIS(II))
         EPS2(IM1)= EPS2(IM2)
         VIST(IM1)= VIST(IM2)

         VELC1    = A11*BOUNU(NP) + A12*BOUNV(NP) + A13*BOUNW(NP)
         UKM1     = RK1*U(II) + RK2*U(II2) ! From the domain
         VKM1     = RK1*V(II) + RK2*V(II2) ! From the domain
         WKM1     = RK1*W(II) + RK2*W(II2) ! From the domain
         VELC2    = A21*UKM1 + A22*VKM1 + A23*WKM1
         VELC3    = A31*UKM1 + A32*VKM1 + A33*WKM1

         UU       = A11*VELC1 + A21*VELC2 + A31*VELC3
         VV       = A12*VELC1 + A22*VELC2 + A32*VELC3
         WW       = A13*VELC1 + A23*VELC2 + A33*VELC3
         U(IM1)   = UU
         U(IM2)   = UU
         V(IM1)   = VV
         V(IM2)   = VV
         W(IM1)   = WW
         W(IM2)   = WW

         RMASS    = RK1*RO(II) *(A11*U(II)  + A12*V(II)  + A13*W(II))
     +            + RK2*RO(II2)*(A11*U(II2) + A12*V(II2) + A13*W(II2))
         IF(-RMASS*IDIR > 0.) RO(IM1) = RMASS/VELOC ! Outflow
         RO(IM2)  = RO(IM1)
          
         RM(IM1)  = U(IM1)*RO(IM1)
         RM(IM2)  = U(IM2)*RO(IM2)
         RN(IM1)  = V(IM1)*RO(IM1)
         RN(IM2)  = V(IM2)*RO(IM2)
         RW(IM1)  = W(IM1)*RO(IM1)
         RW(IM2)  = W(IM2)*RO(IM2)
         
         TEMP(IM1)= RK1*TEMP(II) + RK2*TEMP(II2)
         TEMP(IM2)= TEMP(IM1)

         IF(ISTATE >= 6) THEN
         PDIFF(IM1) = RK1*PDIFF(II)+RK2*PDIFF(II2)
         PDIFF(IM2) = PDIFF(IM1)
         P(IM1)     = RK1*P(II)+RK2*P(II2)
         P(IM2)     = P(IM1)
         ELSE IF(ISTATE <= 5) THEN
         P(IM1)     = RGAS*TEMP(IM1)*RO(IM1)
         P(IM2)     = P(IM1)
         PDIFF(IM1) = P(IM1) - FRSPRE
         PDIFF(IM2) = PDIFF(IM1)
         ENDIF

         CALL EFPRO(E(IM1),P(IM1),RO(IM1),1,ISTATE,GAMMA,0.,E0REF,T0REF)
         CALL CONSVE(E(IM1),U(IM1),V(IM1),W(IM1),RK(IM1),RO(IM1),XC(IM1)
     +   ,YC(IM1),ZC(IM1),1,ITURB,FI,MAXSB,NSCAL)
         E(IM2)   = E(IM1)
         
      ELSE

         WRITE(*,'(A,A,I2)') ' No such outlet type yet, contact the',
     +   ' system administrator. Type =',INLOUT
         WRITE(*,*) 'Exiting in OUTLET...'
         STOP
      ENDIF ! INLOUT == 2

********************************************************************************

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         IF(-VELOC*IDIR  > 0.) THEN ! Outflow
         RK(IM1)     = RK1*RK(II) + RK2*RK(II2)
         RK(IM2)     = RK(IM1)
         REPS(IM1)   = RK1*REPS(II) + RK2*REPS(II2)
         REPS(IM2)   = REPS(IM1)
         E(IM1)      = E(IM1) + RK(IM1)
         E(IM2)      = E(IM2) + RK(IM2)
         IF(ITURB == 9) THEN
           RNUT(IM1) = REPS(IM1)/RO(IM1) ; RNUT(IM2) = RNUT(IM1)
         ENDIF

         ELSE
C ...Stay as it is
         ENDIF
      ENDIF

      IF (TRANSL) THEN ! Intermittency variables
         IF(-VELOC*IDIR  > 0.) THEN ! Outflow
         TRM(IM1)%G  = RK1*TRM(II)%G + RK2*TRM(II2)%G
         TRM(IM2)%G  = TRM(IM1)%G
         TRM(IM1)%RET= RK1*TRM(II)%RET + RK2*TRM(II2)%RET
         TRM(IM2)%RET= TRM(IM1)%RET
         ELSE
C ...Stay as it is
         ENDIF
      ENDIF ! Intermittency

      IF (IPRESC == 1) THEN ! Pressure gradient
         IF(-VELOC*IDIR  > 0.) THEN ! Outflow
         PRC(IM1)%DPDX  = RK1*PRC(II)%DPDX + RK2*PRC(II2)%DPDX
         PRC(IM2)%DPDX  = PRC(IM1)%DPDX
         PRC(IM1)%DPDY  = RK1*PRC(II)%DPDY + RK2*PRC(II2)%DPDY
         PRC(IM2)%DPDY  = PRC(IM1)%DPDY
         PRC(IM1)%DPDZ  = RK1*PRC(II)%DPDZ + RK2*PRC(II2)%DPDZ
         PRC(IM2)%DPDZ  = PRC(IM1)%DPDZ
         ELSE
C ...Stay as it is
         ENDIF
      ENDIF ! Pressure gradient

      IF (NSCAL >= 1) THEN
      DO  560 NS = 1,NSCAL
         FI(IM1,NS)  = RK1*FI(II,NS) + RK2*FI(II2,NS)
         FI(IM2,NS)  = FI(IM1,NS)
560   CONTINUE
      ENDIF

      IF (ISTRES >= 1) THEN
      DO  570 NS = 1,5
         BIJ(IM1,NS)  = RK1*BIJ(II,NS) + RK2*BIJ(II2,NS)
         BIJ(IM2,NS)  = BIJ(IM1,NS) 
570   CONTINUE
      ENDIF

      IF (ITURB >= 10.AND.ITURB <= 19 .OR. ISTRES >= 1) THEN
C ... REYNOLDS STRESSES ARE TAKEN FROM COMPUTATIONAL DOMAIN IF GATSKI-TYPE ASM IS USED
         S11(IM1,1)  = S11(II,1)
         S11(IM2,1)  = S11(IM1,1)
         S11(IM1,2)  = S11(II,2)
         S11(IM2,2)  = S11(IM1,2)
         S11(IM1,3)  = S11(II,3)
         S11(IM2,3)  = S11(IM1,3)
         S11(IM1,4)  = S11(II,4)
         S11(IM2,4)  = S11(IM1,4)
         S11(IM1,5)  = S11(II,5)
         S11(IM2,5)  = S11(IM1,5)
         S11(IM1,6)  = S11(II,6)
         S11(IM2,6)  = S11(IM1,6)
      ENDIF
********************************************************************************

      ENDDO
      ENDDO

1000  CONTINUE

C ... WRITE OUTPUT DIST FILE MOVED TO BOUFLO


      RETURN
      END SUBROUTINE OUTLET
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUFLO(F3R,F3RM,F3RN,F3RW,F3E,P,U,V,W,C,E,VOL,A3,A3X,
     2  A3Y,A3Z,OHMI,SS11,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     3  PDI,KX1,KX2,KY1,KY2,PR,PRT,RO,RM,RN,RW,
     4  CX,CY,CZ,DX,DY,DZ,CMX,CMY,CMZ,TOM,IFLUX,M,M2,IN,JN,KN,XC,YC,ZC,
     5  RKLIM,IDIR,F3RK,F3EPS,RK,REPS,UROT,TEMP,CH,FRSDEN,FRSPRE,FRSVEL,
     6  FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,IDIS,
     7  XCP,YCP,ZCP,UWALL,VWALL,WWALL,CPWALL,TWALL,QWALL,QWFRIC,
     8  HFLUX,TAUW1,TAUW2,TAUWX,TAUWY,TAUWZ,SURLE,DSURLE,WMFLUX,POROS,
     9  WHSTAG,WTEMP,RBK,RSDIRX,RSDIRY,RSDIRZ,SURFX,SURFY,SURFZ,
     +  MAX11,XCO,YCO,ZCO,QMFIN,QMEIN,QMFOUT,QMEOUT,QMXIN,QMXOUT,
     1  QMYIN,QMYOUT,QMZIN,QMZOUT,QVFIN,QVFOUT,INTEM,INTET,N,
     2  INLOUT,XMOM,YMOM,ZMOM,CMXIN,CMYIN,CMZIN,
     3  CMXOUT,CMYOUT,CMZOUT,INLETS,nsp,CONVL,PRO,VAR,MULPHL,
     4  QGFIN,QGEIN,QGFOUT,QGEOUT,ICYCLE,IPRINT,NSCAL,ISTRES,NPHASE,
     5  TRANSL,TRM,BOUNG,BOUNRET,BIJ,MAXSB,MAXEB,ICMAX,SURFPX,SURFPY,
     6  SURFPZ,QLFIN,QLEIN,QLFOUT,QLEOUT)

      USE NS3CO, ONLY : LN, GX, GY, GZ, GROUND
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IDIR,IDIS,ITURB,MAX11,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     2  KX1,KX2,KY1,KY2,IFLUX,M,IN,JN,KN,ISTR,JSTR,KSTR,NFL,IHFL,IHT,
     3  KA,I,J,NN,NR,NI,N1,II,IJ,II2,IF1,IF2,NH,NP,IZZ,IF3,II3,IA,N,
     4  INTEM,INTET,LSTRID,INLOUT,M2,IM1,IM2,nsp,ICYCLE,IPRINT,NS,
     5  NSCAL,IFILE,ISTRES,NPHASE,MAXSB,MAXEB,ICMAX

      REAL :: F3R(*),F3RW(*),F3E(*),P(*),V(*),U(*),C(*),RO(*),A3(*),
     2  A3Y(*),A3Z(*),A3X(*),OHMI(*),VIS(*),EPS2(*),VIST(*),PDI(*),
     3  F3RN(*),W(*),VOL(*),F3RM(*),E(*),SS11(MAX11,*),
     4  RK(*),REPS(*),F3RK(*),F3EPS(*),UROT(*),TEMP(*),CH(*),
     5  RM(*),RN(*),RW(*),UWALL(*),VWALL(*),WWALL(*),CPWALL(*),TWALL(*),
     +  QWALL(*),QWFRIC(*),HFLUX(*),TAUW1(*),TAUW2(*),
     1  TAUWX(*),TAUWY(*),TAUWZ(*),SURLE(*),DSURLE(*),
     2  WMFLUX(*),POROS(*),WHSTAG(*),WTEMP(*),RBK(*),
     3  RSDIRX(*),RSDIRY(*),RSDIRZ(*),FI(MAXSB,MAX(1,NSCAL)),
     4  SURFX(*),SURFY(*),SURFZ(*),BOUNG(*),BOUNRET(*),BIJ(MAXEB,*),
     5  SURFPX(*),SURFPY(*),SURFPZ(*)

      REAL :: PR,PRT,TOM,FRSDEN,FRSVEL,T0,DIFPRE,
     2  RC1,RC2,RM1,RM2,ROMFAC,REPFAC,PRS,YMR,
     3  YPG,RSUM,EPS,TSURF,VISLIM,CHLIM,ROLIM,TEMPB,SURF,SURG,EPS2B,
     4  VISB,SURFE,VISTB,CHSB,YPRO1,YPRO2,RKLIM,FRSPRE,FRSVIS,YP1,YP2,
     5  YP4,ROM,POM,EKIN,RMM,RNN,RWW,QMFIN,QMFOUT,QMEIN,QMEOUT,
     6  QMXIN,QMXOUT,QMYIN,QMYOUT,QMZIN,QMZOUT,F2V,F3V,F4V,F5V,F6V,PB,
     7  PB1,PB4,X,Y,Z,XD,YD,ZD,BX,BY,BZ,DELTAP,
     8  CMXIN,CMYIN,CMZIN,CMXOUT,CMYOUT,CMZOUT,CX,CY,CZ,
     9  CMX,CMY,CMZ,DX,DY,DZ,PDM,QVFIN,QVFOUT,
     1  QGFIN,QGEIN,QGFOUT,QGEOUT,QLFIN,QLEIN,QLFOUT,QLEOUT

      REAL :: XCO(*), YCO(*), ZCO(*), XC(*), YC(*), ZC(*),
     2                    XCP(*),YCP(*),ZCP(*)
      REAL :: XMOM,YMOM,ZMOM

      LOGICAL :: INLETS, CONVL, MULPHL, TRANSL

      CHARACTER(LEN=3) ::  NBLC


      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY),    TARGET :: TRM(*)

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      RC1     = 1.5
      RC2     = RC1 - 1.

C ... Calculate balances on the inlets and outlets

      DO 1000 J = KY1,KY2 ! These should not be extended

      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1000 I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IF2 = IF1 + LSTRID                  ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

         ROM     = .5*(RO(II) + RO(II-LSTRID))
         POM     = .5*(P(II)  + P(II-LSTRID))
         PDM     = .5*(PDI(II)+ PDI(II-LSTRID))
         RMM     = .5*(RM(II) + RM(II-LSTRID))
         RNN     = .5*(RN(II) + RN(II-LSTRID))
         RWW     = .5*(RW(II) + RW(II-LSTRID))
          
         EKIN    = .5*(RMM**2 + RNN**2 + RWW**2)/ROM

      IF(INLOUT <= 0) THEN ! A direct boundary, where patch data exists
         X       = XCP(NP) - GROUND
         Y       = YCP(NP) - GROUND
         Z       = ZCP(NP) - GROUND
      ELSE ! Carbage
         X       = 1.E6
         Y       = 1.E6
         Z       = 1.E6
      ENDIF

C ... VECTOR FOR THE MOMENTUM. CENTERPOINT COORDINATES ARE UTILIZED

         XD      = X - XMOM
         YD      = Y - YMOM
         ZD      = Z - ZMOM

         DELTAP  = FRSDEN*(GX*X + GY*Y+ GZ*Z)

         BX      = -IDIR*(F3RM(IF1)-A3(IF1)*A3X(IF1)*( -DELTAP)) 
         BY      = -IDIR*(F3RN(IF1)-A3(IF1)*A3Y(IF1)*( -DELTAP))
         BZ      = -IDIR*(F3RW(IF1)-A3(IF1)*A3Z(IF1)*( -DELTAP))

      IF(M == 1) THEN ! Store the non-cumulative coefficients

         CMX     = CMX + BZ*YD - BY*ZD
         CMY     = CMY + BX*ZD - BZ*XD
         CMZ     = CMZ + BY*XD - BX*YD

         CX      = CX  + BX
         CY      = CY  + BY
         CZ      = CZ  + BZ

         DX      = DX  - IDIR*(A3(IF1)*A3X(IF1)*(POM + DELTAP))
         DY      = DY  - IDIR*(A3(IF1)*A3Y(IF1)*(POM + DELTAP))
         DZ      = DZ  - IDIR*(A3(IF1)*A3Z(IF1)*(POM + DELTAP))

      END IF ! M == 1

      IF(M2 == 1) THEN ! Integrate total fluxes in and out

         IF(INLETS .AND. CONVL) THEN ! Inlets
          QMFIN   = QMFIN - IDIR*F3R(IF1)  ! FLUX-BASED MASS FLOW
          QVFIN   = QVFIN - IDIR*F3R(IF1)/ROM
          QMXIN   = QMXIN + BX !-IDIR*F3RM(IF1)
          QMYIN   = QMYIN + BY !-IDIR*F3RN(IF1)
          QMZIN   = QMZIN + BZ !-IDIR*F3RW(IF1)
          CMXIN   = CMXIN + BZ*YD - BY*ZD
          CMYIN   = CMYIN + BX*ZD - BZ*XD
          CMZIN   = CMZIN + BY*XD - BX*YD
c           QMEIN      = QMEIN + F3E(IF1) ! Unused total energy flux
          QMEIN   = QMEIN - IDIR*F3R(IF1)*(EKIN + POM)/ROM
          IF(MULPHL) THEN
             QGFIN   = QGFIN - IDIR*VAR(IF1)%FRO(2)
             QGEIN   = QGEIN - IDIR*VAR(IF1)%FE(2)
             QLFIN   = QLFIN - IDIR*VAR(IF1)%FRO(1)
             QLEIN   = QLEIN - IDIR*VAR(IF1)%FE(1)
          ENDIF
         ELSE IF(.NOT.INLETS .AND. CONVL) THEN              ! Outlets
          QMFOUT  = QMFOUT - IDIR*F3R(IF1)  ! FLUX-BASED MASS FLOW
          QVFOUT  = QVFOUT - IDIR*F3R(IF1)/ROM  
          QMXOUT  = QMXOUT + BX !-IDIR*F3RM(IF1)
          QMYOUT  = QMYOUT + BY !-IDIR*F3RN(IF1)
          QMZOUT  = QMZOUT + BZ !-IDIR*F3RW(IF1)
          CMXOUT  = CMXOUT + BZ*YD - BY*ZD
          CMYOUT  = CMYOUT + BX*ZD - BZ*XD
          CMZOUT  = CMZOUT + BY*XD - BX*YD
          QMEOUT  = QMEOUT - IDIR*F3R(IF1)*(EKIN + POM)/ROM
          IF(MULPHL) THEN
             QGFOUT  = QGFOUT - IDIR*VAR(IF1)%FRO(2)
             QGEOUT  = QGEOUT - IDIR*VAR(IF1)%FE(2)
             QLFOUT  = QLFOUT - IDIR*VAR(IF1)%FRO(1)
             QLEOUT  = QLEOUT - IDIR*VAR(IF1)%FE(1)
c          write(666,*) qgfout,qgeout,qlfout,qleout
          ENDIF
         ENDIF

         IF(INLOUT == 0) THEN ! Store i/o data
          WMFLUX(NP)= F3R(IF1) /(A3(IF1) + 1.E-20)
          SURFX(NP) = BX/(A3(IF1) + 1.E-20)
          SURFY(NP) = BY/(A3(IF1) + 1.E-20)
          SURFZ(NP) = BZ/(A3(IF1) + 1.E-20)
c          PB4       = MAX(-.99*FRSPRE, RC1*PDI(II) - RC2*PDI(II2))
          PB4       = MAX(-.99*FRSPRE, PDM)
          CPWALL(NP)= PB4
          TWALL(NP) = RC1*TEMP(II) - RC2*TEMP(II2)
         ENDIF ! Store i/o data

      ENDIF ! M2 == 1

1000  CONTINUE


      IF(ICYCLE == IPRINT .AND. M2 == 1 .OR.
     &   ICYCLE >= ICMAX  .AND. M2 == 1) THEN

          
      CALL NUMCH3(NBLC,N)
      OPEN(51,FILE='DIST_CONS'//NBLC//'.OUT',FORM='FORMATTED')
      OPEN(52,FILE='DIST_PRIM'//NBLC//'.OUT',FORM='FORMATTED')

      IFILE = 0 ! This a recommended value
      WRITE(51,*) NSCAL
      WRITE(51,*) KX2-KX1+1,KY2-KY1+1
      WRITE(52,*) NSCAL,IFILE
      WRITE(52,*) KX2-KX1+1,KY2-KY1+1

      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) RO(IM1),RM(IM1),RN(IM1),RW(IM1),E(IM1)
         WRITE(52,*) P(IM1),U(IM1),V(IM1),W(IM1),TEMP(IM1)
      ENDDO
      ENDDO

C ... WRITE TURBULENCE VARIABLES (K-EPSILON)
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) RK(IM1),REPS(IM1)
         WRITE(52,*) RK(IM1),REPS(IM1)
      ENDDO
      ENDDO
      ENDIF

C ... WRITE SCALAR VALUES
      IF(NSCAL >= 1) THEN
      DO 900 NS = 1,NSCAL
      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) FI(IM1,NS)
         WRITE(52,*) FI(IM1,NS)
      ENDDO
      ENDDO
 900  CONTINUE
      ENDIF

C ... WRITE EARSM VALUES
      IF(ISTRES >= 1) THEN
      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) BIJ(IM1,1),BIJ(IM1,2),BIJ(IM1,3),BIJ(IM1,4),
     +   BIJ(IM1,5)
         WRITE(52,*) BIJ(IM1,1),BIJ(IM1,2),BIJ(IM1,3),BIJ(IM1,4),
     +   BIJ(IM1,5)
      ENDDO
      ENDDO
      ENDIF

C ... WRITE MULTIPHASE VALUES
      IF(NPHASE > 1) THEN
      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) VAR(IM1)%ALFA(1),VAR(IM1)%ALFA(2)
         WRITE(52,*) VAR(IM1)%ALFA(1),VAR(IM1)%ALFA(2)
      ENDDO
      ENDDO
      ENDIF

      IF (TRANSL) THEN ! Intermittency variables
      DO J = KY1,KY2 ! These should not be extended
      IJ   = (JN+J-1)*JSTR + KA
      DO I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         IM1     = II - LSTRID            ! First ghost cell
         WRITE(51,*) TRM(IM1)%G,TRM(IM1)%RET
         WRITE(52,*) TRM(IM1)%G,TRM(IM1)%RET
      ENDDO
      ENDDO
      ENDIF

      CLOSE(51) ; CLOSE(52)
      ENDIF ! ICYCLE == IPRINT

      RETURN
      END SUBROUTINE BOUFLO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUMIR(F3R,F3RM,F3RN,F3RW,F3E,P,U,V,W,C,E,VOL,A3,A3X,
     2  A3Y,A3Z,OHMI,SS11,VIS,EPS2,VIST,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     3  PDI,KX1,KX2,KY1,KY2,PR,PRT,RO,RM,RN,RW,
     4  CX,CY,CZ,DX,DY,DZ,CMX,CMY,CMZ,TOM,IFLUX,M,M2,IN,JN,KN,XC,YC,ZC,
     5  RKLIM,IDIR,F3RK,F3EPS,RK,REPS,UROT,TEMP,CH,FRSDEN,FRSPRE,FRSVEL,
     6  FRSVIS,T0,DIFPRE,ITURB,ISTR,JSTR,KSTR,IDIS,
     7  XCP,YCP,ZCP,UWALL,VWALL,WWALL,CPWALL,TWALL,QWALL,QWFRIC,
     8  HFLUX,TAUW1,TAUW2,TAUWX,TAUWY,TAUWZ,SURLE,DSURLE,WMFLUX,POROS,
     9  WHSTAG,WTEMP,RBK,RSDIRX,RSDIRY,RSDIRZ,SURFX,SURFY,SURFZ,
     +  MAX11,XCO,YCO,ZCO,QMFIN,QMEIN,QMFOUT,QMEOUT,QMXIN,QMXOUT,
     1  QMYIN,QMYOUT,QMZIN,QMZOUT,QVFIN,QVFOUT,INTEM,INTET,N,
     2  INLOUT,XMOM,YMOM,ZMOM,CMXIN,CMYIN,CMZIN,
     3  CMXOUT,CMYOUT,CMZOUT,INLETS,nsp,CONVL,PRO,VAR,MULPHL,
     4  QGFIN,QGEIN,QGFOUT,QGEOUT,IDIMP,JDIMP,SURFPX,SURFPY,SURFPZ)

      USE NS3CO, ONLY : LN, GX, GY, GZ, GROUND
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IDIR,IDIS,ITURB,MAX11,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     2  KX1,KX2,KY1,KY2,IFLUX,M,IN,JN,KN,ISTR,JSTR,KSTR,NFL,IHFL,IHT,
     3  KA,I,J,NN,NR,NI,N1,II,IJ,II2,IF1,IF2,NH,NP,IZZ,IF3,II3,IA,N,
     4  INTEM,INTET,LSTRID,INLOUT,M2,IM1,IM2,nsp,IDIMP,JDIMP

      REAL :: F3R(*),F3RW(*),F3E(*),P(*),V(*),U(*),C(*),RO(*),A3(*),
     2  A3Y(*),A3Z(*),A3X(*),OHMI(*),VIS(*),EPS2(*),VIST(*),PDI(*),
     3  F3RN(*),W(*),VOL(*),F3RM(*),E(*),SS11(MAX11,*),
     4  RK(*),REPS(*),F3RK(*),F3EPS(*),UROT(*),TEMP(*),CH(*),
     5  RM(*),RN(*),RW(*),UWALL(*),VWALL(*),WWALL(*),CPWALL(*),TWALL(*),
     +  QWALL(*),QWFRIC(*),HFLUX(*),TAUW1(*),TAUW2(*),
     1  TAUWX(*),TAUWY(*),TAUWZ(*),SURLE(*),DSURLE(*),
     2  WMFLUX(*),POROS(*),WHSTAG(*),WTEMP(*),RBK(*),
     3  RSDIRX(*),RSDIRY(*),RSDIRZ(*),
     4  SURFX(*),SURFY(*),SURFZ(*),SURFPX(*),SURFPY(*),SURFPZ(*)

      REAL :: PR,PRT,TOM,FRSDEN,FRSVEL,T0,DIFPRE,
     2  RC1,RC2,RM1,RM2,ROMFAC,REPFAC,PRS,YMR,
     3  YPG,RSUM,EPS,TSURF,VISLIM,CHLIM,ROLIM,TEMPB,SURF,SURG,EPS2B,
     4  VISB,SURFE,VISTB,CHSB,YPRO1,YPRO2,RKLIM,FRSPRE,FRSVIS,YP1,YP2,
     5  YP4,ROM,POM,EKIN,RMM,RNN,RWW,QMFIN,QMFOUT,QMEIN,QMEOUT,
     6  QMXIN,QMXOUT,QMYIN,QMYOUT,QMZIN,QMZOUT,F2V,F3V,F4V,F5V,F6V,PB,
     7  PB1,PB4,X,Y,Z,XD,YD,ZD,BX,BY,BZ,DELTAP,
     8  CMXIN,CMYIN,CMZIN,CMXOUT,CMYOUT,CMZOUT,CX,CY,CZ,
     9  CMX,CMY,CMZ,DX,DY,DZ,PDM,QVFIN,QVFOUT,
     1  QGFIN,QGEIN,QGFOUT,QGEOUT

      REAL :: XCO(*), YCO(*), ZCO(*), XC(*), YC(*), ZC(*),
     2                    XCP(*), YCP(*), ZCP(*)

      REAL :: XMOM,YMOM,ZMOM

      LOGICAL :: INLETS, CONVL, MULPHL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+KBOT-1)*KSTR

      RC1     = 1.5
      RC2     = RC1 - 1.

C ... Calculate something on the mirror surfaces

      DO 1000 J = KY1,KY2 ! These should not be extended

      IJ      = (JN+J-1)*JSTR + KA
      NN      = (J-KY1)*(KX2-KX1+1) - KX1 + 1
*      NR = (LN+J-KY1)*(IDIMP+2*LN) - KX1 + LN + 1
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
      N1      = (JN+J-1)*ISTRID + IN
      DO 1000 I = KX1,KX2
         II      = (IN+I-1)*ISTR + IJ + 1 ! CELL INDEX (STRIDE=ISTR)
         II2     = II + LSTRID            ! 2. comp cell
         IF1 = II  + NFL                     ! FLUX INDEX
         IF2 = IF1 + LSTRID                  ! 2. face
         NI      = I  + N1                ! PATCH VELOCITY INDEX
         NP  = I   + NR                      ! Patch index with LN ghost cells

         ROM     = .5*(RO(II) + RO(II-LSTRID))
         POM     = .5*(P(II)  + P(II-LSTRID))
         PDM     = .5*(PDI(II)+ PDI(II-LSTRID))
         RMM     = .5*(RM(II) + RM(II-LSTRID))
         RNN     = .5*(RN(II) + RN(II-LSTRID))
         RWW     = .5*(RW(II) + RW(II-LSTRID))
          
         EKIN    = .5*(RMM**2 + RNN**2 + RWW**2)/ROM

      IF(M2 == 1) THEN ! Store the primitive variables

         IF(INLOUT == 0) THEN ! Store i/o data
          WMFLUX(NP)= F3R(IF1) /(A3(IF1) + 1.E-20)
          SURFX(NP) = RMM/ROM
          SURFY(NP) = RNN/ROM
          SURFZ(NP) = RWW/ROM
c          PB4       = MAX(-.99*FRSPRE, RC1*PDI(II) - RC2*PDI(II2))
          PB4       = MAX(-.99*FRSPRE, PDM)
          CPWALL(NP)= PB4
          TWALL(NP) = RC1*TEMP(II) - RC2*TEMP(II2)
         ENDIF ! Store i/o data

      ENDIF ! M2 == 1

1000  CONTINUE

      RETURN
      END SUBROUTINE BOUMIR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C ******************** SECOND-MOMENT CLOSURE FLUX  *****************

      SUBROUTINE FLUXRE(FRO,FRM,FRN,FRW,FE,FRK,FEPS,VOL,A,RO,RM,RN,RW,
     2 RK,REPS,P,U,V,W,E,C,IMAX,JMAX,KMAX,ICYCLE,INTEM,INTET,IDI,VIS,
     3 EPS2,VIST,PDI,A2XA,A2YA,A2ZA,PRLAM,PRT,
     4 UROT,TEMP,CH,DRDH,DRDP,ISTATE,GAMMA,E0REF,T0REF,
     5 FRSDEN,FRSPRE,FRSVEL,HAT1,HAT2,HAT3,HAT4,
     6 FUU,FUV,FUW,FVV,FVW,FWW,RUU,RUV,RUW,RVV,RVW,RWW,
     7 PSIGSC,PSIGS2,NSCAL,MAXSB,KMAXP1,KSTR,IDIR,ZZZ,MAW,MAXW,RKSI,
     8 INCHIML)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),RK(*),REPS(*),FRK(*)
     4  ,FEPS(*),UROT(*),TEMP(*),CH(*),DRDH(*),DRDP(*),HAT1(*)
     5  ,HAT2(*),HAT3(*),HAT4(*),PDI(*)
     6  ,FUU(*),FUV(*),FUW(*),FVV(*),FVW(*),FWW(*)
     7  ,RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),PSIGSC(*),PSIGS2(*)
     8  ,RKSI(*)

      REAL, POINTER ::   
     +     ROIP(:),RNIP(:),RMIP(:),EIP(:),RWIP(:),RKIP(:),
     +     ROIM(:),RNIM(:),RMIM(:),EIM(:),RWIM(:),RKIM(:),
     +      PMI(:), PPI(:),REIM(:),REIP(:),ERL(:), ELR(:),DEDPH(:),
     +     DEDRH(:),
     +     RUUP(:),RUVP(:),RUWP(:),RVVP(:),RVWP(:),RWWP(:),
     +     RUUM(:),RUVM(:),RUWM(:),RVVM(:),RVWM(:),RWWM(:)
      REAL, TARGET ::  ZZZ(MAXW)
      LOGICAL INCHIML

      ROIP => ZZZ( 0*MAW+1: 1*MAW);RNIP => ZZZ( 1*MAW+1: 2*MAW)
      RMIP => ZZZ( 2*MAW+1: 3*MAW);EIP  => ZZZ( 3*MAW+1: 4*MAW)
      ROIM => ZZZ( 4*MAW+1: 5*MAW);RNIM => ZZZ( 5*MAW+1: 6*MAW) 
      RMIM => ZZZ( 6*MAW+1: 7*MAW);EIM  => ZZZ( 7*MAW+1: 8*MAW) 
      RWIM => ZZZ( 8*MAW+1: 9*MAW);RWIP => ZZZ( 9*MAW+1:10*MAW)
      PMI  => ZZZ(10*MAW+1:11*MAW); PPI => ZZZ(11*MAW+1:12*MAW)
      RKIM => ZZZ(12*MAW+1:13*MAW);RKIP => ZZZ(13*MAW+1:14*MAW) 
      REIM => ZZZ(14*MAW+1:15*MAW);REIP => ZZZ(15*MAW+1:16*MAW)
      ERL  => ZZZ(16*MAW+1:17*MAW); ELR => ZZZ(17*MAW+1:18*MAW)
      DEDRH=> ZZZ(18*MAW+1:19*MAW);DEDPH=> ZZZ(19*MAW+1:20*MAW)
      RUUP => ZZZ(20*MAW+1:21*MAW);RUUM => ZZZ(21*MAW+1:22*MAW)
      RUVP => ZZZ(22*MAW+1:23*MAW);RUVM => ZZZ(23*MAW+1:24*MAW)
      RUWP => ZZZ(24*MAW+1:25*MAW);RUWM => ZZZ(25*MAW+1:26*MAW)
      RVVP => ZZZ(26*MAW+1:27*MAW);RVVM => ZZZ(27*MAW+1:28*MAW)
      RVWP => ZZZ(28*MAW+1:29*MAW);RVWM => ZZZ(29*MAW+1:30*MAW)
      RWWP => ZZZ(30*MAW+1:31*MAW);RWWM => ZZZ(31*MAW+1:32*MAW)


C ... 13 AND 14 ARE NOT USED YET PPR 3.11.1994

      PSI33   = 1./1.3
      EPS     = 1.E-10

C ... SECOND-ORDER UPWIND IS FORCED FOR TURBULENCE QUANTITIES
C              RK3   = 0.
C              RK4   = 2.
C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW
C

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES
C
      CALL XXTRAP(  RO,ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   U,RMIP,RMIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   V,RNIP,RNIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   W,RWIP,RWIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(   P, PPI, PMI,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)
      CALL XXTRAP(  RK,RKIP,RKIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP(REPS,REIP,REIM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RUU,RUUP,RUUM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RUV,RUVP,RUVM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RUW,RUWP,RUWM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RVV,RVVP,RVVM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RVW,RVWP,RVWM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)
      CALL XXTRAP( RWW,RWWP,RWWM,IJTRID,IL,KL,IA,INTET,RKSI,INCHIML)

      DO  500 IG= 1,IJTRID
         I       = IA + IG
         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         PPI(IG)     = MAX(0.001*FRSPRE,PPI(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         PMI(IG)     = MAX(0.001*FRSPRE,PMI(IG))
500   CONTINUE

C ... SPECIFIC INTERNAL ENERGY AS A FUNCTION OF PRESSURE AND DENSITY

      CALL EFPRO(EIP,PPI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(EIM,PMI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ELR,PPI,ROIM,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)
      CALL EFPRO(ERL,PMI,ROIP,IJTRID,ISTATE,GAMMA,0.,E0REF,T0REF)

C ... DERIVATIVES ACCORDING TO GLAISTER

      I1     = IA + 1
      CALL DERIV(RO(I1),P(I1),DEDPH,DEDRH,DRDP(I1),DRDH(I1),
     + EIP,EIM,ELR,ERL,PPI,PMI,ROIP,ROIM,IJTRID,IL)

      DO  600 IG= 1,IJTRID
         I       = IA + IG

         RMIPI    = RMIP(IG)
         RNIPI    = RNIP(IG)
         RWIPI    = RWIP(IG)
         RMIMI    = RMIM(IG)
         RNIMI    = RNIM(IG)
         RWIMI    = RWIM(IG)

C ... CONTRAVARIANT MASS FLUXES (RMIP,UIP) AND PRESSURES

         RKIP(IG)    = MAX(1.E-12,RKIP(IG))
         REIP(IG)    = MAX(1.E-12,REIP(IG))
         RKIM(IG)    = MAX(1.E-12,RKIM(IG))
         REIM(IG)    = MAX(1.E-12,REIM(IG))

C ... LIMIT NORMAL REYNOLDS STRESSES

         U11      = RUUP(IG)
         U12      = RUVP(IG)
         U13      = RUWP(IG)
         U22      = RVVP(IG)
         U23      = RVWP(IG)
         U33      = RWWP(IG)

         V11      = MAX(EPS,RUUM(IG))
         V12      = RUVM(IG)
         V13      = RUWM(IG)
         V22      = MAX(EPS,RVVM(IG))
         V23      = RVWM(IG)
         V33      = MAX(EPS,RWWM(IG))

C ... FROM THE PRIMITIVE TO THE CONTRAVARIANT VARIABLES

         A11      = A2XA(I)
         A12      = A2YA(I)
         A13      = A2ZA(I)

         CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

*         A21      = A2ZA(I) - A2YA(I)
*         A22      = A2XA(I) - A2ZA(I)
*         A23      = A2YA(I) - A2XA(I)
*         PALPO    = A21**2 + A22**2 + A23**2
*         IF(PALPO < 1.E-3) THEN
*            A21   = A2ZA(I) - A2YA(I)
*            A22   = A2XA(I) + A2ZA(I)
*            A23   =-A2YA(I) - A2XA(I)
*            PALPO = A21**2 + A22**2 + A23**2
*         ENDIF
*         PALPO    = 1./SQRT(PALPO)
*         A21      = A21*PALPO
*         A22      = A22*PALPO
*         A23      = A23*PALPO
*         A31      = A2YA(I)*A23 - A2ZA(I)*A22
*         A32      = A2ZA(I)*A21 - A2XA(I)*A23
*         A33      = A2XA(I)*A22 - A2YA(I)*A21

         RMIP(IG) = A11*RMIPI + A12*RNIPI + A13*RWIPI
         RNIP(IG) = A21*RMIPI + A22*RNIPI + A23*RWIPI
         RWIP(IG) = A31*RMIPI + A32*RNIPI + A33*RWIPI
         RMIM(IG) = A11*RMIMI + A12*RNIMI + A13*RWIMI
         RNIM(IG) = A21*RMIMI + A22*RNIMI + A23*RWIMI
         RWIM(IG) = A31*RMIMI + A32*RNIMI + A33*RWIMI

C ... ROTATE REYNOLDS STRESS TENSOR 

         B11 = A11*U11 + A12*U12 + A13*U13
         B12 = A11*U12 + A12*U22 + A13*U23
         B13 = A11*U13 + A12*U23 + A13*U33
         B21 = A21*U11 + A22*U12 + A23*U13
         B22 = A21*U12 + A22*U22 + A23*U23
         B23 = A21*U13 + A22*U23 + A23*U33
         B31 = A31*U11 + A32*U12 + A33*U13
         B32 = A31*U12 + A32*U22 + A33*U23
         B33 = A31*U13 + A32*U23 + A33*U33
               
         F11 = A11*V11 + A12*V12 + A13*V13
         F12 = A11*V12 + A12*V22 + A13*V23
         F13 = A11*V13 + A12*V23 + A13*V33
         F21 = A21*V11 + A22*V12 + A23*V13
         F22 = A21*V12 + A22*V22 + A23*V23
         F23 = A21*V13 + A22*V23 + A23*V33
         F31 = A31*V11 + A32*V12 + A33*V13
         F32 = A31*V12 + A32*V22 + A33*V23
         F33 = A31*V13 + A32*V23 + A33*V33
               
         C11 = B11*A11 + B12*A12 + B13*A13
         C12 = B11*A21 + B12*A22 + B13*A23
         C13 = B11*A31 + B12*A32 + B13*A33
         C22 = B21*A21 + B22*A22 + B23*A23
         C23 = B21*A31 + B22*A32 + B23*A33
         C33 = B31*A31 + B32*A32 + B33*A33
                                          
         G11 = F11*A11 + F12*A12 + F13*A13
         G12 = F11*A21 + F12*A22 + F13*A23
         G13 = F11*A31 + F12*A32 + F13*A33
         G22 = F21*A21 + F22*A22 + F23*A23
         G23 = F21*A31 + F22*A32 + F23*A33
         G33 = F31*A31 + F32*A32 + F33*A33


         RUUP(IG)    = MAX(EPS,C11)
         RUVP(IG)    = C12
         RUWP(IG)    = C13
         RVVP(IG)    = MAX(EPS,C22)
         RVWP(IG)    = C23
         RWWP(IG)    = MAX(EPS,C33)

         RUUM(IG)    = MAX(EPS,G11)
         RUVM(IG)    = G12
         RUWM(IG)    = G13
         RVVM(IG)    = MAX(EPS,G22)
         RVWM(IG)    = G23
         RWWM(IG)    = MAX(EPS,G33)

C*******************************************************************
C ... ROE'S FLUX SPLITTING                                         *
C*******************************************************************

C ... THERMODYNAMIC VARIABLES

      DPDE   =  1./DEDPH(IG)
      DPDRO  = -DEDRH(IG)/DEDPH(IG)

      YPROR  = 1./ROIP(IG)
      UR     = RMIP(IG)
      VR     = RNIP(IG)
      WR     = RWIP(IG)
      RKR    = RKIP(IG)*YPROR
      RUUR   = RUUP(IG)*YPROR
      RUVR   = RUVP(IG)*YPROR
      RUWR   = RUWP(IG)*YPROR
      RVVR   = RVVP(IG)*YPROR
      RVWR   = RVWP(IG)*YPROR
      RWWR   = RWWP(IG)*YPROR


      YPROL  = 1./ROIM(IG)
      UL     = RMIM(IG)
      VL     = RNIM(IG)
      WL     = RWIM(IG)
      RKL    = RKIM(IG)*YPROL
      RUUL   = RUUM(IG)*YPROL
      RUVL   = RUVM(IG)*YPROL
      RUWL   = RUWM(IG)*YPROL
      RVVL   = RVVM(IG)*YPROL
      RVWL   = RVWM(IG)*YPROL
      RWWL   = RWWM(IG)*YPROL

C ... TOTAL INTERNAL ENERGY

      EIR    = EIP(IG) +.5*(UR**2 + VR**2 + WR**2 +RUUR+RVVR+RWWR)
      EIL    = EIM(IG) +.5*(UL**2 + VL**2 + WL**2 +RUUL+RVVL+RWWL)

      PR     = PPI(IG)
      PKR    = .6667*RKIP(IG)
      EPSR   = REIP(IG)*YPROR
      HR     = EIR + PR*YPROR
      ER     = EIP(IG)
      ROHAT  = SQRT(ROIP(IG)*ROIM(IG))
      PRUR   = RUUP(IG)
      PRVR   = RUVP(IG)
      PRWR   = RUWP(IG)

      PL     = PMI(IG)
      PKL    = .6667*RKIM(IG)
      EPSL   = REIM(IG)*YPROL
      HL     = EIL + YPROL*PL
      EL     = EIM(IG)
      PRUL   = RUUM(IG)
      PRVL   = RUVM(IG)
      PRWL   = RUWM(IG)

      SQRL   = 1./(ROIM(IG) + ROHAT)
      CQL    = ROIM(IG)*SQRL
      CQR    = 1.- CQL

      UHAT   = CQR*UR  + CQL*UL
      VHAT   = CQR*VR  + CQL*VL
      WHAT   = CQR*WR  + CQL*WL
      HHAT   = CQR*HR  + CQL*HL
      EHAT   = CQR*ER  + CQL*EL
      RKHAT  = CQR*RKR + CQL*RKL
      EPSHAT = CQR*EPSR+ CQL*EPSL
      RUUHAT = CQR*RUUR+ CQL*RUUL
      RUVHAT = CQR*RUVR+ CQL*RUVL
      RUWHAT = CQR*RUWR+ CQL*RUWL
      RVVHAT = CQR*RVVR+ CQL*RVVL
      RVWHAT = CQR*RVWR+ CQL*RVWL
      RWWHAT = CQR*RWWR+ CQL*RWWL
      PHAT   = MAX(ROHAT*(HHAT - EHAT -.5*(UHAT**2+VHAT**2+WHAT**2)
     +         - RKHAT),0.001*FRSPRE)
C *********************************************************
C     GAMMA   = 1.4
C     GM1     = GAMMA - 1.
C     CHAT    = SQRT(GAMMA*(PHAT/ROHAT+.6667*RKHAT))
C     PDPDE   = 1./(ROHAT*GM1)
C *********************************************************
C      APPP   = (PHAT/ROHAT - ROHAT*DEDRH(IG) +
C     +        .6667*RKHAT*(1./(ROHAT+EPS)+DEDPH(IG)))
C     +        /(ROHAT*DEDPH(IG))
      APPP   = (DPDE*PHAT/ROHAT**2 + DPDRO + 3*RUUHAT)
      CMIN   = MIN(C(I),C(I-IL))
      CMAX   = MAX(C(I),C(I-IL))
c     APPP   = MIN(CMAX**2,APPP)
c     APPP   = MAX(CMIN**2,APPP)    ! EVALUOITAVA
      CHAT   = SQRT(MAX(APPP,EPS))
      PDPDE  = DEDPH(IG)

C ... LIMITERS FOR RATIOS OF U1J OVER SQRT(U11)

      DELUV  = RUVR - RUVL
      DELUW  = RUWR - RUWL
      SU11   = SQRT(MAX(RUUHAT,1.E-5))
      VVLIM  = 2.*SU11
      WWLIM  = 2.*SU11
      RU2Y1  = SIGN(MIN(ABS(RUVHAT)/SU11,SQRT(RVVHAT)),RUVHAT)
      RU3Y1  = SIGN(MIN(ABS(RUWHAT)/SU11,SQRT(RWWHAT)),RUWHAT)
      DU2Y1  = SIGN(MIN(ABS(DELUV)/SU11,VVLIM),DELUV)
      DU3Y1  = SIGN(MIN(ABS(DELUW)/SU11,WWLIM),DELUW)

C ... MULTIPLIERS (CHARACTERISTIC VARIABLES)

      DRO    = ROIP(IG)-ROIM(IG)
      DU     = UR - UL
      DV     = VR - VL
      DW     = WR - WL
      DK     = RKR - RKL
C      DP     = PR - PL + .6667*(RKHAT*DRO + ROHAT*DK)
      DE     = ER - EL
      DUU    = RUUR - RUUL
      PCHAT2 = 1./CHAT**2
      PCUUT2 = 1./(CHAT**2 - RUUHAT)
      PRO2   = PHAT/ROHAT**2
      APU1   = (.5*(DPDE*PRO2-DPDRO)*DRO - DPDE*DE - ROHAT*DUU)*PCUUT2
      APU2   = 2.*((-DPDE*PRO2+DPDRO-RUUHAT)*DRO +
     +     2.*DPDE*DE + 2.*ROHAT*DUU)*PCHAT2
      ALFA1  = ((DPDE*PRO2+2.*RUUHAT)*DRO - DPDE*DE - ROHAT*DUU)*PCHAT2
      ALFA2  = .5*((DPDRO+RUUHAT)*DRO + ROHAT*CHAT*DU + DPDE*DE + 
     +     ROHAT*DUU)*PCHAT2
      ALFA3  = APU1*RU2Y1 + ROHAT*RUVHAT*PCUUT2*DU-.5*ROHAT*DV
     +    +.5*ROHAT*DU2Y1
      ALFA4  = APU1*RU3Y1 + ROHAT*RUWHAT*PCUUT2*DU-.5*ROHAT*DW
     +    +.5*ROHAT*DU3Y1
      ALFA5  = ALFA2  - ROHAT/CHAT*DU
      ALFA6  = ROHAT*(RKR - RKL)
      ALFA7  = ROHAT*(EPSR - EPSL)
      ALFA8  = 2.*(-RUUHAT*(DPDRO+RUUHAT)*DRO - DPDE*RUUHAT*DE +
     +     ROHAT*(.5*CHAT**2-RUUHAT)*DUU)*PCHAT2
      ALFA9  = ALFA3 - 2.*ROHAT*RUVHAT*PCUUT2*DU + ROHAT*DV
      ALFA10 = ALFA4 - 2.*ROHAT*RUWHAT*PCUUT2*DU + ROHAT*DW
      ALFA11 = APU2*RU2Y1**2 - 2.*ROHAT*RU2Y1*DU2Y1 + ROHAT*(RVVR-RVVL)
      ALFA12 = APU2*RU2Y1*RU3Y1 - ROHAT*RU3Y1*DU2Y1 - 
     +                            ROHAT*RU2Y1*DU3Y1 + ROHAT*(RVWR-RVWL)
      ALFA13 = APU2*RU3Y1**2 - 2.*ROHAT*RU3Y1*DU3Y1 + ROHAT*(RWWR-RWWL)

C ... ABSOLUTE VALUES OF THE CHARACTERISTIC SPEEDS

      UC     = UROT(I)
      RL1    = ABS(UHAT-UC)
      RL2    = ABS(UHAT-UC + CHAT)
      RL3    = ABS(UHAT-UC-SU11)
      RL4    = RL3
      RL5    = ABS(UHAT-UC - CHAT)
      RL6    = RL1
      RL7    = RL1
      RL8    = RL1
      RL9    = ABS(UHAT-UC+SU11)
      RL10   = RL9
      RL11   = RL1
      RL12   = RL1
      RL13   = RL1

      ALF1   = RL1*ALFA1
      ALF2   = RL2*ALFA2
      ALF3   = RL3*ALFA3
      ALF4   = RL4*ALFA4
      ALF5   = RL5*ALFA5
      ALF6   = RL6*ALFA6
      ALF7   = RL7*ALFA7
      ALF8   = RL8*ALFA8
      ALF9   = RL9*ALFA9
      ALF10  = RL10*ALFA10
      ALF11  = RL11*ALFA11
      ALF12  = RL12*ALFA12
      ALF13  = RL13*ALFA13

C ... CONVECTIVE FLUXES

      ROUCR  = ROIP(IG)*(UR - UC)
      ROUCL  = ROIM(IG)*(UL - UC)
      DAMP   =  ALF2 + ALF5 + ALF1
C ... AVERAGED FLUXES
      F1PI   = ROUCR + ROUCL 
      F2PI   = ROUCR*UR+PR+PRUR + ROUCL*UL+PL+PRUL 
      F3PI   = ROUCR*VR   +PRVR + ROUCL*VL   +PRVL
      F4PI   = ROUCR*WR   +PRWR + ROUCL*WL   +PRWL
      F5PI   = ROUCR*(HR+RUUR)  + VL*RUVL + WL*RUWL +  
     +         ROUCL*(HL+RUUL)  + VR*RUVR + WR*RUWR
     +         + UC*(PR+PRUR + PL+PRUL) 
      F6PI   = .5*A(I)*(ROUCR*RKR + ROUCL*RKL - DAMP*RKHAT - ALF6)
      F7PI   = .5*A(I)*(ROUCR*EPSR+ROUCL*EPSL - DAMP*EPSHAT- ALF7)
      F8PI   = ROUCR*RUUR+ROUCL*RUUL
      F9PI   = ROUCR*RUVR+ROUCL*RUVL
      F10PI  = ROUCR*RUWR+ROUCL*RUWL
      F11PI  = ROUCR*RVVR+ROUCL*RVVL
      F12PI  = ROUCR*RVWR+ROUCL*RVWL
      F13PI  = ROUCR*RWWR+ROUCL*RWWL

C ... DAMPING TERMS

      F1PIG  = .5*A(I)*(F1PI - DAMP)
      F2PIG  = .5*A(I)*(F2PI - (DAMP*UHAT+(ALF2-ALF5)*CHAT))
      F3PIG  = .5*A(I)*(F3PI - (DAMP*VHAT+
     +     2.*CHAT*RUVHAT*PCUUT2*(ALF2-ALF5) - ALF3 + ALF9))
      F4PIG  = .5*A(I)*(F4PI - (DAMP*WHAT+
     +     2.*CHAT*RUWHAT*PCUUT2*(ALF2-ALF5) - ALF4 + ALF10))
      F5PIG  = .5*A(I)*(F5PI - (DAMP*HHAT - 
     +     ROHAT*(CHAT**2-2.*RUUHAT)/DPDE*ALF1 +(ALF2-ALF5)*UHAT*CHAT+
     +     (ALF2+ALF5)*(RUUHAT + 2.*PCUUT2*
     +     (RUVHAT**2+RUWHAT**2+CHAT*VHAT*RUVHAT+CHAT*WHAT*RUWHAT)) +
     +     ALF3*(RU2Y1-VHAT) + ALF4*(RU3Y1-WHAT) + 
     +     ALF8*(.5-ROHAT/DPDE) + ALF9*(RU2Y1+VHAT)+ALF10*(RU3Y1+WHAT)+
     +     ALF11*.5 + ALF13*.5))
      F6PIG  = F6PI
      F7PIG  = F7PI
      F8PIG  = .5*A(I)*(F8PI - ((DAMP+2.*(ALF2+ALF5))*RUUHAT + 
     +     ALF8))
      F9PIG  = .5*A(I)*(F9PI - (2.*CHAT**2*PCUUT2*RUVHAT*(ALF2+ALF5) +
     +    SU11*(ALF3 + ALF9)))
      F10PIG = .5*A(I)*(F10PI- (2.*CHAT**2*PCUUT2*RUWHAT*(ALF2+ALF5) +
     +    SU11*(ALF4 + ALF10)))
      F11PIG = .5*A(I)*(F11PI - (DAMP*RVVHAT+
     +     4.*PCUUT2*RUVHAT**2*    (ALF2+ALF5) +
     +     2.*RU2Y1*(ALF3+ALF9)                     + ALF11))
      F12PIG = .5*A(I)*(F12PI - (DAMP*RVWHAT+
     +     4.*PCUUT2*RUVHAT*RUWHAT*(ALF2+ALF5) +
     +        RU3Y1*(ALF3+ALF9) + RU2Y1*(ALF4+ALF10)+ ALF12))
      F13PIG = .5*A(I)*(F13PI - (DAMP*RWWHAT+
     +     4.*PCUUT2*RUWHAT**2*    (ALF2+ALF5) +
     +     2.*RU3Y1*(ALF4+ALF10)                    + ALF13))
C*    HATS ARE SAVED FOR SCALAR TRANSPORT SUBROUTINE PPR 10.2.
      HAT1(I) = DAMP
      HAT2(I) = ROUCR
      HAT3(I) = ROUCL
      HAT4(I) = RL1

c      F1P  = .5*A(I)*( - DAMP)
c      F2P  = .5*A(I)*( - (DAMP*UHAT+(ALF2-ALF5)*CHAT))
c      F3P  = .5*A(I)*( - (DAMP*VHAT+
c     +     2.*CHAT*RUVHAT*PCUUT2*(ALF2-ALF5) - ALF3
c     +     + ALF9))
c      F4P  = .5*A(I)*( - (DAMP*WHAT+
c     +     2.*CHAT*RUWHAT*PCUUT2*(ALF2-ALF5) - ALF4
c     +     + ALF10))
c      F5P  = .5*A(I)*( - (DAMP*HHAT - 
c     +     ROHAT*(CHAT**2-2.*RUUHAT)/DPDE*ALF1 +(ALF2-ALF5)*UHAT*CHAT+
c     +     (ALF2+ALF5)*(RUUHAT + 2.*PCUUT2*
c     +     (RUVHAT**2+RUWHAT**2+CHAT*VHAT*RUVHAT+CHAT*WHAT*RUWHAT)) +
c     +     ALF3*(RU2Y1-VHAT) + ALF4*(RU3Y1-WHAT) + 
c     +     ALF8*(.5-ROHAT/DPDE) + ALF9*(RU2Y1+VHAT)+ALF10*(RU3Y1+WHAT)+
c     +     ALF11*.5 + ALF13*.5))
c      F6P  = F6PI
c      F7P  = F7PI
c      F8P  = .5*A(I)*(- ((DAMP+2.*(ALF2+ALF5))*RUUHAT + 
c     +     ALF8))
c      F9P  = .5*A(I)*(- (2.*CHAT**2*PCUUT2*RUVHAT*(ALF2+ALF5) +
c     +    SU11*ALF3 +  SU11*ALF9))
c      F10P = .5*A(I)*(- (2.*CHAT**2*PCUUT2*RUWHAT*(ALF2+ALF5) +
c     +    SU11*ALF4 +  SU11*ALF10))
c      F11P = .5*A(I)*(- (DAMP*RVVHAT+
c     +     4.*PCUUT2*RUVHAT**2*    (ALF2+ALF5) +
c     +     2.*RU2Y1*(ALF3+ALF9)                     + ALF11))
c      F12P = .5*A(I)*(- (DAMP*RVWHAT+
c     +     4.*PCUUT2*RUVHAT*RUWHAT*(ALF2+ALF5) +
c     +        RU3Y1*(ALF3+ALF9) + RU2Y1*(ALF4+ALF10)+ ALF12))
c      F13P = .5*A(I)*(- (DAMP*RWWHAT+
c     +     4.*PCUUT2*RUWHAT**2*    (ALF2+ALF5) +
c     +     2.*RU3Y1*(ALF4+ALF10)                    + ALF13))

c      if (abs(f1pi*.5*a(i)) <= abs(f1p) )
c     +     write(*,*) 'f1 ',f1pi*.5*a(i),f1p,f1pi*.5*a(i)/f1p
c      if (abs(f2pi*.5*a(i)) <= abs(f2p) )
c     +     write(*,*) 'f2 ',f2pi*.5*a(i),f2p,f2pi*.5*a(i)/f2p
c      if (abs(f3pi*.5*a(i)) <= abs(f3p) )
c     +     write(*,*) 'f3 ',f3pi*.5*a(i),f3p,f3pi*.5*a(i)/f3p
c      if (abs(f4pi*.5*a(i)) <= abs(f4p) )
c     +     write(*,*) 'f4 ',f4pi*.5*a(i),f4p,f4pi*.5*a(i)/f4p
c      if (abs(f5pi*.5*a(i)) <= abs(f5p) )
c     +     write(*,*) 'f5 ',f5pi*.5*a(i),f5p,f5pi*.5*a(i)/f5p
c      if (abs(f8pi*.5*a(i)) <= abs(f8p) )
c     +     write(*,*) 'f8 ',f8pi*.5*a(i),f8p,f8pi*.5*a(i)/f8p
c      if (abs(f9pi*.5*a(i)) <= abs(f9p) )
c     +     write(*,*) 'f9 ',f9pi*.5*a(i),f9p,f9pi*.5*a(i)/f9p
c      if (abs(f10pi*.5*a(i)) <= abs(f10p) )
c     +     write(*,*) 'f10',f10pi*.5*a(i),f10p,f10pi*.5*a(i)/f10p
c      if (abs(f11pi*.5*a(i)) <= abs(f11p) )
c     +     write(*,*) 'f11',f11pi*.5*a(i),f11p,f11pi*.5*a(i)/f11p
c      if (abs(f12pi*.5*a(i)) <= abs(f12p) )
c     +     write(*,*) 'f12',f12pi*.5*a(i),f12p,f12pi*.5*a(i)/f12p
c      if (abs(f13pi*.5*a(i)) <= abs(f13p)) 
c     +     write(*,*) 'f13',f13pi*.5*a(i),f13p,f13pi*.5*a(i)/f13p


C ... ROTATE THE COORDINATE SYSTEM
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES

      F2PM    = A11*F2PIG + A21*F3PIG + A31*F4PIG
      F3PM    = A12*F2PIG + A22*F3PIG + A32*F4PIG
      F4PM    = A13*F2PIG + A23*F3PIG + A33*F4PIG

C ... ROTATE REYNOLDS STRESS TENSOR 
      C11     = F8PIG
      C12     = F9PIG
      C13     = F10PIG
      C22     = F11PIG
      C23     = F12PIG
      C33     = F13PIG
C ... STRESS TENSOR IS SYMMETRIC
      C21 = C12
      C31 = C13
      C32 = C23
      D11 = A11*C11 + A21*C21 + A31*C31
      D12 = A11*C12 + A21*C22 + A31*C32
      D13 = A11*C13 + A21*C23 + A31*C33
      D21 = A12*C11 + A22*C21 + A32*C31
      D22 = A12*C12 + A22*C22 + A32*C32
      D23 = A12*C13 + A22*C23 + A32*C33
      D31 = A13*C11 + A23*C21 + A33*C31
      D32 = A13*C12 + A23*C22 + A33*C32
      D33 = A13*C13 + A23*C23 + A33*C33

      E11 = D11*A11 + D12*A21 + D13*A31
      E12 = D11*A12 + D12*A22 + D13*A32
      E13 = D11*A13 + D12*A23 + D13*A33
      E22 = D21*A12 + D22*A22 + D23*A32
      E23 = D21*A13 + D22*A23 + D23*A33
      E33 = D31*A13 + D32*A23 + D33*A33

C
C ... FLUXES
C
      FRO(I)  = F1PIG  + FRO(I)
      FRM(I)  = F2PM   + FRM(I)
      FRN(I)  = F3PM   + FRN(I)
      FRW(I)  = F4PM   + FRW(I)
      FE(I)   = F5PIG  + FE(I)
      FRK(I)  = F6PIG  + FRK(I)
      FEPS(I) = F7PIG  + FEPS(I)

      FUU(I)  = E11  + FUU(I)
      FUV(I)  = E12  + FUV(I)
      FUW(I)  = E13  + FUW(I)
      FVV(I)  = E22  + FVV(I)
      FVW(I)  = E23  + FVW(I)
      FWW(I)  = E33  + FWW(I)
600   CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE FLUXRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXVV(FRO,FRM,FRN,FRW,FE,VOL,A,RO,RM,RN,RW,P,U,V,W,E,
     2 C,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,EPS2,VIST,A2XA,A2YA,
     3 A2ZA,PRLAM,PRT,UROT,TEMP,CH,
     4 DRDH,DRDP,ISTATE,GAMMA,FRSDEN,FRSPRE,HAT1,HAT2,HAT3,HAT4,
     5 KSTR,IDIR)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: A(*),RO(*),RN(*),RW(*),U(*),E(*),FRO(*),FRN(*),FE(*)
     2  ,FRM(*),V(*),RM(*),C(*),A2XA(*),A2YA(*),VIS(*),EPS2(*),VIST(*)
     3  ,VOL(*),A2ZA(*),FRW(*),W(*),P(*),UROT(*),TEMP(*),CH(*)
     4  ,DRDP(*),DRDH(*),HAT1(*),HAT2(*),HAT3(*),HAT4(*)
      LOGICAL :: INVIS

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID
      IA      = ILL*KN
 
      PRS     = PRLAM/PRT
      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

C
C ... THIN-LAYER VISCOUS FLUXES ARE CALCULATED WITH REYNOLDS STRESSES
C
      IF(IDI /= 0) THEN
      DO 900 KG = 1,KMAXID
         IA      = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = IA + (JG+JN-1)*ISTRID + IN
         DO 850 IG = 1,IMAXID
            I      = IA2 + IG
            YPD2    = A(I)/(VOL(I) + VOL(I-IL))
            SURG    = -A(I)*YPD2
            SURF    = SURG*(VIS(I) + VIS(I-IL))             ! LAMINAR
            SURFE    = SURG*((1.+PRS*VIST(I)/VIS(I))*CH(I) + 
     +      (1.+PRS*VIST(I-IL)/VIS(I-IL))* CH(I-IL))        ! BOUSSINESQ
            
            DUC     = .333333*(A2XA(I)*(U(I) - U(I-IL))   +
     +           A2YA(I)*(V(I) - V(I-IL)) + A2ZA(I)*(W(I) - W(I-IL)))

            F2VIG   = SURF*(U(I)-U(I-IL) + A2XA(I)*DUC)
            F3VIG   = SURF*(V(I)-V(I-IL) + A2YA(I)*DUC)
            F4VIG   = SURF*(W(I)-W(I-IL) + A2ZA(I)*DUC)
            F5VIG   = .5*((U(I)+U(I-IL))*F2VIG+(V(I)+V(I-IL))*F3VIG+
     +           (W(I)+W(I-IL))*F4VIG) + SURFE*(TEMP(I) - TEMP(I-IL))

            FRO(I)  = 0.
            FRM(I)  = F2VIG
            FRN(I)  = F3VIG
            FRW(I)  = F4VIG
            FE(I)   = FE(I) + F5VIG
C                     =====
C ... POSSIBLE DIFFUSION OF TURBULENT KINECTIC ENERGY PPR 10.8.1995
C
 850     CONTINUE
 900  CONTINUE
      ENDIF

      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
         FRO(I)  = 0.
         FRM(I)  = 0.
         FRN(I)  = 0.
         FRW(I)  = 0.
         FE(I)   = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLUXVV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXIS(FRO,FRM,FRN,FRW,FE,RO,P,U,V,W,E,
     2     RUU,RUV,RUW,RVV,RVW,RWW,
     3     VOL,A,AXA,AYA,AZA,INTERK,IDI,IMAX,JMAX,KMAX,
     4     KMAXP1,KSTR)
      
      USE NS3CO, ONLY : IN, JN, KN
      
C ... This subroutine corrects fluxes of the k-epsilon model by
C ... anisotrophy of the Reynold's stresses. PPR 20.4.1995

      DIMENSION :: FRO(*),FRM(*),FRN(*),FRW(*),FE(*),
     2     RO(*),P(*),U(*),V(*),W(*),E(*),
     4     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     5     VOL(*),A(*),AXA(*),AYA(*),AZA(*)

      EPS     = 1.E-10
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN

      IL      = KSTR
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)

      DO 1001 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
      DO 851 IG  = 1,IJTRID
         I       = IA + IG
         SURG    =  .5*A(I)
         U2      = U(I) + U(I-IL)
         V2      = V(I) + V(I-IL)
         W2      = W(I) + W(I-IL)
         FUU     = RUU(I) + RUU(I-IL)
         FUV     = RUV(I) + RUV(I-IL)
         FUW     = RUW(I) + RUW(I-IL)
         FVV     = RVV(I) + RVV(I-IL)
         FVW     = RVW(I) + RVW(I-IL)
         FWW     = RWW(I) + RWW(I-IL)
         XC      = -.3333333*(FUU + FVV + FWW)
         F1M     = SURG*(AXA(I)*(FUU + XC) + AYA(I)*FUV + AZA(I)*FUW)
         F1N     = SURG*(AYA(I)*(FVV + XC) + AXA(I)*FUV + AZA(I)*FVW)
         F1W     = SURG*(AZA(I)*(FWW + XC) + AXA(I)*FUW + AYA(I)*FVW)
         FRM(I)  = FRM(I) + F1M 
         FRN(I)  = FRN(I) + F1N 
         FRW(I)  = FRW(I) + F1W
         FE(I)   = FE(I)  + .5*(U2*F1M + V2*F1N + W2*F1W)
 851  CONTINUE
      
 1001 CONTINUE

      RETURN
      END SUBROUTINE FLUXIS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ***************** BOUNDARY ROUTINES **********************************

      SUBROUTINE FLUXIN(RO,U,V,W,E,A1,A1XA,A1YA,A1ZA,XC,YC,ZC,IMAX,
     + JMAX,KMAX,IN,JN,KN,IST,INITC)

      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN
      REAL :: RO(*),U(*),V(*),W(*),E(*),A1(*),A1XA(*),A1YA(*),A1ZA(*)
      REAL :: XC(*), YC(*), ZC(*)
C
C ... INITIALIZATION BASED ON MASS BALANCE IN A CHANNEL (CONSTANT RHO)
C     TOIMIIKS TAMA? EI, AINAKAAN KUNNOLLISESTI
C ... assumed constant density and changed to primitivi variables 6.1.00 PR
C
      EPS     = 1.E-20

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      ILPO    = (1-IST)/2*IL
      IL      = IL*IST

      IF(IST == 1) THEN
        KALKU = 1
        KOPPU = KMAX
      ELSE
        KALKU = KMAX
        KOPPU = 1
      ENDIF
      K       = KALKU
      KA      = (KN+K-1)*ISTRID*JSTRID

      AOS     = 0.
      ANI     = EPS
      RMASST  = 0.

      DO 700 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 700 I = 1,IMAX
      KK      = JJ + I

      RNX     = (XC(KK+IL) - XC(KK-IL))*.5*IST
      RNY     = (YC(KK+IL) - YC(KK-IL))*.5*IST
      RNZ     = (ZC(KK+IL) - ZC(KK-IL))*.5*IST
      RLEN    = 1./SQRT(RNX**2 + RNY**2 + RNZ**2 + EPS)
      RNX     = RNX*RLEN
      RNY     = RNY*RLEN
      RNZ     = RNZ*RLEN

C ... APPROXIMATE FLOW AREAS AND MASS FLOW AT THE CHANNEL INLET

      AX1     =  A1(KK+ILPO)* (RNX*A1XA(KK+ILPO) +
     +           RNY*A1YA(KK+ILPO) + RNZ*A1ZA(KK+ILPO))
      AX2     =  A1(KK+ILPO+IL)*(RNX*A1XA(KK+ILPO+IL) +
     +           RNY*A1YA(KK+ILPO+IL) + RNZ*A1ZA(KK+ILPO+IL))
      AOS     =  AOS + A1(KK+ILPO)* (RNX*A1XA(KK+ILPO) +
     +           RNY*A1YA(KK+ILPO) + RNZ*A1ZA(KK+ILPO))
      ANI     =  ANI + A1(KK+ILPO+IL)*(RNX*A1XA(KK+ILPO+IL) +
     +           RNY*A1YA(KK+ILPO+IL) + RNZ*A1ZA(KK+ILPO+IL))
      AREA    = .5*(AX1+AX1*AX1/AX2)
      RMAS2   = (A1XA(KK+ILPO)*U(KK-IL)
     +          +A1YA(KK+ILPO)*V(KK-IL) + A1ZA(KK+ILPO)*W(KK-IL))
      RMASST  = RMASST + RMAS2*AREA

 700  CONTINUE
C      SUHDE   = AOS/ANI

      DO 2000 K = KALKU,KOPPU,IST
      KA      = (KN+K-1)*ISTRID*JSTRID

      AOS     = 0.
      ANI     = EPS

      DO 1500 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 1500 I = 1,IMAX
      KK      = JJ + I
      RNX     = (XC(KK+IL) - XC(KK-IL))*.5*IST
      RNY     = (YC(KK+IL) - YC(KK-IL))*.5*IST
      RNZ     = (ZC(KK+IL) - ZC(KK-IL))*.5*IST
      RLEN    = 1./SQRT(RNX**2 + RNY**2 + RNZ**2 + EPS)
      RNX     = RNX*RLEN
      RNY     = RNY*RLEN
      RNZ     = RNZ*RLEN


C ... APPROXIMATE FLOW AREAS OF THE CHANNEL

      AOS     =  AOS + A1(KK+ILPO)* (RNX*A1XA(KK+ILPO) +
     +           RNY*A1YA(KK+ILPO) + RNZ*A1ZA(KK+ILPO))
      ANI     =  ANI + A1(KK+ILPO+IL)*(RNX*A1XA(KK+ILPO+IL) +
     +           RNY*A1YA(KK+ILPO+IL) + RNZ*A1ZA(KK+ILPO+IL))
 1500 CONTINUE
      AREAL   = (AOS+ANI)*.5
      RMASSF  = RMASST/AREAL

      DO 2000 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 2000 I = 1,IMAX
      KK      = JJ + I

      RNX     = (XC(KK+IL) - XC(KK-IL))*.5*IST
      RNY     = (YC(KK+IL) - YC(KK-IL))*.5*IST
      RNZ     = (ZC(KK+IL) - ZC(KK-IL))*.5*IST
      RLEN    = 1./SQRT(RNX**2 + RNY**2 + RNZ**2 + EPS)
      RNX     = RNX*RLEN
      RNY     = RNY*RLEN
      RNZ     = RNZ*RLEN


C ... MOMENTUM COMPONENTS AND ENERGY

      E(KK)   = E(KK) -.5*(U(KK)**2 + V(KK)**2 + W(KK)**2)*RO(KK)
      U(KK)  = RNX*RMASSF
      V(KK)  = RNY*RMASSF
      W(KK)  = RNZ*RMASSF
      E(KK)   = E(KK) +.5*(U(KK)**2 + V(KK)**2 + W(KK)**2)*RO(KK)
2000  CONTINUE

      RETURN
      END SUBROUTINE FLUXIN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLUXP(FRM,FRN,FRW,RO,RM,RN,RW,P,U,V,W,UROT,A,AX,AY,AZ,
     &     VOL,RK,EPS2,VIS,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,IL,
     &     KADD,IDI)
C     ****************

C     Flux for axisymmetric flow-case.  by Petri Kaurinkoski 10.3.99
C     Increment in DIM=1 direction is IL. Kinetic energy of turbulence
C     added 11.12.2003 by TSii

      IMPLICIT NONE

      REAL :: FRM(*),FRN(*),FRW(*),P(*),A(*),AX(*),AY(*),AZ(*),VOL(*)
      REAL :: RO(*),RM(*),RN(*),RW(*),U(*),V(*),W(*),EPS2(*),VIS(*)
      REAL :: RK(*)
      REAL :: PAM,PAP,RKP,UM,VM,WM,UP,VP,WP,VFLOM,VFLOP,PERRO,UROT(*)
      REAL :: PER3,YPD2,SURG,SURF,DUC,F2V,F3V,F4V,RMCM,RMCP
      INTEGER :: IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,IL
      INTEGER :: IJK,K1,IJ,J1,I1,II,IP,IM
      INTEGER :: ISTRID,JSTRID,IJSLB,JMAXP1,I,KG,IG,IA,KADD,IDI,IJMAX

C     K1 = K-1
C     J1 = J-1
C     I1 = I-1

      IJMAX = IMAX*JMAX
      PER3   = 1.0/3.0

      DO IJK= 1,IJMAX*KMAX
         K1 = (IJK-1)/IJMAX
         IJ = IJK - K1*IJMAX
         J1 = (IJ-1)/IMAX
         I1 = (IJ-1) - J1*IMAX
         II = (K1+KN)*KSTR + (J1+JN)*JSTR + (I1+IN)*ISTR +1
         IP = II + IL
         IM = II - IL

         UM = 0.5*(RM(IM)+RM(II))
         VM = 0.5*(RN(IM)+RN(II))
         WM = 0.5*(RW(IM)+RW(II))

         UP = 0.5*(RM(IP)+RM(II))
         VP = 0.5*(RN(IP)+RN(II))
         WP = 0.5*(RW(IP)+RW(II))

         PERRO = 1.0/RO(II)

         RMCM  = AX(II)*UM + AY(II)*VM + AZ(II)*WM
         RMCP  = AX(IP)*UP + AY(IP)*VP + AZ(IP)*WP

         VFLOM = A(II)*(RMCM*PERRO-UROT(II))
         VFLOP = A(IP)*(RMCP*PERRO-UROT(IP))

         RKP     = 2.*PER3*RK(II)
         PAM     = (P(II)+RKP)*A(II)
         FRM(II) = PAM*AX(II) + VFLOM*UM
         FRN(II) = PAM*AY(II) + VFLOM*VM
         FRW(II) = PAM*AZ(II) + VFLOM*WM

         PAP     = (P(II)+RKP)*A(IP)
         FRM(IP) = PAP*AX(IP) + VFLOP*UP
         FRN(IP) = PAP*AY(IP) + VFLOP*VP
         FRW(IP) = PAP*AZ(IP) + VFLOP*WP
      END DO

      IF (IDI > 0) THEN        ! Viscous fluxes

         JSTRID = JMAX + 2*JN
         ISTRID = IMAX + 2*IN

         JMAXP1 = JMAX + 1
         IJSLB  = JMAXP1*ISTRID

         DO KG = 1,KMAX+KADD
            IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
            DO IG=1,IJSLB
               I    = IA + IG
               YPD2 = A(I)/(VOL(I) + VOL(I-IL))
               SURG = -A(I)*YPD2
               SURF = SURG*(EPS2(I)*VIS(I) + EPS2(I-IL)*VIS(I-IL))
               DUC  = ( AX(I)*(U(I) - U(I-IL)) +
     &              AY(I)*(V(I) - V(I-IL)) +
     &              AZ(I)*(W(I) - W(I-IL)) )*PER3

               F2V = SURF*(U(I) - U(I-IL) + AX(I)*DUC)
               F3V = SURF*(V(I) - V(I-IL) + AY(I)*DUC)
               F4V = SURF*(W(I) - W(I-IL) + AZ(I)*DUC)

               FRM(I) = FRM(I) + F2V
               FRN(I) = FRN(I) + F3V
               FRW(I) = FRW(I) + F4V
            END DO
         END DO

      END IF

      RETURN
      END SUBROUTINE FLUXP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUDEF(NBL,IBWALL,FR,FRM,FRN,FRW,FE,RO,RM,RN,RW,E,
     +     TEMP,U,V,W,P,PDIFF,FRSPRE,ISTATE,GAMMA,RGAS,
     +     IMAX,JMAX,KMAX,VIS,EPS2,VIST,IN,JN,KN,K1,K2,ITURB,FRK,FEPS,
     +     RK,REPS,FI,FFI,BIJ,S11,XC,YC,ZC,ZZZ,MAXSB,MAXEB,NSCAL,BOUNDF,
     +     IUPDAT,IOLD,IBTYPE,ISTRES,MAXSS,FRSTEM,IUPPTL)

      USE CONSTANTS, ONLY : PR,PRT,VISU0,EXPSU,TSU0,
     +     E0REF,T0REF,PII,EPS

      DIMENSION :: RO(*),RN(*),RW(*),E(*),FR(*),FRN(*),FE(*),
     +     FRM(*),RM(*),VIS(*),EPS2(*),VIST(*),FRW(*),FRK(*),FEPS(*),
     +     RK(*),REPS(*),FI(MAXSB,*),FFI(MAXSB,*),BIJ(MAXEB,*),
     +     S11(MAXSS,*),ZZZ(6+NSCAL,*),
     +     TEMP(*),U(*),V(*),W(*),PDIFF(*),P(*)

      INTEGER :: IN, JN, KN, IUPPTL
      REAL    :: FRSPRE,FRSTEM
      REAL    :: E1(1)
      REAL    :: XC(*), YC(*), ZC(*)

      CHARACTER(LEN=80) :: BOUNDF


C
C   INDEX SYSTEM (LOOP INDECES IN PARANTHESIS):
C  ____________________________          ____________________________
C  |          |    |          |          |          |     |         |
C  |    K2    | K1 |    1     |    OR    |   KMAX   | K1  | K2      |
C  | (KBOUND) | (K)| (KFIELD) |          | (KFIELD) | (K) | (KBOUND)|
C  |          |    |          |          |          |     |         |
C  ----------------------------          ----------------------------
C ... VARIABLES ARE READ FROM A TABLE FOR THE INFLOW BOUNDARY. K-EPSILON
C
C ... PULTATUT ARVOT
C ... TAMA ALIOHJELMA ON RUUVATTU PERAVIRTAUSAMMUKSEN KASITTELYA VARTEN
C ... KOKEILLAAN ALUKSI RUUVATA PERASUIHKUN ARVOT. PPR 9.7
C ... PITAA KAYTTAA DIMENSIOLLISIA MUUTTUJIA
C ... OHJELMA LUKEE KAIKKI MUUTTUJAT
C ... DISTRI ALIOHJELMAN AVULLA TIEDOSTOSTA DISTXY MISSA X ON LOHKO JA
C ... Y ON SEINA. PPR 16.3

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN

      IL      = ISTRID*JSTRID
      IND     = K2 - K1
      IA      = (K1+KN-1)*ISTRID*JSTRID + JN*ISTRID

      IF(ISTATE /= 1) THEN
         WRITE(45,*) 'BCs works only with air'
         WRITE(45,*) 'Contact Patrik'
         WRITE(45,*) 'BCs works with water too. ISTATE =',ISTATE
      ENDIF

      IF (IUPDAT == 0) THEN
      CALL DISTFI(NBL,IBWALL,FR,FRM,FRN,FRW,FE,FRK,FEPS,FFI,S11,XC,YC,
     +  ZC,IMAX,JMAX,ITURB,MAXSB,MAXEB,NSCAL,BOUNDF,IBTYPE,ZZZ,K1,K2,
     +  ISTRES,IFILE,IN,JN,KN)

      DO 600 J = 1,JMAX
      IG = IA + (J-1)*ISTRID + IN
      IG2= (J-1)*IMAX
      DO 600 I = 1,IMAX
         K = I + IG
         N2= I + IG2
         KFIELP  = K  - 2*IND*IL
         KFIELD  = K  - IND*IL
         KBOUND  = K  + IND*IL

         IF(IFILE == 0) THEN
            P(K)      = FR(N2) 
            P(KBOUND) = P(K)
            TEMP(K)   = FE(N2)
            TEMP(KBOUND) = TEMP(K)
         ELSE             
            RO(K)      = FR(N2)
            RO(KBOUND) = RO(K)
            E(K)       = FE(N2)
            E(KBOUND)  = E(K)
         ENDIF

         RM(K)      = FRM(N2)
         RM(KBOUND) = RM(K)
         RN(K)      = FRN(N2)
         RN(KBOUND) = RN(K)
         RW(K)      = FRW(N2)
         RW(KBOUND) = RW(K)

C ... IF NOT K-EPSILON MODEL FLOW IS LAMINAR

         EPS2(K)    = 1.
         EPS2(KBOUND) = EPS2(K)
         VIST(K)    = 0.
         VIST(KBOUND) = VIST(K)

         IF(IOLD <= 0) THEN
            IF(IFILE == 0) THEN
              P(KFIELD)  = P(K)
              RO(KFIELP) = P(K)
              TEMP(KFIELD) = TEMP(K)
              TEMP(KFIELP) = TEMP(K)
            ELSE
              RO(KFIELD) = RO(K)
              RO(KFIELP) = RO(K)
              E(KFIELD)  = E(K)
              E(KFIELP)  = E(K)
            ENDIF

            RM(KFIELD) = RM(K)
            RM(KFIELP) = RM(K)
            RN(KFIELD) = RN(K)
            RN(KFIELP) = RN(K)
            RW(KFIELD) = RW(K)
            RW(KFIELP) = RW(K)
            EPS2(KFIELD) = 1.
            EPS2(KFIELP) = 1.

         ENDIF

 600  CONTINUE

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 700 J = 1,JMAX
      IG = IA + (J-1)*ISTRID + IN
      IG2= (J-1)*IMAX
      DO 700 I = 1,IMAX
         K = I + IG
         N2= I + IG2
         KFIELP  = K  - 2*IND*IL
         KFIELD  = K  - IND*IL
         KBOUND  = K  + IND*IL
         RK(K)        = FRK(N2)
         RK(KBOUND)   = RK(K)
         REPS(K)      = FEPS(N2)
         REPS(KBOUND) = REPS(K)
         EPS2(K)      = 1.+ 0.09*RK(K)**2/(REPS(K)+EPS)/(VIS(K)+EPS)
         EPS2(KBOUND) = EPS2(K)
         VIST(K)      = 0.09*RK(K)**2/(REPS(K)+EPS)
         VIST(KBOUND) = VIST(K)
      
         IF(IOLD <= 0) THEN
         RK(KFIELD)   = RK(K)
         RK(KFIELP)   = RK(K)
         REPS(KFIELD) = REPS(K)
         REPS(KFIELP) = REPS(K)
         EPS2(KFIELD) = EPS2(K)
         EPS2(KFIELP) = EPS2(K)
         VIST(KFIELD) = VIST(K)
         VIST(KFIELP) = VIST(K)
         ENDIF
 700  CONTINUE

      DO 800 NS= 1,NSCAL
      DO 800 J = 1,JMAX
      IG = IA + (J-1)*ISTRID + IN
      IG2= (J-1)*IMAX
      DO 800 I = 1,IMAX
         K = I + IG
         N2= I + IG2
         KFIELP  = K  - 2*IND*IL
         KFIELD  = K  - IND*IL
         KBOUND  = K  + IND*IL
         FI(K,NS)     = FFI(N2,NS)
         FI(KBOUND,NS)= FI(K,NS)
         IF(IOLD <= 0) THEN
         FI(KFIELD,NS)= FI(K,NS)
         FI(KFIELP,NS)= FI(K,NS)
         ENDIF
 800  CONTINUE

      IF(ISTRES > 0) THEN
      DO 900 NS= 1,6
      DO 900 J = 1,JMAX
      IG = IA + (J-1)*ISTRID + IN
      IG2= (J-1)*IMAX
      DO 900 I = 1,IMAX
         K = I + IG
         N2= I + IG2
         KFIELP  = K  - 2*IND*IL
         KFIELD  = K  - IND*IL
         KBOUND  = K  + IND*IL
         BIJ(K,NS)     = S11(N2,NS)
         BIJ(KBOUND,NS)= BIJ(K,NS)
         IF(IOLD <= 0) THEN
         BIJ(KFIELD,NS)= BIJ(K,NS)
         BIJ(KFIELP,NS)= BIJ(K,NS)
         ENDIF
 900  CONTINUE
      ENDIF

      ENDIF ! ITURB >= 3

C ... Specify either conservative or primitive variables depending
C     on the instruction given in DIST-files

      DO J = 1,JMAX

      IG = IA + (J-1)*ISTRID + IN
      IG2= (J-1)*IMAX

      DO I = 1,IMAX

      K = I + IG
      N2= I + IG2
      KFIELP  = K  - 2*IND*IL
      KFIELD  = K  - IND*IL
      KBOUND  = K  + IND*IL

      IF(IFILE == 0) THEN ! Primitive variables were given

         CALL ROFPT(TEMP(K),P(K),RO(K),1,ISTATE,RGAS,GAMMA,FRSDEN,
     +     FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)

         U(K) = RM(K)/RO(K)
         V(K) = RN(K)/RO(K)
         W(K) = RW(K)/RO(K)
         U(KBOUND) = U(K)
         V(KBOUND) = V(K)
         W(KBOUND) = W(K)
     
         CALL EFPT(E1(1),P,TEMP,1,ISTATE,RGAS,GAMMA,FRSPRE,FRSSIE,
     +     E0REF,T0REF)

         IF(ITURB >= 3 .AND. ITURB /= 8) E1(1) = E1(1) + RK(K)/RO(K)
         E(K)     = RO(K)*(E1(1)+.5*(U(K)**2+V(K)**2+W(K)**2))
         E(KBOUND)= E(K)

         PDIFF(K) = P(K) - FRSPRE
         PDIFF(KBOUND) = PDIFF(K)

         IF(IOLD <= 0) THEN
            E(KFIELD)     = E(K)
            PDIFF(KFIELD) = PDIFF(K)
            U(KFIELD)     = U(K)
            V(KFIELD)     = V(K)
            W(KFIELD)     = W(K)

            E(KFIELP)     = E(K)
            PDIFF(KFIELP) = PDIFF(K)
            U(KFIELP)     = U(K)
            V(KFIELP)     = V(K)
            W(KFIELP)     = W(K)
         ENDIF
     
      ELSE ! Conservative variables were given
       
C ... temporary. Better would be that the BCs are given in t, rho*u and p

         U(K) = RM(K)/RO(K)
         V(K) = RN(K)/RO(K)
         W(K) = RW(K)/RO(K)

         IF(ISTATE >= 8) THEN
            WRITE(*,*) 'ISTATE >= 8 requires currently IFILE = 0'
            STOP
         ENDIF
         
         EIN  = E(K)/RO(K) - .5*(U(K)**2+V(K)**2+W(K)**2)
         IF(ITURB >= 3 .AND. ITURB /= 8) EIN = EIN - RK(K)/RO(K)
         E1(1)= EIN
         PDIFF(K) = (GAMMA-1.)*RO(K)*EIN - FRSPRE
         P(K)     = (GAMMA-1.)*RO(K)*EIN
c        TEMP(K)  = (PDIFF(K)+FRSPRE)/(RO(K)*RGAS)
         CALL TFROE(RO(K),E1(1),TEMP(K),1,ISTATE,RGAS,GAMMA,FRSPRE,
     +      E0REF,T0REF)
         IF(ISTATE >= 6) THEN ! viritys maximus
            PDIFF(K) = 0.
c           TEMP(K) = 293.+(EIN - 83.77E+3)/4180. 
            CALL TFROE(RO(K),E1(1),TEMP(K),1,ISTATE,RGAS,GAMMA,FRSPRE,
     +      E0REF,T0REF)
         ENDIF

         TEMP(KBOUND)  = TEMP(K)
         U(KBOUND)     = U(K)
         V(KBOUND)     = V(K)
         W(KBOUND)     = W(K)
         PDIFF(KBOUND) = PDIFF(K)
         P(KBOUND)     = P(K)

         IF(IOLD <= 0) THEN
            TEMP(KFIELD)  = TEMP(K)
            U(KFIELD)     = U(K)
            V(KFIELD)     = V(K)
            W(KFIELD)     = W(K)
            PDIFF(KFIELD) = PDIFF(K)
            P(KFIELD)     = P(K)

            TEMP(KFIELP)  = TEMP(K)
            U(KFIELP)     = U(K)
            V(KFIELP)     = V(K)
            W(KFIELP)     = W(K)
            PDIFF(KFIELP) = PDIFF(K)
            P(KFIELP)     = P(K)
         ENDIF
         ENDIF ! IFILE == 0

      ENDDO
      ENDDO


      ENDIF                     ! IUPDAT ==  0

      RETURN
      END SUBROUTINE BOUDEF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DISTFI(NBL,IBWALL,F1R,F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,
     +     F1FI,S11,XC,YC,ZC,IMAX,JMAX,ITURB,MAXSB,MAXEB,NSCAL,BOUNDF,
     +     IBTYPE,ZZZ,K1,K2,ISTRES,IFILE,IN,JN,KN)

      REAL              :: F1FI(MAXSB,MAX(1,NSCAL)), F1R(*), F1RM(*),
     +                     F1RN(*), F1RW(*), F1E(*), F1RK(*), F1EPS(*),
     +                     ZZZ(6+NSCAL,*), S11(MAXEB,*)
      REAL              :: XC(*), YC(*), ZC(*)
      CHARACTER(LEN=6)  :: NAME
      CHARACTER(LEN=80) :: BOUNDF, LINE1
      LOGICAL           :: THERE
      INTEGER           :: IN, JN, KN
C
C ... READ IN THE BOUNDARY VALUES FOR SCALAR VARIABLES PPR 10.3
C ... IMAX0 IS DIMENSION AT THE TOP LEVEL AND IMAX IS DIMENSION
C ... AT CURRENT LEVEL
C
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      
      IL      = ISTRID*JSTRID
      IND     = K2 - K1
      IA2      = (K1+KN-1)*ISTRID*JSTRID + JN*ISTRID

      INQUIRE(FILE=BOUNDF,EXIST=THERE)
      IF(.NOT. THERE) THEN
         WRITE(*,*) 'Cannot find file:'
         WRITE(*,'(A50)') BOUNDF
         WRITE(*,*) NBL,IBWALL
         WRITE(*,*) ' Exiting...'
         STOP
      ENDIF

      OPEN(51,FILE=BOUNDF,STATUS='OLD',FORM='FORMATTED')

C ... INITIALIZATION FOR AUXLIARY VARIABLES

      DO 100 N = 1,IMAX*JMAX
         F1R(N)  = 0.
         F1RM(N) = 0.
         F1RN(N) = 0.
         F1RW(N) = 0.
         F1E(N)  = 0.
 100  CONTINUE

C ... INITIALIZATION FOR TURBULENCE VARIABLES

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 200 N = 1,IMAX*JMAX
         F1RK(N)  = 0.
         F1EPS(N) = 0.
 200  CONTINUE
      ENDIF

C ... INITIALIZATION FOR SCALARS

      DO 300 NS = 1,NSCAL
      DO 300 N  = 1,IMAX*JMAX
         F1FI(N,NS) = 0.
 300  CONTINUE

C ... INITIALIZATION FOR EARSM ANISOTROPY COMPONENTS

      IF(ISTRES > 0) THEN
      DO 350 NS = 1,6
      DO 350 N  = 1,IMAX*JMAX
         S11(N,NS) = 0.
 350  CONTINUE
      ENDIF

C ... Start to read in distributions

      IF(IBTYPE <= 10) THEN

      IFILE1 = 1 ! Old default value

c     READ(51,*) NSCAL1

      READ(51,'(1A80)') LINE1

	        READ(LINE1,*,END=2013,ERR=2011) NSCAL1,IFILE 
       	 GO TO 2012
2011     WRITE(*,*) 'The following card is broken or erronous:'
	        WRITE(*,*) LINE1
	        STOP
2013     IFILE = IFILE1 ! With old files, the old default value
2012  CONTINUE     ! This line succesfully read

      IF (NSCAL1 < NSCAL) THEN
         WRITE(*,*) 'WRONG NUMBER OF SCALARS AT BOUNDARIES AT FILE',
     +        BOUNDF
         WRITE(*,*) 'EXITING...'
*         CALL EXIT
         STOP
      ENDIF

      READ(51,*) IMAX0,JMAX0
      IF ((MOD(IMAX0,IMAX) /= 0) .OR. (MOD(JMAX0,JMAX) /= 0)) THEN
         WRITE(*,*) 'In inlet file:',imax0,jmax0
         WRITE(*,*) 'In INPUT file:',imax ,jmax
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',BOUNDF
         WRITE(*,*) 'EXITING...'
*         CALL EXIT
         STOP
      ENDIF
      IF (IMAX0/IMAX /= JMAX0/JMAX .AND. IMAX0 /= 1 .AND. JMAX0 /= 1)
     +     THEN
         WRITE(*,*) 'In inlet file:',imax0,jmax0
         WRITE(*,*) 'In INPUT file:',imax ,jmax
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',BOUNDF
         WRITE(*,*) 'EXITING...'
*         CALL EXIT
         STOP
      ENDIF

      NDIVI = INT(IMAX0/IMAX)
      NDIVJ = INT(JMAX0/JMAX)
      IG1 = IMAX0*JMAX0

      ISTRID  = IMAX + 2*IN
      NTOT0   = IMAX0*JMAX0
      NTWO    = NDIVI*NDIVJ
      NONE1I  = NDIVI-1
      NONE1J  = NDIVJ-1


C ... READ PRIMITIVE VALUES IN.

      DO 400 J = 1,JMAX
      IG = (J-1)*IMAX
      DO 400 J2= 0,NONE1J
         DO 400 I = 1,IMAX
         IA = IG + I
         DO 400 I2= 0,NONE1I
            READ(51,*) APU1,APU2,APU3,APU4,APU5
            F1R(IA)  = F1R(IA)  + APU1
            F1RM(IA) = F1RM(IA) + APU2
            F1RN(IA) = F1RN(IA) + APU3
            F1RW(IA) = F1RW(IA) + APU4
            F1E(IA)  = F1E(IA)  + APU5
 400  CONTINUE
      DO 600 N  = 1,IMAX*JMAX
         F1R(N)  = F1R(N)/NTWO
         F1RM(N) = F1RM(N)/NTWO
         F1RN(N) = F1RN(N)/NTWO
         F1RW(N) = F1RW(N)/NTWO
         F1E(N)  = F1E(N)/NTWO
 600  CONTINUE

C ... READ TURBULENCE VARIABLES (K-EPSILON)

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 700 J = 1,JMAX
      IG = (J-1)*IMAX
      DO 700 J2= 0,NONE1J
         DO 700 I = 1,IMAX
         IA = IG + I
         DO 700 I2= 0,NONE1I
            READ(51,*) APU1,APU2
            F1RK(IA)  = F1RK(IA)  + APU1
            F1EPS(IA) = F1EPS(IA) + APU2
 700  CONTINUE
      DO 800 N  = 1,IMAX*JMAX
         F1RK(N)  = F1RK(N)/NTWO
         F1EPS(N) = F1EPS(N)/NTWO
 800  CONTINUE

C ... READ SCALAR VALUES IN.

      DO 900 NS = 1,NSCAL
      DO 900 J = 1,JMAX
      IG = (J-1)*IMAX
      DO 900 J2= 0,NONE1J
         DO 900 I = 1,IMAX
         IA = IG + I
         DO 900 I2= 0,NONE1I
            READ(51,*) APU
            F1FI(IA,NS) = F1FI(IA,NS) + APU
 900  CONTINUE

      DO 1000 NS = 1,NSCAL
      DO 1000 N  = 1,IMAX*JMAX
         F1FI(N,NS) = F1FI(N,NS)/NTWO
1000  CONTINUE

      IF (ISTRES > 0) THEN ! b_ij distributions
      DO 1100 J = 1,JMAX
      IG = (J-1)*IMAX
      DO 1100 J2= 0,NONE1J
         DO 1100 I = 1,IMAX
         IA = IG + I
         DO 1100 I2= 0,NONE1I
            READ(51,*) APU1,APU2,APU3,APU4,APU5
            S11(IA,1) = S11(IA,1) + APU1
            S11(IA,2) = S11(IA,2) + APU2
            S11(IA,3) = S11(IA,3) + APU3
            S11(IA,4) = S11(IA,4) + APU4
            S11(IA,5) = S11(IA,5) + APU5
            S11(IA,6) =-S11(IA,1) - S11(IA,4)
 1100  CONTINUE

      DO 1200 NS = 1,6
      DO 1200 N  = 1,IMAX*JMAX
         S11(N,NS) = S11(N,NS)/NTWO
1200  CONTINUE
      ENDIF ! ISTRES > 0

      ENDIF ! ITURB >= 3

C ... Use free-stream values for a bulk flow

      ELSE                      ! (IBTYPE <= 10)

         READ(51,*)  ! explanation
         READ(51,*)  ! explanation
         READ(51,*)  ! explanation
         READ(51,*)  ! explanation
         READ(51,*)  ! explanation
         READ(51,*) IIMAX,IX,IVEL,NSCAL1
C ... IIMAX = number of points
C ... IX    : 1= function of x, 2=function of y,....
C ... IVEL  : 1= velocity x-direction, ...
C ... NSCAL1: number of scalar in file
         IF (NSCAL1 < NSCAL) THEN
            WRITE(*,*) 'WRONG NUMBER OF SCALARS AT BOUNDARIES AT FILE',
     +           BOUNDF
            WRITE(*,*) 'EXITING...'
*            CALL EXIT
            STOP
         ENDIF
         READ(51,*) 
C ... skaling explanation: XSCALE,XGRIG,ROSCALE,USCALE,e_inskale,REFL
         READ(51,*) XSCALE,XGRIG,ROSCAL,USCAL,EINSKA,REFLEN
         READ(51,*) 
C ... x location, density, velocity, internal, energy,k,epsilon,scalars
         DO I = 1,IIMAX
            IF(ITURB >= 3 .AND. ITURB /= 8)
     +      READ(51,*) (ZZZ(NS,I),NS=1,6+NSCAL)
            IF(ITURB < 3 .OR. ITURB == 8)
     +      READ(51,*) (ZZZ(NS,I),NS=1,4+NSCAL)
         ENDDO
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO J = 1,JMAX
         IG = (J-1)*IMAX
         IGG = IA2 + (J-1)*ISTRID + IN
         DO I = 1,IMAX
            K  = I + IGG
            IA = IG + I
            XAVE = (XC(K)-XGRIG)/XSCALE
            IF(IX == 2) XAVE = (YC(K)-XGRIG)/XSCALE
            IF(IX == 3) XAVE = (ZC(K)-XGRIG)/XSCALE
            DO KK3 = 1,IIMAX-1
               IF(ZZZ(1,KK3) <= XAVE .AND. ZZZ(1,KK3+1) > XAVE) THEN
                   XRES2 = (XAVE - ZZZ(1,KK3))/(ZZZ(1,KK3+1)-ZZZ(1,KK3))
                   XRES1 = 1. - XRES2 
                   F1R(IA)  = XRES1*ZZZ(2,KK3) + XRES2*ZZZ(2,KK3+1)
                   F1RM(IA) = XRES1*ZZZ(3,KK3) + XRES2*ZZZ(3,KK3+1)
                   F1E(IA)  = XRES1*ZZZ(4,KK3) + XRES2*ZZZ(4,KK3+1)
                   F1RK(IA) = XRES1*ZZZ(5,KK3) + XRES2*ZZZ(5,KK3+1)
                   F1EPS(IA)= XRES1*ZZZ(6,KK3) + XRES2*ZZZ(6,KK3+1)
                   DO NS = 1,NSCAL
                      F1FI(IA,NS) = XRES1*ZZZ(6+NS,KK3) + 
     +                     XRES2*ZZZ(6+NS,KK3+1)
                   ENDDO
                   GOTO 313
                ENDIF
             ENDDO
             WRITE(*,*) 'DISTFI1: no hit',I,J,XAVE
             WRITE(*,*) BOUNDF
             STOP
 313         CONTINUE
          ENDDO
       ENDDO
      ELSE                      ! ITURB >= 3
      DO J = 1,JMAX
         IGG = IA2 + (J-1)*ISTRID + IN
         IG = (J-1)*IMAX
         DO I = 1,IMAX
            K  = I + IGG
            IA = IG + I
            XAVE = (XC(K)-XGRIG)/XSCALE
            IF(IX == 2) XAVE = (YC(K)-XGRIG)/XSCALE
            IF(IX == 3) XAVE = (ZC(K)-XGRIG)/XSCALE
            DO KK3 = 1,IIMAX-1
               IF(ZZZ(1,KK3) <= XAVE .AND. ZZZ(1,KK3+1) > XAVE) THEN
                   XRES2 = (XAVE - ZZZ(1,KK3))/(ZZZ(1,KK3+1)-ZZZ(1,KK3))
                   XRES1 = 1. - XRES2 
                   F1R(IA)  = XRES1*ZZZ(2,KK3) + XRES2*ZZZ(2,KK3+1)
                   F1RM(IA) = XRES1*ZZZ(3,KK3) + XRES2*ZZZ(3,KK3+1)
                   F1E(IA)  = XRES1*ZZZ(4,KK3) + XRES2*ZZZ(4,KK3+1)
                   DO NS = 1,NSCAL
                      F1FI(IA,NS) = XRES1*ZZZ(4+NS,KK3) + 
     +                     XRES2*ZZZ(4+NS,KK3+1)
                   ENDDO
                   GOTO 413
                ENDIF
             ENDDO
             WRITE(*,*) 'DISTFI2: no hit',I,J,XAVE
             WRITE(*,*) BOUNDF
             STOP
 413         CONTINUE
          ENDDO
       ENDDO
      ENDIF                     ! ITURB >= 3
C ... skaling explanation: XSCALE,XGRIG,ROSCALE,USCALE,e_inskale,REFL
C       READ(51,*) XSCALE,XGRIG,ROSCAL,USCAL,EINSKA,REFLEN
      DO N  = 1,IMAX*JMAX
         F1R(N)  = F1R(N) * ROSCAL
         F1RM(N) = F1R(N)*F1RM(N)*USCAL
         F1RN(N) = 0.
         F1RW(N) = 0.
c         F1E(N)  = F1R(N)*F1E(N)*EINSKA + .5*(F1RM(N)**2)/F1R(N)
C ... seem to work better. ppr 15.4.99
         F1E(N)  = F1R(N)*F1E(N)*EINSKA !+ .5*(F1RM(N)**2)/F1R(N)
         DO NS = 1,NSCAL
            F1FI(N,NS) = F1R(N)*F1FI(IA,NS)
         ENDDO
      ENDDO
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         DO N  = 1,IMAX*JMAX
            F1RK(N)  = F1R(N)*F1RK(N)*USCAL**2
            F1EPS(N) = F1R(N)*F1EPS(N)*USCAL**3/REFLEN
C            F1E(N)   =  F1E(N) + F1RK(N)/F1R(N)
         ENDDO
      ENDIF
      IF(IVEL == 2) THEN
         DO N  = 1,IMAX*JMAX
            F1RN(N) = F1RM(N)
            F1RM(N) = 0.
            F1RW(N) = 0.
         ENDDO
      ELSEIF(IVEL == 3) THEN
         DO N  = 1,IMAX*JMAX
            F1RW(N) = F1RM(N)
            F1RM(N) = 0.
            F1RN(N) = 0.
         ENDDO
      ENDIF
      IF(ITURB >= 20) THEN
         DO NS = 1,NSCAL
            DO N  = 1,IMAX*JMAX
               F1FI(N,NS) = F1FI(N,NS)*USCAL**2
            ENDDO
         ENDDO
      ENDIF
      ENDIF                     ! IBTYPE >= 10
      CLOSE(51)

      RETURN
      END SUBROUTINE DISTFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

C ... Calculate two perpendicular face vectors (A2 and A3) based on the 
C ... face normal vector A1. Vector A2 is obtained from the vector 
C ... product between normal vector A1 and vector (1,0,0), (0,1,0) or
C ... (0,0,1) depending on which of the absolute values of the normal 
C ... vector components is smallest. A3 is obtained as a vector product  
C ... between A2 and A1. 

      IMPLICIT NONE

      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33, PALPO, AM
      INTEGER :: IDIR

      AM   = ABS(A11)
      IDIR = 1
      IF(ABS(A12) < AM) THEN
         AM = ABS(A12)
         IDIR = 2
      ENDIF
      IF(ABS(A13) < AM) THEN
         IDIR = 3
      ENDIF

      SELECT CASE(IDIR)
         CASE (1)
            A21 =  0.0
            A22 = -A13
            A23 =  A12
         CASE (2)
            A21 =  A13
            A22 =  0.0
            A23 = -A11
         CASE (3)
            A21 = -A12
            A22 =  A11
            A23 =  0.0
      END SELECT
 
      PALPO = 1./SQRT(A21**2 + A22**2 + A23**2)

      A21 = A21*PALPO 
      A22 = A22*PALPO 
      A23 = A23*PALPO 

      A31 = A12*A23 - A13*A22
      A32 = A13*A21 - A11*A23
      A33 = A11*A22 - A12*A21

      RETURN
      END SUBROUTINE SURFVECTS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
