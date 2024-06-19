C ********************************************************************
C     MULTIGRID ROUTINES                                             *
C ********************************************************************

C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANSV(VOL2,VOL,IMAX2,JMAX2,KMAX2,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: VOL2(*),VOL(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VOLUMES AT THE BOUNDARIES FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID

      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN
C ... ISTRID AND JSTRID ARE THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR
C ... TRANSFER OF THE VOLUMES VARIABLES FOR THE BOUNDARIES
      DO 1000 K = 0,KMAX2+1,KMAX2+1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1000 J = 0,JMAX2+1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2+1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJPK   = IJK   + ISTRID + 1
      IJPK    = IJK   + ISTRID
      IPJK    = IJK   + 1
      IJKP    = IJK   + IMJM2             !PPR
      IPJPKP  = IPJPK + IMJM2             !PPR
      IJPKP   = IJPK  + IMJM2             !PPR
      IPJKP   = IPJK  + IMJM2             !PPR

      VOL2(N) = VOL(IPJPK) + VOL(IJK) + VOL(IPJK) + VOL(IJPK)
     2        + VOL(IPJPKP) + VOL(IJKP) + VOL(IPJKP) + VOL(IJPKP)
1000  CONTINUE
      DO 1100 K = 0,KMAX2+1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1100 J = 0,JMAX2+1,JMAX2+1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1100 I = 0,IMAX2+1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJPK   = IJK   + ISTRID + 1
      IJPK    = IJK   + ISTRID
      IPJK    = IJK   + 1
      IJKP    = IJK   + IMJM2             !PPR
      IPJPKP  = IPJPK + IMJM2             !PPR
      IJPKP   = IJPK  + IMJM2             !PPR
      IPJKP   = IPJK  + IMJM2             !PPR

      VOL2(N) = VOL(IPJPK) + VOL(IJK) + VOL(IPJK) + VOL(IJPK)
     2        + VOL(IPJPKP) + VOL(IJKP) + VOL(IPJKP) + VOL(IJPKP)
1100  CONTINUE
      DO 1200 K = 0,KMAX2+1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1200 J = 0,JMAX2+1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1200 I = 0,IMAX2+1,IMAX2+1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJPK   = IJK   + ISTRID + 1
      IJPK    = IJK   + ISTRID
      IPJK    = IJK   + 1
      IJKP    = IJK   + IMJM2             !PPR
      IPJPKP  = IPJPK + IMJM2             !PPR
      IJPKP   = IJPK  + IMJM2             !PPR
      IPJKP   = IPJK  + IMJM2             !PPR

      VOL2(N) = VOL(IPJPK) + VOL(IJK) + VOL(IPJK) + VOL(IJPK)
     2        + VOL(IPJPKP) + VOL(IJKP) + VOL(IPJKP) + VOL(IJPKP)
1200  CONTINUE

C ... EXTRA BOUNDARY LAYER
      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = J*IMAX22 + IN + (JN -1)*IMAX22  + KA
2100  VOL2(II-1)  = VOL2(II)

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = J*IMAX22 + (JN -1)*IMAX22  + KA + IMAX22-1
      II      = J*IMAX22 + (JN -1)*IMAX22  + KA + IMAX22-IN+1
2200  VOL2(II+1)  = VOL2(II)

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
2300  VOL2(N+IMAX22)  = VOL2(N)

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
2400  VOL2(N-IMAX22)  = VOL2(N)

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
2500  VOL2(N-IMCJMC)  = VOL2(N)

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
2600  VOL2(N+IMCJMC)  = VOL2(N)

      RETURN
      END SUBROUTINE TRANSV
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANSF(TEMP2,U2,V2,W2,PD2,TOLD2,UOLD2,VOLD2,
     1     WOLD2,POLD2,ROP2H,RMP2H,RNP2H,RWP2H,EP2H,TEMP,U,V,W,PD,
     2     DRO,DM,DN,DW,DE,VOL2,VOL,EPS22,VIST2,EPS2,VIST,
     3     IMAX2,JMAX2,KMAX2,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: TEMP2(*),U2(*),V2(*),TOLD2(*),PD2(*),UOLD2(*),VOLD2(*),
     2 ROP2H(*),RMP2H(*),RNP2H(*),EP2H(*),TEMP(*),U(*),
     3 V(*),DRO(*),DM(*),DN(*),DE(*),VOL2(*),VOL(*),
     4 W2(*),WOLD2(*),RWP2H(*),W(*),DW(*),PD(*),
     4 EPS22(*),VIST2(*),EPS2(*),VIST(*),POLD2(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR
      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                        !PPR
      IPJPKP  = IPJPK + IMJM2                        !PPR
      IJPKP   = IJPK  + IMJM2                        !PPR
      IPJKP   = IPJK  + IMJM2                        !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES

      TEMP2(N)  = (VOL(IPJPK)*TEMP(IPJPK)  + VOL(IJK)  *TEMP(IJK)  +
     2          VOL(IPJK)  *TEMP(IPJK)   + VOL(IJPK) *TEMP(IJPK) +
     3          VOL(IPJPKP)*TEMP(IPJPKP) + VOL(IJKP) *TEMP(IJKP) +
     4          VOL(IPJKP) *TEMP(IPJKP)  + VOL(IJPKP)*TEMP(IJPKP))*PVOL
      U2(N)  = (VOL(IPJPK)*U(IPJPK)  + VOL(IJK)  *U(IJK)  +
     2          VOL(IPJK)  *U(IPJK)   + VOL(IJPK) *U(IJPK) +
     3          VOL(IPJPKP)*U(IPJPKP) + VOL(IJKP) *U(IJKP) +
     4          VOL(IPJKP) *U(IPJKP)  + VOL(IJPKP)*U(IJPKP))*PVOL
      V2(N)  = (VOL(IPJPK)*V(IPJPK)  + VOL(IJK)  *V(IJK)  +
     2          VOL(IPJK)  *V(IPJK)   + VOL(IJPK) *V(IJPK) +
     3          VOL(IPJPKP)*V(IPJPKP) + VOL(IJKP) *V(IJKP) +
     4          VOL(IPJKP) *V(IPJKP)  + VOL(IJPKP)*V(IJPKP))*PVOL
      W2(N)  = (VOL(IPJPK)*W(IPJPK)  + VOL(IJK)  *W(IJK)  +
     2          VOL(IPJK)  *W(IPJK)   + VOL(IJPK) *W(IJPK) +
     3          VOL(IPJPKP)*W(IPJPKP) + VOL(IJKP) *W(IJKP) +
     4          VOL(IPJKP) *W(IPJKP)  + VOL(IJPKP)*W(IJPKP))*PVOL
      PD2(N)  = (VOL(IPJPK)* PD(IPJPK) + VOL(IJK)  * PD(IJK)  +
     2          VOL(IPJK)  * PD(IPJK)  + VOL(IJPK) * PD(IJPK) +
     3          VOL(IPJPKP)* PD(IPJPKP)+ VOL(IJKP) * PD(IJKP) +
     4          VOL(IPJKP) * PD(IPJKP) + VOL(IJPKP)* PD(IJPKP))*PVOL
      EPS22(N)= (VOL(IPJPK)*EPS2(IPJPK)+ VOL(IJK)  *EPS2(IJK) +
     2          VOL(IPJK)  *EPS2(IPJK) + VOL(IJPK) *EPS2(IJPK)+
     3          VOL(IPJPKP)*EPS2(IPJPKP)+VOL(IJKP) *EPS2(IJKP)+
     4          VOL(IPJKP) *EPS2(IPJKP)+ VOL(IJPKP)*EPS2(IJPKP))*PVOL
      VIST2(N)= (VOL(IPJPK)*VIST(IPJPK)+ VOL(IJK)  *VIST(IJK) +
     2          VOL(IPJK)  *VIST(IPJK) + VOL(IJPK) *VIST(IJPK)+
     3          VOL(IPJPKP)*VIST(IPJPKP)+VOL(IJKP) *VIST(IJKP)+
     4          VOL(IPJKP) *VIST(IPJKP)+ VOL(IJPKP)*VIST(IJPKP))*PVOL

C ... TRANSFERRED RESIDUAL IS STORED AS A FORCING FUNCTION

      ROP2H(N)= (VOL(IPJPK)*DRO(IPJPK)  + VOL(IJK)  *DRO(IJK)  +
     2          VOL(IPJK)  *DRO(IPJK)   + VOL(IJPK) *DRO(IJPK) +
     3          VOL(IPJPKP)*DRO(IPJPKP) + VOL(IJKP) *DRO(IJKP) +
     4          VOL(IPJKP) *DRO(IPJKP)  + VOL(IJPKP)*DRO(IJPKP))*PVOL
      RMP2H(N)= (VOL(IPJPK)*DM (IPJPK)  + VOL(IJK)  *DM (IJK)  +
     2          VOL(IPJK)  *DM (IPJK)   + VOL(IJPK) *DM (IJPK) +
     3          VOL(IPJPKP)*DM (IPJPKP) + VOL(IJKP) *DM (IJKP) +
     4          VOL(IPJKP) *DM (IPJKP)  + VOL(IJPKP)*DM (IJPKP))*PVOL
      RNP2H(N)= (VOL(IPJPK)*DN (IPJPK)  + VOL(IJK)  *DN (IJK)  +
     2          VOL(IPJK)  *DN (IPJK)   + VOL(IJPK) *DN (IJPK) +
     3          VOL(IPJPKP)*DN (IPJPKP) + VOL(IJKP) *DN (IJKP) +
     4          VOL(IPJKP) *DN (IPJKP)  + VOL(IJPKP)*DN (IJPKP))*PVOL
      RWP2H(N)= (VOL(IPJPK)*DW (IPJPK)  + VOL(IJK)  *DW (IJK)  +
     2          VOL(IPJK)  *DW (IPJK)   + VOL(IJPK) *DW (IJPK) +
     3          VOL(IPJPKP)*DW (IPJPKP) + VOL(IJKP) *DW (IJKP) +
     4          VOL(IPJKP) *DW (IPJKP)  + VOL(IJPKP)*DW (IJPKP))*PVOL
      EP2H(N) = (VOL(IPJPK)*DE (IPJPK)  + VOL(IJK)  *DE (IJK)  +
     2          VOL(IPJK)  *DE (IPJK)   + VOL(IJPK) *DE (IJPK) +
     3          VOL(IPJPKP)*DE (IPJPKP) + VOL(IJKP) *DE (IJKP) +
     4          VOL(IPJKP) *DE (IPJKP)  + VOL(IJPKP)*DE (IJPKP))*PVOL
1000  CONTINUE

      DO 2000 K = 1,KMAX2
      KA      = (KN+K-1)*IMCJMC
      DO 2000 J = 1,JMAX2
      II      = (JN+J-1)*IMAX22 + IN + KA
      DO 2000 I = 1,IMAX2
      N       = I + II
      TOLD2(N) = TEMP2(N)
      UOLD2(N) = U2(N)
      VOLD2(N) = V2(N)
      WOLD2(N) = W2(N)
      POLD2(N) = PD2(N)
2000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      ROP2H(II) = 0.
      RMP2H(II) = 0.
      RNP2H(II) = 0.
      RWP2H(II) = 0.
      EP2H(II)  = 0.
      ROP2H(II-1) = 0.
      RMP2H(II-1) = 0.
      RNP2H(II-1) = 0.
      RWP2H(II-1) = 0.
      EP2H(II-1)  = 0.
      TEMP2(II-1) = TEMP2(II)
      U2(II-1)    = U2(II)
      V2(II-1)    = V2(II)
      W2(II-1)    = W2(II)
      PD2(II-1)   = PD2(II)
      EPS22(II-1) = EPS22(II)
      VIST2(II-1) = VIST2(II)
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      ROP2H(II) = 0.
      RMP2H(II) = 0.
      RNP2H(II) = 0.
      RWP2H(II) = 0.
      EP2H(II)  = 0.
      ROP2H(II+1) = 0.
      RMP2H(II+1) = 0.
      RNP2H(II+1) = 0.
      RWP2H(II+1) = 0.
      EP2H(II+1)  = 0.
      TEMP2(II+1) = TEMP2(II)
      U2(II+1)    = U2(II)
      V2(II+1)    = V2(II)
      W2(II+1)    = W2(II)
      PD2(II+1)   = PD2(II)
      EPS22(II+1) = EPS22(II)
      VIST2(II+1) = VIST2(II)
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      RNP2H(N)= 0.
      RWP2H(N)= 0.
      EP2H(N) = 0.
      ROP2H(N+IMAX22) = 0.
      RMP2H(N+IMAX22) = 0.
      RNP2H(N+IMAX22) = 0.
      RWP2H(N+IMAX22) = 0.
      EP2H(N+IMAX22)  = 0.
      TEMP2(N+IMAX22) = TEMP2(N)
      U2(N+IMAX22)    = U2(N)
      V2(N+IMAX22)    = V2(N)
      W2(N+IMAX22)    = W2(N)
      PD2(N+IMAX22)   = PD2(N)
      EPS22(N+IMAX22) = EPS22(N)
      VIST2(N+IMAX22) = VIST2(N)
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      RNP2H(N)= 0.
      RWP2H(N)= 0.
      EP2H(N) = 0.
      ROP2H(N-IMAX22) = 0.
      RMP2H(N-IMAX22) = 0.
      RNP2H(N-IMAX22) = 0.
      RWP2H(N-IMAX22) = 0.
      EP2H(N-IMAX22)  = 0.
      TEMP2(N-IMAX22) = TEMP2(N)
      U2(N-IMAX22)    = U2(N)
      V2(N-IMAX22)    = V2(N)
      W2(N-IMAX22)    = W2(N)
      PD2(N-IMAX22)   = PD2(N)
      EPS22(N-IMAX22) = EPS22(N)
      VIST2(N-IMAX22) = VIST2(N)
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      RNP2H(N)= 0.
      RWP2H(N)= 0.
      EP2H(N) = 0.
      ROP2H(N-IMCJMC) = 0.
      RMP2H(N-IMCJMC) = 0.
      RNP2H(N-IMCJMC) = 0.
      RWP2H(N-IMCJMC) = 0.
      EP2H(N-IMCJMC)  = 0.
      TEMP2(N-IMCJMC) = TEMP2(N)
      U2(N-IMCJMC)    = U2(N)
      V2(N-IMCJMC)    = V2(N)
      W2(N-IMCJMC)    = W2(N)
      PD2(N-IMCJMC)   = PD2(N)
      EPS22(N-IMCJMC) = EPS22(N)
      VIST2(N-IMCJMC) = VIST2(N)
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      RNP2H(N)= 0.
      RWP2H(N)= 0.
      EP2H(N) = 0.
      ROP2H(N+IMCJMC) = 0.
      RMP2H(N+IMCJMC) = 0.
      RNP2H(N+IMCJMC) = 0.
      RWP2H(N+IMCJMC) = 0.
      EP2H(N+IMCJMC)  = 0.
      TEMP2(N+IMCJMC) = TEMP2(N)
      U2(N+IMCJMC)    = U2(N)
      V2(N+IMCJMC)    = V2(N)
      W2(N+IMCJMC)    = W2(N)
      PD2(N+IMCJMC)   = PD2(N)
      EPS22(N+IMCJMC) = EPS22(N)
      VIST2(N+IMCJMC) = VIST2(N)
2600  CONTINUE

      RETURN
      END SUBROUTINE TRANSF
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANKE(RO2,RM2,ROLD2,RMOLD2,ROP2H,RMP2H,RO,RM,DRO,DM,
     2 VOL2,VOL,EPS22,VIST2,EPS2,VIST,IMAX2,JMAX2,KMAX2,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RO2(*),RM2(*),ROLD2(*),RMOLD2(*),ROP2H(*),RMP2H(*),
     3 RO(*),RM(*),DRO(*),DM(*),VOL2(*),VOL(*),
     4 EPS22(*),VIST2(*),EPS2(*),VIST(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2. THIS IS USED FOR K AND EPSILON
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR

      IMCJMC  = IMAX22*JMAX22
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2               !PPR
      IPJPKP  = IPJPK + IMJM2               !PPR
      IJPKP   = IJPK  + IMJM2               !PPR
      IPJKP   = IPJK  + IMJM2               !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES AND TURBULENCE VISCOSITIES

      RO2(N)  = (VOL(IPJPK)*RO(IPJPK) + VOL(IJK)*RO(IJK) +
     2          VOL(IPJK)*RO(IPJK) + VOL(IJPK)*RO(IJPK) +
     3          VOL(IPJPKP)*RO(IPJPKP) + VOL(IJKP)*RO(IJKP) +
     4          VOL(IPJKP)*RO(IPJKP) + VOL(IJPKP)*RO(IJPKP))*PVOL
      RM2(N)  = (VOL(IPJPK)*RM(IPJPK) + VOL(IJK)*RM(IJK) +
     2          VOL(IPJK)*RM(IPJK) + VOL(IJPK)*RM(IJPK) +
     3          VOL(IPJPKP)*RM(IPJPKP) + VOL(IJKP)*RM(IJKP) +
     4          VOL(IPJKP)*RM(IPJKP) + VOL(IJPKP)*RM(IJPKP))*PVOL
      EPS22(N)= (VOL(IPJPK)*EPS2(IPJPK) + VOL(IJK)*EPS2(IJK) +
     2          VOL(IPJK)*EPS2(IPJK) + VOL(IJPK)*EPS2(IJPK) +
     3          VOL(IPJPKP)*EPS2(IPJPKP) + VOL(IJKP)*EPS2(IJKP) +
     4          VOL(IPJKP)*EPS2(IPJKP) + VOL(IJPKP)*EPS2(IJPKP))*PVOL
      VIST2(N)= (VOL(IPJPK)*VIST(IPJPK) + VOL(IJK)*VIST(IJK) +
     2          VOL(IPJK)*VIST(IPJK) + VOL(IJPK)*VIST(IJPK) +
     3          VOL(IPJPKP)*VIST(IPJPKP) + VOL(IJKP)*VIST(IJKP) +
     4          VOL(IPJKP)*VIST(IPJKP) + VOL(IJPKP)*VIST(IJPKP))*PVOL

C ... TRANSFERRED RESIDUAL IS STORED AS A FORCING FUNCTION
      ROP2H(N)= (VOL(IPJPK)*DRO(IPJPK) + VOL(IJK)*DRO(IJK) +
     2          VOL(IPJK)*DRO(IPJK) + VOL(IJPK)*DRO(IJPK) +
     3          VOL(IPJPKP)*DRO(IPJPKP) + VOL(IJKP)*DRO(IJKP) +
     4          VOL(IPJKP)*DRO(IPJKP) + VOL(IJPKP)*DRO(IJPKP))*PVOL
      RMP2H(N)= (VOL(IPJPK)*DM (IPJPK) + VOL(IJK)*DM (IJK) +
     2          VOL(IPJK)*DM (IPJK) + VOL(IJPK)*DM (IJPK) +
     3          VOL(IPJPKP)*DM (IPJPKP) + VOL(IJKP)*DM (IJKP) +
     4          VOL(IPJKP)*DM (IPJKP) + VOL(IJPKP)*DM (IJPKP))*PVOL
1000  CONTINUE
      DO 2000 K = 1,KMAX2
      KA      = (KN+K-1)*IMCJMC
      DO 2000 J = 1,JMAX2
      II      = (JN+J-1)*IMAX22 + IN + KA
      DO 2000 I = 1,IMAX2
      N       = I + II
      ROLD2(N)  = RO2(N)
      RMOLD2(N) = RM2(N)
2000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      ROP2H(II) = 0.
      RMP2H(II) = 0.
      ROP2H(II-1) = 0.
      RMP2H(II-1) = 0.
      RO2(II-1)   = RO2(II)
      RM2(II-1)   = RM2(II)
      EPS22(II-1) = EPS22(II)
      VIST2(II-1) = VIST2(II)
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      ROP2H(II) = 0.
      RMP2H(II) = 0.
      ROP2H(II+1) = 0.
      RMP2H(II+1) = 0.
      RO2(II+1)   = RO2(II)
      RM2(II+1)   = RM2(II)
      EPS22(II+1) = EPS22(II)
      VIST2(II+1) = VIST2(II)
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      ROP2H(N+IMAX22) = 0.
      RMP2H(N+IMAX22) = 0.
      RO2(N+IMAX22)   = RO2(N)
      RM2(N+IMAX22)   = RM2(N)
      EPS22(N+IMAX22) = EPS22(N)
      VIST2(N+IMAX22) = VIST2(N)
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      ROP2H(N-IMAX22) = 0.
      RMP2H(N-IMAX22) = 0.
      RO2(N-IMAX22)   = RO2(N)
      RM2(N-IMAX22)   = RM2(N)
      EPS22(N-IMAX22) = EPS22(N)
      VIST2(N-IMAX22) = VIST2(N)
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      ROP2H(N-IMCJMC) = 0.
      RMP2H(N-IMCJMC) = 0.
      RO2(N-IMCJMC)   = RO2(N)
      RM2(N-IMCJMC)   = RM2(N)
      EPS22(N-IMCJMC) = EPS22(N)
      VIST2(N-IMCJMC) = VIST2(N)
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      ROP2H(N)= 0.
      RMP2H(N)= 0.
      ROP2H(N+IMCJMC) = 0.
      RMP2H(N+IMCJMC) = 0.
      RO2(N+IMCJMC)   = RO2(N)
      RM2(N+IMCJMC)   = RM2(N)
      EPS22(N+IMCJMC) = EPS22(N)
      VIST2(N+IMCJMC) = VIST2(N)
2600  CONTINUE

      RETURN
      END SUBROUTINE TRANKE
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANMF(PRO2,VAR2,PRO,VAR,P2H,VOL2,VOL,IMAX2,JMAX2,
     2 KMAX2,KMAX,NPHASE)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      REAL :: VOL2(*),VOL(*)

      TYPE(PROPERTIES)       :: PRO(*),PRO2(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*),VAR2(*)
      TYPE(MGRID_VARIABLES)  :: P2H(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22
      NTOT    = IMAX22*JMAX22*KMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR

      DO IPHASE = 1,NPHASE

C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR
      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                        !PPR
      IPJPKP  = IPJPK + IMJM2                        !PPR
      IJPKP   = IJPK  + IMJM2                        !PPR
      IPJKP   = IPJK  + IMJM2                        !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES

      PRO2(N)%DTEMP(IPHASE)  = 
     1          (VOL(IPJPK) *PRO(IPJPK)%DTEMP(IPHASE)  + 
     2           VOL(IJK)   *PRO(IJK)%DTEMP(IPHASE)    +
     3           VOL(IPJK)  *PRO(IPJK)%DTEMP(IPHASE)   + 
     4           VOL(IJPK)  *PRO(IJPK)%DTEMP(IPHASE)   +
     5           VOL(IPJPKP)*PRO(IPJPKP)%DTEMP(IPHASE) + 
     6           VOL(IJKP)  *PRO(IJKP)%DTEMP(IPHASE)   +
     7           VOL(IPJKP) *PRO(IPJKP)%DTEMP(IPHASE)  + 
     8           VOL(IJPKP) *PRO(IJPKP)%DTEMP(IPHASE))*PVOL
      VAR2(N)%ALFA(IPHASE)  = 
     1          (VOL(IPJPK) *VAR(IPJPK)%ALFA(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%ALFA(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%ALFA(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%ALFA(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%ALFA(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%ALFA(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%ALFA(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%ALFA(IPHASE))*PVOL
      VAR2(N)%X(IPHASE)     = 
     1          (VOL(IPJPK) *VAR(IPJPK)%X(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%X(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%X(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%X(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%X(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%X(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%X(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%X(IPHASE))*PVOL
      VAR2(N)%EVAPR(IPHASE) = 
     1          (VOL(IPJPK) *VAR(IPJPK)%EVAPR(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%EVAPR(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%EVAPR(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%EVAPR(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%EVAPR(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%EVAPR(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%EVAPR(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%EVAPR(IPHASE))*PVOL
      VAR2(N)%EVAPH(IPHASE) = 
     1          (VOL(IPJPK) *VAR(IPJPK)%EVAPH(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%EVAPH(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%EVAPH(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%EVAPH(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%EVAPH(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%EVAPH(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%EVAPH(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%EVAPH(IPHASE))*PVOL
      VAR2(N)%EQL(IPHASE) = 
     1          (VOL(IPJPK) *VAR(IPJPK)%EQL(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%EQL(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%EQL(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%EQL(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%EQL(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%EQL(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%EQL(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%EQL(IPHASE))*PVOL
      PRO2(N)%QIF(IPHASE) = 
     1          (VOL(IPJPK) *PRO(IPJPK)%QIF(IPHASE)  + 
     2           VOL(IJK)   *PRO(IJK)%QIF(IPHASE)    +
     3           VOL(IPJK)  *PRO(IPJK)%QIF(IPHASE)   + 
     4           VOL(IJPK)  *PRO(IJPK)%QIF(IPHASE)   +
     5           VOL(IPJPKP)*PRO(IPJPKP)%QIF(IPHASE) + 
     6           VOL(IJKP)  *PRO(IJKP)%QIF(IPHASE)   +
     7           VOL(IPJKP) *PRO(IPJKP)%QIF(IPHASE)  + 
     8           VOL(IJPKP) *PRO(IJPKP)%QIF(IPHASE))*PVOL

c       write(866,*) i,j,k,VAR(IPJPK)%ALFA(IPHASE),
c     2           VAR(IJK)%ALFA(IPHASE)    ,
c     3           VAR(IPJK)%ALFA(IPHASE)   , 
c     4           VAR(IJPK)%ALFA(IPHASE)   ,
c     5           VAR(IPJPKP)%ALFA(IPHASE) , 
c     6           VAR(IJKP)%ALFA(IPHASE)   ,
c     7           VAR(IPJKP)%ALFA(IPHASE)  , 
c     8           VAR(IJPKP)%ALFA(IPHASE)

C ... TRANSFERRED RESIDUAL IS STORED AS A FORCING FUNCTION

      P2H(N)%ROP2H(IPHASE)= 
     1          (VOL(IPJPK) *VAR(IPJPK)%DX(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%DX(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%DX(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%DX(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%DX(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%DX(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%DX(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%DX(IPHASE)) * PVOL
      P2H(N)%EP2H(IPHASE)= 
     1          (VOL(IPJPK) *VAR(IPJPK)%DE(IPHASE)  + 
     2           VOL(IJK)   *VAR(IJK)%DE(IPHASE)    +
     3           VOL(IPJK)  *VAR(IPJK)%DE(IPHASE)   + 
     4           VOL(IJPK)  *VAR(IJPK)%DE(IPHASE)   +
     5           VOL(IPJPKP)*VAR(IPJPKP)%DE(IPHASE) + 
     6           VOL(IJKP)  *VAR(IJKP)%DE(IPHASE)   +
     7           VOL(IPJKP) *VAR(IPJKP)%DE(IPHASE)  + 
     8           VOL(IJPKP) *VAR(IJPKP)%DE(IPHASE)) * PVOL
1000  CONTINUE

      DO 2000 K = 1,KMAX2
      KA      = (KN+K-1)*IMCJMC
      DO 2000 J = 1,JMAX2
      II      = (JN+J-1)*IMAX22 + IN + KA
      DO 2000 I = 1,IMAX2
      N       = I + II
      VAR2(N)%DTOLD(IPHASE) = PRO2(N)%DTEMP(IPHASE)
      VAR2(N)%XOLD(IPHASE)  = VAR2(N)%X(IPHASE)
2000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      P2H(II)%ROP2H(IPHASE)   = 0.
      P2H(II)%EP2H(IPHASE)    = 0.
      P2H(II-1)%ROP2H(IPHASE) = 0.
      P2H(II-1)%EP2H(IPHASE)  = 0.
      PRO2(II-1)%DTEMP(IPHASE)= PRO2(II)%DTEMP(IPHASE)
      VAR2(II-1)%ALFA(IPHASE) = VAR2(II)%ALFA(IPHASE)
      VAR2(II-1)%X(IPHASE)    = VAR2(II)%X(IPHASE)
      VAR2(II-1)%EVAPR(IPHASE)= VAR2(II)%EVAPR(IPHASE)
      VAR2(II-1)%EVAPH(IPHASE)= VAR2(II)%EVAPH(IPHASE)
      PRO2(II-1)%QIF(IPHASE)  = PRO2(II)%QIF(IPHASE)
      VAR2(II-1)%EQL(IPHASE)  = VAR2(II)%EQL(IPHASE)
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      P2H(II)%ROP2H(IPHASE)   = 0.
      P2H(II)%EP2H(IPHASE)    = 0.
      P2H(II+1)%ROP2H(IPHASE) = 0.
      P2H(II+1)%EP2H(IPHASE)  = 0.
      PRO2(II+1)%DTEMP(IPHASE) = PRO2(II)%DTEMP(IPHASE)
      VAR2(II+1)%ALFA(IPHASE) = VAR2(II)%ALFA(IPHASE)
      VAR2(II+1)%X(IPHASE)    = VAR2(II)%X(IPHASE)
      VAR2(II+1)%EVAPR(IPHASE)= VAR2(II)%EVAPR(IPHASE)
      VAR2(II+1)%EVAPH(IPHASE)= VAR2(II)%EVAPH(IPHASE)
      PRO2(II+1)%QIF(IPHASE)  = PRO2(II)%QIF(IPHASE)
      VAR2(II+1)%EQL(IPHASE)  = VAR2(II)%EQL(IPHASE)
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      P2H(N)%ROP2H(IPHASE)= 0.
      P2H(N)%EP2H(IPHASE) = 0.
      P2H(N+IMAX22)%ROP2H(IPHASE)= 0.
      P2H(N+IMAX22)%EP2H(IPHASE) = 0.
      PRO2(N+IMAX22)%DTEMP(IPHASE)= PRO2(N)%DTEMP(IPHASE)
      VAR2(N+IMAX22)%ALFA(IPHASE) = VAR2(N)%ALFA(IPHASE)
      VAR2(N+IMAX22)%X(IPHASE)    = VAR2(N)%X(IPHASE)
      VAR2(N+IMAX22)%EVAPR(IPHASE)= VAR2(N)%EVAPR(IPHASE)
      VAR2(N+IMAX22)%EVAPH(IPHASE)= VAR2(N)%EVAPH(IPHASE)
      PRO2(N+IMAX22)%QIF(IPHASE)  = PRO2(N)%QIF(IPHASE)
      VAR2(N+IMAX22)%EQL(IPHASE)  = VAR2(N)%EQL(IPHASE)
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      P2H(N)%ROP2H(IPHASE)= 0.
      P2H(N)%EP2H(IPHASE) = 0.
      P2H(N-IMAX22)%ROP2H(IPHASE)= 0.
      P2H(N-IMAX22)%EP2H(IPHASE) = 0.
      PRO2(N-IMAX22)%DTEMP(IPHASE)= PRO2(N)%DTEMP(IPHASE)
      VAR2(N-IMAX22)%ALFA(IPHASE) = VAR2(N)%ALFA(IPHASE)
      VAR2(N-IMAX22)%X(IPHASE)    = VAR2(N)%X(IPHASE)
      VAR2(N-IMAX22)%EVAPR(IPHASE)= VAR2(N)%EVAPR(IPHASE)
      VAR2(N-IMAX22)%EVAPH(IPHASE)= VAR2(N)%EVAPH(IPHASE)
      PRO2(N-IMAX22)%QIF(IPHASE)  = PRO2(N)%QIF(IPHASE)
      VAR2(N-IMAX22)%EQL(IPHASE)  = VAR2(N)%EQL(IPHASE)
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      P2H(N)%ROP2H(IPHASE)= 0.
      P2H(N)%EP2H(IPHASE) = 0.
      P2H(N-IMCJMC)%ROP2H(IPHASE)= 0.
      P2H(N-IMCJMC)%EP2H(IPHASE) = 0.
      PRO2(N-IMCJMC)%DTEMP(IPHASE)= PRO2(N)%DTEMP(IPHASE)
      VAR2(N-IMCJMC)%ALFA(IPHASE) = VAR2(N)%ALFA(IPHASE)
      VAR2(N-IMCJMC)%X(IPHASE)    = VAR2(N)%X(IPHASE)
      VAR2(N-IMCJMC)%EVAPR(IPHASE)= VAR2(N)%EVAPR(IPHASE)
      VAR2(N-IMCJMC)%EVAPH(IPHASE)= VAR2(N)%EVAPH(IPHASE)
      PRO2(N-IMCJMC)%QIF(IPHASE)  = PRO2(N)%QIF(IPHASE)
      VAR2(N-IMCJMC)%EQL(IPHASE)  = VAR2(N)%EQL(IPHASE)
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      P2H(N)%ROP2H(IPHASE)= 0.
      P2H(N)%EP2H(IPHASE) = 0.
      P2H(N+IMCJMC)%ROP2H(IPHASE)= 0.
      P2H(N+IMCJMC)%EP2H(IPHASE) = 0.
      PRO2(N+IMCJMC)%DTEMP(IPHASE)= PRO2(N)%DTEMP(IPHASE)
      VAR2(N+IMCJMC)%ALFA(IPHASE) = VAR2(N)%ALFA(IPHASE)
      VAR2(N+IMCJMC)%X(IPHASE)    = VAR2(N)%X(IPHASE)
      VAR2(N+IMCJMC)%EVAPR(IPHASE)= VAR2(N)%EVAPR(IPHASE)
      VAR2(N+IMCJMC)%EVAPH(IPHASE)= VAR2(N)%EVAPH(IPHASE)
      PRO2(N+IMCJMC)%QIF(IPHASE)  = PRO2(N)%QIF(IPHASE)
      VAR2(N+IMCJMC)%EQL(IPHASE)  = VAR2(N)%EQL(IPHASE)
2600  CONTINUE

      PRO2(1:NTOT)%TEMP(IPHASE) = PRO2(1:NTOT)%DTEMP(IPHASE)
      VAR2(1:NTOT)%TOLD(IPHASE) = VAR2(1:NTOT)%DTOLD(IPHASE)

      ENDDO
      RETURN
      END SUBROUTINE TRANMF
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANTR(TRM2,TRM,P2HTRM,
     2                  VOL2,VOL,IMAX2,JMAX2,KMAX2,KMAX)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      REAL :: VOL2(*),VOL(*)

      TYPE(INTERMITTENCY)  TRM2(*),TRM(*)
      TYPE(MGRID_INTERMIT) P2HTRM(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR


C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR
      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                        !PPR
      IPJPKP  = IPJPK + IMJM2                        !PPR
      IJPKP   = IJPK  + IMJM2                        !PPR
      IPJKP   = IPJK  + IMJM2                        !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES

      TRM2(N)%G  = 
     1          (VOL(IPJPK) *TRM(IPJPK)%G  + 
     2           VOL(IJK)   *TRM(IJK)%G    +
     3           VOL(IPJK)  *TRM(IPJK)%G   + 
     4           VOL(IJPK)  *TRM(IJPK)%G   +
     5           VOL(IPJPKP)*TRM(IPJPKP)%G + 
     6           VOL(IJKP)  *TRM(IJKP)%G   +
     7           VOL(IPJKP) *TRM(IPJKP)%G  + 
     8           VOL(IJPKP) *TRM(IJPKP)%G)*PVOL
      TRM2(N)%RET  = 
     1          (VOL(IPJPK) *TRM(IPJPK)%RET  + 
     2           VOL(IJK)   *TRM(IJK)%RET    +
     3           VOL(IPJK)  *TRM(IPJK)%RET   + 
     4           VOL(IJPK)  *TRM(IJPK)%RET   +
     5           VOL(IPJPKP)*TRM(IPJPKP)%RET + 
     6           VOL(IJKP)  *TRM(IJKP)%RET   +
     7           VOL(IPJKP) *TRM(IPJKP)%RET  + 
     8           VOL(IJPKP) *TRM(IJPKP)%RET)*PVOL

C ... TRANSFERRED RESIDUAL IS STORED AS A FORCING FUNCTION

      P2HTRM(N)%PG= 
     1          (VOL(IPJPK) *TRM(IPJPK)%DG  + 
     2           VOL(IJK)   *TRM(IJK)%DG    +
     3           VOL(IPJK)  *TRM(IPJK)%DG   + 
     4           VOL(IJPK)  *TRM(IJPK)%DG   +
     5           VOL(IPJPKP)*TRM(IPJPKP)%DG + 
     6           VOL(IJKP)  *TRM(IJKP)%DG   +
     7           VOL(IPJKP) *TRM(IPJKP)%DG  + 
     8           VOL(IJPKP) *TRM(IJPKP)%DG) * PVOL

      P2HTRM(N)%PRET= 
     1          (VOL(IPJPK) *TRM(IPJPK)%DRET  + 
     2           VOL(IJK)   *TRM(IJK)%DRET    +
     3           VOL(IPJK)  *TRM(IPJK)%DRET   + 
     4           VOL(IJPK)  *TRM(IJPK)%DRET   +
     5           VOL(IPJPKP)*TRM(IPJPKP)%DRET + 
     6           VOL(IJKP)  *TRM(IJKP)%DRET   +
     7           VOL(IPJKP) *TRM(IPJKP)%DRET  + 
     8           VOL(IJPKP) *TRM(IJPKP)%DRET) * PVOL
1000  CONTINUE

      DO 2000 K = 1,KMAX2
      KA      = (KN+K-1)*IMCJMC
      DO 2000 J = 1,JMAX2
      II      = (JN+J-1)*IMAX22 + IN + KA
      DO 2000 I = 1,IMAX2
      N       = I + II
      TRM2(N)%GOLD   = TRM2(N)%G
      TRM2(N)%RETOLD = TRM2(N)%RET
2000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      P2HTRM(II)%PG     = 0.
      P2HTRM(II)%PRET   = 0.
      P2HTRM(II-1)%PG   = 0.
      P2HTRM(II-1)%PRET = 0.
      TRM2(II-1)%G      = TRM2(II)%G
      TRM2(II-1)%RET    = TRM2(II)%RET
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      P2HTRM(II)%PG     = 0.
      P2HTRM(II)%PRET   = 0.
      P2HTRM(II+1)%PG   = 0.
      P2HTRM(II+1)%PRET = 0.
      TRM2(II+1)%G      = TRM2(II)%G
      TRM2(II+1)%RET    = TRM2(II)%RET
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      P2HTRM(N)%PG          = 0.
      P2HTRM(N)%PRET        = 0.
      P2HTRM(N+IMAX22)%PG   = 0.
      P2HTRM(N+IMAX22)%PRET = 0.
      TRM2(N+IMAX22)%G      = TRM2(N)%G
      TRM2(N+IMAX22)%RET    = TRM2(N)%RET
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      P2HTRM(N)%PG          = 0.
      P2HTRM(N)%PRET        = 0.
      P2HTRM(N-IMAX22)%PG   = 0.
      P2HTRM(N-IMAX22)%PRET = 0.
      TRM2(N-IMAX22)%G      = TRM2(N)%G
      TRM2(N-IMAX22)%RET    = TRM2(N)%RET
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      P2HTRM(N)%PG          = 0.
      P2HTRM(N)%PRET        = 0.
      P2HTRM(N-IMCJMC)%PG   = 0.
      P2HTRM(N-IMCJMC)%PRET = 0.
      TRM2(N-IMCJMC)%G      = TRM2(N)%G
      TRM2(N-IMCJMC)%RET    = TRM2(N)%RET
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      P2HTRM(N)%PG          = 0.
      P2HTRM(N)%PRET        = 0.
      P2HTRM(N+IMCJMC)%PG   = 0.
      P2HTRM(N+IMCJMC)%PRET = 0.
      TRM2(N+IMCJMC)%G      = TRM2(N)%G
      TRM2(N+IMCJMC)%RET    = TRM2(N)%RET
2600  CONTINUE
      RETURN
      END SUBROUTINE TRANTR
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANSC(FI2,FIOLD2,FIP2H,FI,DFI,VOL2,VOL,
     2  IMAX2,JMAX2,KMAX2,KMAX,MAXSB,MAXS2,NSCAL)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: FI2(MAXSB,MAX(1,NSCAL)),FIOLD2(MAXSB,MAX(1,NSCAL)),
     2 FIP2H(MAXS2,MAX(1,NSCAL)),
     3 FI(MAXSB,MAX(1,NSCAL)),DFI(MAXSB,MAX(1,NSCAL)),VOL2(*),VOL(*)

C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2. THIS IS USED FOR SCALAR EQ PPR 1.3
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR

      DO 3000 NS = 1,NSCAL
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                     !PPR
      IPJPKP  = IPJPK + IMJM2                     !PPR
      IJPKP   = IJPK  + IMJM2                     !PPR
      IPJKP   = IPJK  + IMJM2                     !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES AND TURBULENCE VISCOSITIES

      FI2(N,NS)  = (VOL(IPJPK)*FI(IPJPK,NS) + VOL(IJK)*FI(IJK,NS) +
     2        VOL(IPJK)*FI(IPJK,NS) + VOL(IJPK)*FI(IJPK,NS) +
     3        VOL(IPJPKP)*FI(IPJPKP,NS) + VOL(IJKP)*FI(IJKP,NS) +
     4        VOL(IPJKP)*FI(IPJKP,NS) + VOL(IJPKP)*FI(IJPKP,NS))*PVOL

C ... TRANSFERRED RESIDUAL IS STORED AS A FORCING FUNCTION
      FIP2H(N,NS)= (VOL(IPJPK)*DFI(IPJPK,NS) + VOL(IJK)*DFI(IJK,NS) +
     2      VOL(IPJK)*DFI(IPJK,NS) + VOL(IJPK)*DFI(IJPK,NS) +
     3      VOL(IPJPKP)*DFI(IPJPKP,NS) + VOL(IJKP)*DFI(IJKP,NS) +
     4      VOL(IPJKP)*DFI(IPJKP,NS) + VOL(IJPKP)*DFI(IJPKP,NS))*PVOL
1000  CONTINUE
      DO 2000 K = 1,KMAX2
      KA      = (KN+K-1)*IMCJMC
      DO 2000 J = 1,JMAX2
      II      = (JN+J-1)*IMAX22 + IN + KA
      DO 2000 I = 1,IMAX2
      N       = I + II
      FIOLD2(N,NS) = FI2(N,NS)
2000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      FIP2H(II,NS)   = 0.
      FIP2H(II-1,NS) = 0.
      FI2(II-1,NS)   = FI2(II,NS)
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      FIP2H(II,NS)   = 0.
      FIP2H(II+1,NS) = 0.
      FI2(II+1,NS)   = FI2(II,NS)
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      FIP2H(N,NS)        = 0.
      FIP2H(N+IMAX22,NS) = 0.
      FI2(N+IMAX22,NS)   = FI2(N,NS)
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      FIP2H(N,NS)= 0.
      FIP2H(N-IMAX22,NS) = 0.
      FI2(N-IMAX22,NS)   = FI2(N,NS)
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      FIP2H(N,NS)= 0.
      FIP2H(N-IMCJMC,NS) = 0.
      FI2(N-IMCJMC,NS)   = FI2(N,NS)
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      FIP2H(N,NS)= 0.
      FIP2H(N+IMCJMC,NS) = 0.
      FI2(N+IMCJMC,NS)   = FI2(N,NS)
2600  CONTINUE
3000  CONTINUE

      RETURN
      END SUBROUTINE TRANSC
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANSO(FI2,FI,VOL2,VOL,
     2  IMAX2,JMAX2,KMAX2,KMAX,MAXSB,NSCAL)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: FI2(MAXSB,MAX(1,NSCAL)),FI(MAXSB,MAX(1,NSCAL)),
     &        VOL2(*),VOL(*)
C
C ... THIS SUBROUTINE TRANSFERS THE VARIABLES AND THE RESIDUALS FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2. THIS IS USED FOR SCALAR EQ PPR 1.3
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR

      DO 3000 NS = 1,NSCAL
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                     !PPR
      IPJPKP  = IPJPK + IMJM2                     !PPR
      IJPKP   = IJPK  + IMJM2                     !PPR
      IPJKP   = IPJK  + IMJM2                     !PPR
      PVOL    = 1./VOL2(N)
C ... TRANSFER OF THE DEPENDENT VARIABLES AND TURBULENCE VISCOSITIES

      FI2(N,NS)  = (VOL(IPJPK)*FI(IPJPK,NS) + VOL(IJK)*FI(IJK,NS) +
     2        VOL(IPJK)*  FI(IPJK,NS)  + VOL(IJPK)* FI(IJPK,NS) +
     3        VOL(IPJPKP)*FI(IPJPKP,NS)+ VOL(IJKP)* FI(IJKP,NS) +
     4        VOL(IPJKP)* FI(IPJKP,NS) + VOL(IJPKP)*FI(IJPKP,NS))*PVOL

1000  CONTINUE

C ... BOUNDARY CONDITIONS

      DO 2100 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2100 J = -1,JMAX2+2
      II      = (JN+J-1)*IMAX22 + IN + KA
      FI2(II-1,NS)   = FI2(II,NS)
2100  CONTINUE

      DO 2200 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      DO 2200 J = -1,JMAX2+2
*      II      = (JN+J-1)*IMAX22 + KA + IMAX22-1
      II      = (JN+J-1)*IMAX22 + KA + IMAX22-IN+1
      FI2(II+1,NS)   = FI2(II,NS)
2200  CONTINUE

      DO 2300 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JMAX2+JN-1)*IMAX22 + IMAX22 + IN + KA
      DO 2300 I = -1,IMAX2+2
      N       = II + I
      FI2(N+IMAX22,NS)   = FI2(N,NS)
2300  CONTINUE

      DO 2400 K = -1,KMAX2+2
      KA      = (KN+K-1)*IMCJMC
      II      = (JN-1)*IMAX22  + IN + KA
      DO 2400 I = -1,IMAX2+2
      N       = II + I
      FI2(N-IMAX22,NS)   = FI2(N,NS)
2400  CONTINUE

      DO 2500 J = -1,JMAX2+2
      JA      = IMCJMC
*      II      = (J+JN-1)*IMAX22  + IN + JA
      II      = (J+JN-1)*IMAX22  + IN + (KN-1)*JA
      DO 2500 I = -1,IMAX2+2
      N       = II + I
      FI2(N-IMCJMC,NS)   = FI2(N,NS)
2500  CONTINUE

      DO 2600 J = -1,JMAX2+2
*      JA      = (KMAX22-2)*IMCJMC
      JA      = (KMAX22-KN)*IMCJMC
      II      = (J+JN-1)*IMAX22  + IN + JA
      DO 2600 I = -1,IMAX2+2
      N       = II + I
      FI2(N+IMCJMC,NS)   = FI2(N,NS)
2600  CONTINUE
3000  CONTINUE

      RETURN
      END SUBROUTINE TRANSO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANTI(RO2,RM2,RN2,RW2,E2,RO,RM,RN,RW,E,VOL2,VOL,
     2 IMAX2,JMAX2,KMAX2,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RO2(*),RM2(*),RN2(*),E2(*),RO(*),RM(*),RN(*),E(*),
     2 VOL2(*),VOL(*),RW2(*),RW(*)
C
C ... THIS SUBROUTINE TRANSFERS THE OLD LEVEL OF VARIABLES FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE
C ...  GRID IS IMAX2 X JMAX2 X KMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1

      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR
      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                        !PPR
      IPJPKP  = IPJPK + IMJM2                        !PPR
      IJPKP   = IJPK  + IMJM2                        !PPR
      IPJKP   = IPJK  + IMJM2                        !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES

      RO2(N)  = (VOL(IPJPK)*RO(IPJPK) + VOL(IJK)*RO(IJK) +
     2          VOL(IPJK)*RO(IPJK) + VOL(IJPK)*RO(IJPK) +
     3          VOL(IPJPKP)*RO(IPJPKP) + VOL(IJKP)*RO(IJKP) +
     4          VOL(IPJKP)*RO(IPJKP) + VOL(IJPKP)*RO(IJPKP))*PVOL
      RM2(N)  = (VOL(IPJPK)*RM(IPJPK) + VOL(IJK)*RM(IJK) +
     2          VOL(IPJK)*RM(IPJK) + VOL(IJPK)*RM(IJPK) +
     3          VOL(IPJPKP)*RM(IPJPKP) + VOL(IJKP)*RM(IJKP) +
     4          VOL(IPJKP)*RM(IPJKP) + VOL(IJPKP)*RM(IJPKP))*PVOL
      RN2(N)  = (VOL(IPJPK)*RN(IPJPK) + VOL(IJK)*RN(IJK) +
     2          VOL(IPJK)*RN(IPJK) + VOL(IJPK)*RN(IJPK) +
     3          VOL(IPJPKP)*RN(IPJPKP) + VOL(IJKP)*RN(IJKP) +
     4          VOL(IPJKP)*RN(IPJKP) + VOL(IJPKP)*RN(IJPKP))*PVOL
      RW2(N)  = (VOL(IPJPK)*RW(IPJPK) + VOL(IJK)*RW(IJK) +
     2          VOL(IPJK)*RW(IPJK) + VOL(IJPK)*RW(IJPK) +
     3          VOL(IPJPKP)*RW(IPJPKP) + VOL(IJKP)*RW(IJKP) +
     4          VOL(IPJKP)*RW(IPJKP) + VOL(IJPKP)*RW(IJPKP))*PVOL
      E2(N)   = (VOL(IPJPK)* E(IPJPK) + VOL(IJK)* E(IJK) +
     2          VOL(IPJK)* E(IPJK) + VOL(IJPK)* E(IJPK) +
     3          VOL(IPJPKP)* E(IPJPKP) + VOL(IJKP)* E(IJKP) +
     4          VOL(IPJKP)* E(IPJKP) + VOL(IJPKP)* E(IJPKP))*PVOL
1000  CONTINUE

      RETURN
      END SUBROUTINE TRANTI
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANTK(RO2,RM2,RO,RM,VOL2,VOL,IMAX2,JMAX2,KMAX2,
     2 KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RO2(*),RM2(*),RO(*),RM(*),VOL2(*),VOL(*)
C
C ... THIS SUBROUTINE TRANSFERS THE OLD LEVEL OF VARIABLES FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE
C ...  GRID IS IMAX2 X JMAX2 X KMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1

      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR
      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                        !PPR
      IPJPKP  = IPJPK + IMJM2                        !PPR
      IJPKP   = IJPK  + IMJM2                        !PPR
      IPJKP   = IPJK  + IMJM2                        !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE DEPENDENT VARIABLES ( K AND EPSILON)

      RO2(N)  = (VOL(IPJPK)*RO(IPJPK) + VOL(IJK)*RO(IJK) +
     2          VOL(IPJK)*RO(IPJK) + VOL(IJPK)*RO(IJPK) +
     3          VOL(IPJPKP)*RO(IPJPKP) + VOL(IJKP)*RO(IJKP) +
     4          VOL(IPJKP)*RO(IPJKP) + VOL(IJPKP)*RO(IJPKP))*PVOL
      RM2(N)  = (VOL(IPJPK)*RM(IPJPK) + VOL(IJK)*RM(IJK) +
     2          VOL(IPJK)*RM(IPJK) + VOL(IJPK)*RM(IJPK) +
     3          VOL(IPJPKP)*RM(IPJPKP) + VOL(IJKP)*RM(IJKP) +
     4          VOL(IPJKP)*RM(IPJKP) + VOL(IJPKP)*RM(IJPKP))*PVOL
1000  CONTINUE

      RETURN
      END SUBROUTINE TRANTK
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TRANTS(FI2,FI3,FI2D,FI3D,VOL2,VOL,IMAX2,JMAX2,
     2  KMAX2,KMAX,MAXSB,NSCAL)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: FI2(MAXSB,MAX(1,NSCAL)),FI3(MAXSB,MAX(1,NSCAL)),
     2        FI2D(MAXSB,MAX(1,NSCAL)),
     3        FI3D(MAXSB,MAX(1,NSCAL)),VOL2(*),VOL(*)
C
C ... THIS SUBROUTINE TRANSFERS THE OLD LEVELS OF VARIABLES FROM
C ... THE FINE GRID TO THE COARSE GRID. DIMENSION OF THE COARSE GRID
C ... IS IMAX2 X JMAX2 X KMAX2. THIS IS USED FOR SCALAR EQ PPR 1.3
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      KMAX22  = KMAX2 + 2*KN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      IMJM2   = IMJM                              !PPR
      KIND    = 2                                 !PPR
C ... IMJM2,KIND MAKES 2D CALCULATION POSSIBLE    !PPR
      IF (KMAX == 1) THEN                       !PPR
         IMJM2 = 0                                !PPR
         KIND  = 1                                !PPR
      ENDIF                                       !PPR

      DO 3000 NS = 1,NSCAL
C ... OLI -2
C     DO 1000 K = 0,KMAX22-1
      DO 1000 K = 0,KMAX2 +1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM             !PPR

      DO 1000 J = 0,JMAX2 +1
      IIC     = (JN+J-1)*IMAX22   + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 0,IMAX2 +1
      N       = I     + IIC
      IJK     = 2*I   + IIF - 1
      IPJK    = IJK   + 1
      IJPK    = IJK   + ISTRID
      IPJPK   = IJK   + ISTRID + 1
      IJKP    = IJK   + IMJM2                     !PPR
      IPJPKP  = IPJPK + IMJM2                     !PPR
      IJPKP   = IJPK  + IMJM2                     !PPR
      IPJKP   = IPJK  + IMJM2                     !PPR
      PVOL    = 1./VOL2(N)

C ... TRANSFER OF THE SCALARS

      FI2(N,NS)  = (VOL(IPJPK)*FI2D(IPJPK,NS)+VOL(IJK)*FI2D(IJK,NS) +
     2        VOL(IPJK)*FI2D(IPJK,NS) + VOL(IJPK)*FI2D(IJPK,NS) +
     3        VOL(IPJPKP)*FI2D(IPJPKP,NS) + VOL(IJKP)*FI2D(IJKP,NS) +
     4    VOL(IPJKP)*FI2D(IPJKP,NS) + VOL(IJPKP)*FI2D(IJPKP,NS))*PVOL

      FI3(N,NS)  = (VOL(IPJPK)*FI3D(IPJPK,NS)+VOL(IJK)*FI3D(IJK,NS) +
     2        VOL(IPJK)*FI3D(IPJK,NS) + VOL(IJPK)*FI3D(IJPK,NS) +
     3        VOL(IPJPKP)*FI3D(IPJPKP,NS) + VOL(IJKP)*FI3D(IJKP,NS) +
     4    VOL(IPJKP)*FI3D(IPJKP,NS) + VOL(IJPKP)*FI3D(IJPKP,NS))*PVOL
1000  CONTINUE
3000  CONTINUE

      RETURN
      END SUBROUTINE TRANTS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE FORCFU(ROP2H,RMP2H,RNP2H,RWP2H,EP2H,VAR,P2H,DRO,DM,
     2 DN,DW,DE,RKP2H,EPSP2H,DRK,DEPS,FIP2H,DFI,NTOT,MAXSB,MAXS2,NSCAL,
     2 ITURB,MULPHL,NPHASE,TRANSL,TRM,P2HTRM,MULPHC)
      USE TYPE_ARRAYS

      INTEGER :: NPHASE, IPHASE
      REAL    :: ROP2H(*),RMP2H(*),RNP2H(*),EP2H(*),DRO(*),DM(*),DN(*),
     &           DE(*),RWP2H(*),DW(*),RKP2H(*),EPSP2H(*),DRK(*),DEPS(*),
     &           FIP2H(MAXS2,MAX(1,NSCAL)),DFI(MAXSB,MAX(1,NSCAL))
      LOGICAL :: MULPHL,TRANSL
      CHARACTER(10)          :: MULPHC


      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(MGRID_VARIABLES)  :: P2H(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(MGRID_INTERMIT)   :: P2HTRM(*)

C
C ... CALCULATE THE FORCING FUNCTION FROM THE TRANSFERRED RESIDUAL AND
C ... THE NEW RESIDUAL
C
      IF(ITURB <= 2 .OR. ITURB == 8) THEN
      DO 1000 N = 1,NTOT
         ROP2H(N) = ROP2H(N) - DRO(N)
         RMP2H(N) = RMP2H(N) - DM(N)
         RNP2H(N) = RNP2H(N) - DN(N)
         RWP2H(N) = RWP2H(N) - DW(N)
         EP2H(N)  = EP2H(N)  - DE(N)
 1000 CONTINUE
      ELSEIF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 1100 N = 1,NTOT
         ROP2H(N) = ROP2H(N) - DRO(N)
         RMP2H(N) = RMP2H(N) - DM(N)
         RNP2H(N) = RNP2H(N) - DN(N)
         RWP2H(N) = RWP2H(N) - DW(N)
         EP2H(N)  = EP2H(N)  - DE(N)
         RKP2H(N) = RKP2H(N) - DRK(N)
         EPSP2H(N)= EPSP2H(N)- DEPS(N)
 1100 CONTINUE
      ENDIF

      IF(MULPHL) THEN
      DO IPHASE = 1,NPHASE
      DO 1300 N = 1,NTOT
         P2H(N)%ROP2H(IPHASE) = P2H(N)%ROP2H(IPHASE) - VAR(N)%DX(IPHASE)
         P2H(N)%EP2H(IPHASE)  = P2H(N)%EP2H(IPHASE)  - VAR(N)%DE(IPHASE)
         IF(MULPHC == 'MULTI') THEN
         P2H(N)%RMP2H(IPHASE) = P2H(N)%RMP2H(IPHASE) - VAR(N)%DM(IPHASE)
         P2H(N)%RNP2H(IPHASE) = P2H(N)%RNP2H(IPHASE) - VAR(N)%DN(IPHASE)
         P2H(N)%RWP2H(IPHASE) = P2H(N)%RWP2H(IPHASE) - VAR(N)%DW(IPHASE)
         ENDIF
 1300 CONTINUE
      ENDDO
      ENDIF

      IF(TRANSL) THEN
      DO 1400 N = 1,NTOT
         P2HTRM(N)%PG   = P2HTRM(N)%PG   - TRM(N)%DG
         P2HTRM(N)%PRET = P2HTRM(N)%PRET - TRM(N)%DRET
 1400 CONTINUE
      ENDIF

      IF(NSCAL > 0) THEN
      DO 1200 NS = 1,NSCAL
      DO 1200 N  = 1,NTOT
         FIP2H(N,NS) = FIP2H(N,NS) - DFI(N,NS)
 1200 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FORCFU
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE ADDRES(ROP2H,RMP2H,RNP2H,RWP2H,EP2H,VAR,P2H,DRO,DM,DN,
     2 DW,DE,RKP2H,EPSP2H,DRK,DEPS,FIP2H,DFI,NTOT,MAXSB,MAXS2,NSCAL,
     3 ITURB,MULPHL,NPHASE,TRANSL,TRM,P2HTRM,MULPHC)
      USE TYPE_ARRAYS

      INTEGER :: NPHASE, IPHASE
      REAL    :: ROP2H(*),RMP2H(*),RNP2H(*),EP2H(*),DRO(*),DM(*),DN(*),
     &           DE(*),RWP2H(*),DW(*),RKP2H(*),EPSP2H(*),DRK(*),DEPS(*),
     &           FIP2H(MAXS2,MAX(1,NSCAL)),DFI(MAXSB,MAX(1,NSCAL))
      LOGICAL :: MULPHL,TRANSL

      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(MGRID_VARIABLES)  :: P2H(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(MGRID_INTERMIT)   :: P2HTRM(*)
      CHARACTER(10)          :: MULPHC


C
C ... CALCULATE THE NEW RESIDUAL FROM THE FORCING FUNCTION AND THE
C ... THE OLD RESIDUAL
C
      IF(ITURB <= 2 .OR. ITURB == 8) THEN
      DO 1000 N = 1,NTOT
         DRO(N)  = ROP2H(N) + DRO(N)
         DM(N)   = RMP2H(N) + DM(N)
         DN(N)   = RNP2H(N) + DN(N)
         DW(N)   = RWP2H(N) + DW(N)
         DE(N)   = EP2H(N)  + DE(N)
 1000 CONTINUE
      ELSEIF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 1100 N = 1,NTOT
         DRO(N)  = ROP2H(N) + DRO(N)
         DM(N)   = RMP2H(N) + DM(N)
         DN(N)   = RNP2H(N) + DN(N)
         DW(N)   = RWP2H(N) + DW(N)
         DE(N)   = EP2H(N)  + DE(N)
         DRK(N)  = RKP2H(N) + DRK(N)
         DEPS(N) = EPSP2H(N)+ DEPS(N)
 1100 CONTINUE
      ENDIF

      IF(MULPHL) THEN
      DO IPHASE = 1,NPHASE
      DO 1300 N = 1,NTOT
         VAR(N)%DX(IPHASE)  = P2H(N)%ROP2H(IPHASE) + VAR(N)%DX(IPHASE)
         VAR(N)%DE(IPHASE)  = P2H(N)%EP2H(IPHASE)  + VAR(N)%DE(IPHASE)
         IF(MULPHC == 'MULTI') THEN
         VAR(N)%DM(IPHASE)  = P2H(N)%RMP2H(IPHASE) + VAR(N)%DM(IPHASE)
         VAR(N)%DN(IPHASE)  = P2H(N)%RNP2H(IPHASE) + VAR(N)%DN(IPHASE)
         VAR(N)%DW(IPHASE)  = P2H(N)%RWP2H(IPHASE) + VAR(N)%DW(IPHASE)
         ENDIF
 1300 CONTINUE
      ENDDO
      ENDIF

      IF(TRANSL) THEN
      DO 1400 N = 1,NTOT
         TRM(N)%DG   = P2HTRM(N)%PG   + TRM(N)%DG
         TRM(N)%DRET = P2HTRM(N)%PRET + TRM(N)%DRET
 1400 CONTINUE
      ENDIF

      IF(NSCAL > 0) THEN
      DO 1200 NS = 1,NSCAL
      DO 1200 N = 1,NTOT
         DFI(N,NS)   = FIP2H(N,NS) + DFI(N,NS)
 1200 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE ADDRES
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERP(TEMP2,U2,V2,W2,P2,PD2,TOLD2,UOLD2,VOLD2,
     2 WOLD2,POLD2,TEMP,U,V,W,P,PD,IMAX2,JMAX2,KMAX2,KMAX,
     3 MIB,MIT,MJB,MJT,MKB,MKT)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: TEMP2(*),U2(*),V2(*),W2(*),TOLD2(*),UOLD2(*),
     2 VOLD2(*),WOLD2(*),TEMP(*),U(*),V(*),W(*),
     3 P(*),P2(*),POLD2(*),PD2(*),PD(*)

C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C ... TO THE FINE GRID. DIMENSION OF THE COARSE GRID IS IMAX2 X JMAX2
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      DO 1000 K = 1,KMAX
C      DO 1000 K = MKB,2*KMAX2 + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = 1,2*JMAX2
C      DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = 1,2*IMAX2
C      DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC     = (I+1)/2
      IJK    = I + IIF
      N      = IC + IIC
      DTOK   = TEMP2(N)  - TOLD2(N)
      DUK    =  U2(N)    - UOLD2(N)
      DVK    =  V2(N)    - VOLD2(N)
      DWK    =  W2(N)    - WOLD2(N)
      DPK    = PD2(N)    - POLD2(N)
      
      TEMP(IJK) = TEMP(IJK) + DTOK
      U(IJK)    =    U(IJK) + DUK
      V(IJK)    =    V(IJK) + DVK
      W(IJK)    =    W(IJK) + DWK
      P(IJK)    =    P(IJK) + DPK
      PD(IJK)   =   PD(IJK) + DPK ! APPLIED FOR TWO VARIABLES

1000  CONTINUE

      RETURN
      END SUBROUTINE INTERP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERM(VAR2,PRO2,VAR,PRO,IMAX2,JMAX2,KMAX2,KMAX,
     2 MIB,MIT,MJB,MJT,MKB,MKT,NPHASE)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      INTEGER :: NPHASE, IPHASE

      TYPE(MPHASE_VARIABLES) :: VAR2(*),VAR(*)
      TYPE(PROPERTIES)       :: PRO2(*),PRO(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C     LIMITATION LIKE FOR TURBULENCE QUANTITIES IS OPTIONAL (COMMENTS)
C
      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22
      
      EPS     = 1.E-10
      CC1     = 0.05 ! First guess 11.11.2010 does not work
      CC2     = CC1 + EPS
      CCX1    = 5.E-2
      CCX2    = CCX1 + EPS

      DO IPHASE = 1,NPHASE

      DO 1000 K = 1,KMAX
C      DO 1000 K = MKB,2*KMAX2 + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = 1,2*JMAX2
C      DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = 1,2*IMAX2
C      DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC      = (I+1)/2
      IJK     = I + IIF
      N       = IC + IIC
      DROK    = VAR2(N)%X(IPHASE)   - VAR2(N)%XOLD(IPHASE)
      ALPO    = ABS(DROK)
c      if(i == 1 .and. j ==  13 .and.k == 1)
c     + write(688,*) iphase,VAR(IJK)%X(IPHASE),VAR2(N)%X(IPHASE)

c       write(666,*) iphase,i,j,VAR2(N)%X(IPHASE),VAR2(N)%XOLD(IPHASE),
c     +  VAR(IJK)%X(IPHASE),drok
      VAR(IJK)%X(IPHASE) = VAR(IJK)%X(IPHASE) + 
     + CCX1*VAR(IJK)%X(IPHASE)*DROK/(CCX2*VAR(IJK)%X(IPHASE) + ALPO+EPS)
c      if(i == 1 .and. j ==  13 .and.k == 1)
c     + write(688,*) iphase,VAR(IJK)%X(IPHASE),DROK,ALPO,
c     + CCX1*VAR(IJK)%X(IPHASE)*DROK/(CCX2*VAR(IJK)%X(IPHASE) + ALPO+EPS)
!     + CC1*VAR(IJK)%X(IPHASE)*DROK/(CC2*VAR(IJK)%X(IPHASE) + ALPO)
c       write(688,*) iphase,i,j,VAR(IJK)%X(IPHASE),drok       
      DROK    = PRO2(N)%DTEMP(IPHASE) - VAR2(N)%DTOLD(IPHASE)
      ALPO    = ABS(DROK)
      PRO(IJK)%DTEMP(IPHASE) = PRO(IJK)%DTEMP(IPHASE) +
     + CC1*PRO(IJK)%DTEMP(IPHASE)*DROK/(CC2*PRO(IJK)%DTEMP(IPHASE)+ALPO)
1000  CONTINUE
      ENDDO
      RETURN
      END SUBROUTINE INTERM
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERTR(TRM2,TRM,IMAX2,JMAX2,KMAX2,KMAX,
     2                   MIB,MIT,MJB,MJT,MKB,MKT)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      TYPE(INTERMITTENCY)  TRM2(*),TRM(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C     FOR THE INTERMITTENCY VARIABLES
C
      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN

      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22
      
      EPS     = 1.E-10
      CCG1    = 0.05       !max change of G
      CCG2    = CCG1
      CCR1    = 0.05       !max change of RET
      CCR2    = CCR1

      DO 1000 K = MKB,KMAX + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC      = (I+1)/2
      IJK     = I + IIF
      N       = IC + IIC

      DROK    = TRM2(N)%G   - TRM2(N)%GOLD
c      TRM(IJK)%G   = TRM(IJK)%G + DROK
      ALPO    = ABS(DROK)
      TRM(IJK)%G   = TRM(IJK)%G + 
     &               CCG1*TRM(IJK)%G*DROK/
     &               (CCG2*TRM(IJK)%G + ALPO + EPS)
      
      DROK    = TRM2(N)%RET - TRM2(N)%RETOLD
c      TRM(IJK)%RET = TRM(IJK)%RET + DROK
      ALPO    = ABS(DROK)
      TRM(IJK)%RET = TRM(IJK)%RET + 
     &               CCR1*TRM(IJK)%RET*DROK/
     &               (CCR2*TRM(IJK)%RET + ALPO + EPS)

1000  CONTINUE
      RETURN
      END SUBROUTINE INTERTR
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERV(RO2,ROLD2,RO,ROAUX,IMAX2,JMAX2,KMAX2,KMAX,
     2 ISTRID,JSTRID,IMAX22,JMAX22,IN,JN,KN,MIB,MIT,MJB,MJT,MKB,MKT)

      REAL :: RO2(*),ROLD2(*),RO(*),ROAUX(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C     LIMITATION LIKE FOR TURBULENCE QUANTITIES IS OPTIONAL (COMMENTS)
C
      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22

      DO 1000 K = 1,KMAX
C     DO 1000 K = MKB,2*KMAX2 + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = 1,2*JMAX2
C     DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = 1,2*IMAX2
C     DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC      = (I+1)/2
      IJK     = I + IIF
      N       = IC + IIC
      DROK    = RO2(N)  - ROLD2(N)
      RO(IJK) = RO(IJK) + DROK       ! THE SAME CORRECTION CAN BE 
      ROAUX(IJK) = ROAUX(IJK) + DROK ! APPLIED FOR TWO VARIABLES
1000  CONTINUE

      RETURN
      END SUBROUTINE INTERV
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERO(RO2,ROLD2,RO,IMAX2,JMAX2,KMAX2,KMAX,ISTRID,
     2 JSTRID,IMAX22,JMAX22,IN,JN,KN,CC1,MIB,MIT,MJB,MJT,MKB,MKT)

      REAL :: RO2(*),ROLD2(*),RO(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C
      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22
      EPS     = 1.E-7
      CC2     = CC1 + EPS

      DO 1000 K = MKB,KMAX + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC      = (I+1)/2
      IJK     = I + IIF
      N       = IC + IIC

C ... TRANSFER OF THE CORRECTION

      DROK    = RO2(N) - ROLD2(N)
      ALPO    = ABS(DROK)
      RO(IJK) = RO(IJK) + CC1*RO(IJK)*DROK/(CC2*RO(IJK)+ ALPO)
c     RO(IJK) = RO(IJK) + DROK/(1.+ ALPO/(CC1*RO(IJK)))
c     RO(IJK) = RO(IJK) + SIGN(MIN(CC1*RO(IJK),ALPO),  DROK)
1000  CONTINUE

      RETURN
      END SUBROUTINE INTERO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERT(RO2,RM2,ROLD2,RMOLD2,RO,RM,IMAX2,JMAX2,KMAX2,
     + KMAX,CMGK,CMGEPS,MIB,MIT,MJB,MJT,MKB,MKT)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RO2(*),RM2(*),ROLD2(*),RMOLD2(*),RO(*),RM(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C ... TO THE FINE GRID. DIMENSION OF THE COARSE GRID IS IMAX2 X JMAX2
C ... TWO VARIABLES, I.E. K AND EPSILON ARE TRANSFERRED
C

C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN

      CALL INTERO(RO2,ROLD2,RO,IMAX2,JMAX2,KMAX2,KMAX,ISTRID,JSTRID,
     2 IMAX22,JMAX22,IN,JN,KN,CMGK,MIB,MIT,MJB,MJT,MKB,MKT)
      CALL INTERO(RM2,RMOLD2,RM,IMAX2,JMAX2,KMAX2,KMAX,ISTRID,JSTRID,
     2 IMAX22,JMAX22,IN,JN,KN,CMGEPS,MIB,MIT,MJB,MJT,MKB,MKT)

      RETURN
      END SUBROUTINE INTERT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERS(FI2,FIOLD2,FI,IMAX2,JMAX2,KMAX2,KMAX,MIB,MIT,
     2 MJB,MJT,MKB,MKT,MAXSB,NSCAL,CC1,CC2)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: FI2(MAXSB,MAX(1,NSCAL)),FIOLD2(MAXSB,MAX(1,NSCAL)),
     &        FI(MAXSB,MAX(1,NSCAL))
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C ... TO THE FINE GRID. DIMENSION OF THE COARSE GRID IS IMAX2 X JMAX2
C ... SCALAR VARIABLES ARE TRANSFERRED PPR 1.3.
C
C ... IMAX22 IS THE STRIDE ON ETA-DIRECTION ON THE COARSE GRID
C ... ISTRID IS THE SAME ON THE FINE GRID

      ISTRID  = 2*IMAX2 + 2*IN
      JSTRID  = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN

      DO 1000 NS = 1,NSCAL
c      CALL INTERV(FI2(1,NS),FIOLD2(1,NS),FI(1,NS),IMAX2,JMAX2,KMAX2,
c     2 KMAX,ISTRID,JSTRID,IMAX22,JMAX22,IN,JN,KN,
c     3 MIB,MIT,MJB,MJT,MKB,MKT)
c      CALL INTERO(FI2(1,NS),FIOLD2(1,NS),FI(1,NS),IMAX2,JMAX2,KMAX2,
c     2 KMAX,ISTRID,JSTRID,IMAX22,JMAX22,IN,JN,KN,0.2,
c     3 MIB,MIT,MJB,MJT,MKB,MKT)
      CALL INTERE(FI2(1,NS),FIOLD2(1,NS),FI(1,NS),IMAX2,JMAX2,KMAX2,
     2 KMAX,ISTRID,JSTRID,IMAX22,JMAX22,IN,JN,KN,CC1,CC2,
     3 MIB,MIT,MJB,MJT,MKB,MKT)

1000  CONTINUE

      RETURN
      END SUBROUTINE INTERS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE INTERE(RO2,ROLD2,RO,IMAX2,JMAX2,KMAX2,KMAX,ISTRID,
     2 JSTRID,IMAX22,JMAX22,IN,JN,KN,CC1,CC2,MIB,MIT,MJB,MJT,MKB,MKT)

      REAL :: RO2(*), ROLD2(*), RO(*)
C
C ... THIS SUBROUTINE TRANSFERS THE CORRECTION FROM THE COARSE GRID
C
      IMJM    = ISTRID*JSTRID
      IMCJMC  = IMAX22*JMAX22
      EPS     = 1.E-10

      DO 1000 K = MKB,KMAX + MKT
      KC      = (K+1)/2
      KKC     = (KN+KC-1)*IMCJMC + IN
      KKF     = (KN+K -1)*IMJM   + IN
      DO 1000 J = MJB,2*JMAX2 + MJT
      JC      = (J+1)/2
      IIC     = (JN+JC-1)*IMAX22 + KKC
      IIF     = (JN+J -1)*ISTRID + KKF
      DO 1000 I = MIB,2*IMAX2 + MIT

C ... TRANSFER OF THE CORRECTION

      IC      = (I+1)/2
      IJK     = I + IIF
      N       = IC + IIC

C ... TRANSFER OF THE CORRECTION

      DROK    = RO2(N) - ROLD2(N)
      FAPU    = CC1*ABS(RO(IJK))+CC2*ABS(DROK)
      ALPO    = FAPU*DROK/(FAPU+ABS(DROK) +EPS)
      RO(IJK) = RO(IJK) + ALPO

1000  CONTINUE

      RETURN
      END SUBROUTINE INTERE
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE ININIT(RO2,RO,IMAX2,JMAX2,KMAX2,IMAX1,JMAX1,KMAX1,
     2 IN,JN,KN)

      REAL :: RO2(*),RO(*)
C
C ... THIS SUBROUTINE TRANSFERS THE INITIAL VALUES
C ... FROM THE COARSE GRID TO THE FINEST GRID
C
      ISTRID  = IMAX1 + 2*IN
      JSTRID  = JMAX1 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTRID*JSTRID
      IMJM2   = IMJM
      IMCJMC  = IMAX22*JMAX22
      KMAX3   = KMAX2          !PPR
      IF (KMAX1 == 1) THEN   !PPR
         KMAX3 = 1             !PPR
         IMJM2 = 0             !PPR
      ENDIF                    !PPR
      DO 1000 K = 1,KMAX3      !PPR
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+2*K-2)*IMJM
      DO 1000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 1,IMAX2
      N       = I + IIC
      IJK     = 2*I + IIF - 1
      IPJPK   = IJK + ISTRID + 1
      IJPK    = IJK + ISTRID
      IPJK    = IJK + 1
      IJKP    = IJK + IMJM2
      IPJPKP  = IPJPK + IMJM2
      IJPKP   = IJPK  + IMJM2
      IPJKP   = IPJK  + IMJM2

C ... TRANSFER OF THE VALUES
      DROK       = RO2(N)
      RO(IJK)    = DROK
      RO(IPJPK)  = DROK
      RO(IJPK)   = DROK
      RO(IPJK)   = DROK
      RO(IJKP)   = DROK
      RO(IPJPKP) = DROK
      RO(IJPKP)  = DROK
1000  RO(IPJKP)  = DROK

      RETURN
      END SUBROUTINE ININIT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE ININITm(RO2,RO,IMAX2,JMAX2,KMAX2,IMAX1,JMAX1,KMAX1,
     2 IN,JN,KN)

      REAL :: RO2(*), RO(*)
C
C ... THIS SUBROUTINE TRANSFERS THE INITIAL VALUES
C ... FROM THE COARSE GRID TO THE FINEST GRID
C
c     JOM 12.3.08
c       write(670,*) '---------'
      ISTRID  = IMAX1 + 2*IN
      JSTRID  = JMAX1 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTRID*JSTRID
      IMJM2   = IMJM
      IMCJMC  = IMAX22*JMAX22
      KMAX3   = KMAX2          !PPR
      IF (KMAX1 == 1) THEN   !PPR
         KMAX3 = 1             !PPR
         IMJM2 = 0             !PPR
      ENDIF                    !PPR
      DO 1000 K = 1,KMAX3      !PPR
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+2*K-2)*IMJM
      DO 1000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTRID + IN + KKF
      DO 1000 I = 1,IMAX2
      N       = I + IIC
      IJK     = 2*I + IIF - 1
      IPJPK   = IJK + ISTRID + 1
      IJPK    = IJK + ISTRID
      IPJK    = IJK + 1
      IJKP    = IJK + IMJM2
      IPJPKP  = IPJPK + IMJM2
      IJPKP   = IJPK  + IMJM2
      IPJKP   = IPJK  + IMJM2
C ... TRANSFER OF THE VALUES
      DROK       = RO2(N)
      RO(IJK)    = DROK
      RO(IPJPK)  = DROK
      RO(IJPK)   = DROK
      RO(IPJK)   = DROK
      RO(IJKP)   = DROK
      RO(IPJPKP) = DROK
      RO(IJPKP)  = DROK
c     JOM 12.3.08
c      write(670,*) i,j,k,n,ro2(n),ijk,ro(ijk)
1000  RO(IPJKP)  = DROK

      RETURN
      END SUBROUTINE ININITm
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
C ... Standard manipulation subroutines from res3c.f
      SUBROUTINE SORT50(   YC,ZC,ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: YC(*), ZC(*), APU(*)
C
C ... FROM IJK TO JKI FOR THE CALCULATION OF ROTATIONAL FLUXES
C
CCC   CALL SORTZZ(XC  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(YC  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(ZC  ,ISTRID,JSTRID,KSTRID,NN,APU)
      RETURN
      END SUBROUTINE SORT50
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT51(XC,YC,ZC,ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: XC(*),YC(*),ZC(*)
      REAL, ALLOCATABLE :: APUD(:)
      REAL :: APU(*) 
C
C ... FROM IJK TO JKI FOR THE CALCULATION OF ROTATIONAL FLUXES
C
      ALLOCATE (APUD(ISTRID*JSTRID*KSTRID),STAT=IERRCODE)
      CALL SORTZD(XC  ,ISTRID,JSTRID,KSTRID,NN,APUD)
      CALL SORTZD(YC  ,ISTRID,JSTRID,KSTRID,NN,APUD)
      CALL SORTZD(ZC  ,ISTRID,JSTRID,KSTRID,NN,APUD)
      DEALLOCATE (APUD)

      RETURN
      END SUBROUTINE SORT51
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT38(A1,A1XA,A1YA,A1ZA,D1,ISTRID,JSTRID,KSTRID,
     1 NN,APU)

      REAL :: A1(*),A1XA(*),A1YA(*),A1ZA(*),D1(*),APU(*)
C
C ... FROM XI-ETA TO ETA-XI FOR THE CALCULATION OF FLUXES
C
C     CALL SORTZZ(VOL ,ISTRID,JSTRID,KSTRID,NN,APU) ISOMPI HAH HAH HAA!

      CALL SORTZZ(A1  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1XA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1YA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(D1  ,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT38
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT40(VOL,A1,A1XA,A1YA,A1ZA,D1,ISTRID,JSTRID,KSTRID,
     1 NN,APU)

      REAL :: A1(*),A1XA(*),A1YA(*),A1ZA(*),VOL(*),D1(*),APU(*)
C
C ... FROM XI-ETA TO ETA-XI FOR THE CALCULATION OF FLUXES
C
      CALL SORTZZ(VOL ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1XA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1YA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(D1  ,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT40
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT39(A1,A1XA,A1YA,A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: A1(*),A1XA(*),A1YA(*),A1ZA(*),APU(*)
C
C ... FROM IJK TO JKI FOR THE CALCULATION OF FLUXES
C
      CALL SORTZZ(A1  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1XA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1YA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT39
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT41(VOL,A1,A1XA,A1YA,A1ZA,ISTRID,JSTRID,KSTRID,NN,
     +     APU)

      REAL :: A1(*),A1XA(*),A1YA(*),A1ZA(*),VOL(*),APU(*)
C
C ... FROM IJK TO JKI FOR THE CALCULATION OF FLUXES
C
      CALL SORTZZ(VOL ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1XA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1YA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT41
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT42(P,U,V,W,C,RO,RM,RN,RW,E,VIS,EPS2,VIST,OHMI,
     2 ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: VIS(*),EPS2(*),W(*),P(*),
     2 U(*),V(*),C(*),RO(*),RM(*),RN(*),RW(*),E(*),VIST(*),OHMI(*),
     3 APU(*)
C
C ... FROM XI-ETA TO ETA-XI FOR THE CALCULATION OF FLUXES
C
      CALL SORTZZ(P   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(U   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(V   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(W   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(C   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(RO  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(RM  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(RN  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(RW  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(E   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(VIS ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(EPS2,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(VIST,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(OHMI,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT42
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT2(DTL,DRO,DM,DN,DW,DE,U,V,W,C,RO,VIS,VIST,
     2 ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: DTL(*),DRO(*),DM(*),DN(*),DW(*),DE(*),VIS(*),
     2 U(*),V(*),W(*),C(*),RO(*),VIST(*),APU(*)
C
C ... FROM XI-ETA TO ETA-XI
C
      CALL SORTZZ(DTL ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(DRO ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(DM  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(DN  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(DW  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(DE  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(U   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(V   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(W   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(C   ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(RO  ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(VIS ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(VIST,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT2
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORT43(FWC,FWX,FWY,FWZ,A1XA,A1YA,A1ZA,
     +     ISTRID,JSTRID,KSTRID,NN,APU)

      REAL :: FWX(*),FWY(*),FWZ(*),A1XA(*),A1YA(*),A1ZA(*),
     +             FWC(*),APU(*)
C
C ... FROM IJK TO JKI FOR THE CALCULATION OF MULTI-PHASE WALL FORCE
C
      CALL SORTZZ(FWC ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(FWX ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(FWY ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(FWZ ,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1XA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1YA,ISTRID,JSTRID,KSTRID,NN,APU)
      CALL SORTZZ(A1ZA,ISTRID,JSTRID,KSTRID,NN,APU)

      RETURN
      END SUBROUTINE SORT43
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORTSS(FI,ISTRID,JSTRID,KSTRID,NN,MAX,NSCAL,APU)

      REAL :: FI(MAX,MAX0(1,NSCAL)), APU(*)
C
C ... FROM XI-ETA TO ETA-XI
C
C ... FOR SCALAR FUNCTIONS PPR 11.2
      DO 10 NS = 1,NSCAL
      CALL SORTZZ(FI(1,NS),ISTRID,JSTRID,KSTRID,NN,APU)

 10   CONTINUE

      RETURN
      END SUBROUTINE SORTSS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
C ... sort subroutine for SGI about 20% faster than for cray
      SUBROUTINE SORTZZ(PHIY,IM,JM,KM,NN,APU)

      REAL :: PHIY(*), APU(*)
C
C ... ARRAYS ARE SORTED FROM IJK TO JKI
C
      IF(NN == 1) THEN
           DO 1000 J = 1,JM*KM
           DO 1000 I = 1,IM
           IB        = I - IM
           IA        = (I-1)*JM*KM
1000       APU(J+IA) = PHIY(J*IM + IB)
      ELSE IF(NN == 2) THEN
           IJM       = IM*JM
           DO 1005 K = 1,IJM
           KB        = K - IJM
           KA        = (K-1)*KM
           DO 1005 I = 1,KM
1005       APU(I+KA) = PHIY(I*IJM + KB)
      ENDIF

      PHIY(1:IM*JM*KM) = APU(1:IM*JM*KM)
c      DO 1010 I = 1,IM*JM*KM
c1010  PHIY(I) = APU(I)

      RETURN
      END SUBROUTINE SORTZZ
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SORTZD(PHIY,IM,JM,KM,NN,APU)

      REAL :: PHIY(*), APU(*)
C
C ... ARRAYS ARE SORTED FROM IJK TO JKI
C
      IF(NN == 1) THEN
           DO 1000 J = 1,JM*KM
           DO 1000 I = 1,IM
           IB        = I - IM
           IA        = (I-1)*JM*KM
1000       APU(J+IA) = PHIY(J*IM + IB)
      ELSE IF(NN == 2) THEN
           IJM       = IM*JM
           DO 1005 K = 1,IJM
           KB        = K - IJM
           KA        = (K-1)*KM
           DO 1005 I = 1,KM
1005       APU(I+KA) = PHIY(I*IJM + KB)
      ENDIF

      PHIY(1:IM*JM*KM) = APU(1:IM*JM*KM)
c      DO 1010 I = 1,IM*JM*KM
c1010  PHIY(I) = APU(I)

      RETURN
      END SUBROUTINE SORTZD
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
C ... sort subroutine for cray C90 about 50% faster than above
      SUBROUTINE SORTZC(PHIY,IM,JM,KM,NN,APU)

      REAL :: PHIY(*), APU(*)

C
C ... ARRAYS ARE SORTED FROM IJK TO JKI
C
      IF(NN == 1) THEN
           DO 1000 I = 1,IM
           IB        = I - IM
           IA        = (I-1)*JM*KM
           DO 1000 J = 1,JM*KM
1000       APU(J+IA) = PHIY(J*IM + IB)
      ELSE IF(NN == 2) THEN
           IJM       = IM*JM
           DO 1005 I = 1,KM
           DO 1005 K = 1,IJM
           KB        = K - IJM
           KA        = (K-1)*KM
1005       APU(I+KA) = PHIY(I*IJM + KB)
      ENDIF

      DO 1010 I = 1,IM*JM*KM
1010  PHIY(I) = APU(I)

      RETURN
      END SUBROUTINE SORTZC
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE IND3C(IMAX2,JMAX2,KMAX2,NTOT2,IT2,IL2,IK2,
     2 IMAX,JMAX,KMAX,IT,IL,IK)

      USE NS3CO, ONLY : IN, JN, KN

C
C ... GRID SIZE IDEXES FOR THE NEXT COARSER GRID LEVEL
C

      IF(IMAX2 < 1) IMAX2 = 1
      IF(JMAX2 < 1) JMAX2 = 1
      IF(KMAX2 < 1) KMAX2 = 1
      NTOT2   = (IMAX2+2*IN)*(JMAX2+2*JN)*(KMAX2+2*KN)

      IT2     = (IT+1)/2
      IK2     = (IK+1)/2
      IL2     = IL

      RETURN
      END SUBROUTINE IND3C
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE IND3F2(JBOT2,JTOP2,JBOT,JTOP)
C
C ... SOLID WALL IDEXES FOR THE NEXT COARSER GRID LEVEL
C

      JBOT2 = 0
      JTOP2 = 0
      IF(JBOT /= 0)  JBOT2   = (JBOT  )/2 + 1
      IF(JTOP /= 0)  JTOP2   = (JTOP-1)/2 + 1

      RETURN
      END SUBROUTINE IND3F2
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SETPAT(ICON,NPATCH,NBCS,LEVEL,NBLOCK,
     +     KMAX,MGM,NPPV,IGRID,NPROCE,IPRO)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*),NPATCH(*),KMAX(MGM,*),IGRID(*),NPROCE(*)
C
C ... TRANSFER OF PATCH DATA ON THE APPLIED COARSER LEVEL AND THE
C ... NUMBER OF PATCHES / BLOCK
C
      DO N = 1,NBLOCK
         NPATCH(N) = 0
      ENDDO

      DO L = 1,NBCS
         IVOL = ICON(2,L)
         ICON(23,L) = IVOL
         NPATCH(IVOL) = NPATCH(IVOL) + 1
      ENDDO

      DO N = 1,NBLOCK
         IF(NPATCH(N) > NPPV) THEN
            WRITE(*,*) 'Too many patches in block',N
            WRITE(*,*) 'Change parameter NPPV in NS3CO at least',
     +           NPATCH(N),' (current=',NPPV,')'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
      ENDDO

9500  FORMAT(3X,' WARNING: ',
     +  'Global block',I4,' is freely moving. IMOV was changed to',
     +  ' IMOV=4 in patch',I4)
      WRITE(45,*)


      DO M = 2,LEVEL  ! IF M = 1, NO TRANSFER
         DO L = 1,NBCS
C ... kuulunee tahan ??
            IVOL      = ICON(2,L)

            ICON(4,L) = (ICON(4,L)+1)/2
            ICON(5,L) = (ICON(5,L)+1)/2
            ICON(6,L) = (ICON(6,L)+1)/2
            ICON(7,L) = (ICON(7,L)+1)/2

            ICON(10,L)= (ICON(10,L)+1)/2
            ICON(11,L)= (ICON(11,L)+1)/2
            ICON(12,L)= (ICON(12,L)+1)/2
            ICON(13,L)= (ICON(13,L)+1)/2


            IF(ICON(1,L) == 1) THEN                          !PPR
               IF(KMAX(1,IVOL) == 1 .AND. ICON(3,L) /= 3 .AND. !PPR
     +              ICON(3,L) /= 6) THEN                       !PPR
                  ICON(7,L)  = MAX0(ICON(13,L),1)              !PPR
                  ICON(13,L) = MAX0(ICON(13,L),1)              !PPR
               ENDIF                                           !PPR
            ENDIF

         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SETPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE IPATCH(ICON,IC,NPATCH,MGRID,NBCS,NBLOCK,IAPU,M,N,
     + MGM,KMAX,RCO1,RCO2)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*),NPATCH(*),MGRID(*),IC(MGM,*),IAPU(IC9,*),
     +           KMAX(MGM,*)
      REAL    :: RCO1(IC9,*),RCO2(IC9,*)

C
C ... TRANSFER OF PATCH DATA ON A COARSER LEVELS
C

      NP = 0

      DO 1900 L = 1,NBCS

         IVOL = IAPU(2,L)

      IF(N == IVOL) THEN

         NP = NP + 1

         IF(M > 1) THEN

            ICON(1,NP) = IAPU(1,L)
            ICON(2,NP) = IAPU(2,L)
            ICON(3,NP) = IAPU(3,L)
            
            ICON(4,NP) = (IAPU(4,L)+1)/2
            ICON(5,NP) = (IAPU(5,L)+1)/2
            ICON(6,NP) = (IAPU(6,L)+1)/2
            ICON(7,NP) = (IAPU(7,L)+1)/2

            ICON(8,NP) = IAPU(8,L)
            ICON(9,NP) = IAPU(9,L)
            
            ICON(10,NP)= (IAPU(10,L)+1)/2
            ICON(11,NP)= (IAPU(11,L)+1)/2
            ICON(12,NP)= (IAPU(12,L)+1)/2
            ICON(13,NP)= (IAPU(13,L)+1)/2

            DO II = 14,IC9
               ICON(II,NP)= IAPU(II,L)
            ENDDO

            DO NCO = 1,IC9
               RCO1(NCO,NP) = RCO2(NCO,L)
            ENDDO

         ENDIF

      ENDIF  ! N == IVOL

      IF(NP /= 0) THEN
      IF(ICON(1,NP) == 1) THEN                        !PPR
      IF(KMAX(1,IVOL) == 1 .AND. ICON(3,NP) /= 3 .AND.  !PPR
     +   ICON(3,NP) /= 6) THEN
ccc         write(6,*)   MAX0(ICON(13,NP),2) !PPR
         ICON(7,NP)  = MAX0(ICON(13,NP),1)              !PPR Kokeilee 1
         ICON(13,NP) = MAX0(ICON(13,NP),1)              !PPR
      ENDIF                                             !PPR
      ENDIF                                             !TSii
      ENDIF ! BOUNDARY BLOCK PATCHES


      IF(N < NBLOCK) IC(M,N+1) = NPATCH(N)*IC9 + IC(M,N)
 1900 CONTINUE

      RETURN
      END SUBROUTINE IPATCH
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SETBOT(ICON,NPATCH,IMAX,JMAX,KMAX,IBOT,ITOP,
     + JBOT,JTOP,KBOT,KTOP,IVOL)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*)
C
C ... SET UP THE SOME OF THE TRADITIONAL WALL INDECES FROM ICON ARRAY
C
      IBOT = 0
      ITOP = 0
      JBOT = 0
      JTOP = 0
      KBOT = 0
      KTOP = 0

      NSOLID  = 0
      DO 100  NP = 1,NPATCH
         ITY = ICON(1,NP)
         IF(ITY == 8 .OR. ITY == 9 .OR. ITY == 10 .OR. ITY == 15 .OR.
     +      ITY == 16) NSOLID = NSOLID+1
 100  CONTINUE

      IF(NSOLID  >= 1) WRITE(45,*) NSOLID,
     +     ' SOLID WALL INDICES FOR BLOCK',IVOL
      IF(NSOLID == 0) RETURN

      DO 1000 NP = 1,NPATCH
C ... SOLID WALLS AND TEMPORARILY SINGULARITY AND MORE TEMPORARILY FREE SURFACE
      ITY = ICON(1,NP)
      IF(ITY == 8 .OR. ITY == 9 .OR. ITY == 10 .OR. ITY == 13 .OR. ! Mersu
     +   ITY == 15 .OR. ITY == 16) THEN
         IF(ICON(3,NP) == 1) THEN
         IBOT = 1
         ENDIF

         IF(ICON(3,NP) == 2) THEN
         JBOT = 1
         ENDIF

         IF(ICON(3,NP) == 3) THEN
         KBOT = 1
         ENDIF

         IF(ICON(3,NP) == 4) THEN
         ITOP = IMAX                 !
         ENDIF

         IF(ICON(3,NP) == 5) THEN
         JTOP = JMAX          !+ 1
         ENDIF

         IF(ICON(3,NP) == 6) THEN
         KTOP = KMAX         !+ 1
         ENDIF

      ENDIF

1000  CONTINUE

      RETURN
      END SUBROUTINE SETBOT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SETKCP(ICON,NPATCH,IMAX,JMAX,KMAX,
     + ICP,JCP,KCP,IJMASK,IVOL)

      USE NS3CO, ONLY : IN, JN, KN, IC9

      INTEGER :: ICON(IC9,*),ICP(*),JCP(*),KCP(*),IJMASK(*)
C
C ... SET UP KCP FROM ICON ARRAY 
C
      IF(IVOL == 1) THEN ! only print once
      WRITE(45,*) 'SETUP ICP, JCP AND KCP FROM THE ICON ARRAY'
      WRITE(45,*) 'WALL INDICES (IX1,IX2 ETC.) ARE NOT UPDATED'
      ENDIF

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      NSOLID  = 0

      DO 5100 J = 1,JMAX
      IA      = (   J-1)*ISTRID + IN
      DO 5100 I = 1,IMAX
5100  IJMASK(I+IA)  = 1

      DO 100  NP = 1,NPATCH
         ITY = ICON(1,NP)
         IF(ITY == 8 .OR. ITY == 9 .OR. ITY == 10 .OR. ITY == 15 .OR.
     +      ITY == 16)  NSOLID = NSOLID+1
 100     CONTINUE
      IF(NSOLID  >= 1) WRITE(45,*) NSOLID,
     +        ' SOLID WALL INDICES FOR BLOCK',IVOL
      IF(NSOLID == 0) RETURN

      DO 5000 NP = 1,NPATCH

C ... SOLID WALLS AND TEMPORARILY SINGULARITIES
      ITY = ICON(1,NP)

      IF(ITY == 8 .OR. ITY == 9 .OR. ITY == 10 .OR. ITY == 15 .OR.
     +   ITY == 16)  THEN

         IF(ICON(3,NP) == 1 .OR. ICON(3,NP) == 4) THEN
         IY1  = ICON(4,NP)
         IY2  = ICON(5,NP)
         IZ1  = ICON(6,NP)
         IZ2  = ICON(7,NP)
         WRITE(45,9101) ICON(3,NP),IY1,IY2,IZ1,IZ2
9101     FORMAT(' FACE =',I4,' IY1 =',I4,' IY2 =',I4,' IZ1 =',I4,
     +   ' IZ2 =',I4)
           DO 1000 J = IZ1,IZ2
           DO 1000 I = IY1,IY2
           N       = (J-1+JN)*JSTRID + I  + JN
1000       ICP(N)  = 1
         ENDIF

         IF(ICON(3,NP) == 2  .OR. ICON(3,NP) == 5) THEN
         JX1  = ICON(4,NP)
         JX2  = ICON(5,NP)
         JZ1  = ICON(6,NP)
         JZ2  = ICON(7,NP)
         WRITE(45,9102) ICON(3,NP),JX1,JX2,JZ1,JZ2
9102     FORMAT(' FACE =',I4,' JX1 =',I4,' JX2 =',I4,' JZ1 =',I4,
     +   ' JZ2 =',I4)
           DO 2000 J = JX1,JX2
           DO 2000 I = JZ1,JZ2
           N       = (J-1+JN)*KSTRID + I  + KN
2000       JCP(N)  = 1
         ENDIF

         IF(ICON(3,NP) == 3 .OR. ICON(3,NP) == 6) THEN
         KX1  = ICON(4,NP)
         KX2  = ICON(5,NP)
         KY1  = ICON(6,NP)
         KY2  = ICON(7,NP)
         WRITE(45,9103) ICON(3,NP),KX1,KX2,KY1,KY2
9103     FORMAT(' FACE =',I4,' KX1 =',I4,' KX2 =',I4,' KY1 =',I4,
     +   ' KY2 =',I4)
           DO 3000 J = KY1,KY2
           DO 3000 I = KX1,KX2
           N       = (J-1+JN)*ISTRID + I  + IN
3000       KCP(N)  = 1
         ENDIF

      ENDIF       !      FOR WALLS AND SINGULARITIES
5000  CONTINUE

      WRITE(45,*)
      WRITE(45,*)
      RETURN
      END SUBROUTINE SETKCP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE IND3W3(IW,ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,M,N,
     & MGRID,NBLOCK,NB,NPPV,ITURB,ISTRES,NSCAL,ITIMES,IWMAX,IPMAX,
     & MULPHL,NPHASES,TRANSL,TWO_FLUIDL)
C
C ... POINTERS FOR PATCHES IN WORK ARRAY
C
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,M,N,MGRID,NBLOCK,NB,
     &        NPPV,ITURB,ISTRES,NSCAL,ITIMES,IWMAX,IPMAX,NPHASES
      INTEGER :: IWOLD, IOLOC, IOGLO, L, ITYPE
      INTEGER :: IA, IM, JA, JM, IA2, IM2, JA2, JM2
      INTEGER :: IW(NPPV,NB,*), ICON(IC9,*)
      LOGICAL :: MULPHL, TRANSL, TWO_FLUIDL

C ... CHECK OUT LENGHT OF THE IO PATHC

      CALL IWLENG(ITIMES,ITURB,NSCAL,ISTRES,MULPHL,NPHASES,TRANSL,
     &   TWO_FLUIDL)

      IWOLD = IW(1,N,M)
      IOLOC = 2*IN*2*JN
      IOGLO = 2*IN*2*JN

      IF(ICON(1,1) ==  1 .OR. ICON(1,1) ==  6 .OR. ICON(1,1) == 11 .OR.
     &   ICON(1,1) == 14 .OR. ICON(1,1) == 15) THEN

         L = 1

         CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2) ! Extend patch
                  
         IOLOC = (IM-IA+1)*(JM-JA+1)
         IOGLO = (IM2-IA2+1)*(JM2-JA2+1)

         IF(ICON(18,1) /= 1) IOLOC = MAX0(IOLOC,IOGLO)
      ENDIF

      IF(M == 1) THEN

      DO 16 L = 2,NPATCH

         ITYPE = ICON(1,L)
         
         IF(ITYPE ==  1 .OR. ITYPE ==  6 .OR. ITYPE == 11 .OR.
     &      ITYPE == 14 .OR. ITYPE == 15) THEN

            CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2) ! Extend patch

            IW(L,N,M) = IWOLD + ITIMES*2*IOLOC
            IPMAX     = MAX(IPMAX,ITIMES*2*IOLOC)
            IOLOC = (IM-IA+1)*(JM-JA+1)
            IOGLO = (IM2-IA2+1)*(JM2-JA2+1)
            IF(ICON(18,L) /= 1) IOLOC = MAX0(IOLOC,IOGLO)
            IWOLD     = IW(L,N,M)
            IWMAX     = IW(L,N,M) + ITIMES*2*IOLOC
         ENDIF
 16   CONTINUE

      IW(1,N+1,1) = IWOLD + ITIMES*2*IOLOC
      IWMAX       = IWOLD + ITIMES*2*IOLOC
      IPMAX       = MAX(IPMAX,ITIMES*2*IOLOC)

      ELSE

      DO 17 L = 1,NPATCH

         ITYPE = ICON(1,L)

         IF(ITYPE ==  1 .OR. ITYPE ==  6 .OR. ITYPE == 11 .OR.
     &      ITYPE == 14 .OR. ITYPE == 15) THEN
            IW(L,N,M) = IW(L,N,M-1) 
         ENDIF

 17   CONTINUE

      IW(1,N,M+1) = IW(1,N,M)

      ENDIF

      RETURN
      END SUBROUTINE IND3W3
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE PRIICO(ICON,NBCS)

      USE NS3CO, ONLY : IC9

      INTEGER :: ICON(IC9,*)

      IBEDEL = 0
      WRITE(45,1102)
      WRITE(45,1105)
      WRITE(45,1104)
 1102 FORMAT(//4X,'BOUNDARY CONDITION PATCHES')
 1105 FORMAT(//4X,' 1  = Connectivity     2  = External     ',
     &            ' 3  = Inlet',
     &        /4X,' 4  = Mirror           5  = Outlet       ',
     &            ' 6  = Cyclic',
     &        /4X,' 7  = Singularity      8  = Solid        ',
     &            ' 9  = Rotating solid',
     &        /4X,' 10 = Moving solid    11  = Sliding      ',
     &            ' 12 = Chimera',
     &        /4X,' 13 = Free surface    14  = Mixing       ',
     &            ' 15 = Coupling with solid',
     &        /4X,' 16 = Slip            17  = Circulating inlet')
 1103 FORMAT(//4X,'LOCAL BLOCK NO. ',I4,3X,'GLOBAL BLOCK NO. ',I4,//,
     &         '   PTC    1    2    3    4    5    6    7    8    9',
     &         '   10   11   12   13   14',/,
     &         '  =================================================',
     &         '=========================')
 1107 FORMAT(//'   PTC   15   16   17   18   19   20   21   22   23',
     &         '   24   25   26   27   28',/,
     &         '  =================================================',
     &         '=========================')
 1104 FORMAT(//2X,
     &' 1 = Boundary condition type              '/,2X,
     &' 2 = Process local patch (or block)       '/,2X,
     &' 3 = Face number                          '/,2X,
     &'4-7= x,y low,up own                       '/,2X,
     &' 8 = number of connective block           '/,2X,
     &' 9 = face of the c. patch or MOV-type (0,1,2,3,4)'/,2X,
     &'10-13=x,y low,up in connective patch      '/,2X,
     &'14 = Orientation                          '/,2X,
     &'15 = connective patch n. (block l.)       '/,2X,
     &'16 = connective patch n. (proces l.)      '/,2X,
     &     '17 = number of connectivity process  '/,2X,
     &     '18 = connectivity method (1=normal,2)'/,2X,
     &     '19 = Multigrid levels in connective p'/,2X,
     &'20 = H. flux t.(0=adi,1=temp,2=h.f. and 1x read from files)',
     &      'or CON/PER switch'/,2X,
     &     '21 = Order of connection             '/,2X,
     &     '22 = nonmatching (0=norm. 1=nomatch) '/,2X,
     &     '23 = Local block number              '/,2X,
     &     '24 = Global block number             '/,2X,
     &     '25 = Global patch number             '/,2X,
     &     '26 = Actuator disc number            '/,2X,
     &     '27 = Inlet/outlet type               '/,2X,
     &     '28 = Connective global block number  ')
       
      ILOBC = 0
      DO 17 I=1,NBCS
        IF(ICON(24,I) > IBEDEL) THEN
           ILOBC = ILOBC + 1
          WRITE(45,1103) ILOBC,ICON(24,I)
          IBEDEL = ICON(24,I)
          NOLD   = I
        ENDIF
        WRITE(45,111) I,(ICON(J,I),J=1,14)
        IF(I /= NBCS) THEN
        IF(ICON(24,I+1) > IBEDEL) THEN
           WRITE(45,1107)
           DO  133 II = NOLD,I
              WRITE(45,111) II,(ICON(J,II),J=15,28)
 133       CONTINUE
        ENDIF
        ENDIF
 111    FORMAT(1X,15I5)
 17   CONTINUE
      WRITE(45,1107)
      DO  134 II = NOLD,NBCS
         WRITE(45,111) II,(ICON(J,II),J=15,28)
 134  CONTINUE
      WRITE(45,111)
      WRITE(45,111)
      WRITE(45,111)

      RETURN
      END SUBROUTINE PRIICO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE PRIPAT(ICON,IC,NPATCH,MGRID,NBCS,NBLOCK,M,MGM,IC11)

      INTEGER :: ICON(IC11,*), NPATCH(*), MGRID(*), IC(MGM,*)
C
C ... TRANSFER OF PATCH DATA ON A COARSER LEVELS
C
      MGLEV     = 0

C ... CHECK WETHER THE CURRENT GRID LEVEL (M) IS ACTIVE IN ANY BLOCK

      DO 1900 N = 1,NBLOCK
      IF(M <= MGRID(N)) MGLEV = MGLEV + 1
1900  CONTINUE

C ... REPLACE THE BLOCK NUMBER BY A GLOBAL PATCH NUMBER

      NP         = 0
      DO 2100 N  = 1,NBLOCK
      DO 2100 L  = 1,NPATCH(N)
      NP         = NP + 1
      ICON(2,NP) = NP
2100  CONTINUE

C
C ... Write out the boundary condition patches for checking purposes
C

       IF(MGLEV > 0) THEN
          IF(M == 1) WRITE(45,1101)
          IF(M > 1) WRITE(45,1100) M
          WRITE(45,1102)
 1100 FORMAT(//4X,'BOUNDARY CONDITION PATCHES ON LEVEL ',I3)
 1101 FORMAT(//4X,'BOUNDARY CONDITION PATCHES ON THE FIRST ',
     &            'GRID LEVEL:')
 1102 FORMAT(//4X,' 1  = Connectivity     2  = External     ',
     &            ' 3  = Inlet',
     &        /4X,' 4  = Mirror           5  = Outlet       ',
     &            ' 6  = Cyclic',
     &        /4X,' 7  = Singularity      8  = Solid        ',
     &            ' 9  = Rotating solid',
     &        /4X,' 10 = Moving solid    11  = Sliding      ',
     &            ' 12 = Chimera',
     &        /4X,' 13 = Free surface    14  = Mixing       ',
     &            ' 15 = Coupling with solid',
     &        /4X,' 16 = Slip            17  = Circulating inlet')
c 1103 FORMAT(//4X,'BLOCK NO. ',I4,' LEVEL = ',I2,'  NUMBER OF PATCHES =
c     &',I4,'        CONNECTIVITY'//,
c     &         4X,'BC   BLK FACE  XLO  XUP  YLO  YUP  BLK FACE  ',
c     &'XLO  XUP  YLO  YUP EMP LPN GPN PRO TY MGS PER'/)
 1103 FORMAT(//4X,'BLOCK NO.',I6,' LEVEL = ',I2,'  NUMBER OF PATCHES =
     &',I4,'        CONNECTIVITY'//,
     &         4X,'BC   BLK FACE XLO  XUP  YLO  YUP    BLK  FACE ',
     &'XLO  XUP  YLO  YUP EMP LPN   GPN PRO  TY MGS PER'/)
c 1104 FORMAT(/4x,'STARTING ADRESS IC(',I1,',',I3,')',I6)
 1104 FORMAT(/4x,'STARTING ADDRESS IC(',I1,',',I6,')',I8)
      IALKU   = 1
      DO 17 N = 1,NBLOCK
      IF(M <= MGRID(N)) THEN
        WRITE(45,1103) N,M,NPATCH(N)
      DO 16 I = IALKU,IALKU+NPATCH(N)-1
        WRITE(45,111) (ICON(J,I),J=1,20)
c 111    FORMAT(4X,I3,I4,1X,I4,4I5,2X,I3,2X,I3,4I5,I4,I4,I4,5I4)
 111    FORMAT(3X,I3,I6,I4,4I5,1X,I6,I5,4I5,I4,I4,I6,5I4)
 16   CONTINUE
        IALKU = IALKU + NPATCH(N)
        WRITE(45,1104) M,N,IC(M,N)
      ENDIF
 17   CONTINUE

      ENDIF
      WRITE(45,*)
      WRITE(45,*)

      RETURN
      END SUBROUTINE PRIPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE YESNO(STARTC,STARTL)

      CHARACTER(LEN=80) :: STARTC
      LOGICAL :: STARTL

      STARTL = STARTC == 'YES' .OR. STARTC == 'yes' .OR. STARTC == 'Yes'
     &.OR. STARTC == 'START'.OR. STARTC == 'start'.OR. STARTC == 'Start'

      RETURN
      END SUBROUTINE YESNO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE YESNO_STARTL(STARTC,STARTL)

      CHARACTER(LEN=80) :: STARTC
      LOGICAL :: STARTL

      STARTL = 
     &   STARTC == 'START' .OR. STARTC == 'start' .OR. STARTC == 'Start'

      RETURN
      END SUBROUTINE YESNO_STARTL
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE IJKPAI(LL,IMAX,JMAX,KMAX,M,N,L)

C ... Solve I, J, K from the cell index LL

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, ISTRID, JSTRID, IL, LL, LAPU, L, N, M
      
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      LAPU   = LL

      L      = (LAPU-1)/IL - KN + 1
      LAPU   = LAPU - (L+KN-1)*IL

      N      = (LAPU-1)/ISTRID -JN + 1
      LAPU   = LAPU - (N+JN-1)*ISTRID

      M      = LAPU - IN

      RETURN
      END SUBROUTINE IJKPAI
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUMCH6(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFER NUMBER ICYCLE TO THE 6 DIGIT CHARACTER FILEE

      CHARACTER(LEN=6) :: FILEE

      IAPU  = ICYCLE
      I7    = INT(IAPU/1000000)
      IAPU  = IAPU - I7*100000
      I6    = INT(IAPU/100000)
      IAPU  = IAPU - I6*100000
      I5    = INT(IAPU/10000)
      IAPU  = IAPU - I5*10000
      I4    = INT(IAPU/1000)
      IAPU  = IAPU - I4*1000
      I3    = INT(IAPU/100)
      IAPU  = IAPU - I3*100
      I2    = INT(IAPU/10)
      IAPU  = IAPU - I2*10
      I1    = INT(IAPU)
      FILEE = CHAR(I6+48)//CHAR(I5+48)//CHAR(I4+48)//CHAR(I3+48)//
     +     CHAR(I2+48)//CHAR(I1+48)

      RETURN
      END SUBROUTINE NUMCH6
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUMCHA(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFER NUMBER ICYCLE TO THE 4 DIGIT CHARACTER FILEE

      CHARACTER(LEN=5) :: FILEE

      IAPU  = ICYCLE
      I6    = INT(IAPU/100000)
      IAPU  = IAPU - I6*100000
      I5    = INT(IAPU/10000)
      IAPU  = IAPU - I5*10000
      I4    = INT(IAPU/1000)
      IAPU  = IAPU - I4*1000
      I3    = INT(IAPU/100)
      IAPU  = IAPU - I3*100
      I2    = INT(IAPU/10)
      IAPU  = IAPU - I2*10
      I1    = INT(IAPU)
      FILEE = CHAR(I5+48)//CHAR(I4+48)//CHAR(I3+48)//CHAR(I2+48)//
     +     CHAR(I1+48)

      RETURN
      END SUBROUTINE NUMCHA
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUMCH4(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFER NUMBER ICYCLE TO THE 4 DIGIT CHARACTER FILEE

      CHARACTER(LEN=4) :: FILEE

      IAPU  = ICYCLE
      I6    = INT(IAPU/100000)
      IAPU  = IAPU - I6*100000
      I5    = INT(IAPU/10000)
      IAPU  = IAPU - I5*10000
      I4    = INT(IAPU/1000)
      IAPU  = IAPU - I4*1000
      I3    = INT(IAPU/100)
      IAPU  = IAPU - I3*100
      I2    = INT(IAPU/10)
      IAPU  = IAPU - I2*10
      I1    = INT(IAPU)
      FILEE = CHAR(I4+48)//CHAR(I3+48)//CHAR(I2+48)//CHAR(I1+48)

      RETURN
      END SUBROUTINE NUMCH4
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUMCH3(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFER NUMBER ICYCLE TO THE 3 DIGIT CHARACTER FILEE

      CHARACTER(LEN=3) :: FILEE

      IAPU  = ICYCLE
      I6    = INT(IAPU/100000)
      IAPU  = IAPU - I6*100000
      I5    = INT(IAPU/10000)
      IAPU  = IAPU - I5*10000
      I4    = INT(IAPU/1000)
      IAPU  = IAPU - I4*1000
      I3    = INT(IAPU/100)
      IAPU  = IAPU - I3*100
      I2    = INT(IAPU/10)
      IAPU  = IAPU - I2*10
      I1    = INT(IAPU)
      FILEE = CHAR(I3+48)//CHAR(I2+48)//CHAR(I1+48)

      RETURN
      END SUBROUTINE NUMCH3
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUM2ABC(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFERS TWO DIGIT NUMBER TO TWO ASCII CHARACTERS

      CHARACTER(LEN=2) :: FILEE

      IAPU  = ICYCLE

      I2    = INT(IAPU/10)
      IAPU  = IAPU - I2*10
      I1    = INT(IAPU)
      FILEE = CHAR(I2+97)//CHAR(I1+97)

      RETURN
      END SUBROUTINE NUM2ABC
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NUMCH1(FILEE,ICYCLE)

C ... THIS SUBROUTINE TRANSFER NUMBER ICYCLE TO THE 1 DIGIT CHARACTER FILEE

      CHARACTER(LEN=1) :: FILEE

      IAPU  = ICYCLE
      I1    = INT(IAPU)
      FILEE = CHAR(I1+48)

      RETURN
      END SUBROUTINE NUMCH1
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SETCON(JLOC2,JTRA2,APP2,JLOC,JTRA,APP,MREQ1)

      INTEGER :: JLOC2(*), JTRA2(*), JLOC(*), JTRA(*)
      REAL    :: APP2(*), APP(*)

      DO I = 1,MREQ1
         JLOC(I) = JLOC2(I)
         JTRA(I) = JTRA2(I)
         APP(I)  = APP2(I)
      ENDDO

      RETURN
      END SUBROUTINE SETCON
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE DOMAW(MAXW,IMAX,JMAX,IN,JN,MAW)

*     USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: MAXW, IMAX, JMAX, IN, JN, MAW

C ... optimum of T3E primary cache. Must be changed in MODIMS also if changed 

*      MAW = (IMAX+4)*(JMAX+4)! tata voisi optimoida
      MAW = (IMAX+2*IN)*(JMAX+2*JN)

      RETURN
      END SUBROUTINE DOMAW
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE MODIMS(MS)

C ... optimum of T3E primary cache. Must be changed in DOMAW also if changed 

      NTIM = INT((MS-1)/127)
      MS   = NTIM*127 + MAX0(MOD(MS-1,127),2*NTIM-1) + 1
         
      RETURN
      END SUBROUTINE MODIMS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE NATRAP(DRO,IMAX,JMAX,KMAX,IN,JN,KN,
     &                  NB,M,MERKKI,IREPEA,NREPEA)

C ... THIS SUBROUTINE FINDS NANS

      USE INTEGERS,  ONLY : IPRO

      IMPLICIT NONE

      INTEGER :: IREPEA,NREPEA,IPAIVA,IMONTH,IHOUR,IMINUT,ISECON,IFILE,
     +   ISTRID,IMAX,JMAX,KMAX,IN,JN,KN,IL,J2,JNAN,I,J,K,N,M,NB

      INTEGER :: VALUES(8)

      REAL              :: DRO(*)

      CHARACTER(LEN=5)  :: MERKKI
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      CHARACTER(LEN=3)  :: CN

C ... Control the number of NATRAP calls
c        write(6,*) irepea,nrepea,IPRO
      IF(IREPEA >= NREPEA) THEN

         WRITE(*,*) ' I am going to stop because of nans. Process=',IPRO
         WRITE(13,*)' Too many nan values found in NATRAP. Exiting...'
         WRITE(4,*) ' Too many nan values found in NATRAP. Exiting...'
         WRITE(45,*)' Too many nan values found in NATRAP. Exiting...'

         CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

         WRITE(13,*)
         WRITE(13,'(A,I2,A1,I1,A1,I4,A,I2,A1,I2,A1,I2,A)') 
     + '  Simulation is finished ',
     +    VALUES(3),'.',VALUES(2),'.',VALUES(1),' at ',
     +    VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
         
         WRITE(13,'(A)')    '  Hermit : I will be resting now... '
         WRITE(13,'(A,I3)') '  This message is from process ',IPRO
         WRITE(13,'(A)')    '  (Sorry about this trouble chief). '

c         IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR) ! Needs PARALLEL

         WRITE(13,*) 'END.'
         CLOSE(13)
*         CLOSE(IFILE)
         STOP

      ENDIF ! IREPEA >= NREPEA

      ISTRID = IMAX + 2*IN
*      IL     = ISTRID*(JMAX+4)
      IL     = ISTRID*(JMAX+2*JN)
      J2     = 0
      JNAN   = 0

      DO K = -1,KMAX+2
      DO J = -1,JMAX+2
      DO I = -1,IMAX+2
         N = 1 + (I-1+IN) + (J-1+JN)*ISTRID +  (K-1+KN)*IL
         IF(DRO(N) > 0. .or. DRO(N) <= 0.) THEN
c     WRITE(*,111) I,J,K,IMAX,JMAX,KMAX,NB,M,DRO(N),MERKKI
c            j2= j2+1 ! useless
            CONTINUE
         ELSE
            JNAN = JNAN + 1
            IF(JNAN <= 10 .AND. IREPEA <= 1) THEN
               WRITE(*,111) I,J,K,IMAX,JMAX,KMAX,NB,M,DRO(N),MERKKI
            ELSE
               IREPEA = IREPEA + 1
               IF(IREPEA == 1) THEN
               IFILE = 8900+IPRO
               CALL NUMCH3(CN,IPRO-1)
               OPEN(IFILE,FILE= 'NANS.'//CN,STATUS='UNKNOWN',
     +         FORM='FORMATTED')
               ENDIF
               WRITE(IFILE,111) I,J,K,IMAX,JMAX,KMAX,NB,M,DRO(N),MERKKI
            ENDIF
               
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      IF(JNAN > 10) THEN
      WRITE(*,*) ' More than 10 NANS found. Rest is in file NAN.8900+',
     +     JNAN
      ENDIF
 111  FORMAT("nan at ",3I3,"(",3I4,") block=",I2," M=",I1," value=",
     +     E16.8," vari=",A5)

      RETURN
      END SUBROUTINE NATRAP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE FLUXME(FLUXTY,IFLUX,ITURB,PARALLEL)

      USE MPI
      CHARACTER(LEN=80) :: FLUXTY
      LOGICAL :: PARALLEL
CCC TPS 25.2.2008      
      INTEGER :: ERRORCODE
CCC  

      IF(FLUXTY == 'van leer' .OR. FLUXTY == 'VAN LEER' .OR.
     + FLUXTY == 'van Leer' .OR. FLUXTY == 'Van Leer') THEN ! VAN LEER
         IFLUX = 0
      ELSEIF(FLUXTY == 'roe' .OR. FLUXTY == 'ROE' .OR. FLUXTY == 'Roe') 
     +        THEN
         IFLUX = 1
      ELSE IF(FLUXTY == 'wps' .OR. FLUXTY == 'WPS' .OR. FLUXTY == 'Wps'
     + ) THEN !  WAVE/PARTICLE SPLITTING SCHEME
         IFLUX = 2
      ELSE IF(FLUXTY == 'ausm'.OR.FLUXTY == 'AUSM' .OR. FLUXTY == 'Ausm'
     + ) THEN ! WAVE/PARTICLE SPLITTING SCHEME USING AUSM IMPLEMENTATION
         IFLUX = 3
      ELSE IF(FLUXTY == 'inco'.OR.FLUXTY == 'INCO' .OR. FLUXTY == 'Inco'
     + ) THEN ! APPROXIMATIVE INCOMPRESSIBLE FLUX. PSEUDOCOMPRESSIBILITY
         IFLUX = 4
      ELSE IF(FLUXTY == 'roersm'.OR.FLUXTY == 'ROERSM' .OR. 
     + FLUXTY == 'RoeRSM') THEN ! APPROXIMATIVE Roe for RSM
         IFLUX = 5
         IF(ITURB < 20) THEN
            WRITE(*,*) 'You have not activated RSM. You cannot have '
     +           ,'flux-splitting ROERSM'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
      ELSE IF(FLUXTY == 'MPHASE'.OR.FLUXTY == 'Mphase' .OR. 
     +        FLUXTY == 'mphase'
     + ) THEN ! TWO-PHASE FLOW WITH EQUAL VELOCITIES
         IFLUX = 6
      ELSE IF(FLUXTY == 'KT'.OR.FLUXTY == 'Kt' .OR. 
     +        FLUXTY == 'kt'.OR.FLUXTY == 'TADMOR'
     + ) THEN ! Kurganov-Tadmor flux
         IFLUX = 7
      ELSE IF(FLUXTY == 'H-CUSP'.OR.FLUXTY == 'H-cusp' .OR. 
     +        FLUXTY == 'h-cusp'.OR.FLUXTY == 'H-Cusp'
     + ) THEN ! Kurganov-Tadmor flux
         IFLUX = 8
      ELSE
         WRITE( 4,9160)
         WRITE( 6,9160)
         WRITE(44,9160)
         WRITE(45,9160)
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
*         CALL EXIT
         STOP
      ENDIF
 9160 FORMAT('  FLUX-TYPE  AMBIGIOUS. TRY AGAIN')

      RETURN
      END SUBROUTINE FLUXME
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C


