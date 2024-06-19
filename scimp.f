C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOCDTF(DFI,DTL,NTOT,DTG,MAXSB,NSCAL)

      DIMENSION :: DFI(MAXSB,MAX(1,NSCAL)),DTL(*)
C
C ... USE THE LOCAL TIME-STEP BEFORE THE IMPLICIT SWEEP
C
c      PDT     = 1./DTG
      DO 1100 NS= 1,NSCAL
      DO 1100 K = 1,NTOT
      DT1       = DTL(K)
      DFI(K,NS) = DT1*DFI(K,NS)
1100  CONTINUE

      RETURN
      END SUBROUTINE LOCDTF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DERISC(DFI,DTL,DT,NTOT,MAXSB,NSCAL)

      DIMENSION :: DFI(MAXSB,MAX(1,NSCAL)),DTL(*)

      EPS     = 1.E-20

      DO 1000 NS= 1,NSCAL
      DO 1000 N = 1,NTOT
      PVOL      = 1./(1.+1.5*DTL(N)/DT)
      DFI(N,NS) = DFI(N,NS) *PVOL
1000  CONTINUE

      RETURN
      END SUBROUTINE DERISC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SOUREY(DUU,UU,E,SUU,VOL,DTL,NTOT,RO,VIS,MAXSB)

      REAL :: DUU(MAXSB,6), UU(MAXSB,6), SUU(MAXSB,6), VOL(*), DTL(*),
     + RO(*),VIS(*),E(*)

c ... CALCULATE SOURCE TERM IN IMPLICIT STEP
      EPS     = 1.E-20
      CR      = 5.

      DO 1000 N  = 1,NTOT
         PVOL      = DTL(N)
         XKIN      = (UU(N,1) + UU(N,4) + UU(N,6))/3.
         UUN1      = CR/(UU(N,1)+EPS)
         UUN4      = CR/(UU(N,4)+EPS)
         UUN6      = CR/(UU(N,6)+EPS)
         UUN2      = CR/(ABS(UU(N,2))+EPS+.1*XKIN)
         UUN3      = CR/(ABS(UU(N,3))+EPS+.1*XKIN)
         UUN5      = CR/(ABS(UU(N,5))+EPS+.1*XKIN)
         ELIM      = CR/(ABS(0.1*E(N)-XKIN)+EPS)
         SOUR1    = PVOL*ABS(SUU(N,1))*MAX(UUN1,ELIM)
         SOUR2    = PVOL*ABS(SUU(N,2))*MAX(UUN2,ELIM)
         SOUR3    = PVOL*ABS(SUU(N,3))*MAX(UUN3,ELIM)
         SOUR4    = PVOL*ABS(SUU(N,4))*MAX(UUN4,ELIM)
         SOUR5    = PVOL*ABS(SUU(N,5))*MAX(UUN5,ELIM)
         SOUR6    = PVOL*ABS(SUU(N,6))*MAX(UUN6,ELIM)
         DUU(N,1) = DUU(N,1)/(1. + SOUR1) 
         DUU(N,2) = DUU(N,2)/(1. + SOUR2) 
         DUU(N,3) = DUU(N,3)/(1. + SOUR3) 
         DUU(N,4) = DUU(N,4)/(1. + SOUR4) 
         DUU(N,5) = DUU(N,5)/(1. + SOUR5) 
         DUU(N,6) = DUU(N,6)/(1. + SOUR6) 
1000  CONTINUE

      RETURN
      END SUBROUTINE SOUREY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CORRFI(DFI,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 RO,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 PR,PRT,IN,JN,KN,NSCAL,MAXSB,BDFI,MAW,MAXW,TRCASE)


      REAL :: DFI(MAXSB,MAX(1,NSCAL)), VIST(*), VIS(*), DTL(*), A2(*),
     2 VOL(*),U(*),V(*),RO(*),A2X(*),A2Y(*),A2Z(*),W(*)
      REAL :: BDFI(MAW)
      REAL :: CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1,DIFF
      INTEGER :: TRCASE
C
C ... CORRECTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C

      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

      CDIFF  = 2.
      C1     = 0.5
      EPS    = 1.E-6

C ... Intermittency diffusion coefficients
      IF(TRCASE == 1) THEN  !G
         CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
         SIGOT = 1.
      ELSEIF(TRCASE == 2) THEN  !RET
         CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
         SIGF  = 1.
      ELSE
         SIGOT = 1.
         SIGF  = 1.
      ENDIF

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 NS = 1,NSCAL

C ... STARTING FLUX
      DO 900 I = 1,JMAX*ISTRID
 900  BDFI(I) = 0.

      DO 1000 K = 1,KMAX

      IF(K <= IDI) THEN
         CDIFF = 2.
      ELSE
         CDIFF = 0.
      ENDIF

      JJ      = (K+KN-1)*IL + JN*ISTRID

      DO 1100 I = IN+1,JMAX*ISTRID,ISTEP

      N       = JJ + I
      T       = 0.
      RL1     = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N)
      RNUM = .5*(VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
      DIFF    = (  (VIST(N)   /SIGF + VIS(N)   )*SIGOT 
     2           + (VIST(N+IL)/SIGF + VIS(N+IL))*SIGOT  )
      T       = CDIFF*A2(N+IL)/RNUM*DIFF
      X1      = MAX(RL1,0.) + T
      X1A     = ABS(RL1) + T

      DT1     = THETA*DTL(N)
CC    DT1     = DTG
      ALPO    = A2(N)*BDFI(I)     ! Antakee mersu
      IF(ABS(ALPO) < 1.E-30) ALPO = 0.
      ALKU    = DT1*ALPO
      YDFI    = ALKU + VOL(N)*DFI(N,NS)
C      SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
      SEPSA   = 0.
      YDFI    = YDFI/(VOL(N) + DT1*(X1A*A2(N+IL)+SEPSA))
      IF(ABS(YDFI) < 1.E-20) YDFI = 0.
      DFI(N,NS) = YDFI
1100  BDFI(I) = X1*YDFI

1000  CONTINUE

      RETURN
      END SUBROUTINE CORRFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PREDFI(DFI,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 RO,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 PR,PRT,IN,JN,KN,NSCAL,MAXSB,BDFI,MAW,MAXW,TRCASE)

      REAL :: DFI(MAXSB,MAX(1,NSCAL)), VIST(*), VIS(*), DTL(*), A2(*),
     2 VOL(*),U(*),V(*),RO(*),A2X(*),A2Y(*),A2Z(*),W(*)
      REAL :: BDFI(MAW)
      REAL :: CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1,DIFF
      INTEGER :: TRCASE
C
C ... PREDICTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C

      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID
      CDIFF  = 2.
      C1     = 0.5
      EPS    = 1.E-6

C ... Intermittency diffusion coefficients
      IF(TRCASE == 1) THEN  !G
         CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
         SIGOT = 1.
      ELSEIF(TRCASE == 2) THEN  !RET
         CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
         SIGF  = 1.
      ELSE
         SIGOT = 1.
         SIGF  = 1.
      ENDIF

      ISTEP = 1

      IF(IMAX == 1) ISTEP = 2*IN+1
     
      DO 1000 NS = 1,NSCAL
C ... STARTING FLUX
      DO I = 1,JMAX*ISTRID
         BDFI(I) = 0.
      ENDDO

      KLOP    = 1

      DO 1000 K = KMAX,KLOP,-1

      JJ      = (K+KN-1)*IL + JN*ISTRID
      FF      = .5
      IF(K == 1) FF = 0.25
      IF(K <= IDI) THEN
         CDIFF= 2.
      ELSE
         CDIFF= 0.
      ENDIF

      DO 1100 I = IN+1,JMAX*ISTRID,ISTEP

      N       = JJ + I
      T       = 0.
      RL1     = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N)
      RNUM    = FF*(VOL(N)+VOL(N-IL))*(RO(N)+RO(N-IL))
      DIFF    = ( (VIST(N)   /SIGF + VIS(N)   )*SIGOT 
     2          + (VIST(N-IL)/SIGF + VIS(N-IL))*SIGOT )
      T       = CDIFF*A2(N)/RNUM*DIFF
      X1      = MAX(-RL1,0.) + T
      X1A     = ABS(RL1) + T
      DT1     = THETA*DTL(N)
      YDFI    = DT1*A2(N+IL)*BDFI(I) + VOL(N)*DFI(N,NS)
C      SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
      SEPSA   = 0.
      YDFI    = YDFI/(VOL(N) + DT1*(X1A*A2(N)+SEPSA))
      DFI(N,NS) = YDFI
      BDFI(I) = X1*YDFI
1100  CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE PREDFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DIAGFI(DFI,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 IMAX,JMAX,KMAX,THETA,IN,JN,KN,NSCAL,MAXSB)

      REAL :: DFI(MAXSB,MAX(1,NSCAL)), DTL(*), A2(*),
     2        VOL(*), U(*), V(*), A2X(*), A2Y(*), A2Z(*), W(*)
C
C ... MULTPLICATION WITH DIAGONAL MATRIX IN DD-ADI
C
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

      CDIFF  = 2.
      C1     = 0.5
      EPS    = 1.E-6

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 NS = 1,NSCAL

      DO 1000 K = 1,KMAX

      JJ  = (K+KN-1)*IL + JN*ISTRID

      DO 1300 I = IN+1,JMAX*ISTRID,ISTEP

      N       = JJ + I
      YDFI    = DFI(N,NS)
      DT1     = THETA*DTL(N)

      UCON    = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N)
      RL1N    = MAX(- UCON,0.)

      UCON    = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N)
      RL1P    = MAX(UCON,0.)

      DT1     = DT1/VOL(N)
C      SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
      SEPSA   = 0.
      YDFI    = YDFI*(1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N + SEPSA))
      DFI(N,NS) = YDFI
1300  CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE DIAGFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COPRFI(DFI,FI,RO,DRO,NTOT,NSCAL,MAXSB)

      REAL :: DRO(*), RO(*), DFI(MAXSB,MAX(1,NSCAL)), 
     2        FI(MAXSB,MAX(1,NSCAL))
C
C ... FROM THE CONSERVATIVE TO THE PRIMITIVE VARIABLES
C ... THIS SUBROUTINE IS USED FOR SCALAR EQ. PPR 24.2
C
      DO 1100 NS = 1,NSCAL
      DO 1000 I = 1,NTOT
      YPRO    = 1./RO(I)
      DFI(I,NS) = YPRO*(-FI(I,NS)*YPRO*DRO(I) + DFI(I,NS))
 1000 CONTINUE
 1100 CONTINUE

      RETURN
      END SUBROUTINE COPRFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRCOFI(DFI,FI,RO,DRO,NTOT,NSCAL,MAXSB)

      REAL :: DRO(*), RO(*), DFI(MAXSB,MAX(1,NSCAL)), 
     2        FI(MAXSB,MAX(1,NSCAL))
C
C ... FROM THE PRIMITIVE TO THE CONSERVATIVE VARIABLES 
C ... THIS SUBROUTINE IS USED FOR SCALAR EQ. 24.2
C
      DO 1000 NS = 1,NSCAL
      DO 1000 I = 1,NTOT
      YPRO    = 1./RO(I)
      DFI(I,NS) = YPRO*FI(I,NS)*DRO(I) + RO(I)*DFI(I,NS)
1000  CONTINUE

      RETURN
      END SUBROUTINE PRCOFI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C*******************************************************
C*                                                     *
C*    RSM IMPLICID SWEEP SUBROUTINES                   *
C*                                                     *
C*******************************************************

      SUBROUTINE COPRRE(DRO,DM,DN,DW,DE,DRK,DEPS,DFI,RO,E,U,V,W,
     +                  RK,REPS,FI,NTOT,MAXSB)

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),RO(*),U(*),V(*),W(*),E(*),
     +         DRK(*),DEPS(*),RK(*),REPS(*),DFI(MAXSB,6),FI(MAXSB,6)
C
C ... FROM THE CONSERVATIVE TO THE PRIMITIVE VARIABLES
C ... THIS SUBROUTINE IS USED FOR RSM
C
      DO 1000 I = 1,NTOT
      BETA     = .5*(U(I)**2 + V(I)**2 + W(I)**2)
      YPRO     = 1./RO(I)
      DMI      = DM(I)
      DNI      = DN(I)
      DWI      = DW(I)
      DKI      = DRK(I)
      DM(I)    = YPRO*(-U(I)*DRO(I) + DMI)
      DN(I)    = YPRO*(-V(I)*DRO(I) + DNI)
      DW(I)    = YPRO*(-W(I)*DRO(I) + DWI)
      DRK(I)   = YPRO*(-RK(I)*YPRO*DRO(I) + DKI)
      DEPS(I)  = YPRO*(-REPS(I)*YPRO*DRO(I) + DEPS(I))
      DE(I)    = YPRO*(- U(I)*DMI - V(I)*DNI - W(I)*DWI + DE(I) +
     2         (-E(I)*YPRO + 2.*BETA + RK(I)*YPRO)*DRO(I) - DKI)
      DFI(I,1) = YPRO*(-FI(I,1)*YPRO*DRO(I) + DFI(I,1))
      DFI(I,2) = YPRO*(-FI(I,2)*YPRO*DRO(I) + DFI(I,2))
      DFI(I,3) = YPRO*(-FI(I,3)*YPRO*DRO(I) + DFI(I,3))
      DFI(I,4) = YPRO*(-FI(I,4)*YPRO*DRO(I) + DFI(I,4))
      DFI(I,5) = YPRO*(-FI(I,5)*YPRO*DRO(I) + DFI(I,5))
      DFI(I,6) = YPRO*(-FI(I,6)*YPRO*DRO(I) + DFI(I,6))
1000  CONTINUE

      RETURN
      END SUBROUTINE COPRRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRCORE(DRO,DM,DN,DW,DE,DRK,DEPS,DFI,RO,E,P,PDIFF,U,V,W,
     +     C,RK,REPS,FI,XC,YC,ZC,FRSVEL,FRSDEN,NTOT,MAXSB)

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),RO(*),U(*),V(*),W(*),E(*),
     +        DRK(*),DEPS(*),RK(*),REPS(*),P(*),PDIFF(*),C(*),
     +        XC(*),YC(*),ZC(*),DFI(MAXSB,6),FI(MAXSB,6)
C
C ... FROM THE PRIMITIVE TO THE CONSERVATIVE VARIABLES
C ... THIS SUBROUTINE IS USED FOR RSM
C
      DO 1000 I = 1,NTOT
      YPRO    = 1./RO(I)
      BETA    = .5*(U(I)**2 + V(I)**2 + W(I)**2) + 
     +     .5*(FI(I,1) + FI(I,4) + FI(I,6))*YPRO
      DMI     = DM(I)
      DNI     = DN(I)
      DWI     = DW(I)
      DKI     = DRK(I)
      DKIR    = .5*(DFI(I,1) + DFI(I,4) + DFI(I,6))
      DM(I)   = U(I)*DRO(I) + RO(I)*DMI
      DN(I)   = V(I)*DRO(I) + RO(I)*DNI
      DW(I)   = W(I)*DRO(I) + RO(I)*DWI
      DRK(I)  = YPRO*RK(I)*DRO(I) + RO(I)*DKI
      DEPS(I) = YPRO*REPS(I)*DRO(I) + RO(I)*DEPS(I)
      DE(I)   = E(I)/RO(I)*DRO(I) + RO(I)*(U(I)*DMI +
     2          V(I)*DNI + W(I)*DWI + DKIR + DE(I))
      DFI(I,1) = YPRO*FI(I,1)*DRO(I) + RO(I)*DFI(I,1)
      DFI(I,2) = YPRO*FI(I,2)*DRO(I) + RO(I)*DFI(I,2)
      DFI(I,3) = YPRO*FI(I,3)*DRO(I) + RO(I)*DFI(I,3)
      DFI(I,4) = YPRO*FI(I,4)*DRO(I) + RO(I)*DFI(I,4)
      DFI(I,5) = YPRO*FI(I,5)*DRO(I) + RO(I)*DFI(I,5)
      DFI(I,6) = YPRO*FI(I,6)*DRO(I) + RO(I)*DFI(I,6)
1000  CONTINUE

      RETURN
      END SUBROUTINE PRCORE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PREDRE(DRO,DM,DN,DW,DE,DRK,DEPS,DUU,DUV,DUW,DVV,DVW,
     2     DWW,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,C,RO,P,RK,REPS,SRK,SEPS,
     3     RUU,RUV,RUW,RVV,RVW,RWW,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,
     4     VIST,DRDP,DRDH,PR,PRT,IN,JN,KN,UROT,ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),A2(*),
     2     VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3     W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),
     4     DUU(*),DUV(*),DUW(*),DVV(*),DVW(*),DWW(*),
     5     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     6     DRDP(*),DRDH(*),P(*)
      REAL, POINTER :: WDRO(:),WDM(:),WDN(:),WDW(:),WDE(:),WDK(:),
     +     WDEP(:),WDUU(:),WDUV(:),WDUW(:),WDVV(:),WDVW(:),WDWW(:)
      REAL, TARGET ::  ZZZ(MAXW)

      WDRO => ZZZ( 0*MAW+1: 1*MAW);WDM => ZZZ( 1*MAW+1: 2*MAW)
      WDN  => ZZZ( 2*MAW+1: 3*MAW);WDW => ZZZ( 3*MAW+1: 4*MAW)
      WDE  => ZZZ( 4*MAW+1: 5*MAW);WDK => ZZZ( 5*MAW+1: 6*MAW)
      WDEP => ZZZ( 6*MAW+1: 7*MAW);WDUU=> ZZZ( 7*MAW+1: 8*MAW)
      WDUV => ZZZ( 8*MAW+1: 9*MAW);WDUW=> ZZZ( 9*MAW+1:10*MAW)
      WDVV => ZZZ(10*MAW+1:11*MAW);WDVW=> ZZZ(11*MAW+1:12*MAW)
      WDWW => ZZZ(12*MAW+1:13*MAW)

C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... PREDICTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C
C ... ORDER OF EIGENVALUES
C ...        RO U V W E RK REPS UU UV UW VW WW
C ...   X    1  2 5 5 3  1   1   1 4  4  1  1
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ISTRID  = IMAX + 2*IN
      IL      = JSTRID*ISTRID

      CDIFF   = 2.
      C1      = 0.5
      EPS     = 1.E-6
C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
         WDRO(I) = 0.
         WDM(I)  = 0.
         WDN(I)  = 0.
         WDW(I)  = 0.
         WDE(I)  = 0.
         WDK(I)  = 0.
         WDEP(I) = 0.
         WDUU(I) = 0.
         WDUV(I) = 0.
         WDUW(I) = 0.
         WDVV(I) = 0.
         WDVW(I) = 0.
         WDWW(I) = 0.
 900  CONTINUE

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      KLOP = 1

      DO 1000 K = KMAX,KLOP,-1

      JJ = (K+KN-1)*IL + JN*ISTRID

      FF      = .5
      IF(K == 1) FF = 0.25
      F1 = 0.
      IF(K == 1) F1 = 1.

      DO 1400 I = IN+1,JMAX*ISTRID,ISTEP

         N       = JJ + I
         RUUN    = A2X(N)**2*RUU(N) + A2Y(N)**2*RVV(N) + 
     +        A2Z(N)**2*RWW(N) - 2.*(A2X(N)*A2Y(N)*RUV(N) + 
     +        A2X(N)*A2Z(N)*RUW(N) + A2Y(N)*A2Z(N)*RVW(N))
         UUN     = SQRT(MAX(RUUN/RO(N),1.E-20))
         DT1     = THETA*DTL(N)
CC       DT1     = DTG
         WIRO = DT1*A2(N+IL)*WDRO(I) + VOL(N)*DRO(N)
         WIM = DT1*A2(N+IL)*WDM(I)  + VOL(N)*DM(N)
         WIN = DT1*A2(N+IL)*WDN(I)  + VOL(N)*DN(N)
         WIW = DT1*A2(N+IL)*WDW(I)  + VOL(N)*DW(N)
         WIE = DT1*A2(N+IL)*WDE(I)  + VOL(N)*DE(N)
         WIK = DT1*A2(N+IL)*WDK(I)  + VOL(N)*DRK(N)
         WIEP = DT1*A2(N+IL)*WDEP(I) + VOL(N)*DEPS(N)
         WIUU = DT1*A2(N+IL)*WDUU(I) + VOL(N)*DUU(N)
         WIUV = DT1*A2(N+IL)*WDUV(I) + VOL(N)*DUV(N)
         WIUW = DT1*A2(N+IL)*WDUW(I) + VOL(N)*DUW(N)
         WIVV = DT1*A2(N+IL)*WDVV(I) + VOL(N)*DVV(N)
         WIVW = DT1*A2(N+IL)*WDVW(I) + VOL(N)*DVW(N)
         WIWW = DT1*A2(N+IL)*WDWW(I) + VOL(N)*DWW(N)
         DT   = DT1
1100  CONTINUE

      CALL SYIRE
     1 (YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     2  WIRO,WIM,WIN,WIW,WIE,WIK,WIEP,WIUU,WIUV,WIUW,WIVV,WIVW,WIWW,
     3     C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4     A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

         UC      = UROT(N)
         RL1     = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
         RL2     = RL1  + C(N)
         RL3     = RL1  - C(N)
         RL4     = RL1  + UUN
         RL5     = RL1  - UUN
         RNUM    = FF*(VOL(N)+VOL(N-IL)+F1)*(RO(N)+RO(N-IL))
         T       = CDIFF*A2(N)/RNUM*
     2        (VIST(N)+VIS(N) + VIST(N-IL)+VIS(N-IL))
         
         X1      = MAX(-RL1,0.) + T
         X2      = MAX(-RL2,0.) + T
         X3      = MAX(-RL3,0.) + T
         X4      = MAX(-RL4,0.) + T
         X5      = MAX(-RL5,0.) + T
         X1A     = ABS(RL1) + T
         X2A     = ABS(RL2) + T
         X3A     = ABS(RL3) + T
         X4A     = ABS(RL4) + T
         X5A     = ABS(RL5) + T

C        SRKA    = ABS(SRK(N))/(RK(N)+EPS)*C1*0.
C        SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
         SRKA    = 0.
         SEPSA   = 0.
         YDRO    = YDRO/(VOL(N) + DT*X1A*A2(N))
         YDM     = YDM /(VOL(N) + DT*X2A*A2(N))
         YDN     = YDN /(VOL(N) + DT*X5A*A2(N))
         YDW     = YDW /(VOL(N) + DT*X5A*A2(N))
         YDE     = YDE /(VOL(N) + DT*X3A*A2(N))
         YDK     = YDK /(VOL(N) + DT*(X1A*A2(N)+SRKA))
         YDEP    = YDEP/(VOL(N) + DT*(X1A*A2(N)+SEPSA))
         YDUU    = YDUU/(VOL(N) + DT*X1A*A2(N))
         YDUV    = YDUV/(VOL(N) + DT*X4A*A2(N))
         YDUW    = YDUW/(VOL(N) + DT*X4A*A2(N))
         YDVV    = YDVV/(VOL(N) + DT*X1A*A2(N))
         YDVW    = YDVW/(VOL(N) + DT*X1A*A2(N))
         YDWW    = YDWW/(VOL(N) + DT*X1A*A2(N))


      CALL SYM1RE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),
     1        DUU(N),DUV(N),DUW(N),DVV(N),DVW(N),DWW(N),
     2  YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4        A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

         YDRO = X1*YDRO
         YDM  = X2*YDM
         YDN  = X5*YDN
         YDW  = X5*YDW
         YDE  = X3*YDE
         YDK  = X1*YDK
         YDEP = X1*YDEP
         YDUU = X1*YDUU
         YDUV = X4*YDUV
         YDUW = X4*YDUW
         YDVV = X1*YDVV
         YDVW = X1*YDVW
         YDWW = X1*YDWW

      CALL SYM1RE
     1 (WIRO,WIM,WIN,WIW,WIE,WIK,WIEP,WIUU,WIUV,WIUW,WIVV,WIVW,WIWW,
     2  YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4        A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

         WDRO(I) = WIRO
         WDM(I)  = WIM 
         WDN(I)  = WIN 
         WDW(I)  = WIW 
         WDE(I)  = WIE 
         WDK(I)  = WIK 
         WDEP(I) = WIEP
         WDUU(I) = WIUU
         WDUV(I) = WIUV
         WDUW(I) = WIUW
         WDVV(I) = WIVV
         WDVW(I) = WIVW
         WDWW(I) = WIWW
 1400 CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE PREDRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DIAGRE(DRO,DM,DN,DW,DE,DRK,DEPS,DUU,DUV,DUW,DVV,DVW,
     2     DWW,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,C,RO,P,RK,REPS,SRK,SEPS,
     3     RUU,RUV,RUW,RVV,RVW,RWW,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,
     4     VIST,DRDP,DRDH,PR,PRT,IN,JN,KN,UROT,ZZ2,MA2,MA2W)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),A2(*),
     2     VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3     W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),
     4     DUU(*),DUV(*),DUW(*),DVV(*),DVW(*),DWW(*),
     5     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     6     DRDP(*),DRDH(*),P(*)

C ... CONTRAVARIANTTI MUOTO ON HIEMAN KAKSMIELINEN PPR 23.5.95
C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... MULTPLICATION WITH DIAGONAL MATRIX IN DD-ADI
C
C ... ORDER OF EIGENVALUES
C ...        RO U V W E RK REPS UU UV UW VW WW
C ...   X    1  2 5 5 3  1   1   1 4  4  1  1
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

      CDIFF  = 2.
      C1     = 0.5
      EPS    = 1.E-6

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 K = 1,KMAX

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1100 I = IN+1,JMAX*ISTRID,ISTEP

         N       = JJ + I
         RUUN    = A2X(N+IL)**2*RUU(N) + A2Y(N+IL)**2*RVV(N) + 
     +        A2Z(N+IL)**2*RWW(N) - 2.*(A2X(N+IL)*A2Y(N+IL)*RUV(N) + 
     +        A2X(N+IL)*A2Z(N+IL)*RUW(N) + A2Y(N+IL)*A2Z(N+IL)*RVW(N))
         UUN  = SQRT(MAX(RUUN/RO(N),1.E-20))
         DT   = THETA*DTL(N)

      CALL SYIRE
     1 (YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     2        DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),
     +        DUU(N),DUV(N),DUW(N),DVV(N),DVW(N),DWW(N),
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4        A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

         UC      = UROT(N+IL)
         UCON    = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
         RL1N    = MAX(- UCON,0.)
         RL2N    = MAX(-(UCON+C(N)), 0.)
         RL3N    = MAX(-(UCON-C(N)), 0.)
         RL4N    = MAX(-(UCON+UUN), 0.)
         RL5N    = MAX(-(UCON-UUN), 0.)
         
         UC      = UROT(N)
         RUUN    = A2X(N)**2*RUU(N) + A2Y(N)**2*RVV(N) + 
     +        A2Z(N)**2*RWW(N) - 2.*(A2X(N)*A2Y(N)*RUV(N) + 
     +        A2X(N)*A2Z(N)*RUW(N) + A2Y(N)*A2Z(N)*RVW(N))
         UUP     = SQRT(MAX(RUUN/RO(N),1.E-20))
         UCON    = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
         RL1P    = MAX(UCON,0.)
         RL2P    = MAX(UCON+C(N), 0.)
         RL3P    = MAX(UCON-C(N), 0.)
         RL4P    = MAX(UCON+UUP, 0.)
         RL5P    = MAX(UCON-UUP, 0.)
         
         DT1     = DT/VOL(N)
C        SRKA    = ABS(SRK(N))/(RK(N)+EPS)*C1*0.
C        SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
         SRKA    = 0.
         SEPSA   = 0.
         YDRO = YDRO*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N))
         YDM  = YDM *(1.+DT1*(A2(N)*RL2P+A2(N+IL)*RL2N))
         YDN  = YDN *(1.+DT1*(A2(N)*RL5P+A2(N+IL)*RL5N))
         YDW  = YDW *(1.+DT1*(A2(N)*RL5P+A2(N+IL)*RL5N))
         YDE  = YDE *(1.+DT1*(A2(N)*RL3P+A2(N+IL)*RL3N))
         YDK  = YDK *(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N + SRKA))
         YDEP = YDEP*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N + SEPSA))
         YDUU = YDUU*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N))
         YDUV = YDUV*(1.+DT1*(A2(N)*RL4P+A2(N+IL)*RL4N))
         YDUW = YDUW*(1.+DT1*(A2(N)*RL4P+A2(N+IL)*RL4N))
         YDVV = YDVV*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N))
         YDVW = YDVW*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N))
         YDWW = YDWW*(1.+DT1*(A2(N)*RL1P+A2(N+IL)*RL1N))

      CALL SYM1RE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),
     +        DUU(N),DUV(N),DUW(N),DVV(N),DVW(N),DWW(N),
     2  YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4        A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

1100  CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE DIAGRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CORRRE(DRO,DM,DN,DW,DE,DRK,DEPS,DUU,DUV,DUW,DVV,DVW,
     2     DWW,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,C,RO,P,RK,REPS,SRK,SEPS,
     3     RUU,RUV,RUW,RVV,RVW,RWW,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,
     4     VIST,DRDP,DRDH,PR,PRT,IN,JN,KN,UROT,ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),A2(*),
     2 VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),
     4     DUU(*),DUV(*),DUW(*),DVV(*),DVW(*),DWW(*),
     5     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     6     DRDP(*),DRDH(*),P(*)

      REAL, POINTER :: WDRO(:),WDM(:),WDN(:),WDW(:),WDE(:),WDK(:),
     +     WDEP(:),WDUU(:),WDUV(:),WDUW(:),WDVV(:),WDVW(:),WDWW(:)
      REAL, TARGET ::  ZZZ(MAXW)

      WDRO => ZZZ( 0*MAW+1: 1*MAW);WDM => ZZZ( 1*MAW+1: 2*MAW)
      WDN  => ZZZ( 2*MAW+1: 3*MAW);WDW => ZZZ( 3*MAW+1: 4*MAW)
      WDE  => ZZZ( 4*MAW+1: 5*MAW);WDK => ZZZ( 5*MAW+1: 6*MAW)
      WDEP => ZZZ( 6*MAW+1: 7*MAW);WDUU=> ZZZ( 7*MAW+1: 8*MAW)
      WDUV => ZZZ( 8*MAW+1: 9*MAW);WDUW=> ZZZ( 9*MAW+1:10*MAW)
      WDVV => ZZZ(10*MAW+1:11*MAW);WDVW=> ZZZ(11*MAW+1:12*MAW)
      WDWW => ZZZ(12*MAW+1:13*MAW)

C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... CORRECTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C
C ... ORDER OF EIGENVALUES
C ...        RO U V W E RK REPS UU UV UW VW WW
C ...   X    1  2 5 5 3  1   1   1 4  4  1  1

      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ISTRID  = IMAX + 2*IN
      IL      = JSTRID*ISTRID

      CDIFF   = 2.
      C1      = 0.5
      EPS     = 1.E-6

C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
         WDRO(I) = 0.
         WDM(I)  = 0.
         WDN(I)  = 0.
         WDW(I)  = 0.
         WDE(I)  = 0.
         WDK(I)  = 0.
         WDEP(I) = 0.
         WDUU(I) = 0.
         WDUV(I) = 0.
         WDUW(I) = 0.
         WDVV(I) = 0.
         WDVW(I) = 0.
         WDWW(I) = 0.
 900  CONTINUE

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 K = 1,KMAX

      IF(K <= IDI) THEN
         CDIFF = 2.
      ELSE
         CDIFF = 0.
      ENDIF

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1400 I = IN+1,JMAX*ISTRID,ISTEP

         N       = JJ + I
         DT1     = THETA*DTL(N)
CC       DT1     = DTG
         RUUN    = A2X(N+IL)**2*RUU(N) + A2Y(N+IL)**2*RVV(N) + 
     +        A2Z(N+IL)**2*RWW(N) - 2.*(A2X(N+IL)*A2Y(N+IL)*RUV(N) + 
     +        A2X(N+IL)*A2Z(N+IL)*RUW(N) + A2Y(N+IL)*A2Z(N+IL)*RVW(N))
         UUN  = SQRT(MAX(RUUN/RO(N),1.E-20))
         WIRO = DT1*A2(N)*WDRO(I) + VOL(N)*DRO(N)
         WIM  = DT1*A2(N)*WDM(I)  + VOL(N)*DM(N)
         WIN  = DT1*A2(N)*WDN(I)  + VOL(N)*DN(N)
         WIW  = DT1*A2(N)*WDW(I)  + VOL(N)*DW(N)
         WIE  = DT1*A2(N)*WDE(I)  + VOL(N)*DE(N)
         WIK  = DT1*A2(N)*WDK(I)  + VOL(N)*DRK(N)
         WIEP = DT1*A2(N)*WDEP(I) + VOL(N)*DEPS(N)
         WIUU = DT1*A2(N)*WDUU(I) + VOL(N)*DUU(N)
         WIUV = DT1*A2(N)*WDUV(I) + VOL(N)*DUV(N)
         WIUW = DT1*A2(N)*WDUW(I) + VOL(N)*DUW(N)
         WIVV = DT1*A2(N)*WDVV(I) + VOL(N)*DVV(N)
         WIVW = DT1*A2(N)*WDVW(I) + VOL(N)*DVW(N)
         WIWW = DT1*A2(N)*WDWW(I) + VOL(N)*DWW(N)
         DT   = DT1

      CALL SYIRE
     1 (YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     2  WIRO,WIM,WIN,WIW,WIE,WIK,WIEP,WIUU,WIUV,WIUW,WIVV,WIVW,WIWW,
     3  C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4  A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)


         UC      = UROT(N+IL)
         RL1     = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
         RL2     = RL1  + C(N)
         RL3     = RL1  - C(N)
         RL4     = RL1  + UUN
         RL5     = RL1  - UUN
         RNUM    = .5*(VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
         T       = CDIFF*A2(N+IL)/RNUM*
     2        (VIST(N)+VIS(N) + VIST(N+IL)+VIS(N+IL))
         
         X1      = MAX(RL1,0.) + T
         X2      = MAX(RL2,0.) + T
         X3      = MAX(RL3,0.) + T
         X4      = MAX(RL4,0.) + T
         X5      = MAX(RL5,0.) + T
         X1A     = ABS(RL1) + T
         X2A     = ABS(RL2) + T
         X3A     = ABS(RL3) + T
         X4A     = ABS(RL4) + T
         X5A     = ABS(RL5) + T
C        SRKA    = ABS(SRK(N))/(RK(N)+EPS)*C1*0.
C        SEPSA   = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
         SRKA    = 0.
         SEPSA   = 0.
         YDRO    = YDRO/(VOL(N) + DT*X1A*A2(N+IL))
         YDM     = YDM/ (VOL(N) + DT*X2A*A2(N+IL))
         YDN     = YDN/ (VOL(N) + DT*X5A*A2(N+IL))
         YDW     = YDW/ (VOL(N) + DT*X5A*A2(N+IL))
         YDE     = YDE/ (VOL(N) + DT*X3A*A2(N+IL))
         YDK     = YDK/ (VOL(N) + DT*(X1A*A2(N+IL)+SRKA))
         YDEP    = YDEP/(VOL(N) + DT*(X1A*A2(N+IL)+SEPSA))
         YDUU    = YDUU/(VOL(N) + DT*X1A*A2(N+IL))
         YDUV    = YDUV/(VOL(N) + DT*X4A*A2(N+IL))
         YDUW    = YDUW/(VOL(N) + DT*X4A*A2(N+IL))
         YDVV    = YDVV/(VOL(N) + DT*X1A*A2(N+IL))
         YDVW    = YDVW/(VOL(N) + DT*X1A*A2(N+IL))
         YDWW    = YDWW/(VOL(N) + DT*X1A*A2(N+IL))

      CALL SYM1RE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),
     +        DUU(N),DUV(N),DUW(N),DVV(N),DVW(N),DWW(N),YDRO,
     2        YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     +        A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

         
         YDRO = X1*YDRO
         YDM  = X2*YDM
         YDN  = X5*YDN
         YDW  = X5*YDW
         YDE  = X3*YDE
         YDK  = X1*YDK
         YDEP = X1*YDEP
         YDUU = X1*YDUU
         YDUV = X4*YDUV
         YDUW = X4*YDUW
         YDVV = X1*YDVV
         YDVW = X1*YDVW
         YDWW = X1*YDWW

      CALL SYM1RE
     1 (WIRO,WIM,WIN,WIW,WIE,WIK,WIEP,WIUU,WIUV,WIUW,WIVV,WIVW,WIWW,
     2  YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,YDUU,YDUV,YDUW,YDVV,YDVW,YDWW,
     3        C(N),RO(N),P(N),DRDP(N),DRDH(N),UUN,
     4        A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

         WDRO(I) = WIRO
         WDM(I)  = WIM 
         WDN(I)  = WIN 
         WDW(I)  = WIW 
         WDE(I)  = WIE 
         WDK(I)  = WIK 
         WDEP(I) = WIEP
         WDUU(I) = WIUU
         WDUV(I) = WIUV
         WDUW(I) = WIUW
         WDVV(I) = WIVV
         WDVW(I) = WIVW
         WDWW(I) = WIWW
 1400 CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE CORRRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SYM1RE(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,
     2                  W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,
     3     C,RO,P,DRDP,DRDH,UUN,A2X,A2Y,A2Z,A21,A22,A23,A31,A32,A33)
C
C ... MULTIPLICATION OF Y BY SYM1 (INDEPENDENT OF THE STATE EQUATION)
C
      EPS = 1.E-10
C
C ... FROM THE CHARACTERISTIC TO THE CONTRAVARIANT VARIABLES
C
C ... THERMODYNAMIC VARIABLES
      DEDP1   = -DRDP/DRDH - 1./RO
      DEDR1   =       1./DRDH + P/RO**2
      DPDE    =  1./DEDP1
      DPDRO   = -DEDR1/DEDP1
      SUUN    = UUN**2

      PC2     = .5/C**2
      PC1     = RO*(C**2-2.*SUUN)*2.*PC2
      PROC    = .5/(RO*C)
      DPPU    = DPDRO + SUUN
      X1      = W1 + PC2*(W2 + W5) - PC1/(DPPU + EPS)*W8
      X2I     = PROC*(W2 - W5)
c      ap3     = w9/(UUN+EPS)
c      ap4     = w10/(UUN+EPS)
c      if (abs(ap3) >= .01) write(*,*) 'x3',i,ap3,UUN,w9
c      if (abs(ap4) >= .01) write(*,*) 'x4',i,ap4,UUN,w10
      X3I     = .5*(W3 +  W9/(UUN+EPS))
      X4I     = .5*(W4 + W10/(UUN+EPS))
      X5   = 1./DPDE*(-(DPDRO+SUUN)*W1 +
     +     (.5 -PC2*(DPDRO+3.*SUUN))*(W2+W5))
      X6   = W6
      X7   = W7
      C11     = 2.*PC2*SUUN/(RO + EPS)*(W2 + W5 - 
     +     2.*RO*W8) + W8
      C12     = .5*(-UUN*W3 + W9)
      C13     = .5*(-UUN*W4 + W10)
      C22     = W11
      C23     = W12
      C33     = W13
C ... STRESS TENSOR IS SYMMETRIC
      C21     = C12
      C31     = C13
      C32     = C23

C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES
C

      A11 = A2X
      A12 = A2Y
      A13 = A2Z

      X2   = A11*X2I + A21*X3I + A31*X4I
      X3   = A12*X2I + A22*X3I + A32*X4I
      X4   = A13*X2I + A23*X3I + A33*X4I

      D11     =   A11*C11 - A21*C21 - A31*C31
      D12     = - A11*C12 + A21*C22 - A31*C32
      D13     = - A11*C13 - A21*C23 + A31*C33
      D21     =   A12*C11 - A22*C21 - A32*C31
      D22     = - A12*C12 + A22*C22 - A32*C32
      D23     = - A12*C13 - A22*C23 + A32*C33
      D31     =   A13*C11 - A23*C21 - A33*C31
      D32     = - A13*C12 + A23*C22 - A33*C32
      D33     = - A13*C13 - A23*C23 + A33*C33

      X8   =   D11*A11 + D12*A21 + D13*A31
      X9   = -(D11*A12 + D12*A22 + D13*A32)
      X10  = -(D11*A13 + D12*A23 + D13*A33)
      X11  =   D21*A12 + D22*A22 + D23*A32
      X12  = -(D21*A13 + D22*A23 + D23*A33)
      X13  =   D31*A13 + D32*A23 + D33*A33

1000  CONTINUE

      RETURN
      END SUBROUTINE SYM1RE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SYIRE(W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,
     2                 X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,
     3     C,RO,P,DRDP,DRDH,UUN,A2X,A2Y,A2Z,A21,A22,A23,A31,A32,A33)
C
C ... MULTIPLICATION OF Y BY SY (INDEPENDENT OF THE STATE EQUATION)
C
      EPS = 1.E-10

      Z1      = X1
      Z5      = X5
      Z6      = X6
      Z7      = X7
C
C ... FROM THE PRIMITIVE (RO,U,V,W,P) TO THE CONTRAVARIANT VARIABLES
C
      A11      = A2X
      A12      = A2Y
      A13      = A2Z

      CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

*      A21      = A2Z - A2Y
*      A22      = A2X - A2Z
*      A23      = A2Y - A2X
*      PALPO    = 1./SQRT(A21**2 + A22**2 + A23**2)
*      A21      = A21*PALPO
*      A22      = A22*PALPO
*      A23      = A23*PALPO
*      A31      = A2Y*A23 - A2Z*A22
*      A32      = A2Z*A21 - A2X*A23
*      A33      = A2X*A22 - A2Y*A21

      Z2       = A11*X2 + A12*X3 + A13*X4
      Z3       = A21*X2 + A22*X3 + A23*X4
      Z4       = A31*X2 + A32*X3 + A33*X4

C ... ROTATE REYNOLDS STRESS TENSOR
      
      U11 = X8
      U12 = X9
      U13 = X10
      U22 = X11
      U23 = X12
      U33 = X13
      B11 =   A11*U11 - A12*U12 - A13*U13
      B12 = - A11*U12 + A12*U22 - A13*U23
      B13 = - A11*U13 - A12*U23 + A13*U33
      B21 =   A21*U11 - A22*U12 - A23*U13
      B22 = - A21*U12 + A22*U22 - A23*U23
      B23 = - A21*U13 - A22*U23 + A23*U33
      B31 =   A31*U11 - A32*U12 - A33*U13
      B32 = - A31*U12 + A32*U22 - A33*U23
      B33 = - A31*U13 - A32*U23 + A33*U33
      
      Z8  =   B11*A11 + B12*A12 + B13*A13
      Z9  = -(B11*A21 + B12*A22 + B13*A23)
      Z10 = -(B11*A31 + B12*A32 + B13*A33)
      Z11 =   B21*A21 + B22*A22 + B23*A23
      Z12 = -(B21*A31 + B22*A32 + B23*A33)
      Z13 =   B31*A31 + B32*A32 + B33*A33

C
C ... FROM THE CONTRAVARIANT TO THE CHARACTERISTIC VARIABLES
C
C ... THERMODYNAMIC VARIABLES
      DEDP1   = -DRDP/DRDH - 1./RO
      DEDR1   =       1./DRDH + P/RO**2
      DPDE    =  1./DEDP1
      DPDRO   = -DEDR1/DEDP1
      SUUN    = UUN**2
      DPPU    = DPDRO + SUUN
      YPCRO   = (DPDRO+3.*SUUN)/(C**2*(DPPU + EPS))
      PC3     = 1./(C**2-2.*SUUN)
      ROC     = RO*C
      TKIN    = DPPU*Z1 + DPDE*Z5 + RO*Z8
      W1   = Z1 -  YPCRO*TKIN + RO/(DPPU + EPS)*Z8
      W2   = ROC*Z2 + TKIN
      W3   = Z3 -  Z9/(UUN+EPS)
      W4   = Z4 - Z10/(UUN+EPS)
      W5   = -ROC*Z2 + TKIN
      W6   = Z6
      W7   = Z7
      W8   = -2.*SUUN/RO*PC3*TKIN + PC3*C**2*Z8
      W9   = Z9  + UUN*Z3
      W10  = Z10 + UUN*Z4
      W11  = Z11
      W12  = Z12
      W13  = Z13

1000  CONTINUE

      RETURN
      END SUBROUTINE SYIRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
