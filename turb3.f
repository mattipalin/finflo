C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C						
      SUBROUTINE TURBBL(EPS2,VIST,VIS,OHMI,U,V,W,RO,D3,Y,F,
     2 PR,PRT,IMAX,JMAX,KMAX,IN,JN,KN,KBOT,KTOP,KCP,M,TURLIM,
     3 UBK,VBK,WBK,UTK,VTK,WTK,nbl,ZZZ,MAW,MAXW,RKSI)

      DIMENSION :: U(*),RO(*),EPS2(*),VIST(*),V(*),VIS(*),W(*),
     2             D3(*),KCP(*),OHMI(*),Y(*),F(*),
     3             UBK(*),VBK(*),WBK(*),UTK(*),VTK(*),WTK(*),RKSI(*)

      REAL, POINTER ::  UMAX1(:),UMIN1(:),FMAX1(:),SCALE(:),YMAX1(:)
      REAL, TARGET  ::  ZZZ(MAXW)

      UMAX1=> ZZZ( 0*MAW+1: 1*MAW);UMIN1=> ZZZ( 1*MAW+1: 2*MAW)
      FMAX1=> ZZZ( 2*MAW+1: 3*MAW);SCALE=> ZZZ( 3*MAW+1: 4*MAW)
      YMAX1=> ZZZ( 4*MAW+1: 5*MAW)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      IA1     = (KN-1)*IL + JN*ISTRID
      IA2     = IA1 + 1
      IA3     = IA1 + JMAX*ISTRID
      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
C ... TEST FOR KUPP
      IF(KBOT == 0) KUPP = 0
      EPS = 1.E-10
      YKOE    = 1.
      YMAXP   = 0.
      FMAXP   = 0.
      PRS     = PR/PRT
C
C ... CALCULATE TURBULENT VISCOSITIES BASED ON THE BALDWIN-LOMAX MODEL
C
C ... FIND MAXIMUM VALUES OF UTOT AND F(Y)
      DO 300 IG =1,JMAX*ISTRID
         FMAX1(IG) = 0.
         YMAX1(IG) = 0.
         UMAX1(IG) = 0.
         UMIN1(IG) = 1.E10
 300  CONTINUE

C  -- TURBULENCE CAN BE EVALUATED ONLY NEAR THE WALLS ----------------|
C     IF (KCP(IG) /= 0) THEN
C     ENDIF
C  -------------------------------------------------------------------|

      DO 100 IG  = 1,JMAX*ISTRID
C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...
C ... CALCULATE F(Y) AND HIGHT OF THE LAYER Y(X,Y,Z)
      L          = IG + IA1 + IL
      IGW        = IG + JN*ISTRID
      Y(L)       = .5*D3(L)
      SCALE(IG)  = SQRT(RO(L)*OHMI(L)/VIS(L))
      YEXP       = MAX(-Y(L)*SCALE(IG)/26.,-40.)
      F(L)       = Y(L)*OHMI(L)*(1.-EXP(YEXP))
      IF(F(L) >  FMAX1(IG)) THEN
          FMAX1(IG) = F(L)
          YMAX1(IG) = Y(L)
C          FKOE      = 0.90*FMAX1(IG)
      ENDIF
      UREL       = U(L) - UBK(IGW)
      VREL       = V(L) - VBK(IGW)
      WREL       = W(L) - WBK(IGW)
      UTOT       = SQRT(UREL**2 + VREL**2 + WREL**2)
      IF(UTOT < UMIN1(IG)) UMIN1(IG) = UTOT
      IF(UTOT > UMAX1(IG)) UMAX1(IG) = UTOT
100   CONTINUE

      DO 9002 KG = KLOW+1,KUPP
C ... DIRECT THE COMPILER TO IGNORE DATA DEPENDENCIES (-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO 4000 IG = 1,JMAX*ISTRID
      L       = KG*IL + IA1 + IG
      IGW     = IG + JN*ISTRID
      Y(L)    = Y(L-IL) + .5*(D3(L)+D3(L-IL))
C ... CALCULATE F(Y) AND TOTAL VELOCITY
      YEXP    = MAX(-Y(L)*SCALE(IG)/26.,-40.)
      F(L)    = Y(L)*OHMI(L)*(1.-EXP(YEXP))
C ... FIND MAXIMUM VALUES OF UTOT AND F(Y)
      IF(F(L) >  FMAX1(IG)) THEN
          FMAX1(IG) = F(L)
          YMAX1(IG) = Y(L)
      ENDIF
      UREL       = U(L) - UBK(IGW)
      VREL       = V(L) - VBK(IGW)
      WREL       = W(L) - WBK(IGW)
      UTOT       = SQRT(UREL**2 + VREL**2 + WREL**2)
      IF(UTOT < UMIN1(IG)) UMIN1(IG) = UTOT
      IF(UTOT > UMAX1(IG)) UMAX1(IG) = UTOT
4000  CONTINUE
9002  CONTINUE

C ... PROFILE
      DO 5000 IG = 1,JMAX*ISTRID
      IA      = IA1 + IG
      IGW     = IG  + JN*ISTRID
      IF(KCP(IGW) == 1) UMIN1(IG) = 0.
5000  CONTINUE

C ... CALCULATE THE OUTER LAYER VISCOSITY
      Y0  = 0.
      DO 9000 KG = KLOW,KUPP
      DO 6000 IG = 1,JMAX*ISTRID
      L       = KG*IL + IA1 + IG
      IGW     = IG + JN*ISTRID
      XPOT    = MIN(.3*Y(L)/(YMAX1(IG)+EPS),30.)
      FKLEB   = 1./(1.+ 5.5*XPOT**6)
      CWK     = .25*(1.+ 60.*Y0)
      FWAKE   = CWK*YMAX1(IG)*(UMAX1(IG)-UMIN1(IG))**2/(FMAX1(IG)+EPS)
      FWAKE   = MIN(FWAKE,YMAX1(IG)*FMAX1(IG))
      RMUIL   = 0.16*RO(L)*F(L)**2/(OHMI(L)+EPS)
      RMUOL   = .0168*1.6*RO(L)*FWAKE*FKLEB
C
C ... CALCULATE THE MULTIPLIERS FOR TURBULENT FLOW
C
      RMU     = MIN(RMUIL,RMUOL)

      if(kcp(igw) /= 0) then
         RMYT    = RMU/VIS(L)
      else
      RMYT    = 0.
      endif
C ... MULTIPLIER IN MOMENTUM EQUATIONS IS ADDED TO PREVIOUS VALUES
c     EPS2(L) = 1.+ RMYT
      EPS2(L) = 1.+ SQRT(RMYT**2 + (EPS2(L)-1.)**2)
      EPS2(L) = MIN(TURLIM,EPS2(L))
C ... ENERGY EQUATION
c     EPS4(L) = 1.+ RMYT*PRS
      VIST(L) = SQRT((RMYT*VIS(L))**2 + VIST(L)**2)
      VIST(L) = MIN(TURLIM*VIS(L),VIST(L))
6000  CONTINUE
9000  CONTINUE

C ********************************************************************
C ... BOUNDARY NODE  REFLECTED SIDE .....                            *
C ********************************************************************
      IF (KTOP /= 0) THEN
C ... FIND MAXIMUM VALUES OF UTOT AND F(Y)
      DO 1200 IG = 1,JMAX*ISTRID
         FMAX1(IG) = 0.
         YMAX1(IG) = 0.
         UMAX1(IG) = 0.
         UMIN1(IG) = 1.E10
 1200 CONTINUE

      DO 1000 IG = 1,JMAX*ISTRID
      L          = IG + IA1 + IL*KTOP
      IGW        = IG + JN*ISTRID
      Y(L)       = .5*D3(L)
      SCALE(IG)  = SQRT(RO(L)*OHMI(L)/VIS(L))
      YEXP       = MAX(-Y(L)*SCALE(IG)/26.,-40.)
      F(L)       = Y(L)*OHMI(L)*(1.-EXP(YEXP))
      IF(F(L) >  FMAX1(IG)) THEN
          FMAX1(IG) = F(L)
          YMAX1(IG) = Y(L)
      ENDIF
      UREL       = U(L) - UTK(IGW)
      VREL       = V(L) - VTK(IGW)
      WREL       = W(L) - WTK(IGW)
      UTOT       = SQRT(UREL**2 + VREL**2 + WREL**2)
      IF(UTOT < UMIN1(IG)) UMIN1(IG) = UTOT
      IF(UTOT > UMAX1(IG)) UMAX1(IG) = UTOT
1000  CONTINUE

      DO 9001 KG = KTOP-1,KUPP+1,-1
C ... DIRECT THE COMPILER TO IGNORE DATA DEPENDENCIES (-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO 4050 IG = 1,JMAX*ISTRID
      L       = KG*IL + IA1 + IG
      IGW     = IG + JN*ISTRID
      Y(L)    = Y(L+IL) + .5*(D3(L)+D3(L+IL))
C ... CALCULATE F(Y) AND TOTAL VELOCITY
      YEXP    = MAX(-Y(L)*SCALE(IG)/26.,-40.)
      F(L)    = Y(L)*OHMI(L)*(1.-EXP(YEXP))
2050  CONTINUE
C ... FIND MAXIMUM VALUES OF UTOT AND F(Y)
      IF(F(L) >  FMAX1(IG)) THEN
          FMAX1(IG) = F(L)
          YMAX1(IG) = Y(L)
      ENDIF
      UREL       = U(L) - UTK(IGW)
      VREL       = V(L) - VTK(IGW)
      WREL       = W(L) - WTK(IGW)
      UTOT       = SQRT(UREL**2 + VREL**2 + WREL**2)
      IF(UTOT < UMIN1(IG)) UMIN1(IG) = UTOT
      IF(UTOT > UMAX1(IG)) UMAX1(IG) = UTOT
4050  CONTINUE
9001  CONTINUE

C ... PROFILE
      DO 9011 IG = 1,JMAX*ISTRID
      IA      = IA1 + IG
      IGW     = IG + JN*ISTRID
      IF(KCP(IGW) == 1) UMIN1(IG) = 0.
9011  CONTINUE

C ... CALCULATE THE OUTER LAYER VISCOSITY
      Y0  = 0.
      DO 9012 KG = KTOP,KUPP+1,-1
      DO 6050 IG = 1,JMAX*ISTRID
      L   = KG*IL + IA1 + IG
      IF(RKSI(L) <= 1.) THEN
      IGW     = IG + JN*ISTRID
      XPOT    = MIN(.3*Y(L)/(YMAX1(IG)+EPS),30.)
      FKLEB   = 1./(1.+ 5.5*XPOT**6)
      CWK     = .25*(1.+ 60.*Y0)
      FWAKE   = CWK*YMAX1(IG)*(UMAX1(IG)-UMIN1(IG))**2/(FMAX1(IG)+EPS)
      FWAKE   = MIN(FWAKE,YMAX1(IG)*FMAX1(IG))
      RMUIL   = 0.16*RO(L)*F(L)**2/(OHMI(L)+EPS)
      RMUOL   = .0168*1.6*RO(L)*FWAKE*FKLEB

C
C ... CALCULATE THE MULTIPLIERS FOR TURBULENT FLOW
C
      RMU     = MIN(RMUIL,RMUOL)
c      RMYT    = RMU/VIS(L)
      if(kcp(igw) /= 0) then
      RMYT    = RMU/VIS(L)
      else
      RMYT    = 0.
      endif
C ... MULTIPLIER IN MOMENTUM EQUATIONS IS ADDED TO PREVIOUS VALUES
c     EPS2(L) = 1.+ RMYT
      EPS2(L) = 1.+ SQRT(RMYT**2 + (EPS2(L)-1.)**2)
      EPS2(L) = MIN(TURLIM,EPS2(L))
C ... ENERGY EQUATION
c     EPS4(L) = 1.+ RMYT*PRS
      VIST(L) = SQRT((RMYT*VIS(L))**2 + VIST(L)**2)
      VIST(L) = MIN(TURLIM*VIS(L),VIST(L))
      ENDIF
6050  CONTINUE
9012  CONTINUE
C ... END OF TOP SURFACE
      ENDIF

      RETURN
      END SUBROUTINE TURBBL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C						
      SUBROUTINE TURBCS(EPS2,VIST,VIS,OHMI,U,V,W,RO,D3,F,Y,
     2 PR,PRT,IMAX,JMAX,KMAX,IN,JN,KN,KBOT,KTOP,KCP,M,TURLIM,
     3 UBK,VBK,WBK,UTK,VTK,WTK,nbl,ZZZ,MAW,MAXW,RKSI)

      DIMENSION :: U(*),RO(*),EPS2(*),VIST(*),V(*),VIS(*),W(*),
     3 D3(*),KCP(*),OHMI(*),Y(*),F(*),
     3 UBK(*),VBK(*),WBK(*),UTK(*),VTK(*),WTK(*),RKSI(*)

      INTEGER, ALLOCATABLE :: KMIN(:)  ! varsinainen viritys

      REAL, POINTER :: UMAX(:),UMIN(:),FMAX(:),SCALE(:),YMAX(:),
     +                 UDEL(:)
      REAL, TARGET  :: ZZZ(MAXW)

      UMAX => ZZZ( 0*MAW+1: 1*MAW);UMIN => ZZZ( 1*MAW+1: 2*MAW)
      FMAX => ZZZ( 2*MAW+1: 3*MAW);SCALE=> ZZZ( 3*MAW+1: 4*MAW)
      YMAX => ZZZ( 4*MAW+1: 5*MAW);UDEL => ZZZ( 5*MAW+1: 6*MAW)

      ALLOCATE(KMIN(MAW))

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      IA1     = (KN-1)*IL + JN*ISTRID
      IA2     = IA1 + 1
      IA3     = IA1 + JMAX*ISTRID
      IA4     = IL + IA1

      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      EPS = 1.E-10
      YKOE    = 1.
      PRS = PR/PRT
C
C ... CALCULATE TURBULENT VISCOSITIES BASED ON THE CEBECI-SMITH MODEL
C ... MODIFIED BY STOCK AND HAASE

C ... SET ZEROS FOR FMAX AND UMAX
      DO 100 IG   = 1,JMAX*ISTRID
         FMAX(IG) = 0.
         UMAX(IG) = 0.
         UMIN(IG) = 1.E10
 100  CONTINUE

C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...
C ... CALCULATE F(Y) AND HIGHT OF THE LAYER Y(X,Y,Z)
C ... FIND FMAX, YMAX, UMAX AND UMIN
      DO 110 IG = 1,JMAX*ISTRID
         L        = IA4 + IG
         Y(L)     = .5*D3(L)
         IGW      = IG + JN*ISTRID
         SCALE(IG)= SQRT(RO(L)*OHMI(L)/VIS(L))
         F(L)     = Y(L)*OHMI(L)*(1.-EXP(-Y(L)*SCALE(IG)/26.))
         UTOT = SQRT(U(L)**2 + V(L)**2 + W(L)**2)
         IF(F(L) > FMAX(IG)) THEN
            FMAX(IG)  = F(L)
            YMAX(IG)  = Y(L)
C           FKOE      =.90*FMAX
         ENDIF
         IF(UTOT < UMIN(IG)) THEN
            UMIN(IG) = UTOT
            KMIN(IG) = 1
         ENDIF
         IF(UTOT > UMAX(IG)) UMAX(IG) = UTOT
 110  CONTINUE

      DO 150 KG = 2,KUPP
C ... DIRECT THE COMPILER TO IGNORE DATA DEPENDENCIES (-pfa)
C*$*ASSERTDO(CONCURRENT)
         DO 150 IG = 1,JMAX*ISTRID
            L      = KG*IL + IA1 + IG
            IGW    = IG + JN*ISTRID
            Y(L)   = Y(L-IL) + .5*(D3(L)+D3(L-IL))
            F(L)   = Y(L)*OHMI(L)*(1.-EXP(-Y(L)*SCALE(IG)/26.))
            UTOT = SQRT(U(L)**2 + V(L)**2 + W(L)**2)
            IF(F(L) > FMAX(IG)) THEN
               FMAX(IG)  = F(L)
               YMAX(IG)  = Y(L)
C              FKOE      =.90*FMAX
            ENDIF

            IF(UTOT < UMIN(IG)) THEN
               UMIN(IG) = UTOT
               KMIN(IG) = KG
            ENDIF
            IF(UTOT > UMAX(IG)) UMAX(IG) = UTOT
 150  CONTINUE


C ... PROFILE
      DO 200 IG = 1,JMAX*ISTRID
         IGW    = IG + JN*ISTRID
         IF(KCP(IGW) == 1) THEN
            UMIN(IG) = 0.
            KMIN(IG) = 1
         ENDIF
 200  CONTINUE

C ... CALCULATE BOUNDARY LAYER THICKNESS AND DISPLACEMENT THICKNESS
      DO 400 IG = 1,JMAX*ISTRID
         L       = IA4 + IG
         IGW     = IG + JN*ISTRID
         DELTA   = 1.936*YMAX(IG)
         KDEL    = 1
         UDEL(IG)= SQRT(U(L)**2 + V(L)**2 + W(L)**2)
         DO 300 KG = 1,KUPP
            L          = KG*IL + IA1 + IG
            IF(Y(L) < DELTA) THEN
               KDEL    = KG
               UDEL(IG)= SQRT(U(L)**2 + V(L)**2 + W(L)**2)
            ENDIF
 300     CONTINUE
         DELTAX  = 0.
         PUDEL   = 1./(UDEL(IG) + EPS)
         DO 350 KG = KMIN(IG),KDEL
            L      = KG*IL + IA1 + IG
            DELTAX = DELTAX + MAX(1.-SQRT(U(L)**2+V(L)**2+
     +           W(L)**2)*PUDEL,0.)*D3(L)
 350     CONTINUE
         UDEL(IG)  = DELTAX*UDEL(IG)
 400  CONTINUE

C ... CALCULATE INNER AND OUTER LAYER VISCOSITY
      DO 600 KG = 1,KUPP
         DO 600 IG = 1,JMAX*ISTRID
            L      = KG*IL + IA1+IG
            IF(RKSI(L) <= 1.) THEN
            IGW    = IG + JN*ISTRID
            DELTA  = 1.936*YMAX(IG)
            XPOT   = MIN(Y(L)/(DELTA+EPS),30.)
            FKLEB  = 1./(1.+ 5.5*XPOT**6)
C            FKLEB  = 1./(1.+ 5.5*(Y(L)/(DELTA+EPS))**6)
C           FWAKE  = .25*YMAX(IG)*(UMAX(IG)-UMIN(IG))**2/(FMAX(IG)+EPS
C           FWAKE  = MIN(FWAKE,YMAX(IG)*FMAX(IG))
            RMUO   = .0168*RO(L)*UDEL(IG)*FKLEB
            RMUI   = 0.16*RO(L)*F(L)**2/(OHMI(L)+EPS)
C ... CALCULATE THE MULTIPLIERS FOR TURBULENT FLOW
            RMU    = MIN(RMUI,RMUO)
            RMYT   = RMU/VIS(L)
C ... MULTIPLIER IN MOMENTUM EQUATIONS
            EPS2(L)= 1.+ RMYT
C ... ENERGY EQUATION
            VIST(L)= RMU
            ENDIF
 600  CONTINUE
      DEALLOCATE(KMIN)

      RETURN
      END SUBROUTINE TURBCS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE TURBSC(EPS2,VIST,VIS,OHMI,U,V,W,RO,D3,
     2 PR,PRT,IMAX,JMAX,KMAX,IN,JN,KN,KBOT,KTOP,KCP,M,TURLIM,
     3 ZZZ,MAW,MAXW,RKSI)

      DIMENSION :: U(*),RO(*),EPS2(*),VIST(*),V(*),VIS(*),W(*),
     3             D3(*),KCP(*),OHMI(*),RKSI(*)

      REAL, POINTER :: Y(:),F(:),RMUI(:),RMUO(:),UTOT(:),UTAN(:)
      REAL, TARGET  :: ZZZ(MAXW)

      Y    => ZZZ( 0*MAW+1: 1*MAW);F    => ZZZ( 1*MAW+1: 2*MAW)
      RMUI => ZZZ( 2*MAW+1: 3*MAW);RMUO => ZZZ( 3*MAW+1: 4*MAW)
      UTOT => ZZZ( 4*MAW+1: 5*MAW);UTAN => ZZZ( 5*MAW+1: 6*MAW)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KMAX/2
      IF(KBOT == 0) KUPP = 0
      FMAX    = 0.
      EPS     = 1.E-10
C
C ... CALCULATE TURBULENT VISCOSITIES BASED ON SCHETZ (=SC) MODEL
C ... MODEL IS SUITABLE ONLY FOR JETS
C

C ... CALCULATE THE INITIAL RADIUS OF THE JET (PULTATTU, I=3, K=16)
      RJET    = 0.
      IA      = (KN-1)*IL + JN*ISTRID + 3
      DO 4500 KG = 1,8
         L       = KG*IL + IA
         RJET    = RJET + D3(L)
 4500 CONTINUE
      RJET   = 2.6395E-3
      DO 9000 IG = 1,JMAX*ISTRID
         IA      = (KN-1)*IL + JN*ISTRID + IG
         KINF    = KUPP*IL + IA
         DELTAX2 = 0.
         YD      = -.5*D3(IL+IA)
         UDEL    = RO(KINF)*SQRT(U(KINF)**2 + V(KINF)**2 + W(KINF)**2)
         
C ... CALCULATE THE DISPLACEMENT THICKNESS
         DO 4600 KG = 1,KUPP
            L       = KG*IL + IA
            YD      = YD + D3(L)
            UTO     = SQRT(U(L)**2 + V(L)**2 + W(L)**2)
            DELTAX2 = DELTAX2 + MAX(RO(L)*UTO-UDEL,0.)*6.2832*YD*D3(L)
 4600    CONTINUE
         DELTAX2 = DELTAX2/(3.14159*(UDEL + EPS))

C ... CALCULATE THE OUTER LAYER VISCOSITY
         DO 5000 KG = 1,KUPP
            L       = KG*IL + IA
            FKLEB   = 1.
            RMUO(KG)= .018*UDEL*DELTAX2/RJET*FKLEB *.5
 5000    CONTINUE

C
C ... CALCULATE THE MULTIPLIERS FOR TURBULENT FLOW
C
         PRS     = PR/PRT
         DO 6000 KG = 1,KUPP
            L   = KG*IL + IA
            IF(RKSI(L) <= 1.) THEN
            RMU     = RMUO(KG)
            RMYT    = RMU/VIS(L)
            EPS2(L) = 1.+ RMYT
C ... ENERGY EQUATION
            VIST(L) = RMU
            ENDIF
 6000    CONTINUE
         
 9000 CONTINUE

      RETURN
      END SUBROUTINE TURBSC
C
C ----------------------------------------------------------------------
C --- Two-equation Models ----------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TURBCO(ITURB,JRDIS,JRPRE,C1,C2,C3,C21,CMU,CTA,PSIGK,
     +     PSIGE,AA1,ETA0)

      USE CONSTANTS, ONLY : CB1, CB2, SRNU, CV1, CW1, CW2, CW3

      AA1  = 0.
      C1   = 0.
      C2   = 0.
      C21  = 0.
      CC   = 0.
      CMU  = 0.
      CTA  = 0.
      ETA0 = 0.
      PSIGE = 0.
      PSIGK = 0.
      
C ... SPECIFY THE COEFFICIENTS APPLIED IN 2-EQUATION MODELS

      IF(ITURB == 3 .OR. ITURB >= 21 .OR. ITURB == 11) THEN
C ... MODIFIED CHIEN. IT IS ALSO USED WITH REYNOLDS STRESS MODEL IN
C ... ENERGY EQUATION.
         C1      = 1.44
         C2      = 1.92
         PSIGK   = 1./1.
         PSIGE   = 1./1.3
         C21     = 0.22
         CMU     = 0.09
         IF(JRDIS == 2) THEN
            C1 = 1.39
            C2 = 1.83
         ENDIF
      ENDIF

      IF(ITURB == 4) THEN
C ... MODIFIED CHIEN + MENTER'S LIMITATION OF TURBULENT VISCOSITY
         C1      = 1.44
         C2      = 1.92
         PSIGK   = 1./1.
         PSIGE   = 1./1.3
         C21     = 0.22
         CMU     = 0.09
         AA1     = 0.31
      ENDIF

      IF(ITURB == 5) THEN
C ... RNG K-EPSILON
         C1      = 1.42
         C2      = 1.68
         PSIGK   = 1./1.39
         PSIGE   = 1./1.39
         C21     = 0.22
         CMU     = 0.085
         AA1     = 0.31
         ETA0    = SQRT((C2 - 1.)/(CMU*(C1 - 1.)))
      ENDIF

      IF(ITURB == 6) THEN 
C ... CROSS-DIFFUSION
         XKAPPA  = 0.4
         C2      = 1.92
c         C1      = 1.17!=C2 - XKAPPA**2/SQRT(CMU)*PSIGE XKAPPA = 0.41
C     ei limitysta lahdetermissa. silloin nama hyvat
C     lahtee oletuksesta psige = 1./.8
         C1      = 1.25!=C2 - XKAPPA**2/SQRT(CMU)*PSIGE XKAPPA = 0.4
         C3      = -1.13  
         PSIGE   = 1./.8
C     lahtee oletuksesta psige = 1./.75
         C1      = 1.173!=C2 - XKAPPA**2/SQRT(CMU)*PSIGE XKAPPA = 0.41
         C3      = -1.26  
         PSIGE   = 1./.75


         PSIGK   = 1./1.
         C21     = 0.
         CMU     = 0.09
         AA1     = 0.
      ENDIF

      IF(ITURB == 7) THEN
C ...  K-OMEGA TESTS USING EPSILON-EQUATION
         C1      = 1.44
         C2      = 1.92
         PSIGK   = 1./.5
         PSIGE   = 1./.5
         C21     = 0.22
         CMU     = 0.09
         AA1     = 0.31
      ENDIF

      IF(ITURB == 8) THEN
C ... Smagorinsky
         CMU     = 0.2
      ENDIF

      IF(ITURB == 9) THEN
C ... Spalart-Allmaras
         CB1     = 0.1335 ; CB2 = 0.622
         CV1     = 7.1
         SRNU    = 2./3.   ! Kokeili *.001
         CKAPPA  = 0.41
         CW1     = CB1/CKAPPA**2 + (1.+CB2)/SRNU ; CW2 = 0.3 ; CW3 = 2.0
      ENDIF

      IF(ITURB == 10) THEN
C ... Speziale et al. 1995. Explicit ASM model
         C2      = 1.83
         PSIGK   = 1./1.
         PSIGE   = 1./1.3
         C21     = 0.22
         CMU     = 0.088
         CTA     = 12.5
         C1      = 1.39 !=C2 - KAPPA**2/SQRT(CMU)*PSIGE
      ENDIF

      RETURN
      END SUBROUTINE TURBCO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLVD(VIS,OHMI,RO,D3,IMAX,JMAX,KMAX,IN,JN,KN,
     +     KBOT,KTOP,KCP,LAMIN,F,Y,F2,ITURB,
     +     KX1,KX2,KY1,KY2,IDIR2,SCALE)

      REAL    :: RO(*),VIS(*),D3(*),OHMI(*),F(*),Y(*),F2(*)
      INTEGER :: KCP(*)
      REAL    :: SCALE(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      IF(LAMIN == 2) RETURN

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      IB      = (KN-1)*IL + JN*ISTRID

      RC1     = 1.0
      RC2     = 0.0

C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF
C *************************************************************

C ... CALCULATE THE VAN DRIEST WALL DISSIPATION TERMS

C *********************************************************************
C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...

      KA        = (KN+KBOT2-1)*IL
      DO 800 J  = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO 800 I  = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C     -------- THE EFFECT OF THE WALLS ----------------------------------
         OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L+IDIR2*IL))
               SCALE(IG)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY BASED y+
               Y(L)     = .5*D3(L) ! Could be replaced by wall distances
               YEXP     = MAX(-Y(L)*SCALE(IG)/A,-30.)
               F(L)     = (1.- EXP(YEXP))*F(L)
c              F2(L)    = MAX(1.,F2(L)) ! Unused
 800     CONTINUE

      DO 2000 KG = KZ1,KZ2,IDIR2
      KA        = (KN+KG-1)*IL
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO 2000 J = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO 2000 I = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C  -------- THE EFFECT OF THE WALLS ----------------------------------
C ...       CALCULATE F(Y)
               Y(L)     = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
               YEXP     = MAX(-Y(L)*SCALE(IG)/A,-30.)
               F(L)     = (1.- EXP(YEXP))*F(L)
C              F2(L)    = MAX(1.,F2(L)) ! Unused
 2000    CONTINUE

C ************** WALL CORRECTION CALCULATED ***************************

      RETURN
      END SUBROUTINE WALLVD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLKE(VIS,OHMI,RO,D3,IMAX,JMAX,KMAX,IN,JN,KN,
     +     KBOT,KTOP,KCP,LAMIN,F,Y,F2,ITURB,
     +     KX1,KX2,KY1,KY2,IDIR2,SCALE)

      REAL    :: RO(*),VIS(*),D3(*),OHMI(*),F(*),Y(*),F2(*)
      INTEGER :: KCP(*)
      REAL    :: SCALE(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      IF(LAMIN == 2) RETURN

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      IB      = (KN-1)*IL + JN*ISTRID

      RC1     = 1.0
      RC2     = 0.0

C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF
C *************************************************************


C
C ... CALCULATE THE WALL DISSIPATION TERMS IN THE K-EPSILON MODEL
C

C ... INITIAL VALUES FOR WALL FUNCTION COULD BE GIVEN HERE
C     IF(LAMIN == 1) THEN
C     DO 750 IG =1,IL*(KMAX+KN*2)
C       F(IG)  = 1.
C750  CONTINUE
C     ENDIF
C *********************************************************************
C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...
      IF (ITURB /= 10) THEN ! Speziale's model does not need wall functions
      KA        = (KN+KBOT2-1)*IL
      DO 800 J  = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO 800 I  = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C     -------- THE EFFECT OF THE WALLS ----------------------------------
         OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L+IDIR2*IL))
               SCALE(IG)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY BASED
C              APU      = SQRT(REPS(L)/(VIS(L)))   ! KOLMOGOROV'S
C              SCALE(IG)= SQRT(RO(L)/VIS(L)*APU)   ! SCALE
               Y(L)     = .5*D3(L)
               YEXP     = MAX(-Y(L)*SCALE(IG)*.0115,-30.)
               F(L)     = (1.- EXP(YEXP))*F(L)
               F2(L)    = MAX(1.,F2(L))
 800     CONTINUE

      DO 2000 KG = KZ1,KZ2,IDIR2
      KA        = (KN+KG-1)*IL
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO 2000 J = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO 2000 I = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C  -------- THE EFFECT OF THE WALLS ----------------------------------
C ...       CALCULATE F(Y)
               Y(L)     = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
               YEXP     = MAX(-Y(L)*SCALE(IG)*0.0115,-30.)
               F(L)     = (1.- EXP(YEXP))*F(L)
               F2(L)    = MAX(1.,F2(L))
 2000    CONTINUE
      ENDIF ! IF (ITURB /= 10)
C ************** WALL CORRECTION CALCULATED ***************************

      RETURN
      END SUBROUTINE WALLKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EPS2KE(EPS2,VIST,VIS,S,OHMI,RO,RK,REPS,C,IMAX,JMAX,
     + KMAX,IN,JN,KN,M,F,F2,DDEPS,ITURB,TURLIM,RKLIM,EPSLIM,CMU,AA1,
     + RKSI,VOL)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,M,ISTRID,JSTRID,IL,ITURB,
     +           I,J,K,KA,KK,JJ,L

      REAL :: EPS,EPSLIM,RMUFRS,RKLIM,XMA,EPD,RMU,CMU,RMULO,
     +        TURBLE,TURLIM,RMU1,AA1,RMU0,RMYT,C1,C2,C3,C4,C5,GG,
     +        ALF1,ALF2,ALF3,SQR2,EPSI,RKPE,ETA,XSI,CMUT

      REAL :: RO(*),EPS2(*),VIST(*),VIS(*),OHMI(*),RK(*),REPS(*),
     +        F(*),F2(*),DDEPS(*),C(*),S(*),RKSI(*),VOL(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      EPS     = EPSLIM/196.047                    ! ANTAKEE MERSU
      RMUFRS  = 0.09*RKLIM**2/EPSLIM

C
C ... CALCULATE TURBULENT VISCOSITIES BASED ON THE K-EPSILON MODEL
C
      IF(M == 1) THEN

C ************** VISCOSITY ON THE FIRST GRID LEVEL ********************

       IF(ITURB == 3 .OR. ITURB >= 20 .OR. ITURB == 11 .OR. 
     +    ITURB == 5 .OR. ITURB == 7) THEN !CHIEN'S MODEL
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7000 K = 1,KMAX
        KK      = K*IL + KA
        DO 7000 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7000 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU     = CMU*F(L)*RK(L)**2/(REPS(L)+EPS)
            RMYT    = MAX(0.,RMU/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*REPS(L)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE)*MIN(1.*F(L),1.)
C ...                                     kopin ekaan arvoon johtoreunalla
C ...                                     saattaa tulla ei fys ilman tata
C ...                                     PPR 18.11.98
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7000  CONTINUE
      ENDIF

      IF(ITURB == 4) THEN   ! CHIEN WITH MENTER'S LIMITATION
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7100 K = 1,KMAX
        KK      = K*IL + KA
        DO 7100 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7100 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            RMU0    = CMU*F(L)*RK(L)**2/(REPS(L)+EPS)
            RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*REPS(L)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
c      if(k == 1 .and. j == 14 .and. i == 1)then
c      write(667,*) 'RMU= ',RMU
c      write(667,*) 'RMYT= ',RMYT
c      write(667,*) 'RMULO,RMUFRS= ',RMULO,RMUFRS
c      write(667,*) 'TURBLE= ',TURBLE
c      write(667,*) 'F(L)= ',F(L)
c      write(667,*) 'EPS2= ',EPS2(L)
c      write(667,*) 'EPS2OR= ',1.+RMYT
c      write(667,*) 'reps,epslim= ',REPS(L),EPSLIM
c      endif
7100  CONTINUE
      ENDIF

      IF(ITURB == 55) THEN   ! RNG K-EPSILON. CURRENT VERSION DOESN'T WORK
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7200 K = 1,KMAX
        KK      = K*IL + KA
        DO 7200 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7200 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            RMU0    = CMU*RK(L)**2/(REPS(L)+DDEPS(L)+EPS)
            RMU0    = RMU0*(1.+ 2.*SQRT(VIS(L)/(RMU0+EPS)))
C           RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*REPS(L)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7200  CONTINUE
      ENDIF

      IF(ITURB == 7) THEN    ! K-OMEGA TESTS USING EPSILON-EQUATION
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO K = 1,KMAX
        KK      = K*IL + KA
        DO J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            RMU0    = CMU*F(L)*RK(L)**2/(REPS(L)+EPS)
            RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*REPS(L)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
        ENDDO; ENDDO; ENDDO ! ANtakee mersu
      ENDIF

      IF(ITURB == 8) THEN    ! Smagorinsky (zero-equation model)
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7300 K = 1,KMAX
        KK      = K*IL + KA
        DO 7300 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7300 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU0    = RO(L)*F(L)*(CMU*VOL(L)**.33)**2*S(L)
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
            EPS2(L) = MIN(TURLIM,EPS2(L))
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            VIST(L) = MAX(0.,RMU0)
            ENDIF
7300  CONTINUE
      ENDIF
C **********************************************************************
      ENDIF

       IF(ITURB == 10) THEN  ! Speziale
C ... LRR mallin vakiot

      C1 = 3.
      C2 = 0.8
      C3 = 1.75
      C4 = 1.31
      C5 = 2.

C ... SSG mallin vakiot

      C1 = 6.8
      C2 = 0.36
      C3 = 1.25
      C4 = 0.4
      C5 = 1.88

      GG      = 1./(C1/2. + C5 - 1.)
      ALF1    = (4./3.-C2)*GG/2.
      ALF2    = (2. - C3)**2*GG**2/4.
      ALF3    = (2. - C4)**2*GG**2/4.
      SQR2    = 1./SQRT(2.)
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7400 K = 1,KMAX
        KK      = K*IL + KA
        DO 7400 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7400 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
         EPSI    = REPS(L) + DDEPS(L)
         RKPE    = RK(L)/(EPSI+EPS)
         ETA     = ALF2*(S(L)*RKPE)**2
         XSI     = ALF3*(SQR2*OHMI(L)*RKPE)**2
         CMUT    = 3.*(1. + ETA)/(3. + ETA + 6.*XSI*ETA + 6.*XSI)*ALF1
            RMU     = CMUT*RK(L)*RKPE
            RMYT    = MAX(0.,RMU/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*REPS(L)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
            VIST(L) = MAX(0.,RMU)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7400  CONTINUE
      ENDIF
 8000 CONTINUE
      RETURN
      END SUBROUTINE EPS2KE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLF(FWC,FWX,FWY,FWZ,A3XA,A3YA,A3ZA,D3,IMAX,JMAX,KMAX,
     +     IN,JN,KN,KBOT,KTOP,LAMIN,Y,ITURB,
     +     KX1,KX2,KY1,KY2,IDIR2,IWFORCE)

      USE CONSTANTS, ONLY : EPS10

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IDIR2,IN,JN,KN,KBOT,KTOP,LAMIN,ITURB,
     +           KX1,KX2,KY1,KY2,ISTRID,JSTRID,IL,KLOW,KUPP,
     +           IB,KBOT2,KZ1,KZ2,KA,J,IJ,I,N1,IG,L,IGW,KG,IWFORCE

      REAL :: RC1,RC2,CW,RKF,DBUBBLE,UR
      REAL :: FWC(*),FWX(*),FWY(*),FWZ(*),A3XA(*),A3YA(*),A3ZA(*),
     +        D3(*),Y(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      IF(LAMIN == 2) RETURN

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      IB      = (KN-1)*IL + JN*ISTRID

      RC1     = 1.0
      RC2     = 0.0

C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF

C *********************************************************************
C ... Wall force (1). Hosokawa et al. (2002) towards all walls
C *********************************************************************

C .... Artificial parameters

c      DBUBBLE = 1.E-3

C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...

      KA        = (KN+KBOT2-1)*IL
      DO J  = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO I  = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C     -------- THE EFFECT OF THE WALLS ----------------------------------
         Y(L)    = .5*D3(L)
         RKF     = FWC(L)/Y(L)**2
         FWX(L)  = FWX(L) + IDIR2*RKF*A3XA(L)
         FWY(L)  = FWY(L) + IDIR2*RKF*A3YA(L)
         FWZ(L)  = FWZ(L) + IDIR2*RKF*A3ZA(L)
c         write(666,*) i,j,kbot,fwc(l),fwx(l),fwy(l),fwz(l)
      ENDDO; ENDDO

      DO KG = KZ1,KZ2,IDIR2
      KA        = (KN+KG-1)*IL
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO J = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO I = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C  -------- THE EFFECT OF THE WALLS ----------------------------------
         Y(L)     = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
         RKF     = FWC(L)/Y(L)**2
         FWX(L)  = FWX(L) + IDIR2*RKF*A3XA(L)
         FWY(L)  = FWY(L) + IDIR2*RKF*A3YA(L)
         FWZ(L)  = FWZ(L) + IDIR2*RKF*A3ZA(L)
c         write(666,*) i,j,kg,fwc(l),fwx(l),fwy(l),fwz(l)

      ENDDO; ENDDO; ENDDO

      RETURN
      END SUBROUTINE WALLF
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLMF(VIS,OHMI,RO,D3,IMAX,JMAX,KMAX,IN,JN,KN,
     +     KBOT,KTOP,LAMIN,F,Y,ITURB,
     +     KX1,KX2,KY1,KY2,IDIR2,SCALE)

      USE CONSTANTS, ONLY : EPS10

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,KBOT,KTOP,LAMIN,ITURB,
     +           KX1,KX2,KY1,KY2,IDIR2,ISTRID,JSTRID,IL,KLOW,KUPP,
     +           IB,KBOT2,KZ1,KZ2,KA,J,IJ,I,N1,IG,L,IGW,KG

      REAL :: RC1,RC2,CS,ALFAG,DBUBBLE,UR,OHMIB,YEXP
      REAL :: RO(*),VIS(*),D3(*),OHMI(*),F(*),Y(*)
      REAL :: SCALE(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      IF(LAMIN == 2) RETURN

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      IB      = (KN-1)*IL + JN*ISTRID

      RC1     = 1.0
      RC2     = 0.0

C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF

C *********************************************************************
C ... Bubble induced turbulence. Sato (1981) simplified
C *********************************************************************

C .... Artificial parameters

      CS      = 0.6 ! Sato parameter
      ALFAG   = 1.
      DBUBBLE = 1.E-3
      UR      = 0.2

C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...

      KA        = (KN+KBOT2-1)*IL
      DO J  = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO I  = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C     -------- THE EFFECT OF THE WALLS ----------------------------------
            OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L+IDIR2*IL))
            SCALE(IG)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY BASED
            Y(L)     = .5*D3(L)
            YEXP     = MAX(-Y(L)*SCALE(IG)/16.,-30.)
            IF(F(L) <= 0.) F(L) = 1.
            F(L)     = (1.- EXP(YEXP))**2*F(L)
      ENDDO; ENDDO

      DO KG = KZ1,KZ2,IDIR2
      KA        = (KN+KG-1)*IL
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO J = KY1,KY2
      IJ        = (JN+J-1)*ISTRID + KA
      N1        = (JN+J-1)*ISTRID + IN 
      DO I = KX1,KX2
         IG     = I + (J-1)*ISTRID
         L      = (IN+I-1)*1 + IJ + 1
         IGW    = I + N1
C  -------- THE EFFECT OF THE WALLS ----------------------------------
               Y(L)     = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
               YEXP     = MAX(-Y(L)*SCALE(IG)/16.,-30.)
               IF(F(L) <= 0.) F(L) = 1.
               F(L)     = (1.- EXP(YEXP))**2*F(L)
      ENDDO; ENDDO; ENDDO

      RETURN
      END SUBROUTINE WALLMF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EPS2CO(EPS2,VIST,VIS,S,OHMI,RO,RK,REPS,C,IMAX,JMAX,
     + KMAX,IN,JN,KN,M,F,F2,DDEPS,ITURB,TURLIM,RKLIM,EPSLIM,CMU,AA1,
     + RKSI)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,M,ISTRID,JSTRID,IL,ITURB,
     +           I,J,K,KA,KK,JJ,L

      REAL :: EPS,EPSLIM,RMUFRS,RKLIM,XMA,EPD,RMU,CMU,RMULO,
     +        TURBLE,TURLIM,RMU1,AA1,RMU0,RMYT

      REAL :: RO(*),EPS2(*),VIST(*),VIS(*),OHMI(*),RK(*),REPS(*),
     +        F(*),F2(*),DDEPS(*),C(*),S(*),RKSI(*)

C ... SAME AS EPS2KE. ONLY COMPRESSIBLE EFFECTS HAVE BEEN INCLUDED

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      EPS     = EPSLIM/196.047                    ! ANTAKEE MERSU
      RMUFRS  = 0.09*RKLIM**2/EPSLIM

C
C ... CALCULATE TURBULENT VISCOSITIES BASED ON THE K-EPSILON MODEL
C
      IF(M == 1) THEN

C ************** VISCOSITY ON THE FIRST GRID LEVEL ********************

       IF(ITURB == 3 .OR. ITURB >= 20 .OR. ITURB == 11 .OR. 
     +    ITURB == 5 .OR. ITURB == 7) THEN !CHIEN'S MODEL
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7000 K = 1,KMAX
        KK      = K*IL + KA
        DO 7000 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7000 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            RMU     = CMU*F(L)*RK(L)**2/(REPS(L)+EPS+EPD)
            RMYT    = MAX(0.,RMU/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*(REPS(L)+EPD)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7000  CONTINUE
      ENDIF

      IF(ITURB == 4) THEN   ! CHIEN WITH MENTER'S LIMITATION
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7100 K = 1,KMAX
        KK      = K*IL + KA
        DO 7100 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7100 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            RMU0    = CMU*F(L)*RK(L)**2/(REPS(L)+EPD+EPS)
            RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*(REPS(L)+EPD)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7100  CONTINUE
      ENDIF

      IF(ITURB == 55) THEN   ! RNG K-EPSILON. CURRENT VERSION DOESN'T WORK
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7200 K = 1,KMAX
        KK      = K*IL + KA
        DO 7200 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7200 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            RMU0    = CMU*RK(L)**2/(REPS(L)+DDEPS(L)+EPD+EPS)
            RMU0    = RMU0*(1.+ 2.*SQRT(VIS(L)/(RMU0+EPS)))
C           RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*(REPS(L)+EPD)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7200  CONTINUE
      ENDIF

      IF(ITURB == 7) THEN    ! K-OMEGA TESTS USING EPSILON-EQUATION
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7300 K = 1,KMAX
        KK      = K*IL + KA
        DO 7300 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7300 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            RMU1    = AA1*RK(L)/(OHMI(L)+EPS)
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            RMU0    = CMU*F(L)*RK(L)**2/(REPS(L)+EPD+EPS)
            RMU0    = RMU0/MAX(1.,RMU0/RMU1*F2(L))
            RMYT    = MAX(0.,RMU0/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*(REPS(L)+EPD)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
C           EPS2(L) = MIN(TURLIM,EPS2(L))
            VIST(L) = MAX(0.,RMU0)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7300  CONTINUE
      ENDIF
C **********************************************************************
      ENDIF

       IF(ITURB == 10) THEN  ! Speziale
        KA      = (KN-1)*IL + (JN-1)*ISTRID + IN
        DO 7400 K = 1,KMAX
        KK      = K*IL + KA
        DO 7400 J = 1,JMAX
        JJ      = J*ISTRID + KK
        DO 7400 I = 1,IMAX
            L       = JJ + I
            IF(RKSI(L) <= 1.) THEN
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            RMU     = CMU*RK(L)**2/(REPS(L)+EPS+EPD)
            RMYT    = MAX(0.,RMU/VIS(L))
            EPS2(L) = 1.+ RMYT
C ... TURLIM IS MODIFIED ON THE EDGE OF THE BOUNDARY LAYER (LOW REPS)
            RMULO   = 10.*RMUFRS/VIS(L)
            TURBLE  = TURLIM*(1.-1.1/(1.+.1*(REPS(L)+EPD)/EPSLIM))
            TURBLE  = MAX(RMULO,TURBLE) 
            EPS2(L) = MIN(TURBLE+1.,EPS2(L)) 
            VIST(L) = MAX(0.,RMU)
            VIST(L) = MIN(TURBLE*VIS(L),VIST(L))
            ENDIF
7400  CONTINUE
      ENDIF
 8000 CONTINUE
      RETURN
      END SUBROUTINE EPS2CO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TURBKE(EPS2,VIST,VIS,OHMI,U,V,W,RO,D3,RK,REPS,DDEPS,
     +     SRK,SEPS,PTUR,VOL,A3,A3X,A3Y,A3Z,PR,PRT,IMAX,JMAX,KMAX,
     +     IN,JN,KN,KBOT,KTOP,KCP,LAMIN,NBL,M,EPS5,F,Y,ITURB,
     +     TURLIM,RKLIM,EPSLIM,CMU,C1,ISTR,JSTR,KSTR,
     +     KX1,KX2,KY1,KY2,IDIR2,SCALE)

      REAL :: U(*),RO(*),EPS2(*),VIST(*),V(*),VIS(*),W(*),
     +        OHMI(*),VOL(*),A3(*),A3X(*),A3Y(*),A3Z(*),D3(*),DDEPS(*),
     +        RK(*),REPS(*),SRK(*),SEPS(*),PTUR(*),EPS5(*),F(*),Y(*)

      INTEGER :: KCP(*)

      REAL SCALE(*)
C
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = KSTR

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0
      IB      = (KN-1)*IL + JN*ISTRID
      IB0     = (KN-1)*IL +  0*ISTRID

      EPS     = EPSLIM/196.047                    ! ANTAKEE MERSU
      RC1     = 1.0
      RC2     = 0.0


C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF
C *************************************************************

C
C ... CALCULATE WALL SOURCE FOR K-EPSILON.
C ... CALCULATION OF TURBULENT VISCOSITIES CAN BE ACTIVATED
C ... CHIEN'S MODEL IS APPLIED
C
C ... SET INITIAL VALUES FOR WALL FUNCTION.(IF DIRECTIONAL TURBULENCE)
C     IF(LAMIN == 1) THEN
C     DO 750 IG = 1,IL*(KMAX + 2*KN)
C        F(IG)  = 1.
C750  CONTINUE
C     ENDIF
C *********************************************************************
C ... BOUNDARY NODE     NORMAL SIDE, INCREASING OR DECREASING (IDIR2)K ...
      KA        = (KN+KBOT2-1)*KSTR
      DO 800 J  = KY1,KY2
         IJ        = (JN+J-1)*JSTR + KA
         N1        = (JN+J-1)*ISTRID + IN 
      DO 800 I  = KX1,KX2
            IG     = I + (J-1)*ISTRID
            L      = (IN+I-1)*ISTR + IJ + 1
            IGW    = I + N1
            Y(L)   = .5*D3(L)
          OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L+IDIR2*IL))
          SCALE(IG)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY
C              RETUR    = RK(L)**2/(VIS(L)*(REPS(L)+1.E-16))
C              YEXP     = MAX(-3.4/(1.+.02*RETUR)**2,-30.)
C              F(L)     = EXP(YEXP)
C              YEXP     = MAX(-Y(L)*SCALE(IG)*.0115,-30.)
C              F(L)     = (1.- EXP(YEXP))

               YEXP2    = MAX(-Y(L)*SCALE(IG)*0.5,-30.)
               YPRO     = 1./RO(L)
               PYN2     = 1./Y(L)**2
               VAPU     = 2.*VIS(L)*YPRO*PYN2
               SRK(L)   = SRK(L)  - VAPU*RK(L)
               DDEPS(L) = DDEPS(L)+ VAPU*RK(L)
               SEPS(L)  = SEPS(L) - VAPU*REPS(L)*EXP(YEXP2)
 800     CONTINUE


      DO 2000 KG = KZ1,KZ2,IDIR2
         KA        = (KN+KG-1)*KSTR
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
      DO 2000 J = KY1,KY2
         IJ        = (JN+J-1)*JSTR + KA
         N1        = (JN+J-1)*ISTRID + IN 
      DO 2000 I = KX1,KX2
            IG     = I + (J-1)*ISTRID
            L      = (IN+I-1)*ISTR + IJ + 1
            IGW    = I + N1
C     ----- THE EFFECT OF THE WALLS -----------------------------------
                Y(L)    = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
                YEXP2   = MAX(-Y(L)*SCALE(IG)*0.5,-30.)
                YPRO    = 1./RO(L)
                PYN2    = 1./Y(L)**2
                VAPU    = 2.*VIS(L)*YPRO*PYN2
                SRK(L)  = SRK(L)  - VAPU*RK(L)
                DDEPS(L) = DDEPS(L)+ VAPU*RK(L)
                SEPS(L) = SEPS(L) - VAPU*REPS(L)*EXP(YEXP2)
 2000    CONTINUE

      RETURN
      END SUBROUTINE TURBKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TURBSP(VIS,RO,D3,RK,REPS,SRK,SEPS,VOL,A3,A3X,A3Y,A3Z,
     +     Y,IMAX,JMAX,KMAX,IN,JN,KN,KBOT,KTOP,KCP,CTA,ISTR,JSTR,KSTR,
     +     KX1,KX2,KY1,KY2,IDIR2)

      REAL :: RO(*),VIS(*),VOL(*),A3(*),A3X(*),A3Y(*),A3Z(*),D3(*),
     +        RK(*),REPS(*),SRK(*),SEPS(*),Y(*)

      INTEGER :: KCP(*)

C ... Calculate wall dependent source for Speziales ASM model
C ...  (ITURB=10)
      ISTRID  = IMAX + 2*IN
      IL      = KSTR

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 0

C ... WALL FUNCTIONS ACT ONLY TO THE HAFL OF THE CHANNEL**********
      KBOT2 = KBOT
      KZ1   = KLOW+1            !
      KZ2   = KUPP

      IF(IDIR2 == -1) THEN
         KBOT2 = KTOP
         KZ1   = KTOP-1
         KZ2   = KUPP+1
      ENDIF
C *************************************************************

C ... WALL FUNCTIONS ACT HOLE THE CHANNEL************************
      KBOT2 = 1
      KZ1   = 2
      KZ2   = KMAX
      IF(IDIR2 == -1) THEN
         KBOT2 = KMAX
         KZ1   = KMAX-1
         KZ2   = 1
      ENDIF
C *************************************************************

C *********************************************************************
C ... BOUNDARY NODE     NORMAL SIDE, INCREASING K ...
      KA        = (KN+KBOT2-1)*KSTR
      DO 800 J  = KY1,KY2
         IJ        = (JN+J-1)*JSTR + KA
         N1        = (JN+J-1)*ISTRID + IN 
      DO 800 I  = KX1,KX2
            L      = (IN+I-1)*ISTR + IJ + 1
            IGW    = I + N1
C     ----- THE EFFECT OF THE WALLS -----------------------------------
               Y(L)     = .5*D3(L)
               RY       = SQRT(RK(L)*RO(L))*Y(L)/VIS(L)
               YEXP     = MAX(-RY/CTA,-30.)
               F2       = 1. - EXP(YEXP)
               SEPS(L)  = SEPS(L) * F2  
 800     CONTINUE


      DO 2000 KG = KZ1,KZ2,IDIR2
         KA        = (KN+KG-1)*KSTR
C ...    NEXT LINE TELLS COMPILER TO IGNORE DATA DEPENDENCIES(-pfa)
C*$*ASSERTDO(CONCURRENT)
         DO 2000 J = 1,JMAX
         IJ        = (JN+J-1)*JSTR + KA
         N1        = (JN+J-1)*ISTRID + IN 
         DO 2000 I = 1,IMAX
            L      = (IN+I-1)*ISTR + IJ + 1
            IGW    = I + N1
C     ----- THE EFFECT OF THE WALLS -----------------------------------
         Y(L)     = Y(L-IDIR2*IL) + .5*(D3(L)+D3(L-IDIR2*IL))
               RY       = SQRT(RK(L)*RO(L))*Y(L)/VIS(L)
               YEXP     = MAX(-RY/CTA,-30.)
               F2       = 1. - EXP(YEXP)
               SEPS(L)  = SEPS(L)  * F2 
 2000    CONTINUE

      RETURN
      END SUBROUTINE TURBSP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SOURKE(SRK,SEPS,PTUR,VIS,OHMI,U,V,W,RO,RK,REPS,C,DDEPS,
     +     S,EPSTRU,VOL,PR,PRT,IMAX,JMAX,KMAX,IN,JN,KN,LAMIN,M,
     +     C1,C2,C21,ETA0,ITURB,IEPSMA,RI,IROTCO,TRANSL,TRM,NGL)

      USE TYPE_ARRAYS

      REAL :: U(*),RO(*),V(*),VIS(*),W(*),OHMI(*),VOL(*),C(*),RK(*),
     + REPS(*),SRK(*),SEPS(*),PTUR(*),S(*),EPSTRU(*),DDEPS(*),RI(*)

      REAL :: YEXP,APUF,RETUR

      LOGICAL :: TRANSL

      TYPE(INTERMITTENCY) :: TRM(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      EPS     = 1.E-20

      IF(IROTCO >= 1) THEN 
         RCOR=1.0
      ELSE
         RCOR=0.0
      ENDIF
C
C ... COMPUTE WALL-INDEPENDENT PRODUCTION TERMS FOR K-EPSILON MODEL
C
C ... Chien's 2-eq. model
      IF(ITURB /= 10) THEN
         DO 9000 KG = 1,KMAX
         IA      = (KG+KN-1)*IL + JN*ISTRID
         DO 9000 IG = 1,JMAX*ISTRID
            L       = IA + IG

C ... RICHARDSON NUMBER IS BOUNDED IN THE LOWER LIMIT
            RI(L)    = MAX(RI(L),-1.0E+10)
C ... TRUE EPSILON BASED ON EPSILON(WORM)
            EPSTRU(L) = REPS(L) + DDEPS(L) + EPS
            ETA     = MIN(S(L)*RK(L)/EPSTRU(L),100.)
            S(L)    = ETA
            
            RETUR   = MIN(RK(L)**2/(VIS(L)*(REPS(L)+EPS)),33.0)
            YEXP    = -.02778*RETUR**2

C ... Intermittency variables
            IF(TRANSL) THEN
*              SRK(L)  = SRK(L) + TRM(L)%G*PTUR(L) 
*     &                - MIN(MAX(TRM(L)%G, 0.1), 1.0)*REPS(L)
              SRK(L)  = SRK(L) + TRM(L)%GEFF*PTUR(L) 
     &                - MIN(MAX(TRM(L)%GEFF, 0.1), 1.0)*REPS(L)
            ELSE
              SRK(L)  = SRK(L) + PTUR(L) - REPS(L) 
            ENDIF 

C ... In Reynolds stress calculation PTUR can be negative.
C ... In dissipation equation this is not good thing. PPR 28.2.97
C ... After calculating cylinder, got a feeling that PTUR may be
C ... negative in dissipation eq. PPR 24.6.97 PPR
C           PTUR2   = MAX(0.,PTUR(L))
            PTUR2   = PTUR(L)
C ... The following does not work with the intel compiler with '-xW'
c           APUF    = C21*EXP(YEXP)
            APUF    = C21*2.7182818**YEXP
            SEPS(L) = SEPS(L) + REPS(L)/(RK(L)+EPS)*(C1*PTUR2
     +      - C2*(1.-0.2*RI(L)*RCOR)*REPS(L)*(1.- APUF))
c            T = MAX(RK(L)/REPS(L),SQRT(VIS(L)/REPS(L))) ! muutta tulosta
c     SEPS(L) = SEPS(L) + (C1*PTUR2 ! joskus
c     +      - C2*(1.-0.2*RI(L)*RCOR)*REPS(L)*(1.- C21*EXP(YEXP)))/T
9000  CONTINUE
C ... Speziale  et al. explicit ASM model
      ELSE
         DO 9100 KG = 1,KMAX
         IA      = (KG+KN-1)*IL + JN*ISTRID
         DO 9100 IG = 1,JMAX*ISTRID
            L       = IA + IG

C ... TRUE EPSILON BASED ON EPSILON(WORM)
            EPSTRU(L) = REPS(L) + DDEPS(L) + EPS
            ETA     = MIN(S(L)*RK(L)/EPSTRU(L),100.)
            S(L)    = ETA
            SRK(L)  = SRK(L) + PTUR(L) - REPS(L) 
            F2      = SEPS(L)
C ... In Reynolds stress calculation PTUR can be negative.
C ... In dissipation equation this is not good thing. PPR 28.2.97
            PTUR2   = MAX(0.,PTUR(L))
            SEPS(L) = REPS(L)/(RK(L)+EPS)*(C1*PTUR2 - C2*F2*REPS(L))
 9100    CONTINUE
      ENDIF

      IF(IEPSMA == 1) THEN
C ... COMPRESSIBLE DISSIPATION
      DO KG = 1,KMAX
         IA    = (KG+KN-1)*IL + JN*ISTRID
         DO IG = 1,JMAX*ISTRID
            L       = IA + IG
            XMA = 2.*RK(L)/RO(L)/C(L)**2 ! TURBULENT MACH NUMBER
            EPD = XMA*REPS(L)            ! COMPRESSIBLE DISSIPATION
            SRK(L) = SRK(L) - EPD
         ENDDO
      ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SOURKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SOURNG(SRK,SEPS,PTUR,VIS,OHMI,U,V,W,RO,RK,REPS,DDEPS,S,
     + EPSTRU,VOL,PR,PRT,IMAX,JMAX,KMAX,IN,JN,KN,LAMIN,M,
     + C1,C2,C21,ETA0)

      REAL :: U(*),RO(*),V(*),VIS(*),W(*),OHMI(*),VOL(*),RK(*),REPS(*),
     +        SRK(*),SEPS(*),PTUR(*),S(*),EPSTRU(*),DDEPS(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      EPS     = 1.E-20
C
C ... COMPUTE RNG- PRODUCTION TERM FOR THE EPSILON EQUATION
C
          DO 9000 KG = 1,KMAX
          IA      = (KG+KN-1)*IL + JN*ISTRID
          DO 9000 IG = 1,JMAX*ISTRID
          L       = IA + IG
          SEPS(L) = SEPS(L) + S(L)*(1.- S(L)/ETA0) / (1.+ .012*S(L)**3)*
     +              PTUR(L)*REPS(L)/(RK(L)+EPS)
9000  CONTINUE

      RETURN
      END SUBROUTINE SOURNG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SOURCD(PTUR,RK,REPS,DDEPS,EPS2,VIST,VIS,SRK,SEPS,
     +     RO,A1,A2,A3,VOL,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     C1,C2,C3,CMU,IMAX,JMAX,KMAX,IN,JN,KN,ITURB,M)

      REAL :: A1(*),A2(*),A1X(*),A1Y(*),A2X(*),A2Y(*),A3(*),A1Z(*),
     +        A2Z(*),A3X(*),A3Y(*),A3Z(*),VOL(*),
     +        PTUR(*),RK(*),REPS(*),EPS2(*),VIS(*),DDEPS(*),VIST(*),
     +        SRK(*),SEPS(*),RO(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
C ... CROSS-DIFFUSION MODEL FOR CHIEN KE
      DO 5000 K = 1,KMAX
      IA  = (KN+K-1)*IL + JN*ISTRID

      DO 1000 JJ = 1,JMAX*ISTRID
         L  = IA + JJ

         RKPR =   RK(L)/RO(L)
         REPR = REPS(L)/RK(L)
C ... XI-DIRECTION
         R1MP = RKPR +   RK(L-1)/RO(L-1)
         E1MP = REPR + REPS(L-1)/RK(L-1)
         R1PP = RKPR +   RK(L+1)/RO(L+1)
         E1PP = REPR + REPS(L+1)/RK(L+1)
C ... ETA-DIRECTION
         R2MP = RKPR +   RK(L-ISTRID)/RO(L-ISTRID)
         E2MP = REPR + REPS(L-ISTRID)/RK(L-ISTRID)
         R2PP = RKPR +   RK(L+ISTRID)/RO(L+ISTRID)
         E2PP = REPR + REPS(L+ISTRID)/RK(L+ISTRID)
C ... ZETA DIRECTION
         R3MP = RKPR +   RK(L-IL)/RO(L-IL)
         E3MP = REPR + REPS(L-IL)/RK(L-IL)
         R3PP = RKPR +   RK(L+IL)/RO(L+IL)
         E3PP = REPR + REPS(L+IL)/RK(L+IL)

         PVOL    = .5/VOL(L)
         DKDXI   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)    *R3PP - A3(L)*A3X(L)*R3MP +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*R2PP - A2(L)*A2X(L)*R2MP +
     3        A1(L+1)     *A1X(L+1)     *R1PP - A1(L)*A1X(L)*R1MP)
         DKDYI     = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)    *R3PP - A3(L)*A3Y(L)*R3MP +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*R2PP - A2(L)*A2Y(L)*R2MP +
     3        A1(L+1)     *A1Y(L+1)     *R1PP - A1(L)*A1Y(L)*R1MP)
         DKDZI     = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)    *R3PP - A3(L)*A3Z(L)*R3MP +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*R2PP - A2(L)*A2Z(L)*R2MP +
     3        A1(L+1)     *A1Z(L+1)     *R1PP - A1(L)*A1Z(L)*R1MP)
         DEDXI     = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)    *E3PP - A3(L)*A3X(L)*E3MP +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*E2PP - A2(L)*A2X(L)*E2MP +
     3        A1(L+1)     *A1X(L+1)     *E1PP - A1(L)*A1X(L)*E1MP)
         DEDYI   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)    *E3PP - A3(L)*A3Y(L)*E3MP +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*E2PP - A2(L)*A2Y(L)*E2MP +
     3        A1(L+1)     *A1Y(L+1)     *E1PP - A1(L)*A1Y(L)*E1MP)
         DEDZI     = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)    *E3PP - A3(L)*A3Z(L)*E3MP +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*E2PP - A2(L)*A2Z(L)*E2MP +
     3        A1(L+1)     *A1Z(L+1)     *E1PP - A1(L)*A1Z(L)*E1MP)

         SUM = DKDXI*DEDXI + DKDYI*DEDYI + DKDZI*DEDZI
         SUM = MIN(0.,SUM)
c         SEPS(L) = SEPS(L) + C3*VIST(L)*RK(L)/RO(L)*SUM
c         SEPS(L) = SEPS(L) + C3*VIST(L)*SUM
         SEPS(L) = SEPS(L) + C3*CMU*RK(L)**2/(REPS(L)+DDEPS(L))*SUM
c ... tama vaihtoehdoista tuntu toimivan parhaiten
c         xpu = C3*CMU*RK(L)**2/(REPS(L)+DDEPS(L))*SUM
c         call ijkpai(l,48,32,1,ll,mm,nn)
c         if(ll == 40 .and. nn == 1) then
c            write(77,*) mm,xpu,DKDYI,DEDYI
c         endif
 1000 CONTINUE
 5000 CONTINUE

      RETURN
      END SUBROUTINE SOURCD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE VORFUN(OHMI,OMX,OMY,OMZ,S,S11,S12,S13,S22,S23,S33,
     2 PTUR,U,V,W,RK,REPS,DDEPS,EPS2,VIST,VIS,A1,A2,A3,VOL,A1X,A1Y,
     3 A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,D1,D2,D3,XC,YC,ZC,IMAX,JMAX,KMAX,
     4 IN,JN,KN,ITURB,M,STRESL,ZZZ,PHITUR,MAXW,IDERI,BIJ,MAXEB,ISTRES)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ITURB,M,MAXW,MAW,I,J,K,L,JJ,IL,
     2           KA,IA,LLL,MMM,NNN,ISTRID,JSTRID,KSTRID,MAXEB

      INTEGER :: IALA,IYLA,JALA,JYLA,KALA,KYLA,N,IDERI,ISTRES

      INTEGER, DIMENSION(6) :: IS

      REAL :: U(*),VOL(*),V(*),A1(*),A2(*),A1X(*),A1Y(*),A2X(*),
     2 A2Y(*),W(*),A3(*),A1Z(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),OHMI(*),
     3 OMX(*),OMY(*),OMZ(*),PTUR(*),RK(*),REPS(*),EPS2(*),VIS(*),S(*),
     4 S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),DDEPS(*),VIST(*),
     5 D1(*),D2(*),D3(*),BIJ(MAXEB,*),PHITUR(*)

      REAL :: XC(*), YC(*), ZC(*)

      REAL, POINTER :: DUDY(:),DUDZ(:),DVDX(:),DVDZ(:),DWDY(:),DWDX(:),
     2      DUDX(:),DVDY(:),DWDZ(:)
      REAL, TARGET ::  ZZZ(MAXW)

      REAL :: U1M,V1M,W1M,U1P,V1P,W1P,U2M,V2M,W2M,U2P,V2P,W2P,U3M,V3M,
     2     W3M,U3P,V3P,W3P,PVOL,DUDXI,DUDYI,DUDZI,DVDXI,DVDYI,DVDZI,
     3     DWDXI,DWDYI,DWDZI,XTRACE,DIV,R2O3,ALF1,ALF2,ALF3,ALF4,ALF5,
     4     C1,C2,C3,C4,C5,GG,TERM1,TERM2,TERM3,DAMP,DAMP2,SQR2,EPS,
     5     XSI,CMUT,VMUT,ETA,EPSI,RKPE,RKPE2,P3,RMUT,DIVVP3

      REAL :: R11,R22,R33,R12,R13,R23,SIJ11,SIJ12,SIJ13,SIJ22,SIJ23,
     2     SIJ33,UIJ11,UIJ12,UIJ13,UIJ21,UIJ22,UIJ23,UIJ31,UIJ32,
     3     UIJ33,HUU,HUV,HUW,HVV,HVW,HWW

      REAL :: W1,W2,W3,DX,DY,DZ,BETA,ALFA1,ALFA2,ALFA3

      LOGICAL :: STRESL

*      CALL DOMAW(MAXW,IMAX,JMAX,MAW)
*      CALL DOMAW(MAXW,IMAX,JMAX,IN,JN,MAW)

*      DUDY => ZZZ( 0*MAW+1: 1*MAW);DUDZ=> ZZZ( 1*MAW+1: 2*MAW)
*      DVDX => ZZZ( 2*MAW+1: 3*MAW);DVDZ=> ZZZ( 3*MAW+1: 4*MAW)
*      DWDX => ZZZ( 4*MAW+1: 5*MAW);DWDY=> ZZZ( 4*MAW+1: 5*MAW)
*      DUDX => ZZZ( 6*MAW+1: 7*MAW);DVDY=> ZZZ( 6*MAW+1: 7*MAW)
*      DWDZ => ZZZ( 8*MAW+1: 9*MAW)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID

      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghostcells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF

      R2O3    = 2./3.
      P3      = 1.0/3.0      

C ... Constants for a Speziale's model (ITURB=10) AIAA Journal, Vol.33, No. 10.

      ALF1    = .375
      ALF2    = .116
      ALF3    = .108

      SQR2    = 1./SQRT(2.)
      EPS     = 1.E-20

C ... LRR mallin vakiot

      C1 = 3.
      C2 = 0.8
      C3 = 1.75
      C4 = 1.31
      C5 = 2.

C ... SSG mallin vakiot

      C1 = 6.8
      C2 = 0.36
      C3 = 1.25
      C4 = 0.4
      C5 = 1.88

      GG    = 1./(C1/2. + C5 - 1.)
      ALF1  = (4./3.-C2)*GG/2.
      ALF2  = (2. - C3)**2*GG**2/4.
      ALF3  = (2. - C4)**2*GG**2/4.
      ALF4  = (2. - C4)*GG/2.
      ALF5  = (2. - C3)*GG

      DO L = 1,ISTRID*JSTRID*KSTRID
         S(L) = 0.
      ENDDO


C ... Three different approaches are available for calculating gradients:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach


      DO K = KALA,KYLA
         KA = (KN+K-1)*IL
         DO J = JALA,JYLA
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I = IALA,IYLA
               L = JJ + I

               SELECT CASE(IDERI)

               CASE(1)  ! Gauss approach weighted by distances

C ... XI-DIRECTION

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2  = 2.0 - W1

                  U1M = A1(L)  *(W2*U(L) + W1*U(L-1))
                  V1M = A1(L)  *(W2*V(L) + W1*V(L-1))
                  W1M = A1(L)  *(W2*W(L) + W1*W(L-1))

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) )
                  W2  = 2.0 - W1

                  U1P = A1(L+1)*(W2*U(L) + W1*U(L+1))
                  V1P = A1(L+1)*(W2*V(L) + W1*V(L+1))
                  W1P = A1(L+1)*(W2*W(L) + W1*W(L+1))

C ... ETA-DIRECTION

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) )
                  W2  = 2.0 - W1

                  U2M = A2(L)       *(W2*U(L) + W1*U(L-ISTRID))
                  V2M = A2(L)       *(W2*V(L) + W1*V(L-ISTRID))
                  W2M = A2(L)       *(W2*W(L) + W1*W(L-ISTRID))

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) )
                  W2  = 2.0 - W1

                  U2P = A2(L+ISTRID)*(W2*U(L) + W1*U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(W2*V(L) + W1*V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W2*W(L) + W1*W(L+ISTRID))

C ... ZETA DIRECTION

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) )
                  W2  = 2.0 - W1

                  U3M = A3(L)   *(W2*U(L) + W1*U(L-IL))
                  V3M = A3(L)   *(W2*V(L) + W1*V(L-IL))
                  W3M = A3(L)   *(W2*W(L) + W1*W(L-IL))

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) )
                  W2  = 2.0 - W1

                  U3P = A3(L+IL)*(W2*U(L) + W1*U(L+IL))
                  V3P = A3(L+IL)*(W2*V(L) + W1*V(L+IL))
                  W3P = A3(L+IL)*(W2*W(L) + W1*W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)

               CASE(2)  ! Simple Gaussian derivative calculation

C ... XI-DIRECTION

                  U1M = A1(L)  *(U(L) + U(L-1))
                  V1M = A1(L)  *(V(L) + V(L-1))
                  W1M = A1(L)  *(W(L) + W(L-1))
                  U1P = A1(L+1)*(U(L) + U(L+1))
                  V1P = A1(L+1)*(V(L) + V(L+1))
                  W1P = A1(L+1)*(W(L) + W(L+1))

C ... ETA-DIRECTION

                  U2M = A2(L)       *(U(L) + U(L-ISTRID))
                  V2M = A2(L)       *(V(L) + V(L-ISTRID))
                  W2M = A2(L)       *(W(L) + W(L-ISTRID))
                  U2P = A2(L+ISTRID)*(U(L) + U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(V(L) + V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W(L) + W(L+ISTRID))

C ... ZETA DIRECTION

                  U3M = A3(L)   *(U(L) + U(L-IL))
                  V3M = A3(L)   *(V(L) + V(L-IL))
                  W3M = A3(L)   *(W(L) + W(L-IL))
                  U3P = A3(L+IL)*(U(L) + U(L+IL))
                  V3P = A3(L+IL)*(V(L) + V(L+IL))
                  W3P = A3(L+IL)*(W(L) + W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)

               CASE(3)  ! Least-squares approach

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001

                  IS(1) = -1
                  IS(2) = -ISTRID
                  IS(3) = -IL
                  IS(4) =  1
                  IS(5) =  ISTRID
                  IS(6) =  IL

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO

                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))

                  BETA = (R12*R23 - R13*R22) / (R11*R22)

                  DUDXI = 0.0
                  DUDYI = 0.0
                  DUDZI = 0.0
                  DVDXI = 0.0
                  DVDYI = 0.0
                  DVDZI = 0.0
                  DWDXI = 0.0
                  DWDYI = 0.0
                  DWDZI = 0.0

                  DO N = 1,6
                     DX    = XC(L + IS(N)) - XC(L)
                     DY    = YC(L + IS(N)) - YC(L)
                     DZ    = ZC(L + IS(N)) - ZC(L)
                     ALFA1 =  DX / R11**2
                     ALFA2 = (DY - R12*DX/R11) / R22**2
                     ALFA3 = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1    = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2    = ALFA2 - R23*ALFA3/R22
                     W3    = ALFA3
                     DUDXI = DUDXI + W1*(U(L+IS(N))-U(L))
                     DUDYI = DUDYI + W2*(U(L+IS(N))-U(L))
                     DUDZI = DUDZI + W3*(U(L+IS(N))-U(L))
                     DVDXI = DVDXI + W1*(V(L+IS(N))-V(L))
                     DVDYI = DVDYI + W2*(V(L+IS(N))-V(L))
                     DVDZI = DVDZI + W3*(V(L+IS(N))-V(L))
                     DWDXI = DWDXI + W1*(W(L+IS(N))-W(L))
                     DWDYI = DWDYI + W2*(W(L+IS(N))-W(L))
                     DWDZI = DWDZI + W3*(W(L+IS(N))-W(L))
                  ENDDO

               END SELECT  ! IDERI

                  OMZ(L)  = .5*(DUDYI - DVDXI)
                  OMY(L)  = .5*(DUDZI - DWDXI)
                  OMX(L)  = .5*(DVDZI - DWDYI)
                  OHMI(L) = 2.*SQRT(OMX(L)**2 + OMY(L)**2 + OMZ(L)**2)

                  IF(STRESL) THEN

                     UIJ11  = DUDXI
                     UIJ12  = DUDYI
                     UIJ13  = DUDZI
         
                     UIJ21  = DVDXI
                     UIJ22  = DVDYI
                     UIJ23  = DVDZI
         
                     UIJ31  = DWDXI
                     UIJ32  = DWDYI
                     UIJ33  = DWDZI
         
                     DIV    = R2O3*(DUDXI + DVDYI + DWDZI)
                     SIJ11  = UIJ11 + UIJ11 - DIV
                     SIJ12  = UIJ12 + UIJ21
                     SIJ13  = UIJ13 + UIJ31
                     SIJ22  = UIJ22 + UIJ22 - DIV
                     SIJ23  = UIJ23 + UIJ32
                     SIJ33  = UIJ33 + UIJ33 - DIV
     
C ... Calculate strain
     
                     S(L)   = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +                      + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +                      + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
                     S(L)   = SQRT(S(L))

                     S11(L) = .5*SIJ11
                     S12(L) = .5*SIJ12
                     S13(L) = .5*SIJ13
                     S22(L) = .5*SIJ22
                     S23(L) = .5*SIJ23
                     S33(L) = .5*SIJ33

                  ENDIF

                  IF(ITURB == 10 .AND. STRESL) THEN ! Speziale's model

                     EPSI    = REPS(L) + DDEPS(L)
                     RKPE    = RK(L)/(EPSI+EPS)
                     RKPE2   = RK(L)**2/(EPSI+EPS)
                     ETA     = ALF2*(S(L)*RKPE)**2
                     XSI     = ALF3*(SQR2*OHMI(L)*RKPE)**2
                     CMUT    = 3.*(1. + ETA)/(3. + ETA 
     +                       + 6.*XSI*ETA + 6.*XSI)*ALF1
                     VMUT    = CMUT*RKPE2
                     TERM1   = 2.     *VMUT
                     TERM2   = 2.*ALF4*VMUT*RKPE
                     TERM3   = 2.*ALF5*VMUT*RKPE         

                     R11     = R2O3*RK(L) 
     +                       - (TERM1*  S11(L) 
     +                       +  TERM2*( S12(L)*OMZ(L) - S13(L)*OMY(L))
     +                       -  TERM3*( S11(L)**2 + S12(L)**2 
     +                       +          S13(L)**2 - S(L)**2/3.))

                     R22     = R2O3*RK(L) 
     +                       - (TERM1*  S22(L)
     +                       +  TERM2*(-S12(L)*OMZ(L) + S23(L)*OMX(L))
     +                       -  TERM3*( S12(L)**2 + S22(L)**2 
     +                       +          S23(L)**2 - S(L)**2/3.))

                     R33     = R2O3*RK(L) 
     +                       - (TERM1*  S33(L)
     +                       +  TERM2*( S13(L)*OMY(L) - S23(L)*OMX(L))
     +                       -  TERM3*( S13(L)**2 + S23(L)**2 
     +                       +          S33(L)**2 - S(L)**2/3.))

                     R12     =            
     +                       - (TERM1*S12(L) 
     +                       +  TERM2*.5*(-S11(L)*OMZ(L)+S13(L)*OMX(L) 
     +                       +             S22(L)*OMZ(L)-S23(L)*OMY(L))
     +                       -  TERM3*    (S11(L)*S12(L) 
     +                       +             S12(L)*S22(L) 
     +                       +             S13(L)*S23(L)))

                     R13     =            
     +                       - (TERM1*S13(L) 
     +                       +  TERM2*.5*( S11(L)*OMY(L)-S12(L)*OMX(L) 
     +                       +             S23(L)*OMZ(L)-S33(L)*OMY(L))
     +                       -  TERM3*    (S11(L)*S13(L) 
     +                       +             S12(L)*S23(L) 
     +                       +             S13(L)*S33(L)))

                     R23     =            
     +                       - (TERM1*S23(L) 
     +                       +  TERM2*.5*( S12(L)*OMY(L)-S22(L)*OMX(L) 
     +                       -             S13(L)*OMZ(L)-S33(L)*OMX(L))
     +                       -  TERM3*    (S12(L)*S13(L) 
     +                       +             S22(L)*S23(L) 
     +                       +             S23(L)*S33(L)))

                    PTUR(L) = -(
     +                           R11*DUDXI+R12*(DUDYI+DVDXI)+
     +                           R22*DVDYI+R23*(DVDZI+DWDYI)+
     +                           R33*DWDZI+R13*(DWDXI+DUDZI))

C ... Menter's modification for the maximum production
         
                     PTUR(L) = MIN(PTUR(L),20.*EPSI)

                     S11(L)  = R11
                     S12(L)  = R12
                     S13(L)  = R13
                     S22(L)  = R22
                     S23(L)  = R23
                     S33(L)  = R33

                  ENDIF

                  IF(ITURB == 11 .AND. STRESL) THEN ! Something else

                     EPSI    = REPS(L) + DDEPS(L)
                     RKPE    = RK(L)/(EPSI+EPS)
                     RKPE2   = RK(L)**2/(EPSI+EPS)
                     ETA     = (  .5*ALF3/ALF1*   S(L)*RKPE)**2
                     XSI     = (SQR2*ALF2/ALF1*OHMI(L)*RKPE)**2
                     DAMP2   = 3.*(1. + ETA)/(3. + ETA 
     +                       + 6.*XSI*ETA + 6.*XSI)*RKPE
                     DAMP    = 2.*VIST(L)/ALF1/RK(L) ! Chien's model

                     TERM1   = DAMP*ALF1*RK(L)
                     TERM2   = DAMP*ALF2*RKPE2
                     TERM3   = DAMP*ALF3*RKPE2

                     R11     = R2O3*RK(L) 
     +                       - (TERM1*  S11(L) 
     +                       +  TERM2*( S12(L)*OMZ(L) - S13(L)*OMY(L))
     +                       -  TERM3*( S11(L)**2 + S12(L)**2 
     +                       +          S13(L)**2 - S(L)**2/3.))

                     R22     = R2O3*RK(L) 
     +                       - (TERM1*  S22(L)
     +                       +  TERM2*(-S12(L)*OMZ(L) + S23(L)*OMX(L))
     +                       -  TERM3*( S12(L)**2 + S22(L)**2 
     +                       +          S23(L)**2 - S(L)**2/3.))

                     R33     = R2O3*RK(L) 
     +                       - (TERM1*  S33(L)
     +                       +  TERM2*( S13(L)*OMY(L) - S23(L)*OMX(L))
     +                       -  TERM3*( S13(L)**2 + S23(L)**2 
     +                       +          S33(L)**2 - S(L)**2/3.))

                     R12     =            
     +                       - (TERM1*S12(L) 
     +                       +  TERM2*.5*(-S11(L)*OMZ(L)+S13(L)*OMX(L) 
     +                       +             S22(L)*OMZ(L)-S23(L)*OMY(L))
     +                       -  TERM3*    (S11(L)*S12(L) 
     +                       +             S12(L)*S22(L) 
     +                       +             S13(L)*S23(L)))

                     R13     =            
     +                       - (TERM1*S13(L) 
     +                       +  TERM2*.5*( S11(L)*OMY(L)-S12(L)*OMX(L) 
     +                       +             S23(L)*OMZ(L)-S33(L)*OMY(L))
     +                       -  TERM3*    (S11(L)*S13(L) 
     +                       +             S12(L)*S23(L) 
     +                       +             S13(L)*S33(L)))

                     R23     =            
     +                       - (TERM1*S23(L) 
     +                       +  TERM2*.5*( S12(L)*OMY(L)-S22(L)*OMX(L) 
     +                       -             S13(L)*OMZ(L)-S33(L)*OMX(L))
     +                       -  TERM3*    (S12(L)*S13(L) 
     +                       +             S22(L)*S23(L) 
     +                       +             S23(L)*S33(L)))

                     PTUR(L) = -(
     +                           R11*DUDXI+R12*(DUDYI+DVDXI)+
     +                           R22*DVDYI+R23*(DVDZI+DWDYI)+
     +                           R33*DWDZI+R13*(DWDXI+DUDZI))

C ... Menter's modification for the maximum production
         
                     PTUR(L) = MIN(PTUR(L),20.*EPSI)

                     S11(L)  = R11
                     S12(L)  = R12
                     S13(L)  = R13
                     S22(L)  = R22
                     S23(L)  = R23
                     S33(L)  = R33

                  ENDIF

C ... EARSM of Wallin and Johansson

                  IF(ISTRES >= 1) THEN
                  
                     RMUT    = MAX(EPS2(L) - 1.0 , EPS)*VIS(L)
c                    DIVVP3  = P3*(S11(L) + S22(L) + S33(L))
                     DIVVP3  = 0.   ! Was already taken into account

C     Half of the Reynolds stresses (with rho) 0.5*rho*u_i*u_j

                     HUU     = MAX((RK(L)*(P3 + BIJ(L,1))
     +                       - RMUT*(S11(L) - DIVVP3)) , 0.0)
                     HVV     = MAX((RK(L)*(P3 + BIJ(L,4))
     +                       - RMUT*(S22(L) - DIVVP3)) , 0.0)
                     HWW     = MAX((RK(L)*(P3 - BIJ(L,1) - BIJ(L,4))
     +                       - RMUT*(S33(L) - DIVVP3)) , 0.0)
                     HUV     = RK(L)*BIJ(L,2)  - RMUT*S12(L)
                     HUW     = RK(L)*BIJ(L,3)  - RMUT*S13(L) 
                     HVW     = RK(L)*BIJ(L,5)  - RMUT*S23(L)

C     Production of turbulence

                     PTUR(L) =-2.0*(HUU*S11(L)+HVV*S22(L)+HWW*S33(L)
     +                       + 2.0*(HUV*S12(L)+HUW*S13(L)+HVW*S23(L)))
            
C     Forced laminarity by limited production if PHITUR is zero
                     PTUR(L) = PTUR(L)*PHITUR(L)

                     S11(L)  = 2.*HUU
                     S12(L)  = 2.*HUV
                     S13(L)  = 2.*HUW
                     S22(L)  = 2.*HVV
                     S23(L)  = 2.*HVW
                     S33(L)  = 2.*HWW

                  ENDIF

            ENDDO
         ENDDO
      ENDDO

C **********************************************************************

C ... Finally change the tensor components into vector components

      DO K = KALA,KYLA
         KA = (KN+K-1)*IL
         DO J = JALA,JYLA
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I = IALA,IYLA
               L = JJ + I
               OMX(L) = -2.*OMX(L) ! Sign changed 28.6.2019
               OMY(L) =  2.*OMY(L)
               OMZ(L) = -2.*OMZ(L)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE VORFUN 
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DERFUN(OHMI,W12,W13,W23,S11,S12,S13,S22,S23,S33,
     +     U,V,W,A1,A2,A3,VOL,D1,D2,D3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     +     A3X,A3Y,A3Z,XC,YC,ZC,IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,
     +     STRAIN,VELLAP)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: U(*),VOL(*),V(*),A1(*),A2(*),A1X(*),A1Y(*),A2X(*),A2Y(*),
     2 W(*),A3(*),A1Z(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),OHMI(*),W12(*),
     3 W13(*),W23(*),S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),
     4 D1(*),D2(*),D3(*),STRAIN(*),VELLAP(*)

      REAL :: R2O3,SQR2,U1M,V1M,W1M,U1P,V1P,W1P,U2M,V2M,W2M,U2P,V2P,W2P,
     2 U3M,V3M,W3M,U3P,V3P,W3P,PVOL,DUDXI,DUDYI,DUDZI,DVDXI,DVDYI,
     3 DVDZI,DWDXI,DWDYI,DWDZI,DIV,SIJ11,SIJ22,SIJ33,SIJ12,SIJ13,SIJ23,
     4 D2UDX2,D2VDX2,D2WDX2

      REAL :: W1,W2,W3,R11,R12,R22,R13,R23,R33,DX,DY,DZ,BETA,
     2        ALFA1,ALFA2,ALFA3

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     2           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA

      INTEGER, DIMENSION(6) :: IS

C ... Calculate velocity gradients and store them as vorticity and strain

C ... Three different approaches are available:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID
      
      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghostcells are included on the first level
c         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
c         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF

      R2O3    = 2./3.
      SQR2    = SQRT(2.)

      SELECT CASE(IDERI)

      CASE(1)  ! Gauss approach weighted by distances

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2  = 2.0 - W1

                  U1M = A1(L)  *(W2*U(L) + W1*U(L-1))
                  V1M = A1(L)  *(W2*V(L) + W1*V(L-1))
                  W1M = A1(L)  *(W2*W(L) + W1*W(L-1))

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) )
                  W2  = 2.0 - W1

                  U1P = A1(L+1)*(W2*U(L) + W1*U(L+1))
                  V1P = A1(L+1)*(W2*V(L) + W1*V(L+1))
                  W1P = A1(L+1)*(W2*W(L) + W1*W(L+1))

C ... ETA-DIRECTION

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) )
                  W2  = 2.0 - W1

                  U2M = A2(L)       *(W2*U(L) + W1*U(L-ISTRID))
                  V2M = A2(L)       *(W2*V(L) + W1*V(L-ISTRID))
                  W2M = A2(L)       *(W2*W(L) + W1*W(L-ISTRID))

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) )
                  W2  = 2.0 - W1

                  U2P = A2(L+ISTRID)*(W2*U(L) + W1*U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(W2*V(L) + W1*V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W2*W(L) + W1*W(L+ISTRID))

C ... ZETA DIRECTION

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) )
                  W2  = 2.0 - W1

                  U3M = A3(L)   *(W2*U(L) + W1*U(L-IL))
                  V3M = A3(L)   *(W2*V(L) + W1*V(L-IL))
                  W3M = A3(L)   *(W2*W(L) + W1*W(L-IL))

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) )
                  W2  = 2.0 - W1

                  U3P = A3(L+IL)*(W2*U(L) + W1*U(L+IL))
                  V3P = A3(L+IL)*(W2*V(L) + W1*V(L+IL))
                  W3P = A3(L+IL)*(W2*W(L) + W1*W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)
              
C ... CALCULATE VORTICITIES

                  W12(L)  = .5*(DUDYI - DVDXI)
                  W13(L)  = .5*(DUDZI - DWDXI)
                  W23(L)  = .5*(DVDZI - DWDYI)
                  OHMI(L) = 2.*SQRT(W12(L)**2 + W13(L)**2 + W23(L)**2)

                  SIJ11   = .5*(DUDXI + DUDXI) ! - DIV) obs!
                  SIJ22   = .5*(DVDYI + DVDYI) ! - DIV)
                  SIJ33   = .5*(DWDZI + DWDZI) ! - DIV)
                  SIJ12   = .5*(DUDYI + DVDXI)
                  SIJ13   = .5*(DUDZI + DWDXI)
                  SIJ23   = .5*(DVDZI + DWDYI)
C     
C ... CALCULATE STRAIN IF NEEDED
C     
      STRAIN(L)  = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +           + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +           + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
      STRAIN(L)  = SQRT(2.*STRAIN(L))

                  S11(L)= SIJ11
                  S12(L)= SIJ12
                  S13(L)= SIJ13
                  S22(L)= SIJ22
                  S23(L)= SIJ23
                  S33(L)= SIJ33

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(2) ! Simple Gaussian derivative calculation

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  U1M = A1(L)  *(U(L) + U(L-1))
                  V1M = A1(L)  *(V(L) + V(L-1))
                  W1M = A1(L)  *(W(L) + W(L-1))
                  U1P = A1(L+1)*(U(L) + U(L+1))
                  V1P = A1(L+1)*(V(L) + V(L+1))
                  W1P = A1(L+1)*(W(L) + W(L+1))

C ... ETA-DIRECTION

                  U2M = A2(L)       *(U(L) + U(L-ISTRID))
                  V2M = A2(L)       *(V(L) + V(L-ISTRID))
                  W2M = A2(L)       *(W(L) + W(L-ISTRID))
                  U2P = A2(L+ISTRID)*(U(L) + U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(V(L) + V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W(L) + W(L+ISTRID))

C ... ZETA DIRECTION

                  U3M = A3(L)   *(U(L) + U(L-IL))
                  V3M = A3(L)   *(V(L) + V(L-IL))
                  W3M = A3(L)   *(W(L) + W(L-IL))
                  U3P = A3(L+IL)*(U(L) + U(L+IL))
                  V3P = A3(L+IL)*(V(L) + V(L+IL))
                  W3P = A3(L+IL)*(W(L) + W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)
              
C ... CALCULATE VORTICITIES

                  W12(L)  = .5*(DUDYI - DVDXI)
                  W13(L)  = .5*(DUDZI - DWDXI)
                  W23(L)  = .5*(DVDZI - DWDYI)
                  OHMI(L) = 2.*SQRT(W12(L)**2 + W13(L)**2 + W23(L)**2)

                  SIJ11   = .5*(DUDXI + DUDXI) ! - DIV) obs!
                  SIJ22   = .5*(DVDYI + DVDYI) ! - DIV)
                  SIJ33   = .5*(DWDZI + DWDZI) ! - DIV)
                  SIJ12   = .5*(DUDYI + DVDXI)
                  SIJ13   = .5*(DUDZI + DWDXI)
                  SIJ23   = .5*(DVDZI + DWDYI)
C     
C ... CALCULATE STRAIN IF NEEDED
C     
      STRAIN(L)  = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +      + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +      + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
      STRAIN(L)  = SQRT(2.*STRAIN(L))

                  S11(L)= SIJ11
                  S12(L)= SIJ12
                  S13(L)= SIJ13
                  S22(L)= SIJ22
                  S23(L)= SIJ23
                  S33(L)= SIJ33

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(3) ! Least-squares approach

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001

         IS(1) = -1
         IS(2) = -ISTRID
         IS(3) = -IL
         IS(4) =  1
         IS(5) =  ISTRID
         IS(6) =  IL

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO

                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))

                  BETA = (R12*R23 - R13*R22) / (R11*R22)

                  DUDXI = 0.0
                  DUDYI = 0.0
                  DUDZI = 0.0
                  DVDXI = 0.0
                  DVDYI = 0.0
                  DVDZI = 0.0
                  DWDXI = 0.0
                  DWDYI = 0.0
                  DWDZI = 0.0

                  DO N = 1,6
                     DX    = XC(L + IS(N)) - XC(L)
                     DY    = YC(L + IS(N)) - YC(L)
                     DZ    = ZC(L + IS(N)) - ZC(L)
                     ALFA1 =  DX / R11**2
                     ALFA2 = (DY - R12*DX/R11) / R22**2
                     ALFA3 = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1    = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2    = ALFA2 - R23*ALFA3/R22
                     W3    = ALFA3
                     DUDXI = DUDXI + W1*(U(L+IS(N))-U(L))
                     DUDYI = DUDYI + W2*(U(L+IS(N))-U(L))
                     DUDZI = DUDZI + W3*(U(L+IS(N))-U(L))
                     DVDXI = DVDXI + W1*(V(L+IS(N))-V(L))
                     DVDYI = DVDYI + W2*(V(L+IS(N))-V(L))
                     DVDZI = DVDZI + W3*(V(L+IS(N))-V(L))
                     DWDXI = DWDXI + W1*(W(L+IS(N))-W(L))
                     DWDYI = DWDYI + W2*(W(L+IS(N))-W(L))
                     DWDZI = DWDZI + W3*(W(L+IS(N))-W(L))
                  ENDDO
              
C ... CALCULATE VORTICITIES

                  W12(L)  = .5*(DUDYI - DVDXI)
                  W13(L)  = .5*(DUDZI - DWDXI)
                  W23(L)  = .5*(DVDZI - DWDYI)
                  OHMI(L) = 2.*SQRT(W12(L)**2 + W13(L)**2 + W23(L)**2)

                  SIJ11   = .5*(DUDXI + DUDXI) ! - DIV) obs!
                  SIJ22   = .5*(DVDYI + DVDYI) ! - DIV)
                  SIJ33   = .5*(DWDZI + DWDZI) ! - DIV)
                  SIJ12   = .5*(DUDYI + DVDXI)
                  SIJ13   = .5*(DUDZI + DWDXI)
                  SIJ23   = .5*(DVDZI + DWDYI)
C     
C ... CALCULATE STRAIN
C     
      STRAIN(L)  = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +      + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +      + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
      STRAIN(L)  = SQRT(2.*STRAIN(L))

                  S11(L) = SIJ11
                  S12(L) = SIJ12
                  S13(L) = SIJ13
                  S22(L) = SIJ22
                  S23(L) = SIJ23
                  S33(L) = SIJ33

               ENDDO
            ENDDO
         ENDDO

      END SELECT

C **********************************************************************

C ... Calculate the Laplacian of the velocity vector
C ... Thin-layer approximation!

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  U1M = A1(L)**2  *(U(L) - U(L-1))
                  V1M = A1(L)**2  *(V(L) - V(L-1))
                  W1M = A1(L)**2  *(W(L) - W(L-1))
                  U1P = A1(L+1)**2*(U(L) - U(L+1))
                  V1P = A1(L+1)**2*(V(L) - V(L+1))
                  W1P = A1(L+1)**2*(W(L) - W(L+1))

C ... ETA-DIRECTION

                  U2M = A2(L)**2       *(U(L) - U(L-ISTRID))
                  V2M = A2(L)**2       *(V(L) - V(L-ISTRID))
                  W2M = A2(L)**2       *(W(L) - W(L-ISTRID))
                  U2P = A2(L+ISTRID)**2*(U(L) - U(L+ISTRID))
                  V2P = A2(L+ISTRID)**2*(V(L) - V(L+ISTRID))
                  W2P = A2(L+ISTRID)**2*(W(L) - W(L+ISTRID))

C ... ZETA DIRECTION

                  U3M = A3(L)**2   *(U(L) - U(L-IL))
                  V3M = A3(L)**2   *(V(L) - V(L-IL))
                  W3M = A3(L)**2   *(W(L) - W(L-IL))
                  U3P = A3(L+IL)**2*(U(L) - U(L+IL))
                  V3P = A3(L+IL)**2*(V(L) - V(L+IL))
                  W3P = A3(L+IL)**2*(W(L) - W(L+IL))

                  PVOL  = 1./VOL(L)**2 ! These are not averages !

                  D2UDX2 =-PVOL*(U3P - U3M +
     2                           U2P - U2M +
     3                           U1P - U1M)
                  D2VDX2 =-PVOL*(V3P - V3M +
     2                           V2P - V2M +
     3                           V1P - V1M)
                  D2WDX2 =-PVOL*(W3P - W3M +
     2                           W2P - W2M +
     3                           W1P - W1M)  
                  VELLAP(L) = SQRT(D2UDX2**2 + D2VDX2**2 + D2WDX2**2) 

               ENDDO
            ENDDO
         ENDDO



C **********************************************************************


C ... On the coarse levels extend the values to the ghostcells
C ... This affects only the convergence rate and is not very important

      IF(M >= 2) THEN
         CALL EXTEND(OHMI,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(W12,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(W13,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(W23,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S11,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S12,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S13,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S22,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S23,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(S33,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(VELLAP,IMAX,JMAX,KMAX,IN,JN,KN)
      ENDIF      

      RETURN
      END SUBROUTINE DERFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TIJFUN(TIJ,RO,U,V,W,A1,A2,A3,VOL,D1,D2,D3,
     +   A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     +   IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,NGL)

C ... Lighthill's tensor gradient TIJ(*)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: U(*),VOL(*),V(*),A1(*),A2(*),A1X(*),A1Y(*),A2X(*),A2Y(*),
     &        W(*),A3(*),A1Z(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     &        TIJ(*),RO(*),D1(*),D2(*),D3(*)

      REAL, DIMENSION(:), ALLOCATABLE :: RUU, RUV, RUW, RVV, RVW, RWW 
      REAL, DIMENSION(:), ALLOCATABLE :: SIGMA1, SIGMA2, SIGMA3 

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     &           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA,
     &           NSIZE,NGL ! NGL is needed for printing

C ... Three different approaches are available:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID
      NSIZE  = IL*KSTRID 
      
      ALLOCATE(RUU(NSIZE))
      ALLOCATE(RUV(NSIZE))
      ALLOCATE(RUW(NSIZE))
      ALLOCATE(RVV(NSIZE))
      ALLOCATE(RVW(NSIZE))
      ALLOCATE(RWW(NSIZE))
      ALLOCATE(SIGMA1(NSIZE))
      ALLOCATE(SIGMA2(NSIZE))
      ALLOCATE(SIGMA3(NSIZE))

      SIGMA1 = 0.0
      SIGMA2 = 0.0
      SIGMA3 = 0.0

      DO L=1,NSIZE
         RUU(L) = RO(L)*U(L)*U(L)
         RUV(L) = RO(L)*U(L)*V(L)
         RUW(L) = RO(L)*U(L)*W(L) 
         RVV(L) = RO(L)*V(L)*V(L) 
         RVW(L) = RO(L)*V(L)*W(L) 
         RWW(L) = RO(L)*W(L)*W(L) 
      ENDDO

      KALA = 1
      KYLA = KMAX
      JALA = 1
      JYLA = JMAX
      IALA = 1
      IYLA = IMAX

      IF(M == 1) THEN  ! The ghostcells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF

      CALL NABLA_DOT_V(SIGMA1,RUU,RUV,RUW,A1,A2,A3,VOL,D1,D2,D3,
     &     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     &     IALA,IYLA,JALA,JYLA,KALA,KYLA,ISTRID,IL,IN,JN,KN,IDERI)

      CALL NABLA_DOT_V(SIGMA2,RUV,RVV,RVW,A1,A2,A3,VOL,D1,D2,D3,
     &     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     &     IALA,IYLA,JALA,JYLA,KALA,KYLA,ISTRID,IL,IN,JN,KN,IDERI)

      CALL NABLA_DOT_V(SIGMA3,RUW,RVW,RWW,A1,A2,A3,VOL,D1,D2,D3,
     &     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     &     IALA,IYLA,JALA,JYLA,KALA,KYLA,ISTRID,IL,IN,JN,KN,IDERI)

C ... Compute the gradient

      KALA = 1
      KYLA = KMAX
      JALA = 1
      JYLA = JMAX
      IALA = 1
      IYLA = IMAX

      CALL NABLA_DOT_V(TIJ,SIGMA1,SIGMA2,SIGMA3,A1,A2,A3,VOL,D1,D2,D3,
     &     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     &     IALA,IYLA,JALA,JYLA,KALA,KYLA,ISTRID,IL,IN,JN,KN,IDERI)

C ... On the coarse levels extend the values to the ghostcells

      IF(M >= 2) CALL EXTEND(TIJ,IMAX,JMAX,KMAX,IN,JN,KN)

      DEALLOCATE(RUU,RUV,RUW,RVV,RVW,RWW,SIGMA1,SIGMA2,SIGMA3)

      RETURN
      END SUBROUTINE TIJFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NABLA_DOT_V(SIGMA,U,V,W,A1,A2,A3,VOL,D1,D2,D3,
     &           A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,
     &           IALA,IYLA,JALA,JYLA,KALA,KYLA,ISTRID,IL,IN,JN,KN,IDERI)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: SIGMA(*),U(*),V(*),W(*),VOL(*),A1(*),A2(*),A3(*),
     & A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     & D1(*),D2(*),D3(*)

      REAL :: DUDXI,DUDYI,DUDZI,DVDXI,DVDYI,DVDZI,DWDXI,DWDYI,DWDZI,
     &        U1M,V1M,W1M,U1P,V1P,W1P,U2M,V2M,W2M,U2P,V2P,W2P,
     &        U3M,V3M,W3M,U3P,V3P,W3P,PVOL 

      REAL :: W1,W2,W3,R11,R12,R22,R13,R23,R33,DX,DY,DZ,BETA,
     &        ALFA1,ALFA2,ALFA3,R2O3,SQR2

      INTEGER :: I,J,K,KA,JJ,IDERI,IN,JN,KN,L,N,ISTRID,IL,
     &           IALA,IYLA,JALA,JYLA,KALA,KYLA
    
      INTEGER, DIMENSION(6) :: IS


      R2O3    = 2./3.
      SQR2    = SQRT(2.)


      SELECT CASE(IDERI)

      CASE(1)  ! Gauss approach weighted by distances

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2  = 2.0 - W1

                  U1M = A1(L)  *(W2*U(L) + W1*U(L-1))
                  V1M = A1(L)  *(W2*V(L) + W1*V(L-1))
                  W1M = A1(L)  *(W2*W(L) + W1*W(L-1))
                  
                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) )
                  W2  = 2.0 - W1
                  
                  U1P = A1(L+1)*(W2*U(L) + W1*U(L+1))
                  V1P = A1(L+1)*(W2*V(L) + W1*V(L+1))
                  W1P = A1(L+1)*(W2*W(L) + W1*W(L+1))

C ... ETA-DIRECTION

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) )
                  W2  = 2.0 - W1
                  
                  U2M = A2(L)       *(W2*U(L) + W1*U(L-ISTRID))
                  V2M = A2(L)       *(W2*V(L) + W1*V(L-ISTRID))
                  W2M = A2(L)       *(W2*W(L) + W1*W(L-ISTRID))
                  
                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) )
                  W2  = 2.0 - W1
                  
                  U2P = A2(L+ISTRID)*(W2*U(L) + W1*U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(W2*V(L) + W1*V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W2*W(L) + W1*W(L+ISTRID))

C ... ZETA DIRECTION

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) )
                  W2  = 2.0 - W1

                  U3M = A3(L)   *(W2*U(L) + W1*U(L-IL))
                  V3M = A3(L)   *(W2*V(L) + W1*V(L-IL))
                  W3M = A3(L)   *(W2*W(L) + W1*W(L-IL))

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) )
                  W2  = 2.0 - W1

                  U3P = A3(L+IL)*(W2*U(L) + W1*U(L+IL))
                  V3P = A3(L+IL)*(W2*V(L) + W1*V(L+IL))
                  W3P = A3(L+IL)*(W2*W(L) + W1*W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     &                 A3X(L+IL)    *U3P - A3X(L)*U3M +
     &                 A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     &                 A1X(L+1)     *U1P - A1X(L)*U1M)
                  DVDYI = PVOL*(
     &                 A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     &                 A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     &                 A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DWDZI = PVOL*(
     &                 A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     &                 A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     &                 A1Z(L+1)     *W1P - A1Z(L)*W1M)

                  SIGMA(L) = DUDXI + DVDYI + DWDZI

               ENDDO
            ENDDO
         ENDDO

      CASE(2) ! Simple Gaussian derivative calculation

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  U1M = A1(L)  *(U(L) + U(L-1))
                  V1M = A1(L)  *(V(L) + V(L-1))
                  W1M = A1(L)  *(W(L) + W(L-1))
                  U1P = A1(L+1)*(U(L) + U(L+1))
                  V1P = A1(L+1)*(V(L) + V(L+1))
                  W1P = A1(L+1)*(W(L) + W(L+1))

C ... ETA-DIRECTION

                  U2M = A2(L)       *(U(L) + U(L-ISTRID))
                  V2M = A2(L)       *(V(L) + V(L-ISTRID))
                  W2M = A2(L)       *(W(L) + W(L-ISTRID))
                  U2P = A2(L+ISTRID)*(U(L) + U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(V(L) + V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W(L) + W(L+ISTRID))

C ... ZETA DIRECTION

                  U3M = A3(L)   *(U(L) + U(L-IL))
                  V3M = A3(L)   *(V(L) + V(L-IL))
                  W3M = A3(L)   *(W(L) + W(L-IL))
                  U3P = A3(L+IL)*(U(L) + U(L+IL))
                  V3P = A3(L+IL)*(V(L) + V(L+IL))
                  W3P = A3(L+IL)*(W(L) + W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     &                 A3X(L+IL)    *U3P - A3X(L)*U3M +
     &                 A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     &                 A1X(L+1)     *U1P - A1X(L)*U1M)
                  DVDYI = PVOL*(
     &                 A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     &                 A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     &                 A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DWDZI = PVOL*(
     &                 A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     &                 A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     &                 A1Z(L+1)     *W1P - A1Z(L)*W1M)

                  SIGMA(L) = DUDXI + DVDYI + DWDZI

               ENDDO
            ENDDO
         ENDDO

      CASE(3) ! Least-squares approach

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

                  IS(1) = -1
                  IS(2) = -ISTRID
                  IS(3) = -IL
                  IS(4) =  1
                  IS(5) =  ISTRID
                  IS(6) =  IL

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO
                  
                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))
                  
                  BETA = (R12*R23 - R13*R22) / (R11*R22)
                  
                  DUDXI = 0.0
                  DVDYI = 0.0
                  DWDZI = 0.0
                  
                  DO N = 1,6
                     DX    = XC(L + IS(N)) - XC(L)
                     DY    = YC(L + IS(N)) - YC(L)
                     DZ    = ZC(L + IS(N)) - ZC(L)
                     ALFA1 =  DX / R11**2
                     ALFA2 = (DY - R12*DX/R11) / R22**2
                     ALFA3 = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1    = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2    = ALFA2 - R23*ALFA3/R22
                     W3    = ALFA3
                     DUDXI = DUDXI + W1*(U(L+IS(N))-U(L))
                     DVDYI = DVDYI + W2*(V(L+IS(N))-V(L))
                     DWDZI = DWDZI + W3*(W(L+IS(N))-W(L))
                  ENDDO

                  SIGMA(L) = DUDXI + DVDYI + DWDZI

               ENDDO
            ENDDO
         ENDDO
                  
      END SELECT

      RETURN
      END SUBROUTINE NABLA_DOT_V
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXTEND(OHMI,IMAX,JMAX,KMAX,IN,JN,KN)

C ... This subroutine extends the values into the first ghost cells
C ... using a firts-order extrapolation. Could be replaced by a 2nd 
C ... order method.

      IMPLICIT NONE

      REAL    :: OHMI(*)
      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,KSTRID,IL
      INTEGER :: I,J,K,L,M,KA,JJ

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID


C ... K=0 wall

      DO K=0,0
         KA = (KN+K-1)*IL
         DO J=1,JMAX
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=1,IMAX
               L = JJ + I
               M = L + IL
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   


C ... K=KMAX wall

      DO K=KMAX+1,KMAX+1
         KA = (KN+K-1)*IL
         DO J=1,JMAX
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=1,IMAX
               L = JJ + I
               M = L - IL
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   


C ... J=0 wall

      DO K=1,KMAX
         KA = (KN+K-1)*IL
         DO J=0,0
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=1,IMAX
               L = JJ + I
               M = L + ISTRID
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   


C ... J=JMAX+1 wall

      DO K=1,KMAX
         KA = (KN+K-1)*IL
         DO J=JMAX+1,JMAX+1
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=1,IMAX
               L = JJ + I
               M = L - ISTRID
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   


C ... I=0 wall

      DO K=1,KMAX
         KA = (KN+K-1)*IL
         DO J=1,JMAX
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=0,0
               L = JJ + I
               M = L + 1
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   


C ... I=IMAX+1 wall

      DO K=1,KMAX
         KA = (KN+K-1)*IL
         DO J=1,JMAX
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I=IMAX+1,IMAX+1
               L = JJ + I
               M = L - 1
               OHMI(L) = OHMI(M)
            ENDDO
         ENDDO
      ENDDO   

      RETURN
      END SUBROUTINE EXTEND
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TEMDER(TEMP,DTDX,DTDY,DTDZ,A1,A2,A3,VOL,D1,D2,D3,
     + A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,IMAX,JMAX,KMAX,
     + IN,JN,KN,IDERI,M)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: TEMP(*),DTDX(*),DTDY(*),DTDZ(*),A1(*),A2(*),A3(*),VOL(*),
     2        A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),
     3        A3Z(*),D1(*),D2(*),D3(*)

      REAL :: T1MP,T1PP,T2MP,T2PP,T3MP,T3PP,PVOL

      REAL :: W1,W2,W3,R11,R12,R22,R13,R23,R33,DX,DY,DZ,BETA,
     2        ALFA1,ALFA2,ALFA3

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     2           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA

      INTEGER, DIMENSION(6) :: IS
     
C ... Calculate a gradient for a single variable

C ... Three different approaches are available:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach
 
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID

      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghost cells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF


      SELECT CASE(IDERI)

      CASE(1)  ! Gauss approach weighted by distances

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  W1   = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2   = 2.0 - W1
                  T1MP = W2*TEMP(L) + W1*TEMP(L-1)
                  W1   = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) ) 
                  W2   = 2.0 - W1
                  T1PP = W2*TEMP(L) + W1*TEMP(L+1)

C ... ETA-DIRECTION

                  W1   = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) ) 
                  W2   = 2.0 - W1
                  T2MP = W2*TEMP(L) + W1*TEMP(L-ISTRID)
                  W1   = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) ) 
                  W2   = 2.0 - W1
                  T2PP = W2*TEMP(L) + W1*TEMP(L+ISTRID)

C ... ZETA DIRECTION

                  W1   = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) ) 
                  W2   = 2.0 - W1
                  T3MP = W2*TEMP(L) + W1*TEMP(L-IL)
                  W1   = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) ) 
                  W2   = 2.0 - W1
                  T3PP = W2*TEMP(L) + W1*TEMP(L+IL)

                  PVOL = .5/VOL(L)

                  DTDX(L) = PVOL*(
     1             A3(L+IL)    *A3X(L+IL)    *T3PP - A3(L)*A3X(L)*T3MP +
     2             A2(L+ISTRID)*A2X(L+ISTRID)*T2PP - A2(L)*A2X(L)*T2MP +
     3             A1(L+1)     *A1X(L+1)     *T1PP - A1(L)*A1X(L)*T1MP)
                  DTDY(L) = PVOL*(
     1             A3(L+IL)    *A3Y(L+IL)    *T3PP - A3(L)*A3Y(L)*T3MP +
     2             A2(L+ISTRID)*A2Y(L+ISTRID)*T2PP - A2(L)*A2Y(L)*T2MP +
     3             A1(L+1)     *A1Y(L+1)     *T1PP - A1(L)*A1Y(L)*T1MP)
                  DTDZ(L)     = PVOL*(
     1             A3(L+IL)    *A3Z(L+IL)    *T3PP - A3(L)*A3Z(L)*T3MP +
     2             A2(L+ISTRID)*A2Z(L+ISTRID)*T2PP - A2(L)*A2Z(L)*T2MP +
     3             A1(L+1)     *A1Z(L+1)     *T1PP - A1(L)*A1Z(L)*T1MP)

               ENDDO
            ENDDO
         ENDDO

C ******************************************************************

      CASE(2) ! Simple Gaussian derivative calculation

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  T1MP = TEMP(L) + TEMP(L-1)
                  T1PP = TEMP(L) + TEMP(L+1)

C ... ETA-DIRECTION

                  T2MP = TEMP(L) + TEMP(L-ISTRID)
                  T2PP = TEMP(L) + TEMP(L+ISTRID)

C ... ZETA DIRECTION

                  T3MP = TEMP(L) + TEMP(L-IL)
                  T3PP = TEMP(L) + TEMP(L+IL)

                  PVOL = .5/VOL(L)

                  DTDX(L) = PVOL*(
     1             A3(L+IL)    *A3X(L+IL)    *T3PP - A3(L)*A3X(L)*T3MP +
     2             A2(L+ISTRID)*A2X(L+ISTRID)*T2PP - A2(L)*A2X(L)*T2MP +
     3             A1(L+1)     *A1X(L+1)     *T1PP - A1(L)*A1X(L)*T1MP)
                  DTDY(L) = PVOL*(
     1             A3(L+IL)    *A3Y(L+IL)    *T3PP - A3(L)*A3Y(L)*T3MP +
     2             A2(L+ISTRID)*A2Y(L+ISTRID)*T2PP - A2(L)*A2Y(L)*T2MP +
     3             A1(L+1)     *A1Y(L+1)     *T1PP - A1(L)*A1Y(L)*T1MP)
                  DTDZ(L) = PVOL*(
     1             A3(L+IL)    *A3Z(L+IL)    *T3PP - A3(L)*A3Z(L)*T3MP +
     2             A2(L+ISTRID)*A2Z(L+ISTRID)*T2PP - A2(L)*A2Z(L)*T2MP +
     3             A1(L+1)     *A1Z(L+1)     *T1PP - A1(L)*A1Z(L)*T1MP)

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(3) ! Least Squares method

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001
         
         IS(1) = -1
         IS(2) = -ISTRID
         IS(3) = -IL
         IS(4) =  1
         IS(5) =  ISTRID
         IS(6) =  IL


         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO

                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))

                  BETA = (R12*R23 - R13*R22) / (R11*R22)

                  DTDX(L) = 0.0
                  DTDY(L) = 0.0
                  DTDZ(L) = 0.0

                  DO N = 1,6
                     DX      = XC(L + IS(N)) - XC(L)
                     DY      = YC(L + IS(N)) - YC(L)
                     DZ      = ZC(L + IS(N)) - ZC(L)
                     ALFA1   =  DX / R11**2
                     ALFA2   = (DY - R12*DX/R11) / R22**2
                     ALFA3   = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1      = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2      = ALFA2 - R23*ALFA3/R22
                     W3      = ALFA3
                     DTDX(L) = DTDX(L) + W1*(TEMP(L+IS(N))-TEMP(L))
                     DTDY(L) = DTDY(L) + W2*(TEMP(L+IS(N))-TEMP(L))
                     DTDZ(L) = DTDZ(L) + W3*(TEMP(L+IS(N))-TEMP(L))

                  ENDDO

               ENDDO
            ENDDO
         ENDDO

      END SELECT


C ... On the coarse levels extend the values to the ghostcells
C ... This affects only the convergence rate and is not very important

      IF(M >= 2) THEN
         CALL EXTEND(DTDX,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDY,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDZ,IMAX,JMAX,KMAX,IN,JN,KN)
      ENDIF      

      RETURN
      END SUBROUTINE TEMDER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TEMDEP(TEMP,DTDX,DTDY,DTDZ,A1,A2,A3,VOL,D1,D2,D3,
     + A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,IMAX,JMAX,KMAX,
     + IN,JN,KN,IDERI,M)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)
        
      real :: sxi,syi,szi,slen,alfa,frki,frkj,dkdxi,dkdyi,dkdzi,suhde

      REAL :: TEMP(*),DTDX(*),DTDY(*),DTDZ(*),A1(*),A2(*),A3(*),VOL(*),
     2        A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),
     3        A3Z(*),D1(*),D2(*),D3(*)

      REAL :: T1MP,T1PP,T2MP,T2PP,T3MP,T3PP,PVOL

      REAL :: W1,W2,W3,R11,R12,R22,R13,R23,R33,DX,DY,DZ,BETA,
     2        ALFA1,ALFA2,ALFA3

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     2           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA

      INTEGER, DIMENSION(6) :: IS

C ...  Temporary routine with printing option
     
C ... Calculate a gradient for a single variable

C ... Three different approaches are available:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach
 
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID

      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghost cells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF


      SELECT CASE(IDERI)

      CASE(1)  ! Gauss approach weighted by distances

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  W1   = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2   = 2.0 - W1
                  T1MP = W2*TEMP(L) + W1*TEMP(L-1)
                  W1   = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) ) 
                  W2   = 2.0 - W1
                  T1PP = W2*TEMP(L) + W1*TEMP(L+1)

C ... ETA-DIRECTION

                  W1   = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) ) 
                  W2   = 2.0 - W1
                  T2MP = W2*TEMP(L) + W1*TEMP(L-ISTRID)
                  W1   = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) ) 
                  W2   = 2.0 - W1
                  T2PP = W2*TEMP(L) + W1*TEMP(L+ISTRID)

C ... ZETA DIRECTION

                  W1   = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) ) 
                  W2   = 2.0 - W1
                  T3MP = W2*TEMP(L) + W1*TEMP(L-IL)
                  W1   = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) ) 
                  W2   = 2.0 - W1
                  T3PP = W2*TEMP(L) + W1*TEMP(L+IL)

                  PVOL = .5/VOL(L)

                  DTDX(L) = PVOL*(
     1             A3(L+IL)    *A3X(L+IL)    *T3PP - A3(L)*A3X(L)*T3MP +
     2             A2(L+ISTRID)*A2X(L+ISTRID)*T2PP - A2(L)*A2X(L)*T2MP +
     3             A1(L+1)     *A1X(L+1)     *T1PP - A1(L)*A1X(L)*T1MP)
                  DTDY(L) = PVOL*(
     1             A3(L+IL)    *A3Y(L+IL)    *T3PP - A3(L)*A3Y(L)*T3MP +
     2             A2(L+ISTRID)*A2Y(L+ISTRID)*T2PP - A2(L)*A2Y(L)*T2MP +
     3             A1(L+1)     *A1Y(L+1)     *T1PP - A1(L)*A1Y(L)*T1MP)
                  DTDZ(L)     = PVOL*(
     1             A3(L+IL)    *A3Z(L+IL)    *T3PP - A3(L)*A3Z(L)*T3MP +
     2             A2(L+ISTRID)*A2Z(L+ISTRID)*T2PP - A2(L)*A2Z(L)*T2MP +
     3             A1(L+1)     *A1Z(L+1)     *T1PP - A1(L)*A1Z(L)*T1MP)

               ENDDO
            ENDDO
         ENDDO

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I
            SXI    = XC(L) - XC(L-ISTRID)
            SYI    = YC(L) - YC(L-ISTRID)
            SZI    = ZC(L) - ZC(L-ISTRID)
            SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
            SXI    = SXI/SLEN
            SYI    = SYI/SLEN
            SZI    = SZI/SLEN
            ALFA   = SXI*A2X(L) + SYI*A2Y(L) + SZI*A2Z(L)
            frki =  -(TEMP(L) - TEMP(L-ISTRID))/SLEN 
            DKDXI   = 0.5*(DTDX(L) + DTDX(L-ISTRID))
            DKDYI   = 0.5*(DTDY(L) + DTDY(L-ISTRID))
            DKDZI   = 0.5*(DTDZ(L) + DTDZ(L-ISTRID))

            frkj = -(DKDXI*SXI + DKDYI*SYI + DKDZI*SZI)
            suhde = frki/(frkj + 1.e-21)
c            write(201,9876) i,j,k,l,frki,frkj,suhde
c            write(202,9876) i,j,k,l,dkdxi,dkdyi,dkdzi
c            write(203,9876) i,j,k,l,sxi,syi,szi,slen
c            write(203,9876) i,j,k,l,dtdy(l),dtdy(l-istrid),
c     2                   dtdz(l),dtdz(l-istrid)
9876        FORMAT(4I6,4E12.5)

               ENDDO
            ENDDO
         ENDDO

C ******************************************************************

      CASE(2) ! Simple Gaussian derivative calculation

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  T1MP = TEMP(L) + TEMP(L-1)
                  T1PP = TEMP(L) + TEMP(L+1)

C ... ETA-DIRECTION

                  T2MP = TEMP(L) + TEMP(L-ISTRID)
                  T2PP = TEMP(L) + TEMP(L+ISTRID)

C ... ZETA DIRECTION

                  T3MP = TEMP(L) + TEMP(L-IL)
                  T3PP = TEMP(L) + TEMP(L+IL)

                  PVOL = .5/VOL(L)

                  DTDX(L) = PVOL*(
     1             A3(L+IL)    *A3X(L+IL)    *T3PP - A3(L)*A3X(L)*T3MP +
     2             A2(L+ISTRID)*A2X(L+ISTRID)*T2PP - A2(L)*A2X(L)*T2MP +
     3             A1(L+1)     *A1X(L+1)     *T1PP - A1(L)*A1X(L)*T1MP)
                  DTDY(L) = PVOL*(
     1             A3(L+IL)    *A3Y(L+IL)    *T3PP - A3(L)*A3Y(L)*T3MP +
     2             A2(L+ISTRID)*A2Y(L+ISTRID)*T2PP - A2(L)*A2Y(L)*T2MP +
     3             A1(L+1)     *A1Y(L+1)     *T1PP - A1(L)*A1Y(L)*T1MP)
                  DTDZ(L) = PVOL*(
     1             A3(L+IL)    *A3Z(L+IL)    *T3PP - A3(L)*A3Z(L)*T3MP +
     2             A2(L+ISTRID)*A2Z(L+ISTRID)*T2PP - A2(L)*A2Z(L)*T2MP +
     3             A1(L+1)     *A1Z(L+1)     *T1PP - A1(L)*A1Z(L)*T1MP)

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(3) ! Least Squares method

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001
         
         IS(1) = -1
         IS(2) = -ISTRID
         IS(3) = -IL
         IS(4) =  1
         IS(5) =  ISTRID
         IS(6) =  IL


         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO

                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))

                  BETA = (R12*R23 - R13*R22) / (R11*R22)

                  DTDX(L) = 0.0
                  DTDY(L) = 0.0
                  DTDZ(L) = 0.0

                  DO N = 1,6
                     DX      = XC(L + IS(N)) - XC(L)
                     DY      = YC(L + IS(N)) - YC(L)
                     DZ      = ZC(L + IS(N)) - ZC(L)
                     ALFA1   =  DX / R11**2
                     ALFA2   = (DY - R12*DX/R11) / R22**2
                     ALFA3   = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1      = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2      = ALFA2 - R23*ALFA3/R22
                     W3      = ALFA3
                     DTDX(L) = DTDX(L) + W1*(TEMP(L+IS(N))-TEMP(L))
                     DTDY(L) = DTDY(L) + W2*(TEMP(L+IS(N))-TEMP(L))
                     DTDZ(L) = DTDZ(L) + W3*(TEMP(L+IS(N))-TEMP(L))
                  ENDDO

               ENDDO
            ENDDO
         ENDDO

      END SELECT


C ... On the coarse levels extend the values to the ghostcells
C ... This affects only the convergence rate and is not very important

      IF(M >= 2) THEN
         CALL EXTEND(DTDX,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDY,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDZ,IMAX,JMAX,KMAX,IN,JN,KN)
      ENDIF      

      RETURN
      END SUBROUTINE TEMDEP


      SUBROUTINE PDER(TEMP,DTDX,DTDY,DTDZ,A1,A2,A3,VOL,D1,D2,D3,
     + A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,XC,YC,ZC,IMAX,JMAX,KMAX,
     + IN,JN,KN,IDERI,M,INTERI,INTERJ,INTERK,RKSI,INCHIML,NGL)

      IMPLICIT NONE


      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     2           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA

      INTEGER :: INTER,INTERI,INTERJ,INTERK,INTE2,IT0,IT1,IT2,NGL

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: TEMP(*),DTDX(*),DTDY(*),DTDZ(*),A1(*),A2(*),A3(*),VOL(*),
     2        A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),
     3        A3Z(*),D1(*),D2(*),D3(*),RKSI(*)

      REAL :: T1MP,T1PP,T2MP,T2PP,T3MP,T3PP,PVOL,TL,TR,
     2        RK1,RK2,R1,R2,EPS,DIFF0,DIFF1,DIFF2,DIFF12,PHIBI,PHIFI

      LOGICAL :: INCHIML

      REAL RKK1(5)
      DATA RKK1/0.,1.,1.3333333,2.,1.5/

C ... INTER     Method
C ... 1         SECOND-ORDER UPWIND
C ... 2         SECOND-ORDER UPWIND BIASED
C ... 3         THIRD-ORDER UPWIND BIASED
C ... 4         CENTRAL
C ... 5         QUICK

      EPS     = 1.E-10
 
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID

      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghost cells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF

      DO K = KALA,KYLA
         KA = (KN+K-1)*IL
         DO J = JALA,JYLA
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I = IALA,IYLA
               L = JJ + I

C ... XI-DIRECTION

      IT0     = L - 2
      IT1     = L - 1
      IT2     = L + 1
      INTER   = IABS(INTERI)
      INTE2   = MIN0(5,MAX0(1,INTER))
      RK1     = RKK1(INTE2)
      RK2     = 2.-RK1
      R1      = .25*RK1
      R2      = .25*RK2

      IF (INTERI <= 0 .AND. INTERI >= -5) THEN ! WITHOUT LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            TL      = TEMP(L)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1) + (DIFF1*R1 + DIFF0*R2)
            T1MP    = TL + TR
c       write(779,*)'t1mp',i,j,k,t1mp,tl,tr,temp(l),temp(it1),r1,r2,inte2

            DIFF0   = TEMP(IT1+1) - TEMP(IT0+1)
            DIFF1   = TEMP(L+1)   - TEMP(IT1+1)
            DIFF2   = TEMP(IT2+1) - TEMP(L+1)
            TL      = TEMP(L+1)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1+1) + (DIFF1*R1 + DIFF0*R2)
            T1PP    = TL + TR
c       write(779,*) 't1pp',i,j,k,t1pp,tl,tr,temp(l+1),temp(it1+1)

c       ENDIF

      ELSE IF(INTERI <= 4) THEN ! FLUX LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL     = TEMP(L)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR     = TEMP(IT1) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            T1MP   = TL + TR
            DIFF0   = TEMP(IT1+1) - TEMP(IT0+1)
            DIFF1   = TEMP(L+1)   - TEMP(IT1+1)
            DIFF2   = TEMP(IT2+1) - TEMP(L+1)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL     = TEMP(L+1)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR     = TEMP(IT1+1) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            T1PP   = TL + TR

      ENDIF

C ... ETA-DIRECTION

      IT0     = L - 2*ISTRID
      IT1     = L - ISTRID
      IT2     = L + ISTRID
      INTER   = IABS(INTERJ)
      INTE2   = MIN0(5,MAX0(1,INTER))
      RK1     = RKK1(INTE2)
      RK2     = 2.-RK1
      R1      = .25*RK1
      R2      = .25*RK2

      IF (INTERJ <= 0 .AND. INTERJ >= -5) THEN ! WITHOUT LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            TL      = TEMP(L)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1) + (DIFF1*R1 + DIFF0*R2)
            T2MP    = TL + TR
            DIFF0   = TEMP(IT1+ISTRID) - TEMP(IT0+ISTRID)
            DIFF1   = TEMP(L+ISTRID)   - TEMP(IT1+ISTRID)
            DIFF2   = TEMP(IT2+ISTRID) - TEMP(L+ISTRID)
            TL      = TEMP(L+ISTRID)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1+ISTRID) + (DIFF1*R1 + DIFF0*R2)
c      if(ngl == 1 .and. j == 16) TL = TR

            T2PP    = TL + TR
c       ENDIF

      ELSE IF(INTERJ <= 4) THEN ! FLUX LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL       = TEMP(L)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR     = TEMP(IT1) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            T2MP    = TL + TR
            DIFF0   = TEMP(IT1+ISTRID) - TEMP(IT0+ISTRID)
            DIFF1   = TEMP(L+ISTRID)   - TEMP(IT1+ISTRID)
            DIFF2   = TEMP(IT2+ISTRID) - TEMP(L+ISTRID)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL       = TEMP(L+ISTRID)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR       = TEMP(IT1+ISTRID) + PHIBI*(DIFF1*R1 + DIFF0*R2)

            T2PP     = TL + TR
      ENDIF

C ... ZETA-DIRECTION

      IT0     = L - 2*IL
      IT1     = L - IL
      IT2     = L + IL
      INTER   = IABS(INTERK)
      INTE2   = MIN0(5,MAX0(1,INTER))
      RK1     = RKK1(INTE2)
      RK2     = 2.-RK1
      R1      = .25*RK1
      R2      = .25*RK2

      IF (INTERK <= 0 .AND. INTERK >= -5) THEN ! WITHOUT LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            TL      = TEMP(L)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1) + (DIFF1*R1 + DIFF0*R2)
            T3MP    = TL + TR
            DIFF0   = TEMP(IT1+IL) - TEMP(IT0+IL)
            DIFF1   = TEMP(L+IL)   - TEMP(IT1+IL)
            DIFF2   = TEMP(IT2+IL) - TEMP(L+IL)
            TL      = TEMP(L+IL)   - (DIFF1*R1 + DIFF2*R2)
            TR      = TEMP(IT1+IL) + (DIFF1*R1 + DIFF0*R2)
            T3PP    = TL + TR
c       ENDIF

      ELSE IF(INTERK <= 4) THEN ! FLUX LIMITATION

c         IF((RKSI(I) /= 0 .AND. RKSI(IT1) == 0.  .OR. ! First-order
c     +       RKSI(I) == 0 .AND. RKSI(IT1) /= 0.) .AND. INCHIML) THEN
c            ROIP(IG)= RO(I)  
c            ROIM(IG)= RO(IT1)
c         ELSE
            DIFF0   = TEMP(IT1) - TEMP(IT0)
            DIFF1   = TEMP(L)   - TEMP(IT1)
            DIFF2   = TEMP(IT2) - TEMP(L)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL       = TEMP(L)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR       = TEMP(IT1) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            T3MP     = TL + TR
            DIFF0   = TEMP(IT1+IL) - TEMP(IT0+IL)
            DIFF1   = TEMP(L+IL)   - TEMP(IT1+IL)
            DIFF2   = TEMP(IT2+IL) - TEMP(L+IL)
            DIFF12  = DIFF1**2 + EPS
C ... LIMITERS (VAN ALBADA)
            PHIBI    = (DIFF12 + DIFF1*DIFF0)/(DIFF12 + DIFF0**2)
            PHIFI    = (DIFF12 + DIFF1*DIFF2)/(DIFF12 + DIFF2**2)
C ... EXTRAPOLATING
            TL       = TEMP(L+IL)   - PHIFI*(DIFF1*R1 + DIFF2*R2)
            TR       = TEMP(IT1+IL) + PHIBI*(DIFF1*R1 + DIFF0*R2)
            T3PP     = TL + TR
      ENDIF

                  PVOL = .5/VOL(L)

                  DTDX(L) = PVOL*(
     1             A3(L+IL)    *A3X(L+IL)    *T3PP - A3(L)*A3X(L)*T3MP +
     2             A2(L+ISTRID)*A2X(L+ISTRID)*T2PP - A2(L)*A2X(L)*T2MP +
     3             A1(L+1)     *A1X(L+1)     *T1PP - A1(L)*A1X(L)*T1MP)
                  DTDY(L) = PVOL*(
     1             A3(L+IL)    *A3Y(L+IL)    *T3PP - A3(L)*A3Y(L)*T3MP +
     2             A2(L+ISTRID)*A2Y(L+ISTRID)*T2PP - A2(L)*A2Y(L)*T2MP +
     3             A1(L+1)     *A1Y(L+1)     *T1PP - A1(L)*A1Y(L)*T1MP)
                  DTDZ(L) = PVOL*(
     1             A3(L+IL)    *A3Z(L+IL)    *T3PP - A3(L)*A3Z(L)*T3MP +
     2             A2(L+ISTRID)*A2Z(L+ISTRID)*T2PP - A2(L)*A2Z(L)*T2MP +
     3             A1(L+1)     *A1Z(L+1)     *T1PP - A1(L)*A1Z(L)*T1MP)
c      if(ngl == 1 .and.j >= 15) then
c      write(7778,*) i,j,k,vol(l),dtdx(l),dtdy(l),dtdz(l)
c      write(7778,*) t1mp,t1pp,t2mp,t2pp,t3mp,t3pp
c      write(7778,*) temp(l),temp(l-1),temp(l-istrid),temp(l-il)
c      write(7778,*) temp(l),temp(l+1),temp(l+istrid),temp(l+il)
c      endif
c       if(abs(dtdx(l)) > 1.e8) write(7000,*) ngl,i,j,k,dtdx(l)
               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

C ... On the coarse levels extend the values to the ghostcells
C ... This affects only the convergence rate and is not very important

      IF(M >= 2) THEN
         CALL EXTEND(DTDX,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDY,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(DTDZ,IMAX,JMAX,KMAX,IN,JN,KN)
      ENDIF      

      RETURN
      END SUBROUTINE PDER


      SUBROUTINE NABLA_CROSS_V(NABLAXV,NABLAYV,NABLAZV,D1,D2,D3,
     + VOL,A1,A2,A3,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,U,V,W,
     + A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      REAL :: EPS,W1,W2,U3P,V3P,W3P
      INTEGER :: L,ISTR,JSTR,IL,I,J,K,II,JJ,IS
      REAL :: NABLAXV(*),NABLAYV(*),NABLAZV(*),!OHMI(*),
     + D1(*),D2(*),D3(*),VOL(*),
     + A1(*),A2(*),A3(*),U(*),V(*),W(*),A1X(*),A1Y(*),
     + A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)
    
C
C ... This should calculate a rotor of a vector
C ... directly from Sum_j (n_j x V_j)S_j
C

      NABLAXV(1:NTOT) = 0.; NABLAYV(1:NTOT) = 0.; NABLAZV(1:NTOT) = 0.
     
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
      IS   = ISTR ! Antakee mersu
      EPS  = 1.E-20

C ... K-direction face by face into two cells simultaneously

      DO K = 0,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

         L   = I + II
         W1  = D3(L) / ( D3(L) + D3(L+IL) )
         W2  = 1.0 - W1
c         w1=0.5; w2 = 0.5

         U3P = A3(L+IL)*(W2*U(L) + W1*U(L+IL))
         V3P = A3(L+IL)*(W2*V(L) + W1*V(L+IL))
         W3P = A3(L+IL)*(W2*W(L) + W1*W(L+IL))

         NABLAXV(L)    = NABLAXV(L)    + A3Y(L+IL)*W3P - A3Z(L+IL)*V3P
         NABLAXV(L+IL) = NABLAXV(L+IL) - A3Y(L+IL)*W3P + A3Z(L+IL)*V3P

         NABLAYV(L)    = NABLAYV(L)    + A3Z(L+IL)*U3P - A3X(L+IL)*W3P
         NABLAYV(L+IL) = NABLAYV(L+IL) - A3Z(L+IL)*U3P + A3X(L+IL)*W3P

         NABLAZV(L)    = NABLAZV(L)    + A3X(L+IL)*V3P - A3Y(L+IL)*U3P
         NABLAZV(L+IL) = NABLAZV(L+IL) - A3X(L+IL)*V3P + A3Y(L+IL)*U3P

      ENDDO ; ENDDO ; ENDDO

C ... J-direction similarly

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

         L   = I + II
         W1  = D2(L) / ( D2(L) + D2(L+IS) ) 
         W2  = 1.0 - W1
c         w1=0.5; w2 = 0.5

         U3P = A2(L+IS)*(W2*U(L) + W1*U(L+IS))
         V3P = A2(L+IS)*(W2*V(L) + W1*V(L+IS))
         W3P = A2(L+IS)*(W2*W(L) + W1*W(L+IS))

         NABLAXV(L)    = NABLAXV(L)    + A2Y(L+IS)*W3P - A2Z(L+IS)*V3P
         NABLAXV(L+IS) = NABLAXV(L+IS) - A2Y(L+IS)*W3P + A2Z(L+IS)*V3P

         NABLAYV(L)    = NABLAYV(L)    + A2Z(L+IS)*U3P - A2X(L+IS)*W3P
         NABLAYV(L+IS) = NABLAYV(L+IS) - A2Z(L+IS)*U3P + A2X(L+IS)*W3P

         NABLAZV(L)    = NABLAZV(L)    + A2X(L+IS)*V3P - A2Y(L+IS)*U3P
         NABLAZV(L+IS) = NABLAZV(L+IS) - A2X(L+IS)*V3P + A2Y(L+IS)*U3P

      ENDDO ; ENDDO ; ENDDO

C ... Finally I-direction

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX

         L   = I + II
         W1  = D1(L) / ( D1(L) + D1(L+1) ) 
         W2  = 1.0 - W1
c         w1=0.5; w2 = 0.5

         U3P = A1(L+1)*(W2*U(L) + W1*U(L+1))
         V3P = A1(L+1)*(W2*V(L) + W1*V(L+1))
         W3P = A1(L+1)*(W2*W(L) + W1*W(L+1))

         NABLAXV(L)   = NABLAXV(L)   + A1Y(L+1)*W3P - A1Z(L+1)*V3P
         NABLAXV(L+1) = NABLAXV(L+1) - A1Y(L+1)*W3P + A1Z(L+1)*V3P

         NABLAYV(L)   = NABLAYV(L)   + A1Z(L+1)*U3P - A1X(L+1)*W3P
         NABLAYV(L+1) = NABLAYV(L+1) - A1Z(L+1)*U3P + A1X(L+1)*W3P

         NABLAZV(L)   = NABLAZV(L)   + A1X(L+1)*V3P - A1Y(L+1)*U3P
         NABLAZV(L+1) = NABLAZV(L+1) - A1X(L+1)*V3P + A1Y(L+1)*U3P

      ENDDO ; ENDDO ; ENDDO

C ... Finally divide by the cell volume

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

         L   = I + II
         NABLAXV(L) = NABLAXV(L)/(VOL(L)+EPS)
         NABLAYV(L) = NABLAYV(L)/(VOL(L)+EPS)
         NABLAZV(L) = NABLAZV(L)/(VOL(L)+EPS)

C ... These are calculated eelsewhere (VORFUN, DERFUN, TURFUN)
C         OHMI(L) = NABLAXV(L)**2 + NABLAYV(L)**2 + NABLAZV(L)**2
C         OHMI(L) = SQRT(ohmi(L))

      ENDDO; ENDDO; ENDDO

C ... These are probably not necessary

C         CALL EXTEND(OHMI,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(NABLAXV,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(NABLAYV,IMAX,JMAX,KMAX,IN,JN,KN)
         CALL EXTEND(NABLAZV,IMAX,JMAX,KMAX,IN,JN,KN)

      RETURN
      END SUBROUTINE NABLA_CROSS_V


      SUBROUTINE GATSKI(OHMI,W12,W13,W23,S11,S12,S13,S22,S23,S33,
     +     UU,UV,UW,VV,VW,WW,PTUR,U,V,W,RK,REPS,DDEPS,EPS2,VIST,VIS,
     +     OMEGA,IMAX,JMAX,KMAX,IN,JN,KN,ITURB,M)

      REAL :: OHMI(*),UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),
     +        W12(*),W13(*),W23(*),PTUR(*),RK(*),REPS(*),EPS2(*),VIS(*),
     +        S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),DDEPS(*),VIST(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      R2O3    = 2./3.

C ... Constants for a Speziale's model (ITURB=10) AIAA Journal, Vol.33, No. 10.

      ALF1    = .375
      ALF2    = .116
      ALF3    = .108

      SQR2    = 1./SQRT(2.)
      EPS     = 1.E-20

C ... LRR mallin vakiot

      C1 = 3.
      C2 = 0.8
      C3 = 1.75
      C4 = 1.31
      C5 = 2.

C ... SSG mallin vakiot

      C1 = 6.8
      C2 = 0.36
      C3 = 1.25
      C4 = 0.4
      C5 = 1.88

      GG      = 1./(C1/2. + C5 - 1.)
      ALF1    = (4./3.-C2)*GG/2.
      ALF2    = (2. - C3)**2*GG**2/4.
      ALF3    = (2. - C4)**2*GG**2/4.
      ALF4    = (2. - C4)*GG/2.
      ALF5    = (2. - C3)*GG

      DO 5000 K = 1,KMAX
      IA  = (KN+K-1)*IL + JN*ISTRID

      IF(ITURB == 10) THEN
      DO 4000 JJ = 1,JMAX
      DO 4000 II = 1,IMAX
         L       = IA + (JJ-1)*ISTRID + (II+IN-1) + 1

         DIV  = (S11(L)+S22(L)+S33(L))/3.

         S11L = S11(L)-DIV
         S22L = S22(L)-DIV
         S33L = S33(L)-DIV


         W13L = W13(L)
         W12L = W12(L) + 2./(C4-2.)*OMEGA ! vasta 26.1.2000
         W23L = W23(L) 
         W31L = -W13L
         W21L = -W12L
         W32L = -W23L

         SL  = S11L*S11L     + S12(L)*S12(L) + S13(L)*S13(L)
     +       + S12(L)*S12(L) + S22L*S22L     + S23(L)*S23(L)
     +       + S13(L)*S13(L) + S23(L)*S23(L) + S33L*S33L
         SL   = SQRT(SL)


         EPSI    = REPS(L) + DDEPS(L)
         RKPE    = RK(L)/(EPSI+EPS)
         RKPE2   = RK(L)**2/(EPSI+EPS)
         ETA     = ALF2*(SL*RKPE)**2
         XSI     = ALF3*(SQR2*OHMI(L)*RKPE)**2
c         XSI     = ALF3*(SL*RKPE)**2
         CMUT    = 3.*(1. + ETA)/(3. + ETA + 6.*XSI*ETA + 6.*XSI)*ALF1
C ... MODIFICATIONS....
c         CMUT    = ALF1*(3.*(1. + ETA) +   .2*(ETA**3 + XSI**3))/ 
c     +        (3. + ETA +6.*XSI*ETA + 6.*XSI + ETA**3 + XSI**3)
         VMUT    = CMUT*RKPE2
         call ijkpai(l,imax,jmax,kmax,mmm,nnn,lll)
c         if(mmm ==  40 .and. lll == 1 .and.jmax == 32) then
c            write(77,313) nnn,vmut/vis(l),cmut,reps(l),ddeps(l)
c 313        format(I5,6E14.5)
c            write(88,*) nnn,rkpe2,eta,xsi
c         endif
c        SKEP    = SQRT(2.0)*SL*RKPE
c        OKEP    = SQRT(0.5)*OHMI(L)*RKPE
c        CMUT2    = 1.0/(9.0+1.4*SKEP)         
c        GG      = 5.6/(20.0+CMUT2*(SKEP**2+OKEP**2)) !2.5,2.6,3.0/9.0
c        ALF4    = (2. - C4)*GG/2.
c        ALF5    = (2. - C3)*GG



         TERM1   = 2.     *VMUT
         TERM2   = 2.*ALF4*VMUT*RKPE
         TERM3   = 2.*ALF5*VMUT*RKPE         

         R11     = R2O3*RK(L) - (TERM1*S11L +
     +        2.*TERM2*(S12(L)*W21L + S13(L)*W31L) -
     +        TERM3*(S11L**2 + S12(L)**2 + S13(L)**2 - SL**2/3.))
         R22     = R2O3*RK(L) - (TERM1*S22L +
     +        2.*TERM2*(S12(L)*W12(L) + S23(L)*W32L) -
     +        TERM3*(S12(L)**2 + S22L**2 + S23(L)**2 - SL**2/3.))
         R33     = R2O3*RK(L) - (TERM1*S33L +
     +        2.*TERM2*(S13(L)*W13(L) + S23(L)*W23(L)) -
     +        TERM3*(S13(L)**2 + S23(L)**2 + S33L**2 - SL**2/3.))
         R12     =            - (TERM1*S12(L) + TERM2*
     +        (S11L*W12(L)+S13(L)*W32L + S22L*W21L+S23(L)*W31L)-
     +        TERM3*(S11L*S12(L) + S12(L)*S22L + S13(L)*S23(L)))
         R13     =            - (TERM1*S13(L) + TERM2*
     +        (S11L*W13(L)+S12(L)*W23(L) + S23(L)*W21L+S33L*W31L)-
     +        TERM3*(S11L*S13(L) + S12(L)*S23(L) + S13(L)*S33L))
         R23     =            - (TERM1*S23(L) + TERM2*
     +        (S12(L)*W13(L)+S22L*W23(L) + S13(L)*W12(L)+S33L*W23(L))-
     +        TERM3*(S12(L)*S13(L) + S22L*S23(L) + S23(L)*S33L))

C ... ADD HOC TRANSITION
c          CALL IJKPAI(L,IMAX,JMAX,KMAX,LL,MM,NN)
c          IF(LL <= 22) THEN
c             R12 = 0.
c          ENDIF

         DUDXI   = S11(L)
         DVDYI   = S22(L)
         DWDZI   = S33(L)
         DUDYI   = S12(L) + W12(L)
         DUDZI   = S13(L) + W13(L)
         DVDXI   = S12(L) - W12(L)
         DVDZI   = S23(L) + W23(L)
         DWDXI   = S13(L) - W13(L)
         DWDYI   = S23(L) - W23(L)
         PTUR(L) = -(
     +        R11*DUDXI+R12*(DUDYI+DVDXI)+
     +        R22*DVDYI+R23*(DVDZI+DWDYI)+
     +        R33*DWDZI+R13*(DWDXI+DUDZI))

C ... MENTER'S MODIFICATION FOR THE MAXIMUM PRODUCTION
         
         PTUR(L) = MIN(PTUR(L),20.*EPSI)
 1234    FORMAT(1I6,5E14.5)

         UU(L) = R11
         UV(L) = R12
         UW(L) = R13
         VV(L) = R22
         VW(L) = R23
         WW(L) = R33
 4000 CONTINUE
      ENDIF
      IF(ITURB == 11) THEN
      DO 4500 JJ = 1,JMAX*ISTRID
         L       = IA + JJ

         DIV  = (S11(L)+S22(L)+S33(L))/3.

         S11L = S11(L)-DIV
         S22L = S22(L)-DIV
         S33L = S33(L)-DIV

         W31L = -W13(L) 
         W21L = -W12(L) 
         W32L = -W23(L) 

         SL  = S11L*S11L     + S12(L)*S12(L) + S13(L)*S13(L)
     +       + S12(L)*S12(L) + S22L*S22L     + S23(L)*S23(L)
     +       + S13(L)*S13(L) + S23(L)*S23(L) + S33L*S33L
         SL   = SQRT(SL)

         EPSI    = REPS(L) + DDEPS(L)
         RKPE    = RK(L)/(EPSI+EPS)
         RKPE2   = RK(L)**2/(EPSI+EPS)
         ETA     = (  .5*ALF3/ALF1*   SL*RKPE)**2
         XSI     = (SQR2*ALF2/ALF1*OHMI(L)*RKPE)**2
         DAMP2   = 3.*(1. + ETA)/(3. + ETA + 6.*XSI*ETA + 6.*XSI)*RKPE
         DAMP    = 2.*VIST(L)/ALF1/RK(L) !chien's model
         TERM1   = DAMP*ALF1*RK(L)
         TERM2   = DAMP*ALF2*RKPE2
         TERM3   = DAMP*ALF3*RKPE2

         R11     = R2O3*RK(L) - (TERM1*S11L +
     +        2.*TERM2*(S12(L)*W21L + S13(L)*W31L) -
     +        TERM3*(S11L**2 + S12(L)**2 + S13(L)**2 - SL**2/3.))
         R22     = R2O3*RK(L) - (TERM1*S22L +
     +        2.*TERM2*(S12(L)*W12(L) + S23(L)*W32L) -
     +        TERM3*(S12(L)**2 + S22L**2 + S23(L)**2 - SL**2/3.))
         R33     = R2O3*RK(L) - (TERM1*S33L +
     +        2.*TERM2*(S13(L)*W13(L) + S23(L)*W23(L)) -
     +        TERM3*(S13(L)**2 + S23(L)**2 + S33L**2 - SL**2/3.))
         R12     =            - (TERM1*S12(L) + TERM2*
     +        (S11L*W12(L)+S13(L)*W32L + S22L*W21L+S23(L)*W31L)-
     +        TERM3*(S11L*S12(L) + S12(L)*S22L + S13(L)*S23(L)))
         R13     =            - (TERM1*S13(L) + TERM2*
     +        (S11L*W13(L)+S12(L)*W23(L) + S23(L)*W21L+S33L*W31L)-
     +        TERM3*(S11L*S13(L) + S12(L)*S23(L) + S13(L)*S33L))
         R23     =            - (TERM1*S23(L) + TERM2*
     +        (S12(L)*W13(L)+S22L*W23(L) + S13(L)*W12(L)+S33L*W23(L))-
     +        TERM3*(S12(L)*S13(L) + S22L*S23(L) + S23(L)*S33L))

C ... ADD HOC TRANSITION
c          CALL IJKPAI(L,IMAX,JMAX,KMAX,LL,MM,NN)
c          IF(LL <= 22) THEN
c             R12 = 0.
c          ENDIF

         DUDXI   = S11(L)
         DVDYI   = S22(L)
         DWDZI   = S33(L)
         DUDYI   = S12(L) + W12(L)
         DUDZI   = S13(L) + W13(L)
         DVDXI   = S12(L) - W12(L)
         DVDZI   = S23(L) + W23(L)
         DWDXI   = S13(L) - W13(L)
         DWDYI   = S23(L) - W23(L)
         PTUR(L) = -(
     +        R11*DUDXI+R12*(DUDYI+DVDXI)+
     +        R22*DVDYI+R23*(DVDZI+DWDYI)+
     +        R33*DWDZI+R13*(DWDXI+DUDZI))

C ... MENTER'S MODIFICATION FOR THE MAXIMUM PRODUCTION
         
         PTUR(L) = MIN(PTUR(L),20.*EPSI)
 1235    FORMAT(1I6,5E14.5)

         UU(L) = R11
         UV(L) = R12
         UW(L) = R13
         VV(L) = R22
         VW(L) = R23
         WW(L) = R33
 4500 CONTINUE
      ENDIF
5000  CONTINUE

      RETURN
      END SUBROUTINE GATSKI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRODKE(OHMI,W12,W13,W23,S11,S12,S13,S22,S23,S33,
     +     PTUR,RK,REPS,DDEPS,TRUEPS,EPS2,VIST,VIS,RST11,RST22,
     +     RST33,RST12,RST13,RST23,VTRAN,IMAX,JMAX,KMAX,IN,JN,KN,
     +     M,ITURB,KATOL)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,M,ITURB,K,L,JJ,IA,
     +           ISTRID,JSTRID,KSTRID,IL,mmm,nnn,kkk

      REAL :: OHMI(*),W12(*),W13(*),W23(*),PTUR(*),RK(*),REPS(*),
     +        EPS2(*),VIS(*),S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),
     +        DDEPS(*),VIST(*),TRUEPS(*),RST11(*),RST22(*),RST33(*),
     +        RST12(*),RST13(*),RST23(*),VTRAN(*)

      REAL :: R2O3,DUDXI,DUDYI,DUDZI,DVDXI,DVDYI,DVDZI,
     +        DWDXI,DWDYI,DWDZI,
     +        DIV,VISCO,TRANS,pturl,STRAIN

      LOGICAL :: KATOL

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      R2O3    = 2./3.

      DO 5000 K = 1,KMAX
      IA  = (KN+K-1)*IL + JN*ISTRID

      DO 2000 JJ = 1,JMAX*ISTRID
         L       = IA + JJ

C ... Velocity derivatives are stored as strain rate and vorticity tensors

         DUDXI   = S11(L)
         DVDYI   = S22(L)
         DWDZI   = S33(L)
         DUDYI   = S12(L) + W12(L)
         DUDZI   = S13(L) + W13(L)
         DVDXI   = S12(L) - W12(L)
         DVDZI   = S23(L) + W23(L)
         DWDXI   = S13(L) - W13(L)
         DWDYI   = S23(L) - W23(L)

         STRAIN  = S11(L)*S11(L) + S12(L)*S12(L) + S13(L)*S13(L)
     +           + S12(L)*S12(L) + S22(L)*S22(L) + S23(L)*S23(L)
     +           + S13(L)*S13(L) + S23(L)*S23(L) + S33(L)*S33(L)
         STRAIN  = SQRT(2.*STRAIN)
         
         TRANS   = 0.

      IF(ITURB < 10) THEN ! An eddy viscosity model         

         DIV     = DUDXI + DVDYI + DWDZI
         VISCO   = (EPS2(L)-1.+TRANS)*VIS(L)
     
C ... Calculate the production of turbulence
     
         IF(KATOL) THEN ! Kato-Launder modeification

         PTUR(L) = VISCO*STRAIN*OHMI(L)

         ELSE

         PTUR(L) = VISCO*(2.*(DUDXI**2 + DVDYI**2 + DWDZI**2)
     1        + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2        + (DVDZI + DWDYI)**2 - R2O3*DIV**2) - R2O3*RK(L)*DIV

         ENDIF ! KATOL

         
C        PTUR(L) = VISCO*OHMI(L)**2       ! PRACTICALLY THE SAME
c        PTUR(L) = MIN(PTUR(L),PTURL)   ! PROPABLY USELESS (Antakee mersu)    

C ... Menter's modification for the maximum production
         
c        PTUR(L) = MIN(PTUR(L),20.*REPS(L))   ! ORIGINALLY 20.
         PTUR(L) = MIN(PTUR(L),20.*TRUEPS(L)) ! ORIGINALLY 20.
c            call ijkpai(l,imax,jmax,kmax,mmm,nnn,kkk)! You will need this

      ELSE IF(ITURB >= 21 .AND. ITURB <= 23) THEN ! A Reynolds stress model. 

C ... Obs! The Reynolds stresses RSTij are stored into FI or S11 arrays

         PTUR(L) = -(
     +    RST11(L)*DUDXI + RST12(L)*(DUDYI+DVDXI)+
     +    RST22(L)*DVDYI + RST23(L)*(DVDZI+DWDYI)+
     +    RST33(L)*DWDZI + RST13(L)*(DWDXI+DUDZI))

C ... This is the same as for eddy-viscosity models (for testing)

c        PTUR(L) = VISCO*(
c     +    (2.*S11(L) - R2O3*DIV)*DUDXI + 2.*S12(L)*(DUDYI+DVDXI)+
c     +    (2.*S22(L) - R2O3*DIV)*DVDYI + 2.*S23(L)*(DVDZI+DWDYI)+
c     +    (2.*S33(L) - R2O3*DIV)*DWDZI + 2.*S13(L)*(DWDXI+DUDZI))-
c     +     R2O3*RK(L)*DIV 
      ENDIF

C ... Transition via VTRAN

      PTUR(L) = PTUR(L)*VTRAN(L)

 2000 CONTINUE

 5000 CONTINUE

      RETURN
      END SUBROUTINE PRODKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TURFUN(S11,S12,S13,S22,S23,S33,U,V,W,RK,REPS,DDEPS,
     2     EPS2,VIST,VIS,PROK,PROV,A1,A2,A3,VOL,A1X,A1Y,A1Z,A2X,A2Y,
     3     A2Z,A3X,A3Y,A3Z,D1,D2,D3,XC,YC,ZC,IMAX,JMAX,KMAX,IN,JN,KN,
     4     ITURB,M,ZZZ,VTRAN,MAXW,IDERI,NGL)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,IDERI,M,IALA,IYLA,JALA,JYLA,
     2           KALA,KYLA,JJ,I,J,K,L,IA,IL,ISTRID,JSTRID,KSTRID,N,KA,
     3           ITURB,MAXW,MAW,NGL

      INTEGER, DIMENSION(6) :: IS

      REAL :: XC(*), YC(*), ZC(*)

      REAL :: U(*),VOL(*),V(*),A1(*),A2(*),A1X(*),A1Y(*),A2X(*),A2Y(*),
     2        W(*),A3(*),A1Z(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),RK(*),
     3        REPS(*),EPS2(*),VIS(*),PROK(*),PROV(*),D1(*),D2(*),D3(*),
     4        S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),DDEPS(*),
     5        VIST(*),VTRAN(*)

      REAL, POINTER :: DUDY(:),DUDZ(:),DVDX(:),DVDZ(:),DWDY(:),DWDX(:),
     2     DUDX(:),DVDY(:),DWDZ(:)

      REAL, TARGET :: ZZZ(MAXW)

      REAL :: W1,W2,W3,R11,R12,R22,R13,R23,R33,DX,DY,DZ,BETA,
     2        ALFA1,ALFA2,ALFA3

      REAL :: U1M,V1M,W1M,U1P,V1P,W1P,U2M,V2M,W2M,U2P,V2P,W2P,U3M,V3M,
     2        W3M,U3P,V3P,W3P,PVOL,DUDXI,DUDYI,DUDZI,DVDXI,DVDYI,DVDZI,
     3        DWDXI,DWDYI,DWDZI,XTRACE,DIV,R2O3

*      CALL DOMAW(MAXW,IMAX,JMAX,MAW)
*      CALL DOMAW(MAXW,IMAX,JMAX,IN,JN,MAW)

*      DUDY => ZZZ( 0*MAW+1: 1*MAW);DUDZ=> ZZZ( 1*MAW+1: 2*MAW)
*      DVDX => ZZZ( 2*MAW+1: 3*MAW);DVDZ=> ZZZ( 3*MAW+1: 4*MAW)
*      DWDX => ZZZ( 4*MAW+1: 5*MAW);DWDY=> ZZZ( 5*MAW+1: 6*MAW)
*      DUDX => ZZZ( 6*MAW+1: 7*MAW);DVDY=> ZZZ( 7*MAW+1: 8*MAW)
*      DWDZ => ZZZ( 8*MAW+1: 9*MAW)


C ... THIS SUBROUTINE CALCULATES REYNOLDS STRESSES FOR POST PROCESSING
C ... FOR 2 EQ AND ALGEBRAIC (zero eq.) TURBULENCE MODELS

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      
      KALA   = 1
      KYLA   = KMAX
      JALA   = 1
      JYLA   = JMAX
      IALA   = 1
      IYLA   = IMAX

      IF(M == 1) THEN  ! The ghostcells are included on the first level
         IF(KMAX > 1) THEN
            KALA = 0
            KYLA = KMAX + 1
         ENDIF
         JALA = 0
         JYLA = JMAX + 1
         IALA = 0
         IYLA = IMAX + 1
      ENDIF

      R2O3    = 2./3.


C ... Three different approaches are available for calculating gradients:

C ... IDERI = 1: Gauss approach weighted by distances
C ... IDERI = 2: Gauss approach
C ... IDERI = 3: Least-squares approach


      SELECT CASE(IDERI)

      CASE(1)  ! Gauss approach weighted by distances

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L-1) ) ) 
                  W2  = 2.0 - W1

                  U1M = A1(L)  *(W2*U(L) + W1*U(L-1))
                  V1M = A1(L)  *(W2*V(L) + W1*V(L-1))
                  W1M = A1(L)  *(W2*W(L) + W1*W(L-1))

                  W1  = 2.0*( D1(L) / ( D1(L) + D1(L+1) ) )
                  W2  = 2.0 - W1

                  U1P = A1(L+1)*(W2*U(L) + W1*U(L+1))
                  V1P = A1(L+1)*(W2*V(L) + W1*V(L+1))
                  W1P = A1(L+1)*(W2*W(L) + W1*W(L+1))

C ... ETA-DIRECTION

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L-ISTRID) ) )
                  W2  = 2.0 - W1

                  U2M = A2(L)       *(W2*U(L) + W1*U(L-ISTRID))
                  V2M = A2(L)       *(W2*V(L) + W1*V(L-ISTRID))
                  W2M = A2(L)       *(W2*W(L) + W1*W(L-ISTRID))

                  W1  = 2.0*( D2(L) / ( D2(L) + D2(L+ISTRID) ) )
                  W2  = 2.0 - W1

                  U2P = A2(L+ISTRID)*(W2*U(L) + W1*U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(W2*V(L) + W1*V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W2*W(L) + W1*W(L+ISTRID))

C ... ZETA DIRECTION

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L-IL) ) )
                  W2  = 2.0 - W1

                  U3M = A3(L)   *(W2*U(L) + W1*U(L-IL))
                  V3M = A3(L)   *(W2*V(L) + W1*V(L-IL))
                  W3M = A3(L)   *(W2*W(L) + W1*W(L-IL))

                  W1  = 2.0*( D3(L) / ( D3(L) + D3(L+IL) ) )
                  W2  = 2.0 - W1

                  U3P = A3(L+IL)*(W2*U(L) + W1*U(L+IL))
                  V3P = A3(L+IL)*(W2*V(L) + W1*V(L+IL))
                  W3P = A3(L+IL)*(W2*W(L) + W1*W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)

                  IF(ITURB >=3 .AND. ITURB < 8) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE) + R2O3*RK(L)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE) + R2O3*RK(L)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE) + R2O3*RK(L)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)
C ... make correction so that trace of Reynolds stresses is equal to RK.
                     S33(L)  = 2.*RK(L) - S11(L) - S22(L) ! Why?

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
      
                  ELSEIF(ITURB < 3 .OR. ITURB==8 .OR. ITURB==9) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROK(L) = VIST(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
                     PROK(L) = PROK(L) * VTRAN(L) ! Transition
                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(2) ! Simple Gaussian derivative calculation

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... XI-DIRECTION

                  U1M = A1(L)  *(U(L) + U(L-1))
                  V1M = A1(L)  *(V(L) + V(L-1))
                  W1M = A1(L)  *(W(L) + W(L-1))
                  U1P = A1(L+1)*(U(L) + U(L+1))
                  V1P = A1(L+1)*(V(L) + V(L+1))
                  W1P = A1(L+1)*(W(L) + W(L+1))

C ... ETA-DIRECTION

                  U2M = A2(L)       *(U(L) + U(L-ISTRID))
                  V2M = A2(L)       *(V(L) + V(L-ISTRID))
                  W2M = A2(L)       *(W(L) + W(L-ISTRID))
                  U2P = A2(L+ISTRID)*(U(L) + U(L+ISTRID))
                  V2P = A2(L+ISTRID)*(V(L) + V(L+ISTRID))
                  W2P = A2(L+ISTRID)*(W(L) + W(L+ISTRID))

C ... ZETA DIRECTION

                  U3M = A3(L)   *(U(L) + U(L-IL))
                  V3M = A3(L)   *(V(L) + V(L-IL))
                  W3M = A3(L)   *(W(L) + W(L-IL))
                  U3P = A3(L+IL)*(U(L) + U(L+IL))
                  V3P = A3(L+IL)*(V(L) + V(L+IL))
                  W3P = A3(L+IL)*(W(L) + W(L+IL))

                  PVOL  = .5/VOL(L) ! .5 is because of the averages above

                  DUDXI = PVOL*(
     1              A3X(L+IL)    *U3P - A3X(L)*U3M +
     2              A2X(L+ISTRID)*U2P - A2X(L)*U2M +
     3              A1X(L+1)     *U1P - A1X(L)*U1M)
                  DUDYI = PVOL*(
     1              A3Y(L+IL)    *U3P - A3Y(L)*U3M +
     2              A2Y(L+ISTRID)*U2P - A2Y(L)*U2M +
     3              A1Y(L+1)     *U1P - A1Y(L)*U1M)
                  DUDZI = PVOL*(
     1              A3Z(L+IL)    *U3P - A3Z(L)*U3M +
     2              A2Z(L+ISTRID)*U2P - A2Z(L)*U2M +
     3              A1Z(L+1)     *U1P - A1Z(L)*U1M)
                  DVDXI = PVOL*(
     1              A3X(L+IL)    *V3P - A3X(L)*V3M +
     2              A2X(L+ISTRID)*V2P - A2X(L)*V2M +
     3              A1X(L+1)     *V1P - A1X(L)*V1M)
                  DVDYI = PVOL*(
     1              A3Y(L+IL)    *V3P - A3Y(L)*V3M +
     2              A2Y(L+ISTRID)*V2P - A2Y(L)*V2M +
     3              A1Y(L+1)     *V1P - A1Y(L)*V1M)
                  DVDZI = PVOL*(
     1              A3Z(L+IL)    *V3P - A3Z(L)*V3M +
     2              A2Z(L+ISTRID)*V2P - A2Z(L)*V2M +
     3              A1Z(L+1)     *V1P - A1Z(L)*V1M)
                  DWDXI = PVOL*(
     1              A3X(L+IL)    *W3P - A3X(L)*W3M +
     2              A2X(L+ISTRID)*W2P - A2X(L)*W2M +
     3              A1X(L+1)     *W1P - A1X(L)*W1M)
                  DWDYI = PVOL*(
     1              A3Y(L+IL)    *W3P - A3Y(L)*W3M +
     2              A2Y(L+ISTRID)*W2P - A2Y(L)*W2M +
     3              A1Y(L+1)     *W1P - A1Y(L)*W1M)
                  DWDZI = PVOL*(
     1              A3Z(L+IL)    *W3P - A3Z(L)*W3M +
     2              A2Z(L+ISTRID)*W2P - A2Z(L)*W2M +
     3              A1Z(L+1)     *W1P - A1Z(L)*W1M)

                  IF(ITURB >=3 .AND. ITURB < 8) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE) + R2O3*RK(L)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE) + R2O3*RK(L)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE) + R2O3*RK(L)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)
C ... make correction so that trace of Reynolds stresses is equal to RK.
                     S33(L)  = 2.*RK(L) - S11(L) - S22(L) ! Why?

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
      
                  ELSEIF(ITURB < 3 .OR. ITURB==8 .OR. ITURB==9) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROK(L) = VIST(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
                     PROK(L) = PROK(L) * VTRAN(L) ! Transition
                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

C **********************************************************************

      CASE(3) ! Least-squares approach

C ... Computational Fluid Dynamics: Principles and Applications p. 163
C ... J. Blazek, Elsevier 2001

         IS(1) = -1
         IS(2) = -ISTRID
         IS(3) = -IL
         IS(4) =  1
         IS(5) =  ISTRID
         IS(6) =  IL

         DO K = KALA,KYLA
            KA = (KN+K-1)*IL
            DO J = JALA,JYLA
               JJ = (JN+J-1)*ISTRID + IN + KA 
               DO I = IALA,IYLA
                  L = JJ + I

C ... Compute the six entries in the upper triangular matrix.
C ... These entries should be pre-computed and stored for each 
C ... node due to efficiency reasons.

                  R11 = 0.0
                  R12 = 0.0
                  R22 = 0.0
                  R13 = 0.0
                  R23 = 0.0
                  R33 = 0.0

                  DO N = 1,6
                     DX  = XC(L + IS(N)) - XC(L)
                     DY  = YC(L + IS(N)) - YC(L)
                     DZ  = ZC(L + IS(N)) - ZC(L)
                     R11 = R11 + DX*DX
                     R12 = R12 + DX*DY
                     R22 = R22 + DY*DY
                     R13 = R13 + DX*DZ
                     R23 = R23 + DY*DZ
                     R33 = R33 + DZ*DZ
                  ENDDO

                  R11 = SQRT(R11)
                  R12 = R12 / R11
                  R22 = SQRT(R22-R12**2)
                  R13 = R13 / R11
                  R23 = (R23 - R12*R13) / R22
                  R33 = SQRT(R33-(R13**2 + R23**2))

                  BETA = (R12*R23 - R13*R22) / (R11*R22)

                  DUDXI = 0.0
                  DUDYI = 0.0
                  DUDZI = 0.0
                  DVDXI = 0.0
                  DVDYI = 0.0
                  DVDZI = 0.0
                  DWDXI = 0.0
                  DWDYI = 0.0
                  DWDZI = 0.0

                  DO N = 1,6
                     DX    = XC(L + IS(N)) - XC(L)
                     DY    = YC(L + IS(N)) - YC(L)
                     DZ    = ZC(L + IS(N)) - ZC(L)
                     ALFA1 =  DX / R11**2
                     ALFA2 = (DY - R12*DX/R11) / R22**2
                     ALFA3 = (DZ - R23*DY/R22 + BETA*DX) / R33**2
                     W1    = ALFA1 - R12*ALFA2/R11 + BETA*ALFA3
                     W2    = ALFA2 - R23*ALFA3/R22
                     W3    = ALFA3
                     DUDXI = DUDXI + W1*(U(L+IS(N))-U(L))
                     DUDYI = DUDYI + W2*(U(L+IS(N))-U(L))
                     DUDZI = DUDZI + W3*(U(L+IS(N))-U(L))
                     DVDXI = DVDXI + W1*(V(L+IS(N))-V(L))
                     DVDYI = DVDYI + W2*(V(L+IS(N))-V(L))
                     DVDZI = DVDZI + W3*(V(L+IS(N))-V(L))
                     DWDXI = DWDXI + W1*(W(L+IS(N))-W(L))
                     DWDYI = DWDYI + W2*(W(L+IS(N))-W(L))
                     DWDZI = DWDZI + W3*(W(L+IS(N))-W(L))
                  ENDDO

                  IF(ITURB >=3 .AND. ITURB < 8) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE) + R2O3*RK(L)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE) + R2O3*RK(L)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE) + R2O3*RK(L)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)
C ... make correction so that trace of Reynolds stresses is equal to RK.
                     S33(L)  = 2.*RK(L) - S11(L) - S22(L) ! Why?

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
      
                  ELSEIF(ITURB < 3 .OR. ITURB==8 .OR. ITURB==9) THEN

C ... Calculate Reynolds stresses from Boussinesq approximation

                     XTRACE  = (DUDXI + DVDYI + DWDZI)/3.
                     S11(L)  = -VIST(L)*2.*(DUDXI-XTRACE)
                     S22(L)  = -VIST(L)*2.*(DVDYI-XTRACE)
                     S33(L)  = -VIST(L)*2.*(DWDZI-XTRACE)
                     S12(L)  = -VIST(L)*(DUDYI + DVDXI)
                     S13(L)  = -VIST(L)*(DUDZI + DWDXI)
                     S23(L)  = -VIST(L)*(DWDYI + DVDZI)

                     DIV     = DUDXI + DVDYI + DWDZI

                     PROK(L) = VIST(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)
                     PROK(L) = PROK(L) * VTRAN(L) ! Transition
                     PROV(L) = VIS(L)*(2.*(DUDXI**2+DVDYI**2+DWDZI**2)
     1                       + (DUDYI + DVDXI)**2 + (DUDZI + DWDXI)**2
     2                       + (DVDZI + DWDYI)**2 - R2O3*DIV**2)

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

      END SELECT
C **********************************************************************
*     CALL EXTEND ?????

      RETURN
      END SUBROUTINE TURFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REFLEC(U,V,W,X,Y,Z,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,
     + JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,INRE,JNRE,KNRE,
     + GAMMA,E0REF,T0REF,FRSDEN,
     + FRSPRE,NBL,NPATCH,ICON,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,
     + EPS2,TEMP,C,DRDP,DRDH,RK,REPS,DDEPS,PRO,VAR,RGAS,T0,ITURB,
     + ISTATE,IB,MAXB,RCON,UWALL,VWALL,WWALL,TWALL,IHF,TLOLIM,TUPLIM,
     + SOLUTION_TYPE,MULPHL,REFLECL,TRANSL,TRM,
     + A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,PRC,IPRESC)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : LN, IC9

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     & IDI1,IDI2,IDI3,IN,JN,KN,INRE,JNRE,KNRE,
     & NBL,NPATCH,ITURB,ISTATE,IB,MAXB,ISTRID,JSTRID,
     & KSTRID,IL,IP,IBC,IPL,IFACE,I1,I2,J1,J2,KX1,KX2,KY1,KY2,IDIR,
*     & IST2,ISTR,JSTR,KSTR,ISTR3,JSTR3,KSTR3,I,J,K,ID,JST2,IUP,IAP,
     & ISTR,JSTR,KSTR,ISTR3,JSTR3,KSTR3,I,J,K,ID,
     & IHFL,IHEAT,KA3,IJ,N1,NR,L,NP,LB,KK,IC,LT,IPH,NTOT,mmm,nnn,lll,
     & LC,IPRESC

      INTEGER :: IHF(*), ICON(IC9,*)

      REAL :: U(*),V(*),W(*),
     &        RO(*),P(*),PDIFF(*),RM(*),RN(*),RW(*),E(*),VIS(*),
     &        VIST(*),CH(*),EPS2(*),TEMP(*),RK(*),REPS(*),DDEPS(*),C(*),
     &        DRDP(*),DRDH(*),RCON(IC9,*),UWALL(*),
     &        VWALL(*),WWALL(*),TWALL(*),
     &    A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      REAL :: GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,RGAS,XAUX,TSURF,
     &        VISB,EPSURF,SQRK,USURF2,TLOLIM,TUPLIM,T0,A11,A12,A13,A21,
     &        A22,A23,A31,A32,A33   ,UN,VN,WN     

      REAL :: Y(*), Z(*), X(*)

      REAL :: EINT(1),PSURF(1),ROSURF(1)
C      REAL, ALLOCATABLE :: zzz(:) ! You may need this

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(PRE_COR)          :: PRC(*)

      CHARACTER(LEN=10) :: SOLUTION_TYPE

      LOGICAL :: MULPHL,REFLECL,TRANSL

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      NTOT    = IL*KSTRID
C      ALLOCATE (zzz(NTOT))
      XAUX    = 0.
      IF (ITURB == 10) XAUX = 1.

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
    
      IF(IBC >= 7 .AND. IBC <= 10 .OR. IBC == 15 .OR. IBC == 3) THEN

      IPL   = ICON(2,IP)  ! proces local patch number
      IFACE = ICON(3,IP)
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)

C ... So-called standard wall-indeces for the patch data (mersu)

         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
      IF(IFACE == 2 .OR. IFACE == 5) THEN
         KY1  = ICON(4,IP)
         KY2  = ICON(5,IP)
         KX1  = ICON(6,IP)
         KX2  = ICON(7,IP)
      ENDIF
      IDIR    = 1
*      IST2    = 1
      IF(IFACE >= 4) IDIR = -1

C ... XI-DIRECTION

      IF(IFACE == 1 .OR. IFACE == 4) THEN

         IN = JNRE
         JN = KNRE
         KN = INRE 

*         ISTR = ISTRID  ! Not used
*         JSTR = IL
         KSTR = 1

         ISTR3 = ISTRID
         JSTR3 = IL
         KSTR3 = 1

         K    = IBOT
         IF(IBC == 3) K = 1
         ID   = IDI1
*         JST2 = JSTRID
*         IUP  = 1
*         IAP  = 1
         IF(IFACE == 4) THEN
            K   = ITOP
*            IUP  = 10
            IF(IBC == 7 .OR. IBC == 3) K = IMAX
         ENDIF

C ... ETA-DIRECTION

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

*         IN = INRE
*         JN = KNRE
*         KN = JNRE 
         IN = KNRE
         JN = INRE
         KN = JNRE 

*         ISTR = IL  ! Not used
*         JSTR = 1
         KSTR = ISTRID

         ISTR3 = IL
         JSTR3 = 1
         KSTR3 = ISTRID

         K    = JBOT
         IF(IBC ==3) K = 1
         ID   = IDI2
*         IST2 = KSTRID
*         JST2 = 1
*         IUP  = 4
*         IAP  = 2
         IF(IFACE == 5) THEN
            K   = JTOP
*            IUP  = 13
            IF(IBC == 7 .OR. IBC==3) K = JMAX
         ENDIF

C ... ZETA DIRECTION

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         IN = INRE
         JN = JNRE
         KN = KNRE 

*         ISTR = 1       ! Not used
*         JSTR = ISTRID
         KSTR = IL

         ISTR3 = 1
         JSTR3 = ISTRID
         KSTR3 = IL

         K = KBOT
         IF(IBC == 3) K = 1
         ID   = IDI3
*         JST2 = ISTRID
*         IUP  = 7
*         IAP  = 3
         IF(IFACE == 6) THEN
            K   = KTOP
*            IUP  = 16
            IF(IBC == 7 .OR. IBC == 3) K = KMAX
         ENDIF
      ENDIF
         
      IF(IBC == 7 .AND. IFACE <= 3) K = 1

      TSURF   = 0.
      IHFL    = 0
      IHEAT   = ICON(20,IP)

      IF(IHEAT == 1) TSURF = RCON(5,IP)
      IF(IHEAT == 1) IHFL  = 1
      IF(IHEAT == 3) TSURF = T0
      IF(IHEAT == 3) IHFL  = 0 ! Prevents a jump in the second cell

      KA3      = (KN+K-1)*KSTR3

C ... END OF SET-UP OF INDECES. FIRSTLY PRESSURES AND VELOCITIES

      IF(SOLUTION_TYPE /= 'SOLID') THEN

      IF(IBC == 3) THEN ! Extrapolate inlet pressures at the inlet

         DO J = KY1,KY2
            IJ        = (JN+J-1)*JSTR3 + KA3 
            DO I = KX1,KX2
               L      = 1 + (IN+I-1)*ISTR3 + IJ ! Cell index
               LB     = L - IDIR*KSTR
               LT     = L + IDIR*KSTR
               LC     = L - 2*IDIR*KSTR
               P(LB)  = P(L)
               PDIFF(LB) = PDIFF(L)
               P(LC)  = P(L)
               PDIFF(LC) = PDIFF(L)
            ENDDO 
         ENDDO

      ENDIF                     ! IBC ==  3 Inlet

C ... Firstly THERMODYNAMICS. (This was moved here 20.3.2017 from the end).

      IF(IBC >= 7 .AND. IBC <= 10) THEN
         DO 1200 J = KY1,KY2
         IJ        = (JN+J-1)*JSTR3 + KA3 
         DO 1200 I = KX1,KX2
            L      = 1 + (IN+I-1)*ISTR3 + IJ  ! Cell index
            LB     = L - IDIR*KSTR
            LT     = L + IDIR*KSTR
            LC     = L - 2*IDIR*KSTR
            TEMP(LB) = (1-IHFL)*(1.*TEMP(L) - 0.*TEMP(LT))+
c            TEMP(LB) = (1-IHFL)*(2.*TEMP(L) - TEMP(LT))+
     +        IHFL*(2.*TSURF-TEMP(L))
            TEMP(LB) = MIN(TEMP(LB),TUPLIM)
            TEMP(LB) = MAX(TEMP(LB),TLOLIM)
            TEMP(LC) = TEMP(LB)
            VIS(LB)  = VIS(L)
            CH(LB)   = CH(L)
            C(LB)    = C(L)
            DRDP(LB) = DRDP(L)
            DRDH(LB) = DRDH(L)
            IF(MULPHL) THEN
               PRO(LB)  = PRO(L)
               VAR(LB)  = VAR(L)
            ENDIF ! MULPHL

	      IF(TRANSL) THEN  !Intermittency variables
	         TRM(LB)  = TRM(L)	
	      ENDIF ! TRANSL
 1200    CONTINUE
      ENDIF                     ! IBC /= 15

      IF(IBC == 15 .AND. MULPHL) THEN ! Coupling with solid block
         DO 1300 J = KY1,KY2
         IJ        = (JN+J-1)*JSTR3 + KA3 
         DO 1300 I = KX1,KX2
            L      = 1 + (IN+I-1)*ISTR3 + IJ  ! Cell index
            LB     = L - IDIR*KSTR
            LT     = L + IDIR*KSTR
            VAR(LB)= VAR(L)
 1300 CONTINUE
      ENDIF                     ! IBC == 15

C ... Next velocities (antakee uuden mallinen mersu, 20.3.2017)

      IF(((ID /= 0 .AND. IFACE <= 3).OR.(ID > K .AND. IFACE >= 4))
     +        .AND. IBC /= 7 .AND. IBC /= 3) THEN ! FOR VISCOUS SURFACES

C ... A separate loop for the wall velocities (Aah NP, mersujen mersu)

         DO 500 J = KY1,KY2
         IJ       = (JN+J-1)*JSTR3 + KA3 
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO 500 I = KX1,KX2
            L       = 1 + (IN+I-1)*ISTR3 + IJ  ! Cell index
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR            
            LC      = L - 2*IDIR*KSTR          ! Second ghost cell
            NP = I  + NR                      ! Patch index with LN ghost cells

            IF(.NOT.REFLECL) THEN
               U(LB)   = - U(L) + 2.*UWALL(NP)
               V(LB)   = - V(L) + 2.*VWALL(NP)
               W(LB)   = - W(L) + 2.*WWALL(NP)
               U(LC)   =   U(LB)
               V(LC)   =   V(LB)
               W(LC)   =   W(LB)
               IF(SOLUTION_TYPE == 'MULTI') THEN
                  DO IPH = 1,NPHASES
                  VAR(LB)%U(IPH)   = - VAR(L)%U(IPH) + 2.*UWALL(NP)
                  VAR(LB)%V(IPH)   = - VAR(L)%V(IPH) + 2.*VWALL(NP)
                  VAR(LB)%W(IPH)   = - VAR(L)%W(IPH) + 2.*WWALL(NP)
                  VAR(LC)%U(IPH)   =   VAR(LB)%U(IPH)
                  VAR(LC)%V(IPH)   =   VAR(LB)%V(IPH)
                  VAR(LC)%W(IPH)   =   VAR(LB)%W(IPH)
                  ENDDO
               ENDIF

            ELSE IF(REFLECL) THEN

               U(LB)  = 2.*U(L) - 1.*U(LT) ! Try 1 and 0 in desporate cases
               V(LB)  = 2.*V(L) - 1.*V(LT)
               W(LB)  = 2.*W(L) - 1.*W(LT)
               U(LC)  =    U(LB)
               V(LC)  =    V(LB)
               W(LC)  =    W(LB)
               IF(SOLUTION_TYPE == 'MULTI') THEN
                  DO IPH = 1,NPHASES
                  VAR(LB)%U(IPH)  = 2.*VAR(L)%U(IPH) - 1.*VAR(LT)%U(IPH)
                  VAR(LB)%V(IPH)  = 2.*VAR(L)%V(IPH) - 1.*VAR(LT)%V(IPH)
                  VAR(LB)%W(IPH)  = 2.*VAR(L)%W(IPH) - 1.*VAR(LT)%W(IPH)
                  VAR(LC)%U(IPH)  =   VAR(LB)%U(IPH)
                  VAR(LC)%V(IPH)  =   VAR(LB)%V(IPH)
                  VAR(LC)%W(IPH)  =   VAR(LB)%W(IPH)
                  ENDDO
               ENDIF
            ENDIF

            IF(ISTATE /= 10) THEN ! Then realistic equation of state
               RO(LB)  =  MAX(0.001*FRSDEN,1.*RO(L)-0.*RO(LT))
               P(LB)   =  MAX(0.001*FRSPRE,2.* P(L)- P(LT))
               PDIFF(LB)= MAX(-0.999*FRSPRE,2.* PDIFF(L)- PDIFF(LT))
            ELSE ! Fully incompressible
               RO(LB)  =  1.*RO(L)-0.*RO(LT)
               P(LB)   =  2.* P(L)- P(LT)
               PDIFF(LB)= 2.* PDIFF(L)- PDIFF(LT)
            ENDIF

            EPS2(LB)= 1.       
            VIST(LB)= 0.  ! -VIST(L) shows error better? 2.2.19
            RM(LB)  = RO(LB)*U(LB)
            RN(LB)  = RO(LB)*V(LB)
            RW(LB)  = RO(LB)*W(LB)
            E(LB)   = E(L)

            RM(LC)  = RM(LB)
            RN(LC)  = RN(LB)
            RW(LC)  = RW(LB)
            E(LC)   = E(LB)
            EPS2(LC)= EPS2(LB)
            VIST(LC)= VIST(LB)
            PDIFF(LC)  = PDIFF(LB)
            P(LC)   = P(LB)
            RO(LC)  = RO(LB)

            IF(IPRESC == 1) THEN
            PRC(LB) = PRC(L)
            PRC(LC) = PRC(LB)
            ENDIF

            IF(MULPHL) THEN
               PRO(LB)  = PRO(L)
c               VAR(LB)  = VAR(L)
            ENDIF ! MULPHL

          IF(ITURB >= 3 .AND. ITURB /= 8) THEN
C ... in Speziales model distances should be transferred (antakee mersu):
c            DIST    = 1./D2(L)
c            SQRK    =IDIR*(9.*SQRT(RK(L)/RO(L))-
c     +                       SQRT(RK(LT)/RO(LT)))*DIST/3.
            SQRK    = 0.
            VISB    = VIS(L)
            EPSURF  = XAUX*2.*VISB*SQRK**2 ! = zero at the moment
            RK(LB)  = -RK(L)
            IF(ITURB /= 6) THEN
            REPS(LB)= -REPS(L) + 2.*EPSURF
            ELSE
            REPS(LB)= MAX(0.,2.*REPS(L) - REPS(LT))
            ENDIF
            DDEPS(LB) = MAX(0.,1.*DDEPS(L) - 0.*DDEPS(LT))
          ENDIF ! Turbulence quantities
 500    CONTINUE

      ELSEIF(IBC /= 3) THEN !INVISCID SURFACES OR SINGULARITIES (IBC=7)

         DO 1150 J = KY1,KY2
         IJ        = (JN+J-1)*JSTR3 + KA3 
         DO 1150 I = KX1,KX2
            L      = 1 + (IN+I-1)*ISTR3 + IJ  ! Cell index
            LB     = L - IDIR*KSTR
            LT     = L + IDIR*KSTR
            IF(ISTATE /= 10) THEN ! Then realistic equation of state
               RO(LB)  =  MAX(0.001*FRSDEN,1.*RO(L)-0.*RO(LT))
               P(LB)   =  MAX(0.001*FRSPRE,2.* P(L)- P(LT))
               PDIFF(LB)= MAX(-0.999*FRSPRE,2.* PDIFF(L)- PDIFF(LT))
            ELSE ! Fully incompressible
               RO(LB)  =  1.*RO(L)-0.*RO(LT)
               P(LB)   =  2.* P(L)- P(LT)
               PDIFF(LB)= 2.* PDIFF(L)- PDIFF(LT)
            ENDIF

C ... This was modified (January 2015). Testaamati

            IF(IFACE == 1) THEN
               A11 = A1X(L)
               A12 = A1Y(L)
               A13 = A1Z(L)
            ELSEIF(IFACE == 2) THEN
               A11 = A2X(L)
               A12 = A2Y(L)
               A13 = A2Z(L)
            ELSEIF(IFACE == 3) THEN
               A11 = A1X(L)
               A12 = A1Y(L)
               A13 = A1Z(L)
            ELSEIF(IFACE == 4) THEN
               A11 = A1X(L+KSTR)
               A12 = A1Y(L+KSTR)
               A13 = A1Z(L+KSTR)
            ELSEIF(IFACE == 5) THEN
               A11 = A2X(L+KSTR)
               A12 = A2Y(L+KSTR)
               A13 = A2Z(L+KSTR)
            ELSEIF(IFACE == 6) THEN
               A11 = A1X(L+KSTR)
               A12 = A1Y(L+KSTR)
               A13 = A1Z(L+KSTR)
            ENDIF

            CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

            UN = A11*U(L) + A12*V(L) + A13*W(L)
            VN = A21*U(L) + A22*V(L) + A23*W(L)
            WN = A31*U(L) + A32*V(L) + A33*W(L)
            UN = -UN

            U(LB) = A11*UN + A21*VN + A31*WN
            V(LB) = A12*UN + A22*VN + A32*WN
            W(LB) = A13*UN + A23*VN + A33*WN

            IF(IPRESC == 1) THEN

               UN = A11*PRC(L)%DPDX + A12*PRC(L)%DPDY + A13*PRC(L)%DPDZ
               VN = A21*PRC(L)%DPDX + A22*PRC(L)%DPDY + A23*PRC(L)%DPDZ
               WN = A31*PRC(L)%DPDX + A32*PRC(L)%DPDY + A33*PRC(L)%DPDZ
               UN = -UN

               PRC(LB)%DPDX = A11*UN + A21*VN + A31*WN
               PRC(LB)%DPDY = A12*UN + A22*VN + A32*WN
               PRC(LB)%DPDZ = A13*UN + A23*VN + A33*WN

            ENDIF

c            call ijkpai(lb,imax,jmax,kmax,mmm,nnn,lll)! You will need this
            IF(SOLUTION_TYPE == 'MULTI') THEN

            DO IPH = 1,NPHASES

              UN = A11*VAR(L)%U(IPH)+A12*VAR(L)%V(IPH)+A13*VAR(L)%W(IPH)
              VN = A21*VAR(L)%U(IPH)+A22*VAR(L)%V(IPH)+A23*VAR(L)%W(IPH)
              WN = A31*VAR(L)%U(IPH)+A32*VAR(L)%V(IPH)+A33*VAR(L)%W(IPH)
              UN = -UN

              VAR(LB)%U(IPH) = A11*UN + A21*VN + A31*WN
              VAR(LB)%V(IPH) = A12*UN + A22*VN + A32*WN
              VAR(LB)%W(IPH) = A13*UN + A23*VN + A33*WN

            ENDDO

            ENDIF ! MULTI

c            U(LB)  = 2.*U(L) - 1.*U(LT) ! oli 2 ja -1 1.10.2008-21.1.2015
c            V(LB)  = 2.*V(L) - 1.*V(LT)
c            W(LB)  = 2.*W(L) - 1.*W(LT)

            EPS2(LB) = 1.*EPS2(L) - 0.*EPS2(LT)
            VIST(LB) = 1.*VIST(L) - 0.*VIST(LT) ! tänne

            IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Singularity

               RK(LB)    = MAX(0.,1.*RK(L)    - 0.*RK(LT))  
               REPS(LB)  = MAX(0.,1.*REPS(L)  - 0.*REPS(LT))
               DDEPS(LB) = MAX(0.,1.*DDEPS(L) - 0.*DDEPS(LT))
               EPS2(LB)  = 1.*EPS2(L) - 0.*EPS2(LT)
               VIST(LB)  = 1.*VIST(L) - 0.*VIST(LT)

            ENDIF
            
 1150    CONTINUE

      ENDIF ! OF INVISCID SURFACES OR SINGULARITIES (IBC=7)


      ELSE IF(SOLUTION_TYPE == 'SOLID') THEN ! A solid block


      IF(((ID /= 0 .AND. IFACE <= 3).OR.(ID > K .AND. IFACE >= 4))
     +        .AND. IBC /= 7) THEN ! FOR VISCOUS SURFACES

         DO J = KY1,KY2

            IJ = (JN+J-1)*JSTR3 + KA3 
            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 

            DO I = KX1,KX2

               L  = 1 + (IN+I-1)*ISTR3 + IJ ! Cell index
               LB = L - IDIR*KSTR
               LT = L + IDIR*KSTR
               NP = I + NR      ! Patch index with LN ghost cells

               U(LB)  = - U(L) + 2.*UWALL(NP)
               V(LB)  = - V(L) + 2.*VWALL(NP)
               W(LB)  = - W(L) + 2.*WWALL(NP)
               RM(LB) = RO(LB)*U(LB)
               RN(LB) = RO(LB)*V(LB)
               RW(LB) = RO(LB)*W(LB)
               E(LB)  = E(L)

               IF(ITURB >= 3 .AND. ITURB /= 8) THEN
                  RK(LB)    = 0.
                  REPS(LB)  = 0.
                  DDEPS(LB) = 0.
               ENDIF            ! OF VISCOUS SURFACES

            ENDDO
         ENDDO

      ELSE     ! INVISCID SURFACES OR SINGULARITIES (IBC=7)

         DO J = KY1,KY2
            IJ = (JN+J-1)*JSTR3 + KA3 
            DO I = KX1,KX2
               L        = 1 + (IN+I-1)*ISTR3 + IJ ! Cell index
               LB       = L - IDIR*KSTR
               LT       = L + IDIR*KSTR

               RO(LB)   = MAX(0.001*FRSDEN,2.*RO(L)-RO(LT))
               P(LB)    = 0.
               V(LB)    = 0.
               W(LB)    = 0.
               EPS2(LB) = 0.
               VIST(LB) = 0.
               IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Singularity
                  RK(LB)    = 0.  
                  REPS(LB)  = 0.
                  DDEPS(LB) = 0.
               ENDIF
            ENDDO
         ENDDO
      ENDIF ! OF INVISCID SURFACES OR SINGULARITIES (IBC=7)

C ... THEN THERMODYNAMICS

      IF(IBC /= 15) THEN
         DO J = KY1,KY2
            IJ = (JN+J-1)*JSTR3 + KA3 
            DO I = KX1,KX2

               L  = 1 + (IN+I-1)*ISTR3 + IJ ! Cell index
               LB = L - IDIR*KSTR
               LT = L + IDIR*KSTR
               LC = L - 2*IDIR*KSTR

               TEMP(LB) = (1-IHFL)*(2.*TEMP(L) - TEMP(LT))+
     +              IHFL*(2.*TSURF-TEMP(L))
               TEMP(LB) = MIN(TEMP(LB),TUPLIM)
               TEMP(LB) = MAX(TEMP(LB),TLOLIM)
               TEMP(LT) = TEMP(LB)
               IF(MULPHL) THEN
                 DO IPH = 1,NPHASES
                    PRO(LB)%TEMP(IPH) = (1-IHFL)*(2.*PRO(L)%TEMP(IPH) - 
     +                   PRO(LT)%TEMP(IPH))+
     +                   IHFL*(2.*TSURF-PRO(L)%TEMP(IPH))
                    PRO(LB)%TEMP(IPH) = MIN(PRO(LB)%TEMP(IPH),TUPLIM)
                    PRO(LB)%TEMP(IPH) = MAX(PRO(LB)%TEMP(IPH),TLOLIM)
                    PRO(LB)%DTEMP(IPH) = (1-IHFL)*(2.*PRO(L)%DTEMP(IPH)- 
     +                   PRO(LT)%DTEMP(IPH))+
     +                   IHFL*(2.*TSURF-PRO(L)%DTEMP(IPH))
                    PRO(LB)%DTEMP(IPH) = MIN(PRO(LB)%DTEMP(IPH),TUPLIM)
                    PRO(LB)%DTEMP(IPH) = MAX(PRO(LB)%DTEMP(IPH),TLOLIM)
                 ENDDO
               ENDIF
               IF(TRANSL) THEN  !Intermittency variables
                  TRM(LB)  = TRM(L)
               ENDIF
               CH(LB)   = CH(L)
            ENDDO
         ENDDO
      ENDIF                     ! IBC /= 15

      ENDIF                     ! SOLUTION_TYPE /= 'SOLID'
      ENDIF                     ! IBC >= 7 .AND. IBC <= 15

7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE REFLEC


      SUBROUTINE REFCEN(U,V,W,VAR,TRM,IMAX,JMAX,KMAX,IBOT,ITOP,
     + JBOT,JTOP,KBOT,KTOP,INRE,JNRE,KNRE,FRSDEN,FRSPRE,NPATCH,ICON,
     + RO,P,PDIFF,E,RM,RN,RW,RK,REPS,ITURB,ISTATE,SOLUTION_TYPE,TRANSL)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9,RKLIM,EPSLIM
      
      IMPLICIT NONE
      
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, INRE, JNRE, KNRE, IPH,
     1           ISTRID,JSTRID,KSTRID,IL,IP,IBC,IFACE,I1,I2,J1,J2,IDIR,
     2           ISTR,JSTR,KSTR,I,J,K,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     3           NPATCH,KK,L,LB,LT,ISTATE,ITURB
      INTEGER :: ICON(IC9,*)

      REAL    :: FRSDEN,FRSPRE
      REAL    :: U(*),V(*),W(*),
     1           RO(*),P(*),RM(*),RN(*),RW(*),E(*),
     2             RK(*),REPS(*),PDIFF(*)

      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      CHARACTER(LEN=10) :: SOLUTION_TYPE

      LOGICAL :: TRANSL


C ... THIS SUBROUTINE IS FOR CENTRAL DIFFERENCE CALCULATION NEAR WALLS

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
      IF(IBC >= 7 .AND. IBC <= 10) THEN
      IFACE = ICON(3,IP)
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)
      IDIR    = 1
      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN

         IN = JNRE
         JN = KNRE
         KN = INRE

         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = IBOT
         IF(IFACE == 4) THEN
            K   = ITOP
            IF(IBC == 7) K = IMAX
         ENDIF
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

         IN = INRE
         JN = KNRE
         KN = JNRE

         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         K    = JBOT
         IF(IFACE == 5) THEN
            K   = JTOP
            IF(IBC == 7) K = JMAX
         ENDIF
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         IN = INRE
         JN = JNRE
         KN = KNRE

         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K = KBOT
         IF(IFACE == 6) THEN
            K   = KTOP
            IF(IBC == 7) K = KMAX
         ENDIF
      ENDIF
      IF(IBC == 7 .AND. IFACE <= 3) K = 1

C ... END OF SET-UP OF INDECES. FIRSTLY PRESSURES AND VELOCITIES

         DO 1100 J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO 1100 I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            IF(ISTATE /= 10) THEN ! Then realistic equation of state
               RO(LB)  =  MAX(0.001*FRSDEN,2.*RO(L)-RO(LT))
               P(LB)   =  MAX(0.001*FRSPRE,2.* P(L)- P(LT))
               PDIFF(LB)= MAX(-0.999*FRSPRE,2.* PDIFF(L)- PDIFF(LT))
            ELSE ! Fully incompressible
               RO(LB)  =  2.*RO(L)-RO(LT)
               P(LB)   =  2.* P(L)- P(LT)
               PDIFF(LB)= 2.* PDIFF(L)- PDIFF(LT)
            ENDIF
            E(LB)  = 2.*E(L)  -  E(LT)
            U(LB)  = 2.*U(L)  -  U(LT)
            V(LB)  = 2.*V(L)  -  V(LT)
            W(LB)  = 2.*W(L)  -  W(LT)
            RM(LB) = 2.*RM(L) - RM(LT)
            RN(LB) = 2.*RN(L) - RN(LT)
            RW(LB) = 2.*RW(L) - RW(LT)
            IF(SOLUTION_TYPE == 'MULTI') THEN
               DO IPH = 1,NPHASES
               VAR(LB)%U(IPH)  = 2.*VAR(L)%U(IPH) - VAR(LT)%U(IPH)
               VAR(LB)%V(IPH)  = 2.*VAR(L)%V(IPH) - VAR(LT)%V(IPH)
               VAR(LB)%W(IPH)  = 2.*VAR(L)%W(IPH) - VAR(LT)%W(IPH)
               ENDDO
            ENDIF

 1100    CONTINUE

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         DO 1150 J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO 1150 I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            RK(LB)   =  MAX(RKLIM,2.*RK(L)    - 1.*RK(LT))
            REPS(LB) =  MAX(EPSLIM,2.*REPS(L)  - 1.*REPS(LT))
C            DDEPS(LB)=  1.*DDEPS(L) - 0.*DDEPS(LT) ei ollut vanhaskaan
	      IF(TRANSL) THEN  !Intermittency variables
	         TRM(LB)  = TRM(L)	
	      ENDIF ! TRANSL
 1150    CONTINUE
          ENDIF

      ENDIF                     ! IBC >= 7 .AND. IBC <= 10
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE REFCEN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REFCOR(U,V,W,X,Y,Z,IMAX,JMAX,KMAX,IN,JN,KN,
     + FRSDEN,FRSPRE,NBL,NPATCH,
     + RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,DRDH,
     + RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRANSL,TRM,MULPHL,ITURB,NSCAL,
     + MAXSB,MAXEB)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ITURB,NSCAL,MAXSB,MAXEB,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,
     +           NBL,I,J,K,ISTR,JSTR,KSTR,IL,ISTRID,JSTRID,KSTRID

      REAL :: Y(*), Z(*), X(*)

      REAL :: U(*),V(*),W(*),
     3        RO(*),P(*),PDIFF(*),RM(*),RN(*),RW(*),E(*),VIS(*),
     4        VIST(*),CH(*),EPS2(*),TEMP(*),RK(*),REPS(*),DDEPS(*),C(*),
     5        DRDP(*),DRDH(*),FI(MAXSB,MAX(1,NSCAL)),BIJ(MAXEB,*)

      REAL :: FRSPRE,FRSDEN

      LOGICAL :: MULPHL,TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      CHARACTER(LEN=10) :: SOLUTION_TYPE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      ISTR   = 1
      JSTR   = ISTRID
      KSTR   = IL

C ************************************************************************
C     FILL EMPTY CORNERS
C ********************************************************************
C ... ZETA-DIRECTION

      CALL REFVAR(U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,
     +     C,DRDP,DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,MULPHL,
     +     TRANSL,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,NSCAL,MAXSB,
     +     MAXEB)

C ********************************************************************
C ... XI-DIRECTION

      CALL REFVAR(U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,
     +     C,DRDP,DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,MULPHL,
     +     TRANSL,JMAX,KMAX,IMAX,JSTR,KSTR,ISTR,JN,KN,IN,NSCAL,MAXSB,
     +     MAXEB)

C ********************************************************************
C ... ETA DIRECTION

      CALL REFVAR(U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,
     +     C,DRDP,DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,MULPHL,
     +     TRANSL,KMAX,IMAX,JMAX,KSTR,ISTR,JSTR,KN,IN,JN,NSCAL,MAXSB,
     +     MAXEB)

C ********************************************************************

C ... CORNER POINT (1,1,1)

      I       = 0
      J       = 0
      K       = 0
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (IMAX,1,1)

      I       = IMAX+1
      J       = 0
      K       = 0
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (1,JMAX,1)

      I       = 0
      J       = JMAX+1
      K       = 0
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (IMAX,JMAX,1)

      I       = IMAX+1
      J       = JMAX+1
      K       = 0
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (1,1,KMAX)

      I       = 0
      J       = 0
      K       = KMAX+1
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (IMAX,1,KMAX)

      I       = IMAX+1
      J       = 0
      K       = KMAX+1
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (1,JMAX,KMAX)

      I       = 0
      J       = JMAX+1
      K       = KMAX+1
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

C ... CORNER POINT (IMAX,JMAX,KMAX)

      I       = IMAX+1
      J       = JMAX+1
      K       = KMAX+1
       CALL EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     + DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,
     + MULPHL,TRANSL)

      RETURN
      END SUBROUTINE REFCOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REFVAR(U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,
     +     TEMP,C,DRDP,DRDH,RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,
     +     MULPHL,TRANSL,
     +     IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,NSCAL,MAXSB,MAXEB)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: NSCAL,MAXSB,MAEB,IMAX,JMAX,KMAX,IN,JN,KN,
     &           ISTR,JSTR,KSTR,ITURB,MAXEB,I1,J1,K1,I2,J2,K2,
     &           K,KK,IG,JG,NS,LLL2,LLU2,LUL2,LUU2,LLL,LLU,LUL,LUU

      REAL :: U(*),V(*),W(*),RO(*),P(*),PDIFF(*),RM(*),RN(*),RW(*),
     +        E(*),VIS(*),VIST(*),CH(*),EPS2(*),TEMP(*),RK(*),REPS(*),
     +        DDEPS(*),C(*),DRDP(*),DRDH(*),
     +        FI(MAXSB,MAX(1,NSCAL)),BIJ(MAXEB,*)

      LOGICAL :: MULPHL, TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      I1  = 1
      J1  = 1
      K1  = 1
      I2  = IMAX
      J2  = JMAX
      K2  = KMAX
C ... (IG,JG) = (1,1) --------- REAL GHOSTCELL
C ... (IG,JG) = (2,1) --------- CORNER GHOSTCELL
C ... (IG,JG) = (1,2) --------- CORNER GHOSTCELL
C ... (IG,JG) = (2,2) --------- BEHIND REAL GHOSTCELL
      DO 1000 K = K1,K2
         KK     = (K-1+KN)*KSTR
         LLL2   = 1+(IN+I1-1   )*ISTR+(JN+J1-1   )*JSTR+KK
         LLU2   = 1+(IN+I1-1   )*ISTR+(JN+J2-1   )*JSTR+KK
         LUL2   = 1+(IN+I2-1   )*ISTR+(JN+J1-1   )*JSTR+KK
         LUU2   = 1+(IN+I2-1   )*ISTR+(JN+J2-1   )*JSTR+KK
      DO 1000 IG = 1,2
      DO 1000 JG = 1,2
         LLL    = 1+(IN+I1-1-IG)*ISTR+(JN+J1-1-JG)*JSTR+KK
         LLU    = 1+(IN+I1-1-IG)*ISTR+(JN+J2-1+JG)*JSTR+KK
         LUL    = 1+(IN+I2-1+IG)*ISTR+(JN+J1-1-JG)*JSTR+KK
         LUU    = 1+(IN+I2-1+IG)*ISTR+(JN+J2-1+JG)*JSTR+KK
C ... LOWER-LOWER
         RO(LLL)   = RO(LLL2)
         RM(LLL)   = RM(LLL2)
         RN(LLL)   = RN(LLL2)
         RW(LLL)   = RW(LLL2)
         E(LLL)    = E(LLL2)
         U(LLL)    = U(LLL2)
         V(LLL)    = V(LLL2)
         W(LLL)    = W(LLL2)
         P(LLL)    = P(LLL2)
         PDIFF(LLL)= PDIFF(LLL2)
         EPS2(LLL) = EPS2(LLL2)
         VIST(LLL) = VIST(LLL2)
         TEMP(LLL) = TEMP(LLL2)
         VIS(LLL)  = VIS(LLL2)
         CH(LLL)   = CH(LLL2)
         C(LLL)    = C(LLL2)
         DRDP(LLL) = DRDP(LLL2)
         DRDH(LLL) = DRDH(LLL2)

C ... LOWER-UPPER
         RO(LLU)   = RO(LLU2)
         RM(LLU)   = RM(LLU2)
         RN(LLU)   = RN(LLU2)
         RW(LLU)   = RW(LLU2)
         E(LLU)    = E(LLU2)
         U(LLU)    = U(LLU2)
         V(LLU)    = V(LLU2)
         W(LLU)    = W(LLU2)
         P(LLU)    = P(LLU2)
         PDIFF(LLU)= PDIFF(LLU2)
         EPS2(LLU) = EPS2(LLU2)
         VIST(LLU) = VIST(LLU2)
         TEMP(LLU) = TEMP(LLU2)
         VIS(LLU)  = VIS(LLU2)
         CH(LLU)   = CH(LLU2)
         C(LLU)    = C(LLU2)
         DRDP(LLU) = DRDP(LLU2)
         DRDH(LLU) = DRDH(LLU2)

C ... UPPER-LOWER
         RO(LUL)   = RO(LUL2)
         RM(LUL)   = RM(LUL2)
         RN(LUL)   = RN(LUL2)
         RW(LUL)   = RW(LUL2)
         E(LUL)    = E(LUL2)
         U(LUL)    = U(LUL2)
         V(LUL)    = V(LUL2)
         W(LUL)    = W(LUL2)
         P(LUL)    = P(LUL2)
         PDIFF(LUL)= PDIFF(LUL2)
         EPS2(LUL) = EPS2(LUL2)
         VIST(LUL) = VIST(LUL2)
         TEMP(LUL) = TEMP(LUL2)
         VIS(LUL)  = VIS(LUL2)
         CH(LUL)   = CH(LUL2)
         C(LUL)    = C(LUL2)
         DRDP(LUL) = DRDP(LUL2)
         DRDH(LUL) = DRDH(LUL2)

C ... UPPER-UPPER
         RO(LUU)   = RO(LUU2)
         RM(LUU)   = RM(LUU2)
         RN(LUU)   = RN(LUU2)
         RW(LUU)   = RW(LUU2)
         E(LUU)    = E(LUU2)
         U(LUU)    = U(LUU2)
         V(LUU)    = V(LUU2)
         W(LUU)    = W(LUU2)
         P(LUU)    = P(LUU2)
         PDIFF(LUU)= PDIFF(LUU2)
         EPS2(LUU) = EPS2(LUU2)
         VIST(LUU) = VIST(LUU2)
         TEMP(LUU) = TEMP(LUU2)
         VIS(LUU)  = VIS(LUU2)
         CH(LUU)   = CH(LUU2)
         C(LUU)    = C(LUU2)
         DRDP(LUU) = DRDP(LUU2)
         DRDH(LUU) = DRDH(LUU2)

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            RK(LLL)   = RK(LLL2)  
            REPS(LLL) = REPS(LLL2)
            DDEPS(LLL)= DDEPS(LLL2)

            RK(LLU)   = RK(LLU2)  
            REPS(LLU) = REPS(LLU2)
            DDEPS(LLU)= DDEPS(LLU2)

            RK(LUL)   = RK(LUL2)  
            REPS(LUL) = REPS(LUL2)
            DDEPS(LUL)= DDEPS(LUL2)

            RK(LUU)   = RK(LUU2)  
            REPS(LUU) = REPS(LUU2)
            DDEPS(LUU)= DDEPS(LUU2)
         ENDIF
         IF(NSCAL >= 1) THEN
            DO 900 NS = 1,NSCAL
               FI(LLL,NS) = FI(LLL2,NS)
               FI(LLU,NS) = FI(LLU2,NS)
               FI(LUL,NS) = FI(LUL2,NS)
               FI(LUU,NS) = FI(LUU2,NS)
 900        CONTINUE
         ENDIF
         IF(MAXEB > 1) THEN
            DO 910 NS = 1,6
               BIJ(LLL,NS) = BIJ(LLL2,NS)
               BIJ(LLU,NS) = BIJ(LLU2,NS)
               BIJ(LUL,NS) = BIJ(LUL2,NS)
               BIJ(LUU,NS) = BIJ(LUU2,NS)
 910        CONTINUE
         ENDIF
         IF(MULPHL) THEN
c            DO 920 NS = 1,NPHASES
c               PRO(LLL)%RO(NS)   = PRO(LLL2)%RO(NS)
c               PRO(LLU)%RO(NS)   = PRO(LLU2)%RO(NS)
c               PRO(LUL)%RO(NS)   = PRO(LUL2)%RO(NS)
c               PRO(LUU)%RO(NS)   = PRO(LUU2)%RO(NS)
c               PRO(LLL)%E(NS)    = PRO(LLL2)%E(NS)
c               PRO(LLU)%E(NS)    = PRO(LLU2)%E(NS)
c               PRO(LUL)%E(NS)    = PRO(LUL2)%E(NS)
c               PRO(LUU)%E(NS)    = PRO(LUU2)%E(NS)
c               PRO(LLL)%TEMP(NS) = PRO(LLL2)%TEMP(NS)
c               PRO(LLU)%TEMP(NS) = PRO(LLU2)%TEMP(NS)
c               PRO(LUL)%TEMP(NS) = PRO(LUL2)%TEMP(NS)
c               PRO(LUU)%TEMP(NS) = PRO(LUU2)%TEMP(NS)
c               VAR(LLL)%X(NS)    = VAR(LLL2)%X(NS)
c               VAR(LLU)%X(NS)    = VAR(LLU2)%X(NS)
c               VAR(LUL)%X(NS)    = VAR(LUL2)%X(NS)
c               VAR(LUU)%X(NS)    = VAR(LUU2)%X(NS)
c               VAR(LLL)%ALFA(NS) = VAR(LLL2)%ALFA(NS)
c               VAR(LLU)%ALFA(NS) = VAR(LLU2)%ALFA(NS)
c               VAR(LUL)%ALFA(NS) = VAR(LUL2)%ALFA(NS)
c               VAR(LUU)%ALFA(NS) = VAR(LUU2)%ALFA(NS)
c               PRO(LLL)%VIS(NS)  = PRO(LLL2)%VIS(NS)
c               PRO(LLU)%VIS(NS)  = PRO(LLU2)%VIS(NS)
c               PRO(LUL)%VIS(NS)  = PRO(LUL2)%VIS(NS)
cc               PRO(LUU)%VIS(NS)  = PRO(LUU2)%VIS(NS)
c               PRO(LLL)%C(NS)    = PRO(LLL2)%C(NS)
c               PRO(LLU)%C(NS)    = PRO(LLU2)%C(NS)
c               PRO(LUL)%C(NS)    = PRO(LUL2)%C(NS)
c               PRO(LUU)%C(NS)    = PRO(LUU2)%C(NS)
c               PRO(LLL)%DRODP(NS)= PRO(LLL2)%DRODP(NS)
c               PRO(LLU)%DRODP(NS)= PRO(LLU2)%DRODP(NS)
c               PRO(LUL)%DRODP(NS)= PRO(LUL2)%DRODP(NS)
c               PRO(LUU)%DRODP(NS)= PRO(LUU2)%DRODP(NS)
c               PRO(LLL)%DRODH(NS)= PRO(LLL2)%DRODH(NS)
c               PRO(LLU)%DRODH(NS)= PRO(LLU2)%DRODH(NS)
c               PRO(LUL)%DRODH(NS)= PRO(LUL2)%DRODH(NS)
c               PRO(LUU)%DRODH(NS)= PRO(LUU2)%DRODH(NS)
C ... Näinkö reflectoi osuuskauppaväki?
               PRO(LLL) = PRO(LLL2)
               PRO(LLU) = PRO(LLU2)
               PRO(LUL) = PRO(LUL2)
               PRO(LUU) = PRO(LUU2)
               VAR(LLL) = VAR(LLL2)
               VAR(LLU) = VAR(LLU2)
               VAR(LUL) = VAR(LUL2)
               VAR(LUU) = VAR(LUU2)
c 920        CONTINUE
         ENDIF ! MULPHL

         IF(TRANSL) THEN ! Intermittency variables
               TRM(LLL) = TRM(LLL2)
               TRM(LLU) = TRM(LLU2)
               TRM(LUL) = TRM(LUL2)
               TRM(LUU) = TRM(LUU2)
         ENDIF       
    
 1000 CONTINUE

      RETURN
      END SUBROUTINE REFVAR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C	
      SUBROUTINE EXTCOR(I,J,K,IMAX,JMAX,KMAX,IN,JN,KN,
     + U,V,W,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,DRDH,
     + RK,REPS,DDEPS,FI,BIJ,PRO,VAR,TRM,ITURB,NSCAL,MAXSB,MAXEB,MULPHL,
     + TRANSL)
      
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ITURB,NSCAL,MAXSB,MAXEB,I,J,K,IMAX,JMAX,KMAX,
     &           IN,JN,KN,IL,IDIR,JDIR,KDIR,ISTRID,JSTRID,
     &           IBASE,JBASE,KBASE,LBASE,I1,J1,K1,L,NS

      REAL :: U(*),V(*),W(*),RO(*),P(*),PDIFF(*),RM(*),RN(*),RW(*),
     &        E(*),VIS(*),VIST(*),CH(*),EPS2(*),TEMP(*),RK(*),REPS(*),
     &        DDEPS(*),C(*),DRDP(*),DRDH(*),FI(MAXSB,MAX(1,NSCAL)),
     &        BIJ(MAXEB,*)

      LOGICAL :: MULPHL, TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
      IDIR    = -1
      JDIR    = -1
      KDIR    = -1
      IF (IMAX+1 == I) IDIR = 1
      IF (JMAX+1 == J) JDIR = 1
      IF (KMAX+1 == K) KDIR = 1
      IBASE   = I - IDIR
      JBASE   = J - JDIR
      KBASE   = K - KDIR
         
      LBASE   = 1+(IN+IBASE-1) + (JN+JBASE-1)*ISTRID + 
     +     (KN+KBASE-1)*IL
      DO 1000 K1 = 1,2
      DO 1000 J1 = 1,2
      DO 1000 I1 = 1,2
      L   = 1+(IN+IBASE-1+IDIR*I1) + (JN+JBASE-1+JDIR*J1)*ISTRID + 
     +     (KN+KBASE-1+KDIR*K1)*IL

          RO(L)   = RO(LBASE)
          RM(L)   = RM(LBASE)
          RN(L)   = RN(LBASE)
          RW(L)   = RW(LBASE)
          E(L)    = E(LBASE)
          U(L)    = U(LBASE)
          V(L)    = V(LBASE)
          W(L)    = W(LBASE)
          P(L)    = P(LBASE)
          PDIFF(L)= PDIFF(LBASE)
          EPS2(L) = EPS2(LBASE)
          VIST(L) = VIST(LBASE)
          TEMP(L) = TEMP(LBASE)
          VIS(L)  = VIS(LBASE)
          CH(L)   = CH(LBASE)
          C(L)    = C(LBASE)
          DRDP(L) = DRDP(LBASE)
          DRDH(L) = DRDH(LBASE)
          IF(ITURB >= 3 .AND. ITURB /= 8) THEN
             RK(L)   = RK(LBASE)  
             REPS(L) = REPS(LBASE)
             DDEPS(L)= DDEPS(LBASE)
          ENDIF
          IF(NSCAL >= 1) THEN
             DO 900 NS  = 1,NSCAL
               FI(L,NS) = FI(LBASE,NS)
 900        CONTINUE
         ENDIF
          IF(MAXEB > 1) THEN
             DO 910 NS  = 1,6
               BIJ(L,NS) = BIJ(LBASE,NS)
 910        CONTINUE
         ENDIF
          IF(MULPHL) THEN
c             DO 910 NS  = 1,6
               PRO(L) = PRO(LBASE)
               VAR(L) = VAR(LBASE)
c 910        CONTINUE
          ENDIF
          IF(TRANSL) THEN  !Intermittency variables
               TRM(L) = TRM(LBASE)
          ENDIF
 1000 CONTINUE

      RETURN
      END SUBROUTINE EXTCOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE RICHAR(OHMI,W12,W13,W23,S11,S12,S13,S22,S23,S33,LRICH,
     +     U,V,W,RK,REPS,DDEPS,A1,A2,A3,VOL,A1X,A1Y,A1Z,
     +     A2X,A2Y,A2Z,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IN,JN,KN,RI)

C ... THIS SUBROUTINE CALCULATES THE RICHARSON-LIKE NUMBER 
C ... Andrei Khodak and Charles Hirsch. Second Order Non-Linear k-e
C ... Models with Explicit Effect of curvature and Rotation. Proceedings
C ... of the Third ECCOMAS Computational Fluid Dynamics Conference, 9-13.
C ... September 1996, Paris, France. 
C ... See also: Muistio CFD/TERMO-17b-97

      REAL :: U(*),VOL(*),V(*),A1(*),A2(*),A1X(*),A1Y(*),A2X(*),
     +        A2Y(*),W(*),A3(*),A1Z(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +        OMX,OMY,OMZ,RK(*),REPS(*),DDEPS(*),S,RI(*),
     +        OHMI(*),W12(*),W13(*),W23(*),S11(*),S12(*),S13(*),
     +        S22(*),S23(*),S33(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      R2O3    = 2./3.
      SQR2    = SQRT(2.)
      EPS     = 1.E-10

      DO 5000 K = 1,KMAX
         IA  = (KN+K-1)*IL + JN*ISTRID

         IF(LRICH == 1) THEN  ! ANTTI
         DO JJ = 1,JMAX*ISTRID
            L  = IA + JJ

            DIV     = 2.*(S11(L) + S22(L) + S33(L))/3.
            SIJ11   = (S11(L) - DIV)
            SIJ22   = (S22(L) - DIV)
            SIJ33   = (S33(L) - DIV)
            SIJ12   = S12(L)
            SIJ13   = S13(L)
            SIJ23   = S23(L)

            S  = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +           + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +           + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
            S  = SQRT(2.*S)
C ... CALCULATE RICHARDSON NUMBER
            RAPU    = MAX(OHMI(L),EPS)/(MAX(S,EPS))
            RI(L) = RAPU*(RAPU - 1.0) 
            RI(L) = MIN(RI(L), 1.0E6)
         ENDDO
         ELSEIF(LRICH == 2) THEN  ! HIRSCH
         DO JJ = 1,JMAX*ISTRID
            L  = IA + JJ

            DIV     = 2.*(S11(L) + S22(L) + S33(L))/3.
            SIJ11   = (S11(L) - DIV)
            SIJ22   = (S22(L) - DIV)
            SIJ33   = (S33(L) - DIV)
            SIJ12   = S12(L)
            SIJ13   = S13(L)
            SIJ23   = S23(L)

            S  = SIJ11*SIJ11 + SIJ12*SIJ12 + SIJ13*SIJ13
     +           + SIJ12*SIJ12 + SIJ22*SIJ22 + SIJ23*SIJ23
     +           + SIJ13*SIJ13 + SIJ23*SIJ23 + SIJ33*SIJ33
            S  = SQRT(2.*S)
C ... CALCULATE RICHARDSON NUMBER
            T       = RK(L) / (REPS(L)+DDEPS(L)+EPS) !TURBULENCE TIME-SCALE
            RI(L) = -1. * T**2 * OHMI(L) * (S - OHMI(L))
         ENDDO
         ENDIF
5000  CONTINUE

      RETURN
      END SUBROUTINE RICHAR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE PERCHE(IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     +     ALPHA,BETA,NPATCH,ICON,XCO,YCO,ZCO,N,PERCHL,MULPHC)

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: ICON(IC9,*)

      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, INRE, JNRE, KNRE

      REAL :: XCO(*), YCO(*), ZCO(*)

      LOGICAL :: PERCHL

      CHARACTER(*) :: MULPHC

C
C ... check periodic patches if they are inlet or outlets

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      U10     = COS(ALPHA)*COS(BETA)
      V10     = SIN(ALPHA)
      W10     = COS(ALPHA)*SIN(BETA)

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
         IBC     = ICON(1,IP)
      IF(IBC == 1.AND.ICON(20,IP) == 2 .OR.  !?
     +   IBC == 1.AND.ICON(20,IP) == 3 .OR.  !?
     +   IBC == 1.AND.ICON(20,IP) == 1 .OR. 
     +   IBC == 17) THEN
         PERCHL = .TRUE.
      IFACE   = ICON(3,IP)
      I1    = ICON(4,IP)
      I2    = ICON(5,IP) + 1
      J1    = ICON(6,IP)
      J2    = ICON(7,IP) + 1
      IDIR  = 1
      K     = 1
      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         IF(IFACE == 4) K   = IMAX+1
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = INRE
         JN   = KNRE
         KN   = JNRE
         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         IF(IFACE == 5) K   = JMAX+1
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         IF(IFACE == 6) K   = KMAX+1
      ENDIF
         
      AX = 0.
      AY = 0.
      AZ = 0.
      NN = 0
      DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO I = I1,I2
            NN = NN + 1
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR

            AX      = AX + XCO(LT)-XCO(L)
            AY      = AY + YCO(LT)-YCO(L)
            AZ      = AZ + ZCO(LT)-ZCO(L)
         ENDDO
      ENDDO
      SUM = SQRT(AX**2 + AY**2 + AZ**2)
      AX = AX/SUM
      AY = AY/SUM
      AZ = AZ/SUM

      DIR = U10*AX + V10*AY + W10*AZ

      IF(DIR >=   .5) ICON(20,IP) = 2 ! inlet
      IF(DIR <=  -.5) ICON(20,IP) = 3 ! outlet

      NBG = ICON(24,IP)

      IF(DIR >=   .5) THEN
         WRITE(45,313) NBG,IP,' inlet'
      ELSEIF(DIR <=  -.5) THEN
         WRITE(45,313) NBG,IP,'outlet'
      ELSE
         WRITE(*,314)  NBG,IP
         WRITE(13,314) NBG,IP
         WRITE(45,314) NBG,IP
         WRITE(45,*) DIR
      ENDIF
        
 313  FORMAT('  Block',I3,' pathc',I3,' is periodic ',A6)
 314  FORMAT('  Block',I3,' pathc',I3,' is periodic, but unknown type.',
     + /,'  Side periodicity is assumed, check the situation.')

      ENDIF                     ! IBC == 1 .AND. ICON(20,IP) /= 0
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE PERCHE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE PERNE1(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,
     +     PRO,VAR,TRM,XC,YC,ZC,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     +     FRSDEN,FRSVEL,ALPHA,BETA,NBL,NPATCH,ICON,
     +     A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     ITURB,MAXB,BAL,INLRC,OUTRC,MULPHC)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +           NBL,NPATCH,ITURB,
     +           MAXB,ISTRID,JSTRID,KSTRID,NTOT,IL,L,IP,IBC,IFACE,
     +           I1,I2,J1,J2,IDIR,ISTR,JSTR,KSTR,IUP,IAP,I2M,J2M,I,J,K,
     +           KK,LB,LT,LA,IOUT,INLRC,OUTRC,INLRCS

      REAL :: TEMP(*),RO(*),P(*),PDIFF(*),U(*),V(*),W(*),E(*),
     +        RK(*),REPS(*),BAL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),
     +        A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      REAL :: FRSDEN,FRSVEL,ALPHA,BETA,U10,V10,W10,XMASS,XMAST,PAVE,
     +        PAVG,AREA,EAVE,EAVC,EIN,E1,E2,UFACE,VFACE,WFACE

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      REAL :: XC(*), YC(*), ZC(*)
      
      REAL, ALLOCATABLE :: AA(:,:),AAXA(:,:)

      CHARACTER(*) :: MULPHC

C ... Periodic inlet and outlets. First the global values are calculated 
C ... over all patches. PPR 2.10.98. Arrangement of A1 and A1X changed 2.2.04
C ... Two-phase missing

C ... BAL is balance table
C ... what                       INLET I  OUTLET I  Variable
C ... mass flux                     1        8       XMASS
C ... mass flux                     2        9       XMAST
C ... P*A in computational domain   3       10       PAVE             
C ... P*A in ghost cells            4       11       PAVG    
C ... area                          5       12       AREA 
C ... T*A                           6       13       EAVE 
C ... T*A in computational domain   7       14       EAVC             

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
      IL      = ISTRID*JSTRID

      U10     = FRSVEL*COS(ALPHA)*COS(BETA)
      V10     = FRSVEL*SIN(ALPHA)
      W10     = FRSVEL*COS(ALPHA)*SIN(BETA)

C ... Antakee mersu

      ALLOCATE(AA(NTOT,3),AAXA(NTOT,9))

      DO L = 1,NTOT
         AA(L,1)   = A1(L)
         AA(L,2)   = A2(L)
         AA(L,3)   = A3(L)
         AAXA(L,1) = A1X(L)
         AAXA(L,2) = A1Y(L)
         AAXA(L,3) = A1Z(L)
         AAXA(L,4) = A2X(L)
         AAXA(L,5) = A2Y(L)
         AAXA(L,6) = A2Z(L)
         AAXA(L,7) = A3X(L)
         AAXA(L,8) = A3Y(L)
         AAXA(L,9) = A3Z(L)
      ENDDO

C ... Obsolate block above, AAX could be replaced by a patch-sized array

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
      IF(IBC == 1.AND.ICON(20,IP) >= 2 .OR. IBC == 17) THEN ! inlet or
      IFACE   = ICON(3,IP)                                  ! outlet
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)
      IDIR  = 1
      K     = 1
      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         IUP  = 1
         IAP  = 1
         I2M  = MIN0(I2,JMAX)
         J2M  = MIN0(J2,KMAX)
         IF(IFACE == 4) THEN
            K = IMAX
            IF(IBC == 7) K = IMAX
         ENDIF
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = INRE
         JN   = KNRE
         KN   = JNRE
         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         IUP  = 4
         IAP  = 2
         I2M  = MIN0(I2,IMAX)
         J2M  = MIN0(J2,KMAX)
         IF(IFACE == 5) THEN
            K = JMAX
            IF(IBC == 7) K = JMAX
         ENDIF
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         IUP  = 7
         IAP  = 3
         I2M  = MIN0(I2,IMAX)
         J2M  = MIN0(J2,JMAX)
         IF(IFACE == 6) THEN
            K = KMAX
            IF(IBC == 7) K = KMAX
         ENDIF
      ENDIF
         
      XMASS = 0. ! mass flux
      XMAST = 0. ! target mass
      PAVE  = 0. ! average pressure in computational domain
      PAVG  = 0. ! average pressure in ghost cells
      AREA  = 0. ! area
      EAVE  = 0. ! average internal energy
      EAVC  = 0. ! average internal in computational domain

      DO J = MAX0(1,J1),J2M
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO I = MAX0(1,I1),I2M
            L       = 1 + (IN+I-1)*ISTR + KK ! 1. computational cell
            LB      = L - IDIR*KSTR          ! 1. ghost cell
            LT      = L + IDIR*KSTR          ! 2. computational cell
            LA      = L + (1-IDIR)*KSTR/2    ! cell face

C ... Simple central difference, the same should be applied in INFPER

      IF(ICON(20,IP) == 2) INLRCS = INLRC
      IF(ICON(20,IP) == 3) INLRCS = OUTRC

      SELECT CASE(INLRCS)

c      CASE(0,1,2) ! Central
c            UFACE   = U(L)
c            VFACE   = V(L)
c            WFACE   = W(L)
      CASE(0,1,2) ! Central ! Miksei tämä toimi???? Kait se toimii...
            UFACE   = .5*(U(L)+U(LB))
            VFACE   = .5*(V(L)+V(LB))
            WFACE   = .5*(W(L)+W(LB))
      CASE(3) ! 1st order
            UFACE   = U(LB)
            VFACE   = V(LB)
            WFACE   = W(LB)
      END SELECT
c      write (6000,*)UFACE,U(L),U(LB)

            XMASS = XMASS + RO(LA)*AA(LA,IAP)*(AAXA(LA,IUP)*UFACE+ 
     +           AAXA(LA,IUP+1)*VFACE+AAXA(LA,IUP+2)*WFACE)

            XMAST = XMAST + AA(LA,IAP)*FRSDEN*(AAXA(LA,IUP)*U10+ 
     +           AAXA(LA,IUP+1)*V10+AAXA(LA,IUP+2)*W10)
            AREA = AREA + AA(LA,IAP)
            PAVE = PAVE + AA(LA,IAP)*(2.*PDIFF(L)- PDIFF(LT))
            PAVG = PAVG + AA(LA,IAP)*PDIFF(LB)

            EIN  = TEMP(LB)
            EAVE = EAVE + AA(LA,IAP)*EIN

            E1   = TEMP(L)
            E2   = TEMP(LT)

            EAVC = EAVC + AA(LA,IAP)*(2.*E1-E2)

         ENDDO
      ENDDO

      PAVE = PAVE/AREA
      PAVG = PAVG/AREA
      EAVE = EAVE/AREA
      EAVC = EAVC/AREA

      IOUT  = 0
      IF(ICON(20,IP) == 3) IOUT = 7
      BAL(1+IOUT) = BAL(1+IOUT) + XMASS
      BAL(2+IOUT) = BAL(2+IOUT) + XMAST   
      BAL(3+IOUT) = BAL(3+IOUT) + PAVE*AREA   
      BAL(4+IOUT) = BAL(4+IOUT) + PAVG*AREA      
      BAL(5+IOUT) = BAL(5+IOUT) + AREA     
      BAL(6+IOUT) = BAL(6+IOUT) + EAVE*AREA      
      BAL(7+IOUT) = BAL(7+IOUT) + EAVC*AREA      

      ENDIF                     ! IBC >= 1 .AND. ICON....
7000  CONTINUE ! LOOP OVER THE PATCHES

      DEALLOCATE (AA,AAXA)
      RETURN
      END SUBROUTINE PERNE1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE PERNE2(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,PRO,VAR,
     +     TRM,XC,YC,ZC,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     +     GAMMA,E0,T0,FRSDEN,FRSPRE,FRSVEL,FRSSIE,FRSTEM,
     +     A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     NBL,NPATCH,ICON,ITURB,MAXB,BAL,INLRC,OUTRC,MULPHC)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9, GX, GY, GZ, GROUND

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +           NBL,NPATCH,ITURB,
     +           MAXB,ISTRID,JSTRID,KSTRID,NTOT,IL,L,IP,IBC,IFACE,
     +           I1,I2,J1,J2,IDIR,ISTR,JSTR,KSTR,IUP,IAP,IOUT,I,J,K,KK,
     +           LB,LB2,LT,LA,INLRC,OUTRC,INLRCS

      REAL :: RO(*),TEMP(*),P(*),PDIFF(*),U(*),V(*),W(*),E(*),RK(*),
     +        REPS(*),BAL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      REAL :: GAMMA,E0,T0,FRSDEN,FRSVEL,FRSSIE,FRSTEM,   
     +        XMASS,XMAST,AREA,PAVE,PAVG,EAVE,EAVC,EIN,EI2,REL,DELTAP,
     +        FRSPRE,UTARGET

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      REAL :: XC(*), YC(*), ZC(*)

      REAL, ALLOCATABLE :: AA(:,:),AAXA(:,:)

      CHARACTER(*) :: MULPHC

C ... Periodic inlet and outlets. Pressure, massflux and internal
C ... energy is manipulated. PPR 21.9.98. Two-phase missing!

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
      IL      = ISTRID*JSTRID

C ... Antakee mersu

      ALLOCATE(AA(NTOT,3),AAXA(NTOT,9))

      DO L = 1,NTOT
         AA(L,1)   = A1(L)
         AA(L,2)   = A2(L)
         AA(L,3)   = A3(L)
         AAXA(L,1) = A1X(L)
         AAXA(L,2) = A1Y(L)
         AAXA(L,3) = A1Z(L)
         AAXA(L,4) = A2X(L)
         AAXA(L,5) = A2Y(L)
         AAXA(L,6) = A2Z(L)
         AAXA(L,7) = A3X(L)
         AAXA(L,8) = A3Y(L)
         AAXA(L,9) = A3Z(L)
      ENDDO

C ... Obsolate block above, AAX could be replaced by a patch-sized array

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
      IF(IBC == 1.AND.ICON(20,IP) >= 2 .OR. IBC == 17) THEN ! inlet or
      IFACE   = ICON(3,IP)                                  !  outlet
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)
      IDIR    = 1
      K       = 1
      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         IUP  = 1
         IAP  = 1
         IF(IFACE == 4) THEN
            K   = IMAX
         ENDIF
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = INRE
         JN   = KNRE
         KN   = JNRE
         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         IUP  = 4
         IAP  = 2
         IF(IFACE == 5) THEN
            K   = JMAX
         ENDIF
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         IUP  = 7
         IAP  = 3
         IF(IFACE == 6) THEN
            K   = KMAX
         ENDIF
      ENDIF
         
      IOUT  = 0
      IF(ICON(20,IP) == 3) IOUT = 7

      XMASS = BAL(1+IOUT)       ! mass flux
      XMAST = BAL(2+IOUT)       ! target mass
      AREA  = BAL(5+IOUT)       ! area
      PAVE  = BAL(3+IOUT)/AREA  ! average pressure in computational domain
      PAVG  = BAL(4+IOUT)/AREA  ! average pressure in ghost cells
      EAVE  = BAL(6+IOUT)/AREA  ! average temperature
      EAVC  = BAL(7+IOUT)/AREA  ! average temperature in computational domain

C ... Central-difference correction (check this if discretization is changed)

c      REL   = (XMAST+2.*(XMAST-XMASS))/XMASS ! XMAST/XMASS ! First-order upw

c      REL   = XMAST/XMASS ! First-order upw

      IF(ICON(20,IP) == 2) INLRCS = INLRC
      IF(ICON(20,IP) == 3) INLRCS = OUTRC

      SELECT CASE(INLRCS)

      CASE(0,1,2) ! Central
            REL   = XMAST/XMASS 
c            REL   = (XMASS+2.*(XMAST-XMASS))/XMASS ! Toimiiko??? Ei
      CASE(3) ! 1st order
            REL   = XMAST/XMASS 
      END SELECT


      IF(ICON(20,IP) == 2 .OR. IBC == 17) THEN ! inlet
      DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO I = I1,I2
            L    = 1 + (IN+I-1)*ISTR + KK ! 1. computational cell
            LB   = L - IDIR*KSTR          ! 1. ghost cell
            LB2  = L - 2*IDIR*KSTR        ! 2. ghost cell
            LT   = L + IDIR*KSTR          ! 2. computational cell
            LA   = L + (1-IDIR)*KSTR/2    ! cell face

C ... Pressure must be extrapolated, temperature is adjusted to the inlet
C ... value and velocities to the target mass flow.

            EIN      = TEMP(LB)
            TEMP(LB) = EIN - EAVE + FRSTEM
            EI2      = TEMP(LB2)
            TEMP(LB2)= EI2 - EAVE + FRSTEM
            
c            UTARGET= U(LB)*REL
            U(LB)  = U(LB )*REL
c            write(667,*) nbl,xmass,xmast,rel,u(lb)
            V(LB)  = V(LB )*REL
            W(LB)  = W(LB )*REL
            U(LB2) = U(LB2)*REL
            V(LB2) = V(LB2)*REL
            W(LB2) = W(LB2)*REL
            PDIFF(LB)  = PDIFF(LB)  + (PAVE-PAVG)
            PDIFF(LB2) = PDIFF(LB2) + (PAVE-PAVG)
            DELTAP = FRSDEN*(GX*(XC(LB) - GROUND) +
     &                       GY*(YC(LB) - GROUND) +
     &                       GZ*(ZC(LB) - GROUND))
            P(LB)      = PDIFF(LB) + FRSPRE + DELTAP
            DELTAP = FRSDEN*(GX*(XC(LB2) - GROUND) +
     &                       GY*(YC(LB2) - GROUND) +
     &                       GZ*(ZC(LB2) - GROUND))
            P(LB2)      = PDIFF(LB2) + FRSPRE + DELTAP

C ... Turbulence variables are not changed !

         ENDDO
      ENDDO

      ELSEIF(ICON(20,IP) == 3) THEN ! outlet

      DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK ! 1. computational cell
            LB      = L - IDIR*KSTR          ! 1. ghost cell
            LB2     = L - 2*IDIR*KSTR        ! 2. ghost cell
            LT      = L + IDIR*KSTR          ! 2. computational cell
            LA      = L + (1-IDIR)*KSTR/2    ! cell face

            PDIFF(LB)  = PDIFF(LB)  -PAVG
            PDIFF(LB2) = PDIFF(LB2) -PAVG
            DELTAP = FRSDEN*(GX*(XC(LB) - GROUND) +
     &                       GY*(YC(LB) - GROUND) +
     &                       GZ*(ZC(LB) - GROUND))
            P(LB)      = PDIFF(LB) + FRSPRE + DELTAP
            DELTAP = FRSDEN*(GX*(XC(LB2) - GROUND) +
     &                       GY*(YC(LB2) - GROUND) +
     &                       GZ*(ZC(LB2) - GROUND))
            P(LB2)      = PDIFF(LB2) + FRSPRE + DELTAP

            TEMP(LB ) = 2.*TEMP(L ) - TEMP(LT)
            TEMP(LB2) = 2.*TEMP(LB) - TEMP(L )

            EIN  = TEMP(LB)
c            TEMP(LB) = EIN - EAVE + EAVC
            EI2  = TEMP(LB2)
c            TEMP(LB2) = EI2 - EAVE + EAVC

         ENDDO
      ENDDO
      ENDIF                     ! ICON(20,...

      ENDIF                     ! IBC >= 1 .AND. ICON....
7000  CONTINUE ! LOOP OVER THE PATCHES

      DEALLOCATE (AA,AAXA)
 
      RETURN
      END SUBROUTINE PERNE2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE PERNE3(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,PRO,VAR,
     +     TRM,XC,YC,ZC,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     +     GAMMA,E0,T0,FRSDEN,FRSPRE,FRSVEL,FRSSIE,FRSTEM,
     +     A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     NBL,NPATCH,ICON,ITURB,MAXB,BAL,VIST,EPS2,MULPHL,TRANSL,
     +     MULPHC)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9, GX, GY, GZ, GROUND

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +           NBL,NPATCH,ITURB,
     +           MAXB,ISTRID,JSTRID,KSTRID,NTOT,IL,L,IP,IBC,IFACE,
     +           I1,I2,J1,J2,IDIR,ISTR,JSTR,KSTR,IUP,IAP,I,J,K,KK,
     +           LB,LB2,LT,LA,KU,KUP,LU,LU2

      REAL :: RO(*),TEMP(*),P(*),PDIFF(*),U(*),V(*),W(*),E(*),RK(*),
     +        REPS(*),BAL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),VIST(*),EPS2(*)

      REAL :: GAMMA,E0,T0,FRSDEN,FRSVEL,FRSSIE,FRSTEM,   
     +        XMASS,XMAST,AREA,PAVE,PAVG,EAVE,EAVC,EIN,EI2,REL,DELTAP,
     +        FRSPRE

      REAL :: XC(*), YC(*), ZC(*)

      LOGICAL :: MULPHL, TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)   
      TYPE(INTERMITTENCY)    :: TRM(*) 

      CHARACTER(*) :: MULPHC

c      REAL, ALLOCATABLE :: AA(:,:),AAXA(:,:)


C ... Periodic inlet and outlets. Pressure, massflux and internal
C ... energy is manipulated. PPR 21.9.98
C ... This routine is only for the circulating inlet (INP)

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
      IL      = ISTRID*JSTRID

C ... LOOP OVER THE PATCHES OF THE BLOCK 
          
      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)

      IF(IBC == 17) THEN ! inlet only periodic
      IFACE   = ICON(3,IP)
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)
      IDIR    = 1
      K       = 1

      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         IUP  = 1
         IAP  = 1
         KUP  = IMAX
         IF(IFACE == 4) THEN
            K   = IMAX
            KUP = 1
            IF(IBC == 7) K = IMAX
         ENDIF
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = INRE
         JN   = KNRE
         KN   = JNRE
         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         IUP  = 4
         IAP  = 2
         KUP  = JMAX
         IF(IFACE == 5) THEN
            K   = JMAX
            KUP = 1
            IF(IBC == 7) K = JMAX
         ENDIF
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         IUP  = 7
         IAP  = 3
         KUP  = KMAX
         IF(IFACE == 6) THEN
            K   = KMAX
            KUP = 1
            IF(IBC == 7) K = KMAX
         ENDIF
      ENDIF
       

      IF(ICON(20,IP) == 2) THEN ! inlet extended as in connections
      DO J = J1-1,J2+1
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         KU      = (JN+J-1)*JSTR + (KN+KUP-1)*KSTR
         DO I = I1-1,I2+1
            L      = 1 + (IN+I-1)*ISTR + KK ! 1. computational cell
            LB     = L - IDIR*KSTR          ! 1. ghost cell
            LB2    = L - 2*IDIR*KSTR        ! 2. ghost cell
            LU     = 1 + (IN+I-1)*ISTR + KU ! Last computational cell
            LU2    = LU- IDIR*KSTR          ! 2. last comp. cell

            LT     = L + IDIR*KSTR          ! 2. computational cell
            LA     = L + (1-IDIR)*KSTR/2    ! cell face

            U(LB)  = U(LU )
            V(LB)  = V(LU )
            W(LB)  = W(LU )
            U(LB2) = U(LU2)
            V(LB2) = V(LU2)
            W(LB2) = W(LU2)

            TEMP(LB)   = TEMP(LU)
            TEMP(LB2)  = TEMP(LU2)
            RO(LB)     = RO(LU)
            RO(LB2)    = RO(LU2)

            RK(LB)     = RK(LU)
            RK(LB2)    = RK(LU2)
            REPS(LB)   = REPS(LU)
            REPS(LB2)  = REPS(LU2)
            VIST(LB)   = VIST(LU)
            VIST(LB2)  = VIST(LU2)
            EPS2(LB)   = EPS2(LU)
            EPS2(LB2)  = EPS2(LU2)
          
            PDIFF(LB)  = PDIFF(LU)
            PDIFF(LB2) = PDIFF(LU2)
            DELTAP = FRSDEN*(GX*(XC(LB) - GROUND) +
     &                       GY*(YC(LB) - GROUND) +
     &                       GZ*(ZC(LB) - GROUND))
            P(LB)      = PDIFF(LB) + FRSPRE + DELTAP
            DELTAP = FRSDEN*(GX*(XC(LB2) - GROUND) +
     &                       GY*(YC(LB2) - GROUND) +
     &                       GZ*(ZC(LB2) - GROUND))
            P(LB2)      = PDIFF(LB2) + FRSPRE + DELTAP

            IF(MULPHL) THEN ! Multiphase variables ! Two-fluid missing
              PRO(LB)  = PRO(LU)
              PRO(LB2) = PRO(LU2)
            ENDIF ! MULPHL

            IF(TRANSL) THEN ! Intermittency variables
              TRM(LB)  = TRM(LU)
              TRM(LB2) = TRM(LU2)
            ENDIF ! TRANSL

         ENDDO
      ENDDO
      ENDIF                     ! ICON(20,...

      ENDIF                     ! IBC >= 17 .AND. ICON....
7000  CONTINUE ! LOOP OVER THE PATCHES
 
      RETURN
      END SUBROUTINE PERNE3
C
C ----------------------------------------------------------------------
C --- Omega equation source terms by A. Hellsten -----------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE FUN1KO(FUN1,CROSSD,VIS,RO,RK,ROMEGA,DISTW,TRANSL,
     &         IMAX,JMAX,KMAX,IN,JN,KN,BSTAR,SIGOM2,RKLIM,CHLREF,RE,NGL)

C     Calculate the blending function FUN1 of Menter's SST model.
C                                               AH 31.1.1997

      REAL :: FUN1(*),CROSSD(*),VIS(*),RO(*),RK(*),ROMEGA(*),DISTW(*)
      REAL :: RY, F3
      LOGICAL :: TRANSL

      ISTR    = 1 
      JSTR    = IMAX + 2*IN
      KSTR    = JSTR*(JMAX + 2*JN)

      EPS     = 1.E-20
      SIGO24  = 4.*SIGOM2
      CDMAX   = 0.0
      SQRTRE  = SQRT(RE)

CC     New lo-limiter for the cross-diffusion term in ARG14 (AH 1.9.1997):
C      FREECO = RKLIM/CHLREF**2
C      CDMIN  = 5E6*FREECO
CC     This doesn't wor always properly, let's have something else... 

C     Again new lo-limiter for the cross-diffusion term in ARG14. 
C     The maximum value of CROSSD is first sought in order to find
C     the order of magnitude (AH 5.9.1997):
      LBEG = KN*KSTR + JN*JSTR + IN*ISTR + 1
      LEND = (KMAX-1+KN)*KSTR + (JMAX-1+JN)*JSTR + (IMAX-1+IN)*ISTR + 1
      DO  L = LBEG, LEND
         CDMAX = MAX(CROSSD(L),CDMAX)
      END DO
      CDMIN = 1E-8*CDMAX
      DO K  = 1,KMAX
         IA = (K+KN-1)*KSTR + JN*JSTR + IN*ISTR
         DO IJ = 1,IMAX*JMAX
            J1= (IJ-1)/IMAX       ! J1 = J-1
            I = IJ - J1*IMAX
            L = IA + J1*JSTR + I

            PROL    = 1./RO(L)
            RKL     = MAX(RK(L)*PROL, 0.0)	! Matti Palin 25.09.2017 - lisasin MAX()
            OML     = ROMEGA(L)*PROL
            DIS2    = DISTW(L)**2
            CDKOM   = MAX(CROSSD(L),CDMIN) 
            ARG11   = SQRT(RKL)/(MAX(EPS,BSTAR*OML*DISTW(L)))
            ARG12   = 500.*VIS(L)/(MAX(EPS,ROMEGA(L)*DIS2))
            ARG13   = MAX(ARG11,ARG12)
            ARG14   = SIGO24*RK(L)/(MAX(EPS,CDKOM*DIS2))
            ARG1    = MIN(ARG13,ARG14)
            FUN1(L) = TANH(ARG1**4)

            IF(TRANSL) THEN  ! Intermittency variables
              RY      = RO(L)*DISTW(L)*SQRT(RK(L))/VIS(L)
              F3      = EXP(-(RY/120.)**8)
              FUN1(L) = MAX( FUN1(L), F3 )
            ENDIF
            DISLIM  = .40559E-7*SQRTRE ! Based on the Hornet leading edge
            IF(DISTW(L) < DISLIM) FUN1(L) = 1. ! Otherwise fails
         END DO
      END DO

      RETURN
      END SUBROUTINE FUN1KO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE EPS2KO(EPS2,VIST,VIS,OHMI,RO,RK,ROMEGA,DISTW,
     &     TURLIM,FRSMUT,IMAX,JMAX,KMAX,IN,JN,KN,
     &     PR,PRT,BSTAR,A1KLEB,KOVER,RKSI)

C     Calculate the eddy viscosity of Menter's sst model at the 
C     finest grid level.
C                                               AH 31.1.1997
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,
     +           IL,K,IA,IJ,J1,I,L,KOVER,mmm,nnn,lll

      REAL :: EPS,SST,BSL,PROL,RKL,DIS2,ARG21,AUXVI,ARG22,RMUBAS,OML,
     +        BSTAR,ARG2,ARG4,FUN2,FUN3,A1KLEB,RMUTRB,RMUMEN,FRSMUT,
     +        PR,PRT,TURLIM,RMUTPV

      REAL :: EPS2(*),VIS(*),OHMI(*),RO(*),RK(*),ROMEGA(*),
     &        DISTW(*),VIST(*),RKSI(*)

      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      IL     = ISTRID*JSTRID

      EPS    = 1.E-20
cc      PRS    = PR/PRT
cc      EP4LIM = 1.0 + (TURLIM -1.0)*PRS

      IF (KOVER < 5) THEN
C     SST-model selected:
         SST = 1.0
         BSL = 0.0
      ELSE
C     BSL-model selected:
         SST = 0.0
         BSL = 1.0
      END IF

C     Update eddy viscosities at the finest grid level.
C     FUN2 makes the SST-limitation passive in free 
C     shear flows and FUN3 makes it passive in viscous sublayer.


      DO K  = 1,KMAX
         IA      = (K+KN-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX       ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*ISTRID + I

            IF(RKSI(L) <= 1) THEN ! Added for Chimera 7.1.2011
c            call ijkpai(l,imax,jmax,kmax,mmm,nnn,lll)
            PROL    = 1./RO(L)
            RKL     = MAX(RK(L)*PROL,0.0)	! Matti Palin 25.09.2017 - lisasin MAX()
            OML     = ROMEGA(L)*PROL
            DIS2    = DISTW(L)**2

            ARG21   = 2.*SQRT(RKL)/(MAX(EPS,BSTAR*OML*DISTW(L)))
            AUXVI   = VIS(L)/(MAX(EPS,ROMEGA(L)*DIS2))
            ARG22   = 500.*AUXVI
            ARG2    = MAX(ARG21,ARG22)
            ARG4    = 150.*AUXVI

            FUN2    = TANH(ARG2**2)
            FUN3    = 1. - TANH(ARG4**4)
  
            RMUBAS  = RK(L)/OML 
            RMUMEN  = A1KLEB*RK(L)/(MAX(EPS,OHMI(L)*FUN2*FUN3))
            RMUTRB  = MIN(RMUBAS,RMUMEN)
            RMUTRB  = RMUTRB*SST + RMUBAS*BSL
C           This is an important lo-limiter  AH 5.5.1997
            RMUTRB  = MAX(RMUTRB,FRSMUT)            
            RMUTRB  = MIN(RMUTRB,TURLIM*VIS(L))            

            VIST(L)  = RMUTRB
            RMUTPV  = RMUTRB/VIS(L)
            EPS2(L) = 1.0 + RMUTPV

            ENDIF
              
         END DO
      END DO

      RETURN
      END SUBROUTINE EPS2KO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE EPS2MF(EPS2,VIST,VIS,OHMI,RO,RK,ROMEGA,DISTW,FMF,VAR,
     &     TURLIM,FRSMUT,IMAX,JMAX,KMAX,IN,JN,KN,
     &     PR,PRT,BSTAR,A1KLEB,KOVER,RKSI)

C     Calculate the bubble induced turbulence (Sato model)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : TWO_FLUIDL
      USE CONSTANTS, ONLY : EPS10


      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,
     +           IL,K,IA,IJ,J1,I,L,KOVER,ii,jj,kk

      REAL :: EPS,SST,BSL,PROL,RKL,DIS2,ARG21,AUXVI,ARG22,RMUBAS,
     +        OML,BSTAR,ARG2,ARG4,FUN2,FUN3,A1KLEB,RMUTRB,RMUMEN,
     +        FRSMUT,PR,PRT,TURLIM,RMUTPV,CS,ALFAG,DBUBBLE,UR,RMUS

      REAL :: EPS2(*),VIS(*),OHMI(*),RO(*),RK(*),ROMEGA(*),
     &        DISTW(*),VIST(*),RKSI(*),FMF(*)

      TYPE(MPHASE_VARIABLES) VAR(*)

      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      IL     = ISTRID*JSTRID

      EPS    = 1.E-20
      CS     = 0.6

      DO K  = 1,KMAX
         IA      = (K+KN-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX       ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*ISTRID + I

            IF(RKSI(L) <= 1) THEN ! Added for Chimera 7.1.2011
            ALFAG   = MIN(VAR(L)%ALFA(1),VAR(L)%ALFA(2))
*           DBUBBLE = VAR(L)%BUBDIA(1)                       ! 2012 version
            DBUBBLE = MIN(VAR(L)%BUBDIA(1),VAR(L)%BUBDIA(2)) ! Currently
            DBUBBLE = MAX(DBUBBLE,EPS10) ! 

            IF(TWO_FLUIDL) THEN ! Velocity difference
               UR  = (VAR(L)%U(1) - VAR(L)%U(2))**2 +
     +               (VAR(L)%V(1) - VAR(L)%V(2))**2 +
     +               (VAR(L)%W(1) - VAR(L)%W(2))**2
               UR  = SQRT(UR) + 1.E-3
            ELSE ! By tuning for cavitation
               UR      = 0.2
               DBUBBLE = 1.E-3
            ENDIF ! TWO_FLUIDL

            RMUS    = CS*ALFAG*MIN(DISTW(L),DBUBBLE)*UR*RO(L)*FMF(L)
C           CALL IJKPAI(L,IMAX,JMAX,KMAX,II,JJ,KK)           
c       write(6666,*) ii,jj,rmus,fmf(l)

C           This is an important lo-limiter  AH 5.5.1997
            RMUTRB  = VIST(L) + RMUS
            RMUTRB  = MAX(RMUTRB,FRSMUT)            
            RMUTRB  = MIN(RMUTRB,TURLIM*VIS(L))            

            VIST(L) = RMUTRB
            RMUTPV  = RMUTRB/VIS(L)

            EPS2(L) = 1.0 + RMUTPV

            ENDIF
              
         END DO
      END DO

      RETURN
      END SUBROUTINE EPS2MF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE SOURKO(SRK,SOMEGA,PTUR,EPS2,VIS,RO,RK,ROMEGA,FUN1,
     &     RICH,D1,D2,D3,NGL,IMAX,JMAX,KMAX,IN,JN,KN,KOVER,
     &     BSTAR,CKAPPA,BETA1,SIGOM1,BETA2,SIGOM2,ZETA2,CMU,
     &     TURDESL,TTS,TUR_FRS_SOURCE,TRANSL,TRM,T2,STRAIN,
     &     VELLAP,QSAS,OHMI,DISTW)

C     Calculate the source terms for k and omega in Menter's SST model.
C                                                        AH 31.1.1997
      USE CONSTANTS, ONLY   : EPS10
      USE MAIN_ARRAYS, ONLY : BLKS
      USE TYPE_ARRAYS

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IMAX,JMAX,KMAX,IN,JN,KN,KOVER
      REAL, INTENT(IN)    :: BSTAR,CKAPPA,BETA1,SIGOM1
      REAL, INTENT(IN)    :: BETA2,SIGOM2,ZETA2,CMU
      REAL, DIMENSION(*)  :: SRK,SOMEGA,PTUR,VIS,RO,RK,ROMEGA,EPS2,QSAS
      REAL, DIMENSION(*)  :: FUN1,RICH,D1,D2,D3,TTS,T2,STRAIN,VELLAP
      REAL, DIMENSION(*)  :: OHMI,DISTW

      REAL, PARAMETER :: BERO=3.6, EPS=1.0E-10
      REAL, PARAMETER :: cdes1=0.78, cdes2=0.61

      REAL :: PBSTAR,PERSBS,ROTCO,ONEMF1,BETA,SIGOM,GAMKOM
      REAL :: RMUTRB,FCO,OMEPRO,OMEDIS,LVK3D,FSAS,T1,TURLEN,QSASL
      REAL :: prol,cdes,stke,ransle,lesle,desle,RKLIM,EPSLIM,TRUEOMEGA
      REAL :: RD,FD

      INTEGER :: ISTRID,JSTRID,IL,K,IA,IJ,J1,I,L,NGL

      LOGICAL :: TURDESL,TUR_FRS_SOURCE,TRANSL

      TYPE(INTERMITTENCY) :: TRM(*)

c      TYPE(BLOCKS) :: BLKS(*)

      PBSTAR = 1.0/BSTAR
      PERSBS = 1.0/SQRT(BSTAR)
      FSAS   = 1.25 ! For tuning the SST-SAS model
      FSAS   = 1.00 ! vmv
      
C     The rotation and curvature correction  (Added Nov 27 1997 and
C     updated Mar 12, 1998 by AH):
      IF ((KOVER == 0).OR.(KOVER == 5) .AND..NOT.TURDESL) THEN
         ROTCO = 0.0
      ELSE IF((KOVER == 1).OR.(KOVER == 6) .AND. .NOT.TURDESL) THEN
         ROTCO = 1.0
      ELSE IF(TURDESL) THEN ! No rotational correction in DES
         ROTCO = 0.
      END IF

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

C ... Free-stream values

      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      QSASL   = 0.

      DO K  = 1,KMAX
         IA = (K+KN-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX            ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*ISTRID + I

C           Blended model coefficients:
            ONEMF1    = MAX(0.0,(1.0-FUN1(L)))  
            BETA      = FUN1(L)*BETA1  + ONEMF1*BETA2
            SIGOM     = FUN1(L)*SIGOM1 + ONEMF1*SIGOM2          
            GAMKOM    = BETA*PBSTAR - SIGOM*CKAPPA**2*PERSBS
            cdes      = fun1(l)*cdes1  + onemf1*cdes2

C           Eddy viscosity:
            RMUTRB    = MAX(EPS,EPS2(L)-1.0)*VIS(L)

C           The rotation correction  (Added Nov 27 1997 by AH):
            FCO       = MAX(0.0,(1.0 - ROTCO)) 
     &                + MAX(0.0,ROTCO/(1.0 + BERO*RICH(L)))

            IF(TUR_FRS_SOURCE) THEN ! Preserve free-stream values
               TRUEOMEGA = MAX(0.,ROMEGA(L)-EPSLIM)
            ELSE
               TRUEOMEGA = ROMEGA(L)
            ENDIF

            PROL   = 1.0/RO(L)

C           Source term additions:

            IF(.NOT.TURDESL) THEN ! RANS source terms  

               QSASL = 0.      ! This line should be useless 
               SOMEGA(L) = ONEMF1*SOMEGA(L)

               IF(TRANSL) THEN          ! Intermittency variables
*                  SRK(L)    = TRM(L)%G*PTUR(L) - 
*     &          MIN(MAX(TRM(L)%G,0.1),1.0)*BSTAR*RK(L)*TRUEOMEGA*PROL
                  SRK(L)    = TRM(L)%GEFF*PTUR(L) - 
     &          MIN(MAX(TRM(L)%GEFF,0.1),1.0)*BSTAR*RK(L)*TRUEOMEGA*PROL
               ELSE
                  SRK(L)    = PTUR(L) - BSTAR*RK(L)*TRUEOMEGA*PROL
               ENDIF 

            ELSE  ! DES or SAS sources  

C     DES, length scales, AH 25.2.2008.
            STKE   = SQRT(MAX(PROL*RK(L),0.0))		! Matti Palin 25.09.2017 - lisasin MAX()
            RANSLE = STKE/(BSTAR*PROL*ROMEGA(L)) + 1.0E-7
            LESLE  = MAX(D1(L),D2(L),D3(L))
            IF(KOVER == 0 .OR. KOVER == 5) THEN ! The basic model
              LESLE  = CDES*LESLE + 1.0E-7
              DESLE  = MIN(RANSLE,LESLE)
              SRK(L) = PTUR(L) - RK(L)*STKE/DESLE
            ELSEIF(KOVER == 2 .OR. KOVER == 7) THEN ! DDES
              LESLE = RANSLE-CDES*LESLE
C              DESLE = RANSLE - FUN1(L)1.414213*MAX(0.,LESLE)
              RD    = EPS2(L)*VIS(L)*1.414214 /
     &       (RO(L)*(CKAPPA*DISTW(L))**2*SQRT(STRAIN(L)**2+OHMI(L)**2))
              FD    = 1 - TANH((20.*RD)**3)
              DESLE = RANSLE - FD*MAX(0.,LESLE)
              SRK(L)= PTUR(L) - RK(L)*STKE/DESLE
            ELSEIF(KOVER == 3 .OR. KOVER == 8) THEN ! SST-SAS (8=BSL_SAS?!)
c vmv         LVK3D = CKAPPA*STRAIN(L)/(VELLAP(L)+EPS10)
              LVK3D = CKAPPA*ABS(STRAIN(L)/(VELLAP(L)+EPS10)) ! vmv
              TURLEN= SQRT(RK(L))/(CMU*ROMEGA(L)+EPS10)
c             T1    = ZETA2*CKAPPA*STRAIN(L)**2*(TURLEN/LVK3D)**2
c ...          vmv 4.3.2018:
              T1    = RO(L)*ZETA2*STRAIN(L)**2*(TURLEN/LVK3D)**2
              QSASL = FSAS*RO(L)*MAX(0.,T1-T2(L))
              SOMEGA(L) = ONEMF1*SOMEGA(L) + QSASL
              QSAS(L)   = QSASL
            ELSE
              STOP 'KOVER in SOURKO'
            ENDIF

C ... Store length scale as 'time-scale'
            TTS(L)  = DESLE ! RANS TTS is calculated in subroutine TTIMES
C     DES-dissipation, AH 25.2.2008.

C           
            IF(TRANSL) THEN ! Intermittency variables (not tested)
*              SRK(L)    = TRM(L)%G*PTUR(L) - MIN(MAX(TRM(L)%G,0.1),1.0)
*     &                                  *RK(L)*STKE/DESLE
              SRK(L)    = TRM(L)%GEFF*PTUR(L) - 
     &                    MIN(MAX(TRM(L)%GEFF,0.1),1.0)*RK(L)*STKE/DESLE
            ENDIF           ! Intermittency in DES

            ENDIF ! RANS or DES

            OMEPRO    = GAMKOM*RO(L)*PTUR(L)/(RMUTRB*ROMEGA(L))
            OMEDIS    = FCO*BETA*PROL*TRUEOMEGA !ROMEGA(L)
            SOMEGA(L) = SOMEGA(L) + ROMEGA(L)*(OMEPRO - OMEDIS)

         END DO
      END DO
      
      RETURN
      END SUBROUTINE SOURKO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE CROSKO(CROSSD,RO,RK,ROMEGA,VOL,
     &     A1,A1X,A1Y,A1Z,A2,A2X,A2Y,A2Z,A3,A3X,A3Y,A3Z,
     &     IMAX,JMAX,KMAX,IN,JN,KN,SIGOM2,SIGPHI,T2,KOVER)
C
C     Computes the cross-diffusion term of the k-omega SST-model. 
C     (ITURB=6).
C                                           AH 7.2.1997
C                                                     
      USE CONSTANTS, ONLY : EPS10
      REAL :: CROSSD(*),RO(*),RK(*),ROMEGA(*),VOL(*),T2(*),
     &        A1(*),A1X(*),A1Y(*),A1Z(*),
     &        A2(*),A2X(*),A2Y(*),A2Z(*),
     &        A3(*),A3X(*),A3Y(*),A3Z(*)

      SIGO22  = 2.0*SIGOM2

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
      NTOT    = IL*KSTRID

      ISTR    = 1
      JSTR    = ISTRID
      KSTR    = IL

C     Convert ro*k and ro*omega to k and omega temporarily.
      DO  L = 1, NTOT
         PROL      = 1./RO(L)
         RK(L)     = RK(L)*PROL
         ROMEGA(L) = ROMEGA(L)*PROL
      END DO

      DO K  = 1,KMAX
         IA = (K+KN-1)*KSTR + JN*JSTR + IN*ISTR + 1
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX            ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*JSTR + (I-1)*ISTR

C           Xi-direction:

            RK1M = A1(L)     *(RK(L)     + RK(L-ISTR))
            RK1P = A1(L+ISTR)*(RK(L)     + RK(L+ISTR))
            OM1M = A1(L)     *(ROMEGA(L) + ROMEGA(L-ISTR))
            OM1P = A1(L+ISTR)*(ROMEGA(L) + ROMEGA(L+ISTR))

C           Eta-direction:

            RK2M = A2(L)     *(RK(L)     + RK(L-JSTR))
            RK2P = A2(L+JSTR)*(RK(L)     + RK(L+JSTR))
            OM2M = A2(L)     *(ROMEGA(L) + ROMEGA(L-JSTR))
            OM2P = A2(L+JSTR)*(ROMEGA(L) + ROMEGA(L+JSTR))

C           Zeta-direction:

            RK3M = A3(L)     *(RK(L)     + RK(L-KSTR))
            RK3P = A3(L+KSTR)*(RK(L)     + RK(L+KSTR))
            OM3M = A3(L)     *(ROMEGA(L) + ROMEGA(L-KSTR))
            OM3P = A3(L+KSTR)*(ROMEGA(L) + ROMEGA(L+KSTR))

            PVOL  = 0.5/VOL(L)

            DKDX = PVOL*(
     &             A3X(L+KSTR)*RK3P - A3X(L)*RK3M
     &           + A2X(L+JSTR)*RK2P - A2X(L)*RK2M
     &           + A1X(L+ISTR)*RK1P - A1X(L)*RK1M)
            DKDY = PVOL*(
     &             A3Y(L+KSTR)*RK3P - A3Y(L)*RK3M
     &           + A2Y(L+JSTR)*RK2P - A2Y(L)*RK2M
     &           + A1Y(L+ISTR)*RK1P - A1Y(L)*RK1M)
            DKDZ = PVOL*(
     &             A3Z(L+KSTR)*RK3P - A3Z(L)*RK3M
     &           + A2Z(L+JSTR)*RK2P - A2Z(L)*RK2M
     &           + A1Z(L+ISTR)*RK1P - A1Z(L)*RK1M)

            DODX = PVOL*(
     &             A3X(L+KSTR)*OM3P - A3X(L)*OM3M
     &           + A2X(L+JSTR)*OM2P - A2X(L)*OM2M
     &           + A1X(L+ISTR)*OM1P - A1X(L)*OM1M)
            DODY = PVOL*(
     &             A3Y(L+KSTR)*OM3P - A3Y(L)*OM3M
     &           + A2Y(L+JSTR)*OM2P - A2Y(L)*OM2M
     &           + A1Y(L+ISTR)*OM1P - A1Y(L)*OM1M)
            DODZ = PVOL*(
     &             A3Z(L+KSTR)*OM3P - A3Z(L)*OM3M
     &           + A2Z(L+JSTR)*OM2P - A2Z(L)*OM2M
     &           + A1Z(L+ISTR)*OM1P - A1Z(L)*OM1M)

            DIVDOT    = DKDX*DODX + DKDY*DODY + DKDZ*DODZ
            DIVDK     = DKDX*DKDX + DKDY*DKDY + DKDZ*DKDZ
            DIVDOM    = DODX*DODX + DODY*DODY + DODZ*DODZ
            CROSSD(L) = SIGO22*RO(L)*DIVDOT/ROMEGA(L)

c            IF(KOVER == 3 .OR. KOVER == 8) THEN ! SST-SAS source
               T2(L) = 2.*RK(L)/SIGPHI * 
     &         MAX(DIVDK/(RK(L)**2+EPS10), DIVDOM/(ROMEGA(L)**2+EPS10))
c            ENDIF

         END DO
      END DO

C     Convert k and omega back to ro*k and ro*omega:
      DO L = 1, NTOT
         RK(L)     = RK(L)*RO(L)
         ROMEGA(L) = ROMEGA(L)*RO(L)
      END DO

      RETURN
      END SUBROUTINE CROSKO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE MAKE_EPS(RO,RK,REPS,SRK,EPSTRU,ITURB,NTOT)

      INTEGER :: ITURB,NTOT,L
      REAL , DIMENSION(NTOT) :: RO,RK,REPS,SRK,EPSTRU

C     Create true-epsilon based on srk and epsilon-tilde or k and omega

      IF (ITURB == 6) THEN
         EPSTRU = 0.09*REPS*RK/RO
      ELSE IF (ITURB == 3) THEN
         EPSTRU = REPS + SRK
      ELSE
         EPSTRU = REPS  
      END IF

      RETURN
      END SUBROUTINE MAKE_EPS


      SUBROUTINE REFLKO(RK,ROMEGA,RO,SIJ,VIS,RBK,IMAX,JMAX,KMAX,
     &     IN,JN,KN,RSURF,ICON,IHF,NPATCH,IB,MAX11,DISTW)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*),NPATCH,IB,IHF(*)
      REAL    :: RK(*),ROMEGA(*),RO(*),VIS(*),SIJ(MAX11,6),
     &           RBK(*),DISTW(*)

c      REAL RSURF(IB,6)
      REAL :: RSURF(1,6) ! pinnan karheus aerokoodissa 
C
C     Reflect RK and ROMEGA behind the solid/moving, solid/rotating 
C     or solid patches.
C                                              AH  10.2.1997
C     Massively re-written 4.2.1999, PK. RBK added 26.9.2003, TSii

      ISTR   = 1
      JSTR   = IMAX + 2*IN
      KSTR   =(JMAX + 2*JN)*JSTR

      DO LP=1,NPATCH

         ITYPE  = ICON(1,LP)
         IF2    = IHF(ICON(2,LP))

         IF (ITYPE < 8 .OR. ITYPE > 10) THEN
            CYCLE
         END IF
         IW     = ICON(3,LP)
         IF ( (IW == 1) .OR. (IW == 4) ) THEN

            IMI = ICON(4,LP)
            IMA = ICON(5,LP)
            JMI = ICON(6,LP)
            JMA = ICON(7,LP)
            IMX = JMAX          ! Fastest index on the wall

            IF (IW <= 3) THEN
               KW   = 1
               IDIR = 1
            ELSE
               KW   = IMAX
               IDIR = -1
            END IF

            CALL REFL1W_KO(IMI,IMA,JMI,JMA,KW,IMX,JSTR,KSTR,ISTR,JN,KN,
     &         IN,IDIR,SIJ,VIS,RO,RK,ROMEGA,RBK(IF2),RSURF(1,IW),MAX11,
     &         DISTW)

         ELSE IF ( (IW == 2) .OR. (IW == 5) ) THEN

            JMI = ICON(4,LP)
            JMA = ICON(5,LP)
            IMI = ICON(6,LP)
            IMA = ICON(7,LP)
            IMX = KMAX          ! Fastest index on the wall

            IF (IW <= 3) THEN
               KW   = 1
               IDIR = 1
            ELSE
               KW   = JMAX
               IDIR = -1
            END IF
            
            CALL REFL1W_KO(IMI,IMA,JMI,JMA,KW,IMX,KSTR,ISTR,JSTR,KN,IN,
     &         JN,IDIR,SIJ,VIS,RO,RK,ROMEGA,RBK(IF2),RSURF(1,IW),MAX11,
     &         DISTW)

         ELSE IF ( (IW == 3) .OR. (IW == 6) ) THEN

            IMI = ICON(4,LP)
            IMA = ICON(5,LP)
            JMI = ICON(6,LP)
            JMA = ICON(7,LP)
            IMX = IMAX          ! Fastest index on the wall

            IF (IW <= 3) THEN
               KW   = 1
               IDIR = 1
            ELSE
               KW   = KMAX
               IDIR = -1
            END IF

            CALL REFL1W_KO(IMI,IMA,JMI,JMA,KW,IMX,ISTR,JSTR,KSTR,IN,JN,
     &         KN,IDIR,SIJ,VIS,RO,RK,ROMEGA,RBK(IF2),RSURF(1,IW),MAX11,
     &         DISTW)

         END IF
      END DO                    ! Patch loop.

      RETURN
      END SUBROUTINE REFLKO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C							
      SUBROUTINE REFL1W_KO(IMI,IMA,JMI,JMA,KW,IMAX,ISTR,JSTR,KSTR,
     &     IN,JN,KN,IDIR,SIJ,VIS,RO,RK,ROMEGA,RBK,RSURF,MAX11,DISTW)

      USE NS3CO, ONLY : LN
      
      REAL :: VIS(*),RO(*),RK(*),ROMEGA(*),RSURF(*),SIJ(MAX11,6),RBK(*),
     &        DISTW(*)
      REAL :: RC1,RC2,EPS
      REAL :: XMULTI(6)

      DATA XMULTI /1.,2.,2.,1.,2.,1./

      RC1  = 1.5 
      RC2  = 1.0 - RC1
      EPS  = 1E-20

      JSTRID = IMAX+2*IN        ! Increment in the slower surface direction
      L1     = IDIR*KSTR        ! K-increment to move into the block

      KK = (KW - 1 + KN)*KSTR + (JN-1)*JSTR + (IN-1)*ISTR + 1
      KJ = (JN-1)*JSTRID + IN
      DO J = JMI,JMA
         II = KK + J*JSTR
         JJ = KJ + J*JSTRID
         NR = (LN+J-JMI)*(IMA-IMI+1+2*LN) - IMI + LN + 1
         DO I = IMI,IMA
            L  = II + I*ISTR    ! Cell just above the wall
            KC = JJ + I
            NP = I  + NR                      ! Patch index with LN ghost cells

C ... Changed to simple extrapolation

            RK(L-L1)     = -RK(L) ! 2.0*RK(L)     - RK(L+L1)
            ROMEGA(L-L1) = 2.0*ROMEGA(L) - ROMEGA(L+L1)
         END DO
      END DO

      RETURN
      END SUBROUTINE REFL1W_KO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE COEFKO(ITURB,A1KLEB,BSTAR,CKAPPA,
     &     BETA1,SIGRK1,SIGOM1,BETA2,SIGRK2,SIGOM2,
     &     SIGPHI,ZETA2,KOVER)

C     Set values for the model constants of Menter's k-omega SST model.
C                                                   AH 3.2.1997

      BETA1  = 0.
      BETA2  = 0.
      BSTAR  = 0.
      CKAPPA = 0.
      SIGPHI = 1.
      ZETA2  = 0.97 ! 0.


      IF (ITURB == 6) THEN
C        Common constants + set 1:
         A1KLEB = 0.31
         BSTAR  = 0.09
         CKAPPA = 0.41
         BETA1  = 0.075
         IF (KOVER < 5) THEN
C     SST-versions:
            SIGRK1 = 0.85
         ELSE
C     BSL-versions:
            SIGRK1 = 0.5
         END IF
         SIGOM1 = 0.5
C        Set 2:
         BETA2  = 0.0828
         SIGRK2 = 1.0
         SIGOM2 = 0.856
         IF(KOVER == 3 .OR. KOVER == 8) THEN ! SAS
         SIGPHI = 1.
         ZETA2   = 0.97
         ENDIF
      END IF

      RETURN
      END SUBROUTINE COEFKO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE SPALART_ALLMARAS(EPS2,VIST,RO,VIS,REPS,CV1,IMAX,JMAX,
     +    KMAX,IN,JN,KN,TURLIM,RKSI)

C     Calculate the turbulent viscosity for the Spalart-Allmaras model
C                                                        TSii 1.7.2009
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,IL,IA,IJ,J1,I,K,L

      REAL :: FV1,CV1,RKHI,TURLIM
      REAL :: EPS2(*),VIST(*),VIS(*),RO(*),REPS(*),RKSI(*)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      DO K  = 1,KMAX
         IA = (K+KN-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX            ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*ISTRID + I
            IF(RKSI(L) <= 1) THEN ! Added for Chimera 7.1.2011
            RKHI    = REPS(L)/VIS(L)
            FV1     = RKHI**3/(RKHI**3+CV1**3)
            VIST(L) = REPS(L)*FV1
            VIST(L) = MIN(TURLIM*VIS(L),VIST(L))
            EPS2(L) = 1.+ VIST(L)/VIS(L)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SPALART_ALLMARAS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE SOURSA(PTUR,RK,REPS,DDEPS,EPS2,VIST,VIS,SRK,SEPS,RO,
     + A1,A2,A3,VOL,DISTW,OHMI,S11,S12,S13,S22,S23,S33,DREPSDX,DREPSDY,
     + DREPSDZ,DEPSDX,DEPSDY,DEPSDZ,
     + CB1,CB2,SRNU,CV1,CW1,CW2,CW3,IMAX,JMAX,KMAX,IN,JN,KN,ITURB,M,
     + D1,D2,D3,TURDESL,KOVER,NGL)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,IN,JN,KN,ITURB,
     +           M,IL,K,IA,IJ,J1,I,L,KOVER,NGL

      REAL :: CB1,CB2,SRNU,CV1,CW1,CW2,CW3,EPS,STRAIN,S,RKHI,CKAPPA,FV1,
     +        FV2,SHAT,D2K2,FW,YRNU,DVELDXI,DELTA
      REAL :: LESLE ! Mersu
      REAL :: PTUR(*),RK(*),REPS(*),DDEPS(*),EPS2(*),VIST(*),VIS(*),
     +        SRK(*),SEPS(*),RO(*),A1(*),A2(*),A3(*),VOL(*),OHMI(*),
     +        S11(*),S12(*),S13(*),S22(*),S23(*),S33(*),DISTW(*),
     +        DREPSDX(*),DREPSDY(*),DREPSDZ(*),DEPSDX(*),DEPSDY(*),
     +        DEPSDZ(*),D1(*),D2(*),D3(*)

      REAL :: CW3G,G,R

      LOGICAL :: TURDESL

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      EPS     = 1.E-20
      CKAPPA  = 0.41
   
C ... COMPUTE WALL-INDEPENDENT PRODUCTION TERMS FOR SPALART_ALLMARAS MODEL
C ... TURDESL AND KOVER DETERMINE THE DES MODEL TO BE USED

      IF(ITURB == 9) THEN
      DO K  = 1,KMAX
         IA = (K+KN-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX            ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*ISTRID + I

C ... Length-scale according to RANS or DES

            IF(.NOT. TURDESL) THEN
               LESLE = DISTW(L)
            ELSE ! DEtached eddy models
               IF(KOVER == 0) THEN
               DELTA = MAX(D1(L),D2(L),D3(L),EPS)
               LESLE = MIN(DISTW(L),0.65*DELTA)
               ELSEIF(KOVER == 1) THEN ! Options
                  STOP 'KOVER in SOURSA'
               ENDIF
            ENDIF

C ... Strain
            STRAIN  = S11(L)*S11(L) + S12(L)*S12(L) + S13(L)*S13(L)
     +              + S12(L)*S12(L) + S22(L)*S22(L) + S23(L)*S23(L)
     +              + S13(L)*S13(L) + S23(L)*S23(L) + S33(L)*S33(L)
            STRAIN  = SQRT(2.*STRAIN)
            S       = OHMI(L) !+ 2.*MIN(0.,STRAIN-OHMI(L)) ! Aktivoi?
C ... Production
            D2K2    = (CKAPPA*LESLE)**2 + EPS
            RKHI    = REPS(L)/VIS(L)
            FV1     = RKHI**3/(RKHI**3+CV1**3)
            FV2     = 1.- RKHI/(1+RKHI*FV1)
            SHAT    = S + REPS(L)/RO(L)*FV2/D2K2
            PTUR(L) = CB1*SHAT*REPS(L)
C ... Derivative
            DVELDXI = DEPSDX(L)*DREPSDX(L) + DEPSDY(L)*DREPSDY(L) 
     +              + DEPSDZ(L)*DREPSDZ(L)
            DVELDXI = DVELDXI/SRNU*CB2
C ... Destruction
            R       = REPS(L)/(MAX(1.E-10,RO(L)*SHAT*D2K2))
            G       = R + CW2*(R**6-R)
            CW3G    = CW3/G
            IF(CW3G < 1.E-6) CW3G = 0.
            CW3G    = Cw3G**6
            FW      = ((1+CW3**6)/(1.+CW3G))**.166667
            YRNU    = CW1*FW/RO(L)*(REPS(L)/LESLE)**2

C ... Total sources

            SEPS(L) = PTUR(L) + DVELDXI - YRNU
            SRK(L)  = YRNU ! koe 0.

         ENDDO ; ENDDO ! Mersu

      ELSE
         WRITE(*,*) 'No such option in Spart-Allmaras. Exiting...'
         STOP
      ENDIF

      RETURN
      END SUBROUTINE SOURSA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C						
