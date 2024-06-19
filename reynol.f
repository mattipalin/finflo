C     ***************************************
C     *       REYNOLDS STRESS SUBROUTINES   *
C     ***************************************



C*********************************************************
C*                                                       *
C*    SOURSE TERMS IN THE SECOND-MOMENT CLOSURE          *
C*                                                       *
C*********************************************************


      SUBROUTINE REPROD(RO,U,V,W,RK,REPS,DDEPS,UU,UV,UW,VV,VW,WW,VIS,
     +     SDRXX,SDRXY,SDRXZ,SDRYY,SDRYZ,SDRZZ,OHMI,
     +     FWLL,D3,VOL,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IMAXRE,
     +     JMAXRE,KMAXRE,JRDIS,JRPRE,IN,JN,KN,KBOT,KTOP,KCP,LAMIN,
     +     F,Y,SOURL,ISTR,JSTR,KSTR,ZZZ)

      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     + UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),VIS(*),
     + SDRXX(*),SDRXY(*),SDRXZ(*),SDRYY(*),SDRYZ(*),SDRZZ(*),
     + D3(*),VOL(*),A3(*),A3X(*),A3Y(*),A3Z(*),FWLL(*),OHMI(*),
     + KCP(*),F(*),Y(*),ZZZ(*)

      LOGICAL SOURL


C ... THIS SUBROUTINE CALCULATES WALL CORRECTIONS FOR
C ... SHIMA'S AND SSG MODEL IS APPLIED. 7.8.97 PR

      IF(JRDIS == 1 .AND. JRPRE == 1) THEN  !SHIMAS MODEL
         CALL RESHIM(RO,U,V,W,RK,REPS,DDEPS,UU,UV,UW,VV,VW,WW,VIS,
     +        SDRXX,SDRXY,SDRXZ,SDRYY,SDRYZ,SDRZZ,OHMI,
     +        FWLL,D3,VOL,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IMAXRE,
     +        JMAXRE,KMAXRE,IN,JN,KN,KBOT,KTOP,KCP,LAMIN,
     +        F,Y,SOURL,ISTR,JSTR,KSTR)

      ELSEIF(JRDIS == 1 .OR. JRPRE == 1) THEN
         WRITE(*,*) 'Non valid compination of JRDIS and JRPRE=',
     +        JRDIS,JRPRE
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      IF(JRPRE == 2) THEN  !SSG MODEL VELOCITY-PRESSURE CORRELATION
         CALL REDSSG(RO,U,V,W,RK,REPS,DDEPS,UU,UV,UW,VV,VW,WW,VIS,
     +        SDRXX,SDRXY,SDRXZ,SDRYY,SDRYZ,SDRZZ,OHMI,
     +        FWLL,D3,VOL,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IMAXRE,
     +        JMAXRE,KMAXRE,IN,JN,KN,KBOT,KTOP,KCP,LAMIN,
     +        F,Y,SOURL,ISTR,JSTR,KSTR,ZZZ)
      ENDIF
         
      RETURN
      END


      SUBROUTINE RESHIM(RO,U,V,W,RK,REPS,DDEPS,UU,UV,UW,VV,VW,WW,VIS,
     +     SDRXX,SDRXY,SDRXZ,SDRYY,SDRYZ,SDRZZ,OHMI,
     +     FWLL,D3,VOL,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IMAXRE,
     +     JMAXRE,KMAXRE,IN,JN,KN,KBOT,KTOP,KCP,LAMIN,
     +     F,Y,SOURL,ISTR,JSTR,KSTR)

      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     + UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),VIS(*),
     + SDRXX(*),SDRXY(*),SDRXZ(*),SDRYY(*),SDRYZ(*),SDRZZ(*),
     + D3(*),VOL(*),A3(*),A3X(*),A3Y(*),A3Z(*),FWLL(*),OHMI(*),
     + KCP(*),F(*),Y(*)

      LOGICAL SOURL

C ... THIS SUBROUTINE CALCULATES WALL CORRECTIONS FOR
C ... SHIMA'S MODEL 16.8.93 PR

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = KSTR

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 1
C ... GIVE VALUES FOR CONSTANTS
      WCONST  = 0.015

C ... SET ZEROS FOR WALL CORRECTION
      MAXX1    = ISTRID*JSTRID*KSTRID
      DO 10 L = 1,MAXX1
         Y(L)  = 1000.
 10      F(L)  = 0.

C *********************************************************************
C ... DETERMINE WALL CORRECTION
      IF (KBOT /=  0) THEN
         IA        = (KN+KBOT-1)*KSTR
         DO 800 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 800 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW    = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = .5*D3(L)
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
               F(L)    = EXP(YEXP)
               SDRXX(L)= (1.-F(L))*SDRXX(L)
            ENDIF
 800     CONTINUE


         DO 2000 KG = KLOW+1,KUPP
         IA        = (KN+KG-1)*KSTR
         DO 2000 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 2000 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = Y(L-IL)  + .5*(D3(L)+D3(L-IL))
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
               F(L)    = EXP(YEXP)
               SDRXX(L)= (1.-F(L))*SDRXX(L)
            ENDIF
 2000    CONTINUE
      ENDIF
C ********************************************************************
C ... BOUNDARY NODE  REFLECTED SIDE .....                            *
C ********************************************************************
      IF (KTOP /= 0) THEN
         IA        = (KTOP-1+KN)*KSTR
         DO 2020 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 2020 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = .5*D3(L)
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
               F(L)    = EXP(YEXP)
               SDRXX(L)= (1.-F(L))*SDRXX(L)
            ENDIF
 2020    CONTINUE

         DO 3050 KG = KTOP-1,KUPP+1,-1
         IA        = (KN+KG-1)*KSTR
         DO 3050 J  = 1,JMAX
         IA2       = IA + (J-1+JN)*JSTR
         DO 3050 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+i + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = Y(L+IL)  + .5*(D3(L)+D3(L+IL))
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
               F(L)    = EXP(YEXP)
               SDRXX(L)= (1.-F(L))*SDRXX(L)
            ENDIF
 3050    CONTINUE
      ENDIF
C ************** WALL CORRECTION CALCULATED ***************************
C *********************************************************************


      RETURN
      END


      SUBROUTINE REDSSG(RO,U,V,W,RK,REPS,DDEPS,UU,UV,UW,VV,VW,WW,VIS,
     +     SDRXX,SDRXY,SDRXZ,SDRYY,SDRYZ,SDRZZ,OHMI,
     +     FWLL,D3,VOL,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,IMAXRE,
     +     JMAXRE,KMAXRE,IN,JN,KN,KBOT,KTOP,KCP,LAMIN,
     +     F,Y,SOURL,ISTR,JSTR,KSTR,SCALE)

      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     + UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),VIS(*),
     + SDRXX(*),SDRXY(*),SDRXZ(*),SDRYY(*),SDRYZ(*),SDRZZ(*),
     + D3(*),VOL(*),A3(*),A3X(*),A3Y(*),A3Z(*),FWLL(*),OHMI(*),
     + KCP(*),F(*),Y(*)

      REAL SCALE(*)
      LOGICAL SOURL

C ... THIS SUBROUTINE CALCULATES THE PRESSURE-VELOCITY CORRELATION 
C ... TERM IN THE REYNOLDS-STRESS MODEL
C ... AND WALL CORRECTIONS IF
C ... SSG MODEL IS USED 4.4.1996 PR
C ... LOW-REYNOLDS NUMBER REGION IS BASED ON SHIMAS MODEL

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = KSTR

      KLOW    = MAX0(KBOT,1)
      KUPP    = KMAX
      IF(KTOP /= 0) KUPP = KTOP/2
      IF(KBOT == 0) KUPP = 1

C ... SHIMA
      WCONST  = 0.015
      FDIV    = 100.


C ... proposed by Aksoy et al.
c      WCONST  = 0.015

      RC1     = 1.0
      RC2     = 0.0

C ... SET ZEROS FOR WALL CORRECTION
      MAXX1    = ISTRID*JSTRID*KSTRID
      DO 10 L = 1,MAXX1
         Y(L)  = 1000.
 10      F(L)  = 0.

         NP = 4
C *********************************************************************
C ... DETERMINE WALL CORRECTION
      IF (KBOT /=  0) THEN
         IA        = (KN+KBOT-1)*KSTR
         DO 800 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 800 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW    = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = .5*D3(L)
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
c               YEXP    = -(RK(L)**2/REPS(L)/VIS(L)/FDIV)**2
               F(L)    = EXP(YEXP)
               OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L+IL))
               SCALE(IGW)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY
               YPLUS   = SCALE(IGW)*Y(L)
               SDRXZ(L) = YPLUS/((YPLUS/SDRXZ(L))**NP+1.)**(1./NP)
            ENDIF
 800     CONTINUE


         DO 2000 KG = KLOW+1,KUPP
         IA        = (KN+KG-1)*KSTR
         DO 2000 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 2000 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = Y(L-IL)  + .5*(D3(L)+D3(L-IL))
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
c               YEXP    = -(RK(L)**2/REPS(L)/VIS(L)/FDIV)**2
               F(L)    = EXP(YEXP)
               YPLUS   = SCALE(IGW)*Y(L)
               SDRXZ(L) = YPLUS/((YPLUS/SDRXZ(L))**NP+1.)**(1./NP)

            ENDIF
 2000    CONTINUE
      ENDIF
C ********************************************************************
C ... BOUNDARY NODE  REFLECTED SIDE .....                            *
C ********************************************************************
      IF (KTOP /= 0) THEN
         IA        = (KTOP-1+KN)*KSTR
         DO 2020 J  = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 2020 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+I + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = .5*D3(L)
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
c               YEXP    = -(RK(L)**2/REPS(L)/VIS(L)/FDIV)**2
               F(L)    = EXP(YEXP)
               OHMIB    = MAX(.05,RC1*OHMI(L) - RC2*OHMI(L-IL))
               SCALE(IGW)= SQRT(RO(L)*OHMIB/VIS(L)) ! VORTICITY
               YPLUS   = SCALE(IGW)*Y(L)
               SDRXZ(L) = YPLUS/((YPLUS/SDRXZ(L))**NP+1.)**(1./NP)
            ENDIF
 2020    CONTINUE

         DO 3050 KG = KTOP-1,KUPP+1,-1
         IA        = (KN+KG-1)*KSTR
         DO 3050 J  = 1,JMAX
         IA2       = IA + (J-1+JN)*JSTR
         DO 3050 I  = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
            IGW     = IN+i + (JN+J-1)*ISTRID 
            IF (KCP(IGW) /= 0) THEN
               Y(L)    = Y(L+IL)  + .5*(D3(L)+D3(L+IL))
               YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
c               YEXP    = -(RK(L)**2/REPS(L)/VIS(L)/FDIV)**2
               F(L)    = EXP(YEXP)
               YPLUS   = SCALE(IGW)*Y(L)
               SDRXZ(L) = YPLUS/((YPLUS/SDRXZ(L))**NP+1.)**(1./NP)
            ENDIF
 3050    CONTINUE
      ENDIF
C ************** WALL CORRECTION CALCULATED ***************************
C *********************************************************************


C ... this must be commented to take every direction in to account
C ... PPR 24.6.97
C ... CALCULATE DERIVATIVES BY USING THIN LAYER APPROXIMATION IN K-
C ... DIRECTION

         DO 9000 KG = 1,KMAX
         IA        = (KN+KG-1)*KSTR
         DO 9000 J = 1,JMAX
         IA2       = IA + (JN+J-1)*JSTR
         DO 9000 I = 1,IMAX
            L      = 1 + (IN+I-1)*ISTR + IA2
C ... COEFFICIENTS. CLOSE WALL SHIMA, FAR AWAY SSG
c            YEXP    = -(WCONST*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2 
            YEXP    = -(0.011*Y(L)/VIS(L))**4*(RK(L)*RO(L))**2
C ... TASSA OLI VIRHE. ALLA OLEVA ON DIMENSIOLLINEN JA SE OLI VIRITETTY
C ... KANAVA LASKULLE. NYT KERROINTA KORJATTIIN KOSKA RO = .7361 KANAVASSA
C ... PPR 21.3.97
c            YEXP    = -(0.02*Y(L) *RO(L)/VIS(L))**4*RK(L)**2
            F1    = EXP(YEXP)
            F2    = 1. - F1 
            FWLL(L) = F2*FWLL(L)
            SDRXX(L) = SDRXX(L)*(1.-F(L))
            SDRXY(L) = Y(L)/((Y(L)/SDRXY(L))**NP+1.)**(1./NP)
               
 9000    CONTINUE
      RETURN
      END


      SUBROUTINE SOURRE(OHMI,W12,S11,RO,U,V,W,RK,REPS,DDEPS,C,VIS,
     +     SDRXX,DFI,IMAX,JMAX,KMAX,IN,JN,KN,JRDIS,JRPRE,SOURL,
     +     UU,PRO,PI,DIS,FWLL,MAXSB,MAX11)

      REAL SDRXX(MAXSB,6),DFI(MAXSB,6),UU(MAXSB,6),PRO(MAXSB,6),
     +     PI(MAXSB,6),DIS(MAXSB,6),
     +     OHMI(*),W12(MAX11,3),S11(MAX11,6),RO(*),U(*),V(*),W(*),
     +     RK(*),REPS(*),DDEPS(*),C(*),FWLL(*),VIS(*)
      LOGICAL SOURL

C ... JRDIS   DISSIPATION MODEL               EX/TL INSTALLATION DATE  BY
C     ====================================================================
C       1     SHIMAS DISSIPATION                 EX      16. 8.1993    PPR
C       2     SSG    DISSIPATION                 EX            1996    PPR
C     ====================================================================

C ... JRPRE   REDISTRIPUTION MODEL             EX/TL INSTALLATION DATE  BY
C     ====================================================================
C       1     SHIMAS REDISTRIBUTION              EX      16. 8.1993    PPR
C       2     SSG    REDISTRIBUTION              EX            1996    PPR
C     ====================================================================

C ... THIS SUBROUTINE CALCULATES THE DISSIPATION THE REDISTRIBUTION
C ... AND THE PRODUCTION TERM IN THE REYNOLDS-STRESS MODEL 16.8.93 PR
C
C ... REVISITED 28.4.1994 FOR SOME MODIFICATIONS FOR WALL TREATMENTS PPR
C ... REVISITED 7.9.97 for exact calculation of sources.

      IF(JRDIS == 1 .AND. JRPRE == 1) THEN  !SHIMAS MODEL
       CALL SOUSHI(OHMI,W12(1,1),W12(1,2),W12(1,3),S11(1,1),S11(1,2),
     +        S11(1,3),S11(1,4),S11(1,5),S11(1,6),
     +        RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +        IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +         UU(1,1) ,UU(1,2), UU(1,3), UU(1,4), UU(1,5), UU(1,6),
     +        PRO(1,1),PRO(1,2),PRO(1,3),PRO(1,4),PRO(1,5),PRO(1,6),
     +         PI(1,1) ,PI(1,2) ,PI(1,3) ,PI(1,4) ,PI(1,5) ,PI(1,6),
     +        DIS(1,1),DIS(1,2),DIS(1,3),DIS(1,4),DIS(1,5),DIS(1,6),
     +        FWLL)
      ELSEIF(JRDIS == 1 .OR. JRPRE == 1) THEN
         WRITE(*,*) 'Non valid compination of JRDIS and JRPRE=',
     +        JRDIS,JRPRE
         WRITE(*,*) 'Exiting ...'
         STOP

      ELSEIF(JRPRE == 2) THEN  !SSG MODEL VELOCITY-PRESSURE CORRELATION
       CALL SOUSSG(OHMI,W12(1,1),W12(1,2),W12(1,3),S11(1,1),S11(1,2),
     +        S11(1,3),S11(1,4),S11(1,5),S11(1,6),
     +        RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +        IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +         UU(1,1) ,UU(1,2), UU(1,3), UU(1,4), UU(1,5), UU(1,6),
     +        PRO(1,1),PRO(1,2),PRO(1,3),PRO(1,4),PRO(1,5),PRO(1,6),
     +         PI(1,1) ,PI(1,2) ,PI(1,3) ,PI(1,4) ,PI(1,5) ,PI(1,6),
     +        DIS(1,1),DIS(1,2),DIS(1,3),DIS(1,4),DIS(1,5),DIS(1,6),
     +        FWLL)
      ELSEIF(JRPRE == 3) THEN  !SSG MODEL with elliptic correction
       CALL SSGELL(OHMI,W12(1,1),W12(1,2),W12(1,3),S11(1,1),S11(1,2),
     +        S11(1,3),S11(1,4),S11(1,5),S11(1,6),
     +        RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +        IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +         UU(1,1) ,UU(1,2), UU(1,3), UU(1,4), UU(1,5), UU(1,6),
     +        PRO(1,1),PRO(1,2),PRO(1,3),PRO(1,4),PRO(1,5),PRO(1,6),
     +         PI(1,1) ,PI(1,2) ,PI(1,3) ,PI(1,4) ,PI(1,5) ,PI(1,6),
     +        DIS(1,1),DIS(1,2),DIS(1,3),DIS(1,4),DIS(1,5),DIS(1,6),
     +        FWLL)
      ELSE
         WRITE(*,*) 'Non valid compination of JRDIS and JRPRE=',
     +        JRDIS,JRPRE
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      RETURN
      END


      SUBROUTINE SOUSSG(OHMI,W12I,W13I,W23I,S11I,S12I,S13I,S22I,S23I,
     +     S33I,RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +     IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +     UU,UV,UW,VV,VW,WW,
     +     PXX,PXY,PXZ,PYY,PYZ,PZZ,
     +     PIXX,PIXY,PIXZ,PIYY,PIYZ,PIZZ,
     +     DISXX,DISXY,DISXZ,DISYY,DISYZ,DISZZ,FWLL)
      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     +     UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),SDRXX(MAXSB,6),
     +     DFI(MAXSB,6),FWLL(*),C(*),VIS(*),
     +     PXX(*),PXY(*),PXZ(*),PYY(*),PYZ(*),PZZ(*),
     +     PIXX(*), PIXY(*), PIXZ(*), PIYY(*), PIYZ(*), PIZZ(*),
     +     DISXX(*),DISXY(*),DISXZ(*),DISYY(*),DISYZ(*),DISZZ(*),
     +     OHMI(*),W12I(*),W13I(*),W23I(*),S11I(*),S12I(*),S13I(*),
     +     S22I(*),S23I(*),S33I(*)
      LOGICAL SOURL

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

C ... GIVE VALUES FOR CONSTANTS
      EPS     = 1.E-10
      R13     = 1./3.
      R23     = 2./3.
      C1S      = 3.4
      C1STS    = 1.80
      C2S      = 4.2
      C3S      = 0.8
      C3STS    = 1.30
      C4S      = 1.25
      C5S      = 0.40

C ... Launder-Reece-Rodi
      C1L      = 3.0
      C1STL    = 0.
      C2L      = 0.
      C3L      = 0.8
      C3STL    = 0.
      C4L      = 1.75
      C5L      = 1.31
C ... end


C ... omat
      ALPHA   = .2
      BETA    = .03  ! OR -.03
      GAMMA   = .08

C ... SHIMA
      ALPHA   = .45
      BETA    = .0   ! OR -.03
      GAMMA   = .04  ! 0.055 antaa parempia tuloksia???? 9.11.99
      GAMMA   = .055  ! 0.055 antaa parempia tuloksia???? 9.11.99


C ... proposed by Aksoy et al.
c      ALPHA   = .29
c      BETA    = 0.
c      GAMMA   = .065
c      WCONST  = 0.015

C ... SSG VELOCITY-PRESSURE CORRELATION (WITH SHIMA)

      IB      = (KN-1)*IL + JN*ISTRID
      DO 1200 KG = 1,KMAX
         IA      = KG*IL + IB
         DO 1200 IG = 1,JMAX*ISTRID
            L       = IA + IG
C ... switching to shimas model close to the surface
            F2   = FWLL(L)
            F1   = 1. - F2
            C1    = C1L  *F1 + C1S  *F2
            C1ST  = C1STL*F1 + C1STS*F2
            C2    = C2L  *F1 + C2S  *F2
            C3    = C3L  *F1 + C3S  *F2
            C3ST  = C3STL*F1 + C3STS*F2
            C4    = C4L  *F1 + C4S  *F2
            C5    = C5L  *F1 + C5S  *F2

            FL = 1.- SDRXX(L,1)
            YDIST = SDRXX(L,2)
            YPLUS = SDRXX(L,3)
               
            REP2 = REPS(L) + DDEPS(L)
C ... ANISOTROPY OF REYNOLDS STRESS
C ... ANISOTROPY OF REYNOLDS STRESS
            YP2K = 1./(UU(L)+VV(L)+WW(L)+EPS)
            B11 = YP2K*UU(L) - R13
            B22 = YP2K*VV(L) - R13
C            B33 = YP2K*UU(L,6) - R13 ! bad values for limiting cases
            B33 = - B11 - B22
            B12 = YP2K*UV(L)
            B13 = YP2K*UW(L)
            B23 = YP2K*VW(L)
            B21 = B12
            B31 = B13
            B32 = B23
            PIB = B11*B11 + B22*B22 + B33*B33 + 2.*(
     +            B12*B12 + B13*B13 + B23*B23)

C ... STRAIN RATE
            DIV = (S11I(L) + S22I(L) + S33I(L))*.33333333
            S11 = S11I(L) - DIV
            S22 = S22I(L) - DIV
            S33 = S33I(L) - DIV
            S12 = S12I(L)
            S13 = S13I(L)
            S23 = S23I(L)
            S21 = S12
            S31 = S13
            S32 = S23

C ... ROTATION TENSOR
            W12 = W12I(L)
            W13 = W13I(L)
            W23 = W23I(L)
            W21 = -W12
            W31 = -W13
            W32 = -W23

C ... DERIVATIVES
            DUX = S11I(L)
            DVY = S22I(L)
            DWZ = S33I(L)
            DUY = S12I(L) + W12I(L)
            DUZ = S13I(L) + W13I(L)
            DVX = S12I(L) - W12I(L)
            DVZ = S23I(L) + W23I(L)
            DWX = S13I(L) - W13I(L)
            DWY = S23I(L) - W23I(L)

C ... CALCULATE PRODUCTION TERM

            P2XX  = -2.*(UU(L)*DUX+UV(L)*DUY+UW(L)*DUZ)
            P2YY  = -2.*(UV(L)*DVX+VV(L)*DVY+VW(L)*DVZ)
            P2ZZ  = -2.*(UW(L)*DWX+VW(L)*DWY+WW(L)*DWZ)
            P2XY  = -(UU(L)*DVX+UV(L)*DVY+UW(L)*DVZ+
     +                UV(L)*DUX+VV(L)*DUY+VW(L)*DUZ)
            P2XZ  = -(UU(L)*DWX+UV(L)*DWY+UW(L)*DWZ+
     +                UW(L)*DUX+VW(L)*DUY+WW(L)*DUZ)
            P2YZ  = -(UV(L)*DWX+VV(L)*DWY+VW(L)*DWZ+
     +                UW(L)*DVX+VW(L)*DVY+WW(L)*DVZ)
            
C ... THESE TERMS ARE NEEDED IN REDISTRIBUTION TERM (PHI2)

            DXX   = -2.*(UU(L)*DUX+UV(L)*DVX+UW(L)*DWX)
            DYY   = -2.*(UV(L)*DUY+VV(L)*DVY+VW(L)*DWY)
            DZZ   = -2.*(UW(L)*DUZ+VW(L)*DVZ+WW(L)*DWZ)
            DXY   = -(UU(L)*DUY+UV(L)*DVY+UW(L)*DWY+
     +                UV(L)*DUX+VV(L)*DVX+VW(L)*DWX)
            DXZ   = -(UU(L)*DUZ+UV(L)*DVZ+UW(L)*DWZ+
     +                UW(L)*DUX+VW(L)*DVX+WW(L)*DWX)
            DYZ   = -(UV(L)*DUZ+VV(L)*DVZ+VW(L)*DWZ+
     +                UW(L)*DUY+VW(L)*DVY+WW(L)*DWY)
            
            P  =  .5*(P2XX + P2YY + P2ZZ)

C ... PHI 1 AND DISSIPATION - THE WALL CORECTION

C ... CSTAR IS not (vrt Shima) ONLY EFFECT OF WALL AND DISSIPATION ON WALL
C ... P has an effect to the CTART
            C1OVER = C1 + C1ST*MAX(P,0.)/(REPS(L)+EPS)
            CSTAR1 = (C1OVER - (C1OVER-2.)*FL)*REP2
            CSTAR2 =  C2*(1. - FL)*REPS(L)

            AXX1=-CSTAR1*B11+CSTAR2*(B11*B11+B12*B12+B13*B13 - R13*PIB)
            AYY1=-CSTAR1*B22+CSTAR2*(B12*B12+B22*B22+B23*B23 - R13*PIB)
            AZZ1=-CSTAR1*B33+CSTAR2*(B13*B13+B23*B23+B33*B33 - R13*PIB)
            AXY1=-CSTAR1*B12+CSTAR2*(B11*B12+B12*B22+B13*B23)         
            AXZ1=-CSTAR1*B13+CSTAR2*(B11*B13+B12*B23+B13*B33)          
            AYZ1=-CSTAR1*B23+CSTAR2*(B12*B13+B22*B23+B23*B33)          


C ... PHI-WALL

            CO1   = ALPHA*FL 
            CO2   = BETA*FL 
c            CO3   = 2.*GAMMA*FL*RK(L)
            CO3   = 2.*GAMMA*FL*.5*(UU(L)+VV(L)+WW(L))

            AXX2 = CO1*(P2XX-R23*P)+CO2*(DXX-R23*P)+CO3*S11
            AYY2 = CO1*(P2YY-R23*P)+CO2*(DYY-R23*P)+CO3*S22
            AZZ2 = CO1*(P2ZZ-R23*P)+CO2*(DZZ-R23*P)+CO3*S33

            AXY2 = CO1* P2XY      + CO2*DXY   +     CO3*S12
            AXZ2 = CO1* P2XZ      + CO2*DXZ   +     CO3*S13
            AYZ2 = CO1* P2YZ      + CO2*DYZ   +     CO3*S23

C ... PHI2
            BS11 = B11*S11 + B12*S12 + B13*S13
            BS12 = B11*S12 + B12*S22 + B13*S23
            BS13 = B11*S13 + B12*S23 + B13*S33
            BS21 = B12*S11 + B22*S12 + B23*S13
            BS22 = B12*S12 + B22*S22 + B23*S23
            BS23 = B12*S13 + B22*S23 + B23*S33
            BS31 = B13*S11 + B23*S12 + B33*S13
            BS32 = B13*S12 + B23*S22 + B33*S23
            BS33 = B13*S13 + B23*S23 + B33*S33

            BSSS = BS11 + BS22 + BS33

            BW11 =           B12*W12 + B13*W13
            BW12 = B11*W21 +           B13*W23
            BW13 = B11*W31 + B12*W32
            BW21 =           B22*W12 + B23*W13
            BW22 = B21*W21 +           B23*W23
            BW23 = B21*W31 + B22*W32
            BW31 =           B32*W12 + B33*W13
            BW32 = B31*W21 +           B33*W23
            BW33 = B31*W31 + B32*W32

C           CO3   = (C3-C3ST/2.*SQRT(PIB + EPS))*RK(L)
C ... VIRHE.MISTAKOHAN TULI ***
            CO3   = (C3-C3ST   *SQRT(PIB + EPS))*RK(L)
            CO4   = C4*RK(L)
            CO5   = C5*RK(L)

            AXX2 = AXX2+CO3*S11+CO4*(2.*BS11-R23*BSSS)+CO5*2.*BW11
            AYY2 = AYY2+CO3*S22+CO4*(2.*BS22-R23*BSSS)+CO5*2.*BW22
            AZZ2 = AZZ2+CO3*S33+CO4*(2.*BS33-R23*BSSS)+CO5*2.*BW33
            AXY2 = AXY2+CO3*S12+CO4*(BS12+BS21)+CO5*(BW12+BW21)
            AXZ2 = AXZ2+CO3*S13+CO4*(BS13+BS31)+CO5*(BW13+BW31)
            AYZ2 = AYZ2+CO3*S23+CO4*(BS23+BS32)+CO5*(BW23+BW32)

c            AYY2z = CO4*(2.*BS22-R23*BSSS)+CO5*2.*BW22
c            AYY22 = CO4*(2.*BS22-R23*BSSS)
c            AYY23 = CO5*2.*BW22
c            if(j == 30) write(77,'(10F11.2)') 10000.*y(l),
c     +           2.*BS22-R23*BSSS,2.*BW22,(2.*BS22-R23*BSSS)/(2.*BW22)

C ... UPDATE SOURCE TERM

            AUXXX = -R23*REP2
            SDRXX(L,1) = AXX1 + AXX2 + P2XX + AUXXX
            SDRXX(L,4) = AYY1 + AYY2 + P2YY + AUXXX
            SDRXX(L,6) = AZZ1 + AZZ2 + P2ZZ + AUXXX
            SDRXX(L,2) = AXY1 + AXY2 + P2XY
            SDRXX(L,3) = AXZ1 + AXZ2 + P2XZ
            SDRXX(L,5) = AYZ1 + AYZ2 + P2YZ
            
C ... UPDATE REDISTRIBUTION TERM. NOT USED IN CALCULATION

            CSTAR1 = C1OVER*(1. - FL)*REP2

            AXX1=-CSTAR1*B11+CSTAR2*(B11*B11+B12*B12+B13*B13 - R13*PIB)
            AYY1=-CSTAR1*B22+CSTAR2*(B12*B12+B22*B22+B23*B23 - R13*PIB)
            AZZ1=-CSTAR1*B33+CSTAR2*(B13*B13+B23*B23+B33*B33 - R13*PIB)
            AXY1=-CSTAR1*B12+CSTAR2*(B11*B12+B12*B22+B13*B23)         
            AXZ1=-CSTAR1*B13+CSTAR2*(B11*B13+B12*B23+B13*B33)          
            AYZ1=-CSTAR1*B23+CSTAR2*(B12*B13+B22*B23+B23*B33)          


            PIXX(L) = AXX1 + AXX2
            PIYY(L) = AYY1 + AYY2
            PIZZ(L) = AZZ1 + AZZ2
            PIXY(L) = AXY1 + AXY2    
            PIXZ(L) = AXZ1 + AXZ2  
            PIYZ(L) = AYZ1 + AYZ2     

C ... UPDATE DISSIPATION TERM. NOT USED IN CALCULATION
            CSTAR1 = -2.*FL*REP2

            DISXX(L) = CSTAR1*B11 + AUXXX
            DISYY(L) = CSTAR1*B22 + AUXXX
            DISZZ(L) = CSTAR1*B33 + AUXXX
            DISXY(L) = CSTAR1*B12
            DISXZ(L) = CSTAR1*B13
            DISYZ(L) = CSTAR1*B23

C ... UPDATE PRODUCTION TERM. NOT USED IN CALCULATION
            PXX(L)   = P2XX 
            PXY(L)   = P2XY 
            PXZ(L)   = P2XZ 
            PYY(L)   = P2YY 
            PYZ(L)   = P2YZ 
            PZZ(L)   = P2ZZ 

 1200 CONTINUE
      RETURN
      END



      SUBROUTINE SSGELL(OHMI,W12I,W13I,W23I,S11I,S12I,S13I,S22I,S23I,
     +     S33I,RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +     IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +     UU,UV,UW,VV,VW,WW,
     +     PXX,PXY,PXZ,PYY,PYZ,PZZ,
     +     PIXX,PIXY,PIXZ,PIYY,PIYZ,PIZZ,
     +     DISXX,DISXY,DISXZ,DISYY,DISYZ,DISZZ,FWLL)
      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     +     UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),SDRXX(MAXSB,6),
     +     DFI(MAXSB,6),FWLL(*),C(*),VIS(*),
     +     PXX(*),PXY(*),PXZ(*),PYY(*),PYZ(*),PZZ(*),
     +     PIXX(*), PIXY(*), PIXZ(*), PIYY(*), PIYZ(*), PIZZ(*),
     +     DISXX(*),DISXY(*),DISXZ(*),DISYY(*),DISYZ(*),DISZZ(*),
     +     OHMI(*),W12I(*),W13I(*),W23I(*),S11I(*),S12I(*),S13I(*),
     +     S22I(*),S23I(*),S33I(*)
      LOGICAL SOURL

C ... SSG VELOCITY-PRESSURE CORRELATION (WITH Durbin's elliptic correction)

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

C ... GIVE VALUES FOR CONSTANTS
      EPS     = 1.E-10
      R13     = 1./3.
      R23     = 2./3.
      C1      = 3.4
      C1ST    = 1.80
      C2      = 4.2
      C3      = 0.8
      C3ST    = 1.30
      C4      = 1.25
      C5      = 0.40


      IB      = (KN-1)*IL + JN*ISTRID
      DO 1200 KG = 1,KMAX
         IA      = KG*IL + IB
         DO 1200 IG = 1,JMAX*ISTRID
            L       = IA + IG

            REP2 = REPS(L) + DDEPS(L)

C ... ANISOTROPY OF REYNOLDS STRESS
            YP2K = 1./(UU(L)+VV(L)+WW(L)+EPS)
            B11 = YP2K*UU(L) - R13
            B22 = YP2K*VV(L) - R13
C            B33 = YP2K*UU(L,6) - R13 ! bad values for limiting cases
            B33 = - B11 - B22
            B12 = YP2K*UV(L)
            B13 = YP2K*UW(L)
            B23 = YP2K*VW(L)
            B21 = B12
            B31 = B13
            B32 = B23
            PIB = B11*B11 + B22*B22 + B33*B33 + 2.*(
     +            B12*B12 + B13*B13 + B23*B23)

C ... STRAIN RATE
            DIV = (S11I(L) + S22I(L) + S33I(L))*.33333333
            S11 = S11I(L) - DIV
            S22 = S22I(L) - DIV
            S33 = S33I(L) - DIV
            S12 = S12I(L)
            S13 = S13I(L)
            S23 = S23I(L)
            S21 = S12
            S31 = S13
            S32 = S23

C ... ROTATION TENSOR
            W12 = W12I(L)
            W13 = W13I(L)
            W23 = W23I(L)
            W21 = -W12
            W31 = -W13
            W32 = -W23

C ... DERIVATIVES
            DUX = S11I(L)
            DVY = S22I(L)
            DWZ = S33I(L)
            DUY = S12I(L) + W12I(L)
            DUZ = S13I(L) + W13I(L)
            DVX = S12I(L) - W12I(L)
            DVZ = S23I(L) + W23I(L)
            DWX = S13I(L) - W13I(L)
            DWY = S23I(L) - W23I(L)

C ... CALCULATE PRODUCTION TERM

            P2XX  = -2.*(UU(L)*DUX+UV(L)*DUY+UW(L)*DUZ)
            P2YY  = -2.*(UV(L)*DVX+VV(L)*DVY+VW(L)*DVZ)
            P2ZZ  = -2.*(UW(L)*DWX+VW(L)*DWY+WW(L)*DWZ)
            P2XY  = -(UU(L)*DVX+UV(L)*DVY+UW(L)*DVZ+
     +                UV(L)*DUX+VV(L)*DUY+VW(L)*DUZ)
            P2XZ  = -(UU(L)*DWX+UV(L)*DWY+UW(L)*DWZ+
     +                UW(L)*DUX+VW(L)*DUY+WW(L)*DUZ)
            P2YZ  = -(UV(L)*DWX+VV(L)*DWY+VW(L)*DWZ+
     +                UW(L)*DVX+VW(L)*DVY+WW(L)*DVZ)
            
            P  =  .5*(P2XX + P2YY + P2ZZ)

C ... PHI 1
            TIME   = SQRT((RK(L)/(REP2+EPS))**2 + 36.*VIS(L)/(REP2+EPS))
            CSTAR1 = (C1 + C1ST*P/REP2)*RK(L)/TIME
            CSTAR2 = C2*RK(L)/TIME

            AXX1=-CSTAR1*B11+CSTAR2*(B11*B11+B12*B12+B13*B13 - R13*PIB)
            AYY1=-CSTAR1*B22+CSTAR2*(B12*B12+B22*B22+B23*B23 - R13*PIB)
            AZZ1=-CSTAR1*B33+CSTAR2*(B13*B13+B23*B23+B33*B33 - R13*PIB)
            AXY1=-CSTAR1*B12+CSTAR2*(B11*B12+B12*B22+B13*B23)         
            AXZ1=-CSTAR1*B13+CSTAR2*(B11*B13+B12*B23+B13*B33)          
            AYZ1=-CSTAR1*B23+CSTAR2*(B12*B13+B22*B23+B23*B33)          

C ... PHI2
            BS11 = B11*S11 + B12*S12 + B13*S13
            BS12 = B11*S12 + B12*S22 + B13*S23
            BS13 = B11*S13 + B12*S23 + B13*S33
            BS21 = B12*S11 + B22*S12 + B23*S13
            BS22 = B12*S12 + B22*S22 + B23*S23
            BS23 = B12*S13 + B22*S23 + B23*S33
            BS31 = B13*S11 + B23*S12 + B33*S13
            BS32 = B13*S12 + B23*S22 + B33*S23
            BS33 = B13*S13 + B23*S23 + B33*S33

            BSSS = BS11 + BS22 + BS33

            BW11 =           B12*W12 + B13*W13
            BW12 = B11*W21 +           B13*W23
            BW13 = B11*W31 + B12*W32
            BW21 =           B22*W12 + B23*W13
            BW22 = B21*W21 +           B23*W23
            BW23 = B21*W31 + B22*W32
            BW31 =           B32*W12 + B33*W13
            BW32 = B31*W21 +           B33*W23
            BW33 = B31*W31 + B32*W32

            CO3   = (C3-C3ST*SQRT(PIB + EPS))*RK(L)
            CO4   = C4*RK(L)
            CO5   = C5*RK(L)

            AXX2 = CO3*S11+CO4*(2.*BS11-R23*BSSS)+CO5*2.*BW11
            AYY2 = CO3*S22+CO4*(2.*BS22-R23*BSSS)+CO5*2.*BW22
            AZZ2 = CO3*S33+CO4*(2.*BS33-R23*BSSS)+CO5*2.*BW33
            AXY2 = CO3*S12+CO4*(BS12+BS21)+CO5*(BW12+BW21)
            AXZ2 = CO3*S13+CO4*(BS13+BS31)+CO5*(BW13+BW31)
            AYZ2 = CO3*S23+CO4*(BS23+BS32)+CO5*(BW23+BW32)


c            call ijkpai(l,imax,jmax,kmax,mmm,nnn,kkk)
c            if(mmm == 40) write(87,*) nnn,b12,u(l),s12


C ... DISSIPATION
            E11 = -UU(L)/(RK(L)+EPS)*REP2
            E22 = -VV(L)/(RK(L)+EPS)*REP2
            E33 = -WW(L)/(RK(L)+EPS)*REP2
            E12 = -UV(L)/(RK(L)+EPS)*REP2
            E13 = -UW(L)/(RK(L)+EPS)*REP2
            E23 = -VW(L)/(RK(L)+EPS)*REP2

C ... UPDATE SOURCE TERM
            SDRXX(L,1) = P2XX + E11
            SDRXX(L,4) = P2YY + E22
            SDRXX(L,6) = P2ZZ + E33
            SDRXX(L,2) = P2XY + E12
            SDRXX(L,3) = P2XZ + E13
            SDRXX(L,5) = P2YZ + E23
            
C ... UPDATE REDISTRIBUTION TERM for elliptic correction

            PIXX(L) = AXX1 + AXX2
            PIYY(L) = AYY1 + AYY2
            PIZZ(L) = AZZ1 + AZZ2
            PIXY(L) = AXY1 + AXY2    
            PIXZ(L) = AXZ1 + AXZ2  
            PIYZ(L) = AYZ1 + AYZ2     

C ... UPDATE DISSIPATION TERM. NOT USED IN CALCULATION

            DISXX(L) = E11
            DISYY(L) = E22
            DISZZ(L) = E33
            DISXY(L) = E12
            DISXZ(L) = E13
            DISYZ(L) = E23

C ... UPDATE PRODUCTION TERM. NOT USED IN CALCULATION
            PXX(L)   = P2XX 
            PXY(L)   = P2XY 
            PXZ(L)   = P2XZ 
            PYY(L)   = P2YY 
            PYZ(L)   = P2YZ 
            PZZ(L)   = P2ZZ 

 1200 CONTINUE
      RETURN
      END

      SUBROUTINE ELLICO(HAT1,HAT2,HAT3,HAT4,F1RM,F1RN,F1RW,RK,REPS,
     +     DDEPS,VIS,RO,FI,SIJ,FL,A1,A2,A3,VOL,IMAX,JMAX,KMAX,IN,JN,KN,
     +     MAXSB,MAX11)

      REAL HAT1(*),HAT2(*),HAT3(*),HAT4(*),RK(*),REPS(*),VIS(*),
     +     A1(*),A2(*),A3(*),VOL(*),F1RM(*),F1RN(*),F1RW(*),RO(*),
     +     FI(MAXSB,6),FL(*),SIJ(MAX11,6),DDEPS(*)

C ... This subroutine solves coefficients for elliptic smoothing of the
C ... velocity-pressure strain term. PPR 5.9.97
C ... HAT1 is diagonal matrix    
C ... HAT2 is i-direction      F1RM is right side
C ... HAT3 is j-direction      F1RN is right side
C ... HAT4 is k-direction      F1RW is right side

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      NTOT   = ISTRID*JSTRID*KSTRID
      IL     = ISTRID*JSTRID

      R13 = 1./3.
      EPS = 1.E-8
      CL  = 0.4
      CNU = 30.
      CL  = 0.2
      CNU = 120.
      CL2  = 0.23
      CNU2 = 130.
      CL  = 0.4
      CNU = 142.
      CL  = 0.24
      CL  = 0.3
      CNU = 400.
      CNU = .24/CL*400.
      DO I = 1,NTOT
         HAT1(I) = 0.
         HAT2(I) = 0.
         HAT3(I) = 0.
         HAT4(I) = 0.
         F1RM(I) = 0.
         F1RN(I) = 0.
         F1RW(I) = 0.
      ENDDO
c      if(imax == 48) write(77,*) 
      DO K = 1,KMAX
      KA   = (KN+K-1)*IL
      DO J = 1,JMAX
      JJ   = (JN+J-1)*ISTRID + IN + KA
      DO I = 1,IMAX
         L    = JJ + I
         PVOL = 1./VOL(L)
         YP2K = 1./(FI(L,1)+FI(L,4)+FI(L,6)+EPS)
         B11 = YP2K*FI(L,1) - R13
         B22 = YP2K*FI(L,4) - R13
C     B33 = YP2K*FI(L,6) - R13 ! bad values for limiting cases
         B33 = - B11 - B22
         B12 = YP2K*FI(L,2)
         B13 = YP2K*FI(L,3)
         B23 = YP2K*FI(L,5)
         B21 = B12
         B31 = B13
         B32 = B23
         PIB = B11*B11 + B22*B22 + B33*B33 + 2.*(
     +        B12*B12 + B13*B13 + B23*B23)
         PIII=B11**3+B22**3+B33**3+3.* B11*(B12**2+B13**2)+
     +        3.* B22*(B12**2+B23**2)+3.*B33*(B13**2 + B23**2) + 
     +        6.* B12*B23*B13
         DIV  = (SIJ(L,1)+SIJ(L,4)+SIJ(L,6))/3.

         S11L = SIJ(L,1)-DIV
         S22L = SIJ(L,4)-DIV
         S33L = SIJ(L,6)-DIV
         S12L = SIJ(L,2)
         S13L = SIJ(L,3)
         S23L = SIJ(L,5)

         SL  = S11L*S11L + S12L*S12L + S13L*S13L
     +       + S12L*S12L + S22L*S22L + S23L*S23L
     +       + S13L*S13L + S23L*S23L + S33L*S33L

         REP2 = REPS(L) + DDEPS(L)
         TIME   = SQRT((RK(L)/(REP2+EPS))**2 + 36.*VIS(L)/(REP2+EPS))
C         PIB  = .03*TIME*SQRT(SL)

         XLE  = PIB*CL2**2*(RK(L)**3/REP2**2 + 
     +        CNU2**2*VIS(L)*SQRT(VIS(L)/REP2))/RO(L)
         FL(L)   = SQRT(XLE)
        
         FL1     = SQRT(PIB*4.)*CL*SQRT(RK(L)**3/REP2**2/RO(L))
         FL2     = CL*CNU*SQRT(VIS(L)*SQRT(VIS(L)/REP2)/RO(L))*
     +        EXP(-((1.+8.*PIII)/(0.1+4.*PIB))**2)
         FL(L) = sqrt(FL1**2 + FL2**2)
         if(i == 40) write(77,*) j,FL(L),FL1,FL2,SQRT(XLE)
         XLE = FL(L)**2*PVOL

         BBI  = XLE*2.*A1(L       )**2/(VOL(L) + VOL(L-1     ))
         BBJ  = XLE*2.*A2(L       )**2/(VOL(L) + VOL(L-ISTRID))
         BBK  = XLE*2.*A3(L       )**2/(VOL(L) + VOL(L-IL    ))
         AAI  = XLE*2.*A1(L+1     )**2/(VOL(L) + VOL(L+1     ))
         AAJ  = XLE*2.*A2(L+ISTRID)**2/(VOL(L) + VOL(L+ISTRID))
         AAK  = XLE*2.*A3(L+IL    )**2/(VOL(L) + VOL(L+IL    ))

c         BBI  = 0.
c         BBK  = 0.
c         AAI  = 0.
c         AAK  = 0.

         HAT2(L) = BBI
         HAT3(L) = BBJ
         HAT4(L) = BBK
         F1RM(L) = AAI
         F1RN(L) = AAJ
         F1RW(L) = AAK
         DD      = -BBI - BBJ - BBK - AAI - AAJ - AAK
         HAT1(L) = DD - 1.
c         write(*,*) i,bbi,xle,A1(L)**2/(VOL(L) + VOL(L-1))
C     This is different for  ^^  normal Poisson
      ENDDO

c      if(imax == 48) write(77,*) 
      ENDDO
      ENDDO
C ... boundary conditions must be given in different systems.


      RETURN
      END


      SUBROUTINE ELLBOU(HAT,F1RM,PI,ICM,IMAX,JMAX,KMAX,IDI1,IDI2,IDI3,
     +     IN,JN,KN,NBL,NPATCH,ICON,RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,
     +     EPS2,VOL,A1,A1X,RK,REPS,DDEPS,ITURB,ISTATE,MAXB,OHMI,U,V,W,
     +     SIJ,MAX11,FI)

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: HAT(MAXB,*),ICON(IC9,*),RO(*),P(*),PDIFF(*),
     +     RM(*),RN(*),RW(*),E(*),VIS(*),VIST(*),EPS2(*),
     +     RK(*),REPS(*),DDEPS(*),VOL(*),A1(MAXB,*),A1X(MAXB,*),
     +     PI(MAXB,6),F1RM(MAXB,*),OHMI(*),U(*),V(*),W(*),SIJ(MAX11,6),
     +     FI(MAXB,6)
C
      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      XAUX    = 0.
      EPS = 1.E-10

      UUB =  10.*7.4E-5         !tangential component
      VVB = -20.*7.4E-5         !normal component
      WWB =  10.*7.4E-5         !tangential component
      UVB =   4.*5.0E-6         !shear stress component
      UVB =   8.*1.4E-3/4.      ! B/y+, y+=4 shear stress component
      UWB =   0.*5.0E-6         !shear stress component ?????
      VWB =   0.*5.0E-6         !shear stress component ?????


C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
      IFACE   = ICON(3,IP)
      I1    = ICON(4,IP)
      I2    = ICON(5,IP)
      J1    = ICON(6,IP)
      J2    = ICON(7,IP)
      IDIR    = 1
      IST2    = 1
      IF(IFACE >= 4) IDIR = -1
C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         ID   = IDI1
         JST2 = JSTRID
         IA   = 1
         K1  = IMAX
         K2  = IMAX
C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         ISTR = 1
         JSTR = IL
         KSTR = ISTRID
         ID   = IDI2
         IST2 = KSTRID
         JST2 = 1
         IA   = 2
         K1   = JMAX
         K2   = JMAX
C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         ID   = IDI3
         JST2 = ISTRID
         IA   = 3
         K1   = KMAX
         K2   = KMAX
      ENDIF
      IF(IFACE <= 3) K1 = 1
      IF(IFACE > 3) K2 = 1

C ... END OF SET-UP OF INDECES.
         DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K1-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
C ... TAA ON TODENNAKOISESTI TURHAAN. PITAA TARKISTAA....
            HAT(LB,2)  = 0.
            HAT(LB,3)  = 0.
            HAT(LB,4)  = 0.
            F1RM(LB,1) = 0.
            F1RM(LB,2) = 0.
            F1RM(LB,3) = 0.
            PI(LB,1)   = 0.
            PI(LB,2)   = 0.
            PI(LB,3)   = 0.
            PI(LB,4)   = 0.
            PI(LB,5)   = 0.
            PI(LB,6)   = 0.
        ENDDO
         ENDDO

      IF(.NOT.(((ID /= 0 .AND. IFACE <= 3).OR.(ID > K1 .AND. 
     +     IFACE >= 4)) .AND. 
     +     (IBC >= 8 .AND. IBC <= 10))) THEN ! FOR not VISCOUS SURFACES

         
         IF(IFACE <= 3) THEN
         DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K1-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            F1RM(LB,IA)  =  HAT(L,IA+1)
            HAT(LB,1)    = -HAT(L,IA+1)
            if(iface == 5) then
            HAT(LB,1)    =  HAT(L,IA+1)
            PI(LB,1)   = -2.*PI(L,1)*HAT(LB,1)/(1.E-8+RK(L))
            PI(LB,2)   = -2.*PI(L,2)*HAT(LB,1)/(1.E-8+RK(L))
            PI(LB,3)   = -2.*PI(L,3)*HAT(LB,1)/(1.E-8+RK(L))
            PI(LB,4)   = -2.*PI(L,4)*HAT(LB,1)/(1.E-8+RK(L))
            PI(LB,5)   = -2.*PI(L,5)*HAT(LB,1)/(1.E-8+RK(L))
            PI(LB,6)   = -2.*PI(L,6)*HAT(LB,1)/(1.E-8+RK(L))
            endif
         ENDDO
         ENDDO
         ELSE
         DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K1-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            HAT(LB,IA+1) =  F1RM(L,IA)
            HAT(LB,1)    = -F1RM(L,IA)
            if(iface == 5) then
             HAT(LB,1)    =  F1RM(L,IA)
            PI(LB,1)   = 2.*PI(L,1)*HAT(LB,1) /(1.E-8+RK(L))
            PI(LB,2)   = 2.*PI(L,2)*HAT(LB,1) /(1.E-8+RK(L))
            PI(LB,3)   = 2.*PI(L,3)*HAT(LB,1) /(1.E-8+RK(L))
            PI(LB,4)   = 2.*PI(L,4)*HAT(LB,1) /(1.E-8+RK(L))
            PI(LB,5)   = 2.*PI(L,5)*HAT(LB,1) /(1.E-8+RK(L))
            PI(LB,6)   = 2.*PI(L,6)*HAT(LB,1) /(1.E-8+RK(L))
            endif
         ENDDO
         ENDDO
         ENDIF
      ELSE
         IF(IFACE <= 3) THEN
         DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K1-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            F1RM(LB,IA)  =  HAT(L,IA+1)
            HAT(LB,1)    =  HAT(L,IA+1)
            AX      = A1X(L,(IA-1)*3 + 1)
            AY      = A1X(L,(IA-1)*3 + 2)
            AZ      = A1X(L,(IA-1)*3 + 3)
            YPRO    = 1./RO(L)
            OHMIB   = OHMI(L)
            SCALE   = SQRT(RO(L)*OHMIB/VIS(L)) ! UTAU/NU
            UTAU    = SCALE*VIS(L)*YPRO
            VISB    = VIS(L)
            XNU     = VISB*YPRO

            SIG12    = -IDIR*(SIGN(.5,FI(L,2)) + SIGN(.5,FI(L,2)-EPS))
            SIG13    = -IDIR*(SIGN(.5,FI(L,3)) + SIGN(.5,FI(L,3)-EPS))
            SIG23    = -IDIR*(SIGN(.5,FI(L,5)) + SIGN(.5,FI(L,5)-EPS))

            
            SCAL  = SCALE**4*UTAU**2

            XHALF = 0.5
            XI22  = 0.
            XI12  = 0.
            YP   = 0.
            DIV  = (SIJ(L,1)+SIJ(L,4)+SIJ(L,6))/3.
            S11L = SIJ(L,1)-DIV
            S22L = SIJ(L,4)-DIV
            S33L = SIJ(L,6)-DIV
            S12L = SIJ(L,2)
            S13L = SIJ(L,3)
            S23L = SIJ(L,5)

            SL  = S11L*S11L + S12L*S12L + S13L*S13L
     +           + S12L*S12L + S22L*S22L + S23L*S23L
     +           + S13L*S13L + S23L*S23L + S33L*S33L
            SLEGHT = 6.*SQRT(VIS(L)/(SQRT(SL)*RO(L)))
            
            DO K = K1,K2
               L2  = 1 + (IN+I-1)*ISTR + (JN+J-1)*JSTR + (KN+K-1)*KSTR

               dist= (2.*VOL(L2))/(A1(L2,IA)+A1(L2+IDIR*KSTR,IA))
               YP = YP + XHALF*DIST
               REP2 = REPS(L2) + DDEPS(L2)
               VLEGHT = 1.*SQRT(VIS(L2)**3/REP2)/RO(L2)
               YPLE = SQRT(YP**2+VLEGHT)
c               YPLE = SQRT(YP**2+0.65E-4**2)
c               YPLE = SQRT(YP**2+3.5*SLEGHT**2)
               XI22 = XI22 + FI(L2,4)/YPLE**5*DIST/RO(L2)
               XI12 = XI12 + FI(L2,2)/YPLE**5*DIST/RO(L2)
               XHALF = 1.
c               sc88 = SCALE**4*UTAU**2
c               if(i == 40) write(94,313) yp*scale,xi22/sc88,
c     +            FI(L2,4)/YPLE**5*DIST/sc88/RO(L2),sqrt(vleght)*scale,
c     +              SLEGHT*scale,yple*scale
 313           format(10(1x,E10.3))
            ENDDO
            vvw4 = -20.*.28*XI22
            uvw4 = -20.*.38*XI12
c            vvw4 = -200.*.28*XI22
c            uvw4 = -200.*.38*XI12


            dist= (2.*VOL(L))/(A1(L,IA)+A1(LT,IA))
            vvw1 = -20.*fi( l,4)/(dist/2.)**4*YPRO
            vvw2 = -20.*fi(lt,4)/(dist/2.*3.)**4*YPRO
            vvw3 = -20.*(fi(lt,4)-fi( l,4))/(dist)**4*YPRO
            SUFRVV = VVW4

            uvw1 = -20.*fi( l,2)/(dist/2.)**4*YPRO
            uvw2 = -20.*fi(lt,2)/(dist/2.*3.)**4*YPRO
            uvw3 = -20.*(fi(lt,2)-fi( l,2))/(dist)**4*YPRO
            SUFRUV = UVW4
            SUFRUW = SIG13*UWB*SCAL
            SUFRVW = SIG23*VWB*SCAL

            SUFRUU = -.5* SUFRVV
            SUFRWW = -.5* SUFRVV
c            SUFRUU = 0.
c            SUFRWW = 0.
c            SUFRUV = SUFRVV/(fi(lt,4)/fi(lt,2))

            ASUR  = VISB**2/(REPS(L) + DDEPS(L))*YPRO
            TIME  = sqrt(VISB/(REPS(L) + DDEPS(L)))
c             if(i2 >= 40) then
c               write(87,*) i,UVB*SCAL,UVW1,UVW3,UVW4
c               write(89,*) i,VVB*SCAL,VVW1,VVW3,VVW4
c               write(90,*) i,1./(VVW1*ASUR),-20.*TIME
c               write(89,*) i,UUMAX,YUMAX,IUUK,FI(iumax,4)
c               write(*,*) i,SUFRVV/SUFRUV,fi( l,4)/fi( l,2),
c     +              fi(lt,4)/fi(lt,2)
c            endif
            EPSSUR  = VISB**2/(REPS(L) + DDEPS(L))*YPRO*2.*HAT(L,IA+1)
            PI(LB,1)   = SUFRUU*EPSSUR
            PI(LB,2)   = SUFRUV*EPSSUR
            PI(LB,3)   = SUFRUW*EPSSUR
            PI(LB,4)   = SUFRVV*EPSSUR
            PI(LB,5)   = SUFRVW*EPSSUR
            PI(LB,6)   = SUFRWW*EPSSUR

         ENDDO
         ENDDO
         ELSE
         DO J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K1-1)*KSTR
         DO I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            HAT(LB,IA+1) =  F1RM(L,IA)
            HAT(LB,1)    =  F1RM(L,IA)
            AX      = A1X(LB,(IA-1)*3 + 1)
            AY      = A1X(LB,(IA-1)*3 + 2)
            AZ      = A1X(LB,(IA-1)*3 + 3)
            YPRO    = 1./RO(L)

            OHMIB   = OHMI(L)
            SCALE   = SQRT(RO(L)*OHMIB/VIS(L)) ! UTAU/NU
            UTAU    = SCALE*VIS(L)*YPRO
            VISB    = VIS(L)
            XNU     = VISB*YPRO

            SIG12    = -IDIR*(SIGN(.5,FI(L,2)) + SIGN(.5,FI(L,2)-EPS))
            SIG13    = -IDIR*(SIGN(.5,FI(L,3)) + SIGN(.5,FI(L,3)-EPS))
            SIG23    = -IDIR*(SIGN(.5,FI(L,5)) + SIGN(.5,FI(L,5)-EPS))
            
            SCAL  = SCALE**4*UTAU**2

            dist= (2.*VOL(L))/(A1(L,IA)+A1(LT,IA))
            vvw1 = -5.*fi( l,4)/(dist/2.)**4*YPRO
            vvw2 = -5.*fi(lt,4)/(dist/2.*3.)**4*YPRO
            SUFRVV = MAX(VVB*SCAL,VVW1)
            uvw1 = -2.*fi( l,2)/(dist/2.)**4*YPRO
            uvw2 = -2.*fi(lt,2)/(dist/2.*3.)**4*YPRO
            SUFRUV = MIN(SIG12*UVB*SCAL,UVW1)
            SUFRUW = SIG13*UWB*SCAL
            SUFRVW = SIG23*VWB*SCAL

            SUFRUU = -.5* SUFRVV
            SUFRWW = -.5* SUFRVV


         ENDDO
         ENDDO
         ENDIF

      ENDIF
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END


      SUBROUTINE SOUSHI(OHMI,W12I,W13I,W23I,S11I,S12I,S13I,S22I,S23I,
     +     S33I,RO,U,V,W,RK,REPS,C,DDEPS,SDRXX,DFI,VIS,
     +     IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,JRDIS,JRPRE,SOURL,
     +     UU,UV,UW,VV,VW,WW,
     +     PXX,PXY,PXZ,PYY,PYZ,PZZ,
     + PIXX,PIXY,PIXZ,PIYY,PIYZ,PIZZ,
     + DISXX,DISXY,DISXZ,DISYY,DISYZ,DISZZ,FWLL)
      DIMENSION RO(*),U(*),V(*),W(*),RK(*),REPS(*),DDEPS(*),
     +     UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),SDRXX(MAXSB,6),
     +     DFI(MAXSB,6),FWLL(*),C(*),VIS(*),
     +     PXX(*),PXY(*),PXZ(*),PYY(*),PYZ(*),PZZ(*),
     +  PIXX(*), PIXY(*), PIXZ(*), PIYY(*), PIYZ(*), PIZZ(*),
     +     DISXX(*),DISXY(*),DISXZ(*),DISYY(*),DISYZ(*),DISZZ(*),
     +     OHMI(*),W12I(*),W13I(*),W23I(*),S11I(*),S12I(*),S13I(*),
     +     S22I(*),S23I(*),S33I(*)
      LOGICAL SOURL

C ... THIS SUBROUTINE CALCULATES THE DISSIPATION THE REDISTRIBUTION
C ... AND THE PRODUCTION TERM IN THE REYNOLDS-STRESS MODEL
C ... AND WALL CORRECTIONS FOR
C ... SHIMA'S MODEL IS APPLIED. 16.8.93 PR
C
C ... REVISITED 28.4.1994 FOR SOME MODIFICATIONS FOR WALL TREATMENTS PPR
C ... REVISITED 7.8.97 FOR EXACT DERIVATIVES PPR

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

C ... GIVE VALUES FOR CONSTANTS
      EPS     = 1.E-10
      R23     = 2./3.
      C1      = 1.5
      C2      = .4
      ALPHA   = .45
      BETA    = .0   ! OR -.03
      GAMMA   = .08  ! talla laskettu vanhat kanavat
      GAMMA   = .06

      CP1     = (C2+8.)/11.
      CP2     = (30.*C2-2.)/55.
      CP3     = (8.*C2-2.)/11.

      IB      = (KN-1)*IL + JN*ISTRID
      DO 1200 KG = 1,KMAX
         IA      = KG*IL + IB
         DO 1200 IG = 1,JMAX*ISTRID
            L       = IA + IG

            FL = 1.- SDRXX(L,1)

C ... CORRECT DISSIPATION TERM PPR 18.5.1994
            REP2 = REPS(L) + DDEPS(L)

C ... DERIVATIVES
            DUX = S11I(L)
            DVY = S22I(L)
            DWZ = S33I(L)
            DUY = S12I(L) + W12I(L)
            DUZ = S13I(L) + W13I(L)
            DVX = S12I(L) - W12I(L)
            DVZ = S23I(L) + W23I(L)
            DWX = S13I(L) - W13I(L)
            DWY = S23I(L) - W23I(L)

C ... CALCULATE PRODUCTION TERM

            P2XX  = -2.*(UU(L)*DUX+UV(L)*DUY+UW(L)*DUZ)
            P2YY  = -2.*(UV(L)*DVX+VV(L)*DVY+VW(L)*DVZ)
            P2ZZ  = -2.*(UW(L)*DWX+VW(L)*DWY+WW(L)*DWZ)
            P2XY  = -(UU(L)*DVX+UV(L)*DVY+UW(L)*DVZ+
     +                UV(L)*DUX+VV(L)*DUY+VW(L)*DUZ)
            P2XZ  = -(UU(L)*DWX+UV(L)*DWY+UW(L)*DWZ+
     +                UW(L)*DUX+VW(L)*DUY+WW(L)*DUZ)
            P2YZ  = -(UV(L)*DWX+VV(L)*DWY+VW(L)*DWZ+
     +                UW(L)*DVX+VW(L)*DVY+WW(L)*DVZ)
            
C ... THESE TERMS ARE NEEDED IN REDISTRIBUTION TERM (PHI2)

            DXX   = -2.*(UU(L)*DUX+UV(L)*DVX+UW(L)*DWX)
            DYY   = -2.*(UV(L)*DUY+VV(L)*DVY+VW(L)*DWY)
            DZZ   = -2.*(UW(L)*DUZ+VW(L)*DVZ+WW(L)*DWZ)
            DXY   = -(UU(L)*DUY+UV(L)*DVY+UW(L)*DWY+
     +                UV(L)*DUX+VV(L)*DVX+VW(L)*DWX)
            DXZ   = -(UU(L)*DUZ+UV(L)*DVZ+UW(L)*DWZ+
     +                UW(L)*DUX+VW(L)*DVX+WW(L)*DWX)
            DYZ   = -(UV(L)*DUZ+VV(L)*DVZ+VW(L)*DWZ+
     +                UW(L)*DUY+VW(L)*DVY+WW(L)*DWY)
            
            P  =  .5*(P2XX + P2YY + P2ZZ)

C ... PHI 1 AND DISSIPATION - THE WALL CORECTION

C ... CSTAR IS ONLY EFFECT OF WALL AND DISSIPATION ON WALL

            CSTAR1 =  (-C1*REPS(L)+REP2)*FL   + C1*REPS(L)
            CSTAR2 =  R23*FL*REPS(L)*(-C1+1.) + R23*C1*REPS(L)
            AXX1 = -CSTAR1*UU(L)/(RK(L)+EPS) + CSTAR2
            AYY1 = -CSTAR1*VV(L)/(RK(L)+EPS) + CSTAR2
            AZZ1 = -CSTAR1*WW(L)/(RK(L)+EPS) + CSTAR2
            AXY1 = -CSTAR1*UV(L)/(RK(L)+EPS)
            AXZ1 = -CSTAR1*UW(L)/(RK(L)+EPS)
            AYZ1 = -CSTAR1*VW(L)/(RK(L)+EPS)

C ... PHI2 AND PHI-WALL
            CO1   = ALPHA*FL-CP1
            CO2   = (GAMMA*FL-CP2)*RK(L)
            CO3   = BETA*FL-CP3
            AXX2 = (CO1*(P2XX - R23*P) + CO2*2.*DUX + CO3*(DXX-R23*P))
            AYY2 = (CO1*(P2YY - R23*P) + CO2*2.*DVY + CO3*(DYY-R23*P)) 
            AZZ2 = (CO1*(P2ZZ - R23*P) + CO2*2.*DWZ + CO3*(DZZ-R23*P))
               AXY2 = CO1*P2XY + CO2*(DUY + DVX) + CO3*DXY
C ... KOKEILUA
c            AXY2 = -CP1*P2XY - CP2*RK(L)*(DUY + DVX) -CP3*DXY
C ... KOKEILU BLOKKI LOPPU
            AXZ2 = CO1*P2XZ + CO2*(DUZ + DWX) + CO3*DXZ
            AYZ2 = CO1*P2YZ + CO2*(DVZ + DWY) + CO3*DYZ


C ... UPDATE SOURCE TERM
            AUXXX      = -R23*REPS(L)
            SDRXX(L,1) = AXX1+AXX2+P2XX + AUXXX
            SDRXX(L,4) = AYY1+AYY2+P2YY + AUXXX
            SDRXX(L,6) = AZZ1+AZZ2+P2ZZ + AUXXX
            SDRXX(L,2) = AXY1+AXY2+P2XY
            SDRXX(L,3) = AXZ1+AXZ2+P2XZ
            SDRXX(L,5) = AYZ1+AYZ2+P2YZ

C ... UPDATE REDISTRIBUTION TERM. NOT USED IN CALCULATION
            CSTAR1  = C1*(1.-FL)*REPS(L)/(RK(L)+EPS)
            CSTAR2  = R23*C1*(1.-FL)*REPS(L)
            PIXX(L) = -CSTAR1*UU(L) + CSTAR2 + AXX2
            PIYY(L) = -CSTAR1*VV(L) + CSTAR2 + AYY2
            PIZZ(L) = -CSTAR1*WW(L) + CSTAR2 + AZZ2
            PIXY(L) = -CSTAR1*UV(L)         + AXY2
            PIXZ(L) = -CSTAR1*UW(L)         + AXZ2
            PIYZ(L) = -CSTAR1*VW(L)         + AYZ2

C ... UPDATE DISSIPATION TERM. NOT USED IN CALCULATION
            CSTAR1   = -FL*REP2
            CSTAR2   =  R23*FL*REPS(L)
            DISXX(L) = CSTAR1*UU(L)/(RK(L)+EPS) + CSTAR2 + AUXXX
            DISYY(L) = CSTAR1*VV(L)/(RK(L)+EPS) + CSTAR2 + AUXXX
            DISZZ(L) = CSTAR1*WW(L)/(RK(L)+EPS) + CSTAR2 + AUXXX
            DISXY(L) = CSTAR1*UV(L)/(RK(L)+EPS)         
            DISXZ(L) = CSTAR1*UW(L)/(RK(L)+EPS)         
            DISYZ(L) = CSTAR1*VW(L)/(RK(L)+EPS)         

C ... UPDATE PRODUCTION TERM. NOT USED IN CALCULATION
            PXX(L)   = P2XX
            PYY(L)   = P2YY
            PZZ(L)   = P2ZZ
            PXY(L)   = P2XY
            PXZ(L)   = P2XZ
            PYZ(L)   = P2YZ

 1200 CONTINUE
      RETURN
      END




C*******************************************************************
C*                                                                 *
C*            DIFFUSION SUBROUTINES FOR RSM                        *
C*                                                                 *
C*******************************************************************

      SUBROUTINE DIFF(RO,RK,REPS,VIS,EPS2,VOL,A,AX,AY,AZ,D,
     +     RUU,RUV,RUW,RVV,RVW,RWW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     DXX,DXY,DXZ,DYY,DYZ,DZZ,VXX,VXY,VXZ,VYY,VYZ,VZZ,
     +     IDI,PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,SOURL,
     +     KSTR,IDIR,ICYCLE,ICMAX)

C     THIS SUBROUTINE DETERMINES DIFFUSION TERM IN SECOND-MOMENT CLOS.

      USE NS3CO, ONLY : IN, JN, KN

      LOGICAL SOURL,SOUR2
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A(*),AX(*),AY(*),AZ(*),D(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     FUU(*),FUV(*),FUW(*),FVV(*),FVW(*),FWW(*),
     +     DXX(*),DXY(*),DXZ(*),DYY(*),DYZ(*),DZZ(*),
     +     VXX(*),VXY(*),VXZ(*),VYY(*),VYZ(*),VZZ(*),
     +     PSIGSC(*),PSIGS2(*)

C ... JRDIF   MODEL                            EX/TL INSTALLATION DATE  BY
C     ====================================================================
C       1     SCALAR DIFFUSION MODEL (TS&PPR)    TL       5. 5.1994    PPR
C       2     Daly and Harlow (1970)             TL      13. 4.1995    PPR
C       3     Hanjalic and Launder (1972)        TL      ??????????    ???
C       4     SCALAR DIFFUSION MODEL (TS&PPR)    EX      26. 5.1995    PPR
C       5     Daly and Harlow (1970)             EX      25. 4.1995    PPR
C       6     Hanjalic and Launder (1972)        EX      ??????????    ???
C     ====================================================================

      IF (JRDIF <= 0 .OR. JRDIF > 6) THEN
         WRITE(*,*) 'Unleagal diffusion option in Sub. DIFF. '//
     +  'JRDIF = ',JRDIF
         WRITE(*,*) 'Exiting...'
         STOP
      ENDIF

      SOUR2 = .FALSE.
      IF(SOURL .AND. ICYCLE >= ICMAX) SOUR2=.TRUE.

      GOTO(1,2,3,4,5,6),JRDIF

 1    CALL DIFSCA(RO,RK,REPS,VIS,EPS2,VOL,A,AX,AY,AZ,D,
     +     RUU,RUV,RUW,RVV,RVW,RWW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     DXX,DXY,DXZ,DYY,DYZ,DZZ,VXX,VXY,VXZ,VYY,VYZ,VZZ,
     +     IDI,PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,SOUR2,
     +     KSTR,IDIR)
      GOTO 100

 2    CALL DIFDLY(RO,RK,REPS,VIS,EPS2,VOL,A,AX,AY,AZ,D,
     +     RUU,RUV,RUW,RVV,RVW,RWW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     DXX,DXY,DXZ,DYY,DYZ,DZZ,VXX,VXY,VXZ,VYY,VYZ,VZZ,
     +     IDI,PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,SOUR2,
     +     KSTR,IDIR)
      GOTO 100

 3    WRITE(*,*) 'Thinlayer diffusion by Hanjalic and Launder is not '//
     +           'yet implemented (26.4.1995)'
      WRITE(*,*) 'Exiting...'
      STOP

 4    CONTINUE
C ... IS DONE IN SUBROUTINE DIFFEX
      GOTO 100 

 5    CONTINUE
C ... IS DONE IN SUBROUTINE DIFFEX
      GOTO 100 

 6    CONTINUE
C ... IS DONE IN SUBROUTINE DIFFEX
      GOTO 100 

 100  CONTINUE
      RETURN
      END


      SUBROUTINE DIFEX(RO,RK,REPS,VIS,EPS2,VOL,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     RUU,RUV,RUW,RVV,RVW,RWW,DRUU,DUU,
     +     TUU,VUU,
     +     PSIGSC,PSIGS2,IMAX,JMAX,KMAX,IN,JN,KN,
     +     D11X,D11Y,D11Z,JRDIF,SOURL)

C ... THIS SUBROUTINE DETERMINES EXACT DIFFUSIONS IN SECOND-MOMENT CLOS.

      LOGICAL SOURL
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +     A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     DUU(*),DRUU(*),TUU(*),VUU(*),
     +     D11X(*),D11Y(*),D11Z(*)

C ... JRDIF   MODEL                            EX/TL INSTALLATION DATE  BY
C     ====================================================================
C       1     SCALAR DIFFUSION MODEL (TS&PPR)    TL       5. 5.1994    PPR
C       2     Daly and Harlow (1970)             TL      13. 4.1995    PPR
C       3     Hanjalic and Launder (1972)        TL      ??????????    ???
C       4     SCALAR DIFFUSION MODEL (TS&PPR)    EX      26. 4.1995    PPR
C       5     Daly and Harlow (1970)             EX      25. 4.1995    PPR
C       6     Hanjalic and Launder (1972)        EX      ??????????    ???
C     ====================================================================

      IF (JRDIF <= 0 .OR. JRDIF > 6) THEN
         WRITE(*,*) 'Unleagal diffusion option in Sub. DIFF. '//
     +  'JRDIF = ',JRDIF
         WRITE(*,*) 'Exiting...'
         STOP
      ENDIF
      GOTO(1,2,3,4,5,6),JRDIF

 1    CONTINUE
C ... IS DONE IN SUBROUTINE DIFF
      GOTO 100 

 2    CONTINUE
C ... IS DONE IN SUBROUTINE DIFF
      GOTO 100 

 3    CONTINUE
C ... IS DONE IN SUBROUTINE DIFF
      GOTO 100 

 4    CALL DIFEXS(RO,RK,REPS,VIS,EPS2,VOL,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     RUU,RUV,RUW,RVV,RVW,RWW,DRUU,DUU,
     +     TUU,VUU,
     +     PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,
     +     D11X,D11Y,D11Z,SOURL)
      GOTO 100

 5    CALL DIFEXH(RO,RK,REPS,VIS,EPS2,VOL,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     RUU,RUV,RUW,RVV,RVW,RWW,DRUU,DUU,
     +     TUU,VUU,
     +     PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,
     +     D11X,D11Y,D11Z,SOURL)
      GOTO 100

 6    WRITE(*,*) 'Exact diffusion by Hanjalic and Launder is not yet'//
     +           'implemented (26.4.1995)'
      WRITE(*,*) 'Exiting...'
      STOP

 100  CONTINUE

      RETURN
      END

C**************************************************
C*    DIFFUSION BASED ON THINLAYER APPROXIMATION  *
C**************************************************

      SUBROUTINE DIFSCA(RO,RK,REPS,VIS,EPS2,VOL,A,AX,AY,AZ,D,
     +     RUU,RUV,RUW,RVV,RVW,RWW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     DXX,DXY,DXZ,DYY,DYZ,DZZ,VXX,VXY,VXZ,VYY,VYZ,VZZ,
     +     IDI,PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,SOURL,
     +     KSTR,IDIR)

C ... THIS SUBROUTINE DETERMINES DIFFUSION TERM IN SECOND-MOMENT CLOS.
C ... THIS IS SO CALLED SCALAR DIFFUSION MODEL
      
      LOGICAL SOURL,INVIS
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A(*),AX(*),AY(*),AZ(*),D(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     FUU(*),FUV(*),FUW(*),FVV(*),FVW(*),FWW(*),
     +     DXX(*),DXY(*),DXZ(*),DYY(*),DYZ(*),DZZ(*),
     +     VXX(*),VXY(*),VXZ(*),VYY(*),VYZ(*),VZZ(*),
     +     PSIGSC(*),PSIGS2(*)

      EPS     = 1.E-10

C ... SCHMIDTS NUMBERS FOR REYNOLDS STRESSES

      PUU1  = 1./PSIGSC(1)
      PUU2  = 1./PSIGS2(1)
      PUV1  = 1./PSIGSC(2)
      PUV2  = 1./PSIGS2(2)
      PUW1  = 1./PSIGSC(3)
      PUW2  = 1./PSIGS2(3)
      PVV1  = 1./PSIGSC(4)
      PVV2  = 1./PSIGS2(4)
      PVW1  = 1./PSIGSC(5)
      PVW2  = 1./PSIGS2(5)
      PWW1  = 1./PSIGSC(6)
      PWW2  = 1./PSIGS2(6)
C
C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (2ND MOMENT CL. TURB.)
C

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ILL     = ISTRID*JSTRID

      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

C ... MOLECYLAR VISCOSITY (FOR PRINTING)
C     ==================

      IF(IDI /= 0) THEN
      IF (SOURL) THEN
      DO 900 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 850 IG = 1,IMAXID
         I      = IA2 + IG
         YPD2    =  A(I)/(VOL(I) + VOL(I-IL))
         SURG    = -A(I)*YPD2
         SURF    = SURG*(VIS(I) + VIS(I-IL))
         YPRO1   = 1./RO(I)
         YPRO2   = 1./RO(I-IL)
         
         FUU(I)  = (PUU2*SURF)*(RUU(I)*YPRO1 - RUU(I-IL)*YPRO2)
         FUV(I)  = (PUV2*SURF)*(RUV(I)*YPRO1 - RUV(I-IL)*YPRO2)
         FUW(I)  = (PUW2*SURF)*(RUW(I)*YPRO1 - RUW(I-IL)*YPRO2)
         FVV(I)  = (PVV2*SURF)*(RVV(I)*YPRO1 - RVV(I-IL)*YPRO2)
         FVW(I)  = (PVW2*SURF)*(RVW(I)*YPRO1 - RVW(I-IL)*YPRO2)
         FWW(I)  = (PWW2*SURF)*(RWW(I)*YPRO1 - RWW(I-IL)*YPRO2)
 850  CONTINUE
 900  CONTINUE

      DO 1200 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 1100 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 1100 IG = 1,IMAXID
         N      = IA2 + IG
         DT1    = 1./VOL(N)
         VXX(N) = VXX(N)-(FUU(N+IL) - FUU(N))*DT1
         VXY(N) = VXY(N)-(FUV(N+IL) - FUV(N))*DT1
         VXZ(N) = VXZ(N)-(FUW(N+IL) - FUW(N))*DT1
         VYY(N) = VYY(N)-(FVV(N+IL) - FVV(N))*DT1
         VYZ(N) = VYZ(N)-(FVW(N+IL) - FVW(N))*DT1
         VZZ(N) = VZZ(N)-(FWW(N+IL) - FWW(N))*DT1
 1100 CONTINUE
 1200 CONTINUE

C ... TURBULENT DIFFUSION (FOR PRINTING)
C     ===================
      DO 2000 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 1850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 1850 IG = 1,IMAXID
         I      = IA2 + IG

C ... SCHMIDT'S NUMBER  = PSIUU = 1/PSIUU1
C ... SURFU1 IS TURBULENT AND SURFU2 IS MOLECYLAR VISCOSITY

         YPD2   = A(I)/(VOL(I) + VOL(I-IL))
         SURG   = -A(I)*YPD2
         SURF   = SURG*((EPS2(I)-1.)*VIS(I) + (EPS2(I-IL)-1.)*VIS(I-IL))
         
         YPRO1  = 1./RO(I)
         YPRO2  = 1./RO(I-IL)
         
         FUU(I) = (PUU1*SURF)*(RUU(I)*YPRO1 - RUU(I-IL)*YPRO2)
         FUV(I) = (PUV1*SURF)*(RUV(I)*YPRO1 - RUV(I-IL)*YPRO2)
         FUW(I) = (PUW1*SURF)*(RUW(I)*YPRO1 - RUW(I-IL)*YPRO2)
         FVV(I) = (PVV1*SURF)*(RVV(I)*YPRO1 - RVV(I-IL)*YPRO2)
         FVW(I) = (PVW1*SURF)*(RVW(I)*YPRO1 - RVW(I-IL)*YPRO2)
         FWW(I) = (PWW1*SURF)*(RWW(I)*YPRO1 - RWW(I-IL)*YPRO2)
 1850 CONTINUE
 2000 CONTINUE
      
      DO 2200 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 2100 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 2100 IG = 1,IMAXID
         N      = IA2 + IG
         DT1    = 1./VOL(N)
         DXX(N) = DXX(N)-(FUU(N+IL) - FUU(N))*DT1
         DXY(N) = DXY(N)-(FUV(N+IL) - FUV(N))*DT1
         DXZ(N) = DXZ(N)-(FUW(N+IL) - FUW(N))*DT1
         DYY(N) = DYY(N)-(FVV(N+IL) - FVV(N))*DT1
         DYZ(N) = DYZ(N)-(FVW(N+IL) - FVW(N))*DT1
         DZZ(N) = DZZ(N)-(FWW(N+IL) - FWW(N))*DT1
 2100 CONTINUE
 2200 CONTINUE
      ENDIF

C ... TURBULENT AND MOLECYLAR DIFFUSION
C     =================================
      DO 4000 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 4850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 4850 IG = 1,IMAXID
         I      = IA2 + IG
C ... SCHMIDT'S NUMBER  = PSIUU = 1/PSIUU1
C ... SURFU1 IS TURBULENT AND SURFU2 IS MOLECYLAR VISCOSITY

         YPD2   = A(I)/(VOL(I) + VOL(I-IL))
         SURG   = -A(I)*YPD2
         SURF   = SURG*((EPS2(I)-1.)*VIS(I) + (EPS2(I-IL)-1.)*VIS(I-IL))
         SURV   = SURG*(VIS(I) + VIS(I-IL))
         
         YPRO1  = 1./RO(I)
         YPRO2  = 1./RO(I-IL)
         
         FUU(I) = (PUU1*SURF+PUU2*SURV)*(RUU(I)*YPRO1-RUU(I-IL)*YPRO2)
         FUV(I) = (PUV1*SURF+PUV2*SURV)*(RUV(I)*YPRO1-RUV(I-IL)*YPRO2)
         FUW(I) = (PUW1*SURF+PUW2*SURV)*(RUW(I)*YPRO1-RUW(I-IL)*YPRO2)
         FVV(I) = (PVV1*SURF+PVV2*SURV)*(RVV(I)*YPRO1-RVV(I-IL)*YPRO2)
         FVW(I) = (PVW1*SURF+PVW2*SURV)*(RVW(I)*YPRO1-RVW(I-IL)*YPRO2)
         FWW(I) = (PWW1*SURF+PWW2*SURV)*(RWW(I)*YPRO1-RWW(I-IL)*YPRO2)
 4850 CONTINUE
 4000 CONTINUE
      ENDIF

C ... LAMINAR FLOW
      IF(INVIS) THEN
      DO 3000 KG = KMINID,KMAX
      IA      = (KN+KG-1)*ILL
      DO 2950 JG = JMINID,JMAX
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 2950 IG = IMINID,IMAX
         I       = IA2 + IG
         FUU(I)  = 0.
         FUV(I)  = 0.
         FUW(I)  = 0.
         FVV(I)  = 0.
         FVW(I)  = 0.
         FWW(I)  = 0.
 2950 CONTINUE
      IF(SOURL) THEN
      DO 2960 JG = JMINID,JMAX
      IA2        = IA + (JN+JG-1)*ISTRID + IN
      DO 2960 IG = IMINID,IMAX
         I       = IA2 + IG
         DXX(I)  = 0.
         DXY(I)  = 0.
         DXZ(I)  = 0.
         DYY(I)  = 0.
         DYZ(I)  = 0.
         DZZ(I)  = 0.
         VXX(I)  = 0.
         VXY(I)  = 0.
         VXZ(I)  = 0.
         VYY(I)  = 0.
         VYZ(I)  = 0.
         VZZ(I)  = 0.
 2960 CONTINUE
      ENDIF
 3000 CONTINUE
      ENDIF

      RETURN
      END

      SUBROUTINE DIFDLY(RO,RK,REPS,VIS,EPS2,VOL,A,AX,AY,AZ,D,
     +     RUU,RUV,RUW,RVV,RVW,RWW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     DXX,DXY,DXZ,DYY,DYZ,DZZ,VXX,VXY,VXZ,VYY,VYZ,VZZ,
     +     IDI,PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,SOURL,
     +     KSTR,IDIR)

C ... THIS SUBROUTINE DETERMINES DIFFUSION TERM IN SECOND-MOMENT CLOS.
C ... THIS IS Daly and Harlow model (1970)

      LOGICAL SOURL,INVIS
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A(*),AX(*),AY(*),AZ(*),D(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     FUU(*),FUV(*),FUW(*),FVV(*),FVW(*),FWW(*),
     +     DXX(*),DXY(*),DXZ(*),DYY(*),DYZ(*),DZZ(*),
     +     VXX(*),VXY(*),VXZ(*),VYY(*),VYZ(*),VZZ(*),
     +     PSIGSC(*),PSIGS2(*)

      EPS     = 1.E-10

      CS      = .22 
C ... values between .2 and .25 have been found in literature
C
C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (2ND MOMENT CL. TURB.)
C
C ... SCHMIDTS NUMBERS FOR REYNOLDS STRESSES

      PUU2  = 1./PSIGS2(1)
      PUV2  = 1./PSIGS2(2)
      PUW2  = 1./PSIGS2(3)
      PVV2  = 1./PSIGS2(4)
      PVW2  = 1./PSIGS2(5)
      PWW2  = 1./PSIGS2(6)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ILL     = ISTRID*JSTRID

      IL      = KSTR

      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

C ... MOLECYLAR VISCOSITY (FOR PRINTING)
C     ==================
      IF(IDI /= 0) THEN
      IF (SOURL) THEN
      DO 900 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 850 IG = 1,IMAXID
         I      = IA2 + IG
         YPD2    =  A(I)/(VOL(I) + VOL(I-IL))
         SURG    = -A(I)*YPD2
         SURF    = SURG*(VIS(I) + VIS(I-IL))
         YPRO1   = 1./RO(I)
         YPRO2   = 1./RO(I-IL)
         
         FUU(I)  = (PUU2*SURF)*(RUU(I)*YPRO1 - RUU(I-IL)*YPRO2)
         FUV(I)  = (PUV2*SURF)*(RUV(I)*YPRO1 - RUV(I-IL)*YPRO2)
         FUW(I)  = (PUW2*SURF)*(RUW(I)*YPRO1 - RUW(I-IL)*YPRO2)
         FVV(I)  = (PVV2*SURF)*(RVV(I)*YPRO1 - RVV(I-IL)*YPRO2)
         FVW(I)  = (PVW2*SURF)*(RVW(I)*YPRO1 - RVW(I-IL)*YPRO2)
         FWW(I)  = (PWW2*SURF)*(RWW(I)*YPRO1 - RWW(I-IL)*YPRO2)
 850  CONTINUE
 900  CONTINUE

      DO 1200 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 1100 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 1100 IG = 1,IMAXID
         N      = IA2 + IG
         DT1    = 1./VOL(N)
         VXX(N) = VXX(N)-(FUU(N+IL) - FUU(N))*DT1
         VXY(N) = VXY(N)-(FUV(N+IL) - FUV(N))*DT1
         VXZ(N) = VXZ(N)-(FUW(N+IL) - FUW(N))*DT1
         VYY(N) = VYY(N)-(FVV(N+IL) - FVV(N))*DT1
         VYZ(N) = VYZ(N)-(FVW(N+IL) - FVW(N))*DT1
         VZZ(N) = VZZ(N)-(FWW(N+IL) - FWW(N))*DT1
 1100 CONTINUE
 1200 CONTINUE

C ... TURBULENT DIFFUSION (FOR PRINTING)
C     ===================
      DO 2000 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 1850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 1850 IG = 1,IMAXID
         I      = IA2 + IG

         YPD2   = A(I)/(VOL(I) + VOL(I-IL))
         SURG   = -A(I)*YPD2
         CORR   = CS*(RK(I)+ RK(I-IL))/(REPS(I) + REPS(I-IL))
         DDIF   = (
     +    AX(I)**2*(RUU(I)+RUU(I-IL))+2.*AX(I)*AY(I)*(RUV(I)+RUV(I-IL))+
     +    AY(I)**2*(RVV(I)+RVV(I-IL))+2.*AY(I)*AZ(I)*(RVW(I)+RVW(I-IL))+
     +    AZ(I)**2*(RWW(I)+RWW(I-IL))+2.*AX(I)*AZ(I)*(RUW(I)+RUW(I-IL)))
         SURF   = SURG*CORR*DDIF

         YPRO1  = 1./RO(I)
         YPRO2  = 1./RO(I-IL)
         
         FUU(I) = SURF*(RUU(I)*YPRO1 - RUU(I-IL)*YPRO2)
         FUV(I) = SURF*(RUV(I)*YPRO1 - RUV(I-IL)*YPRO2)
         FUW(I) = SURF*(RUW(I)*YPRO1 - RUW(I-IL)*YPRO2)
         FVV(I) = SURF*(RVV(I)*YPRO1 - RVV(I-IL)*YPRO2)
         FVW(I) = SURF*(RVW(I)*YPRO1 - RVW(I-IL)*YPRO2)
         FWW(I) = SURF*(RWW(I)*YPRO1 - RWW(I-IL)*YPRO2)
 1850 CONTINUE
 2000 CONTINUE
      
      DO 2200 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 2100 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 2100 IG = 1,IMAXID
         N      = IA2 + IG
         DT1    = 1./VOL(N)
         DXX(N) = DXX(N)-(FUU(N+IL) - FUU(N))*DT1
         DXY(N) = DXY(N)-(FUV(N+IL) - FUV(N))*DT1
         DXZ(N) = DXZ(N)-(FUW(N+IL) - FUW(N))*DT1
         DYY(N) = DYY(N)-(FVV(N+IL) - FVV(N))*DT1
         DYZ(N) = DYZ(N)-(FVW(N+IL) - FVW(N))*DT1
         DZZ(N) = DZZ(N)-(FWW(N+IL) - FWW(N))*DT1
 2100 CONTINUE
 2200 CONTINUE
      ENDIF


C ... TURBULENT AND MOLECYLAR DIFFUSION
C     =================================
      DO 4000 KG = 1,KMAXID
      IA      = (KN+KG-1)*ILL
      DO 4850 JG = 1,JMAXID
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 4850 IG = 1,IMAXID
         I      = IA2 + IG

C ... SCHMIDT'S NUMBER  = PSIUU = 1/PSIUU1

         YPD2   = A(I)/(VOL(I) + VOL(I-IL))
         SURG   = -A(I)*YPD2
         CORR   = CS*(RK(I) + RK(I-IL))/(REPS(I) + REPS(I-IL))
         DDIF   = (
     +    AX(I)**2*(RUU(I)+RUU(I-IL))+2.*AX(I)*AY(I)*(RUV(I)+RUV(I-IL))+
     +    AY(I)**2*(RVV(I)+RVV(I-IL))+2.*AY(I)*AZ(I)*(RVW(I)+RVW(I-IL))+
     +    AZ(I)**2*(RWW(I)+RWW(I-IL))+2.*AX(I)*AZ(I)*(RUW(I)+RUW(I-IL)))
         SURF   = SURG*CORR*DDIF
         SURV   = SURG*(VIS(I) + VIS(I-IL))
         
         YPRO1  = 1./RO(I)
         YPRO2  = 1./RO(I-IL)
         
         FUU(I) = (SURF+PUU2*SURV)*(RUU(I)*YPRO1-RUU(I-IL)*YPRO2)
         FUV(I) = (SURF+PUV2*SURV)*(RUV(I)*YPRO1-RUV(I-IL)*YPRO2)
         FUW(I) = (SURF+PUW2*SURV)*(RUW(I)*YPRO1-RUW(I-IL)*YPRO2)
         FVV(I) = (SURF+PVV2*SURV)*(RVV(I)*YPRO1-RVV(I-IL)*YPRO2)
         FVW(I) = (SURF+PVW2*SURV)*(RVW(I)*YPRO1-RVW(I-IL)*YPRO2)
         FWW(I) = (SURF+PWW2*SURV)*(RWW(I)*YPRO1-RWW(I-IL)*YPRO2)
 4850 CONTINUE
 4000 CONTINUE
      ENDIF

C ... LAMINAR FLOW
      IF(INVIS) THEN
      DO 3000 KG = KMINID,KMAX
      IA      = (KN+KG-1)*ILL
      DO 2950 JG = JMINID,JMAX
      IA2       = IA + (JN+JG-1)*ISTRID + IN
      DO 2950 IG = IMINID,IMAX
         I       = IA2 + IG
         FUU(I)  = 0.
         FUV(I)  = 0.
         FUW(I)  = 0.
         FVV(I)  = 0.
         FVW(I)  = 0.
         FWW(I)  = 0.
 2950 CONTINUE
      IF(SOURL) THEN
      DO 2960 JG = JMINID,JMAX
      IA2        = IA + (JN+JG-1)*ISTRID + IN
      DO 2960 IG = IMINID,IMAX
         I       = IA2 + IG
         DXX(I)  = 0.
         DXY(I)  = 0.
         DXZ(I)  = 0.
         DYY(I)  = 0.
         DYZ(I)  = 0.
         DZZ(I)  = 0.
         VXX(I)  = 0.
         VXY(I)  = 0.
         VXZ(I)  = 0.
         VYY(I)  = 0.
         VYZ(I)  = 0.
         VZZ(I)  = 0.
 2960 CONTINUE
      ENDIF
 3000 CONTINUE
      ENDIF

      RETURN
      END

C***********************************************
C*    EXACT DIFFUSION SUBROUTINES              *
C***********************************************

      SUBROUTINE DIFEXS(RO,RK,REPS,VIS,EPS2,VOL,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     RUU,RUV,RUW,RVV,RVW,RWW,DRUU,DUU,
     +     TUU,VUU,
     +     PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,
     +     D11X,D11Y,D11Z,SOURL)

C ... THIS SUBROUTINE DETERMINES DIFFUSION TERM IN SECOND-MOMENT CLOS.
C ... THIS IS Scalar diffusion model. Thin layer approximation is
C ... not applied.

      LOGICAL SOURL
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +     A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     DUU(*),DRUU(*),TUU(*),VUU(*),
     +     D11X(*),D11Y(*),D11Z(*)

      EPS     = 1.E-10

C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (2ND MOMENT CL. TURB.)
C
C ... SCHMIDTS NUMBERS FOR REYNOLDS STRESSES

      PUU1  = 1./PSIGSC
      PUU2  = 1./PSIGS2

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      NTOT    = ISTRID*JSTRID*KSTRID

      KMAXP1  = KMAX + 1
      IL      = ISTRID*JSTRID

      DO 5000 KG = 1,KMAX
      IA  = (KN+KG-1)*IL + JN*ISTRID
      DO 1000 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         YPRO1   = 1./RO(L)
         YPRO2   = 1./RO(L-1)
         YPRO3   = 1./RO(L+1)
         U11MX   = YPRO1*DRUU(L) + YPRO2*DRUU(L-1)
         U11PX   = YPRO1*DRUU(L) + YPRO3*DRUU(L+1)
C ... ETA-DIRECTION
         YPRO2   = 1./RO(L-ISTRID)
         YPRO3   = 1./RO(L+ISTRID)
         U11MY   = YPRO1*DRUU(L) + YPRO2*DRUU(L-ISTRID)
         U11PY   = YPRO1*DRUU(L) + YPRO3*DRUU(L+ISTRID)
C ... ZETA DIRECTION
         YPRO2   = 1./RO(L-IL)
         YPRO3   = 1./RO(L+IL)
         U11MZ   = YPRO1*DRUU(L) + YPRO2*DRUU(L-IL)
         U11PZ   = YPRO1*DRUU(L) + YPRO3*DRUU(L+IL)
                 
         PVOL    = .5/VOL(L)
         D11XI   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    U11PZ - A3(L)*A3X(L)*U11MZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*U11PY - A2(L)*A2X(L)*U11MY +
     3        A1(L+1)     *A1X(L+1)*     U11PX - A1(L)*A1X(L)*U11MX)
         D11YI   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    U11PZ - A3(L)*A3Y(L)*U11MZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*U11PY - A2(L)*A2Y(L)*U11MY +
     3        A1(L+1)     *A1Y(L+1)*     U11PX - A1(L)*A1Y(L)*U11MX)
         D11ZI   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    U11PZ - A3(L)*A3Z(L)*U11MZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*U11PY - A2(L)*A2Z(L)*U11MY +
     3        A1(L+1)     *A1Z(L+1)*     U11PX - A1(L)*A1Z(L)*U11MX)

         D11X(L) = D11XI
         D11Y(L) = D11YI
         D11Z(L) = D11ZI
 1000 CONTINUE
 5000 CONTINUE

      CALL REFDER(D11X,D11Y,D11Z,IMAX,JMAX,KMAX,IN,JN,KN)

C ... MOLECYLAR VISCOSITY
C     ===================
      DO 1200 KG = 1,KMAX
      IA        = (KN+KG-1)*IL + JN*ISTRID
      DO 1100 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         ICENT   = L
         ILEFT   = L - 1
         IRIGH   = L + 1
         D11MXX  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYX  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZX  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXX  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYX  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZX  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

C ... ETA-DIRECTION
         ICENT   = L
         ILEFT   = L - ISTRID
         IRIGH   = L + ISTRID
         D11MXY  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYY  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZY  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXY  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYY  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZY  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

C ... ZETA DIRECTION
         ICENT   = L
         ILEFT   = L - IL
         IRIGH   = L + IL
         D11MXZ  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYZ  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZZ  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXZ  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYZ  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZZ  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

         PVOL    = .5/VOL(L)
         D11X2   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    D11PXZ - A3(L)*A3X(L)*D11MXZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*D11PXY - A2(L)*A2X(L)*D11MXY +
     3        A1(L+1)     *A1X(L+1)*     D11PXX - A1(L)*A1X(L)*D11MXX)
         D11Y2   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    D11PYZ - A3(L)*A3Y(L)*D11MYZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*D11PYY - A2(L)*A2Y(L)*D11MYY +
     3        A1(L+1)     *A1Y(L+1)*     D11PYX - A1(L)*A1Y(L)*D11MYX)
         D11Z2   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    D11PZZ - A3(L)*A3Z(L)*D11MZZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*D11PZY - A2(L)*A2Z(L)*D11MZY +
     3        A1(L+1)     *A1Z(L+1)*     D11PZX - A1(L)*A1Z(L)*D11MZX)

         DUU(L)  = PUU2*(D11X2 + D11Y2 + D11Z2)
 1100 CONTINUE
 1200 CONTINUE

C ... MOLECYLAR VISCOSITY (FOR PRINTING)
C     ===================
      IF (SOURL) THEN
      DO 1300 N = 1,NTOT
         VUU(N) = DUU(N)
 1300 CONTINUE
      ENDIF

C ... TURBULENT DIFFUSION
C     ===================
      DO 1500 KG = 1,KMAX
      IA        = (KN+KG-1)*IL + JN*ISTRID
      DO 1400 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         ICENT   = L
         ILEFT   = L - 1
         IRIGH   = L + 1
         VICENT  = (EPS2(ICENT)-1.)*VIS(ICENT)
         VILEFT  = (EPS2(ILEFT)-1.)*VIS(ILEFT)
         VIRIGH  = (EPS2(IRIGH)-1.)*VIS(IRIGH)
         D11MXX  = VICENT*D11X(ICENT) + VILEFT*D11X(ILEFT)
         D11MYX  = VICENT*D11Y(ICENT) + VILEFT*D11Y(ILEFT)
         D11MZX  = VICENT*D11Z(ICENT) + VILEFT*D11Z(ILEFT)

         D11PXX  = VICENT*D11X(ICENT) + VIRIGH*D11X(IRIGH)
         D11PYX  = VICENT*D11Y(ICENT) + VIRIGH*D11Y(IRIGH)
         D11PZX  = VICENT*D11Z(ICENT) + VIRIGH*D11Z(IRIGH)

C ... ETA-DIRECTION
         ILEFT   = L - ISTRID
         IRIGH   = L + ISTRID
         VILEFT  = (EPS2(ILEFT)-1.)*VIS(ILEFT)
         VIRIGH  = (EPS2(IRIGH)-1.)*VIS(IRIGH)
         D11MXY  = VICENT*D11X(ICENT) + VILEFT*D11X(ILEFT)
         D11MYY  = VICENT*D11Y(ICENT) + VILEFT*D11Y(ILEFT)
         D11MZY  = VICENT*D11Z(ICENT) + VILEFT*D11Z(ILEFT)

         D11PXY  = VICENT*D11X(ICENT) + VIRIGH*D11X(IRIGH)
         D11PYY  = VICENT*D11Y(ICENT) + VIRIGH*D11Y(IRIGH)
         D11PZY  = VICENT*D11Z(ICENT) + VIRIGH*D11Z(IRIGH)

C ... ZETA DIRECTION
         ILEFT   = L - IL
         IRIGH   = L + IL
         VILEFT  = (EPS2(ILEFT)-1.)*VIS(ILEFT)
         VIRIGH  = (EPS2(IRIGH)-1.)*VIS(IRIGH)
         D11MXZ  = VICENT*D11X(ICENT) + VILEFT*D11X(ILEFT)
         D11MYZ  = VICENT*D11Y(ICENT) + VILEFT*D11Y(ILEFT)
         D11MZZ  = VICENT*D11Z(ICENT) + VILEFT*D11Z(ILEFT)

         D11PXZ  = VICENT*D11X(ICENT) + VIRIGH*D11X(IRIGH)
         D11PYZ  = VICENT*D11Y(ICENT) + VIRIGH*D11Y(IRIGH)
         D11PZZ  = VICENT*D11Z(ICENT) + VIRIGH*D11Z(IRIGH)

         PVOL    = .5/VOL(L)
         D11X2   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    D11PXZ - A3(L)*A3X(L)*D11MXZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*D11PXY - A2(L)*A2X(L)*D11MXY +
     3        A1(L+1)     *A1X(L+1)*     D11PXX - A1(L)*A1X(L)*D11MXX)
         D11Y2   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    D11PYZ - A3(L)*A3Y(L)*D11MYZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*D11PYY - A2(L)*A2Y(L)*D11MYY +
     3        A1(L+1)     *A1Y(L+1)*     D11PYX - A1(L)*A1Y(L)*D11MYX)
         D11Z2   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    D11PZZ - A3(L)*A3Z(L)*D11MZZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*D11PZY - A2(L)*A2Z(L)*D11MZY +
     3        A1(L+1)     *A1Z(L+1)*     D11PZX - A1(L)*A1Z(L)*D11MZX)
         DUU(L)  = DUU(L) + PUU1*(D11X2 + D11Y2 + D11Z2)

 1400 CONTINUE
 1500 CONTINUE

C ... TURBULENT DIFFUSION (FOR PRINTING)
C     ===================

      IF (SOURL) THEN
      DO 1600 N = 1,NTOT
         TUU(N) = (DUU(N) - VUU(N))
 1600 CONTINUE
      ENDIF

      
      RETURN
      END


      SUBROUTINE DIFEXH(RO,RK,REPS,VIS,EPS2,VOL,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     RUU,RUV,RUW,RVV,RVW,RWW,DRUU,DUU,
     +     TUU,VUU,
     +     PSIGSC,PSIGS2,JRDIF,IMAX,JMAX,KMAX,IN,JN,KN,
     +     D11X,D11Y,D11Z,SOURL)

C ... THIS SUBROUTINE DETERMINES DIFFUSION TERM IN SECOND-MOMENT CLOS.
C ... THIS IS Daly and Harlow model (1970). Thin layer approximation is
C ... not applied.

      USE CHARACTERS
cc      INCLUDE 'NS3CH.C'

      LOGICAL SOURL
      DIMENSION  RO(*),RK(*),REPS(*),VIS(*),EPS2(*),
     +     VOL(*),A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +     A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +     RUU(*),RUV(*),RUW(*),RVV(*),RVW(*),RWW(*),
     +     DUU(*),DRUU(*),TUU(*),VUU(*),
     +     D11X(*),D11Y(*),D11Z(*)

      

      EPS     = 1.E-10

      CS      = .22 
C ... values between .2 and .25 have been found in literature
C
C ... VISCOUS FLUXES FOR THREE-DIMENSIONAL FLOW (2ND MOMENT CL. TURB.)
C
C ... SCHMIDTS NUMBERS FOR REYNOLDS STRESSES

      PUU2  = 1./PSIGS2

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      NTOT    = ISTRID*JSTRID*KSTRID

      KMAXP1  = KMAX + 1
      IL      = ISTRID*JSTRID

      DO 5000 KG = 1,KMAX
      IA  = (KN+KG-1)*IL + JN*ISTRID
      DO 1000 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         YPRO1   = 1./RO(L)
         YPRO2   = 1./RO(L-1)
         YPRO3   = 1./RO(L+1)
         U11MX   = YPRO1*DRUU(L) + YPRO2*DRUU(L-1)
         U11PX   = YPRO1*DRUU(L) + YPRO3*DRUU(L+1)
C ... ETA-DIRECTION
         YPRO2   = 1./RO(L-ISTRID)
         YPRO3   = 1./RO(L+ISTRID)
         U11MY   = YPRO1*DRUU(L) + YPRO2*DRUU(L-ISTRID)
         U11PY   = YPRO1*DRUU(L) + YPRO3*DRUU(L+ISTRID)
C ... ZETA DIRECTION
         YPRO2   = 1./RO(L-IL)
         YPRO3   = 1./RO(L+IL)
         U11MZ   = YPRO1*DRUU(L) + YPRO2*DRUU(L-IL)
         U11PZ   = YPRO1*DRUU(L) + YPRO3*DRUU(L+IL)
                 
         PVOL    = .5/VOL(L)
         D11XI   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    U11PZ - A3(L)*A3X(L)*U11MZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*U11PY - A2(L)*A2X(L)*U11MY +
     3        A1(L+1)     *A1X(L+1)*     U11PX - A1(L)*A1X(L)*U11MX)
         D11YI   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    U11PZ - A3(L)*A3Y(L)*U11MZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*U11PY - A2(L)*A2Y(L)*U11MY +
     3        A1(L+1)     *A1Y(L+1)*     U11PX - A1(L)*A1Y(L)*U11MX)
         D11ZI   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    U11PZ - A3(L)*A3Z(L)*U11MZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*U11PY - A2(L)*A2Z(L)*U11MY +
     3        A1(L+1)     *A1Z(L+1)*     U11PX - A1(L)*A1Z(L)*U11MX)

         D11X(L) = D11XI
         D11Y(L) = D11YI
         D11Z(L) = D11ZI
 1000 CONTINUE
 5000 CONTINUE

      CALL REFDER(D11X,D11Y,D11Z,IMAX,JMAX,KMAX,IN,JN,KN)

C ... MOLECYLAR VISCOSITY
C     ===================
      DO 1200 KG = 1,KMAX
      IA        = (KN+KG-1)*IL + JN*ISTRID
      DO 1100 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         ICENT   = L
         ILEFT   = L - 1
         IRIGH   = L + 1
         D11MXX  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYX  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZX  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXX  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYX  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZX  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

C ... ETA-DIRECTION
         ICENT   = L
         ILEFT   = L - ISTRID
         IRIGH   = L + ISTRID
         D11MXY  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYY  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZY  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXY  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYY  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZY  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

C ... ZETA DIRECTION
         ICENT   = L
         ILEFT   = L - IL
         IRIGH   = L + IL
         D11MXZ  = VIS(ICENT)*D11X(ICENT) + VIS(ILEFT)*D11X(ILEFT)
         D11MYZ  = VIS(ICENT)*D11Y(ICENT) + VIS(ILEFT)*D11Y(ILEFT)
         D11MZZ  = VIS(ICENT)*D11Z(ICENT) + VIS(ILEFT)*D11Z(ILEFT)

         D11PXZ  = VIS(ICENT)*D11X(ICENT) + VIS(IRIGH)*D11X(IRIGH)
         D11PYZ  = VIS(ICENT)*D11Y(ICENT) + VIS(IRIGH)*D11Y(IRIGH)
         D11PZZ  = VIS(ICENT)*D11Z(ICENT) + VIS(IRIGH)*D11Z(IRIGH)

         PVOL    = .5/VOL(L)
         D11X2   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    D11PXZ - A3(L)*A3X(L)*D11MXZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*D11PXY - A2(L)*A2X(L)*D11MXY +
     3        A1(L+1)     *A1X(L+1)*     D11PXX - A1(L)*A1X(L)*D11MXX)
         D11Y2   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    D11PYZ - A3(L)*A3Y(L)*D11MYZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*D11PYY - A2(L)*A2Y(L)*D11MYY +
     3        A1(L+1)     *A1Y(L+1)*     D11PYX - A1(L)*A1Y(L)*D11MYX)
         D11Z2   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    D11PZZ - A3(L)*A3Z(L)*D11MZZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*D11PZY - A2(L)*A2Z(L)*D11MZY +
     3        A1(L+1)     *A1Z(L+1)*     D11PZX - A1(L)*A1Z(L)*D11MZX)

         DUU(L)  = PUU2*(D11X2 + D11Y2 + D11Z2)
 1100 CONTINUE
 1200 CONTINUE

C ... MOLECYLAR VISCOSITY (FOR PRINTING)
C     ===================
      IF (SOURL) THEN
      DO 1300 N = 1,NTOT
         VUU(N) = DUU(N)
 1300 CONTINUE
      ENDIF

C ... TURBULENT DIFFUSION
C     ===================
      DO 1500 KG = 1,KMAX
      IA        = (KN+KG-1)*IL + JN*ISTRID
      DO 1400 IG = 1,JMAX*ISTRID
         L       = IA + IG

C ... XI-DIRECTION
         ICE     = L
         ILE     = L - 1
         IRI     = L + 1
         REC     = CS*RK(ICE)/(REPS(ICE)+EPS)
         REL     = CS*RK(ILE)/(REPS(ILE)+EPS)
         RER     = CS*RK(IRI)/(REPS(IRI)+EPS)
         D11MXX  = REC*RUU(ICE)*D11X(ICE) + REL*RUU(ILE)*D11X(ILE) +
     +             REC*RUV(ICE)*D11Y(ICE) + REL*RUV(ILE)*D11Y(ILE) +
     +             REC*RUW(ICE)*D11Z(ICE) + REL*RUW(ILE)*D11Z(ILE)
         D11MYX  = REC*RUV(ICE)*D11X(ICE) + REL*RUV(ILE)*D11X(ILE) +
     +             REC*RVV(ICE)*D11Y(ICE) + REL*RVV(ILE)*D11Y(ILE) +
     +             REC*RVW(ICE)*D11Z(ICE) + REL*RVW(ILE)*D11Z(ILE)
         D11MZX  = REC*RUW(ICE)*D11X(ICE) + REL*RUW(ILE)*D11X(ILE) +
     +             REC*RVW(ICE)*D11Y(ICE) + REL*RVW(ILE)*D11Y(ILE) +
     +             REC*RWW(ICE)*D11Z(ICE) + REL*RWW(ILE)*D11Z(ILE)

         D11PXX  = REC*RUU(ICE)*D11X(ICE) + RER*RUU(IRI)*D11X(IRI) +
     +             REC*RUV(ICE)*D11Y(ICE) + RER*RUV(IRI)*D11Y(IRI) +
     +             REC*RUW(ICE)*D11Z(ICE) + RER*RUW(IRI)*D11Z(IRI)
         D11PYX  = REC*RUV(ICE)*D11X(ICE) + RER*RUV(IRI)*D11X(IRI) +
     +             REC*RVV(ICE)*D11Y(ICE) + RER*RVV(IRI)*D11Y(IRI) +
     +             REC*RVW(ICE)*D11Z(ICE) + RER*RVW(IRI)*D11Z(IRI)
         D11PZX  = REC*RUW(ICE)*D11X(ICE) + RER*RUW(IRI)*D11X(IRI) +
     +             REC*RVW(ICE)*D11Y(ICE) + RER*RVW(IRI)*D11Y(IRI) +
     +             REC*RWW(ICE)*D11Z(ICE) + RER*RWW(IRI)*D11Z(IRI)

C ... ETA-DIRECTION
         ICE     = L
         ILE     = L - ISTRID
         IRI     = L + ISTRID
         REC     = CS*RK(ICE)/(REPS(ICE)+EPS)
         REL     = CS*RK(ILE)/(REPS(ILE)+EPS)
         RER     = CS*RK(IRI)/(REPS(IRI)+EPS)
         D11MXY  = REC*RUU(ICE)*D11X(ICE) + REL*RUU(ILE)*D11X(ILE) +
     +             REC*RUV(ICE)*D11Y(ICE) + REL*RUV(ILE)*D11Y(ILE) +
     +             REC*RUW(ICE)*D11Z(ICE) + REL*RUW(ILE)*D11Z(ILE)
         D11MYY  = REC*RUV(ICE)*D11X(ICE) + REL*RUV(ILE)*D11X(ILE) +
     +             REC*RVV(ICE)*D11Y(ICE) + REL*RVV(ILE)*D11Y(ILE) +
     +             REC*RVW(ICE)*D11Z(ICE) + REL*RVW(ILE)*D11Z(ILE)
         D11MZY  = REC*RUW(ICE)*D11X(ICE) + REL*RUW(ILE)*D11X(ILE) +
     +             REC*RVW(ICE)*D11Y(ICE) + REL*RVW(ILE)*D11Y(ILE) +
     +             REC*RWW(ICE)*D11Z(ICE) + REL*RWW(ILE)*D11Z(ILE)

         D11PXY  = REC*RUU(ICE)*D11X(ICE) + RER*RUU(IRI)*D11X(IRI) +
     +             REC*RUV(ICE)*D11Y(ICE) + RER*RUV(IRI)*D11Y(IRI) +
     +             REC*RUW(ICE)*D11Z(ICE) + RER*RUW(IRI)*D11Z(IRI)
         D11PYY  = REC*RUV(ICE)*D11X(ICE) + RER*RUV(IRI)*D11X(IRI) +
     +             REC*RVV(ICE)*D11Y(ICE) + RER*RVV(IRI)*D11Y(IRI) +
     +             REC*RVW(ICE)*D11Z(ICE) + RER*RVW(IRI)*D11Z(IRI)
         D11PZY  = REC*RUW(ICE)*D11X(ICE) + RER*RUW(IRI)*D11X(IRI) +
     +             REC*RVW(ICE)*D11Y(ICE) + RER*RVW(IRI)*D11Y(IRI) +
     +             REC*RWW(ICE)*D11Z(ICE) + RER*RWW(IRI)*D11Z(IRI)

C ... ZETA-DIRECTION
         ICE     = L
         ILE     = L - IL
         IRI     = L + IL
         REC     = CS*RK(ICE)/(REPS(ICE)+EPS)
         REL     = CS*RK(ILE)/(REPS(ILE)+EPS)
         RER     = CS*RK(IRI)/(REPS(IRI)+EPS)
         D11MXZ  = REC*RUU(ICE)*D11X(ICE) + REL*RUU(ILE)*D11X(ILE) +
     +             REC*RUV(ICE)*D11Y(ICE) + REL*RUV(ILE)*D11Y(ILE) +
     +             REC*RUW(ICE)*D11Z(ICE) + REL*RUW(ILE)*D11Z(ILE)
         D11MYZ  = REC*RUV(ICE)*D11X(ICE) + REL*RUV(ILE)*D11X(ILE) +
     +             REC*RVV(ICE)*D11Y(ICE) + REL*RVV(ILE)*D11Y(ILE) +
     +             REC*RVW(ICE)*D11Z(ICE) + REL*RVW(ILE)*D11Z(ILE)
         D11MZZ  = REC*RUW(ICE)*D11X(ICE) + REL*RUW(ILE)*D11X(ILE) +
     +             REC*RVW(ICE)*D11Y(ICE) + REL*RVW(ILE)*D11Y(ILE) +
     +             REC*RWW(ICE)*D11Z(ICE) + REL*RWW(ILE)*D11Z(ILE)

         D11PXZ  = REC*RUU(ICE)*D11X(ICE) + RER*RUU(IRI)*D11X(IRI) +
     +             REC*RUV(ICE)*D11Y(ICE) + RER*RUV(IRI)*D11Y(IRI) +
     +             REC*RUW(ICE)*D11Z(ICE) + RER*RUW(IRI)*D11Z(IRI)
         D11PYZ  = REC*RUV(ICE)*D11X(ICE) + RER*RUV(IRI)*D11X(IRI) +
     +             REC*RVV(ICE)*D11Y(ICE) + RER*RVV(IRI)*D11Y(IRI) +
     +             REC*RVW(ICE)*D11Z(ICE) + RER*RVW(IRI)*D11Z(IRI)
         D11PZZ  = REC*RUW(ICE)*D11X(ICE) + RER*RUW(IRI)*D11X(IRI) +
     +             REC*RVW(ICE)*D11Y(ICE) + RER*RVW(IRI)*D11Y(IRI) +
     +             REC*RWW(ICE)*D11Z(ICE) + RER*RWW(IRI)*D11Z(IRI)

         PVOL    = .5/VOL(L)
         D11X2   = PVOL*(
     1        A3(L+IL)    *A3X(L+IL)*    D11PXZ - A3(L)*A3X(L)*D11MXZ +
     2        A2(L+ISTRID)*A2X(L+ISTRID)*D11PXY - A2(L)*A2X(L)*D11MXY +
     3        A1(L+1)     *A1X(L+1)*     D11PXX - A1(L)*A1X(L)*D11MXX)
         D11Y2   = PVOL*(
     1        A3(L+IL)    *A3Y(L+IL)*    D11PYZ - A3(L)*A3Y(L)*D11MYZ +
     2        A2(L+ISTRID)*A2Y(L+ISTRID)*D11PYY - A2(L)*A2Y(L)*D11MYY +
     3        A1(L+1)     *A1Y(L+1)*     D11PYX - A1(L)*A1Y(L)*D11MYX)
         D11Z2   = PVOL*(
     1        A3(L+IL)    *A3Z(L+IL)*    D11PZZ - A3(L)*A3Z(L)*D11MZZ +
     2        A2(L+ISTRID)*A2Z(L+ISTRID)*D11PZY - A2(L)*A2Z(L)*D11MZY +
     3        A1(L+1)     *A1Z(L+1)*     D11PZX - A1(L)*A1Z(L)*D11MZX)
         DUU(L)  = DUU(L) + (D11X2 + D11Y2 + D11Z2)

 1400 CONTINUE
 1500 CONTINUE

C ... TURBULENT DIFFUSION (FOR PRINTING)
C     ===================

      IF (SOURL) THEN
      DO 1600 N = 1,NTOT
         TUU(N) = (DUU(N) - VUU(N))
 1600 CONTINUE
      ENDIF

      
      RETURN
      END

C********************************
C*    END OF DIFFUSION BLOCK    *
C********************************

      subroutine testi(rk,REPS,DDEPS,uu,sfi,pro,dis,pi,ro,vis,
     +     imax,jmax,kmax,maxsb)

      USE NS3CO, ONLY : IN, JN, KN

      real rk(*),REPS(*),DDEPS(*),ro(*),vis(*),
     +     uu(maxsb,6),sfi(maxsb,6),pro(maxsb,6),dis(maxsb,6),
     +     pi(maxsb,6)

      istrid = imax + 2*in
      jstrid = jmax + 2*jn
      IL      = ISTRID*JSTRID
      xtot = JMAX*ISTRID*kmax
      xtot2 = imax*jmax*kmax
      apu5 = 0.

      write(21,*) 'lahdetermien summat ovat:'
      do 900 ns = 1,6
      sum = 0.
      DO 9000 KG = 1,KMAX
         IA      = (KN+KG-1)*IL + JN*ISTRID
         DO 9000 IG = 1,JMAX*ISTRID
            L       = IA + IG
            apu1 = sfi(l,ns)
            apu2 = apu1-pro(l,ns)
            apu3 = apu2-dis(l,ns)
            apu4 = apu3-pi(l,ns)
            sum = sum +apu4
9000  continue
      write(21,*) ns,sum/xtot
900   continue

      xsum  = 0.
      DO 800 KG = 1,KMAX
         IA      = (KN+KG-1)*IL + JN*ISTRID
         DO 800 IG = 1,JMAX*ISTRID
            L       = IA + IG
            apu1 = .5*(uu(l,1)+uu(l,4)+uu(l,6))
            apu2 = rk(l)-apu1
            xsum = xsum +apu2
 800  continue
      write(21,*)'Avarage differece in k. energy of tur',xsum/xtot
      apu6 = 0.
      apu7 = 0.
      DO 300 KG = 1,KMAX
         IA      = (KN+KG-1)*IL + JN*ISTRID
         DO 300 JG = 1,JMAX
         DO 300 IG = 1,IMAX
            L       = IA + (JG-1)*ISTRID + (IG+IN-1) + 1
            apu5 = apu5 + .5*(dis(l,1)+dis(l,4)+dis(l,6))+
     +           reps(l)+ddeps(L)
            sumpi = pi(l,1) + pi(l,4) + pi(l,6)
            su2pi = pi(l,1)**2 + pi(l,4)**2 + pi(l,6)**2
            apu6 = apu6 + sumpi**2
            apu7 = apu7 + su2pi
c            write(88,*) ig,jg,sumpi,su2pi,sumpi/(su2pi+1.E-4)
c            uk   = sqrt(rk(l)/ro(l))
c            ueps = (vis(l)*reps(l)/ro(l)**2)**(1./4.)
c            if(ig == 30) write(88,*) jg,uk,ueps,uk/ueps
 300  continue
      WRITE(21,*) 'Avarage Difference in dissipation',apu5/XTOT2
      WRITE(21,*) 'L2 norm of total velocity-pressure gradient',
     +     SQRT(apu6/XTOT2),SQRT(apu6/XTOT2)/SQRT(apu7/3./XTOT2)
      return
      end

      SUBROUTINE ADDUUK(RK,UU,IMAX,JMAX,KMAX,MAXSB)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RK(*), UU(MAXSB,6)

C ... CALCULATE KINETIC ENERGY FROM REYNOLDS STRESSES

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      NTOT   = ISTRID*JSTRID*KSTRID

      DO 9000 L = 1,NTOT
         RK(L)  = .5*(UU(L,1)+UU(L,4)+UU(L,6))
 9000 CONTINUE

      RETURN
      END SUBROUTINE ADDUUK


      SUBROUTINE ADDKUU(RK,UU,IMAX,JMAX,KMAX,MAXSB)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: RK(*), UU(MAXSB,6)

C ... CALCULATE REYNOLDS STRESSES FROM KINETIC ENERGY

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      NTOT   = ISTRID*JSTRID*KSTRID
      EPS    = 1.E-10

      DO 9000 L = 1,NTOT
         XTOTAL  = RK(L)/(.5*(UU(L,1)+UU(L,4)+UU(L,6)+EPS))
         UU(L,1) = UU(L,1)*XTOTAL
         UU(L,4) = UU(L,4)*XTOTAL
         UU(L,6) = UU(L,6)*XTOTAL
 9000 CONTINUE

      RETURN
      END SUBROUTINE ADDKUU


      SUBROUTINE REFLUU(UU,UV,UW,VV,VW,WW,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,
     + JTOP,KBOT,KTOP,INRE,JNRE,KNRE,NPATCH,ICON,NBL)

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: UU(*), UV(*), UW(*), VV(*), VW(*), WW(*), ICON(IC9,*)

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

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
         
      IF(IBC == 7.AND.IFACE <= 3) K = 1

C ... END OF SET-UP OF INDECES. FIRSTLY PRESSURES AND VELOCITIES

      DO 1100 J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO 1100 I = I1,I2
            L       = 1 + (IN+I-1)*ISTR + KK
            LB      = L - IDIR*KSTR
            LT      = L + IDIR*KSTR
            UU(LB)  = 2.*UU(L) - UU(LT)
            UV(LB)  = 2.*UV(L) - UV(LT)
            UW(LB)  = 2.*UW(L) - UW(LT)
            VV(LB)  = 2.*VV(L) - VV(LT)
            VW(LB)  = 2.*VW(L) - VW(LT)
            WW(LB)  = 2.*WW(L) - WW(LT)
 1100 CONTINUE
      ENDIF                     ! IBC >= 7 .AND. IBC <= 10
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE REFLUU
				
				
      SUBROUTINE REFDER(UU,UV,UW,IMAX,JMAX,KMAX,IN,JN,KN)

      DIMENSION UU(*),UV(*),UW(*)

C ... THIS SUBROUTINE REFLECT DERIVATIVES OF REYNOLDS STRESSES
C ... IN THE GHOST SHELLS.

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID
C ... XI-DIRECTION
      I   = 1
      DO 1000 K = 1,KMAX
      KK      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
         II      = (JN+J-1)*ISTRID + IN + KK
         L       = I + II
         UU(L-1) = 2.*UU(L) - UU(L+1)
         UV(L-1) = 2.*UV(L) - UV(L+1)
         UW(L-1) = 2.*UW(L) - UW(L+1)
         UU(L-2) = 3.*UU(L) - 2.*UU(L+1)
         UV(L-2) = 3.*UV(L) - 2.*UV(L+1)
         UW(L-2) = 3.*UW(L) - 2.*UW(L+1)
 1000 CONTINUE

      I   = IMAX
      DO 2000 K = 1,KMAX
      KK      = (KN+K-1)*IL
      DO 2000 J = 1,JMAX
         II      = (JN+J-1)*ISTRID + IN + KK
         IC      = (KN+K-1)*JSTRID + J + JN
         L       = I + II
         UU(L+1) = 2.*UU(L) - UU(L-1)
         UV(L+1) = 2.*UV(L) - UV(L-1)
         UW(L+1) = 2.*UW(L) - UW(L-1)
         UU(L+2) = 3.*UU(L) - 2.*UU(L-1)
         UV(L+2) = 3.*UV(L) - 2.*UV(L-1)
         UW(L+2) = 3.*UW(L) - 2.*UW(L-1)
 2000 CONTINUE

C ... ETA-DIRECTION
      J   = 1
      DO 3000 K = 1,KMAX
      II      = (JN+J-1)*ISTRID + IN + KK
      DO 3000 I = 1,IMAX
         L       = I + II
         UU(L-ISTRID) = 2.*UU(L) - UU(L+ISTRID)
         UV(L-ISTRID) = 2.*UV(L) - UV(L+ISTRID)
         UW(L-ISTRID) = 2.*UW(L) - UW(L+ISTRID)
         UU(L-2*ISTRID) = 3.*UU(L) - 2.*UU(L+ISTRID)
         UV(L-2*ISTRID) = 3.*UV(L) - 2.*UV(L+ISTRID)
         UW(L-2*ISTRID) = 3.*UW(L) - 2.*UW(L+ISTRID)
3000  CONTINUE

      J   = JMAX
      DO 4000 K = 1,KMAX
      II      = (JN+J-1)*ISTRID + IN + KK
      DO 4000 I = 1,IMAX
         L       = I + II
         UU(L+ISTRID) = 2.*UU(L) - UU(L-ISTRID)
         UV(L+ISTRID) = 2.*UV(L) - UV(L-ISTRID)
         UW(L+ISTRID) = 2.*UW(L) - UW(L-ISTRID)
         UU(L+2*ISTRID) = 3.*UU(L) - 2.*UU(L-ISTRID)
         UV(L+2*ISTRID) = 3.*UV(L) - 2.*UV(L-ISTRID)
         UW(L+2*ISTRID) = 3.*UW(L) - 2.*UW(L-ISTRID)
 4000 CONTINUE

C ... ZETA DIRECTION
      K = 1
      KK      = (KN+K-1)*IL
      DO 5000 J = 1,JMAX
      II      = (JN+J-1)*ISTRID + IN + KK
      DO 5000 I = 1,IMAX
         L       = I + II
         UU(L-IL) = 2.*UU(L) - UU(L+IL)
         UV(L-IL) = 2.*UV(L) - UV(L+IL)
         UW(L-IL) = 2.*UW(L) - UW(L+IL)
         UU(L-2*IL) = 3.*UU(L) - 2.*UU(L+IL)
         UV(L-2*IL) = 3.*UV(L) - 2.*UV(L+IL)
         UW(L-2*IL) = 3.*UW(L) - 2.*UW(L+IL)
 5000 CONTINUE

C ... ZETA DIRECTION
      K = KMAX
      KK      = (KN+K-1)*IL
      DO 6000 J = 1,JMAX
      II      = (JN+J-1)*ISTRID + IN + KK
      DO 6000 I = 1,IMAX
         L       = I + II
         UU(L+IL) = 2.*UU(L) - UU(L-IL)
         UV(L+IL) = 2.*UV(L) - UV(L-IL)
         UW(L+IL) = 2.*UW(L) - UW(L-IL)
         UU(L+2*IL) = 3.*UU(L) - 2.*UU(L-IL)
         UV(L+2*IL) = 3.*UV(L) - 2.*UV(L-IL)
         UW(L+2*IL) = 3.*UW(L) - 2.*UW(L-IL)
 6000 CONTINUE
 6001 CONTINUE

      RETURN
      END SUBROUTINE REFDER


      SUBROUTINE TURRSM(NBL,M,U,V,W,RO,A1,A1XA,A1YA,A1ZA,
     + A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,XC,YC,ZC,VOL,D1,D2,
     + D3,F1R,F1RM,F1RN,F1RW,F1E,VIS,EPS2,VIST,OHMI,ICP,JCP,KCP,
     + IMAX,JMAX,KMAX,ICYCLE,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     + IDI1,IDI2,IDI3,IT,IL,IK,IPRINT,PR,PRT,
     + UBI,VBI,WBI,UBJ,VBJ,WBJ,UBK,VBK,WBK,UTI,VTI,WTI,
     + UTJ,VTJ,WTJ,UTK,VTK,WTK,
     + RK,REPS,DDEPS,PTUR,TURLIM,ITURB,FI,F1FI,MAXSB,NSCAL,STRESL,
     + ZZZ,MAXW,IDER)

C ... LAST LINE IS FOR SCALARS PPR 11.2
CIBM  INCLUDE (NS3CH)
CIBM  INCLUDE (WORK3)

      USE CHARACTERS
      USE NS3CO, ONLY : IN, JN, KN
      
      LOGICAL :: STRESL, JEP

      DIMENSION :: U(*),V(*),W(*),RO(*),A1(*),A1XA(*),A1YA(*),A1ZA(*),
     2 A2(*),A2XA(*),A2YA(*),A2ZA(*),A3(*),A3XA(*),A3YA(*),A3ZA(*),
     3 VOL(*),F1R(*),F1RM(*),F1RN(*),F1E(*),F1RW(*),VIS(*),EPS2(*),
     4 VIST(*),ICP(*),JCP(*),KCP(*),OHMI(*),ZC(*),YC(*),XC(*),
     5 D1(*),D2(*),D3(*),PTUR(*),RK(*),REPS(*),DDEPS(*),
     6 UBI(*),VBI(*),WBI(*),UBJ(*),VBJ(*),WBJ(*),UBK(*),VBK(*),WBK(*),
     7 UTI(*),VTI(*),WTI(*),UTJ(*),VTJ(*),WTJ(*),UTK(*),VTK(*),WTK(*),
     8 FI(MAXSB,MAX(1,NSCAL)),F1FI(MAXSB,MAX(1,NSCAL))

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      NTOT    = ISTRID*JSTRID*KSTRID
      JEP     = .TRUE.
      EPS     = 1.E-10
      RELAX   = .7

C ... ADJUST CELLS WITH FRICTION TO THE GRID LEVEL USED

      IDM1    = (IDI1-1)/2**(M-1) + 1
      IDM2    = (IDI2-1)/2**(M-1) + 1
      IDM3    = (IDI3-1)/2**(M-1) + 1


      PRS     = PR/PRT

C **********************************************************************
C                                                                      *
C ... CALCULATION OF TWO-EQUATION TURBULENCE MODELS                    *
C     CALCULATE TURBULENCE VISCOSITIES                                  *
C                                                                      *
C **********************************************************************

C ... VORTICITY, ITS COMPONENTS (=F1RM,F1RN,F1RW) AND STRAIN (=F1R)

c      CALL REFLEG(U,V,W,XC,YC,ZC,OMEGA,IMAX,JMAX,KMAX,IBOT,ITOP,
c     2 JBOT,JTOP,KBOT,KTOP,ICP,JCP,KCP,IDM1,IDM2,IDM3,IN,JN,KN,
c     3 IROTB,IROTT,JROTB,JROTT,KROTB,KROTT,
c     4 UBI,VBI,WBI,UBJ,VBJ,WBJ,UBK,VBK,WBK,UTI,VTI,WTI,
c     5 UTJ,VTJ,WTJ,UTK,VTK,WTK,nbl)
C ... VORFUN calculates vector components, not tensor components.
C ... This should be changed accordingly, but TURRSM is never called!
      CALL VORFUN(OHMI,F1RM,F1RN,F1RW,F1R,F1FI(1,1),F1FI(1,2),F1FI(1,3),
     2 F1FI(1,4),F1FI(1,5),F1FI(1,6),PTUR,U,V,W,RK,REPS,DDEPS,EPS2,VIST,
     2 VIS,A1,A2,A3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     3 D1,D2,D3,IMAX,JMAX,KMAX,IN,JN,KN,ITURB,M,JEP,ZZZ,MAXW,IDER)
      CALL REFDER(F1FI(1,2),F1FI(1,3),F1FI(1,5),IMAX,JMAX,KMAX,IN,JN,KN)

c      DO 9 L = 1,ISTRID*JSTRID*KSTRID
c ... venymaan perustuvat: ehka pitaisi olla neliolliset summat???
c         SXY     = SQRT(FI(L,2)**2+FI(L,3)**2+FI(L,5)**2)/
c     +        SQRT(F1FI(L,2)**2+F1FI(L,3)**2+F1FI(L,5)**2+EPS)
c         RMYT    = SXY/VIS(L)
c         EPS2(L) = RELAX*EPS2(L) + (1.-RELAX)*(1. + RMYT)
c         VIST(L) = RELAX*VIST(L) + (1.-RELAX)*SXY
c         EPS2(L) = MIN(TURLIM,EPS2(L))
c         VIST(L) = MIN(TURLIM*VIS(L),VIST(L))
c 9    CONTINUE
      DO 9 L = 1,ISTRID*JSTRID*KSTRID
c ... k-epsilon mallin tapainen suora yhteys
         AUXI    = (FI(L,1)+FI(L,4)+FI(L,6))**1.5*SQRT(3.*FI(L,4))
         SXY     = .09*AUXI/(4.*REPS(L)*RO(L))
         RMYT    = SXY/VIS(L)
         EPS2(L) = 1. + RMYT
         VIST(L) = SXY
         EPS2(L) = MIN(TURLIM,EPS2(L))
         VIST(L) = MIN(TURLIM*VIS(L),VIST(L))
 9    CONTINUE
 
      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' AFTER BOUSSINESQ'
         CALL PRINYS(3,EPS2C,EPS2, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF

      RETURN
      END SUBROUTINE TURRSM


