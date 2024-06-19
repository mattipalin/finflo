C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_FREVOL(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,
     +     MAXSB,NSCAL,ITURB,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     VOL,XC,YC,ZC,XVL,YVL,ZVL,UROT,VROT,WROT,
     +     ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,XVEL,FRSPRE,WH,
     +     XXI,YXI,XETA,YETA,WFS,IOLD,ICYCLE,ICMAX,A3,DTL,TE,
     *     XHULL,YHULL,ZHULL,
     +     LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,INWH,
     +     JFIRST,NFSD,DTWMAX,DWMV,ICFST,
     +     XCO,YCO,ZCO,UU,VV,WW,PRN,IHULL,FRSDEN,
     +     DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH)

      USE NS3CO, ONLY : IC9, GX, GY, GZ
      
      IMPLICIT NONE

      INTEGER :: MAXSB,NSCAL,ITURB,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,
     +           LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,
     +           KTOPGR,INWH,KBEGIN,KSCAL,ISTR,JSTR,KSTR,L,IA,IB,IM,
     +           JA,JM,KM,IFACE

      INTEGER :: ICON(IC9,*),IOLD,ICYCLE,ICMAX,TE,JFIRST,NFSD,ICFST

      INTEGER :: IHULL(*)

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),VOL(*),
     +        XXI(*),YXI(*),XETA(*),YETA(*),WFS(*),A3(*),DTL(*),
     +        XVL(*),YVL(*),ZVL(*),UROT(*),VROT(*),WROT(*),
     +        XHULL(*),YHULL(*),ZHULL(*),WH(*),UU(*),VV(*),WW(*)

      REAL :: DTWMAX,DWMV,XVEL,FRSPRE,FRSDEN,DWMAX,
     +        FLODWH,WHMAX,WHMIN,SUMDWH

      REAL :: XC(*), YC(*), ZC(*), XCO(*), YCO(*), ZCO(*)

      CHARACTER(LEN=3) :: PRN
C
C ... MIRROR THE VECTORS ON THE REQUIRED  BOUNDARIES OF THE BLOCK
C ... this is for free surface
C ... XVEL is 1 for velocities and -1 for vorticities

C ... IF REYNOLDS STRESSES THERE IS NO NEED TO free surface THEM

      KBEGIN = 1
      KSCAL = NSCAL
      IF(NSCAL /= 0) THEN
         IF(ITURB >= 20) THEN
            KBEGIN = 7
            KSCAL = NSCAL - 6
         ENDIF
      ENDIF

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 13) THEN
C ------------------------------------------------------------------
      IB = ICON(8,L)
      IA = ICON(4,L) - 1 ! Free-surface patches are extended  
      IM = ICON(5,L) + 1
      JA = ICON(6,L) - 1
      JM = ICON(7,L) + 1
      KM = 1
      IFACE = ICON(3,L)
      IF(IFACE == 1 .OR. IFACE == 4) THEN
       IF(IFACE == 4) KM = IMAX
       CALL FS_FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A1X,A1Y,A1Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,JN,KN,IN,IFACE,XVEL,FRSPRE,
     +        WH,XXI,YXI,XETA,YETA,WFS,IOLD,ICYCLE,ICMAX,A3,DTL,IB,TE,
     +        XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,
     +        INWH,JFIRST,NFSD,DTWMAX,DWMV,ICFST,UROT,
     +        VROT,WROT,UU,VV,WW,PRN,FRSDEN,DWMAX,FLODWH,
     +        WHMAX,WHMIN,SUMDWH)

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
       IF(IFACE == 5) KM = JMAX
       CALL FS_FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A2X,A2Y,A2Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,ISTR,KSTR,JSTR,IN,KN,JN,IFACE,XVEL,FRSPRE,
     +        WH,XXI,YXI,XETA,YETA,WFS,IOLD,ICYCLE,ICMAX,A3,DTL,IB,TE,
     +        XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,
     +        INWH,JFIRST,NFSD,DTWMAX,DWMV,ICFST,UROT,
     +        VROT,WROT,UU,VV,WW,PRN,FRSDEN,DWMAX,FLODWH,
     +        WHMAX,WHMIN,SUMDWH)

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
       IF(IFACE == 6) KM = KMAX
       CALL FS_FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A3X,A3Y,A3Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IN,JN,KN,IFACE,XVEL,FRSPRE,
     +        WH,XXI,YXI,XETA,YETA,WFS,IOLD,ICYCLE,ICMAX,A3,DTL,IB,TE,
     +        XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,
     +        INWH,JFIRST,NFSD,DTWMAX,DWMV,ICFST,UROT,
     +        VROT,WROT,UU,VV,WW,PRN,FRSDEN,DWMAX,FLODWH,
     +        WHMAX,WHMIN,SUMDWH)
      ENDIF
      
      ELSEIF(ICON(1,L) == 24) THEN ! Level-set method is not active
C ------------------------------------------------------------------
      IA      = ICON(4,L) - 1 ! Free-surface patches are extended
      IM      = ICON(5,L) + 1
      JA      = ICON(6,L) - 1
      JM      = ICON(7,L) + 1
      KM      = 1
      IFACE = ICON(3,L)
         IF(IFACE == 1 .OR. IFACE == 4) THEN
          IF(IFACE == 4) KM = IMAX
          IF(GX == 0) THEN
          WRITE(*,*)
     *'DIRECTION OF GRAVITY (I) DOES NOT MATCH WITH F.S. BOUNDARY (1-4)'
          STOP
          ENDIF
c         CALL FRESU2(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
c     +        MAXSB,A1X,A1Y,A1Z,
c     +        IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,IFACE,XVEL,
c     +        XXI,YXI,XETA,YETA,WH,WFS,A3,DTL,XC,YC,ZC,TE,IHULL)

         ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
          IF(IFACE == 5) KM = JMAX
          IF(GY == 0) THEN
          WRITE(*,*)
     *'DIRECTION OF GRAVITY (J) DOES NOT MATCH WITH F.S. BOUNDARY (2-5)'
          STOP
          ENDIF
c         CALL FRESU2(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
c     +        MAXSB,A2X,A2Y,A2Z,
c     +        IA,IM,JA,JM,KM,ISTR,KSTR,JSTR,IFACE,XVEL,
c     +        XXI,YXI,XETA,YETA,WH,WFS,A3,DTL,YC,XC,ZC,TE,IHULL)

         ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
          IF(IFACE == 6) KM = KMAX
          IF(GZ == 0) THEN
          WRITE(*,*)
     *'DIRECTION OF GRAVITY (K) DOES NOT MATCH WITH F.S. BOUNDARY (3-6)'
          STOP
          ENDIF
c         CALL FRESU2(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
c     +        MAXSB,A3X,A3Y,A3Z,
c     +        IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IFACE,XVEL,
c     +        XXI,YXI,XETA,YETA,WH,WFS,A3,DTL,ZC,XC,YC,TE,IHULL)
         ENDIF


C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS
1000  CONTINUE
      RETURN
      END SUBROUTINE FS_FREVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FOUVOL(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,MAXSB,
     +     NSCAL,ITURB,A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     VOL,XC,YC,ZC,
     +     ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,XVEL,FRSPRE,WH)

      USE NS3CO, ONLY : IC9
       
      IMPLICIT NONE

      INTEGER :: MAXSB,NSCAL,ITURB,NPATCH,IMAX,JMAX,KMAX,
     +           IN,JN,KN,KBEGIN,KSCAL,ISTR,JSTR,KSTR,L,
     +           IA,IB,IM,JA,JM,KM,IFACE
      INTEGER :: ICON(IC9,*)

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),A1(*),A2(*),A3(*),
     +        A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +        VOL(*),WH(*)
      REAL :: XVEL,FRSPRE

      REAL :: XC(*),YC(*),ZC(*)
      
C
C ... MIRROR THE VECTORS ON THE REQUIRED  BOUNDARIES OF THE BLOCK
C ... this is for free surface
C ... XVEL is 1 for velocities and -1 for vorticities

C ... IF REYNOLDS STRESSES THERE IS NO NEED TO free surface THEM
      KBEGIN = 1
      KSCAL = NSCAL
      IF(NSCAL /= 0) THEN
         IF(ITURB >= 20) THEN
            KBEGIN = 7
            KSCAL = NSCAL - 6
         ENDIF
      ENDIF

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR
      DO 1000 L = 1,NPATCH
C ------------------------------------------------------------------
c      IF(ICON(1,L) == 14) THEN
      IF(ICON(1,L) == 1555) THEN   !ASC0299
C ------------------------------------------------------------------
      IB      = ICON(8,L)
      IA      = ICON(4,L) - 1 ! Free-surface patches are extended
      IM      = ICON(5,L) + 1
      JA      = ICON(6,L) - 1
      JM      = ICON(7,L) + 1
      KM      = 1
      IFACE = ICON(3,L)
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IF(IFACE == 4) KM = IMAX
       CALL WAVEOUT(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A1,A2,A3,A1X,A1Y,A1Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,JN,KN,IN,
     +        IFACE,XVEL,FRSPRE,WH)

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IF(IFACE == 5) KM = JMAX
         CALL WAVEOUT(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A1,A2,A3,A2X,A2Y,A2Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,ISTR,KSTR,JSTR,IN,KN,JN,
     +        IFACE,XVEL,FRSPRE,WH)

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IF(IFACE == 6) KM = KMAX
         CALL WAVEOUT(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        KSCAL,ITURB,MAXSB,A1,A2,A3,A3X,A3Y,A3Z,VOL,XC,YC,ZC,
     +        IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IN,JN,KN,
     +        IFACE,XVEL,FRSPRE,WH)
      ENDIF


C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS

1000  CONTINUE

      RETURN
      END SUBROUTINE FOUVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,
     +     KSCAL,ITURB,MAXSB,NX,NY,NZ,VOL,XC,YC,ZC,
     +     I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IN,JN,KN,IWALL,XVEL,FRSPRE,WH,
     +     XXI,YXI,XETA,YETA,WFS,IOLD,ICYCLE,ICMAX,A3,DTL,IB,TE,
     +     XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,
     +     INWH,JFIRST,NFSD,DTWMAX,DWMV,ICFST,
     +     UROT,VROT,WROT,UU,VV,WW,PRN,FRSDEN,
     +     DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH)          !ASC0798_0200
C
C     TESTING OF FREESURFACE APPROACHES (WRITTEN BY TOM SUNDELL)
C
      USE NS3CO, ONLY : G0

      IMPLICIT NONE

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,
     +           KDIR,KOFF,KB,KB1,KB2

      INTEGER :: IT,IT1,IT2,JL,IL,IOLD,ICYCLE,ICMAX,IDWMAX,
     +           JDWMAX,NFSD,IB,TE,ICFST,JFIRST,TR,KSCAL,ITURB,MAXSB,
     +           LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,
     +           KTOPGR,INWH,KT1,KT2,I,J,K,N,IC1,IC2,NS,IN,JN,KN

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),NX(*), NY(*), NZ(*),
     +        WH(*),VOL(*),
     +        XXI(*),YXI(*),XETA(*),YETA(*),WFS(*),A3(*),DTL(*),
     +        UROT(*),VROT(*),WROT(*),UU(*),VV(*),WW(*),
     +        XHULL(*),YHULL(*),ZHULL(*)
      REAL :: DWH,UCO,VCO,DETJ,DWMAX,WHN,WHS,WHE,WHW,
     +        WHNN,WHNNN,WHWW,WHWWW,WHSS,WHSSS,FLODWH,Q,D,DX1,DY,
     +        DY1,ALFA,DISSCAL,DTW,DTWMAX,US,VS,WS,TEN,
     +        DWMV,XVEL,FRSPRE,FRSDEN,PIEN,GRADEN,RATIO,AAA1,WHEE,DX,
     +        DTMAX,RLIM,WHMAX,WHMIN,SUMDWH
      REAL :: XC(*), YC(*), ZC(*)

      CHARACTER(LEN=3) :: PRN

      DISSCAL = 0.00  ! scale-factor od dissipation in Farmee's method
      PIEN    = 1E-7
      GRADEN  = G0*FRSDEN
      TR      = 0       ! Initialization??
       
C..    Free surface difference operators (at least NFSD=4 works). Note ! simlified
C..    treatment of others (may cause problems in separated flow etc.)
C.....NFSD=1=>Farmer, NFSD=2=>3rd.order upwind/cent.diff, NFSD=3=>4 pnt. upwind
C.....NFSD=4=>3rd.order upwind in both directions
      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine FRESUR !'
      ENDIF
       
      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR
      KT1 = (KW          - 1 + KN)*KSTR
      KT2 = (KW -   KDIR - 1 + KN)*KSTR

C ... If first cycle amd solution from lower grid-leve exists
C.... read WH-file and init. WH-array
         IF (ICYCLE == 0 .AND. IOLD == -1) THEN
         PRINT *, 'INIT. OF WH FROM THE PREVIOUS GRID-LEVEL'
         OPEN(UNIT=399,FILE='WH',STATUS='OLD',ERR=3)
         DO J=J1+1,(J2-1)/2
          DO 4 I=I1+1,(I2-1)/2 
           READ(399,*) WH(I+1+((I2-1)/2+2)*J)
4         CONTINUE
         ENDDO
3        CONTINUE
         CLOSE(399)
         DO 1 J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO 2 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            K=INT((I+1)/2)
            N=INT((J+1)/2)
            
           
           WFS(I+1+(I2+1)*J)=WH(K+1+((I2-1)/2+2)*N)
          
           PDIFF(IT1) = 2.*WFS(I+1+(I2+1)*J)*9810.-PDIFF(IC1)
           PDIFF(IT2) = 2.*PDIFF(IT1)-PDIFF(IC1)
           

2       CONTINUE
1       CONTINUE
C...    init. end
        ENDIF

c         DO  J=J1+1,J2-1
c         JL = (J-1+JN)*JSTR
c         DO  I=I1+1,I2-1
c         write(icycle,*) i,j,WH(I+1+(I2+1)*J)
c         ENDDO
c         ENDDO
c         write(icycle,*) '......................'
       
C......   If firs cycle init WH (INWH=0)          
         IF (INWH == 0) THEN
         DO 10 J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO 20 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2 

            WH(I+1+(I2+1)*J)=(1.5*PDIFF(IT1)-0.5*PDIFF(IT2))/9810.
        
20      CONTINUE
10      CONTINUE
        ELSE IF (INWH == 2) THEN   !add from this line ...
         DO J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2 

           RATIO=(ZC(IC1)-ZC(IC2))/(ZC(IT1)-ZC(IC1))
           AAA1=(ZC(IT1)-ZC(IC1))/(1.+RATIO)
           WH(I+1+(I2+1)*J)=ZC(IT1)-AAA1
c           if(j == j1+1)write(*,'(2I4,6f10.3)')
c          write(666,*)'sab',
c     *     i,j,WH(I+1+(I2+1)*J),ZC(IC2),ZC(IC1),ZC(IT1),RATIO,AAA1
        ENDDO
        ENDDO                          !... to this line 
        ENDIF 
        IF (TE == 0)   GOTO 401  !addition of "Surface tension" 
        DO 100 J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO 200 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2


            RO(IT1) =    RO(IC1)
            E (IT1) =     E(IC1)
            EPS2(IT1) =  EPS2(IC1)
            VIST(IT1) = VIST(IC1)

            RO   (IT2) =    RO(IC2)
            E    (IT2) =     E(IC2)
            EPS2 (IT2) =  EPS2(IC1)
            VIST (IT2) =  VIST(IC1)

            

CCCCCCCCCCCCCCCCCCCCCCCCC++++++++++CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       EVALUATION OF FREE SURFAVE ELEVATION BY THE KINEMATIC
C       BOUNDARY CONDITIO
C       By: Tom Sundell
CCCCCCCCCCCCCCCCCCCCCCCCC++++++++++CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c...    Velocities at the free durface

c       US=U(IC1)-(U(IC2)-U(IC1))*.5*(ZC(IC1)-ZC(IT1))
c     &     /(ZC(IC1)-ZC(IC2))
c       VS=V(IC1)-(V(IC2)-V(IC1))*.5*(ZC(IC1)-ZC(IT1))
c     &     /(ZC(IC1)-ZC(IC2))
c       WS=W(IC1)-(W(IC2)-W(IC1))*.5*(ZC(IC1)-ZC(IT1))
c     &    /(ZC(IC1)-ZC(IC2))

C...    Zero gradient extrapolation ov velocities
        US=(9.*U(IC1)-U(IC2))/8.
        VS=(9.*V(IC1)-V(IC2))/8.
        WS=(9.*W(IC1)-W(IC2))/8.
       
C...  Evaluate contravariant velocities at free surface
         DETJ=XXI(I+1+(I2+1)*J)*YETA(I+1+(I2+1)*J)-
     &     XETA(I+1+(I2+1)*J)*YXI(I+1+(I2+1)*J)
         UCO= US/RO(IC1)*YETA(I+1+(I2+1)*J)/DETJ-  
     &       VS/RO(IC1)*XETA(I+1+(I2+1)*J)/DETJ
         VCO= -US/RO(IC1)*YXI(I+1+(I2+1)*J)/DETJ+
     &       VS/RO(IC1)*XXI(I+1+(I2+1)*J)/DETJ 
c         IF (UCO < 0 .and. J < J2) THEN
c         PRINT *, '*********  WARNING !!!! ************'
c         PRINT *, 'WAVE BREAKING AT NODE:',I,J
c         PRINT *, '************************************'
c         ENDIF

c...    Neighbouring wave heights to north, south, etc.

         WHN=WH(I+1+(I2+1)*J-1)
         WHNN=WH(I+1+(I2+1)*J-2)
         WHNNN=WH(I+1+(I2+1)*J-3)
         WHS=WH(I+1+(I2+1)*J+1)
         WHE=WH(I+1+(I2+1)*J-(I2-I1+1))
         WHW=WH(I+1+(I2+1)*J+(I2-I1+1))
         WHSS=WH(I+1+(I2+1)*J+2)
         WHSSS=WH(I+1+(I2+1)*J+3)
         WHWW=WH(I+1+(I2+1)*J+2*(I2-I1+1))
         WHWWW=WH(I+1+(I2+1)*J+3*(I2-I1+1))
C         WHEE=WH(I+1+(I2+1)*J-2*(I2-I1+1))
         WHEE=WH(MAX(1,I+1+(I2+1)*J-2*(I2-I1+1))) ! Huom, menee nollan alle
C..      Handle the Imin,Imax,Jmin and Jmax cells 
C..      Using ghos cell pressure values to tak into account 
C..      boundary conditions on other patches. 
         IF (J == J1+1) THEN
            WHE=(PDIFF(IC1-JSTR)+PDIFF(IT1-JSTR))/(GRADEN*2.)
            WHEE=(PDIFF(IC1-2*JSTR)+PDIFF(IT1-2*JSTR))/(GRADEN*2.)
            IF (WHEE == 0. .AND. WHE /= 0.) THEN
            WHEE=2*WHE-WH(I+1+(I2+1)*J)
            ENDIF
          ELSEIF (J == J2-1) THEN
            WHW =(PDIFF(IC1+JSTR)+PDIFF(IT1+JSTR))/(GRADEN*2.)
            WHWW=(PDIFF(IC1+2*JSTR)+PDIFF(IT1+2*JSTR))/(GRADEN*2.)
            IF (WHEE == 0. .AND. WHE /= 0.) THEN
            WHWW=2*WHW-WH(I+1+(I2+1)*J)
            ENDIF
         ELSEIF (J == J1+2) THEN
            WHEE=(PDIFF(IC1-2*JSTR)+PDIFF(IT1-2*JSTR))/(GRADEN*2.)
         ELSEIF (J == J2-2) THEN
            WHWW=(PDIFF(IC1+2*JSTR)+PDIFF(IT1+2*JSTR))/(GRADEN*2.)
         ENDIF  ! ASC23.09.98
         IF (I == I1+1) THEN
            WHN =(PDIFF(IC1-ISTR)+PDIFF(IT1-ISTR))/(GRADEN*2.)
            WHNN=(PDIFF(IC1-2*ISTR)+PDIFF(IT1-2*ISTR))/(GRADEN*2.)
c            WHN =WH(I+1+(I2+1)*J)  ! ASC23.09.98
c            WHNN=WHS               ! ASC23.09.98
         ELSEIF (I == I2-1) THEN
         WHS =(PDIFF(IC1+ISTR)+PDIFF(IT1+ISTR))/(GRADEN*2.)
            WHSS=(PDIFF(IC1+2*ISTR)+PDIFF(IT1+2*ISTR))/(GRADEN*2.)
c           WHS =WH(I+1+(I2+1)*J)   ! ASC23.09.98
c           WHSS=WHN                ! ASC23.09.98
         ELSEIF (I == I1+2) THEN
c            WHNN=WHN               ! ASC23.09.98
         ELSEIF (I == I2-2) THEN
c            WHSS=WHS               ! ASC23.09.98
         ENDIF
C         WHEE=WH(I+1+(I2+1)*J-2*(I2-I1+1))
         WHEE=WH(MAX(1,I+1+(I2+1)*J-2*(I2-I1+1))) ! Huom, menee nollan alle

CC.... Evaluate Q=>WFS (velocity used for time-integration of WH) 
 
       IF (NFSD == 1) THEN  ! FARMER
       
           Q= WS/RO(IC1)-
     &        0.5*(WHS-WHN)*UCO
     &       -0.5*(WHW-WHE)*VCO

         
C.... Artificcial dissipation
C.... X-dir
         
        ALFA=DISSCAL*
     *  (ABS( U(IC1+ISTR)/RO(IC1+ISTR)*YETA(I+1+(I2+1)*J+1)
     &  -V(IC1+ISTR)/RO(IC1+ISTR)*XETA(I+1+(I2+1)*J+1))+
     &   ABS( U(IC1)/RO(IC1)*YETA(I+1+(I2+1)*J)-
     &   V(IC1)/RO(IC1)*XETA(I+1+(I2+1)*J)))    
     
        DX  =  ALFA*((WH(I+1+(I2+1)*J) -2*WHS +WHSS)-
     &          (WHN - 2*WH(I+1+(I2+1)*J) + WHS))
        DX1 =  ALFA*((WHS - 2*WHSS + WHSSS)-
     &            (WH(I+1+(I2+1)*J) -2*WHS + WHSS))
   
C.... Y-dir
        ALFA=DISSCAL*
     &  (ABS(V(IC1+JSTR)/RO(IC1+JSTR)* XXI(I+1+(I2+1)*J+1)
     &  -U(IC1+ISTR)/RO(IC1+ISTR)*YXI(I+1+(I2+1)*J+1))+
     &   ABS(U(IC1)/RO(IC1)*XXI(I+1+(I2+1)*J)-
     &   U(IC1)/RO(IC1)*YXI(I+1+(I2+1)*J))) 
   
        DY  =  ALFA*((WH(I+1+(I2+1)*J) -2*WHW +WHWW)-
     &          (WHE - 2*WH(I+1+(I2+1)*J) + WHW))
        DY1 =  ALFA*((WHW - 2*WHWW + WHWWW)-
     &            (WH(I+1+(I2+1)*J) -2*WHW + WHWW))

      D   =  DX1 - DX + DY1 - DY

          ELSEIF (NFSD == 2) THEN  ! 3TH-ORDER SEMISYM/CENT.
            IF (UCO > 0.) THEN
            Q=WS/RO(IC1)-
     &      (WHNN/6.-WHN+WH(I+1+(I2+1)*J)/2.+WHS/3.)*UCO
     &      -0.5*(WHW-WHE)*VCO
          ELSE 
             Q=WS/RO(IC1)+
     &      (WHSS/6.-WHS+WH(I+1+(I2+1)*J)/2.+WHN/3.)*UCO
     &      -0.5*(WHW-WHE)*VCO
          ENDIF
c ....testind dis.term in transversdir
         ALFA=DISSCAL*
     &  (ABS(V(IC1+JSTR)/RO(IC1+JSTR)* XXI(I+1+(I2+1)*J+1)
     &  -U(IC1+ISTR)/RO(IC1+ISTR)*YXI(I+1+(I2+1)*J+1))+
     &   ABS(U(IC1)/RO(IC1)*XXI(I+1+(I2+1)*J)-
     &   U(IC1)/RO(IC1)*YXI(I+1+(I2+1)*J))) 
   
        DY  =  ALFA*((WH(I+1+(I2+1)*J) -2*WHW +WHWW)-
     &          (WHE - 2*WH(I+1+(I2+1)*J) + WHW))
        DY1 =  ALFA*((WHW - 2*WHWW + WHWWW)-
     &            (WH(I+1+(I2+1)*J) -2*WHW + WHWW))

      D   =  DX1 - DX + DY1 - DY
c....testing end
      ELSEIF (NFSD == 3) THEN  ! 4PNT. UPWIND
          Q=WS/RO(IC1)-
     &   (-1./6.)*(-WHNNN+6.*WHNN-15.*WHN+10.*WH(I+1+(I2+1)*J)/2.)*UCO
     &    -0.5*(WHW-WHE)*VCO
      ELSEIF (NFSD == 4) THEN  ! 3TH-ORDER SEMISYM.

      IF (UCO < 0.) THEN
              IF (VCO  < 0.) THEN
               Q=WS/RO(IC1)+
     &         (WHSS/6.-WHS+WH(I+1+(I2+1)*J)/2.+WHN/3.)*UCO
     &         +(WHWW/6.-WHW+WH(I+1+(I2+1)*J)/2.+WHE/3.)*VCO
              ELSE
               Q=WS/RO(IC1)+
     &         (WHSS/6.-WHS+WH(I+1+(I2+1)*J)/2.+WHN/3.)*UCO
     &         -(WHEE/6.-WHE+WH(I+1+(I2+1)*J)/2.+WHW/3.)*VCO
       ENDIF
              ELSE
       IF (VCO  < 0.) THEN
               Q=WS/RO(IC1)-
C     &         (WHNN/2.-2.*WHN+1.5*WH(I+1+(I2+1)*J))*UCO ! pure upwinding
     &         (WHNN/6.-WHN+WH(I+1+(I2+1)*J)/2.+WHS/3.)*UCO
     &         +(WHWW/6.-WHW+WH(I+1+(I2+1)*J)/2.+WHE/3.)*VCO
              ELSE
               Q=WS/RO(IC1)-
C     &         (WHNN/2.-2.*WHN+1.5*WH(I+1+(I2+1)*J))*UCO ! pure upwinding
     &         (WHNN/6.-WHN+WH(I+1+(I2+1)*J)/2.+WHS/3.)*UCO
     &         -(WHEE/6.-WHE+WH(I+1+(I2+1)*J)/2.+WHW/3.)*VCO
       ENDIF
       ENDIF

       ENDIF
c      ENDIF
c         Q  = Q*RO(IC1)  ! Velocities are used directly

           WFS(I+1+(I2+1)*J)=Q
c      ENDIF
200     CONTINUE
100     CONTINUE
c...    Check fs. time-step
          DTW=DTWMAX
cc           IF(ITURB == 0) THEN ! Euler Time Stepping  !!
cc          DO 110 J=J1+1,J2-1
cc          DO 120 I=I1+1,I2-1
cc            IF (DTW > ABS(DWMV/WFS(I+1+(I2+1)*J))) THEN
cc              DTW=ABS(DWMV/WFS(I+1+(I2+1)*J))
cc          ENDIF

cc120   CONTINUE
cc110   CONTINUE
cc        ENDIF !!
        DO 400 J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO 300 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
C           DTW=DTL(IC1)
C           IF (J == J1) THEN
C           DWH=WFS(I+1+(I2+1)*(J+1))*DTW
C          ELSE
c           IF(ITURB > 0) THEN ! Viscous Time Stepping !
c           PRINT *, 'Viscous Timesteping'

C...    Evaluate time step (a funtion of ICYCLE/ICFST and DTWMAX
C...    Specified in the INPUT-file

      IF (ICYCLE < ICFST) THEN
             DTW=DTWMAX*(REAL(ICYCLE)/REAL(ICFST))*
     &       (((U(IC1)**2+V(IC1)**2+W(IC1)**2)/RO(IC1)**2)**(-0.5))
      ELSEIF (ICYCLE >= ICFST) THEN
             DTW=DTWMAX*
     &       (((U(IC1)**2+V(IC1)**2+W(IC1)**2)/RO(IC1)**2)**(-0.5))
      ELSE
              DTW=0.
      ENDIF
             
C...    If last cylce no free-srfa ce update
CCCC           IF (ICYCLE == ICMAX) DTW=0.    ASC25.01.99        
C...    Change of wave height DWH
           DWH=WFS(I+1+(I2+1)*J)*DTW
c           IF(ITURB > 0 .and. J >= JFIRST) THEN ! Update WH 

C...    Limit maximum change of WH with DWMV (specified in INPUT-file)
            IF (ABS(DWH) > DWMV) THEN
               IF (DWH < 0.) THEN
                  DWH=-DWMV
               ELSE
                  DWH=DWMV
               ENDIF    
            ENDIF
         
C ... kokeilee
 
c      DWH = 0. !DWH/(1.+20.*ABS(DWH/(1.E-4+WH(I+1+(I2+1)*J))))

         WH(I+1+(I2+1)*J)=WH(I+1+(I2+1)*J)+DWH
    
C...   Monitoring convergence
      IF (J > JFIRST .AND. J < J2 .AND. 
     &    I > I1 .AND. I < I2) THEN
         IF (ABS(DWH) > ABS(DWMAX)) THEN
           DWMAX=DWH
           IDWMAX=I
           JDWMAX=J
           DTMAX=DTW
         ENDIF
           WHMAX = MAX(WHMAX,WH(I+1+(I2+1)*J))
           WHMIN = MIN(WHMIN,WH(I+1+(I2+1)*J))
      ENDIF
      FLODWH=FLODWH+DWH*A3(IC1)
      SUMDWH=SUMDWH+ABS(DWH)
300   CONTINUE
400   CONTINUE
c         PRINT *, 'WH:',DWMAX,IDWMAX,JDWMAX,DTMAX  
      
      IF (ITURB > 0 ) THEN
C... Extrapolate the innermost wave heights 
        
         CALL LEASTSQR(WH,XC,YC,I1,I2,J1,J2,JFIRST,
     &                 ISTR,JSTR,JN,IN,KT1,KB1,TR,ICYCLE)
      ENDIF
      
C... Ltrans takes care of buttock flow stern transoms not ready June-98 TS !
       CALL LTRANS(XC,YC,ZC,WH,
     *    XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,
     &    KTOPGR,I1,I2,J1,J2,JFIRST,ISTR,JSTR,JN,IN,KT1,KB1,TR) 

401      DO 405 J=J1+1,J2-1
         JL = (J-1+JN)*JSTR
         DO 404 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            
C       update dynamic boundary condition
C       extrapolating pressure to the free surface 
c       "surface tension" to smooth the free surface (only j-direction)... 

c             WHE=WH(I+1+(I2+1)*J-(I2-I1+1))
c             WHW=WH(I+1+(I2+1)*J+(I2-I1+1))
c             DX1=0.5*(YETA(I+1+(I2+1)*J)+YETA(I+1+(I2+1)*J-(I2-I1+1))) !!!YETAcd
c             DX2=0.5*(YETA(I+1+(I2+1)*J)+YETA(I+1+(I2+1)*J+(I2-I1+1)))
c             DY=WHE+DX1*(WHW-WHE)/(DX1+DX2)-WH(I+1+(I2+1)*J)
c             IF (ABS(DY) < PIEN) THEN 
c              JAK=ATAN(1./PIEN)
c             ELSE
c              JAK =ATAN(ABS(0.5*(dx1+dx2)/DY))
c             ENDIF
c             IF (ABS(SIN(3.141593-2*JAK)) < PIEN) THEN
c              RAD=(dx1+dx2)/PIEN
c             ELSE
c             RAD=(dx1+dx2)/(2*SIN(3.141593-2*JAK))
c             ENDIF
c             IF (DY > 0.) THEN
c              RAD=-RAD
c             ENDIF
            
c             IF (TE == 1 ) THEN
c              TEN =0.
c             ELSEIF (J >= JFIRST) THEN
c              TEN=0.035/RAD
c             IF (ABS(TEN) > 30.) THEN
c             TEN=30.*TEN/ABS(TEN)
c              ENDIF
c             ENDIF
            TEN=0.
            RLIM=0.
c...surface tension end 
c pressures on the free surface ....

            PDIFF(IT1) =(2.*(WH(I+1+(I2+1)*J)*9810.+TEN)-PDIFF(IC1))*
     &                   (1.-RLIM) + PDIFF(IT1)*RLIM
            PDIFF(IT2) = (2.*PDIFF(IT1)-PDIFF(IC1))*(1.-RLIM)
     &                   + PDIFF(IT2)*RLIM
cc       write(687,*) i,j,I+1+(I2+1)*J,it1,it2,ic1,WH(I+1+(I2+1)*J)
C   update mometum boundary condition OBS. for some reason U but 
c   returns RM, constant density assumed
C   Flat surface approximation (inviscid) 
            U(IT1) = U(IC1)
            V(IT1) = V(IC1)
            W(IT1) = W(IC1)
            U(IT2) = U(IC2)
            V(IT2) = V(IC2)
            W(IT2) = W(IC2)
C ... Above momentum components (ineffective)
            UU(IT1) = UU(IC1)
            VV(IT1) = VV(IC1)
            WW(IT1) = WW(IC1)
            UU(IT2) = UU(IC2)
            VV(IT2) = VV(IC2)
            WW(IT2) = WW(IC2)
c            write(570,*) i,j,ic1,it1,pdiff(it1),pdiff(ic1),
c     +      pdiff(it2),wh(I+1+(I2+1)*J)
C...   Limiting the maximum Z-velocity to secure convergence in the begining of Hamburger computations. 
C...   Add hoc added June-98 should be removed in future TS. 
c            IF (ITURB == 0) THEN !HAMBURGER
c              IF (ABS(W(IC1)) > 
c     &            ABS(0.3*(U(IC1)**2+V(IC1)**2))**0.5) THEN
c     	        PRINT *, 'limitter',I,J
c	       ENDIF    	   !PRINT		
c	      IF (W(IC1) <  0.) THEN
c              W(IT2) = MAX(W(IC1),-ABS(0.3*(U(IC1)**2+V(IC1)**2)**0.5)) 
c  	      W(IT1) = MAX(W(IC1),-ABS(0.3*(U(IC1)**2+V(IC1)**2)**0.5))
c              ELSE
c		 W(IT2) = MIN(W(IC1),0.3*(U(IC1)**2+V(IC1)**2)**0.5) 
c  	         W(IT1) = MIN(W(IC1),0.3*(U(IC1)**2+V(IC1)**2)**0.5)
c 	      ENDIF
c            ENDIF
C... Hamburger end
404	 CONTINUE
405	 CONTINUE
	  
           



CCCCCCCCCCCCCCCCCCCCCCCCC*********CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (ICYCLE < ICMAX.OR.MOD(ICYCLE,100) == 0) THEN        
C...     Write WH and WFS to file
c        OPEN(UNIT=399,FILE='WH'//PRN,STATUS='UNKNOWN')
        OPEN(UNIT=399,FILE='WH',STATUS='UNKNOWN')
         DO 410 J=J1+1,J2-1
            JL = (J-1+JN)*JSTR
         DO 420 I=I1+1,I2-1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IC1 = 1 + IL + KT1
         WRITE(399,*) WH(I+1+(I2+1)*J),0.5*(XC(IC1)+XC(IT1))
420    CONTINUE
410    CONTINUE
       write(399,*) '____________________________________________'
       CLOSE(UNIT=399)
c       OPEN(UNIT=110,FILE='WFS',STATUS='UNKNOWN')
c         DO 411 J=J1+1,J2-1
c         DO 421 I=I1+1,I2-1
c         WRITE(110,*) WFS(I+1+(I2+1)*J)
c421    CONTINUE
c411    CONTINUE
c       CLOSE(UNIT=110)
           ENDIF       
           
CCCCCCCCCCCCCCCCCCCCCCCCC********CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C...   Mirror variables
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 600 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 500 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            RK  (IT1) =    RK(IC1)
            REPS(IT1) =  REPS(IC1)
            RK  (IT2) =    RK(IC2)
            REPS(IT2) =  REPS(IC2)
 500     CONTINUE
 600  CONTINUE
      ENDIF

      IF(KSCAL >= 1) THEN
      DO NS = 1,KSCAL
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            FI(IT1,NS) =  FI(IC1,NS)
            FI(IT2,NS) =  FI(IC2,NS)
         ENDDO
      ENDDO
      ENDDO
      ENDIF
      RETURN
      END SUBROUTINE FS_FRESUR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_COVOL(ICON,XCO,YCO,ZCO,IMAX,JMAX,KMAX,M,NPATCH,
     &  IN,JN,KN,XXI,YXI,XETA,YETA,XC,YC,ZC,IG,IC,JF,MGM,NBLOCK)

C ... Handle the computing of dx/dxi, dy/dxi, dx/deta, dy/deta

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: IN,JN,KN,MGM,NBLOCK,N,IG1,IC1,
     &           IF1,ISTR,JSTR,KSTR,L,IA,IM,JA,JM,KM,IFACE,M,
     &           IA2,IM2,JA2,JM2

      INTEGER :: ICON(*),NPATCH(*)

      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),
     &           IG(MGM,*),IC(MGM,*),JF(MGM,*)

      REAL    :: XXI(*),YXI(*),XETA(*),YETA(*)

      REAL    :: XCO(*), YCO(*), ZCO(*), XC(*), YC(*), ZC(*)

      DO N = 1,NBLOCK
         
         IG1  = IG(1,N)
         IC1  = IC(1,N)
         IF1  = JF(1,N)
         ISTR = 1
         JSTR =  IMAX(1,N)+2*IN
         KSTR = (JMAX(1,N)+2*JN)*JSTR

      DO L = 1,NPATCH(N)

         IF(ICON(IC1+(L-1)*IC9) == 13) THEN

C ... Extended free surface patch

            CALL PATCHE(ICON(IC1),L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

            KM    = 1  
            IFACE = ICON(3+(L-1)*IC9)
    
            IF(IFACE == 1 .OR. IFACE == 4) THEN
               IF(IFACE == 4) KM = IMAX(1,N)
               CALL FRSCO(IFACE,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,
     &         XCO(IG1),YCO(IG1),ZCO(IG1),XXI(IF1),YXI(IF1),XETA(IF1),
     &         YETA(IF1),XC(IG1),YC(IG1),ZC(IG1))
            ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
               IF(IFACE == 5) KM = JMAX(1,N)
               CALL FRSCO(IFACE,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,
     &         XCO(IG1),YCO(IG1),ZCO(IG1),XXI(IF1),YXI(IF1),XETA(IF1),
     &         YETA(IF1),XC(IG1),YC(IG1),ZC(IG1))
            ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
               IF(IFACE == 6) KM = KMAX(1,N)
               CALL FRSCO(IFACE,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,
     &         XCO(IG1),YCO(IG1),ZCO(IG1),XXI(IF1),YXI(IF1),XETA(IF1),
     &         YETA(IF1),XC(IG1),YC(IG1),ZC(IG1))
            ENDIF
         ENDIF

      ENDDO

      ENDDO

      RETURN
      END SUBROUTINE FS_COVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FRSCO(IFACE,I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,
     & XCO,YCO,ZCO,XXI,YXI,XETA,YETA,XC,YC,ZC)

C ... Compute dx/dxi, dy/dxi, dx/deta, dy/deta

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IFACE,KW,ISTR,JSTR,KSTR,I1,I2,J1,J2,KDIR,
     &           KOFF,KB,KB1,KB2,KT1,KT2,I,J,JL,IL,IC1,IT1

      REAL    :: XXI(*),YXI(*),XETA(*),YETA(*)
      REAL    :: XCO(*),YCO(*),ZCO(*),XC(*),YC(*),ZC(*)

      IF ((IFACE >= 1) .AND. (IFACE <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IFACE >= 4) .AND. (IFACE <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IFACE
         STOP 'Illegal wall specification for subroutine RSCO !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR
      KT1 = (KW          - 1 + KN)*KSTR
      KT2 = (KW -   KDIR - 1 + KN)*KSTR
      JSTR=JSTR
      ISTR=ISTR
      DO 100 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IC1 = 1 + IL + KT1
            IT1 = 1 + IL + KB1
 
c           XXI(I+1+(I2+1)*J)=
c     &       0.25*((XC(IT1+ISTR)+XC(IC1+ISTR)) - 
c     &             (XC(IT1-ISTR)+XC(IC1-ISTR)))
c           YXI(I+1+(I2+1)*J)=
c     &       0.25*((YC(IT1+ISTR)+YC(IC1+ISTR)) - 
c     &	            (YC(IT1-ISTR)+YC(IC1-ISTR)))
c           XETA(I+1+(I2+1)*J)=
c     &       0.25*((XC(IT1 +JSTR)+XC(IC1+JSTR)) - 
c     &	            (XC(IT1-JSTR)+XC(IC1-JSTR)))
c           YETA(I+1+(I2+1)*J)=
c     &       0.25*((YC(IT1+JSTR)+YC(IC1+JSTR)) - 
c     &	             YC((IT1-JSTR)+YC(IC1-JSTR)))

C Hei: kuka tän muutoksen on tehnyt - oletko varma, että xco... päivitetään
c jossei niin nää pysyy vakioina koko ajan

           XXI(I+1+(I2+1)*J)= 0.5*(XCO(IC1+ISTR)-XCO(IC1) + 
     &		           XCO(IC1+ISTR+JSTR)-XCO(IC1+JSTR))
           YXI(I+1+(I2+1)*J)= 0.5*(YCO(IC1+ISTR)-YCO(IC1) + 
     &                YCO(IC1+ISTR+JSTR)-YCO(IC1+JSTR))
	          XETA(I+1+(I2+1)*J)=0.5*(XCO(IC1+JSTR)-XCO(IC1) + 
     &	              XCO(IC1+ISTR+JSTR)-XCO(IC1+ISTR))
	          YETA(I+1+(I2+1)*J)=0.5*(YCO(IC1+JSTR)-YCO(IC1) + 
     &	                   YCO(IC1+ISTR+JSTR)-YCO(IC1+ISTR))
200	CONTINUE
100	CONTINUE

      RETURN
      END SUBROUTINE FRSCO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C	
      SUBROUTINE HYDROSTAT(CS,CMY,TRIM,SINK,FRSVEL,FRSDEN,AREF,CHLREF)
     
      USE NS3CO, ONLY : G0

      IMPLICIT NONE

      REAL :: CS,CMY,TRIM,SINK,FRSVEL,FRSDEN,AREF,CHLREF,DSINK,DTRIM

C...... Only for wigley testing computation of Iwl and Awl should be added
      DSINK=(CS*0.5*AREF*FRSVEL**2)/(0.9375*G0)
      DTRIM=ATAN(CMY*0.5*AREF*FRSVEL**2*CHLREF/(3.369*G0))	!3.369 ?
      TRIM=TRIM+DTRIM
      SINK=SINK+DSINK
C        PRINT *,  DTRIM,TRIM,DSINK,SINK

      RETURN
      END SUBROUTINE HYDROSTAT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LEASTSQR(WH,XC,YC,I1,I2,J1,J2,JFIRST,
     &                    ISTR,JSTR,JN,IN,KT1,KB1,TR,ICYCLE)

      IMPLICIT NONE

      INTEGER :: I1,I2,J1,J2,JFIRST,ISTR,JSTR,JN,IN,KT1,KB1,I,J,K,
     &           IL,IT1,IC1,KL,NL,JL,ICYCLE

C     Least square fitting (line)
C     Tom Sundell 20.05.1997

      REAL :: SX,SXX,SY,SXY,APU1,APU2,S,X
      REAL :: WH(*), XC(*), YC(*)

      INTEGER :: TR

         IF (JFIRST == 0) GOTO 10
         J=JFIRST-1 !TMikkola 090500
1         JL = (J-1+JN)*JSTR
         DO I=I1,I2
2          IL  = JL + (I-1+IN)*ISTR
3          IT1 = 1 + IL + KB1
4          IC1 = 1 + IL + KT1
           X=0   			!length along the line
           S=0
           SX=0
  	   SY=0 
      SXX=0
      SXY=0
c......compute the length along the 'line'
          DO  K=J1,JFIRST-1
            KL = (K-1+JN)*JSTR		
            NL  = KL + (I-1+IN)*ISTR	
            IT1 = 1 + NL + KB1		
            IC1 = 1 + NL + KT1		
            X=SQRT((XC(IC1)-XC(IC1-JSTR))**2+
     &       (YC(IC1)-YC(IC1-JSTR))**2)+X
          ENDDO
c........evaluate A and B
        
           DO K=1,5 
            X=SQRT((XC(IC1+K*JSTR)-XC(IC1+(K-1)*JSTR))**2+
     &              (YC(IC1+K*JSTR)-YC(IC1+(K-1)*JSTR))**2)+X
            S=1+S
            SX=X+SX
            SXX=X**2+SXX
	           SY=WH(I+1+(I2+1)*(J+K))+SY
	           SXY=X*WH(I+1+(I2+1)*(J+K))+SXY
           ENDDO
           apu1=(SXX*SY-SX*SXY)/(S*SXX-SX**2)
           apu2=(S*SXY-SX*SY)/(S*SXX-SX**2)
ccc.......extrapolatin
	  X=0     				!  X-count init. 
          
cc          IF (TR == 1 .AND. I == 70) THEN
cc          PRINT *, 'no leastsq',TR,I
cc          ELSE
          DO  K=J1,JFIRST-1
            KL = (K-1+JN)*JSTR		!1
            NL  = KL + (I-1+IN)*ISTR	!2
           IT1 = 1 + NL + KB1		!3
           IC1 = 1 + NL + KT1		!4
           X=SQRT((XC(IC1)-XC(IC1-JSTR))**2+
     &      (YC(IC1)-YC(IC1-JSTR))**2)+X
           WH(I+1+(I2+1)*K)=apu1+apu2*X
         ENDDO
cc	ENDIF
        ENDDO
10      CONTINUE

      RETURN
      END SUBROUTINE LEASTSQR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WAVEOUT(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,
     +     KSCAL,ITURB,MAXSB,A1,A2,A3,NX,NY,NZ,VOL,XC,YC,ZC,
     +     I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IN,JN,KN,IWALL,XVEL,FRSPRE,WH)

C
C 	TESTING OF FREESURFACE APPROACHES (WRITTEN BY TOM SUNDEL)
C

      IMPLICIT NONE

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF,
     +           KB,KB1,KB2,IN,JN,KN
      INTEGER :: IT,IT1,IT2,JL,IL,
     +           KSCAL,ITURB,MAXSB,I,J,KT1,KT2,IC1,IC2,NS
    
      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),NX(*), NY(*), NZ(*),WH(*),
     +        VOL(*),A1(*),A2(*),A3(*)

      REAL :: FRSPRE, XVEL

      REAL :: XC(*),YC(*),ZC(*)
  
      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine FS_FRESUR !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR
      KT1 = (KW          - 1 + KN)*KSTR
      KT2 = (KW -   KDIR - 1 + KN)*KSTR

c      WN=1
C ... ja initialisoidaan......WH voidaan poistaa jatkossa
         

         DO 10 J=J1-1,J2+1
         JL = (J-1+JN)*JSTR
         DO 20 I=I1-1,I2+1
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2 
c           WH(I+1+(I2+1)*J)=(1.5*PDIFF(IT1)-0.5*PDIFF(IT2))/9810.
          
          
          
20	CONTINUE
10	CONTINUE
	
        DO 100 J=J1+1,J2
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            
            RO(IT1) =    RO(IC1)
            E (IT1) =     E(IC1)
            EPS2(IT1) =  EPS2(IC1)
            VIST(IT1) = VIST(IC1)
	           U(IT1)    = U(IC1)
	           V(IT1)    = V(IC1)
	           W(IT1)    = W(IC1)
            PDIFF(IT1) = PDIFF(IC1)
            RO   (IT2) =    RO(IC2)
            E    (IT2) =     E(IC2)
            EPS2 (IT2) =  EPS2(IC1)
            VIST (IT2) =  VIST(IC1)
	           U(IT2)    = U(IC2)
	           V(IT2)    = V(IC2)
	           W(IT2)    = W(IC2)
            PDIFF(IT2) =PDIFF(IC2)
200     CONTINUE
100     CONTINUE

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
        DO 600 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 500 I=I1,I2
           IL  = JL + (I-1+IN)*ISTR
           IT1 = 1 + IL + KB1
           IT2 = 1 + IL + KB2
           IC1 = 1 + IL + KT1
           IC2 = 1 + IL + KT2

            RK  (IT1) =    RK(IC1)
            REPS(IT1) =  REPS(IC1)
            RK  (IT2) =    RK(IC2)
            REPS(IT2) =  REPS(IC2)
 500    CONTINUE
 600   CONTINUE
      ENDIF
      IF(KSCAL >= 1) THEN
         DO NS = 1,KSCAL
          DO J=J1,J2
           JL = (J-1+JN)*JSTR
           DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            FI(IT1,NS) =  FI(IC1,NS)
            FI(IT2,NS) =  FI(IC2,NS)
           ENDDO
          ENDDO
         ENDDO
        ENDIF

      RETURN
      END SUBROUTINE WAVEOUT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LTRANS(XC,YC,ZC,WH,
     & XHULL,YHULL,ZHULL,LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,
     & KTOPGR,I1,I2,J1,J2,JFIRST,ISTR,JSTR,JN,IN,KT1,KB1,TR)

      IMPLICIT NONE

      INTEGER :: I1,I2,J1,J2,JFIRST,ISTR,JSTR,JN,IN,KT1,KB1,TR,
     &           LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,
     &           KTOPGR,JF,JJ,ITR,L1M0,L1MM,IUP,ILO,IUP2,ILO2,I,KL,
     &           NL,IT1,IC1,J

      REAL :: WH(*),XHULL(*),YHULL(*),ZHULL(*)

      REAL :: XLI,TGFC,X,WH2,WHH,ZH1,ZH2,ZHH,XH1,XH2,XL1,XL2,YH1,
     &        YH2,YL1,YL2,TGHULL,DAU1,DAU2,XHF,AU,AUX,TGWAVE,TGWH1,
     &        ZL1,ZL2

      REAL :: XC(*), YC(*), ZC(*)

C      IF (JFIRST == 0) GOTO 10
         JF=JFIRST 
         JJ=J1+1
         TR=0
         ITR=0
         TGFC=0.9
          IF (JFIRST == 0) THEN
            JF=J1+1
          ENDIF
      L1M0=(LHULL+1)*(MMHUL+1)*MHULL
      L1MM=(LHULL+1)*(MMHUL+1)*(MHULL-1)

         IUP=IBOTGR+1
         ILO=ITOPGR
         IUP2=IBOTGR+1
         ILO2=ITOPGR

         DO 10 I=IBOTGR+1,ITOPGR
           XLI=0.
           X=0.
      	   WH2=WH(I+1+(I2+1)*JJ+(I2-I1+1))
c   extrapolate WH to hull surface...
           WHH=1.5*WH(I+1+(I2+1)*JJ)-0.5*WH2
c   
           ZH1=ZHULL((I-IBOTGR)+L1M0)
           ZH2=ZHULL((I-IBOTGR+1)+L1M0)
c100400           ZHH=0.5*(ZH1+ZH2)
           ZHH=MIN(ZH1,ZH2)

             IF (WHH >= ZHH) THEN
           ZL1=ZH1-ZHULL((I-IBOTGR)+L1MM)
           ZL2=ZH2-ZHULL((I-IBOTGR+1)+L1MM)   !ASC0798
           XH1=XHULL((I-IBOTGR)+L1M0)
           XH2=XHULL((I-IBOTGR+1)+L1M0)
           XL1=XH1-XHULL((I-IBOTGR)+L1MM)
           XL2=XH2-XHULL((I-IBOTGR+1)+L1MM)   !ASC0798
           YH1=YHULL((I-IBOTGR)+L1M0)
           YH2=YHULL((I-IBOTGR+1)+L1M0)
           YL1=YH1-YHULL((I-IBOTGR)+L1MM)
           YL2=YH2-YHULL((I-IBOTGR+1)+L1MM)   !ASC0798
         TGHULL=TGFC*(ZL1+ZL2)/SQRT((YL1+YL2)*(YL1+YL2)+(XL1+XL2)
     &*(XL1+XL2))
         IF(TGHULL < 0.)THEN
           TGHULL=TGHULL/TGFC/TGFC
             IF(I <= ILO2)ILO2=I
             IF(I >= IUP2)IUP2=I
             ITR=1
         ENDIF

c             WRITE(*,1973) 'Using dry transom model. WHH=',
c     &WHH,' ZHH=',ZHH,' TGHULL=',TGHULL,' I=',I
             IF(I <= ILO)ILO=I
             IF(I >= IUP)IUP=I
             TR=1
                 KL = ((J1+1)-1+JN)*JSTR
                 NL = KL+(I-1+IN)*ISTR
                 IT1 = 1+NL+KB1
                 IC1 = 1+NL+KT1
C
           DAU1=SQRT((XC(IC1)-XC(IC1+JSTR))**2+
     &          (YC(IC1)-YC(IC1+JSTR))**2+
     &          (ZC(IC1)-ZC(IC1+JSTR))**2)
           DAU2=SQRT((XC(IC1+JSTR*2)-XC(IC1+JSTR))**2+
     &          (YC(IC1+JSTR*2)-YC(IC1+JSTR))**2+
     &          (ZC(IC1+JSTR*2)-ZC(IC1+JSTR))**2)
           XHF=DAU1*DAU1/(DAU1+DAU2)
           DAU1=SQRT((XC(IC1)-XC(IC1+JSTR))**2+
     &          (YC(IC1)-YC(IC1+JSTR))**2)
           DAU2=SQRT((XC(IC1+JSTR*2)-XC(IC1+JSTR))**2+
     &          (YC(IC1+JSTR*2)-YC(IC1+JSTR))**2)
           AU=DAU1*DAU1/(DAU1+DAU2)
C
	      XLI=XHF         !ASC0798
               DO J=J1+2,JF-1   
                 KL = (J-1+JN)*JSTR
                 NL = KL+(I-1+IN)*ISTR
                 IT1 = 1+NL+KB1
                 IC1 = 1+NL+KT1
                 XLI=SQRT((XC(IC1)-XC(IC1-JSTR))**2+
     &           (YC(IC1)-YC(IC1-JSTR))**2+
     &           (ZC(IC1)-ZC(IC1-JSTR))**2)+XLI
              ENDDO
c	      X=-XHF
              J=J1+1
                KL = (J-1+JN)*JSTR
                NL = KL+(I-1+IN)*ISTR
                IT1 = 1+NL+KB1
                IC1 = 1+NL+KT1
c                X=SQRT((XC(IC1)-XC(IC1-JSTR))**2+
c     &          (YC(IC1)-YC(IC1-JSTR))**2+
c     &          (ZC(IC1)-ZC(IC1-JSTR))**2)+X
	      X=XHF
                AUX=SQRT((YC(IC1)-YC(IC1-JSTR))**2+
     &          (XC(IC1)-XC(IC1-JSTR))**2)
                TGWAVE=(ZC(IC1)-ZC(IC1-JSTR))/AUX
                WH(I+1+(I2+1)*J)=ZHH+(X/XLI)*
     &          (WH(I+1+(I2+1)*(JF-1))-ZHH)
                TGWH1=(WH(I+1+(I2+1)*J)-ZHH)/AU
      IF(TGWH1 > TGHULL)WH(I+1+(I2+1)*J)=ZHH+TGHULL*AU
c           WRITE(*,*)X/XLI,' TGWAVE',TGWAVE,' TGWH',TGWH1
              DO J=J1+2,JF-1   
                KL = (J-1+JN)*JSTR
                NL = KL+(I-1+IN)*ISTR
                IT1 = 1+NL+KB1
                IC1 = 1+NL+KT1
                X=SQRT((XC(IC1)-XC(IC1-JSTR))**2+
     &          (YC(IC1)-YC(IC1-JSTR))**2+
     &          (ZC(IC1)-ZC(IC1-JSTR))**2)+X
                AUX=SQRT((YC(IC1)-YC(IC1-JSTR))**2+
     &          (XC(IC1)-XC(IC1-JSTR))**2)
                TGWAVE=(ZC(IC1)-ZC(IC1-JSTR))/AUX
                WH(I+1+(I2+1)*J)=ZHH+(X/XLI)*
     &          (WH(I+1+(I2+1)*(JF-1))-ZHH)
                TGWH1=(WH(I+1+(I2+1)*J)-WH(I+1+(I2+1)*(J-1)))/AUX
      IF(TGWH1 > TGHULL)WH(I+1+(I2+1)*J)=WH(I+1+(I2+1)*(J-1))+
     *TGHULL*AUX
cc           WRITE(6,*)X/XLI,WH(I+1+(I2+1)*J),WH(I+1+(I2+1)*JF),ZHH
c           WRITE(6,*)X/XLI,' TGWAVE',TGWAVE,' TGWH',TGWH1
              ENDDO
              DO J=JF,JF+5
                KL = (J-1+JN)*JSTR
                NL = KL+(I-1+IN)*ISTR
                IT1 = 1+NL+KB1
                IC1 = 1+NL+KT1
                AUX=SQRT((YC(IC1)-YC(IC1-JSTR))**2+
     &          (XC(IC1)-XC(IC1-JSTR))**2)
                TGWAVE=(ZC(IC1)-ZC(IC1-JSTR))/AUX
                TGWH1=(WH(I+1+(I2+1)*J)-WH(I+1+(I2+1)*(J-1)))/AUX
      IF(TGWH1 > TGHULL)WH(I+1+(I2+1)*J)=WH(I+1+(I2+1)*(J-1))+
     *TGHULL*AUX
c           WRITE(6,*)(J-JF)*1.,' TGWAVE',TGWAVE,' TGWH',TGWH1
              ENDDO
            ENDIF       
 10       CONTINUE
    
       IF(TR == 1) THEN
          WRITE(*,1974) '    Using dry transom model from I=',
     *    ILO,' to I=',IUP
          IF(ITR == 1)
     *    WRITE(*,1974) '    Hull tangent  < 0  between   I=',ILO2,
     *    ' to I=',IUP2
       ENDIF
 1974  FORMAT(15X,(A),X,I3,(A),I3,(A),I3,(A),I3)
 1973  FORMAT((A),X,F6.4,(A),X,F6.4,(A),X,F6.4,(A),I3)

      RETURN      
      END SUBROUTINE LTRANS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETWAVEH(WAVEH,WH,IHF,ICON,NPATCH,IMAX,
     +     JMAX,KMAX,IN,JN,KN,ICYCLE,IPRINT)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,
     +           IA,IM,JA,JM,IF2,KM,ILMO,IFACE,L,ICYCLE,ISTRID,
     +           JSTRID,KSTRID,IPRINT

      INTEGER :: IHF(*),ICON(IC9,*)

      REAL :: WH(*)

      REAL :: WAVEH(*)

      ISTR   = 1
      JSTR   = IMAX +2*IN
      KSTR   = (JMAX+2*JN)*JSTR
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      
      DO L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 13) THEN
C ------------------------------------------------------------------

      IF2   = IHF(ICON(2,L))       ! Starting address of this patch
      IFACE = ICON(3,L)
      IA    = ICON(4,L) !- 1        ! Indeces are extended for FS_SET 
      IM    = ICON(5,L) !+ 1        ! (mersu)
      JA    = ICON(6,L) !- 1
      JM    = ICON(7,L) !+ 1
      IF(IFACE == 2 .OR. IFACE == 5) THEN
        JA  = ICON(4,L) !- 1
        JM  = ICON(5,L) !+ 1
        IA  = ICON(6,L) !- 1
        IM  = ICON(7,L) !+ 1
      ENDIF

      KM   = 1
      ILMO = 1
      IF(IFACE > 3) ILMO = -1
       
      IF(IFACE == 1 .OR. IFACE == 4) THEN
        IF(IFACE == 4) KM = IMAX
        CALL FS_SET(WAVEH(IF2),WH,IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,
     +       JSTRID,KSTRID,ISTRID,JN,KN,IN,IFACE,ILMO,ICYCLE,IPRINT)
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
        IF(IFACE == 5) KM = JMAX
        CALL FS_SET(WAVEH(IF2),WH,IA,IM,JA,JM,KM,ISTR,KSTR,JSTR,
     +       KSTRID,ISTRID,JSTRID,KN,IN,JN,IFACE,ILMO,ICYCLE,IPRINT)
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
        IF(IFACE == 6) KM = KMAX
        CALL FS_SET(WAVEH(IF2),WH,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,
     +       ISTRID,JSTRID,KSTRID,IN,JN,KN,IFACE,ILMO,ICYCLE,IPRINT)
      ENDIF
      ENDIF ! ICON(1,L) == 13 (Free surface
      ENDDO

      RETURN
      END SUBROUTINE FS_SETWAVEH
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETWAVES(WAVEH,WAVES,IHF,ICON,NPATCH,IMAX,
     +     JMAX,KMAX,IN,JN,KN,ICYCLE,IPRINT)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,
     +           IA,IM,JA,JM,IF2,KM,ILMO,IFACE,L,ICYCLE,ISTRID,
     +           JSTRID,KSTRID,IPRINT

      INTEGER :: IHF(*), ICON(IC9,*)

      REAL :: WAVEH(*), WAVES(*)

      ISTR   = 1
      JSTR   = IMAX +2*IN
      KSTR   = (JMAX+2*JN)*JSTR
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      
      DO L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 13) THEN
C ------------------------------------------------------------------

      IF2   = IHF(ICON(2,L))       ! Starting address of this patch
      IFACE = ICON(3,L)
      IA    = ICON(4,L)         
      IM    = ICON(5,L)         
      JA    = ICON(6,L) 
      JM    = ICON(7,L) 

      IF(IFACE == 2 .OR. IFACE == 5) THEN
        JA  = ICON(4,L)
        JM  = ICON(5,L) 
        IA  = ICON(6,L) 
        IM  = ICON(7,L) 
      ENDIF

      KM   = 1
      ILMO = 1
      IF(IFACE > 3) ILMO = -1
          
      IF(IFACE == 1 .OR. IFACE == 4) THEN
        IF(IFACE == 4) KM = IMAX
        CALL VERPAT(WAVES(IF2),WAVEH(IF2),IA,IM,JA,JM)
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
        IF(IFACE == 5) KM = JMAX
        CALL VERPAT(WAVES(IF2),WAVEH(IF2),IA,IM,JA,JM)
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
        IF(IFACE == 6) KM = KMAX
        CALL VERPAT(WAVES(IF2),WAVEH(IF2),IA,IM,JA,JM)
      ENDIF
      ENDIF ! ICON(1,L) == 13 (Free surface
      ENDDO

      RETURN
      END SUBROUTINE FS_SETWAVES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SET(WAVEH,WH,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,
     +     ISTRID,JSTRID,KSTRID,IN,JN,KN,IFACE,ILMO,ICYCLE,IPRINT)

      USE NS3CO, ONLY : LN
      
C ... Set a values from a 1-ghostcell array to ordinary patch arrays

      IMPLICIT NONE

      INTEGER :: IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IFACE,ILMO,ICYCLE,
     +           LSTRID,NFL,I,J,NN,NR,N1,IN,JN,KN,
     +           ISTRID,NI,NH,NP,JSTRID,KSTRID,IPRINT

      REAL :: WAVEH(*)

      REAL :: WH(*)

      LSTRID  = KSTR*ILMO
      IF(ILMO ==  1) NFL = 0
      IF(ILMO == -1) NFL = -LSTRID

*      DO J = JA,JM
      DO J = JA-1,JM+1
      NN = (J-JA)*(IM-IA+1) - IA + 1
      NR = (LN+J-JA)*(IM-IA+1+2*LN) - IA + LN + 1
      N1 = (JN+J-1)*ISTRID + IN
*      DO I = IA,IM
      DO I = IA-1,IM+1
         NI = I + N1                         ! PATCH VELOCITY INDEX
         NH = I + NN                         ! New patch index
         NP = I + NR                         ! Patch index with LN ghost cells
         WAVEH(NI) = WH(NH)  ! + WAVEH(NI))*.5
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FS_SET
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETWAVEHBACK(WAVEH,WH,IHF,ICON,NPATCH,IMAX,
     +     JMAX,KMAX,IN,JN,KN,ICYCLE)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,
     +           IA,IM,JA,JM,IF2,KM,ILMO,IFACE,L,ICYCLE,ISTRID,
     +           JSTRID,KSTRID

      INTEGER :: IHF(*),ICON(IC9,*)

      REAL :: WH(*)

      REAL :: WAVEH(*)

      ISTR   = 1
      JSTR   = IMAX +2*IN
      KSTR   = (JMAX+2*JN)*JSTR
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      
      DO L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 13) THEN
C ------------------------------------------------------------------

      IF2   = IHF(ICON(2,L))       ! Starting address of this patch
      IFACE = ICON(3,L)
      IA    = ICON(4,L) !- 1        ! Indeces are lowered for FS_SETBACK
      IM    = ICON(5,L) !+ 1        ! (mersu)
      JA    = ICON(6,L) !- 1
      JM    = ICON(7,L) !+ 1
      IF(IFACE == 2 .OR. IFACE == 5) THEN
        JA  = ICON(4,L) !- 1
        JM  = ICON(5,L) !+ 1
        IA  = ICON(6,L) !- 1
        IM  = ICON(7,L) !+ 1
      ENDIF

      KM   = 1
      ILMO = 1
      IF(IFACE > 3) ILMO = -1
     
      IF(IFACE == 1 .OR. IFACE == 4) THEN
        IF(IFACE == 4) KM = IMAX
        CALL FS_SETBACK(WAVEH(IF2),WH,IA,IM,JA,JM,KM,JSTR,
     +       KSTR,ISTR,JSTRID,KSTRID,ISTRID,JN,KN,IN,IFACE,ILMO,ICYCLE)
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
        IF(IFACE == 5) KM = JMAX
        CALL FS_SETBACK(WAVEH(IF2),WH,IA,IM,JA,JM,KM,ISTR,
     +       KSTR,JSTR,KSTRID,ISTRID,JSTRID,KN,IN,JN,IFACE,ILMO,ICYCLE)
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
        IF(IFACE == 6) KM = KMAX
        CALL FS_SETBACK(WAVEH(IF2),WH,IA,IM,JA,JM,KM,ISTR,
     +       JSTR,KSTR,ISTRID,JSTRID,KSTRID,IN,JN,KN,IFACE,ILMO,ICYCLE)
      ENDIF
      ENDIF ! ICON(1,L) == 13 (Free surface
      ENDDO

      RETURN
      END SUBROUTINE FS_SETWAVEHBACK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETBACK(WAVEH,WH,IA,IM,JA,JM,KM,ISTR,JSTR,
     +     KSTR,ISTRID,JSTRID,KSTRID,IN,JN,KN,IFACE,ILMO,ICYCLE)

C ... Set a values from ordinary patch arrays to a 1-ghostcell array

      IMPLICIT NONE

      INTEGER :: IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IFACE,ILMO,ICYCLE,
     +           LSTRID,I,J,NN,N1,IN,JN,KN,ISTRID,
     +           NI,NH,JSTRID,KSTRID

      REAL :: WAVEH(*)

      REAL :: WH(*)

      LSTRID  = KSTR*ILMO
*      DO J = JA,JM
      DO J = JA-1,JM+1
      NN      = (J-JA)*(IM-IA+1) - IA + 1
      N1      = (JN+J-1)*ISTRID + IN
*      DO I = IA,IM
      DO I = IA-1,IM+1
         NI     = I  + N1                  ! PATCH VELOCITY INDEX
         NH     = I  + NN                  ! New patch index
         WH(NH) = WAVEH(NI) 
c         write(689,*) i,j,nh,ni,wh(nh)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FS_SETBACK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETWAVEHVALUE(WAVEH,WHEIGHT,IHF,ICON,NPATCH,IMAX,
     +     JMAX,KMAX,IN,JN,KN,ICYCLE,NGL,XCP,YCP,ZCP,XBULB,YBULB,ABULB)

      USE NS3CO, ONLY : LN, IC9

      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,
     +           IA,IM,JA,JM,IF2,KM,ILMO,IFACE,L,ICYCLE,ISTRID,
     +           JSTRID,KSTRID,NGL

      INTEGER :: IHF(*),ICON(IC9,*)

      REAL :: WHEIGHT,XBULB,YBULB,ABULB

      REAL :: XCP(*),YCP(*),ZCP(*)

      REAL :: WAVEH(*)

      ISTR    = 1
      JSTR    = IMAX +2*IN
      KSTR    = (JMAX+2*JN)*JSTR
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      
      DO L = 1,NPATCH

         IF(ICON(1,L) == 13) THEN

            IF2   = IHF(ICON(2,L)) ! Starting address of this patch
            IFACE = ICON(3,L)
            IA    = ICON(4,L) !- 1    ! Indeces are lowered for FS_SETBACK
            IM    = ICON(5,L) !+ 1    ! (mersu)
            JA    = ICON(6,L) !- 1
            JM    = ICON(7,L) !+ 1
            KM    = KMAX
            ILMO  = -1

            CALL FS_SETVALUE(WAVEH(IF2),WHEIGHT,IA,IM,JA,JM,KM,
     +           ISTR,KSTR,ISTR,JSTRID,KSTRID,ISTRID,IN,JN,KN,
     +           IFACE,ILMO,ICYCLE,NGL,XCP(IF2),YCP(IF2),ZCP(IF2),
     +           XBULB,YBULB,ABULB)
         ENDIF                  ! ICON(1,L) == 13 (Free surface

      ENDDO

      RETURN
      END SUBROUTINE FS_SETWAVEHVALUE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_SETVALUE(WAVEH,WHEIGHT,IA,IM,JA,JM,KM,ISTR,JSTR,
     +     KSTR,ISTRID,JSTRID,KSTRID,IN,JN,KN,IFACE,ILMO,ICYCLE,NGL,
     +     XCP,YCP,ZCP,XBULB,YBULB,ABULB)

      USE NS3CO, ONLY : LN

C ... Set a values from ordinary patch arrays to a 1-ghostcell array

      IMPLICIT NONE

      INTEGER :: IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IFACE,ILMO,ICYCLE,
     +           I,J,NR,IN,JN,KN,ISTRID,
     +           NP,JSTRID,KSTRID,NGL

      REAL :: WAVEH(*)

      REAL :: WHEIGHT, COR, XBULB, YBULB, ABULB

      REAL :: XCP(*), YCP(*), ZCP(*)

      DO J = JA-2,JM+2
         NR = (LN+J-JA)*(IM-IA+1+2*LN) - IA + LN + 1
         DO I = IA-2,IM+2
            NP = I + NR  ! Patch index with LN ghost cells
            COR = ABULB**(SQRT((XCP(NP)-XBULB)**2+(YCP(NP)-YBULB)**2)) 
            WAVEH(NP) = WHEIGHT*COR
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FS_SETVALUE
C
C ----------------------------------------------------------------------
C -------------- Free surface by ESa and PR 31.1.2000 ------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE WSHAPE(XINI,YINI,ZINI,  
     &           XNEW,YNEW,ZNEW,IMAX,JMAX,KMAX,
     &           IHF,ICON,NPATCH,WAVES,ICYCLE,IPRINT,IREPEA,NREPEA,NGL)

C ... Subroutine modifies block face shapes and redistributes
C ... the volume grid nodes along the grid lines.

      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      REAL, DIMENSION(*):: XINI, YINI, ZINI
      REAL, DIMENSION(*):: XNEW, YNEW, ZNEW

      REAL, DIMENSION(*):: WAVES

      INTEGER :: ICON(IC9,NPATCH), IHF(*),IREPEA(*), NREPEA(*)

      INTEGER :: IMAX, JMAX, KMAX, I, L, IPL, ICYCLE, IPRINT
      INTEGER :: JA, JM, KA, KM, IFACE, NPATCH, NNODE, NGL

      NNODE = (IMAX+2*IN)*(JMAX+2*JN)*(KMAX+2*KN)

      DO I=1,NNODE
         XNEW(I) = XINI(I)
         YNEW(I) = YINI(I)
         ZNEW(I) = ZINI(I)
      ENDDO

      DO L = 1,NPATCH

         IF(ICON(1,L) == 13) THEN  ! Free surface

            IPL   = ICON(2,L)      ! Proces local patch number
            IFACE = ICON(3,L)
            JA    = ICON(4,L)
            JM    = ICON(5,L)
            KA    = ICON(6,L)
            KM    = ICON(7,L)
            IF(IFACE == 2 .OR. IFACE == 5) THEN
               JA    = ICON(6,L)
               JM    = ICON(7,L)
               KA    = ICON(4,L)
               KM    = ICON(5,L)
            ENDIF
              
            CALL REDIST(XNEW,YNEW,ZNEW,XINI,YINI,ZINI,
     &           IMAX,JMAX,KMAX,IFACE,JA,JM,KA,KM,
     &           WAVES,IHF(IPL),ICYCLE,IPRINT,IREPEA(1),NREPEA(1),NGL)

         ENDIF

      ENDDO

      IREPEA(1) = IREPEA(1) + 1

      RETURN           
      END SUBROUTINE WSHAPE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REDIST(XNEW,YNEW,ZNEW,XINI,YINI,ZINI,
     &           IMAX,JMAX,KMAX,IFACE,JA,JM,KA,KM,
     &           WAVES,IHF,ICYCLE,IPRINT,IPRIFS,NREPEA1,NGL)

C ... Redistribution of the grid nodes along grid lines

      USE NS3CO, ONLY : IN, JN, KN, LN, GX, GY, GZ, GROUND
      USE CONSTANTS, ONLY : PII

      IMPLICIT NONE

      INTEGER, PARAMETER :: MNODES = 1025  ! Maximum number of points per line

      REAL, DIMENSION(*)      :: XNEW, YNEW, ZNEW,
     &                                       XINI, YINI, ZINI, WAVES
      REAL, DIMENSION(MNODES) :: XN, YN, ZN, RR, SS

      REAL :: WLEVEL,RLEN, SLEN, TLEN
      REAL :: DR, DS, RATIO, A, B, U, SF, SFACC
      REAL :: DS1, DS2

      REAL :: AX, AY, AZ, BX, BY, BZ

      INTEGER, DIMENSION(MNODES) :: INODES
      INTEGER :: IPRIFS, INTO, NREPEA1

      LOGICAL :: FOUND

      INTEGER :: IMAX, JMAX, KMAX, ISTRID, JSTRID, KSTRID
      INTEGER :: ILE, IFACE, ISTEP, ISLAB 
      INTEGER :: IS, IE, JA, JM, KA, KM 
      INTEGER :: IPRINT, I, J, K, N, IP, NN, IQ
      INTEGER :: L1, L2, NP, NH, ICYCLE, IHF, NGL

      REAL :: VP, VS, GX1, GY1, GZ1

      INTO   = 1
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ILE    = ISTRID*JSTRID

      IF(MAX(ISTRID,JSTRID,KSTRID) > MNODES) THEN
         WRITE(*,*) 'REDIST: CONTACT PATRIK. EXITING ...'
         STOP
      ENDIF

      GX1 = GX/SQRT(GX**2+GY**2+GZ**2)
      GY1 = GY/SQRT(GX**2+GY**2+GZ**2)
      GZ1 = GZ/SQRT(GX**2+GY**2+GZ**2)

      SELECT CASE(IFACE)
         CASE(1)
            ISTEP = -1
            ISLAB = -1
            IS    =  IMAX + 1
            IE    =  1
         CASE(2)
            ISTEP = -1
            ISLAB = -ISTRID
            IS    =  JMAX + 1
            IE    =  1
         CASE(3)
            ISTEP = -1
            ISLAB = -ISTRID*JSTRID
            IS    =  KMAX + 1
            IE    =  1
         CASE(4)
            ISTEP =  1
            ISLAB =  1
            IS    =  1
            IE    =  IMAX + 1
         CASE(5)
            ISTEP =  1
            ISLAB =  ISTRID
            IS    =  1
            IE    =  JMAX + 1
         CASE(6)
            ISTEP =  1
            ISLAB =  ISTRID*JSTRID
            IS    =  1
            IE    =  KMAX + 1
      END SELECT


      IF(ICYCLE+1 == IPRINT) THEN  ! Mersu

         WRITE(46,*)      
         WRITE(46,*) ' From WSHAPE: Surface height. ICYCLE = ',ICYCLE
         WRITE(46,'(" imin,imax,jmin,jmax",4I5)') JA,JM,KA,KM

         DO K = KA,KM+1
            DO J = JA,JM+1

               NN = (LN+K-KA)*(JM-JA+1+2*LN) - JA + LN + IHF
               NH = NN + J

               WRITE(46,'(2I4,2I6,G14.6)') J,K,NH-IHF+1,NH,WAVES(NH)
            ENDDO
         ENDDO
         WRITE(46,*)
      ENDIF

      DO K = KA,KM+1
         DO J = JA,JM+1

           NN = (LN+K-KA)*(JM-JA+1+2*LN) - JA + LN + IHF
           NH = NN + J 

C ... Total length of the original grid line and the length after
C ... the modification.

            WLEVEL = WAVES(NH) + GROUND
            RLEN   = 0.0
            SLEN   = 0.0
            NP     = 1
            RR(NP) = RLEN
            FOUND  = .FALSE.
            I      = IE

            SELECT CASE(IFACE)
               CASE(1,4)
                  L1 = I+IN+(J-1+JN)*ISTRID+(K-1+KN)*ILE
               CASE(2,5)
                  L1 = K+IN+(I-1+JN)*ISTRID+(J-1+KN)*ILE
               CASE(3,6)
                  L1 = J+IN+(K-1+JN)*ISTRID+(I-1+KN)*ILE
            END SELECT

         VP = -(GX1*XNEW(L1) + GY1*YNEW(L1) + GZ1*ZNEW(L1))

C ... Check if the free surface lies higher than the background grid.
         
         IF(VP < WLEVEL) THEN 
            IF(IPRIFS < NREPEA1) THEN 
               WRITE(13,11) ICYCLE,IE,J,K,VP,WLEVEL
            ELSE IF(IPRIFS == NREPEA1) THEN
               IF(INTO == 1) WRITE(13,12) IPRIFS
               INTO = INTO + 1
            ENDIF

            XNEW(L1) = XNEW(L1) - GX1*(WAVES(NH) + GROUND - VP)
            YNEW(L1) = YNEW(L1) - GY1*(WAVES(NH) + GROUND - VP)
            ZNEW(L1) = ZNEW(L1) - GZ1*(WAVES(NH) + GROUND - VP)
              
            WLEVEL   = WAVES(NH) + GROUND

         ENDIF

 11      FORMAT('  Pinta vuotaa yli',I6,1X,3I5,2X,5G11.4)
 12      FORMAT('  Pinta vuosi yli - repeated',I3,' times. Silence')

            DO I=IS,IE-ISTEP,ISTEP

               SELECT CASE(IFACE)
                  CASE(1,4)
                     L1 = I+IN+(J-1+JN)*ISTRID+(K-1+KN)*ILE
                  CASE(2,5)
                     L1 = K+IN+(I-1+JN)*ISTRID+(J-1+KN)*ILE
                  CASE(3,6)
                     L1 = J+IN+(K-1+JN)*ISTRID+(I-1+KN)*ILE
               END SELECT
               
               L2 = L1 + ISLAB

               INODES(NP)   = L1
               INODES(NP+1) = L2

               NP = NP + 1

               DR = SQRT((XINI(L1)-XINI(L2))**2
     &                 + (YINI(L1)-YINI(L2))**2
     &                 + (ZINI(L1)-ZINI(L2))**2)
               DS = SQRT((XNEW(L1)-XNEW(L2))**2
     &                 + (YNEW(L1)-YNEW(L2))**2
     &                 + (ZNEW(L1)-ZNEW(L2))**2)

               VS = -(GX1*XNEW(L1) + GY1*YNEW(L1) + GZ1*ZNEW(L1))
               VP = -(GX1*XNEW(L2) + GY1*YNEW(L2) + GZ1*ZNEW(L2))

C ... Find the intersection point of the grid line and the free surface.
               
               IF(.NOT.FOUND) THEN

                  IF((VS-WLEVEL)*(VP-WLEVEL) <= 0.0) THEN
                     DS = (WLEVEL-VS)/(VP-VS)*DS
                     FOUND = .TRUE.
                  ENDIF

               ELSE
                  DS = 0.0
               ENDIF
               
               RLEN   = RLEN + DR
               SLEN   = SLEN + DS
               RR(NP) = RLEN

            ENDDO

            
            RATIO  = SLEN/RLEN

            IF(NP <= 5) THEN  ! Old clustering (linear scaling)

               DO IP=1,NP
                  RR(IP) = RATIO*RR(IP)
               ENDDO

            ELSE  ! New clustering (TANH with first cell size in both ends)

               DS1 = RR(2)/SLEN
               DS2 = (RR(NP)-RR(NP-1))/SLEN

C ... TANH distribution [0,1]

               SF = 1.0
               SFACC = 1.0e-8
               N  = NP - 1
               A  = SQRT(DS2)/SQRT(DS1)
               B  = 1.0/(N*SQRT(DS2*DS1))

               CALL SFACTOR(DS1,DS2,B,SF,SFACC)

               SS(1) = 0.0

               DO I=0,N

                  U = 0.5*(1.0+(TANH(SF*(REAL(I)/REAL(N)-0.5)))/
     &                         (TANH(0.5*SF)))    

                  SS(I+1) = U/(A+(1.0-A)*U)

               ENDDO

               RR(1) = 0.0
               DO I=2,NP
                  RR(I) = SS(I)*SLEN
               ENDDO

            ENDIF


C ... Redistribution of the nodes along a grid line (original)

            TLEN = 0.0
            IP   = 1
            IQ   = 1

C ... Just in case we can't find the last point on the free surface.

            XN(NP) = XNEW(INODES(NP))           
            YN(NP) = YNEW(INODES(NP))           
            ZN(NP) = ZNEW(INODES(NP))           
            
            
 33         CONTINUE

            IP = IP + 1

            IF(IP > NP) GOTO 35

               L1 = INODES(IP-1)
               L2 = INODES(IP)
            
               DS = SQRT((XNEW(L1)-XNEW(L2))**2
     &                + (YNEW(L1)-YNEW(L2))**2
     &                + (ZNEW(L1)-ZNEW(L2))**2)
               TLEN = TLEN + DS
 
 34            CONTINUE

               IQ = IQ + 1

            IF(IQ > NP) GOTO 35
               IF(TLEN >= RR(IQ)) THEN
                  RATIO  = (RR(IQ)-(TLEN-DS))/DS
                  XN(IQ) = XNEW(L1) + RATIO*(XNEW(L2)-XNEW(L1)) 
                  YN(IQ) = YNEW(L1) + RATIO*(YNEW(L2)-YNEW(L1)) 
                  ZN(IQ) = ZNEW(L1) + RATIO*(ZNEW(L2)-ZNEW(L1)) 
                  GOTO 34
               ELSE
                  IQ = IQ - 1 
                  GOTO 33
               ENDIF

 35         CONTINUE

            DO IP=2,NP
               L1 = INODES(IP)
               XNEW(L1) = XN(IP)
               YNEW(L1) = YN(IP)
               ZNEW(L1) = ZN(IP)
            ENDDO


C ... A new idea: project the distribution linearly from a line
C ... connecting the end points of the new grid line. Don't remove  
C ... the original code above. The intersection point of the free
C ... surface and the grid line is utilized below.
    
            L1 = INODES(NP)
            L2 = INODES(1)

            BX = XNEW(L1)-XNEW(L2)
            BY = YNEW(L1)-YNEW(L2)
            BZ = ZNEW(L1)-ZNEW(L2)

            SLEN = SQRT(BX**2 + BY**2 + BZ**2)

            RR(1) = 0.0

            DO I=2,NP
               RR(I) = SS(I)*SLEN
            ENDDO

            TLEN = 0.0
            IP   = 1
            IQ   = 1
            
 37         CONTINUE

            IP = IP + 1

            IF(IP > NP) GOTO 39
 
               L1 = INODES(IP-1)
               L2 = INODES(IP)

               AX = XINI(L2)-XINI(L1)
               AY = YINI(L2)-YINI(L1)
               AZ = ZINI(L2)-ZINI(L1)

               DS = (AX*BX+AY*BY+AZ*BZ)/SLEN

               TLEN = TLEN + DS

 38            CONTINUE

               IQ = IQ + 1

            IF(IQ > NP) GOTO 39
               IF(TLEN >= RR(IQ)) THEN
                  RATIO  = (RR(IQ)-(TLEN-DS))/DS
                  XN(IQ) = XINI(L1) + RATIO*(XINI(L2)-XINI(L1)) 
                  YN(IQ) = YINI(L1) + RATIO*(YINI(L2)-YINI(L1)) 
                  ZN(IQ) = ZINI(L1) + RATIO*(ZINI(L2)-ZINI(L1))
                  GOTO 38
               ELSE
                  IQ = IQ - 1
                  GOTO 37
               ENDIF

 39         CONTINUE

            DO IP=2,NP
               L1 = INODES(IP)
               XNEW(L1) = XN(IP)
               YNEW(L1) = YN(IP)
               ZNEW(L1) = ZN(IP)
            ENDDO
           
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE REDIST
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SFACTOR(X1,X2,B,SF,SFACC)

C ... Subroutine solves the streching factor SF with accuracy SFACC.

      IMPLICIT NONE

      INTEGER, PARAMETER :: JMAX = 50
      INTEGER :: J

      REAL :: X1, X2, B, SF, SFACC, F, FMID, DX, XMID

      FMID = SINH(X2)/X2-B
      F    = SINH(X1)/X1-B

      IF(F < 0.0) THEN
        SF = X1
        DX = X2-X1
      ELSE
        SF = X2
        DX = X1-X2
      ENDIF

      DO J=1,JMAX

        DX   = DX*0.5
        XMID = SF+DX
        FMID = SINH(XMID)/XMID-B
        IF(FMID < 0.0) SF = XMID
        IF(ABS(DX) < SFACC .OR. FMID == 0.0) RETURN

      ENDDO

      STOP 'TANH node distribution failed.'

      END SUBROUTINE SFACTOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
