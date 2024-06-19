C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PARTICLE(OSKU,XCG,YCG,ZCG,
     & PSIR,THETAR,PHIR,PSIM,IGR,TDBL,DTB,TRMODE)

                      
C ... Original version (version 1.0)  of the "PARTICLE PROGRAM" 
C ... has been developed in October 1999. This version of the 
C ... "PARTICLE PROGRAM" has been developed by editing 
C ... the version 1.0 in 2007. 
 
      USE CONSTANTS   , ONLY : DEG2RAD 
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : M,S,C1,B,IX,IY,IZ,IXY,IXZ,IYZ,DAMPN,
     & DAMPT,XEA,YEA,ZEA,UPAR,VPAR,WPAR,FIIPAR,THETAPAR,PSIIPAR,PPAR,
     & QPAR,RPAR,TOL,TOL0,H1,HMIN,XEND,TDEL,XT,TSTEER,VLEKO,AKONE,BKONE,
     & UP,VP,WP,PP,QP,RP,FIIP,THETAP,PSIIP,
     & XP,YP,ZP,IPATH,
     &     particleQuaternion0,
     &     particleQuaternion1,
     &     particleQuaternion2,
     &     particleQuaternion3
      USE NS3CO , ONLY : TIMEL, USE_QUATERNIONS_LOGICAL
      
      IMPLICIT NONE

      INTEGER :: ICYCLE,IPRINT,IGR,NGRIFL,TRMODE
      INTEGER :: I, IANS, IO, IARGC

      TYPE (SIXDOF) :: OSKU(*)

      REAL :: XCG,YCG,ZCG,PSIR,THETAR,PHIR,PSIM,TMPXYZ
      REAL :: ROTX,ROTY,ROTZ,VXINI,VYINI,VZINI,TDBL,DTB,TDBLDTB


C ... Select "Free movement" (=6-dof), "Forced movement", "Ejector movement" 
C ... or 2-dof calculation

      IF (TRMODE == 6) THEN   ! Current time step is "Free movement"
         GOTO 100              
      ENDIF                     ! IF (TRMODE == 6) ...

      IF (TRMODE == 61) THEN  ! Current time step is "Forced movement"
         TDBLDTB = TDBL+DTB
         CALL TRAJECTORY_FORMAT(TDBLDTB,DTB,XCG,YCG,ZCG,PSIR,THETAR,
     &        PHIR,PSIM,OSKU,IGR,TIMEL,TRMODE,0.,.FALSE.)
         IF (TRMODE == 6) THEN ! Next time step is "Free movement"
            GOTO 100             ! New values are calculated by RUNGE()
         ELSEIF (TRMODE == 61) THEN ! Next time step is "Forced movement"
            GOTO 200            ! New values comes from TRAJECTORY_FORMAT
         ENDIF                  ! Next time step is "Free movement" ...  
      ENDIF                     ! IF (TRMODE == 61) ...

      IF (TRMODE == 62) THEN
         TDBLDTB = TDBL+DTB               
         CALL TRAJECTORY_FORMAT(TDBLDTB,DTB,XCG,YCG,ZCG,PSIR,THETAR,
     &           PHIR,PSIM,OSKU,IGR,TIMEL,TRMODE,0.,.FALSE.)
         GOTO 100
      ENDIF                     ! IF (TRMODE == 62) ...

      IF (TRMODE == 22) THEN   ! Current time step is 2-dof calculation
         GOTO 100               
      ENDIF                     ! IF (TRMODE == 22) ...

C ... Continue

 100  CONTINUE

C ... Update values for the trajectory simulation

C ... IPATH is the path number, which will be calculated in this program

      IPATH  = IGR

C ... Read the initial conditions for the path

      M     = OSKU(IPATH)%RMASS
      S     = OSKU(IPATH)%AREF 
      C1    = OSKU(IPATH)%CHLREF 
      B     = OSKU(IPATH)%CHLREF
      IX    = OSKU(IPATH)%IX
      IY    = OSKU(IPATH)%IY
      IZ    = OSKU(IPATH)%IZ
      IXY   = OSKU(IPATH)%IXY
      IXZ   = OSKU(IPATH)%IXZ
      IYZ   = OSKU(IPATH)%IYZ

      DAMPN  = OSKU(IPATH)%DAMPN
      DAMPT  = OSKU(IPATH)%DAMPT
      TOL    = OSKU(IPATH)%TOL 
        TOL0 = TOL 
      H1     = OSKU(IPATH)%H1 
      HMIN   = OSKU(IPATH)%HMIN 

      XEND      = TDBL
      TDEL      = DTB
         
      VLEKO     = OSKU(IPATH)%SPEED
      AKONE     = OSKU(IPATH)%ALPHA
      BKONE     = OSKU(IPATH)%BETA

C ... Read the initial conditions for the path (which are given 
C ... in the FINFLO coordinate system) and tranform them to the
C ... flight mechanics coordinate system

      ROTX   = OSKU(IPATH)%ROTX
      ROTY   = OSKU(IPATH)%ROTY
      ROTZ   = OSKU(IPATH)%ROTZ
      VXINI  = OSKU(IPATH)%VX
      VYINI  = OSKU(IPATH)%VY
      VZINI  = OSKU(IPATH)%VZ
      VYINI = VYINI - VLEKO*SIN(AKONE*DEG2RAD)
      VXINI = VXINI - VLEKO*COS(AKONE*DEG2RAD)

      XEA    = -XCG
      TMPXYZ = -YCG
      YEA    = -ZCG
      ZEA    = -TMPXYZ

      FIIPAR   = PHIR
      THETAPAR = THETAR
      PSIIPAR  = PSIR

      IF (USE_QUATERNIONS_LOGICAL) THEN
      CALL EULQUA(PSIIPAR,THETAPAR,FIIPAR,
     &   particleQuaternion0,particleQuaternion1, ! Matti Palin October 2020
     &   particleQuaternion2,particleQuaternion3,
     &   1)
      ENDIF  ! IF (USE_QUATERNIONS_LOGICAL) THEN

      IF (IPATH == 31) THEN
         CALL FINEUL(PSIIPAR,FIIPAR,THETAPAR,
     &        VXINI,VYINI,VZINI,UPAR,VPAR,WPAR,2)
         CALL FINEUL(PSIIPAR,FIIPAR,THETAPAR,
     &        ROTX,ROTY,ROTZ,PPAR,QPAR,RPAR,2)
      ELSEIF (IPATH /= 31) THEN
         CALL FINEUL(PSIIPAR,THETAPAR,FIIPAR,
     &        VXINI,VYINI,VZINI,UPAR,VPAR,WPAR,1)
         CALL FINEUL(PSIIPAR,THETAPAR,FIIPAR,
     &        ROTX,ROTY,ROTZ,PPAR,QPAR,RPAR,1)
      ENDIF

C ... Call "RUNGE" subroutine (in "Free movement" or "Ejector movement"
C ... calculation)

      CALL RUNGE()

C ... Continue and then return

 200  CONTINUE
      
      RETURN
      END SUBROUTINE PARTICLE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     &DERIVS,IND)


C ... Runge-Kutta driver with adaptive stepsize control. Integrate the
C     NVAR starting values YSTART from X1 to X2 with accuracy EPS.
C     H1 should be set as a guessed first stepsize, HMIN as the minimum
C     allowed stepsize (can be zero). On output NOK and NBAD are the
C     number of good and bad (but retried and fixed) steps taken, and
C     YSTART is replaced by values at the end of the integration
C     interval. The subroutine DERIVS (FCN) contains the equations of
C ... motion.

      USE PATH
      IMPLICIT NONE

      INTEGER :: NOK,NBAD,I,NVAR,NSTP,IND,MAXSTP
 
      REAL :: X,X1,H,H1,X2,HDID,XSAV,EPS,HNEXT,HMIN,TINY
      REAL, DIMENSION(NVAR) :: YSTART
      REAL, DIMENSION(16)   :: YSCAL,Y,DYDX

      EXTERNAL DERIVS

      KMAX   = 0
      X      = X1
      H      = SIGN(H1,X2-X1)
      NOK    = 0
      NBAD   = 0
      KOUNT  = 0
      XXP    = 0. ; YYP = 0. ; HDID = 0. ; YSCAL = 0. ; Y=0.; DYDX=0.
      MAXSTP =10000
      TINY   = 1.E-10

      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE

      XSAV = X-DXSAV*2.0
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
12      CONTINUE

        IF(KMAX > 0)THEN
          IF(ABS(X-XSAV) > ABS(DXSAV)) THEN
            IF(KOUNT < KMAX-1)THEN
              KOUNT=KOUNT+1
              XXP(KOUNT)=X
              DO 13 I=1,NVAR
                YYP(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF

        IF((X+H-X2)*(X+H-X1) > 0.0) H=X2-X

        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)       

	WRITE (*,*) "After RKQC: DYDX( 4)= ", DYDX(4)
	WRITE (*,*) "After RKQC: DYDX( 5)= ", DYDX(5)
	WRITE (*,*) "After RKQC: DYDX( 6)= ", DYDX(6)
	WRITE (*,*) "After RKQC: DYDX( 7)= ", DYDX(7)
	WRITE (*,*) "After RKQC: DYDX( 8)= ", DYDX(8)
	WRITE (*,*) "After RKQC: DYDX( 9)= ", DYDX(9)
        
        IF(HDID == H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF

        IF((X-X2)*(X2-X1) >= 0.0)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX /= 0)THEN
            KOUNT=KOUNT+1
            XXP(KOUNT)=X
            DO 15 I=1,NVAR
              YYP(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN

        ENDIF
        IND = 2
        IF(ABS(HNEXT) < HMIN) THEN
          IND = -3
          RETURN
        ENDIF
        H = HNEXT
16    CONTINUE
      STOP 'Too many steps.'

      RETURN
      END SUBROUTINE ODEINT
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)


C ... Fifth-order Runge-Kutta step with monitoring of local truncation
C     error to ensure accuracy and adjust stepsize. Input are the
C     dependent variable vector Y of length N and its derivative DYDX
C     at the starting value of the independent variable X. Also input
C     are the stepsize to be attempted HTRY, the required accuracy EPS,
C     and the vector YSCAL against which the error is scaled. On output,
C     Y and X are replaced by their new values, HDID is the stepsize
C     which was actually accomplished, and HNEXT is the estimate next
C     step. The subroutine DERIVS (FCN) contains the equations of
C ... motion.

      IMPLICIT NONE

      INTEGER :: I,N

      REAL :: PGROW,PSHRNK,XSAV,X,H,HTRY,HH,ERRMAX,EPS,
     &        SAFETY,HDID,ERRCON,HNEXT,FCOR
      REAL, DIMENSION(N)  :: Y,DYDX,YSCAL
      REAL, DIMENSION(16) :: YTEMP,YSAV,DYSAV

      EXTERNAL DERIVS

      PGROW  = -0.20
      PSHRNK = -0.25
      XSAV   = X
      FCOR   = 0.0666666667
      SAFETY = 0.9
      ERRCON = 1.0E-4

      DO 11 I=1,N
        YSAV(I)  = Y(I)
        DYSAV(I) = DYDX(I)
11    CONTINUE
      H  = HTRY
1     HH = 0.5*H

      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)

      X  = XSAV+HH

      CALL DERIVS(X,YTEMP,DYDX)
 
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)

      X  = XSAV+H
      IF(X == XSAV) STOP 'Stepsize not significant in RKQC.'

      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.0

      DO 12 I=1,N
        YTEMP(I) = Y(I)-YTEMP(I)
        ERRMAX   = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE

      ERRMAX = ERRMAX/EPS

      IF(ERRMAX > 1.0) THEN
        H = SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID = H
        IF(ERRMAX > ERRCON)THEN
          HNEXT = SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT = 4.*H
        ENDIF
      ENDIF

      DO 13 I=1,N
        Y(I) = Y(I)+YTEMP(I)*FCOR
13    CONTINUE

      RETURN
      END SUBROUTINE RKQC
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)


C ... Given values for N variables Y and their derivatives DYDX known
C     at X, use the fourth-order Runge-Kutta method to advance the
C     solution over an interval H and return the incremented variables
C     as YOUT. This subroutine calls the subroutine DERIVS (FCN),
C ... which contains the equations of motion.

      IMPLICIT NONE

      INTEGER :: I, N

      REAL, DIMENSION(16) :: Y, DYDX, YOUT, YT, DYT, DYM
      REAL  :: X, H, HH, H6, XH

      EXTERNAL DERIVS

      HH = H*0.5
      H6 = H/6.
      XH = X+HH

      DO 11 I=1,N               ! First step
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE

      CALL DERIVS(XH,YT,DYT)    ! Second step

      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE

      CALL DERIVS(XH,YT,DYM)    ! Third step

      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE

      CALL DERIVS(X+H,YT,DYT)   ! Fourth step

      DO 14 I=1,N           ! Accumulate increments with proper weights
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE

      RETURN
      END SUBROUTINE RK4
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
      SUBROUTINE FCN(X,Y,YPRIME)


C ... Equations of motion

      USE NS3CO     , ONLY   : GX, GY, GZ
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT

      IMPLICIT NONE

      INTEGER :: ICTL

      REAL, DIMENSION(12) :: Y, YPRIME
      REAL :: X
      REAL :: COSALF, SINALF, COSBET, SINBET      
      REAL :: DML, DMD, DMY, DLL, DLD, DLY, DNL, DND, DNY, DM, DL, DN 
      REAL :: GYROX, GYROY, GYROZ, GYROXP, GYROYP, GYROZP
      REAL :: CROL, CM, CN
      REAL :: CXK, CYK, CZK
      REAL :: REY, MACH, MA
      REAL :: Y1F, Y2F, Y3F, Y1, Y2, Y3, V0, V1, X0, Y0, Z0
      REAL :: FX, FY, FZ, MX, MY, MZ, GXX, GYY, GZZ 
      REAL :: EX, EY, EZ, EMX ,EMY, EMZ
      REAL :: TX, TY, TZ, TMX ,TMY, TMZ
      REAL :: HX, HY, HZ, HMX, HMY, HMZ

      INTRINSIC SIN 

C ... Aerodynamic forces and moments, gravity, ejector forces and moments,
C ... thrust forces and moments in the particle body coordinate

C ... Initialize helicopter forces and moments

      HX = 0. ; HY = 0. ; HZ = 0. ; HMX = 0. ; HMY = 0. ; HMZ = 0.

      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%CX,
     &            OSKU(IPATH)%CY,OSKU(IPATH)%CZ,FX,FY,FZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%CMX,
     &            OSKU(IPATH)%CMY,OSKU(IPATH)%CMZ,MX,MY,MZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),GX,GY,GZ,
     &            GXX,GYY,GZZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%FXE,
     &            OSKU(IPATH)%FYE,OSKU(IPATH)%FZE,EX,EY,EZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%MXE,
     &            OSKU(IPATH)%MYE,OSKU(IPATH)%MZE,EMX,EMY,EMZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%FXT,
     &            OSKU(IPATH)%FYT,OSKU(IPATH)%FZT,TX,TY,TZ,1)
      CALL FINEUL(Y(9),Y(8),Y(7),OSKU(IPATH)%MXT,
     &            OSKU(IPATH)%MYT,OSKU(IPATH)%MZT,TMX,TMY,TMZ,1)

C ... Helicopter case is handled in a different way
      IF (TRMODE(IPATH) == 31) THEN
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%CX,
     &        OSKU(IPATH)%CY,OSKU(IPATH)%CZ,FX,FY,FZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%CMX,
     &        OSKU(IPATH)%CMY,OSKU(IPATH)%CMZ,MX,MY,MZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),GX,GY,GZ,
     &        GXX,GYY,GZZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%FXT,
     &        OSKU(IPATH)%FYT,OSKU(IPATH)%FZT,TX,TY,TZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%MXT,
     &        OSKU(IPATH)%MYT,OSKU(IPATH)%MZT,TMX,TMY,TMZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%FXE,
     &        OSKU(IPATH)%FYE,OSKU(IPATH)%FZE,EX,EY,EZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%MXE,
     &        OSKU(IPATH)%MYE,OSKU(IPATH)%MZE,EMX,EMY,EMZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%FXH,
     &        OSKU(IPATH)%FYH,OSKU(IPATH)%FZH,HX,HY,HZ,2)
         CALL FINEUL(Y(9),Y(7),Y(8),OSKU(IPATH)%MXH,
     &        OSKU(IPATH)%MYH,OSKU(IPATH)%MZH,HMX,HMY,HMZ,2)
C ... Gravity forces come via TX,TY,TZ
         GXX = 0.0
         GYY = 0.0
         GZZ = 0.0
C ... Centrifugal terms are edited, because YPIRME terms 
C ... contain part of these centrifugal terms
         EX = EX - (-Y(5)*Y(3)+Y(6)*Y(2))*M
         EY = EY - (-Y(6)*Y(1) + Y(4)*Y(3))*M
         EZ = EZ - (-Y(4)*Y(2)+Y(5)*Y(1))*M
         EMX = EMX - (DL- IXY*(Y(4)*Y(6)-QP)+IXZ*(Y(4)*Y(5)+RP)
     &         + IYZ*(Y(5)**2-Y(6)**2)-Y(5)*Y(6)*(IZ-IY)
     &         - GYROXP+Y(6)*GYROY-Y(5)*GYROZ)
         EMY = EMY - (DM + IXY*(Y(5)*Y(6)+PP)
     &         - IXZ*(Y(4)**2-Y(6)**2)-IYZ*(Y(4)*Y(5)-RP)
     &         - Y(4)*Y(6)*(IX-IZ)-GYROYP-Y(6)*GYROX+Y(4)*GYROZ)
         EMZ = EMZ - (DN - IXY*(Y(5)**2-Y(4)**2)-IXZ*(Y(5)*Y(6)-PP)
     &         + IYZ*(Y(4)*Y(6)+QP)-Y(4)*Y(5)*(IY-IX)
     &         - GYROZP+Y(5)*GYROX-Y(4)*GYROY)
      ENDIF ! IF (TRMODE(IPATH) == 31) THEN
C ... Helicopter definitions end

C ... Set some values to zero
 
      GYROX  = 0.0
      GYROY  = 0.0
      GYROZ  = 0.0
      GYROXP = 0.0
      GYROYP = 0.0
      GYROZP = 0.0

      DML    = 0.0
      DMD    = 0.0
      DMY    = 0.0
      DLL    = 0.0
      DLD    = 0.0
      DLY    = 0.0
      DNL    = 0.0
      DND    = 0.0
      DNY    = 0.0
 
      DM     = DML + DMD + DMY
      DL     = DLL + DLD + DLY
      DN     = DNL + DND + DNY
 
C ... Euler's equations:
C ... =================
 
C ... Forces equations
 
      YPRIME(1) = FX/M + GXX + EX/M + TX/M + HX/M 
     &     - Y(5)*Y(3) + Y(6)*Y(2)

      YPRIME(2) = FY/M + GYY + EY/M + TY/M + HY/M
     &     - Y(6)*Y(1) + Y(4)*Y(3)

      YPRIME(3) = FZ/M + GZZ + EZ/M + TZ/M + HZ/M
     &     - Y(4)*Y(2) + Y(5)*Y(1)

C ... Moment equations
 
      YPRIME(4) = (MX+EMX+TMX+HMX+DL
     &          - IXY*(Y(4)*Y(6)-QP)+IXZ*(Y(4)*Y(5)+RP)
     &          + IYZ*(Y(5)**2-Y(6)**2)-Y(5)*Y(6)*(IZ-IY)
     &          - GYROXP+Y(6)*GYROY-Y(5)*GYROZ)/IX

      YPRIME(5) = (MY+EMY+TMY+HMY+DM+IXY*(Y(5)*Y(6)+PP)
     &          - IXZ*(Y(4)**2-Y(6)**2)-IYZ*(Y(4)*Y(5)-RP)
     &          - Y(4)*Y(6)*(IX-IZ)-GYROYP-Y(6)*GYROX+Y(4)*GYROZ)/IY

      YPRIME(6) = (MZ+EMZ+TMZ+HMZ+DN
     &          - IXY*(Y(5)**2-Y(4)**2)-IXZ*(Y(5)*Y(6)-PP)
     &          + IYZ*(Y(4)*Y(6)+QP)-Y(4)*Y(5)*(IY-IX)
     &          - GYROZP+Y(5)*GYROX-Y(4)*GYROY)/IZ

C ... Kinematic connections

      IF(ABS(ABS(Y(8))-1.570796) <= 0.0017) THEN
         ICTL = 1
         YPRIME(7) = 0.0
         YPRIME(8) = Y(5)*COS(Y(7))-Y(6)*SIN(Y(7))
         YPRIME(9) = 0.0
      ELSE
         ICTL = 0
         YPRIME(7) = Y(4)+Y(5)*SIN(Y(7))*TAN(Y(8))
     &             + Y(6)*COS(Y(7))*TAN(Y(8))
         YPRIME(8) = Y(5)*COS(Y(7))-Y(6)*SIN(Y(7))
         YPRIME(9) = Y(5)*SIN(Y(7))/COS(Y(8))
     &             + Y(6)*COS(Y(7))/COS(Y(8))
      ENDIF
 
C ... Location equations

      YPRIME(10) = Y(1)*COS(Y(8))*COS(Y(9))
     &           + Y(2)*(SIN(Y(7))*SIN(Y(8))*COS(Y(9))
     &           - COS(Y(7))*SIN(Y(9)))
     &           + Y(3)*(COS(Y(7))*SIN(Y(8))*COS(Y(9))
     &           + SIN(Y(7))*SIN(Y(9)))
 
      YPRIME(11) = Y(1)*COS(Y(8))*SIN(Y(9))
     &           + Y(2)*(SIN(Y(7))*SIN(Y(8))*SIN(Y(9))
     &           + COS(Y(7))*COS(Y(9)))
     &           + Y(3)*(COS(Y(7))*SIN(Y(8))*SIN(Y(9))
     &           - SIN(Y(7))*COS(Y(9)))
 
      YPRIME(12) = Y(1)*SIN(Y(8))-Y(2)*SIN(Y(7))*COS(Y(8))
     &           - Y(3)*COS(Y(7))*COS(Y(8))

C ... Define the time derivatives

      UP     = YPRIME(1)
      VP     = YPRIME(2)
      WP     = YPRIME(3)
      PP     = YPRIME(4)
      QP     = YPRIME(5)
      RP     = YPRIME(6)
      FIIP   = YPRIME(7)
      THETAP = YPRIME(8)
      PSIIP  = YPRIME(9)
      XP     = YPRIME(10)
      YP     = YPRIME(11)
      ZP     = YPRIME(12)

      IF(ABS(YPRIME(1)) <= 1E-15) UP = 0.0
      IF(ABS(YPRIME(2)) <= 1E-15) VP = 0.0
      IF(ABS(YPRIME(3)) <= 1E-15) WP = 0.0
      IF(ABS(YPRIME(4)) <= 1E-15) PP = 0.0
      IF(ABS(YPRIME(5)) <= 1E-15) QP = 0.0
      IF(ABS(YPRIME(6)) <= 1E-15) RP = 0.0

      RETURN
      END SUBROUTINE FCN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RUNGE

C ... Integration of the equations of motion

      USE CONSTANTS , ONLY : DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT
      USE NS3CO , ONLY : ROTORW, USE_QUATERNIONS_LOGICAL
 
      IMPLICIT NONE

      INTEGER, PARAMETER :: NPMAX = 50000
      INTEGER :: N, IER, NPRINT
      INTEGER :: NOK, NBAD, IND, IND0, NP

      REAL, DIMENSION(16) :: Y
      REAL :: Y1, Y2, Y3, Y4, Y5, Y6
      REAL :: RBLADE

      EXTERNAL FCN, FCNQuaternion

C ... Initialize something

      NPRINT = 1
      NP     = 0
      N      = 16
      IND    = 2
      IND0   = IND
      XT     = 0.0

C ... Initialize the time derivatives (only PP, QP, RP are needed)

      UP     = 0.0
      VP     = 0.0 
      WP     = 0.0 
      PP     = 0.0
      QP     = 0.0
      RP     = 0.0 
      FIIP   = 0.0 
      THETAP = 0.0 
      PSIIP  = 0.0 
      XP     = 0.0 
      YP     = 0.0 
      ZP     = 0.0

C ... Initialize the variables to be integrated in time
      
      Y(1)   = UPAR
      Y(2)   = VPAR
      Y(3)   = WPAR
      Y(4)   = PPAR
      Y(5)   = QPAR
      Y(6)   = RPAR
      Y(7)   = FIIPAR
      Y(8)   = THETAPAR
      Y(9)   = PSIIPAR
      Y(10)  = XEA
      Y(11)  = YEA
      Y(12)  = ZEA
      Y(13)  = particleQuaternion0  ! The quaternion components
      Y(14)  = particleQuaternion1  ! Added November 2020 by Matti Palin.
      Y(15)  = particleQuaternion2
      Y(16)  = particleQuaternion3

  200 CONTINUE

   90 IF(NPRINT == 0) GO TO 100

C ... VX,VY,VZ,ROTX,ROTY,ROTZ at t = XEND

      IF (TRMODE(IPATH) == 31) THEN ! helicopter case
         CALL EULFIN(Y(9),Y(7),Y(8),Y(1),Y(2),Y(3),Y1,Y2,Y3,2)
         CALL EULFIN(Y(9),Y(7),Y(8),Y(4),Y(5),Y(6),Y4,Y5,Y6,2)
      ELSEIF (TRMODE(IPATH) /= 31) THEN ! normal case
         CALL EULFIN(Y(9),Y(8),Y(7),Y(1),Y(2),Y(3),Y1,Y2,Y3,1)
         CALL EULFIN(Y(9),Y(8),Y(7),Y(4),Y(5),Y(6),Y4,Y5,Y6,1)
      ENDIF

  100 CONTINUE
      
C ... Next time step
      TNEXT = XEND + TDEL
      nok = 0 ; nbad = 0

      
  300 CONTINUE

C ... Integration of the equations of motion using Runge-Kutta method

      XT = TNEXT - TDEL   
  
      IF(USE_QUATERNIONS_LOGICAL) THEN
           CALL ODEINT(Y,N,XT,TNEXT,TOL,H1,HMIN,NOK,NBAD,
     & FCNQuaternion,IND)

           particleQuaternion0 = Y(13)	! The quaternion components
           particleQuaternion1 = Y(14)	! Added January 4, 2019 by Matti Palin.
           particleQuaternion2 = Y(15)
           particleQuaternion3 = Y(16)

         CALL QUAEUL(			! Convert quaternions back to Euler angles
     &    Y(9),Y(8), Y(7),
     &    particleQuaternion0,
     &    particleQuaternion1,
     &    particleQuaternion2,
     &    particleQuaternion3, 1)
      ELSE
           CALL ODEINT(Y,N,XT,TNEXT,TOL,H1,HMIN,NOK,NBAD,FCN,IND)
      ENDIF

C ... Happy end? In unhappy case increase the tolerance and try again. 

      IF(IND /= -3) GO TO 250

C ... Increase the tolerance

      TOL = TOL*10.
      IF(TOL >= 0.1) GO TO 350
      IER = 0
      IND = IND0
      WRITE(6,1500) TOL
 1500 FORMAT('  !!! Tolerance increased, new value',1X,E10.1)

C ... Try again
      GO TO 300

  350 CONTINUE

C ... Cannot continue without exceeding maximum tolerance 
      WRITE(6,1600)
 1600 FORMAT('  !!! Maximum tolerance 0.1 exeeded. STOP.')

      GO TO 50

  250 CONTINUE

C ... Restore the original value of tolerance
      TOL = TOL0

C ... Write information to the screen if time derivatives of the 
C ... yaw angle and the roll angle is set to zero

      IF(ABS(ABS(Y(8))-1.570796) <= 0.0017) WRITE(6,1700) XEND
 1700 FORMAT(/,' !!! Time derivatives of the yaw angle and the roll',/,
     &         '     angle is set to zero at the timelevel t=',F7.3,'s')

      IF(IND > 0) GO TO 400

      WRITE(6,1800)
 1800 FORMAT('  !!! Execution of ODEINT failed.')

C ... Exit
  400 IF(IND < 0) GO TO 50

C ... VX,VY,VZ,ROTX,ROTY,ROTZ at t = TNEXT

      IF (TRMODE(IPATH) == 31) THEN ! helicopter case
         CALL EULFIN(Y(9),Y(7),Y(8),Y(1),Y(2),Y(3),Y1,Y2,Y3,2)
         CALL EULFIN(Y(9),Y(7),Y(8),Y(4),Y(5),Y(6),Y4,Y5,Y6,2)
      ELSEIF (TRMODE(IPATH) /= 31) THEN ! normal case
         CALL EULFIN(Y(9),Y(8),Y(7),Y(1),Y(2),Y(3),Y1,Y2,Y3,1)
         CALL EULFIN(Y(9),Y(8),Y(7),Y(4),Y(5),Y(6),Y4,Y5,Y6,1)
      ENDIF
  

C ... Update FLIGHT and SIXDOF modules

      IF (TRMODE(IPATH) == 6 .OR. TRMODE(IPATH) == 61 .OR. 
     &     TRMODE(IPATH) == 62) THEN ! 6-dof

      XCG(IPATH) = -Y(10)+TDEL*(VLEKO*COS(AKONE*DEG2RAD))
      YCG(IPATH) = Y(12)+TDEL*(VLEKO*SIN(AKONE*DEG2RAD))
      ZCG(IPATH) = -Y(11)

      OSKU(IPATH)%VX = Y1+(VLEKO*COS(AKONE*DEG2RAD))
      OSKU(IPATH)%VY = Y2+(VLEKO*SIN(AKONE*DEG2RAD))
      OSKU(IPATH)%VZ = Y3

      PHIR(IPATH)   = Y(7)
      THETAR(IPATH) = Y(8)
      PSIR(IPATH)   = Y(9)

      OSKU(IPATH)%ROTX = Y4
      OSKU(IPATH)%ROTY = Y5
      OSKU(IPATH)%ROTZ = Y6

      ELSEIF (TRMODE(IPATH) == 22) THEN ! 2-dof

      XCG(IPATH) = -Y(10)
      YCG(IPATH) = Y(12)
      ZCG(IPATH) = -Y(11)

      OSKU(IPATH)%VX = Y1
      OSKU(IPATH)%VY = Y2
      OSKU(IPATH)%VZ = Y3

      PHIR(IPATH)   = Y(7)
      THETAR(IPATH) = Y(8)
      PSIR(IPATH)   = Y(9)
       
      OSKU(IPATH)%ROTX = Y4
      OSKU(IPATH)%ROTY = Y5
      OSKU(IPATH)%ROTZ = Y6

      TRIMA(IPATH) = PSIR(IPATH)
      SINK(IPATH)  = DRAUGHT(IPATH)+(ZCGIS(IPATH)-ZCG(IPATH))

      ELSEIF (TRMODE(IPATH) == 31) THEN ! Helicopter rotor

      RBLADE=SQRT(XCG(IPATH)*XCG(IPATH)
     &           +ZCG(IPATH)*ZCG(IPATH))
      PSIM(IPATH)= PSIM(IPATH)-ROTORW*TDEL ! updates shaft angle

      XCG(IPATH) = RBLADE*SIN(PSIM(IPATH)) ! updates translations
      YCG(IPATH) = YCG(IPATH)
      ZCG(IPATH) = -RBLADE*COS(PSIM(IPATH))

      OSKU(IPATH)%VX = -ROTORW*RBLADE*COS(PSIM(IPATH)) ! updates speeds
      OSKU(IPATH)%VY = OSKU(IPATH)%VY
      OSKU(IPATH)%VZ = -ROTORW*RBLADE*SIN(PSIM(IPATH))


c      XCG(IPATH) = -Y(10) ! updates translations
c      YCG(IPATH) = Y(12)
c      ZCG(IPATH) = -Y(11)

c      OSKU(IPATH)%VX = Y1 ! updates speeds
c      OSKU(IPATH)%VY = Y2 ! naa ei toimi en kasita miksi
c      OSKU(IPATH)%VZ = Y3

      PHIR(IPATH)      = Y(7) ! updates  attitude angles
      THETAR(IPATH)    = Y(8)
      PSIR(IPATH)      = Y(9)

      OSKU(IPATH)%ROTX = Y4 ! updates rotational speeds
      OSKU(IPATH)%ROTY = Y5
      OSKU(IPATH)%ROTZ = Y6


CC      ROTPP(IPATH) = YPRIME(4) ! updates rotational accelerations
CC      ROTQP(IPATH) = YPRIME(5)
CC      ROTRP(IPATH) = YPRIME(6)

      ENDIF ! (TRMODE(IPATH) >= 31 .AND. <= 39) THEN

   50 CONTINUE

      RETURN
      END SUBROUTINE RUNGE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FINEUL(PSI,THETA,PHI,XO,YO,ZO,XONEW,YONEW,ZONEW,ICASE)


C ... FINFLO (PLOT3D) coordinates -> Local body coordinates. 

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE
      REAL :: PSI, THETA, PHI
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: S1, S2, S3, C1, C2, C3
      REAL :: XO, YO, ZO, XN, YN, ZN
      REAL :: XFM, YFM, ZFM, XONEW, YONEW, ZONEW     


C ... Transformation matrix

      S1 = SIN(PHI)
      S2 = SIN(THETA)
      S3 = SIN(PSI)
      C1 = COS(PHI)
      C2 = COS(THETA)
      C3 = COS(PSI)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A21 =  C2*S3
      A31 = -S2
      A12 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A32 =  C2*S1
      A13 =  C1*C3*S2 + S1*S3
      A23 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A21 =  C3*S1*S2 + C1*S3
      A31 = -C2*S1
      A12 = -C2*S3 
      A22 =  C2*C3
      A32 =  S2
      A13 =  C3*S1 + C1*S2*S3
      A23 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A21 =  S2
      A31 = -C2*S3
      A12 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A32 =  C3*S1 + C1*S2*S3
      A13 =  C3*S1*S2 + C1*S3
      A23 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A21 =  C1*C3*S2 + S1*S3
      A31 =  C1*S2*S3 - C3*S1
      A12 = -S2
      A22 =  C2*C3
      A32 =  C2*S3
      A13 =  C2*S1
      A23 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A21 =  C3*S1 + C1*S2*S3
      A31 =  S1*S3 - C1*C3*S2
      A12 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A32 =  C3*S1*S2 + C1*S3
      A13 =  S2
      A23 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A21 =  C2*S1
      A31 =  C3*S1*S2 - C1*S3
      A12 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A32 =  C1*C3*S2 + S1*S3
      A13 =  C2*S3
      A23 = -S2
      A33 =  C2*C3

      END SELECT


C ... From FINFLO coordinate system to the flight mechanics 
C ... coordinate system.

      XFM = -XO
      YFM = -ZO 
      ZFM = -YO

C ... Rotations according to the Euler angles.

      XN = A11*XFM + A21*YFM + A31*ZFM
      YN = A12*XFM + A22*YFM + A32*ZFM
      ZN = A13*XFM + A23*YFM + A33*ZFM

      XONEW = XN
      YONEW = YN
      ZONEW = ZN

      RETURN
      END SUBROUTINE FINEUL
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE EULFIN(PSI,THETA,PHI,XO,YO,ZO,XONEW,YONEW,ZONEW,ICASE)


C ... Local body coordinates -> FINFLO (PLOT3D) coordinates.

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE

      REAL :: PSI, THETA, PHI
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XO, YO, ZO, XN, YN, ZN
      REAL :: XONEW, YONEW, ZONEW     
      REAL :: S1, S2, S3, C1, C2, C3

      
C ... Transformation matrix

      S1 = SIN(PHI)
      S2 = SIN(THETA)
      S3 = SIN(PSI)
      C1 = COS(PHI)
      C2 = COS(THETA)
      C3 = COS(PSI)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A12 =  C2*S3
      A13 = -S2
      A21 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A23 =  C2*S1
      A31 =  C1*C3*S2 + S1*S3
      A32 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A12 =  C3*S1*S2 + C1*S3
      A13 = -C2*S1
      A21 = -C2*S3 
      A22 =  C2*C3
      A23 =  S2
      A31 =  C3*S1 + C1*S2*S3
      A32 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A12 =  S2
      A13 = -C2*S3
      A21 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A23 =  C3*S1 + C1*S2*S3
      A31 =  C3*S1*S2 + C1*S3
      A32 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A12 =  C1*C3*S2 + S1*S3
      A13 =  C1*S2*S3 - C3*S1
      A21 = -S2
      A22 =  C2*C3
      A23 =  C2*S3
      A31 =  C2*S1
      A32 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A12 =  C3*S1 + C1*S2*S3
      A13 =  S1*S3 - C1*C3*S2
      A21 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A23 =  C3*S1*S2 + C1*S3
      A31 =  S2
      A32 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A12 =  C2*S1
      A13 =  C3*S1*S2 - C1*S3
      A21 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A23 =  C1*C3*S2 + S1*S3
      A31 =  C2*S3
      A32 = -S2
      A33 =  C2*C3

      END SELECT


C ... Rotations according to the Euler angles.

      XN = A11*XO + A21*YO + A31*ZO
      YN = A12*XO + A22*YO + A32*ZO
      ZN = A13*XO + A23*YO + A33*ZO


C ... Global directions back to the FINFLO directions.

      XONEW = -XN
      YONEW = -ZN
      ZONEW = -YN

      RETURN
      END SUBROUTINE EULFIN
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE PARTICLE_FORCES(CX,CY,CZ,CMX,CMY,CMZ,ICON,
     &     XCGP,YCGP,ZCGP,CBTWF,IBC1,IBC2,IBC3,IGR,IGRID)

C     This subroutine has been developed from "AUXFOR" subroutine.
C     This subroutine calculates the forces excerted on various
C     sub-systems based on patch-oriented data.
C
C     INPUT:
C        CXB,CYB,CZB,CMXB,CMYB,CMZB, = Forces and moments excerted
C        ALPHA                       = Angle of attack
C        BETA                        = Sideslip angle
C        ICON                        = BC patch data
C        NBCS                        = Number of BCs
C        APATCH                      = Patch area
C                                      Flying Object 
C        XCGP,YCGP,ZCGP              = Flying Object location in the 
C                                      FINFLO coordinate system
C        XMOM,YMOM,ZMOM              = Moment refrence point
C        AREF,CHLREF,REFPRE          = Reference area,length and pressure
C
C
C        CBTWF                       = Integer, used when coefficients
C                                      are changed to forces and moments 
C                                      CBTW=1 -> coefficients->forces+moments
C        IBC1,IBC2,IBC3              = Boundrary condition type
C        IGR                         = Flying object number
C        IGRID                       = IGR value for every block
C
C
C     OUTPUT:
C        CX,CY,CZ,CMX,CMY,CMZ        = Forces and moments in the 
C                                      FINFLO coordinate system 

      USE NS3CO, ONLY : IC9,ALPHA,BETA,PARALLEL,AREF,CHLREF,REFPRE,
     &     XMOM,YMOM,ZMOM,PARALLEL
      USE MPI

      USE MAIN_ARRAYS , ONLY : CXB,CYB,CZB,CMXB,CMYB,CMZB,APATCH
      USE INTEGERS, ONLY : NBCS,IPRO
      
      INTEGER :: ICON(IC9,*)

      INTEGER :: ITIMES,CBTWF,IBC1,IBC2,IBC3,IGR,IGRID(*)
          
      REAL :: CXS,CYS,CZS,CMXS,CMYS,CMZS,CLS,CDS,CSS,AAS
      REAL :: CX,CY,CZ,CMX,CMY,CMZ
      
      REAL :: XCGP,YCGP,ZCGP

      CHARACTER(LEN=1) :: GNAME

      LOGICAL :: GROUP, GROUPX

c      DO J=1,52    ! Scan through all alphabetics

         
         GROUP = .FALSE.

         CXS   = 0.0
         CYS   = 0.0
         CZS   = 0.0
         CMXS  = 0.0
         CMYS  = 0.0
         CMZS  = 0.0
         CLS   = 0.0
         CDS   = 0.0
         CSS   = 0.0
         AAS   = 0.0


         DO I=1,NBCS  ! Scan through all patches
            
            IF (ICON(1,I) == IBC1 .OR. ICON(1,I) == IBC2) THEN
                                          !IBC1,IBC2 surface types
               IF (IBC1 == 1 .OR. IBC2 == 1) THEN ! ACT disk
                  
                  IF (ICON(20,I) == IBC3 .AND. ! ACT disk => ICON(20,I)=8  
     &                 IGRID(ICON(24,I)) == IGR) THEN ! IGRID number
                     GROUP = .TRUE.              
                     CXS   =  CXS  +  CXB(I)
                     CYS   =  CYS  +  CYB(I)
                     CZS   =  CZS  +  CZB(I)
                     CMXS  =  CMXS +  CMXB(I)
                     CMYS  =  CMYS +  CMYB(I)
                     CMZS  =  CMZS +  CMZB(I)
                     AAS   =  AAS  +  APATCH(I)
                  ENDIF

               ELSE ! NO ACT disk
                         
                  IF (IGRID(ICON(24,I)) == IGR) THEN ! IGRID number
                     GROUP = .TRUE.
                     CXS   =  CXS  +  CXB(I)
                     CYS   =  CYS  +  CYB(I)
                     CZS   =  CZS  +  CZB(I)
                     CMXS  =  CMXS +  CMXB(I)
                     CMYS  =  CMYS +  CMYB(I)
                     CMZS  =  CMZS +  CMZB(I)
                     AAS   =  AAS  +  APATCH(I)
                  ENDIF

               ENDIF
                  
            ENDIF

         ENDDO

         CLS = (-CXS*SIN(ALPHA) + CYS*COS(ALPHA))
         CDS = CXS*COS(ALPHA)*COS(BETA) + CYS*SIN(ALPHA)*COS(BETA)
     &        +CZS*SIN(BETA)
         CSS = CXS*COS(ALPHA)*SIN(BETA) + CYS*SIN(ALPHA)*SIN(BETA)
     &        +CZS*COS(BETA)

C ... Final forces and moments are formed by using 2 If Statements 
         
         IF (GROUP .AND. .NOT. PARALLEL) THEN ! If Statement #1, not MPI
               CX  = CXS
               CY  = CYS
               CZ  = CZS
               CMX = CMXS
               CMY = CMYS
               CMZ = CMZS

         ENDIF
         
         IF (PARALLEL) THEN         ! If Statement #2, yes MPI
            CALL MPI_REDUCE(GROUP,GROUPX,1,MPI_LOGICAL,MPI_LOR,0,
     &           MPI_COMM_WORLD,IERR)
            IF (IPRO == 1) GROUP = GROUPX
            CALL MPI_BCAST(GROUP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

            IF (GROUP) THEN
               CALL MPI_REDUCE(CLS, CLSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CDS, CDSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CSS, CSSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CXS, CXSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CYS, CYSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CZS, CZSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(DXS, DXSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(DYS, DYSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(DZS, DZSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(QTS, QTSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(QWS, QWSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(QHS, QHSM, 1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CMXS,CMXSM,1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CMYS,CMYSM,1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(CMZS,CMZSM,1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(TOS,TOSM,  1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(AAS,AASM,  1,MPI_REAL8,MPI_SUM,0,
     &              MPI_COMM_WORLD,IERR)
            ENDIF

            IF (GROUP .AND. IPRO == 1) THEN 
               CX  = CXSM
               CY  = CYSM
               CZ  = CZSM
               CMX = CMXSM
               CMY = CMYSM
               CMZ = CMZSM
            ENDIF               ! IF(GROUP .AND. IPRO == 1) THEN 

         ENDIF
      
c      ENDDO ! j = 1,52

C ... Change moments from the moment refrence point 
C ... to the CG of the Flying Object

      IF (GROUP .AND. IPRO == 1) THEN
         IF (CBTWF == 1) THEN ! Coefficients -> Forces and Moments
            CMX = CMX + CZ*(YMOM-YCGP)/CHLREF - CY*(ZMOM-ZCGP)/CHLREF
            CMY = CMY - CZ*(XMOM-XCGP)/CHLREF + CX*(ZMOM-ZCGP)/CHLREF
            CMZ = CMZ + CY*(XMOM-XCGP)/CHLREF - CX*(YMOM-YCGP)/CHLREF

            CX  = CX*REFPRE*AREF
            CY  = CY*REFPRE*AREF
            CZ  = CZ*REFPRE*AREF
            CMX = CMX*REFPRE*AREF*CHLREF
            CMY = CMY*REFPRE*AREF*CHLREF
            CMZ = CMZ*REFPRE*AREF*CHLREF
         ELSEIF (CBTWF == 2) THEN ! Forces and Moments -> Forces and Moments
            CMX = CMX + CZ*(YMOM-YCGP) - CY*(ZMOM-ZCGP)
            CMY = CMY - CZ*(XMOM-XCGP) + CX*(ZMOM-ZCGP)
            CMZ = CMZ + CY*(XMOM-XCGP) - CX*(YMOM-YCGP)

            CX  = CX
            CY  = CY         
            CZ  = CZ
            CMX = CMX
            CMY = CMY
            CMZ = CMZ
         ELSEIF (CBTWF == 3) THEN
            CONTINUE
         ENDIF
      ENDIF                     ! IF(GROUP .AND. IPRO == 1) THEN

      RETURN
      END SUBROUTINE PARTICLE_FORCES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLAPPOS(XCO,YCO,ZCO,XGRIG,YGRIG,
     &                   ZGRIG,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,N,TDBL,DTB,NGPTS)

      USE CONSTANTS, ONLY : DEG2RAD

      IMPLICIT NONE

      INTEGER :: IBLOCK,IMAX,JMAX,KMAX,IN,JN,KN,N,IG1,IHEADER,I,NGPTS

      REAL, DIMENSION(*) :: XCO,YCO,ZCO,XGRIG,YGRIG,ZGRIG
      REAL :: TDBL,DTB,TIME
      REAL :: OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,
     &        ROTANG

      CHARACTER(LEN=235) :: HEADER

      LOGICAL :: THERE

C ... Update an original grid to a new grid
      DO I = 1,NGPTS
         XCO(I) = XGRIG(I)
         YCO(I) = YGRIG(I)
         ZCO(I) = ZGRIG(I)
      ENDDO

C ... Check if "FLAPPOS.dat" file exists and then open channel 500 or return
      INQUIRE(FILE='FLAPPOS.dat',EXIST=THERE)
      IF(THERE) THEN          
         OPEN(500, FILE='FLAPPOS.dat',
     &         STATUS='UNKNOWN', FORM='FORMATTED') 
      ELSE 
         RETURN
      ENDIF

C ... Skip headerlines
      DO IHEADER=1,3
         READ(500,'(A235)') HEADER
      ENDDO

C ... Search a correct row (=correct timelevel and correct block)
 100  CONTINUE
      READ(500,*,END=300) TIME,IBLOCK,OMEGAX,OMEGAY,OMEGAZ,
     &     CENAX,CENAY,CENAZ,ROTANG
C ... Rotation angle degs. -> rads        
         ROTANG = ROTANG*DEG2RAD
C ... Flap movement for a block, which is defined in "FLAPPOS.dat" file
      IF ((ABS(TIME-TDBL))/DTB <= 0.33333 .AND. N == IBLOCK) THEN
         CALL ROTGRI(XCO,YCO,ZCO,XGRIG,
     &        YGRIG,ZGRIG,IMAX,JMAX,KMAX,
     &        ROTANG,IN,JN,KN,OMEGAX,OMEGAY,OMEGAZ,
     &        CENAX,CENAY,CENAZ)
      ELSE 
         GOTO 100               
      ENDIF

C ... Continue, format of "FLAPPOS.dat" file, close cahnnel 500
 300  CONTINUE
 400  FORMAT(F6.4,I7,7F8.3)
      CLOSE(500)
      
      RETURN
      END SUBROUTINE FLAPPOS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_WRITE(TDBL,DTB,XCGI,YCGI,ZCGI,
     &           XCG,YCG,ZCG,PSIR,THETAR,PHIR,PSIM,
     &           OSKU,IGR,TIMEL,TRMODE,ROTA1,ROTB1,CONEA)


      USE CONSTANTS   , ONLY : RAD2DEG
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE NS3CO       , ONLY : ICYCLE,ICYOLD,IOLD
      USE BLADE_VARIABLES , ONLY : DOTZE,DOTBETA,DOTTHETA,SHAFT,
     &                             THETACOL,THCYCLON,THCYCLAT,QTH,
     &                             QTH_A,QBE_A
      USE FLIGHT      , ONLY : FLYOBJ,ADV,FDSP,UTSP,THRUST,TORQUE,IFA,
     &     RTMSP,FDSP,FXTSP,QFACT,XFAKEPNEW,YFAKEPNEW,ZFAKEPNEW

      INTEGER :: IGR,IHEADER,I,ITRAJEC,TRMODE,TRMODE2

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      REAL :: PSIR,THETAR,PHIR,PSIM,XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &        P,Q,R,ROTA1,ROTB1,CONEA,TX,TY,TZ,TDBL,DTB,TIME

      CHARACTER(LEN=188) :: TRAJECF
      CHARACTER(LEN=188), DIMENSION(10000) :: TRAJECTF
      CHARACTER(LEN=437) :: TRAJECM
      CHARACTER(LEN=437), DIMENSION(10000) :: TRAJECTM
      CHARACTER(LEN=437) :: HEADER
      CHARACTER(LEN=3)   :: IGRNUM

      LOGICAL :: TIMEL,THERE

C ... Check if "TRAJECTORY_IGR" files exist, open channels 400+IGls R and 500+IGR
      CALL NUMCH3(IGRNUM,IGR)

      INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.dat',EXIST=THERE)
      IF(THERE) THEN          
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &        STATUS='UNKNOWN', FORM='FORMATTED')
      ELSE 
         WRITE(*,*)' '
         WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.dat',
     &                ' was not found'
         WRITE(*,*)' in TRAJECTORY_WRITE Subroutine.' 
         WRITE(*,*)' '
         WRITE(*,*)' Exiting...'
         WRITE(*,*)' '
         STOP
      ENDIF

C ... Binary file, channel 500+IGR flying object and ships only
      INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',EXIST=THERE)
      IF(THERE) THEN          
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &        STATUS='UNKNOWN', FORM='UNFORMATTED')
      ELSE 
         WRITE(*,*)' '
         WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.BIN',
     &        ' was not found'
         WRITE(*,*)' in TRAJECTORY_WRITE Subroutine.' 
         WRITE(*,*)' '
         WRITE(*,*)' Exiting...'
         WRITE(*,*)' '
         STOP
      ENDIF

      
C ... Skip headerlines
      DO IHEADER=1,16
         READ(400+IGR,'(A307)') HEADER
      ENDDO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... NON-TIME ACCURARTE CALCULATION
      IF (.NOT. TIMEL) THEN     
C ... TRMODE is set to "61" and ejector forces and moments to "0" if NOT
C ... ship/helicopter rotor case
c         IF (TRMODE /= 10.AND.TRMODE < 21 .OR. TRMODE  > 49) THEN
         IF (FLYOBJ) THEN
            TRMODE = 61
            OSKU(IGR)%FXE = 0
            OSKU(IGR)%FYE = 0
            OSKU(IGR)%FZE = 0
            OSKU(IGR)%MXE = 0
            OSKU(IGR)%MYE = 0
            OSKU(IGR)%MZE = 0
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Write values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2000)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
C ... Write values to TRAJECTORY_IGR.BIN file
            WRITE(500+IGR)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
            
C ... Ship case
         ELSEIF (TRMODE >= 21 .AND. TRMODE <= 29) THEN
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Search a correct row (TRAJECTORY_IGR.dat file)
 40         CONTINUE
            READ(400+IGR,'(1X,I10)',END=50) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 50
            ELSE
               GOTO 40
            ENDIF
 50         CONTINUE
            BACKSPACE(400+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2500)ICYCLE,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &           XFAKEPNEW(IGR),YFAKEPNEW(IGR),ZFAKEPNEW(IGR)
C ... Search a correct row (TRAJECTORY_IGR.BIN file)
 52         CONTINUE
            READ(500+IGR,END=54) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 54
            ELSE
               GOTO 52
            ENDIF
 54         CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            IF (IOLD < 0) THEN
               CLOSE(500+IGR)
               OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &              STATUS='UNKNOWN', FORM='UNFORMATTED')
            ENDIF
            WRITE(500+IGR)ICYCLE,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &           XFAKEPNEW(IGR),YFAKEPNEW(IGR),ZFAKEPNEW(IGR)

C ... Helicopter rotor case
         ELSEIF (TRMODE >= 31 .AND. TRMODE <= 39) THEN

C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,PHIR,THETAR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,2)
C ... Search a correct row (TRAJECTORY_IGR.dat)
 70        CONTINUE
            READ(400+IGR,'(1X,I10)',END=80) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 80
            ELSE
               GOTO 70
            ENDIF
 80         CONTINUE
            BACKSPACE(400+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,3000)ICYCLE,PSIM*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &           OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &           OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
C ... Search a correct row (TRAJECTORY_IGR.BIN)
 82         CONTINUE
            READ(500+IGR,END=83) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 83
            ELSE
               GOTO 82
            ENDIF
 83         CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to TRAJECTORY_IGR.BIN file
            WRITE(500+IGR)ICYCLE,PSIM*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &           OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &           OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
           
C ... Actuator disk case
         ELSEIF (TRMODE >= 40 .AND. TRMODE <= 49) THEN
C ... Search a correct row TRAJECTORY.dat file
 85         CONTINUE
            READ(400+IGR,'(1X,I10)',END=88) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 88
            ELSE
               GOTO 85
            ENDIF
 88         CONTINUE
            BACKSPACE(400+IGR)
C ... Search a correct row TRAJECTORY.BIN file
 89         CONTINUE
            READ(500+IGR,END=91) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 91
            ELSE
               GOTO 89
            ENDIF
 91         CONTINUE
            BACKSPACE(500+IGR)    
C ... Write calculated values to TRAJECTORY_IGR.dat and TRAJECTORY.BIN file
            IF(TRMODE == 40) THEN ! Patria's helicopt
               IF(IFA(IGR) == 12 .OR. IFA(IGR) == 13)THEN
                  IF ((ICYCLE-ICYOLD) < 1) THEN
                     THETACOL(IGR)=THETACOL(IGR)+0.0
                     THCYCLON(IGR)=THCYCLON(IGR)+0.0
                     THCYCLAT(IGR)=THCYCLAT(IGR)+0.0
                  ELSE
                     THETACOL(IGR)=THETACOL(IGR)+QTH(IGR)
                     THCYCLON(IGR)=THCYCLON(IGR)+QTH_A(IGR)
                     THCYCLAT(IGR)=THCYCLAT(IGR)+QBE_A(IGR)
                  ENDIF
               ENDIF
               CALL EULFIN(0.0,ROTB1,ROTA1,
     &              0.0,0.0,-SHAFT(IGR)/ABS(SHAFT(IGR)),
     &              TX,TY,TZ,2)
               TORQUE(IGR) = (OSKU(IGR)%MXT*TX+OSKU(IGR)%MYT*TY
     &              +OSKU(IGR)%MZT*TZ)
               FXTSP(IGR) = (OSKU(IGR)%FXT*TX+OSKU(IGR)%FYT*TY
     &              +OSKU(IGR)%FZT*TZ)
               WRITE(400+IGR,3100)ICYCLE,XCG,YCG,ZCG,
     &              ROTA1*RAD2DEG,ROTB1*RAD2DEG,
     &              CONEA*RAD2DEG,THETACOL(IGR)*RAD2DEG,
     &              THCYCLON(IGR)*RAD2DEG,THCYCLAT(IGR)*RAD2DEG,
     &              OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &              OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &              THRUST(IGR),FXTSP(IGR),TORQUE(IGR),QFACT(IGR),
     &              OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT
               WRITE(500+IGR)ICYCLE,XCG,YCG,ZCG,
     &              ROTA1*RAD2DEG,ROTB1*RAD2DEG,
     &              CONEA*RAD2DEG,THETACOL(IGR)*RAD2DEG,
     &              THCYCLON(IGR)*RAD2DEG,THCYCLAT(IGR)*RAD2DEG,
     &              OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &              OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &              THRUST(IGR),FXTSP(IGR),TORQUE(IGR),QFACT(IGR),
     &              OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT
               RTMSP(IGR) = OSKU(IGR)%MZT/FXTSP(IGR)
               FDSP(IGR)  = OSKU(IGR)%MXT/FXTSP(IGR)
            ELSE
               WRITE(400+IGR,3101)ICYCLE,XCG,YCG,ZCG, ! TRAJECTORY.dat
     &              ROTA1*RAD2DEG,ROTB1*RAD2DEG,CONEA*RAD2DEG,
     &              OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &              OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &              THRUST(IGR),TORQUE(IGR),ADV(IGR),UTSP(IGR),
     &              SHAFT(IGR),FDSP(IGR),OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT
               WRITE(500+IGR)ICYCLE,XCG,YCG,ZCG, ! TRAJECTORY.BIN
     &              ROTA1*RAD2DEG,ROTB1*RAD2DEG,CONEA*RAD2DEG,
     &              OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &              OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &              THRUST(IGR),TORQUE(IGR),ADV(IGR),UTSP(IGR),
     &              SHAFT(IGR),FDSP(IGR),OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT
            ENDIF


         ELSEIF (TRMODE == 10) THEN
C ... Search a correct row (TRAJECTORY_IGR.dat)
 92         CONTINUE
            READ(400+IGR,'(1X,I10)',END=93) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 93
            ELSE
               GOTO 92
            ENDIF
 93         CONTINUE
            BACKSPACE(400+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,3110)ICYCLE,SHAFT(IGR)*RAD2DEG,
     &           -PHIR*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
C ... Search a correct row (TRAJECTORY_IGR.BIN)
 95         CONTINUE
            READ(500+IGR,END=96) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 96
            ELSE
               GOTO 95
            ENDIF
 96         CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to TRAJECTORY_IGR.BIN file
            WRITE(500+IGR)ICYCLE,SHAFT(IGR)*RAD2DEG,
     &           -PHIR*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT    
            
         ENDIF                  ! (TRMODE > 10.AND.TRMODE < 21 ...


      END IF                    ! IF (.NOT. TIMEL) ...

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... TIME ACCURATE CALCULATION
C ... Three cases: FREE MOVEMENT (TRMODE = 6), FORCED MOVEMENT (TRMODE = 61)
C ... or EJECTOR MOVEMENT (TRMODE = 62)
      IF (TIMEL) THEN                                   



C ... FREE MOVEMENT
         IF (TRMODE == 6) THEN
C ... Ejector forces and moments are set to "0"
            OSKU(IGR)%FXE = 0
            OSKU(IGR)%FYE = 0
            OSKU(IGR)%FZE = 0
            OSKU(IGR)%MXE = 0
            OSKU(IGR)%MYE = 0
            OSKU(IGR)%MZE = 0
C ... Search the correct row = the last row 
 100        CONTINUE                                                     
            READ(400+IGR,'(11X,F10.6)',END=200)
            GOTO 100
C ... Continue
 200        CONTINUE
            BACKSPACE(400+IGR)
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2000)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
         ENDIF                  ! IF (TRMODE == 6


C ... FORCED MOVEMENT                                                      
         IF (TRMODE == 61) THEN 
C ... Ejector forces and moments are set to "0"
            OSKU(IGR)%FXE = 0
            OSKU(IGR)%FYE = 0
            OSKU(IGR)%FZE = 0
            OSKU(IGR)%MXE = 0
            OSKU(IGR)%MYE = 0
            OSKU(IGR)%MZE = 0
C ... Search the correct row 
 300        CONTINUE
            READ(400+IGR,'(11X,F10.6)',END=500) TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN 
               GOTO 400               
            ELSE 
               GOTO 300               
            ENDIF               ! IF ((ABS(TIME-TDBL))/DTB <= 0.33333) ...
C ... Continue
 400        CONTINUE
            BACKSPACE(400+IGR)
 500        CONTINUE
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Calculate the amount of rows in the forced movement phase
            ITRAJEC = 0
 600        CONTINUE
            READ(400+IGR,'(1X,I10)',END=700) TRMODE2
            BACKSPACE(400+IGR)
            ITRAJEC = ITRAJEC+1
            READ(400+IGR,'(A190)') TRAJECF ! Read the trajectory data and
            TRAJECTF(ITRAJEC) = TRAJECF ! transfer it to the TRAJECTF-table
            IF (ITRAJEC >= 10000) THEN ! 10000 rows is maximum...
               WRITE(*,*)'Too many steps in the forced movement phase'
               WRITE(*,*)'maximum is 10000. STOP'
               STOP
            ENDIF
            IF (TRMODE2 == 6) THEN 
               GOTO 800               
            ELSE 
               GOTO 600               
            ENDIF               ! (TRMODE2 == 6) THEN
C ... Continue
 700        CONTINUE
            BACKSPACE(400+IGR)
 800        CONTINUE
            BACKSPACE(400+IGR)
C ... The cursor is moved to the correct row by using BACKSPACE
            DO I = 1,ITRAJEC-1
               BACKSPACE(400+IGR)
            ENDDO      
C ... Write calculated values to the TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2000)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
C ... Write information from the TRAJECTF-table back to the TRAJECTORY-file
            DO I = 2,ITRAJEC
               WRITE(400+IGR,'(A190)') TRAJECTF(I)
            ENDDO      
C ... ENDIF TRMODE == 61
         ENDIF                  ! (TRMODE == 61) THEN 


C ... EJECTOR MOVEMENT                                                     
         IF (TRMODE == 62) THEN 
C ... Search the correct row 
 900        CONTINUE
            READ(400+IGR,'(11X,F10.6)',END=1100) TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN 
               GOTO 1000               
            ELSE 
               GOTO 900               
            ENDIF               ! IF ((ABS(TIME-TDBL))/DTB <= 0.33333) ...
C ... Continue
 1000       CONTINUE
            BACKSPACE(400+IGR)
C ... Update ejector forces and moments to correct timelevel
            READ(400+IGR,'(A237,6F12.3)')HEADER,
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE 
c     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
c     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
            BACKSPACE(400+IGR)
 1100       CONTINUE
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Calculate the amount of rows in the ejector movement phase
            ITRAJEC = 0
 1200       CONTINUE
            READ(400+IGR,'(1X,I10)',END=1300) TRMODE2
            BACKSPACE(400+IGR)
            ITRAJEC = ITRAJEC+1
            READ(400+IGR,'(A439)') TRAJECM ! Read the trajectory data and
            TRAJECTM(ITRAJEC) = TRAJECM ! transfer it to the TRAJECTM-table
            IF (ITRAJEC >= 10000) THEN ! 10000 rows is maximum...
               WRITE(*,*)'Too many steps in the ejector movement'
               WRITE(*,*)'phase maximum is 10000. STOP'
               STOP
            ENDIF
            IF (TRMODE2 == 6) THEN 
               GOTO 1400               
            ELSE 
               GOTO 1200               
            ENDIF               ! (TRMODE2 == 6) THEN
C ... Continue
 1300       CONTINUE
            BACKSPACE(400+IGR)
 1400       CONTINUE
            BACKSPACE(400+IGR)
C ... The cursor is moved to the correct row by using BACKSPACE
            DO I = 1,ITRAJEC-1
               BACKSPACE(400+IGR)
            ENDDO      
C ... Write calculated values to the TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2000)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
C ... Write information from the TRAJECTM-table back to the TRAJECTORY-file
            DO I = 2,ITRAJEC
               WRITE(400+IGR,'(A439)') TRAJECTM(I)
            ENDDO      
C ... ENDIF TRMODE == 62
         ENDIF                  ! (TRMODE == 62) THEN 

         
         IF (TRMODE == 6 .OR. TRMODE == 61 .OR. TRMODE == 62) THEN
C ... Search the correct row ( the last row) in the TRAJECTORY_IGR.BIN file 
 1420       CONTINUE                                                     
            READ(500+IGR,END=1430)
            GOTO 1420
C ... Continue
 1430       CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to the TRAJECTORY_IGR.BIN file            
            WRITE(500+IGR)TRMODE,TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT            
         ENDIF
         


C ... HELICOPTER ROTOR CASE
         IF (TRMODE >= 31 .AND. TRMODE <= 39) THEN
C ... Search the correct row = the last row (TRAJECTORY.dat file)
 1500       CONTINUE                                                     
            READ(400+IGR,'(11X,F8.4)',END=1600)
            GOTO 1500
C ... Continue
 1600       CONTINUE
C ... Search the correct row = the last row (TRAJECTORY.BIN file)
 1610       CONTINUE                                                     
            READ(500+IGR,END=1620)
            GOTO 1610
C ... Continue
 1620       CONTINUE            
            BACKSPACE(400+IGR)
            BACKSPACE(500+IGR)
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,PHIR,THETAR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,2)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,3200)TDBL,PSIM*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &           OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
C ... Write calculated values to TRAJECTORY_IGR.BIN file
            WRITE(500+IGR)TDBL,PSIM*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &           OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &           OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &           OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE            
         ENDIF                  ! IF (TRMODE >= 31 .AND.

         
C ... SHIP CASE
         IF (TRMODE >= 21 .AND. TRMODE <= 29) THEN
C ... ROTX,ROTY,ROTZ FINFLO coordinate system -> Euler coordinate system
            CALL FINEUL(PSIR,THETAR,PHIR,OSKU(IGR)%ROTX,
     &           OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,P,Q,R,1)
C ... Search a correct row (TRAJECTORY_IGR.dat file)
 1640       CONTINUE
            READ(400+IGR,'(1X,E10.3)',END=1650) TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN 
               GOTO 1650
            ELSE
               GOTO 1640
            ENDIF
 1650       CONTINUE
            BACKSPACE(400+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,2550)TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
C ... Search a correct row (TRAJECTORY_IGR.BIN file)
 1653       CONTINUE
            READ(500+IGR,END=1655) TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
               GOTO 1655
            ELSE
               GOTO 1653
            ENDIF
 1655       CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to TRAJECTORY_IGR.BIN file
            WRITE(500+IGR)TDBL,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,R,
     &           OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &           OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT

C ... SHIP PROPULSION CASE 
         ELSEIF (TRMODE == 10) THEN
C ... Search a correct row (TRAJECTORY_IGR.dat)
 1660       CONTINUE
            READ(400+IGR,'(F10.6)',END=1670) TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN 
               GOTO 1670
            ELSE
               GOTO 1660
            ENDIF
 1670       CONTINUE
            BACKSPACE(400+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(400+IGR,3115)TDBL,SHAFT(IGR)*RAD2DEG,
     &           -PHIR*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT     
C            PSIM = PHIR         ! store PHIR to PSIM

C ... Search a correct row (TRAJECTORY_IGR.BIN)
 1680       CONTINUE
            READ(500+IGR,END=1690)  TIME
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN 
               GOTO 1690
            ELSE
               GOTO 1680
            ENDIF
 1690       CONTINUE
            BACKSPACE(500+IGR)
C ... Write calculated values to TRAJECTORY_IGR.dat file
            WRITE(500+IGR)TDBL,SHAFT(IGR)*RAD2DEG,
     &           -PHIR*RAD2DEG,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &           OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &           OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT     
            PSIM = PHIR ! store PHIR to PSIM
            
         ENDIF                  ! IF (TRMODE >= 21

C ... ENDIF
      ENDIF                     ! IF (TIMEL) THEN 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... Writing format      
 2000  FORMAT(1X,I10,F10.6,12F12.4,18E12.4) 
 2500  FORMAT(1X,I10,8F15.8,15F21.8) ! ship format
 2550  FORMAT(1X,F10.6,8F15.8,12F21.8) ! ship format
 3000  FORMAT(1X,I8,F11.4,12F12.4,18E12.4) ! helicopter format
 3100  FORMAT(1X,I10,9E13.5,12E15.7) ! actuator disk format in Patria's case
 3101  FORMAT(1X,I10,6E13.5,14E15.7) ! actuator disk format
 3110  FORMAT(1X,I10,17E13.5) ! ship propulsion format
 3115  FORMAT(1X,F10.6,17E13.5) ! ship propulsion format (ta)
 3200  FORMAT(1X,F8.4,F12.4,12F12.4,18E12.4) ! helicopter format

C ... Close channels (400+IGR) and (500+IGR)
      CLOSE(400+IGR)
      CLOSE(500+IGR)
      
      RETURN
      END SUBROUTINE TRAJECTORY_WRITE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT(TDBL,DTB,XCG,YCG,ZCG,PSIR,THETAR,
     &           PHIR,PSIM,OSKU,IGR,TIMEL,TRMODE,DRAUGHT,STARTSTEP)


      USE NS3CO       , ONLY : ALPHA,BETA,AGAMMA,ABANK,IOLD
      USE CONSTANTS   , ONLY : RAD2DEG, DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : SHIP,FLYOBJ

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE

      REAL :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,ROTX,ROTY,ROTZ,
     &        DRAUGHT,TDBL,DTB,TIME,P,Q,R

      CHARACTER(LEN=235) :: HEADER
      CHARACTER(LEN=3)   :: IGRNUM

      LOGICAL :: THERE, TIMEL, STARTSTEP

      ! January 2021: In order to make this subroutine more readable,
      !  I started to separate the main functions into their own
      !  subroutines. The new subroutine names start with TRAJECTORY_FORMAT_
      !  and they are placed after this one. I don't promise that they work 100%
      !  but I'm working on it ... Matti Palin

C ... THIS SUBROUTINE CREATES OR READS TRAJECTORY_IGR.dat FILE.
C ... FILE IS CREATED WHEN SIMULATION IS NON-TIME ACCURATE SIMULATION.
C ... IN TIME ACCURATE SIMULATION EARLIER CREATED TRAJECTORY FILE IS READ.
C ... BOTH CASES HAVE OWN IF STATEMENT (IF STATEMENT #1 AND #2).

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... IF STATEMENT #1
C ... Non-time accurate calculation

      IF (.NOT.TIMEL) THEN

         IF (FLYOBJ) THEN ! Flying object 
         ! Create a new TRAJECTORY_IGR.dat file, into the flying object format
         CALL NUMCH3(IGRNUM,IGR)

         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,*)' '  
         WRITE(400+IGR,*)'# FLYING OBJECT NUMBER ',IGR       
 5       FORMAT(1X,A,F7.2)
         WRITE(400+IGR,5)'# A/C aoa   [deg]              ',
     &   ALPHA*RAD2DEG
         WRITE(400+IGR,5)'# A/C beta  [deg]              ',
     &   BETA*RAD2DEG
         WRITE(400+IGR,5)'# A/C GAMMA [deg]              ',
     &   AGAMMA*RAD2DEG
         WRITE(400+IGR,5)'# A/C BANK  [deg]              ',
     &   ABANK*RAD2DEG
         WRITE(400+IGR,10)
 10      FORMAT(' '/,
     &      ' # TRAJECTORY MODE is type of the movement: FREE',
     &      ' MOVEMENT=6, FORCED MOVEMENT=61 (movements are given by',
     &      ' the user) and EJECTOR MOVEMENT=62 (FXE...MZE are',
     &      ' given by the user).'/,
     &      ' # x,y,z are locations of the Flying Object in the',
     &      ' FINFLO coordinate system.'/,
     &      ' # PSII,THETA,PHI are Euler angles of the',
     &      ' Flying Object, i.e. angles between the fixed coordinate',
     &      ' system (based on the FINFLO coordinate system)',
     &      ' and the particle body coordinate system.'/,
     &      ' # Vx,Vy,Vz are relative velocities of',
     &      ' the Flying Object compared to the A/C in the',
     &      ' FINFLO coordinate system.', 
     &      ' P,Q,R are angular velocities in the particle',
     &      ' body coordinate system.'/,
     &      ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &      ' moments in the FINFLO coordinate system.'
     &      ' Fxe,Fye,Fze,Mxe,Mye,Mze are ejector forces and moments'
     &      ' in the FINFLO coordinate system.'/,
     &      ' # Fxt,Fyt,Fzt,Mxt,Myt,Mzt are thrust forces and moments'
     &      ' in the FINFLO coordinate system.'/,/,
     &      ' # TRAJECTORY    time       x           y           z',
     &      '          PSII       THETA        PHI',
     &      '          Vx         Vy          Vz',
     &      '           P           Q           R',
     &      '           Fx          Fy          Fz',
     &      '           Mx          My          Mz',
     &      '          Fxe         Fye         Fze',
     &      '         Mxe         Mye         Mze',
     &      '         Fxt         Fyt         Fzt',
     &      '         Mxt         Myt         Mzt',/
     &      ' # MODE          [s]       [m]         [m]         [m]',
     &      '        [deg]       [deg]       [deg]',
     &      '        [m/s]      [m/s]       [m/s]',
     &      '       [rad/s]     [rad/s]     [rad/s]',
     &      '      [N]         [N]         [N]',
     &      '          [Nm]        [Nm]        [Nm]',
     &      '         [N]         [N]         [N]',
     &      '         [Nm]        [Nm]        [Nm]',
     &      '        [N]         [N]         [N]',
     &      '         [Nm]        [Nm]        [Nm]')

         ! Create a new TRAJECTORY_IGR.BIN file
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED')

         ENDIF !  FLYOBJ

	 SELECT CASE (TRMODE)

	 CASE(21:29)
	    CALL TRAJECTORY_FORMAT_SHIP1(OSKU,IGRNUM,IGR,
     &      TRMODE,XCG,YXG,ZCG,P,Q,R,PHIR,THETAR,PSIR,DRAUGHT)

	 CASE(31:39)		! Helicopter rotor 
	    CALL TRAJECTORY_FORMAT_HELICOPTERROTOR1(OSKU,
     &      IGR,P,Q,R,PHIR,THETAR,PSIR,PSIM,XCG,YCG,ZCG)

	 CASE(40:49)		! Actuator disk
            IF (IOLD < 1) THEN
	    CALL TRAJECTORY_FORMAT_ACTUATORDISC1(OSKU,IGRNUM,IGR,TRMODE)
            ELSEIF (IOLD >= 1) THEN
	    CALL TRAJECTORY_FORMAT_ACTUATORDISC2(OSKU,IGRNUM,IGR,TRMODE)
            ENDIF

         CASE(10) ! Ship propulsion 
            IF (IOLD < 1) THEN  ! Create a new TRAJECTORY_IGR.dat file
                                   ! into the helicopter rotor format
               CALL NUMCH3(IGRNUM,IGR)
               OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &            STATUS='UNKNOWN', FORM='FORMATTED')
               WRITE(400+IGR,*)''
               WRITE(400+IGR,*)''
               WRITE(400+IGR,'(A,I12)')' # PROPEL NUMBER         ',IGR
               WRITE(400+IGR,*)''
               WRITE(400+IGR,*)''
               WRITE(400+IGR,47)

 47            FORMAT(' '/,
     &            ' # OMEGA is rotation velocity',
     &            ' of the propel.'/,
     &            ' # ROTANG is rotation angle',
     &            ' of the propel.'/,
     &            ' # x,y,z are locations of the rotation point',
     &            ' in the FINFLO coordinate system.'/,
     &            ' # PSIR,THETAR,RHIR are attitude angles of the',
     &            ' rotation axel.'/,
     &            ' # PSIR is pitch angle pos. dir. bow down,',
     &            ' THETAR is yaw angle pos. dir. bow right,',
     &            ' PHIR is roll angle = negative ROTANG angle.'/,
     &            ' # Vx,Vy,Vz is the vector of rotation',
     &            ' in the FINFLO coordinate system.'/,
     &            ' # Fx,Fy,Fz,Mx,My,Mz are hydrodynamic forces and',
     &            ' moments in the FINFLO coordinate system.'/,/,
     &            ' # ICYCLE      OMEGA       ROTANG',
     &            '       x            y            z',
     &            '            PSIR        THETAR        PHIR',
     &            '        Vx           Vy           Vz',
     &            '           Fx           Fy           Fz',
     &            '           Mx           My           Mz',/
     &            ' #             [1/s]       [deg]       [m]',
     &            '          [m]          [m]',
     &            '          [deg]        [deg]        [deg]',
     &            '        [ ]          [ ]          [ ]',
     &            '          [N]          [N]          [N]',
     &            '          [Nm]         [Nm]         [Nm]')
		! Creates TRAJECTORY.BIN file
               OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &              STATUS='REPLACE', FORM='UNFORMATTED')
            ENDIF ! IF (IOLD < 1) THEN

	 END SELECT

C ... ENDIF IF STATEMENT #1
      ENDIF                     ! IF (.NOT.TIMEL ...

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... IF STATEMENT #2

      !Time accurate calculation 
      IF (TIMEL) THEN

         ! Check and select the correct part part of TIMEL-statement
         IF(TRMODE >= 31 .AND. TRMODE <= 39) THEN
            CALL TRAJECTORY_FORMAT_HELICOPTERROTOR2(OSKU,IGR,
     &           PSIR,THETAR,PHIR,PSIM,XCG,YCG,ZCG,TDBL,DTB)
            GOTO 600
         ELSEIF(TRMODE >= 21 .AND. TRMODE <= 29) THEN
	    CALL TRAJECTORY_FORMAT_SHIP2(OSKU,IGR,
     &           P,Q,R,XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,
     &           DRAUGHT,TDBL,DTB,TRMODE)
            GOTO 600
         ELSEIF(TRMODE == 10) THEN  ! SHIP PROPULSION CALCULATION
            CALL TRAJECTORY_FORMAT_SHIPPROPULSION(OSKU,
     &           IGR,TDBL,DTB,XCG,YCG,ZCG,P,Q,R,
     &           PHIR,THETAR,PSIR,PSIM)
            GOTO 600
         ELSE
            GOTO 50
         ENDIF


       
 50      CONTINUE ! TRAJECTORY CALCULATION
         ! Check if "TRAJECTORY_IGR" files exist and open channels 400+IGR 500+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.dat',EXIST=THERE)
         IF(THERE) THEN          
            OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.dat',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &         STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF
         
         ! Skip headerlines
         DO IHEADER=1,16
            READ(400+IGR,'(A235)') HEADER
         ENDDO

         ! Search a correct row (=correct timelevel)
 100     CONTINUE
         READ(400+IGR,'(1X,I10,F10.6)',END=200) TRMODE,TIME
         IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
            GOTO 300               
         ELSE 
            GOTO 100               
         ENDIF

	! Program will be stopped if the correct timelevel is not found.
	! In FORCED CALCULATION and EJECTOR CALCULATION next time step is 
	! set to be FREE MOVEMENT.

 200     CONTINUE
         IF (TRMODE == 6) THEN  ! FREE  MOVEMENT calculation 
            CONTINUE                     
            IF ((ABS(TIME-TDBL))/DTB > 0.33333) THEN
               WRITE(*,*)''
               WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &                   ' able to find timelevel',TDBL
               WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.dat'
               WRITE(*,*)''
               WRITE(*,*)' Exiting...'
               WRITE(*,*)''
               STOP
            ENDIF
         ELSEIF (TRMODE == 61.OR. TRMODE == 62) THEN ! FORCED or EJECTOR 
            IF ((ABS(TIME-(TDBL-DTB)))/DTB > 0.33333)THEN! MOVEMENT calculation
               WRITE(*,*)''
               WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &                   ' able to find timelevel',TDBL
               WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.dat'
               WRITE(*,*)''
               WRITE(*,*)' Exiting...'
               WRITE(*,*)''
               STOP
            ENDIF               ! IF ((ABS(TIME-TDBL))/DTB > 0.33333) ...
            BACKSPACE(400+IGR)
            BACKSPACE(400+IGR)
            READ(400+IGR,*)TRMODE,TIME,XCG,YCG,ZCG,
     &         PSIR,THETAR,PHIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &         OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &         OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 

            TRMODE        = 6

            GOTO 600
         ENDIF                  ! IF (TRMODE == 6 .OR.  ...

 300     CONTINUE
         BACKSPACE(400+IGR)
         ! Read values from the TRAJECTORY_IGR.dat file

         IF (TRMODE == 6) THEN      ! FREE MOVEMENT calculation
            READ(400+IGR,*)TRMODE,TIME,XCG,YCG,ZCG,
     &         PSIR,THETAR,PHIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &         OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &         OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT  

         ELSEIF (TRMODE == 61) THEN    ! FORCED MOVEMENT calculation
            READ(400+IGR,*)TRMODE,TIME,XCG,YCG,ZCG,
     &         PSIR,THETAR,PHIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R

C ... EJECTOR MOVEMENT calculation TRMODE is set to "62" 
C ... because of BACKSPACE

         ELSEIF (TRMODE == 62) THEN           
            BACKSPACE(400+IGR)
            READ(400+IGR,*)TRMODE,TIME,XCG,YCG,ZCG,
     &         PSIR,THETAR,PHIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &         OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &         OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT

            TRMODE=62		! Why is TRMODE set to 62 here, if we test for this earlier? terv. MPalin

         ENDIF                  ! IF (TRMODE == 6

         IF (STARTSTEP) THEN
            IF (TRMODE == 6 .OR. TRMODE == 61 .OR. TRMODE == 62) THEN
 330           CONTINUE
               ! Data is read from binary file
               READ(500+IGR)TRMODE,TIME,XCG,YCG,ZCG,  
     &            PSIR,THETAR,PHIR,
     &            OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &            OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &            OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &            OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &            OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &            OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &            OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
               IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
                  BACKSPACE(500+IGR)
                  WRITE(500+IGR)TRMODE,TIME,XCG,YCG,ZCG,  
     &            PSIR,THETAR,PHIR,
     &            OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &            OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &            OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &            OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &            OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &            OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &            OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
                  CONTINUE               
               ELSE 
                  GOTO 330               
               ENDIF
            ENDIF
         ENDIF

	 ! Write BACKSPACE and previously read data one more time, 
	 ! this is needed when simulation is started by using RSTART file 

540      FORMAT(1X,I10,F10.6,12F12.4,18E12.4) 

         IF (TRMODE == 6) THEN
            BACKSPACE(400+IGR)
            WRITE(400+IGR,540)
     &         TRMODE,TDBL,XCG,YCG,ZCG,
     &         PSIR,THETAR,PHIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,P,Q,R,
     &         OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &         OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE, 
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT 
         ENDIF ! IF (TRMODE == 6

600      CONTINUE

         ! Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
         IF (IOLD < 0 .AND. .NOT. SHIP 
     &  .OR. IOLD > 0 .AND. .NOT. SHIP) THEN
            PSIR   = PSIR*DEG2RAD
            THETAR = THETAR*DEG2RAD
            PHIR   = PHIR*DEG2RAD
            PSIM   = PSIM*DEG2RAD
            IF (TRMODE >= 31 .AND. TRMODE <= 39) THEN
               CALL EULFIN(PSIR,PHIR,THETAR,P,Q,R,
     &              ROTX,ROTY,ROTZ,2)
            ELSE
               CALL EULFIN(PSIR,THETAR,PHIR,P,Q,R,
     &              ROTX,ROTY,ROTZ,1)
            ENDIF
            OSKU(IGR)%ROTX = ROTX 
            OSKU(IGR)%ROTY = ROTY 
            OSKU(IGR)%ROTZ = ROTZ 
         ENDIF

C ... ENDIF IF STATEMENT #2
      ENDIF ! IF (TIMEL...

      !Close channels 400+IGR and 500+IGR
      CLOSE(400+IGR)
      CLOSE(500+IGR)
      
      RETURN
      END SUBROUTINE TRAJECTORY_FORMAT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_SHIP1(OSKU,IGRNUM,IGR,
     & TRMODE,XCG,YXG,ZCG,P,Q,R,PHIR,THETAR,PSIR,DRAUGHT)

      USE NS3CO       , ONLY : RE,REFVEL,CHLREF,
     &                         GX, GY, GZ, IOLD, ICYCLE
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : SINK,TRIMA,ZCGIS,ACTDISK
      USE CONSTANTS   , ONLY : DEG2RAD
      
      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*)

      INTEGER :: IGR,TRMODE,TRMODE2
      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,ROTX,ROTY,ROTZ,
     &           DRAUGHT,P,Q,R,A13,A14,A15
      REAL    :: XCG2,YCG2,ZCG2,PSIR2,VX2,VY2,VZ2,R2
      CHARACTER(LEN=3) :: IGRNUM
      LOGICAL :: THERE,TIMEL


 12   FORMAT('#'/,
     &   '# x,y,z are locations of the c.g. in the',
     &   ' FINFLO coordinate system.'/,
     &   '# TRIM is TRIM angle in the FINFLO'
     &   ' coordinate system (positive TRIM = bow down).'/,
     &   '# Vx,Vy,Vz are velocities of',
     &   ' the ship in the FINFLO coordinate system.'/, 
     &   '# Movements towards z direction represent',
     &   ' SINK direction (negative z = positive SINK).'/,
     &   '# TRIM vel is trimming velocity in the particle',
     &   ' body coordinate system.'/,
     &   '# Fx,Fy,Fz,Mx,My,Mz are hydrodynamic forces and'
     &   ' moments in the FINFLO coordinate system.'/,
     &   '# Fxt,Fyt,Fzt,Mxt,Myt,Mzt are thrust forces and'
     &   ' moments in the FINFLO coordinate system.'
     &   ' XFakep,YFakep,ZFakep are the location of'
     &   ' the fake propeller in the FINFLO'
     &   ' coordinate system.'/,'#'/,
     &   '#     ICYCLE      x              y',
     &   '              z             TRIM', 
     &   '           Vx             Vy',
     &   '             Vz             TRIM vel'
     &   '             Fx                   Fy',
     &   '                   Fz',
     &   '                    Mx',
     &   '                  My',          
     &   '                   Mz'
     &   '                  Fxt                  Fyt',
     &   '                  Fzt',
     &   '                   Mxt',
     &   '                 Myt',          
     &   '                  Mzt',
     &   '                   XFakep',
     &   '               YFakep',          
     &   '               ZFakep',/
     &   '#                [m]            [m]',
     &   '            [m]            [deg] ',
     &   '         [m/s]          [m/s]',
     &   '          [m/s]          [rad/s]',
     &   '              [N]                  [N]',
     &   '                  [N]',
     &   '                  [Nm]',
     &   '                [Nm]',        
     &   '                 [Nm]'
     &   '                 [N]                  [N]',
     &   '                  [N]',
     &   '                  [Nm]',
     &   '                [Nm]',        
     &   '                 [Nm]',
     &   '                   [m]',
     &   '                  [m]',        
     &   '                  [m]')  


      IF (IOLD == 0) THEN  
         ! Create a new TRAJECTORY_IGR.dat file for the ship case
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,'(A1)')'#'
         WRITE(400+IGR,'(A,I12)')'# SHIP NUMBER          ',IGR        
         WRITE(400+IGR,'(A,E12.5)')
     &      '# Reynolds number               ',RE
         WRITE(400+IGR,'(A,F12.6)')
     &      '# Froude number                 ',
     &      REFVEL/SQRT(SQRT(GX*GX+GY*GY+GZ*GZ)*CHLREF)
         WRITE(400+IGR,'(A,F12.6,A)')
     &      '# Draught                       ',DRAUGHT,' m'
         WRITE(400+IGR,12)

             
         ! Create a new TRAJECTORY_IGR.BIN file
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &           STATUS='REPLACE', FORM='UNFORMATTED')

         ELSEIF (IOLD >= 1 .AND. TRMODE == 22 .OR.
     &           IOLD >= 1 .AND. TRMODE == 23 .OR.
     &           IOLD >= 1 .AND. TRMODE == 21 .AND. ACTDISK) THEN
            ! Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
            CALL NUMCH3(IGRNUM,IGR)
            INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &         EXIST=THERE)
            IF(THERE) THEN          
               OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &         STATUS='UNKNOWN', FORM='UNFORMATTED')
            ELSE 
               WRITE(*,*)' '
               WRITE(*,*)' Warning:',
     &                   'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                   ' was not found'
               WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
               WRITE(*,*)' '
               WRITE(*,*)' Exiting...'
               WRITE(*,*)' '
               STOP
            ENDIF

            ! ... Skip headerlines
            !DO IHEADER=1,16
            !    READ(400+IGR,'(A235)') HEADER
            !ENDDO
            ! ... Search a correct row

 13         CONTINUE
            READ(500+IGR,END=14) TRMODE2
            IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
               GOTO 15
            ELSE
               GOTO 13
            ENDIF
 14         CONTINUE
            WRITE(*,*)''
            WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &                ' able to find cycle',ICYCLE
            WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.BIN'
            WRITE(*,*)''
            WRITE(*,*)' Exiting...'
            WRITE(*,*)''
            STOP
 15         CONTINUE 
            BACKSPACE(500+IGR)
            ! Read data in if not actdisk
            IF(TRMODE == 22 .OR. TRMODE == 23) THEN
               READ(500+IGR)TRMODE2,XCG,YCG,ZCG,PSIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &         R,OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &         A13,A14,A15

               ! Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
               PSIR   = PSIR*DEG2RAD ! Angles deg -> rad
               CALL EULFIN(PSIR,THETAR,PHIR,P,Q,R,
     &              ROTX,ROTY,ROTZ,1)
               OSKU(IGR)%ROTX = ROTX 
               OSKU(IGR)%ROTY = ROTY 
               OSKU(IGR)%ROTZ = ROTZ
               IF(TRMODE == 23) THEN ! FIXED SHIP
                   TRIMA(IGR) = PSIR
                   SINK(IGR)  = DRAUGHT+(ZCGIS(IGR)-ZCG)
               ENDIF
    
            ENDIF !TRMODE == 22 .OR ...
               ! Read data in if actdisk
               IF(TRMODE == 21) THEN
                  READ(500+IGR)TRMODE2,XCG2,YCG2,ZCG2,
     &               PSIR2,VX2,VY2,VZ2,R2,
     &               OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &               OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &               OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &               OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &               A13,A14,A15
               ENDIF

            ELSEIF (IOLD == -1) THEN
               ! Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
               CALL NUMCH3(IGRNUM,IGR)
               INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &                 EXIST=THERE)                  
               IF(THERE) THEN
                  OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &               STATUS='UNKNOWN', FORM='UNFORMATTED')
               ELSE 
                  WRITE(*,*)' '
                  WRITE(*,*)' Warning:',
     &                      'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                      ' was not found'
                  WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
                  WRITE(*,*)' '
                  WRITE(*,*)' Exiting...'
                  WRITE(*,*)' '
                  STOP
               ENDIF
	       ! Skip headerlines
	       !DO IHEADER=1,16
	       !   READ(400+IGR,'(A235)') HEADER
	       !ENDDO
	       ! Search the correct row = the last row 
 16            CONTINUE                                                     
               READ(500+IGR,END=17) TRMODE2
               GOTO 16
 17            CONTINUE 
               BACKSPACE(500+IGR)
               BACKSPACE(500+IGR)

               !Read data in
               READ(500+IGR)TRMODE2,XCG,YCG,ZCG,PSIR,
     &            OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &            R,OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &            OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &            OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &            OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     &            A13,A14,A15

               !Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
               PSIR   = PSIR*DEG2RAD
               CALL EULFIN(PSIR,THETAR,PHIR,P,Q,R,
     &              ROTX,ROTY,ROTZ,1)
               OSKU(IGR)%ROTX = ROTX 
               OSKU(IGR)%ROTY = ROTY 
               OSKU(IGR)%ROTZ = ROTZ
               TRIMA(IGR)     = PSIR
               SINK(IGR)      = DRAUGHT+(ZCGIS(IGR)-ZCG)

               ! Close channel 400+IGR
               ! CLOSE(400+IGR)

               !Create a new TRAJECTORY_IGR.dat file for the ship case
               OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &              STATUS='UNKNOWN', FORM='FORMATTED')
               WRITE(400+IGR,'(A1)')'#'
               WRITE(400+IGR,'(A,I12)')'# SHIP NUMBER          ',IGR        
               WRITE(400+IGR,'(A,E12.5)')
     &              '# Reynolds number               ',RE
               WRITE(400+IGR,'(A,F12.6)')
     &              '# Froude number                 ',
     &              REFVEL/SQRT(SQRT(GX*GX+GY*GY+GZ*GZ)*CHLREF)
               WRITE(400+IGR,'(A,F12.6,A)')
     &              '# Draught                       ',DRAUGHT,' m'
               WRITE(400+IGR,12)

            ENDIF               ! IF (IOLD == 0) THEN

      END SUBROUTINE TRAJECTORY_FORMAT_SHIP1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_SHIP2(OSKU,IGR,
     & P,Q,R,XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,
     & DRAUGHT,TDBL,DTB,TRMODE)

      USE NS3CO       , ONLY : RE,REFVEL,CHLREF,
     &                         GX, GY, GZ, IOLD
      USE CONSTANTS   , ONLY : DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : SINK,TRIMA,ZCGIS

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE,TRMODE2
      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,
     &           ROTX,ROTY,ROTZ,DRAUGHT,TDBL,DTB,TIME,P,Q,R

      CHARACTER(LEN=3) :: IGRNUM
      LOGICAL :: THERE

 481  FORMAT('#'/,
     &   '# x,y,z are locations of the c.g. in the',
     &   ' FINFLO coordinate system.'/,
     &   '# TRIM is TRIM angle in the FINFLO'
     &   ' coordinate system (positive TRIM = bow down).'/,
     &   '# Vx,Vy,Vz are velocities of',
     &   ' the ship in the FINFLO coordinate system.'/, 
     &   '# Movements towards z direction represent',
     &   ' SINK direction (negative z = positive SINK).'/,
     &   '# TRIM vel is trimming velocity in the particle',
     &   ' body coordinate system.'/,
     &   '# Fx,Fy,Fz,Mx,My,Mz are hydrodynamic forces and'
     &   ' moments in the FINFLO coordinate system.'/,
     &   '# Fxt,Fyt,Fzt,Mxt,Myt,Mzt are thrust forces and'
     &   ' moments in the FINFLO coordinate system.'/,'#'/,
     &   '#     t           x              y',
     &   '              z             TRIM', 
     &   '           Vx             Vy',
     &   '             Vz             TRIM vel'
     &   '             Fx                   Fy',
     &   '                   Fz',
     &   '                    Mx',
     &   '                  My',          
     &   '                   Mz'
     &   '                  Fxt                  Fyt',
     &   '                  Fzt',
     &   '                   Mxt',
     &   '                 Myt',          
     &   '                  Mzt',/
     &   '#    [s]         [m]            [m]',
     &   '            [m]            [deg] ',
     &   '         [m/s]          [m/s]',
     &   '          [m/s]          [rad/s]',
     &   '              [N]                  [N]',
     &   '                  [N]',
     &   '                  [Nm]',
     &   '                [Nm]',        
     &   '                 [Nm]'
     &   '                 [N]                  [N]',
     &   '                  [N]',
     &   '                  [Nm]',
     &   '                [Nm]',        
     &   '                 [Nm]')

      IF (IOLD == 0) THEN  
         ! Create a new TRAJECTORY_IGR.dat file for the ship case
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,'(A1)')'#'
         WRITE(400+IGR,'(A,I12)')'# SHIP NUMBER          ',IGR        
         WRITE(400+IGR,'(A,E12.5)')
     &      '# Reynolds number               ',RE
         WRITE(400+IGR,'(A,F12.6)')
     &      '# Froude number                 ',
     &      REFVEL/SQRT(SQRT(GX*GX+GY*GY+GZ*GZ)*CHLREF)
         WRITE(400+IGR,'(A,F12.6,A)')
     &      '# Draught                       ',DRAUGHT,' m'
         WRITE(400+IGR,481)


         ! Create a new TRAJECTORY_IGR.BIN file              
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='REPLACE', FORM='UNFORMATTED')
         
      ELSEIF (IOLD >= 1) THEN
         ! Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &         STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

	 ! ... Skip headerlines
	 !DO IHEADER=1,16
	 !   READ(400+IGR,'(A235)') HEADER
	 !ENDDO

         !Last row is selected if timelevel is zero
         IF (TDBL/DTB <= 0.33333) THEN
 484        CONTINUE
            READ(500+IGR,END=485)TRMODE2,XCG,YCG,ZCG,PSIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      R,OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &      OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
            GOTO 484
 485        CONTINUE
            GOTO 498
         ENDIF
         ! ... Search a correct row
         IF (TRMODE == 23) THEN
 491        CONTINUE
            READ(500+IGR,END=492) TIME

            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
               GOTO 493
            ELSE
               GOTO 491
            ENDIF

 492        CONTINUE
            WRITE(*,*)''
            WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &                ' able to find timelevel',T
            WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.BIN'
            WRITE(*,*)''
            WRITE(*,*)' Exiting...'
            WRITE(*,*)''
            STOP
 493        CONTINUE 
            BACKSPACE(500+IGR)
            ! ... Read values from the TRAJECTORY_IGR.dat file
            ! ... SHIP ROTOR calculation
            ! ... Read data in
            READ(500+IGR)TIME,XCG,YCG,ZCG,PSIR,
     &         OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &         R,OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &         OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT

            ! Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
            PSIR   = PSIR*DEG2RAD ! Angles deg -> rad
            CALL EULFIN(PSIR,THETAR,PHIR,P,Q,R,
     &           ROTX,ROTY,ROTZ,1)
            OSKU(IGR)%ROTX = ROTX 
            OSKU(IGR)%ROTY = ROTY 
            OSKU(IGR)%ROTZ = ROTZ
            TRIMA(IGR)     = PSIR
            SINK(IGR)      = DRAUGHT+(ZCGIS(IGR)-ZCG)
         ENDIF	! IF (TRMODE == 23) THEN 
                  
      ELSEIF (IOLD == -1) THEN
         ! Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &   EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN ',
     &                'was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF
         ! ... Skip headerlines
         !  DO IHEADER=1,16
         !     READ(400+IGR,'(A235)') HEADER
         !  ENDDO

         ! Search the correct row = the last row 
 496     CONTINUE                                                     
         READ(500+IGR,END=497) TIME
         GOTO 496
 497     CONTINUE 
         BACKSPACE(500+IGR)
         BACKSPACE(500+IGR)

         ! Read data in
         READ(500+IGR)TIME,XCG,YCG,ZCG,PSIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      R,OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &      OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
 498     CONTINUE
         ! Create a new TRAJECTORY.BIN file                  
         CLOSE(500+IGR)
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='REPLACE', FORM='UNFORMATTED')

         ! Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
         PSIR   = PSIR*DEG2RAD ! Angles deg -> rad
         CALL EULFIN(PSIR,THETAR,PHIR,P,Q,R,
     &        ROTX,ROTY,ROTZ,1)
         OSKU(IGR)%ROTX = ROTX 
         OSKU(IGR)%ROTY = ROTY 
         OSKU(IGR)%ROTZ = ROTZ
         TRIMA(IGR)     = PSIR
         SINK(IGR)      = DRAUGHT+(ZCGIS(IGR)-ZCG)
         ! Close channel 400+IGR
         !            CLOSE(400+IGR)
         ! Create a new TRAJECTORY_IGR.dat file for the ship case
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &   STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,'(A1)')'#'
         WRITE(400+IGR,'(A,I12)')'# SHIP NUMBER          ',IGR        
         WRITE(400+IGR,'(A,E12.5)')
     &   '# Reynolds number               ',RE
         WRITE(400+IGR,'(A,F12.6)')
     &   '# Froude number                 ',
     &   REFVEL/SQRT(SQRT(GX*GX+GY*GY+GZ*GZ)*CHLREF)
         WRITE(400+IGR,'(A,F12.6,A)')
     &   '# Draught                       ',DRAUGHT,' m'
         WRITE(400+IGR,481)
      ENDIF     ! IF (IOLD == 0) THEN

      END SUBROUTINE TRAJECTORY_FORMAT_SHIP2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_SHIPPROPULSION(OSKU,
     &           IGR,TDBL,DTB,XCG,YCG,ZCG,P,Q,R,
     &           PHIR,THETAR,PSIR,PSIM)

      USE NS3CO       , ONLY : IOLD
      USE CONSTANTS   , ONLY : DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE2
      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,
     &           TDBL,DTB,TIME,R,P
      CHARACTER(LEN=3) :: IGRNUM
      LOGICAL :: THERE

      IF (IOLD < 1 .OR. IOLD >= 1 .AND. 
     &   ABS(TDBL/DTB) <= 0.33333) THEN 
         !Check if "TRAJECTORY_IGR" file exists
         ! and then open channel 500+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      EXIST=THERE)
         IF(THERE) THEN
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         !Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         !Search the correct row = the last row 
 501     CONTINUE                                                     
         READ(500+IGR,END=502)
         GOTO 501
 502     CONTINUE 
         BACKSPACE(500+IGR)
         BACKSPACE(500+IGR)

         ! Read data in
         IF (IOLD < 0) THEN
            READ(500+IGR)TIME,R,P,XCG,YCG,ZCG,PSIR,THETAR,
     &      PHIR,OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &      OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
         ELSE
            READ(500+IGR)TRMODE2,R,P,XCG,YCG,ZCG,PSIR,THETAR,
     &      PHIR,OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &      OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT
         ENDIF
         PSIR   = PSIR*DEG2RAD ! update angles
         THETAR = THETAR*DEG2RAD
         PHIR   = PHIR*DEG2RAD                 
         PSIM   = PHIR

         ! Creates a new TRAJECOTORY.BIN file
         CLOSE(500+IGR)
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   STATUS='REPLACE', FORM='UNFORMATTED')                  
         !Creates a new TRAJECTORY_IGR.dat file, into the ship propel format
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &   STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,*)''
         WRITE(400+IGR,*)''
         WRITE(400+IGR,'(A,I12)')' # PROPEL NUMBER         ',IGR
         WRITE(400+IGR,*)''
         WRITE(400+IGR,*)''
         WRITE(400+IGR,505)
 505     FORMAT(' '/,
     &      ' # OMEGA is rotation velocity',
     &      ' of the propel.'/,
     &      ' # ROTANG is rotation angle',
     &      ' of the propel.'/,
     &      ' # x,y,z are locations of the rotation point',
     &      ' in the FINFLO coordinate system.'/,
     &      ' # PSIR,THETAR,RHIR are attitude angles of the',
     &      ' rotation axel.'/,
     &      ' # PSIR is pitch angle pos. dir. bow down,',
     &      ' THETAR is yaw angle pos. dir. bow right,',
     &      ' PHIR is roll angle = negative ROTANG angle.'/,
     &      ' # Vx,Vy,Vz is the vector of rotation',
     &      ' in the FINFLO coordinate system.'/,
     &      ' # Fx,Fy,Fz,Mx,My,Mz are hydrodynamic forces and',
     &      ' moments in the FINFLO coordinate system.'/,/,
     &      ' #   t         OMEGA       ROTANG',
     &      '       x            y            z',
     &      '            PSIR        THETAR        PHIR',
     &      '        Vx           Vy           Vz',
     &      '           Fx           Fy           Fz',
     &      '           Mx           My           Mz',/
     &      ' #  [s]        [1/s]       [deg]       [m]',
     &      '          [m]          [m]',
     &      '          [deg]        [deg]        [deg]',
     &      '        [ ]          [ ]          [ ]',
     &      '          [N]          [N]          [N]',
     &      '          [Nm]         [Nm]         [Nm]')

      ELSEIF (IOLD >= 1 .AND. ABS(TDBL/DTB) > 0.33333) THEN 
         !Check if "TRAJECTORY_IGR" file exists and then open channel 500+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         !Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         !Search a correct row
 506     CONTINUE
         READ(500+IGR,END=507) TIME
         IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
            GOTO 508
         ELSE
            GOTO 506
         ENDIF

 507     CONTINUE
         WRITE(*,*)''
         WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &             ' able to find timelevel',TDBL
         WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.dat'
         WRITE(*,*)''
         WRITE(*,*)' Exiting...'
         WRITE(*,*)''
         STOP
 508     CONTINUE 
         BACKSPACE(500+IGR)
         !Read data in
         READ(500+IGR)TIME,R,P,XCG,YCG,ZCG,PSIR,THETAR,
     &      PHIR,OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &      OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT

         PSIR   = PSIR*DEG2RAD ! update angles
         THETAR = THETAR*DEG2RAD
         PHIR   = PHIR*DEG2RAD
         PSIM   = PHIR

      ENDIF !IF (IOLD < 1 .OR. IOLD >= 1...   

      !Writing formats
      !540  FORMAT(1X,I10,F10.6,12F12.4,18E12.4)  ! Not used?
      !550  FORMAT(1X,F8.4,F12.4,12F12.4,18E12.4) 

      END SUBROUTINE TRAJECTORY_FORMAT_SHIPPROPULSION
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_HELICOPTERROTOR1(OSKU,
     &IGR,P,Q,R,PHIR,THETAR,PSIR,PSIM,XCG,YCG,ZCG)

      USE NS3CO       , ONLY : ALPHA,BETA,IOLD, ICYCLE
      USE CONSTANTS   , ONLY : RAD2DEG, DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF

      USE BLADE_VARIABLES , ONLY : DOTZE,DOTBETA,DOTTHETA

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE,TRMODE2
      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,
     &           PSIM,ROTX,ROTY,ROTZ,P,Q,R

      CHARACTER(LEN=3) :: IGRNUM
      LOGICAL :: THERE

 20   FORMAT(' '/,
     &   ' # PSII(shaft) is shaft angle in the flight',
     &   ' mechanic coordinate system.'/,
     &   ' # x,y,z are locations of the rotation point',
     &   ' in the FINFLO coordinate system.'/,
     &   ' # PSII,THETA,BETA are attitude angles',
     &   ' of the blade, i.e.',
     &   ' angles between the fixed coordinate',
     &   ' system (=FINFLO coordinate system)',
     &   ' and the particle body',
     &   ' coordinate system (Rotation',
     &   ' order: PSII,BETA,THETA).'/,
     &   ' # Vx,Vy,Vz are relative velocities of',
     &   ' the blade in the FINFLO coordinate system.', 
     &   ' DOTZETA,DOTBETA,DOTTHETA are angular velocities',
     &   ' in the helicopter coordinate system.'/,
     &   ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &   ' moments in the FINFLO coordinate system.'
     &   ' Fx_hin,Fy_hin,Fz_hin,Mx_hin,My_hin,Mz_hin'
     &   ' are spring and damper forces'
     &   ' and moments in the FINFLO coordinate system.'/,
     &   ' # Fx_cen,Fy_cen,Fz_cen,Mx_cen,My_cen,Mz_cen'
     &   ' are centrifugal'
     &   ' forces and moments in'
     &   ' the FINFLO coordinate system.'/,/,
     &   ' # t [s]     PSII(shaft)',
     &   '    x           y           z',
     &   '         PSII       THETA        PHI',
     &   '          Vx         Vy          Vz',
     &   '        DOTZETA     DOTBETA    DOTTHETA',
     &   '        Fx          Fy          Fz',
     &   '           Mx          My          Mz',
     &   '       Fx_hin      Fy_hin      Fz_hin',
     &   '       Mx_hin      My_hin      Mz_hin',
     &   '     Fx_cen      Fy_cen      Fz_cen',
     &   '       Mx_cen      My_cen      Mz_cen',/
     &   ' # or Cycle   [deg]        [m]',
     &   '         [m]         [m]',
     &   '       [deg]       [deg]       [deg]',
     &   '        [m/s]      [m/s]       [m/s]',
     &   '       [rad/s]     [rad/s]     [rad/s]',
     &   '      [N]         [N]         [N]',
     &   '          [Nm]        [Nm]        [Nm]',
     &   '         [N]         [N]         [N]',
     &   '         [Nm]        [Nm]        [Nm]',
     &   '        [N]         [N]         [N]',
     &   '         [Nm]        [Nm]        [Nm]')


      IF (IOLD == 0) THEN 

         !Create a new TRAJECTORY_IGR.dat file, into the helicopter rotor format
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &   STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,*)' '
         WRITE(400+IGR,*)'# BLADE NUMBER         ',IGR        
 19      FORMAT(1X,A,F7.2)
         WRITE(400+IGR,19)'# H/C aoa   [deg]              ',
     &   ALPHA*RAD2DEG
         WRITE(400+IGR,19)'# H/C beta  [deg]              ',
     &   BETA*RAD2DEG
         WRITE(400+IGR,*)''
         WRITE(400+IGR,*)''
         WRITE(400+IGR,20)

         !Create a new TRAJECTORY_IGR.BIN file
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   STATUS='UNKNOWN', FORM='UNFORMATTED')
               
      ELSEIF (IOLD >= 1) THEN
         !Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         ! ... Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         !Search a correct row
 23      CONTINUE
         READ(500+IGR,END=24) TRMODE2

         IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
           GOTO 25
         ELSE
           GOTO 23
         ENDIF
 24      CONTINUE
         WRITE(*,*)''
         WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &             ' able to find cycle',ICYCLE
         WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.BIN'
         WRITE(*,*)''
         WRITE(*,*)' Exiting...'
         WRITE(*,*)''
         STOP
 25      CONTINUE 
         BACKSPACE(500+IGR)

         !Read data in
         READ(500+IGR)TRMODE2,PSIM,
     &      XCG,YCG,ZCG,PSIR,THETAR,PHIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &      OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &      OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &      OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &      OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &      OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE

         !Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
         PSIR   = PSIR*DEG2RAD
         THETAR = THETAR*DEG2RAD
         PHIR   = PHIR*DEG2RAD
         PSIM   = PSIM*DEG2RAD
         CALL EULFIN(PSIR,PHIR,THETAR,P,Q,R,
     &        ROTX,ROTY,ROTZ,2)
         OSKU(IGR)%ROTX = ROTX 
         OSKU(IGR)%ROTY = ROTY 
         OSKU(IGR)%ROTZ = ROTZ 

      ELSEIF (IOLD == -1) THEN
         !Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      EXIST=THERE)
         IF(THERE) THEN          
           OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         ! ... Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         !Search the correct row = the last row 
 26      CONTINUE                                                     
         READ(500+IGR,END=27)
         GOTO 26
 27      CONTINUE 
         BACKSPACE(500+IGR)
         BACKSPACE(500+IGR)

         ! Read data in
         READ(500+IGR)TRMODE2,PSIM,
     &        XCG,YCG,ZCG,PSIR,THETAR,PHIR,
     &        OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &        DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &        OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &        OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,     
     &        OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &        OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &        OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &        OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE

         !Euler coordinate system (P,Q,R) -> FINFLO coordinate system 
         PSIR   = PSIR*DEG2RAD ! Angles deg -> rad
         THETAR = THETAR*DEG2RAD
         PHIR   = PHIR*DEG2RAD
         PSIM   = PSIM*DEG2RAD
         CALL EULFIN(PSIR,PHIR,THETAR,P,Q,R,
     &        ROTX,ROTY,ROTZ,2)
         OSKU(IGR)%ROTX = ROTX 
         OSKU(IGR)%ROTY = ROTY 
         OSKU(IGR)%ROTZ = ROTZ 
         !Creates a new TRAJECTORY.BIN file
         CLOSE(500+IGR)
         OPEN(500+IGR,FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='REPLACE', FORM='UNFORMATTED')
               
         !Create a new TRAJECTORY_IGR.dat file for the helicopter case
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED')
         !Create a new TRAJECTORY_IGR.dat file, into the helicopter rotor format
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &      STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,*)' '
         WRITE(400+IGR,*)'# BLADE NUMBER         ',IGR        
 29      FORMAT(1X,A,F7.2)
         WRITE(400+IGR,29)'# H/C aoa   [deg]              ',
     &      ALPHA*RAD2DEG
         WRITE(400+IGR,29)'# H/C beta  [deg]              ',
     &      BETA*RAD2DEG
         WRITE(400+IGR,*)''
         WRITE(400+IGR,*)''
         WRITE(400+IGR,20)
      ENDIF ! IF (IOLD == 0) THEN

      END SUBROUTINE TRAJECTORY_FORMAT_HELICOPTERROTOR1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_HELICOPTERROTOR2(OSKU,IGR,
     & PSIR,THETAR,PHIR,PSIM,XCG,YCG,ZCG,TDBL,DTB)

      USE NS3CO       , ONLY : ALPHA,BETA,IOLD
      USE CONSTANTS   , ONLY : RAD2DEG
      USE TYPE_ARRAYS , ONLY : SIXDOF

      USE BLADE_VARIABLES , ONLY : DOTZE,DOTBETA,DOTTHETA

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE2

      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,TDBL,DTB,TIME

      CHARACTER(LEN=235) :: HEADER
      CHARACTER(LEN=3)   :: IGRNUM
      LOGICAL :: THERE

C ... HELICOPTER ROTOR CALCULATION

 550  FORMAT(1X,F8.4,F12.4,12F12.4,18E12.4) 

 405  FORMAT(' '/,
     &   ' # PSII(shaft) is shaft angle in the flight',
     &   ' mechanic coordinate system.'/,
     &   ' # x,y,z are locations of the rotation point',
     &   ' in the FINFLO coordinate system.'/,
     &   ' # PSII,THETA,BETA are attitude angles',
     &   ' of the blade, i.e.',
     &   ' angles between the fixed coordinate',
     &   ' system (=FINFLO coordinate system)',
     &   ' and the particle body',
     &   ' coordinate system (Rotation',
     &   ' order: PSII,BETA,THETA).'/,
     &   ' # Vx,Vy,Vz are relative velocities of',
     &   ' the blade in the FINFLO coordinate system.', 
     &   ' DOTZETA,DOTBETA,DOTTHETA are angular velocities',
     &   ' in the helicopter coordinate system.'
     &   ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &   ' moments in the FINFLO coordinate system.'
     &   ' Fx_hin,Fy_hin,Fz_hin,Mx_hin,My_hin,Mz_hin'
     &   ' are spring and damper forces'
     &   ' and moments in the FINFLO coordinate system.'/,
     &   ' # Fx_cen,Fy_cen,Fz_cen,Mx_cen,My_cen,Mz_cen'
     &   ' are centrifugal'
     &   ' forces and moments in'
     &   ' the FINFLO coordinate system.'/,/,
     &   ' # t [s]     PSII(shaft)',
     &   '   x           y           z',
     &   '          PSII       THETA        PHI',
     &   '          Vx         Vy          Vz',
     &   '        DOTZETA     DOTBETA    DOTTHETA',
     &   '        Fx          Fy          Fz',
     &   '           Mx          My          Mz',
     &   '       Fx_hin      Fy_hin      Fz_hin',
     &   '       Mx_hin      My_hin      Mz_hin',
     &   '     Fx_cen      Fy_cen      Fz_cen',
     &   '       Mx_cen      My_cen      Mz_cen',/
     &   ' # or Cycle   [deg]       [m]',
     &   '         [m]         [m]',
     &   '        [deg]       [deg]       [deg]',
     &   '        [m/s]      [m/s]       [m/s]',
     &   '       [rad/s]     [rad/s]     [rad/s]',
     &   '      [N]         [N]         [N]',
     &   '          [Nm]        [Nm]        [Nm]',
     &   '         [N]         [N]         [N]',
     &   '         [Nm]        [Nm]        [Nm]',
     &   '        [N]         [N]         [N]',
     &   '         [Nm]        [Nm]        [Nm]')

      IF (IOLD == 0) THEN
         !Create a new TRAJECTORY_IGR.dat file, into the helicopter rotor format
         CALL NUMCH3(IGRNUM,IGR)
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &        STATUS='UNKNOWN', FORM='FORMATTED')
         WRITE(400+IGR,*)' '
         WRITE(400+IGR,*)'# BLADE NUMBER         ',IGR        
 404     FORMAT(1X,A,F7.2)
         WRITE(400+IGR,404)'# H/C aoa   [deg]              ',
     &        ALPHA*RAD2DEG
         WRITE(400+IGR,404)'# H/C beta  [deg]              ',
     &        BETA*RAD2DEG
         WRITE(400+IGR,*)''
         WRITE(400+IGR,*)''
         WRITE(400+IGR,405)

         !Creates a new TRAJECTORY_IGR.BIN file               
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   STATUS='REPLACE', FORM='UNFORMATTED') 

      !If IOLD is not zero then old TRAJECTORY file is read and manipulated
      ELSEIF (IOLD /= 0) THEN
         !Check if "TRAJECTORY_IGR" file exists and then open channel 500+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF
         ! ... Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         ! Last row is selected if timelevel is zero


         IF (TDBL/DTB <= 0.33333) THEN
 410        CONTINUE
            READ(500+IGR,END=420) TRMODE2

            GOTO 410
 420        BACKSPACE(500+IGR)
            GOTO 450
         ENDIF

         ! Search a correct row (=correct timelevel)
 430     CONTINUE
         READ(500+IGR,END=440) TIME
         IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
            GOTO 450               
         ELSE 
            GOTO 430               
         ENDIF
 440     CONTINUE   

                                   
         IF ((ABS(TIME-TDBL))/DTB > 0.33333) THEN
            WRITE(*,*)''
            WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &           ' able to find timelevel',TDBL
            WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.BIN'
            WRITE(*,*)''
            WRITE(*,*)' Exiting...'
            WRITE(*,*)''
            STOP
         ENDIF      

 450     CONTINUE

         ! Read values from the TRAJECTORY_IGR.dat file
         ! HELICOPTER ROTOR calculation
         BACKSPACE(500+IGR)
         IF (IOLD > 0 .AND. TDBL/DTB <= 0.33333) THEN
            READ(500+IGR)TRMODE2,PSIM,XCG,YCG,ZCG,
     &      PSIR,THETAR,PHIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &      OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &      OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH,
     &      OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &      OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
         ELSE
            READ(500+IGR)TIME,PSIM,XCG,YCG,ZCG,
     &      PSIR,THETAR,PHIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &      OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &      OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH,
     &      OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &      OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE   
         ENDIF

         !Closes and opens TRAJECTORY file if time level is zero
         IF (TDBL/DTB <= 0.33333) THEN
            CLOSE(400+IGR)
            CLOSE(500+IGR)
            OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &           STATUS='UNKNOWN', FORM='FORMATTED')
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &           STATUS='REPLACE', FORM='UNFORMATTED') 
            DO IHEADER=1,16
               READ(400+IGR,'(A235)') HEADER
            ENDDO

            !Write data to the TRAJECTORY.dat file TRAJECTORY.BIN file
            WRITE(400+IGR,550)TDBL,PSIM,XCG,YCG,ZCG,
     &      PSIR,THETAR,PHIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &      OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &      OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &      OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &      OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
            WRITE(500+IGR)TDBL,PSIM,XCG,YCG,ZCG,
     &      PSIR,THETAR,PHIR,
     &      OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &      DOTZE(IGR),DOTBETA(IGR),DOTTHETA(IGR),
     &      OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     &      OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     &      OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     &      OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH, 
     &      OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     &      OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE
         ENDIF
      ENDIF ! (IF (IOLD == 0
      END SUBROUTINE TRAJECTORY_FORMAT_HELICOPTERROTOR2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_ACTUATORDISC1(OSKU,IGRNUM,IGR,TRMODE)

      USE NS3CO       , ONLY : IOLD, ICYCLE
      USE CONSTANTS   , ONLY : DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : ROTA1,ROTB1,CONEA,
     &                         ADV,FDSP,UTSP,THRUST,
     &                         TORQUE,RTMSP,FXTSP,QFACT
      USE BLADE_VARIABLES,ONLY : SHAFT,THETACOL,THCYCLON,THCYCLAT

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE
      REAL    :: XCG,YCG,ZCG,A1,TX,TY,TZ,P,Q,R
      REAL    :: A11,A12
      CHARACTER(LEN=235) :: HEADER
      CHARACTER(LEN=3)   ::  IGRNUM
      LOGICAL :: THERE

 40   FORMAT(' #'/,
     &   ' # x,y,z are locations of the actuator disk in the',
     &   ' FINFLO coordinate system.'/,
     &   ' # A0 is longitudinal angle of the actuator disk.',
     &   ' B0 is lateral angle of the actuator disk.',
     &   ' CON is conic angle of the the actuator disk.'/,
     &   ' # A0,BO and CON are defined in the FINFLO input.'/,
     &   ' # THETACOL, THETACYCLON and THETACYCLAT are THETA',
     &   ' angles iterated by the blade elment theory.'/,
     &   ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic/hydrodynamic'
     &   ' forces and moments of the actuator disk'
     &   ' in the FINFLO coordinate system.'/,
     &   ' # T_c is target thrust for the the actuator',
     &   ' disk in the blade element iteration.',
     &   ' T_c is defined in the FINFLO input.',
     &   ' Thrust is thrust of the actuator disk.'/,
     &   ' # Torque is torque of the actuator',
     &   ' disk. QFACT is a multiplier for tangetial',
     &   ' force components. DAMPN and DAMPT are factors in',
     &   ' the PD controller.'/,
     &   ' # '/,
     &   ' #    ICYCLE      x            y',
     &   '            z          A0', 
     &   '          B0           CON',
     &   '        THETACOL     THETACYCLON',
     &   '  THETACYCLAT',
     &   '     Fx             Fy',
     &   '             Fz',
     &   '             Mx',
     &   '             My',          
     &   '             Mz',
     &   '            T_c',
     &   '            Thrust'
     &   '          Torque'
     &   '          QFACT',
     &   '          DAMPN',
     &   '          DAMPT'/
     &   ' #               [m]          [m]',
     &   '          [m]        [deg] ',
     &   '      [deg]        [deg]',
     &   '        [deg]        [deg]',
     &   '        [deg]',
     &   '          [N]            [N]',
     &   '            [N]',
     &   '           [Nm]',
     &   '           [Nm]',        
     &   '           [Nm]',
     &   '           [N]',
     &   '            [N]',
     &   '            [Nm]',
     &   '                                        [N/deg]')

 41   FORMAT(' #'/,' #'/,
     &   ' # x,y,z are locations of the actuator disk in the',
     &   ' FINFLO coordinate system.'/,
     &   ' # A0 is longitudinal angle of the actuator disk.',
     &   ' B0 is longitudinal angle of the actuator disk.',
     &   ' CON is conic angle of the the actuator disk.'/,
     &   ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic/hydrodynamic'
     &   ' forces and moments of the actuator disk'
     &   ' in the FINFLO coordinate system.'/,
     &   ' # T is thrust of the the actuator disk.',
     &   ' Q is torque of the the actuator disk.'/,
     &   ' # J is advance ratio of the actuator disk. UT is',
     &   ' average flow velocity through the actuator disk.'/,
     &   ' # N is propeller rotation speed. FD is skin',
     &   ' friction corretion force.'/,
     &   ' #'/,
     &   ' #    ICYCLE      x            y',
     &   '            z          A0', 
     &   '          B0           CON',
     &   '           Fx             Fy',
     &   '             Fz',
     &   '             Mx',
     &   '             My',          
     &   '             Mz',
     &   '             T',
     &   '              Q',
     &   '              J',
     &   '              UT',
     &   '             N',
     &   '              FD',
     &   '             DAMPN',
     &   '          DAMPT'/
     &   ' #               [m]          [m]',
     &   '          [m]        [deg] ',
     &   '      [deg]        [deg]',
     &   '          [N]            [N]',
     &   '            [N]',
     &   '           [Nm]',
     &   '           [Nm]',        
     &   '           [Nm]',
     &   '           [N]',        
     &   '            [Nm]',
     &   '               ',        
     &   '           [m/s]',
     &   '          [1/s]',        
     &   '          [N]',
     &   '               ',
     &   '               ')


      IF (IOLD < 0) THEN
         !Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
         CALL NUMCH3(IGRNUM,IGR)
         INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   EXIST=THERE)
         IF(THERE) THEN          
            OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &      STATUS='UNKNOWN', FORM='UNFORMATTED') 
         ELSE 
            WRITE(*,*)' '
            WRITE(*,*)' Warning:',
     &                'File TRAJECTORY_'//IGRNUM//'.BIN',
     &                ' was not found'
            WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
            WRITE(*,*)' '
            WRITE(*,*)' Exiting...'
            WRITE(*,*)' '
            STOP
         ENDIF

         ! ... Skip headerlines
         !DO IHEADER=1,16
         !   READ(400+IGR,'(A235)') HEADER
         !ENDDO

         !Data in from the last line
 31      CONTINUE                                                     
         READ(500+IGR,END=32)
         GOTO 31
 32      CONTINUE 
         BACKSPACE(500+IGR)
         BACKSPACE(500+IGR)
         IF(TRMODE == 40) THEN ! Patria's helicopt.
            READ(500+IGR)ICYCLE,XCG,YCG,ZCG,
     &         ROTA1(IGR),ROTB1(IGR),CONEA(IGR),THETACOL(IGR),
     &         THCYCLON(IGR),THCYCLAT(IGR),
     &         OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &         OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,     
     &         THRUST(IGR),A1,TORQUE(IGR),A11,A12
            CALL EULFIN(0.0,ROTB1(IGR),ROTA1(IGR),
     &         0.0,0.0,-SHAFT(IGR)/ABS(SHAFT(IGR)),
     &         TX,TY,TZ,2)

            FXTSP(IGR) = (OSKU(IGR)%FXT*TX+OSKU(IGR)%FYT*TY
     &      +OSKU(IGR)%FZT*TZ)

            RTMSP(IGR) = OSKU(IGR)%MZT/FXTSP(IGR)
            FDSP(IGR)  = OSKU(IGR)%MXT/FXTSP(IGR)
         ELSE
            READ(500+IGR)ICYCLE,XCG,YCG,ZCG,
     &      ROTA1(IGR),ROTB1(IGR),CONEA(IGR),OSKU(IGR)%FXT,
     &      OSKU(IGR)%FYT,OSKU(IGR)%FZT,OSKU(IGR)%MXT,
     &      OSKU(IGR)%MYT,OSKU(IGR)%MZT,THRUST(IGR),     
     &      TORQUE(IGR),ADV(IGR),UTSP(IGR),SHAFT(IGR),
     &      FDSP(IGR),OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT 
         ENDIF
         THETACOL(IGR) = THETACOL(IGR)*DEG2RAD
         THCYCLON(IGR) = THCYCLON(IGR)*DEG2RAD
         THCYCLAT(IGR) = THCYCLAT(IGR)*DEG2RAD
         ROTA1(IGR) = ROTA1(IGR)*DEG2RAD
         ROTB1(IGR) = ROTB1(IGR)*DEG2RAD
         CONEA(IGR) = CONEA(IGR)*DEG2RAD
         CLOSE(500+IGR)
      ENDIF

      !Create a new TRAJECTORY_IGR.dat file, into the actuator disk format
      CALL NUMCH3(IGRNUM,IGR)
      OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(400+IGR,*)'#'
      WRITE(400+IGR,*)'#'
      WRITE(400+IGR,*)'#' 
      WRITE(400+IGR,*)'# ACTUATOR DISK NUMBER         ',IGR        
      WRITE(400+IGR,*)'#'

      IF(TRMODE == 40) THEN
         WRITE(400+IGR,40)
      ELSE
         WRITE(400+IGR,41)
      ENDIF

      !Creates TRAJECTORY.BIN file
      OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   STATUS='REPLACE', FORM='UNFORMATTED')
               
      END SUBROUTINE TRAJECTORY_FORMAT_ACTUATORDISC1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_FORMAT_ACTUATORDISC2(OSKU,IGRNUM,IGR,TRMODE)

      USE NS3CO       , ONLY : ICYCLE
      USE CONSTANTS   , ONLY : DEG2RAD
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT      , ONLY : ROTA1,ROTB1,CONEA,
     &                         ADV,FDSP,UTSP,IFA,THRUST,
     &                         TORQUE,RTMSP,FXTSP
      USE BLADE_VARIABLES , ONLY : SHAFT,THETACOL,THCYCLON,THCYCLAT

      TYPE(SIXDOF), DIMENSION(19) :: OSKU(*) 

      INTEGER :: IGR,TRMODE,TRMODE2

      REAL    :: XCG,YCG,ZCG,PHIR,THETAR,PSIR,PSIM,
     &           A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,TX,TY,TZ,
     &           TDBL,DTB
      REAL    :: A11,A12

      CHARACTER(LEN=3)   ::  IGRNUM

      LOGICAL :: THERE

      !Check if "TRAJECTORY_IGR" file exists and then open channel 400+IGR
      CALL NUMCH3(IGRNUM,IGR)
      INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   EXIST=THERE)
      IF(THERE) THEN          
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.BIN',
     &   STATUS='UNKNOWN', FORM='UNFORMATTED') 
      ELSE 
         WRITE(*,*)' '
         WRITE(*,*)' Warning:',
     &             'File TRAJECTORY_'//IGRNUM//'.BIN',
     &             ' was not found'
         WRITE(*,*)' in TRAJECTORY_FORMAT Subroutine.' 
         WRITE(*,*)' '
         WRITE(*,*)' Exiting...'
         WRITE(*,*)' '
         STOP
      ENDIF

      ! ... Skip headerlines
      !DO IHEADER=1,16
      !   READ(400+IGR,'(A235)') HEADER
      !ENDDO

      !Search a correct row
 43   CONTINUE
      READ(500+IGR,END=44) TRMODE2
      IF (ABS(TRMODE2-ICYCLE) <= 1E-6) THEN
         GOTO 45
      ELSE
         GOTO 43
      ENDIF
 44      CONTINUE
      WRITE(*,*)''
      WRITE(*,*)' Subroutine TRAJECTORY_FORMAT is not',
     &          ' able to find cycle',ICYCLE
      WRITE(*,*)' from file: TRAJECTORY_'//IGRNUM//'.BIN'
      WRITE(*,*)''
      WRITE(*,*)' Exiting...'
      WRITE(*,*)''
      STOP
 45   CONTINUE
      BACKSPACE(500+IGR)

      !Read data in
      IF(TRMODE == 40)THEN ! Patria's helicopt.
         READ(500+IGR)ICYCLE,XCG,YCG,ZCG,A1,A2,A3,
     &                A4,A5,A6,
     &                OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &                OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,     
     &                A7,A8,A9,A10,A11,A12
         CALL EULFIN(0.0,ROTB1(IGR),ROTA1(IGR),
     &               0.0,0.0,-SHAFT(IGR)/ABS(SHAFT(IGR)),
     &               TX,TY,TZ,2)

         FXTSP(IGR) = (OSKU(IGR)%FXT*TX+OSKU(IGR)%FYT*TY
     &   +OSKU(IGR)%FZT*TZ)

         RTMSP(IGR) = OSKU(IGR)%MZT/FXTSP(IGR)
         FDSP(IGR)  = OSKU(IGR)%MXT/FXTSP(IGR)
      ELSE
         READ(500+IGR)ICYCLE,XCG,YCG,ZCG,
     &                ROTA1(IGR),ROTB1(IGR),CONEA(IGR),
     &                OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     &                OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,     
     &                A1,A2,A3,A4,A5,A6,A11,A12 
         ROTA1(IGR) = ROTA1(IGR)*DEG2RAD
         ROTB1(IGR) = ROTB1(IGR)*DEG2RAD
         CONEA(IGR) = CONEA(IGR)*DEG2RAD
         THETACOL(IGR) = THETACOL(IGR)*DEG2RAD
      ENDIF

      IF (IFA(IGR) == 12) THEN ! Patria's helicopt. 
         THETACOL(IGR) = A4*DEG2RAD
      ELSEIF (IFA(IGR) == 13) THEN !Patria's helicopt. 
         THETACOL(IGR) = A4*DEG2RAD
         THCYCLON(IGR) = A5*DEG2RAD
         THCYCLAT(IGR) = A6*DEG2RAD
      ELSEIF (IFA(IGR) == 14) THEN ! Patria's helicopt.           
         ROTA1(IGR)    = A1*DEG2RAD
         ROTB1(IGR)    = A2*DEG2RAD
         CONEA(IGR)    = A3*DEG2RAD
         THETACOL(IGR) = A4*DEG2RAD
         THCYCLON(IGR) = A5*DEG2RAD
         THCYCLAT(IGR) = A6*DEG2RAD
      ELSEIF(IFA(IGR) > 7 .AND. IFA(IGR) < 11) THEN
         THRUST(IGR) = A1
         TORQUE(IGR) = A2
         ADV(IGR)    = A3
         UTSP(IGR)   = A4
         SHAFT(IGR)  = A5
         FDSP(IGR)   = A6
         OSKU(IGR)%DAMPN = A11
         OSKU(IGR)%DAMPT = A12
      ENDIF
      END SUBROUTINE TRAJECTORY_FORMAT_ACTUATORDISC2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_CONVERTER(IGR,TRMODE,TDBL,DTB)

      USE CONSTANTS      , ONLY : PII, DEG2RAD, RAD2DEG
      USE NS3CO          , ONLY : ROTORW,TIMEL

      INTEGER :: IGR,IHEADER,TRMODE,CYCLES
      REAL    :: PSIR,THETAR,PHIR,PSIM,PSIRH,THETARH,PHIRH
      REAL    :: TIME,XCG,YCG,ZCG,VX,VY,VZ,P,Q,R,CX,CY,CZ,
     &           MXT,MYT,MZT,CMX,CMY,CMZ,FXE,FYE,FZE,MXE,MYE,MZE,
     &           FXT,FYT,FZT,FXH,FYH,FZH,MXH,MYH,MZH
      REAL    :: U,V,W,CXE,CYE,CZE,CMXE,CMYE,CMZE,FXEE,FYEE,FZEE,
     &           MXEE,MYEE,MZEE,FXTE,FYTE,FZTE,MXTE,MYTE,MZTE,
     &           FXHE,FYHE,FZHE,MXHE,MYHE,MZHE,ROTX,ROTY,ROTZ,
     &           DOTZE,DOTBETA,DOTTHETA,TDBL,DTB
      CHARACTER(LEN=307) :: HEADER
      CHARACTER(LEN=3)   :: IGRNUM
      LOGICAL :: THERE

C ... This routine works only external store separation and Helicopter case
      IF (TRMODE == 6 .OR. TRMODE == 61 .OR. TRMODE == 62) THEN
         GOTO 100
      ELSEIF (TRMODE >= 31 .AND. TRMODE <= 39) THEN
         GOTO 400
      ELSE
         GOTO 1000
      ENDIF

C ... Continue EXTERNAL STORE SEPARATION
 100     CONTINUE
C ... Create IGR index
         CALL NUMCH3(IGRNUM,IGR)
C ... Open TRAJECTORY_pb file and write headerlines
         OPEN(500+IGR, FILE='TRAJECTORY_pb_'//IGRNUM//'.dat',
     &         STATUS='UNKNOWN', FORM='FORMATTED') 
         WRITE(500+IGR,120)
 120      FORMAT(' '/,
     &     ' # x,y,z are locations of the Flying Object in the',
     &     ' FINFLO coordinate system.'/,
     &     ' # PSII,THETA,PHI are Euler angles of the',
     &     ' Flying Object, i.e. angles between the fixed coordinate',
     &     ' system (based on the FINFLO coordinate system)',
     &     ' and the particle body coordinate system.'/,
     &     ' # u,v,w are relative velocities of',
     &     ' the Flying Object compared to the A/C in the',
     &     ' particle body coordinate system.', 
     &     ' P,Q,R are angular velocities in the particle',
     &     ' body coordinate system.'/,
     &     ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &     ' moments in the particle body coordinate system.'
     &     ' Fxe,Fye,Fze,Mxe,Mye,Mze are ejector forces and moments'
     &     ' in the particle body coordinate system.'/,
     &     ' # Fxt,Fyt,Fzt,Mxt,Myt,Mzt are thrust forces and moments'
     &     ' in the particle body coordinate system.'/,/,
     &     ' #   time          x               y               z',
     &     '               PSII           THETA            PHI',
     &     '             u               v               w',
     &     '               P               Q               R',
     &     '            Fx            Fy            Fz',
     &     '             Mx            My            Mz',
     &     '            Fxe           Fye           Fze',
     &     '           Mxe           Mye           Mze',
     &     '           Fxt           Fyt           Fzt',
     &     '          Mxt           Myt           Mzt',/
     &     ' #   [s]          [m]             [m]             [m]',
     &     '             [deg]           [deg]           [deg]',
     &     '           [m/s]           [m/s]           [m/s]',
     &     '         [rad/s]         [rad/s]         [rad/s]',
     &     '         [N]           [N]           [N]',
     &     '           [Nm]          [Nm]          [Nm]',
     &     '           [N]           [N]           [N]',
     &     '           [Nm]          [Nm]          [Nm]',
     &     '          [N]           [N]           [N]',
     &     '          [Nm]          [Nm]          [Nm]')

C ... Open TRAJECTORY file
      INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.dat',EXIST=THERE)
      IF(THERE) THEN          
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &         STATUS='UNKNOWN', FORM='FORMATTED')  
      ELSE 
         WRITE(*,*)' '
         WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.dat',
     &                ' was not found'
         WRITE(*,*)' in TRAJECTORY_CONVERT Subroutine.' 
         WRITE(*,*)' '
         WRITE(*,*)' Exiting...'
         WRITE(*,*)' '
         STOP
      ENDIF
C ... Skip headerlines (TRAJECTORY file)
      DO IHEADER=1,16
         READ(400+IGR,'(A307)') HEADER
      ENDDO
C ... Continue by reading data (TRAJECTORY file) 
 200  CONTINUE
C ... Read data (TRAJECTORY file)
      READ(400+IGR,*,END=800)TRMODE,TIME,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,
     &     VX,VY,VZ,P,Q,R,
     &     CX,CY,CZ,CMX,CMY,CMZ,
     &     FXE,FYE,FZE,MXE,MYE,MZE, 
     &     FXT,FYT,FZT,MXT,MYT,MZT 
C ... Attitude angles in radians
      PSIR   = PSIR*DEG2RAD
      THETAR = THETAR*DEG2RAD
      PHIR   = PHIR*DEG2RAD
C ... Data from the FINFLO coord to the particle coord
            CALL FINEUL(PSIR,THETAR,PHIR,VX,VY,VZ,U,V,W,1)
            CALL FINEUL(PSIR,THETAR,PHIR,CX,CY,CZ,CXE,CYE,CZE,1)
            CALL FINEUL(PSIR,THETAR,PHIR,CMX,CMY,CMZ,CMXE,CMYE,CMZE,1)
            CALL FINEUL(PSIR,THETAR,PHIR,FXE,FYE,FZE,FXEE,FYEE,FZEE,1)
            CALL FINEUL(PSIR,THETAR,PHIR,MXE,MYE,MZE,MXEE,MYEE,MZEE,1)
            CALL FINEUL(PSIR,THETAR,PHIR,FXT,FYT,FZT,FXTE,FYTE,FZTE,1)
            CALL FINEUL(PSIR,THETAR,PHIR,MXT,MYT,MZT,MXTE,MYTE,MZTE,1)
C ... Write particle body data to the (TRAJECTORY_pb file)
            WRITE(500+IGR,300)TIME,XCG,YCG,ZCG,
     &           PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &           U,V,W,P,Q,R,CXE,CYE,CZE,CMXE,CMYE,CMZE,FXEE,FYEE,FZEE,
     &           MXEE,MYEE,MZEE,FXTE,FYTE,FZTE,MXTE,MYTE,MZTE
 300        FORMAT(1X,F9.5,12E16.8,18E14.6)
C ... Return reading phase or goto 800
            IF ((ABS(TIME-TDBL))/DTB <= 0.33333) THEN
               GOTO 800
            ELSE
               GOTO 200
            ENDIF
         
C ... Continue HELICOPTER CASE
 400  CONTINUE
C ... Create IGR index
      CALL NUMCH3(IGRNUM,IGR)
C ... Open TRAJECTORY_pb file and write headerlines
      OPEN(500+IGR, FILE='TRAJECTORY_pb_'//IGRNUM//'.dat',
     &     STATUS='UNKNOWN', FORM='FORMATTED') 
      IF (TIMEL) THEN
         WRITE(500+IGR,420)
      ELSE
         WRITE(500+IGR,430)
      ENDIF
 420  FORMAT(' '/,
     &     ' # SHAFT is shaft angle in the helicopter coordinate',
     &     ' system. x,y,z are locations of the rotation point in the',
     &     ' FINFLO coordinate system.'/,
     &     ' # PSII,THETA,BETA are attitude angles of the blade i.e.',
     &     ' angles between the fixed coordinate system',
     &     ' (=FINFLO coordinate system) and the particle body',
     &     ' coordinate system (Rotation order: PSII,BETA,THETA).'/,
     &     ' # ksi,beta,theta are attitude angles of the',
     &     ' blade in the helicopter coordinate system.'/,
     &     ' # u,v,w are velocities of',
     &     ' the blade in the rotating coordinate system.', 
     &     ' zeta_dot,theta_dot,beta_dot are angular velocities',
     &     ' in the helicopter coordinate system.'/,
     &     ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &     ' moments in the helicopter coordinate system.'
     &     ' Fx_hin,Fy_hin,Fz_hin,Mx_hin,My_hin,Mz_hin are'
     &     ' spring and damper forces'
     &     ' and moments in the helicopter coordinate system.'/,
     &     ' # Fx_cen,Fy_cen,Fz_cen,Mx_cen,My_cen,Mz_cen are'
     &     ' centrifugal forces and moments in the  helicopter'
     &     ' coordinate system.'/,/,
     &     ' # time       SHAFT         x           y           z',
     &     '         PSII       THETA        PHI',
     &     '          zeta       beta        theta',
     &     '         u          v           w ',
     &     '       zeta_dot   theta_dot    beta_dot'
     &     '      Fx          Fy          Fz',
     &     '           Mx          My          Mz',
     &     '        Fx_hin      Fy_hin      Fz_hin',
     &     '      Mx_hin      My_hin      Mz_hin',
     &     '     Fx_cen      Fy_cen      Fz_cen',
     &     '       Mx_cen      My_cen      Mz_cen',/
     &     ' # [s]        [deg]        [m]         [m]         [m]',
     &     '       [deg]       [deg]       [deg]',
     &     '        [deg]       [deg]       [deg]',
     &     '        [m/s]      [m/s]       [m/s]',
     &     '      [rad/s]     [rad/s]     [rad/s]',
     &     '      [N]         [N]         [N]',
     &     '          [Nm]        [Nm]        [Nm]',
     &     '         [N]         [N]         [N]',
     &     '         [Nm]        [Nm]        [Nm]',
     &     '        [N]         [N]         [N]',
     &     '         [Nm]        [Nm]        [Nm]')
 430  FORMAT(' '/,
     &     ' # SHAFT is shaft angle in the helicopter coordinate',
     &     ' system. x,y,z are locations of the rotation point in the',
     &     ' FINFLO coordinate system.'/,
     &     ' # PSII,THETA,BETA are attitude angles of the blade i.e.',
     &     ' angles between the fixed coordinate system',
     &     ' (=FINFLO coordinate system) and the particle body',
     &     ' coordinate system (Rotation order: PSII,BETA,THETA).'/,
     &     ' # ksi,beta,theta are attitude angles of the',
     &     ' blade in the helicopter coordinate system.'/,
     &     ' # u,v,w are velocities of',
     &     ' the blade in the rotating coordinate system.', 
     &     ' zeta_dot,theta_dot,beta_dot are angular velocities',
     &     ' in the helicopter coordinate system.'/,
     &     ' # Fx,Fy,Fz,Mx,My,Mz are aerodynamic forces and',
     &     ' moments in the helicopter coordinate system.'
     &     ' Fx_hin,Fy_hin,Fz_hin,Mx_hin,My_hin,Mz_hin are'
     &     ' spring and damper forces'
     &     ' and moments in the helicopter coordinate system.'/,
     &     ' # Fx_cen,Fy_cen,Fz_cen,Mx_cen,My_cen,Mz_cen are'
     &     ' centrifugal forces and moments in the helicopter'
     &     ' coordinate system.'/,/,
     &     ' # Cycles     SHAFT         x           y           z',
     &     '         PSII       THETA        PHI',
     &     '          zeta       beta        theta',
     &     '         u          v           w ',
     &     '       zeta_dot   theta_dot    beta_dot'
     &     '      Fx          Fy          Fz',
     &     '           Mx          My          Mz',
     &     '        Fx_hin      Fy_hin      Fz_hin',
     &     '      Mx_hin      My_hin      Mz_hin',
     &     '     Fx_cen      Fy_cen      Fz_cen',
     &     '       Mx_cen      My_cen      Mz_cen',/
     &     ' #            [deg]        [m]         [m]         [m]',
     &     '       [deg]       [deg]       [deg]',
     &     '        [deg]       [deg]       [deg]',
     &     '        [m/s]      [m/s]       [m/s]',
     &     '      [rad/s]     [rad/s]     [rad/s]',
     &     '      [N]         [N]         [N]',
     &     '          [Nm]        [Nm]        [Nm]',
     &     '         [N]         [N]         [N]',
     &     '         [Nm]        [Nm]        [Nm]',
     &     '        [N]         [N]         [N]',
     &     '         [Nm]        [Nm]        [Nm]')
C ... Open TRAJECTORY file
      INQUIRE(FILE='TRAJECTORY_'//IGRNUM//'.dat',EXIST=THERE)
      IF(THERE) THEN          
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &         STATUS='UNKNOWN', FORM='FORMATTED')  
      ELSE 
         WRITE(*,*)' '
         WRITE(*,*)' Warning: File TRAJECTORY_'//IGRNUM//'.dat',
     &                ' was not found'
         WRITE(*,*)' in TRAJECTORY_CONVERT Subroutine.' 
         WRITE(*,*)' '
         WRITE(*,*)' Exiting...'
         WRITE(*,*)' '
         STOP
      ENDIF
C ... Skip headerlines (TRAJECTORY file)
      DO IHEADER=1,16
         READ(400+IGR,'(A307)') HEADER
      ENDDO
C ... Continue by reading data (TRAJECTORY file) 
 500  CONTINUE
C ... Read data (TRAJECTORY file)
      IF (.NOT. TIMEL) THEN
         READ(400+IGR,*,END=800)CYCLES,PSIM,XCG,YCG,ZCG,
     &        PSIR,THETAR,PHIR,
     &        VX,VY,VZ,DOTZE,DOTBETA,DOTTHETA,
     &        CX,CY,CZ,CMX,CMY,CMZ,
     &        FXH,FYH,FZH,MXH,MYH,MZH, 
     &        FXE,FYE,FZE,MXE,MYE,MZE 
      ELSEIF (TIMEL) THEN
         READ(400+IGR,*,END=800)TIME,PSIM,XCG,YCG,ZCG,
     &        PSIR,THETAR,PHIR,
     &        VX,VY,VZ,DOTZE,DOTBETA,DOTTHETA,
     &        CX,CY,CZ,CMX,CMY,CMZ,
     &        FXH,FYH,FZH,MXH,MYH,MZH, 
     &        FXE,FYE,FZE,MXE,MYE,MZE
      ENDIF

C ... Attitude angles in radians
      PSIM   = PSIM*DEG2RAD
      PSIR   = PSIR*DEG2RAD
      THETAR = THETAR*DEG2RAD
      PHIR   = PHIR*DEG2RAD
C ... Data from the FINFLO coord to the particle coord
c      CALL FINEUL(PSIR,PHIR,THETAR,VX,VY,VZ,U,V,W,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,CX,CY,CZ,CXE,CYE,CZE,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,CMX,CMY,CMZ,CMXE,CMYE,CMZE,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,FXE,FYE,FZE,FXEE,FYEE,FZEE,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,MXE,MYE,MZE,MXEE,MYEE,MZEE,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,FXH,FYH,FZH,FXHE,FYHE,FZHE,2)
c      CALL FINEUL(PSIR,PHIR,THETAR,MXH,MYH,MZH,MXHE,MYHE,MZHE,2)
C ... Data from the particle coord to the FINFLO coord
      CALL EULFIN(PSIR,PHIR,THETAR,P,Q,R,ROTX,ROTY,ROTZ,2)
C ... Translation velocities to the helicopter system  
      U = -COS(PSIM)*VX-SIN(PSIM)*VZ
      V = SIN(PSIM)*VX-COS(PSIM)*VZ
      W = -VY
C ... Aerodynamic forces and moments to the helicopter system  
      CXE = -COS(PSIR)*CX-SIN(PSIR)*CZ
      CYE = SIN(PSIR)*CX-COS(PSIR)*CZ
      CZE = -CY
      CMXE = -COS(PSIR)*CMX-SIN(PSIR)*CMZ
      CMYE = SIN(PSIR)*CMX-COS(PSIR)*CMZ
      CMZE = -CMY
C ... Hinge forces and moments to the helicopter system  
      FXHE = -COS(PSIR)*FXH-SIN(PSIR)*FZH
      FYHE = SIN(PSIR)*FXH-COS(PSIR)*FZH
      FZHE = -FYH
      MXHE = -COS(PSIR)*MXH-SIN(PSIR)*MZH
      MYHE = SIN(PSIR)*MXH-COS(PSIR)*MZH
      MZHE = -MYH
C ... Centrifugal forces and moments to the helicopter system  
      FXEE = -COS(PSIR)*FXE-SIN(PSIR)*FZE
      FYEE = SIN(PSIR)*FXE-COS(PSIR)*FZE
      FZEE = -FYE
      MXEE = -COS(PSIR)*MXE-SIN(PSIR)*MZE
      MYEE = SIN(PSIR)*MXE-COS(PSIR)*MZE
      MZEE = -MYE
C ... Angular velocities to the helicopter system      
c      DOTBETA  = COS(PSIR)*ROTX+SIN(PSIR)*ROTZ
c      DOTTHETA = SIN(PSIR)*ROTX-COS(PSIR)*ROTZ
c      DOTZE    = ROTORW-ROTY
C ... Angles to the helicotper system
      PSIRH   = PSIR-PSIM
      THETARH = THETAR
      PHIRH   = -PHIR
      SHAFT   = -PSIM+PII/2

C ... Write particle body data to the (TRAJECTORY_pb file)
      IF (.NOT. TIMEL) THEN
         WRITE(500+IGR,600)CYCLES,SHAFT*RAD2DEG,XCG,YCG,ZCG,
     &        PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &        PSIRH*RAD2DEG,PHIRH*RAD2DEG,THETARH*RAD2DEG,
     &        U,V,W,DOTZE,DOTTHETA,DOTBETA,CXE,CYE,CZE,CMXE,CMYE,CMZE,
     &        FXHE,FYHE,FZHE,MXHE,MYHE,MZHE,
     &        FXEE,FYEE,FZEE,MXEE,MYEE,MZEE
 600     FORMAT(1X,I8,F11.4,15F12.4,18E12.4)
      ELSEIF (TIMEL) THEN
         WRITE(500+IGR,601)TIME,SHAFT*RAD2DEG,XCG,YCG,ZCG,
     &        PSIR*RAD2DEG,THETAR*RAD2DEG,PHIR*RAD2DEG,
     &        PSIRH*RAD2DEG,PHIRH*RAD2DEG,THETARH*RAD2DEG,
     &        U,V,W,DOTZE,DOTTHETA,DOTBETA,CXE,CYE,CZE,CMXE,CMYE,CMZE,
     &        FXHE,FYHE,FZHE,MXHE,MYHE,MZHE,
     &        FXEE,FYEE,FZEE,MXEE,MYEE,MZEE
 601     FORMAT(1X,F8.4,F12.4,15F12.4,18E12.4)
      ENDIF
C ... Return reading phase
      GOTO 500

C ... Continue and close channels
 800  CONTINUE
      CLOSE(400+IGR)
      CLOSE(500+IGR)

C ... Continue
 1000 CONTINUE
      
      RETURN
      END SUBROUTINE TRAJECTORY_CONVERTER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRAJECTORY_BACKUP(IGR)

      IMPLICIT NONE

      INTEGER :: IGR
      CHARACTER(LEN=3) :: IGRNUM
      CHARACTER(LEN=1000) :: RIVI

C ... Create IGR index
         CALL NUMCH3(IGRNUM,IGR)
C ... Open TRAJECTORY file         
         OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &         STATUS='UNKNOWN', FORM='FORMATTED')  
C ... Open TRAJECTORY back up file and write headerlines
         OPEN(500+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat.BAK',
     &         STATUS='UNKNOWN', FORM='FORMATTED') 
C ... Continue
 100  CONTINUE
C ... Read and write data
      READ(400+IGR,'(A)',END=200)RIVI
      WRITE(500+IGR,'(A)')TRIM(RIVI)
C ... Return reading phase
      GOTO 100
C ... Continue and close channels
 200  CONTINUE
      CLOSE(400+IGR)
      CLOSE(500+IGR)

      END SUBROUTINE TRAJECTORY_BACKUP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTORPOS(XCO,YCO,ZCO,XGRIG,YGRIG,
     &                   ZGRIG,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,N,T,DT,NGPTS)

      IMPLICIT NONE

      INTEGER :: IBLOCK,IMAX,JMAX,KMAX,IN,JN,KN,N,IG1,IHEADER,I,NGPTS

      REAL, DIMENSION(*) :: XCO,YCO,ZCO,XGRIG,YGRIG,ZGRIG
      REAL :: TIME,OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,T,DT,
     &        ROTANG

      CHARACTER(LEN=235) :: HEADER

      LOGICAL :: THERE

C ... Update an original grid to a new grid
      DO I = 1,NGPTS
         XCO(I) = XGRIG(I)
         YCO(I) = YGRIG(I)
         ZCO(I) = ZGRIG(I)
      ENDDO

cc      CALL TRAJECTORY_FORMAT(T,DT,XCG(IGR),YCG(IGR),ZCG(IGR),
cc     +     PSIR(IGR),THETAR(IGR),PHIR(IGR),
cc     +     OSKU,IGR,TIMEL,TRMODE(IGR),IGRID(A,1))
      
      RETURN
      END SUBROUTINE ROTORPOS
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE ROTOR_FORCES(M,XCG,YCG,ZCG,
     +     ROTX,ROTY,ROTZ,VX,VY,VZ,
     +     CX,CY,CZ,
     +     CMX,CMY,CMZ,
     +     FXE,FYE,FZE,
     +     MXE,MYE,MZE,
     +     FXT,FYT,FZT,
     +     MXT,MYT,MZT,
     +     FXH,FYH,FZH,
     +     MXH,MYH,MZH,
     +     PSI,THETA,PHI,RCG,
     +     IX,IY,IZ,IXY,IXZ,IYZ,
     +     ROTORW,NBLOCK,TRMODE,IGR)

      USE FLIGHT , ONLY :  PSIM
      USE CONSTANTS , ONLY : DEG2RAD
      USE MPI
      USE NS3CO     , ONLY : TIMEL, T, DT, GX, GY, GZ, GROUND
    
      IMPLICIT NONE

      INTEGER :: I,NBLOCK,TRMODE,IGR,IHEADER,IERR,IPRO

      REAL :: M,ROTORW,ROTORF,RCG,RBLADE,
     +     CX,CY,CZ,CMX,CMY,CMZ,
     +     FXE,FYE,FZE,MXE,MYE,MZE,
     +     FXEE,FYEE,FZEE,MXEE,MYEE,MZEE,
     +     FXTE,FYTE,FZTE,MXTE,MYTE,MZTE,FXT,FYT,FZT,MXT,MYT,MZT,
     +     FXGE,FYGE,FZGE,MXGE,MYGE,MZGE,FXG,FYG,FZG,MXG,MYG,MZG,
     +     XCGBI,ZCGBI,XCGHI,ZCGHI,
     +     FXHE,FYHE,FZHE,MXHE,MYHE,MZHE,FXH,FYH,FZH,MXH,MYH,MZH,
     +     ROTX,ROTY,ROTZ,VX,VY,VZ,
     +     TIME
      REAL :: PSI,THETA,PHI,ROTU,ROTV,ROTW,ROTP,ROTQ,ROTR,
     +        CXE,CYE,CZE,CMXE,CMYE,CMZE,
     +        GXX,GYY,GZZ,FX,FY,FZ,MX,MY,MZ,
     +        XCG,YCG,ZCG,XCGB,YCGB,ZCGB
      REAL :: PDOT,QDOT,RDOT,IX,IY,IZ,IXY,IXZ,IYZ,QADD
      REAL :: THETACOL,THCYCLAT,THCYCLON,KDELTA3

      CHARACTER(LEN=235) :: HEADER
      CHARACTER(LEN=128) :: VALI
      CHARACTER(LEN=3)   :: IGRNUM

      LOGICAL :: PARALLEL


C ... This subprogram is performed only in the 6-dof helicopter calculation
      IF (TRMODE /= 31) THEN
         RETURN                
      ENDIF

C ... Velocities of the blade from the FINFLO coord to the particle body coord
      CALL FINEUL(PSI,PHI,THETA,VX,VY,VZ,ROTU,ROTV,ROTW,2)
      CALL FINEUL(PSI,PHI,THETA,ROTX,ROTY,ROTZ,ROTP,ROTQ,ROTR,2)
C ... Forces and moments of the blade from the FINFLO coord 
C ... to the particle body coord
      CALL FINEUL(PSI,PHI,THETA,CX,CY,CZ,CXE,CYE,CZE,2)
      CALL FINEUL(PSI,PHI,THETA,CMX,CMY,CMZ,CMXE,CMYE,CMZE,2)  

C ... Angular accelerations
C ... If time is zero
      IF (T/DT <= 0.33333) THEN
         PDOT = 0.
         QDOT = 0.
         RDOT = 0.
         GOTO 400
      ELSE
         CONTINUE
      ENDIF
C ... These procedures are needed if time is not zero
C ... Open "TRAJECTORY_IGR" file 
      CALL NUMCH3(IGRNUM,IGR)
      OPEN(400+IGR, FILE='TRAJECTORY_'//IGRNUM//'.dat',
     &     STATUS='UNKNOWN', FORM='FORMATTED') 
C ... Skip headerlines
         DO IHEADER=1,16
            READ(400+IGR,'(A235)') HEADER
         ENDDO
C ... Search a correct row (=correct timelevel)
 100     CONTINUE
         READ(400+IGR,'(1X,F8.4)') TIME
         IF ((ABS(TIME-T+DT))/DT <= 0.33333) THEN
            GOTO 300               
         ELSE 
            GOTO 100               
         ENDIF
 300     CONTINUE
         BACKSPACE(400+IGR)
         READ(400+IGR,'(128X,3F12.4)') PDOT,QDOT,RDOT ! current timelevel
         PDOT = (ROTP-PDOT)/DT ! accelerations
         QDOT = (ROTQ-QDOT)/DT ! accelerations
         RDOT = (ROTR-RDOT)/DT ! accelerations         
         CLOSE(400+IGR)
 400     CONTINUE

C ... Gravity forces and moments in the particle body coordinate system
      CALL FINEUL(PSI,PHI,THETA,0.0,GY,0.0,GXX,GYY,GZZ,2)
      FXTE  = GXX*M
      FYTE  = GYY*M
      FZTE  = GZZ*M
      MXTE  = 0.
      MYTE  = 0. 
      MZTE  = 0.
C ... Gravity data to the FINFLO coordinate system
      CALL EULFIN(PSI,PHI,THETA,FXTE,FYTE,FZTE,FX,FY,FZ,2)
      CALL EULFIN(PSI,PHI,THETA,MXTE,MYTE,MZTE,
     +     MX,MY,MZ,2)
      FXT   = FX
      FYT   = FY
      FZT   = FZ
      MXT   = MX
      MYT   = MY
      MZT   = MZ

C ... Centrifugal forces and moments in the particle body coordinate system
      FXEE = (-ROTQ*(ROTP*RCG))*M
     &     +(-ROTQ*ROTW+ROTR*ROTV)*M
      FYEE = (-ROTR*(-ROTR*RCG)
     &     +ROTP*(ROTP*RCG))*M+(-ROTR*ROTU+ROTP*ROTW)*M
      FZEE = (ROTQ*(-ROTR*RCG))*M+(-ROTP*ROTV+ROTQ*ROTU)*M
      MXEE = (-ROTP*ROTV+ROTQ*ROTU)*M*RCG
     &     +(-IXY*(ROTP*ROTR-QDOT)+IXZ*(ROTP*ROTQ+RDOT)
     &     +IYZ*(ROTQ**2-ROTR**2)-ROTQ*ROTR*(IZ-IY))
      MYEE = 0.+(IXY*(ROTQ*ROTR+PDOT)-IXZ*(ROTP**2-ROTR**2)
     &       -IYZ*(ROTP*ROTQ-RDOT)-ROTP*ROTR*(IX-IZ))
      MZEE = -(-ROTQ*ROTW+ROTR*ROTV)*M*RCG
     &     +(-IXY*(ROTQ**2-ROTP**2)-IXZ*(ROTQ*ROTR-PDOT)
     &     +IYZ*(ROTP*ROTR+QDOT)-ROTP*ROTQ*(IY-IX))
C ... Centrifugal data to the FINFLO coordinate system
      CALL EULFIN(PSI,PHI,THETA,FXEE,FYEE,FZEE,FX,FY,FZ,2)
      CALL EULFIN(PSI,PHI,THETA,MXEE,MYEE,MZEE,
     +     MX,MY,MZ,2)
      FXE = FX
      FYE = FY
      FZE = FZ
      MXE = MX
      MYE = MY
      MZE = MZ

C ... Hinge forces in the particle body coordinate system
      FXHE = -CXE-FXTE-FXEE
      FYHE = -CYE-FYTE-FYEE
      FZHE = -CZE-FZTE-FZEE
C ... These values must be changed to the input data
      THETACOL = 4.0
      THCYCLAT = 0.0
      THCYCLON = 0.0
      KDELTA3  = 0.0
      THETACOL = DEG2RAD*THETACOL
      THCYCLON = DEG2RAD*THCYCLON
      THCYCLAT = DEG2RAD*THCYCLAT
C ... This angular velocity is used to stabilize theta angle      
      QADD = -ROTQ+ROTR*TAN(PHI)
     &     -((THETA- 
     &     THETACOL+THCYCLAT*SIN(PSIM(IGR))
     &     +THCYCLON*COS(PSIM(IGR))-KDELTA3*(-PHI))
     &     /DT)/COS(PHI)          
C ... Hinge moments in the particle body coordinate system
      MYHE = -MYTE-CMYE-MYEE+(QADD/DT)*IY 
      MXHE = -3600.0*PHI            ! elastomeerin momentti SUBROUTINE HINGE_MOMENT(L) would be better here
      MZHE = 3600.0*(PSIM(IGR)-PSI) ! elastomeerin momentti SUBROUTINE HINGE_MOMENT(L) would be better here
C ... Hinge data to the FINFLO coordinate system
      CALL EULFIN(PSI,PHI,THETA,FXHE,FYHE,FZHE,FX,FY,FZ,2)
      CALL EULFIN(PSI,PHI,THETA,MXHE,MYHE,MZHE,MX,MY,MZ,2)
      FXH = FX
      FYH = FY
      FZH = FZ
      MXH = MX
      MYH = MY
      MZH = MZ
      MYH = (MYH - 56800.0*(PSIM(IGR)-PSI)-1170*(ROTY-ROTORW)) ! sylinterien momentti SUBROUTINE HINGE_MOMENT(L) would be better here 
      RETURN
      END SUBROUTINE ROTOR_FORCES
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE BLADE_CYCLES(ICYCLEG,TIMEG,NGRIFL)

C ... Printing routine for the 2-dof model
C ... Similar to the emulator

      USE BLADE_CONSTANTS
      USE BLADE_VARIABLES
      USE FLIGHT,   ONLY : PSIM,PHIR,PSIR,THETAR,XCG,YCG,ZCG,OSKU,TRMODE
      USE NS3CO,    ONLY : TIMEL,ROTORW
      USE INTEGERS, ONLY : IREPEA,NREPEA
      USE CONSTANTS,ONLY : PII

      IMPLICIT NONE

      INTEGER :: ICYCLEG,L,ICMAX,NGRIFL
      REAL    :: TIMEG,XCGGG,ZCGGG,APU1

C ... Initial values
 
c      ICYCLE = 0
c      TIME   = -DT
c      IPR    = 1

1000  CONTINUE


      ICMAX = DTB/DTB



      IF(ICYCLEG == 1) THEN
         WRITE(4,'(/1X,A,I5,2A)')'In BLADES ',ICMAX,' iteration cycles',
     &   ' are made / flow solution step'
      ENDIF

      IF(ICMAX < 1 .AND. IREPEA(6) < NREPEA(6)) THEN
         WRITE(13,'(2A/2A,I1)')'  Rotor time step is larger than the ',
     +          ' flow solution step.','  The latter is being applied',
     +          ' for BLADE ',L  
      ENDIF

C ... Transfer flight mechanic data to the helicopter variables
      DO L = 1,NGRIFL
      BETAH(L)   = -PHIR(L)
      ZETAH(L)   = PSIR(L) - PSIM(L)
      THETA(L)   = THETAR(L)
      SHAFT(L)   = -PSIM(L)+PII/2
C      DOTBETA(L) = COS(PSIR(L))*OSKU(L)%ROTX+SIN(PSIR(L))*OSKU(L)%ROTZ
C      DOTZE(L)   = ROTORW-OSKU(L)%ROTY
      Y(1,L)     = ZETAH(L)
      Y(2,L)     = BETAH(L)
      Y(3,L)     = DOTZE(L)
      Y(4,L)     = DOTBETA(L)
      YOLD       = Y 
      ENDDO

      ICMAX = MAX(1,ICMAX)
      DO ICYCLE = 1,ICMAX ! Internal steps    ! Tarvitaanko tata do-luuppia?
      DO L = 1,NGRIFL

      TIMEB(L)   = TIMEB(L) + DTB  ! onko naa jo vanhoja juttuja?
      ICYCLEB(L) = ICYCLEB(L) + 1 ! onko naa jo vanhoja juttuja?

C ... Solve the equations at time + dt
         
c      IF(IORD <= 2) THEN
c         CALL SOLVER
c      ELSE IF(IORD == 4) THEN

C ... New rotor angle in time accurate scheme
       IF(TIMEL .AND. TDBL/DTB >= 0.33333) THEN
          PSIM(L) = PSIM(L) - DTB*ROTORW          
          SHAFT(L)= -PSIM(L)+PII/2
       ENDIF

C ... Call Blade Solver routine
      CALL BLADE_SOLVER_RK4(L)
c      ELSE
c         STOP 'Solution option is not available'
c      ENDIF

cc_jil 15.3.2010      XCG(L)    = COS(PSIM(L))*RNIVEL(L)
cc_jil 15.3.2010      YCG(L)    = 0.
cc_jil 15.3.2010      ZCG(L)    =-SIN(PSIM(L))*RNIVEL(L)




c         beta(l) = 0.; zeta(L) = 0. ! paikallaan
c      ZETAH(L) = -10.

c      write(800+L,*) TIMEB(L),L,xcg(l),ycg(l),zcg(l),psim(l),rnivel(l),
c     & omega,zeta(l),beta(l)
c      write(800+l,*) y(1,L),y(2,l),y(3,l),y(4,l),yold(1,l),yold(2,l),
c     & yold(3,l),yold(4,l)

c       write(601,*) l,iprb(l),icycleb(l),icycleg
ccc      IF(IPRB(L) == ICYCLEB(L) .AND..NOT.TIMEL) THEN ! PRINT
ccc         CALL BLADE_OUTPUT(ICYCLEG,TIMEG,L)
ccc         IPRB(L) = IPRB(L) + IPRINT ! Printing interval
ccc      ENDIF

cc_jil 15.3.2010    PSIR(L) = ZETAH(L) - PSIM(L) ! Negative in flight mechanics

C ... Angles to the helicopter coordinate system
      IF (TRMODE(L) == 37) THEN
         ZETAH(L)   = ZETAHCOL(L)
         BETAH(L)   = BETAHCOL(L) - BECYCLAT(L)*SIN(SHAFT(L))
     &        - BECYCLON(L)*COS(SHAFT(L))
         DOTBETA(L)= ((BETAHCOL(L) -BECYCLAT(L)*SIN(SHAFT(L)+DTB*ROTORW)
     &        - BECYCLON(L)*COS(SHAFT(L)+DTB*ROTORW)) - BETAH(L))/DTB
         Y(1,L)     = ZETAH(L)
         Y(2,L)     = BETAH(L)
         Y(3,L)     = DOTZE(L)
         Y(4,L)     = DOTBETA(L)
         CALL HINGE_MOMENT(L) ! updates hinge moments
         CALL AERODYN_MOMENT(L) ! updates aerodyn moments
         CALL BLADE_DERIV(L) ! updates centrifugal moments
      ELSEIF (TRMODE(L) == 38) THEN
         ZETAH(L)  = MIN(.25*PII,Y(1,L));ZETAH(L)=MAX(-.25*PII,ZETAH(L))
         BETAH(L)  = Y(2,L)
         DOTZE(L)  = Y(3,L)  
         DOTBETA(L)= Y(4,L)
      ENDIF

C ... Angles and dottheta to the flight mechanics system
      PSIR(L)      = ZETAH(L) + PSIM(L)
      PHIR(L)      = -BETAH(L) 
      THETA(L)     = THETACOL(L) - THCYCLAT(L)*SIN(SHAFT(L))
     &     - THCYCLON(L)*COS(SHAFT(L))-KDELTA3(L)*BETAH(L)
      DOTTHETA(L)  = ((THETACOL(L) -THCYCLAT(L)*SIN(SHAFT(L)+DTB*ROTORW)
     &     - THCYCLON(L)*COS(SHAFT(L)+DTB*ROTORW) - KDELTA3(L)*BETAH(L)) 
     &     - THETA(L))/DTB
      THETAR(L)    = THETA(L)

C ... Blade position, translation and angular velocity to the FINFLO coord.
      XCG(L)     = SIN(PSIM(L))*RHINGE(L)
      YCG(L)     = 0.
      ZCG(L)     = -COS(PSIM(L))*RHINGE(L)
      OSKU(L)%VX = -COS(PSIM(L))*RHINGE(L)*ROTORW
      OSKU(L)%VY = 0.
      OSKU(L)%VZ = -SIN(PSIM(L))*RHINGE(L)*ROTORW
      OSKU(L)%ROTX = COS(PSIR(L))*DOTBETA(L)+SIN(PSIR(L))*DOTTHETA(L)
      OSKU(L)%ROTY = -DOTZE(L)+ROTORW
      OSKU(L)%ROTZ = SIN(PSIR(L))*DOTBETA(L)-COS(PSIR(L))*DOTTHETA(L)


C ... paikallaan edellä!

      ENDDO ; ENDDO

ccc      IF(TIMEL) THEN
ccc      DO L = 1,NGRIFL
ccc         CALL BLADE_OUTPUT(ICYCLEG,TIMEG,L)
ccc      ENDDO
ccc      ENDIF ! TIMEL

      END SUBROUTINE BLADE_CYCLES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BLADE_OUTPUT(ICYCLEG,TIMEG,L)

C ... Printing routine for the 2-dof model
C ... Similar to the emulator

      USE BLADE_CONSTANTS
      USE BLADE_VARIABLES
      USE FLIGHT, ONLY : PSIM,PSIR,THETAR,PHIR,OSKU
      USE NS3CO,    ONLY : TIMEL,ROTORW

      IMPLICIT NONE

      INTEGER :: ICYCLEG,L,N,IC,ID
      REAL    :: ANGLE,RXO,RYO,RZO,RX,RY,RZ,APU1,APU2,BETAC,ZETAC
      REAL    :: TIMEG

C ... Initial values
 
c      ICYCLE = 0
c      TIME   = -DT
c      IPR    = 1

c1000  CONTINUE

c      TIME    = TIME + DT
c      ICYCLE  = ICYCLE + 1
c      YOLD = Y
c        write(6,*) 'print', time,timeg,icycle,icycleg,ipr,iprint
C ... Solve the equations at time + dt
         
c      IF(IORD <= 2) THEN
c         CALL SOLVER
c      ELSE IF(IORD == 4) THEN
c         CALL SOLVER_RK4
c      ELSE
c         STOP 'Solution option is not available'
c      ENDIF
        
c      PSI     = PSI + DT*OMEGA

c      DO L = 1,NGRIFL
c        write(6,*) icycleg,timeg,l
C ... Absolute blade angle, blades 2-L lead blade 1: 
cc_jil 15.3.2010      ANGLE     = PSIM(L)! + (L-1)*SHIFT
       ANGLE     = -PSIM(L)! + (L-1)*SHIFT
c      PSIM(L)   = ANGLE
      
C ... Blade pitch angle control
      THETA(L)    = THETACOL(L) - THCYCLAT(L)*SIN(ANGLE) 
     &      - THCYCLON(L)*COS(ANGLE) - KDELTA3(L)*BETAH(L)
      THETAR(L)   = THETA(L)

c      IF(IPR == ICYCLE) THEN ! PRINT
c      TIME   = TIMEB(L)

      IF(.NOT. TIMEL) THEN
         TIME = ICYCLEB(L)
      ELSE
         TIME = TIMEG
      ENDIF

      IC     = 1000*L+13                                   ! File OUTDER
      ID     = 1000*L+113
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,9(1X,D13.6))") TIME,ANGLE*RAD2DEG,(F(N,L),N=1,4)
      WRITE(ID,"(F8.5,9(1X,D13.6))") TIME,ANGLE*RAD2DEG,(F(N,L),N=1,4)
      ELSE
      WRITE(IC,"(F10.0,9(1X,D13.6))") TIME,ANGLE*RAD2DEG,(F(N,L),N=1,4)
      ENDIF
      IC     = 1000*L+12                                   ! File OUTPUT
      ID     = 1000*L+112
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,6(1X,D13.6))") TIME,ANGLE*RAD2DEG,
     &  (Y(N,L)*RAD2DEG,N=1,2),(Y(N,L),N=3,4),THETA(L)*RAD2DEG
      WRITE(ID,"(F8.5,6(1X,D13.6))") TIME,ANGLE*RAD2DEG,
     &  (Y(N,L)*RAD2DEG,N=1,2),(Y(N,L),N=3,4),THETA(L)*RAD2DEG
      ELSE
      WRITE(IC,"(F10.0,6(1X,D13.6))") TIME,ANGLE*RAD2DEG,
     &  (Y(N,L)*RAD2DEG,N=1,2),(Y(N,L),N=3,4),THETA(L)*RAD2DEG
      ENDIF
      RXO   = 0.   ; RYO  = 1.  ; RZO  = 0. ; APU1 = THETA(L) !0. ;
      CALL EULFIND(ZETAH(L)-ANGLE,-BETAH(L),THETA(L),RXO,RYO,RZO,
     &   RX,RY,RZ,2)
cc      write(532,*) l,osku(l)%cx,osku(l)%cy,osku(l)%cz
c      call eulfind(ZETAH(L)-ANGLE,-BETAH(L),THETA(L),
c     & 1.0,0.0,0.0,rx,ry,rz,2)
c      write(701,*) l,rad2deg*(zeta(l)-angle),rx,ry,rz

      CALL EULFIND(ZETAH(L)-ANGLE,-BETAH(L),THETA(L),RXO,RYO,RZO,
     &   RX,RY,RZ,2)

c       write(900+l,*)ZETAH(L)-ANGLE,BETAH(L),THETA(L),
c     & RXO,RYO,RZO,RX,RY,RZ
c      write(464+l,*)zeta(l),angle,beta(l),theta(l),apu1,rxo,ryo,
c     & rzo,rx,ry,rz
C ... ZETAH(L)-ANGLE,BETAH(L) correspond to Phi,Theta in AC flight mechanics
      BETAC = ASIN(-RZ)*RAD2DEG   ! Beta from the position vector
      APU2  = 1.
      APU1  = MIN(APU2,-SIN(ANGLE)*RX - COS(ANGLE)*RY)
      APU2  =-1.
      APU1  = MAX(APU1,APU2)
      ZETAC = -ASIN(APU1)*RAD2DEG  ! ZETA from the position vector
      IC     = 1000*L+21                                   ! File OUTPOSI  
      ID     = 1000*L+121                                  ! File OUTPOSI  
      IF(TIMEL) THEN     
      WRITE(IC,"(F8.5,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      WRITE(ID,"(F8.5,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ELSE
      WRITE(IC,"(F10.0,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ENDIF
      RXO  = 0.   ; RYO  = 1.     ; RZO  = 0. ; APU1 = THETA(L) !0. ;
      CALL EULFIND(ZETAH(L)-ANGLE,THETA(L),-BETAH(L),RXO,RYO,RZO,
     &   RX,RY,RZ,2)
      IC     = 1000*L+22                                   ! File OUTPOSJ  
      ID     = 1000*L+122                                  ! File OUTPOSJ  
      IF(TIMEL) THEN    
      WRITE(IC,"(F8.5,1X,4(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      WRITE(ID,"(F8.5,1X,4(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ELSE
      WRITE(IC,"(F10.0,1X,4(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ENDIF
      RXO  = 0.   ; RYO  = 0.     ; RZO  = 1. ;  APU1 = THETA(L) !0. ;
      CALL EULFIND(ZETAH(L)-ANGLE,THETA(L),-BETAH(L),RXO,RYO,RZO,
     &   RX,RY,RZ,2)
      IC     = 1000*L+23                                   ! File OUTPOSK
      IC     = 1000*L+123                                  ! File OUTPOSK
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      WRITE(ID,"(F8.5,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ELSE
      WRITE(IC,"(F10.0,4(1X,D13.6))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & RX,RY,RZ
      ENDIF
      IC     = 1000*L+14                                    ! File OUTFOR
      ID     = 1000*L+114                                   ! File OUTFOR
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,1X,9(D13.6,1X))")TIME,ZETAH(L)*RAD2DEG,FHINGEX(L),
     &   FHINGEY(L),FHINGEI(L),FHINGEJ(L),FHINGEK(L),RMHINGEI(L),
     &   RMHINGEJ(L),RMHINGEK(L)
      WRITE(ID,"(F8.5,1X,9(D13.6,1X))")TIME,ZETAH(L)*RAD2DEG,FHINGEX(L),
     &   FHINGEY(L),FHINGEI(L),FHINGEJ(L),FHINGEK(L),RMHINGEI(L),
     &   RMHINGEJ(L),RMHINGEK(L)
      ELSE
      WRITE(IC,"(F10.0,9(1X,D13.6))") TIME,ZETAH(L)*RAD2DEG,FHINGEX(L),
     &   FHINGEY(L),FHINGEI(L),FHINGEJ(L),FHINGEK(L),RMHINGEI(L),
     &   RMHINGEJ(L),RMHINGEK(L)
      ENDIF
      IC     = 1000*L+15                                    ! File OUTDAM
      ID     = 1000*L+115                                   ! File OUTDAM
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QZE_S1(L),
     &   QZE_S2(L),QZE_D(L),QBE_S1(L),QBE_S2(L),QBE_D(L),QZE(L),QBE(L),
     &   QZE_A(L),QBE_A(L)
      WRITE(ID,"(F8.5,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QZE_S1(L),
     &   QZE_S2(L),QZE_D(L),QBE_S1(L),QBE_S2(L),QBE_D(L),QZE(L),QBE(L),
     &   QZE_A(L),QBE_A(L)
      ELSE
      WRITE(IC,"(F10.0,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QZE_S1(L),
     &   QZE_S2(L),QZE_D(L),QBE_S1(L),QBE_S2(L),QBE_D(L),QZE(L),QBE(L),
     &   QZE_A(L),QBE_A(L)
      ENDIF
      IC     = 1000*L+17                                  ! file OUTAERO
      ID     = 1000*L+117                                 ! file OUTAERO
      RXO = OSKU(L)%CMX
      RYO = OSKU(L)%CMY
      RZO = OSKU(L)%CMZ
c      CALL FINEULD(ZETAH(L)-ANGLE,-BETAH(L),THETA(L),RXO,RYO,RZO,
c     &   RX,RY,RZ,2)
      CALL FINEULD(ZETAH(L)-ANGLE,0.0,0.0,RXO,RYO,RZO,
     &   RX,RY,RZ,2)
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QBE_A(L),
     &   QZE_A(L),RX,RY,RZ
      WRITE(ID,"(F8.5,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QBE_A(L),
     &   QZE_A(L),RX,RY,RZ
      ELSE
      WRITE(IC,"(F10.0,11(1X,D13.6))") TIME,ANGLE*RAD2DEG,QBE_A(L),
     &   QZE_A(L),RX,RY,RZ
      ENDIF
      IC     = 1000*L+16                                   ! File OUTANG
      ID     = 1000*L+116                                  ! File OUTANG
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,9(1X,D13.6))")TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & THETA(L)*RAD2DEG,ZETAH(L)*RAD2DEG,BETAH(L)*RAD2DEG,ZETAC,BETAC,
     & (THETA(L)-Y(4,L)/ROTORW)*RAD2DEG
      WRITE(ID,"(F8.5,9(1X,D13.6))")TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & THETA(L)*RAD2DEG,ZETAH(L)*RAD2DEG,BETAH(L)*RAD2DEG,ZETAC,BETAC,
     & (THETA(L)-Y(4,L)/ROTORW)*RAD2DEG
      ELSE
      WRITE(IC,"(F10.0,19(1X,D13.6))")TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     & THETA(L)*RAD2DEG,ZETAH(L)*RAD2DEG,BETAH(L)*RAD2DEG,ZETAC,BETAC,
     & (THETA(L)-Y(4,L)/ROTORW)*RAD2DEG
      ENDIF
      IC     = 1000*L+18                                   ! File OUTOSKU
      ID     = 1000*L+118                                  ! File OUTOSKU
      RXO = OSKU(L)%CX
      RYO = OSKU(L)%CY
      RZO = OSKU(L)%CZ
      RX  = OSKU(L)%CMX
      RY  = OSKU(L)%CMY
      RZ  = OSKU(L)%CMZ
      IF(TIMEL) THEN
      WRITE(IC,"(F8.5,1X,11(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     &   RXO,RYO,RZO,RX,RY,RZ
      WRITE(ID,"(F8.5,1X,11(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     &   RXO,RYO,RZO,RX,RY,RZ
      ELSE
      WRITE(IC,"(F10.0,1X,11(D13.6,1X))") TIME,(ZETAH(L)-ANGLE)*RAD2DEG,
     &   RXO,RYO,RZO,RX,RY,RZ
      ENDIF

c      ENDIF ! IPRINT

c      ENDDO
         
c      IF(IPR == ICYCLE) IPR = IPR + IPRINT ! Printing interval

c      IF(TIME < TMAX) GO TO 1000

      END SUBROUTINE BLADE_OUTPUT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE HINGE_MOMENT(L)

      USE BLADE_CONSTANTS

      USE BLADE_VARIABLES, ONLY : RMHINGEI,RMHINGEJ,RMHINGEK,RHINGE,
     & Y,RD,RKZE,RKBE,RKD,QZE_S1,QZE_S2,QBE_S1,QBE_S2,
     & QZE_D,QBE_D,QZE,QBE,CD,QTH,RKTH,THETA,RCZE,RCBE

      USE FLIGHT,          ONLY : NGRIFL,OSKU 
      USE NS3CO,           ONLY : ROTORW

      IMPLICIT NONE

      INTEGER :: L,LP1,LM1

      REAL :: APU1
CCCC HUOMIO HUOMIO! TAMA ALIOHJELMA VOI LAAHATA/OLLA ASKELEEN JALJESSA/EDESSA
CCCC PITANEE SELVITTAA JOSKUS ...
C ... Inertial moments due to the hinge offset

      RMHINGEI(L) = 0.
      RMHINGEJ(L) = OSKU(L)%RMASS*RHINGE(L)*(OSKU(L)%RCG
     &     -RHINGE(L))*ROTORW**2
     &         * SIN(Y(2,L))*COS(Y(1,L))
      RMHINGEK(L) = OSKU(L)%RMASS*RHINGE(L)*(OSKU(L)%RCG-RHINGE(L))
     &         * ROTORW**2*COS(Y(2,L))*SIN(Y(1,L))
c       write(454+l,*) rmhingej(l),rmhingek(l),rm(l),rhinge(l),
c     & omega,y(1,L),y(2,l)
 
C ... Moments caused by the elastomer and the lag dampers

C ... Elastomer:
      QZE_S1(L) = - RKZE(L)*Y(1,L) - RCZE(L)*Y(3,L)
      QBE_S1(L) = - RKBE(L)*Y(2,L) - RCBE(L)*Y(4,L)

C ... Two dampers for each blade:
      CALL INDECES(L,LP1,LM1,NGRIFL)
      APU1     = 0.5*(RD(L)-RHINGE(L))/(RD(L)+1.E-10)

C ... Stiffnesses:                           
      QZE_S2(L) =  0.5*RKD(L)*(RD(L)-RHINGE(L))**2*(Y(1,LP1)-2.*Y(1,L)
     &           + Y(1,LM1) + .5*(Y(2,LP1)**2+2.*Y(2,L)**2+Y(2,LM1)**2)
     &           - APU1*((Y(2,LP1)-Y(2,L))**2 + (Y(2,L)-Y(2,LM1))**2))

      QBE_S2(L) = -0.5*RKD(L)*(RD(L)-RHINGE(L))**2*((Y(1,LP1)-Y(1,L)
     &           +0.5*(Y(2,LP1)**2+Y(2,L)**2)-APU1*(Y(2,LP1)-Y(2,L))**2)
     &           * (Y(2,LP1)-Y(2,L)) 
     &           + (Y(1,LM1)-Y(1,L)+.5*(Y(2,LM1)**2+Y(2,L)**2)
     &           - APU1*(Y(2,LM1)-Y(2,L))**2)*(Y(2,LM1)-Y(2,L)))

C ... Damping moments:
      QZE_D(L) =  0.5*CD(L)*(RD(L)-RHINGE(L))**2*(Y(3,LP1)-2.*Y(3,L)
     &           + Y(3,LM1))
      QBE_D(L) = -0.5*CD(L)*(RD(L)-RHINGE(L))*((Y(3,LP1)-Y(3,L))
     &           * (Y(2,LP1)-Y(2,L))+(Y(3,LM1)-Y(3,L))*(Y(2,LM1)
     &           - Y(2,L)))

C ... Add together
      QZE(L)   = QZE_S1(L) + QZE_S2(L) + QZE_D(L)
      QBE(L)   = QBE_S1(L) + QBE_S2(L) + QBE_D(L)
      QTH(L)   = -RKTH(L)*THETA(L)

      RETURN
      END SUBROUTINE HINGE_MOMENT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AERODYN_MOMENT(L)

      USE BLADE_CONSTANTS
      USE BLADE_VARIABLES, ONLY : QZE_A,QBE_A,QTH_A,DTB
      USE FLIGHT, ONLY          : PSIR,OSKU
      USE NS3CO, ONLY           : TIMEL,ROTORW

      IMPLICIT NONE
      INTEGER :: L
      REAL    :: RX,RY,RZ,RX0,RY0,RZ0

C ... Moment components from the FINFLO coord sys to the helicopter coord sys
      RX0 = OSKU(L)%CMX
      RY0 = OSKU(L)%CMY
      RZ0 = OSKU(L)%CMZ
      IF(TIMEL) THEN
         CALL FINEUL(PSIR(L)+DTB*ROTORW,0.0,0.0,
     +        RX0,RY0,RZ0,RX,RY,RZ,2)
      ELSEIF(.NOT. TIMEL) THEN
         CALL FINEUL(PSIR(L),0.0,0.0,
     +        RX0,RY0,RZ0,RX,RY,RZ,2)
      ENDIF

C ... Primitive quasi-steady aerodynamic moments for NH90 (ISA, SL)
C ... Moments are behind comment (c_priaero)
C ... Moments can be used by removing comments
c_priaero      THETAL   = THETA(L) + 10.8*DEG2RAD ! This shifting must be here
c_priaero      QZE_A(L) = (ROTORW-Y(3,L))**2*(619.*(THETAL-Y(4,L)/ROTORW)**2
c_priaero     &           - 117.*(THETAL-Y(4,L)/ROTORW)+8.5)
c_priaero      QBE_A(L) = (1230.*(THETAL-Y(4,L)/ROTORW)-116)*ROTORW**2
c_priaero      RZ       = QZE_A(L)
c_priaero      RX       = -QBE_A(L)
c_priaero      RY       = 0.

C ... Store moment components to the helicopter positions
      QZE_A(L) = RZ
      QBE_A(L) =-RX
      QTH_A(L) = RY

      RETURN
      END SUBROUTINE AERODYN_MOMENT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BLADE_SOLVER_RK4(L)

      USE BLADE_CONSTANTS

      USE BLADE_VARIABLES, ONLY : PSI,Y,YOLD,F,DTB

      IMPLICIT NONE

      INTEGER :: N,L

      REAL :: PP(4,4),QQ(4,4),RR(4,4),SS(4,4),YY(4,4)

c     write(511,*) l,psim(l)
C ... Solution using 4th-order Runge-Kutta

c      DO L = 1,NGRIFL

C      APU1       = (L-1)*SHIFT
C      ANGLE      = -PSI - APU1 - .5*DTB*OMEGA
C      FHINGEX(L) = -OMEGA**2*RHINGE(L)*COS(ANGLE)*RM(L)
C      FHINGEY(L) = -OMEGA**2*RHINGE(L)*SIN(ANGLE)*RM(L)
C      FHINGEZ(L) = 0. ! Unused forces for comparison to 3dof
C      CALL FINEUL(Y(1,L)+ANGLE,THETA(L),Y(2,L),FHINGEX(L),FHINGEY(L),
C     &     FHINGEZ(L),FHINGEI(L),FHINGEJ(L),FHINGEK(L),2)

C ... First step
      CALL HINGE_MOMENT(L)
      CALL AERODYN_MOMENT(L)
      CALL BLADE_DERIV(L)
      DO N = 1,4
      YY(N,L) = YOLD(N,L) + .5*DTB*F(N,L)
      PP(N,L) = F(N,L)
      ENDDO

C      APU1       = (L-1)*SHIFT
C      ANGLE      = -PSI - APU1 - .5*DTB*OMEGA
C      FHINGEX(L) = -OMEGA**2*RHINGE(L)*COS(ANGLE)*RM(L)
C      FHINGEY(L) = -OMEGA**2*RHINGE(L)*SIN(ANGLE)*RM(L)
C      FHINGEZ(L) = 0. !; FHINGEX(L) = 0. ; FHINGEY(L) = 0. 
C      CALL FINEUL(YY(1,L)+ANGLE,THETA(L),YY(2,L),FHINGEX(L),FHINGEY(L),
C     &     FHINGEZ(L),FHINGEI(L),FHINGEJ(L),FHINGEK(L),2)

C ... Second step
      CALL HINGE_MOMENT(L)
      CALL AERODYN_MOMENT(L)
      CALL BLADE_DERIV(L)
      DO N = 1,4
      YY(N,L) = YOLD(N,L) + .5*DTB*F(N,L)
      QQ(N,L) = F(N,L)
      ENDDO

C ... Third step
      CALL HINGE_MOMENT(L)
      CALL AERODYN_MOMENT(L)
      CALL BLADE_DERIV(L)
      DO N = 1,4
      YY(N,L) = YOLD(N,L) + .5*DTB*F(N,L)
      RR(N,L) = F(N,L)
      ENDDO

C ... Fourth step
      CALL HINGE_MOMENT(L)
      CALL AERODYN_MOMENT(L)
      CALL BLADE_DERIV(L)
      DO N = 1,4
      YY(N,L) = YOLD(N,L) + .5*DTB*F(N,L) ! Useless
      SS(N,L) = F(N,L)
      ENDDO

C ... Final values
      DO N = 1,4
      Y(N,L) = YOLD(N,L) + DTB*(PP(N,L)+2.*QQ(N,L)+2.*RR(N,L)+SS(N,L))/6.
      ENDDO

c      ENDDO

      DO N = 1,4
      YOLD(N,L)  = Y(N,L)
      ENDDO

      RETURN
      END SUBROUTINE BLADE_SOLVER_RK4
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BLADE_DERIV(L)

C ... Gives rates of zeta, beta, zetadot and betadot in vector F

      USE BLADE_CONSTANTS

      USE BLADE_VARIABLES, ONLY : F,Y,RMHINGEI,RMHINGEJ,RMHINGEK,
     &   QZE,QBE,RIZE,RIBE,DISSZE,DISSBE,QZE_A,QBE_A,QTH

      USE CONSTANTS, ONLY : PII

      USE NS3CO, ONLY : TIMEL,ROTORW
      
      USE FLIGHT, ONLY : OSKU,PSIR

      IMPLICIT NONE

      INTEGER :: L
      REAL    :: APUB,APUZ,DZ1,DB1

C ... Calculate the derivatives

C ... Artificial damping
      IF(.NOT. TIMEL) THEN
      DZ1    = DISSZE(L)
      DB1    = DISSBE(L)
      APUZ   = MAX(1.E-2,.25*PII-ABS(Y(1,L))) ! Try to stop the movement
      APUB   = MAX(1.E-2,.25*PII-ABS(Y(2,L))) ! at large angles
      ELSE
         DZ1    = 0.
         DB1    = 0.
         APUZ   = 1.e5! Mita nama on? TIMEL-laskussa laitettu arvoiksi 1, etta
         APUB   = 1.e5! paastaa PATRIAN 5.6.2009 muistion tuloksiin
      ENDIF           ! NH90:n roottorimallinnuksen tukitoimia... (kuvat 8-13)

C ... Calculate the derivatives
      F(1,L) = Y(3,L)
      F(2,L) = Y(4,L)
      F(3,L) = (-RIZE(L)*2.*SIN(Y(2,L))*COS(Y(2,L))*(ROTORW-Y(3,L))
     &       * Y(4,L)-RMHINGEK(L)+ QZE(L)+QZE_A(L))
     &       /(RIZE(L)*COS(Y(2,L))**2)
     &       - DZ1*Y(3,L)!*ABS(Y(3,L))    ! Damping, alternatives<+
     &       - Y(3,L)/APUZ**2             ! Hold your horses, angle too large
      F(4,L) = (-RIBE(L)*SIN(Y(2,L))*COS(Y(2,L))*(ROTORW-Y(3,L))**2     
     &        - RMHINGEJ(L) + QBE(L)+QBE_A(L)) / RIBE(L)
     &       - DB1*Y(4,L)!*ABS(Y(4,L))
     &       - Y(4,L)/APUB**2             ! Hold your horses, angle too large

C ... Centrifugal data to the FINFLO coordinate system
      OSKU(L)%FXE = SIN(PSIR(L))*ROTORW**2*OSKU(L)%RCG*OSKU(L)%RMASS
      OSKU(L)%FYE = 0.0
      OSKU(L)%FZE = -COS(PSIR(L))*ROTORW**2*OSKU(L)%RCG*OSKU(L)%RMASS
      OSKU(L)%MXE = COS(PSIR(L))*
     &             (-RIBE(L)*SIN(Y(2,L))*COS(Y(2,L))*(ROTORW-Y(3,L))**2     
     &             -RMHINGEJ(L))
      OSKU(L)%MYE =-(-RIZE(L)*2.*SIN(Y(2,L))*COS(Y(2,L))*(ROTORW-Y(3,L))
     &       * Y(4,L)-RMHINGEK(L)) ! pitaa olla tuo eka miinus
      OSKU(L)%MZE = SIN(PSIR(L))*
     &             (-RIBE(L)*SIN(Y(2,L))*COS(Y(2,L))*(ROTORW-Y(3,L))**2     
     &             -RMHINGEJ(L))

C ... Damper data to the FINFLO coordinate system
      OSKU(L)%FXH = 0
      OSKU(L)%FYH = 0
      OSKU(L)%FZH = 0
      OSKU(L)%MXH = COS(PSIR(L))*(QBE(L)-(DB1*Y(4,L)+(Y(4,L)/APUB**2))*
     &     RIBE(L))+SIN(PSIR(L))*QTH(L)
      OSKU(L)%MYH = -QZE(L)+(DZ1*Y(3,L)+(Y(3,L)/APUZ**2))*
     &     (RIZE(L)*COS(Y(2,L))**2) 
      OSKU(L)%MZH = SIN(PSIR(L))*(QBE(L)-(DB1*Y(4,L)+(Y(4,L)/APUB**2))*
     &     RIBE(L))-COS(PSIR(L))*QTH(L)

      RETURN
      END SUBROUTINE BLADE_DERIV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INDECES(L,LP1,LM1,NGRIFL)

C ... This subroutine determines L+1 and L-1 indeces

      IMPLICIT NONE

      INTEGER :: L, LP1, LM1, NGRIFL

      IF(NGRIFL <= 2) THEN
         WRITE(*,*)'The number of blades is too small for damping model'
         STOP 'INDECES'
      ENDIF

      IF(L == 1) THEN
         LP1 = L + 1
         LM1 = NGRIFL
      ELSE IF(L == NGRIFL) THEN
         LP1 = 1
         LM1 = L - 1
      ELSE
         LP1 = L + 1
         LM1 = L - 1
      ENDIF

      RETURN
      END SUBROUTINE INDECES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FINEULD(PSI,THETA,PHI,XO,YO,ZO,XONEW,YONEW,ZONEW,ICASE)

C ... FINFLO (PLOT3D) coordinates -> Local body coordinates. 

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE
      REAL :: PSI, THETA, PHI
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: S1, S2, S3, C1, C2, C3
      REAL :: XO, YO, ZO, XN, YN, ZN
      REAL :: XFM, YFM, ZFM, XONEW, YONEW, ZONEW     

C ... Transformation matrix.

      S1 = SIN(PHI)
      S2 = SIN(THETA)
      S3 = SIN(PSI)
      C1 = COS(PHI)
      C2 = COS(THETA)
      C3 = COS(PSI)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A21 =  C2*S3
      A31 = -S2
      A12 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A32 =  C2*S1
      A13 =  C1*C3*S2 + S1*S3
      A23 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A21 =  C3*S1*S2 + C1*S3
      A31 = -C2*S1
      A12 = -C2*S3 
      A22 =  C2*C3
      A32 =  S2
      A13 =  C3*S1 + C1*S2*S3
      A23 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A21 =  S2
      A31 = -C2*S3
      A12 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A32 =  C3*S1 + C1*S2*S3
      A13 =  C3*S1*S2 + C1*S3
      A23 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A21 =  C1*C3*S2 + S1*S3
      A31 =  C1*S2*S3 - C3*S1
      A12 = -S2
      A22 =  C2*C3
      A32 =  C2*S3
      A13 =  C2*S1
      A23 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A21 =  C3*S1 + C1*S2*S3
      A31 =  S1*S3 - C1*C3*S2
      A12 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A32 =  C3*S1*S2 + C1*S3
      A13 =  S2
      A23 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A21 =  C2*S1
      A31 =  C3*S1*S2 - C1*S3
      A12 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A32 =  C1*C3*S2 + S1*S3
      A13 =  C2*S3
      A23 = -S2
      A33 =  C2*C3

      END SELECT


C ... From FINFLO coordinate system to the flight mechanics 
C ... coordinate system.

      XFM = -XO
      YFM = -ZO 
      ZFM = -YO

C ... Rotations according to the Euler angles.

      XN = A11*XFM + A21*YFM + A31*ZFM
      YN = A12*XFM + A22*YFM + A32*ZFM
      ZN = A13*XFM + A23*YFM + A33*ZFM

      XONEW = XN
      YONEW = YN
      ZONEW = ZN

      RETURN
      END SUBROUTINE FINEULD
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE EULFIND(PSI,THETA,PHI,XO,YO,ZO,XONEW,YONEW,ZONEW,ICASE)


C ... Local body coordinates -> FINFLO (PLOT3D) coordinates.

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE

      REAL :: PSI, THETA, PHI
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XO, YO, ZO, XN, YN, ZN
      REAL :: XONEW, YONEW, ZONEW     
      REAL :: S1, S2, S3, C1, C2, C3

      S1 = SIN(PHI)
      S2 = SIN(THETA)
      S3 = SIN(PSI)
      C1 = COS(PHI)
      C2 = COS(THETA)
      C3 = COS(PSI)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A12 =  C2*S3
      A13 = -S2
      A21 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A23 =  C2*S1
      A31 =  C1*C3*S2 + S1*S3
      A32 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A12 =  C3*S1*S2 + C1*S3
      A13 = -C2*S1
      A21 = -C2*S3 
      A22 =  C2*C3
      A23 =  S2
      A31 =  C3*S1 + C1*S2*S3
      A32 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A12 =  S2
      A13 = -C2*S3
      A21 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A23 =  C3*S1 + C1*S2*S3
      A31 =  C3*S1*S2 + C1*S3
      A32 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A12 =  C1*C3*S2 + S1*S3
      A13 =  C1*S2*S3 - C3*S1
      A21 = -S2
      A22 =  C2*C3
      A23 =  C2*S3
      A31 =  C2*S1
      A32 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A12 =  C3*S1 + C1*S2*S3
      A13 =  S1*S3 - C1*C3*S2
      A21 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A23 =  C3*S1*S2 + C1*S3
      A31 =  S2
      A32 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A12 =  C2*S1
      A13 =  C3*S1*S2 - C1*S3
      A21 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A23 =  C1*C3*S2 + S1*S3
      A31 =  C2*S3
      A32 = -S2
      A33 =  C2*C3

      END SELECT


C ... Rotations according to the Euler angles.

      XN = A11*XO + A21*YO + A31*ZO
      YN = A12*XO + A22*YO + A32*ZO
      ZN = A13*XO + A23*YO + A33*ZO


C ... Global directions back to the FINFLO directions.

      XONEW = -XN
      YONEW = -ZN
      ZONEW = -YN

      RETURN
      END SUBROUTINE EULFIND
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE CONEDISK(X,Y,Z,XNEW,YNEW,ZNEW,NGPTS,
     &                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,CAC)

C ... This subroutine moves, rotates, and changes the cone angle of
C ... an actuator disk block. The manipulation is done in three phases.

      IMPLICIT NONE

      REAL, DIMENSION(*) :: X, Y, Z, XNEW, YNEW, ZNEW  
      REAL :: ANGLE 
      REAL :: X1, Y1, Z1, X2, Y2, Z2, CAC
      REAL :: X3, Y3, Z3, X4, Y4, Z4
      REAL :: PX, PY, PZ, PD
      REAL :: QX, QY, QZ, QD
      REAL :: RX2, RY2, RZ2, DR2
      REAL :: A, B, C, COSTHE, SINTHE
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XW, YW, ZW

      INTEGER :: N, NGPTS

C ........................................................................
C ... Phase 1: Move and rotate the disk so that its center point
C ... lies in 0,0,0 and so that its rotation axis is y-axis. 
C ........................................................................

      PX  = X2 - X1
      PY  = Y2 - Y1
      PZ  = Z2 - Z1
      PD  = SQRT(PX*PX+PY*PY+PZ*PZ)
      PX  = PX/PD
      PY  = PY/PD
      PZ  = PZ/PD
 
      QX  = 0.0
      QY  = 1.0
      QZ  = 0.0

      RX2 = PY*QZ - PZ*QY
      RY2 = PZ*QX - PX*QZ
      RZ2 = PX*QY - PY*QX
      DR2 = SQRT(RX2*RX2+RY2*RY2+RZ2*RZ2)

      IF(DR2 > 0.0) THEN
         A = RX2/DR2
         B = RY2/DR2
         C = RZ2/DR2
      ELSE
         A = 1.0
         B = 0.0
         C = 0.0
      ENDIF

      ANGLE = ACOS(PX*QX + PY*QY + PZ*QZ)

      ANGLE  = -ANGLE

      COSTHE = COS(ANGLE)
      SINTHE = SIN(ANGLE)

C ... Rotation matrix

      A11 = A*A*(1.0-COSTHE) +   COSTHE
      A12 = B*A*(1.0-COSTHE) + C*SINTHE
      A13 = C*A*(1.0-COSTHE) - B*SINTHE
      A21 = A*B*(1.0-COSTHE) - C*SINTHE
      A22 = B*B*(1.0-COSTHE) +   COSTHE
      A23 = C*B*(1.0-COSTHE) + A*SINTHE
      A31 = A*C*(1.0-COSTHE) + B*SINTHE
      A32 = B*C*(1.0-COSTHE) - A*SINTHE
      A33 = C*C*(1.0-COSTHE) +   COSTHE


      DO N = 1,NGPTS

         XNEW(N) = A11*(X(N)-X1) 
     &           + A12*(Y(N)-Y1) 
     &           + A13*(Z(N)-Z1)
         YNEW(N) = A21*(X(N)-X1) 
     &           + A22*(Y(N)-Y1) 
     &           + A23*(Z(N)-Z1)
         ZNEW(N) = A31*(X(N)-X1) 
     &           + A32*(Y(N)-Y1) 
     &           + A33*(Z(N)-Z1)

      ENDDO

C ........................................................................
C ... Phase 2: Manipulate the actuator disk cone angle. 
C ........................................................................

      DO N = 1,NGPTS

         XW = XNEW(N)
         YW = YNEW(N)
         ZW = ZNEW(N)

         IF(XW /= 0.0 .OR. ZW /=0.0) THEN

            PX  = XW
            PY  = YW
            PZ  = ZW
            PD  = SQRT(PX*PX+PY*PY+PZ*PZ)
            PX  = PX/PD
            PY  = PY/PD
            PZ  = PZ/PD
 
            QX  = XW
            QY  = YW - 0.1
            QZ  = ZW
            QD  = SQRT(QX*QX+QY*QY+QZ*QZ)
            QX  = QX/QD
            QY  = QY/QD
            QZ  = QZ/QD

            RX2 = PY*QZ - PZ*QY
            RY2 = PZ*QX - PX*QZ
            RZ2 = PX*QY - PY*QX
            DR2 = SQRT(RX2*RX2+RY2*RY2+RZ2*RZ2)

            IF(DR2 > 0.0) THEN
               A = RX2/DR2
               B = RY2/DR2
               C = RZ2/DR2
            ELSE
               A = 1.0
               B = 0.0
               C = 0.0
            ENDIF

            ANGLE = CAC

            COSTHE = COS(ANGLE)
            SINTHE = SIN(ANGLE)

C ... Rotation matrix

            A11 = A*A*(1.0-COSTHE) +   COSTHE
            A12 = B*A*(1.0-COSTHE) + C*SINTHE
            A13 = C*A*(1.0-COSTHE) - B*SINTHE
            A21 = A*B*(1.0-COSTHE) - C*SINTHE
            A22 = B*B*(1.0-COSTHE) +   COSTHE
            A23 = C*B*(1.0-COSTHE) + A*SINTHE
            A31 = A*C*(1.0-COSTHE) + B*SINTHE
            A32 = B*C*(1.0-COSTHE) - A*SINTHE
            A33 = C*C*(1.0-COSTHE) +   COSTHE

            XNEW(N) = A11*XW 
     &              + A12*YW 
     &              + A13*ZW
            YNEW(N) = A21*XW 
     &              + A22*YW 
     &              + A23*ZW
            ZNEW(N) = A31*XW 
     &              + A32*YW 
     &              + A33*ZW

         ENDIF

      ENDDO

C ........................................................................
C ... Phase 3: Move and rotate the disk to the desired tilt and location.
C ........................................................................

      PX  = X4 - X3
      PY  = Y4 - Y3
      PZ  = Z4 - Z3
      PD  = SQRT(PX*PX+PY*PY+PZ*PZ)
      PX  = PX/PD
      PY  = PY/PD
      PZ  = PZ/PD
 
      QX  = 0.0
      QY  = 1.0
      QZ  = 0.0

      RX2 = PY*QZ - PZ*QY
      RY2 = PZ*QX - PX*QZ
      RZ2 = PX*QY - PY*QX
      DR2 = SQRT(RX2*RX2+RY2*RY2+RZ2*RZ2)

      IF(DR2 > 0.0) THEN
         A = RX2/DR2
         B = RY2/DR2
         C = RZ2/DR2
      ELSE
         A = 1.0
         B = 0.0
         C = 0.0
      ENDIF

      ANGLE = ACOS(PX*QX + PY*QY + PZ*QZ)

      COSTHE = COS(ANGLE)
      SINTHE = SIN(ANGLE)

C ... Rotation matrix

      A11 = A*A*(1.0-COSTHE) +   COSTHE
      A12 = B*A*(1.0-COSTHE) + C*SINTHE
      A13 = C*A*(1.0-COSTHE) - B*SINTHE
      A21 = A*B*(1.0-COSTHE) - C*SINTHE
      A22 = B*B*(1.0-COSTHE) +   COSTHE
      A23 = C*B*(1.0-COSTHE) + A*SINTHE
      A31 = A*C*(1.0-COSTHE) + B*SINTHE
      A32 = B*C*(1.0-COSTHE) - A*SINTHE
      A33 = C*C*(1.0-COSTHE) +   COSTHE


      DO N = 1,NGPTS

         XW = XNEW(N)
         YW = YNEW(N)
         ZW = ZNEW(N)

         XNEW(N) = A11*XW 
     &           + A12*YW 
     &           + A13*ZW + X3
         YNEW(N) = A21*XW 
     &           + A22*YW 
     &           + A23*ZW + Y3
         ZNEW(N) = A31*XW 
     &           + A32*YW 
     &           + A33*ZW + Z3

      ENDDO

      RETURN
      END SUBROUTINE CONEDISK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TQIHM(IACTU,ICYCLE)

      USE CONSTANTS, ONLY : PII
      USE NS3CO,     ONLY : FRSVEL, ICYOLD, FRSDEN
      USE FLIGHT,    ONLY : RIN,ROUT,THRUST,TORQUE,RTMSP,FDSP,UTSP,
     &                      FXSP,FXTSP,ADV,ROTA1,ROTB1
      USE BLADE_VARIABLES , ONLY : SHAFT,QTH,QTH_A,QBE_A,QZE_A

C ... This routine has been written by Miklos Lakatos.



C ... Program for solving the torque which corresponds the given thrust 
C ... value,usng the Interval halving method
      IMPLICIT NONE
C ... INPUT PARAMETERS from *.INPUT file
      REAL:: R1,R2,RTM,FD,UT
      CHARACTER(25):: POWFILE,COMMENTS 
      CHARACTER(LEN=3) :: IGRNUM
C ... INPUT PARAMETERS from TRAJECTORY_001.DAT file
      REAL:: FX,FXT
C ... Parameters to be computed
      REAL:: A0,VA,WFT,THDF,T,Q,DM
	 
C ... INETRVAL HALVING METHOD for defining the torque corresponding 
C ... to thrust target value from propeller open water curves
      INTEGER ::I, K, L, M, X, IMAX, NR, IACTU, ICYCLE
      REAL :: JA, JB, T1, T2, T3, MAXER,NOUT
      REAL :: JOUT, KTOUT, KQOUT, ETA0OUT
      REAL :: JJ(3), N(3), QQ(3), ETA0(3), TT(3)
      REAL :: KTI(100), KQI(100)
      REAL :: KT(100), KQ(100), J(100)

      REAL ::  TVECX,TVECY,TVECZ

      LOGICAL :: THERE

C***********************************************************************
cc	OPEN(8,FILE='INPUT.dat')
	
cc 	READ(8,*) AREF		! Ship wetted surface[m^2]
cc	READ(8,*) CHLREF	! Ship length LOS [m]
C ... AREF, CHLREF are not needed.
cc	READ(8,*) FRSVEL	! Free surface (water)speed [m/s]
cc	READ(8,*) FRSDEN	! Free surface (water) density [kg/m^3]
cc	READ(8,*) RIN		! Prop. hub (Act.Disk inner) radius [m]
cc	READ(8,*) ROUT		! Prop. (Act.Disk outer) radius [m]
cc	READ(8,*) RTM		! Resistance in Towing [N]
cc	READ(8,*) FD		! Skin Friction Correction Force[N]
cc	READ(8,*) UT		! Total inflo velocity at PP [m/s]
cc	READ(8,*) FX		! Resistance in Self propulsion [N]
cc	READ(8,*) FXT		! Thrust applyed on to Act.Disc
cc	READ(8,*) POWFILE	! File with Propeller Open Water Curves
cc	CLOSE(8)
C ... FRSVEL,FRSDEN,RIN,ROUT are to be taken from *.INPUT file
C ... RTM,FD,need to be added to  *.INPUT file and taken from there.




C ... Propeller open water curves 
c	OPEN(5,FILE=POWFILE)
c	OPEN(5,FILE='POW_p2411.dat')
        CALL NUMCH3(IGRNUM,IACTU)
        INQUIRE(FILE='OPENWATER_'//IGRNUM//'.dat',
     &       EXIST=THERE)
        IF(THERE) THEN
           IF (ICYCLE-ICYOLD < 2) THEN
              WRITE(*,*)'File OPENWATER_'//IGRNUM//'.dat',
     &             ' is found. Self propulsion calculation ',
     &             'will be performed.'
           ENDIF          
           OPEN(300+IACTU, FILE='OPENWATER_'//IGRNUM//'.dat', ! Read data in 
     &          STATUS='UNKNOWN', FORM='FORMATTED')
           READ(300+IACTU,*) NR
           READ(300+IACTU,*) COMMENTS
           READ(300+IACTU,*) COMMENTS
           DO I=1, NR
              READ(300+IACTU,*) J(I), KT(I), KQ(I) 
           END DO       
           CLOSE(300+IACTU)
        ELSE 
           IF (ICYCLE-ICYOLD < 2) THEN
              WRITE(*,*)'File OPENWATER_'//IGRNUM//'.dat',
     &             ' is not found. Exiting.'
              STOP
           ENDIF
           RETURN
        ENDIF

C ... Data is stored to correct place
        FX  = ABS(FXSP(IACTU))
        FXT = ABS(FXTSP(IACTU))
        R2  = ROUT(IACTU)
        R1  = RIN(IACTU)
        RTM = RTMSP(IACTU)
        FD  = FDSP(IACTU)
        UT  = UTSP(IACTU)
ccc        write(*,*)'e',FX,FXT,FD,UT

C ... Propeller disk are in model scale
	DM=2*R2
	A0= PII*R2**2-PII*R1**2
	VA=UT-FXT/(2.*UT*FRSDEN*A0) ! Advance velocity model scale
	WFT=(FRSVEL-VA)/FRSVEL
	THDF=(FXT-RTM+FD)/FXT
	T=FX-FD			! Target thrust
ccc        write(*,*)'e',FRSVEL,A0,R2,R1

C ... Thrust vector changed to the actuator disc coord. 
         CALL FINEUL(0.0,ROTA1(IACTU),ROTB1(IACTU),
     +        -T,0.0,0.0,TVECZ,TVECY,TVECX,1)

cc         write(*,*)TVECX,TVECY,TVECZ
         T = TVECX

C ... Iteration
	IMAX = 200
	MAXER =1.E-5
	
	JA=J(1)
	JB=J(NR)

    	DO M = 1,IMAX
	   JJ(1)=JA+(JB-JA)*0.25
	   JJ(2)=JA+(JB-JA)*0.5
	   JJ(3)=JA+(JB-JA)*0.75
	   DO K=1,3
	      DO L=1,NR
		 IF (JJ(K) == J(L)) THEN
		    KTI(K) = KT(L)				
		    KQI(K) = KQ(L)
		    GOTO 100
		 ELSE IF (JJ(K)==J(L+1)) THEN
		    KTI(K) = KT(L+1)
		    KQI(K) = KQ(L+1)
		    GOTO 100
		 ELSE IF (JJ(K)>J(L) .AND. JJ(K)<J(L+1)) THEN
		    KTI(K) = KT(L)+(KT(L+1)-KT(L))/(J(L+1)-J(L))*
     +		    (JJ(K)-J(L))
		    KQI(K) = KQ(L)+(KQ(L+1)-KQ(L))/(J(L+1)-J(L))*
     +              (JJ(K)-J(L))
		    GOTO 100
		 END IF 
	      END DO       
       
 100	      N(K) = VA/(JJ(K)*DM)
	      TT(K) = KTI(K)*1000.*N(K)**2*DM**4
	      QQ(K) = KQI(K)/10.*1000.*N(K)**2*DM**5
	      ETA0(K) = (JJ(K)*KTI(K)*10.)/(2.*PII*KQI(K))
	      
	      T1 = ABS(TT(1)-T)
	      IF (K >= 2) THEN
		 T2 = ABS(TT(2)-T)
	      END IF
	      IF (K >= 3) THEN
		 T3 = ABS(TT(3)-T)
	      END IF
	   END DO  
	   
	   IF (T1 <= T2 .AND. T2 <= T3) THEN
	      JB = JJ(2)
	   ELSE IF (T1 >= T2 .AND. T2 >= T3) THEN
	      JA = JJ(2)
	   ELSE IF (T1 >= T2 .AND. T2 <= T3) THEN
	      JA = JJ(1)
	      JB = JJ(3)
	   END IF
	   
	   IF (T1 < MAXER .OR. T2 < MAXER .OR. T3 < MAXER)THEN
	      GOTO 200
	   END IF
	END DO

 200	IF (T1 <= T2 .AND. T2 <= T3) THEN
	   X=1
    	ELSE IF (T1 >= T2 .AND. T2 >= T3) THEN
	   X=3
    	ELSE IF (T1 >= T2 .AND. T2 <= T3) THEN
	   X=2
    	ELSE IF (T1 <= T2 .AND. T2 >= T3) THEN
	   X=3
    	END IF

C ... Final values from IHM
	SHAFT(IACTU)  = N(X)
    	ADV(IACTU)    = JJ(X)
    	KTOUT         = KTI(X)
    	KQOUT         = KQI(X)/10.
    	THRUST(IACTU) = TT(X)
    	TORQUE(IACTU) = QQ(X)
    	ETA0OUT       = ETA0(X)

        QTH(IACTU)    = SHAFT(IACTU)
        QTH_A(IACTU)  = ADV(IACTU)
        QBE_A(IACTU)  = THRUST(IACTU)
        QZE_A(IACTU)  = TORQUE(IACTU)
        
ccc      write(*,*)'j',THRUST(IACTU),TORQUE(IACTU),ADV(IACTU),SHAFT(IACTU)
CCC	WRITE(*,FMT='(6A9)')'CHLREF','AREF','FRSDEN','R1','R2','A0'
CCC	WRITE(*,FMT='(6F9.3)')CHLREF,AREF,FRSDEN,R1,R2,A0
CCC	WRITE(*,*)
CCC	WRITE(*,FMT='(7A7)')'FRSVEL','RTM','FD','FX','THRUST','THDF','ITER'
CCC	WRITE(*,FMT='(6F7.3,I7)')FRSVEL,RTM,FD,FX,THRUST(IACTU),THDF,M
CCC	WRITE(*,*)
CCC	WRITE(*,FMT='(4A7)')'FRSVEL','UT','VA','WFT'
CCC	WRITE(*,FMT='(4F7.3)') FRSVEL,UT,VA,WFT
CCC	WRITE(*,*)
CCC	WRITE(*,FMT='(9A7)')'FRSVEL','N','J','KT','KQ','THRUST','TORQUE','ETA0'
CCC	WRITE(*,FMT='(8F7.3,I7)')FRSVEL,NOUT,JOUT,KTOUT,KQOUT,THRUST(IACTU),
CCC	1 TORQUE(IACTU),ETA0OUT

CCC	OPEN(4,FILE='OUTPUT.dat')
CCC	WRITE(4,FMT='(6A9)')'CHLREF','AREF','FRSDEN','R1','R2','A0'
CCC	WRITE(4,FMT='(6F9.3)')CHLREF,AREF,FRSDEN,R1,R2,A0
CCC	WRITE(4,*)
CCC	WRITE(4,FMT='(7A7)')'FRSVEL','RTM','FD','FX','THRUST','THDF','ITER'
CCC	WRITE(4,FMT='(6F7.3,I7)')FRSVEL,RTM,FD,FX,THRUST(IACTU),THDF,M
CCC	WRITE(4,*)
CCC	WRITE(4,FMT='(4A7)')'FRSVEL','UT','VA','WFT'
CCC	WRITE(4,FMT='(4F7.3)') FRSVEL,UT,VA,WFT
CCC	WRITE(4,*)
CCC	WRITE(4,FMT='(9A7)')'FRSVEL','N','J','KT','KQ','THRUST','TORQUE',
CCC	1 'ETA0'
CCC	WRITE(4,FMT='(8F7.3,I7)')FRSVEL,NOUT,JOUT,KTOUT,KQOUT,THRUST(IACTU),
CCC	1 TORQUE(IACTU),ETA0OUT

CCC	CLOSE(4)


      RETURN	
      END SUBROUTINE TQIHM
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE SET_THETAS(IFA,IACTU)

      USE CONSTANTS, ONLY : DEG2RAD
      USE FLIGHT,    ONLY : NGRIFL,OSKU,FXTSP,THRUST,ROUT,
     &     RTMSP,FDSP,ROTA1,ROTB1
      USE BLADE_VARIABLES , ONLY : QTH,QTH_A,QBE_A,SHAFT


C ... This routine is based on Jaakko Hoffren's model sept 2014.

      IMPLICIT NONE
      REAL :: TCUR,TERR,DTHETACOL,XERR,DTHCYCLON,ZERR,DTHCYCLAT
      REAL :: TX,TY,TZ
      INTEGER :: IACTU,IFA


C ... Current Thrust vector

      CALL EULFIN(0.0,ROTB1(IACTU),ROTA1(IACTU),
     &     0.0,0.0,-SHAFT(IACTU)/ABS(SHAFT(IACTU)),TX,TY,TZ,2)
      TCUR = (OSKU(IACTU)%FXT*TX+OSKU(IACTU)%FYT*TY+
     &     OSKU(IACTU)%FZT*TZ)


C ... New THETACOL 

      TERR = TCUR+OSKU(IACTU)%DAMPN*(TCUR-FXTSP(IACTU))-
     &     THRUST(IACTU) 

      DTHETACOL=MAX(-0.1,MIN((-0.4*TERR/OSKU(IACTU)%DAMPT),0.1))
      DTHETACOL=DTHETACOL*DEG2RAD
      QTH(IACTU)=DTHETACOL


C ... New THCYCLON 

      XERR = OSKU(IACTU)%MZT/TCUR+OSKU(IACTU)%DAMPN*
     &     (OSKU(IACTU)%MZT/TCUR-RTMSP(IACTU))      

      DTHCYCLON=MAX(-0.1,MIN((-4*XERR/ROUT(IACTU)),0.1))
      DTHCYCLON=DTHCYCLON*DEG2RAD
      QTH_A(IACTU)=DTHCYCLON
      IF (IFA == 12) THEN 
         QTH_A(IACTU)= 0.0
      ENDIF


C ... New THCYCLAT

      ZERR = OSKU(IACTU)%MXT/TCUR+OSKU(IACTU)%DAMPN*
     &     (OSKU(IACTU)%MXT/TCUR-FDSP(IACTU))

      DTHCYCLAT=MAX(-0.1,MIN((-4*ZERR/ROUT(IACTU)),0.1))
      DTHCYCLAT=DTHCYCLAT*DEG2RAD
      QBE_A(IACTU)=DTHCYCLAT
      IF (IFA == 12) THEN 
         QBE_A(IACTU)= 0.0
      ENDIF

 
      RETURN	
      END SUBROUTINE SET_THETAS
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE BLELEMENT(PSIL,RVECLEN,R1,R2,UM,VM,WM,COM,
     1     cl,cd,alphai,VE,IACTU,NX,NY,NZ,PX,PY,PZ,RX,RY,RZ)

C ... This subroutine calculates lift and drag coeff for the blade element.
C ... Equitations are written by Jaakko Hoffren (1.9.2014).
C ... Programming work is done by Finflo Ltd


      USE CONSTANTS, ONLY : PII,RAD2DEG,DEG2RAD,EPS6

      USE FLIGHT, ONLY : ROTA1,ROTB1,CONEA,NBLADE,OSKU

      USE BLADE_VARIABLES, ONLY : THETACOL,THCYCLON,THCYCLAT,SHAFT,ETIP,
     +     CBLADE,QTH,QTH_A,QBE_A

      USE NS3CO, ONLY : ICYCLE,ICYOLD

      IMPLICIT NONE

      INTEGER :: IACTU

      REAL :: PSIL,RVECLEN,R1,R2,UM,VM,WM,RY,RX,RZ,COM,cl,cd,alphai,VE,
     +     NX,NY,NZ,PX,PY,PZ,VT,alphae,MA,cla,cd0,cls,alphas,FOUT,EFOUT,
     +     FIN,EFIN,R,ai,THET_0,THET_A,THET_B,RY1,RX1,RZ1,R1PRAN,R2PRAN

      REAL :: PSI,VECX,VECY,VECZ


C ... Rotation angle

      PSI = PSIL


C ... Coordinate settings (Euler vs Helicopter ACT disk): x=-p , y=r , z=-n
C ... Unit vector n

      CALL EULFIN(0.0,-SIN(PSI)*CONEA(IACTU)+ROTB1(IACTU),
     +     -COS(PSI)*CONEA(IACTU)+ROTA1(IACTU),
     +     0.0,0.0,-SHAFT(IACTU)/ABS(SHAFT(IACTU)),
     +     VECX,VECY,VECZ,2)
      NX = VECX
      NY = VECY
      NZ = VECZ


C ... Unit vector r

      RX1 = RX/RVECLEN
      RY1 = RY/RVECLEN
      RZ1 = RZ/RVECLEN


C ... Unit vector p

      PX = NZ*RY1-NY*RZ1
      PY = NX*RZ1-NZ*RX1
      PZ = NY*RX1-NX*RY1

         
C ... Rotor RPM is changed to the FINFLO coord. sys. (ROTX,ROTY,ROTZ)

      CALL EULFIN(0.0,ROTB1(IACTU),ROTA1(IACTU),
     +     0.0,0.0,-SHAFT(IACTU)*RAD2DEG,
     +     VECX,VECY,VECZ,2)
      OSKU(IACTU)%ROTX = VECX
      OSKU(IACTU)%ROTY = VECY
      OSKU(IACTU)%ROTZ = VECZ


C ... Tangential velocity component

      VT = UM*PX+VM*PY+WM*PZ+
     +     (OSKU(IACTU)%ROTZ*RY-OSKU(IACTU)%ROTY*RZ)*PX+
     +     (OSKU(IACTU)%ROTX*RZ-OSKU(IACTU)%ROTZ*RX)*PY+
     +     (OSKU(IACTU)%ROTY*RX-OSKU(IACTU)%ROTX*RY)*PZ


C ... Flow field angle compared to the disk plane

      alphai = -atan((UM*NX+VM*NY+WM*NZ)/(VT))


C ... Effective angle of attack

c      alphae = THETACOL(IACTU)-ROTA1(IACTU)*COS(PSI)-
c     +     ROTB1(IACTU)*SIN(PSI)-(0.75-RVECLEN/R2)*ETIP(IACTU)
c     +     -DBLE(alphai)

      IF ((ICYCLE-ICYOLD) < 1) THEN
         THET_0 = THETACOL(IACTU)
         THET_A = THCYCLON(IACTU)
         THET_B = THCYCLAT(IACTU)
      ELSE
         THET_0 = THETACOL(IACTU)+QTH(IACTU)
         THET_A = THCYCLON(IACTU)+QTH_A(IACTU)
         THET_B = THCYCLAT(IACTU)+QBE_A(IACTU)
      ENDIF

      alphae = THET_0-THET_A*COS(PSI)-THET_B*SIN(PSI)-
     +     (0.75-RVECLEN/R2)*ETIP(IACTU)-alphai


C ... Effective flow field velocity
      VE = SQRT(VT**2+(UM*NX+VM*NY+WM*NZ)**2)


C ... Local Mach number

      MA = VE/COM


C ... Prandtl's tip loss correction

      ai     = MAX(alphai,0.0349066)
      R2PRAN = MAX((R2-RVECLEN),EPS6)
      R1PRAN = MAX((RVECLEN-R1),EPS6)
      EFOUT  = (-NBLADE(IACTU)*R2PRAN)/(2*RVECLEN*SIN(ai))
      FOUT   = 2/PII*ACOS(EXP(EFOUT))
      EFIN   = (-NBLADE(IACTU)*R1PRAN)/(2*RVECLEN*SIN(ai))
      FIN    = 2/PII*ACOS(EXP(EFIN))
      alphae = alphae*FOUT*FIN
ccc      cl = cl*F


C ... Aerodynamic coefficient for the profile

      IF (MA <= 0.82) THEN ! subsonic low aoa
         cla = 6/SQRT(1-MA**2)
         cd0 = 0.008
      ELSE ! transonic low aoa
         cla = 4.
         cd0 = 0.008+MIN((3.1*(MA-0.82)**2),0.07)
      ENDIF

      cl  = cla*alphae
      cls = MAX(1.5*COS(0.5*PII*MA**1.5),0.1)
      cd  = cd0+0.03*(cl/cls)**4

      IF (ABS(cl) > ABS(cls)) THEN ! stall
         alphas = cls/cla
         cl = MAX((cls-9*(alphae-alphas)),2*sin(alphae)*cos(alphae))
         
         IF (alphae < 0 ) THEN
            cl =MIN((-cls+9*(-alphae-alphas)),2*sin(alphae)*cos(alphae))
         ENDIF

         IF(ABS(alphae) < 20.*DEG2RAD) THEN
            cd = cd0+0.03+(ABS(alphae)-alphas)/(20.*DEG2RAD-alphas)*
     +           (1.8*(sin(20.*DEG2RAD))**2-cd0-0.03)
         ELSE
            cd = MAX(1.8*(sin(alphae))**2,0.02)
         ENDIF

      ENDIF ! IF (ABS(cl) > ABS(cls)) ...






C ... Files for checking purpose


c          WRITE(707,*)REAL(PSI*57.29577951),RX/RVECLEN,RY/RVECLEN,
c     +     RZ/RVECLEN,RVECLEN
c          WRITE(777,*)REAL(PSI*57.29577951),REAL(NX),REAL(NY),REAL(NZ),
c     +         REAL(sqrt(NX**2+NY**2+NZ**2)-1)
c          WRITE(787,*)REAL(PSI*57.29577951),REAL(PX),REAL(PY),REAL(PZ),
c     +         REAL(sqrt(PX**2+PY**2+PZ**2)-1)
c      IF (PSI >= 1.50 .AND. PSI <= 1.58) THEN
c          WRITE(797,*)RVECLEN,alphae*57.29577951,
c     +         alphai*57.29577951,(THETACOL(IACTU)-
c     +        SIN(PSI)*THCYCLAT(IACTU)-COS(PSI)*THCYCLON(IACTU)-
c     +        (0.75-DBLE(RVECLEN/R2))*ETIP(IACTU))*57.29577951
c          WRITE(727,*)ETIP(IACTU)*57.2957795,CBLADE(IACTU),NBLADE(IACTU)
c          WRITE(737,*)RVECLEN,MA,SQRT(UM**2+VM**2+WM**2)
c          WRITE(757,*)RVECLEN,alphae*57.29577951,UM,VM,WM
c          WRITE(767,*)RVECLEN,alphae*57.29577951,cl,cd,cla
c          WRITE(747,*)RVECLEN,alphae*57.29577951,VT,VE
c          WRITE(717,*)RVECLEN,OSKU(IACTU)%ROTX,OSKU(IACTU)%ROTY,
c     +         OSKU(IACTU)%ROTZ
c          write(555,*)RVECLEN,EFOUT,(EXP(EFOUT)),FOUT,alphai
c          write(556,*)RVECLEN,EFIN,(EXP(EFIN)),FOUT*FIN,alphai
c          write(557,*)OSKU(IACTU)%ROTX,OSKU(IACTU)%ROTY,OSKU(IACTU)%ROTZ
c          write(558,*)PSI*57.29577951,RX,RY,RZ


c       ENDIF
       
      RETURN	
      END SUBROUTINE BLELEMENT
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE FAKEPROPELLER(IGRFAKEP)

      USE FLIGHT,    ONLY : OSKU,XCG,YCG,ZCG,PSIR,DRAUGHTI,DRAUGHT,
     &     XFAKEP,YFAKEP,ZFAKEP,ROTB1FAKEP,IFAKEP,
     &     XFAKEPNEW,YFAKEPNEW,ZFAKEPNEW
      USE CONSTANTS, ONLY : DEG2RAD

      IMPLICIT NONE
      
      INTEGER :: IGRFAKEP
      REAL, DIMENSION(3) :: FVEC, RVEC, MFAKEP
      REAL :: XPROP, YPROP, ZPROP, ROTB1PROP,
     &        XCGAFAKEP, YCGAFAKEP, ZCGAFAKEP, ZCGDFAKEP, ZERO

C ... Calculate new location and angle of the fake propeller
      ZERO = 0.0
      ZCGDFAKEP = ZFAKEP(IGRFAKEP)+
     &     DRAUGHTI(IGRFAKEP)-DRAUGHT(IGRFAKEP) ! add draught
      CALL EULFIN(PSIR(IGRFAKEP),ZERO,ZERO, ! rotation
     &     -(XFAKEP(IGRFAKEP)-XCG(IGRFAKEP)),
     &     -(ZCGDFAKEP-ZCG(IGRFAKEP)),
     &     -(YFAKEP(IGRFAKEP)-YCG(IGRFAKEP)),XCGAFAKEP,
     &     YCGAFAKEP,ZCGAFAKEP,1) 
      XPROP     = XCGAFAKEP+XCG(IGRFAKEP) ! new locations
      YPROP     = YCGAFAKEP+YCG(IGRFAKEP)
      ZPROP     = ZCGAFAKEP+ZCG(IGRFAKEP)
      ROTB1PROP = ROTB1FAKEP(IGRFAKEP)-PSIR(IGRFAKEP)
      XFAKEPNEW(IGRFAKEP)=XPROP ! data is stored
      YFAKEPNEW(IGRFAKEP)=YPROP
      ZFAKEPNEW(IGRFAKEP)=ZPROP

C ... Fake propeller forces in X and Z direction
      OSKU(IGRFAKEP)%FXT = -OSKU(IGRFAKEP)%CX
      OSKU(IGRFAKEP)%FZT = -tan(-ROTB1PROP*DEG2RAD)*
     &     OSKU(IGRFAKEP)%CX

C ... Fake propelle moment in Y direction
      FVEC(1) = OSKU(IGRFAKEP)%FXT  ! Force components
      FVEC(2) = OSKU(IGRFAKEP)%FYT
      FVEC(3) = OSKU(IGRFAKEP)%FZT
      RVEC(1) = XCG(IGRFAKEP)-XPROP ! Shaft componnets
      RVEC(2) = YCG(IGRFAKEP)-YPROP
      RVEC(3) = ZCG(IGRFAKEP)-ZPROP               
      CALL CROSS_PRODUCT(FVEC,RVEC,MFAKEP) ! Cross product
      OSKU(IGRFAKEP)%MYT = MFAKEP(2)  ! MY moment

      RETURN	
      END SUBROUTINE FAKEPROPELLER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FCNQuaternion(X,Y,YPRIME)
C ... This is the "quaternion" version of FCN. So basically all the derivatives are
C ... calculated in the same way as in FCN, but instead of Euler angles, now the
C ... orientation quaternion is used.
C ... April 2020, Matti Palin. 

      USE NS3CO     , ONLY   : GX, GY, GZ
      USE CONSTANTS , ONLY   : PII
      USE TYPE_ARRAYS , ONLY : SIXDOF
      USE FLIGHT
      IMPLICIT NONE

      INTEGER :: ICTL
      REAL, DIMENSION(16) :: Y, YPRIME
      REAL :: X
      REAL :: COSALF, SINALF, COSBET, SINBET      
      REAL :: DML, DMD, DMY, DLL, DLD, DLY, DNL, DND, DNY, DM, DL, DN 
      REAL :: GYROX, GYROY, GYROZ, GYROXP, GYROYP, GYROZP
      REAL :: CROL, CM, CN
      REAL :: CXK, CYK, CZK
      REAL :: REY, MACH, MA
      REAL :: Y1F, Y2F, Y3F, Y1, Y2, Y3, V0, V1, X0, Y0, Z0
      REAL :: FX, FY, FZ, MX, MY, MZ, GXX, GYY, GZZ 
      REAL :: EX, EY, EZ, EMX ,EMY, EMZ
      REAL :: TX, TY, TZ, TMX ,TMY, TMZ
      REAL :: HX, HY, HZ, HMX, HMY, HMZ
      REAL :: OmegaAuxX, OmegaAuxY, OmegaAuxZ

      INTRINSIC SIN 

C ------- PART I : convert forces and moments to the body coordinate system  ----- 

C ... Initialize helicopter forces and moments
      HX = 0. ; HY = 0. ; HZ = 0. ; HMX = 0. ; HMY = 0. ; HMZ = 0.

      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%CX,
     &     OSKU(IPATH)%CY,OSKU(IPATH)%CZ,FX,FY,FZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%CMX,
     &     OSKU(IPATH)%CMY,OSKU(IPATH)%CMZ,MX,MY,MZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),GX,GY,GZ,
     &     GXX,GYY,GZZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%FXE,
     &     OSKU(IPATH)%FYE,OSKU(IPATH)%FZE,EX,EY,EZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%MXE,
     &     OSKU(IPATH)%MYE,OSKU(IPATH)%MZE,EMX,EMY,EMZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%FXT,
     &     OSKU(IPATH)%FYT,OSKU(IPATH)%FZT,TX,TY,TZ)
      CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%MXT,
     &     OSKU(IPATH)%MYT,OSKU(IPATH)%MZT,TMX,TMY,TMZ)

C ... Helicopter case is handled in a different way
      IF (TRMODE(IPATH) == 31) THEN
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%CX,
     &        OSKU(IPATH)%CY,OSKU(IPATH)%CZ,FX,FY,FZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%CMX,
     &        OSKU(IPATH)%CMY,OSKU(IPATH)%CMZ,MX,MY,MZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),GX,GY,GZ,GXX,GYY,GZZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%FXT,
     &        OSKU(IPATH)%FYT,OSKU(IPATH)%FZT,TX,TY,TZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%MXT,
     &        OSKU(IPATH)%MYT,OSKU(IPATH)%MZT,TMX,TMY,TMZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%FXE,
     &        OSKU(IPATH)%FYE,OSKU(IPATH)%FZE,EX,EY,EZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%MXE,
     &        OSKU(IPATH)%MYE,OSKU(IPATH)%MZE,EMX,EMY,EMZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%FXH,
     &        OSKU(IPATH)%FYH,OSKU(IPATH)%FZH,HX,HY,HZ)
         CALL FINEULQua(Y(13),Y(14),Y(15),Y(16),OSKU(IPATH)%MXH,
     &        OSKU(IPATH)%MYH,OSKU(IPATH)%MZH,HMX,HMY,HMZ)
C ... Gravity forces come via TX,TY,TZ
         GXX = 0.0
         GYY = 0.0
         GZZ = 0.0
C ... Centrifugal terms are edited, because YPRIME terms 
C ... contain part of these centrifugal terms
         EX = EX - (-Y(5)*Y(3)+Y(6)*Y(2))*M
         EY = EY - (-Y(6)*Y(1) + Y(4)*Y(3))*M
         EZ = EZ - (-Y(4)*Y(2)+Y(5)*Y(1))*M
         EMX = EMX - (DL- IXY*(Y(4)*Y(6)-QP)+IXZ*(Y(4)*Y(5)+RP)
     &         + IYZ*(Y(5)**2-Y(6)**2)-Y(5)*Y(6)*(IZ-IY)
     &         - GYROXP+Y(6)*GYROY-Y(5)*GYROZ)
         EMY = EMY - (DM + IXY*(Y(5)*Y(6)+PP)
     &         - IXZ*(Y(4)**2-Y(6)**2)-IYZ*(Y(4)*Y(5)-RP)
     &         - Y(4)*Y(6)*(IX-IZ)-GYROYP-Y(6)*GYROX+Y(4)*GYROZ)
         EMZ = EMZ - (DN - IXY*(Y(5)**2-Y(4)**2)-IXZ*(Y(5)*Y(6)-PP)
     &         + IYZ*(Y(4)*Y(6)+QP)-Y(4)*Y(5)*(IY-IX)
     &         - GYROZP+Y(5)*GYROX-Y(4)*GYROY)
      ENDIF ! IF (TRMODE(IPATH) == 31) THEN



C ------- PART II : Initialize a number of values  ----- 
      OmegaAuxX = 0.0
      OmegaAuxY = 0.0
      OmegaAuxZ = 0.0
      GYROX  = 0.0
      GYROY  = 0.0
      GYROZ  = 0.0
      GYROXP = 0.0
      GYROYP = 0.0
      GYROZP = 0.0
      DML    = 0.0
      DMD    = 0.0
      DMY    = 0.0
      DLL    = 0.0
      DLD    = 0.0
      DLY    = 0.0
      DNL    = 0.0
      DND    = 0.0
      DNY    = 0.0
      DM     = 0.0 !DML + DMD + DMY
      DL     = 0.0 !DLL + DLD + DLY
      DN     = 0.0 !DNL + DND + DNY



C ------- PART III : Euler equations  ----- 
C ... Force equations
 
      YPRIME(1) = FX/M + GXX + EX/M + TX/M + HX/M 
     &     - Y(5)*Y(3) + Y(6)*Y(2)

      YPRIME(2) = FY/M + GYY + EY/M + TY/M + HY/M
     &     - Y(6)*Y(1) + Y(4)*Y(3)

      YPRIME(3) = FZ/M + GZZ + EZ/M + TZ/M + HZ/M
     &     - Y(4)*Y(2) + Y(5)*Y(1)

C ... Moment equations
      YPRIME(4) = (MX+EMX+TMX+HMX+DL
     &          - IXY*(Y(4)*Y(6)-QP)+IXZ*(Y(4)*Y(5)+RP)
     &          + IYZ*(Y(5)**2-Y(6)**2)-Y(5)*Y(6)*(IZ-IY)
     &          - GYROXP+Y(6)*GYROY-Y(5)*GYROZ)/IX

      YPRIME(5) = (MY+EMY+TMY+HMY+DM+IXY*(Y(5)*Y(6)+PP)
     &          - IXZ*(Y(4)**2-Y(6)**2)-IYZ*(Y(4)*Y(5)-RP)
     &          - Y(4)*Y(6)*(IX-IZ)-GYROYP-Y(6)*GYROX+Y(4)*GYROZ)/IY

      YPRIME(6) = (MZ+EMZ+TMZ+HMZ+DN
     &          - IXY*(Y(5)**2-Y(4)**2)-IXZ*(Y(5)*Y(6)-PP)
     &          + IYZ*(Y(4)*Y(6)+QP)-Y(4)*Y(5)*(IY-IX)
     &          - GYROZP+Y(5)*GYROX-Y(4)*GYROY)/IZ

C ... Quaternion time derivative

      CALL EULFINQua(Y(13),Y(14),Y(15),Y(16),Y(4),Y(5),Y(6),
     & OmegaAuxX, OmegaAuxY, OmegaAuxZ)

      YPRIME(13) = 0.5*(-OmegaAuxX*Y(14) - OmegaAuxY*Y(15)
     &                -  OmegaAuxZ*Y(16))
      YPRIME(14) = 0.5*( OmegaAuxX*Y(13) - OmegaAuxZ*Y(15)
     &                 + OmegaAuxY*Y(16))
      YPRIME(15) = 0.5*( OmegaAuxY*Y(13) + OmegaAuxZ*Y(14)
     &                -  OmegaAuxX*Y(16))
      YPRIME(16) = 0.5*( OmegaAuxZ*Y(13) - OmegaAuxY*Y(14)
     &                +  OmegaAuxX*Y(15))
      ! If you want to check that this was calculated correctly,
      ! check that the lengh of YPRIME(13)---YPRIME(16) equals 
      ! half of the length of Y(4)---Y(6)

C .... Location equations

      YPRIME(10) = 
     & (Y(13)**2 + Y(14)**2 - Y(15)**2 - Y(16)**2) * Y(1)
     &          + 2.0*(Y(14)*Y(15) - Y(13)*Y(16)) * Y(2)
     &          + 2.0*(Y(14)*Y(16) + Y(13)*Y(15)) * Y(3)

      YPRIME(11) =
     &            2.0*(Y(14)*Y(15) + Y(13)*Y(16)) * Y(1)
     & +(Y(13)**2 - Y(14)**2 + Y(15)**2 - Y(16)**2)* Y(2)
     &          + 2.0*(Y(15)*Y(16) - Y(13)*Y(14)) * Y(3)

      YPRIME(12) =
     &           -2.0*(Y(14)*Y(16) - Y(13)*Y(15)) * Y(1)
     &          - 2.0*(Y(15)*Y(16) + Y(13)*Y(14)) * Y(2)
     & -(Y(13)**2 - Y(14)**2 - Y(15)**2 + Y(16)**2)* Y(3)

C ... Define the time derivatives
      UP     = YPRIME(1)
      VP     = YPRIME(2)
      WP     = YPRIME(3)
      PP     = YPRIME(4)
      QP     = YPRIME(5)
      RP     = YPRIME(6)
      FIIP   = YPRIME(7)
      THETAP = YPRIME(8)
      PSIIP  = YPRIME(9)
      XP     = YPRIME(10)
      YP     = YPRIME(11)
      ZP     = YPRIME(12)

      IF(ABS(YPRIME(1)) <= 1E-15) UP = 0.0
      IF(ABS(YPRIME(2)) <= 1E-15) VP = 0.0
      IF(ABS(YPRIME(3)) <= 1E-15) WP = 0.0
      IF(ABS(YPRIME(4)) <= 1E-15) PP = 0.0
      IF(ABS(YPRIME(5)) <= 1E-15) QP = 0.0
      IF(ABS(YPRIME(6)) <= 1E-15) RP = 0.0


      RETURN
      END SUBROUTINE FCNQuaternion
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FINEULQua(Q0,Q1,Q2,Q3,XO,YO,ZO,XONEW,YONEW,ZONEW)
C ... FINFLO (PLOT3D) coordinates -> Local body coordinates.
C ... But this time using the quaternion representation instead of Euler angles!
C ... So we take the orientation in terms of the orientation quaternion com-
C ... ponents Q0, Q1, Q2 and Q3. Then the transformation is calculated similarly
C ... to as in FINEUL. Notice that we do not need to know the "ICASE" now.
C ... Reference: Finflo Report F-115, Equation (5.14)
C ... Matti Palin, April 2020

      USE CONSTANTS
      
      IMPLICIT NONE

      REAL :: Q0, Q1, Q2, Q3
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XO, YO, ZO, XN, YN, ZN
      REAL :: XFM, YFM, ZFM, XONEW, YONEW, ZONEW

C ... Initialize matrix elements
      A11 = 0.0; A12 = 0.0; A13 = 0.0
      A21 = 0.0; A22 = 0.0; A23 = 0.0
      A31 = 0.0; A32 = 0.0; A33 = 0.0
     
C ... Transformation matrix. Alternative formulations are in comments.
      A11 = Q0**2 + Q1**2 - Q2**2 - Q3**2
      A12 = 2.0*(Q1*Q2 - Q0*Q3)
      A13 = 2.0*(Q1*Q3 + Q0*Q2)

      A21 = 4.0*Q1*Q2 - A12			!A21 = 2.0*(Q1*Q2 + Q0*Q3)
      A22 = Q0**2 - Q1**2 + Q2**2 - Q3**2
      A23 = 2.0*(Q2*Q3 - Q0*Q1)

      A31 = 4.0*(Q1*Q3) - A13		!A31 = 2.0*(Q1*Q3 - Q0*Q2)
      A32 = 4.0*(Q2*Q3) - A23		!A32 = 2.0*(Q2*Q3 + Q0*Q1)
      A33 = Q0**2 - Q1**2 - Q2**2 + Q3**2

C ... From FINFLO coordinate system to the flight mechanics 
C ... coordinate system.
      XFM = -XO
      YFM = -ZO 
      ZFM = -YO

C ... Round small values to zero in order to prevent errors with machine precision
      IF (ABS(A11) < EPS12) A11 = 0.0
      IF (ABS(A12) < EPS12) A12 = 0.0
      IF (ABS(A13) < EPS12) A13 = 0.0

      IF (ABS(A21) < EPS12) A21 = 0.0
      IF (ABS(A22) < EPS12) A22 = 0.0
      IF (ABS(A23) < EPS12) A23 = 0.0

      IF (ABS(A31) < EPS12) A31 = 0.0
      IF (ABS(A32) < EPS12) A32 = 0.0
      IF (ABS(A33) < EPS12) A33 = 0.0

C ... Rotations according to the orientation quaternion.
      XN = A11*XFM + A21*YFM + A31*ZFM
      YN = A12*XFM + A22*YFM + A32*ZFM
      ZN = A13*XFM + A23*YFM + A33*ZFM

      XONEW = XN		! Matti Palin - testi!
      YONEW = YN
      ZONEW = ZN

      RETURN
      END SUBROUTINE FINEULQua
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EULFINQua(Q0,Q1,Q2,Q3,XO,YO,ZO,XONEW,YONEW,ZONEW)

      USE CONSTANTS
      
      IMPLICIT NONE

      REAL :: Q0, Q1, Q2, Q3
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XO, YO, ZO
      REAL :: XFM, YFM, ZFM, XONEW, YONEW, ZONEW   

C ... Transformation matrix.
      A11 = Q0**2 + Q1**2 - Q2**2 - Q3**2
      A21 = 2.0*(Q1*Q2 - Q0*Q3)
      A31 = 2.0*(Q1*Q3 + Q0*Q2)

      A12 = 2.0*(Q1*Q2 + Q0*Q3)
      A22 = Q0**2 - Q1**2 + Q2**2 - Q3**2
      A32 = 2.0*(Q2*Q3 - Q0*Q1)

      A13 = 2.0*(Q1*Q3 - Q0*Q2)
      A23 = 2.0*(Q2*Q3 + Q0*Q1)
      A33 = Q0**2 - Q1**2 - Q2**2 + Q3**2

C ... Round small values to zero in order to prevent errors with machine precision
      IF (ABS(A11) < EPS12) A11 = 0.0
      IF (ABS(A12) < EPS12) A12 = 0.0
      IF (ABS(A13) < EPS12) A13 = 0.0

      IF (ABS(A21) < EPS12) A21 = 0.0
      IF (ABS(A22) < EPS12) A22 = 0.0
      IF (ABS(A23) < EPS12) A23 = 0.0

      IF (ABS(A31) < EPS12) A31 = 0.0
      IF (ABS(A32) < EPS12) A32 = 0.0
      IF (ABS(A33) < EPS12) A33 = 0.0

C ... From FINFLO coordinate system to the flight mechanics 
C ... coordinate system.
      XFM = XO
      YFM = YO 
      ZFM = ZO


C ... Rotations according to the orientation quaternion.
      XONEW = A11*XFM + A21*YFM + A31*ZFM		! Matti Palin - testi
      YONEW = A12*XFM + A22*YFM + A32*ZFM
      ZONEW = A13*XFM + A23*YFM + A33*ZFM

      RETURN
      END SUBROUTINE EULFINQua
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EULQUA(PSI,THETA,PHI,Q0,Q1,Q2,Q3,ICASE)

C ... Converts Euler angles to quaternion representation
C ... Based on Finflo report F-115, Table 5.1. Matti Palin, July 2019.

C ... Q0, Q1, Q2 and Q3 are the quaternion components and they form a unit quaternion
C ... that represents an orientation in three dimensions. The conversion
C ... between Euler angles and quaternions can be derived by multiplying
C ... an orientation quaternion from the left by simple rotation quaternions
C ... that correspond to rotations around the axes in three dimensions.
C ... This idea is also presented by P.H. Zipfel.

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE
      REAL :: PSI, THETA, PHI
      REAL :: Q0, Q1, Q2, Q3
      REAL :: S1, S2, S3, C1, C2, C3

      S1 = SIN(PHI  /2.0)
      S2 = SIN(THETA/2.0)
      S3 = SIN(PSI  /2.0)

      C1 = COS(PHI  /2.0)
      C2 = COS(THETA/2.0)
      C3 = COS(PSI  /2.0)

      SELECT CASE(ICASE)

      CASE(1)     ! Rotation case xyz

      Q0 = C3*C2*C1 + S3*S2*S1
      Q1 = C3*C2*S1 - S3*S2*C1
      Q2 = S3*C2*S1 + C3*S2*C1
      Q3 = S3*C2*C1 - C3*S2*S1

      CASE(2)     ! Rotation case yxz

      Q0 = C3*C2*C1 - S3*S2*S1
      Q1 = C3*S2*C1 - S3*C2*S1
      Q2 = S3*S2*C1 + C3*C2*S1
      Q3 = S3*C2*C1 + C3*S2*S1

      CASE(3)     ! Rotation case xzy

      Q0 = C3*C2*C1 - S3*S2*S1
      Q1 = C3*C3*S1 + S3*S2*C1
      Q2 = S3*C2*C1 + C3*S2*S1
      Q3 = C3*S2*C1 - S3*C2*S1

      CASE(4)     ! Rotation case yzx

      Q0 = C3*C2*C1 + S3*S2*S1
      Q1 = S3*C2*C1 - C3*S2*S1
      Q2 = C3*C2*S1 - S3*S2*C1
      Q3 = S3*C2*S1 + C3*S2*C1

      CASE(5)     ! Rotation case zyx

      Q0 = C3*C2*C1 - S3*S2*S1
      Q1 = S3*C2*C1 + C3*S2*S1
      Q2 = C3*S2*C1 - S3*C2*S1
      Q3 = S3*S2*C1 + C3*C2*S1

      CASE(6)     ! Rotation case zxy

      Q0 = C3*C2*C1 + S3*S2*S1
      Q1 = C3*S2*C1 + S3*C2*S1
      Q2 = S3*C2*C1 - C3*S2*S1
      Q3 = C3*C2*S1 - S3*S2*C1

      END SELECT

C ... Convention requires that Q0 be non-negative, and therefore it is possible
C ... that the sign has to be changed. The sign of a quaternion does not affect
C ... its use otherwise.

      IF (Q0<0) THEN
      Q0 = -Q0
      Q1 = -Q1
      Q2 = -Q2
      Q3 = -Q3
      END IF

      END SUBROUTINE EULQUA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE QUAEUL(PSI,THETA,PHI,Q0,Q1,Q2,Q3,ICASE)
C ... Converts a unit quaternion (an orientation representation)
C ... into Euler angles. Reference: Finflo report F-115.
C ... This is the inverse calculation of EULQUA, and the formulas
C ... can be extracted from those equations by inspection.
C ... In reality they were found using a trial-and-error method :)
C ... Matti Palin, July 2019.

      USE CONSTANTS

      IMPLICIT NONE

      INTEGER :: ICASE
      REAL :: PSI, THETA, PHI
      REAL :: Q0, Q1, Q2, Q3

      SELECT CASE(ICASE)

      CASE(1)     ! Rotation case xyz

      PHI   = ATAN2( 2*(Q0*Q1 + Q2*Q3)  , Q0**2 - Q1**2 - Q2**2 + Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q2 - Q1*Q3), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q3 + Q1*Q2)  , Q0**2 + Q1**2 - Q2**2 - Q3**2)

      CASE(2)     ! Rotation case yxz

      PHI   = ATAN2( 2*(Q0*Q2 - Q1*Q3)  , Q0**2 - Q1**2 - Q2**2 + Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q1 + Q2*Q3), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q3 - Q1*Q2)  , Q0**2 - Q1**2 + Q2**2 - Q3**2)

      CASE(3)     ! Rotation case xzy

      PHI   = ATAN2( 2*(Q0*Q1 - Q2*Q3)  , Q0**2 - Q1**2 + Q2**2 - Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q3 + Q1*Q2), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q2 - Q1*Q3)  , Q0**2 + Q1**2 - Q2**2 - Q3**2)

      CASE(4)     ! Rotation case yzx

      PHI   = ATAN2( 2*(Q0*Q2 + Q1*Q3)  , Q0**2 + Q1**2 - Q2**2 - Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q3 - Q1*Q2), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q1 + Q2*Q3)  , Q0**2 - Q1**2 + Q2**2 - Q3**2)

      CASE(5)     ! Rotation case zyx

      PHI   = ATAN2( 2*(Q0*Q3 - Q1*Q2)  , Q0**2 + Q1**2 - Q2**2 - Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q2 + Q1*Q3), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q1 - Q2*Q3)  , Q0**2 - Q1**2 - Q2**2 + Q3**2)


      CASE(6)     ! Rotation case zxy

      PHI   = ATAN2( 2*(Q0*Q3 + Q1*Q2)  , Q0**2 - Q1**2 + Q2**2 - Q3**2)
      THETA = ASIN(  MAX(MIN(2*(Q0*Q1 - Q2*Q3), 1.0),-1.0) )
      PSI   = ATAN2( 2*(Q0*Q2 + Q1*Q3)  , Q0**2 - Q1**2 - Q2**2 + Q3**2)

      END SELECT

      END SUBROUTINE QUAEUL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      function roundnominal(num)

      REAL :: roundnominal, num, absnum
      INTEGER :: magnitude
      INTEGER(KIND=8) :: largeint, multiplier

      if (num == 0.0) then
            roundnominal = 0.0
            return
      end if
      absnum = ABS(num)
      magnitude = FLOOR(LOG10(absnum))
      multiplier = 10.0**(10-magnitude)
      largeint = NINT(absnum*multiplier)
      roundnominal = SIGN(largeint/multiplier,INT(num,8))

C      WRITE (*,'(A20,E28.15)') "ROUNDN     num = ",num	! Matti Palin, September 2020
C      WRITE (*,'(A20,I12)') "ROUNDNmagnitude= ",magnitude
C      WRITE (*,'(A20,E28.15)') "ROUNDmultiplier= ",multiplier
C      WRITE (*,'(A20,E28.15)') "ROUND      res = ",roundnominal
C      WRITE (*,*) " -----"

      return
      end function
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
