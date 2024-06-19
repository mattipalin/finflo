C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPS(NBL,M,LUSGS,RO,RM,RN,RW,E,P,PDIFF,U,V,W,C,CP,
     + DRDP,DRDH,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,
     + D1,D2,D3,DTL,DT,GAMMA,FRSDEN,FRSVEL,FRSPRE,DRO,DM,DN,DW,DE,
     + VOL,VIS,EPS2,VIST,ITURB,KOVER,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,
     + IT,IL,IK,CFL,THETA,IDI1,IDI2,IDI3,NTOT,PR,PRT,OHMI,
     + RK,REPS,DRK,DEPS,SRK,SEPS,PTUR,XC,YC,ZC,CSIMPS,
     + ICP,JCP,KCP,LAMIN,JRDIS,JRPRE,UROT,VROT,WROT,
     + OMEGA,OMEX,OMEY,OMEZ,TURLIM,
     + RKLIM,EPSLIM,FI,DFI,SFI,MAXSB,NSCAL,JRIMP,TIMEL,IPRESC,PSEUCO,
     + RKMAX,ABU,ABV,ABW,DPCORR,A1PP,A1RM,A1RN,A1RW,RMASS,ICON,NPATCH,
     + IFLUX,INTERI,INTERJ,INTERK,ITERM,TOLER,ROFOR,RMFOR,RNFOR,RWFOR,
     + EFOR,PDFOR,RKFOR,REFOR,RKSI,FIFOR,NCHIM,
     + ZZZ,MAXW,MA2,IDIS,XMASSB,PRO,VAR,BLKS,MULPHL,
     + TRANSL,TRM,CDIFF,CDIFFT)

      USE CHARACTERS
      USE TYPE_ARRAYS
      USE INTEGERS, ONLY : IREPEA,NREPEA
      USE NS3CO, ONLY : IN, JN, KN, IC9

      INTEGER :: IPHASE

      REAL :: RO(*),RM(*),RN(*),RW(*),E(*),P(*),U(*),V(*),W(*),C(*),
     +     CP(*),A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),
     +     A2XA(*),A2YA(*),A1ZA(*),DRO(*),DM(*),DN(*),DE(*),DW(*),
     +     VOL(*),DTL(*),VIS(*),EPS2(*),VIST(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),RK(*),REPS(*),DRK(*),DEPS(*),
     +     SRK(*),SEPS(*),PTUR(*),OHMI(*),D1(*),D2(*),D3(*),
     +     UROT(*),VROT(*),WROT(*),DRDP(*),DRDH(*),
     +     FI(MAXSB,MAX(1,NSCAL)),DFI(MAXSB,MAX(1,NSCAL)),
     +     SFI(MAXSB,MAX(1,NSCAL)),PDIFF(*),
     +     ABU(*),ABV(*),ABW(*),A1RM(*),A1RN(*),A1RW(*),DPCORR(*),
     +     A1PP(*),RMASS(*),
     +     ROFOR(*),RMFOR(*),RNFOR(*),RWFOR(*),EFOR(*),PDFOR(*),
     +     RKFOR(*),REFOR(*),RKSI(*),FIFOR(MAXSB,MAX(1,NSCAL)),ZZZ(*),
     +     XMASSB(4)

      REAL :: CDIFF, CDIFFT

      REAL, ALLOCATABLE :: DXI(:),XI(:),EVAPR(:),ROG(:),
     +                     DGI(:),GII(:),DRETI(:),RETI(:)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(BLOCKS)           :: BLKS(*)
      TYPE(INTERMITTENCY)     :: TRM(*)

      REAL :: XC(*), YC(*), ZC(*)

C ... LAST LINES ARE FOR THE PRESSURE CORRECTION TEST. APU MAY BE DANGEROUS!
      LOGICAL :: TIMEL,MULPHL,TRANSL

      INTEGER :: ICP(*), JCP(*), KCP(*), ICON(IC9,*)

      DTURB = 1. ! MAGNIFICATION OF THE TURBULENCE QUATITIES

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN

      ISTR = 1
      JSTR = (IMAX + 2*IN)
      KSTR = (JMAX + 2*JN)*JSTR

C ... ADJUST CELLS WITH FRICTION TO THE GRID LEVEL USED

      IDM1    = (IDI1-1)/2**(M-1) + 1
      IDM2    = (IDI2-1)/2**(M-1) + 1
      IDM3    = (IDI3-1)/2**(M-1) + 1

      LAMIN1  = LAMIN/100
      LAMIN2  = (LAMIN - LAMIN1*100)/10
      LAMIN3  =  LAMIN - LAMIN1*100 - LAMIN2*10
      IF (ITURB <= 23) JRIMP = 0
      KSCAL   = NSCAL
      IF (ITURB >= 24 .AND. JRIMP == 1) KSCAL = NSCAL - 6
      KBEGIN  = NSCAL - KSCAL + 1 ! BEGINING OF SCALARS

      IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT') THEN
         ALLOCATE(DXI(NTOT),XI(NTOT),EVAPR(NTOT),ROG(NTOT))
         DXI(1:NTOT) = VAR(1:NTOT)%DX(2)
      ELSE
         ALLOCATE(DXI(1))
      ENDIF ! CAVIT

      IF(TRANSL) THEN
         ALLOCATE(GII(NTOT),RETI(NTOT),DGI(NTOT),DRETI(NTOT))
         DGI(1:NTOT)   = TRM(1:NTOT)%DG
         DRETI(1:NTOT) = TRM(1:NTOT)%DRET
      ELSE
         ALLOCATE(GII(1),RETI(1),DGI(1),DRETI(1))
      ENDIF ! TRANSL

C
C ... IMPLICIT SWEEPS
C
C-------------------------------------------------------------------
       
      IF (NCHIM == 0) THEN
C ... LOCAL TIME-STEPPING. DIFF. STEP: PVISC   = 1.8    !*8.**(M-1)

      IF(LUSGS == 0) PVISC = 1.8     !  TESTED 1.5.1990
      IF(LUSGS == 1) PVISC = 1.8     !  WILL BE TESTED SOON

      CALL TIME3(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDM1,IDM2,IDM3,M,TIMEL,
     + PVISC,ABV,1,FRSPRE,CDIFF)
      ENDIF

C ... Total mass imbalance in the block
      IF(M == 1) CALL SUMDRO(XMASSB(1),DRO,VOL,IMAX,JMAX,KMAX)

      IF(ITURB <= 2 .OR. ITURB == 8) THEN
         CALL LOCDT(DRO,DM,DN,DW,DE,DTL,NTOT,DT)
      ELSE IF(ITURB >= 3 .AND. ITURB/=8) THEN
         CALL LOCDTK(DRO,DM,DN,DW,DE,DRK,DEPS,DTL,NTOT,DT,DTURB,
     +   BLKS(NBL)%SOLUTION_TYPE,DXI,TRANSL,DGI,DRETI)
      ENDIF
      IF(NSCAL >= 1) THEN
         CALL LOCDTF(DFI,DTL,NTOT,DT,MAXSB,NSCAL)
      ENDIF

C ... Total mass imbalance after the local time stepping
      IF(M == 1) CALL SUMDRO(XMASSB(2),DRO,VOL,IMAX,JMAX,KMAX)

      IF(ICYCLE == IPRINT) THEN
C ... nan trap
         CALL NATRAP(DTL,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dTL ',
     +   IREPEA(8),NREPEA(8))
         CALL NATRAP(DRO,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dRO ',
     +   IREPEA(8),NREPEA(8))
         CALL NATRAP(DM,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dRM ',
     +   IREPEA(8),NREPEA(8))
         CALL NATRAP(DN,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dRN ',
     +   IREPEA(8),NREPEA(8))
         CALL NATRAP(DW,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dRW ',
     +   IREPEA(8),NREPEA(8))
         CALL NATRAP(DE ,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dE  ',
     +   IREPEA(8),NREPEA(8))
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER LOCAL TIME STEPPING'
         CALL PRINYS(3,TIC,  DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DMC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DNC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DEC,  DE,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT')
     +   CALL PRINYS(3,CHAR_VAR(11,2),DXI,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +   NBL,M)
         CALL PRINYS(3,UROTC,UROT,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,VROTC,VROT,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,WROTC,WROT,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)

         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF

         IF(TRANSL) THEN
           CALL PRINYS(3,DGC,DGI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
           CALL PRINYS(3,DRETC,DRETI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF		

C ... FOR SCALAR EQ PPR 24.2
         DO 8180 NS = 1,NSCAL
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NBL,M)
 8180    CONTINUE
      ENDIF

C ... EXPLICIT WITH A SMALL COURANT NUMBER

      IF(CFL  < .05) RETURN

C-----SOURCES-----------------------------------------------------------

      IF(OMEGA /= 0. .AND. .NOT.TIMEL) THEN
         CALL ROTDIA(DM,DN,DW,DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT)
         IF(ITURB >= 21) THEN
            CALL ROTRSM(DFI(1,1),DFI(1,2),DFI(1,3),DFI(1,4),DFI(1,5),
     +           DFI(1,6),DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT) 
         ENDIF
         IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)'                         '
            WRITE(3,*)' RESULTS AFTER ROTATIONAL CORRECTION'
            CALL PRINYS(3,DMC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DNC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF
      ENDIF

      IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT') THEN
         XI(1:NTOT)    = VAR(1:NTOT)%X(2)
         EVAPR(1:NTOT) = VAR(1:NTOT)%EVAPR(2)
         ROG(1:NTOT)   = PRO(1:NTOT)%RO(2)
c         CALL EVAPIM(DXI,RO,DTL,VAR(1:NTOT)%X(2),VAR(1:NTOT)%EVAPR(2),
c     +   NTOT)
cccccccc         CALL EVAPIM(DXI,RO,ROG,DTL,XI,EVAPR,NTOT)
         IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)'                         '
            WRITE(3,*)' RESULTS AFTER THE EVAP LINEARIZATION'
            CALL PRINYS(3,CHAR_VAR(11,2),DXI,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +      NBL,M)
         ENDIF
      ENDIF ! CAVIT

      IF(TRANSL) THEN
         GII(1:NTOT)   = TRM(1:NTOT)%G
         RETI(1:NTOT)  = TRM(1:NTOT)%RET
      ENDIF

C ... LINEARIZATION OF THE k-epsilon SOURCE IS PREVENTED IN LU-SGS
C ... DUE TO THE STRONG INTERNAL DISSIPATION

      IF(IDIS ==  1 .AND. LUSGS. EQ. 0) THEN
         CALL TURBCO(ITURB,JRDIS,JRPRE,CE1,CE2,C3,C21,CMU,CTA,PSIGK,
     +        PSIGE,AA1,ETA0)
         IF(ITURB /= 8 .OR. ITURB /= 9)
     +   CALL SOURIM(DRK,DEPS,RK,REPS,E,SRK,SEPS,PTUR,VOL,DTL,
     +        NTOT,M,RO,VIS,CSIMPS,RKMAX,DTURB,CE1,CE2)
         IF(ITURB == 9)
     +   CALL SPALIM(DRK,DEPS,RK,REPS,E,SRK,SEPS,PTUR,VOL,DTL,
     +        NTOT,M,RO,VIS,CSIMPS,RKMAX,DTURB,CE1,CE2,EPSLIM)
      ELSEIF(IDIS == 2) THEN
         CALL IMSOKO(DRK,DEPS,RK,REPS,PTUR,EPS2,VIS,DTL,RO,U,V,W,E,
     &        GAMMA,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,CSIMPS,
     &        KOVER)
         IF(TRANSL) THEN
              CALL IMSOTR(TRM,DGI,DRETI,DTL,NTOT,
     &           IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN)
         ENDIF
      ENDIF

      IF(ITURB == 9) DRK(1:NTOT) = 0. ! In Spalart-Allmaras model

      IF (ITURB >= 23) THEN
         CALL SOUREY(DFI,FI,E,SFI,VOL,DTL,NTOT,RO,VIS,MAXSB)
      ENDIF

C ... COMMENT THE FOLLOWING IF-BLOCK FOR TEST-STAGE 2
c      IF(TIMEL) THEN          
c          CALL DERIMP(DRO,DM,DN,DW,DE,DTL,DT,NTOT)
c          IF(ITURB >= 3) CALL DERIKE(DRK,DEPS,DTL,DT,NTOT)
c          IF(NSCAL > 0) CALL DERISC(DFI,DTL,DT,NTOT,MAXSB,NSCAL)
c      ENDIF
         IF(TRANSL) THEN
            IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)
            WRITE(3,*) ' Results after transition-source linearization'
              CALL PRINYS(3,DGC,DGI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
              CALL PRINYS(3,DRETC,DRETI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,
     &        NBL,M)
            ENDIF
         ENDIF		

C-----------------------------------------------------------------------

C ... FROM THE CONSERVATIVE TO THE PRIMITIVE VARIABLES (IN STATE EQ.)

      IF (KSCAL >= 1) THEN
          CALL COPRFI(DFI(1,KBEGIN),FI(1,KBEGIN),RO,DRO,NTOT,KSCAL,
     +    MAXSB)
      ENDIF
      IF(ITURB <= 2 .OR. ITURB == 8) THEN
         CALL COPR(DRO,DM,DN,DW,DE,RO,E,P,PDIFF,U,V,W,C,DRDP,DRDH,
     +   FRSVEL,NTOT)
      ELSE IF(JRIMP == 0) THEN
         CALL COPRKE(DRO,DM,DN,DW,DE,DRK,DEPS,RO,E,P,PDIFF,U,V,W,C,RK,
     +   REPS,DRDP,DRDH,FRSVEL,BLKS(NBL)%SOLUTION_TYPE,PRO,VAR,DXI,DTL,
     +   NTOT,TRANSL,GII,RETI,DGI,DRETI,PDFOR,RKSI,VOL,NBL,IMAX,JMAX,
     +   KMAX)
      ELSE IF(JRIMP == 1) THEN
         CALL COPRRE(DRO,DM,DN,DW,DE,DRK,DEPS,DFI,RO,E,U,V,W,
     +   RK,REPS,FI,NTOT,MAXSB)
      ELSE
         WRITE(*,*) 'Unlegal ITURB value. ITURB=',ITURB,'. Exiting...'
         STOP
      ENDIF

C ... Total change in pressure inside the block
      IF(M == 1) CALL SUMDRO(XMASSB(3),DE,VOL,IMAX,JMAX,KMAX)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' AFTER UPDATING SOURCE TERM AND TO THE PRIMITIVE'
         CALL PRINYS(3,DHC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DUC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DVC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DPC,  DE,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)

         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF

         IF(TRANSL) THEN
           CALL PRINYS(3,DGC,DGI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
           CALL PRINYS(3,DRETC,DRETI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF		

         IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT')
     +   CALL PRINYS(3,CHAR_VAR(11,2),DXI,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +   NBL,M)
         DO NS = 1,NSCAL
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NBL,M)
         ENDDO
      ENDIF

C **********************************************************************
      IF(LUSGS == 0) THEN   !  DDADI METHOD
C **********************************************************************
C
C ... THE FIRST SWEEP IS MADE IN K-DIRECTION
C
      CALL IMPSWE(DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,U,V,W,
     +     C,RO,P,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,DRK,DEPS,
     +     RK,REPS,SRK,SEPS,DRDP,DRDH,PR,PRT,KCP,IN,JN,KN,IPRESC,
     +     WROT,DFI,FI,MAXSB,NSCAL,KSCAL,KBEGIN,ITURB,0,JRIMP,
     +     DTURB,ZZZ,MAXW,ABU,DXI,BLKS(NBL)%SOLUTION_TYPE,
     +     TRANSL,DGI,DRETI,NBL,M)
             
      IF(ICYCLE == IPRINT) THEN
           WRITE(3,*) 
           WRITE(3,*) 'after k-dir'
           CALL PRINYS(3,DPC,  DE,  IT,IL,2,IK,JMAX,KMAX,IMAX,NBL,M)
           CALL PRINYS(3,DUC,  DM,  IT,IL,2,IK,JMAX,KMAX,IMAX,NBL,M)
           CALL PRINYS(3,DVC,  DN,  IT,IL,2,IK,JMAX,KMAX,IMAX,NBL,M)
           CALL PRINYS(3,DWC,  DW,  IT,IL,2,IK,JMAX,KMAX,IMAX,NBL,M)
      ENDIF
C
C ... THE SECOND SWEEP IS MADE IN I-DIRECTION IN THE CASE OF THE
C ... BIDIAGONAL METHOD. THE ARRAYS ARE SORTED FROM IJK TO JKI
C ... TO OBTAIN A STRIDE = 1. 
C
      CALL SORT38(A1,A1XA,A1YA,A1ZA,UROT,ISTRID,JSTRID,KSTRID,1,ZZZ)
      CALL IMPSWE(DRO,DM,DN,DW,DE,A1,VOL,A1XA,A1YA,A1ZA,DTL,U,V,W,
     +     C,RO,P,JMAX,KMAX,IMAX,THETA,DT,IDM1,VIS,VIST,DRK,DEPS,
     +     RK,REPS,SRK,SEPS,DRDP,DRDH,PR,PRT,KCP,JN,KN,IN,IPRESC,
     +     UROT,DFI,FI,MAXSB,NSCAL,KSCAL,KBEGIN,ITURB,1,JRIMP,
     +     DTURB,ZZZ,MAXW,ABU,DXI,BLKS(NBL)%SOLUTION_TYPE,
     +     TRANSL,DGI,DRETI,NBL,M)
      IF(ICYCLE == IPRINT) THEN
           WRITE(3,*) 
           WRITE(3,*) 'after i-dir'
           CALL PRINYS(3,DPC,  DE,  IT,IL,1,IK,KMAX,IMAX,JMAX,NBL,M)
           CALL PRINYS(3,DUC,  DM,  IT,IL,1,IK,KMAX,IMAX,JMAX,NBL,M)
           CALL PRINYS(3,DVC,  DN,  IT,IL,1,IK,KMAX,IMAX,JMAX,NBL,M)
           CALL PRINYS(3,DWC,  DW,  IT,IL,1,IK,KMAX,IMAX,JMAX,NBL,M)
      ENDIF
C
C ... THE THIRD SWEEP IS MADE IN THE J-DIRECTION FOR THE
C ... BIDIAGONAL SOLUTION. 
C
      CALL SORT38(A2,A2XA,A2YA,A2ZA,VROT,ISTRID,JSTRID,KSTRID,2,ZZZ)
      CALL IMPSWE(DRO,DM,DN,DW,DE,A2,VOL,A2XA,A2YA,A2ZA,DTL,U,V,W,
     +     C,RO,P,KMAX,IMAX,JMAX,THETA,DT,IDM2,VIS,VIST,DRK,DEPS,
     +     RK,REPS,SRK,SEPS,DRDP,DRDH,PR,PRT,KCP,KN,IN,JN,IPRESC,
     +     VROT,DFI,FI,MAXSB,NSCAL,KSCAL,KBEGIN,ITURB,2,JRIMP,
     +     DTURB,ZZZ,MAXW,ABU,DXI,BLKS(NBL)%SOLUTION_TYPE,
     +     TRANSL,DGI,DRETI,NBL,M)
      IF(ICYCLE == IPRINT) THEN
           WRITE(3,*) 
           WRITE(3,*) 'after j-dir'
           CALL PRINYS(3,DPC,  DE,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF

C **********************************************************************
      ELSE IF(LUSGS == 1) THEN      ! LU-SGS METHOD
C **********************************************************************
         IF(IMAX >= 996) STOP 'IROW is too small'
         CALL DOMAW(MAXW,IMAX,JMAX,IN,JN,MAW)
         IF(ITURB <= 2 .OR. ITURB == 8) THEN
            CALL PREDGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,IMAX,
     +           JMAX,KMAX,THETA,DT,IDM1,IDM2,IDM3,VIS,VIST,PR,
     +           PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,ZZZ,MAW,MAXW)
            CALL DIAGGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,IMAX,
     +           JMAX,KMAX,THETA,DT,IDM1,IDM2,IDM3,VIS,VIST,PR,
     +           PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,ZZZ,MAW,MAXW)
            CALL CORRGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,IMAX,
     +           JMAX,KMAX,THETA,DT,IDM1,IDM2,IDM3,VIS,VIST,PR,
     +           PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,ZZZ,MAW,MAXW)
          ELSEIF(JRIMP == 0) THEN
            CALL NATRAP(DEPS,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' depe  ',
     +           IREPEA(8),NREPEA(8))
            CALL PRKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,DRDP,DRDH,
     +           IMAX,JMAX,KMAX,THETA,DT,IDM1,IDM2,IDM3,VIS,VIST,
     +           DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,IN,JN,KN,UROT,VROT,
     +           WROT,RKSI,ZZZ,MAW,MAXW,TRANSL,DGI,DRETI)
            IF(ICYCLE == IPRINT) THEN
            WRITE(3,*) 
            WRITE(3,*) 'after PREDGS'
            CALL PRINYS(3,DPC,  DE,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DUC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DVC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DGC, DGI,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            ENDIF
            CALL NATRAP(DEPS,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dep1  ',
     +           IREPEA(8),NREPEA(8))
            CALL DIKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,DRDP,DRDH,
     +           IMAX,JMAX,KMAX,THETA,VIS,VIST,
     +           DRK,DEPS,RK,REPS,SRK,SEPS,IN,JN,KN,UROT,VROT,WROT,
     +           RKSI,TRANSL,DGI,DRETI)
            IF(ICYCLE == IPRINT) THEN
            WRITE(3,*) 
            WRITE(3,*) 'after DIKEGS'
            CALL PRINYS(3,DGC, DGI,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            ENDIF
            CALL NATRAP(DEPS,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dep2  ',
     +           IREPEA(8),NREPEA(8))
            CALL COKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,
     +           A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,DRDP,DRDH,
     +           IMAX,JMAX,KMAX,THETA,DT,IDM1,IDM2,IDM3,VIS,VIST,
     +           DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,IN,JN,KN,UROT,VROT,
     +           WROT,RKSI,ZZZ,MAW,MAXW,TRANSL,DGI,DRETI)
            CALL NATRAP(DEPS,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dep3  ',
     +           IREPEA(8),NREPEA(8))
            IF(ICYCLE == IPRINT) THEN
            WRITE(3,*) 
            WRITE(3,*) 'after COKEGS'
            CALL PRINYS(3,DGC, DGI,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            ENDIF


      ENDIF
C **********************************************************************
      ENDIF ! LUSGS == 0
C **********************************************************************
        
      IF(IPRESC == 1) THEN
         IF(ICYCLE == IPRINT) THEN
            WRITE(92,*)' RESULTS BEFORE PRESSURE CORRECTION'
            CALL PRINYS(92,DUC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(92,DNC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(92,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(92,DHC,  DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(92,DPC,  DE,  IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF

C ... IS THIS THE CORRECT PLACE FOR THE PRESSURE CORRECTION?
c        write(6,*) icycle,iprint,interk,iflux,iterm ! oikeellisuustsekki
      CALL PRECOR(NBL,M,RO,RM,RN,RW,E,P,PDIFF,U,V,W,C,DM,DN,DW,DE,DRDP,
     + A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,VOL,VIS,
     + VIST,DTL,ABU,ABV,ABW,DPCORR,RMASS,A1PP,A1RM,A1RN,A1RW,ICP,
     + JCP,KCP,IMAX,JMAX,KMAX,ICYCLE,INTERI,INTERJ,INTERK,JBOT,JTOP,
     + IDM1,IDM2,IDM3,IT,IL,IK,IPRINT,ICON,NPATCH,IBOT,ITOP,KBOT,KTOP,
     + IFLUX,UROT,VROT,WROT,ITURB,FRSDEN,FRSPRE,ITERM,TOLER,ZZZ,MAXW)       

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' PRESSURE CHANGE AFTER THE PRESSURE CORRECTION'
         CALL PRINYS(3,DPC,DE,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF

      ENDIF  ! IPRESC == 1

C **********************************************************************

C ... FROM THE PRIMITIVE TO THE CONSERVATIVE VARIABLES 
C ... is removed concerning the main variables 
C ...(primitive variables DeltaP and DeltaT are used)

C ... Maximum change in pressure after the implicit sweeps
      IF(M == 1) CALL SUMDRO(XMASSB(4),DE,VOL,IMAX,JMAX,KMAX)

C ... Change from dh to dtemp

      CALL DIVV12(DRO,DRO,CP,NTOT)

C ... Temporary change from dtemp to phase-variables

      IF(MULPHL) THEN
      DO IPHASE = 1,NPHASES
c      CALL SETVIN(VAR(1:NTOT)%DTEMP(IPHASE),DRO,IMAX,JMAX,KMAX,IN,JN,KN)
c      CALL SETVIN(VAR(1:NTOT)%DTEMP(IPHASE),DRO,IMAX,JMAX,KMAX,IN,JN,KN)
      VAR(1:NTOT)%DTEMP(IPHASE)=DRO(1:NTOT)
      ENDDO
      ENDIF

C ... Put the turbulence and scalar quantities to conservative form

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL MULV12(DRK, RO,NTOT)
         CALL MULV12(DEPS,RO,NTOT)
      ENDIF
      IF(ITURB == 9) DRK(1:NTOT) = 0.
      IF(TRANSL) THEN
         CALL MULV12(DGI,  RO,NTOT)
         CALL MULV12(DRETI,RO,NTOT)
      ENDIF
      IF(NSCAL >= 1) THEN
         DO NS = KBEGIN,KBEGIN + KSCAL - 1
            CALL MULV12(DFI(1,NS),RO,NTOT)
         ENDDO
      ENDIF

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER THE IMPLICIT SWEEPS'
         CALL PRINYS(3,DTEMPC,DRO,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
         ZZZ(1:NTOT) = VAR(1:NTOT)%DTEMP(IPHASE)
         CALL PRINYS(3,CHAR_VAR(4,IPHASE),ZZZ(1),
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDDO
         ENDIF
         CALL PRINYS(3,DUC,    DM,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DVC,    DN,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DWC,    DW,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DPC,    DE,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF ! ITURB >= 0

         IF(TRANSL) THEN
           CALL PRINYS(3,DGC,DGI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
           CALL PRINYS(3,DRETC,DRETI(1),IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         ENDIF		

         IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT') THEN
            CALL PRINYS(3,CHAR_VAR(11,2),DXI,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +      NBL,M)
            DO IPHASE = 1,NPHASES
            ZZZ(1:NTOT) = VAR(1:NTOT)%DTEMP(IPHASE)
            CALL PRINYS(3,CHAR_PH(2,IPHASE),ZZZ(1),IT,IL,0,IK,
     +      IMAX,JMAX,KMAX,NBL,M)
            ENDDO
         ENDIF
C ... FOR SCALAR EQ PPR 24.2
         IF (NSCAL > 0) THEN
            DO 8181 NS = 1,NSCAL
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NBL,M)
 8181       CONTINUE
         ENDIF ! NSCAL > 0
      ENDIF ! ICYCLE == IPRINT

      IF(BLKS(NBL)%SOLUTION_TYPE == 'CAVIT') THEN
         VAR(1:NTOT)%DX(1) =-DXI(1:NTOT)
         VAR(1:NTOT)%DX(2) = DXI(1:NTOT)
         DO IPHASE = 1,NPHASES
         VAR(1:NTOT)%DTEMP(IPHASE)=0.
         ENDDO
         DEALLOCATE(DXI,XI,EVAPR,ROG)
      ENDIF

      IF(TRANSL) THEN
         TRM(1:NTOT)%DG   = DGI(1:NTOT)
         TRM(1:NTOT)%DRET = DRETI(1:NTOT)
         DEALLOCATE(GII,RETI,DGI,DRETI)
      ENDIF

      RETURN
      END SUBROUTINE IMPS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPSWE(DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +     U,V,W,C,RO,P,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,DRK,
     +     DEPS,RK,REPS,SRK,SEPS,DRDP,DRDH,PR,PRT,KCP,IN,JN,KN,IPRESC,
     +     WROT,DFI,FI,MAXSB,NSCAL,KSCAL,KBEGIN,ITURB,IDIR,JRIMP,
     +     DTURB,ZZZ,MAXW,ABU,DXI,SOLUTION_TYPE,TRANSL,DGI,DRETI,NBL,M)

      USE INTEGERS, ONLY : IREPEA,NREPEA

      REAL ::RO(*),P(*),U(*),V(*),C(*),DRO(*),DM(*),DN(*),DE(*),VOL(*),
     +     DTL(*),VIS(*),VIST(*),W(*),DW(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),RK(*),REPS(*),DRK(*),DEPS(*),
     +     WROT(*),DFI(MAXSB,NSCAL),FI(MAXSB,NSCAL),
     +     SRK(*),SEPS(*),DRDP(*),DRDH(*),ZZZ(*),ABU(*),DXI(*),
     +     DGI(*),DRETI(*)

      INTEGER :: KCP(*)
      LOGICAL :: TRANSL
      CHARACTER(LEN=10) :: SOLUTION_TYPE

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      NTOT    = ISTRID*JSTRID*KSTRID

      CALL DOMAW(MAXW,IMAX,JMAX,IN,JN,MAW)

      IF(KMAX > 1 ) THEN
      IF(IPRESC /= 1)THEN
         IF(ITURB <= 2 .OR. ITURB == 8) THEN
            CALL PRED (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,
     +           VIST,PR,PRT,IN,JN,KN,WROT,ZZZ,MAW,MAXW)
            CALL DIAG (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,
     +           VIST,PR,PRT,IN,JN,KN,WROT)
            CALL CORR (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,
     +           VIST,PR,PRT,IN,JN,KN,WROT,ZZZ,MAW,MAXW)
         ELSEIF(JRIMP == 0) THEN
            CALL PREDKE (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DT,IDM3,
     +           VIS,VIST,DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,
     +           IN,JN,KN,WROT,DTURB,ZZZ,MAW,MAXW)       
                 CALL NATRAP(DE,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dp1  ',
     +           IREPEA(8),NREPEA(8))

            CALL DIAGKE (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DT,IDM3,
     +           VIS,VIST,DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,
     +           IN,JN,KN,WROT,DTURB,ZZZ,MAW,MAXW)
                 CALL NATRAP(DE,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dp2  ',
     +           IREPEA(8),NREPEA(8))
            CALL CORRKE (DRO,DM,DN,DW,DE,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DT,IDM3,
     +           VIS,VIST,DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,
     +           IN,JN,KN,WROT,DTURB,ZZZ,MAW,MAXW)
                 CALL NATRAP(DE,IMAX,JMAX,KMAX,IN,JN,KN,NBL,M,' dp3  ',
     +           IREPEA(8),NREPEA(8))
         ELSEIF(JRIMP == 1) THEN
            CALL PREDRE(DRO,DM,DN,DW,DE,DRK,DEPS,
     +           DFI(1,1),DFI(1,2),DFI(1,3),DFI(1,4),DFI(1,5),DFI(1,6),
     +           A3,VOL,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,RK,REPS,SRK,
     +           SEPS,FI(1,1),FI(1,2),FI(1,3),FI(1,4),FI(1,5),FI(1,6),
     +           IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,DRDP,DRDH,
     +           PR,PRT,IN,JN,KN,WROT,ZZZ,MAW,MAXW)
            CALL DIAGRE(DRO,DM,DN,DW,DE,DRK,DEPS,
     +           DFI(1,1),DFI(1,2),DFI(1,3),DFI(1,4),DFI(1,5),DFI(1,6),
     +           A3,VOL,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,RK,REPS,SRK,
     +           SEPS,FI(1,1),FI(1,2),FI(1,3),FI(1,4),FI(1,5),FI(1,6),
     +           IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,DRDP,DRDH,
     +           PR,PRT,IN,JN,KN,WROT,ZZZ,MAW,MAXW)
            CALL CORRRE(DRO,DM,DN,DW,DE,DRK,DEPS,
     +           DFI(1,1),DFI(1,2),DFI(1,3),DFI(1,4),DFI(1,5),DFI(1,6),
     +           A3,VOL,A3XA,A3YA,A3ZA,DTL,U,V,W,C,RO,P,RK,REPS,SRK,
     +           SEPS,FI(1,1),FI(1,2),FI(1,3),FI(1,4),FI(1,5),FI(1,6),
     +           IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,DRDP,DRDH,
     +           PR,PRT,IN,JN,KN,WROT,ZZZ,MAW,MAXW)
         ENDIF
C  ...   FOR TRANSITION VARIABLES
         IF(TRANSL) THEN
            CALL PREDFI (DGI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,1) 
            CALL DIAGFI (DGI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DGI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,1)

            CALL PREDFI (DRETI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,2)
            CALL DIAGFI (DRETI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DRETI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,2)
         ENDIF
      ENDIF ! IPRESC /= 1

C ... PRESSURE CORRECTION TEST

      IF(IPRESC == 1) THEN

            CALL PREDFI (DM,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DM,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DM,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,
     +           VIST,PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)

            CALL PREDFI (DN,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DN,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DN,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)

            CALL PREDFI (DW,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DW,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DW,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)

            CALL PREDFI (DRO,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DRO,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,NTOT)
            CALL CORRFI (DRO,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,NTOT,ZZZ,MAW,MAXW,0)
      ENDIF ! IPRESC == 1

C ... FOR SCALARS PPR 24.2
      IF(KSCAL >= 1) THEN
            CALL PREDFI (DFI(1,KBEGIN),A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,KSCAL,MAXSB,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DFI(1,KBEGIN),A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,KSCAL,MAXSB)
            CALL CORRFI (DFI(1,KBEGIN),A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,IDM3,VIS,VIST,
     +           PR,PRT,IN,JN,KN,KSCAL,MAXSB,ZZZ,MAW,MAXW,0)
      ENDIF
C ... FOR CAVITATING FLOW TSII 28.6       
          
      IF(SOLUTION_TYPE == 'CAVIT') THEN       
            CALL PREDFI (DXI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,0,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,MAXW,ZZZ,MAW,MAXW,0)
            CALL DIAGFI (DXI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,IMAX,JMAX,KMAX,THETA,IN,JN,KN,1,MAXW)
            CALL CORRFI (DXI,A3,VOL,A3XA,A3YA,A3ZA,DTL,
     +           U,V,W,RO,IMAX,JMAX,KMAX,THETA,DT,0,VIS,VIST,
     +           PR,PRT,IN,JN,KN,1,MAXW,ZZZ,MAW,MAXW,0)
      ENDIF

      ENDIF ! KMAX > 1
C
C ... FROM ??? TO IJK
C
1000  CONTINUE

      IF(IDIR /= 0) THEN
        CALL SORT38(A3,A3XA,A3YA,A3ZA,WROT,ISTRID,JSTRID,KSTRID,3-IDIR,
     +        ZZZ)
      ENDIF
C ... SORT FOR NEXT SWEEP IJK TO JKI
      CALL SORTZZ(VOL,ISTRID,JSTRID,KSTRID,1,ZZZ)
C      CALL SORTZZ(P,ISTRID,JSTRID,KSTRID,1,ZZZ)
      CALL SORT2(DTL,DRO,DM,DN,DW,DE,U,V,W,C,RO,VIS,VIST,
     2     ISTRID,JSTRID,KSTRID,1,ZZZ)
      CALL SORT50(DRDP,DRDH,ISTRID,JSTRID,KSTRID,1,ZZZ)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL SORT39(DRK,DEPS,RK,REPS,ISTRID,JSTRID,KSTRID,1,ZZZ)
C         CALL SORT50(SRK,SEPS,ISTRID,JSTRID,KSTRID,1,ZZZ)
      ENDIF
      IF (NSCAL > 0) THEN
         CALL SORTSS(DFI,ISTRID,JSTRID,KSTRID,1,MAXSB,NSCAL,ZZZ)
         CALL SORTSS( FI,ISTRID,JSTRID,KSTRID,1,MAXSB,NSCAL,ZZZ)
      ENDIF

      IF(SOLUTION_TYPE == 'CAVIT') THEN
         CALL SORTZZ(DXI,ISTRID,JSTRID,KSTRID,1,ZZZ)
      ENDIF

      IF(TRANSL) THEN
         CALL SORTZZ(DGI,  ISTRID,JSTRID,KSTRID,1,ZZZ)
         CALL SORTZZ(DRETI,ISTRID,JSTRID,KSTRID,1,ZZZ)
      ENDIF

      RETURN
      END SUBROUTINE IMPSWE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOCDT(DRO,DM,DN,DW,DE,DTL,NTOT,DTG)

      DIMENSION :: DRO(*), DM(*), DN(*), DE(*), DTL(*), DW(*)
C
C ... USE THE LOCAL TIME-STEP BEFORE THE IMPLICIT SWEEP
C
c      PDT     = 1./DTG
      DO 1100 K = 1,NTOT
      DT1     = DTL(K)
      DRO(K ) = DT1*DRO(K)
      DM(K )  = DT1*DM(K)
      DN(K )  = DT1*DN(K)
      DW(K )  = DT1*DW(K)
      DE(K )  = DT1*DE(K)
1100  CONTINUE

      RETURN
      END SUBROUTINE LOCDT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOCDTK(DRO,DM,DN,DW,DE,DRK,DEPS,DTL,NTOT,DTG,DTURB,
     + SOLUTION_TYPE,DXI,TRANSL,DGI,DRETI)

      REAL :: DRO(*),DM(*),DN(*),DE(*),DTL(*),DW(*),DRK(*),DEPS(*),
     +        DXI(*),DGI(*),DRETI(*)

      LOGICAL :: TRANSL

      CHARACTER(LEN=10) :: SOLUTION_TYPE
C
C ... USE THE LOCAL TIME-STEP BEFORE THE IMPLICIT SWEEP
C
c      PDT     = 1./DTG

      DO 1100 K = 1,NTOT
      DT1     = DTL(K)
      DRO(K)  = DT1*DRO(K)
      DM(K)   = DT1*DM(K)
      DN(K)   = DT1*DN(K)
      DW(K)   = DT1*DW(K)
      DE(K)   = DT1*DE(K)
      DRK(K)  = DTURB*DT1*DRK(K)
      DEPS(K) = DTURB*DT1*DEPS(K)
      IF(SOLUTION_TYPE == 'CAVIT') THEN
         DXI(K) = DT1*DXI(K)
      ENDIF
      IF(TRANSL) THEN
         DGI(K)   = DT1*DGI(K)
         DRETI(K) = DT1*DRETI(K)
      ENDIF
1100  CONTINUE

      RETURN
      END SUBROUTINE LOCDTK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DERIMP(DRO,DM,DN,DW,DE,DTL,DT,NTOT)

      REAL :: DRO(*), DM(*), DN(*), DW(*), DE(*), DTL(*)

      EPS     = 1.E-20

      DO 1000 N = 1,NTOT
      PVOL    = 1./(1.+1.5*DTL(N)/DT)
      DRO(N)  = DRO(N)*PVOL
      DM(N)   = DM(N) *PVOL
      DN(N)   = DN(N) *PVOL
      DW(N)   = DW(N) *PVOL
      DE(N)   = DE(N) *PVOL
1000  CONTINUE

      RETURN
      END SUBROUTINE DERIMP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DERIKE(DRK,DEPS,DTL,DT,NTOT)

      REAL :: DRK(*),DEPS(*),DTL(*)

      EPS     = 1.E-20

      DO 1000 N = 1,NTOT
      PVOL    = 1./(1.+1.5*DTL(N)/DT)
      DRK(N)  = DRK(N) *PVOL
      DEPS(N) = DEPS(N)*PVOL
1000  CONTINUE

      RETURN
      END SUBROUTINE DERIKE

C *******************************************************************
C     ROUTINES FOR ALGEBRAIC MODELS                                 *
C *******************************************************************
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COPR(DRO,DM,DN,DW,DE,RO,E,P,PDIFF,U,V,W,C,DRDP,DRDH,
     + FRSVEL,NTOT)

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),RO(*),U(*),V(*),W(*),
     + P(*),E(*),DRDP(*),DRDH(*),C(*),PDIFF(*)
C
C ... FROM THE CONSERVATIVE TO THE PRIMITIVE VARIABLES (dh,du,dv,dw,dp)
C ... OUTPUT: DRO = DH AND DE = DP
C

      DO 1000 I = 1,NTOT

      VEL2    = (U(I)**2 + V(I)**2 + W(I)**2)
      BETA    = .5*VEL2
      DRODH   = DRDH(I)   
      DRODP   = DRDP(I) 

      DRODPM  = 1./DRODP
      YPRO    = 1./RO(I)

      DRI     = DRO(I)
      DMI     = DM(I)
      DNI     = DN(I)
      DWI     = DW(I)
      DRE     = DE(I)
      DKIN    = U(I)*DMI + V(I)*DNI + W(I)*DWI

      DM(I)   = YPRO*(-U(I)*DRI + DMI)
      DN(I)   = YPRO*(-V(I)*DRI + DNI)
      DW(I)   = YPRO*(-W(I)*DRI + DWI)

      HI      = YPRO*(E(I) + P(I)) 
      DEI     = DRE - (HI-2*BETA)*DRI - DKIN
      DPI     = C(I)**2*(DRI - DRODH*YPRO*DEI)
      DHI     = (DRI - DRODP*DPI)/DRODH
      DRO(I)  = DHI
      DE(I)   = DPI

1000  CONTINUE

      RETURN
      END SUBROUTINE COPR


C *********************************************************************
C     K-EPSILON ROUTINES                                              *
C *********************************************************************
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SOURIM(DRK,DEPS,RK,REPS,E,SRK,SEPS,PTUR,VOL,DTL,NTOT,
     + M,RO,VIS,C1,RKMAX,DTURB,CE1,CE2)

      REAL :: DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),VOL(*),DTL(*),
     + PTUR(*),RO(*),VIS(*),E(*)

      EPS     = 1.E-20
      RKLIM   = 0.

      DO 1000 N = 1,NTOT
      PVOL    = DTURB*DTL(N)
      RKN     = 1./(RK(N)+EPS)
      ELIM1   = 1./(ABS(0.1*E(N)-RK(N))+EPS)
      ELIM2   = 1./(ABS(RKLIM-RK(N))+EPS)
      ELIM    = MAX(ELIM1,ELIM2)
      SOURCE  = C1*PVOL*ABS(PTUR(N))*MAX(RKN,ELIM)
      DRK(N)  = DRK(N)/(1. + 2.*PVOL*REPS(N)*RKN + 2.*SOURCE )
c      DEPS(N) = (DEPS(N) + 0.*2.*1.92*PVOL*(REPS(N)*RKN)**2*DRK(N))
      DEPS(N) = DEPS(N)
     +          /(1.+ 2.*CE2*PVOL*REPS(N)*RKN + CE1*SOURCE)
1000  CONTINUE

      RETURN
      END SUBROUTINE SOURIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMSOKO(DRK,DOM,RK,ROMEGA,PTUR,EPS2,VIS,DTL,RO,U,V,W,
     &     E,GAMMA,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN,
     &     CSIMPS,KOVER)

C     Treat the source terms of the k-omega models.
C                                                   AH 11.2.1997

      REAL :: DRK(*),DOM(*),RK(*),ROMEGA(*),PTUR(*),EPS2(*),VIS(*),
     &     DTL(*),RO(*),U(*),V(*),W(*),E(*)

      EPS  = 1.0E-10

C     Model coefficients:
      CALL COEFKO(6,A1KLEB,BSTAR,CKAPPA,
     &     BETA1,SIGRK1,SIGOM1,BETA2,SIGRK2,SIGOM2,
     &     SIGPHI,ZETA2,KOVER)
      GAMKOM = BETA1/BSTAR - SIGOM1*CKAPPA**2/(SQRT(BSTAR))
      CKMAX  = 1.0/(2.0*CSIMPS)
      COMAX  = 1.0/(2.0*CSIMPS)

      DO K = 1,KMAX
         IA = (KN+K-1)*KSTR + JN*JSTR + IN*ISTR
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX    ! J1 = J-1
            I  = IJ - J1*IMAX
            L  = IA + J1*JSTR + I*ISTR
C           Maximum allowable changes:
C            EMEANF = P(L)/(GAMMA-1.) 
C     &             + 0.5*RO(L)*(U(L)**2 + V(L)**2 + W(L)**2)
            EMEANF = E(L) - RK(L)
            DRKMAX = MAX(CKMAX*MIN(RK(L),ABS(0.1*EMEANF-0.9*RK(L))),
     &           EPS)
            DOMMAX = COMAX*ROMEGA(L)
            PROL   = 1./RO(L)
            PROL2  = 2.0*PROL
            ATPROD = ABS(PTUR(L))
            RNUTRB = MAX((EPS2(L)-1.0),EPS)*VIS(L)*PROL
            DRKDEN = 1.0 + DTL(L)*(ATPROD/DRKMAX+PROL*BSTAR*ROMEGA(L))
            DOMDEN = 1.0 + DTL(L)*(GAMKOM*ATPROD/(RNUTRB*DOMMAX) 
     &             + 2.*BETA1*ROMEGA(L))
C    &             + PROL2*BETA1*ROMEGA(L)) ! Aeron mersu

            DRK(L) = DRK(L)/DRKDEN
            DOM(L) = DOM(L)/DOMDEN
         END DO
      END DO

      RETURN
      END SUBROUTINE IMSOKO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMSOTR(TRM,DGI,DRETI,DTL,NTOT,
     &           IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN)
C ... Source term linearization for the transition equations
C     when IPRESC = 0

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: NTOT,N,K,IA,IJ,J1,I,
     &           IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN

      REAL :: DGI(*),DRETI(*),DTL(*),EPSLIM,DGN,DRETN

      TYPE(INTERMITTENCY) TRM(*)

      EPSLIM  = 1.0E-10
      DO K = 1,KMAX
         IA = (KN+K-1)*KSTR + JN*JSTR + IN*ISTR
         DO IJ = 1,IMAX*JMAX
            J1 = (IJ-1)/IMAX
            I  = IJ - J1*IMAX
            N  = IA + J1*JSTR + I*ISTR
            DGN      = 1.0 + DTL(N)*(
     &             ABS( TRM(N)%SG/(2.*0.1+EPSLIM) )
     &           + ABS( TRM(N)%LG ) )
            DRETN    = 1.0 + DTL(N)*(
     &             ABS( TRM(N)%SRET/(2000.*0.1+EPSLIM) )
     &           + ABS( TRM(N)%LRET ) )
            DGI(N)   = DGI(N)  /DGN
            DRETI(N) = DRETI(N)/DRETN
         ENDDO
      ENDDO
      END SUBROUTINE IMSOTR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SPALIM(DRK,DEPS,RK,REPS,E,SRK,SEPS,PTUR,VOL,DTL,NTOT,
     + M,RO,VIS,C1,RKMAX,DTURB,CE1,CE2,EPSLIM)

      REAL :: DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),VOL(*),DTL(*),
     + PTUR(*),RO(*),VIS(*),E(*)

      EPS     = 1.E-20
      RKLIM   = 0.

      DO 1000 N = 1,NTOT
      PVOL    = DTURB*DTL(N)
      ELIM2   = 1./(ABS(EPSLIM-REPS(N))+EPS)
      DRK(N)  = 0.
      DEPS(N) = DEPS(N)
     + /(1.+ PVOL*(2.*ABS(SRK(N))/(REPS(N)+EPS) + ABS(PTUR(N))*ELIM2))
1000  CONTINUE

      RETURN
      END SUBROUTINE SPALIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EVAPIM(DXI,RO,ROG,DTL,XI,EVAPG,NTOT)

      IMPLICIT NONE

      INTEGER :: NTOT,N
      REAL    :: EPS,XLIM,ALIM
      REAL    :: DXI(*),RO(*),DTL(*),XI(*),EVAPG(*),ROG(*)

      EPS     = 1.E-9
      ALIM    = 0.001

      DO 1000 N = 1,NTOT
c      XLIM    = 0.001*XI(N) + EPS
      XLIM    = ALIM*ROG(N)/RO(N) + EPS
      DXI(N)  = DXI(N)/(1.+DTL(N)*ABS(EVAPG(N))/(RO(N)*XLIM))
1000  CONTINUE

      RETURN
      END SUBROUTINE EVAPIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COPRKE(DRO,DM,DN,DW,DE,DRK,DEPS,RO,E,P,PDIFF,U,V,W,C,
     + RK,REPS,DRDP,DRDH,FRSVEL,SOLUTION_TYPE,PRO,VAR,DXI,DTL,NTOT,
     + TRANSL,GII,RETI,DGI,DRETI,PDFOR,RKSI,VOL,NGL,IMAX,JMAX,KMAX)
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: NTOT,I,NGL,IMAX,JMAX,KMAX
      REAL :: FRSVEL,VEL2,DRODH,DRODP,YPRO,RKI,DRI,DMI,DNI,DWI,DRE,
     + DKI,DKIN,HI,DEI,DPI,DHI,DPMAX,RMUL,ALIM,EPS,XLIM,DXII,DXIA,DPIA
      
      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),RO(*),U(*),V(*),W(*),E(*),
     + DRK(*),DEPS(*),RK(*),REPS(*),DRDP(*),DRDH(*),P(*),C(*),PDIFF(*),
     + DXI(*),DTL(*),GII(*),RETI(*),DGI(*),DRETI(*),PDFOR(*),RKSI(*),
     + VOL(*)

      LOGICAL :: TRANSL

      CHARACTER(LEN=10) :: SOLUTION_TYPE

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
C
C ... FROM THE CONSERVATIVE TO THE PRIMITIVE VARIABLES (dh,dV,dp,dk,deps)
C ...                                              (dro,drv,de,drk,dreps)
C ... THIS SUBROUTINE IS USED FOR K-EPSILON
C ... OUTPUT: DRO = DH AND DE = DP
C
      EPS     = 1.E-9
      ALIM    = 1.E-2!0.001 ! viimeisin
        
      DO 1000 I = 1,NTOT

      VEL2    = (U(I)**2 + V(I)**2 + W(I)**2)
c      BETA    = .5*VEL2
      DRODH   = DRDH(I)   
      DRODP   = DRDP(I) 

      YPRO    = 1./RO(I)
      RKI     = RK(I)*YPRO

      DRI     = DRO(I)
      DMI     = DM(I)
      DNI     = DN(I)
      DWI     = DW(I)
      DRE     = DE(I)
      DKI     = DRK(I)
      DKIN    = U(I)*DMI + V(I)*DNI + W(I)*DWI + DKI

      DM(I)   = YPRO*(-U(I)*DRO(I) + DMI)
      DN(I)   = YPRO*(-V(I)*DRO(I) + DNI)
      DW(I)   = YPRO*(-W(I)*DRO(I) + DWI)
      DRK(I)  = YPRO*(-RKI *DRO(I) + DKI)
      DEPS(I) = YPRO*(-REPS(I)*YPRO*DRO(I) + DEPS(I))

      IF(TRANSL) THEN
         DGI(I)   = YPRO*(-GII(I) *DRO(I)*YPRO + DGI(I)  )
         DRETI(I) = YPRO*(-RETI(I)*DRO(I)*YPRO + DRETI(I))  
      ENDIF

      IF(SOLUTION_TYPE == 'FLUID') THEN ! Checked 7.2.2018
         HI      = YPRO*(E(I) + P(I)) 
         DEI     = DRE - (HI-VEL2-RKI)*DRI - DKIN
         DPI     = C(I)**2*(DRI - DRODH*YPRO*DEI)
         DPIA    = DPI

C         IF(RKSI(I)  > 0. .AND. DPI /=  0) THEN
Cc          call ijkpai(I,imax,jmax,kmax,ii,jj,kk)
Cc          DPIA = DPI
C           DPI  = RKSI(I)*(PDFOR(I)-PDIFF(I))
CC         ENDIF
         DHI     = (DRI - DRODP*DPI)/DRODH
      ELSE IF(SOLUTION_TYPE == 'CAVIT') THEN
         DXIA    = DXI(I)
         XLIM    = ALIM*MAX(1.E-3,MIN(VAR(I)%ALFA(2),1.-VAR(I)%ALFA(2)))
     +             *PRO(I)%RO(2)/RO(I) + EPS
         XLIM    = MIN(0.0001,XLIM) !0.00001
         DXI(I)  = YPRO*(DXI(I) - VAR(I)%X(2)*DRI)
         DXII    = DXI(I)*RO(I)
         DPMAX  = 1.*ABS(PDIFF(I)-PRO(I)%DPSAT) + .01
         DXI(I)  = DXI(I)/(1.+DTL(I)*ABS(VAR(I)%EVAPR(2))/(RO(I)*XLIM))
         DPMAX  = .1*ABS(PDIFF(I)-PRO(I)%DPSAT) + 1.E-2 ! oli .25 3.1.2008

        RMUL    = 1./(RO(I)*C(I)**2) + ABS(VAR(I)%EVAPH(1))*DTL(I)/!c0111 oli evapr(2) 18.1.2008
     +             DPMAX*(1./PRO(I)%RO(2) - 1./PRO(I)%RO(1)) ! y.o. väärin!
        RMUL  = 1./RMUL
          DPI     = RMUL*(DRI  /PRO(I)%RO(1) ! 2809
     +             +  DXIA ! DXII  ! 2809
     +             *(1./PRO(I)%RO(2) - 1./PRO(I)%RO(1)))
         DHI     = 0.
      ENDIF

      DRO(I)  = DHI
      DE(I)   = DPI
      DE(I)   = DPI + .667*DKI
1000  CONTINUE

      RETURN
      END SUBROUTINE COPRKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPFRE(NGL,M,PRC,
     + IMAX,JMAX,KMAX,INRE,JNRE,KNRE,NTOT,ICYCLE,IPRINT,
     + IT,IL,IK,NPATCH,ICON,IHF,SURLE,DSURLE,VOL,RO,WAVEH,IWAVEB,
     + DWMAX,FLODWH,WHMAX,WHMIN,A1,A2,A3,XCP,YCP,ZCP,FREDIN,PDIFF,
     + POLD,D1,D2,D3,DTL,CHLREF,SUMDWH,CFLFRE,CFL,CFLL,DTWMIN,DTWMAX,
     + DWMV,U,V,W,F1RM,F1RN,F1RW,FPRINTL,IGRID)

      USE CHARACTERS
      USE TYPE_ARRAYS
      USE NS3CO, ONLY : LN, IC9

      IMPLICIT NONE

      INTEGER :: NGL,M,IMAX,JMAX,KMAX,NTOT,ICYCLE,IPRINT,IT,IL,IK,L,
     +   NPATCH,IF2,IA,IM,JA,JM,IFACE,I,J,K,IPL,IPG,NP,NR,
     +   ISTRID,JSTRID,KSTRID,ISTR,JSTR,KSTR,KA,IJ,IP,
     +   ISLAB,JSLAB,IS,ILPO,KOKO,MGRIDA,MCYCLE,ITERMA,ITERHA,
     +   ICONV3,IMAXSL,JMAXSL,KMAXSL,NFL,IF1,IDWMAX,JDWMAX,IGRID,
     +   IN,JN,KN,INRE,JNRE,KNRE

      INTEGER :: ICON(IC9,*),IHF(*),IWAVEB(*)

      REAL :: DWMAX,SUMDWH,WHMAX,WHMIN,DWH,FREDIF,DISTI,DISTJ,DISTK,
     +   DISTIP,DISTJP,DISTKP,APUA,FREDIN,DTLOC,CHLREF,FLODWH,
     +   CFLFRE,CFL,CFLL,DTWMIN,DTWMAX,DWMV,CFLMUL

      REAL :: SURLE(*),DSURLE(*),VOL(*),RO(*),A1(*),A2(*),A3(*),
     +   PDIFF(*),POLD(*),D1(*),D2(*),D3(*),DTL(*),
     +   U(*),V(*),W(*),F1RM(*),F1RN(*),F1RW(*)

      REAL :: WAVEH(*),XCP(*),YCP(*),ZCP(*)

      TYPE(PRE_COR) :: PRC(*)

      REAL, ALLOCATABLE :: AE(:),AW(:),AN(:),AS(:),AT(:),AD(:),AP(:),
     +   DFRE(:),APU(:)

      LOGICAL :: FPRINTL

C ... A routine for the patchwise free-surface solution

      IN = INRE
      JN = JNRE
      KN = KNRE

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      KSTRID  = ISTRID*JSTRID

      FREDIF  = 1.*FREDIN             ! Add implicit damping ??
      CFLMUL  = CFLFRE/CFL            ! Currently from local DT
      IF(M > 1) CFLMUL = CFLFRE/CFLL

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE FREE-SURFACE SOLUTION'
         WRITE(3,*)'                         '
         F1RM(1:NTOT) =  PRC(1:NTOT)%F1R                    
         F1RN(1:NTOT) =  PRC(1:NTOT)%F2R                    
         F1RW(1:NTOT) =  PRC(1:NTOT)%F3R                    
         CALL PRINYS(3,F1RC,F1RM(1), IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         CALL PRINYS(3,F2RC,F1RN(1), IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         CALL PRINYS(3,F3RC,F1RW(1), IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
      ENDIF
       
      DO IP = 1,NPATCH

C ------------------------------------------------------------------
      IF(ICON(1,IP) == 13) THEN
C ------------------------------------------------------------------

      IFACE   = ICON(3,IP)
      IPL     = ICON(2,IP)            ! Proces local patch number
      IPG     = ICON(25,IP)           ! Global patch number (unused)
      IF2     = IHF(IPL)              ! Starting address of this patch

      IA      = ICON(4,IP)         
      IM      = ICON(5,IP)         
      JA      = ICON(6,IP) 
      JM      = ICON(7,IP) 

      IF(IFACE == 2 .OR. IFACE == 5) THEN
        JA    = ICON(4,IP)
        JM    = ICON(5,IP) 
        IA    = ICON(6,IP) 
        IM    = ICON(7,IP) 
      ENDIF


C ... Slab size

      IMAXSL  = IM - IA + 1
      JMAXSL  = JM - JA + 1
      KMAXSL  = 1

      ISLAB   = IMAXSL + 2*IN
      JSLAB   = JMAXSL + 2*JN
      KOKO    = (1+2*KN)*ISLAB*JSLAB
         
      ALLOCATE (AP(KOKO),AE(KOKO),AW(KOKO),AN(KOKO),AS(KOKO),
     +   AT(KOKO),AD(KOKO),DFRE(KOKO),APU(KOKO))

      AP = 1.
      AE = 0.
      AW = 0.
      AN = 0.
      AS = 0.
      AT = 0.
      AD = 0.
      DFRE = 0.

C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = KSTRID
         KSTR = 1
         K    = 1
         NFL  = 0 
         IF(IFACE == 4) K   = IMAX 
         IF(IFACE == 4) NFL = 1

C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = KSTRID
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         NFL  = 0
         IF(IFACE == 5) K   = JMAX
         IF(IFACE == 5) NFL = ISTRID

C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = KSTRID
         K    = 1
         NFL  = 0
         IF(IFACE == 6) K   = KMAX
         IF(IFACE == 6) NFL = KSTRID 
      ENDIF

      KA      = (KN+K-1)*KSTR
        
***********************************************************************
C ... Start to establish the implicit stage

      DO J = JA,JM

         IJ   = (JN+J-1)*JSTR + KA
*         NR   = (JN+J-JA)*ISLAB - IA + IN + IF2
         NR   = (LN+J-JA)*(IM-IA+1+2*LN) - IA + LN + IF2
         ILPO = KN*ISLAB*JSLAB + (JN-JA+J)*ISLAB + IN-IA+1

         DO I = IA,IM
           L     = 1  + (IN+I-1)*ISTR + IJ ! CELL INDEX (STRIDE=ISTR)
           NP  = I + NR                       ! Patch index with LN ghost cells
           IF1 = L + NFL                      ! Face and flux index
           IS    = ILPO + I                ! Local patch index

           DFRE(IS)  = -DSURLE(NP)/(VOL(L)*RO(L)) !nyt miinus AMG:lle
           APU(IS)   = DFRE(IS)
           DSURLE(NP)= DFRE(IS)                   !Explicit change (overwritten)
 
           IF(IFACE == 1 .OR. IFACE == 4) THEN
              DISTI  = (D2(L-ISTRID)+D2(L)) + 1.E-20
              DISTIP = (D2(L+ISTRID)+D2(L)) + 1.E-20
              DISTJ  = (D3(L-KSTRID)+D3(L)) + 1.E-20
              DISTJP = (D3(L+KSTRID)+D3(L)) + 1.E-20
              DISTK  = A3(L)       *(RO(L)+RO(L-KSTRID))*FREDIF/DISTJ
              DISTKP = A3(L+KSTRID)*(RO(L)+RO(L+KSTRID))*FREDIF/DISTJP
              DISTJ  = A2(L)       *(RO(L)+RO(L-ISTRID))*FREDIF/DISTI
              DISTJP = A2(L+ISTRID)*(RO(L)+RO(L+ISTRID))*FREDIF/DISTIP
              DISTI  = 0.
              DISTIP = 0.
           ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
              DISTI  = (D3(L-KSTRID)+D3(L)) + 1.E-20
              DISTIP = (D3(L+KSTRID)+D3(L)) + 1.E-20
              DISTJ  = (D1(L-1)     +D1(L)) + 1.E-20
              DISTJP = (D1(L+1)     +D1(L)) + 1.E-20
              DISTK  = A3(L)       *(RO(L)+RO(L-KSTRID))*FREDIF/DISTI
              DISTKP = A3(L+KSTRID)*(RO(L)+RO(L+KSTRID))*FREDIF/DISTIP
              DISTI  = A1(L)       *(RO(L)+RO(L-1))     *FREDIF/DISTJ
              DISTIP = A1(L+1)     *(RO(L)+RO(L+1))     *FREDIF/DISTJP
              DISTJ  = 0.
              DISTJP = 0.
           ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
              DISTI  = (D1(L-1)     +D1(L)) + 1.E-20
              DISTIP = (D1(L+1)     +D1(L)) + 1.E-20
              DISTJ  = (D2(L-ISTRID)+D2(L)) + 1.E-20
              DISTJP = (D2(L+ISTRID)+D2(L)) + 1.E-20
              DISTI  = A1(L)       *(RO(L)+RO(L-1))     *FREDIF/DISTI
              DISTIP = A1(L+1)     *(RO(L)+RO(L+1))     *FREDIF/DISTIP
              DISTJ  = A2(L)       *(RO(L)+RO(L-ISTRID))*FREDIF/DISTJ
              DISTJP = A2(L+ISTRID)*(RO(L)+RO(L+ISTRID))*FREDIF/DISTJP
              DISTK  = 0.
              DISTKP = 0.
           ENDIF

           AN(IS) = (-MAX(0.,-PRC(L+ISTRID)%F2R))
           AS(IS) = (-MAX(0., PRC(L)%F2R)       )
           AE(IS) = (-MAX(0.,-PRC(L+1)%F1R)     )
           AW(IS) = (-MAX(0., PRC(L)%F1R)       )
           AT(IS) = (-MAX(0.,-PRC(L+KSTRID)%F3R))
           AD(IS) = (-MAX(0., PRC(L)%F3R)       )
           APUA   =  -AN(IS) -AS(IS)-AE(IS)-AW(IS)-AT(IS)-AD(IS)

c           VELOC  = SQRT(U(L)**2+V(L)**2+W(L)**2+1.E-20)
c           DTLOC  = 0.3/VELOC
c            go to 110

C ... Boundaries in KSI-direction
      
            IF(I == IA) THEN
            IF(IWAVEB(L-ISTR) == 4  .OR. IWAVEB(L-ISTR) == 5  .OR. 
     +         IWAVEB(L-ISTR) == 8  .OR. IWAVEB(L-ISTR) == 9  .OR. 
     +         IWAVEB(L-ISTR) == 10 .OR. IWAVEB(L-ISTR) == 15)
     +         THEN
               AE(IS-1)=-1.
            ELSE
               AE(IS-1)= 0.
            ENDIF
            ENDIF
            IF(I == IM) THEN
            IF(IWAVEB(L+ISTR) == 4  .OR. IWAVEB(L+ISTR) == 5  .OR. 
     +         IWAVEB(L+ISTR) == 8  .OR. IWAVEB(L+ISTR) == 9  .OR. 
     +         IWAVEB(L+ISTR) == 10 .OR. IWAVEB(L+ISTR) == 15)
     +         THEN
               AW(IS+1)=-1.
            ELSE
               AW(IS+1)= 0.
            ENDIF
            ENDIF

C ... Boundaries in ETA-direction

            IF(J == JA) THEN
            IF(IWAVEB(L-JSTR) == 4  .OR. IWAVEB(L-JSTR) == 5  .OR. 
     +         IWAVEB(L-JSTR) == 8  .OR. IWAVEB(L-JSTR) == 9  .OR. 
     +         IWAVEB(L-JSTR) == 10 .OR. IWAVEB(L-JSTR) == 15)
     +         THEN
               AN(IS-ISLAB)=-1.
            ELSE
               AN(IS-ISLAB)= 0.
            ENDIF
            ENDIF
            IF(J == JM) THEN
            IF(IWAVEB(L+JSTR) == 4  .OR. IWAVEB(L+JSTR) == 5  .OR. 
     +         IWAVEB(L+JSTR) == 8  .OR. IWAVEB(L+JSTR) == 9  .OR. 
     +         IWAVEB(L+JSTR) == 10 .OR. IWAVEB(L+JSTR) == 15)
     +         THEN
               AS(IS+ISLAB)=-1.
            ELSE
               AS(IS+ISLAB)= 0.
            ENDIF
            ENDIF
c110    continue

           DTLOC  = CFLMUL*DTL(L)
           DTLOC  = MAX(DTLOC,DTWMIN*CHLREF*.25) ! Calibrated to CHLREF = 4 m
           DTLOC  = MIN(DTLOC,DTWMAX*CHLREF*.25) ! Calibrated to CHLREF = 4 m

           APUA   = APUA  + VOL(L)*RO(L)/DTLOC !0001 Not divided by zero
           AP(IS) = APUA  !* 20. !+ ABS(DFRE(IS))/0.0001 Under-relaxation

           AN(IS) = (AN(IS) - DISTJP)/(VOL(L)*RO(L))
           AS(IS) = (AS(IS) - DISTJ )/(VOL(L)*RO(L))
           AE(IS) = (AE(IS) - DISTI )/(VOL(L)*RO(L))
           AW(IS) = (AW(IS) - DISTIP)/(VOL(L)*RO(L))
           AT(IS) = (AT(IS) - DISTKP)/(VOL(L)*RO(L))
           AD(IS) = (AD(IS) - DISTK )/(VOL(L)*RO(L))
           AP(IS) = (AP(IS) + DISTJP+DISTJ+DISTIP+DISTI+DISTKP+DISTK)/
     +              (VOL(L)*RO(L))
         ENDDO
      ENDDO

***********************************************************************

C ... Jacobian computed. Solution with an algebraic multigrid

      MCYCLE = 1
      ITERMA = 10
      ITERHA = 0
      MGRIDA = 1
      ICONV3 = 0

      DO J = JA-1,JM+1
         IJ   = (JN+J-1)*JSTR + KA
*         NR = (JN+J-JA)*ISLAB - IA + IN + IF2
         NR = (LN+J-JA)*(IM-IA+1+2*LN) - IA + LN + IF2
         ILPO = KN*ISLAB*JSLAB + (JN-JA+J)*ISLAB + IN-IA+1
         DO I = IA-1,IM+1
           L     = 1  + (IN+I-1)*ISTR + IJ ! CELL INDEX (STRIDE=ISTR)
           NP  = I + NR                       ! Patch index with LN ghost cells
           IF1 = L + NFL                      ! Face and flux index
           IS    = ILPO + I                ! Local patch index
           AN(IS)= AN(IS)/AP(IS)
           AS(IS)= AS(IS)/AP(IS)
           AE(IS)= AE(IS)/AP(IS)
           AW(IS)= AW(IS)/AP(IS)
           AT(IS)= AT(IS)/AP(IS)
           AD(IS)= AD(IS)/AP(IS)
           DFRE(IS) = DFRE(IS)/AP(IS)
           AP(IS)= 1.
         ENDDO
      ENDDO

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(2956+NGL,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DFRE,IMAXSL,JMAXSL,KMAXSL,
     +   IN,JN,KN,2956+NGL)
         CALL PRINXY(3,DSURFC,DFRE,0.,ICYCLE,IMAXSL,JMAXSL,1,NGL,M)
      ENDIF

      CALL SOLAMG(DFRE,DFRE,AP,AN,AS,AE,AW,AT,AD,IMAXSL,JMAXSL,KMAXSL,
     + IN,JN,KN,MCYCLE,ITERMA,ITERHA,MGRIDA,ICONV3,0,NGL,ICYCLE)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*) ' AFTER SOLAMG:'
         CALL PRINXY(3,DSURFC,DFRE,0.,ICYCLE,IMAXSL,JMAXSL,1,NGL,M)
      ENDIF

      DO J = JA-1,JM+1
         IJ   = (JN+J-1)*JSTR + KA
*         NR = (JN+J-JA)*(IM-IA+1+2*IN) - IA + IN + IF2
         NR = (LN+J-JA)*(IM-IA+1+2*LN) - IA + LN + IF2
         ILPO = KN*ISLAB*JSLAB + (JN-JA+J)*ISLAB + IN-IA+1
         DO I = IA-1,IM+1
            NP  = I + NR                      ! Patch index with LN ghost cells
            IS    = ILPO + I
            L     = 1  + (IN+I-1)*ISTR + IJ ! CELL INDEX (STRIDE=ISTR)
            IF1 = L + NFL                     ! Face and flux index
c         write(667,*) i,j,np-if2+1,is,dfre(is)

c            DWH        = DFRE(IS)/(1.+ABS(DFRE(IS)/0.0005))
            DWH        = DFRE(IS)/(1.+ABS(DFRE(IS)/DWMV))
            IF (IGRID == 6) THEN ! jaadyttaa aaltopinnan
               DWH = 0.0           ! onko nain hyva?
            ENDIF
            DSURLE(NP) = DWH !  DFRE(IS) ! Implicit change 
            WAVEH(NP)  = WAVEH(NP) + DWH
            DFRE(IS)   = DWH 

	    IF(ABS(DWH) > ABS(DWMAX)) THEN
               DWMAX  = DWH
               IDWMAX = I
	       JDWMAX = J
            ENDIF

            WHMAX  = MAX(WHMAX,WAVEH(NP))
            WHMIN  = MIN(WHMIN,WAVEH(NP))
            SUMDWH = SUMDWH + ABS(DWH)
            IF(IFACE == 1 .OR. IFACE == 4) THEN
               FLODWH = FLODWH + DWH*A1(IF1)
            ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
               FLODWH = FLODWH + DWH*A2(IF1)
            ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
               FLODWH = FLODWH + DWH*A3(IF1)
            ENDIF

         ENDDO
      ENDDO
      IF(IPRINT == ICYCLE) THEN
         WRITE(3,*) ' FINAL VALUES:'
          CALL PRINXY(3,DSURFC,DFRE,0.,ICYCLE,IMAXSL,JMAXSL,1,NGL,M)
      ENDIF

      DEALLOCATE (AP,AE,AW,AN,AS,AT,AD,DFRE)
        
      ENDIF ! ICON(1,IP) == 13 (Free surface)

      ENDDO ! NPATCH
     
      RETURN
      END SUBROUTINE IMPFRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPSOL(NBL,M,CFM,SMAX,TEMP,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DE,RO,CH,CP,
     + VOL,VIS,EPS2,VIST,ITURB,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,PR,PRT,
     + XC,YC,ZC,TIMEL,MCYCLE,MGRIDA,ITERMA,ITERHA,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,A,B,C,D,DTEMP,FPRINTL,ICONV3) 

      USE CHARACTERS
      USE NS3CO, ONLY : IN, JN, KN
      USE CONSTANTS, ONLY : EPS

      IMPLICIT NONE

      INTEGER ::
     +     NBL,M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,ITERMA,MCYCLE,ITERHA,MGRIDA,ICONV3,IPOST,ITERMAX

      REAL :: TEMP(*),DTEMP(*),RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),
     +        A1(*),A1XA(*),A1YA(*),A1ZA(*),
     +        A2ZA(*),A2(*),A2XA(*),A2YA(*),
     +        A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +        D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),
     +        AT(*),AD(*),RESP(*),A(*),B(*),C(*),D(*)

      REAL :: XC(*),YC(*),ZC(*)

      REAL :: DT,CFM,PR,PRT,SMAX

C ... Solution of the heat conduction equation

      LOGICAL :: TIMEL,FPRINTL
        
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      KSTR = ISTR*JSTR

      DTEMP(1:NTOT) = 0.

C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG

C     ICONV3  = 0 ! No convergence in AMG
      IPOST   = 0 ! No post smoothing

C ... Define the local time step
          
      CALL TIMEHS(DTL,CH,RO,CP,D1,D2,D3,SMAX,IMAX,JMAX,KMAX,IN,JN,KN)

C ... Form the matrix for the implicit stage

      DE(1:NTOT) = -DE(1:NTOT)*VOL(1:NTOT) ! In this solution system 

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE SUBTRACTING VOLUMES'
         CALL PRINYS(3,DEC,      DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF
          
      CALL AAMATG(AP,AN,AS,AE,AW,AT,AD,D1,D2,D3,RO,CP,VOL,CH,
     +     A1,A2,A3,DTL,DE,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,DEC,      DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,CPC,      CP, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,ROC,      RO, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,CHC,      CH, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,APC,      AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(950+NBL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(950+NBL+(M-1)*1000,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DE,
     +   IMAX,JMAX,KMAX,IN,JN,KN,950+NBL+(M-1)*1000)
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DE(1:NTOT) = DE(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT)  = 0.


C ... Solve using an algebraic multigrid
        ITERMAX = 40 !; MGRIDA = 3

      CALL SOLAMG(DTEMP,DE,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
c     + IN,JN,KN,MCYCLE,ITERMA,ITERHA,MGRIDA,ICONV3,IPOST,NBL,ICYCLE)
c     + IN,JN,KN,4,1,5,MGRIDA,5,0,NBL,ICYCLE) ! Oli MGRID
     + IN,JN,KN,ITERMAX,1,5,MGRIDA,ICONV3,0,NBL,ICYCLE)

C ... Set the residual to zero

      CALL SETVAL(DE,0.,NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,TIC,     DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
         CALL PRINYS(3,DEC,      DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NBL,M)
      ENDIF

      RETURN
      END

