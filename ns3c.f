C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NS3D
      
      USE MPI

      USE CHARACTERS

      USE CONSTANTS,   ONLY : EPS

      USE INTEGERS

      USE MAIN_ARRAYS

      USE TYPE_ARRAYS

      USE NS3CO             

      USE FLIGHT , ONLY : NGRIFL,XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,PSIM,OSKU,TRMODE,DRAUGHT,SINK,TRIMA,
     &     ROTA1,ROTB1,CONEA,MVSHIP,SHIPWAVE,HMVBLADE,
     &     SHIPPROP,IGRSLAVE,FLYOBJ,SHIP,ACTDISK,FXSP,FXTSP,RTMSP,
     &     UTSP,ADV,THRUST,TORQUE,FDSP,IFA

      USE BLADE_VARIABLES , ONLY : SHAFT,THETACOL,THCYCLON,THCYCLAT,
     &     QTH,QTH_A,QBE_A,DTB,TMAXB,TDBL,QZE_A,QZE_S1,QBE_S1,QZE_D

      USE MPCCIVARS, ONLY : transfered, strongcoupling

      IMPLICIT NONE

      INTEGER :: IGR     

      CHARACTER(LEN=4) :: RUN
      CHARACTER(LEN=1) :: CLEVW

      LOGICAL :: TIME2, FIRSTC, FIRSTOUTP, THERE, MASTER

      INTEGER :: IERR, KOKO1, I, N, NGL,
     &           IG1, ICYCTI, IRUN, NTO3, NSC, IC1, IF1, IP1,
     &           IT1, IGRSLA, IGN, IPHASE, M, II, JJ, KK, ICHARP

      REAL    :: DREL, AANGLE, CLDRO, HEADT, EFFIC,
     &           XMUL1,XMUL2,QMEIN2,QMEOU2,DELTAP,GTOT,FRSG,FRSRET

c      REAL (DBLREAL) :: TDBL


C***********************************************************************
      MASTER = IPRO == 1 

      IUPDAT = 0              ! Indicates boundary updating
      DT    = DTB ; T = TDBL  ! Enough?

C***********************************************************************
C     START OF CALCULATION                                             *
C***********************************************************************
C     PROGRAM NOZZLE

C ... BIDIAGONAL LU-FACTORED METHOD FOR TWO-DIMENSIONAL VISCOUS FLOW
C ... 2.2.1988  . UPDATED 5.12.1988. THE CORRESPONDING INVISCID CODE
C ... IS LU2.     REVISED VERSION 16.3.1989. 3D VERSION NS3 31.3.1989
C ... TURBULENCE IS INCLUDED 7.2.1990. THIS IS VERSION 14.07.1991
C ... DEVELOPMENT OF THE NEW TURBULENCE MODELS STARTED 26.11.1992
C ... K-EPSILON 2.2.1993. SCALARS 1.4.1994. NEW BOUNDARY TREATMENT 20.5
C ... VERSION 16.9.1994. SO-CALLED '31.12.94' VERSION 24.2.1995
C ... 'BLOCK LOOPS' ARE MARKED FOR PARALLEL EXECUTION 24.12.1995
C ... SO-CALLED 'ORIGINAL VERSION' 3.2.1997. Free-surface model 1.6.2000
C ... Linux versions at Finflo Ltd. started from no. 7.1 (1.1.2002)
C ... Program FINFLO-7.2 (27.4.2004). Multi-phase modeling 2.2.2006
C ... FINFLO-8.0 started 1.3.2006. FINFLO-8.1 1.4.2008.
C ... FINFLO 8.2.0 started 2.2.2009

      FIRSTC    = .TRUE.   ! First cycle after start or restart
      FIRSTOUTP = .TRUE.   ! First output to MOVIE files
      NPHASE    = NPHASES

C ... INP is devided into two parts for dynamic allocations for JTRA etc.
C***********************************************************************
      IF(FRESUL) CALL INP_HULL ! Requires HULL.BIN

      CALL INP
      
C **********************************************************************
C ... Make a dynamic memory allocation for APP,JTRA and JLOC
      CALL ALLAJJ
      
C***********************************************************************
      CALL INPGRI
C***********************************************************************
      
      CALL FS_VERPAT(1,NBLOCK,1)

      CLOSE(2)
      CLOSE(41)
      CLOSE(42)
          
C ... NEW PROBLEM!
       
      IF(FRESUL .AND. IOLD > 0) THEN
         CALL NUMCH1(CLEVW,LEVEL)
         INQUIRE(FILE='WH.DAT'//CLEVW,EXIST=THERE)
         IF(THERE) THEN         ! Open the free-surface file
           OPEN(21,FILE='WH.DAT'//CLEVW,STATUS='UNKNOWN',
     +     FORM='FORMATTED')
         ELSE 
           WRITE(*,*)  ' Warning: File WH.DAT was not found.'
           WRITE(45,*) ' Warning: File WH.DAT was not found.'
           WRITE(4 ,*) ' Warning: File WH.DAT was not found.'
           IF(IFSBC == 1)THEN ! Extrapolate from pressures
             INWH = 0
           ELSE                 ! Assume zero wave height
             INWH = 1
           ENDIF
         ENDIF
      ENDIF

C *******************************************************************
C ... CREATE/READ TRAJECTORY_IGR.DAT FILE IF NGRIFL > 0  
C *******************************************************************
              
      IF (NGRIFL > 0 .AND. MASTER) THEN
         DO IGR = 1,NGRIFL     
c           IF(TRMODE(IGR) /= 38) THEN
            CALL TRAJECTORY_FORMAT(TDBL,DTB,XCG(IGR),YCG(IGR),ZCG(IGR),
     +           PSIR(IGR),THETAR(IGR),PHIR(IGR),PSIM(IGR),
     +           OSKU,IGR,TIMEL,TRMODE(IGR),DRAUGHT(IGR),.TRUE.)
c           ENDIF
            IGRSLA = IGRSLAVE(IGR)
            IF (SHIP .AND. ACTDISK .AND. IGRSLA /= 0) THEN 
               FXSP(IGR)  = OSKU(IGRSLA)%CX  ! Self propulsion updates
               FXTSP(IGR) = OSKU(IGR)%FXT
            ENDIF
         ENDDO
      ENDIF ! IF (NGRIFL > 0 ...) 
      
C ... Trajectory data is send to slave prosessors
      IF(PARALLEL .AND. ACTDISK) THEN 
         IF (SHIP) THEN
            CALL MPI_BCAST(FXSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FXTSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FDSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THRUST,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(TORQUE,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(UTSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ADV,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)  
            CALL MPI_BCAST(SHAFT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR) 
            CALL MPI_BCAST(OSKU%DAMPN,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%DAMPT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR) ! Self propulsion updates
         ELSE ! Patria's ACT disk
            CALL MPI_BCAST(FXTSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(RTMSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FDSP,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ROTA1,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ROTB1,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CONEA,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THETACOL,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THCYCLON,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THCYCLAT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THRUST,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%FXT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%FYT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%FZT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%MXT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%MYT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(OSKU%MZT,NGRIFL+10,MPI_REAL8,0,
     +           MPI_COMM_WORLD,IERR)
         ENDIF
      ENDIF
      DO IGR = 1,NGRIFL
         IF (PARALLEL .AND. IFA(IGR) == 12 .OR. ! Patria's ACT disk
     +        PARALLEL .AND. IFA(IGR) == 13) THEN
            QTH(IGR)   = -1.E10
            QTH_A(IGR) = -1.E10
            QBE_A(IGR) = -1.E10
         ELSEIF (PARALLEL .AND. IFA(IGR) == 8 .OR. ! TQIHM
     +           PARALLEL .AND. IFA(IGR) == 10) THEN
            QTH(IGR)    = 0.0
            QZE_A(IGR)  = 0.0
            QTH_A(IGR)  = 0.0
            QBE_A(IGR)  = 0.0
            QZE_S1(IGR) = 0.0
            QBE_S1(IGR) = 0.0
            QZE_D(IGR)  = 0.0
         ENDIF
         
      ENDDO
      IF (PARALLEL .AND. TIMEL .AND. NGRIFL > 0) THEN
         CALL MPI_BCAST(OSKU,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETAR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIM,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
      ENDIF
      IF(PARALLEL .AND. SHIP .AND. .NOT. TIMEL) THEN
         CALL MPI_BCAST(OSKU,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETAR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIM,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
      ENDIF         
C ... Unexpected flow field initialization needed when IN, JN and KN are 
C ... greater than the number of active ghost cell slabs. Antakee Audi.

      RO(1:MAXB)   = 1.0
      P(1:MAXB)    = 1.0
      TEMP(1:MAXB) = 1.0
     
C *******************************************************************
      IF (IOLD <= 0) CALL INIT
C *******************************************************************
 
C ... Read the solution part of the RSTART-file and NEWTIM
C ... Some initialization inside as IOLD <= 0 (sigh)
 
      CALL READRS(TIME2,DREL)

C ... Leikataan READRS-ohjelman loppuun. Definitions needed there ???

      DO N = 1,NBLOCK

         M   = 1
         IC1 = IC(M,N)
         IG1 = IG(M,N)
         IP1 = 1
         IT1 = 1
         IF(MULPHL) IP1 = IG1
         IF(TRANSL) IT1 = IG1
         NGL = NPROCE(1+N,IPRO)     ! Global block number

         FRSDEN  = BLKS(NGL)%FRSDEN
         FRSPRE  = BLKS(NGL)%FRSPRE
         FRSVEL  = BLKS(NGL)%FRSVEL
         FRSVIS  = BLKS(NGL)%FRSVIS
         FRSTEM  = BLKS(NGL)%FRSTEM
         FRSSIE  = BLKS(NGL)%FRSSIE
         CHLREF  = BLKS(NGL)%CHLREF
         RKLIM   = BLKS(NGL)%RKLIM
         EPSLIM  = BLKS(NGL)%EPSLIM
         FRSMUT  = BLKS(NGL)%FRSMUT
         SOLTEM  = BLKS(NGL)%SOLTEM
         FRSG    = BLKS(NGL)%FRSG
         FRSRET  = BLKS(NGL)%FRSRET

         CALL GHOSTEXT(FRSVEL,ALPHA,BETA, ! Update values in EXT ghost cells
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,
     +        XC(IG1),YC(IG1),ZC(IG1),
     +        NPATCH(N),ICON(IC1),FRSDEN,FRSPRE,FRSSIE,FRSTEM,FRSVIS,
     +        RKLIM,EPSLIM,FRSMUT,SOLTEM,ITURB,NGL,
     +        RO(IG1),P(IG1),PDIFF(IG1),RM(IG1),RN(IG1),RW(IG1),
     +        U(IG1),V(IG1),W(IG1),E(IG1),TEMP(IG1),VIS(IG1),RK(IG1),
     +        REPS(IG1),VIST(IG1),EPS2(IG1),RNUT(IG1),BLKS,PRO(IP1),
     +        VAR(IP1),MULPHL,INITC(N),FRSG,FRSRET,TRANSL,TRM(IT1),
     +        CENAX(N),CENAY(N),CENAZ(N))

      ENDDO

C ... READRS this far

      IF(FRESUL) INQUIRE(FILE='WH.DAT'//CLEVW,EXIST=THERE)
      IF(FRESUL .AND. THERE .AND. IOLD >= 1) THEN

*         CALL READWH(IPRO,PARALLEL,.FALSE.)
         CALL READWH(.FALSE.)

         IF(IFSBC == 1) THEN ! Use the WH array
           DO N = 1,NBLOCK
           IF1  = JF(1,N)
           IC1  = IC(1,N)
           CALL FS_SETWAVEHBACK(WAVEH,WH(IF1),IHF(1,1),ICON(IC1),
     +     NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ICYCLE)
           INWH = 1
           ENDDO
         ENDIF ! IFSBC == 1
      ENDIF    ! FRESUL
         
      IF(MASTER) THEN
         CLOSE(20)
         CLOSE(21)
      ENDIF
      
      IF(COORL .AND. IOLD /= 0) CALL GRIORI     ! Initialize XGRI from XORI
      
      IF(IOLD >= 1) THEN

ccc         IF(COORL) CALL GRIORI ! Initialize XGRI from XORI

         DO I = 1,MAXB

            DELTAP  = FRSDEN*(GX*(XC(I) - GROUND) + GY*(YC(I) - GROUND)
     +              +         GZ*(ZC(I) - GROUND))
            P(I)    = PDIFF(I) + FRSPRE + DELTAP
            IF(MULPHL) THEN
               PRO(I)%TSAT = TEMP(I)
            ENDIF

         ENDDO

      ENDIF

c      IF(IOLD <= 0) THEN ! Read viscosity. For tests
c       OPEN(927,FILE='EPS.FMT',FORM='FORMATTED')
c        do I=1,NTOT(1,1)
c       READ(927,*)I,EPS2(I),VIST(I)
c       ENDDO
c       ENDIF
         
C***********************************************************************
      IF(NCHIM > 0 .AND. IOLD < 0) THEN 
         CALL CHIMER(1)
         CALL CHIMUP(1)
      ENDIF
C***********************************************************************

C ... Initialize (free surface), Reynolds stresses and boundary blocks

C***********************************************************************
      CALL INPOST(DREL,TIME2)
C***********************************************************************
         
C ... INITIALIZATION OF REYNOLDS STRESSES

C***********************************************************************
      IF (ITURB >= 21 .AND. ITURB < 23) CALL BOUSS
C***********************************************************************

C***********************************************************************

      CALL SUBIN(0,FIRSTC)
      CALL BOUNDI(1)

C***********************************************************************

      CALL CPUTIM

C **********************************************************************
C     TIME INTEGRATION IS STARTED
C     TDBL = T

      CALL MPCCI_INITCOUPLING
      CALL MODALFSI_INITCOUPLING

900   CONTINUE
C **********************************************************************
 
      IF(TIMEL) THEN

C         T     = T + DT
         TDBL   = TDBL + DTB
         T      = TDBL
         ITIMES = ITIMES + 1
         ICYCLE = 0
C        IPRINT = 1  ! IF ACCORDING TO CYCLES INSIDE THE TIME-STEP

C ... The following line is for droplet evaporation model. The inlet 
C ... patches should be in a same processor (output EVAPMF in kg/s)
C ... Subroutine DROPLET should make one time step or several steps
C ... inside DT using constant pressure PJET. This only for subsonic
C ... cases and evaporation model must take care of that. Evaporation
C ... module should be put to the last position in charac.f

c        CALL DROPLET(T,DT,PJET)
c        CALL BASE_BLEED

        DO N = 1,NBLOCK
            CALL STATIM(1,N,2,.TRUE.) ! Should be removed?
        ENDDO

         IF(ITIMES == IPRINT) THEN
            WRITE(3,*)
            WRITE(3,*) ' TIME = ',REAL(T,4),
     +                 ' DT = ',REAL(DT,4),
     +                 ' STEPS = ',REAL(ITIMES,4)
            WRITE(3,*) ' DROMAX FROM THE PREVIOUS CYCLE = ',
     +                   REAL(DROMAX,4)
            DO N = 1,NBLOCK
               NGL     = NPROCE(1+N,IPRO) !Global block number
               IG1 = IG(1,N)
               CALL PRINYS(3,ROC,RO(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +              IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,'  ROLE2 ',ROLE2(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,'  ROLE3 ',ROLE3(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               IF(MULPHL) THEN
               IGN  = IG1 + NTOT(1,N) - 1
               DO IPHASE = 1,BLKS(NGL)%NPHASE
               F1RM(IG1:IGN) = VAR(IG1:IGN)%AROLE2(IPHASE)
               F1RN(IG1:IGN) = VAR(IG1:IGN)%AROLE3(IPHASE)
               CALL PRINYS(3,' AROLE2 ',F1RM(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' AROLE3 ',F1RN(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               ENDDO ! IPHASE
               ENDIF
            ENDDO
         ENDIF


         CALL NEWTIM

         CALL MPCCI_DOTRANSFER(FIRSTC)
         CALL MODALFSI_DOTRANSFER
       
C ... GENERATE A NEW GRID AT THE BEGINNING OF THE ITERATION LOOP
            
C***********************************************************************
         IF (FRESUL) CALL SURGRI(1)
         IF (COORL)  CALL NEWROT
C***********************************************************************
     
      ENDIF  ! TIMEL
      
C **********************************************************************
C     ITERATION IS STARTED
C **********************************************************************

1000  CONTINUE

      IF(.NOT.TIMEL .OR. TIMEL .AND. StrongCoupling /= 0) THEN

         CALL MPCCI_DOTRANSFER(FIRSTC)
         CALL MODALFSI_DOTRANSFER

C ... GENERATE A NEW GRID AT THE BEGINNING OF THE ITERATION LOOP
      
C***********************************************************************
         IF (FRESUL) CALL SURGRI(1)
         IF (COORL)  CALL NEWGRI
C***********************************************************************

      ENDIF  ! .NOT.TIMEL
      
      ICYCLE  = ICYCLE + 1
      ICYTOT  = ICYTOT + 1

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT) THEN
          WRITE(3,*)
          WRITE(3,*) ' CYCLE = ',ICYCLE,' DT = ',REAL(DT,4)
          WRITE(3,*) ' DROMAX FROM THE PREVIOUS CYCLE = ',REAL(DROMAX,4)
      ENDIF

      IF(TIMEL .AND. ICYCLE == IPRINT) THEN ! OUTPUT WITHIN TIME-STEP
C     IF(TIMEL .AND. ITIMES == IPRINT .AND. ICYCLE == 1) THEN 
          WRITE(3,*)
          WRITE(3,*) ' CYCLE = ',ICYCLE,' DT = ',REAL(DT,4)
          WRITE(3,*) ' DROMAX FROM THE PREVIOUS CYCLE = ',REAL(DROMAX,4)
      ENDIF

C ... Initialize global forces and moments

      XMASSB(1) = 0.
      XMASSB(2) = 0.
      XMASSB(3) = 0.
      XMASSB(4) = 0.
      QMFIN     = 0.
      QMEIN     = 0.
      QVFIN     = 0.
      QMFOUT    = 0.
      QVFOUT    = 0.
      QMEOUT    = 0.
      QMXIN     = 0.
      QMYIN     = 0.
      QMZIN     = 0.
      QMXOUT    = 0.
      QMYOUT    = 0.
      QMZOUT    = 0.
      CMXIN     = 0.
      CMYIN     = 0.
      CMZIN     = 0.
      CMXOUT    = 0.
      CMYOUT    = 0.
      CMZOUT    = 0.
      IF(MULPHL) THEN
         QGFIN   = 0.
         QGFOUT  = 0.
         QGEIN   = 0.
         QGEOUT  = 0.
         QLFIN   = 0.
         QLFOUT  = 0.
         QLEIN   = 0.
         QLEOUT  = 0.
      ENDIF

C ... CHIMER subroutine is visited every 1.cycle in time accurate calculation 

      IF (TIMEL) THEN
         IF (ICYCLE > 1) THEN
            FIRSTC =.FALSE.
         ELSE
            FIRSTC =.TRUE.
         ENDIF
      ENDIF
      
C ######################################################################
      CALL CYCLES
C ######################################################################
   
C***********************************************************************
C ... ITERATION CYCLE IS COMPLETED. UPDATE BOUNDARIES AND SUBIN ARRAYS
C **********************************************************************

C ... START TIME TICK FOR COMMUNICATION IN MPI
      IF(PARALLEL) WTBCST = MPI_WTIME()

C***********************************************************************
      CALL BOUNDI(1)
C***********************************************************************
*         CALL MPCCI_DOTRANSFER(FIRSTC)
C***********************************************************************
     
      IF(NCHIM > 0) THEN ! Make Chimera interpolation

         IF(.NOT.TIMEL) THEN

            IF(IOLD >= 0  .AND. FIRSTC 
     +      .OR. MPCCIL   .AND. transfered > 0
     +      .OR. HMVBLADE 
     +      .OR. MVSHIP   .AND. MOD(ICYCLE,10) == 0 
     +      .OR. SHIPWAVE .AND. MOD(ICYCLE,10) == 0)
     +           CALL CHIMER(1)
         ELSEIF(TIMEL) THEN

            IF(FIRSTC .OR. MPCCIL .AND. transfered > 0) CALL CHIMER(1)   

         ENDIF

         
         CALL CHIMUP(1)

         
         IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)
            WRITE(3,*) ' From CHIMUP'
            DO N = 1,NBLOCK
               NGL     = NPROCE(1+N,IPRO) !Global block number
               IG1 = IG(1,N)
               CALL PRINYS(3,ROC,RO(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +              IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
            ENDDO
         ENDIF

      ENDIF ! NCHIM > 0

C***********************************************************************

C ... END TIME TICK FOR COMMUNICATION IN MPI
      IF(PARALLEL) WTBCEN = MPI_WTIME()

      ICYCTI = ICYTOT

C***********************************************************************
      CALL SUBIN(ICYCTI,FIRSTC)
C***********************************************************************

      IUPDAT = 1
      INWH = 1        !H.Piippo

C ... END OF ITERATION CYCLE

      IF(TIMEL) THEN 
                     
C ... IF YOU WANT TO HAVE OUTPUT, TO MAKE MOVIE OR TO FOLLOW
C ... CONVERGENCE DURING THE TIME-STEP, COMMENT THE MARKED LINES

         DTOLD  = DT                   ! Hope this works    
         ICYCTI = ITIMES               ! COMMENT THIS
         IF(ABS(DROMAX) < DROLIM .AND. ICYCLE >= 2) GO TO 1900
         IF(ICYCLE >= ICMAX) GO TO 1900
         GO TO 1000                    ! COMMENT THIS 
      ENDIF ! TIMEL
1900  CONTINUE

C **********************************************************************
C ... PROCESS POSSIBLE OUTPUT/(ITERATION CYCLE OR TIME-STEP)
C **********************************************************************

C **********************************************************************
      CALL CYCOUT(ICYCTI,FIRSTC,FIRSTOUTP)
C **********************************************************************

      IF(STARTL .AND. MOD(ICYCLE,KRRINT) == 0) THEN ! Write Rstart
         CALL CALLRS(IPRO,RUN,1)
      ENDIF ! STARTL

C ... Terminate iteration

      IF(ABS(DROMAX) < DROLIM .AND. ICYCLE >= 3) GO TO 2000
      IF(ICYCLE >= ICMAX) GO TO 2000

C ... START INTERRUPT ROUTINE ===============================
C ... INTERUPT EXCECUTION?
C ... IN PARALLEL COMPUTING ONLY MASTER PROCESS READS RUN-FILE

      IF(CONVL) THEN

         CALL INTERRUPT(IPRO,IRUN)
         IF(IRUN == 2) RETURN

C ... Write Rstart during the interrupt 

         IF(IRUN == 3)THEN
            CALL CALLRS(IPRO,RUN,1)
         ENDIF ! IRUN == 3

         IF(IRUN < 1 .OR. IRUN > 3) GO TO 2000

      ENDIF ! CONVL

C ... END INTERRUPT ROUTINE ===============================

      FIRSTC = .FALSE.

      GO TO 1000  

C ********************************************************************
C     OUTPUT/ITERATION CYCLE OR TIME-STEP IS FINISHED                *
C ********************************************************************

2000  CONTINUE

c      IF(TIMEL) THEN
c         CALL MPCCI_DOTRANSFER(FIRSTC)
c         CALL MODALFSI_DOTRANSFER
c      ENDIF
      
C **********************************************************************
C ..  CONVERGENCE MONITOR ACCORDING TO TIME-STEPS
C **********************************************************************

      IF(TIMEL) THEN
         IF(MASTER) THEN
            WRITE(6,9993) T,ICYCTI
            WRITE(7,9993) T,ICYCTI
            IF(PARALLEL) WRITE(11,9993) T,ICYCTI
         ENDIF
         CALL CONMON(ITIMES,DROMAX,CL,CD,CMZ,IXERR,INERR,IMACH,
     +        IMAX(1,INERR),JMAX(1,INERR),IPRO,NPROCE,MBPRO,
     +        NPRO,CONVL)  
         IF(MASTER) THEN
            IF(COORL) THEN
               AANGLE = -1000000.
               DO N = 1,NBLOCG
                  IF(OMEGA(N) /= 0.) AANGLE = MAX(AANGLE,ROTANG(N))
               ENDDO
               WRITE(6,9995) AANGLE*180./3.1415927
               WRITE(7,9995) AANGLE*180./3.1415927
               IF(PARALLEL) WRITE(11,9995) AANGLE*180./3.1415927
             ELSE
               WRITE(6,9994)
               WRITE(7,9994)
               IF(PARALLEL) WRITE(11,9994)
            ENDIF
         ENDIF
 9993 FORMAT(' --TIME=',E14.5,' seconds-------','ICYCTI=',I7,
     +       '--------------------------')
 9994 FORMAT(' --------------------------------------------------',
     +       '--------------------------')
 9995 FORMAT(' ---------------------------------ANGLE =',F14.5,
     +       ' degrees--------------')
         CLDRO     = ABS(DROMAX)+1.E-15
************************************************************************

C **********************************************************************
C ... UPDATE TBDIAG
C **********************************************************************

      NTO3 = 0
      DO 8195 N = 1,NBLOCK
         NTO3 = NTO3 + IMAX(1,N)*JMAX(1,N)*KMAX(1,N)
 8195 CONTINUE

      CLDRO     = ABS(DROMAX)+1.E-15
      QMEIN2    = QMEIN /(QMFIN +1.E-20)
      QMEOU2    = QMEOUT/(QMFOUT+1.E-20)
      EFFIC     = 100.*ABS(QMFIN*(QMEOU2-QMEIN2)/(TOM+1.E-20)) ! ADD HOC
      HEADT     = (QMEOU2 - QMEIN2)/G0
      IBOUT     = 1

      CALL IBDWRI(59,IBOUT,ITIMES,T,CLDRO,RESI,'TBDIAG',PRN,CD,CL,
     +     CS,CMX,CMY,CMZ,IMACH,NPROCE,MBPRO,NPRO,NTO3,IPRO,
     +     NRESI,VOLTOT,CONVL,QMFIN,QMFOUT,QMEIN,QMEOUT,HEADT,EFFIC,
     +     TOM,ZZZ,JET,QVFIN,QVFOUT,DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH,
     +     SINK,TRIMA,QGFIN,QGFOUT,QGEIN,QGEOUT)

C ... Update total forces in time accurate calculation
      
      IF(GROUP) THEN
      CALL AUXFOR(156,CXB,CYB,CZB,CMXB,CMYB,CMZB,DXB,DYB,DZB,QTB,QWB,
     +     QHB,TOMEGA,ALPHA,BETA,ICON,BOUNDF,NBCS,IPRO,3,T,
     +     ITIMES,APATCH,XMOM,YMOM,ZMOM,REFPRE,AREF,CHLREF)
      ENDIF

C***********************************************************************
C ... Update "TRAJECTORY_IGR.DAT" file in time accurate calculation
C***********************************************************************
      IF (NGRIFL > 0 .AND. TIMEL) THEN
      DO IGR = 1,NGRIFL ! particle_forces loop
         IGRSLA = IGRSLAVE(IGR)
      CALL PARTICLE_FORCES(OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     +     OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,ICON,
     +     XCG(IGR),YCG(IGR),ZCG(IGR),2,10,10,0,IGR,IGRID(:,3)) ! MOV surfaces
      IF(INOUTL .AND. FLYOBJ) THEN ! thrust for the flying objects
      CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +    OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +    XCG(IGR),YCG(IGR),ZCG(IGR),2,3,5,0,IGR,IGRID(:,3)) ! INL and OUT surf 
      ENDIF
      IF(SHIPPROP .AND. IGRSLA /= 0) THEN ! SHIP PROPULSION
      CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +    OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +    XCG(IGR),YCG(IGR),ZCG(IGR),2,9,9,0,IGR,IGRID(:,3))  ! ROT surfaces
      CALL PARTICLE_FORCES(OSKU(IGRSLA)%FXT,OSKU(IGRSLA)%FYT,
     +     OSKU(IGRSLA)%FZT,OSKU(IGRSLA)%MXT,OSKU(IGRSLA)%MYT,
     +     OSKU(IGRSLA)%MZT,ICON,XCG(IGRSLA),YCG(IGRSLA),ZCG(IGRSLA),
     +     2,9,9,0,IGR,IGRID(:,3)) ! ROT surfaces, slave IGRID
      ENDIF
      ENDDO

      DO IGR = 1,NGRIFL ! trajectory_write loop
      IF (MASTER) THEN
         CALL ROTOR_FORCES(OSKU(IGR)%RMASS,XCG(IGR),YCG(IGR),ZCG(IGR),
     +        OSKU(IGR)%ROTX,OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,
     +        OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,! helikopterin roottorin 
     +        OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,! keskihaku ja gravitaatio
     +        OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ, ! tama on 6dof systeemi
     +        OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE, ! ja tata ei kayteta
     +        OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE,
     +        OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +        OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     +        OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     +        OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH,
     +        PSIR(IGR),THETAR(IGR),PHIR(IGR),OSKU(IGR)%RCG,
     +        OSKU(IGR)%IX,OSKU(IGR)%IY,OSKU(IGR)%IZ,
     +        OSKU(IGR)%IXY,OSKU(IGR)%IXZ,OSKU(IGR)%IYZ,
     +        ROTORW,NBLOCK,TRMODE(IGR),IGR)
         CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     +        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     +        PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,TRMODE(IGR),
     +        ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
         CALL TRAJECTORY_CONVERTER(IGR,TRMODE(IGR),TDBL,DTB)
         CALL TRAJECTORY_BACKUP(IGR)
      ENDIF
      ENDDO

      ENDIF ! IF (NGRIFL > 0 .AND. TIMEL) THEN

C***********************************************************************
C ... Calculate averaged values
C***********************************************************************

      XMUL1 = TAVER/(TAVER+DTB)
      XMUL2 = DTB  /(TAVER+DTB)

      DO I = 1,MAXB ! korjattu
         ROAV1(I) = XMUL1*ROAV1(I) + XMUL2*RO(I)
         RMAV1(I) = XMUL1*RMAV1(I) + XMUL2*RM(I)
         RNAV1(I) = XMUL1*RNAV1(I) + XMUL2*RN(I)
         RWAV1(I) = XMUL1*RWAV1(I) + XMUL2*RW(I)
         EAV1(I)  = XMUL1*EAV1(I)  + XMUL2*E(I)
         PAV1(I)  = XMUL1*PAV1(I)  + XMUL2*P(I)
         TAV1(I)  = XMUL1*TAV1(I)  + XMUL2*TEMP(I)
      ENDDO
      DO I = 1,MAXB
         ROAV2(I) = XMUL1*ROAV2(I) + XMUL2*RO(I)**2
         RMAV2(I) = XMUL1*RMAV2(I) + XMUL2*RM(I)**2/(RO(I)+EPS)
         RNAV2(I) = XMUL1*RNAV2(I) + XMUL2*RN(I)**2/(RO(I)+EPS)
         RWAV2(I) = XMUL1*RWAV2(I) + XMUL2*RW(I)**2/(RO(I)+EPS)
         EAV2 (I) = XMUL1*EAV2 (I) + XMUL2* E(I)**2/(RO(I)+EPS)
         RUVAV(I) = XMUL1*RUVAV(I) + XMUL2*RM(I)*RN(I)/(RO(I)+EPS)
         RUWAV(I) = XMUL1*RUWAV(I) + XMUL2*RM(I)*RW(I)/(RO(I)+EPS)
         RVWAV(I) = XMUL1*RVWAV(I) + XMUL2*RN(I)*RW(I)/(RO(I)+EPS)

         RMAV3(I) = RMAV2(I) - RMAV1(I)**2/(ROAV1(I)+EPS)
         RNAV3(I) = RNAV2(I) - RNAV1(I)**2/(ROAV1(I)+EPS)
         RWAV3(I) = RWAV2(I) - RWAV1(I)**2/(ROAV1(I)+EPS)
      ENDDO

      IF(ITURB >= 3 .AND.ITURB /= 8) THEN
      DO I = 1,MAXB     ! Ottiko t
         RKAV1(I)  = XMUL1*RKAV1(I) + XMUL2*RK(I)
         RKAV2(I)  = XMUL1*RKAV2(I) + XMUL2*RK(I)**2
         EPSAV1(I) = XMUL1*EPSAV1(I) + XMUL2*REPS(I)
         EPSAV2(I) = XMUL1*EPSAV2(I) + XMUL2*REPS(I)**2
      ENDDO
      ENDIF

      DO NSC = 1,NSCAL
      DO I   = 1,MAXSB  ! Multiple loop
         FIAV1(I,NSC) = XMUL1*FIAV1(I,NSC) + XMUL2*FI(I,NSC)
         FIAV2(I,NSC) = XMUL1*FIAV2(I,NSC) + XMUL2*FI(I,NSC)**2
      ENDDO
      ENDDO  ! Antakee mersu

      TAVER = TAVER + DT
      NAVER = NAVER + 1

C ... START INTERRUPT ROUTINE ===============================
C ... INTERUPT EXCECUTION?
C ... IN PARALLEL COMPUTING ONLY MASTER PROCESS READS RUN-FILE

      IF(CONVL) THEN

         CALL INTERRUPT(IPRO,IRUN)
         IF(IRUN == 2) RETURN

C ... Write time-accurate restart in interrupt

         IF(IRUN == 3)THEN
            CALL CALLRS(IPRO,RUN,2)
         ENDIF ! IRUN == 3

        IF(IRUN < 1 .OR. IRUN > 3) GO TO 2013

      ENDIF ! CONVL

C ... END INTERRUPT ROUTINE ===============================

      IF((TDBL+DTB/2.) < TMAX) GO TO 900

      ENDIF ! TIMEL

C ********************************************************************
C     TIME INTEGRATION IS FINISHED                                   *
C ********************************************************************

 2013 CONTINUE

      CALL MPCCI_EXITCOUPLING
      CALL MODALFSI_EXITCOUPLING

C ... CPU TIME

      CALL CPU_TIME(RTIME)

C ... WRITE ON THE RESTART FILE AND CLOSE FILES 

      CALL CALLRS(IPRO,RUN,3)

      CLOSE(7)                           ! MONITOR
      CLOSE(44)                          ! MCHECK
      CLOSE(46)                          ! GRIPUT
      CLOSE(90)                          ! MASSB
      IF(MASTER) CLOSE(59)               ! TBDIAG
      IF(MASTER) CLOSE(9)                ! IBDIAG
      IF(MASTER .AND. GROUP) CLOSE(155)  ! IBDIAG_FGS
      IF(PRN == '001') CLOSE(11)         ! MONITOR
      WRITE(45,*) 'NS3D: MONITOR, MCHECK... succesfully closed.'
      
C **********************************************************************
      CALL OUTP(1,FIRSTOUTP)  ! Output at the end of simulation
C **********************************************************************
       
      CLOSE(3)

      DEALLOCATE (ZZZ,STAT=IERR)
         IF(IERR > 0) WRITE(45,*) 'NS3D: ZZZ succesfully deallocated.'
      DEALLOCATE(JLOC,JTRA,APP,STAT=IERR)
         IF(IERR > 0) WRITE(45,*) 'JLOC... succesfully deallocated.'
      IF(NCHIM > 0) THEN 
         DEALLOCATE(SKINX,SKINY,SKINZ,SKINXN,SKINYN,SKINZN,
     &   IBSKIN,ITSKIN,STAT = IERR)
         IF(IERR > 0) WRITE(45,*) 'SKINX... succesfully deallocated.'
      ENDIF
      DEALLOCATE (FTR2,STAT=IERR)
         IF(IERR > 0) WRITE(45,*) 'FTR2... succesfully deallocated.'

      RETURN
      END SUBROUTINE NS3D
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CALLRS(IPRO,RUN,ISELEC)

      USE INTEGERS, ONLY : MGM,MAXSS,MAXSB,IBF

      USE NS3CO

      USE MAIN_ARRAYS

      USE BLADE_VARIABLES, ONLY : TDBL,DTB

      IMPLICIT NONE

      INTEGER :: IPRO,ISELEC,N,NGL,IGR
      REAL    :: RC
      LOGICAL :: MASTER
      CHARACTER(LEN=6) :: RSTNUM
      CHARACTER(LEN=4) :: RUN
      CHARACTER(LEN=1) :: CLEVW

************************************************************************
      CHARACTER(LEN=80) :: COMMAND
************************************************************************

      MASTER = IPRO == 1


      IF(STARTL) THEN                    ! Start to write the restart

      IF(MASTER .AND. ISELEC == 1) THEN  ! Open the appropriate files

         CALL NUMCH6(RSTNUM,ICYCLE)
         CALL NUMCH4(RUN,LEVEL)

************************************************************************
C ... Make backups of the force group files. Note that under Windows the 
C ... copy command is copy.

         IF(GROUP) THEN
            COMMAND = 'cp IBDIAG_FGS IBDIAG_FGS_'//RSTNUM//'.'//RUN(4:4)
            COMMAND = TRIM(COMMAND)
            CALL EXECUTE_COMMAND_LINE(COMMAND)
         ENDIF
         IF(GROUP .AND. TIMEL) THEN
            COMMAND = 'cp TBDIAG_FGS TBDIAG_FGS_'//RSTNUM//'.'//RUN(4:4)
            COMMAND = TRIM(COMMAND)
            CALL EXECUTE_COMMAND_LINE(COMMAND)
         ENDIF

************************************************************************

         OPEN(19,FILE='Ibdiag'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')
         IF(TIMEL)  OPEN (69,FILE='Tbdiag'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(20,FILE='Rstart'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')

         WRITE(20)

     &         TDBL,    DTB,     THETA,   H0,      DROMAX,   
     &         CD,      CL,      CS,      
     &         CX,      CY,      CZ,      CMX,     CMY,     CMZ,   
     &         RMUFRS,  FRSMUT,  TOM,    
     &         FRSVIS,  FRSSIE,  T0,      ANGLE,   
     &         VOLTOT,  QMFIN,   QMEIN,   QMFOUT,  QMEOUT,  RKMAX,  
     &         TAVER,   NBLOCG,  ROTANG(1:NBLOCG),        

     &         ICYCLE,  IPRINT,  IXERR, 
     &         INERR,   IMACH,   NBLOCK,  ICYOLD, 
     &         IUPDAT,  NSCAL,   JSCAL,   IBOUT,   ITIMES,  ICYTOT,
     &         NAVER

         IF(FRESUL) OPEN (21,FILE='Wh.dat'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='FORMATTED')
         IF(COORL)  OPEN(28,FILE='Xyd.bin'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')

      ELSEIF(MASTER .AND. ISELEC == 2) THEN

         CALL NUMCH6(RSTNUM,ITIMES)
         CALL NUMCH4(RUN,LEVEL)
         OPEN(19,FILE='IbdiagTIME'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')
         IF(TIMEL)  OPEN(69,FILE='TbdiagTIME'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(20,FILE='RstartTIME'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')

         WRITE(20)

     &         TDBL,    DTB,     THETA,   H0,      DROMAX,   
     &         CD,      CL,      CS,      
     &         CX,      CY,      CZ,      CMX,     CMY,     CMZ,   
     &         RMUFRS,  FRSMUT,  TOM,    
     &         FRSVIS,  FRSSIE,  T0,      ANGLE,   
     &         VOLTOT,  QMFIN,   QMEIN,   QMFOUT,  QMEOUT,  RKMAX,  
     &         TAVER,   NBLOCG,  ROTANG(1:NBLOCG),        

     &         ICYCLE,  IPRINT,  IXERR, 
     &         INERR,   IMACH,   NBLOCK,  ICYOLD, 
     &         IUPDAT,  NSCAL,   JSCAL,   IBOUT,   ITIMES,  ICYTOT,
     &         NAVER

         IF(FRESUL) OPEN(21,FILE='Wh.datTIME'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='FORMATTED')
         IF(COORL) OPEN(28,FILE='Xyd.binTIME'//RSTNUM//'.'//RUN(4:4),
     +        STATUS='UNKNOWN',FORM='UNFORMATTED')

      ELSEIF(MASTER .AND. ISELEC == 3) THEN

         OPEN(19,FILE='IBDIAG.BAK',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(69,FILE='TBDIAG.BAK',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(20,FILE='RSTART',STATUS='UNKNOWN',FORM='UNFORMATTED')

         WRITE(20)

     &         TDBL,    DTB,     THETA,   H0,      DROMAX,   
     &         CD,      CL,      CS,      
     &         CX,      CY,      CZ,      CMX,     CMY,     CMZ,   
     &         RMUFRS,  FRSMUT,  TOM,    
     &         FRSVIS,  FRSSIE,  T0,      ANGLE,   
     &         VOLTOT,  QMFIN,   QMEIN,   QMFOUT,  QMEOUT,  RKMAX,  
     &         TAVER,   NBLOCG,  ROTANG(1:NBLOCG),        

     &         ICYCLE,  IPRINT,  IXERR, 
     &         INERR,   IMACH,   NBLOCK,  ICYOLD, 
     &         IUPDAT,  NSCAL,   JSCAL,   IBOUT,   ITIMES,  ICYTOT,
     &         NAVER

         IF(FRESUL) THEN 
         CALL NUMCH1(CLEVW,LEVEL)
         OPEN(21,FILE='WH.DAT'//CLEVW,STATUS='UNKNOWN',FORM='FORMATTED')
         ENDIF
         OPEN(28,FILE='XYD.BIN',FORM='UNFORMATTED')

      ENDIF

      CALL WRITRS
       
      IF(FRESUL) CALL WRITWH

      IF(COORL .OR. ISELEC == 3)  CALL WRIXYD

      ELSEIF (MASTER) THEN
         WRITE(*,*)
         WRITE(*,*) '   WARNING: RSTART-FILE WAS NOT WRITTEN'
         WRITE(13,*)
         WRITE(13,*)'   WARNING: RSTART-FILE WAS NOT WRITTEN'

      ENDIF ! STARTL


ccc      DO N = 1,NBLOCK

ccc      NGL  = NPROCE(N+1,IPRO)

ccc      IF(IGRID(NGL,1) == 38 .OR. IGRID(NGL,1) == 39) THEN ! Helicopter blade
ccc      IGR  = IGRID(NGL,3)
ccc      CALL BLADE_OUTPUT(ICYCLE,T,IGR)
ccc      ENDIF
ccc      ENDDO

C ... Write a back-up file for IBDIAG and TBDIAG

      IF(MASTER) THEN
         IF (.NOT. TIMEL) THEN 
            CALL WTCBAK(9,19,ICYCLE)
         ELSEIF (TIMEL) THEN
            CALL WTCBAK(9,19,ICYTOT)
            CALL WTCBAK(59,69,ICYCLE)
         ENDIF     ! IF (.NOT. TIMEL) ...
      ENDIF

      IF(MASTER) THEN
         CLOSE(19)                               ! IBDIAG.BAK
         CLOSE(20)                               ! RSTART
         IF(COORL .OR.ISELEC == 3)  CLOSE(28)    ! XYD.BIN
         IF(FRESUL) CLOSE(21)                    ! WH.DAT
         IF(TIMEL)  CLOSE(69)                    ! TBDIAG.BAK
      ENDIF

      RETURN
      END SUBROUTINE CALLRS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INTERRUPT(IPRO,IRUN)

      USE NS3CO, ONLY : ICMAX, PARALLEL

      USE MPI

      IMPLICIT NONE

      INTEGER :: IPRO,IRUN,IERR
      LOGICAL :: MASTER

      CHARACTER(LEN=4) :: RUN

      MASTER = IPRO == 1 

      IF(MASTER) THEN
         OPEN(UNIT=1,FILE='RUN',STATUS='UNKNOWN')
         READ(1,9991) RUN
         READ(1,*) ICMAX
         CLOSE(1)
         IRUN = 0
         IF(RUN == 'run ') IRUN = 1
         IF(RUN == 'exit') IRUN = 2
         IF(RUN == 'rsta') IRUN = 3
         IF(IRUN == 3) THEN
            OPEN(UNIT=1,FILE='RUN',STATUS='UNKNOWN')
            WRITE(1,'(A4)') 'run'
            WRITE(1,'(I7,3X,A5)') ICMAX,'ICMAX'
            CLOSE(1)
         ENDIF
      ENDIF

 9991 FORMAT(A4)
 
      IF(PARALLEL) THEN 
         CALL MPI_BCAST( IRUN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ICMAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      ENDIF

      RETURN
      END SUBROUTINE INTERRUPT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ALLAJJ

      USE MPI

      USE MAIN_ARRAYS, ONLY : APP,JTRA,JLOC,F1RN,F1RW,F1E,IWAPP,JET

      USE NS3CO, ONLY : PARALLEL, NBLOCK

      USE INTEGERS, ONLY : MBITS, NPRO, MAXB, IPRO

      IMPLICIT NONE

      INTEGER :: MREQ1, IERRC1, IERRC2, IERRC3, ERRORCODE, IERR,
     &           MREQ1B

      MREQ1  = IWAPP(1,NBLOCK+1)

      IERRC1 = 3
      IERRC2 = 3
      IERRC3 = 3

      ALLOCATE( APP(MREQ1),STAT=IERRC1)
      ALLOCATE(JTRA(MREQ1),STAT=IERRC2)
      ALLOCATE(JLOC(MREQ1),STAT=IERRC3)

      IF(IERRC1 > 0 .OR. IERRC2 > 0 .OR. IERRC3 > 0) THEN
         WRITE(*,1022)  IPRO, REAL(MREQ1,8)/1048576.
         WRITE(45,1022) IPRO, REAL(MREQ1,8)/1048576.
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

      MREQ1B = MREQ1

      IF(PARALLEL) CALL MPI_REDUCE(MREQ1,MREQ1B,1,MPI_INTEGER,MPI_SUM,0,
     +     MPI_COMM_WORLD,IERR)

      IF(PARALLEL) CALL MPI_BCAST(MREQ1B,1,MPI_INTEGER,0,
     +     MPI_COMM_WORLD,IERR)

      IF(MREQ1B > NPRO) THEN

         CALL SETCON(JET(1),JET(MAXB/2+1),F1E,JLOC,JTRA,APP,MREQ1)
         CALL ZEROZZ(F1RN,MREQ1)
         CALL ZEROZZ(F1RW,MREQ1)
         CALL ZEROZZ(F1E,MREQ1)

      ENDIF

1022  FORMAT(
     +/'INP :Not enough memory for nonmatching connectionin process ',
     +     I3,'. Desired ',F6.2,'MB. aborting...'/)

      RETURN
      END SUBROUTINE ALLAJJ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C          
      SUBROUTINE ALLCHI

C ... Skin elements must be recognized for some Chimera routines.

      USE MPI
      USE NS3CO
      USE MAIN_ARRAYS, ONLY : IPSKIN, IBSKIN, ITSKIN, 
     &   SKINX, SKINY, SKINZ, SKINXN, SKINYN, SKINZN, ICON
      USE INTEGERS

      IMPLICIT NONE

      INTEGER :: IERR, IERR1, IERR2, IERR3, IERR4, IERR5, ERRORCODE

      ALLOCATE(IPSKIN(NPRO))

      CALL SKINEL(ICON,NBCS,IPRO,NPRO,NSKIN,IPSKIN)

      IF(ALLOCATED(SKINX)) DEALLOCATE(SKINX,SKINY,SKINZ,SKINXN,
     &                                SKINYN,SKINZN,IBSKIN,ITSKIN)

      ALLOCATE(SKINX(4,NSKIN), STAT=IERR1)
      ALLOCATE(SKINY(4,NSKIN), STAT=IERR2)
      ALLOCATE(SKINZ(4,NSKIN), STAT=IERR3)
      ALLOCATE(SKINXN(9,NSKIN),STAT=IERR1)
      ALLOCATE(SKINYN(9,NSKIN),STAT=IERR2)
      ALLOCATE(SKINZN(9,NSKIN),STAT=IERR3)
      ALLOCATE(IBSKIN(NSKIN),  STAT=IERR4)
      ALLOCATE(ITSKIN(NSKIN),  STAT=IERR5)

      IF(IERR1 > 0 .OR. IERR2 > 0 .OR. IERR3 > 0 .OR. 
     &   IERR4 > 0 .OR. IERR5 > 0) THEN
         WRITE(*,*) ' Not enough memory for skin elements ',NSKIN
         WRITE(*,*) ' Exiting...'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
         ELSE
         WRITE(45,*) 'In ALLCHI: IPSKIN, SKINX... allocated.'
      ENDIF

      CALL SKINXY(ICON,NBCS,NBLOCK,NBLOCG,MGM,SKINX,SKINY,SKINZ,
     &            SKINXN,SKINYN,SKINZN,ITSKIN,IBSKIN,NSKIN,IPSKIN,
     &            IPRO,NPRO,1)

      DEALLOCATE(IPSKIN)

      RETURN
      END SUBROUTINE ALLCHI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CPUTIM

      USE MPI
      USE NS3CO

      IMPLICIT NONE

      INTEGER, DIMENSION(8) :: VALUES
      INTEGER :: IERR

      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

C ... For bench marking
      IF (.NOT.CONVL.AND.PARALLEL)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR) 

C ... START UP TIME IN MPI-TIME CHECKING
      IF(PARALLEL) WTIMST = MPI_WTIME()
      WCAL = 0.0
      WCOM = 0.0
      IF(PARALLEL) WRITE(45,'("Starting time",E14.4)') REAL (WTIMST)
      IF(PARALLEL) WRITE(45,*) 
     +     'ICYCLE  TTOTAL(s) COMTIM(s) COMT./TTOTAL COMT./CALC'

C ... Monitor CPU-time. Activate 'CALL SECOND' if does not work

      CALL CPU_TIME(RTIME1)
      
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

      IF(.NOT. TIMEL) THEN

      IF(IOLD <= 0) THEN ! A new problem
         WRITE(13,'(A,I2,A1,I2,A1,I2,A)') '  First cycle starts at ',
     +   VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ELSE ! Continuing an old one
         WRITE(13,'(A,I7,A4,I2,A1,I2,A1,I2,A)')'  Starting from cycle ',
     +   ICYCLE,' at ',VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ENDIF

      ELSE IF(TIMEL) THEN

      IF(IOLD <= 0) THEN ! A new problem
         WRITE(13,'(A,A,F6.3,A,I2,A1,I2,A1,I2,A)') '  Time-accurate ',
     +   'simulation starts from T = ',T,' at ',
     +   VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ELSE ! Continuing an old one
         WRITE(13,'(A,A,F6.3,A,I2,A1,I2,A1,I2,A)') '  Time-accurate ',
     +   'simulation continues from T = ',T,' at ',
     +   VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ENDIF

      ENDIF ! .NOT. TIMEL

      IF(IDRXX == 1) THEN
         WRITE(13,*)    ' DROMAX is followed in convergence history'
         WRITE(4,'(/A)')' DROMAX is followed in convergence history'
      ELSEIF(IDRXX == 2) THEN
         WRITE(13,*) ' DUMAX is followed in convergence history'
      ELSEIF(IDRXX == 3) THEN
         WRITE(13,*) ' DVMAX is followed in convergence history'
      ELSEIF(IDRXX == 4) THEN
         WRITE(13,*) ' DWMAX is followed in convergence history'
      ELSEIF(IDRXX == 5) THEN
         WRITE(13,*) ' DTMAX is followed in convergence history'
      ELSEIF(IDRXX == 6) THEN
         WRITE(13,*)    ' DPMAX is followed in convergence history'
         WRITE(4,'(/A)')' DPMAX is followed in convergence history'
      ELSEIF(IDRXX == 7) THEN
         WRITE(13,*) ' DRKMAX is followed in convergence history'
      ELSEIF(IDRXX == 8) THEN
         WRITE(13,*) ' DEPMAX is followed in convergence history'
      ELSEIF(IDRXX >= 9) THEN
         WRITE(13,*) IDRXX-8,' is followed in convergence history'
      ENDIF

      RETURN
      END SUBROUTINE CPUTIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYCLES

      USE MPI

      USE CHARACTERS

      USE INTEGERS, ONLY : IPRO, MAXSB, IB, MAXB, MAXSS, MAXS2, NBCS

      USE FLIGHT,   ONLY : NGRIFL, OSKU, XCG, YCG, ZCG, ACTDISK,
     &                     IGRSLAVE, SHIPPROP, FLYOBJ, FXSP, FXTSP,
     &                     ISCALE_ACT,RTMSP,FDSP,ROTA1,ROTB1,CONEA,
     &                     THRUST,IFA,ADV,TORQUE,UTSP

      USE BLADE_VARIABLES, ONLY :  THETACOL,THCYCLON,THCYCLAT,QTH,
     &     QTH_A,QBE_A,SHAFT,QZE_S1,QBE_S1,QZE_A,QZE_D

      USE MAIN_ARRAYS

      USE NS3CO

      USE TYPE_ARRAYS , ONLY : SIXDOF

      IMPLICIT NONE

      INTEGER :: N,NGL,M,I,IG1,IQ1,IF1,IH1,IHQ1,IHT1,IR1,IC1,IP1,IG2,
     &    IG3,IF2,IH2,IHQ2,IHT2,IQ2,IQ3,IC2,IR2,IM2,ISS,NOREY,KSTATE,
     &    IGN,IHP1,IHP2,IM3,IPHASE,ICHARP,III,IHN,IGR,IERR,IGRSLA,
     &    IT1,IT2,IT3,ITN,ITP1,ITP2,ITPN,IGM,IPC2,II,JJ,KK

      REAL :: REPSLIM,FRSG,FRSRET,FXTAPU,FYTAPU,FZTAPU,
     &        MXTAPU,MYTAPU,MZTAPU
      REAL :: EKAQTH,EKAQTH_A,EKAQBE_A,EKAQZE_S1,EKAQBE_S1,
     &        EKAQZE_A,EKAQZE_D,
     &        TOKAQTH,TOKAQTH_A,TOKAQBE_A,TOKAQZE_S1,TOKAQBE_S1,
     &        TOKAQZE_A,TOKAQZE_D

      LOGICAL :: AMGDIVG,AMGDIVL,MASTER

      AMGDIVL = .FALSE. ; AMGDIVG = .FALSE.
      MASTER = IPRO == 1 
      ISCALE_ACT = 0 ! Initialize the actuator force-scaling (antakee Mosse)

C **********************************************************************
C ... BLOCK LOOP 1 = ITERATION CYCLE BEGINS
C **********************************************************************

      DO 8050 N = 1,NBLOCK

      IG1     = IG(1,N)
      IP1     = 1
      IGN     = 1

      NGL     = NPROCE(N+1,IPRO)

      IF(MULPHL) THEN
         IP1  = IG1
         IGN  = IG1 + NTOT(1,N) - 1
      ENDIF
      IF(TRANSL) THEN
         IT1  = IG1
         ITN  = IG1 + NTOT(1,N) - 1
      ENDIF
      
C ... Blockwise properties

      PR      = BLKS(NGL)%PR
      PRT     = BLKS(NGL)%PRT
      IFLUX   = BLKS(NGL)%IFLUX
      TURLIM  = BLKS(NGL)%TURLIM
      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSVEL  = BLKS(NGL)%FRSVEL
      FRSVIS  = BLKS(NGL)%FRSVIS
      T0      = BLKS(NGL)%T0
      DIFPRE  = BLKS(NGL)%DIFPRE
      CHLREF  = BLKS(NGL)%CHLREF
      IPRESC  = BLKS(NGL)%IPRESC
      FRSG    = BLKS(NGL)%FRSG
      FRSRET  = BLKS(NGL)%FRSRET


C ... NEW PRESSURES AND VELOCITIES AT THE FINE GRID LEVEL

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,'(A,I4)') ' Iteration starting for block',NGL
         WRITE(3,'(A)') ' Pressure before statim (first level)'
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,UROTC,UROT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF

C ######################################################################
      CALL STATIM(1,N,0,.TRUE.)
C ######################################################################

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' Pressure after statim (first level)'
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF

C ... TURBULENCE IS ALSO EVALUATED AT THE FINE GRID LEVEL
C ... IN THE PREVIOUS VERSION THIS WAS BEFORE THE RESIDU-CALL 

C***********************************************************************
      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'MULTI'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

      CALL TURBDI(1,N,IPRINT)

      ENDIF
C***********************************************************************
c       REWIND(927)
c        do I=1,NTOT(1,1)
c       READ(927,*)I,EPS2(I),VIST(I)
c       ENDDO

C ... FINEST (FIRST) GRID LEVEL

C *********************************************************************

C ... CALCULATE THE RESIDUAL (EXPLICIT STEP)
      
C ######################################################################
      CALL RESIDI(1,N,1,IPRINT)
C ######################################################################
     
C ... PERFORM THE TIME INTEGRATION (IMPLICIT IF CFL > .05)
     
      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' Pressure before the implicit stage (first level)'
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF
      
C ######################################################################
      CALL IMPSDI(1,N,1,AMGDIVL)
C ######################################################################

C ... UPDATE  VALUES AFTER THE FIRST STEP

      IG1   = IG(1,N)
      IQ1   = IQ(1,N)
      IF1   = JF(1,N)

      IP1 = 1
      IT1 = 1
      IF(MULPHL) IP1 = IG1          ! Beginning of Two-phase arrays  
      IF(TRANSL) IT1 = IG1

      NOREY = 1

      IF (JSCAL > 0) THEN
         NOREY = NSCAL - JSCAL + 1  ! Begining address of updated scalars
      ENDIF

      CALL FIN(TEMP(IG1),U(IG1),V(IG1),W(IG1),PDIFF(IG1),
     +     DRO(IG1),DM(IG1),DN(IG1),DW(IG1),DE(IG1),
     +     TOLD(IG1),UOLD(IG1),VOLD(IG1),WOLD(IG1),POLD(IG1),
     +     RK(IG1),REPS(IG1),DRK(IG1),DEPS(IG1),
     +     RKOLD(IG1),EPSOLD(IG1),RKLIM,EPSLIM,PRO(IP1),VAR(IP1),
     +     FI(IQ1,NOREY),DFI(IQ1,NOREY),FIOLD(IQ1,NOREY),RLOLIM,UPPLIM,
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),MAXSB,JSCAL,NSCAL,ITURB,IPRESC,
     +     TURLIM,VIS(IG1),TRANSL,TRM(IT1),BLKS,NGL,ICYCLE)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' Pressure after updating the first level'
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF

C ... scalar type updating


C ... Update and limit the pressure and temperature

      IF(MULPHL) THEN
        DO I = IG1,IGN
        PRO(I)%TSAT = MIN(PRO(I)%TSAT,BLKS(NGL)%TUPLIM)
        PRO(I)%TSAT = MAX(PRO(I)%TSAT,BLKS(NGL)%TLOLIM)
        DO IPHASE   = 1,BLKS(NGL)%NPHASE
        PRO(I)%TEMP(IPHASE) = MIN(PRO(I)%TEMP(IPHASE),BLKS(NGL)%TUPLIM)
        PRO(I)%DTEMP(IPHASE)= MIN(PRO(I)%DTEMP(IPHASE),BLKS(NGL)%TUPLIM)
        VAR(I)%X(IPHASE)    = MIN(VAR(I)%X(IPHASE),   1.)
        VAR(I)%ALFA(IPHASE) = MIN(VAR(I)%ALFA(IPHASE),1.)
        PRO(I)%DTEMP(IPHASE)= MAX(PRO(I)%DTEMP(IPHASE),BLKS(NGL)%TLOLIM)
        VAR(I)%X(IPHASE)    = MAX(VAR(I)%X(IPHASE),   0.)
        VAR(I)%ALFA(IPHASE) = MAX(VAR(I)%ALFA(IPHASE),0.)
        ENDDO
        ENDDO
      ENDIF ! MULPHL
      CALL UPLIM(TEMP(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     BLKS(NGL)%TUPLIM)
      CALL LOLIM(TEMP(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     BLKS(NGL)%TLOLIM)
      CALL UPLIM(U(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     5000.)
      CALL LOLIM(U(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     -5000.)
      CALL UPLIM(V(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     5000.)
      CALL LOLIM(V(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     -5000.)
      CALL UPLIM(W(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     5000.)
      CALL LOLIM(W(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     -5000.)
      CALL PDTOP(P(IG1),PDIFF(IG1),XC(IG1),YC(IG1),ZC(IG1),FRSDEN,
     +     FRSPRE,NTOT(1,N))
      IF(ITURB /= 8 .OR. ITURB /= 9) THEN ! 
      CALL LOLIM(REPS(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),IN,JN,KN,EPSLIM)
      ELSE IF(ITURB == 9) THEN! Spalart-Allmaras near the walls (res3c)??
      CALL LOLIM(REPS(IG1),IMAX(1,N),JMAX(1,N),
c     +     KMAX(1,N),IN,JN,KN,0.1*EPSLIM)
     +     KMAX(1,N),IN,JN,KN,EPSLIM)
      REPSLIM = 2.*TURLIM*FRSVIS
      CALL UPLIM(REPS(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),IN,JN,KN,REPSLIM)
      ENDIF ! ITURB /= 9

      IF (TRANSL) THEN
      DO I = IG1,IGN
         TRM(I)%G   = MIN(TRM(I)%G,     5.)
         TRM(I)%RET = MIN(TRM(I)%RET,3000.)
      ENDDO
      ENDIF

      KSTATE  = JSTATE(NGL,1) ! Single-phase flow

      IF(KSTATE /= 10) THEN ! LIMIT THE PRESSURE
         CALL LOLIMP(P(IG1),PDIFF(IG1),FRSDEN,FRSPRE,BLKS(NGL)%DFRSPRE,
     +    BLKS(NGL)%DFRSTEM,XC(IG1),YC(IG1),ZC(IG1),
     +    IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
      ENDIF
      
      IF (ITURB >= 24) THEN
        CALL ADDUUK(RK(IG1),FI(IG1,1),IMAX(1,N),JMAX(1,N),
     +  KMAX(1,N),MAXSB)
      ENDIF

C *********************************************************************
C ... START OF MULTIGRID CYCLE?
C**********************************************************************
         
      DO 7000 M = 2,MGRID(N)

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
         WRITE(3,*)
         WRITE(3,*) ' Recalculating the residual on level',M-1
      ENDIF

C ... NEW PRESSURES AND VELOCITIES AT THE PREVIOUS GRID LEVEL

C     CALL BOUNDI(M-1)

C ... CALCULATE VELOCITIES AND SPECIFIC INTERNAL ENERGY
C ######################################################################
      CALL STATIM(M-1,N,0,.TRUE.) 
C ######################################################################

C ######################################################################
      CALL RESIDI(M-1,N,M,IPRINT)
C ######################################################################

      IG1  = IG(M-1,N)
      IGM  = IG1 + NTOT(M-1,N) - 1
      IF1  = JF(M-1,N)
      IH1  = IH(M-1,N)
      IHQ1 = IH(M-1,N)
      IHT1 = IH(M-1,N)
      IP1  = IG(M-1,N)
      IHP1 = IH(M-1,N)
      IT1  = IG(M-1,N)
      ITP1 = IH(M-1,N)

      IF(.NOT.MULPHL) IP1  = 1
      IF(.NOT.MULPHL) IHP1 = 1
      IF(.NOT.TRANSL) IT1  = 1
      IF(.NOT.TRANSL)  ITP1 = 1
      IF(ITURB < 3 .OR. ITURB/=8)   IHT1 = 1
      IF(NSCAL == 0)  IHQ1 = 1

      IR1  = IR(M-1,N)
      IQ1  = IQ(M-1,N)
      IC1  = IC(M-1,N)
      IG2  = IG(M,N)
      IGN  = IG2 + NTOT(M,N) - 1
      IG3  = IG(M+1,N)  ! Unused
      IF2  = JF(M,N)
      IH2  = IH(M,N)
      IHQ2 = IH(M,N)
      IHT2 = IH(M,N)
      IM2  = IG(M,N)
      IHP2 = IH(M,N)
      IHN  = IHP2 + NTOT(M,N) - 1
      IT2  = IG(M,N)
      ITP2    = IH(M,N)
      ITPN    = ITP2 + NTOT(M,N) - 1

      IF(.NOT.MULPHL) IM2  = 1
      IF(.NOT.MULPHL) IHP2 = 1
      IF(.NOT.TRANSL) IT2  = 1
      IF(.NOT.TRANSL)  ITP2 = 1
      IF(ITURB < 3 .OR. ITURB/=8)   IHT2 = 1
      IF(NSCAL == 0)  IHQ2 = 1

      IQ2  = IQ(M,N)

C ... ADD THE NEW RESIDUAL AND THE FORCING FUNCTION

      IF(M > 2) THEN
      
      CALL ADDRES(ROP2H(IH1),RMP2H(IH1),RNP2H(IH1),RWP2H(IH1),EP2H(IH1),
     +     VAR(IP1),P2H(IHP1),DRO(IG1),DM(IG1),DN(IG1),DW(IG1),DE(IG1),
     +     RKP2H(IH1),EPSP2H(IH1),DRK(IG1),DEPS(IG1),FIP2H(IHQ1,1),
     +     DFI(IQ1,1),NTOT(M-1,N),MAXSB,MAXS2,NSCAL,ITURB,MULPHL,
     +     BLKS(NGL)%NPHASE,TRANSL,TRM(IT1),P2HTRM(ITP1),
     +     BLKS(NGL)%SOLUTION_TYPE)

      ENDIF
         
C ... TRANSFER OF DATA TO THE COARSE GRID LEVEL
      CALL TRANSF(TEMP(IG2),U(IG2),V(IG2),W(IG2),PDIFF(IG2),
     1      TOLD(IG2),UOLD(IG2),VOLD(IG2),WOLD(IG2),POLD(IG2),
     2      ROP2H(IH2),RMP2H(IH2),RNP2H(IH2),RWP2H(IH2),EP2H(IH2),
     3      TEMP(IG1),U(IG1),V(IG1),W(IG1),PDIFF(IG1),
     4      DRO(IG1),DM(IG1),DN(IG1),DW(IG1),DE(IG1),VOL(IG2),VOL(IG1),
     3      EPS2(IG2),VIST(IG2),EPS2(IG1),VIST(IG1),
     5      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
       CALL TRANKE(RK(IG2),REPS(IG2),RKOLD(IG2),EPSOLD(IG2),
     1      RKP2H(IH2),EPSP2H(IH2),RK(IG1),REPS(IG1),
     2      DRK(IG1),DEPS(IG1),VOL(IG2),VOL(IG1),
     3      EPS2(IG2),VIST(IG2),EPS2(IG1),VIST(IG1),
     4      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))
      ENDIF

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
         WRITE(3,*)
         WRITE(3,*) ' Transferred residuals to level',M
         CALL PRINYS(3,DROC,ROP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEC, EP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEPSC, EPSP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF

      IF(MULPHL) THEN
       CALL TRANMF(PRO(IM2),VAR(IM2),PRO(IP1),VAR(IP1),P2H(IHP2),
     2      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     3      KMAX(M-1,N),BLKS(NGL)%NPHASE)
         IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     +   TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
         WRITE(3,*)
         WRITE(3,*) ' Transferred variables to level',M
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         DO I = IG2,IGN
         F1RM(I) = VAR(I)%ALFA(IPHASE)
         F1RN(I) = VAR(I)%X(IPHASE)
         ENDDO
         F1RW(IHP2:IHN) = P2H(IHP2:IHN)%ROP2H(IPHASE)
         CALL PRINYS(3,CHAR_VAR(1,ICHARP),F1RM(IG2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RN(IG2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),F1RW(IHP2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         DO I = IG2,IGN
         F1RM(I) = PRO(I)%TEMP(IPHASE)
         F1RN(I) = VAR(I)%EVAPR(IPHASE)
         ENDDO
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(IG2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,CHAR_VAR(3,ICHARP),F1RN(IG2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         ENDDO
         ENDIF
         WRITE(3,*)
         WRITE(3,*) ' End of transferred variables to level',M
         WRITE(3,*)
      ENDIF ! MULPHL

      IF(TRANSL) THEN
         CALL TRANTR(TRM(IT2),TRM(IT1),P2HTRM(ITP2),
     1   VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))
      ENDIF 

      IF(NSCAL > 0) THEN
       CALL TRANSC(FI(IQ2,1),FIOLD(IQ2,1),FIP2H(IHQ2,1),FI(IQ1,1),
     1      DFI(IQ1,1),VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),
     2      KMAX(M,N),KMAX(M-1,N),NTOT(M-1,N),NTOT(M,N),NSCAL)
c    2      KMAX(M,N),KMAX(M-1,N),MAXSB,MAXS2,NSCAL) ! testaamati
      ENDIF 
      IF(ITURB >= 10 .AND. ITURB <= 19) THEN
       CALL TRANSO(S11(IG2,1),S11(IG1,1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXSS,6)
       CALL TRANSO(PTUR(IG2),PTUR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N), MAXB,1)
      ENDIF 

      IF(NCHIM > 0) THEN
       CALL TRANSO(ROFOR(IG2),ROFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO(RMFOR(IG2),RMFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO(RNFOR(IG2),RNFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO(RWFOR(IG2),RWFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO( EFOR(IG2), EFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO(PDFOR(IG2),PDFOR(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL TRANSO(RKSI(IG2),RKSI(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
         WRITE(3,*)
         WRITE(3,*) ' Transferred rksis to level',M
         CALL PRINYS(3,RKSIC,RKSI(IG2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF
      ENDIF

      IF(NCHIM > 0) THEN

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
       CALL  TRANSO(RKFOR(IG2),RKFOR(IG1),VOL(IG2),VOL(IG1),
     &         IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
C ... Before 14.1.2003 epsilon was not transferred!
       CALL  TRANSO(REFOR(IG2),REFOR(IG1),VOL(IG2),VOL(IG1),
     &         IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
      ENDIF

      IF(TRANSL) THEN
       DO I = IG1,IGM
         F1RM(I) = TRM(I)%GFOR
         F1RN(I) = TRM(I)%RETFOR
       ENDDO
       CALL  TRANSO(F1RM(IG2),F1RM(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       CALL  TRANSO(F1RN(IG2),F1RN(IG1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
       DO I = IG2,IGN ! Antakee mersy
         TRM(I)%GFOR   = F1RM(I)
         TRM(I)%RETFOR = F1RN(I)
       ENDDO
      ENDIF

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
            F1RM(IG1:IGM) = MPFOR(IG1:IGM)%ALFAFOR(IPHASE)
            F1RN(IG1:IGM) = MPFOR(IG1:IGM)%XFOR(IPHASE)
            F1RW(IG1:IGM) = MPFOR(IG1:IGM)%DTEMPFOR(IPHASE)
            CALL TRANSO(F1RM(IG2),F1RM(IG1),VOL(IG2),VOL(IG1),
     &           IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
            CALL TRANSO(F1RN(IG2),F1RN(IG1),VOL(IG2),VOL(IG1),
     &           IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
            CALL TRANSO(F1RW(IG2),F1RN(IG1),VOL(IG2),VOL(IG1),
     &           IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,1)
            MPFOR(IG2:IGN)%ALFAFOR(IPHASE)  = F1RM(IG2:IGN)
            MPFOR(IG2:IGN)%XFOR(IPHASE)     = F1RN(IG2:IGN)
            MPFOR(IG2:IGN)%DTEMPFOR(IPHASE) = F1RW(IG2:IGN)
         ENDDO
      ENDIF

      IF(NSCAL > 0) THEN
       CALL  TRANSO(FIFOR(IG2,1),FIFOR(IG1,1),VOL(IG2),VOL(IG1),
     1      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MAXB,NSCAL)
      ENDIF
      ENDIF ! NCHIM > 0

cccc  CALL BOUNDI(M)

C ... PRESSURES AND VELOCITIES ON THE COARSE GRID LEVEL

C ######################################################################
      CALL STATIM(M,N,2,.TRUE.)
C ######################################################################

C ... RESIDUAL ON THE COARSE GRID LEVEL (FIRST ORDER)
 
C ######################################################################
      CALL RESIDI(M,N,M,IPRINT)
C ######################################################################

C ... FORCING FUNCTION
        
      CALL FORCFU(ROP2H(IH2),RMP2H(IH2),RNP2H(IH2),RWP2H(IH2),EP2H(IH2),
     +     VAR(IM2),P2H(IHP2),DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +     RKP2H(IH2),EPSP2H(IH2),DRK(IG2),DEPS(IG2),FIP2H(IHQ2,1),
     +     DFI(IQ2,1),NTOT(M,N),MAXSB,MAXS2,NSCAL,ITURB,MULPHL,
     +     BLKS(NGL)%NPHASE,TRANSL,TRM(IT2),P2HTRM(ITP2),
     +     BLKS(NGL)%SOLUTION_TYPE)

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
         WRITE(3,*)
         WRITE(3,*) ' After CALL FORCFU on level',M
         CALL PRINYS(3, DEPSC,EPSP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF

C ... ADD THE FORCING FUNCTION AND THE RESIDUAL TOGETHER

      CALL ADDRES(ROP2H(IH2),RMP2H(IH2),RNP2H(IH2),RWP2H(IH2),EP2H(IH2),
     +     VAR(IM2),P2H(IHP2),DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +     RKP2H(IH2),EPSP2H(IH2),DRK(IG2),DEPS(IG2),FIP2H(IHQ2,1),
     +     DFI(IQ2,1),NTOT(M,N),MAXSB,MAXS2,NSCAL,ITURB,MULPHL,
     +     BLKS(NGL)%NPHASE,TRANSL,TRM(IT2),P2HTRM(ITP2),
     +     BLKS(NGL)%SOLUTION_TYPE)

C ######################################################################

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
         WRITE(3,*)
         WRITE(3,*) ' Forcing function on level',M
         CALL PRINYS(3,DROC,ROP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEC, EP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEPSC,EPSP2H(IH2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(IHP2:IHN) = P2H(IHP2:IHN)%ROP2H(IPHASE)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),F1RM(IHP2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         ENDDO
         ENDIF
         WRITE(3,*) ' Final residual on level',M
         CALL PRINYS(3,DROC,DRO(IG2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEC, DE(IG2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, DEPSC,DEPS(IG2),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(IG2:IGN) = VAR(IG2:IGN)%DX(IPHASE)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),F1RM(IG2),
     +   IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         ENDDO
         ENDIF
      ENDIF

C ... IMPLICIT SWEEP ON THE LOWER LEVEL. COURANT NUMBER = CFLL

C ######################################################################
      CALL IMPSDI(M,N,M,AMGDIVL)
C ######################################################################

C ... UPDATE THE VALUES ON THE NEXT COARSER GRID

      NOREY = 1
      IF(JSCAL > 0) THEN
      NOREY = NSCAL - JSCAL + 1  ! BEGINING ADDRESS OF UPPDATED SCALARS
      ENDIF

      CALL FIN(TEMP(IG2),U(IG2),V(IG2),W(IG2),PDIFF(IG2),
     +     DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +     TOLD(IG2),UOLD(IG2),VOLD(IG2),WOLD(IG2),POLD(IG2),
     +     RK(IG2),REPS(IG2),DRK(IG2),DEPS(IG2),
     +     RKOLD(IG2),EPSOLD(IG2),RKLIM,EPSLIM,PRO(IM2),VAR(IM2),
     +     FI(IQ2,NOREY),DFI(IQ2,NOREY),FIOLD(IQ2,NOREY),RLOLIM,UPPLIM,
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),MAXSB,JSCAL,NSCAL,ITURB,IPRESC,
     +     TURLIM,VIS(IG2),TRANSL,TRM(IT2),BLKS,NGL,ICYCLE)

C ... Update (?) and limit the pressure and temperature

ccc      CALL PDTOP(P(IG2),PDIFF(IG2),XC(IG2),YC(IG2),ZC(IG2),FRSDEN,
ccc     +     FRSPRE,NTOT(M,N))
      CALL UPLIM(TEMP(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,TUPLIM)
      CALL LOLIM(TEMP(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,TLOLIM)
      IF(MULPHL) THEN
      DO IPHASE = 1,BLKS(NGL)%NPHASE
      F1RM(IG2:IGN) = PRO(IG2:IGN)%TEMP(IPHASE)
      F1RN(IG2:IGN) = VAR(IG2:IGN)%X(IPHASE)
      F1RW(IG2:IGN) = VAR(IG2:IGN)%ALFA(IPHASE)
      CALL UPLIM(F1RM(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,BLKS(NGL)%TUPLIM)
      CALL UPLIM(F1RN(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,1.)
      CALL UPLIM(F1RW(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,1.)
      CALL LOLIM(F1RM(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,BLKS(NGL)%TLOLIM)
      CALL LOLIM(F1RN(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,0.)
      CALL LOLIM(F1RW(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,0.)
      PRO(IG2:IGN)%TEMP(IPHASE) = F1RM(IG2:IGN)
      F1RM(IG2:IGN) = PRO(IG2:IGN)%DTEMP(IPHASE)
      CALL UPLIM(F1RM(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,BLKS(NGL)%TUPLIM)
      CALL LOLIM(F1RM(IG2),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,BLKS(NGL)%TLOLIM)
      PRO(IG2:IGN)%DTEMP(IPHASE) = F1RM(IG2:IGN)
      ENDDO
      ENDIF

      IF(TRANSL) THEN
      DO I = IG2,IGN
        TRM(I)%G    = MIN(TRM(I)%G,     2.)
        TRM(I)%RET  = MIN(TRM(I)%RET,2000.)
      ENDDO
      ENDIF

c      IF(IPRESC /= 1) THEN ! LIMIT THE PRESSURE
         CALL LOLIMP(P(IG2),PDIFF(IG2),FRSDEN,FRSPRE,BLKS(NGL)%DFRSPRE,
     +    BLKS(NGL)%DFRSTEM,XC(IG2),YC(IG2),ZC(IG2),
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN)
c      ENDIF !Varo
     
      IF (ITURB >= 24) THEN
        CALL ADDUUK(RK(IG2),FI(IQ2,1),IMAX(M,N),JMAX(M,N),
     +  KMAX(M,N),MAXSB)
      ENDIF

c      CALL STATIM(M,N,.TRUE.) ! siirretty
      
7000  CONTINUE

C
C ... INTERPOLATION BACK TO THE FINER GRID LEVELS
C
       
      DO 8000 M = MGRID(N),2,-1

      IG3     = IG(M,N)
      IG2     = IG(M-1,N)
      IGN     = IG2 + NTOT(M-1,N) - 1
      IQ3     = IQ(M,N)
      IQ2     = IQ(M-1,N)
      IM3     = IG(M,N)
      IM2     = IG(M-1,N)
      IT3     = IG(M,N)
      IT2     = IG(M-1,N)
           
      CALL INTERP(TEMP(IG3),U(IG3),V(IG3),W(IG3),P(IG3),
     1     PDIFF(IG3),TOLD(IG3),UOLD(IG3),VOLD(IG3),WOLD(IG3),
     2     POLD(IG3),TEMP(IG2),U(IG2),V(IG2),W(IG2),
     3     P(IG2),PDIFF(IG2),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     4     KMAX(M-1,N),MIB(N),MIT(N),MJB(N),MJT(N),MKB(N),MKT(N))
         
      IF(MULPHL) THEN
      CALL INTERM(VAR(IG3),PRO(IG3),VAR(IG2),PRO(IG2),IMAX(M,N),
     2     JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MIB(N),MIT(N),MJB(N),
     3     MJT(N),MKB(N),MKT(N),BLKS(NGL)%NPHASE)
      IF(ICYCLE == IPRINT) THEN
         ICHARP = BLKS(NGL)%ICHAR(1)
         WRITE(3,*)' RESULTS AFTER INTERM'
         F1RM(1:NTOT(M-1,N)) = VAR(1:NTOT(M-1,N))%X(1)
         CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RM(1),
     +   IT(M-1,N),IL(M-1,N),0,IK(M-1,N),
     +   IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),NGL,M-1)
         ICHARP = BLKS(NGL)%ICHAR(2)
         WRITE(3,*)' RESULTS AFTER INTERM'
         F1RM(1:NTOT(M-1,N)) = VAR(1:NTOT(M-1,N))%X(2)
         CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RM(1),
     +   IT(M-1,N),IL(M-1,N),0,IK(M-1,N),
     +   IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),NGL,M-1)
         ICHARP = BLKS(NGL)%ICHAR(1)
         WRITE(3,*)' RESULTS AFTER INTERM'
         F1RM(1:NTOT(M-1,N)) = PRO(1:NTOT(M-1,N))%DTEMP(1)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(M-1,N),IL(M-1,N),0,IK(M-1,N),
     +   IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),NGL,M-1)
         ICHARP = BLKS(NGL)%ICHAR(2)
         WRITE(3,*)' RESULTS AFTER INTERM'
         F1RM(1:NTOT(M-1,N)) = PRO(1:NTOT(M-1,N))%DTEMP(2)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(M-1,N),IL(M-1,N),0,IK(M-1,N),
     +   IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),NGL,M-1)
      ENDIF
      ENDIF

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL INTERT(RK(IG3),REPS(IG3),RKOLD(IG3),EPSOLD(IG3),RK(IG2),
     1     REPS(IG2),IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),
     2     CMGK,CMGEPS,MIB(N),MIT(N),MJB(N),MJT(N),MKB(N),MKT(N))
       
      CALL UPLIM(TEMP(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,TUPLIM)
      CALL LOLIM(TEMP(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,BLKS(NGL)%TLOLIM)
      CALL LOLIM(RK(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,RKLIM)
      IF(ITURB /= 8 .OR. ITURB /= 9) THEN ! 
      CALL LOLIM(REPS(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,EPSLIM)
      ELSE IF(ITURB == 9) THEN! Spalart-Allmaras near the walls (res3c)??
      CALL LOLIM(REPS(IG2),IMAX(M-1,N),JMAX(M-1,N),
c     +     KMAX(M-1,N),IN,JN,KN,0.1*EPSLIM)
     +     KMAX(M-1,N),IN,JN,KN,EPSLIM)
      REPSLIM = 2.*TURLIM*FRSVIS
      CALL UPLIM(REPS(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,REPSLIM)
      ENDIF ! ITURB == 9
      ENDIF ! ITURB >= 3
          
      IF(MULPHL) THEN
      DO IPHASE = 1,BLKS(NGL)%NPHASE
      F1RM(IG2:IGN) = PRO(IG2:IGN)%TEMP(IPHASE)
      F1RN(IG2:IGN) = VAR(IG2:IGN)%X(IPHASE)
      F1RW(IG2:IGN) = VAR(IG2:IGN)%ALFA(IPHASE)
      CALL UPLIM(F1RM(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,BLKS(NGL)%TUPLIM)
      CALL UPLIM(F1RN(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,1.)
      CALL UPLIM(F1RW(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,1.)
      CALL LOLIM(F1RM(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,BLKS(NGL)%TLOLIM)
      CALL LOLIM(F1RN(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,0.)
      CALL LOLIM(F1RW(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,0.)
      PRO(IG2:IGN)%TEMP(IPHASE) = F1RM(IG2:IGN)
      VAR(IG2:IGN)%X(IPHASE)    = F1RN(IG2:IGN)
      VAR(IG2:IGN)%ALFA(IPHASE) = F1RW(IG2:IGN)
      F1RM(IG2:IGN) = PRO(IG2:IGN)%DTEMP(IPHASE)
      CALL UPLIM(F1RM(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,BLKS(NGL)%TUPLIM)
      CALL LOLIM(F1RM(IG2),IMAX(M-1,N),JMAX(M-1,N),
     +     KMAX(M-1,N),IN,JN,KN,BLKS(NGL)%TLOLIM)
      PRO(IG2:IGN)%DTEMP(IPHASE) = F1RM(IG2:IGN)

      ENDDO
      ENDIF

      IF(TRANSL) THEN
      CALL INTERTR(TRM(IT3),TRM(IT2),
     2     IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),
     3     MIB(N),MIT(N),MJB(N),MJT(N),MKB(N),MKT(N))
      DO I = IG2,IGN
        TRM(I)%G    = MIN(TRM(I)%G,     2.)
        TRM(I)%G    = MAX(TRM(I)%G,     0.)
        TRM(I)%RET  = MIN(TRM(I)%RET,2000.)
        TRM(I)%RET  = MAX(TRM(I)%RET,   0.)
      ENDDO
      ENDIF

c      IF(IPRESC /= 1) THEN ! LIMIT THE PRESSURE
         CALL LOLIMP(P(IG2),PDIFF(IG2),FRSDEN,FRSPRE,BLKS(NGL)%DFRSPRE,
     +    BLKS(NGL)%DFRSTEM,XC(IG2),YC(IG2),ZC(IG2),
     +    IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),IN,JN,KN)
c      ENDIF ! Varo

      IF(JSCAL > 0) THEN
      NOREY = NSCAL - JSCAL + 1  ! BEGINING ADDRESS OF UPPDATED SCALARS
      CALL INTERS(FI(IQ3,NOREY),FIOLD(IQ3,NOREY),FI(IQ2,NOREY),
     1     IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N),MIB(N),MIT(N),
     2     MJB(N),MJT(N),MKB(N),MKT(N),MAXSB,JSCAL,CC1,CC2)

      CALL LOLISC(FI(IQ2,1),IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),
     +      IN,JN,KN,RLOLIM,UPPLIM,MAXSB,NSCAL)
      ENDIF
     
 8000 CONTINUE

C **********************************************************************
C ... END OF MULTIGRID CYCLE
C **********************************************************************

      IG2     = IG(1,N)
      IF2     = JF(1,N)
      IH2     = IH(1,N)
      IQ2     = IQ(1,N)
      IC2     = IC(1,N)
      IR2     = IR(1,N)        
      ISS     = 1
      IM2     = 1
      IT2     = 1
      IPC2    = 1
      IF (STRESL) ISS = IG2
      IF (MULPHL) IM2 = IG2
      IF (TRANSL) IT2 = IG2
      IF(IPRESC == 1) IPC2 = IG2
c ... testing
      IG1  = IG(1,N)
C ######################################################################
c      CALL STATIM(1,N,.TRUE.) ! testing
C ######################################################################

C ... UPDATE VELOCITIES AND MOMENTUM BELOW THE SOLID SURFACES

      NGL     = NPROCE(1+N,IPRO) !Global block number
      KSTATE  = JSTATE(NGL,1)

      CALL REFLEC(U(IG2),V(IG2),W(IG2),XC(IG2),YC(IG2),ZC(IG2),
     1 IMAX(1,N),JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),
     2 JTOP(1,N),KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     3 GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,NGL,NPATCH(N),ICON(IC2),
     4 RO(IG2),P(IG2),PDIFF(IG2),RM(IG2),RN(IG2),RW(IG2),E(IG2),
     5 VIS(IG2),VIST(IG2),CH(IG2),EPS2(IG2),TEMP(IG2),C(IG2),DRDP(IG2),
     6 DRDH(IG2),RK(IR2),REPS(IR2),DDEPS(IR2),PRO(IM2),VAR(IM2),
     7 RGAS,T0,ITURB,KSTATE,IB,MAXB,RCON(IC2),UWALL(1),VWALL(1),
     8 WWALL(1),TWALL(1),IHF(1,1),TLOLIM,TUPLIM,BLKS(NGL)%SOLUTION_TYPE,
     9 MULPHL,REFLECL,TRANSL,TRM(IT2),
     1 A1XA(IG2),A1YA(IG2),A1ZA(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     2 A2XA(IG2),A2YA(IG2),A2ZA(IG2),PRC(IPC2),IPRESC)

C ... UPDATE KINETIC ENERGY OF TURBULENCE BY REYNOLDS STRESSES

      IF (ITURB >= 24) THEN
        CALL ADDUUK(RK(IG2),FI(IQ2,1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +  MAXSB)
      ENDIF

C ... TURBULENCE IS ALSO EVALUATED AT THE FINE GRID LEVEL
C ... IN THE PREVIOUS VERSION THIS WAS BEFORE THE RESIDU-CALL 

C **********************************************************************
c      CALL TURBDI(1,N,IPRINT)
C **********************************************************************

8050   CONTINUE
c **********************************************************************
C ... END OF BLOCK LOOP 1. RETURN TO MAIN PROGRAM
C **********************************************************************

C ... Check troubles in algebraic multigrid of pressure correction

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(AMGDIVL,AMGDIVG,1,MPI_LOGICAL,MPI_LOR,
     &                         MPI_COMM_WORLD,IERR)
      ELSE
         AMGDIVG = AMGDIVL
      ENDIF

      IF(MASTER .AND. AMGDIVG) THEN
       WRITE(*,*) ' Severe troubles in a convergence of pressure',
     +   ' correction. See RUN.LOG or with mpi: grep AMG RUN.LOG*'
         WRITE(*,*) ' I am finishing..'
         WRITE(13,*) ' Severe troubles in a convergence of pressure',
     +   ' correction. See above or with mpi: grep AMG RUN.LOG*'
         WRITE(13,*) ' I am finishing..'
         STOP
      ENDIF

c ... onko tama hyva kohta/pitaako tehda joka kierroksella (tutkitaan)
C ... Calls PARTICLE_FORCES routine between iteration loops

      IF (NGRIFL > 0 .AND.  .NOT.TIMEL) THEN       ! not time-accurate

c      IF (IGRID(NGL,1) >= 21 .AND. IGRID(NGL,1) <= 29 .OR.! Ship case or
c     +    IGRID(NGL,1) == 37 .OR. IGRID(NGL,1) == 38.OR.! ! 2-dof blade or
c     +    IGRID(NGL,1) >= 40 .AND. IGRID(NGL,1) < 50) THEN ! Actuator disk

         FXTAPU=0.0 ! These are needed when ACTDISK
         FYTAPU=0.0
         FZTAPU=0.0
         MXTAPU=0.0
         MYTAPU=0.0
         MZTAPU=0.0
            
      DO IGR = 1,NGRIFL
         IGRSLA = IGRSLAVE(IGR)
         
      CALL PARTICLE_FORCES(OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     +     OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,ICON,
     +     XCG(IGR),YCG(IGR),ZCG(IGR),2,10,10,0,IGR,IGRID(:,3))  ! MOV surfaces
      IF(INOUTL .AND. FLYOBJ) THEN ! thrust for the flying objects
         CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +        OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +        XCG(IGR),YCG(IGR),ZCG(IGR),2,3,5,0,IGR,IGRID(:,3))
      ENDIF                     ! INL and OUT surfaces
      
      IF(ACTDISK) THEN           
         CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +        OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +        XCG(IGR),YCG(IGR),ZCG(IGR),3,1,1,8,IGR,IGRID(:,3))! ACT disk surf
         IF(IGRSLA /= 0) THEN 
            CALL PARTICLE_FORCES(OSKU(IGRSLA)%FXT,OSKU(IGRSLA)%FYT,
     +           OSKU(IGRSLA)%FZT,OSKU(IGRSLA)%MXT,OSKU(IGRSLA)%MYT,
     +           OSKU(IGRSLA)%MZT,ICON,XCG(IGRSLA),YCG(IGRSLA),
     +           ZCG(IGRSLA),2,1,1,8,IGR,IGRID(:,3))! ACT disk surf, slave IGRID
            FXTAPU=FXTAPU+OSKU(IGRSLA)%FXT 
            FYTAPU=FYTAPU+OSKU(IGRSLA)%FYT ! FXTAPU ... MZTAPU are needed
            FZTAPU=FZTAPU+OSKU(IGRSLA)%FZT ! to summarize if more than one
            MXTAPU=MXTAPU+OSKU(IGRSLA)%MXT ! ACT disk are in use.
            MYTAPU=MYTAPU+OSKU(IGRSLA)%MYT
            MZTAPU=MZTAPU+OSKU(IGRSLA)%MZT
            OSKU(IGRSLA)%FXT=FXTAPU 
            OSKU(IGRSLA)%FYT=FYTAPU ! Final values (FXTAPU ... MZTAPU)
            OSKU(IGRSLA)%FZT=FZTAPU ! are stored in OSKU tables
            OSKU(IGRSLA)%MXT=MXTAPU
            OSKU(IGRSLA)%MYT=MYTAPU
            OSKU(IGRSLA)%MZT=MZTAPU            
            IF(IPRO == 1) THEN
               FXSP(IGR)  = OSKU(IGRSLA)%CX
               FXTSP(IGR) = OSKU(IGR)%FXT     ! Self propulsion updates
            ENDIF
            IF(PARALLEL) THEN
               CALL MPI_BCAST(FXSP,NGRIFL+10,MPI_REAL8,0,
     +              MPI_COMM_WORLD,IERR)
               CALL MPI_BCAST(FXTSP,NGRIFL+10,MPI_REAL8,0,
     +              MPI_COMM_WORLD,IERR) ! Self propulsion updates
            ENDIF
         ENDIF
         IF(PARALLEL .AND. IFA(IGR) == 8 .OR. ! MPI commands transfer data
     +        PARALLEL .AND. IFA(IGR) == 10 .OR. ! from slaves to the master
     +        PARALLEL .AND. IFA(IGR) == 12 .OR. 
     +        PARALLEL .AND. IFA(IGR) == 13) THEN 
            IF(IFA(IGR) == 12 .OR. IFA(IGR) == 13) THEN ! Patria's ACT disk
               EKAQTH    = QTH(IGR)
               EKAQTH_A  = QTH_A(IGR)
               EKAQBE_A  = QBE_A(IGR)
               CALL MPI_REDUCE(EKAQTH,TOKAQTH,1,MPI_REAL8,
     +              MPI_MAX,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQTH_A,TOKAQTH_A,1,
     +              MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQBE_A,TOKAQBE_A,1,
     +              MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,IERR)
            ELSEIF(IFA(IGR) == 8 .OR. IFA(IGR) == 10) THEN ! TQIHM
               EKAQTH    = QTH(IGR)
               EKAQTH_A  = QTH_A(IGR)
               EKAQBE_A  = QBE_A(IGR)
               EKAQZE_A  = QZE_A(IGR)            
               EKAQZE_S1 = QZE_S1(IGR)
               EKAQBE_S1 = QBE_S1(IGR)
               EKAQZE_D  = QZE_D(IGR)
               CALL MPI_REDUCE(EKAQTH,TOKAQTH,1,MPI_REAL8,
     +              MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQTH_A,TOKAQTH_A,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQBE_A,TOKAQBE_A,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)  
               CALL MPI_REDUCE(EKAQZE_A,TOKAQZE_A,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQZE_S1,TOKAQZE_S1,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQBE_S1,TOKAQBE_S1,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(EKAQZE_D,TOKAQZE_D,1,
     +              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
            ENDIF
            IF (IPRO == 1) THEN ! Data is stored when master proces
               IF (IFA(IGR) == 12 .OR. IFA(IGR) == 13) THEN ! Patria's ACT disk
                  QTH(IGR)   = TOKAQTH
                  QTH_A(IGR) = TOKAQTH_A
                  QBE_A(IGR) = TOKAQBE_A
               ELSEIF (IFA(IGR) == 8 .OR. IFA(IGR) == 10) THEN ! TQIHM
                  SHAFT(IGR)  = TOKAQTH*0.5 
                  ADV(IGR)    = TOKAQTH_A*0.5
                  THRUST(IGR) = TOKAQBE_A*0.5
                  TORQUE(IGR) = TOKAQZE_A*0.5
                  UTSP(IGR)   = TOKAQZE_D*0.5
                  OSKU(IGR)%DAMPT = TOKAQZE_S1*0.5
                  OSKU(IGR)%DAMPN = TOKAQBE_S1*0.5
               ENDIF
            ENDIF
         ENDIF ! PARALLEL .AND. IFA(IGR) == 8 .OR.  
      ENDIF                     ! IF(ACTDISK) THEN
      
      IF(SHIPPROP .AND. IGRSLA /= 0) THEN ! SHIP PROPULSION
      CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +     OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +     XCG(IGR),YCG(IGR),ZCG(IGR),2,9,9,0,IGR,IGRID(:,3)) ! ROT surfaces
      CALL PARTICLE_FORCES(OSKU(IGRSLA)%FXT,OSKU(IGRSLA)%FYT,
     +     OSKU(IGRSLA)%FZT,OSKU(IGRSLA)%MXT,OSKU(IGRSLA)%MYT,
     +     OSKU(IGRSLA)%MZT,ICON,XCG(IGRSLA),YCG(IGRSLA),ZCG(IGRSLA),
     +     2,9,9,0,IGR,IGRID(:,3)) ! ROT surfaces, slave IGRID
      ENDIF

      ENDDO 

c      ENDIF
      ENDIF ! NGRIFL > 0 .AND. .NOT.TIMEL
      RETURN
      END SUBROUTINE CYCLES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RESIDI(M,N,M2,IPR)

      USE INTEGERS,    ONLY : MAXEB,MAX11,MAXW,IPRO,MAXSB,MAXSS,IBF

      USE MAIN_ARRAYS

      USE NS3CO

      IMPLICIT NONE

      INTEGER M,N,M2,IPR,IS1,N1,NGL,IPRI,IG11,IRR2,IGL,IQL,
     &    IRL,IGC,ICH,ICQ,IE2,INTI,INTJ,INTK,IDERI,IG2,IF2,
     &    IH2,IR2,IQ2,IC2,ISS,IM2,IPC2,KSTATE,iii,LSTATE(3),
     &    IFR2,IC2M,IT2,IMCH

      REAL CFM      

      IS1 = 1
      DO 1000 N1 = 1,N-1
1000  IS1 = IS1 + NPATCH(N1)               ! CORRECT ADDRESS FOR PATCH-ARRAYS

      IF(M == 1) THEN
      DO 2000 N1 = IS1,IS1+NPATCH(N)-1
             CXB(N1)   = 0.
             CYB(N1)   = 0.
             CZB(N1)   = 0.
             DXB(N1)   = 0.
             DYB(N1)   = 0.
             DZB(N1)   = 0.
             CMXB(N1)  = 0.
             CMYB(N1)  = 0.
             CMZB(N1)  = 0.
             TOMEGA(N1)= 0.
             QTB(N1)   = 0.
             QWB(N1)   = 0.
             QHB(N1)   = 0.
             BUX(N1)   = 0.
             BUY(N1)   = 0.
             BUZ(N1)   = 0.
             VMXB(N1)  = 0.
             VMYB(N1)  = 0.
             VMZB(N1)  = 0.
2000  CONTINUE
      ENDIF

      NGL     = NPROCE(1+N,IPRO)     ! Global block number
      KSTATE  = JSTATE(NGL,1)        ! For single phase
      LSTATE(1:3) = JSTATE(NGL,1:3)  ! For resike routines
c      if(ngl == 3 .and. mulphl) kstate = jstate(ngl,2) ! viritys testikeisiin

      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSTEM  = BLKS(NGL)%FRSTEM
      FRSSIE  = BLKS(NGL)%FRSSIE
      ARTSSP  = BLKS(NGL)%ARTSSP
      FRSSSP  = BLKS(NGL)%FRSSSP
      FRSVEL  = BLKS(NGL)%FRSVEL
      TLOLIM  = BLKS(NGL)%TLOLIM
      TUPLIM  = BLKS(NGL)%TUPLIM
      IPRESC  = BLKS(NGL)%IPRESC

      CFM     = CFL
      IF(M >= 2) CFM = CFLL

      IPRI    = IPR
      IG2     = IG(M,N)

      IG11    = IG(M,N) -  IG(1,N) + 1

      IF2     = JF(M,N)
      IH2     = IH(M,N)
      IC2     = IC(M,N)
      IR2     = IR(M,N)
      IQ2     = IQ(M,N)

      IRR2    = 1
      ISS     = 1
      IGL     = 1
      IQL     = 1
      IRL     = 1
      IGC     = 1
      ICH     = 1
      IMCH    = 1
      ICQ     = 1
      IE2     = 1
      IM2     = 1
      IPC2    = 1
      IFR2    = 1
      IT2     = 1

      IF (ITURB >= 21)IRR2 = IQ2
      IF (STRESL)     ISS  = IG2
      IF (TIMEL) THEN
         IGL     = IG2
         IQL     = IQ2
         IRL     = IR2
      ENDIF

      IF (COORL)        IGC = IG2
      IF(MULPHL)        IM2 = IG2
      IF(MULPHL .AND. CHIMEL) IMCH = IG2
      IF (CHIMEL)       ICH = IG2 ! Kai pelaa?
*      IF (NCHIM > 0)   ICH = IG2
      IF(ISTRES > 0) IE2 = IG2
      IF(IPRESC == 1) IPC2= IG2
      IF(FRESUL)        IFR2 = IG2
      IF(FRESUL)        IPC2 = IG2
      IF(NCHIM > 0 .AND. NSCAL > 0) ICQ = IQ2
      IF(TRANSL)        IT2 = IG2

      INTI    = INTERI(N)
      INTJ    = INTERJ(N)
      INTK    = INTERK(N)
      IDERI   = IDER(N)
      IF(M >= 2) THEN
         INTI = 55
         INTJ = 55
         INTK = 55
      ENDIF

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID' .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'MULTI' .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

       IF(ITURB >= 3 .AND. ITURB /= 8 .AND. ITURB /= 9) ! Two-equation model
     + CALL RESIKE(N,M,M2,NGL,ISSB(1,NGL),RO(IG2),RM(IG2),RN(IG2),
     +    RW(IG2),E(IG2),P(IG2),PDIFF(IG2),U(IG2),V(IG2),W(IG2),C(IG2),
     +    CH(IG2),CP(IG2),A1(IG2),D1(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +    A2(IG2),D2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),D3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),DISTW(IG2),
     +    XCO(IG2),YCO(IG2),ZCO(IG2),DTL(IG2),
     +    UWALL(1),VWALL(1),WWALL(1),CPWALL(1),TWALL(1),QWALL(1),
     +    HFLUX(1),QWFRIC(1),TAUW1(1),TAUW2(1),SURFX(1),SURFY(1),
     +    SURFZ(1),TAUWX(1),TAUWY(1),TAUWZ(1),SURLE(1),DSURLE(1),
     +    WMFLUX(1),POROS(1),RMLOSS(1),WHSTAG(1),WTEMP(1),RSDIRX(1),
     +    RSDIRY(1),RSDIRZ(1),RBK(1),XCP(1),YCP(1),ZCP(1),IHF(1,M),
     +    DT,GAMMA,RGAS,DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +    VOL(IG2),F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),TEMP(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),
     +    ALPHA,ICP(IF2),JCP(IF2),KCP(IF2),IJMASK(IF2),
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),INTI,INTJ,INTK,IDERI,
     +    JBOT(M,N),JTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N),
     +    IK(M,N),IPRI,ICMAX,MPRINT,ICYCLE,ITIMES,T,PR,PRT,LUSGS,CFM,
     +    CXB(IS1),CYB(IS1),CZB(IS1),CMXB(IS1),CMYB(IS1),CMZB(IS1),
     +    DXB(IS1),DYB(IS1),DZB(IS1),QTB(IS1),QWB(IS1),QHB(IS1),
     +    TOMEGA(IS1),OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    ICON(IC2),NPATCH(N),XMOM,
     +    YMOM,ZMOM,IBOT(M,N),ITOP(M,N),KBOT(M,N),KTOP(M,N),LAMIN(N),
     +    IFLUX,UBI(IF2),VBI(IF2),WBI(IF2),UBJ(IF2),VBJ(IF2),WBJ(IF2),
     +    UBK(IF2),VBK(IF2),WBK(IF2),UTI(IF2),VTI(IF2),WTI(IF2),
     +    UTJ(IF2),VTJ(IF2),WTJ(IF2),UTK(IF2),VTK(IF2),WTK(IF2),
     +    UROT(IG2),VROT(IG2),WROT(IG2),RCON(IC2),
     +    RK(IG2),REPS(IG2),DDEPS(IG2),DRK(IG2),DEPS(IG2),SRK(IG2),
     +    SEPS(IG2),PTUR(IG2),FUN1(IG2),F1RK(IG2),F1EPS(IG2),TTS(IG2),
     +    BIJ(IE2,1),WIR(IE2,1),MAXEB,ISTRES,
     +    TURLIM,RKLIM,EPSLIM,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     +    CHLREF,DRDH(IG2),DRDP(IG2),KSTATE,LSTATE,ITURB,IDIS,KOVER,
     +    E0REF,T0REF,PSIGSC,PSIGS2,
     +    FI(IQ2,1),F1FI(IQ2,1),DFI(IQ2,1),SFI(IQ2,1),MAXSB,NSCAL,KSCAL,
     +    HAT1(IG2),HAT2(IG2),HAT3(IG2),HAT4(IG2),RLOLIM,PROD(IRR2),
     +    SPI(IRR2),DIF(IRR2),DIS(IRR2),VVIS(IRR2),FWLL(IRR2),
     +    JRDIF,JRDIS,JRPRE,S11(ISS,1),MAXSS,STRESL,SOURL,TIMEL,GRAVIL,
     +    IEPSMA,
     +    ROLE2(IGL),RMLE2(IGL),RNLE2(IGL),RWLE2(IGL),ELE2(IGL),
     +    RKLE2(IRL),EPSLE2(IRL),FILE2(IQL,1),
     +    XLE2(IGC),YLE2(IGC),ZLE2(IGC),
     +    ROLE3(IGL),RMLE3(IGL),RNLE3(IGL),RWLE3(IGL),ELE3(IGL),
     +    RKLE3(IRL),EPSLE3(IRL),FILE3(IQL,1),
     +    XLE3(IGC),YLE3(IGC),ZLE3(IGC),W12(IG11),SIJ(IG11),GRADT(IG11),
     +    GRADK(IG11),GREPS(IG11),MAX11,
     +    NCHIM,ROFOR(ICH),RMFOR(ICH),RNFOR(ICH),RWFOR(ICH),EFOR(ICH),
     +    PDFOR(ICH),RKFOR(ICH),REFOR(ICH),MPFOR(IMCH),
     +    RKSI(ICH),FIFOR(ICQ,1),
     +    QMFIN,QVFIN,QMEIN,QMFOUT,QVFOUT,QMEOUT,QMXIN,QMXOUT,QMYIN,
     +    QMYOUT,QMZIN,QMZOUT,CMXIN,CMYIN,CMZIN,CMXOUT,CMYOUT,CMZOUT,
     +    IROTCO,FULLNS,IPRESC,ZZZ,MAXW,BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,
     +    BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,
     +    BOUNA2,IBF,PRO(IM2),VAR(IM2),BLKS,MULPHL,PRC(IPC2),TURCOR,
     +    XXTRAL,FRESUL,IFSBC,F1H(IFR2),WAVEH(1),IWAVEB(IFR2),FREDIF,
     +    IGRID(NGL,2),BUX(IS1),BUY(IS1),BUZ(IS1),VMXB(IS1),VMYB(IS1),
     +    VMZB(IS1),FRSSIE,VTRAN(IG2),INTERTU,RE,RNUT(IR2),VORT(IG2,1),
     +    SHEAR(IG2,1),STRAIN(IG2),PSEUCO,RJK2,RJK4,TURDESL,INCHIML,
     +    ENTROPY_FIX,QGFIN,QGEIN,QGFOUT,QGEOUT,
     +    TUR_FRS_SOURCE,TRANSL,TRM(IT2),BOUNG,BOUNRET,PLE2(IGL),
     +    PLE3(IGL),IDIMSG,JDIMSG,SURFA,SURFT,SURF2X,SURF2Y,SURF2Z,
     +    SURFMX,SURFMY,SURFMZ,CDIFF,VELLAP(IG2),QSAS(IG2),IUPPT(NGL),
     +    BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2,SURFPX(1),
     +    SURFPY(1),SURFPZ(1),QLFIN,QLEIN,QLFOUT,QLEOUT,UTAUM(1),
     +    EPSOLD(IG2))
C ... 2. LAST LINE CONTAINS INTERPOLATED VALUES FROM THE CHIMERA TO THE 
C     BACKGROUND

       IF(ITURB < 3 .OR. ITURB == 8 .OR. ITURB == 9)  ! Algebraic turbulence model
     + CALL RESIDU(N,M,M2,NGL,ISSB(1,NGL),RO(IG2),RM(IG2),RN(IG2),
     +    RW(IG2),
     +    E(IG2),P(IG2),PDIFF(IG2),U(IG2),V(IG2),W(IG2),C(IG2),CH(IG2),
     +    CP(IG2),A1(IG2),D1(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +    A2(IG2),D2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),D3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),DISTW(IG2),
     +    XCO(IG2),YCO(IG2),ZCO(IG2),DTL(IG2),
     +    UWALL(1),VWALL(1),WWALL(1),CPWALL(1),TWALL(1),QWALL(1),
     +    HFLUX(1),QWFRIC(1),TAUW1(1),TAUW2(1),SURFX(1),SURFY(1),
     +    SURFZ(1),TAUWX(1),TAUWY(1),TAUWZ(1),SURLE(1),DSURLE(1),
     +    WMFLUX(1),POROS(1),WHSTAG(1),WTEMP(1),RSDIRX(1),RSDIRY(1),
     +    RSDIRZ(1),RBK(1),XCP(1),YCP(1),ZCP(1),IHF(1,M),
     +    DT,GAMMA,RGAS,E0REF,T0REF,
     +    DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),VOL(IG2),
     +    F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),TEMP(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),
     +    ALPHA,ICP(IF2),JCP(IF2),KCP(IF2),IJMASK(IF2),
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),INTI,INTJ,INTK,IDERI,
     +    JBOT(M,N),JTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N),
     +    IK(M,N),IPRI,ICMAX,MPRINT,ICYCLE,ITIMES,T,PR,PRT,LUSGS,CFM,
     +    CXB(IS1),CYB(IS1),CZB(IS1),CMXB(IS1),CMYB(IS1),CMZB(IS1),
     +    DXB(IS1),DYB(IS1),DZB(IS1),QTB(IS1),QWB(IS1),QHB(IS1),
     +    TOMEGA(IS1),OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    ICON(IC2),NPATCH(N),XMOM,
     +    YMOM,ZMOM,IBOT(M,N),ITOP(M,N),KBOT(M,N),KTOP(M,N),LAMIN(N),
     +    IFLUX,UBI(IF2),VBI(IF2),WBI(IF2),UBJ(IF2),VBJ(IF2),WBJ(IF2),
     +    UBK(IF2),VBK(IF2),WBK(IF2),UTI(IF2),VTI(IF2),WTI(IF2),
     +    UTJ(IF2),VTJ(IF2),WTJ(IF2),UTK(IF2),VTK(IF2),WTK(IF2),
     +    UROT(IG2),VROT(IG2),WROT(IG2),RCON(IC2),
     +    RK(IR2),REPS(IR2),DDEPS(IR2),DRK(IR2),DEPS(IR2),SRK(IR2),
     +    SEPS(IR2),PTUR(IR2),F1RK(IG2),F1EPS(IG2),
     +    TURLIM,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     +    DRDH(IG2),DRDP(IG2),KSTATE,JSTATE,ITURB,PSIGSC,PSIGS2,
     +    FI(IQ2,1),F1FI(IQ2,1),DFI(IQ2,1),SFI(IQ2,1),MAXSB,NSCAL,
     +    HAT1(IG2),HAT2(IG2),HAT3(IG2),HAT4(IG2),RLOLIM,TIMEL,GRAVIL,
     +    ROLE2(IGL),RMLE2(IGL),RNLE2(IGL),RWLE2(IGL),ELE2(IGL),
     +    RKLE2(IRL),EPSLE2(IRL),FILE2(IQL,1),
     +    XLE2(IGC),YLE2(IGC),ZLE2(IGC),
     +    ROLE3(IGL),RMLE3(IGL),RNLE3(IGL),RWLE3(IGL),ELE3(IGL),
     +    RKLE3(IRL),EPSLE3(IRL),FILE3(IQL,1),
     +    XLE3(IGC),YLE3(IGC),ZLE3(IGC),W12(IG11),SIJ(IG11),GRADT(IG11),
     +    GRADK(IG11),GREPS(IG11),MAX11,
     +    NCHIM,ROFOR(ICH),RMFOR(ICH),RNFOR(ICH),RWFOR(ICH),EFOR(ICH),
     +    PDFOR(ICH),RKSI(ICH),FIFOR(ICQ,1),RKFOR(ICH),REFOR(ICH),
     +    MPFOR(IMCH),FULLNS,QMFIN,QVFIN,QMEIN,
     +    QMFOUT,QVFOUT,QMEOUT,QMXIN,QMXOUT,QMYIN,QMYOUT,QMZIN,QMZOUT,
     +    CMXIN,CMYIN,CMZIN,CMXOUT,CMYOUT,CMZOUT,IPRESC,ZZZ,MAXW,BOUNMF,
     +    BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,
     +    BOUNEP,BOUNFI,BOUNBI,BOUNA1,BOUNA2,IBF,PRO(IM2),VAR(IM2),
     +    BLKS,MULPHL,PRC(IPC2),XXTRAL,FRESUL,IFSBC,F1H(IFR2),WAVEH(1),
     +    IWAVEB(IFR2),FREDIF,IGRID(NGL,2),BUX(IS1),BUY(IS1),BUZ(IS1),
     +    VMXB(IS1),VMYB(IS1),VMZB(IS1),FRSSIE,RNUT(IR2),VORT(IG2,1),
     +    SHEAR(IG2,1),STRAIN(IG2),PSEUCO,RJK2,RJK4,TURDESL,KOVER,
     +    INCHIML,ENTROPY_FIX,QGFIN,QGEIN,QGFOUT,QGEOUT,TRANSL,TRM(IT2),
     +    BOUNG,BOUNRET,PLE2(IGL),PLE3(IGL),IDIMSG,JDIMSG,CDIFF,
     +    VELLAP(IG2),QSAS(IG2),IUPPT(NGL),
     +    BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2,SURFPX(1),
     +    SURFPY(1),SURFPZ(1),QLFIN,QLEIN,QLFOUT,QLEOUT,UTAUM(1),
     +    EPSOLD(IG2))

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'SOLID') THEN
  
       CALL RESIHS(N,M,M2,NGL,ISSB(1,NGL),RO(IG2),RM(IG2),RN(IG2),
     +    RW(IG2),E(IG2),P(IG2),PDIFF(IG2),U(IG2),V(IG2),W(IG2),C(IG2),
     +    CH(IG2),A1(IG2),D1(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +    A2(IG2),D2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),D3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),DISTW(IG2),
     +    XCO(IG2),YCO(IG2),ZCO(IG2),DTL(IG2),
     +    UWALL(1),VWALL(1),WWALL(1),CPWALL(1),TWALL(1),QWALL(1),
     +    HFLUX(1),QWFRIC(1),TAUW1(1),TAUW2(1),SURFX(1),SURFY(1),
     +    SURFZ(1),TAUWX(1),TAUWY(1),TAUWZ(1),SURLE(1),DSURLE(1),
     +    WMFLUX(1),POROS(1),WHSTAG(1),WTEMP(1),RSDIRX(1),RSDIRY(1),
     +    RSDIRZ(1),RBK(1),XCP(1),YCP(1),ZCP(1),IHF(1,M),
     +    DT,GAMMA,RGAS,DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +    VOL(IG2),F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),TEMP(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),ALPHA,
     +    ICP(IF2),JCP(IF2),KCP(IF2),IJMASK(IF2),
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),INTI,INTJ,INTK,IDERI,
     +    JBOT(M,N),JTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N),
     +    IK(M,N),IPRI,ICMAX,MPRINT,ICYCLE,ITIMES,T,PR,PRT,LUSGS,CFM,
     +    CXB(IS1),CYB(IS1),CZB(IS1),CMXB(IS1),CMYB(IS1),CMZB(IS1),
     +    DXB(IS1),DYB(IS1),DZB(IS1),QTB(IS1),QWB(IS1),QHB(IS1),
     +    TOMEGA(IS1),OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    ICON(IC2),NPATCH(N),XMOM,
     +    YMOM,ZMOM,IBOT(M,N),ITOP(M,N),KBOT(M,N),KTOP(M,N),LAMIN(N),
     +    IFLUX,UBI(IF2),VBI(IF2),WBI(IF2),UBJ(IF2),VBJ(IF2),WBJ(IF2),
     +    UBK(IF2),VBK(IF2),WBK(IF2),UTI(IF2),VTI(IF2),WTI(IF2),
     +    UTJ(IF2),VTJ(IF2),WTJ(IF2),UTK(IF2),VTK(IF2),WTK(IF2),
     +    UROT(IG2),VROT(IG2),WROT(IG2),RCON(IC2),
     +    RK(IG2),REPS(IG2),DDEPS(IG2),DRK(IG2),DEPS(IG2),SRK(IG2),
     +    SEPS(IG2),PTUR(IG2),FUN1(IG2),F1RK(IG2),F1EPS(IG2),TTS(IG2),
     +    BIJ(IE2,1),WIR(IE2,1),MAXEB,ISTRES,
     +    TURLIM,RKLIM,EPSLIM,FRSDEN,FRSPRE,FRSVEL,FRSVIS,T0,DIFPRE,
     +    CHLREF,DRDH(IG2),DRDP(IG2),KSTATE,ITURB,IDIS,KOVER,
     +    E0REF,T0REF,PSIGSC,PSIGS2,
     +    FI(IQ2,1),F1FI(IQ2,1),DFI(IQ2,1),SFI(IQ2,1),MAXSB,NSCAL,KSCAL,
     +    HAT1(IG2),HAT2(IG2),HAT3(IG2),HAT4(IG2),RLOLIM,PROD(IRR2),
     +    SPI(IRR2),DIF(IRR2),DIS(IRR2),VVIS(IRR2),FWLL(IRR2),
     +    JRDIF,JRDIS,JRPRE,S11(ISS,1),MAXSS,STRESL,SOURL,TIMEL,GRAVIL,
     +    IEPSMA,
     +    ROLE2(IGL),RMLE2(IGL),RNLE2(IGL),RWLE2(IGL),ELE2(IGL),
     +    RKLE2(IRL),EPSLE2(IRL),FILE2(IQL,1),
     +    XLE2(IGC),YLE2(IGC),ZLE2(IGC),
     +    ROLE3(IGL),RMLE3(IGL),RNLE3(IGL),RWLE3(IGL),ELE3(IGL),
     +    RKLE3(IRL),EPSLE3(IRL),FILE3(IQL,1),
     +    XLE3(IGC),YLE3(IGC),ZLE3(IGC),W12(IG11),SIJ(IG11),GRADT(IG11),
     +    GRADK(IG11),GREPS(IG11),MAX11,
     +    NCHIM,ROFOR(ICH),RMFOR(ICH),RNFOR(ICH),RWFOR(ICH),EFOR(ICH),
     +    PDFOR(ICH),RKFOR(ICH),REFOR(ICH),RKSI(ICH),FIFOR(ICQ,1),
     +    QMFIN,QMEIN,QMFOUT,QMEOUT,QMXIN,QMXOUT,QMYIN,QMYOUT,QMZIN,
     +    QMZOUT,CMXIN,CMYIN,CMZIN,CMXOUT,CMYOUT,CMZOUT,
     +    IROTCO,FULLNS,IPRESC,ZZZ,MAXW,BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,
     +    BOUNE,BOUNP,BOUNPD,BOUNT,BOUNRK,BOUNEP,BOUNFI,BOUNBI,IBF,
     +    PRO(IM2),VAR(IM2),BLKS,MULPHL,RLDIST(IG2))
      ENDIF

C ... Store the explicit residuals
       
      IF(M == 1) THEN

      CALL RESL2(VOL(IG2),DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +     DFI(IQ2,1),DRK(IR2),DEPS(IR2),VAR(IM2),RESI(1,N),RKSI(ICH),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),ITURB,NSCAL,MAXSB,NRESI,
     +     RKLIM,EPSLIM,MULPHL,NGL,M,M2,PRC(IPC2))
      ENDIF

      RETURN
      END SUBROUTINE RESIDI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPSDI(M,N,M2,AMGDIVL)

      USE MAIN_ARRAYS

      USE NS3CO

      USE CHARACTERS, ONLY: DROC

      USE INTEGERS,    ONLY : MAXW,MAW,IPRO,MAXSB,ITERAD,ITERAC,MCYCAM,
     &    IBF

      IMPLICIT NONE

      INTEGER M,N,NGL,IGL,IQL,IRL,IGC,ICH,ICQ,INTI,INTJ,INTK,
     &    IG2,IF2,IH2,IR2,IQ2,IC2,IM2,IPC2,ii,M2,IP1,IGN,IT2

      REAL CFM

      LOGICAL :: AMGDIVL

C      Antakee mersu. Values from chim.f to RKSIDI, now test with zeroes
      CALL ZEROZZ(F1RK,NTOT(M,N))

C ... A CALLING ROUTINE FOR THE IMPLICIT STAGE
C ... THE SAME ROUTINE IS USED FOR NEW TURBULENCE MODELS

      NGL     = NPROCE(1+N,IPRO) !Global block number
       
      IPRESC  = BLKS(NGL)%IPRESC
      CFM     = CFL
      IF (M >= 2) CFM = CFLL
       
      IG2     = IG(M,N)
      IF2     = JF(M,N)
      IH2     = IH(M,N)
      IR2     = IR(M,N)
      IQ2     = IQ(M,N)
      IC2     = IC(M,N)
      IM2     = 1
      IPC2    = 1
      IGL     = 1
      IQL     = 1
      IRL     = 1
      IGC     = 1
      ICH     = 1
      ICQ     = 1   
      IT2     = 1

      IF(MULPHL) IM2 = IG2  

      IF (TIMEL) THEN
         IGL     = IG2
         IQL     = IQ2
         IRL     = IR2
      ENDIF

      IF (COORL)      IGC = IG2
      IF(IPRESC == 1) IPC2= IG2
      IF (NCHIM > 0) ICH = IG2
      IF(NCHIM > 0 .AND. NSCAL > 0) ICQ = IQ2
      IF(TRANSL)      IT2 = IG2


      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN
      
      IF(BLKS(NGL)%IPRESC == 1) THEN

       CALL IMPSEG(NGL,M,CFM,SMAX,TEMP(IG2),A1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),D1(IG2),D2(IG2),D3(IG2),
     +    DTL(IG2),DT,DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),DRK(IR2),
     +    DEPS(IR2),RO(IG2),U(IG2),V(IG2),W(IG2),E(IG2),C(IG2),RK(IR2),
     +    REPS(IR2),CH(IG2),CP(IG2),PRO(IM2),VAR(IM2),VOL(IG2),VIS(IG2),
     +    EPS2(IG2),VIST(IG2),ITURB,JRDIS,JRPRE,IMAX(M,N),JMAX(M,N),
     +    KMAX(M,N),IBOT(M,N),ITOP(M,N),JBOT(M,N),JTOP(M,N),KBOT(M,N),
     +    KTOP(M,N),ICYCLE,IPRINT,IT(M,N),IL(M,N),IK(M,N),IDI1(N),
     +    IDI2(N),IDI3(N),NTOT(M,N),PR,PRT,XC(IG2),YC(IG2),ZC(IG2),
     +    TIMEL,MCYCAM,MGRIDA(N),ITERAD,ITERAC,NCHIM,ZZZ,HAT1,HAT3,HAT4,
     +    F1R,F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,PRC(IPC2),LUSGS,P(IG2),
     +    OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    UROT(IG2),VROT(IG2),WROT(IG2),IDI1(N),IDI2(N),IDI3(N),CFL,
     +    FRSDEN,FRSVEL,RKMAX,CSIMPS,SRK(IR2),SEPS(IR2),PTUR(IR2),KOVER,
     +    PDIFF(IG2),DRDP(IG2),DRDH(IG2),ALFAP,FI(IQ2,1),DFI(IQ2,1),
     +    SFI(IQ2,1),MAXSB,NSCAL,MULPHL,BLKS, 
     +    NPATCH(N),ICON(IC2),IHF(1,M),XCP(1),YCP(1),ZCP(1),
     +    XCO(IG2),YCO(IG2),ZCO(IG2),FPRINTL,FRSPRE,XFC(IG2),YFC(IG2),
     +    ZFC(IG2),SDI(IG2),AMGDIVL,ICONV3,RKSI(IG2),PDFOR(ICH),
     +    PSIGSC,PSIGS2,TRANSL,TRM(IT2),BLKS(NGL)%IFLUX,CDIFF,CDIFFT,
     +    FREDIF)

      ELSE

      CALL IMPS(NGL,M,LUSGS,RO(IG2),RM(IG2),RN(IG2),RW(IG2),E(IG2),
     +    P(IG2),PDIFF(IG2),U(IG2),V(IG2),W(IG2),C(IG2),CP(IG2),
     +    DRDP(IG2),DRDH(IG2),A1(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +    A2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3(IG2),A3XA(IG2),
     +    A3YA(IG2),A3ZA(IG2),D1(IG2),D2(IG2),D3(IG2),DTL(IG2),DT,GAMMA,
     +    FRSDEN,FRSVEL,FRSPRE,DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     +    VOL(IG2),VIS(IG2),EPS2(IG2),VIST(IG2),ITURB,KOVER,IMAX(M,N),
     +    JMAX(M,N),KMAX(M,N),IBOT(M,N),ITOP(M,N),JBOT(M,N),JTOP(M,N),
     +    KBOT(M,N),KTOP(M,N),ICYCLE,IPRINT,IT(M,N),IL(M,N),IK(M,N),CFM,
     +    THETA,IDI1(N),IDI2(N),IDI3(N),NTOT(M,N),PR,PRT,OHMI(IG2),
     +    RK(IR2),REPS(IR2),DRK(IR2),DEPS(IR2),
     +    SRK(IR2),SEPS(IR2),PTUR(IR2),XC(IG2),YC(IG2),ZC(IG2),CSIMPS,
     +    ICP(IF2),JCP(IF2),KCP(IF2),LAMIN(N),JRDIS,JRPRE,
     +    UROT(IG2),VROT(IG2),WROT(IG2),
     +    OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    TURLIM,RKLIM,EPSLIM,
     +    FI(IQ2,1),DFI(IQ2,1),SFI(IQ2,1),MAXSB,NSCAL,JRIMP,TIMEL,
     +    IPRESC,PSEUCO,RKMAX,HAT1(IG2),
     +    HAT2(IG2),HAT3(IG2),HAT4(IG2),F1R(IG2),F1RM(IG2),F1RN(IG2),
     +    F1RW(IG2),F1E(IG2),ICON(IC2),NPATCH(N),IFLUX,INTI,INTJ,INTK,
     +    ITERMA,TOLER,ROFOR(ICH),RMFOR(ICH),RNFOR(ICH),RWFOR(ICH),
     +    EFOR(ICH),PDFOR(ICH),RKFOR(ICH),REFOR(ICH),RKSI(ICH),
     +    FIFOR(ICQ,1),NCHIM,ZZZ,MAXW,MAW,IDIS,XMASSB,PRO(IM2),VAR(IM2),
     +    BLKS,MULPHL,TRANSL,TRM(IT2),CDIFF,CDIFFT)
      ENDIF

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
c     +BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

       CALL IMPSMF(NGL,M,CFM,SMAX,TEMP(IG2),A1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),D1(IG2),D2(IG2),D3(IG2),
     +    DTL(IG2),DT,DRO(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),DRK(IR2),
     +    DEPS(IR2),RO(IG2),U(IG2),V(IG2),W(IG2),E(IG2),C(IG2),RK(IR2),
     +    REPS(IR2),CH(IG2),CP(IG2),PRO(IM2),VAR(IM2),VOL(IG2),VIS(IG2),
     +    EPS2(IG2),VIST(IG2),ITURB,JRDIS,JRPRE,IMAX(M,N),JMAX(M,N),
     +    KMAX(M,N),IBOT(M,N),ITOP(M,N),JBOT(M,N),JTOP(M,N),KBOT(M,N),
     +    KTOP(M,N),ICYCLE,IPRINT,IT(M,N),IL(M,N),IK(M,N),IDI1(N),
     +    IDI2(N),IDI3(N),NTOT(M,N),PR,PRT,XC(IG2),YC(IG2),ZC(IG2),
     +    TIMEL,MCYCAM,MGRIDA(N),ITERAD,ITERAC,NCHIM,ZZZ,HAT1,HAT3,HAT4,
     +    F1R,F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,PRC(IPC2),LUSGS,P(IG2),
     +    OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +    UROT(IG2),VROT(IG2),WROT(IG2),IDI1(N),IDI2(N),IDI3(N),CFL,
     +    FRSDEN,FRSVEL,RKMAX,CSIMPS,SRK(IR2),SEPS(IR2),PTUR(IR2),KOVER,
     +    PDIFF(IG2),DRDP(IG2),DRDH(IG2),ALFAP,FI(IQ2,1),DFI(IQ2,1),
     +    SFI(IQ2,1),MAXSB,NSCAL,MULPHL,BLKS,
     +    NPATCH(N),ICON(IC2),IHF(1,M),XCP(1),YCP(1),ZCP(1),
     +    XCO(IG2),YCO(IG2),ZCO(IG2),FPRINTL,FRSPRE,XFC(IG2),YFC(IG2),
     +    ZFC(IG2),SDI(IG2),AMGDIVL,ICONV3,RKSI(IG2),PDFOR(ICH),
     +    PSIGSC,PSIGS2,TRANSL,TRM(IT2),BLKS(NGL)%IFLUX,CDIFF,CDIFFT,
     +    FREDIF)

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'SOLID') THEN

       CALL IMPSOL(NGL,M,CFM,SMAX,TEMP(IG2),A1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),D1(IG2),D2(IG2),D3(IG2),
     +    DTL(IG2),DT,DE(IG2),RO(IG2),CH(IG2),CP(IG2),VOL(IG2),VIS(IG2),
     +    EPS2(IG2),VIST(IG2),ITURB,IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +    IBOT(M,N),ITOP(M,N),JBOT(M,N),JTOP(M,N),KBOT(M,N),KTOP(M,N),
     +    ICYCLE,IPRINT,IT(M,N),IL(M,N),IK(M,N),IDI1(N),IDI2(N),IDI3(N),
     +    NTOT(M,N),PR,PRT,XC(IG2),YC(IG2),ZC(IG2),TIMEL,MCYCLE,
     +    MGRIDA(N),ITERMA,ITERHA,NCHIM,ZZZ,HAT1,HAT2,HAT3,HAT4,F1R,
     +    F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,DRO(IG2),FPRINTL,ICONV3) 

      ELSE

         WRITE(*,*) ' SOLUTION TYPE =',BLKS(NGL)%SOLUTION_TYPE,' does',
     +              ' not exist. Exiting..'
         WRITE(4,*) ' SOLUTION TYPE =',BLKS(NGL)%SOLUTION_TYPE,' does',
     +              ' not exist. Exiting..'
         WRITE(13,*)' SOLUTION TYPE =',BLKS(NGL)%SOLUTION_TYPE,' does',
     +              ' not exist. Exiting..'
         STOP

      ENDIF
    
      IF(FRESUL .AND. M2 == 1 .AND. IFSBC == 4 .AND. IGRID(NGL,2) >= 5
     +     .AND. IGRID(NGL,2) <= 9)THEN
         IF(IFSBC == 1) SURLE(1:IBF) = WAVEH(1:IBF)  
         CALL IMPFRE(NGL,M,PRC(IG2),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,
     +        NTOT(M,N),ICYCLE,IPRINT,IT(M,N),IL(M,N),IK(M,N),NPATCH(N),
     +        ICON(IC2),IHF(1,M),SURLE(1),DSURLE(1),VOL(IG2),RO(IG2),
     +        WAVEH(1),IWAVEB(IG2),DWMAX,FLODWH,WHMAX,WHMIN,A1(IG2),
     +        A2(IG2),A3(IG2),XCP(1),YCP(1),ZCP(1),FREDIF,PDIFF(IG2),
     +        POLD(IG2),D1(IG2),D2(IG2),D3(IG2),DTL(IG2),CHLREF,SUMDWH,
     +        CFLFRE,CFL,CFLL,DTWMIN,DTWMAX,DWMV,U(IG2),V(IG2),W(IG2),
     +        F1RM(IG2),F1RN(IG2),F1RW(IG2),FPRINTL,IGRID(NGL,2))
      ENDIF

      RETURN
      END SUBROUTINE IMPSDI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TURBDI(M,N,IPR)

      USE INTEGERS,    ONLY : MAXW,MAXEB,MAX11,IPRO,MAXSB,MAXSS

      USE MAIN_ARRAYS

      USE NS3CO

      IMPLICIT NONE

      INTEGER M,N,IPR,IPRI,NGL,IG11,IRR2,IE2,IDERI,INTI,INTJ,INTK,
     &    IG2,IF2,IH2,IQ2,IC2,IR2,ISS,IT2,IM2

      IPRI    = IPRINT

C ... A CALLING ROUTINE FOR THE EVALUATION OF TURBULENT VISCOSITIES

      NGL     = NPROCE(1+N,IPRO) !Global block number
      IG2     = IG(M,N)
      IF2     = JF(M,N)
      IH2     = IH(M,N)
      IQ2     = IQ(M,N)
      IC2     = IC(M,N)
      IR2     = IR(M,N)        !  NOT USED
      IC2     = IC(M,N)
      ISS     = 1
      IG11    = IG(M,N) -  IG(1,N) + 1   ! Oli NBMAX
      IT2     = 1
      IF (STRESL) ISS = IG2
      IRR2    = 1
      IF (ITURB >= 21) IRR2 = IQ2
      IE2     = 1
      IF(ISTRES > 0) IE2 = IG2
      IDERI   = IDER(N)
      IF(TRANSL) IT2 = IG2
      IM2     = 1
      IF(MULPHL) IM2 = IG2  

      IF(ITURB < 3 .AND. ITURB /= 0) THEN
       CALL TURBML(NGL,M,RO(IG2),
     +    U(IG2),V(IG2),W(IG2),C(IG2),CH(IG2),A1(IG2),D1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),D2(IG2),A2XA(IG2),A2YA(IG2),
     +    A2ZA(IG2),A3(IG2),D3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),DT,GAMMA,DRO(IG2),
     +    DM(IG2),DN(IG2),DW(IG2),DE(IG2),VOL(IG2),   ! 34
     +    F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),TEMP(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),
     +    ALPHA,ICP(IF2),JCP(IF2),KCP(IF2),IJMASK(IF2),   ! 49
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),ICYCLE,INTI,INTJ,INTK,IDERI,
     +    JBOT(M,N),JTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N), ! 65
     +    IK(M,N),IPRI,PR,PRT,CXB(N),CYB(N),CZB(N),CMXB(N),
     +    CMYB(N),CMZB(N),DXB(N),DYB(N),DZB(N),TOMEGA(N),XMOM,
     +    YMOM,ZMOM,IBOT(M,N),ITOP(M,N),KBOT(M,N),KTOP(M,N),LAMIN(N),  ! 87
     +    IFLUX,UBI(IF2),VBI(IF2),WBI(IF2),UBJ(IF2),VBJ(IF2),WBJ(IF2),
     +    UBK(IF2),VBK(IF2),WBK(IF2),UTI(IF2),VTI(IF2),WTI(IF2),
     +    UTJ(IF2),VTJ(IF2),WTJ(IF2),UTK(IF2),VTK(IF2),WTK(IF2),
     +    UROT(IG2),VROT(IG2),WROT(IG2),W12(IG11),SIJ(IG11),MAX11,
     +    TURLIM,FRSDEN,FRSPRE,FRSVEL,T0,
     +    DRDH(IG2),DRDP(IG2),ITURB,PSIGSC,PSIGS2,
     +    FI(IQ2,1),F1FI(IQ2,1),DFI(IQ2,1),SFI(IQ2,1),MAXSB,NSCAL,
     +    HAT1(IG2),HAT2(IG2),HAT3(IG2),HAT4(IG2),S11(ISS,1),MAXSS,
     +    STRESL,ZZZ(1),MAXW,VORT(IG2,1),SHEAR(IG2,1),STRAIN(IG2),
     +    RKSI(IG2),VELLAP(IG2))

      ELSEIF(ITURB >= 3) THEN ! .AND. ITURB /= 8) THEN
       
       CALL TURBEQ(NGL,M,RO(IG2),
     +    U(IG2),V(IG2),W(IG2),C(IG2),CH(IG2),A1(IG2),D1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),D2(IG2),A2XA(IG2),A2YA(IG2),
     +    A2ZA(IG2),A3(IG2),D3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),DISTW(IG2),DT,GAMMA,DRO(IG2),
     +    DM(IG2),DN(IG2),DW(IG2),DE(IG2),VOL(IG2),
     +    F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),TEMP(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),
     +    ALPHA,ICP(IF2),JCP(IF2),KCP(IF2),IJMASK(IF2),
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),ICYCLE,INTI,INTJ,INTK,IDERI,
     +    JBOT(M,N),JTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N),
     +    IK(M,N),IPRI,PR,PRT,CXB(N),CYB(N),CZB(N),CMXB(N),
     +    CMYB(N),CMZB(N),DXB(N),DYB(N),DZB(N),TOMEGA(N),
     +    OMEGA(NPROCE(1+N,IPRO)),XMOM,
     +    YMOM,ZMOM,IBOT(M,N),ITOP(M,N),KBOT(M,N),KTOP(M,N),LAMIN(N),
     +    IFLUX,
     +    UROT(IG2),VROT(IG2),WROT(IG2),W12(IG11),SIJ(IG11),MAX11,
     +    RK(IG2),REPS(IG2),DDEPS(IG2),DRK(IG2),DEPS(IG2),SRK(IG2),
     +    SEPS(IG2),PTUR(IG2),F1RK(IG2),F1EPS(IG2),TTS(IG2),
     +    BIJ(IE2,1),WIR(IE2,1),MAXEB,ISTRES,
     +    TURLIM,FRSVIS,RKLIM,EPSLIM,FRSDEN,FRSPRE,FRSVEL,T0,
     +    DRDH(IG2),DRDP(IG2),ITURB,IDIS,KOVER,IEPSMA,PSIGSC,
     +    PSIGS2,FI(IQ2,1),MAXSB,NSCAL,KSCAL,JRDIS,JRPRE,
     +    HAT1(IG2),HAT2(IG2),HAT3(IG2),HAT4(IG2),
     +    S11(ISS,1),MAXSS,STRESL,NPATCH(N),ICON(IC2),IROTCO,ZZZ,MAXW,
     +    VTRAN(IG2),RNUT(IG2),VORT(IG2,1),SHEAR(IG2,1),STRAIN(IG2),
     +    RKSI(IG2),TURDESL,TRANSL,TRM(IT2),MULPHL,TUR_MULPHL,
     +    PRO(IM2),VAR(IM2),VELLAP(IG2),QSAS(IG2),FUN1(IG2),BLKS)

c      ELSEIF(ITURB >= 24) THEN  !ei toimi profiililla ollenkaan
      ELSEIF(ITURB >= 2400) THEN ! MITA TEKEE?
       CALL TURRSM(NGL,M,U(IG2),V(IG2),W(IG2),RO(IG2),A1(IG2),A1XA(IG2),
     +    A1YA(IG2),A1ZA(IG2),A2(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     +    A3(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +    XC(IG2),YC(IG2),ZC(IG2),VOL(IG2),D1(IG2),D2(IG2),D3(IG2),
     +    F1R(IG2),F1RM(IG2),F1RN(IG2),F1RW(IG2),F1E(IG2),
     +    VIS(IG2),EPS2(IG2),VIST(IG2),OHMI(IG2),
     +    ICP(IF2),JCP(IF2),KCP(IF2),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +    ICYCLE,IBOT(M,N),ITOP(M,N),JBOT(M,N),JTOP(M,N),KBOT(M,N),
     +    KTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IT(M,N),IL(M,N),IK(M,N),
     +    IPRI,PR,PRT,
     +    UBI(IF2),VBI(IF2),WBI(IF2),UBJ(IF2),VBJ(IF2),WBJ(IF2),
     +    UBK(IF2),VBK(IF2),WBK(IF2),UTI(IF2),VTI(IF2),WTI(IF2),
     +    UTJ(IF2),VTJ(IF2),WTJ(IF2),UTK(IF2),VTK(IF2),WTK(IF2),
     +    RK(IG2),REPS(IG2),DDEPS(IG2),PTUR(IG2),TURLIM,ITURB,FI(IQ2,1),
     +    F1FI(IQ2,1),MAXSB,NSCAL,STRESL,ZZZ,MAXW,IDER(N))

      ELSEIF(ITURB == 0) THEN
         CALL ADJUIN(EPS2(IG2),1.,IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN)
         CALL ADJUIN(VIST(IG2),0.,IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN)
      ENDIF

      RETURN
      END SUBROUTINE TURBDI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYCOUT(ICYCTI,FIRSTC,FIRSTOUTP)

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : IPRO,NMOV,MAXB

      USE MAIN_ARRAYS

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: N,M,IG1,IGN,ICYCTI,NGL,NS,IC77,ICHARP,IPHASE,I
      INTEGER :: IERR, NMOVG
      LOGICAL :: FIRSTC, FIRSTOUTP

C **********************************************************************
C ... BLOCK LOOP 1 = WRITE ON THE OUTPUT FILE
C **********************************************************************

      IF(.NOT.TIMEL .AND. ICYCLE == IPRINT .OR.
     + TIMEL .AND. MOD(ICYCLE,IPRINT) == 0) THEN ! ACCORDING TO CYCLES
C    + TIMEL .AND. ITIMES == IPRINT)        THEN ! ACCORDING TO STEPS
      DO 8190 N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1 = IG(1,N)
         IGN    = IG1 + NTOT(1,N)
         WRITE(3,*)
         WRITE(3,*) ' Iteration cycle',ICYCLE,' is finished'
         CALL PRINYS(3,PC,P(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CC,C(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,DRDPC,DRDP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,DRDHC,DRDH(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CPC,CP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(BLKS(NGL)%SOLUTION_TYPE /= 'MULTI') THEN
         CALL PRINYS(3,UC,U(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,VC,V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,WC,W(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,RM(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RNC,RN(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RWC,RW(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSEIF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
         DO IPHASE = 1,1 !BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         ZZZ(1:NTOT(1,N)) = VAR(IG1:IGN)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),ZZZ(1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         ZZZ(1:NTOT(1,N)) = VAR(IG1:IGN)%V(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),ZZZ(1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         ZZZ(1:NTOT(1,N)) = VAR(IG1:IGN)%W(IPHASE)
         CALL PRINYS(3,CHAR_VAR(35,ICHARP),ZZZ(1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF
         CALL PRINYS(3,ROC,RO(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,EC,E(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(ITURB >= 3) THEN
              CALL PRINYS(3,RKC,RK(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              IF(ITURB /= 8 .AND. ITURB /= 9) THEN
              CALL PRINYS(3,REPSC,REPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              ELSE IF (ITURB == 9) THEN
              CALL PRINYS(3,RNUTC,REPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              ENDIF
C             CALL PRINYS(3,EPS2C,EPS2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
C    +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF
         IF(MULPHL) THEN
         F1RM(IG1:IGN) = PRO(IG1:IGN)%PSAT
         F1RN(IG1:IGN) = PRO(IG1:IGN)%TSAT
         CALL PRINYS(3,PSATC,F1RM(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,TSATC,F1RN(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         DO I = IG1,IGN
         F1R(I)  = VAR(I)%EVAPR(IPHASE)
         F1RM(I) = VAR(I)%ALFA(IPHASE)
         F1RN(I) = PRO(I)%DTEMP(IPHASE)
         F1RW(I) = PRO(I)%HSAT(IPHASE)
         F1E(I)  = PRO(I)%E(IPHASE)
         F1RK(I) = PRO(I)%TAUF(IPHASE)
         ENDDO
         CALL PRINYS(3,CHAR_VAR(3,ICHARP),F1R(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_VAR(1,ICHARP),F1RM(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_PH(2,ICHARP), F1RN(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_PH(9,ICHARP), F1RW(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_PH(12,ICHARP),F1E(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_PH(13,ICHARP),F1RK(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDDO
      ENDIF ! MULPHL
C ... FOR SCALAR EQ PPR 14.2

              DO 8180 NS = 1,NSCAL
              CALL PRINYS(3,FIC(NS),FI(IG1,NS),IT(1,N),IL(1,N),0,
     +        IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
8180          CONTINUE

              CALL PRINYS(3,DROC,DRO(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              IF(BLKS(NGL)%SOLUTION_TYPE /= 'MULTI') THEN
              CALL PRINYS(3,DMC,DM(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              CALL PRINYS(3,DNC,DN(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              CALL PRINYS(3,DWC,DW(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              ELSEIF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
              DO IPHASE = 1,BLKS(NGL)%NPHASE
              ICHARP = BLKS(NGL)%ICHAR(IPHASE)
              ZZZ(IG1:IGN) = VAR(IG1:IGN)%DM(IPHASE)
              CALL PRINYS(3,CHAR_VAR(30,ICHARP),ZZZ(IG1),IT(1,N),
     +        IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              ZZZ(IG1:IGN) = VAR(IG1:IGN)%DN(IPHASE)
              CALL PRINYS(3,CHAR_VAR(31,ICHARP),ZZZ(IG1),IT(1,N),
     +        IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              ZZZ(IG1:IGN) = VAR(IG1:IGN)%DW(IPHASE)
              CALL PRINYS(3,CHAR_VAR(32,ICHARP),ZZZ(IG1),IT(1,N),
     +        IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF
              CALL PRINYS(3,DEC,DE(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
              CALL PRINYS(3,DKC,DRK(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
              CALL PRINYS(3,DEPSC,DEPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF
              DO 8181 NS = 1,NSCAL
              CALL PRINYS(3,DFIC(NS),DFI(IG1,NS),IT(1,N),IL(1,N),0,
     +        IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
8181       CONTINUE
8190  CONTINUE

C ... END OF BLOCK LOOP 1
           IPRINT = IPRINT + KP  ! SO-CALLED NORMAL OUTPUT,
      ENDIF

C **********************************************************************
C ... WRITE ITERATION HISTORY ON PLOT FILES
C **********************************************************************

      IF(.NOT.TIMEL) IC77 = ICYCLE
      IF(TIMEL)      IC77 = ITIMES ! SO-CALLED NORMAL OUTPUT
       
C ... MOVIE FILES

      IF(PARALLEL) THEN 
         CALL MPI_ALLREDUCE(NMOV,NMOVG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)
      ELSE
         NMOVG = NMOV
      ENDIF

 
C ... Create MOVIE directory for animations if it does not exist

      IF(NMOVG >= 1) THEN
         CALL EXECUTE_COMMAND_LINE("mkdir -p MOVIE")
      ENDIF
     
C **********************************************************************
*      IF(MOD(IC77,MPRINT) == 0 .AND. NMOVG >= 1) CALL MOVIE
      IF(MOD(IC77,MPRINT) == 0 .AND. NMOVG >= 1) THEN
         CALL OUTP(2,FIRSTOUTP)
         FIRSTOUTP = .FALSE.
      ENDIF
C **********************************************************************

      RETURN
      END SUBROUTINE CYCOUT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INP

      USE MPI

      USE CHARACTERS

      USE TYPE_ARRAYS

      USE CONSTANTS,   ONLY : EPS,PII

      USE INTEGERS,    ONLY : MAXB,NB,NBGG,IPRO,MBPRO,NPRO,MGM,
     &    MAX11,NBCS,NPPV,MAXW,MAW,NCON,MAXMP

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,RO,RM,RN,RW,E,DRO,DM,
     &    DN,DW,DE,TOLD,UOLD,VOLD,WOLD,POLD,OHMI,U,V,W,P,PDIFF,
     &    EPS2,VIST,C,TEMP,F1R,F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,HAT1,
     &    HAT2,HAT3,HAT4,VIS,DTL,CP,CH,DRDP,DRDH,RKSI,BLANK,W12,
     &    GRADT,GRADK,GREPS,SIJ,IW,IC,JF,NTOT,ZZZ,XCO,YCO,ZCO,IG,
     &    NPROCE,IMAXG,JMAXG,KMAXG,ICON,NPATCH,MGRID,AAA,XCOL,XCOG,
     &    NLOCAL,RCON,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,ICP,JCP,KCP,
     &    IJMASK,IL,IK,IT,IWAPP,NCPAT,NPNUM,OMEGA,JET,PRO,VAR,XGRI,
     &    YGRI,ZGRI,IGRID,RNUT,FTR2,ROLD

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: I,M,N,IC1,IC2,NBG,MSTP,ISTP,JSTP,KSTP,ISTRID,JSTRID,
     &    KSTRID,L,IS,IG1,IG2,IERR,ITME,IWMAX,IPMAX,NCG,NTOT1,
     &    KU1,KU2,KU3,KU4,KU5,KU6,IMEM,IERRCODE,ERRORCODE,III

C ... THIS SUBROUTINE READS GRID, MAKES INDEX MANIPULATIONS AND
C ... CREATES CONNECTION ARRAYS

      INTEGER IX(NB),IY(NB),IZ(NB),IXG(NBGG),IYG(NBGG),IZG(NBGG),
     +     IWAXIS(6)
      LOGICAL PERCH2

C **********************************************************************
C ... READ, SEND AND RECEIVE THE GRID WHEN MPI IS USED
  
      IF(IOLD > 0 .AND. IPRO == 1) THEN
         OPEN(28,FILE='XYD.BIN',FORM='UNFORMATTED')
      ENDIF
         
      IF(IOLD <= 0) THEN
      IF(PARALLEL) THEN 
         CALL GRMPI(XCO,YCO,ZCO,IX,IY,IZ,IXG,IYG,IZG,LEVEL,
     +        GRILEN,IG,IPRO,NPROCE,MBPRO,NPRO,NBLOCK,MGM)
      ELSE
C ... READ GRID IN WHEN IOTY=2 AND MULTIPROCESSOR ARE USED OR ONE PROCESS
C ... MODE IS IN USE
         CALL GRIMP(XCO,YCO,ZCO,IX,IY,IZ,LEVEL,GRILEN,
     +        IG,IPRO,NPROCE,MBPRO,NPRO,NBLOCK,MGM,IN,JN,KN)
      ENDIF ! PARALLEL
      
      ELSE  ! IOLD > 0

         CALL READXY(IPRO,NBLOCG,NLOCAL,NPNUM,NPROCE,MBPRO,
     +     IMAXG,JMAXG,KMAXG,IMAX,JMAX,KMAX,MGM,NBLOCK,IG,
     +     IXG,IYG,IZG,IX,IY,IZ,XCO,YCO,ZCO)
      ENDIF
       
      IF(IOLD > 0 .AND. IPRO == 1) CLOSE(28) ! XYD.BIN

C **********************************************************************

C ... Next loop re-initializes the memory from density to end of MAXB 
C ... size arrays.
      
C     DO 2101 I=1,MREQ-(IP-1)
C        RO(I)  = 0.0

C 2101 CONTINUE ! Sorry Patrik, but let's do this explicitly...

      DO 2101 I = 1,MAXB
      RO(I)     = 0.
      RM(I)     = 0.
      RN(I)     = 0.
      RW(I)     = 0.
      E(I)      = 0.
      DRO(I)    = 0.
      DM(I)     = 0.
      DN(I)     = 0.
      DW(I)     = 0.
      DE(I)     = 0.
      TOLD(I)   = 0.
      UOLD(I)   = 0.
      VOLD(I)   = 0.
      WOLD(I)   = 0.
      POLD(I)   = 0.
      ROLD(I)   = 0.
      OHMI(I)   = 0.
      U(I)      = 0.
      V(I)      = 0.
      W(I)      = 0.
      P(I)      = 0.
      PDIFF(I)  = 0.
      EPS2(I)   = 0.
      VIST(I)   = 0.
      C(I)      = 0.
      TEMP(I)   = 0.
      F1R(I)    = 0.
      F1RM(I)   = 0.
      F1RN(I)   = 0.
      F1RW(I)   = 0.
      F1E(I)    = 0.
      F1RK(I)   = 0.
      F1EPS(I)  = 0.
      HAT1(I)   = 0.
      HAT2(I)   = 0.
      HAT3(I)   = 0.
      HAT4(I)   = 0.
      VIS(I)    = 0.
      DTL(I)    = 0.
      CP(I)     = 0.
      CH(I)     = 0.
      DRDP(I)   = 0.
      DRDH(I)   = 0.
      RKSI(I)   = 0.
      BLANK(I)  = 0.
      RNUT(I)   = 0.
 2101 CONTINUE

      DO 2102 I = 1,3*MAX11
      W12(I)  = 0.
      GRADT(I)= 0.
      GRADK(I)= 0.
      GREPS(I)= 0.
 2102 CONTINUE

      DO 2103 I = 1,6*MAX11
      SIJ(I)  = 0.
 2103 CONTINUE

C **********************************************************************

C ... ARRAY POINTERS FOR THE BEGINNING

      IW(1,1,1)  = 1
      IC(1,1)    = 1

C ... ADJUST ICON ARRAY ON THE GRID LEVEL USED IN THE CALCULATION
C ... (PATCH LOOP INSIDE)

      CALL SETPAT(ICON,NPATCH,NBCS,LEVEL,NBLOCK,KMAX,MGM,NPPV, !PPR
     +            IGRID(:,1),NPROCE,IPRO)

C ... ADJUST ICON ARRAY ON THE COARSE GRID LEVELS AND ADDRESSES
C ... (IC ARRAY) BLOCK LOOP 1 BEGINS

      DO 1949 M = 1,MGM
      DO 1948 N = 1,NBLOCK
         IC1 = IC(M,N)
         IF(M == 1) IC2 = IC(1,1)
         IF(M >  1) IC2 = IC(M-1,1)

         CALL IPATCH(ICON(IC1),IC,NPATCH,MGRID,NBCS,NBLOCK,ICON(IC2),
     +               M,N,MGM,KMAX,RCON(IC1),RCON(IC2))

1948  CONTINUE
      IF(M < MGM) IC(M+1,1) = IC(M,NBLOCK) + NPATCH(NBLOCK)*IC9
1949  CONTINUE
      
C ... END OF BLOCK LOOP 1

C ... PRINT OUT ICON ARRAY TO THE MEMORY FILE

      DO 1946 M = 1,MGM
         IC1  = IC(M,1)
         CALL PRIPAT(ICON(IC1),IC,NPATCH,MGRID,NBCS,NBLOCK,M,MGM,IC9)
1946  CONTINUE

C **** USAGE OF ANCIENT SOLID WALL ARRAYS ******************************

      DO 1950 N = 1,NBLOCK   ! BLOCK LOOP 2 BEGINS
      IC1     = IC(1,N)
      CALL SETBOT(ICON(IC1),NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     2 IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),KTOP(1,N),N)
1950  CONTINUE               ! END OF BLOCK LOOP 2

C **********************************************************************
C ... CHECK PROBLEM DIMENSIONS. BLOCK LOOP 3 BEGINS
            
      DO 1100 N = 1,NBLOCK
      NBG = NPROCE(1+N,IPRO)
      WRITE(45,9876) IX(N),IY(N),IZ(N),N
      MSTP      = 2**(LEVEL-1)
      ISTP      = MSTP
      JSTP      = MSTP
      KSTP      = MSTP
      IF(IX(N) == 3) ISTP = 1
      IF(IY(N) == 3) JSTP = 1
      IF(IZ(N) <= 3) KSTP = 1   !PPR
      IF(N <= NBLOCK) THEN
         IF(((IX(N)-1)/ISTP) /= IMAX(1,N) .OR.
     +        ((IY(N)-1)/JSTP) /= JMAX(1,N) .OR.
     +        ((IZ(N)-1)/KSTP) /= KMAX(1,N)) THEN
            WRITE (*,4000)
            WRITE (*,6001) NBG
            WRITE (*,6002) IMAX(1,N),JMAX(1,N),KMAX(1,N)
            WRITE (*,6003) (IX(N)-1)/ISTP,(IY(N)-1)/JSTP,(IZ(N)-1)/KSTP
            STOP
      ENDIF
      ENDIF
1100  CONTINUE

C ... END OF BLOCK LOOP 3
C ***********************************************************************
 4000 FORMAT('THE DIMENSIONS OF THE GRID FILE AND THE INPUT FILE DO NOT'
     &     ,' MATCH.'/'CHECK YOUR FILES AND RERUN THE CODE!'/
     &     'ABORTING ...')
 6001 FORMAT('  Block No. ',I3)
 6002 FORMAT('  Dimensions from the INPUT file: (',
     +          I3,',',I3,',',I3,')')
 6003 FORMAT('  Dimensions from the grid  file: (',
     +          I3,',',I3,',',I3,')')
 9876 FORMAT('IMAXP1 = ',I4,' JMAXP1 = ',I4,' KMAXP1 = ',I4,' N = ',I4)

C *****************************************************************

C***********************************************************************
      IF(.NOT.TIMEL .AND. IOLD > 0) THEN
        IPRINT = ICYCLE + KP
      ENDIF
      IF(TIMEL) THEN
        IPRINT = KP   ! KPRINT is number of printed cycle
        IBOUT  = 1
      ENDIF
         
C ... BLOCK LOOP 4 BEGINS

      DO 7600 N = 1,NBLOCK

      ISTRID     = IMAX(1,N) + 2*IN
      JSTRID     = JMAX(1,N) + 2*JN
      KSTRID     = KMAX(1,N) + 2*KN
      NTOT(1,N)  = ISTRID*JSTRID*KSTRID

      IC1        = IC(1,N)
      L          = JF(1,N)

      CALL SETKCP(ICON(IC1),NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     3     ICP(L),JCP(L),KCP(L),IJMASK(L),N)

C ... INDEXES FOR THE COARSE GRID LEVELS
      
      DO 7000 M = 2,MGRID(N)
      IC1     = IC(M,N)
      CALL IND3C(IMAX(M,N),JMAX(M,N),KMAX(M,N),NTOT(M,N),IT(M,N),
     2   IL(M,N),IK(M,N),IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),IT(M-1,N),
     3   IL(M-1,N),IK(M-1,N))
      CALL IND3F2(IBOT(M,N),ITOP(M,N),IBOT(M-1,N),ITOP(M-1,N))
      CALL IND3F2(JBOT(M,N),JTOP(M,N),JBOT(M-1,N),JTOP(M-1,N))
7000  CALL IND3F2(KBOT(M,N),KTOP(M,N),KBOT(M-1,N),KTOP(M-1,N))

      DO 7500 M = 2,MGRID(N)
         IC1     = IC(M,N)
         L       = JF(M,N)
      CALL SETKCP(ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     3 ICP(L),JCP(L),KCP(L),IJMASK(L),N)
7500  CONTINUE
7600  CONTINUE

C ... END OF BLOCK LOOP 4

C *******************************************************************
C ... EXTRAPOLATE CORNER POINT COORDINATES AND FIND PATCH CORNERS
      
      IS        = 1
      DO 7650 N = 1,NBLOCK              !  BLOCK LOOP 5 BEGINS
      IG1    = IG(1,N)
      IC1    = IC(1,N)
         CALL GRIDEX(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(1,N),JMAX(1,N),
     +        KMAX(1,N),ICON(IC1),NPATCH(N),XCOL(1,1,IS),IN,JN,KN)
         IS  = IS + NPATCH(N)         
7650  CONTINUE                          !  END OF BLOCK LOOP 5
C *******************************************************************

C ... Solve the patch orientations in case of connectivity BC or
C ... sliding BC and the cyclicity matrices in case of cyclic BC
            
      IF(IOLD <= 0) THEN
      DO I=1,NBLOCG
         F1RM(I) = 0.0
      ENDDO
      ELSE
      DO I=1,NBLOCG
         F1RM(I) = ROTANG(I)
      ENDDO
      ENDIF
               
      WRITE(45,*)
      WRITE(45,'(A)') 'Grid angles before PATCHO:'
      WRITE(45,*)
      WRITE(45,9470)
      WRITE(45,9475)
      WRITE(45,9480)

      DO 5010 NCG = 1,NBLOCG    
      WRITE(45,9500) NCG,OMEGA(NCG),ROTANG(NCG)*180./PII,ROTANG(NCG),
     +   ROTAT(NCG)*180./PII,ROTAT(NCG)
5010  CONTINUE
      WRITE(45,*)

c9470  FORMAT('BLOCK','  OMEGA (1/s) ','  Current angle (ROTANG)      ',
c     +       '     Initial angle  (ROTAT)    ')
9470  FORMAT(' BLOCK','  OMEGA (1/s) ','  Current angle (ROTANG)      ',
     +       '     Initial angle  (ROTAT)    ')
c9475  FORMAT(22X,'(degrees)',4X,'(radians)',11X,'(degrees)',5X,
c     +           '(radians)')
9475  FORMAT(23X,'(degrees)',4X,'(radians)',11X,'(degrees)',5X,
     +           '(radians)')
9480  FORMAT(80('='))
c9500  FORMAT(I3,1X,F10.2,6X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)
9500  FORMAT(I6,1X,F10.2,4X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)

      CALL PATCHO(F1RN,F1R,XCOL,XCOG,AAA,NBCS,ICON,NPATCH,F1RM,NPROCE,
     &     IPRO,MBPRO,CHLREF)

C ... Boundary condition patches are extended on all grid levels. 

*      DO N = 1,NBLOCK
*         DO M = 1,MGRID(N)
*            IC1 = IC(M,N)
*            CALL RISEC2(ICON(IC1),NPATCH(N),IMAXG,JMAXG,KMAXG,
*     &                  NPROCE,MBPRO,NPRO,IPRO,M)
*         ENDDO
*      ENDDO

C *******************************************************************
C ... SET UP IW -ARRAY

      IWMAX   = 0
      IPMAX   = 0
      DO 7643 N = 1,NBLOCK
         DO 7643 M = 1,MGRID(N)
         IC1     = IC(M,N)
         CALL IND3W3(IW,ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),
     +        KMAX(M,N),IN,JN,KN,M,N,MGRID(N),NBLOCK,NB,NPPV,ITURB,
     +        ISTRES,NSCAL,ITME,IWMAX,IPMAX,MULPHL,NPHASES,TRANSL,
     +        TWO_FLUIDL)
 7643 CONTINUE

      WRITE(45,*)'*****************************************************'
      WRITE(45,*)' number of group connective variables is:',ITME
      WRITE(45,*)'*****************************************************'
      WRITE(45,*)' ARRAY POINTERS FOR IW:'
      WRITE(45,*)'IW(I,NBLOCK,MULTIGRID) '
      DO 7653 N = 1,NBLOCK
         DO M = 1,MGRID(N)
            WRITE(45,312)N,M,(IW(I,N,M),I=1,NPATCH(N))
         ENDDO
 7653 CONTINUE
 312  FORMAT(' IW(I,',I4,',',I2') =',50I8) 

      WRITE(45,*) 'IWMAX = ',IWMAX, ' MAXW =',MAXW
      MAXW = MAX0(MAXW,IWMAX)
C **********************************************************************

      WRITE(45,*) 'NBCS*2 = ',NBCS*2, ' MAXW =',MAXW
      WRITE(45,*)
      MAXW  = MAX0(MAXW,NBCS*2)
C **********************************************************************

      WRITE(45,*)'ARRAY SIZE (MAXW)= ',MAXW
      WRITE(45,*)'Maximum message length for this process is ',IPMAX
      WRITE(45,*)
*      MAXW = MAX0(MAXW,NPRO*(NRESI+9)/2*43)
      MAXW = MAX0(MAXW,NPRO*(NRESI+14)/2*43)
      WRITE(45,*)'Derired MAXW and MAW sizes are ',MAXW,MAW
      WRITE(45,*) 'Allocating ZZZ for MAXW = ',MAXW

      ALLOCATE(ZZZ(MAXW))
      ZZZ = 0.

      WRITE(45,*) 'ZZZ allocated. MAXW = ', MAXW,'MAW = ',MAW,
     +     ' MAXW/46 =',MAXW/46 ! FLUMPH usage
C **********************************************************************
      IF(IPMAX > MAXB) THEN

      WRITE(45,*)'Array size (IPMAX)= ',IPMAX
      WRITE(45,*)'Array size of F1R (MAXB)= ',MAXB,' => increased',
     + ' to IPMAX for connections'
      WRITE(45,*)

      WRITE(45,*) 'Rellocating F1R ...'

      DEALLOCATE(F1R)
      ALLOCATE(F1R(IPMAX))
      F1R = 0.
      
      WRITE(45,*) 'F1R reallocated. IPMAX = ', IPMAX
      WRITE(45,*)
      ENDIF
C **********************************************************************

      IF(PARALLEL .AND. MULPHL) THEN

         IMEM = 0

         DO NCG = 1,NBLOCG

            IMEM = MAX(IMEM,(IMAXG(NCG)+2*IN)
     &                     *(JMAXG(NCG)+2*JN)
     &                     *(KMAXG(NCG)+2*KN))
         ENDDO

         ALLOCATE(FTR2(IMEM),STAT=IERRCODE)

         IF (IERRCODE > 0) THEN
            WRITE(*,1022) IPRO,REAL(4*IMEM,4)/1048576.
            IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
 1022       FORMAT(/'MAIN :      Not enough memory in process ',I3,
     +     '. Desired ',F6.2,'MB. aborting...'/)
         ELSE
         WRITE(45,*) 'FTR2 allocated. IMEM = ', IMEM
         WRITE(45,*)
         ENDIF

         FTR2 = 0.

      ENDIF

      CALL SUPERC(ICON(IC(1,1)),NBCS,NCON,NCPAT)

      DO 7644 M = 1,1
      IC1  = IC(M,1)
      CALL PRIPAT(ICON(IC1),IC,NPATCH,MGRID,NBCS,NBLOCK,M,MGM,IC9)
7644  CONTINUE    !  MIKA ERO ALLA OLEVAAN?
C
C      WRITE(45,*) '  FINAL FORM OF ICON'
C      WRITE(45,*)
C      CALL PRIICO(ICON,NBCS)
             
C ... Check connections matching each other. If not make a transpolation
C ... table. This is for non matching boundaries

C ... F1RN = JLOC
C ... F1RW = JTRA
C ... F1E  = APP
          
      CALL CONEX(XCO,YCO,ZCO,F1R,F1RM,JET(1),JET(MAXB/2+1),F1E,
     &           1,1,NBLOCK,MAXB)

      IF (IWAPP(1,NBLOCK+1) /= 1) THEN
      WRITE(45,*)' ARRAY nonmatching connective POINTERS FOR IWAPP:'
      DO N = 1,NBLOCK+1
         WRITE(45,311) N,(IWAPP(I,N),I=1,6)
      ENDDO
      ENDIF
 311  FORMAT(' IWAPP(I,',I2,') =',6I7) 

C ... make check- periodic boundary is inlet/outlet or side.
C ... ICON(20,N) 
C ... = 0 = normal connection
C ... = 1 = side periodic
C ... = 2 = inlet periodic
C ... = 3 = outlet periodic
      
      PERCHL = .FALSE.
      PERCH2 = .FALSE.

      DO N = 1,NBLOCK
         
         IG2 = IG(1,N)
         IC2 = IC(1,N)
         CALL PERCHE(IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ALPHA,BETA,
     +        NPATCH(N),ICON(IC2),XCO(IG2),YCO(IG2),ZCO(IG2),N,PERCHL)
      ENDDO
      
      PERCH2 = PERCHL
      IF(PARALLEL) CALL MPI_ALLREDUCE(PERCH2,PERCHL,1,MPI_LOGICAL,
     &         MPI_LOR,MPI_COMM_WORLD,IERR)

      IF(PERCHL) THEN
          WRITE(45,*) ' There are periodic patches'
          WRITE(4,*)  ' There are periodic patches'
          WRITE(13,*) ' There are periodic patches'
          WRITE(45,*) 
      ENDIF

      WRITE(45,*) '  FINAL FORM OF ICON'
      WRITE(45,*)
      CALL PRIICO(ICON,NBCS)

      RETURN
      END SUBROUTINE INP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INPGRI

      USE MPI

      USE CHARACTERS

      USE CONSTANTS,   ONLY : EPS

      USE INTEGERS,    ONLY : MAXB,NMOV,IPRO,MGM,MAXFS,IREPEA,NREPEA,
     &                        NBCS

      USE GM

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,MOVPO,IMAXM,JMAXM,
     &    KMAXM,IMINM,JMINM,KMINM,IT,IL,IK,XCO,YCO,ZCO,HFLUX,BLANK2,
     &    MGRID,IG,IC,NPATCH,ICON,IHF,VOL,MOV,NPROCE,JF,XC,YC,ZC,
     &    XXI,YXI,XETA,YETA,NTOT,ICONH,ICOGH,IHULL,JHULL,XHULL,YHULL,
     &    ZHULL,LHULL,MMHUL,MHULL,IBOTGR,JBOTGR,KBOTGR,ITOPGR,JTOPGR,
     &    KTOPGR,XI,ETA,ZETA,ZZTOP,ICNH,IBTGR,ITPGR,JBTGR,JTPGR,BLANK,
     &    XSP,YSP,ZSP,WH,XGRI,YGRI,ZGRI,XCOL,PDIFF,IGRID,XORI,YORI,ZORI,
     &    XCP,YCP,ZCP,WAVEH,OMEGA,OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,
     &    CENAZ,UROTCP,VROTCP,WROTCP


      USE NS3CO           

      USE FLIGHT, ONLY : XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,PSIRI,THETARI,PHIRI,DRAUGHTI,DRAUGHT,
     &     XCGIS,YCGIS,ZCGIS,ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI

      USE BLADE_VARIABLES, ONLY : TDBL,DTB

      IMPLICIT NONE

      INTEGER :: ERRORCODE,IERR,N,M,IG1,IC1,I1,I2,IB,NMOVO,J1,J2,K,I,L,
     &           ILE4,NBG,NP,ITYPE,IFACE,NTP2,NM,
     &           IIMIN,JJMIN,KKMIN,IIMAX,JJMAX,KKMAX,II,JJ,KK,JSTR,KSTR,
     &           NTYPE,JH1,JG1,NLVL,IG2,NGL,ISTRID,JSTRID,
     &           ILE,KA,J,NTOT1,KU1,KU2,KU3,KU4,KU5,KU6,IS,IGR,
     &           IP,iiirc,KSTRID,III,IF1


C ... THIS SUBROUTINE MAKES GRID AT THE BEGINNING OF CALCULATION

C **********************************************************************

C ... EXTRAPOLATE AND CONNECT THE CORNER POINT COORDINATES

      IS = 1
                    
      DO N = 1,NBLOCK

         IG1  = IG(1,N)
         IC1  = IC(1,N)

         CALL GRIDEX(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(1,N),JMAX(1,N),
     &         KMAX(1,N),ICON(IC1),NPATCH(N),XCOL(1,1,IS),IN,JN,KN)

         IS  = IS + NPATCH(N) 

      ENDDO        

      CALL MIRC2(XCO,YCO,ZCO,1,1,NBLOCK)
        
      CALL CONECX(XCO,1,1,NBLOCK,1)
      CALL CONECX(YCO,1,1,NBLOCK,1)
      CALL CONECX(ZCO,1,1,NBLOCK,1)

C ... Store the generated grid and initialize the computational grid
      
      IF (COORL .AND. IOLD <= 0) THEN
           
      DO N = 1,NBLOCK
      NTOT1= (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IG1  = IG(1,N)
      IC1 = IC(1,N) 
      M   = 1
      NBG = NPROCE(1+N,IPRO)

C ... Store the grid as it is at the beginning of simulation

      DO L = 1,NTOT1
         I       = L + IG1 - 1
         XORI(I) = XCO(I)
         YORI(I) = YCO(I)
         ZORI(I) = ZCO(I)
         XGRI(I) = XCO(I)
         YGRI(I) = YCO(I)
         ZGRI(I) = ZCO(I)
      ENDDO
        
C **********************************************************************
C ... XGRI (and also XCO) can be put to a new position here (e.g. sink)
C **********************************************************************

      IGR = IGRID(NBG,3)                               ! Particle number

      IF(IGRID(NBG,1) >= 21 .AND. 
     &     IGRID(NBG,1) <= 29) THEN ! Sink and trim
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Set draught of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        0.0,0.0,DRAUGHTI(IGR)-DRAUGHT(IGR),0.0,
     &        0.0,0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Set trim angle of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGIS(IGR),YCGIS(IGR),ZCGIS(IGR),
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),0.0,
     &        0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) == 10) THEN ! Propel in ship case
C ... Set location(s) and attitude(s) of the propel(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the propel(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 11 .AND. 
     &     IGRID(NBG,1) <= 19) THEN ! Flying object
C ... Set angle(s) of control surface(s) to (X,Y,Z) ORI-table
         CALL FLAPPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &        IN,JN,KN,NBG,TDBL,DTB,NTOT(1,N))
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 31 .AND. 
     &     IGRID(NBG,1) <= 39) THEN ! Helicopter rotor blades
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, PSIRI(IGR),PHIRI(IGR),
     &        THETARI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),PHIR(IGR),
     &        THETAR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 40 .AND. 
     &     IGRID(NBG,1) <= 48) THEN ! Actuator disk
C ... Set location(s) and attitude(s) of the actuator disk(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, 0.0,ROTB1I(IGR),
     &        ROTA1I(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Deformation of the the actuator disk(s) 
         IF(IGRID(NBG,2) == 15) THEN
         CALL CONEDISK(XGRI(IG1),YGRI(IG1),ZGRI(IG1),
     &     XGRI(IG1),YGRI(IG1),ZGRI(IG1),NTOT(1,N),
     &     0.0,0.0,0.0,0.0,1.0,0.0,
     &     0.0,0.0,0.0,0.0,1.0,0.0,
     &     CONEA(IGR)-CONEAI(IGR))
         ENDIF
C ... Update location(s) and rotation angle(s) of the actuator disk(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),0.0,ROTB1(IGR),
     &        ROTA1(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) == 49) THEN ! Actuator disk in ship case
C ... Set location(s) and attitude(s) of the actuator disk(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, 0.0,ROTA1I(IGR),
     &        ROTB1I(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Deformation of the the actuator disk(s) 
c         IF(IGRID(NBG,2) == 15) THEN
c            CALL CONEDISK(XGRI(IG1),YGRI(IG1),ZGRI(IG1),
c     &           XGRI(IG1),YGRI(IG1),ZGRI(IG1),NTOT(1,N),
c     &           0,0,0,0,1.0,0,0,0,0,0,1.0,0,CONEAI(IGR)-CONEA(IGR))
c         ENDIF
C ... Update location(s) and rotation angle(s) of the actuator disk(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),0.0,ROTA1(IGR),
     &        ROTB1(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)
      ENDIF ! IF(IGRID(NBG,1) >= 21 .AND ...

      ENDDO ! NBLOCK
      ENDIF ! COORL .AND. IOLD <= 0
      
C ... ROTATIONAL AXIS AND THE PLACE AND VELOCITY (Global and obsolate)

      OMEX = 1.                 ! rotational axis
      OMEY = 0.
      OMEZ = 0.
      CENX = 0.                 ! center of the rotation
      CENY = 0.
      CENZ = 0.
c      GVEX = 0.                 ! grid velocity
c      GVEY = 0.
c      GVEZ = 0.

      NLVL = 2**(LEVEL-1)

      IF(ABS(OMEX**2 + OMEY**2 + OMEZ**2 - 1.) >= 1.E-4) THEN
         WRITE(*,*) 'OMEX, OMEY and OMEZ is not lenght of unity'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

C*******************CHECK****** For what purpose?
       IF (IOLD > 0.AND.IFSBC == 11) THEN

       OPEN(UNIT=997,FILE='GRID.BIN'//PRN,FORM='UNFORMATTED')

       DO N=1,NBLOCK
         IG1=IG(1,N)-1
         READ(997)(XCO(IG1+K),YCO(IG1+K),ZCO(IG1+K),K=1,NTOT(1,N))
       ENDDO
       CLOSE (997)

       ENDIF
C*******************CHECK****** For what purpose was this?

C ... Calculate the cell centerpoint coordinates and velocities

      DO N = 1,NBLOCK
      IG1  = IG(1,N)
      NGL  = NPROCE(1+N,IPRO)

      CALL CELLCP(XCO(IG1),YCO(IG1),ZCO(IG1),XC(IG1),YC(IG1),ZC(IG1),
     +            OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +            CENAX(NGL),CENAY(NGL),CENAZ(NGL),UROTCP(IG1),
     +            VROTCP(IG1),WROTCP(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +            IN,JN,KN,PDIFF(IG1),FRESUL,FRSDEN)
      ENDDO

C ... CALCULATE THE GEOMETRIC GRID PROPERTIES
C **********************************************************************

      IF(IPRINT-1 == ICYCLE) WRITE(46,*)
      IF(IPRINT-1 == ICYCLE) WRITE(46,*) ' CALL FROM INPGRI'
      WRITE(45,*)
      WRITE(45,*) 'GRIVAL CALL FROM INPGRI'
c      WRITE(46,*)
c      WRITE(46,*) 'GRIVAL CALL FROM INPGRI'
c      WRITE(45,*) '-----------------------'
      
      CALL GRIVAL(0,.FALSE.,'INPGRI') 

        DO I=1,NBLOCK

         IC1  = IC(1,I)
         IG1  = IG(1,I)
         NGL  = NPROCE(1+I,IPRO)
c            write(6,*) ngl,i,IPRO
c      IF(NGL == 69) THEN
      ISTRID   = IMAX(1,I) + 2*IN
      JSTRID   = JMAX(1,I) + 2*JN
      KSTRID   = ISTRID*JSTRID

      DO K = 1,KMAX(1,I)
      KA      = (KN+K-1)*KSTRID
      DO J = 1,JMAX(1,I)
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO  III = 1,IMAX(1,I)
         KK      = JJ + III + IG1
         BLANK2(KK) = BLANK(KK)
c         WRITE(670,*) iii,j,k,blank(kk)
      ENDDO ; ENDDO ; ENDDO
c      ENDIF
       enddo

      IREPEA(3) = IREPEA(3) + 1

C **********************************************************************

C ... Initialize XHULL with surface type 9

      IF(IFSBC /= 3) THEN

         XSP(1:MAXFS) = XCO(1:MAXFS)
         YSP(1:MAXFS) = YCO(1:MAXFS)
         ZSP(1:MAXFS) = ZCO(1:MAXFS)
        
         DO I=1,NBLOCK

         IC1  = IC(1,I)
         IG1  = IG(1,I)
         NGL  = NPROCE(1+I,IPRO)

         IF(IGRID(NGL,2) == 9) THEN ! Original FINFLOSHIP treatment
                                      
         DO L = 1,NPATCH(I)
           IF(ICON(IC1+(L-1)*IC9) == 13) THEN
           JH1  = ICONH(I)
           IF(JH1 /= 0 .AND. JH1 /= NBLOCG+1)THEN 
           JG1 = ICOGH(JH1)
           CALL INIHUL(                          !IHULL(JG1),JHULL(JG1),
     &     XHULL(JG1),YHULL(JG1),ZHULL(JG1),
     &     LHULL(JH1),MMHUL(JH1),MHULL(JH1),XSP(IG1),YSP(IG1),ZSP(IG1),
     &     IMAX(1,I),JMAX(1,I),KMAX(1,I),IBOTGR(I),ITOPGR(I),JBOTGR(I),
     &     JTOPGR(I),KBOTGR(I),KTOPGR(I),XI(IG1),ETA(IG1),ZETA(IG1),
     &     ZZTOP(IG1),IGLOBL,NGLOBL,NLVL,IOLD,ICYCLE,ICMAX,GML,CZ,CMY,
     &     PRN)
           ENDIF
           GOTO 2102
           ENDIF
         ENDDO
 2102    CONTINUE

         ENDIF ! IGRID(NGL,2) == 9
         ENDDO

c         XCO(1:MAXFS) = XSP(1:MAXFS)
c         YCO(1:MAXFS) = YSP(1:MAXFS)
c         ZCO(1:MAXFS) = ZSP(1:MAXFS)
C        
      IF(IFSBC == 1) CALL FS_COVOL(ICON,XCO,YCO,ZCO,IMAX,JMAX,KMAX,M,
     & NPATCH,IN,JN,KN,XXI,YXI,XETA,YETA,XC,YC,ZC,IG,IC,JF,MGM,NBLOCK)

      ELSEIF(IFSBC == 3) THEN
    
      DO N   = 1,NBLOCK
      DO M   = 1,MGRID(N)
         IG1 = IG(M,N)
         IC1 = IC(M,N)
C ... CALCULATE SURFACE HEIGHT
      CALL SURHEI(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),IN,JN,KN,HFLUX,N,M,NPATCH(N),ICON(IC1),IHF(1,M),
     +     VOL(IG1),IB,MAXB)
      ENDDO
      ENDDO

      ENDIF ! IFSBC == 1

      WRITE(45,*)
      WRITE(45,*)'END OF INPUT TREATMENT' 
      WRITE(45,*)'*****************************************************'
      WRITE(45,*)

C ... DIRECTIONS OF MOVIE FILES (all mirror walls and solids)

C ... MOV(N) = 1 same patches as OUTPUT directions
C ... MOV(N) = 2 mirror, fresurface and solid + 1 patches
C ... MOV(N) = 3 whole block
C ... if more is specified rememger to change sub. INPUT in NMOV calculation
      
      NMOVO    = NMOV
      NMOV     = 0
      MOVPO(1) = 1

      DO 1960 N = 1,NBLOCK
         IF(MOV(N) == 1) THEN ! OUTPUT directive
            NMOV = NMOV + 1
            MOVPO(NMOV) = N
            IMAXM(NMOV) = IMAX(1,N) + 1
            JMAXM(NMOV) = JMAX(1,N) + 1
            KMAXM(NMOV) = KMAX(1,N) + 1
            IMINM(NMOV) = 1
            JMINM(NMOV) = 1
            KMINM(NMOV) = 1
            IF(IL(1,N) == 1) IMAXM(NMOV) = IK(1,N)
            IF(IL(1,N) == 2) JMAXM(NMOV) = IK(1,N)
            IF(IL(1,N) == 3) KMAXM(NMOV) = IK(1,N)
            IF(IL(1,N) == 1) IMINM(NMOV) = IK(1,N)
            IF(IL(1,N) == 2) JMINM(NMOV) = IK(1,N)
            IF(IL(1,N) == 3) KMINM(NMOV) = IK(1,N)
         ELSEIF(MOV(N) == 2) THEN ! mirrors and solid+1
            DO NTYPE = 6,10  ! arrange order so that mirrors first, then solids
               NTP2  = NTYPE
               IF(NTYPE == 6) NTP2 = 4
               IF(NTYPE == 7) NTP2 = 13
               IC1     = IC(1,N) - 1
                  DO NP = 1,NPATCH(N)
                     ITYPE = ICON(IC1+(NP-1)*IC9 + 1)
                     IF(ITYPE == NTP2) THEN
                        IFACE = ICON(IC1+(NP-1)*IC9 + 3)
                        I1    = ICON(IC1+(NP-1)*IC9 + 4)
                        I2    = ICON(IC1+(NP-1)*IC9 + 5)
                        J1    = ICON(IC1+(NP-1)*IC9 + 6)
                        J2    = ICON(IC1+(NP-1)*IC9 + 7)
                        NMOV  = NMOV + 1
                        MOVPO(NMOV) = N
                        IMAXM(NMOV) = MIN0(IMAX(1,N),I2) + 1
                        JMAXM(NMOV) = MIN0(JMAX(1,N),I2) + 1
                        KMAXM(NMOV) = MIN0(KMAX(1,N),J2) + 1
                        IMINM(NMOV) = MAX0(1,I1)
                        JMINM(NMOV) = MAX0(1,I1)
                        KMINM(NMOV) = MAX0(1,J1)
                        ILE4 = 1
                        IF(ITYPE >= 8 .AND. ITYPE <= 10) ILE4 = 2 ! BCP

                        IF(IFACE == 1) IMAXM(NMOV) = ILE4
                        IF(IFACE == 2) JMAXM(NMOV) = ILE4
                        IF(IFACE == 3) KMAXM(NMOV) = ILE4

                        IF(IFACE == 4) IMINM(NMOV) = IMAX(1,N) + 2-ILE4
                        IF(IFACE == 5) JMINM(NMOV) = JMAX(1,N) + 2-ILE4
                        IF(IFACE == 6) KMINM(NMOV) = KMAX(1,N) + 2-ILE4

                        IF(IFACE == 1) IMINM(NMOV) = IMAXM(NMOV)
                        IF(IFACE == 2) JMINM(NMOV) = JMAXM(NMOV)
                        IF(IFACE == 3) KMINM(NMOV) = KMAXM(NMOV)
      
                        IF(IFACE == 4) IMAXM(NMOV) = IMINM(NMOV)
                        IF(IFACE == 5) JMAXM(NMOV) = JMINM(NMOV)
                        IF(IFACE == 6) KMAXM(NMOV) = KMINM(NMOV)
                        IF(IFACE == 3 .OR. IFACE == 6) THEN
                           JMAXM(NMOV) = MIN0(JMAX(1,N),J2) + 1
                           JMINM(NMOV) = MAX0(1,J1)
                        ENDIF
                     ENDIF
                  ENDDO
            ENDDO               !DO NTYPE = 7,10
         ELSEIF(MOV(N) == 3) THEN ! whole block
            NMOV = NMOV + 1
            MOVPO(NMOV) = N
            IMAXM(NMOV) = IMAX(1,N) + 1
            JMAXM(NMOV) = JMAX(1,N) + 1
            KMAXM(NMOV) = KMAX(1,N) + 1
            IMINM(NMOV) = 1
            JMINM(NMOV) = 1
            KMINM(NMOV) = 1
         ELSEIF(MOV(N) /= 0) THEN
            NBG = NPROCE(1+N,IPRO)
            WRITE(*,*) 'Impossible MOV value in block',NBG,
     +           '. MOV=',MOV(N),N
            IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
         ENDIF
 1960 CONTINUE
      IF(NMOV /= NMOVO) THEN
         WRITE(*,*) 'Different number of allocated NMOV and defined???'
         WRITE(*,*) 'NMOVO,NMOV',NMOVO,NMOV
      ENDIF

      GOTO 1980
* 1970 WRITE(*,*) 'If you want to make movie files you have to make '//
*     +     'subdirectory MOVIE'
      WRITE(*,*) 'Exiting...'
      IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
      STOP
 1980 CONTINUE

      IF (FRESUL) THEN
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO)
         IF1  = JF(1,N)
         IC1  = IC(1,N)
         CALL FS_SETWAVEHVALUE(WAVEH,WHEIGHT,IHF(1,1),ICON(IC1),
     +        NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +        ICYCLE,NGL,XCP,YCP,ZCP,XBULB,YBULB,ABULB)
      ENDDO
      ENDIF

      RETURN
      END SUBROUTINE INPGRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRIORI      

      USE CONSTANTS,   ONLY : EPS

      USE INTEGERS,    ONLY : IPRO,IREPEA,NREPEA

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,XCO,YCO,ZCO,MGRID,IG,
     &    NPROCE,NTOT,XC,YC,ZC,XGRI,YGRI,ZGRI,XORI,YORI,ZORI,IGRID

      USE NS3CO           

      USE FLIGHT,     ONLY  : XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,PSIRI,THETARI,PHIRI,DRAUGHTI,DRAUGHT,
     &     XCGIS,YCGIS,ZCGIS,ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI

      USE BLADE_VARIABLES, ONLY : TDBL,DTB

      IMPLICIT NONE

      INTEGER :: N, IG1, IGR, NTOT1, NBG
           
      DO N  = 1,NBLOCK
      NTOT1 = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IG1   = IG(1,N)
      NBG   = NPROCE(1+N,IPRO)
    
C **********************************************************************
C ... XGRI (and XCO) can be put to a new position in restart (e.g. sink)
C **********************************************************************

      IGR = IGRID(NBG,3)  ! Particle number  

      IF(IGRID(NBG,1) <= 9) THEN                       ! Use the original 

         XGRI(IG1:IG1+NTOT1-1) = XORI(IG1:IG1+NTOT1-1)
         YGRI(IG1:IG1+NTOT1-1) = YORI(IG1:IG1+NTOT1-1)  ! pitaako alustaa
         ZGRI(IG1:IG1+NTOT1-1) = ZORI(IG1:IG1+NTOT1-1)  ! aina kun 
                                                        ! IGRID(NBG,1) <= 10
                                                        ! ROT:lla ainakin
                                                        ! =IGRID(NBG,1)=1

      ELSEIF(IGRID(NBG,1) >= 21 .AND. IGRID(NBG,1) <= 29) THEN ! Sink and trim
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Set draught of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        0.0,0.0,DRAUGHTI(IGR)-DRAUGHT(IGR),0.0,
     &        0.0,0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Set trim angle of the moving grid(s)
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGIS(IGR),YCGIS(IGR),ZCGIS(IGR),
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),0.0,
     &        0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) == 10 .AND. .NOT. TIMEL) THEN ! Propel in ship case
C ... Set location(s) and attitude(s) of the propel(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the propel(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 11 .AND. IGRID(NBG,1) <= 19 ! Flying object,
     &        .AND. .NOT. TIMEL) THEN  ! time-accurate newpos is done later 
C ... Set angle(s) of control surface(s) to (X,Y,Z) ORI-table
         CALL FLAPPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &        IN,JN,KN,NBG,TDBL,DTB,NTOT(1,N))
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 31 .AND. 
     &        IGRID(NBG,1) <= 39 .AND. .NOT. TIMEL) THEN ! Helicopter blade
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),PHIRI(IGR),
     &        THETARI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),PHIR(IGR),
     &        THETAR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 40 .AND. 
     &     IGRID(NBG,1) <= 48) THEN ! Actuator disk
C ... Set location(s) and attitude(s) of the actuator disk(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, 0.0,ROTB1I(IGR),
     &        ROTA1I(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Deformation of the the actuator disk(s) 
         IF(IGRID(NBG,2) == 15) THEN
         CALL CONEDISK(XGRI(IG1),YGRI(IG1),ZGRI(IG1),
     &     XGRI(IG1),YGRI(IG1),ZGRI(IG1),NTOT(1,N),
     &     0.0,0.0,0.0,0.0,1.0,0.0,
     &     0.0,0.0,0.0,0.0,1.0,0.0,
     &     CONEA(IGR)-CONEAI(IGR))
         ENDIF
C ... Update location(s) and rotation angle(s) of the actuator disk(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),0.0,ROTB1(IGR),
     &        ROTA1(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) == 49) THEN ! Actuator disk in ship case
C ... Set location(s) and attitude(s) of the actuator disk(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, 0.0,ROTA1I(IGR),
     &        ROTB1I(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Deformation of the the actuator disk(s) 
c         IF(IGRID(NBG,2) == 15) THEN
c            CALL CONEDISK(XGRI(IG1),YGRI(IG1),ZGRI(IG1),
c     &           XGRI(IG1),YGRI(IG1),ZGRI(IG1),NTOT(1,N),
c     &           0,0,0,0,1.0,0,0,0,0,0,1.0,0,CONEAI(IGR)-CONEA(IGR))
c         ENDIF
C ... Update location(s) and rotation angle(s) of the actuator disk(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),0.0,ROTA1(IGR),
     &        ROTB1(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ENDIF ! IF(IGRID(NBG,1) <= 9) THEN
      ENDDO ! NBLOCK

      END SUBROUTINE GRIORI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INIT         

      USE MPI

      USE CHARACTERS

      USE CONSTANTS,   ONLY : EPS,PII

      USE INTEGERS,    ONLY : MAXB,MGM,IPRO,MAXSB,MAXEB,MBPRO

      USE TYPE_ARRAYS

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,IG,NTOT,INITC,RO,P,
     &    PDIFF,RM,RN,RW,E,U,V,W,TEMP,VIS,RK,REPS,VIST,EPS2,FI,
     &    RKSI,ROFOR,RMFOR,RNFOR,RWFOR,EFOR,PDFOR,RKFOR,REFOR,FIFOR,
     &    XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,ZZZ,TTS,A1XA,A1YA,A1ZA,A2XA,
     &    A2YA,A2ZA,A3XA,A3YA,A3ZA,F1R,F1RM,BIJ,XC,YC,ZC,JLOC,JTRA,
     &    APP,XCO,YCO,ZCO,ROLE2,RMLE2,RNLE2,RWLE2,ELE2,ROLE3,RMLE3,
     &    RNLE3,RWLE3,ELE3,RKLE2,EPSLE2,RKLE3,EPSLE3,FILE2,FILE3,
     &    IGRID,NPROCE,OMEGA,UROT,VROT,WROT,IT,IL,IK,A1,A2,A3,
     &    IQ,IR,JF,IC,PRO,VAR,BLKS,JSTATE,F1RN,RNUT,TRM,MPFOR,
     &    ROAV1,RMAV1,RNAV1,RWAV1,RKAV1,EPSAV1,EAV1,PAV1,TAV1,
     &    ROAV2,RMAV2,RNAV2,RWAV2,RKAV2,EPSAV2,EAV2,RUVAV,RUWAV,RVWAV,
     &    IUPPT,CENAX,CENAY,CENAZ,PRC

      USE NS3CO            

      USE INTEGERS, ONLY : MAXTI

      USE FLIGHT,      ONLY : NGRIFL,OSKU,XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR

      USE BLADE_VARIABLES, ONLY : TDBL
          
      IMPLICIT NONE

      REAL :: C1,C2,C3,C21,CMU,CTA,P11GK,P11GE,AA1,
     &        ETA0,U10,V10,W10,UU10,VV10,WW10,FRSMUD,APU,DELTAP,
     &        UIROT,VIROT,WIROT,GTOT,RKHI,FV1,FRSG,FRSRET,GM1,H
 
      REAL, ALLOCATABLE, DIMENSION(:) :: FRSDEN_L, FRSPRE_L, FRSSIE_L, 
     &                                   FRSTEM_L, FRSVIS_L

      REAL :: DFRSTEM

      INTEGER :: N,ISTRID,JSTRID,KSTRID,IG1,I,J,K,L,KA,KK,NN,NS,NGL,
     &           JJ,NTOT1,LEVEP1,ISTFLX,INITCC,ISTFL2,IERR,IG2,IR2,IQ2,
     &           IF2,IC2,M,IGN,IGT,ICHARP,IPHASE,IGR,III

      CHARACTER(LEN=3) :: IGRNUM

      CHARACTER(LEN=1) :: CLEVR, CLEVW

      LOGICAL :: THERE

      REAL :: RADCG, RADGP

      REAL, EXTERNAL :: ATMOSRHO, ATMOSP, ATMOST

      
C ... INITIALIZATION FOR NOZZLE FLOW ALGORITHM
C ... OUTPUT : BOUNDARY CONDITIONS AND INITIAL CONDITIONS

      T       = 0.
      TDBL    = 0.
      INERR   = 1
      ICYCLE  = 0
      ICYTOT  = 0
      ITIMES  = 0
      IPRINT  = 1
      IF(.NOT.TIMEL) IBOUT  = 3
      DROMAX  = 1.


      CALL TURBCO(ITURB,JRDIS,JRPRE,C1,C2,C3,C21,CMU,CTA,P11GK,
     +     P11GE,AA1,ETA0)

      
C ... BLOCK LOOP 1 BEGINS

      DO 2100 N = 1,NBLOCK

      ISTRID  = IMAX(1,N) + 2*IN
      JSTRID  = JMAX(1,N) + 2*JN
      KSTRID  = ISTRID*JSTRID
      IG1     = IG(1,N) - 1
      NGL     = NPROCE(1+N,IPRO)

      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSVEL  = BLKS(NGL)%FRSVEL
      FRSVIS  = BLKS(NGL)%FRSVIS
      FRSTEM  = BLKS(NGL)%FRSTEM
      TEMINI  = BLKS(NGL)%TEMINI
      FRSSIE  = BLKS(NGL)%FRSSIE
      SIEINI  = BLKS(NGL)%SIEINI
      CHLREF  = BLKS(NGL)%CHLREF
      RMUINI  = BLKS(NGL)%RMUINI
      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      FRSMUT  = BLKS(NGL)%FRSMUT
      SOLTEM  = BLKS(NGL)%SOLTEM
      DFRSTEM = BLKS(NGL)%DFRSTEM
      FRSG    = BLKS(NGL)%FRSG
      FRSRET  = BLKS(NGL)%FRSRET
      REFVEL  = BLKS(NGL)%REFVEL
      GVEX    = BLKS(NGL)%GVEX
      GVEY    = BLKS(NGL)%GVEY
      GVEZ    = BLKS(NGL)%GVEZ

      GM1 = GAMMA - 1.0

      ALLOCATE(FRSDEN_L(NTOT(1,N)),FRSPRE_L(NTOT(1,N)),
     &   FRSSIE_L(NTOT(1,N)),FRSTEM_L(NTOT(1,N)),FRSVIS_L(NTOT(1,N)))

      GTOT = SQRT(GX**2 + GY**2 + GZ**2)

      DO L=1,NTOT(1,N)

         IF(ALTITUDE >= GROUND) THEN  ! Standard atmosphere

            IF(INITC(N) == 2) THEN  ! Symmetric pull-up

               RADCG = SQRT((CENAX(N)-XMOM)**2
     &                    + (CENAY(N)-YMOM)**2)
               RADGP = SQRT((CENAX(N)-XC(L+IG1))**2
     &                    + (CENAY(N)-YC(L+IG1))**2)

               H = -(GX*XMOM + GY*YMOM + GZ*ZMOM 
     &             + GTOT*GROUND)/GTOT + (RADCG-RADGP)

            ELSE

*               H = -(GX*XC(L+IG1) + GY*YC(L+IG1) + GZ*ZC(L+IG1) 
*     &             + GTOT*GROUND)/GTOT

               H = -(GX*(XC(L+IG1)-XMOM) +
     &               GY*(YC(L+IG1)-YMOM) +
     &               GZ*(ZC(L+IG1)-ZMOM) + 
     &               GTOT*GROUND)/GTOT

            ENDIF
            
            FRSDEN_L(L) = ATMOSRHO(H)                           
            FRSPRE_L(L) = ATMOSP(H)
            FRSTEM_L(L) = ATMOST(H)
            FRSSIE_L(L) = -E0REF + RGAS/GM1*(FRSTEM-T0REF)
            FRSVIS_L(L) = VISU0*FRSTEM**EXPSU/(FRSTEM+TSU0)

         ELSE

            FRSDEN_L(L) = FRSDEN                           
            FRSPRE_L(L) = FRSPRE
            FRSTEM_L(L) = FRSTEM
            FRSSIE_L(L) = FRSSIE
            FRSVIS_L(L) = FRSVIS 

         ENDIF

      ENDDO

C ... INITIAL (AND BOUNDARY) CONDITIONS

      U10 = FRSVEL*COS(ALPHA)*COS(BETA)
      V10 = FRSVEL*SIN(ALPHA)
      W10 = FRSVEL*COS(ALPHA)*SIN(BETA)

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'MULTI'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

      IF(INITC(N) /= 2) THEN
         UU10 = U10
         VV10 = V10
         WW10 = W10
      ELSE IF(INITC(N) == 2) THEN  ! A pull-up case
         UU10 = 0. 
         VV10 = 0.
         WW10 = 0.
      ENDIF

      DO 2000 L = 1,NTOT(1,N)

         I = IG1 + L

         FRSDEN = FRSDEN_L(L)                          
         FRSPRE = FRSPRE_L(L)
         FRSTEM = FRSTEM_L(L)
         FRSSIE = FRSSIE_L(L)
         FRSVIS = FRSVIS_L(L)

         RO(I)    = FRSDEN
         P(I)     = FRSPRE
*         PDIFF(I) = 0.
         PDIFF(I) = FRSPRE - BLKS(NGL)%FRSPRE
         RM(I)    = FRSDEN*UU10
         RN(I)    = FRSDEN*VV10
         RW(I)    = FRSDEN*WW10
         U(I)     = UU10
         V(I)     = VV10
         W(I)     = WW10
         E(I)     = FRSDEN*(FRSSIE +.5*(UU10**2 + VV10**2 + WW10**2))
         TEMP(I)  = FRSTEM
         VIS(I)   = FRSVIS  ! First estimate for viscosity

2000  CONTINUE

      IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
            DO L = 1,NTOT(1,N)
               I    = IG1 + L
               PRO(I)%TEMP(IPHASE)  = FRSTEM
               PRO(I)%DTEMP(IPHASE) = DFRSTEM
               VAR(I)%ALFA(IPHASE)  = BLKS(NGL)%FRSALFA(IPHASE)
               VAR(I)%X(IPHASE)     = BLKS(NGL)%FRSX(IPHASE)
               PRO(I)%TSAT          = FRSTEM !  Temporary initialization
            ENDDO
         ENDDO 
      ENDIF
        
      IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN ! Two-fluid model

         DO IPHASE = 1,BLKS(NGL)%NPHASE
            DO L = 1,NTOT(1,N)
               I = IG1 + L
               VAR(I)%U(IPHASE) = UU10            ! Currently these
               VAR(I)%V(IPHASE) = VV10
               VAR(I)%W(IPHASE) = WW10
            ENDDO
         ENDDO 

      ENDIF

   
      IF(TRANSL) THEN  ! Intermittency variables
         DO L = 1,NTOT(1,N)
            I = IG1 + L
            TRM(I)%G   = FRSG
            TRM(I)%RET = FRSRET
         ENDDO 
      ENDIF

      DO K = 1,KMAX(1,N)
         KA      = (KN+K-1)*KSTRID
         DO J = 1,JMAX(1,N)
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO I = 1,IMAX(1,N)
            KK        = JJ + I + IG1

            FRSDEN = FRSDEN_L(JJ+I)                          
            FRSPRE = FRSPRE_L(JJ+I)
            FRSTEM = FRSTEM_L(JJ+I)
            FRSSIE = FRSSIE_L(JJ+I)
            FRSVIS = FRSVIS_L(JJ+I)

            TEMP(KK) = TEMINI
            IF(ALTITUDE >= GROUND) TEMP(KK) = FRSTEM

            IF(MULPHL) THEN
               DO IPHASE = 1,BLKS(NGL)%NPHASE
               PRO(KK)%TEMP(IPHASE)  = TEMINI
               PRO(KK)%DTEMP(IPHASE) = TEMINI
               ENDDO
            ENDIF

            E(KK) = FRSDEN*(SIEINI+.5*(UU10**2 + VV10**2 + WW10**2)) 

            ENDDO
         ENDDO
      ENDDO

      FRSDEN = BLKS(NGL)%FRSDEN 
      FRSPRE = BLKS(NGL)%FRSPRE 
      FRSTEM = BLKS(NGL)%FRSTEM 
      FRSSIE = BLKS(NGL)%FRSSIE 
      FRSVIS = BLKS(NGL)%FRSVIS 

      IF(ITURB >= 3) THEN        ! K-EPSILON
      FRSMUD = 1.0E-3            ! Default value for free-stream mut/mu

      IF ((TURBLE < EPS).OR.(RMUINI/FRSVIS < EPS)) THEN
         IF(IDIS == 2) THEN      ! k-omega model
            FRSMUT = FRSMUD*FRSVIS
            RMUINI = FRSMUT
            FRSEPS = FRSDEN*REFVEL/CHLREF ! Changed on May 25th 2009
            FRSRK  = FRSMUT*FRSEPS/FRSDEN
         ELSEIF(IDIS == 1) THEN  ! k-epsilon model
            FRSMUT = FRSMUD*FRSVIS
            RMUINI = FRSMUT
            FRSRK  = RMUINI*1.0*REFVEL/CHLREF
            FRSEPS = CMU*FRSRK**2/RMUINI
         ENDIF
      ELSE
         IF(IDIS == 1) THEN      ! k-epsilon model
            FRSRK  = 1.5*FRSDEN*(TURBLE*REFVEL)**2 ! NEW 17.11.95, but
            FRSEPS = CMU*FRSRK**2/RMUINI           ! changed 25.5.2009
         ELSEIF(IDIS == 2) THEN  ! k-omega model
            FRSRK  = 1.5*FRSDEN*(TURBLE*REFVEL)**2
            FRSEPS = FRSDEN*FRSRK/RMUINI
      ENDIF
         IF(ITURB == 5) THEN     ! RNG K-EPSILON
            APU    = (1+SQRT(1.+RMUINI))**2
            FRSEPS = CMU*APU*FRSRK**2/RMUINI
      ENDIF
      IF(FRSRK < RKLIM) WRITE(4,*) 'Warning, RKIN smaller than RKLIM'
      IF(FRSEPS < EPSLIM) WRITE(4,*)'Warning, EPSIN smaller th. EPSLIM'
      ENDIF

C ... Added 20.3.2003
      FRSRK  = MAX(FRSRK,RKLIM)
      FRSEPS = MAX(FRSEPS,EPSLIM)

C ... Free-stream values are specified here

      IF(ITURB /= 9) THEN ! Two-equation models

      DO 2010 L = 1,NTOT(1,N)

         I = IG1 + L

         FRSDEN = FRSDEN_L(L)                          
         FRSPRE = FRSPRE_L(L)
         FRSTEM = FRSTEM_L(L)
         FRSSIE = FRSSIE_L(L)
         FRSVIS = FRSVIS_L(L)

         RK(I)   = RKLIM   ! FRSRK
         REPS(I) = EPSLIM  ! FRSEPS
         VIST(I) = FRSMUT  ! RMUINI
         EPS2(I) = 1. + FRSMUT/FRSVIS  ! RMUINI/FRSVIS
         E(I)    = E(I) + RK(I)

 2010 CONTINUE

      ELSE IF(ITURB == 9) THEN  ! Spalart-Allmaras

         RKHI = 1.

         DO III  = 1,20
            RKHI = SQRT(SQRT(FRSMUT/FRSVIS*(RKHI**3+7.1**3)))
         ENDDO ! Iteration also in main.f (mersu)

         DO L = 1,NTOT(1,N)

         I = IG1 + L

         FRSDEN = FRSDEN_L(L)                          
         FRSPRE = FRSPRE_L(L)
         FRSTEM = FRSTEM_L(L)
         FRSSIE = FRSSIE_L(L)
         FRSVIS = FRSVIS_L(L)

         REPS(I) = RKHI*FRSVIS          ! Value for RNUT*RO
         RNUT(I) = RKHI*FRSVIS/RO(I)
         VIST(I) = FRSMUT ! RMUINI
         EPS2(I) = 1. + FRSMUT/FRSVIS   ! RMUINI/FRSVIS
         RK(I)   = 0.
         ENDDO

      ENDIF

C ... Initial conditions are given here

      IF(TRANSL) THEN ! Intermittency variables 

         DO K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO I = 1,IMAX(1,N)
                  KK = JJ + I + IG1
                  TRM(KK)%G   = FRSG
                  TRM(KK)%RET = FRSRET
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF(ITURB /= 8 .OR. ITURB /= 9) THEN ! Two-equation models

         DO 2020 K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO 2020 J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO 2020 I = 1,IMAX(1,N)
                  KK      = JJ + I + IG1

                  FRSDEN = FRSDEN_L(JJ+I)                          
                  FRSPRE = FRSPRE_L(JJ+I)
                  FRSTEM = FRSTEM_L(JJ+I)
                  FRSSIE = FRSSIE_L(JJ+I)
                  FRSVIS = FRSVIS_L(JJ+I)

                  RK(KK)   = FRSRK
                  REPS(KK) = FRSEPS
                  VIST(KK) = RMUINI
                  EPS2(KK) = 1. + RMUINI/FRSVIS
                  E(KK)    = E(KK) + RK(KK)

2020  CONTINUE

      ELSE IF(ITURB == 9) THEN  ! Spalart-Allmaras

         RKHI = 1.

         DO III  = 1,20
            RKHI = SQRT(SQRT(RMUINI/FRSVIS*(RKHI**3+7.1**3)))
         ENDDO

         DO K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO I = 1,IMAX(1,N)
                  KK       = JJ + I + IG1

                  FRSDEN = FRSDEN_L(JJ+I)                          
                  FRSPRE = FRSPRE_L(JJ+I)
                  FRSTEM = FRSTEM_L(JJ+I)
                  FRSSIE = FRSSIE_L(JJ+I)
                  FRSVIS = FRSVIS_L(JJ+I)

                  REPS(KK) = RKHI*FRSVIS
                  RNUT(KK) = REPS(KK)/RO(KK) ! Value for RNUT*RO
                  VIST(KK) = RMUINI
                  EPS2(KK) = 1. + RMUINI/FRSVIS
                  RK(KK)   = 0.

               ENDDO
            ENDDO
         ENDDO

      ENDIF


      ENDIF    ! K-EPSILON

      IF(IPRO == 1) THEN
         WRITE(4,*)
         WRITE(4,'(2A,I4)') ' INLET (OR FREE STREAM) VALUES AFTER INIT',
     +               ' FOR BLOCK ',NGL
         WRITE(4,*) '**************************************************'
         WRITE(4,*) ' FRSDEN =',REAL(FRSDEN,4),' kg/m3'
         WRITE(4,*) ' FRSPRE =',REAL(FRSPRE,4),' N/m2'
         WRITE(4,*) ' FRSVEL =',REAL(FRSVEL,4),' m/s'
         WRITE(4,*) ' REFVEL =',REAL(REFVEL,4),' m/s'
         WRITE(4,*) ' GVEX, GVEY, GVEZ =',
     +                REAL(GVEX,4),REAL(GVEY,4),REAL(GVEZ,4),' m/s'
         WRITE(4,*) ' FRSTEM =',REAL(FRSTEM,4),' K'
         WRITE(4,*) ' TOTAL ENERGY(no turbulence)', 
     +        REAL(FRSDEN*(SIEINI +.5*(UU10**2 + VV10**2 + WW10**2)),4),
     +        ' J/kg'
         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            WRITE(4,*) ' RMUINI =',REAL(RMUINI,4), ' kg/ms'
            WRITE(4,*) ' RHO K  = ',REAL(FRSRK,4), ' J/m3'
            IF(ITURB /= 6) THEN
            WRITE(4,*) ' RHO EPSILON = ',REAL(FRSEPS,4),' J/m3s'
      ELSE
            WRITE(4,*) ' RHO OMEGA = ',REAL(FRSEPS,4),' kg/m3s'
            ENDIF
         ENDIF
      ENDIF

      IF(ITURB >= 20) THEN ! RSTM
      DO 2030 L = 1,NTOT(1,N)
         I      = IG1 + L
         FI(I,1) = 2./3.*FRSRK
         FI(I,4) = 2./3.*FRSRK
         FI(I,6) = 2./3.*FRSRK
         FI(I,2) = 0.
         FI(I,3) = 0.
         FI(I,5) = 0.
 2030 CONTINUE
      ENDIF    ! RSTM

      IF(INITC(N) == 0 .OR. INITC(N) == 5) THEN

C ... INITIAL CONDITION OF VELOCITY IS PUT TO ZERO

         DO 2060 K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO 2060 J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO 2060 I = 1,IMAX(1,N)
                  KK      = JJ + I + IG1

                  FRSDEN = FRSDEN_L(JJ+I)                          
                  FRSPRE = FRSPRE_L(JJ+I)
                  FRSTEM = FRSTEM_L(JJ+I)
                  FRSSIE = FRSSIE_L(JJ+I)
                  FRSVIS = FRSVIS_L(JJ+I)
                  
                  RM(KK)   = 0.
                  RN(KK)   = 0.
                  RW(KK)   = 0.
                  U(KK)    = 0.
                  V(KK)    = 0.
                  W(KK)    = 0.
                  TEMP(KK) = TEMINI
                  E(KK)    = FRSDEN*SIEINI

2060  CONTINUE

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN

         DO 2070 K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO 2070 J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO 2070 I = 1,IMAX(1,N)
                  KK      = JJ + I + IG1

                  FRSDEN = FRSDEN_L(JJ+I)                          
                  FRSPRE = FRSPRE_L(JJ+I)
                  FRSTEM = FRSTEM_L(JJ+I)
                  FRSSIE = FRSSIE_L(JJ+I)
                  FRSVIS = FRSVIS_L(JJ+I)

                  E(KK)  = FRSDEN*(SIEINI + RK(KK))

2070  CONTINUE

      ENDIF  ! ITURB >= 3
      ENDIF  ! INITC(N)

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'SOLID') THEN

      DO L = 1,NTOT(1,N)

         I = IG1 + L

         FRSDEN = FRSDEN_L(L)                          
         FRSPRE = FRSPRE_L(L)
         FRSTEM = FRSTEM_L(L)
         FRSSIE = FRSSIE_L(L)
         FRSVIS = FRSVIS_L(L)

         TEMP(I)  = SOLTEM
         P(I)     = FRSPRE
         PDIFF(I) = 0.
         U(I)     = 0.
         V(I)     = 0.
         W(I)     = 0.
         E(I)     = 0. 
         EPS2(I)  = 1.
         VIST(I)  = 0.

         IF(MULPHL) THEN 
            DO IPHASE = 1,NPHASES
               PRO(I)%TEMP(IPHASE) = SOLTEM
               VAR(I)%ALFA(IPHASE) = BLKS(NGL)%FRSALFA(IPHASE)
               VAR(I)%X(IPHASE)    = BLKS(NGL)%FRSX(IPHASE)
            ENDDO
         ENDIF

C     IF(TRANSL) THEN
C     Solid ehdot tanne.
C     ENDIF

      ENDDO

      CALL ROFPT(TEMP(IG1+1),P(IG1+1),RO(IG1+1),NTOT(1,N),JSTATE(NGL,1),
     +     RGAS,GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPT(NGL))
      
      ENDIF

      DEALLOCATE(FRSDEN_L,FRSPRE_L,FRSSIE_L,FRSTEM_L,FRSVIS_L)

2100  CONTINUE

C ... END OF BLOCK LOOP 1

C ... BLOCK LOOP 2 BEGINS

      DO 3100 N = 1,NBLOCK
      ISTRID   = IMAX(1,N) + 2*IN
      JSTRID   = JMAX(1,N) + 2*JN
      KSTRID   = ISTRID*JSTRID
      IG1     = IG(1,N) - 1

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID' .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'MULTI' .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

      DO 3000 K = 1,KMAX(1,N)
      KA      = (KN+K-1)*KSTRID
      DO 3000 J = 1,JMAX(1,N)
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 3000 I = 1,IMAX(1,N)
      KK      = JJ + I + IG1
      EPS2(KK) = 1.
      VIST(KK) = 0.
 3000 CONTINUE

      IF(GRAVIL) THEN ! HYDROSTATIC PRESSURE DIFFERENCE

      DO 2005 L = 1,NTOT(1,N)
      I        = IG1 + L
      DELTAP   = RO(I)*(GX*(XC(I) - GROUND) + GY*(YC(I) - GROUND)
     +          +       GZ*(ZC(I) - GROUND))
      P(I)     = P(I) + DELTAP
C...  PDIFF WITHOUT GRAVIATION LEVEL
      E(I)     = E(I)
2005  CONTINUE
      ENDIF   ! GRAVIL
      ENDIF   ! FLUID and MULTI
 3100 CONTINUE! END OF BLOCK LOOP 2

C ... If chimera is used

      DO I = 1,MAXB
         RKSI(I) = 0.
      ENDDO

      IF(NCHIM > 0) THEN

         DO I = 1,MAXB
            ROFOR(I) = RO(I)
            RMFOR(I) = 0.
            RNFOR(I) = 0.
            RWFOR(I) = 0.
            EFOR(I)  = E(I)
            PDFOR(I) = PDIFF(I)
         ENDDO

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            DO I = 1,MAXB
               RKFOR(I) = 0.
               REFOR(I) = 0.
            ENDDO
         ENDIF

         IF(TRANSL) THEN
            DO I = 1,MAXB
               TRM(I)%GFOR   = 0.
               TRM(I)%RETFOR = 0.
            ENDDO
         ENDIF

         IF(MULPHL) THEN
            DO IPHASE = 1,NPHASE
               DO I = 1,MAXB
                  MPFOR(I)%ALFAFOR(IPHASE)   = 0.
                  MPFOR(I)%XFOR(IPHASE)      = 0.
                  MPFOR(I)%DTEMPFOR(IPHASE)  = 0.
               ENDDO
            ENDDO
         ENDIF

         IF(NSCAL > 0) THEN
            DO NS = 1,NSCAL
            DO I = 1,MAXSB
               FIFOR(I,1) = 0.
            ENDDO
            ENDDO
      ENDIF

      ENDIF  ! NCHIM > 0

C ... GENERATE A NEW GRID AT THE BEGINNING OF THE CALCULATION

      IF (COORL) THEN
      
      DO 3300 N = 1,NBLOCK
      NTOT1   = NTOT(1,N)
      IG1     = IG(1,N)
            DO 3300 L  = 1,NTOT1
               I       = L + IG1 - 1

               XLE2(I) = XCO(I)
               YLE2(I) = YCO(I)
               ZLE2(I) = ZCO(I)

               XLE3(I) = XLE2(I)
               YLE3(I) = YLE2(I)
               ZLE3(I) = ZLE2(I)
 3300    CONTINUE

      ENDIF ! COORL

C ... INITIALIZATION FROM THE PREVIOUS GRID LEVEL. OPEN AN OLD COMPUT FILE

      IF(IOLD < 0) THEN

      LEVEP1 = LEVEL + 1
      CALL NUMCH1(CLEVR,LEVEP1)
      OPEN(16,FILE='COMPUT.'//CLEVR,STATUS='UNKNOWN',FORM='UNFORMATTED')
       
      IF(FRESUL) THEN
         INQUIRE(FILE='WH.DAT'//CLEVR,EXIST=THERE)
         IF(THERE) THEN         ! Open the free-surface file
           OPEN(21,FILE='WH.DAT'//CLEVR,STATUS='UNKNOWN',
     +     FORM='FORMATTED')
         ELSE 
           WRITE(*,*)  ' Warning: File WH.DAT was not found.'
           WRITE(45,*) ' Warning: File WH.DAT was not found.'
           WRITE(4 ,*) ' Warning: File WH.DAT was not found.'
         ENDIF
      ENDIF

      CALL INITV
      
      ENDIF

C ######################################################################
      IF (COORL .AND. IOLD <= 0)  CALL NEWROT
C ######################################################################


C ... SET INITIAL VALUES AT THE BOUNDARIES

      WRITE(3,*) 'menossa boundiin 1. kerran'

      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         IGN    = IG1 - 1 + NTOT(1,N)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(ITURB == 9) THEN
         CALL PRINYS(3,NUTC,RNUT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF

         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP    = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(1:NTOT(1,N)) = PRO(IG1:IGN)%TEMP(IPHASE)
         F1RN(1:NTOT(1,N)) = VAR(IG1:IGN)%X(IPHASE)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RN(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF

         CALL PRINYS(3,PC,   P(IG1),   IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,ROC,  RO(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,  RM(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDDO

C***********************************************************************
      CALL BOUNDI(1)
C***********************************************************************
          
      WRITE(3,*) 'Values after the first CALL BOUND'

      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         IGN    = IG1 - 1 + NTOT(1,N)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,UC,U(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(1:NTOT(1,N)) = PRO(IG1:IGN)%TEMP(IPHASE)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         F1RM(1:NTOT(1,N)) = VAR(IG1:IGN)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),F1RM(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF

         CALL PRINYS(3,PC,   P(IG1),   IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,ROC,  RO(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,  RM(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDDO

C ... SET INITIAL VALUES FOR INTERNAL FLOWS

      IF(IOLD >= 0) THEN

C **********************************************************************
C ... BLOCK lOOP 3 BEGINS
      
      DO 5000 NN = 1,NBLOCG

         ISTFLX = 0
         IGT    = IG(1,NBLOCK+1) ! Maximum array size

      DO 4900 N = 1,NBLOCK

      IG1    = IG(1,N)
      NGL    = NPROCE(1+N,IPRO)

      FRSDEN = BLKS(NGL)%FRSDEN
      FRSPRE = BLKS(NGL)%FRSPRE
      FRSVEL = BLKS(NGL)%FRSVEL
      FRSVIS = BLKS(NGL)%FRSVIS
      FRSTEM = BLKS(NGL)%FRSTEM
      TEMINI = BLKS(NGL)%TEMINI
      FRSSIE = BLKS(NGL)%FRSSIE
      SIEINI = BLKS(NGL)%SIEINI
      CHLREF = BLKS(NGL)%CHLREF
      RMUINI = BLKS(NGL)%RMUINI
      RKLIM  = BLKS(NGL)%RKLIM
      EPSLIM = BLKS(NGL)%EPSLIM
      FRSMUT = BLKS(NGL)%FRSMUT
      SOLTEM = BLKS(NGL)%SOLTEM

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'MULTI'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

         IF(INITC(N) > 10) THEN
         INITCC = INITC(N)
         IF(INITC(N) >= 20) INITCC = INITC(N) - 10
      CALL SETFLX(RO(IG1),U(IG1),V(IG1),W(IG1),E(IG1),XC(IG1),
     + YC(IG1),ZC(IG1),A1(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     + A2(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3(IG1),A3XA(IG1),
     + A3YA(IG1),A3ZA(IG1),INITCC,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + ZZZ)
      ISTFLX = 1
      ENDIF ! INITC(N) > 10

      ENDIF ! FLUID and MULTI
 4900 CONTINUE

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(ISTFLX,ISTFL2,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,IERR)
         ISTFLX = ISTFL2
      ENDIF

      IF(ISTFLX == 1) THEN
      CALL MIR(TEMP,U,V,W,E,RO,VIST,P,PDIFF,RK,REPS,TTS,FI,PRO,VAR,
     +     TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     +     1,NBLOCK,1.,0,PRC) ! huomaa ro
      CALL CYC(U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     +1,NBLOCK,F1R,F1RM,VAR)
c      IF(TWO_FLUIDL) THEN
c      VAR(1:MAXB)%U(IPHASE) = U(1:MAXB)
c      VAR(1:MAXB)%V(IPHASE) = V(1:MAXB)
c      VAR(1:MAXB)%W(IPHASE) = W(1:MAXB)
c      ENDIF
c      ENDDO
      IF(IFSBC == 3) THEN ! Boundary treatment without grid movement
        CALL FRE(1,1,NBLOCK,1.,0) ! huomaa ro
      ENDIF
      CALL CONEC(U,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(V,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(W,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)

      CALL CONEC(RO,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(E,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(EPS2,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(VIST,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(   P,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(PDIFF,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(TEMP, F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
         CALL CONEC(PRO(1:IGT)%TEMP(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,0)
         CALL CONEC(VAR(1:IGT)%ALFA(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,0)
         CALL CONEC(VAR(1:IGT)%X(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,0)
         PRO(1:IGT)%DTEMP(IPHASE) = PRO(1:IGT)%TEMP(IPHASE)
         IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
         CALL CONEC(VAR(1:IGT)%U(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(VAR(1:IGT)%V(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(VAR(1:IGT)%W(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         ENDIF
         ENDDO
      ENDIF

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL CONEC(RK,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(REPS,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      ENDIF

      IF(TRANSL) THEN           !Intermittency variables
           CALL CONEC(TRM(1:IGT)%G,
     +                        F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
           CALL CONEC(TRM(1:IGT)%RET,
     +                        F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
        ENDIF

C ... SCALAR EQUATIONS
      DO 4100 NS = 1,NSCAL
         CALL CONEC(FI(1,NS),F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
4100  CONTINUE
      ENDIF ! ISTFLX == 1

 5000 CONTINUE
      WRITE(3,*) 'Temperature ennen 2. boundia'
      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(ITURB == 9) THEN
             RNUT(IG1:IG1+NTOT(1,N)-1) = REPS(IG1:IG1+NTOT(1,N)-1)/
     +       RO(IG1:IG1+NTOT(1,N)-1)
             CALL PRINYS(3,RNUTC,REPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +       IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
             CALL PRINYS(3,NUTC,RNUT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +       IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF ! ITURB .EQ: 9
      ENDDO
      WRITE(3,*) 'menee boundiin toisen kerran',nblock
C***********************************************************************
      CALL BOUNDI(1)
C***********************************************************************

      WRITE(3,*) 'Values after the second CALL BOUND'
      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         IGN    = IG1 - 1 + NTOT(1,N)
         
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(1:NTOT(1,N)) = PRO(IG1:IGN)%TEMP(IPHASE)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF

         CALL PRINYS(3,PC,   P(IG1),   IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,ROC,  RO(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,  RM(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDDO
       
C ... END OF BLOCK LOOP 3

C ... GRAVITATION ARE UTILIZED AS INITIAL CONDITIONS
C ... BLOCK LOOP 4 BEGINS

      DO 5200 N = 1,NBLOCK

      IF(OMEGA(NPROCE(1+N,IPRO)) /= 0. .AND. INITC(N) >= 5
     +        .AND. INITC(N) <= 20) THEN

      ISTRID  = IMAX(1,N) + 2*IN
      JSTRID  = JMAX(1,N) + 2*JN
      KSTRID  = ISTRID*JSTRID
      NGL     = NPROCE(1+N,IPRO)
      IG1     = IG(1,N) - 1
    
      DO 5100 K = 1,KMAX(1,N)
      KA      = (KN+K-1)*KSTRID
      DO 5100 J = 1,JMAX(1,N)
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 5100 I = 1,IMAX(1,N)
      KK      = JJ + I + IG1
      UIROT   = A1XA(KK)*UROT(KK) + A1XA(KK+1)     *UROT(KK+1)      +
     +          A2XA(KK)*VROT(KK) + A2XA(KK+ISTRID)*VROT(KK+ISTRID) +
     +          A3XA(KK)*WROT(KK) + A3XA(KK+KSTRID)*WROT(KK+KSTRID)
      VIROT   = A1YA(KK)*UROT(KK) + A1YA(KK+1)     *UROT(KK+1)      +
     +          A2YA(KK)*VROT(KK) + A2YA(KK+ISTRID)*VROT(KK+ISTRID) +
     +          A3YA(KK)*WROT(KK) + A3YA(KK+KSTRID)*WROT(KK+KSTRID)
      WIROT   = A1ZA(KK)*UROT(KK) + A1ZA(KK+1)     *UROT(KK+1)      +
     +          A2ZA(KK)*VROT(KK) + A2ZA(KK+ISTRID)*VROT(KK+ISTRID) +
     +          A3ZA(KK)*WROT(KK) + A3ZA(KK+KSTRID)*WROT(KK+KSTRID)
      UIROT   = .5*UIROT
      VIROT   = .5*VIROT
      WIROT   = .5*WIROT
      U(KK)   = U(KK) + UIROT
      V(KK)   = V(KK) + VIROT
      W(KK)   = W(KK) + WIROT
      E(KK)   = E(KK)  + .5*RO(KK)*(UIROT**2 + VIROT**2 + WIROT**2)
      IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
      DO IPHASE = 1,NPHASE
         VAR(KK)%U(IPHASE) = VAR(KK)%U(IPHASE) + UIROT
         VAR(KK)%V(IPHASE) = VAR(KK)%V(IPHASE) + VIROT
         VAR(KK)%W(IPHASE) = VAR(KK)%W(IPHASE) + WIROT
      ENDDO
      ENDIF
5100  CONTINUE
      ENDIF   ! OMEGA /= 0. AND. INITC(N) >= 5

5200  CONTINUE
C ... END OF BLOCK LOOP 4
C **********************************************************************
      ENDIF   ! IOLD >= 0

      WRITE(3,*) 'Values before CALL STATIM'
      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         IGN    = IG1 - 1 + NTOT(1,N)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,PC,   P(IG1),   IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,ROC,  RO(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,  RM(IG1),  IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         F1RM(1:NTOT(1,N)) = PRO(IG1:IGN)%TEMP(IPHASE)
         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RM(1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF

      ENDDO

C ... NEW PRESSURES AND VELOCITIES AT THE FINE GRID LEVEL
C ... BLOCK LOOP 5 BEGINS

      M       = 1
      DO 4800 N = 1,NBLOCK
      IG2     = IG(M,N)
      IR2     = IR(M,N)
      IQ2     = IQ(M,N)
      IF2     = JF(M,N)
      IC2     = IC(M,N)
            
C ... CALCULATE VELOCITIES AND INTERNAL ENERGY 

C***********************************************************************
      CALL STATIM(1,N,2,.TRUE.)
C***********************************************************************

4800  CONTINUE
C ... END OF BLOCK LOOP 5

c     IF(TIMEL .OR. COORL) THEN ! INITIALIZE THE SECOND TIME-LEVEL
      IF(TIMEL) THEN ! INITIALIZE THE SECOND TIME-LEVEL
      
      CALL STARTI(NTOT,IG,NBLOCK,TIMEL,COORL,ITURB,NSCAL,MGM,
     +  IGRID(:,1),RO,RM,RN,RW,E,RK,REPS,FI,XCO,YCO,ZCO,
     +  ROLE2,RMLE2,RNLE2,RWLE2,ELE2,RKLE2,EPSLE2,FILE2,
     +  ROLE3,RMLE3,RNLE3,RWLE3,ELE3,RKLE3,EPSLE3,FILE3,
     +  XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,BLKS,VAR,MAXSB,IPRO,NPROCE,MBPRO,
     +  TWO_FLUIDL)

          DO I = 1,MAXTI ! And initialize the averages
             ROAV1(I)  = RO(I)
             RMAV1(I)  = RM(I)
             RNAV1(I)  = RN(I)
             RWAV1(I)  = RW(I)
             RKAV1(I)  = RK(I)
             EPSAV1(I) = REPS(I)
             EAV1(I)   = E(I)
             PAV1(I)   = P(I)
             TAV1(I)   = TEMP(I)

             ROAV2(I)  = RO(I)**2
             RMAV2(I)  = RM(I)**2
             RNAV2(I)  = RN(I)**2/(RO(I)+EPS)
             RWAV2(I)  = RW(I)**2/(RO(I)+EPS)
             EAV2(I)   = E(I)**2/(RO(I)+EPS)
             RUVAV(I)  = RM(I)*RN(I)/(RO(I)+EPS)
             RUWAV(I)  = RM(I)*RW(I)/(RO(I)+EPS)
             RVWAV(I)  = RN(I)*RW(I)/(RO(I)+EPS)
          ENDDO

      ENDIF ! TIMEL .OR. COORL

      RETURN
      END SUBROUTINE INIT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INITV

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : MAXB,IPRO,MGM,MAXSB,MBPRO,NPRO

      USE MAIN_ARRAYS, ONLY : IG,IMAX,JMAX,KMAX,F1R,F1RM,F1RN,
     &    F1RW,F1E,HAT1,F1RK,F1EPS,F1FI,NPROCE,ZZZ,P,TEMP,U,V,
     &    W,PDIFF,RK,REPS,FI,NTOT,BLKS,PRO,VAR,WAVEH,WH,JF,IC,
     &    IHF,ICON,NPATCH,HAT2,TRM
            
      USE NS3CO

      USE TYPE_ARRAYS

      IMPLICIT NONE
      
      INTEGER :: N,NGL,I,IG1,IG2,IGN,IGN2,IGT,NS,IPHASE,IC1,IF1,LEVEP1,
     +           KMAX2,NTO2

      CHARACTER(LEN=1) :: CLEVW, CLEVR

      LOGICAL :: THERE

cc      TYPE(PROPERTIES) PRO(*)
cc      TYPE(MPHASE_VARIABLES) VAR(*)

C ... INITIALIZATION  FROM THE PREVIOUS GRID LEVEL

C ... IF MPI IS USED MASTER PROCES READ COMPUT AND SEND IT TO SLAVES
    
      IF(PARALLEL .AND. IPRO == 1) THEN
         CALL COMRMP(NBLOCK,IG,IMAX,JMAX,KMAX,IN,JN,KN,MGM,F1R,
     +        F1RM,F1RN,F1RW,F1E,HAT1,F1RK,F1EPS,F1FI,MAXSB,NSCAL,ITURB,
     +        IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,VAR,0,
     +        16,MULPHL,BLKS,HAT2,TRANSL,TRM)
      ELSE
C ... ORDINARY ONE PROCESSOR MODE READING COMPUT OR WITH MPI SLAVES
C ... RECEIVE COMPUT FROM MASTER
         CALL COMSLA(NBLOCK,IG,IMAX,JMAX,KMAX,IN,JN,KN,MGM,F1R,
     +        F1RM,F1RN,F1RW,F1E,HAT1,F1RK,F1EPS,F1FI,MAXSB,NSCAL,ITURB,
     +        IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,VAR,0,
     +        16,MULPHL,BLKS,HAT2,TRANSL,TRM)
      ENDIF
C ... END OF BLOCK LOOP 1 INSIDE COMRMP AND COMSLA

      DO 3000 N = 1,NBLOCK

         NGL    = NPROCE(1+N,IPRO)
         IG2    = IG(2,N)
         IG1    = IG(1,N)
         IGN    = IG(1,N) + NTOT(1,N) - 1
         KMAX2  = MAX(1,KMAX(1,N)/2) ! Was a trap without the multigrid
         NTO2   = (IMAX(2,N)+2*IN)*(JMAX(2,N)+2*JN)*(KMAX2+2*KN)
         IGN2   = IG(1,N) + NTO2 - 1
         CALL ININIT( F1R(IG1),TEMP(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL ININIT(F1RM(IG1),U(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL ININIT(F1RN(IG1),V(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL ININIT(F1RW(IG1),W(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL ININIT( F1E(IG1),PDIFF(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
C ... turha
         CALL ININIT(HAT1(IG1),P(IG1),IMAX(2,N),JMAX(2,N),KMAX(2,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL ININIT( F1RK(IG1),RK(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            CALL ININIT(F1EPS(IG1),REPS(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         ENDIF

C ... For intermittency variables
         IF(TRANSL) THEN
            F1RM(IG1:IGN2)    = TRM(IG1:IGN2)%G
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            TRM(IG1:IGN)%G   = F1RN(IG1:IGN)

            F1RM(IG1:IGN2)  = TRM(IG1:IGN2)%RET
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            TRM(IG1:IGN)%RET = F1RN(IG1:IGN)
         ENDIF

C ... For multiphase modeling
         IF(MULPHL) THEN
            DO IPHASE = 1,BLKS(NGL)%NPHASE
            F1RM(IG1:IGN2) = VAR(IG1:IGN2)%F1R(IPHASE) ! Temperature (SP)
            F1RN(IG1:IGN)  = PRO(IG1:IGN)%TEMP(IPHASE)
c            CALL ININITm(F1RM(IG1),PRO(IG1:IGN)%TEMP(IPHASE),
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            PRO(IG1:IGN)%TEMP(IPHASE) = F1RN(IG1:IGN)
            PRO(IG1:IGN)%DTEMP(IPHASE) = F1RN(IG1:IGN)

            F1RM(IG1:IGN2) = VAR(IG1:IGN2)%F2R(IPHASE) ! Mass fraction
            F1RN(IG1:IGN)  = VAR(IG1:IGN)%X(IPHASE)
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            VAR(IG1:IGN)%X(IPHASE) = F1RN(IG1:IGN)

            IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
            F1RM(IG1:IGN2) = VAR(IG1:IGN2)%FRM(IPHASE) ! u-velocity
            F1RN(IG1:IGN)  = VAR(IG1:IGN)%U(IPHASE)
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            VAR(IG1:IGN)%U(IPHASE) = F1RN(IG1:IGN)

            F1RM(IG1:IGN2) = VAR(IG1:IGN2)%FRN(IPHASE) ! v-velocity
            F1RN(IG1:IGN)  = VAR(IG1:IGN)%V(IPHASE)
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            VAR(IG1:IGN)%V(IPHASE) = F1RN(IG1:IGN)

            F1RM(IG1:IGN2) = VAR(IG1:IGN2)%FRW(IPHASE) ! w-velocity
            F1RN(IG1:IGN)  = VAR(IG1:IGN)%W(IPHASE)
            CALL ININIT(F1RM(IG1),F1RN(IG1),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
            VAR(IG1:IGN)%W(IPHASE) = F1RN(IG1:IGN)
            ENDIF
            ENDDO
         ENDIF
C ... FOR SCALAR EQ. PPR 1.3
         DO NS = 1,NSCAL
            CALL ININIT(F1FI(1,NS),FI(IG1,NS),IMAX(2,N),JMAX(2,N),
     +           KMAX(2,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)
         ENDDO
         DO I = 1,NTOT(1,N)
            P(I+IG1-1) = BLKS(NGL)%FRSPRE + PDIFF(I+IG1-1)
         ENDDO
      
 3000 CONTINUE

      WRITE(4,'(/,2X,A)')  'Initialization performed from a previous gri
     +d level'
      WRITE(4,'(2X,A/)')   '============================================
     +======='
      WRITE(45,'(/,2X,A)') 'Initialization from a previous grid level'
      WRITE(45,'(2X,A/)')  '========================================='
      WRITE(13,'(/,2X,A)') 'Initialization from a previous grid level'
         
      IF(FRESUL) THEN

         LEVEP1 = LEVEL + 1
         CALL NUMCH1(CLEVR,LEVEP1)
       
         INQUIRE(FILE='WH.DAT'//CLEVR,EXIST=THERE)

         IF(THERE) THEN         ! Read the free-surface file

*         CALL READWH(IPRO,PARALLEL,.TRUE.)
         CALL READWH(.TRUE.)
    
         ELSE 
           WRITE(*,*)  ' No initialization for the free surface.'
           WRITE(45,*) ' No initialization for the free surface.'
           WRITE(4 ,*) ' No initialization for the free surface.'
         ENDIF
    
         IF(IFSBC == 1) THEN ! Use the WH array
           DO N = 1,NBLOCK
           IF1  = JF(1,N)
           IC1  = IC(1,N)
           CALL FS_SETWAVEHBACK(WAVEH,WH(IF1),IHF(1,1),ICON(IC1),
     +     NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ICYCLE)
           INWH = 1
           ENDDO
         ENDIF
      CLOSE(21) ! The coarse grid file

      CALL NUMCH1(CLEVW,LEVEL)
      OPEN(21,FILE='WH.DAT'//CLEVW,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF

c       CALL WRITWH
      REWIND(16)
      RETURN
      END SUBROUTINE INITV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRIVAL(IALL,LDISTAN,CALSUB)

      USE MPI

      USE CHARACTERS
      
      USE INTEGERS, ONLY : MAXB,MAXEB,MAXSB,IBF,IPRO,MGM,NB,
     &    NPRO,MAXSS,IREPEA,NREPEA,NBGG

      USE TYPE_ARRAYS
      
      USE MAIN_ARRAYS, ONLY : VOL,A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,
     &    A2YA,A2ZA,A3XA,A3YA,A3ZA,XC,YC,ZC,UROT,VROT,WROT,TWALL,
     &    IMAX,JMAX,KMAX,MGRID,IG,D1,D2,D3,UWALL,VWALL,WWALL,HFLUX,
     &    XCP,YCP,ZCP,F1R,F1RM,F1RN,NORMAL,NSPB,NSPP,WMFLUX,POROS,
     &    RMLOSS,WHSTAG,WTEMP,RSDIRX,RSDIRY,RSDIRZ,RBK,DISTW,BLANK,NTOT,
     &    ICON,IMAXG,JMAXG,KMAXG,NPATCH,NLOCAL,NPNUM,LOCDIS,JET,JLOC,
     &    JTRA,APP,NPROCE,IGRID,IC,XCO,YCO,ZCO,VOLN,OMEGAX,OMEGAY,
     &    OMEGAZ,OMEGA,CENAX,CENAY,CENAZ,XCOR,DRO,DM,DN,DW,DE,IT,IL,
     &    IK,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,RCON,IHF,
     &    JF,NSOLPA,IDIMS,JDIMS,KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,ISTRS,
     &    JSTRS,KSTRS,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNMF,
     &    BOUNT,BOUNP,BOUNPD,BOUNRK,BOUNEP,BOUNA1,BOUNA2,BOUNFI,BOUNBI,
     &    INITC,JSTATE,NSPG,PRO,VAR,BLKS,APATCH,U,V,W,EPS2,VIST,RO,TEMP,
     &    P,PDIFF,PTUR,E,RM,RN,RW,RK,REPS,FI,S11,BIJ,TIJ,XGRI,YGRI,ZGRI,
     &    XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,WTRAN,VTRAN,XFC,YFC,ZFC,SDI,
     &    BOUNG,BOUNRET,TRM,IDIMSG,JDIMSG,SURFA,SURFT,
     &    UROTCP,VROTCP,WROTCP,BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,
     &    IUPPT,BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2
 
      USE NS3CO

      USE FLIGHT, ONLY : ACTDISK,MVSHIP,IFA,IFR,IFT

      IMPLICIT NONE

      REAL    :: APU(NBGG), FRSALFA(3), SKEWAVE_MIN, SKEW_MIN, VOLNA

      INTEGER :: IALL,M,N,IG1,IC1,NGL,ERRORCODE,IERR,IF1,KSTATE,
     &           IG0,IP,LLL,IG2,IG3,KU1,KU2,KU3,KU4,KU5,KU6,IGN,ITYPE,
     &           IBTYPE,NIOS,IPL,IPG,IFACE,IY1,III,iiirc,IMAXP1,JMAXP1,
     &           KMAXP1

      LOGICAL :: NEGVOL,LDISTAN,CHANGE_GRID,INOUTLX

      LOGICAL :: WDISTS_NEEDED, WDISTS_FORCED, MASTER, SKEW_WARNING
           
      CHARACTER(6) :: CALSUB

      REAL, ALLOCATABLE, DIMENSION(:) :: BXMINW, BXMAXW, BYMINW, BYMAXW, 
     &                                   BZMINW, BZMAXW

      NEGVOL       = .FALSE.
      CHANGE_GRID  = .FALSE.
      SKEW_WARNING = .FALSE.

      MASTER = IPRO == 1

      BXMIN(1:NBLOCG) = 0.0
      BXMAX(1:NBLOCG) = 0.0
      BYMIN(1:NBLOCG) = 0.0
      BYMAX(1:NBLOCG) = 0.0
      BZMIN(1:NBLOCG) = 0.0
      BZMAX(1:NBLOCG) = 0.0

C ... BLOCK LOOP 1 BEGINS

      DO 2100 N = 1,NBLOCK

      NGL = NPROCE(1+N,IPRO)

      IF(IGRID(NGL,1) >= 1 .AND. IGRID(NGL,1) <= 9) THEN
         CHANGE_GRID = .TRUE.              ! Rotation              
      ELSE IF(IGRID(NGL,1) >= 11 .AND. IGRID(NGL,1) <= 19 .AND.
     &        TIMEL) THEN                  ! Flying objects
         CHANGE_GRID = .TRUE. 
      ELSE IF(IGRID(NGL,1) >= 21 .AND. IGRID(NGL,1) <= 29) THEN  
         CHANGE_GRID = .TRUE.              ! Moving ships
      ELSE IF(IGRID(NGL,1) >= 31 .AND. IGRID(NGL,1) <= 39) THEN 
         CHANGE_GRID = .TRUE.              ! Helicopter blades  
      ELSE IF(IGRID(NGL,2) >= 0) THEN
         CHANGE_GRID = .TRUE.              ! Grid deformation    
      ENDIF

      IF(IALL == 0 .OR. CHANGE_GRID) THEN
         IG1    = IG(1,N)
         IC1    = IC(1,N)
         IY1    = 1
         IF(COORL) IY1 = IG1

C *** CALCULATE THE GRID PROPERTIES ************************************

         CALL GRID1B(N,A1XA(IG1),A1YA(IG1),A1ZA(IG1),A1(IG1),
     1        A2XA(IG1),A2YA(IG1),A2ZA(IG1),A2(IG1),
     2        A3XA(IG1),A3YA(IG1),A3ZA(IG1),A3(IG1),
     3        VOL(IG1),XC(IG1),YC(IG1),ZC(IG1),UROT(IG1),VROT(IG1),
     4        WROT(IG1),XCO(IG1),YCO(IG1),ZCO(IG1),
C ... On this line the possible face centerpoint array
     5        IMAX(1,N),JMAX(1,N),KMAX(1,N),IG(1,N),MGRID(N),
     6        LEVEL,ICON(IC1),NPATCH(N),GRILEN,VOLN(N),OMEGAX(NGL),
     7        OMEGAY(NGL),OMEGAZ(NGL),CENAX(NGL),CENAY(NGL),CENAZ(NGL),
     8        XCOR(1,N),NEGVOL,PRN,IREPEA,NREPEA,
     9        XLE2(IY1),YLE2(IY1),ZLE2(IY1),XLE3(IY1),YLE3(IY1),
     1        ZLE3(IY1),IGRID(NGL,1),IGRID(NGL,2),DT,DTOLD,TIMEL,
     2        XFC(IG1),YFC(IG1),ZFC(IG1))

         CALL BOUNDX(XC(IG1),YC(IG1),ZC(IG1),IMAX(1,N),JMAX(1,N),
     1        KMAX(1,N),1) !     CENTERPOINTS FOR THE FREE-STREAM VOLUMES

*         WRITE(45,*)
*         WRITE(45,'(1X,A,I4)') 'SKEWNESS VALUES FOR GLOBAL BLOCK:',NGL
*         WRITE(45,*)
*         WRITE(45,*) 'SKEWMIN(I)  SKEWMIN(J)  SKEWMIN(K)    SKEWAVE(I)',
*     +               '  SKEWAVE(J)  SKEWAVE(K)'
*         WRITE(45,*) '================================================',
*     +               '========================'

         CALL SKEWNESS(XFC(IG1),YFC(IG1),ZFC(IG1),XC(IG1),YC(IG1),
     +        ZC(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),
     +        A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),M,NGL,BLKS,SDI(IG1))
         CALL BOUNDV(VOL(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N)) ! POSITIVITY
         CALL D1D2D3(N,IMAX,JMAX,KMAX,MGRID,IG,A1,A2,A3,VOL,D1,D2,D3,
     +        MGM)
      ENDIF ! IALL
 2100 CONTINUE
      
C ... END OF BLOCK LOOP 1

C ... Stop if negative volumes are found

      IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR) ! Wait ...

      IF(NEGVOL) THEN
         WRITE(*,*) 'FROM ',CALSUB,': Negative volumes found. See file',
     +   ' NEGVOL. Exiting...'
         WRITE(45,*)'FROM ',CALSUB,': Negative volumes found. See file',
     +   ' NEGVOL. Exiting...'
         WRITE(46,*)'FROM ',CALSUB,': Negative volumes found. See file',
     +   ' NEGVOL. Exiting...'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF     
     
      IF(NEGVL .AND. MASTER .AND. ICYCLE <= ICYOLD) THEN
         WRITE(*,*)''
         WRITE(*,*)'ATTENTION: Negative volumes allowed!'
         WRITE(*,*)'VOLUMES',NEGV,'... 1.E-11 are changed to 1.E-11'
         WRITE(*,*)''
         WRITE(45,*)'ATTENTION: Negative volumes allowed!'
         WRITE(45,*)'Volumes',NEGV,'... 1.E-11 -> 1.E-11'
         WRITE(46,*)'ATTENTION: Negative volumes allowed!'
         WRITE(46,*)'Volumes',NEGV,'... 1.E-11 -> 1.E-11'
      ENDIF
     
C ... Analyze skewness values

      IF(PARALLEL) THEN

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWI,APU(1),NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWI = APU(1:NBLOCG)

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWJ,APU,NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWJ = APU(1:NBLOCG)

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWK,APU,NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWK = APU(1:NBLOCG)

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWIMIN,APU,NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWIMIN = APU(1:NBLOCG)

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWJMIN,APU,NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWJMIN = APU(1:NBLOCG)

         APU = 0.0
         CALL MPI_REDUCE(BLKS%SKEWKMIN,APU,NBLOCG,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         IF(MASTER) BLKS(1:NBLOCG)%SKEWKMIN = APU(1:NBLOCG)

         APU = 0.0

      ENDIF

      IF(MASTER) THEN ! Output to MEMORY file.

         WRITE(45,'(3/)')
         IF(PARALLEL) THEN
            WRITE(45,*) 'GLOBAL GRID SKEWNESS ANALYSIS:'
         ELSE
            WRITE(45,*) 'GRID SKEWNESS ANALYSIS:'
         ENDIF
         WRITE(45,*)
         IF(PARALLEL) WRITE(45,*) 'Global'
         WRITE(45,*) 'Block ID   SKEWIAVE  SKEWJAVE  SKEWKAVE  ',
     &                          'SKEWIMIN  SKEWJMIN  SKEWKMIN  '   
         WRITE(45,*) '=======================================',
     &                        '=============================='   
         DO III = 1,NBLOCG
            WRITE(45,'(I6,4X,6F10.6)') III,
     &      BLKS(III)%SKEWI,   BLKS(III)%SKEWJ,   BLKS(III)%SKEWK,
     &      BLKS(III)%SKEWIMIN,BLKS(III)%SKEWJMIN,BLKS(III)%SKEWKMIN
            IF(IALL == 0) THEN
            SKEWAVE_MIN = MIN(BLKS(III)%SKEWI,BLKS(III)%SKEWJ,
     &                     BLKS(III)%SKEWK)
            SKEW_MIN    = MIN(BLKS(III)%SKEWIMIN,BLKS(III)%SKEWJMIN,
     &                     BLKS(III)%SKEWKMIN)
            IF(SKEWAVE_MIN < 0.5 .AND. SKEWAVE_MIN >= 0.3) THEN
               WRITE(13,'(2A,I3,A)') '  An average skewness value is',
     &         ' low for block',III,'. See MEMORY-files.'
            ENDIF
            IF(SKEWAVE_MIN < 0.3) THEN
               WRITE(13,'(2A,I3,A/A)') '  An average skewness value is',
     &         '  very low for block',III,'. See MEMORY-files.',
     &         '  A grid-generation course is recommended'
               SKEW_WARNING = .TRUE.
            ENDIF
            IF(SKEW_MIN < 0.1 .AND. SKEW_MIN >= 0.05) THEN
               WRITE(13,'(2A,I3,A)') '  A minimum skewness value is',
     &         ' low for block',III,'. See MEMORY-files.'
            ENDIF
            IF(SKEW_MIN < 0.05 .AND. SKEW_MIN >= 0.01) THEN
               WRITE(13,'(2A,I3,A)') '  A minimum skewness value is',
     &         ' very low for block',III,'. See MEMORY-files.'
               SKEW_WARNING = .TRUE.
            ENDIF
            IF(SKEW_MIN < 0.01) THEN
               WRITE(13,'(2A,I3,A)') '  A minimum skewness value is',
     &         ' almost catasthrophic for block',III,
     &         '. See MEMORY-files.'
               SKEW_WARNING = .TRUE.
            ENDIF
            ENDIF ! IALL
         ENDDO
         WRITE(45,'(3/)')
      ENDIF ! MASTER

      IF(SKEW_WARNING) WRITE(*,*)' GRID1B :  There are skewness',
     &                           ' warnings in RUN.LOG(001) -file.'


C ... BLOCK LOOP 2 BEGINS
          
      IF(IALL == 0) NSPTOT  = 0

      DO 2200 N = 1,NBLOCK

      NGL    = NPROCE(1+N,IPRO)

      IF(IGRID(NGL,1) >= 1 .AND. IGRID(NGL,1) <= 9) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 11 .AND. IGRID(NGL,1) <= 19 .AND.
     +        TIMEL) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 21 .AND. IGRID(NGL,1) <= 29) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 31 .AND. IGRID(NGL,1) <= 39) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,2) >= 0) THEN
         CHANGE_GRID = .TRUE.
      ENDIF  ! onko eo. maarittelyt oiken tassa kohtaa?
     
      IF(IALL == 0 .OR. CHANGE_GRID) THEN

      DO 3210 M = 1,MGRID(N)

      NGL     = NPROCE(1+N,IPRO) ! A global block number
      IG1     = IG(M,N)
      IF1     = JF(M,N)
      IC1     = IC(M,N)
      IY1     = 1
      IF(COORL) IY1 = IG1
      KSTATE  = JSTATE(NGL,1)
       
      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSVEL  = BLKS(NGL)%FRSVEL
      FRSVIS  = BLKS(NGL)%FRSVIS
      FRSTEM  = BLKS(NGL)%FRSTEM
      TEMINI  = BLKS(NGL)%TEMINI
      FRSSIE  = BLKS(NGL)%FRSSIE
      SIEINI  = BLKS(NGL)%SIEINI
      CHLREF  = BLKS(NGL)%CHLREF
      RMUINI  = BLKS(NGL)%RMUINI
      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      FRSMUT  = BLKS(NGL)%FRSMUT
      SOLTEM  = BLKS(NGL)%SOLTEM
      GVEX    = BLKS(NGL)%GVEX
      GVEY    = BLKS(NGL)%GVEY
      GVEZ    = BLKS(NGL)%GVEZ
      
      DO LLL  = 1,NPHASES ! Antakee mersu
      FRSALFA(LLL) = BLKS(NGL)%FRSALFA(LLL)
      ENDDO
       
      NPHASE  = BLKS(NGL)%NPHASE

C ... LUMP CORNER POINTS TO THE LOWER GRID LEVELS

      IF(M >= 2) THEN
         IG0 = IG(M-1,N)
         CALL LUMXCO(XCO(IG1),YCO(IG1),ZCO(IG1),XCO(IG0),YCO(IG0),
     +        ZCO(IG0),IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))
      ENDIF

C ... Rotating reference frame

      IF(.NOT. TIMEL) THEN
         SELECT CASE (IGRID(NGL,1))
            CASE (1,10,21,22,23,24,25,26,27,28,29,31,37,38)
               CALL ROTFUN(UROT(IG1),VROT(IG1),WROT(IG1),NTOT(M,N),
     +              OMEGA(NPROCE(1+N,IPRO)))
            CASE DEFAULT
               CALL ROTFUN(UROT(IG1),VROT(IG1),WROT(IG1),NTOT(M,N),0.0)
         END SELECT
      ENDIF

      IF(TIMEL) THEN
         SELECT CASE (IGRID(NGL,1))
            CASE (1,10,31)
               CALL ROTFUN(UROT(IG1),VROT(IG1),WROT(IG1),NTOT(M,N),
     +              OMEGA(NPROCE(1+N,IPRO)))
            CASE (0,21,22)
               IF(IGRID(NGL,2) == 0) THEN
               CALL ROTFUN(UROT(IG1),VROT(IG1),WROT(IG1),NTOT(M,N),0.0)
               ENDIF   
         END SELECT
      ENDIF

      
      IF (IPRINT == ICYCLE) THEN
       CALL PRINYS(3,UROTC,UROT(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF

C ... The grid is moving
        
      IF(SQRT(GVEX**2 + GVEY**2 + GVEZ**2) >= 1.E-5) THEN
         CALL VELFUN(UROT(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     +        GVEX,GVEY,GVEZ,NTOT(M,N))
         CALL VELFUN(VROT(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),
     +        GVEX,GVEY,GVEZ,NTOT(M,N))
         CALL VELFUN(WROT(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +        GVEX,GVEY,GVEZ,NTOT(M,N))
      ENDIF
      
      IF(IGRID(NGL,1) /= 0.) THEN
       
      CALL CHECK(UROT(IG1),VROT(IG1),WROT(IG1),IMAX(M,N),JMAX(M,N),
     + KMAX(M,N),A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),
     + A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),A1(IG1),
     + A2(IG1),A3(IG1),DRO(IG1),DE(IG1),DM(IG1),DN(IG1))

      IF(IPRINT == ICYCLE) THEN
       CALL PRINYS(3,SURFBXC,DRO(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
       CALL PRINYS(3,SURFBYC,DE(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
       CALL PRINYS(3,SURFBZC,DM(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
       CALL PRINYS(3,UROTBC,DN(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF
      ENDIF

C ... CALCULATE SURFACE VELOCITIES

      CALL SURPRO(UWALL,VWALL,WWALL,
     +     XCP,YCP,ZCP,XCO(IG1),YCO(IG1),ZCO(IG1),
     +     A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),
     +     A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +     A1(IG1),A2(IG1),A3(IG1),
     +     OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     +     CENAX(NGL),CENAY(NGL),CENAZ(NGL),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +     IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,N,M,
     +     NPATCH(N),ICON(IC1),RCON(IC1),IHF(1,M),
     +     NSPTOT,NSOLPA(N),NORMAL,IDIMS,JDIMS,KX1S,KX2S,
     +     KY1S,KY2S,KZ1S,KZ2S,ISTRS,JSTRS,KSTRS,NSPB,NSPP,
     +     XC(IG1),YC(IG1),ZC(IG1),NSPG,APATCH,
     +     IREPEA,NREPEA,XLE2(IY1),YLE2(IY1),ZLE2(IY1),XLE3(IY1),
     +     YLE3(IY1),ZLE3(IY1),DT,TIMEL,IGRID(NGL,1),POROS,
     +     RMLOSS,IDIMSG,JDIMSG,BOUNDN,IGRID(NGL,3),SURFA,SURFT,
     +     IALL,IPRO,DTOLD,GVEX,GVEY,GVEZ)


      IF(IALL == 0 .AND. M == 1 .AND. NSOLPA(N) > 0) THEN
         NSPTOT = NSPTOT + NSOLPA(N)
         IF(IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*) 'Total number of solids =',NSPTOT
         WRITE(45,*) '===================================='
         ENDIF ! IREPEA(3) <= 10
      ENDIF

      CALL SURHTP(UWALL,TWALL,HFLUX,WMFLUX,POROS,
     +     WHSTAG,WTEMP,RSDIRX,RSDIRY,RSDIRZ,RBK,
     +     A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),
     +     A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),IBOT(M,N),ITOP(M,N),
     +     JBOT(M,N),JTOP(M,N),KBOT(M,N),KTOP(M,N),
     +     IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,N,M,
     +     NPATCH(N),ICON(IC1),RCON(IC1),IHF(1,M),
     +     NSOLPA(N),ISTRS,JSTRS,KSTRS,NSPB,
     +     KSTATE,MAXB,IREPEA,NREPEA,WTRAN,ngl)

      CALL BOUNDP(UWALL,VWALL,WWALL,TWALL,HFLUX,WMFLUX,POROS,
     +     WHSTAG,WTEMP,RBK,
     +     XCP,YCP,ZCP,XCO(IG1),YCO(IG1),ZCO(IG1),VOL(IG1),
     +     A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),
     +     A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +     A1(IG1),A2(IG1),A3(IG1),UROT(IG1),VROT(IG1),WROT(IG1),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),IBOT(M,N),ITOP(M,N),
     +     JBOT(M,N),JTOP(M,N),KBOT(M,N),KTOP(M,N),
     +     IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,NGL,N,M,
     +     NPATCH(N),ICON(IC1),RCON(IC1),IHF(1,M),
     +     NSPTOT,NSOLPA(N),NORMAL,IDIMS,JDIMS,KX1S,KX2S,
     +     KY1S,KY2S,KZ1S,KZ2S,ISTRS,JSTRS,KSTRS,NSPB,NSPP,
     +     KSTATE,JSTATE,NBLOCG,MAXB,XC(IG1),YC(IG1),ZC(IG1),
     +     BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNMF,BOUNT,BOUNP,BOUNRK,
     +     BOUNEP,FRSVEL,ALPHA,BETA,FRSDEN,FRSSIE,FRSTEM,FRSPRE,RKLIM,
     +     EPSLIM,ITURB,INITC(N),BOUNDN,BOUNPD,BOUNA1,BOUNA2,BOUNFI,
     +     BOUNBI,LEVEL,NSCAL,MAXEB,MAXSB,IBF,ISTRES,RGAS,GAMMA,E0REF,
     +     T0REF,MULPHL,NPHASE,FRSALFA,APATCH,IREPEA,NREPEA,
     +     BOUNG,BOUNRET,TRANSL,BLKS(NGL)%FRSTUR,BLKS(NGL)%FRSVIS,
     +     IUPPT(NGL),BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2,
     +     BLKS(NGL)%FRADEN)

3210  CONTINUE
      ENDIF
2200  CONTINUE

      CLOSE(74)

C ... END OF BLOCK LOOP 2

C ... Update boundaries

      CALL CONEPA_SP(WTRAN,'SOL')

C ... Update transition data in the volumes

      DO 4200 N = 1,NBLOCK !  Block loop 3

      NGL    = NPROCE(1+N,IPRO)

      IF(IGRID(NGL,1) >= 1 .AND. IGRID(NGL,1) <= 9) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 11 .AND. IGRID(NGL,1) <= 19 .AND. 
     +        TIMEL) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 21 .AND. IGRID(NGL,1) <= 29) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,1) >= 31 .AND. IGRID(NGL,1) <= 39) THEN
         CHANGE_GRID = .TRUE.
      ELSE IF(IGRID(NGL,2) >= 0) THEN
         CHANGE_GRID = .TRUE.
      ENDIF  ! onko eo. maarittelyt oikein? Probably yes.
     
      IF(IALL == 0 .OR. CHANGE_GRID) THEN

      DO 4210 M = 1,MGRID(N)

      IG1     = IG(M,N)
      IC1     = IC(M,N)

      CALL SURTRA(VTRAN(IG1),WTRAN,NPATCH(N),ICON(IC1),IHF(1,M),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN)

4210  CONTINUE

      ENDIF ! IALL == 0 .OR. CHANGE_GRID

4200  CONTINUE ! End of block loop 3

C ... Start of a solid patch loop

      IF(NSPTOT > 0 .AND. IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*)
         WRITE(45,*) 'Information for the solid-like patches'
         WRITE(45,*) '======================================'
         WRITE(45,9878)
         WRITE(45,9460)
          
         DO 2300 IP = 1,NSPTOT

         WRITE(45,9877) IP,IDIMS(IP),JDIMS(IP),KX1S(IP),KX2S(IP),
     &        KY1S(IP),KY2S(IP),KZ1S(IP),KZ2S(IP),ISTRS(IP),JSTRS(IP),
     &        KSTRS(IP),NORMAL(IP),NSPB(IP),NSPP(IP),NSPG(IP),
c     &        APATCH(IP)
     &        APATCH(NSPP(IP))
2300     CONTINUE
      ENDIF ! NSPTOT > 0

9460  FORMAT(3X,97('='))
9877  FORMAT(1X,9I5,2I7,2I5,1X,3I5,E12.4)
9878  FORMAT(/4X,'ISP IDIM JDIM  XLO  XUP  YLO  YUP  ZLO  ZUP ',
     &' ISTR   JSTR KSTR NORM  BLOCK PATCH G.PATCH  AREA')

C ... Start of a inlet/outlet patch loop


      IF(IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*)
         WRITE(45,*) 'Information for the inlet/outlet patches'
         WRITE(45,*) '========================================'
         WRITE(45,9888)
         WRITE(45,9470)

         NIOS   = 0
           
         DO 2400 N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IC1    = IC(1,N) - 1
         DO IP  = 1,NPATCH(N)
         ITYPE  = ICON(IC1+(IP-1)*IC9 + 1)   ! Boundary type
         IBTYPE = ICON(IC1+(IP-1)*IC9 + 8)   ! Inlet/outlet type
         IPL    = ICON(IC1+(IP-1)*IC9 + 2)   ! Local patch number
         IPG    = ICON(IC1+(IP-1)*IC9 + 25)  ! Global patch number
         IFACE  = ICON(IC1+(IP-1)*IC9 + 3)   ! Block face

         IF(ITYPE == 3 .OR. ITYPE == 5) THEN ! BCP
         NIOS   = NIOS + 1
         WRITE(45,9887) IP,ITYPE,IBTYPE,N,NGL,IFACE,IPL,IPG,
     &        APATCH(IPL)
         ENDIF
         ENDDO
2400     CONTINUE


         IF(NIOS == 0) THEN 
             WRITE(45,*)
             WRITE(45,*) 'No inlet/outlet surfaces found'
             INOUTL = .FALSE.
         ELSE
             WRITE(45,*)
             WRITE(45,9471)NIOS
             WRITE(4,*)
             WRITE(4,9471) NIOS
             WRITE(4,*)
             WRITE(4,*) 'Information for the inlet/outlet patches'
             WRITE(4,*) '========================================'
             WRITE(4,9888)
             WRITE(4,9470)
             INOUTL = .TRUE.

             DO 2500 N = 1,NBLOCK
             NGL    = NPROCE(1+N,IPRO)
             IC1    = IC(1,N) - 1
             DO IP  = 1,NPATCH(N)
             ITYPE  = ICON(IC1+(IP-1)*IC9 + 1) ! Boundary type
             IBTYPE = ICON(IC1+(IP-1)*IC9 + 8) ! Inlet/outlet type
             IPL    = ICON(IC1+(IP-1)*IC9 + 2) ! Local patch number
             IPG    = ICON(IC1+(IP-1)*IC9 + 25)! Global patch number
             IFACE  = ICON(IC1+(IP-1)*IC9 + 3) ! Block face
             IF(ITYPE == 3 .OR. ITYPE == 5) THEN ! BCP
             WRITE(4,9887) IP,ITYPE,IBTYPE,N,NGL,IFACE,IPL,IPG, ! BCP
     &        APATCH(IPL)
             ENDIF
             ENDDO

2500         CONTINUE
        ENDIF ! NIOS == 0
      ENDIF ! IREPEA(3) <= 10

      IF (PARALLEL) THEN         ! Knowledge about INOUTL is broadcasted
            CALL MPI_REDUCE(INOUTL,INOUTLX,1,MPI_LOGICAL,MPI_LOR,0,
     &           MPI_COMM_WORLD,IERR)
            IF (IPRO == 1) INOUTL = INOUTLX
            CALL MPI_BCAST(INOUTL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      ENDIF ! PARALLEL

9470  FORMAT(3X,75('='))
9471  FORMAT(I4,' inlet/outlet surfaces found')
9887  FORMAT(1X,I5,7(I5,3X),E12.4) 
9888  FORMAT(/4X,'IP  ITYPE BC-TYPE  BLOCK    NGL     FACE   PATCH',
     &'   GLOBAL    AREA')
      WRITE(45,*)
     

C ... Calculate the wall distances.

      WDISTS_NEEDED = LDISTAN .OR.
     &                ITURB == 6 .OR.
     &                ITURB == 9 .OR. 
     &                NCHIM > 0
      
      WDISTS_FORCED = LDISTAN .OR. TIMEL .OR. IALL /= 0

C ... Cancel the wall distance calculation in next situations:
      
      IF(.NOT.LDISTL .AND. IREPEA(2) >= 1) WDISTS_NEEDED = .FALSE.

      IF(ACTDISK .AND. .NOT. MVSHIP .AND. IREPEA(2) >= 1)
     &   WDISTS_NEEDED = .FALSE.

      IF (IALL == 0 .AND. COORL .AND. IOLD <= 0 .OR.
     &    IALL == 0 .AND. COORL .AND. TIMEL) WDISTS_NEEDED = .FALSE.

*      WDISTS_NEEDED = MPCCIL .OR. MODALFSIL
*      WDISTS_FORCED = MPCCIL .OR. MODALFSIL

      IF (WDISTS_NEEDED) THEN
*         write(*,*) CALSUB
         CALL ALLCHI   ! Collect skin elements.

         IF(.NOT.TRUE_DISTL) THEN
         CALL DISTANCES(XC,YC,ZC,DISTW,F1R,F1RM,F1RN,IMAX,JMAX,KMAX,
     &        NTOT,IG,MGRID,IN,JN,KN,ICON,IC9,NPATCH,NB,MGM,NBLOCK,
     &        LEVEL,30,IPRO-1,IMAXG,JMAXG,KMAXG,NLOCAL,NBLOCG,
     &        NPNUM,NPRO,LOCDIS,JET,IREPEA,NREPEA,WDISTS_FORCED)
         ELSE
           CALL TRUE_DISTANCES(30,WDISTS_FORCED)
         ENDIF

      ENDIF  ! WDISTS_NEEDED


C **** GEOMETRY ON THE BOUNDARIES IS CONNECTED **********************

      DO 8000 LLL = 1,2
      DO 8000 M = 1,1
C ... Does not connect although M is increased
         CALL CONEC(VOL,   F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL CONEC(D1,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL CONEC(D2,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL CONEC(D3,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL CONEC(DISTW, F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL CONEC(VTRAN, F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
         CALL MIRC(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A1,A2,A3,A1XA,A1YA,
     +        A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,1,NBLOCK) 
c         CALL CONEC(XC,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,1)
c         CALL CONEC(YC,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,1)
c         CALL CONEC(ZC,    F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,1)
C ... Pellaako ? TSii 9.5.2003

c      IF(NBLOCK > 1) THEN
         CALL CONECX(XC,M,1,NBLOCK,1)
         CALL CONECX(YC,M,1,NBLOCK,1)
         CALL CONECX(ZC,M,1,NBLOCK,1)
c      ENDIF
8000  CONTINUE

C ... ROTATE GEOMETRIC VALUES IN SLD > SLIMES 
C ... ADDED 1.9.2005 ---------------------------

      CALL SLD(U,V,W,EPS2,VIST,RO,TEMP,P,PDIFF,PTUR,E,RM,RN,RW,RK,REPS,
     +     FI,S11,BIJ,TIJ,PRO,VAR,TRM,F1R,F1RM,F1RN,MAXB,MAXEB,MAXSB,
     +     MAXSS,XCO,YCO,ZCO,M,1,NBLOCK,ISTRES,D1,D2,D3,VOL,DM,DN,DW,3)

C ... BLOCK LOOP 3 BEGINS

      DO 8101 N = 1,NBLOCK
      IF(IALL == 0 .OR. IGRID((NPROCE(1+N,IPRO)),1) /= 0) THEN
      NGL    = NPROCE(1+N,IPRO)
      IMAXP1 = IMAX(1,N) + 1
      JMAXP1 = JMAX(1,N) + 1
      KMAXP1 = KMAX(1,N) + 1
      DO 8100 M = 1,MGRID(N)
         IG2     = IG(M,N)
         IG3     = IG(M+1,N)
         IF(M < MGRID(N)) THEN

            CALL BOUNDX(XC(IG3),YC(IG3),ZC(IG3),IMAX(M+1,N),JMAX(M+1,N),
     +      KMAX(M+1,N),M+1) !  CENTERPOINTS FOR THE FREE-STREAM VOLUMES
            CALL TRANSV(VOL(IG3),VOL(IG2),IMAX(M+1,N),JMAX(M+1,N),
     +      KMAX(M+1,N),KMAX(M,N))
c            CALL SKEWNESS(XFC(IG3),YFC(IG3),ZFC(IG3),XC(IG3),YC(IG3),
c     +      ZC(IG3),A1XA(IG3),A1YA(IG3),A1ZA(IG3),A2XA(IG3),
c     +      A2YA(IG3),A2ZA(IG3),A3XA(IG3),A3YA(IG3),A3ZA(IG3),
c     +      IMAX(M+1,N),JMAX(M+1,N),KMAX(M+1,N),M+1,NGL,BLKS,SDI(IG3))
         ENDIF
 8100 CONTINUE
      ENDIF

      IG2     = IG(1,N)
      WRITE(45,*) 
      WRITE(45,*) ' Check for negative volumes after connec ',NGL
      CALL CHECKV(VOL(IG2),IMAX(1,N),JMAX(1,N),KMAX(1,N),N,VOLNA,
     2             NEGVOL,PRN,IREPEA,NREPEA)

 8101 CONTINUE

C ... END OF BLOCK LOOP 3

C ... Block bounding box (2nd ghost cell center points included)

      DO N = 1,NBLOCK
         IG1 = IG(1,N)
         NGL = NPROCE(1+N,IPRO)
         CALL BOUNDS(IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &               XC(IG1),YC(IG1),ZC(IG1),
     &               BXMIN(NGL),BXMAX(NGL),
     &               BYMIN(NGL),BYMAX(NGL),
     &               BZMIN(NGL),BZMAX(NGL),IN,JN,KN)
      ENDDO

C ... Collect block bounding boxes

      IF(PARALLEL) THEN

         ALLOCATE(BXMINW(NBLOCG))
         ALLOCATE(BXMAXW(NBLOCG))
         ALLOCATE(BYMINW(NBLOCG))
         ALLOCATE(BYMAXW(NBLOCG))
         ALLOCATE(BZMINW(NBLOCG))
         ALLOCATE(BZMAXW(NBLOCG))

         BXMINW(1:NBLOCG) = BXMIN(1:NBLOCG)
         BXMAXW(1:NBLOCG) = BXMAX(1:NBLOCG)
         BYMINW(1:NBLOCG) = BYMIN(1:NBLOCG)
         BYMAXW(1:NBLOCG) = BYMAX(1:NBLOCG)
         BZMINW(1:NBLOCG) = BZMIN(1:NBLOCG)
         BZMAXW(1:NBLOCG) = BZMAX(1:NBLOCG)

         BXMIN(1:NBLOCG) = 0.0
         BXMAX(1:NBLOCG) = 0.0
         BYMIN(1:NBLOCG) = 0.0
         BYMAX(1:NBLOCG) = 0.0
         BZMIN(1:NBLOCG) = 0.0
         BZMAX(1:NBLOCG) = 0.0

         CALL MPI_ALLREDUCE(BXMINW,BXMIN,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(BXMAXW,BXMAX,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(BYMINW,BYMIN,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(BYMAXW,BYMAX,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(BZMINW,BZMIN,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(BZMAXW,BZMAX,NBLOCG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

         DEALLOCATE(BXMINW,BXMAXW,BYMINW,BYMAXW,BZMINW,BZMAXW)

      ENDIF

C *********************************************************************

      IF(AREF < 1.E-12 ) THEN
         WRITE(*,*) 'AREF must be calculated in post-processing'
         WRITE(*,*) 'Exiting ...'
         STOP
C         WRITE(*,*) 'SURFACE AREA USED IN CALCULATIONS = ',AREF
C         WRITE(*,*)
      ENDIF  ! AREF < 1.E-12 

      IF (IPRINT-1 <= ICYCLE) THEN

C ... BLOCK LOOP 4 BEGINS. PRINT QUANTITIES RELATED TO THE GEOMETRY
       
      WRITE(46,*)
      WRITE(46,*) ' ICYCLE = ',ICYCLE
      DO 9100 N = 1,NBLOCK
      NGL     = NPROCE(1+N,IPRO)
      IF(IALL == 0 .OR. IGRID(NGL,1) /= 0 .OR. IGRID(NGL,2) /= 0) THEN
      DO 9000 M = 1,MGRID(N)

      NGL     = NPROCE(1+N,IPRO) !Global block number
      IG2     = IG(M,N)
      KU1     = IT(M,N)
      KU2     = IL(M,N)
      KU3     = IK(M,N)
      KU4     = IMAX(M,N)
      KU5     = JMAX(M,N)
      KU6     = KMAX(M,N)
      IGN     = IG2 - 1 + NTOT(M,N)
      IF(MULPHL) THEN ! Antakee mersu, kokeilee
      PRO(IG2:IGN)%ALPO(1) = A1XA(IG2:IGN)
      ENDIF

      IF(KU2  <= 3) THEN
      WRITE(46,*)
      WRITE(46,*)'FROM ',CALSUB,', IN GRIVAL(4):'
      ENDIF
      CALL PRINYD(46,XCOC,XCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,YCOC,YCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,ZCOC,ZCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      IF(COORL) THEN
      IF(KU2 <= 3) THEN 
      WRITE(46,*)
      WRITE(46,*)'FROM ',CALSUB,', IN GRIVAL(4): ORIGINAL CORNER-POINTS'
      ENDIF
      CALL PRINYD(46,XORC,XGRI(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,YORC,YGRI(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,ZORC,ZGRI(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      ENDIF
      CALL PRINYD(46,XCC ,XC(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,YCC ,YC(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,ZCC ,ZC(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)

      CALL PRINYD(46,'XLE2    ',XLE2(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,'YLE2    ',YLE2(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,'ZLE2    ',ZLE2(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,'XLE3    ',XLE3(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,'YLE3    ',YLE3(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,'ZLE3    ',ZLE3(IG2),
     &            KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)

      CALL PRINYS(46,VOLC,VOL(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A1C ,A1(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A2C ,A2(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A3C ,A3(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A1XC,A1XA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A1YC,A1YA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A1ZC,A1ZA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A2XC,A2XA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A2YC,A2YA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A2ZC,A2ZA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A3XC,A3XA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A3YC,A3YA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,A3ZC,A3ZA(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,D1C ,D1(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,D2C ,D2(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,D3C ,D3(IG2)  ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,DISTWC,DISTW(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,VTRANC,VTRAN(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
CESA      CALL PRINYS(46,DISTWC,LOCDIS(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
c      IF(OMEGA(NPROCE(1+N,IPRO)) /= 0) THEN
      CALL PRINYS(46,UROTC,UROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,VROTC,VROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,WROTC,WROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,XSKC,XFC(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,YSKC,YFC(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,ZSKC,ZFC(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
c      ENDIF

9000  CONTINUE
      ENDIF
9100  CONTINUE

C ... END OF BLOCK LOOP 4
C ********************************************************************
      ENDIF  ! (IPRINT == ICYCLE) 

C ... BLOCK LOOP 5 BEGINS. PRINT QUANTITIES RELATED TO THE GEOMETRY
C ... IF GRID is Moving

      IF(TIMEL) THEN
      WRITE(46,*)
      WRITE(46,*) '  IN GRIVAL: TIME =',T   
      WRITE(46,*)
      DO 9150 N = 1,NBLOCK
      NGL     = NPROCE(1+N,IPRO)
      IF(IGRID(NGL,1) /= 0 .OR. IGRID(NGL,2) /= 0) THEN
      DO 9050 M = 1,MGRID(N)

      NGL     = NPROCE(1+N,IPRO) !Global block number
      IG2     = IG(M,N)
      KU1     = IT(M,N)
      KU2     = IL(M,N)
      KU3     = IK(M,N)
      KU4     = IMAX(M,N)
      KU5     = JMAX(M,N)
      KU6     = KMAX(M,N)

      CALL PRINYS(46,VOLC,VOL(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      IF(KU2 <= 3) THEN
      WRITE(46,*)
      WRITE(46,*) 'FROM ',CALSUB,', IN GRIVAL(TIMEL):'
      ENDIF
      CALL PRINYD(46,XCOC,XCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,YCOC,YCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYD(46,ZCOC,ZCO(IG2) ,KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      IF(OMEGA(NPROCE(1+N,IPRO)) /= 0 .OR. COORL) THEN ! Rotating velocities
      CALL PRINYS(46,UROTC,UROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,VROTC,VROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,WROTC,WROT(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      CALL PRINYS(46,UROTBC, DN(IG2),KU1,KU2,0,KU3,KU4,KU5,KU6,NGL,M)
      ENDIF

9050  CONTINUE
      ENDIF
9150  CONTINUE
      ENDIF  ! IF (TIMEL)
C ... END OF BLOCK LOOP 5

C ********************************************************************
      IF(IREPEA(3) <= NREPEA(3)) THEN
      WRITE(45,*)
      VOLTOT = 0.
      DO 9200 N = 1,NBLOCK
         WRITE(45,111) N,VOLN(N)
         VOLTOT  = VOLTOT + VOLN(N) ! TOTAL VOLUME
 9200 CONTINUE
      WRITE(45,112) VOLTOT
 111  FORMAT(' Volume of block ',I4,' is ',E13.5)
 112  FORMAT(' ====================================='/
     &       ' Total volume    ',4X,' is ',E13.5)
      ENDIF ! IREPEA(3) <= 10

      RETURN
      END SUBROUTINE GRIVAL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INPOST(DREL,TIME2)

      USE CHARACTERS

      USE TYPE_ARRAYS

      USE INTEGERS,    ONLY : MAXB,MAXEB,MAXSB,IB,MGM,IPRO,MAXW,
     &    MAXSS,MBPRO

      USE MAIN_ARRAYS, ONLY : TEMP,U,V,W,E,C,EPS2,VIST,P,PDIFF,RK,
     &    REPS,FI,BIJ,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     &    XC,YC,ZC,ROLE2,RMLE2,RNLE2,RWLE2,ELE2,RKLE2,EPSLE2,
     &    RKLE3,EPSLE3,ROLE3,RMLE3,RNLE3,RWLE3,ELE3,FILE2,FILE3,
     &    ZZZ,F1R,F1RM,JLOC,JTRA,APP,IG,RO,RM,RN,RW,NTOT,IGRID,
     &    XCO,YCO,ZCO,XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,IR,IQ,JF,IC,
     &    VIS,CP,CH,DRDP,DRDH,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     &    KBOT,KTOP,IDI1,IDI2,IDI3,NPATCH,ICON,RCON,UWALL,VWALL,
     &    WWALL,IHF,OHMI,DM,DN,DW,DE,S11,VOL,A1,A2,A3,D1,D2,D3,IDER,
     &    NPROCE,DDEPS,PTUR,UROT,VROT,WROT,IT,IL,IK,OMEGA,JSTATE,VTRAN,
     &    BLKS,VAR,F1RN,F1RK,PLE2,PLE3

      USE NS3CO,       ONLY : IN,JN,KN,REANEW,ITURB,STRESL,TIMEL,
     &    PARALLEL,GX,GY,GZ,GROUND,FRSDEN,NBLOCK,FRSPRE,
     &    ISTRES,GAMMA,PR,PRT,RGAS,VISU0,EXPSU,TSU0,E0REF,FRSTEM,
     &    FRSVEL,FRSSSP,IPRESC,PSEUCO,NSCAL,JRIMP,T0,TLOLIM,
     &    TUPLIM,IOLD,COORL,T0REF,IPRINT,ICYCLE,IDRXX,IFSBC,MULPHL,
     &    TWO_FLUIDL

      IMPLICIT NONE

      INTEGER :: N,IG2,IR2,IQ2,IF2,IC2,IE2,IQQ2,M,I,L,NS,IG1,NGL,
     &           IPHASE,NPHASE,IGN,ICHARP

      REAL    :: DREL

      LOGICAL :: TIME2   

C *** INITIALIZATION OF STATE VARIABLES AFTER RESTART ******************

C ... compute the free surface (in IFSBC=3 only)

      IF(IFSBC == 3) THEN
        CALL FRE(1,1,NBLOCK,1.,1)
      ENDIF

      IF(REANEW .AND. ITURB >= 10 .AND. ITURB <= 19) THEN 
      WRITE(45,*) 'Starting INPOST'

C ... MAKE A FIRST GUESS FOR THE REYNOLDS STRESSES WITH ALGEBRAIC MODELS

      DO 8400 N = 1,NBLOCK   !  BLOCK LOOP 10 BEGINS
      
      NGL     = NPROCE(1+N,IPRO) !Global block number
      IG2     = IG(1,N)
      IR2     = IR(1,N)
      IQ2     = IQ(1,N)
      IF2     = JF(1,N)
      IC2     = IC(1,N)
      IE2     = 1
      IF(ISTRES >= 1) IE2 = IG2
      IQQ2    = 1
      IF(STRESL) IQQ2 = IG2

C ... Firstly the equation of state

      CALL STATIM(1,N,2,.FALSE.)

C ... Then vorticity components (= DM, DN, DW), strain (= DE) and
C ... the Reynolds stresses (=S11)

      CALL VORFUN(OHMI(IG2),DM(IG2),DN(IG2),DW(IG2),DE(IG2),
     1 S11(IQQ2,1),S11(IQQ2,2),S11(IQQ2,3),S11(IQQ2,4),
     2 S11(IQQ2,5),S11(IQQ2,6),PTUR(IG2),U(IG2),V(IG2),
     3 W(IG2),RK(IR2),REPS(IR2),DDEPS(IR2),EPS2(IG2),VIST(IG2),VIS(IG2),
     4 A1(IG2),A2(IG2),A3(IG2),VOL(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     5 A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     6 D1(IG2),D2(IG2),D3(IG2),XC(IG2),YC(IG2),ZC(IG2),IMAX(1,N),
     7 JMAX(1,N),KMAX(1,N),IN,JN,KN,ITURB,1,STRESL,ZZZ,VTRAN(IR2),MAXW,
     8 IDER(N),BIJ(IE2,1),MAXEB,ISTRES)
 8400 CONTINUE

      CALL MIRR(S11(1,1),S11(1,2),S11(1,3),S11(1,4),
     +     S11(1,5),S11(1,6),A1XA,A1YA,A1ZA,
     +     A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,1,NBLOCK,MAXSS)
      DO 900 NS = 1,6
         CALL CONEC(S11(1,NS),F1R,F1RM,JLOC,JTRA,APP,
     +   1,1,NBLOCK,0,0)
 900  CONTINUE
      ENDIF  ! reynolds stress guest

      IF(IOLD > 0) THEN

      DO 8300 N = 1,NBLOCK   !  BLOCK LOOP 10 BEGINS (for all blocks...)
c      DO 8300 M = 1,MGRID(N)! is this needed? PPR 29.10.00
      DO 8300 M = 1,1 

      NGL     = NPROCE(N+1,IPRO)
      IG2     = IG(M,N)
      IR2     = IR(M,N)
      IQ2     = IQ(M,N)
      IF2     = JF(1,N)
      IC2     = IC(1,N)
      NPHASE  = BLKS(NGL)%NPHASE
         
      CALL STATIM(M,N,2,.FALSE.)

C ... CALCULATE VELOCITIES AND INTERNAL ENERGY FOR THE BOUNDARY BLOCKS

      IF(TIMEL .AND. ABS(DREL - 1.) >=  1.E-4) THEN

         IF(IPRO == 1) WRITE(*,313) DREL
         WRITE(45,313) DREL
 313     FORMAT(3X,
     & 'Time step has been changed. Old time levels are changed by ',
     &        F10.5)

         DO I=1,NTOT(M,N)
            L = IG2+I-1
            ROLE2(L) = DREL*ROLE2(L) + (1.-DREL)*RO(L)
            RMLE2(L) = DREL*RMLE2(L) + (1.-DREL)*RM(L)
            RNLE2(L) = DREL*RNLE2(L) + (1.-DREL)*RN(L)
            RWLE2(L) = DREL*RWLE2(L) + (1.-DREL)*RW(L)
            ELE2(L)  = DREL*ELE2(L)  + (1.-DREL)*E(L)
         ENDDO

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         DO I=1,NTOT(M,N)
            L = IG2+I-1
            RKLE2(L) = DREL*RKLE2(L) + (1.-DREL)*RK(L)
            EPSLE2(L)= DREL*EPSLE2(L)+ (1.-DREL)*REPS(L)
         ENDDO
         ENDIF

         IF(NPHASE > 1) THEN 
         DO IPHASE = 1,NPHASE
         DO I=1,NTOT(M,N)
            L = IG2+I-1
            VAR(L)%AROLE2(IPHASE) = DREL*VAR(L)%AROLE2(IPHASE) +
     &                          (1.-DREL)*VAR(L)%ARO(IPHASE)
            VAR(L)%ARELE2(IPHASE) = DREL*VAR(L)%ARELE2(IPHASE) +
     &                          (1.-DREL)*VAR(L)%ARE(IPHASE)
         ENDDO
         ENDDO
         IF(TWO_FLUIDL) THEN
         DO IPHASE = 1,NPHASE
         DO I=1,NTOT(M,N)
            L = IG2+I-1
            VAR(L)%ARMLE2(IPHASE) = DREL*VAR(L)%ARMLE2(IPHASE) +
     &                          (1.-DREL)*VAR(L)%ARM(IPHASE)
            VAR(L)%ARNLE2(IPHASE) = DREL*VAR(L)%ARNLE2(IPHASE) +
     &                          (1.-DREL)*VAR(L)%ARN(IPHASE)
            VAR(L)%ARWLE2(IPHASE) = DREL*VAR(L)%ARWLE2(IPHASE) +
     &                          (1.-DREL)*VAR(L)%ARW(IPHASE)
         ENDDO
         ENDDO
         ENDIF ! TWO_FLUIDL
         ENDIF ! NPHASE > 1

         IF (NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19))THEN
            DO NS = 1,NSCAL
               DO I=1,NTOT(M,N)
                  L = IG2+I-1
                  FILE2(L,NS) = DREL*FILE2(L,NS)+(1.-DREL)*FI(L,NS)
               ENDDO
            ENDDO
      ENDIF
      ENDIF ! DREL TIME STEP CHANGE

      IF(TIMEL) THEN

C ... First order approximation for the third time level

         DO I=1,NTOT(M,N)
            L = IG2+I-1
            RKLE3(L) = 2.*RKLE2(L) - RK(L)
            EPSLE3(L)= 2.*EPSLE2(L)- REPS(L)
            ROLE3(L) = 2.*ROLE2(L) - RO(L)
            RMLE3(L) = 2.*RMLE2(L) - RM(L)
            RNLE3(L) = 2.*RNLE2(L) - RN(L)
            RWLE3(L) = 2.*RWLE2(L) - RW(L)
            ELE3(L)  = 2.*ELE2(L)  - E(L)
            PLE3(L)  = PDIFF(L)
            PLE2(L)  = PDIFF(L)
         ENDDO
         IF(NPHASE > 1) THEN 
         DO IPHASE = 1,NPHASE
         DO I=1,NTOT(M,N)
            L = IG2+I-1
            VAR(L)%AROLE3(IPHASE) = 2.*VAR(L)%AROLE2(IPHASE)    -
     &                                VAR(L)%ARO(IPHASE)
            VAR(L)%ARELE3(IPHASE) = 2.*VAR(L)%ARELE2(IPHASE) -
     &                                VAR(L)%ARE(IPHASE)
         ENDDO
         ENDDO ! NPHASE
         IF(TWO_FLUIDL) THEN
         DO IPHASE = 1,NPHASE
         DO I=1,NTOT(M,N)
            L = IG2+I-1
            VAR(L)%ARMLE3(IPHASE) = 2.*VAR(L)%ARMLE2(IPHASE)    -
     &                                VAR(L)%ARM(IPHASE)
            VAR(L)%ARNLE3(IPHASE) = 2.*VAR(L)%ARNLE2(IPHASE)    -
     &                                VAR(L)%ARN(IPHASE)
            VAR(L)%ARWLE3(IPHASE) = 2.*VAR(L)%ARWLE2(IPHASE)    -
     &                                VAR(L)%ARW(IPHASE)
         ENDDO
         ENDDO
         ENDIF ! TWO_FLUIDL

         ENDIF ! NPHASE > 1
         IF (NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19))THEN
            DO NS = 1,NSCAL
               IF(ABS(DREL - 1.) >=  1.E-4) THEN
                  DO I=1,NTOT(M,N)
                     L = IG2+I-1
                     FILE2(L,NS) = DREL*FILE2(L,NS)+(1.-DREL)*FI(L,NS)
                  ENDDO
               ENDIF
               DO I=1,NTOT(M,N)
                  L = IG2+I-1
                  FILE3(L,NS) = 2.*FILE2(L,NS) - FI(L,NS)
               ENDDO
            ENDDO
         ENDIF
      ENDIF ! TIMEL

      
 8300 CONTINUE                        ! END OF BLOCK LOOP 10
       
      IF(TIMEL .AND. .NOT. TIME2) THEN 
      CALL STARTI(NTOT,IG,NBLOCK,TIMEL,COORL,ITURB,NSCAL,MGM,
     +  IGRID(:,1),RO,RM,RN,RW,E,RK,REPS,FI,XCO,YCO,ZCO,
     +  ROLE2,RMLE2,RNLE2,RWLE2,ELE2,RKLE2,EPSLE2,FILE2,
     +  ROLE3,RMLE3,RNLE3,RWLE3,ELE3,RKLE3,EPSLE3,FILE3,
     +  XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,BLKS,VAR,MAXSB,IPRO,NPROCE,MBPRO,
     +  TWO_FLUIDL)
      ENDIF

      ENDIF !What is this? REANEW ?
 
      IF (IPRINT == (ICYCLE+1)) THEN
      WRITE(3,*) ' INITIAL VALUES'

      DO 3000 N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         IF(OMEGA(NPROCE(1+N,IPRO)) /= 0.) THEN
         CALL PRINYS(3,UROTC,UROT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,VROTC,VROT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,WROTC,WROT(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         CALL PRINYS(3,ROC,RO(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3, PC, P(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,TEMPC,TEMP(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,RM(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RNC,RN(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RWC,RW(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3, EC, E(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3, CC, C(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(TIMEL) THEN
         WRITE(3,*)
         WRITE(3,*)' SECOND TIMEL-LEVEL'
         CALL PRINYS(3,ROC,ROLE2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RMC,RMLE2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RNC,RNLE2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,RWC,RWLE2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(MULPHL) THEN
               IGN  = IG1 + NTOT(1,N) - 1
               DO IPHASE = 1,BLKS(NGL)%NPHASE
               ICHARP = BLKS(NGL)%ICHAR(IPHASE)
               F1RK(IG1:IGN) = VAR(IG1:IGN)%ARO(IPHASE)
               F1RM(IG1:IGN) = VAR(IG1:IGN)%AROLE2(IPHASE)
               F1RN(IG1:IGN) = VAR(IG1:IGN)%AROLE3(IPHASE)
               CALL PRINYS(3,CHAR_VAR(16,ICHARP),F1RK(IG1),IT(1,N),
     +            IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' AROLE2 ',F1RM(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' AROLE3 ',F1RN(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               ENDDO ! IPHASE
               IF(TWO_FLUIDL) THEN
               DO IPHASE = 1,BLKS(NGL)%NPHASE
               F1RK(IG1:IGN) = VAR(IG1:IGN)%ARM(IPHASE)
               F1RM(IG1:IGN) = VAR(IG1:IGN)%ARMLE2(IPHASE)
               F1RN(IG1:IGN) = VAR(IG1:IGN)%ARMLE3(IPHASE)
               CALL PRINYS(3,'  ARM   ',F1RK(IG1),IT(1,N),
     +            IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARMLE2 ',F1RM(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARMLE3 ',F1RN(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

               F1RK(IG1:IGN) = VAR(IG1:IGN)%ARN(IPHASE)
               F1RM(IG1:IGN) = VAR(IG1:IGN)%ARNLE2(IPHASE)
               F1RN(IG1:IGN) = VAR(IG1:IGN)%ARNLE3(IPHASE)
               CALL PRINYS(3,'  ARN   ',F1RK(IG1),IT(1,N),
     +            IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARNLE2 ',F1RM(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARNLE3 ',F1RN(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

               F1RK(IG1:IGN) = VAR(IG1:IGN)%ARW(IPHASE)
               F1RM(IG1:IGN) = VAR(IG1:IGN)%ARWLE2(IPHASE)
               F1RN(IG1:IGN) = VAR(IG1:IGN)%ARWLE3(IPHASE)
               CALL PRINYS(3,'  ARW   ',F1RK(IG1),IT(1,N),
     +            IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARWLE2 ',F1RM(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               CALL PRINYS(3,' ARWLE3 ',F1RN(IG1),IT(1,N),IL(1,N),0,
     +              IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               ENDDO ! IPHASE
               ENDIF ! TWO_FLUIDL
          ENDIF
           
         ENDIF
         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL PRINYS(3,RKC,RK(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(ITURB /= 8 .OR. ITURB /= 9)
     +   CALL PRINYS(3,REPSC,REPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         IF(ITURB == 9)
     +   CALL PRINYS(3,RNUTC,REPS(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,EPS2C,EPS2(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         IF(ITURB >= 10 .AND. ITURB <= 19) THEN
         DO 8171 NS = 1,6
c         IG3 = IG1 + (NS-1)*MAXSS
         CALL PRINYS(3,UUC(NS),S11(IG1,NS),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
 8171    CONTINUE
         ENDIF
 3000 CONTINUE
      ENDIF

C ... ONLY MASTER WRITES THE ITERATION HISTORY ON THE SCREEN
        
      IF(IPRO == 1 .OR. .NOT. PARALLEL) THEN
      WRITE(*,*)
      IF(IDRXX == 1) THEN
         WRITE(6,111) 'DROMAX','******'
      ELSEIF(IDRXX == 2) THEN
         WRITE(6,111) 'DUMAX ','******'
      ELSEIF(IDRXX == 3) THEN
         WRITE(6,111) 'DVMAX ','******'
      ELSEIF(IDRXX == 4) THEN
         WRITE(6,111) 'DWMAX ','******'
      ELSEIF(IDRXX == 5) THEN
         WRITE(6,111) 'DTMAX ','***** '
      ELSEIF(IDRXX == 6) THEN
         WRITE(6,111) 'DPMAX ','***** '
      ELSEIF(IDRXX == 7) THEN
         WRITE(6,111) 'DRKMAX','******'
      ELSEIF(IDRXX == 8) THEN
         WRITE(6,111) 'DEPMAX','******'
      ELSEIF(IDRXX >= 9) THEN
         WRITE(6,112) IDRXX-8,'******'
      ENDIF
      ENDIF
 111  FORMAT('  ICYCLE      ',A6,'       N     I   J   K     NSPSC',
     +       '       CD           CL'/
     +       '  ******      ',A6,'      ***   *** *** ***    *****',
     +       '    ********     ********')
 112  FORMAT('  ICYCLE      DF',I1,'MAX       N     I   J   K     NSPSC'
     +      ,'       CD           CL'/
     +       '  ******      ',A6,'      ***   *** *** ***    *****',
     +       '    ********     ********')
C***********************************************************************

      RETURN
      END SUBROUTINE INPOST
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NEWROT

      USE MPI

      USE INTEGERS,    ONLY : MAXB,MBPRO,IPRO,IREPEA,NREPEA,NBCS

      USE MAIN_ARRAYS, ONLY : IGRID,XCO,YCO,ZCO,NPROCE,HFLUX,
     &     NTOT,OMEGA,XGRI,YGRI,ZGRI,XORI,YORI,ZORI,IG,IMAX,JMAX,KMAX,
     &     OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,AMPL,OSCLS,
     &     IC,ICON,MGRID,
     &     XC,YC,ZC,IHF,NPATCH,WAVES,XCOL,CXB,CYB,CZB,CMXB,CMYB,CMZB,
     &     APATCH,XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3   

      USE TYPE_ARRAYS, ONLY : SIXDOF

      USE FLIGHT,      ONLY : XCGI,YCGI,ZCGI,XCG,YCG,
     &     ZCG,PSIR,THETAR,PHIR,PSIM,NGRIFL,OSKU,TRMODE,
     &     PSIRI,THETARI,PHIRI

      USE BLADE_VARIABLES, ONLY : ZETAH,TDBL,DTB

      USE NS3CO,       ONLY : IN,JN,KN,ROTANG,OMEX,OMEY,OMEZ,
     &     CENX,CENY,CENZ,GVEX,GVEY,GVEZ,NBLOCK,TIMEL,T,IOLD,
     &     IPRINT,ICYCLE,ITIMES,NBLOCG,ROTAT,IPRIFS,LDISTL,ALPHA,
     &     BETA,PARALLEL,DT,PARALLEL,XMOM,YMOM,ZMOM,TROT

      USE CHARACTERS,  ONLY : BOUNDF

      USE CONSTANTS, ONLY : PII

      IMPLICIT NONE

      INTEGER :: N,NTOT1,IG1,IC1,M,NBG,III,IS,IGR,IERR,
     &           ICASE,ICASEX,I,L

      CHARACTER(LEN=80) :: RIVI

C ... ROTATIONAL AXIS AND THE PLACE AND VELOCITY. Could be needed 
C ... different axis for different blocks?

      OMEX = 1.                 ! rotational axis
      OMEY = 0.
      OMEZ = 0.
      CENX = 0.                 ! center of the rotation
      CENY = 0.
      CENZ = 0.
c      GVEX = 0.                 ! grid velocity
c      GVEY = 0.
c      GVEZ = 0.
      
C***********************************************************************
C ... Update the position of grid blocks
C***********************************************************************

C ... CHANGE THE GRID ACCORDING TO ROTATING ANGLE
C ... Note that in parallel calculation this is in global mode

      DO N = 1,NBLOCG           ! Block loop '0' begins (mersu)
         IF(IGRID(N,1) >= 1.AND.IGRID(N,1) <= 10.AND.OMEGA(N) /= 0.)THEN
            IF(TIMEL .AND. IGRID(N,1) < 10) THEN ! ADVANCE ROTANG 

               IF(T > TROT .AND. OSCLS(N) == 1) THEN ! Sinus shape oscillation 
                  ROTANG(N) = ROTAT(N) + SIN(OMEGA(N)*(T+DT-TROT))*
     &                       (AMPL(N)/2.*PII/180.)
               ELSE
                  ROTANG(N) = ROTAT(N) + OMEGA(N)*(T-TROT) 
               ENDIF

            ELSEIF(TIMEL .AND. IGRID(N,1) == 10) THEN ! SHIP PROPULSION
               ROTANG(N) = -PSIM(IGRID(N,3))+ OMEGA(N)*DT   
               PHIR(IGRID(N,3)) = -ROTANG(N)
            ELSEIF(.NOT. TIMEL .AND. IOLD > 0) THEN  
               WRITE(*,*) 'KUINKA PYGRITETAAN JOS ON EI AJAN SUHT TARK'
            ENDIF
         ENDIF
      ENDDO
        
ccc      DO N = 1,NBLOCK           ! BLOCK LOOP 1 BEGINS (works without MPI)?

      ICASE = 0

      DO NBG = 1,NBLOCG           ! BLOCK LOOP 1 BEGINS (works without MPI)?

ccc        NBG  = NPROCE(1+N,IPRO)
      
         IF(IGRID(NBG,1) >= 11 .AND. IGRID(NBG,1) < 19) THEN ! 6-dof
            ICASE  = MAX(ICASE,1)
         ELSEIF (IGRID(NBG,1) == 22) THEN ! 2-dof movement for a ship
            ICASE  = MAX(ICASE,1)
         ELSEIF(IGRID(NBG,1) == 37 .OR. IGRID(NBG,1) == 38) THEN
            ICASE  = MAX(ICASE,2)
         ELSEIF(IGRID(NBG,1) == 39) THEN
            ICASE  = MAX(ICASE,3)
         ELSEIF(IGRID(NBG,1) >= 31 .AND. IGRID(NBG,1) < 37) THEN
            ICASE  = MAX(ICASE,1)
         ENDIF

      ENDDO


ccc      IF (PARALLEL) THEN         ! Knowledge about ICASE is broadcasted
ccc            CALL MPI_REDUCE(ICASE,ICASEX,1,MPI_INTEGER,MPI_MAX,0,
ccc     +      MPI_COMM_WORLD,IERR)
ccc            IF (IPRO == 1) ICASE = ICASEX
ccc            CALL MPI_BCAST(ICASE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
ccc      ENDIF ! PARALLEL

      WRITE(45,*)

C ... A new position of the flying object is calculated in the 
C ... PARTICLE-routine. This is done when TIMEL is true 
      
      IF(TIMEL .AND. IPRO == 1) THEN

      SELECT CASE(ICASE)

      CASE(1)

         DO IGR = 1,NGRIFL
C ... Calling routine of the trajectory 'IGR'
            IF(TRMODE(IGR) /= 37 .OR. TRMODE(IGR) /= 38) THEN
           CALL PARTICLE(OSKU,XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     &              THETAR(IGR),PHIR(IGR),PSIM(IGR),IGR,TDBL-DTB,DTB,
     &              TRMODE(IGR))
            ENDIF
         ENDDO

      CASE(2)

cc         DO IGR = 1,NGRIFL
C ... Does not work yet
c         write(805+igr,*) OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
c     +   OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ
cc            CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),ZCGI(IGR),
cc     +           XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
cc     +           PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,TRMODE(IGR))
cc         ENDDO         

C ... Call blade_cycles         
c         IF (T/DT >= 0.33333) THEN
            CALL BLADE_CYCLES(ICYCLE,TDBL,NGRIFL) ! Blade loop inside
            IREPEA(6) = IREPEA(6) + 1 ! For time-step warning
c         ENDIF

      END SELECT

      ENDIF                     ! IF(TIMEL .AND. IPRO == 1

C ... In a parallel mode send the required information

      IF (PARALLEL .AND. TIMEL .AND. NGRIFL > 0) THEN
         CALL MPI_BCAST(OSKU%CX,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CY,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CZ,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CMX,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CMY,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CMZ,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR) 
         CALL MPI_BCAST(YCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETAR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIM,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZETAH,NGRIFL,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
C tarvitaanko n\E4it??????
c         CALL MPI_BCAST(BETA,NGRIFL,MPI_REAL8,0,
c     +        MPI_COMM_WORLD,IERR)
      ENDIF
         
C***********************************************************************
C ... Apply grid movements
C***********************************************************************

      IS = 1
       
      DO 2000 N = 1,NBLOCK      ! BLOCK LOOP 2 BEGINS
          
         IG1 = IG(1,N)
         IC1 = IC(1,N) 
         M   = 1
         NBG = NPROCE(1+N,IPRO)

C ... Rotating blocks for turbomachinery

      IF(IGRID(NBG,1) >= 1.AND.IGRID(NBG,1) <= 9.AND.ROTANG(NBG) /= 0.)
     +   THEN                                     ! Grid rotation
         NTOT1   = NTOT(1,N)                      ! IGRID(N,1)=1...10
         CALL ROTGRI(XCO(IG1),YCO(IG1),ZCO(IG1),XGRI(IG1),
     +        YGRI(IG1),ZGRI(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +        ROTANG(NBG),IN,JN,KN,OMEGAX(NBG),OMEGAY(NBG),OMEGAZ(NBG),
     +        CENAX(NBG),CENAY(NBG),CENAZ(NBG))
         WRITE(45,333) N,1000.*T,ROTANG(NBG)*180./3.1415926
 333     FORMAT(' NEW GRID POSITION FOR BLOCK',I4,' ROTANG(',E10.3,
     +        ' ms) = ',F10.4,' degrees')

C ... Ship propel

      ELSEIF(IGRID(NBG,1) == 10 .AND. TIMEL) THEN ! Propel in ship case
         NTOT1   = NTOT(1,N) 
         IGR = IGRID(NBG,3) 
C ... Set location(s) and attitude(s) of the propel(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the propel(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

C ... Flying 6-DoF blocks

      ELSEIF(IGRID(NBG,1) >= 11 .AND. 
     &       IGRID(NBG,1) < 20 .AND. TIMEL) THEN ! TRAJECTORY Calculation
         NTOT1   = NTOT(1,N) 
         IGR = IGRID(NBG,3)                            ! Particle number
c         CALL FLPATH(OSKU,XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),   ! turha?
c     &        THETAR(IGR),PHIR(IGR),ICYCLE,IPRINT,NGRIFL,T,DT,TRMODE)

C ... A rigid grid block movement using IREPEA(4) for 
C ... print check. The new position of the flying object is calculated 
C ... in the DO-loop above (DO IGR = 11,10+NGRIFL).
C ... Update angle(s) of control surface(s) to (X,Y,Z)ORI-table
C ... A flap movement by the FLAPPOS program  

C ... Set angle(s) of control surface(s) to (X,Y,Z) ORI-table
         CALL FLAPPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &        IN,JN,KN,NBG,TDBL,DTB,NTOT(1,N))
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)
       
C ... Floating 6-DoF blocks

      ELSEIF(IGRID(NBG,1) == 22 .AND. TIMEL
     &    .OR. IGRID(NBG,1) >= 24 .AND. IGRID(NBG,1) < 30 .AND. TIMEL) 
     &        THEN              ! Ship trajectories
         WRITE(45,*) 'Ships are not yet moving. Exiting...'
         WRITE(13,*) 'Ships are not yet moving. Exiting...'
         WRITE(*,*)  'Ships are not yet moving. Exiting...'
         STOP

C ... Blocks for rotors

      ELSEIF(IGRID(NBG,1) >= 31 .AND. 
     &       IGRID(NBG,1) <= 39 .AND. TIMEL) THEN    ! Helicopter rotors
    
         NTOT1 = NTOT(1,N)
         IGR   = IGRID(NBG,3)                        ! Particle number
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),PHIRI(IGR),
     &        THETARI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),PHIR(IGR),
     &        THETAR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) /= 0) THEN ! tama on vahan outo testi, joo on

         WRITE(45,*) 'Nblock = ',NBG,' IGRID(NBG,1) = ',IGRID(NBG,1),
     &              'ROTANG = ',ROTANG(NBG)
         WRITE(45,*) 'Grid was not changed ...'
c        STOP

      ENDIF ! Grid rotation

C ... Update the cell center point coordinates (ESa 4.2.1997)
          
      CALL GRIDEX(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),ICON(IC1),NPATCH(N),XCOL(1,1,IS),IN,JN,KN)
      IS  = IS + NPATCH(N)         

 2000 CONTINUE                  ! END OF BLOCK LOOP 2

C***********************************************************************
C ... Next perform grid deformation
C***********************************************************************
  
      IS      = 1

      DO 3000 N = 1,NBLOCK      ! BLOCK LOOP 3 BEGINS

         IG1 = IG(1,N)
         IC1 = IC(1,N) 
         M   = 1
         NBG = NPROCE(1+N,IPRO)

C ... A free-surface movement (types 5...8) using IREPEA(1) for print check

      IF(IGRID(NBG,2) >= 5 .AND. IGRID(NBG,2) < 9) THEN  ! A free surface

         CALL WSHAPE(XGRI(IG1),YGRI(IG1),ZGRI(IG1),  
     &        XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     &        IHF(1,M),ICON(IC1),NPATCH(N),WAVES,ICYCLE,
     &        IPRINT,IREPEA,NREPEA,NBG)

      IF (TIMEL) THEN !XLE tables are updated in time-accurate simulation
         NTOT1   = NTOT(1,N) 
         DO L  = 1,NTOT1
            I       = L + IG1 - 1
            XLE2(I) = XCO(I)
            YLE2(I) = YCO(I)
            ZLE2(I) = ZCO(I)
            XLE3(I) = XLE2(I)
            YLE3(I) = YLE2(I)
            ZLE3(I) = ZLE2(I)
         ENDDO
      ENDIF ! TIMEL


      ELSEIF(IGRID(NBG,2) == 20) THEN  ! Does not work yet for blades

         IGR = IGRID(NBG,3)                       ! Number of movement type

         CALL BLADE_POS(OSKU,XCGI(IGR),YCGI(IGR),ZCGI(IGR),XCG(IGR),
     &        YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,NGRIFL) ! Parameters nonsense

         CALL GRIDEF(XCO(IG1),YCO(IG1),ZCO(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     &        PHIR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4))
     
C ... Deform the Hornet grid

      ELSEIF(IGRID(NBG,2) == 21) THEN                           ! ylarajaa
         WRITE(*,*) 'Nblock = ',NBG,' IGRID(NBG,2)=',IGRID(NBG,2) ! voitaneen 
         WRITE(*,*) 'Deformation is not yet possible'             ! nostaa
         WRITE(*,*) 'Exiting ...'                                 ! joskus
         STOP

C     CALL CHAGRI   ! Simo Malmi's program


C ... Add deformation on the helicopter blade

      ELSEIF(IGRID(NBG,2) == 32) THEN ! Deformation of a blade ! ylarajaa
         WRITE(*,*) 'Nblock = ',NBG,' IGRID(NBG,2)=',IGRID(NBG,2) ! voitaneen 
         WRITE(*,*) 'Deformation is not yet possible'             ! nostaa
         WRITE(*,*) 'Exiting ...'                                 ! joskus
         STOP

C     CALL BLADE_DEF   ! Antero Miettinen's program

      ENDIF ! A free surface

C ... Update the cell center point coordinates (ESa 4.2.1997)
          
      CALL GRIDEX(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),ICON(IC1),NPATCH(N),XCOL(1,1,IS),IN,JN,KN)
      IS  = IS + NPATCH(N)         

 3000 CONTINUE                  ! END OF BLOCK LOOP 3

C ... CONNECT THE CORNER POINT COORDINATES
       
      CALL MIRC2(XCO,YCO,ZCO,1,1,NBLOCK)

      CALL CONECX(XCO,1,1,NBLOCK,1)
      CALL CONECX(YCO,1,1,NBLOCK,1)
      CALL CONECX(ZCO,1,1,NBLOCK,1)

C ... EXTRAPOLATION AND CONNECTION IS NOT NEEDED 

C***********************************************************************

      IF(ICYCLE-1 == IPRINT) WRITE(45,*)
      IF(IREPEA(3) < NREPEA(3)) THEN ! IREPEA(3) for NEWROT
        WRITE(45,*)
        WRITE(45,*) 'GRIVAL CALL FROM NEWROT:'
        WRITE(45,'(1X,24("-"))')
      ELSE IF(IREPEA(3) == NREPEA(3)) THEN
        WRITE(45,*)
        WRITE(45,*) 'GRIVAL CALL FROM NEWROT:'
        RIVI ="(' THIS HAS BEEN REPEATED',I3,' TIMES. SILENCE')"
        WRITE(45,'(1X,24("="))')
        WRITE(45,RIVI) IREPEA(3)
        WRITE(45,'(1X,40("-"))')
      ENDIF

      CALL GRIVAL(1,LDISTL,'NEWROT')

      IREPEA(3) = IREPEA(3) + 1
C***********************************************************************

C ... EXTRAPOLATE CENTER POINTS FOR SLIDING PATCHES

      DO 2200 N = 1,NBLOCK
      NBG = NPROCE(1+N,IPRO)
      IF(IGRID(NBG,1) >= 1 .AND. IGRID(NBG,1) < 9) THEN ! Rotating
      DO 3210 M = 1,MGRID(N)                             
         IG1       = IG(M,N)
         IC1     = IC(M,N)
         CALL SLICOR(XC(IG1),YC(IG1),ZC(IG1),IMAX(M,N),JMAX(M,N),
     +        KMAX(M,N),IN,JN,KN,N,NPATCH(N),ICON(IC1),ROTANG,NPROCE,
     +        IPRO,MBPRO)
 3210 CONTINUE
      ENDIF
2200  CONTINUE

      RETURN
      END SUBROUTINE NEWROT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NEWGRI

      USE FLIGHT,     ONLY : NGRIFL,OSKU,XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,PSIRI,THETARI,PHIRI,TRMODE,PSIM,DRAUGHTI,
     &     DRAUGHT,XCGIS,YCGIS,ZCGIS,DAMPC1,DAMPC2,
     &     ROTA1I,ROTB1I,CONEAI,ROTA1,ROTB1,CONEA,MVSHIP,ROTB1S,
     &     ACTDISK,IGRSLAVE,RTMSP,FDSP,FXTSP,THRUST,TORQUE,IFA,RGML,
     &     IFAKEP
      USE INTEGERS,   ONLY : MAXB,MBPRO,IPRO,IREPEA,NREPEA,MGM,
     &    MAXFS
      USE MAIN_ARRAYS,ONLY : IGRID,XCO,YCO,ZCO,NPROCE,HFLUX,
     &    NTOT,OMEGA,XGRI,YGRI,ZGRI,IG,IMAX,JMAX,KMAX,OMEGAX,
     &    OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,IC,ICON,MGRID,XC,
     &    YC,ZC,IHF,NPATCH,WAVES,IT,IK,IL,WAVEH,
     &    XCOL,JF,XXI,YXI,XETA,YETA,XSP,YSP,ZSP,IBOTGR,JBOTGR,
     &    KBOTGR,ITOPGR,JTOPGR,KTOPGR,XHULL,YHULL,ZHULL,LHULL,
     &    MHULL,MMHUL,ICONH,ICOGH,XI,ETA,ZETA,WH,PDIFF,ZZTOP,
     &    XVL,YVL,ZVL,XORI,YORI,ZORI,SIXDOF,UROTCP,VROTCP,WROTCP
      USE NS3CO,      ONLY : IN,JN,KN,ROTANG,OMEX,OMEY,OMEZ,
     &    CENX,CENY,CENZ,GVEX,GVEY,GVEZ,NBLOCK,TIMEL,T,IOLD,
     &    IPRINT,ICYCLE,NBLOCG,NBLOCG,ROTAT,IPRIFS,IFSBC,LEVEL,
     &    FRSVEL,ICMAX,AREF,CHLREF,REFPRE,GX,GY,GZ,CX,CY,CZ,INDXTS,
     &    CMY,XMOM,YMOM,ZMOM,NGLOBL,FRESUL,IGLOBL,FRSDEN,
     &    ITIMES,DT,PARALLEL,COORL,LDISTL,IBOUT,ICYOLD
      USE CHARACTERS            
      USE TYPE_ARRAYS,ONLY : SIXDOF
      USE TRSINK
      USE GM
      USE BLADE_VARIABLES, ONLY : ICYCLEB,THETACOL,THCYCLON,THCYCLAT,
     &     QTH,QTH_A,QBE_A,TDBL,DTB
      USE CONSTANTS, ONLY : RAD2DEG,DEG2RAD
      USE MPI

      IMPLICIT NONE

      INTEGER :: N, NTOT1, IG1, IC1, M, NBG, iii,
     &           KU1, KU2, KU3, KU4, KU5, KU6,
     &           IS, IF1, NLVL, ISTRID, JSTRID, L, ILE, INDXFS,
     &           JH1, JG1, IGR, ICASE, ICASEX, IERR

      LOGICAL :: GOBACK, MVBLADE, LDISTAN
      CHARACTER(LEN=80) :: RIVI

      REAL, DIMENSION(19) :: ROTB1A, ZCGD, XCGA, YCGA, ZCGA

      REAL :: ZERO
      
      GOBACK  = .TRUE.
      LDISTAN = COORL .AND. LDISTL  ! Grid is changing
      MVBLADE = .FALSE.
      ICASE  = 0

C***********************************************************************
C ... Apply grid movements in steady-state cases (not time-accurate)
C***********************************************************************

      DO N = 1,NBLOCK           ! BLOCK LOOP 1 BEGINS

        NBG  = NPROCE(1+N,IPRO)
        ! CASES-RAKENNE HIEMAN HUONO, KUN VOI OLLA KAKSI (TAI 3) IF-TESTI TOTTA
        ! SAMASSA AJOSSA (ESIM. LAIVA + ACT-LEVY) ... Thus ??
        IF(IGRID(NBG,1) == 22) THEN ! 2-dof ship
           IF(.NOT. ACTDISK) THEN
            GOBACK = .FALSE.
            ICASE  = 1
           ELSEIF(ACTDISK) THEN
              GOBACK = .FALSE.                             
              ICASE  = 6                                   
           ENDIF
         ELSEIF(IGRID(NBG,1) == 37 .OR. IGRID(NBG,1) == 38) THEN
            MVBLADE = .TRUE.
            ICASE  = 2
         ELSEIF(IGRID(NBG,1) == 39) THEN
            MVBLADE = .TRUE.
            ICASE  = 3
         ELSEIF(IGRID(NBG,1) >= 31 .AND. IGRID(NBG,1) < 37) THEN
            ICASE  = 4
         ELSEIF (IGRID(NBG,1) == 21 .OR. IGRID(NBG,1) == 23) THEN
            GOBACK = .FALSE.
         ELSEIF(IGRID(NBG,2) >= 5 .AND. IGRID(NBG,2) <= 9) THEN
            GOBACK = .FALSE.
         ELSEIF(IGRID(NBG,1) >= 40 .AND. IGRID(NBG,1) < 50) THEN
            ICASE  = 5
         ENDIF

      ENDDO

      IF (ACTDISK .AND. .NOT. MVSHIP) THEN ! Onpa hieno CASES-rakenne
         ICASE = 5 ! PolarsternII jee jee !!!!!!
      ENDIF
      
      IF (PARALLEL) THEN         ! Knowledge about ICASE is broadcasted
            CALL MPI_REDUCE(ICASE,ICASEX,1,MPI_INTEGER,MPI_MAX,0,
     +      MPI_COMM_WORLD,IERR)
            IF (IPRO == 1) ICASE = ICASEX
            CALL MPI_BCAST(ICASE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      ENDIF ! PARALLEL

        
      IF(IPRO == 1) THEN ! Calculate new position

      SELECT CASE(ICASE)

         CASE(1) ! A ship in a quasi-steady case

         DO IGR = 1,NGRIFL
C ... Sting forces and moments
            OSKU(IGR)%FXE   = -OSKU(IGR)%CX-OSKU(IGR)%FXT-
     &           GX*OSKU(IGR)%RMASS
            IF(IFAKEP(IGR) == 1) THEN ! Fake propeller is active
               OSKU(IGR)%FXE   = 0*(-OSKU(IGR)%CX)-0*OSKU(IGR)%FXT-
     &              GX*OSKU(IGR)%RMASS
            ENDIF
            OSKU(IGR)%FYE   = -OSKU(IGR)%CY-OSKU(IGR)%FYT-
     &           GY*OSKU(IGR)%RMASS
            OSKU(IGR)%MXE   = -OSKU(IGR)%CMX-OSKU(IGR)%MXT
            OSKU(IGR)%MZE   = -OSKU(IGR)%CMZ-OSKU(IGR)%MZT
C ... Damp forces and moments           
            OSKU(IGR)%FZE   = DAMPC1(IGR)*OSKU(IGR)%RMASS*
     &           SQRT(ABS(GZ/CHLREF))*OSKU(IGR)%VZ
            OSKU(IGR)%MYE   = DAMPC2(IGR)*OSKU(IGR)%IZ*
     &           SQRT(ABS(GZ/CHLREF))*OSKU(IGR)%ROTY
cc            +                   !EQ. 4.49 Laivan kelluvuus
cc     &           OSKU(IGR)%CZ*RGML(IGR)*SIN(PSIR(IGR))! ja vakavuus
C ... Fake propeller generates X-force, Z-force and Y-moment
            IF (IFAKEP(IGR) == 1) THEN
               CALL FAKEPROPELLER(IGR)
            ENDIF
C ... Calls PARTICLE routine
            CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),
     +           ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     +           THETAR(IGR),PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,
     +           TRMODE(IGR),ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
            CALL PARTICLE(OSKU,XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     +           THETAR(IGR),PHIR(IGR),PSIM(IGR),IGR,0.0,DTB,
     +           TRMODE(IGR))
         ENDDO

         CASE(2) ! 2-Dof helicopter blade

            DO IGR = 1,NGRIFL
C ... Does not work yet
c         write(800+igr,*) OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
c     +  OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ

C ... Calls writing routines
               CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),
     +              ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     +              THETAR(IGR),PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,
     +              TRMODE(IGR),ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
               IF(MOD(ICYCLE,10*IBOUT) == 0) THEN
                  CALL TRAJECTORY_CONVERTER(IGR,TRMODE(IGR),TDBL,DTB)
               ENDIF
            ENDDO

C ... Calls BLADE routine
            CALL BLADE_CYCLES(ICYCLE,TDBL,NGRIFL) ! Blade loop inside
            IREPEA(6) = IREPEA(6) + 1          ! For time-step warning

         CASE(3) ! 3-Dof helicopter blade
            WRITE(13,*) ' Blade calculation with 3-Dof model does not ',
     +               'work yet. Change to option 38. Now exiting..'
            STOP '  Stop in blade calculation. See RUN.LOG, exiting...'

cc         CASE(4) ! Reduced 6-Dof helicopter blade

cc            WRITE(13,*) ' Blade calculation with 6-Dof model does not',
cc     +              ' work yet. Change to option 38. Now exiting..'
cc            STOP '  Stop in blade calculation. See RUN.LOG, exiting...'

         CASE(5) ! Actuator disk
c            IF (ICYCLE > ICYOLD) THEN
               DO IGR = 1,NGRIFL
C     ... Calls writing routines
                  CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),
     +                 ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),
     +                 PSIR(IGR),THETAR(IGR),PHIR(IGR),PSIM(IGR),
     +                 OSKU,IGR,TIMEL,TRMODE(IGR),
     +                 ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
                  IF (PARALLEL .AND. IFA(IGR) == 12 .OR. ! Patria's ACT disk 
     +                 PARALLEL .AND. IFA(IGR) == 13) THEN
                     QTH(IGR)   = -1.E10
                     QTH_A(IGR) = -1.E10
                     QBE_A(IGR) = -1.E10
                  ENDIF
               ENDDO
c            ENDIF

         CASE(6) ! A ship in a quasi-steady case + act disk

         DO IGR = 1,NGRIFL
            IF(IGRSLAVE(IGR) == 0) THEN
C ... Sting forces and moments
               OSKU(IGR)%FXE = -OSKU(IGR)%CX-OSKU(IGR)%FXT-
     &              GX*OSKU(IGR)%RMASS
               OSKU(IGR)%FYE = -OSKU(IGR)%CY-OSKU(IGR)%FYT-
     &              GY*OSKU(IGR)%RMASS
               OSKU(IGR)%MXE = -OSKU(IGR)%CMX-OSKU(IGR)%MXT
               OSKU(IGR)%MZE = -OSKU(IGR)%CMZ-OSKU(IGR)%MZT
C ... Damp forces and moments           
               OSKU(IGR)%FZE = DAMPC1(IGR)*OSKU(IGR)%RMASS*
     &              SQRT(ABS(GZ/CHLREF))*OSKU(IGR)%VZ
               OSKU(IGR)%MYE = DAMPC2(IGR)*OSKU(IGR)%IZ*
     &              SQRT(ABS(GZ/CHLREF))*OSKU(IGR)%ROTY
C ... Writes data to the TRAJECTORY file
               CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),
     &              ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     &              THETAR(IGR),PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,
     &              TRMODE(IGR),ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
C ... Calls PARTICLE routine
               CALL PARTICLE(OSKU,XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     &              THETAR(IGR),PHIR(IGR),PSIM(IGR),IGR,
     &              0.0,DTB,TRMODE(IGR))
            ELSEIF(IGRSLAVE(IGR) /= 0) THEN
C ... Writes data to the TRAJECTORY file
               CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),
     &              ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),
     &              THETAR(IGR),PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,
     &              TRMODE(IGR),ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
C ... Calculate new location and angle of the actuator disk
               ZERO = 0.0
               ZCGD(IGR) = ZCGIS(IGR)+
     &              DRAUGHTI(IGRSLAVE(IGR))-DRAUGHT(IGRSLAVE(IGR)) ! add draught
               CALL EULFIN(PSIR(IGRSLAVE(IGR)),ZERO,ZERO, ! rotation
     &              -(XCGIS(IGR)-XCG(IGRSLAVE(IGR))),
     &              -(ZCGD(IGR)-ZCG(IGRSLAVE(IGR))),
     &              -(YCGIS(IGR)-YCG(IGRSLAVE(IGR))),XCGA(IGR),
     &              YCGA(IGR),ZCGA(IGR),1) 
               XCG(IGR)    = XCGA(IGR)+XCG(IGRSLAVE(IGR)) ! new locations
               YCG(IGR)    = YCGA(IGR)+YCG(IGRSLAVE(IGR))
               ZCG(IGR)    = ZCGA(IGR)+ZCG(IGRSLAVE(IGR))
               ROTB1A(IGR) = ROTB1S(IGR)-PSIR(IGRSLAVE(IGR))
               ROTB1(IGR)  = ROTB1A(IGR)
            ENDIF
         ENDDO

      END SELECT

      ENDIF  ! IPRO === 1

C***********************************************************************
C ... Next update grid and perform a grid deformation if required
C***********************************************************************

      IF (PARALLEL .AND. NGRIFL > 0) THEN ! Position update
         CALL MPI_BCAST(XCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR) 
         CALL MPI_BCAST(YCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCG,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETAR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIR,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIM,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTA1,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTB1,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CONEA,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%FXT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%FYT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%FZT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%MXT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%MYT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%MZT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FDSP,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RTMSP,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FXTSP,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETACOL,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THCYCLON,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THCYCLAT,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QTH,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QTH_A,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QBE_A,NGRIFL+10,MPI_REAL8,0,
     +        MPI_COMM_WORLD,IERR)
      ENDIF

      IS      = 1
      INDXFS  = 0

      DO 2000 N = 1,NBLOCK      ! BLOCK LOOP 1 BEGINS

         M   = 1
         IG1  = IG(M,N) 
         IC1  = IC(M,N) 
         IF1  = JF(1,N)
         NBG  = NPROCE(1+N,IPRO)
         NTOT1   = NTOT(1,N)
       
         NLVL   = 2**(LEVEL-1)
         istrid = imax(1,n) + 2*IN
         jstrid = jmax(1,n) + 2*jN
         ILE    = ISTRID*JSTRID

         KU1  = IT(M,N)
         KU2  = IL(M,N)
         KU3  = IK(M,N)
         KU4  = IMAX(M,N)
         KU5  = JMAX(M,N)
         KU6  = KMAX(M,N)
       
C***********************************************************************
C ... Update positions in quasi-steady cases
C***********************************************************************

      IGR    = IGRID(NBG,3)
 
      IF (IGRID(NBG,1) == 22) THEN ! Ship position is updated 

c      GOBACK = .FALSE.

C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0,0.0,0.0,PSIRI(IGR),THETARI(IGR),
     &        PHIRI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Set draught of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0,0.0,0.0,
     &        0.0,0.0,DRAUGHTI(IGR)-DRAUGHT(IGR),0.0,
     &        0.0,0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),XCGIS(IGR),YCGIS(IGR),ZCGIS(IGR),
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),0.0,
     &        0.0,ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) >= 37 .AND. IGRID(NBG,1) <= 39) THEN !Helicopter blade

c         GOBACK = .FALSE.
C ... Set location(s) and attitude(s) of the moving grid(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, PSIRI(IGR),PHIRI(IGR),
     &        THETARI(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
C ... Update location(s) and rotation angle(s) of the moving grid(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),PHIR(IGR),
     &        THETAR(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),2)
c      endif
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ELSEIF(IGRID(NBG,1) == 49 .AND. IGRID(NBG,4) /= 0 
     &        .AND. MVSHIP) THEN ! ACT disk in Ship case
C ... Set location(s) and attitude(s) of the actuator disk(s) to "zero position"
         CALL REVPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XORI(IG1),YORI(IG1),
     &        ZORI(IG1),NTOT(1,N),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &        0.0, 0.0, 0.0, 0.0,ROTA1I(IGR),
     &        ROTB1I(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
C ... Deformation of the the actuator disk(s)
c         IF(IGRID(NBG,2) == 15) THEN     t\E4m\E4 ei alkanut toimimaan en tied\E4 
c            CALL CONEDISK(XGRI(IG1),YGRI(IG1),ZGRI(IG1),              miksi?
c     &           XGRI(IG1),YGRI(IG1),ZGRI(IG1),NTOT(1,N),
c     &         0,0,0,0,1.0,0,0,0,0,0,1.0,0,CONEAI(IGR)-CONEA(IGR))
c         ENDIF
C ... Update location(s) and rotation angle(s) of the actuator disk(s) 
         CALL NEWPOS(XGRI(IG1),YGRI(IG1),ZGRI(IG1),XGRI(IG1),YGRI(IG1),
     &        ZGRI(IG1),NTOT(1,N),0.0, 0.0, 0.0,
     &        XCG(IGR),YCG(IGR),ZCG(IGR),0.0,ROTA1(IGR),
     &        ROTB1(IGR),ICYCLE,IPRINT,IREPEA(4),NREPEA(4),1)
         XCO(IG1:IG1+NTOT1-1) = XGRI(IG1:IG1+NTOT1-1)
         YCO(IG1:IG1+NTOT1-1) = YGRI(IG1:IG1+NTOT1-1)
         ZCO(IG1:IG1+NTOT1-1) = ZGRI(IG1:IG1+NTOT1-1)

      ENDIF

C***********************************************************************
C ... Deform grids in quasi-steady cases (Here helicopter blade)
C***********************************************************************

      IF(IGRID(NBG,2) >= 5 .AND. IGRID(NBG,2) < 9) THEN   ! WSHAPE free-surface
   
c      GOBACK = .FALSE.

      CALL WSHAPE(XGRI(IG1),YGRI(IG1),ZGRI(IG1),  
     &     XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     &     IHF(1,M),ICON(IC1),NPATCH(N),WAVES,ICYCLE,
     &     IPRINT,IREPEA,NREPEA,NBG)

      ELSEIF(IGRID(NBG,2) == 9) THEN ! GRID3 deformation model (single pre.)

c      GOBACK = .FALSE.

C *** FREE SURFACE, GRID3: Grid updating by Lehtimaki (object file libgrid.a)

       XSP(1:MAXFS) = XCO(1:MAXFS)
       YSP(1:MAXFS) = YCO(1:MAXFS)
       ZSP(1:MAXFS) = ZCO(1:MAXFS)
        
       osink  = sink
       otrim  = trimm
       INDXFS = 1
       JH1    = ICONH(N)
       JG1    = ICOGH(JH1)

       CALL GRID3(XSP(IG1),YSP(IG1),ZSP(IG1),
     & IMAX(1,N),JMAX(1,N),KMAX(1,N),XI(IG1),ETA(IG1),ZETA(IG1),
     & WH(IF1),XHULL(JG1),YHULL(JG1),ZHULL(JG1),
     & LHULL(JH1),MMHUL(JH1),MHULL(JH1),
     & IBOTGR(N),ITOPGR(N),JBOTGR(N),JTOPGR(N),KBOTGR(N),KTOPGR(N),
     & PDIFF(IG1),ZZTOP(IG1),IGLOBL,NGLOBL,CX,CY,CZ,CMY,XMOM,YMOM,ZMOM,
     & FRSVEL,GZ,GML,NLVL,XVL(IG1),YVL(IG1),ZVL(IG1),ICYCLE,ICMAX,
     $ AREF,CHLREF,REFPRE,INDXFS,N,IPRO)
       
       XCO(1:MAXFS) = XSP(1:MAXFS)
       YCO(1:MAXFS) = YSP(1:MAXFS)
       ZCO(1:MAXFS) = ZSP(1:MAXFS)

      ELSEIF(IGRID(NBG,2) > 20) THEN 
      GOBACK= .FALSE.
         WRITE(*,*) 'Nblock = ',NBG,' IGRID(NBG,2)=',IGRID(NBG,2)
         WRITE(*,*) 'Deformation is not yet possible'
         WRITE(*,*) 'Exiting ...'
         STOP

C     CALL CHAGRI 
C     ***********

C***********************************************************************
      ENDIF ! Wave (or blade) deformation
C***********************************************************************


C ... EXTRAPOLATE CORNER POINT COORDINATES AND FIND PATCH CORNERS

      CALL GRIDEX(XCO(IG1),YCO(IG1),ZCO(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),ICON(IC1),NPATCH(N),XCOL(1,1,IS),IN,JN,KN)
      IS  = IS + NPATCH(N)         

 2000 CONTINUE                  ! END OF BLOCK LOOP 2

      IF(INDXFS == 1) THEN ! Probably ineffective
       CALL FS_COVOL(ICON,XCO,YCO,ZCO,IMAX,JMAX,KMAX,M,NPATCH,IN,JN,KN,
     &    XXI,YXI,XETA,YETA,XC,YC,ZC,IG,IC,JF,MGM,NBLOCK)
      ENDIF

      IF(GOBACK) RETURN
**********************************************************************

C ... CONNECT THE CORNER POINT COORDINATES
      
      CALL MIRC2(XCO,YCO,ZCO,1,1,NBLOCK)

      CALL CONECX(XCO,1,1,NBLOCK,1)
      CALL CONECX(YCO,1,1,NBLOCK,1)
      CALL CONECX(ZCO,1,1,NBLOCK,1)
        
      DO 3000 N = 1,NBLOCK      ! BLOCK LOOP 1 BEGINS

      IG1  = IG(1,N)
      M    = 1
      NBG  = NPROCE(1+N,IPRO)

      KU1  = IT(M,N)
      KU2  = IL(M,N)
      KU3  = IK(M,N)
      KU4  = IMAX(M,N)
      KU5  = JMAX(M,N)
      KU6  = KMAX(M,N)

C ... Update the cell center point coordinates (ESa 4.2.1997)

c      CALL CELLCP(XCO(IG1),YCO(IG1),ZCO(IG1),XC(IG1),YC(IG1),ZC(IG1),
c     +            IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,PDIFF(IG1),
c     +            FRESUL,FRSDEN)

      CALL CELLCP(XCO(IG1),YCO(IG1),ZCO(IG1),XC(IG1),YC(IG1),ZC(IG1),
     +            OMEGA(NBG),OMEGAX(NBG),OMEGAY(NBG),OMEGAZ(NBG),
     +            CENAX(NBG),CENAY(NBG),CENAZ(NBG),UROTCP(IG1),
     +            VROTCP(IG1),WROTCP(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +            IN,JN,KN,PDIFF(IG1),FRESUL,FRSDEN)
      
      IF(ICYCLE == IPRINT) THEN
      WRITE(3,*)
      WRITE(3,*) 'IN NEWGRI:'
      CALL PRINYD(3,XCOC,XCO(IG1),KU1,KU2,0,KU3,KU4,KU5,KU6,NBG,M)
      CALL PRINYD(3,YCOC,YCO(IG1),KU1,KU2,0,KU3,KU4,KU5,KU6,NBG,M)
      CALL PRINYD(3,ZCOC,ZCO(IG1),KU1,KU2,0,KU3,KU4,KU5,KU6,NBG,M)
      ENDIF
3000  CONTINUE
       
C ... EXTRAPOLATION AND CONNECTION IS NOT NEEDED

C***********************************************************************
c      IF(IPRINT-1 == ICYCLE) WRITE(46,*)
c      IF(IPRINT-1 == ICYCLE) WRITE(46,*) ' CALL FROM NEWGRI'

      IF(IREPEA(3) < NREPEA(3)) THEN ! IREPEA(3) for NEWGRI
        WRITE(45,*)
        WRITE(45,*) 'GRIVAL CALL FROM NEWGRI:'
        WRITE(45,'(1X,24("-"))')
      ELSE IF(IREPEA(3) == NREPEA(3)) THEN
        WRITE(45,*)
        WRITE(45,*) 'GRIVAL CALL FROM NEWGRI:'
        RIVI ="(' THIS HAS BEEN REPEATED',I3,' TIMES. SILENCE')"
        WRITE(45,'(1X,24("="))')
        WRITE(45,RIVI) IREPEA(3)
        WRITE(45,'(1X,40("-"))')
      ENDIF
      
      CALL GRIVAL(1,LDISTAN,'NEWGRI')
       
      IREPEA(3) = IREPEA(3) + 1

      RETURN
      END SUBROUTINE NEWGRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SURGRI(M)

      USE MAIN_ARRAYS, ONLY : TEMP,U,V,W,E,EPS2,VIST,P,PDIFF,RK,
     +     REPS,FI,BIJ,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,
     +     A3ZA,XC,YC,ZC,XCO,YCO,ZCO,IT,IL,IK,IMAX,JMAX,KMAX,A1,
     +     A2,A3,WH,IBOTGR,ITOPGR,KTOPGR,LHULL,MMHUL,MHULL,VOL,
     +     DTL,UROT,VROT,WROT,ICONH,ICOGH,IHULL,WAVEH,RO,RM,RN,RW,
     +     IG,JF,IC,IHF,NPATCH,ICON,IGRID,NPROCE

      USE NS3CO,       ONLY : GX,GY,GZ,GROUND,FRSDEN,IFSBC,INWH,
     +     NBLOCK,JFIRST,ICFST,ITURB,IOLD,ICMAX,ICYCLE,DTWMAX,
     +     WHMAX,WHMIN,FLODWH,NFSD,DWMAX,DWMV,PRN,IPRINT,SUMDWH,
     +     IN, JN, KN

      USE INTEGERS,    ONLY : MAXEB,MAXSB,IPRO

      USE CHARACTERS
 
      IMPLICIT NONE

      INTEGER :: M, N, IG1, IF1, IC1, NGL, III


      IF(IFSBC == 1) THEN ! Original mod 

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*) 'INWH=',inwh
       DO N = 1,NBLOCK
       IG1 = IG(1,N)
       CALL PRINYS(3,pdifC,pdiff(IG1),IT(1,N),
     + IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),N,1)
       ENDDO
      ENDIF
         
      CALL FS_FRE(RO,RM,RN,RW,E,EPS2,VIST,RK,REPS,FI,MAXSB,
     +     A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +     1,NBLOCK,1.,VOL,XC,YC,ZC,A3,DTL,1,
     +     LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR,INWH,
     +     UROT,VROT,WROT,XCO,YCO,ZCO,ICONH,ICOGH,IHULL,PRN,
     +     NFSD,ICYCLE,ICFST,DTWMAX,DWMV,ITURB,JFIRST,IOLD,ICMAX,
     +     DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH)
      IF(ICYCLE == IPRINT) WRITE(3,*) 'fs_fre'
        
c      CALL MIR(RO,RM,RN,RW,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,MAXSB,
c     +     A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
c     +     1,NBLOCK,1.,PRC)
      CALL FOU(RO,RM,RN,RW,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,MAXSB,
     +     A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     +     1,NBLOCK,1.,VOL,XC,YC,ZC,WH)

      IF(ICYCLE == IPRINT) THEN
       DO N = 1,NBLOCK
       IG1 = IG(1,N)
       CALL PRINYS(3,pdifC,pdiff(IG1),IT(1,N),
     + IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),N,1)
       ENDDO
      ENDIF

      ELSE IF(IFSBC == 2) THEN
         STOP 'SURGRI: Level-set does not exist man! Change IFSBC'
      ELSE IF(IFSBC == 3) THEN 
         CALL FRE(M,1,NBLOCK,1.,1) ! Only here surface calculation ISUR = 1
      ELSE IF(IFSBC == 44) THEN 
         CALL FRE(M,1,NBLOCK,1.,2) ! Only here surface calculation ISUR = 2
      ENDIF

C ... Connect the wave heights (MIR, SOL and CON)
       
      CALL CONEPA(WAVEH,'FRE')

C ... Use the WH array for IFSBC = 1 or for Hatchback

      IF(IFSBC == 1) THEN 
        DO N = 1,NBLOCK
        IF1  = JF(1,N)
        IC1  = IC(1,N)
        CALL FS_SETWAVEHBACK(WAVEH,WH(IF1),IHF(1,1),ICON(IC1),
     +   NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ICYCLE)
        ENDDO
      ELSE IF(IFSBC == 4) THEN
        DO N = 1,NBLOCK
        NGL  = NPROCE(1+N,IPRO)
        IF(IGRID(NGL,2) == 9) THEN
        IF1  = JF(1,N)
        IC1  = IC(1,N)
        CALL FS_SETWAVEHBACK(WAVEH,WH(IF1),IHF(1,1),ICON(IC1),
     +   NPATCH(N),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ICYCLE)
        ENDIF ! IGRID(NGL,2) == 9
        ENDDO
      ENDIF

C ... Update the nodal height values (array WAVES)

      CALL FS_VERPAT(1,NBLOCK,M)

      RETURN
      END SUBROUTINE SURGRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUSS

      USE INTEGERS,    ONLY : MAXB,MAXW,MAXSB,IB,IPRO
      USE MAIN_ARRAYS, ONLY : RO,RM,RN,RW,E,U,V,W,IMAX,JMAX,KMAX,
     +    ZZZ,F1R,F1RM,JLOC,JTRA,APP,RLOLIM,UPPLIM,IG,IR,JF,IH,IC,
     +    XC,YC,ZC,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,IDI1,IDI2,
     +    IDI3,ICON,NPATCH,P,PDIFF,VIS,VIST,CH,EPS2,TEMP,C,DRDP,
     +    DRDH,RK,REPS,DDEPS,RCON,UWALL,VWALL,WWALL,FI,HAT1,HAT2,
     +    A1,A2,A3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +    D1,D2,D3,IDER,IHF,JSTATE,NPROCE,TWALL,BLKS,PRO,VAR,VTRAN,TRM,
     +    PRC
      USE NS3CO,       ONLY : IN,JN,KN,NBLOCK,ICYCLE,GAMMA,E0REF,
     +    T0REF,FRSDEN,FRSPRE,RGAS,T0,ITURB,NSCAL,TLOLIM,TUPLIM,
     +    ISTATE,MULPHL,REFLECL,TRANSL

      IMPLICIT NONE 

      INTEGER :: M, IPR, N, NS, IPRI, IG1, IG2, IR2, IF2, IH2, IC2,
     &           NGL, KSTATE, IM2, IT2, IPC2, IPRESC       

C ... Makes an approximation for the Reynolds stresses stored as scalars.
                        
      M   = 1
      IPR = 1
                        
      DO 500 N = 1,NBLOCK
      NGL     = NPROCE(1+N,IPRO) !Global block number
      IPRESC  = BLKS(NGL)%IPRESC
      KSTATE = JSTATE(NGL,1) 
      IPRI    = ICYCLE
      IG2     = IG(M,N)
      IR2     = IR(M,N)
      IF2     = JF(M,N)
      IH2     = IH(M,N)
      IC2     = IC(M,N)
      IM2     = 1
      IT2     = 1
      IPC2    = 1
      IF(MULPHL) IM2 = IG2
      IF(TRANSL) IT2 = IG2
      IF(IPRESC == 1) IPC2 = IG2
                        
      CALL NEWVEL(RO(IG2),RM(IG2),RN(IG2),RW(IG2),U(IG2),V(IG2),
     + W(IG2),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN)

      CALL REFLEC(U(IG2),V(IG2),W(IG2),XC(IG2),YC(IG2),ZC(IG2),
     1 IMAX(1,N),JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),
     2 JTOP(1,N),KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     3 GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,NGL,NPATCH(N),ICON(IC2),
     4 RO(IG2),P(IG2),PDIFF(IG2),RM(IG2),RN(IG2),RW(IG2),E(IG2),
     5 VIS(IG2),VIST(IG2),CH(IG2),EPS2(IG2),TEMP(IG2),C(IG2),DRDP(IG2),
     6 DRDH(IG2),RK(IR2),REPS(IR2),DDEPS(IR2),PRO(IM2),VAR(IM2),
     7 RGAS,T0,ITURB,KSTATE,IB,MAXB,RCON(IC2),UWALL(1),VWALL(1),
     8 WWALL(1),TWALL(1),IHF(1,M),TLOLIM,TUPLIM,BLKS(NGL)%SOLUTION_TYPE,
     9 MULPHL,REFLECL,TRANSL,TRM(IT2),
     1 A1XA(IG2),A1YA(IG2),A1ZA(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),
     2 A2XA(IG2),A2YA(IG2),A2ZA(IG2),PRC(IPC2),IPRESC)

      CALL TURFUN(FI(IG2,1),FI(IG2,2),FI(IG2,3),FI(IG2,4),FI(IG2,5),
     1 FI(IG2,6),U(IG2),V(IG2),W(IG2),RK(IR2),REPS(IR2),DDEPS(IR2),
*     1 EPS2(IG2),VIST(IG2),VIS(IG2),HAT1(IG1),HAT2(IG1),
     1 EPS2(IG2),VIST(IG2),VIS(IG2),HAT1(IG2),HAT2(IG2),
     2 A1(IG2),A2(IG2),A3(IG2),VOL(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     3 A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     4 D1(IG2),D2(IG2),D3(IG2),XC(IG2),YC(IG2),ZC(IG2),
     5 IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ITURB,1,ZZZ,VTRAN(IR2),
     6 MAXW,IDER(N),NGL)

      CALL LOLISC(FI(IG2,1),IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,
     +   RLOLIM,UPPLIM,MAXSB,NSCAL)
 500  CONTINUE

      DO 900 NS = 1,NSCAL
         CALL CONEC(FI(1,NS),F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
 900  CONTINUE

C ... END BOUSSINESQ-APPROXIMATION

      RETURN
      END SUBROUTINE BOUSS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONEPA(WAVEH,SURCHA)

      USE NS3CO,       ONLY : NBLOCK, IN, JN, KN
      USE INTEGERS,    ONLY : MAXB,IREPEA,NREPEA,IPRO
      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,NPATCH,ICON,IHF,IG,JF,IC,
     +   NPROCE,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,IT,IK,IL

      IMPLICIT NONE

      INTEGER :: M,N,IG1,IC1,IF1,NGL,ku1,ku2,ku3,ku4,ku5,ku6,iii
      REAL, ALLOCATABLE  :: DPU1(:)
      REAL, DIMENSION(*) :: WAVEH
      CHARACTER(LEN=3) :: SURCHA

      ALLOCATE(DPU1(MAXB))

      DPU1(1:MAXB) = 0.

      M = 1  ! Only densiest grid level

C ... Put everything into a block array
      
      DO N = 1,NBLOCK

         NGL = NPROCE(1+N,IPRO) ! Global block number
         IG1 = IG(M,N)
         IF1 = JF(M,N)
         IC1 = IC(M,N)

         CALL PATVOL(WAVEH,DPU1(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN,NGL,M,NPATCH(N),ICON(IC1),IHF,MAXB,
     +        IREPEA,NREPEA,SURCHA,.FALSE.)

         KU1 = IT(M,N)
         KU2 = IL(M,N)
         KU3 = IK(M,N)
         KU4 = IMAX(M,N)
         KU5 = JMAX(M,N)
         KU6 = KMAX(M,N)

c      CALL IPRINT0(DPU1,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,466)   
c      CALL PRINYD(466,'dpu',DPU1(IG1),KU1,KU2,0,KU3,KU4,KU5,KU6,N,M)

      ENDDO

C ... Mirror (IBC = 4)

      DO N = 1,NBLOCK

         IG1 = IG(1,N)
         IC1 = IC(1,N)

      CALL MIRDTE(DPU1(IG1),4,IHF(1,1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),KTOP(1,N),M,
     + IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,ICON(IC1),NPATCH(N),10,MAXB,1,2) 

      ENDDO

C ... Connection twice in order to fill in the empty corners

      CALL CONECX(DPU1,M,1,NBLOCK,0)
      CALL CONECX(DPU1,M,1,NBLOCK,0)

C ... Replacement back to patch arrays

      DO N = 1,NBLOCK

         NGL = NPROCE(1+N,IPRO) ! Global block number
         IG1 = IG(M,N)
         IF1 = JF(M,N)
         IC1 = IC(M,N)

         CALL PATVOL(WAVEH,DPU1(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN,NGL,M,NPATCH(N),ICON(IC1),IHF,MAXB,
     +        IREPEA,NREPEA,SURCHA,.TRUE.)
         
      ENDDO

      RETURN
      END SUBROUTINE CONEPA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONEPA_SP(WAVEH,SURCHA)

      USE NS3CO,       ONLY : NBLOCK, IN, JN, KN
      USE INTEGERS,    ONLY : MAXB,IREPEA,NREPEA,IPRO
      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,NPATCH,ICON,IHF,IG,JF,IC,
     +   NPROCE,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,IT,IK,IL,
     +   F1R,F1RM,JTRA,JLOC,APP

      IMPLICIT NONE

      INTEGER :: M,N,IG1,IC1,IF1,NGL,ku1,ku2,ku3,ku4,ku5,ku6,iii
      REAL, ALLOCATABLE  :: DPU1(:)
      REAL, DIMENSION(*) :: WAVEH
      CHARACTER(LEN=3)   :: SURCHA

      ALLOCATE(DPU1(MAXB))

      DPU1(1:MAXB) = 0.

      M = 1  ! Only densiest grid level

C ... Put everything into a block array

      DO N = 1,NBLOCK

         NGL = NPROCE(1+N,IPRO) !Global block number
         IG1 = IG(M,N)
         IF1 = JF(M,N)
         IC1 = IC(M,N)

         CALL PATVOL_SP(WAVEH,DPU1(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN,NGL,M,NPATCH(N),ICON(IC1),IHF,MAXB,
     +        IREPEA,NREPEA,SURCHA,.FALSE.)
         KU1 = IT(M,N)
         KU2 = IL(M,N)
         KU3 = IK(M,N)
         KU4 = IMAX(M,N)
         KU5 = JMAX(M,N)
         KU6 = KMAX(M,N)

c      CALL IPRINT0(DPU1,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,466)   
c      CALL PRINYD(466,'dpu',DPU1(IG1),KU1,KU2,0,KU3,KU4,KU5,KU6,N,M)

      ENDDO

C ... Mirror (IBC = 4)

      DO N = 1,NBLOCK

         IG1 = IG(1,N)
         IC1 = IC(1,N)

      CALL MIRDTE_SP(DPU1(IG1),4,IHF(1,1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),KTOP(1,N),M,
     + IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,ICON(IC1),NPATCH(N),10,MAXB,1,2) 

      ENDDO

C ... Connection twice in order to fill in the empty corners

      CALL CONEC(DPU1,F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)
      CALL CONEC(DPU1,F1R,F1RM,JLOC,JTRA,APP,M,1,NBLOCK,0,0)

C ... Replacement back to patch arrays

      DO N = 1,NBLOCK

         NGL = NPROCE(1+N,IPRO) !Global block number
         IG1 = IG(M,N)
         IF1 = JF(M,N)
         IC1 = IC(M,N)

         CALL PATVOL_SP(WAVEH,DPU1(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +        IN,JN,KN,NGL,M,NPATCH(N),ICON(IC1),IHF,MAXB,
     +        IREPEA,NREPEA,SURCHA,.TRUE.)

      ENDDO

      RETURN
      END SUBROUTINE CONEPA_SP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONEC(VOL,APU,APU2,JLOC,JTRA,APP,M,N1,N2,IEXT,MPE)

      USE INTEGERS,    ONLY : NPPV,NCON,IPRO
      USE MAIN_ARRAYS, ONLY : ICON,MGRID,ITAG,NCPAT,ZZZ,NPATCH,IW,
     & IG,IC,IWAPP,JAPP,IMAX,JMAX,KMAX
      USE NS3CO           

      IMPLICIT NONE

      REAL    :: VOL(*), APU(*), APU2(*), APP(*)
      INTEGER :: JLOC(*), JTRA(*)
      INTEGER :: M, N1, N2, IEXT, MPE, LC, LPO, L, N, IG1, IC1, NN

C
C ... THIS SUBROUTINE FIRSTLY EXPORTS ALL THE PATCHES
C
C ... IEXT = 0 IF CELL CENTERED DATA
C ... IEXT = 1 IF CORNER DATA
C ... MPE  = 0 IF NORMAL CONNECTIVITY AND CYCLIC AND MIRROR AND 
C ...             FREE SURFACE DATA
C ... MPE  = 1 IF ONLY NORMAL CONNECTIVITY (no cyclic and periodic)
C ... MPE  = 2 IF ONLY CYCLIC              (no normal)

C ... COUNTER FOR THE CORRESPONDANCE     ! BLOCK LOOP 1

      DO 7800 N = 1,NBLOCK
      DO 7800 L = 1,NPPV
7800  ITAG(L,N)   = 0

      LC = NCPAT(1)
      DO 1000 LPO = 1,NCON
         L      = NCPAT(LPO)        ! PATCH NUMBER
         N      = ICON((L-1)*IC9+23) ! BLOCK NUMBER
         IF(M <= MGRID(N) .AND. N >= N1 .AND. N <= N2) THEN
         LC = L
         IF(N > 1) THEN
            DO NN = 1,N-1
               LC = LC - NPATCH(NN) ! BLOCK LOCAL PATCH NUMBER
            ENDDO
         ENDIF

         IG1    = IG(M,N)
         IC1    = IC(M,1)
         CALL EXPVOL(VOL(IG1),ZZZ,APU,ICON(IC1),NPATCH(N),IW(LC,N,M),
     +   IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG,IN,JN,KN,NPPV,IEXT,M,MPE,
     +   IPRO,N,L)
      ENDIF
1000  CONTINUE

      DO 7950 N = N1,N2
         IC1     = IC(M,N)
         CALL ITAGMO(ITAG,ICON(IC1),NPATCH(N),NPPV,N)
 7950 CONTINUE

C ... AND THEN IMPORTS ALL THE PATCHES.
C ... PROCESSES SHOULD BE SYNCHRONIZED BEFORE DATA IMPORT.(= EOBL)
C ... IMPORT COULD BE PERFORMED SOMEWHERE ELSE BEFORE FLUX EVALUATION
C ... ITAG ACTIVATES THE BLOCK, WHERE RECEIVE IS NEEDED

      DO 8000 N = 1,NBLOCK               ! BLOCK LOOP 3
      IF(M <=  MGRID(N)) THEN
      IG1       = IG(M,N)
      IC1       = IC(M,N)
        CALL IMPVOL(VOL(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +  IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,IEXT,IPRO,
     +  MPE,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N))
      ENDIF
8000  CONTINUE

      RETURN
      END SUBROUTINE CONEC
C
C ----------------------------------------------------------------------
C --- Subroutine for handling connective boundary conditions -----------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONECX(VOL,M,N1,N2,IEXT)

      USE INTEGERS,    ONLY : NPPV,NCON,IPRO,MAXB,MAXW
      USE MAIN_ARRAYS, ONLY : ICON,MGRID,ITAG,NCPAT,NPATCH,IW,
     & IG,IC,IWAPP,IMAX,JMAX,KMAX
      USE NS3CO
         
      IMPLICIT NONE

      INTEGER :: M,N1,N2,IERR,KOKO1,IEXT,MPE,N,L,LC,LPO,NN,IG1,IC1

      REAL    :: BKOKO1

      REAL :: VOL(*)
      REAL, ALLOCATABLE :: DPU1(:), DPU2(:), DPU3(:)

C ... THIS SUBROUTINE FIRSTLY EXPORTS ALL THE PATCHES FOR CORNER POINTS
C
C ... IEXT = 0 IF CELL CENTERED DATA
C ... IEXT = 1 IF CORNER DATA

      ALLOCATE(DPU1(MAXB), DPU2(MAXB), DPU3(MAXW),STAT=IERR)
      
cc      MREQ1  = IWAPP(1,NBLOCK+1)

      IF(ICYCLE <= 1) THEN ! Print to the MEMORY file
      WRITE(45,*)
      IF(IERR == 0) THEN
        KOKO1  = 2*MAXB + MAXW
        BKOKO1 = (REAL(KOKO1,4)*8.)/1048576.
        WRITE(45,9050) KOKO1,BKOKO1
9050    FORMAT('In CONECX: DPU1, DPU2 and DPU3 succesfully allo',
     +  'cated.'/'Size (DP) = ',I8,' = ',F6.2,' Mbytes')
        ELSE
        WRITE(45,*) 'Not enough memory in CONECX. Exiting...'
        STOP
      ENDIF
      ENDIF !(ICYCLE <= 1)

c      IEXT = 1 ! From parameter list
c      MPE  = 0 ! Unused

C ... COUNTER FOR THE CORRESPONDANCE     ! BLOCK LOOP 1

      DO 7800 N = 1,NBLOCK
      DO 7800 L = 1,NPPV
7800  ITAG(L,N)   = 0
       
      LC = NCPAT(1)

      DO 1000 LPO = 1,NCON
         L      = NCPAT(LPO)        ! PATCH NUMBER
         N      = ICON((L-1)*IC9+23) ! BLOCK NUMBER
         IF(M <= MGRID(N) .AND. N >= N1 .AND. N <= N2) THEN
         LC = L
         IF(N > 1) THEN
            DO NN = 1,N-1
               LC = LC - NPATCH(NN) ! BLOCK LOCAL PATCH NUMBER
            ENDDO
         ENDIF

         IG1    = IG(M,N)
         IC1    = IC(M,1)
         CALL EXPVOD(VOL(IG1),DPU3,DPU1,ICON(IC1),NPATCH(N),IW(LC,N,M),
     +   IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG,IN,JN,KN,NPPV,M,IPRO,N,L,
     +   IEXT)
         ENDIF ! M <= MGRID(N) .AND. N >= N1 .AND. N <= N2 
1000  CONTINUE
      
      DO 7950 N = N1,N2
         IC1     = IC(M,N)
         CALL ITAGMO(ITAG,ICON(IC1),NPATCH(N),NPPV,N)
 7950 CONTINUE

C ... AND THEN IMPORTS ALL THE PATCHES.
C ... PROCESSES SHOULD BE SYNCHRONIZED BEFORE DATA IMPORT.(= EOBL)
C ... IMPORT COULD BE PERFORMED SOMEWHERE ELSE BEFORE FLUX EVALUATION
C ... ITAG ACTIVATES THE BLOCK, WHERE RECEIVE IS NEEDED

      DO 8000 N = 1,NBLOCK               ! BLOCK LOOP 3
      IF(M <=  MGRID(N)) THEN
        IG1     = IG(M,N)
        IC1     = IC(M,N)
        CALL IMPVOD(VOL(IG1),DPU3,DPU1,DPU2,ICON(IC1),NPATCH(N),
     +  IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,
     +  KN,NPPV,IPRO,N,M,IEXT)
      END IF
8000  CONTINUE
      DEALLOCATE (DPU1, DPU2, DPU3,STAT=IERR)

      IF(ICYCLE <= 1) THEN ! Print to the MEMORY file
      WRITE(45,9100) IERR
9100  FORMAT('In CONECX: DPU1... deallocated. IERR = ',I6)
      ENDIF
      IF(IERR /= 0) STOP

      RETURN
      END SUBROUTINE CONECX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONEX(XCO,YCO,ZCO,APU,APU2,JLOC,JTRA,APP,M,N1,N2,MAXB)

      USE INTEGERS,    ONLY : NPPV,NCON,IPRO,MAXW,NB
      USE MAIN_ARRAYS, ONLY : ICON,ITAG,IWAPP,JAPP,NCPAT,NPATCH,IW,
     & MGRID,IG,IC,IMAX,JMAX,KMAX
      USE NS3CO          

      IMPLICIT NONE

      INTEGER :: M, N1, N2, MAXB, IERR1, KOKO1, KOKO2,
     &           N, L, LC, LPO, NN, IG1, IC1, I, LL
      REAL :: BKOKO1,BKOKO2
      REAL :: APU(MAXB/2), APU2(MAXB/2)
      REAL :: APP(*)
      REAL :: XCO(MAXB), YCO(MAXB), ZCO(MAXB)
      REAL, ALLOCATABLE :: DPU1(:),DPU2(:),DPU3(:)
      INTEGER :: JLOC(*), JTRA(*)
      LOGICAL :: NOMATC

C ... THIS SUBROUTINE FIRSTLY EXPORTS ALL THE PATCHES

C ... JAPP is number of connective relation in patch (magnitude)
C ... IWAPP is the starting number for patch connection relations (pointter)

      NOMATC = .FALSE.
      ALLOCATE (DPU1(MAXB),DPU2(MAXB),DPU3(MAXW),STAT=IERR1)

      IF(ICYCLE <= 1) THEN ! Print to the MEMORY file
      WRITE(45,*)
      IF(IERR1 == 0) THEN
        KOKO1  = 2*MAXB + MAXW
        BKOKO1 = REAL(KOKO1,4)*8./1048576.
        KOKO2  = 2*MAXB
        BKOKO2 = REAL(KOKO2,4)*2./1048576.
        WRITE(45,9050) KOKO1,BKOKO1
9050    FORMAT('In CONEX: DPU1, DPU2 and DPU3 succesfully allo',
     +  'cated.'/'Size (DP) = ',I8,' = ',F6.2,' Mbytes')
        WRITE(45,9060) KOKO2,BKOKO2
9060    FORMAT('In CONEX: JLOC and JTRA succesfully allocated.'/
     + 'Size (INTEGER) = ',I8' = ',F6.2,' Mbytes')
        ELSE
        WRITE(45,*) 'Not enough memory in CONEX. Exiting...'
        STOP
      ENDIF
      ENDIF ! ICYCLE <= 1

C ... COUNTER FOR THE CORRESPONDANCE     ! BLOCK LOOP 1

      DO 7800 N = 1,NBLOCK
      DO 7800 L = 1,NPPV
7800  ITAG(L,N)   = 0

C ... do connection only for the multigrid level 1. All patc data are 
C ... transmited, but only desired are used.

      IF(M /= 1) WRITE(*,*) 'You are in trouble'
      LC = NCPAT(1)
      DO 1000 LPO = 1,NCON
         L = NCPAT(LPO)          ! PATCH NUMBER
         N = ICON((L-1)*IC9+23)  ! BLOCK NUMBER
         LC = L
         IF(N > 1) THEN
            DO NN = 1,N-1
               LC = LC - NPATCH(NN)  ! BLOCK LOCAL PATCH NUMBER
            ENDDO
      ENDIF
      
         IG1 = IG(M,N)
         IC1 = IC(M,1)
         CALL EXPVXX(XCO(IG1),YCO(IG1),ZCO(IG1),DPU3,DPU1,ICON(IC1),
     +        NPATCH(N),IW(LC,N,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG,
     +        IN,JN,KN,NPPV,M,IPRO,N,L)
1000  CONTINUE

      DO 7950 N = N1,N2
         IC1 = IC(M,N)
         CALL ITAGMO(ITAG,ICON(IC1),NPATCH(N),NPPV,N)
 7950 CONTINUE

C ... AND THEN Check if  ALL THE PATCHES fit to each other.
C ... ITAG ACTIVATES THE BLOCK, WHERE RECEIVE IS NEEDED

      DO 8000 N = 1,NBLOCK               ! BLOCK LOOP 3

      IF(M <= MGRID(N)) THEN

        IG1 = IG(M,N)
        IC1 = IC(M,N)

        CALL IMPVXX(XCO(IG1),YCO(IG1),ZCO(IG1),DPU3,DPU1,DPU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,M,N,NOMATC)

      ENDIF

8000  CONTINUE

      IF(NOMATC) WRITE(*,*) 
     + '  NBG PATCH  # of cell #con        face a    face b  DIR'

C ... Weight the fluxes over block connection in case of non-patched
C ... connection.  

C ... initialize the pointer array and the magnitude array
      
      DO N = 1,NB
      DO I = 1,NPPV
         JAPP(I,N)  = 1
         IWAPP(I,N) = 1
      ENDDO
      ENDDO

      DO 8100 N = 1,NBLOCK               ! BLOCK LOOP 4
      IF(M <=  MGRID(N)) THEN
        IG1     = IG(M,N)
        IC1     = IC(M,N)
        CALL APPRAI(XCO(IG1),YCO(IG1),ZCO(IG1),DPU3,APU2,DPU2,JLOC,JTRA,
     +       APP,ICON(IC1),JAPP(1,N),IWAPP,
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,MAXB)
        DO NN = N+1,NBLOCK+1
           DO LL = 1,NPPV
              IWAPP(LL,NN) = IWAPP(1,N+1)
           ENDDO
        ENDDO
      ENDIF
8100  CONTINUE

      DEALLOCATE (DPU1,DPU2,DPU3,STAT=IERR1)

      IF(ICYCLE <= 1) THEN ! Print to the MEMORY file
      WRITE(45,*)
      WRITE(45,9070) IERR1
9070  FORMAT('In CONEX: DPU1, DPU2 ... deallocated. IERR1 = ',I6)
      WRITE(45,*)
      ENDIF ! ICYCLE <= 1

      RETURN
      END SUBROUTINE CONEX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONNEC(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,FI,
     +S11,BIJ,PRO,VAR,TRM,BLKS,MAXEB,MAXSB,MAXSS,APU,APU2,JLOC,JTRA,APP,
     +M,N1,N2,MULPHC)

      USE TYPE_ARRAYS
      USE INTEGERS,    ONLY : NPPV,NCON,IPRO,MAXW,NB
      USE MAIN_ARRAYS, ONLY : ICON,ITAG,IWAPP,JAPP,NCPAT,NPATCH,IW,
     + MGRID,NTOT,IG,IC,IR,IQ,IMAX,JMAX,KMAX,NPROCE,ZZZ,F1RN,F1RW,F1E
      USE NS3CO            

      IMPLICIT NONE

      INTEGER M,N1,N2,MAXEB,MAXSB,MAXSS,N,L,LC,LPO,NN,IG1,IS1,IR1,IQ1,
     +     IC1,IP1,IGN,NEX,NS,NTOT1,NGL,IPHASE,I,IT1

      REAL :: RO(*), RM(*), RN(*), RW(*), E(*), EPS2(*), VIST(*),
     +        P(*), PDIFF(*), RK(*), REPS(*), FI(MAXSB,*), S11(MAXSS,*),
     +        APU(*), APU2(*), APP(*), BIJ(MAXEB,*), CH(*)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(BLOCKS)           :: BLKS(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      INTEGER :: JLOC(*), JTRA(*)

      CHARACTER(*) :: MULPHC

C ... This subroutine firstly export the main variables for the patches

C ... Counter for the correspondance     ! Block loop 1
   
      DO 7800 N = 1,NBLOCK
      DO 7800 L = 1,NPPV
7800  ITAG(L,N)   = 0

C ... do connection only for the multigrid level 1. All patc data are 
C ... transmited, but only desired are used.

      IF(M /= 1) WRITE(*,*) 'Something wrong in connections'
      LC = NCPAT(1)

      DO 1000 LPO = 1,NCON          ! Total number of connectins
         L      = NCPAT(LPO)        ! PATCH NUMBER
         N      = ICON((L-1)*IC9+23) ! BLOCK NUMBER
         LC = L
         IF(N > 1) THEN
            DO NN = 1,N-1
               LC = LC - NPATCH(NN) ! BLOCK LOCAL PATCH NUMBER
            ENDDO
      ENDIF

         IG1    = IG(M,N)
         IS1    = IG(M,N)
         IF(NSCAL == 0) IS1 = 1
         IR1    = 1
         IF(ITURB >= 10.AND.ITURB <= 19 .OR. ISTRES > 0) IR1 = IR(M,N)
         IQ1    = IQ(M,N)
         IC1    = IC(M,1)
         IP1    = 1
         IT1    = 1
         IF(MULPHL) IP1 = IG1
         IF(TRANSL) IT1 = IG1

         CALL EXPTOT(RO(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),EPS2(IG1),
     +      VIST(IG1),CH(IG1),P(IG1),PDIFF(IG1),RK(IG1),REPS(IG1),ITURB,
     +      ISTRES,FI(IS1,1),S11(IR1,1),BIJ(IR1,1),PRO(IP1),VAR(IP1),
     +     TRM(IT1),MAXEB,MAXSB,MAXSS,NSCAL,ZZZ,APU,ICON(IC1),NPATCH(N),
     +      IW(LC,N,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG,IN,JN,KN,
     +     NPPV,M,IPRO,N,L,MULPHL,TRANSL,F1RN,MULPHC,TWO_FLUIDL)
 1000 CONTINUE

      DO 7950 N = N1,N2
         IC1     = IC(M,N)
         CALL ITAGMO(ITAG,ICON(IC1),NPATCH(N),NPPV,N)
 7950 CONTINUE

C ... AND THEN IMPORTS THE MAIN ARRAYS FOR ALL THE PATCHES.
C ... PROCESSES SHOULD BE SYNCHRONIZED BEFORE DATA IMPORT.(= EOBL)
C ... IMPORT COULD BE PERFORMED SOMEWHERE ELSE BEFORE FLUX EVALUATION
C ... ITAG ACTIVATES THE BLOCK, WHERE RECEIVE IS NEEDED

      DO 8000 N = 1,NBLOCK               ! BLOCK LOOP 3
      IF(M <=  MGRID(N)) THEN
        IG1     = IG(M,N)
        IGN     = IG1 + NTOT(M,N) - 1
        IR1     = IR(M,N)
        IQ1     = IQ(M,N)
        IC1     = IC(M,N)
        NTOT1   = NTOT(M,N)
        NGL     = NPROCE(1+N,IPRO) !Global block number
      CALL IMPVO2(RO(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),0)
        CALL IMPVO2(RM(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),1)
        CALL IMPVO2(RN(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),2)
        CALL IMPVO2(RW(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),3)
        CALL IMPVO2(E(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),4)
        CALL IMPVO2(EPS2(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),5)
        CALL IMPVO2(VIST(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),6)
        CALL IMPVO2(CH(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),7)
        CALL IMPVO2(P(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),IW(1,1,M),
     +       IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,NPPV,
     +       IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),8)
        CALL IMPVO2(PDIFF(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,
     +       JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),9)
        NEX = 10
C ... K-E MODEL
         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
        CALL IMPVO2(RK(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),10)
        CALL IMPVO2(REPS(IG1),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),11)
        NEX = 12
         ENDIF
C ... SCALARS
      IF(NSCAL >= 1) THEN
      DO 720 NS = 1,NSCAL
        CALL IMPVO2(FI(IG1,NS),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),NEX)
      NEX = NEX + 1
 720  CONTINUE
         ENDIF
C ... ARSM
      IF(ITURB >= 10.AND.ITURB <= 19 .OR.ISTRES > 0) THEN
      DO NS = 1,6
        CALL IMPVO2(S11(IG1,NS),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),NEX)
      NEX = NEX + 1
      ENDDO
      DO NS = 1,5
        CALL IMPVO2(BIJ(IG1,NS),ZZZ,APU,APU2,ICON(IC1),NPATCH(N),
     +       IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),ITAG(1,N),IN,JN,KN,
     +       NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),JAPP(1,N),NEX)
      NEX = NEX + 1
      ENDDO
      ENDIF ! (ITURB >= 10.AND.ITURB <= 19 .OR.ISTRES > 0)
C ... Multiphase variables
      IF(MULPHL) THEN
        DO IPHASE = 1,NPHASES

        DO I = 1,NTOT1
        F1RN(I) = PRO(IG1+I-1)%TEMP(IPHASE) ! Single-precision at the moment
        F1RW(I) = VAR(IG1+I-1)%ALFA(IPHASE)
        F1E(I)  = VAR(IG1+I-1)%X(IPHASE)
        ENDDO

        CALL IMPVO2(F1RN,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1
        CALL IMPVO2(F1RW,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1
        CALL IMPVO2(F1E,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1

        DO I = 1,NTOT1
        PRO(IG1+I-1)%TEMP(IPHASE) = F1RN(I)
        VAR(IG1+I-1)%ALFA(IPHASE) = F1RW(I)
        VAR(IG1+I-1)%X(IPHASE)    = F1E(I)
        ENDDO

        IF(TWO_FLUIDL) THEN

        DO I = 1,NTOT1
        F1RN(I) = VAR(IG1+I-1)%U(IPHASE)
        F1RW(I) = VAR(IG1+I-1)%V(IPHASE)
        F1E(I)  = VAR(IG1+I-1)%W(IPHASE)
        ENDDO

        CALL IMPVO2(F1RN,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1
        CALL IMPVO2(F1RW,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1
        CALL IMPVO2(F1E,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1

        DO I = 1,NTOT1
        VAR(IG1+I-1)%U(IPHASE) = F1RN(I)
        VAR(IG1+I-1)%V(IPHASE) = F1RW(I)
        VAR(IG1+I-1)%W(IPHASE) = F1E(I)
        ENDDO

        ENDIF ! MULPHC

       ENDDO ! NPHASES
      ENDIF ! MULPHL
C ... Intermittency variables
      IF(TRANSL) THEN

        DO I = 1,NTOT1
           F1RN(I) = TRM(IG1+I-1)%G
           F1RW(I) = TRM(IG1+I-1)%RET
        ENDDO

        CALL IMPVO2(F1RN,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1
        CALL IMPVO2(F1RW,ZZZ,APU,APU2,ICON(IC1),
     +       NPATCH(N),IW(1,1,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +       ITAG(1,N),IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP(1,N),
     +       JAPP(1,N),NEX)
        NEX = NEX + 1

        DO I = 1,NTOT1
           TRM(IG1+I-1)%G   = F1RN(I)
           TRM(IG1+I-1)%RET = F1RW(I)
        ENDDO

      ENDIF ! TRANSL

      ENDIF   ! (M <=  MGRID(N))
8000  CONTINUE

      RETURN
      END SUBROUTINE CONNEC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERNEC(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,PRO,VAR,TRM,
     +                  XC,YC,ZC,MAXB,M,N1,N2,VIST,EPS2)

      USE TYPE_ARRAYS
      USE MPI
      USE INTEGERS,    ONLY : NPPV,NCON,IPRO,MAXW,NB 
      USE MAIN_ARRAYS, ONLY : ICON,ITAG,IWAPP,JAPP,NCPAT,NPATCH,IW,
     & MGRID,IG,IC,IR,IQ,JF,IMAX,JMAX,KMAX,ZZZ,NPROCE,A1,A2,A3,A1XA,
     & A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,BLKS

      USE NS3CO            

      IMPLICIT NONE

      INTEGER :: MAXB, M, N1, N2, I, N, IG2, IR2, IQ2, IF2, IC2, IM2,
     &           ISS, IERR, NGL, IT2

      REAL :: RO(*),TEMP(*),U(*),V(*),W(*),E(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),VIST(*),EPS2(*)
      REAL :: BAL(14),BALS(14)

      REAL :: XC(*), YC(*), ZC(*)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

C ... BAL is balance table
C ... what                       INLET I  OUTLET I  Variable
C ... mass flux                     1        8       XMASS
C ... mass flux                     2        9       XMAST
C ... P*A in computational domain   3       10       PAVE             
C ... P*A in ghost cells            4       11       PAVG    
C ... area                          5       12       AREA 
C ... e*A                           6       13       EAVE 
C ... e*A in computational domain   7       14       EAVC             

      IF(PERCHL) THEN

      DO I = 1,14
         BAL(I) = 0.
      ENDDO

      IF(M /= 1) STOP 'PERNEC M /= 1'

      DO N = N1,N2
      IG2     = IG(1,N)
      IR2     = IR(1,N)
      IQ2     = IQ(1,N)
      IF2     = JF(1,N)
      IC2     = IC(1,N)
      ISS     = 1
      IM2     = 1
      IT2     = 1
      IF (STRESL) ISS = IG2
      IF (MULPHL) IM2 = IG2
      IF(TRANSL)  IT2 = IG2

      NGL     = NPROCE(1+N,IPRO) !Global block number
             
C ... For circulating inlet only
             
      CALL PERNE3(RO(IG2),TEMP(IG2),U(IG2),V(IG2),W(IG2),E(IG2),
     +     P(IG2),PDIFF(IG2),RK(IR2),REPS(IR2),PRO(IM2),VAR(IM2),
     +     TRM(IT2),XC(IG2),YC(IG2),ZC(IG2),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,FRSVEL,FRSSIE,FRSTEM,
     +     A1(IG2),A2(IG2),A3(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +     A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +     N,NPATCH(N),ICON(IC2),ITURB,MAXB,BAL,VIST(IG2),EPS2(IG2),
     +     MULPHL,TRANSL,BLKS(NGL)%SOLUTION_TYPE)

      DO I = 1,14
         BALS(I) = 0.
      ENDDO

C ... Update periodic boundaries

      CALL PERNE1(RO(IG2),TEMP(IG2),U(IG2),V(IG2),W(IG2),E(IG2),
     +     P(IG2),PDIFF(IG2),RK(IR2),REPS(IR2),PRO(IM2),VAR(IM2),
     +     TRM(IT2),XC(IG2),YC(IG2),ZC(IG2),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),IN,JN,KN,FRSDEN,FRSVEL,ALPHA,BETA,N,NPATCH(N),
     +     ICON(IC2),A1(IG2),A2(IG2),A3(IG2),A1XA(IG2),A1YA(IG2),
     +     A1ZA(IG2),A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3XA(IG2),A3YA(IG2),
     +     A3ZA(IG2),ITURB,MAXB,BALS,BLKS(NGL)%INLRC,BLKS(NGL)%OUTRC,
     +     BLKS(NGL)%SOLUTION_TYPE)
      DO I = 1,14
         BAL(I) = BAL(I) + BALS(I)
      ENDDO
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(BAL,BALS,14,MPI_REAL8,MPI_SUM,
     &         MPI_COMM_WORLD,IERR)
         DO I = 1,14
            BAL(I) = BALS(I)
         ENDDO
      ENDIF

      IF(IPRO == 1) THEN ! Mass flows in and out (all in BOUFLO 14.2.2011)
         QMFIN  = QMFIN  !- BAL(1)
         QMFOUT = QMFOUT !+ BAL(8)
         QMEIN  = QMEIN  !- BAL(7)
         QMEOUT = QMEOUT !+ BAL(14)
      ENDIF

      DO N = N1,N2
      IG2     = IG(1,N)
      IR2     = IR(1,N)
      IQ2     = IQ(1,N)
      IF2     = JF(1,N)
      IC2     = IC(1,N)
      ISS     = 1
      IM2     = 1
      IF (STRESL) ISS = IG2
      IF (MULPHL) IM2 = IG2
      IF (TRANSL) IT2 = IG2

C ... Update periodic boundaries

      CALL PERNE2(RO(IG2),TEMP(IG2),U(IG2),V(IG2),W(IG2),E(IG2),
     +     P(IG2),PDIFF(IG2),RK(IR2),REPS(IR2),PRO(IM2),VAR(IM2),
     +     TRM(IT2),XC(IG2),YC(IG2),ZC(IG2),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +     GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,FRSVEL,FRSSIE,FRSTEM,
     +     A1(IG2),A2(IG2),A3(IG2),A1XA(IG2),A1YA(IG2),A1ZA(IG2),
     +     A2XA(IG2),A2YA(IG2),A2ZA(IG2),A3XA(IG2),A3YA(IG2),A3ZA(IG2),
     +     N,NPATCH(N),ICON(IC2),ITURB,MAXB,BAL,BLKS(NGL)%INLRC,
     +     BLKS(NGL)%OUTRC,BLKS(NGL)%SOLUTION_TYPE)
      ENDDO

      ENDIF
 
      RETURN
      END SUBROUTINE PERNEC
C
C ----------------------------------------------------------------------
C --- Subroutine for handling mirror boundary conditions ------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,
     +     PRO,VAR,TRM,MAXSB,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     M,N1,N2,XVEL,ISCL,PRC)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IR,IQ,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,IWAVEB,BLKS,NPROCE,NTOT
      USE INTEGERS, ONLY : IPRO
      USE TYPE_ARRAYS
      USE NS3CO           

      IMPLICIT NONE

      INTEGER :: MAXSB,M,N1,N2,ISCL,N,IG1,IR1,IQ1,IC1,IP1,IF1,IT1,NGL,
     +           IPC1,I,II,JJ,KK

      REAL :: XVEL

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),TTS(*),FI(MAXSB,*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(PRE_COR)          :: PRC(*)

C
C ... THIS SUBROUTINE MIRRORS ALL THE PATCHES
C
      DO N = N1,N2           ! BLOCK LOOP 1

         IF(M <= MGRID(N)) THEN

            IG1 = IG(M,N)
            IR1 = IR(M,N)
            IQ1 = IQ(M,N)
            IC1 = IC(M,N)
            IP1 = 1
            IF1 = 1
            IT1 = 1
            IPC1= 1
            NGL = NPROCE(N+1,IPRO)
            IF(BLKS(NGL)%IPRESC == 1) IPC1= IG1

            IF(MULPHL) IP1 = IG1
            IF(FRESUL) IF1 = IG1
            IF(TRANSL) IT1 = IG1

            CALL MIRVOL(RO(IG1),U(IG1),V(IG1),W(IG1),E(IG1),EPS2(IG1),
     +           VIST(IG1),P(IG1),PDIFF(IG1),RK(IR1),REPS(IR1),TTS(IR1),
     +           FI(IQ1,1),PRO(IP1),VAR(IP1),TRM(IT1),MAXSB,NSCAL,ITURB,
     +           MULPHL,TRANSL,A1X(IG1),A1Y(IG1),A1Z(IG1),
     +           A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     +           ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +           IN,JN,KN,XVEL,ISCL,IWAVEB(IF1),FRESUL,M,
     +           BLKS(NGL)%SOLUTION_TYPE,PRC(IPC1),BLKS(NGL)%IPRESC)

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE MIR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRC(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,M,N1,N2)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,IW
      USE NS3CO            

      IMPLICIT NONE

      INTEGER :: M,N1,N2,N,IG1,IC1

      REAL :: VOL(*),D1(*),D2(*),D3(*),DISTW(*),BLANK(*),
     +        A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)
      REAL :: XC(*),YC(*),ZC(*)

C
C ... THIS SUBROUTINE MIRRORS COORDINATES AT THE CENTERPOINTS
C ... AND ALSO FOR A FREE SURFACE
C

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IC1     = IC(M,N)

      CALL MIRCOR(XC(IG1),YC(IG1),ZC(IG1),VOL(IG1),D1(IG1),D2(IG1),
     +     D3(IG1),DISTW(IG1),BLANK(IG1),A1(IG1),A2(IG1),A3(IG1),
     +     A1X(IG1),A1Y(IG1),A1Z(IG1),A2X(IG1),A2Y(IG1),A2Z(IG1),
     +     A3X(IG1),A3Y(IG1),A3Z(IG1),
     +     ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +     IN,JN,KN,M)
      ENDIF
7900  CONTINUE

      RETURN
      END SUBROUTINE MIRC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRC2(XCO,YCO,ZCO,M,N1,N2)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,IW
      USE NS3CO

      IMPLICIT NONE

      INTEGER :: M, N1, N2, N, IG1, IC1
      REAL :: XCO(*), YCO(*), ZCO(*)

C
C ... THIS SUBROUTINE MIRRORS CORNER POINTS OF ALL THE PATCHES
C ... AND ALSO FOR A FREE SURFACE

      DO N = N1,N2           ! BLOCK LOOP 1

         IF(M <= MGRID(N)) THEN

            IG1 = IG(M,N)
            IC1 = IC(M,N)

            CALL MIRCCC(XCO(IG1),YCO(IG1),ZCO(IG1),ICON(IC1),NPATCH(N),
     +           IW(1,N,M),IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,M)
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE MIRC2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRR(UU,UV,UW,VV,VW,WW,A1X,A1Y,A1Z,
     +     A2X,A2Y,A2Z,A3X,A3Y,A3Z,M,N1,N2,MAXSB)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,F1R,F1RM,F1RN,F1RW,F1E,F1RK,NTOT
      USE NS3CO

      IMPLICIT NONE

      INTEGER :: M, N1, N2, MAXSB, N, IG1, IC1, ii
      REAL :: UU(*) ,VV(*) ,WW(*)
      REAL :: UV(*) ,UW(*) ,VW(*)
      REAL :: A1X(*),A1Y(*),A1Z(*)
      REAL :: A2X(*),A2Y(*),A2Z(*)
      REAL :: A3X(*),A3Y(*),A3Z(*)

C
C ... THIS SUBROUTINE MIRRORS THE REYNOLDS STRESSES FOR ALL THE PATCHES
C
      DO N = N1,N2           ! BLOCK LOOP 1

         IF(M <= MGRID(N)) THEN

            IG1 = IG(M,N)
            IC1 = IC(M,N)

            CALL MIRREL(UU(IG1),UV(IG1),UW(IG1),VV(IG1),VW(IG1),WW(IG1),
     +           F1R,F1RM,F1RN,F1RW,F1E,F1RK,A1X(IG1),A1Y(IG1),A1Z(IG1),
     +           A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     +           ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +           IN,JN,KN,MAXSB,M) 

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE MIRR
C
C ----------------------------------------------------------------------
C --- Subroutine for handling cyclic boundary conditions ------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYC(U,V,W,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,M,N1,N2,
     +     F1R,F1RM,VAR)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,JLOC,JTRA,APP,AAA,IW,IT,IL,IK
      USE NS3CO
      USE INTEGERS,    ONLY : MAXB

      USE TYPE_ARRAYS         

      IMPLICIT NONE

      INTEGER M,N1,N2,MAXSB,N,IG1,IC1,IS,IPHASE
      REAL U(*),V(*),W(*),A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),
     + A3X(*),A3Y(*),A3Z(*),F1R(*),F1RM(*)
      TYPE(MPHASE_VARIABLES) VAR(*)

 
C ... CYCLIC BCS FOR ALL PATCHES

      DO IPHASE = 1,NPHASE

      IF(TWO_FLUIDL) THEN
      U(1:MAXB) = VAR(1:MAXB)%U(IPHASE)
      V(1:MAXB) = VAR(1:MAXB)%V(IPHASE)
      W(1:MAXB) = VAR(1:MAXB)%W(IPHASE)
      ENDIF

      CALL CONEC(U,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)
      CALL CONEC(V,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)
      CALL CONEC(W,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)

      IS = 1                     ! CORRECT ADDRESS FOR AAA-ARRAY
      DO 7800 N = 1,N1-1
7800  IS = IS + NPATCH(N)

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IC1     = IC(M,N)

      CALL CYCVOL(U(IG1),V(IG1),W(IG1),A1X(IG1),A1Y(IG1),A1Z(IG1),
     + A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     + AAA(1,IS),ICON(IC1),NPATCH(N),IW(1,N,M),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,M,N,IT(M,N),IL(M,N),
     + IK(M,N))  ! FOR PRINTING. POIS JOSKUS
      ENDIF
      IS = IS + NPATCH(N)
7900  CONTINUE

      IF(TWO_FLUIDL) THEN
      VAR(1:MAXB)%U(IPHASE) = U(1:MAXB)
      VAR(1:MAXB)%V(IPHASE) = V(1:MAXB)
      VAR(1:MAXB)%W(IPHASE) = W(1:MAXB)
      ENDIF

      ENDDO  ! NPHASE

      RETURN
      END SUBROUTINE CYC
C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface boundary conditions ---------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYCV(U,V,W,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,M,N1,N2,
     +     F1R,F1RM,VAR)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     + NPATCH,JLOC,JTRA,APP,AAA,IW,IT,IL,IK
      USE NS3CO
      USE INTEGERS,    ONLY : MAXB

      USE TYPE_ARRAYS         

      IMPLICIT NONE

      INTEGER M,N1,N2,MAXSB,N,IG1,IC1,IS,IPHASE
      REAL U(*),V(*),W(*),A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),
     + A3X(*),A3Y(*),A3Z(*),F1R(*),F1RM(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
 
C ... CYCLIC BCS FOR ALL PATCHES FOR OTHER VARIABLES THAN U, V and W

      CALL CONEC(U,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)
      CALL CONEC(V,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)
      CALL CONEC(W,F1R,F1RM,JLOC,JTRA,APP,M,N1,N2,0,2)

      IS = 1                     ! CORRECT ADDRESS FOR AAA-ARRAY
      DO 7800 N = 1,N1-1
7800  IS = IS + NPATCH(N)

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IC1     = IC(M,N)

      CALL CYCVOL(U(IG1),V(IG1),W(IG1),A1X(IG1),A1Y(IG1),A1Z(IG1),
     + A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     + AAA(1,IS),ICON(IC1),NPATCH(N),IW(1,N,M),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,M,N,IT(M,N),IL(M,N),
     + IK(M,N))  ! FOR PRINTING. POIS JOSKUS
      ENDIF
      IS = IS + NPATCH(N)
7900  CONTINUE

      RETURN
      END SUBROUTINE CYCV
C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface boundary conditions ---------
C ----------------------------------------------------------------------
C
      SUBROUTINE FRE(M,N1,N2,XVEL,ISUR)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,NPATCH,
     +     JLOC,JTRA,APP,AAA,IW,IT,IL,IK,IR,IQ,HFLUX,XCO,YCO,ZCO,
     +     TEMP,IHF,VOL,A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     +     A3XA,A3YA,A3ZA,NPROCE,IGRID,XC,YC,ZC,RO,U,V,W,E,EPS2,
     +     VIST,P,PDIFF,RK,REPS,FI,BIJ,WAVEH,POLD,DSURLE,IWAVEB

      USE NS3CO,       ONLY : IN,JN,KN,ISTRES,NSCAL,ITURB,FRSDEN,
     +     FRSVEL,CHLREF,ICYCLE,CFL,FRSPRE,GX,GY,GZ,GROUND

      USE INTEGERS,    ONLY : IB,MAXB,IPRO,MAXSB,MAXEB         

      IMPLICIT NONE

      INTEGER M,N1,N2,ISUR,N,IG1,IR1,IQ1,IC1,IE1,NGL 

      REAL XVEL

C ... THIS SUBROUTINE HANDLES FREE SURFACE PATCHES

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IR1     = IR(M,N)
      IQ1     = IQ(M,N)
      IC1     = IC(M,N)
      IE1     = 1
      IF(ISTRES > 0) IE1 = IG1
      NGL     = NPROCE(1+N,IPRO)
        
      CALL FREVOL(RO(IG1),U(IG1),V(IG1),W(IG1),E(IG1),EPS2(IG1),
     +     VIST(IG1),P(IG1),PDIFF(IG1),RK(IR1),REPS(IR1),FI(IQ1,1),
     +     BIJ(IE1,1),MAXEB,MAXSB,NSCAL,ITURB,A1XA(IG1),A1YA(IG1),
     +     A1ZA(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),
     +     A3ZA(IG1),XC,YC,ZC,FRSDEN,WAVEH(1),IHF(1,M),
     +     M,ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,
     +     XVEL,FRSPRE,ISTRES,POLD(IG1),DSURLE(1),IWAVEB(IG1),NGL)

      IF(ISUR == 1) THEN ! Update the free-surface here instead of SURGRI
      CALL SURSET(XCO(IG1),YCO(IG1),ZCO(IG1),XC(IG1),YC(IG1),ZC(IG1),
     +     U(IG1),V(IG1),W(IG1),TEMP(IG1),PDIFF(IG1),
     +     FRSDEN,FRSVEL,CHLREF,ICYCLE,IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +     IN,JN,KN,HFLUX,N,M,NPATCH(N),ICON(IC1),IHF(1,M),
     +     VOL(IG1),A1(IG1),A2(IG1),A3(IG1),A1XA(IG1),A1YA(IG1),
     +     A1ZA(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),
     +     A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     +     IB,MAXB,IGRID((NPROCE(1+N,IPRO)),2),CFL)
      ENDIF
      ENDIF
7900  CONTINUE

      RETURN
      END SUBROUTINE FRE
C
C ----------------------------------------------------------------------
C --- Subroutine for handling sliding mesh boundary conditions ---------
C ----------------------------------------------------------------------
C
      SUBROUTINE SLD(U,V,W,EPS2,VIST,RO,TEMP,P,PDIFF,PTUR,E,RM,RN,RW,
     +     RK,REPS,FI,S11,BIJ,TIJ,PRO,VAR,TRM,F1R,F1RM,F1RN,
     +     MAXB,MAXEB,MAXSB,MAXSS,
     +     XCO,YCO,ZCO,M,N1,N2,ISTRES,D1,D2,D3,VOL,DM,DN,DW,ISLIM)

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IC,ICON,IMAX,JMAX,KMAX,
     +     NPATCH,JLOC,JTRA,APP,AAA,IW,IT,IL,IK,IR,IQ,JF
      USE NS3CO,       ONLY : STRESL,NSCAL,ITURB,TIMEL,MULPHL,TRANSL
      USE INTEGERS,    ONLY : MAXW,MBPRO,IPRO        
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: M,N1,N2,MAXB,MAXEB,MAXSB,MAXSS,ISTRES,N,IG1,IR1,IQ1,
     &           IF1,IC1,IE1,IQQ1,IP1,ISLIM,IT1
     
      REAL :: U(*),V(*),W(*),EPS2(*),VIST(*),RO(*),TEMP(*),
     &        P(*),PDIFF(*),
     &        PTUR(*),E(*),RM(*),RN(*),RW(*),F1R(*),F1RM(*),F1RN(*),
     &        RK(*),REPS(*),FI(MAXSB,*),S11(MAXSS,*),BIJ(MAXEB,*),
     &        TIJ(*),D1(*),D2(*),D3(*),VOL(*),DM(*),DN(*),DW(*)

      REAL :: XCO(*), YCO(*), ZCO(*)

      TYPE(PROPERTIES)       PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(INTERMITTENCY)    TRM(*)

      DO 7900 N = N1,N2      

      IF(M <= MGRID(N)) THEN

         IG1 = IG(1,N)
         IR1 = IR(1,N)
         IQ1 = IQ(1,N)
         IF1 = JF(1,N)
         IC1 = IC(1,N)
         IE1 = 1

         IF(ISTRES > 0) IE1 = IG1
         IQQ1 = 1

         IF(STRESL) IQQ1 = IG1
         IP1  = 1

         IF(MULPHL) IP1 = IG1
         IT1  = 1

         IF(TRANSL) IT1 = IG1

         CALL SLIMES(U(IG1),V(IG1),W(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +       N,NPATCH(N),ICON(IC1),RO(IG1),TEMP(IG1),P(IG1),PDIFF(IG1),
     +       PTUR(IG1),E(IG1),RM(IG1),RN(IG1),RW(IG1),RK(IG1),REPS(IG1),
     +       EPS2(IG1),VIST(IG1),FI(IQ1,1),S11(IQQ1,1),BIJ(IE1,1),
     +       TIJ(IG1),PRO(IP1),VAR(IP1),TRM(IT1),MAXB,MAXEB,MAXSB,
     +       MAXSS,NSCAL,XCO(IG1),YCO(IG1),ZCO(IG1),M,ITURB,
     +       IPRO,MBPRO,TIMEL,MAXW,ISTRES,D1(IG1),D2(IG1),D3(IG1),
     +       VOL(IG1),DM(IG1),DN(IG1),DW(IG1),ISLIM,MULPHL,TRANSL)

         ENDIF

7900  CONTINUE

      RETURN
      END SUBROUTINE SLD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NEWTIM

      USE TYPE_ARRAYS

      USE NS3CO, ONLY : NBLOCK, TIMEL, COORL, ITURB, NSCAL, RKLIM,
     &    EPSLIM, IPRESC, IN, JN, KN, TWO_FLUIDL

      USE MAIN_ARRAYS, ONLY :  NTOT, IG, IMAX, JMAX, KMAX, MGRID, IGRID,
     &    RO, RM, RN, RW, E, RK, REPS, FI, TEMP, U, V, W, PDIFF,
     &    TOLD, UOLD, VOLD, WOLD, POLD, RKOLD, EPSOLD, FIOLD,
     &    XCO, YCO, ZCO,
     &    ROLE2, RMLE2, RNLE2, RWLE2, ELE2, RKLE2, EPSLE2, FILE2,
     &    ROLE3, RMLE3, RNLE3, RWLE3, ELE3, RKLE3, EPSLE3, FILE3,
     &    XLE2, YLE2, ZLE2, XLE3, YLE3, ZLE3, 
     &    VOL, RLOLIM, UPPLIM, NPROCE, BLKS, VAR, F1R, F1RM, PLE2, PLE3

      USE INTEGERS, ONLY : MGM, MAXSB, IPRO, MBPRO 

      IMPLICIT NONE

      INTEGER :: N, NTOT1, IG1, ISTRID, JSTRID, KSTRID, 
     &           NGL, I, J, K, KA, JJ, KK, NS, L,
     &           NBG, M, IG2, IPHASE, NPHASE, NTOT2, IGN1, IGN2
      REAL :: R1, R2, R3, R1R, R2R, R3R, RON, RMN, RNN, RWN, EN, 
     &        RKN, REN, RONA, RONE, PN

C
C ... INITIALIZE ARRAYS FOR PREVIOUS TIME LEVELS
C
      
C ... zero order initial guess
      R1 = 1.
      R2 = 0.
      R3 = 0.
C ... first order initial guess
c      R1 = 2.
c      R2 = -1.
c      R3 = 0.
C ... second order initial guess is always unstable!
c      R1 = 3.
c      R2 = -3.
c      R3 = 1.
c      R1  = 1.5
c      R2  = -.5
c      R3  = 0.

      R1R = R1
      R2R = R2
      R3R = R3

c      IF(IPRESC >= 1) THEN ! Vanha paikka
c         R1R = 1.
c         R2R = 0.
c         R3R = 0.
c      ENDIF

      DO 9000 N = 1,NBLOCK
      NTOT1   = NTOT(1,N)
      IG1     = IG(1,N)
      ISTRID  = IMAX(1,N) + 2*IN
      JSTRID  = JMAX(1,N) + 2*JN
      KSTRID  = ISTRID*JSTRID
      NGL     = NPROCE(N+1,IPRO)

C ... Blockwise properties

      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      IF(BLKS(NGL)%IPRESC >= 1) THEN
         R1R = 1.
         R2R = 0.
         R3R = 0.
      ENDIF

      NPHASE = BLKS(NGL)%NPHASE

      IF(TIMEL) THEN ! IN CASE ONLY COORL IS ACTIVATED
         DO K = 1,KMAX(1,N)
            KA      = (KN+K-1)*KSTRID
            DO J = 1,JMAX(1,N)
               JJ      = (JN+J-1)*ISTRID + IN + KA
               DO I = 1,IMAX(1,N)
                  KK      = JJ + I + IG1 - 1
                  RON = R1R*RO(KK)+ R2R*ROLE2(KK)+ R3R*ROLE3(KK)
                  RMN = R1*RM(KK) + R2*RMLE2(KK) + R3*RMLE3(KK)
                  RNN = R1*RN(KK) + R2*RNLE2(KK) + R3*RNLE3(KK)
                  RWN = R1*RW(KK) + R2*RWLE2(KK) + R3*RWLE3(KK)
                  EN  = R1*E(KK)  + R2*ELE2(KK)  + R3*ELE3(KK)
                  PN  = R1*PDIFF(KK) + R2*PLE2(KK)  + R3*PLE3(KK)
                  ROLE3(KK) = ROLE2(KK)
                  RMLE3(KK) = RMLE2(KK)
                  RNLE3(KK) = RNLE2(KK)
                  RWLE3(KK) = RWLE2(KK)
                  ELE3(KK)  = ELE2(KK)
                  PLE3(KK)  = PLE2(KK)
                  ROLE2(KK) = RO(KK)
                  RMLE2(KK) = RM(KK)
                  RNLE2(KK) = RN(KK)
                  RWLE2(KK) = RW(KK)
                  ELE2(KK)  = E(KK)
                  PLE2(KK)  = PDIFF(KK)
                  RO(KK)    = RON
                  RM(KK)    = RMN
                  RN(KK)    = RNN
                  RW(KK)    = RWN
                  E(KK)     = EN
                  IF(NPHASE > 1) THEN ! (Else meaningless)
                  DO IPHASE = 1,NPHASE
                  RONA      = R1R * VAR(KK)%ARO(IPHASE)    +
     &                        R2R * VAR(KK)%AROLE2(IPHASE) +
     &                        R3R * VAR(KK)%AROLE3(IPHASE)
                  RONE      = R1R * VAR(KK)%ARE(IPHASE)    +
     &                        R2R * VAR(KK)%ARELE2(IPHASE) +
     &                        R3R * VAR(KK)%ARELE3(IPHASE)
                  VAR(KK)%AROLE3(IPHASE) = VAR(KK)%AROLE2(IPHASE)
                  VAR(KK)%AROLE2(IPHASE) = VAR(KK)%ARO(IPHASE)
                  VAR(KK)%ARO(IPHASE)    = RONA
                  VAR(KK)%ARELE3(IPHASE) = VAR(KK)%ARELE2(IPHASE)
                  VAR(KK)%ARELE2(IPHASE) = VAR(KK)%ARE(IPHASE)
                  VAR(KK)%ARE(IPHASE)    = RONE
                  ENDDO
                  IF(TWO_FLUIDL) THEN
                  DO IPHASE = 1,NPHASE
                  RMN = R1*VAR(KK)%ARM(IPHASE) + 
     &            R2*VAR(KK)%ARMLE2(IPHASE) + R3*VAR(KK)%ARMLE3(IPHASE)
                  RNN = R1*VAR(KK)%ARN(IPHASE) + 
     &            R2*VAR(KK)%ARNLE2(IPHASE) + R3*VAR(KK)%ARNLE3(IPHASE)
                  RWN = R1*VAR(KK)%ARW(IPHASE) + 
     &            R2*VAR(KK)%ARWLE2(IPHASE) + R3*VAR(KK)%ARWLE3(IPHASE)
                  VAR(KK)%ARMLE3(IPHASE) = VAR(KK)%ARMLE2(IPHASE)
                  VAR(KK)%ARMLE2(IPHASE) = VAR(KK)%ARM(IPHASE)
                  VAR(KK)%ARM(IPHASE)    = RMN
                  VAR(KK)%ARNLE3(IPHASE) = VAR(KK)%ARNLE2(IPHASE)
                  VAR(KK)%ARNLE2(IPHASE) = VAR(KK)%ARN(IPHASE)
                  VAR(KK)%ARN(IPHASE)    = RNN
                  VAR(KK)%ARWLE3(IPHASE) = VAR(KK)%ARWLE2(IPHASE)
                  VAR(KK)%ARWLE2(IPHASE) = VAR(KK)%ARW(IPHASE)
                  VAR(KK)%ARW(IPHASE)    = RWN
                  ENDDO ! IPHASE
                  ENDIF ! TWO_FLUIDL
                  ENDIF ! NPHASE
                  U(KK)     = RMN/RON ! velocities are guest
                  V(KK)     = RNN/RON ! velocities are guest
                  W(KK)     = RWN/RON ! velocities are guest
                  TOLD(KK)  = TEMP(KK) 
                  UOLD(KK)  = U(KK) 
                  VOLD(KK)  = V(KK) 
                  WOLD(KK)  = W(KK) 
                  POLD(KK)  = PDIFF(KK)
C ... but for pressure and temperature old values are used for now 
C ... PPR 30.12.99
               ENDDO
            ENDDO
         ENDDO

c         DO 3000 L   = 1,NTOT1
c            KK        = L + IG1 - 1
c            TOLD(KK)  = TEMP(KK) 
c            UOLD(KK)  = U(KK) 
c            VOLD(KK)  = V(KK) 
c            WOLD(KK)  = W(KK) 
c            POLD(KK)  = PDIFF(KK)
c 3000    CONTINUE

         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            DO K = 1,KMAX(1,N)
               KA      = (KN+K-1)*KSTRID
               DO J = 1,JMAX(1,N)
                  JJ      = (JN+J-1)*ISTRID + IN + KA
                  DO I = 1,IMAX(1,N)
                     KK      = JJ + I + IG1 - 1
                     RKN = R1*RK(KK)   + R2* RKLE2(KK) + R3* RKLE3(KK)
                     REN = R1*REPS(KK) + R2*EPSLE2(KK) + R3*EPSLE3(KK)
                     RKLE3(KK)  = RKLE2(KK)
                     EPSLE3(KK) = EPSLE2(KK)
                     RKLE2(KK)  = RK(KK)
                     EPSLE2(KK) = REPS(KK)
                     RK(KK)     = RKN
                     REPS(KK)   = REN
                     RKOLD(I)   = RKN
                     EPSOLD(I)  = REN
C ... And now also for the turbulence quantities!!
                  ENDDO
               ENDDO
            ENDDO
            
c            DO 3100 L    = 1,NTOT1
c               I         = L + IG1 - 1
c               RKOLD(I)  = RKN
c               EPSOLD(I) = REN
c 3100       CONTINUE

            CALL LOLIM(RK(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +           RKLIM)
            CALL LOLIM(REPS(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +           EPSLIM)
         ENDIF                  ! ITURB
         
         IF (NSCAL /= 0) THEN
            DO 3200 NS     = 1,NSCAL
            DO K = 1,KMAX(1,N)
               KA      = (KN+K-1)*KSTRID
               DO J = 1,JMAX(1,N)
                  JJ      = (JN+J-1)*ISTRID + IN + KA
                  DO I = 1,IMAX(1,N)
                     KK      = JJ + I + IG1 - 1
           RKN = R1*FI(KK,NS)   + R2*FILE2(KK,NS) + R3*FILE3(KK,NS)
                     FILE3(KK,NS) = FILE2(KK,NS)
                     FILE2(KK,NS) = FI(KK,NS)
                     FI(KK,NS)    = RKN
                  ENDDO
               ENDDO
            ENDDO
            DO 3200 L      = 1,NTOT1
               I           = L + IG1 - 1
               FIOLD(I,NS) = FI(I,NS)
 3200       CONTINUE
            CALL LOLISC(FI,IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +           RLOLIM,UPPLIM,MAXSB,NSCAL)
         ENDIF                  ! SCALARS
      ENDIF                     ! TIMEL

      NBG     = NPROCE(1+N,IPRO)               ! global block number
      IF (COORL .AND. IGRID(NBG,1) >= 11) THEN ! not with rotation grid
         DO L  = 1,NTOT1
            I       = L + IG1 - 1
            XLE3(I) = XLE2(I)
            YLE3(I) = YLE2(I)
            ZLE3(I) = ZLE2(I)
            
            XLE2(I) = XCO(I)
            YLE2(I) = YCO(I)
            ZLE2(I) = ZCO(I)
         ENDDO

      ELSEIF(COORL .AND. IGRID(NBG,1) < 11) THEN ! Is this needed??

         DO L  = 1,NTOT1
            I       = L + IG1 - 1
            XLE3(I) = XLE2(I)
            YLE3(I) = YLE2(I)
            ZLE3(I) = ZLE2(I)

            XLE2(I) = XCO(I)
            YLE2(I) = YCO(I)
            ZLE2(I) = ZCO(I)
         ENDDO
      ENDIF ! COORL


C ... update for multigrid acceleration
      IF(TIMEL) THEN ! IN CASE ONLY COORL IS ACTIVATED
      DO 7000 M = 2,MGRID(N)
      IG1     = IG(M-1,N)
      IG2     = IG(M,N)
      NTOT2   = NTOT(M,N)
      IGN1    = IG1 + NTOT1 - 1
      IGN2    = IG2 + NTOT2 - 1
      CALL TRANTI(ROLE2(IG2),RMLE2(IG2),RNLE2(IG2),RWLE2(IG2),ELE2(IG2),
     1            ROLE2(IG1),RMLE2(IG1),RNLE2(IG1),RWLE2(IG1),ELE2(IG1),
     2     VOL(IG2),VOL(IG1),
     3      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))
      CALL TRANTI(ROLE3(IG2),RMLE3(IG2),RNLE3(IG2),RWLE3(IG2),ELE3(IG2),
     1            ROLE3(IG1),RMLE3(IG1),RNLE3(IG1),RWLE3(IG1),ELE3(IG1),
     2     VOL(IG2),VOL(IG1),
     3      IMAX(M,N),JMAX(M,N),KMAX(M,N),KMAX(M-1,N))
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
       CALL TRANTK(RKLE2(IG2),EPSLE2(IG2),RKLE2(IG1),EPSLE2(IG1),
     1      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     2      KMAX(M-1,N))
       CALL TRANTK(RKLE3(IG2),EPSLE3(IG2),RKLE3(IG1),EPSLE3(IG1),
     1      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     2      KMAX(M-1,N))
      ENDIF
      IF(NPHASE > 1) THEN ! Multiphase variables
      DO IPHASE = 1,NPHASE
       F1R(IG1:IGN1)  = VAR(IG1:IGN1)%AROLE2(IPHASE)
       F1RM(IG1:IGN1) = VAR(IG1:IGN1)%AROLE3(IPHASE)
       CALL TRANTS(F1R(IG2),F1RM(IG2),F1R(IG1),F1RM(IG1),
     1      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     2      KMAX(M-1,N),MAXSB,NSCAL)
       VAR(IG2:IGN2)%AROLE2(IPHASE) = F1R(IG2:IGN2)
       VAR(IG2:IGN2)%AROLE3(IPHASE) = F1RM(IG2:IGN2)
       F1R(IG1:IGN1)  = VAR(IG1:IGN1)%ARELE2(IPHASE)
       F1RM(IG1:IGN1) = VAR(IG1:IGN1)%ARELE3(IPHASE)
       CALL TRANTS(F1R(IG2),F1RM(IG2),F1R(IG1),F1RM(IG1),
     1      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     2      KMAX(M-1,N),MAXSB,NSCAL)
       VAR(IG2:IGN2)%ARELE2(IPHASE) = F1R(IG2:IGN2)
       VAR(IG2:IGN2)%ARELE3(IPHASE) = F1RM(IG2:IGN2)
      ENDDO
      ENDIF
      IF(NSCAL > 0) THEN
       CALL TRANTS(FILE2(IG2,1),FILE3(IG2,1),FILE2(IG1,1),FILE3(IG1,1),
     1      VOL(IG2),VOL(IG1),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     2      KMAX(M-1,N),MAXSB,NSCAL)
      ENDIF
 7000 CONTINUE                  ! END OF MULTIGRID LOOP
      ENDIF  ! (TIMEL)
 9000 CONTINUE                  ! END OF BLOCK LOOP 1
      RETURN
      END SUBROUTINE NEWTIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ######################################################################
      SUBROUTINE SUBIN(ICYCTI,FIRSTC)
C ######################################################################

      USE MPI
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS10
      USE TYPE_ARRAYS
      USE INTEGERS,    ONLY : IPRO,MBPRO,NPRO,MAXSB,MAXW,NBCS
      USE MAIN_ARRAYS, ONLY : DRO,RO,TEMP,TOLD,RK,DRK,RKOLD,DEPS,REPS,
     &   EPSOLD,DFI,FI,FIOLD,DM,U,UOLD,DN,V,VOLD,DW,W,WOLD,DE,P,PDIFF,
     &   POLD,E,F1R,F1E,DRDP,DRDH,RESI,ZZZ,NPROCE,JET,NTOT,IG,CP,VOLN,
     &   IMAX,JMAX,KMAX,NPATCH,ICON,IR,IQ,C,TOMEGA,CXB,CYB,CZB,CMXB,
     &   CMYB,CMZB,XC,YC,ZC,RM,RN,RW,PTUR,VOL,IT,IL,IK,BLKS,PRO,VAR,
     &   F1RM,F1RN,F1RW,A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     &   A3XA,A3YA,A3ZA,HAT3,DDEPS,EPS2,VIST,VIS,D1,D2,D3,HAT1,HAT2,
     &   IDER,VTRAN,RMAV3,RNAV3,RWAV3,TRM,ROLD,
     &   APATCH,DXB,DYB,DZB,QTB,QWB,QHB
      USE NS3CO
      USE FLIGHT,      ONLY : SINK,TRIMA

      IMPLICIT NONE

      INTEGER :: N,NTOT1,IG1,L,I,J,K,NS,IMACH1,ISTRID,JSTRID,KSTRID,
     & KSTR,III,II9,ISD,ISU,N1,ICON1,IR1,IQ1,NTO3,NGL,JJ,ICYCTI,ICY,
     & IPHASE,IGL,IP1,IIND,JIND,KIND,IT1

      REAL :: CIII,APUMA,QMEIN2,QMEOU2,EFFIC,HEADT,CLDRO,VOLTOB,CXIO,
     & CYIO,CZIO,CMXIO,CMYIO,CMZIO         

      REAL :: TTOTAL, COMTIM

      CHARACTER (LEN=10) :: SOLUTION_TYPE

      LOGICAL :: FIRSTC

C ... SUBINITIALIZE ARRAYS

      DROMAX = 0.
      CD     = 0.
      CL     = 0.
      CS     = 0.
      CX     = 0.
      CY     = 0.
      CZ     = 0.
      CMX    = 0.
      CMY    = 0.
      CMZ    = 0.
      TOM    = 0.
      CXIO   = 0.
      CYIO   = 0.
      CZIO   = 0.
      CMXIO  = 0.
      CMYIO  = 0.
      CMZIO  = 0.

      IMACH  = 0
      IXERR  = 1
      DROMAX = 1.E-20        

C ... BLOCK LOOP 1 BEGINS

      DO 9000 N = 1,NBLOCK

      NTOT1 = NTOT(1,N)
      IG1   = IG(1,N)
      NGL   = NPROCE(1+N,IPRO)
      SOLUTION_TYPE = BLKS(NGL)%SOLUTION_TYPE

      IF(ICYCLE == 0) GO TO 999

      IF(ITURB >= 3 .AND. ITURB <= 23 .AND. ITURB /= 8) THEN
         rkmax = .5*blks(ngl)%frsden*blks(ngl)%refvel**2 * 2.
         DO 1600 L = 1,NTOT1
         I       = L + IG1 - 1

C ... THIS SHOULD BE TRADITIONAL

C ... KINETIC ENERGY OF TURBULENCE CANNOT EXCEED THE TOTAL ENERGY

         RK(I)   = MIN(RK(I),MAX(0.1*E(I),RKMAX))
         RK(I)   = MIN(RK(I),RKMAX)
         DRK(I)  = RK(I) - RKOLD(I)
 1600    DEPS(I) = REPS(I) - EPSOLD(I)
      ENDIF

      IF (ITURB >= 24) THEN
      DO 1800 L  =1,NTOT1
         I       = L + IG1 - 1
         RK(I)   = .5*(FI(I,1) + FI(I,4) + FI(I,6))
         DRK(I)  = RK(I) - RKOLD(I)
         DEPS(I) = REPS(I) - EPSOLD(I)
 1800 CONTINUE
      ENDIF

      DO 1500 L = 1,NTOT1

         I = L + IG1 - 1

         DM(I)   = RO(I)*(U(I) - UOLD(I))
         DN(I)   = RO(I)*(V(I) - VOLD(I))
         DW(I)   = RO(I)*(W(I) - WOLD(I))
         DE(I)   = PDIFF(I)    - POLD(I)
         HAT1(I) = TEMP(I)     - TOLD(I)
         
         IF(ROLD(I) > EPS10) THEN
            DRO(I)  = RO(I) - ROLD(I)
         ELSE
            DRO(I)  = DRDP(I)*DE(I) + DRDH(I)*CP(I)*(TEMP(I)-TOLD(I))
         ENDIF

 1500 CONTINUE

      IF(TRANSL) THEN ! Intermittency variables
      DO 1650 L = 1,NTOT1
         I        = L + IG1 - 1
         TRM(I)%DG   = TRM(I)%G   - TRM(I)%GOLD
         TRM(I)%DRET = TRM(I)%RET - TRM(I)%RETOLD
 1650 CONTINUE
      ENDIF

C ... FOR SCALAR EQUATIONS PPR 14.2

      IF (NSCAL /= 0) THEN
      DO 1700 NS = 1,NSCAL
      DO 1700 L  = 1,NTOT1
         I           = L + IG1 - 1
         DFI(I,NS)   = FI(I,NS) - FIOLD(I,NS)
 1700 CONTINUE
      ENDIF

      IF(CONVL) THEN            ! Global test (mersu)

      NGL     = NPROCE(1+N,IPRO)! Global block number

      IF(BLKS(NGL)%CONVL) THEN  ! Blockwise test
      IF (IDRXX == 1) THEN
         CALL DROMXX(DRO(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 2) THEN
         CALL DROMXX( DM(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 3) THEN
         CALL DROMXX( DN(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 4) THEN
         CALL DROMXX( DW(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 5) THEN
         CALL DROMXX(HAT1(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 6) THEN
         CALL DROMXX( DE(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 7) THEN
         CALL DROMXX(DRK(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX == 8) THEN
         CALL DROMXX(DEPS(IG1),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))
      ELSEIF (IDRXX >= 9) THEN
         CALL DROMXX(DFI(IG1,IDRXX-8),NTOT1,N,IXERR,INERR,DROMAX,
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N))

      ENDIF ! IDRXX == 1

      IMACH1 = 0 ! Old supersonic cell counter with the ghost cells

      DO L = IG1,IG1+NTOT1-1
         CIII    = C(L)**2
         APUMA   = U(L)**2+V(L)**2+W(L)**2
         IF(APUMA > CIII+1.E-7) IMACH1 = IMACH1 + 1
      ENDDO

      ISTRID = IMAX(1,N) + 2*IN
      JSTRID = JMAX(1,N) + 2*JN
      KSTR   = ISTRID*JSTRID

C ... New supersonic cell counter without the ghost cells

      DO K = 1,KMAX(1,N)
      III  = (K+1)*KSTR + IG1 - 1
         DO J = 1,JMAX(1,N)
         II9  = III + (J+1)*ISTRID + IN
           DO I = 1,IMAX(1,N)
           L = I + II9
           CIII    = C(L)**2
           APUMA   = U(L)**2+V(L)**2+W(L)**2
           IF(APUMA > CIII+1.E-7) IMACH = IMACH + 1
           ENDDO
         ENDDO
      ENDDO

      ENDIF ! BLKS(NGL)%CONVL

      ELSE
         DROMAX = 1.
         IXERR  = 1
         INERR  = 1
         IMACH  = 0

      ENDIF ! CONVL

C ... FIRST CYCLE STARTS FROM HERE

 999  CONTINUE
      DO 2000 L = 1,NTOT1
         I         = L + IG1 - 1
         TOLD(I)   = TEMP(I) 
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
         VAR(I)%TOLD(IPHASE)  = PRO(I)%TEMP(IPHASE)
         VAR(I)%DTOLD(IPHASE) = PRO(I)%DTEMP(IPHASE)
         VAR(I)%XOLD(IPHASE)  = VAR(I)%X(IPHASE)
         IF(SOLUTION_TYPE == 'MULTI') THEN ! Multi-phase model
         VAR(I)%UOLD(IPHASE)  = VAR(I)%U(IPHASE)
         VAR(I)%VOLD(IPHASE)  = VAR(I)%V(IPHASE)
         VAR(I)%WOLD(IPHASE)  = VAR(I)%W(IPHASE)
         ENDIF
         ENDDO
         ENDIF
         UOLD(I)   = U(I) 
         VOLD(I)   = V(I) 
         WOLD(I)   = W(I) 
         POLD(I)   = PDIFF(I)
         ROLD(I)   = RO(I)
 2000 CONTINUE

      IF (NSCAL /= 0) THEN
         DO 2700 NS = 1,NSCAL
         DO 2700 L =1,NTOT1
         I         = L + IG1 - 1
         FIOLD(I,NS) = FI(I,NS)
 2700    CONTINUE
      ENDIF

      IF(TRANSL) THEN
         DO L = 1,NTOT1
            I         = L + IG1 - 1
            TRM(I)%GOLD   = TRM(I)%G
            TRM(I)%RETOLD = TRM(I)%RET
         ENDDO
      ENDIF

* 1001 IF(ITURB >= 3 .AND. ITURB /= 8) THEN

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         DO L = 1,NTOT1
            I         = L + IG1 - 1
            RKOLD(I)  = RK(I)
            EPSOLD(I) = REPS(I)
         ENDDO
      ENDIF

C ... END OF BLOCK LOOP 1

9000  CONTINUE

C ... NORMALIZE DROMAX FOR TIME ACCURATE CALCULATION

      IF(TIMEL) THEN
         IF(IDRXX == 1) THEN
            DROMAX = DROMAX/FRSDEN
         ELSEIF(IDRXX <= 4) THEN
            DROMAX = DROMAX/(REFVEL)
         ELSEIF(IDRXX == 5) THEN
            DROMAX = DROMAX/(FRSTEM)
c         ELSEIF(IDRXX == 6) THEN
c            DROMAX = DROMAX/(.5*FRSDEN*FRSVEL**2)
         ENDIF
      ENDIF
       

C **********************************************************************
C ... FOLLOW THE HISTORY OF FORCE COEFFICIENTS
C **********************************************************************

c      IF(ICYCTI /= 0 .OR. IOLD == 0) THEN
      IF(ICYCTI /= 0) THEN

      ISD     = 1

      DO 8200 N = 1,NBLOCK                 ! BLOCK LOOP 2 BEGINS

      NGL     = NPROCE(1+N,IPRO)           ! Global block number

      IF(BLKS(NGL)%CONVL) THEN             ! Blockwise test

      ISU     = ISD + NPATCH(N)

      DO 8199 N1 = ISD,ISU-1

         IF(GROUP) THEN         ! only uppercase alphabetics are monitored

         ICON1  = ICON((N1-1)*IC9+1) ! BC - type
         IF(ICON1 == 8 .OR. ICON1 == 9 .OR. ICON1 == 10) THEN
            DO J=1,26        ! Scan through uppercase alphabetics
               IF(BOUNDF(N1)(J:J) /= ' ') THEN
                  TOM     = TOM + TOMEGA(N1)
                  CX      = CX  + CXB(N1)
                  CY      = CY  + CYB(N1)
                  CZ      = CZ  + CZB(N1)
                  CMX     = CMX + CMXB(N1)
                  CMY     = CMY + CMYB(N1)
                  CMZ     = CMZ + CMZB(N1)
                  GOTO 8199
               ENDIF
            ENDDO ! alphabetics
         ENDIF

         ELSE                   ! all force patches are monitored

         TOM     = TOM + TOMEGA(N1)
         CX      = CX  + CXB(N1)
         CY      = CY  + CYB(N1)
         CZ      = CZ  + CZB(N1)
         CMX     = CMX + CMXB(N1)
         CMY     = CMY + CMYB(N1)
         CMZ     = CMZ + CMZB(N1)
         ENDIF ! (GROUP)
8199  CONTINUE
      ISD     = ISU
      ENDIF ! BLKS(NGL)%CONVL
8200  CONTINUE                             ! END OF BLOCK LOOP 2


      IF(.NOT.GROUP) THEN
         CXIO    = QMXIN+QMXOUT
         CYIO    = QMYIN+QMYOUT
         CZIO    = QMZIN+QMZOUT
         CMXIO   = CMXIN+CMXOUT
         CMYIO   = CMYIN+CMYOUT
         CMZIO   = CMZIN+CMZOUT
      ENDIF

      IF(ICYCTI /= 0) THEN
         CX  = (CX-CXIO)/(AREF*REFPRE)
         CY  = (CY-CYIO)/(AREF*REFPRE)
         CZ  = (CZ-CZIO)/(AREF*REFPRE)
         CMX = (CMX-CMXIO)/(AREF*CHLREF*REFPRE)
         CMY = (CMY-CMYIO)/(AREF*CHLREF*REFPRE)
         CMZ = (CMZ-CMZIO)/(AREF*CHLREF*REFPRE)
         CL  = (-CX*SIN(ALPHA) + CY*COS(ALPHA))
         CD  = CX*COS(ALPHA)*COS(BETA) + CY*SIN(ALPHA)*COS(BETA)
     +        -CZ*SIN(BETA)
         CS  = CX*COS(ALPHA)*SIN(BETA) + CY*SIN(ALPHA)*SIN(BETA)
     +        +CZ*COS(BETA)
      ENDIF

C ... is anyone looking these variables during iteration?? PPR 29.9.98

      QMEIN2 = QMEIN /(QMFIN +1.E-20)
      QMEOU2 = QMEOUT/(QMFOUT+1.E-20)
      EFFIC  = 100.*ABS(QMFIN*(QMEOU2-QMEIN2)/(TOM+1.E-20)) ! ADD HOC
      HEADT  = (QMEOU2 - QMEIN2)/G0


C ... DO PERFOMANCE CHECK FOR MPI

      IF(PARALLEL) THEN
         WTIMEN = MPI_WTIME()
         TTOTAL = WTIMEN - WTIMST   ! TOTAL         TIME IN ITERATION SWEEP
         COMTIM = WTBCEN - WTBCST   ! COMMUNICATION TIME IN ITERATION SWEEP
         WCAL   = WCAL + TTOTAL     ! TOTAL TIME
         WCOM   = WCOM + COMTIM     ! TOTAL COMMUNICATION TIME
         WTIMST = WTIMEN
         WRITE(45,1111) ICYCLE,TTOTAL,COMTIM,COMTIM/TTOTAL,
     +        COMTIM/(TTOTAL-COMTIM)
         IF(.NOT.CONVL) DROMAX = COMTIM/TTOTAL
      ENDIF
 1111 FORMAT(I6,1X,6F10.3)

C ... CONVERGENCE MONITOR (ON TERMINAL AND MONITOR FILES)

      IF(ICYCTI == 0) THEN
         DROMAX = 0.0
         CL     = 0.0
         CD     = 0.0
      ENDIF

      CALL CONMON(ICYCLE,DROMAX,CL,CD,CMZ,IXERR,INERR,IMACH,
     +     IMAX(1,INERR),JMAX(1,INERR),IPRO,NPROCE,MBPRO,
     +     NPRO,CONVL)

c      CALL WHMON(ICYCLE,DWMAX,SUMDWH,IPRO,NPROCE,MBPRO,
c     +     NPRO,CONVL)

C ... SEND THE FORCE COEFFICIENTS TO THE MONITORING PROCESSOR

      DO 8300 N = 1,NBLOCK                 ! BLOCK LOOP 3 BEGINS

      DO I = 1,30
         RESI(I,N) = 0.
      ENDDO

      NGL     = NPROCE(1+N,IPRO)           ! Global block number

      IF(BLKS(NGL)%CONVL) THEN             ! Blockwise test

      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSVEL  = BLKS(NGL)%FRSVEL
      FRSVIS  = BLKS(NGL)%FRSVIS
      FRSTEM  = BLKS(NGL)%FRSTEM
      TEMINI  = BLKS(NGL)%TEMINI
      FRSSIE  = BLKS(NGL)%FRSSIE
      SIEINI  = BLKS(NGL)%SIEINI
      CHLREF  = BLKS(NGL)%CHLREF
      RMUINI  = BLKS(NGL)%RMUINI
      RKLIM   = BLKS(NGL)%RKLIM
      EPSLIM  = BLKS(NGL)%EPSLIM
      FRSMUT  = BLKS(NGL)%FRSMUT
      SOLTEM  = BLKS(NGL)%SOLTEM

C ... L2-NORM OF THE RESIDUALS


      IG1     = IG(1,N)
      IR1     = IR(1,N)
      IQ1     = IQ(1,N)
      IGL     = 1
      IF (TIMEL) THEN
         IGL     = IG1
      ENDIF
      IP1     = 1
      IF(MULPHL) THEN
         IP1  = IG1
       ENDIF
      IT1 = 1
      IF(TRANSL) THEN
         IT1 = IG1
      ENDIF


      CALL SETV12(F1E(IG1),E(IG1),NTOT(1,N))  ! STORE TOTAL ENERGY
      CALL PRIMVE(U(IG1),V(IG1),W(IG1),RK(IG1),RO(IG1),F1E(IG1),
     + XC(IG1),YC(IG1),ZC(IG1),NTOT(1,N),FRSDEN,ITURB)
      CALL MULV12(F1E(IG1),RO(IG1),NTOT(1,N)) ! CREATE INT.ENERGY/VOLUME

C ... Viscous dissipation

      CALL TURFUN(F1R(IG1),F1RM(IG1),F1RN(IG1),F1RW(IG1),HAT3(IG1),
     1 F1E(IG1),U(IG1),V(IG1),W(IG1),RK(IR1),REPS(IR1),DDEPS(IR1),
     1 EPS2(IG1),VIST(IG1),VIS(IG1),HAT1(IG1),HAT2(IG1),
     2 A1(IG1),A2(IG1),A3(IG1),VOL(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     3 A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     4 D1(IG1),D2(IG1),D3(IG1),XC(IG1),YC(IG1),ZC(IG1),
     5 IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ITURB,1,ZZZ,VTRAN(IR1),
     6 MAXW,IDER(N),NGL)

      CALL RESL3(RO(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),FI(IQ1,1),
     +     RK(IR1),REPS(IR1),PTUR(IR1),VOL(IG1),P(IG1),PDIFF(IG1),
     +     VAR(IP1),F1E(IG1),
     +     TEMP(IG1),DRO(IG1),DM(IG1),DN(IG1),DW(IG1),DE(IG1),
     +     DFI(IQ1,1),DRK(IR1),DEPS(IR1),RESI(1,N),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),ITURB,NSCAL,MAXSB,NRESI,
     +     RKLIM,EPSLIM,HAT2(IG1),DDEPS(IG1),RMAV3(IGL),RNAV3(IGL),
     +     RWAV3(IGL),TIMEL,MULPHL,NGL,TRANSL,TRM(IT1))
      ENDIF ! BLKS(NGL)%CONVL

8300  CONTINUE                             ! END OF BLOCK LOOP 3

      VOLTOB = 0.
      NTO3 = 0
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) ! Global block number
         IF(BLKS(NGL)%CONVL) THEN   ! Blockwise test
         NTO3 = NTO3 + IMAX(1,N)*JMAX(1,N)*KMAX(1,N)
         VOLTOB  = VOLTOB + VOLN(N) ! TOTAL VOLUME
         ENDIF ! BLKS(NGL)%CONVL
      ENDDO

C ... Into the OUTPUT-file

      IF(ICYCLE == IPRINT .OR. (.NOT.TIMEL .AND.ICYCLE >= ICMAX))THEN
         WRITE(3,*)
         WRITE(3,*) ' Residuals in different blocks'
         DO L = 1,7           ! only values up to epsilon
            WRITE(3,*)
            DO N = 1,NBLOCK
               NGL     = NPROCE(1+N,IPRO) !Global block number
               IF(BLKS(NGL)%CONVL) THEN   ! Blockwise test
               WRITE(3,8301) L,NGL,SQRT(RESI(L,N))/REAL(NTO3,4)
               ENDIF ! BLKS(NGL)%CONVL
            ENDDO
         ENDDO
         L = 49
            WRITE(3,*)
            DO N = 1,NBLOCK
               NGL     = NPROCE(1+N,IPRO) !Global block number
               IF(BLKS(NGL)%CONVL) THEN   ! Blockwise test
               WRITE(3,8301) L,NGL,SQRT(RESI(L,N))/REAL(NTO3,4)
               ENDIF ! BLKS(NGL)%CONVL
            ENDDO
      ENDIF
 8301 FORMAT('Residual ',I2,' in block ',I3,' is ',G15.5)

      DO 8195 L = 1,NRESI
      DO 8195 N = 2,NBLOCK
         NGL    = NPROCE(1+N,IPRO) ! Global block number
         IF(BLKS(NGL)%CONVL) THEN  ! Blockwise test
         RESI(L,1) = RESI(L,1) + RESI(L,N)
         ENDIF ! BLKS(NGL)%CONVL
 8195 CONTINUE

C ... L2-NORMS
      DO 8196 L = 1,21
         RESI(L,1) = SQRT(RESI(L,1))/REAL(NTO3,4)
 8196 CONTINUE
C ... AVERAGED VALUES
      DO 8197 L = 22,23
         RESI(L,1) = RESI(L,1)/VOLTOB
 8197 CONTINUE
      CLDRO     = ABS(DROMAX)+1.E-15
C ... L2-NORMS OF THE EXPLICIT RESIDUALS
      DO 8198 L = 31,37
         RESI(L,1) = SQRT(RESI(L,1))/REAL(NTO3,4)
 8198 CONTINUE
         RESI(39,1)= SQRT(RESI(39,1))/REAL(NTO3,4)
      DO L = 41,43 ! Nobody knows what are 44, 45
         RESI(L,1) = SQRT(RESI(L,1))/REAL(NTO3,4)
      ENDDO

C ... L2-NORMS OF THE INTERMITTENCY VARIABLES
      DO 8299 L = 49,50
         RESI(L,1) = SQRT(RESI(L,1))/REAL(NTO3,4)
 8299 CONTINUE

C ... L2-NORMS OF THE EXPLICIT RESIDUALS (multi phase)
      DO L = 62,67
         RESI(L,1) = SQRT(RESI(L,1))/REAL(NTO3,4)
      ENDDO

************************************************************************
C ... UPDATE IBDIAG-FILE !

      IF (TIMEL) THEN
         ICY = ICYTOT
      ELSE
         ICY = ICYCTI
      ENDIF

      CALL IBDWRI(9,IBOUT,ICY,T,CLDRO,RESI,'IBDIAG',PRN,CD,CL,CS,
     +     CMX,CMY,CMZ,IMACH,NPROCE,MBPRO,NPRO,NTO3,IPRO,NRESI,
     +     VOLTOB,CONVL,QMFIN,QMFOUT,QMEIN,QMEOUT,HEADT,EFFIC,TOM,ZZZ,
     +     JET,QVFIN,QVFOUT,DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH,
     +     SINK,TRIMA,QGFIN,QGFOUT,QGEIN,QGEOUT)

      CALL AUXFOR(155,CXB,CYB,CZB,CMXB,CMYB,CMZB,DXB,DYB,DZB,QTB,QWB,
     +     QHB,TOMEGA,ALPHA,BETA,ICON,BOUNDF,NBCS,IPRO,3,T,
     +     ICY,APATCH,XMOM,YMOM,ZMOM,REFPRE,AREF,CHLREF)

      IF(IPRO == 1) THEN
      WRITE(90,'(I7,4E15.5)') ICYCTI,(XMASSB(JJ),JJ=1,4)
      ENDIF

************************************************************************
      ENDIF ! CONVERGENCE HISTORY AND FORCE COEFFICIENTS
************************************************************************
       
C ... Initialize convergence variables for a free surface

      WHMAX   = 0. !WHMAX   =-99999. 8.1.2010 jil
      WHMIN   = 0. !WHMIN   =99999. 8.1.2010 jil
      DWMAX   = 0.
      SUMDWH  = 0.
      FLODWH  = 0.

      RETURN
      END SUBROUTINE SUBIN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ######################################################################
      SUBROUTINE STATIM(M,N,INN,REFLEL)
C ######################################################################

C ... A CALLING ROUTINE FOR THE EQUATION OF STATE AND TIME-STEP SIZE

      USE CHARACTERS
      USE INTEGERS
      USE MAIN_ARRAYS
      USE NS3CO          
 
      IMPLICIT NONE

      INTEGER :: M,N,NGL,IG1,IGN,IR1,IQ1,IC1,IF1,KSTATE,
     +           INN,II,ig11,ICHARP,IPHASE,IP1,iii,IT1,IDERI,IPC1,
     +           INTI,INTJ,INTK

      LOGICAL REFLEL         
       
      NGL     = NPROCE(1+N,IPRO) !Global block number
      IPRESC  = BLKS(NGL)%IPRESC
      NPHASE  = BLKS(NGL)%NPHASE

      IG1     = IG(M,N)
      IR1     = IR(M,N)
      IQ1     = IQ(M,N)
      IC1     = IC(M,N)
      IF1     = JF(M,N)
      IP1     = 1
      IGN     = 1
      IT1     = 1
      IPC1    = 1
      IF(MULPHL) THEN
         IP1  = IG1
         IGN  = IG(M,N) + NTOT(M,N) - 1
      ENDIF
      IF(TRANSL) THEN
         IT1  = IG1
         IGN  = IG(M,N) + NTOT(M,N) - 1
      ENDIF
      IF(IPRESC == 1) IPC1 = IG1

      KSTATE  = JSTATE(NGL,1)  
      IDERI   = IDER(N)
      INTI    = INTERI(N)
      INTJ    = INTERJ(N)
      INTK    = INTERK(N)

      FRSPRE  = BLKS(NGL)%FRSPRE
      FRSDEN  = BLKS(NGL)%FRSDEN
      FRSTEM  = BLKS(NGL)%FRSTEM
      FRSSIE  = BLKS(NGL)%FRSSIE
      ARTSSP  = BLKS(NGL)%ARTSSP
      FRSSSP  = BLKS(NGL)%FRSSSP
      FRSVEL  = BLKS(NGL)%FRSVEL
      TLOLIM  = BLKS(NGL)%TLOLIM
      TUPLIM  = BLKS(NGL)%TUPLIM

      IF(BLKS(NGL)%SOLUTION_TYPE /= 'SOLID') THEN

      IF (ITURB >= 24) THEN
        CALL ADDUUK(RK(IG1),FI(IG1,1),IMAX(M,N),JMAX(M,N),
     +  KMAX(M,N),MAXSB)
      ENDIF

C ... Calculate pressure from pressure difference and gravity

      CALL PDTOP(P(IG1),PDIFF(IG1),XC(IG1),YC(IG1),ZC(IG1),FRSDEN,
     +     FRSPRE,NTOT(M,N))

      END IF ! BLKS(NGL)%SOLUTION_TYPE /= 'SOLID'

C ... Calculate density, specific internal energy, density derivatives,
C ... viscosity, specific heat and thermal conductivity as a f(p,T)
C *********************************************************************

      IF(BLKS(NGL)%SOLUTION_TYPE /= 'MULTI'  .AND.
     +   BLKS(NGL)%SOLUTION_TYPE /= 'CAVIT') THEN

C ... Calculate the equation of state for a single-phase flow
C *********************************************************************
      
      CALL STATE(VIS(IG1),P(IG1),PDIFF(IG1),RO(IG1),E(IG1),TEMP(IG1),
     +     CP(IG1),CH(IG1),DRDP(IG1),DRDH(IG1),GAMMA,PR,PRT,NTOT(M,N),
     +     ITURB,KSTATE,FRSPRE,RGAS,VISU0,EXPSU,TSU0,E0REF,T0REF,
     +    IMAX(M,N),JMAX(M,N),KMAX(M,N),INN,RMULTV,FRSDEN,FRSVIS,FRSTEM,
     +    FRSSIE,F1RM(IG1),IUPPT(NGL))

      IF(MULPHL .AND. BLKS(NGL)%SOLUTION_TYPE == 'SOLID') THEN
        DO IPHASE = 1,NPHASES
        PRO(IG1:IGN)%TEMP(IPHASE) = TEMP(IG1:IGN)
        ENDDO
      ENDIF

      IF(ICYCLE == IPRINT) THEN
      WRITE(3,*) 'STATIM output:'
      CALL PRINYS(3,DRDPC,DRDP(IG1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF

      IPHASE = 1
C ... Modify DRODP according to PSEUCO        
      IF(BLKS(NGL)%IPRESC == 2)
     +CALL PREUPM(RO(IG1),PRO(IP1),DRDP(IG1),DRDH(IG1),U(IG1),V(IG1),
     +     W(IG1),FRSPRE,FRSDEN,ARTSSP,FRSSSP,NTOT(M,N),NGL,KSTATE,
     +     IPRESC,PSEUCO,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,IPHASE,
     +     MULPHL,INN)
C ... CALCULATE SOUND SPEED AND TRANSFORM E TO TOTAL ENERGY
      CALL CONSVA(RO(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),C(IG1),
     +     U(IG1),V(IG1),W(IG1),RK(IR1),P(IG1),DRDP(IG1),DRDH(IG1),
     +     XC(IG1),YC(IG1),ZC(IG1),PRO(IP1),VAR(IP1),0.  ,FRSVEL,
     +     NTOT(M,N),ITURB,JRIMP,FI(IQ1,1),MAXSB,NSCAL,IPRESC,PSEUCO,
     +     BLKS(NGL)%SOLUTION_TYPE,1,NGL,RNUT(IG1),REPS(IG1))
C *********************************************************************
      END IF ! BLKS(NGL)%SOLUTION_TYPE /= 'MULTI' .OR. 'CAVIT'

      IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI'  .OR.
     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

C ... Calculate the equation of state for the multiphase variables
C *********************************************************************

C      vikapisto?
c      F1EPS(1:NTOT(M,N))= PLE2(IP1:IGN)+FRSPRE ! Old pressure without g

      DO IPHASE = 1,NPHASE

      CALL STATEM(PRO(IP1),P(IG1),PDIFF(IG1),GAMMA,PR,PRT,NTOT(M,N),
     +     ITURB,BLKS(NGL)%ISTATE(IPHASE),FRSPRE,RGAS,
     +     VISU0,EXPSU,TSU0,E0REF,T0REF,IMAX(M,N),JMAX(M,N),
     +     KMAX(M,N),INN,IPHASE,BLKS(NGL)%FRADEN(IPHASE),RMULTV)

      IF(ICYCLE == IPRINT) THEN
      WRITE(3,*) 'STATIM output2:'
      ICHARP = BLKS(NGL)%ICHAR(IPHASE)
      F1RM(1:NTOT(M,N)) = PRO(IG1:IGN)%DRODP(IPHASE)
      F1RN(1:NTOT(M,N)) = PRO(IG1:IGN)%RO(IPHASE)
      CALL PRINYS(3,CHAR_PH(3,ICHARP),F1RM(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_PH(1,ICHARP),F1RN(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,ROC,RO(IG1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDIF
       
      ENDDO ! IPHASE

      CALL MFDENS(PRO(IP1),VAR(IP1),PDIFF(IG1),P(IG1),RO(IG1),TEMP(IG1),
     +    DRDP(IG1),DRDH(IG1),VIS(IG1),CH(IG1),CP(IG1),BLKS(NGL)%FRSTEM,
     +    NPHASE,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,INN,ICYCLE,
     +    BLKS(NGL)%SOLUTION_TYPE,BLKS(NGL)%FRSPRE,BLKS(NGL)%DFRSPRE,
     +    BLKS(NGL)%DFRSTEM,BLKS(NGL)%EVAP_TYPE,KSTATE,
     +    BLKS(NGL)%DFRSDEN,XC(IG1),YC(IG1),ZC(IG1),NGL,IREPEA,NREPEA)

      DO IPHASE = 1,NPHASE
C     Modify DRODP according to PSEUCO
      IF(BLKS(NGL)%IPRESC == 2)
     +CALL PREUPM(RO(IG1),PRO(IP1),DRDP(IG1),DRDH(IG1),U(IG1),V(IG1),
     +     W(IG1),FRSPRE,FRSDEN,ARTSSP,FRSSSP,NTOT(M,N),NGL,KSTATE,
     +     IPRESC,PSEUCO,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,IPHASE,
     +     MULPHL,INN)
      ENDDO ! IPHASE

C ... CALCULATE SOUND SPEED AND TRANSFORM E TO TOTAL ENERGY
      CALL CONSVA(RO(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),C(IG1),
     +     U(IG1),V(IG1),W(IG1),RK(IR1),P(IG1),DRDP(IG1),DRDH(IG1),
     +     XC(IG1),YC(IG1),ZC(IG1),PRO(IP1),VAR(IP1),0.  ,FRSVEL,
     +     NTOT(M,N),ITURB,JRIMP,FI(IQ1,1),MAXSB,NSCAL,IPRESC,PSEUCO,
     +     BLKS(NGL)%SOLUTION_TYPE,NPHASE,NGL,RNUT(IG1),REPS(IG1))
      
C ... Calculate bubble diameter and corresponding variables
      IF(MULPHL) THEN
         CALL INTPRO(PRO(IP1),VAR(IP1),BLKS,IMAX(M,N),JMAX(M,N),
     +   KMAX(M,N),NGL,M)
      ENDIF
      
C *********************************************************************
      END IF ! BLKS(NGL)%SOLUTION_TYPE == 'MULTI' .OR. 'CAVIT'

C ... REFLECT THE VELOCITIES UNDER THE SURFACES

      IF(REFLEL) THEN 
      CALL REFLEC(U(IG1),V(IG1),W(IG1),XC(IG1),YC(IG1),ZC(IG1),
     1 IMAX(M,N),JMAX(M,N),KMAX(M,N),IBOT(M,N),ITOP(M,N),JBOT(M,N),
     2 JTOP(M,N),KBOT(M,N),KTOP(M,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     3 GAMMA,E0REF,T0REF,FRSDEN,FRSPRE,NGL,NPATCH(N),ICON(IC1),
     4 RO(IG1),P(IG1),PDIFF(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),
     5 VIS(IG1),VIST(IG1),CH(IG1),EPS2(IG1),TEMP(IG1),C(IG1),DRDP(IG1),
     6 DRDH(IG1),RK(IR1),REPS(IR1),DDEPS(IR1),PRO(IP1),VAR(IP1),
     7 RGAS,T0,ITURB,KSTATE,IB,MAXB,RCON(IC1),UWALL(1),VWALL(1),
     8 WWALL(1),TWALL(1),IHF(1,M),TLOLIM,TUPLIM,BLKS(NGL)%SOLUTION_TYPE,
     9 MULPHL,REFLECL,TRANSL,TRM(IT1),
     1 A1XA(IG1),A1YA(IG1),A1ZA(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),
     2 A2XA(IG1),A2YA(IG1),A2ZA(IG1),PRC(IPC1),IPRESC)
      ENDIF
       
c      WRITE(3,*) 'STATIM after REFLEC'
c      CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
c     + IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
       
c      IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI'  .OR.
c     +   BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

      IF(IPRESC == 1) THEN

      F1RM(IG1:IGN) = PRC(IG1:IGN)%DPDX ! Start from the current values
      F1RN(IG1:IGN) = PRC(IG1:IGN)%DPDY
      F1RW(IG1:IGN) = PRC(IG1:IGN)%DPDZ

c      CALL PDER(PDIFF(IG1),F1RM(IG1),F1RN(IG1),F1RW(IG1),A1(IG1),
c     +     A2(IG1),A3(IG1),VOL(IG1),D1(IG1),D2(IG1),D3(IG1),A1XA(IG1),
c     +     A1YA(IG1),A1ZA(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),
c     +     A3YA(IG1),A3ZA(IG1),XC(IG1),YC(IG1),ZC(IG1),IMAX(M,N),
c     +     JMAX(M,N),KMAX(M,N),IN,JN,KN,IDERI,M,
cc     +     INTI,INTJ,INTK,RKSI(IG1),INCHIML,NGL)
c     +     -3.,-3.,-3.,RKSI(IG1),INCHIML)

c      CALL TEMDEP(PDIFF(IG1),F1RM(IG1),F1RN(IG1),F1RW(IG1),A1(IG1),
c     +     A2(IG1),A3(IG1),VOL(IG1),D1(IG1),D2(IG1),D3(IG1),A1XA(IG1),
c     +     A1YA(IG1),A1ZA(IG1),A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),
c     +     A3YA(IG1),A3ZA(IG1),XC(IG1),YC(IG1),ZC(IG1),IMAX(M,N),
c     +     JMAX(M,N),KMAX(M,N),IN,JN,KN,IDERI,M)

c      PRC(IG1:IGN)%DPDX = F1RM(IG1:IGN)
c      PRC(IG1:IGN)%DPDY = F1RN(IG1:IGN)
c      PRC(IG1:IGN)%DPDZ = F1RW(IG1:IGN)

      END IF ! IPRESC == 1

      IF(ICYCLE == IPRINT) THEN

      WRITE(3,*) 'STATIM after REFLEC'
c      CALL PRINYS(3,UC,U(IG1),IT(1,N),IL(1,N),0,IK(1,N),
c     + IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      CALL PRINYS(3,UC,U(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     + IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)

      IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
      DO IPHASE = 1,BLKS(NGL)%NPHASE
      ICHARP = BLKS(NGL)%ICHAR(IPHASE)
      F1RM(1:NTOT(M,N)) = VAR(IP1:IGN)%U(IPHASE)
      CALL PRINYS(3,CHAR_VAR(33,ICHARP),F1RM(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      ENDDO
      ENDIF                 

      IF(MULPHL) THEN

      DO IPHASE = 1,BLKS(NGL)%NPHASE
      ICHARP = BLKS(NGL)%ICHAR(IPHASE)
      WRITE(3,*) 'STATIM output3:'
      F1RM(1:NTOT(M,N)) = VAR(IP1:IGN)%ALFA(IPHASE)
      F1RN(1:NTOT(M,N)) = VAR(IP1:IGN)%X(IPHASE)
      F1RW(1:NTOT(M,N)) = PRO(IP1:IGN)%DTEMP(IPHASE)
      F1E(1:NTOT(M,N))  = PRO(IP1:IGN)%DRODP(IPHASE)
      F1RK(1:NTOT(M,N)) = VAR(IP1:IGN)%ARO(IPHASE)
      F1EPS(1:NTOT(M,N))= VAR(IP1:IGN)%ARE(IPHASE)
      CALL PRINYS(3,CHAR_VAR(1,ICHARP),F1RM(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RN(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_PH(2,ICHARP),F1RW(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,ROC,RO(IG1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_PH(3,ICHARP),F1E(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_VAR(16,ICHARP),F1RK(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      CALL PRINYS(3,CHAR_VAR(17,ICHARP),F1EPS(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      F1EPS(1:NTOT(M,N))= PRO(IP1:IGN)%DPSDT
      CALL PRINYS(3,DPSDTC,F1EPS(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
      F1RK(1:NTOT(M,N)) = PRO(IP1:IGN)%E(IPHASE)
      CALL PRINYS(3,CHAR_PH(12,ICHARP),F1RK(1),
     + IT(M,N),IL(M,N),0,IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)

      ENDDO
      ENDIF ! MULPHL
      ENDIF

      RETURN
      END SUBROUTINE STATIM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ######################################################################
      SUBROUTINE BOUNDI(M)
C ######################################################################

C ... A CALLING ROUTINE FOR UPDATING THE BOUNDARY CONDITIONS.

      USE CHARACTERS

      USE INTEGERS,    ONLY : IPRO,MAXSB,MAXEB,MAXSS,MAXB,IB

      USE MAIN_ARRAYS, ONLY : U,V,W,XC,YC,ZC,IMAX,JMAX,KMAX,NPATCH,
     &   RO,P,PDIFF,RM,RN,RW,E,VIS,VIST,CH,EPS2,TEMP,C,DRDP,DRDH,RK,
     &   REPS,DDEPS,FI,BIJ,TIJ,TTS,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,
     &   A3YA,A3ZA,F1R,F1RM,F1RN,S11,JLOC,JTRA,APP,XCO,YCO,ZCO,JET,
     &   ZZZ,IG,IR,IQ,JF,IC,INTERI,INTERJ,INTERK,NPROCE,A1,A2,A3,
     &   VOL,D1,D2,D3,DRO,DM,DN,DW,DE,UROT,VROT,WROT,F1RW,F1E,OHMI,
     &   ICP,JCP,KCP,IJMASK,IDI1,IDI2,IDI3,IBOT,JBOT,KBOT,ITOP,JTOP,
     &   KTOP,IT,IL,IK,CXB,CYB,CZB,CMXB,CMYB,CMZB,LAMIN,ICON,F1RK,
     &   F1EPS,F1FI,NTOT,CP,JSTATE,PRO,VAR,BLKS,PTUR,
     &   WH,XXI,YXI,XETA,YETA,WFS,DTL,XHULL,YHULL,ZHULL,MHULL,MMHUL,
     &   IBOTGR,JBOTGR,KBOTGR,XVL,YVL,ZVL,ICONH,ICOGH,IHULL,LHULL,
     &   ITOPGR,JTOPGR,KTOPGR,WAVEH,IHF,IWAPP,TRM,QWALL,PRC

      USE NS3CO,       ONLY : IN,JN,KN,NBLOCK,ISTRES,FRSDEN,FRSPRE,
     &   ITURB,NSCAL,GX,GY,GZ,GROUND,QMFIN,QMFOUT,QMEIN,QMEOUT,QMXIN,
     &   QMXOUT,QMYIN,QMYOUT,QMZIN,QMZOUT,STRESL,GAMMA,E0REF,
     &   IPRINT,DT,T0REF,ALPHA,ICYCLE,PR,PRT,PRN,FRSVEL,RKLIM,EPSLIM,
     &   TURLIM,IUPDAT,IOLD,RGAS,VISU0,EXPSU,TSU0,JRIMP,IPRESC,PSEUCO,
     &   FRSTEM,FRSSSP,ISTATE,MULPHL,IFSBC,INWH,NFSD,ICFST,
     &   DTWMAX,DWMV,JFIRST,ICMAX,FRESUL,DWMAX,WHMAX,WHMIN,NPHASE,
     &   TRANSL,TWO_FLUIDL

      IMPLICIT NONE

      INTEGER :: M,IG1,IR1,IQ1,IF1,IC1,IE1,N,IG2,IR2,IQ2,IC2,
     &   IF2,IE2,IQQ2,IP1,IGN,INTI,INTJ,INTK,NGL,NS,KSTATE,III,IPHASE,
     &   IT1,ICHARP

      REAL, ALLOCATABLE :: DPU1(:)

      IF(MULPHL) THEN
         ALLOCATE(DPU1(MAXB))
      ENDIF

C***********************************************************************
C ... UPDATE INTERNAL BOUNDARIES
C***********************************************************************

      DO N = 1,NBLOCK

         IG1 = IG(1,N)
         IR1 = IR(1,N)
         IQ1 = IQ(1,N)
         IF1 = JF(1,N)
         IC1 = IC(1,N)
         IE1 = 1
         IP1 = 1
         IT1 = 1

         IF(ISTRES > 0) IE1 = IG1
         IF(MULPHL)     IP1 = IG1
         IF(TRANSL)     IT1 = IG1

         IGN = IG1 + NTOT(1,N) - 1
         NGL = NPROCE(1+N,IPRO)!Global block number

      
C ... REFLECTION OF VARIABLES INTO THE EMPTY CORNERS

         CALL REFCOR(U(IG1),V(IG1),W(IG1),XC(IG1),YC(IG1),ZC(IG1),
     &               IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     &               FRSDEN,FRSPRE,N,NPATCH(N),RO(IG1),P(IG1),
     &               PDIFF(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),
     &               VIS(IG1),VIST(IG1),CH(IG1),EPS2(IG1),TEMP(IG1),
     &               C(IG1),DRDP(IG1),DRDH(IG1),RK(IR1),REPS(IR1),
     &               DDEPS(IR1),FI(IQ1,1),BIJ(IE1,1),PRO(IP1),VAR(IP1),
     &               TRANSL,TRM(IT1),MULPHL,ITURB,NSCAL,MAXSB,MAXEB)

      ENDDO

C***********************************************************************
C
C ... FIRSTLY CONNECT ALL THE PATCHES TOGETHER (BLOCK LOOP INSIDE)
C
C***********************************************************************

      CALL MIR(TEMP,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,PRO,VAR,
     +     TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +     1,NBLOCK,1.,0,PRC)
                    
      IF(ICYCLE == IPRINT) THEN

      IF(MULPHL) THEN
         IPHASE = 1
         F1RM(1:MAXB) = VAR(1:MAXB)%V(IPHASE)
         F1RN(1:MAXB) = VAR(1:MAXB)%W(IPHASE)
      ENDIF
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         IF(IL(1,N) >= 1 .AND. IL(1,N) <= 3) THEN
         WRITE(3,*)
         WRITE(3,*) ' Before first CYC'
         IF(TWO_FLUIDL) THEN
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),F1RM(IG1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSE
         CALL PRINYS(3,VC, V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         ENDIF
      ENDDO

      ENDIF ! ICYCLE == IPRINT

                    
      CALL CYC(U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +1,NBLOCK,F1R,F1RM,VAR)

      CALL MIR(TEMP,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,PRO,VAR,
     +     TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +     1,NBLOCK,1.,0,PRC) ! This is done twice

       IF(ICYCLE == IPRINT) THEN

      IF(MULPHL) THEN
         IPHASE = 1
         F1RM(1:MAXB) = VAR(1:MAXB)%V(IPHASE)
         F1RN(1:MAXB) = VAR(1:MAXB)%W(IPHASE)
      ENDIF
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         IF(IL(1,N) >= 1 .AND. IL(1,N) <= 3) THEN
         WRITE(3,*)
         WRITE(3,*) ' Before second CYC'
         IF(TWO_FLUIDL) THEN
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),F1RM(IG1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSE
         CALL PRINYS(3,VC, V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         ENDIF
      ENDDO

      ENDIF ! ICYCLE == IPRINT

      CALL CYC(U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +1,NBLOCK,F1R,F1RM,VAR)  ! Also this, because CYC should dominate

      IF(IPRESC == 1) THEN
         F1RN(1:MAXB) = PRC(1:MAXB)%DPDX
         F1RW(1:MAXB) = PRC(1:MAXB)%DPDY
         F1E(1:MAXB)  = PRC(1:MAXB)%DPDZ
      CALL CYCV(F1RN,F1RW,F1E,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,
     + A3ZA,M,1,NBLOCK,F1R,F1RM,VAR)
         PRC(1:MAXB)%DPDX = F1RN(1:MAXB)
         PRC(1:MAXB)%DPDY = F1RW(1:MAXB)
         PRC(1:MAXB)%DPDZ = F1E(1:MAXB)
      ENDIF

      IF(ICYCLE == IPRINT) THEN

      IF(TWO_FLUIDL) THEN
         IPHASE = 1
         F1RM(1:MAXB) = VAR(1:MAXB)%V(IPHASE)
      ENDIF
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         IF(IL(1,N) >= 1 .AND. IL(1,N) <= 3) THEN ! Check
         WRITE(3,*)
         WRITE(3,*) ' After second CYC'
         IF(TWO_FLUIDL) THEN
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),F1RM(IG1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSE
         CALL PRINYS(3,VC, V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         ENDIF
      ENDDO

      ENDIF ! ICYCLE == IPRINT

      IF(FRESUL) THEN

      IF(IFSBC == 1) THEN ! Original mod (free-stream only) 

         CALL FOU(RO,RM,RN,RW,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,MAXSB,
     +     A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     +     1,NBLOCK,1.,VOL,XC,YC,ZC,WH)

         IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' FOU from BOUNDI'
         DO N = 1,NBLOCK
           IG1 = IG(1,N)
           CALL PRINYS(3,pdifC,pdiff(IG1),IT(1,N),
     +     IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDDO
         ENDIF

      ELSE IF(IFSBC == 2) THEN

         STOP 'In BOUNDI: Level-set does not work man! Change IFSBC'

      ELSE IF(IFSBC == 3) THEN ! Probably does not worka

         CALL FRE(M,1,NBLOCK,1.,1)  ! ESa 12.11.2000
*                 M,1,NBLOCK,1.,0)  ! Only BC:s 1.10.06

      ELSE IF(IFSBC == 4) THEN ! Recommended free-surface conditions
      
         CALL FRE(M,1,NBLOCK,1.,0)

      ENDIF ! IFSBC == 1

      ENDIF ! FRESUL

C ... SCALAR EQUATIONS IF Reynolds stress model is applied

      IF (ITURB >= 21) THEN
         CALL MIRR(FI(1,1),FI(1,2),FI(1,3),FI(1,4),FI(1,5),FI(1,6),
     +        A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +        M,1,NBLOCK,MAXSB)
      ENDIF

C ... Old ASM versions

      IF (ITURB >= 10 .AND. ITURB <= 19) THEN

         CALL MIRR(S11(1,1),S11(1,2),S11(1,3),S11(1,4),
     +        S11(1,5),S11(1,6),A1XA,A1YA,A1ZA,A2XA,
     +        A2YA,A2ZA,A3XA,A3YA,A3ZA,M,1,NBLOCK,MAXSS)
      ENDIF

C ... New EARSM versions

      IF (ISTRES >= 1) THEN

         CALL MIRR(BIJ(1,1),BIJ(1,2),BIJ(1,3),BIJ(1,4),BIJ(1,5),
     +        BIJ(1,6),A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +        M,1,NBLOCK,MAXEB)

C ... Currently S11 needs a VORFUN call to be updated

         CALL MIRR(S11(1,1),S11(1,2),S11(1,3),S11(1,4),
     +        S11(1,5),S11(1,6),A1XA,A1YA,A1ZA,A2XA,
     +        A2YA,A2ZA,A3XA,A3YA,A3ZA,M,1,NBLOCK,MAXSS)
      ENDIF

C ... Calculate the material properties before connections

         IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' STATIM from BOUNDI'
         ENDIF

      DO N = 1,NBLOCK
         CALL STATIM(1,N,2,.TRUE.)
      ENDDO

C ... Connection after mirror and cycle

      CALL CONNEC(TEMP,U,V,W,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,FI,S11,BIJ,
     + PRO,VAR,TRM,BLKS,MAXEB,MAXSB,MAXSS,F1R,F1RM,JLOC,JTRA,APP,M,1,
     + NBLOCK,BLKS(NGL)%SOLUTION_TYPE)

      DO N = 1,NBLOCK
         CALL STATIM(1,N,2,.TRUE.)
      ENDDO

C ... Pressure gradient after STATIM

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*) ' AFTER CONNEC in BOUNDI'
         DO N = 1,NBLOCK
           IG1 = IG(1,N)
           CALL PRINYS(3,ROC,RO(IG1),IT(1,N),
     +     IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
           CALL PRINYS(3,EC,E(IG1),IT(1,N),
     +     IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
        ENDDO
      ENDIF

C ... Connection of a solid-type patches

      CALL CONEPA_SP(QWALL,'SOL')

      CALL CYC(U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +1,NBLOCK,F1R,F1RM,VAR)! Also this, because CYC should dominate

C ... Update double-precision variables

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
         DPU1(1:MAXB) = VAR(1:MAXB)%DTEMP(IPHASE)
         CALL CONECX(DPU1,1,1,NBLOCK,0)
         VAR(1:MAXB)%DTEMP(IPHASE) = DPU1(1:MAXB)
         ENDDO
      ENDIF
       
C ... Update sliding surfaces

      CALL SLD(U,V,W,EPS2,VIST,RO,TEMP,P,PDIFF,PTUR,E,RM,RN,RW,RK,REPS,
     +     FI,S11,BIJ,TIJ,PRO,VAR,TRM,F1R,F1RM,F1RN,MAXB,MAXEB,MAXSB,
     +     MAXSS,XCO,YCO,ZCO,M,1,NBLOCK,ISTRES,D1,D2,D3,VOL,DM,DN,DW,1)

C ... Update periodic surfaces

      CALL PERNEC(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,PRO,VAR,TRM,XC,YC,ZC,
     +     MAXB,1,1,NBLOCK,VIST,EPS2)
 
      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,'(A)') ' Boundaries have been updated'

      DO N = 1,NBLOCK
         NGL    = NPROCE(1+N,IPRO)
         IG1    = IG(1,N)
         CALL PDTOP(P(IG1),PDIFF(IG1),XC(IG1),YC(IG1),ZC(IG1),FRSDEN,
     +     FRSPRE,NTOT(1,N))

         CALL PRINYS(3,PDIFC,PDIFF(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +        IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,UC,U(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDDO                    

      IF(MULPHL) THEN
         IPHASE = 1
         F1RM(1:MAXB) = VAR(1:MAXB)%DM(IPHASE)
c         F1RM(1:MAXB) = VAR(1:MAXB)%V(IPHASE) ! 34
         F1RN(1:MAXB) = VAR(1:MAXB)%W(IPHASE)
      ENDIF
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         IF(IL(1,N) >= 1 .AND. IL(1,N) <= 3) THEN
         WRITE(3,*)
         WRITE(3,*) ' After SLD, MIR and CYC'
         IF(TWO_FLUIDL) THEN
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),F1RM(IG1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSE
         CALL PRINYS(3,VC, V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
      ENDIF
      ENDDO

      ENDIF ! ICYCLE == IPRINT


      CALL MIR(TEMP,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,PRO,VAR,
     +     TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,M,
     +     1,NBLOCK,1.,0,PRC)

      RETURN
      END SUBROUTINE BOUNDI
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
      SUBROUTINE GHOSTEXT(FRSVEL,ALPHA,BETA,IMAX,JMAX,KMAX,IN,JN,KN,
     +     XC,YC,ZC,NPATCH,ICONC,
     +     FRSDEN,FRSPRE,FRSSIE,FRSTEM,FRSVIS,RKLIM,
     +     EPSLIM,FRSMUT,SOLTEM,ITURB,NGL,RO,P,PDIFF,RM,RN,RW,U,V,W,E,
     +     TEMP,VIS,RK,REPS,VIST,EPS2,RNUT,BLKS,PRO,
     +     VAR,MULPHL,INITC,FRSG,FRSRET,TRANSL,TRM,CENAX,CENAY,CENAZ)

C ... This subroutine is created by copying and editing subroutine DONORS

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9, GROUND, GX, GY,GZ, G0, ALTITUDE,
     &     GAMMA, RGAS, E0REF, T0REF, VISU0, EXPSU, TSU0,
     &     XMOM, YMOM, ZMOM
      USE CONSTANTS, ONLY : PII

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,L,NPATCH,IFACE
      INTEGER :: IA,IM,JA,JM,IFIRST,ILAST,JFIRST,JLAST,KFIRST,KLAST
      INTEGER :: K,KA,J,JJ,I,KR,ITURB,III,NGL,IPHASE,INITC

      INTEGER, DIMENSION(IC9,*) :: ICONC

      REAL    :: FRSVEL,ALPHA,BETA,U10,V10,W10
      REAL    :: FRSDEN,FRSPRE,FRSSIE,FRSTEM,FRSVIS,FRSG,FRSRET
      REAL    :: FRSDEN_S,FRSPRE_S,FRSSIE_S,FRSTEM_S,FRSVIS_S
      REAL    :: RKLIM,EPSLIM,FRSMUT,SOLTEM,RKHI
      REAL    :: GTOT,H,GM1

      REAL    :: RADCG, RADGP, CENAX, CENAY, CENAZ

      REAL    :: RO(*),P(*),PDIFF(*),RM(*),RN(*),RW(*),U(*),V(*),W(*),
     &           E(*),TEMP(*),VIS(*),RK(*),REPS(*),VIST(*),EPS2(*),
     &           RNUT(*)

      REAL :: XC(*), YC(*), ZC(*)

      LOGICAL :: MULPHL,TRANSL

      REAL, EXTERNAL :: ATMOSRHO, ATMOSP, ATMOST

      TYPE(BLOCKS)           :: BLKS(*)
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      GM1 = GAMMA - 1.0

      FRSDEN_S = FRSDEN
      FRSPRE_S = FRSPRE
      FRSTEM_S = FRSTEM
      FRSSIE_S = FRSSIE
      FRSVIS_S = FRSVIS

      U10 = FRSVEL*COS(ALPHA)*COS(BETA)  ! New velocities
      V10 = FRSVEL*SIN(ALPHA)
      W10 = FRSVEL*COS(ALPHA)*SIN(BETA)

      IF(INITC == 2) THEN  ! Pull-up simulation
         U10 = 0.0
         V10 = 0.0
         W10 = 0.0
      ENDIF

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN

C ... Update the free-stream boundary conditions (EXT)

      DO L = 1,NPATCH

         IF(ICONC(1,L) == 2) THEN

            IFACE = ICONC(3,L)
            IA    = ICONC(4,L)
            IM    = ICONC(5,L)
            JA    = ICONC(6,L)
            JM    = ICONC(7,L)

            IF(IFACE == 1) THEN
               IFIRST = -1
               ILAST  =  0
               JFIRST = IA
               JLAST  = IM
               KFIRST = JA
               KLAST  = JM
            ELSEIF(IFACE == 2) THEN
               IFIRST = IA
               ILAST  = IM
               JFIRST = -1
               JLAST  =  0
               KFIRST = JA
               KLAST  = JM
            ELSEIF(IFACE == 3) THEN
               IFIRST = IA
               ILAST  = IM
               JFIRST = JA
               JLAST  = JM
               KFIRST = -1
               KLAST  =  0
            ELSEIF(IFACE == 4) THEN
               IFIRST = IMAX + 1
               ILAST  = IMAX + 2
               JFIRST = IA
               JLAST  = IM
               KFIRST = JA
               KLAST  = JM
            ELSEIF(IFACE == 5) THEN
               IFIRST = IA
               ILAST  = IM
               JFIRST = JMAX + 1
               JLAST  = JMAX + 2
               KFIRST = JA
               KLAST  = JM
            ELSEIF(IFACE == 6) THEN
               IFIRST = IA
               ILAST  = IM
               JFIRST = JA
               JLAST  = JM
               KFIRST = KMAX + 1
               KLAST  = KMAX + 2
            ENDIF


            DO K = KFIRST,KLAST
               KA = (KN+K-1)*ISTRID*JSTRID
               DO J = JFIRST,JLAST
                  JJ = (JN+J-1)*ISTRID + IN + KA 
                  DO I = IFIRST,ILAST
                     KR = JJ + I

                     IF(ALTITUDE >= GROUND) THEN  ! Flight altitude is given

                        GTOT   = SQRT(GX**2 + GY**2 + GZ**2)

                        IF(INITC == 2) THEN ! Symmetric pull-up

                           RADCG = SQRT((CENAX-XMOM)**2
     &                           + (CENAY-YMOM)**2)
                           RADGP = SQRT((CENAX-XC(KR))**2
     &                           + (CENAY-YC(KR))**2)

                           H  =  -(GX*XMOM + GY*YMOM + GZ*ZMOM 
     &                           + GTOT*GROUND)/GTOT + (RADCG-RADGP)

                        ELSE

*                           H = -(GX*XC(KR) + GY*YC(KR) + GZ*ZC(KR)
*     &                         + GTOT*GROUND)/GTOT

                           H = -(GX*(XC(KR)-XMOM) +
     &                           GY*(YC(KR)-YMOM) +
     &                           GZ*(ZC(KR)-ZMOM) + 
     &                           GTOT*GROUND)/GTOT

                        ENDIF

                        FRSDEN = ATMOSRHO(H)                           
                        FRSPRE = ATMOSP(H)
                        FRSTEM = ATMOST(H)
                        FRSSIE = -E0REF + RGAS/GM1*(FRSTEM-T0REF)
                        FRSVIS = VISU0*FRSTEM**EXPSU/(FRSTEM+TSU0)

                     ENDIF

                     RO(KR)    = FRSDEN      
                     P(KR)     = FRSPRE
*                     PDIFF(KR) = 0.
                     PDIFF(KR) = FRSPRE - FRSPRE_S
                     RM(KR)    = FRSDEN*U10
                     RN(KR)    = FRSDEN*V10
                     RW(KR)    = FRSDEN*W10
                     U(KR)     = U10
                     V(KR)     = V10
                     W(KR)     = W10
                     E(KR)     = FRSDEN*(FRSSIE
     &                         + .5*(U10**2+V10**2+W10**2))
                     TEMP(KR)  = FRSTEM
                     VIS(KR)   = FRSVIS

C ... Free-stream turbulence values are specified here

                     IF(ITURB >= 3 .AND. ITURB /= 8) THEN

                        IF(ITURB /= 9) THEN   ! Two-equation models

                           RK(KR)   = RKLIM   ! FRSRK
                           REPS(KR) = EPSLIM  ! FRSEPS
                           VIST(KR) = FRSMUT  ! RMUINI
                           EPS2(KR) = 1. + FRSMUT/FRSVIS  ! RMUINI/FRSVIS
                           E(KR)    = E(KR) + RK(KR)

                        ELSE IF(ITURB == 9) THEN  ! Spalart-Allmaras

                           RKHI = 1.

                           DO III = 1,20
                              RKHI = SQRT(SQRT(FRSMUT/FRSVIS
     &                             *(RKHI**3+7.1**3)))
                           ENDDO  ! Iteration also in main.f (mersu)

                           REPS(KR) = RKHI*FRSVIS ! Value for RNUT*RO
                           RNUT(KR) = RKHI*FRSVIS/RO(KR)
                           VIST(KR) = FRSMUT ! RMUINI
                           EPS2(KR) = 1. + FRSMUT/FRSVIS ! RMUINI/FRSVIS
                           RK(KR)   = 0.

                        ENDIF   ! ITURB >= 3

                     ENDIF      ! Two-equation models

                     IF(BLKS(NGL)%SOLUTION_TYPE == 'SOLID') THEN
                        TEMP(KR) = SOLTEM
                     ENDIF

                     IF(MULPHL) THEN

                     DO IPHASE = 1,BLKS(NGL)%NPHASE
                        PRO(KR)%TEMP(IPHASE) = FRSTEM
                        PRO(KR)%DTEMP(IPHASE)= BLKS(NGL)%DFRSTEM
                        VAR(KR)%ALFA(IPHASE) = BLKS(NGL)%FRSALFA(IPHASE)
                        VAR(KR)%X(IPHASE)    = BLKS(NGL)%FRSX(IPHASE)
                        PRO(KR)%TSAT         = FRSTEM  ! Temporary init
                     ENDDO 

                     ENDIF  ! MULPHL

                     IF(TRANSL) THEN  ! Intermittency variables
                        TRM(KR)%G   = FRSG
                        TRM(KR)%RET = FRSRET
                     ENDIF

                  ENDDO  ! I
               ENDDO     ! J
            ENDDO        ! K

         ENDIF  ! EXT

      ENDDO

      FRSDEN = FRSDEN_S
      FRSPRE = FRSPRE_S
      FRSTEM = FRSTEM_S
      FRSSIE = FRSSIE_S
      FRSVIS = FRSVIS_S

      RETURN
      END SUBROUTINE GHOSTEXT


C  ATMOS FORTRAN ---------------------------------------------------------
C
C     International Standard Atmosphere
C
C  -- Kari Renko -- 28.12.1989 -------------------------------------------
C
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
      REAL FUNCTION ATMOSA(H)

C     Returns the speed of sound at a given height (ISO atmosphere)

      REAL :: H
      REAL :: GAMMA, R0

      GAMMA  = 1.4
      R0     = 287.05287

      ATMOSA = SQRT(GAMMA*R0*ATMOST(H))

      RETURN
      END FUNCTION ATMOSA
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
      REAL FUNCTION ATMOSP(H)

C     Returns static pressure at a given height (ISO atmosphere)
C     (H <= 80 km)

      USE NS3CO, ONLY : G0, RGAS

      REAL :: H
      REAL :: ALPHA, TB, HB, RHO0, R0, A0

      IF (H > 80.E3) STOP

      R0 = RGAS

      IF (H <= 11.0E3) THEN
        PB    = 1.01325E5
        ALPHA = -6.5E-3
        TB    = 288.15
        HB    = 0.
      ENDIF
      IF (H > 11.0E3 .AND. H <= 20.0E3) THEN
        PB    = 2.26320E4
        ALPHA = 0.0
        TB    = 216.65
        HB    = 11.0E3
      ENDIF
      IF (H > 20.0E3 .AND. H <= 32.0E3) THEN
        PB    = 5.47488E3
        ALPHA = 1.0E-3
        TB    = 216.65
        HB    = 20.0E3
      ENDIF
      IF (H > 32.0E3 .AND. H < 47.0E3) THEN
        PB    = 8.68016E2
        ALPHA = 2.8E-3
        TB    = 228.65
        HB    = 32.0E3
      ENDIF
      IF (H > 47.0E3 .AND. H <= 51.0E3) THEN
        PB    = 1.10906E2
        ALPHA = 0.0
        TB    = 270.65
        HB    = 47.0E3
      ENDIF
      IF (H > 51.0E3 .AND. H <= 71.0E3) THEN
        PB    = 6.69385E1
        ALPHA = -2.80E-3
        TB    = 270.65
        HB    = 51.0E3
      ENDIF
      IF (H > 71.0E3 .AND. H <= 80.0E3) THEN
        PB    = 3.95639
        ALPHA = -2.0E-3
        TB    = 214.65
        HB    = 71.0E3
      ENDIF

      IF (ALPHA /= 0.0) THEN
        ATMOSP = PB * (1 + (ALPHA/TB)*(H-HB))**(-(G0/(ALPHA*R0)))
      ELSE
        ATMOSP = PB * EXP(-(G0/(R0*TB))*(H-HB))
      ENDIF

      RETURN
      END FUNCTION ATMOSP
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
      REAL FUNCTION ATMOSRHO(H)

C     Returns density a given height (ISO atmosphere)
C     (H <= 80 km)

      USE NS3CO, ONLY : G0, RGAS

      REAL :: H
      REAL :: ALPHA, TB, HB, RHO0, R0, A0

      IF (H > 80.E3) STOP

      R0 = RGAS

      IF (H <= 11.0E3) THEN
        RHOB  = 1.225
        ALPHA = -6.5E-3
        TB    = 288.15
        HB    = 0.
      ENDIF
      IF (H > 11.0E3 .AND. H <= 20.0E3) THEN
        RHOB  = 0.363918
        ALPHA = 0.0
        TB    = 216.65
        HB    = 11.0E3
      ENDIF
      IF (H > 20.0E3 .AND. H <= 32.0E3) THEN
        RHOB  = 8.80347E-2
        ALPHA = 1.0E-3
        TB    = 216.65
        HB    = 20.0E3
      ENDIF
      IF (H > 32.0E3 .AND. H <= 47.0E3) THEN
        RHOB  = 1.32250E-2
        ALPHA = 2.8E-3
        TB    = 228.65
        HB    = 32.0E3
      ENDIF
      IF (H > 47.0E3 .AND. H <= 51.0E3) THEN
        RHOB  = 1.42753E-3
        ALPHA = 0.0
        TB    = 270.65
        HB    = 47.0E3
      ENDIF
      IF (H > 51.0E3 .AND. H <= 71.0E3) THEN
        RHOB  = 8.61601E-4
        ALPHA = -2.80E-3
        TB    = 270.65
        HB    = 51.0E3
      ENDIF
      IF (H > 71.0E3 .AND. H <= 80.0E3) THEN
        RHOB  = 6.42106E-5
        ALPHA = -2.0E-3
        TB    = 214.65
        HB    = 71.0E3
      ENDIF

      IF (ALPHA /= 0.0) THEN
        ATMOSRHO = RHOB * (1 + (ALPHA/TB)*(H-HB))**(-(1+G0/(ALPHA*R0)))
      ELSE
        ATMOSRHO = RHOB * EXP(-(G0/(R0*TB))*(H-HB))
      ENDIF

      RETURN
      END FUNCTION ATMOSRHO
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
      REAL FUNCTION ATMOST(H)

C     Returns static temperature at a given height (ISO atmosphere)
C     (H <= 80 km)

      REAL :: H
      REAL :: ALPHA, TB, HB

      IF (H > 80.E3) STOP

      IF (H <= 11.0E3) THEN
        ALPHA = -6.5E-3
        TB    = 288.15
        HB    = 0.
      ENDIF
      IF (H > 11.0E3 .AND. H <= 20.0E3) THEN
        ALPHA = 0.0
        TB    = 216.65
        HB    = 11.0E3
      ENDIF
      IF (H > 20.0E3 .AND. H <= 32.0E3) THEN
        ALPHA = 1.0E-3
        TB    = 216.65
        HB    = 20.0E3
      ENDIF
      IF (H > 32.0E3 .AND. H <= 47.0E3) THEN
        ALPHA = 2.8E-3
        TB    = 228.65
        HB    = 32.0E3
      ENDIF
      IF (H > 47.0E3 .AND. H <= 51.0E3) THEN
        ALPHA = 0.0
        TB    = 270.65
        HB    = 47.0E3
      ENDIF
      IF (H > 51.0E3 .AND. H <= 71.0E3) THEN
        ALPHA = -2.80E-3
        TB    = 270.65
        HB    = 51.0E3
      ENDIF
      IF (H > 71.0E3 .AND. H <= 80.0E3) THEN
        ALPHA = -2.0E-3
        TB    = 214.65
        HB    = 71.0E3
      ENDIF

      ATMOST  = TB + ALPHA*(H-HB)

      RETURN
      END FUNCTION ATMOST
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C      
