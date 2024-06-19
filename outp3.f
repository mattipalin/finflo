C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OPENER(CN)
C     *****************

      CHARACTER(LEN=3)  :: CN
      CHARACTER(LEN=50) :: INPUTC
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      INTEGER :: VALUES(8)
      
      IF(CN == '   ' .OR. CN == '001') THEN  ! This is master process
         INPUTC(1:50) = ' '
C     Pick the file name from command-line for INPUT, if one is available
         I0 = IARGC()
         IF (I0 > 0) THEN
            DO I=1,I0
               CALL GETARG(I,INPUTC)
               IF(INPUTC(1:1) /= '-') GOTO 100
            ENDDO
         ENDIF
         INPUTC = 'INPUT'
 100     CONTINUE
         
         OPEN(  2,FILE=  INPUTC,STATUS='UNKNOWN')
c        OPEN( 16,FILE='COMPUT',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN( 17,FILE= 'LOKKI',STATUS='UNKNOWN',FORM='FORMATTED')
         OPEN( 43,FILE='SCALAR',STATUS='UNKNOWN')
         OPEN(  9,FILE='IBDIAG',STATUS='UNKNOWN',FORM='UNFORMATTED')
c         OPEN( 19,FILE='IBDIAG.BAK',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN( 59,FILE='TBDIAG',STATUS='UNKNOWN',FORM='UNFORMATTED')
c         OPEN( 69,FILE='TBDIAG.BAK',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN( 90,FILE= 'MASSB',STATUS='UNKNOWN',FORM='FORMATTED')
      ENDIF ! Master process only
         
      OPEN(13,FILE='RUN.LOG'//CN,STATUS='UNKNOWN',FORM='FORMATTED')
c      OPEN(13,FILE='RUN.LOG'//CN,STATUS='UNKNOWN',FORM='FORMATTED',
c     &                           BUFFERED="NO")
      WRITE(13,*)

	
      OPEN(3, FILE= 'OUTPUT'//CN,STATUS='UNKNOWN')
      OPEN(4, FILE= 'FORCES'//CN,STATUS='UNKNOWN')
      OPEN(7, FILE='MONITOR'//CN,STATUS='UNKNOWN')
C     OPEN(21,FILE=    'KOE'//CN,STATUS='UNKNOWN')
C     OPEN(39,FILE= 'Cmdata'//CN,STATUS='UNKNOWN',FORM=  'FORMATTED')
      OPEN(44,FILE= 'MCHECK'//CN,STATUS='UNKNOWN')
      OPEN(45,FILE= 'MEMORY'//CN,STATUS='UNKNOWN')
      OPEN(46,FILE= 'GRIPUT'//CN,STATUS='UNKNOWN')
c      OPEN(999,FILE='GRID.BIN'//CN,FORM='UNFORMATTED') !TKM 27.4.2001 99->999
      IF(CN == '001') THEN
         OPEN(11, FILE='MONITOR',STATUS='UNKNOWN')
      ENDIF
          
      RETURN
      END SUBROUTINE OPENER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE HEADING(ICHAN)

C ... Print out the heading and version number

      WRITE(ICHAN,*)
      WRITE(ICHAN,*)
      WRITE(ICHAN,*)
     &'          FFFFFF   III     N   N     FFFFFF   L        OOOO     '
      WRITE(ICHAN,*)
     &'          F         I      NN  N     F        L       O    O    '
      WRITE(ICHAN,*)
     &'          FFFF      I      N N N     FFFF     L       O    O    '
      WRITE(ICHAN,*)
     &'          F         I      N  NN     F        L       O    O    '
      WRITE(ICHAN,*)
     &'          F        III     N   N     F        LLLLL    OOOO     '
      WRITE(ICHAN,*)
     &'                                                                '
      WRITE(ICHAN,*)
      WRITE(ICHAN,*)
     &'                          Release r001 -> r822                  '
      WRITE(ICHAN,*)
      WRITE(ICHAN,*)
     &'                  Copyright © 2001 - 2022 Elomatic Ltd.         '
      WRITE(ICHAN,*)
      WRITE(ICHAN,*)

      WRITE(ICHAN,*)'**************************************************'
     +     ,'***********************'
      RETURN
      END SUBROUTINE HEADING
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRINYS(LAITE,F1C,F1,IXX,IL,IDIR,IK2,IMAX,JMAX,KMAX,N,M)

      USE INTEGERS
      USE NS3CO, ONLY : IN, JN, KN 

      DIMENSION        :: F1(*)
      CHARACTER(LEN=8) :: F1C
      CHARACTER(LEN=5) :: CDI(3)
      CHARACTER(LEN=6) :: CDJ(3)

      INTEGER :: INR, JNR, KNR
       

      DATA CDI   /' I = ',' J = ',' K = '/
      DATA CDJ/', J = ',', K = ',', I = '/
C
C ... FORMATTED OUTPUT OF ARRAY F1 IN Y-DIRECTION, X=IX

      IF(IL <= 0 .OR. IL >= 4) RETURN
      IF(IDIR  < 0 .OR. IDIR  >= 3) THEN
         WRITE(*,*) 'Illigal printing option in prinys IDIR',IDIR,N
         WRITE(*,*) F1C
         WRITE(*,*) 'Exiting...'
         STOP
      ENDIF

      SELECT CASE(IDIR)
      CASE(0)
        INR = IN
        JNR = JN
        KNR = KN
      CASE(1)
        INR = KN
        JNR = IN
        KNR = JN
      CASE(2)
        INR = JN
        JNR = KN
        KNR = IN
      END SELECT

      IX      = MAX(IXX,-1)
      WRITE(LAITE,*)
      NMAXX   = IX + NUMCOL - 1 ! Default is five colmuns of output
      ILA     = IL + IDIR
      GO TO(1,2,3,1,2),ILA

 1    NMAX    = JMAX + 2
      NMAY    = KMAX + 1
      INI     = JNR
      JNI     = KNR
      KNI     = INR
      IK      = MIN0(IMAX+2,IK2)
      ISTEP   = IMAX + 2*INR
      JSTEP   = (IMAX + 2*INR)*(JMAX + 2*JNR)
      KSTEP   = 1

      NMAXX   = MIN0(NMAX,NMAXX)
      GOTO 100

 2    NMAX    = KMAX + 2
      NMAY    = IMAX + 1
      INI     = KNR
      JNI     = INR
      KNI     = JNR
      IK      = MIN0(JMAX+2,IK2)
      ISTEP   = (IMAX + 2*INR)*(JMAX + 2*JNR)
      JSTEP   = 1
      KSTEP   = IMAX + 2*INR

      NMAXX   = MIN0(NMAX,NMAXX)
      GOTO 100

 3    NMAX    = IMAX + 2
      NMAY    = JMAX + 1
      INI     = INR
      JNI     = JNR
      KNI     = KNR
      IK      = MIN0(KMAX+2,IK2)
      ISTEP   = 1
      JSTEP   = IMAX + 2*INR
      KSTEP   = (IMAX + 2*INR)*(JMAX + 2*JNR)

      NMAXX   = MIN0(NMAX,NMAXX)


 100  CONTINUE
      WRITE(LAITE,9102) N,M,CDI(IL),IK,CDJ(IL),IX,NMAXX
 9102 FORMAT(2X,'BLOCK = ',I4,' LEVEL = ',I3,10X,A5,I3,A6,I3,
     +       ',',I3)
      WRITE(LAITE,9100) (F1C,K=IX,NMAXX)
      WRITE(LAITE,*)
 9100 FORMAT(2X,'NODE',1X,10(1X,A8,4X))

      DO 1000 J = -1,NMAY+1
      ILK     =  (JNI+J-1)*JSTEP + (INI-1)*ISTEP + (KNI+IK-1)*KSTEP + 1
      WRITE(LAITE,9200) J,(F1(ILK+I*ISTEP),I=IX,NMAXX)
1000  CONTINUE

 9200 FORMAT(2X,I4,10(1X,E12.5))
 9201 FORMAT(2X,I4,10(1X,E15.8))

      RETURN
      END SUBROUTINE PRINYS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRINYD(LAITE,F1C,F1,IXX,IL,IDIR,IK2,IMAX,JMAX,KMAX,N,M)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: F1(*)

      CHARACTER(LEN=8) :: F1C
      CHARACTER(LEN=5) :: CDI(3)
      CHARACTER(LEN=6) :: CDJ(3)

      DATA CDI   /' I = ',' J = ',' K = '/
      DATA CDJ/', J = ',', K = ',', I = '/
C
C ... FORMATTED OUTPUT OF ARRAY F1 IN Y-DIRECTION, X=IX
C
      IF(IL <= 0 .OR. IL >= 4) RETURN
      IF(IDIR  < 0 .OR. IDIR  >= 3) THEN
         WRITE(*,*) 'Illigal printing option in prinys IDIR',IDIR
         WRITE(*,*) 'Exiting...'
         STOP
      ENDIF

      IX      = MAX(IXX,-1)
      WRITE(LAITE,*)
      NMAXX   = IX + 4
      ILA     = IL + IDIR
      GO TO(1,2,3,1,2),ILA

 1    NMAX    = JMAX + 2
      NMAY    = KMAX + 1
      INI     = JN
      JNI     = KN
      KNI     = IN
      IK      = MIN0(IMAX+2,IK2)
      ISTEP   = IMAX + 2*IN
      JSTEP   = (IMAX + 2*IN)*(JMAX + 2*JN)
      KSTEP   = 1

      NMAXX   = MIN0(NMAX,NMAXX)
      GOTO 100

 2    NMAX    = KMAX + 2
      NMAY    = IMAX + 1
      INI     = KN
      JNI     = IN
      KNI     = JN
      IK      = MIN0(JMAX+2,IK2)
      ISTEP   = (IMAX + 2*IN)*(JMAX + 2*JN)
      JSTEP   = 1
      KSTEP   = IMAX + 2*IN

      NMAXX   = MIN0(NMAX,NMAXX)
      GOTO 100

 3    NMAX    = IMAX + 2
      NMAY    = JMAX + 1
      INI     = IN
      JNI     = JN
      KNI     = KN
      IK      = MIN0(KMAX+2,IK2)
      ISTEP   = 1
      JSTEP   = IMAX + 2*IN
      KSTEP   = (IMAX + 2*IN)*(JMAX + 2*JN)

      NMAXX   = MIN0(NMAX,NMAXX)


 100  CONTINUE
      WRITE(LAITE,9102) N,M,CDI(IL),IK,CDJ(IL),IX,NMAXX
 9102 FORMAT(2X,'BLOCK = ',I4,' LEVEL = ',I3,10X,A5,I3,A6,I3,
     +       ',',I3)
      WRITE(LAITE,9100) (F1C,K=IX,NMAXX)
      WRITE(LAITE,*)
 9100 FORMAT(2X,'NODE',10(1X,A8,4X))

      DO 1000 J = -1,NMAY+1
      ILK     =  (JNI+J-1)*JSTEP + (INI-1)*ISTEP + (KNI+IK-1)*KSTEP + 1
      WRITE(LAITE,9200) J,(REAL(F1(ILK+I*ISTEP),4),I=IX,NMAXX)
1000  CONTINUE

 9200 FORMAT(2X,I4,10(1X,E12.5))
 9201 FORMAT(2X,I4,10(1X,E15.8))
      RETURN
      END SUBROUTINE PRINYD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRINXY(LAITE,F1C,F1,TIME,ICYCLE,IMAX,JMAX,IK,N,M)
C
C ... FORMATTED OUTPUT OF TWO-DIMENSIONAL ARRAY F1
C

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION F1(*)
      CHARACTER*8 F1C

      IMAXP2  = IMAX + 2*IN
      JMAXP2  = JMAX + 2*JN

      WRITE(LAITE,*)
      WRITE(LAITE,9100) F1C,TIME,ICYCLE
      WRITE(LAITE,*)'---------------------------------------------------
     2--------'
      WRITE(LAITE,*)
      WRITE(LAITE,9101) N,M
9100  FORMAT(2X,'FIELD VALUES OF ',A8,2X,'TIME = ',E10.3,' ICYCLE =',I5)
9101  FORMAT(2X,'BLOCK = ',I4,' LEVEL = ',I3)
9200  FORMAT(    7(1X,E9.2))
9300  FORMAT(2X,'IY = ',I3)
      DO 1000 J = 1,JMAX + 1
      IALKU   = (JN+J-1)*IMAXP2 + 1 + IN + (KN+IK-1)*IMAXP2*JMAXP2
      WRITE(LAITE,9300) J
      WRITE(LAITE,9200) (F1(ILK),ILK = IALKU,IALKU+IMAX)

1000  CONTINUE

      RETURN
      END SUBROUTINE PRINXY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE STATE_REPORT(NBLOCK,NPROCE,NBLOCG,NPRO,JSTATE,IPRO,
     2 MULPHL,IPHASE,IPRESC)

      IMPLICIT NONE

      INTEGER :: NBLOCK,NBLOCG,NPRO,IPRO,ISTATE,IPHASE,IPRESC,N,NGL
      INTEGER :: NPROCE(NBLOCG+1,*),JSTATE(NBLOCG,*)
      LOGICAL :: MULPHL
       
      DO N   = 1,NBLOCK
      NGL    = NPROCE(N+1,IPRO)
      ISTATE = JSTATE(NGL,IPHASE)

      SELECT CASE(ISTATE)

      CASE(0)                        !  This phase does not exist
      CASE(1)                        !  PERFECT GAS
         WRITE(4,7167) N,NGL,ISTATE
      CASE(2)                        !  Chemically Reacting Equilibrium Air
         WRITE(4,7167) N,NGL,ISTATE
      CASE(3)                        !  Perfect gas, variable Cp (air)
         WRITE(4,7168) N,NGL,ISTATE
      CASE(4)                        !  SEMI-PERFECT GAS PPR 12.5.1995
         WRITE(4,7170) N,NGL,ISTATE
      CASE(5)
         WRITE(*,*)' WRONG EQUATION OF STATE'
         WRITE(*,*) 'ISTATE = 5 not in use anymore.'
         WRITE(*,*) 'Use ISTATE = 1.  Exiting ...'
         STOP
      CASE(6)                        ! ... WATER joo joo
         WRITE(4,7169) N,NGL,ISTATE
         IF(IPRESC == 0) THEN
          WRITE(6,*)' Do not possible work without pseudo compressible'
         ENDIF
      CASE(7)                        !  WATER joo joo
         WRITE(4,7173) N,NGL,ISTATE
         IF(IPRESC == 0) THEN
          WRITE(*,*)' Do not possible work without pseudo compressible'
         ENDIF
      CASE(8)                        !  Water P = 0.5 ... 1.5 bar
         WRITE(4,7174) N,NGL,ISTATE
         IF(IPRESC == 0) THEN
          WRITE(*,*)' Do not possible work without pseudo compressible'
         ENDIF
      CASE(9)                        !  Steam P = 0.5 ... 1.5 bar
         WRITE(4,7175) N,NGL,ISTATE
      CASE(101)                      !  Iron (cast) T = 0 ... 1000 K
         WRITE(4,7267) N,NGL,ISTATE
      END SELECT

      ENDDO

 7167 FORMAT(3I5,' PERFECT GAS AS A FUNCTION OF RHO AND E')
 7168 FORMAT(3I5,' CHEMICALLY REACTING EQUILIBRIUM AIR')
 7169 FORMAT(3I5,' WATER AS A FUNCTION OF P AND T')
 7170 FORMAT(3I5,' SEMI PERFECT GAS')
 7171 FORMAT(3I5,' PERFECT GAS A FUNCTION OF P AND E')
 7172 FORMAT(3I5,' EQUATION OF STATE DOES NOT EXIST'//'  Exiting ...')
 7173 FORMAT(3I5,' WATER AS A FUNCTION OF P AND T. MODIFIED VISCOSITY!')
 7174 FORMAT(3I5,' WATER AS A FUNCTION OF P AND T. RANGE 0.5 - 1.5 BAR')
 7175 FORMAT(3I5,' STEAM AS A FUNCTION OF P AND T. RANGE 0.5 - 1.5 BAR')

 7267 FORMAT(3I5,' CAST IRON AS A FUNCTION T. RANGE 0 - 1000 K')

      END SUBROUTINE STATE_REPORT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
      SUBROUTINE BLOCK_REPORT(NBLOCK,NPROCE,NBLOCG,IPRO,BLKS,
     &                        FRSDEN,FRSVEL,FRSTEM,FRSPRE,FRSEG,FRSMUT,
     &                        RMUINI,T0,TLOLIM,TUPLIM,RKLIM,FRSLEN,
     &                        EPSLIM,TEMINI,ITURB,JSTATE,FRSSSP,ARTSSP,
     &                        FRSVIS,FRSSIE,SIEINI,ROTAT,ROTANG,OMEGA,
     &                        OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,
     &                        MGRID,INTERI,INTERJ,INTERK,MGRIDA,
     &                        MULPHL,IGRID,REFVEL,IPRESC,GROUND,
     &                        CAVLEV,TWO_FLUIDL)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : PII,RAD2DEG
      USE FLIGHT   , ONLY : OSKU

      IMPLICIT NONE

      INTEGER :: NBLOCK,NBLOCG,NPRO,IPRO,N,NGL,I,ITURB,MPCASE,IPRESC,
     &        IEVAP
      INTEGER,DIMENSION(NBLOCG+1,*) :: NPROCE, IGRID
      INTEGER,DIMENSION(NBLOCG,  *) :: JSTATE
      INTEGER,DIMENSION(*) :: MGRID,INTERI,INTERJ,INTERK,MGRIDA
      REAL :: FRSE,FRSDEN,FRSVEL,FRSTEM,FRSPRE,FRSEG,FRSMUT,RMUINI,T0,
     &        TLOLIM,TUPLIM,RKLIM,FRSLEN,EPSLIM,TEMINI,FRSSSP,ARTSSP,
     &        FRSVIS,FRSSIE,SIEINI,PSAT,CAVNO,FRSMU,FRSKE,RKPKE,RMACH,
     &        FRSMUA,REFVEL,GROUND,CAVLEV
      REAL,DIMENSION(*) :: ROTAT,ROTANG,OMEGA,OMEGAX,OMEGAY,OMEGAZ,
     &                     CENAX,CENAY,CENAZ
      REAL :: DRSPRE,DRSTEM,DPSAT
      CHARACTER(LEN=3) :: CHRLN3
      CHARACTER(LEN=6) :: CHRLN6
      CHARACTER(LEN=7) :: FRSALFA1
      CHARACTER(LEN=30):: MATERIAL1
      CHARACTER(LEN=12):: STATE1
      CHARACTER(LEN=8),DIMENSION(8) :: FLUXTYB = (/
     &    "Van Leer","  Roe   ","  WPS   ","  AUSM  ",
     &    "  INCO  ","Roe RSM ","M-PHASE ","Kurganov" /)
      CHARACTER(LEN=14),DIMENSION(4) :: SOL_T = (/
     &    "Density based ","Pressure corr.",
     &    "Preconditioned","Conduction    " /)
      LOGICAL :: MULPHL,NOTICE,TWO_FLUIDL
      TYPE(BLOCKS),DIMENSION(*) :: BLKS

! ... FORMAT statements

 9090 FORMAT(84("=")/)
 9091 FORMAT(84("="))
 9092 FORMAT(84("-"))
c 9101 FORMAT(1X,A3,1X,3F9.4,4(1X,E11.5))
c 9103 FORMAT(1X,A3,1X,3F9.3,4(1X,E11.5))
c 9105 FORMAT(1X,A3,2F8.2,1X,E11.5,1X,A11,4(1X,E11.5))
 9101 FORMAT(2X,A6,1X,3F9.2,4(1X,E11.4))
 9103 FORMAT(2X,A6,3F9.2,4(1X,E11.4))
 9105 FORMAT(2X,A6,1X,2F9.2,1X,E11.4,1X,A11,4(1X,E11.4))
c 9106 FORMAT(1X,I3,1X,2F9.5,1X,4(E11.5,1X),7X,L1)
c 9107 FORMAT(I3,1X,F10.2,1X,3(1X,F8.5),2X,3(1X,E11.5))
c 9109 FORMAT(I3,1X,F10.2,6X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)
c 9111 FORMAT(I3,5X,A10,1X,I2," : ",3I4," <=> ",3A10,1X,A7)
 9106 FORMAT(I5,2E10.3,1X,4(E11.4,1X),7X,L1)          ! F changed to E in the second entry. Matti Palin 19.01.2021
 9107 FORMAT(I5,F9.2,1X,3(1X,F9.2),2X,3(1X,E11.4))
 9109 FORMAT(I5,F9.2,6X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)
 9111 FORMAT(I5,3X,A10,1X,I2," : ",3I4," <=> ",3A10,1X,A7)
 9011 FORMAT(I5,3X,A10,1X,I2," : ", I4,"       <=>       ", A10,5X,A7)
 9012 FORMAT(I5,3X,A10,1X,I2," : ",2I4,"    <=>    ",2A10,5X,A7)
 9013 FORMAT(I5,3X,A10,1X,I2," : ",A12," <=> ",A30,1X,A7)
 9112 FORMAT(I5,6(3X,I5))
 9113 FORMAT(I5,A9,"(",I1,") ",3I5,3X,I5,1X,I5,3X,A14,I2,2X,I2,1X,I2,
     &     1X,I2,2X,I2)
 9114 FORMAT(I5,F8.3,3X,F8.6,3X,I1,8X,L1,7X,L1,3X,I5,1X,I5,1X,I5,1X,
     & I5,2X,L5)
 9115 FORMAT(I5,2X,I3,1X,F10.3,2X,E12.5,2X,E12.5,2X,E12.5,2(1X,F9.2))
 9116 FORMAT(I5,A9,"(",I1,") ",3I5,4X,4F8.2)

 9450 FORMAT("   BLOCK    FRSDEN   REFVEL   FRSTEM   FRSPRE",
     &       "      FRSENE      FRSMUT      RMUINI")
 9460 FORMAT("   BLOCK    T0      TLOLIM   TUPLIM     RKLIM",
     &       "     FRSLEN     ", A8 ,"     TEMINI")
 9470 FORMAT("   BLOCK    FRSSSP  ARTSSP   PLOLIM       PUPLIM",
     &       "     FRSVIS      FRSSIE      SIEINI")
 9510 FORMAT("GRAVITY VECTOR: [  ",3(F7.4,2X),"].")
 9540 FORMAT("BLOCK  SOLUTION     NUMBER OF PHASES ",
     &       "       STATE MATERIALS           FRSALFA(1)")
 9541 FORMAT("BLOCK    TURMUL  DRAGMUL LIFTMUL DISPER   VMASS  WFORCE")
c     &       "       STATE MATERIALS           FRSALFA(1)")
 9550 FORMAT("BLOCK  OMEGA (1/s)   VECTOR OF ROTATION          ",
     &       "         CENTER POINT         ")
 9560 FORMAT("BLOCK  EVAP.TYPE  CAVITATION NUMBER ",
     &       "SATURATION PRESSURE   FRSTEM       TEMINI")
 9561 FORMAT("BLOCK EVAP.TY CAV.NO. ",
     &       "SAT. PRESSURE    FRSTEM       FRSPRE",
     &       "       CHLCAV   REFCAV")
9462  FORMAT('BLOCK','    Tu(FRS) ',' MUT(FRS)',' TURLIM ',
     + '     RKLIM/K.e. ','K.e.(FRS)[J/m3] ','MUT(max)','     KATOL')
9463  FORMAT(/'DIMENSIONLESS VARIABLES:'/)
9464  FORMAT('Warning: As Tu is lower than 0.01, free-stream eddy-',
     + 'viscosity should be smaller than',/'MUT(max). In this case k-e',
     + ' model may occasionally produce non-physical (high) values'/
     + 'for the eddy-viscosity.')
 9570 FORMAT("BLOCK  OMEGA (1/s)   CURRENT ANGLE  (ROTANG)     ",
     &       "     INITIAL ANGLE  (ROTAT)    ")
 9575 FORMAT(25X,"(deg)",7X,"(rad)",16X,"(deg)",7X,"(rad)")
 9580 FORMAT('BLOCK','  Flux type ',' Interpolation(IJK) ',
     +  ' MGRID MGRIDA   Solution  LUSGS  IGRID IDIFF')
 9582 FORMAT('BLOCK',' Flux type  ','Void interpolation(IJK)',
     +  '      FREDIF(IJK)  ','     RINTF  ')
 9576 FORMAT('BLOCK','  Flux type ',' Interpolation(IJK) ',
     +  ' MGRID MGRIDA   Solution  LUSGS  IGRID IDIFF')
 9581 FORMAT('BLOCK','  ALFAP      ALFAU  PREVEL SKEWCORR COMPCORR',
     +  ' INLRC OUTRC ICERO ITERMP FLUXCORR')

! ... Executables begin

      WRITE(4,"(/A/36('*')/)") "INFORMATION ABOUT BLOCK PROPERTIES :"

c      IF(.NOT.MULPHL) THEN

      WRITE(4,9450)
      WRITE(4,9091)
      WRITE(4,9101) "Def",FRSDEN,REFVEL,FRSTEM,
     &    FRSPRE,FRSEG,FRSMUT,RMUINI
      WRITE(4,9092)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       FRSE = BLKS(NGL)%FRSDEN*(BLKS(NGL)%FRSSIE+.5*BLKS(NGL)%FRSVEL**2)
c       WRITE(CHRLN3,"(I3)") NGL
c       WRITE(CHRLN3,"(I6)") NGL
c       WRITE(4,9101) CHRLN3,BLKS(NGL)%FRSDEN,BLKS(NGL)%REFVEL,
       WRITE(CHRLN6,"(I6)") NGL
       WRITE(4,9101) CHRLN6,BLKS(NGL)%FRSDEN,BLKS(NGL)%REFVEL,
     &     BLKS(NGL)%FRSTEM,BLKS(NGL)%FRSPRE,FRSE,
     &     BLKS(NGL)%FRSMUT,BLKS(NGL)%RMUINI
      END DO
      WRITE(4,9090)

c      ENDIF ! .NOT.MULPHL


      SELECT CASE(ITURB)
      CASE(9)   ; WRITE(4,9460) " RNUTLIM"
      CASE(6:8) ; WRITE(4,9460) "OMEGALIM"
      CASE(3:5) ; WRITE(4,9460) " EPSLIM "
      END SELECT
      WRITE(4,9091)
      WRITE(4,9103) "Def",T0,TLOLIM,TUPLIM,RKLIM,FRSLEN,EPSLIM,TEMINI
      WRITE(4,9092)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
c       WRITE(CHRLN3,"(I3)") NGL
       WRITE(CHRLN6,"(I6)") NGL
c       WRITE(4,9103) CHRLN3,BLKS(NGL)%T0,   BLKS(NGL)%TLOLIM,
       WRITE(4,9103) CHRLN6,BLKS(NGL)%T0,   BLKS(NGL)%TLOLIM,
     &     BLKS(NGL)%TUPLIM,BLKS(NGL)%RKLIM,BLKS(NGL)%FRSLEN,
     &     BLKS(NGL)%EPSLIM,BLKS(NGL)%TEMINI
      END DO

      SELECT CASE(ITURB)
      CASE(3:5)
         WRITE(4,9463) ; WRITE(4,9462) ; WRITE(4,9091)
         NOTICE = .FALSE.
         DO N   = 1,NBLOCK
         NGL    = NPROCE(N+1,IPRO)
         FRSMU  = BLKS(NGL)%FRSMUT/BLKS(NGL)%FRSVIS
         FRSKE  = .5*BLKS(NGL)%FRSDEN*BLKS(NGL)%REFVEL**2
         RKPKE  = BLKS(NGL)%RKLIM/FRSKE
         RMACH  = BLKS(NGL)%REFVEL/BLKS(NGL)%FRSSSP
         FRSMUA = 1.E5*BLKS(NGL)%FRSTUR**2
         WRITE(4,9106) NGL,BLKS(NGL)%FRSTUR,FRSMU,
     +   BLKS(NGL)%TURLIM,RKPKE,FRSKE,FRSMUA,BLKS(NGL)%KATOL
         IF(BLKS(NGL)%FRSTUR <= 1.E-2 .AND. FRSMU > FRSMUA) 
     +   NOTICE = .TRUE.
         ENDDO
      CASE(6:7,9)  !Following should be k-omega based !!!!
         WRITE(4,9463) ; WRITE(4,9462) ; WRITE(4,9091)
         NOTICE = .FALSE.
         DO N   = 1,NBLOCK
         NGL    = NPROCE(N+1,IPRO)
         FRSMU  = BLKS(NGL)%FRSMUT/BLKS(NGL)%FRSVIS
         FRSKE  = .5*BLKS(NGL)%FRSDEN*BLKS(NGL)%REFVEL**2
         RKPKE  = BLKS(NGL)%RKLIM/FRSKE
         RMACH  = BLKS(NGL)%REFVEL/BLKS(NGL)%FRSSSP
         FRSMUA = 1.E5*BLKS(NGL)%FRSTUR**2
         WRITE(4,9106) NGL,BLKS(NGL)%FRSTUR,FRSMU,
     +   BLKS(NGL)%TURLIM,RKPKE,FRSKE,FRSMUA,BLKS(NGL)%KATOL
         ENDDO
      END SELECT
      IF(NOTICE) THEN
         WRITE(4,*)
         WRITE(4,9464)
         NOTICE = .FALSE.
      ENDIF

      WRITE(4,9090)

      WRITE(4,9470)
      WRITE(4,9091)
      WRITE(4,9105) "Def",FRSSSP,ARTSSP,0.01*FRSPRE,"  Unused   ",
     &    FRSVIS,FRSSIE,SIEINI
      WRITE(4,9092)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       FRSE = BLKS(NGL)%FRSDEN*(BLKS(NGL)%FRSSIE+.5*BLKS(NGL)%FRSVEL**2)
*       WRITE(CHRLN3,"(I3)") NGL
*       WRITE(4,9105) CHRLN3,BLKS(NGL)%FRSSSP,BLKS(NGL)%ARTSSP,
       WRITE(CHRLN6,"(I6)") NGL
       WRITE(4,9105) CHRLN6,BLKS(NGL)%FRSSSP,BLKS(NGL)%ARTSSP,
     &     0.01*BLKS(NGL)%FRSPRE,"  Unused   ",BLKS(NGL)%FRSVIS,
     &     BLKS(NGL)%FRSSIE,BLKS(NGL)%SIEINI
      END DO
      
      WRITE(4,9090)

      WRITE(4,9540)
      WRITE(4,9091)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       IF (MULPHL) THEN
        WRITE(FRSALFA1,"(F7.5)") BLKS(NGL)%FRSALFA(1)
        WRITE(MATERIAL1,"(3(A10))")(BLKS(NGL)%MATERIAL(I),I=1,NPHASES)
        WRITE(STATE1,"(3(I4))")(JSTATE(NGL,I),I=1,NPHASES)
       ELSE
        FRSALFA1 = "  N/A  "
        WRITE(MATERIAL1,"(3(A10))")(BLKS(NGL)%MATERIAL(I),I=1,NPHASES)
        WRITE(STATE1,"(3(I4))")(JSTATE(NGL,I),I=1,NPHASES)
       END IF
       WRITE(4,9013) NGL,BLKS(NGL)%SOLUTION_TYPE,BLKS(NGL)%NPHASE,
     &     STATE1,MATERIAL1,FRSALFA1
      END DO

      IF(TWO_FLUIDL) THEN

      WRITE(4,9090)

      WRITE(4,9541)
      WRITE(4,9091)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       WRITE(4,9112) NGL,BLKS(NGL)%TURMUL,BLKS(NGL)%DRAGMUL,
     & BLKS(NGL)%LIFTMUL,BLKS(NGL)%DISPERMUL,BLKS(NGL)%VMASS,
     & BLKS(NGL)%WFORCEMUL
      END DO

      ENDIF ! TWO_FLUIDL

      WRITE(4,9090)

      IF (MULPHL) THEN

       WRITE(4,9561)
       WRITE(4,9091)

!       OPEN(1974,FILE="MPCASES")
         
       DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)

        MPCASE = BLKS(NGL)%MPCASE
        IEVAP  = BLKS(NGL)%EVAP_TYPE

        IF(IEVAP >= 3 .AND. IEVAP /= 6 .AND. MPCASE == 1) THEN
            WRITE(*,*)' MPCASE should be ge 2 as temperature based'
            WRITE(*,*)' mass transfer model is applied. Exiting..'
            WRITE(13,*)' MPCASE should be ge 2 as temperature based'
            WRITE(13,*)' mass transfer model is applied. Exiting..'
            STOP
         ENDIF

        SELECT CASE(MPCASE)
        CASE(0) ! Use existing FRSTEM and FRSDEN
         CAVNO  = BLKS(NGL)%CAVNO
         BLKS(NGL)%DFRSTEM = BLKS(NGL)%FRSTEM
         BLKS(NGL)%DFRSPRE = BLKS(NGL)%FRSPRE
         BLKS(NGL)%DFRSDEN = BLKS(NGL)%FRSDEN
         DRSTEM = BLKS(NGL)%FRSTEM
         DRSPRE = BLKS(NGL)%FRSPRE
         CALL DPSATLOG(DPSAT,DRSTEM,JSTATE(NGL,1))
         DPSAT  = DPSAT - DRSPRE
        CASE(1) ! Determine FRSTEM and FRSDEN
!         CALL FCAVNO(BLKS(NGL)%CAVNO, BLKS(NGL)%FRSTEM,
!     &               BLKS(NGL)%FRSPRE,BLKS(NGL)%FRSDEN,
!     &               BLKS(NGL)%FRSVEL,BLKS(NGL)%DFRSDEN,
!     &               BLKS(NGL)%DFRSTEM,JSTATE(NGL,1),NGL)
         IF(BLKS(NGL)%CAVNO < 1.E-7) THEN
            WRITE(*,*)' You have given a zero cavitation number.'
            WRITE(*,*)' Free stream values cannot be found. Exiting..'
            WRITE(13,*)' You have given a zero for cavitation number.'
            WRITE(13,*)' Free stream values cannot be found. Exit.'
            STOP
         ENDIF
         CALL FTEMPS(BLKS(NGL)%CAVNO, BLKS(NGL)%FRSTEM,
     &               BLKS(NGL)%FRSPRE,BLKS(NGL)%FRSDEN,
     &               BLKS(NGL)%REFVEL,BLKS(NGL)%DFRSDEN,
     &               BLKS(NGL)%DFRSTEM,JSTATE(NGL,1),NGL,
     &               GROUND,CAVLEV)
         CAVNO  = BLKS(NGL)%CAVNO
         DRSTEM = BLKS(NGL)%DFRSTEM
         DRSPRE = BLKS(NGL)%FRSPRE
         CALL DPSATLOG(DPSAT,DRSTEM,JSTATE(NGL,1))
         DPSAT  = DPSAT - DRSPRE
!         CALL PSATF(PSAT,BLKS(NGL)%FRSTEM,8)
!         CALL DPSATF(DPSAT,DRSTEM,DRSPRE,8)
!         CAVNO = -2.0*DPSAT/(BLKS(NGL)%DFRSDEN*BLKS(NGL)%FRSVEL**2)
!         CAVNO = 2.0*REAL(DRSPRE-DPSAT)
!     &       /(BLKS(NGL)%FRSDEN*BLKS(NGL)%FRSVEL**2)
        CASE(2) ! Determine saturation pressure and cavitation #
         CALL FSIGMA(BLKS(NGL)%CAVNO,  BLKS(NGL)%FRSTEM,
     &               BLKS(NGL)%FRSPRE, BLKS(NGL)%FRSDEN,
     &               BLKS(NGL)%REFVEL, BLKS(NGL)%DFRSDEN,
     &               BLKS(NGL)%DFRSTEM,DPSAT,JSTATE(NGL,1),NGL,
     &               GROUND,CAVLEV)
         CAVNO  = BLKS(NGL)%CAVNO
         DRSTEM = BLKS(NGL)%DFRSTEM
         DRSPRE = BLKS(NGL)%FRSPRE
         DPSAT  = DPSAT - DRSPRE
        CASE(3) ! Determine FRSPRE and FRSDEN
         IF(BLKS(NGL)%CAVNO < 1.E-7) THEN
            WRITE(*,*)' You have given a zero cavitation no. MPCASE?'
            WRITE(*,*)' Free stream valueses cannot be found. Exiting..'
            WRITE(13,*)' You have given a zero for cavitation number.'
            WRITE(13,*)' Free stream valueses cannot be found. Exit.'
            STOP
         ENDIF
         CALL FPRESS(BLKS(NGL)%CAVNO,  BLKS(NGL)%FRSTEM,
     &               BLKS(NGL)%FRSPRE, BLKS(NGL)%FRSDEN,
     &               BLKS(NGL)%REFVEL, BLKS(NGL)%DFRSDEN,
     &               BLKS(NGL)%DFRSTEM,DPSAT,JSTATE(NGL,1),NGL,
     &               GROUND,CAVLEV)
         CAVNO  = BLKS(NGL)%CAVNO
         FRSPRE = BLKS(NGL)%FRSPRE
         DRSTEM = BLKS(NGL)%DFRSTEM
         DRSPRE = BLKS(NGL)%FRSPRE
         DPSAT  = DPSAT - DRSPRE
        CASE(4) ! Determine FRSTEM and FRSDEN as temperature based
         IF(BLKS(NGL)%CAVNO < 1.E-7) THEN
            WRITE(*,*)' You have given a zero cavitation no. MPCASE?'
            WRITE(*,*)' Free stream valueses cannot be found. Exiting..'
            WRITE(13,*)' You have given a zero for cavitation number.'
            WRITE(13,*)' Free stream valueses cannot be found. Exit.'
            STOP
         ENDIF
         CALL FTEMPT(BLKS(NGL)%CAVNO,  BLKS(NGL)%FRSTEM,
     &               BLKS(NGL)%FRSPRE, BLKS(NGL)%FRSDEN,
     &               BLKS(NGL)%REFVEL, BLKS(NGL)%DFRSDEN,
     &               BLKS(NGL)%DFRSTEM,DPSAT,JSTATE(NGL,1),NGL,
     &               GROUND,CAVLEV)
         CAVNO  = BLKS(NGL)%CAVNO
         FRSPRE = BLKS(NGL)%FRSPRE
         DRSTEM = BLKS(NGL)%DFRSTEM
         DRSPRE = BLKS(NGL)%FRSPRE
         DPSAT  = DPSAT - DRSPRE
        END SELECT

!        IF (BLKS(NGL)%CAVNO > 0.) THEN ! Determine FRSTEM and FRSDEN
!        ELSE ! Determine saturation pressure and cavitation number
!        END IF

        BLKS(NGL)%CAVNO   = CAVNO
        BLKS(NGL)%DFRSPRE = BLKS(NGL)%FRSPRE
        IF (BLKS(NGL)%SOLUTION_TYPE == "CAVIT" .OR.
     &      BLKS(NGL)%SOLUTION_TYPE == "MULTI") THEN
         BLKS(NGL)%TEMINI = BLKS(NGL)%FRSTEM
         WRITE(4,9115) NGL,BLKS(NGL)%EVAP_TYPE,CAVNO,DPSAT,
     &       BLKS(NGL)%FRSTEM,BLKS(NGL)%FRSPRE,BLKS(NGL)%CHLCAV,
     &       BLKS(NGL)%REFCAV
        ELSE
         WRITE(4,9115) NGL,BLKS(NGL)%EVAP_TYPE,CAVNO,DPSAT + DRSPRE,
     &       BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,BLKS(NGL)%CHLCAV,
     &       BLKS(NGL)%REFCAV
        END IF

!        WRITE(1974,"(/2I5,2E15.6/8E15.6/)")
!     &      NGL,JSTATE(NGL,1),BLKS(NGL)%CAVNO,BLKS(NGL)%TEMINI,
!     &      BLKS(NGL)%FRSTEM,BLKS(NGL)%DFRSTEM,BLKS(NGL)%FRSPRE,
!     &      BLKS(NGL)%FRSDEN,BLKS(NGL)%DFRSDEN,BLKS(NGL)%FRSVEL

       END DO

       WRITE(4,9090)

!       CLOSE(1974)

       SELECT CASE(MPCASE)
       CASE(1) ; WRITE(4,"(2(/A))")
     &     "Free-stream temperature has been, at least in one block,",
     &     "re-evaluated from the cavitation number "//
     &     "and pressure."
       CASE(2) ; WRITE(4,"(2(/A))")
     &     "Cavitation number has been, at least in one block,",
     &     "re-evaluated from the free-stream temperature "//
     &     "and pressure."
       CASE(3) ; WRITE(4,"(2(/A))")
     &     "Free-stream pressure has been, at least in one block,",
     &     "re-evaluated from the free-stream temperature "//
     &     "and cavitation number."
       CASE(4) ; WRITE(4,"(2(/A))")
     &     "Free-stream temperature has been, at least in one block,",
     &     "re-evaluated from the cavitation number. "
       END SELECT

      WRITE(4,*)
      WRITE(4,9450)
      WRITE(4,9091)
      WRITE(4,9101) "Def",FRSDEN,REFVEL,FRSTEM,
     &    FRSPRE,FRSEG,FRSMUT,RMUINI
      WRITE(4,9092)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       FRSE = BLKS(NGL)%FRSDEN*(BLKS(NGL)%FRSSIE+.5*BLKS(NGL)%FRSVEL**2)
c       WRITE(CHRLN3,"(I3)") NGL
c       WRITE(4,9101) CHRLN3,BLKS(NGL)%FRSDEN,BLKS(NGL)%REFVEL,
       WRITE(CHRLN6,"(I6)") NGL
       WRITE(4,9101) CHRLN6,BLKS(NGL)%FRSDEN,BLKS(NGL)%REFVEL,
     &     BLKS(NGL)%FRSTEM,BLKS(NGL)%FRSPRE,FRSE,
     &     BLKS(NGL)%FRSMUT,BLKS(NGL)%RMUINI
      END DO
      WRITE(4,9090)


      END IF
         
      WRITE(4,"(//A/31('*')/)") "INFORMATION ABOUT GRID ROTATION"

      WRITE(4,9550)
      WRITE(4,9091)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       WRITE(4,9107) NGL,OMEGA(NGL),
     &     OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),
     &     CENAX(NGL),CENAY(NGL),CENAZ(NGL)
       IF(IGRID(NGL,1) == 10) THEN
          OSKU(IGRID(NGL,3))%VX =OMEGAX(NGL)
          OSKU(IGRID(NGL,3))%VY =OMEGAY(NGL)
          OSKU(IGRID(NGL,3))%VZ =OMEGAZ(NGL)
       ENDIF
      END DO

      WRITE(4,9090)

      WRITE(4,9570)
      WRITE(4,9575)
      WRITE(4,9091)
      DO N = 2,NBLOCK+1 ; NGL = NPROCE(N,IPRO)
       WRITE(4,9109) NGL,OMEGA(NGL),
     &     ROTANG(NGL)*RAD2DEG,ROTANG(NGL),
     &      ROTAT(NGL)*RAD2DEG, ROTAT(NGL)
      END DO

      WRITE(4,9090)

      WRITE(4,"(/A/26('*')/)") "INFORMATION ABOUT NUMERICS"

      WRITE(4,9580)
      WRITE(4,9091)
      DO N = 1,NBLOCK ; NGL = NPROCE(N+1,IPRO)
       WRITE(4,9113) NGL,FLUXTYB(BLKS(NGL)%IFLUX+1),BLKS(NGL)%IFLUX,
     &     INTERI(N),INTERJ(N),INTERK(N),MGRID(N),MGRIDA(N),
     &     SOL_T(BLKS(NGL)%IPRESC+1),BLKS(NGL)%LUSGS,IGRID(NGL,1),
     &     IGRID(NGL,2),IGRID(NGL,3),BLKS(NGL)%IDIFF
      END DO

      WRITE(4,9090)

      IF(IPRESC /= 0) THEN

      WRITE(4,9581)
      WRITE(4,9091)
      DO N = 1,NBLOCK ; NGL = NPROCE(N+1,IPRO)
       WRITE(4,9114)NGL,BLKS(NGL)%ALFAP,BLKS(NGL)%ALFAU,BLKS(NGL)%PREVEL
     & ,BLKS(NGL)%SKEWCORR,BLKS(NGL)%COMPCORR,BLKS(NGL)%INLRC,
     &  BLKS(NGL)%OUTRC,BLKS(NGL)%ICERO,BLKS(NGL)%ITERMP,
     &  BLKS(NGL)%FLUXCORR
      END DO

      WRITE(4,9090)

      ENDIF

      IF(MULPHL) THEN

      WRITE(4,9582)
      WRITE(4,9091)
      DO N = 1,NBLOCK ; NGL = NPROCE(N+1,IPRO)
       WRITE(4,9116)NGL,FLUXTYB(BLKS(NGL)%IFLUX+1),BLKS(NGL)%IFLUX,
     &     BLKS(NGL)%INTERI,BLKS(NGL)%INTERJ,BLKS(NGL)%INTERK,
     &     BLKS(NGL)%FREDIFI,BLKS(NGL)%FREDIFJ,BLKS(NGL)%FREDIFK,
     &     BLKS(NGL)%RINTF
      END DO

      WRITE(4,9090)

      ENDIF

      END SUBROUTINE BLOCK_REPORT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLIGHT_REPORT(OSKU,XCG,YCG,ZCG,PSIR,THETAR,PHIR,
     & NGRIFL,XCGI,YCGI,ZCGI,PSIRI,THETARI,PHIRI)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : PII, RAD2DEG
      USE NS3CO, ONLY : ROTORW,REFVEL
      USE INTEGERS, ONLY : IPRO
      USE FLIGHT , ONLY : RIN,ROUT,THRUST,TORQUE,ADV,ROTA1,ROTB1,ROTA1I,
     &     ROTB1I,CONEA,CONEAI,ACTUA,VTIP,CDBLADE,CFBLADE,SIGMA,
     &     IFA,IFT,IFR,NSERIES,TRIMA,TRIMAI,DRAUGHT,
     &     DRAUGHTI,DAMPC1,DAMPC2,IGRSLAVE

      USE BLADE_VARIABLES, ONLY : SHAFT,RHINGE,RCG,RIBE,RIZE,THETACOL,
     &     THCYCLON,THCYCLAT,BETAHCOL,BECYCLON,BECYCLAT,ZETAHCOL,
     &     KDELTA3,RKZE,RKBE,RD,RKD,CD,RCZE,RCBE

      IMPLICIT NONE

      INTEGER :: IGR,NGRIFL

      TYPE(SIXDOF) :: OSKU(*)
    
      REAL, DIMENSION(*) :: XCG,YCG,ZCG,XCGI,YCGI,ZCGI,
     &     PSIR,THETAR,PHIR,PSIRI,THETARI,PHIRI
      
      IF(IPRO == 1) THEN 
         GOTO 9000
      ELSE
         RETURN
      ENDIF
     
 9000 CONTINUE

      WRITE(4,*)
      WRITE(4,*) 'INFORMATION ABOUT THE MOVING OBJECTS'
      WRITE(4,*) '************************************'

      DO IGR = 1,NGRIFL

      IF (OSKU(IGR)%CHAR == 'F') THEN
9650    FORMAT(/,' FLYING OBJECTS NO.',I2/,/,
     & ' MASS  = ',F9.3, '[kg]','  AREF   = ',F9.3, '[m²]',
     & '   CHLREF  = ',F9.3, '[m]'/,/,         
     & ' IX    = ',F9.3, '[Nm]','  IY     = ',F9.3, '[Nm]',
     & '   IZ      = ',F9.3, '[Nm]'/,
     & ' IXY   = ',F9.3, '[Nm]','  IXZ    = ',F9.3, '[Nm]',
     & '   IYZ     = ',F9.3, '[Nm]'/,
     & ' DAMPN = ',F9.5,'      DAMPT  = ',F9.5/,/,
     & ' TOL   = ',F9.5,'      H1     = ',F9.5, '[s]',
     & '    HMIN    = ',F9.5, '[s]'/
     & ' TDEL  = ',F9.5, '[s]','   TSTEER (=tmax+dt) = ',F9.5, '[s]'/,/,
     & ' SPEED = ',F9.3, '[m/s]',' ALPHAF = ',F9.3, '[deg]',
     & ' BETAF    = ',F9.3, '[deg]'/,
     & ' GAMMAF= ',F9.3, '[deg]',' BANKF  = ',F9.3, '[deg]'/,/,
     & ' XCG   = ',F9.3,'[m]','   YCG    = ',F9.3,'[m]',
     & '   ZCG      = ',F9.3,'[m]'/,
     & ' VX    = ',F9.3, '[m/s]',' VY     = ',F9.3, '[m/s]',
     & ' VZ       = ',F9.3, '[m/s]'/,/,
     & ' PSIR (att -y) = ',F9.3, '[deg]',
     & '   THETAR (att -z) = ',F9.3, '[deg]'/,
     & ' PHIR (att -x) = ',F9.3, '[deg]',
     & '   ROTX            = ',F9.3, '[rad/s]'/,
     & ' ROTY          = ',F9.3, '[rad/s]',
     & ' ROTZ            = ',F9.3, '[rad/s]'/,/)

        WRITE(4,9650)IGR,OSKU(IGR)%RMASS,OSKU(IGR)%AREF,
     &   OSKU(IGR)%CHLREF,
     &   OSKU(IGR)%IX,OSKU(IGR)%IY,OSKU(IGR)%IZ,
     &   OSKU(IGR)%IXY,OSKU(IGR)%IXZ,OSKU(IGR)%IYZ,
     &   OSKU(IGR)%DAMPN,OSKU(IGR)%DAMPT, 
     &   OSKU(IGR)%TOL,OSKU(IGR)%H1,OSKU(IGR)%HMIN,
     &   OSKU(IGR)%TDEL,OSKU(IGR)%TSTEER,   
     &   OSKU(IGR)%SPEED,OSKU(IGR)%ALPHA,
     &   OSKU(IGR)%BETA,OSKU(IGR)%GAMMA,OSKU(IGR)%BANK,
     &   XCG(IGR),YCG(IGR),ZCG(IGR),
     &   OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     &   PSIR(IGR)*RAD2DEG,THETAR(IGR)*RAD2DEG,PHIR(IGR)*RAD2DEG,      
     &   OSKU(IGR)%ROTX,OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ        

      ELSEIF (OSKU(IGR)%CHAR == 'A') THEN
9660    FORMAT(/,' ACTUATOR DISK NO.',I2/,/,
     & ' THRUST = ',F9.2, '[N]','   TORQUE  = ',F9.2, '[Nm]',
     & '  ADV     = ',F9.5, ' ',/ 
     & ' RIN    = ',F9.3, '[m]',         
     & '   ROUT    = ',F9.3, '[m]','   AREA    = ',F9.3, '[m²]'/,
     & ' VTIP   = ',F9.3, '[m/s]',' CFBLADE = ',F9.5, '  ',
     & '    CDBLADE = ',F9.5, '  ' /,
     & ' SIGMA  = ',F9.5, '  '/,/,
     & ' XCGI   = ',F9.3, '[m]','   YCGI    = ',F9.3, '[m]',
     & '   ZCGI    = ',F9.3, '[m]'/,
     & ' XCG    = ',F9.3, '[m]','   YCG     = ',F9.3, '[m]',
     & '   ZCG     = ',F9.3, '[m]'/,
     & ' ROTA1I = ',F9.3, '[deg]',' ROTB1I  = ',F9.3, '[deg]',
     & ' CONEAI  = ',F9.3, '[deg]'/,
     & ' ROTA1  = ',F9.3, '[deg]',' ROTB1   = ',F9.3, '[deg]',
     & ' CONEA   = ',F9.3, '[deg]'/,/,
     & ' IFA     = ',I2,'                IFT  = ',I2,
     & '                IFR  = ',I2,/
     & ' NSERIES = ',I2,/,/)

        WRITE(4,9660)IGR,THRUST(IGR),TORQUE(IGR),ADV(IGR),
     &       RIN(IGR),ROUT(IGR),
     &       PII*(ROUT(IGR)**2-RIN(IGR)**2),
     &       VTIP(IGR),CFBLADE(IGR),CDBLADE(IGR),SIGMA(IGR),
     &       XCGI(IGR),YCGI(IGR),ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),
     &       ROTA1I(IGR)*RAD2DEG,ROTB1I(IGR)*RAD2DEG,
     &       CONEAI(IGR)*RAD2DEG,ROTA1(IGR)*RAD2DEG,
     &       ROTB1(IGR)*RAD2DEG,CONEA(IGR)*RAD2DEG,
     &       IFA(IGR),IFT(IGR),IFR(IGR),NSERIES(IGR)

      ELSEIF (OSKU(IGR)%CHAR == 'B') THEN
9661    FORMAT(/,' ACTUATOR DISK NO.',I2
     & ',   CONNECTED TO SHIP NO. ',I2/,/,
     & ' THRUST = ',F9.2, '[N]','   TORQUE  = ',F9.2, '[Nm]',
     & '  ADV     = ',F9.5, ' ',/ 
     & ' RIN    = ',F9.3, '[m]',         
     & '   ROUT    = ',F9.3, '[m]','   AREA    = ',F9.3, '[m²]'/,
     & ' VTIP   = ',F9.3, '[m/s]',' CFBLADE = ',F9.5, '  ',
     & '    CDBLADE = ',F9.5, '  ' /,
     & ' SIGMA  = ',F9.5, '  '/,/,
     & ' XCGI   = ',F9.3, '[m]','   YCGI    = ',F9.3, '[m]',
     & '   ZCGI    = ',F9.3, '[m]'/,
     & ' XCG    = ',F9.3, '[m]','   YCG     = ',F9.3, '[m]',
     & '   ZCG     = ',F9.3, '[m]'/,
     & ' ROTA1I = ',F9.3, '[deg]',' ROTB1I  = ',F9.3, '[deg]',
     & ' CONEAI  = ',F9.3, '[deg]'/,
     & ' ROTA1  = ',F9.3, '[deg]',' ROTB1   = ',F9.3, '[deg]',
     & ' CONEA   = ',F9.3, '[deg]'/,/,
     & ' IFA     = ',I2,'                IFT  = ',I2,
     & '                IFR  = ',I2,/
     & ' NSERIES = ',I2,/,/)

        WRITE(4,9661)IGR,IGRSLAVE(IGR),THRUST(IGR),
     &       TORQUE(IGR),ADV(IGR),RIN(IGR),ROUT(IGR),
     &       PII*(ROUT(IGR)**2-RIN(IGR)**2),
     &       VTIP(IGR),CFBLADE(IGR),CDBLADE(IGR),SIGMA(IGR),
     &       XCGI(IGR),YCGI(IGR),ZCGI(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),
     &       ROTA1I(IGR)*RAD2DEG,ROTB1I(IGR)*RAD2DEG,
     &       CONEAI(IGR)*RAD2DEG,ROTA1(IGR)*RAD2DEG,
     &       ROTB1(IGR)*RAD2DEG,CONEA(IGR)*RAD2DEG,
     &       IFA(IGR),IFT(IGR),IFR(IGR),NSERIES(IGR)

      ELSEIF (OSKU(IGR)%CHAR == 'P') THEN
 9662    FORMAT(/,' SHIP PROPEL NO.',I2
     & ',   CONNECTED TO SHIP NO. ',I2/,/,
     & ' CENAX  = ',F9.5, '[m]',' CENYX   = ',F9.5, '[m]',
     & ' CENAZ  = ',F9.5, '[m]',/ 
     & ' OMEGAX = ',F9.5, '[ ]',' OMEGAY  = ',F9.5, '[ ]',
     & ' OMEGAZ = ',F9.5, '[ ]',/ 
     & ' OMEGA  = ',F9.5, '[1/s]',/,/)

        WRITE(4,9662)IGR,IGRSLAVE(IGR),XCG(IGR),YCG(IGR),ZCG(IGR),
     &        OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,SHAFT(IGR)

      ELSEIF (OSKU(IGR)%CHAR == 'S') THEN
9670    FORMAT(/,' SHIP NO.',I2/,/,
     & ' MASS = ',F9.3, '[kg]','  AREF  = ',F9.3, '[m²]',
     & '   CHLREF  = ',F9.3, '[m]'/,         
     & ' IX   = ',F9.3, '[Nm]',
     & '              IY (trim direction) = ',F9.3, '[Nm]',/
     & ' IZ   = ',F9.3, '[Nm]','  IXY   = ',F9.3, '[Nm]',/
     & ' IXZ  = ',F9.3, '[Nm]','  IYZ   = ',F9.3, '[Nm]',/,/
     & ' XCG  = ',F9.3,'[m]','   YCG   = ',F9.3,'[m]',/
     & ' ZCG (sink direction)        = ',F9.3,'[m]',/
     & ' XCGI    = ',F9.3,'[m]',' YCGI = ',F9.3,'[m]',/
     & ' ZCGI (sink direction)       = ',F9.3,'[m]',/,/
     & ' TRIM    = ',F9.5,'[deg]',' TRIMAI   = ',F9.5,'[deg]',/
     & ' DRAUGHT = ',F9.3,'[m]','   DRAUGHTI = ',F9.3,'[m]',/
     & ' DAMPC1  = ',F9.3,'   ','   DAMPC2   = ',F9.3,''/,/)

        WRITE(4,9670)IGR,OSKU(IGR)%RMASS,OSKU(IGR)%AREF,
     &   OSKU(IGR)%CHLREF,
     &   OSKU(IGR)%IX,OSKU(IGR)%IY,OSKU(IGR)%IZ,
     &   OSKU(IGR)%IXY,OSKU(IGR)%IXZ,OSKU(IGR)%IYZ,
     &   XCG(IGR),YCG(IGR),ZCG(IGR),XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     &   TRIMA(IGR)*RAD2DEG,TRIMAI(IGR)*RAD2DEG,
     &   DRAUGHT(IGR),DRAUGHTI(IGR),DAMPC1(IGR),DAMPC2(IGR)

      ELSEIF (OSKU(IGR)%CHAR == 'R') THEN
9680    FORMAT(/,' BLADE NO.',I2/,/,
     & ' XCG    = ',F9.3,'[m]','   YCG   = ',F9.3,'[m]',
     & '   ZCG    = ',F9.3,'[m]',/
     & ' XCGI   = ',F9.3,'[m]','   YCGI  = ',F9.3,'[m]',
     & '   ZCGI   = ',F9.3,'[m]',/
     & ' PSIRI  = ',F9.3,'[deg]',' SHAFT = ',F9.3,'[deg]',
     & ' ROTORW = ',F9.3,'[rad/s]',/
     & ' REFVEL = ',F9.3,'[m/s]',' MASS  = ',F9.3,'[kg]',
     & '  RHINGE = ',F9.3,'[m]',/,/
     & ' RCG      = ',F9.3,'[m]','   RIZE     = ',F9.3,'[kgm²]',
     & ' RIBE     = ',F9.3,'[kgm²]',/
     & ' THETACOL = ',F9.3,'[deg]',' THCYCLON = ',F9.3,'[deg]',
     & '  THCYCLAT = ',F9.3,'[deg]',/
     & ' BETAHCOL = ',F9.3,'[deg]',' BECYCLON = ',F9.3,'[deg]',
     & '  BECYCLAT = ',F9.3,'[deg]',/
     & ' ZETAHCOL = ',F9.3,'[deg]',/,/
     & ' KDELTA3 = ',F9.3,'','      RKZE = ',F9.3,'[Nm/rad]',
     & ' RKBE = ',F9.3,'[Nm/rad]',/
     & ' RD      = ',F9.3,'[m]','   RKD  = ',F9.3,'[N/m]',
     & '    CD   = ',F9.3,'[N/(m/s)]',/
     & ' RCZE    = ',F9.3,'[Nm/(rad/s)]',
     &   '                   RCBE = ',F9.3,'[Nm/(rad/s)]',/,/)

        WRITE(4,9680)IGR,XCG(IGR),YCG(IGR),ZCG(IGR),
     &   XCGI(IGR),YCGI(IGR),ZCGI(IGR),PSIRI(IGR)*RAD2DEG,
     &   SHAFT(IGR)*RAD2DEG,ROTORW,REFVEL,OSKU(IGR)%RMASS,
     &   RHINGE(IGR),RCG(IGR),RIZE(IGR),RIBE(IGR),
     &   THETACOL(IGR)*RAD2DEG,THCYCLON(IGR)*RAD2DEG,
     &   THCYCLAT(IGR)*RAD2DEG,BETAHCOL(IGR)*RAD2DEG,
     &   BECYCLON(IGR)*RAD2DEG,BECYCLAT(IGR)*RAD2DEG,
     &   ZETAHCOL(IGR)*RAD2DEG,KDELTA3(IGR),RKZE(IGR),RKBE(IGR),
     &   RD(IGR),RKD(IGR),CD(IGR),RCZE(IGR),RCBE(IGR)


      ENDIF

      ENDDO

      RETURN
      END SUBROUTINE FLIGHT_REPORT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE READRS(TIME2,DREL)

      USE MPI

      USE TYPE_ARRAYS

      USE CONSTANTS,   ONLY : EPS

      USE NS3CO, ONLY : IOLD, NBLOCK, NBLOCG, IN, JN, KN, DT, DTOLD,
     &    TIMEL, COORL, PARALLEL, CHIMEL, REANEW, ITURB, IDIS, NSCAL,
     &    FRSPRE, TAVER, NAVER, MULPHL, FRESUL, TRANSL, !LONG_STARTL,
     &    LONG_TRANSL, AVER_STARTL, T, TWO_FLUIDL, START_CAVL

      USE MAIN_ARRAYS, ONLY : NLOCAL, NPNUM, IMAX, JMAX, KMAX,
     &    IMAXG, JMAXG, KMAXG, IG, NTOT, IGRID, U, V, W, RK, REPS, 
     &    DDEPS, FI, S11, VIST, EPS2, OHMI, PTUR, P, PDIFF, TEMP, 
     &    ROLE2, RMLE2, RNLE2, RWLE2, ELE2, RKLE2, EPSLE2, FILE2,
     &    ROLE3, RMLE3, RNLE3, RWLE3, ELE3, RKLE3, EPSLE3, FILE3,
     &    ROAV1, RMAV1, RNAV1, RWAV1, EAV1, RKAV1, EPSAV1, FIAV1,
     &    ROAV2, RMAV2, RNAV2, RWAV2, EAV2, RKAV2, EPSAV2, FIAV2,
     &    RUVAV, RUWAV, RVWAV, RMAV3, RNAV3,RWAV3,   PAV1,  TAV1,
     &    XCO, YCO, ZCO, XLE2, YLE2, ZLE2, XORI, YORI, ZORI,
     &    ROFOR, RMFOR, RNFOR, RWFOR, EFOR, PDFOR, RKSI, 
     &    RKFOR, REFOR, FIFOR, PRO, VAR, WAVEH, F1R, FTR2, TRM, PRC,
     &    BLKS,  F1RM,  F1RN,  F1RW, MPFOR

      USE INTEGERS, ONLY : MGM, IBF, MAXSS, MAXSB, MAXB, 
     &    MAXTI, MAXTS, IPRO


      IMPLICIT NONE

      REAL, ALLOCATABLE :: SBUF(:)
      REAL, ALLOCATABLE :: DBUF(:)

      INTEGER :: IERR, NCG, N, NROOT, NR1, IG1, IALKU, NTRFOR,
     +           NNMAX, I, NS, K, ISUM, J, IMEM, IERRCODE, IRK, MMAX, 
     +           ERRORCODE, IGN, ISUM2, STATUS(MPI_STATUS_SIZE),IPRESC,
     +           IPHASE

      REAL :: DREL

      LOGICAL :: TIME2, COOR2, CHIME2, MASTER

C     Read the solution part of the RSTART-file
C ... tama pitaisi periaatteessa kirjoittaa uudestaan. En viitsi. 30.12.99

      MASTER = IPRO == 1
      
      IF (IOLD > 0) THEN ! Read RSTART
***********************************************************************

      IF(MASTER) THEN
         
         READ(20) TIME2, COOR2, CHIME2
         IF(TIME2 .AND. .NOT.TIMEL) WRITE(*,*) 
     +' Odd! Starting steady calculation from time accurate RSTART.'
         IF(COOR2 .AND. .NOT.COORL) WRITE(*,*)'Odd! Starting ',
     +'nondeforming grid calculation from deforming grid RSTART.'
         IF(CHIME2 .AND. .NOT.CHIMEL) WRITE(*,*) 
     +' Odd! Starting nonchimera calculation from chimera RSTART.'
         IF(.NOT. CHIME2 .AND. CHIMEL) WRITE(*,*) 
     +        ' Chimeara calculation is starded using normal solution'
      ENDIF

      IF(PARALLEL) THEN
         CALL MPI_BCAST( TIME2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST( COOR2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CHIME2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      ENDIF

      DO 1040 NCG = 1,NBLOCG    ! Loop over the global blocks

      N     = NLOCAL(NCG)       ! local block number
      NROOT = NPNUM(NCG)        ! process number
      NR1   = NROOT - 1         ! process number used by MPI
      IPRESC = BLKS(NCG)%IPRESC

      IF(MASTER .AND. N /= -1) THEN ! Read master's blocks
 
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         READ(20) (   U(IG1+I),I=1,NNMAX)
         READ(20) (   V(IG1+I),I=1,NNMAX)
         READ(20) (   W(IG1+I),I=1,NNMAX)
         READ(20) (EPS2(IG1+I),I=1,NNMAX)
         READ(20) (VIST(IG1+I),I=1,NNMAX)
         READ(20) (OHMI(IG1+I),I=1,NNMAX)
         READ(20) (  RK(IG1+I),I=1,NNMAX)
         READ(20) (REPS(IG1+I),I=1,NNMAX)
         READ(20) (PTUR(IG1+I),I=1,NNMAX)
         READ(20)(DDEPS(IG1+I),I=1,NNMAX)
         READ(20)(PDIFF(IG1+I),I=1,NNMAX)
         READ(20)( TEMP(IG1+I),I=1,NNMAX)
         
         IF(IPRESC == 1 .OR. FRESUL) THEN
            READ(20) (PRC(IG1+I)%FIR,I=1,NNMAX)
            READ(20) (PRC(IG1+I)%FJR,I=1,NNMAX)
            READ(20) (PRC(IG1+I)%FKR,I=1,NNMAX)
         ELSE ! You have to read these although unused
            READ(20) (F1RM(IG1+I),I=1,NNMAX)
            READ(20) (F1RN(IG1+I),I=1,NNMAX)
            READ(20) (F1RW(IG1+I),I=1,NNMAX)
            WRITE(4,'(2X,2A)')'Mass flows',
     +      ' were read in READRS, but not used.'
            WRITE(13,'(2X,2A)')'Mass flows',
     +      ' were read in READRS, but not used.'
         ENDIF

         IF(TRANSL .AND. .NOT.LONG_TRANSL) THEN ! Intermittency variables
            READ(20) (TRM(IG1+I)%G,  I=1,NNMAX)
            READ(20) (TRM(IG1+I)%RET,I=1,NNMAX)
            WRITE(4,'(2X,2A,I4)')'Transition',
     +      ' variables were read in READRS for block',NCG
            WRITE(13,'(2X,2A,I4)')'Transition',
     +      ' variables were read in READRS for block',NCG
         ELSE IF (TRANSL .AND. LONG_TRANSL) THEN! These are initialized here
            DO I = 1,NNMAX
               TRM(IG1+I)%G   = BLKS(NCG)%FRSG
               TRM(IG1+I)%RET = BLKS(NCG)%FRSRET
            ENDDO
            WRITE(4,'(/,2X,2A)')'Initialization of transition',
     +      ' variables was made in READRS'
            WRITE(13,'(/,2X,2A)')'Initialization of transition',
     +      ' variables was made in READRS'
         ELSE IF (TRANSL) THEN
            WRITE(4,'(/,2X,2A)')'I am stymied of transition',
     +      '  initialization. This should be impossible in READRS?'
            WRITE(13,'(/,2X,2A)')'I am stymied of transition',
     +      '  initialization. This should be impossible in READRS?'
         ENDIF

C ... Scalars or Reynolds stresses
         IF(NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19))THEN
            DO NS = 1,NSCAL
               READ(20) (FI(IG1+I,NS),I=1,NNMAX)
            ENDDO
         ENDIF

         IF(.NOT.REANEW .AND. ITURB >= 10 .AND. ITURB <= 19) THEN
            DO NS = 1,6
               READ(20) (S11(IG1+I,NS),I=1,NNMAX)
            ENDDO
         ENDIF

         IF(MULPHL) THEN ! Multiphase variables
            DO NS = 1,NPHASES
               READ(20) (PRO(IG1+I)%TEMP(NS),I=1,NNMAX)
               READ(20) (VAR(IG1+I)%ALFA(NS),I=1,NNMAX)
               READ(20) (VAR(IG1+I)%X(NS),   I=1,NNMAX)
               READ(20) (PRO(IG1+I)%QIF(NS), I=1,NNMAX)
               IF(BLKS(NCG)%SOLUTION_TYPE == 'MULTI') THEN
               IF(.NOT. START_CAVL) THEN
               READ(20) (VAR(IG1+I)%U(NS),   I=1,NNMAX)
               READ(20) (VAR(IG1+I)%V(NS),   I=1,NNMAX)
               READ(20) (VAR(IG1+I)%W(NS),   I=1,NNMAX)
               ELSEIF(START_CAVL) THEN
               DO I = 1,NNMAX ! Single processor restart from CAVIT
               VAR(IG1+I)%U(NS) = U(IG1+I)
               VAR(IG1+I)%V(NS) = V(IG1+I)
               VAR(IG1+I)%W(NS) = W(IG1+I)
               ENDDO
             ENDIF ! START_CAVL
             ENDIF ! MULTI
             IF(BLKS(NCG)%SOLUTION_TYPE  == 'CAVIT'.AND.START_CAVL) THEN
               READ(20) (   U(IG1+I),I=1,NNMAX)
               READ(20) (   V(IG1+I),I=1,NNMAX)
               READ(20) (   W(IG1+I),I=1,NNMAX)
               ENDIF
            ENDDO

            DO NS = 1,NPHASES
               READ(20) (PRO(IG1+I)%DTEMP(NS),I=1,NNMAX)
c              Activate for old RSTART-files
c               READ(20) (PRO(IG1+I)%QIF(NS),  I=1,NNMAX) 
            ENDDO
         ENDIF

      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN ! Receive from the master

      IG1   = IG(1,N)
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IGN   = IG1 + NNMAX - 1
      K = NCG

      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
      
      CALL MPI_RECV(   U(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(   V(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(   W(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(EPS2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(VIST(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(OHMI(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(  RK(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(REPS(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(PTUR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      CALL MPI_RECV(DDEPS(IG1),NNMAX,MPI_REAL8,0,K,
     +     MPI_COMM_WORLD,STATUS,IERR)
      CALL MPI_RECV(PDIFF(IG1),NNMAX,MPI_REAL8,0,K,
     +     MPI_COMM_WORLD,STATUS,IERR)
      CALL MPI_RECV(TEMP(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)

!      IF(LONG_STARTL .AND. (IPRESC == 1 .OR. FRESUL)) THEN ! Antakee mersu
      CALL MPI_RECV(F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      IF(IPRESC == 1 .OR. FRESUL) PRC(IG1:IGN)%FIR = F1R(IG1:IGN)
      CALL MPI_RECV(F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      IF(IPRESC == 1 .OR. FRESUL) PRC(IG1:IGN)%FJR = F1R(IG1:IGN)
      CALL MPI_RECV(F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,STATUS,
     +     IERR)
      IF(IPRESC == 1 .OR. FRESUL) PRC(IG1:IGN)%FKR = F1R(IG1:IGN)
!      ENDIF

      IF(TRANSL .AND. .NOT.LONG_TRANSL) THEN
         CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)
         TRM(IG1:IGN)%G = F1R(IG1:IGN)
         CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)
         TRM(IG1:IGN)%RET = F1R(IG1:IGN)
      ENDIF

      IF(NSCAL > 0.AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
         DO NS = 1,NSCAL
            CALL MPI_RECV(FI(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF

      IF(.NOT.REANEW .AND. ITURB >= 10 .AND. ITURB <= 19) THEN
         DO NS = 1,6
            CALL MPI_RECV(S11(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF

      IF(MULPHL) THEN
         DO NS = 1,NPHASES
C TPS 12.9.08
c prosenttimerkit poistettu MPI_RECV kutsusta aputaulukon avulla
c      CALL MPI_RECV(PRO(IG1:IGN)%TEMP(NS),NNMAX,MPI_REAL8,0,K,
c     +     MPI_COMM_WORLD,STATUS,IERR)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            PRO(IG1:IGN)%TEMP(NS) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            VAR(IG1:IGN)%ALFA(NS) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            VAR(IG1:IGN)%X(NS) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            PRO(IG1:IGN)%QIF(NS) = F1R(IG1:IGN)
            IF(TWO_FLUIDL) THEN
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            VAR(IG1:IGN)%U(NS) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            VAR(IG1:IGN)%V(NS) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            VAR(IG1:IGN)%W(NS) = F1R(IG1:IGN)
            ENDIF

         ENDDO

         ALLOCATE(DBUF(NNMAX))

         DO NS = 1,NPHASES
            CALL MPI_RECV(DBUF(1:NNMAX),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            PRO(IG1:IGN)%DTEMP(NS) = DBUF(1:NNMAX)
         ENDDO

         DEALLOCATE(DBUF)

      ENDIF


      ELSEIF(MASTER .AND. N == -1) THEN ! Read and send to the slaves
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)
         ALLOCATE(SBUF(NNMAX))

         ISUM = 0
         IF(.NOT.REANEW .AND.ITURB >= 10 .AND. ITURB <= 19) ISUM = 6
         IF(NSCAL > 0.AND. (ITURB >= 23 .OR. ITURB <= 19))
     +        ISUM = ISUM + NSCAL
         IF(TRANSL .AND. .NOT.LONG_TRANSL) ISUM = ISUM + 2
         IF(MULPHL) ISUM = ISUM + 4*NPHASES
         IF(TWO_FLUIDL) ISUM = ISUM + 3*NPHASES

!         IALKU = 12
         IALKU = 15
!     IF(LONG_STARTL .AND. (IPRESC == 1 .OR. FRESUL)) IALKU = 15         
       
         DO I = 1,IALKU+ISUM
            READ(20) (SBUF(J),J=1,NNMAX)
            CALL MPI_SSEND(SBUF,NNMAX,MPI_REAL8,NR1,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO

         DEALLOCATE(SBUF)

         IF(MULPHL) THEN
         ALLOCATE(DBUF(NNMAX))      
         DO I = 1,NPHASES
            READ(20) (DBUF(J),J=1,NNMAX)
            CALL MPI_SSEND(DBUF,NNMAX,MPI_REAL8,NR1,K,
     +           MPI_COMM_WORLD,IERR)
         ENDDO

         DEALLOCATE(DBUF)

         ENDIF ! MULPHL
      ENDIF                     ! (MASTER .AND. N /= -1)

      IF(N /= -1) THEN ! Evaluate the local pressure (Miksköhän??)
         DO 10 I = 1,NNMAX
            P(IG1+I) = PDIFF(IG1+I) + FRSPRE
 10      CONTINUE
      ENDIF

 1040 CONTINUE

************************************************************************

      IF (TRANSL .AND. LONG_TRANSL) THEN  ! These are initialized here
         DO N=1,NBLOCK
            IG1   = IG(1,N)
            NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
            IGN   = IG1 + NNMAX - 1
            DO I = IG1,IGN
               TRM(I)%G   = BLKS(N)%FRSG
               TRM(I)%RET = BLKS(N)%FRSRET
            ENDDO
         ENDDO
         WRITE(4,'(/,2X,2A)')'Initialization of transition',
     +        ' variables was made in READRS'
         WRITE(13,'(/,2X,2A)')'Initialization of transition',
     +        ' variables was made in READRS'
      ELSE
         WRITE(4,'(2X,2A)')'Transition variables were not initialized'
     +        ,'  in READRS.'
        WRITE(13,'(/,2X,2A)')'Transition variables were not initialized'
     +        ,'  in READRS. I assume you know what you are doing. :-}'
      ENDIF

************************************************************************

      IF(TIMEL .AND. TIME2) THEN ! Continue time-accurate simulation
      
      DO 1050 NCG = 1,NBLOCG     ! Loop over the global blocks

      N     = NLOCAL(NCG)        ! local block number
      NROOT = NPNUM(NCG)         ! process number
      NR1   = NROOT - 1          ! process number used by MPI
      
      IF(MASTER .AND. N /= -1) THEN ! Read master's own blocks

         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         IGN   = IG1 + NNMAX - 1

         READ(20)( ROLE2(IG1+I),I=1,NNMAX)
         READ(20)( RMLE2(IG1+I),I=1,NNMAX)
         READ(20)( RNLE2(IG1+I),I=1,NNMAX)
         READ(20)( RWLE2(IG1+I),I=1,NNMAX)
         READ(20)( RKLE2(IG1+I),I=1,NNMAX)
         READ(20)(EPSLE2(IG1+I),I=1,NNMAX)
         READ(20)(  ELE2(IG1+I),I=1,NNMAX)

         IF(MULPHL) THEN
            DO NS = 1,NPHASES
               READ(20)(F1R(IG1+I),I=1,NNMAX)
               VAR(IG1+1:IGN+1)%AROLE2(NS) = F1R(IG1+1:IGN+1)
               READ(20)(F1R(IG1+I),I=1,NNMAX)
               VAR(IG1+1:IGN+1)%ARELE2(NS) = F1R(IG1+1:IGN+1)
            ENDDO
         ENDIF ! MULPHL

         IF(NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
            READ(20) ((FILE2(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF

         READ(20)( ROAV1(IG1+I),I=1,NNMAX)
         READ(20)( RMAV1(IG1+I),I=1,NNMAX)
         READ(20)( RNAV1(IG1+I),I=1,NNMAX)
         READ(20)( RWAV1(IG1+I),I=1,NNMAX)
         READ(20)( RKAV1(IG1+I),I=1,NNMAX)
         READ(20)(EPSAV1(IG1+I),I=1,NNMAX)
         READ(20)(  EAV1(IG1+I),I=1,NNMAX)

         IF(NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
            READ(20) ((FIAV1(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF

         READ(20)( ROAV2(IG1+I),I=1,NNMAX)
         READ(20)( RMAV2(IG1+I),I=1,NNMAX)
         READ(20)( RNAV2(IG1+I),I=1,NNMAX)
         READ(20)( RWAV2(IG1+I),I=1,NNMAX)
         READ(20)(  PAV1(IG1+I),I=1,NNMAX)  ! Here, antakee mersu
         READ(20)(  TAV1(IG1+I),I=1,NNMAX)
*         READ(20)( RKAV2(IG1+I),I=1,NNMAX) ! Unused k and epsilon
*         READ(20)(EPSAV2(IG1+I),I=1,NNMAX)
         READ(20)(  EAV2(IG1+I),I=1,NNMAX)
         READ(20)( RUVAV(IG1+I),I=1,NNMAX)
         READ(20)( RUWAV(IG1+I),I=1,NNMAX)
         READ(20)( RVWAV(IG1+I),I=1,NNMAX)

         IF(NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19))THEN
            READ(20) ((FIAV2(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
      ENDIF                     ! (MASTER .AND. N /= -1)
      
      IF(.NOT.MASTER .AND. N /= -1) THEN ! Receive from the master

      IG1   = IG(1,N)
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IGN   = IG1 + NNMAX - 1
      K = NCG

      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
      
      CALL MPI_RECV( ROLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RMLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RNLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RWLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RKLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV(EPSLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV(  ELE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)

      IF(MULPHL) THEN
         DO NS = 1,NPHASES
         CALL MPI_RECV(F1R(IG1),NNMAX,MPI_REAL8,0,K,
     +                          MPI_COMM_WORLD,STATUS,IERR)
         VAR(IG1:IGN)%AROLE2(NS) = F1R(IG1:IGN)
         CALL MPI_RECV(F1R(IG1),NNMAX,MPI_REAL8,0,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)                       
         VAR(IG1:IGN)%ARELE2(NS) = F1R(IG1:IGN)
         ENDDO
      ENDIF ! MULPHL

      IF (NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
         DO NS = 1,NSCAL
            CALL MPI_RECV(FILE2(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF

      CALL MPI_RECV( ROAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RMAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RNAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RWAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RKAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV(EPSAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV(  EAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      IF (NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
         DO NS = 1,NSCAL
            CALL MPI_RECV(FIAV1(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF

      CALL MPI_RECV( ROAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RMAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RNAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV( RWAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV(  PAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)                                           
      CALL MPI_RECV(  TAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV(  EAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV( RUVAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV( RUWAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)
      CALL MPI_RECV( RVWAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +     STATUS,IERR)

      IF (NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) THEN
         DO NS = 1,NSCAL
            CALL MPI_RECV(FIAV2(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF

      ELSEIF(MASTER .AND. N == -1) THEN  ! Read and send to slaves

         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,MPI_COMM_WORLD,STATUS,
     +        IERR)

         ALLOCATE(SBUF(NNMAX))
         ISUM = 0

         IF(NSCAL > 0 .AND. (ITURB >= 23 .OR. ITURB <= 19)) ISUM = NSCAL

         ISUM2 = 0
         IF(MULPHL) ISUM2 = 2*NPHASES
     
         DO I = 1,3*(7+ISUM)+3+ISUM2  ! Mersu
            READ(20) (SBUF(J),J=1,NNMAX)
            CALL MPI_SSEND(SBUF,NNMAX,MPI_REAL8,NR1,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO

         DEALLOCATE(SBUF)

      ENDIF                     ! (MASTER .AND. N /= -1)

      DREL = DT/DTOLD ! time step change is done in INPOST

1050  CONTINUE

      ENDIF                     ! TIMEL

C ... End time accurate data
***************************************************************************

      IF(COORL .AND. COOR2) THEN ! Coordinates are changing

      DO 1060 NCG = 1,NBLOCG         ! Loop over the global blocks

      N     = NLOCAL(NCG)       ! local block number
      NROOT = NPNUM(NCG)        ! process number
      NR1   = NROOT - 1         ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN  ! Master's own blocks
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         READ(20)( XLE2(IG1+I),I=1,NNMAX)
         READ(20)( YLE2(IG1+I),I=1,NNMAX)
         READ(20)( ZLE2(IG1+I),I=1,NNMAX)
         READ(20)( XORI(IG1+I),I=1,NNMAX)
         READ(20)( YORI(IG1+I),I=1,NNMAX)
         READ(20)( ZORI(IG1+I),I=1,NNMAX)
      ENDIF                     ! (MASTER .AND. N /= -1)
      
      IF(.NOT.MASTER .AND. N /= -1) THEN  ! Slaves receive from the master
         IG1   = IG(1,N)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         K = NCG

         CALL MPI_SEND(NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)

         CALL MPI_RECV(XLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(YLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(ZLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(XORI(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(YORI(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(ZORI(IG1),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)

      ELSEIF(MASTER .AND. N == -1) THEN ! Master reads and sends to the slaves
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,MPI_COMM_WORLD,STATUS,
     +        IERR)

         ALLOCATE(DBUF(NNMAX))

         DO I = 1,6
            READ(20) (DBUF(J),J=1,NNMAX)
            CALL MPI_SSEND(DBUF,NNMAX,MPI_REAL8,NR1,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO

         DEALLOCATE(DBUF)

      ENDIF                     ! (MASTER .AND. N /= -1)
1060  CONTINUE                  ! NCG = 1,NBLOCG   
      ENDIF                     ! COORL


      IF(COORL .AND. .NOT.COOR2) THEN  ! Coordinates are start to change
         
      DO 1065 NCG = 1,NBLOCG           ! Loop over the global blocks

      N     = NLOCAL(NCG)              ! local block number
      NROOT = NPNUM(NCG)               ! process number
      NR1   = NROOT - 1                ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN    ! Master's own blocks

         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)

         DO I = 1,NNMAX
            XLE2(IG1+I) = XCO(IG1+I)
            YLE2(IG1+I) = YCO(IG1+I)
            ZLE2(IG1+I) = ZCO(IG1+I)
            XORI(IG1+I) = XCO(IG1+I)
            YORI(IG1+I) = YCO(IG1+I)
            ZORI(IG1+I) = ZCO(IG1+I)
         ENDDO

      ENDIF                     ! (MASTER .AND. N /= -1)
      
      IF(.NOT.MASTER .AND. N /= -1) THEN  ! Slaves

         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         K = NCG

         DO I = 1,NNMAX
            XLE2(IG1+I) = XCO(IG1+I)
            YLE2(IG1+I) = YCO(IG1+I)
            ZLE2(IG1+I) = ZCO(IG1+I)
            XORI(IG1+I) = XCO(IG1+I)
            YORI(IG1+I) = YCO(IG1+I)
            ZORI(IG1+I) = ZCO(IG1+I)
         ENDDO

      ENDIF                     ! (MASTER .AND. N /= -1)

1065  CONTINUE                  ! NCG = 1,NBLOCG   

      ENDIF                     ! COORL started from .NOT.COORL

C ... End coordinates changed
****************************************************************************

*      IF(CHIMEL .AND. CHIME2 .AND. COORL) THEN ! Chimera blocks
      IF(CHIMEL .AND. CHIME2) THEN   ! Chimera blocks
      
      DO 1070 NCG = 1,NBLOCG         ! Loop over the global blocks

      N     = NLOCAL(NCG)            ! local block number
      NROOT = NPNUM(NCG)             ! process number
      NR1   = NROOT - 1              ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN  ! Master's own blocks

         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)

         READ(20)( ROFOR(IG1+I),I=1,NNMAX)
         READ(20)( RMFOR(IG1+I),I=1,NNMAX)
         READ(20)( RNFOR(IG1+I),I=1,NNMAX)
         READ(20)( RWFOR(IG1+I),I=1,NNMAX)
         READ(20)(  EFOR(IG1+I),I=1,NNMAX)
         READ(20)( PDFOR(IG1+I),I=1,NNMAX)
         READ(20)( RKFOR(IG1+I),I=1,NNMAX)
         READ(20)( REFOR(IG1+I),I=1,NNMAX)
         READ(20)(  RKSI(IG1+I),I=1,NNMAX)

         IF(TRANSL .AND. .NOT.LONG_TRANSL) THEN
            READ(20) (TRM(IG1+I)%GFOR,  I=1,NNMAX)
            READ(20) (TRM(IG1+I)%RETFOR,I=1,NNMAX)
         ENDIF

         IF(MULPHL) THEN
            DO IPHASE = 1,NPHASES
               READ(20) (MPFOR(IG1+I)%ALFAFOR(IPHASE),  I=1,NNMAX)
               READ(20) (MPFOR(IG1+I)%XFOR(IPHASE),     I=1,NNMAX)
               READ(20) (MPFOR(IG1+I)%DTEMPFOR(IPHASE), I=1,NNMAX)
            ENDDO
         ENDIF

         IF (NSCAL > 0) THEN
             DO NS = 1,NSCAL
                READ(20) (FIFOR(IG1+I,NS),I=1,NNMAX)
             ENDDO
         ENDIF
	 
      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN  ! Slaves receive
      
      IG1   = IG(1,N)      
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      K = NCG

      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)

      CALL MPI_RECV( ROFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( RMFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( RNFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( RWFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV(  EFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( PDFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( RKFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV( REFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
      CALL MPI_RECV(  RKSI(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +        STATUS,IERR)

      IF(TRANSL .AND. .NOT.LONG_TRANSL) THEN
         CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         TRM(IG1:IGN)%GFOR = F1R(IG1:IGN)
         CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,STATUS,IERR)
         TRM(IG1:IGN)%RETFOR = F1R(IG1:IGN)
      ENDIF

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,STATUS,IERR)
                 MPFOR(IG1:IGN)%ALFAFOR(IPHASE) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,STATUS,IERR)
                 MPFOR(IG1:IGN)%XFOR(IPHASE) = F1R(IG1:IGN)
            CALL MPI_RECV(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,STATUS,IERR)
                 MPFOR(IG1:IGN)%DTEMPFOR(IPHASE) = F1R(IG1:IGN)
         ENDDO
      ENDIF

      IF (NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_RECV(FIFOR(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDDO
      ENDIF
      
      ELSEIF(MASTER .AND. N == -1) THEN  ! Master reads and sends
      
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)

         ALLOCATE(SBUF(NNMAX))

         NTRFOR = 0
         IF(TRANSL .AND. .NOT.LONG_TRANSL) NTRFOR = 2      

         DO I = 1,9+NTRFOR+NSCAL
            READ(20) (SBUF(J),J=1,NNMAX)
            CALL MPI_SSEND(SBUF,NNMAX,MPI_REAL8,NR1,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO

         DEALLOCATE(SBUF)

      ENDIF                     ! (MASTER .AND. N /= -1)

1070  CONTINUE                  ! NCG = 1,NBLOCG   

      ENDIF                     ! CHIMEL
      
C ... End Chimera blocks

**********************************************************************      
      ENDIF                     ! IOLD > 0
**********************************************************************      

C ... Do something (probably usefull) anyway

      IF(TIMEL) THEN

         IF(AVER_STARTL) THEN ! Initialize averaging

           WRITE(4,'(/2A,E12.5)') ' Averaging was started from',
     &     ' TIME =',T

           TAVER = 0.
           NAVER = 0.

         ELSE ! Continue averaging

           WRITE(4,'(/2A,E12.5,A,E12.5)')' Averaging is continued from',
     &     ' time =',T-TAVER,'. Current time =',T

         ENDIF

         DO I = 1,MAXTI
           RMAV3(I) = RMAV2(I) - RMAV1(I)**2/(ROAV1(I)+EPS)
           RNAV3(I) = RNAV2(I) - RNAV1(I)**2/(ROAV1(I)+EPS)
           RWAV3(I) = RWAV2(I) - RWAV1(I)**2/(ROAV1(I)+EPS)
         ENDDO

      ENDIF ! TIMEL

      RETURN
      END SUBROUTINE READRS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITRO(IOLD,NBLOCK,IMAX,JMAX,KMAX,IG,MGM,IN,JN,KN,
     +  NTOT,TIMEL,COORL,CHIMEL,REANEW,ITURB,NSCAL,MAXSS,MAXSB,
     +  RO,RM,RN,RW,E,RK,REPS,DDEPS,FI,S11,VIST,EPS2,OHMI,PTUR,P,PDIFF,
     +  ROLE2,RMLE2,RNLE2,RWLE2,ELE2,RKLE2,EPSLE2,FILE2,
     +  XCO,YCO,ZCO,XLE2,YLE2,ZLE2,
     +  ROFOR,RMFOR,RNFOR,RWFOR,EFOR,PDFOR,RKSI,RKFOR,REFOR,FIFOR)

      DIMENSION :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),IG(MGM,*)

      REAL :: RO(*),RM(*),RN(*),RW(*),E(*),RK(*),REPS(*),FI(MAXSB,*),
     +        ROLE2(*),RMLE2(*),RNLE2(*),RWLE2(*),
     +        ELE2(*),RKLE2(*),EPSLE2(*),FILE2(MAXSB,*),S11(MAXSS,*),
     +        VIST(*),EPS2(*),PTUR(*),P(*),
     +        PDIFF(*),DDEPS(*),OHMI(*),
     +        ROFOR(*),RMFOR(*),RNFOR(*),RWFOR(*),EFOR(*),PDFOR(*),
     +        RKSI(*),RKFOR(*),REFOR(*),FIFOR(MAXSB,*)

      REAL :: XCO(*),YCO(*),ZCO(*),XLE2(*),YLE2(*),ZLE2(*)

      LOGICAL :: TIMEL,COORL,REANEW,CHIMEL


      WRITE(20) TIMEL,COORL,CHIMEL
      DO 1040 N=1,NBLOCK
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20) (RO(IG1+I),I=1,NNMAX)
         WRITE(20) (  RM(IG1+I),I=1,NNMAX)
         WRITE(20) (  RN(IG1+I),I=1,NNMAX)
         WRITE(20) (  RW(IG1+I),I=1,NNMAX)
         WRITE(20) (   E(IG1+I),I=1,NNMAX)
         WRITE(20) (EPS2(IG1+I),I=1,NNMAX)
         WRITE(20) (VIST(IG1+I),I=1,NNMAX)
         WRITE(20) (OHMI(IG1+I),I=1,NNMAX)
         WRITE(20) (  RK(IG1+I),I=1,NNMAX)
         WRITE(20) (REPS(IG1+I),I=1,NNMAX)
         WRITE(20) (PTUR(IG1+I),I=1,NNMAX)
         WRITE(20)(DDEPS(IG1+I),I=1,NNMAX)
         WRITE(20)(PDIFF(IG1+I),I=1,NNMAX)
         IF (NSCAL > 0) THEN
            WRITE(20) ((FI(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
         IF (ITURB >= 10 .AND. ITURB <= 19) THEN
            WRITE(20) ((S11(IG1+I,NS),I=1,NNMAX),NS=1,6)
         ENDIF
 1040 CONTINUE

      IF(TIMEL) THEN
      DO 1050 N=1,NBLOCK
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20)( ROLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RMLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RNLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RWLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RKLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(EPSLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(  ELE2(IG1+I),I=1,NNMAX)
         IF (NSCAL > 0) THEN
            WRITE(20) ((FILE2(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
 1050 CONTINUE
      ENDIF           ! TIMEL
      IF(COORL) THEN
      DO 1060 N=1,NBLOCK
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20)( XLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( YLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( ZLE2(IG1+I),I=1,NNMAX)
 1060 CONTINUE
      ENDIF           ! COORL
      IF(CHIMEL) THEN
      DO 1070 N=1,NBLOCK
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20)( ROFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RMFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RNFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RWFOR(IG1+I),I=1,NNMAX)
         WRITE(20)(  EFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( PDFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RKFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( REFOR(IG1+I),I=1,NNMAX)
         WRITE(20)(  RKSI(IG1+I),I=1,NNMAX)
         IF (NSCAL > 0) THEN
            WRITE(20) ((FIFOR(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
 1070 CONTINUE
      ENDIF           ! CHIMEL

      RETURN
      END
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITRS

      USE MPI
      USE TYPE_ARRAYS

      USE NS3CO, ONLY : IOLD, NBLOCK, NBLOCG, IN, JN, KN, TIMEL, 
     &    COORL,CHIMEL, PARALLEL, REANEW, ITURB, NSCAL, MULPHL, 
     &    FRESUL, TRANSL,TWO_FLUIDL

      USE MAIN_ARRAYS, ONLY : NLOCAL, NPNUM, IMAX, JMAX, KMAX,
     &    IMAXG, JMAXG, KMAXG, IG, NTOT, U, V, W, RK, REPS, 
     &    DDEPS, FI, S11, VIST, EPS2, OHMI, PTUR, P, PDIFF, TEMP,
     &    ROLE2, RMLE2, RNLE2, RWLE2, ELE2, RKLE2, EPSLE2, FILE2,
     &    ROAV1, RMAV1, RNAV1, RWAV1, EAV1, RKAV1, EPSAV1, FIAV1,
     &    ROAV2, RMAV2, RNAV2, RWAV2, EAV2, RKAV2, EPSAV2, FIAV2,
     &    RUVAV, RUWAV, RVWAV,  PAV1, TAV1,
     &    XCO, YCO, ZCO, XLE2, YLE2, ZLE2, XORI, YORI, ZORI,
     &    ROFOR, RMFOR, RNFOR, RWFOR, EFOR, PDFOR, RKSI, 
     &    RKFOR, REFOR, FIFOR, PRO,VAR,WAVEH, F1R, FTR2, TRM, PRC,
     &    BLKS, MPFOR

      USE INTEGERS, ONLY : IPRO, MGM, IBF, MAXSS, MAXSB, MAXB 

      IMPLICIT NONE

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      INTEGER :: NCG, N, NROOT, NR1, IG1, NNMAX, I, NS, K, IERR, NTRFOR,
     &           ISUM, J, IMEM, IERRCODE,ERRORCODE, IGN, IPRESC, IALKU,
     &           IPHASE

      REAL, ALLOCATABLE :: SBUF(:)
      REAL, ALLOCATABLE :: DBUF(:)

      LOGICAL :: MASTER

      MASTER = IPRO == 1

      IF(MASTER) WRITE(20) TIMEL, COORL, CHIMEL

      DO 1040 NCG = 1,NBLOCG    ! Loop over the global blocks

      N     = NLOCAL(NCG)       ! local block number
      NROOT = NPNUM(NCG)        ! process number
      NR1   = NROOT - 1         ! process number used by MPI
      IPRESC = BLKS(NCG)%IPRESC

      IF(MASTER .AND. N /= -1) THEN  ! Write master's blocks
      
         IG1   = IG(1,N)-IG(1,1)

         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)

         WRITE(20) (   U(IG1+I),I=1,NNMAX)
         WRITE(20) (   V(IG1+I),I=1,NNMAX)
         WRITE(20) (   W(IG1+I),I=1,NNMAX)
         WRITE(20) (EPS2(IG1+I),I=1,NNMAX)
         WRITE(20) (VIST(IG1+I),I=1,NNMAX)
         WRITE(20) (OHMI(IG1+I),I=1,NNMAX)
         WRITE(20) (  RK(IG1+I),I=1,NNMAX)
         WRITE(20) (REPS(IG1+I),I=1,NNMAX)
         WRITE(20) (PTUR(IG1+I),I=1,NNMAX)
         WRITE(20)(DDEPS(IG1+I),I=1,NNMAX)
         WRITE(20)(PDIFF(IG1+I),I=1,NNMAX)
         WRITE(20) (TEMP(IG1+I),I=1,NNMAX)
         
         IF(IPRESC == 1 .OR. FRESUL) THEN
            WRITE(20) (PRC(IG1+I)%FIR,I=1,NNMAX)
            WRITE(20) (PRC(IG1+I)%FJR,I=1,NNMAX)
            WRITE(20) (PRC(IG1+I)%FKR,I=1,NNMAX)
         ELSE ! PRC is not allocated, let us use these
            WRITE(20)(DDEPS(IG1+I),I=1,NNMAX)
            WRITE(20)(PDIFF(IG1+I),I=1,NNMAX)
            WRITE(20) (TEMP(IG1+I),I=1,NNMAX)
         ENDIF


         IF(TRANSL) THEN  ! Intermittency variables
            WRITE(20) (TRM(IG1+I)%G,   I=1,NNMAX)
            WRITE(20) (TRM(IG1+I)%RET, I=1,NNMAX)
         ENDIF

         IF(NSCAL > 0) THEN  ! Scalars
            DO NS = 1,NSCAL
               WRITE(20) (FI(IG1+I,NS),I=1,NNMAX)
            ENDDO
         ENDIF

         IF (ITURB >= 10 .AND. ITURB <= 19) THEN ! Reynolds stresses
            DO NS = 1,6
               WRITE(20) (S11(IG1+I,NS),I=1,NNMAX)
            ENDDO
         ENDIF
	 
         IF (MULPHL) THEN  ! Multiphase variables
            DO NS = 1,NPHASES
               WRITE(20) (PRO(IG1+I)%TEMP(NS), I=1,NNMAX)
               WRITE(20) (VAR(IG1+I)%ALFA(NS), I=1,NNMAX)
               WRITE(20) (VAR(IG1+I)%X(NS),    I=1,NNMAX)
               WRITE(20) (PRO(IG1+I)%QIF(NS),  I=1,NNMAX)
               IF(BLKS(NCG)%SOLUTION_TYPE == 'MULTI') THEN
               WRITE(20) (VAR(IG1+I)%U(NS),   I=1,NNMAX)
               WRITE(20) (VAR(IG1+I)%V(NS),   I=1,NNMAX)
               WRITE(20) (VAR(IG1+I)%W(NS),   I=1,NNMAX)
               ENDIF

            ENDDO

            DO NS = 1,NPHASES
               WRITE(20) (PRO(IG1+I)%DTEMP(NS),I=1,NNMAX)
c              Activate to write RSTART for old versions               
c              WRITE(20) (PRO(IG1+I)%QIF(NS),  I=1,NNMAX)
            ENDDO
         ENDIF

      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN ! Send to the master
      
      IG1   = IG(1,N)
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IGN   = IG1 + NNMAX - 1

      K = NCG
      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
      
      CALL MPI_SSEND(    U(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(    V(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(    W(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND (EPS2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( VIST(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( OHMI(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(   RK(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( REPS(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( PTUR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(DDEPS(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(PDIFF(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( TEMP(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)

!      IF(IPRESC == 1 .OR. FRESUL) THEN
      IF(IPRESC == 1 .OR. FRESUL) F1R(IG1:IGN) = PRC(IG1:IGN)%FIR
      CALL MPI_SSEND( F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      IF(IPRESC == 1 .OR. FRESUL) F1R(IG1:IGN) = PRC(IG1:IGN)%FJR
      CALL MPI_SSEND( F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      IF(IPRESC == 1 .OR. FRESUL) F1R(IG1:IGN) = PRC(IG1:IGN)%FKR
      CALL MPI_SSEND( F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
!      ENDIF

      IF(TRANSL) THEN
         F1R(IG1:IGN) = TRM(IG1:IGN)%G
         CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         F1R(IG1:IGN) = TRM(IG1:IGN)%RET
         CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
      ENDIF

      IF(NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_SSEND(FI(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF

      IF(ITURB >= 10 .AND. ITURB <= 19) THEN
         DO NS = 1,6
            CALL MPI_SSEND(S11(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF

      IF (MULPHL) THEN
         DO NS = 1,NPHASES
            F1R(IG1:IGN) = PRO(IG1:IGN)%TEMP(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = VAR(IG1:IGN)%ALFA(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = VAR(IG1:IGN)%X(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = PRO(IG1:IGN)%QIF(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            IF(TWO_FLUIDL) THEN
            F1R(IG1:IGN) = VAR(IG1:IGN)%U(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = VAR(IG1:IGN)%V(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = VAR(IG1:IGN)%W(NS)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                    MPI_COMM_WORLD,IERR)
            ENDIF

         ENDDO

      DO NS = 1,NPHASES
         FTR2(1:NNMAX) = PRO(IG1:IGN)%DTEMP(NS)
         CALL MPI_SSEND(FTR2(1:NNMAX),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
      ENDDO

      ENDIF

      
      ELSEIF(MASTER .AND. N == -1) THEN ! Receive from the slaves and write

         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,MPI_COMM_WORLD,STATUS,
     +        IERR)
         ALLOCATE(SBUF(NNMAX))
         ISUM = 0
         IF(ITURB >= 10 .AND. ITURB <= 19) ISUM = 6
         IF(TRANSL) ISUM = ISUM + 2
         IF(MULPHL) ISUM = ISUM + 4*NPHASES
         IF(TWO_FLUIDL) ISUM = ISUM + 3*NPHASES
!         IALKU = 12
         IALKU = 15
!         IF(IPRESC == 1 .OR. FRESUL) IALKU = 15
         DO I = 1,IALKU+NSCAL+ISUM ! Long list always
            CALL MPI_RECV(SBUF,NNMAX,MPI_REAL8,NR1,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            WRITE(20) (SBUF(J),J=1,NNMAX)
         ENDDO
         DEALLOCATE(SBUF)
           
         IF(MULPHL) THEN
            ALLOCATE(DBUF(NNMAX))
            DO I = 1,NPHASES
               CALL MPI_RECV(DBUF,NNMAX,MPI_REAL8,NR1,K,
     +                       MPI_COMM_WORLD,STATUS,IERR)
               WRITE(20) (DBUF(J),J=1,NNMAX)
            ENDDO
         DEALLOCATE(DBUF)
         ENDIF  ! MULPHL

      ENDIF  ! (MASTER .AND. N /= -1)

 1040 CONTINUE

*********************************************************************

      IF(TIMEL) THEN                 ! Store time-dependent data
      
      DO NCG = 1,NBLOCG              ! Loop over the global blocks
  
      N     = NLOCAL(NCG)            ! local block number
      NROOT = NPNUM(NCG)             ! process number
      NR1   = NROOT - 1              ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN  ! Write master's own blocks
      
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         IGN   = IG1 + NNMAX - 1

         WRITE(20)( ROLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RMLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RNLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RWLE2(IG1+I),I=1,NNMAX)
         WRITE(20)( RKLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(EPSLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(  ELE2(IG1+I),I=1,NNMAX)

         IF(MULPHL) THEN
            DO NS = 1,NPHASES
               F1R(IG1+1:IGN+1) = VAR(IG1+1:IGN+1)%AROLE2(NS)
               WRITE(20)(F1R(IG1+I),I=1,NNMAX)
               F1R(IG1+1:IGN+1) = VAR(IG1+1:IGN+1)%ARELE2(NS)
               WRITE(20)(F1R(IG1+I),I=1,NNMAX)
            ENDDO
         ENDIF  ! MULPHL

         IF (NSCAL > 0) THEN
            WRITE(20) ((FILE2(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF

         WRITE(20)( ROAV1(IG1+I),I=1,NNMAX)
         WRITE(20)( RMAV1(IG1+I),I=1,NNMAX)
         WRITE(20)( RNAV1(IG1+I),I=1,NNMAX)
         WRITE(20)( RWAV1(IG1+I),I=1,NNMAX)
         WRITE(20)( RKAV1(IG1+I),I=1,NNMAX)
         WRITE(20)(EPSAV1(IG1+I),I=1,NNMAX)
         WRITE(20)(  EAV1(IG1+I),I=1,NNMAX)
         IF(NSCAL > 0) THEN
            WRITE(20) ((FIAV1(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
         WRITE(20)( ROAV2(IG1+I),I=1,NNMAX)
         WRITE(20)( RMAV2(IG1+I),I=1,NNMAX)
         WRITE(20)( RNAV2(IG1+I),I=1,NNMAX)
         WRITE(20)( RWAV2(IG1+I),I=1,NNMAX)
         WRITE(20)(  PAV1(IG1+I),I=1,NNMAX)
         WRITE(20)(  TAV1(IG1+I),I=1,NNMAX)
*         WRITE(20)( RKAV2(IG1+I),I=1,NNMAX)
*         WRITE(20)(EPSAV2(IG1+I),I=1,NNMAX)
         WRITE(20)(  EAV2(IG1+I),I=1,NNMAX)
         WRITE(20)( RUVAV(IG1+I),I=1,NNMAX)
         WRITE(20)( RUWAV(IG1+I),I=1,NNMAX)
         WRITE(20)( RVWAV(IG1+I),I=1,NNMAX)
         IF (NSCAL > 0) THEN
            WRITE(20) ((FIAV2(IG1+I,NS),I=1,NNMAX),NS=1,NSCAL)
         ENDIF
	 
      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN  ! Send to the master
      
      IG1   = IG(1,N)
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      IGN   = IG1 + NNMAX - 1
      K = NCG

      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
      
      CALL MPI_SSEND( ROLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND( RMLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND( RNLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND( RWLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND( RKLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND(EPSLE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)
      CALL MPI_SSEND(  ELE2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                IERR)

      IF(MULPHL) THEN
         DO NS = 1,NPHASES
         F1R(IG1:IGN) = VAR(IG1:IGN)%AROLE2(NS)
         CALL MPI_SSEND(F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                  IERR)
         F1R(IG1:IGN) = VAR(IG1:IGN)%ARELE2(NS)
         CALL MPI_SSEND(F1R(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +                  IERR)
         ENDDO
      ENDIF  ! MULPHL

      IF(NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_SSEND(FILE2(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF

      CALL MPI_SSEND(ROAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RMAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RNAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RWAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RKAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(EPSAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,
     +               IERR)
      CALL MPI_SSEND( EAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)

      IF(NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_SSEND(FIAV1(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF

      CALL MPI_SSEND(ROAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RMAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RNAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RWAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( PAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( TAV1(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( EAV2(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RUVAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RUWAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RVWAV(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)

      IF(NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_SSEND(FIAV2(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF
      
      ELSEIF(MASTER .AND. N == -1) THEN  ! Receive and write
      
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,MPI_COMM_WORLD,STATUS,
     +        IERR)
         ALLOCATE(SBUF(NNMAX))
         ISUM = 0

         IF(MULPHL) ISUM = 2*NPHASES

         DO I = 1,3*(7+NSCAL)+3+ISUM
             CALL MPI_RECV(SBUF,NNMAX,MPI_REAL8,NR1,K,MPI_COMM_WORLD,
     +           STATUS,IERR)
            WRITE(20) (SBUF(J),J=1,NNMAX)
         ENDDO

         DEALLOCATE(SBUF)

      ENDIF  ! (MASTER .AND. N /= -1)

      ENDDO  ! NCG = 1,NBLOCG   

      ENDIF  ! TIMEL

C ... End time accurate data
***************************************************************************

      IF(COORL) THEN                 ! Coordinates are changing

      DO NCG = 1,NBLOCG              ! Loop over the global blocks

      N     = NLOCAL(NCG)            ! local block number
      NROOT = NPNUM(NCG)             ! process number
      NR1   = NROOT - 1              ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN  ! Master's own blocks
      
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20)(XLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(YLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(ZLE2(IG1+I),I=1,NNMAX)
         WRITE(20)(XORI(IG1+I),I=1,NNMAX)
         WRITE(20)(YORI(IG1+I),I=1,NNMAX)
         WRITE(20)(ZORI(IG1+I),I=1,NNMAX)
	 
      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN  ! Slaves send to the master
      
         IG1   = IG(1,N)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         K = NCG
         CALL MPI_SEND(NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(XLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(YLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(ZLE2(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(XORI(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(YORI(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_SSEND(ZORI(IG1),NNMAX,MPI_REAL8,0,K,
     +                  MPI_COMM_WORLD,IERR)
     
      ELSEIF(MASTER .AND. N == -1) THEN  ! Master receives and writes
      
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)

         ALLOCATE(DBUF(NNMAX))

         DO I = 1,6
            CALL MPI_RECV(DBUF,NNMAX,MPI_REAL8,NR1,K,
     +           MPI_COMM_WORLD,STATUS,IERR)
            WRITE(20) (DBUF(J),J=1,NNMAX)
         ENDDO

         DEALLOCATE(DBUF)

      ENDIF  ! (MASTER .AND. N /= -1)

      ENDDO  ! NCG = 1,NBLOCG   

      ENDIF  ! COORL

C ... End coordinates changed
****************************************************************************

      IF(CHIMEL) THEN                ! Chimera blocks
      
      DO NCG = 1,NBLOCG              ! Loop over the global blocks

      N     = NLOCAL(NCG)            ! local block number
      NROOT = NPNUM(NCG)             ! process number
      NR1   = NROOT - 1              ! process number used by MPI

      IF(MASTER .AND. N /= -1) THEN  ! Master's own blocks
      
         IG1   = IG(1,N)-IG(1,1)
         NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
         WRITE(20)( ROFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RMFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RNFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RWFOR(IG1+I),I=1,NNMAX)
         WRITE(20)(  EFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( PDFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( RKFOR(IG1+I),I=1,NNMAX)
         WRITE(20)( REFOR(IG1+I),I=1,NNMAX)
         WRITE(20)(  RKSI(IG1+I),I=1,NNMAX)

         IF(TRANSL) THEN  ! Intermittency variables
            WRITE(20) (TRM(IG1+I)%GFOR,   I=1,NNMAX)
            WRITE(20) (TRM(IG1+I)%RETFOR, I=1,NNMAX)
         ENDIF

         IF(MULPHL) THEN 
            DO IPHASE = 1,NPHASES
               WRITE(20) (MPFOR(IG1+I)%ALFAFOR(IPHASE),  I=1,NNMAX)
               WRITE(20) (MPFOR(IG1+I)%XFOR(IPHASE),     I=1,NNMAX)
               WRITE(20) (MPFOR(IG1+I)%DTEMPFOR(IPHASE), I=1,NNMAX)
            ENDDO
         ENDIF

         IF(NSCAL > 0) THEN
            DO NS = 1,NSCAL
               WRITE(20) (FIFOR(IG1+I,NS),I=1,NNMAX)
            ENDDO
         ENDIF
	 
      ELSEIF(.NOT.MASTER .AND. N /= -1) THEN ! Slaves send to the master
      
      IG1   = IG(1,N)
      NNMAX = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)
      K = NCG
      CALL MPI_SEND(  NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(ROFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RMFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RNFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RWFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( EFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(PDFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(RKFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND(REFOR(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)
      CALL MPI_SSEND( RKSI(IG1),NNMAX,MPI_REAL8,0,K,MPI_COMM_WORLD,IERR)

      IF(TRANSL) THEN
         F1R(IG1:IGN) = TRM(IG1:IGN)%GFOR
         CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,IERR)
         F1R(IG1:IGN) = TRM(IG1:IGN)%RETFOR
         CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +        MPI_COMM_WORLD,IERR)
      ENDIF

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
            F1R(IG1:IGN) = MPFOR(IG1:IGN)%ALFAFOR(IPHASE)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = MPFOR(IG1:IGN)%XFOR(IPHASE)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
            F1R(IG1:IGN) = MPFOR(IG1:IGN)%DTEMPFOR(IPHASE)
            CALL MPI_SSEND(F1R(IG1:IGN),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF

      IF(NSCAL > 0) THEN
         DO NS = 1,NSCAL
            CALL MPI_SSEND(FIFOR(IG1,NS),NNMAX,MPI_REAL8,0,K,
     +                     MPI_COMM_WORLD,IERR)
         ENDDO
      ENDIF
      
      ELSEIF(MASTER .AND. N == -1) THEN  ! Master receives and writes
      
         K = NCG
         CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,
     +                 MPI_COMM_WORLD,STATUS,IERR)
         ALLOCATE(SBUF(NNMAX))

         NTRFOR = 0
         IF(TRANSL) NTRFOR = 2

         DO I = 1,9+NTRFOR+NSCAL
            CALL MPI_RECV(SBUF,NNMAX,MPI_REAL8,NR1,K,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            WRITE(20) (SBUF(J),J=1,NNMAX)
         ENDDO

         DEALLOCATE(SBUF)

      ENDIF  ! (MASTER .AND. N /= -1)

      ENDDO  ! NCG = 1,NBLOCG   

      ENDIF  ! CHIMEL

C ... End Chimera blocks
**********************************************************************      

      RETURN
      END SUBROUTINE WRITRS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITWH

      USE MPI

      USE INTEGERS, ONLY    : IPRO, NBCS
      USE NS3CO, ONLY       : PARALLEL, IN, JN, KN, LN, IC9
      USE MAIN_ARRAYS, ONLY : WAVEH, ICON, IHF

      IMPLICIT NONE

      INTEGER :: IP,IPW,IPL,IPG,IFACE,IBC,KX1,KX2,KY1,KY2,NR,NP,NCG,
     +           I,J,II,IERR,ISTRID,JSTRID,NNMAX,NBCSG,ITSME,OWNER

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      REAL, ALLOCATABLE, DIMENSION(:) :: WAVES

      LOGICAL :: MASTER, FOUNDL, FOUNDG


      MASTER = IPRO == 1

      NBCSG  = NBCS

      IF(PARALLEL) CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,
     &                                MPI_SUM,MPI_COMM_WORLD,IERR)

 
      DO IPW = 1,NBCSG              ! Global patch loop

         FOUNDL = .FALSE.
         ITSME  = 0

         DO IP = 1,NBCS             ! Local patch loop

         II  = 0                    ! Free-surface pointer
         IBC = ICON(IC9*(IP-1)+ 1)  ! Patch type
         IPG = ICON(IC9*(IP-1)+25)  ! Global patch number

         IF(IPW == IPG .AND. IBC == 13) THEN

            FOUNDL = .TRUE.

            IPL    = ICON(IC9*(IP-1)+ 2)  ! Proces local patch number
            IFACE  = ICON(IC9*(IP-1)+ 3)  ! Face number
            NCG    = ICON(IC9*(IP-1)+24)  ! Global block number
 
            IF(PARALLEL) CALL MPI_COMM_RANK(MPI_COMM_WORLD,ITSME,IERR)
         
C ... Patch indeces without extensions

            IF(IFACE == 2 .OR. IFACE == 5) THEN
               KY1 = ICON(IC9*(IP-1)+4) !+ 1
               KY2 = ICON(IC9*(IP-1)+5) !- 1
               KX1 = ICON(IC9*(IP-1)+6) !+ 1
               KX2 = ICON(IC9*(IP-1)+7) !- 1
            ELSE
               KX1 = ICON(IC9*(IP-1)+4) !+ 1
               KX2 = ICON(IC9*(IP-1)+5) !- 1
               KY1 = ICON(IC9*(IP-1)+6) !+ 1
               KY2 = ICON(IC9*(IP-1)+7) !- 1
            ENDIF

*            ISTRID = KX2-KX1+1 + 2*IN
*            JSTRID = KY2-KY1+1 + 2*JN
            ISTRID = KX2-KX1+1 + 2*LN
            JSTRID = KY2-KY1+1 + 2*LN
            NNMAX  = ISTRID*JSTRID

            ALLOCATE (WAVES(NNMAX),STAT=IERR)

*            DO J = KY1-JN,KY2+JN
            DO J = KY1-2,KY2+2
               NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL,1)
*               DO I = KX1-IN,KX2+IN
               DO I = KX1-2,KX2+2
                  NP =  I + NR                ! Patch index with LN ghost cells
                  II = II + 1 
*                  WRITE(212,*)I,J,NP-ihf(ipl)+1,NP,WAVEH(NP)
                  WAVES(II) = WAVEH(NP)
               ENDDO
            ENDDO

         ENDIF                  ! IBC == 13 (FRE)

      ENDDO                     ! Local patch loop


      IF(PARALLEL) THEN

         CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,
     &                      MPI_LOR,MPI_COMM_WORLD,IERR)

         IF(FOUNDG) CALL MPI_ALLREDUCE(ITSME,OWNER,1,MPI_INTEGER,
     &                                 MPI_SUM,MPI_COMM_WORLD,IERR)

         IF(MASTER .AND. FOUNDL) THEN
            WRITE(21,*)  NNMAX,NCG         
            WRITE(21,*) (WAVES(II),II=1,NNMAX)
         ELSEIF(.NOT.MASTER .AND. FOUNDL) THEN
            CALL MPI_SEND(NCG,1,MPI_INTEGER,0,IPW,
     +                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(NNMAX,1,MPI_INTEGER,0,IPW,
     +                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(WAVES,NNMAX,MPI_REAL8,0,IPW,
     +                    MPI_COMM_WORLD,IERR)
         ELSEIF(MASTER .AND. .NOT.FOUNDL .AND. FOUNDG) THEN
            CALL MPI_RECV(NCG,1,MPI_INTEGER,OWNER,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            CALL MPI_RECV(NNMAX,1,MPI_INTEGER,OWNER,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            ALLOCATE (WAVES(NNMAX),STAT=IERR)
            CALL MPI_RECV(WAVES,NNMAX,MPI_REAL8,OWNER,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            WRITE(21,*) NNMAX,NCG         
            WRITE(21,*) (WAVES(II),II=1,NNMAX)         
            DEALLOCATE (WAVES)
         ENDIF

      ELSE                      ! Single process run

         IF(FOUNDL) THEN
            WRITE(21,*)  NNMAX,NCG         
            WRITE(21,*) (WAVES(II),II=1,NNMAX)
         ENDIF

      ENDIF

      IF(ALLOCATED(WAVES)) DEALLOCATE(WAVES)

      ENDDO                     ! Global patch loop   

      END SUBROUTINE WRITWH
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE READWH(INITVL)

      USE MPI

      USE INTEGERS, ONLY    : IPRO, NBCS
      USE NS3CO, ONLY       : PARALLEL, IN, JN, KN, LN, IC9
      USE MAIN_ARRAYS, ONLY : WAVEH, ICON, IHF, NPNUM

      IMPLICIT NONE

      INTEGER :: IP,IPW,IPL,IPG,IFACE,IBC,KX1,KX2,KY1,KY2,NR,NP,NCG,
     +           I,J,II,IERR,ISTRID,JSTRID,NNMAX,NBCSG,TARGET,
     +           NP1,NPJ,NPJ1

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      REAL, ALLOCATABLE, DIMENSION(:) :: WAVES

      LOGICAL :: MASTER, FOUNDL, FOUNDG, INITVL

      MASTER = IPRO == 1

      NBCSG  = NBCS

      IF(PARALLEL) CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,
     &                                MPI_SUM,MPI_COMM_WORLD,IERR)


      DO IPW = 1,NBCSG            ! Global patch loop

         FOUNDL = .FALSE.

         DO IP = 1,NBCS           ! Local patch loop

         IBC = ICON(IC9*(IP-1)+ 1) ! Patch type
         IPG = ICON(IC9*(IP-1)+25) ! Global patch number

         IF(IPW == IPG .AND. IBC == 13) THEN

            FOUNDL = .TRUE.

            IPL    = ICON(IC9*(IP-1)+ 2) ! Proces local patch number
            IFACE  = ICON(IC9*(IP-1)+ 3) ! Face number
            NCG    = ICON(IC9*(IP-1)+24) ! Global block number

C ... Patch indeces without extensions

            IF(IFACE == 2 .OR. IFACE == 5) THEN
               KY1 = ICON(IC9*(IP-1)+4) !+ 1
               KY2 = ICON(IC9*(IP-1)+5) !- 1
               KX1 = ICON(IC9*(IP-1)+6) !+ 1
               KX2 = ICON(IC9*(IP-1)+7) !- 1
               IF(INITVL) THEN
                  KY1 = (KY1+1)/2
                  KY2 =  KY2/2
                  KX1 = (KX1+1)/2
                  KX2 =  KX2/2
               ENDIF
            ELSE
               KX1 = ICON(IC9*(IP-1)+4) !+ 1
               KX2 = ICON(IC9*(IP-1)+5) !- 1
               KY1 = ICON(IC9*(IP-1)+6) !+ 1
               KY2 = ICON(IC9*(IP-1)+7) !- 1
               IF(INITVL) THEN
                  KX1 = (KX1+1)/2
                  KX2 =  KX2/2
                  KY1 = (KY1+1)/2
                  KY2 =  KY2/2
               ENDIF
            ENDIF

*            ISTRID = KX2-KX1+1 + 2*IN
*            JSTRID = KY2-KY1+1 + 2*JN
            ISTRID = KX2-KX1+1 + 2*LN
            JSTRID = KY2-KY1+1 + 2*LN

            NNMAX  = ISTRID*JSTRID

         ENDIF                  ! IBC == 13 (FRE) 

      ENDDO                     ! Local patch loop


      IF(PARALLEL) THEN

         CALL MPI_REDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,
     &                   MPI_LOR,0,MPI_COMM_WORLD,IERR)

         IF(MASTER .AND. FOUNDL) THEN
            READ(21,*)  NNMAX,NCG         
            ALLOCATE (WAVES(NNMAX),STAT=IERR)            
            READ(21,*) (WAVES(II),II=1,NNMAX)
         ELSEIF(MASTER .AND. .NOT.FOUNDL .AND. FOUNDG) THEN
            READ(21,*)  NNMAX,NCG         
            ALLOCATE (WAVES(NNMAX),STAT=IERR)            
            READ(21,*) (WAVES(II),II=1,NNMAX)
            TARGET = NPNUM(NCG) - 1
            CALL MPI_SEND(NCG,1,MPI_INTEGER,TARGET,IPW,
     +                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(NNMAX,1,MPI_INTEGER,TARGET,IPW,
     +                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(WAVES,NNMAX,MPI_REAL8,TARGET,IPW,
     +                    MPI_COMM_WORLD,IERR)
         ELSEIF(.NOT.MASTER .AND. FOUNDL) THEN
            CALL MPI_RECV(NCG,1,MPI_INTEGER,0,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            CALL MPI_RECV(NNMAX,1,MPI_INTEGER,0,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
            ALLOCATE (WAVES(NNMAX),STAT=IERR)
            CALL MPI_RECV(WAVES,NNMAX,MPI_REAL8,0,IPW,
     +                    MPI_COMM_WORLD,STATUS,IERR)
         ENDIF

      ELSE                      ! Single process run

         IF(FOUNDL) THEN
            READ(21,*)  NNMAX,NCG         
            ALLOCATE (WAVES(NNMAX),STAT=IERR)            
            READ(21,*) (WAVES(II),II=1,NNMAX)
         ENDIF

      ENDIF

      IF(FOUNDL .AND. .NOT.INITVL) THEN

         II = 0
*         DO J = KY1-JN,KY2+JN
         DO J = KY1-2,KY2+2
            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL,1)
*            DO I = KX1-IN,KX2+IN
            DO I = KX1-2,KX2+2
               NP = I + NR                    ! Patch index with LN ghost cells
               II = II + 1
               WAVEH(NP) = WAVES(II)
            ENDDO
c            IF(INITVL) II = II + (JN-1)*2
         ENDDO

      ELSEIF(FOUNDL .AND. INITVL) THEN

*         II = 1 + KX2-KX1+1+2*IN
         II = 1 + KX2-KX1+1+2*LN
         DO J = KY1-1,KY2+1
*            NR = (JN+2*J-2*KY1)*(2*KX2-2*KX1+2+2*IN)-2*KX1+IN+IHF(IPL,1)
            NR = (LN+2*J-2*KY1)*(2*KX2-2*KX1+2+2*LN)-2*KX1+LN+IHF(IPL,1)
            DO I = KX1-1,KX2+1
               NP   = 2*I + NR                ! Patch index with LN ghost cells
               NP1  = NP  + 1
*               NPJ  = NP  + 2*KX2-2*KX1+2+2*IN
               NPJ  = NP  + 2*KX2-2*KX1+2+2*LN
               NPJ1 = NPJ + 1
               II   = II  + 1
               WAVEH(NP)   = WAVES(II)
               WAVEH(NP1)  = WAVES(II)
               WAVEH(NPJ)  = WAVES(II)
               WAVEH(NPJ1) = WAVES(II)
            ENDDO
*            II = II + (JN-1)*2
            II = II + (LN-1)*2
         ENDDO

      ENDIF

      IF(ALLOCATED(WAVES)) DEALLOCATE(WAVES)

      ENDDO                     ! Global patch loop   

      END SUBROUTINE READWH
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EOF(ICHAN,ICYCLE)

C ... THIS SUBROUTINE WILL FIND THE END OF A FILE

      IMPLICIT NONE

      INTEGER :: ICHAN, ICYCLE
      REAL    :: APU
      LOGICAL :: END_OF_FILE

      END_OF_FILE = .FALSE.

      DO WHILE (.NOT.END_OF_FILE)
         READ(ICHAN,*,END=20,ERR=20) APU
      ENDDO

 20   CONTINUE

      BACKSPACE ICHAN

      RETURN
      END SUBROUTINE EOF      
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WTCFOR(ICHAN,ICYCLE)

C ... Correctly continue writing the iteration history of force groups, 
C ... i.e., scan the IBDIAG_FGS file until the line containing ICYCLE
C ... is found and continue writing there.

      IMPLICIT NONE

      REAL    :: APU 
      INTEGER :: ITERA, ICHAN, ICYCLE

      ITERA = 0

      DO WHILE (ITERA <= ICYCLE)
         READ(ICHAN,*,END=10,ERR=20) APU,ITERA
      ENDDO

 10   CONTINUE

      BACKSPACE ICHAN

 20   CONTINUE

      RETURN
      END SUBROUTINE WTCFOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOCATE_LINE(ICHAN,ICYCLE)

C ... Locate line containing current cycle in an ASCII file.

      IMPLICIT NONE

      INTEGER :: ITERA, ICHAN, ICYCLE

      ITERA = 0

      DO WHILE (ITERA <= ICYCLE)
         READ(ICHAN,*,END=20,ERR=20) ITERA
      ENDDO

 20   BACKSPACE ICHAN

* 20   CONTINUE

      RETURN
      END SUBROUTINE LOCATE_LINE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WTC(ICHAN,ICYCLE)

C ... Correctly continue writing the iteration history, i.e.,
C ... scan the IBDIAG file until the line containing ICYCLE is found
C ... and continue writing there.
C ... File changed to binary-file 20.2.1998, PK.

      IMPLICIT NONE
      CHARACTER(LEN=79) :: NAME
      INTEGER :: MAJOR, MINOR, ITERA, ICHAN, ICYCLE, IDRXX,ITURB,NRESI,
     & KSCAL
      REAL(KIND=4) :: RMACH,ALPHA,RE,T

      ITERA = 0
      READ(ICHAN) MAJOR,MINOR   ! IBDIAG v. Major.Minor
      READ(ICHAN) NAME

      READ(ICHAN) RMACH,ALPHA,RE,T,ITURB,KSCAL,NRESI !,IDRXX ! Ei haittaa!

      DO WHILE (ITERA <= ICYCLE)
         READ(ICHAN,END=10,ERR=20) ITERA
      ENDDO

 10   CONTINUE

      BACKSPACE ICHAN

 20   CONTINUE

      RETURN
      END SUBROUTINE WTC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WTCFGS(ICHAN,ICYCLE)

C ... Correctly continue writing the iteration history of force groups, 
C ... i.e., scan the IBDIAG_FGS file until the line containing ICYCLE
C ... is found and continue writing there.

      IMPLICIT NONE

      INTEGER :: ITERA, ICHAN, ICYCLE

      ITERA = 0

      DO WHILE (ITERA <= ICYCLE)
         READ(ICHAN,END=10,ERR=20) ITERA
      ENDDO

 10   CONTINUE

      BACKSPACE ICHAN

 20   CONTINUE

      RETURN
      END SUBROUTINE WTCFGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WTCBAK(ICHAN,ICHANB,ICYCLE)

C     Correctly continue writing the iteration history, i.e.,
C     scan the IBDIAG file until the line containing ICYCLE is found
C     and continue writing there.
C     File changed to binary-file 20.2.1998, PK.

      IMPLICIT NONE
      CHARACTER(LEN=79) :: NAME
      INTEGER :: MAJOR,MINOR,ITERA,ICHAN,ICYCLE,ICHANB,I,IQ,ITURB,
     2           KSCAL,NRESI
      REAL(KIND=4), ALLOCATABLE :: QQ1(:),QQ2(:)
      REAL(KIND=4) :: RMACH,ALPHA,RE,T

      REWIND(ICHAN)
      REWIND(ICHANB) ! Kokeilu 31.12.03

      ITERA = 0

      READ(ICHAN) MAJOR,MINOR   ! IBDIAG v. Major.Minor
      READ(ICHAN) NAME
      READ(ICHAN) RMACH,ALPHA,RE,T,ITURB,KSCAL,NRESI

      WRITE(ICHANB) MAJOR,MINOR   ! IBDIAG v. Major.Minor
      WRITE(ICHANB) NAME
      WRITE(ICHANB) RMACH,ALPHA,RE,T,ITURB,KSCAL,NRESI

      ALLOCATE (QQ1(7),QQ2(NRESI+7))

      DO WHILE (ITERA < ICYCLE)
         READ (ICHAN,END=10,ERR=10)  ITERA,(QQ1(I),I=1,7),IQ,
     &        (QQ2(I),I=1,NRESI+7)
         WRITE(ICHANB) ITERA,(QQ1(I),I=1,7),IQ,
     &        (QQ2(I),I=1,NRESI+7)
      ENDDO

      DEALLOCATE (QQ1, QQ2)
      RETURN
10    CONTINUE
      WRITE(45,*) 'OBS! IBDIAG ENDED IN BACK-UP PROCESS WTCBAK'
      DEALLOCATE (QQ1, QQ2)

      RETURN
      END SUBROUTINE WTCBAK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITE_BNDFILE(NBLOCK,IBC,NBCS,BCFILE,IPTRFG,
     &                         LEVEL,ENSBND_FILE,FVBND_FILE)

C ... This subroutine writes boundary information for EnSight and
C ... for Fieldview.

      USE NS3CO, ONLY : IC9, PLOT3D, ENSIGHT

      IMPLICIT NONE

      INTEGER :: IBC(IC9,*),IPTRFG(52), NORMAL,NBG,ITYPE,NBLOCK,
     +           NBCS,LEVEL,NDIV,I,J,K,IMOV,IM,IX,JM,JX,KM,KX

      CHARACTER(LEN=80) :: BCFILE(*)
      CHARACTER(LEN=80) :: ENSBND_FILE, FVBND_FILE

      LOGICAL :: HAVRES, INLOUT, SOLIDS, FRESUL, FRESU

      IF(ENSIGHT) OPEN(33,FILE = ENSBND_FILE)
      IF(PLOT3D)  OPEN(66,FILE = FVBND_FILE)

      IF(ENSIGHT) WRITE(33,664)
      IF(ENSIGHT) WRITE(33,666)
      IF(PLOT3D)  WRITE(66,665)
      IF(PLOT3D)  WRITE(66,666)

      K = 0
      IPTRFG(1:52) = 0

      DO J=1,52  ! Uppercase and lowercase alphabetics / Force groups
         DO I=1,NBCS

         INLOUT = (IBC(1,I) >= 3.AND.  IBC(1,I) <= 5   ! Inlet/outlet or MIR
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 1   ! Periodic i/o surface
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 7   ! Minor loss
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 8   ! Actuator disc
     +       .OR. IBC(1,I) == 17)                      ! Circulating inlet
         SOLIDS = (IBC(1,I) >= 8 .AND. IBC(1,I) <= 10
     +        .OR. IBC(1,I) == 15)                     ! Solid surface
         FRESU  = (IBC(1,I) == 13)                     ! Free surface

         IF(INLOUT .OR. SOLIDS .OR. FRESU) THEN       ! A solid or i/o surface 
               IF(BCFILE(I)(J:J) /= ' ') THEN
                  IF(J <= 26 .AND. ENSIGHT) WRITE(33,667) CHAR(J+64) 
                  IF(J > 26 .AND. ENSIGHT) WRITE(33,667) CHAR(J+70)
                  IF(J <= 26 .AND. PLOT3D)  WRITE(66,667) CHAR(J+64) 
                  IF(J > 26 .AND. PLOT3D)  WRITE(66,667) CHAR(J+70)
                  K = K + 1
                  IPTRFG(J) = K 
                  GOTO 975
               ENDIF
            ENDIF
         ENDDO
 975  CONTINUE
      ENDDO
      
      IF(ENSIGHT) WRITE(33,668)
      IF(PLOT3D)  WRITE(66,668)

C ... Check this if boundary types are increased in XYZ.BIN.fvbnd (15 options)

 664  FORMAT('ENSBND 1.0')
 665  FORMAT('FVBND 1 4')
 666  FORMAT('Solid wall'/
     +     'Rotating wall'/
     +     'Moving wall(1)'/
     +     'Moving wall(2)'/
     +     'Moving wall(3)'/
     +     'Solidplus1'/
     +     'Symmetry'/
     +     'Inlet'/
     +     'Outlet'/
     +     'Solid-fluid coupling'/
     +     'Free surface'/
     +     'Arbitrarily moving wall(4)'/
     +     'Slip boundary'/
     +     'Periodic i/o'/
     +     'Circulating inlet'/
     +     'Minor pressure loss'/
     +     'Actuator disc')
 667  FORMAT('Zone: Force group ',A1)
 668  FORMAT('BOUNDARIES')


      HAVRES = .TRUE.
      NDIV = 2**(LEVEL-1)

       DO I = 1,NBCS

         NBG  = IBC(24,I)                              ! Global block number
         IMOV = IBC(9,I)                               ! Movement type

         INLOUT =(IBC(1,I) == 3 .OR.  IBC(1,I)  == 5   ! Inlet or outlet
     +       .OR. IBC(1,I) == 4                        ! MIR
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 1   ! Periodic i/o surface
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 7   ! Minor loss
     +       .OR. IBC(1,I) == 1 .AND. IBC(20,I) == 8   ! Actuator disc
     +       .OR. IBC(1,I) == 17)                      ! Circulating inlet

         SOLIDS = (IBC(1,I) >= 8 .AND. IBC(1,I) <= 10
     +        .OR. IBC(1,I) == 15)                     ! Solid surface

         FRESUL = (IBC(1,I) == 13)                     ! Free surface

         IF(INLOUT .OR. SOLIDS .OR. FRESUL) THEN       ! A solid or i/o surface 

         IF(IBC(24,I) <= NBLOCK) THEN

            IF(IBC(1,I) ==  1) ITYPE = 14
            IF(IBC(1,I) ==  3) ITYPE =  8
            IF(IBC(1,I) ==  4) ITYPE =  7
            IF(IBC(1,I) ==  5) ITYPE =  9
            IF(IBC(1,I) ==  8) ITYPE =  1
            IF(IBC(1,I) ==  9) ITYPE =  2
            IF(IBC(1,I) == 13) ITYPE = 11
            IF(IBC(1,I) == 15) ITYPE = 10
            IF(IBC(1,I) == 16) ITYPE =  3
            IF(IBC(1,I) == 17) ITYPE = 15
            IF(IBC(1,I) ==  1 .AND. IBC(20,I) == 7) ITYPE = 16
            IF(IBC(1,I) ==  1 .AND. IBC(20,I) == 8) ITYPE = 17
            IF(IBC(1,I) == 10 .AND. IMOV == 1)      ITYPE =  3
            IF(IBC(1,I) == 10 .AND. IMOV == 2)      ITYPE =  4
            IF(IBC(1,I) == 10 .AND. IMOV == 3)      ITYPE =  5
            IF(IBC(1,I) == 10 .AND. IMOV == 4)      ITYPE = 12

            IM = (IBC(4,I)-1)/NDIV + 1
            IX = (IBC(5,I)-0)/NDIV + 1
            JM = (IBC(6,I)-1)/NDIV + 1
            JX = (IBC(7,I)-0)/NDIV + 1
            KM = (IBC(6,I)-1)/NDIV + 1
            KX = MAX0((IBC(7,I)-0)/NDIV+1,2) ! 2D case

            IF(IBC(3,I) == 1 .OR. IBC(3,I) == 4) THEN
               JM = (IBC(4,I)-1)/NDIV + 1
               JX = (IBC(5,I)-0)/NDIV + 1
            ENDIF
           
            IF(IBC(3,I) == 1) THEN
               IM = 1
               IX = 1
               NORMAL = 1
            ELSEIF(IBC(3,I) == 2) THEN
               JM = 1
               JX = 1
               NORMAL = 1
            ELSEIF(IBC(3,I) == 3) THEN
               KM = 1
               KX = 1
               NORMAL = 1
            ELSEIF(IBC(3,I) == 4) THEN
               IM = IBC(21,I)/NDIV + 1
               IX = IBC(21,I)/NDIV + 1
               NORMAL = -1
            ELSEIF(IBC(3,I) == 5) THEN
               JM = IBC(21,I)/NDIV + 1
               JX = IBC(21,I)/NDIV + 1
               NORMAL = -1
            ELSEIF(IBC(3,I) == 6) THEN
               KM = MAX(1,IBC(21,I)/NDIV) + 1
               KX = MAX(1,IBC(21,I)/NDIV) + 1
               NORMAL = -1
            ENDIF

            IF(ENSIGHT) WRITE(33,'(8I5)')
     &      ITYPE,IBC(24,I),IM,IX,JM,JX,KM,KX
            IF(PLOT3D)  WRITE(66,'(8I5,L5,I5)')
     &      ITYPE,IBC(24,I),IM,IX,JM,JX,KM,KX,HAVRES,NORMAL


C ... Force group patches. NOTE: IPTRFG(J)+17, where 17 is the number of 
C ... previous boundary types in XYZ.BIN.fvbnd.
C ... Check this if boundary types are increased in XYZ.BIN.fvbnd

            DO J=1,52  ! Uppercase and lowercase alphabetics 
               IF(BCFILE(I)(J:J) /= ' ' .AND. IPTRFG(J) /= 0) THEN
                  IF(ENSIGHT) WRITE(33,'(8I5)') IPTRFG(J)+17,
     &            IBC(24,I),IM,IX,JM,JX,KM,KX
                  IF(PLOT3D)  WRITE(66,'(8I5,L5,I5)') IPTRFG(J)+17,
     &            IBC(24,I),IM,IX,JM,JX,KM,KX,.FALSE.,NORMAL
               ENDIF
            ENDDO

C ... Solid + 1 patches.

            IF(SOLIDS) THEN 

               ITYPE = 6            

               IM = (IBC(4,I)-1)/NDIV + 1
               IX = (IBC(5,I)-0)/NDIV + 1
               JM = (IBC(6,I)-1)/NDIV + 1
               JX = (IBC(7,I)-0)/NDIV + 1
               KM = (IBC(6,I)-1)/NDIV + 1
               KX = MAX0((IBC(7,I)-0)/NDIV+1,2) ! 2D case

               IF(IBC(3,I) == 1 .OR. IBC(3,I) == 4) THEN
                  JM = (IBC(4,I)-1)/NDIV + 1
                  JX = (IBC(5,I)-0)/NDIV + 1
               ENDIF
           
               IF(IBC(3,I) == 1) THEN
                  IM = 2
                  IX = 2
                  NORMAL = 1
               ELSEIF(IBC(3,I) == 2) THEN
                  JM = 2
                  JX = 2
                  NORMAL = 1
               ELSEIF(IBC(3,I) == 3) THEN
                  KM = 2
                  KX = 2
                  NORMAL = 1
               ELSEIF(IBC(3,I) == 4) THEN
                  IM = IBC(21,I)/NDIV
                  IX = IBC(21,I)/NDIV
                  NORMAL = -1
               ELSEIF(IBC(3,I) == 5) THEN
                  JM = IBC(21,I)/NDIV
                  JX = IBC(21,I)/NDIV
                  NORMAL = -1
               ELSEIF(IBC(3,I) == 6) THEN
                  KM = IBC(21,I)/NDIV
                  KX = IBC(21,I)/NDIV
                  NORMAL = -1
               ENDIF

               IF(ENSIGHT) WRITE(33,'(8I5)')
     &         ITYPE,IBC(24,I),IM,IX,JM,JX,KM,KX
               IF(PLOT3D)  WRITE(66,'(8I5,L5,I5)')
     &         ITYPE,IBC(24,I),IM,IX,JM,JX,KM,KX,.FALSE.,NORMAL

            ENDIF

         ENDIF
         ENDIF ! Solid or MIR or i/o surface
      ENDDO

      IF(ENSIGHT) CLOSE(33)
      IF(PLOT3D)  CLOSE(66)

      RETURN
      END SUBROUTINE WRITE_BNDFILE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOADD(NB,IDMN,JDMN,X,Y,Z,FX,FY,FZ,FFX,FFY,FFZ,TEMP,A,
     &                 AREF,CHLREF,REFPRE)
                                                               
C ... This subroutine determines the load distributions for 
C ... selected combinations of surface pathes.

      USE CONSTANTS,   ONLY : EPS 

      IMPLICIT NONE   

      REAL, PARAMETER :: TL=1.0E-3

      INTEGER, PARAMETER :: NSLOT=100
      INTEGER, DIMENSION(NSLOT,1,1,2) :: LBOXES
      INTEGER, DIMENSION(1,NSLOT,1,2) :: MBOXES
      INTEGER, DIMENSION(1,1,NSLOT,2) :: NBOXES

      INTEGER, DIMENSION(*)           :: IDMN, JDMN

      INTEGER :: NNODE, NOFEL
      INTEGER :: I, J, K, NI, NJ, NK, IB, NB, IERRCODE, IC, IL, IE, IM 
      INTEGER :: NCB, LC, LBS, LBE
      INTEGER :: ILOW, IUPP, JLOW, JUPP, IP, NPOINTS, I1, J1, K1, ISKIN

      INTEGER, DIMENSION(:), ALLOCATABLE :: KSC, NR

      REAL, DIMENSION(*) :: X, Y, Z
      REAL, DIMENSION(*) :: FX, FY, FZ, FFX, FFY, FFZ, TEMP, A
      REAL, DIMENSION(:,:), ALLOCATABLE :: SKINX, SKINY, SKINZ
      REAL :: AREF, CHLREF, REFPRE

      REAL :: XI,YI,ZI, X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3, X4,Y4,Z4, 
     &                  X5,Y5,Z5, XN,YN,ZN, DD, DDMIN

      REAL :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      REAL :: DX, DY, DZ, DXR, DYR, DZR

      REAL, DIMENSION(NSLOT,2) :: FXX, FYX, FZX
      REAL, DIMENSION(NSLOT,2) :: FXY, FYY, FZY
      REAL, DIMENSION(NSLOT,2) :: FXZ, FYZ, FZZ

      REAL, DIMENSION(NSLOT,2) :: FFXX, FFYX, FFZX
      REAL, DIMENSION(NSLOT,2) :: FFXY, FFYY, FFZY
      REAL, DIMENSION(NSLOT,2) :: FFXZ, FFYZ, FFZZ

      REAL :: MINFXX, MAXFXX, MINFYX, MAXFYX, MINFZX, MAXFZX
      REAL :: MINFXY, MAXFXY, MINFYY, MAXFYY, MINFZY, MAXFZY
      REAL :: MINFXZ, MAXFXZ, MINFYZ, MAXFYZ, MINFZZ, MAXFZZ

      REAL, DIMENSION(3,5) :: A1, B1, C1
      REAL                 :: ASHARE, ASHARE1, ASHARE2
      REAL                 :: APROJ,  APROJ1,  APROJ2, AELEM
  
      REAL :: XC, YC, ZC
      REAL :: CX, CY, CZ
      REAL :: CXM, CYM, CZM, XCP, YCP, ZCP
      REAL :: QKINS
      REAL :: SUMX, SUMY, SUMZ
      REAL :: FRICSUMX, FRICSUMY, FRICSUMZ

      LOGICAL :: FLATX, FLATY, FLATZ


      QKINS = REFPRE*AREF


      NNODE = 0
      NOFEL = 0

      DO IB=1,NB
         NNODE = NNODE + (IDMN(IB)+1)*(JDMN(IB)+1)
         NOFEL = NOFEL + IDMN(IB)*JDMN(IB)
      ENDDO

C ......................................................................

C ... Allocate memory for element corner point arrays

      ALLOCATE(SKINX(4,NOFEL),STAT=IERRCODE)
      ALLOCATE(SKINY(4,NOFEL),STAT=IERRCODE)
      ALLOCATE(SKINZ(4,NOFEL),STAT=IERRCODE)

C ... Not enough memory ? 

      IF(IERRCODE /= 0) THEN                                          
         WRITE(*,*)'LOADD:  Not enough memory, aborting ...'   
         STOP
      ENDIF                                                           

C ......................................................................

C ... Fill the element corner point arrays

      IL = 0
      IM = 0
      DO IB=1,NB    
         NI = IDMN(IB)+1
         NJ = JDMN(IB)+1
         DO J=1,NJ-1
            DO I=1,NI-1
               IC = IL + I + (J-1)*NI       ! Node point number
               IE = IM + I + (J-1)*(NI-1)   ! Element number
               SKINX(1,IE) = X(IC)
               SKINY(1,IE) = Y(IC)
               SKINZ(1,IE) = Z(IC)
               SKINX(2,IE) = X(IC+1)
               SKINY(2,IE) = Y(IC+1)
               SKINZ(2,IE) = Z(IC+1)
               SKINX(3,IE) = X(IC+NI+1)
               SKINY(3,IE) = Y(IC+NI+1)
               SKINZ(3,IE) = Z(IC+NI+1)
               SKINX(4,IE) = X(IC+NI)
               SKINY(4,IE) = Y(IC+NI)
               SKINZ(4,IE) = Z(IC+NI)
            ENDDO
         ENDDO
         IL = IL + (IDMN(IB)+1)*(JDMN(IB)+1)
         IM = IM +  IDMN(IB)*JDMN(IB)
      ENDDO

C ......................................................................


C ... INTEGRATE THE LOAD DISTRIBUTIONS


C ... Force distributions as function of x.

C ... Group the surface elements according to their locations in space.

      XMIN = MINVAL(SKINX)
      XMAX = MAXVAL(SKINX)
      YMIN = MINVAL(SKINY) 
      YMAX = MAXVAL(SKINY)
      ZMIN = MINVAL(SKINZ) 
      ZMAX = MAXVAL(SKINZ)

      DX = (XMAX-XMIN)/NSLOT+EPS
      DY = (YMAX-YMIN)+EPS
      DZ = (ZMAX-ZMIN)+EPS

      FLATX = DX < MIN(DY,DZ)/1000000.
      FLATY = DY < MIN(DX,DZ)/1000000.
      FLATZ = DZ < MIN(DX,DY)/1000000.

      IF(FLATX) THEN
         DXR = 1.0/1000000.
      ELSE
         DXR = 1.0/DX
      ENDIF

      IF(FLATY) THEN
         DYR = 1.0/1000000.
      ELSE
         DYR = 1.0/DY
      ENDIF

      IF(FLATZ) THEN
         DZR = 1.0/1000000.
      ELSE
         DZR = 1.0/DZ
      ENDIF

      ALLOCATE(KSC(1))
      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           LBOXES,NSLOT,1,1,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)
      DEALLOCATE(KSC)

      ALLOCATE(KSC(NCB),STAT=IERRCODE)

      IF(IERRCODE /= 0) THEN                                          
         WRITE(*,*)'LOADD:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           LBOXES,NSLOT,1,1,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)


      FXX = 0.0
      FYX = 0.0
      FZX = 0.0

      DO I1 = 1,NSLOT

C ... Load point location

         XC = XMIN + (I1-1)*DX + 0.5*DX

         FXX(I1,1)  = XC
         FYX(I1,1)  = XC
         FZX(I1,1)  = XC
         FXX(I1,2)  = 0.0
         FYX(I1,2)  = 0.0
         FZX(I1,2)  = 0.0
         FFXX(I1,2) = 0.0
         FFYX(I1,2) = 0.0
         FFZX(I1,2) = 0.0

C ... Load slot corner points (y constant plane)

         A1(1,1) = XMIN + (I1-1)*DX
         A1(2,1) = YMIN
         A1(3,1) = ZMIN
         A1(1,2) = XMIN + (I1-1)*DX
         A1(2,2) = YMIN
         A1(3,2) = ZMAX
         A1(1,3) = XMIN + I1*DX
         A1(2,3) = YMIN
         A1(3,3) = ZMAX
         A1(1,4) = XMIN + I1*DX
         A1(2,4) = YMIN
         A1(3,4) = ZMIN

C ... Load slot corner points (z constant plane)

         B1(1,1) = XMIN + (I1-1)*DX
         B1(2,1) = YMIN
         B1(3,1) = ZMIN
         B1(1,2) = XMIN + (I1-1)*DX
         B1(2,2) = YMAX
         B1(3,2) = ZMIN
         B1(1,3) = XMIN + I1*DX
         B1(2,3) = YMAX
         B1(3,3) = ZMIN
         B1(1,4) = XMIN + I1*DX
         B1(2,4) = YMIN
         B1(3,4) = ZMIN

         DO J1 = 1,1
            DO K1 = 1,1

               LBS = LBOXES(I1,J1,K1,2)                          
               LBE = LBOXES(I1,J1,K1,2) + LBOXES(I1,J1,K1,1) - 1

               DO LC=LBS,LBE

                  ISKIN = KSC(LC)

C ... Corner point coordinates of a surface element which lies at least 
C ... partly inside the load slot.

                  C1(1,1) = SKINX(1,ISKIN)
                  C1(2,1) = SKINY(1,ISKIN)
                  C1(3,1) = SKINZ(1,ISKIN)
                  C1(1,2) = SKINX(2,ISKIN)
                  C1(2,2) = SKINY(2,ISKIN)
                  C1(3,2) = SKINZ(2,ISKIN)
                  C1(1,3) = SKINX(3,ISKIN)
                  C1(2,3) = SKINY(3,ISKIN)
                  C1(3,3) = SKINZ(3,ISKIN)
                  C1(1,4) = SKINX(4,ISKIN)
                  C1(2,4) = SKINY(4,ISKIN)
                  C1(3,4) = SKINZ(4,ISKIN)

                  CALL CAREA(A1,C1,AELEM,APROJ,ASHARE)

                  ASHARE1 = ASHARE
                  APROJ1  = APROJ

                  CALL CAREA(B1,C1,AELEM,APROJ,ASHARE)

                  ASHARE2 = ASHARE
                  APROJ2  = APROJ

                  ASHARE = ASHARE1

                  IF(APROJ2 > APROJ1) ASHARE = ASHARE2

                  FXX(I1,2) = FXX(I1,2)+ASHARE*AELEM*FX(ISKIN)
                  FYX(I1,2) = FYX(I1,2)+ASHARE*AELEM*FY(ISKIN)
                  FZX(I1,2) = FZX(I1,2)+ASHARE*AELEM*FZ(ISKIN)

                  FFXX(I1,2) = FFXX(I1,2)+ASHARE*AELEM*FFX(ISKIN)
                  FFYX(I1,2) = FFYX(I1,2)+ASHARE*AELEM*FFY(ISKIN)
                  FFZX(I1,2) = FFZX(I1,2)+ASHARE*AELEM*FFZ(ISKIN)

               ENDDO

            ENDDO
         ENDDO

         FXX(I1,2) = FXX(I1,2)/DX
         FYX(I1,2) = FYX(I1,2)/DX
         FZX(I1,2) = FZX(I1,2)/DX

         FFXX(I1,2) = FFXX(I1,2)/DX
         FFYX(I1,2) = FFYX(I1,2)/DX
         FFZX(I1,2) = FFZX(I1,2)/DX

      ENDDO


      CX  = 0.0      
      CY  = 0.0      
      CZ  = 0.0
      CYM = 0.0
      CZM = 0.0
      MINFXX =  1.0E+30      
      MAXFXX = -1.0E+30      
      MINFYX =  1.0E+30      
      MAXFYX = -1.0E+30      
      MINFZX =  1.0E+30      
      MAXFZX = -1.0E+30      

      DO I1=1,NSLOT

         CX  = CX  + FXX(I1,2)*DX/QKINS
         CY  = CY  + FYX(I1,2)*DX/QKINS
         CZ  = CZ  + FZX(I1,2)*DX/QKINS
         CYM = CYM + (FYX(I1,1)-XMIN)*FYX(I1,2)*DX/QKINS
         CZM = CZM + (FZX(I1,1)-XMIN)*FZX(I1,2)*DX/QKINS
         IF(FXX(I1,2) < MINFXX) MINFXX = FXX(I1,2)
         IF(FXX(I1,2) > MAXFXX) MAXFXX = FXX(I1,2)
         IF(FYX(I1,2) < MINFYX) MINFYX = FYX(I1,2)
         IF(FYX(I1,2) > MAXFYX) MAXFYX = FYX(I1,2)
         IF(FZX(I1,2) < MINFZX) MINFZX = FZX(I1,2)
         IF(FZX(I1,2) > MAXFZX) MAXFZX = FZX(I1,2)

      ENDDO

      
      WRITE(400,*)
      WRITE(400,*) "Total loads integrated from next distributions ",
     &             "for checking purposes:"
      WRITE(400,*)
      IF(FLATX) THEN
         WRITE(400,*) "CY =      N/A       ",
     &               "XCP =      N/A        (load center)"     
         WRITE(400,*) "CZ =      N/A       ",
     &               "XCP =      N/A        (load center)" 
      ELSE
         WRITE(400,*) "CX =",REAL(CX,4)
         IF(CY == 0.0) THEN
         WRITE(400,*) "CY =",REAL(CY,4),
     &               "XCP =      N/A        (load center)"     
         ELSE
            WRITE(400,*) "CY =",REAL(CY,4),"XCP = ",
     &                          REAL(XMIN+CYM/CY,4), "(load center)" 
         ENDIF
         IF(CZ == 0.0) THEN
         WRITE(400,*) "CZ =",REAL(CZ,4),
     &               "XCP =      N/A        (load center)"     
         ELSE    
            WRITE(400,*) "CZ =",REAL(CZ,4),"XCP = ",
     &                          REAL(XMIN+CZM/CZ,4), "(load center)" 
         ENDIF
      ENDIF   
      WRITE(400,*)
      WRITE(400,*) "Minimum and maximum values for plotting routines:"
      WRITE(400,*)
      WRITE(400,*) "Min x  = ", REAL(FYX(1,1),4),
     &         "    Max x  = ", REAL(FYX(NSLOT,1),4) 
      IF(FLATX) THEN
         WRITE(400,*) "Min Fx =        N/A          Max Fx =        N/A" 
         WRITE(400,*) "Min Fy =        N/A          Max Fy =        N/A" 
         WRITE(400,*) "Min Fz =        N/A          Max Fz =        N/A"
      ELSE
         WRITE(400,*) "Min Fx = ", REAL(MINFXX,4),
     &            "    Max Fx = ", REAL(MAXFXX,4) 
         WRITE(400,*) "Min Fy = ", REAL(MINFYX,4),
     &            "    Max Fy = ", REAL(MAXFYX,4) 
         WRITE(400,*) "Min Fz = ", REAL(MINFZX,4),
     &            "    Max Fz = ", REAL(MAXFZX,4)
      ENDIF 
      WRITE(400,*)
      WRITE(400,*)
     &     "          i     x[m]  ",                                     
     &     "            Fx(x)[N/m]        Fy(x)[N/m]        Fz(x)[N/m]",
     &     "        sumFx(x)[N]       sumFy(x)[N]       sumFz(x)[N]" ,
     &     "       fx(x)[N/m]        fy(x)[N/m]        fz(x)[N/m]",
     &     "        sumfx(x)[N]       sumfy(x)[N]       sumfz(x)[N]"
      WRITE(400,*)

      DO I1=1,NSLOT
         IF(FLATX) THEN
            WRITE(400,'(I12,1E18.8,216A)') I1, REAL(FYX(I1,1),4), 
     &                 '      N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A'
         ELSE
            IF(I1 > 1) THEN
               SUMX = SUMX + DX*FXX(I1,2) 
               SUMY = SUMY + DX*FYX(I1,2)
               SUMZ = SUMZ + DX*FZX(I1,2)
               FRICSUMX = FRICSUMX + DX*FFXX(I1,2) 
               FRICSUMY = FRICSUMY + DX*FFYX(I1,2)
               FRICSUMZ = FRICSUMZ + DX*FFZX(I1,2)
            ELSE
               SUMX = DX*FXX(I1,2)
               SUMY = DX*FYX(I1,2)
               SUMZ = DX*FZX(I1,2)
               FRICSUMX = DX*FFXX(I1,2)
               FRICSUMY = DX*FFYX(I1,2)
               FRICSUMZ = DX*FFZX(I1,2)
            ENDIF
            WRITE(400,'(I12,13E18.8)')        I1, REAL(FYX(I1,1),4), 
     &      REAL(FXX(I1,2),4), REAL(FYX(I1,2),4), REAL(FZX(I1,2),4),
     &      REAL(SUMX,4),      REAL(SUMY,4),      REAL(SUMZ,4),
     &      REAL(FFXX(I1,2),4),REAL(FFYX(I1,2),4),REAL(FFZX(I1,2),4),
     &      REAL(FRICSUMX,4),  REAL(FRICSUMY,4),  REAL(FRICSUMZ,4)
         ENDIF
      ENDDO

      DEALLOCATE(KSC)



C ... Force distributions as function of y.

C ... Group the surface elements according to their locations in space.

      XMIN = MINVAL(SKINX)
      XMAX = MAXVAL(SKINX)
      YMIN = MINVAL(SKINY) 
      YMAX = MAXVAL(SKINY)
      ZMIN = MINVAL(SKINZ) 
      ZMAX = MAXVAL(SKINZ)

      DX = (XMAX-XMIN)+EPS
      DY = (YMAX-YMIN)/NSLOT+EPS
      DZ = (ZMAX-ZMIN)+EPS

      FLATX = DX < MIN(DY,DZ)/1000000.
      FLATY = DY < MIN(DX,DZ)/1000000.
      FLATZ = DZ < MIN(DX,DY)/1000000.

      IF(FLATX) THEN
         DXR = 1.0/1000000.
      ELSE
         DXR = 1.0/DX
      ENDIF

      IF(FLATY) THEN
         DYR = 1.0/1000000.
      ELSE
         DYR = 1.0/DY
      ENDIF

      IF(FLATZ) THEN
         DZR = 1.0/1000000.
      ELSE
         DZR = 1.0/DZ
      ENDIF

      ALLOCATE(KSC(1))
      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           MBOXES,1,NSLOT,1,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)
      DEALLOCATE(KSC)

      ALLOCATE(KSC(NCB),STAT=IERRCODE)

      IF(IERRCODE /= 0) THEN                                          
         WRITE(*,*)'LOADD:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           MBOXES,1,NSLOT,1,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)


      FXY = 0.0
      FYY = 0.0
      FZY = 0.0

      DO J1 = 1,NSLOT

C ... Load point location

         YC = YMIN + (J1-1)*DY + 0.5*DY

         FXY(J1,1)  = YC
         FYY(J1,1)  = YC
         FZY(J1,1)  = YC
         FXY(J1,2)  = 0.0
         FYY(J1,2)  = 0.0
         FZY(J1,2)  = 0.0
         FFXY(J1,2) = 0.0
         FFYY(J1,2) = 0.0
         FFZY(J1,2) = 0.0

C ... Load slot corner points (x constant plane)

         A1(1,1) = XMIN
         A1(2,1) = YMIN + (J1-1)*DY
         A1(3,1) = ZMIN
         A1(1,2) = XMIN
         A1(2,2) = YMIN + (J1-1)*DY
         A1(3,2) = ZMAX
         A1(1,3) = XMIN
         A1(2,3) = YMIN + J1*DY
         A1(3,3) = ZMAX
         A1(1,4) = XMIN
         A1(2,4) = YMIN + J1*DY
         A1(3,4) = ZMIN

C ... Load slot corner points (z constant plane)

         B1(1,1) = XMIN
         B1(2,1) = YMIN + (J1-1)*DY
         B1(3,1) = ZMIN
         B1(1,2) = XMAX
         B1(2,2) = YMIN + (J1-1)*DY
         B1(3,2) = ZMIN
         B1(1,3) = XMAX
         B1(2,3) = YMIN + J1*DY
         B1(3,3) = ZMIN
         B1(1,4) = XMIN
         B1(2,4) = YMIN + J1*DY
         B1(3,4) = ZMIN

         DO I1 = 1,1
            DO K1 = 1,1

               LBS = MBOXES(I1,J1,K1,2)                          
               LBE = MBOXES(I1,J1,K1,2) + MBOXES(I1,J1,K1,1) - 1

               DO LC=LBS,LBE

                  ISKIN = KSC(LC)

C ... Corner point coordinates of a surface element which lies at least 
C ... partly inside the load slot.

                  C1(1,1) = SKINX(1,ISKIN)
                  C1(2,1) = SKINY(1,ISKIN)
                  C1(3,1) = SKINZ(1,ISKIN)
                  C1(1,2) = SKINX(2,ISKIN)
                  C1(2,2) = SKINY(2,ISKIN)
                  C1(3,2) = SKINZ(2,ISKIN)
                  C1(1,3) = SKINX(3,ISKIN)
                  C1(2,3) = SKINY(3,ISKIN)
                  C1(3,3) = SKINZ(3,ISKIN)
                  C1(1,4) = SKINX(4,ISKIN)
                  C1(2,4) = SKINY(4,ISKIN)
                  C1(3,4) = SKINZ(4,ISKIN)

                  CALL CAREA(A1,C1,AELEM,APROJ,ASHARE)

                  ASHARE1 = ASHARE
                  APROJ1  = APROJ

                  CALL CAREA(B1,C1,AELEM,APROJ,ASHARE)

                  ASHARE2 = ASHARE
                  APROJ2  = APROJ

                  ASHARE = ASHARE1

                  IF(APROJ2 > APROJ1) ASHARE = ASHARE2

                  FXY(J1,2) = FXY(J1,2)+ASHARE*AELEM*FX(ISKIN)
                  FYY(J1,2) = FYY(J1,2)+ASHARE*AELEM*FY(ISKIN)
                  FZY(J1,2) = FZY(J1,2)+ASHARE*AELEM*FZ(ISKIN)

                  FFXY(J1,2) = FFXY(J1,2)+ASHARE*AELEM*FFX(ISKIN)
                  FFYY(J1,2) = FFYY(J1,2)+ASHARE*AELEM*FFY(ISKIN)
                  FFZY(J1,2) = FFZY(J1,2)+ASHARE*AELEM*FFZ(ISKIN)


               ENDDO

            ENDDO
         ENDDO
         
         FXY(J1,2) = FXY(J1,2)/DY
         FYY(J1,2) = FYY(J1,2)/DY
         FZY(J1,2) = FZY(J1,2)/DY
         
         FFXY(J1,2) = FFXY(J1,2)/DY
         FFYY(J1,2) = FFYY(J1,2)/DY
         FFZY(J1,2) = FFZY(J1,2)/DY

      ENDDO


      CX  = 0.0      
      CY  = 0.0      
      CZ  = 0.0
      CXM = 0.0
      CZM = 0.0      
      MINFXY =  1.0E+30      
      MAXFXY = -1.0E+30      
      MINFYY =  1.0E+30      
      MAXFYY = -1.0E+30      
      MINFZY =  1.0E+30      
      MAXFZY = -1.0E+30      

      DO J1=1,NSLOT

         CX  = CX  + FXY(J1,2)*DY/QKINS
         CY  = CY  + FYY(J1,2)*DY/QKINS
         CZ  = CZ  + FZY(J1,2)*DY/QKINS
         CXM = CXM + (FXY(J1,1)-YMIN)*FXY(J1,2)*DY/QKINS
         CZM = CZM + (FZY(J1,1)-YMIN)*FZY(J1,2)*DY/QKINS
         IF(FXY(J1,2) < MINFXY) MINFXY = FXY(J1,2)
         IF(FXY(J1,2) > MAXFXY) MAXFXY = FXY(J1,2)
         IF(FYY(J1,2) < MINFYY) MINFYY = FYY(J1,2)
         IF(FYY(J1,2) > MAXFYY) MAXFYY = FYY(J1,2)
         IF(FZY(J1,2) < MINFZY) MINFZY = FZY(J1,2)
         IF(FZY(J1,2) > MAXFZY) MAXFZY = FZY(J1,2)

      ENDDO


      WRITE(400,*)
      WRITE(400,*)
      WRITE(400,*) "Total loads integrated from next distributions ",
     &             "for checking purposes:"
      WRITE(400,*)
      IF(FLATY) THEN
         WRITE(400,*) "CX =      N/A       ",
     &               "YCP =      N/A        (load center)"     
         WRITE(400,*) "CZ =      N/A       ",
     &               "YCP =      N/A        (load center)" 
      ELSE
         IF(CX == 0.0) THEN
         WRITE(400,*) "CX =", REAL(CX,4),
     &               "YCP =      N/A        (load center)"     
         ELSE
            WRITE(400,*) "CX =", REAL(CX,4),
     &                  "YCP = ",REAL(YMIN+CXM/CX,4), "(load center)" 
         ENDIF
         WRITE(400,*) "CY =", REAL(CY,4)
         IF(CZ == 0.0) THEN
         WRITE(400,*) "CZ =", REAL(CZ,4),
     &               "YCP =      N/A        (load center)" 
         ELSE   
            WRITE(400,*) "CZ =", REAL(CZ,4),
     &                  "YCP = ", REAL(YMIN+CZM/CZ,4), "(load center)" 
         ENDIF
      ENDIF   
      WRITE(400,*)
      WRITE(400,*) "Minimum and maximum values for plotting routines:"
      WRITE(400,*)
      WRITE(400,*) "Min y  = ", REAL(FXY(1,1),4),
     &         "    Max y  = ", REAL(FXY(NSLOT,1),4) 
      IF(FLATY) THEN
         WRITE(400,*) "Min Fx =        N/A          Max Fx =        N/A" 
         WRITE(400,*) "Min Fy =        N/A          Max Fy =        N/A"
         WRITE(400,*) "Min Fz =        N/A          Max Fz =        N/A"
      ELSE
         WRITE(400,*) "Min Fx = ", REAL(MINFXY,4),
     &            "    Max Fx = ", REAL(MAXFXY,4) 
         WRITE(400,*) "Min Fy = ", REAL(MINFYY,4),
     &            "    Max Fy = ", REAL(MAXFYY,4) 
         WRITE(400,*) "Min Fz = ", REAL(MINFZY,4),
     &            "    Max Fz = ", REAL(MAXFZY,4) 
      ENDIF
      WRITE(400,*)
      WRITE(400,*)
     &     "          i     y[m]  ",
     &     "            Fx(y)[N/m]        Fy(y)[N/m]        Fz(y)[N/m]",
     &     "        sumFx(y)[N]       sumFy(y)[N]       sumFz(y)[N]" ,
     &     "       fx(y)[N/m]        fy(y)[N/m]        fz(y)[N/m]",
     &     "        sumfx(y)[N]       sumfy(y)[N]       sumfz(y)[N]"
      WRITE(400,*)

      DO J1=1,NSLOT
         IF(FLATY) THEN
            WRITE(400,'(I12,1E18.8,216A)') I1, REAL(FZY(I1,1),4), 
     &                 '      N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A'
         ELSE
            IF(J1 > 1) THEN
               SUMX = SUMX + DY*FXY(J1,2) 
               SUMY = SUMY + DY*FYY(J1,2)
               SUMZ = SUMZ + DY*FZY(J1,2)
               FRICSUMX = FRICSUMX + DY*FFXY(J1,2) 
               FRICSUMY = FRICSUMY + DY*FFYY(J1,2)
               FRICSUMZ = FRICSUMZ + DY*FFZY(J1,2)
            ELSE
               SUMX = DY*FXY(J1,2)
               SUMY = DY*FYY(J1,2)
               SUMZ = DY*FZY(J1,2)
               FRICSUMX = DY*FFXY(J1,2)
               FRICSUMY = DY*FFYY(J1,2)
               FRICSUMZ = DY*FFZY(J1,2)
            ENDIF
            WRITE(400,'(I12,13E18.8)')        J1, REAL(FZY(J1,1),4), 
     &      REAL(FXY(J1,2),4), REAL(FYY(J1,2),4), REAL(FZY(J1,2),4),
     &      REAL(SUMX,4),      REAL(SUMY,4),      REAL(SUMZ,4),
     &      REAL(FFXY(J1,2),4),REAL(FFYY(J1,2),4),REAL(FFZY(J1,2),4),
     &      REAL(FRICSUMX,4),  REAL(FRICSUMY,4),  REAL(FRICSUMZ,4)
         ENDIF
      ENDDO

      DEALLOCATE(KSC)



C ... Force distributions as function of z.

C ... Group the surface elements according to their locations in space.

      XMIN = MINVAL(SKINX)
      XMAX = MAXVAL(SKINX)
      YMIN = MINVAL(SKINY) 
      YMAX = MAXVAL(SKINY)
      ZMIN = MINVAL(SKINZ) 
      ZMAX = MAXVAL(SKINZ)

      DX = (XMAX-XMIN)+EPS
      DY = (YMAX-YMIN)+EPS
      DZ = (ZMAX-ZMIN)/NSLOT+EPS

      FLATX = DX < MIN(DY,DZ)/1000000.
      FLATY = DY < MIN(DX,DZ)/1000000.
      FLATZ = DZ < MIN(DX,DY)/1000000.

      IF(FLATX) THEN
         DXR = 1.0/1000000.
      ELSE
         DXR = 1.0/DX
      ENDIF

      IF(FLATY) THEN
         DYR = 1.0/1000000.
      ELSE
         DYR = 1.0/DY
      ENDIF

      IF(FLATZ) THEN
         DZR = 1.0/1000000.
      ELSE
         DZR = 1.0/DZ
      ENDIF

      ALLOCATE(KSC(1))
      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           NBOXES,1,1,NSLOT,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)
      DEALLOCATE(KSC)

      ALLOCATE(KSC(NCB),STAT=IERRCODE)

      IF(IERRCODE /= 0) THEN                                          
         WRITE(*,*)'LOADD:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL(NOFEL,SKINX,SKINY,SKINZ,KSC,NCB,
     &           NBOXES,1,1,NSLOT,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)


      FXZ = 0.0
      FYZ = 0.0
      FZZ = 0.0

      DO K1 = 1,NSLOT

C ... Load point location

         ZC = ZMIN + (K1-1)*DZ + 0.5*DZ

         FXZ(K1,1)  = ZC
         FYZ(K1,1)  = ZC
         FZZ(K1,1)  = ZC
         FXZ(K1,2)  = 0.0
         FYZ(K1,2)  = 0.0
         FZZ(K1,2)  = 0.0
         FFXZ(K1,2) = 0.0
         FFYZ(K1,2) = 0.0
         FFZZ(K1,2) = 0.0

C ... Load slot corner points (x constant plane)

         A1(1,1) = XMIN
         A1(2,1) = YMIN
         A1(3,1) = ZMIN + (K1-1)*DZ
         A1(1,2) = XMIN
         A1(2,2) = YMAX
         A1(3,2) = ZMIN + (K1-1)*DZ
         A1(1,3) = XMIN
         A1(2,3) = YMAX
         A1(3,3) = ZMIN + K1*DZ
         A1(1,4) = XMIN
         A1(2,4) = YMIN
         A1(3,4) = ZMIN + K1*DZ

C ... Load slot corner points (y constant plane)

         B1(1,1) = XMIN
         B1(2,1) = YMIN
         B1(3,1) = ZMIN + (K1-1)*DZ
         B1(1,2) = XMAX
         B1(2,2) = YMIN
         B1(3,2) = ZMIN + (K1-1)*DZ
         B1(1,3) = XMAX
         B1(2,3) = YMIN
         B1(3,3) = ZMIN + K1*DZ
         B1(1,4) = XMIN
         B1(2,4) = YMIN
         B1(3,4) = ZMIN + K1*DZ

         DO I1 = 1,1
            DO J1 = 1,1

               LBS = NBOXES(I1,J1,K1,2)                          
               LBE = NBOXES(I1,J1,K1,2) + NBOXES(I1,J1,K1,1) - 1

               DO LC=LBS,LBE

                  ISKIN = KSC(LC)

C ... Corner point coordinates of a surface element which lies at least 
C ... partly inside the load slot.

                  C1(1,1) = SKINX(1,ISKIN)
                  C1(2,1) = SKINY(1,ISKIN)
                  C1(3,1) = SKINZ(1,ISKIN)
                  C1(1,2) = SKINX(2,ISKIN)
                  C1(2,2) = SKINY(2,ISKIN)
                  C1(3,2) = SKINZ(2,ISKIN)
                  C1(1,3) = SKINX(3,ISKIN)
                  C1(2,3) = SKINY(3,ISKIN)
                  C1(3,3) = SKINZ(3,ISKIN)
                  C1(1,4) = SKINX(4,ISKIN)
                  C1(2,4) = SKINY(4,ISKIN)
                  C1(3,4) = SKINZ(4,ISKIN)

                  CALL CAREA(A1,C1,AELEM,APROJ,ASHARE)

                  ASHARE1 = ASHARE
                  APROJ1  = APROJ

                  CALL CAREA(B1,C1,AELEM,APROJ,ASHARE)

                  ASHARE2 = ASHARE
                  APROJ2  = APROJ

                  ASHARE = ASHARE1

                  IF(APROJ2 > APROJ1) ASHARE = ASHARE2

                  FXZ(K1,2) = FXZ(K1,2)+ASHARE*AELEM*FX(ISKIN)
                  FYZ(K1,2) = FYZ(K1,2)+ASHARE*AELEM*FY(ISKIN)
                  FZZ(K1,2) = FZZ(K1,2)+ASHARE*AELEM*FZ(ISKIN)

                  FFXZ(K1,2) = FFXZ(K1,2)+ASHARE*AELEM*FFX(ISKIN)
                  FFYZ(K1,2) = FFYZ(K1,2)+ASHARE*AELEM*FFY(ISKIN)
                  FFZZ(K1,2) = FFZZ(K1,2)+ASHARE*AELEM*FFZ(ISKIN)

               ENDDO

            ENDDO
         ENDDO

         FXZ(K1,2) = FXZ(K1,2)/DZ
         FYZ(K1,2) = FYZ(K1,2)/DZ
         FZZ(K1,2) = FZZ(K1,2)/DZ

         FFXZ(K1,2) = FFXZ(K1,2)/DZ
         FFYZ(K1,2) = FFYZ(K1,2)/DZ
         FFZZ(K1,2) = FFZZ(K1,2)/DZ

      ENDDO


      CX  = 0.0      
      CY  = 0.0
      CZ  = 0.0
      CXM = 0.0
      CYM = 0.0      
      MINFXZ =  1.0E+30      
      MAXFXZ = -1.0E+30      
      MINFYZ =  1.0E+30      
      MAXFYZ = -1.0E+30      
      MINFZZ =  1.0E+30      
      MAXFZZ = -1.0E+30      

      DO K1=1,NSLOT

         CX  = CX  + FXZ(K1,2)*DZ/QKINS
         CY  = CY  + FYZ(K1,2)*DZ/QKINS
         CZ  = CZ  + FZZ(K1,2)*DZ/QKINS
         CXM = CXM + (FXZ(K1,1)-ZMIN)*FXZ(K1,2)*DZ/QKINS
         CYM = CYM + (FYZ(K1,1)-ZMIN)*FYZ(K1,2)*DZ/QKINS
         IF(FXZ(K1,2) < MINFXZ) MINFXZ = FXZ(K1,2)
         IF(FXZ(K1,2) > MAXFXZ) MAXFXZ = FXZ(K1,2)
         IF(FYZ(K1,2) < MINFYZ) MINFYZ = FYZ(K1,2)
         IF(FYZ(K1,2) > MAXFYZ) MAXFYZ = FYZ(K1,2)
         IF(FZZ(K1,2) < MINFZZ) MINFZZ = FZZ(K1,2)
         IF(FZZ(K1,2) > MAXFZZ) MAXFZZ = FZZ(K1,2)

      ENDDO


      WRITE(400,*)
      WRITE(400,*)
      WRITE(400,*) "Total loads integrated from next distributions ",
     &             "for checking purposes:"
      WRITE(400,*)
      IF(FLATZ) THEN
         WRITE(400,*) "CX =      N/A       ",
     &               "ZCP =      N/A        (load center)"     
         WRITE(400,*) "CY =      N/A       ",
     &               "ZCP =      N/A        (load center)" 
      ELSE
         IF(CX == 0.0) THEN
         WRITE(400,*) "CX =", REAL(CX,4),   
     &               "ZCP =      N/A        (load center)"     
         ELSE
            WRITE(400,*) "CX =", REAL(CX,4),
     &                  "ZCP = ",REAL(ZMIN+CXM/CX,4), "(load center)" 
         ENDIF
         IF(CY == 0.0) THEN
         WRITE(400,*) "CY =", REAL(CY,4),
     &               "ZCP =      N/A        (load center)"     
         ELSE    
            WRITE(400,*) "CY =", REAL(CY,4),
     &                  "ZCP = ",REAL(ZMIN+CYM/CY,4), "(load center)" 
         ENDIF
         WRITE(400,*) "CZ =", REAL(CZ,4)
      ENDIF    
      WRITE(400,*)
      WRITE(400,*) "Minimum and maximum values for plotting routines:"
      WRITE(400,*)
      WRITE(400,*) "Min z  = ", REAL(FXZ(1,1),4),
     &         "    Max z  = ", REAL(FXZ(NSLOT,1),4)
      IF(FLATZ) THEN 
         WRITE(400,*) "Min Fx =        N/A          Max Fx =        N/A" 
         WRITE(400,*) "Min Fy =        N/A          Max Fy =        N/A"
         WRITE(400,*) "Min Fz =        N/A          Max Fz =        N/A"
      ELSE
         WRITE(400,*) "Min Fx = ", REAL(MINFXZ,4),
     &            "    Max Fx = ", REAL(MAXFXZ,4) 
         WRITE(400,*) "Min Fy = ", REAL(MINFYZ,4),
     &            "    Max Fy = ", REAL(MAXFYZ,4) 
         WRITE(400,*) "Min Fz = ", REAL(MINFZZ,4),
     &            "    Max Fz = ", REAL(MAXFZZ,4) 
      ENDIF
      WRITE(400,*)
      WRITE(400,*)
     &     "          i     z[m]  ",
     &     "            Fx(z)[N/m]        Fy(z)[N/m]        Fz(z)[N/m]",
     &     "        sumFx(z)[N]       sumFy(z)[N]       sumFz(z)[N]" ,
     &     "       fx(z)[N/m]        fy(z)[N/m]        fz(z)[N/m]",
     &     "        sumfx(z)[N]       sumfy(z)[N]       sumfz(z)[N]"
      WRITE(400,*)

      DO K1=1,NSLOT
         IF(FLATZ) THEN
            WRITE(400,'(I12,1E18.8,216A)') I1, REAL(FXZ(I1,1),4), 
     &                 '      N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A',
     &        '               N/A               N/A               N/A'
         ELSE
            IF(K1 > 1) THEN
               SUMX = SUMX + DZ*FXZ(K1,2) 
               SUMY = SUMY + DZ*FYZ(K1,2)
               SUMZ = SUMZ + DZ*FZZ(K1,2)
               FRICSUMX = FRICSUMX + DZ*FFXZ(K1,2) 
               FRICSUMY = FRICSUMY + DZ*FFYZ(K1,2)
               FRICSUMZ = FRICSUMZ + DZ*FFZZ(K1,2)
            ELSE
               SUMX = DZ*FXZ(K1,2)
               SUMY = DZ*FYZ(K1,2)
               SUMZ = DZ*FZZ(K1,2)
               FRICSUMX = DZ*FFXZ(K1,2)
               FRICSUMY = DZ*FFYZ(K1,2)
               FRICSUMZ = DZ*FFZZ(K1,2)
            ENDIF
            WRITE(400,'(I12,13E18.8)')        K1, REAL(FXZ(K1,1),4), 
     &      REAL(FXZ(K1,2),4), REAL(FYZ(K1,2),4), REAL(FZZ(K1,2),4),
     &      REAL(SUMX,4),      REAL(SUMY,4),      REAL(SUMZ,4),
     &      REAL(FFXZ(K1,2),4),REAL(FFYZ(K1,2),4),REAL(FFZZ(K1,2),4),
     &      REAL(FRICSUMX,4),  REAL(FRICSUMY,4),  REAL(FRICSUMZ,4)
         ENDIF
      ENDDO

      DEALLOCATE(KSC)
      DEALLOCATE(SKINX,SKINY,SKINZ)

      RETURN                                                              
      END SUBROUTINE LOADD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CAREA(A,E,AELEM,APROJ,ASHARE)      

C ... This subroutine computes the common area of two overlaping
C ... polygons. The second of the polygons (E, a surface element)
C ... will be projected onto the first one (A, load slot plane). 
C ... APROJ is the projected area and ASHARE is the share of the
C ... projected area which lies inside the slot (first polygon). 

      IMPLICIT REAL (A-H,O-Z)

      REAL :: A(3,5),B(3,5),C(4,25),D(3,4),E(3,5),F(3,3),G(4),
     +        XCD,YCD,ZCD,D1,D2,D3,D4

      REAL :: APROJ, ASHARE, XCA,YCA,ZCA,XCB,YCB,ZCB,RA,RB

      LOGICAL :: SORTED
     
      EPS = 1.0E-10


C ... Normal of polygon A 

      CX = A(1,3) - A(1,1)
      CY = A(2,3) - A(2,1)
      CZ = A(3,3) - A(3,1)

      DX = A(1,4) - A(1,2)
      DY = A(2,4) - A(2,2)
      DZ = A(3,4) - A(3,2)

      ANX = CY*DZ - CZ*DY
      ANY = CZ*DX - CX*DZ
      ANZ = CX*DY - CY*DX
      AA  = SQRT(ANX**2 + ANY**2 + ANZ**2)
      AA  = MAX(AA,EPS)
      ANX = ANX/AA
      ANY = ANY/AA
      ANZ = ANZ/AA


C ... Surface elements area

      AX = E(1,3) - E(1,1)
      AY = E(2,3) - E(2,1)
      AZ = E(3,3) - E(3,1)

      BX = E(1,4) - E(1,2)
      BY = E(2,4) - E(2,2)
      BZ = E(3,4) - E(3,2)

      AELEM = 0.5*SQRT((AY*BZ-AZ*BY)**2
     &               + (AZ*BX-AX*BZ)**2
     &               + (AX*BY-AY*BX)**2)

C ... Project the corner points of the surface element onto the
C ... load slot plane

      DO K = 1,4
         S = ANX*(A(1,1)-E(1,K))+ANY*(A(2,1)-E(2,K))+ANZ*(A(3,1)-E(3,K))
         B(1,K) = E(1,K) + S*ANX
         B(2,K) = E(2,K) + S*ANY
         B(3,K) = E(3,K) + S*ANZ
      ENDDO


C ... Area of the projection

      AX = B(1,3) - B(1,1)
      AY = B(2,3) - B(2,1)
      AZ = B(3,3) - B(3,1)

      BX = B(1,4) - B(1,2)
      BY = B(2,4) - B(2,2)
      BZ = B(3,4) - B(3,2)

      APROJ = 0.5*SQRT((AY*BZ-AZ*BY)**2
     &               + (AZ*BX-AX*BZ)**2
     &               + (AX*BY-AY*BX)**2)

C ... Close the polygons

      A(1,5) = A(1,1)
      A(2,5) = A(2,1)
      A(3,5) = A(3,1)

      B(1,5) = B(1,1)
      B(2,5) = B(2,1)
      B(3,5) = B(3,1)


      NPOINT = 0

C ... Check first if the corner(s) of polygon A lie inside
C ... polygon B and vice versa.

      DO N=1,4
         DO K=1,4 
            L = K + 1 
            AX = B(1,L) - B(1,K)
            AY = B(2,L) - B(2,K)
            AZ = B(3,L) - B(3,K)
            BX = A(1,N) - B(1,K)
            BY = A(2,N) - B(2,K)
            BZ = A(3,N) - B(3,K)
            SX = AY*BZ  - AZ*BY
            SY = AZ*BX  - AX*BZ
            SZ = AX*BY  - AY*BX
            G(K) = ANX*SX + ANY*SY + ANZ*SZ
         ENDDO
         GMIN = MIN(G(1),G(2),G(3),G(4))
         GMAX = MAX(G(1),G(2),G(3),G(4))
         IF(GMIN*GMAX >= 0.0) THEN
            NPOINT = NPOINT + 1
            C(1,NPOINT) = A(1,N)
            C(2,NPOINT) = A(2,N) 
            C(3,NPOINT) = A(3,N)
         ENDIF 
      ENDDO

      DO N=1,4
         DO K=1,4 
            L = K + 1 
            AX = A(1,L) - A(1,K)
            AY = A(2,L) - A(2,K)
            AZ = A(3,L) - A(3,K)
            BX = B(1,N) - A(1,K)
            BY = B(2,N) - A(2,K)
            BZ = B(3,N) - A(3,K)
            SX = AY*BZ  - AZ*BY
            SY = AZ*BX  - AX*BZ
            SZ = AX*BY  - AY*BX
            G(K) = ANX*SX + ANY*SY + ANZ*SZ
         ENDDO
         GMIN = MIN(G(1),G(2),G(3),G(4))
         GMAX = MAX(G(1),G(2),G(3),G(4))
         IF(GMIN*GMAX >= 0.0) THEN
          NPOINT = NPOINT + 1
            C(1,NPOINT) = B(1,N)
            C(2,NPOINT) = B(2,N) 
            C(3,NPOINT) = B(3,N)
        ENDIF 
      ENDDO



C ... Search next locations where the polygon edges intersect each other.

      DO N=1,4
        M = N + 1
         DO L=1,4
          K = L + 1

            AX = A(1,L) - A(1,K)                
            AY = A(2,L) - A(2,K)                
            AZ = A(3,L) - A(3,K)                

            BX = B(1,M) - B(1,N)                
            BY = B(2,M) - B(2,N)                
            BZ = B(3,M) - B(3,N)                

            DX = B(1,M) - A(1,K)                
            DY = B(2,M) - A(2,K)
            DZ = B(3,M) - A(3,K)

            UX = DY*BZ  - DZ*BY
            UY = DZ*BX  - DX*BZ
            UZ = DX*BY  - DY*BX
            UL = SQRT(UX**2 + UY**2 + UZ**2)

            VX = AY*BZ  - AZ*BY
            VY = AZ*BX  - AX*BZ
            VZ = AX*BY  - AY*BX
            VL = SQRT(VX**2 + VY**2 + VZ**2)

            PX = B(1,N) - A(1,K) 
            PY = B(2,N) - A(2,K) 
            PZ = B(3,N) - A(3,K) 

            QX = A(1,L) - B(1,N) 
            QY = A(2,L) - B(2,N) 
            QZ = A(3,L) - B(3,N) 

            RX = B(1,M) - A(1,L) 
            RY = B(2,M) - A(2,L) 
            RZ = B(3,M) - A(3,L) 

            SX = A(1,K) - B(1,M) 
            SY = A(2,K) - B(2,M) 
            SZ = A(3,K) - B(3,M) 

            PRQX = PY*QZ - PZ*QY
            PRQY = PZ*QX - PX*QZ
            PRQZ = PX*QY - PY*QX

            QRRX = QY*RZ - QZ*RY
            QRRY = QZ*RX - QX*RZ
            QRRZ = QX*RY - QY*RX

            RRSX = RY*SZ - RZ*SY
            RRSY = RZ*SX - RX*SZ
            RRSZ = RX*SY - RY*SX

            SRPX = SY*PZ - SZ*PY
            SRPY = SZ*PX - SX*PZ
            SRPZ = SX*PY - SY*PX

            G(1) = ANX*PRQX + ANY*PRQY + ANZ*PRQZ           
            G(2) = ANX*QRRX + ANY*QRRY + ANZ*QRRZ           
            G(3) = ANX*RRSX + ANY*RRSY + ANZ*RRSZ           
            G(4) = ANX*SRPX + ANY*SRPY + ANZ*SRPZ           

            GMIN = MIN(G(1),G(2),G(3),G(4))
            GMAX = MAX(G(1),G(2),G(3),G(4))

            IF(GMIN*GMAX > 0.0) THEN
               IF(VL /= 0.0) THEN
                  S = UL/VL
                  IF(S > 0.0 .AND. S < 1.0) THEN
                     NPOINT = NPOINT + 1
                     C(1,NPOINT) = A(1,K) + S*AX 
                     C(2,NPOINT) = A(2,K) + S*AY 
                     C(3,NPOINT) = A(3,K) + S*AZ 
                  ENDIF
               ENDIF
            ENDIF                
         ENDDO
      ENDDO


C ... Search the points which were found twice or even more often

      DO IP=1,NPOINT-1
         DO JP=IP+1,NPOINT
            DD = SQRT((C(1,JP)-C(1,IP))**2 + (C(2,JP)-C(2,IP))**2
     &                                      + (C(3,JP)-C(3,IP))**2)
            IF(DD < EPS) THEN
               C(1,JP) = 99999.0 
               C(2,JP) = 99999.0
               C(3,JP) = 99999.0
          ENDIF 
        ENDDO
      ENDDO

      KP = 0
      DO IP=1,NPOINT
         IF(C(1,IP) < 99998.0) THEN 
            KP = KP + 1
            C(1,KP) = C(1,IP)
            C(2,KP) = C(2,IP)
            C(3,KP) = C(3,IP)
         ENDIF
      ENDDO

      NPOINT = KP

      IF(NPOINT < 3) THEN
         ASHARE = 0.0
         RETURN
      ENDIF


C ... Compute the center of the overlaping area.

      XCP = 0.0
      YCP = 0.0
      ZCP = 0.0
      DO NP = 1,NPOINT
         XCP = XCP + C(1,NP)
         YCP = YCP + C(2,NP)
         ZCP = ZCP + C(3,NP)
      ENDDO
      XCP = XCP/NPOINT
      YCP = YCP/NPOINT 
      ZCP = ZCP/NPOINT


C ... Reorder the polygon corner points

      CX = C(1,1) - XCP  ! Vector from center to point 1 (vector a)
      CY = C(2,1) - YCP
      CZ = C(3,1) - ZCP
      CC = SQRT(CX**2 + CY**2 + CZ**2)
      C(4,1) = 0.0

      DO I = 2,NPOINT
         DX = C(1,I) - XCP ! Vector from center to point i (vector b)
         DY = C(2,I) - YCP
         DZ = C(3,I) - ZCP
         DD = SQRT(DX**2 + DY**2 + DZ**2)
         HX = CY*DZ-CZ*DY
         HY = CZ*DX-CX*DZ
         HZ = CX*DY-CY*DX
         HH = ANX*HX + ANY*HY + ANZ*HZ
         AA = MAX((CX*DX+CY*DY+CZ*DZ)/(CC*DD), -1.0)
         AA = MIN(AA, 1.0)
         C(4,I) = SIGN(ACOS(AA),HH)  ! Angle between vectors a and b
      ENDDO


C ... Bubble sort (use angle as the sorting key)

      I = 1
 10   CONTINUE
         I = I + 1
         SORTED = .TRUE.
         DO J=NPOINT,I,-1
            IF(C(4,J-1) > C(4,J)) THEN
               SORTED   = .FALSE.
               XAPU     = C(1,J-1)
               YAPU     = C(2,J-1)
               ZAPU     = C(3,J-1)
               AAPU     = C(4,J-1)
               C(1,J-1) = C(1,J)
               C(2,J-1) = C(2,J)
               C(3,J-1) = C(3,J)
               C(4,J-1) = C(4,J)
               C(1,J)   = XAPU 
               C(2,J)   = YAPU 
               C(3,J)   = ZAPU 
               C(4,J)   = AAPU 
            ENDIF
         ENDDO
         IF(SORTED) GOTO 20
      GOTO 10
 20   CONTINUE


C ... Compute the area of the polygon (overlaping area)

      ASHARE  = 0.0

      F(1,1)  = C(1,1)
      F(2,1)  = C(2,1)
      F(3,1)  = C(3,1)

      DO I = 2,NPOINT-1

         F(1,2) = C(1,I)
         F(2,2) = C(2,I)
         F(3,2) = C(3,I)
         F(1,3) = C(1,I+1)
         F(2,3) = C(2,I+1)
         F(3,3) = C(3,I+1)

         AX = F(1,2) - F(1,1)
         AY = F(2,2) - F(2,1)
         AZ = F(3,2) - F(3,1)

         BX = F(1,3) - F(1,1)
         BY = F(2,3) - F(2,1)
         BZ = F(3,3) - F(3,1)
      
         ASHARE = ASHARE + 0.5*SQRT((AY*BZ-AZ*BY)**2
     &                            + (AZ*BX-AX*BZ)**2
     &                            + (AX*BY-AY*BX)**2)

      ENDDO

      
      IF(APROJ > 0.0) THEN
         ASHARE = ASHARE/APROJ
      ELSE
         ASHARE = 0.0
      ENDIF

      RETURN
      END SUBROUTINE CAREA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTP(OUTPUT_LOCATION,FIRSTOUTP)

C ... OUTPUT FOR THE THREE-DIMENSIONAL SOLVER

      USE MPI

      USE CHARACTERS

      USE CONSTANTS,   ONLY : PII,EPS

      USE INTEGERS,    ONLY : IPRO,NPRO,NBCS,MAXB,IB,MAXW,MAXEB,NMOV,
     &                        MAXSS,MGM,MBPRO,MAXSB,IBF,IREPEA,NREPEA

      USE MAIN_ARRAYS

      USE NS3CO

      USE TYPE_ARRAYS, ONLY : SIXDOF
     
      USE FLIGHT , ONLY : NGRIFL,OSKU,XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIR,THETAR,PHIR,PSIM,TRMODE,MVSHIP,ROTA1,ROTB1,CONEA,
     &     ACTDISK,FLYOBJ,IFAKEP,SHIP

      USE BLADE_VARIABLES , ONLY : TDBL,DTB

      IMPLICIT NONE

      INTEGER :: N,INSCA,ICY22,NCELLS,NBLCOM,NGL,NTOTIN,IG1,IF1,IR1,IC1,
     &    ISU,N1,IBTYPE,INC2,L,ISTRID,JSTRID,KSTRID,IQ1,IE1,IQQ1,NN,
     &    NS,IHP,I,KMAXXP,NSP,KX1,KX2,KY1,KY2,KZ1,KZ2,ISTR,JSTR,KSTR,
     &    ISP,IWALL,J,K,IT1,IL1,IK1,IM1,JM1,KM1,NMAX,IMAXP1,JMAXP1,
     &    KMAXP1,ILE,INODE,IPT,IMEM,NCG,IERR,KOKO1,NROOT,NR1,NNMAX,
     &    II,JJ,KK,NBLOCM,IERRCODE,ERRORCODE,ISD,IGN,IGT,IPHASE,IGR,
     &    ICHARP,III,NTOT1,RC,OUTPUT_LOCATION,IUNIT,ITIMES_DUMMY,NMOVG,
     &    NOFSTEPS, ISTART, INCREMENT

      INTEGER STATUS(MPI_STATUS_SIZE)

      REAL :: WCALI,WCOMI,RTIMEI,RTIMEC,TOMTOT,CXTOT,CYTOT,CZTOT,
     &    DXTOT,DYTOT,DZTOT,CMXTOT,CMYTOT,CMZTOT,PAVE,TAVE,ROAVE,EAVE,
     &    EINAVE,UUAVE,RKAVE,REAVE,RWAVE,SUP,TUP,SUR,SUI,SUU,SUK,
     &    DX,DY,DZ,QT,QW,QH,CLB,CDB,CSB,CLP,CDP,VX,VY,VZ,CLV,CDV,VXTOT,
     &    VYTOT,VZTOT,RF1,RF2,QMEIN2,QMEOU2,HEADT,EFFIC,TOMTOG,VISAVE,
     &    PROAVE,EPSAVE,EPWAVE,SUV,SUPRO,SUEPS,SUEPW,BKOKO1,
     &    SUE,CXIO,CYIO,CZIO,CLBIO,CDBIO,CSBIO,CMXIO,CMYIO,CMZIO,
     &    CDBTOT,CLBTOT,CSBTOT,CDPTOT,CDVTOT,CLPTOT,CLVTOT,
     &    RX,RY,RZ,CLR,CDR,RXTOT,RYTOT,RZTOT,CDRTOT,CLRTOT,
     &    VMX,VMY,VMZ,VMXTOT,VMYTOT,VMZTOT,NOSTEX,NOSTEY,NOSTEZ,
     &    T_DUMMY,SUA,ROGAVE,SUG,EVAPGE

      LOGICAL :: MASTER, SIMULATION_END, FIRSTOUTP

      CHARACTER(LEN=1) :: CLEVW

      MASTER = IPRO == 1

      SIMULATION_END = OUTPUT_LOCATION == 1
      
      INSCA = 1

      IF(ITURB >= 6 .AND. ITURB /= 8) INSCA = 7


      IF(SIMULATION_END) THEN  ! Write only at the end of simulation

      IF(PARALLEL) THEN
         ICY22 = (ICYTOT - ICYOLD)
         WCALI = WCAL/ICY22
         WCOMI = WCOM/ICY22
         WRITE(45,*)
         WRITE(45,*) 'Total times in proces per iteration sweep in ',
     +        'process ',IPRO
         WRITE(45,*) '  TTIME/cycle     TCOM/cycle   commucation/total',
     +        ' communication/calcul.'
         WRITE(45,1111) REAL(WCALI,4),REAL(WCOMI,4),REAL(WCOM/WCAL,4),
     +        REAL(WCOM/(WCAL-WCOM),4)
      ENDIF

 1111 FORMAT(1X,6E15.5)

      RTIMEI = (RTIME-RTIME1)/(ICYTOT - ICYOLD)

      NCELLS = 0
      NBLOCM = NBLOCK
      DO 1 N = 1,NBLOCM
1     NCELLS = NCELLS + IMAX(1,N)*JMAX(1,N)*KMAX(1,N)

      RTIMEC = RTIMEI/NCELLS*1000000
      WRITE(3,*)
      WRITE(3,*) ' ITERATION IS FINISHED'
      IF(TIMEL) WRITE(3,*)' TIME = ',T,' DT = ',DT,' STEPS = ',ITIMES
      IF(.NOT.TIMEL) WRITE(3,*) ' CYCLE = ',ICYCLE
      WRITE(3,*) ' DROMAX = ',DROMAX
c     WRITE(4,11) RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS

      IF(TIMEL .AND. IPRO == 1) WRITE(*,12) ITIMES,T
 12   FORMAT(4X,'STEPS                    =',I9,3X,'TIME   =',E12.5)
C ... PRINT THE CONVERGENCE HISTORY ON THE SCREEN
      CALL PRIRES(RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS,CD,CL,CMZ,
     +     IPRO,NPRO,WCALI,WCOMI)
      IF(TIMEL .AND. IPRO == 1) THEN
         WRITE(*,4040) REAL(TAVER,4),NAVER
         WRITE(4,4040) REAL(TAVER,4),NAVER
 4040    FORMAT(' Averaging is based on ',E11.4,
     +        ' seconds of simulation with ',I4,' time steps.'/
     +        ' Results are in files: RES.AVE and FUN.AVE.')
         ENDIF

C *** NON-STANDARD OUTPUT **********************************************
           
      WRITE(4,*)
      WRITE(4,*)'------ AT THE END BLOCK BY BLOCK -----------------'
      WRITE(4,*)'--------------------------------------------------'

      TOMTOT  = 0.
      CXTOT   = 0.
      CYTOT   = 0.
      CZTOT   = 0.
      DXTOT   = 0.
      DYTOT   = 0.
      DZTOT   = 0.
      RXTOT   = 0.
      RYTOT   = 0.
      RZTOT   = 0.
      CMXTOT  = 0.
      CMYTOT  = 0.
      CMZTOT  = 0.
      VMXTOT  = 0.
      VMYTOT  = 0.
      VMZTOT  = 0.
      PAVE    = 0.
      TAVE    = 0.
      ROAVE   = 0.
      ROGAVE  = 0.
      EVAPGE  = 0.
      EAVE    = 0.
      EINAVE  = 0.
      UUAVE   = 0.
      RKAVE   = 0.
      REAVE   = 0. ! dissipation tilde
      RWAVE   = 0. ! wall effects
      SUP     = 0.
      TUP     = 0.
      SUR     = 0.
      SUA     = 0.
      SUG     = 0.
      SUE     = 0.
      SUI     = 0.
      SUU     = 0.
      SUK     = 0.

      ISD     = 1

      RF1     = AREF*REFPRE
      RF2     = AREF*REFPRE*CHLREF

C ... Here the methods of the integrated parameters should be described

      WRITE(4,5000)
 5000 FORMAT(//'   DESCRIPTION OF THE CALCULATED VARIABLES'/
     +         '   ======================================='/
     +  ' These might be out of the date, so most recent modification'/
     +     ,' can be seen in source code ns3c.f (subroutine OUTP) '/
     + ' Averaged pressure is calculated T= Sum(p_i*V_i)/Sum(V_i)'/
     + ' Averaged temperature is calculated T= Sum(T_i*V_i)/Sum(V_i)'/
     +     ' Total mass            = Sum(ro_i*V_i)'/
     +     ' Total energy          = Sum(E_i*V_i)'/
     +     ' Total internal energy = Sum(e_i*V_i)'/
     +     ' Total kinetic energy  = Sum(.5*(u^2+v^2+w^2)_i*ro_i*V_i)'/
     +     ' Total turbulent energy= Sum(k_i*ro_i*V_i)'/)
C ... NON-STANDARD BLOCK LOOP BEGINS
C     OPEN(8, FILE= 'PLOTPI'//PRN,STATUS='UNKNOWN',FORM='FORMATTED')

C ... Non-standard block loop begins

      DO 1000 N = 1,NBLOCK

      NGL     = NPROCE(1+N,IPRO) !Global block number
      NTOTIN  = IMAX(1,N)*JMAX(1,N)*KMAX(1,N)
      IG1     = IG(1,N)
      IGN     = IG1 + NTOT(1,N)
      IF1     = JF(1,N)
      IR1     = IR(1,N)
      IC1     = IC(1,N)
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

C ... CALCULATE THE AVERAGE STATE IN EACH BLOCK

      CALL SETV12(F1E(IG1),E(IG1),NTOT(1,N))  ! STORE TOTAL ENERGY
      CALL PRIMVE(U(IG1),V(IG1),W(IG1),RK(IG1),RO(IG1),F1E(IG1),
     + XC(IG1),YC(IG1),ZC(IG1),NTOT(1,N),FRSDEN,ITURB)
      CALL MULV12(F1E(IG1),RO(IG1),NTOT(1,N)) ! CREATE INT.ENERGY/VOLUME

      CALL SUMFU(SUP,   P(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      CALL SUMFU(TUP,TEMP(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      CALL SUMFU(SUR,  RO(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      IF(MULPHL .OR. CAVITL) THEN
      F1R(IG1:IGN)   = VAR(IG1:IGN)%ALFA(2)
      F1RM(IG1:IGN)  = PRO(IG1:IGN)%RO(1)
      CALL SUMFA(SUA,F1R(IG1),F1RM(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      ! write(6,*) sua
      F1R(IG1:IGN)   = VAR(IG1:IGN)%ARO(2)
      CALL SUMFU(SUA,F1R(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      F1R(IG1:IGN)   = VAR(IG1:IGN)%EVAPR(2)
      CALL SUMFU(SUG,F1R(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)

      ENDIF
      CALL SUMFU(SUE,   E(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      CALL SUMFU(SUI, F1E(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      CALL SUUFU(SUU, RM(IG1),RN(IG1),RW(IG1),RO(IG1),VOL(IG1),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),0.)
        IF(ITURB >= 3 .AND. ITURB /= 8) THEN
       CALL SUMFU(SUK,  RK(IR1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),RKLIM)
      ENDIF

      PAVE    = SUP + PAVE  
      TAVE    = TUP + TAVE  
      ROAVE   = SUR + ROAVE 
      ROGAVE  = SUA + ROGAVE
      EVAPGE  = EVAPGE + SUG
      EAVE    = SUE + EAVE  
      EINAVE  = SUI + EINAVE
      UUAVE   = SUU + UUAVE
      RKAVE   = SUK + RKAVE

      WRITE(4,*)
      WRITE(4,*) 'BLOCK = ',NGL
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED BALANCES --------------------'
      WRITE(4,*)
      WRITE(4,*) 'AVERAGE PRESSURE    = ',REAL(SUP/VOLN(N),4), ' N/m2'
      WRITE(4,*) 'AVERAGE TEMPERATURE = ',REAL(TUP/VOLN(N),4),' K'
      WRITE(4,*) 'TOTAL MASS          = ',REAL(SUR,4),       ' kg'
      IF(MULPHL) THEN
      WRITE(4,*) 'TOTAL GAS MASS      = ',REAL(SUA,4),      ' kg'
      WRITE(4,*) 'TOTAL EVAPORATION   = ',REAL(SUG,4),      ' kg/s'
      ENDIF
      WRITE(4,*) 'TOTAL ENERGY        = ',REAL(SUE,4),   ' kg/s'
      WRITE(4,*) 'TOTAL INT. ENERGY   = ',REAL(SUI,4),      ' J'
      WRITE(4,*) 'TOTAL KINETIC ENERGY= ',REAL(SUU,4),      ' J'
      IF(ITURB >= 3 .AND. ITURB /= 8)
     +WRITE(4,*) 'TOTAL TURB. ENERGY  = ',REAL(SUK,4),      ' J'
      WRITE(4,*)
      WRITE(4,*)'------------- INDIVIDUAL PATHCES -----------------'
c9700  FORMAT(' PATCH',3X,'CX',9X,'CY',9X,'CZ',8X,'CMX',8X,'CMY',
c     + 8X,'CMZ',6X,'TOMEGA')
9700  FORMAT(' PATCH',4X,'CX',9X,'CY',9X,'CZ',8X,'CMX',8X,'CMY',
     + 8X,'CMZ',6X,'TOMEGA')
      WRITE(4,9700)
      WRITE(4,*)
c9701  FORMAT(I3,2X,7(E10.4,1X))
9701  FORMAT(I5,1X,7(E10.3,1X))

      ISU     = ISD + NPATCH(N)
      TOM     = 0.
      CX      = 0.
      CY      = 0.
      CZ      = 0.
      DX      = 0.
      DY      = 0.
      DZ      = 0.
      RX      = 0.
      RY      = 0.
      RZ      = 0.
      QT      = 0.
      QW      = 0.
      QH      = 0.
      CMX     = 0.
      CMY     = 0.
      CMZ     = 0.
      VMX     = 0.
      VMY     = 0.
      VMZ     = 0.      

      DO 8199 N1 = ISD,ISU-1
      CXB(N1) = CXB(N1)/(AREF*REFPRE)
      CYB(N1) = CYB(N1)/(AREF*REFPRE)
      CZB(N1) = CZB(N1)/(AREF*REFPRE)
      DXB(N1) = DXB(N1)/(AREF*REFPRE)
      DYB(N1) = DYB(N1)/(AREF*REFPRE)
      DZB(N1) = DZB(N1)/(AREF*REFPRE)
      BUX(N1) = BUX(N1)/(AREF*REFPRE)
      BUY(N1) = BUY(N1)/(AREF*REFPRE)
      BUZ(N1) = BUZ(N1)/(AREF*REFPRE)
      CMXB(N1)= CMXB(N1)/(AREF*CHLREF*REFPRE)
      CMYB(N1)= CMYB(N1)/(AREF*CHLREF*REFPRE)
      CMZB(N1)= CMZB(N1)/(AREF*CHLREF*REFPRE)
      VMXB(N1)= VMXB(N1)/(AREF*CHLREF*REFPRE)
      VMYB(N1)= VMYB(N1)/(AREF*CHLREF*REFPRE)
      VMZB(N1)= VMZB(N1)/(AREF*CHLREF*REFPRE)

      INC2    = 1+(N1-ISD)*IC9+IC1-1
      IBTYPE  = ICON(INC2)

      IF(IBTYPE >= 8 .AND. IBTYPE <= 10) THEN ! BCP
      INC2 = 24+(N1-ISD)*IC9+IC1
      WRITE(4,9701) ICON(INC2),CXB(N1),CYB(N1),CZB(N1),
     + CMXB(N1),CMYB(N1),CMZB(N1),TOMEGA(N1)
      ENDIF
     
      TOM     = TOM + TOMEGA(N1)
      CX      = CX  + CXB(N1)
      CY      = CY  + CYB(N1)
      CZ      = CZ  + CZB(N1)
      DX      = DX  + DXB(N1)
      DY      = DY  + DYB(N1)
      DZ      = DZ  + DZB(N1)
      RX      = RX  + BUX(N1)      
      RY      = RY  + BUY(N1)  
      RZ      = RZ  + BUZ(N1)  
      QT      = QT  + QTB(N1)
      QW      = QW  + QWB(N1)
      QH      = QH  + QHB(N1)
      CMX     = CMX + CMXB(N1)
      CMY     = CMY + CMYB(N1)
      CMZ     = CMZ + CMZB(N1)
      VMX     = VMX + VMXB(N1)
      VMY     = VMY + VMYB(N1)
      VMZ     = VMZ + VMZB(N1)

8199  CONTINUE

      WRITE(4,*)
      WRITE(4,*)'-----PRESSURE AND VISCOUS FORCES BASED FOR '//
     +     'INDIVIDUAL PATHCES ------'
      WRITE(4,9702)
c9702  FORMAT(' PATCH',3X,'DX',9X,'DY',9X,'DZ',9X,'VX',9X,'VY',9X,'VZ')
9702  FORMAT(' PATCH',4X,'DX',9X,'DY',9X,'DZ',9X,'VX',9X,'VY',9X,'VZ')
      WRITE(4,*)

      DO 8200 N1 = ISD,ISU-1
      INC2   = 1+(N1-ISD)*IC9+IC1-1
      IBTYPE = ICON(INC2)
      IF(IBTYPE >= 8 .AND. IBTYPE <= 10) THEN ! BCP
      INC2 = 24+(N1-ISD)*IC9+IC1
      WRITE(4,9701) ICON(INC2),DXB(N1),DYB(N1),DZB(N1),
     +       CXB(N1)- DXB(N1),CYB(N1)- DYB(N1),CZB(N1)- DZB(N1)
      ENDIF
8200  CONTINUE

      WRITE(4,*)
      WRITE(4,*)'------- BUOYANCY FORCES AND VISCOUS MOMENTS -------'
      WRITE(4,9703)
c9703  FORMAT(' PATCH',2X,'BUX',8X,'BUY',8X,'BYZ',8X,'MVX',8X,'MVY',
c     + 8X,'MVZ',6X,'TOMEGA')
9703  FORMAT(' PATCH',3X,'BUX',8X,'BUY',8X,'BYZ',8X,'MVX',8X,'MVY',
     + 8X,'MVZ',6X,'TOMEGA')
      WRITE(4,*)

      DO 8201 N1 = ISD,ISU-1
      INC2       = 1+(N1-ISD)*IC9+IC1-1
      IBTYPE     = ICON(INC2)
      IF(IBTYPE >= 8 .AND. IBTYPE <= 10) THEN ! BCP
      INC2 = 24+(N1-ISD)*IC9+IC1
      WRITE(4,9701) ICON(INC2),BUX(N1),BUY(N1),BUZ(N1),
     +       VMXB(N1),VMYB(N1),VMZB(N1),TOMEGA(N1)
      ENDIF
8201  CONTINUE

      WRITE(4,*)
      WRITE(4,*)'----- HEAT FLUXES OR HEAT TRANSFER COEFFICIENTS ---'
      WRITE(4,9704)
c9704  FORMAT(' PATCH',3X,'Qtot',6X,'Qwall',7X,'Qheat')
9704  FORMAT(' PATCH',4X,'Qtot',6X,'Qwall',7X,'Qheat')
      WRITE(4,*)

      DO 8202 N1 = ISD,ISU-1
      INC2       = 1+(N1-ISD)*IC9+IC1-1
      IBTYPE     = ICON(INC2)
      IF(IBTYPE >= 8 .AND. IBTYPE <= 10 .OR. IBTYPE == 15) THEN ! BCP
      INC2 = 24+(N1-ISD)*IC9+IC1
      WRITE(4,9701) ICON(INC2),QTB(N1),QWB(N1),QHB(N1)
      ENDIF
8202  CONTINUE

      WRITE(4,*)

      ISD     = ISU  ! The patch number in the next block

      CXTOT   = CXTOT  + CX
      CYTOT   = CYTOT  + CY
      CZTOT   = CZTOT  + CZ
      DXTOT   = DXTOT  + DX
      DYTOT   = DYTOT  + DY
      DZTOT   = DZTOT  + DZ
      RXTOT   = RXTOT  + RX
      RYTOT   = RYTOT  + RY
      RZTOT   = RZTOT  + RZ
      CMXTOT  = CMXTOT + CMX
      CMYTOT  = CMYTOT + CMY
      CMZTOT  = CMZTOT + CMZ
      VMXTOT  = VMXTOT + VMX
      VMYTOT  = VMYTOT + VMY
      VMZTOT  = VMZTOT + VMZ 
      TOMTOT  = TOMTOT + TOM

      CLB     = (-CX*SIN(ALPHA) + CY*COS(ALPHA))
      CDB     = CX*COS(ALPHA)*COS(BETA) + CY*SIN(ALPHA)*COS(BETA)
     +         +CZ*SIN(BETA)
      CSB     = CX*COS(ALPHA)*SIN(BETA) + CY*SIN(ALPHA)*SIN(BETA)
     +         +CZ*COS(BETA)
      CLP     = (-DX*SIN(ALPHA) + DY*COS(ALPHA))
      CDP     = DX*COS(ALPHA)*COS(BETA) + DY*SIN(ALPHA)*COS(BETA)
     +         +DZ*SIN(BETA)
      CLR     = (-RX*SIN(ALPHA) + RY*COS(ALPHA))
      CDR     = RX*COS(ALPHA)*COS(BETA) + RY*SIN(ALPHA)*COS(BETA)
     +         +RZ*SIN(BETA)
      VX      = CX  - DX
      VY      = CY  - DY
      VZ      = CZ  - DZ
      CLV     = CLB - CLP
      CDV     = CDB - CDP

      WRITE(4,*)
     +'CX = ',REAL(CX,4),' CY = ',REAL(CY,4),' CZ = ',REAL(CZ,4)
      WRITE(4,*)
     +'DX = ',REAL(DX,4),' DY = ',REAL(DY,4),' DZ = ',REAL(DZ,4)
      WRITE(4,*)
     +'VX = ',REAL(VX,4),' VY = ',REAL(VY,4),' VZ = ',REAL(VZ,4)
      WRITE(4,*)
     +'RX = ',REAL(RX,4),' RY = ',REAL(RY,4),' RZ = ',REAL(RZ,4)
      WRITE(4,*)
     +'MX = ',REAL(CMX,4),' MY = ',REAL(CMY,4),' MZ = ',REAL(CMZ,4)
      WRITE(4,*)
     +'MVX= ',REAL(VMX,4),' MVY= ',REAL(VMY,4),' MVZ= ',REAL(VMZ,4)
      WRITE(4,*) 'CD  = ',REAL(CDB,4),'  CL  = ',REAL(CLB,4)
      WRITE(4,*) 'CDP = ',REAL(CDP,4),'  CLP = ',REAL(CLP,4)
      WRITE(4,*) 'CDV = ',REAL(CDV,4),'  CLV = ',REAL(CLV,4)
      WRITE(4,*) 'CLR = ',REAL(CLR,4),'  CDR = ',REAL(CDR,4)
      WRITE(4,*) 'TOMEGA =',REAL(TOM,4),' W'

1000  CONTINUE   !  END OF NON-STANDARD BLOCK LOOP
         
C      IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C      IF(IPRO == 1) THEN
C         CALL COLPLO(ZZZ,MAXW,SPLIT,NOB,IDOB,JDOB,KDOB,NSSB,ISSB,NB,
C     +        NPROCE,MBPRO,NPRO)
C      ENDIF
C ... voimat loppu (ESa teki loppuun 31.12.2004: FVSRF)

C ... WRITE THE TOTAL COEFFICIENTS ON FORCES FILE

      PAVE     = PAVE/VOLTOT
      TAVE     = TAVE/VOLTOT
      VXTOT    = CXTOT - DXTOT
      VYTOT    = CYTOT - DYTOT
      VZTOT    = CZTOT - DZTOT
      
      CLBTOT  = (-CXTOT*SIN(ALPHA) + CYTOT*COS(ALPHA))
      CDBTOT  = CXTOT*COS(ALPHA)*COS(BETA) + CYTOT*SIN(ALPHA)*COS(BETA)
     +        + CZTOT*SIN(BETA)
      CSBTOT  = CXTOT*COS(ALPHA)*SIN(BETA) + CYTOT*SIN(ALPHA)*SIN(BETA)
     +        + CZTOT*COS(BETA)
      CLPTOT  = (-DX*SIN(ALPHA) + DY*COS(ALPHA))
      CDPTOT  = DXTOT*COS(ALPHA)*COS(BETA) + DYTOT*SIN(ALPHA)*COS(BETA)
     +         +DZTOT*SIN(BETA)
      CLVTOT  = CLBTOT - CLPTOT
      CDVTOT  = CDBTOT - CDPTOT
      CLRTOT  = (-RX*SIN(ALPHA) + RY*COS(ALPHA))
      CDRTOT  = RXTOT*COS(ALPHA)*COS(BETA) + RYTOT*SIN(ALPHA)*COS(BETA)
     +         +RZTOT*SIN(BETA)

      CXIO    = (QMXIN+QMXOUT)/RF1
      CYIO    = (QMYIN+QMYOUT)/RF1
      CZIO    = (QMZIN+QMZOUT)/RF1

      CMXIO   = (CMXIN+CMXOUT)/RF2
      CMYIO   = (CMYIN+CMYOUT)/RF2
      CMZIO   = (CMZIN+CMZOUT)/RF2

      CLBIO   = (-CXIO*SIN(ALPHA) + CYIO*COS(ALPHA))
      CDBIO   = CXIO*COS(ALPHA)*COS(BETA) + CYIO*SIN(ALPHA)*COS(BETA)
     +         +CZIO*SIN(BETA)
      CSBIO   = CXIO*COS(ALPHA)*SIN(BETA) + CYIO*SIN(ALPHA)*SIN(BETA)
     +         +CZIO*COS(BETA)

      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED COEFFICIENTS ----------------
     +------------'
      WRITE(4,*)
      WRITE(4,*)  'CD = ',REAL(CDBTOT,4),' CL  = ',REAL(CLBTOT,4),
     +           ' CS = ',REAL(CSBTOT,4)
      WRITE(4,*) 'CDp = ',REAL(CDPTOT,4),' CDv = ',REAL(CDVTOT,4),
     +          ' CDr = ',REAL(CDRTOT,4)
      WRITE(4,*) 'CLp = ',REAL(CLPTOT,4),' CLv = ',REAL(CLVTOT,4),
     +          ' CLr = ',REAL(CLRTOT,4)
      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED COEFFICIENTS FOR SOLIDS -----
     +------------'
      WRITE(4,*)
      WRITE(4,*) 'CD  = ',REAL(CDBTOT-CDBIO,4),
     +          ' CL  = ',REAL(CLBTOT-CLBIO,4),
     +          ' CS = ' ,REAL(CSBTOT-CDBIO,4)
      WRITE(4,*) 'CDp = ',REAL(CDPTOT-CDBIO,4),
     +          ' CDv = ',REAL(CDVTOT,4),' CDr = ',REAL(CDRTOT,4)
      WRITE(4,*) 'CLp = ',REAL(CLPTOT-CLBIO,4),' CLv = ',REAL(CLVTOT,4),
     +          ' CLr = ',REAL(CLRTOT,4)
      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED COEFFICIENTS FOR INLETS AND O
     +UTLETS -------'
      WRITE(4,*)
      WRITE(4,*) 'CD = ',REAL(CDBIO,4),' CL = ',REAL(CLBIO,4),
     +          ' CS = ',REAL(CSBIO,4)
      WRITE(4,*) 'MX = ',REAL(CMXIO,4),' MY = ',REAL(CMYIO,4),
     +          ' MZ = ',REAL(CMZIO,4)
      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED FORCES ----------------------
     +----------'
      WRITE(4,*)
      WRITE(4,*) 'CX = ',REAL(CXTOT,4),' CY = ', REAL(CYTOT,4),
     +          ' CZ = ',REAL(CZTOT,4)
      WRITE(4,*) 'DX = ',REAL(DXTOT,4),' DY = ', REAL(DYTOT,4),
     +          ' DZ = ',REAL(DZTOT,4)
      WRITE(4,*) 'VX = ',REAL(VXTOT,4),' VY = ', REAL(VYTOT,4),
     +          ' VZ = ',REAL(VZTOT,4)
      WRITE(4,*) 'RX = ',REAL(RXTOT,4),' RY = ', REAL(RYTOT,4),
     +          ' RZ = ',REAL(RZTOT,4)
      WRITE(4,*) 'MX = ',REAL(CMXTOT,4),' MY = ',REAL(CMYTOT,4),
     +          ' MZ = ',REAL(CMZTOT,4)
      WRITE(4,*) 'MVX= ',REAL(VMXTOT,4),' MVY= ',REAL(VMYTOT,4),
     +          ' MVZ= ',REAL(VMZTOT,4)
      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED FORCES IN SI-UNITS-----------
     +-----------------'

      WRITE(4,*)
      WRITE(4,*) 'CX = ',REAL(CXTOT*RF1,4),' N',' CY = ',
     +REAL(CYTOT*RF1,4),' N',' CZ = ',REAL(CZTOT*RF1,4),' N'
      WRITE(4,*) 'DX = ',REAL(DXTOT*RF1,4),' N',' DY = ',
     +REAL(DYTOT*RF1,4),' N',' DZ = ',REAL(DZTOT*RF1,4),' N'
      WRITE(4,*) 'VX = ',REAL(VXTOT*RF1,4),' N',
     +' VY = ',REAL(VYTOT*RF1,4),' N',' VZ = ',REAL(VZTOT*RF1,4),' N'
      WRITE(4,*) 'RX = ',REAL(RXTOT*RF1,4),' N',
     +' RY = ',REAL(RYTOT*RF1,4),' N',' RZ = ',REAL(RZTOT*RF1,4),' N'
      WRITE(4,*) 'MX = ',REAL(CMXTOT*RF2,4),'Nm',
     +' MY = ',REAL(CMYTOT*RF2,4),'Nm',' MZ = ',REAL(CMZTOT*RF2,4),'Nm'
      WRITE(4,*) 'MVX= ',REAL(VMXTOT*RF2,4),'Nm',
     +' MVY= ',REAL(VMYTOT*RF2,4),'Nm',' MVZ= ',REAL(VMZTOT*RF2,4),'Nm'
      WRITE(4,*)
      WRITE(4,*) 'ABS.TOMEGA =',REAL(ABS(TOMTOT),4),' W'
      IF (GRAVIL) THEN
         NOSTEX=RXTOT*RF1
         NOSTEY=RYTOT*RF1
         NOSTEZ=RZTOT*RF1
         IF (ABS(GX) < 1.0E-5) THEN
            NOSTEX = 0.
         ENDIF 
         IF (ABS(GY) < 1.0E-5) THEN
            NOSTEY = 0.
         ENDIF 
         IF (ABS(GZ) < 1.0E-5) THEN
            NOSTEZ = 0.
         ENDIF
         WRITE(4,*) 'VOLUME INSIDE SOLID PATCHES =',
     +        SQRT(NOSTEX**2+NOSTEY**2+NOSTEZ**2) /
     +        (SQRT(GX**2+GY**2+GZ**2)*FRSDEN),'m^3'
         WRITE(4,*) '(VOLUME IS BASED ON HYDROSTATIC FORCE)'
       ENDIF                    ! IF (GRAVIL)
      WRITE(4,*)
      WRITE(4,*)
      WRITE(4,*) 'THE FOLLOWING VALUES ARE ACCURATE. MORE ABOUT ',
     +'FLUX-BASED BALANCES'
      WRITE(4,*) 'MAY BE FOUND IN MCHECK-FILE'
      WRITE(4,*)'------------------------------------------------------
     +----------------------'
      WRITE(4,*) 'MASS FLUX IN   =',REAL(QMFIN,4),' kg/s',
     + '  MASS FLUX OUT   =',REAL(QMFOUT,4),'kg/s'
      IF(MULPHL) THEN
       WRITE(4,*) 'GAS M.FLUX IN  =',REAL(QGFIN,4),' kg/s',
     + '  GAS M.FLUX OUT  =',REAL(QGFOUT,4),'kg/s'
       WRITE(4,*) 'GAS ENERGY FLUX IN =',REAL(QGEIN,4),'W',
     + '  GAS ENERGY FLUX OUT =',REAL(QGEOUT,4),'W'
       WRITE(4,*) 'LIQ M.FLUX IN  =',REAL(QLFIN,4),' kg/s',
     + '  LIQ M.FLUX OUT  =',REAL(QLFOUT,4),'kg/s'
       WRITE(4,*) 'LIQ ENERGY FLUX IN =',REAL(QLEIN,4),'W',
     + '  LIQ ENERGY FLUX OUT =',REAL(QLEOUT,4),'W'
      ENDIF
      WRITE(4,*) 'MECH.ENERGY IN =',REAL(QMEIN,4),' W', 
     + '     MECH.ENERGY OUT =',REAL(QMEOUT,4),'W'
      WRITE(4,*) 'X-MOMENTUM IN = ',REAL(QMXIN,4),' N', 
     + '     X-MOMENTUM OUT =',REAL(QMXOUT,4),' N'
      WRITE(4,*) 'Y-MOMENTUM IN = ',REAL(QMYIN,4),' N', 
     + '     Y-MOMENTUM OUT =',REAL(QMYOUT,4),' N'
      WRITE(4,*) 'Z-MOMENTUM IN = ',REAL(QMZIN,4),' N', 
     + '     Z-MOMENTUM OUT =',REAL(QMZOUT,4),' N'
      QMEIN2    = QMEIN /(QMFIN +1.E-20)
      QMEOU2    = QMEOUT/(QMFOUT+1.E-20)
      WRITE(4,*) 'AVERAGE SPECIFIC MECH.ENERGY IN  =',
     +            REAL(QMEIN2,4),' J/kg'
      WRITE(4,*) 'AVERAGE SPECIFIC MECH.ENERGY OUT =',
     +            REAL(QMEOU2,4),' J/kg'
      IF(TOMTOT /= 0.) THEN
      EFFIC = 100.*ABS(QMFIN*(QMEOU2-QMEIN2)/TOMTOT)
      HEADT = (QMEOU2 - QMEIN2)/G0
      WRITE(4,*) 'ESTIMATED HEAD                   =',
     +            REAL(HEADT,4),' m'
      WRITE(4,*) 'ESTIMATED EFFICIENCY             =',
     +            REAL(EFFIC,4),' %'
      ENDIF
      WRITE(4,*)

C ... CALCULATE AUXILIARLY FORCES.

      WRITE(4,*)'------------- INTEGRATED SUBFORCES -------------------
     +-----------'
      CALL AUXFOR(155,CXB,CYB,CZB,CMXB,CMYB,CMZB,DXB,DYB,DZB,QTB,QWB,
     +     QHB,TOMEGA,ALPHA,BETA,ICON,BOUNDF,NBCS,IPRO,1,0.,0,
     +     APATCH,XMOM,YMOM,ZMOM,REFPRE,AREF,CHLREF)
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED BALANCES --------------------'
      WRITE(4,*)
      WRITE(4,*) 'AVERAGE PRESSURE    = ',REAL(PAVE,4),  ' N/m2'
      WRITE(4,*) 'AVERAGE TEMPERATURE = ',REAL(TAVE,4),  ' K'
      WRITE(4,*) 'TOTAL MASS          = ',REAL(ROAVE,4), ' kg'
      WRITE(4,*) 'TOTAL GAS MASS      = ',REAL(ROGAVE,4),' kg'
      WRITE(4,*) 'TOTAL EVAPORATION   = ',REAL(EVAPGE,4),' kg/s'
      WRITE(4,*) 'TOTAL ENERGY        = ',REAL(EAVE,4),  ' J'
      WRITE(4,*) 'TOTAL INT. ENERGY   = ',REAL(EINAVE,4),' J'
      WRITE(4,*) 'TOTAL KINETIC ENERGY= ',REAL(UUAVE,4), ' J'
      IF(ITURB >= 3 .AND. ITURB /= 8)
     +WRITE(4,*) 'TOTAL TURB. ENERGY  = ',REAL(RKAVE,4), ' J'

C ... Update "TRAJECTORY_IGR.DAT" file in non-time accurate calculation

      IF (NGRIFL > 0 .AND. .NOT. TIMEL) THEN
      DO IGR = 1,NGRIFL
      IF(.NOT.MVSHIP) THEN   
      CALL PARTICLE_FORCES(OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     +     OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,ICON,
     +     XCG(IGR),YCG(IGR),ZCG(IGR),1,10,10,0,IGR,IGRID(:,3)) ! MOV surfaces
      ENDIF
      IF(SHIP .AND. IFAKEP(IGR) == 1)THEN
         CALL FAKEPROPELLER(IGR) !Fake propeller: X-force, Z-force and Y-moment
      ENDIF
      IF(INOUTL .AND. FLYOBJ) THEN ! thrust for the flying objects
      CALL PARTICLE_FORCES(OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +     OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,ICON,
     +     XCG(IGR),YCG(IGR),ZCG(IGR),1,3,5,0,IGR,IGRID(:,3))
      ENDIF                                           ! INL and OUT surfaces
C      IF(ACTDISK) THEN !    ! Act surfaces ! tassa ei tarvita particle_forces
c      CALL PARTICLE_FORCES(CXB,CYB,CZB,CMXB,CMYB,CMZB,
c     +     ALPHA,BETA,ICON,BOUNDF,NBCS,PARALLEL,IPRO,
c     +     APATCH,OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
c     +     OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
c     +     XCG(IGR),YCG(IGR),ZCG(IGR),XMOM,YMOM,ZMOM,
c     +     AREF,CHLREF,REFPRE,3,1,1,8,IGR,IGRID(:,3))
C      ENDIF                                           ! Act surfaces
      IF (IPRO == 1) THEN
c ... tassa kohtaa lasketaan helikopterin roottorin keskihaku ja gravitaatio
         CALL ROTOR_FORCES(OSKU(IGR)%RMASS,XCG(IGR),YCG(IGR),ZCG(IGR),
     +        OSKU(IGR)%ROTX,OSKU(IGR)%ROTY,OSKU(IGR)%ROTZ,
     +        OSKU(IGR)%VX,OSKU(IGR)%VY,OSKU(IGR)%VZ,
     +        OSKU(IGR)%CX,OSKU(IGR)%CY,OSKU(IGR)%CZ,
     +        OSKU(IGR)%CMX,OSKU(IGR)%CMY,OSKU(IGR)%CMZ,
     +        OSKU(IGR)%FXE,OSKU(IGR)%FYE,OSKU(IGR)%FZE,
     +        OSKU(IGR)%MXE,OSKU(IGR)%MYE,OSKU(IGR)%MZE,
     +        OSKU(IGR)%FXT,OSKU(IGR)%FYT,OSKU(IGR)%FZT,
     +        OSKU(IGR)%MXT,OSKU(IGR)%MYT,OSKU(IGR)%MZT,
     +        OSKU(IGR)%FXH,OSKU(IGR)%FYH,OSKU(IGR)%FZH,
     +        OSKU(IGR)%MXH,OSKU(IGR)%MYH,OSKU(IGR)%MZH,
     +        PSIR(IGR),THETAR(IGR),PHIR(IGR),OSKU(IGR)%RCG,
     +        OSKU(IGR)%IX,OSKU(IGR)%IY,OSKU(IGR)%IZ,
     +        OSKU(IGR)%IXY,OSKU(IGR)%IXZ,OSKU(IGR)%IYZ,
     +        ROTORW,GY,NBLOCK,TRMODE(IGR),IGR)
         CALL TRAJECTORY_WRITE(TDBL,DTB,XCGI(IGR),YCGI(IGR),ZCGI(IGR),
     +        XCG(IGR),YCG(IGR),ZCG(IGR),PSIR(IGR),THETAR(IGR),
     +        PHIR(IGR),PSIM(IGR),OSKU,IGR,TIMEL,TRMODE(IGR),
     +        ROTA1(IGR),ROTB1(IGR),CONEA(IGR))
         CALL TRAJECTORY_CONVERTER(IGR,TRMODE(IGR),TDBL,DTB)
         CALL TRAJECTORY_BACKUP(IGR)
      ENDIF
      ENDDO
      ENDIF

      IF(PARALLEL) CALL TOTBAL(PAVE,TAVE,ROAVE,EAVE,EINAVE,UUAVE,RKAVE,
     +     VOLTOT,ITURB,IPRO,QMFIN,QMFOUT,QMEIN,QMEOUT,TOMTOT,
     +     TOMTOG,QGFIN,QGFOUT,QGEIN,QGEOUT,MULPHL,ROGAVE,EVAPGE)

      DO 1997 L = 1,MAXB
         DM(L) = 0.
         DN(L) = 0.
         DW(L) = 0.
 1997 CONTINUE
  
C **********************************************************************
C ... CALCULATE THE COMPONENTS OF THE VORTICITY VECTOR
C ... BLOCK LOOP 1 BEGINS

C ... RE-EVALUATION OF PRESSURE (ESPECIALLY FOR EMPTY CORNERS)

      IF(ICYCLE == IPRINT) THEN
          IF(MULPHL) THEN
          IPHASE = 1
          F1RM(1:MAXB) = VAR(1:MAXB)%V(IPHASE)
          F1RN(1:MAXB) = VAR(1:MAXB)%W(IPHASE)
          ENDIF
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO) !Global block number
         IG1     = IG(1,N)
         CALL STATIM(1,N,2,.TRUE.)
         CALL PRINYS(3,EC, E(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

         IF(IL(1,N) >= 1 .AND. IL(1,N) <= 3) THEN
         WRITE(3,*)
         WRITE(3,*) ' Before MIR and CYC'
         IF(MULPHL) THEN
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),F1RN(1),IT(1,N),IL(1,N),
     +        0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ELSE
         CALL PRINYS(3,VC, V(IG1),IT(1,N),IL(1,N),0,IK(1,N),
     +   IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         ENDIF
         ENDIF
      ENDDO
      ENDIF  ! ICYCLE == IPRINT


      ENDIF  ! Write only at the end of simulation
      

C ... These where added to get updated velocity components

      CALL MIR(STRAIN,U,V,W,TEMP,C,VIS,CH,PDIFF,DE,DDEPS,TTS,FI,PRO,
     +     VAR,TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +     1,1,NBLOCK,1.,0,PRC)

      CALL CYC(U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     + 1,NBLOCK,F1R,F1RM,VAR)

c     IF(TWO_FLUIDL) THEN
c      VAR(1:MAXB)%U(IPHASE) = U(1:MAXB)
c      VAR(1:MAXB)%V(IPHASE) = V(1:MAXB)
c      VAR(1:MAXB)%W(IPHASE) = W(1:MAXB)
c      ENDIF
c      ENDDO


C ... Do these have some effect here?

      CALL CONEPA_SP(QWALL,'SOL')
      CALL CONEPA_SP(CPWALL,'SOL')

C ... Non standard way for computing space requirement.
      NN = IG(1,NBLOCK) + NTOT(1,NBLOCK) - 1

      ALLOCATE(RS1(NN))
      ALLOCATE(RS2(NN))
      ALLOCATE(RS3(NN))
      ALLOCATE(RS4(NN))
      ALLOCATE(RS5(NN))
      ALLOCATE(RS6(NN))

      ALLOCATE(CVV1(NN))
      ALLOCATE(CVV2(NN))
      ALLOCATE(CVV3(NN))
      ALLOCATE(CVV4(NN))
      ALLOCATE(CVV5(NN))
      ALLOCATE(CVV6(NN))
      ALLOCATE(CVV7(NN))
      ALLOCATE(CVV8(NN))
      ALLOCATE(CVV9(NN))

      ALLOCATE(GVX(NN))
      ALLOCATE(GVY(NN))
      ALLOCATE(GVZ(NN))

      ALLOCATE(VISTND(NN))

      ALLOCATE(UREL(NN))
      ALLOCATE(VREL(NN))
      ALLOCATE(WREL(NN))

      ALLOCATE(RMREL(NN))
      ALLOCATE(RNREL(NN))
      ALLOCATE(RWREL(NN))
      ALLOCATE(EREL(NN))

      IF(MULPHL) THEN
         ALLOCATE(MPUREL(NPHASE,NN))
         ALLOCATE(MPVREL(NPHASE,NN))
         ALLOCATE(MPWREL(NPHASE,NN))
      ELSE
         ALLOCATE(MPUREL(NPHASE,1))
         ALLOCATE(MPVREL(NPHASE,1))
         ALLOCATE(MPWREL(NPHASE,1))
      ENDIF
      
      DO 1998 N = 1,NBLOCK

      NGL     = NPROCE(1+N,IPRO) !Global block number
      IG1     = IG(1,N)
      IR1     = IR(1,N)
      IE1     = 1

      IF(ISTRES >= 1) IE1 = IG1
      IQQ1    = 1
      IF(STRESL) IQQ1 = IG1

C ... Vorticity and its components, strain (DE) and S11 with the ASM's
      
      CALL VORFUN(OHMI(IG1),DM(IG1),DN(IG1),DW(IG1),DE(IG1),
     1 S11(IQQ1,1),S11(IQQ1,2),S11(IQQ1,3),S11(IQQ1,4),
     2 S11(IQQ1,5),S11(IQQ1,6),PTUR(IG1),U(IG1),V(IG1),
     3 W(IG1),RK(IR1),REPS(IR1),DDEPS(IR1),EPS2(IG1),VIST(IG1),VIS(IG1),
     4 A1(IG1),A2(IG1),A3(IG1),VOL(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     5 A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     6 D1(IG1),D2(IG1),D3(IG1),XC(IG1),YC(IG1),ZC(IG1),IMAX(1,N),
     7 JMAX(1,N),KMAX(1,N),IN,JN,KN,ITURB,1,STRESL,ZZZ,VTRAN(IR1),MAXW,
     8 IDER(N),BIJ(IE1,1),MAXEB,ISTRES)

C ... HAT1 is production of turbulence for algebraic eddy-viscosity models
C     HAT2 is viscous dissipation, RS1...RS6 are the Reynolds stresses

      CALL TURFUN(RS1(IG1),RS2(IG1),RS3(IG1),RS4(IG1),RS5(IG1),RS6(IG1),
     1 U(IG1),V(IG1),W(IG1),RK(IR1),REPS(IR1),DDEPS(IR1),
     1 EPS2(IG1),VIST(IG1),VIS(IG1),HAT1(IG1),HAT2(IG1),
     2 A1(IG1),A2(IG1),A3(IG1),VOL(IG1),A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     3 A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     4 D1(IG1),D2(IG1),D3(IG1),XC(IG1),YC(IG1),ZC(IG1),
     5 IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,ITURB,1,ZZZ,VTRAN(IR1),
     6 MAXW,IDER(N))

      CALL DERFUN(OHMI(IG1),VORT(IG1,1),VORT(IG1,2),VORT(IG1,3),
     1 SHEAR(IG1,1),
     1 SHEAR(IG1,2),SHEAR(IG1,3),SHEAR(IG1,4),SHEAR(IG1,5),SHEAR(IG1,6),
     1 U(IG1),V(IG1),W(IG1),
     2 A1(IG1),A2(IG1),A3(IG1),VOL(IG1),D1(IG1),D2(IG1),D3(IG1),
     2 A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     3 A2XA(IG1),A2YA(IG1),A2ZA(IG1),A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     4 XC(IG1),YC(IG1),ZC(IG1),
     5 IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,IDER(N),1,STRAIN(IG1),
     6 VELLAP(IG1))

C ... Calculate Lighthill's tensor gradient

      CALL TIJFUN(TIJ(IG1),RO(IG1),U(IG1),V(IG1),W(IG1),
     + A1(IG1),A2(IG1),A3(IG1),VOL(IG1),
     + D1(IG1),D2(IG1),D3(IG1),
     + A1XA(IG1),A1YA(IG1),A1ZA(IG1),
     + A2XA(IG1),A2YA(IG1),A2ZA(IG1),
     + A3XA(IG1),A3YA(IG1),A3ZA(IG1),
     + XC(IG1),YC(IG1),ZC(IG1),
     + IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,IDER(N),1,NGL)

C ... Nondimensionalize turbulent viscosity for visualization and 
C ... change the sign in the vorticity vector (mersu)

      DO NN = 1,NTOT(1,N)
         N1 = NN + IG(1,N) - 1
         VISTND(N1) = VIST(N1)/(VIS(N1)+1.E-10)
      ENDDO

C ... and set the strains into reynolds stresses

      IF(ITURB >= 10 .AND. ITURB <= 19 .OR. ISTRES == 1) THEN

      DO NN = 1,NTOT(1,N)
         N1      = NN + IG(1,N) - 1
         RS1(N1) = S11(N1,1)    
         RS2(N1) = S11(N1,2) 
         RS3(N1) = S11(N1,3)
         RS4(N1) = S11(N1,4)
         RS5(N1) = S11(N1,5)
         RS6(N1) = S11(N1,6)
      ENDDO

      ELSEIF(ITURB >= 20) THEN

      DO NN = 1,NTOT(1,N)
         N1      = NN + IG(1,N) - 1
         RS1(N1) = FI(N1,1)
         RS2(N1) = FI(N1,2)
         RS3(N1) = FI(N1,3)
         RS4(N1) = FI(N1,4)
         RS5(N1) = FI(N1,5)
         RS6(N1) = FI(N1,6)
      ENDDO

      ENDIF ! ITURB >= 10 .AND. ITURB <= 19
 
 1998 CONTINUE

C ... END OF BLOCK LOOP 1
            
C ... BLOCK LOOP (INSIDE) FOR BOUNDARY VORTICITIES AND TURB.VISCOSITIES

      CALL CONEC(OHMI,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(TIJ,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(STRAIN,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL MIR(OHMI,DM,DN,DW,TEMP,C,VIS,CH,PDIFF,DE,DDEPS,TTS,FI,PRO,
     +     VAR,TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +     1,1,NBLOCK,-1.,0,PRC)
      CALL MIR(STRAIN,DM,DN,DW,TEMP,C,VIS,CH,PDIFF,DE,DDEPS,TTS,FI,PRO,
     +     VAR,TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +     1,1,NBLOCK,-1.,0,PRC)
      CALL MIR(TIJ,DM,DN,DW,TEMP,C,VIS,CH,PDIFF,DE,DDEPS,TTS,FI,PRO,
     +     VAR,TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +     1,1,NBLOCK,-1.,0,PRC)
      CALL MIRR(SHEAR(1,1),SHEAR(1,2),SHEAR(1,3),SHEAR(1,4),
     +     SHEAR(1,5),SHEAR(1,6),A1XA,A1YA,A1ZA,
     +     A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,1,NBLOCK,MAXB)
      IF(ISTRES >= 1)
     +     CALL MIRR(S11(1,1),S11(1,2),S11(1,3),S11(1,4),
     +     S11(1,5),S11(1,6),A1XA,A1YA,A1ZA,
     +     A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,1,NBLOCK,MAXEB)
      CALL CYCV(DM,DN,DW,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     + 1,NBLOCK,F1R,F1RM,VAR)
C ... CH IS NOT MIRRORED IN POST PROCESSING. PPR 29.5.96
C ... HAT4 is used temporary for pressure difference
      IF(IFSBC == 3) THEN
        CALL FRE(1,1,NBLOCK,-1.,0)
      ENDIF
      CALL CONEC(DM,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(DN,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(DW,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
C --- END OF VORTICITY LOOP -----------------------------------------
      CALL CONEC(P,    F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(PDIFF,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
C --- END OF VORTICITY LOOP -----------------------------------------

      CALL MIR(RO,RM,RN,RW,E,EPS2,VIST,P,FUN1,RK,REPS,PTUR,FI,PRO,VAR,
     +     TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     +     1,NBLOCK,1.,0,PRC)

      CALL CYCV(RM,RN,RW,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,1,
     + 1,NBLOCK,F1R,F1RM,VAR)
      IF(IFSBC == 3) THEN
       CALL FRE(1,1,NBLOCK,1.,0) 
      ENDIF

      CALL CONEC(RM,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(RN,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(RW,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)

      CALL CONEC(RO,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(E,    F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)

      CALL CONEC( EPS2,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC( VIST,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(    P,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(PDIFF,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(RLDIST,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)


C ... CONNECT THERMODYNAMIC VARIABLES

      CALL CONEC(TEMP,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1) ! no cyc and per
      CALL CONEC(   C,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC( VIS,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(  CH,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,1)
      CALL CONEC(TEMP,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,2) ! only cyc ==>
      CALL CONEC(   C,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,2) ! per is not con
      CALL CONEC( VIS,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,2)
      CALL CONEC(  CH,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,2)

C ... Turbulence quantities

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL CONEC(RK,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(REPS, F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(DE,   F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(DDEPS,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(FUN1, F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(TTS,  F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(PTUR, F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      ENDIF

C ... Intermittency variables
             
      IF (TRANSL) THEN
         IGT  = IG(1,NBLOCK+1)
         F1RW(1:IGT) = TRM(1:IGT)%G
         CALL CONEC(F1RW(1), F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         TRM(1:IGT)%G = F1RW(1:IGT)
         F1RW(1:IGT)  = TRM(1:IGT)%RET
         CALL CONEC(F1RW(1), F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         TRM(1:IGT)%RET = F1RW(1:IGT)
      ENDIF

C ... Multiphase quantities
            
      IF (MULPHL) THEN
         IGT  = IG(1,NBLOCK+1)
         DO IPHASE = 1,NPHASES
         F1RW(1:IGT) = PRO(1:IGT)%TEMP(IPHASE)
         CALL CONEC(F1RW(1), F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         PRO(1:IGT)%TEMP(IPHASE) = F1RW(1:IGT)
         F1RW(1:IGT) = PRO(1:IGT)%DTEMP(IPHASE)
         CALL CONEC(F1RW(1), F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         PRO(1:IGT)%DTEMP(IPHASE) = F1RW(1:IGT)
         F1RW(1:IGT) = VAR(1:IGT)%ALFA(IPHASE)
         CALL CONEC(F1RW(1), F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         VAR(1:IGT)%ALFA(IPHASE) = F1RW(1:IGT)
         F1RW(1:IGT) = VAR(1:IGT)%EVAPR(IPHASE)
         CALL CONEC(F1RW(1),F1R,F1RM,JLOC,JTRA,APP,
     &   1,1,NBLOCK,0,0)
         VAR(1:IGT)%EVAPR(IPHASE) = F1RW(1:IGT)
         IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
         CALL CONEC(VAR(1:IGT)%U(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(VAR(1:IGT)%V(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(VAR(1:IGT)%W(IPHASE),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FITOT(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FITOT(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FITOT(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FD(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FD(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FD(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FL(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FL(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FL(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FW(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FW(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FW(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FVM(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FVM(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FVM(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FTD(1),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FTD(2),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         CALL CONEC(PRO(1:IGT)%FTD(3),F1R,F1RM,JLOC,JTRA,APP,1,1,
     +   NBLOCK,0,1)
         ENDIF

         ENDDO
      ENDIF
         
C ... SCALAR EQUATIONS

      IF(NSCAL >= 1) THEN
         DO NS = 1,NSCAL
            CALL CONEC(FI(1,NS),F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         ENDDO
      ENDIF

C ... REYNOLDS STRESSES

      DO IHP = 1,2 ! MUST BE DONE BECAUSE TOLD ARE NOT OLD VARIABLES
C ... HAT3 on nykyaan vw 30.12.99
C --- MANU

*         CALL MIRR(TOLD,UOLD,VOLD,WOLD,HAT3,POLD,
*     +        A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
*     +        1,1,NBLOCK,MAXB)

         CALL MIRR(RS1,RS2,RS3,RS4,RS5,RS6,
     +        A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +        1,1,NBLOCK,MAXB)

      IF (ITURB >= 10 .OR. ISTRES == 1) THEN
         CALL CONEC(RS1,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(RS2,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(RS3,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(RS4,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(RS5,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         CALL CONEC(RS6,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      ENDIF

      ENDDO ! IHP = 1,2

      IF(STRESL) THEN
         DO NS = 1,6
            CALL CONEC(S11(1,NS),F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
         ENDDO
      ENDIF
              
C ... Update sliding surfaces 

      CALL SLD(U,V,W,EPS2,VIST,RO,TEMP,P,PDIFF,PTUR,E,RM,RN,RW,RK,REPS,
     +     FI,S11,BIJ,TIJ,PRO,VAR,TRM,F1R,F1RM,F1RN,MAXB,MAXEB,MAXSB,
     +     MAXSS,XCO,YCO,ZCO,1,1,NBLOCK,ISTRES,D1,D2,D3,VOL,DM,DN,DW,2)

      CALL PERNEC(RO,TEMP,U,V,W,E,P,PDIFF,RK,REPS,PRO,VAR,TRM,XC,YC,ZC,
     +     MAXB,1,1,NBLOCK,VIST,EPS2)

C ... REFLECTIONS ONCE MORE. BLOCK LOOP 2 BEGINS

      DO  998 N = 1,NBLOCK

      IG1 = IG(1,N)
      IF1 = JF(1,N)
      IC1 = IC(1,N)
      IR1 = IR(1,N)
      IQ1 = IQ(1,N)
      IGN = IG1 + NTOT(1,N) - 1

      CALL VIEWTE(TEMP(IG1),TWALL(1),IHF(1,1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Temperature
         CALL PDTOP(P(IG1),PDIFF(IG1),XC(IG1),YC(IG1),ZC(IG1),FRSDEN,
     +     FRSPRE,NTOT(1,N))

      IF(MULPHL) THEN
      DO IPHASE = 1,NPHASES
      F1RW(IG1:IGN) = PRO(IG1:IGN)%TEMP(IPHASE)
      CALL VIEWTE(F1RW(IG1),TWALL(1),IHF(1,1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Phase temperature
      PRO(IG1:IGN)%TEMP(IPHASE) = F1RW(IG1:IGN)
      F1RW(IG1:IGN) = PRO(IG1:IGN)%DTEMP(IPHASE) ! Voi mersu
      CALL VIEWTE(F1RW(IG1),TWALL(1),IHF(1,1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Phase temperature
      PRO(IG1:IGN)%DTEMP(IPHASE) = F1RW(IG1:IGN)
      ENDDO
      ENDIF

      CALL VIEWRF(DM(IG1),DN(IG1),DW(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     2 IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),KTOP(1,N),
     3 IN,JN,KN,ICON(IC1),NPATCH(N))

C ... THIS IS FOR SCALARS PPR 7.3

      IF (NSCAL >= 1) THEN
      CALL VIEWSC(FI(IG1,1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),ITURB,MAXSB,NSCAL,0)
      ENDIF

      IF (STRESL) THEN
      CALL VIEWSC(S11(IG1,1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),ITURB,MAXSS,6,0) 
      ENDIF
         
      CALL VIEWSC(OHMI(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Vorticity
      CALL VIEWSC(TIJ(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Lighthill's tensor gradient
      CALL VIEWSC(STRAIN(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,2) ! Shear strain

      CALL VIEWSC(RS1(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! uu
      CALL VIEWSC(RS2(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! uv
      CALL VIEWSC(RS3(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! uw
      CALL VIEWSC(RS4(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! vv
      CALL VIEWSC(RS5(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! vw
      CALL VIEWSC(RS6(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0) ! ww, 6 Reynolds stresses!

      CALL VIEWSC(HAT3(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),10,MAXB,1,0)
      CALL VIEWSC(FUN1(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),3,MAXB,1,3)
      CALL VIEWSC(TTS(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),3,MAXB,1,3)
      CALL VIEWSC(RLDIST(IG1),IMAX(1,N),
     2 JMAX(1,N),KMAX(1,N),IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     3 KBOT(1,N),KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     4 ICON(IC1),NPATCH(N),3,MAXB,1,3)
 998  CONTINUE
           
C ... TEMPin tilalla oli VIST
      CALL MIR(HAT1,RM,RN,RW,HAT2,EPS2,TEMP,P,PDIFF,RK,REPS,TTS,FI,PRO,
     +     VAR,TRM,MAXSB,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     +     1,1,NBLOCK,1.,0,PRC)

C ... HAT1 is production of turbulence for eddy-viscosity models
C ... HAT2 is viscous dissipation

      CALL CONEC(HAT1,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(HAT2,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)
      CALL CONEC(TIJ,F1R,F1RM,JLOC,JTRA,APP,1,1,NBLOCK,0,0)

      DO N = 1,NBLOCK
     
         IG1 = IG(1,N)
         IC1 = IC(1,N)

         CALL VIEWSC(HAT1(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &     IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),
     &     KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     &     ICON(IC1),NPATCH(N),3,MAXB,1,0)
         CALL VIEWSC(HAT2(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &     IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),
     &     KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     &     ICON(IC1),NPATCH(N),3,MAXB,1,0)

      ENDDO

C ... WRITE THE COMPUT-FILE. BLOCK LOOP 3 BEGINS

      IF(IOLD /= 2 .AND. IOLD /= -1 .AND. LEVEL > 1) THEN

C ... Open a new COMPUT file

      CALL NUMCH1(CLEVW,LEVEL)

      IF(IPRO == 1) THEN
      OPEN(15,FILE='COMPUT.'//CLEVW,STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(45,*) 
      WRITE(45,*) 'COMPUT file opened on level',LEVEL
      ENDIF

C ... In parallel simulation only the master process writes the COMPUT file.

      IF(PARALLEL) THEN

         CALL COMPIW(NBLOCK,IG,JF,IR,IMAX,JMAX,KMAX,IN,JN,KN,
     +   MGM,TEMP,U,V,W,PDIFF,P,RK,REPS,FI,MAXSB,NSCAL,ITURB,
     +   IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,0,15,MULPHL,BLKS,
     +   F1RW,VAR,PRO,TRANSL,TRM)

      ELSE

      DO 1980 N = 1,NBLOCK

         NGL = NPROCE(1+N,IPRO) !Global block number
         IG1 = IG(1,N)
         IF1 = JF(1,N)
         IC1 = IC(1,N)
         IR1 = IR(1,N)
         IGN = IG1 + NTOT(1,N) - 1 ! lis. 7.7.2007

         
      CALL WRITER(15,0,TEMP(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)
      CALL WRITER(15,0,U(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)
      CALL WRITER(15,0,V(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)
      CALL WRITER(15,0,W(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)
      CALL WRITER(15,0,PDIFF(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)
C ... turha
      CALL WRITER(15,0,P(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN)

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         CALL WRITER(15,0,RK(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +        IN,JN,KN)
         CALL WRITER(15,0,REPS(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +        IN,JN,KN)
      ENDIF  ! ITURB >= 3

C ... For intermittency variables

      IF(TRANSL) THEN
c         WRITE(3,*) 'writeriin meno'
         F1R(IG1:IGN)  = TRM(IG1:IGN)%G
         F1RM(IG1:IGN) = TRM(IG1:IGN)%RET
         CALL WRITER(15,0,F1R(IG1), ZZZ,IMAX(1,N),
     +   JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL WRITER(15,0,F1RM(IG1),ZZZ,IMAX(1,N),
     +   JMAX(1,N),KMAX(1,N),IN,JN,KN)
      ENDIF ! TRANSL

C     ... FOR SCALAR EQ. PPR 1.3.94

      DO NS = 1,NSCAL
         CALL WRITER(15,0,FI(IG1,NS),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +   IN,JN,KN)
      ENDDO

C     ... For multiphase modeling (not in COMPIW)

      IF(MULPHL) THEN

      DO IPHASE = 1,BLKS(NGL)%NPHASE

         ICHARP    = BLKS(NGL)%ICHAR(IPHASE)
c         WRITE(3,*) 'writeriin meno'
         F1R(IG1:IGN)  = PRO(IG1:IGN)%DTEMP(IPHASE)
         F1RM(IG1:IGN) = VAR(IG1:IGN)%X(IPHASE)

         IF(SIMULATION_END) THEN 

         CALL PRINYS(3,CHAR_PH(2,ICHARP),F1R(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         CALL PRINYS(3,CHAR_VAR(2,ICHARP),F1RM(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)

         ENDIF

         CALL WRITER(15,0,F1R(IG1), ZZZ,IMAX(1,N),
     +   JMAX(1,N),KMAX(1,N),IN,JN,KN)
         CALL WRITER(15,0,F1RM(IG1),ZZZ,IMAX(1,N),
     +   JMAX(1,N),KMAX(1,N),IN,JN,KN)

         IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN

            F1RM(IG1:IGN) = VAR(IG1:IGN)%U(IPHASE)
            F1RN(IG1:IGN) = VAR(IG1:IGN)%V(IPHASE)
            F1RW(IG1:IGN) = VAR(IG1:IGN)%W(IPHASE)
            CALL WRITER(15,0,F1RM(IG1), ZZZ,IMAX(1,N),
     +           JMAX(1,N),KMAX(1,N),IN,JN,KN)
            CALL WRITER(15,0,F1RN(IG1), ZZZ,IMAX(1,N),
     +           JMAX(1,N),KMAX(1,N),IN,JN,KN)
            CALL WRITER(15,0,F1RW(IG1), ZZZ,IMAX(1,N),
     +           JMAX(1,N),KMAX(1,N),IN,JN,KN)

         ENDIF                  ! MULTI

      ENDDO

      ENDIF                     ! MULPHL

 1980 CONTINUE
          
C ... END OF BLOCK LOOP 3                                             *

      ENDIF ! PARALLEL
      ENDIF ! IOLD /= 2 .AND. IOLD /= -1 .AND. LEVEL > 1
      
      IF(ITURB >= 3 .AND. ITURB /= 8 .AND. ITURB /= 9) THEN
         DO I = 1,MAXB          ! save PTUR/REPS in F1R for printing
            F1R(I) = PTUR(I)/(REPS(I) + DDEPS(I) + 1.E-8)
         ENDDO
      ELSEIF(ITURB == 9) THEN
         DO I = 1,MAXB          ! save PTUR/SRK in F1R for printing
            F1R(I) = PTUR(I)/(SRK(I) + 1.E-8)
         ENDDO
      ELSE
         DO I = 1,MAXB          ! no PTUR/REPS with algebraic models
            F1R(I) = 1.
         ENDDO
      ENDIF
           
C **********************************************************************
C ... Calculate turbulent dissipation, production and viscous dissipation
C **********************************************************************

      VISAVE  = 0.
      PROAVE  = 0.
      EPSAVE  = 0.
      EPWAVE  = 0.
        
      DO N = 1,NBLOCK
         NGL     = NPROCE(1+N,IPRO)  ! Global block number
         NTOTIN  = IMAX(1,N)*JMAX(1,N)*KMAX(1,N)
         IG1     = IG(1,N)
         IF1     = JF(1,N)
         IR1     = IR(1,N)
         IC1     = IC(1,N)

C ... CALCULATE THE AVERAGE STATE IN EACH BLOCK

C ... HAT1 is production of turbulence for algebraic models
C ... HAT2 is viscous dissipation

      CALL SUMFU(SUV,HAT2(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)

      IF(ITURB < 3) THEN ! Algebraic model

      CALL SUMFU(SUPRO,HAT1(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      SUEPS = SUPRO
      SUEPW = 0.

      ELSE 

      CALL SUMFU(SUPRO,PTUR(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),0.)
      IF(ITURB == 6) THEN ! k-omega
           CALL SUEFU(SUEPS, REPS(IG1),RK(IG1),RO(IG1),VOL(IG1),
     +     IMAX(1,N),JMAX(1,N),KMAX(1,N),EPSLIM)
      ELSE ! k-epsilon
           CALL SUMFU(SUEPS,REPS(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     +     KMAX(1,N),EPSLIM)
      ENDIF
      CALL SUMFU(SUEPW,DDEPS(IG1),VOL(IG1),IMAX(1,N),JMAX(1,N),
     + KMAX(1,N),0.)

      ENDIF ! Turbulence model type

      VISAVE  = VISAVE + SUV
      PROAVE  = PROAVE + SUPRO
      EPSAVE  = EPSAVE + SUEPS
      EPWAVE  = EPWAVE + SUEPW
          

      IF(SIMULATION_END) THEN  ! Write only at the end of simulation

      WRITE(4,*)
      WRITE(4,*) 'BLOCK = ',NGL
      WRITE(4,*)
      WRITE(4,*)'------------- INTEGRATED DISSIPATIONS ----------------'
      WRITE(4,*)
      WRITE(4,*) 'Viscous dissipation = ',REAL(SUV,4),  ' W'
      WRITE(4,*) 'Production power    = ',REAL(SUPRO,4),' W'
      WRITE(4,*) 'Dissipation tilde   = ',REAL(SUEPS,4),' W'
      WRITE(4,*) 'Dissipation wall    = ',REAL(SUEPW,4),' W'
      WRITE(4,*) 'Turb. dissipation   = ',REAL((SUEPS+SUEPW),4),' W'
      WRITE(4,*) 'Total dissipation   = ',REAL((SUEPS+SUEPW+SUV),4),' W'
      WRITE(4,*) 'Ratios: P/epsilon   = ',
     +            REAL(SUPRO/(SUEPS+SUEPW+1.E-20),4)
      WRITE(4,*) 'viscous/turbu diss. = ',
     +            REAL(SUV/(SUEPS+SUEPW+1.E-20),4)

      ENDIF  ! Write only at the end of simulation

      ENDDO  ! dissipations

      IF(SIMULATION_END) THEN  ! Write only at the end of simulation

      WRITE(4,*)
      WRITE(4,*)'------------- TOTAL DISSIPATIONS --------------------'
      WRITE(4,*)
      WRITE(4,*) 'Viscous dissipation = ',REAL(VISAVE,4), ' W'
      WRITE(4,*) 'Production power    = ',REAL(PROAVE,4),' W'
      WRITE(4,*) 'Dissipation tilde   = ',REAL(EPSAVE,4),' W'
      WRITE(4,*) 'Dissipation wall    = ',REAL(EPWAVE,4),' W'
      WRITE(4,*) 'Turb. dissipation   = ',REAL((EPSAVE+EPWAVE),4),' W'
      WRITE(4,*) 'Total dissipation   = ',
     +            REAL((EPSAVE+EPWAVE+VISAVE),4),' W'
      WRITE(4,*) 'Ratios: P/epsilon   = ',
     +            REAL(PROAVE/(EPSAVE+EPWAVE+1.E-20),4)
      WRITE(4,*) 'Viscous/turbu diss. = ',
     +            REAL(VISAVE/(EPSAVE+EPWAVE+1.E-20),4)
      IF(ABS(TOMTOT) >= 1.E-5) THEN
         WRITE(4,*) 'TOMEGA/T. dissi.    = ',REAL(-TOMTOT/
     +        (EPSAVE+EPWAVE+VISAVE+1.E-20),4)
      ENDIF

         IF(ABS(TOMTOT) >= 1.E-5) THEN
            WRITE(4,*) 'TOMEGA/T. dissi.    = ',REAL(-TOMTOT/
     &           (EPSAVE+EPWAVE+VISAVE+1.E-20),4)
         ENDIF

         IF(PARALLEL) CALL TOTDIS(VISAVE,PROAVE,EPSAVE,EPWAVE,
     &                            ITURB,IPRO,TOMTOG)

      ENDIF  ! SIMULATION_END


C ... Some surface properties before data is stored (y+,

      F1RK = 0.
       
      CALL FVPRO
       
C **********************************************************************
C ... Write a function file for visualization of the surface variables
C ... This must be done before changing to cell-vertex format ('REDUCE')
C ... Surface patch distributions: XYZBND.FMT, RESBND.FMT etc.
C **********************************************************************
         
      IF(FRESUL) SURLE(1:IBF) = WAVEH(1:IBF)         

      CALL FVSRF(IN,JN,KN,SIMULATION_END)
      
C**********************************************************************

C ... Collect and write post-processing Plot3D files XYZ.BIN, 
C ... XYZ_orig.BIN, RES.BIN and RES.AVE. 

C***********************************************************************

      IF(MASTER .AND. ENSIGHT) THEN

         IF(SIMULATION_END) THEN
            OPEN(1130,FILE='finflo.case',FORM='FORMATTED')
            WRITE(1130,'(A)') '# EnSight Gold case file from FINFLO run' 
            WRITE(1130,'(A)') '' 
            WRITE(1130,'(A)') 'FORMAT' 
            WRITE(1130,'(A)') 'type: ensight gold' 
            WRITE(1130,'(A)') ''
            WRITE(1130,'(A)') 'GEOMETRY'
         ENDIF

         IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
            OPEN(1130,FILE='MOVIE/finflo.case',FORM='FORMATTED')
            WRITE(1130,'(A)') '# EnSight Gold case file from FINFLO run' 
            WRITE(1130,'(A)') '' 
            WRITE(1130,'(A)') 'FORMAT' 
            WRITE(1130,'(A)') 'type: ensight gold' 
            WRITE(1130,'(A)') '' 
            WRITE(1130,'(A)') 'GEOMETRY' 
         ENDIF

      ENDIF

      CALL WRIP3D(SIMULATION_END,FIRSTOUTP)
      
C**********************************************************************

C ... Calculate true epsilon (or omega) and reduce

C***********************************************************************

      CALL OUTTUR

C***********************************************************************

C ... Reduce scalars to corner point format and write the function file

C***********************************************************************

C      CALL OUTSCL('SCL',FI) - obsolete

C***********************************************************************

C ... Write a function file for visualization the rest of the variables

C***********************************************************************

      CALL FVFUN(SIMULATION_END,FIRSTOUTP)

      IF(MASTER .AND. ENSIGHT .AND. SIMULATION_END) CLOSE(1130) 

C ... Save the output times of MOVIE files
      
      IF(MASTER .AND. .NOT.SIMULATION_END) THEN

         IUNIT = 287
         OPEN(IUNIT,FILE='MOVIE/OUTP_TIMES',STATUS='UNKNOWN')

         ITIMES_DUMMY = 0
         IF(.NOT.TIMEL) ITIMES = ICYCLE

         DO WHILE (ITIMES_DUMMY < ITIMES)
            READ(IUNIT,*,END=10,ERR=20) T_DUMMY,ITIMES_DUMMY
         ENDDO

 10      CONTINUE

         BACKSPACE IUNIT

 20      CONTINUE

         IF(TIMEL) THEN
            WRITE(IUNIT,'(F16.10,3X,I8)') T,ITIMES
         ELSE
            WRITE(IUNIT,'(F16.10,3X,I8)') REAL(ICYCLE,4),ICYCLE
         ENDIF
            
         CLOSE(IUNIT)

      ENDIF


C ... Final touch of the EnSight case file if needed

C ... Did we create MOVIE files ?

      IF(PARALLEL) THEN 
         CALL MPI_ALLREDUCE(NMOV,NMOVG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)
      ELSE
         NMOVG = NMOV
      ENDIF

C ... If yes ...
      
      IF(MASTER.AND.ENSIGHT.AND.NMOVG >= 1.AND.SIMULATION_END) THEN

         IUNIT = 287
         OPEN(IUNIT,FILE='MOVIE/OUTP_TIMES',STATUS='OLD')

         OPEN(1130,FILE='MOVIE/finflo.case',FORM='FORMATTED')

         NOFSTEPS = 0

C ... Read the time steps
         
 5       CONTINUE
         READ(IUNIT,*,END=11) T_DUMMY,ITIMES_DUMMY
         NOFSTEPS = NOFSTEPS + 1
         IF(NOFSTEPS == 1) ISTART = ITIMES_DUMMY
         IF(NOFSTEPS == 2) INCREMENT = ITIMES_DUMMY - ISTART
         GOTO 5

 11      CONTINUE
         REWIND IUNIT

C ... Scroll the MOVIE/finflo.case file to the end
         
 6       CONTINUE
         READ(1130,'(A)',END=13)
         GOTO 6
 13      CONTINUE
         BACKSPACE 1130
         
C ... Write the time step history to the end of the MOVIE/finflo.case file
         
         WRITE(1130,'(A)') ''
         WRITE(1130,'(A)')  'TIME'
         WRITE(1130,'(A)')  'time set: 1'
         WRITE(1130,'(A,I6)') 'number of steps:',NOFSTEPS
         WRITE(1130,'(A,I6)') 'filename start number:',ISTART
         WRITE(1130,'(A,I6)') 'filename increment:',INCREMENT
         WRITE(1130,'(A)')  'time values:'

         DO N=1,NOFSTEPS
            READ(IUNIT,*) T_DUMMY,ITIMES_DUMMY
            WRITE(1130,'(F16.6)') T_DUMMY
         ENDDO

         CLOSE(IUNIT)
         CLOSE(1130)
     
      ENDIF

       
C***********************************************************************

C ... PRINTING SOURCE TERMS.. TEMPORARY - obsolete
C
C      IF (ITURB >= 21) THEN
C         NSCAL = 6
C         CALL OUTSCL('PRO',PROD)
C         CALL OUTSCL('SPI',SPI)
C         CALL OUTSCL('DIS',DIS)
C         CALL OUTSCL('DIF',DIF)
C         CALL OUTSCL('VVS',VVIS)
C      ENDIF

      DEALLOCATE(RS1,RS2,RS3,RS4,RS5,RS6)
      DEALLOCATE(CVV1,CVV2,CVV3,CVV4,CVV5,CVV6,CVV7,CVV8,CVV9)
      DEALLOCATE(GVX,GVY,GVZ,VISTND)
      DEALLOCATE(UREL,VREL,WREL,RMREL,RNREL,RWREL,EREL)
      DEALLOCATE(MPUREL,MPVREL,MPWREL)

      RETURN
      END SUBROUTINE OUTP
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRIP3D(SIMULATION_END,FIRSTOUTP)

C ... Collect (in parallel run) and write flow field for post-processing.
C ... In Plot3D format XYZ.BIN, XYZ_orig.BIN, RES.BIN and RES.AVE.
C ... In EnSight Gold format XYZ.BIN.geo etc.

      USE MPI

      USE CHARACTERS

      USE CONSTANTS, ONLY : EPS

      USE INTEGERS,  ONLY : IPRO

      USE MAIN_ARRAYS

      USE NS3CO


      IMPLICIT NONE

      INTEGER :: N, NGL, IG1, IF1, IC1, IR1, IT1, IL1, IK1, IPT, 
     &           IM1, JM1, KM1, NS, INSCA, NMAX, I, J, K, L, IGN,
     &           IMAXP1, JMAXP1, KMAXP1, ISTRID, JSTRID, ILE, INODE,
     &           IPHASE, NROOT, NR1, IERR, II, NCELLS, IMX, JMX, KMX,
     &           NTRAIN, NSPOINTS, IPB, LL, MM, NG, IG3, IG4

      INTEGER STATUS(MPI_STATUS_SIZE)

      LOGICAL :: MASTER, SLAVE, SIMULATION_END, FIRSTOUTP

      CHARACTER(LEN=35) :: VNAME
      CHARACTER(LEN=80) :: FNAME
      CHARACTER(LEN=6)  :: FILEE
      CHARACTER(LEN=80) :: COMMAND
      
      REAL, ALLOCATABLE, DIMENSION(:) :: TRAINW, BLANKW

      MASTER = IPRO == 1 
      SLAVE = .NOT.MASTER    
      INSCA  = 1
      
      IF(TIMEL) THEN
         CALL NUMCH6(FILEE,ITIMES)
      ELSE
         CALL NUMCH6(FILEE,ICYCLE)
      ENDIF

      
      IF(ITURB >= 6 .AND. ITURB /= 8) INSCA = 7


      DO 2000 N = 1,NBLOCK

      NGL     = NPROCE(1+N,IPRO)  ! Global block number
      IG1     = IG(1,N)
      IF1     = JF(1,N)
      IC1     = IC(1,N)
      IR1     = IR(1,N)

      IF (COORL) THEN
         IG3  = IG1
      ELSE
         IG3 = 1
      ENDIF

      IF (MULPHL) THEN
         IG4  = IG1
      ELSE
         IG4 = 1
      ENDIF

      IGN = IG1 + NTOT(1,N) - 1 ! lis. 5.1.2010
 
C **********************************************************************
C *** FROM THE ABSOLUTE TO THE RELATIVE COORDINATES ********************

      IT1  = IT(1,N)
      IL1  = IL(1,N)
      IK1  = IK(1,N)
      IM1  = IMAX(1,N)
      JM1  = JMAX(1,N)
      KM1  = KMAX(1,N)

      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation
         WRITE(3,*)
         WRITE(3,*) ' BEFORE THE ABSOLUTE TO THE RELATIVE COORDINATES'
         CALL PRINYS(3,   UC,   U(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      ENDIF  ! Write only at the end of the simulation

C ... Initialize the relative variables
      
      UREL(IG1:IGN) = U(IG1:IGN)
      VREL(IG1:IGN) = V(IG1:IGN)
      WREL(IG1:IGN) = W(IG1:IGN)
      
      RMREL(IG1:IGN) = RM(IG1:IGN)
      RNREL(IG1:IGN) = RN(IG1:IGN)
      RWREL(IG1:IGN) = RW(IG1:IGN)
       EREL(IG1:IGN) =  E(IG1:IGN)

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
            MPUREL(IPHASE,IG1:IGN) = VAR(IG1:IGN)%U(IPHASE)
            MPVREL(IPHASE,IG1:IGN) = VAR(IG1:IGN)%V(IPHASE)
            MPWREL(IPHASE,IG1:IGN) = VAR(IG1:IGN)%W(IPHASE)
         ENDDO
      ENDIF

      
      IF(IROTVE(N) > 0) THEN
        CALL CARROT(UREL(IG1),VREL(IG1),WREL(IG1),
     +  XC(IG1),YC(IG1),ZC(IG1),
     +  OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),CENAX(NGL),
     +  CENAY(NGL),CENAZ(NGL),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +  VAR(IG4),MPUREL(1,IG4),MPVREL(1,IG4),MPWREL(1,IG4),
     +  BLKS(NGL)%SOLUTION_TYPE,NPHASE,NGL)
      ENDIF

      IF(IROTVE(N) == 2) THEN
        CALL RMROT(RO(IG1),RMREL(IG1),RNREL(IG1),RWREL(IG1),
     +  UREL(IG1),VREL(IG1),WREL(IG1),EREL(IG1),
     +  IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     +  VAR(IG4),PRO(IG4),BLKS(NGL)%SOLUTION_TYPE,NPHASE)
      ENDIF

      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation

      CALL PRINYS(3,   PC,    P(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,PDIFC,PDIFF(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,TEMPC, TEMP(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      WRITE(3,*)
      WRITE(3,*) ' AFTER THE ABSOLUTE TO THE RELATIVE COORDINATES'
      CALL PRINYS(3,  RMC,RMREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   UC, UREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   CC,    C(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  RNC,RNREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   VC, VREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  RWC,RWREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   WC, WREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3,  ROC,   RO(IG1),IT1,IL1,0,  0,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  ROC,   RO(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   EC, EREL(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3, DROC,  DRO(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3, VISC,  VIS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,EPS2C, EPS2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,VISTC, VIST(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
 
      IF(TWO_FLUIDL) THEN
        F1R(IG1:IGN) = PRO(IG1:IGN)%FL(1)
         CALL PRINYS(3,FLXC,F1R(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         F1R(IG1:IGN) = PRO(IG1:IGN)%FL(2)
         CALL PRINYS(3,FLYC,F1R(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
         F1R(IG1:IGN) = PRO(IG1:IGN)%FL(3)
         CALL PRINYS(3,FLZC,F1R(IG1),
     +   IT(1,N),IL(1,N),0,IK(1,N),IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
      ENDIF
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL PRINYS(3,  RKC,   RK(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3,  DKC,  DRK(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3,   DC,DDEPS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(ITURB /= 8 .OR. ITURB /= 9)
     +CALL PRINYS(3,REPSC, REPS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(ITURB == 9)
     +CALL PRINYS(3,RNUTC, REPS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3,DEPSC, DEPS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      IF(N <= NBLOCK)
     +CALL PRINYS(3,STRAIC,STRAIN(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,PTURC, PTUR(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,PTUREC, F1R(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,FUN1C, FUN1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,TTSC,   TTS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,TIJC,   TIJ(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      ENDIF

      IF (ITURB == 6) THEN ! .AND. KOVER == 3) THEN
        CALL PRINYS(3,PTURC,VELLAP(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
        CALL PRINYS(3,PTURC,  QSAS(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      ENDIF

C ... SCALAR FUNCTIONS PPR 23.3
      DO 1930 NS = INSCA,NSCAL
      CALL PRINYS(3,FIC(NS),FI(IG1,NS),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
 1930 CONTINUE

      IF (STRESL) THEN ! Muuta tama ed. kaltaiseksi
c ... B_ij can be printed here
      DO 1933 NS = 1,6
c      IG3     = IG1 + (NS-1)*MAXSS
      CALL PRINYS(3,S11C(NS),BIJ(IG1,NS),IT1,IL1,0,IK1,IM1,JM1,KM1,
     +  NGL,1)
1933  CONTINUE
      DO 1931 NS = 1,6
c       IG3     = IG1 + (NS-1)*MAXSS
      CALL PRINYS(3,S11C(NS),S11(IG1,NS),IT1,IL1,0,IK1,IM1,JM1,KM1,
     +  NGL,1)
 1931 CONTINUE
      ENDIF

C ... Reynolds stresses
           
*      CALL PRINYS(3,UUC(1),TOLD(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
*      CALL PRINYS(3,UUC(2),UOLD(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
*      CALL PRINYS(3,UUC(3),VOLD(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
*      CALL PRINYS(3,UUC(4),WOLD(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
*      CALL PRINYS(3,UUC(5),HAT3(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
*      CALL PRINYS(3,UUC(6),POLD(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
           
      CALL PRINYS(3,UUC(1),RS1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(2),RS2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(3),RS3(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(4),RS4(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(5),RS5(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(6),RS6(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)

      ENDIF  ! Write only at the end of the simulation
      
      
      IF(TIMEL) THEN ! Calculate Reynolds stresses from fluctuations

      NMAX = (IM1+2*IN)*(JM1+2*JN)*(KM1+2*KN)

      DO I = IG1,IG1+NMAX-1 !  change quadratic to reynolds stresses
         RMAV2(I) = RMAV2(I) - RMAV1(I)**2/(ROAV1(I)+EPS)
         RNAV2(I) = RNAV2(I) - RNAV1(I)**2/(ROAV1(I)+EPS)
         RWAV2(I) = RWAV2(I) - RWAV1(I)**2/(ROAV1(I)+EPS)
         RUVAV(I) = RUVAV(I) - RMAV1(I)*RNAV1(I)/(ROAV1(I)+EPS)
         RUWAV(I) = RUWAV(I) - RMAV1(I)*RWAV1(I)/(ROAV1(I)+EPS)
         RVWAV(I) = RVWAV(I) - RNAV1(I)*RWAV1(I)/(ROAV1(I)+EPS)
      ENDDO


      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation
      
      WRITE(3,*)
      WRITE(3,*) '  Time averaged values for timeperiod and steps',
     +              REAL(TAVER,4),NAVER

      CALL PRINYS(3,  ROC,ROAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  RMC,RMAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  RNC,RNAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,  RWC,RWAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   EC, EAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,   PC, PAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,TEMPC, TAV1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(1),RMAV2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(4),RNAV2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(6),RWAV2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(2),RUVAV(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(3),RUWAV(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,UUC(5),RVWAV(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)

      ENDIF  ! Write only at the end of the simulation

      ENDIF  ! TIMEL


      IF(SIMULATION_END) THEN  ! Write only at the of the simulation 

      IF(ITURB < 3) THEN
      CALL PRINYS(3,PROKC,  HAT1(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      ELSE
      CALL PRINYS(3,PTURC,  PTUR(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      ENDIF
      CALL PRINYS(3,PROVC,  HAT2(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)

      ENDIF  ! Write only at the end of the simulation
   
      
      IMAXP1  = IMAX(1,N) + 1
      JMAXP1  = JMAX(1,N) + 1
      KMAXP1  = KMAX(1,N) + 1
      ISTRID  = IMAX(1,N) + 2*IN
      JSTRID  = JMAX(1,N) + 2*JN
      ILE     = ISTRID*JSTRID

************************************************************************

      IF(N <= NBLOCK) THEN
   
      IF(SIMULATION_END) THEN  ! Write only at the of the simulation 

      CALL PRINYS(3,OHMIC,OHMI(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,OHMIXC, DM(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,OHMIYC, DN(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,OHMIZC, DW(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,STRAINC,STRAIN(IG1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(1),SHEAR(IG1,1),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(2),SHEAR(IG1,2),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(3),SHEAR(IG1,3),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(4),SHEAR(IG1,4),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(5),SHEAR(IG1,5),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)
      CALL PRINYS(3,SRC(6),SHEAR(IG1,6),IT1,IL1,0,IK1,IM1,JM1,KM1,NGL,1)

      ENDIF  ! Write only at the end of the simulation 

      
C     **********************VISUALIZATION*******************************
C ... FROM CELL CENTRED TO CELL VERTEX SYSTEM                          *
C     ******************************************************************

*      CALL REDUCE(EPS2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
*      CALL REDUCE(RLDIST(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))

C ... Reduce grid velocities

      GVEX = BLKS(NGL)%GVEX
      GVEY = BLKS(NGL)%GVEY
      GVEZ = BLKS(NGL)%GVEZ

      CALL REDUCE_SRF(GVX(IG1),GVY(IG1),GVZ(IG1),XCO(IG1),YCO(IG1),
     + ZCO(IG1),XLE2(IG3),YLE2(IG3),ZLE2(IG3),
     + XLE3(IG3),YLE3(IG3),ZLE3(IG3),GVEX,GVEY,GVEZ,
     + OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),CENAX(NGL),
     + CENAY(NGL),CENAZ(NGL),DT,DTOLD,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     + IN,JN,KN,TIMEL,IGRID(NGL,1),NGL)


*      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
* 
*         DO NS = 1,6
*            CALL REDUCE(SHEAR(IG1,NS),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
*         ENDDO
*
*      ENDIF

C ... Needed at node points?
*      CALL REDUCE(HAT1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))

      ENDIF  ! N <= NBLOCK

2000  CONTINUE  ! End of a very long block loop

************************************************************************
      
      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation
         IF(PARALLEL .AND. MASTER) THEN
            WRITE(*,*)
            WRITE(*,*)' Collecting and writing XYZ.BIN  ...'
            IF(COORL) THEN
               WRITE(*,*)
               WRITE(*,*)' Collecting and writing XYZ_orig.BIN  ...'
            ENDIF
         ENDIF
      ENDIF  ! Write only at the end of the simulation

      DO NG = 1,NBLOCG

         N     = NLOCAL(NG)     ! local block number
         NROOT = NPNUM(NG)      ! process number
         NR1   = NROOT - 1      ! process number used by MPI
         
         NCELLS   = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
         NSPOINTS = (IMAXG(NG)+1)*(JMAXG(NG)+1)*(KMAXG(NG)+1)
         NTRAIN   = 4*NSPOINTS
         IF(COORL) NTRAIN = NTRAIN + 3*NSPOINTS

         IF(MASTER) ALLOCATE(TRAINW(NTRAIN)) 


         IF(N /= -1) THEN

            IG1 = IG(1,N)
            IC1 = IC(1,N)

            ALLOCATE(BLANKW(NCELLS)) 


            IF(NCHIM > 0) THEN  ! Add IBLANK vector into XYZ.BIN

C ... Trim the blanking vector

            BLANKW(1:NCELLS) = BLANK(IG1:IG1+NCELLS-1)

            CALL TBLANK(BLANKW,IMAX(1,N),JMAX(1,N),KMAX(1,N),LEVEL)

C ... Blanking from cell center points to grid points (cell vertex)

            CALL REDUCE(BLANKW,ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))

            NMAX = (IMAX(1,N) + 1)*(JMAX(1,N) + 1)*(KMAX(1,N) + 1)

            DO INODE = 1,NMAX
               IF(BLANKW(INODE) < 0.0) BLANKW(INODE) = 0.0
               IF(BLANKW(INODE) > 0.0) BLANKW(INODE) = 1.0
            ENDDO

C ... Mark the grid points on solid surfaces ('wall' points) in order 
C ... to avoid blanking them.

            CALL WALLGP(BLANKW,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &           NPATCH(N),ICON(IC1))

            ENDIF  ! NCHIM > 0
************************************************************************

            IF(SLAVE) ALLOCATE(TRAINW(NTRAIN)) 

            IMAXP1 = IMAX(1,N) + 1
            JMAXP1 = JMAX(1,N) + 1
            KMAXP1 = KMAX(1,N) + 1
            ISTRID = IMAX(1,N) + 2*IN
            JSTRID = JMAX(1,N) + 2*JN
            ILE    = ISTRID*JSTRID
            IPB    = 0
            
            DO K = 1,KMAXP1
               DO J = 1,JMAXP1
                  DO I = 1,IMAXP1

                     IPB = IPB + 1
                     MM  = IG1+I+IN-1 + (J+JN-1)*ISTRID + (K+KN-1)*ILE

                     TRAINW(IPB             ) = XCO(MM)      
                     TRAINW(IPB + 1*NSPOINTS) = YCO(MM)      
                     TRAINW(IPB + 2*NSPOINTS) = ZCO(MM)      
                     TRAINW(IPB + 3*NSPOINTS) = BLANKW(IPB)

                     IF(COORL) THEN
                        TRAINW(IPB + 4*NSPOINTS) = XGRI(MM)      
                        TRAINW(IPB + 5*NSPOINTS) = YGRI(MM)      
                        TRAINW(IPB + 6*NSPOINTS) = ZGRI(MM)      
                     ENDIF
      
                  ENDDO            
               ENDDO            
            ENDDO            

            DEALLOCATE(BLANKW)
            
            IF(SLAVE) THEN
               CALL MPI_SEND(TRAINW,NTRAIN,MPI_REAL8,0,
     &              NG,MPI_COMM_WORLD,IERR)
               DEALLOCATE(TRAINW)
            ENDIF

         ELSEIF(MASTER .AND. N == -1) THEN

            CALL MPI_RECV(TRAINW,NTRAIN,MPI_REAL8,NR1,NG,
     &           MPI_COMM_WORLD,STATUS,IERR)

         ENDIF                  ! (N /= -1)

  
C ... Finally write the grid points in PLOT3D or EnSight Gold format 
C ... for visualization and other post-processing purposes. 

         IF(MASTER) THEN

            IF(SIMULATION_END) THEN
               FNAME='XYZ.BIN'
            ELSE
               FNAME='MOVIE/XYZ'//FILEE//'.BIN'
            ENDIF
            FNAME=TRIM(FNAME)

            IF(PLOT3D) CALL WritePlot3dGeometry(FNAME,27,
     &           NCHIM,NG,NBLOCG,IMAXG,JMAXG,KMAXG,
     &           TRAINW(1),TRAINW(NSPOINTS+1),TRAINW(2*NSPOINTS+1),
     &           TRAINW(3*NSPOINTS+1))

            IF(PLOT3D .AND. .NOT.SIMULATION_END) THEN
               FNAME='MOVIE/XYZ'//FILEE//'.BIN.fvbnd'
               FNAME=TRIM(FNAME)
               CALL CopyAsciiFile('XYZ.BIN.fvbnd',FNAME)
            ENDIF
            
            IF(SIMULATION_END) THEN
               FNAME='XYZ.BIN.geo'
            ELSE
               FNAME='MOVIE/XYZ'//FILEE//'.BIN.geo'
            ENDIF
            FNAME = TRIM(FNAME)

            IF(ENSIGHT) THEN

               CALL WriteEnSightGeometry(FNAME,1127,
     &              NCHIM,NG,NBLOCG,
     &              IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &              TRAINW(1),TRAINW(NSPOINTS+1),
     &              TRAINW(2*NSPOINTS+1),TRAINW(3*NSPOINTS+1))

               IF(NG == 1) THEN

                  IF(SIMULATION_END) THEN
                     WRITE(1130,'(A)')'model: XYZ.BIN.geo'
                     WRITE(1130,'(A)')'boundary: XYZ.BIN.ensbnd'
                  ENDIF

*                  IF(.NOT.SIMULATION_END) THEN
                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
*                     FNAME='MOVIE/XYZ'//FILEE//'.BIN.ensbnd'
                     FNAME='MOVIE/XYZ.BIN.ensbnd'
                     FNAME=TRIM(FNAME)
                     CALL CopyAsciiFile('XYZ.BIN.ensbnd',FNAME)
                  ENDIF

                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
                     WRITE(1130,'(A)')'model: XYZ******.BIN.geo'
*                     WRITE(1130,'(A)')'boundary: XYZ******.BIN.ensbnd'
                     WRITE(1130,'(A)')'boundary: XYZ.BIN.ensbnd'
                  ENDIF

               ENDIF

            ENDIF 

            IF(SIMULATION_END) THEN
               FNAME='XYZ_orig.BIN'
            ELSE
               FNAME='MOVIE/XYZ'//FILEE//'_orig.BIN'
            ENDIF
            FNAME = TRIM(FNAME)

            IF(COORL) THEN
               IF(PLOT3D) CALL WritePlot3dGeometry(FNAME,32,0,
     &              NG,NBLOCG,IMAXG,JMAXG,KMAXG,
     &              TRAINW(4*NSPOINTS+1),TRAINW(5*NSPOINTS+1),
     &              TRAINW(6*NSPOINTS+1),TRAINW(3*NSPOINTS+1))
            ENDIF

            DEALLOCATE(TRAINW)

         ENDIF 

      ENDDO                     ! Global block loop


************************************************************************

C ... XYZ.BIN has been written, start to store the results

C **********************************************************************

      IF(MASTER .AND. ENSIGHT) THEN

         IF(SIMULATION_END .OR. .NOT.SIMULATION_END.AND.FIRSTOUTP) THEN
            WRITE(1130,'(A)') ''
            WRITE(1130,'(A)') 'VARIABLE'
         ENDIF

      ENDIF

      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation
         IF(PARALLEL .AND. MASTER) THEN
            WRITE(*,*)
            WRITE(*,*)' Collecting and writing RES.BIN  ...'
         ENDIF
      ENDIF  ! Write only at the end of the simulation

      DO NG = 1,NBLOCG

         N     = NLOCAL(NG)     ! local block number
         NROOT = NPNUM(NG)      ! process number
         NR1   = NROOT - 1      ! process number used by MPI

         NCELLS = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
         NTRAIN = 5*NCELLS

         IF(MASTER) ALLOCATE(TRAINW(NTRAIN)) 


         IF(N /= -1) THEN

            IG1 = IG(1,N)

            IMX = IMAX(1,N)
            JMX = JMAX(1,N)
            KMX = KMAX(1,N)

            IF(SLAVE) ALLOCATE(TRAINW(NTRAIN))

            TRAINW(0*NCELLS+1:1*NCELLS) =    RO(IG1:IG1+NCELLS-1)
            TRAINW(1*NCELLS+1:2*NCELLS) = RMREL(IG1:IG1+NCELLS-1)
            TRAINW(2*NCELLS+1:3*NCELLS) = RNREL(IG1:IG1+NCELLS-1)
            TRAINW(3*NCELLS+1:4*NCELLS) = RWREL(IG1:IG1+NCELLS-1)
            TRAINW(4*NCELLS+1:5*NCELLS) =  EREL(IG1:IG1+NCELLS-1)

            DO LL = 0,4
               CALL REDUCE(TRAINW(LL*NCELLS+1),ZZZ,IMX,JMX,KMX)
            ENDDO

            IF(SLAVE) THEN
               CALL MPI_SEND(TRAINW,NTRAIN,MPI_REAL8,0,
     &              NG,MPI_COMM_WORLD,IERR)
               DEALLOCATE(TRAINW)
            ENDIF

         ELSEIF(MASTER .AND. N == -1) THEN

            CALL MPI_RECV(TRAINW,NTRAIN,MPI_REAL8,NR1,NG,
     &           MPI_COMM_WORLD,STATUS,IERR)

         ENDIF                  ! (N /= -1)


         IF(MASTER) THEN

            IF(SIMULATION_END) THEN
               FNAME='RES.BIN'
            ELSE
               FNAME='MOVIE/RES'//FILEE//'.BIN'
            ENDIF
            FNAME = TRIM(FNAME)

            IF(PLOT3D) CALL WritePlot3dSolution(FNAME,17,NCHIM,
     &           NG,NBLOCG,IMAXG,JMAXG,KMAXG,RMACH,ALPHA,RE,T,ICYCLE,
     &           TIMEL,TRAINW(1),TRAINW(NCELLS+1),
     &           TRAINW(2*NCELLS+1),TRAINW(3*NCELLS+1),
     &           TRAINW(4*NCELLS+1))

            IF(ENSIGHT) THEN

               IF(SIMULATION_END) THEN
                  FNAME = 'RES.BIN.rho'
               ELSE
                  FNAME = 'MOVIE/RES'//FILEE//'.BIN.rho'
               ENDIF
                  
               CALL WriteEnsightScalar(FNAME,1131,
     &              NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &              TRAINW(1))
 
               VNAME = 'Density'

               IF(NG == 1) THEN

                  IF(SIMULATION_END) THEN
                     WRITE(1130,'(A)')
     &               'scalar per node: '//VNAME//'RES.BIN.rho'
                  ENDIF

                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
                     WRITE(1130,'(A)')
     &               'scalar per node: '//VNAME//'RES******.BIN.rho'
                  ENDIF

               ENDIF

               
               IF(SIMULATION_END) THEN
                  FNAME = 'RES.BIN.E'
               ELSE
                  FNAME = 'MOVIE/RES'//FILEE//'.BIN.E'
               ENDIF                  
               
               CALL WriteEnsightScalar(FNAME,1132,
     &              NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &              TRAINW(4*NSPOINTS+1))

               VNAME = 'Total_E'

               IF(NG == 1) THEN

                  IF(SIMULATION_END) THEN
                     WRITE(1130,'(A)')
     &               'scalar per node: '//VNAME//'RES.BIN.E'
                  ENDIF

                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
                     WRITE(1130,'(A)')
     &               'scalar per node: '//VNAME//'RES******.BIN.E'
                  ENDIF

               ENDIF

 
* If you need both velocity and momentum in EnSight Gold format remove '*'s.
*
*               IF(SIMULATION_END) THEN
*                  FNAME = 'RES.BIN.mom'
*               ELSE
*                  FNAME = 'MOVIE/RES'//FILEE//'.BIN.mom'
*               ENDIF                  
*               
*               CALL WriteEnsightVector(FNAME,1133,
*     &              NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
*     &              TRAINW(1*NSPOINTS+1),
*     &              TRAINW(2*NSPOINTS+1),
*     &              TRAINW(3*NSPOINTS+1))
*
*               VNAME = 'Momentum'
*
*               IF(NG == 1) THEN
*
*                  IF(SIMULATION_END) THEN
*                     WRITE(1130,'(A)')
*     &               'vector per node: '//VNAME//'RES.BIN.mom'
*                  ENDIF
*
*                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
*                     WRITE(1130,'(A)')
*     &               'vector per node: '//VNAME//'RES******.BIN.mom'
*                  ENDIF
*
*               ENDIF

               
               L = 0
               DO K=1,KMAXG(NG)+1
                  DO J=1,JMAXG(NG)+1
                     DO I=1,IMAXG(NG)+1
                        L  = L + 1
                        F1RM(L) = TRAINW(  NSPOINTS+L)/TRAINW(L)
                        F1RN(L) = TRAINW(2*NSPOINTS+L)/TRAINW(L)  
                        F1RW(L) = TRAINW(3*NSPOINTS+L)/TRAINW(L)  
                     ENDDO
                  ENDDO
               ENDDO
 

               IF(SIMULATION_END) THEN
                  FNAME = 'RES.BIN.vel'
               ELSE
                  FNAME = 'MOVIE/RES'//FILEE//'.BIN.vel'
               ENDIF                  

               CALL WriteEnsightVector(FNAME,1134,
     &              NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &              F1RM,F1RN,F1RW)

               VNAME = 'Velocity'

               IF(NG == 1) THEN

                  IF(SIMULATION_END) THEN
                     WRITE(1130,'(A)')
     &               'vector per node: '//VNAME//'RES.BIN.vel'
                  ENDIF

                  IF(.NOT.SIMULATION_END .AND. FIRSTOUTP) THEN
                     WRITE(1130,'(A)')
     &               'vector per node: '//VNAME//'RES******.BIN.vel'
                  ENDIF

               ENDIF

            ENDIF

            DEALLOCATE(TRAINW)

         ENDIF 

      ENDDO                     ! Global block loop


C ... Time averaged values                         

      IF (TIMEL .AND. SIMULATION_END) THEN

         IF(PARALLEL .AND. MASTER) THEN
            WRITE(*,*)
            WRITE(*,*)' Collecting and writing RES.AVE  ...'
         ENDIF

         DO NG  = 1,NBLOCG

            N     = NLOCAL(NG)  ! local block number
            NROOT = NPNUM(NG)   ! process number
            NR1   = NROOT - 1   ! process number used by MPI

            NSPOINTS = (IMAXG(NG)+1)*(JMAXG(NG)+1)*(KMAXG(NG)+1)
            NTRAIN   = 5*NSPOINTS

            IF(MASTER) ALLOCATE(TRAINW(NTRAIN)) 

            IF(N /= -1) THEN

               IG1 = IG(1,N)

               CALL REDUCE(ROAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RMAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RNAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RWAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE( EAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE( PAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE( TAV1(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(ROAV2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RMAV2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RNAV2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RWAV2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE( EAV2(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RUVAV(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RUWAV(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
               CALL REDUCE(RVWAV(IG1),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))

               IF(SLAVE) ALLOCATE(TRAINW(NTRAIN)) 

               IMAXP1 = IMAX(1,N) + 1
               JMAXP1 = JMAX(1,N) + 1
               KMAXP1 = KMAX(1,N) + 1
               IPB    = 0

               DO K = 1,KMAXP1
               DO J = 1,JMAXP1
               DO I = 1,IMAXP1

                  IPB = IPB + 1
                  LL  = IG1-1+I+IMAXP1*(J-1)+IMAXP1*JMAXP1*(K-1)

                  TRAINW(IPB             ) = ROAV1(LL)
                  TRAINW(IPB + 1*NSPOINTS) = RMAV1(LL)
                  TRAINW(IPB + 2*NSPOINTS) = RNAV1(LL)
                  TRAINW(IPB + 3*NSPOINTS) = RWAV1(LL)
                  TRAINW(IPB + 4*NSPOINTS) =  EAV1(LL)
      
               ENDDO            
               ENDDO            
               ENDDO            

               IF(SLAVE) THEN
                  CALL MPI_SEND(TRAINW,NTRAIN,MPI_REAL8,0,
     &                          NG,MPI_COMM_WORLD,IERR)
                  DEALLOCATE(TRAINW)
               ENDIF

            ELSEIF(MASTER .AND. N == -1) THEN

               CALL MPI_RECV(TRAINW,NTRAIN,MPI_REAL8,NR1,NG,
     &         MPI_COMM_WORLD,STATUS,IERR)

            ENDIF ! (N /= -1)
           

            IF(MASTER) THEN

               IF(PLOT3D) CALL WritePlot3dSolution('RES.AVE',62,NCHIM,
     &              NG,NBLOCG,IMAXG,JMAXG,KMAXG,RMACH,ALPHA,RE,
     &              TAVER,ICYCLE,TIMEL,
     &              TRAINW(1),TRAINW(NSPOINTS+1),
     &              TRAINW(2*NSPOINTS+1),TRAINW(3*NSPOINTS+1),
     &              TRAINW(4*NSPOINTS+1))

               IF(ENSIGHT) THEN
 
                  VNAME = 'Density_average'
                  CALL WriteEnsightScalar('RES.AVE.rho',1135,
     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                 TRAINW(1))
                  IF(NG == 1) WRITE(1130,'(A)') 
     &               'scalar per node: '//VNAME//'RES.AVE.rho' 
 
                  VNAME = 'Total_E_average'
                  CALL WriteEnsightScalar('RES.AVE.E',1136,
     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                 TRAINW(4*NSPOINTS+1))
                  IF(NG == 1) WRITE(1130,'(A)') 
     &               'scalar per node: '//VNAME//'RES.AVE.E' 
 
* If you need the momentum vector average in EnSight Gold format at the
* end of the simulation remove '*'s.               
*
*                  VNAME = 'Momentum_average'
*                  CALL WriteEnsightVector('RES.AVE.mom',1137,
*     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
*     &                 TRAINW(  NSPOINTS+1),
*     &                 TRAINW(2*NSPOINTS+1),
*     &                 TRAINW(3*NSPOINTS+1))
*                  IF(NG == 1) WRITE(1130,'(A)') 
*     &               'vector per node: '//VNAME//'RES.AVE.mom' 

                  L = 0
                  DO K=1,KMAXG(NG)+1
                     DO J=1,JMAXG(NG)+1
                        DO I=1,IMAXG(NG)+1
                           L  = L + 1
                           F1RM(L) = TRAINW(  NSPOINTS+L)/TRAINW(L)
                           F1RN(L) = TRAINW(2*NSPOINTS+L)/TRAINW(L)  
                           F1RW(L) = TRAINW(3*NSPOINTS+L)/TRAINW(L)  
                        ENDDO
                     ENDDO
                  ENDDO
 
                  VNAME = 'Velocity_average'
                  CALL WriteEnsightVector('RES.AVE.vel',1138,
     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                 F1RM,F1RN,F1RW)

                  IF(NG == 1) WRITE(1130,'(A)') 
     &               'vector per node: '//VNAME//'RES.AVE.vel' 

               ENDIF

               DEALLOCATE(TRAINW)

            ENDIF 

         ENDDO  ! Global block loop
      
      ENDIF  ! TIMEL

      RETURN
      END SUBROUTINE WRIP3D
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRIXYD

      USE MPI

      USE CHARACTERS

      USE CONSTANTS, ONLY : PII,EPS
 
      USE INTEGERS,  ONLY : IPRO

      USE NS3CO,     ONLY : PARALLEL, NBLOCG, NBLOCK, ROTANG, ROTAT,
     &                      IN, JN, KN

      USE MAIN_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMEM,N,KOKO,KOKO1,I,J,K,II,JJ,IERRCODE,ERRORCODE,IERR,
     +   NROOT,NR1,NCG,IG1,IMAXP1,JMAXP1,KMAXP1,NNMAX,ILE,ISTRID,JSTRID 

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      REAL    :: BKOKO1,RC

      REAL, ALLOCATABLE, DIMENSION(:) :: FTR1

      LOGICAL :: MASTER

      MASTER = IPRO == 1

      IF(MASTER) THEN
         WRITE(28) NBLOCG
         WRITE(28) (IMAXG(N)+1,JMAXG(N)+1,KMAXG(N)+1,N=1,NBLOCG)
      ENDIF

      IMEM = 0
      
      IF(MASTER) THEN
         DO NCG = 1,NBLOCG
         IMEM = MAX0(IMEM,(IMAXG(NCG)+1)*(JMAXG(NCG)+1)*(KMAXG(NCG)+1))
         ENDDO
      ELSE
         DO N = 1,NBLOCK
         IMEM = MAX0(IMEM,(IMAX(1,N)+1)*(JMAX(1,N)+1)*(KMAX(1,N)+1))
         ENDDO
      ENDIF

      IMEM = 3*IMEM

      ALLOCATE(FTR1(IMEM),STAT=IERRCODE)
      WRITE(45,*)

      IF (IERRCODE > 0) THEN
            WRITE(*,1022)  IPRO,REAL(4*IMEM,4)/1048576.
            WRITE(45,1022) IPRO,REAL(4*IMEM,4)/1048576.
            IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
 1022 FORMAT(/'NS3D :      Not enough memory in process ',I3,
     +     '. Desired ',F6.2,'MB. for FTR1. Aborting...'/)
      ELSE
            KOKO1  = IMEM
            BKOKO1 = (REAL(KOKO1)*4.)/1048576.
            WRITE(45,9050) KOKO1,BKOKO1
9050  FORMAT('NS3D: FTR1 allocated for post-processing',
     +  /'Size (REAL) = ',I8,' = ',F6.2,' Mbytes')

      ENDIF

      WRITE(45,*)
      WRITE(45,'(A)') 'Grid angles after writing XYD.BIN:'
      WRITE(45,*)
      WRITE(45,9470)
      WRITE(45,9475)
      WRITE(45,9480)

      DO NCG = 1,NBLOCG

      N     = NLOCAL(NCG)       ! local block number
      NROOT = NPNUM(NCG)        ! process number
      NR1   = NROOT - 1         ! process number used by MPI

      NNMAX = (IMAXG(NCG)+1)*(JMAXG(NCG)+1)*(KMAXG(NCG)+1)

      IF(N /= -1) THEN

         IG1     = IG(1,N)
         IMAXP1  = IMAX(1,N) + 1
         JMAXP1  = JMAX(1,N) + 1
         KMAXP1  = KMAX(1,N) + 1
         ISTRID  = IMAX(1,N) + 2*IN
         JSTRID  = JMAX(1,N) + 2*JN
         ILE     = ISTRID*JSTRID

         DO K=1,KMAXP1
            DO J=1,JMAXP1
               DO I=1,IMAXP1
                  II = IG1+I+IN-1 + (J+JN-1)*ISTRID + (K+KN-1)*ILE
                  JJ = I + (J-1)*IMAXP1 + (K-1)*IMAXP1*JMAXP1
                  FTR1(JJ        ) = XCO(II)
                  FTR1(JJ+  NNMAX) = YCO(II)
                  FTR1(JJ+2*NNMAX) = ZCO(II)
               ENDDO
            ENDDO
         ENDDO

         IF(.NOT.MASTER) THEN
            K = NCG
            CALL MPI_SEND(FTR1,3*NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,IERR)
         ENDIF

      ELSEIF(N == -1 .AND. MASTER) THEN

         K = NCG
         CALL MPI_RECV(FTR1,3*NNMAX,MPI_REAL8,NR1,K,
     +        MPI_COMM_WORLD,STATUS,IERR)

      ENDIF 

      IF(MASTER) WRITE(28) (FTR1(I),I=1,3*NNMAX)
              
      WRITE(45,9500) NCG,OMEGA(NCG),ROTANG(NCG)*180./PII,ROTANG(NCG),
     +   ROTAT(NCG)*180./PII,ROTAT(NCG)

c9470  FORMAT('BLOCK','  OMEGA (1/s) ','  Current angle (ROTANG)      ',
c     +       '     Initial angle  (ROTAT)    ')
c9475  FORMAT(22X,'(degrees)',4X,'(radians)',11X,'(degrees)',5X,
c     +           '(radians)')
c9480  FORMAT(80('='))
c9500  FORMAT(I3,1X,F10.2,6X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)

9470  FORMAT(' BLOCK','  OMEGA (1/s) ','  Current angle (ROTANG)      ',
     +       '     Initial angle  (ROTAT)    ')
9475  FORMAT(23X,'(degrees)',4X,'(radians)',11X,'(degrees)',5X,
     +           '(radians)')
9480  FORMAT(80('='))
9500  FORMAT(I6,1X,F10.2,4X,F11.5,1X,F11.5,10X,F11.5,1X,F11.5)

      ENDDO
        
      IF(MASTER) CLOSE(28)

      DEALLOCATE(FTR1,STAT=IERRCODE)

      IF(IERRCODE <= 0) THEN
         WRITE(45,9070) 
      ELSE
         WRITE(45,9080) IERRCODE
      ENDIF

9070  FORMAT('FTR1 deallocated')
9080  FORMAT('Problems with FTR1, continuing until the end. IERR = ',I6)

      RETURN
      END SUBROUTINE WRIXYD
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FVFUN(SIMULATION_END,FIRSTOUTP)

      USE MPI

      USE CHARACTERS

      USE INTEGERS, ONLY : IPRO

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,F1RM,F1RN,F1RW,XCO,
     &    YCO,ZCO,RM,RN,RW,U,V,W,TOLD,UOLD,VOLD,WOLD,HAT3,POLD,
     &    IG,NTOT,S11,PDIFF,P,DE,TEMP,C,VIS,CH,DM,DN,DW,VIST,RK,
     &    REPS,DDEPS,PTUR,FUN1,TTS,RKSI,HAT1,DISTW,HAT2,PRO,VAR,
     &    VTRAN, IMAXG, JMAXG, KMAXG, NPNUM, NLOCAL, STRAIN, OHMI,
     &    RMAV2,RNAV2,RWAV2,RUVAV,RUWAV,RVWAV,RMAV1,RMAV3,TRM,
     &    VELLAP,QSAS,RLDIST,TIJ,RO,A1,A2,A3,VOL,D1,D2,D3,IDER,
     &    A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,XC,YC,ZC,
     &    UROT,VROT,WROT,IT,IL,IK,NPROCE,PAV1,TAV1,ZZZ,BLKS,
     &    RS1,RS2,RS3,RS4,RS5,RS6,CVV1,CVV2,CVV3,CVV4,CVV5,CVV6,
     &    CVV7,CVV8,CVV9,GVX,GVY,GVZ,MPUREL,MPVREL,MPWREL

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: N, IMAXP1, JMAXP1, KMAXP1, IGA, ISTRID, ILE, NISTR, 
     &           NJSTR, NKSTR, NN, N1, NOVAR, NMAX, L, RC, NGL,
     &           IG1, IGN, NTRAIN

      LOGICAL :: MASTER, SLAVE, SIMULATION_END

      REAL, ALLOCATABLE, DIMENSION(:) :: TRAINW

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: NSPOINTS, NCELLS, IPB, LL, MM, NG, IERR  
      INTEGER :: NROOT, NR1, ICHAN, ICHAN1, ICHAN2, ICHAN3
      INTEGER :: IMX, JMX, KMX, I, J, K, IPHASE

      INTEGER(KIND=8) :: II, I1, I2, J1, J2

      INTEGER, DIMENSION(NOVARMAX) :: FI, JUMP, JF, JL

      LOGICAL, DIMENSION(NOVARMAX) :: F

      LOGICAL ::FIRSTOUTP

      CHARACTER(LEN=35), DIMENSION(NOVARMAX) :: VNAME
      CHARACTER(LEN=35) :: VANAME
      CHARACTER(LEN=80) :: FNAME
      CHARACTER(LEN=80) :: ENSNAME
      CHARACTER(LEN=2)  :: FABC
      CHARACTER(LEN=3)  :: FNUM
      CHARACTER(LEN=6)  :: FILEE


      MASTER = IPRO == 1
      SLAVE  = .NOT.MASTER

      
C ... Initialize the output variable selector
         
      FI = OUTVAR
      
      IF(.NOT.MULPHL) THEN
         FI(42:47) = 0            ! Data not available
         IF(.NOT.TWO_FLUIDL) THEN
            FI(48:71) = 0         ! Data not available
         ENDIF
      ENDIF

      IF(.NOT.TRANSL) THEN
         FI(72:73) = 0            ! Data not available
      ENDIF
      
C ... Generate logical function selector

      F = .FALSE.

      DO L = 1,NOVARMAX
         IF(FI(L) == 1)                           F(L) = .TRUE.
         IF(FI(L) == 2 .AND. SIMULATION_END)      F(L) = .TRUE.
         IF(FI(L) == 3 .AND. .NOT.SIMULATION_END) F(L) = .TRUE.
         FI(L) = 0; IF(F(L)) FI(L) = 1
      ENDDO

      NOVAR = SUM(FI)  ! Number of variables to be output for post-processing

      IF(NOVAR == 0) RETURN

      
      IF(TIMEL) THEN
         CALL NUMCH6(FILEE,ITIMES)
      ELSE
         CALL NUMCH6(FILEE,ICYCLE)
      ENDIF


      DO 500 N = 1,NBLOCK

         IMAXP1 = IMAX(1,N) + 1
         JMAXP1 = JMAX(1,N) + 1
         KMAXP1 = KMAX(1,N) + 1
         IGA    = IG(1,N)
         ISTRID = IMAX(1,N) + 2*IN
         ILE    = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)
         NISTR  = 1
         NJSTR  = IMAXP1
         NKSTR  = IMAXP1*JMAXP1
         NMAX   = IMAXP1*JMAXP1*KMAXP1

         NGL    = NPROCE(1+N,IPRO) ! Someday you need this

C ... Calculate covariant vectors to temporary arrays

         CALL COVEC(CVV1(IGA),CVV2(IGA),CVV3(IGA),
     +              XCO(IGA),YCO(IGA),ZCO(IGA),IMAXP1,JMAXP1,KMAXP1,
     +              1,ISTRID,ILE,NISTR,NJSTR,NKSTR)
         CALL COVEC(CVV5(IGA),CVV6(IGA),CVV4(IGA),
     +              YCO(IGA),ZCO(IGA),XCO(IGA),JMAXP1,KMAXP1,IMAXP1,
     +              ISTRID,ILE,1,NJSTR,NKSTR,NISTR)
         CALL COVEC(CVV9(IGA),CVV7(IGA),CVV8(IGA),
     +              ZCO(IGA),XCO(IGA),YCO(IGA),KMAXP1,IMAXP1,JMAXP1,
     +              ILE,1,ISTRID,NKSTR,NISTR,NJSTR)

500   CONTINUE
      

      IF(MASTER .AND. PLOT3D) THEN

         FNAME = 'FUN.BIN'
         IF(.NOT.SIMULATION_END) FNAME = 'MOVIE/FUN'//FILEE//'.BIN'
         FNAME = TRIM(FNAME)
         OPEN(70,FILE=FNAME,FORM='UNFORMATTED')
         FNAME = 'FUN.nam'
         IF(.NOT.SIMULATION_END) FNAME = 'MOVIE/FUN'//FILEE//'.nam'
         FNAME = TRIM(FNAME)
         OPEN(71,FILE=FNAME,FORM='FORMATTED')

      ENDIF


C ... Generate variable name table
      
      IF(MASTER) THEN
         
         VNAME( 1) = 'Pressure_difference'
         VNAME( 2) = 'Temperature'
         VNAME( 3) = 'Pressure_coefficient'
         VNAME( 4) = 'Pressure'
         VNAME( 5) = 'Strain?'
         VNAME( 6) = 'Speed_of_sound'
         VNAME( 7) = 'Viscosity'
         VNAME( 8) = 'Heat_conductivity'
         VNAME( 9) = 'Vorticity_(x-dir);_Vorticity'
         VNAME(10) = 'Vorticity_(y-dir)'
         VNAME(11) = 'Vorticity_(z-dir)'
         VNAME(12) = 'Eddy_viscosity'
         VNAME(13) = 'Kinetic_energy_of_turbulence'
         VNAME(14) = 'Variable_epsilon'
         VNAME(15) = 'True_epsilon'
         VNAME(16) = 'Production_of_turbulence'
         VNAME(17) = "Menter's_blending_function"
         VNAME(18) = 'Turbulence_time-scale'
         VNAME(19) = 'Wall_distance'
         VNAME(20) = 'Viscous_dissipation'
         VNAME(21) = 'ru''u'''
         VNAME(22) = 'ru''v'''
         VNAME(23) = 'ru''w'''
         VNAME(24) = 'rv''v'''
         VNAME(25) = 'rv''w'''
         VNAME(26) = 'rw''w'''
         VNAME(27) = 'ncix;_Cov._vector_in_i-dir.'
         VNAME(28) = 'nciy'
         VNAME(29) = 'nciz'
         VNAME(30) = 'ncjx;_Cov._vector_in_j-dir.'
         VNAME(31) = 'ncjy'
         VNAME(32) = 'ncjz'
         VNAME(33) = 'nckx;_Cov._vector_in_k-dir.'
         VNAME(34) = 'ncky'
         VNAME(35) = 'nckz'
         VNAME(36) = 'Rksi'
         VNAME(37) = 'VTRAN'
         VNAME(38) = "Lighthill's_tensor_gradient"
         VNAME(39) = 'gvex;_Grid_motion_velocity_vector'
         VNAME(40) = 'gvey'
         VNAME(41) = 'gvez'
         VNAME(42) = 'Evaporation_rate'
         VNAME(43) = 'Gas_volume_fraction'
         VNAME(44) = 'Liquid_temperature'
         VNAME(45) = 'Gas_temperature'
         VNAME(46) = 'Saturation_temperature'
         VNAME(47) = 'Saturation_pressure'
         VNAME(48) = 'UL;_Liquid_velocity'
         VNAME(49) = 'VL'
         VNAME(50) = 'WL'
         VNAME(51) = 'UG;_Gas_velocity'
         VNAME(52) = 'VG'
         VNAME(53) = 'WG'
         VNAME(54) = 'FITX;_Interfacial_force'
         VNAME(55) = 'FITY'
         VNAME(56) = 'FITZ'
         VNAME(57) = 'FDX;_Interfacial_drag_force'
         VNAME(58) = 'FDY'
         VNAME(59) = 'FDZ'
         VNAME(60) = 'FLX;_Interfacial_lift_force'
         VNAME(61) = 'FLY'
         VNAME(62) = 'FLZ'
         VNAME(63) = 'FTDX;_Turbulent_dispersion'
         VNAME(64) = 'FTDY'
         VNAME(65) = 'FTDZ'
         VNAME(66) = 'FVMX;_Virtual_mass'
         VNAME(67) = 'FVMY'
         VNAME(68) = 'FVMZ'
         VNAME(69) = 'FWX;_Wall_force'
         VNAME(70) = 'FWY'
         VNAME(71) = 'FWZ'
         VNAME(72) = 'Intermittency'
         VNAME(73) = 'Re_mt'
         VNAME(74) = 'Q-criterion'

      ENDIF
      

C ... Write selected function names to the FUN.nam file. Some post-processing
C ... programs don't like underscores in function names, but EnSight requires
C ... them.
      
      IF(MASTER .AND. PLOT3D) THEN

         DO L = 1,NOVARMAX
            IF(F(L)) THEN
               VANAME = VNAME(L)
               DO N=1,35
                  IF(VANAME(N:N) == '_') VANAME(N:N) = ' '
               ENDDO
               WRITE(71,"(A)") TRIM(VANAME)
            ENDIF
         ENDDO

      ENDIF


      IF(MASTER .AND. PLOT3D) THEN
         WRITE(70) NBLOCG
         WRITE(70)(IMAXG(N)+1,JMAXG(N)+1,KMAXG(N)+1,NOVAR, N=1,NBLOCG)
      ENDIF

      
      IF(SIMULATION_END) THEN  ! Write only at the end of the simulation

         IF(PARALLEL .AND. MASTER) THEN
            WRITE(*,*)
            WRITE(*,*)' Collecting and writing FUN.BIN  ...'
         ENDIF

      ENDIF  ! Write only at the end of the simulation


      DO NG = 1,NBLOCG
         
         N     = NLOCAL(NG)     ! local block number
         NROOT = NPNUM(NG)      ! process number
         NR1   = NROOT - 1      ! process number used by MPI

         NSPOINTS = (IMAXG(NG)+1)*(JMAXG(NG)+1)*(KMAXG(NG)+1)
         NCELLS   = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
         NTRAIN   = NOVAR*NCELLS

C ... Initialize the TRAINW pointer jump array (cells)

         DO L=1,NOVARMAX
            JUMP(L) = (SUM(FI(1:L))-1)*NCELLS
            JF(L)   = JUMP(L)+1
            JL(L)   = JUMP(L)+NCELLS
         ENDDO

         IF(MASTER) ALLOCATE(TRAINW(NTRAIN)) 


         IF(N /= -1) THEN

         IF(SLAVE) ALLOCATE(TRAINW(NTRAIN)) 

         IG1 = IG(1,N)
         IGN = IG1 + NCELLS - 1

         
         IF(MULPHL) THEN        ! Multiphase variables

            DO I = IG1,IGN
               VAR(I)%F1R(1) = PRO(I)%TSAT
               VAR(I)%F1R(2) = MIN(VAR(I)%BUBDIA(1),VAR(I)%BUBDIA(2))
C ...          PRO(I)%DPSAT  ! Varo muutosta
            ENDDO
                 
            IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN

C ... Total interfacial force (or something else)

            DO I = IG1,IGN
               PRO(I)%FITOT(1) = PRO(I)%EO
               PRO(I)%FITOT(2) = MIN(VAR(I)%BUBDIA(1),VAR(I)%BUBDIA(2))
               PRO(I)%FITOT(3) = VAR(I)%BUBDIA(2)
            ENDDO

            ENDIF               ! MULTI
         
            DO IPHASE = 1,BLKS(NGL)%NPHASE

               DO I = IG1,IGN
                  PRO(I)%TEMP(IPHASE) = PRO(I)%DTEMP(IPHASE)
               ENDDO

            ENDDO

         ENDIF                  ! MULPHL


         IF(F( 1)) TRAINW(JF( 1):JL( 1)) = PDIFF(IG1:IGN)
         IF(F( 2)) TRAINW(JF( 2):JL( 2)) = TEMP(IG1:IGN)
         IF(F( 3)) TRAINW(JF( 3):JL( 3)) = PDIFF(IG1:IGN)/REFPRE
         IF(F( 4)) TRAINW(JF( 4):JL( 4)) = P(IG1:IGN)
         IF(F( 5)) TRAINW(JF( 5):JL( 5)) = STRAIN(IG1:IGN)
         IF(F( 6)) TRAINW(JF( 6):JL( 6)) = C(IG1:IGN)
         IF(F( 7)) TRAINW(JF( 7):JL( 7)) = VIS(IG1:IGN)
         IF(F( 8)) TRAINW(JF( 8):JL( 8)) = CH(IG1:IGN)
         IF(F( 9)) TRAINW(JF( 9):JL( 9)) = DM(IG1:IGN)
         IF(F(10)) TRAINW(JF(10):JL(10)) = DN(IG1:IGN)
         IF(F(11)) TRAINW(JF(11):JL(11)) = DW(IG1:IGN)
         IF(F(12)) TRAINW(JF(12):JL(12)) = VIST(IG1:IGN)
         IF(F(13)) TRAINW(JF(13):JL(13)) = RK(IG1:IGN)
         IF(F(14)) TRAINW(JF(14):JL(14)) = REPS(IG1:IGN)
         IF(F(15)) TRAINW(JF(15):JL(15)) = DDEPS(IG1:IGN)
         IF(F(16)) TRAINW(JF(16):JL(16)) = PTUR(IG1:IGN)
         IF(F(17)) TRAINW(JF(17):JL(17)) = FUN1(IG1:IGN)
         IF(F(18)) TRAINW(JF(18):JL(18)) = TTS(IG1:IGN)
         IF(F(19)) TRAINW(JF(19):JL(19)) = DISTW(IG1:IGN)
         IF(F(20)) TRAINW(JF(20):JL(20)) = HAT2(IG1:IGN)
         IF(F(21)) TRAINW(JF(21):JL(21)) = RS1(IG1:IGN)
         IF(F(22)) TRAINW(JF(22):JL(22)) = RS2(IG1:IGN)
         IF(F(23)) TRAINW(JF(23):JL(23)) = RS3(IG1:IGN)
         IF(F(24)) TRAINW(JF(24):JL(24)) = RS4(IG1:IGN)
         IF(F(25)) TRAINW(JF(25):JL(25)) = RS5(IG1:IGN)
         IF(F(26)) TRAINW(JF(26):JL(26)) = RS6(IG1:IGN)
         IF(F(27)) TRAINW(JF(27):JL(27)) = CVV1(IG1:IGN)
         IF(F(28)) TRAINW(JF(28):JL(28)) = CVV2(IG1:IGN)
         IF(F(29)) TRAINW(JF(29):JL(29)) = CVV3(IG1:IGN)
         IF(F(30)) TRAINW(JF(30):JL(30)) = CVV4(IG1:IGN)
         IF(F(31)) TRAINW(JF(31):JL(31)) = CVV5(IG1:IGN)
         IF(F(32)) TRAINW(JF(32):JL(32)) = CVV6(IG1:IGN)
         IF(F(33)) TRAINW(JF(33):JL(33)) = CVV7(IG1:IGN)
         IF(F(34)) TRAINW(JF(34):JL(34)) = CVV8(IG1:IGN)
         IF(F(35)) TRAINW(JF(35):JL(35)) = CVV9(IG1:IGN)
         IF(F(36)) TRAINW(JF(36):JL(36)) = RKSI(IG1:IGN)
         IF(F(37)) TRAINW(JF(37):JL(37)) = VTRAN(IG1:IGN)
         IF(F(38)) TRAINW(JF(38):JL(38)) = TIJ(IG1:IGN)
         IF(F(39)) TRAINW(JF(39):JL(39)) = GVX(IG1:IGN)
         IF(F(40)) TRAINW(JF(40):JL(40)) = GVY(IG1:IGN)
         IF(F(41)) TRAINW(JF(41):JL(41)) = GVZ(IG1:IGN)
         IF(F(42)) TRAINW(JF(42):JL(42)) = VAR(IG1:IGN)%EVAPR(2)
         IF(F(43)) TRAINW(JF(43):JL(43)) = VAR(IG1:IGN)%ALFA(2)
         IF(F(44)) TRAINW(JF(44):JL(44)) = PRO(IG1:IGN)%TEMP(1) ! Mersu
         IF(F(45)) TRAINW(JF(45):JL(45)) = PRO(IG1:IGN)%TEMP(2)
         IF(F(46)) TRAINW(JF(46):JL(46)) = VAR(IG1:IGN)%F1R(1)
         IF(F(47)) TRAINW(JF(47):JL(47)) = VAR(IG1:IGN)%F1R(2)  ! Kamalaa
         IF(F(48)) TRAINW(JF(48):JL(48)) = MPUREL(1,IG1:IGN)
         IF(F(49)) TRAINW(JF(49):JL(49)) = MPVREL(1,IG1:IGN)
         IF(F(50)) TRAINW(JF(50):JL(50)) = MPWREL(1,IG1:IGN)
         IF(F(51)) TRAINW(JF(51):JL(51)) = MPUREL(2,IG1:IGN)
         IF(F(52)) TRAINW(JF(52):JL(52)) = MPVREL(2,IG1:IGN)
         IF(F(53)) TRAINW(JF(53):JL(53)) = MPWREL(2,IG1:IGN)
         IF(F(54)) TRAINW(JF(54):JL(54)) = PRO(IG1:IGN)%FITOT(1)
         IF(F(55)) TRAINW(JF(55):JL(55)) = PRO(IG1:IGN)%FITOT(2)
         IF(F(56)) TRAINW(JF(56):JL(56)) = PRO(IG1:IGN)%FITOT(3)
         IF(F(57)) TRAINW(JF(57):JL(57)) = PRO(IG1:IGN)%FD(1)/REFPRE
         IF(F(58)) TRAINW(JF(58):JL(58)) = PRO(IG1:IGN)%FD(2)/REFPRE
         IF(F(59)) TRAINW(JF(59):JL(59)) = PRO(IG1:IGN)%FD(3)/REFPRE
         IF(F(60)) TRAINW(JF(60):JL(60)) = PRO(IG1:IGN)%FL(1)/REFPRE
         IF(F(61)) TRAINW(JF(61):JL(61)) = PRO(IG1:IGN)%FL(2)/REFPRE
         IF(F(62)) TRAINW(JF(62):JL(62)) = PRO(IG1:IGN)%FL(3)/REFPRE
         IF(F(63)) TRAINW(JF(63):JL(63)) = PRO(IG1:IGN)%FTD(1)/REFPRE
         IF(F(64)) TRAINW(JF(64):JL(64)) = PRO(IG1:IGN)%FTD(2)/REFPRE
         IF(F(65)) TRAINW(JF(65):JL(65)) = PRO(IG1:IGN)%FTD(3)/REFPRE
         IF(F(66)) TRAINW(JF(66):JL(66)) = PRO(IG1:IGN)%FVM(1)/REFPRE
         IF(F(67)) TRAINW(JF(67):JL(67)) = PRO(IG1:IGN)%FVM(2)/REFPRE
         IF(F(68)) TRAINW(JF(68):JL(68)) = PRO(IG1:IGN)%FVM(3)/REFPRE
         IF(F(69)) TRAINW(JF(69):JL(69)) = PRO(IG1:IGN)%FW(1)/REFPRE
         IF(F(70)) TRAINW(JF(70):JL(70)) = PRO(IG1:IGN)%FW(2)/REFPRE
         IF(F(71)) TRAINW(JF(71):JL(71)) = PRO(IG1:IGN)%FW(3)/REFPRE
         IF(F(72)) TRAINW(JF(72):JL(72)) = TRM(IG1:IGN)%G
         IF(F(73)) TRAINW(JF(73):JL(73)) = TRM(IG1:IGN)%RET
         IF(F(74)) TRAINW(JF(74):JL(74)) = 0.25*OHMI(IG1:IGN)**2
     &                                   - 0.50*STRAIN(IG1:IGN)**2

         IMX = IMAX(1,N)
         JMX = JMAX(1,N)
         KMX = KMAX(1,N)

C ... From cell center points to node points
         
         DO LL = 0,NOVAR-1
            CALL REDUCE(TRAINW(LL*NCELLS+1),ZZZ,IMX,JMX,KMX)
         ENDDO

C ... The contravariant vectors were calculated in node points already
C ... before the REDUCE call so they must be reset.         
         
         IF(F(27)) TRAINW(JF(27):JL(27)) = CVV1(IG1:IGN)
         IF(F(28)) TRAINW(JF(28):JL(28)) = CVV2(IG1:IGN)
         IF(F(29)) TRAINW(JF(29):JL(29)) = CVV3(IG1:IGN)
         IF(F(30)) TRAINW(JF(30):JL(30)) = CVV4(IG1:IGN)
         IF(F(31)) TRAINW(JF(31):JL(31)) = CVV5(IG1:IGN)
         IF(F(32)) TRAINW(JF(32):JL(32)) = CVV6(IG1:IGN)
         IF(F(33)) TRAINW(JF(33):JL(33)) = CVV7(IG1:IGN)
         IF(F(34)) TRAINW(JF(34):JL(34)) = CVV8(IG1:IGN)
         IF(F(35)) TRAINW(JF(35):JL(35)) = CVV9(IG1:IGN)

C ... The grid velocities were calculated in node points already
C ... before the REDUCE call so they must be reset.         
         
         IF(F(39)) TRAINW(JF(39):JL(39)) = GVX(IG1:IGN)
         IF(F(40)) TRAINW(JF(40):JL(40)) = GVY(IG1:IGN)
         IF(F(41)) TRAINW(JF(41):JL(41)) = GVZ(IG1:IGN)

C ... Pack the TRAINW vector.
         
         DO LL = 0,NOVAR-1
            I1 = LL*NSPOINTS+1
            I2 = I1+NSPOINTS-1
            J1 = LL*NCELLS
            DO II = I1,I2
               J1 = J1+1
               TRAINW(II) = TRAINW(J1)
            ENDDO
          ENDDO

         
C ... Initialize the TRAINW pointer jump array (grid points)

         DO L=1,NOVARMAX
            JUMP(L) = (SUM(FI(1:L))-1)*NSPOINTS
         ENDDO

         IF(SLAVE) THEN
            CALL MPI_SEND(TRAINW,NOVAR*NSPOINTS,MPI_REAL8,0,
     &           NG,MPI_COMM_WORLD,IERR)
            DEALLOCATE(TRAINW)
         ENDIF

         ELSEIF(MASTER .AND. N == -1) THEN

            CALL MPI_RECV(TRAINW,NOVAR*NSPOINTS,MPI_REAL8,NR1,NG,
     &           MPI_COMM_WORLD,STATUS,IERR)

         ENDIF                     ! (N /= -1)

            IF(MASTER) THEN

               IF(PLOT3D) WRITE(70)
     &                   (REAL(TRAINW(II),4), II=1,NOVAR*NSPOINTS)
 
               IF(ENSIGHT) THEN

                  DO L = 1,NOVARMAX

                     IF(F(L)) THEN

                        ICHAN = 1138+L
                        CALL NUM2ABC(FABC,L)

                        IF(SIMULATION_END) THEN
                           FNAME   = 'FUN.BIN.'//FABC
                           ENSNAME = 'scalar per node: '//
     &                     VNAME(L)//'FUN.BIN.'//FABC
                        ELSE
                           FNAME   = 'MOVIE/FUN'//FILEE//'.BIN.'//FABC
                           ENSNAME = 'scalar per node: '//
     &                     VNAME(L)//'FUN******.BIN.'//FABC
                        ENDIF

                        FNAME = TRIM(FNAME)
                        
                        CALL WriteEnsightScalar(FNAME,ICHAN,NG,NBLOCG,
     &                       IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                       TRAINW(JUMP(L)+1))

                        IF(NG == 1 .AND. (SIMULATION_END.OR.FIRSTOUTP))
     &                     WRITE(1130,'(A)') ENSNAME
       
                     ENDIF

                  ENDDO

* If you need Reynolds stress tensor in EnSight Gold format remove '*'s.
*
*                  IF(SUM(FI(21:26)) == 6) THEN ! Reynolds stresses
*
*                     ICHAN = ICHAN + 1
*
*                        IF(SIMULATION_END) THEN
*                           FNAME = 'FUN.BIN.Rey'
*                           ENSNAME = 'tensor symm per node: '//
*     &                     'Reynolds_stresses             '//FNAME
*                        ELSE
*                           FNAME = 'MOVIE/FUN'//FILEE//'.BIN.Rey'
*                           ENSNAME = 'tensor symm per node: '//
*     &                               'Reynolds_stresses             '//
*     &                               'FUN******.BIN.Rey'
*                        ENDIF
*
*                        FNAME = TRIM(FNAME)
* 
*                        CALL WriteEnsightTensor(FNAME,ICHAN,NG,NBLOCG,
*     &                       IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
*     &                       TRAINW(JUMP(21)+1),TRAINW(JUMP(24)+1),
*     &                       TRAINW(JUMP(26)+1),TRAINW(JUMP(22)+1),
*     &                       TRAINW(JUMP(23)+1),TRAINW(JUMP(25)+1))
*
*                        IF(NG == 1 .AND. (SIMULATION_END.OR.FIRSTOUTP))
*     &                     WRITE(1130,'(A)') ENSNAME
*                    
*                  ENDIF

               ENDIF

               DEALLOCATE(TRAINW)

            ENDIF


      ENDDO ! Global block loop

      IF(MASTER .AND. PLOT3D) CLOSE(70)
      IF(MASTER .AND. PLOT3D) CLOSE(71)

C **********************************************************************

C ... Time averaged values

C **********************************************************************

      IF (TIMEL .AND. SIMULATION_END) THEN

      IF(MASTER.AND.PLOT3D) OPEN(63,FILE='FUN.AVE',FORM='UNFORMATTED')
      IF(MASTER.AND.PLOT3D) OPEN(64,FILE='FUN_AVE.nam',FORM='FORMATTED')

      NOVAR = 8

      IF(.NOT.PARALLEL .AND. PLOT3D) THEN
         WRITE(63) NBLOCK
         WRITE(63) (IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,NOVAR,
     +              N = 1,NBLOCK)
      ELSE
         IF(MASTER .AND. PLOT3D) THEN
            WRITE(63) NBLOCG
            WRITE(63)(IMAXG(N)+1,JMAXG(N)+1,KMAXG(N)+1,NOVAR,N=1,NBLOCG)
         ENDIF
      ENDIF


      IF(MASTER .AND. PLOT3D) THEN
         WRITE(64,"(A)") 'Pressure'
         WRITE(64,"(A)") 'Temperature'
         WRITE(64,"(A)") 'ru''u'''
         WRITE(64,"(A)") 'ru''v'''
         WRITE(64,"(A)") 'ru''w'''
         WRITE(64,"(A)") 'rv''v'''
         WRITE(64,"(A)") 'rv''w'''
         WRITE(64,"(A)") 'rw''w'''
      ENDIF

      IF(.NOT.PARALLEL) THEN

         DO N = 1,NBLOCK

            NMAX = (IMAX(1,N) + 1)*(JMAX(1,N) + 1)*(KMAX(1,N) + 1)
            IG1  = IG(1,N) - 1

            IF(PLOT3D)
     +      WRITE(63) (REAL( PAV1(IG1+L),4),L=1,NMAX),
     +                (REAL( TAV1(IG1+L),4),L=1,NMAX),
     +                (REAL(RMAV2(IG1+L),4),L=1,NMAX),
     +                (REAL(RUVAV(IG1+L),4),L=1,NMAX),
     +                (REAL(RUWAV(IG1+L),4),L=1,NMAX),
     +                (REAL(RNAV2(IG1+L),4),L=1,NMAX),
     +                (REAL(RVWAV(IG1+L),4),L=1,NMAX),
     +                (REAL(RWAV2(IG1+L),4),L=1,NMAX)

            IF(ENSIGHT) THEN
              
               IF(N == 1) THEN
                  ICHAN1 = ICHAN + 1
                  ICHAN2 = ICHAN + 2
                  ICHAN3 = ICHAN + 3
               ENDIF
                  
               VANAME = 'Pressure_average'
               CALL WriteEnsightScalar('FUN.AVE.p',ICHAN1,
     &              N,NBLOCK,IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,
     &              PAV1(IG1+1))
               IF(N == 1) WRITE(1130,'(A)') 
     &            'scalar per node: '//VANAME//'FUN.AVE.p' 

               VANAME = 'Temperature_average'
               CALL WriteEnsightScalar('FUN.AVE.T',ICHAN2,
     &              N,NBLOCK,IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,
     &              TAV1(IG1+1))
               IF(N == 1) WRITE(1130,'(A)') 
     &            'scalar per node: '//VANAME//'FUN.AVE.T' 


* If you need the Reynolds stress average tensor in EnSight Gold format at the
* end of the simulation remove '*'s.               
*
*               FNAME = TRIM('FUN.AVE.Rey')
*               ENSNAME = 'tensor symm per node: '//
*     &                   'Reynolds stress averages      '//FNAME
*               CALL WriteEnsightTensor(FNAME,ICHAN3,N,NBLOCK,
*     &              IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,
*     &              RMAV2(IG1+1),RNAV2(IG1+1),RWAV2(IG1+1),
*     &              RUVAV(IG1+1),RUWAV(IG1+1),RVWAV(IG1+1))
*
*               IF(N == 1) WRITE(1130,'(A)') ENSNAME
               
            ENDIF

         ENDDO

      ELSE  ! Parallel

         IF(MASTER) THEN
            WRITE(*,*)
            WRITE(*,*)' Collecting and writing FUN.AVE  ...'
         ENDIF

         DO NG  = 1,NBLOCG

            N     = NLOCAL(NG)  ! local block number
            NROOT = NPNUM(NG)   ! process number
            NR1   = NROOT - 1   ! process number used by MPI

            NSPOINTS = (IMAXG(NG)+1)*(JMAXG(NG)+1)*(KMAXG(NG)+1)
            NTRAIN   = 8*NSPOINTS

            IF(MASTER) ALLOCATE(TRAINW(NTRAIN)) 


            IF(N /= -1) THEN

               IF(.NOT.MASTER) ALLOCATE(TRAINW(NTRAIN)) 

               IG1    = IG(1,N)
               IMAXP1 = IMAX(1,N) + 1
               JMAXP1 = JMAX(1,N) + 1
               KMAXP1 = KMAX(1,N) + 1
               IPB    = 0

               DO K = 1,KMAXP1
                  DO J = 1,JMAXP1
                     DO I = 1,IMAXP1

                        IPB = IPB + 1
                        LL  = IG1-1+I+IMAXP1*(J-1)+IMAXP1*JMAXP1*(K-1)

                        TRAINW(IPB             ) =  PAV1(LL)
                        TRAINW(IPB + 1*NSPOINTS) =  TAV1(LL)
                        TRAINW(IPB + 2*NSPOINTS) = RMAV2(LL)
                        TRAINW(IPB + 3*NSPOINTS) = RUVAV(LL)
                        TRAINW(IPB + 4*NSPOINTS) = RUWAV(LL)
                        TRAINW(IPB + 5*NSPOINTS) = RNAV2(LL)
                        TRAINW(IPB + 6*NSPOINTS) = RVWAV(LL)
                        TRAINW(IPB + 7*NSPOINTS) = RWAV2(LL)
      
                     ENDDO            
                  ENDDO            
               ENDDO            

               IF(.NOT.MASTER) THEN
                  CALL MPI_SEND(TRAINW,NTRAIN,MPI_REAL8,0,
     &                          NG,MPI_COMM_WORLD,IERR)
                  DEALLOCATE(TRAINW)
               ENDIF

            ELSEIF(MASTER .AND. N == -1) THEN

               CALL MPI_RECV(TRAINW,NTRAIN,MPI_REAL8,NR1,NG,
     &         MPI_COMM_WORLD,STATUS,IERR)

            ENDIF ! (N /= -1)


            IF(MASTER) THEN
               
               IF(PLOT3D) WRITE(63) (REAL(TRAINW(II),4), II=1,NTRAIN)

               IF(ENSIGHT) THEN
            
                  IF(NG == 1) THEN
                     ICHAN1 = ICHAN + 1
                     ICHAN2 = ICHAN + 2
                     ICHAN3 = ICHAN + 3
                  ENDIF
                  
                  VANAME = 'Pressure_average'
                  CALL WriteEnsightScalar('FUN.AVE.p',ICHAN1,
     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                 TRAINW(1))
                  IF(NG == 1) WRITE(1130,'(A)') 
     &                 'scalar per node: '//VANAME//'FUN.AVE.p' 
 
                  VANAME = 'Temperature_average'
                  CALL WriteEnsightScalar('FUN.AVE.T',ICHAN2,
     &                 NG,NBLOCG,IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
     &                 TRAINW(1+NSPOINTS))
                  IF(NG == 1) WRITE(1130,'(A)') 
     &                 'scalar per node: '//VANAME//'FUN.AVE.T' 

* If you need the Reynolds stress average tensor in EnSight Gold format at the
* end of the simulation remove '*'s.               
*
*                  FNAME = TRIM('FUN.AVE.Rey')
*                  ENSNAME = 'tensor symm per node: '//
*     &                 'Reynolds stress averages      '//FNAME
*                  CALL WriteEnsightTensor(FNAME,ICHAN3,NG,NBLOCG,
*     &                 IMAXG(NG)+1,JMAXG(NG)+1,KMAXG(NG)+1,
*     &                 TRAINW(1+2*NSPOINTS),TRAINW(1+5*NSPOINTS),
*     &                 TRAINW(1+7*NSPOINTS),TRAINW(1+3*NSPOINTS),
*     &                 TRAINW(1+4*NSPOINTS),TRAINW(1+6*NSPOINTS))
*
*                  IF(NG == 1) WRITE(1130,'(A)') ENSNAME
                   
               ENDIF
             
               DEALLOCATE(TRAINW)

            ENDIF 

         ENDDO ! Global block loop

      ENDIF
      
      IF(MASTER .AND. PLOT3D) CLOSE(63)
      IF(MASTER .AND. PLOT3D) CLOSE(64)

      ENDIF  ! TIMEL .AND. SMULATION_END


      RETURN
      END SUBROUTINE FVFUN
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FVPRO

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : IPRO

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,NSPB,NSPP,ICON,
     &    IG,IHF,KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,F1RK,ISTRS,JSTRS,
     &    KSTRS,OHMI,XCP,YCP,ZCP,XC,YC,ZC,D1,D2,D3,RO,VIS,BLKS,
     &    NPROCE,UTAUM 

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: NSP,N,NGL,M,ISP,IWALL,IG1,IF1,ISTR,JSTR,KSTR,
     &           KX1,KX2,KY1,KY2,KZ1,KZ2,IDIR

C ... Calculates properties to be stored in patch function files

      DO 1000 NSP = 1,NSPTOT

         N   = NSPB(NSP)        ! Local block number
         NGL = NPROCE(1+N,IPRO) ! Global block number

         IF(BLKS(NGL)%SOLUTION_TYPE /= 'SOLID') THEN
        
            ISP   = NSPP(NSP)   ! Patch number
            IWALL = ICON((ISP-1)*IC9+ 3) ! WALL NUMBER
            M     = 1           ! Level

            IG1 = IG(1,N)
            IF1 = IHF(ISP,1)

            ISTR = ISTRS(NSP)
            JSTR = JSTRS(NSP)
            KSTR = KSTRS(NSP)

            KX1 = KX1S(NSP)
            KX2 = KX2S(NSP)       
            KY1 = KY1S(NSP)
            KY2 = KY2S(NSP)       

            IF(IWALL == 2 .OR. IWALL == 5) THEN
               KY1 = KX1S(NSP)
               KY2 = KX2S(NSP)
               KX1 = KY1S(NSP)
               KX2 = KY2S(NSP)
            ENDIF
            
            KZ1 = KZ1S(NSP)
            KZ2 = KZ2S(NSP)

            IF(IWALL >= 4) THEN       
               KZ1 = KZ1S(NSP) - 1
               KZ2 = KZ2S(NSP) - 1
            ENDIF

C ... Calculate the y+ values into F1RK array      

            IDIR = -1

            IF(IWALL == 1 .OR. IWALL == 4) THEN
               CALL YPLUS1(F1RK(IF1),OHMI(IG1),
     &              D1(IG1),RO(IG1),VIS(IG1),ISTR,JSTR,KSTR,JN,KN,IN,
     &              KX1,KX2,KY1,KY2,KZ1,KZ2,ISP,IDIR,UTAUM(IF1))
            ELSE IF(IWALL == 2 .OR. IWALL == 5) THEN
               CALL YPLUS1(F1RK(IF1),OHMI(IG1),
     &              D2(IG1),RO(IG1),VIS(IG1),ISTR,JSTR,KSTR,KN,IN,JN,
     &              KX1,KX2,KY1,KY2,KZ1,KZ2,ISP,IDIR,UTAUM(IF1))
            ELSE IF(IWALL == 3 .OR. IWALL == 6) THEN
               CALL YPLUS1(F1RK(IF1),OHMI(IG1),
     &              D3(IG1),RO(IG1),VIS(IG1),ISTR,JSTR,KSTR,IN,JN,KN,
     &              KX1,KX2,KY1,KY2,KZ1,KZ2,ISP,IDIR,UTAUM(IF1))
            ENDIF

            CALL EXTPAT(F1RK(IF1),KX1,KX2,KY1,KY2)
         ENDIF                  ! Not a solid block

 1000 CONTINUE

      RETURN
      END SUBROUTINE FVPRO
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FVSRF(INRE,JNRE,KNRE,SIMULATION_END)

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : IPRO, NBCS, IREPEA, NREPEA

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,IDIMS,JDIMS,NSPB,NSPP,
     &    ICON,IG,IHF,KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,CPWALL,QWALL,
     &    QWFRIC,TWALL,TAUW1,TAUW2,TAUWX,TAUWY,TAUWZ,SURFX,SURFY,
     &    SURFZ,HFLUX,SURLE,DSURLE,ZZZ,F1RK,ISTRS,JSTRS,KSTRS,A1,
     &    A2,A3,WMFLUX,BOUNR,BOUNU,BOUNV,BOUNW,BOUNT,BOUNP,BOUNE,
     &    BOUNRK,BOUNEP,BOUNMF,NSPG,XCO,YCO,ZCO,WTRAN,SURFPX,SURFPY,
     &    SURFPZ,
     &    CPWALL_,QWALL_, QWFRIC_,TWALL_,  TAUW1_,  TAUW2_, 
     &    TAUWX_, TAUWY_, TAUWZ_, SURFX_,  SURFY_,  SURFZ_, 
     &    HFLUX_, SURLE_, DSURLE_,F1RK_,   WTRAN_,  WMFLUX_,
     &    BOUNR_, BOUNT_, BOUNP_, BOUNRK_, BOUNEP_, BOUNMF_,
     &    BOUNU_, BOUNV_, BOUNW_, BOUNE_

      USE NS3CO, ONLY : IC9, LN, REFPRE, PARALLEL, AREF, CHLREF,
     &     GROUP, NSPTOT, ICYCLE, IPRINT, ITIMES, TIMEL, MPRINT

      USE FLIGHT, ONLY : XCGI,YCGI,ZCGI,XCG,YCG,ZCG,PSIR,THETAR,PHIR,
     &     FLYOBJ

      USE FORCE_GROUPS

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SOLIDL

      INTEGER :: NOVAR,KMAXXP,NWOOD,NSP,N,M,ISP,IG1,IF1,KX1,KX2,KY1,
     &    KY2,KZ1,KZ2,I,J,K,KX1W,KX2W,KY1W,KY2W,KZ1W,KZ2W,L,LA,
     &    IPA,IPB,IPC,IPT,IN,JN,KN,INRE,JNRE,KNRE,IC77,
     &    ISTR,JSTR,KSTR,IMAXP1,JMAXP1,IWALL,IPSTR,JPSTR,NN,IERR,II

      INTEGER :: NSG, IWGRID, IGGRID, IWSOLU, IGSOLU, RC

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      LOGICAL :: MASTER, FOUNDL, FOUNDG, FOUND, SIMULATION_END
      LOGICAL :: BEEN_HERE_ALREADY

      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMSW, JDIMSW, ITYPEW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMSG, JDIMSG
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMST, JDIMST
      INTEGER :: IPG, NBCSG, NSPTOTG, NSPTOTT, NSPI, IG1G, ITYPE 
      INTEGER :: NSPOINTS, NSPOINTC, NSPOINTSG
      INTEGER :: NTRAIN, ISLAVE, IDIMT, JDIMT, ITYPET, OWNER, ITSME

      REAL, ALLOCATABLE, DIMENSION(:) :: XCOW, YCOW, ZCOW
      REAL, ALLOCATABLE, DIMENSION(:) :: XCOP, YCOP, ZCOP
      REAL, ALLOCATABLE, DIMENSION(:) :: XCOT, YCOT, ZCOT
      REAL, ALLOCATABLE, DIMENSION(:) :: XPATRIA,YPATRIA 
      REAL, ALLOCATABLE, DIMENSION(:) :: ZPATRIA 

      REAL, ALLOCATABLE, DIMENSION(:) :: SURFXW, SURFYW, SURFZW
      REAL, ALLOCATABLE, DIMENSION(:) :: FRICXW, FRICYW, FRICZW
      REAL, ALLOCATABLE, DIMENSION(:) :: TWALLW, AREAW
      REAL, ALLOCATABLE, DIMENSION(:) :: SURFXP, SURFYP, SURFZP
      REAL, ALLOCATABLE, DIMENSION(:) :: FRICXP, FRICYP, FRICZP
      REAL, ALLOCATABLE, DIMENSION(:) :: FRICXT, FRICYT, FRICZT
      REAL, ALLOCATABLE, DIMENSION(:) :: TWALLP, AREAP     

      REAL, ALLOCATABLE, DIMENSION(:) :: FXWPATRIA, FYWPATRIA, FZWPATRIA      
      
      REAL, ALLOCATABLE, DIMENSION(:) :: SURFXT, SURFYT, SURFZT
      REAL, ALLOCATABLE, DIMENSION(:) :: TWALLT, AREAT

      REAL, ALLOCATABLE, DIMENSION(:) :: TRAINFV, TRAINLI
      REAL, ALLOCATABLE, DIMENSION(:) :: TRAINFW, TRAINLW

      REAL :: SURFXW_D, SURFYW_D, SURFZW_D
      REAL :: SURFXW_APU, SURFYW_APU, SURFZW_APU
      
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: LOADGG
      CHARACTER(LEN=80) :: LOADGT
      CHARACTER(LEN=1)  :: GNAME
      CHARACTER(LEN=6)  :: FILEE


      MASTER = IPRO == 1
      
      IF(TIMEL) THEN
         IC77 = ITIMES
         CALL NUMCH6(FILEE,ITIMES)
      ELSE
         IC77 = ICYCLE
         CALL NUMCH6(FILEE,ICYCLE)
      ENDIF
      
C ... Firstly write the XYZBND.FMT Woodpecker file (without IBLANK)

      IF(PARALLEL) THEN  ! XYZBND.FMT data in a parallel run

C ... Total number of boundary condition patches from all processes (NBCSG)
         CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)
         
C ... Total number of selected surface patches from all processes (NSPTOTG)
         CALL MPI_ALLREDUCE(NSPTOT,NSPTOTG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)

         NSPI = 0

         ALLOCATE(IDIMSW(NSPTOTG),JDIMSW(NSPTOTG),ITYPEW(NSPTOTG))
         ALLOCATE(IDIMSG(NSPTOTG),JDIMSG(NSPTOTG))
         ALLOCATE(IDIMST(NSPTOTG),JDIMST(NSPTOTG))
         ALLOCATE(LOADGG(NSPTOTG))

         LOADGG = ' '


C ... Send the patch dimensions to all processes

         DO IPG = 1,NBCSG  ! Global patch loop

            IF(NSPI == NSPTOTG) CYCLE  ! All patches found

            FOUNDL = .FALSE.
            IDIMT  = 0
            JDIMT  = 0
            ITYPET = 0
            LOADGT = ' '
            ITSME  = 0

            DO NSP = 1,NSPTOT  ! Local patch loop

               IF(IPG == NSPG(NSP)) THEN

                  FOUNDL = .TRUE.

                  IDIMT  = IDIMS(NSP)
                  JDIMT  = JDIMS(NSP)
                  ITYPET = ICON((NSPP(NSP)-1)*IC9+1)
                  LOADGT = BOUNDF(NSPP(NSP))(1:80)

                  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ITSME, IERR)

               ENDIF

            ENDDO 

            CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,MPI_LOR,
     &                         MPI_COMM_WORLD,IERR)

            IF(FOUNDG) THEN
               NSPI = NSPI + 1
               CALL MPI_ALLREDUCE(IDIMT,IDIMSW(NSPI),1,
     &              MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
               CALL MPI_ALLREDUCE(JDIMT,JDIMSW(NSPI),1,
     &              MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
               CALL MPI_ALLREDUCE(ITYPET,ITYPEW(NSPI),1,
     &              MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
               CALL MPI_ALLREDUCE(ITSME,OWNER,1,
     &              MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
               IF(SOLIDL(ITYPEW(NSPI))) THEN ! Only solid like patches
                  LOADGG(NSPI) = LOADGT
                  CALL MPI_BCAST(LOADGG(NSPI),80,MPI_CHARACTER,OWNER,
     &                           MPI_COMM_WORLD,IERR)
               ENDIF
            ENDIF

         ENDDO


C ... Transfer the surface patch coordinates to the master process

         IG1G      = 1
         NSPI      = 0
         NSPOINTSG = 0 

         DO ISP = 1,NSPTOTG
            NSPOINTSG = NSPOINTSG + (IDIMSW(ISP)+1)*(JDIMSW(ISP)+1)
         ENDDO

         ALLOCATE(XCOW(NSPOINTSG),YCOW(NSPOINTSG),ZCOW(NSPOINTSG))

         XCOW = 0.0; YCOW = 0.0; ZCOW = 0.0 

         ALLOCATE(XCOP(NSPOINTSG),YCOP(NSPOINTSG),ZCOP(NSPOINTSG))

         XCOP = 0.0; YCOP = 0.0; ZCOP = 0.0 

         
         DO IPG = 1,NBCSG ! Global patch loop

C ... We know the size of the patch to be found next even when it is not 
C ... in this process.

            IF(NSPI == NSPTOTG) CYCLE  ! All patches found
 
            NSPOINTS = (IDIMSW(NSPI+1)+1)*(JDIMSW(NSPI+1)+1)
            ALLOCATE(XCOT(NSPOINTS),YCOT(NSPOINTS),ZCOT(NSPOINTS))
            XCOT = 0.0; YCOT = 0.0; ZCOT = 0.0 

            FOUNDL = .FALSE.

            DO NSP = 1,NSPTOT

               IF(IPG == NSPG(NSP)) THEN ! Global patch number

                  FOUNDL = .TRUE.

                  N     = NSPB(NSP) ! Local block number
                  IG1   = IG(1,N)

                  KX1   = KX1S(NSP)
                  KX2   = KX2S(NSP) + 1
                  KY1   = KY1S(NSP)
                  KY2   = KY2S(NSP) + 1
                  KZ1   = KZ1S(NSP)
                  KZ2   = KZ2S(NSP)

                  ISTR  = ISTRS(NSP)
                  JSTR  = JSTRS(NSP)
                  KSTR  = KSTRS(NSP)

                  ISP   = NSPP(NSP)              ! Local patch number
                  IWALL = ICON((ISP-1)*IC9 + 3)  ! Wall number

                  IF(IWALL == 2 .OR. IWALL == 5) THEN
                     ISTR = JSTRS(NSP)
                     JSTR = ISTRS(NSP)
                  ENDIF

                  IF(IWALL == 1 .OR. IWALL == 4) THEN
                     IN = JNRE
                     JN = KNRE
                     KN = INRE
                  ELSEIF(IWALL == 2 .OR. IWALL == 5) THEN
                     IN = INRE
                     JN = KNRE
                     KN = JNRE
                  ELSEIF(IWALL == 3 .OR. IWALL == 6) THEN
                     IN = INRE
                     JN = JNRE
                     KN = KNRE
                  ENDIF
              
                  IPT = 0
                  DO K=KZ1,KZ2
                     DO J=KY1,KY2
                        DO I=KX1,KX2
                           L = (I+IN-1)*ISTR+(J+JN-1)*JSTR+(K+KN-1)*KSTR
                           IPT = IPT + 1
                           XCOT(IPT) = XCO(IG1+L)
                           YCOT(IPT) = YCO(IG1+L)
                           ZCOT(IPT) = ZCO(IG1+L)
                        ENDDO
                     ENDDO
                  ENDDO
                
               ENDIF

            ENDDO 

            CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,MPI_LOR,
     &                            MPI_COMM_WORLD,IERR)
            IF(FOUNDG) THEN
               NSPI = NSPI + 1
               CALL MPI_REDUCE(XCOT,XCOW(IG1G),NSPOINTS,
     &              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(YCOT,YCOW(IG1G),NSPOINTS,
     &              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(ZCOT,ZCOW(IG1G),NSPOINTS,
     &              MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
               IG1G = IG1G + NSPOINTS
            ENDIF

            DEALLOCATE(XCOT,YCOT,ZCOT)

         ENDDO

      ELSE  ! XYZBND.FMT data in a single processor run


C ... Store the surface data into work arrays for later usage 
C ... instead of the original direct output.
       
         NSPTOTG = NSPTOT
         ALLOCATE(IDIMSW(NSPTOTG),JDIMSW(NSPTOTG),ITYPEW(NSPTOTG))
         ALLOCATE(IDIMSG(NSPTOTG),JDIMSG(NSPTOTG))
         ALLOCATE(IDIMST(NSPTOTG),JDIMST(NSPTOTG))
         ALLOCATE(LOADGG(NSPTOTG))

         LOADGG = ' '

         DO NSP = 1,NSPTOTG
            IDIMSW(NSP) = IDIMS(NSP)
            JDIMSW(NSP) = JDIMS(NSP)
            ITYPEW(NSP) = ICON((NSPP(NSP)-1)*IC9+1)
            IF(SOLIDL(ITYPEW(NSP)))  ! Only solid like patches
     &      LOADGG(NSP)(1:80) = BOUNDF(NSPP(NSP))(1:80)
         ENDDO
        
         NSPOINTSG = 0 

         DO ISP = 1,NSPTOTG
            NSPOINTSG = NSPOINTSG + (IDIMSW(ISP)+1)*(JDIMSW(ISP)+1)
         ENDDO

         ALLOCATE(XCOW(NSPOINTSG),YCOW(NSPOINTSG),ZCOW(NSPOINTSG))
         ALLOCATE(XCOP(NSPOINTSG),YCOP(NSPOINTSG),ZCOP(NSPOINTSG))

         IPT = 0

         DO NSP = 1,NSPTOT
            N     = NSPB(NSP)  ! Block number
            IG1   = IG(1,N)

            KX1   = KX1S(NSP)
            KX2   = KX2S(NSP) + 1
            KY1   = KY1S(NSP)
            KY2   = KY2S(NSP) + 1
            KZ1   = KZ1S(NSP)
            KZ2   = KZ2S(NSP)

            ISTR  = ISTRS(NSP)
            JSTR  = JSTRS(NSP)
            KSTR  = KSTRS(NSP)

            ISP   = NSPP(NSP)             ! Patch number
            IWALL = ICON((ISP-1)*IC9+3)   ! WALL NUMBER

            IF(IWALL == 2 .OR. IWALL == 5) THEN
               ISTR = JSTRS(NSP)
               JSTR = ISTRS(NSP)
            ENDIF

            IF(IWALL == 1 .OR. IWALL == 4) THEN
               IN = JNRE
               JN = KNRE
               KN = INRE
            ELSEIF(IWALL == 2 .OR. IWALL == 5) THEN
               IN = INRE
               JN = KNRE
               KN = JNRE
            ELSEIF(IWALL == 3 .OR. IWALL == 6) THEN
               IN = INRE
               JN = JNRE
               KN = KNRE
            ENDIF

         
C ... Save the data instead of direct writing.


            DO K=KZ1,KZ2
               DO J=KY1,KY2
                  DO I=KX1,KX2
                     L = (I+IN-1)*ISTR+(J+JN-1)*JSTR+(K+KN-1)*KSTR
                     IPT = IPT + 1
                     XCOW(IPT) = XCO(IG1+L)
                     YCOW(IPT) = YCO(IG1+L)
                     ZCOW(IPT) = ZCO(IG1+L)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO

      ENDIF  ! XYZBND.FMT surface data ready



C ... Write XYZBND.FMT

      IF(MASTER) THEN
         
         OPEN(29,FILE='XYZBND.FMT',FORM='FORMATTED')

C ... Select solid like pathes to be written to the surface load files.

         NSPTOTT = 0
         DO NSPI=1,NSPTOTG
            IF(SOLIDL(ITYPEW(NSPI))) THEN  ! Only solid like patches
               NSPTOTT = NSPTOTT + 1
               IDIMST(NSPTOTT) = IDIMSW(NSPI)
               JDIMST(NSPTOTT) = JDIMSW(NSPI)
            ENDIF
         ENDDO
 
         WRITE(29,*) NSPTOTT         
         KMAXXP = 1
         WRITE(29,*) (IDIMST(NSPI)+1,JDIMST(NSPI)+1,KMAXXP,
     &                NSPI=1,NSPTOTT)
         IG1G = 1

         DO NSPI=1,NSPTOTG
            II = (IDIMSW(NSPI)+1)*(JDIMSW(NSPI)+1)
            IF(SOLIDL(ITYPEW(NSPI))) ! Only solid like patches            
     &           WRITE(29,*)
     &           (XCOW(I),I=IG1G,IG1G+II-1), 
     &           (YCOW(I),I=IG1G,IG1G+II-1), 
     &           (ZCOW(I),I=IG1G,IG1G+II-1)
            IG1G = IG1G + II           

         ENDDO

         CLOSE(29)

      ENDIF

C ... XYZBND.FMT written


C ... Collect and write the patch function files

      IF(MASTER) THEN

         IF(SIMULATION_END) THEN
            OPEN(31,FILE = 'RESBND.FMT.nam',    FORM = 'FORMATTED')
            OPEN(72,FILE = 'RES.BIN.fvsrf',     FORM = 'UNFORMATTED')
            OPEN(73,FILE = 'RES.BIN.fvsrf.nam', FORM = 'FORMATTED')
            OPEN(74,FILE = 'LINER.BIN',         FORM = 'UNFORMATTED')
            OPEN(75,FILE = 'LINER.FMT',         FORM = 'FORMATTED')
            OPEN(76,FILE = 'LINER.nam',         FORM = 'FORMATTED')
            OPEN(400,FILE= 'LOADDISTS',         FORM = 'FORMATTED',
     &            STATUS = 'UNKNOWN',           RECL = 250)
         ELSE
            OPEN(72,FILE = 'MOVIE/RES'//FILEE//'.BIN.fvsrf',
     &              FORM = 'UNFORMATTED')
            OPEN(73,FILE = 'MOVIE/RES'//FILEE//'.BIN.fvsrf.nam',
     &              FORM = 'FORMATTED')
            OPEN(400,FILE= 'MOVIE/LOADDISTS'//FILEE//'.txt',
     &               STATUS='UNKNOWN',FORM='FORMATTED',RECL=250)
         ENDIF

      ENDIF
              

      KMAXXP = 1
      NOVAR  = 28
      NWOOD  = 5

      IF(MASTER .AND. SIMULATION_END) THEN
         WRITE(72)   NSPTOTG         
         WRITE(72)  (IDIMSW(NSPI),  JDIMSW(NSPI),  NOVAR,NSPI=1,NSPTOTG)
         WRITE(74)   NSPTOTG
         WRITE(74)  (IDIMSW(NSPI)+1,JDIMSW(NSPI)+1,NOVAR,NSPI=1,NSPTOTG)
         WRITE(75,*) NSPTOTG
         WRITE(75,*)(IDIMSW(NSPI)+1,JDIMSW(NSPI)+1,NOVAR,NSPI=1,NSPTOTG)
      ENDIF

      IF(MASTER .AND. .NOT.SIMULATION_END) THEN
         WRITE(72)   NSPTOTG         
         WRITE(72)  (IDIMSW(NSPI),  JDIMSW(NSPI),  NOVAR,NSPI=1,NSPTOTG)
      ENDIF

        
      IG1G      = 1
      NSPI      = 0
      NSPOINTSG = 0 

      DO ISP = 1,NSPTOTG
         NSPOINTSG = NSPOINTSG + IDIMSW(ISP)*JDIMSW(ISP)
      ENDDO

      ALLOCATE(SURFXW(NSPOINTSG),SURFYW(NSPOINTSG),SURFZW(NSPOINTSG),
     &         FRICXW(NSPOINTSG),FRICYW(NSPOINTSG),FRICZW(NSPOINTSG),
     &         TWALLW(NSPOINTSG),AREAW(NSPOINTSG))

      SURFXW = 0.0; SURFYW = 0.0; SURFZW = 0.0
      FRICXW = 0.0; FRICYW = 0.0; FRICZW = 0.0
      TWALLW = 0.0; AREAW  = 0.0

      ALLOCATE(SURFXP(NSPOINTSG),SURFYP(NSPOINTSG),SURFZP(NSPOINTSG),
     &         FRICXP(NSPOINTSG),FRICYP(NSPOINTSG),FRICZP(NSPOINTSG),
     &         TWALLP(NSPOINTSG),AREAP(NSPOINTSG))


      IF(.NOT.PARALLEL) NBCSG = NBCS
      
      DO 999 IPG = 1,NBCSG          ! Global patch loop

C ... We know the size of the patch to be found next even when it  
C ... is not in this process.

         IF(NSPI == NSPTOTG) CYCLE  ! All patches found

         NSPOINTS = IDIMSW(NSPI+1)*JDIMSW(NSPI+1)
         NTRAIN   = NOVAR*NSPOINTS

         ALLOCATE(TRAINFV(NTRAIN),TRAINFW(NTRAIN)) 

         TRAINFV = 0.0; TRAINFW = 0.0  

         FOUNDL = .FALSE.

         IPB = 0

         DO 996 NSP = 1,NSPTOT
            
            IF(IPG /= NSPG(NSP)) CYCLE

            FOUNDL = .TRUE.

            N      = NSPB(NSP)              ! Block number
            ISP    = NSPP(NSP)              ! Patch number
            IWALL  = ICON((ISP-1)*IC9 + 3)  ! Wall number

            IG1    = IG(1,N)
            IF1    = IHF(ISP,1)

C ... In the case of the j-direction the data is stored in in jki-order (sigh)

            IF(IWALL == 2 .OR. IWALL == 5) THEN
               KY1   = KX1S(NSP)
               KY2   = KX2S(NSP)
               KX1   = KY1S(NSP)
               KX2   = KY2S(NSP)
               IPSTR = KX2-KX1 + 1 + 2*LN
               JPSTR = 1
               NN    = (LN-KY1)*IPSTR - KX1 + LN + IF1
            ENDIF               ! Patch index with LN ghost cells

C ... but is printed for Fieldview in a jik order

            KX1 = KX1S(NSP)
            KX2 = KX2S(NSP)
            KY1 = KY1S(NSP)
            KY2 = KY2S(NSP) 

C ... and the order in i and k-directions is 'normal'

            IF(IWALL /= 2 .AND. IWALL /= 5) THEN
               IPSTR = 1
               JPSTR = KX2-KX1 + 1 + 2*LN
               NN    = (LN-KY1)*JPSTR - KX1 + LN + IF1
            ENDIF               ! Patch index with LN ghost cells

C ... Collect the FIELDVIEW function file for the patches

            DO J=KY1,KY2
               DO I=KX1,KX2

                  IPB = IPB + 1
                  L   = I*IPSTR+J*JPSTR+NN

                  TRAINFV(IPB              ) = CPWALL(L)/REFPRE
                  TRAINFV(IPB +  1*NSPOINTS) = QWALL(L)
                  TRAINFV(IPB +  2*NSPOINTS) = QWFRIC(L)
                  TRAINFV(IPB +  3*NSPOINTS) = TWALL(L)
                  TRAINFV(IPB +  4*NSPOINTS) = TAUW1(L)/REFPRE
                  TRAINFV(IPB +  5*NSPOINTS) = TAUW2(L)/REFPRE
                  TRAINFV(IPB +  6*NSPOINTS) = TAUWX(L)/REFPRE
                  TRAINFV(IPB +  7*NSPOINTS) = TAUWY(L)/REFPRE
                  TRAINFV(IPB +  8*NSPOINTS) = TAUWZ(L)/REFPRE
                  TRAINFV(IPB +  9*NSPOINTS) = SURFX(L)
                  TRAINFV(IPB + 10*NSPOINTS) = SURFY(L)
                  TRAINFV(IPB + 11*NSPOINTS) = SURFZ(L)
                  TRAINFV(IPB + 12*NSPOINTS) = HFLUX(L)
                  TRAINFV(IPB + 13*NSPOINTS) = SURLE(L)
                  TRAINFV(IPB + 14*NSPOINTS) = DSURLE(L)
                  TRAINFV(IPB + 15*NSPOINTS) = F1RK(L)
                  TRAINFV(IPB + 16*NSPOINTS) = WMFLUX(L)
                  TRAINFV(IPB + 17*NSPOINTS) = BOUNR(L)
                  TRAINFV(IPB + 18*NSPOINTS) = BOUNT(L)
                  TRAINFV(IPB + 19*NSPOINTS) = BOUNP(L)
                  TRAINFV(IPB + 20*NSPOINTS) = BOUNRK(L)
                  TRAINFV(IPB + 21*NSPOINTS) = BOUNEP(L)
                  TRAINFV(IPB + 22*NSPOINTS) = BOUNMF(L)
                  TRAINFV(IPB + 23*NSPOINTS) = BOUNU(L)
                  TRAINFV(IPB + 24*NSPOINTS) = BOUNV(L)
                  TRAINFV(IPB + 25*NSPOINTS) = BOUNW(L)
                  TRAINFV(IPB + 26*NSPOINTS) = BOUNE(L)
                  TRAINFV(IPB + 27*NSPOINTS) = WTRAN(L)
               ENDDO
            ENDDO

 996     CONTINUE


C ... Write one FIELDVIEW function file patch

         IF(PARALLEL) THEN
            CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,MPI_LOR,
     &                         MPI_COMM_WORLD,IERR)
         ELSE
            FOUNDG = FOUNDL
         ENDIF
         IF(FOUNDG) THEN
            IF(PARALLEL) THEN
               CALL MPI_REDUCE(TRAINFV,TRAINFW,NTRAIN,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
            ELSE
               TRAINFW = TRAINFV 
            ENDIF
            IF(MASTER) WRITE(72) (REAL(TRAINFW(I),4),I=1,NTRAIN)
         ENDIF

         DEALLOCATE(TRAINFV, TRAINFW)
      
C ... One FIELDVIEW function file patch written


         ALLOCATE(SURFXT(NSPOINTS),SURFYT(NSPOINTS),SURFZT(NSPOINTS),
     &            FRICXT(NSPOINTS),FRICYT(NSPOINTS),FRICZT(NSPOINTS),
     &            TWALLT(NSPOINTS),AREAT (NSPOINTS))

         SURFXT = 0.0; SURFYT = 0.0; SURFZT = 0.0
         FRICXT = 0.0; FRICYT = 0.0; FRICZT = 0.0
         TWALLT = 0.0; AREAT  = 0.0

         FOUNDL = .FALSE.

         IPA = 0
         IPT = 0

         DO 997 NSP = 1,NSPTOT

            IF(IPG /= NSPG(NSP)) CYCLE

            FOUNDL = .TRUE.

            KX1W = KX1S(NSP)
            KX2W = KX2S(NSP)
            KY1W = KY1S(NSP)
            KY2W = KY2S(NSP)
            KZ1W = KZ1S(NSP)
            KZ2W = KZ2S(NSP)

            ISTR = ISTRS(NSP)
            JSTR = JSTRS(NSP)
            KSTR = KSTRS(NSP)


            IF(IWALL == 2 .OR. IWALL == 5) THEN
               ISTR = JSTRS(NSP)
               JSTR = ISTRS(NSP)
            ENDIF

            IF(IWALL == 1 .OR. IWALL == 4) THEN
               IN = JNRE
               JN = KNRE
               KN = INRE
            ELSEIF(IWALL == 2 .OR. IWALL == 5) THEN
               IN = INRE
               JN = KNRE
               KN = JNRE
            ELSEIF(IWALL == 3 .OR. IWALL == 6) THEN
               IN = INRE
               JN = JNRE
               KN = KNRE
            ENDIF
      
            DO J=KY1,KY2
               DO I=KX1,KX2
                  IPT = IPT + 1
                  L   = I*IPSTR+J*JPSTR+NN
                  SURFXT(IPT) = SURFX(L)
                  SURFYT(IPT) = SURFY(L)
                  SURFZT(IPT) = SURFZ(L)
                  FRICXT(IPT) = TAUWX(L)
                  FRICYT(IPT) = TAUWY(L)
                  FRICZT(IPT) = TAUWZ(L)
                  TWALLT(IPT) = TWALL(L)
               ENDDO
            ENDDO

            DO K=KZ1W,KZ2W
               DO J=KY1W,KY2W
                  DO I=KX1W,KX2W
                     IPA = IPA + 1
                     LA=(I+IN-1)*ISTR+(J+JN-1)*JSTR+(K+KN-1)*KSTR
                     LA = IG1 + LA
                     IF(IWALL == 1 .OR. IWALL == 4) AREAT(IPA) = A1(LA)
                     IF(IWALL == 2 .OR. IWALL == 5) AREAT(IPA) = A2(LA)
                     IF(IWALL == 3 .OR. IWALL == 6) AREAT(IPA) = A3(LA)
                  ENDDO
               ENDDO
            ENDDO

 997     CONTINUE

         IF(PARALLEL) THEN
            CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,MPI_LOR,
     &                         MPI_COMM_WORLD,IERR)
         ELSE
            FOUNDG = FOUNDL
         ENDIF
               
         IF(FOUNDG) THEN
            IF(PARALLEL) THEN
               CALL MPI_REDUCE(SURFXT,SURFXW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(SURFYT,SURFYW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(SURFZT,SURFZW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(FRICXT,FRICXW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(FRICYT,FRICYW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(FRICZT,FRICZW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(TWALLT,TWALLW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
               CALL MPI_REDUCE(AREAT,AREAW(IG1G),NSPOINTS,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
            ELSE
               SURFXW(IG1G:IG1G+NSPOINTS-1) = SURFXT(1:NSPOINTS)
               SURFYW(IG1G:IG1G+NSPOINTS-1) = SURFYT(1:NSPOINTS)
               SURFZW(IG1G:IG1G+NSPOINTS-1) = SURFZT(1:NSPOINTS)
               FRICXW(IG1G:IG1G+NSPOINTS-1) = FRICXT(1:NSPOINTS)
               FRICYW(IG1G:IG1G+NSPOINTS-1) = FRICYT(1:NSPOINTS)
               FRICZW(IG1G:IG1G+NSPOINTS-1) = FRICZT(1:NSPOINTS)
               TWALLW(IG1G:IG1G+NSPOINTS-1) = TWALLT(1:NSPOINTS)
                AREAW(IG1G:IG1G+NSPOINTS-1) =  AREAT(1:NSPOINTS)
            ENDIF
            IG1G = IG1G + NSPOINTS
         ENDIF

         DEALLOCATE(SURFXT,SURFYT,SURFZT,TWALLT,AREAT)
         DEALLOCATE(FRICXT,FRICYT,FRICZT)



C *** Liner files **********************************************************

C ... In 'PLOTPI-format' the indeces run differently (in j-direction as jki)

         NSPOINTC = (IDIMSW(NSPI+1)+1)*(JDIMSW(NSPI+1)+1)
         NTRAIN   =  NOVAR*NSPOINTC
         ALLOCATE(TRAINLI(NTRAIN), TRAINLW(NTRAIN))
         TRAINLI = 0.0; TRAINLW = 0.0 


         FOUNDL = .FALSE.

         IPC = 0


C ... Temporary arrays for output

         CPWALL_ = CPWALL
         QWALL_  = QWALL
         QWFRIC_ = QWFRIC
         TWALL_  = TWALL
         TAUW1_  = TAUW1
         TAUW2_  = TAUW2
         TAUWX_  = TAUWX
         TAUWY_  = TAUWY
         TAUWZ_  = TAUWZ
         SURFX_  = SURFX
         SURFY_  = SURFY
         SURFZ_  = SURFZ
         HFLUX_  = HFLUX
         SURLE_  = SURLE
         DSURLE_ = DSURLE
         F1RK_   = F1RK
         WMFLUX_ = WMFLUX
         BOUNR_  = BOUNR
         BOUNT_  = BOUNT
         BOUNP_  = BOUNP
         BOUNRK_ = BOUNRK
         BOUNEP_ = BOUNEP
         BOUNMF_ = BOUNMF
         BOUNU_  = BOUNU
         BOUNV_  = BOUNV
         BOUNW_  = BOUNW
         BOUNE_  = BOUNE
         WTRAN_  = WTRAN

                  
         DO 998 NSP = 1,NSPTOT

            IF(IPG /= NSPG(NSP)) CYCLE

            FOUNDL = .TRUE.

            N     = NSPB(NSP)   ! Block number
            ISP   = NSPP(NSP)   ! Patch number
            IWALL = ICON((ISP-1)*IC9 + 3) ! WALL NUMBER

            IG1 = IG(1,N)
            IF1 = IHF(ISP,1)

            KX1 = KX1S(NSP)
            KX2 = KX2S(NSP)
            KY1 = KY1S(NSP)
            KY2 = KY2S(NSP) 
            
            IF(IWALL == 2 .OR. IWALL == 5) THEN
               KY1 = KX1S(NSP)
               KY2 = KX2S(NSP)
               KX1 = KY1S(NSP)
               KX2 = KY2S(NSP)
            ENDIF

C ... Transfer into cell-vertex format (use temporary arrays)

            CALL REDPAT(CPWALL_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( QWALL_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT(QWFRIC_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TWALL_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TAUW1_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TAUW2_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TAUWX_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TAUWY_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( TAUWZ_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( SURFX_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( SURFY_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( SURFZ_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( HFLUX_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( SURLE_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT(DSURLE_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT(  F1RK_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT(WMFLUX_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            CALL REDPAT( BOUNR_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( BOUNT_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( BOUNP_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT(BOUNRK_(IF1),ZZZ,KX1,KX2,KY1,KY2) 
            CALL REDPAT(BOUNEP_(IF1),ZZZ,KX1,KX2,KY1,KY2) 
            CALL REDPAT(BOUNMF_(IF1),ZZZ,KX1,KX2,KY1,KY2) 
            CALL REDPAT( BOUNU_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( BOUNV_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( BOUNW_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( BOUNE_(IF1),ZZZ,KX1,KX2,KY1,KY2)  
            CALL REDPAT( WTRAN_(IF1),ZZZ,KX1,KX2,KY1,KY2)
            
            IF(IWALL == 2 .OR. IWALL == 5) THEN
               IPSTR  = KX2-KX1 + 1 + 1
               JPSTR  = 1
            ENDIF

C ... but is printed for Liner in a jik order

            KX1 = KX1S(NSP)
            KX2 = KX2S(NSP)       
            KY1 = KY1S(NSP)
            KY2 = KY2S(NSP)       

C ... and the order in i and k-directions is 'normal'

            IF(IWALL /= 2 .AND. IWALL /= 5) THEN
               IPSTR  = 1
               JPSTR  = KX2-KX1 + 1 + 1
            ENDIF
            
            IMAXP1 = KX2 - KX1 + 2
            JMAXP1 = KY2 - KY1 + 2
            NN     = -IPSTR-JPSTR + IF1


            DO J=1,JMAXP1
               DO I=1,IMAXP1
                  IPC = IPC + 1
                  L   = I*IPSTR+J*JPSTR+NN
                  
                  TRAINLI(IPC              ) = CPWALL_(L)/REFPRE
                  TRAINLI(IPC +  1*NSPOINTC) = QWALL_(L)
                  TRAINLI(IPC +  2*NSPOINTC) = QWFRIC_(L)
                  TRAINLI(IPC +  3*NSPOINTC) = TWALL_(L)
                  TRAINLI(IPC +  4*NSPOINTC) = TAUW1_(L)/REFPRE
                  TRAINLI(IPC +  5*NSPOINTC) = TAUW2_(L)/REFPRE
                  TRAINLI(IPC +  6*NSPOINTC) = TAUWX_(L)/REFPRE
                  TRAINLI(IPC +  7*NSPOINTC) = TAUWY_(L)/REFPRE
                  TRAINLI(IPC +  8*NSPOINTC) = TAUWZ_(L)/REFPRE
                  TRAINLI(IPC +  9*NSPOINTC) = SURFX_(L)
                  TRAINLI(IPC + 10*NSPOINTC) = SURFY_(L)
                  TRAINLI(IPC + 11*NSPOINTC) = SURFZ_(L)
                  TRAINLI(IPC + 12*NSPOINTC) = HFLUX_(L)
                  TRAINLI(IPC + 13*NSPOINTC) = SURLE_(L)
                  TRAINLI(IPC + 14*NSPOINTC) = DSURLE_(L)
                  TRAINLI(IPC + 15*NSPOINTC) = F1RK_(L)
                  TRAINLI(IPC + 16*NSPOINTC) = WMFLUX_(L)
                  TRAINLI(IPC + 17*NSPOINTC) = BOUNR_(L)
                  TRAINLI(IPC + 18*NSPOINTC) = BOUNT_(L)
                  TRAINLI(IPC + 19*NSPOINTC) = BOUNP_(L)
                  TRAINLI(IPC + 20*NSPOINTC) = BOUNRK_(L)
                  TRAINLI(IPC + 21*NSPOINTC) = BOUNEP_(L)
                  TRAINLI(IPC + 22*NSPOINTC) = BOUNMF_(L)
                  TRAINLI(IPC + 23*NSPOINTC) = BOUNU_(L)
                  TRAINLI(IPC + 24*NSPOINTC) = BOUNV_(L)
                  TRAINLI(IPC + 25*NSPOINTC) = BOUNW_(L)
                  TRAINLI(IPC + 26*NSPOINTC) = BOUNE_(L)
                  TRAINLI(IPC + 27*NSPOINTC) = WTRAN_(L)
               ENDDO
            ENDDO

 998     CONTINUE


C ... Write one Liner file patch

         IF(PARALLEL) THEN
            CALL MPI_ALLREDUCE(FOUNDL,FOUNDG,1,MPI_LOGICAL,MPI_LOR,
     &                         MPI_COMM_WORLD,IERR)
         ELSE
            FOUNDG = FOUNDL
         ENDIF

         IF(FOUNDG) THEN
            NSPI = NSPI + 1
            IF(PARALLEL) THEN
               CALL MPI_REDUCE(TRAINLI,TRAINLW,NTRAIN,MPI_REAL8,
     &                         MPI_SUM,0,MPI_COMM_WORLD,IERR)
            ELSE
               TRAINLW = TRAINLI
            ENDIF
            IF(MASTER .AND. SIMULATION_END) THEN
               WRITE(74) (REAL(TRAINLW(I),4),I=1,NTRAIN)
               WRITE(75,9600) (TRAINLW(I),I=1,NTRAIN)
            ENDIF
         ENDIF

C ... One Liner file patch written

         DEALLOCATE(TRAINLI,TRAINLW)

 999  CONTINUE

 9600 FORMAT(5E15.6)


C ... Write RESBND.FMT

      IF(MASTER) THEN

         IF(SIMULATION_END) THEN ! Write only at the end of simulation
         
            OPEN(30,FILE='RESBND.FMT',FORM='FORMATTED')

            WRITE(30,*) NSPTOTT         
            WRITE(30,*) (IDIMST(NSPI),JDIMST(NSPI),NWOOD,NSPI=1,NSPTOTT)

            IG1G = 1

            DO NSPI=1,NSPTOTG

               II = IDIMSW(NSPI)*JDIMSW(NSPI)

               IF(SOLIDL(ITYPEW(NSPI))) ! Only solid like patches
     &              WRITE(30,9600)
     &              (SURFXW(I),I=IG1G,IG1G+II-1), 
     &              (SURFYW(I),I=IG1G,IG1G+II-1), 
     &              (SURFZW(I),I=IG1G,IG1G+II-1), 
     &              (TWALLW(I),I=IG1G,IG1G+II-1), 
     &              (AREAW(I), I=IG1G,IG1G+II-1) 
               IG1G = IG1G + II           

            ENDDO

            CLOSE(30)

         ENDIF  ! Write only at the end of simulation

      ENDIF

C ... RESBND.FMT written



C ... Compute the load distributions for the force groups

      IF(.NOT.PARALLEL .OR. (PARALLEL .AND. MASTER .AND. GROUP)) THEN

         WRITE(400,*)
         WRITE(400,*) "LOAD DISTRIBUTIONS FOR SELECTED FORCE GROUPS"
         WRITE(400,*) "============================================"
         WRITE(400,*)
         WRITE(400,*) "Free stream dynamic pressure: ",
     &                 REAL(REFPRE,4), "N/m^2"
         WRITE(400,*) "Reference area:   ", REAL(AREF,4), "m^2"
         WRITE(400,*) "Reference length: ", REAL(CHLREF,4), "m"

         DO J=1,52
     
            FOUND = .FALSE.

            IF(J <= 26) GNAME = CHAR(J+64) ! Uppercase characters
            IF(J > 26) GNAME = CHAR(J+70) ! Lowercase characters

            NSG    = 0
            IWGRID = 1
            IGGRID = 1
            IWSOLU = 1
            IGSOLU = 1
             
            DO I=1,NSPTOTG

               IF(LOADGG(I)(J:J) /= ' ') THEN

                  FOUND = .TRUE.
                  NSG = NSG + 1

                  IDIMSG(NSG) = IDIMSW(I)
                  JDIMSG(NSG) = JDIMSW(I)

                  XCOP(IGGRID:IGGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)
     &          = XCOW(IWGRID:IWGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)
                  YCOP(IGGRID:IGGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)
     &          = YCOW(IWGRID:IWGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)
                  ZCOP(IGGRID:IGGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)
     &          = ZCOW(IWGRID:IWGRID+(IDIMSW(I)+1)*(JDIMSW(I)+1)-1)

                  SURFXP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = SURFXW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
                  SURFYP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = SURFYW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
                  SURFZP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = SURFZW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)

                  FRICXP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = FRICXW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
                  FRICYP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = FRICYW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
                  FRICZP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = FRICZW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)

                  TWALLP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          = TWALLW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
                   AREAP(IGSOLU:IGSOLU+IDIMSW(I)*JDIMSW(I)-1)
     &          =  AREAW(IWSOLU:IWSOLU+IDIMSW(I)*JDIMSW(I)-1)
         
                  IGGRID = IGGRID + (IDIMSW(I)+1)*(JDIMSW(I)+1)
                  IGSOLU = IGSOLU +  IDIMSW(I)   * JDIMSW(I)

               ENDIF
              
               IWGRID = IWGRID + (IDIMSW(I)+1)*(JDIMSW(I)+1)
               IWSOLU = IWSOLU +  IDIMSW(I)   * JDIMSW(I)

            ENDDO

            IF(FOUND) THEN
 
             WRITE(400,*)
             WRITE(400,*)
             WRITE(400,*) "########################################",
     &                    "##############################"
             WRITE(400,321) GNAME,FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24)
             WRITE(400,*) "########################################",
     &                    "##############################"
             CALL LOADD(NSG,IDIMSG,JDIMSG,XCOP,YCOP,ZCOP,
     &                  SURFXP,SURFYP,SURFZP,FRICXP,FRICYP,FRICZP,
     &                  TWALLP,AREAP,AREF,CHLREF,REFPRE)

            ENDIF

         ENDDO

         CLOSE(400)

      ENDIF

 321  FORMAT(' Force group: ',1A1,2X,A24)


      DEALLOCATE(IDIMSW,JDIMSW,ITYPEW)        
      DEALLOCATE(IDIMSG,JDIMSG)        
      DEALLOCATE(IDIMST,JDIMST)        
      DEALLOCATE(XCOW,YCOW,ZCOW)       
      DEALLOCATE(SURFXW,SURFYW,SURFZW) 
      DEALLOCATE(FRICXW,FRICYW,FRICZW) 
      DEALLOCATE(TWALLW,AREAW)         
      DEALLOCATE(XCOP,YCOP,ZCOP)       
      DEALLOCATE(SURFXP,SURFYP,SURFZP) 
      DEALLOCATE(FRICXP,FRICYP,FRICZP) 
      DEALLOCATE(TWALLP,AREAP)         
      DEALLOCATE(LOADGG)               


      IF(MASTER) THEN

      WRITE(73,"(A)") 'Cp_wall'
      WRITE(73,"(A)") 'Q_wall'
      WRITE(73,"(A)") 'Qf_wall'
      WRITE(73,"(A)") 'T_wall'
      WRITE(73,"(A)") 'Cf_1'
      WRITE(73,"(A)") 'Cf_2'
      WRITE(73,"(A)") 'Cf_x;Cf_wall'
      WRITE(73,"(A)") 'Cf_y'
      WRITE(73,"(A)") 'Cf_z'
      WRITE(73,"(A)") 'Fx_wall (u with MIR);F_wall'
      WRITE(73,"(A)") 'Fy_wall (v with MIR)'
      WRITE(73,"(A)") 'Fz_wall (w with MIR)'
      WRITE(73,"(A)") 'Qt_wall'
      WRITE(73,"(A)") 'Surf1'
      WRITE(73,"(A)") 'Surf2'
      WRITE(73,"(A)") 'y+(1)'
      WRITE(73,"(A)") 'Mflux'
      WRITE(73,"(A)") 'Dist_density'
      WRITE(73,"(A)") 'Dist_temp'
      WRITE(73,"(A)") 'Dist_pre'
      WRITE(73,"(A)") 'Dist_rk'
      WRITE(73,"(A)") 'Dist_re'
      WRITE(73,"(A)") 'Dist_Mass_flux'
      WRITE(73,"(A)") 'Dist_Momentum_x'
      WRITE(73,"(A)") 'Dist_Momentum_y'
      WRITE(73,"(A)") 'Dist_Momentum_z'
      WRITE(73,"(A)") 'Dist_Energy'
      WRITE(73,"(A)") 'Wtran'

      IF(SIMULATION_END) THEN  ! Write only at the end of simulation

      WRITE(76,"(A)") 'Cp_wall'
      WRITE(76,"(A)") 'Q_wall'
      WRITE(76,"(A)") 'Qf_wall'
      WRITE(76,"(A)") 'T_wall'
      WRITE(76,"(A)") 'Cf_1'
      WRITE(76,"(A)") 'Cf_2'
      WRITE(76,"(A)") 'Cf_x;Cf_wall'
      WRITE(76,"(A)") 'Cf_y'
      WRITE(76,"(A)") 'Cf_z'
      WRITE(76,"(A)") 'Fx_wall (u with MIR);F_wall'
      WRITE(76,"(A)") 'Fy_wall (v with MIR)'
      WRITE(76,"(A)") 'Fz_wall (w with MIR)'
      WRITE(76,"(A)") 'Qt_wall'
      WRITE(76,"(A)") 'Surf1'
      WRITE(76,"(A)") 'Surf2'
      WRITE(76,"(A)") 'y+(1)'
      WRITE(76,"(A)") 'Mflux'
      WRITE(76,"(A)") 'Dist_density'
      WRITE(76,"(A)") 'Dist_temp'
      WRITE(76,"(A)") 'Dist_pre'
      WRITE(76,"(A)") 'Dist_rk'
      WRITE(76,"(A)") 'Dist_re'
      WRITE(76,"(A)") 'Dist_Mass_flux'
      WRITE(76,"(A)") 'Dist_Momentum_x'
      WRITE(76,"(A)") 'Dist_Momentum_y'
      WRITE(76,"(A)") 'Dist_Momentum_z'
      WRITE(76,"(A)") 'Dist_Energy'
      WRITE(76,"(A)") 'Wtran'
      
      WRITE(31,"(A)") 'Fx_wall;F_wall'
      WRITE(31,"(A)") 'Fy_wall'
      WRITE(31,"(A)") 'Fz_wall'
      WRITE(31,"(A)") 'T_wall'
      WRITE(31,"(A)") 'Area'

      ENDIF  ! Write only at the end of simulation

      CLOSE(31)
      CLOSE(72)
      CLOSE(73)
      CLOSE(74)
      CLOSE(75)
      CLOSE(76)

      ENDIF


      RETURN
      END SUBROUTINE FVSRF
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      LOGICAL FUNCTION SOLIDL(ITYPE)

C ... Function for selecting solid like surface patches. Current choices
C ... are SOL, ROT, MOV, HTS and SLI. Note that you must not activate
C ... INL and OUT (ITYPEs 3 and 5) here.  

      IMPLICIT NONE

      INTEGER :: ITYPE
     
      SOLIDL = (ITYPE == 8  .OR. ITYPE == 9 .OR. ITYPE == 10 .OR.
     &          ITYPE == 15 .OR. ITYPE == 16) ! BCP

      RETURN
      END FUNCTION SOLIDL
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTTUR

      USE CHARACTERS

      USE INTEGERS, ONLY : MAXB

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,IG,RO,RK,REPS,DDEPS,
     &    FUN1,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,IDI1,IDI2,IDI3,ICON,
     &    NPATCH,F1R,TTS,IC

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: N,M,I1,J1,K1,L,IG1,IC1,IMAXP1,JMAXP1,KMAXP1,NMAX

C ... Modification of turbulence variables for output

      DO N = 1,NBLOCK

         M   = 1
         I1  = IMAX(M,N)
         J1  = JMAX(M,N)
         K1  = KMAX(M,N)
         IG1 = IG(M,N)
         IC1 = IC(1,N)

         IF(IDIS == 2) THEN  ! Reduce wall distance and blending function
        
            CALL VIEWSC(FUN1(IG1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &           IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),KBOT(1,N),
     &           KTOP(1,N),IDI1(N),IDI2(N),IDI3(N),IN,JN,KN,
     &           ICON(IC1),NPATCH(N),3,MAXB,1,0)

*            CALL REDUCE(FUN1(IG1),F1R(IG1),I1,J1,K1)

         ENDIF

C ... Reduce turbulence time scale

*         CALL REDUCE(TTS(IG1),F1R(IG1),I1,J1,K1)

         IG1    = IG(1,N)   - 1
*         IMAXP1 = IMAX(1,N) + 1
*         JMAXP1 = JMAX(1,N) + 1
*         KMAXP1 = KMAX(1,N) + 1
*         NMAX   = IMAXP1*JMAXP1*KMAXP1
         NMAX   = (IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)*(KMAX(1,N)+2*KN)

C ... CALCULATE TRUE EPSILON TO SRK

         IF (IDIS == 2) THEN    ! k-omega

            DO L = 1,NMAX
               DDEPS(IG1+L) = 0.09*REPS(IG1+L)*RK(IG1+L)/RO(IG1+L)
            ENDDO

         ELSE IF (IDIS == 1) THEN

            DO L = 1,NMAX       ! k-epsilon
               DDEPS(IG1+L) = DDEPS(IG1+L) + REPS(IG1+L)
            ENDDO

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE OUTTUR
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTSCL(FILEE,FI)

C ... This subroutine is obsolete.

      USE CHARACTERS

      USE INTEGERS,    ONLY : MAXB,MAXSB,IPRO

      USE TYPE_ARRAYS

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,IC,IBOT,JBOT,KBOT,
     &    ITOP,JTOP,KTOP,ICON,NPATCH,NPROCE,IT,IL,IK,IG,RO,RM,RN,
     &    RW,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,A1XA,A1YA,A1ZA,A2XA,
     &    A2YA,A2ZA,A3XA,A3YA,A3ZA,F1R,F1RM,JLOC,JTRA,APP,ZZZ,
     &    PRO,VAR,TRM

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: III,L,N,NS,IG1,IC1,NGL,IG3,IMAXP1,JMAXP1,KMAXP1,NMAX

      CHARACTER(LEN=3) ::  FILEE

      REAL :: FI(MAXSB,*)

      TYPE(PRE_COR) :: PRC

      IF (NSCAL == 0) RETURN

      OPEN(24,FILE=FILEE//'.BIN'//PRN,FORM='UNFORMATTED')
      IF(IPRO == 1) OPEN(71,FILE='SCL.nam',FORM='FORMATTED')

      WRITE(24) NBLOCK
      
C ... SCALAR OUTPUT FOR THE THREE-DIMENSIONAL SOLVER PPR 7.3

      IF (FILEE /= 'SCL') THEN

         DO III = 1,2           ! Try to fill the ghost cells correctly
            DO NS = 1,NSCAL
               CALL MIR(RO,RM,RN,RW,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,
     &              PRO,VAR,TRM,MAXSB,
     &              A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,
     &              A3YA,A3ZA,1,1,NBLOCK,1.,1,PRC)
               CALL CONEC(FI(1,NS),F1R,F1RM,JLOC,JTRA,APP,1,1,
     &              NBLOCK,0,0)
            ENDDO
         ENDDO

         DO N = 1,NBLOCK

            IG1     = IG(1,N)
            IC1     = IC(1,N)

            CALL VIEWSC(FI(IG1,1),IMAX(1,N),JMAX(1,N),KMAX(1,N),
     &           IBOT(1,N),ITOP(1,N),JBOT(1,N),JTOP(1,N),
     &           KBOT(1,N),KTOP(1,N),
     &           0,0,0,IN,JN,KN,ICON(IC1),NPATCH(N),0,MAXSB,NSCAL,0)

            NGL = NPROCE(1+N,IPRO) ! Global block number

            IF(N <= NBLOCK) THEN

               DO NS = 1,NSCAL
                  IG3 = IG1 + (NS-1)*MAXB
                  CALL PRINYS(22,FILEE//'no= '//CHAR(NS+48),FI(IG1,NS),
     &                 IT(1,N),IL(1,N),0,IK(1,N),
     &                 IMAX(1,N),JMAX(1,N),KMAX(1,N),NGL,1)
               ENDDO

            ENDIF

         ENDDO

      ENDIF
      
      DO NS = 1,NSCAL
         DO N  = 1,NBLOCK      
            IG1 = IG(1,N)
            CALL REDUCE(FI(IG1,NS),ZZZ,IMAX(1,N),JMAX(1,N),KMAX(1,N))
         ENDDO
      ENDDO

C ... Function file for scalars (not tested 24.9.02)

      IF(ITURB > 19 .AND. FILEE == 'SCL') THEN
         WRITE(24) (IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,NSCAL+1-7,
     +   N = 1,NBLOCK)
      ELSEIF(NSCAL /= 0 .AND. FILEE == 'SCL') THEN
         WRITE(24) (IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,NSCAL,
     +   N = 1,NBLOCK)
      ENDIF

      IF (FILEE /= 'SCL') THEN
         WRITE(24) (IMAX(1,N)+1,JMAX(1,N)+1,KMAX(1,N)+1,6,N = 1,NBLOCK)
      ENDIF

      DO N = 1,NBLOCK

         IG1    = IG(1,N)   - 1
         IMAXP1 = IMAX(1,N) + 1
         JMAXP1 = JMAX(1,N) + 1
         KMAXP1 = KMAX(1,N) + 1
         NMAX   = IMAXP1*JMAXP1*KMAXP1

         IF(ITURB > 19 .AND. FILEE == 'SCL') THEN
            WRITE(24) ((FI(IG1+L,NS),L=1,NMAX),NS=7,NSCAL)
         ELSEIF(NSCAL /= 0 .AND. FILEE == 'SCL') THEN
            WRITE(24) ((FI(IG1+L,NS),L=1,NMAX),NS=1,NSCAL)
         ENDIF

         IF (FILEE /= 'SCL') THEN
            WRITE(24) ((FI(IG1+L,NS),L=1,NMAX),NS=1,6)
         ENDIF
      
      ENDDO

      CLOSE(24)

      RETURN
      END SUBROUTINE OUTSCL
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WriteEnSightGeometry(name,iunit,nchim,ib,nb,
     &           idim,jdim,kdim,x,y,z,blank)

C ... Write EnSight Gold geometry file.

      implicit none

      character(len=*)  :: name
      character(len=80) :: buffer
      character(len=5)  :: block_number
      integer :: i, n, ib, nb, nchim, iunit, idim, jdim, kdim
      real, dimension(*) :: x, y, z, blank

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         buffer = "Fortran Binary"
         write(iunit) buffer
         buffer = "Geometry file from FINFLO output"
         write(iunit) buffer
         buffer = "EnSight Gold"
         write(iunit) buffer
         buffer = "node id off"
         write(iunit) buffer
         buffer = "element id off"
         write(iunit) buffer

      endif

      buffer = "part"
      write(iunit) buffer

      write(iunit) ib

      call numcha(block_number,ib)

      buffer = "FINFLO block #"//block_number
      write(iunit) buffer

      if(nchim > 0) then
         buffer = "block iblanked"
         write(iunit) buffer
      else
         buffer = "block"
         write(iunit) buffer
      endif      

c ... write block dimensions

      write(iunit) idim, jdim, kdim

      n = idim*jdim*kdim

c ... write node coordinates

      write(iunit) (real(X(I),4),i=1,n)
      write(iunit) (real(Y(I),4),i=1,n)
      write(iunit) (real(Z(I),4),i=1,n)

c ... write iblanking

      if(nchim > 0) write(iunit) (NINT(blank(i)),i=1,n)

      if(ib == nb) close(iunit)

      return 
      end subroutine WriteEnSightGeometry		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WriteEnSightScalar(name,iunit,ib,nb,
     &           idim,jdim,kdim,s)

C ... Write EnSight Gold per_node scalar file.

      implicit none

      character(len=*)  :: name
      character(len=80) :: buffer
      integer :: i, n, ib, nb, iunit, idim, jdim, kdim
      real, dimension(*) :: s

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         buffer = "EnSight Gold per_node scalar file"
         write(iunit) buffer

      endif

      buffer = "part"
      write(iunit) buffer

      write(iunit) ib

      buffer = "block"
      write(iunit) buffer

      n = idim*jdim*kdim

      write(iunit) (real(s(i),4),i=1,n)

      if(ib == nb) close(iunit)

      return 
      end subroutine WriteEnSightScalar		
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WriteEnSightVector(name,iunit,ib,nb,
     &           idim,jdim,kdim,vx,vy,vz)

C ... Write EnSight Gold per_node vector file.

      implicit none

      character(len=*)  :: name
      character(len=80) :: buffer
      integer :: i, n, ib, nb, iunit, idim, jdim, kdim
      real, dimension(*) :: vx, vy, vz

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         buffer = "EnSight Gold per_node vector file"
         write(iunit) buffer

      endif

      buffer = "part"
      write(iunit) buffer

      write(iunit) ib

      buffer = "block"
      write(iunit) buffer

      n = idim*jdim*kdim

      write(iunit) (real(vx(i),4),i=1,n)
      write(iunit) (real(vy(i),4),i=1,n)
      write(iunit) (real(vz(i),4),i=1,n)

      if(ib == nb) close(iunit)

      return 
      end subroutine WriteEnSightVector		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WriteEnSightTensor(name,iunit,ib,nb,
     &           idim,jdim,kdim,v11,v22,v33,v12,v13,v23)

C ... Write EnSight Gold per_node symmetric tensor file.

      implicit none

      character(len=*)  :: name
      character(len=80) :: buffer
      integer :: i, n, ib, nb, iunit, idim, jdim, kdim
      real, dimension(*) :: v11, v22, v33, v12, v13, v23 

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         buffer = "EnSight Gold per_node symmetric tensor file"
         write(iunit) buffer

      endif

      buffer = "part"
      write(iunit) buffer

      write(iunit) ib

      buffer = "block"
      write(iunit) buffer

      n = idim*jdim*kdim

      write(iunit) (real(v11(i),4),i=1,n)
      write(iunit) (real(v22(i),4),i=1,n)
      write(iunit) (real(v33(i),4),i=1,n)
      write(iunit) (real(v12(i),4),i=1,n)
      write(iunit) (real(v13(i),4),i=1,n)
      write(iunit) (real(v23(i),4),i=1,n)

      if(ib == nb) close(iunit)

      return 
      end subroutine WriteEnSightTensor		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WriteEnSightTensor9(name,iunit,ib,nb,
     &           idim,jdim,kdim,v11,v12,v13,v21,v22,v23,v31,v32,v33)

C ... Write EnSight Gold per_node asymmetric tensor file.

      implicit none

      character(len=*)  :: name
      character(len=80) :: buffer
      integer :: i, n, ib, nb, iunit, idim, jdim, kdim
      real, dimension(*) :: v11, v12, v13, v21, v22, v23, v31, v32, v33

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         buffer = "EnSight Gold per_node asymmetric tensor file"
         write(iunit) buffer

      endif

      buffer = "part"
      write(iunit) buffer

      write(iunit) ib

      buffer = "block"
      write(iunit) buffer

      n = idim*jdim*kdim

      write(iunit) (real(v11(i),4),i=1,n)
      write(iunit) (real(v12(i),4),i=1,n)
      write(iunit) (real(v13(i),4),i=1,n)
      write(iunit) (real(v21(i),4),i=1,n)
      write(iunit) (real(v22(i),4),i=1,n)
      write(iunit) (real(v23(i),4),i=1,n)
      write(iunit) (real(v31(i),4),i=1,n)
      write(iunit) (real(v32(i),4),i=1,n)
      write(iunit) (real(v33(i),4),i=1,n)

      if(ib == nb) close(iunit)

      return 
      end subroutine WriteEnSightTensor9		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WritePlot3dGeometry(name,iunit,nchim,ib,nb,
     &           idim,jdim,kdim,x,y,z,blank)

C ... Write Plot3d geometry file.

      implicit none

      character(len=*) :: name
      integer :: i, n, ib, nb, nchim, iunit
      integer, dimension(*) :: idim, jdim, kdim
      real, dimension(*) :: x, y, z, blank

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         write(iunit) nb
         write(iunit) (idim(i)+1,jdim(i)+1,kdim(i)+1,i=1,nb)

      endif

      n = (idim(ib)+1)*(jdim(ib)+1)*(kdim(ib)+1)

c ... Write node coordinates and optional blanking

      if(nchim > 0) then
         write(iunit) (real(X(I),4),I=1,n),(real(Y(I),4),I=1,n),
     &                (real(Z(I),4),I=1,n),(NINT(blank(i)),i=1,n)
      else
         write(iunit) (real(X(I),4),I=1,n),(real(Y(I),4),I=1,n),
     &                (real(Z(I),4),I=1,n)
      endif

      if(ib == nb) close(iunit)

      return 
      end subroutine WritePlot3dGeometry		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine WritePlot3dSolution(name,iunit,nchim,ib,nb,
     &           idim,jdim,kdim,rmach,alpha,re,t,icycle,timel,
     &           s1,s2,s3,s4,s5)

C ... Write Plot3d solution file.

      implicit none

      character(len=*) :: name
      integer :: i, n, ib, nb, nchim, iunit, icycle
      real :: rmach, alpha, re, t
      logical :: timel
      integer, dimension(*) :: idim, jdim, kdim
      real, dimension(*) :: s1, s2, s3, s4, s5

      if(ib == 1) then

         open(iunit,file=name,status='unknown',form='unformatted')

         write(iunit) nb
         write(iunit) (idim(i)+1,jdim(i)+1,kdim(i)+1,i=1,nb)

      endif

      if (timel) then   
         write(iunit) REAL(rmach,4),REAL(alpha,4),REAL(re,4),REAL(t,4)
      else
         write(iunit) REAL(rmach,4),REAL(alpha,4),REAL(re,4),
     +                REAL(icycle,4)
      endif                

      n = (idim(ib)+1)*(jdim(ib)+1)*(kdim(ib)+1)

c ... write solution

      write(iunit) (real(s1(i),4),i=1,n),(real(s2(i),4),i=1,n),
     &             (real(s3(i),4),i=1,n),(real(s4(i),4),i=1,n),
     &             (real(s5(i),4),i=1,n) 

      if(ib == nb) close(iunit)

      return 
      end subroutine WritePlot3dSolution		
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      subroutine CopyAsciiFile(SourceFile,TargetFile)

C ... This subroutine copies an ASCII file from one location to another.

      integer :: unitin, unitou
      
      character(len=80) :: Line
      character(len=*)  :: SourceFile, TargetFile

      unitin = 285
      unitou = 286

      open(unitin,file=SourceFile,status='old',form='formatted')
      open(unitou,file=TargetFile,status='unknown',form='formatted')

 1    continue
      read(unitin,'(A)',END=2) Line
      write(unitou,'(A)') Line
      goto 1

 2    continue

      close(285)
      close(286)
      
      return
      end subroutine CopyAsciiFile
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

