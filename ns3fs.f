C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface boundary conditions ---------
C --- by Farmer's approach
C ----------------------------------------------------------------------
C
      SUBROUTINE FS_FRE(RO,RM,RN,RW,E,EPS2,VIST,RK,REPS,FI,
     +     MAXSB,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     M,N1,N2,XVEL,VOL,XC,YC,ZC,A3,DTL,TE,
     +     LHULL,MMHUL,MHULL,IBOTGR,ITOPGR,KTOPGR, !for transom
     +     INWH,UROT,VROT,WROT,XCO,YCO,ZCO,ICONH,ICOGH,IHULL,
     +     PRN,NFSD,ICYCLE,ICFST,DTWMAX,DWMV,ITURB,JFIRST,
     +     IOLD,ICMAX,DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH)
     
      USE MAIN_ARRAYS, ONLY : MGRID,IG,IR,IQ,IC,JF,IHF,ICON,NPATCH,
     +     IMAX,JMAX,KMAX,WH,XXI,YXI,XETA,YETA,WFS,XHULL,YHULL,ZHULL,
     +     XVL,YVL,ZVL,P,PDIFF,U,V,W,WAVEH,WAVES,IGRID

      USE NS3CO, ONLY : FRSPRE, GROUND, FRSDEN, NSCAL, IPRINT,
     &                  GX, GY, GZ, IN, JN, KN

      IMPLICIT NONE

      INTEGER :: MAXSB,M,N1,N2,N,IG1,IR1,IQ1,IC1,IF1,JH1,JG1,ITURB,
     +           IOLD,INWH,ICYCLE,ICMAX,JFIRST,NFSD,ICFST,iii
      INTEGER :: TE
      INTEGER :: LHULL(*),MMHUL(*),MHULL(*)
      INTEGER :: IBOTGR(*),ITOPGR(*),KTOPGR(*)
      INTEGER :: ICONH(*),ICOGH(*),IHULL(*)

      REAL :: RO(*),RM(*),RN(*),RW(*),E(*),EPS2(*),VIST(*),
     +        RK(*),REPS(*),FI(MAXSB,*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +        VOL(*),UROT(*),VROT(*),WROT(*),!WH(*),
     +        A3(*),DTL(*)

      REAL :: XVEL,DTWMAX,DWMV,DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH

      REAL :: XCO(*),YCO(*),ZCO(*),XC(*),YC(*),ZC(*)
     
      CHARACTER(LEN=3) :: PRN
C
C ... THIS SUBROUTINE HANDLES FREE SURFACE PATCHES
C

C ... Initialize convergence variables

c      WHMAX   = -99999.
c      WHMIN   = 99999.
c      DWMAX   = 0.
c      SUMDWH  = 0.

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IR1     = IR(M,N)
      IQ1     = IQ(M,N)
      IC1     = IC(M,N)
      IF1     = JF(M,N)		!ASC090300
      JH1     = ICONH(N)
      JG1     = ICOGH(JH1)

      CALL FS_FREVOL(RO(IG1),RM(IG1),RN(IG1),RW(IG1),E(IG1),EPS2(IG1),
     +     VIST(IG1),P(IG1),PDIFF(IG1),RK(IR1),REPS(IR1),FI(IQ1,1),
     +     MAXSB,NSCAL,ITURB,A1X(IG1),A1Y(IG1),A1Z(IG1),
     +     A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     +     VOL(IG1),XC(IG1),YC(IG1),ZC(IG1),XVL(IG1),YVL(IG1),ZVL(IG1),
     +     UROT(IG1),VROT(IG1),WROT(IG1),
     +     ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +     IN,JN,KN,XVEL,FRSPRE,WH(IF1),XXI(IF1),YXI(IF1),
     +     XETA(IF1),YETA(IF1),WFS(IF1),IOLD,ICYCLE,ICMAX,
     +     A3(IG1),DTL,TE,XHULL(JG1),YHULL(JG1),ZHULL(JG1),LHULL(JH1),
     +     MMHUL(JH1),MHULL(JH1),IBOTGR(N),ITOPGR(N),KTOPGR(N),
     +     INWH,JFIRST,NFSD,DTWMAX,DWMV,ICFST,
     +     XCO(IG1),YCO(IG1),ZCO(IG1),U(IG1),V(IG1),
     +     W(IG1),PRN,IHULL(IG1),FRSDEN,DWMAX,FLODWH,WHMAX,WHMIN,SUMDWH)
      ENDIF

C ... Put WH array into double-precision array WAVEH

      CALL FS_SETWAVEH(WAVEH,WH(IF1),IHF(1,M),ICON(IC1),NPATCH(N),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,ICYCLE,IPRINT)

7900  CONTINUE

      CLOSE(99)
      RETURN
      END SUBROUTINE FS_FRE


      SUBROUTINE FS_VERPAT(N1,N2,M)

      USE NS3CO, ONLY : ICYCLE, IPRINT, IN, JN, KN

      USE MAIN_ARRAYS, ONLY : NPATCH,IC,IHF,IMAX,JMAX,KMAX,WAVEH,
     +   WAVES,ICON

      IMPLICIT NONE

      INTEGER :: N1,N2,M,N,IC1,III

C ... THIS SUBROUTINE SETS FREE SURFACE PATCHES INTO CELL-VERTEX FORM

      DO N = N1,N2           ! BLOCK LOOP 1

      IC1 = IC(M,N)
      CALL FS_SETWAVES(WAVEH,WAVES,IHF(1,M),ICON(IC1),NPATCH(N),
     +     IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,ICYCLE,IPRINT)

      ENDDO

      RETURN
      END SUBROUTINE FS_VERPAT

C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface boundary conditions ---------
C ---A simplified Outflow condition
C ----------------------------------------------------------------------
C
      SUBROUTINE FOU(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,MAXSB,
     +     A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     M,N1,N2,XVEL,VOL,XC,YC,ZC,WH)
     
      USE NS3CO

      USE MAIN_ARRAYS, ONLY : MGRID,IG,IR,IQ,IC,ICON,NPATCH,IMAX,
     +     JMAX,KMAX

      IMPLICIT NONE

      INTEGER :: MAXSB,M,N1,N2,N,IG1,IR1,IQ1,IC1

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +        VOL(*),WH(*),A1(*),A2(*),A3(*)

      REAL :: XVEL

      REAL ::  XC(*),YC(*),ZC(*)

C
C ... THIS SUBROUTINE HANDLES FREE SURFACE PATCHES
C

      DO 7900 N = N1,N2           ! BLOCK LOOP 1
      IF(M <= MGRID(N)) THEN
      IG1     = IG(M,N)
      IR1     = IR(M,N)
      IQ1     = IQ(M,N)
      IC1     = IC(M,N)
      CALL FOUVOL(RO(IG1),U(IG1),V(IG1),W(IG1),E(IG1),EPS2(IG1),
     +     VIST(IG1),P(IG1),PDIFF(IG1),RK(IR1),REPS(IR1),FI(IQ1,1),
     +     MAXSB,NSCAL,ITURB,A1(IG1),A2(IG1),A3(IG1),A1X(IG1),
     +     A1Y(IG1),A1Z(IG1),
     +     A2X(IG1),A2Y(IG1),A2Z(IG1),A3X(IG1),A3Y(IG1),A3Z(IG1),
     +     VOL(IG1),XC(IG1),YC(IG1),ZC(IG1),
     +     ICON(IC1),NPATCH(N),IMAX(M,N),JMAX(M,N),KMAX(M,N),
     +     IN,JN,KN,XVEL,FRSPRE,WH) 
      ENDIF
7900  CONTINUE

      RETURN
      END SUBROUTINE FOU

C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface memory allocation------------
C ----------------------------------------------------------------------
C

C ***********************************************************************
      SUBROUTINE MEMHUL(MAXHB,LEVEL,MGRID)
C DETERMINE THE MEMORY SIZE MAXHB FOR XHULL,YHULL AND ZHULL

      IMPLICIT NONE

      INTEGER MAXHB,LEVEL,NMLEV,NBH,I,J,LHULL,MHULL,NHULL

      INTEGER :: DEALLOCSTAT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDHULL,JDHULL,KDHULL
      INTEGER MGRID(*)
      LOGICAL THERE

C NUMBER OF MULTIGRID LEVELS (WE DO NOT NEED HULL FOR THE COARSER LEVELS)
      
      INQUIRE(FILE='HULL.BIN',EXIST=THERE)
      IF(.NOT.THERE) THEN
          WRITE(*,*) 'File HULL.BIN does not exist musjo. Exiting..'
          STOP
      ENDIF

      NMLEV=1

      MAXHB=1
      OPEN(171,FILE='HULL.BIN',STATUS='OLD',FORM='UNFORMATTED',ERR=3)
      READ(171)NBH
      ALLOCATE(IDHULL(NBH),JDHULL(NBH),KDHULL(NBH)) 
      READ(171)(IDHULL(I),JDHULL(I),KDHULL(I),I=1,NBH)
      CLOSE(171)
      DO I=1,NBH
C     NMLEV=MGRID(I)
      DO J=LEVEL,LEVEL+NMLEV-1
      LHULL=(IDHULL(I)-1)/(2**(J-1))
      MHULL=(JDHULL(I)-1)/(2**(J-1))
      NHULL=(KDHULL(I)-1)/(2**(J-1))
C NO GHOST CELLS OR ANY OTHER INTERMEDIATE CELLS I RECKON
      MAXHB=MAXHB + (LHULL+1)*(MHULL+1)*(NHULL+1)
      END DO
      END DO
      DEALLOCATE(IDHULL,JDHULL,KDHULL,STAT=DEALLOCSTAT)
3     RETURN
      END SUBROUTINE MEMHUL
C
C
C ----------------------------------------------------------------------
C --- Subroutine for handling free surface grid             ------------
C ----------------------------------------------------------------------
C
C***********************************************************************
      SUBROUTINE RDHULL(XHULL,YHULL,ZHULL,LHULL1,MHULL1,NHULL1,NLVL)
C READ HULL CARTESIAN COORDINATES

      IMPLICIT NONE

      INTEGER LHULL1,MHULL1,NHULL1,NLVL,ISKIP,IPAS,LHULL,MHULL,NHULL,
     +     IND,JND,KND,I,J,K,IPASS

      REAL XHULL(*)
      REAL YHULL(*)
      REAL ZHULL(*)
      REAL DUMMY

      ISKIP=NLVL
      IPASS=ISKIP-1
      LHULL=LHULL1/ISKIP
      MHULL=MHULL1/ISKIP
      NHULL=NHULL1/ISKIP
      READ(171)
     $ XHULL(1),((DUMMY,IND=1,IPASS),XHULL(1+I),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ XHULL(1+J*(LHULL+1)),((DUMMY,IND=1,IPASS),
     $ XHULL(1+I+J*(LHULL+1)),I=1,LHULL),J=1,MHULL),     
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ DUMMY,((DUMMY,IND=1,IPASS),
     $ DUMMY,I=1,LHULL),J=1,MHULL),KND=1,IPASS),     
     $ XHULL(1+K*(LHULL+1)*(MHULL+1)),
     $ ((DUMMY,IND=1,IPASS),XHULL(1+I+K*(LHULL+1)*(MHULL+1)),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ XHULL(1+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),((DUMMY,IND=1,IPASS),
     $ XHULL(1+I+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),
     $ I=1,LHULL),J=1,MHULL),K=1,NHULL),     
     $ YHULL(1),((DUMMY,IND=1,IPASS),YHULL(1+I),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ YHULL(1+J*(LHULL+1)),((DUMMY,IND=1,IPASS),
     $ YHULL(1+I+J*(LHULL+1)),I=1,LHULL),J=1,MHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ DUMMY,((DUMMY,IND=1,IPASS),
     $ DUMMY,I=1,LHULL),J=1,MHULL),KND=1,IPASS),
     $ YHULL(1+K*(LHULL+1)*(MHULL+1)),
     $ ((DUMMY,IND=1,IPASS),YHULL(1+I+K*(LHULL+1)*(MHULL+1)),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ YHULL(1+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),((DUMMY,IND=1,IPASS),
     $ YHULL(1+I+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),
     $ I=1,LHULL),J=1,MHULL),K=1,NHULL),
     $ ZHULL(1),((DUMMY,IND=1,IPASS),ZHULL(1+I),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ ZHULL(1+J*(LHULL+1)),((DUMMY,IND=1,IPASS),
     $ ZHULL(1+I+J*(LHULL+1)),I=1,LHULL),J=1,MHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ DUMMY,((DUMMY,IND=1,IPASS),
     $ DUMMY,I=1,LHULL),J=1,MHULL),KND=1,IPASS),
     $ ZHULL(1+K*(LHULL+1)*(MHULL+1)),
     $ ((DUMMY,IND=1,IPASS),ZHULL(1+I+K*(LHULL+1)*(MHULL+1)),I=1,LHULL),
     $ ((DUMMY,((DUMMY,IND=1,IPASS),DUMMY,I=1,LHULL),JND=1,IPASS),
     $ ZHULL(1+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),((DUMMY,IND=1,IPASS),
     $ ZHULL(1+I+J*(LHULL+1)+K*(LHULL+1)*(MHULL+1)),
     $ I=1,LHULL),J=1,MHULL),K=1,NHULL)

      WRITE(45,'(4X,A,A,3I4)')'HULL.BIN succesfully read. LHULL, ',
     $  'MHULL, NHULL =',LHULL,MHULL,NHULL

      RETURN
      END SUBROUTINE RDHULL
      
C ######################################################################
      SUBROUTINE INP_HULL
C ######################################################################

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : MAXB,NBCS,IPRO,MAXFS

      USE MAIN_ARRAYS, ONLY : IMAX,JMAX,KMAX,LHULL,MMHUL,MHULL,XHULL,
     $ YHULL,ZHULL,ICONH,IBTGR,JBTGR,KBTGR,ITOPGR,JTOPGR,KTOPGR,IBOTGR,
     $ JBOTGR,KBOTGR,ICOGH,ICON,NPROCE,ITPGR,JTPGR,KTPGR,ICNH,IGRID

      USE NS3CO

      IMPLICIT NONE

      INTEGER NBH,L,IFSBCX,IERR,NLVL,I,NB,IEDL,IED2,IHU,KOKO

      LOGICAL READER


C**********   Reading HULL.BIN   *****************************************
C ASC030299

C       IFSBC=0
c
c       DO L = 1,NBCS     !NPATCH(1)
c        IF(ICON(1+(L-1)*25) == 13) THEN
c           IF(IFSBC /= 2)THEN
c           IFSBC=1
c           ELSE IF(IFSBC /= 0)THEN
c           WRITE(*,*)'STOP: FRE and LVS boundary conditions cannot be '
c     *     ,'used at the same time'
c           ENDIF
c        ELSE IF(ICON(1+(L-1)*25) == 14) THEN
c           IF(IFSBC /= 1)THEN
c           IFSBC=2
c           ELSE IF(IFSBC /= 0)THEN
c           WRITE(*,*)'STOP: FRE and LVS boundary conditions cannot be '
c     *     ,'used at the same time'
c           ENDIF
c        ENDIF
c       ENDDO

       READER = .FALSE.
       
       DO I = 1,NBLOCK
         IF(IGRID(I,2) == 9) READER = .TRUE. ! taa on nyt
       ENDDO                                   ! deformaatiota vai?

        IF(PARALLEL) THEN
           CALL MPI_ALLREDUCE(IFSBC,IFSBCX,1,MPI_INTEGER,MPI_MAX,
     &     MPI_COMM_WORLD,IERR)
           IFSBC = IFSBCX
        ENDIF

C ASC030299

      IF(READER) THEN

      NLVL=2**(LEVEL-1)
                    
      OPEN(171,FILE='HULL.BIN',FORM='UNFORMATTED')

      READ(171)NBH
       
      READ(171)(LHULL(I),MMHUL(I),MHULL(I),I=1,NBH)
       
         IF(NBH > NBLOCG)THEN
          WRITE(*,*)'Not allowed more blocks in HULL.BIN than NBLOCG'
          STOP
         ENDIF

C FOR 2D HULL ALLOW ALL IJ- AND IK- AND JK-SURFACE HULLS
C one of these should be = 1

      DO I=1,NBH
      IF(LHULL(I) == 1) THEN
        LHULL(I)=MHULL(I)
        MHULL(I)=MMHUL(I)
      ELSE IF(MHULL(I) == 1) THEN
        MHULL(I)=MMHUL(I)
        MMHUL(I)=1
      ELSE IF(MMHUL(I) == 1) THEN
      ELSE
      ENDIF
      ENDDO

      IEDL = 0
      KOKO = 0

      DO I = 1,NBH
      KOKO = KOKO + LHULL(I)*MMHUL(I)*MHULL(I) + 100 ! Mersu
      ENDDO

      DEALLOCATE (XHULL,YHULL,ZHULL)
      ALLOCATE(XHULL(KOKO),YHULL(KOKO),ZHULL(KOKO))
      WRITE(45,*)
      WRITE(45,*)' In INP_HULL XHULL; YHULL, ZHULL reallocated',
     $           ' Size =',KOKO

      IF(IPRO == 1) THEN

      DO I=1,NBH
       CALL RDHULL(XHULL(IEDL+1),YHULL(IEDL+1),ZHULL(IEDL+1),
     $  LHULL(I)-1,MMHUL(I)-1,MHULL(I)-1,NLVL)
        IEDL = IEDL + ((LHULL(I)-1)/NLVL+1)*
     $  ((MHULL(I)-1)/NLVL+1)*((MMHUL(I)-1)/NLVL+1)
      
      ENDDO
C1
      ENDIF ! IPRO == 1
      
   
      IF(PARALLEL) THEN
         CALL MPI_BCAST(XHULL,KOKO,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YHULL,KOKO,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZHULL,KOKO,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      ENDIF

       
      IED2= 0

      DO I=1,NBH

      IF(MOD(LHULL(I)-1,NLVL) /= 0..OR.MOD(MHULL(I)-1,NLVL) /= 0.)THEN
      WRITE(*,*)
     *'WARNING HULL.BIN IS NOT SUFFICIENTLY EVEN FOR THIS GRID LEVEL'
      STOP
      ENDIF

      LHULL(I)=(LHULL(I)-1)/NLVL
      MHULL(I)=(MHULL(I)-1)/NLVL
      MMHUL(I)=(MMHUL(I)-1)/NLVL
C2
      CLOSE(171)

      DO IHU=1+IED2,IED2+(LHULL(I)+1)*(MMHUL(I)+1)*(MHULL(I)+1)
        XHULL(IHU)= XHULL(IHU)*GRILEN
        YHULL(IHU)= YHULL(IHU)*GRILEN
        ZHULL(IHU)= ZHULL(IHU)*GRILEN
      ENDDO

      ICOGH(I)= IED2+1
      IED2 = IED2 + (LHULL(I)+1)*(MMHUL(I)+1)*(MHULL(I)+1)
     
      ENDDO

c      WRITE(*,*)'GRILEN=', GRILEN

      ICOGH(NBLOCG+1)=1

      DO I=1,NBLOCK
      IBOTGR(I)=IBTGR(NPROCE(1+I,IPRO))
      JBOTGR(I)=JBTGR(NPROCE(1+I,IPRO))
      KBOTGR(I)=KBTGR(NPROCE(1+I,IPRO))
      ITOPGR(I)=ITPGR(NPROCE(1+I,IPRO))
      JTOPGR(I)=JTPGR(NPROCE(1+I,IPRO))
      KTOPGR(I)=KTPGR(NPROCE(1+I,IPRO))
      ICONH(I) =ICNH(NPROCE(1+I,IPRO))

c      READ(199,*)
c      READ(199,*)
c      READ(199,*)
c      READ(199,*) INDXTS
      REWIND 199

      ENDDO

      ELSE ! IF READER

       DO I=1,NBLOCK
       ICOGH(I)=1
       ICONH(I)=1
       ENDDO
      
      ENDIF ! (IFSBC == 1)
c         do i=1,nbh
c         write(6,*) 'no',IPRO,i,iconh(i),icogh(i),nblock
c         enddo

      END SUBROUTINE INP_HULL
