C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE CHIMER(M)

C ... Calling routine for Chimera interpolation

      USE MPI

      USE CHARACTERS

      USE INTEGERS, ONLY : IPRO, MAXB
      
      USE MAIN_ARRAYS

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: N,M,I,NC3,NCG,NC,N3,NG,IG1,IR1,IC1,IG2,IR2,IC2,
     &           IMAXW,JMAXW,KMAXW,NPATCW,NCHIMW,IG1W,NGL,III,NOFDONORS

      REAL :: RKSIV

C ... Allocatable arrays for data to be received from other processes
C ... in case of parallel simulation.

      REAL, ALLOCATABLE :: BLANKW(:), VOLW(:), DISTWW(:),RKSIW(:),
     &                     WGH1W(:), WGH2W(:), WGH3W(:), WGH4W(:)

      REAL, ALLOCATABLE :: XCW(:), YCW(:), ZCW(:)

      INTEGER, ALLOCATABLE :: ICONW(:), LOCDISW(:), INTPW(:),
     &                        II1W(:), II2W(:), II3W(:), II4W(:)

      INTEGER, DIMENSION(5) :: MSG
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: st
      INTEGER :: pid, tag, rc
      INTEGER :: S, LI, IPHASE

      LOGICAL :: NOBBOXHIT
             
C ... Chi value in the fortified algorithm.
      RKSIV = 1.0E+10

C ... Initialize the fortify vectors

      DO I=1,MAXB

         ROFOR(I) = 0.0
         RMFOR(I) = 0.0
         RNFOR(I) = 0.0
         RWFOR(I) = 0.0
         EFOR(I)  = 0.0
         PDFOR(I) = 0.0
         RKFOR(I) = 0.0
         REFOR(I) = 0.0
         RKSI(I)  = 0.0

         IF(TRANSL) THEN
            TRM(I)%GFOR   = 0.0
            TRM(I)%RETFOR = 0.0
         ENDIF

         IF(MULPHL) THEN
            DO IPHASE = 1,NPHASES
               MPFOR(I)%ALFAFOR(IPHASE)  = 0.0
               MPFOR(I)%XFOR(IPHASE)     = 0.0
               MPFOR(I)%DTEMPFOR(IPHASE) = 0.0
            ENDDO
         ENDIF

         IF(.NOT.TRUE_DISTL) BLANK(I) = 1.0  

      ENDDO


C ... Initialize the Chimera donor cell information

      DO I=1,MAXB
         II1(I)  = 0
         II2(I)  = 0
         II3(I)  = 0
         II4(I)  = 0
         INTP(I) = 0
         WGH1(I) = 0.0
         WGH2(I) = 0.0
         WGH3(I) = 0.0
         WGH4(I) = 0.0
      ENDDO


      DO NC3=1,NBLOCG                 ! Chimera block loop

         NCG = NCHOR(NC3)             ! global block number
         NC  = NLOCAL(NCG)            ! local block number
         IF(NCHIMT(NCG) < -1) CYCLE   ! Skip Chimera interpolation

      DO N3=1,NBLOCG                  ! Background block loop

         NG  = NCHOR(N3)              ! global block number
         N   = NLOCAL(NG)             ! local block number

         IF(NG == NCG) CYCLE

         IF(NOBBOXHIT(NG,NCG)) CYCLE  ! Block bounding boxes don't hit

         IF(NEARBLOCK((NG-1)*NBLOCG + NCG) < 1 .AND.
     &   (NCHIMT(NG) /= -1 .AND. NCHIMT(NCG) /= -1)) CYCLE

      IF((NCHIMT(NG) >= -1 .AND. NCHIMT(NG) < NCHIMT(NCG)) .OR.
     &   (NCHIMT(NG) == -1 .AND. NCHIMT(NCG) == -1)) THEN
 
         IF(N /= -1) THEN
            IG1 = IG(M,N) 
            IR1 = IR(M,N)
            IC1 = IC(M,N) 
         ENDIF

         IF(NC /= -1) THEN
            IG2 = IG(M,NC) 
            IR2 = IR(M,NC)
            IC2 = IC(M,NC) 
         ENDIF


      IF(N /= -1 .AND. NC /= -1) THEN 


C ... Background block and Chimera block are both local in a parallel run
C ... or this is a single processor run

      CALL DONORS(ITURB,LEVEL,NCHIMT,TRUE_DISTL,
     &   NG,IMAX(M,N),JMAX(M,N),KMAX(M,N),IN,JN,KN,
     &   XC(IG1),YC(IG1),ZC(IG1),BLANK(IG1),DISTW(IG1),LOCDIS(IG1),
     &   VOL(IG1),ICON(IC1),NPATCH(N),M,
     &   NCG,IMAX(M,NC),JMAX(M,NC),KMAX(M,NC),
     &   XC(IG2),YC(IG2),ZC(IG2),BLANK(IG2),DISTW(IG2),LOCDIS(IG2),
     &   VOL(IG2),ICON(IC2),NPATCH(NC),
     &   RKSI(IG1),RKSI(IG2),RKSIV,
     &   II1(IG1),II2(IG1),II3(IG1),II4(IG1),INTP(IG1),
     &   II1(IG2),II2(IG2),II3(IG2),II4(IG2),INTP(IG2),
     &   WGH1(IG1),WGH2(IG1),WGH3(IG1),WGH4(IG1),IG1,
     &   WGH1(IG2),WGH2(IG2),WGH3(IG2),WGH4(IG2),IG2,NOFDONORS)

************************************************************************
      IF(NOFDONORS == 0) NEARBLOCK((NG-1)*NBLOCG + NCG) = 0
************************************************************************

C ... If only the Chimera block is local, then receive the background 
C ... block arrays and call the Chimera routines. If only the background 
C ... block is local, then send the background arrays to the Chimera 
C ... block process.

      ELSEIF(NC /= -1) THEN  ! Chimera block is local

      pid    = NPNUM(NG) - 1 
      tag    = NG
      CALL MPI_RECV(MSG,5,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      IMAXW  = MSG(1)
      JMAXW  = MSG(2)
      KMAXW  = MSG(3)
      NPATCW = MSG(4)
      NCHIMW = MSG(5)

      S = (IMAXW + 2*IN)*(JMAXW + 2*JN)*(KMAXW + 2*KN)

      ALLOCATE(XCW(S),YCW(S),ZCW(S),BLANKW(S),DISTWW(S),LOCDISW(S))
      ALLOCATE(RKSIW(S),VOLW(S))
      ALLOCATE(ICONW(IC9*NPATCW))
      ALLOCATE(II1W(S),II2W(S),II3W(S),II4W(S),INTPW(S))
      ALLOCATE(WGH1W(S),WGH2W(S),WGH3W(S),WGH4W(S))

      CALL MPI_RECV(XCW,   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(YCW,   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(ZCW,   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(BLANKW,S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(VOLW,  S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(RKSIW, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)

      CALL MPI_RECV(DISTWW, S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(LOCDISW,S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      CALL MPI_RECV(II1W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II2W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II3W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II4W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(INTPW,S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH1W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH2W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH3W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH4W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(IG1W, 1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      LI = IC9*NPATCW
      CALL MPI_RECV(ICONW,LI,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      CALL DONORS(ITURB,LEVEL,NCHIMT,TRUE_DISTL,
     &   NG,IMAXW,JMAXW,KMAXW,IN,JN,KN,
     &   XCW,YCW,ZCW,BLANKW,DISTWW,LOCDISW,
     &   VOLW,ICONW,NPATCW,M,
     &   NCG,IMAX(M,NC),JMAX(M,NC),KMAX(M,NC),
     &   XC(IG2),YC(IG2),ZC(IG2),BLANK(IG2),DISTW(IG2),LOCDIS(IG2),
     &   VOL(IG2),ICON(IC2),NPATCH(NC),
     &   RKSIW,RKSI(IG2),RKSIV,
     &   II1W,II2W,II3W,II4W,INTPW,
     &   II1(IG2),II2(IG2),II3(IG2),II4(IG2),INTP(IG2),
     &   WGH1W,WGH2W,WGH3W,WGH4W,IG1W,
     &   WGH1(IG2),WGH2(IG2),WGH3(IG2),WGH4(IG2),IG2,NOFDONORS)

************************************************************************
      IF(NOFDONORS == 0) NEARBLOCK((NG-1)*NBLOCG + NCG) = 0
      CALL MPI_SEND(NOFDONORS,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
************************************************************************
       
C ... Send the fortify arrays and the blanking information back to the 
C ... background block process and deallocate the working arrays.

      IF(NOFDONORS > 0) THEN

         CALL MPI_SEND(RKSIW, S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(BLANKW,S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(II1W,  S,MPI_INTEGER, pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(II2W,  S,MPI_INTEGER, pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(II3W,  S,MPI_INTEGER, pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(II4W,  S,MPI_INTEGER, pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(INTPW, S,MPI_INTEGER, pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(WGH1W, S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(WGH2W, S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(WGH3W, S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(WGH4W, S,MPI_REAL8,    pid,tag,MPI_COMM_WORLD,rc)

      ENDIF


      DEALLOCATE(XCW,YCW,ZCW,BLANKW,DISTWW,LOCDISW)
      DEALLOCATE(RKSIW,VOLW)
      DEALLOCATE(ICONW)
      DEALLOCATE(II1W,II2W,II3W,II4W,INTPW)
      DEALLOCATE(WGH1W,WGH2W,WGH3W,WGH4W)


      ELSEIF(N /= -1) THEN  ! Background block is local
      
      pid    = NPNUM(NCG) - 1 
      tag    = NG
      MSG(1) = IMAX(M,N)
      MSG(2) = JMAX(M,N)
      MSG(3) = KMAX(M,N)
      MSG(4) = NPATCH(N)
      MSG(5) = NCHIMT(NG)
      CALL MPI_SEND(MSG,5,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)

      S = (IMAX(M,N) + 2*IN)*(JMAX(M,N) + 2*JN)*(KMAX(M,N) + 2*KN)

      CALL MPI_SEND(XC(IG1),   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,rc)
      CALL MPI_SEND(YC(IG1),   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,rc)
      CALL MPI_SEND(ZC(IG1),   S,MPI_REAL8,pid,tag,
     &              MPI_COMM_WORLD,rc)
      CALL MPI_SEND(BLANK(IG1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(VOL(IG1),  S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(RKSI(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(DISTW(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(LOCDIS(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II1(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II2(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II3(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II4(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(INTP(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH1(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH2(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH3(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH4(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(IG1,      1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)

      LI = IC9*NPATCH(N)
      CALL MPI_SEND(ICON(IC1),LI,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)


C ... Receive the updated fortify arrays and the blanking information
C ... from the Chimera block process. 

************************************************************************
      CALL MPI_RECV(NOFDONORS,1,MPI_INTEGER,pid,tag,
     &              MPI_COMM_WORLD,st,rc)
      IF(NOFDONORS == 0) NEARBLOCK((NG-1)*NBLOCG + NCG) = 0
************************************************************************
 
      IF(NOFDONORS > 0) THEN

      CALL MPI_RECV(RKSI(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(BLANK(IG1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II1(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II2(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II3(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II4(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(INTP(IG1),S,MPI_INTEGER,pid,tag,
     &                                            MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH1(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH2(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH3(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH4(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)

      ENDIF

      ENDIF  ! Chimera block or background block not local

      ENDIF  ! Chimera order


      ENDDO  ! Background block loop

      ENDDO  ! Chimera block loop


C ... Try to fix isolated cells.

      DO N = 1,NBLOCK
         IG1 = IG(1,N)
         NGL = NPROCE(1+N,IPRO)  ! Global block number
         CALL FIXKSI(RKSI(IG1),WGH1(IG1),
     &   IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,RKSIV,NGL,IPRO)
      ENDDO

      RETURN
      END SUBROUTINE CHIMER
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE CHIMUP(M)

C ... This subroutine updates the flow field connections 
C ... in the overlapping grids using the precomputed donor 
C ... cell information (subroutine CHIMER). 

      USE MPI

      USE CHARACTERS

      USE INTEGERS,    ONLY : IPRO, MAXB

      USE MAIN_ARRAYS, ONLY : NTOT,IG,ROFOR,RMFOR,RNFOR,RWFOR,EFOR,
     &    RKFOR,REFOR,TRM,MPFOR,VAR,PRO,NEARBLOCK,
     &    RO,RM,RN,RW,E,PDIFF,RK,REPS,PDFOR,TEMP,
     &    U,V,W,INTP,RKSI,WGH1,WGH2,WGH3,WGH4,II1,II2,II3,II4,NPROCE,XC,
     &    YC,ZC,BLANK,OMEGA,OMEGAX,OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,
     &    IMAX,JMAX,KMAX,IT,IL,IK,IR,NCHOR,NLOCAL,IC,NCHIMT,
     &    NPATCH,NPNUM,VIST,EPS2,BLANK2,BLKS,P,F1RM,F1RN,ZZZ

      USE NS3CO

      USE TYPE_ARRAYS , ONLY : BLOCKS, NPHASES

      IMPLICIT NONE

      INTEGER :: M,NN,I,J,IERRORCODE,IERR,NGL,IG1,IR1,IG2,IR2
      INTEGER :: N,NG,N3,NC,NCG,NC3,IG1W
      INTEGER :: II1I, II2I, II3I, II4I

C ... Allocatable arrays for data to be received from other processes
C ... in case of parallel simulation.

      REAL, ALLOCATABLE :: PDIFFW(:),TEMPW(:),
     &      ROW(:),RMW(:),RNW(:),RWW(:),EW(:),RKW(:),REPSW(:),
     &      GW(:),RETW(:),VISTW(:),EPS2W(:),
     &      WGH1W(:),WGH2W(:),WGH3W(:),WGH4W(:),
     &      RKSIW(:), TRAIN(:), ALFAW(:,:), XW(:,:), MPTRAIN(:)
      REAL, ALLOCATABLE :: DTEMPW(:,:), MPTRAINDP(:)
      REAL, ALLOCATABLE :: BUFFERDP(:)
 
      INTEGER, ALLOCATABLE :: II1W(:),II2W(:),II3W(:),II4W(:),INTPW(:)

      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: st
      INTEGER :: pid, tag, rc, S
      INTEGER :: INTPI, INTPB, NSEATS, NS, IPHASE
      INTEGER :: L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12
      INTEGER :: LMP, IMP

      LOGICAL :: MPIB, MPIC, NOBBOXHIT

           
      DO NC3 = 1,NBLOCG               ! Chimera block loop

         NCG = NCHOR(NC3)             ! global block number
         NC  = NLOCAL(NCG)            ! local block number
         IF(NCHIMT(NCG) < -1) CYCLE   ! Skip Chimera interpolation

      DO N3 = 1,NBLOCG                ! Background block loop

         NG  = NCHOR(N3)              ! global block number
         N   = NLOCAL(NG)             ! local block number

      IF(NG == NCG) CYCLE

      IF(NOBBOXHIT(NG,NCG)) CYCLE  ! Block bounding boxes don't hit

      IF(NEARBLOCK((NG-1)*NBLOCG + NCG) < 1 .AND.
     &(NCHIMT(NG) /= -1 .AND. NCHIMT(NCG) /= -1)) CYCLE

      IF((NCHIMT(NG) >= -1 .AND. NCHIMT(NG) < NCHIMT(NCG)) .OR.
     &   (NCHIMT(NG) == -1 .AND. NCHIMT(NCG) == -1)) THEN
 
         IF(N /= -1) THEN
            IG1 = IG(M,N) 
            IR1 = IR(M,N)
         ENDIF

         IF(NC /= -1) THEN
            IG2 = IG(M,NC) 
            IR2 = IR(M,NC)
         ENDIF


      IF(N /= -1 .AND. NC /= -1) THEN   

C ... Background block and Chimera block are both local in a parallel run
C ... or this is a single processor run

C ... First interpolate from the background block to the Chimera block 

      DO NN = 1,NTOT(M,NC)

         I = NN - 1 + IG2      

         INTPI = MOD(INTP(I),10)

         IF(INTPI == 2 .OR. RKSI(I) > 5.0E9) THEN

         INTPB = (INTP(I)-INTPI)/10

         IF(INTPB /= NG) CYCLE

         II1I = II1(I)
         II2I = II2(I)
         II3I = II3(I)
         II4I = II4(I)

         ROFOR(I) = WGH1(I)*RO(II1I)+WGH2(I)*RO(II2I)
     &            + WGH3(I)*RO(II3I)+WGH4(I)*RO(II4I) 

         RMFOR(I) = WGH1(I)*RM(II1I)+WGH2(I)*RM(II2I)
     &            + WGH3(I)*RM(II3I)+WGH4(I)*RM(II4I) 

         RNFOR(I) = WGH1(I)*RN(II1I)+WGH2(I)*RN(II2I)
     &            + WGH3(I)*RN(II3I)+WGH4(I)*RN(II4I) 

         RWFOR(I) = WGH1(I)*RW(II1I)+WGH2(I)*RW(II2I)
     &            + WGH3(I)*RW(II3I)+WGH4(I)*RW(II4I) 

          EFOR(I) = WGH1(I)*E(II1I)+WGH2(I)*E(II2I)
     &            + WGH3(I)*E(II3I)+WGH4(I)*E(II4I) 
       
         PDFOR(I) = WGH1(I)*PDIFF(II1I)+WGH2(I)*PDIFF(II2I)
     &            + WGH3(I)*PDIFF(II3I)+WGH4(I)*PDIFF(II4I)

          VIST(I) = WGH1(I)*VIST(II1I)+WGH2(I)*VIST(II2I)
     &            + WGH3(I)*VIST(II3I)+WGH4(I)*VIST(II4I) 
       
          EPS2(I) = WGH1(I)*EPS2(II1I)+WGH2(I)*EPS2(II2I)
     &            + WGH3(I)*EPS2(II3I)+WGH4(I)*EPS2(II4I)

         IF(ITURB >= 3 .AND. ITURB /=8) THEN

            RKFOR(I) = WGH1(I)*RK(II1I)+WGH2(I)*RK(II2I)
     &               + WGH3(I)*RK(II3I)+WGH4(I)*RK(II4I) 

            REFOR(I) = WGH1(I)*REPS(II1I)+WGH2(I)*REPS(II2I)
     &               + WGH3(I)*REPS(II3I)+WGH4(I)*REPS(II4I)

         ENDIF

         IF(TRANSL) THEN

            TRM(I)%GFOR   = WGH1(I)*TRM(II1I)%G+WGH2(I)*TRM(II2I)%G
     &                    + WGH3(I)*TRM(II3I)%G+WGH4(I)*TRM(II4I)%G 

            TRM(I)%RETFOR = WGH1(I)*TRM(II1I)%RET+WGH2(I)*TRM(II2I)%RET
     &                    + WGH3(I)*TRM(II3I)%RET+WGH4(I)*TRM(II4I)%RET

         ENDIF

         IF(MULPHL) THEN

            DO IPHASE = 1,NPHASES        

            MPFOR(I)%ALFAFOR(IPHASE)  = WGH1(I)*VAR(II1I)%ALFA(IPHASE)
     &                                + WGH2(I)*VAR(II2I)%ALFA(IPHASE)
     &                                + WGH3(I)*VAR(II3I)%ALFA(IPHASE)
     &                                + WGH4(I)*VAR(II4I)%ALFA(IPHASE)

            MPFOR(I)%XFOR(IPHASE)     = WGH1(I)*VAR(II1I)%X(IPHASE)
     &                                + WGH2(I)*VAR(II2I)%X(IPHASE)
     &                                + WGH3(I)*VAR(II3I)%X(IPHASE)
     &                                + WGH4(I)*VAR(II4I)%X(IPHASE)

            MPFOR(I)%DTEMPFOR(IPHASE) = WGH1(I)*PRO(II1I)%DTEMP(IPHASE)
     &                                + WGH2(I)*PRO(II2I)%DTEMP(IPHASE)
     &                                + WGH3(I)*PRO(II3I)%DTEMP(IPHASE)
     &                                + WGH4(I)*PRO(II4I)%DTEMP(IPHASE)

            ENDDO

         ENDIF

            IF(INTPI == 2) THEN ! This is a Chimera ghost cell

               RO(I) = ROFOR(I)
               RM(I) = RMFOR(I)
               RN(I) = RNFOR(I)
               RW(I) = RWFOR(I)
                E(I) = EFOR(I)
            PDIFF(I) = PDFOR(I)
      
             TEMP(I) = WGH1(I)*TEMP(II1I)+WGH2(I)*TEMP(II2I)
     &               + WGH3(I)*TEMP(II3I)+WGH4(I)*TEMP(II4I) 

                U(I) = RM(I)/RO(I)
                V(I) = RN(I)/RO(I)
                W(I) = RW(I)/RO(I)
      
           IF(ITURB >= 3 .AND. ITURB/=8) THEN
               RK(I) = RKFOR(I)
             REPS(I) = REFOR(I)
           ENDIF

           IF(TRANSL) THEN
              TRM(I)%G   = TRM(I)%GFOR
              TRM(I)%RET = TRM(I)%RETFOR
           ENDIF

           IF(MULPHL) THEN

              DO IPHASE = 1,NPHASES        

                 VAR(I)%ALFA(IPHASE)  = MPFOR(I)%ALFAFOR(IPHASE)
                 VAR(I)%X(IPHASE)     = MPFOR(I)%XFOR(IPHASE)
                 PRO(I)%DTEMP(IPHASE) = MPFOR(I)%DTEMPFOR(IPHASE)

              ENDDO

           ENDIF

           ENDIF

         ENDIF

      ENDDO

C ... and then interpolate from the Chimera block to the background block.

      DO NN = 1,NTOT(M,N)

         I = IG1 + NN - 1       

         IF(RKSI(I) > 5.0E9) THEN

         INTPI = MOD(INTP(I),10)
         INTPB = (INTP(I)-INTPI)/10

         IF(INTPB /= NCG) CYCLE

         II1I = II1(I)
         II2I = II2(I)
         II3I = II3(I)
         II4I = II4(I)

         ROFOR(I) = WGH1(I)*RO(II1I)+WGH2(I)*RO(II2I)
     &            + WGH3(I)*RO(II3I)+WGH4(I)*RO(II4I) 

         RMFOR(I) = WGH1(I)*RM(II1I)+WGH2(I)*RM(II2I)
     &            + WGH3(I)*RM(II3I)+WGH4(I)*RM(II4I) 

         RNFOR(I) = WGH1(I)*RN(II1I)+WGH2(I)*RN(II2I)
     &            + WGH3(I)*RN(II3I)+WGH4(I)*RN(II4I) 

         RWFOR(I) = WGH1(I)*RW(II1I)+WGH2(I)*RW(II2I)
     &            + WGH3(I)*RW(II3I)+WGH4(I)*RW(II4I) 

          EFOR(I) = WGH1(I)*E(II1I)+WGH2(I)*E(II2I)
     &            + WGH3(I)*E(II3I)+WGH4(I)*E(II4I) 
       
         PDFOR(I) = WGH1(I)*PDIFF(II1I)+WGH2(I)*PDIFF(II2I)
     &            + WGH3(I)*PDIFF(II3I)+WGH4(I)*PDIFF(II4I)

          VIST(I) = WGH1(I)*VIST(II1I)+WGH2(I)*VIST(II2I)
     &            + WGH3(I)*VIST(II3I)+WGH4(I)*VIST(II4I) 
       
          EPS2(I) = WGH1(I)*EPS2(II1I)+WGH2(I)*EPS2(II2I)
     &            + WGH3(I)*EPS2(II3I)+WGH4(I)*EPS2(II4I)

         IF(ITURB >= 3 .AND. ITURB /=8) THEN

            RKFOR(I) = WGH1(I)*RK(II1I)+WGH2(I)*RK(II2I)
     &               + WGH3(I)*RK(II3I)+WGH4(I)*RK(II4I) 

            REFOR(I) = WGH1(I)*REPS(II1I)+WGH2(I)*REPS(II2I)
     &               + WGH3(I)*REPS(II3I)+WGH4(I)*REPS(II4I)

         ENDIF

         IF(TRANSL) THEN

            TRM(I)%GFOR   = WGH1(I)*TRM(II1I)%G+WGH2(I)*TRM(II2I)%G
     &                    + WGH3(I)*TRM(II3I)%G+WGH4(I)*TRM(II4I)%G 

            TRM(I)%RETFOR = WGH1(I)*TRM(II1I)%RET+WGH2(I)*TRM(II2I)%RET
     &                    + WGH3(I)*TRM(II3I)%RET+WGH4(I)*TRM(II4I)%RET

         ENDIF

         IF(MULPHL) THEN

            DO IPHASE = 1,NPHASES        

            MPFOR(I)%ALFAFOR(IPHASE)  = WGH1(I)*VAR(II1I)%ALFA(IPHASE)
     &                                + WGH2(I)*VAR(II2I)%ALFA(IPHASE)
     &                                + WGH3(I)*VAR(II3I)%ALFA(IPHASE)
     &                                + WGH4(I)*VAR(II4I)%ALFA(IPHASE)

            MPFOR(I)%XFOR(IPHASE)     = WGH1(I)*VAR(II1I)%X(IPHASE)
     &                                + WGH2(I)*VAR(II2I)%X(IPHASE)
     &                                + WGH3(I)*VAR(II3I)%X(IPHASE)
     &                                + WGH4(I)*VAR(II4I)%X(IPHASE)

            MPFOR(I)%DTEMPFOR(IPHASE) = WGH1(I)*PRO(II1I)%DTEMP(IPHASE)
     &                                + WGH2(I)*PRO(II2I)%DTEMP(IPHASE)
     &                                + WGH3(I)*PRO(II3I)%DTEMP(IPHASE)
     &                                + WGH4(I)*PRO(II4I)%DTEMP(IPHASE)

            ENDDO

         ENDIF

         ENDIF

      ENDDO


      ELSEIF(NC /= -1) THEN  ! Only the Chimera block is local

C ... If only the Chimera block is local, then receive the background 
C ... block arrays and execute the interpolation routines. If only the 
C ... background block is local, then send the background arrays to the 
C ... Chimera block process.
 
      pid    = NPNUM(NG) - 1 
      tag    = NG

      MPIC = .FALSE.

      NSEATS = 0
      DO NN = 1,NTOT(M,NC)
         I  = IG2 + NN - 1
         INTPI = MOD(INTP(I),10)
         IF(INTPI == 2 .OR. RKSI(I) > 5.0E9) THEN  
            INTPB = (INTP(I)-INTPI)/10
            IF(INTPB == NG) NSEATS = NSEATS + 1
         ENDIF
      ENDDO

      MPIC = NSEATS > 0

      CALL MPI_SEND(MPIC,1,MPI_LOGICAL,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_RECV(MPIB,1,MPI_LOGICAL,pid,tag,MPI_COMM_WORLD,st,rc)


      IF(MPIC) THEN  ! Chimera block wants data from background

      CALL MPI_RECV(S,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      ALLOCATE(ROW(S),RMW(S),RNW(S),RWW(S),EW(S))
      ALLOCATE(PDIFFW(S),TEMPW(S),VISTW(S),EPS2W(S))
   
      IF(ITURB >= 3 .AND. ITURB /=8) THEN
         ALLOCATE(RKW(S),REPSW(S))
      ENDIF
   
      IF(TRANSL) THEN
         ALLOCATE(GW(S),RETW(S))
      ENDIF
   
      IF(MULPHL) THEN
         ALLOCATE(ALFAW(NPHASES,S),XW(NPHASES,S),DTEMPW(NPHASES,S))
         ALLOCATE(BUFFERDP(S))
      ENDIF

      CALL MPI_RECV(ROW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(RMW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(RNW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(RWW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(EW,    S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(PDIFFW,S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(TEMPW, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(VISTW, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(EPS2W, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(IG1W,  1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      IF(ITURB >= 3 .AND. ITURB /=8) THEN
         CALL MPI_RECV(RKW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
         CALL MPI_RECV(REPSW, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      ENDIF

      IF(TRANSL) THEN
         CALL MPI_RECV(GW,   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
         CALL MPI_RECV(RETW, S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
      ENDIF

      IF(MULPHL) THEN
         DO IPHASE = 1,NPHASES
            CALL MPI_RECV(F1RM,S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
               ALFAW(IPHASE,1:S) = F1RM(1:S) 
            CALL MPI_RECV(F1RM,S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
               XW(IPHASE,1:S) = F1RM(1:S) 
            CALL MPI_RECV(BUFFERDP,S,
     &               MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
               DTEMPW(IPHASE,1:S) = BUFFERDP(1:S) 
         ENDDO
      ENDIF


C ... First interpolate from the background block to the Chimera block 

      DO NN = 1,NTOT(M,NC)

         I = IG2 + NN - 1
      
         INTPI = MOD(INTP(I),10)

         IF(INTPI == 2 .OR. RKSI(I) > 5.0E9) THEN

         INTPB = (INTP(I)-INTPI)/10

         IF(INTPB /= NG) CYCLE

         II1I = II1(I) - IG1W + 1
         II2I = II2(I) - IG1W + 1
         II3I = II3(I) - IG1W + 1
         II4I = II4(I) - IG1W + 1

         ROFOR(I) = WGH1(I)*ROW(II1I)+WGH2(I)*ROW(II2I)
     &            + WGH3(I)*ROW(II3I)+WGH4(I)*ROW(II4I) 

         RMFOR(I) = WGH1(I)*RMW(II1I)+WGH2(I)*RMW(II2I)
     &            + WGH3(I)*RMW(II3I)+WGH4(I)*RMW(II4I) 

         RNFOR(I) = WGH1(I)*RNW(II1I)+WGH2(I)*RNW(II2I)
     &            + WGH3(I)*RNW(II3I)+WGH4(I)*RNW(II4I) 

         RWFOR(I) = WGH1(I)*RWW(II1I)+WGH2(I)*RWW(II2I)
     &            + WGH3(I)*RWW(II3I)+WGH4(I)*RWW(II4I) 

          EFOR(I) = WGH1(I)*EW(II1I)+WGH2(I)*EW(II2I)
     &            + WGH3(I)*EW(II3I)+WGH4(I)*EW(II4I) 
       
         PDFOR(I) = WGH1(I)*PDIFFW(II1I)+WGH2(I)*PDIFFW(II2I)
     &            + WGH3(I)*PDIFFW(II3I)+WGH4(I)*PDIFFW(II4I)

          VIST(I) = WGH1(I)*VISTW(II1I)+WGH2(I)*VISTW(II2I)
     &            + WGH3(I)*VISTW(II3I)+WGH4(I)*VISTW(II4I) 
       
          EPS2(I) = WGH1(I)*EPS2W(II1I)+WGH2(I)*EPS2W(II2I)
     &            + WGH3(I)*EPS2W(II3I)+WGH4(I)*EPS2W(II4I)

         IF(ITURB >= 3 .AND. ITURB /=8) THEN

            RKFOR(I) = WGH1(I)*RKW(II1I)+WGH2(I)*RKW(II2I)
     &               + WGH3(I)*RKW(II3I)+WGH4(I)*RKW(II4I) 

            REFOR(I) = WGH1(I)*REPSW(II1I)+WGH2(I)*REPSW(II2I)
     &               + WGH3(I)*REPSW(II3I)+WGH4(I)*REPSW(II4I)

         ENDIF

         IF(TRANSL) THEN

            TRM(I)%GFOR   = WGH1(I)*GW(II1I)+WGH2(I)*GW(II2I)
     &                    + WGH3(I)*GW(II3I)+WGH4(I)*GW(II4I) 

            TRM(I)%RETFOR = WGH1(I)*RETW(II1I)+WGH2(I)*RETW(II2I)
     &                    + WGH3(I)*RETW(II3I)+WGH4(I)*RETW(II4I)

         ENDIF

         IF(MULPHL) THEN

            DO IPHASE = 1,NPHASES

               MPFOR(I)%ALFAFOR  = WGH1(I)*ALFAW(IPHASE,II1I)
     &                           + WGH2(I)*ALFAW(IPHASE,II2I)
     &                           + WGH3(I)*ALFAW(IPHASE,II3I)
     &                           + WGH4(I)*ALFAW(IPHASE,II4I) 

               MPFOR(I)%XFOR     = WGH1(I)*XW(IPHASE,II1I)
     &                           + WGH2(I)*XW(IPHASE,II2I)
     &                           + WGH3(I)*XW(IPHASE,II3I)
     &                           + WGH4(I)*XW(IPHASE,II4I) 

               MPFOR(I)%DTEMPFOR = WGH1(I)*DTEMPW(IPHASE,II1I)
     &                           + WGH2(I)*DTEMPW(IPHASE,II2I)
     &                           + WGH3(I)*DTEMPW(IPHASE,II3I)
     &                           + WGH4(I)*DTEMPW(IPHASE,II4I) 

            ENDDO

         ENDIF

            IF(INTPI == 2) THEN ! This is a Chimera ghost cell         
         
               RO(I) = ROFOR(I)
               RM(I) = RMFOR(I)
               RN(I) = RNFOR(I)
               RW(I) = RWFOR(I)
                E(I) = EFOR(I)
            PDIFF(I) = PDFOR(I)
      
             TEMP(I) = WGH1(I)*TEMPW(II1I)+WGH2(I)*TEMPW(II2I)
     &               + WGH3(I)*TEMPW(II3I)+WGH4(I)*TEMPW(II4I) 

                U(I) = RM(I)/RO(I)
                V(I) = RN(I)/RO(I)
                W(I) = RW(I)/RO(I)
      
           IF(ITURB >= 3 .AND. ITURB /=8) THEN
               RK(I) = RKFOR(I)
             REPS(I) = REFOR(I)
           ENDIF
      
           IF(TRANSL) THEN
              TRM(I)%G   = TRM(I)%GFOR
              TRM(I)%RET = TRM(I)%RETFOR
           ENDIF

           IF(MULPHL) THEN

              DO IPHASE = 1,NPHASES        

                 VAR(I)%ALFA(IPHASE)  = MPFOR(I)%ALFAFOR(IPHASE)
                 VAR(I)%X(IPHASE)     = MPFOR(I)%XFOR(IPHASE)
                 PRO(I)%DTEMP(IPHASE) = MPFOR(I)%DTEMPFOR(IPHASE)

              ENDDO

           ENDIF

           ENDIF

         ENDIF

      ENDDO

      DEALLOCATE(ROW,RMW,RNW,RWW,EW,PDIFFW,TEMPW,VISTW,EPS2W)

      IF(ITURB >= 3 .AND. ITURB /=8) THEN
         DEALLOCATE(RKW,REPSW)
      ENDIF

      IF(TRANSL) THEN
         DEALLOCATE(GW,RETW)
      ENDIF

      IF(MULPHL) THEN
         DEALLOCATE(ALFAW,XW,DTEMPW,BUFFERDP)
      ENDIF

      ENDIF ! MPIC


C ... and then interpolate from the Chimera block to the background block.


      IF(MPIB) THEN  ! Background block wants data from Chimera block 

      CALL MPI_RECV(S,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(NSEATS,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)

      NS = 12*NSEATS

      ALLOCATE(RKSIW(S),TRAIN(NS))
      IF(MULPHL) THEN
         ALLOCATE(MPTRAIN(2*NSEATS*NPHASES),MPTRAINDP(NSEATS*NPHASES))
      ENDIF
      ALLOCATE(II1W(S),II2W(S),II3W(S),II4W(S),INTPW(S))
      ALLOCATE(WGH1W(S),WGH2W(S),WGH3W(S),WGH4W(S))

      CALL MPI_RECV(RKSIW,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II1W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II2W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II3W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(II4W, S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(INTPW,S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH1W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH2W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH3W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)
      CALL MPI_RECV(WGH4W,S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,st,rc)


      L1  = 0
      LMP = 0

      DO I = 1,S
      
         IF(RKSIW(I) > 5.0E9) THEN

         IF(II1W(I) == 0) CYCLE

         INTPI = MOD(INTPW(I),10)
         INTPB = (INTPW(I)-INTPI)/10

         IF(INTPB /= NCG) CYCLE

         II1I = II1W(I)
         II2I = II2W(I)
         II3I = II3W(I)
         II4I = II4W(I)

         L1  = L1  + 1
         L2  = L1  + NSEATS
         L3  = L2  + NSEATS
         L4  = L3  + NSEATS
         L5  = L4  + NSEATS
         L6  = L5  + NSEATS
         L7  = L6  + NSEATS
         L8  = L7  + NSEATS
         L9  = L8  + NSEATS
         L10 = L9  + NSEATS
         L11 = L10 + NSEATS
         L12 = L11 + NSEATS
         LMP = LMP + 1


         TRAIN(L1) = WGH1W(I)*RO(II1I)+WGH2W(I)*RO(II2I)
     &             + WGH3W(I)*RO(II3I)+WGH4W(I)*RO(II4I) 

         TRAIN(L2) = WGH1W(I)*RM(II1I)+WGH2W(I)*RM(II2I)
     &             + WGH3W(I)*RM(II3I)+WGH4W(I)*RM(II4I) 

         TRAIN(L3) = WGH1W(I)*RN(II1I)+WGH2W(I)*RN(II2I)
     &             + WGH3W(I)*RN(II3I)+WGH4W(I)*RN(II4I) 

         TRAIN(L4) = WGH1W(I)*RW(II1I)+WGH2W(I)*RW(II2I)
     &             + WGH3W(I)*RW(II3I)+WGH4W(I)*RW(II4I) 

         TRAIN(L5) = WGH1W(I)*E(II1I)+WGH2W(I)*E(II2I)
     &             + WGH3W(I)*E(II3I)+WGH4W(I)*E(II4I) 
       
         TRAIN(L6) = WGH1W(I)*PDIFF(II1I)+WGH2W(I)*PDIFF(II2I)
     &             + WGH3W(I)*PDIFF(II3I)+WGH4W(I)*PDIFF(II4I)

         TRAIN(L7) = WGH1W(I)*VIST(II1I)+WGH2W(I)*VIST(II2I)
     &             + WGH3W(I)*VIST(II3I)+WGH4W(I)*VIST(II4I) 
       
         TRAIN(L8) = WGH1W(I)*EPS2(II1I)+WGH2W(I)*EPS2(II2I)
     &             + WGH3W(I)*EPS2(II3I)+WGH4W(I)*EPS2(II4I)

         IF(ITURB >= 3 .AND. ITURB /=8) THEN

            TRAIN(L9)  = WGH1W(I)*RK(II1I)+WGH2W(I)*RK(II2I)
     &                 + WGH3W(I)*RK(II3I)+WGH4W(I)*RK(II4I) 

            TRAIN(L10) = WGH1W(I)*REPS(II1I)+WGH2W(I)*REPS(II2I)
     &                 + WGH3W(I)*REPS(II3I)+WGH4W(I)*REPS(II4I)

         ENDIF

         IF(TRANSL) THEN

            TRAIN(L11) = WGH1W(I)*TRM(II1I)%G+WGH2W(I)*TRM(II2I)%G
     &                 + WGH3W(I)*TRM(II3I)%G+WGH4W(I)*TRM(II4I)%G 

            TRAIN(L12) = WGH1W(I)*TRM(II1I)%RET+WGH2W(I)*TRM(II2I)%RET
     &                 + WGH3W(I)*TRM(II3I)%RET+WGH4W(I)*TRM(II4I)%RET

         ENDIF

         IF(MULPHL) THEN

            DO IPHASE = 1,NPHASES

               IMP = LMP + (IPHASE-1)*NSEATS

               MPTRAIN(IMP) = WGH1W(I)*VAR(II1I)%ALFA(IPHASE)
     &                      + WGH2W(I)*VAR(II2I)%ALFA(IPHASE)
     &                      + WGH3W(I)*VAR(II3I)%ALFA(IPHASE)
     &                      + WGH4W(I)*VAR(II4I)%ALFA(IPHASE) 

               IMP = LMP + (IPHASE-0)*NSEATS

               MPTRAIN(IMP) = WGH1W(I)*VAR(II1I)%X(IPHASE)
     &                      + WGH2W(I)*VAR(II2I)%X(IPHASE)
     &                      + WGH3W(I)*VAR(II3I)%X(IPHASE)
     &                      + WGH4W(I)*VAR(II4I)%X(IPHASE) 

               IMP = LMP + (IPHASE-1)*NSEATS

               MPTRAINDP(IMP) = WGH1W(I)*PRO(II1I)%DTEMP(IPHASE)
     &                        + WGH2W(I)*PRO(II2I)%DTEMP(IPHASE)
     &                        + WGH3W(I)*PRO(II3I)%DTEMP(IPHASE)
     &                        + WGH4W(I)*PRO(II4I)%DTEMP(IPHASE)

            ENDDO

         ENDIF

         ENDIF

      ENDDO

       
C ... Send the fortify arrays back to the background block process 
C ... and deallocate the working arrays.

      CALL MPI_SEND(TRAIN,NS,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      IF(MULPHL) THEN
         CALL MPI_SEND(TRAIN,2*NPHASES*NSEATS,MPI_REAL8,
     &        pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(MPTRAINDP,NPHASES*NSEATS,MPI_REAL8,
     &        pid,tag,MPI_COMM_WORLD,rc)
      ENDIF

      DEALLOCATE(TRAIN,RKSIW)
      IF(MULPHL) DEALLOCATE(MPTRAIN,MPTRAINDP)
      DEALLOCATE(II1W,II2W,II3W,II4W,INTPW)
      DEALLOCATE(WGH1W,WGH2W,WGH3W,WGH4W)

      ENDIF ! MPIB


      ELSEIF(N /= -1) THEN  ! Only the background block is local

      pid    = NPNUM(NCG) - 1 
      tag    = NG

      NSEATS = 0
      MPIB   = .FALSE.
      
      DO NN = 1,NTOT(M,N)
         I  = IG1 + NN - 1 
         INTPI = MOD(INTP(I),10)
         INTPB = (INTP(I)-INTPI)/10
         IF(INTPB == NCG .AND. RKSI(I) > 5.0E9) THEN
            NSEATS = NSEATS + 1
         ENDIF  
      ENDDO

      MPIB = NSEATS > 0

      CALL MPI_SEND(MPIB,1,MPI_LOGICAL,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_RECV(MPIC,1,MPI_LOGICAL,pid,tag,MPI_COMM_WORLD,st,rc)


      IF(MPIC) THEN ! Chimera block wants data from background

      S = NTOT(M,N)

      CALL MPI_SEND(S,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)

      CALL MPI_SEND(RO(IG1),   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(RM(IG1),   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(RN(IG1),   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(RW(IG1),   S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(E(IG1),    S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(PDIFF(IG1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(TEMP(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(VIST(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(EPS2(IG1), S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(IG1,      1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
                                                            
      IF(ITURB >= 3 .AND. ITURB /=8) THEN                                 
         CALL MPI_SEND(RK(IR1),  S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(REPS(IR1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      ENDIF
                                                            
      IF(TRANSL) THEN
         DO I = IR1,IR1+S-1
            F1RM(I) = TRM(I)%G
            F1RN(I) = TRM(I)%RET
         ENDDO
         CALL MPI_SEND(F1RM(IR1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
         CALL MPI_SEND(F1RN(IR1),S,MPI_REAL8,pid,tag,MPI_COMM_WORLD,rc)
      ENDIF
                                                            
      IF(MULPHL) THEN

         DO IPHASE = 1,NPHASES

            DO I = IG1,IG1+S-1
               F1RM(I) = VAR(I)%ALFA(IPHASE)
            ENDDO
            CALL MPI_SEND(F1RM(IG1),S,MPI_REAL8,
     &                    pid,tag,MPI_COMM_WORLD,rc)

            DO I = IG1,IG1+S-1
               F1RM(I) = VAR(I)%X(IPHASE)
            ENDDO
            CALL MPI_SEND(F1RM(IG1),S,MPI_REAL8,
     &                    pid,tag,MPI_COMM_WORLD,rc)

            DO I = IG1,IG1+S-1
               F1RM(I) = PRO(I)%DTEMP(IPHASE)
            ENDDO
            CALL MPI_SEND(F1RM(IG1),S,MPI_REAL8,
     &                    pid,tag,MPI_COMM_WORLD,rc)

         ENDDO

      ENDIF

      ENDIF ! MPIC

      
      IF(MPIB) THEN  ! Background block wants data from Chimera block 

      S = NTOT(M,N)

      CALL MPI_SEND(S,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(NSEATS,1,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)

      CALL MPI_SEND(RKSI(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II1(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II2(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II3(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(II4(IG1), S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(INTP(IG1),S,MPI_INTEGER,pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH1(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH2(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH3(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)
      CALL MPI_SEND(WGH4(IG1),S,MPI_REAL8,   pid,tag,MPI_COMM_WORLD,rc)


C ... Receive the updated fortify arrays from the Chimera block process. 

      NS = 12*NSEATS
      ALLOCATE(TRAIN(NS))
      IF(MULPHL) THEN
         ALLOCATE(MPTRAIN(2*NSEATS*NPHASES),MPTRAINDP(NSEATS*NPHASES))
      ENDIF
 
      CALL MPI_RECV(TRAIN,NS,MPI_REAL8,pid,tag,MPI_COMM_WORLD,st,rc)
    
      IF(MULPHL) THEN
         CALL MPI_RECV(MPTRAIN,2*NPHASES*NSEATS,MPI_REAL8,
     &                 pid,tag,MPI_COMM_WORLD,st,rc)
         CALL MPI_RECV(MPTRAINDP,NPHASES*NSEATS,MPI_REAL8,
     &                 pid,tag,MPI_COMM_WORLD,st,rc)
      ENDIF

      L1  = 0
      LMP = 0

      DO NN = 1,NTOT(M,N)

         I = IG1 + NN - 1
         J = IR1 + NN - 1
      
         INTPI = MOD(INTP(I),10)

         IF(RKSI(I) > 5.0E9) THEN

         INTPB = (INTP(I)-INTPI)/10

         IF(INTPB /= NCG) CYCLE

         L1  = L1  + 1
         L2  = L1  + NSEATS
         L3  = L2  + NSEATS
         L4  = L3  + NSEATS
         L5  = L4  + NSEATS
         L6  = L5  + NSEATS
         L7  = L6  + NSEATS
         L8  = L7  + NSEATS
         L9  = L8  + NSEATS
         L10 = L9  + NSEATS
         L11 = L10 + NSEATS
         L12 = L11 + NSEATS
         LMP = LMP + 1

         ROFOR(I) = TRAIN(L1)
         RMFOR(I) = TRAIN(L2) 
         RNFOR(I) = TRAIN(L3) 
         RWFOR(I) = TRAIN(L4) 
          EFOR(I) = TRAIN(L5) 
         PDFOR(I) = TRAIN(L6) 
          VIST(I) = TRAIN(L7) 
          EPS2(I) = TRAIN(L8) 

         IF(ITURB >= 3 .AND. ITURB /=8) THEN
            RKFOR(J) = TRAIN(L9)
            REFOR(J) = TRAIN(L10)
         ENDIF

         IF(TRANSL) THEN
            TRM(J)%GFOR   = TRAIN(L11)
            TRM(J)%RETFOR = TRAIN(L12)
         ENDIF

         IF(MULPHL) THEN

            DO IPHASE = 1,NPHASES 

               IMP = LMP + (IPHASE-1)*NSEATS
               MPFOR(I)%ALFAFOR(IPHASE)  = MPTRAIN(IMP)

               IMP = LMP + (IPHASE-0)*NSEATS
               MPFOR(I)%XFOR(IPHASE)     = MPTRAIN(IMP)

               IMP = LMP + (IPHASE-1)*NSEATS
               MPFOR(I)%DTEMPFOR(IPHASE) = MPTRAINDP(IMP)

            ENDDO

         ENDIF

         ENDIF

      ENDDO

      DEALLOCATE(TRAIN)
      IF(MULPHL) DEALLOCATE(MPTRAIN,MPTRAINDP)

      ENDIF ! MPIB


      ENDIF  ! Chimera block or background block not local

      ENDIF  ! Chimera priority

      ENDDO  ! Background block loop

      ENDDO  ! Chimera block loop


************************************************************************
C ... Set the flow field slowly to rest or to grid speed in cells 
C ... that lie inside the structure.

* This should be further tested. 
*
       DO N = 1,NBLOCK
 
         NGL = NPROCE(1+N,IPRO) !Global block number
         IG1 = IG(1,N)

         FRSDEN  = BLKS(NGL)%FRSDEN
         FRSPRE  = BLKS(NGL)%FRSPRE
         FRSVEL  = BLKS(NGL)%FRSVEL
         FRSVIS  = BLKS(NGL)%FRSVIS
         FRSSIE  = BLKS(NGL)%FRSSIE
         RKLIM   = BLKS(NGL)%RKLIM
         EPSLIM  = BLKS(NGL)%EPSLIM
         FRSMUT  = BLKS(NGL)%FRSMUT
 
c         CALL ZEROVELO(U(IG1),V(IG1),W(IG1),XC(IG1),YC(IG1),ZC(IG1),
c     +   INTP(IG1),BLANK(IG1),
c     +   OMEGA(NGL),OMEGAX(NGL),OMEGAY(NGL),OMEGAZ(NGL),CENAX(NGL),
c     +   CENAY(NGL),CENAZ(NGL),IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
c     +   NGL)
         
         IF (BLKS(NGL)%ZEROV) THEN
            CALL ZEROVAR(U(IG1),V(IG1),W(IG1),XC(IG1),YC(IG1),ZC(IG1),
     +           INTP(IG1),BLANK(IG1),OMEGA(NGL),OMEGAX(NGL),
     +           OMEGAY(NGL),OMEGAZ(NGL),CENAX(NGL),CENAY(NGL),
     +           CENAZ(NGL),
     +           VIST(IG1),EPS2(IG1),RO(IG1),P(IG1),PDIFF(IG1),
     +           TEMP(IG1),RK(IG1),
     +           REPS(IG1),FRSDEN,FRSPRE,FRSSIE,RKLIM,EPSLIM,FRSMUT,
     +           FRSVIS,FRSTEM,IMAX(1,N),JMAX(1,N),KMAX(1,N),
     +           IN,JN,KN,NGL)
         ENDIF
 
       ENDDO
************************************************************************

 
C ... Print fortify functions

      DO N = 1,NBLOCK
         IG1 = IG(M,N)      
         IR1 = IR(M,N)      
         NGL = NPROCE(1+N,IPRO)  ! Global block number

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)
         WRITE(3,*)' RESULTS FORTIFY FUNCTIONS'
         CALL PRINYS(3, RKSIC, RKSI(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)         
         CALL PRINYS(3,ROFORC,ROFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,RMFORC,RMFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,RNFORC,RNFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,RWFORC,RWFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3, EFORC, EFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,PDFORC,PDFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)

         IF(ITURB >= 3 .AND. ITURB /=8) THEN
         CALL PRINYS(3,RKFORC,RKFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         CALL PRINYS(3,REFORC,REFOR(IG1),IT(M,N),IL(M,N),0,IK(M,N),
     +        IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         ENDIF

         IF(TRANSL) THEN
            ZZZ(1:NTOT(M,N)) = TRM(IG1:IG1+NTOT(M,N)-1)%GFOR
            CALL PRINYS(3,GFORC,ZZZ(1),IT(M,N),IL(M,N),0,
     +           IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
            ZZZ(1:NTOT(M,N)) = TRM(IG1:IG1+NTOT(M,N)-1)%RETFOR
            CALL PRINYS(3,RETFORC,ZZZ(1),IT(M,N),IL(M,N),0,
     +           IK(M,N),IMAX(M,N),JMAX(M,N),KMAX(M,N),NGL,M)
         ENDIF

      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CHIMUP
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE DONORS(ITURB,LEVEL,NCHIMT,TRUE_DISTL,
     &     NG,IMAX,JMAX,KMAX,IN,JN,KN,X,Y,Z,BLANK,DISTW,LOCDIS,
     &     VOL,ICON,NPATCH,M,
     &     NCG,IMAXC,JMAXC,KMAXC,XC,YC,ZC,BLANKC,DISTWC,LOCDISC,
     &     VOLC,ICONC,NPATCC,
     &     RKSI,RKSIC,RKSIV,
     &     II1, II2, II3, II4, INTP, 
     &     II1C,II2C,II3C,II4C,INTPC,
     &     WGH1, WGH2, WGH3, WGH4, IG1,
     &     WGH1C,WGH2C,WGH3C,WGH4C,IG1C,NOFDONORS)

C ... This subroutine 'interpolates' the flow field properties 
C ... between two grids. In fact, the interpolation is done in subroutine
C ... CHIMUP. This subroutine only searches the dominating cells and saves
C ... the corresponding weight factors.
C ...  
C ... The 'interpolation' is done in three phases. 
C ... In the first phase the values from the background grid are 
C ... interpolated to the ghost cells of the Chimera grid.
C ... In the second phase the values from the background grid are 
C ... interpolated to the Chimera grid and in the final
C ... phase the solution is interpolated from the Chimera grid to
C ... the background grid. The domination criterion is based on the 
C ... local nearest wall distance. The search process needed in the
C ... interpolation is accelerated by sharing the starting grid
C ... cells into small groups according to their locations in
C ... space. Note that these 'cells' are not actual cells, but 
C ... cells formed by the center points of the actual cells. 

      USE NS3CO, ONLY : IC9
      USE MAIN_ARRAYS, ONLY : NLOCAL, BXMIN, BXMAX, BYMIN, BYMAX,
     &                                BZMIN, BZMAX
      
      IMPLICIT NONE

      INTEGER, ALLOCATABLE :: IBOXES(:,:,:,:)  ! Cell sorting grid
      INTEGER, ALLOCATABLE :: KSC(:)           ! Cell pointer vector
      INTEGER, ALLOCATABLE :: NEARWL(:)        ! Wall vicinity vector

      INTEGER, DIMENSION(IC9,*) :: ICON, ICONC
      INTEGER, DIMENSION(8) :: II

      INTEGER :: NXBOX, NYBOX, NZBOX
      INTEGER :: IN, JN, KN, ITURB, LEVEL, NCB
      INTEGER :: IMAX, JMAX, KMAX, IMAXC, JMAXC, KMAXC
      INTEGER :: M, N, NC, NPATCH, NPATCC, NCHN, NCHNC
      INTEGER :: ISTRID, JSTRID, KSTRID, ISTRIC, JSTRIC, KSTRIC
      INTEGER :: IFIRST, ILAST, JFIRST, JLAST, KFIRST, KLAST
      INTEGER :: IA, IM, JA, JM, IA2, IM2, JA2, JM2
      INTEGER :: I, J, K, L, I1, J1, K1, LBS, LBE, KP, KR
      INTEGER :: IC, KA, JJ, JC, KC, LC, IG1, IG1C, NG, NCG
      INTEGER :: INTPKR, INTPCKR, IFACE, IERRCODE, LCE, NOFDONORS

      REAL, DIMENSION(*) :: 
     &     RKSI,RKSIC,BLANK,BLANKC,
     &     VOL,VOLC,DISTW,DISTWC

      REAL, DIMENSION(*) :: 
     &     WGH1,WGH2,WGH3,WGH4,WGH1C,WGH2C,WGH3C,WGH4C

      INTEGER, DIMENSION(*) :: LOCDIS, LOCDISC, NCHIMT,
     &                         II1,  II2,  II3,  II4,  INTP,   
     &                         II1C, II2C, II3C, II4C, INTPC

      REAL, DIMENSION(3) :: XX
      REAL, DIMENSION(4) :: A0

      REAL    :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      REAL    :: DX, DY, DZ, DXR, DYR, DZR
      REAL    :: RKSIV, VTETRA, VRATIO, DCHIM

      REAL, DIMENSION(*) :: X,Y,Z,XC,YC,ZC

      LOGICAL :: ITISIN, TWODIM, ALLCORNERS, DOMINATE, LDCHIM,
     &           TRUE_DISTL, CHIMGC, INSIDE, connected

      NOFDONORS = 0

      NC    = NLOCAL(NCG)
      NCHN  = NCHIMT(NG)
      NCHNC = NCHIMT(NCG)

C ... Cell volume ratio used as a 'strong' blanking criteria.

C ... Sopiva arvo on tapauskohtainen ja hilataso (LEVEL) vaikuttaa siihen. 
C ... Periaatteessa: harvempi hila -> suurempi VRATIO.
C ... Arvo vaikuttaa vain jalkikasittelya varten tehtavaan blankkaukseen,
C ... ei itse virtausratkaisuun. 
C ... 20.10.2000: Nyt kiinteiden pintojen sisaan jaavat pisteet 
C ... voidaan poistaa/poistetaan tarkasti, joten 'voimakas blankkaus' on
C ... tarpeeton. VRATIOnkaan arvoa ei siis tarvitse pohtia. Toistaiseksi
C ... 'voimakas blankkaus' on kuitenkin aktivoituna, koska se joissakin
C ... tilanteissa auttaa 'siistimpaan' blankkaukseen.

      VRATIO = 10.0

C ... Two-dimensional case ?

      TWODIM = (KMAX == 1 .AND. KMAXC == 1)


C ... Allocate memory for the cell sorting grid. Use the size of the
C ... chimera block as an 'optimum' for the sorting grid size. 

*      NXBOX = 100
*      NYBOX = 100
*      NZBOX = 100

      NXBOX = IMAXC
      NYBOX = JMAXC
      NZBOX = KMAXC

      ALLOCATE(IBOXES(NXBOX,NYBOX,NZBOX,2),STAT=IERRCODE)

C ... Not enough memory ?                                                      

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'DONORS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    


C ... Maximum dimensions of the Chimera grid

*      CALL BOUNDS(IMAXC,JMAXC,KMAXC,XC,YC,ZC,
*     &            XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,IN,JN,KN)

      XMIN = BXMIN(NCG)
      XMAX = BXMAX(NCG)
      YMIN = BYMIN(NCG)
      YMAX = BYMAX(NCG)
      ZMIN = BZMIN(NCG)
      ZMAX = BZMAX(NCG)


C ... Box size (cells will be grouped into these 'boxes').
C ... Stretch the sorting grid slightly.

      DX  = XMAX-XMIN
      DY  = YMAX-YMIN
      DZ  = ZMAX-ZMIN

      XMIN = XMIN - 0.05*DX
      XMAX = XMAX + 0.05*DX
      YMIN = YMIN - 0.05*DY
      YMAX = YMAX + 0.05*DY
      ZMIN = ZMIN - 0.05*DZ
      ZMAX = ZMAX + 0.05*DZ

      DX  = (XMAX-XMIN)/NXBOX
      DY  = (YMAX-YMIN)/NYBOX
      DZ  = (ZMAX-ZMIN)/NZBOX

      DXR = 1.0/DX
      DYR = 1.0/DY
      DZR = 1.0/DZ


C ... Calculate the total length of the cell pointer vector in order
C ... to allocate enough memory.

      ALLOCATE(KSC(1))

      CALL BOXIT(IMAX,JMAX,KMAX,IN,JN,KN,X,Y,Z,KSC,NCB,
     &           IBOXES,NXBOX,NYBOX,NZBOX,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)

C ... Allocate memory for the cell pointer vector KSC

      DEALLOCATE(KSC)
      ALLOCATE(KSC(NCB),STAT=IERRCODE)

C ... Not enough memory ?                                                      

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'DONORS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    


C ... Share the cells into groups ('boxes')

      CALL BOXIT(IMAX,JMAX,KMAX,IN,JN,KN,X,Y,Z,KSC,NCB,
     &           IBOXES,NXBOX,NYBOX,NZBOX,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)

C ######################################################################
C
C Interpolate from the background grid to the Chimera's ghost cells.   
C
C ######################################################################

      ISTRID = IMAX  + 2*IN
      JSTRID = JMAX  + 2*JN
      KSTRID = KMAX  + 2*KN

      ISTRIC = IMAXC + 2*IN
      JSTRIC = JMAXC + 2*JN
      KSTRIC = KMAXC + 2*KN

      DCHIM  = 1.0E+20
      LDCHIM = .FALSE.

      LCE    = 0


      connected = .false.

      if(nchn == -1 .and. nchnc == -1) then
         do l = 1,NPATCC
            IF(ICONC(1,L) == 1 .and. ICONC(28,L) == NG) then
               connected = .true.
            endif
         enddo
      endif


      DO 1000 L = 1,NPATCC

      IF(ICONC(1,L) == 12) THEN    ! Boundary condition type CHI

C ... Extend the chimera patch
         
      CALL PATCHE(ICONC,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IFACE = ICONC(3,L)

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
         IFIRST = IMAXC + 1
         ILAST  = IMAXC + 2
         JFIRST = IA
         JLAST  = IM
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 5) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = JMAXC + 1
         JLAST  = JMAXC + 2
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 6) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = JA
         JLAST  = JM
         KFIRST = KMAXC + 1
         KLAST  = KMAXC + 2
      ENDIF

      DO KC = KFIRST,KLAST
         KA = (KN+KC-1)*ISTRIC*JSTRIC
      DO JC = JFIRST,JLAST
         JJ = (JN+JC-1)*ISTRIC + IN + KA 
      DO IC = IFIRST,ILAST

      KR = JJ + IC

      IF(MOD(INTPC(KR),10) /= 2) INTPC(KR) = 10*NG + 3

C ... An estimation for the chimera block size
      IF(LOCDISC(KR) == NC .AND. DISTWC(KR) < DCHIM) THEN
         DCHIM = DISTWC(KR)
         LDCHIM = .TRUE.
      ENDIF

      XX(1) = XC(KR) 
      XX(2) = YC(KR) 
      XX(3) = ZC(KR)

C ... In a two-dimensional case we can choose any location in the
C ... z-direction as long as the point certainly lies within the
C ... boundaries. By this we try to avoid a situation where the
C ... target point would lie exactly on the face of the donor cell.  

      IF(TWODIM) XX(3) = 0.51*(ZMIN+ZMAX)

C ... Find the right 'box'

      I1 = MIN(INT((XX(1)-XMIN)*DXR) + 1,NXBOX) 
      J1 = MIN(INT((XX(2)-YMIN)*DYR) + 1,NYBOX)
      K1 = MIN(INT((XX(3)-ZMIN)*DZR) + 1,NZBOX)
      I1 = MAX(I1,1)
      J1 = MAX(J1,1)
      K1 = MAX(K1,1)


C ... The cell pointer range in this 'box'

      LBS = IBOXES(I1,J1,K1,2)                             ! Start
      LBE = IBOXES(I1,J1,K1,2) + IBOXES(I1,J1,K1,1) - 1    ! End

      DO LC = LBS,LBE

         KP = KSC(LC)

         CALL 
     &   ISITIN(XX,A0,KP,II,ISTRID,JSTRID,KSTRID,X,Y,Z,ITISIN,VTETRA)
         IF(ITISIN) THEN
            LCE = LC
            GOTO 313
         ENDIF

      ENDDO


      KP = 0

 313  CONTINUE

         IF(KP /= 0) THEN

         INSIDE = BLANK(II(1)) == -8.0 
     &       .OR. BLANK(II(2)) == -8.0
     &       .OR. BLANK(II(3)) == -8.0
     &       .OR. BLANK(II(4)) == -8.0 
            IF(.NOT.INSIDE .and..not.connected) THEN
                NOFDONORS = NOFDONORS + 1
                 II1C(KR) = II(1) + IG1 - 1             
                 II2C(KR) = II(2) + IG1 - 1             
                 II3C(KR) = II(3) + IG1 - 1             
                 II4C(KR) = II(4) + IG1 - 1
                INTPC(KR) = 10*NG + 2
                WGH1C(KR) = A0(1)             
                WGH2C(KR) = A0(2)             
                WGH3C(KR) = A0(3)             
                WGH4C(KR) = A0(4)
               BLANKC(KR) = 1.0
            ENDIF

         ENDIF

      ENDDO ! i-loop (the chimera block)
      ENDDO ! j-loop (the chimera block)
      ENDDO ! k-loop (the chimera block)

      ENDIF

 1000 CONTINUE

      IF(NCHN == -1) THEN
         DEALLOCATE(KSC)
         RETURN
      ENDIF

C ######################################################################
C
C Interpolate from the background grid to the Chimera grid.   
C
C ######################################################################

C ... Allocate and create the wall vicinity vector

      ALLOCATE(NEARWL(ISTRID*JSTRID*KSTRID),STAT=IERRCODE)
      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'DONORS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    
      CALL WALLVI(IMAX,JMAX,KMAX,ICON,NPATCH,NEARWL,TWODIM,IN,JN,KN)


      LCE   = 0

      DO KC = 0,KMAXC+1
         KA = (KN+KC-1)*ISTRIC*JSTRIC
      DO JC = 0,JMAXC+1
         JJ = (JN+JC-1)*ISTRIC + IN + KA 
      DO IC = 0,IMAXC+1

      KR = JJ + IC

      IF(MOD(INTPC(KR),10) == 0) INTPC(KR) = 10*NG + 4 

      XX(1) = XC(KR) 
      XX(2) = YC(KR) 
      XX(3) = ZC(KR)
      IF(TWODIM) XX(3) = 0.51*(ZMIN+ZMAX)

C ... Find the right 'box'

      I1 = MIN(INT((XX(1)-XMIN)*DXR) + 1,NXBOX) 
      J1 = MIN(INT((XX(2)-YMIN)*DYR) + 1,NYBOX)
      K1 = MIN(INT((XX(3)-ZMIN)*DZR) + 1,NZBOX)
      I1 = MAX(I1,1)
      J1 = MAX(J1,1)
      K1 = MAX(K1,1)


C ... The cell pointer range in this 'box'

      LBS = IBOXES(I1,J1,K1,2)                             ! Start
      LBE = IBOXES(I1,J1,K1,2) + IBOXES(I1,J1,K1,1) - 1    ! End

      DO LC = LBS,LBE

         KP = KSC(LC)

         CALL 
     &   ISITIN(XX,A0,KP,II,ISTRID,JSTRID,KSTRID,X,Y,Z,ITISIN,VTETRA)
         IF(ITISIN) THEN
            LCE = LC
            GOTO 314
         ENDIF

      ENDDO

      KP = 0

 314  CONTINUE

         INTPCKR = MOD(INTPC(KR),10)
         IF(KP /= 0 .AND. INTPCKR /= 2 .AND. RKSIC(KR) == 0.0) THEN

            INSIDE = BLANK(II(1)) == -8.0 
     &          .OR. BLANK(II(2)) == -8.0
     &          .OR. BLANK(II(3)) == -8.0
     &          .OR. BLANK(II(4)) == -8.0

            IF(.NOT.INSIDE) THEN
               NOFDONORS = NOFDONORS + 1
                II1C(KR) = II(1) + IG1 - 1             
                II2C(KR) = II(2) + IG1 - 1             
                II3C(KR) = II(3) + IG1 - 1             
                II4C(KR) = II(4) + IG1 - 1
               INTPC(KR) = 10*NG + 1
               WGH1C(KR) = A0(1)             
               WGH2C(KR) = A0(2)             
               WGH3C(KR) = A0(3)             
               WGH4C(KR) = A0(4)             
            ENDIF

            IF(TRUE_DISTL) THEN

                DOMINATE =  NCHIMT(LOCDIS(II(1))) == NCHN
     &                .AND. NCHIMT(LOCDIS(II(2))) == NCHN
     &                .AND. NCHIMT(LOCDIS(II(3))) == NCHN
     &                .AND. NCHIMT(LOCDIS(II(4))) == NCHN
     &                .AND. NCHIMT(LOCDISC(KR))   == NCHN

* This criterion was better in some 2D tests
*                DOMINATE = (NCHIMT(LOCDIS(II(1))) == NCHN
*     &                 .OR. NCHIMT(LOCDIS(II(2))) == NCHN
*     &                 .OR. NCHIMT(LOCDIS(II(3))) == NCHN
*     &                 .OR. NCHIMT(LOCDIS(II(4))) == NCHN)
*     &                .AND. NCHIMT(LOCDISC(KR))   == NCHN

            ELSE

              ALLCORNERS = NCHIMT(LOCDIS(II(1))) == NCHN
     &               .AND. NCHIMT(LOCDIS(II(2))) == NCHN
     &               .AND. NCHIMT(LOCDIS(II(3))) == NCHN
     &               .AND. NCHIMT(LOCDIS(II(4))) == NCHN
                DOMINATE = (ALLCORNERS
     &               .AND. NCHIMT(LOCDISC(KR))   == NCHN)
     &               .OR.  (ALLCORNERS
     &               .AND. NCHIMT(LOCDISC(KR))   /= NCHN
     &               .AND. DISTWC(KR) > MAX(DISTW(II(1)),DISTW(II(2)),
     &                                      DISTW(II(3)),DISTW(II(4))))

            ENDIF


            IF(DOMINATE .AND. .NOT.INSIDE .AND. NCHNC < 900) THEN
               RKSIC(KR)  = RKSIV
               BLANKC(KR) = 0.0
            ENDIF

*            IF (NCHIMT(LOCDIS(II(1))) /= NCHNC
*     &    .AND. NCHIMT(LOCDIS(II(2))) /= NCHNC
*     &    .AND. NCHIMT(LOCDIS(II(3))) /= NCHNC
*     &    .AND. NCHIMT(LOCDIS(II(4))) /= NCHNC) BLANKC(KR) = 0.0

            IF(.NOT.TRUE_DISTL) THEN
            IF(NEARWL(KP) == 1 .AND. VOLC(KR)/VOL(KP) > VRATIO) THEN 
               BLANKC(KR) = -8.0 ! Strong blanking near walls
            ENDIF
            ENDIF

         ENDIF     ! KP/=0

      ENDDO ! i-loop (the chimera block)
      ENDDO ! j-loop (the chimera block)
      ENDDO ! k-loop (the chimera block)

      DEALLOCATE(NEARWL)


C ######################################################################
C
C Interpolate from the Chimera grid to the background grid.   
C
C ######################################################################

C ... Calculate the total length of the cell pointer vector in order
C ... to allocate enough memory.


      CALL BOXIT(IMAXC,JMAXC,KMAXC,IN,JN,KN,XC,YC,ZC,KSC,NCB,
     &           IBOXES,NXBOX,NYBOX,NZBOX,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)


C ... Allocate memory for the cell pointer vector KSC

      DEALLOCATE(KSC)
      ALLOCATE(KSC(NCB),STAT=IERRCODE)
    
C ... Not enough memory ?                                                      

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'DONORS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    


C ... Share the cells into groups ('boxes')

      CALL BOXIT(IMAXC,JMAXC,KMAXC,IN,JN,KN,XC,YC,ZC,KSC,NCB,
     &           IBOXES,NXBOX,NYBOX,NZBOX,
     &           XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)

C ... Allocate and create the wall vicinity vector

      ALLOCATE(NEARWL(ISTRIC*JSTRIC*KSTRIC),STAT=IERRCODE)
      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'DONORS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    
      CALL WALLVI(IMAXC,JMAXC,KMAXC,ICONC,NPATCC,NEARWL,TWODIM,IN,JN,KN)

      LCE    = 0


      DO K  = 0,KMAX+1
         KA = (KN+K-1)*ISTRID*JSTRID
      DO J  = 0,JMAX+1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 0,IMAX+1

      KR = JJ + I

      IF(TRUE_DISTL) THEN

         IF(MOD(INTP(KR),10) == 0) INTP(KR) = 10*NCG + 5

      ELSE

          IF(MOD(INTP(KR),10) == 0
     &   .AND. NCHIMT(LOCDIS(KR)) > NCHN
     &   .AND. LDCHIM    
     &   .AND. DISTW(KR) < 0.25*DCHIM) INTP(KR) = 10*NCG + 5 

      ENDIF  

      XX(1) = X(KR) 
      XX(2) = Y(KR) 
      XX(3) = Z(KR)

      IF(TWODIM) XX(3) = 0.51*(ZMIN+ZMAX)

      IF(XX(1) < XMIN) CYCLE
      IF(XX(2) < YMIN) CYCLE
      IF(XX(3) < ZMIN) CYCLE
      IF(XX(1) > XMAX) CYCLE
      IF(XX(2) > YMAX) CYCLE
      IF(XX(3) > ZMAX) CYCLE

C ... Find the right 'box'

      I1 = MIN(INT((XX(1)-XMIN)*DXR) + 1,NXBOX) 
      J1 = MIN(INT((XX(2)-YMIN)*DYR) + 1,NYBOX)
      K1 = MIN(INT((XX(3)-ZMIN)*DZR) + 1,NZBOX)
      I1 = MAX(I1,1)
      J1 = MAX(J1,1)
      K1 = MAX(K1,1)


C ... The cell pointer range in this 'box'

      LBS = IBOXES(I1,J1,K1,2)                             ! Start
      LBE = IBOXES(I1,J1,K1,2) + IBOXES(I1,J1,K1,1) - 1    ! End

      DO LC = LBS,LBE

         KP = KSC(LC)

         CALL 
     &   ISITIN(XX,A0,KP,II,ISTRIC,JSTRIC,KSTRIC,XC,YC,ZC,ITISIN,VTETRA)
         IF(ITISIN) THEN
            LCE = LC
            GOTO 315
         ENDIF

      ENDDO

      KP = 0

 315  CONTINUE

         INTPKR = MOD(INTP(KR),10)
         IF(KP /= 0 .AND. INTPKR /= 2 .AND. RKSI(KR) == 0.0) THEN

            INSIDE = BLANKC(II(1)) == -8.0 
     &          .OR. BLANKC(II(2)) == -8.0
     &          .OR. BLANKC(II(3)) == -8.0
     &          .OR. BLANKC(II(4)) == -8.0

            IF(.NOT.INSIDE) THEN
              NOFDONORS = NOFDONORS + 1
                II1(KR) = II(1) + IG1C - 1             
                II2(KR) = II(2) + IG1C - 1             
                II3(KR) = II(3) + IG1C - 1             
                II4(KR) = II(4) + IG1C - 1
               INTP(KR) = 10*NCG + 1
               WGH1(KR) = A0(1)             
               WGH2(KR) = A0(2)             
               WGH3(KR) = A0(3)             
               WGH4(KR) = A0(4)             
            ENDIF  

            IF(TRUE_DISTL) THEN

              DOMINATE =  NCHIMT(LOCDISC(II(1))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(2))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(3))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(4))) == NCHNC
     &              .AND. NCHIMT(LOCDIS(KR))     == NCHNC

* This criterion was better in some 2D tests
*              DOMINATE = (NCHIMT(LOCDISC(II(1))) == NCHNC
*     &               .OR. NCHIMT(LOCDISC(II(2))) == NCHNC
*     &               .OR. NCHIMT(LOCDISC(II(3))) == NCHNC
*     &               .OR. NCHIMT(LOCDISC(II(4))) == NCHNC)
*     &              .AND. NCHIMT(LOCDIS(KR))     == NCHNC

            ELSE

            ALLCORNERS =  NCHIMT(LOCDISC(II(1))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(2))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(3))) == NCHNC
     &              .AND. NCHIMT(LOCDISC(II(4))) == NCHNC
              DOMINATE = (ALLCORNERS
     &              .AND. NCHIMT(LOCDIS(KR))     == NCHNC)
     &              .OR. (ALLCORNERS
     &              .AND. NCHIMT(LOCDIS(KR))     /= NCHNC
     &              .AND. DISTW(KR) > MAX(DISTWC(II(1)),DISTWC(II(2)),
     &                                    DISTWC(II(3)),DISTWC(II(4))))

            ENDIF


            CHIMGC = MOD(INTPC(II(1)),10) == 2
     &          .OR. MOD(INTPC(II(2)),10) == 2
     &          .OR. MOD(INTPC(II(3)),10) == 2
     &          .OR. MOD(INTPC(II(4)),10) == 2

            IF((DOMINATE .AND. .NOT.INSIDE 
     &                   .AND. .NOT.CHIMGC) .OR. NCHNC >= 900) THEN
               RKSI(KR)  = RKSIV
               BLANK(KR) = 0.0
            ENDIF
     
*            IF(.NOT.CHIMGC) THEN 
*                IF (NCHIMT(LOCDISC(II(1))) /= NCHN
*     &        .AND. NCHIMT(LOCDISC(II(2))) /= NCHN
*     &        .AND. NCHIMT(LOCDISC(II(3))) /= NCHN
*     &        .AND. NCHIMT(LOCDISC(II(4))) /= NCHN) BLANK(KR) = 0.0
*            ENDIF

            IF(.NOT.TRUE_DISTL) THEN
            IF(NEARWL(KP) == 1 .AND. VOL(KR)/VOLC(KP) > VRATIO) THEN 
               BLANK(KR) = -8.0  ! Strong blanking near walls
            ENDIF
            ENDIF

         ENDIF     ! KP/=0

      ENDDO ! i-loop (the background block)
      ENDDO ! j-loop (the background block)
      ENDDO ! k-loop (the background block)

      DEALLOCATE(KSC,NEARWL,IBOXES)

      RETURN
      END SUBROUTINE DONORS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ISITIN(X0,A0,K,II,NX,NY,NZ,X,Y,Z,ITISIN,VOL)

C ... If X0 lies inside cell K, then ITISIN = .TRUE.

      IMPLICIT NONE

      INTEGER          :: K, NX, NY, NZ, IT, II(8)
      REAL             :: VOL, A0(4), X0(3)
      REAL    :: X(*), Y(*), Z(*)
      LOGICAL          :: ITISIN
            
*      II(1) = K
*      II(2) = K + 1
*      II(3) = K     + NX
*      II(4) = K + 1 + NX
*      II(5) = K          + NX*NY
*      II(6) = K + 1      + NX*NY
*      II(7) = K     + NX + NX*NY
*      II(8) = K + 1 + NX + NX*NY

      II(1) = K
      II(2) = K + 1
      II(3) = K     + NX
      II(4) = II(3) + 1
      II(5) = K + NX*NY
      II(6) = II(5) + 1
      II(7) = II(5) + NX
      II(8) = II(7) + 1
      
      ITISIN = .FALSE.
      
      IF(X0(1) < MIN(X(II(1)),X(II(2)),X(II(3)),X(II(4)),
     &               X(II(5)),X(II(6)),X(II(7)),X(II(8)))) RETURN
      IF(X0(1) > MAX(X(II(1)),X(II(2)),X(II(3)),X(II(4)),
     &               X(II(5)),X(II(6)),X(II(7)),X(II(8)))) RETURN
      IF(X0(2) < MIN(Y(II(1)),Y(II(2)),Y(II(3)),Y(II(4)),
     &               Y(II(5)),Y(II(6)),Y(II(7)),Y(II(8)))) RETURN
      IF(X0(2) > MAX(Y(II(1)),Y(II(2)),Y(II(3)),Y(II(4)),
     &               Y(II(5)),Y(II(6)),Y(II(7)),Y(II(8)))) RETURN
      IF(X0(3) < MIN(Z(II(1)),Z(II(2)),Z(II(3)),Z(II(4)),
     &               Z(II(5)),Z(II(6)),Z(II(7)),Z(II(8)))) RETURN
      IF(X0(3) > MAX(Z(II(1)),Z(II(2)),Z(II(3)),Z(II(4)),
     &               Z(II(5)),Z(II(6)),Z(II(7)),Z(II(8)))) RETURN
      
      ITISIN = .TRUE.
      
*      DO IT = 1,6
      DO IT = 1,18

         CALL NODNUM(K,IT,II,NX,NY)
     
         CALL TETVOL(X0,A0,VOL,II,NX*NY*NZ,X,Y,Z)
      
         IF(VOL > 0.0) THEN
            IF(MIN(A0(1),A0(2),A0(3),A0(4)) >= 0.0) RETURN
         END IF

      ENDDO
      
      ITISIN = .FALSE.
      
      RETURN
      END SUBROUTINE ISITIN
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE NODNUM(K,IT,II,NX,NY)

C ... Calculate node numbers: 
C ... cell corners/tetrahedron corners/face corners

      IMPLICIT NONE

      INTEGER :: K, IT, NX, NY, II(8)  
*      INTEGER :: GLOBAC(8), LOCALC(4,6)
      INTEGER :: GLOBAC(8), LOCALC(4,18)
      
*      DATA LOCALC/1,2,5,3, 2,6,5,3, 3,6,5,7,
*     &            3,8,7,5, 3,4,8,5, 1,4,3,5/
      DATA LOCALC/1,2,5,3, 2,6,5,3, 3,6,5,7,
     &            3,8,7,5, 3,4,8,5, 1,4,3,5,
     &            5,6,1,7, 6,2,1,7, 7,2,1,3,
     &            7,4,3,1, 7,8,4,1, 5,8,7,1,
     &            8,5,4,6, 5,1,4,6, 6,1,4,2,
     &            6,3,2,4, 6,7,3,4, 8,7,6,4/
      
*      GLOBAC(1) = 0
*      GLOBAC(2) = 1
*      GLOBAC(3) = 1 + NX
*      GLOBAC(4) =   + NX
*      GLOBAC(5) =        + NX*NY
*      GLOBAC(6) = 1      + NX*NY
*      GLOBAC(7) = 1 + NX + NX*NY
*      GLOBAC(8) =   + NX + NX*NY

      GLOBAC(1) = 0
      GLOBAC(2) = 1
      GLOBAC(3) = NX + 1
      GLOBAC(4) = NX
      GLOBAC(5) = NX*NY
      GLOBAC(6) = GLOBAC(5) + 1
      GLOBAC(7) = GLOBAC(6) + NX
      GLOBAC(8) = GLOBAC(7) - 1
      
      II(1) = K + GLOBAC(LOCALC(1,IT))
      II(2) = K + GLOBAC(LOCALC(2,IT))
      II(3) = K + GLOBAC(LOCALC(3,IT))
      II(4) = K + GLOBAC(LOCALC(4,IT))
      
      RETURN
      END SUBROUTINE NODNUM
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE TETVOL(X0,A0,VOL,II,N,X,Y,Z)

C ... Compute tetrahedron volume and volume coordinates for point X0
C ... inside tetrahedron.      

      IMPLICIT NONE

      INTEGER          :: N, II(4)
      REAL             :: X0(3), A0(4), VOL 
      REAL    :: X(*), Y(*), Z(*)      
      REAL    :: V0, A01, A02, A03, A04, DETD3
      REAL    :: A(3), B(3), C(3), D(3)
      

      A(1) = X(II(2))-X(II(1))
      A(2) = Y(II(2))-Y(II(1))
      A(3) = Z(II(2))-Z(II(1))

      B(1) = X(II(3))-X(II(1))
      B(2) = Y(II(3))-Y(II(1))
      B(3) = Z(II(3))-Z(II(1))

      C(1) = X(II(4))-X(II(1))
      C(2) = Y(II(4))-Y(II(1))
      C(3) = Z(II(4))-Z(II(1))

      V0   = C(1)*(A(2)*B(3)-B(2)*A(3)) +
     &       C(2)*(A(3)*B(1)-B(3)*A(1)) +
     &       C(3)*(A(1)*B(2)-B(1)*A(2))
      
      IF(ABS(V0) <= 1.0E-30) THEN
         VOL = 0.0
         RETURN
      ENDIF
      
      D(1)  = X0(1)-X(II(1))
      D(2)  = X0(2)-Y(II(1))
      D(3)  = X0(3)-Z(II(1))
      
      A02   =(D(1)*(B(2)*C(3)-C(2)*B(3)) +
     &        D(2)*(B(3)*C(1)-C(3)*B(1)) +
     &        D(3)*(B(1)*C(2)-C(1)*B(2)))/V0
      A03   =(D(1)*(C(2)*A(3)-A(2)*C(3)) +
     &        D(2)*(C(3)*A(1)-A(3)*C(1)) +
     &        D(3)*(C(1)*A(2)-A(1)*C(2)))/V0
      A04   =(D(1)*(A(2)*B(3)-B(2)*A(3)) +
     &        D(2)*(A(3)*B(1)-B(3)*A(1)) +
     &        D(3)*(A(1)*B(2)-B(1)*A(2)))/V0
      
      A01   = 1.0-A02-A03-A04
      
      VOL   = ABS(V0/6.0)
      
      A0(1) = A01
      A0(2) = A02
      A0(3) = A03
      A0(4) = A04
      
      RETURN
      END SUBROUTINE TETVOL
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE BOXIT(IMAX,JMAX,KMAX,IN,JN,KN,X,Y,Z,KSC,NCB,
     &                 IBOXES,NXBOX,NYBOX,NZBOX,
     &                 XMIN,YMIN,ZMIN,DXR,DYR,DZR,MT)

C ... This subroutine shares the cells into groups in order to 
C ... accelerate the search process. If MT=1, only the space
C ... requirement is calculated. In case MT=2, the cells are
C ... actually shared.

      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, NXBOX, NYBOX, NZBOX, MT
      INTEGER :: NP, KSLABI, NCB, ICP
      INTEGER :: I, J, K, L, I1, J1, K1, L1, M
      INTEGER :: IBS, IBE, JBS, JBE, KBS, KBE
      INTEGER, DIMENSION(8) :: IB, JB, KB, N

      INTEGER, DIMENSION(*) :: KSC
      INTEGER, DIMENSION(NXBOX,NYBOX,NZBOX,2) :: IBOXES

      REAL, DIMENSION(*) :: X, Y, Z
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR

      IF(MT == 1) THEN     ! Initialize the box pointer array
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,1) = 0
                  IBOXES(I1,J1,K1,2) = 0
               ENDDO
            ENDDO
         ENDDO


         NCB = 0  ! Total number of cells in groups. 
                  ! Note that one cell can belong to 
                  ! several groups.
      ENDIF


      IF(MT == 2) THEN     ! Starting pointer for each group ('box')
         NP = 1
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,2) = NP
                  NP = NP + IBOXES(I1,J1,K1,1)
                  IBOXES(I1,J1,K1,1) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDIF


      KSLABI = (IMAX+2*IN)*(JMAX+2*JN)

*      DO K=-1,KMAX+1
*         DO J=-1,JMAX+1
*            DO I=-1,IMAX+1

      DO K=0,KMAX
         DO J=0,JMAX
            DO I=0,IMAX

               N(1)  = I+IN+(J-1+JN)*(IMAX+2*IN)+(K-1+KN)*KSLABI
               N(2)  = N(1) + 1
               N(3)  = N(1) + (IMAX+2*IN) + 1
               N(4)  = N(3) - 1
               N(5)  = N(1) + KSLABI
               N(6)  = N(2) + KSLABI 
               N(7)  = N(3) + KSLABI 
               N(8)  = N(4) + KSLABI

               DO L=1,8
                  IB(L) = INT((X(N(L))-XMIN)*DXR) + 1
                  JB(L) = INT((Y(N(L))-YMIN)*DYR) + 1
                  KB(L) = INT((Z(N(L))-ZMIN)*DZR) + 1
               ENDDO

               IBS = MAX(MINVAL(IB),1)
               IBE = MIN(MAXVAL(IB),NXBOX)
               JBS = MAX(MINVAL(JB),1)
               JBE = MIN(MAXVAL(JB),NYBOX)
               KBS = MAX(MINVAL(KB),1)
               KBE = MIN(MAXVAL(KB),NZBOX)

               DO I1=IBS,IBE
                  DO J1=JBS,JBE
                     DO K1=KBS,KBE
                        IF(MT == 1) THEN
                           IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                           NCB = NCB + 1
                        ELSE
                           ICP = IBOXES(I1,J1,K1,1)+IBOXES(I1,J1,K1,2)
                           KSC(ICP) = N(1)
                           IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE BOXIT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOXEL(NSKIN,SKINX,SKINY,SKINZ,KSC,NCB,
     &                 IBOXES,NXBOX,NYBOX,NZBOX,
     &                 XMIN,YMIN,ZMIN,DXR,DYR,DZR,MT)

C ... This subroutine shares the skin elements into groups in order to 
C ... accelerate the search process. If MT=1, only the space
C ... requirement is calculated. In case MT=2, the elements are
C ... actually shared.

      IMPLICIT NONE

      INTEGER :: NSKIN, NXBOX, NYBOX, NZBOX
      INTEGER, DIMENSION(NXBOX,NYBOX,NZBOX,2) :: IBOXES
      INTEGER, DIMENSION(*) :: KSC
      INTEGER :: MT, I1, J1, K1, N, ISKIN, NCB,ICP
      INTEGER :: IB1, IB2, IB3, IB4
      INTEGER :: JB1, JB2, JB3, JB4
      INTEGER :: KB1, KB2, KB3, KB4
      INTEGER :: IBS, IBE, JBS, JBE, KBS, KBE

      REAL, DIMENSION(4,NSKIN) :: SKINX, SKINY, SKINZ
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR

      IF(MT == 1) THEN     ! Initialize the box pointer array
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,1) = 0
                  IBOXES(I1,J1,K1,2) = 0
               ENDDO
            ENDDO
         ENDDO

         NCB = 0  ! Total number of elements in groups. Note that one 
                  ! element can belong to several groups.
      ENDIF


      IF(MT == 2) THEN     ! Starting pointer for each group ('box')
         N = 1
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,2) = N
                  N = N + IBOXES(I1,J1,K1,1)
                  IBOXES(I1,J1,K1,1) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDIF


      DO ISKIN=1,NSKIN
         
         IB1 = INT((SKINX(1,ISKIN)-XMIN)*DXR) + 1
         IB2 = INT((SKINX(2,ISKIN)-XMIN)*DXR) + 1
         IB3 = INT((SKINX(3,ISKIN)-XMIN)*DXR) + 1
         IB4 = INT((SKINX(4,ISKIN)-XMIN)*DXR) + 1
         
         JB1 = INT((SKINY(1,ISKIN)-YMIN)*DYR) + 1
         JB2 = INT((SKINY(2,ISKIN)-YMIN)*DYR) + 1
         JB3 = INT((SKINY(3,ISKIN)-YMIN)*DYR) + 1
         JB4 = INT((SKINY(4,ISKIN)-YMIN)*DYR) + 1
         
         KB1 = INT((SKINZ(1,ISKIN)-ZMIN)*DZR) + 1
         KB2 = INT((SKINZ(2,ISKIN)-ZMIN)*DZR) + 1
         KB3 = INT((SKINZ(3,ISKIN)-ZMIN)*DZR) + 1
         KB4 = INT((SKINZ(4,ISKIN)-ZMIN)*DZR) + 1
         
         IBS = MIN(IB1,IB2,IB3,IB4)
         IBE = MAX(IB1,IB2,IB3,IB4)
         JBS = MIN(JB1,JB2,JB3,JB4)
         JBE = MAX(JB1,JB2,JB3,JB4)
         KBS = MIN(KB1,KB2,KB3,KB4)
         KBE = MAX(KB1,KB2,KB3,KB4)
         IBS = MAX(IBS,1)
         IBE = MIN(IBE,NXBOX)
         JBS = MAX(JBS,1)
         JBE = MIN(JBE,NYBOX)
         KBS = MAX(KBS,1)
         KBE = MIN(KBE,NZBOX)
         
         DO I1=IBS,IBE
            DO J1=JBS,JBE
               DO K1=KBS,KBE
                  IF(MT == 1) THEN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                     NCB = NCB + 1
                  ELSE
                     ICP = IBOXES(I1,J1,K1,1)+IBOXES(I1,J1,K1,2)
                     KSC(ICP) = ISKIN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE BOXEL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUNDS(IMACHI,JMACHI,KMACHI,XCHI,YCHI,ZCHI,
     +                  XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,IN,JN,KN)

C ... This subroutine computes XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
C ... values of cell center point coordinates in a block.

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XCHI, YCHI, ZCHI
      REAL    :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, EPS
      INTEGER :: IMACHI, JMACHI, KMACHI, IN, JN, KN
      INTEGER :: ISTRID, JSTRID, KSTRID, I, J, K, KA, JJ, NODE
     

      EPS  =  1.0E-4

      XMIN =  1E30                                                      
      XMAX = -1E30                                                      
      YMIN =  1E30                                                      
      YMAX = -1E30                                                      
      ZMIN =  1E30                                                      
      ZMAX = -1E30                                                      
   
      ISTRID = IMACHI + 2*IN
      JSTRID = JMACHI + 2*JN
      KSTRID = KMACHI + 2*KN

      DO K = -1,KMACHI+2
         KA = (KN+K-1)*ISTRID*JSTRID
         DO J = -1,JMACHI+2
            JJ = (JN+J-1)*ISTRID + IN + KA 
            DO I = -1,IMACHI+2
               NODE = JJ + I 
               XMIN = MIN(XMIN,XCHI(NODE))
               XMAX = MAX(XMAX,XCHI(NODE))
               YMIN = MIN(YMIN,YCHI(NODE))
               YMAX = MAX(YMAX,YCHI(NODE))
               ZMIN = MIN(ZMIN,ZCHI(NODE))
               ZMAX = MAX(ZMAX,ZCHI(NODE))
            ENDDO
         ENDDO
      ENDDO

      IF(XMAX-XMIN < EPS) THEN
         XMIN = XMIN - EPS
         XMAX = XMAX + EPS
      ENDIF
      IF(YMAX-YMIN < EPS) THEN
         YMIN = YMIN - EPS
         YMAX = YMAX + EPS
      ENDIF
      IF(ZMAX-ZMIN < EPS) THEN
         ZMIN = ZMIN - EPS
         ZMAX = ZMAX + EPS
      ENDIF

      RETURN
      END SUBROUTINE BOUNDS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALLVI(IMAX,JMAX,KMAX,ICON,NPATCH,NEARWL,TWODIM,
     &                  IN,JN,KN)

C ... This subroutine marks cells in the vicinity of solid walls.

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*), NEARWL(*)
      INTEGER :: IN, JN, KN
      LOGICAl :: TWODIM

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN

      DO I = 1,ISTRID*JSTRID*KSTRID
         NEARWL(I) = 0
      ENDDO

      DO 1000 L = 1,NPATCH

C ... Check steady solid walls, rotating walls, and moving walls

      IF(ICON(1,L) == 8 .OR. ICON(1,L) == 9 .OR. ICON(1,L) == 10) THEN  

      IFACE = ICON(3,L)
      IA    = ICON(4,L)
      IM    = ICON(5,L)
      JA    = ICON(6,L)
      JM    = ICON(7,L)

      IF(IFACE == 1) THEN
         IFIRST = -1
         ILAST  = IMAX/2
         JFIRST = IA
         JLAST  = IM
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 2) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = -1
         JLAST  = JMAX/2
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 3) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = JA
         JLAST  = JM
         KFIRST = -1
         KLAST  = KMAX/2
      ELSEIF(IFACE == 4) THEN
         IFIRST = IMAX / 2 + 1
         ILAST  = IMAX + 2
         JFIRST = IA
         JLAST  = IM
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 5) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = JMAX / 2 + 1
         JLAST  = JMAX + 2
         KFIRST = JA
         KLAST  = JM
      ELSEIF(IFACE == 6) THEN
         IFIRST = IA
         ILAST  = IM
         JFIRST = JA
         JLAST  = JM
         KFIRST = KMAX / 2 + 1
         KLAST  = KMAX + 2
      ENDIF

      IF(TWODIM) THEN
         KFIRST = -1
         KLAST  =  3 
      ENDIF

      DO KC = KFIRST,KLAST
         KA = (KN+KC-1)*ISTRID*JSTRID
      DO JC = JFIRST,JLAST
         JJ = (JN+JC-1)*ISTRID + IN + KA 
      DO IC = IFIRST,ILAST

         KR = JJ + IC
         NEARWL(KR) = 1

      ENDDO
      ENDDO
      ENDDO

      ENDIF

 1000 CONTINUE

      END SUBROUTINE WALLVI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TBLANK(BLANK,IMAX,JMAX,KMAX,LEVEL)

C ... This subroutine manipulates the blanking vector before
C ... it will be written into XYZ.BIN.

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, IMAXP1, JMAXP1, KMAXP1
      INTEGER :: ISTRID, JSTRID, ILE, LEVEL, NEXTENT
      INTEGER :: I, J, K, L, KA, KR, JJ, NMAX
      INTEGER :: INEXT, JNEXT, KNEXT, INEXTP1, JNEXTP1, KNEXTP1
      INTEGER :: IPREV, JPREV, KPREV, IPREVM1, JPREVM1, KPREVM1

      REAL :: BLANK(*)


      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      ILE    = ISTRID*JSTRID

      NMAX = (IMAX + 2*IN)*(JMAX + 2*JN)*(KMAX + 2*KN)

C ... Trim the blanking vector, i.e. extend the visible area sligthly.
C ... The purpose of this step is to avoid ugly holes in the 
C ... post-processing.

      IF(LEVEL == 1) THEN
*         NEXTENT = 2  ! Does not improve the post-processing in all cases.
         NEXTENT = 1
      ELSE
         NEXTENT = 1
      ENDIF
      
      DO L=1,NEXTENT

C ... Increasing I-direction

      DO K  = 1,KMAX
         KA = (KN+K-1)*ILE
      DO J  = 1,JMAX
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 1,IMAX

      KR      = JJ + I

      INEXT   = KR + 1
      INEXTP1 = KR + 2
      JNEXT   = KR + ISTRID
      JNEXTP1 = KR + 2*ISTRID
      KNEXT   = KR + ILE
      KNEXTP1 = KR + 2*ILE

      IF(BLANK(INEXT) == 1.0 .AND. BLANK(INEXTP1) == 1.0 .AND.
     &   I /= IMAX .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


C ... Decreasing I-direction

      DO K  = KMAX,1,-1
         KA = (KN+K-1)*ILE
      DO J  = JMAX,1,-1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = IMAX,1,-1

      KR      = JJ + I

      IPREV   = KR - 1
      IPREVM1 = KR - 2
      JPREV   = KR - ISTRID
      JPREVM1 = KR - 2*ISTRID
      KPREV   = KR - ILE
      KPREVM1 = KR - 2*ILE

      IF(BLANK(IPREV) == 1.0 .AND. BLANK(IPREVM1) == 1.0 .AND. 
     &   I /= 1 .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


C ... Increasing J-direction

      DO K  = 1,KMAX
         KA = (KN+K-1)*ILE
      DO J  = 1,JMAX
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 1,IMAX

      KR      = JJ + I

      INEXT   = KR + 1
      INEXTP1 = KR + 2
      JNEXT   = KR + ISTRID
      JNEXTP1 = KR + 2*ISTRID
      KNEXT   = KR + ILE
      KNEXTP1 = KR + 2*ILE

      IF(BLANK(JNEXT) == 1.0 .AND. BLANK(JNEXTP1) == 1.0 .AND.
     &   J /= JMAX .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


C ... Decreasing J-direction

      DO K  = KMAX,1,-1
         KA = (KN+K-1)*ILE
      DO J  = JMAX,1,-1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = IMAX,1,-1

      KR      = JJ + I

      IPREV   = KR - 1
      IPREVM1 = KR - 2
      JPREV   = KR - ISTRID
      JPREVM1 = KR - 2*ISTRID
      KPREV   = KR - ILE
      KPREVM1 = KR - 2*ILE

      IF(BLANK(JPREV) == 1.0 .AND. BLANK(JPREVM1) == 1.0 .AND. 
     &   J /= 1 .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


C ... Increasing K-direction

      DO K  = 1,KMAX
         KA = (KN+K-1)*ILE
      DO J  = 1,JMAX
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 1,IMAX

      KR      = JJ + I

      INEXT   = KR + 1
      INEXTP1 = KR + 2
      JNEXT   = KR + ISTRID
      JNEXTP1 = KR + 2*ISTRID
      KNEXT   = KR + ILE
      KNEXTP1 = KR + 2*ILE

      IF(BLANK(KNEXT) == 1.0 .AND. BLANK(KNEXTP1) == 1.0 .AND. 
     &   K /= KMAX .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


C ... Decreasing K-direction

      DO K  = KMAX,1,-1
         KA = (KN+K-1)*ILE
      DO J  = JMAX,1,-1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = IMAX,1,-1

      KR      = JJ + I

      IPREV   = KR - 1
      IPREVM1 = KR - 2
      JPREV   = KR - ISTRID
      JPREVM1 = KR - 2*ISTRID
      KPREV   = KR - ILE
      KPREVM1 = KR - 2*ILE

      IF(BLANK(KPREV) == 1.0 .AND. BLANK(KPREVM1) == 1.0 .AND. 
     &   K /= 1 .AND. BLANK(KR) >= 0.0) BLANK(KR) = 1.0

      ENDDO
      ENDDO
      ENDDO


      ENDDO   ! L=1,2


      RETURN
      END SUBROUTINE TBLANK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE WALLGP(BLANK,IMAX,JMAX,KMAX,NPATCH,ICON)
     
C ... This subroutine marks the grid points on solid surfaces
C ... ('wall' points) in the blanking vector (iblank = 2) in order 
C ... to avoid blanking them.

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, NPATCH
      INTEGER :: IPATCH, IFACE, IA, IM, JA, JM
      INTEGER :: I, J, K, ISTRID, JSTRID, ILE, IP
      INTEGER :: IFIRST, ILAST, JFIRST, JLAST, KFIRST, KLAST

      INTEGER, DIMENSION(IC9,NPATCH) :: ICON

      REAL, DIMENSION(*) :: BLANK

      ISTRID = IMAX + 1
      JSTRID = JMAX + 1
      ILE    = ISTRID*JSTRID


      DO IPATCH=1,NPATCH

         IF(ICON(1,IPATCH) ==  1 .AND. ICON(20,IPATCH) == 8 .OR. ! ACT
     &      ICON(1,IPATCH) ==  8 .OR.     ! SOL
     &      ICON(1,IPATCH) ==  9 .OR.     ! ROT
     &      ICON(1,IPATCH) == 10) THEN    ! MOV

            IFACE = ICON(3,IPATCH)
            IA    = ICON(4,IPATCH)
            IM    = ICON(5,IPATCH) + 1
            JA    = ICON(6,IPATCH)
            JM    = ICON(7,IPATCH) + 1

            SELECT CASE(IFACE)
               CASE(1)
                  IFIRST = 1
                  ILAST  = 1
                  JFIRST = IA
                  JLAST  = IM
                  KFIRST = JA
                  KLAST  = JM
               CASE(2)
                  IFIRST = IA
                  ILAST  = IM
                  JFIRST = 1
                  JLAST  = 1
                  KFIRST = JA
                  KLAST  = JM
               CASE(3)
                  IFIRST = IA
                  ILAST  = IM
                  JFIRST = JA
                  JLAST  = JM
                  KFIRST = 1
                  KLAST  = 1
               CASE(4)
                  IFIRST = IMAX+1
                  ILAST  = IMAX+1
                  JFIRST = IA
                  JLAST  = IM
                  KFIRST = JA
                  KLAST  = JM
               CASE(5)
                  IFIRST = IA
                  ILAST  = IM
                  JFIRST = JMAX+1
                  JLAST  = JMAX+1
                  KFIRST = JA
                  KLAST  = JM
               CASE(6)
                  IFIRST = IA
                  ILAST  = IM
                  JFIRST = JA
                  JLAST  = JM
                  KFIRST = KMAX+1
                  KLAST  = KMAX+1
            END SELECT


            DO K = KFIRST,KLAST
               DO J = JFIRST,JLAST
                  DO I = IFIRST,ILAST
                     IP = I + (J-1)*ISTRID + (K-1)*ILE
                     BLANK(IP) = 2.0
                  ENDDO
               ENDDO
            ENDDO

         ENDIF   ! SOL, ROT or MOV

      ENDDO   ! Patch loop

      RETURN
      END SUBROUTINE WALLGP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IPOINT
     &(XI,YI,ZI,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,TOL,ICODE)

C ... This subroutine calculates the intersection point (XI,YI,ZI) of a
C     line defined by two points (X1,Y1,Z1) and (X2,Y2,Z2) and a plane
C     defined by three points (X3,Y3,Z3), (X4,Y4,Z4) and (X5,Y5,Z5).
C     Subroutine also checks whether the point lies within the triangle
C     or not.
C
C     The formulas for the line and the plane are:
C
C               XI - X1     YI - Y1     ZI - Z1
C               -------  =  -------  =  -------
C               X2 - X1     Y2 - Y1     Z2 - Z1
C
C              | XI - X3    YI - Y3    ZI - Z3 |
C              | X4 - X3    Y4 - Y3    Z4 - Z3 | = 0
C              | X5 - X3    Y5 - Y3    Z5 - Z3 |
C
C     Return codes:
C
C       ICODE = 0:  The line and the plane have an intersection point
C                   which lies within the triangle determined by the
C                   points (X3,Y3,Z3), (X4,Y4,Z4), and (X5,Y5,Z5).
C                   The points (X1,Y1,Z1) and (X2,Y2,Z2) are located
C                   on opposite sides of the plane.
C
C       ICODE = 1:  The line and the plane have an intersection point
C                   which lies within the triangle determined by the
C                   points (X3,Y3,Z3), (X4,Y4,Z4), and (X5,Y5,Z5).
C                   The points (X1,Y1,Z1) and (X2,Y2,Z2) are located
C                   on the same side of the plane.
C
C       ICODE = 2:  The line and the plane have an intersection point but
C                   it does not lie within the triangle determined by
C                   the points (X3,Y3,Z3), (X4,Y4,Z4), and (X5,Y5,Z5).
C
C       ICODE = 3:  Points (X1,Y1,Z1) and (X2,Y2,Z2) are one and the
C                   same point, i.e., they do not determine a line.
C
C       ICODE = 4:  At least two of the points (X3,Y3,Z3), (X4,Y4,Z4),
C                   and (X5,Y5,Z5) are one and the same point, i.e.,
C                   the points do not determine a plane.
C
C       ICODE = 5:  The line is parallel with the plane, i.e., the line
C                   does not intersect the plane or the line lies on the
C ...               plane.

      IMPLICIT REAL (A-H,O-Z)

      PARAMETER (EPS=1.0E-12)

C ... Normal of the plane
      XH = (Y4-Y3)*(Z5-Z3) - (Z4-Z3)*(Y5-Y3)   
      YH = (Z4-Z3)*(X5-X3) - (X4-X3)*(Z5-Z3)  
      ZH = (X4-X3)*(Y5-Y3) - (Y4-Y3)*(X5-X3) 

      CE = XH*(X3-X1)+YH*(Y3-Y1)+ZH*(Z3-Z1)
      CD = XH*(X2-X1)+YH*(Y2-Y1)+ZH*(Z2-Z1)

C ... Check the validity of the geometry
      IF(ABS(XH*(X2-X1) + YH*(Y2-Y1) + ZH*(Z2-Z1)) < EPS) THEN
        ICODE = 5
        D1 = SQRT((X4-X3)*(X4-X3) + (Y4-Y3)*(Y4-Y3) + (Z4-Z3)*(Z4-Z3))
        D2 = SQRT((X5-X3)*(X5-X3) + (Y5-Y3)*(Y5-Y3) + (Z5-Z3)*(Z5-Z3))
        D3 = SQRT((X5-X4)*(X5-X4) + (Y5-Y4)*(Y5-Y4) + (Z5-Z4)*(Z5-Z4))
        D4 = SQRT((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1) + (Z2-Z1)*(Z2-Z1))
        IF(D1 < EPS .OR. D2 < EPS .OR. D3 < EPS) ICODE = 4
        IF(D4 < EPS) ICODE = 3
        RETURN
      END IF

C ... Intersection point of the line and the plane
      XI = X1 + CE/CD*(X2-X1)
      YI = Y1 + CE/CD*(Y2-Y1)
      ZI = Z1 + CE/CD*(Z2-Z1)


C ... Check if the line intersects the triangle

C ... In order to avoid problems caused by a numerical inaccuracy we
C     use in the checking process an intersection point which is 
C ... temporarily moved slightly towards the centre of the triangle. 

      XG   = (X3+X4+X5)/3.0
      YG   = (Y3+Y4+Y5)/3.0
      ZG   = (Z3+Z4+Z5)/3.0
 
      XJ   = XG + (1.0-TOL/100.)*(XI-XG)
      YJ   = YG + (1.0-TOL/100.)*(YI-YG)
      ZJ   = ZG + (1.0-TOL/100.)*(ZI-ZG)

      FX34 = (Y3-YJ)*(Z4-ZJ) - (Z3-ZJ)*(Y4-YJ)
      FY34 = (Z3-ZJ)*(X4-XJ) - (X3-XJ)*(Z4-ZJ)
      FZ34 = (X3-XJ)*(Y4-YJ) - (Y3-YJ)*(X4-XJ)
      FX45 = (Y4-YJ)*(Z5-ZJ) - (Z4-ZJ)*(Y5-YJ)
      FY45 = (Z4-ZJ)*(X5-XJ) - (X4-XJ)*(Z5-ZJ)
      FZ45 = (X4-XJ)*(Y5-YJ) - (Y4-YJ)*(X5-XJ)
      FX53 = (Y5-YJ)*(Z3-ZJ) - (Z5-ZJ)*(Y3-YJ)
      FY53 = (Z5-ZJ)*(X3-XJ) - (X5-XJ)*(Z3-ZJ)
      FZ53 = (X5-XJ)*(Y3-YJ) - (Y5-YJ)*(X3-XJ)

      A    = (X2-X1)*FX34 + (Y2-Y1)*FY34 + (Z2-Z1)*FZ34 
      B    = (X2-X1)*FX45 + (Y2-Y1)*FY45 + (Z2-Z1)*FZ45 
      C    = (X2-X1)*FX53 + (Y2-Y1)*FY53 + (Z2-Z1)*FZ53

      IF(A*B < 0.0) THEN
        ICODE = 2
      ELSEIF(A*C < 0.0) THEN
        ICODE = 2
      ELSE
        ICODE = 1
      ENDIF


C ... Check the location of points (X1,Y1,Z1) and (X2,Y2,Z2)

      IF(ICODE == 1) THEN
        ICODE = 0
        XA = X1 - XI
        YA = Y1 - YI
        ZA = Z1 - ZI 
        XB = X2 - XI
        YB = Y2 - YI
        ZB = Z2 - ZI
        IF((XA*XB+YA*YB+ZA*ZB) > 0.0) ICODE = 1
      END IF

      RETURN
      END SUBROUTINE IPOINT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SKINEL(ICON,NBCS,IPRO,NPRO,NSKIN,IPSKIN)

C ... Find all surface elements. Mostly these are SOL elements but
C ... also ROT and MOV elements must be included. Some extra BCs are 
C ... needed for closing the skin or for getting the wall distance 
C ... related Chimera domination right.

      USE MPI
      USE NS3CO, ONLY : NBLOCK, IC9, PARALLEL

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SKINTY

      INTEGER :: N, IPRO, NPRO, NBCS, NSK, NSKIN, IP, IERR 
      INTEGER :: ILOW, IUPP, JLOW, JUPP
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWSKIN 
      INTEGER, DIMENSION(IC9,NBCS)        :: ICON
      INTEGER, DIMENSION(NPRO)         :: IPSKIN
      LOGICAL, DIMENSION(:), ALLOCATABLE :: SOLIDS
      
      LOGICAL :: MASTER

      ALLOCATE(IWSKIN(NPRO))  ! Work space
      ALLOCATE(SOLIDS(NBLOCK)) 

      MASTER = IPRO == 1

      SOLIDS = .FALSE.

      DO IP=1,NBCS         ! Are there solids in block N
         N = ICON(23,IP)
         IF(ICON(1,IP) ==  8 .OR.
     &      ICON(1,IP) ==  9 .OR.
     &      ICON(1,IP) == 10) SOLIDS(N) = .TRUE.
      ENDDO

      DO IP=1,NPRO
         IPSKIN(IP) = 0
         IWSKIN(IP) = 0
      ENDDO

      NSK  = 0

      DO IP=1,NBCS

         N = ICON(23,IP)  ! Local block number

C ... Select skin elements to be used in the wall distance calculation.

         IF(SKINTY(ICON,NBCS,SOLIDS,NBLOCK,N,IP)) THEN

            ILOW = ICON(4,IP)
            IUPP = ICON(5,IP)
            JLOW = ICON(6,IP)
            JUPP = ICON(7,IP)

            NSK  = NSK + (IUPP-ILOW+1)*(JUPP-JLOW+1)

         ENDIF

      ENDDO

      NSKIN   = NSK
      IWSKIN(IPRO) = NSKIN


C ... Share the information

      IF(PARALLEL) THEN

         CALL MPI_ALLREDUCE(NSK, NSKIN,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(IWSKIN,IPSKIN,NPRO,MPI_INTEGER,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)

C ... Revise the starting pointers for skin element groups and
C ... share this information with all processes.

         IF(MASTER) THEN
            DO IP=2,NPRO
               IPSKIN(IP) = IPSKIN(IP) + IPSKIN(IP-1)
            ENDDO
            DO IP=NPRO,2,-1
               IPSKIN(IP) = IPSKIN(IP-1)
            ENDDO
            IPSKIN(1) = 0
         ENDIF

         CALL MPI_BCAST(IPSKIN,NPRO,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      ENDIF

      DEALLOCATE(IWSKIN)

      RETURN
      END SUBROUTINE SKINEL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SKINXY(ICON,NBCS,NBLOCK,NBLOCG,
     &                  MGM,SKINX,SKINY,SKINZ,
     &                  SKINXN,SKINYN,SKINZN,ITSKIN,
     &                  IBSKIN,NSKIN,IPSKIN,IPRO,NPRO,M)

C ... This subroutine finds the surface element corner point coordinates
C ... and saves them in SKINX, SKINY and SKINZ. The surface element normals
C ... are saved in SKINXN, SKINYN and SKINZN. Total of nine normals will
C ... be calculated for each element, i.e. surface normal, normals
C ... at corner points (average from four surrounding elements) and normals
C ... for edges (average from two adjacent elements). 

C ... Local element (IELE) and corner point (IC) numbering.

C     IC(14) ---- IC(13) ---- IC(12) ---- IC(11)
C       |           |           |           |
C       |     8     |     7     |     6     |     
C       |           |           |           |
C     IC(15) ----  IC(4) ----  IC(3) ---- IC(10)
C       |           |           |           |
C       |     9     |     1     |     5     |     
C       |           |           |           |
C     IC(16) ----  IC(1) ----  IC(2) ----  IC(9)
C       |           |           |           |
C       |     2     |     3     |     4     |     
C       |           |           |           |
C     IC(5) ----  IC(6) ----  IC(7) ----  IC(8)


      USE MPI
      USE NS3CO, ONLY : IN, JN, KN, IC9, PARALLEL

      USE MAIN_ARRAYS, ONLY : NPATCH, IG, XCO, YCO, ZCO,
     &                        IMAX, JMAX, KMAX 

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SKINTY

      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: STATUS
      INTEGER :: NBLOCK, NBLOCG, NBCS, NSK, NSKIN
      INTEGER :: N, NG, IBL, IBG, IFACE, IP, ISTRID, JSTRID, KSTRID, ILE
      INTEGER :: IXLO, IXUP, IYLO, IYUP, IXLOT, IXUPT, IYLOT, IYUPT
      INTEGER :: IDIM, JDIM, KDIM, IPRO, NPRO, IP1
      INTEGER :: IPTR, ICORNER, INORMAL, IELE, IERR, IG1, M, MGM, IK, IA
      INTEGER :: IC(16)
      INTEGER :: IDC(4,9)

      INTEGER, DIMENSION(:), ALLOCATABLE :: IVSKIN, IWSKIN 
      INTEGER, DIMENSION(NPRO)  :: IPSKIN
      INTEGER, DIMENSION(NSKIN)   :: IBSKIN, ITSKIN
      INTEGER, DIMENSION(IC9,NBCS) :: ICON

      REAL, DIMENSION(16) :: X, Y, Z
      REAL, DIMENSION(9)  :: ANX, ANY, ANZ
      REAL :: CX, CY, CZ, DX, DY, DZ, AA

      REAL, DIMENSION(4,NSKIN) :: SKINX, SKINY, SKINZ 
      REAL, DIMENSION(9,NSKIN) :: SKINXN, SKINYN, SKINZN 
      REAL, DIMENSION(:,:), ALLOCATABLE :: 
     &      SKINXW, SKINYW, SKINZW, SKINXNW, SKINYNW, SKINZNW 

      LOGICAL, DIMENSION(:), ALLOCATABLE :: SOLIDS

      LOGICAL :: MASTER

      DATA IDC/1,2,3,4,    5,6,1,16,  6,7,2,1,    7,8,9,2,    2,9,10,3,
     &         3,10,11,12, 4,3,12,13, 15,4,13,14, 16,1,4,15/


      ALLOCATE(SKINXW (4,NSKIN))   ! Work space
      ALLOCATE(SKINYW (4,NSKIN))   ! Work space
      ALLOCATE(SKINZW (4,NSKIN))   ! Work space
      ALLOCATE(SKINXNW(9,NSKIN))   ! Work space
      ALLOCATE(SKINYNW(9,NSKIN))   ! Work space
      ALLOCATE(SKINZNW(9,NSKIN))   ! Work space
      ALLOCATE(IVSKIN(NSKIN))      ! Work space
      ALLOCATE(IWSKIN(NSKIN))      ! Work space
      ALLOCATE(SOLIDS(NBLOCK)) 

      MASTER = IPRO == 1

      SOLIDS = .FALSE.

      DO IP=1,NBCS         ! Are there solids in block N
         N = ICON(23,IP)
         IF(ICON(1,IP) ==  8 .OR.
     &      ICON(1,IP) ==  9 .OR.
     &      ICON(1,IP) == 10) SOLIDS(N) = .TRUE.
      ENDDO


      DO IP=1,NSKIN

         IBSKIN(IP) = 0
         IWSKIN(IP) = 0
         ITSKIN(IP) = 0
         IVSKIN(IP) = 0

         DO ICORNER=1,4
            SKINX (ICORNER,IP)  = 0.0
            SKINY (ICORNER,IP)  = 0.0
            SKINZ (ICORNER,IP)  = 0.0
            SKINXW(ICORNER,IP)  = 0.0
            SKINYW(ICORNER,IP)  = 0.0
            SKINZW(ICORNER,IP)  = 0.0
         ENDDO

         DO INORMAL=1,9
            SKINXN (INORMAL,IP) = 0.0
            SKINYN (INORMAL,IP) = 0.0
            SKINZN (INORMAL,IP) = 0.0
            SKINXNW(INORMAL,IP) = 0.0
            SKINYNW(INORMAL,IP) = 0.0
            SKINZNW(INORMAL,IP) = 0.0
         ENDDO

      ENDDO

      IPTR = IPSKIN(IPRO) 

      NSK  = 0 

      IP1 = 1

      DO N=1,NBLOCK

         IDIM   = IMAX(M,N)
         JDIM   = JMAX(M,N)
         KDIM   = KMAX(M,N)

         ISTRID = IDIM + 2*IN
         JSTRID = JDIM + 2*JN
         KSTRID = KDIM + 2*KN

         ILE    = ISTRID*JSTRID

         IG1    = IG(M,N)

         DO IP=IP1,IP1+NPATCH(N)-1

C ... Select skin elements to be used in the wall distance calculation.

         IF(SKINTY(ICON,NBCS,SOLIDS,NBLOCK,N,IP)) THEN

            IBL = ICON(23,IP)   ! Local block number
            IBG = ICON(24,IP)   ! Global block number

            IFACE = ICON(3,IP)

            IF((ICON(1,IP) == 3 .OR. ICON(1,IP) == 5) .AND. IFACE > 3)
     &          IFACE = IFACE - 3
            IF((ICON(1,IP) == 3 .OR. ICON(1,IP) == 5) .AND. IFACE < 4)
     &          IFACE = IFACE + 3

            IXLOT = ICON(4,IP)
            IXUPT = ICON(5,IP)
            IYLOT = ICON(6,IP)
            IYUPT = ICON(7,IP)

            DO IXLO = IXLOT,IXUPT
               IXUP = IXLO + 1
            DO IYLO = IYLOT,IYUPT
               IYUP = IYLO + 1

               SELECT CASE(IFACE)
               CASE(1)
                  IC(1)  = 1+IN + (IXLO-1+JN)*ISTRID + (IYLO-1+KN)*ILE
                  IC(2)  = 1+IN + (IXLO-1+JN)*ISTRID + (IYUP-1+KN)*ILE
                  IC(3)  = 1+IN + (IXUP-1+JN)*ISTRID + (IYUP-1+KN)*ILE
                  IC(4)  = 1+IN + (IXUP-1+JN)*ISTRID + (IYLO-1+KN)*ILE
                  IC(5)  = IC(1) - ISTRID - ILE
                  IC(6)  = IC(1) - ISTRID
                  IC(7)  = IC(2) - ISTRID
                  IC(8)  = IC(2) - ISTRID + ILE
                  IC(9)  = IC(2)          + ILE
                  IC(10) = IC(3)          + ILE
                  IC(11) = IC(3) + ISTRID + ILE
                  IC(12) = IC(3) + ISTRID
                  IC(13) = IC(4) + ISTRID
                  IC(14) = IC(4) + ISTRID - ILE
                  IC(15) = IC(4)          - ILE
                  IC(16) = IC(1)          - ILE 
               CASE(2)
                  IC(1)  = IXLO+IN + (1-1+JN)*ISTRID + (IYLO-1+KN)*ILE
                  IC(2)  = IXUP+IN + (1-1+JN)*ISTRID + (IYLO-1+KN)*ILE
                  IC(3)  = IXUP+IN + (1-1+JN)*ISTRID + (IYUP-1+KN)*ILE
                  IC(4)  = IXLO+IN + (1-1+JN)*ISTRID + (IYUP-1+KN)*ILE
                  IC(5)  = IC(1) - 1      - ILE
                  IC(6)  = IC(1)          - ILE
                  IC(7)  = IC(2)          - ILE
                  IC(8)  = IC(2) + 1      - ILE
                  IC(9)  = IC(2) + 1
                  IC(10) = IC(3) + 1
                  IC(11) = IC(3) + 1      + ILE
                  IC(12) = IC(3)          + ILE
                  IC(13) = IC(4)          + ILE
                  IC(14) = IC(4) - 1      + ILE
                  IC(15) = IC(4) - 1
                  IC(16) = IC(1) - 1
               CASE(3)
                  IC(1)  = IXLO+IN + (IYLO-1+JN)*ISTRID + (1-1+KN)*ILE
                  IC(2)  = IXLO+IN + (IYUP-1+JN)*ISTRID + (1-1+KN)*ILE
                  IC(3)  = IXUP+IN + (IYUP-1+JN)*ISTRID + (1-1+KN)*ILE
                  IC(4)  = IXUP+IN + (IYLO-1+JN)*ISTRID + (1-1+KN)*ILE
                  IC(5)  = IC(1) - 1 - ISTRID
                  IC(6)  = IC(1) - 1
                  IC(7)  = IC(2) - 1
                  IC(8)  = IC(2) - 1 + ISTRID
                  IC(9)  = IC(2)     + ISTRID       
                  IC(10) = IC(3)     + ISTRID         
                  IC(11) = IC(3) + 1 + ISTRID
                  IC(12) = IC(3) + 1
                  IC(13) = IC(4) + 1
                  IC(14) = IC(4) + 1 - ISTRID 
                  IC(15) = IC(4)     - ISTRID          
                  IC(16) = IC(1)     - ISTRID         
               CASE(4)
                  IC(1)  = IDIM+1+IN+(IXLO-1+JN)*ISTRID+(IYLO-1+KN)*ILE
                  IC(2)  = IDIM+1+IN+(IXUP-1+JN)*ISTRID+(IYLO-1+KN)*ILE
                  IC(3)  = IDIM+1+IN+(IXUP-1+JN)*ISTRID+(IYUP-1+KN)*ILE
                  IC(4)  = IDIM+1+IN+(IXLO-1+JN)*ISTRID+(IYUP-1+KN)*ILE
                  IC(5)  = IC(1) - ISTRID - ILE
                  IC(6)  = IC(1)          - ILE
                  IC(7)  = IC(2)          - ILE
                  IC(8)  = IC(2) + ISTRID - ILE
                  IC(9)  = IC(2) + ISTRID
                  IC(10) = IC(3) + ISTRID
                  IC(11) = IC(3) + ISTRID + ILE
                  IC(12) = IC(3)          + ILE
                  IC(13) = IC(4)          + ILE
                  IC(14) = IC(4) - ISTRID + ILE
                  IC(15) = IC(4) - ISTRID
                  IC(16) = IC(1) - ISTRID 
               CASE(5)
                  IC(1)  = IXLO+IN+(JDIM+1-1+JN)*ISTRID+(IYLO-1+KN)*ILE
                  IC(2)  = IXLO+IN+(JDIM+1-1+JN)*ISTRID+(IYUP-1+KN)*ILE
                  IC(3)  = IXUP+IN+(JDIM+1-1+JN)*ISTRID+(IYUP-1+KN)*ILE
                  IC(4)  = IXUP+IN+(JDIM+1-1+JN)*ISTRID+(IYLO-1+KN)*ILE
                  IC(5)  = IC(1) - 1      - ILE
                  IC(6)  = IC(1) - 1
                  IC(7)  = IC(2) - 1
                  IC(8)  = IC(2) - 1      + ILE
                  IC(9)  = IC(2)          + ILE
                  IC(10) = IC(3)          + ILE
                  IC(11) = IC(3) + 1      + ILE
                  IC(12) = IC(3) + 1
                  IC(13) = IC(4) + 1
                  IC(14) = IC(4) + 1      - ILE
                  IC(15) = IC(4)          - ILE
                  IC(16) = IC(1)          - ILE
               CASE(6)
                  IC(1)  = IXLO+IN+(IYLO-1+JN)*ISTRID+(KDIM+1-1+KN)*ILE
                  IC(2)  = IXUP+IN+(IYLO-1+JN)*ISTRID+(KDIM+1-1+KN)*ILE
                  IC(3)  = IXUP+IN+(IYUP-1+JN)*ISTRID+(KDIM+1-1+KN)*ILE
                  IC(4)  = IXLO+IN+(IYUP-1+JN)*ISTRID+(KDIM+1-1+KN)*ILE
                  IC(5)  = IC(1) - 1 - ISTRID
                  IC(6)  = IC(1)     - ISTRID 
                  IC(7)  = IC(2)     - ISTRID  
                  IC(8)  = IC(2) + 1 - ISTRID
                  IC(9)  = IC(2) + 1     
                  IC(10) = IC(3) + 1       
                  IC(11) = IC(3) + 1 + ISTRID
                  IC(12) = IC(3)     + ISTRID 
                  IC(13) = IC(4)     + ISTRID 
                  IC(14) = IC(4) - 1 + ISTRID 
                  IC(15) = IC(4) - 1         
                  IC(16) = IC(1) - 1         
               END SELECT

               NSK = NSK + 1

C ... Corner points

               IA = IG1 - 1 

               DO ICORNER=1,4
                  SKINX(ICORNER,IPTR+NSK) = XCO(IC(ICORNER)+IA)
                  SKINY(ICORNER,IPTR+NSK) = YCO(IC(ICORNER)+IA)
                  SKINZ(ICORNER,IPTR+NSK) = ZCO(IC(ICORNER)+IA)
               ENDDO                             

C ... Normals

               DO IELE=1,9

                  CX = XCO(IC(IDC(4,IELE))+IA) - XCO(IC(IDC(2,IELE))+IA)
                  CY = YCO(IC(IDC(4,IELE))+IA) - YCO(IC(IDC(2,IELE))+IA)
                  CZ = ZCO(IC(IDC(4,IELE))+IA) - ZCO(IC(IDC(2,IELE))+IA)

                  DX = XCO(IC(IDC(3,IELE))+IA) - XCO(IC(IDC(1,IELE))+IA)
                  DY = YCO(IC(IDC(3,IELE))+IA) - YCO(IC(IDC(1,IELE))+IA)
                  DZ = ZCO(IC(IDC(3,IELE))+IA) - ZCO(IC(IDC(1,IELE))+IA)

                  ANX(IELE) = CY*DZ - CZ*DY
                  ANY(IELE) = CZ*DX - CX*DZ
                  ANZ(IELE) = CX*DY - CY*DX
                  AA  = SQRT(ANX(IELE)**2 + ANY(IELE)**2 + ANZ(IELE)**2)
                  ANX(IELE) = ANX(IELE)/AA
                  ANY(IELE) = ANY(IELE)/AA
                  ANZ(IELE) = ANZ(IELE)/AA
            
               ENDDO


C ... Skin element surface normal

               SKINXN(1,IPTR+NSK) = ANX(1)
               SKINYN(1,IPTR+NSK) = ANY(1)
               SKINZN(1,IPTR+NSK) = ANZ(1)


C ... Skin element normal on edge 1-2

               CX = ANX(1) + ANX(3)
               CY = ANY(1) + ANY(3)
               CZ = ANZ(1) + ANZ(3)
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(6,IPTR+NSK) = CX/AA
               SKINYN(6,IPTR+NSK) = CY/AA 
               SKINZN(6,IPTR+NSK) = CZ/AA 


C ... Skin element normal on edge 2-3

               CX = ANX(1) + ANX(5)
               CY = ANY(1) + ANY(5)
               CZ = ANZ(1) + ANZ(5)
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(7,IPTR+NSK) = CX/AA
               SKINYN(7,IPTR+NSK) = CY/AA 
               SKINZN(7,IPTR+NSK) = CZ/AA 


C ... Skin element normal on edge 3-4

               CX = ANX(1) + ANX(7)
               CY = ANY(1) + ANY(7)
               CZ = ANZ(1) + ANZ(7)
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(8,IPTR+NSK) = CX/AA
               SKINYN(8,IPTR+NSK) = CY/AA 
               SKINZN(8,IPTR+NSK) = CZ/AA 


C ... Skin element normal on edge 4-1

               CX = ANX(1) + ANX(9)
               CY = ANY(1) + ANY(9)
               CZ = ANZ(1) + ANZ(9)
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(9,IPTR+NSK) = CX/AA
               SKINYN(9,IPTR+NSK) = CY/AA 
               SKINZN(9,IPTR+NSK) = CZ/AA 


C ... Surface normals at element corners are currently based on adjacent 
C ... edge normals. Originally, the corner normals were calculated based on 
C ... surrounding four element surface normals. However, this way leads to 
C ... erroneous normals at four solid patch corners. This is because of an
C ... inexact ghost cell geometry at corners.

C ... Skin element normal at corner 1

               CX = SKINXN(9,IPTR+NSK) + SKINXN(6,IPTR+NSK)
               CY = SKINYN(9,IPTR+NSK) + SKINYN(6,IPTR+NSK) 
               CZ = SKINZN(9,IPTR+NSK) + SKINZN(6,IPTR+NSK) 
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(2,IPTR+NSK) = CX/AA
               SKINYN(2,IPTR+NSK) = CY/AA 
               SKINZN(2,IPTR+NSK) = CZ/AA 


C ... Skin element normal at corner 2

               CX = SKINXN(6,IPTR+NSK) + SKINXN(7,IPTR+NSK)
               CY = SKINYN(6,IPTR+NSK) + SKINYN(7,IPTR+NSK) 
               CZ = SKINZN(6,IPTR+NSK) + SKINZN(7,IPTR+NSK) 
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(3,IPTR+NSK) = CX/AA
               SKINYN(3,IPTR+NSK) = CY/AA 
               SKINZN(3,IPTR+NSK) = CZ/AA 


C ... Skin element normal at corner 3

               CX = SKINXN(7,IPTR+NSK) + SKINXN(8,IPTR+NSK)
               CY = SKINYN(7,IPTR+NSK) + SKINYN(8,IPTR+NSK) 
               CZ = SKINZN(7,IPTR+NSK) + SKINZN(8,IPTR+NSK) 
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(4,IPTR+NSK) = CX/AA
               SKINYN(4,IPTR+NSK) = CY/AA 
               SKINZN(4,IPTR+NSK) = CZ/AA 


C ... Skin element normal at corner 4

               CX = SKINXN(8,IPTR+NSK) + SKINXN(9,IPTR+NSK)
               CY = SKINYN(8,IPTR+NSK) + SKINYN(9,IPTR+NSK) 
               CZ = SKINZN(8,IPTR+NSK) + SKINZN(9,IPTR+NSK) 
               AA = SQRT(CX**2 + CY**2 + CZ**2)
               SKINXN(5,IPTR+NSK) = CX/AA
               SKINYN(5,IPTR+NSK) = CY/AA 
               SKINZN(5,IPTR+NSK) = CZ/AA 


C ... Skin element belongs globally to block ...

               IBSKIN(IPTR+NSK)  = IBG 

C ... Skin element type

               ITSKIN(IPTR+NSK)  = ICON(1,IP)
               IF(ICON(1,IP) == 1 .AND. ICON(20,IP) == 7)  ! CML 
     &            ITSKIN(IPTR+NSK) = 18
               IF(ICON(1,IP) == 1 .AND. ICON(20,IP) == 8)  ! ACT
     &            ITSKIN(IPTR+NSK) = 19


            ENDDO
            ENDDO

 
            ENDIF  ! SKINTY
         
         ENDDO  ! Patch loop

         IP1 = IP1 + NPATCH(N)
 
      ENDDO  ! Block loop


C ... The usage of the MPI_REDUCE command to combine the
C ... skin point vectors from different processes is clumsy.
C ... The MPI_GATHERV would probably do the same thing more
C ... efficiently. Especially, when the number of processes is large.

      IF(PARALLEL) THEN

         CALL MPI_REDUCE(IBSKIN,IWSKIN,NSKIN,MPI_INTEGER,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(ITSKIN,IVSKIN,NSKIN,MPI_INTEGER,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)

C ... Skin element cornerpoints

         CALL MPI_REDUCE(SKINX,SKINXW,4*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SKINY,SKINYW,4*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SKINZ,SKINZW,4*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)

         IF(MASTER) THEN
            IK = 0
            DO NG = 1,NBLOCG
               DO IP=1,NSKIN
                  IF(IWSKIN(IP) == NG) THEN
                     IK = IK + 1
                     DO ICORNER=1,4
                        SKINX(ICORNER,IK) = SKINXW(ICORNER,IP)
                        SKINY(ICORNER,IK) = SKINYW(ICORNER,IP)
                        SKINZ(ICORNER,IK) = SKINZW(ICORNER,IP)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         CALL MPI_BCAST(SKINX,4*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SKINY,4*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SKINZ,4*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)

C ... Skin element normals (surface, corners, edges)

         CALL MPI_REDUCE(SKINXN,SKINXNW,9*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SKINYN,SKINYNW,9*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SKINZN,SKINZNW,9*NSKIN,MPI_REAL8,
     &        MPI_SUM,0,MPI_COMM_WORLD,IERR)

         IF(MASTER) THEN
            IK = 0
            DO NG = 1,NBLOCG
               DO IP=1,NSKIN
                  IF(IWSKIN(IP) == NG) THEN
                     IK = IK + 1
                     DO INORMAL=1,9
                        SKINXN(INORMAL,IK) = SKINXNW(INORMAL,IP)
                        SKINYN(INORMAL,IK) = SKINYNW(INORMAL,IP)
                        SKINZN(INORMAL,IK) = SKINZNW(INORMAL,IP)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         CALL MPI_BCAST(SKINXN,9*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SKINYN,9*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SKINZN,9*NSKIN,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)

C ... Skin element belongs to block ...

         IF(MASTER) THEN
            IK = 0
            DO NG = 1,NBLOCG
               DO IP=1,NSKIN
                  IF(IWSKIN(IP) == NG) THEN
                     IK = IK + 1
                     IBSKIN(IK) = IWSKIN(IP)
                     ITSKIN(IK) = IVSKIN(IP)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         CALL MPI_BCAST(IBSKIN,NSKIN,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITSKIN,NSKIN,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      ENDIF

      DEALLOCATE(IWSKIN,IVSKIN)
      DEALLOCATE(SKINXW,SKINYW,SKINZW,SKINXNW,SKINYNW,SKINZNW)

      RETURN
      END SUBROUTINE SKINXY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C ... Select solid type skin elements for wall distance calculations.
C ... Basically, you should add only solid BC types here. However, some
C ... extra BCs are needed for closing the skin or for getting the wall
C ... distance related Chimera domination right.

      LOGICAL FUNCTION SKINTY(ICON,NBCS,SOLIDS,NBLOCK,N,IP)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: NBCS, NBLOCK, N, IP
      INTEGER, DIMENSION(IC9,NBCS) :: ICON
      LOGICAL, DIMENSION(NBLOCK)  :: SOLIDS

      IF(ICON(1,IP) ==  8  .OR.                ! SOL
     &   ICON(1,IP) == 10  .OR.                ! MOV
     &   ICON(1,IP) ==  3  .AND. SOLIDS(N) .OR.! INL
     &   ICON(1,IP) ==  5  .AND. SOLIDS(N) .OR.! OUT
c     &   ICON(1,IP) == 13  .OR.                ! FRE
     &   ICON(1,IP) == 15  .OR.                ! HTS
c     &   ICON(1,IP) == 16  .OR.               ! SLI
     &   ICON(1,IP) ==  9  .OR.                ! ROT
     &   ICON(1,IP) ==  1  .AND.
     &   ICON(20,IP) == 7 .OR.                                  ! CML
     &   ICON(1,IP) ==  1  .AND.
     &   ICON(20,IP) == 8) THEN                                 ! ACT
         SKINTY = .TRUE.
      ELSE
         SKINTY = .FALSE.
      ENDIF 

      RETURN
      END FUNCTION SKINTY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FIXKSI(RKSI,WGH1,IMAX,JMAX,KMAX,IN,JN,KN,RKSIV,NGL,
     + IPRO)

C ... This subroutine updates the RKSI values (fortify pointers)
C ... in isolated cells.  

      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, ISTRID, JSTRID, ILE, I, J, K, IP
      INTEGER :: IN, JN, KN, KA, JJ, INTPI, INTPB
      INTEGER :: KSI0, KSI1, NGL, IPRO
 
      REAL :: RKSI(*),WGH1(*)
      REAL :: RKSIV

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      ILE    = ISTRID*JSTRID

      IF(KMAX == 1) THEN ! Two-dimensional case

      DO K  = 1,1
         KA = (KN+K-1)*ISTRID*JSTRID
      DO J  = 0,JMAX+1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 0,IMAX+1

         IP = JJ + I

         IF(RKSI(IP-1)           == 0.0 .AND. RKSI(IP+1)      == 0.0 
     &   .AND. RKSI(IP-ISTRID)   == 0.0 .AND. RKSI(IP+ISTRID) == 0.0 
     &   .AND. RKSI(IP)  > 0.0)                       RKSI(IP) = 0.0

         IF(RKSI(IP-1)            > 0.0 .AND. RKSI(IP+1)       > 0.0 
     &   .AND. RKSI(IP-ISTRID)    > 0.0 .AND. RKSI(IP+ISTRID)  > 0.0 
     &   .AND. RKSI(IP) == 0.0 .AND. WGH1(IP) > 0.0) RKSI(IP) = RKSIV

      ENDDO
      ENDDO
      ENDDO

      ELSE

      DO K  = 0,KMAX+1
         KA = (KN+K-1)*ISTRID*JSTRID
      DO J  = 0,JMAX+1
         JJ = (JN+J-1)*ISTRID + IN + KA 
      DO I  = 0,IMAX+1

         IP = JJ + I
         IF(RKSI(IP-1)            == 0.0 .AND. RKSI(IP+1)      == 0.0 
     &   .AND. RKSI(IP-ISTRID)    == 0.0 .AND. RKSI(IP+ISTRID) == 0.0 
     &   .AND. RKSI(IP-ILE)       == 0.0 .AND. RKSI(IP+ILE)    == 0.0 
     &   .AND. RKSI(IP)  > 0.0)                       RKSI(IP)  = 0.0

         IF(RKSI(IP-1)             > 0.0 .AND. RKSI(IP+1)      > 0.0 
     &   .AND. RKSI(IP-ISTRID)     > 0.0 .AND. RKSI(IP+ISTRID) > 0.0 
     &   .AND. RKSI(IP-ILE)        > 0.0 .AND. RKSI(IP+ILE)    > 0.0 
     &   .AND. RKSI(IP) == 0.0 .AND. WGH1(IP) > 0.0) RKSI(IP) = RKSIV

      ENDDO
      ENDDO
      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE FIXKSI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ZEROVELO(U,V,W,XC,YC,ZC,INTP,BLANK,OMEGA,A,B,C,
     + CENX,CENY,CENZ,IMAX,JMAX,KMAX,IN,JN,KN,NGL)

C ... Set the flow field slowly to rest or to grid speed in cells 
C ... that lie inside the structure.

      IMPLICIT NONE

      REAL, DIMENSION(*) :: U, V, W, BLANK
      INTEGER, DIMENSION(*) :: INTP
      REAL, DIMENSION(*) :: XC, YC, ZC
      REAL    :: OMEGA, A, B, C, CENX, CENY, CENZ, DX, DY, DZ
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, NGL
      INTEGER :: ISTRID, JSTRID, KSTRID, I, J, K, II, JJ, KK, INTPII

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = ISTRID*JSTRID

*      DO K = 0,KMAX+1
      DO K = -1,KMAX+2
         KK = (KN+K-1)*KSTRID
*         DO J = 0,JMAX+1
         DO J = -1,JMAX+2
            JJ = (JN+J-1)*ISTRID + IN + KK
*            DO I = 0,IMAX+1
            DO I = -1,IMAX+2

               II = I + JJ

               INTPII = MOD(INTP(II),10)
cc              if(ngl == 69) write(666,*)  i,j,k,intpii,blank(ii)

*               IF(INTPII==3 .OR. INTPII==4 .OR. INTPII==5) THEN
               IF(ABS(BLANK(II)+8.0) < 1.E-3) THEN

                  DX = XC(II)-CENX
                  DY = YC(II)-CENY
                  DZ = ZC(II)-CENZ

                  U(II)=OMEGA*(B*DZ-C*DY)+0.95*(U(II)-OMEGA*(B*DZ-C*DY))
                  V(II)=OMEGA*(C*DX-A*DZ)+0.95*(V(II)-OMEGA*(C*DX-A*DZ))
                  W(II)=OMEGA*(A*DY-B*DX)+0.95*(W(II)-OMEGA*(A*DY-B*DX))

c                  BLANK(II) = 0.0

               ENDIF

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ZEROVELO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ZEROVAR(U,V,W,XC,YC,ZC,INTP,BLANK,OMEGA,A,B,C,
     + CENX,CENY,CENZ,VIST,EPS2,RO,P,PDIFF,TEMP,RK,REPS,FRSDEN,FRSPRE,
     + FRSSIE,RKLIM,EPSLIM,FRSMUT,FRSVIS,FRSTEM,
     + IMAX,JMAX,KMAX,IN,JN,KN,NGL)

C ... Set all flow-field values slowly to background, rest or grid-speed 
C ... values in cells that lie inside the structure.

      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, NGL
      INTEGER :: ISTRID, JSTRID, KSTRID, I, J, K, II, JJ, KK, INTPII
      INTEGER, DIMENSION(*) :: INTP
      REAL, DIMENSION(*) :: U, V, W, BLANK, VIST,
     +            EPS2, RO, P, PDIFF,TEMP, RK, REPS
     +           
      REAL, DIMENSION(*) :: XC, YC, ZC
      REAL    :: OMEGA, A, B, C, CENX, CENY, CENZ, DX, DY, DZ,
     +           UII, VII, WII, FRSDEN, FRSPRE, FRSSIE, RKLIM,
     +           EPSLIM, FRSMUT, FRSVIS, EPS2II, ROII, FRSTEM

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = ISTRID*JSTRID

*      DO K = 1,KMAX
      DO K = -1,KMAX+2
         KK = (KN+K-1)*KSTRID
*         DO J = 1,JMAX
         DO J = -1,JMAX+2
            JJ = (JN+J-1)*ISTRID + IN + KK
*            DO I = 1,IMAX
            DO I = -1,IMAX+2

               II = I + JJ

               INTPII = MOD(INTP(II),10)
cc              if(ngl == 69) write(666,*)  i,j,k,intpii,blank(ii)

*               IF(INTPII==3 .OR. INTPII==4 .OR. INTPII==5) THEN
C ... BLANK2 is hopefully temporary 

               IF(ABS(BLANK(II)+8.0) < 1.E-3) THEN

                  DX = XC(II)-CENX
                  DY = YC(II)-CENY
                  DZ = ZC(II)-CENZ

                  UII=OMEGA*(B*DZ-C*DY)+0.95*(U(II)-OMEGA*(B*DZ-C*DY))
                  VII=OMEGA*(C*DX-A*DZ)+0.95*(V(II)-OMEGA*(C*DX-A*DZ))
                  WII=OMEGA*(A*DY-B*DX)+0.95*(W(II)-OMEGA*(A*DY-B*DX))
                  ROII      = FRSDEN+0.95*(RO(II)-FRSDEN)

C ... Ideas above were left behind (mersu), temperature was added

                  P(II)     = FRSPRE + 0.95*(P(II)-FRSPRE)
                  PDIFF(II) = P(II) - FRSPRE

                  TEMP(II)  = FRSTEM + 0.95*(TEMP(II) - FRSTEM)

                  VIST(II)  = FRSMUT+ 0.95*(VIST(II)-FRSMUT)
                  EPS2(II)  = 1. + FRSMUT/FRSVIS
c                  write(1000+ngl,*) P(II),TEMP(II),EPS2(II)
c                  EPS2(II)  = EPS2II+ 0.95*(EPS2(II)-EPS2II)

               ENDIF

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ZEROVAR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      LOGICAL FUNCTION NOBBOXHIT(NG,NCG)

C ... Check if the bounding boxes of two blocks hit or not.

      USE MAIN_ARRAYS, ONLY : BXMIN, BXMAX, BYMIN, BYMAX, BZMIN, BZMAX

      IMPLICIT NONE

      INTEGER :: NG, NCG

      NOBBOXHIT = .FALSE.

      IF(BXMAX(NG) < BXMIN(NCG)) NOBBOXHIT = .TRUE.
      IF(BXMIN(NG) > BXMAX(NCG)) NOBBOXHIT = .TRUE.

      IF(BYMAX(NG) < BYMIN(NCG)) NOBBOXHIT = .TRUE.
      IF(BYMIN(NG) > BYMAX(NCG)) NOBBOXHIT = .TRUE.

      IF(BZMAX(NG) < BZMIN(NCG)) NOBBOXHIT = .TRUE.
      IF(BZMIN(NG) > BZMAX(NCG)) NOBBOXHIT = .TRUE.

      RETURN
      END FUNCTION NOBBOXHIT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

