C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ITAGMO(ITAG,ICON,NPATCH,NPPV,N)

C ... THIS SUBROUTINE MAKES ITAG SO THAT COMMUCATION WORKS FOR MOST
C ... CASES. MULTIGRID AND ONLY THE PARTLY CONNECTED GROUP OF PATCHES
C ... WILL MAKE DIFFICULTIES. SO ALL BLOCKS MUST BE CONNECTED SAME TIME!!
C ... DO NOT APPLY TO THE BC-BLOCKS IF THEY ARE IN SAME PROCESS THEY
C ... ARE CONNECTED. PPR 28.5.96

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*), ITAG(NPPV,*)

      DO 1000 L = 1,NPATCH
         IVOL  = ICON(8,L)
         LP    = ICON(15,L)
         IF(ICON(18,L) >= 2 ) ITAG(L,N) = 1
 1000 CONTINUE

      RETURN
      END SUBROUTINE ITAGMO

C     ****************************************************
C     * MPI handing subroutines                          *
C     ****************************************************

      SUBROUTINE MPRO1(IPRO,PRN,NPRO,IOTY)

C ... READ MULTIPROCESSOR MANAGER FILE IN AND MAKE SOME CHECKS

      USE MPI

      USE NS3CO, ONLY : PARALLEL

      IMPLICIT NONE

      CHARACTER(LEN=3)  :: PRN
      CHARACTER(LEN=80) :: EPSF
      LOGICAL :: THERE, THER2
      INTEGER :: I0, I, IPRO, NPRO, IOTY, MYID, IERR, NUMPROCS, IARGC
      INTEGER :: ERRORCODE
      
C     Pick the file name from command-line, if one is available

      IOTY = 0

      I0 = IARGC()
      IF (I0 > 0) THEN
         DO I=1,I0
            CALL GETARG(I,EPSF)
            IF (EPSF(1:2) == '-p' .or. EPSF(1:2) == '-a') THEN ! parallel run
               IOTY = 4
               GOTO 111
            ENDIF
         ENDDO

 111     CONTINUE
      ENDIF
       
      IF(IOTY == 0) THEN  ! A single or multiprocessor run without arguments
         THERE = .FALSE.
c        INQUIRE(FILE='PROCES',EXIST=THERE) ! Purpose unknown, commented
         IF (THERE) THEN
            OPEN(68,FILE='PROCES',FORM='FORMATTED')
            READ(68,*) NPRO,IOTY
            IF(IOTY /= 4 .AND. NPRO == 1) PRN = '   '
            CLOSE(68)
         ELSE             ! A single-processor run if the file does not exist
            IPRO   = 1
            PRN    = '   '
            NPRO = 1
            IOTY   = 1
         ENDIF
      ENDIF
          
      PARALLEL = .FALSE.
      IF(IOTY == 3 .OR. IOTY == 4) THEN
         PARALLEL = .TRUE.
         CALL MPI_INIT(IERR)
         CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
         CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
         IF(IOTY == 4) NPRO = NUMPROCS
         IF(IOTY == 4) IOTY = 3
         IF(NPRO /= NUMPROCS) THEN
            IF(myid == 0) WRITE(*,313) numprocs, NPRO
            CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
         ENDIF
 313     FORMAT('Number of prosessors in MPIRUN command (',I3,')'/
     +        'do not match the number in file PROCES (',I3,').'/
     +        'Exiting ...')
         IPRO = MYID + 1
         CALL NUMCH3(PRN,IPRO)
      ENDIF
      GOTO 200
        
 100  WRITE(*,*) 'Error in reading file PROCES'
      WRITE(*,*) 'Exiting ...'
      STOP

 200  CONTINUE

      RETURN
      END SUBROUTINE MPRO1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MPRO2(NPROCE,MBPRO)

C ... READ MULTIPROCESSOR MANAGER FILE IN AND MAKE SOME CHECKS

      IMPLICIT NONE

      INTEGER :: MBPRO,NBLOCK,NPRO,IOTY,N,J,NP,INUM,
     +           IERR,MPI_COMM_WORLD
      INTEGER :: NPROCE(MBPRO+1,*)

      REAL :: ERRORCODE

      LOGICAL :: THER2

C ... MULTIPROCESSOR RUN

      NBLOCK = 0
      OPEN(68,FILE='PROCES',FORM='FORMATTED')
      READ(68,*) NPRO,IOTY

      IF(IOTY /= 4) THEN
      DO 10 N = 1,NPRO
         READ(68,*,ERR = 100) NPROCE(1,N)
         READ(68,*,ERR = 100) (NPROCE(J,N),J=2,NPROCE(1,N)+1)
         NBLOCK = NBLOCK + NPROCE(1,N)
 10   CONTINUE
 
      DO 20 N = 1,NPRO
         DO 20 J = 3,NPROCE(1,N)+1
            IF (NPROCE(J-1,N) >= NPROCE(J,N)) THEN
               WRITE(*,*) 'Blocks are not in right order in file PROCES'
               WRITE(*,*) N,J,NPROCE(J-1,N), NPROCE(J,N)
               WRITE(*,*) 'Exiting ...'
               STOP
            ENDIF
 20   CONTINUE

      CLOSE(68)
      THER2 = .FALSE.

      DO 30 N = 1,NBLOCK
         INUM = 0
         DO NP = 1,NPRO
            DO J = 1,NPROCE(1,NP)
               IF(N == NPROCE(J+1,NP)) INUM = INUM + 1
            ENDDO
         ENDDO
         IF(INUM /= 1) WRITE(*,*) 'Block',N,' is defined',INUM,
     +        ' times in PROCES-file'
         IF(INUM /= 1) THER2 = .TRUE.
30    CONTINUE

      IF(NBLOCK /= MBPRO) THEN
         WRITE(*,*)'Different number of blocks in INPUT and PROCES file'
     +        ,NBLOCK,MBPRO
         WRITE(*,*) 'Exiting'
         CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF
         
      IF(THER2) STOP
      ELSE
         IOTY = 3
      ENDIF
      GOTO 200

 100  WRITE(*,*) 'Error in reading file PROCES'
      WRITE(*,*) 'Exiting ...'
      STOP

 200  CONTINUE

      RETURN
      END SUBROUTINE MPRO2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MKPROC(NPROCE,MBPRO,NPRO,ICON,NBCS,NBLOCK,IX,IY,IZ,
     &            INPARA,INPAR,LEVEL)

C ... This subroutine shares the blocks for parallel processing.

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: NPROCE(MBPRO+1,*), ICON(IC9,*), INPARA(INPAR,*)
      INTEGER   :: IX(*), IY(*), IZ(*)
      LOGICAL   :: FORWARD, TURN, SWAP, SORTED

      INTEGER, ALLOCATABLE :: NBCSF(:,:)
      INTEGER, ALLOCATABLE :: NOFBPP(:)     
      INTEGER, ALLOCATABLE :: JAKO(:)     
      INTEGER, ALLOCATABLE :: JAKOPTR(:)     
      INTEGER, ALLOCATABLE :: IPFIRST(:)     
      INTEGER, ALLOCATABLE :: IPLAST(:)     
      INTEGER, ALLOCATABLE :: NCELL(:,:)     
      INTEGER, ALLOCATABLE :: MCELL(:,:)     
      INTEGER, ALLOCATABLE :: NCELLZ(:)     
      INTEGER, ALLOCATABLE :: NCELLT(:)     

      ALLOCATE(NBCSF(6,NBLOCK+1))
      ALLOCATE(NOFBPP(NBLOCK+1))
      ALLOCATE(JAKO(NBLOCK+1))
      ALLOCATE(JAKOPTR(NBLOCK+1))
      ALLOCATE(IPFIRST(NBLOCK+1))
      ALLOCATE(IPLAST(NBLOCK+1))
      ALLOCATE(NCELL(3,NBLOCK+1))
      ALLOCATE(MCELL(3,NBLOCK+1))
      ALLOCATE(NCELLZ(0:NBLOCK+1))
      ALLOCATE(NCELLT(NBLOCK+1))

      OPEN(68,FILE='PROCES',FORM='FORMATTED',STATUS='UNKNOWN')

      NC   = NBLOCK
      NCPU = NPRO

      WRITE(45,*)
      WRITE(45,*) 'Allocating memory for NBCSF etc. in MKPROC'

C ......................................................................

C ... Calculate the number of INLET/OUTLET blocks

      NIN = 0


C ... Calculate the block sizes

      DO IB = 1,NC

         NCELL(1,IB) = IB
         MCELL(1,IB) = IB
         NCELL(2,IB) = 0
         MCELL(2,IB) = 0

         DO IM = LEVEL,LEVEL+INPARA(16,IB)-1

            IDIV = 2**(IM-1)
            IXLE = MAX0((IX(IB)-1)/IDIV,2)
            IYLE = MAX0((IY(IB)-1)/IDIV,2)
            IZLE = MAX0((IZ(IB)-1)/IDIV,2)
            NELE =  IXLE   * IYLE   * IZLE     ! Ghost cells excluded
            NILE = (IXLE+2)*(IYLE+2)*(IZLE+2)  ! First ghost cells included 
            MCELL(2,IB) = MCELL(2,IB) + NELE
            NCELL(2,IB) = NCELL(2,IB) + NILE 

         ENDDO

      ENDDO

C ... Calculate the number of BC patches on each block face

      DO IB=1,NC
         DO IFACE=1,6
            K = 0
            DO IP=1,NBCS
               IF(ICON(2,IP) == IB .AND. ICON(3,IP) == IFACE) THEN 
                  K = K + 1
               ENDIF  
            ENDDO
            NBCSF(IFACE,IB) = K
         ENDDO
      ENDDO

C ... BC patch number limits in each block

      IPFIRST(1) = 1
      IPLAST(1)  = 0
      DO IFACE=1,6
         IPLAST(1)  = IPLAST(1) + NBCSF(IFACE,1)
      ENDDO
      DO IB=2,NC
         IP1 = 0
         IP2 = 0
         DO IFACE=1,6
            IP1 = IP1 + NBCSF(IFACE,IB-1)
            IP2 = IP2 + NBCSF(IFACE,IB)
         ENDDO
         IPFIRST(IB) = IPFIRST(IB-1) + IP1
         IPLAST(IB)  = IPLAST(IB-1)  + IP2
      ENDDO

C ... Sort the blocks according to the sizes
      
      I = 1
 641  CONTINUE
        I = I + 1
        SORTED = .TRUE.
        DO J=NC-NIN,I,-1
          IF(NCELL(2,J-1) < NCELL(2,J)) THEN
            SORTED       = .FALSE.
            KAPU         = NCELL(1,J-1)
            LAPU         = NCELL(2,J-1)
            NCELL(1,J-1) = NCELL(1,J)
            NCELL(2,J-1) = NCELL(2,J)
            NCELL(1,J)   = KAPU
            NCELL(2,J)   = LAPU
            KAPU         = MCELL(1,J-1)
            LAPU         = MCELL(2,J-1)
            MCELL(1,J-1) = MCELL(1,J)
            MCELL(2,J-1) = MCELL(2,J)
            MCELL(1,J)   = KAPU
            MCELL(2,J)   = LAPU
          ENDIF
        ENDDO
        IF(SORTED) GOTO 642
      GOTO 641
      
 642  CONTINUE

      DO J = 1,NC
         MCELL(3,J) = MCELL(1,J)
         NCELL(3,J) = NCELL(1,J)
      ENDDO

 742  CONTINUE
      DO J = 1,NC
         MCELL(1,J) = MCELL(3,J)
         NCELL(1,J) = NCELL(3,J)
      ENDDO


      N1 = 0
      DO IB = 1,NC-NIN
        DO IFACE = 1,6
           N1 = N1 + NBCSF(IFACE,IB)
        ENDDO
      ENDDO

      N2 = 0
      DO IB = 1,NC
        DO IFACE = 1,6
           N2 = N2 + NBCSF(IFACE,IB)
        ENDDO
      ENDDO


C ... Share the blocks to processes

      DO I = 1,NBLOCK
         NOFBPP(I) = 0
      ENDDO

      DO I=1,NCPU
         JAKOPTR(I) = 1
         NCELLZ(I)  = 0
         NCELLT(I)  = 0
      ENDDO

      IB      = 0
      I       = 1
      IMAX    = 1
      FORWARD = .FALSE.

 841  CONTINUE

      TURN = .FALSE.
      
      IF((I == 1 .AND. .NOT.FORWARD) .OR.
     &   (I > 1 .AND. .NOT.FORWARD .AND. 
     &    NCELLZ(I-1) > NCELLZ(I))) THEN
         TURN    = .TRUE.
         FORWARD = .TRUE.
      ELSEIF((I == NCPU .AND. FORWARD) .OR.
     &       (I < NCPU .AND. FORWARD .AND. 
     &        NCELLZ(I+1) > NCELLZ(I))) THEN
         TURN    = .TRUE.
         FORWARD = .FALSE.
      ENDIF

      IF(TURN) THEN
         I = I
      ELSE
         IF(FORWARD) THEN
            I = I + 1
            IF(I > IMAX) IMAX = I
         ELSE
            I = I - 1
         ENDIF
      ENDIF

      IF(IB < NC-NIN) THEN
         IB = IB + 1
         NOFBPP(I) = NOFBPP(I) + 1

C ... Choose from similar blocks that one, which has connections to
C ... blocks already in this process.

         IBSIZE    = NCELL(2,IB)
         ID        = IB - 1
         IS        = IB
         PRIORITY  = 0.0
         SWAP      = .FALSE.

 123     CONTINUE

         ID        = ID + 1
         IF(NCELL(2,ID) /= IBSIZE .OR. ID > NC-NIN) GOTO 124
         IF(NOFBPP(I) > 1) THEN
            IPSIZE = 0
            DO IC = JAKOPTR(I)-NOFBPP(I)+1,JAKOPTR(I)-1
               DO IP = IPFIRST(NCELL(1,ID)),IPLAST(NCELL(1,ID))
                  IF(ICON(1,IP) == 1 .AND. 
     &               ICON(8,IP) == JAKO(IC)) THEN
                     IPSIZE = IPSIZE + (ICON(5,IP) - ICON(4,IP) + 1)
     &                               * (ICON(7,IP) - ICON(6,IP) + 1)
     &                               + 250 ! Overhead due to MPI
                     SWAP = .TRUE.
                  ENDIF
               ENDDO
            ENDDO
            IF(IPSIZE > PRIORITY) THEN
               PRIORITY = IPSIZE
               IS = ID
            ENDIF
         ENDIF
         GOTO 123
 124     CONTINUE

         IF(SWAP) THEN
            ID          = IS
            IAPU        = NCELL(1,ID)
            JAPU        = NCELL(2,ID)
            KAPU        = MCELL(1,ID)
            LAPU        = MCELL(2,ID)
            NCELL(1,ID) = NCELL(1,IB)
            NCELL(2,ID) = NCELL(2,IB)
            MCELL(1,ID) = MCELL(1,IB)
            MCELL(2,ID) = MCELL(2,IB)
            NCELL(1,IB) = IAPU
            NCELL(2,IB) = JAPU
            MCELL(1,IB) = KAPU
            MCELL(2,IB) = LAPU
         ENDIF

         NCELLZ(I) = NCELLZ(I) + NCELL(2,IB) ! First ghost cells included
         NCELLT(I) = NCELLT(I) + MCELL(2,IB) ! Ghost cells excluded

         DO J=NC,JAKOPTR(I)+1,-1
            JAKO(J) = JAKO(J-1)
         ENDDO
         DO J=I+1,NCPU
            JAKOPTR(J) = JAKOPTR(J) + 1
         ENDDO

         JAKO(JAKOPTR(I)) = NCELL(1,IB)
         JAKOPTR(I) = JAKOPTR(I) + 1

C ... Each INLET-/OUTLET-block must locate in the same process as the
C ... block where it is connected.
 
         DO IP = N1+1,N2
            IF((ICON(1,IP) == 1 .OR. ICON(1,IP) == 11) .AND. 
     &         ICON(8,IP) == NCELL(1,IB)) THEN
               NOFBPP(I) = NOFBPP(I) + 1

               DO J=NC,JAKOPTR(I)+1,-1
                  JAKO(J) = JAKO(J-1)
               ENDDO
               DO J=I+1,NCPU
                  JAKOPTR(J) = JAKOPTR(J) + 1
               ENDDO

               JAKO(JAKOPTR(I)) = ICON(2,IP)
               JAKOPTR(I) = JAKOPTR(I) + 1
               NCELLT(I)  = NCELLT(I)  + MCELL(2,ICON(2,IP)) 
            ENDIF
         ENDDO
      ELSE
         GOTO 843
      ENDIF

      GOTO 841

 843  CONTINUE

      IF(IMAX < NCPU) NCPU = IMAX-1
      K = NCPU


C ... Sort jako in increasing order in every process

      N = 0
      DO I = 1,K

 646     SORTED = .TRUE.
         DO J=N+2,N+NOFBPP(I)
            IF(JAKO(J) <= JAKO(J-1)) THEN
               IAPU      = JAKO(J-1)
               JAKO(J-1) = JAKO(J)
               JAKO(J)   = IAPU
               SORTED = .FALSE.
            ENDIF
         ENDDO
         IF(SORTED) GOTO 647
         GOTO 646
 647     N = N + NOFBPP(I)
      ENDDO
      
C ... Fill the NPROCE array

      N = 0
      DO I = 1,NPRO
         NPROCE(1,I) = NOFBPP(I)
         DO J = 2,NOFBPP(I)+1
            NPROCE(J,I) = JAKO(N+J-1)
         ENDDO
         N = N + NOFBPP(I)
      ENDDO

C ... Write the PROCES-file

      N = 0
      WRITE(68,109) K, 3
      DO I = 1,K
         WRITE(68,110) NOFBPP(I),I,NCELLT(I)
         WRITE(68,'(500I6)') (JAKO(J),J=N+1,N+NOFBPP(I))
         N = N + NOFBPP(I)
      ENDDO

 109  FORMAT(2I4,   '    ! TOTALNUMBER OF PROCESSES   DEFAULT IOTYPE')
 110  FORMAT(1I4,4X,'    ! BLOCKS IN PROCESSOR NUMBER',I4,':',I7,
     &              ' cells')

      CLOSE(68)

      DEALLOCATE(NBCSF)
      DEALLOCATE(NOFBPP)
      DEALLOCATE(JAKO)
      DEALLOCATE(JAKOPTR)
      DEALLOCATE(IPFIRST)
      DEALLOCATE(IPLAST)
      DEALLOCATE(NCELL)
      DEALLOCATE(MCELL)
      DEALLOCATE(NCELLZ)
      DEALLOCATE(NCELLT)
      WRITE(45,*) 'Deallocating memory for NBCSF etc. in MKPROC'
      WRITE(45,*)

      RETURN
      END SUBROUTINE MKPROC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MPICHE(IPRO)

      USE MPI

      LOGICAL :: FLAGI

C ... make some MPI enviromental checks and interrupt run if 
C ... something is wrong. Could be more. PPR 28.5.96

      CALL MPI_ATTR_GET(MPI_COMM_WORLD,MPI_TAG_UB,ITGI,FLAGI,IERROR)
      IF(.NOT.FLAGI) THEN
         WRITE(*,*) 'MPI_TAG_UB is not something wrong'
         GOTO 1000
      ENDIF
C ... CHECK UPPER LIMIT OF MPI_TAG_UB
      IF(ITGI < 32767) THEN
         WRITE(*,*) 'MPI_TAG_UB might be too small',ITGI
         GOTO 1000
      ENDIF

      CALL MPI_ATTR_GET(MPI_COMM_WORLD,MPI_IO,IOMP,FLAGI,IERROR)
*     CALL MPI_COMM_GET_ATTR(MPI_COMM_WORLD,MPI_IO,IOMP,FLAGI,IERROR)
      IF(.NOT.FLAGI) THEN
         WRITE(*,*) 'MPI_IO is not something wrong'
         GOTO 1000
      ENDIF
C ... CHECK IF ALL PROCESES ARE ABLE TO COMMUNICATE
*      IF(MPI_ANY_SOURCE /= IOMP) THEN
*         WRITE(*,*) 'PROCES',IPRO,' IS NOT ABLE TO COMMUNICATE'
*         GOTO 1000
*      ENDIF

      WRITE(45,*)
      WRITE(45,*) 'MPI_TAG_UB (tag upper limit)   is ',ITGI
      WRITE(45,*) 'MPI_WTICK (accuracy in timing) is ',MPI_WTICK(), ' s'
      WRITE(45,*)
      GOTO 2000

 1000 WRITE(*,*) ' Exiting in proces',IPRO
      CALL MPI_ABORT(MPI_COMM_WORLD,IERRORCODE,IERR)
      STOP

 2000 CONTINUE

      RETURN
      END SUBROUTINE MPICHE

C     ****************************************************************
C     *        COMMUNICATION SUBROUTINES                             *
C     ****************************************************************

      SUBROUTINE SEFILE(ZZZ,JPRO,IP,NTOT,IOTY)

      USE MPI

      REAL :: ZZZ(*)

C ... SEND VECTOR ZZZ TO PACHT IP IN PROSESS JPRO BY USING MPI
      IF(IOTY == 3) THEN
         CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,JPRO-1,IP,MPI_COMM_WORLD,
     +        IERR)
       ELSE
         WRITE(*,*) 'None existing iotype=',IOTY
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      RETURN
      END SUBROUTINE SEFILE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REFILE(ZZZ,JPRO,IFRO,IP,NTOT,IOTY)

      USE MPI

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL :: ZZZ(*)
      
C ... RECEIVE VECTOR ZZZ FROM PACHT IP IN PROSESS JPRO  BY USING MPI
      IF(IOTY == 3) THEN
         CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,IFRO-1,IP,MPI_COMM_WORLD,
     +        STATUS,IERR)
      ELSE
         WRITE(*,*) 'None existing iotype=',IOTY
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      RETURN
      END SUBROUTINE REFILE

C     ***************************************************************
C     *       START UP SUBROUTINES                                  *
C     ***************************************************************

      SUBROUTINE MPITIN(NBLOCK,INPARA,INPAR,NPROCE,MBPRO,NPRO) 

      USE MPI

      INTEGER :: NPROCE(MBPRO+1,*),INPARA(INPAR,*)!,IAA(1000)
      LOGICAL :: CHECK

      IF(INPAR > 1000) STOP 'TASSA 11'

C ... SEND INPUT FROM MASTER TO THE SLAVE PROCESSES BY USING MPI
C ... These inputs are not changed from processor to processor
C ... just picking up right blocks. 
C ... Subroutine MPITOU is a communication partner with this subroutine
C ... INPARA contains all block data.

      WRITE(45,*)
      WRITE(45,*) 'MPITIN: sending INPUT'

      DO IPRO = 1,NPRO

         DO NBL = 1,NPROCE(1,IPRO)

            NBG = NPROCE(NBL+1,IPRO)

            IF(IPRO /= 1) THEN 
               CALL MPI_SEND(INPARA(1,NBG),INPAR,MPI_INTEGER,IPRO-1,
     +              900+NBL,MPI_COMM_WORLD,IERR)
            ENDIF

         ENDDO

      ENDDO

      WRITE(45,*) 'MPITIN: INPUT sent'
      WRITE(45,*)

      RETURN
      END SUBROUTINE MPITIN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MPITOU(NBLOCK,IPRO,NPROCE,MBPRO,NPRO,IMAX,JMAX,KMAX,
     +     INTERI,INTERJ,INTERK,IDER,LAMIN,INITC,IT,IL,IK,IDI1,IDI2,
     +     IDI3,MOV,MGRID,MIB,MIT,MJB,MJT,MKB,MKT,IROTVE,MGM,
     *     INPAR,INPARA) 

      USE MPI

      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),
     +           INTERI(*),INTERJ(*),INTERK(*),IDER(*),LAMIN(*),
     +           INITC(*),IT(MGM,*),IL(MGM,*),IK(MGM,*),
     +           IDI1(*),IDI2(*),IDI3(*),
     +           MGRID(*),MIB(*),MIT(*),MJB(*),MJT(*),MKB(*),MKT(*),
     +           IROTVE(*),IAA(INPAR),MOV(*),STATUS(MPI_STATUS_SIZE),
     +           INPARA(INPAR,*),NPROCE(MBPRO+1,*)

      LOGICAL :: CHECK,CHECK2

C ... RECEIVE INPUT FROM MASTER TO THE SLAVE AND ALSO TO THE MASTER
C ... PROCESS BY USING MPI

C ... Subroutines MPITIN and MPITOU rewritten so that the master 
C ... process does not send messages to itself. ESa Mar 3, 2007. 

C ... Subroutine MPITIN is communication partner to this subroutine.

      DO 30 NBL = 1,NBLOCK

         IF(IPRO /= 1) THEN
         CALL MPI_RECV(IAA,INPAR,MPI_INTEGER,0,900+NBL,MPI_COMM_WORLD,
     +        STATUS,IERR)
         IMAX(1,NBL) = IAA(1) 
         JMAX(1,NBL) = IAA(2) 
         KMAX(1,NBL) = IAA(3) 
         INTERI(NBL) = IAA(4) 
         INTERJ(NBL) = IAA(5) 
         INTERK(NBL) = IAA(6)
         IDER(NBL)   = IAA(24) ! Could be changed some day..
         LAMIN(NBL)  = IAA(7) 
         INITC(NBL)  = IAA(8) 
         IT(1,NBL)   = IAA(9) 
         IL(1,NBL)   = IAA(10)
         IK(1,NBL)   = IAA(11)
         IDI1(NBL)   = IAA(12)
         IDI2(NBL)   = IAA(13)
         IDI3(NBL)   = IAA(14)
         MOV(NBL)    = IAA(15)
         MGRID(NBL)  = IAA(16)
         MIB(NBL)    = IAA(17)
         MIT(NBL)    = IAA(18)
         MJB(NBL)    = IAA(19)
         MJT(NBL)    = IAA(20)
         MKB(NBL)    = IAA(21)
         MKT(NBL)    = IAA(22)
         IROTVE(NBL) = IAA(23)
         ELSE
         NBG    = NPROCE(NBL+1,IPRO)
         IMAX(1,NBL) = INPARA(1,NBG) 
         JMAX(1,NBL) = INPARA(2,NBG) 
         KMAX(1,NBL) = INPARA(3,NBG) 
         INTERI(NBL) = INPARA(4,NBG) 
         INTERJ(NBL) = INPARA(5,NBG) 
         INTERK(NBL) = INPARA(6,NBG)
         IDER(NBL)   = INPARA(24,NBG) ! Could be changed some day..
         LAMIN(NBL)  = INPARA(7,NBG) 
         INITC(NBL)  = INPARA(8,NBG) 
         IT(1,NBL)   = INPARA(9,NBG) 
         IL(1,NBL)   = INPARA(10,NBG)
         IK(1,NBL)   = INPARA(11,NBG)
         IDI1(NBL)   = INPARA(12,NBG)
         IDI2(NBL)   = INPARA(13,NBG)
         IDI3(NBL)   = INPARA(14,NBG)
         MOV(NBL)    = INPARA(15,NBG)
         MGRID(NBL)  = INPARA(16,NBG)
         MIB(NBL)    = INPARA(17,NBG)
         MIT(NBL)    = INPARA(18,NBG)
         MJB(NBL)    = INPARA(19,NBG)
         MJT(NBL)    = INPARA(20,NBG)
         MKB(NBL)    = INPARA(21,NBG)
         MKT(NBL)    = INPARA(22,NBG)
         IROTVE(NBL) = INPARA(23,NBG)
         ENDIF
 30   CONTINUE


      RETURN
      END SUBROUTINE MPITOU
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INPUTO(NBLOCK,IMAX,JMAX,KMAX,INTERI,INTERJ,INTERK,
     +     IDER,LAMIN,INITC,IT,IL,IK,IDI1,IDI2,IDI3,MOV,MGRID,
     +     MIB,MIT,MJB,MJT,MKB,MKT,IROTVE,MGM,INPARA,INPAR) 

      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),
     +           INTERI(*),INTERJ(*),INTERK(*),IDER(*),LAMIN(*),
     +           INITC(*),IT(MGM,*),IL(MGM,*),IK(MGM,*),
     +           IDI1(*),IDI2(*),IDI3(*),
     +     MGRID(*),MIB(*),MIT(*),MJB(*),MJT(*),MKB(*),MKT(*),
     +     IROTVE(*),IAA(23),MOV(*),INPARA(INPAR,*)
      LOGICAL :: CHECK,CHECK2

C ... PUT INPUT TO THE NORMAL MODE

      DO 30 NBL = 1,NBLOCK

         IMAX(1,NBL) = INPARA(1 ,NBL) 
         JMAX(1,NBL) = INPARA(2 ,NBL) 
         KMAX(1,NBL) = INPARA(3 ,NBL) 
         INTERI(NBL) = INPARA(4 ,NBL) 
         INTERJ(NBL) = INPARA(5 ,NBL) 
         INTERK(NBL) = INPARA(6 ,NBL)
c         IDER(NBL)   = 2 ! Constant in the past ...
         IDER(NBL)   = INPARA(24,NBL) ! Numbering cuold be changed some day..
         LAMIN(NBL)  = INPARA(7 ,NBL) 
         INITC(NBL)  = INPARA(8 ,NBL) 
         IT(1,NBL)   = INPARA(9 ,NBL) 
         IL(1,NBL)   = INPARA(10,NBL)
         IK(1,NBL)   = INPARA(11,NBL)
         IDI1(NBL)   = INPARA(12,NBL)
         IDI2(NBL)   = INPARA(13,NBL)
         IDI3(NBL)   = INPARA(14,NBL)
         MOV(NBL)    = INPARA(15,NBL)
         MGRID(NBL)  = INPARA(16,NBL)
         MIB(NBL)    = INPARA(17,NBL)
         MIT(NBL)    = INPARA(18,NBL)
         MJB(NBL)    = INPARA(19,NBL)
         MJT(NBL)    = INPARA(20,NBL)
         MKB(NBL)    = INPARA(21,NBL)
         MKT(NBL)    = INPARA(22,NBL)
         IROTVE(NBL) = INPARA(23,NBL)
 30   CONTINUE

      RETURN
      END SUBROUTINE INPUTO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MPITBC(ICON,RCON,NBCS,IBC,RBC,NPROCE,MBPRO,
     +                  NPRO,NBLOCK,BCFILE,IOTY,BCAPU,BNFILE,BNAPU)

C ... This subroutine splits boundary conditions (ICON array) for
C ... multi-processing.

      USE MPI
      USE NS3CO, ONLY : IC9

      INTEGER           :: ICON(IC9,*), IBC(IC9,*), NPROCE(MBPRO+1,*)
      CHARACTER(LEN=80) :: BCFILE(*), BCAPU(*), BNFILE(*), BNAPU(*)
      LOGICAL           :: CHECK, MPICH, MPICH2
      REAL              :: RCON(IC9,*), RBC(IC9,*)

C ... SEND ICON FILE FROM THE MASTER TO THE SLAVES BY USING MPI
C ... this is quite complex structure. Subroutine does division between
C ... different blocks. 

      WRITE(45,*)
      WRITE(45,*) 'MPITBC: sending BC to slaves'

      DO 1000 IPRO = NPRO,1,-1
      MPICH = .FALSE.
      IOINT = 1     ! default internal communication type
      IOOUT = IOTY  ! default outternal communication type

      DO 5 NP = 1,NBCS
         IF(ICON(18,NP) == 3) MPICH =.TRUE.
 5    CONTINUE

      IF(MPICH) THEN
         CALL MPI_INITIALIZED(MPICH2,IERR)
         IF(.NOT. MPICH2) CALL MPI_INIT(IERR)
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      ENDIF
c      IF(NPRO == 1) RETURN

C ... CHECK THAT ALL BLOCKS ARE FOUND IN NPROCE         
      DO 10 N = 1,NBLOCK
         CHECK = .FALSE.
         DO 20 J = 1,NPRO
         DO 20 I = 2,NPROCE(1,J)+1
            IF (CHECK .AND. (N == NPROCE(I,J))) STOP 'Multiple blocks'
            IF (N == NPROCE(I,J)) CHECK = .TRUE.
 20      CONTINUE
         IF(.NOT. CHECK) THEN
            WRITE(*,*) 'Cannot find block ',N,' in PROCES'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
 10   CONTINUE


C ... TAKE LOCAL BC DATA FOR PROCESS IPRO

      NPL = 0
      DO 100 NPG = 1,NBCS
         DO 30 J = 1,NPROCE(1,IPRO)
            IF(ICON(2,NPG) == NPROCE(J+1,IPRO)) THEN
               NPL = NPL + 1

               DO 29 IT = 1,IC9
                  IBC(IT,NPL) = ICON(IT,NPG)
                  RBC(IT,NPL) = RCON(IT,NPG)
 29            CONTINUE
C ...
            ENDIF
 30      CONTINUE
 100  CONTINUE
      NPCSL  = NPL
      WRITE(45,*) '  Number of patches in this prosess(IBC):',NPCSL,IPRO

C ... CHECK CONNECTION, CYCLIC BOUNDARIES, MIRRORS
      DO 200 NPL = 1,NPCSL
         ITYP = IBC(1,NPL)
         IF (ITYP == 1 .OR. ITYP == 6 .OR. ITYP == 11 .OR. 
     +       ITYP == 4 .OR. ITYP  == 13 .OR. ITYP == 14 .OR.
     +       ITYP  == 15) THEN
C ... FIND NEW LOCAL AND GLOBAL PATCH NUMPERS, 
C ... BLOCK NUMBER AND PROSESS NUMBER WHERE TO CONNECT.
            DO 150 NPR = 1,NPRO
            DO 150 J = 1,NPROCE(1,NPR)
               IF(IBC(8,NPL) == NPROCE(J+1,NPR)) THEN
                  MPROCE = NPR
                  LOCABL = J
                  GOTO 155
               ENDIF
 150        CONTINUE
C ... MPROCE IS THE NUMBER OF PROCESS TO CONNECT
C ... LOCALBL IS THE NUMBER OF LOCAL BLOCK TO CONNECT
 155        LOCAPA = 0
            NP     = 0
 160        NP = NP + 1
            DO 170 J = 1,NPROCE(1,MPROCE)
               IF(ICON(2,NP) == NPROCE(J+1,MPROCE)) THEN
                  LOCAPA = LOCAPA + 1
                  GOTO 175
               ENDIF
 170        CONTINUE
 175        IF(IBC(16,NPL) == NP) GOTO 180
            GOTO 160
 180        CONTINUE
C ... LOCAPA IS THE LOCAL PATCH NUMBER IN CONNECTIVE PROCESS

            IBC( 8,NPL) = LOCABL
            IBC(16,NPL) = LOCAPA
            IBC(16,NPL) = LOCAPA
            IBC(17,NPL) = MPROCE
 
            IF(IBC(17,NPL) /= IPRO .AND. IBC(18,NPL)  == 1) THEN
               IBC(18,NPL) = IOOUT
            ENDIF

            IF(IBC(18,NPL) == 2) MPICH = .TRUE.

            IF(IBC(1,NPL) == 4 .OR. IBC(1,NPL) == 13) THEN
               IBC(18,NPL) = 0
               IBC(17,NPL) = 0
            ENDIF
         ENDIF !(ITYP == 1 .OR. ITYP == 6 .....)
 200  CONTINUE

C ... DROP BLK IN LOCAL MODE
      DO 300 NPL = 1,NPCSL
         DO 270 J = 1,NPROCE(1,IPRO)
            IF(IBC(2,NPL) == NPROCE(J+1,IPRO)) THEN
               IBC(2,NPL) = J
               GOTO 300
            ENDIF
 270     CONTINUE
 300  CONTINUE
      DO 35 IP = 1,NPCSL
*         IF (IBC(1,IP) == 3 .OR. IBC(1,IP) == 5) THEN
           BCAPU(IP) = BCFILE(IBC(25,IP))
           BNAPU(IP) = BNFILE(IBC(25,IP)) 
*         ENDIF
 35   CONTINUE

      WRITE(45,*) 'MPITBC: sending BC to PROCES',IPRO
      IF(IPRO /= 1) THEN
      CALL MPI_SEND(NPCSL,1,MPI_INTEGER,IPRO-1,100+IPRO,MPI_COMM_WORLD,
     +     IERR)
      CALL MPI_SEND(IBC(1,1),IC9*NPCSL,MPI_INTEGER,IPRO-1,300+IPRO,
     +     MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(RBC(1,1),IC9*NPCSL,MPI_REAL8,IPRO-1,1300+IPRO,
     +     MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(BCAPU(1),80*NPCSL,MPI_CHARACTER,IPRO-1,500+IPRO,
     +     MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(BNAPU(1),80*NPCSL,MPI_CHARACTER,IPRO-1,1500+IPRO,
     +     MPI_COMM_WORLD,IERR)
      ELSE ! master proces

         NBCS = NPCSL

         DO NPL = 1,NBCS

            DO I = 1,IC9
               ICON(I,NPL) = IBC(I,NPL) 
               RCON(I,NPL) = RBC(I,NPL)
            ENDDO

            BCFILE(NPL) = BCAPU(NPL)
            BNFILE(NPL) = BNAPU(NPL)

         ENDDO
      ENDIF

      WRITE(45,*) 'MPITBC: sent    BC to PROCES',IPRO
      

      IBEDEL = 0
      WRITE(45,1102) IPRO
      WRITE(45,1105)
 1102 FORMAT(//4X,'BOUNDARY CONDITION PATCHES ON THE FINEST ',
     &            'GRID LEVEL FOR PROCESS:',I4)
 1105 FORMAT(//4X,' 1  = Connectivity     2  = External     ',
     &            ' 3  = Inlet',
     &        /4X,' 4  = Mirror           5  = Outlet       ',
     &            ' 6  = Cyclic',
     &        /4X,' 7  = Singularity      8  = Solid        ',
     &            ' 9  = Rotating solid',
     &        /4X,' 10 = Moving solid    11  = Sliding      ',
     &            ' 12 = Chimera',
     &        /4X,' 13 = Free surface    14  = Mixing       ',
     &            ' 15 = Coupling with solid 16  = SLIP     ',
     &        /4X ' 17 = Circulating inlet')
c 1103 FORMAT(//4X,'BLOCK NO. ',I4,'                     CONNECTIVITY'//,
c     &         4X,'BC BLK FACE XLO XUP  YLO  YUP BLK FACE ',
c     &'XLO XUP YLO  YUP EMP GPN GPN PRO TY MGS PER'/)
 1103 FORMAT(//4X,'BLOCK NO.',I6,'                     CONNECTIVITY'//,
     &         4X,'BC   BLK FACE XLO  XUP  YLO  YUP    BLK  FACE ',
     &'XLO  XUP  YLO  YUP EMP GPN   GPN PRO  TY MGS PER'/)

      DO 17 I=1,NPCSL
        IF(IBC(2,I) > IBEDEL) THEN
          WRITE(45,1103) IBC(2,I)
          IBEDEL = IBC(2,I)
        ENDIF
        WRITE(45,111) (IBC(J,I),J=1,20)
c 111    FORMAT(4X,3I3,4I5,2I3,4I5,I3,I5,I4,2I3,4I4)
 111    FORMAT(3X,I3,I6,I4,4I5,1X,I6,I5,4I5,I4,I4,I6,5I4)
 17   CONTINUE
      WRITE(45,111)
      WRITE(45,111)
      WRITE(45,111)
 1000 CONTINUE !IPRO = NPRO,1,-1
      WRITE(45,*) 'MPITBC: done'
      RETURN
      END SUBROUTINE MPITBC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRIMP(XCO,YCO,ZCO,IX,IY,IZ,LEVEL,GRILEN,
     +           IG,IPRO,NPROCE,MBPRO,NPRO,NBLOCK,MGM,IN,JN,KN)

      IMPLICIT NONE

      INTEGER :: NGRID,IN,JN,KN,IERRCODE,LEVEL,MNDSTP,NPRO,
     &        IPRO,MBPRO,NBLOCK,MGM,N,IG1,INDSTP,JNDSTP,KNDSTP,
     &        IMAXP1,JMAXP1,KMAXP1,IMAX,JMAX,KMAX,ISTRID,JSTRID,
     &        KSTRID,ILE,NTOT,I,K,KA,J,JJ,KK,IL,IL2

      INTEGER :: IX(*), IY(*), IZ(*), IG(MGM,*)
      INTEGER :: NPROCE(MBPRO+1,*)

      REAL :: XCO(*), YCO(*), ZCO(*)
      REAL, ALLOCATABLE :: X(:),Y(:),Z(:)

      REAL :: GRILEN, BKOKO

      MNDSTP = 2**(LEVEL-1)

C *** READ GRID FILE ***************************************************

      READ (41) NBLOCK
      READ (41) (IX(N),IY(N),IZ(N),N=1,NBLOCK)
      IF (NPRO /= 1) THEN
         WRITE(*,*) 'Something wrong in GRIMP. NPRO=',NPRO
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      DO 900 N = 1,NBLOCK
         IG1    = IG(1,N) - 1
         INDSTP = MNDSTP
         JNDSTP = MNDSTP
         KNDSTP = MNDSTP
         IF(IX(N) == 3) INDSTP = 1
         IF(IY(N) == 3) JNDSTP = 1
         IF(IZ(N) <= 3) KNDSTP = 1
         IMAXP1 = (IX(N)-1)/INDSTP+1
         JMAXP1 = (IY(N)-1)/JNDSTP+1
         KMAXP1 = (IZ(N)-1)/KNDSTP+1

         IMAX    = IMAXP1 - 1
         JMAX    = JMAXP1 - 1
         KMAX    = KMAXP1 - 1
         ISTRID  = IMAX + 2*IN
         JSTRID  = JMAX + 2*JN
         KSTRID  = KMAX + 2*KN
         ILE     = ISTRID*JSTRID
         NTOT    = IMAXP1*JMAXP1*KMAXP1
         NGRID   = IX(N)*IY(N)*IZ(N)

         ALLOCATE (X(NGRID),Y(NGRID),Z(NGRID),STAT=IERRCODE)
         IF(IERRCODE /= 0) THEN
           WRITE(*,*) 'Could not allocate space to read in the grid',
     &     IERRCODE
           WRITE(45,*)
           WRITE(45,*) 'Could not allocate space to read in the grid',
     &     IERRCODE
           CALL EXIT ! I made this
           ELSE
           WRITE(45,*)
           BKOKO = REAL(3*NGRID)*8./1048576.
           WRITE(45,9060) N,3*NGRID,BKOKO
9060       FORMAT('In GRIMP: Allocated temporary arrays X, Y and Z for ' 
     &     'block',I3,' Total size (DP) = ',I8,' = ',F6.2,'MB')
         ENDIF
         
C ... READ GRID AND TRANFER IT INTO THE DESIRED LEVEL

         READ(41)(X(I),I=1,NGRID),(Y(I),I=1,NGRID),(Z(I),I=1,NGRID)

         DO 130 K = 1,KMAXP1
         KA    = (KN+K-1)*ILE
         DO 120 J = 1,JMAXP1
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO 110 I = 1,IMAXP1
               IL =1+((I-1)*INDSTP)+((J-1)*JNDSTP)*IX(N)+
     &               ((K-1)*KNDSTP)*IY(N)*IX(N)
               IL2=1+((I-1)       )+((J-1)       )*IMAXP1+
     &               ((K-1)       )*IMAXP1*JMAXP1
               X(IL2) = X(IL)
               Y(IL2) = Y(IL)
               Z(IL2) = Z(IL)
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE

      DO 170 K = 1,KMAXP1
         KA    = (KN+K-1)*ILE
         DO 160 J = 1,JMAXP1
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO 150 I = 1,IMAXP1
               KK      = JJ + I
               IL2=1+((I-1)       )+((J-1)       )*IMAXP1+
     &               ((K-1)       )*IMAXP1*JMAXP1
               XCO(KK+IG1) = GRILEN*X(IL2)
               YCO(KK+IG1) = GRILEN*Y(IL2)
               ZCO(KK+IG1) = GRILEN*Z(IL2)
 150        CONTINUE
 160     CONTINUE
 170  CONTINUE

      DEALLOCATE (X,Y,Z,STAT=IERRCODE)
      WRITE(45,9070)
9070  FORMAT('X, Y and Z deallocated in GRIMP')
 900  CONTINUE

      RETURN
      END SUBROUTINE GRIMP      
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRIDIN(XCO,YCO,ZCO,IX,IY,IZ,LEVEL,GRILEN,X,Y,Z,
     &     NB,IPRO)

      USE MPI
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      REAL :: X(*),Y(*),Z(*),XCO(*),YCO(*),ZCO(*)

      REAL :: GRILEN

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IX,IY,IZ,LEVEL,NB,IPRO,INDSTP,JNDSTP,KNDSTP,
     &        MNDSTP,IMAX,JMAX,KMAX,IMAXP1,JMAXP1,KMAXP1,ILE,NTOT,
     &        ISTRID,JSTRID,KSTRID,I,J,K,III,IL2,JJ,KK,KA,IERR

C ... READ ONE BLOCK OF THE GRID IN

      MNDSTP = 2**(LEVEL-1)
      INDSTP = MNDSTP
      JNDSTP = MNDSTP
      KNDSTP = MNDSTP
      IF(IX == 3) INDSTP = 1
      IF(IY == 3) JNDSTP = 1
      IF(IZ <= 3) KNDSTP = 1
      IMAXP1 = (IX-1)/INDSTP+1
      JMAXP1 = (IY-1)/JNDSTP+1
      KMAXP1 = (IZ-1)/KNDSTP+1

      IMAX    = IMAXP1 - 1
      JMAX    = JMAXP1 - 1
      KMAX    = KMAXP1 - 1
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ILE     = ISTRID*JSTRID
      NTOT    = IMAXP1*JMAXP1*KMAXP1

C ... RECEIVE GRID FROM MASTER PROCESS IF MPI IS USED

      IF(IPRO /= 1) THEN
         CALL MPI_RECV(X(1),NTOT,MPI_REAL8,0,110+NB,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(Y(1),NTOT,MPI_REAL8,0,120+NB,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(Z(1),NTOT,MPI_REAL8,0,130+NB,
     +        MPI_COMM_WORLD,STATUS,IERR)
      ELSE
         DO III = 1,NTOT
            X(III)   = XCO(III)
            Y(III)   = YCO(III)
            Z(III)   = ZCO(III)
            XCO(III) = 0.0
            YCO(III) = 0.0
            ZCO(III) = 0.0
         ENDDO
      ENDIF

      DO 150 K = 1,KMAXP1
         KA    = (KN+K-1)*ILE
         DO 160 J = 1,JMAXP1
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO 170 I = 1,IMAXP1
               KK      = JJ + I
               IL2=1+((I-1)       )+((J-1)       )*IMAXP1+
     &               ((K-1)       )*IMAXP1*JMAXP1
               XCO(KK) = GRILEN*X(IL2)
               YCO(KK) = GRILEN*Y(IL2)
               ZCO(KK) = GRILEN*Z(IL2)
 170        CONTINUE
 160     CONTINUE
 150  CONTINUE

      RETURN
      END SUBROUTINE GRIDIN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRMPI(XCO,YCO,ZCO,IX,IY,IZ,IXG,IYG,IZG,LEVEL,
     +     GRILEN,IG,IPRO,NPROCE,MBPRO,NPRO,NBLOCK,MGM)

      USE MPI

      IMPLICIT NONE

      INTEGER :: IERRCODE,NGRID,NB,NBLOCG,IPROG,N,NBG,NBL,LEVEL,
     +        IPRO,MBPRO,NPRO,NBLOCK,MGM,IPROL,IERR,IG1

      INTEGER :: IX(*),IY(*),IZ(*),IG(MGM,*)
      INTEGER :: NPROCE(MBPRO+1,*),STATUS(MPI_STATUS_SIZE),
     +        IXG(*),IYG(*),IZG(*)

      REAL :: XCO(*),YCO(*),ZCO(*)
      REAL, ALLOCATABLE :: X(:),Y(:),Z(:)

      REAL :: GRILEN, BKOKO

C ... READ, SEND AND RECEIVE THE GRID WHEN MPI IS USED

      IF(IPRO == 1) THEN

C *** READ GRID FILE IN A MASTER PROCESS *******************

      READ (41) NBLOCG
      READ (41) (IXG(N),IYG(N),IZG(N),N=1,NBLOCG)

      DO 110 NB = 1,NBLOCG

C ... FIND THE NUMBER OF THE PROCESS

         DO 100 IPROG = 1,NPRO
            DO 100 N = 1,NPROCE(1,IPROG)
               NBG = NPROCE(N+1,IPROG)
               IF(NBG == NB) THEN
                  IPROL = IPROG ! process number
                  NBL   = N     ! process local block number
               ENDIF
 100     CONTINUE

C ... SEND THE GRID AND DIMENSIONS TO THE SLAVES

         NGRID  = IXG(NB)*IYG(NB)*IZG(NB)
         ALLOCATE (X(NGRID),Y(NGRID),Z(NGRID),STAT=IERRCODE)
         IF(IERRCODE /= 0) THEN
           WRITE(*,*) 'Could not allocate space to read in the grid',
     &     ' in subroutine GRMPI',IERRCODE
           WRITE(45,*)
           WRITE(45,*) 'Could not allocate space to read in the grid',
     &     IERRCODE
           CALL EXIT ! I made this
         ELSE
           WRITE(45,*)
           BKOKO = REAL(3*NGRID)*8./1048576.
           WRITE(45,9060) NB,3*NGRID,BKOKO,IPROL
9060       FORMAT('In GRMPI: Allocated temporary arrays X, Y and Z for ' 
     &     'block',I3,/' Total size (DP) = ',I8,' = ',F6.2,'MB',
     &     ', Process no.',I3)
         ENDIF

         IF(IPROL /= 1) THEN
         CALL MPI_SEND(IXG(NB),1,MPI_INTEGER,IPROL-1,140+NBL,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(IYG(NB),1,MPI_INTEGER,IPROL-1,150+NBL,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(IZG(NB),1,MPI_INTEGER,IPROL-1,160+NBL,
     +        MPI_COMM_WORLD,IERR)
         ENDIF
         CALL GRMPIS(IXG(NB),IYG(NB),IZG(NB),LEVEL,X,Y,Z,IPROL,NBL,
     +        XCO,YCO,ZCO,IG,MGM) ! Read and send the grid
         DEALLOCATE (X,Y,Z)
         WRITE(45,*)'Read and send the grid. X, Y, and Z deallocated'     

 110  CONTINUE

      ENDIF !(IPRO == 1) ! The grid has been read into array 'XYZ'

      DO 200 N = 1,NPROCE(1,IPRO)

         IF(IPRO /= 1) THEN
         CALL MPI_RECV(IX(N),1,MPI_INTEGER,0,140+N,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(IY(N),1,MPI_INTEGER,0,150+N,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(IZ(N),1,MPI_INTEGER,0,160+N,
     +        MPI_COMM_WORLD,STATUS,IERR)
         ELSE
              NBG = NPROCE(N+1,IPRO)
              IX(N) = IXG(NBG)
              IY(N) = IYG(NBG)
              IZ(N) = IZG(NBG)
         ENDIF

         NGRID  = IX(N)*IY(N)*IZ(N)
         ALLOCATE (X(NGRID),Y(NGRID),Z(NGRID),STAT=IERRCODE)
         IF(IERRCODE > 0) THEN
           WRITE(*,*) 'Could not allocate space to transfer the grid',
     &     ' in subroutine GRMPI',IERRCODE
           WRITE(45,*)
           WRITE(45,*) 'Could not allocate space to read in the grid',
     &     IERRCODE
           CALL EXIT ! I made this
           ELSE
           WRITE(45,*)
           BKOKO = REAL(3*NGRID)*8./1048576.
           WRITE(45,9060) N,3*NGRID,BKOKO,IPRO
         ENDIF

         IG1    = IG(1,N)
         CALL GRIDIN(XCO(IG1),YCO(IG1),ZCO(IG1),IX(N),IY(N),IZ(N),
     +               LEVEL,GRILEN,X,Y,Z,N,IPRO)
         DEALLOCATE (X,Y,Z)
         WRITE(45,*)'The grid is received. X, Y, and Z are deallocated'     
       
 200  CONTINUE

      IF(NBLOCK /= NPROCE(1,IPRO)) THEN
         WRITE(*,*) 'Something is wrong with multiprocessing???',
     +        NBLOCK,NPROCE(1,IPRO),IPRO
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

C ... BETTER READING ROUTINE FOR A GRID MUST BE WRITTEN. PPR 13.10.95

      RETURN
      END SUBROUTINE GRMPI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRMPIS(IX,IY,IZ,LEVEL,X,Y,Z,IPROL,NBL,
     +           XCO,YCO,ZCO,IG,MGM)
C     *****************

      USE MPI
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IPROL,NBL,MGM,LEVEL,INDSTP,JNDSTP,KNDSTP,
     &        MNDSTP,IMAX,JMAX,KMAX,IMAXP1,JMAXP1,KMAXP1,ILE,NTOT,
     &        ISTRID,JSTRID,KSTRID,I,J,K,II,JJ,KK,III,IL,IL2,KA,
     &        IX,IY,IZ,IERR,IG1
      INTEGER :: IG(MGM,*)

      REAL :: X(*),Y(*),Z(*),XCO(*),YCO(*),ZCO(*)

C ... READ AND SEND THE GRID TO THE SLAVE PROCESSES (IPROL)

      MNDSTP = 2**(LEVEL-1)
      INDSTP = MNDSTP
      JNDSTP = MNDSTP
      KNDSTP = MNDSTP

      IF(IX == 3) INDSTP = 1
      IF(IY == 3) JNDSTP = 1
      IF(IZ <= 3) KNDSTP = 1

      READ(41) (X(I),I=1,IX*IY*IZ),
     &         (Y(I),I=1,IX*IY*IZ),
     &         (Z(I),I=1,IX*IY*IZ)

      IMAXP1 = (IX-1)/INDSTP+1
      JMAXP1 = (IY-1)/JNDSTP+1
      KMAXP1 = (IZ-1)/KNDSTP+1

      IMAX    = IMAXP1 - 1
      JMAX    = JMAXP1 - 1
      KMAX    = KMAXP1 - 1
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ILE     = ISTRID*JSTRID
      NTOT    = IMAXP1*JMAXP1*KMAXP1

      DO 150 K = 1,KMAXP1
         KA    = (KN+K-1)*ILE
         DO 160 J = 1,JMAXP1
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO 170 I = 1,IMAXP1
               IL =1+((I-1)*INDSTP)+((J-1)*JNDSTP)*IX+
     &               ((K-1)*KNDSTP)*IY*IX
               IL2=1+((I-1)       )+((J-1)       )*IMAXP1+
     &               ((K-1)       )*IMAXP1*JMAXP1
               KK      = JJ + I
               X(IL2) = X(IL)
               Y(IL2) = Y(IL)
               Z(IL2) = Z(IL)
 170        CONTINUE
 160     CONTINUE
 150  CONTINUE

      IF(IPROL /= 1) THEN
      CALL MPI_SEND(X(1),NTOT,MPI_REAL8,IPROL-1,110+NBL,
     +        MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(Y(1),NTOT,MPI_REAL8,IPROL-1,120+NBL,
     +        MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(Z(1),NTOT,MPI_REAL8,IPROL-1,130+NBL,
     +        MPI_COMM_WORLD,IERR)
      ELSE
         DO III = 1,NTOT
            IG1 = IG(1,NBL) + III - 1
            XCO(IG1) = X(III)
            YCO(IG1) = Y(III)
            ZCO(IG1) = Z(III)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GRMPIS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE READXY(IPRO,NBLOCG,NLOCAL,NPNUM,NPROCE,MBPRO,
     +     IMAXG,JMAXG,KMAXG,IMAX,JMAX,KMAX,MGM,NBLOCK,IG,
     +     IXG,IYG,IZG,IX,IY,IZ,XCO,YCO,ZCO)

      USE MPI
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER :: IPRO,NBLOCG,MBPRO,MGM,NCG,N,NROOT,NR1,
     +     IERR,IERRCODE,ERRORCODE,NBLOCK,N7,IM11,JM11,KM11,IMEM,IG1,
     +     IMAXP1,JMAXP1,KMAXP1,ISTRID,JSTRID,ILE,NNMAX,I,J,K,II,JJ 

      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),IG(MGM,*),
     +     IMAXG(*),JMAXG(*),KMAXG(*),NLOCAL(*),NPNUM(*),
     +     IXG(*),IYG(*),IZG(*),IX(*),IY(*),IZ(*),NPROCE(MBPRO+1,*)

      REAL, ALLOCATABLE :: FTR1(:)
      REAL :: XCO(*), YCO(*), ZCO(*)

      LOGICAL :: MASTER, SLAVE

C ... This subroutine reads grid from file XYD.BIN at restart

C *** READ GRID FILE IN A MASTER PROCESS *******************

      MASTER = IPRO == 1
      SLAVE  = .NOT.MASTER 
      
      IF(MASTER) THEN
         READ (41) NBLOCG
         READ (41) (IXG(N),IYG(N),IZG(N),N=1,NBLOCG)
      ENDIF

      DO 110 NCG = 1,NBLOCG

         N     = NLOCAL(NCG)    ! local block number
         NROOT = NPNUM(NCG)     ! process number
         NR1   = NROOT - 1      ! process number used by MPI

C ... Send the grid dimensions to the slaves

      IF(MASTER .AND. N == -1) THEN

      CALL MPI_SEND(IXG(NCG),1,MPI_INTEGER,NR1,NCG,MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(IYG(NCG),1,MPI_INTEGER,NR1,NCG,MPI_COMM_WORLD,IERR)
      CALL MPI_SEND(IZG(NCG),1,MPI_INTEGER,NR1,NCG,MPI_COMM_WORLD,IERR)

      ELSEIF(MASTER .AND. N /= -1) THEN

         IX(N) = IXG(NCG)
         IY(N) = IYG(NCG)
         IZ(N) = IZG(NCG)

      ELSEIF(SLAVE .AND. N /= -1) THEN  ! Receive the grid dimensions

         CALL MPI_RECV(IX(N),1,MPI_INTEGER,0,NCG,
     +        MPI_COMM_WORLD,STATUS,IERR)    
         CALL MPI_RECV(IY(N),1,MPI_INTEGER,0,NCG,
     +        MPI_COMM_WORLD,STATUS,IERR)    
         CALL MPI_RECV(IZ(N),1,MPI_INTEGER,0,NCG,
     +        MPI_COMM_WORLD,STATUS,IERR)

      ENDIF  ! MASTER .AND. N == -1

 110  CONTINUE
 
      IF(NBLOCK /= NPROCE(1,IPRO)) THEN
         WRITE(*,*) 'Something is wrong with multiprocessing???',
     +        NBLOCK,NPROCE(1,IPRO),IPRO
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

C ... BETTER READING ROUTINE FOR A GRID MUST BE WRITTEN. PPR 13.10.95
C ... BETTER READING ROUTINE FOR A GRID IS DONE. PPR 21.7.98

      IF(MASTER) THEN
         READ(28) N7
         READ(28) (IM11,JM11,KM11,N=1,N7)
         IF(N7 /= NBLOCG) THEN
            WRITE(*,*) ' Non matching XYD.BIN and INPUT'
            WRITE(*,*) ' Exiting ...'
            STOP
         ENDIF
      ENDIF

C ... Find the largest block

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

      IF (IERRCODE > 0) THEN
            WRITE(*,1022) IPRO,REAL(4*IMEM)/1048576.
            WRITE(45,1022) IPRO,REAL(4*IMEM)/1048576.
            CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
 1022 FORMAT(/'READXY :      Not enough memory in process ',I3,
     +     '. Desired ',F6.2,'MB. aborting...'/)
      ELSE
            WRITE(45,9070) IMEM
9070        FORMAT(4X,'FTR1 allocated in READXY. Size(DP)=',I8)
      ENDIF

      DO NCG = 1,NBLOCG

      N     = NLOCAL(NCG)       ! local block number
      NROOT = NPNUM(NCG)        ! process number
      NR1   = NROOT - 1         ! process number used by MPI

      IF(N /= -1) THEN  ! The block exists in this process

         IG1    = IG(1,N)
         IMAXP1 = IMAX(1,N) + 1
         JMAXP1 = JMAX(1,N) + 1
         KMAXP1 = KMAX(1,N) + 1
         ISTRID = IMAX(1,N) + 2*IN
         JSTRID = JMAX(1,N) + 2*JN
         ILE    = ISTRID*JSTRID
         NNMAX  = IMAXP1*JMAXP1*KMAXP1

         IF(SLAVE) THEN

            K = NCG
            CALL MPI_SEND(NNMAX,1,MPI_INTEGER,0,K,MPI_COMM_WORLD,IERR)
            CALL MPI_RECV(FTR1,3*NNMAX,MPI_REAL8,0,K,
     +           MPI_COMM_WORLD,STATUS,IERR)

         ELSE

            READ(28) (FTR1(I),I=1,3*NNMAX)

         ENDIF  ! SLAVE       

      ELSEIF(N == -1 .AND. MASTER) THEN  ! Master receives

       K = NCG
       CALL MPI_RECV(NNMAX,1,MPI_INTEGER,NR1,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
       READ(28) (FTR1(I),I=1,3*NNMAX)
        CALL MPI_SEND(FTR1,3*NNMAX,MPI_REAL8,NR1,K,
     +        MPI_COMM_WORLD,IERR)

      ENDIF                     !Master received data

      IF(N /= -1) THEN          ! The block exists in thi sprocess

         DO K=1,KMAXP1
         DO J=1,JMAXP1
         DO I=1,IMAXP1
            II = IG1+I+IN-1 + (J+JN-1)*ISTRID + (K+KN-1)*ILE
            JJ = I + (J-1)*IMAXP1 + (K-1)*IMAXP1*JMAXP1
            XCO(II) = FTR1(JJ        )
            YCO(II) = FTR1(JJ+  NNMAX)
            ZCO(II) = FTR1(JJ+2*NNMAX)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      ENDDO

      DEALLOCATE(FTR1)
      WRITE(45,*) 'FTR1 deallocated in READXY'

      RETURN
      END SUBROUTINE READXY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PRIME(MREQ,IPRO,NPRO,MBITS)

C ... Print memory requirement into the memory file and on the screen

      USE MPI

      USE NS3CO, ONLY : MAXPRO, PARALLEL

      IMPLICIT NONE

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      INTEGER :: IPRO, NPRO, IERR, NP, K

      INTEGER(KIND=8), DIMENSION(MAXPRO) :: IAA
      INTEGER(KIND=8) :: MREQ,MREQL,MREQB,MREQBL,MREQG,MREQBG,MBITS

      LOGICAL :: MASTER

      MREQB = MREQ*MBITS

      WRITE (45,1200) MREQ, MREQB, REAL(MREQB)/1048576.

      IF(PARALLEL) THEN

         IF(NPRO >= MAXPRO) THEN
            WRITE(*,*) 'MAXPRO must be changed from',MAXPRO,
     &        ' to at least',NPRO
            WRITE(*,*) 'in subroutine PRIME'
            WRITE(*,*) 'Aborting'
            STOP
         ENDIF

         IAA(1) = MREQ
           
         MASTER = IPRO == 1

         IF(.NOT.MASTER) THEN
            K = IPRO + 900
            CALL MPI_SEND(IAA(1),1,MPI_INTEGER8,0,K,MPI_COMM_WORLD,IERR)
         ELSEIF(MASTER) THEN
            DO NP = 2,NPRO
              K = NP + 900
              CALL MPI_RECV(MREQ,1,MPI_INTEGER8,NP-1,K,MPI_COMM_WORLD,
     +        STATUS,IERR)
              IAA(NP) = MREQ
            ENDDO
         ENDIF

         IF(.NOT.MASTER) GOTO 1300
            
         WRITE(*,1000)
         WRITE(45,1000)

         MREQG  = 0
         MREQBG = 0

         DO 200 NP = 1,NPRO
            
            MREQL  = IAA(NP)
            MREQBL = MREQL*MBITS
            MREQG  = MREQG + MREQL
            MREQBG = MREQBG+ MREQBL
c            WRITE (*,1100)  NP,MREQL,MREQBL,REAL(MREQBL)/1048576.
c            WRITE (45,1100) NP,MREQL,MREQBL,REAL(MREQBL)/1048576.
            WRITE (*,1100)  NP,REAL(MREQBL,8)/1048576.
            WRITE (45,1100) NP,REAL(MREQBL,8)/1048576.
 200     CONTINUE

c         WRITE(*,1150)  MREQG,MREQBG,REAL(MREQBG)/1048576.
c         WRITE(45,1150) MREQG,MREQBG,REAL(MREQBG)/1048576.
         WRITE(*,1150)  REAL(MREQBG,8)/1048576.
         WRITE(45,1150) REAL(MREQBG,8)/1048576.

c 1000 FORMAT(/,'     Memory requirements:'/
c     &         '     Process    words       bytes        MB'/
c     &         '     ===========================================')
c 1100 FORMAT(5X,I5,I12,I13,F10.3)
c 1150 FORMAT('     ==========================================='/
c     &  5X,'Total',I12,I13,F10.3)

 1000 FORMAT(/,'     Memory requirements:'/
     &         '     Process    MBytes'/
     &         '     ====================')
 1100 FORMAT(5X,I6,F12.3)
 1150 FORMAT('     ===================='/
     &  6X,'Total',F12.3)
    
      ELSE

*         IF(ABS(MREQB) > 1000000000) THEN
*            WRITE (*,1210) MREQ,MREQB,REAL(MREQB)/1073741824.
*         ELSE
*            WRITE (*,1200) MREQ,MREQB,REAL(MREQB)/1048576.
*         ENDIF

         IF(ABS(MREQB) > 1000000000) THEN
            WRITE (*,1410) REAL(MREQB,8)/1073741824.
         ELSE
            WRITE (*,1400) REAL(MREQB,8)/1048576.
         ENDIF

      ENDIF

 1300 CONTINUE

c 1200 FORMAT(/9X,'Memory requirement is',I11,' words of main memory'/
c     &        9X,'or equivalently',I17,' bytes of main memory (='
c     &        ,F7.2,'MB).')
c 1210 FORMAT(/9X,'Memory requirement is',I11,' words of main memory'/
c     &        9X,'or equivalently',I17,' bytes of main memory (='
c     &        ,F7.2,'GB).')
 1200 FORMAT(/9X,'Memory requirement is',I11,' words of main memory'/
     &        9X,'or equivalently',I17,' bytes of main memory (='
     &        ,F11.2,'MB).')
 1210 FORMAT(/9X,'Memory requirement is',I11,' words of main memory'/
     &        9X,'or equivalently',I17,' bytes of main memory (='
     &        ,F11.2,'GB).')

 1400 FORMAT(/2X,'Memory requirement is ',F7.2,' Mbytes')
 1410 FORMAT(/2X,'Memory requirement is ',F7.2,' Gbytes')

      RETURN
      END SUBROUTINE PRIME
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE READMP(LAITE,IFORM,RO,APUP2,IMAX,JMAX,KMAX,IN,JN,KN,
     &                  IPRO,NTG,N5)

C ... THIS SUBROUTINE READ COMPUT FILE AND ALSO RECEIVE DATA FROM
C ... MASTER PROCES

      USE MPI

      USE NS3CO, ONLY : PARALLEL

       
      IMPLICIT NONE
      INTEGER :: LAITE, IFORM, IMAX, JMAX, KMAX, IN, JN, KN, IPRO, 
     &           NTG, N5, ISTRID, JSTRID, ISTRIC, JSTRIC, NTOT,
     &           IERR, N, K, KK, KH, J, JJ, JH, I, IJ, IH
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: RO(*),APUP2(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      ISTRIC  = IMAX + IN
      JSTRIC  = JMAX + JN
      NTOT = (IMAX+2)*(JMAX+2)*(KMAX+2)
      IF(PARALLEL .AND. IPRO /= 1) THEN
C ... MPI IS USED AND PROCES IS SLAVE. RECEIVE DATA FROM MASTER
         CALL MPI_RECV(APUP2(1),NTOT,MPI_REAL8,0,NTG+N5,
     +        MPI_COMM_WORLD,STATUS,IERR)
      ELSE
C ... MASTER PROCES READ COMPUT ITSELF
         IF(IFORM == 0) THEN
            READ(LAITE)  (APUP2(N),N=1,NTOT)
         ELSE
            READ(LAITE,*)(APUP2(N),N=1,NTOT)
         ENDIF
      ENDIF
      N       = 0
      DO 3620 K = 1,KMAX
      KK      = (KN+K-1)*ISTRID*JSTRID
      KH      = (KN+K-2)*ISTRIC*JSTRIC
      DO 3620 J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KK
      JH      = (JN+J-2)*ISTRIC + IN + KH - 1
      DO 3620 I = 1,IMAX
      IJ      = I + JJ
      IH      = I + JH
      N       = N + 1
3620  RO(IJ)  = APUP2(IH)
      RETURN
      END SUBROUTINE READMP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COMRMP(NBLOCK,IG,IMAX,JMAX,KMAX,IN,JN,KN,MGM,
     +   F1R,F1RM,F1RN,F1RW,F1E,HAT1,F1RK,F1EPS,F1FI,MAXSB,NSCAL,ITURB,
     +   IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,VAR,IFORM,
     +   LAITE,MULPHL,BLKS,HAT2,TRANSL,TRM)

C ... READ COMPUT FILE BY MASTER PROCES AND SEND IT TO SLAVES
C ... MASTER PROCESS DATA IS NOT SEND BY MPI

      USE MPI
      USE TYPE_ARRAYS
      USE NS3CO, ONLY : PARALLEL

      IMPLICIT NONE

      TYPE(MPHASE_VARIABLES)  :: VAR(*)
      TYPE (BLOCKS)           :: BLKS(*)
      TYPE(INTERMITTENCY)     :: TRM(*)
      INTEGER :: NBLOCK,IN,JN,KN,MGM,MAXSB,NSCAL,ITURB,IPRO,
     +           MBPRO,NPRO,NBLOCG,IFORM,LAITE,N,NTG,NP,
     +           LB,IG1,KMAX2,NS,IS1,N6,NTO1,IGN,NGL,IPHASE,NTOT,
     +           IERR,N5,N4,IERRORCODE,IMUL
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: NPROCE(MBPRO+1,*)
      INTEGER :: IMAX(MGM,*), JMAX(MGM,*), KMAX(MGM,*), IG(MGM,*)
      LOGICAL :: MULPHL, TRANSL
      REAL    :: ZZZ(*),F1R(*),F1RM(*),F1RN(*),F1RW(*),F1E(*),F1RK(*),
     +           F1EPS(*),F1FI(*),HAT1(*),HAT2(*)
       
C ... SWEEP OVER GLOBAL BLOCK
      DO 1000 N = 1,NBLOCG              ! GLOBAL BLOCK LOOP
      NTG  = N*100 + 570                ! TAG NUMBER     
      DO 900 NP = 1,NPRO              ! PROCESSOR LOOP
      DO 800 LB=1,NPROCE(1,NP)          ! LOCAL BLOCK LOOP
cc      NGL  = NPROCE(LB+1,NP)
      IF(N == NPROCE(LB+1,NP)) THEN
      IF(NP == 1) THEN 
C ... MASTER PROCES OWNS THIS BLOCK AND IT IS READ ORDINARY WAY
         IG1     = IG(1,LB)
         NTG    = 0
         KMAX2 = MAX0(1,KMAX(2,LB))
         CALL READMP(LAITE,IFORM, F1R(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,1)
         CALL READMP(LAITE,IFORM,F1RM(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,2)
         CALL READMP(LAITE,IFORM,F1RN(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,3)
         CALL READMP(LAITE,IFORM,F1RW(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,4)
         CALL READMP(LAITE,IFORM, F1E(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,5)
         CALL READMP(LAITE,IFORM,HAT1(IG1),ZZZ,IMAX(2,LB),JMAX(2,LB),
     +        KMAX2,IN,JN,KN,IPRO,NTG,6)
         N6  = 6
         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            N6 = N6 + 1
            CALL READMP(LAITE,IFORM, F1RK(IG1),ZZZ,IMAX(2,LB),
     +           JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            N6 = N6 + 1
            CALL READMP(LAITE,IFORM,F1EPS(IG1),ZZZ,IMAX(2,LB),
     +           JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         ENDIF

C ... Intermittency variables
         IF(TRANSL) THEN
            NTO1 = (IMAX(2,LB)+2*IN)*(JMAX(2,LB)+2*JN)*(KMAX(2,LB)+2*KN)
            IG1  = IG(1,LB)
            IGN  = IG1 + NTO1 - 1
            NGL  = NPROCE(LB+1,NP)
            N6   = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            TRM(IG1:IGN)%G = HAT2(IG1:IGN)
            N6   = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            TRM(IG1:IGN)%RET = HAT2(IG1:IGN)
         ENDIF

C ... FOR SCALAR EQ. PPR 1.3
         DO 1900 NS = 1,NSCAL
            IS1 = IG1 + (NS-1)*MAXSB
            CALL READMP(LAITE,IFORM,F1FI(IS1),ZZZ,IMAX(2,LB),
     +           JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6+NS)
 1900    CONTINUE

C ... Multiphase variables
         IF(MULPHL) THEN
            N6   = N6 + NSCAL
            NTO1 = (IMAX(2,LB)+2*IN)*(JMAX(2,LB)+2*JN)*(KMAX(2,LB)+2*KN)
            IG1  = IG(1,LB)
            IGN  = IG1 + NTO1 - 1
            NGL  = NPROCE(LB+1,NP)
            DO IPHASE = 1,BLKS(NGL)%NPHASE
            N6   = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            VAR(IG1:IGN)%F1R(IPHASE) = HAT2(IG1:IGN)
            N6   = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            VAR(IG1:IGN)%F2R(IPHASE) = HAT2(IG1:IGN)

            IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
            N6 = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            VAR(IG1:IGN)%FRM(IPHASE) = HAT2(IG1:IGN) ! U-velocity
            N6 = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            VAR(IG1:IGN)%FRN(IPHASE) = HAT2(IG1:IGN) ! V-velocity
            N6 = N6 + 1
            CALL READMP(LAITE,IFORM, HAT2(IG1),ZZZ,
     +           IMAX(2,LB),JMAX(2,LB),KMAX2,IN,JN,KN,IPRO,NTG,N6)
            VAR(IG1:IGN)%FRW(IPHASE) = HAT2(IG1:IGN) ! W-velocity
            ENDIF ! MULTI

            ENDDO
         ENDIF
         GOTO 1000

      ELSE                      ! SLAVE BLOCK
C ... SLAVE NP OWNS THIS BLOCK AND IT IS READ HERE AND SEND TO THE SLAVE
C
C ... RECEIVE SIZE OF THE DESIRED VECTOR
         CALL MPI_RECV(NTOT,1,MPI_INTEGER,NP-1,NTG,MPI_COMM_WORLD,
     +        STATUS,IERR)
C ... MAIN VARIABLES
         DO 1950 N5 = 1,6
            IF(IFORM == 0) THEN
               READ(LAITE)   (ZZZ(N4),N4=1,NTOT)
            ELSE
               READ(LAITE,*) (ZZZ(N4),N4=1,NTOT)
            ENDIF
            CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +           MPI_COMM_WORLD,IERR)
            N6 = 6
 1950    CONTINUE
C ... KE VARIABLES
         IF(ITURB >= 3 .AND. ITURB /= 8) THEN
            DO 1960 N5 = N6+1,N6+2
               IF(IFORM == 0) THEN
                  READ(LAITE)   (ZZZ(N4),N4=1,NTOT)
               ELSE
                  READ(LAITE,*) (ZZZ(N4),N4=1,NTOT)
               ENDIF
               CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +              MPI_COMM_WORLD,IERR)
 1960       CONTINUE
            N6 = N6 + 2
         ENDIF
C ... Intermittency variables
         IF(TRANSL) THEN
            DO 1965 N5 = N6+1,N6+2
               IF(IFORM == 0) THEN
                  READ(LAITE)   (ZZZ(N4),N4=1,NTOT)
               ELSE
                  READ(LAITE,*) (ZZZ(N4),N4=1,NTOT)
               ENDIF
               CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +              MPI_COMM_WORLD,IERR)
 1965       CONTINUE
            N6 = N6 + 2
         ENDIF ! TRANSL
C ... SCALARS
         IF(NSCAL >= 1) THEN
            DO 1970 N5 = 1,NSCAL
               IF(IFORM == 0) THEN
                  READ(LAITE)   (ZZZ(N4),N4=1,NTOT)
               ELSE
                  READ(LAITE,*) (ZZZ(N4),N4=1,NTOT)
               ENDIF
               CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5+N6,
     +              MPI_COMM_WORLD,IERR)
 1970       CONTINUE
         ENDIF
C ...MULTIPHASE
         IF(MULPHL) THEN
            N6   = NTG+NSCAL+N6
            NGL  = NPROCE(LB+1,NP)
            IMUL = 2
            IF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') IMUL = 5
            DO 1990 IPHASE = 1,BLKS(NGL)%NPHASE
            DO 1980 N5 = 1,IMUL
               N6 = N6 + 1
               IF(IFORM == 0) THEN
                  READ(LAITE)   (ZZZ(N4),N4=1,NTOT)
               ELSE
                  READ(LAITE,*) (ZZZ(N4),N4=1,NTOT)
               ENDIF
               CALL MPI_SEND(ZZZ(1),NTOT,MPI_REAL8,NP-1,N6,
     +              MPI_COMM_WORLD,IERR)
 1980       CONTINUE
 1990       CONTINUE
         ENDIF
         GOTO 1000
      ENDIF                     !   (NP == 1):ELSE
      ENDIF                     !   (N == NPROCE(LB+1,NP))
 800  CONTINUE
 900  CONTINUE
      WRITE(*,*) 'No hits in PROSEC. Something wrong in COMRMP'
      WRITE(*,*) 'Contact Patrik Rautaheimo'
      WRITE(*,*) 'Exiting ...'
      CALL MPI_ABORT(MPI_COMM_WORLD,IERRORCODE,IERR)
      STOP
 1000 CONTINUE

      RETURN
      END SUBROUTINE COMRMP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COMSLA(NBLOCK,IG,IMAX,JMAX,KMAX,IN,JN,KN,MGM,
     +    F1R,F1RM,F1RN,F1RW,F1E,HAT1,F1RK,F1EPS,F1FI,MAXSB,NSCAL,ITURB,
     +    IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,VAR,IFORM,
     +    LAITE,MULPHL,BLKS,HAT2,TRANSL,TRM)
      
C ... RECEIVE COMPUT FILE FROM MASTER PROCES OR READ COMPUT IN ONE PROCES
C ... MODE

      USE MPI
      USE TYPE_ARRAYS
      USE NS3CO, ONLY : PARALLEL

      IMPLICIT NONE

      INTEGER :: NBLOCK,IN,JN,KN,MGM,MAXSB,NSCAL,ITURB,IPRO,
     +   MBPRO,NPRO,NBLOCG,IFORM,LAITE,N,IG1,IG2,KMAX2,
     +   NTOT,NTO2,IGN,NG,NTG,IERR,NS,IS1,N6,IPHASE
      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),IG(MGM,*)
      INTEGER :: NPROCE(MBPRO+1,*)
      TYPE (MPHASE_VARIABLES) :: VAR(*)
      TYPE (BLOCKS)           :: BLKS(*)
      TYPE (INTERMITTENCY)    :: TRM(*)
      LOGICAL :: MULPHL, TRANSL
      REAL    :: ZZZ(*),F1R(*),F1RM(*),F1RN(*),F1RW(*),F1E(*),F1RK(*),
     +     F1EPS(*),F1FI(*),HAT1(*),HAT2(*)

C ... BLOCK LOOP 1 BEGINS
      DO 2000 N = 1,NBLOCK
      IG2     = IG(2,N)
      IG1     = IG(1,N)
      KMAX2 = MAX0(1,KMAX(2,N))
         NTOT = (IMAX(2,N)+2*IN)*(JMAX(2,N)+2*JN)*(KMAX2+2*KN)
         NTO2 = (IMAX(2,N)+2)*(JMAX(2,N)+2)*(KMAX2+2)
         IGN    = IG1 + NTOT - 1
      IF(PARALLEL) THEN
         NG  = NPROCE(N+1,IPRO)         ! GLOBAL BLOCK NUMBER
         NTG  = NG*100 + 570            ! TAG NUMBER     
         CALL MPI_SEND(NTO2,1,MPI_INTEGER,0,NTG,MPI_COMM_WORLD,
     +        IERR)
      ELSE
         NG  = N
         NTG = 0
      ENDIF

      CALL READMP(LAITE,IFORM, F1R(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,1)
      CALL READMP(LAITE,IFORM,F1RM(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,2)
      CALL READMP(LAITE,IFORM,F1RN(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,3)
      CALL READMP(LAITE,IFORM,F1RW(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,4)
      CALL READMP(LAITE,IFORM, F1E(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,5)
      CALL READMP(LAITE,IFORM,HAT1(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,6)
      N6 = 6
      
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM, F1RK(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +        KMAX2,IN,JN,KN,IPRO,NTG,N6)
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,F1EPS(IG1),ZZZ,IMAX(2,N),JMAX(2,N),
     +        KMAX2,IN,JN,KN,IPRO,NTG,N6)
      ENDIF

C ... Intermittency variables
      IF(TRANSL) THEN
         N6   = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         TRM(IG1:IGN)%G = HAT2(IG1:IGN)
         N6   = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         TRM(IG1:IGN)%RET = HAT2(IG1:IGN)
      ENDIF ! TRANSL
C ... FOR SCALAR EQ. PPR 1.3
      DO 1900 NS = 1,NSCAL
         IS1 = IG1 + (NS-1)*MAXSB
         CALL READMP(LAITE,IFORM,F1FI(IS1),ZZZ,IMAX(2,N),JMAX(2,N),
     +     KMAX2,IN,JN,KN,IPRO,NTG,N6+NS)
 1900 CONTINUE
C ... MULTIPHASE
      IF(MULPHL) THEN
         N6 = N6+NSCAL !NS
         DO 1990 IPHASE = 1,BLKS(NG)%NPHASE
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         VAR(IG1:IGN)%F1R(IPHASE) = HAT2(IG1:IGN) ! Temperature
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         VAR(IG1:IGN)%F2R(IPHASE) = HAT2(IG1:IGN) ! Mass fraction
         IF(BLKS(NG)%SOLUTION_TYPE == 'MULTI') THEN
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         VAR(IG1:IGN)%FRM(IPHASE) = HAT2(IG1:IGN) ! U-velocity
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         VAR(IG1:IGN)%FRN(IPHASE) = HAT2(IG1:IGN) ! V-velocity
         N6 = N6 + 1
         CALL READMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(2,N),
     +     JMAX(2,N),KMAX2,IN,JN,KN,IPRO,NTG,N6)
         VAR(IG1:IGN)%FRW(IPHASE) = HAT2(IG1:IGN) ! W-velocity
         ENDIF ! MULTI
 1990    CONTINUE
      ENDIF ! MULPHL
 2000 CONTINUE

      RETURN
      END SUBROUTINE COMSLA


C     ***************************************************************
C     *      CONVERGENCE HISTORY SUBROUTINES                        *
C     ***************************************************************

      SUBROUTINE WHMON(ICYCLE,DWMAX,SUMDWH,IPRO,NPROCE,
     + MBPRO,NPRO,CONVL)
C ... WRITE THE CONVERGENCE MONITORING DATA TO TERMINAL AND FILES.
C ... NOW ALSO MULTIPROCESSING IS AVAILABLE BY MPI

      USE MPI

      USE NS3CO, ONLY : MAXPRO, PARALLEL

      LOGICAL :: CONVL
      INTEGER :: ERRORCODE
      INTEGER :: NPROCE(MBPRO+1,*),IAA(6*MAXPRO)
      REAL    :: XAA(3*MAXPRO), XAAW(3*MAXPRO)

      IF(NPRO >= MAXPRO) THEN
         WRITE(*,*) 'MAXPRO must be changed from',MAXPRO,' to at least',
     +        NPRO
         WRITE(*,*) 'in subroutine WHMON'
         WRITE(*,*) 'Aborting'
         STOP
      ENDIF
      DROL = DWMAX
      IF(PARALLEL .AND. CONVL) THEN
         IAA    = 0 ! For a later use?
         XAA(1) = DROL
         XAA(2) = SUMDWH
         XAA(3) = 0.
         XAAW(1) = XAA(1) 
         XAAW(2) = XAA(2) 
         XAAW(3) = XAA(3) 
         CALL MPI_GATHER(XAAW(1),3,MPI_REAL8,XAA(1),3,MPI_REAL8,0,
     +     MPI_COMM_WORLD,IERR)
         IF(IPRO /= 1) GOTO 1100
         DO 100 NP = 2,NPRO
            IND = (NP-1)*6
            INR = (NP-1)*3
            IF(IAA(IND+1) /= ICYCLE) THEN
               WRITE(*,*) 'Processes are not synchronoused.'
               WRITE(*,*) 'ICYCLE is ',ICYCLE,' in master and'
               WRITE(*,*) 'ICYCLE is ',IAA(IND+1),' in slave no:',NP
               WRITE(*,*) 'Exiting ...'
               CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
               STOP
            ENDIF
            IF(ABS(DROL) < ABS(XAA(INR+1))) THEN
               DROL   = XAA(INR+1)
               IDROMX = IAA(IND+3)
               JDROMX = IAA(IND+4)
               KDROMX = IAA(IND+5)
            ENDIF
            SUMDWH    = SUMDWH + XAA(INR+2)
 100     CONTINUE
      ENDIF ! PARALLEL .AND. CONVL

 1100 CONTINUE
      IF(PARALLEL) THEN
         CALL MPI_BCAST(DROL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         DWMAX = DROL
      ENDIF

 1000 FORMAT(1X,I6,E15.5,2X,I4,I6,I4,I4,1X,I8,3X,F10.6,3X,F10.6)
 2000 FORMAT(1X,I6,E15.5,2X,I4,I6,I4,I4,1X,I8,1X,E12.5,3X,F10.6)
 3000 FORMAT(1X,I6,E15.5,2X,I4,I6,I4,I4,1X,I8,3X,F10.6,1X,E12.5)
 4000 FORMAT(1X,I6,E15.5,2X,I4,I6,I4,I4,1X,I8,1X,E12.5,1X,E12.5)

      RETURN
      END SUBROUTINE WHMON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONMON(ICYCLE,DROMAX,CL,CD,CM,IXERR,INERR,IMACH,IMAX,
     + JMAX,IPRO,NPROCE,MBPRO,NPRO,CONVL)

C ... WRITE THE CONVERGENCE MONITORING DATA TO TERMINAL AND FILES.
C ... NOW ALSO MULTIPROCESSING IS AVAILABLE BY MPI

      USE MPI
      USE NS3CO, ONLY : IN, JN, KN, MAXPRO, PARALLEL

      LOGICAL :: CONVL
      INTEGER :: ERRORCODE
      INTEGER :: NPROCE(MBPRO+1,*),IAA(6*MAXPRO),IAAW(6*MAXPRO)
      REAL :: XAA(3*MAXPRO), XAAW(3*MAXPRO)

      IF(NPRO >= MAXPRO) THEN
         WRITE(*,*) 'MAXPRO must be changed from',MAXPRO,' to at least',
     +        NPRO
         WRITE(*,*) 'in subroutine CONMON'
         WRITE(*,*) 'Aborting'
         STOP
      ENDIF

      IM  = IMACH ! Global number of supersonic cells
      CDL = CD
      CLL = CL
      INER2 = NPROCE(INERR+1,IPRO)
      DROL = DROMAX

      IAPU=(IMAX+2*IN)
      JAPU=(JMAX+2*JN)*IAPU
      KDROMX=IXERR/JAPU+1-KN
      JDROMX=(IXERR-(KDROMX-1+KN)*JAPU)/IAPU+1-JN
      IDROMX=IXERR-(KDROMX-1+KN)*JAPU-(JDROMX-1+JN)*IAPU-IN
      IF(IDROMX == -2) THEN
           IDROMX = IMAX + IN
           JDROMX = JDROMX - 1
      ENDIF

      ABSCD = ABS(CDL)
      ABSCL = ABS(CLL)
      IF(ABSCD < 100. .AND. ABSCL < 100.) THEN
      WRITE(7,1000) ICYCLE,DROL,INERR,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD >= 100. .AND. ABSCL < 100.) THEN
      WRITE(7,2000) ICYCLE,DROL,INERR,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD < 100. .AND. ABSCL >= 100.) THEN
      WRITE(7,3000) ICYCLE,DROL,INERR,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE
      WRITE(7,4000) ICYCLE,DROL,INERR,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ENDIF
      IF(PARALLEL .AND. CONVL) THEN
         IAA(1) = ICYCLE
         IAA(2) = NPROCE(INERR+1,IPRO)
         IAA(3) = IDROMX
         IAA(4) = JDROMX
         IAA(5) = KDROMX
         IAA(6) = IM
         XAA(1) = DROL
         XAA(2) = CDL
         XAA(3) = CLL
         IAAW(1) = IAA(1)
         IAAW(2) = IAA(2)
         IAAW(3) = IAA(3)
         IAAW(4) = IAA(4)
         IAAW(5) = IAA(5)
         IAAW(6) = IAA(6)
         XAAW(1) = XAA(1)
         XAAW(2) = XAA(2)
         XAAW(3) = XAA(3)
         CALL MPI_GATHER(IAAW(1),6,MPI_INTEGER,IAA(1),6,MPI_INTEGER,0,
     +     MPI_COMM_WORLD,IERR) 
         CALL MPI_GATHER(XAAW(1),3,MPI_REAL8,XAA(1),3,MPI_REAL8,0,
     +     MPI_COMM_WORLD,IERR)
         IF(IPRO /= 1) GOTO 1100
         DO 100 NP = 2,NPRO
            IND = (NP-1)*6
            INR = (NP-1)*3
            IF(IAA(IND+1) /= ICYCLE) THEN
               WRITE(*,*) 'Processes are not synchronoused.'
               WRITE(*,*) 'ICYCLE is ',ICYCLE,' in master and'
               WRITE(*,*) 'ICYCLE is ',IAA(IND+1),' in slave no:',NP
               WRITE(*,*) 'Exiting ...'
               CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
               STOP
            ENDIF
            IF(ABS(DROL) < ABS(XAA(INR+1))) THEN
               DROL = XAA(INR+1)
               IDROMX = IAA(IND+3)
               JDROMX = IAA(IND+4)
               KDROMX = IAA(IND+5)
               INER2 = IAA(IND+2)
            ENDIF
            CDL    = CDL + XAA(INR+2)
            CLL    = CLL + XAA(INR+3)
            IM = IM + IAA(IND+6)
 100     CONTINUE
         ABSCD = ABS(CDL)
         ABSCL = ABS(CLL)
      IF(ABSCD < 100. .AND. ABSCL < 100.) THEN
      WRITE(11,1000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD >= 100. .AND. ABSCL < 100.) THEN
      WRITE(11,2000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD < 100. .AND. ABSCL >= 100.) THEN
      WRITE(11,3000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE
      WRITE(11,4000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ENDIF
      ENDIF ! PARALLEL .AND. CONVL
      IF (IPRO /= 1) GOTO 1100
      ABSCD = ABS(CDL)
      ABSCL = ABS(CLL)
         
      IF(ABSCD < 100. .AND. ABSCL < 100.) THEN
      WRITE(6,1000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD >= 100. .AND. ABSCL < 100.) THEN
      WRITE(6,2000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE IF(ABSCD < 100. .AND. ABSCL >= 100.) THEN
      WRITE(6,3000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ELSE
      WRITE(6,4000) ICYCLE,DROL,INER2,IDROMX,JDROMX,KDROMX,IM,CDL,CLL
      ENDIF
 
 1100 CONTINUE
      IF(PARALLEL) THEN
         CALL MPI_BCAST(DROL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         DROMAX = DROL
      ENDIF

 1000 FORMAT(1X,I7,E15.5,2X,I4,I6,I4,I4,1X,I8,3X,F10.6,3X,F10.6)
 2000 FORMAT(1X,I7,E15.5,2X,I4,I6,I4,I4,1X,I8,1X,E12.5,3X,F10.6)
 3000 FORMAT(1X,I7,E15.5,2X,I4,I6,I4,I4,1X,I8,3X,F10.6,1X,E12.5)
 4000 FORMAT(1X,I7,E15.5,2X,I4,I6,I4,I4,1X,I8,1X,E12.5,1X,E12.5)

      RETURN
      END SUBROUTINE CONMON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IBDWRI(LAITE,IBOUT,ICYCTI,T,CLDRO,RESI,FILEE,
     +     PRN,CD,CL,CS,CMX,CMY,CMZ,IMACH,NPROCE,MBPRO,NPRO,
     +     NTO3,IPRO,NRESI,VOLTOT,CONVL,QMFIN,QMFOUT,QMEIN,QMEOUT,
     +     HEADT,EFFIC,TOM,RESM,IESM,QVFIN,QVFOUT,DWMAX,FLODWH,
     +     WHMAX,WHMIN,SUMDWH,SINK,TRIM,QGFIN,QGFOUT,QGEIN,QGEOUT)

C ... WRITE THE CONVERGENCE MONITORING DATA TO TERMINAL AND FILES.
C ... NOW ALSO MULTIPROCESSING IS AVAILABLE BY MPI

      USE MPI
      USE CONSTANTS, ONLY : RAD2DEG
      USE NS3CO, ONLY : G0, PARALLEL

      CHARACTER(LEN=6) :: FILEE
      CHARACTER(LEN=3) :: PRN

      REAL    :: RESI(NRESI,*), RESM(*)      ! *ESM-files have a large 
      REAL    :: SINK(*), TRIM(*)            ! buffer... (ZZZ/2)
      REAL    :: RESMW(1000),HEADTT

      INTEGER :: ERRORCODE
      INTEGER :: NPROCE(MBPRO+1,*),IESM(*)   ! Whatch the size of this

      INTEGER :: IESMW(1000)
      INTEGER :: IRESI
      LOGICAL :: CONVL

      CLDROL = CLDRO  
                       ! Paivitetaan IBDIAG joka kolmas kierros (esim.)
      IBOUTP = IBOUT   ! Luetaan INPUTista tai sijoitus ohjelman alussa

      IF(.NOT.PARALLEL) THEN
         QMEIN2     = QMEIN /(QMFIN +1.E-20)
         QMEOU2     = QMEOUT/(QMFOUT+1.E-20)
         RESI(51,1) = QGFOUT
         RESI(52,1) = QGEOUT
         RESI(53,1) = SINK(1)
         RESI(54,1) = TRIM(1)*RAD2DEG
         RESI(55,1) = SUMDWH
         RESI(56,1) = WHMAX
         RESI(57,1) = WHMIN
         RESI(58,1) = ABS(DWMAX)
         RESI(59,1) = FLODWH
         RESI(60,1) = QVFIN
         RESI(61,1) = QVFOUT

*         WRITE(LAITE) ICYCTI,T,CD,CL,CS,CMX,CMY,CMZ,IMACH,HEADT,TOM,
*     +      QMFIN,QMFOUT,QMEIN2,QMEOU2,CLDRO,(RESI(JJ,1),JJ=1,NRESI)

         WRITE(LAITE) ICYCTI,REAL(T,4),REAL(CD,4),REAL(CL,4),
     +        REAL(CS,4),REAL(CMX,4),REAL(CMY,4),REAL(CMZ,4),IMACH,
     +        REAL(HEADT,4),REAL(TOM,4),REAL(QMFIN,4),REAL(QMFOUT,4),
     +        REAL(QMEIN2,4),REAL(QMEOU2,4),REAL(CLDRO,4),
     +       (REAL(RESI(JJ,1),4),JJ=1,NRESI)

         IF(CONVL) THEN
            IF (ICYCTI >= 100)  IBOUTP = 10*IBOUT
            IF (ICYCTI >= 1000) IBOUTP = 1000000
            IF(MOD(ICYCTI,IBOUTP) == 0) THEN
               CLOSE(LAITE)
               OPEN(LAITE,FILE=FILEE//PRN,STATUS='UNKNOWN',
     +              FORM='UNFORMATTED')
               CALL WTC(LAITE,ICYCTI)
            ENDIF
         ENDIF
      ENDIF

      IF(PARALLEL .AND. CONVL) THEN
         DO 50 JJ = 1,NRESI
            RESM(JJ) = RESI(JJ,1)
 50      CONTINUE

C ... AVERAGED  VALUES 
         DO 55 JJ = 22,23
            RESM(JJ) = RESI(JJ,1)*VOLTOT
 55      CONTINUE

         RESM(51) = QGFOUT !Gas mass flux out
         RESM(52) = QGEOUT !Gas energy flux out 
         RESM(53) = SINK(1)!Sink of the first ship
         RESM(54) = TRIM(1)*RAD2DEG!Trim of the first ship
         RESM(55) = SUMDWH !L1-norm of the wave height residual
         RESM(56) = WHMAX  !Maximum wave height
         RESM(57) = WHMIN  !Minimum wave height
         RESM(58) = ABS(DWMAX)  !Maximum change in wave height
         RESM(59) = FLODWH !Surface total volumetric flow
         RESM(60) = QVFIN  !Volume flux in
         RESM(61) = QVFOUT !Volume flux out 
C ... End of mersu
         RESM(NRESI+1)  = T      !TIME LEVEL
         RESM(NRESI+2)  = CLDRO  !MAXIMUN CHANGE        
         RESM(NRESI+3)  = CD     !DRAG COEFFICIENT
         RESM(NRESI+4)  = CL     !LIFT COEFFICIENT
         RESM(NRESI+5)  = CS     !SIDE FORCE COEFFICIENT
         RESM(NRESI+6)  = CMX    !MOMENTUM COEFFICIENT
         RESM(NRESI+7)  = CMY    !MOMENTUM COEFFICIENT
         RESM(NRESI+8)  = CMZ    !MOMENTUM COEFFICIENT
         RESM(NRESI+9)  = VOLTOT !TOTAL VOLUME OF PROCESS
         RESM(NRESI+10) = QMFIN  !TOTAL INLET FLOW OF PROCESS
         RESM(NRESI+11) = QMFOUT !TOTAL OUTLET FLOW OF PROCESS
         RESM(NRESI+12) = QMEIN /(QMFIN +1.E-20) !TOTAL INLET ENERGY FLOW 
         RESM(NRESI+13) = QMEOUT/(QMFOUT+1.E-20) !TOTAL OUTLET ENERGY FLOW 
         RESM(NRESI+14) = TOM    !TORQUE x OMEGA = P
         IESM(1) = ICYCTI        !CYCLE
         IESM(2) = IMACH         !NUMBER OF SUP. CELLS
         IESM(3) = NTO3          !TOTAL NUMBER OF CELLS IN ONE PROCESS
         IESMW(1) = IESM(1)
         IESMW(2) = IESM(2)
         IESMW(3) = IESM(3)
         DO IRESI=1,NRESI+14
            RESMW(IRESI) = RESM(IRESI)
         ENDDO

         CALL MPI_GATHER(IESMW,3,MPI_INTEGER,IESM,3,MPI_INTEGER,0,
     +        MPI_COMM_WORLD,IERR) 
         CALL MPI_GATHER(RESMW,NRESI+14,MPI_REAL8,RESM,NRESI+14,
     +        MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(HEADTT,HEADT,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)

         IF(IPRO /= 1) GOTO 1100

         DO JJ = 1,21
            RESM(JJ) = (RESM(JJ)*REAL(NTO3))**2
         ENDDO ! Mersu
         DO JJ = 31,37
            RESM(JJ) = (RESM(JJ)*REAL(NTO3))**2
         ENDDO
            RESM(39) = (RESM(39)*REAL(NTO3))**2
         DO JJ = 41,45
            RESM(JJ) = (RESM(JJ)*REAL(NTO3))**2
         ENDDO
         DO JJ = 49,50
            RESM(JJ) = (RESM(JJ)*REAL(NTO3))**2
         ENDDO
         DO JJ = 62,67
            RESM(JJ) = (RESM(JJ)*REAL(NTO3))**2
         ENDDO


         DO 100 NP = 2,NPRO
            IND = (NP-1)*3
            INR = (NP-1)*(NRESI+14)
            IF(IESM(IND+1) /= ICYCTI) THEN
               WRITE(*,*) 'Processes are not synchronoused.'
               WRITE(*,*) 'ICYCTI is ',ICYCTI,' in master and'
               WRITE(*,*) 'ICYCTI is ',IESM(IND+1),' in slave no:',NP
               WRITE(*,*) 'Exiting ...'
               CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
               STOP
            ENDIF
C ... SUM UP THE L2 NORMS
            DO 60 JJ = 1,21
               RESM(JJ) = RESM(JJ) + (RESM(JJ+INR)*REAL(IESM(IND+3)))**2
 60         CONTINUE
C ... SUM UP TOTAL AND AVERAGED VALUES
            DO 65 JJ = 22,30 
               RESM(JJ) = RESM(JJ) + RESM(JJ+INR)
 65         CONTINUE
            DO JJ = 46,48 
               RESM(JJ) = RESM(JJ) + RESM(JJ+INR)
            ENDDO
            RESM(38)    = RESM(38) + RESM(38+INR) ! Viscous dissipation
            RESM(40)    = RESM(40) + RESM(40+INR) ! LES kinetic energy
C ... SUM UP THE L2 NORMS
            DO 70 JJ = 31,37
               RESM(JJ) = RESM(JJ) + (RESM(JJ+INR)*REAL(IESM(IND+3)))**2
 70         CONTINUE
C ... Two-phase flows
               RESM(39) = RESM(39) + (RESM(39+INR)*REAL(IESM(IND+3)))**2
            DO JJ = 41,45
               RESM(JJ) = RESM(JJ) + (RESM(JJ+INR)*REAL(IESM(IND+3)))**2
            ENDDO
C ... Re_t transition model
            DO 71 JJ = 49,50
               RESM(JJ) = RESM(JJ) + (RESM(JJ+INR)*REAL(IESM(IND+3)))**2
 71         CONTINUE
C ... Two-fluid momentum residuals
            DO JJ = 62,67
               RESM(JJ) = RESM(JJ) + (RESM(JJ+INR)*REAL(IESM(IND+3)))**2
            ENDDO

C ... SUM UP TOTAL AND AVERAGED VALUES
*           DO 75 JJ = NRESI-2,NRESI !oletus etta loput resissa ovat naita.
            DO 75 JJ = 59,61
               RESM(JJ) = RESM(JJ) + RESM(JJ+INR)
 75         CONTINUE
            RESM(55) = RESM(55) + RESM(55+INR)
C ... FIND MAXIMUMs
            RESM(NRESI+2) = MAX(RESM(NRESI+2),RESM(NRESI+2+INR)) ! CLDRO
            RESM(58)      = MAX(RESM(58),RESM(58+INR))           ! DWH
            RESM(56)      = MAX(RESM(56),RESM(56+INR))           ! WH
C ... AND MINIMUMs
            RESM(57)      = MIN(RESM(57),RESM(57+INR))           ! WH
C ... SUM UP FORCE COEFFICIENTS AND INLET&OUTLET FLOWS
            DO I4 = 3,14
               RESM(NRESI+I4) = RESM(NRESI+I4) + RESM(NRESI+I4+INR) 
            ENDDO

            IF(FILEE == 'IBDIAG') THEN
            IESM(2)  = IESM(2)  + IESM(2+IND) !NUMBER OF SUP. CELLS
            ENDIF
            IESM(3)  = IESM(3)  + IESM(3+IND) !TOTAL NUMBER OF CELLS
 100     CONTINUE
C ... NORMALIZED L2-NORMS
         DO 105 JJ = 1,21
            RESM(JJ) = SQRT(RESM(JJ))/REAL(IESM(3))
 105      CONTINUE
         DO 106 JJ = 31,37
            RESM(JJ) = SQRT(RESM(JJ))/REAL(IESM(3))
 106      CONTINUE
            RESM(39) = SQRT(RESM(39))/REAL(IESM(3))
         DO JJ = 41,45
            RESM(JJ) = SQRT(RESM(JJ))/REAL(IESM(3))
         ENDDO
         DO 107 JJ = 49,50
            RESM(JJ) = SQRT(RESM(JJ))/REAL(IESM(3))
 107      CONTINUE

         DO JJ = 62,67
            RESM(JJ) = SQRT(RESM(JJ))/REAL(IESM(3))
         ENDDO
C ... Remove tiny values
         IF(RESM(62) <= 2.E-8) RESM(62) = RESM(63)
         IF(RESM(64) <= 2.E-8) RESM(64) = RESM(65)
         IF(RESM(66) <= 2.E-8) RESM(66) = RESM(67)

C ... AVERAGED VALUES
      DO 8197 JJ = 22,23
         RESM(JJ)= RESM(JJ)/RESM(NRESI+9)
 8197 CONTINUE

cc      QMEIN2    = RESM(12) /(RESM(10)+1.E-20)
cc      QMEOU2    = RESM(13) /(RESM(11)+1.E-20)
cc      EFFICT    = 100.*ABS(RESM(10)*(QMEOU2-QMEIN2)/(RESM(14)+1.E-20))
c      HEADTT    = (RESM(13)-RESM(12))/(RESM(10)+1.E-20)/G0  ! ADD HOC
cc      HEADTT    = (RESM(13)/(RESM(11)+1.E-20) -
cc     +             RESM(12)/(RESM(10)+1.E-20) )/G0  ! ADD HOC
                       ! Paivitetaan IBDIAG joka kolmas kierros (esim.)
      IBOUTP = IBOUT   ! Luetaan INPUTista tai sijoitus ohjelman alussa
      IF (ICYCTI >= 100)  IBOUTP = 10*IBOUT
      IF (ICYCTI >= 1000) IBOUTP = 1000000
cc      WRITE(LAITE) ICYCTI,RESM(NRESI+1),(RESM(JJ),JJ=NRESI+3,NRESI+8),
cc     +     IESM(2),HEADTT,RESM(NRESI+14),(RESM(JJ),JJ=NRESI+10,
cc     +     NRESI+13),RESM(NRESI+2),(RESM(JJ),JJ=1,NRESI)
! pistetään kirjoituskomento näin 16.5.2012 JIL
*      WRITE(LAITE) ICYCTI,RESM(NRESI+1),(RESM(JJ),JJ=NRESI+3,NRESI+8),
*     +     IESM(2),HEADTT,RESM(NRESI+14),(RESM(JJ),JJ=NRESI+10,
*     +     NRESI+11),(RESMW(JJ),JJ=NRESI+12,NRESI+13),
*     +     RESM(NRESI+2),(RESM(JJ),JJ=1,NRESI)

      WRITE(LAITE) ICYCTI,REAL(RESM(NRESI+1),4),
     +     (REAL(RESM(JJ),4),JJ=NRESI+3,NRESI+8),
     +      IESM(2),REAL(HEADTT,4),REAL(RESM(NRESI+14),4),
     +     (REAL(RESM(JJ),4),JJ=NRESI+10,NRESI+11),
     +     (REAL(RESMW(JJ),4),JJ=NRESI+12,NRESI+13),
     +      REAL(RESM(NRESI+2),4),(REAL(RESM(JJ),4),JJ=1,NRESI)

      IF(MOD(ICYCTI,IBOUTP) == 0) THEN
        CLOSE(LAITE)
        OPEN(LAITE,FILE=FILEE,STATUS='UNKNOWN',FORM='UNFORMATTED')
        CALL WTC(LAITE,ICYCTI)
      ENDIF
      ENDIF                     !(PARALLEL .AND. CONVL)

 1100 CONTINUE

      RETURN
      END SUBROUTINE IBDWRI

C     ***************************************************************
C     *      END OF THE ITERATION                                   *
C     ***************************************************************

      SUBROUTINE PRIRES(RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS,CD,CL,CMZ,
     +     IPRO,NPRO,WCALI,WCOMI)

C ... PRINT THE RESULT ON THE SCREEN

      USE MPI

      USE NS3CO, ONLY : MAXPRO, PARALLEL
      
      REAL    :: XAA(8*MAXPRO), XAAW(8*MAXPRO)
      INTEGER :: IAA(MAXPRO), IAAW(MAXPRO)

      IF(PARALLEL) THEN
         IF(NPRO >= MAXPRO) THEN
         WRITE(*,*) 'MAXPRO must be changed from',MAXPRO,' to at least',
     +        NPRO
         WRITE(*,*) 'in subroutine PRIRES'
         WRITE(*,*) 'Aborting'
         STOP
         ENDIF
         IAA(1) = NCELLS
         XAA(1) = RTIME
         XAA(2) = RTIMEI
         XAA(3) = RTIMEC
         XAA(4) = CD
         XAA(5) = CL
         XAA(6) = CMZ
         XAA(7) = WCALI
         XAA(8) = WCOMI
         IAAW(1) = IAA(1)
         XAAW(1) = XAA(1)
         XAAW(2) = XAA(2)
         XAAW(3) = XAA(3)
         XAAW(4) = XAA(4)
         XAAW(5) = XAA(5)
         XAAW(6) = XAA(6)
         XAAW(7) = XAA(7)
         XAAW(8) = XAA(8)
         CALL MPI_GATHER(IAAW(1),1,MPI_INTEGER,IAA(1),1,MPI_INTEGER,0,
     +     MPI_COMM_WORLD,IERR) 
         CALL MPI_GATHER(XAAW(1),8,MPI_REAL8,XAA(1),8,MPI_REAL8,0,
     +     MPI_COMM_WORLD,IERR) 
         IF(IPRO /= 1) GOTO 1300
         RTG = 0.
         RIG = 0.
         RCG = 0.
         CDG = 0.
         CLG = 0.
         CMG = 0.
         WCAG = 0.
         WCOG = 0.
         ICE = 0
         WRITE(*,*)
         WRITE(13,*)
         WRITE(*,1001)
         WRITE(13,1001)
         WRITE(4,*)
         WRITE(4,1000)
         WRITE(45,*)
         WRITE(45,1001)
 1000    FORMAT
     + ('Process   CPU-time/cycle  Number of Cells'/
     +  '     CPU-time   CPU-time/cycle/cell(mus)  CD  ',
     +  '        CL             CMZ'/
     +  '=============================================================='
     +  ,'================')
 1001    FORMAT
     + ('Process   CPU-time/cycle   Number of Cells'/
     +  '     CPU-time  CPU-time /cycle/cell(mus)  Wtime',
     +  '/c        Wcom/c         Wcom/Wtim'/					!
     +  '=============================================================='
     +  ,'================')
         DO 200 NP = 1,NPRO
            INR   = (NP-1)*8
            RTL   = XAA(1+INR)
            RIL   = XAA(2+INR)
            RCL   = XAA(3+INR)
            CDL   = XAA(4+INR)
            CLL   = XAA(5+INR)
            CML   = XAA(6+INR)
            WCALL = XAA(7+INR)
            WCOML = XAA(8+INR)
            RTG   = RTG + RTL
            RIG   = RIG + RIL
            RCG   = RCG + RCL
            CDG   = CDG + CDL
            CLG   = CLG + CLL
            CMG   = CMG + CML
            WCAG  = WCAG + WCALL
            WCOG  = WCOG + WCOML
            ICE   = ICE + IAA(NP)
C ... PRINT FORCE COEFFICIENT IN SCREEN
C           WRITE(*,1100) NP,RTL,RIL,RCL,IAA(NP),CDL,CLL,CML
C ... PRINT MPI PERFORMANCE DATA IN SCREEN
            WRITE(*,1101) NP,RTL,RIL,RCL,IAA(NP),WCALL,WCOML,WCOML/WCALL
            WRITE(13,1101)NP,RTL,RIL,RCL,IAA(NP),WCALL,WCOML,WCOML/WCALL
            WRITE(45,1101)NP,RTL,RIL,RCL,IAA(NP),WCALL,WCOML,WCOML/WCALL
            WRITE(4,1100) NP,RTL,RIL,RCL,IAA(NP),CDL,CLL,CML
 200     CONTINUE
C ... PRINT FORCE COEFFICIENT IN SCREEN
C         WRITE(*,1150) RTG,RIG,RIG/ICE*1.E6,ICE,CDG,CLG,CMG
C ... PRINT MPI PERFORMANCE DATA IN SCREEN
         WRITE(*,1151) RTG,RIG,RIG/ICE*1.E6,ICE,WCAG,WCOG,WCOG/WCAG
         WRITE(13,1151)RTG,RIG,RIG/ICE*1.E6,ICE,WCAG,WCOG,WCOG/WCAG
         WRITE(6,*)
         WRITE(6,*) ' CD = ',REAL(CDG,4),
     +            '   CL = ',REAL(CLG,4),
     +           '   CMZ = ',REAL(CMG,4)
         WRITE(13,*)
         WRITE(13,*) ' CD = ',REAL(CDG,4),
     +             '   CL = ',REAL(CLG,4),
     +            '   CMZ = ',REAL(CMG,4)


 1100    FORMAT(I3,F9.2,F9.2,1X,F9.2,I8,1X,3G14.7)
 1101    FORMAT(I3,F9.2,F9.2,1X,F9.2,I8,1X,3(1X,G14.7))			!
 1150    FORMAT(
     +  '=============================================================='
     +  ,'================'/
     +  'Sum',F9.2,F9.2,1X,F9.2,I8,1X,3G14.7)
 1151    FORMAT(
     +  '=============================================================='
     +  ,'================'/
     +  'Sum',F9.2,F9.2,1X,F9.2,I8,1X,3(1X,G14.7))
      ELSE
         WRITE(6,*)
         IF(ICYTOT <= 999999) THEN
         WRITE(6,11) RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS
         ELSE
         WRITE(6,12) RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS
         ENDIF
         WRITE(6,*) ' CD = ',REAL(CD,4),
     +            '   CL = ',REAL(CL,4),
     +           '   CMZ = ',REAL(CMZ,4)
         WRITE(13,*)
         IF(ICYTOT <= 999999) THEN
         WRITE(13,11)RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS
         ELSE
         WRITE(13,12) RTIME,ICYTOT,RTIMEI,RTIMEC,NCELLS
         ENDIF
         WRITE(13,*)' CD = ',REAL(CD,4),
     +            '   CL = ',REAL(CL,4),
     +           '   CMZ = ',REAL(CMZ,4)
      ENDIF

 11   FORMAT(4X,'CPU-TIME                 =',F9.2,3X,'CYCLES =',I6/
     &     4X,'CPU-TIME/CYCLE           =',F9.2/
     &     4X,'CPU-TIME/CYCLE/CELL(mus) =',F9.2,3X,
     &     'NUMBER OF CELLS =',I8)
 12   FORMAT(4X,'CPU-TIME                 =',F9.2,3X,'CYCLES =',I8/
     &     4X,'CPU-TIME/CYCLE           =',F9.2/
     &     4X,'CPU-TIME/CYCLE/CELL(mus) =',F9.2,3X,
     &     'NUMBER OF CELLS =',I8)


 1300 CONTINUE

      RETURN
      END SUBROUTINE PRIRES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COMPIW(NBLOCK,IG,JF,IR,IMAX,JMAX,KMAX,IN,JN,KN,MGM,
     +        RO,RM,RN,RW,E,PDIFF,RK,REPS,FI,MAXSB,NSCAL,ITURB,
     +        IPRO,NPROCE,MBPRO,NPRO,NBLOCG,ZZZ,IFORM,LAITE,
     +        MULPHL,BLKS,HAT2,VAR,PRO,TRANSL,TRM)

C ... WRITE COMPUT FILE BY MASTER PROCES.
C ... THIS IS DONE SO THAT NO BUFFER IS NEEDED IN MPI

      USE MPI
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: ERRORCODE,IPRO,N,NBLOCG,NTG,NP,NPRO,LB,IF1,IG1,IR1,
     +           IGN,KMAX2,NTOT,NTOT2,IN,JN,KN,LAITE,IFORM,ITURB,NS,
     +           NSCAL,IPHASE,IERR,N4,N5,NG,MGM,MAXSB,MBPRO,NBLOCK,N6
      INTEGER :: NPROCE(MBPRO+1,*)
      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),JF(MGM,*),
     +           IG(MGM,*),IR(MGM,*)
      REAL :: RO(*),RM(*),RN(*),RW(*),E(*),RK(*),REPS(*),FI(MAXSB,*),
     +     ZZZ(*),PDIFF(*),HAT2(*)
      LOGICAL :: MULPHL,TRANSL
      TYPE(MPHASE_VARIABLES)  :: VAR(*)
      TYPE(PROPERTIES)        :: PRO(*)
      TYPE (BLOCKS)           :: BLKS(*)
      TYPE(INTERMITTENCY)     :: TRM(*)

C ... MASTER PROCES
      IF (IPRO == 1) THEN
         DO 1100 N = 1,NBLOCG        ! GLOBAL BLOCK LOOP (Miksi N)?
         NTG  = N*100 + 570          ! TAG NUMBER     
         DO 900 NP = 1,NPRO        ! PROCESSOR LOOP
         DO 800 LB=1,NPROCE(1,NP)    ! LOCAL BLOCK LOOP
            IF(N == NPROCE(LB+1,NP)) THEN
            IF(NP == 1) THEN       ! MASTER BLOCK
C ... THIS IS OWN BLOCK AND WILL BE WRITTEN IN ORDINARY WAY
               IG1   = IG(1,LB)
               IF1   = JF(1,LB)
               IR1   = IR(1,LB)
               KMAX2 = MAX(1,KMAX(1,LB))
               NTOT2 = (IMAX(1,LB)+2*IN)*(JMAX(1,LB)+2*JN)*(KMAX2+2*KN)
               IGN   = IG1 + NTOT2 - 1   

               CALL  WRITER(LAITE,IFORM,RO(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               CALL  WRITER(LAITE,IFORM,RM(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               CALL  WRITER(LAITE,IFORM,RN(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               CALL  WRITER(LAITE,IFORM,RW(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               CALL  WRITER(LAITE,IFORM, E(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               CALL  WRITER(LAITE,IFORM,PDIFF(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               IF(ITURB >= 3 .AND. ITURB /= 8) THEN
                  CALL  WRITER(LAITE,IFORM,RK(IG1),ZZZ,IMAX(1,LB),
     +                 JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                  CALL  WRITER(LAITE,IFORM,REPS(IG1),ZZZ,IMAX(1,LB),
     +                 JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               ENDIF
               IF(TRANSL) THEN
                     HAT2(IG1:IGN) = TRM(IG1:IGN)%G 
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     HAT2(IG1:IGN) = TRM(IG1:IGN)%RET 
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
               ENDIF ! TRANSL
               DO 1900 NS = 1,NSCAL
                  CALL WRITER(LAITE,IFORM,FI(IG1,NS),ZZZ,IMAX(1,LB),
     +                 JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
 1900          CONTINUE
               IF(MULPHL) THEN
                  DO IPHASE = 1,BLKS(N)%NPHASE
                     HAT2(IG1:IGN) = PRO(IG1:IGN)%DTEMP(IPHASE) 
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     HAT2(IG1:IGN) = VAR(IG1:IGN)%X(IPHASE) 
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     IF(BLKS(N)%SOLUTION_TYPE == 'MULTI') THEN
                     HAT2(IG1:IGN) = VAR(IG1:IGN)%U(IPHASE) ! u-velocity
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     HAT2(IG1:IGN) = VAR(IG1:IGN)%V(IPHASE) ! v-velocity
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     HAT2(IG1:IGN) = VAR(IG1:IGN)%W(IPHASE) ! w-velocity
                     CALL WRITER(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +                    JMAX(1,LB),KMAX(1,LB),IN,JN,KN)
                     ENDIF ! MULTI
                  ENDDO
               ENDIF ! MULPHL
               GOTO 1000
            ELSE                ! SLAVE BLOCK
C ... SLAVE BLOCK IS RECVEIDED FROM SLAVE AND WRITTEN IN COMPUT
               CALL MPI_RECV(NTOT,1,MPI_INTEGER,NP-1,NTG,
     +              MPI_COMM_WORLD,STATUS,IERR)
               DO 1950 N5 = 1,6
                  CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                 MPI_COMM_WORLD,STATUS,IERR)
                  IF(IFORM == 0) THEN
                     WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                  ELSE
                     WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                  ENDIF
 1950          CONTINUE
               N6 = 6
               IF(ITURB >= 3 .AND. ITURB /= 8) THEN
                  DO 1960 N5 = N6+1,N6+2
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
 1960             CONTINUE
               N6 = N6 + 2
               ENDIF
               IF(TRANSL) THEN ! Intermittency variables
                  DO 1965 N5 = N6+1,N6+2
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
 1965             CONTINUE
               N6 = N6 + 2
               ENDIF
               IF(NSCAL >= 1) THEN
                  DO 1970 N5 = 1,NSCAL
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5+N6,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
 1970             CONTINUE
                  N6 = N6 + NSCAL
               ENDIF
               IF(MULPHL) THEN
                  N5 = N6
                  DO IPHASE = 1,BLKS(N)%NPHASE
                     N5 = N5 + 1
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
                     N5 = N5 + 1
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
                     IF(BLKS(N)%SOLUTION_TYPE == 'MULTI') THEN
                     N5 = N5 + 1
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
                     N5 = N5 + 1
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
                     N5 = N5 + 1
                     CALL MPI_RECV(ZZZ(1),NTOT,MPI_REAL8,NP-1,NTG+N5,
     +                    MPI_COMM_WORLD,STATUS,IERR)
                     IF(IFORM == 0) THEN
                        WRITE(LAITE)   (ZZZ(N4),N4=1,NTOT)
                     ELSE
                        WRITE(LAITE,*) (ZZZ(N4),N4=1,NTOT)
                     ENDIF
                     ENDIF ! MULTI
                  ENDDO
                  N6 = N5 ! Antakee mersu
               ENDIF ! MULPHL   
               GOTO 1000
            ENDIF
            ENDIF   !   (N == NPROCE(LB+1,NP))
 800     CONTINUE
 900     CONTINUE
         WRITE(*,*) 'No hits in PROSEC. Something wrong'
         WRITE(*,*) 'Contact Patrik Rautaheimo'
         WRITE(*,*) 'Exiting ...'
         CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
 1000    CONTINUE
 1100    CONTINUE

      ELSE  !(IPRO == 1)
C ... SLAVE PROSESSES SEND THE DATA TO THE MASTER PROSES
         
         DO 2800 LB=1,NPROCE(1,IPRO)     ! LOCAL BLOCK LOOP
            NG    = NPROCE(LB+1,IPRO)    ! GLOBAL BLOCK NUMBER ! Mersu
            NTG   = NG*100 + 570         ! TAG NUMBER     
            IG1   = IG(1,LB)
            IF1   = JF(1,LB)
            IR1   = IR(1,LB)
            NTOT  = (IMAX(1,LB)+2)*(JMAX(1,LB)+2)*(KMAX(1,LB)+2)
            KMAX2 = MAX(1,KMAX(1,LB))
            NTOT2 = (IMAX(1,LB)+2*IN)*(JMAX(1,LB)+2*JN)*(KMAX2+2*KN)
            IGN   = IG1 + NTOT2 - 1   

            CALL MPI_SEND(NTOT,1,MPI_INTEGER,0,NTG,MPI_COMM_WORLD,
     +           IERR)
            CALL  WRITMP(LAITE,IFORM,RO(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+1)
            CALL  WRITMP(LAITE,IFORM,RM(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+2)
            CALL  WRITMP(LAITE,IFORM,RN(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+3)
            CALL  WRITMP(LAITE,IFORM,RW(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+4)
            CALL  WRITMP(LAITE,IFORM, E(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+5)
            CALL  WRITMP(LAITE,IFORM,PDIFF(IG1),ZZZ,IMAX(1,LB),
     +           JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+6)
            N6 = 6
            IF(ITURB >= 3 .AND. ITURB /= 8) THEN
               N6 = N6 + 1
               CALL  WRITMP(LAITE,IFORM,RK(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N6)
               N6 = N6 + 1
               CALL  WRITMP(LAITE,IFORM,REPS(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N6)
            ENDIF
            IF(TRANSL) THEN ! Intermittency variables
               N6 = N6 + 1
               HAT2(IG1:IGN) = TRM(IG1:IGN)%G
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N6)
               N6 = N6 + 1
               HAT2(IG1:IGN) = TRM(IG1:IGN)%RET
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N6)
            ENDIF
            DO 2900 NS = 1,NSCAL
               CALL WRITMP(LAITE,IFORM,FI(IG1,NS),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N6+NS)
 2900       CONTINUE
            IF(MULPHL) THEN
            N5 = NSCAL + N6
            DO IPHASE = 1,BLKS(NG)%NPHASE
               N5 = N5 + 1
               HAT2(IG1:IGN) = PRO(IG1:IGN)%TEMP(IPHASE)
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N5)
               N5 = N5 + 1
               HAT2(IG1:IGN) = VAR(IG1:IGN)%X(IPHASE)
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N5)
               IF(BLKS(NG)%SOLUTION_TYPE == 'MULTI') THEN
               N5 = N5 + 1
               HAT2(IG1:IGN) = VAR(IG1:IGN)%U(IPHASE)
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N5)
               N5 = N5 + 1
               HAT2(IG1:IGN) = VAR(IG1:IGN)%V(IPHASE)
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N5)
               N5 = N5 + 1
               HAT2(IG1:IGN) = VAR(IG1:IGN)%W(IPHASE)
               CALL  WRITMP(LAITE,IFORM,HAT2(IG1),ZZZ,IMAX(1,LB),
     +              JMAX(1,LB),KMAX(1,LB),IN,JN,KN,NTG+N5)
               ENDIF ! MULTI
            ENDDO
            N6 = N5 ! Mersu
            ENDIF ! MULPHL
 2800    CONTINUE
      ENDIF     !(ELSE .... IPRO == 1)

      RETURN
      END SUBROUTINE COMPIW
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITMP(LAITE,IFORM,RO,APUP2,IMAX,JMAX,KMAX,IN,JN,KN,
     +     NTG)

      USE MPI

      REAL :: RO(*),APUP2(*)

C ... small helping subroutine for COMPIW

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      N       = 0
      DO 3620 K = 0,KMAX + 1
      KK      = (KN+K-1)*ISTRID*JSTRID
      DO 3620 J = 0,JMAX + 1
      DO 3620 I = 0,IMAX + 1
         IJ      = (JN+J-1)*ISTRID + I + IN + KK
         N       = N + 1
         APUP2(N)= RO(IJ)
 3620 CONTINUE
      NTOT = (IMAX+2)*(JMAX+2)*(KMAX+2)
      CALL MPI_SSEND(APUP2(1),NTOT,MPI_REAL8,0,NTG,MPI_COMM_WORLD,
     +        IERR)
      RETURN
      END SUBROUTINE WRITMP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONORD(ICON,NBCS,NPROCE,MBPRO,NPRO,IOTY,ITURB,ISTRES,
     +     NSCAL,IMAX,JMAX,KMAX,MGM,MGRID,LEVEL,MULPHL,NPHASES,TRANSL,
     +     TWO_FLUIDL)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: NBCS,NPRO,MBPRO,IOTY,ITURB,ISTRES,NSCAL,MGM,LEVEL
      INTEGER :: IERR1,IERR2,IP,IGBL,IPRO,J,IL1,ITYPE,IGCO,NTOT,N,IBLO,
     +           MG,ILE,NN,NMPI,ITIME,NP,MMIN,IP1,JOBTO,ILENG,ICM,
     +           NPHASES
      INTEGER :: ICON(IC9,*),NPROCE(MBPRO+1,*),IMAX(MGM,*),JMAX(MGM,*),
     +     KMAX(MGM,*),MGRID(MGM,*)
      LOGICAL :: SORTED,MULPHL,TRANSL,TWO_FLUIDL

      INTEGER, ALLOCATABLE :: IPOINT(:),ICMPI(:),IJOB(:),IPROL(:)
      
      ALLOCATE(IPOINT(NPRO),ICMPI(NPRO),IJOB(NPRO),STAT=IERR1)
      ALLOCATE(IPROL(NBCS),STAT=IERR2)

      WRITE(45,*)
      IF(IERR1 == 0 .AND. IERR2 == 0) THEN
         WRITE(45,*) 'In CONORD: IPOINT etc. allocated'
      ELSE
         WRITE(45,*) 'In CONORD: Troubles in allocation'
      ENDIF

C ... this subroutine decide what order should communication be done.


C ... FIND LOCAL PROSESSOR (IPROL)
      DO IP = 1,NBCS
         IGBL  = ICON(24,IP)   ! GLOBAL BLOCK NUMBER
         IPROL(IP) = 1
         DO IPRO = 1,NPRO
            DO J = 1,NPROCE(1,IPRO)
               IF(IGBL == NPROCE(J+1,IPRO)) IPROL(IP) = IPRO !LOCAL PROCESS
            ENDDO
         ENDDO
      ENDDO
      IL1  = 2**(LEVEL-1)

C ... SET CONNECTIONS

      DO IP = 1,NBCS

      ITYPE = ICON(1,IP)
      IF(ITYPE ==  1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR. 
     &   ITYPE == 14 .OR. ITYPE == 15) THEN
      IGCO  = ICON(8,IP)   ! CONNECTIVITY BLOCK NUMBER

C ... FIND CONNECTIVITY PROCESSOR (ICON(17,IP))
      DO IPRO = 1,NPRO
         DO J = 1,NPROCE(1,IPRO)
            IF(IGCO == NPROCE(J+1,IPRO)) ICON(17,IP) = IPRO
         ENDDO
      ENDDO
      IF(IOTY == 3 .AND. IPROL(IP) /= ICON(17,IP)) ICON(18,IP) = 3
      ENDIF
      ENDDO

C ... do work variable. It has large number for prosesses that have
C ... much work and small for smaller
      DO IPRO = 1,NPRO
         NTOT = 0
         DO J = 1,NPROCE(1,IPRO)
            N    = NPROCE(J+1,IPRO)
C ... give walls work for extra calculation. Work is  equal to one slab.
            DO IP = 1,NBCS
               IBLO  = ICON(24,IP) ! GLOBAL BLOCK NUMBER
               IF(IBLO == N) THEN
               ITYPE = ICON(1,IP)
               IF(ITYPE == 8 .OR. ITYPE == 9 .OR. ITYPE == 10) THEN
                  NTOT = (ICON(5,IP)-ICON(4,IP)+1)*
     +                 (ICON(7,IP)-ICON(6,IP)+1)/IL1**2 + NTOT
               ENDIF
               ENDIF
            ENDDO
C ... end walls
            DO MG = 1,MGRID(1,N)
               ILE  = 2**(MG+LEVEL-2)
               NTOT = (IMAX(1,N)/ILE+4)*(JMAX(1,N)/ILE+2)*
     +              MAX0(KMAX(1,N)/ILE,1) + NTOT
            ENDDO
         ENDDO
         IJOB(IPRO) = NTOT
      ENDDO


C ... begin job scheduling

C ... idea of this procedure is that all the processes can and will 
C ... communicate parallel to each other and no buffering is needed. 
C ... Communication are devided for different time levels.
C ... Main goal is that every process make one patch communication in each
C ... time level. Subroutine do not take acount different size of patches
C ... or in the other words: every patches are assumed to be equal size. 
C ... Always should be possible to organized communications so 
C ... that there is only one more time level than maximum 
C ... number of communicating patches in one process.

C ... idea is continued so that less working blocks are done first and
C ... more working last. This is because less working generally are first
C ... in communications.

C ... SET ZEROS
      NN = 1
      NMPI = 0
      ITIME = 0
      DO IP = 1,NBCS
      ITYPE = ICON(1,IP)
      IF(ITYPE ==  1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &   ITYPE == 14 .OR. ITYPE == 15) ICON(21,IP) = 0
      ENDDO

      WRITE(45,*)
      WRITE(45,*)
C ... IPOINT CHECKS THAT CONNECTIONS ARE IN BALANCE
 10   DO NP = 1,NPRO
         IPOINT(NP) = 0
      ENDDO
      SORTED = .FALSE.
c 112  FORMAT('============MPI CONNECTIONS BY TIME LEVEL ',I4,
c     +     '==============='/
c     +   '  Blocks      patches       proceses   size     total sizes')
 112  FORMAT('=========== MPI CONNECTIONS BY TIME LEVEL ',I4,
     +     ' =============='/
     +   '    Blocks          patches         proceses   size ',
     +   '    total sizes')
      WRITE(45,112) ITIME+1

C ... FIND MINIMUN WORK PATCH
 31   MMIN = 1000000000
      IP1 = -1
      DO IP = 1,NBCS
      ITYPE = ICON(1,IP)
      IF(ITYPE ==  1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &   ITYPE == 14 .OR. ITYPE == 15) THEN
         NP    = ICON(17,IP)   ! CONNECTIVE PROCESS NUMBER
            IF(ICON(21,IP) == 0) THEN
            IF((IPOINT(NP) == 0 .AND. IPOINT(IPROL(IP)) == 0) .OR. 
     +           (NP == IPROL(IP))) THEN
               JOBTO = IJOB(IPROL(IP)) + IJOB(NP)
               IF(JOBTO < MMIN) IP1 = IP
               IF(JOBTO < MMIN) MMIN = JOBTO
            ENDIF
         ENDIF
         ENDIF
      ENDDO
C      DO IP = 1,NBCS
      IP = IP1
      IF(IP1 == -1) GOTO 113
      NP    = ICON(17,IP)       ! CONNECTIVE PROCESS NUMBER
      IGBL  = ICON(24,IP)       ! GLOBAL BLOCK NUMBER

C ... for mpi
      IF(NP /= IPROL(IP)) THEN
         IPOINT(NP)    = 1      ! THIS ENSURES THAT NO CONNECTIONS
         IPOINT(IPROL(IP)) = 1  ! IS MADE IN SIDE SAME STEP ???
         SORTED = .TRUE.
      ENDIF
      ICON(21,IP) = NN
      IF(ICON(21,ICON(16,IP)) /= 0) THEN
         WRITE(*,*) 'Error in patch',IP
         WRITE(*,*) IP,ICON(16,IP),ITIME,NP,IPROL(IP),
     +        ICON(21,ICON(16,IP))
         STOP 'ERROR HERE'
      ENDIF
      ICON(21,ICON(16,IP)) = NN + 1
      NN = NN + 2
      IF(ICON(18,IP) == 3) THEN
         NMPI = NMPI + 1
         WRITE(45,111) ICON(24,IP),ICON(8,IP),IP,ICON(16,IP),IPROL(IP),
     +      NP,(ICON(5,IP)-ICON(4,IP)+5)*(ICON(7,IP)-ICON(6,IP)+5)/
     +      IL1**2,IJOB(IPROL(IP)),IJOB(NP),IJOB(IPROL(IP))+IJOB(NP)
      ENDIF
      
      GOTO 31
c 111  FORMAT(I3,'<==>',I3,2X,I4,'<==>',I4,2X,I3,'<==>',I3,2X,I5,
c     +     2X,I7,'+',I7,'=',I8)
 111  FORMAT(I5,' <=>',I5,2X,I6,' <=>',I6,2X,I3,' <=>',I3,2X,I5,
     +     2X,I7,'+',I7,' =',I8)
 113  ITIME = ITIME + 1
      IF(SORTED) GOTO 10
C ... end job scheduling

      WRITE(45,*) 
      WRITE(45,*) 
      WRITE(45,*) (NN-1)/2,' pairs of connections where found.'
      WRITE(45,*) NMPI,' pairs of connections by MPI.'
      WRITE(45,*) 'Connections where made in ',ITIME-1,
     +     ' different time level.'
      WRITE(45,*) 'Best possible value is value of maximum MPI',
     +     ' connections per one process'

      CALL IWLENG(ILENG,ITURB,NSCAL,ISTRES,MULPHL,NPHASES,TRANSL,
     +   TWO_FLUIDL)

C ... COUNT MPI CONNECTIONS PER PROCESSOR (JUST FOR CHECKING)
      DO NP = 1,NPRO
         IPOINT(NP) = 0
         ICMPI(NP)  = 0
      ENDDO

      DO IP = 1,NBCS
      ITYPE = ICON(1,IP)
      IF(ITYPE ==  1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &   ITYPE == 14 .OR. ITYPE == 15) THEN
      IF(ICON(18,IP) == 3) THEN
         IPOINT(IPROL(IP)) = IPOINT(IPROL(IP)) + 1
         ICMPI(IPROL(IP))  = ICMPI(IPROL(IP)) + ILENG*2*
     +        (ICON(5,IP)-ICON(4,IP)+5)*(ICON(7,IP)-ICON(6,IP)+5)
      ENDIF
      ENDIF
      ENDDO

      WRITE(45,*) 
     +'       Process  MPI conn.  Changed variables/iteration (1.level)'
      DO NP = 1,NPRO
         WRITE(45,*) NP,IPOINT(NP),ICMPI(NP)
      ENDDO
      WRITE(45,*)

C ... Calculate approximate size of MPI connections
      ICM = 0
      DO IP = 1,NBCS
         ITYPE = ICON(1,IP)
         IF(ICON(18,IP) == 3) ICM = ICM + ILENG*2*
     +        (ICON(5,IP)-ICON(4,IP)+5)*(ICON(7,IP)-ICON(6,IP)+5)
      ENDDO

      WRITE(45,*) 
     + 'Length of the total variables (1.level) needed to exchange in'
      WRITE(45,*) 'every iteration sweep is roughly equal to',ICM

      DEALLOCATE(IPOINT,ICMPI,IJOB,IPROL)
      WRITE(45,*) 'In CONORD: IPOINT etc. deallocated.'

      RETURN
      END SUBROUTINE CONORD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NLOCAS(NLOCAL,NPNUM,NBLOCG,NBLOCK,NPROCE,MBPRO,IPRO,
     +     NPRO)

      USE MPI

      INTEGER :: ERRORCODE
      INTEGER :: NPROCE(MBPRO+1,*),NLOCAL(*),NPNUM(*)

C ... make transfer table from global to local blocks
C ... should have been done much earlier.

      DO N = 1,NBLOCG           !initialization
         NLOCAL(N) = -1
         NPNUM(N) = 123456
      ENDDO

      DO 110 NB = 1,NBLOCG
C ... FIND THE NUMBER OF PROCESS
         DO 100 IPROG = 1,NPRO
            DO 100 N = 1,NPROCE(1,IPROG)
               NBG = NPROCE(N+1,IPROG)
               IF(NBG == NB) THEN
                  IPROL = IPROG ! process number
                  NBL   = N     ! process local block number
                  NPNUM(NB) = IPROL
                  IF(IPRO == IPROG) NLOCAL(NB) = NBL
               ENDIF
 100     CONTINUE
 110  CONTINUE
      DO N = 1,NPROCE(1,IPRO)
         NBG = NPROCE(N+1,IPRO)
         N2  = NLOCAL(NBG)
         IF(N2 /= N) THEN
            WRITE(*,*) ' NLOCAL is not well defined',N2,N,NBG,IPRO
            WRITE(*,*) ' Exiting ...'
            CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
            STOP
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE NLOCAS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SUPCHI(NCHIMT,NBLOCK,NCHIM,NCHOR)

C ... put indexes so that NCHOR gives NCHIMT in order of size
C ... NCHIMT(NCHOR(1..NCHIM))

      IMPLICIT NONE

      INTEGER, DIMENSION(*) :: NCHIMT, NCHOR
      INTEGER :: N, NBLOCK, NCHIM, NCHIM2, N1, N2, NAPU
      LOGICAL :: THERE

      NCHIM2 = 0
      
      DO N = 1,NBLOCK
         IF(NCHIMT(N) >= 1) NCHIM2 = NCHIM2 + 1
      ENDDO

      NCHIM = NCHIM2
      WRITE(45,222) NCHIM
      IF(NCHIM > 1) THEN 
         WRITE(13,222) NCHIM
      ELSEIF(NCHIM == 1) THEN
         WRITE(13,223) NCHIM
      ENDIF

 222  FORMAT(/'  I found',I4,' Chimera blocks'/)
 223  FORMAT(/'  I found',I4,' Chimera block'/)


C ... Initialisation

      DO N = 1,NBLOCK
         NCHOR(N) = N
      ENDDO


 111  THERE = .FALSE.           ! Bubble sorting
      DO N = 1,NBLOCK-1
         N1 = NCHOR(N)
         N2 = NCHOR(N+1)

         IF(NCHIMT(N1) > NCHIMT(N2)) THEN
            THERE     = .TRUE.
            NAPU       = NCHOR(N)
            NCHOR(N)   = NCHOR(N+1)
            NCHOR(N+1) = NAPU
         ENDIF

      ENDDO
      IF(THERE) GOTO 111

c      DO N = 1,NCHIM         ! tiputetaan nollat pois (eli ei chimmit)
c 11      N1 = NCHOR(N)   ei tehda enaa
c         IF(NCHIMT(N1) <= 0) THEN
c            DO NN = N,NBLOCK-1
c               NCHOR(NN) = NCHOR(NN+1)
c            ENDDO
c            GOTO 11
c         ENDIF
c      ENDDO

*      IF(NCHIM > 0) THEN
*      WRITE(*,444) 
* 444  FORMAT('   Chim order          Block    level of block') 
*      DO N = 1,NBLOCK
*         IF(NCHIMT(NCHOR(N)) /= 0 .AND. NCHIMT(NCHOR(N)) /= -2)
*     +        WRITE(*,*) N,NCHOR(N),NCHIMT(NCHOR(N))
*      ENDDO
*      ENDIF

      RETURN
      END SUBROUTINE SUPCHI
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE CHSEND(FTRA,ZZZ,INUM,NTR,IPMIN,IPMAX,MAXW,NROOT,N7)

      USE MPI

      INTEGER :: ERRORCODE
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: FTRA(NTR,*),ZZZ(*)

      IPTOT = IPMAX - IPMIN + 1

      IF(IPTOT*INUM > MAXW) THEN
         WRITE(*,*) ' ZZZ size is too small. Change algorithm'
         WRITE(*,*) ' Exiting ...'
         CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
      ENDIF

c      write(*,*) 'se',inum*iptot,ntr,ipmin,ipmax,nroot
      DO N = 1,INUM
         DO I = IPMIN,IPMAX
            IZ = I - IPMIN + 1 + (N-1)*IPTOT
            ZZZ(IZ) = FTRA(I,N+3)
         ENDDO
      ENDDO

      CALL MPI_SEND(ZZZ,IPTOT*INUM,MPI_REAL8,NROOT,N7,MPI_COMM_WORLD,
     +     IERR)

      RETURN
      END SUBROUTINE CHSEND
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE CHRECV(FTRA,ZZZ,INUM,NTR,IPMIN,IPMAX,MAXW,NROO1,N7)

      USE MPI

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: FTRA(NTR,*),ZZZ(*)

      IPTOT = IPMAX - IPMIN + 1

c      write(*,*) 're',inum*iptot,ntr,ipmin,ipmax,nroo1
      CALL MPI_RECV(ZZZ,IPTOT*INUM,MPI_REAL8,NROO1,N7,MPI_COMM_WORLD,
     +     STATUS,IERR)

      DO N = 1,INUM
         DO I = IPMIN,IPMAX
            IZ = I - IPMIN + 1 + (N-1)*IPTOT
            FTRA(I,N+3) = ZZZ(IZ)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CHRECV
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SUPERC(ICON,NBCS,NCON,NCPAT)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*),NCPAT(*),IPOINT(NBCS)
      LOGICAL :: SORTED

C ... SOLVE ORDER OF SUPER CONNECTION. THIS ENSURES THAT NO
C ... BUFFERING IS NEEDED FOR MPI IMPLEMENTATION.
C ... It means that different processors communicates in same order
C ... It uses the size of the ICON(21,I) variable for each connection.
C ... ICON(21,I) is build in CONORD subroutine.

C ... Alustetaan sorttaustaulukko

      do i=1,NBCS
        ipoint(i) = i
      enddo

C ... Lasketaan nollasta poikkeavien alkioiden lukumaara
      NONOLLA = 0
      DO I=1,NBCS
        IF(ICON(21,I) /= 0) NONOLLA = NONOLLA + 1
      ENDDO

      i = 1
 641  continue
      i = i + 1
      sorted = .true.
      do j=NBCS,i,-1
         if(ICON(21,ipoint(j-1)) > ICON(21,ipoint(j))) then
            sorted      = .false.
            lapu        = ipoint(j-1)
            ipoint(j-1) = ipoint(j)
            ipoint(j)   = lapu
         endif
      enddo
      if(sorted) goto 642
      goto 641

 642    continue

C ... Nollat ovat taulukon alussa, joten skipataan ne ...
      do i=NBCS-nonolla+1,NBCS
c         write(*,*) i-(NBCS-nonolla),ipoint(i),icon(21,IPOINT(i))
         NCPAT(I-(NBCS-NONOLLA)) = IPOINT(i)
      enddo

      NCON = NONOLLA
      WRITE(45,*)
      WRITE(45,*) 'Connections',NCON
      WRITE(45,*) 'Connection order is following'

      
      WRITE(45,*) (NCPAT(I),I=1,NCON)
      WRITE(45,*)

      RETURN
      END SUBROUTINE SUPERC
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE AUXFOR(IFG,CXB,CYB,CZB,CMXB,CMYB,CMZB,DXB,DYB,DZB,
     &     QTB,QWB,QHB,TOMEGA,ALPHA,BETA,ICON,BOUNDF,NBCS,IPRO,
     &     IMETH,T,ITIMES,APATCH,XMOMG,YMOMG,ZMOMG,REFPRE,AREF,CHLREF)
C
C     This subroutine calculates the forces excerted on various
C     sub-systems based on patch-oriented data.
C
C     INPUT:
C        CXB,CYB,CZB,CMXB,CMYB,CMZB, = Forces and moments excerted
C        DXB,DYB,DZB                   on the patches
C        QTB,QWB,QHB,                = Corresponding heat fluxes
C        TOMEGA                      = Power of moving solids
C        ALPHA                       = Angle of attack
C        BETA                        = Sideslip angle
C        ICON                        = BC patch data
C        BOUNDF                      = Force group identifiers
C        NBCS                        = Number of BCs
C        IMETH                       = 1 for end of calculation
C                                    = 2 for time accurate monitoring
C        T                           = time for IMETH = 2
C        ITIMES                      = time cycle for IMETH = 2
C        APATCH                      = Patch area
C
C     OUTPUT:
C        I/O unit 4 (FORCES)
C
C     Petri Kaurinkoski, (C) 2.5.1995
C
C     Subroutine completely rewritten 4.3.1996 by Esa Salminen.
C     Now the force groups are defined in the BC-file and not in 
C     separate AUXFOR-file(s). 
C
C     MPI support added  Oct 1, 1996  ESa

      USE MPI
      USE FORCE_GROUPS
      USE NS3CO, ONLY : IC9, PARALLEL
      USE MAIN_ARRAYS, ONLY : IGRID
      USE FLIGHT, ONLY : NGRIFL, XCG, YCG, ZCG, XCGI, YCGI, ZCGI,
     &                   PSIR, THETAR, PHIR, PSIRI, THETARI, PHIRI

      IMPLICIT NONE

      REAL :: CXB(*),CYB(*),CZB(*),CMXB(*),CMYB(*),CMZB(*),APATCH(*),
     &        DXB(*),DYB(*),DZB(*),QTB(*),QWB(*),QHB(*),TOMEGA(*),
     &        CXS,CYS,CZS,DXS,DYS,DZS,CMXS,CMYS,CMZS,CLS,CDS,CSS,AAS,
     &        CAS,DAS,REFPRE,AREF,CHLREF,C,QTS,QWS,QHS,TOS,CMA,
     &        CXT,CYT,CZT,CMXT,CMYT,CMZT,ALPHA,BETA,T,
     &        CXSM,CYSM,CZSM,DXSM,DYSM,DZSM,
     &        QTSM,QWSM,QHSM,CMXSM,CMYSM,CMZSM,TOSM,AASM,CMAM,CASM,DASM,
     &        QA, QAC 

      INTEGER :: ICON(IC9,*)
      INTEGER :: I, J, IPRO, NBCS, IMETH, IERR, ITIMES
      INTEGER :: IACTU, IFLY, JFLY, IFG 

      REAL    :: XMOMG,  YMOMG, ZMOMG
      REAL    :: XMOMRJ, YMOMRJ, ZMOMRJ
      REAL    :: XMOMAJ, YMOMAJ, ZMOMAJ
      REAL    :: AX, AY, AZ, AA
      REAL    :: BX, BY, BZ, BA, BXAPU, BYAPU, BZAPU
      REAL    :: CX, CY, CZ, CA, CXAPU, CYAPU, CZAPU

      CHARACTER(LEN=80) :: BOUNDF(*)
      CHARACTER(LEN=1)  :: GNAME
      LOGICAL           :: GROUP, GROUPX, MASTER

      MASTER = IPRO == 1 
      
      C   = CHLREF
      QA  = REFPRE*AREF
      QAC = QA*C
      
      DO J=1,52    ! Scan through all alphabetics

         GROUP = .FALSE.

         CXS  = 0.0
         CYS  = 0.0
         CZS  = 0.0
         DXS  = 0.0
         DYS  = 0.0
         DZS  = 0.0
         QTS  = 0.0
         QWS  = 0.0
         QHS  = 0.0
         CMXS = 0.0
         CMYS = 0.0
         CMZS = 0.0
         CLS  = 0.0
         CDS  = 0.0
         CSS  = 0.0
         TOS  = 0.0
         AAS  = 0.0

         XMOMRJ = XMOMR(J)
         YMOMRJ = YMOMR(J)
         ZMOMRJ = ZMOMR(J)
         XMOMAJ = XMOMA(J)
         YMOMAJ = YMOMA(J)
         ZMOMAJ = ZMOMA(J)
      
         IF(J <= 26) GNAME = CHAR(J+64)   ! Uppercase characters
         IF(J >  26) GNAME = CHAR(J+70)   ! Lowercase characters

         DO I=1,NBCS

            IF(ICON(1,I) == 8  .OR.  ICON(1,I)  ==  9 
     &    .OR. ICON(1,I) == 10 .OR.  ICON(1,I)  == 15
     &    .OR. ICON(1,I) == 1  .AND. ICON(20,I) ==  7
     &    .OR. ICON(1,I) == 1  .AND. ICON(20,I) ==  8
     &    .OR. ICON(1,I) == 3  .OR.  ICON(1,I)  ==  5) THEN

               IF(BOUNDF(I)(J:J) /= ' ') THEN
                  GROUP = .TRUE.

                  CXT  =  CXB(I)
                  CYT  =  CYB(I)
                  CZT  =  CZB(I)

                  CXS  =  CXS  +  CXT
                  CYS  =  CYS  +  CYT
                  CZS  =  CZS  +  CZT

                  DXS  =  DXS  +  DXB(I)
                  DYS  =  DYS  +  DYB(I)
                  DZS  =  DZS  +  DZB(I)

                  QTS  =  QTS  +  QTB(I)
                  QWS  =  QWS  +  QWB(I)
                  QHS  =  QHS  +  QHB(I)

                  CMXT =  CMXB(I)
                  CMYT =  CMYB(I)
                  CMZT =  CMZB(I)

                  IFLY = IGRID(ICON(24,I),1)

                  IF(IFLY > 10 .AND. IFLY < 20) THEN

                     JFLY = IGRID(ICON(24,I),3)

C ... Vectors pointing from the center of gravity to two points on an
C ... arbitrary moment reference axis, e.g. flap hinge line 

                     BX = XMOMR(J) - XCGI(JFLY)
                     BY = YMOMR(J) - YCGI(JFLY)
                     BZ = ZMOMR(J) - ZCGI(JFLY)
                     CX = XMOMA(J) - XCGI(JFLY)
                     CY = YMOMA(J) - YCGI(JFLY)
                     CZ = ZMOMA(J) - ZCGI(JFLY)

                     CALL FINEUL(PSIRI,THETARI,PHIRI,BX,BY,BZ,
     &                    BXAPU,BYAPU,BZAPU,1)
                     CALL FINEUL(PSIRI,THETARI,PHIRI,CX,CY,CZ,
     &                    CXAPU,CYAPU,CZAPU,1)

C ... Rotate the arbitrary moment reference axis (two points) according 
C ... to the flying object position 

                     CALL EULFIN(PSIR,THETAR,PHIR,BXAPU,BYAPU,BZAPU,
     &                           XMOMRJ,YMOMRJ,ZMOMRJ,1)
                     XMOMRJ = XMOMRJ + XCG(JFLY) 
                     YMOMRJ = YMOMRJ + YCG(JFLY) 
                     ZMOMRJ = ZMOMRJ + ZCG(JFLY) 

                     CALL EULFIN(PSIR,THETAR,PHIR,CXAPU,CYAPU,CZAPU,
     &                           XMOMAJ,YMOMAJ,ZMOMAJ,1)
                     XMOMAJ = XMOMAJ + XCG(JFLY) 
                     YMOMAJ = YMOMAJ + YCG(JFLY) 
                     ZMOMAJ = ZMOMAJ + ZCG(JFLY) 
                     
                  ENDIF
                  
C ... Actuator disk moment reference point from local to global.

                  IF(ICON(1,I) == 1 .AND. ICON(20,I) == 8) THEN
                     IACTU = IGRID(ICON(24,I),3)
                     IF(IMETH == 1 .OR. IMETH == 2) THEN
                        CMXT = CMXT + CZT*(YCG(IACTU)-YMOMG)/C 
     &                       - CYT*(ZCG(IACTU)-ZMOMG)/C
                        CMYT = CMYT - CZT*(XCG(IACTU)-XMOMG)/C 
     &                       + CXT*(ZCG(IACTU)-ZMOMG)/C
                        CMZT = CMZT + CYT*(XCG(IACTU)-XMOMG)/C 
     &                       - CXT*(YCG(IACTU)-YMOMG)/C
                     ELSEIF (IMETH == 3) THEN
                        CMXT = CMXT + CZT*(YCG(IACTU)-YMOMG) 
     &                       - CYT*(ZCG(IACTU)-ZMOMG)
                        CMYT = CMYT - CZT*(XCG(IACTU)-XMOMG) 
     &                       + CXT*(ZCG(IACTU)-ZMOMG)
                        CMZT = CMZT + CYT*(XCG(IACTU)-XMOMG) 
     &                       - CXT*(YCG(IACTU)-YMOMG)
                     ENDIF
                  ENDIF

                  CMXS = CMXS + CMXT
                  CMYS = CMYS + CMYT
                  CMZS = CMZS + CMZT

                  TOS  = TOS  +  TOMEGA(I)
                  AAS  = AAS  +  APATCH(I)

               ENDIF

            ENDIF

         ENDDO  ! DO I=1,NBCS


         IF(PARALLEL) THEN

            CALL MPI_REDUCE(GROUP,GROUPX,1,MPI_LOGICAL,MPI_LOR,0,
     &           MPI_COMM_WORLD,IERR)
            IF(MASTER) GROUP = GROUPX
            CALL MPI_BCAST(GROUP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

            IF(GROUP) THEN

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

               IF(MASTER) THEN
                  CXS  = CXSM
                  CYS  = CYSM
                  CZS  = CZSM
                  DXS  = DXSM
                  DYS  = DYSM
                  DZS  = DZSM
                  QTS  = QTSM
                  QWS  = QWSM
                  QHS  = QHSM
                  CMXS = CMXSM
                  CMYS = CMYSM
                  CMZS = CMZSM
                  TOS  = TOSM
                  AAS  = AASM
               ENDIF

            ENDIF

         ENDIF  ! PARALLEL


         IF(MASTER) THEN

         CLS = (-CXS*SIN(ALPHA)           + CYS*COS(ALPHA))
         CDS =   CXS*COS(ALPHA)*COS(BETA) + CYS*SIN(ALPHA)*COS(BETA)
     &         + CZS*SIN(BETA)
         CSS =   CXS*COS(ALPHA)*SIN(BETA) + CYS*SIN(ALPHA)*SIN(BETA)
     &         + CZS*COS(BETA)

C ... Moment/force refrence axis

          AX = XMOMAJ - XMOMRJ
          AY = YMOMAJ - YMOMRJ
          AZ = ZMOMAJ - ZMOMRJ

          AA = SQRT(AX**2 + AY**2 + AZ**2)

          IF(AA > 0.0) THEN
             AX = AX/AA
             AY = AY/AA
             AZ = AZ/AA
          ELSE
             AX = 1.0
             AY = 0.0
             AZ = 0.0
          ENDIF

          
         IF(GROUP .AND. IMETH == 1) THEN

C ... Move global moment reference point MOMG to local MOMR (force group)
C ... Note that loads are non-dimensional when IMETH==1.

            CMXS = CMXS + CZS*(YMOMG-YMOMRJ)/C - CYS*(ZMOMG-ZMOMRJ)/C
            CMYS = CMYS - CZS*(XMOMG-XMOMRJ)/C + CXS*(ZMOMG-ZMOMRJ)/C
            CMZS = CMZS + CYS*(XMOMG-XMOMRJ)/C - CXS*(YMOMG-YMOMRJ)/C

            CMA  = AX*CMXS + AY*CMYS + AZ*CMZS  
            CAS  = AX*CXS  + AY*CYS  + AZ*CZS  
            DAS  = AX*DXS  + AY*DYS  + AZ*DZS  

            IF(PARALLEL) THEN
               WRITE(4,*)
               WRITE(4,*)
               WRITE(4,*)'     *************************************'
               WRITE(4,*)'     TOTAL GROUP FORCES FROM ALL PROCESSES:'
            ENDIF

            WRITE(4,300) GNAME,FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24)

c         IF(XMOMAJ < -4.0E20) THEN
c         IF(XMOMAJ < -4.0E10) THEN  ! Matti Palin 25.09.2017
         IF(XMOMAJ < -4.0D10) THEN   ! Esa Oct 25, 2018

            WRITE(4,310) CLS,CDS,CSS,
     &                   CXS,CYS,CZS,
     &                   DXS,DYS,DZS,
     &                   CXS-DXS,CYS-DYS,CZS-DZS,
     &                   CMXS,CMYS,CMZS,
     &                   XMOMRJ,YMOMRJ,ZMOMRJ,
     &                   CXS*QA,CYS*QA,CZS*QA,
     &                   DXS*QA,DYS*QA,DZS*QA,
     &                  (CXS-DXS)*QA,(CYS-DYS)*QA,(CZS-DZS)*QA,
     &                   CMXS*QAC,CMYS*QAC,CMZS*QAC,
     &                   QTS,QWS,QHS,
     &                   QTS/AAS,QWS/AAS,QHS/AAS,
     &                   TOS,AAS

         ELSE 

            WRITE(4,311) CLS,CDS,CSS,
     &                   CXS,CYS,CZS,
     &                   DXS,DYS,DZS,
     &                   CXS-DXS,CYS-DYS,CZS-DZS,
     &                   CAS,DAS,CAS-DAS,
     &                   CMXS,CMYS,CMZS,
     &                   XMOMRJ,YMOMRJ,ZMOMRJ,
     &                   CMA,
     &                   XMOMAJ,YMOMAJ,ZMOMAJ,
     &                   CXS*QA,CYS*QA,CZS*QA,
     &                   DXS*QA,DYS*QA,DZS*QA,
     &                  (CXS-DXS)*QA,(CYS-DYS)*QA,(CZS-DZS)*QA,
     &                   CAS*QA,DAS*QA,(CAS-DAS)*QA,
     &                   CMXS*QAC,CMYS*QAC,CMZS*QAC,
     &                   QTS,QWS,QHS,
     &                   QTS/AAS,QWS/AAS,QHS/AAS,
     &                   TOS,AAS

         ENDIF

         IF(PARALLEL) THEN
            WRITE(4,*)'     *************************************'
            WRITE(4,*)
            WRITE(4,*)
         ENDIF

         ELSEIF(GROUP .AND. IMETH == 2) THEN

            WRITE(102+J,111) T,ITIMES,CLS,CDS,CSS,CXS,CYS,CZS,
     &           DXS,DYS,DZS,CMXS,CMYS,CMZS,QTS,QWS,QHS,
     &           QTS/AAS,QWS/AAS,QHS/AAS,TOS,AAS

         ELSEIF(GROUP .AND. IMETH == 3) THEN
            
C ... Move global moment reference point MOMG to local MOMR (force group)
C ... Note that loads are dimensional when IMETH==3.

            CMXS = CMXS + CZS*(YMOMG-YMOMRJ) - CYS*(ZMOMG-ZMOMRJ)
            CMYS = CMYS - CZS*(XMOMG-XMOMRJ) + CXS*(ZMOMG-ZMOMRJ)
            CMZS = CMZS + CYS*(XMOMG-XMOMRJ) - CXS*(YMOMG-YMOMRJ)

            WRITE(IFG) ITIMES,J,
     &                 FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24),REAL(T,4), 
     &                 REAL(CXS,4),REAL(CYS,4),REAL(CZS,4),
     &                 REAL(DXS,4),REAL(DYS,4),REAL(DZS,4),
     &                 REAL(CMXS,4),REAL(CMYS,4),REAL(CMZS,4),
     &                 REAL(REFPRE,4),REAL(AREF,4),REAL(CHLREF,4),
     &                 REAL(XMOMRJ,4),REAL(YMOMRJ,4),REAL(ZMOMRJ,4),
     &                 REAL(XMOMAJ,4),REAL(YMOMAJ,4),REAL(ZMOMAJ,4)

         ENDIF

         ENDIF  ! MASTER

      ENDDO  ! DO J = 1,52


 111  FORMAT(G11.4,I6,20(1X,G11.4))

 300  FORMAT(/6X,'Force group: ',1A1,2X,A24)

 310  FORMAT(/6X,'CL  = ',E12.5,5X,'CD  = ',E12.5,5X,'CS  = ',E12.5/
     &        6X,'CX  = ',E12.5,5X,'CY  = ',E12.5,5X,'CZ  = ',E12.5/
     &        6X,'DX  = ',E12.5,5X,'DY  = ',E12.5,5X,'DZ  = ',E12.5/
     &        6X,'VX  = ',E12.5,5X,'VY  = ',E12.5,5X,'VZ  = ',E12.5/
     &        6X,'CMX = ',E12.5,5X,'CMY = ',E12.5,5X,'CMZ = ',E12.5/
     &        6X,'XMR = ',E12.5,' m   ','YMR = ',E12.5,' m   ',
     &           'ZMR = ',E12.5,' m   '/
     &        6X,'FX  = ',E11.4,' N    ','FY  = ',E11.4,' N    ',
     &           'FZ  = ',E11.4,' N    '/
     &        6X,'FXp = ',E11.4,' N    ','FYp = ',E11.4,' N    ',
     &           'FZp = ',E11.4,' N    '/
     &        6X,'FXf = ',E11.4,' N    ','FYf = ',E11.4,' N    ',
     &           'FZf = ',E11.4,' N    '/
     &        6X,'FMX = ',E11.4,' Nm   ','FMY = ',E11.4,' Nm   ',
     &           'FMZ = ',E11.4,' Nm   '/
     &        6X,'QT  = ',E11.4,' W    ','QW  = ',E11.4,' W    ',
     &           'QH  = ',E11.4,' W    '/
     &        6X,'Q"T = ',E11.4,' W/m2 ','Q"W = ',E11.4,' W/m2 ',
     &           'Q"H = ',E11.4,' W/m2 '/
     &        6X,'TOM = ',E11.4,' W    ','AREA= ',E11.4,' m2')

 311  FORMAT(/6X,'CL  = ',E12.5,5X,'CD  = ',E12.5,5X,'CS  = ',E12.5/
     &        6X,'CX  = ',E12.5,5X,'CY  = ',E12.5,5X,'CZ  = ',E12.5/
     &        6X,'DX  = ',E12.5,5X,'DY  = ',E12.5,5X,'DZ  = ',E12.5/
     &        6X,'VX  = ',E12.5,5X,'VY  = ',E12.5,5X,'VZ  = ',E12.5/
     &        6X,'CA  = ',E12.5,5X,'DA  = ',E12.5,5X,'VA  = ',E12.5/
     &        6X,'CMX = ',E12.5,5X,'CMY = ',E12.5,5X,'CMZ = ',E12.5/
     &        6X,'XMR = ',E12.5,' m   ','YMR = ',E12.5,' m   ',
     &           'ZMR = ',E12.5,' m   '/
     &        6X,'CMA = ',E12.5/
     &        6X,'XMA = ',E12.5,' m   ','YMA = ',E12.5,' m   ',
     &           'ZMA = ',E12.5,' m   '/
     &        6X,'FX  = ',E11.4,' N    ','FY  = ',E11.4,' N    ',
     &           'FZ  = ',E11.4,' N    '/
     &        6X,'FXp = ',E11.4,' N    ','FYp = ',E11.4,' N    ',
     &           'FZp = ',E11.4,' N    '/
     &        6X,'FXf = ',E11.4,' N    ','FYf = ',E11.4,' N    ',
     &           'FZf = ',E11.4,' N    '/
     &        6X,'FA  = ',E11.4,' N    ','FAp = ',E11.4,' N    ',
     &           'FAf = ',E11.4,' N    '/
     &        6X,'FMX = ',E11.4,' Nm   ','FMY = ',E11.4,' Nm   ',
     &           'FMZ = ',E11.4,' Nm   '/
     &        6X,'QT  = ',E11.4,' W    ','QW  = ',E11.4,' W    ',
     &           'QH  = ',E11.4,' W    '/
     &        6X,'Q"T = ',E11.4,' W/m2 ','Q"W = ',E11.4,' W/m2 ',
     &           'Q"H = ',E11.4,' W/m2 '/
     &        6X,'TOM = ',E11.4,' W    ','AREA= ',E11.4,' m2')

      RETURN
      END SUBROUTINE AUXFOR
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TOTBAL(PAVE,TAVE,ROAVE,EAVE,EINAVE,UUAVE,RKAVE,
     +     VOLTOT,ITURB,IPRO,QMFIN,QMFOUT,QMEIN,QMEOUT,TOMTOT,
     +     TOMTOG,QGFIN,QGFOUT,QGEIN,QGEOUT,MULPHL,ROGAVE,EVAPGE)
C
C     This subroutine calculates the integrated balances over
C     prosessors
C
C     OUTPUT:
C        I/O unit 4 (FORCES001)
C

      USE MPI

      IMPLICIT NONE

      REAL :: PAVE, ROAVE, EAVE, EINAVE, UUAVE, RKAVE, QMFIN,
     &        TOMTOT, TOMTOG, QGFIN, QGFOUT, QGEIN, 
     &        QGEOUT, EFFIC, QMEIG, QMEOUG, QGFOUG, QGEOUG, QGEIG,
     &        QGFIG, QMFOUG, QMFIG, RKGLO, VGLO, ROGLO, EGLO, EINGLO,
     &        UUGLO, PGLO, TAVE, TGLO, VOLTOT, QMFOUT, QMEIN, QMEOUT,
     &        ROGAVE,EVAPGE, ROGGLO,EVAGLO    
      INTEGER :: ITURB, IPRO, IERR

      LOGICAL :: MULPHL, MASTER

      MASTER = IPRO == 1

      CALL MPI_REDUCE(PAVE*VOLTOT,PGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(TAVE*VOLTOT,TGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(VOLTOT,VGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(ROAVE,ROGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(EAVE,EGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(EINAVE,EINGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(UUAVE,UUGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL MPI_REDUCE(RKAVE,RKGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      ENDIF
      IF(MULPHL) THEN
      CALL MPI_REDUCE(ROGAVE,ROGGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(EVAPGE,EVAGLO,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      ENDIF


      IF(MASTER) THEN
         
      WRITE(4,*)
      WRITE(4,*)'------ INTEGRATED BALANCES OVER PROCESSORS -----------'
      WRITE(4,*) 'AVERAGE PRESSURE    = ',REAL(PGLO/VGLO,4),  ' N/m2'
      WRITE(4,*) 'AVERAGE TEMPERATURE = ',REAL(TGLO/VGLO,4),  ' K'
      WRITE(4,*) 'TOTAL MASS          = ',REAL(ROGLO,4), ' kg'
      IF(MULPHL) THEN
      WRITE(4,*) 'TOTAL GAS MASS      = ',REAL(ROGGLO,4),' kg'
      WRITE(4,*) 'TOTAL EVAPORATION   = ',REAL(EVAGLO,4),' kg/s'
      ENDIF
      WRITE(4,*) 'TOTAL ENERGY        = ',REAL(EGLO,4),  ' J'
      WRITE(4,*) 'TOTAL INT. ENERGY   = ',REAL(EINGLO,4),' J'
      WRITE(4,*) 'TOTAL KINETIC ENERGY= ',REAL(UUGLO,4), ' J'
      IF(ITURB >= 3 .AND. ITURB /= 8)
     +WRITE(4,*) 'TOTAL TURB. ENERGY  = ',REAL(RKGLO,4), ' J'

      ENDIF

      RETURN
      END SUBROUTINE TOTBAL
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE TOTDIS(VISAVE,PROAVE,EPSAVE,EPWAVE,ITURB,IPRO,
     +     TOMTOG)
C
C     This subroutine calculates the integrated dissipations over
C     prosessors
C
C     OUTPUT:
C        I/O unit 4 (FORCES001)
C

      USE MPI

      CALL MPI_REDUCE(VISAVE,VISAVG,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(PROAVE,PROAVG,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(EPSAVE,EPSAVG,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(EPWAVE,EPWAVG,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)

      IF(IPRO == 1) THEN
      WRITE(4,*)
      WRITE(4,*)'------------- TOTAL DISSIPATIONS OVER PROCESSORS -----'
      WRITE(4,*)
      WRITE(4,*) 'Viscous dissipation = ',VISAVG, ' W'
      WRITE(4,*) 'Production power    = ',PROAVG,' W'
      WRITE(4,*) 'Dissipation tilde   = ',EPSAVG,' W'
      WRITE(4,*) 'Dissipation wall    = ',EPWAVG,' W'
      WRITE(4,*) 'Turb. dissipation   = ',(EPSAVG+EPWAVG),' W'
      WRITE(4,*) 'Total dissipation   = ',(EPSAVG+EPWAVG+VISAVG),' W'
      WRITE(4,*) 'Ratios: P/epsilon   = ',PROAVG/(EPSAVG+EPWAVG+1.E-20)
      WRITE(4,*) 'viscous/turbu diss. = ',VISAVG/(EPSAVG+EPWAVG+1.E-20)
      WRITE(4,*)
      IF(ABS(TOMTOG) >= 1.E-5) THEN
         WRITE(4,*) 'TOMEGA/T. dissi.    = ',-TOMTOG/
     +        (EPSAVG+EPWAVG+VISAVG+1.E-20)
      ENDIF
      ENDIF

      RETURN
      END SUBROUTINE TOTDIS
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE COLPLO(APU2,MAXW,SPLIT,NOB,IDOB,JDOB,KDOB,NSSB,ISSB,NB,
     +                  NPROCE,MBPRO,NPRO)

      USE NS3CO, ONLY : PARALLEL
      
      REAL :: APU2(MAXW)
      DIMENSION :: IDOB(NB),JDOB(NB),KDOB(NB),NSSB(NB),ISSB(7,NB)
      INTEGER :: NPROCE(MBPRO+1,*)
      LOGICAL :: SPLIT
      CHARACTER(LEN=25) :: NAME, TMPSTR
      CHARACTER(LEN=7)  :: TESTER(24)

C ... piti tulla ohjelma joka keraa PLOTPIt mutta loppu puhti

      DATA TESTER/
     &     ' TAUBI ',' TAUBJ ',' TAUBK ',' TAUTI ',' TAUTJ ',' TAUTK ',
     &     ' TAU2BI',' TAU2BJ',' TAU2BK',' TAU2TI',' TAU2TJ',' TAU2TK',
     &     ' QAUBI ',' QAUBJ ',' QAUBK ',' QAUTI ',' QAUTJ ',' QAUTK ',
     &     ' QAU2BI',' QAU2BJ',' QAU2BK',' QAU2TI',' QAU2TJ',' QAU2TK'/

      IF(.NOT.SPLIT .AND. .NOT. PARALLEL) RETURN

C ... Check that MAXW is big enough
      NMAX = 0
      DO NBL = 1,NOB
         NPX  = JDOB(NBL)*KDOB(NBL)
         NPY  = IDOB(NBL)*KDOB(NBL)
         NPZ  = IDOB(NBL)*JDOB(NBL)
         NSI  = 1
         NSJ  = NSI + NPX*8
         NSK  = NSJ + NPY*8
         NEND = NSK + NPZ*8
         NMAX = MAX0(NMAX,NEND)
         IF(NEND > MAXW) THEN
            WRITE(*,*) 'Cannot collect PLOTPIs. Too small MAXW.'
            WRITE(*,*) 'Should be at least ',NMAX
            RETURN
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE COLPLO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C



