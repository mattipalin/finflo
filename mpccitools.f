C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MPCCI_INITCOUPLING

C ... Calling routine for MpCCI initialization.

      USE MPCCIVARS

      USE INTEGERS,    ONLY : IPRO

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: IERR
 
      LOGICAL :: MASTER

      
      MASTER = IPRO == 1

      IF(MPCCIL) THEN

         CALL SIZE_FOR_MPCCI

         ALLOCATE(fnodecoor(fnnodes,3))
         ALLOCATE(fnodecoor_orig(fnnodes,3))
         ALLOCATE(fnodedisp(3*fnnodes))
         ALLOCATE(felemnodes(fnelems,4))
         ALLOCATE(felemforce(fnelems,3))
         ALLOCATE(felempres(fnelems))
         ALLOCATE(felemtemp(fnelems))
         ALLOCATE(felemarea(fnelems))
         ALLOCATE(fnodechimt(fnnodes))

C ... Original wetted surface for the grid point distance calculation
C ... in FINFLO coordinate system.
         
         CALL MESH_FOR_MPCCI(3,IN,JN,KN)  ! Original wetted surface 

         fnodecoor_orig = fnodecoor       ! Save the original wetted surface

         IF(RMESHT == 2) CALL GP_DISTANCES(30,.FALSE.)

C ... Original wetted surface in the FEM model coordinate system.
         
         CALL MESH_FOR_MPCCI(1,IN,JN,KN)  ! Original wetted surface 

         fnodecoor_orig = fnodecoor       ! Save the original wetted surface
         
         CALL MESH_FOR_MPCCI(2,IN,JN,KN)  ! Deformed data for the coupled code

         IF(MASTER) call initcoupling(t,dt,strongcoupling)

         IF(PARALLEL) THEN
            CALL MPI_BCAST(strongcoupling,1,MPI_INTEGER,0,
     &           MPI_COMM_WORLD,IERR)
         ENDIF

      ENDIF

      RETURN     
      END SUBROUTINE MPCCI_INITCOUPLING
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MPCCI_DOTRANSFER(FIRSTC)

C ... Calling routine for MpCCI data transfer actions.

      USE MPI

      USE MPCCIVARS

      USE INTEGERS,    ONLY : IPRO

      USE NS3CO

      IMPLICIT NONE

      INTEGER :: I, IERR, IS
      LOGICAL :: MASTER, FIRSTC

      REAL :: apu1, apu2, apu3

      MASTER = IPRO == 1

      IF(MPCCIL) THEN

         CALL MESH_FOR_MPCCI(2,IN,JN,KN)

         CALL DATA_FOR_MPCCI(IN,JN,KN)
         
C .. Only the master process exchanges data with the coupling code
       
         IF(MASTER) call dotransfer(t,dt,transfered,strongcoupling)

C ... Send all displacements to all nodes in a parallel simulation

         IF(PARALLEL) THEN
            CALL MPI_BCAST(transfered,1,MPI_INTEGER,0,
     &           MPI_COMM_WORLD,IERR)
         ENDIF

         IF(transfered == 0) RETURN

         IF(PARALLEL) THEN
            CALL MPI_BCAST(fnodedisp,3*fnnodes,MPI_REAL8,0,
     &           MPI_COMM_WORLD,IERR)            
         ENDIF

C ... Rotate the false FEM coordinate system back to the FINFLO system

         SELECT CASE(FEMXYZ)

         CASE(-1)  ! CFD coordinates xyz correspond to xz-y in FEM model
         
            do i=1,fnnodes
            
               apu1 = fnodecoor(i,1)
               apu2 = fnodecoor(i,2)
               apu3 = fnodecoor(i,3)

               fnodecoor(i,1) =  apu1 
               fnodecoor(i,2) =  apu3        
               fnodecoor(i,3) = -apu2

            enddo

         CASE(0)  ! Same coordinate system in CFD and FEM
         
            do i=1,fnnodes
            
               apu1 = fnodecoor(i,1)
               apu2 = fnodecoor(i,2)
               apu3 = fnodecoor(i,3)

               fnodecoor(i,1) = apu1 
               fnodecoor(i,2) = apu2        
               fnodecoor(i,3) = apu3

            enddo

         CASE(1)  ! CFD coordinates xyz correspond to yzx in FEM model (HN)
         
            do i=1,fnnodes
            
               apu1 = fnodecoor(i,1)
               apu2 = fnodecoor(i,2)
               apu3 = fnodecoor(i,3)

               fnodecoor(i,1) = apu2 
               fnodecoor(i,2) = apu3        
               fnodecoor(i,3) = apu1

            enddo

         CASE(2)  ! CFD coordinates xyz correspond to zxy in FEM model
         
            do i=1,fnnodes
            
               apu1 = fnodecoor(i,1)
               apu2 = fnodecoor(i,2)
               apu3 = fnodecoor(i,3)

               fnodecoor(i,1) = apu3 
               fnodecoor(i,2) = apu1        
               fnodecoor(i,3) = apu2

            enddo

         END SELECT
         
C ... Deform the grid

         IF(.NOT.FIRSTC) CALL REMESH

      ENDIF

      RETURN     
      END SUBROUTINE MPCCI_DOTRANSFER
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MPCCI_EXITCOUPLING

C ... Calling routine for MpCCI exit.

      USE MPCCIVARS

      USE INTEGERS,    ONLY : IPRO

      USE NS3CO

      IMPLICIT NONE

      LOGICAL :: MASTER


      MASTER = IPRO == 1

      IF(MPCCIL) THEN

         IF(MASTER) call exitcoupling()

         DEALLOCATE(fnodecoor)
         DEALLOCATE(fnodecoor_orig)
         DEALLOCATE(fnodedisp)
         DEALLOCATE(felemforce)
         DEALLOCATE(felemnodes)
         DEALLOCATE(felempres)
         DEALLOCATE(felemtemp)
         DEALLOCATE(felemarea)
         DEALLOCATE(fnodechimt)

      ENDIF

      RETURN     
      END SUBROUTINE MPCCI_EXITCOUPLING
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE adapterOutput(msg)

C ... Print output to the screen

        character msg*(*)

        WRITE(*,'(A)') msg
c       important: Flush output buffer
        call FLUSH(6)
      END SUBROUTINE adapterOutput
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE error(msg)

      character msg*(*)

      WRITE(*,'(2A)') 'ERROR: ',msg

      call FLUSH(6)

      END SUBROUTINE error
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE getSurfaceID(surfaceID,name)

C ...
C ...  Get ID of coupling component
C ...

      IMPLICIT NONE
      character :: name*(*)
      integer   :: i, length, surfaceID

*        include 'data.h'
c        WRITE(*,*) 'getSurfaceID, name is "',name,'"'
c        call Flush(6)
c        surfaceID = MAXFOUNDATIONS + 1
*        WRITE(*,*) 'getSurfaceID, between'
*        call Flush(6)
*        length = LEN(name)
*        WRITE(*,*) 'getSurfaceID, length is "',length,'"'
*        call Flush(6)
*        DO i=1, fndnumber
*          IF(fnamelen(i) == length) THEN
*            IF(fname(i) == name) THEN
*              surfaceID = i
*            END IF
*          END IF
*        END DO
*        IF(surfaceID > MAXFOUNDATIONS) THEN
*          call error('Coupling surface not found!!!')
*        END IF

      surfaceID = 1
      WRITE(*,*) 'The ID is "',surfaceID,'"'
      call Flush(6)

      RETURN
      END SUBROUTINE getSurfaceID
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE getInfo(fnd,nnodes,nelems,dim)

      USE MPCCIVARS

      INTEGER :: fnd, nnodes, nelems, dim

      WRITE(*,*) 'getInfo, ID is "',fnd,'"'
      call Flush(6)
      WRITE(*,*) 'getting nnodes'
      call Flush(6)
      nnodes = fnnodes
      WRITE(*,*) 'nnodes',nnodes
      call Flush(6)
      nelems = fnelems
      WRITE(*,*) 'nelems',nelems
      call Flush(6)
      dim = 3
      WRITE(*,*) 'finished getInfo'
      call Flush(6)

      RETURN
      END SUBROUTINE getInfo
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE getMesh(fnd,nodeNumbers,nodeCoords,elemNodes)

C ... Get surface mesh data for MpCCI from FINFLO's MPCCIVARS module.

      USE MPCCIVARS

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: nodeCoords
      INTEGER, DIMENSION(*) :: nodeNumbers, elemNodes
      INTEGER :: fnd, n, i
      write(*,*) 
      WRITE(*,*) 'getMesh ID is "',fnd,'"'
      call Flush(6)

      IF(MPCCIL) THEN

         n=0
         DO i=1,fnnodes
            nodeNumbers(n+1)  = i
            nodeCoords(3*n+1) = fnodecoor(i,1)
            nodeCoords(3*n+2) = fnodecoor(i,2)
            nodeCoords(3*n+3) = fnodecoor(i,3)
            n = n+1
         END DO
      
         n=0
         
         DO i=1,fnelems
            elemNodes(4*n+1) = felemnodes(i,1)
            elemNodes(4*n+2) = felemnodes(i,2)
            elemNodes(4*n+3) = felemnodes(i,3)
            elemNodes(4*n+4) = felemnodes(i,4)
            n = n+1
         END DO

      ENDIF

      WRITE(*,*) 'getMesh finished'
      
      call Flush(6)
      
      RETURN
      END SUBROUTINE getMesh
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE getPressure(fnd,values)

C ... Get pressure difference for MpCCI from FINFLO's MPCCIVARS module.

      USE MPCCIVARS

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values
      INTEGER :: fnd, n, i

      n=0

      DO i=1,fnelems
         n = n+1
         values(n) = felempres(i)
      END DO
      
      RETURN
      END SUBROUTINE getPressure
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE getWallTemp(fnd,values)

C ... Get wall temperature for MpCCI from FINFLO's MPCCIVARS module.

      USE MPCCIVARS

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values
      INTEGER :: fnd, n, i

      n = 0

      DO i=1,fnelems
         n = n+1
         values(n) = felemtemp(i)
      END DO
      
      RETURN
      END SUBROUTINE getWallTemp
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE getWallForce(fnd,values)

C ... Get element forces for MpCCI from FINFLO's MPCCIVARS module.

      USE MPCCIVARS

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values
      INTEGER :: fnd, n, i

      n=0

      DO i=1,fnelems
         values(3*n+1) = felemforce(i,1)
         values(3*n+2) = felemforce(i,2)
         values(3*n+3) = felemforce(i,3)
         n = n+1
      END DO

      RETURN
      END SUBROUTINE getWallForce
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE getphysicaltime(values)

C ... Get FINFLO's physical time.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      values(1) = T
      
      RETURN
      END SUBROUTINE getphysicaltime
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE gettimestepsize(values)

C ... Get FINFLO's time step size.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      values(1) = DT
      
      RETURN
      END SUBROUTINE gettimestepsize
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE gettimestepcount(values)

C ... Get FINFLO's time step count.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      values(1) = DBLE(ICYTOT)
      
      RETURN
      END SUBROUTINE gettimestepcount
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE getreferencepressure(values)

C ... Get FINFLO's reference pressure.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      values(1) = REFPRE
      
      RETURN
      END SUBROUTINE getreferencepressure
C     
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE putNPosition(fnd,values)

C ... Put displacements from coupled code to FINFLO's MPCCIVARS module.

      USE MPCCIVARS

      USE NS3CO

      IMPLICIT NONE

      integer :: fnd, n, d, i
      
      REAL, DIMENSION(*) :: values
      
      n = 0

      SELECT CASE(FEMXYZ)  ! If displacements relative to the previous shape
                           ! are needed, replace fnodecoor_orig by fnodecoor.
      
      CASE(-1)  ! CFD coordinates xyz correspond to xz-y in FEM model

      DO i=1,fnnodes

         fnodedisp(3*n+1) =   values(3*n+1) - fnodecoor_orig(i,1)
         fnodedisp(3*n+2) =   values(3*n+3) - fnodecoor_orig(i,3)
         fnodedisp(3*n+3) = -(values(3*n+2) - fnodecoor_orig(i,2))

         n = n + 1

      ENDDO
      
      CASE(0)  ! Same coordinate system in CFD and FEM

      DO i=1,fnnodes

         fnodedisp(3*n+1) = values(3*n+1) - fnodecoor_orig(i,1)
         fnodedisp(3*n+2) = values(3*n+2) - fnodecoor_orig(i,2)
         fnodedisp(3*n+3) = values(3*n+3) - fnodecoor_orig(i,3)

         n = n + 1

      ENDDO

      CASE(1)  ! CFD coordinates xyz correspond to yzx in FEM model (HN)

      DO i=1,fnnodes

         fnodedisp(3*n+1) = values(3*n+2) - fnodecoor_orig(i,2)
         fnodedisp(3*n+2) = values(3*n+3) - fnodecoor_orig(i,3)
         fnodedisp(3*n+3) = values(3*n+1) - fnodecoor_orig(i,1)

         n = n + 1

      ENDDO

      CASE(2)  ! CFD coordinates xyz correspond to zxy in FEM model

      DO i=1,fnnodes

         fnodedisp(3*n+1) = values(3*n+3) - fnodecoor_orig(i,3)
         fnodedisp(3*n+2) = values(3*n+1) - fnodecoor_orig(i,1)
         fnodedisp(3*n+3) = values(3*n+2) - fnodecoor_orig(i,2)

         n = n + 1

      ENDDO

      END SELECT

      RETURN
      END SUBROUTINE putNPosition
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE putphysicaltime(values)

C ... Physical time from coupling code.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      T = values(1)
      
      RETURN
      END SUBROUTINE putphysicaltime
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE puttimestepsize(values)

C ... Time step size from coupling code.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      DT = values(1)
      
      RETURN
      END SUBROUTINE puttimestepsize
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE puttimestepcount(values)

C ... Timemstep count from coupling code.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      ICYTOT = NINT(values(1))
      
      RETURN
      END SUBROUTINE puttimestepcount
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 
      SUBROUTINE putreferencepressure(values)

C ... Reference pressure from coupling code.

      USE NS3CO

      IMPLICIT NONE

      REAL, DIMENSION(*) :: values

      REFPRE = values(1)
      
      RETURN
      END SUBROUTINE putreferencepressure
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C 


