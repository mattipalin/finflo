C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MODALFSI_INITCOUPLING

C ... Calling routine for Modal FSI initialization.

      USE MPCCIVARS

      USE NS3CO

      USE INTEGERS, ONLY : IPRO

      IMPLICIT NONE

      LOGICAL :: MASTER

      MASTER = IPRO == 1
      
      IF (MODALFSIL) THEN

************************************************************************
*      if(master) then
*      OPEN(999,FILE='modal-forces.jans.dat_orig',
*     &        STATUS='OLD',FORM='FORMATTED')
*      endif
************************************************************************

         eta = 0.0
         etadot = 0.0
         etaddot = 0.0
         
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

         CALL READ_JOINTS
         CALL READ_MODE_SHAPES

C ... Original wetted surface for the grid point distance calculation
C ... in FINFLO coordinate system.
         
         CALL MESH_FOR_MPCCI(3,IN,JN,KN)  ! Original wetted surface 

         fnodecoor_orig = fnodecoor       ! Save the original wetted surface

         IF(RMESHT == 2) CALL GP_DISTANCES(30,.FALSE.)

      ENDIF

************************************************************************
C ... Map modes to the CFD grid 
      if(master) CALL MAP_MODES
************************************************************************

      RETURN     
      END SUBROUTINE MODALFSI_INITCOUPLING
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MODALFSI_DOTRANSFER

C ... Calling routine for Modal FSI actions.

      USE MPCCIVARS
      USE NS3CO

      IMPLICIT NONE

      INTEGER :: I
      
      IF (MODALFSIL) THEN
         
         CALL MESH_FOR_MPCCI(2,IN,JN,KN)  ! Deformed wetted surface

         CALL DATA_FOR_MPCCI(IN,JN,KN)

         CALL MODALTRANSFER

         CALL REMESH
         
      ENDIF

      RETURN     
      END SUBROUTINE MODALFSI_DOTRANSFER
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE MODALFSI_EXITCOUPLING

C ... Calling routine for Modal FSI exit.

      USE MPCCIVARS

      USE NS3CO

      USE INTEGERS, ONLY : IPRO

      IMPLICIT NONE

      LOGICAL :: MASTER

      MASTER = IPRO == 1

      IF (MODALFSIL) THEN

         DEALLOCATE(fnodecoor)
         DEALLOCATE(fnodecoor_orig)
         DEALLOCATE(fnodedisp)
         DEALLOCATE(felemforce)
         DEALLOCATE(felemnodes)
         DEALLOCATE(felempres)
         DEALLOCATE(felemtemp)
         DEALLOCATE(felemarea)
         DEALLOCATE(fnodechimt)

         DEALLOCATE(jointxyz)
         DEALLOCATE(jointxyzo)
         DEALLOCATE(jointdisp)

         IF (MASTER) THEN
************************************************************************
            DEALLOCATE(shapefelem)
************************************************************************
            DEALLOCATE(shapexyz)
            DEALLOCATE(shapedydx)
            DEALLOCATE(shapedxdz)
            DEALLOCATE(jointxyzn)
            DEALLOCATE(jointf)
         ENDIF

      ENDIF

      RETURN     
      END SUBROUTINE MODALFSI_EXITCOUPLING
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE READ_JOINTS

C ... This subroutine must be tailored to fit the modal joints file
C ... or the modal joints file must be tailored to fit the read
C ... statements in this subroutine.      

      USE INTEGERS, ONLY : IPRO
      USE NS3CO, ONLY : PARALLEL
      USE MPCCIVARS
      USE MPI

      IMPLICIT NONE
      
      LOGICAL :: MASTER
      CHARACTER(LEN=80) :: LINE
      INTEGER :: I, IJOINT, IERR
      REAL, DIMENSION(3) :: CGJOINTS

      MASTER = IPRO == 1

      IF (MASTER) THEN

         OPEN(789,FILE=JOINTS_FILE,STATUS='OLD',FORM='FORMATTED')

C ... Skip comment lines starting with '#' in the beginning of the file.
         
 100     CONTINUE
         READ(789,'(A80)') LINE
         IF (LINE(1:1) == '#') GOTO 100
         BACKSPACE 789

C ... Read number of joints.
         
         READ(789,*) IDIMJOINTS, JDIMJOINTS

         NJOINTS = IDIMJOINTS*JDIMJOINTS

         ALLOCATE(jointxyzn(3,NJOINTS),STAT=IERR) 
         ALLOCATE(jointf(3,NJOINTS),STAT=IERR) 

      ENDIF
 
      IF (PARALLEL) CALL MPI_BCAST(NJOINTS,1,MPI_INTEGER,0,
     &                   MPI_COMM_WORLD,IERR)
      IF (PARALLEL) CALL MPI_BCAST(IDIMJOINTS,1,MPI_INTEGER,0,
     &                   MPI_COMM_WORLD,IERR)
      IF (PARALLEL) CALL MPI_BCAST(JDIMJOINTS,1,MPI_INTEGER,0,
     &                   MPI_COMM_WORLD,IERR)

      ALLOCATE(jointxyz(3,NJOINTS),STAT=IERR)
      ALLOCATE(jointxyzo(3,NJOINTS),STAT=IERR)
      ALLOCATE(jointdisp(3,NJOINTS),STAT=IERR)

C ... Read joint locations.

      IF (MASTER) THEN

         DO IJOINT = 1,NJOINTS
            READ(789,*) I, JOINTXYZ(1,I), JOINTXYZ(2,I), JOINTXYZ(3,I)
        ENDDO

************************************************************************
* Save original modal model for plotting.        
*               FNAME='SI_MODE0.P3D'
*               FNAME=TRIM(FNAME)
*               jntx = jointxyz(1,:)
*               jnty = jointxyz(2,:)
*               jntz = jointxyz(3,:)
*            CALL WritePlot3dGeometry(FNAME,888,
*     &           0,1,1,idimjoints-1,jdimjoints-1,1-1,
*     &           jntx,jnty,jntz,jntz)
************************************************************************

        jointxyzo = jointxyz

      ENDIF
     
      IF (PARALLEL) THEN
         CALL MPI_BCAST(JOINTXYZ,3*NJOINTS,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JOINTXYZO,3*NJOINTS,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
      ENDIF
         
      IF (MASTER) CLOSE(789)
    
      RETURN
      END SUBROUTINE READ_JOINTS
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE READ_MODE_SHAPES

C ... This subroutine must be tailored to fit the mode shapes file
C ... or the mode shapes file must be tailored to fit the read
C ... statements in this subroutine.      

      USE INTEGERS, ONLY : IPRO
      USE NS3CO, ONLY : PARALLEL
      USE MPCCIVARS
      USE MPI

      IMPLICIT NONE
      
      LOGICAL :: MASTER
      CHARACTER(LEN=80) :: LINE
      INTEGER :: I, IMODE, IJOINT, IERR
      real :: mpf, memf, mew

*      real :: ak(10)
*      ak = 0.0
      
      MASTER = IPRO == 1

      IF (MASTER) THEN

         OPEN(790,FILE=MODE_SHAPES_FILE,STATUS='OLD',FORM='FORMATTED')

C ... Read number of frequencies, i.e. modes.
         
         READ(790,*) NFREQ

      ENDIF

      IF (PARALLEL) CALL MPI_BCAST(NFREQ,1,MPI_INTEGER,0,
     &                   MPI_COMM_WORLD,IERR)

      IF (MASTER) THEN
     
************************************************************************
         ALLOCATE(SHAPEFELEM(3,NFREQ*FNELEMS),STAT=IERR) 
************************************************************************
         ALLOCATE(SHAPEXYZ(3,NFREQ*NJOINTS),STAT=IERR) 
         ALLOCATE(SHAPEDYDX(NFREQ*NJOINTS),STAT=IERR) 
         ALLOCATE(SHAPEDXDZ(NFREQ*NJOINTS),STAT=IERR) 
         
C ... Read the natural frequencies.

         DO IMODE = 1,NFREQ
            READ(790,*) FREQ(IMODE)
         ENDDO

         IMODE = 0
         
 100     CONTINUE

         IMODE = IMODE + 1
         
C ... Skip comment lines starting with '#' above each mode.
         
 200     CONTINUE

         READ(790,'(A80)') LINE
         IF (LINE(1:1) == '#') THEN
            GOTO 200
         ENDIF
         BACKSPACE 790

************************************************************************

* Modal Effective Mass Fraction (AGARD 445.6)
      if(imode == 1) memf = 5.082378D-01
      if(imode == 2) memf = 1.050280D-01
      if(imode == 3) memf = 9.780459D-02
      if(imode == 4) memf = 5.224878D-02
      if(imode == 5) memf = 2.917382D-02
      if(imode == 6) memf = 2.989627D-02
      if(imode == 7) memf = 4.834736D-05
      if(imode == 8) memf = 3.433454D-03

* Modal Participation Factor (AGARD 445.6)
      if(imode == 1) mpf = -9.089388D-01
      if(imode == 2) mpf =  4.131937D-01
      if(imode == 3) mpf = -3.987317D-01  
      if(imode == 4) mpf =  2.914333D-01
      if(imode == 5) mpf = -2.177699D-01
      if(imode == 6) mpf =  2.204498D-01
      if(imode == 7) mpf = -8.865180D-03
      if(imode == 8) mpf =  7.470794D-02

* Modal Effective Mass = Modal Effective Weight (AGARD 445.6)     
      if(imode == 1) mew =  8.261697D-01
      if(imode == 2) mew =  1.707290D-01
      if(imode == 3) mew =  1.589870D-01
      if(imode == 4) mew =  8.493338D-02
      if(imode == 5) mew =  4.742371D-02
      if(imode == 6) mew =  4.859810D-02
      if(imode == 7) mew =  7.859141D-05
      if(imode == 8) mew =  5.581276D-03
       
************************************************************************

C ... Read one mode shape.

         DO IJOINT = 1,NJOINTS

            READ(790,*) I, SHAPEXYZ(1,(IMODE-1)*NJOINTS+I),
     &                     SHAPEXYZ(2,(IMODE-1)*NJOINTS+I),
     &                     SHAPEXYZ(3,(IMODE-1)*NJOINTS+I),
     &                     SHAPEDYDX((IMODE-1)*NJOINTS+I),
     &                     SHAPEDXDZ((IMODE-1)*NJOINTS+I)

************************************************************************
            SHAPEXYZ(2,(IMODE-1)*NJOINTS+I) =
     &      SHAPEXYZ(2,(IMODE-1)*NJOINTS+I)
************************************************************************

         ENDDO
         
************************************************************************
* Save mode shapes for plotting.
*               CALL NUMCH1(FILEE,IMODE)
*               FNAME='SI_MODE'//FILEE//'.P3D'
*               FNAME=TRIM(FNAME)
*               IE = (IMODE-1)*NJOINTS + 1
*               IV = IMODE*NJOINTS
*               jntx = jointxyz(1,:) + SHAPEXYZ(1,IE:IV)
*               jnty = jointxyz(2,:) + SHAPEXYZ(2,IE:IV)
*               jntz = jointxyz(3,:) + SHAPEXYZ(3,IE:IV)
*            CALL WritePlot3dGeometry(FNAME,888,
*     &           0,1,1,idimjoints-1,jdimjoints-1,1-1,
*     &           jntx,jnty,jntz,jntz)
************************************************************************

         IF (IMODE < NFREQ) GOTO 100

         CLOSE(790)
         
      ENDIF

      RETURN
      END SUBROUTINE READ_MODE_SHAPES
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MODALTRANSFER

C ... Transfer data between FINFLO and the modal structural model. 

      USE NS3CO
      USE MPI
      USE MPCCIVARS
      USE INTEGERS, ONLY : IPRO
      USE CONSTANTS, ONLY : PII

      IMPLICIT NONE

      INTEGER :: IMODE, INODE, IS, MODES, IERR
      LOGICAL :: MASTER
      real    :: tim

      MASTER = IPRO == 1

      MODES = ACTIVEFREQS
      
      IF (MASTER) THEN

         FOMEGA = 2.0*PII*FREQ
        
         CALL MAP_CFD_LOADS

**************************************************************************
*         read(999,*) tim,fmx(1,1),fmx(2,1),fmx(3,1),fmx(4,1)
*************************************************************************
         
         DO IMODE = 1,MODES

            CALL MODAL_LOADS(IMODE)

*            CALL SOLVE(IMODE)
            CALL IMPLICITSOLVE(IMODE)

         ENDDO

******** Modal forces **************************************************
         write(681,*) REAL(T,4),
     &                REAL(FM(1,1),4),
     &                REAL(FM(2,1),4),
     &                REAL(FM(3,1),4),
     &                REAL(FM(4,1),4)
************************************************************************
         write(*,'(F6.3,4F12.3)') REAL(T,4),
     &                REAL(FM(1,1),4),
     &                REAL(FM(2,1),4),
     &                REAL(FM(3,1),4),
     &                REAL(FM(4,1),4)
         write(*,'(F6.3,4F12.3)') REAL(T,4),
     &                REAL(FMX(1,1),4),
     &                REAL(FMX(2,1),4),
     &                REAL(FMX(3,1),4),
     &                REAL(FMX(4,1),4)
******** Etat **********************************************************
         write(679,*) REAL(T,4),
     &                REAL(ETA(1,1),4),
     &                REAL(ETA(2,1),4),
     &                REAL(ETA(3,1),4),
     &                REAL(ETA(4,1),4)
************************************************************************

         CALL MODAL_POSITIONS

      ENDIF

      IF (PARALLEL) THEN
         CALL MPI_BCAST(JOINTXYZ,3*NJOINTS,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JOINTDISP,3*NJOINTS,MPI_REAL8,0,
     &        MPI_COMM_WORLD,IERR)
      ENDIF

      RETURN
      END SUBROUTINE MODALTRANSFER
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MAP_CFD_LOADS

C ... This subroutine projects the modal structural model joint points
C ... onto the CFD model surface in order to extrapolate the flow
C ... solution loads.
      
      USE MPCCIVARS

      IMPLICIT NONE
      
      INTEGER :: I, J, N, IJOINT, IMELEM, NMELEMS
      INTEGER :: J1, J2, J3, J4, ME1, ME2, ME3, ME4

      REAL, DIMENSION(3)   :: V1, V2, V3, CFDTF, MELTF, JOINTSTF
      REAL, DIMENSION(3,5) :: A1, B1, E1
      REAL :: AELEM, APROJ, ASHARE
      REAL, DIMENSION(3,(IDIMJOINTS-1)*(JDIMJOINTS-1)) :: MELEMF
      
C ... Scan through all modal model elements and calculate the CFD model
C ... element projections and shares on the modal elements. Element     
C ... voxelization would accelerate the algorithm.
      
      IMELEM = 0
      MELEMF = 0.0
      
      DO J = 1,JDIMJOINTS-1
         DO I = 1,IDIMJOINTS-1

            IMELEM = IMELEM + 1

            J1 = (J-1)*IDIMJOINTS+I
            J2 = J1 + 1
            J3 = J*IDIMJOINTS+I+1
            J4 = J3 - 1
            
C ... Modal surface element corner point coordinates.
               
            A1(:,1) = jointxyz(:,J1)
            A1(:,2) = jointxyz(:,J2)
            A1(:,3) = jointxyz(:,J3)
            A1(:,4) = jointxyz(:,J4)

C ... Scan trhough all CFD surface elements against one modal element.
            
            DO N = 1,fnelems

C ... CFD model surface element corner point coordinates.
               
               E1(:,1) = fnodecoor(felemnodes(N,1),:)
               E1(:,2) = fnodecoor(felemnodes(N,2),:)
               E1(:,3) = fnodecoor(felemnodes(N,3),:)
               E1(:,4) = fnodecoor(felemnodes(N,4),:)

C ... CFD model surface element projection on modal element.

               CALL CAREA(A1,E1,AELEM,APROJ,ASHARE)

               IF(ASHARE > 0.0) THEN
                   MELEMF(:,IMELEM) = MELEMF(:,IMELEM) +
     &                                ASHARE*felemforce(N,:)
               ENDIF

            ENDDO
            
         ENDDO
      ENDDO

      NMELEMS = IMELEM

C ... Element forces to node forces.

      JOINTF = 0.0
      
      DO J = 1,JDIMJOINTS
         DO I = 1,IDIMJOINTS

            J1  = (J-1)*IDIMJOINTS + I
            
            ME1 = (J-2)*(IDIMJOINTS-1) + I - 1
            ME2 = ME1 + 1
            ME3 = ME2 + (IDIMJOINTS-1)
            ME4 = ME3 - 1

            if(I == 1 .and. j == 1) then  ! corner points
               JOINTF(:,J1) = MELEMF(:,ME3)/4.0  
            elseif(I == 1 .and. j == jdimjoints) then
               JOINTF(:,J1) = MELEMF(:,ME2)/4.0  
            elseif(I == idimjoints .and. j == 1) then
               JOINTF(:,J1) = MELEMF(:,ME4)/4.0  
            elseif(I == idimjoints .and. j == jdimjoints) then 
               JOINTF(:,J1) = MELEMF(:,ME1)/4.0  
            elseif(I == 1 .and. j /= 1 .and. j /= jdimjoints) then  ! edges
               JOINTF(:,J1) = (MELEMF(:,ME2) + MELEMF(:,ME3))/4.0  
            elseif(I == idimjoints .and. j/=1 .and. j/=jdimjoints) then 
               JOINTF(:,J1) = (MELEMF(:,ME1) + MELEMF(:,ME4))/4.0  
            elseif(j == 1 .and. i /= 1 .and. i /= idimjoints) then 
               JOINTF(:,J1) = (MELEMF(:,ME3) + MELEMF(:,ME4))/4.0  
            elseif(j == jdimjoints .and. i/=1 .and. i/=idimjoints) then 
               JOINTF(:,J1) = (MELEMF(:,ME1) + MELEMF(:,ME2))/4.0  
            else  ! inner points
               JOINTF(:,J1) = (MELEMF(:,ME1) + MELEMF(:,ME2) +
     &                         MELEMF(:,ME3) + MELEMF(:,ME4))/4.0
            endif

         ENDDO
      ENDDO
      
C ... Calculate the CFD solution load total and the joints load total. 

      CFDTF = 0.0

      DO N = 1,fnelems
         CFDTF = CFDTF + felemforce(N,:)
      ENDDO

      MELTF = 0.0

      DO IMELEM = 1,NMELEMS
         MELTF = MELTF + MELEMF(:,IMELEM)
      ENDDO

      JOINTSTF = 0.0

      DO J1 = 1,NJOINTS
         JOINTSTF(:) = JOINTSTF(:) + JOINTF(:,J1)
      ENDDO

************************************************************************
      write(*,*) REAL(CFDTF(1),4), REAL(MELTF(1),4), REAL(JOINTSTF(1),4)
      write(*,*) REAL(CFDTF(2),4), REAL(MELTF(2),4), REAL(JOINTSTF(2),4)
      write(*,*) REAL(CFDTF(3),4), REAL(MELTF(3),4), REAL(JOINTSTF(3),4)
************************************************************************
      
C ... Scale the joint point forces to give the CFD total load.

*      DO I = 1,NJOINTS
*         jointf(1,I) = jointf(1,I)*CFDTF(1)/JOINTSTF(1)
*         jointf(2,I) = jointf(2,I)*CFDTF(2)/JOINTSTF(2)
*         jointf(3,I) = jointf(3,I)*CFDTF(3)/JOINTSTF(3)
*      ENDDO

      RETURN
      END SUBROUTINE MAP_CFD_LOADS
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MAP_DISPLACEMENTS

C ... This subroutine calculates the modal structural model displacements
C ... at the locations of the CFD model surface (wetted) points.
      
      USE MPCCIVARS
      USE NS3CO
      USE MPI
      USE INTEGERS, ONLY : IPRO
      
      IMPLICIT NONE
      
      INTEGER :: I, J, IP, J1, J2, J3, J4, IT, ICODE, IERR

      REAL, DIMENSION(3) :: S1, S2, S3, S4, S5, V1, V2, V3
      REAL, DIMENSION(3) :: C1, C2, C3, C4, C5
      REAL, DIMENSION(3) :: T1, T2, T3, TC, TN, TR, D, Q, DQ, POINT
      REAL, DIMENSION(3) :: T4, T5, T6, D1, D2, D3

      DO IP = 1,fnnodes

         point = fnodecoor_orig(IP,:)
            
         DO J = 1,JDIMJOINTS-1
            DO I = 1,IDIMJOINTS-1

               J1 = (J-1)*IDIMJOINTS+I
               J2 = J1 + 1
               J3 = J*IDIMJOINTS+I+1
               J4 = J3 - 1
               
C ... Original modal surface element corner point coordinates. 

               S1 = jointxyzo(:,J1)
               S2 = jointxyzo(:,J2)
               S3 = jointxyzo(:,J3)
               S4 = jointxyzo(:,J4)
               S5 = (S1+S2+S3+S4)/4.
           
C ... Original modal surface element normat

               V1 = S3 - S1
               V2 = S2 - S4

               CALL CROSS_PRODUCT(V1,V2,V3)
            
C ... Deformed modal surface element corner point coordinates.

               C1 = jointxyz(:,J1)
               C2 = jointxyz(:,J2)
               C3 = jointxyz(:,J3) 
               C4 = jointxyz(:,J4)
               C5 = (C1+C2+C3+C4)/4.
               
C ... Divide elements into four triangles.
               
               DO IT = 1,4

                  SELECT CASE(IT)

                  CASE(1)
                     T1 = S1
                     T2 = S2
                     T3 = S5
                     D1 = C1
                     D2 = C2
                     D3 = C5
                  CASE(2)
                     T1 = S2
                     T2 = S3
                     T3 = S5
                     D1 = C2
                     D2 = C3
                     D3 = C5
                  CASE(3)
                     T1 = S3
                     T2 = S4
                     T3 = S5
                     D1 = C3
                     D2 = C4
                     D3 = C5
                  CASE(4)
                     T1 = S4
                     T2 = S1
                     T3 = S5
                     D1 = C4
                     D2 = C1
                     D3 = C5

                  END SELECT

               TC = (T1 + T2 + T3) / 3.0  ! Original triangle center.

C ... Enlarge the triangle a little bit in order to ensure successful mapping

               T4 = TC + 1.15*(T1-TC)
               T5 = TC + 1.15*(T2-TC)
               T6 = TC + 1.15*(T3-TC)
                  
               TR = MAX(SUM((T4-TC)**2),SUM((T5-TC)**2),SUM((T6-TC)**2))
               TR = SQRT(TR)
               
C ... Project a wetted surface point onto the modal model surface by
C ... calculating an intersection point of a normal and a plane defined
C ... by a triangle.

                  CALL POINT_PLANE(POINT,T4,T5,T6,V3,TC,TR,D,Q,ICODE)
 
                  IF(ICODE == 1) THEN

C ... Transform point Q -> DQ onto the deformed triangle

                     CALL TRANSFORM_TRIANGLE(Q,T1,T2,T3,DQ,D1,D2,D3)
                     
                     fnodedisp(3*(IP-1)+1) = DQ(1) - Q(1) 
                     fnodedisp(3*(IP-1)+2) = DQ(2) - Q(2) 
                     fnodedisp(3*(IP-1)+3) = DQ(3) - Q(3)

                     GOTO 100

                  ENDIF
                  
               ENDDO

            ENDDO
         ENDDO

 100     CONTINUE

      ENDDO

      RETURN
      END SUBROUTINE MAP_DISPLACEMENTS
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MAP_MODES

C ... This subroutine interpolates the mode shapes to the CFD model
C ... surface element center points.
      
      USE MPCCIVARS
      USE NS3CO
      USE MPI
      USE INTEGERS, ONLY : IPRO
      
      IMPLICIT NONE
      
      INTEGER :: I, J, N, J1, J2, J3, J4, IT, ICODE, IMODE, IERR

      REAL, DIMENSION(3) :: S1, S2, S3, S4, S5, V1, V2, V3
      REAL, DIMENSION(3) :: T1, T2, T3, TC, TN, TR, D, Q, DQ, POINT
      REAL, DIMENSION(3) :: M1, M2, M3, MM1, MM2, MM3, MM4, MM5
      REAL, DIMENSION(3,4) :: E1
      REAL :: AT, A1, A2, A3


      DO IMODE = 1,NFREQ
      
      DO N = 1,fnelems

C ... CFD surface element corner points 
         
         E1(:,1) = fnodecoor_orig(felemnodes(N,1),:)
         E1(:,2) = fnodecoor_orig(felemnodes(N,2),:)
         E1(:,3) = fnodecoor_orig(felemnodes(N,3),:)
         E1(:,4) = fnodecoor_orig(felemnodes(N,4),:)

C ... CFD surface element center point

         point = (E1(:,1)+E1(:,2)+E1(:,3)+E1(:,4))/4.
            
         DO J = 1,JDIMJOINTS-1
            DO I = 1,IDIMJOINTS-1

               J1 = (J-1)*IDIMJOINTS+I
               J2 = J1 + 1
               J3 = J*IDIMJOINTS+I+1
               J4 = J3 - 1
               
C ... Original modal surface element corner point coordinates and
C ... mode shapes               

               S1 = jointxyzo(:,J1)
               S2 = jointxyzo(:,J2)
               S3 = jointxyzo(:,J3)
               S4 = jointxyzo(:,J4)
               S5 = (S1+S2+S3+S4)/4.

               MM1 = shapexyz(:,(IMODE-1)*njoints + J1)
               MM2 = shapexyz(:,(IMODE-1)*njoints + J2)
               MM3 = shapexyz(:,(IMODE-1)*njoints + J3)
               MM4 = shapexyz(:,(IMODE-1)*njoints + J4)
               MM5 = (MM1 + MM2 + MM3 + MM4)/4. 
           
C ... Original surface element normal

               V1 = S3 - S1
               V2 = S2 - S4

               CALL CROSS_PRODUCT(V1,V2,V3)
               
C ... Divide elements into four triangles.
               
               DO IT = 1,4

                  SELECT CASE(IT)

                  CASE(1)
                     T1 = S1
                     T2 = S2
                     T3 = S5
                     M1 = MM1
                     M2 = MM2
                     M3 = MM5
                  CASE(2)
                     T1 = S2
                     T2 = S3
                     T3 = S5
                     M1 = MM2
                     M2 = MM3
                     M3 = MM5
                  CASE(3)
                     T1 = S3
                     T2 = S4
                     T3 = S5
                     M1 = MM3
                     M2 = MM4
                     M3 = MM5
                  CASE(4)
                     T1 = S4
                     T2 = S1
                     T3 = S5
                     M1 = MM4
                     M2 = MM1
                     M3 = MM5

                  END SELECT

               TC = (T1 + T2 + T3) / 3.0  ! Original triangle center.
               
C ... Project a wetted surface point onto the modal model surface by
C ... calculating an intersection point of a normal and a plane defined
C ... by a triangle.
                  
               TR = MAX(SUM((T1-TC)**2),SUM((T2-TC)**2),SUM((T3-TC)**2))
               TR = SQRT(TR)

                  CALL POINT_PLANE(POINT,T1,T2,T3,V3,TC,TR,D,Q,ICODE)
 
                  IF(ICODE == 1) THEN

                     CALL CROSS_PRODUCT(T2-T1,T3-T1,V3)
                     AT = 0.5*SQRT(DOT_PRODUCT(V3,V3))
                     CALL CROSS_PRODUCT(Q-T2,Q-T3,V3)
                     A1 = 0.5*SQRT(DOT_PRODUCT(V3,V3))
                     CALL CROSS_PRODUCT(Q-T1,Q-T3,V3)
                     A2 = 0.5*SQRT(DOT_PRODUCT(V3,V3))
                     CALL CROSS_PRODUCT(Q-T1,Q-T2,V3)
                     A3 = 0.5*SQRT(DOT_PRODUCT(V3,V3))

                     shapefelem(:,(imode-1)*fnelems+N) =
     &               A1/AT*M1 + A2/AT*M2 + A3/AT*M3
                     
                     GOTO 100

                  ENDIF
                  
               ENDDO

            ENDDO
         ENDDO

 100     CONTINUE

      ENDDO

      ENDDO
      
      RETURN
      END SUBROUTINE MAP_MODES
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE TRANSFORM_TRIANGLE(P,A,B,C,PX,AX,BX,CX)

C ... This subroutine transforms triangle point P location to another 
C ... deformed triangle (PX).
      
      IMPLICIT NONE

      REAL, DIMENSION(3) :: P, A, B, C, PX, AX, BX, CX
      REAL, DIMENSION(3) :: AP, AB, AC, APxAB, APxAC, ACxAB, ABxAC
      REAL               :: U, V
      
      AP = P - A
      AB = B - A
      AC = C - A

      CALL CROSS_PRODUCT(AP,AB,APxAB)
      CALL CROSS_PRODUCT(AP,AC,APxAC)
      CALL CROSS_PRODUCT(AC,AB,ACxAB)
      CALL CROSS_PRODUCT(AB,AC,ABxAC)

      V = DOT_PRODUCT(APxAB,ACxAB)/DOT_PRODUCT(ACxAB,ACxAB)
      U = DOT_PRODUCT(APxAC,ABxAC)/DOT_PRODUCT(ABxAC,ABxAC)
       
      PX = AX + U*(BX-AX) + V*(CX-AX)
      
      RETURN
      END SUBROUTINE TRANSFORM_TRIANGLE
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MODAL_LOADS(IMODE)

      USE MPCCIVARS

      IMPLICIT NONE

      INTEGER :: IMODE, INODE, IS

      REAL, DIMENSION(3) :: VA, VB
      REAL :: memf, mpf, mew
      
      FM(IMODE,1) = 0.0
      FMX(IMODE,1) = 0.0

* Modal Effective Mass Fraction (AGARD 445.6)
      if(imode == 1) memf = 5.082378D-01
      if(imode == 2) memf = 1.050280D-01
      if(imode == 3) memf = 9.780459D-02
      if(imode == 4) memf = 5.224878D-02
      if(imode == 5) memf = 2.917382D-02
      if(imode == 6) memf = 2.989627D-02
      if(imode == 7) memf = 4.834736D-05
      if(imode == 8) memf = 3.433454D-03

* Modal Participation Factor (AGARD 445.6)
      if(imode == 1) mpf = -9.089388D-01
      if(imode == 2) mpf =  4.131937D-01
      if(imode == 3) mpf = -3.987317D-01  
      if(imode == 4) mpf =  2.914333D-01
      if(imode == 5) mpf = -2.177699D-01
      if(imode == 6) mpf =  2.204498D-01
      if(imode == 7) mpf = -8.865180D-03
      if(imode == 8) mpf =  7.470794D-02

* Modal Effective Mass = Modal Effective Weight (AGARD 445.6)      
      if(imode == 1) mew = 8.261697D-01
      if(imode == 2) mew = 1.707290D-01
      if(imode == 3) mew = 1.589870D-01
      if(imode == 4) mew = 8.493338D-02
      if(imode == 5) mew = 4.742371D-02
      if(imode == 6) mew = 4.859810D-02
      if(imode == 7) mew = 7.859141D-05
      if(imode == 8) mew = 5.581276D-03

      
      DO INODE = 1,NJOINTS
         
         IS = (IMODE-1)*NJOINTS + INODE

         VA = SHAPEXYZ(:,IS)
         VB = JOINTF(:,INODE)
                
         FM(IMODE,1) = FM(IMODE,1) + DOT_PRODUCT(VA,VB)

      ENDDO

************************************************************************
C ... Modal loads based on mode shapes at the CFD element center points.
C ... Not used but calculated for comparison with the modal loads based 
C ... on mode shapes at the modal model joint points.

      do inode=1,fnelems

         IS = (IMODE-1)*fnelems + INODE
         VA = SHAPEFELEM(:,IS)
         VB = felemforce(inode,:)

         FMX(IMODE,1) = FMX(IMODE,1) + DOT_PRODUCT(VA,VB)

      enddo
************************************************************************

      RETURN
      END SUBROUTINE MODAL_LOADS
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE SOLVE(IMODE)

C ... Integrate modes shape factors.

      USE NS3CO, ONLY : DT, T
      USE MPCCIVARS

      IMPLICIT NONE

      INTEGER :: IMODE
      REAL    :: A0, AP1, AM1

C ... Initial disturbance
*      if(t <= 0.0003) eta(imode,2) = fomega(imode)*0.00005/2.0
*      if(t <= 0.0025) eta(1,2) = fomega(imode)*0.00001
*      if(imode==1 .and. t <= 0.005) eta(1,2) = fomega(imode)*0.00005
      
      A0  = -2./DT**2 + FOMEGA(IMODE)**2/2.
      AP1 =  1./DT**2 + ZETA(IMODE)*FOMEGA(IMODE)/DT+FOMEGA(IMODE)**2/4.
      AM1 =  1./DT**2 - ZETA(IMODE)*FOMEGA(IMODE)/DT+FOMEGA(IMODE)**2/4. 

      ETA(IMODE,1) = 1./AP1*((FM(IMODE,1) + 2.*FM(IMODE,2) +
     &                        FM(IMODE,3))/4.0 -
     &                    A0*ETA(IMODE,2) - AM1*ETA(IMODE,3))

      ETA(IMODE,3) = ETA(IMODE,2) 
      ETA(IMODE,2) = ETA(IMODE,1)
      FM(IMODE,3)  = FM(IMODE,2)   
      FM(IMODE,2)  = FM(IMODE,1)

      RETURN
      END SUBROUTINE SOLVE
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      
      SUBROUTINE IMPLICITSOLVE(IMODE)

      USE NS3CO, ONLY : DT, T
      USE MPCCIVARS

      IMPLICIT NONE
      
      REAL    :: M, C, K, GAMMA, BETA, EPSILON, ETA_PREV, ETADDOT_PREV 
      INTEGER :: IMODE, ITER, MAX_ITER
      LOGICAL :: CONVERGED
      
      GAMMA = 0.5               ! Newmark gamma parameter
      BETA  = 0.25              ! Newmark beta parameter
      EPSILON = 1.0E-8          ! convergence tolerance
      MAX_ITER = 1000           ! maximum number of iterations
      
      M = 1.0
      C = 2.*ZETA(IMODE)*FOMEGA(IMODE)  
      K = FOMEGA(IMODE)**2
      
! Initial disturbance

*      if(t <= 0.0011) then
*         ETADOT(1,2) =  3.5
*         ETADOT(2,2) =  3.5
*         ETADOT(3,2) = 10.0
*         ETADOT(4,2) = 10.0
*      endif

!     Check for convergence

      CONVERGED = .FALSE.
      ITER = 1
      
      DO WHILE (.NOT.CONVERGED .AND. ITER < MAX_ITER)

!     Store previous displacement (Firstly from the previous time step level)
                  
      ETADDOT_PREV = ETADDOT(IMODE,1)

!     Calculate new acceleration

*        ETADDOT(IMODE,1) = (FM(IMODE,1) - C*ETADOT(IMODE,1)
*     &                                  - K*ETA(IMODE,1)) / M
        ETADDOT(IMODE,1) = (FMX(IMODE,1) - C*ETADOT(IMODE,1)
     &                                   - K*ETA(IMODE,1)) / M

!     Update displacement and velocity

*        ETA(IMODE,1) = ETA(IMODE,2) + DT*ETADOT(IMODE,1)
*     &                              + (0.5-BETA)*DT**2*ETADDOT(IMODE,1)
*     &                              + BETA*DT**2*ETADDOT(IMODE,2)   

        ETA(IMODE,1) = ETA(IMODE,2) + DT*ETADOT(IMODE,2)
     &                              + (0.5-BETA)*DT**2*ETADDOT(IMODE,2)
     &                              + BETA*DT**2*ETADDOT(IMODE,1)   

        ETADOT(IMODE,1) = ETADOT(IMODE,2) 
     &                  + DT*((1.0-GAMMA)*ETADDOT(IMODE,2)
     &                  + GAMMA*ETADDOT(IMODE,1))

!     Calculate error

        IF (ABS(ETADDOT(IMODE,1)-ETADDOT_PREV) < EPSILON) THEN
           CONVERGED = .TRUE.
        END IF

        ITER = ITER + 1

      END DO
         
      ETA(IMODE,2)     = ETA(IMODE,1)
      ETADOT(IMODE,2)  = ETADOT(IMODE,1)
      ETADDOT(IMODE,2) = ETADDOT(IMODE,1)
      FM(IMODE,2)      = FM(IMODE,1)

      RETURN
         
      END SUBROUTINE IMPLICITSOLVE
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE MODAL_POSITIONS

C ... Calculate the joint point displacements and deform the modal model. 

      USE MPCCIVARS
      
      IMPLICIT NONE

      INTEGER :: I, IM, IS
      real    :: tim

************************************************************************
*      READ(999,*) tim,eta(1,1),eta(2,1),eta(3,1),eta(4,1)
*      write(*,*) real(tim,4),real(eta(1,1),4),
*     &           real(eta(2,1),4),real(eta(3,1),4),real(eta(4,1),4)
************************************************************************
      
      jointdisp = 0.0

      DO IM = 1,ACTIVEFREQS

         DO I = 1,njoints

            IS = (IM-1)*njoints + I
            
            jointdisp(1,I) = jointdisp(1,I) + eta(IM,1)*shapexyz(1,IS)
            jointdisp(2,I) = jointdisp(2,I) + eta(IM,1)*shapexyz(2,IS)
            jointdisp(3,I) = jointdisp(3,I) + eta(IM,1)*shapexyz(3,IS)

         ENDDO
*            write(*,*) IM,maxval(abs(jointdisp(2,:)))
      ENDDO

      jointxyz = jointxyzo + jointdisp
      
      RETURN
      END SUBROUTINE MODAL_POSITIONS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SIZE_FOR_MPCCI

C ... Size of data (number of elements and number of nodes) to be 
C ... transferred to the coupled code in multiphysics simulations, e.g.
C ... MpCCI or Modal FSI.
      
      USE MPI
      USE CHARACTERS
      USE FORCE_GROUPS
      USE MPCCIVARS

      USE INTEGERS,    ONLY : NBCS, IPRO

      USE MAIN_ARRAYS, ONLY : ICON, IG, IDIMS, JDIMS, NSPG, NSPP         

      USE NS3CO

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SOLIDL

      INTEGER :: NSP, I, J, K, L, N, M, ISP, IG1, IF1, KX1, KX2, KY1,
     &    KY2, KZ1, KZ2, KX1W, KX2W, KY1W, KY2W, KZ1W, KZ2W, LA,
     &    IPA, IPT, ISTR, JSTR, KSTR,
     &    IWALL, IPSTR, JPSTR, IERR

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMSW, JDIMSW, ITYPEW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMT,  JDIMT,  ITYPET
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPEG, ITYPEGW
      INTEGER :: IPG, NBCSG, NSPTOTG, NSPI

      CHARACTER(LEN=1) :: FGNAME
      

      IF(PARALLEL) THEN 

C ... Total number of boundary condition patches from all processes (NBCSG)
         CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)
         
C ... Total number of selected surface patches from all processes (NSPTOTG)
         CALL MPI_ALLREDUCE(NSPTOT,NSPTOTG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)

      ELSE
         NBCSG   = NBCS
         NSPTOTG = NSPTOT
      ENDIF

      NSPI = 0

      ALLOCATE(IDIMSW(NSPTOTG),JDIMSW(NSPTOTG),ITYPEW(NSPTOTG))
      ALLOCATE(IDIMT(NSPTOTG), JDIMT(NSPTOTG), ITYPET(NSPTOTG))

      IDIMSW = 0; JDIMSW = 0; ITYPEW = 0
      IDIMT  = 0; JDIMT  = 0; ITYPET = 0

      ALLOCATE(ITYPEG(NBCSG),ITYPEGW(NBCSG))

      ITYPEG   = 0; ITYPEGW   = 0

      DO IPG = 1,NBCSG
         DO NSP = 1,NSPTOT
            DO K=1,80
               IF(BOUNDF(NSPP(NSP))(K:K) /= ' ') THEN
                  FGNAME = BOUNDF(NSPP(NSP))(K:K)
                  IF(ICHAR(FGNAME) >= 65 .AND. ICHAR(FGNAME) <= 90)
     &               J = ICHAR(FGNAME) - 64
                  IF(ICHAR(FGNAME) >= 97 .AND. ICHAR(FGNAME) <= 122)
     &               J = ICHAR(FGNAME) - 70
                  IF(ADJUSTL(FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24)) ==
     &               'WETTED') ITYPEG(NSPG(NSP)) = 1             
               ENDIF
            ENDDO   
         ENDDO
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(ITYPEG,ITYPEGW,NBCSG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         ITYPEGW = ITYPEG
      ENDIF

      DEALLOCATE(ITYPEG)
         
C ... Send the patch dimensions to all processes
  
      DO IPG = 1,NBCSG                ! Global patch loop            

         IF(NSPI == NSPTOTG) CYCLE    ! All patches found

         IF(ITYPEGW(IPG) == 0) CYCLE  ! Not wetted patch

         DO NSP = 1,NSPTOT            ! Local patch loop

            IF(IPG == NSPG(NSP)) THEN

               IDIMT(NSPI+1)  = IDIMS(NSP)
               JDIMT(NSPI+1)  = JDIMS(NSP)
               ITYPET(NSPI+1) = ICON((NSPP(NSP)-1)*IC9+1)

            ENDIF

         ENDDO 

         NSPI = NSPI + 1

      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(IDIMT,IDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(JDIMT,JDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(ITYPET,ITYPEW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         IDIMSW = IDIMT
         JDIMSW = JDIMT
         ITYPEW = ITYPET
      ENDIF

      DEALLOCATE(IDIMT,JDIMT,ITYPET)        


C ... Transfer the surface patch coordinates to the master process

      fnelems = 0
      fnnodes = 0

      DO ISP = 1,NSPTOTG

         IF(SOLIDL(ITYPEW(ISP))) THEN  ! Only solid like patches
            fnelems = fnelems +  IDIMSW(ISP)   * JDIMSW(ISP) 
            fnnodes = fnnodes + (IDIMSW(ISP)+1)*(JDIMSW(ISP)+1) 
         ENDIF

      ENDDO

      DEALLOCATE(IDIMSW,JDIMSW,ITYPEW,ITYPEGW)        
 
      RETURN
      END SUBROUTINE SIZE_FOR_MPCCI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MESH_FOR_MPCCI(XTYPE,INRE,JNRE,KNRE)

C ... Collect mesh data for multiphysics simulations, e.g. for MpCCI or
C     Modal FSI
C ... XTYPE = 1: Grid from the beginning of the simulation (XORI,YORI,ZORI)
C ... XTYPE = 2: Use deformed grid in restart (XCO,YCO,ZCO)
C ... XTYPE = 3: Same as XTYPE = 1, but no coordinate transformations
C ... XTYPE = 4: Same as XTYPE = 2, but no coordinate transformations

      USE MPI

      USE MPCCIVARS

      USE CHARACTERS

      USE FORCE_GROUPS

      USE INTEGERS,    ONLY : IPRO, NBCS

      USE MAIN_ARRAYS, ONLY : XCO, YCO, ZCO, XORI, YORI, ZORI, 
     &    ISTRS, JSTRS, KSTRS, NPROCE, NCHIMT, 
     &    ICON, IG, KX1S, KX2S, KY1S, KY2S, KZ1S, KZ2S,
     &    NSPB, NSPG, NSPP, IDIMS, JDIMS         

      USE NS3CO, ONLY : IC9, LN, NSPTOT, PARALLEL, FEMXYZ

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SOLIDL

      INTEGER :: NSP, I, J, K, L, N, M, ISP, IG1, IF1, KX1, KX2, KY1,
     &    KY2, KZ1, KZ2, KX1W, KX2W, KY1W, KY2W, KZ1W, KZ2W, LA,
     &    IPA, IPT, ISTR, JSTR, KSTR, IN, JN, KN, INRE, JNRE, KNRE,
     &    IWALL, IPSTR, JPSTR, NN, IERR, II, XTYPE, NG, CHT, ROTXYZ

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMSW, JDIMSW, ITYPEW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMT,  JDIMT,  ITYPET
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPEG, ITYPEGW
      INTEGER :: IPG, NBCSG, NSPTOTG, NSPI, IG1G
      INTEGER :: NSPOINTS, NSPOINTSG, NSPOINTSS

      REAL, ALLOCATABLE, DIMENSION(:) :: XCOW, YCOW, ZCOW
      REAL, ALLOCATABLE, DIMENSION(:) :: XCOT, YCOT, ZCOT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fnodechimtw, fnodechimtt

      LOGICAL :: MASTER, SLAVE

      CHARACTER(LEN=1) :: FGNAME

      MASTER = IPRO == 1
      SLAVE  = .NOT.MASTER


      IF(PARALLEL) THEN 

C ... Total number of boundary condition patches from all processes (NBCSG)
         CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,IERR)
         
C ... Total number of selected surface patches from all processes (NSPTOTG)
         CALL MPI_ALLREDUCE(NSPTOT,NSPTOTG,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,IERR)

      ELSE
         NBCSG   = NBCS
         NSPTOTG = NSPTOT
      ENDIF

      NSPI = 0

      ALLOCATE(IDIMSW(NSPTOTG),JDIMSW(NSPTOTG),ITYPEW(NSPTOTG))
      ALLOCATE(IDIMT(NSPTOTG), JDIMT(NSPTOTG), ITYPET(NSPTOTG))

      IDIMSW = 0; JDIMSW = 0; ITYPEW = 0
      IDIMT  = 0; JDIMT  = 0; ITYPET = 0
      
      ALLOCATE(ITYPEG(NBCSG),ITYPEGW(NBCSG))

      ITYPEG = 0; ITYPEGW = 0

      DO IPG = 1,NBCSG
         DO NSP = 1,NSPTOT
            DO K=1,80
               IF(BOUNDF(NSPP(NSP))(K:K) /= ' ') THEN
                  FGNAME = BOUNDF(NSPP(NSP))(K:K)
                  IF(ICHAR(FGNAME) >= 65 .AND. ICHAR(FGNAME) <= 90)
     &               J = ICHAR(FGNAME) - 64
                  IF(ICHAR(FGNAME) >= 97 .AND. ICHAR(FGNAME) <= 122)
     &               J = ICHAR(FGNAME) - 70
                  IF(ADJUSTL(FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24)) ==
     &               'WETTED') ITYPEG(NSPG(NSP)) = 1             
               ENDIF
            ENDDO   
         ENDDO
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(ITYPEG,ITYPEGW,NBCSG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         ITYPEGW = ITYPEG
      ENDIF

      DEALLOCATE(ITYPEG)
      
         
C ... Send the patch dimensions to all processes

      DO IPG = 1,NBCSG                ! Global patch loop

         IF(NSPI == NSPTOTG) CYCLE    ! All patches found

         IF(ITYPEGW(IPG) == 0) CYCLE  ! Not wetted patch

         DO NSP = 1,NSPTOT            ! Local patch loop
            
            IF(IPG == NSPG(NSP)) THEN               
               IDIMT(NSPI+1)  = IDIMS(NSP)
               JDIMT(NSPI+1)  = JDIMS(NSP)
               ITYPET(NSPI+1) = ICON((NSPP(NSP)-1)*IC9+1)               
            ENDIF

         ENDDO 

         NSPI = NSPI + 1
            
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(IDIMT,IDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(JDIMT,JDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(ITYPET,ITYPEW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         IDIMSW = IDIMT
         JDIMSW = JDIMT
         ITYPEW = ITYPET
      ENDIF

      DEALLOCATE(IDIMT,JDIMT,ITYPET)        


C ... Transfer the surface patch coordinates to the master process

      IG1G      = 1
      NSPI      = 0
      NSPOINTSG = 0 

      DO ISP = 1,NSPTOTG
         NSPOINTSG = NSPOINTSG + (IDIMSW(ISP)+1)*(JDIMSW(ISP)+1)
      ENDDO

      ALLOCATE(XCOW(NSPOINTSG),YCOW(NSPOINTSG),ZCOW(NSPOINTSG))
      ALLOCATE(XCOT(NSPOINTSG),YCOT(NSPOINTSG),ZCOT(NSPOINTSG))
      ALLOCATE(fnodechimtw(NSPOINTSG), fnodechimtt(NSPOINTSG))

      XCOW = 0.0; YCOW = 0.0; ZCOW = 0.0 
      XCOT = 0.0; YCOT = 0.0; ZCOT = 0.0 
      fnodechimtw = 0; fnodechimtt = 0 

      DO IPG = 1,NBCSG                 ! Global patch loop

C ... We know the size of the patch to be found next even when it is not 
C ... in this process.

         IF(NSPI == NSPTOTG) CYCLE     ! All patches found

         IF(ITYPEGW(IPG) == 0) CYCLE   ! Not wetted patch
 
         NSPOINTS = (IDIMSW(NSPI+1)+1)*(JDIMSW(NSPI+1)+1)

         DO NSP = 1,NSPTOT

            IF(IPG == NSPG(NSP)) THEN  ! Global patch number

               N     = NSPB(NSP)       ! Local block number
               NG    = NPROCE(1+N,IPRO)! Global block number
               CHT   = NCHIMT(NG)      ! Chimera priority

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
                        IF(XTYPE == 1 .OR. XTYPE == 3) THEN
                           XCOT(IG1G+IPT) = XORI(IG1+L)
                           YCOT(IG1G+IPT) = YORI(IG1+L)
                           ZCOT(IG1G+IPT) = ZORI(IG1+L)
                           fnodechimtt(IG1G+IPT) = CHT
                        ELSEIF(XTYPE == 2 .OR. XTYPE == 4) THEN
                           XCOT(IG1G+IPT) = XCO(IG1+L)
                           YCOT(IG1G+IPT) = YCO(IG1+L)
                           ZCOT(IG1G+IPT) = ZCO(IG1+L)
                           fnodechimtt(IG1G+IPT) = CHT
                        ENDIF
                        IPT = IPT + 1
                     ENDDO
                  ENDDO
               ENDDO
                
            ENDIF

         ENDDO 

         NSPI = NSPI + 1
         IG1G = IG1G + NSPOINTS               

      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(XCOT,XCOW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(YCOT,YCOW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(ZCOT,ZCOW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(fnodechimtt,fnodechimtw,NSPOINTSG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         XCOW = XCOT
         YCOW = YCOT
         ZCOW = ZCOT
         fnodechimtw = fnodechimtt
      ENDIF
      
      DEALLOCATE(XCOT,YCOT,ZCOT,fnodechimtt)

C ... Save mesh data for MpCCI

      IG1G = 1
      L    = 0
      M    = 0
      N    = 0

      ROTXYZ = FEMXYZ

      IF(XTYPE == 3 .OR. XTYPE == 4) ROTXYZ = 0  ! No XYZ transformation 

      
      DO NSPI=1,NSPTOTG

         II = (IDIMSW(NSPI)+1)*(JDIMSW(NSPI)+1)

         IF(SOLIDL(ITYPEW(NSPI))) THEN  ! Only solid like patches  

            SELECT CASE(ROTXYZ)

            CASE(-1)  ! CFD coordinates xyz correspond to xz-y in FEM model

               DO I=IG1G,IG1G+II-1
                  L = L + 1
                  fnodecoor(L,1) =  XCOW(I)
                  fnodecoor(L,3) =  YCOW(I)
                  fnodecoor(L,2) = -ZCOW(I)
                  fnodechimt(L)  =  fnodechimtw(I)
               ENDDO

            CASE(0)  ! Same coordinate system in CFD and FEM

               DO I=IG1G,IG1G+II-1
                  L = L + 1
                  fnodecoor(L,1) = XCOW(I)
                  fnodecoor(L,2) = YCOW(I)
                  fnodecoor(L,3) = ZCOW(I)
                  fnodechimt(L)  = fnodechimtw(I) 
               ENDDO

            CASE(1)  ! CFD coordinates xyz correspond to yzx in FEM model (HN)

               DO I=IG1G,IG1G+II-1
                  L = L + 1
                  fnodecoor(L,2) = XCOW(I)
                  fnodecoor(L,3) = YCOW(I)
                  fnodecoor(L,1) = ZCOW(I)
                  fnodechimt(L)  = fnodechimtw(I) 
               ENDDO

            CASE(2)  ! CFD coordinates xyz correspond to zxy in FEM model

               DO I=IG1G,IG1G+II-1
                  L = L + 1
                  fnodecoor(L,3) = XCOW(I)
                  fnodecoor(L,1) = YCOW(I)
                  fnodecoor(L,2) = ZCOW(I)
                  fnodechimt(L)  = fnodechimtw(I)
               ENDDO

            END SELECT
               
            M = L - (IDIMSW(NSPI)+1)*(JDIMSW(NSPI)+1)
                       
            DO J=1,JDIMSW(NSPI)
               M = M + 1
               DO I=1,IDIMSW(NSPI)
                  N = N + 1             ! Element number
                  felemnodes(N,1) = M
                  felemnodes(N,2) = M + 1
                  felemnodes(N,3) = felemnodes(N,2)+IDIMSW(NSPI)+1 
                  felemnodes(N,4) = felemnodes(N,3) - 1
                  M = M + 1             ! Node number 
               ENDDO
            ENDDO
             
         ENDIF

         IG1G = IG1G + II           
         
      ENDDO

C ... MpCCI mesh data ready

      DEALLOCATE(IDIMSW,JDIMSW,ITYPEW,ITYPEGW)        
      DEALLOCATE(XCOW,YCOW,ZCOW)       
      DEALLOCATE(fnodechimtw)       

      RETURN
      END SUBROUTINE MESH_FOR_MPCCI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DATA_FOR_MPCCI(INRE,JNRE,KNRE)

C ... Collect flow data for multiphysics simulations, e.g. for MpCCI or
C ... Modal FSI.      

      USE MPI

      USE MPCCIVARS

      USE CHARACTERS

      USE FORCE_GROUPS

      USE INTEGERS,    ONLY : IPRO, NBCS

      USE MAIN_ARRAYS, ONLY : SURFX, SURFY, SURFZ,
     &    TWALL, CPWALL, A1, A2, A3, ISTRS, JSTRS, KSTRS, IHF,
     &    ICON, IG, KX1S, KX2S, KY1S, KY2S, KZ1S, KZ2S,
     &    NSPB, NSPG, NSPP, IDIMS, JDIMS         

      USE NS3CO, ONLY : IC9, LN, NSPTOT, PARALLEL, ICYCLE, FEMXYZ

      IMPLICIT NONE

      LOGICAL, EXTERNAL :: SOLIDL

      INTEGER :: NSP, I, J, K, L, N, M, ISP, IG1, IF1, KX1, KX2, KY1,
     &    KY2, KZ1, KZ2, KX1W, KX2W, KY1W, KY2W, KZ1W, KZ2W, LA,
     &    IPA, IPT, ISTR, JSTR, KSTR, IN, JN, KN, INRE, JNRE, KNRE,
     &    IWALL, IPSTR, JPSTR, NN, IERR, II

      INTEGER :: STATUS(MPI_STATUS_SIZE)

      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMSW, JDIMSW, ITYPEW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIMT,  JDIMT,  ITYPET
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPEG, ITYPEGW
      INTEGER :: IPG, NBCSG, NSPTOTG, NSPI, IG1G
      INTEGER :: NSPOINTS, NSPOINTSG, NSPOINTSS

      REAL, ALLOCATABLE, DIMENSION(:) :: SURFXW, SURFYW, SURFZW
      REAL, ALLOCATABLE, DIMENSION(:) :: TWALLW, AREAW,  PRESW

      REAL, ALLOCATABLE, DIMENSION(:) :: SURFXT, SURFYT, SURFZT
      REAL, ALLOCATABLE, DIMENSION(:) :: TWALLT, AREAT,  PREST

      LOGICAL :: MASTER, SLAVE

      CHARACTER(LEN=1) :: FGNAME

      MASTER = IPRO == 1
      SLAVE  = .NOT.MASTER

      IF(PARALLEL) THEN 

C ... Total number of boundary condition patches from all processes (NBCSG)
         CALL MPI_ALLREDUCE(NBCS,NBCSG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)
         
C ... Total number of selected surface patches from all processes (NSPTOTG)
         CALL MPI_ALLREDUCE(NSPTOT,NSPTOTG,1,MPI_INTEGER,MPI_SUM,
     &                      MPI_COMM_WORLD,IERR)

      ELSE
         NBCSG   = NBCS
         NSPTOTG = NSPTOT
      ENDIF
         
      NSPI = 0

      ALLOCATE(IDIMSW(NSPTOTG),JDIMSW(NSPTOTG),ITYPEW(NSPTOTG))
      ALLOCATE(IDIMT(NSPTOTG), JDIMT(NSPTOTG), ITYPET(NSPTOTG))

      IDIMSW = 0; JDIMSW = 0; ITYPEW = 0
      IDIMT  = 0; JDIMT  = 0; ITYPET = 0

      ALLOCATE(ITYPEG(NBCSG),ITYPEGW(NBCSG))

      ITYPEG = 0; ITYPEGW = 0

      DO IPG = 1,NBCSG
         DO NSP = 1,NSPTOT
            DO K=1,80
               IF(BOUNDF(NSPP(NSP))(K:K) /= ' ') THEN
                  FGNAME = BOUNDF(NSPP(NSP))(K:K)
                  IF(ICHAR(FGNAME) >= 65 .AND. ICHAR(FGNAME) <= 90)
     &               J = ICHAR(FGNAME) - 64
                  IF(ICHAR(FGNAME) >= 97 .AND. ICHAR(FGNAME) <= 122)
     &               J = ICHAR(FGNAME) - 70
                  IF(ADJUSTL(FORCE_GROUP_FULL_NAME((J-1)*24+1:J*24)) ==
     &               'WETTED') ITYPEG(NSPG(NSP)) = 1             
               ENDIF
            ENDDO   
         ENDDO
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(ITYPEG,ITYPEGW,NBCSG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         ITYPEGW = ITYPEG
      ENDIF

      DEALLOCATE(ITYPEG)

C ... Send the patch dimensions to all processes

      DO IPG = 1,NBCSG                ! Global patch loop

         IF(NSPI == NSPTOTG) CYCLE    ! All patches found

         IF(ITYPEGW(IPG) == 0) CYCLE  ! Not wetted patch

         DO NSP = 1,NSPTOT            ! Local patch loop

            IF(IPG == NSPG(NSP)) THEN
               IDIMT(NSPI+1)  = IDIMS(NSP)
               JDIMT(NSPI+1)  = JDIMS(NSP)
               ITYPET(NSPI+1) = ICON((NSPP(NSP)-1)*IC9+1)
            ENDIF

         ENDDO 
         
         NSPI = NSPI + 1
            
      ENDDO

      IF(PARALLEL) THEN
         CALL MPI_ALLREDUCE(IDIMT,IDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(JDIMT,JDIMSW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
         CALL MPI_ALLREDUCE(ITYPET,ITYPEW,NSPTOTG,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         IDIMSW = IDIMT
         JDIMSW = JDIMT
         ITYPEW = ITYPET
      ENDIF

      DEALLOCATE(IDIMT,JDIMT,ITYPET)        

        
      IG1G      = 1
      NSPI      = 0
      NSPOINTSG = 0 

      DO ISP = 1,NSPTOTG
         NSPOINTSG = NSPOINTSG + IDIMSW(ISP)*JDIMSW(ISP)
      ENDDO
      
      ALLOCATE(SURFXW(NSPOINTSG),SURFYW(NSPOINTSG),
     &         SURFZW(NSPOINTSG),TWALLW(NSPOINTSG),
     &         AREAW(NSPOINTSG),PRESW(NSPOINTSG))
      SURFXW=0.0; SURFYW=0.0; SURFZW=0.0
      TWALLW=0.0; AREAW=0.0;  PRESW=0.0

      ALLOCATE(SURFXT(NSPOINTSG),SURFYT(NSPOINTSG),SURFZT(NSPOINTSG),
     &         TWALLT(NSPOINTSG),AREAT (NSPOINTSG),PREST (NSPOINTSG))

      SURFXT=0.0; SURFYT=0.0; SURFZT=0.0
      TWALLT=0.0; AREAT=0.0;  PREST=0.0


      DO 999 IPG = 1,NBCSG            ! Global patch loop

C ... We know the size of the patch to be found next even when it  
C ... is not in this process.

         IF(NSPI == NSPTOTG) CYCLE    ! All patches found

         IF(ITYPEGW(IPG) == 0) CYCLE  ! Not wetted patch

         NSPOINTS = IDIMSW(NSPI+1)*JDIMSW(NSPI+1)

         IPA = 0
         IPT = 0

         DO 997 NSP = 1,NSPTOT

            IF(IPG /= NSPG(NSP)) CYCLE

            N     = NSPB(NSP)   ! Local block number
            ISP   = NSPP(NSP)   ! Local patch number
            IWALL = ICON((ISP-1)*IC9 + 3)

            IF1 = IHF(ISP,1)
            IG1 = IG(1,N)
            
            IF(IWALL == 2 .OR. IWALL == 5) THEN
               KY1   = KX1S(NSP)
               KY2   = KX2S(NSP)
               KX1   = KY1S(NSP)
               KX2   = KY2S(NSP)
               IPSTR = KX2 - KX1 + 1 + 2*LN
               JPSTR = 1
               NN    = (LN-KY1)*IPSTR - KX1 + LN + IF1
            ENDIF
            
            KX1 = KX1S(NSP)
            KX2 = KX2S(NSP)
            KY1 = KY1S(NSP)
            KY2 = KY2S(NSP) 
            
            IF(IWALL /= 2 .AND. IWALL /= 5) THEN
               JPSTR = KX2 - KX1 + 1 + 2*LN
               IPSTR = 1
               NN    = (LN-KY1)*JPSTR - KX1 + LN + IF1
            ENDIF
            
            DO J=KY1,KY2
               DO I=KX1,KX2
                  L   = I*IPSTR+J*JPSTR+NN
                  SURFXT(IG1G+IPT) = SURFX(L)
                  SURFYT(IG1G+IPT) = SURFY(L)
                  SURFZT(IG1G+IPT) = SURFZ(L)
                  TWALLT(IG1G+IPT) = TWALL(L)
                  PREST(IG1G+IPT)  = CPWALL(L)
                  IPT = IPT + 1
               ENDDO
            ENDDO

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
            
            DO K=KZ1W,KZ2W
               DO J=KY1W,KY2W
                  DO I=KX1W,KX2W
                     LA=(I+IN-1)*ISTR+(J+JN-1)*JSTR+(K+KN-1)*KSTR
                     LA = IG1 + LA
                     IF(IWALL==1.OR.IWALL==4) AREAT(IG1G+IPA)=A1(LA)
                     IF(IWALL==2.OR.IWALL==5) AREAT(IG1G+IPA)=A2(LA)
                     IF(IWALL==3.OR.IWALL==6) AREAT(IG1G+IPA)=A3(LA)
                     IPA = IPA + 1
                  ENDDO
               ENDDO
            ENDDO
            
 997     CONTINUE

         NSPI = NSPI + 1
         IG1G = IG1G + NSPOINTS

 999  CONTINUE

      IF(PARALLEL) THEN
         CALL MPI_REDUCE(SURFXT,SURFXW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SURFYT,SURFYW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SURFZT,SURFZW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(TWALLT,TWALLW,NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(AREAT,AREAW,  NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(PREST,PRESW,  NSPOINTSG,
     &        MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      ELSE
         SURFXW = SURFXT
         SURFYW = SURFYT
         SURFZW = SURFZT
         TWALLW = TWALLT
         AREAW  = AREAT
         PRESW  = PREST
      ENDIF
         
      DEALLOCATE(SURFXT,SURFYT,SURFZT,TWALLT,AREAT,PREST)


      IF(MASTER) THEN
      
      M = 0
      N = 0

      DO NSPI=1,NSPTOTG
         
         II = IDIMSW(NSPI)*JDIMSW(NSPI)

         IF(SOLIDL(ITYPEW(NSPI))) THEN  ! Only solid like patches

            L = 0

            SELECT CASE(FEMXYZ)

            CASE(-1)  ! CFD coordinates xyz correspond to xz-y in FEM model

            DO J=1,JDIMSW(NSPI)
               DO I=1,IDIMSW(NSPI)

                  L = L + 1 
                  N = N + 1

                  felemforce(N,1) =  SURFXW(L+M)*AREAW(L+M)
                  felemforce(N,3) =  SURFYW(L+M)*AREAW(L+M)
                  felemforce(N,2) = -SURFZW(L+M)*AREAW(L+M)

                  felemtemp(N)    = TWALLW(L+M)
                  felempres(N)    = PRESW(L+M)
                  felemarea(N)    = AREAW(L+M)

               ENDDO
            ENDDO

            CASE(0)  ! Same coordinate system in CFD and FEM

            DO J=1,JDIMSW(NSPI)
               DO I=1,IDIMSW(NSPI)

                  L = L + 1 
                  N = N + 1

                  felemforce(N,1) = SURFXW(L+M)*AREAW(L+M)
                  felemforce(N,2) = SURFYW(L+M)*AREAW(L+M)
                  felemforce(N,3) = SURFZW(L+M)*AREAW(L+M)

                  felemtemp(N)    = TWALLW(L+M)
                  felempres(N)    = PRESW(L+M)
                  felemarea(N)    = AREAW(L+M)

               ENDDO
            ENDDO

            CASE(1)  ! CFD coordinates xyz correspond to yzx in FEM model (HN)

            DO J=1,JDIMSW(NSPI)
               DO I=1,IDIMSW(NSPI)

                  L = L + 1 
                  N = N + 1

                  felemforce(N,2) = SURFXW(L+M)*AREAW(L+M)
                  felemforce(N,3) = SURFYW(L+M)*AREAW(L+M)
                  felemforce(N,1) = SURFZW(L+M)*AREAW(L+M)

                  felemtemp(N)    = TWALLW(L+M)
                  felempres(N)    = PRESW(L+M)
                  felemarea(N)    = AREAW(L+M)

               ENDDO
            ENDDO

            CASE(2)  ! CFD coordinates xyz correspond to zxy in FEM model

            DO J=1,JDIMSW(NSPI)
               DO I=1,IDIMSW(NSPI)

                  L = L + 1 
                  N = N + 1

                  felemforce(N,3) = SURFXW(L+M)*AREAW(L+M)
                  felemforce(N,1) = SURFYW(L+M)*AREAW(L+M)
                  felemforce(N,2) = SURFZW(L+M)*AREAW(L+M)

                  felemtemp(N)    = TWALLW(L+M)
                  felempres(N)    = PRESW(L+M)
                  felemarea(N)    = AREAW(L+M)

               ENDDO
            ENDDO

            END SELECT
            
         ENDIF

         M = M + II           
         
      ENDDO

      ENDIF  ! MASTER
      
      DEALLOCATE(IDIMSW,JDIMSW,ITYPEW,ITYPEGW)        
      DEALLOCATE(SURFXW,SURFYW,SURFZW,TWALLW,AREAW,PRESW) 

      RETURN
      END SUBROUTINE DATA_FOR_MPCCI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REMESH

C ... Deform the grid

      USE MPCCIVARS

      USE MAIN_ARRAYS, ONLY : IMAX, JMAX, KMAX, IG, XCO, YCO, ZCO,
     &            XC, YC, ZC, PDIFF, XORI, YORI, ZORI, NPROCE,
     &            SKINX, SKINY, SKINZ, IBSKIN, DISTP, IDP 

      USE NS3CO

      USE INTEGERS, ONLY : IPRO

      USE CONSTANTS, ONLY : PII, EPS6

      IMPLICIT NONE   

      INTEGER, PARAMETER :: ISIZE = 101, JSIZE = 101, KSIZE = 101
*      INTEGER, PARAMETER :: ISIZE = 61, JSIZE = 11, KSIZE = 61
*      INTEGER, PARAMETER :: ISIZE = 181, JSIZE = 61, KSIZE = 181


      REAL, DIMENSION(ISIZE*JSIZE*KSIZE) :: XT, YT, ZT
      REAL, DIMENSION(ISIZE*JSIZE*KSIZE) :: DX, DY, DZ
      INTEGER, DIMENSION((ISIZE-1)*(JSIZE-1)*(KSIZE-1)) :: ND, NC

      REAL, DIMENSION((ISIZE-1)*(JSIZE-1)*(KSIZE-1))
     &               :: PXV, PYV, PZV, DXV, DYV, DZV
      REAL, DIMENSION(:,:), ALLOCATABLE  :: DEFDAT

      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX 
      REAL :: XMIND, YMIND, ZMIND, XMAXD, YMAXD, ZMAXD  
      REAL :: SX, SY, SZ, SS, SC, XTAVE, YTAVE, ZTAVE

      INTEGER :: ISZ, JSZ, KSZ, KF, KL
      INTEGER :: IANS, NB, IB, NNODE, LNODE, INODE
      INTEGER :: NI, NJ, NK, I, J, K, N, IG1, NDEFORM
      INTEGER :: IL, IL1, ILF, ILL 
      INTEGER :: IL2, IL3, IL4, IL5, IL6, IL7, IL8, IL9, IL10,IL11, IL12

      REAL :: D1, D2, D3, D4, D5, D6, RELX, DEFMAX, DEFMAXE, RELXDEFMAX
      REAL :: DEFMAXX, DEFMAXY, DEFMAXZ
      REAL :: XDEFMAX, YDEFMAX, ZDEFMAX
      
      INTEGER, DIMENSION(0:12) :: ILP
      INTEGER :: NEXT, IT

      INTEGER, DIMENSION(8) :: II
      REAL :: X0(3), A0(4), VOL

      LOGICAL :: MASTER, TWODIM

      MASTER = IPRO == 1

      ISZ = ISIZE
      JSZ = JSIZE
      KSZ = KSIZE
      
      IF (MPCCIL) THEN
         NDEFORM = fnnodes
      ELSEIF (MODALFSIL) THEN
         NDEFORM = njoints
      ENDIF

C ... Two-dimensional case ?

      TWODIM = .FALSE.
      DO N = 1,NBLOCK
         IF(KMAX(1,N) == 2) TWODIM = .TRUE.
      ENDDO      

      IF(TWODIM) KSZ = 8

      ALLOCATE(DEFDAT(NDEFORM,6))
      
      defmax  = 0.0
      defmaxe = 0.0
      defmaxx = 0.0
      defmaxy = 0.0
      defmaxz = 0.0

      NDEFORM = 0

      IF (MPCCIL) THEN

         DO I=1,fnnodes

            D1 = fnodecoor(I,1)
            D2 = fnodecoor(I,2)
            D3 = fnodecoor(I,3)
            D4 = fnodedisp(3*(I-1)+1)
            D5 = fnodedisp(3*(I-1)+2)
            D6 = fnodedisp(3*(I-1)+3)

C ... Remove zero displacements in voxel based remeshing

*            if(D4 == 0.0 .and. D5 == 0.0 .and. D6 == 0.0) cycle

            NDEFORM = NDEFORM + 1

            DEFDAT(NDEFORM,1) = D1
            DEFDAT(NDEFORM,2) = D2
            DEFDAT(NDEFORM,3) = D3
            DEFDAT(NDEFORM,4) = D4
            DEFDAT(NDEFORM,5) = D5
            DEFDAT(NDEFORM,6) = D6

*            if(DEFDAT(NDEFORM,3) > 1.5) DEFDAT(NDEFORM,5) = 1.1*D5  ! Hornet
*            if(DEFDAT(NDEFORM,3) > 1.5) DEFDAT(NDEFORM,3) = 1.1*D3  ! Hornet 

            defmax = SQRT(D4**2 + D5**2 + D6**2)
         
            if(defmax > defmaxe) THEN
               defmaxe = defmax
               xdefmax = D1
               ydefmax = D2
               zdefmax = D3
               defmaxx = D4
               defmaxy = D5
               defmaxz = D6
            endif
            
         ENDDO

      ELSEIF (MODALFSIL) THEN

         DO I=1,njoints

            D1 = jointxyzo(1,I)
            D2 = jointxyzo(2,I)
            D3 = jointxyzo(3,I)
            D4 = jointdisp(1,I)
            D5 = jointdisp(2,I)
            D6 = jointdisp(3,I)

            NDEFORM = NDEFORM + 1

            DEFDAT(NDEFORM,1) = D1
            DEFDAT(NDEFORM,2) = D2
            DEFDAT(NDEFORM,3) = D3
            DEFDAT(NDEFORM,4) = D4
            DEFDAT(NDEFORM,5) = D5
            DEFDAT(NDEFORM,6) = D6

            defmax = SQRT(D4**2 + D5**2 + D6**2)
         
            if(defmax > defmaxe) THEN
               defmaxe = defmax
               xdefmax = D1
               ydefmax = D2
               zdefmax = D3
               defmaxx = D4
               defmaxy = D5
               defmaxz = D6
            endif

         ENDDO
         
      ENDIF

*      RELX = 1.0
*      RELX = MIN(REAL(ICYCLE-1000)/REAL(2000),1.0)
      RELX = 1.0

      IF (MASTER) THEN
         RELXDEFMAX = RELX*defmax
************************************************************************
*         if(MODALFSIL) then
            write(678,*) REAL(T,4),REAL(jointdisp(2,221),4),
     &                   REAL(jointdisp(2,231),4)
*         elseif(MPCCIL) then
*            write(678,*) REAL(T,4),REAL(fnodedisp(3*(2290-1)+2),4),
*     &                   REAL(fnodedisp(3*(4314-1)+2),4)
*         endif
************************************************************************
         if(timel) then
            write(*,*) 'ITIMES   ',ITIMES
         else
            write(*,*) 'ICYCLE   ',ICYCLE
         endif
         write(*,*) 'defmaxe, RELX, defmaxe*RELX',
     &               REAL(defmaxe,4),REAL(RELX,4),REAL(RELXDEFMAX,4)
         write(*,*) 'defmaxx, defmaxy, defmaxz',
     &               REAL(defmaxx,4),REAL(defmaxy,4),REAL(defmaxz,4)
         write(*,*) 'xdefmax, ydefmax, zdefmax',
     &               REAL(xdefmax,4),REAL(ydefmax,4),REAL(zdefmax,4)
      ENDIF


************************************************************************
* Relax the deflections if needed ?
************************************************************************
      DO I=1,NDEFORM
         DEFDAT(I,4) = DEFDAT(I,4)*RELX 
         DEFDAT(I,5) = DEFDAT(I,5)*RELX 
         DEFDAT(I,6) = DEFDAT(I,6)*RELX
      ENDDO
************************************************************************


      XMIN =  1.0E+30
      YMIN =  1.0E+30
      ZMIN =  1.0E+30
      XMAX = -1.0E+30
      YMAX = -1.0E+30
      ZMAX = -1.0E+30

      DO I = 1,fnnodes
         IF(fnodecoor_orig(I,1) < XMIN) XMIN = fnodecoor_orig(I,1)
         IF(fnodecoor_orig(I,2) < YMIN) YMIN = fnodecoor_orig(I,2)
         IF(fnodecoor_orig(I,3) < ZMIN) ZMIN = fnodecoor_orig(I,3)
         IF(fnodecoor_orig(I,1) > XMAX) XMAX = fnodecoor_orig(I,1)
         IF(fnodecoor_orig(I,2) > YMAX) YMAX = fnodecoor_orig(I,2)
         IF(fnodecoor_orig(I,3) > ZMAX) ZMAX = fnodecoor_orig(I,3)
      ENDDO


C ... Boundaries of all skin elements, not just the wetted elements. 

*      XMIN = MINVAL(SKINX)
*      XMAX = MAXVAL(SKINX)
*      YMIN = MINVAL(SKINY) 
*      YMAX = MAXVAL(SKINY)
*      ZMIN = MINVAL(SKINZ) 
*      ZMAX = MAXVAL(SKINZ)


C ... Equally spaced rectangular interpolation grid around the known 
C ... deformations. The interpolation grid must be enlarged in all 
C ... directions in order to give space for propagating and fading out 
C ... the deformations. External flow assumed.

      SX = XMAX-XMIN
      SY = YMAX-YMIN
      SZ = ZMAX-ZMIN
      SS = MAX(SX,SY,SZ)
      XMIN = XMIN - SS + EPS6
      XMAX = XMAX + SS + EPS6
      YMIN = YMIN - SS + EPS6
      YMAX = YMAX + SS + EPS6
      ZMIN = ZMIN - SS + EPS6
      ZMAX = ZMAX + SS + EPS6

*      if(master) then
*         write(*,*) xmin, xmax
*         write(*,*) ymin, ymax
*         write(*,*) zmin, zmax
*      endif

      SX   = (XMAX-XMIN)/(ISZ-1)
      SY   = (YMAX-YMIN)/(JSZ-1)
      SZ   = (ZMAX-ZMIN)/(KSZ-1)
       
C ... Generate a rectangular transformation grid and first fill it
C ... with zero transformations.
          
      DO K=1,KSZ
         DO J=1,JSZ
            DO I=1,ISZ
               IL = ISZ*JSZ*(K-1)+ISZ*(J-1)+I
               XT(IL) = XMIN + (I-1)*SX
               YT(IL) = YMIN + (J-1)*SY
               ZT(IL) = ZMIN + (K-1)*SZ
               DX(IL) = 0.0
               DY(IL) = 0.0
               DZ(IL) = 0.0
            ENDDO
         ENDDO
      ENDDO


C ... Zero the average deformations in transformation grid cells.      

      DO I=1,ISZ-1
         DO J=1,JSZ-1
            DO K=1,KSZ-1
               IL = (ISZ-1)*(JSZ-1)*(K-1)+(ISZ-1)*(J-1)+I
               ND(IL)  = 0
               DXV(IL) = 0.0
               DYV(IL) = 0.0
               DZV(IL) = 0.0
            ENDDO
         ENDDO
      ENDDO

      
C ... Search deformations within the transformation grid cells. 

      DO N = 1,NDEFORM

         I  = INT((DEFDAT(N,1)-XMIN)/SX) + 1
         J  = INT((DEFDAT(N,2)-YMIN)/SY) + 1
         K  = INT((DEFDAT(N,3)-ZMIN)/SZ) + 1
         IL = (ISZ-1)*(JSZ-1)*(K-1)+(ISZ-1)*(J-1)+I
         DXV(IL) = DXV(IL) + DEFDAT(N,4)
         DYV(IL) = DYV(IL) + DEFDAT(N,5) 
         DZV(IL) = DZV(IL) + DEFDAT(N,6)
         ND(IL)  = ND(IL)  + 1 
      ENDDO

      
C ... Average deformations within the transformation grid cells. 

      DO IL=1,(ISZ-1)*(JSZ-1)*(KSZ-1)

         IF(ND(IL) > 0) THEN
            DXV(IL) = DXV(IL)/ND(IL)
            DYV(IL) = DYV(IL)/ND(IL)
            DZV(IL) = DZV(IL)/ND(IL) 
            ND(IL) = 1
         ENDIF

      ENDDO

C ... Extrapolate the deformations

      IF(TWODIM) THEN
         KF = 1
         KL = 1
      ELSE
         KF = 3
         KL = KSZ-3
      ENDIF

      NEXT = 2  ! Number of extrapolations, i.e., the number of interpolation
                ! grid cells across which the deformations are extrapolated.

      DO N=1,NEXT  ! Repeat the extrapolation NEXT times

      NC = 0

      DO I=3,ISZ-3
         DO J=3,JSZ-3
            DO K=KF,KL

               IL   = (ISZ-1)*(JSZ-1)*(K-1)+(ISZ-1)*(J-1)+I
               IL1  = IL - 1
               IL2  = IL - 2
               IL3  = IL + 1
               IL4  = IL + 2
               IL5  = IL -   (ISZ-1)
               IL6  = IL - 2*(ISZ-1)
               IL7  = IL +   (ISZ-1)
               IL8  = IL + 2*(ISZ-1)

               IF(.NOT.TWODIM) THEN
                  IL9  = IL -   (ISZ-1)*(JSZ-1)
                  IL10 = IL - 2*(ISZ-1)*(JSZ-1)
                  IL11 = IL +   (ISZ-1)*(JSZ-1)
                  IL12 = IL + 2*(ISZ-1)*(JSZ-1)
               ENDIF

               IF(ND(IL) == 0 .AND. ND(IL1) > 0 .AND. ND(IL2) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL1) - DXV(IL2)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL1) - DYV(IL2)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL1) - DZV(IL2)
                  NC(IL)  = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL3) > 0 .AND. ND(IL4) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL3) - DXV(IL4)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL3) - DYV(IL4)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL3) - DZV(IL4)
                  NC(IL)  = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL5) > 0 .AND. ND(IL6) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL5) - DXV(IL6)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL5) - DYV(IL6)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL5) - DZV(IL6)
                  NC(IL)  = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL7) > 0 .AND. ND(IL8) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL7) - DXV(IL8)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL7) - DYV(IL8)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL7) - DZV(IL8)
                  NC(IL)  = NC(IL) + 1
               ENDIF

               IF(.NOT.TWODIM) THEN
               IF(ND(IL) == 0 .AND. ND(IL9) > 0 .AND.ND(IL10) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL9) - DXV(IL10)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL9) - DYV(IL10)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL9) - DZV(IL10)
                  NC(IL)  = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND.ND(IL11) > 0 .AND.ND(IL12) > 0) THEN
                  DXV(IL) = DXV(IL) + 2.0*DXV(IL11) - DXV(IL12)
                  DYV(IL) = DYV(IL) + 2.0*DYV(IL11) - DYV(IL12)
                  DZV(IL) = DZV(IL) + 2.0*DZV(IL11) - DZV(IL12)
                  NC(IL)  = NC(IL) + 1
               ENDIF
               ENDIF

               IF(NC(IL) > 0) THEN
                  DXV(IL) = DXV(IL) / NC(IL)
                  DYV(IL) = DYV(IL) / NC(IL)
                  DZV(IL) = DZV(IL) / NC(IL)
                  NC(IL)  = 1
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      DO IL=1,(ISZ-1)*(JSZ-1)*(KSZ-1)
         IF(NC(IL) == 1) ND(IL) = 1 
      ENDDO

      ENDDO

      
C ... Extend the deformations

      NEXT = 15   ! Number of extensions, i.e., the number of interpolation
                  ! grid cells across which the deformations are faded out.

      IF(TWODIM) THEN
         KF = 1
         KL = 1
      ELSE
         KF = 3
         KL = KSZ-3
      ENDIF

      DO N=1,NEXT

      NC = 0

      SC = SQRT(0.5*(COS((N-1+0.5)*PII/NEXT)+1.0))

      DO I=3,ISZ-3
         DO J=3,JSZ-3
            DO K=KF,KL

               IL   = (ISZ-1)*(JSZ-1)*(K-1)+(ISZ-1)*(J-1)+I
               IL1  = IL - 1
               IL2  = IL - 2
               IL3  = IL + 1
               IL4  = IL + 2
               IL5  = IL -   (ISZ-1)
               IL6  = IL - 2*(ISZ-1)
               IL7  = IL +   (ISZ-1)
               IL8  = IL + 2*(ISZ-1)

               IF(.NOT.TWODIM) THEN
                  IL9  = IL -   (ISZ-1)*(JSZ-1)
                  IL10 = IL - 2*(ISZ-1)*(JSZ-1)
                  IL11 = IL +   (ISZ-1)*(JSZ-1)
                  IL12 = IL + 2*(ISZ-1)*(JSZ-1)
               ENDIF

               IF(ND(IL) == 0 .AND. ND(IL1) > 0) THEN
                  DXV(IL) = DXV(IL) + DXV(IL1)
                  DYV(IL) = DYV(IL) + DYV(IL1)
                  DZV(IL) = DZV(IL) + DZV(IL1)
                  NC(IL) = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL3) > 0) THEN
                  DXV(IL) = DXV(IL) + DXV(IL3)
                  DYV(IL) = DYV(IL) + DYV(IL3)
                  DZV(IL) = DZV(IL) + DZV(IL3)
                  NC(IL) = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL5) > 0) THEN
                  DXV(IL) = DXV(IL) + DXV(IL5)
                  DYV(IL) = DYV(IL) + DYV(IL5)
                  DZV(IL) = DZV(IL) + DZV(IL5)
                  NC(IL) = NC(IL) + 1
               ENDIF
               IF(ND(IL) == 0 .AND. ND(IL7) > 0) THEN
                  DXV(IL) = DXV(IL) + DXV(IL7)
                  DYV(IL) = DYV(IL) + DYV(IL7)
                  DZV(IL) = DZV(IL) + DZV(IL7)
                  NC(IL) = NC(IL) + 1
               ENDIF

               IF(.NOT.TWODIM) THEN
                  IF(ND(IL) == 0 .AND. ND(IL9) > 0) THEN
                     DXV(IL) = DXV(IL) + DXV(IL9)
                     DYV(IL) = DYV(IL) + DYV(IL9)
                     DZV(IL) = DZV(IL) + DZV(IL9)
                     NC(IL) = NC(IL) + 1
                  ENDIF
                  IF(ND(IL) == 0 .AND.ND(IL11) > 0) THEN
                     DXV(IL) = DXV(IL) + DXV(IL11)
                     DYV(IL) = DYV(IL) + DYV(IL11)
                     DZV(IL) = DZV(IL) + DZV(IL11)
                     NC(IL) = NC(IL) + 1
                  ENDIF

               ENDIF

               IF(NC(IL) > 0) THEN
                  DXV(IL) = DXV(IL) / NC(IL) * SC
                  DYV(IL) = DYV(IL) / NC(IL) * SC
                  DZV(IL) = DZV(IL) / NC(IL) * SC
                  NC(IL) = 1
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      DO IL=1,(ISZ-1)*(JSZ-1)*(KSZ-1)
         IF(NC(IL) == 1) ND(IL) = 1 
      ENDDO

      ENDDO


C ... Average deformations to the transfer grid node points.

      IF(TWODIM) THEN
         KF = 1
         KL = 2
      ELSE
         KF = 2
         KL = KSZ-1
      ENDIF


      DO I=2,ISZ-1
         DO J=2,JSZ-1
            DO K=KF,KL

               IL  = ISZ*JSZ*(K-1)+ISZ*(J-1)+I

               IL1 = (ISZ-1)*(JSZ-1)*(K-2)+(ISZ-1)*(J-2)+I-1
               IL2 = IL1 + 1
               IL3 = IL1 + (ISZ-1)
               IL4 = IL3 + 1

               IF(TWODIM) THEN
                  IL  = ISZ*JSZ*(K-1)+ISZ*(J-1)+I
                  IL1 = (ISZ-1)*(J-2)+I-1
                  IL2 = IL1 + 1
                  IL3 = IL1 + (ISZ-1)
                  IL4 = IL3 + 1
               ENDIF

               IF(.NOT.TWODIM) THEN
                  IL5 = IL1 + (ISZ-1)*(JSZ-1)
                  IL6 = IL2 + (ISZ-1)*(JSZ-1)
                  IL7 = IL3 + (ISZ-1)*(JSZ-1)
                  IL8 = IL4 + (ISZ-1)*(JSZ-1)
               ENDIF

               IF(.NOT.TWODIM) THEN
                  ILF = ND(IL1) + ND(IL2) + ND(IL3) + ND(IL4)
     &                + ND(IL5) + ND(IL6) + ND(IL7) + ND(IL8)
               ELSE
                  ILF = ND(IL1) + ND(IL2) + ND(IL3) + ND(IL4)
               ENDIF

               IF(ILF > 0) THEN

               IF(.NOT.TWODIM) THEN
                  DX(IL) = (DXV(IL1) + DXV(IL2) + DXV(IL3) + DXV(IL4) + 
     &        +  DXV(IL5) + DXV(IL6) + DXV(IL7) + DXV(IL8)) / ILF
                  DY(IL) = (DYV(IL1) + DYV(IL2) + DYV(IL3) + DYV(IL4) + 
     &        +  DYV(IL5) + DYV(IL6) + DYV(IL7) + DYV(IL8)) / ILF
                  DZ(IL) = (DZV(IL1) + DZV(IL2) + DZV(IL3) + DZV(IL4) + 
     &        +  DZV(IL5) + DZV(IL6) + DZV(IL7) + DZV(IL8)) / ILF
               ELSE
                  DX(IL) = (DXV(IL1) + DXV(IL2) 
     &        + DXV(IL3) + DXV(IL4)) / ILF
                  DY(IL) = (DYV(IL1) + DYV(IL2) 
     &        + DYV(IL3) + DYV(IL4)) / ILF 
                  DZ(IL) = (DZV(IL1) + DZV(IL2) 
     &        + DZV(IL3) + DZV(IL4)) / ILF
               ENDIF

               ENDIF

            ENDDO
         ENDDO
      ENDDO

      IF(MODALFSIL) THEN
         CALL MAP_DISPLACEMENTS
      ENDIF

C ... Deform the grid. Two alternative ways to deform the volume grid.
C ... DEFORM is based on voxels and DEFORM_IDP is more accurate and based
C ... on movement of the closest (wetted) surface movement.
      
      DO N=1,NBLOCK

         IG1 = IG(1,N)

         IF(RMESHT == 2) THEN

         CALL DEFORM_IDP(N,IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     &        XCO(IG1),YCO(IG1),ZCO(IG1),XORI(IG1),YORI(IG1),ZORI(IG1),
     &        DISTP(IG1),IDP(IG1),DEFMAXX,DEFMAXY,DEFMAXZ)

         ELSE

         CALL DEFORM_VOX(IMAX(1,N),JMAX(1,N),KMAX(1,N),IN,JN,KN,
     &        XCO(IG1),YCO(IG1),ZCO(IG1),XORI(IG1),YORI(IG1),ZORI(IG1),
     &        XT,YT,ZT,DX,DY,DZ,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     &        SX,SY,SZ,ISZ,JSZ,KSZ)

      ENDIF

************************************************************************
* Example of recording movements of certain grid points
*
*     ElasticFlap movement (two grid points in block 1)
*
*      if(NPROCE(1+N,IPRO) == 1) THEN
*            I = 1
*            J = 1
*            K = 1
*            IL = (I-1+IN)+(J-1+JN)*(IMAX(1,N)+2*IN)
*     &           +(K-1+KN)*(IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)
*          write(*,*) REAL(T,4),REAL(XCO(IG1+IL),4),REAL(XORI(IG1+IL),4),
*     &               -REAL(XCO(IG1+IL)-XORI(IG1+IL),4),
*     &                REAL(YCO(IG1+IL)-YORI(IG1+IL),4),
*     &               -REAL(ZCO(IG1+IL)-ZORI(IG1+IL),4)
*        write(678,*) REAL(T,4),REAL(XCO(IG1+IL),4),REAL(XORI(IG1+IL),4),
*     &               -REAL(XCO(IG1+IL)-XORI(IG1+IL),4),
*     &                REAL(YCO(IG1+IL)-YORI(IG1+IL),4),
*     &               -REAL(ZCO(IG1+IL)-ZORI(IG1+IL),4)
*         ENDIF

*     Second point
         
*      if(NPROCE(1+N,IPRO) == 1) THEN
*            I = 49
*            J = 1
*            K = 1
*            IL = (I-1+IN)+(J-1+JN)*(IMAX(1,N)+2*IN)
*     &           +(K-1+KN)*(IMAX(1,N)+2*IN)*(JMAX(1,N)+2*JN)
*        write(679,*) REAL(T,4),REAL(XCO(IG1+IL),4),REAL(XORI(IG1+IL),4),
*     &               -REAL(XCO(IG1+IL)-XORI(IG1+IL),4),
*     &                REAL(YCO(IG1+IL)-YORI(IG1+IL),4),
*     &               -REAL(ZCO(IG1+IL)-ZORI(IG1+IL),4)
*         ENDIF
************************************************************************

      ENDDO

      DEALLOCATE(DEFDAT)

C ... An ugly fix for an IGRID error
      IF(.NOT.TIMEL) CALL GRIVAL(0,LDISTL,'REMESH')

      RETURN
      END SUBROUTINE REMESH
C
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C ------------------------------------------------------------------------
C
      SUBROUTINE DEFORM_VOX(IMAX,JMAX,KMAX,IN,JN,KN,XCO,YCO,ZCO,
     &   XORI,YORI,ZORI,XT,YT,ZT,DX,DY,DZ,
     &   XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,SX,SY,SZ,ISZ,JSZ,KSZ)

C ... Voxel based grid deformation.
      
C ... This subroutine deforms the volume grid according to the
C ... deformations in a rectangular transformation grid.

      USE NS3CO, ONLY : MPCCIL, MODALFSIL
      
      IMPLICIT NONE

      REAL, DIMENSION(*) :: XT, YT, ZT, DX, DY, DZ  
      REAL, DIMENSION(*) :: XCO, YCO, ZCO 
      REAL, DIMENSION(*) :: XORI, YORI, ZORI
      REAL, DIMENSION(3) :: XX
      REAL, DIMENSION(4) :: A0
      REAL :: SX, SY, SZ, RSX, RSY, RSZ
      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
      INTEGER, DIMENSION (8) :: II
      INTEGER :: NNODE, IK, ISZ, JSZ, KSZ, IS, JS, KS, KP, IT 
      INTEGER :: NPOINTS, IMAX, JMAX, KMAX, IN, JN, KN, I, J, K
      INTEGER :: IMAXP1, JMAXP1, KMAXP1, ISTRID, JSTRID, ILE
      LOGICAL :: PINKUB

      RSX = 1.0/SX
      RSY = 1.0/SY
      RSZ = 1.0/SZ

      NPOINTS = 0

      IMAXP1  = IMAX + 1
      JMAXP1  = JMAX + 1
      KMAXP1  = KMAX + 1
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      ILE     = ISTRID*JSTRID
      NNODE   = IMAXP1*JMAXP1*KMAXP1

*      DO K=0,KMAXP1+1
*         DO J=0,JMAXP1+1
*            DO I=0,IMAXP1+1
      DO K=0,KMAXP1
         DO J=0,JMAXP1
            DO I=0,IMAXP1

               IK = I+IN+(J-1+JN)*ISTRID+(K-1+KN)*ILE  ! grid node pointer

               XX(1) = XCO(IK)
               XX(2) = YCO(IK)
               XX(3) = ZCO(IK)

               IF(XX(1) < XMIN) CYCLE
               IF(XX(1) > XMAX) CYCLE
               IF(XX(2) < YMIN) CYCLE
               IF(XX(2) > YMAX) CYCLE
               IF(XX(3) < ZMIN) CYCLE
               IF(XX(3) > ZMAX) CYCLE

               IS = MIN0(INT((XX(1)-XMIN)*RSX) + 1,ISZ-1) 
               JS = MIN0(INT((XX(2)-YMIN)*RSY) + 1,JSZ-1)
               KS = MIN0(INT((XX(3)-ZMIN)*RSZ) + 1,KSZ-1)

               KP = (KS-1)*ISZ*JSZ + (JS-1)*ISZ + IS

               IT = 0

               IF(PINKUB(XX,KP,ISZ,JSZ,KSZ,XT,YT,ZT)) THEN

                  NPOINTS = NPOINTS + 1

                  CALL PVRINT(XX,A0,KP,IT,II,ISZ,JSZ,KSZ,XT,YT,ZT)

C ... Add displacements to the original grid.
                     
                    XCO(IK) = XORI(IK) + A0(1)*DX(II(1))+A0(2)*DX(II(2))
     &                                 + A0(3)*DX(II(3))+A0(4)*DX(II(4)) 
                    YCO(IK) = YORI(IK) + A0(1)*DY(II(1))+A0(2)*DY(II(2))
     &                                 + A0(3)*DY(II(3))+A0(4)*DY(II(4)) 
                    ZCO(IK) = ZORI(IK) + A0(1)*DZ(II(1))+A0(2)*DZ(II(2))
     &                                 + A0(3)*DZ(II(3))+A0(4)*DZ(II(4)) 
                  
C ... Add displacements to the deformed grid. May be needed one day.
c
c                    XCO(IK) = XCO(IK) + A0(1)*DX(II(1))+A0(2)*DX(II(2))
c     &                                + A0(3)*DX(II(3))+A0(4)*DX(II(4)) 
c                    YCO(IK) = YCO(IK) + A0(1)*DY(II(1))+A0(2)*DY(II(2))
c     &                                + A0(3)*DY(II(3))+A0(4)*DY(II(4)) 
c                    ZCO(IK) = ZCO(IK) + A0(1)*DZ(II(1))+A0(2)*DZ(II(2))
c     &                                + A0(3)*DZ(II(3))+A0(4)*DZ(II(4)) 


               ENDIF
               
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DEFORM_VOX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DEFORM_IDP(N,IMAX,JMAX,KMAX,IN,JN,KN,XCO,YCO,ZCO,
     &           XORI,YORI,ZORI,DISTP,IDP,DXM,DYM,DZM)

C ... Closest surface point based grid deformation.
      
C ... This subroutine deforms the volume grid according to the
C ... wetted surface displacements.

      USE MPCCIVARS, ONLY : fnodedisp
      USE INTEGERS,  ONLY : IPRO
      USE CONSTANTS, ONLY : PII
      USE NS3CO,     ONLY : CHLREF
      USE MAIN_ARRAYS, ONLY : NPROCE, NCHIMT
      
      IMPLICIT NONE

      REAL, DIMENSION(*) :: XCO, YCO, ZCO, XORI, YORI, ZORI, DISTP
      INTEGER :: N, NG, IMAX, JMAX, KMAX, IN, JN, KN, I, J, K, IK, IP
      INTEGER :: IMAXP1, JMAXP1, KMAXP1, ISTRID, JSTRID, ILE
      INTEGER :: IDP(*)
      REAL :: DX, DY, DZ, DXM, DYM, DZM, SCX, SCY, SCZ, FSTART 

      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      ILE    = ISTRID*JSTRID

      NG = NPROCE(1+N,IPRO)  ! Global block number

      DO K=0,KMAXP1
         DO J=0,JMAXP1
            DO I=0,IMAXP1

               IK = I+IN+(J-1+JN)*ISTRID+(K-1+KN)*ILE  ! Grid node
               IP = IDP(IK)                            ! Closest skin point

               DX = fnodedisp(3*(IP-1)+1)
               DY = fnodedisp(3*(IP-1)+2)
               DZ = fnodedisp(3*(IP-1)+3)

C ... Coefficients for fading out the displacements far from the  
C ... deforming surface (wetted surface).

*               FSTART = 0.5*CHLREF
*               
*               SCX = 1.0
*               SCY = 1.0
*               SCZ = 1.0
*               if(DISTP(IK) > FSTART) THEN
*               SCX = .5*(1.0+COS(MIN(DISTP(IK)/(2.0*FSTART)*PII,PII)))
*               SCY = .5*(1.0+COS(MIN(DISTP(IK)/(2.0*FSTART)*PII,PII)))
*               SCZ = .5*(1.0+COS(MIN(DISTP(IK)/(2.0*FSTART)*PII,PII)))
*               endif

C ... Linear fading

               FSTART = 0.3*CHLREF

               SCX = 1.0
               SCY = 1.0
               SCZ = 1.0
               if(DISTP(IK) > FSTART) THEN
                  SCX = MAX(2.0-DISTP(IK)/FSTART,0.0)
                  SCY = MAX(2.0-DISTP(IK)/FSTART,0.0)
                  SCZ = MAX(2.0-DISTP(IK)/FSTART,0.0)
               endif

C ... Add the wetted surface displacements to the volume grid.
               
                  XCO(IK) = XORI(IK) + DX*SCX
                  YCO(IK) = YORI(IK) + DY*SCY
                  ZCO(IK) = ZORI(IK) + DZ*SCZ

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DEFORM_IDP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DEFORM_SRF(XT,YT,ZT,DX,DY,DZ,XMIN,YMIN,ZMIN,
     &                      XMAX,YMAX,ZMAX,SX,SY,SZ,ISZ,JSZ,KSZ)

C ... Voxel based CFD surface grid deformation in case of modal FSI. 
C ... This subroutine is currently not used.
      
C ... This subroutine computes the surface displacements according to the
C ... deformations in a rectangular transformation grid.

      USE MPCCIVARS, ONLY : fnodecoor_orig, fnodedisp, fnnodes

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XT, YT, ZT, DX, DY, DZ  
      REAL, DIMENSION(3) :: XX
      REAL, DIMENSION(4) :: A0
      REAL :: SX, SY, SZ, RSX, RSY, RSZ
      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
      INTEGER, DIMENSION (8) :: II
      INTEGER :: IP, ISZ, JSZ, KSZ, IS, JS, KS, KP, IT 
      LOGICAL :: PINKUB

      RSX = 1.0/SX
      RSY = 1.0/SY
      RSZ = 1.0/SZ

      DO IP = 1,fnnodes

         XX(1) = fnodecoor_orig(IP,1)
         XX(2) = fnodecoor_orig(IP,2)
         XX(3) = fnodecoor_orig(IP,3)

         IF(XX(1) < XMIN) CYCLE
         IF(XX(1) > XMAX) CYCLE
         IF(XX(2) < YMIN) CYCLE
         IF(XX(2) > YMAX) CYCLE
         IF(XX(3) < ZMIN) CYCLE
         IF(XX(3) > ZMAX) CYCLE

         IS = MIN0(INT((XX(1)-XMIN)*RSX) + 1,ISZ-1) 
         JS = MIN0(INT((XX(2)-YMIN)*RSY) + 1,JSZ-1)
         KS = MIN0(INT((XX(3)-ZMIN)*RSZ) + 1,KSZ-1)

         KP = (KS-1)*ISZ*JSZ + (JS-1)*ISZ + IS

         IT = 0

         IF(PINKUB(XX,KP,ISZ,JSZ,KSZ,XT,YT,ZT)) THEN

            CALL PVRINT(XX,A0,KP,IT,II,ISZ,JSZ,KSZ,XT,YT,ZT)

C ... 3D interpolate the displacements.
                     
            fnodedisp(3*(IP-1)+1) = A0(1)*DX(II(1))+A0(2)*DX(II(2))
     &                            + A0(3)*DX(II(3))+A0(4)*DX(II(4)) 
            fnodedisp(3*(IP-1)+2) = A0(1)*DY(II(1))+A0(2)*DY(II(2))
     &                            + A0(3)*DY(II(3))+A0(4)*DY(II(4)) 
            fnodedisp(3*(IP-1)+3) = A0(1)*DZ(II(1))+A0(2)*DZ(II(2))
     &                            + A0(3)*DZ(II(3))+A0(4)*DZ(II(4)) 

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE DEFORM_SRF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PVRINT(X0,A0,K,IT,II,NX,NY,NZ,X,Y,Z)

      INTEGER :: K,IT,NX,NY,NZ
      REAL :: X0(3), A0(4)
      REAL, DIMENSION(NX*NY*NZ) :: X, Y, Z
      
C ... Calculate the weight factors for the tetrahedra corners.    
      
      REAL    :: VOL
      INTEGER :: II(8), K0, IT0
      
      K0    = 0
      IT0   = 0
      II(1) = 0     
      II(2) = 0     
      II(3) = 0     
      II(4) = 0     


      IF(K /= K0 .OR. IT /= IT0) THEN
         IF(IT == 0) THEN
            DO 10 IT=1,18
               CALL KTSNOD(K,IT,0,II,NX,NY)
      
               CALL TETAKO(X0,A0,VOL,II,NX*NY*NZ,X,Y,Z)
      
               IF(VOL > 0.0) THEN
                  IF(MIN(A0(1),A0(2),A0(3),A0(4)) >= 0.0) GO TO 15
               END IF
10          CONTINUE
               STOP'*** PVRINT IT=0'
15          CONTINUE
         ELSE
            CALL KTSNOD(K,IT,0,II,NX,NY)
         END IF
      
         K0=K
         IT0=IT
      END IF
      
      CALL TETAKO(X0,A0,VOL,II,NX*NY*NZ,X,Y,Z)
      
      IF(VOL == 0.0) STOP'Error in PVRINT: VOL = 0'
      
      RETURN
      END SUBROUTINE PVRINT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE KTSNOD(K,IT,IS,II,NX,NY)

C ... Node numbers for cell / tetrahedron / tetrahedron's face  K/IT/IS

      INTEGER :: K,IT,IS,NX,NY,II(8)
      
      INTEGER :: NX0,NY0
      INTEGER :: MK(8),JE(4,18),SS(3,4)
      SAVE MK
      
      DATA NX0,NY0/0,0/
      
*      DATA JE/1,2,4,5, 2,3,4,7, 2,4,5,7, 2,5,6,7, 4,5,7,8/  ! 5 tetrahedra

      DATA JE/1,2,5,3, 2,6,5,3, 3,6,5,7,
     &        3,8,7,5, 3,4,8,5, 1,4,3,5,
     &        5,6,1,7, 6,2,1,7, 7,2,1,3,
     &        7,4,3,1, 7,8,4,1, 5,8,7,1,
     &        8,5,4,6, 5,1,4,6, 6,1,4,2,
     &        6,3,2,4, 6,7,3,4, 8,7,6,4/
      
      DATA SS/1,3,2, 1,2,4, 2,3,4, 1,4,3/
      
      IF(NX /= NX0 .OR. NY /= NY0) THEN
         MK(1) = 0
         MK(2) = 1
         MK(3) = 1 + NX
         MK(4) =   + NX
         MK(5) =        + NX*NY
         MK(6) = 1      + NX*NY
         MK(7) = 1 + NX + NX*NY
         MK(8) =   + NX + NX*NY
         NX0   = NX
         NY0   = NY
      END IF
      
      IF(IS /= 0) THEN
         II(1) = K + MK(JE(SS(1,IS),IT))
         II(2) = K + MK(JE(SS(2,IS),IT))
         II(3) = K + MK(JE(SS(3,IS),IT))
      ELSE IF(IT /= 0) THEN
         II(1) = K + MK(JE(1,IT))
         II(2) = K + MK(JE(2,IT))
         II(3) = K + MK(JE(3,IT))
         II(4) = K + MK(JE(4,IT))
      ELSE IF(K /= 0) THEN
         II(1) = K + MK(1)
         II(2) = K + MK(2)
         II(3) = K + MK(3)
         II(4) = K + MK(4)
         II(5) = K + MK(5)
         II(6) = K + MK(6)
         II(7) = K + MK(7)
         II(8) = K + MK(8)
      ELSE
         STOP'*** KTSNOD'
      END IF
      
      RETURN
      END SUBROUTINE KTSNOD 
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE TETAKO(X0,A0,VOL,II,N,X,Y,Z)
      
C ... Compute the volume coordinates.

      INTEGER            :: N, II(4)

      REAL               :: X0(3), A0(4), VOL
      REAL, DIMENSION(N) :: X, Y, Z

      REAL, DIMENSION(3) :: A, B, C, D
      REAL               :: V0,A01,A02,A03,A04,DETD3
      
      CALL V3DIFF(II(2),II(1),A,N,X,Y,Z)
      CALL V3DIFF(II(3),II(1),B,N,X,Y,Z)
      CALL V3DIFF(II(4),II(1),C,N,X,Y,Z)
      
      V0 = DETD3(A,B,C)
      
      IF(ABS(V0) <= 1.0E-30) THEN
         VOL = 0.0
         RETURN
      ENDIF
      
      D(1)  = X0(1)-X(II(1))
      D(2)  = X0(2)-Y(II(1))
      D(3)  = X0(3)-Z(II(1))
      
      A02   = DETD3(B,C,D)/V0
      A03   = DETD3(C,A,D)/V0
      A04   = DETD3(A,B,D)/V0
      
      A01   = 1.0-A02-A03-A04
      
      VOL   = ABS(V0/6.0)
      
      A0(1) = A01
      A0(2) = A02
      A0(3) = A03
      A0(4) = A04
      
      RETURN
      END SUBROUTINE TETAKO
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE V3DIFF(I1,I2,A,N,X,Y,Z)

      IMPLICIT NONE

      INTEGER                      :: I1,I2,N
      REAL, DIMENSION(N) :: X, Y, Z
      REAL, DIMENSION(3) :: A
      
      A(1) = X(I1) - X(I2)
      A(2) = Y(I1) - Y(I2)
      A(3) = Z(I1) - Z(I2)
      
      RETURN
      END SUBROUTINE V3DIFF
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      REAL FUNCTION DETD3(A,B,C)
      
C                   |  A  |
C ... Determinant   |  B  |
C                   |  C  |

      REAL, DIMENSION(3) :: A, B, C
      
      DETD3 = C(1)*(A(2)*B(3)-B(2)*A(3)) +
     &        C(2)*(A(3)*B(1)-B(3)*A(1)) +
     &        C(3)*(A(1)*B(2)-B(1)*A(2))
      
      RETURN
      END FUNCTION DETD3

C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     

      LOGICAL FUNCTION PINKUB(X0,K,NX,NY,NZ,X,Y,Z)
      
C ... If X0 lies inside cell K, then PINKUB = .TRUE.

      IMPLICIT NONE

      INTEGER                   :: K,NX,NY,NZ, IT
      INTEGER, DIMENSION(8)     :: II
      REAL :: X0(3), A0(4)
      REAL, DIMENSION(NX*NY*NZ) :: X, Y, Z
      REAL :: VOL
      
      CALL KTSNOD(K,0,0,II,NX,NY)

      
      PINKUB = .FALSE.
      
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
      
      PINKUB = .TRUE.

      DO IT=1,18

         CALL KTSNOD(K,IT,0,II,NX,NY)
      
         CALL TETAKO(X0,A0,VOL,II,NX*NY*NZ,X,Y,Z)
      
         IF(VOL > 0.0) THEN
            IF(MINVAL(A0) >= 0.0) RETURN
         END IF

      ENDDO
      
      PINKUB = .FALSE.
      
      RETURN
      END FUNCTION PINKUB
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      
