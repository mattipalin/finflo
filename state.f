C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE STATE(VIS,P,PDIFF,RO,E,TEMP,CP,CH,DRDP,DRDH,GAMMA,
     &                 PR,PRT,NTOT,ITURB,ISTATE,FRSPRE,RGAS,VISU0,
     &                 EXPSU,TSU0,E0REF,T0REF,IMAX,JMAX,KMAX,INN,RMULT,
     &                 FRSDEN,FRSVIS,FRSTEM,FRSSIE,TRGS,IUPPTL)

! ... Calling Routine for the Equation of State
! ... INPUT  : P,T
! ... OUTPUT : P,TEMP,DRDP,DRDH,CH,VIS or RO,TEMP,CH,VIS

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER           :: ERCODE,INN,ISTATE,NTOT,IMAX,JMAX,KMAX,ITURB,
     &                     IUPPTL

      REAL              :: FRSVIS,RGAS,GAMMA,PR,PRT,VISU0,EXPSU,TSU0,
     &                     E0REF,T0REF,RMULT,FRSTEM,FRSSIE,FRSPRE,FRSDEN

      REAL,DIMENSION(*) :: P,RO,E,VIS,TEMP,CP,CH,DRDP,DRDH,PDIFF,TRGS

      REAL,DIMENSION(1) :: H,C,S ! Dummy arrays

      ERCODE = 0

      SELECT CASE(ISTATE)

! **********************************************************************
! ... Perfect Gas
! **********************************************************************
      CASE(1)
       CALL PERFGA(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,VISU0,
     &             EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF,RMULT,TRGS,FRSPRE,
     &             FRSTEM,IUPPTL)

! **********************************************************************
! ... Chemically Reacting Equilibrium Air (Arrays are not defined)!!
! **********************************************************************
      CASE(2)
       CALL TGAS(P,RO,E,H,TEMP,C,S,DRDP,DRDH, 2,NTOT,ERCODE)
       CALL TGAS(P,RO,E,H,TEMP,C,S,DRDP,DRDH,11,NTOT,ERCODE)
       CALL UGAS3(E,RO,VIS,NTOT,ERCODE)
       CALL UGAS4(E,RO,CH,NTOT,ERCODE)
       PRINT"(A/A)","This equation of state do not work right now.",
     &              "Should input: temp, p, but now ro and e."

       DRDH(1:NTOT) = 1./(DRDH(1:NTOT) - P(1:NTOT)/RO(1:NTOT)**2)
       DRDP(1:NTOT) =   -(DRDP(1:NTOT) + 1.       /RO(1:NTOT)   )
     &                  * DRDH(1:NTOT)
       STOP

! **********************************************************************
! ... Perfect Gas with cp(T)
! **********************************************************************
      CASE(3)
       CALL PERFGA2(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
     &              VISU0,EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &              TRGS,FRSPRE,FRSTEM,IUPPTL)

! **********************************************************************
! ... Semi perfect Gas
! **********************************************************************
      CASE(4)
       CALL PERSEM(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
     &             VISU0,EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &             TRGS,FRSPRE,FRSTEM,IUPPTL)

! **********************************************************************
! ... Air modeled as a perfect gas (GAMMA and CP are not constants)
! **********************************************************************
      CASE(5)
       CALL PERAIR(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
     &             VISU0,EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &             TRGS,FRSPRE,FRSTEM,IUPPTL)

! **********************************************************************
! ... Water (old version)
! **********************************************************************
      CASE(6)
       CALL WATER(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,NTOT,
     &            IMAX,JMAX,KMAX,IN,JN,KN,INN,FRSTEM,IUPPTL,
     &             FRSDEN,FRSSIE,FRSVIS,PR)

! **********************************************************************
! ... Water with modified viscosity (modified Reynolds number)
! **********************************************************************
      CASE(7)
       CALL PWATER(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN,RMULT,FRSTEM,IUPPTL)

! **********************************************************************
! ... Water as a function of pressure and temperature
! **********************************************************************
      CASE(8)
       CALL WATER2(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN)

! **********************************************************************
! ... Steam as function of pressure and temperature
! **********************************************************************
      CASE(9)
       CALL STEAM3(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN)

! **********************************************************************
! ... Fully incompressible flow based on the free-stream properties
! **********************************************************************
      CASE(10:11)
       CALL FRSVAL(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &   FRSDEN,FRSVIS,FRSTEM,FRSSIE,RMULT)
! **********************************************************************
! ... Liquid R12 as a function of pressure and temperature
! **********************************************************************
      CASE(12)
       CALL R12LIQ(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN)
! **********************************************************************
! ... Gaseous R12 as a function of pressure and temperature
! **********************************************************************
      CASE(13)
       CALL R12GAS(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Routines for solids
! **********************************************************************
! ... Cast iron
! **********************************************************************
      CASE(101)
       CALL CAST_IRON(RO,E,TEMP,CP,CH,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Routine for constant properties in wall-distance equation
! **********************************************************************
! ... Outputs = 1
! **********************************************************************
      CASE(201)
       CALL DISTPRO(RO,E,TEMP,CP,CH,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN)

      END SELECT

      IF (ERCODE > 0) WRITE(4,"(/9X,I10,A)") ERCODE,
     &                "OVERFLOW ERRORS occured in STATE subroutine!"

      END SUBROUTINE STATE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

! **********************************************************************
! ... STATE EQUATIONS
! **********************************************************************

      SUBROUTINE PERFGA(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
     &  VIS0,EXPSU,T0,NTOT,ERCODE,E0REF,T0REF,RMULT,TRGS,FRSPRE,FRSTEM,
     &  IUPPTL)

! ... Equation of State for a Perfect Gas, 17.1.1994, PMK
! ... Sutherland's Formula is used for viscosity
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS,CP

      IMPLICIT NONE

      INTEGER :: NTOT,I,ERCODE,IUPPTL
      REAL,DIMENSION(*) :: RO,E,P,TEMP,DRDP,DRDH,CP,CH,VIS,TRGS
      REAL :: RGAS,GAMMA,VIS0,EXPSU,T0,PR,E0REF,T0REF,AVEG,
     &        CP2,GM1,GM2,RMULT,FRSPRE,FRSTEM

      GM1 = GAMMA - 1.
      GM2 = GAMMA*GM1
      CP2 = GAMMA/GM1*RGAS

      IF (ABS(EXPSU-1.5) > 1e-5) THEN
       TRGS(1:NTOT) =    TEMP(1:NTOT)*RGAS
       DRDP(1:NTOT) = 1./TRGS(1:NTOT)

       SELECT CASE(IUPPTL)
       CASE(0) 
          RO(1:NTOT)  = P(1:NTOT)*DRDP(1:NTOT)  ! Compressible
          VIS(1:NTOT) = VIS0*TEMP(1:NTOT)**EXPSU/(TEMP(1:NTOT)+T0)*RMULT
       CASE(1)
          RO(1:NTOT)  = FRSPRE/(RGAS*FRSTEM)    ! Incompressible
          VIS(1:NTOT) = VIS0*FRSTEM**EXPSU/(FRSTEM+T0)*RMULT
       CASE(2)
          RO(1:NTOT)  = FRSPRE*DRDP(1:NTOT)     ! Constant pressure
          VIS(1:NTOT) = VIS0*TEMP(1:NTOT)**EXPSU/(TEMP(1:NTOT)+T0)*RMULT
       CASE(3)
          RO(1:NTOT)  = P(1:NTOT)/(RGAS*FRSTEM) ! Constant temperature
          VIS(1:NTOT) = VIS0*FRSTEM**EXPSU/(FRSTEM+T0)*RMULT
       END SELECT

       E(1:NTOT)    =    TRGS(1:NTOT)/GM1
       DRDH(1:NTOT) =     -RO(1:NTOT)/E(1:NTOT)/GAMMA
       CP(1:NTOT)   = CP2
       CH(1:NTOT)     = VIS(1:NTOT)*CP2/PR

      ELSE ! SQRT is possible because of EXPSU = 1.5

       TRGS(1:NTOT) =    TEMP(1:NTOT)*RGAS
       DRDP(1:NTOT) = 1./TRGS(1:NTOT)

       SELECT CASE(IUPPTL)
       CASE(0) 
          RO(1:NTOT)  = P(1:NTOT)*DRDP(1:NTOT)  ! Compressible
       VIS(1:NTOT)  = VIS0*TEMP(1:NTOT)*RMULT
     &               *SQRT(TEMP(1:NTOT))/(TEMP(1:NTOT)+T0)
       CASE(1)
          RO(1:NTOT)  = FRSPRE/(RGAS*FRSTEM)    ! Incompressible
          VIS(1:NTOT) = VIS0*FRSTEM*RMULT*SQRT(FRSTEM)/(FRSTEM+T0)
       CASE(2)
          RO(1:NTOT)  = FRSPRE*DRDP(1:NTOT)     ! Constant pressure
       VIS(1:NTOT)  = VIS0*TEMP(1:NTOT)*RMULT
     &               *SQRT(TEMP(1:NTOT))/(TEMP(1:NTOT)+T0)
       CASE(3)
          RO(1:NTOT)  = P(1:NTOT)/(RGAS*FRSTEM) ! Constant temperature
          VIS(1:NTOT) = VIS0*FRSTEM*RMULT*SQRT(FRSTEM)/(FRSTEM+T0)
       END SELECT

       E(1:NTOT)      = TRGS(1:NTOT)/GM1
       DRDH(1:NTOT)   =-RO(1:NTOT)/E(1:NTOT)/GAMMA
       CP(1:NTOT)   = CP2
       CH(1:NTOT)   = VIS(1:NTOT)*CP2/PR

      END IF

      END SUBROUTINE PERFGA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERFGA2(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,
     &                   PR,VIS0,EXPSU,T0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &                   TRGS,FRSPRE,FRSTEM,IUPPTL)

! ... Equation of State for a Perfect Gas, 17.1.1994, PMK
! ... Sutherland's Formula is used for viscosity
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS,CP

      IMPLICIT NONE

      INTEGER :: NTOT,I,ERCODE,IUPPTL
      REAL,DIMENSION(*) :: RO,E,P,TEMP,DRDP,DRDH,CP,CH,VIS,TRGS
      REAL :: RGAS,GAMMA,VIS0,EXPSU,T0,PR,E0REF,T0REF,AVEG,
     &        CP0,CP1,CP2,CP3,GM1,GM2,TRGAS,RMULT,FRSPRE,FRSTEM

      GM1  = GAMMA - 1.
      GM2  = GAMMA*GM1
      CP2  = GAMMA/GM1*RGAS

      IF (ABS(EXPSU-1.5) > 1e-5) THEN
       DO I = 1,NTOT
        TRGAS   = TEMP(I)*RGAS
        DRDP(I) = 1./TRGAS

       SELECT CASE(IUPPTL)
       CASE(0) 
          RO(I)  = P(I)*DRDP(I)  ! Compressible
          VIS(I) = VIS0*TEMP(I)**EXPSU/(TEMP(I)+T0)*RMULT
       CASE(1)
          RO(I)  = FRSPRE/(RGAS*FRSTEM)    ! Incompressible
          VIS(I) = VIS0*FRSTEM**EXPSU/(FRSTEM+T0)*RMULT
       CASE(2)
          RO(I)  = FRSPRE*DRDP(I)     ! Constant pressure
          VIS(I) = VIS0*TEMP(I)**EXPSU/(TEMP(I)+T0)*RMULT
       CASE(3)
          RO(I)  = P(I)/(RGAS*FRSTEM) ! Constant temperature
          VIS(I) = VIS0*FRSTEM**EXPSU/(FRSTEM+T0)*RMULT
       END SELECT

        CP0     =    .124500e+4
        CP1     = ((-.108980e-6 *TEMP(I)+.333093e-3)*TEMP(I)
     &              -.118565e+0)*TEMP(I)+.101582e+4
        CP3     =    .296940e+0 *TEMP(I)+.703060e+3
        IF (CP0 < MAX(CP1,CP3)) THEN
         E(I)   =     .124500e+4 *TEMP(I)-TRGAS
        ELSE IF (CP1 > CP3) THEN
         E(I)   = (((-.272450e-7 *TEMP(I)+.111031e-3)*TEMP(I)
     &               -.592825e-1)*TEMP(I)+.101582e+4)*TEMP(I)-TRGAS
        ELSE
         E(I)   =   ( .148470e+0 *TEMP(I)+.703060e+3)*TEMP(I)-TRGAS
        END IF

        DRDH(I) = -RO(I)*DRDP(I)*GM1/GAMMA
        CP(I)   = MIN(MAX(CP1,CP3),CP0)
        CH(I)   = VIS(I)*CP(I)/PR
       END DO

      ELSE ! SQRT is possible because of EXPSU = 1.5

       DO I = 1,NTOT
        TRGAS   = TEMP(I)*RGAS
        DRDP(I) = 1./TRGAS

       SELECT CASE(IUPPTL)
       CASE(0) 
          RO(I)  = P(I)*DRDP(I)  ! Compressible
          VIS(I) = VIS0*TEMP(I)*RMULT
     &                * SQRT(TEMP(I))/(TEMP(I)+T0)
       CASE(1)
          RO(I)  = FRSPRE/(RGAS*FRSTEM)    ! Incompressible
          VIS(I) = VIS0*FRSTEM*RMULT*SQRT(FRSTEM)/(FRSTEM+T0)
       CASE(2)
          RO(I)  = FRSPRE*DRDP(I)     ! Constant pressure
          VIS(I) = VIS0*TEMP(I)*RMULT
     &                * SQRT(TEMP(I))/(TEMP(I)+T0)
       CASE(3)
          RO(I)  = P(I)/(RGAS*FRSTEM) ! Constant temperature
          VIS(I) = VIS0*FRSTEM*RMULT*SQRT(FRSTEM)/(FRSTEM+T0)
       END SELECT

        CP0     =    .124500e+4
        CP1     = ((-.108980e-6 *TEMP(I)+.333093e-3)*TEMP(I)
     &              -.118565e+0)*TEMP(I)+.101582e+4
        CP3     =    .296940e+0 *TEMP(I)+.703060e+3
        IF (CP0 < MAX(CP1,CP3)) THEN
         E(I)   =     .124500e+4 *TEMP(I)-TRGAS
        ELSE IF (CP1 > CP3) THEN
         E(I)   = (((-.272450e-7 *TEMP(I)+.111031e-3)*TEMP(I)
     &               -.592825e-1)*TEMP(I)+.101582e+4)*TEMP(I)-TRGAS
        ELSE
         E(I)   =   ( .148470e+0 *TEMP(I)+.703060e+3)*TEMP(I)-TRGAS
        END IF

        DRDH(I) = -RO(I)*DRDP(I)*GM1/GAMMA
        CP(I)   = MIN(MAX(CP1,CP3),CP0)
        CH(I)   = VIS(I)*CP2/PR
       END DO
      END IF

      END SUBROUTINE PERFGA2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERAIR(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,
     &                  PR,VIS0,EXPSU,T0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &                  TRGS,FRSPRE,FRSTEM,IUPPTL)


! ... Equation of State for a Perfect Gas, 17.1.1994, PMK
! ... Sutherland's Formula for viscosity
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS,CP and average GAMMA

      IMPLICIT NONE

      INTEGER :: NTOT,I,ERCODE,IUPPTL
      REAL,DIMENSION(*) :: RO,E,P,TEMP,DRDP,DRDH,CP,CH,VIS,TRGS
      REAL :: RGAS,GAMMA,VIS0,EXPSU,T0,PR,E0REF,T0REF,AVEG,
     &        CP2,GM1,GM2,TRGAS,RMULT,FRSPRE,FRSTEM

      AVEG = 0.

      DO I = 1,NTOT
       TRGAS   = TEMP(I)*RGAS
       CP2     = GAMMA*RGAS/(GAMMA-1.)
       CP(I)   = CP2
       GAMMA   = CP2/(CP2-RGAS)
       GM1     = GAMMA - 1.
       GM2     = GAMMA*GM1
       DRDP(I) = 1./TRGAS
       RO(I)   = P(I)*DRDP(I)
       E(I)    = TRGAS/GM1
       DRDH(I) = -RO(I)/E(I)/GAMMA
*       VIS(I)  = VIS0*TEMP(I)**1.5/(TEMP(I)+T0)*RMULT
       VIS(I)  = VIS0*TEMP(I)*SQRT(TEMP(I))/(TEMP(I)+T0)*RMULT
       CH(I)   = VIS(I)*CP2/PR
       AVEG    = GAMMA + AVEG
      END DO

      GAMMA = AVEG/NTOT

      END SUBROUTINE PERAIR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WATER(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,NTOT,
     &             IMAX,JMAX,KMAX,IN,JN,KN,INN,FRSTEM,IUPPTL,
     &             FRSDEN,FRSSIE,FRSVIS,PR)

! ... Equation of State for water; 8.1.96 KPe
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS

      IMPLICIT NONE

      REAL,DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: A,B,C,D,CP2,CV,T,TVIS,PR,FRSTEM

      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IL,IUPPTL
      REAL :: G11, G21, G12, G22,FRSDEN,FRSSIE,FRSVIS		! Viskositeetin kaava -MPa

!      A    =  .291648393e+1
!      B    =  .259588065e-2
!      C    =  .958787090e+2
      A    =  .2414e-4          ! ESa 24.8.2001
      B    =  .2478e+3
      C    =  .1400e+3
      D    =  .113902e-5
      CP2  =  .4180e+4
      CV   =  .4180e+4          ! EI NYT IHAN OIKEA ARVO
      T    = -.4850e+1
      PR   =  .8000e+1

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

!      DO I = 1,NTOT

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

         RO(I)   = FRSDEN
         E(I)    =  .8377e+5 + (TEMP(I)-.2930e+3)*CV
         DRDP(I) =  .5000e-6
         DRDH(I) = -.4970e-4
         TVIS    = MAX(.27315e+3,MIN(.37000e+3,TEMP(I))) ! ?
         TVIS    = TEMP(I)
         IF(IUPPTL == 1 .OR. IUPPTL == 3) TVIS = FRSTEM

	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta TVIS
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

c	G11 = .2414e-4 * 10.**(.2478e+3/(TVIS-.1400e+3))
c	G21 = MAX(.1e-3,(-.363718e-7*TVIS+.503834e-4)*TVIS
c     &  -.233632e-1        +.369516e+1/TVIS)
c	G12 = MAX(0.,MIN(1.0,0.1 * (310. - TVIS)))
c	G22 = 1.0 - G12

c	VIS(I) =  G12 * G11 + G22 * G21

	! ---------------------------------
         CP(I)   = CP2
!         CH(I)   = VIS(I)*CP(I)/PR
!         CH(I)   =  .556e+0 + .125e-2*(TVIS-.275e+3)
       VIS(I)  = FRSVIS               ! *RMULT 
       CH(I)   = VIS(I)*CP(I)/PR

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE WATER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PWATER(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,NTOT,
     &              IMAX,JMAX,KMAX,IN,JN,KN,INN,RMULT,FRSTEM,IUPPTL)

! ... Equation of State for water; 8.1.96 KPe
! ... This version can be used to modify the Reynolds number using viscosity
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS

      IMPLICIT NONE

      REAL,DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: A,B,C,D,CP2,CV,T,TVIS,RMULT,PR,FRSTEM
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IL,IUPPTL
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

      A    =  .2414e-4          ! ESa 24.8.2001
      B    =  .2478e+3
      C    =  .1400e+3
      D    =  .113902e-5
      CP2  =  .4180e+4
      CV   =  .4180e+4          ! EI NYT IHAN OIKEA ARVO
      T    = -.4850e+1
      PR   =  .8000e+1

!      RMULT= 31.558*SQRT(31.558)!
!       onkalo1250. !5640.7 ! 31.558*SQRT(31.558) ! LNG
!      RMULT = 2000. ! Oli 10k

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

         RO(I)   =  .1000e+4
         E(I)    =  .8377e+5 + (TEMP(I)-.2930e+3)*CV
         DRDP(I) =  .5000e-6
         DRDH(I) = -.4970e-4
         TVIS    = MAX(.27315e+3,MIN(.37000e+3,TEMP(I)))
         TVIS    = TEMP(I)
         IF(IUPPTL == 1 .OR. IUPPTL == 3) TVIS = FRSTEM
	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta TVIS
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(TVIS-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*TVIS+.503834e-4)*TVIS
     &  -.233632e-1        +.369516e+1/TVIS)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - TVIS)))
	G22 = 1 - G12

	VIS(I) =  G12 * G11 + G22 * G21

	! ---------------------------------
         CP(I)   = CP2
         CH(I)   =  .556e+0 + .125e-2*(TVIS-.275e+3)! TSii 7.4.2016

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE PWATER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERSEM(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,
     &                  PR,VIS0,EXPSU,T0,NTOT,ERCODE,E0REF,T0REF,RMULT,
     &                  TRGS,FRSPRE,FRSTEM,IUPPTL)

! ... Equation of State for a Semi-Perfect Gas, 12.5.1995, PPR
! ... INPUT  : RO,E
! ... OUTPUT : RO,TEMP,DRDP,DRDH,CH,VIS

      REAL,DIMENSION(*) :: RO,E,P,TEMP,DRDP,DRDH,CH,CP,VIS,TRGS
      REAL :: RGAS,GAMMA,VIS0,EXPSU,T0,FRSPRE,FRSTEM
      INTEGER :: NTOT,I,ERCODE,IUPPTL

      GM1 = GAMMA - 1.
      GM2 = GAMMA*GM1
      CP2 = GAMMA/GM1*RGAS

!      DO I = 1,NTOT
! ... viritysta
!       if (e(i)+e0ref <= 0.) write(*,*) i,e(i),ro(i),e(i)+e0ref
!       Pkoe    = (GAMMA-1.)/100.*(E(I))*RO(I)
!       Poik    = (GAMMA-1.)*(E(I)+E0REF)*RO(I)
!       if (pkoe >= poik) write(*,*) i,pkoe,poik,pkoe-poik,e(i)
!       P(i)    = MAX(PKOE,POIK)
! ... ylla olevalla voidaan estaa paineen meneminen negatiiviseksi
!      END DO

      P(1:NTOT)    = GM1*(E(1:NTOT)+E0REF)*RO(1:NTOT)
      TEMP(1:NTOT) = MAX(T0REF + P(1:NTOT)/(RO(1:NTOT)*RGAS),.1e+3)
      DRDP(1:NTOT) = 1./(GM1*E(1:NTOT))
      DRDH(1:NTOT) = -P(1:NTOT)/(E(1:NTOT)**2 * GM2)

!      CP(1:NTOT)   = GAMMA*RGAS/GM1
! ... Sutherland's Formula
      VIS(1:NTOT)  = VIS0*TEMP(1:NTOT)**EXPSU/(TEMP(1:NTOT)+T0)*RMULT
      CP(1:NTOT)   = CP2
      CH(1:NTOT)   = VIS(1:NTOT)*CP(1:NTOT)/PR

      END SUBROUTINE PERSEM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FRSVAL(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,FRSDEN,
     &           FRSVIS,FRSTEM,FRSSIE,RMULT)

! ... Equation of State for truly incompressible flows
! ... INPUT  : FRSDEN,FRSTEM,FRSSIE,FRSVIS
! ... OUTPUT : RO,TEMP,DRDP,DRDH,CH,VIS

      REAL,DIMENSION(*) :: RO,E,P,TEMP,DRDP,DRDH,CH,CP,VIS
      REAL :: FRSDEN,FRSSIE,FRSTEM,FRSVIS,RMULT
      INTEGER :: NTOT,I,ERCODE

      DO I = 1,NTOT
       DRDP(I) = 1.E-6    !  1./225. !1.E-6
       RO(I)   = FRSDEN
       E(I)    = FRSSIE
       CP(I)   = FRSSIE/FRSTEM
       DRDH(I) = 1.E-7
       VIS(I)  = FRSVIS               ! *RMULT 
       CH(I)   = VIS(I)*CP(I)/PR
      END DO

      END SUBROUTINE FRSVAL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EFPRO(E,P,RO,NTOT,ISTATE,GAMMA,FRSPRE,E0REF,T0REF)

! ... ENERGY AS A FUNCTION OF PRESSURE AND DENSITY (PERFECT GAS)

      IMPLICIT NONE

      REAL,DIMENSION(*) :: E,P,RO
      REAL :: GAMMA,GM1,FRSPRE,E0REF,T0REF,DMMY
      INTEGER :: NTOT,ISTATE,ERCODE

      ERCODE = 0

      SELECT CASE(ISTATE)

      CASE(8,9,11,12)
       PRINT"(//A/A//)",
     &      "ISTATE >= 8 does not have EFPRO option.",
     &      "Exiting..." ; STOP

      CASE(1,3,5) ! PERFECT GAS
       GM1 = GAMMA - 1.
       E(1:NTOT) = (P(1:NTOT)+FRSPRE)/(GM1*RO(1:NTOT))-E0REF

      CASE(4)     ! SEMI-PERFECT GAS, PPR 12.5.1995
       GM1 = GAMMA - 1. ; E(1:NTOT) = P(1:NTOT)/(GM1*RO(1:NTOT))-E0REF

      CASE(2)     ! Chemically Reacting Equilibrium Air
       CALL TGAS(P,RO,E,DMMY,DMMY,DMMY,DMMY,DMMY,DMMY,10,NTOT,ERCODE)
!       PRINT"(//A//)","No such model of E.o.S. Exiting..." ; STOP

      CASE(6:7) ! WATER (8.1.1996 KPe) & PWATER (13.3.2002 TSii)
       E(1:NTOT) = .8377e+5

      CASE(10) ! Fully incompressible

        E(1:NTOT) = 1.E 6 ! Non sense

c       E(1:NTOT) = FRSSIE

      END SELECT

      END SUBROUTINE EFPRO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROFPT(TEMPG,PG,RO,NTOT,ISTATE,RGAS,
     &           GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)

! ... DENSITY AS A FUNCTION OF PRESSURE AND TEMPERATURE

      IMPLICIT NONE

      INTEGER :: ERCODE,NTOT,ISTATE,I,IUPPTL
      REAL,DIMENSION(*) :: TEMPG,PG,RO
      REAL :: GAMMA,RGAS,FRSPRE,E0REF,T0REF,GM1,FRSDEN,SLP1,SLP2,F21,
     &        F00,F22,R2,ROO1,DRDT1,R2P,F31,F32,R3,ROO10,DRDT10,R3P,
     &        F41,F42,R4,DRDT100,R4P,F51,F52,R5,ROO300,DRDT300,R5P,
     &        F61,F62,R6,ROO2M,DRDT2M,R6P,ROO100,LORE,FRSTEM

      REAL,ALLOCATABLE,DIMENSION(:) :: TEMP, P

      ALLOCATE(TEMP(NTOT),P(NTOT))
         
C ... Cases

      SELECT CASE(IUPPTL)
         CASE(0)
         TEMP(1:NTOT) = TEMPG(1:NTOT) ; P(1:NTOT) = PG(1:NTOT)
         CASE(1)
         TEMP(1:NTOT) = FRSTEM ;        P(1:NTOT) = FRSPRE
         CASE(2)
         TEMP(1:NTOT) = TEMPG(1:NTOT) ; P(1:NTOT) = FRSPRE
         CASE(3)
         TEMP(1:NTOT) = FRSTEM ;        P(1:NTOT) = PG(1:NTOT)
      END SELECT

      ERCODE = 0

      SELECT CASE(ISTATE)

      CASE(1,3:4)                     !  PERFECT GAS & SEMI-PERFECT GAS PPR 12.5.1995

       GM1 = GAMMA - 1. ; RO(1:NTOT) = P(1:NTOT)/(RGAS*TEMP(1:NTOT))

      CASE(2)                         !  Chemically Reacting Equilibrium Air
       PRINT"(//A//)","No such model of E.o.S. Exiting..." ; STOP

      CASE(5)
       PRINT"(//A/A//)","ISTATE = 5 is not in use anymore.",
     &                  "Use ISTATE = 1. Exiting..." ; STOP

      CASE(6:7)                       !  WATER & PWATER joo joo
       RO(1:NTOT) = FRSDEN !.1000e+4

      CASE(8)                         !  Water from 0,5 to 1,5 bar
       RO(1:NTOT) = (-.260031e-2*TEMP(1:NTOT)+.122284e+1)*TEMP(1:NTOT)
     &            + ( .761625e-9*TEMP(1:NTOT)+.196854e-6)*P(1:NTOT)
     &            + ( .864232e+3)

      CASE(9)                         !  Steam from 0,01 to 20 bar

! ... Start of new addition by A. Miettinen, 2007/10/01

      DO I = 1,NTOT
! ... 1kPa
      SLP1 = 0.001
      SLP2 = 1.0/9.0e+3
      F21  =            SLP1*P(I)
      F00  = 0.111111 - SLP2*P(I)
      F22  = 1.0 + F00
      R2   = MAX(0.,MIN(F21,F22))
!      ROO1    =  -.134939e-4*TEMP(I)+.111220e-1
      ROO1    = ( .322363e-7*TEMP(I)-.406131e-4)*TEMP(I)+.165018e-1
      DRDT1   =   .644726e-7*TEMP(I)-.406131e-4
      R2P     = LORE(R2 > 0.0)*(SLP1*LORE(F21 < F22)
     &                        - SLP2*LORE(F21 > F22))

! ... 10kPa
      SLP1 = SLP2
      SLP2 = 1.0/9.0e+4
      F31  = -F00
      F00  = 0.111111 - SLP2*P(I)
      F32  = 1.0 + F00
      R3   = MAX(0.,MIN(F31,F32))
!      ROO10   =  -.119482e-3*TEMP(I)+.103647e+0
      ROO10   = ( .281416e-6*TEMP(I)-.367450e-3)*TEMP(I)+.156294e+0
      DRDT10  =   .562832e-6*TEMP(I)-.367450e-3
      R3P     = LORE(R3 > 0.0)*(SLP1*LORE(F31 < F32)
     &                        - SLP2*LORE(F31 > F32))

! ... 100kPa
      SLP1 = SLP2
      SLP2 = 1.0/2.0e+5
      F41  = -F00
      F00  = 0.500000 - SLP2*P(I)
      F42  = 1.0 + F00
      R4   = MAX(0.,MIN(F41,F42))
!      ROO100  =  -.105576e-2*TEMP(I)+.972238e+0
      ROO100  = ( .240937e-5*TEMP(I)-.331947e-2)*TEMP(I)+.149104e+1
      DRDT100 =   .481874e-5*TEMP(I)-.331947e-2
      R4P     = LORE(R4 > 0.0)*(SLP1*LORE(F41 < F42)
     &                        - SLP2*LORE(F41 > F42))

! ... 300kPa
      SLP1 = SLP2
      SLP2 = 1.0/1.7e+6
      F51  = -F00
      F00  = 0.176471 - SLP2*P(I)
      F52  = 1.0 + F00
      R5   = MAX(0.,MIN(F51,F52))
!      ROO300  =  -.304522e-2*TEMP(I)+.286537e+1
      ROO300  = ( .726955e-5*TEMP(I)-.101508e-1)*TEMP(I)+.457358e+1
      DRDT300 =   .145391e-4*TEMP(I)-.101508e-1
      R5P     = LORE(R5 > 0.0)*(SLP1*LORE(F51 < F52)
     &                        - SLP2*LORE(F51 > F52))

! ... 2Mpa
      SLP1 = SLP2
      SLP2 = 0.0
      F61  = -F00
      F62  = 1.0
      R6   = MAX(0.,MIN(F61,F62))
!      ROO2M   =  -.151110e-1*TEMP(I)+.169803e+2
      ROO2M   = ( .348933e-4*TEMP(I)-.585389e-1)*TEMP(I)+.301263e+2
      DRDT2M  =   .697866e-4*TEMP(I)-.585389e-1
      R6P     = LORE(R6 > 0.0)*(SLP1*LORE(F61 < F62)
     &                        - SLP2*LORE(F61 > F62))

      RO(I) = R2 *ROO1 +R3 *ROO10 +R4 *ROO100 +R5 *ROO300 +R6 *ROO2M

      ENDDO

      CASE(10:11)                     !  Fully incompressible
       RO(1:NTOT) = FRSDEN

      CASE(12)                        !  Liquid R-12
       RO(1:NTOT) = 1017.1

      CASE(13)                        !  Gas R-12
       RO(1:NTOT) = 172.06

      CASE(101)                       !  Cast iron
       RO(1:NTOT) = .7272e+4

      END SELECT

      DEALLOCATE(TEMP,P)

      RETURN

      END SUBROUTINE ROFPT
C
C ......................................................................
C ......................................................................
C ......................................................................
C
      REAL FUNCTION LORE(ONEORZERO)

C ... If ONEORZERO is .TRUE. then return 1.0 else return 0.0.

      LOGICAL :: ONEORZERO

      IF(ONEORZERO) THEN
         LORE = 1.0
      ELSE
         LORE = 0.0
      ENDIF
    
      RETURN
      END FUNCTION LORE
C
C ......................................................................
C ......................................................................
C ......................................................................
C
      SUBROUTINE TFROE(RO,E,TEMP,NTOT,ISTATE,RGAS,
     &                 GAMMA,FRSPRE,E0REF,T0REF)

! ... TEMPERATURE AS A FUNCTION OF DENSITY AND INTERNAL ENERGY
! ... USED IN BOUNDARY CONDITIONS

      IMPLICIT NONE

      REAL,DIMENSION(*) :: TEMP,E,RO
      REAL :: GAMMA,RGAS,FRSPRE,E0REF,T0REF,GM1,CV
      INTEGER :: ERCODE,NTOT,ISTATE

      ERCODE = 0

      SELECT CASE(ISTATE)

      CASE(8:)
       PRINT"(//A/A//)",
     &      "ISTATE >= 8 does not have TFROE option.",
     &      "Exiting..." ; STOP

      CASE(1,3:5)       ! PERFECT GAS & SEMI-PERFECT GAS, PPR 12.5.1995
       GM1 = GAMMA - 1. ; TEMP(1:NTOT) = GM1/RGAS*E(1:NTOT)

      CASE(2)           ! Chemically Reacting Equilibrium Air
       PRINT"(//A//)","No such model of E.o.S. Exiting..." ; STOP

      CASE(6:7)         ! WATER & PWATER
       CV = .4180e+4    ! EI NYT IHAN OIKEA ARVO
       TEMP(1:NTOT) = .2930e+3 + (E(1:NTOT) - .8377e+5)/CV

      END SELECT

      END SUBROUTINE TFROE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EFPT(E,P,TEMP,NTOT,ISTATE,RGAS,
     &                GAMMA,FRSPRE,FRSSIE,E0REF,T0REF)

! ... Internal energy as a function of pressure and temperature
! ... used e.g. in BOTSKE to evaluate the state

      IMPLICIT NONE

      REAL,DIMENSION(*) :: TEMP,E,P
      REAL :: GAMMA,RGAS,FRSPRE,E0REF,T0REF,GM1,CV,FRSSIE
      INTEGER :: ERCODE,NTOT,ISTATE

      ERCODE = 0

      SELECT CASE(ISTATE)

      CASE(1,3:5)                    ! PERFECT GAS & SEMI-PERFECT GAS PPR 12.5.1995
       GM1 = GAMMA - 1. ; E(1:NTOT) = RGAS/GM1*TEMP(1:NTOT)

      CASE(2)                        ! Chemically Reacting Equilibrium Air
       PRINT"(//A//)","No such model for E.o.S. Exiting..." ; STOP

      CASE(6:7)                      ! WATER & PWATER joo joo
       CV = .418e+4                  ! EI NYT IHAN OIKEA ARVO
       E(1:NTOT) = CV*(TEMP(1:NTOT) - .293e+3) + .8377e+5

      CASE(8)                        ! Water P = 0.5 ... 1.5 bar
       E(1:NTOT) = ( .386586e+0*TEMP(1:NTOT)+.392553e+4)*TEMP(1:NTOT)
     &            +( .447332e-5*TEMP(1:NTOT)-.169196e-2)*P(1:NTOT)
     &            +(-.109967e+7)

      CASE(9)                        ! Steam P = 0.5 ... 1.5 bar
       E(1:NTOT) = (-.376298e-1*TEMP(1:NTOT)+.150939e+4)*TEMP(1:NTOT)
     &            +( .360841e-3*TEMP(1:NTOT)-.212573e+0)*P(1:NTOT)
     &            +( .195633e+7)

      CASE(10:11)                    ! Fully incompressible
       E(1:NTOT) = FRSSIE

      CASE(12)                       !  Liquid R-12
C       E(1:NTOT) = 1.2693E5 + 1419.8*(TEMP(1:NTOT) - 359.79) - 
C     &             P(1:NTOT)/1017.1

      E(1:NTOT)  = 1419.8*(TEMP(1:NTOT) - 298.15)
     +                 - P(1:NTOT)/1017.1
      E(1:NTOT)  = MAX(E(1:NTOT),1.)

      CASE(13)                       !  Gas R-12
c       E(1:NTOT) =  2.13E5 + 1290.8*(TEMP(1:NTOT) - 359.79) - 
c     &             P(1:NTOT)/172.06

      E(1:NTOT)  = 1419.8 * (359.79 - 298.15) + 86070. +
     + 1290.8 *(TEMP(1:NTOT) - 359.79) -  P(1:NTOT)/172.06
      E(1:NTOT)  = MAX(E(1:NTOT),1.)

      CASE(101)                      ! Cast iron
       CV = .42e+3 ; E(1:NTOT) = CV*TEMP(1:NTOT)

      END SELECT

      END SUBROUTINE EFPT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EFPTG(E,P,TEMP,HSAT,CP,TSAT,RO,PRO,NTOT,ISTATE,RGAS,
     &                GAMMA,FRSPRE,FRSSIE,E0REF,T0REF,IPHASE)

! ... Internal energy as a function of pressure and temperature
! ... This is tuned with the saturation line in multiphase calculation

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL,DIMENSION(*) :: TEMP,E,P,HSAT,CP,TSAT,RO
      REAL :: GAMMA,RGAS,FRSPRE,E0REF,T0REF,GM1,CV,FRSSIE
      INTEGER :: ERCODE,NTOT,ISTATE,IPHASE

      TYPE(PROPERTIES) PRO(*)

      ERCODE = 0 

      SELECT CASE(ISTATE)

      CASE(1,3:5)                    ! PERFECT GAS & SEMI-PERFECT GAS PPR 12.5.1995
       GM1 = GAMMA - 1. ; E(1:NTOT) = RGAS/GM1*TEMP(1:NTOT)

      CASE(2)                        ! Chemically Reacting Equilibrium Air
       PRINT"(//A//)","No such model for E.o.S. Exiting..." ; STOP

      CASE(6:7)                      ! WATER & PWATER joo joo
       CV = .418e+4                  ! EI NYT IHAN OIKEA ARVO
       E(1:NTOT) = CV*(TEMP(1:NTOT) - .293e+3) + .8377e+5

      CASE(8)                        ! Water P = 0.5 ... 1.5 bar
       E(1:NTOT) = ( .386586e+0*TEMP(1:NTOT)+.392553e+4)*TEMP(1:NTOT)
     &            +( .447332e-5*TEMP(1:NTOT)-.169196e-2)*P(1:NTOT)
     &            +(-.109967e+7)
c      E(1:NTOT)  = PRO(1:NTOT)%HSAT(IPHASE) -
c     &           P(1:NTOT)/PRO(1:NTOT)%RO(IPHASE)+
c     &           PRO(1:NTOT)%CP(IPHASE)*(TEMP(1:NTOT)-PRO(1:NTOT)%TSAT)

c      E(1:NTOT)  = HSAT(1:NTOT) !- P(1:NTOT)/RO(1:NTOT)
c     &             + CP(1:NTOT)*(TEMP(1:NTOT)-TSAT(1:NTOT))

      CASE(9)                        ! Steam P = 0.5 ... 1.5 bar
       E(1:NTOT) = (-.376298e-1*TEMP(1:NTOT)+.150939e+4)*TEMP(1:NTOT)
     &            +( .360841e-3*TEMP(1:NTOT)-.212573e+0)*P(1:NTOT)
     &            +( .195633e+7)
c      E(1:NTOT)  = PRO(1:NTOT)%HSAT(IPHASE) -
c     &           P(1:NTOT)/PRO(1:NTOT)%RO(IPHASE)+
c     &           PRO(1:NTOT)%CP(IPHASE)*(TEMP(1:NTOT)-PRO(1:NTOT)%TSAT)

c      E(1:NTOT)  = HSAT(1:NTOT) !- P(1:NTOT)/RO(1:NTOT)
c     &             + CP(1:NTOT)*(TEMP(1:NTOT)-TSAT(1:NTOT))

      CASE(10:11)                    ! Fully incompressible
       E(1:NTOT) = FRSSIE

      CASE(12)                       !  Liquid R-12
c       E(1:NTOT) = 1.2693E5 + 1419.8*(TEMP(1:NTOT) - 359.79) - 
c     &             P(1:NTOT)/1017.1

      E(1:NTOT)  = 1419.8*(TEMP(1:NTOT) - 298.15)
     +                 - P(1:NTOT)/1017.1
      E(1:NTOT)  = MAX(E(1:NTOT),1.)

      CASE(13)                       !  Gas R-12
c       E(1:NTOT) = 2.13E5 + 1290.8*(TEMP(1:NTOT) - 359.79) - 
c     &             P(1:NTOT)/172.06

      E(1:NTOT)  = 1419.8 * (359.79 - 298.15) + 86070. +
     + 1290.8*(TEMP(1:NTOT) - 359.79)- P(1:NTOT)/172.06
      E(1:NTOT)  = MAX(E(1:NTOT),1.)

      CASE(101)                      ! Cast iron
       CV = .42e+3 ; E(1:NTOT) = CV*TEMP(1:NTOT)

      END SELECT

      END SUBROUTINE EFPTG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE STEAM3(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &                  IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input: T = temperature, [K]
! ...        P = pressure, [Pa]

      IMPLICIT NONE

      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,
     &           II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL
      REAL, DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: YC2,CG,CL,DRDT,DRDPA,CP0,CP1,CP5,F,A1,A2,AA,F1,F5
      REAL :: F00,F11,F12,F21,F22,F31,F32,F41,F42,F51,F52,F61,F62
      REAL :: R1,R2,R3,R4,R5,R6
      REAL :: ROO0,ROO1,ROO10,ROO100,ROO300,ROO2M
      REAL :: DRDT0,DRDT1,DRDT10,DRDT100,DRDT300,DRDT2M
      REAL :: FD11,FD12,FD21,FD22,FD31,FD32,FD41,FD42,FD51,FD52
      REAL :: FD61,FD62,RD1,RD2,RD3,RD4,RD5,RD6
      REAL :: R2P,R3P,R4P,R5P,R6P
      REAL :: SLP1,SLP2
      REAL :: RLAM1,RLAM2,RLAM5,PR,LORE

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

!      DO I = 1,NTOT

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

! ... specific heat of gas, [J/(kgK)]

      F00     = 1e-5*MAX(10.,P(I))
      F12     = 1. - F00
      F21     =      F00
      F00     = .25*(F00 - 1.)
      F22     = 1. - F00
      F31     =      F00
      F32     = 1.
      CP0     = ( .944056e-3*TEMP(I)-.294568e+0)*TEMP(I)+.186841e+4
      CP1     = ( .113393e-1*TEMP(I)-.102209e+2)*TEMP(I)+.427457e+4
      CP5     = ( .362500e-1*TEMP(I)-.360791e+2)*TEMP(I)+.110935e+5
      CP(I)   = + MAX(0.,F12)         *CP0
     &          + MAX(0.,MIN(F21,F22))*CP1
     &          + MAX(0.,MIN(F31,F32))*CP5

! ... density of gas, [kg/m^3] and its derivatives

! ... Start of new addition by A. Miettinen, 2007/10/01

! ... 0kPa ... tahan voi laittaa p=0 arvot, jos on tarvis
CAM      F11 = 1.
CAM      F12 = 1. - P/1000.
CAM      R1  = MAX(0.,MIN(F11,F12))
CAM      ROO0 = 0.
CAM      DRDT0 = 0.   

! ... 1kPa
      SLP1 = 0.001
      SLP2 = 1.0/9.0e+3
      F21  =            SLP1*MAX(10.,P(I))
      F00  = 0.111111 - SLP2*MAX(10.,P(I))
      F22  = 1.0 + F00
      R2   = MAX(0.,MIN(F21,F22))
!      ROO1    =  -.134939e-4*TEMP(I)+.111220e-1
      ROO1    = ( .322363e-7*TEMP(I)-.406131e-4)*TEMP(I)+.165018e-1
      DRDT1   =   .644726e-7*TEMP(I)-.406131e-4
CAM      R2P     = LORE(R2 > 0.0)*(SLP1*LORE(F21 < F22)
CAM     &                        - SLP2*LORE(F21 > F22))

! ... 10kPa
      SLP1 = SLP2
      SLP2 = 1.0/9.0e+4
      F31  = -F00
      F00  = 0.111111 - SLP2*P(I)
      F32  = 1.0 + F00
      R3   = MAX(0.,MIN(F31,F32))
!      ROO10   =  -.119482e-3*TEMP(I)+.103647e+0
      ROO10   = ( .281416e-6*TEMP(I)-.367450e-3)*TEMP(I)+.156294e+0
      DRDT10  =   .562832e-6*TEMP(I)-.367450e-3
CAM      R3P     = LORE(R3 > 0.0)*(SLP1*LORE(F31 < F32)
CAM     &                        - SLP2*LORE(F31 > F32))

! ... 100kPa
      SLP1 = SLP2
      SLP2 = 1.0/2.0e+5
      F41  = -F00
      F00  = 0.500000 - SLP2*P(I)
      F42  = 1.0 + F00
      R4   = MAX(0.,MIN(F41,F42))
!      ROO100  =  -.105576e-2*TEMP(I)+.972238e+0
      ROO100  = ( .240937e-5*TEMP(I)-.331947e-2)*TEMP(I)+.149104e+1
      DRDT100 =   .481874e-5*TEMP(I)-.331947e-2
CAM      R4P     = LORE(R4 > 0.0)*(SLP1*LORE(F41 < F42)
CAM     &                        - SLP2*LORE(F41 > F42))

! ... 300kPa
      SLP1 = SLP2
      SLP2 = 1.0/1.7e+6
      F51  = -F00
      F00  = 0.176471 - SLP2*P(I)
      F52  = 1.0 + F00
      R5   = MAX(0.,MIN(F51,F52))
!      ROO300  =  -.304522e-2*TEMP(I)+.286537e+1
      ROO300  = ( .726955e-5*TEMP(I)-.101508e-1)*TEMP(I)+.457358e+1
      DRDT300 =   .145391e-4*TEMP(I)-.101508e-1
CAM      R5P     = LORE(R5 > 0.0)*(SLP1*LORE(F51 < F52)
CAM     &                        - SLP2*LORE(F51 > F52))

! ... 2Mpa
      SLP1 = SLP2
      SLP2 = 0.0
      F61  = -F00
      F62  = 1.0
      R6   = MAX(0.,MIN(F61,F62))
!      ROO2M   =  -.151110e-1*TEMP(I)+.169803e+2
      ROO2M   = ( .348933e-4*TEMP(I)-.585389e-1)*TEMP(I)+.301263e+2
      DRDT2M  =   .697866e-4*TEMP(I)-.585389e-1
CAM      R6P     = LORE(R6 > 0.0)*(SLP1*LORE(F61 < F62)
CAM     &                        - SLP2*LORE(F61 > F62))

      RO(I) = R2 *ROO1 +R3 *ROO10 +R4 *ROO100 +R5 *ROO300 +R6 *ROO2M

      DRDT  = R2 *DRDT1+R3 *DRDT10+R4 *DRDT100+R5 *DRDT300+R6 *DRDT2M

C ... seuraavan saannonmukaisuuksia voi siistia liukkaammiksi ja likiarvot 
C ... 0.379310345 = 55./145. = 11./29.
C ... 0.210526316 = 20./95.  = 4./19.
      FD11 = + 1.
      FD12 = + 1.1 - P(I)/5000.
      RD1  = MAX(0.,MIN(FD11,FD12))
      FD21 = - 0.1       + P(I)/5000.
      FD22 = + 1.1111111 - P(I)/49500.
      RD2  = MAX(0.,MIN(FD21,FD22))
      FD31 = - 0.111111111 + P(I)/49500.
      FD32 = + 1.379310345 - P(I)/145000.
      RD3  = MAX(0.,MIN(FD31,FD32))
      FD41 = - 0.379310345 + P(I)/145000.
      FD42 = + 1.210526316 - P(I)/950000.
      RD4  = MAX(0.,MIN(FD41,FD42))
      FD51 = -0.210526316 + P(I)/950000.
      FD52 = 1.
      RD5  = MAX(0.,MIN(FD51,FD52))
C ... FD54 = 1. eli asetetaan 300kPa<P<2MPa derivaatan vakioarvo kun P>1150kPa
C ... seuraavat jos arvo asetetaan johonkin tosi isolla paineella
C      FD61 = 
C      FD62 = 
C      RD6  = MAX(0.,MIN(FD61,FD62))

C ... seuraavan avulla pultataan p=0. jokin ROO(T) funktio
      ROO0 = 0.
      DRDPA = + (- .001*ROO0
     &           + .001*ROO1)*RD1
     &        + (- .000111111*ROO1 
     &           + .000111111*ROO10)*RD2
     &        + (- .0000111111*ROO10 
     &           + .0000111111*ROO100)*RD3
     &        + (- .000005*ROO100    
     &           + .000005*ROO300)*RD4
     &        + (- .000000588*ROO300 
     &           + .000000588*ROO2M)*RD5
C     &        + 0.*RD6
C ... tuo viimeinen kommentoitu on paineen suhteen vakiona pysyvän tiheyden syyta
C ... kaytannossa 
C ... kun p<500  Pa  DRDPT=vakio
C ... kun p>1150 kPa DRDPT=vakio
      
CAM      DRDPA = R2P*ROO1 +R3P*ROO10 +R4P*ROO100 +R5P*ROO300 +R6P*ROO2M

! ... End of new addition by A. Miettinen, 2007/10/01

!      RO(I)   = ( .413508e-5*TEMP(I)-.395082e-2)*TEMP(I)+.901454e+0
!     &        + (-.104028e-7*TEMP(I)+.976698e-5)*P(I)

!      DRDT    =   .827016e-5*TEMP(I)-.395082e-2-.104028e-7*P(I)

!      DRDPA   =  -.104028e-7*TEMP(I)+.976698e-5

      DRDH(I) = DRDT/CP(I)
      DRDP(I) = DRDPA - DRDH(I)*(1.+ TEMP(I)*DRDT/RO(I))/RO(I)

      YC2     = DRDH(I)/RO(I) + DRDP(I)
      CG      = 1./SQRT(YC2)

      END DO ; END DO ; END DO

!      END DO

!      PR      = 0. ! An average Prandtl number could be evaluated blockwise

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

! ... dynamic viscosity of gas, [kg/(m s)]

      VIS(I)  = ( .016000e-9*TEMP(I)+.254589e-7)*TEMP(I)+.551790e-6

! ... thermal conductivity of gas, [W/(m K)]

      F00     = .25*(1. - 1e-5*P(I))
      F11     = 1.
      F12     = 1.+F00
      F21     =   -F00
      F22     = 1.
      F1      = MAX(0.,MIN(F11,F12))
      F5      = MAX(0.,MIN(F21,F22))
      RLAM1   = ( .132000e-6*TEMP(I)-.297141e-4)*TEMP(I)+.177873e-1
      RLAM5   = ( .832596e-7*TEMP(I)+.204684e-5)*TEMP(I)+.153213e-1
      CH(I)   = F1*RLAM1 + F5*RLAM5

! ... Average Prandtl number for this loop

!      PR      = ((NTOT-1.)*PR + VIS(I)*CP(I)/CH(I))/NTOT

! ... internal energy of gas, [J/(kg K)]

      E(I)    = (-.376298e-1*TEMP(I)+.150939e+4)*TEMP(I)+.195633e+7
     &        + ( .360841e-3*TEMP(I)-.212573e+0)*P(I)

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE STEAM3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WATER2(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &                  IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      IMPLICIT NONE

      REAL,DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IL
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

!         write(666,*) ii,jj,kk,i
! ... specific heat of liquid, [J/(kgK)]

      CP(I)   = MAX(.4182e+4,
     &          ( .560839e-2*TEMP(I)-.284106e+1)*TEMP(I)+.450010e+4)

! ... density of liquid, [kg/m^3]

      RO(I)   = (-.260031e-2*TEMP(I)+.122284e+1)*TEMP(I)+.864232e+3
     &        + ( .761625e-9*TEMP(I)+.196854e-6)*P(I)

      DRDT    =  -.520062e-2*TEMP(I)+.122284e+1+.761625e-9*P(I)
      DRDPA   =   .761625e-9*TEMP(I)+.196854e-6

      DRDH(I) = DRDT/CP(I)
      DRDP(I) = DRDPA - DRDT*(1.+ TEMP(I)*DRDT/RO(I))/(RO(I)*CP(I))

      YC2     = DRDH(I)/RO(I) + DRDP(I)
      CL      = 1./SQRT(YC2)

      END DO ; END DO ; END DO

!      END DO

!      PR      = 0. ! An average Prandtl number could be evaluated blockwise

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta TEMP(I)
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(TEMP(I)-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*TEMP(I)+.503834e-4)*TEMP(I)
     &  -.233632e-1        +.369516e+1/TEMP(I))
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - TEMP(I))))
	G22 = 1.0 - G12

	VIS(I)  =  G12 * G11 + G22 * G21

	! ---------------------------------

! ... thermal conductivity of liquid, [W/Km]

      F11     = 0.1*(.385e+3-TEMP(I))
      F21     = 1.0-F11
      F12     = MAX(0.,MIN(1.,F11))
      F22     = MAX(0.,MIN(1.,F21))
      RLAM1   = (-.92114e-5*TEMP(I)+.713593e-2)*TEMP(I)-.700822e+0
      RLAM2   =   .68085e+0
      CH(I)   = F12*RLAM1 + F22*RLAM2
      CH(I)   = MAX(.56,CH(I))

! ... Average Prandtl number for this loop

!      PR      = ((NTOT-1.)*PR + VIS(I)*CP(I)/CH(I))/NTOT

! ... internal energy of liquid, [J/(kgK)]

      E(I)   = ( .386586e+0*TEMP(I)+.392553e+4)*TEMP(I)-.109967e+7
     &       + ( .447332e-5*TEMP(I)-.169196e-2)*P(I)

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE WATER2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE R12LIQ(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &                  IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      IMPLICIT NONE

      REAL,DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IL

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

! ... specific heat of liquid R-12, [J/(kgK)]

      CP(I)   = 1419.8

! ... density of liquid R-12, [kg/m^3]

      RO(I)   = 1017.1

      DRDT    = -1.E-4
      DRDPA   =  .5E-6

      DRDH(I) = DRDT
      DRDP(I) = DRDPA

      END DO ; END DO ; END DO

!      END DO

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

! ... dynamic viscosity of liquid, [kg/(sm)]

      VIS(I)  = 8.95E-5

! ... thermal conductivity of liquid, [W/Km]

      CH(I)   = 4.57E-2

! ... Average Prandtl number for this loop

!      PR      = ((NTOT-1.)*PR + VIS(I)*CP(I)/CH(I))/NTOT

! ... internal energy of liquid, [J/(kgK)]

      E(I)   = 1.2693E5 + CP(I)*(TEMP(I) - 359.79) - P(I)/RO(I)
      E(I)  = MAX(E(I),1.)


      END DO ; END DO ; END DO

      END SUBROUTINE R12LIQ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE R12GAS(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
     &                  IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      IMPLICIT NONE

      REAL,DIMENSION(*) :: P,RO,DRDP,DRDH,CP,VIS,CH,TEMP,E
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IL

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

! ... specific heat of saturated gas R-12, [J/(kgK)]

      CP(I)   = 1290.8

! ... density of saturated gas, [kg/m^3]

      RO(I)   = 172.06

      DRDT    = -1.E-4
      DRDPA   =  1.E-6

      DRDH(I) = DRDT
      DRDP(I) = DRDPA

      END DO ; END DO ; END DO

!      END DO

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

! ... dynamic viscosity of saturated gas, [kg/(sm)]

      VIS(I)  = 1.572E-5

! ... thermal conductivity of saturated gas, [W/Km]

      CH(I)   = 1.76E-2

! ... Average Prandtl number for this loop

!      PR      = ((NTOT-1.)*PR + VIS(I)*CP(I)/CH(I))/NTOT

! ... internal energy of liquid, [J/(kgK)]

      E(I)   = 2.13E5 + CP(I)*(TEMP(I) - 359.79) - P(I)/RO(I)
      E(I)  = MAX(E(I),1.)

      END DO ; END DO ; END DO

      END SUBROUTINE R12GAS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CAST_IRON(RO,E,TEMP,CP,CH,NTOT,
     &                     IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input : T = temperature, [K]

      IMPLICIT NONE

      REAL,DIMENSION(*) :: RO,CP,CH,TEMP,E
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN
      INTEGER :: I,II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

! ... specific heat of cast iron, [J/(kgK)]
      CP(I)   = .4200e+3

! ... density, [kg/m^3]
      RO(I)   = .7272e+4

! ... internal energy, [J/(kgK)]  (oli [J/(m^3K)])
      E(I)    = CP(I)*TEMP(I)

! ... thermal conductivity, [W/Km]
      CH(I)   = .51e+2 - .40e-1*TEMP(I)

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE CAST_IRON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DISTPRO(RO,E,TEMP,CP,CH,NTOT,
     &                     IMAX,JMAX,KMAX,IN,JN,KN,INN)

! ... Input : Nothing, outputs = 1

      IMPLICIT NONE

      REAL,DIMENSION(*) :: RO,CP,CH,TEMP,E
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN
      INTEGER :: I,II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

! ... specific heat of cast iron, [J/(kgK)]
      CP(I)   = 1.

! ... density, [kg/m^3]
      RO(I)   = 1.

! ... internal energy, [J/(kgK)]  (oli [J/(m^3K)])
      E(I)    = 1.

! ... thermal conductivity, [W/Km]
      CH(I)   = 1.

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE DISTPRO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PDTOP(P,PDIFF,XC,YC,ZC,FRSDEN,FRSPRE,NTOT)

! ... Calculate the pressure from the pressure difference and gravity
! ... GROUND is now aligned with the main axis (requires GROUNDX, etc.)

      USE NS3CO, ONLY : GX, GY, GZ, GROUND, ALTITUDE, TWO_FLUIDL

      IMPLICIT NONE

      INTEGER :: NTOT
      REAL, DIMENSION(*) :: XC,YC,ZC
      REAL, DIMENSION(*) :: P,PDIFF
      REAL :: FRSDEN,FRSPRE,H

c      IF(ALTITUDE >= GROUND .OR. TWO_FLUIDL) THEN
      IF(ALTITUDE >= GROUND) THEN
         P(1:NTOT) = PDIFF(1:NTOT) + FRSPRE
      ELSE
          P(1:NTOT) = PDIFF(1:NTOT) + FRSPRE
     &              + FRSDEN*(GX*(XC(1:NTOT) - GROUND) ! Mersu
     &              +         GY*(YC(1:NTOT) - GROUND)
     &              +         GZ*(ZC(1:NTOT) - GROUND))    
      ENDIF
c       write(667,*) P(25820),pdiff(25820),FRSPRE,FRSDEN,YC(25820),GROUND
      END SUBROUTINE PDTOP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MFDENS(PRO,VAR,PDIFF,P,RO,TEMP,DRDP,DRDH,VIS,CH,CP,
     &                  FRSTEM,NPHASE,IMAX,JMAX,KMAX,IN,JN,KN,INN,
     &                  ICYCLE,MULPHC,FRSPRE,DFRSPRE,DFRSTEM,IEVAP,
     &                  ISTATE,DFRSDEN,XC,YC,ZC,NGL,IREPEA,NREPEA)

! ... Calculate void fraction, density and temperature for the mixture
! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,INN,IN,JN,KN,NPHASE,ICYCLE,IEVAP,IL,NGL
      INTEGER :: II,JJ,KK,I,J,K,JA,KA,INN1,ISTRID,JSTRID,IPHASE,ISTATE
      INTEGER,DIMENSION(*) :: IREPEA,NREPEA
      TYPE(PROPERTIES),      DIMENSION(*) :: PRO
      TYPE(MPHASE_VARIABLES),DIMENSION(*) :: VAR
      REAL,DIMENSION(*) :: P,RO,TEMP,DRDP,DRDH,VIS,CP,CH,PDIFF
      REAL :: RKOE,FRSTEM,FRSPRE,DPSATDT,GTOT,ALPO
      REAL :: DFRSPRE,DFRSTEM,TEMD,PRESD,TS1,TS2,DPH,DFRSDEN,
     &                 DGROUND,dph2
      REAL, DIMENSION(*) :: XC,YC,ZC
      REAL, PARAMETER :: SMALLVOID = 1.e-14
      CHARACTER(LEN=10) :: MULPHC
      REAL, PARAMETER ::
     &    P1 =  0.154906419868284e+6,
     &    Q1 = -0.101158884083608e+5,
     &    R1 =  0.133000059771258e+3,
     &    S1 =  0.146861326101013e-1,
     &    T1 = -0.170477189812508e+2,
     &    V1 =  0.249899160714079e+2,
     &    W1 = -0.504070497017732e+4
      LOGICAL :: CORNER,INSIDE

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      INN1    = INN-1
        alpo = 0.
      IL     = ISTRID*JSTRID

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

         CORNER = (KK==-1.AND.JJ==-1 .OR. KK==-1.AND.II == -1      .OR.
     &   JJ==-1.AND.II==-1                                         .OR.
     &   KK==-1.AND.JJ==JMAX+2 .OR. KK==-1.AND.II==IMAX+2          .OR.
     &   JJ==-1.AND.II==IMAX+2                                     .OR.
     &   KK==KMAX+2.AND.JJ==-1 .OR. KK==KMAX+2.AND.II == -1        .OR.
     &   JJ==JMAX+2.AND.II==-1                                     .OR.
     &   KK==KMAX+2.AND.JJ==JMAX+2 .OR. KK==KMAX+2.AND.II==IMAX+2  .OR.
     &   JJ==JMAX+2.AND.II==IMAX+2)    

         INSIDE = (JJ >= 1.AND.JJ <= JMAX .AND. KK >= 1.AND.KK <= KMAX
     &   .OR.      II >= 1.AND.II <= IMAX .AND. KK >= 1.AND.KK <= KMAX
     &   .OR.      II >= 1.AND.II <= IMAX .AND. JJ >= 1.AND.JJ <= JMAX)

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      RO(I) = 0.

! ... Nullify small voids
      DO IPHASE = 1,NPHASE
       IF (   VAR(I)%X(IPHASE) < SMALLVOID) VAR(I)%X(IPHASE) = 0.
       IF (1.-VAR(I)%X(IPHASE) < SMALLVOID) VAR(I)%X(IPHASE) = 1.
       RO(I) = RO(I) + VAR(I)%X(IPHASE)/PRO(I)%RO(IPHASE)
      END DO
      RO(I) = 1./RO(I)

      TEMP(I) = 0.
      DRDP(I) = 0.
      DRDH(I) = 0.
      VIS(I)  = 0.
      CP(I)   = 0.
      CH(I)   = 0.

      IF(VAR(I)%ALFA(2) <= 0.5) THEN ! Estimate from the previous value

C ... The following is more accurate for small voids

      ALPO    = VAR(I)%X(2)*RO(I)/PRO(I)%RO(2)
          VAR(I)%ALFA(1) = 1. - ALPO
          VAR(I)%ALFA(2) = ALPO

       ELSE

C ... The following is more accurate for high voids

          ALPO    = VAR(I)%X(1)*RO(I)/PRO(I)%RO(1)
          VAR(I)%ALFA(2) = 1. - ALPO
          VAR(I)%ALFA(1) = ALPO

       ENDIF

      DO IPHASE = 1,NPHASE
       VAR(I)%ARO(IPHASE)  = VAR(I)%X(IPHASE)*RO(I)
       DRDP(I) = DRDP(I)   + VAR(I)%X(IPHASE)*PRO(I)%DRODP(IPHASE)
       DRDH(I) = DRDH(I)   + VAR(I)%X(IPHASE)*PRO(I)%DRODH(IPHASE)
C ... Above DRDH-dependence was changed to 'X' (8.10.2010)
       TEMP(I) = TEMP(I)   + VAR(I)%ALFA(IPHASE)*PRO(I)%DTEMP(IPHASE)
       CP(I)   = CP(I)     + VAR(I)%ALFA(IPHASE)*PRO(I)%CP(IPHASE)
       VIS(I)  = VIS(I)    + VAR(I)%ALFA(IPHASE)*PRO(I)%VIS(IPHASE)
       CH(I)   = CH(I)     + VAR(I)%ALFA(IPHASE)*PRO(I)%CH(IPHASE)
      END DO

      IF (MULPHC == 'CAVIT' .OR. MULPHC == 'MULTI') THEN
c       PRO(I)%PSAT  = (( .305416e+0 *FRSTEM -.283410e+3)*FRSTEM
c     &                  +.878017e+5)*FRSTEM -.907075e+7
c       PRO(I)%DPSAT = (( .305416e+0 *DFRSTEM-.283410e+3)*DFRSTEM
c     &                  +.878017e+5)*DFRSTEM-.907075e+7 -DFRSPRE

      SELECT CASE(IEVAP)

      CASE(0,1,2) ! These are based on a constant temperature

c         DGROUND     = -REAL(GROUND,DP)
c         DPH         = DFRSDEN*(REAL(GX,DP)*(XC(I)-DGROUND) +
c     &                          REAL(GY,DP)*(YC(I)-DGROUND) +
c     &                          REAL(GZ,DP)*(ZC(I)-DGROUND))

         DGROUND     = GROUND  ! *DFRSDEN*GY

         DPH         = DFRSDEN*(GX*(XC(I) - DGROUND) +
     &                          GY*(YC(I) - DGROUND) +
     &                          GZ*(ZC(I) - DGROUND))  ! + DGROUND

c          DPH         = DFRSDEN*(REAL(GX,DP)*XC(I) +
c     &                           REAL(GY,DP)*YC(I) +
c     &                           REAL(GZ,DP)*ZC(I) + 9.81*DGROUND) !+ DGROUND
c      if(ngl == 2 .and. kk == 45) then
c      write(766,*) i,ii,jj,real(DPH),real(dph2),real(DPH+GTOT),-gtot
c      endif

         TEMD        = DFRSTEM
         PRESD       = DFRSPRE + PDIFF(I) + DPH ! Oli -
         PRO(I)%PSAT = EXP((P1/TEMD+Q1)/TEMD + R1 + S1*TEMD + 
     &                 T1*LOG(TEMD))
         PRO(I)%PSAT = MAX(.1e-4,PRO(I)%PSAT)
         PRO(I)%DPSAT= PRO(I)%PSAT - DFRSPRE - DPH ! Oli +
      CASE(3,4,5,6,7)   ! Liquid temperature?
         TEMD        = PRO(I)%DTEMP(1)
         PRESD       = DFRSPRE + PDIFF(I) + DPH

         PRO(I)%PSAT = EXP((P1/TEMD+Q1)/TEMD + R1 + S1*TEMD + 
     &                 T1*LOG(TEMD))
         PRO(I)%PSAT = MAX(.1e-4,PRO(I)%PSAT)
         PRO(I)%DPSAT= PRO(I)%PSAT - DFRSPRE

      CASE(12,13) ! R-12
         PRO(I)%PSAT = DFRSPRE
         PRO(I)%PSAT = MAX(.1e-4,PRO(I)%PSAT)
         PRO(I)%DPSAT= PRO(I)%PSAT - DFRSPRE
      END SELECT


      ELSE ! Two-fluid modelin its infancy?

      PRO(I)%PSAT = (( .305416e+0 *TEMP(I)-.283410e+3)*TEMP(I)
     &                  +.878017e+5)*TEMP(I)-.907075e+7
      PRO(I)%PSAT = MAX(.1e-4,PRO(I)%PSAT)
      PRO(I)%DPSAT= PRO(I)%PSAT - DFRSPRE
      END IF


      SELECT CASE(ISTATE)

      CASE(8,9) ! Water and steam
c      PRO(I)%PSAT = MAX(.1e-4,PRO(I)%PSAT)
c      PRO(I)%TSAT =  .297860e+0*SQRT(P(I))
c     &            +( .070623e-9*P(I)-.202998e-3)*P(I)+.298084e+3
         PRESD       = DFRSPRE + PDIFF(I) + DPH
         IF(PRESD <= -0.01*DFRSPRE .AND. INSIDE) THEN
           IF(IREPEA(11) <= NREPEA(11)) THEN
C           WRITE(*,*) ' Pressure less than -0.01*FRSPRE in block', NGL
           WRITE(13,*) ' Pressure less than zero in block', NGL,II,JJ,KK
           WRITE(13,*) ' Not Exiting from MFDENS'
           IF(IREPEA(11)==NREPEA(11)) 
     &     WRITE(13,*)' Silence from cycle ',ICYCLE 
           IREPEA(11) = IREPEA(11) + 1
           ENDIF
           PRESD = 0.01*DFRSPRE
         ENDIF
         PRO(I)%TSAT = W1/(LOG(PRESD) - V1)
         TS1         = W1/(LOG(PRESD+1.) - V1)
         TS2         = W1/(LOG(PRESD-1.) - V1)
         PRO(I)%DPSDT= 1./((TS1-TS2)/2.)
         PRO(I)%DTSDP= (TS1-TS2)/2.
      CASE(12,13) ! Freon R-12
         PRO(I)%TSAT = 86.73 + 273.15
         PRO(I)%DPSDT= 1.
         PRO(I)%DTSDP= 1.
      END SELECT

      END DO ; END DO ; END DO

      END SUBROUTINE MFDENS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

! **********************************************************************
! ... INITIALIZATION ROUTINES
! **********************************************************************

      SUBROUTINE BCINI1(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &                  SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &                  GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,REFVEL,
     &                  GROUND,ALTITUDE)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON PERFECT GAS MODEL

      IMPLICIT NONE

      REAL :: FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &        SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &        GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,RMUFRS,GM1,
     &        REFVEL,GROUND,ALTITUDE

      REAL, EXTERNAL :: ATMOSRHO, ATMOSP, ATMOST, ATMOSA


      IF(ALTITUDE >= GROUND) THEN  ! Free stream conditions based on 
                                   ! flight altitude

      GM1    = GAMMA-1.
      FRSDEN = ATMOSRHO(ALTITUDE)                           
      FRSPRE = ATMOSP(ALTITUDE)
      FRSTEM = ATMOST(ALTITUDE)
      FRSSSP = ATMOSA(ALTITUDE)
      FRSSIE = -E0REF + RGAS/GM1*(FRSTEM-T0REF)
      FRSVIS = VISU0*FRSTEM**EXPSU/(FRSTEM+TSU0)
      SIEINI = FRSSIE

      RMUFRS = VISU0*FRSTEM**EXPSU/(FRSTEM+TSU0)*RMULT

      IF(FRSVEL == 0. .AND. RMACH /= 0.) THEN
         FRSVEL = RMACH*FRSSSP
      WRITE(4,"(2X,A)") "Free-stream quantities based on given RMACH,"
      ELSE IF (FRSVEL /= 0. .AND. RMACH == 0.) THEN
      RMACH = FRSVEL/FRSSSP
      WRITE(4,"(2X,A)") "Free-stream quantities based on given FRSVEL,"
      ELSE IF (FRSVEL == 0. .AND. RMACH == 0. .AND. REFVEL /= 0.) THEN
      RMACH = REFVEL/FRSSSP
      WRITE(4,"(2X,A)") "Free-stream quantities based on given REFVEL,"
      ELSE
      PRINT"(A/A,2F11.5/A)",
     &    "Problem is overdefined. Only one velocity scale may",
     &    "be specified. RMACH and FRSVEL are ",RMACH,FRSVEL,
     &    "Exiting..." ; STOP
      END IF

      IF(RE == 0. .AND. CHLREF /= 0.) THEN
       IF(FRSVEL /= 0.) THEN
        RE = FRSDEN*FRSVEL*CHLREF/RMUFRS
       ELSE
        RE = FRSDEN*REFVEL*CHLREF/RMUFRS
       ENDIF
       WRITE(4,"(A)") " FRSTEM, FRSDEN AND CHLREF"
      ELSEIF(RE /= 0. .AND. CHLREF == 0.) THEN
       IF(FRSVEL /= 0.) THEN
        CHLREF = RE*RMUFRS/(FRSDEN*FRSVEL)
       ELSE
        CHLREF = RE*RMUFRS/(FRSDEN*REFVEL)
       ENDIF
       WRITE(4,"(A)") " FRSTEM, RE AND FRSDEN"
      ELSE IF (RE == 0. .AND. CHLREF /= 0.) THEN
       IF(FRSVEL /= 0.) THEN
        RE = FRSDEN*FRSVEL*CHLREF/RMUFRS
       ELSE
        RE = FRSDEN*REFVEL*CHLREF/RMUFRS
       ENDIF
       WRITE(4,"(A)") "  FRSTEM, FRSPRE AND CHLREF"
      ELSE
       PRINT"(A/A/A,3F11.5/A)",
     &    "Problem is overdefined. Reynolds number ",
     &    "or the reference length must not be specified.",
     &    "RE and CHLREF are ",RE,CHLREF,
     &    "Exiting..." ; STOP
      END IF

      T0 = FRSTEM*(1.+.5*GM1*RMACH**2)

      ELSE  ! Free stream conditions based on perfect gas model

      IF (FRSDEN == 0. .AND. FRSPRE /= 0.) THEN ! Added 3.1.2006
       FRSDEN = FRSPRE/(RGAS*(FRSTEM-T0REF))
      END IF
      FRSSSP  = SQRT(GAMMA*RGAS*FRSTEM)
! ... Sutherland's Formula
      RMUFRS  = VISU0*FRSTEM**EXPSU/(FRSTEM+TSU0)*RMULT
      IF (FRSVEL == 0. .AND. RMACH /= 0.) THEN
       FRSVEL = RMACH*FRSSSP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given RMACH,"
      ELSE IF (FRSVEL /= 0. .AND. RMACH == 0.) THEN
       RMACH = FRSVEL/FRSSSP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given FRSVEL,"
      ELSE IF (FRSVEL == 0. .AND. RMACH == 0. .AND. REFVEL /= 0.) THEN
       RMACH = REFVEL/FRSSSP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given REFVEL,"
      ELSE
       PRINT"(A/A,2F11.5/A)",
     &    "Problem is overdefined. Only one velocity scale may",
     &    "be specified. RMACH and FRSVEL are ",RMACH,FRSVEL,
     &    "Exiting..." ; STOP
      END IF

      IF (FRSDEN == 0. .AND. RE /= 0. .AND. CHLREF /= 0.) THEN
       IF(FRSVEL /= 0.) THEN
       FRSDEN = RE*RMUFRS/(FRSVEL*CHLREF)
       ELSE
        FRSDEN = RE*RMUFRS/(REFVEL*CHLREF)
       ENDIF
       WRITE(4,"(A)") " FRSTEM, RE AND CHLREF"
      ELSE IF (FRSDEN /= 0. .AND. RE == 0. .AND. CHLREF /= 0.) THEN
       IF(FRSVEL /= 0.) THEN
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       ELSE
        RE     = FRSDEN*REFVEL*CHLREF/RMUFRS
       ENDIF
       WRITE(4,"(A)") " FRSTEM, FRSDEN AND CHLREF"
      ELSE IF (FRSDEN /= 0. .AND. RE /= 0. .AND. CHLREF == 0.) THEN
       IF(FRSVEL /= 0.) THEN
       CHLREF = RE*RMUFRS/(FRSDEN*FRSVEL)
       ELSE
        CHLREF = RE*RMUFRS/(FRSDEN*REFVEL)
       ENDIF
       WRITE(4,"(A)") " FRSTEM, RE AND FRSDEN"
      ELSE IF (FRSDEN == 0. .AND. RE == 0. .AND. CHLREF /= 0.
     &                                     .AND. FRSPRE /= 0.) THEN
       FRSDEN = FRSPRE/(RGAS*(FRSTEM-T0REF))
       IF(FRSVEL /= 0.) THEN
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       ELSE
        RE     = FRSDEN*REFVEL*CHLREF/RMUFRS
       ENDIF
       WRITE(4,"(A)") "  FRSTEM, FRSPRE AND CHLREF"
      ELSE
       PRINT"(A/A/A,3F11.5/A)",
     &    "Problem is overdefined. Density, the Reynolds number ",
     &    "or the reference length must not be specified.",
     &    "RE, FRSDEN and CHLREF are ",RE,FRSDEN,CHLREF,
     &    "Exiting..." ; STOP
      END IF

      GM1    = GAMMA-1.
      FRSPRE = RGAS*(FRSTEM-T0REF)*FRSDEN
      FRSSIE = -E0REF + RGAS/GM1*(FRSTEM-T0REF)
      SIEINI = -E0REF + RGAS/GM1*(TEMINI-T0REF)
      T0     = FRSTEM*(1.+.5*GM1*RMACH**2)
      FRSVIS = RMUFRS

      ENDIF

      END SUBROUTINE BCINI1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BCINI2(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &                  SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &                  GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,REFVEL)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON WATER MODEL
! ... 10.1.96 KPe

      IMPLICIT NONE

      REAL :: FRSVEL,REFVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &        SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &        GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMUFRS,CV,RMULT
      REAL :: A,B,C,D,CP,ERO
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

! ... Kertoimet viskositeetin laskemiseksi
!      A   = .291648393e+1
!      B   = .259588065e-2
!      C   = .958787090e+2

      A   =  .2414e-4            ! ESa 24.8.2001
      B   =  .2478e+3
      C   =  .1400e+3
      D   =  .113902e-5
      ERO = -.485e+1             ! ennen seuraava??? 15-FRSTEM+273.15
      CP  =  .418e+4
      CV  =  .418e+4             ! EI NYT IHAN OIKEA ARVO

      FRSSSP = .14848e+4
!      FRSSIE = .8377e+5
!      SIEINI = .8377e+5
      FRSSIE = CV*(FRSTEM - .293e+3) + .8377e+5
      SIEINI = CV*(FRSTEM - .293e+3) + .8377e+5
!      RMUFRS = EXP(((A*ERO)-(B*ERO**2))/(C+FRSTEM-273.15))*D*FRSDEN
	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta FRSTEM
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(FRSTEM-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*FRSTEM+.503834e-4)*FRSTEM
     &  -.233632e-1        +.369516e+1/FRSTEM)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - FRSTEM)))
	G22 = 1.0 - G12

	RMUFRS  = RMULT * (  G12 * G11 + G22 * G21 )

	! ---------------------------------
      FRSVIS = RMUFRS

      IF (FRSVEL == 0. .AND. RMACH /= 0.) THEN
       FRSVEL = RMACH*FRSSSP
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       T0     = FRSTEM+FRSVEL**2*0.5/CP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given RMACH,"
      ELSE IF (FRSVEL /= 0. .AND. RMACH == 0.) THEN
       RMACH = FRSVEL/FRSSSP
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       T0     = FRSTEM+FRSVEL**2*0.5/CP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given FRSVEL,"
      ELSE IF (FRSVEL == 0. .AND. RMACH == 0. .AND. REFVEL /= 0.) THEN
       RMACH = REFVEL/FRSSSP
       RE     = FRSDEN*REFVEL*CHLREF/RMUFRS
       T0     = FRSTEM+REFVEL**2*0.5/CP      
       WRITE(4,"(2X,A)") "Free-stream quantities based on given REFVEL,"
      ELSE
       PRINT"(A/A,2F11.5/A)",
     &    "Problem may be overdefined. Only one velocity scale may",
     &    "be specified. RMACH and FRSVEL are ",RMACH,FRSVEL,
     &    "Exiting..." ; STOP
      END IF

      WRITE(4,"(A)") "  FRSTEM, FRSPRE AND CHLREF"

      END SUBROUTINE BCINI2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BCINI3(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &              SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &              GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,REFVEL)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON WATER MODEL
! ... 10.1.96 KPe
! ... Can be used to modify the viscosity and Reynolds numbes as ISTATE = 7

      IMPLICIT NONE

      REAL :: FRSVEL,REFVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &        SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &        GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,RMUFRS,CV
      REAL :: A,B,C,D,CP,ERO
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

! ... Kertoimet viskositeetin laskemiseksi
!      A   = .291648393e+1
!      B   = .259588065e-2
!      C   = .958787090e+2

!      RMULT = 31.558*SQRT(31.558) !1250. !5640.7 !31.558*SQRT(31.558)
!      RMULT = 2000.

      A   =  .2414e-4            ! ESa 24.8.2001
      B   =  .2478e+3
      C   =  .1400e+3
      D   =  .113902e-5
      ERO = -.485e+1             ! ennen seuraava??? 15-FRSTEM+273.15
      CP  =  .418e+4
      CV  =  .418e+4             ! EI NYT IHAN OIKEA ARVO

      FRSSSP = .14848e+4
!      FRSSIE = .8377e+5
!      SIEINI = .8377e+5
      FRSSIE = CV*(FRSTEM - .293e+3) + .8377e+5
      SIEINI = CV*(FRSTEM - .293e+3) + .8377e+5
	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta FRSTEM
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(FRSTEM-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*FRSTEM+.503834e-4)*FRSTEM
     &  -.233632e-1        +.369516e+1/FRSTEM)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - FRSTEM)))
	G22 = 1.0 - G12

	RMUFRS  = RMULT*( G12 * G11 + G22 * G21)

	! ---------------------------------
      FRSVIS = RMUFRS

      IF (FRSVEL == 0. .AND. RMACH /= 0.) THEN
       FRSVEL = RMACH*FRSSSP
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       T0     = FRSTEM+FRSVEL**2*0.5/CP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given RMACH,"
      ELSE IF (FRSVEL >= 1.E-6 .AND. RMACH == 0.) THEN
       RMACH = FRSVEL/FRSSSP
       RE     = FRSDEN*FRSVEL*CHLREF/RMUFRS
       T0     = FRSTEM+FRSVEL**2*0.5/CP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given FRSVEL,"
      ELSE
       PRINT"(A/A,2F11.5/A)",
     &    "Problem maybe overdefined. Only one velocity scale may",
     &    "be specified. RMACH and FRSVEL are ",RMACH,FRSVEL!,
!     &    "Exiting..." ; STOP
      RE     = FRSDEN*REFVEL*CHLREF/RMUFRS
      T0     = FRSTEM+REFVEL**2*0.5/CP
      END IF
      FRSDEN = 1.e+3

      WRITE(4,"(A)") "  FRSTEM, FRSPRE AND CHLREF"

      END SUBROUTINE BCINI3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BCINI4(FRSVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,
     &             SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,FRSCP,REFVEL)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON WATER OR STEAM MODEL

      IMPLICIT NONE

      REAL :: FRSVEL,REFVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,SIEINI,
     &        FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,FRSCP

      SIEINI = FRSSIE
      IF (FRSVEL == 0. .AND. RMACH /= 0.) THEN
       FRSVEL = RMACH*FRSSSP
       RE = FRSDEN*FRSVEL*CHLREF/FRSVIS
       T0 = FRSTEM+FRSVEL**2*0.5/FRSCP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given RMACH,"
      ELSE IF (FRSVEL /= 0. .AND. RMACH == 0.) THEN
       RMACH = FRSVEL/FRSSSP
       RE = FRSDEN*FRSVEL*CHLREF/FRSVIS
       T0 = FRSTEM+FRSVEL**2*0.5/FRSCP
       WRITE(4,"(2X,A)") "Free-stream quantities based on given FRSVEL,"
      ELSE
       PRINT"(A/A,2F11.5/A)",
     &    "Problem may be overdefined. Only one velocity scale may",
     &    "be specified. RMACH and FRSVEL are ",RMACH,FRSVEL!,
!     &    "Exiting..." ; STOP
       RE = FRSDEN*REFVEL*CHLREF/FRSVIS
       T0 = FRSTEM+REFVEL**2*0.5/FRSCP
      END IF

      WRITE(4,"(A)") "  FRSTEM, FRSPRE AND CHLREF"

      END SUBROUTINE BCINI4
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BCINI5(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &           SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &           GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,PSEUCO,REFVEL)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON INCOMPRESSIBLE MODEL

      IMPLICIT NONE

      REAL :: FRSVEL,REFVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &        SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,
     &        VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,RMUFRS,GM1,PSEUCO
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

      FRSSIE = 1.E 5
      SIEINI = FRSSIE
      FRSSSP = PSEUCO*REFVEL !FRSVEL
	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta FRSTEM
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(FRSTEM-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*FRSTEM+.503834e-4)*FRSTEM
     &  -.233632e-1        +.369516e+1/FRSTEM)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - FRSTEM)))
	G22 = 1.0 - G12

	FRSVIS  =  RMULT *( G12 * G11 + G22 * G21)

	! ---------------------------------
      T0     = FRSTEM
      IF(ABS(REFVEL) < 1.E-20) REFVEL = FRSVEL

      IF(RE /= 0.) THEN
          FRSVIS = REFVEL*FRSDEN*CHLREF/RE
          WRITE(4,'(2X,2A)') 'Free-stream viscosity was based on the',
     &     'Reynolds number.'
      ELSE
          RE     = REFVEL*FRSDEN*CHLREF/FRSVIS
          WRITE(4,"(2X,2A)") 'Free-stream Reynolds number was based on',
     &  ' the given viscosity.'
      ENDIF

      END SUBROUTINE BCINI5
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BCINI6(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,
     &              FRSSIE,SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,
     &           GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,PSEUCO,REFVEL)

! ... FREE STREAM CONDITIONS ARE EVALUATED. BASED ON INCOMPRESSIBLE MODEL

      IMPLICIT NONE

      REAL :: FRSVEL,REFVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,
     &        SIEINI,FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,
     &        VISU0,EXPSU,TSU0,E0REF,T0REF,RMULT,RMUFRS,GM1,PSEUCO

      FRSSIE = 6.E4
      SIEINI = FRSSIE
      FRSSSP = PSEUCO*REFVEL  ! FRSVEL
c      FRSVIS = RMULT*1.E-5
      FRSVIS = ((0.659E-3*(FRSTEM-273.15-1.0)-0.05076)*
     &     (FRSTEM-273.15-1.0)+1.7688)*1E-3
      RE     = REFVEL*FRSDEN*CHLREF/FRSVIS
         WRITE(4,'(2X,2A)') 'Sea water at constant density'
         WRITE(4,"(2X,2A)") 'Free-stream Reynolds number was based on',
     &        ' the given viscosity.'

      END SUBROUTINE BCINI6
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

! **********************************************************************
! ... old equation of state
! **********************************************************************

      SUBROUTINE PERFG2(VIS,P,RO,E,TEMP,CH,DRDP,DRDH,RGAS,GAMMA,
     &                  PR,VIS0,EXPSU,T0,NTOT,ERCODE,E0REF,T0REF)

! ... Equation of State for a Perfect Gas, 17.1.1994, PMK
! ... INPUT  : RO,E
! ... OUTPUT : P,TEMP,DRDP,DRDH,CH,VIS

      IMPLICIT NONE

      REAL,DIMENSION(*) :: VIS,P,RO,E,TEMP,CH,DRDP,DRDH
      REAL :: RGAS,GAMMA,PR,VIS0,EXPSU,T0,E0REF,T0REF,GM1,GM2,CP
      INTEGER :: NTOT,ERCODE

      GM1 = GAMMA - 1.
      GM2 = GAMMA*GM1
      CP  = GAMMA/GM1*RGAS

      IF (ABS(EXPSU-1.5) > .1e-4) THEN
       P(1:NTOT)    = GM1*(E(1:NTOT)+E0REF)*RO(1:NTOT)
       TEMP(1:NTOT) = T0REF + P(1:NTOT)/(RO(1:NTOT)*RGAS)
       DRDP(1:NTOT) = 1./(GM1*E(1:NTOT))
       DRDH(1:NTOT) = -P(1:NTOT)/(GM2*E(1:NTOT)**2)
! ... Sutherland's Formula
       VIS(1:NTOT)  = VIS0*TEMP(1:NTOT)**EXPSU/(TEMP(1:NTOT)+T0)
       CH(1:NTOT)   = VIS(1:NTOT)*CP/PR
      ELSE
       P(1:NTOT)    = GM1*(E(1:NTOT)+E0REF)*RO(1:NTOT)
       TEMP(1:NTOT) = T0REF + P(1:NTOT)/(RO(1:NTOT)*RGAS)
       DRDP(1:NTOT) = 1./(GM1*E(1:NTOT))
       DRDH(1:NTOT) = -P(1:NTOT)/(GM2*E(1:NTOT)**2)
! ... Sutherland's Formula
       VIS(1:NTOT)  = VIS0*TEMP(1:NTOT)**1.5/(TEMP(1:NTOT)+T0)
       CH(1:NTOT)   = VIS(1:NTOT)*CP/PR
      END IF

      END SUBROUTINE PERFG2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

! **********************************************************************
! ... Multiphase routines
! **********************************************************************

      SUBROUTINE STATEM(PRO,P,PDIFF,GAMMA,PR,PRT,NTOT,ITURB,ISTATE,
     &                  FRSPRE,RGAS,VISU0,EXPSU,TSU0,E0REF,T0REF,
     &                  IMAX,JMAX,KMAX,INN,IPHASE,FRSDEN,RMULT)

! ... Calling Routine for the Equation of State
! ... INPUT  : P,T
! ... OUTPUT : P,TEMP,DRDP,DRDH,CH,VIS or RO,TEMP,CH,VIS

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: INN,IPHASE,NTOT,ITURB,ISTATE,IMAX,JMAX,KMAX
      INTEGER :: ERCODE
      REAL,DIMENSION(*) :: P,PDIFF
      REAL :: GAMMA,PR,PRT,FRSPRE,RGAS,VISU0,EXPSU,TSU0,E0REF,T0REF,
     &        FRSDEN,RMULT
      TYPE(PROPERTIES),DIMENSION(*) :: PRO

      ERCODE = 0

      SELECT CASE(ISTATE)

! **********************************************************************
! ... Perfect Gas
! **********************************************************************
      CASE(1)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                "Select ISTATE = 3 for ideal gas. Exiting..."
       STOP

! **********************************************************************
! ... Chemically Reacting Equilibrium Air
! **********************************************************************
      CASE(2)
!       CALL TGAS(P,RO,E,H,TEMP,A,S,DRDP,DRDH, 2,NTOT,ERCODE)
!       CALL TGAS(P,RO,E,H,TEMP,A,S,DRDP,DRDH,11,NTOT,ERCODE)
!       CALL UGAS3(E,RO,VIS,NTOT,ERCODE)
!       CALL UGAS4(E,RO,CH,NTOT,ERCODE)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                "Exiting..."
       STOP

! **********************************************************************
! ... Perfect Gas with cp(T)
! **********************************************************************
      CASE(3)
       CALL PERFGAM(PRO,P,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE,
     &              RGAS,GAMMA,PR,VISU0,EXPSU,TSU0,ERCODE,E0REF,T0REF)

! **********************************************************************
! ... Semi perfect Gas
! **********************************************************************
      CASE(4)
!       CALL PERSEM(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
!     &             VISU0,EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                "Select ISTATE = 3 for ideal gas. Exiting..."
       STOP

! **********************************************************************
! ... Air modeled as a perfect gas (GAMMA and CP are not constants)
! **********************************************************************
      CASE(5)
!       CALL PERAIR(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,RGAS,GAMMA,PR,
!     &             VISU0,EXPSU,TSU0,NTOT,ERCODE,E0REF,T0REF)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                "Select ISTATE = 3 for ideal gas. Exiting..."
       STOP

! **********************************************************************
! ... Water (old version)
! **********************************************************************
      CASE(6)
       CALL WATERM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE,
     &        RMULT,FRSDEN)

! **********************************************************************
! ... Water with modified viscosity (modified Reynolds number)
! **********************************************************************
      CASE(7)
!       CALL PWATER(VIS,P,RO,E,TEMP,CP,CH,DRDP,DRDH,PR,NTOT,
!     &             IMAX,JMAX,KMAX,IN,JN,KN,INN)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                "Modify viscosity directly. Exiting..."
       STOP

! **********************************************************************
! ... Water as a function of pressure and temperature
! **********************************************************************
      CASE(8)
       CALL WATER2M(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE,
     & FRSDEN)

! **********************************************************************
! ... Steam as function of pressure and temperature
! **********************************************************************
      CASE(9)
       CALL STEAM3M(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE)

      CASE(10:11)
       PRINT"(/A/A/)","No such E.o.S. for a multiphase flow.",
     &                " Exiting..."
       STOP
! **********************************************************************
! ... Liquid R-12 freon as a function of pressure and temperature
! **********************************************************************
      CASE(12)
       CALL R12LIQM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE)

! **********************************************************************
! ... Gas R-12 freon as function of pressure and temperature
! **********************************************************************
      CASE(13)
       CALL R12GASM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE)

! **********************************************************************
! ... Routines for solids
! **********************************************************************
! ... Cast iron
! **********************************************************************
      CASE(101)
       CALL CAST_IRONM(PRO,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE)

      END SELECT

      IF (ERCODE > 0) WRITE(4,"(/9X,I10,A)") ERCODE,
     &    " OVERFLOW ERRORS occured in STATEM subroutine!"

      END SUBROUTINE STATEM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CAST_IRONM(PRO,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE)

! ... Input : T = temperature, [K]

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN,IPHASE
      INTEGER :: I,II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

! ... specific heat of cast iron, [J/(kgK)]
      PRO(I)%CP(IPHASE) = .4200e+3

! ... density, [kg/m^3]
      PRO(I)%RO(IPHASE) = .7272e+4

! ... internal energy, [J/(kgK)]  (oli [J/(m^3K)])
      PRO(I)%E(IPHASE)  = PRO(I)%CP(IPHASE)*PRO(I)%TEMP(IPHASE)

! ... thermal conducitivity, [W/Km]
      PRO(I)%CH(IPHASE) = .51e+2 - .40e-1*PRO(I)%TEMP(IPHASE)

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE CAST_IRONM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERFGAM(PRO,P,NTOT,IMAX,JMAX,KMAX,IN,JN,KN,INN,IPHASE,
     &                   RGAS,GAMMA,PR,VIS0,EXPSU,T0,ERCODE,E0REF,T0REF)

! ... Equation of state for perfect gas

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: PR,RGAS,GAMMA,VIS0,EXPSU,T0,E0REF,T0REF
      REAL :: A,B,C,D,T,TVIS,GM1,CV,CP0,CP1,CP2,CP3,PROTEMP,PROTGAS
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN,IPHASE,ERCODE
      INTEGER :: I,II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL

      GM1  = GAMMA - 1.
      CP2  = GAMMA/GM1*RGAS

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

         PROTEMP = PRO(I)%TEMP(IPHASE)
         PROTGAS = PROTEMP*RGAS

         PRO(I)%RO(IPHASE)    = P(I)/PROTGAS
         CP0 =    .124500e+4
         CP1 = ((-.108980e-6 *PROTEMP + .333093e-3)*PROTEMP
     &           -.118565e+0)*PROTEMP + .101582e+4
         CP3 =    .296940e+0 *PROTEMP + .703060e+3

         IF (CP0 < MAX(CP1,CP3)) THEN
          PRO(I)%E(IPHASE)    =
     &            .124500e+4 *PROTEMP
         ELSE IF (CP1 > CP3) THEN
          PRO(I)%E(IPHASE)    =
     &        (((-.272450e-7 *PROTEMP + .111031e-3)*PROTEMP
     &           -.592825e-1)*PROTEMP + .101582e+4)*PROTEMP
         ELSE
          PRO(I)%E(IPHASE)    =
     &          ( .148470e+0 *PROTEMP + .703060e+3)*PROTEMP
         END IF
         PRO(I)%E(IPHASE)     = PRO(I)%E(IPHASE) - PROTGAS

         PRO(I)%DRODP(IPHASE) = 1./PROTGAS
         PRO(I)%DRODH(IPHASE) = -P(I)*GM1/GAMMA*PRO(I)%DRODP(IPHASE)**2
         PRO(I)%VIS(IPHASE)   = VIS0*PROTEMP*SQRT(PROTEMP)/(PROTEMP+T0)
         PRO(I)%CP(IPHASE)    = MIN(MAX(CP1,CP3),CP0)
         PRO(I)%CH(IPHASE)    = PRO(I)%VIS(IPHASE)*CP2/PR

! ... Artificial values for saturation
         PRO(I)%HSAT(IPHASE)  = PRO(I)%E(IPHASE)
         PRO(I)%DHSDP(IPHASE) = 0.

      END DO ; END DO ; END DO

      END SUBROUTINE PERFGAM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WATERM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,
     &                  IN,JN,KN,INN,IPHASE,RMULT,FRSDEN)

! ... Equation of State for water; 8.1.96 KPe
! ... INPUT  : P,TEMP
! ... OUTPUT : RO,E,DRDP,DRDH,CH,VIS

      USE TYPE_ARRAYS
      
      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: PR,RMULT
      REAL :: A,B,C,D,CP2,CV,T,TVIS,PROTEMP,FRSDEN
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN,IPHASE
      INTEGER :: I,II,JJ,KK,JA,KA,ISTRID,JSTRID,INN1,IL
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

!      A    = .291648393e+1
!      B    = .259588065e-2
!      C    = .958787090e+2
      A    =  .2414e-4          ! ESa 24.8.2001
      B    =  .2478e+3
      C    =  .1400e+3
      D    =  .113902e-5
      CP2  =  .418e+4
      CV   =  .418e+4    ! EI NYT IHAN OIKEA ARVO
      T    = -.485e+1

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

         PROTEMP = PRO(I)%TEMP(IPHASE)

         PRO(I)%RO(IPHASE)    = FRSDEN !1025.8 !????????
         PRO(I)%E(IPHASE)     =  .8377e+5 + (PROTEMP-.2930e+3)*CV
         PRO(I)%DRODP(IPHASE) =  .5000e-6
         PRO(I)%DRODH(IPHASE) = -.4970e-4
         TVIS                 = MAX(.27315e+3,MIN(.37000e+3,PROTEMP))
         TVIS                 = PROTEMP
	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta TVIS
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(TVIS-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*TVIS+.503834e-4)*TVIS
     &  -.233632e-1        +.369516e+1/TVIS)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - TVIS)))
	G22 = 1.0 - G12

!!!!	PRO(I)%VIS(IPHASE)  =  G12 * G11 + G22 * G21
	PRO(I)%VIS(IPHASE)  = RMULT * (  G12 * G11 + G22 * G21 )

	! ---------------------------------



         PRO(I)%CP(IPHASE)    = CP2
         PRO(I)%CH(IPHASE)    =  .556e+0 + .125e-2*(PROTEMP-.275e+3)

! ...  hl_sat(T) from WATERM2 (mersu)
         PRO(I)%HSAT(IPHASE)  =  .127896e+4*SQRT(P(I))
     &      +( .323843e-6*P(I) - .887634e+0)*P(I) + .980548e+5
         PRO(I)%DHSDP(IPHASE) =  .639480e+3/SQRT(P(I))
     &      +  .647686e-6*P(I) - .887634e+0

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE WATERM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE STEAM3M(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,INN,IPHASE)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: PR
      REAL :: YC2,CG,CL,DRDT,DRDPA,CP0,CP1,CP5,F,F1,F5,
     &        A1,A2,AA,RLAM1,RLAM2,RLAM5,PROTEMP,PRORO,PROCP
      REAL :: F00,F11,F12,F21,F22,F31,F32,F41,F42,F51,F52,F61,F62
      REAL :: R1,R2,R3,R4,R5,R6
      REAL :: ROO1,ROO10,ROO100,ROO300,ROO2M
      REAL :: DRDT1,DRDT10,DRDT100,DRDT300,DRDT2M
      REAL :: R2P,R3P,R4P,R5P,R6P
      REAL :: SLP1,SLP2,LORE
      INTEGER :: NTOT,IMAX,JMAX,KMAX,INN,IN,JN,KN,IPHASE
      INTEGER :: I,INN1,II,JJ,KK,JA,KA,ISTRID,JSTRID,IL
      REAL, PARAMETER ::
     &    V1 =  0.249899160714079e+2,
     &    W1 = -0.504070497017732e+4

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

c      IF (MULPHC == "CAVIT") THEN

c      PRO(I)%TEMP(2) = W1/(LOG(P(I)) - V1)

      PROTEMP = PRO(I)%DTEMP(IPHASE)

! ... specific heat of gas, [J/(kgK)]

      F00     = 1e-5*MAX(10.,P(I))
      F12     = 1. - F00
      F21     =      F00
      F00     = .25*(F00 - 1.)
      F22     = 1  - F00
      F31     =      F00
      F32     = 1.
      CP0     = ( .944056e-3*PROTEMP - .294568e+0)*PROTEMP + .186841e+4
      CP1     = ( .113393e-1*PROTEMP - .102209e+2)*PROTEMP + .427457e+4
      CP5     = ( .362500e-1*PROTEMP - .360791e+2)*PROTEMP + .110935e+5
      PRO(I)%CP(IPHASE)   = MAX(0. ,F12) *CP0
     &             + MAX(0.,MIN(F21,F22))*CP1
     &             + MAX(0.,MIN(F31,F32))*CP5

! ... density of gas, [kg/m^3] and its derivatives

! ... Start of new addition by A. Miettinen, 2007/10/01

! ... 1kPa
      SLP1 = 0.001
      SLP2 = 1.0/9.0e+3
      F21  =            SLP1*MAX(10.,P(I))
      F00  = 0.111111 - SLP2*MAX(10.,P(I))
      F22  = 1.0 + F00
      R2   = MAX(0.,MIN(F21,F22))
!      ROO1    = -.134939e-4*PROTEMP + .111220e-1
      ROO1    = ( .322363e-7*PROTEMP - .406131e-4)*PROTEMP + .165018e-1
      DRDT1   =   .644726e-7*PROTEMP - .406131e-4
      R2P     = LORE(R2 > 0.0)*(SLP1*LORE(F21 < F22)
     &                        - SLP2*LORE(F21 > F22))

! ... 10kPa
      SLP1 = SLP2
      SLP2 = 1.0/9.0e+4
      F31  = -F00
      F00  = 0.111111 - SLP2*P(I)
      F32  = 1.0 + F00
      R3   = MAX(0.,MIN(F31,F32))
!      ROO10   = -.119482e-3*PROTEMP + .103647e+0
      ROO10   = ( .281416e-6*PROTEMP - .367450e-3)*PROTEMP + .156294e+0
      DRDT10  =   .562832e-6*PROTEMP - .367450e-3
      R3P     = LORE(R3 > 0.0)*(SLP1*LORE(F31 < F32)
     %                        - SLP2*LORE(F31 > F32))

! ... 100kPa
      SLP1 = SLP2
      SLP2 = 1.0/2.0e+5
      F41  = -F00
      F00  = 0.500000 - SLP2*P(I)
      F42  = 1.0 + F00
      R4   = MAX(0.,MIN(F41,F42))
!      ROO100  = -.105576e-2*PROTEMP + .972238e+0
      ROO100  = ( .240937e-5*PROTEMP - .331947e-2)*PROTEMP + .149104e+1
      DRDT100 =   .481874e-5*PROTEMP - .331947e-2
      R4P     = LORE(R4 > 0.0)*(SLP1*LORE(F41 < F42)
     &                        - SLP2*LORE(F41 > F42))

! ... 300kPa
      SLP1 = SLP2
      SLP2 = 1.0/1.7e+6
      F51  = -F00
      F00  = 0.176471 - SLP2*P(I)
      F52  = 1.0 + F00
      R5   = MAX(0.,MIN(F51,F52))
!      ROO300  = -.304522e-2*PROTEMP + .286537e+1
      ROO300  = ( .726955e-5*PROTEMP - .101508e-1)*PROTEMP + .457358e+1
      DRDT300 =   .145391e-4*PROTEMP - .101508e-1
      R5P     = LORE(R5 > 0.0)*(SLP1*LORE(F51 < F52)
     &                        - SLP2*LORE(F51 > F52))

! ... 2Mpa
      SLP1 = SLP2
      SLP2 = 0.0
      F61  = -F00
      F62  = 1.0
      R6   = MAX(0.,MIN(F61,F62))
!      ROO2M   = -.151110e-1*PROTEMP + .169803e+2
      ROO2M   = ( .348933e-4*PROTEMP - .585389e-1)*PROTEMP + .301263e+2
      DRDT2M  =   .697866e-4*PROTEMP - .585389e-1
      R6P     = LORE(R6 > 0.0)*(SLP1*LORE(F61 < F62)
     &                        - SLP2*LORE(F61 > F62))

      PRO(I)%RO(IPHASE)   =
     &        R2 *ROO1 +R3 *ROO10 +R4 *ROO100 +R5 *ROO300 +R6 *ROO2M

      DRDT  = R2 *DRDT1+R3 *DRDT10+R4 *DRDT100+R5 *DRDT300+R6 *DRDT2M

      DRDPA = R2P*ROO1 +R3P*ROO10 +R4P*ROO100 +R5P*ROO300 +R6P*ROO2M

!      PRO(I)%RO(IPHASE)   =
!     &          ( .413508e-5*PROTEMP - .395082e-2)*PROTEMP + .901454e+0
!     &         +(-.104028e-7*PROTEMP + .976698e-5)*P(I)

!      DRDT    = ( .827016e-5*PROTEMP - .395082e-2) -.104028e-7*P(I)
!      DRDPA   =  -.104028e-7*PROTEMP + .976698e-5

      PRORO   = PRO(I)%RO(IPHASE)
      PROCP   = PRO(I)%CP(IPHASE)

      PRO(I)%DRODH(IPHASE) = DRDT/PROCP
      PRO(I)%DRODP(IPHASE) = DRDPA
     &      - PRO(I)%DRODH(IPHASE)*(1.+PROTEMP*DRDT/PRORO)/PRORO

      YC2   = PRO(I)%DRODH(IPHASE)/PRORO+PRO(I)%DRODP(IPHASE)
      CG    = 1./SQRT(YC2)

      END DO ; END DO ; END DO

!      END DO

!      PR      = 0. ! An average Prandtl number could be evaluated blockwise

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

      PROTEMP = PRO(I)%DTEMP(IPHASE)

! ... dynamic viscosity of gas, [kg/(sm)]

      PRO(I)%VIS(IPHASE) =
     &          ( .016000e-9*PROTEMP + .254589e-7)*PROTEMP + .551790e-6

! ... thermal conductivity of gas, [W/Km]

      F00     = .25*(1. - 1e-5*P(I))
      F11     = 1.
      F12     = 1. + F00
      F21     =    - F00
      F22     = 1.
      F1      = MAX(0.,MIN(F11,F12))
      F5      = MAX(0.,MIN(F21,F22))
      RLAM1   = ( .132000e-6*PROTEMP - .297141e-4)*PROTEMP + .177873e-1
      RLAM5   = ( .832596e-7*PROTEMP + .204684e-5)*PROTEMP + .153213e-1
      PRO(I)%CH(IPHASE) = F1*RLAM1 + F5*RLAM5

! ... Average Prandtl number for this loop

!      PR      = ((NTOT-1.)*PR + PRO(I)%VIS(IPHASE)*PRO(I)%CP(IPHASE)/
!     &           PRO(I)%CH(IPHASE))/NTOT

! ... internal energy of gas, [J/(kgK)]

c      PRO(I)%E(IPHASE) =
c     &          (-.376298e-1*PROTEMP + .150939e+4)*PROTEMP + .195633e+7
c     &         +( .360841e-3*PROTEMP - .212573e+0)*P(I)

      PRO(I)%HSAT(IPHASE)  =  .543794e+3*SQRT(MAX(10.,P(I))) ! Obs
     &   +( .149849e-6*P(I) - .437937e+0)*P(I) + .254522e+7

      PRO(I)%E(IPHASE) = PRO(I)%HSAT(IPHASE)-P(I)/PRO(I)%RO(IPHASE) +
     & PRO(I)%CP(IPHASE)*(PRO(I)%DTEMP(IPHASE)-PRO(I)%TSAT)

      PRO(I)%DHSDP(IPHASE) =  .271897e+3/SQRT(P(I))
     &   +  .299698e-6*P(I) - .437937e+0

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE STEAM3M
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WATER2M(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,INN,IPHASE,FRSDEN)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR,
     &        PTMP,PICP,PIRO,PIDRODH,PIDRODP,F00,FRSDEN
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IPHASE,IL
      REAL :: G11, G21, G12, G22		! Viskositeetin kaava -MPa

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... specific heat of liquid, [J/(kgK)]

      PICP    = MAX(.4182e+4,
     &          ( .560839e-2*PTMP-.284106e+1)*PTMP+.450010e+4)

! ... density of liquid, [kg/m^3]

      PIRO    = (-.260031e-2*PTMP+.122284e+1)*PTMP+.864232e+3
     &        + ( .761625e-9*PTMP+.196854e-6)*P(I)

      DRDT    =  -.520062e-2*PTMP+.122284e+1      +.761625e-9*P(I)
      DRDPA   =   .761625e-9*PTMP+.196854e-6

      F00     = DRDT/PIRO
      PIDRODH = DRDT/PICP
      PIDRODP = DRDPA - F00*(1.+ F00*PTMP)/PICP

      YC2 = PIDRODH/PIRO+PIDRODP
      CL  = 1./SQRT(YC2)

      PRO(I)%CP(IPHASE) = PICP
      PRO(I)%RO(IPHASE) = PIRO
      PRO(I)%DRODH(IPHASE) = PIDRODH
      PRO(I)%DRODP(IPHASE) = PIDRODP

      END DO ; END DO ; END DO

!      END DO

!      PR      = 0. ! An average Prandtl number could be evaluated blockwise

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... dynamic viscosity of liquid, [kg/(sm)]

!      PRO(I)%VIS(IPHASE) = MAX(1e-4,
!     &  (-.363718e-7*PTMP+.503834e-4)*PTMP-.233632e-1+.369516e+1/PTMP)

	! ---- dynamic viscosity of water ----
	! dynaaminen viskositeetti riippuu tassa muuttujasta PTMP
	! Korjannut MPa 16.12.2016. Esan kaavaan (=G11) on sijoitettu arvot
	! A   =  .2414e-4
	! B   =  .2478e+3
	! C   =  .1400e+3

	G11 = .2414e-4 * 10.**(.2478e+3/(PTMP-.1400e+3))
	G21 = MAX(.1e-3,(-.363718e-7*PTMP+.503834e-4)*PTMP
     &  -.233632e-1        +.369516e+1/PTMP)
	G12 = MAX(0.,MIN(1.0,0.1 * (310. - PTMP)))
	G22 = 1 - G12

	PRO(I)%VIS(IPHASE) =  G12 * G11 + G22 * G21
	! ---------------------------------

! ... thermal conductivity of liquid, [W/Km]

      F11     =   .385e+2-.100e+0*PTMP
      F21     =  1. - F11
      F12     = MAX(0.,MIN(1.,F11))
      F22     = MAX(0.,MIN(1.,F21))
      RLAM1   = (-.921140e-5*PTMP+.713593e-2)*PTMP-.700822e+0
      RLAM2   =   .680850e+0
      PRO(I)%CH(IPHASE) = MAX(.56e+0,F12*RLAM1+F22*RLAM2)

! ... Average Prandtl number for this loop

!      PR = ((NTOT-1)*PR+PRO(I)%VIS(IPHASE)*CP(I)/PRO(I)%CH(IPHASE))/NTOT

! ... internal energy of liquid, [J/(kgK)]

      PRO(I)%E(IPHASE) = ( .386586e+0*PTMP+.392553e+4)*PTMP-.109967e+7
     &                 + ( .447332e-5*PTMP-.169196e-2)*P(I)

! ...  hl_sat(T)

      PRO(I)%HSAT(IPHASE)  =  .127896e+4*SQRT(P(I))
     &       +( .323843e-6*P(I) - .887634e+0)*P(I) + .980548e+5
      PRO(I)%DHSDP(IPHASE) =  .639480e+3/SQRT(P(I))
     &       +  .647686e-6*P(I) - .887634e+0

! ... Surface tension of liquid

      PRO(I)%SIGMA = 6.088122774539E-7*PTMP**3 -0.00091967773504*PTMP**2
     &             + 0.233461051381321*PTMP + 68.1605213077499
      PRO(I)%SIGMA = MIN(79.01,MAX(PRO(I)%SIGMA,0.1)) * 0.001

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE WATER2M
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE R12LIQM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,INN,IPHASE)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR,
     &        PTMP,PICP,PIRO,PIDRODH,PIDRODP,F00
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IPHASE,IL

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... specific heat of liquid R-12, [J/(kgK)]


      PRO(I)%CP(IPHASE) = 1419.8
      PRO(I)%RO(IPHASE) = 1017.1
      PRO(I)%DRODH(IPHASE) = -1.E-4
      PRO(I)%DRODP(IPHASE) = .5E-6

      END DO ; END DO ; END DO

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II


      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... dynamic viscosity of liquid, [kg/(sm)]

      PRO(I)%VIS(IPHASE) = 8.95E-5

! ... thermal conductivity of liquid, [W/Km]

      PRO(I)%CH(IPHASE) = 4.57E-2

! ... internal energy of liquid, [J/(kgK)]

c      PRO(I)%E(IPHASE) = 1.2693E5 + PRO(I)%CP(IPHASE)*(PTMP - 359.79)
c     +                 - P(I)/PRO(I)%RO(IPHASE)
      PRO(I)%E(IPHASE) = PRO(I)%CP(IPHASE)*(PTMP - 298.15)
     +                 - P(I)/PRO(I)%RO(IPHASE)
      PRO(I)%E(IPHASE) = MAX(PRO(I)%E(IPHASE),1.)


! ...  hl_sat(T)

      PRO(I)%HSAT(IPHASE)  =  PRO(I)%CP(IPHASE)*(359.79 - 298.15)
      PRO(I)%DHSDP(IPHASE) =  1. ! Oikeasti nolla

! ... Surface tension

      PRO(I)%SIGMA =  0.016 ! 

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE R12LIQM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE R12GASM(PRO,P,PR,NTOT,IMAX,JMAX,KMAX,
     &                   IN,JN,KN,INN,IPHASE)

! ... Input : T = temperature, [K]
! ...         P = pressure, [Pa]

      USE TYPE_ARRAYS

      IMPLICIT NONE

      TYPE(PROPERTIES),DIMENSION(*) :: PRO
      REAL,DIMENSION(*) :: P
      REAL :: YC2,CG,CL,DRDT,DRDPA,F12,F21,F22,F31,F32,CP0,CP1,CP5,
     &        F,A1,A2,AA,F1,F5,F11,RLAM1,RLAM2,RLAM5,PR,
     &        PTMP,PICP,PIRO,PIDRODH,PIDRODP,F00
      INTEGER :: NTOT,I,IMAX,JMAX,KMAX,INN,IN,JN,KN,II,JJ,KK,JA,KA,
     &           ISTRID,JSTRID,INN1,IPHASE,IL

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      INN1   = INN-1

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

!      DO I = 1,NTOT

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... specific heat of gas R-12, [J/(kgK)]


      PRO(I)%CP(IPHASE) = 1290.8
      PRO(I)%RO(IPHASE) = 172.06
      PRO(I)%DRODH(IPHASE) = -1.E-4
      PRO(I)%DRODP(IPHASE) = 1.E-6

      END DO ; END DO ; END DO

c      DO KK   = -INN,KMAX+INN1
c       KA     = (KN+KK)*JSTRID + JN
c       DO JJ  = -INN,JMAX+INN1
c        JA    = (KA+JJ)*ISTRID + IN
c        DO II = -INN1,IMAX+INN
c         I    =  JA+II

      DO KK =  1-INN,KMAX+INN
      KA      = (KN+KK-1)*IL
      DO JJ =  1-INN,JMAX+INN
      JA      = (JN+JJ-1)*ISTRID + IN + KA
      DO II =  1-INN,IMAX+INN
         I   = JA + II

      PTMP    = PRO(I)%DTEMP(IPHASE)

! ... dynamic viscosity of gas, [kg/(sm)]

      PRO(I)%VIS(IPHASE) = 1.572E-5

! ... thermal conductivity of gas, [W/Km]

      PRO(I)%CH(IPHASE) = 1.76E-2

! ... internal energy of gas, [J/(kgK)]

c      PRO(I)%E(IPHASE) = 2.13E5 + PRO(I)%CP(IPHASE)*(PTMP - 359.79)
c     +                 - P(I)/PRO(I)%RO(IPHASE)
c      PRO(I)%E(IPHASE) = MAX(PRO(I)%E(IPHASE),1.)

! ...  hl_sat(T)

      PRO(I)%HSAT(IPHASE)  = PRO(I)%CP(1)*(359.79 - 298.15) +
     +                       86070.
      PRO(I)%E(IPHASE) = PRO(I)%HSAT(IPHASE) + 
     + PRO(I)%CP(IPHASE)*(PTMP - 359.79)- P(I)/PRO(I)%RO(IPHASE)
      PRO(I)%E(IPHASE) = MAX(PRO(I)%E(IPHASE),1.)

      PRO(I)%DHSDP(IPHASE) =  1. ! Oikeasti nolla

      END DO ; END DO ; END DO

!      END DO

      END SUBROUTINE R12GASM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PSATF(PSAT,TEMP,ISTATE)

      IMPLICIT NONE

      REAL :: TEMP,PSAT
      INTEGER :: ERCODE,NTOT,ISTATE,I

! ... Saturation pressure as a function of temperature (ISTATE = 8)

      PSAT = MAX(1e-5,(( .305416e+0 *TEMP-.283410e+3)*TEMP
     &                  +.878017e+5)*TEMP-.907075e+7)

      END SUBROUTINE PSATF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DPSATF(PSAT,TEMP,FRSPRE,ISTATE)

      IMPLICIT NONE

      REAL :: TEMP,PSAT,FRSPRE
      INTEGER :: ERCODE,NTOT,ISTATE

! ... Saturation pressure as a function of temperature (ISTATE = 8)

      PSAT = MAX(1e-5,(( .305416e+0 *TEMP-.283410e+3)*TEMP
     &                  +.878017e+5)*TEMP-.907075e+7)-FRSPRE

      END SUBROUTINE DPSATF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DPSATLOG(PSAT,TEMP,ISTATE)

! ... Saturation pressure as a function of temperature (ISTATE = 8 / 18)

      IMPLICIT NONE

      INTEGER :: ISTATE
      REAL, PARAMETER ::
     &    P1 =  0.154906419868284e+6,
     &    Q1 = -0.101158884083608e+5,
     &    R1 =  0.133000059771258e+3,
     &    S1 =  0.146861326101013e-1,
     &    T1 = -0.170477189812508e+2
      REAL :: TEMP,PSAT

      PSAT = EXP((P1/TEMP + Q1)/TEMP + R1 + S1*TEMP + T1*LOG(TEMP))

      END SUBROUTINE DPSATLOG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DTSATLOG(TSAT,PRES,ISTATE)

! ... Saturation temperature as a function of pressure (ISTATE = 8 / 18)
! ... Used as a initial guess for subroutine FTEMPS

      IMPLICIT NONE

      INTEGER :: ISTATE
      REAL, PARAMETER ::
     &    V1 =  0.249899160714079e+2,
     &    W1 = -0.504070497017732e+4
      REAL :: PRES,TSAT

      TSAT = W1/(LOG(PRES) - V1)

      END SUBROUTINE DTSATLOG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FCAVNO(CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,
     &                  DFRSDEN,DFRSTEM,ISTATE,NGL)

! Unused subroutine

      IMPLICIT NONE

      REAL :: CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL
      REAL :: DFRSDEN,DFRSTEM,DFRSPRE,RO,TEMP,
     &        PSAT,PSATT,PSATP1,PSATM1,DPSAT,RELAX
      INTEGER :: ISTATE,NGL,ICYCLE
      CHARACTER(*),PARAMETER :: TEXT0 = "Error in p_sat = ",
     &    TEXT1 = " Convergence in finding CAVNO is poor in block ",
     &    TEXT2 = " Select state = 8 in finding CAVNO,   CAVNO = 0."

! ... Free-stream temperature as a function of cavitation number
! ... Requires (ISTATE = 8)

      RO      = FRSDEN
      TEMP    = FRSTEM
      RELAX   = 1.
      DFRSPRE = FRSPRE
      IF (ISTATE /= 8) THEN ; PRINT"(A)",TEXT2 ; RETURN ; END IF

      DO ICYCLE = 1,25
       IF (ICYCLE > 20) RELAX = .1
!       PSATT = FRSPRE - CAVNO*.5*RO*FRSVEL**2
       PSATT = -CAVNO*.5*RO*FRSVEL**2
       CALL DPSATF(PSAT,TEMP,DFRSPRE,8)
       CALL DPSATF(PSATP1,TEMP+1.,DFRSPRE,8)
       CALL DPSATF(PSATM1,TEMP-1.,DFRSPRE,8)
       DPSAT = (PSATP1-PSATM1)*.5
       TEMP  = TEMP + (PSATT-PSAT)/DPSAT*RELAX
       RO    = (-.260031e-2*TEMP+.122284e+1)*TEMP+.864232e+3
     &       + ( .761625e-9*TEMP+.196854e-6)*DFRSPRE
c       WRITE(666,*) NGL,PSAT,TEMP,PSATT-PSAT,RO
      END DO

      IF (ABS(PSATT-PSAT) > 1e+3) THEN
       WRITE(*,*) TEXT1,NGL,TEXT0,PSATT-PSAT
       WRITE(4,*) TEXT1,NGL,TEXT0,PSATT-PSAT
      END IF

       FRSDEN = RO ;  FRSTEM = TEMP
      DFRSDEN = RO ; DFRSTEM = TEMP

      END SUBROUTINE FCAVNO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FPRESS(CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,
     &                  DFRSDEN,DFRSTEM,PSAT,ISTATE,NGL,GROUND,CAVLEV)

! ... Free-stream pressure as a function of free-stream temperature
! ... and cavitation number
! ... Requires (ISTATE = 8 / 18)

      IMPLICIT NONE

      INTEGER :: ISTATE,NGL,I
      REAL :: CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,GROUND,CAVLEV,GTOT
      REAL :: DFRSDEN,DFRSTEM,DENS,TEMP,PRES,PSAT,VELS,
     &        FCTR,DELT,PUPP,PMID,PLOW,SLOW,SMID,SUPP,DPH
      CHARACTER(*),PARAMETER :: TEXT0 = "Error in p_sat = ",
     &    TEXT1 = "Convergence in finding FRSPRE is poor in block ",
     &    TEXT2 = "No suitable pressure value is found. Exiting...",
     &    TEXT3 = "Select state = 8 or 18 in finding FRSPRE. Exiting..."

      IF ((ISTATE /= 8) .AND. (ISTATE /= 18)) THEN
       WRITE(6,"(A)") TEXT3 ; STOP
      END IF

      GTOT = 9.80665  ! Antakee mersu
      VELS = 0.5*FRSVEL**2
      PRES = FRSPRE ; TEMP = FRSTEM
      DENS = 1e0/(                       .122993806592151e-02
     &   +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &   +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PRES)
      DPH  = GTOT*(GROUND-CAVLEV)


      CALL DPSATLOG(PSAT,TEMP,ISTATE)

      PLOW = PSAT ; SLOW = -CAVNO
      PUPP = PRES ; SUPP = SLOW + (PRES+DPH*DENS-PSAT)/(DENS*VELS) !DPH

      IF (SLOW*SUPP >= 0.0) THEN
       WRITE(6,"(A)") TEXT2 ; STOP
      END IF

      OPEN(1973,FILE="PresEsts",FORM="FORMATTED")
      WRITE(1973,"(A/)") "Iterating for suitable pressure value"
      WRITE(1973,"(A,T31,F15.3,A)") "Given temperature =",TEMP," K"
      WRITE(1973,"(A,T31,F18.6)") "Given cavitation number =",CAVNO
      WRITE(1973,"(A,T31,F15.3,A)") "Initial pressure =",PRES," Pa"
      WRITE(1973,"(/A/A/2(A,T31,F21.9,A/))")
     &    "Using method of modified regula falsi:","Seed values:",
     &    "Lower limit =",PLOW," Pa","Upper limit =",PUPP," Pa"

      DO I = 1,111
       FCTR = SLOW/(SUPP-SLOW) ; DELT = PUPP-PLOW
       PMID = PLOW - FCTR*DELT
       DENS = 1e0/(                  .122993806592151e-02
     &  +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &  +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PMID)
       SMID = (PMID+DPH*DENS-PSAT)/(DENS*VELS) - CAVNO

      WRITE(1973,"(A,T31,F21.9,A)") "Next estimate =",PMID," Pa"

       IF (SMID < 0.0) THEN
        PLOW = PMID ; SLOW = SMID ; SUPP = -FCTR*SUPP
       ELSE
        PUPP = PMID ; SUPP = SMID ; SLOW = -FCTR*SLOW
       END IF
       IF (ABS(DELT) < 1e-9) EXIT
      END DO

      PRES = PLOW - SLOW*(PUPP+DPH*DENS-PLOW)/(SUPP-SLOW)
      DENS = 1e0/(                   .122993806592151e-02
     &  +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &  +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PRES)

      DENS    = (-.260031e-2*TEMP+.122284e+1)*TEMP+.864232e+3
     &        + ( .761625e-9*TEMP+.196854e-6)*PRES ! Consistent

      WRITE(1973,"(/A,T31,F15.3,A)") "Final pressure value =",PRES," Pa"
      CLOSE(1973)

      FRSDEN = DENS ; DFRSDEN = DENS
      FRSTEM = TEMP ; DFRSTEM = TEMP
      FRSPRE = PRES

      END SUBROUTINE FPRESS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FSIGMA(CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,
     &                  DFRSDEN,DFRSTEM,PSAT,ISTATE,NGL,GROUND,CAVLEV)

! ... Cavitation number as a function of free-stream temperature
! ... Requires (ISTATE = 8 / 18)

      IMPLICIT NONE

      INTEGER :: ISTATE,NGL
      REAL :: CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,GROUND,CAVLEV
      REAL :: DFRSDEN,DFRSTEM,DENS,TEMP,PRES,PSAT,DPH,GTOT
      CHARACTER(*),PARAMETER :: TEXT0 = "Error in p_sat = ",
     &    TEXT1 = "Convergence in finding CAVNO is poor in block ",
     &    TEXT3 = "Select state = 8 or 18 in finding CAVNO. Now CAVNO=0"

      IF ((ISTATE /= 8) .AND. (ISTATE /= 18)) THEN
       WRITE(6,"(A)") TEXT3 ; RETURN
      END IF

      GTOT = 9.80665  ! Antakee mersu
      DPH  = GTOT*(GROUND-CAVLEV)
      PRES = FRSPRE + DPH*FRSDEN ; TEMP = FRSTEM
      DENS = 1e0/(                     .122993806592151e-02
     &    +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &    +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PRES)

!      FRSDEN  = 1e0/(                     .122993806593339e-2
!     &    +( .3463480979514521e-08*FRSTEM-.1792563689396157e-05)*FRSTEM
!     &    +(-.7468375989633927e-14*FRSTEM+.1985224115842062e-11)*FRSPRE)

      CALL DPSATLOG(PSAT,TEMP,ISTATE)

      DENS    = (-.260031e-2*TEMP+.122284e+1)*TEMP+.864232e+3
     &        + ( .761625e-9*TEMP+.196854e-6)*PRES ! Consistent


      FRSDEN = DENS ; DFRSDEN = DENS ; DFRSTEM = TEMP

      CAVNO = 2.0*(FRSPRE-PSAT)/(FRSDEN*FRSVEL**2)

      OPEN(1972,FILE="CavitNum",FORM="FORMATTED")
      WRITE(1972,"(A/)") "Calculating the cavitation number value"
      WRITE(1972,"(A,T31,F15.3,A)") "Given temperature =",TEMP," K"
      WRITE(1972,"(A,T31,F15.3,A)") "Given pressure =",PRES," Pa"
      WRITE(1972,"(A,T31,F18.6)") "Computed cavitation number =",CAVNO
      CLOSE(1972)

      END SUBROUTINE FSIGMA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FTEMPS(CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,
     &                  DFRSDEN,DFRSTEM,ISTATE,NGL,GROUND,CAVLEV)

! ... Free-stream temperature as a function of cavitation number
! ... Requires (ISTATE = 8 / 18)

      IMPLICIT NONE

      INTEGER :: ISTATE,NGL,I
      REAL :: CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,GROUND,CAVLEV,PIRO
      REAL :: DFRSDEN,DFRSTEM,DENS,TEMP,PRES,PSAT,VELS,
     &        FCTR,DELT,TUPP,TMID,TLOW,SLOW,SMID,SUPP,DPH,GTOT
      CHARACTER(*),PARAMETER :: TEXT0 = "Error in p_sat = ",
     &    TEXT1 = "Convergence in finding CAVNO is poor in block ",
     &    TEXT2 = "No suitable temperature value is found. Exiting...",
     &    TEXT3 = "Select state = 8 or 18 in finding CAVNO. Exiting..."

      IF ((ISTATE /= 8) .AND. (ISTATE /= 18)) THEN
       WRITE(6,"(A)") TEXT3 ; STOP
      END IF

      GTOT = 9.80665  ! Antakee mersu
      DFRSDEN = FRSDEN
      DPH  = GTOT*(CAVLEV-GROUND)
      VELS = 0.5*FRSVEL**2
      PRES = FRSPRE + DPH*DFRSDEN
      TEMP = FRSTEM
      DENS = 1e0/(                        .122993806592151e-02
     &    +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &    +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PRES)

      CALL DPSATLOG(PSAT,TEMP,ISTATE)
      CALL DTSATLOG(TUPP,PRES,ISTATE) ; SUPP = -CAVNO
      TLOW = TEMP ; SLOW = SUPP + (PRES-PSAT)/(DENS*VELS)

      IF (SLOW*SUPP >= 0.0) THEN
       WRITE(6,"(A)") TEXT2 ; STOP
      END IF
! tps 19.5.2009
! Added IF (NGL == 1) THEN do not reopen the opened file
      IF (NGL == 1) THEN
      OPEN(1971,FILE="TempEsts",FORM="FORMATTED")
      WRITE(1971,"(A/)") "Iterating for suitable temperature value"
      WRITE(1971,"(A,T31,F15.3,A)") "Given pressure =",PRES," Pa"
      WRITE(1971,"(A,T31,F15.3,A)") "Given refvel =",FRSVEL," m/s"
      WRITE(1971,"(A,T31,F18.6)") "Given cavitation number =",CAVNO
      WRITE(1971,"(A,T31,F15.3,A)") "Initial temperature =",TEMP," K"
      WRITE(1971,"(/A/A/2(A,T31,F21.9,A/))")
     &    "Using method of modified regula falsi:","Seed values:",
     &    "Lower limit =",TLOW," K","Upper limit =",TUPP," K"
      ENDIF

      DO I = 1,111
       FCTR = SLOW/(SUPP-SLOW) ; DELT = TUPP-TLOW
       TMID = TLOW - FCTR*DELT ; CALL DPSATLOG(PSAT,TMID,ISTATE)
       DENS = 1e0/(                  .122993806592151e-02
     &  +( .346348097967690e-08*TMID-.179256368943405e-05)*TMID
     &  +(-.746837644005658e-14*TMID+.198522424218165e-11)*PRES)
       SMID = (PRES-PSAT)/(DENS*VELS) - CAVNO

      IF (NGL == 1) THEN
         WRITE(1971,"(A,T31,2(F21.9,A))") "Next estimate =",TMID," K" ,
     &   REAL(DENS,4)," kg/m^3"
      ENDIF
      IF (SMID > 0.0) THEN
        TLOW = TMID ; SLOW = SMID ; SUPP = -FCTR*SUPP
       ELSE
        TUPP = TMID ; SUPP = SMID ; SLOW = -FCTR*SLOW
       END IF
       IF (ABS(DELT) < 1e-9) EXIT
      END DO

      TEMP = TLOW - SLOW*(TUPP-TLOW)/(SUPP-SLOW)
      DENS = 1e0/(                   .122993806592151e-02
     &  +( .346348097967690e-08*TEMP-.179256368943405e-05)*TEMP
     &  +(-.746837644005658e-14*TEMP+.198522424218165e-11)*PRES)

      IF (NGL == 1) THEN
           WRITE(1971,"(/A,T31,F15.3,A)")
     &    "Final temperature value =",TEMP," K"
      CLOSE(1971)
      ENDIF
      
      CALL DPSATLOG(PSAT,TEMP,ISTATE)

      PIRO    = (-.260031e-2*TEMP+.122284e+1)*TEMP+.864232e+3
     &        + ( .761625e-9*TEMP+.196854e-6)*PRES
 
C      FRSDEN = REAL(DENS) ; DFRSDEN = DENS ! This is not consistent

      FRSDEN = PIRO ; DFRSDEN = PIRO
      FRSTEM = TEMP ; DFRSTEM = TEMP

      END SUBROUTINE FTEMPS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FTEMPT(CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,
     &                  DFRSDEN,DFRSTEM,DPSAT,ISTATE,NGL,GROUND,CAVLEV)

! ... Free-stream temperature as a function of cavitation number
! ... Requires (ISTATE = 8 / 18). Uses only temperatures

      IMPLICIT NONE

      INTEGER :: ISTATE,NGL,I
      REAL :: CAVNO,FRSTEM,FRSPRE,FRSDEN,FRSVEL,DPSAT,GROUND,CAVLEV
      REAL :: DFRSDEN,DFRSTEM,DENS,TEMP,PRES,PSAT,VELS,
     &        FCTR,DELT,TUPP,TMID,TLOW,SLOW,SMID,SUPP,DPH,GTOT
      REAL :: TS1,TS2,DPSDT,DTSAT,DPCAV,TEMD
      CHARACTER(*),PARAMETER :: TEXT0 = "Error in p_sat = ",
     &    TEXT1 = "Convergence in finding CAVNO is poor in block ",
     &    TEXT2 = "No suitable temperature value is found. Exiting...",
     &    TEXT3 = "Select state = 8 or 18 in finding CAVNO. Exiting..."
      REAL, PARAMETER ::
     &    P1 =  0.154906419868284e+6,
     &    Q1 = -0.101158884083608e+5,
     &    R1 =  0.133000059771258e+3,
     &    S1 =  0.146861326101013e-1,
     &    T1 = -0.170477189812508e+2,
     &    V1 =  0.249899160714079e+2,
     &    W1 = -0.504070497017732e+4

      IF ((ISTATE /= 8) .AND. (ISTATE /= 18)) THEN
       WRITE(6,"(A,I2)") TEXT3,ISTATE ; STOP
      END IF
      
      DFRSTEM = FRSTEM ; DPCAV = FRSPRE
      GTOT = 9.80665  ! Antakee mersu
      DPH  = GTOT*(CAVLEV-GROUND)

      DO I = 1,11
      DFRSDEN = (-.260031e-2*DFRSTEM+.122284e+1)*DFRSTEM+.864232e+3
     &        + ( .761625e-9*DFRSTEM+.196854e-6)*(FRSPRE+DPH*DFRSDEN)
      DTSAT       = W1/(LOG(DPCAV) - V1)
      TS1         = W1/(LOG(DPCAV+1.) - V1)
      TS2         = W1/(LOG(DPCAV-1.) - V1)

      DPSDT= 1./((TS1-TS2)/2.)
      DFRSTEM = -CAVNO*.5*DFRSDEN*FRSVEL**2/DPSDT + DTSAT
      ENDDO

      TEMD   = DFRSTEM
      DPSAT  = EXP((P1/TEMD+Q1)/TEMD + R1 + S1*TEMD + T1*LOG(TEMD))
      DPSAT  = MAX(.1e-4,DPSAT)

      FRSTEM = DFRSTEM ; ! FRSDEN = DFRSDEN

      FRSDEN  = (-.260031e-2*TEMD+.122284e+1)*TEMD+.864232e+3
     &        + ( .761625e-9*TEMD+.196854e-6)*(FRSPRE+DPH*DFRSDEN)
      DFRSDEN = FRSDEN

      END SUBROUTINE FTEMPT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

