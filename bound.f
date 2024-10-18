C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CALCIOP(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     &   BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,BOUNA2,BOUNG,BOUNRET,
     &   MAXEB,MAXSB,IF1,NSCAL,ITURB,TRANSL,MULPHL,
     &   X,Y,Z,UROT,KX1,KX2,KY1,KY2,
     &   A1,A1X,A1Y,A1Z,IX,IY,IZ,NGL,NWALL,
     &   XMASS,HTOT,PTOT,P,TEMP,TURBLE,RMUINI,RMUL1,RMUL2,TG,TRET,
     &   LRATUP,IFILE,INPTYP,APU,FRSPRE,FRSTEM,FRSVEL,
     &   FRSDEN,FRSSIE,FRSVIS,IUPPTL)

C ... Calculate INL/OUT boundary value distribution starting from user  
C ... given flow condition parameters (see subroutine READIOP).

      USE NS3CO, ONLY : ISTATE, RMACH, RMULTV, ! Some variables moved
     &                  RGAS, GAMMA, VISU0,    ! to the parameter list.
     &                  EXPSU, TSU0, E0REF, T0REF, ISTRES, IN, JN, KN

      IMPLICIT NONE   

      REAL, DIMENSION(*) :: X, Y, Z

C ... Arrays for grid cell geometries, additional scalars and
C ... Reynolds stress components:

      REAL, DIMENSION(*)  :: A1, A1X, A1Y, A1Z, UROT
      REAL, DIMENSION(NSCAL) :: FISET
      REAL, DIMENSION(5)  :: BIJ

      REAL :: RHO, RHOU, E, RK, REPS, P, TEMP, ARE, LRATUP
      REAL :: XMASS, HTOT, PTOT, TURBLE, RMUINI, RMUL1, RMUL2
      REAL :: XROT, YROT, ZROT, OMX, OMY, OMZ, TG, TRET
      REAL :: FRSPRE, FRSTEM, FRSVEL, FRSDEN, FRSSIE, FRSVIS

      INTEGER :: IX, IY, IZ, KX1, KX2, KY1, KY2
      INTEGER :: NB, IL, NGL, NWALL, INPTYP, IFILE, I 
      INTEGER :: MAXEB, MAXSB, IF1, NSCAL, ITURB, IUPPTL
      
      REAL :: BOUNMF(*),BOUNR(*),BOUNU(*),BOUNV(*),BOUNW(*),BOUNE(*),
     +        BOUNT(*),BOUNP(*),BOUNRK(*),BOUNEP(*),BOUNA1(*),
     +        BOUNA2(*),BOUNFI(MAXSB,*),BOUNBI(MAXEB,*),
     +        BOUNG(*),BOUNRET(*),APU(*)

      LOGICAL :: TRANSL, MULPHL

      CALL AREAS(A1X,A1Y,A1Z,A1,IX,IY,IZ,IN,JN,KN,NWALL,ARE)

C ... Compute the relevant 1-D physical quantities at the boundary

      CALL PHYCON(RHO,RHOU,E,RK,REPS,P,TEMP,ARE,XMASS,HTOT,PTOT,
     +     INPTYP,ISTATE,ITURB,TURBLE,RMUINI,FRSVIS,FRSVEL,
     +     RGAS,GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,APU,FRSPRE,
     +     FRSTEM,IUPPTL,RMULTV)

      IF(LRATUP < 0.1) THEN

C ... Generate a simple constant plug-type flow for the whole area

      CALL OUTPLG(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,
     +     BOUNA2,BOUNG,BOUNRET,MAXEB,MAXSB,
     +     RHO,RHOU,E,RK,REPS,FISET,A1,A1X,A1Y,A1Z,
     +     X,Y,Z,UROT,KX1,KX2,KY1,KY2,IX,IY,IZ,IN,JN,KN,
     +     NSCAL,ITURB,TRANSL,MULPHL,NWALL,NGL,RGAS,GAMMA,
     +     P,TEMP,BIJ,ISTRES,IFILE,ISTATE,FRSDEN,FRSPRE,E0REF,T0REF,
     +     FRSVEL,RMUL1,RMUL2,TG,TRET,FRSTEM,IUPPTL)

      ELSE

C ... Create modified distributions based on developing tube flow

      CALL OUTMOD(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,
     +     BOUNA2,BOUNG,BOUNRET,MAXEB,MAXSB,
     +     RHO,RHOU,E,RK,REPS,FISET,A1,A1X,A1Y,A1Z,
     +     X,Y,Z,UROT,KX1,KX2,KY1,KY2,IX,IY,IZ,IN,JN,KN,
     +     NSCAL,ITURB,TRANSL,MULPHL,NWALL,NGL,RGAS,GAMMA,
     +     P,TEMP,BIJ,ISTRES,IFILE,ISTATE,FRSDEN,FRSPRE,E0REF,T0REF,
     +     FRSVEL,RMUL1,RMUL2,TG,TRET,LRATUP,ARE,FRSVIS,FRSTEM,IUPPTL)

      ENDIF

      RETURN 
      END SUBROUTINE CALCIOP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AREAS(A1X,A1Y,A1Z,A1,IX,IY,IZ,
     &                 INRE,JNRE,KNRE,NWALL,ARECRS)
      
C ... This subroutine determines the block face area and its direction.

      IMPLICIT NONE

      INTEGER :: IX, IY, IZ, NWALL, I, J, K, ISTR, JSTR, KSTR, IA, IG
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, INRE, JNRE, KNRE

      REAL, DIMENSION(*) :: A1X, A1Y, A1Z, A1
      REAL :: AREMEM, ARECRS, SAX, SAY, SAZ

C     Determine the indexing for the block face in question
C     recognizing the ghost cells like in AREVC1

      IF     (NWALL == 1 .OR. NWALL == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         IMAX = IY
         JMAX = IZ
         KMAX = IX
         ISTR = KMAX+2*KN
         JSTR = (KMAX+2*KN)*(IMAX+2*IN)
         KSTR = 1
      ELSEIF (NWALL == 2 .OR. NWALL == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         IMAX = IZ
         JMAX = IX
         KMAX = IY
         ISTR = (JMAX+2*JN)*(KMAX+2*KN)
         JSTR = 1
         KSTR = JMAX+2*JN
      ELSEIF (NWALL == 3 .OR. NWALL == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         IMAX = IX
         JMAX = IY
         KMAX = IZ
         ISTR = 1
         JSTR = IMAX+2*IN
         KSTR = (IMAX+2*IN)*(JMAX+2*JN)
      ENDIF

C     Select between index minima or maxima at the wall to define
C     a global starting address for the cell face areas 

      K = 1
      IF (NWALL > 3) K = KMAX
      IA = 1 + (K-1+KN)*KSTR

C     Accumulate the cell area components for averaging

      AREMEM = 0.
      SAX    = 0.
      SAY    = 0.
      SAZ    = 0.
      DO 20 J = 1,JMAX
         DO 10 I = 1,IMAX
         IG = IA + (I-1+IN)*ISTR + (J-1+JN)*JSTR
         AREMEM = AREMEM + A1(IG)
         SAX    = SAX + (A1X(IG)*A1(IG))
         SAY    = SAY + (A1Y(IG)*A1(IG))
         SAZ    = SAZ + (A1Z(IG)*A1(IG))
 10      CONTINUE
 20   CONTINUE

C ... Calculate the averaged face unit normal by dividing the area
C ... component sums with the accumulated membrane area

      SAX = SAX/AREMEM
      SAY = SAY/AREMEM
      SAZ = SAZ/AREMEM

C ... Integrate the effective cross-section area of the potentially
C ... non-flat boundary using the averaged wall normal

      ARECRS = 0.

      DO 40 J = 1,JMAX
         DO 30 I = 1,IMAX
         IG = IA + (I-1+IN)*ISTR + (J-1+JN)*JSTR
         ARECRS = ARECRS + SAX*(A1X(IG)*A1(IG)) + SAY*(A1Y(IG)*A1(IG))
     &                   + SAZ*(A1Z(IG)*A1(IG))
 30      CONTINUE
 40   CONTINUE

C ... If AREMEM is appreciably larger than ARECRS, the grid on the 
C ... boundary is not flat and the resulting final distributions do
C ... not represent constant flow across the boundary. 

 1000 FORMAT(4X,'Membrane area of the boundary',6X,E10.3/
     +       4X,'Cross-section area of the boundary',1X,E10.3/
     +       4X,'Average direction of normal vector is:'/
     +       10X,'Nx',8X,'Ny',8X,'Nz'/
     +       4X,3F10.3/)

      RETURN
      END SUBROUTINE AREAS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PHYCON(RHO,RHOU,E,RK,REPS,P,TEMP,ARE,XMASS,HTOT,PTOT,
     +           INPTYP,ISTATE,ITURB,TURBLE,RMUINI,FRSVIS,FRSVEL,
     +           RGAS,GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,APP,FRSPRE,
     +           FRSTEM,IUPPTL,RMULT)

C ... Compute the relevant 1-D physical quantities at the boundary.
C ... The flow is perpendicular to the cross-sectional area ARE and 
C ... defined as positive into the grid.

      IMPLICIT NONE

      INTEGER :: INPTYP, ISTATE, ITURB, IUPPTL

      REAL :: RHO, RHOU, E, RK, REPS, P, TEMP
      REAL :: ARE, XMASS, UVEL, HTOT, PTOT
      REAL :: TURBLE, RMUINI, FRSVIS, FRSVEL, FRSPRE, FRSTEM
      REAL :: ABSXM, APU
      REAL :: RMACH, TTOT, RHOTOT,RGAS, GAMMA, GM1
      REAL :: VISU0, EXPSU, TSU0, E0REF, T0REF, RMULT

      REAL :: APP(*)


C ... Convert the directional mass flow into rho*unormal

      XMASS = XMASS/ARE
      ABSXM = ABS(XMASS)

C ... Determine the flow state according to the input parameters 
C ... applying the specified thermodynamic fluid model. The general 
C ... idea is to employ suitable equations of state rho = f(p,T)
C ... and to enable sets of primitive variables (p,u,T +turb) or 
C ... conservative variables (rho,rhou,rhoetot +turb) to be passed 
C ... forward. The mass flow here means always its absolute value. 

      IF(INPTYP == 1) THEN
C     Input as rho*u, htot and p (Ideal and real gas possible):
      CALL BTYPE1(ABSXM,HTOT,P,RHO,RHOU,E,RK,
     +            TURBLE,FRSVEL,GAMMA,ISTATE)
      TEMP  = P/(RGAS*RHO)
      ELSEIF(INPTYP == 2) THEN
C     Input as rho*u, htot and ptot (Ideal gas only):
      CALL BTYPE2(ABSXM,HTOT,PTOT,P,TEMP,RHO,RHOU,E,RK,
     +            TURBLE,RGAS,GAMMA,ISTATE)
      ELSEIF(INPTYP == 3) THEN
C     Input as rho*u, htot and T (Ideal gas only):
      CALL BTYPE3(ABSXM,HTOT,TEMP,P,RHO,RHOU,E,RK,
     +            TURBLE,RGAS,GAMMA,ISTATE)
      ELSE
C     Input as rho*u, p and T (5 ISTATE options possible):
      CALL BTYPE4(ABSXM,P,TEMP,RHO,RHOU,E,RK,TURBLE,FRSVIS,
     +     VISU0,EXPSU,TSU0,E0REF,T0REF,RGAS,GAMMA,ISTATE,
     +     APP,FRSPRE,FRSTEM,IUPPTL,RMULT)
      ENDIF

C ... Take the flow direction into account and determine eps or omega
      RHOU  = ABSXM/XMASS*RHOU
      UVEL  = RHOU/RHO
      REPS  = 0.09*RK**2/(RMUINI*FRSVIS)
      IF (ITURB == 6) THEN
C ... Instead of rho*epsilon, rho*omega is stored in REPS
         REPS = RHO*RK/(RMUINI*FRSVIS)
      END IF

      IF(ISTATE <= 5) THEN
C ... Calculate and report additional flow variables for 
C ... compressible flow to facilitate checking
      RMACH = ABS(UVEL)/SQRT(GAMMA*P/RHO)
      GM1   = GAMMA - 1.
      PTOT  = P*(1. + GM1/2.*RMACH**2)**(GAMMA/GM1)
      RHOTOT= RHO*(1. + GM1/2.*RMACH**2)**(1./GM1)
      TTOT  = TEMP*(1+GM1/2.*RMACH**2)
      ENDIF


      RETURN
      END SUBROUTINE PHYCON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BTYPE1(XMASS,HTOT,P,RHO,RHOU,E,RK,
     +           TURBLE,FRSVEL,GAMMA,ISTATE)

C ... Determine the flow state based on rho*u, htot and p using the
C ... thermodynamic model specified by ISTATE

      IMPLICIT NONE

      INTEGER :: ISTATE,I 

      REAL :: P,RHO,RHOU,E,RK, HTOT, XMASS
      REAL :: UVEL,RHOOLD,ESMA,ESMAOLD
      REAL :: TURBLE,FRSVEL,GAMMA,GM1
      REAL :: R23,R53,APU1,APU2,APU9

      GM1  = GAMMA - 1.
      R53  = 5./3.
      R23  = 2./3.
      APU9 = 0.

C ... Initial guesses for density and specific internal energy
C ... ESMA are based on ideal gas relations. The resulting 
C ... algebraic second-order equation is solved as follows:

      RK   = 1.5*(TURBLE*FRSVEL)**2
      APU1 = 0.5*(XMASS*GM1/P)**2
      APU2 = SQRT(GAMMA**2 + 4.*APU1*(HTOT-R53*RK))
      ESMA = 0.5*(-GAMMA + APU2)/APU1
      RHO  = P/(ESMA*GM1)

      IF(ISTATE == 1) WRITE(*,*) 'Perfect gas relation is used'
      IF(ISTATE == 2) WRITE(*,*) 'Real gas relation is used'

      I = 0

C ... BEGIN ITERATION

 10   RHOOLD  = RHO
      ESMAOLD = ESMA
      I = I + 1
      UVEL = XMASS/RHO
C ... EQUATION OF STATE (Requires bits from state.f and tgas.f to
C ... provide e = f(p,rho)
      CALL EFPRO(ESMA,P,RHO,1,ISTATE,GAMMA,0.,0.,0.)
C ... RELAXATION TO REDUCE OVERSHOOTS:
      ESMA = (ESMA*.1 + ESMAOLD*.9)
      APU1 = (HTOT - ESMA - R53*RK)
      APU2 = SQRT(P**2 + 2.*APU1*XMASS**2)
      RHO  = 0.5*(P + APU2)/APU1
      RK   = 1.5*(TURBLE*UVEL)**2
C     Write and follow the status of density iteration
*      WRITE(8,1313) I,RHO,ESMA,RHO-RHOOLD,RK,UVEL
      IF (I > 100) APU9 = APU9 + ABS(RHO-RHOOLD)
      IF (ABS(RHO-RHOOLD) >= (5.E-6*RHO) .AND. I <= 1000) GOTO 10
C ... ITERATION IS FINISHED
CCC   The iteration does not work convincingly. In a test it got stuck
CCC   in an oscillation after about 20 iterations and then ran 980 
CCC   rounds in vain. Because the original convergence criterion of 
CCC   5.E-7 appears too tight, it was relaxed to 5.E-6. JH 10/2013

      IF (I >= 1000 .AND. APU9/900. >= 5.E-6*RHO) THEN
          WRITE(*,*)'Solution did not converge. Try a re-run'//
     +    'with different input values or re-program the code.'
          WRITE(*,*) 'Exiting...'
          CALL EXIT
      ENDIF

C ... SET THE REMAINING VARIABLES

      RHOU = XMASS
      RK   = RHO*RK
      E    = 0.5*RHO*ESMA + RHOU**2/RHO + RK

1313  FORMAT(I4,1X,E10.5,1X,F14.3,1X,E11.5,1X,2F11.5,1X,E11.5)

      RETURN
      END SUBROUTINE BTYPE1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BTYPE2(XMASS,HTOT,PTOT,P,TEMP,RHO,RHOU,E,RK,
     +           TURBLE,RGAS,GAMMA,ISTATE)

C ... Determine the flow state based on rho*u, htot and ptot. In this
C ... version corresponding to newbound81, only the ideal gas model is 
C ... applied and ISTATE is a dummy. 

      IMPLICIT NONE

      INTEGER :: ISTATE,ISUM 

      REAL :: XMASS,HTOT,PTOT,P,TEMP,RHO,RHOU,E,RK
      REAL :: TURBLE,RGAS,GAMMA,GP1,GM1,GP1GM1,CP
      REAL :: TTOT,HSTAT
      REAL :: FH,FHP,C0,ALPHA

      GP1    = (GAMMA + 1.)
      GM1    = (GAMMA - 1.)
      GP1GM1 = GP1/GM1
      CP     = GAMMA*RGAS/GM1

      TTOT   = HTOT/CP
      C0     = PTOT*(RGAS*TTOT/PTOT)**GAMMA
      ALPHA  = 2.*(GM1/GAMMA/C0)**(2./GM1)

C ... Iterate the static enthalpy

      HSTAT  = HTOT
      ISUM   = 0

 10   FH  = ALPHA*HSTAT**GP1GM1 - ALPHA*HTOT*HSTAT**(2./GM1) + XMASS**2
      FHP   = GP1GM1*ALPHA*HSTAT**(2./GM1) 
     +      - (2./GM1)*ALPHA*HTOT*HSTAT**(-(GAMMA-3.)/GM1)
      HSTAT = HSTAT - FH/FHP
      ISUM  = ISUM + 1

      IF(ABS(FH/HTOT) >=  1.E-5 .AND. ISUM <= 1000) GOTO 10
CC      IF(ABS(FH/HTOT) >=  1.E-5 .AND. ISUM <= 10000) GOTO 10

C ... Primitive variables

      RHO   = (GM1/GAMMA/C0 * HSTAT)**(1./GM1)
      P     = C0*RHO**GAMMA
      TEMP  = P/(RHO*RGAS)

C ... Remaining conservative variables

      RHOU  = XMASS
      RK    = 1.5*(TURBLE*RHOU)**2/RHO
      E     = P/GM1 + 0.5*RHOU**2/RHO + RK   

      RETURN
      END SUBROUTINE BTYPE2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BTYPE3(XMASS,HTOT,TEMP,P,RHO,RHOU,E,RK,
     +           TURBLE,RGAS,GAMMA,ISTATE)

C ... Determine the flow state based on rho*u, htot and T. In this
C ... version corresponding to newbound81, only the ideal gas model is 
C ... applied and ISTATE is a dummy. 

      IMPLICIT NONE

      INTEGER :: ISTATE 

      REAL :: XMASS,HTOT,TEMP,P,RHO,RHOU,E,RK
      REAL :: TURBLE,RGAS,GAMMA,GM1,CP
      REAL :: TTOT,RMACH,C

      GM1    = (GAMMA - 1.)
      CP     = GAMMA*RGAS/GM1

      TTOT   = HTOT/CP
      RMACH  = SQRT(2./GM1*(TTOT/TEMP-1.))
      C      = SQRT(GAMMA*RGAS*TEMP)

C ... Primitive variables
      RHO    = XMASS/(RMACH*C)
      P      = C**2*RHO/GAMMA

C ... Remaining conservative variables
      RHOU  = XMASS
      RK    = 1.5*(TURBLE*RHOU)**2/RHO
      E     = P/GM1 + 0.5*RHOU**2/RHO + RK  

      RETURN
      END SUBROUTINE BTYPE3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BTYPE4(XMASS,P,TEMP,RHO,RHOU,E,RK,TURBLE,FRSVIS,
     +           VISU0,EXPSU,TSU0,E0REF,T0REF,RGAS,GAMMA,ISTATE,
     +           APP,FRSPRE,FRSTEM,IUPPTL,RMULT)

C ... Determine the flow state based on rho*u, p and T. In this
C ... version corresponding to newbound81, ISTATE options 1 and 
C ... 5 to 9 are active.

      IMPLICIT NONE

      INTEGER :: ISTATE,ERCODE, IUPPTL

      REAL :: XMASS,P,TEMP,RHO,RHOU,E,RK
      REAL :: TURBLE,RGAS,GAMMA
      REAL :: FRSVIS,ESMA,CP,CH,DRDP,DRDH,PR
      REAL :: VISU0,EXPSU,TSU0,E0REF,T0REF,FRSPRE,FRSTEM,RMULT

      REAL :: APP(*)

C ... Equation of state to obtain mainly rho and e as f(p,T). 
C ... Requires routines from state.f.
C ... A default Prandtl number is set here

      PR = 0.7

      IF(ISTATE == 1) THEN
         CALL PERAIR(FRSVIS,P,RHO,ESMA,TEMP,CP,CH,DRDP,DRDH,
     &   RGAS,GAMMA,PR,VISU0,EXPSU,TSU0,1,ERCODE,E0REF,T0REF)
      ELSE IF(ISTATE == 5) THEN
         CALL PERFGA(FRSVIS,P,RHO,ESMA,TEMP,CP,CH,DRDP,DRDH,
     &   RGAS,GAMMA,PR,VISU0,EXPSU,TSU0,1,ERCODE,E0REF,T0REF,
     &   RMULT,APP,FRSPRE,FRSTEM,IUPPTL)
      ELSE IF(ISTATE == 6) THEN
         CALL WATER (FRSVIS,P,RHO,ESMA,TEMP,CP,CH,DRDP,DRDH,PR,1)
      ELSE IF(ISTATE == 7) THEN
         CALL WATER (FRSVIS,P,RHO,ESMA,TEMP,CP,CH,DRDP,DRDH,PR,1)
      ELSE IF(ISTATE == 8) THEN 
         CALL WATER2(FRSVIS,P,RHO,ESMA,TEMP,CP,CH,DRDP,DRDH,PR,1)
      ELSE
         STOP
      ENDIF

C ... Remaining conservative variables

      RHOU = XMASS
      RK   = 1.5*(TURBLE*RHOU)**2/RHO
      E    = RHO*ESMA + 0.5*RHOU**2/RHO + RK   

      RETURN
      END SUBROUTINE BTYPE4
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTPLG(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,
     +     BOUNA2,BOUNG,BOUNRET,MAXEB,MAXSB,
     +     RHO,RHOU,E,RK,REPS,FISET,A1,A1X,A1Y,A1Z,
     +     X,Y,Z,UROT,KX1,KX2,KY1,KY2,IX,IY,IZ,INRE,JNRE,KNRE,
     +     NSCAL,ITURB,TRANSL,MULPHL,NWALL,NGL,RGAS,GAMMA,
     +     P,TEMP,BIJ,ISTRES,IFILE,ISTATE,FRSDEN,FRSPRE,E0REF,T0REF,
     +     FRSVEL,RMUL1,RMUL2,TG,TRET,FRSTEM,IUPPTL)

C ... The cell-wise boundary condition result file for FINFLO is 
C ... created assuming a simple constant plug-type flow over the whole
C ... boundary area. Either primitive or conservative primary variables 
C ... are used as specified by IFILE. The final velocity or momentum 
C ... components derived from the relative contravariant momentum
C ... RHOU into the grid are defined in the fixed inertial frame and  
C ... thus contain the effects of grid rotation.

      USE CONSTANTS, ONLY : EPS
      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: IX, IY, IZ, IMAX, JMAX, KMAX, KX1, KX2, KY1, KY2
      INTEGER :: ISTR, JSTR, KSTR
      INTEGER :: NSCAL, ITURB, NWALL, NGL, NS, ISTRES, ISTATE 
      INTEGER :: IG0, NR, IA0, NP, I, J, K, IT1, IT2, IT3, IT4, IFILE
      INTEGER :: IN, JN, KN, INRE, JNRE, KNRE
      INTEGER :: MAXEB, MAXSB, N, IUPPTL

      REAL, DIMENSION(*) :: X, Y, Z
      REAL, DIMENSION(*) :: A1, A1X, A1Y, A1Z, FISET, BIJ
      REAL, DIMENSION(*) :: UROT

      REAL :: RHO, RHOU, E, EL, RK, REPS, RGAS, GAMMA, P, TEMP

      REAL :: XC, YC, ZC, ROIN, UVEL, VVEL, WVEL, FACEVEL, DIRIN,
     &        FRSDEN, FRSPRE, E0REF, T0REF, TU, FRSVEL, RMUL1, RMUL2,
     &        TG, TRET, FRSTEM
      
      REAL :: BOUNMF(*),BOUNR(*),BOUNU(*),BOUNV(*),BOUNW(*),BOUNE(*),
     +     BOUNT(*),BOUNP(*),BOUNRK(*),BOUNEP(*),BOUNA1(*),
     +     BOUNA2(*),BOUNFI(MAXSB,*),BOUNBI(MAXEB,*),BOUNG(*),BOUNRET(*)

      LOGICAL :: MULPHL, TRANSL

      CHARACTER(LEN=12) :: NAME

C ... Determine suitable indexing systems for the three wall types.
C ... STEPs are for the grid nodes that do not recognize ghost cells.
C ... STRs are for the cell face areas that are arranged to take the
C ... ghost cells into account similarly as in AREAS and AREVC1.     

      IF     (NWALL == 1 .OR. NWALL == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         IMAX = IY
         JMAX = IZ
         KMAX = IX
         ISTR = KMAX+2*KN
         JSTR = (KMAX+2*KN)*(IMAX+2*IN)
         KSTR = 1
      ELSEIF (NWALL == 2 .OR. NWALL == 5)THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         IMAX = IZ
         JMAX = IX
         KMAX = IY
         ISTR = (JMAX+2*JN)*(KMAX+2*KN)
         JSTR = 1
         KSTR = JMAX+2*JN
      ELSEIF (NWALL == 3 .OR. NWALL == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         IMAX = IX
         JMAX = IY
         KMAX = IZ
         ISTR = 1
         JSTR = IMAX+2*IN
         KSTR = (IMAX+2*IN)*(JMAX+2*JN)

      ELSE
         CALL EXIT
      ENDIF

C ... Initialization of the main variable

      DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
         BOUNR(N)  = 0.
         BOUNU(N)  = 0.
         BOUNV(N)  = 0.
         BOUNW(N)  = 0.
         BOUNE(N)  = 0.
         BOUNT(N)  = 0.
         BOUNP(N)  = 0.         
         BOUNMF(N) = 0.
         BOUNA1(N) = 0.
         BOUNA2(N) = 0.
      ENDDO

C ... Initialization of the turbulence variables

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNRK(N) = 0.
            BOUNEP(N) = 0.
         ENDDO
      ENDIF

C ... Initialization of the Intermittency variables

      IF (TRANSL) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNG(N)   = 0.
            BOUNRET(N) = 0.
         ENDDO
      ENDIF

C ... Initialization of EARSM anisotropic components

      IF(ISTRES > 0) THEN
         DO NS = 1,6
            DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
               BOUNBI(N,NS) = 0.
            ENDDO
         ENDDO
      ENDIF

C ... Initialization of the multiphase variables

      IF (MULPHL) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNA1(N) = 0.
            BOUNA2(N) = 0.
         ENDDO
      ENDIF

C ... Initialization of scalars

      DO NS = 1,NSCAL
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNFI(N,NS) = 0.
         ENDDO
      ENDDO

C ... Select between index minima or maxima at the wall to define
C ... a global starting address for the cell face areas and
C ... define the direction into the grid
      K     = 1
      DIRIN = 1.
      IF (NWALL > 3) THEN
c         K     = KMAX      ! JH 2013
         K     = KMAX + 1  ! ESa 29.1.2016
         DIRIN = -1.
      ENDIF

      IG0   = 1 + (K-1+KN)*KSTR

CCC      CALL INPFIL(XAU1,XAU2,XAU3,XAU4,XAU5,XAU6,XAU7,X,Y,Z,
CCC     +     IMAX,JMAX,KMAX,IX,IY,IZ)

C ... Automatic naming of the output file to contain 
C ... block and wall numbers

      NAME( 1:12) = 'DIST_B000_F0'

      IF(NGL < 10) THEN
         NAME( 9:9 ) = CHAR(NGL+48)
      ELSEIF(NGL >= 10 .AND. NGL < 100) THEN
         NAME( 8:8 ) = CHAR(NGL/10+48)
         NAME( 9:9 ) = CHAR(MOD(NGL,10)+48)
      ELSE
         NAME( 7:7 ) = CHAR(NGL/100+48)
         NAME( 8:8 ) = CHAR(MOD(NGL,100)/10+48)
         NAME( 9:9 ) = CHAR(MOD(NGL,10)+48)
      ENDIF

      NAME(12:12) = CHAR(48+NWALL)

C ... The loop defining the main variables for each boundary cell:

      DO 400 J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

         DO I = KX1,KX2

         IA0 = IG0 + (I-1+IN)*ISTR + (J-1+JN)*JSTR
         NP  = NR + I                     ! Patch index with LN ghost cells


C ... Compute the velocity components resulting from the grid 
C ... rotation 

         IT1 = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT2 = 1+(I  +IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT3 = 1+(I  +IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR
         IT4 = 1+(I-1+IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR

         FACEVEL = UROT(IA0)

         IF(IFILE == 0) THEN

C ... Add the grid velocities to the relative velocities at the 
C ... boundary to obtain absolute inertial velocities before writing
C ... the primitive variables  

         UVEL = RHOU*A1X(IA0)/(RHO+1.E-20)*DIRIN + A1X(IA0)*FACEVEL
         VVEL = RHOU*A1Y(IA0)/(RHO+1.E-20)*DIRIN + A1Y(IA0)*FACEVEL
         WVEL = RHOU*A1Z(IA0)/(RHO+1.E-20)*DIRIN + A1Z(IA0)*FACEVEL

         BOUNP(NP) = P
         BOUNU(NP) = UVEL
         BOUNV(NP) = VVEL
         BOUNW(NP) = WVEL
         BOUNT(NP) = TEMP

         CALL ROFPT(TEMP,P,ROIN,1,ISTATE,RGAS,
     +        GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
         BOUNR(NP) = ROIN

         BOUNMF(NP) = (A1X(IA0)*BOUNU(NP) + A1Y(IA0)*BOUNV(NP) +
     +                 A1Z(IA0)*BOUNW(NP) - FACEVEL)*BOUNR(NP)

         ELSE IF(IFILE == 1) THEN

C ... Modify the total energy to contain the grid velocities and write
C ... the resulting conservative variables

         EL = E + .5*(-RHOU**2 
     &          + (RHOU*A1X(IA0)*DIRIN + RHO*A1X(IA0)*FACEVEL)**2
     &          + (RHOU*A1Y(IA0)*DIRIN + RHO*A1Y(IA0)*FACEVEL)**2
     &          + (RHOU*A1Z(IA0)*DIRIN + RHO*A1Z(IA0)*FACEVEL)**2)/
     &            (RHO+1.E-20)

         BOUNR(NP) = RHO
         BOUNU(NP) = RHOU*A1X(IA0)*DIRIN + RHO*A1X(IA0)*FACEVEL
         BOUNV(NP) = RHOU*A1Y(IA0)*DIRIN + RHO*A1Y(IA0)*FACEVEL
         BOUNW(NP) = RHOU*A1Z(IA0)*DIRIN + RHO*A1Z(IA0)*FACEVEL
         BOUNE(NP) = EL

         BOUNMF(NP) = (A1X(IA0)*BOUNU(NP) + A1Y(IA0)*BOUNV(NP) +
     +                 A1Z(IA0)*BOUNW(NP) - FACEVEL)*BOUNR(NP)

         ELSE
         STOP

         ENDIF

         ENDDO
 400  CONTINUE
       
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN

C ... Turbulence parameters rho*k and rho*epsilon/rho*omega as
C ... a secondary set of results:

      DO 500 J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

         DO I = KX1,KX2

            NP = NR + I                   ! Patch index with LN ghost cells

            IF(TRANSL) THEN
               BOUNRK(NP)  = RK
               BOUNEP(NP)  = REPS
               BOUNG(NP)   = TG
               BOUNRET(NP) = TRET
               TU = SQRT(2./3.*BOUNRK(NP)/(BOUNR(NP)+EPS))/FRSVEL
               CALL INLETRET(TU,BOUNRET(NP))
            ELSE
               BOUNRK(NP) = RK
               BOUNEP(NP) = REPS
            ENDIF

         ENDDO 

 500  CONTINUE

      IF (ISTRES > 0) THEN ! b_ij distributions

         DO J = KY1,KY2

            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2

               NP = NR + I                ! Patch index with LN ghost cells

               BOUNBI(NP,1) = BIJ(1)
               BOUNBI(NP,2) = BIJ(2)
               BOUNBI(NP,3) = BIJ(3) 
               BOUNBI(NP,4) = BIJ(4) 
               BOUNBI(NP,5) = BIJ(5) 
               BOUNBI(NP,6) =-BOUNBI(NP,1) - BOUNBI(NP,4)

            ENDDO

         ENDDO

      ENDIF ! ISTRES > 0

      ENDIF ! (ITURB >= 3 .AND. ITURB /= 8)


C ... Multiphase variables

      IF(MULPHL) THEN

      DO J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2
               
               NP = NR + I                ! Patch index with LN ghost cells

               BOUNA1(NP) = RMUL1
               BOUNA2(NP) = RMUL2

            ENDDO

         ENDDO

      ENDIF ! MULPHL


C ... Scalar variables

      DO NS = 1,NSCAL

         DO J = KY1,KY2

            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2

               NP = NR + I                ! Patch index with LN ghost cells

               BOUNFI(NP,NS) = FISET(NS)

            ENDDO

         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE OUTPLG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OUTMOD(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNA1,
     +     BOUNA2,BOUNG,BOUNRET,MAXEB,MAXSB,
     +     RHO,RHOU,E,RK,REPS,FISET,A1,A1X,A1Y,A1Z,
     +     X,Y,Z,UROT,KX1,KX2,KY1,KY2,IX,IY,IZ,INRE,JNRE,KNRE,
     +     NSCAL,ITURB,TRANSL,MULPHL,NWALL,NGL,RGAS,GAMMA,
     +     P,TEMP,BIJ,ISTRES,IFILE,ISTATE,FRSDEN,FRSPRE,E0REF,T0REF,
     +     FRSVEL,RMUL1,RMUL2,TG,TRET,LRATUP,ARE,FRSVIS,FRSTEM,IUPPTL)

C ... The cell-wise boundary condition result file for FINFLO is 
C ... created. Approximate formulas for a developing tube flow 
C ... presented by Fred Stern of University of Iowa are applied
C ... to create a reasonably realistic flow distribution across
C ... the boundary surface instead of plug flow. The satic pressure,
C ... temperature and density are kept as constants.

C ... Either primitive or conservative primary variables are used
C ... as specified by IFILE. The final velocity or momentum components
C ... derived from the local relative contravariant momentum into the 
C ... grid are defined in the fixed inertial frame and thus contain
C ... the effects of grid rotation.

      USE CONSTANTS, ONLY : EPS
      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: IX, IY, IZ, IMAX, JMAX, KMAX, KX1, KX2, KY1, KY2
      INTEGER :: ISTR, JSTR, KSTR
      INTEGER :: NSCAL, ITURB, NWALL, NGL, NS, ISTRES, ISTATE 
      INTEGER :: IG0, NR, IA0, NP, I, J, K, IT1, IT2, IT3, IT4, IFILE
      INTEGER :: TYPEF, L, MAXEB, MAXSB, N
      INTEGER :: IN, JN, KN, INRE, JNRE, KNRE, IUPPTL

      REAL, DIMENSION(*) :: X, Y, Z
      REAL, DIMENSION(*) :: UROT
      REAL, DIMENSION(*) :: A1, A1X, A1Y, A1Z, FISET, BIJ

      REAL :: RHO, RHOU, E, EL, RK, REPS, RGAS, GAMMA, P, TEMP
      REAL :: XC, YC, ZC, ROIN, UVEL, VVEL, WVEL, FACEVEL, DIRIN,
     &        FRSDEN, FRSPRE, E0REF, T0REF, TU, FRSVEL, TG, TRET, FRSTEM

      REAL :: LRATUP, LRATCR, ARE, AEFF, DREP, RETUBE, FRSVIS 
      REAL :: DELTA1, DELT1X, RHOUC, NVTURB, WDIST, MFTOT, MFCOR
      REAL :: XW1, YW1, ZW1, WD1, XW2, YW2, ZW2, WD2, RMUL1, RMUL2 
      
      REAL :: BOUNMF(*),BOUNR(*),BOUNU(*),BOUNV(*),BOUNW(*),BOUNE(*),
     +     BOUNT(*),BOUNP(*),BOUNRK(*),BOUNEP(*),BOUNA1(*),
     +     BOUNA2(*),BOUNFI(MAXSB,*),BOUNBI(MAXEB,*),BOUNG(*),BOUNRET(*)

      REAL, ALLOCATABLE :: MFDSTR(:)

      LOGICAL :: MULPHL, TRANSL

      CHARACTER(LEN=12) :: NAME

C ... Determine the flow type based on the mean flow velocity,
C ... a representative tube diameter and the input LRATUP. Then 
C ... determine momentum RHOUC at the tube center for each case.

      DREP   = SQRT(4.*ARE/3.141592654)
      RETUBE = ABS(RHOU)*DREP/FRSVIS
      IF (RETUBE < 2300.) THEN
C ... Laminar flow
      LRATCR = 0.06*RETUBE
          IF (LRATUP > LRATCR) THEN
C ...     Fully developed Hagen-Poiseuille flow
          TYPEF = 1
          RHOUC = 2.*RHOU
          ELSE
          TYPEF = 2
C ...     Iterate the boundary-layer displacement thickness
          DELT1X = 0.
  20      AEFF   = 3.141592654/4.*(DREP - 2.*DELT1X)**2
          DELTA1 = 1.83*LRATUP*DREP/SQRT((ABS(RHOU)*LRATUP*DREP*ARE)/
     &             (FRSVIS*AEFF))
            IF (ABS(DELTA1 - DELT1X)/DREP > 0.00001) THEN
            DELT1X = DELTA1
            GOTO 20
            ENDIF
          RHOUC  = RHOU*ARE/AEFF
          ENDIF
      ELSEIF (RETUBE. GT. 4000.) THEN
C ... Turbulent flow
      NVTURB = 7.
      LRATCR = 4.4*RETUBE**0.166667
          IF (LRATUP > LRATCR) THEN
C ...     Fully developed turbulent tube flow
          TYPEF = 3
          RHOUC  = (1.+1./NVTURB)*(2.+1./NVTURB)/2.*RHOU 
          ELSE
          TYPEF = 4
C ...     Iterate the boundary-layer displacement thickness
          DELT1X = 0.
  40      AEFF   = 3.141592654/4.*(DREP - 2.*DELT1X)**2
          DELTA1 = 0.02*LRATUP*DREP/((ABS(RHOU)*LRATUP*DREP*ARE)/
     &             (FRSVIS*AEFF))**(1./NVTURB)
            IF (ABS(DELTA1 - DELT1X)/DREP > 0.00001) THEN
            DELT1X = DELTA1
            GOTO 40
            ENDIF
          RHOUC  = RHOU*ARE/AEFF 
          ENDIF
      ELSE
C ... Transitional flow
      NVTURB = 4.
      LRATCR = 0.06*RETUBE + (RETUBE - 2300.)/1700.*
     &         (4.4*RETUBE**0.166667 - 0.06*RETUBE)
          IF (LRATUP > LRATCR) THEN
C ...     "Fully developed transitional" (?) tube flow     
          TYPEF = 5
C ...     The velocity distribution is guessed to be something 
C ...     between the laminar and turbulent cases
          RHOUC  = (1.+1./NVTURB)*(2.+1./NVTURB)/2.*RHOU
          ELSE
          TYPEF = 6
C ...     The velocity distribution is simply guessed to be similar 
C ...     as in the laminar case!
          DELT1X = 0.
  60      AEFF   = 3.141592654/4.*(DREP - 2.*DELT1X)**2
          DELTA1 = 1.83*LRATUP*DREP/SQRT((ABS(RHOU)*LRATUP*DREP*ARE)/
     &             (FRSVIS*AEFF))
            IF (ABS(DELTA1 - DELT1X)/DREP > 0.00001) THEN
            DELT1X = DELTA1
            GOTO 60
            ENDIF
          RHOUC  = RHOU*ARE/AEFF
          ENDIF
      ENDIF

C ... Determine suitable indexing systems for the three wall types.
C ... STEPs are for the grid nodes that do not recognize ghost cells.
C ... STRs are for the cell face areas that are arranged to take the
C ... ghost cells into account as in AREAS and AREVC1.     

      IF     (NWALL == 1 .OR. NWALL == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         IMAX = IY
         JMAX = IZ
         KMAX = IX
         ISTR = KMAX+2*KN
         JSTR = (KMAX+2*KN)*(IMAX+2*IN)
         KSTR = 1
      ELSEIF (NWALL == 2 .OR. NWALL == 5)THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         IMAX = IZ
         JMAX = IX
         KMAX = IY
         ISTR = (JMAX+2*JN)*(KMAX+2*KN)
         JSTR = 1
         KSTR = JMAX+2*JN
      ELSEIF (NWALL == 3 .OR. NWALL == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         IMAX = IX
         JMAX = IY
         KMAX = IZ
         ISTR = 1
         JSTR = IMAX+2*IN
         KSTR = (IMAX+2*IN)*(JMAX+2*JN)

      ELSE
         CALL EXIT
      ENDIF

C ... Allocate temporary space for MFDSTR

      ALLOCATE(MFDSTR((IMAX+2*IN)*(JMAX+2*JN)))

C ... Initialization of the main variable

      DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
         BOUNR(N)  = 0.
         BOUNU(N)  = 0.
         BOUNV(N)  = 0.
         BOUNW(N)  = 0.
         BOUNE(N)  = 0.
         BOUNT(N)  = 0.
         BOUNP(N)  = 0.         
         BOUNMF(N) = 0.
         BOUNA1(N) = 0.
         BOUNA2(N) = 0.
      ENDDO

C ... Initialization of the turbulence variables

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNRK(N) = 0.
            BOUNEP(N) = 0.
         ENDDO
      ENDIF

C ... Initialization of the Intermittency variables

      IF (TRANSL) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNG(N)   = 0.0
            BOUNRET(N) = 0.0
         ENDDO
      ENDIF

C ... Initialization of EARSM anisotropic components

      IF(ISTRES > 0) THEN
         DO NS = 1,6
            DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
               BOUNBI(N,NS) = 0.
            ENDDO
         ENDDO
      ENDIF

C ... Initialization of the multiphase variables

      IF (MULPHL) THEN
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNA1(N) = 0.
            BOUNA2(N) = 0.
         ENDDO
      ENDIF

C ... Initialization of scalars

      DO NS = 1,NSCAL
         DO N = 1,(KX2-KX1+1+2*LN)*(KY2-KY1+1+2*LN)
            BOUNFI(N,NS) = 0.
         ENDDO
      ENDDO

C ... Select between index minima or maxima at the wall to define
C ... a global starting address for the cell face areas and
C ... define the direction into the grid
      K     = 1
      DIRIN = 1.
      IF (NWALL > 3) THEN
c         K     = KMAX      ! JH 2013
         K     = KMAX + 1  ! ESa 29.1.2016
         DIRIN = -1.
      ENDIF

      IG0   = 1 + (K-1+KN)*KSTR

C ... Determine a preliminary momentum distribution across the
C ... boundary surface according to the flow type and integrate the
C ... resulting mass flow

      MFTOT = 0.

      DO 100 J = 1,JMAX

         DO 100 I = 1,IMAX

         IA0 = IG0 + (I-1+IN)*ISTR + (J-1+JN)*JSTR

C ... For each boundary surface cell, find the shortest distance 
C ... WDIST to an edge taken to represent a solid tube wall 
         IT1 = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT2 = 1+(I  +IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT3 = 1+(I  +IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR
         IT4 = 1+(I-1+IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR
         XC  = .25*(X(IT1)+X(IT2)+X(IT3)+X(IT4))
         YC  = .25*(Y(IT1)+Y(IT2)+Y(IT3)+Y(IT4))
         ZC  = .25*(Z(IT1)+Z(IT2)+Z(IT3)+Z(IT4))
         WDIST = 1.E+6
         DO L = 1,JMAX
C ...    Two opposite edges in one loop
         IT1 = 1+(     IN)*ISTR+(L-1+JN)*JSTR+(K-1+KN)*KSTR
         IT2 = 1+(     IN)*ISTR+(L  +JN)*JSTR+(K-1+KN)*KSTR
         IT3 = 1+(IMAX+IN)*ISTR+(L-1+JN)*JSTR+(K-1+KN)*KSTR
         IT4 = 1+(IMAX+IN)*ISTR+(L  +JN)*JSTR+(K-1+KN)*KSTR
         XW1 = .5*(X(IT1)+X(IT2))
         YW1 = .5*(Y(IT1)+Y(IT2))
         ZW1 = .5*(Z(IT1)+Z(IT2))
         XW2 = .5*(X(IT3)+X(IT4))
         YW2 = .5*(Y(IT3)+Y(IT4))
         ZW2 = .5*(Z(IT3)+Z(IT4))
C ... The z=0 plane is assumed to be a symmetry plane. Therefore,
C ... the related distances are modified to avoid a boundary layer.
         WD1 = SQRT((XC-XW1)**2+(YC-YW1)**2
     &         +((ABS(ZW1)+1.E-6)/(ABS(ZW1)+1.E-8)*(ZC-ZW1))**2)
         WD2 = SQRT((XC-XW2)**2+(YC-YW2)**2
     &         +((ABS(ZW2)+1.E-6)/(ABS(ZW2)+1.E-8)*(ZC-ZW2))**2)
         WDIST = MIN(WDIST,WD1,WD2)
         ENDDO
         DO L = 1,IMAX
         IT1 = 1+(L-1+IN)*ISTR+(     JN)*JSTR+(K-1+KN)*KSTR
         IT2 = 1+(L  +IN)*ISTR+(     JN)*JSTR+(K-1+KN)*KSTR
         IT3 = 1+(L  +IN)*ISTR+(JMAX+JN)*JSTR+(K-1+KN)*KSTR
         IT4 = 1+(L-1+IN)*ISTR+(JMAX+JN)*JSTR+(K-1+KN)*KSTR
         XW1 = .5*(X(IT1)+X(IT2))
         YW1 = .5*(Y(IT1)+Y(IT2))
         ZW1 = .5*(Z(IT1)+Z(IT2))
         XW2 = .5*(X(IT3)+X(IT4))
         YW2 = .5*(Y(IT3)+Y(IT4))
         ZW2 = .5*(Z(IT3)+Z(IT4))
         WD1 = SQRT((XC-XW1)**2+(YC-YW1)**2
     &         +((ABS(ZW1)+1.E-6)/(ABS(ZW1)+1.E-8)*(ZC-ZW1))**2)
         WD2 = SQRT((XC-XW2)**2+(YC-YW2)**2
     &         +((ABS(ZW2)+1.E-6)/(ABS(ZW2)+1.E-8)*(ZC-ZW2))**2)
         WDIST = MIN(WDIST,WD1,WD2)
         ENDDO

C ... Local flow momenta as functions of the wall distance:
         IF (TYPEF == 1) THEN
C ... Approximated Hagen-Poiseuille:
         MFDSTR((J-1)*IMAX+I) = RHOUC*(1. - ((DREP-2.*WDIST)/DREP)**2)
         ELSEIF (TYPEF == 2 .OR. TYPEF == 6) THEN
C ... Laminar boundary layer of thickness 3*DELTA1:
           IF (WDIST < 3.*DELTA1) THEN
           MFDSTR((J-1)*IMAX+I) = RHOUC*(2.*WDIST/(3.*DELTA1) 
     &                            - (WDIST/(3.*DELTA1))**2)
           ELSE
C ... Constant core flow:
           MFDSTR((J-1)*IMAX+I) = RHOUC
           ENDIF
         ELSEIF (TYPEF == 3 .OR. TYPEF == 5) THEN
C ... Approximated fully-developed turbulent or transitional tube low: 
         MFDSTR((J-1)*IMAX+I) = RHOUC*(1. - ((DREP-2.*WDIST)
     &                          /DREP))**(1./NVTURB)
         ELSEIF (TYPEF == 4) THEN
C ... Turbulent boundary layer of thickness 8*DELTA1:
           IF (WDIST < 8.*DELTA1) THEN
           MFDSTR((J-1)*IMAX+I) = RHOUC*(WDIST
     &                            /(8.*DELTA1))**(1./NVTURB) 
           ELSE
C ... Constant core flow:
           MFDSTR((J-1)*IMAX+I) = RHOUC
           ENDIF
         ENDIF
      MFTOT = MFTOT + A1(IA0)*MFDSTR((J-1)*IMAX+I)
 100  CONTINUE

C ... The ratio ARE*RHOU/MFTOT will be applied to adjust the flow  
C ... distribution to match the specified input mass flow
      MFCOR = ARE*RHOU/MFTOT

C ... Automatic naming of the output file to contain 
C ... block and wall numbers

      NAME( 1:12) = 'DIST_B000_F0'

      IF(NGL < 10) THEN
         NAME( 9:9 ) = CHAR(NGL+48)
      ELSEIF(NGL >= 10 .AND. NGL < 100) THEN
         NAME( 8:8 ) = CHAR(NGL/10+48)
         NAME( 9:9 ) = CHAR(MOD(NGL,10)+48)
      ELSE
         NAME( 7:7 ) = CHAR(NGL/100+48)
         NAME( 8:8 ) = CHAR(MOD(NGL,100)/10+48)
         NAME( 9:9 ) = CHAR(MOD(NGL,10)+48)
      ENDIF

      NAME(12:12) = CHAR(48+NWALL)

C ... Loop for the cell-wise output of the main flow variables

      DO 400 J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

         DO I = KX1,KX2

         IA0 = IG0 + (I-1+IN)*ISTR + (J-1+JN)*JSTR
         NP  = NR + I                     ! Patch index with LN ghost cells

C ... Compute the velocity components resulting from the grid 
C ... rotation

         IT1 = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT2 = 1+(I  +IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         IT3 = 1+(I  +IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR
         IT4 = 1+(I-1+IN)*ISTR+(J  +JN)*JSTR+(K-1+KN)*KSTR

         FACEVEL = UROT(IA0)

         IF(IFILE == 0) THEN

C ... Add the grid velocities to the relative velocities at the 
C ... boundary to obtain absolute inertial velocities before writing
C ... the primitive variables  

         UVEL = MFCOR*MFDSTR((J-1)*IMAX+I)*A1X(IA0)
     &          /(RHO+1.E-20)*DIRIN + A1X(IA0)*FACEVEL
         VVEL = MFCOR*MFDSTR((J-1)*IMAX+I)*A1Y(IA0)
     &          /(RHO+1.E-20)*DIRIN + A1Y(IA0)*FACEVEL
         WVEL = MFCOR*MFDSTR((J-1)*IMAX+I)*A1Z(IA0)
     &          /(RHO+1.E-20)*DIRIN + A1Z(IA0)*FACEVEL

         BOUNP(NP) = P
         BOUNU(NP) = UVEL
         BOUNV(NP) = VVEL
         BOUNW(NP) = WVEL
         BOUNT(NP) = TEMP
         CALL ROFPT(TEMP,P,ROIN,1,ISTATE,RGAS,
     +        GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
         BOUNR(NP) = ROIN

         BOUNMF(NP) = (A1X(IA0)*BOUNU(NP) + A1Y(IA0)*BOUNV(NP) +
     +                 A1Z(IA0)*BOUNW(NP) - FACEVEL)*BOUNR(NP)

         ELSE IF(IFILE == 1) THEN

C ... Modify the total energy to model the velocity ditribution and
C ... to contain the grid velocities before writing the resulting 
C ... conservative variables

         EL = E + .5*(-RHOU**2 
     &          + (MFCOR*MFDSTR((J-1)*IMAX+I)*A1X(IA0)*DIRIN
     &          +  RHO*A1X(IA0)*FACEVEL)**2
     &          + (MFCOR*MFDSTR((J-1)*IMAX+I)*A1Y(IA0)*DIRIN
     &          +  RHO*A1Y(IA0)*FACEVEL)**2
     &          + (MFCOR*MFDSTR((J-1)*IMAX+I)*A1Z(IA0)*DIRIN
     &          +  RHO*A1Z(IA0)*FACEVEL)**2)/(RHO+1.E-20)

         BOUNR(NP) = RHO
         BOUNU(NP) = MFCOR*MFDSTR((J-1)*IMAX+I)*A1X(IA0)*DIRIN
     &             + RHO*A1X(IA0)*FACEVEL
         BOUNV(NP) = MFCOR*MFDSTR((J-1)*IMAX+I)*A1Y(IA0)*DIRIN
     &             + RHO*A1Y(IA0)*FACEVEL
         BOUNW(NP) = MFCOR*MFDSTR((J-1)*IMAX+I)*A1Z(IA0)*DIRIN
     &             + RHO*A1Z(IA0)*FACEVEL
         BOUNE(NP) = EL

         BOUNMF(NP) = (A1X(IA0)*BOUNU(NP) + A1Y(IA0)*BOUNV(NP) +
     +                 A1Z(IA0)*BOUNW(NP) - FACEVEL)*BOUNR(NP)

         ELSE
         STOP
         ENDIF
         ENDDO
 400  CONTINUE
       

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN

C ... Turbulence parameters rho*k and rho*epsilon/rho*omega as
C ... a secondary set of results:

      DO 500 J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

         DO I = KX1,KX2

            NP = NR + I                   ! Patch index with LN ghost cells

            IF (RETUBE < 3000.) THEN

C ... Essentially laminar flow with suppressed turbulence energy

               BOUNRK(NP) = 0.01*RK
               BOUNEP(NP) = REPS

            ELSE

C ... Modify RK by the velocity distribution to obtain at least a 
C ... loose connection to velocity gradients. Leave REPS as it is

               BOUNRK(NP) = 1.7*(MFDSTR((J-1)*IMAX+I)/RHOUC)**0.333
     &              *(1.02-MFDSTR((J-1)*IMAX+I)/RHOUC)**0.333*RK
               BOUNEP(NP) = REPS

            ENDIF

            IF(TRANSL) THEN
               BOUNG(NP)   = TG
               BOUNRET(NP) = TRET
               TU = SQRT(2./3.*BOUNRK(NP)/(BOUNR(NP)+EPS))/FRSVEL
               CALL INLETRET(TU,BOUNRET(NP))
            ELSE
            ENDIF

         ENDDO 

 500  CONTINUE

      IF (ISTRES > 0) THEN ! b_ij distributions

         DO J = KY1,KY2

            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2

               NP = NR + I                ! Patch index with LN ghost cells

               BOUNBI(NP,1) = BIJ(1)
               BOUNBI(NP,2) = BIJ(2)
               BOUNBI(NP,3) = BIJ(3) 
               BOUNBI(NP,4) = BIJ(4) 
               BOUNBI(NP,5) = BIJ(5) 
               BOUNBI(NP,6) =-BOUNBI(NP,1) - BOUNBI(NP,4)


            ENDDO

         ENDDO

      ENDIF ! ISTRES > 0

      ENDIF ! (ITURB >= 3 .AND. ITURB /= 8)


C ... Multiphase varibales

      IF(MULPHL) THEN

      DO J = KY1,KY2

         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2
               
               NP = NR + I                ! Patch index with LN ghost cells

               BOUNA1(NP) = RMUL1
               BOUNA2(NP) = RMUL2

            ENDDO

         ENDDO

      ENDIF ! MULPHL


C ... Scalar variables

      DO NS = 1,NSCAL

         DO J = KY1,KY2

            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1

            DO I = KX1,KX2

               NP = NR + I                ! Patch index with LN ghost cells

               BOUNFI(NP,NS) = FISET(NS)

            ENDDO

         ENDDO

      ENDDO


      DEALLOCATE(MFDSTR)

      RETURN
      END SUBROUTINE OUTMOD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
