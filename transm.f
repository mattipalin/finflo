C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRANSM(TRM,IMAX,JMAX,KMAX,RO,VIS,DISTW,
     &             RK,REPS,U,V,W,DUIDXJ,STRAIN,OHMI,VOL,
     &             FRSVEL,NTOT,NGL)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     The subroutine calculates the transition model source 
C     terms. Transition model is based on Langtry PhD thesis 
C     (2006), see Appendix C in Finflo manual. 
C
C     Created 14.09.2012 by Juho Nikulainen
C     ---------------------------------------------------------

      USE CONSTANTS, ONLY : EPS6, EPS8 
      USE NS3CO, ONLY     : IN, JN, KN
      USE TYPE_ARRAYS 

      IMPLICIT NONE

      INTEGER :: ISTRID,JSTRID,KA,K,KK,J,JJ,I,L,
     &           IMAX,JMAX,KMAX,NTOT,NGL

      REAL, INTENT(IN) :: U(*),V(*),W(*),DUIDXJ(NTOT,9),OHMI(*),
     &                    STRAIN(*),RO(*),VIS(*),DISTW(*),RK(*),
     &                    REPS(*),VOL(*)
      REAL, INTENT(IN) :: FRSVEL

      REAL, ALLOCATABLE :: REOT(:),UT(:),TU(:),REOC(:),FLENG(:)

      REAL :: CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1,
     &        REV,THTBL,DBL,DD,REW,FWAKE,FOT,T,POT,
     &        FON,FON1,FON2,FON3,RT,PG1,EG1,FTURB,PG2,EG2,PEG,
     &        FREAT,GSEP

      TYPE(INTERMITTENCY) TRM(*)

      ALLOCATE (REOT(NTOT),UT(NTOT),TU(NTOT),REOC(NTOT),FLENG(NTOT))

C     %%%%%%% EMPIRICAL TERMS & INITIALIZATION
      CALL TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
      CALL SREOT(U,V,W,VIS,RO,RK,DUIDXJ,FRSVEL,
     &            IMAX,JMAX,KMAX,REOT,UT,TU,NTOT)
      CALL SREOC(TRM,NTOT,REOC)
      CALL SFLENG(TRM,NTOT,FLENG)

      TRM(1:NTOT)%SG   = 0.
      TRM(1:NTOT)%ZG   = 0.
      TRM(1:NTOT)%LG   = 0.
      TRM(1:NTOT)%SRET = 0.
      TRM(1:NTOT)%ZRET = 0.
      TRM(1:NTOT)%LRET = 0.      
      

C     %%%%%%% START LOOPS
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KA     = (KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID + IN
      DO K = 1, KMAX
         KK = K*ISTRID*JSTRID + KA            
         DO J = 1, JMAX
            JJ = J*ISTRID + KK
            DO I = 1, IMAX
               L = JJ + I 

C ... TRANSITION ONSET MOMENTUM THICKNESS TRANSPORT EQ
      THTBL    = TRM(L)%RET*VIS(L)/(RO(L)*UT(L))
      DBL      = 15.*THTBL/2.
      DD       = 50.* OHMI(L)*DISTW(L)*DBL/UT(L)
      REW      = REPS(L)*DISTW(L)**2/(VIS(L)+EPS8)
      FWAKE    = EXP(-(REW/(1.E5))**2)
      FOT      = MIN( MAX( FWAKE*EXP(-(DISTW(L)/DD)**4),
     &             1.-( (TRM(L)%G-1./CE2) / (1.-1./CE2) )**2 ), 1.)
      T        = 500.*VIS(L)/(RO(L)*UT(L)**2)

      TRM(L)%SRET =  COT*RO(L)*REOT(L)   *(1.-FOT)/T   !source
      TRM(L)%ZRET = -COT*RO(L)*TRM(L)%RET*(1.-FOT)/T   !sink
      TRM(L)%LRET = TRM(L)%ZRET / (TRM(L)%RET+EPS8)    !dZG/dRET
C     Total influence
      POT         = COT*RO(L)*(REOT(L)-TRM(L)%RET)*(1.-FOT)/T

C ... STRAIN RATE RE
      REV      = RO(L)*DISTW(L)**2*STRAIN(L)/(VIS(L)+EPS8)

C ... INTERMITTENCY TRANSPORT EQ 
      FON1     = REV/(2.193*REOC(L)+EPS8)
      FON2     = MIN(MAX( FON1, FON1**4 ), 2.)
      RT       = RO(L)*RK(L) / (VIS(L)*REPS(L)+EPS8) 
      FON3     = MAX( 1.-(RT/2.5)**3, 0. )
      FON      = MAX( FON2-FON3, EPS8 )
      FTURB    = EXP(-(RT/4.)**4)    
      PG1      =  CA1*FLENG(L)*RO(L)*STRAIN(L)*SQRT(TRM(L)%G*FON)
      EG1      = -CE1*PG1*TRM(L)%G
      PG2      =  CA2*RO(L)*OHMI(L)*TRM(L)%G*FTURB 
      EG2      = -CE2*PG2*TRM(L)%G

      TRM(L)%SG = PG1 + PG2                  !source
      TRM(L)%ZG = EG1 + EG2                  !sink
      TRM(L)%LG = -3./2.*CE1*PG1 -2.*CE2*PG2 !dZG/dG
C     Total influence
      PEG       = PG1 + EG1 + PG2 + EG2

C ... SEPARATION CHECK FOR EFFECTIVE INTERMITTENCY
      FREAT     = EXP(-(RT/20.)**4)
      GSEP      = MIN(S1*FREAT*MAX(0., 
     &                REV/(3.235*REOC(L)+EPS8)-1.), 2.)*FOT
*      TRM(L)%G  = MAX(TRM(L)%G,GSEP)
      TRM(L)%GEFF  = MAX(TRM(L)%G,GSEP)

C ... ADD SOURCES TO RESIDUALS
coli      TRM(L)%DG   = TRM(L)%DG  + PEG     
coli      TRM(L)%DRET = TRM(L)%DRET+ POT

            ENDDO
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE TRANSM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Transition model coefficients.
C     ---------------------------------------------------------

      IMPLICIT NONE 
      REAL, INTENT(OUT) :: CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1

C     RET transport equation
      COT      = 0.03
      SIGOT    = 2.0      !10. Menter
C     Intermittency transport equation 
      CE1      = 1.0
      CE2      = 50.
      CA1      = 2.0      !0.5 Menter
      CA2      = 0.06     !0.03 Menter
      SIGF     = 1.0
C     Other
      S1       = 2.0      !8.0 Menter 

      RETURN
      END SUBROUTINE TRANVAR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SREOT(U,V,W,VIS,RO,RK,DUIDXJ,FRSVEL,
     &                 IMAX,JMAX,KMAX,REOT,UT,TU,NTOT)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     A subroutine to calculate the momentum thickness
C     Reynolds number at transition onset (REOT),
C     local turbulence intensity (TU) and total velocity (UT).       
C     ---------------------------------------------------------

      USE CONSTANTS, ONLY : EPS6 
      USE NS3CO, ONLY     : IN, JN, KN
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: II,IMAX,JMAX,KMAX,
     &           ISTRID,JSTRID,KA,K,KK,J,JJ,I,L,NTOT
      REAL, INTENT(IN) :: U(*),V(*),W(*),VIS(*),RO(*),RK(*),
     &                    DUIDXJ(NTOT,9)
      REAL, INTENT(IN) :: FRSVEL      
      REAL :: Uxt,Uyt,Uzt,ux,uy,uz,vx,vy,vz,wx,wy,wz,DUDS,
     &        LO,FL,THETA,DIFF,REOTOLD
      REAL, INTENT(OUT) :: REOT(NTOT),UT(NTOT),TU(NTOT)      

      REOT(1:NTOT)  = 0.
      UT(1:NTOT)    = 0.
      TU(1:NTOT)    = 0.


C     %%%%%%% START LOOPS
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KA     = (KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID + IN
      DO K = 1, KMAX
         KK = K*ISTRID*JSTRID + KA            
         DO J = 1, JMAX
            JJ = J*ISTRID + KK
            DO I = 1, IMAX
               L = JJ + I 
      

      UT(L)      = SQRT(U(L)**2+V(L)**2+W(L)**2+EPS6)
      TU(L)      = 100.*SQRT(2./3.*RK(L))/UT(L)
      IF (TU(L) < 0.027) THEN
            TU(L) = 0.027
      ENDIF

C     Derivatives of velocity components: 
      ux = DUIDXJ(L,1); uy = DUIDXJ(L,2); uz = DUIDXJ(L,3)
      vx = DUIDXJ(L,4); vy = DUIDXJ(L,5); vz = DUIDXJ(L,6)
      wx = DUIDXJ(L,7); wy = DUIDXJ(L,8); wz = DUIDXJ(L,9)
      
      Uxt      = 0.5*(2.*U(L)*ux+2.*V(L)*vx+2.*W(L)*wx)/UT(L)
      Uyt      = 0.5*(2.*U(L)*uy+2.*V(L)*vy+2.*W(L)*wy)/UT(L)
      Uzt      = 0.5*(2.*U(L)*uz+2.*V(L)*vz+2.*W(L)*wz)/UT(L)
      DUDS     = U(L)/UT(L)*Uxt + V(L)/UT(L)*Uyt + W(L)/UT(L)*Uzt

C     Initial guess
      IF (TU(L) <= 1.3) THEN 
        REOT(L) = 1173.51-589.428*TU(L)+0.2196/TU(L)**2
      ELSE
        REOT(L) = 331.50*(TU(L)-0.5658)**(-0.671)
      ENDIF            
      THETA    = VIS(L)*REOT(L)/(RO(L)*FRSVEL+EPS6)

C	%%% Jacobi method for solving REOT
      II      = 1
      DIFF    = 1.
      DO WHILE (DIFF > 0.01)
      II      = II + 1
      IF (II == 20) EXIT

      REOTOLD  = REOT(L)
      LO       = THETA**2*DUDS*RO(L)/(VIS(L)+EPS6)   
      IF (LO < -0.1) THEN
            LO = -0.1
      ELSEIF (LO > 0.1) THEN
            LO = 0.1
      ENDIF
      
      IF (LO <= 0) THEN
            FL = 1.-(-12.986*LO-123.66*LO**2-405.689*LO**3)
     &               *EXP(-(TU(L)/1.5)**1.5)
      ELSE
            FL = 1.+0.275*(1.-EXP(-35.*LO))*EXP(-TU(L)/0.5)
      ENDIF      

      IF (TU(L) <= 1.3) THEN 
        REOT(L) = (1173.51-589.428*TU(L)+0.2196/TU(L)**2)*FL
      ELSE
        REOT(L) = 331.50*(TU(L)-0.5658)**(-0.671)*FL
      ENDIF

      DIFF     = ABS( (REOT(L)-REOTOLD)/REOT(L) )
      ENDDO !while loop


C     %%%%%%% Limit REOT
      IF (REOT(L) < 20.) THEN 
            REOT(L) = 20.
      ENDIF 

            ENDDO
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE !SREOT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SREOC(TRM,NTOT,REOC)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     A subroutine to calculate the momentum thickness critical
C     Reynolds number (REOC) for the transition model. 
C     ---------------------------------------------------------

      USE TYPE_ARRAYS

      IMPLICIT NONE
      INTEGER :: NTOT
      REAL, INTENT(OUT) :: REOC(NTOT)      
      TYPE(INTERMITTENCY) TRM(*)

      REOC(1:NTOT) = 
     &   MIN(0.625*TRM(1:NTOT)%RET+62.,TRM(1:NTOT)%RET)


      RETURN
      END SUBROUTINE !SREOC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SFLENG(TRM,NTOT,FLENG)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     A subroutine to calculate the length of the transition
C     zone (FLENG) for the transition model.  
C     ---------------------------------------------------------

      USE TYPE_ARRAYS
     
      IMPLICIT NONE
      INTEGER :: NTOT
      REAL, INTENT(OUT) :: FLENG(NTOT)
      TYPE(INTERMITTENCY) TRM(*)
      
      FLENG(1:NTOT) = 
     &   MIN(0.01*EXP(-0.022*TRM(1:NTOT)%RET+12.)+0.57,300.)


      RETURN
      END SUBROUTINE !SFLENG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE INLETRET(FRSTUR,RETIN)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     A subroutine to calculate the inlet value of RET
C     for the transition model.  
C     ---------------------------------------------------------
     
      IMPLICIT NONE
      REAL :: TU
      REAL, INTENT(IN) :: FRSTUR
      REAL, INTENT(OUT) :: RETIN

	TU	= FRSTUR*100.   !Tu in percents
	IF (TU < 0.027) THEN
            TU = 0.027
	ENDIF

	IF (TU <= 1.3) THEN 
            RETIN = 1173.51-589.428*TU+0.2196/TU**2
	ELSE
            RETIN = 331.50*(TU-0.5658)**(-0.671)
	ENDIF	

	IF (RETIN < 20.) THEN	
            RETIN = 20.
	ENDIF      


      RETURN
      END SUBROUTINE INLETRET
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C





