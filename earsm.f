      SUBROUTINE EARSM1(VIS,VIST,EPS2,EPS4,RK,WIR,WIJ,SIJ,BIJ,
     &     IMAX,JMAX,KMAX,MAXEB,MAX11,TTS,PR,PRT,TURLIM,FRSMUT)
C     *****************

C     A subroutine to calculate the extra anisotropy components and
C     the effective eddy-viscosity coefficient used with the EARSM
C                       23.11.2001 Ville H. / Laboratory of Aerodynamics

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: IMAX,JMAX,KMAX,MAX11,MAXEB
      INTEGER :: ISTRID,JSTRID,IL,KA,K,KK,J,JJ,I,L

      REAL , INTENT(IN) , DIMENSION(*) :: VIS,RK,TTS
      REAL , DIMENSION(*) :: EPS2,VIST
      REAL , INTENT(IN) , DIMENSION(MAX11,6) :: SIJ
      REAL , INTENT(IN) , DIMENSION(MAX11,3) :: WIJ
      REAL , DIMENSION(MAXEB,*) :: BIJ 
      REAL , DIMENSION(MAXEB,3) :: WIR

      REAL :: P3,P6,TT,C1,C1P,C1P3,C1PSQ,PRS,EPS2UP,EPS4UP, EPS4
      REAL :: S11,S22,S33,S12,S13,S23,W12,W13,W23
      REAL :: S11SQ,S22SQ,S33SQ,S12SQ,S13SQ,S23SQ,W12SQ,W13SQ,W23SQ
      REAL :: HTTS,SII,WII,WIIP3,SWWIV,SWWIVTT,SSWWV
      REAL :: TERM3C11,TERM3C12,TERM3C13,TERM3C22,TERM3C23 
      REAL :: TERM4C11,TERM4C12,TERM4C13,TERM4C22,TERM4C23 
      REAL :: TERM6C11,TERM6C12,TERM6C13,TERM6C22,TERM6C23 
      REAL :: TERM9C11,TERM9C12,TERM9C13,TERM9C22,TERM9C23 
      REAL :: P1,P2,SQP2,PM,PMP,SQPM,FACOS,FNC,D,FII1,FII2,FN,FNSQ,Q,PQ
      REAL :: BETA1,BETA3,BETA4,BETA6,BETA9,CMUEFF,PVIS,EPS2LO,EPS4LO
      REAL :: A0,PA0,ALPO
      REAL :: PR,PRT,TURLIM,FRSMUT

C     Often used numbers
      P3     = 1.0/3.0          ! One third
      P6     = 0.5*P3           ! One sixth
      TT     = 2.0*P3           ! Two thirds
      
C     Model coefficient and its variants
      C1     = 1.8
      C1P    = 9.0*(C1 - 1.0)/4.0 ! C1_prime
      C1P3   = C1P*P3        
      C1PSQ  = C1P**2
      A0     = 0.72               ! less the minus sign.
*      PA0    = 1.0/A0
      PA0    = 0. ! is not working yet, WIR has no value

C     Upper limits of EPS2 and EPS4
      PRS    = PR/PRT     
      EPS2UP = 1.0 + TURLIM
      EPS4UP = 1.0 + TURLIM*PRS

C     Constants for the loops over the grid
      ISTRID = IMAX + 2*IN      ! Total number of cells in i-direction
      JSTRID = JMAX + 2*JN      ! Total number of cells in j-direction
      IL     = ISTRID*JSTRID    ! Total number of cells in ij-plane
      KA     = (KN-1)*IL + (JN-1)*ISTRID + IN ! The start address
      
C     Actual loops
      DO K = 1, KMAX
         KK = K*IL + KA            
         DO J = 1, JMAX
            JJ = J*ISTRID + KK
            DO I = 1, IMAX
               L = JJ + I 

C     Half of the time scale
               HTTS   = 0.5*TTS(L) ! unused
C     Strain-rate and vorticity tensor components
               S11    = TTS(L)*SIJ(L,1)
               S22    = TTS(L)*SIJ(L,4)
               S33    = TTS(L)*SIJ(L,6)

               S12    = TTS(L)*SIJ(L,2)
               S13    = TTS(L)*SIJ(L,3)
               S23    = TTS(L)*SIJ(L,5)

c               S12    = HTTS*(DUIDXJ(L,2) + DUIDXJ(L,4))
c               S13    = HTTS*(DUIDXJ(L,3) + DUIDXJ(L,7))
c               S23    = HTTS*(DUIDXJ(L,6) + DUIDXJ(L,8))

C     The last term is the curvature term approximating the convection 
C     of a_ij.  
c               W12    = HTTS*(DUIDXJ(L,2) - DUIDXJ(L,4)) - PA0*WIR(L,3)
c               W13    = HTTS*(DUIDXJ(L,3) - DUIDXJ(L,7)) + PA0*WIR(L,2)
c               W23    = HTTS*(DUIDXJ(L,6) - DUIDXJ(L,8)) - PA0*WIR(L,1)

               W12    = TTS(L)*WIJ(L,1) - PA0*WIR(L,3)
               W13    = TTS(L)*WIJ(L,2) + PA0*WIR(L,2)
               W23    = TTS(L)*WIJ(L,3) - PA0*WIR(L,1)

C     Squares of strain-rate and vorticity tensor components
               S11SQ  = S11*S11
               S22SQ  = S22*S22
               S33SQ  = S33*S33

               S12SQ  = S12*S12
               S13SQ  = S13*S13
               S23SQ  = S23*S23

               W12SQ  = W12*W12
               W13SQ  = W13*W13
               W23SQ  = W23*W23

C     Second invariants of the strain rate and vorticity tensors
               SII    = S11SQ + S22SQ + S33SQ + 2.0*(S12SQ+S13SQ+S23SQ)
               WII    =-2.0*(W12SQ + W13SQ + W23SQ)
               WIIP3  = P3*WII  ! One third of the invariant

C     Third invanriant of the strain rate and vorticity tensors:
               SWWIV  =-S11*(W12SQ + W13SQ) - S22*(W12SQ + W23SQ)
     &                - S33*(W13SQ + W23SQ)
     &                + 2.0*(-S12*W13*W23 + S13*W12*W23 - S23*W12*W13)
               SWWIVTT= TT*SWWIV ! Two thirds of the invariant

C     Fourth invariant of the strain rate and vorticity tensors
               SSWWV  = 2.0*(-(S12*S13 + S22*S23 + S23*S33)*W12*W13
     &                       +(S11*S13 + S12*S23 + S13*S33)*W12*W23
     &                       -(S11*S12 + S12*S22 + S13*S23)*W13*W23)
     &              - (S11SQ + S12SQ + S13SQ)*(W12SQ + W13SQ)
     &              - (S12SQ + S22SQ + S23SQ)*(W12SQ + W23SQ)
     &              - (S13SQ + S23SQ + S33SQ)*(W13SQ + W23SQ)

C     Tensor component terms for beta 3
               TERM3C11 = - W12SQ - W13SQ - WIIP3
               TERM3C22 = - W12SQ - W23SQ - WIIP3
               TERM3C12 = - W13*W23
               TERM3C13 =   W12*W23
               TERM3C23 = - W12*W13

C     Tensor component terms for beta 4
               TERM4C11 =-2.0*(S12*W12 + S13*W13)
               TERM4C22 = 2.0*(S12*W12 - S23*W23)
               TERM4C12 = (S11-S22)*W12       - S23*W13       - S13*W23
               TERM4C13 =     - S23*W12 + (S11-S33)*W13       + S12*W23
               TERM4C23 =       S13*W12       + S12*W13 + (S22-S33)*W23

c     tensor component terms for beta 6
               TERM6C11 = -2.0*((S12*W13 - S13*W12)*W23
     &                        + S11*(W12SQ + W13SQ)) - SWWIVTT - WII*S11
               TERM6C22 = -2.0*((S23*W12 + S12*W23)*W13
     &                        + S22*(W12SQ + W23SQ)) - SWWIVTT - WII*S22
               TERM6C12 = -S12*(2.0*W12SQ + W13SQ + W23SQ)
     &             - (S13*W13-S23*W23)*W12 - (S11+S22)*W13*W23 - WII*S12
               TERM6C13 = -S13*(W12SQ + 2.0*W13SQ + W23SQ)
     &             - (S12*W12+S23*W23)*W13 + (S11+S33)*W12*W23 - WII*S13
               TERM6C23 = -S23*(W12SQ + W13SQ + 2.0*W23SQ)
     &             + (S12*W12-S13*W13)*W23 - (S22+S33)*W12*W13 - WII*S23

C     Tensor component terms for beta 9
               TERM9C11 =-2.0*(( S12*W12 + S13*W13 - S23*W23)*W12SQ
     &                        +( S12*W12 + S13*W13 + S23*W23)*W13SQ
     &                        +( S22-S33)*W12*W13*W23)
               TERM9C22 =-2.0*((-S12*W12 - S13*W13 + S23*W23)*W12SQ
     &                        +(-S12*W12 + S13*W13 + S23*W23)*W23SQ
     &                        +(-S11+S33)*W12*W13*W23)
               TERM9C12 = ((S11-S22)*W12 - 2.0*(S13*W23+S23*W13))*W12SQ
     &                  + ((S11-S33)*W12 - 2.0* S13*W23)         *W13SQ
     &                  + ((S33-S22)*W12 - 2.0* S23*W13)         *W23SQ
               TERM9C13 = ((S11-S22)*W13 + 2.0* S12*W23)         *W12SQ
     &                  + ((S11-S33)*W13 + 2.0*(S12*W23-S23*W12))*W13SQ
     &                  + ((S22-S33)*W13 - 2.0* S23*W12)         *W23SQ
               TERM9C23 = ((S22-S11)*W23 + 2.0* S12*W13)         *W12SQ
     &                  + ((S11-S33)*W23 + 2.0* S13*W12)         *W13SQ
     &                  + ((S22-S33)*W23 + 2.0*(S12*W13+S13*W12))*W23SQ

C     Solution of the third degree equation for N_c
               P1        = (C1PSQ/27.0 + 0.45*SII - TT*WII)*C1P
               P2        = P1**2 - (C1PSQ*P3**2 + 0.9*SII + TT*WII)**3
               IF (P2 >= 0.0) THEN
                  SQP2   = SQRT(P2)
                  PM     = P1 - SQP2
                  PMP    = ABS(PM)**P3 
                  FNC    = C1P3 + (P1+SQP2)**P3 + SIGN(PMP,PM)
               ELSE
                  PM     = P1**2 - P2
                  SQPM   = SQRT(PM)
                  ALPO   = P1/SQPM  ! Required by optimization with Ifort10
                  ALPO   = MIN(ALPO,1.)
                  ALPO   = MAX(ALPO,-1.)
                  FACOS  = P3*ACOS(ALPO) 
                  FNC    = C1P3 + 2.0*(PM**P6)*COS(FACOS) 
               END IF

C     Improvement of the approximation of N
               D         = 20.0*(FNC**4)*(FNC - 0.5*C1P)
     &                     - WII*(10.0*FNC + 15.0*C1P)*FNC**2
     &                     + 10.0*C1P*WII**2
               FII1      = SWWIV**2
               FII2      = SSWWV - 0.5*SII*WII
               FN        = FNC + 162.0*(FII1 + FII2*FNC**2)/D

C     The denominator of the betas
               FNSQ      = FN**2
               Q         = 5.0*P6*(FNSQ - 2.0*WII)*(2.0*FNSQ - WII)
               PQ        = 1.0/Q     

C     Coefficients (betas)
               BETA1     = - PQ*FN*(2.0*FNSQ - 7.0*WII)
               BETA3     = - PQ*12.0*SWWIV/FN
               BETA4     = - PQ*2.0*(FNSQ - 2.0*WII)
               BETA6     = - PQ*6.0*FN
               BETA9     =   PQ*6.0
               
C     Extra anisotropy components.  Note that we use b_ij and Wallin
C     uses a_ij, which is a_ij = 2*b_ij. 
               BIJ(L,1) = 0.5*(BETA3*TERM3C11 + BETA4*TERM4C11
     &              + BETA6*TERM6C11 + BETA9*TERM9C11)
               BIJ(L,2) = 0.5*(BETA3*TERM3C12 + BETA4*TERM4C12
     &              + BETA6*TERM6C12 + BETA9*TERM9C12)
               BIJ(L,3) = 0.5*(BETA3*TERM3C13 + BETA4*TERM4C13
     &              + BETA6*TERM6C13 + BETA9*TERM9C13)
               BIJ(L,4) = 0.5*(BETA3*TERM3C22 + BETA4*TERM4C22
     &              + BETA6*TERM6C22 + BETA9*TERM9C22)
               BIJ(L,5) = 0.5*(BETA3*TERM3C23 + BETA4*TERM4C23
     &              + BETA6*TERM6C23 + BETA9*TERM9C23)
               BIJ(L,6) = -BIJ(L,1) - BIJ(L,4)

        if(imax <  0) then

C     Eps2 and Eps4
               PVIS    = 1.0/VIS(L)
               CMUEFF  = -0.5*(BETA1 + WII*BETA6)
               EPS2(L) = 1.0 + CMUEFF*RK(L)*TTS(L)*PVIS
               EPS4    = 1.0 + CMUEFF*RK(L)*TTS(L)*PVIS*PRS 
               
C     Lo-limiter (minimum)
               EPS2LO  = 1.0 + FRSMUT*PVIS
               EPS4LO  = 1.0 + FRSMUT*PVIS*PRS
               EPS2(L) = MAX(EPS2LO, EPS2(L))
               EPS4    = MAX(EPS4LO, EPS4)

C     Hi-limiter (maximum)
               EPS2(L) = MIN(EPS2UP, EPS2(L))
               EPS4    = MIN(EPS4UP, EPS4)

C    Turbulent viscosity
               VIST(L) = (EPS2(L) - 1.0)*VIS(L)
           endif ! esto tähän ja vain epäisotropia
            END DO              ! End of I-loop
         END DO                 ! End of J-loop
      END DO                    ! End of K-loop

      RETURN
      END


      SUBROUTINE ANISFLUX(IMAX,JMAX,KMAX,JADD,KADD,IL,U,V,W,RK,
     &                    BIJ,A,AX,AY,AZ,FRM,FRN,FRW,FE,MAXEB)

C     A subroutine to add the anisotropic parts of the Reynolds 
C     stresses to the momentum and energy fluxes   
C             13.11.2000  Ville H. / Laboratory of Aerodynamics

      USE NS3CO, ONLY : IN,JN,KN

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: IMAX,JMAX,KMAX,JADD,KADD,IL,MAXEB
      REAL , INTENT(IN) , DIMENSION(*) :: U,V,W,RK,A,AX,AY,AZ
      REAL , INTENT(IN) , DIMENSION(MAXEB,*) :: BIJ
      REAL , DIMENSION(*) :: FRM,FRN,FRW,FE

      INTEGER :: ISTRID,JSTRID,LSTRID,IJSLB,KG,IA,IG,I,IMIL
      REAL :: UCM,VCM,WCM,B11CM,B12CM,B13CM,B22CM,B23CM,B33CM
      REAL :: RKA,F2EX,F3EX,F4EX,F5EX

C     Constants for the loops over the grid
      ISTRID = IMAX + 2*IN      ! Total number of cells in i-direction
      JSTRID = JMAX + 2*JN      ! Total number of cells in j-direction
      LSTRID = ISTRID*JSTRID    ! Total number of cells in ij-plane
      IJSLB  = (JMAX+JADD)*ISTRID ! Number of cells in ij-slab
      
      DO KG = 1,KMAX+KADD
         IA = (KG+KN-1)*LSTRID + JN*ISTRID
         DO  IG = 1,IJSLB
            I     = IA + IG
            IMIL  = I - IL      ! IL = ISTR, JSTR or KSTR 

C     Mean values of velocity components
            UCM  = 0.5*(U(I) + U(IMIL))
            VCM  = 0.5*(V(I) + V(IMIL)) 
            WCM  = 0.5*(W(I) + W(IMIL))

C     Twice the mean values of extra anisotrophy components b_ij
C     Remember that a_ij = 2.0*b_ij
            B11CM = BIJ(I,1) + BIJ(IMIL,1)
            B12CM = BIJ(I,2) + BIJ(IMIL,2)           
            B13CM = BIJ(I,3) + BIJ(IMIL,3)
            B22CM = BIJ(I,4) + BIJ(IMIL,4)
            B23CM = BIJ(I,5) + BIJ(IMIL,5)
            B33CM =-BIJ(I,1) - BIJ(I,4) - BIJ(IMIL,1) - BIJ(IMIL,4)

C     The change of the momentum and energy fluxes
            RKA   = 0.5*(RK(I) + RK(IMIL))*A(I)
            F2EX  = RKA*(AX(I)*B11CM + AY(I)*B12CM + AZ(I)*B13CM)
            F3EX  = RKA*(AX(I)*B12CM + AY(I)*B22CM + AZ(I)*B23CM)
            F4EX  = RKA*(AX(I)*B13CM + AY(I)*B23CM + AZ(I)*B33CM)
            F5EX  = UCM*F2EX + VCM*F3EX + WCM*F4EX

C     Updated momentum and energy fluxes
            FRM(I)  = FRM(I) + F2EX
            FRN(I)  = FRN(I) + F3EX
            FRW(I)  = FRW(I) + F4EX
            FE(I)   = FE(I)  + F5EX
            
         END DO                 
      END DO                    

      RETURN
      END

      SUBROUTINE FLUXA(FRM,FRN,FRW,RO,RM,RN,RW,P,RK,BIJ,UROT,A,AX,AY,AZ,
     &     IMAX,JMAX,KMAX,IL,JADD,KADD,MAXEB,ITURB,ISTRES,m)
C     ****************

C     A subroutine to calculate the momentum fluxes for an axisymmetric 
C     flow-case. Increment in DIM=1 direction is IL.
C              Modified 26.2.2001 Ville H. / Laboratory of Aerodynamics

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      REAL , DIMENSION(*) :: FRM,FRN,FRW
      REAL , INTENT(IN) , DIMENSION(*) :: RO,RM,RN,RW,P,RK
      REAL , INTENT(IN) , DIMENSION(*) :: A,AX,AY,AZ,UROT
      INTEGER , INTENT(IN) :: IMAX,JMAX,KMAX,IL,JADD,KADD
      INTEGER , INTENT(IN) :: ITURB,MAXEB,ISTRES,m

      REAL , PARAMETER :: P23 = 2.0/3.0
      REAL :: PAM,PM,UM,VM,WM,RKM,VFLOM,PERRO,RMCM
      REAL :: B11CM,B12CM,B13CM,B22CM,B23CM,B33CM,RKA,F2EX,F3EX,F4EX
      REAL , INTENT(IN) , DIMENSION(MAXEB,*) :: BIJ
      INTEGER :: II,IM,ISTRID,JSTRID,IJSLB,KG,IG,IA,I1

      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN

      IJSLB  = ISTRID*(JMAX+JADD) - IN + 1
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN ! 2-eq models, combine P and RK
         DO KG = 1,KMAX+KADD
            IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
            DO IG = 1+IN, IJSLB
               II      = IA + IG
               IM      = II - IL
               
               UM      = 0.5*(RM(IM)+RM(II))
               VM      = 0.5*(RN(IM)+RN(II))
               WM      = 0.5*(RW(IM)+RW(II))
               RKM     = 0.5*(RK(IM)+RK(II))
               PM      = 0.5*(P(IM) + P(II))
               PERRO   = 1.0/RO(II)
               RMCM    = AX(II)*UM + AY(II)*VM + AZ(II)*WM
               VFLOM   = A(II)*(RMCM*PERRO-UROT(II))
               PAM     = (PM+P23*RKM)*A(II)
               FRM(II) = FRM(II) + PAM*AX(II) + VFLOM*UM
               FRN(II) = FRN(II) + PAM*AY(II) + VFLOM*VM
c               FRW(II) = PAM*AZ(II) + VFLOM*WM
               FRW(II) = FRW(II) + PAM*AZ(II) + VFLOM*WM
            END DO
         END DO

         IF (ISTRES > 0) THEN ! 2-eq models with EARSM
            DO KG = 1, KMAX+KADD
               IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
               DO IG = 1+IN, IJSLB
                  II      = IA + IG
                  IM      = II - IL
C     Twice the mean values of extra anisotrophy components b_ij
C     (a_ij = 2.0*b_ij)
                  B11CM   = BIJ(II,1) + BIJ(IM,1)
                  B12CM   = BIJ(II,2) + BIJ(IM,2)           
                  B13CM   = BIJ(II,3) + BIJ(IM,3)
                  B22CM   = BIJ(II,4) + BIJ(IM,4)
                  B23CM   = BIJ(II,5) + BIJ(IM,5)
                  B33CM   =-BIJ(II,1) -BIJ(II,4) -BIJ(IM,1) -BIJ(IM,4)
                  
C     The change of the momentum fluxes
                  RKA     = 0.5*(RK(II) + RK(IM))*A(II)
                  F2EX    = RKA*(AX(II)*B11CM+AY(II)*B12CM+AZ(II)*B13CM)
                  F3EX    = RKA*(AX(II)*B12CM+AY(II)*B22CM+AZ(II)*B23CM)
                  F4EX    = RKA*(AX(II)*B13CM+AY(II)*B23CM+AZ(II)*B33CM)
                  
c     Updated momentum fluxes
                  FRM(II) = FRM(II) + F2EX
                  FRN(II) = FRN(II) + F3EX
                  FRW(II) = FRW(II) + F4EX
               END DO
            END DO
         END IF

      ELSE                      ! 0-eq models, use P alone

         DO KG = 1, KMAX+KADD
            IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
            DO IG = 1+IN, IJSLB
               II      = IA + IG
               IM      = II - IL

               UM      = 0.5*(RM(IM)+RM(II))
               VM      = 0.5*(RN(IM)+RN(II))
               WM      = 0.5*(RW(IM)+RW(II))

               PERRO   = 1.0/RO(II)
               RMCM    = AX(II)*UM + AY(II)*VM + AZ(II)*WM
               VFLOM   = A(II)*(RMCM*PERRO-UROT(II))

               PAM     = P(II)*A(II)
               FRM(II) = FRM(II) + PAM*AX(II) + VFLOM*UM
               FRN(II) = FRN(II) + PAM*AY(II) + VFLOM*VM
               FRW(II) = FRW(II) + PAM*AZ(II) + VFLOM*WM
            END DO
         END DO
      END IF

      RETURN
      END

      SUBROUTINE FLUXAMF(FRM,FRN,FRW,VAR,PRO,P,RK,BIJ,UROT,A,AX,AY,
     &     AZ,IMAX,JMAX,KMAX,IL,JADD,KADD,MAXEB,ITURB,ISTRES,PRC)
C     ****************

C     A subroutine to calculate the momentum fluxes for an axisymmetric 
C     flow-case. Increment in DIM=1 direction is IL.
C              Modified 26.2.2001 Ville H. / Laboratory of Aerodynamics
C              Modified for the two-fluid model 8.3.2017 TSii

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      REAL , DIMENSION(*) :: FRM,FRN,FRW
      REAL , INTENT(IN) , DIMENSION(*) :: P,RK
      REAL , INTENT(IN) , DIMENSION(*) :: A,AX,AY,AZ,UROT
      INTEGER , INTENT(IN) :: IMAX,JMAX,KMAX,IL,JADD,KADD
      INTEGER , INTENT(IN) :: ITURB,MAXEB,ISTRES

      REAL , PARAMETER :: P23 = 2.0/3.0
      REAL :: PAM,PM,UM,VM,WM,RKM,VFLOM,RMCM,ALPO
      REAL :: B11CM,B12CM,B13CM,B22CM,B23CM,B33CM,RKA,F2EX,F3EX,F4EX
      REAL :: RMIM,RMII,RNIM,RNII,RWIM,RWII,RM,RN,RW
      REAL , INTENT(IN) , DIMENSION(MAXEB,*) :: BIJ
      INTEGER :: II,IM,ISTRID,JSTRID,IJSLB,KG,IG,IA,IPHASE
      INTEGER :: III,JJJ,KKK
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(PRE_COR) PRC(*)

      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN

      IJSLB  = ISTRID*(JMAX+JADD) - IN + 1
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN ! 2-eq models, combine P and RK
         DO IPHASE = 1,NPHASES
         DO KG = 1,KMAX+KADD
            IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
            DO IG = 1+IN, IJSLB
               II      = IA + IG
               IM      = II - IL
            CALL IJKPAI(II,IMAX,JMAX,KMAX,III,JJJ,KKK)
  
               RMIM    = VAR(IM)%U(IPHASE)*VAR(IM)%ALFA(IPHASE)*
     +                  PRO(IM)%RO(IPHASE)             
               RMII    = VAR(II)%U(IPHASE)*VAR(II)%ALFA(IPHASE)*
     +                  PRO(II)%RO(IPHASE)             
               RNIM    = VAR(IM)%V(IPHASE)*VAR(IM)%ALFA(IPHASE)*
     +                  PRO(IM)%RO(IPHASE)             
               RNII    = VAR(II)%V(IPHASE)*VAR(II)%ALFA(IPHASE)*
     +                  PRO(II)%RO(IPHASE)   
               RWIM    = VAR(IM)%W(IPHASE)*VAR(IM)%ALFA(IPHASE)*
     +                  PRO(IM)%RO(IPHASE)             
               RWII    = VAR(II)%W(IPHASE)*VAR(II)%ALFA(IPHASE)*
     +                  PRO(II)%RO(IPHASE)
         
               RM      = 0.5*(RMIM + RMII)
               RN      = 0.5*(RNIM + RNII)
               RW      = 0.5*(RWIM + RWII)

               UM      = 0.5*(VAR(IM)%U(IPHASE)+VAR(II)%U(IPHASE))
               VM      = 0.5*(VAR(IM)%V(IPHASE)+VAR(II)%V(IPHASE))
               WM      = 0.5*(VAR(IM)%W(IPHASE)+VAR(II)%W(IPHASE))

               RKM     = 0.5*(VAR(IM)%ALFA(IPHASE)*RK(IM) + 
     +                        VAR(II)%ALFA(IPHASE)*RK(II))
               PM      = 0.5*(VAR(IM)%ALFA(IPHASE)*P(IM) + 
     +                        VAR(II)%ALFA(IPHASE)*P(II))
      alpo = 0.5*(VAR(IM)%ALFA(1)*P(IM) + 
     +                        VAR(II)%ALFA(1)*P(II)) +
     +       0.5*(VAR(IM)%ALFA(2)*P(IM) + 
     +                        VAR(II)%ALFA(2)*P(II))
               RMCM    = AX(II)*UM + AY(II)*VM + AZ(II)*WM
               VFLOM   = A(II)*(RMCM - UROT(II))
c               PAM    = (PM+P23*RKM)*A(II) ! Antakee mersu
               PAM     = P23*RKM*A(II)
               VAR(II)%FRM(IPHASE) = VAR(II)%FRM(IPHASE) + 
     +         PAM*AX(II) + VFLOM*RM
               VAR(II)%FRN(IPHASE) = VAR(II)%FRN(IPHASE) + 
     +         PAM*AY(II) + VFLOM*RN
               VAR(II)%FRW(IPHASE) = VAR(II)%FRW(IPHASE) + 
     +         PAM*AZ(II) + VFLOM*RW
               VAR(II)%FA(IPHASE) = A(II)*VAR(II)%ALFA(IPHASE)
               FRM(II) = FRM(II) + VAR(II)%FRM(IPHASE)
               FRN(II) = FRN(II) + VAR(II)%FRN(IPHASE)
               FRW(II) = FRW(II) + VAR(II)%FRW(IPHASE)
               PRC(II)%FPM = A(II)*AX(II)*(P(II)+P(IM))*.5
               PRC(II)%FPN = A(II)*AY(II)*(P(II)+P(IM))*.5
               PRC(II)%FPW = A(II)*AZ(II)*(P(II)+P(IM))*.5
            END DO
         END DO
         END DO ! NPHASES

         IF (ISTRES > 0) THEN ! 2-eq models with EARSM
            DO KG = 1, KMAX+KADD
               IA = ((KG-1+KN)*JSTRID + JN)*ISTRID
               DO IG = 1+IN, IJSLB
                  II      = IA + IG
                  IM      = II - IL
C     Twice the mean values of extra anisotrophy components b_ij
C     (a_ij = 2.0*b_ij)
                  B11CM   = BIJ(II,1) + BIJ(IM,1)
                  B12CM   = BIJ(II,2) + BIJ(IM,2)           
                  B13CM   = BIJ(II,3) + BIJ(IM,3)
                  B22CM   = BIJ(II,4) + BIJ(IM,4)
                  B23CM   = BIJ(II,5) + BIJ(IM,5)
                  B33CM   =-BIJ(II,1) -BIJ(II,4) -BIJ(IM,1) -BIJ(IM,4)
                  
C     The change of the momentum fluxes
                  RKA     = 0.5*(RK(II) + RK(IM))*A(II)
                  F2EX    = RKA*(AX(II)*B11CM+AY(II)*B12CM+AZ(II)*B13CM)
                  F3EX    = RKA*(AX(II)*B12CM+AY(II)*B22CM+AZ(II)*B23CM)
                  F4EX    = RKA*(AX(II)*B13CM+AY(II)*B23CM+AZ(II)*B33CM)
                  
c     Updated momentum fluxes (how to divide these between the phases)?
               STOP 'FLUXAMF: EARSM does not work yet with multiphase' 
               END DO
            END DO
         END IF

      ELSE ! 0-eq models do not work in multiphase flow

         STOP 'FLUXAMF: Use two-equation models with multiphase'

      END IF

      RETURN
      END SUBROUTINE FLUXAMF


      SUBROUTINE PEEKOO(VIS,EPS2,SIJ,BIJ,PHITUR,IMAX,JMAX,KMAX,
     &                  RK,PTUR,MAXEB,MAX11)

C     A subroutine to calculate the production of turbulence with
C     the EARSM 
C                 2.11.2000 Ville H. / Laboratory of Aerodynamics

      USE CONSTANTS, ONLY : EPS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE
      
      INTEGER , INTENT(IN) :: IMAX,JMAX,KMAX,MAXEB,MAX11
      INTEGER :: ISTRID,JSTRID,IL,IXCL,K,IA,IJ,L
*     REAL , INTENT(IN) , DIMENSION(MAXEB,9) :: DUIDXJ
      REAL , INTENT(IN) , DIMENSION(MAX11,6) :: SIJ
      REAL , INTENT(IN) , DIMENSION(MAXEB,*) :: BIJ
      REAL , INTENT(IN) , DIMENSION(*) :: EPS2,VIS,RK,PHITUR
      REAL , DIMENSION(*) :: PTUR
      REAL :: P3,S11,S22,S33,S12,S13,S23
      REAL :: RMUT,DIVVP3,HUU,HUV,HUW,HVV,HVW,HWW

C     One third
      P3     = 1.0/3.0      

C     Constants for the loops over the grid
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID
      IXCL   = JN*ISTRID + IN    

      DO K = 1,KMAX
         IA  = (KN+K-1)*IL
         DO IJ = 1+IXCL,IL-IXCL
            L  = IA + IJ

C     Strain rate tensor components
            S11     = SIJ(L,1)
            S22     = SIJ(L,4)
            S33     = SIJ(L,6)
            
            S12     = SIJ(L,2)
            S13     = SIJ(L,3)
            S23     = SIJ(L,5)
            
C     Effective eddy-viscosity coefficient
            RMUT    = MAX(EPS2(L) - 1.0 , EPS)*VIS(L)
            
C     One third of the divergence of the velocity
            DIVVP3  = P3*(S11 + S22 + S33)

C     Half of the Reynolds stresses (with rho) 0.5*rho*u_i*u_j
            HUU     = MAX((RK(L)*(P3 + BIJ(L,1))
     &           - RMUT*(S11 - DIVVP3)) , 0.0)
            HVV     = MAX((RK(L)*(P3 + BIJ(L,4))
     &           - RMUT*(S22 - DIVVP3)) , 0.0)
            HWW     = MAX((RK(L)*(P3 - BIJ(L,1) - BIJ(L,4))
     &           - RMUT*(S33 - DIVVP3)) , 0.0)
            HUV     = RK(L)*BIJ(L,2)  - RMUT*S12
            HUW     = RK(L)*BIJ(L,3)  - RMUT*S13 
            HVW     = RK(L)*BIJ(L,5)  - RMUT*S23

C     Production of turbulence
            PTUR(L) =-2.0*(HUU*S11 + HVV*S22 + HWW*S33
     &              + 2.0*(HUV*S12 + HUW*S13 + HVW*S23))
            
C     Forced laminarity by limited production if PHITUR (VTRAN) is zero
            PTUR(L) = PTUR(L)*PHITUR(L)
            
         END DO  
      END DO     

      RETURN
      END


      SUBROUTINE TTIMES(RO,VIS,RK,REPS,ITURB,IMAX,JMAX,KMAX,TTS)

C     A subroutine to calculate the turbulent time scale for
C     the EARSM
C            2.11.2000 Ville H. / Laboratory of Aerodynamics
 
      USE CONSTANTS, ONLY : EPS
      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: ITURB,IMAX,JMAX,KMAX
      INTEGER :: ISTRID,JSTRID,IL,KA,K,KK,J,JJ,L,I
      REAL , INTENT(IN) , DIMENSION(*) :: VIS,RO,RK,REPS
      REAL , DIMENSION(*) :: TTS
      REAL :: CT,BSTAR,PREPS,TSTUR,TSVIS,BNUM

C     Model coefficients
      CT     = 6.0
      BSTAR  = 0.09

C     Constants for the loops over the grid
      ISTRID = IMAX + 2*IN      ! Total number of cells in i-direction
      JSTRID = JMAX + 2*JN      ! Total number of cells in j-direction
      IL     = ISTRID*JSTRID    ! Total number of cells in ij-plane
      KA     = (KN-1)*IL + (JN-1)*ISTRID + IN ! The start address

      IF (ITURB == 3) THEN      ! k-epsilon models
         DO K = 0, KMAX+1
            KK = K*IL + KA            
            DO J = 0, JMAX+1
               JJ = J*ISTRID + KK
               DO I = 0, IMAX+1
                  L = JJ + I 
                  PREPS  = 1.0/ABS(REPS(L)+EPS)
                  TSTUR  = RK(L)*PREPS
                  TSVIS  = CT*SQRT(VIS(L)*PREPS)
                  TTS(L) = MAX(TSTUR,TSVIS)
               END DO  
            END DO     
         END DO        
         
      ELSE IF (ITURB == 6) THEN ! k-omega models
         DO K = 0, KMAX+1
            KK = K*IL + KA            
            DO J = 0, JMAX+1
               JJ = J*ISTRID + KK
               DO I = 0, IMAX+1
                  L = JJ + I 
                  BNUM   = 1.0/ABS(BSTAR*REPS(L)+EPS)
                  TSTUR  = RO(L)*BNUM
                  TSVIS  = CT*SQRT(RO(L)*VIS(L)*BNUM/ABS(RK(L)+EPS))
                  TTS(L) = MAX(TSTUR,TSVIS) 
               END DO
            END DO
         END DO 
         
C     ELSE IF (ITURB == X) THEN ! k-tau models
C     k-tau models are not in use, yet

      END IF

      RETURN
      END SUBROUTINE TTIMES


