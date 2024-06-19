C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLSCLD(RO,FFI,FI,HAT1,HAT2,HAT3,HAT4,VOL,A,
     2 PSIGSC,PSIGS2,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 EPS2,A2XA,A2YA,A2ZA,PRL,PRT,
     4 NSCAL,MAXSB,NBL,RLOLIM,FRSDEN,ITURB,
     5 KSTR,IDIR)

C ... SCALAR FLUX

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: RO(*),FFI(MAXSB,MAX(1,NSCAL)),FI(MAXSB,MAX(1,NSCAL)),
     2  HAT1(*),HAT2(*),HAT3(*),HAT4(*),VOL(*),A(*),VIS(*),
     3  EPS2(*),A2XA(*),A2YA(*),A2ZA(*),PSIGSC(*),PSIGS2(*),
     4  RLOLIM(*)

      LOGICAL :: INVIS

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      ILL     = ISTRID*JSTRID

      IL      = KSTR
      CALL SINVIS(IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID,
     +     IMAX,JMAX,KMAX,IDI,INVIS,IDIR)

c      write(*,'(8I4)') nbl,IMINID,IMAXID,JMINID,JMAXID,KMINID,KMAXID
      if (iturb >= 25) write(*,*) 'scalaareja ja rsm. en usko'
      IF(ITURB /= 24 ) THEN
      DO 313 NS = 1,NSCAL
      PSIGE1  = 1./PSIGSC(NS)
      PSIGE2  = 1./PSIGS2(NS)
C ... SCHMIDT'S NUMBER  = PSIGE = 1/PSIGE1
C
C ... VISCOUS FLUXES ARE ADDED
C
      DO 900 KG = 1,KMAXID
         IA      = (KG+KN-1)*ILL
         DO 850 JG = 1,JMAXID
         IA2       = IA + (JG+JN-1)*ISTRID + IN
         DO 850 IG = 1,IMAXID
            I      = IA2 + IG
            YPD2    = A(I)/(VOL(I) + VOL(I-IL))
            SURG    = -A(I)*YPD2
            SURFI1 = SURG*PSIGE1*
     +           ((EPS2(I)-1.)*VIS(I) + (EPS2(I-IL)-1.)*VIS(I-IL))
            SURFI2 = SURG*PSIGE2*(VIS(I) + VIS(I-IL))
            SURFFI = SURFI1 + SURFI2
            YPRO1   = 1./RO(I)
            YPRO2   = 1./RO(I-IL)
            F1VIG   = SURFFI*(FI(I,NS)*YPRO1- FI(I-IL,NS)*YPRO2)
            FFI(I,NS) = F1VIG
 850     CONTINUE
 900  CONTINUE
      
      IF(INVIS) THEN ! AN EXPLICITLY GIVEN INVISCID REGION IS INCLUDED
      DO 1000 KG = KMINID,KMAX+1
      IA      = (KG+KN-1)*ILL
      DO 950 JG = JMINID,JMAX+1
      IA2       = IA + (JG+JN-1)*ISTRID + IN
      DO 950 IG = IMINID,IMAX+1
         I       = IA2 + IG
         FFI(I,NS)  = 0.
 950  CONTINUE
 1000 CONTINUE
      ENDIF
 313  CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE FLSCLD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FLSCAL(RO,FFI,FI,HAT1,HAT2,HAT3,HAT4,VOL,A,
     2 PSIGSC,PSIGS2,IMAX,JMAX,KMAX,ICYCLE,INTEM,IDI,VIS,
     3 EPS2,A2XA,A2YA,A2ZA,PRL,PRT,
     4     NSCAL,MAXSB,NBL,RLOLIM,FRSDEN,ITURB,KMAXP1,KSTR,IDIR,
     5     ZZZ,MAW,MAXW,RKSI,INCHIML)

C ... SCALAR FLUX

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION :: RO(*),FFI(MAXSB,MAX(1,NSCAL)),FI(MAXSB,MAX(1,NSCAL)),
     2  HAT1(*),HAT2(*),HAT3(*),HAT4(*),VOL(*),A(*),VIS(*),
     3  EPS2(*),A2XA(*),A2YA(*),A2ZA(*),PSIGSC(*),PSIGS2(*),
     4  RLOLIM(*),RKSI(*)

      REAL, POINTER :: ROIP(:),ROIM(:),REIM(:),REIP(:),ROHAT(:)
      REAL, TARGET ::  ZZZ(MAXW)
      LOGICAL :: INCHIML

      ROIP => ZZZ( 0*MAW+1: 1*MAW);ROIM => ZZZ( 1*MAW+1: 2*MAW)
      REIM => ZZZ( 2*MAW+1: 3*MAW);REIP => ZZZ( 3*MAW+1: 4*MAW)
      ROHAT=> ZZZ( 4*MAW+1: 5*MAW)

C ... SECOND-ORDER UPWIND IS FORCED FOR SCALAR QUATITIES experiment
c              RK3   = 0.
c              RK4   = 2.
C
C ... FLUXES FOR THREE-DIMENSIONAL FLOW
C

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
*      IJTRID  = ISTRID*(JMAX + JN)
      IJTRID  = ISTRID*(JMAX + 2)
      IL      = KSTR
      KL      = 2*IL

      DO 1000 KG = 1,KMAXP1
*      IA      = (KG+KN-1)*ISTRID*JSTRID + ISTRID
      IA      = (KG+KN-1)*ISTRID*JSTRID + (JN-1)*ISTRID
C
C ... EXTRAPOLATION OF THE DEPENDENT VARIABLES

      CALL XXTRAP(RO,ROIP,ROIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)

      DO IG= 1,IJTRID
         I       = IA + IG
         ROIP(IG)    = MAX(0.001*FRSDEN,ROIP(IG))
         ROIM(IG)    = MAX(0.001*FRSDEN,ROIM(IG))
         ROHAT(IG)   = SQRT(ROIP(IG)*ROIM(IG))
      ENDDO

      DO 313 NS = 1,NSCAL

      CALL XXTRAP(FI(1,NS),REIP,REIM,IJTRID,IL,KL,IA,INTEM,RKSI,INCHIML)


C*******************************************************************
C ... ROE'S FLUX SPLITTING                                         *
C*******************************************************************
      DO 600 IG=1,IJTRID
      I       = IA + IG
      YPROR   = 1./ROIP(IG)
      FIR     = MAX(RLOLIM(NS),REIP(IG))*YPROR

      YPROL   = 1./ROIM(IG)
      FIL     = MAX(RLOLIM(NS),REIM(IG))*YPROL
      SQRL    = 1./(ROIM(IG) + ROHAT(IG))
      CQL     = ROIM(IG)*SQRL
      CQR     = 1.- CQL

      FIHAT  = CQR*FIR+ CQL*FIL
      ALFA7   = ROHAT(IG)*(FIR - FIL)
C ... ABSOLUTE VALUES OF THE CHARACTERISTIC SPEEDS

C*    HATS ARE FROM FLUXKE SUBROUTINE PPR 11.2.
      DAMP  = HAT1(I)
      ROUCR = HAT2(I)
      ROUCL = HAT3(I)
      RL7   = HAT4(I)
C ... CONVECTIVE FLUXES
      F1PIG   = .5*A(I)*(ROUCR*FIR+ROUCL*FIL - DAMP*FIHAT- RL7*ALFA7)

      FFI(I,NS) = FFI(I,NS) + F1PIG

        
 600  CONTINUE
 313  CONTINUE

 1000 CONTINUE

      RETURN
      END SUBROUTINE FLSCAL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOSCAL(F3FI,VOL,A3,A3X,A3Y,
     2  A3Z,VIS,EPS2,PSIGSC,PSIGS2,ISTRID,JSTRID,KSTRID,IDIM,KBOT,
     3  KX1,KX2,KY1,KY2,PR,PRT,RO,
     4  M,IN,JN,KN,IDIR,FI,NSCAL,MAXSB,ITURB,
     5  ISTR,JSTR,KSTR)

      REAL :: RO(*),A3(*),PSIGSC(*),PSIGS2(*),
     2 A3Y(*),A3Z(*),A3X(*),VIS(*),EPS2(*),
     3 VOL(*),FI(MAXSB,MAX(1,NSCAL)),F3FI(MAXSB,MAX(1,NSCAL))
C
C ... LU FACTORED METHOD.BOUNDARY CONDITION FOR THE SOLID BOUNDARY
C ... BOUNDARY FLUXES ARE CALCULATED EXPLICITLY. POINTS ARE NUMBERED
C ... XI-WISE
C

      EPS = 1.E-10
      RC1 = 1.500
      RC2 = RC1 - 1.

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID

      YMR = .1666667

      DO 313 NS = 1,NSCAL
C ... IF REYNOLDS STRESSES ARE APPLIED FISURF WILL BE ZERO PPR 18.5.1994
      IF (ITURB >= 20 .AND. NS <= 6) THEN
         RDIM = 0.
      ELSE
         RDIM = 1.
      ENDIF

      KA = (KN+KBOT-1)*KSTR

C ... VISCOUS FLUXES

      IF(IDIM > KBOT) THEN

      DO 1000 J = KY1,KY2
      IJ      = (JN+J-1)*JSTR + KA
      NN       = (J-1)*(KX2-KX1+1) + 1 - KX1
      DO 1000 I = KX1,KX2
         II      = 1 + (IN+I-1)*ISTR  + IJ
         IF1     = II + NFL
         DIST    = (A3(II)+A3(IF1+LSTRID))/VOL(II)
         SURG    =-DIST*YMR
         EPS2B   = MAX(1.,    RC1*EPS2(II) - RC2*EPS2(II+LSTRID))
         VISB    = MAX(EPS,   RC1*VIS(II)  - RC2*VIS(II+LSTRID))
         FISURF  = RDIM*MAX(0.,RC1*FI(II,NS) - RC2*FI(II+LSTRID,NS))
         SURF    = SURG*EPS2B*VISB/PSIGS2(NS)
         YPRO1   = 1./RO(II)
         YPRO2   = 1./RO(II+LSTRID)
         FIVII    = IDIR*(9.000*FI(II,NS)*YPRO1 - 1.000*
     +        FI(II+LSTRID,NS) * YPRO2 - 8.*FISURF)
         F1V     = SURF*FIVII
         F3FI(IF1,NS)    = A3(IF1)*F1V
1000  CONTINUE

      ELSE

      DO J = KY1,KY2

         IJ = (JN+J-1)*JSTR + KA
         NN = (J-1)*(KX2-KX1+1) + 1 - KX1

         DO I = KX1,KX2

            II  = 1  + (IN+I-1)*ISTR  + IJ
            IF1 = II + NFL
            F3FI(IF1,NS) = 0.

         ENDDO

      ENDDO

      ENDIF

C ... CONVECTIVE FLUXES:  EXTRAPOLATION OF THE PRESSURE

 313  CONTINUE

      RETURN
      END SUBROUTINE BOSCAL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SWEEPS(DFI,F3FI,VOL,IMAX,JMAX,KMAX,KSTRID,IN,JN,KN,DT,
     2                  NSCAL,MAXSB,IDO)

      DIMENSION :: VOL(*),DFI(MAXSB,MAX(1,NSCAL)),
     +             F3FI(MAXSB,MAX(1,NSCAL))
C
C ... EXPLICIT STEP OF THE LU-FACTORED METHOD

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
      IF(IDO == 0) THEN         ! first two  sweeps
      DO NS= 1,NSCAL
      DO K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO I = 1,IMAX
         N       = JJ + I
         DFI(N,NS) = DFI(N,NS)-(F3FI(N+KSTRID,NS)- F3FI(N,NS))
c       CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c      write(680,*) ii,jj,kk,F3FI(N+KSTRID,NS),F3FI(N,NS)
      ENDDO ; ENDDO ; ENDDO ; ENDDO

      ELSEIF(IDO == 1) THEN     ! last sweep

      DO NS= 1,NSCAL
      DO K = 1,KMAX
      KA      = (KN+K-1)*IL
      DO J = 1,JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO I = 1,IMAX
         N       = JJ + I
         DT1     = 1./VOL(N)
         DFI(N,NS) = (DFI(N,NS)-(F3FI(N+KSTRID,NS)- F3FI(N,NS)))*DT1
      ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SWEEPS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FINSC(FI,DFI,FIOLD,IMAX,JMAX,KMAX,IJMASK,
     +                 MAXSB,NSCAL)

      USE NS3CO, ONLY : IN, JN, KN

      DIMENSION FI(MAXSB,MAX(1,NSCAL)),DFI(MAXSB,MAX(1,NSCAL)),
     +          FIOLD(MAXSB,MAX(1,NSCAL)),
     +          IJMASK(*)
C ... ADD DELTA VARIABLES TO FI:S PPR 15.2.

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
C
C ... VALUES AFTER THE INTEGRATION STEP
C
      DO 1000 NS = 1,NSCAL
      DO 1000 K = 1,KMAX
      IA      = (KN+K-1)*IL + JN*ISTRID
      DO 1000 I = 1,JMAX*ISTRID
      N       = IA + I
      IF(IJMASK(I) == 1) THEN
           FI(N,NS) = FIOLD(N,NS)  + DFI(N,NS)
      ENDIF
1000  CONTINUE

      RETURN
      END SUBROUTINE FINSC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOLISC(FI,IMAX,JMAX,KMAX,IN,JN,KN,RLIM,ULIM,MAXSB,
     +                  NSCAL)

      REAL :: FI(MAXSB,MAX(1,NSCAL)), RLIM(*), ULIM(*)
C
C ... LIMIT THE VALUES OF FI GREATER THAN RLIM SC EQ. PPR 4.3
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
      DO 1000 NS = 1,NSCAL
      DO 1000 K  = 1,KMAX
      IA       = (KN+K-1)*IL + JN*ISTRID
      DO 1000 I = 1,JMAX*ISTRID
      N        = IA + I
      FI(N,NS) = MAX(FI(N,NS),RLIM(NS))
      FI(N,NS) = MIN(FI(N,NS),ULIM(NS))
1000  CONTINUE
C ... WARNING AD HOC MODIFICATIONS FOR TRANSITION
c ... possible transition here
c      DO 2000 K = 1,KMAX
c         IA     = (KN+K-1)*IL + JN*ISTRID
c         DO 2000 J = 1,JMAX
c            IJ     = IA + (J-1)*ISTRID  
c            DO 2000 I = 85,164
c               N       = I + IJ
c               FI(N,1) = RLIM(1)
c               FI(N,2) = 0.
c               FI(N,3) = 0.
c               FI(N,4) = RLIM(4)
c               FI(N,5) = 0.
c               FI(N,6) = RLIM(6)
c 2000 CONTINUE

      RETURN
      END SUBROUTINE LOLISC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ADDSSC(DFI,SFI,VOL,IMAX,JMAX,KMAX,KSTRID,
     +                  IN,JN,KN,DT,MAXSB,NSCAL)

      REAL :: VOL(*), DFI(MAXSB,MAX(1,NSCAL)), SFI(MAXSB,MAX(1,NSCAL))
C
C ... ADD SOURCE TERM TO SCALAR EQUATIONS
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
      DO 1000 NS = 1,NSCAL
      DO 1000 K = 1,KMAX
      IA      = (KN+K-1)*IL + JN*ISTRID
      DO 1000 I = 1,JMAX*ISTRID
      N       = IA + I
c      DT1     = DT
      DFI(N,NS)  = DFI(N,NS)  + SFI(N,NS)   !*DT1
1000  CONTINUE

      RETURN
      END SUBROUTINE ADDSSC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TIMDSC(DFI,FI,FILE2,FILE3,VOL,IMAX,JMAX,KMAX,KSTRID,
     +                  IN,JN,KN,DT,MAXSB,NSCAL)

      REAL :: VOL(*), DFI(MAXSB,MAX(1,NSCAL)), FI(MAXSB,MAX(1,NSCAL)),
     +        FILE2(MAXSB,MAX(1,NSCAL)), FILE3(MAXSB,MAX(1,NSCAL))
C
C ... ADD 2-ORDER TIME DERIVATIVE TERM TO SCALAR EQUATIONS (VOL=CONST.)
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      DO 1000 NS = 1,NSCAL
      DO 1000 K = 1,KMAX
      IA      = (KN+K-1)*IL + JN*ISTRID
      DO 1000 I = 1,JMAX*ISTRID
      N       = IA + I
      PDT     = 1./(VOL(N)*DT)  ! TRAP
      DFI(N,NS)  = DFI(N,NS) + PDT *   (-1.5*VOL(N)*FI(N,NS)
     +           + 2.*VOL(N)*FILE2(N,NS) -.5*VOL(N)*FILE3(N,NS))
1000  CONTINUE

      RETURN
      END SUBROUTINE TIMDSC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C


