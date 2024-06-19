C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTDIA(DM,DN,DW,DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT)
C
C ... SOURCE TERM ASSOCIATED WITH ROTATION IS FACTORED OUT
C

      IMPLICIT NONE

      REAL :: R1, R2, R3, R4, APU1, APU, OM1, OM2, OM3
      REAL :: OMEGA, OMEX, OMEY, OMEZ
      REAL :: DM(*), DN(*), DW(*), DTL(*)

      INTEGER :: I, NTOT
      
      IF(OMEX >= .9999) THEN    ! old x-axis rotation

         DO I = 1,NTOT

            R3    = DN(I)
            R4    = DW(I)
            APU1  = DTL(I)*OMEGA
            APU   = 1./(1.+ APU1**2)
            DN(I) = (R3 + APU1*R4)*APU
            DW(I) = (R4 - APU1*R3)*APU

         ENDDO

      ELSE                      ! new varying rotation axis

         DO I = 1,NTOT          ! onkohan vaan oikein??? 15.9.1999 PR

            R1    = DM(I)
            R2    = DN(I)
            R3    = DW(I)
            APU1  = DTL(I)*OMEGA
            APU   = 1./(1.+ APU1**2)
            OM1   = DTL(I)*OMEGA*OMEX
            OM2   = DTL(I)*OMEGA*OMEY
            OM3   = DTL(I)*OMEGA*OMEZ
            DM(I) = (1.+OM1**2)*R1+(OM1*OM2+OM3)*R2+(OM1*OM3-OM2)*R3 
            DN(I) = (1.+OM2**2)*R2+(OM1*OM2-OM3)*R1+(OM2*OM3+OM1)*R3 
            DW(I) = (1.+OM3**2)*R3+(OM1*OM3+OM2)*R1+(OM2*OM3-OM1)*R2 
            DM(I) = DM(I)*APU
            DN(I) = DN(I)*APU
            DW(I) = DW(I)*APU

         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE ROTDIA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTDIAMF(VAR,DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT,NPHASE)
C
C ... SOURCE TERM ASSOCIATED WITH ROTATION IS FACTORED OUT
C

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL :: R1, R2, R3, R4, APU1, APU, OM1, OM2, OM3, DMI, DNI, DWI
      REAL :: OMEGA, OMEX, OMEY, OMEZ
      REAL :: DTL(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      INTEGER :: I, NTOT, IPHASE, NPHASE
      
      IF(OMEX >= .9999) THEN    ! old x-axis rotation

      DO IPHASE = 1,NPHASE

         DO I = 1,NTOT

            R3    = VAR(I)%DN(IPHASE)
            R4    = VAR(I)%DW(IPHASE)
            APU1  = DTL(I)*OMEGA
            APU   = 1./(1.+ APU1**2)
            VAR(I)%DN(IPHASE) = (R3 + APU1*R4)*APU
            VAR(I)%DW(IPHASE) = (R4 - APU1*R3)*APU
         ENDDO

      ENDDO

      ELSE                      ! new varying rotation axis

      DO IPHASE = 1,NPHASE

         DO I = 1,NTOT          ! onkohan vaan oikein??? 15.9.1999 PR

            R1    = VAR(I)%DM(IPHASE)
            R2    = VAR(I)%DN(IPHASE)
            R3    = VAR(I)%DW(IPHASE)
            APU1  = DTL(I)*OMEGA
            APU   = 1./(1.+ APU1**2)
            OM1   = DTL(I)*OMEGA*OMEX
            OM2   = DTL(I)*OMEGA*OMEY
            OM3   = DTL(I)*OMEGA*OMEZ
            DMI   = (1.+OM1**2)*R1+(OM1*OM2+OM3)*R2+(OM1*OM3-OM2)*R3 
            DNI   = (1.+OM2**2)*R2+(OM1*OM2-OM3)*R1+(OM2*OM3+OM1)*R3 
            DWI   = (1.+OM3**2)*R3+(OM1*OM3+OM2)*R1+(OM2*OM3-OM1)*R2 
            VAR(I)%DM(IPHASE) = DMI*APU
            VAR(I)%DN(IPHASE) = DNI*APU
            VAR(I)%DW(IPHASE) = DWI*APU

         ENDDO

      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE ROTDIAMF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTRSM(DUU,DUV,DUW,DVV,DVW,DWW,DTL,OMEGA,
     +     OMEX,OMEY,OMEZ,NTOT)

      REAL :: DUU(*), DUV(*), DUW(*), DVV(*), DVW(*), DWW(*), DTL(*)
C
C ... SOURCE TERM ASSOCIATED WITH ROTATION IS FACTORED OUT in RSM
C
      IF(OMEX <=  .999) THEN
         WRITE(*,*) 'For a steady state Reynolds stresses the rotation'
         WRITE(*,*) 'axis must be x-axis, at least currently'
         STOP
      ENDIF

      DO 1000 I = 1,NTOT
      R1     = DUU(I)
      R2     = DUV(I)
      R3     = DUW(I)
      R4     = DVV(I)
      R5     = DVW(I)
      R6     = DWW(I)
      
      APU1   = DTL(I)*OMEGA
      ALFA   = 1./(1.+APU1**2)
      BETA   = 1./(1+4.*APU1**2)

      DUV(I) = (R2 + APU1*R3)*ALFA
      DUW(I) = (R3 - APU1*R2)*ALFA
      DVV(I) = (R4*(1.+2*APU1**2) + 2.*APU1*R5 + 2.*APU1**2*   R6)*BETA
      DVW(I) = (R4*(-APU1)        +         R5 +    APU1*      R6)*BETA
      DWW(I) = (R4*2.*APU1**2     - 2.*APU1*R5 +(1.+2*APU1**2)*R6)*BETA

1000  CONTINUE

      RETURN
      END SUBROUTINE ROTRSM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTFUN(UROT,VROT,WROT,NTOT,OMEGA)
C
C ... SWIRL VELOCITIES FOR ROTATIONAL FLOWS
C

      IMPLICIT NONE

      REAL :: OMEGA
      REAL :: UROT(*), VROT(*), WROT(*)

      INTEGER :: NN, NTOT
      
      DO NN = 1,NTOT

         UROT(NN) = UROT(NN)*OMEGA
         VROT(NN) = VROT(NN)*OMEGA
         WROT(NN) = WROT(NN)*OMEGA

      ENDDO

      RETURN
      END SUBROUTINE ROTFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE VELFUN(UROT,AX,AY,AZ,GVEX,GVEY,GVEZ,NTOT)

      REAL :: UROT(*), AX(*), AY(*), AZ(*)
C
C ... Give velocity for grid
C

      DO 2000 NN = 1,NTOT
         UROT(NN)= UROT(NN) + GVEX*AX(NN) + GVEY*AY(NN)+ GVEZ*AZ(NN)
 2000 CONTINUE

      RETURN
      END SUBROUTINE VELFUN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CHECK(UROT,VROT,WROT,IMAX,JMAX,
     1 KMAX,A1XA,A1YA,A1ZA,A2XA,A2YA
     2 ,A2ZA,A3XA,A3YA,A3ZA,A1,A2,A3,DRO,DE,DM,DN)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: UROT(*),VROT(*),WROT(*),A1(*),A1XA(*),A1YA(*),A1ZA(*),
     1  A2(*),A2XA(*),A2YA(*),A2ZA(*),A3(*),A3XA(*),A3YA(*),A3ZA(*),
     2  DRO(*),DE(*),DM(*),DN(*)

C ... THE GRID PROPERTIES ARE CHECKED

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = ISTRID*JSTRID
      IMAXP1  = IMAX + 1
      JMAXP1  = JMAX + 1
      KMAXP1  = KMAX + 1

      DO 1000 K = 1,KMAX
      KA  = (KN+K-1)*ISTRID*JSTRID
      DO 1000 J = 1,JMAX
      JJ  = KA + (JN+J-1)*ISTRID + IN
      DO 1000 I = 1,IMAX
      II  = JJ + I
      DRO(II) = A1(II+1)*A1XA(II+1) -  A1(II)*A1XA(II)
     +  + A2(II+ISTRID)*A2XA(II+ISTRID)-A2(II)*A2XA(II)
     +  + A3(II+KSTRID)*A3XA(II+KSTRID)-A3(II)*A3XA(II)
      DE(II) = A1(II+1)*A1YA(II+1) -   A1(II)*A1YA(II)
     +  + A2(II+ISTRID)*A2YA(II+ISTRID)-A2(II)*A2YA(II)
     +  + A3(II+KSTRID)*A3YA(II+KSTRID)-A3(II)*A3YA(II)
      DM(II) = A1(II+1)*A1ZA(II+1) -   A1(II)*A1ZA(II)
     +  + A2(II+ISTRID)*A2ZA(II+ISTRID)-A2(II)*A2ZA(II)
     +  + A3(II+KSTRID)*A3ZA(II+KSTRID)-A3(II)*A3ZA(II)
     +  + DE(II) + DRO(II)
      DN(II) = A1(II+1)*UROT(II+1) -   A1(II)*UROT(II)
     +  + A2(II+ISTRID)*VROT(II+ISTRID)-A2(II)*VROT(II)
     +  + A3(II+KSTRID)*WROT(II+KSTRID)-A3(II)*WROT(II)
 1000 CONTINUE

      RETURN
      END SUBROUTINE CHECK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Some notes on HCHECK
C     averaging :
C     quantity                weighting
C-----------------------------------------
C     pressure          PEE   massflow
C     total pressure    P0A   massflow
C     total temperature T0A   massflow
C     total enthalpy    QEH   massflow
C     density           ROAVE area          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HCHECK(F1R,F1RM,F1RN,F1RW,F1E,F1RK,F1EPS,P,RO,RM,RN,RW,
     +     RK,REPS,A1,A1X,A1Y,A1Z,UROT,IMAX,JMAX,KMAX,IN,JN,KN,
     +     ICYCLE,DIR,ISTR,JSTR,KSTR,ITURB,C,TEMP,GAMMA,XCO,YCO,ZCO,
     +     XC,YC,ZC,CH,VIS,PRO,VAR,PR,IPRESC,FRSPRE,FRSDEN,
     +     ISTATE,ISSB,SOLUTION_TYPE,VOL,PDIFF)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND, G0, TWO_FLUIDL
      
      IMPLICIT NONE

      REAL :: F1R(*),F1E(*),A1(*),P(*),F1RM(*),F1RN(*),F1RW(*),
     + RO(*),RM(*),RN(*),RW(*),RK(*),REPS(*),A1X(*),
     + A1Y(*),A1Z(*),UROT(*),TEMP(*),C(*),F1RK(*),F1EPS(*),
     + CH(*),VIS(*),VOL(*),PDIFF(*)!,XC(*),YC(*),ZC(*)
      REAL :: XCO(*),YCO(*),ZCO(*),XC(*),YC(*),ZC(*)
      REAL :: ROM,RKM,GAMMA
      REAL :: ROAVE, VEL,HSTAT1,HSTAT,HEAD2,XKAVE,TUAVE,XAAVE
      REAL :: RMAKK,RMAKM,VELKK,VELKM,VELAVE,RKFLX,EPSFLX
      REAL :: XAVE,YAVE,ZAVE,XOLD,YOLD,ZOLD,DISTA,DDISTA,PR,FRSPRE
      REAL :: GAMMA2,GM2,GEXP,GEXP2,QQQ,QQM,QRM,QYM,QZM,QEE,QEM,QEMA,
     + PEE,AVFL,QP0,QT0,RMASS1,RMASS2,XKIN1,XKIN2,PK1,PK0,VELK1,
     + VELK0,RMAK1,RMAK0,P01,P00,T01,T00,HEM,HEMA,P0A,TOA,HEAD1,DPEE,
     + RLIP,QK1,CP,TEMPM,TH,QEH,TUSIGN,XKIN0,QQQ1,QQQ2,QRM1,QRM2,
     + RLKE,TEMPM1,TEMPM2,EVAPM1,EVAPMH1,QEE1,QEE2,FRSDEN,
     + DELTAP,QQQTOT,QIFM1,QIFM2,EVAPMH2,EVAPMC,PAVE
      INTEGER :: ITURB,IPRESC,ISTATE,IMAX,JMAX,KMAX,IN,JN,KN,ICYCLE,
     + ISTR,JSTR,KSTR,ISSB,I,J,K,K1,K0,INDEX
      CHARACTER*1 :: DIR
      CHARACTER*10 :: SOLUTION_TYPE
      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)

      WRITE(44,9000) DIR,ICYCLE
      WRITE(44,*) 'mf=massflow, gm=gamma'
      WRITE(44,*) 'HEAD1 = p_ave/(ro_ave*g) + 0.5*ro_ave*(mf_tot/'
     +          // '(ro_ave*A))**2 [p=mf.av.,ro=area av.]'
      WRITE(44,*) 'p_0   = p * (1 + (gm-1)/2 * Ma**2)**(gm/(gm-1)'
     +          // ')              [mf. av.]'
      WRITE(44,*) 'T_0   = T / ((p/p_0)**(gm-1/gm))              '
     +          // '               [mf. av.]'
      WRITE(44,*)
      WRITE(44,9050) DIR
        
      GAMMA2 = 1.4
      GM2    = (GAMMA2-1.)/2.
      GEXP   =  GAMMA2/(GAMMA2-1.)
      GEXP2  = (GAMMA2-1.) / GAMMA2

CJO   THE LOOP FOR FIRST MCHECK-PACKET BEGINS HERE
      DO 1001 K = 1,KMAX+1
         QQQ   = 0.
         QQM   = 0.
         QRM   = 0.
         QYM   = 0.
         QZM   = 0.
         QEE   = 0.
         QEM   = 0.
         QEMA  = 0.
         PEE   = 0.
         AVFL  = 0.
         ROAVE = 0.
         QP0   = 0.
         QT0   = 0.
      DO 1000 J = 1,JMAX
      DO 1000 I = 1,IMAX

         K1     = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         K0     = K1 - KSTR

         QQQ    = QQQ + F1R(K1)  ! FLUX-BASED MASS FLOW
         RMASS1 = 0.
         RMASS1 = A1X(K1)*RM(K0) + A1Y(K1)*RN(K0) + A1Z(K1)*RW(K0)
         RMASS2 = A1X(K1)*RM(K1) + A1Y(K1)*RN(K1) + A1Z(K1)*RW(K1)
         ROM    = RO(K1) + RO(K0)
C ... FLOW BASED ON THE AVERAGE VALUES
         QQM    = QQM + A1(K1)*(RMASS1+RMASS2 - UROT(K1)*ROM)*.5
         IF(.NOT.TWO_FLUIDL) THEN
         QRM    = QRM + F1RM(K1)
         QYM    = QYM + F1RN(K1)
         QZM    = QZM + F1RW(K1)         
         ELSEIF(TWO_FLUIDL) THEN ! Add pressure
            PAVE   = A1(K1)*(PDIFF(K1) + PDIFF(K0))*.5
            QRM    = QRM + VAR(K1)%FRM(1)+VAR(K1)%FRM(2) + A1X(K1)*PAVE
            QYM    = QYM + VAR(K1)%FRN(1)+VAR(K1)%FRN(2) + A1Y(K1)*PAVE
            QZM    = QZM + VAR(K1)%FRW(1)+VAR(K1)%FRW(2) + A1Z(K1)*PAVE
c         write(666,*) k,j,i,qym,f1rn(k1),(PDIFF(K1) + PDIFF(K0))*.5,
c     2   VAR(K1)%FRN(1),VAR(K1)%FRN(2)
         ENDIF
       
         QEE    = QEE + F1E(K1)
         XKIN1  = .5*(RM(K1)**2+RN(K1)**2+RW(K1)**2)/(RO(K1)+1.E-20)
         XKIN0  = .5*(RM(K0)**2+RN(K0)**2+RW(K0)**2)/(RO(K0)+1.E-20)

         PK1    = (P(K1) + XKIN1)/(RO(K1)+1.E-20)
         PK0    = (P(K0) + XKIN0)/(RO(K0)+1.E-20)
C ... insentropic total pressure
         VELK1  = (SQRT(2.*XKIN1)/(RO(K1)+1.E-20))
         VELK0  = (SQRT(2.*XKIN0)/(RO(K0)+1.E-20))
         RMAK1  = VELK1/(C(K1)+1.E-20)
         RMAK0  = VELK0/(C(K0)+1.E-20)
         P01    = P(K1)*(1.+ GM2*RMAK1**2)**GEXP
         P00    = P(K0)*(1.+ GM2*RMAK0**2)**GEXP
         QP0    = QP0 + F1R(K1)*(P01+P00)*.5
C ... isentropic total temperature
         T01    = TEMP(K1)*(1.+ GM2*RMAK1**2)
         T00    = TEMP(K0)*(1.+ GM2*RMAK0**2)
         QT0    = QT0 + F1R(K1)*(T01+T00)*.5
C ... mechanical energy flux ((P+.5*ro*u^2 * volumeflow)
         QEM    = QEM  + F1R(K1)*(P(K1)+P(K0)+XKIN1+XKIN0)/(ROM+1.E-20)
c        QEM    = QEM  + F1R(K1)*(P01+P00 + XKIN1 + XKIN0)/ROM
C ... total pressure p0  (area averaged)
         IF(ISTATE == 0) THEN 
            PK1 = P01           ! isentropic pressure for compressible flows
            PK0 = P00
         ENDIF
         QEMA   = QEMA + A1(K1)*(PK1 + PK0)*.5
C ... static pressure (area averaged)
         PEE    = PEE   + A1(K1)*(P(K1)+P(K0))*.5
c        PEE    = PEE   + F1R(K1)*(P(K1)+P(K0))*.5
C ... area of cross-section
         AVFL   = AVFL + A1(K1)
C ... density (area averaged)
         ROAVE=ROAVE+ROM*0.5*A1(K1)
 1000 CONTINUE
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
         DO J = 1,JMAX
         DO I = 1,IMAX
            K1     = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
            K0     = K1 - KSTR
            ROM    = RO(K1) + RO(K0)
            QEM    = QEM  + 5.*F1RK(K1)/3.
            PEE    = PEE  + F1R(K1)*(RK(K1)+RK(K0))/(ROM*2.+1.E-20)/3.
c           PEE    = PEE   + A1(K1)*(RK(K1)+RK(K0))/ROM*2./3.
         ENDDO
         ENDDO
      ENDIF

c        PEE    = PEE  / (QQQ+1.E-20)
         PEE    = PEE  / (AVFL+1.E-20)
         HEM    = QEM  /((QQQ +1.E-20)*G0)
         HEMA   = QEMA /((AVFL+1.E-20)*G0)
         ROAVE  = ROAVE/(AVFL+1.E-20)
         VEL    = QQQ  /(AVFL*ROAVE+1.E-20)
         P0A    = QP0/(QQQ+1.E-20)
         TOA    = QT0/(QQQ+1.E-20)
         HEAD1  = PEE/(ROAVE*G0+1.E-20)+(VEL*VEL)/(2.0*G0)
C   ... INLETIN ARVOT MUISTIIN ...
c         IF (K <= 1) THEN
c            PEE1 = PEE
c         END IF
         DELTAP  = FRSDEN*(GX*(XC(I) - GROUND) +
     &                     GY*(YC(I) - GROUND) +
     &                     GZ*(ZC(I) - GROUND))

         DPEE    = PEE - FRSPRE - DELTAP

1001  WRITE(44,9100) K+ISSB-1,QQQ,QEE,QRM,QYM,QZM,QEM,AVFL,HEAD1,HEM,
     +        DPEE,TOA,P0A
C            
CJO   THE LOOP FOR THE SECOND (MULTI_PHASE) MCHECK-PACKET BEGINS HERE
C...........................................................................

      IF(SOLUTION_TYPE == 'MULTI' .OR.
     +   SOLUTION_TYPE == 'CAVIT') THEN

      WRITE(44,9200) ('-',INDEX=1,140)
      WRITE(44,9070) DIR

         EVAPMC  = 0.

      DO K = 1,KMAX+1
         QQQ   = 0.
         QQQ1  = 0.
         QQQ2  = 0.
         QQM   = 0.
         QRM1  = 0.
         QRM2  = 0.
         QEE   = 0.
         QEE1  = 0.
         QEE2  = 0.
         TEMPM1= 0.
         TEMPM2= 0.
         EVAPM1= 0.
         EVAPMH1 = 0.
         EVAPMH2 = 0.
         QIFM1 = 0.
         QIFM2 = 0.
         AVFL  = 0.
      DO J = 1,JMAX
      DO I = 1,IMAX
         K1     = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
         K0     = K1 - KSTR
         QQQ    = QQQ  + F1R(K1)         ! FLUX-BASED MASS FLOW
         QQQ1   = QQQ1 + VAR(K1)%FRO(1)  ! PHASE 1
         QQQ2   = QQQ2 + VAR(K1)%FRO(2)  ! PHASE 2
         ROM    = RO(K1) + RO(K0)
C ... FLOW BASED ON THE AVERAGE VALUES
         QQM    = QQM  + F1R(K1)*2./ROM
         QRM1   = QRM1 + A1(K1)*VAR(K1)%ALFA(1)
         QRM2   = QRM2 + A1(K1)*VAR(K1)%ALFA(2)
         QEE    = QEE  + F1E(K1)
         QEE1   = QEE1 + VAR(K1)%FE(1)
         QEE2   = QEE2 + VAR(K1)%FE(2)
C ... temperature
         T01    = TEMP(K1)*(1.+ GM2*RMAK1**2)
         T00    = TEMP(K0)*(1.+ GM2*RMAK0**2)
         QT0    = QT0 + F1R(K1)*(T01+T00)*.5
         TEMPM1 = TEMPM1 + A1(K1)*(PRO(K1)%DTEMP(1)+PRO(K0)%DTEMP(1))*.5
         TEMPM2 = TEMPM2 + A1(K1)*(PRO(K1)%DTEMP(2)+PRO(K0)%DTEMP(2))*.5
         AVFL   = AVFL + A1(K1)
C ... evaporation rate (volume based)
c         EVAPM1 = EVAPM1+.5*A1(K1)*(VAR(K1)%EVAPR(1)+VAR(K0)%EVAPR(1))

         EVAPM1 = EVAPM1+.5*(VOL(K1)*VAR(K1)%EVAPR(1)
     +                     + VOL(K0)*VAR(K0)%EVAPR(1))
         EVAPMH1= EVAPMH1+.5*(VOL(K1)*VAR(K1)%EVAPR(1)*PRO(K1)%HSAT(1)
     +                      + VOL(K0)*VAR(K0)%EVAPR(1)*PRO(K0)%HSAT(1))
         EVAPMH2= EVAPMH2+.5*(VOL(K1)*VAR(K1)%EVAPR(2)*PRO(K1)%HSAT(2)
     +                      + VOL(K0)*VAR(K0)%EVAPR(2)*PRO(K0)%HSAT(2))
         QIFM1  = QIFM1 +.5*(VOL(K1)*PRO(K1)%QIF(1)
     +                     + VOL(K0)*PRO(K0)%QIF(1))
         QIFM2  = QIFM2 +.5*(VOL(K1)*PRO(K1)%QIF(2)
     +                     + VOL(K0)*PRO(K0)%QIF(2))

      ENDDO
      ENDDO ! I and J (mersu)
          
         EVAPMC = EVAPMC + EVAPM1
         QRM1   = QRM1/(AVFL+1.E-20)
         QRM2   = QRM2/(AVFL+1.E-20)
         TEMPM1 = TEMPM1/(AVFL+1.E-20)
         TEMPM2 = TEMPM2/(AVFL+1.E-20)
c         EVAPM1 = EVAPM1/(AVFL+1.E-20)

      WRITE(44,9100) K+ISSB-1,QQQ1,QQQ2,QQM,QEE1,QEE2,QRM1,QRM2,TEMPM1,
     +        TEMPM2,EVAPM1,EVAPMC,EVAPMH1,EVAPMH2,QIFM1,QIFM2
      ENDDO ! K index
      WRITE(44,*) '(*) Evaporation rate is lagging by one cycle '//
     + 'and is zero after restart'
      ENDIF ! SOLUTION_TYPE == MULTI .OR. CAVIT
C    
C ...........................................................................
C     ANOTHER PACKET OF SAME KIND OF DATA. THIS COMES BELOW THE PREVIOUS
C     PACKET.
      WRITE(44,9200) ('-',INDEX=1,140)
      WRITE(44,*) 'mf=massflow'
      WRITE(44,*) 'TU       = TUsign*SQRT(2./3.*ABS(k_ave)) / vel_ave '
     +            //'[vel=area av., k=mf.av.]'
      WRITE(44,*) 'ENTHALPY = e + p/ro +v^2/2 + k                     '
     +            //'[mf.av.]'
      WRITE(44,*) 'TEMP     = (mf_cell*h) / (cp *mf_tot)              '
     +            //'[mf.av.]'
      WRITE(44,*) 'HSTAT    = (p_ave-p_in) / (ro*g) - hstat_inlet     '
     +            //'[mf. av.]'
      WRITE(44,*) 'PDIFF    = p - p_inf                               '
     +            //'[mf. av.]'
      WRITE(44,9065) DIR
      DISTA = 0.
      XOLD  = 0.
      YOLD  = 0.
      ZOLD  = 0.

C ... check the turbulence model
      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
         RLKE=1.0
      ELSE 
         RLKE=0.0
      ENDIF

C ... check if the incompressible solver is used
      IF (IPRESC > 0) THEN
         RLIP=1.0
      ELSE
         RLIP=0.0
      ENDIF

      DO 2001 K = 1,KMAX+1
         QQQ    = 0.
         QEE    = 0.
         QK1    = 0.
         CP     = 0.
         PEE    = 0.
         AVFL   = 0.
         ROAVE  = 0.
         XKAVE   = 0.
         TUAVE  = 0.
         XAAVE  = 0.
         VELAVE = 0.
         RKFLX  = 0.
         EPSFLX = 0.
c1       RFTEMP = 0.
         XAVE   = 0.
         YAVE   = 0.
         ZAVE   = 0.
         TEMPM  = 0.
         QQQTOT = 0.
         DO 2000 J = 1,JMAX
         DO 2000 I = 1,IMAX
            K1   = 1+(I-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
            K0 = K1 - KSTR
            ROM  = RO(K1) + RO(K0)
            RKM  = RK(K1) + RK(K0)
C ... FLOW BASED ON THE AVERAGE VALUES
C ...       mass flux
            QQQ  = QQQ + F1R(K1)
C ...       total energy flux (=mf*(e+p/ro + v^2/2 +k))
            QEE  = QEE + F1E(K1)
C ...       if incompressible solver is used, free-stream pressure will be added
            QEE  = QEE + RLIP * F1R(K1) * FRSPRE / (0.5*ROM)
C ...       massflow * Cp
            CP   = CP + F1R(K1) * CH(K1) * PR / (VIS(K1)+1.E-20)
C ...       flow of kinetic energy and k
            QK1  = QK1 + F1R(K1) * ((RM(K1)**2+RN(K1)**2+RW(K1)**2)/
     +                   (0.5*ROM**2) + RLKE * RKM/ROM)
            TEMPM  = TEMPM + A1(K1)*(TEMP(K1) + TEMP(K0))*.5
c           TEMPM  = TEMPM + F1R(K1)*(TEMP(K1) + TEMP(K0))*.5 ! tai sit ei
C ...       massflow averaged static pressure 
            PEE  = PEE   + A1(K1)*(P(K1)+P(K0))*.5
C           PEE  = PEE   + F1R(K1)*(P(K1)+P(K0))*.5 ! tai sit ei
C ...       flow area
            AVFL = AVFL+ A1(K1)
C ...       Area averaged density
            ROAVE= ROAVE+ ROM*0.5*A1(K1)
C ...       mass flow averaged k
            XKAVE= XKAVE + (RKM/(ROM+1.E-20))*F1R(K1)
C ...        velocities
            VELK1=(SQRT(RM(K1)**2+RN(K1)**2+RW(K1)**2)/
     +           (RO(K1)+1.E-20))
            VELK0=(SQRT(RM(K0)**2+RN(K0)**2+RW(K0)**2)/
     +           (RO(K0)+1.E-20))
C ...       Area-averaged turbulence intensity
            TUAVE  = TUAVE + A1(K1)*RKM/ROM
C ...       local Mach number 
            RMAK1  = VELK1/(C(K1)+1.E-20)
            RMAK0  = VELK0/(C(K0)+1.E-20)
            XAAVE  = XAAVE + 0.5*(RMAK1+RMAK0)*A1(K1)
C ...       average velocity
            VELAVE = VELAVE + 0.5*(VELK1+VELK0)*A1(K1)
C ...       k flux
            RKFLX  = RKFLX  + F1RK(K1)
C ...       epsilon flux
c           EPSFLX = EPSFLX + F1EPS(K1)
C ...       Total massflux accross the surface
            QQQTOT = QQQTOT + ABS(F1R(K1))
c1          RFTEMP = RFTEMP+F1R(K1)*(TEMP(K1)+TEMP(K0))*0.5
C ...       average coordinates
            XAVE   = XAVE + XCO(K1)*A1(K1)
            YAVE   = YAVE + YCO(K1)*A1(K1)
            ZAVE   = ZAVE + ZCO(K1)*A1(K1)
 2000    CONTINUE

C ...    static temperature from the enthalpy
         IF (ISTATE == 6) THEN
C          reference temperature is 0.01C in water internal energy
           TH    = 273.16 + (QEE-QK1) / (CP+1.E-20)
         ELSE
           TH    = (QEE-QK1) / (CP+1.E-20) 
         ENDIF
         TEMPM = TEMPM/(AVFL+1.E-20)
c         TEMPM = TEMPM/(QQQ+1.E-20)
C ...    Total enthalpy
         QEH   = QEE/(QQQ+1.E-20)
C ...    pressure
c         PEE   = PEE/(QQQ+1.E-20)
         PEE   = PEE/(AVFL+1.E-20)
C  ...   density
         ROAVE = ROAVE     / (AVFL+1.E-20)
C ...    k
         XKAVE = XKAVE / (QQQ+1.E-20)
C ...    Mach number
         XAAVE = XAAVE /(AVFL+1.E-20)
         VELAVE= VELAVE /(AVFL+1.E-20)
C ...    average coordinates
         XAVE  = XAVE / (AVFL+1.E-20)
         YAVE  = YAVE / (AVFL+1.E-20)
         ZAVE  = ZAVE / (AVFL+1.E-20)

c1         RFTEMP= RFTEMP / (QQQ+1.E-20)

C ...    THE NEGATIVE K IN GHOSTCELLS IS HANDLED SO THAT THERE TU < 0.
C ...    THIS IS NOT PHYSICALLY POSSIBLE BUT SQUAREROOT OF NEGATIVE
C ...    NUMBER NEVER OCCURS. 

         TUSIGN= SIGN(1.,TUAVE)
c         TUAVE = TUSIGN*SQRT(2./3.*ABS(XKAVE))/(VELAVE+1.E-20)
         TUAVE = TUAVE/(AVFL+1.E-20)
         TUAVE = TUSIGN*SQRT(2./3.*ABS(TUAVE))/(VELAVE+1.E-20)
         
c1         RFTEMP= RFTEMP - RFTEM1
c1    RFTEMP is replaced by TH in the following

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
c        WRITE(44,9100) K+ISSB-1,ROAVE,XAAVE,XKAVE,TUAVE,QEH,
c     +  RKFLX,EPSFLX,TH,XAVE,YAVE,ZAVE
c        WRITE(44,9100) K+ISSB-1,ROAVE,XAAVE,XKAVE,TUAVE,QEH,
c     +  RKFLX,EPSFLX,TEMPM,XAVE,YAVE,ZAVE
        WRITE(44,9100) K+ISSB-1,ROAVE,XAAVE,XKAVE,TUAVE,QEH,
     +  RKFLX,QQQTOT,TEMPM,XAVE,YAVE,ZAVE
      ELSE 
c        WRITE(44,9100) K+ISSB-1,ROAVE,XAAVE,0.,0.,QEH,
c     +  0.,0.,TH,XAVE,YAVE,ZAVE
        WRITE(44,9100) K+ISSB-1,ROAVE,XAAVE,0.,0.,QEH,
     +  0.,0.,TEMPM,XAVE,YAVE,ZAVE
      ENDIF

2001  CONTINUE
      RETURN

CJO ... FORMATS.........
9000  FORMAT('  CHECK IN ',A1,'-DIRECTION AFTER',I7,' CYCLES',/)
9050  FORMAT(3X,A1,4X,'MASS FLUX', 3X,'ENERGY FLUX',2X,'X MOM.FLUX ',2X,
     + 'Y MOM.FLUX ',2X,'Z MOM.FLUX ',2X,
     + 'MECH.E.FLUX',2X,'FLOW AREA',5X,'HEAD1',8X,'TOTHEAD_MF',3X,
     + ' PDIFF  ',5X,'T0 ISENTR',4X,'P0 ISENTR')
9060  FORMAT(3X,A1,6X,' RO      ', 2X,'   MA      ',2X,'    K      ',2X,
     +' TU(u/U_L)',3X,' ENTHALPY  ',2X,'  KFLX  ',5X,' EPSFLX  ',6X,
     +' TEMP ',7X,'XAVE',9X,'YAVE',9X,'ZAVE')
9065  FORMAT(3X,A1,6X,' RO      ', 2X,'   MA      ',2X,'    K      ',2X,
     +' TU(u/U_L)',3X,' ENTHALPY  ',2X,'  KFLX  ',5X,' ABSFLX  ',6X,
     +' TEMP ',7X,'XAVE',9X,'YAVE',9X,'ZAVE')
9070  FORMAT(3X,A1,3X,'MASS FLUX(1)',1X,'MASS FLUX(2)',2X,'VOLUM.FLUX',
     + 2X,'EN. FLUX(1)',2X,'EN. FLUX(2)',4X,'ALPHA(1)',5X,'ALPHA(2)',4X,
     + ' TEMP(1) ',4X,' TEMP(2) ',4X,'EVAP.RATE ',4X,'EVAP.CUM ',4X,
     + ' EVAP*HLS',4X,' EVAP*HGS',4X,'    QIL  ',4X,'   QIG    ')
9100  FORMAT(I4,1X,15(1X,E12.5))
9200  FORMAT(140(A1))
      END SUBROUTINE HCHECK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CARROT(U,V,W,XC,YC,ZC,OMEGA,A,B,C,CENX,CENY,CENZ,
     +     IMAX,JMAX,KMAX,IN,JN,KN,VAR,MPUREL,MPVREL,MPWREL,
     +     SOLUTION_TYPE,NPHASE,NGL)

C ... Absolute velocities -> relative velocities in rotating grids

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL, DIMENSION(*) :: U, V, W
      REAL, DIMENSION(*) :: XC, YC, ZC
      REAL, DIMENSION(NPHASE,*) :: MPUREL, MPVREL, MPWREL
      REAL    :: OMEGA, A, B, C, CENX, CENY, CENZ, DX, DY, DZ
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, IPHASE, NPHASE
      INTEGER :: ISTRID, JSTRID, KSTRID, I, J, K, II, JJ, KK, NGL
      CHARACTER(LEN=10) :: SOLUTION_TYPE

      TYPE(MPHASE_VARIABLES) VAR(*)

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = ISTRID*JSTRID

      IF(SOLUTION_TYPE /= 'MULTI') THEN

      DO K = 0,KMAX+1
         KK = (KN+K-1)*KSTRID
         DO J = 0,JMAX+1
            JJ = (JN+J-1)*ISTRID + IN + KK
            DO I = 0,IMAX+1
               II = I + JJ

               DX = XC(II)-CENX
               DY = YC(II)-CENY
               DZ = ZC(II)-CENZ

               U(II) = U(II) - OMEGA*(B*DZ-C*DY)
               V(II) = V(II) - OMEGA*(C*DX-A*DZ)
               W(II) = W(II) - OMEGA*(A*DY-B*DX)

            ENDDO
         ENDDO
      ENDDO

      ELSEIF(SOLUTION_TYPE == 'MULTI') THEN

      DO IPHASE = 1,NPHASE
      DO K = 0,KMAX+1
         KK = (KN+K-1)*KSTRID
         DO J = 0,JMAX+1
            JJ = (JN+J-1)*ISTRID + IN + KK
            DO I = 0,IMAX+1
               II = I + JJ

               DX = XC(II)-CENX
               DY = YC(II)-CENY
               DZ = ZC(II)-CENZ

               MPUREL(IPHASE,II) = VAR(II)%U(IPHASE) - OMEGA*(B*DZ-C*DY)
               MPVREL(IPHASE,II) = VAR(II)%V(IPHASE) - OMEGA*(B*DZ-C*DY)
               MPWREL(IPHASE,II) = VAR(II)%W(IPHASE) - OMEGA*(B*DZ-C*DY)

            ENDDO
         ENDDO
      ENDDO
      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CARROT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RMROT(RO,RM,RN,RW,U,V,W,E,IMAX,JMAX,KMAX,
     + IN,JN,KN,VAR,PRO,SOLUTION_TYPE,NPHASE)

      USE TYPE_ARRAYS

      REAL :: RO(*), RM(*), RN(*), RW(*), U(*), V(*), W(*), E(*)
      CHARACTER(LEN=10) :: SOLUTION_TYPE

      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(PROPERTIES)       :: PRO(*)

C ... NEW MOMENTUM COMPONENTS AND TOTAL ENERGY (FOR VISUALIZATION)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = ISTRID*JSTRID

      DO 1000 K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO 1000 J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO 1000 I = 0,IMAX+1
      ROAPU   = MAX(RO(I+JJ),1.E-5)
      E(I+JJ) = E(I+JJ) - .5/ROAPU*
     +          (RM(I+JJ)**2+RN(I+JJ)**2+RW(I+JJ)**2)
 1000 CONTINUE

      IF(SOLUTION_TYPE /= 'MULTI') THEN

      DO 2000 K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO 2000 J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO 2000 I = 0,IMAX+1
      E(I+JJ) = E(I+JJ) + .5*RO(I+JJ)*
     +          (U(I+JJ)**2 + V(I+JJ)**2 + W(I+JJ)**2)
      RM(I+JJ)= RO(I+JJ)*U(I+JJ)
      RN(I+JJ)= RO(I+JJ)*V(I+JJ)
      RW(I+JJ)= RO(I+JJ)*W(I+JJ)
 2000 CONTINUE

      ELSEIF(SOLUTION_TYPE == 'MULTI') THEN

      DO K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO I = 0,IMAX+1
      E(I+JJ) = E(I+JJ) + .5*RO(I+JJ)*
     +          (U(I+JJ)**2 + V(I+JJ)**2 + W(I+JJ)**2)
      RM(I+JJ)= 0.
      RN(I+JJ)= 0.
      RW(I+JJ)= 0.
      ENDDO; ENDDO; ENDDO

      DO IPHASE = 1,NPHASE

      DO K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO I = 0,IMAX+1
      RM(I+JJ)= RM(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%U(IPHASE)
      RN(I+JJ)= RN(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%V(IPHASE)
      RW(I+JJ)= RW(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%W(IPHASE)
      ENDDO; ENDDO; ENDDO
     
      ENDDO

      ENDIF ! MULTI

      RETURN
      END SUBROUTINE RMROT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NEWVEL(RO,RM,RN,RW,U,V,W,IMAX,JMAX,KMAX,IN,JN,KN,
     +   VAR,PRO,SOLUTION_TYPE,NPHASE)

      USE TYPE_ARRAYS

      REAL :: RO(*), U(*), V(*), W(*), RM(*), RN(*), RW(*)
      CHARACTER(LEN=10) :: SOLUTION_TYPE

      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(PROPERTIES)       :: PRO(*)

C ... CARTESIAN VELOCITIES ARE CALCULATED FROM THE MOMENTUM COMPONENTS

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = ISTRID*JSTRID

      IF(SOLUTION_TYPE /= 'MULTI') THEN

      DO 1000 K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO 1000 J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO 1000 I = 0,IMAX+1
C ... RECALCULATE VELOCITIES
      RM(I+JJ) = RO(I+JJ)*U(I+JJ)
      RN(I+JJ) = RO(I+JJ)*V(I+JJ)
      RW(I+JJ) = RO(I+JJ)*W(I+JJ)
 1000 CONTINUE

      ELSEIF(SOLUTION_TYPE == 'MULTI') THEN

      DO K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO I = 0,IMAX+1
      RM(I+JJ)= 0.
      RN(I+JJ)= 0.
      RW(I+JJ)= 0.
      ENDDO; ENDDO; ENDDO

      DO IPHASE = 1,NPHASE

      DO K = 0,KMAX+1
      KK      = (KN+K-1)*KSTRID
      DO J = 0,JMAX+1
      JJ      = (JN+J-1)*ISTRID + IN + KK
      DO I = 0,IMAX+1
      RM(I+JJ)= RM(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%U(IPASE)
      RN(I+JJ)= RN(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%V(IPASE)
      RW(I+JJ)= RW(I+JJ) + VAR(I+JJ)%ALFA(IPHASE)*PRO(I+JJ)%RO(IPHASE)*
     +                     VAR(I+JJ)%W(IPASE)
      ENDDO; ENDDO; ENDDO
      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE NEWVEL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUNDX(XC,YC,ZC,IMAX,JMAX,KMAX,M)

      USE NS3CO, ONLY : IN, JN, KN
        
      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)

      INTEGER :: IMAX, JMAX, KMAX, M, I, J, K, IREPEAT
      INTEGER :: ISTRID, JSTRID, IL, KA, JJ, KK
      INTEGER :: KM1, KM2, KP1, KP2

C ... COORDINATES FOR FREE STREAM VOLUMES
C ... IN THE PRESENT VERSION (5.2.1996) ONLY THE SECOND LINE
C ... IN THIS VERSION (19.11.2002) THE FIRST LINE FOR THE COARSE LEVELS

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID

      DO IREPEAT=1,3  ! Repeat three times in order to get corners right

      DO 6000 K = -1,KMAX+2      ! jmax
      KA      = (KN+K-1)*IL
      J       = JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 6000 I = -1,IMAX+2
      KK      = JJ + I
      KM1     = KK - ISTRID
      KM2     = KK - 2*ISTRID
      KP1     = KK + ISTRID
      KP2     = KK + 2*ISTRID

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KP1) = 2.*XC(KK)  - XC(KM1)
      YC(KP1) = 2.*YC(KK)  - YC(KM1)
      ZC(KP1) = 2.*ZC(KK)  - ZC(KM1)
      ENDIF ! M > 1

      XC(KP2) = 2.*XC(KP1) - XC(KK)
      YC(KP2) = 2.*YC(KP1) - YC(KK)
      ZC(KP2) = 2.*ZC(KP1) - ZC(KK)

 6000 CONTINUE

      DO 7000 K = -1,KMAX+2      ! jmin
      KA      = (KN+K-1)*IL
      J       = 1
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 7000 I = -1,IMAX+2
      KK      = JJ + I
      KM1     = KK - ISTRID
      KM2     = KK - 2*ISTRID
      KP1     = KK + ISTRID
      KP2     = KK + 2*ISTRID

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KM1) = 2.*XC(KK)  - XC(KP1)
      YC(KM1) = 2.*YC(KK)  - YC(KP1)
      ZC(KM1) = 2.*ZC(KK)  - ZC(KP1)
      ENDIF ! M > 1

      XC(KM2) = 2.*XC(KM1) - XC(KK)
      YC(KM2) = 2.*YC(KM1) - YC(KK)
      ZC(KM2) = 2.*ZC(KM1) - ZC(KK)
 7000 CONTINUE

      DO 6500 J = -1,JMAX+2      ! kmax
      K       = KMAX
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      DO 6500 I = -1,IMAX+2
      KK      = JJ + I
      KM1     = KK - IL
      KP1     = KK + IL
      KP2     = KK + 2*IL

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KP1) = 2.*XC(KK)  - XC(KM1)
      YC(KP1) = 2.*YC(KK)  - YC(KM1)
      ZC(KP1) = 2.*ZC(KK)  - ZC(KM1)
      ENDIF ! M > 1

      XC(KP2) = 2.*XC(KP1) - XC(KK)
      YC(KP2) = 2.*YC(KP1) - YC(KK)
      ZC(KP2) = 2.*ZC(KP1) - ZC(KK)

 6500 CONTINUE
      DO 7500 J = -1,JMAX+2      ! kmin
      K       = 1
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      DO 7500 I = -1,IMAX+2
      KK      = JJ + I
      KM1     = KK + IL
      KP1     = KK - IL
      KP2     = KK - 2*IL

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KP1) = 2.*XC(KK)  - XC(KM1)
      YC(KP1) = 2.*YC(KK)  - YC(KM1)
      ZC(KP1) = 2.*ZC(KK)  - ZC(KM1)
      ENDIF ! M > 1

      XC(KP2) = 2.*XC(KP1) - XC(KK)
      YC(KP2) = 2.*YC(KP1) - YC(KK)
      ZC(KP2) = 2.*ZC(KP1) - ZC(KK)

 7500 CONTINUE
      DO 8500 K = -1,KMAX+2      ! imax
      I       = IMAX
      DO 8500 J = -1,JMAX+2
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      KK      = JJ + I
      KM1     = KK - 1
      KP1     = KK + 1
      KP2     = KK + 2

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KP1) = 2.*XC(KK)  - XC(KM1)
      YC(KP1) = 2.*YC(KK)  - YC(KM1)
      ZC(KP1) = 2.*ZC(KK)  - ZC(KM1)
      ENDIF ! M > 1

      XC(KP2) = 2.*XC(KP1) - XC(KK)
      YC(KP2) = 2.*YC(KP1) - YC(KK)
      ZC(KP2) = 2.*ZC(KP1) - ZC(KK)

 8500 CONTINUE
      DO 9500 K = -1,KMAX+2      ! imin
      I       = 1
      DO 9500 J = -1,JMAX+2
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      KK      = JJ + I
      KM1     = KK + 1
      KP1     = KK - 1
      KP2     = KK - 2

      IF(M > 1) THEN ! Also the first ghost cell layer
      XC(KP1) = 2.*XC(KK)  - XC(KM1)
      YC(KP1) = 2.*YC(KK)  - YC(KM1)
      ZC(KP1) = 2.*ZC(KK)  - ZC(KM1)
      ENDIF ! M > 1

      XC(KP2) = 2.*XC(KP1) - XC(KK)
      YC(KP2) = 2.*YC(KP1) - YC(KK)
      ZC(KP2) = 2.*ZC(KP1) - ZC(KK)

 9500 CONTINUE

      ENDDO  ! IREPEAT

      RETURN
      END SUBROUTINE BOUNDX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SKEWNESS(XFC,YFC,ZFC,XC,YC,ZC,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     +     A3X,A3Y,A3Z,IMAX,JMAX,KMAX,M,NGL,BLKS,SDI)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IN, JN, KN
        
      IMPLICIT NONE

      TYPE(BLOCKS)   :: BLKS(*)
      TYPE(DIFF_COR) :: SDI(*)

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,KSTRID,IT0,IT5,
     +           I,J,K,IL,II,JJ,M,NGL

c      REAL XFC(*),YFC(*),ZFC(*),XC(*),YC(*),ZC(*),A1X(*),A1Y(*),A1Z(*),
      REAL :: XFC(*),YFC(*),ZFC(*),A1X(*),A1Y(*),A1Z(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      REAL :: SXI,SYI,SZI,SLEN,SKEWX,SKEWY,SKEWZ,SKEW,
     +        SKXMIN,SKYMIN,SKZMIN

      REAL :: XC(*), YC(*), ZC(*)

C ... Dot product between the surface normal and cell connecting vectors

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID

      SKEWX  = 0. ; SKEWY  = 0. ; SKEWZ  = 0.
      SKXMIN = 2. ; SKYMIN = 2. ; SKZMIN = 2.

      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTRID + JJ + IN
      DO I = 1,IMAX + 1
         IT0    = I + II
         IT5    = IT0 - 1
         SXI    =-XC(IT5) + XC(IT0)
         SYI    =-YC(IT5) + YC(IT0)
         SZI    =-ZC(IT5) + ZC(IT0)
         SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
         SXI    = SXI/SLEN
         SYI    = SYI/SLEN
         SZI    = SZI/SLEN
         XFC(IT0) = A1X(IT0)*SXI + A1Y(IT0)*SYI + A1Z(IT0)*SZI
         ZFC(IT0) = MIN(1.,ZFC(IT0))
         SKEWX  = SKEWX + XFC(IT0)
         SKXMIN = MIN(SKXMIN,XFC(IT0))
         SDI(IT0)%S1X = SXI ;  SDI(IT0)%S1Y = SYI ;  SDI(IT0)%S1Z = SZI
         SDI(IT0)%S1  = SLEN 
      ENDDO ; ENDDO ; ENDDO

      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX + 1
      II   = (JN+J-1)*ISTRID + JJ + IN
      DO I = 1,IMAX
         IT0    = I + II
         IT5    = IT0 - ISTRID
         SXI    =-XC(IT5) + XC(IT0)
         SYI    =-YC(IT5) + YC(IT0)
         SZI    =-ZC(IT5) + ZC(IT0)
         SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
         SXI    = SXI/SLEN
         SYI    = SYI/SLEN
         SZI    = SZI/SLEN
         YFC(IT0) = A2X(IT0)*SXI + A2Y(IT0)*SYI + A2Z(IT0)*SZI
         YFC(IT0) = MIN(1.,YFC(IT0))
         SKEWY  = SKEWY + YFC(IT0)
         SKYMIN = MIN(SKYMIN,YFC(IT0))
         SDI(IT0)%S2X = SXI ;  SDI(IT0)%S2Y = SYI ;  SDI(IT0)%S2Z = SZI
         SDI(IT0)%S2  = SLEN 
      ENDDO ; ENDDO ; ENDDO

      DO K = 1,KMAX + 1
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTRID + JJ + IN
      DO I = 1,IMAX
         IT0    = I + II
         IT5    = IT0 - IL
         SXI    =-XC(IT5) + XC(IT0)
         SYI    =-YC(IT5) + YC(IT0)
         SZI    =-ZC(IT5) + ZC(IT0)
         SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
         SXI    = SXI/SLEN
         SYI    = SYI/SLEN
         SZI    = SZI/SLEN
         ZFC(IT0) = A3X(IT0)*SXI + A3Y(IT0)*SYI + A3Z(IT0)*SZI
         ZFC(IT0) = MIN(1.,ZFC(IT0))
         SKEWZ  = SKEWZ + ZFC(IT0)
         SKZMIN = MIN(SKZMIN,ZFC(IT0))
         SDI(IT0)%S3X = SXI ;  SDI(IT0)%S3Y = SYI ;  SDI(IT0)%S3Z = SZI
         SDI(IT0)%S3  = SLEN 
      ENDDO ; ENDDO ; ENDDO

      SKEWX  = SKEWX/((IMAX+1)*JMAX*KMAX)
      SKEWY  = SKEWY/(IMAX*(JMAX+1)*KMAX)
      SKEWZ  = SKEWZ/(IMAX*JMAX*(KMAX+1))

*C ... Print out the skewness values
*
*      WRITE(45,'(2X,3(F7.5,5X),2X,3(F7.5,5X))') SKXMIN,SKYMIN,SKZMIN,
*     + SKEWX,SKEWY,SKEWZ

      BLKS(NGL)%SKEWI = SKEWX
      BLKS(NGL)%SKEWJ = SKEWY
      BLKS(NGL)%SKEWK = SKEWZ
      BLKS(NGL)%SKEWIMIN = SKXMIN
      BLKS(NGL)%SKEWJMIN = SKYMIN
      BLKS(NGL)%SKEWKMIN = SKZMIN

      RETURN
      END SUBROUTINE SKEWNESS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOUNDV(VOL,IMAX,JMAX,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: VOL(*)

C ... VOLUMES IN THE GHOST CELLS ARE FORCED TO BE POSITIVE
C ... VOI VOI

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      DO 6000 K = -1,KMAX+2
      KA      = (KN+K-1)*IL
      J       = JMAX
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 6000 I = -1,IMAX+2
      KK      = JJ + I
      KP1     = KK + ISTRID
      KP2     = KK + 2*ISTRID
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 6000 CONTINUE

      DO 7000 K = -1,KMAX+2
      KA      = (KN+K-1)*IL
      J       = 1
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 7000 I = -1,IMAX+2
      KK      = JJ + I
      KP1     = KK - ISTRID
      KP2     = KK - 2*ISTRID
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 7000 CONTINUE

      DO 6500 J = -1,JMAX+2
      K       = KMAX
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      DO 6500 I = -1,IMAX+2
      KK      = JJ + I
      KP1     = KK + IL
      KP2     = KK + 2*IL
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 6500 CONTINUE

      DO 7500 J = -1,JMAX+2
      K       = 1
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      DO 7500 I = -1,IMAX+2
      KK      = JJ + I
      KP1     = KK - IL
      KP2     = KK - 2*IL
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 7500 CONTINUE

      DO 8500 K = -1,KMAX+2
      I       = IMAX
      DO 8500 J = -1,JMAX+2
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      KK      = JJ + I
      KP1     = KK + 1
      KP2     = KK + 2
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 8500 CONTINUE

      DO 9500 K = -1,KMAX+2
      I       = 1
      DO 9500 J = -1,JMAX+2
      JJ      = (JN+J-1)*ISTRID + IN + (KN+K-1)*IL
      KK      = JJ + I
      KP1     = KK - 1
      KP2     = KK - 2
      VOL(KP1) = ABS(VOL(KP1))
      VOL(KP1) = MAX(VOL(KP1),VOL(KK)*.01)
      VOL(KP2) = 2.*VOL(KP1) - VOL(KK)
      VOL(KP2) = ABS(VOL(KP2))
 9500 CONTINUE

      RETURN
      END SUBROUTINE BOUNDV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SWEEPR(DRO,DM,DN,DW,DE,RO,U,V,W,VOL,IMAX,JMAX,KMAX,
     2 OMEGA,OMEX,OMEY,OMEZ,IN,JN,KN,DT)
      REAL VOL(*),DRO(*),DM(*),DN(*),DE(*),DW(*),RO(*),U(*),V(*),W(*)
C
C ... ADD THE ROTATIONAL ACCELERATION TERM
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      IF(OMEX >= .9999) THEN ! old tradional rotation axis (x)
      DO 1000 K = 1,KMAX
      IA = (KN+K-1)*IL + JN*ISTRID + IN
         DO 1000 IJ = 1,IMAX*JMAX
               J1 = (IJ-1)/IMAX    ! J1 = J-1
               N  = IJ + J1*2*IN + IA
               DN(N)   = DN(N) + OMEGA*RO(N)*W(N)
               DW(N)   = DW(N) - OMEGA*RO(N)*V(N)

 1000    CONTINUE
      ELSE                      ! arbirary axis
      DO K = 1,KMAX
      IA = (KN+K-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
               J1 = (IJ-1)/IMAX    ! J1 = J-1
               N  = IJ + J1*2*IN + IA
               DM(N)   = DM(N)+ OMEGA*RO(N)*( OMEZ*V(N) - OMEY*W(N))
               DN(N)   = DN(N)+ OMEGA*RO(N)*( OMEX*W(N) - OMEZ*U(N))
               DW(N)   = DW(N)+ OMEGA*RO(N)*( OMEY*U(N) - OMEX*V(N))
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SWEEPR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SWEEPRMF(VAR,PRO,VOL,IMAX,JMAX,KMAX,
     2 OMEGA,OMEX,OMEY,OMEZ,IN,JN,KN,DT,NPHASE)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,IL,N,IA,IJ,J1,K,
     +           IPHASE,NPHASE
      REAL    :: UN,VN,WN,OMEGA,OMEX,OMEY,OMEZ,DT
      REAL    :: VOL(*)
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
C
C ... ADD THE ROTATIONAL ACCELERATION TERM
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      DO IPHASE = 1,NPHASE

      IF(OMEX >= .9999) THEN ! old tradional rotation axis (x)
      DO 1000 K = 1,KMAX
      IA = (KN+K-1)*IL + JN*ISTRID + IN
         DO 1000 IJ = 1,IMAX*JMAX
           J1 = (IJ-1)/IMAX    ! J1 = J-1
           N  = IJ + J1*2*IN + IA
           VAR(N)%DN(IPHASE) = VAR(N)%DN(IPHASE) + 
     +     OMEGA*VAR(N)%ALFA(IPHASE)*PRO(N)%RO(IPHASE)*VAR(N)%W(IPHASE)
           VAR(N)%DW(IPHASE) = VAR(N)%DW(IPHASE) - 
     +     OMEGA*VAR(N)%ALFA(IPHASE)*PRO(N)%RO(IPHASE)*VAR(N)%V(IPHASE)
 1000    CONTINUE

      ELSE                      ! arbirary axis

      DO K = 1,KMAX
      IA = (KN+K-1)*IL + JN*ISTRID + IN
         DO IJ = 1,IMAX*JMAX
           J1 = (IJ-1)/IMAX    ! J1 = J-1
           N  = IJ + J1*2*IN + IA
           UN = VAR(N)%U(IPHASE)
           VN = VAR(N)%V(IPHASE)
           WN = VAR(N)%W(IPHASE)
           VAR(N)%DM(IPHASE) = VAR(N)%DM(IPHASE) + 
     +     OMEGA*VAR(N)%ALFA(IPHASE)*PRO(N)%RO(IPHASE)*(OMEZ*VN-OMEY*WN)
           VAR(N)%DW(IPHASE) = VAR(N)%DW(IPHASE) + 
     +     OMEGA*VAR(N)%ALFA(IPHASE)*PRO(N)%RO(IPHASE)*(OMEX*WN-OMEZ*UN)
           VAR(N)%DW(IPHASE) = VAR(N)%DW(IPHASE) + 
     +     OMEGA*VAR(N)%ALFA(IPHASE)*PRO(N)%RO(IPHASE)*(OMEY*UN-OMEX*VN)
         ENDDO
       ENDDO
      ENDIF

      ENDDO

      RETURN
      END SUBROUTINE SWEEPRMF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTRE(UU,UV,UW,VV,VW,WW,DUU,DUV,DUW,DVV,DVW,DWW,
     +     VOL,IMAX,JMAX,KMAX,OMEGA,OMEX,OMEY,OMEZ,IN,JN,KN)
      REAL UU(*),UV(*),UW(*),VV(*),VW(*),WW(*),
     +     DUU(*),DUV(*),DUW(*),DVV(*),DVW(*),DWW(*),VOL(*)
C
C ... ADD THE ROTATIONAL ACCELERATION TERM for Reynolds stresses
C ... x is the rotational axis.
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID

      IF(OMEX <=  .999) THEN
         WRITE(*,*) 'For a steady state Reynolds stresses the rotation'
         WRITE(*,*) 'axis must be x-axis, at least currently'
         STOP
      ENDIF

      DO 1000 K = 1,KMAX
      IA      = (KN+K-1)*IL + JN*ISTRID
      DO 1000 I = 1,JMAX*ISTRID
      N       = IA + I
C      DT1     = OMEGA/VOL(N)
C      DT1     = OMEGA

      DUV(N)  = DUV(N) + UW(N)        *OMEGA
      DUW(N)  = DUW(N) - UV(N)        *OMEGA
      DVV(N)  = DVV(N) + 2.*VW(N)     *OMEGA
      DVW(N)  = DVW(N) + (WW(N)-VV(N))*OMEGA
      DWW(N)  = DWW(N) - 2.*VW(N)     *OMEGA

1000  CONTINUE
      RETURN
      END SUBROUTINE ROTRE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTGRI(XCO,YCO,ZCO,XGRIG,YGRIG,ZGRIG,
     2 IMAX,JMAX,KMAX,ROTANG,IN,JN,KN,OMEX,OMEY,OMEZ,CENX,CENY,CENZ)
 
C ... THIS SUBROUTINE CALCULATES THE NEW GRID POINTS FROM
C ... THE ORIGINAL ONES AND FROM THE ROTATION ANGLE

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XCO,YCO,ZCO,XGRIG,YGRIG,ZGRIG
      REAL :: DROTAN, COSTHE, SINTHE, A, B, C, D
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: DX, DY, DZ, XC, YC, ZC
      REAL             :: ROTANG, OMEX, OMEY, OMEZ, CENX, CENY, CENZ
      INTEGER          :: IMAX, JMAX, KMAX, IN, JN, KN, N, NNODE

C ... A,  B,  C   Rotation axis vector
C ... XC, YC, ZC  Center of rotation
C ... DROTAN      Rotation angle

      DROTAN = -ROTANG
      A  = OMEX
      B  = OMEY
      C  = OMEZ
      XC = CENX
      YC = CENY
      ZC = CENZ
      D  = SQRT(A**2+B**2+C**2)
      A  = A/D
      B  = B/D
      C  = C/D

      COSTHE = COS(DROTAN)
      SINTHE = SIN(DROTAN)

C ... Rotation matrix

      A11 = A*A*(1.-COSTHE) +   COSTHE
      A12 = B*A*(1.-COSTHE) + C*SINTHE
      A13 = C*A*(1.-COSTHE) - B*SINTHE
      A21 = A*B*(1.-COSTHE) - C*SINTHE
      A22 = B*B*(1.-COSTHE) +   COSTHE
      A23 = C*B*(1.-COSTHE) + A*SINTHE
      A31 = A*C*(1.-COSTHE) + B*SINTHE
      A32 = B*C*(1.-COSTHE) - A*SINTHE
      A33 = C*C*(1.-COSTHE) +   COSTHE

      NNODE = (IMAX+2*IN)*(JMAX+2*JN)*(KMAX+2*KN)

      DO N = 1,NNODE

         DX = XGRIG(N)-XC
         DY = YGRIG(N)-YC
         DZ = ZGRIG(N)-ZC

         XCO(N) = A11*DX + A12*DY + A13*DZ + XC
         YCO(N) = A21*DX + A22*DY + A23*DZ + YC
         ZCO(N) = A31*DX + A32*DY + A33*DZ + ZC

      ENDDO

      RETURN
      END SUBROUTINE ROTGRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SLIMES(U,V,W,IMAX,JMAX,KMAX,N,NPATCH,ICON,RO,
     +  TEMP,P,PDIFF,PTUR,E,RM,RN,RW,RK,REPS,EPS2,VIST,FI,S11,BIJ,TIJ,
     +  PRO,VAR,TRM,MAXB,MAXEB,MAXSB,MAXSS,NSCAL,XCO,YCO,ZCO,M,ITURB,
     +  IPRO,MBPRO,TIMEL,MAXW,ISTRES,D1,D2,D3,VOL,
     +  DM,DN,DW,ISLIM,MULPHL,TRANSL)

C ... THIS SUBROUTINE HANDLES THE SLIDING SURFACES

C ... This subroutine manipulates the data in the ghost cells by
C ... shifting it in the circumferential direction. The sliding patches
C ... are connected before entering this subroutine.

C ... It is assumed here that the sliding patch covers the whole cyclic
C ... angle, i.e., the sliding patch is not allowed to be splitted and
C ... there must be only one block in the circumferential direction.

      USE MAIN_ARRAYS, ONLY : NPROCE, OMEGA, OMEGAX, OMEGAY, OMEGAZ,
     &                        CENAX, CENAY, CENAZ
      USE NS3CO,       ONLY : IN, JN, KN, ROTANG, IC9
      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : PII

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,MAXB,MAXEB,MAXSB,MAXSS,NSCAL,
     1 IPRO,MAXW,MBPRO,ISTRES,N,NPATCH,ITURB,ISLIM

      INTEGER :: NTOT,MAW,ISTRID,JSTRID,KSTRID,IL,IBC,IFACE,NBG,
     1 NBGCON,IXLO,IXUP,IYLO,IYUP,NPHI,NLE,NUE,I,J,K,L,IK,NC,IXS,IXE,
     2 IYS,IYE,ISTEP,ISTAR1,ISTAR2,NS,IPHASE,IP,M,
     3 IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2


      REAL :: ROTA1,ROTA2,SLIDE,CYCANG,DPHI,Y1,Y2,Y3,Y4,Y5,Y6,X1,X2,
     1 S1,S2,S3,A,B,C,ANGLEA,ANGLEB,ANGLEC,COSTHE,SINTHE,A11,A12,A13,
     2 A21,A22,A23,A31,A32,A33,B11,B12,B13,B21,B22,B23,B31,B32,B33,
     3 C11,C12,C13,C21,C22,C23,C31,C32,C33,ANGLED

      INTEGER :: ICON(IC9,*)

      REAL :: RAL
      REAL :: XCO(*), YCO(*), ZCO(*)
      REAL, DIMENSION(3) :: RC, RA, C1

      REAL :: RO(*),RM(*),RN(*),RW(*),E(*),
     &        P(*),PDIFF(*),U(*),V(*),W(*),TEMP(*),BIJ(MAXEB,*),PTUR(*),
     &        RK(*),REPS(*),FI(MAXSB,*),S11(MAXSS,*),EPS2(*),VIST(*),
     &        D1(*),D2(*),D3(*),VOL(*),DM(*),DN(*),DW(*),TIJ(*)

      REAL, ALLOCATABLE :: PHI(:), SAFE1(:), SAFE2(:), SAFE3(:),
     &                             SAFE4(:), SAFE5(:), SAFE6(:),
     &                             SA(:),SB(:),SC(:)

      LOGICAL :: CD, ACW
      LOGICAL :: AVER,TIMEL,MULPHL,TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      NTOT = (IMAX+2*IN)*(JMAX+2*JN)*(KMAX+2*KN)

      MAW = MAX((IMAX+2*IN)**2,(JMAX+2*JN)**2,(KMAX+2*KN)**2)      

      ALLOCATE(PHI(MAW),SAFE1(MAW),SAFE2(MAW),SAFE3(MAW),SAFE4(MAW),
     &                  SAFE5(MAW),SAFE6(MAW),SA(MAW),SB(MAW),SC(MAW))

      KSTRID = KMAX + 2*KN
      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      IL     = ISTRID*JSTRID

      DO 7000 IP = 1,NPATCH

      IBC    = ICON(1,IP)  ! PATCH TYPE
      IFACE  = ICON(3,IP)  ! FACE NUMBER

      IF(IBC == 11 .OR. IBC == 14) THEN ! SLIDING OR MIXING SURFACE

      NBG    = NPROCE(1+N,IPRO)                 ! global block number
      NBGCON = NPROCE(1+ICON(8,IP),ICON(17,IP)) ! g. connectivy block n.

      RA(1)  = OMEGAX(NBG)  ! Rotation axis
      RA(2)  = OMEGAY(NBG)
      RA(3)  = OMEGAZ(NBG)
      RAL = SQRT(RA(1)**2 + RA(2)**2 + RA(3)**2)

      IF(OMEGA(NBG) == 0.0) THEN
         RA(1) = OMEGAX(NBGCON) 
         RA(2) = OMEGAY(NBGCON)
         RA(3) = OMEGAZ(NBGCON)
         RAL = SQRT(RA(1)**2 + RA(2)**2 + RA(3)**2)
      ENDIF

      RA = RA/RAL

      ROTA1 = ROTANG(NBG)     ! ANGLE OF THIS BLOCK
      ROTA2 = ROTANG(NBGCON)  ! ANGLE OF CONNECTIVE BLOCK

      SLIDE = ROTA2 - ROTA1    ! SLIDE ANGLE

      AVER = IBC == 14         ! Mixing surface

C ... SLD patch extension

      CALL PATCHE(ICON,IP,M,IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2)

      IXLO = IA1 
      IXUP = IM1 + 1    ! Node number (not cell number)
      IYLO = JA1
      IYUP = JM1 + 1    ! Node number (not cell number)

C ... Compute the cell length distribution on the sliding surface  
C ....in circumferential direction by solving the phi-angles.

      CALL COMPHI(PHI,NPHI,IMAX,JMAX,KMAX,IN,JN,KN,
     +            ISTRID,JSTRID,KSTRID,IL,NLE,NUE,
     +            IXLO,IXUP,IYLO,IYUP,XCO,YCO,ZCO,IFACE,NBG,RA,RC,C1,CD)

C ... For historical reasons we need to define 
C ... the ACW (anticlockwise) boolean. 

      ACW = DOT_PRODUCT(RA,C1) > 0.0

C ... Fix the cyclic 'ghost cells' (copy from the actual angles)
         
      CYCANG = PHI(NPHI-NUE) - PHI(1+NLE)   ! Size of the cyclic angle

      DO I=1,NLE
         PHI(I) = PHI(I+NPHI-NLE-NUE-1) - CYCANG
      ENDDO

      DO I=NPHI-NUE+1,NPHI
         PHI(I) = PHI(I-NPHI+NLE    +2) + CYCANG
      ENDDO

C ... Compute the flux shares between sliding surface elements.

      DPHI = AMOD(SLIDE,CYCANG)
      DPHI = DPHI + CYCANG
      DPHI = AMOD(DPHI,CYCANG)

      DO I=1,NPHI*NPHI
         SA(I) = 0.0
         SB(I) = 0.0
         SC(I) = 0.0
      ENDDO

      DO J=2+NLE,NPHI-NUE     ! Exclude patch extension, i.e. ghost cells

         Y1 = DPHI + PHI(J-1)
         Y2 = DPHI + PHI(J)
         Y3 = Y1 - CYCANG
         Y4 = Y2 - CYCANG
         Y5 = Y1 + CYCANG
         Y6 = Y2 + CYCANG

         DO I=2,NPHI

            X1 = PHI(I-1)
            X2 = PHI(I)
            S1 = MAX(MIN(X2,Y2)-MAX(X1,Y1),0.0)/(X2-X1) 
            S2 = MAX(MIN(X2,Y4)-MAX(X1,Y3),0.0)/(X2-X1)
            S3 = MAX(MIN(X2,Y6)-MAX(X1,Y5),0.0)/(X2-X1)

            IF(ACW) THEN
               IK = (NPHI-1)*(NPHI-1)-((J-2)*(NPHI-1)+I-1)+1
            ELSE
               IK = (J-2)*(NPHI-1)+I-1
            ENDIF

            SA(IK) = S1
            SB(IK) = S2
            SC(IK) = S3

         ENDDO

      ENDDO


C ... Number of cells in circumferential direction
      NC = NPHI - 1

C ... Compute next the rotation matrices. We need three matrices:
C ... One for rotating vectors and tensors in case SLIDE > CYCANG,
C ... second for doing the previous rotation plus swaping
C ... vectors over the cyclic angle and third for doing the previous 
C ... rotation but swapping over the cyclic angle into the opposite 
C ... direction. 

      A = RA(1)
      B = RA(2)
      C = RA(3)

      ANGLEA = DPHI - SLIDE
      IF(SLIDE < 0.0) ANGLEA = ANGLEA - CYCANG
      ANGLEB = ANGLEA -     SIGN(CYCANG,SLIDE)
      ANGLEC = ANGLEA +     SIGN(CYCANG,SLIDE)
      ANGLED = ANGLEA - 2.0*SIGN(CYCANG,SLIDE)

      IF(SLIDE >= 0.0) THEN
         COSTHE = COS(ANGLEA)
         SINTHE = SIN(ANGLEA)
      ELSE
         COSTHE = COS(ANGLEB)
         SINTHE = SIN(ANGLEB)
      ENDIF

      A11 = A*A*(1.-COSTHE) +   COSTHE
      A12 = B*A*(1.-COSTHE) + C*SINTHE
      A13 = C*A*(1.-COSTHE) - B*SINTHE
      A21 = A*B*(1.-COSTHE) - C*SINTHE
      A22 = B*B*(1.-COSTHE) +   COSTHE
      A23 = C*B*(1.-COSTHE) + A*SINTHE
      A31 = A*C*(1.-COSTHE) + B*SINTHE
      A32 = B*C*(1.-COSTHE) - A*SINTHE
      A33 = C*C*(1.-COSTHE) +   COSTHE

      IF(SLIDE >= 0.0) THEN
         COSTHE = COS(ANGLEB)
         SINTHE = SIN(ANGLEB)
      ELSE
         COSTHE = COS(ANGLEA)
         SINTHE = SIN(ANGLEA)
      ENDIF

      B11 = A*A*(1.0-COSTHE) +   COSTHE
      B12 = B*A*(1.0-COSTHE) + C*SINTHE
      B13 = C*A*(1.0-COSTHE) - B*SINTHE
      B21 = A*B*(1.0-COSTHE) - C*SINTHE
      B22 = B*B*(1.0-COSTHE) +   COSTHE
      B23 = C*B*(1.0-COSTHE) + A*SINTHE
      B31 = A*C*(1.0-COSTHE) + B*SINTHE
      B32 = B*C*(1.0-COSTHE) - A*SINTHE
      B33 = C*C*(1.0-COSTHE) +   COSTHE

      IF(SLIDE >= 0.0) THEN
         COSTHE = COS(ANGLEC)
         SINTHE = SIN(ANGLEC)
      ELSE
         COSTHE = COS(ANGLED)
         SINTHE = SIN(ANGLED)
      ENDIF

      C11 = A*A*(1.-COSTHE) +   COSTHE
      C12 = B*A*(1.-COSTHE) + C*SINTHE
      C13 = C*A*(1.-COSTHE) - B*SINTHE
      C21 = A*B*(1.-COSTHE) - C*SINTHE
      C22 = B*B*(1.-COSTHE) +   COSTHE
      C23 = C*B*(1.-COSTHE) + A*SINTHE
      C31 = A*C*(1.-COSTHE) + B*SINTHE
      C32 = B*C*(1.-COSTHE) - A*SINTHE
      C33 = C*C*(1.-COSTHE) +   COSTHE

      IXLO = IA1 
      IXUP = IM1
      IYLO = JA1
      IYUP = JM1
 
      IF(CD) THEN
         IXS = IXLO
         IXE = IXUP
         IYS = IYLO
         IYE = IYUP
      ELSE
         IXS = IYLO
         IXE = IYUP
         IYS = IXLO
         IYE = IXUP
      ENDIF


C ... Starting addresses and step sizes. Two ghost cells in
C ... active direction will be updated.

      DO L=IYS,IYE

      IF(IFACE == 1 .OR. IFACE == 4) THEN

          IF(IFACE == 1) I = 0
          IF(IFACE == 4) I = IMAX+1

          IF(CD) THEN
            J = IXS
            K = L
            ISTEP  = ISTRID
          ELSE
            J = L
            K = IXS
            ISTEP  = IL
          ENDIF  

          ISTAR1 = IN + I + (JN-1+J)*ISTRID + (KN-1+K)*IL
          IF(IFACE == 1) ISTAR2 = ISTAR1 - 1
          IF(IFACE == 4) ISTAR2 = ISTAR1 + 1

        ELSE IF(IFACE == 2 .OR. IFACE == 5) THEN

          IF(IFACE == 2) J = 0
          IF(IFACE == 5) J = JMAX+1

          IF(CD) THEN
            I = IXS
            K = L
            ISTEP  = 1
          ELSE
            I = L
            K = IXS
            ISTEP  = IL
          ENDIF  

          ISTAR1 = IN + I + (JN-1+J)*ISTRID + (KN-1+K)*IL
          IF(IFACE == 2) ISTAR2 = ISTAR1 - ISTRID
          IF(IFACE == 5) ISTAR2 = ISTAR1 + ISTRID
 
        ELSE IF(IFACE == 3 .OR. IFACE == 6) THEN
 
          IF(IFACE == 3) K = 0
          IF(IFACE == 6) K = KMAX+1
  
          IF(CD) THEN
            I = IXS
            J = L
            ISTEP  = 1
          ELSE
            I = L
            J = IXS
            ISTEP  = ISTRID
          ENDIF 

          ISTAR1 = IN + I + (JN-1+J)*ISTRID + (KN-1+K)*IL
          IF(IFACE == 3) ISTAR2  = ISTAR1 - IL
          IF(IFACE == 6) ISTAR2  = ISTAR1 + IL
 
        ENDIF

C ... CALLED IN BOUNDI (ISLIM=1) AND OUTP (ISLIM=2)
       IF(ISLIM == 1 .OR. ISLIM == 2) THEN
        CALL SFTVEC(U,V,W,
     +              A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +              B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +              C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +              ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     +              PHI,AVER,ACW,RA)
        CALL SFTVEC(U,V,W,
     +              A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +              B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +              C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +              ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     +              PHI,AVER,ACW,RA)

        CALL SFTVEC(RM,RN,RW,
     +              A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +              B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +              C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +              ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     +              PHI,AVER,ACW,RA)
        CALL SFTVEC(RM,RN,RW,
     +              A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +              B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +              C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +              ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     +              PHI,AVER,ACW,RA)
       ENDIF ! ISLIM = 1 OR 2

C ... DM,DN,DW ADDED. CALLED IN OUTP: ISLIM = 2
       IF(ISLIM == 2) THEN
         CALL SFTVEC(DM,DN,DW,
     +              A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +              B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +              C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +              ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     +              PHI,AVER,ACW,RA)
       ENDIF ! ISLIM = 2


        IF(NSCAL == 1 .AND. ITURB <= 19) THEN
          DO NS=1,NSCAL
          CALL SFTSCL(FI(1,NS),ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,
     +                PHI,AVER,ACW,RA)
          CALL SFTSCL(FI(1,NS),ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,
     +                PHI,AVER,ACW,RA)
          ENDDO
        ELSEIF(ITURB >= 20 .AND. NSCAL > 6) THEN
          DO NS=7,NSCAL
          CALL SFTSCL(FI(1,NS),ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,
     +                PHI,AVER,ACW,RA)
          CALL SFTSCL(FI(1,NS),ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,
     +                PHI,AVER,ACW,RA)
          ENDDO
        ENDIF

        IF(ITURB >= 10 .AND. ITURB <= 19) THEN !EARSM by Gatski
          CALL SFTTEN(S11(1,1),S11(1,2),S11(1,3),
     +                S11(1,4),S11(1,5),S11(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR1,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
          CALL SFTTEN(S11(1,1),S11(1,2),S11(1,3),
     +                S11(1,4),S11(1,5),S11(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR2,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
        ENDIF

        IF(ISTRES > 0) THEN ! EARSM by Wallin and Johansson
          CALL SFTTEN(S11(1,1),S11(1,2),S11(1,3),
     +                S11(1,4),S11(1,5),S11(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR1,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
          CALL SFTTEN(S11(1,1),S11(1,2),S11(1,3),
     +                S11(1,4),S11(1,5),S11(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR2,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
          CALL SFTTEN(BIJ(1,1),BIJ(1,2),BIJ(1,3),
     +                BIJ(1,4),BIJ(1,5),BIJ(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR1,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
          CALL SFTTEN(BIJ(1,1),BIJ(1,2),BIJ(1,3),
     +                BIJ(1,4),BIJ(1,5),BIJ(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR2,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
        ENDIF

        IF(ITURB >= 20) THEN ! RSM
          CALL SFTTEN(FI(1,1),FI(1,2),FI(1,3),
     +                FI(1,4),FI(1,5),FI(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR1,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
          CALL SFTTEN(FI(1,1),FI(1,2),FI(1,3),
     +                FI(1,4),FI(1,5),FI(1,6),
     +                A11,A12,A13,A21,A22,A23,A31,A32,A33,
     +                B11,B12,B13,B21,B22,B23,B31,B32,B33,
     +                C11,C12,C13,C21,C22,C23,C31,C32,C33,
     +                ISTAR2,ISTEP,NC,SA,SB,SC,
     +                SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     +                PHI,AVER,ACW,RA)
        ENDIF

C ... D1,D2,D3,VOL ADDED.
C ... CALLED IN GRIVAL: ISLIM = 3
       IF(ISLIM == 3) THEN
        CALL SFTSCL(D1,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(D1,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
	
        CALL SFTSCL(D2,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(D2,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(D3,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(D3,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)	
	
        CALL SFTSCL(VOL,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
C        CALL SFTSCL(VOL,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW)	
       ENDIF ! ISLIM = 3 


C ... CALLED IN BOUNDI (ISLIM=1) AND OUTP (ISLIM=2)
       IF(ISLIM == 1 .OR. ISLIM == 2) THEN 
        CALL SFTSCL(TEMP,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(TEMP,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(RO,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(RO,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(P,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(P,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(PDIFF,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &              AVER,ACW,RA)
        CALL SFTSCL(PDIFF,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &              AVER,ACW,RA)

C ... PTUR ADDED
        CALL SFTSCL(PTUR,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(PTUR,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)	
	
        CALL SFTSCL(E,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(E,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(EPS2,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(EPS2,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)

        CALL SFTSCL(VIST,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
        CALL SFTSCL(VIST,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,AVER,ACW,RA)
       ENDIF ! ISLIM = 1 OR 2



        IF(ITURB >= 3 .AND. ITURB /= 8) THEN
           CALL SFTSCL(RK,  ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
           CALL SFTSCL(RK,  ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
           CALL SFTSCL(REPS,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
           CALL SFTSCL(REPS,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
           CALL SFTSCL(TIJ,ISTAR1,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
           CALL SFTSCL(TIJ,ISTAR2,ISTEP,NC,SA,SB,SC,SAFE1,PHI,
     &                 AVER,ACW,RA)
        ENDIF

        IF(MULPHL) THEN ! Multiphase model
          DO IPHASE = 1,NPHASES
          CALL SFTSCL(PRO(1:NTOT)%TEMP(IPHASE),ISTAR1,ISTEP,NC,SA,SB,SC,
     &    SAFE1,PHI,AVER,ACW,RA) ! Explicit use of these array maybe slow!!
          CALL SFTSCL(PRO(1:NTOT)%DTEMP(IPHASE),ISTAR1,ISTEP,NC,SA,SB,
     &    SC,SAFE1,PHI,AVER,ACW,RA)
          CALL SFTSCL(PRO(1:NTOT)%RO(IPHASE),  ISTAR1,ISTEP,NC,SA,SB,SC,
     &    SAFE1,PHI,AVER,ACW,RA)
          CALL SFTSCL(VAR(1:NTOT)%ALFA(IPHASE),ISTAR1,ISTEP,NC,SA,SB,SC,
     &    SAFE1,PHI,AVER,ACW,RA)
          CALL SFTSCL(VAR(1:NTOT)%X(IPHASE),   ISTAR1,ISTEP,NC,SA,SB,SC,
     &    SAFE1,PHI,AVER,ACW,RA)
          ENDDO
        ENDIF ! MULPHL

        IF(TRANSL) THEN ! Intermittency variables
           CALL SFTSCL(TRM(1:NTOT)%G,ISTAR1,ISTEP,NC,SA,SB,SC,
     &          SAFE1,PHI,AVER,ACW,RA)
           CALL SFTSCL(TRM(1:NTOT)%RET,ISTAR1,ISTEP,NC,SA,SB,
     &          SC,SAFE1,PHI,AVER,ACW,RA)
        ENDIF ! TRANSL

      ENDDO

      ENDIF   ! BC type 11 or 14 (= sliding or mixing surfaces)

 7000 CONTINUE

      DEALLOCATE(PHI,SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,SA,SB,SC)

      RETURN
      END SUBROUTINE SLIMES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COMPHI(PHI,NPHI,IMAX,JMAX,KMAX,IN,JN,KN,
     &                  ISTRID,JSTRID,KSTRID,IL,NLE,NUE,IXLO,IXUP,IYLO,
     &                  IYUP,XCO,YCO,ZCO,IFACE,NBG,RA,RC,C1,CD1)

C ... This subroutine computes the phi-angles for nodes on arc in 
C ... circumferential direction. First thing to do is to solve which 
C ... direction is the circumferential direction. 

      IMPLICIT NONE

      REAL, EXTERNAL :: PHIANG

      REAL :: PHI(*)

      INTEGER :: NPHI,IMAX,JMAX,KMAX,IN,JN,KN,ISTRID,JSTRID,KSTRID,IL
      INTEGER :: I, J, K, IPHI, NLE,NUE,IXLO,IXUP,IYLO,IYUP,IFACE,NBG
      INTEGER :: N0, N1, N2, N3, M1, M2, M3
      INTEGER :: IX2STEP, IX4STEP, IY2STEP, IY4STEP

      REAL :: XCO(*),YCO(*),ZCO(*)
      REAL :: X1, X2, R1, R2, VLEN1, VLEN2
      REAL, DIMENSION(3) :: C1, RA, RC
      REAL, DIMENSION(3) :: P0, P1

      LOGICAL :: CD1, CD2

      IX2STEP = (IXLO+IXUP)/2 
      IY2STEP = (IYLO+IYUP)/2
      IX4STEP = MAX(1,(IXLO+IXUP)/4)
      IY4STEP = MAX(1,(IYLO+IYUP)/4)

      IF(IFACE == 1 .OR. IFACE == 4) THEN

         IF(IFACE == 1) I = 1      
         IF(IFACE == 4) I = IMAX + 1      

C ... Which direction is the circumferential direction (CD) ?

         J  = IX2STEP           ! A node somewhere inside the patch
         K  = IY2STEP
         N2 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
         N3 = N2 + ISTRID*IX4STEP
         N1 = N2 - ISTRID*IX4STEP
         M1 = N1 + IL*IY4STEP
         M2 = N2 + IL*IY4STEP
         M3 = N3 + IL*IY4STEP
         
         CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD1)
         
         IF(.NOT.CD1) THEN

            N3 = N2 + IL*IY4STEP
            N1 = N2 - IL*IY4STEP
            M1 = N1 + ISTRID*IX4STEP
            M2 = N2 + ISTRID*IX4STEP
            M3 = N3 + ISTRID*IX4STEP

            CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD2)

         ENDIF


         IF(CD1) THEN

            K  = (IYLO+IYUP)/2
            IPHI = 0
            PHI(1) = 0.0
            J = IXLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO J = IXLO,IXUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(J > IXLO) THEN
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IXLO      ! Size of the patch lower extension
            NUE = IXUP - JMAX - 1 ! Size of the patch upper extension

         ELSEIF(CD2) THEN

            J = (IXLO+IXUP)/2
            IPHI = 0
            PHI(1) = 0.0
            K = IYLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO K = IYLO,IYUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(K > IYLO) THEN 
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IYLO
            NUE = IYUP - KMAX - 1

         ELSE

            WRITE(*,*) 'SLD error in COMPHI: NBG, IFACE = ',NBG,IFACE 

         ENDIF

      ELSE IF(IFACE == 2 .OR. IFACE == 5) THEN

         IF(IFACE == 2) J = 1      
         IF(IFACE == 5) J = JMAX + 1      

C ... Which direction is the circumferential direction (CD) ?

         I  = IX2STEP           ! A node somewhere inside the patch
         K  = IY2STEP
         N2 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL 
         N3 = N2 + 1*IX4STEP
         N1 = N2 - 1*IX4STEP 
         M1 = N1 + IL*IY4STEP 
         M2 = N2 + IL*IY4STEP 
         M3 = N3 + IL*IY4STEP 

         CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD1)

         IF(.NOT.CD1) THEN
            
            N3 = N2 + IL*IY4STEP 
            N1 = N2 - IL*IY4STEP 
            M1 = N1 + 1*IX4STEP 
            M2 = N2 + 1*IX4STEP 
            M3 = N3 + 1*IX4STEP 

            CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD2)

         ENDIF


         IF(CD1) THEN

            K  = (IYLO+IYUP)/2
            IPHI = 0
            PHI(1) = 0.0
            I = IXLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO I = IXLO,IXUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(I > IXLO) THEN
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IXLO
            NUE = IXUP - IMAX - 1

         ELSEIF(CD2) THEN
            
            I = (IXLO+IXUP)/2
            IPHI = 0
            PHI(1) = 0.0
            K = IYLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO K = IYLO,IYUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(K > IYLO) THEN 
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IYLO
            NUE = IYUP - KMAX - 1

         ELSE

            WRITE(*,*) 'SLD error in COMPHI: NBG, IFACE = ',NBG,IFACE 

         ENDIF

      ELSE IF(IFACE == 3 .OR. IFACE == 6) THEN

         IF(IFACE == 3) K = 1      
         IF(IFACE == 6) K = KMAX + 1      

C ... Which direction is the circumferential direction (CD) ?

         I  = IX2STEP           ! A node somewhere inside the patch
         J  = IY2STEP
         N2 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
         N3 = N2 + 1*IX4STEP
         N1 = N2 - 1*IX4STEP
         M1 = N1 + ISTRID*IY4STEP
         M2 = N2 + ISTRID*IY4STEP
         M3 = N3 + ISTRID*IY4STEP

         CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD1)

         IF(.NOT.CD1) THEN

            N3 = N2 + ISTRID*IY4STEP
            N1 = N2 - ISTRID*IY4STEP
            M1 = N1 + 1*IX4STEP
            M2 = N2 + 1*IX4STEP
            M3 = N3 + 1*IX4STEP
           
            CALL SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,RA,C1,RC,CD2)

         ENDIF


         IF(CD1) THEN

            J  = (IYLO+IYUP)/2
            IPHI = 0
            PHI(1) = 0.0
            I = IXLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO I = IXLO,IXUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(I > IXLO) THEN
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IXLO
            NUE = IXUP - IMAX - 1

         ELSEIF(CD2) THEN
            
            I = (IXLO+IXUP)/2
            IPHI = 0
            PHI(1) = 0.0
            J = IYLO
            N0 = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL

            P0(1) = XCO(N0)
            P0(2) = YCO(N0)
            P0(3) = ZCO(N0)

            DO J = IYLO,IYUP
               N1   = IN + I + (J-1+JN)*ISTRID + (K-1+KN)*IL
               P1(1) = XCO(N1)
               P1(2) = YCO(N1)
               P1(3) = ZCO(N1)
               IPHI = IPHI + 1             
               IF(J > IYLO) THEN 
                  PHI(IPHI) = PHI(IPHI-1) + PHIANG(P0,P1,RC) 
                  P0 = P1
               ENDIF
            ENDDO

            NLE = 1 - IYLO
            NUE = IYUP - JMAX - 1

         ELSE

            WRITE(*,*) 'SLD error in COMPHI: NBG, IFACE = ',NBG,IFACE 

         ENDIF

      ENDIF

      NPHI = IPHI

      RETURN
      END SUBROUTINE COMPHI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SLDROTAXIS(N1,N2,N3,M1,M2,M3,XCO,YCO,ZCO,
     &                      RA,C1,Q,CD)

C ... Solve the rotation axis based on three (six) points.

      USE NS3CO, ONLY : CHLREF

      IMPLICIT NONE

      REAL :: XCO(*), YCO(*), ZCO(*)
      REAL :: EPS, VLEN1, VLEN2, VLEN3, VLEN4
      REAL, DIMENSION(3) :: C1, C2, RA, N, Q,
     &                                  P1, P2, P3, P4, P5, P6

      LOGICAL :: CD

      INTEGER :: N1, N2, N3, M1, M2, M3

      EPS = 1.0E-4

      P1(1) = XCO(N1)
      P1(2) = YCO(N1)
      P1(3) = ZCO(N1)

      P2(1) = XCO(N2)
      P2(2) = YCO(N2)
      P2(3) = ZCO(N2)

      P3(1) = XCO(N3)
      P3(2) = YCO(N3)
      P3(3) = ZCO(N3)

      P4(1) = XCO(M1)
      P4(2) = YCO(M1)
      P4(3) = ZCO(M1)

      P5(1) = XCO(M2)
      P5(2) = YCO(M2)
      P5(3) = ZCO(M2)

      P6(1) = XCO(M3)
      P6(2) = YCO(M3)
      P6(3) = ZCO(M3)

      CALL CROSS_PRODUCT(P3-P2,P2-P1,N)

      IF(SQRT(SUM(N**2)) == 0.0) THEN
         CD = .FALSE.
         RETURN
      ENDIF
      
      CALL CROSS_PRODUCT(P2-P1,P3-P1,C1)

      VLEN1 = SQRT(SUM(C1**2))

      C1 = C1/VLEN1

      CD = ABS(DOT_PRODUCT(C1,RA)) > 0.99

      IF(CD) THEN

         CALL CENROT(P1,P2,P3,C1)
         CALL CENROT(P4,P5,P6,C2)

         CALL CROSS_PRODUCT(C2-C1,RA,Q)
         IF(SQRT(SUM(Q**2)) > CHLREF*EPS) CD = .FALSE.  

      ENDIF
      
      IF(CD) CALL CENROT(P1,P2,P3,Q)  ! SLD rotation center.

      RETURN
      END SUBROUTINE SLDROTAXIS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CROSS_PRODUCT(AV,BV,CV)

      IMPLICIT NONE

      REAL, DIMENSION(3), INTENT(IN)  :: AV, BV
      REAL, DIMENSION(3), INTENT(OUT) :: CV
 
      CV(1) = AV(2)*BV(3) - AV(3)*BV(2)
      CV(2) = AV(3)*BV(1) - AV(1)*BV(3)
      CV(3) = AV(1)*BV(2) - BV(1)*AV(2)

      RETURN
      END SUBROUTINE CROSS_PRODUCT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CENROT(P1,P2,P3,C)

      REAL, DIMENSION(3), INTENT(IN)  :: P1, P2, P3
      REAL, DIMENSION(3), INTENT(OUT) :: C
      REAL, DIMENSION(3) :: P4, P5, R1, N
      
      CALL CROSS_PRODUCT(P3-P2,P2-P1,N)
      CALL CROSS_PRODUCT(P2-P1,N,P5)

      R1 = P2 + 0.5*(P3-P2)

      P4 = P1 + 0.5*(P2-P1)

      C = P4 + DOT_PRODUCT((P3-P2),R1-P4)/DOT_PRODUCT((P3-P2),P5)*P5
      
      RETURN
      END SUBROUTINE CENROT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      REAL FUNCTION PHIANG(P0,P1,RC)

      REAL, DIMENSION(3) :: P0, P1, RC, A, B
      REAL :: VLA, VLB

      A = P0 - RC  
      B = P1 - RC  

      VLA = SQRT(A(1)**2 + A(2)**2 + A(3)**2)
      VLB = SQRT(B(1)**2 + B(2)**2 + B(3)**2)

      A = A/VLA
      B = B/VLB

      PHIANG = ACOS(DOT_PRODUCT(A,B))

      RETURN
      END FUNCTION PHIANG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C				
      SUBROUTINE SFTSCL(Q,ISTART,ISTEP,NC,SA,SB,SC,SAFE,PHI,AVER,ACW,RA)

C ... This subroutine shifts scalar data on one ghost cell row in
C ... circumferential direction.

      IMPLICIT NONE

      REAL, DIMENSION(3) :: RA

      REAL :: Q(*), SAFE(*)
      REAL :: SA(*), SB(*), SC(*), PHI(*)
      REAL :: A, B, C, SHARE, AV, PK
      
      INTEGER :: I, J, N, ISTART, ISTEP, NC

      LOGICAL :: AVER, ACW

      A = RA(1)
      B = RA(2)
      C = RA(3)

C ... Save the original data

      N = ISTART-ISTEP
      DO I=1,NC
        N = N + ISTEP
        SAFE(I) = Q(N)
      ENDDO   

C ... Update the row

      N = ISTART-ISTEP
      DO J=1,NC
        N = N + ISTEP
        Q(N) = 0.0
        DO I=1,NC
          SHARE = SA((I-1)*NC+J) + SB((I-1)*NC+J) + SC((I-1)*NC+J)
          Q(N) = Q(N) + SHARE*SAFE(I)
        ENDDO    
      ENDDO

      IF(AVER) THEN
C ... Compute average of the flow quantity
        N  = ISTART-ISTEP
        AV = 0.0
        DO I=1,NC
          N  = N  + ISTEP
          IF(ACW) THEN
            PK = (PHI(NC-I+2)-PHI(NC-I+1))/(PHI(NC+1)-PHI(1))
          ELSE
            PK = (PHI(I+1)-PHI(I))/(PHI(NC+1)-PHI(1))
          ENDIF
          AV = AV + PK*Q(N)
        ENDDO   
        N  = ISTART-ISTEP
        DO I=1,NC
          N = N  + ISTEP
          Q(N) = AV
        ENDDO 
      ENDIF  

      RETURN
      END SUBROUTINE SFTSCL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SFTVEC(Q1,Q2,Q3,
     &                  A11,A12,A13,A21,A22,A23,A31,A32,A33,
     &                  B11,B12,B13,B21,B22,B23,B31,B32,B33,
     &                  C11,C12,C13,C21,C22,C23,C31,C32,C33,
     &                  ISTART,ISTEP,NC,SA,SB,SC,SAFE1,SAFE2,SAFE3,
     &                  PHI,AVER,ACW,RA)

C ... This subroutine shifts vector data on one ghost cell row in
C ... circumferential direction.

      REAL, DIMENSION(3) :: RA

      REAL :: Q1(*), Q2(*), Q3(*), SAFE1(*), SAFE2(*), SAFE3(*) 
      REAL :: SA(*), SB(*), SC(*)

      REAL :: PHI(*)
      LOGICAL :: AVER, ACW

      A = RA(1)
      B = RA(2)
      C = RA(3)

C ... Save the original data

      N = ISTART-ISTEP
      DO I=1,NC
        N = N + ISTEP
        SAFE1(I) = Q1(N)
        SAFE2(I) = Q2(N)
        SAFE3(I) = Q3(N)
      ENDDO

C ... Update the row

      N = ISTART-ISTEP
      DO J=1,NC
        N = N + ISTEP
        Q1(N) = 0.0
        Q2(N) = 0.0
        Q3(N) = 0.0
        DO I=1,NC
          C1 = SAFE1(I)
          C2 = SAFE2(I)
          C3 = SAFE3(I)
          IF(SA((I-1)*NC+J) > 0.0) THEN
            UU = C1
            VV = C2
            WW = C3
            C1 = A11*UU + A21*VV + A31*WW
            C2 = A12*UU + A22*VV + A32*WW 
            C3 = A13*UU + A23*VV + A33*WW
          ELSEIF(SB((I-1)*NC+J) > 0.0) THEN
            UU = C1
            VV = C2
            WW = C3
            C1 = B11*UU + B21*VV + B31*WW
            C2 = B12*UU + B22*VV + B32*WW 
            C3 = B13*UU + B23*VV + B33*WW
          ELSEIF(SC((I-1)*NC+J) > 0.0) THEN
            UU = C1
            VV = C2
            WW = C3
            C1 = C11*UU + C21*VV + C31*WW
            C2 = C12*UU + C22*VV + C32*WW 
            C3 = C13*UU + C23*VV + C33*WW
          ENDIF
          SHARE = SA((I-1)*NC+J) + SB((I-1)*NC+J) + SC((I-1)*NC+J)
          Q1(N) = Q1(N) + SHARE*C1
          Q2(N) = Q2(N) + SHARE*C2
          Q3(N) = Q3(N) + SHARE*C3
        ENDDO    
      ENDDO

      IF(AVER) THEN
C ... Compute average of the flow quantity
        N  = ISTART-ISTEP
        AV1 = 0.0
        AV2 = 0.0
        AV3 = 0.0
        apk = 0.0
        DO I=1,NC
          N  = N  + ISTEP
          IF(ACW) THEN
            ANG = 0.5*(PHI(NC-I+2)+PHI(NC-I+1))
          ELSE
            ANG = 0.5*(PHI(I)+PHI(I+1))
          ENDIF
          COSTHE = COS(ANG)
          SINTHE = SIN(ANG)
          E11 = A*A*(1.-COSTHE) +   COSTHE
          E12 = B*A*(1.-COSTHE) + C*SINTHE
          E13 = C*A*(1.-COSTHE) - B*SINTHE
          E21 = A*B*(1.-COSTHE) - C*SINTHE
          E22 = B*B*(1.-COSTHE) +   COSTHE
          E23 = C*B*(1.-COSTHE) + A*SINTHE
          E31 = A*C*(1.-COSTHE) + B*SINTHE
          E32 = B*C*(1.-COSTHE) - A*SINTHE
          E33 = C*C*(1.-COSTHE) +   COSTHE
          IF(ACW) THEN
            PK = (PHI(NC-I+2)-PHI(NC-I+1))/(PHI(NC+1)-PHI(1))
          ELSE
            PK = (PHI(I+1)-PHI(I))/(PHI(NC+1)-PHI(1))
          ENDIF
          apk = apk + pk
          AV1 = AV1 + PK*(E11*Q1(N) + E12*Q2(N) + E13*Q3(N))
          AV2 = AV2 + PK*(E21*Q1(N) + E22*Q2(N) + E23*Q3(N))
          AV3 = AV3 + PK*(E31*Q1(N) + E32*Q2(N) + E33*Q3(N))
        ENDDO

        N  = ISTART-ISTEP
        DO I=1,NC
          N  = N  + ISTEP
          IF(ACW) THEN
            ANG = 0.5*(PHI(NC-I+2)+PHI(NC-I+1))
          ELSE
            ANG = 0.5*(PHI(I)+PHI(I+1))
          ENDIF
          COSTHE = COS(-ANG)
          SINTHE = SIN(-ANG)
          E11 = A*A*(1.-COSTHE) +   COSTHE
          E12 = B*A*(1.-COSTHE) + C*SINTHE
          E13 = C*A*(1.-COSTHE) - B*SINTHE
          E21 = A*B*(1.-COSTHE) - C*SINTHE
          E22 = B*B*(1.-COSTHE) +   COSTHE
          E23 = C*B*(1.-COSTHE) + A*SINTHE
          E31 = A*C*(1.-COSTHE) + B*SINTHE
          E32 = B*C*(1.-COSTHE) - A*SINTHE
          E33 = C*C*(1.-COSTHE) +   COSTHE
          ap1 = Q1(N)
          ap2 = Q2(N)
          ap3 = Q3(N)
          Q1(N) = E11*AV1 + E12*AV2 + E13*AV3
          Q2(N) = E21*AV1 + E22*AV2 + E23*AV3
          Q3(N) = E31*AV1 + E32*AV2 + E33*AV3
        ENDDO 
      ENDIF  

      RETURN
      END SUBROUTINE SFTVEC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SFTTEN(Q1,Q2,Q3,Q4,Q5,Q6,
     &                  A11,A12,A13,A21,A22,A23,A31,A32,A33,
     &                  B11,B12,B13,B21,B22,B23,B31,B32,B33,
     &                  C11,C12,C13,C21,C22,C23,C31,C32,C33,
     &                  ISTART,ISTEP,NC,SA,SB,SC,
     &                  SAFE1,SAFE2,SAFE3,SAFE4,SAFE5,SAFE6,
     &                  PHI,AVER,ACW,RA)

C ... This subroutine shifts tensor data, like Reynolds stresses,
C ... on one ghost cell row in circumferential direction.

      REAL, DIMENSION(3) :: RA

      REAL :: Q1(*), Q2(*), Q3(*), Q4(*), Q5(*), Q6(*) 
      REAL :: SAFE1(*),SAFE2(*),SAFE3(*),SAFE4(*),SAFE5(*),SAFE6(*)
      REAL :: SA(*), SB(*), SC(*)

      REAL :: PHI(*)
      LOGICAL :: AVER, ACW

      A = RA(1)
      B = RA(2)
      C = RA(3)

C ... Save the original data

      N = ISTART-ISTEP
      DO I=1,NC
        N = N + ISTEP
        SAFE1(I) = Q1(N)
        SAFE2(I) = Q2(N)
        SAFE3(I) = Q3(N)
        SAFE4(I) = Q4(N)
        SAFE5(I) = Q5(N)
        SAFE6(I) = Q6(N)
      ENDDO   

C ... Update the row

      N = ISTART-ISTEP
      DO J=1,NC
        N = N + ISTEP
        Q1(N) = 0.0
        Q2(N) = 0.0
        Q3(N) = 0.0
        Q4(N) = 0.0
        Q5(N) = 0.0
        Q6(N) = 0.0
        DO I=1,NC
          C1 = SAFE1(I)
          C2 = SAFE2(I)
          C3 = SAFE3(I)
          C4 = SAFE4(I)
          C5 = SAFE5(I)
          C6 = SAFE6(I)
          IF(SA((I-1)*NC+J) > 0.0) THEN
            UU  = C1
            UV  = C2
            UW  = C3
            VV  = C4
            VW  = C5
            WW  = C6
            D11 = A11*UU  + A12*UV  + A13*UW
            D12 = A11*UV  + A12*VV  + A13*VW
            D13 = A11*UW  + A12*VW  + A13*WW
            D21 = A21*UU  + A22*UV  + A23*UW
            D22 = A21*UV  + A22*VV  + A23*VW
            D23 = A21*UW  + A22*VW  + A23*WW
            D31 = A31*UU  + A32*UV  + A33*UW
            D32 = A31*UV  + A32*VV  + A33*VW
            D33 = A31*UW  + A32*VW  + A33*WW
            C1  = D11*A11 + D12*A12 + D13*A13
            C2  = D11*A21 + D12*A22 + D13*A23
            C3  = D11*A31 + D12*A32 + D13*A33
            C4  = D21*A21 + D22*A22 + D23*A23
            C5  = D21*A31 + D22*A32 + D23*A33
            C6  = D31*A31 + D32*A32 + D33*A33
          ELSEIF(SB((I-1)*NC+J) > 0.0) THEN
            UU  = C1
            UV  = C2
            UW  = C3
            VV  = C4
            VW  = C5
            WW  = C6
            D11 = B11*UU  + B12*UV  + B13*UW
            D12 = B11*UV  + B12*VV  + B13*VW
            D13 = B11*UW  + B12*VW  + B13*WW
            D21 = B21*UU  + B22*UV  + B23*UW
            D22 = B21*UV  + B22*VV  + B23*VW
            D23 = B21*UW  + B22*VW  + B23*WW
            D31 = B31*UU  + B32*UV  + B33*UW
            D32 = B31*UV  + B32*VV  + B33*VW
            D33 = B31*UW  + B32*VW  + B33*WW
            C1  = D11*A11 + D12*A12 + D13*A13
            C2  = D11*A21 + D12*A22 + D13*A23
            C3  = D11*A31 + D12*A32 + D13*A33
            C4  = D21*A21 + D22*A22 + D23*A23
            C5  = D21*A31 + D22*A32 + D23*A33
            C6  = D31*A31 + D32*A32 + D33*A33
          ELSEIF(SC((I-1)*NC+J) > 0.0) THEN
            UU  = C1
            UV  = C2
            UW  = C3
            VV  = C4
            VW  = C5
            WW  = C6
            D11 = C11*UU  + C12*UV  + C13*UW
            D12 = C11*UV  + C12*VV  + C13*VW
            D13 = C11*UW  + C12*VW  + C13*WW
            D21 = C21*UU  + C22*UV  + C23*UW
            D22 = C21*UV  + C22*VV  + C23*VW
            D23 = C21*UW  + C22*VW  + C23*WW
            D31 = C31*UU  + C32*UV  + C33*UW
            D32 = C31*UV  + C32*VV  + C33*VW
            D33 = C31*UW  + C32*VW  + C33*WW
            C1  = D11*A11 + D12*A12 + D13*A13
            C2  = D11*A21 + D12*A22 + D13*A23
            C3  = D11*A31 + D12*A32 + D13*A33
            C4  = D21*A21 + D22*A22 + D23*A23
            C5  = D21*A31 + D22*A32 + D23*A33
            C6  = D31*A31 + D32*A32 + D33*A33
          ENDIF
          SHARE = SA((I-1)*NC+J) + SB((I-1)*NC+J) + SC((I-1)*NC+J)
          Q1(N) = Q1(N) + SHARE*C1
          Q2(N) = Q2(N) + SHARE*C2
          Q3(N) = Q3(N) + SHARE*C3
          Q4(N) = Q4(N) + SHARE*C4
          Q5(N) = Q5(N) + SHARE*C5
          Q6(N) = Q6(N) + SHARE*C6
        ENDDO    
      ENDDO

      IF(AVER) THEN
C ... Compute average of the flow quantity
        N  = ISTART-ISTEP
        AV1 = 0.0
        AV2 = 0.0
        AV3 = 0.0
        AV4 = 0.0
        AV5 = 0.0
        AV6 = 0.0
        DO I=1,NC
          N  = N  + ISTEP
          IF(ACW) THEN
            ANG = 0.5*(PHI(NC-I+2)+PHI(NC-I+1))
          ELSE
            ANG = 0.5*(PHI(I)+PHI(I+1))
          ENDIF
          COSTHE = COS(ANG)
          SINTHE = SIN(ANG)
          E11 = A*A*(1.-COSTHE) +   COSTHE
          E12 = B*A*(1.-COSTHE) + C*SINTHE
          E13 = C*A*(1.-COSTHE) - B*SINTHE
          E21 = A*B*(1.-COSTHE) - C*SINTHE
          E22 = B*B*(1.-COSTHE) +   COSTHE
          E23 = C*B*(1.-COSTHE) + A*SINTHE
          E31 = A*C*(1.-COSTHE) + B*SINTHE
          E32 = B*C*(1.-COSTHE) - A*SINTHE
          E33 = C*C*(1.-COSTHE) +   COSTHE
          F11 = E11*Q1(N) + E12*Q2(N) + E13*Q3(N)
          F12 = E11*Q2(N) + E12*Q4(N) + E13*Q5(N)
          F13 = E11*Q3(N) + E12*Q5(N) + E13*Q6(N)
          F21 = E21*Q1(N) + E22*Q2(N) + E23*Q3(N)
          F22 = E21*Q2(N) + E22*Q4(N) + E23*Q5(N)
          F23 = E21*Q3(N) + E22*Q5(N) + E23*Q6(N)
          F31 = E31*Q1(N) + E32*Q2(N) + E33*Q3(N)
          F32 = E31*Q2(N) + E32*Q4(N) + E33*Q5(N)
          F33 = E31*Q3(N) + E32*Q5(N) + E33*Q6(N)
          IF(ACW) THEN
            PK = (PHI(NC-I+2)-PHI(NC-I+1))/(PHI(NC+1)-PHI(1))
          ELSE
            PK = (PHI(I+1)-PHI(I))/(PHI(NC+1)-PHI(1))
          ENDIF
          AV1 = AV1 + PK*(F11*E11 + F12*E12 + F13*E13)
          AV2 = AV2 + PK*(F11*E21 + F12*E22 + F13*E23)
          AV3 = AV3 + PK*(F11*E31 + F12*E32 + F13*E33)
          AV4 = AV4 + PK*(F21*E21 + F22*E22 + F23*E23)
          AV5 = AV5 + PK*(F21*E31 + F22*E32 + F23*E33)
          AV6 = AV6 + PK*(F31*E31 + F32*E32 + F33*E33)
        ENDDO
   
        N  = ISTART-ISTEP
        DO I=1,NC
          N  = N  + ISTEP
          IF(ACW) THEN
            ANG = 0.5*(PHI(NC-I+2)+PHI(NC-I+1))
          ELSE
            ANG = 0.5*(PHI(I)+PHI(I+1))
          ENDIF
          COSTHE = COS(-ANG)
          SINTHE = SIN(-ANG)
          E11 = A*A*(1.-COSTHE) +   COSTHE
          E12 = B*A*(1.-COSTHE) + C*SINTHE
          E13 = C*A*(1.-COSTHE) - B*SINTHE
          E21 = A*B*(1.-COSTHE) - C*SINTHE
          E22 = B*B*(1.-COSTHE) +   COSTHE
          E23 = C*B*(1.-COSTHE) + A*SINTHE
          E31 = A*C*(1.-COSTHE) + B*SINTHE
          E32 = B*C*(1.-COSTHE) - A*SINTHE
          E33 = C*C*(1.-COSTHE) +   COSTHE
          F11 = E11*AV1 + E12*AV2 + E13*AV3
          F12 = E11*AV2 + E12*AV4 + E13*AV5
          F13 = E11*AV3 + E12*AV5 + E13*AV6
          F21 = E21*AV1 + E22*AV2 + E23*AV3
          F22 = E21*AV2 + E22*AV4 + E23*AV5
          F23 = E21*AV3 + E22*AV5 + E23*AV6
          F31 = E31*AV1 + E32*AV2 + E33*AV3
          F32 = E31*AV2 + E32*AV4 + E33*AV5
          F33 = E31*AV3 + E32*AV5 + E33*AV6
          Q1(N) = F11*E11 + F12*E12 + F13*E13
          Q2(N) = F11*E21 + F12*E22 + F13*E23
          Q3(N) = F11*E31 + F12*E32 + F13*E33
          Q4(N) = F21*E21 + F22*E22 + F23*E23
          Q5(N) = F21*E31 + F22*E32 + F23*E33
          Q6(N) = F31*E31 + F32*E32 + F33*E33
        ENDDO 
      ENDIF  

      RETURN
      END SUBROUTINE SFTTEN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CELLCP(X,Y,Z,XC,YC,ZC,OMEGA,OMX,OMY,OMZ,CENX,CENY,CENZ,
     & UROTCP,VROTCP,WROTCP,IMAX,JMAX,KMAX,IN,JN,KN,PDIFF,
     &                  FRESUL,FRSDEN)
     
C ... This subroutine updates the cell center point coordinates after
C ... changes in the grid, e.g. rotation. Only one ghost cell in each
C ... direction will be updated. (ESa 4.2.1997)

      USE NS3CO, ONLY : GX, GY, GZ, GROUND
      
      REAL :: FRSDEN,OMEGA,OMX,OMY,OMZ,CENX,CENY,CENZ,XAVE,YAVE,ZAVE

      REAL :: PDIFF(*),UROTCP(*),VROTCP(*),WROTCP(*)

      REAL :: X(*),Y(*),Z(*),XC(*),YC(*),ZC(*)

      REAL :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &        X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,
     &        XCE,YCE,ZCE

      LOGICAL :: FRESUL, INSIDE

      ISTR =  1
      JSTR =  IMAX+2*IN
      KSTR = (IMAX+2*IN)*(JMAX+2*JN)

      DO 100 K = 0,KMAX+1
      DO 110 J = 0,JMAX+1
      DO 120 I = 0,IMAX+1

         INSIDE = (K /= 0 .AND. K /= KMAX+1 .AND. J /= 0 .AND.
     &             J /= JMAX+1 .AND. I /= 0 .AND. I /= IMAX+1)

         IT0 = 1+(I-1  +IN)*ISTR+(J-1  +JN)*JSTR+(K-1  +KN)*KSTR

         IT1 = 1+(I-1  +IN)*ISTR+(J-1  +JN)*JSTR+(K-1  +KN)*KSTR
         IT2 = 1+(I-1  +IN)*ISTR+(J-1  +JN)*JSTR+(K-1+1+KN)*KSTR
         IT3 = 1+(I-1  +IN)*ISTR+(J-1+1+JN)*JSTR+(K-1+1+KN)*KSTR
         IT4 = 1+(I-1  +IN)*ISTR+(J-1+1+JN)*JSTR+(K-1  +KN)*KSTR

         IT5 = 1+(I-1+1+IN)*ISTR+(J-1  +JN)*JSTR+(K-1  +KN)*KSTR
         IT6 = 1+(I-1+1+IN)*ISTR+(J-1  +JN)*JSTR+(K-1+1+KN)*KSTR
         IT7 = 1+(I-1+1+IN)*ISTR+(J-1+1+JN)*JSTR+(K-1+1+KN)*KSTR
         IT8 = 1+(I-1+1+IN)*ISTR+(J-1+1+JN)*JSTR+(K-1  +KN)*KSTR

C ... Cell corner point coordinates

         X1 = X(IT1)
         Y1 = Y(IT1)
         Z1 = Z(IT1)
         X2 = X(IT2)
         Y2 = Y(IT2)
         Z2 = Z(IT2)
         X3 = X(IT3)
         Y3 = Y(IT3)
         Z3 = Z(IT3)
         X4 = X(IT4)
         Y4 = Y(IT4)
         Z4 = Z(IT4)
                    
         X5 = X(IT5)
         Y5 = Y(IT5)
         Z5 = Z(IT5)
         X6 = X(IT6)
         Y6 = Y(IT6)
         Z6 = Z(IT6)
         X7 = X(IT7)
         Y7 = Y(IT7)
         Z7 = Z(IT7)
         X8 = X(IT8)
         Y8 = Y(IT8)
         Z8 = Z(IT8)

         XCE = .125*(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8)
         YCE = .125*(Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8)
         ZCE = .125*(Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8)

c         IF(FRESUL .AND. INSIDE) THEN
c            PDIFF(IT0) = PDIFF(IT0) + FRSDEN*(GX*(XCE-XC(IT0))
c     +                 +   GY*(YCE-YC(IT0)) + GZ*(ZCE-ZC(IT0)))  
c         ENDIF           

         XC(IT0) = XCE
         YC(IT0) = YCE
         ZC(IT0) = ZCE

C ... Velocities in the rotating coordinate system
         ZAVE = ZCE - CENZ 
         YAVE = YCE - CENY
         XAVE = XCE - CENX
         UROTCP(IT0) = OMEGA*(OMY*ZAVE - OMZ*YAVE)
         VROTCP(IT0) = OMEGA*(OMZ*XAVE - OMX*ZAVE)
         WROTCP(IT0) = OMEGA*(OMX*YAVE - OMY*XAVE)
         
 120  CONTINUE
 110  CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE CELLCP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BLADE_POS(OSKU,XCG,YCG,ZCG,
     &        PSIR,THETAR,PHIR,ICYCLE,IPRINT,NGRIFL,T,DT,TRMODE)
C
C ... This subroutine should calclate the position of the blade surface
C ... Preferrably the grid points on a surface should be changed and a
C ... new grid is obtained by deformation in GRIDEF. May require changes
C ... Wave shape is calculöated into 2D arrays and the grid point 
C ... distribution is finally obtained in REDIST. The free surface
C ... algorithm has similar features to the blade vibration.
C ... The parameter list is partially nonsense

      USE TYPE_ARRAYS , ONLY : SIXDOF
      
      IMPLICIT NONE

      INTEGER ICYCLE,IPRINT,IGR,NGRIFL,TRMODE

      TYPE (SIXDOF), DIMENSION(19) :: OSKU

      REAL, DIMENSION(19) :: XCG,YCG,ZCG,
     &        PSIR,THETAR,PHIR
      
      REAL :: T,DT

      RETURN
      END SUBROUTINE BLADE_POS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRIDEF(XNEW,YNEW,ZNEW,XGRIG,YGRIG,ZGRIG,NGPTS,
     &                  XCGI,YCGI,ZCGI,XCG,YCG,ZCG,PSIR,THETAR,
     &                  PHIR,ICYCLE,IPRINT,IREPEA,NREPEA)

C ... This deformates a grid block according to the boundary surface.
C ... Subroutine REDIST does a similar job (can be used as model)

C ... XGRIG, YGRIG, ZGRIG   original grid
C ... XNEW,  YNEW,  ZNEW    manipulated grid
C ... XCGI,  YCGI,  ZCGI    original rotation center (e.g. center of gravity)
C ... XCG,   YCG,   ZCG     new location of the rotation center
C ... PSI,   THETA, PHI     Euler angles (in radians)
C ... Partailly nonsense above, for a model

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XNEW, YNEW, ZNEW 
      REAL, DIMENSION(*) :: XGRIG, YGRIG, ZGRIG 
      REAL :: XCGI, YCGI, ZCGI, XCG, YCG, ZCG
      REAL :: XTMP, YTMP, ZTMP
      REAL :: PSIR, THETAR, PHIR
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XGRIGT, YGRIGT, ZGRIGT 
      REAL :: XCGIT, YCGIT, ZCGIT, XCGT, YCGT, ZCGT

      INTEGER :: I, NGPTS, IREPEA, NREPEA, ICYCLE, IPRINT

      RETURN
      END SUBROUTINE GRIDEF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE NEWPOS(XNEW,YNEW,ZNEW,XGRIG,YGRIG,ZGRIG,NGPTS,
     &                  XCGI,YCGI,ZCGI,XCG,YCG,ZCG,PSIR,THETAR,
     &                  PHIR,ICYCLE,IPRINT,IREPEA,NREPEA,ICASE)

C ... This subroutine rotates and moves a grid block.

C ... XGRIG, YGRIG, ZGRIG   original grid
C ... XNEW,  YNEW,  ZNEW    manipulated grid
C ... XCGI,  YCGI,  ZCGI    original rotation center (e.g. center of gravity)
C ... XCG,   YCG,   ZCG     new location of the rotation center
C ... PSI,   THETA, PHI     Euler angles (in radians)
C ... ICASE                 order of rotations:
C ... 1) PSI, THETA, PHI; 2) PSI, PHI, THETA; 3) THETA, PSI, PHI

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XNEW, YNEW, ZNEW 
      REAL, DIMENSION(*) :: XGRIG, YGRIG, ZGRIG 
      REAL :: XCGI, YCGI, ZCGI, XCG, YCG, ZCG
      REAL :: XTMP, YTMP, ZTMP
      REAL :: PSIR, THETAR, PHIR
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XGRIGT, YGRIGT, ZGRIGT 
      REAL :: XCGIT, YCGIT, ZCGIT, XCGT, YCGT, ZCGT
      REAL :: S1, S2, S3, C1, C2, C3

      INTEGER :: I, NGPTS, IREPEA, NREPEA, ICYCLE, IPRINT, ICASE

      S1 = SIN(PHIR)
      S2 = SIN(THETAR)
      S3 = SIN(PSIR)
      C1 = COS(PHIR)
      C2 = COS(THETAR)
      C3 = COS(PSIR)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A12 =  C2*S3
      A13 = -S2
      A21 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A23 =  C2*S1
      A31 =  C1*C3*S2 + S1*S3
      A32 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A12 =  C3*S1*S2 + C1*S3
      A13 = -C2*S1
      A21 = -C2*S3 
      A22 =  C2*C3
      A23 =  S2
      A31 =  C3*S1 + C1*S2*S3
      A32 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A12 =  S2
      A13 = -C2*S3
      A21 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A23 =  C3*S1 + C1*S2*S3
      A31 =  C3*S1*S2 + C1*S3
      A32 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A12 =  C1*C3*S2 + S1*S3
      A13 =  C1*S2*S3 - C3*S1
      A21 = -S2
      A22 =  C2*C3
      A23 =  C2*S3
      A31 =  C2*S1
      A32 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A12 =  C3*S1 + C1*S2*S3
      A13 =  S1*S3 - C1*C3*S2
      A21 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A23 =  C3*S1*S2 + C1*S3
      A31 =  S2
      A32 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A12 =  C2*S1
      A13 =  C3*S1*S2 - C1*S3
      A21 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A23 =  C1*C3*S2 + S1*S3
      A31 =  C2*S3
      A32 = -S2
      A33 =  C2*C3

      END SELECT


 
      DO I = 1,NGPTS

C ... From FINFLO coordinate system to the 'Euler angle' coordinate system.

         XGRIGT = -XGRIG(I)
         YGRIGT = -ZGRIG(I)
         ZGRIGT = -YGRIG(I)
         XCGIT  = -XCGI
         YCGIT  = -ZCGI
         ZCGIT  = -YCGI
         XCGT   = -XCG
         YCGT   = -ZCG
         ZCGT   = -YCG

C ... Body rotations and movements according to the Euler angles.

         XTMP = A11*(XGRIGT-XCGIT) 
     &        + A21*(YGRIGT-YCGIT) 
     &        + A31*(ZGRIGT-ZCGIT) + XCGT
         YTMP = A12*(XGRIGT-XCGIT) 
     &        + A22*(YGRIGT-YCGIT) 
     &        + A32*(ZGRIGT-ZCGIT) + YCGT
         ZTMP = A13*(XGRIGT-XCGIT) 
     &        + A23*(YGRIGT-YCGIT) 
     &        + A33*(ZGRIGT-ZCGIT) + ZCGT

C ... From 'Euler angle' coordinate system back to the FINFLO 
C ... coordinate system. 

         XNEW(I) = -XTMP
         YNEW(I) = -ZTMP
         ZNEW(I) = -YTMP

      ENDDO

      IREPEA = IREPEA + 1 ! Last sentence for print checking
           
      RETURN
      END SUBROUTINE NEWPOS
C
C ----------------------------------------------------------------------
C ---------------- Edited from NEWPOS by JIl 13.3.2008 -----------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REVPOS(XNEW,YNEW,ZNEW,XGRIG,YGRIG,ZGRIG,NGPTS,
     &                  XCGI,YCGI,ZCGI,XCG,YCG,ZCG,PSIRI,THETARI,
     &                  PHIRI,ICYCLE,IPRINT,IREPEA,NREPEA,ICASE)

C ... This subroutine rotates and moves a grid block, but
C ... the opposite direction compared to NEWPOS

C ... XGRIG, YGRIG, ZGRIG   original grid
C ... XNEW,  YNEW,  ZNEW    manipulated grid
C ... XCGI,  YCGI,  ZCGI    original rotation center (e.g. center of gravity)
C ... XCG,   YCG,   ZCG     new location of the rotation center
C ... PSI,   THETA, PHI     Euler angles (in radians)
C ... ICASE                 order of rotations:
C ... 1) PSI, THETA, PHI; 2) PSI, PHI, THETA; 3) THETA, PSI, PHI

      IMPLICIT NONE

      REAL, DIMENSION(*) :: XNEW, YNEW, ZNEW 
      REAL, DIMENSION(*) :: XGRIG, YGRIG, ZGRIG 
      REAL :: XCGI, YCGI, ZCGI, XCG, YCG, ZCG
      REAL :: XTMP, YTMP, ZTMP
      REAL :: PSIRI, THETARI, PHIRI,PSIR, THETAR, PHIR
      REAL :: A11, A12, A13, A21, A22, A23, A31, A32, A33
      REAL :: XGRIGT, YGRIGT, ZGRIGT 
      REAL :: XCGIT, YCGIT, ZCGIT, XCGT, YCGT, ZCGT
      REAL :: S1, S2, S3, C1, C2, C3

      INTEGER :: I, NGPTS, IREPEA, NREPEA, ICYCLE, IPRINT, ICASE
          
*      PSIR   = PSIRI
*      THETAR = THETARI
*      PHIR   = PHIRI

      S1 = SIN(PHIRI)
      S2 = SIN(THETARI)
      S3 = SIN(PSIRI)
      C1 = COS(PHIRI)
      C2 = COS(THETARI)
      C3 = COS(PSIRI)


      SELECT CASE(ICASE)

      CASE(1)     ! Rotation matrix (xyz)

      A11 =  C2*C3
      A21 =  C2*S3
      A31 = -S2
      A12 =  C3*S1*S2 - C1*S3
      A22 =  C1*C3 + S1*S2*S3
      A32 =  C2*S1
      A13 =  C1*C3*S2 + S1*S3
      A23 =  C1*S2*S3 - C3*S1
      A33 =  C1*C2

      CASE(2)     ! Rotation matrix (yxz)

      A11 =  C1*C3 - S1*S2*S3
      A21 =  C3*S1*S2 + C1*S3
      A31 = -C2*S1
      A12 = -C2*S3 
      A22 =  C2*C3
      A32 =  S2
      A13 =  C3*S1 + C1*S2*S3
      A23 =  S1*S3 - C1*C3*S2
      A33 =  C1*C2

      CASE(3)     ! Rotation matrix (xzy)

      A11 =  C2*C3
      A21 =  S2
      A31 = -C2*S3
      A12 =  S1*S3 - C1*C3*S2
      A22 =  C1*C2
      A32 =  C3*S1 + C1*S2*S3
      A13 =  C3*S1*S2 + C1*S3
      A23 = -C2*S1
      A33 =  C1*C3 - S1*S2*S3

      CASE(4)     ! Rotation matrix (yzx)

      A11 =  C1*C2
      A21 =  C1*C3*S2 + S1*S3
      A31 =  C1*S2*S3 - C3*S1
      A12 = -S2
      A22 =  C2*C3
      A32 =  C2*S3
      A13 =  C2*S1
      A23 =  C3*S1*S2 - C1*S3
      A33 =  C1*C3 + S1*S2*S3

      CASE(5)     ! Rotation matrix (zyx)

      A11 =  C1*C2
      A21 =  C3*S1 + C1*S2*S3
      A31 =  S1*S3 - C1*C3*S2
      A12 = -C2*S1
      A22 =  C1*C3 - S1*S2*S3
      A32 =  C3*S1*S2 + C1*S3
      A13 =  S2
      A23 = -C2*S3
      A33 =  C2*C3

      CASE(6)     ! Rotation matrix (zxy)

      A11 =  C1*C3 + S1*S2*S3
      A21 =  C2*S1
      A31 =  C3*S1*S2 - C1*S3
      A12 =  C1*S2*S3 - C3*S1
      A22 =  C1*C2
      A32 =  C1*C3*S2 + S1*S3
      A13 =  C2*S3
      A23 = -S2
      A33 =  C2*C3

      END SELECT



      DO I = 1,NGPTS

C ... From FINFLO coordinate system to the 'Euler angle' coordinate system.

         XGRIGT = -XGRIG(I)
         YGRIGT = -ZGRIG(I)
         ZGRIGT = -YGRIG(I)
         XCGIT  = -XCGI
         YCGIT  = -ZCGI
         ZCGIT  = -YCGI
         XCGT   = -XCG
         YCGT   = -ZCG
         ZCGT   = -YCG

C ... Body rotations and movements according to the Euler angles.

         XTMP = A11*(XGRIGT-XCGIT) 
     &        + A21*(YGRIGT-YCGIT) 
     &        + A31*(ZGRIGT-ZCGIT) + XCGT
         YTMP = A12*(XGRIGT-XCGIT) 
     &        + A22*(YGRIGT-YCGIT) 
     &        + A32*(ZGRIGT-ZCGIT) + YCGT
         ZTMP = A13*(XGRIGT-XCGIT) 
     &        + A23*(YGRIGT-YCGIT) 
     &        + A33*(ZGRIGT-ZCGIT) + ZCGT

C ... From 'Euler angle' coordinate system back to the FINFLO 
C ... coordinate system. 

         XNEW(I) = -XTMP
         YNEW(I) = -ZTMP
         ZNEW(I) = -YTMP

      ENDDO

      IREPEA = IREPEA + 1 ! Last sentence for print checking
           
      RETURN
      END SUBROUTINE REVPOS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
