C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE COPRPC(DH,DRO,DM,DN,DW,DK,DEPS,RO,E,P,U,V,W,RK,REPS,
     + VOL,PDIFF,DRDP,DRDH,IMAX,JMAX,KMAX,IN,JN,KN,PRO,VAR,NPHASE,
     + MULPHL,DRO1,DFI,FI,KSCAL,MAXSB,TRANSL,TRM,KSTATE,SOLUTION_TYPE)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,I,J,K,IL,II,JJ,
     +     NPHASE,IPH,KSCAL,ISC,MAXSB,KSTATE

      REAL :: DH(*),DRO(*),DM(*),DN(*),DW(*),RO(*),E(*),VOL(*),HTOT,
     +        P(*),U(*),V(*),W(*),DK(*),RK(*),DEPS(*),REPS(*),PDIFF(*),
     +        DRDH(*),DRDP(*),DKIN,apu,DRO1(*),DFI(MAXSB,*),FI(MAXSB,*)

      CHARACTER(*) :: SOLUTION_TYPE

      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(INTERMITTENCY) TRM(*)

      LOGICAL :: MULPHL,TRANSL


      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX
       DM(I+II)  = -(DM(I+II)  -  U(I+II) *DRO(I+II))
       DN(I+II)  = -(DN(I+II)  -  V(I+II) *DRO(I+II))
       DW(I+II)  = -(DW(I+II)  -  W(I+II) *DRO(I+II))

       DK(I+II)  = -(DK(I+II)  -  RK(I+II)*DRO(I+II)/RO(I+II))
       DEPS(I+II)= -(DEPS(I+II)-REPS(I+II)*DRO(I+II)/RO(I+II))

       IF(TRANSL) THEN ! Intermittency variables
          TRM(I+II)%DG   = -(TRM(I+II)%DG   - TRM(I+II)%G   *DRO(I+II))
          TRM(I+II)%DRET = -(TRM(I+II)%DRET - TRM(I+II)%RET *DRO(I+II))
       ENDIF

       DO ISC = 1,KSCAL
         DFI(I+II,ISC)=-(DFI(I+II,ISC)-FI(I+II,ISC)*DRO(I+II)/RO(I+II))
       ENDDO


      IF(SOLUTION_TYPE == 'FLUID') THEN
          HTOT     = (E(I+II) + P(I+II))/RO(I+II)
          DKIN     = DM(I+II)*U(I+II)+DN(I+II)*V(I+II)+DW(I+II)*W(I+II)
     +             + DK(I+II)
          DH(I+II) = -(DH(I+II)-DRO(I+II)*HTOT) +  DKIN
      ELSE IF(SOLUTION_TYPE == 'CAVIT') THEN
          HTOT     = (E(I+II) + P(I+II))/RO(I+II)
          DKIN     = DM(I+II)*U(I+II)+DN(I+II)*V(I+II)+DW(I+II)*W(I+II)
     +             + DK(I+II)
          DH(I+II) = -(DH(I+II)-DRO(I+II)*HTOT) +  DKIN
          DO IPH = 1,NPHASE
           HTOT    = VAR(I+II)%ETOT(IPH) + P(I+II)/PRO(I+II)%RO(IPH)
           DKIN    =(DM(I+II)*U(I+II)+DN(I+II)*V(I+II)+DW(I+II)*W(I+II)
     +             + DK(I+II))*VAR(I+II)%X(IPH)  

           VAR(I+II)%DH(IPH) = -(VAR(I+II)%DH(IPH) -
     +                           VAR(I+II)%DX(IPH) * HTOT) +  DKIN

           VAR(I+II)%DX(IPH) = -(VAR(I+II)%DX(IPH) - VAR(I+II)%X(IPH)
     +                       * DRO(I+II))
          ENDDO ! NPHASE

       ENDIF ! SOLUTION TYPE

      IF(KSTATE /= 10) THEN
       DRO(I+II) = -DRO(I+II) - 1./RO(I+II)*DRDH(I+II)*DH(I+II)
      ELSE
         DRO(I+II) = -DRO(I+II)
      ENDIF
c          DRO1(I+II)= DRO1(I+II) - 
c     +    VAR(I+II)%ALFA(1)/PRO(I+II)%RO(1)*PRO(I+II)%DRODH(1)
c     +    *VAR(I+II)%DH(1) -
c     +    VAR(I+II)%ALFA(2)/PRO(I+II)%RO(2)*PRO(I+II)%DRODH(2)
c     +    *VAR(I+II)%DH(2)
cc      ENDIF

      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE COPRPC
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE IMPSEG(NGL,M,CFM,SMAX,TEMP,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DM,DN,
     + DW,DE,DRK,DEPS,RO,U,V,W,E,C,RK,REPS,CH,CP,PRO,VAR,VOL,VIS,EPS2,
     + VIST,ITURB,JRDIS,JRPRE,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,PR,PRT,
     + XC,YC,ZC,TIMEL,MCYCAM,MGRIDA,ITERAD,ITERAC,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,DIFF,DH,DTEMP,PRC,LUSGS,P,OMEGA,OMEX,OMEY,OMEZ,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,
     + CFL,FRSDEN,FRSVEL,RKMAX,CSIMPS,SRK,SEPS,PTUR,KOVER,
     + PDIFF,DRDP,DRDH,ALFAP,FI,DFI,SFI,MAXSB,NSCAL,MULPHL,BLKS,
     + NPATCH,ICON,IHF,XCP,YCP,ZCP,XCO,YCO,ZCO,FPRINTL,FRSPRE,
     + XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,RKSI,PDFOR,PSIGSC,PSIGS2,
     + TRANSL,TRM,IFLUX,CDIFF,CDIFFT,FREDIF) 

      USE TYPE_ARRAYS
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS, SRNU, EPS10, EPS6, EPS20
      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      INTEGER ICON(IC9,*),IHF(*)

      INTEGER :: M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,I,II,ITERAD,MCYCAM,ITERAC,MGRIDA,ICONV3,IPOST,LUSGS,
     +     IDM1,IDM2,IDM3,NGL,JRDIS,JRPRE,KOVER,ICASE,NS,NSCAL,MAXSB,
     +     MGRID,IPHASE,NPATCH,NPHASE,ICHARP,III,ITERMAX,IEVAP,KSCAL,
     +     KBEGIN,lll,mmm,nnn,IFLUX,KSTATE

      REAL :: TEMP(*),RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),
     +     A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),A2XA(*),A2YA(*),A1ZA(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +     D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),
     +     RESP(*),DIFF(*),DH(*),DTEMP(*),DRO(*),DM(*),DN(*),DW(*),
     +     DRK(*),DEPS(*),U(*),V(*),W(*),E(*),C(*),RK(*),REPS(*),
     +     P(*),UROT(*),VROT(*),WROT(*),SRK(*),SEPS(*),PTUR(*),
     +     FI(MAXSB,MAX(NSCAL,1)),
     ?     DFI(MAXSB,MAX(NSCAL,1)),SFI(MAXSB,MAX(NSCAL,1)),PDIFF(*),
     +     DRDP(*),DRDH(*),XFC(*),YFC(*),ZFC(*),
     +     RKSI(*),PDFOR(*),PSIGSC(*),PSIGS2(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*),
     +        XCP(*),YCP(*),ZCP(*)

      REAL :: DT,CFM,PR,PRT,SMAX,PVISC,CFL,FRSDEN,FRSVEL,CFLL,
     +     CMU,CTA,PSIGK,PSIGE,AA1,ETA0,RKMAX,DTURB,CSIMPS,ALFAP,ALFAU,
     +     ALFAPP,ARTSSP,OMEGA,OMEX,OMEY,OMEZ,PSIGG,FRSPRE,DTM,
     +     APSAT,ALIMIT,APU3,QLIMIT,C1,C2,C21,C3,CDIFF,CDIFFT,CFA,FREDIF

      REAL, ALLOCATABLE :: DP(:),DRO1(:),DTUR(:) ! Mersu, all aux arrays spent
      REAL, ALLOCATABLE :: APP(:),APN(:),APS(:),APE(:),APW(:),APT(:),
     +      APD(:),APU(:),APV(:),DRHOL(:),DRHOG(:),PTUL(:),APUU(:)

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)
      TYPE(INTERMITTENCY) TRM(*)

C ... Solution using a segragated approach

      LOGICAL TIMEL,MULPHL,FPRINTL,MOMCOL,AMGDIVG,TRANSL

      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KSTR    = ISTR*JSTR
      MGRID   = MAX(1,MGRIDA+1-M)
      NPHASE  = BLKS(NGL)%NPHASE

      ALLOCATE (DP(NTOT),DRO1(NTOT),DTUR(NTOT),DRHOL(NTOT),DRHOG(NTOT))   

      ALLOCATE(APP(NTOT),APE(NTOT),APW(NTOT),APN(NTOT),APS(NTOT),
     +         APT(NTOT),APD(NTOT),APU(NTOT),APV(NTOT),APUU(NTOT))      

      ALLOCATE(PTUL(NTOT))     

C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG
c     ICONV3  = 0 ! No convergence monitoring in AMG
    
      IPOST   = 0 ! No post smoothing
      KSCAL   = NSCAL ! No Reynolds stresses in this system
      KBEGIN  = 1
C ... Define the local time step
          
c      CALL TIMEHS(DTL,CH,RO,CP,D1,D2,D3,.01,IMAX,JMAX,KMAX,IN,JN,KN)

      PVISC = 1.8E0*6.**(M-1) ! TESTED 1.5.2008, 6. based on Wigley hull
         
      CFLL   = cfm !  CFL
      ARTSSP = BLKS(NGL)%ARTSSP
      KSTATE = BLKS(NGL)%ISTATE(1)

      IF(.NOT.BLKS(NGL)%COMPCORR) THEN
      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFF)
C ... This does not change much, about CFL=1.5 is possible for Roe
      ELSEIF(BLKS(NGL)%COMPCORR) THEN
      PVISC = 1.8 ! Or does it?
      CALL TIME3(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDM1,IDM2,IDM3,M,TIMEL,
     + PVISC,PTUL,1,FRSPRE,CDIFF)
       ENDIF
       CALL EXTEND(DTL,IMAX,JMAX,KMAX,IN,JN,KN)


C ... To the primitive variables
          
      DO I = 1,NTOT
      DRO(I)   = DRO(I)*VOL(I)
      DM(I)    = DM(I) *VOL(I)
      DN(I)    = DN(I) *VOL(I)
      DW(I)    = DW(I) *VOL(I)
      DH(I)    = DE(I) *VOL(I)
      DTEMP(I) = 0.
      DP(I)    = 0.
      DTUR(I)  = 0.
      ENDDO

c      IF(IPRESC == 1) THEN
      DO I = 1,NTOT
      PRC(I)%DRO   = -PRC(I)%DRO*VOL(I) ! The other DRO is in COPRPC
      ENDDO
c      ENDIF

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO I = 1,NTOT
      DRK(I)   = DRK(I) *VOL(I)
      DEPS(I)  = DEPS(I)*VOL(I)
      ENDDO
      ENDIF

      IF(MULPHL) THEN
         DO I = 1,NTOT
         VAR(I)%DX(1) = VAR(I)%DX(1) * VOL(I)
         VAR(I)%DX(2) = VAR(I)%DX(2) * VOL(I)
         DRHOL(I)     = VAR(I)%DX(1)
         DRHOG(I)     = VAR(I)%DX(2)
         DRO1(I) = -VAR(I)%DX(1)/PRO(I)%RO(1)-VAR(I)%DX(2)/PRO(I)%RO(2)
c     +   - DRO(I)/RO(I) !)*.5  ! April
         VAR(I)%DH(1) = VAR(I)%DE(1) * VOL(I)
         VAR(I)%DH(2) = VAR(I)%DE(2) * VOL(I)
         ENDDO
      ENDIF ! MULPHL
       
      IF(KSCAL > 0) THEN
         DO NS = 1,NSCAL
         DO I  = 1,NTOT
            DFI(I,NS) = DFI(I,NS) * VOL(I)
         ENDDO ; ENDDO
      ENDIF ! KSCAL

      IF(TRANSL) THEN ! Intermittency variables
         DO I = 1,NTOT
            TRM(I)%DG   = TRM(I)%DG   * VOL(I)
            TRM(I)%DRET = TRM(I)%DRET * VOL(I)
         ENDDO
      ENDIF

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE CONSERVATIVE TO PRIMITIVE (COPRPC'
         CALL PRINYS(3,TIC,  DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DX(IPHASE)
         CALL PRINYS(3,CHAR_VAR(13,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DH(IPHASE)
         CALL PRINYS(3,CHAR_VAR(14,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
	   ENDIF ! MULPHL
         CALL PRINYS(3,DMC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DNC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DEC,   DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
C ... FOR SCALAR EQ PPR 24.2
         DO NS = 1,NSCAL
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
         ENDDO
      ENDIF

C-----SOURCES-----------------------------------------------------------

      IF(OMEGA /= 0. .AND. .NOT.TIMEL) THEN
         CALL ROTDIA(DM,DN,DW,DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT)
         IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)'                         '
            WRITE(3,*)' RESULTS AFTER ROTATIONAL CORRECTION'
            CALL PRINYS(3,DMC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DNC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
      ENDIF

C***********************************************************************
C ... From the conservative to the primitive variables
        
      CALL COPRPC(DH,DRO,DM,DN,DW,DRK,DEPS,RO,E,P,U,V,W,RK,
     +  REPS,VOL,PDIFF,DRDP,DRDH,IMAX,JMAX,KMAX,IN,JN,KN,PRO,VAR,
     +  NPHASE,MULPHL,DRO1,DFI(1,KBEGIN),FI(1,KBEGIN),KSCAL,MAXSB,
     +  TRANSL,TRM,KSTATE,BLKS(NGL)%SOLUTION_TYPE)

c      IF (KSCAL >= 1) THEN
c          CALL COPRFI(DFI(1,KBEGIN),FI(1,KBEGIN),RO,DRO,NTOT,KSCAL,
c     +    MAXSB)
c      ENDIF

C ... Print out the mass fluxes of multi-phase flow

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (AFTER COPRPC)'
         WRITE(3,*)'                         '
         CALL PRINYS(3,ROC,  RO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = PRC(1:NTOT)%F1R
         CALL PRINYS(3,F1RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F1R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(5,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRC(1:NTOT)%F2R
         CALL PRINYS(3,F2RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F2R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(6,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRC(1:NTOT)%F3R
         CALL PRINYS(3,F3RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F3R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(7,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC,DRO,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DX(IPHASE)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DH(IPHASE)
         CALL PRINYS(3,CHAR_VAR(14,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRO(1:NTOT)%CP(IPHASE)
         CALL PRINYS(3,CHAR_PH(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
         WRITE(3,*)'                         '
         WRITE(3,*)' IMPSEG: End of multiphase variables'
         WRITE(3,*)'                         '
         ENDIF ! MULPHL

         CALL PRINYS(3,DMC,  DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DNC,  DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,  DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DHC,  DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
      ENDIF

C **********************************************************************
C ... Solve the pressure-velocity coupling
C **********************************************************************

C ... Calculate diffusion coefficient for the momentum equations
          
      CALL DIFFUM(DIFF,VIST,EPS2,VIS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)
      ICASE = BLKS(NGL)%PREVEL

      SELECT CASE(ICASE)

      CASE(1,3) ! Pressures first, a SIMPLEC type solution is possible

      MOMCOL = .true. !.FALSE. ! The effect is not clear, Wigley is more stable?

      CALL PREVEL1(NGL,M,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,
     + A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,DW,DP,DE,RO,
     + U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,VIST,ITURB,JRDIS,JRPRE,
     + IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,ICYCLE,IPRINT,
     + IT,IL,IK,IDI1,IDI2,IDI3,NTOT,XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,
     + ITERAC,NCHIM,AP,AN,AS,AE,AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,
     + KOVER,ALFAP,MULPHL,BLKS,NPATCH,ICON,IHF,XCP,YCP,ZCP,
     + XCO,YCO,ZCO,FPRINTL,FRSPRE,MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,
     + PDIFF,RKSI,PDFOR,APU,TIMEL,DH,IFLUX) 

      CASE(2,4) ! A traditional SIMPLE with velocity correction

      MOMCOL = .false.  ! true Decreases stability?

      CALL PREVEL2(NGL,M,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,
     + A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,DW,DP,DE,RO,
     + U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,VIST,ITURB,JRDIS,JRPRE,
     + IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,ICYCLE,IPRINT,
     + IT,IL,IK,IDI1,IDI2,IDI3,NTOT,XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,
     + ITERAC,NCHIM,AP,AN,AS,AE,AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,
     + KOVER,ALFAP,MULPHL,BLKS,NPATCH,ICON,IHF,XCP,YCP,ZCP,
     + XCO,YCO,ZCO,FPRINTL,FRSPRE,MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,
     + PDIFF,RKSI,PDFOR,APU,TIMEL,DH,IFLUX) 

      END SELECT

      IF(ICYCLE == IPRINT) THEN
         CALL PRINYS(3,DEC,   DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,APC, DIFF, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DUC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DVC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
   
C **********************************************************************
C ... Solve the energy equation(s)
C **********************************************************************

      DO IPHASE = 1,BLKS(NGL)%NPHASE ! Korjaa

      IF(.NOT.MULPHL) THEN

      CALL DIFFUE(DIFF,VIST,CH,CP,PRT,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +     EPS2,VIS)
c         CALL PRINYS(3,DHC,      DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,2,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(956+NGL,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DH,
     +   IMAX,JMAX,KMAX,IN,JN,KN,956+NGL)
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER AAMATF'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
          
      ELSE IF(MULPHL) THEN

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' FOR PHASE :',IPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = PRO(1:NTOT)%CP(IPHASE)
         CALL PRINYS(3,CHAR_PH(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRO(1:NTOT)%CH(IPHASE)
         CALL PRINYS(3,CHAR_PH( 6,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%ETOT(IPHASE) +
     +   P(1:NTOT)/PRO(1:NTOT)%RO(IPHASE)
         CALL PRINYS(3,CHAR_PH(16,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      CALL DIFMFE(DIFF,VIS,EPS2,PRO,VAR,PRT,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

      IEVAP = BLKS(NGL)%EVAP_TYPE

      CALL AAMAMF(AP,AN,AS,AE,AW,AT,AD,PRO,VAR,D1,D2,D3,RO,DRO,VOL,
     +     DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,M,IN,JN,KN,2,IPHASE,
     +     IEVAP,RKSI,NCHIM,ICYCLE,u,v,w,pdiff,DRHOL,DRHOG,APU,TIMEL,DT)

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
        WRITE(960+NGL+(M-1)*1000,*)'               PHASE :',IPHASE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DH,
     +   IMAX,JMAX,KMAX,IN,JN,KN,956+NGL)
      ENDIF

      IF(IPHASE == 1) THEN
         ALIMIT = 10.E-4 ! Edellinen 1.E-4
         QLIMIT = .2E-5

      DO I = 1,NTOT
         DTM    = DTL(I)*.1 + EPS10
         APSAT  = MAX(0.,ALIMIT-VAR(I)%ALFA(IPHASE))*
     +            VOL(I)/DTM * PRO(I)%RO(IPHASE)
         DH(I)  = VAR(I)%DH(IPHASE) 
     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
     +          * PRO(I)%CP(IPHASE)
         AP(I)  = AP(I) + APSAT

c         APSAT  = MAX(0.,QLIMIT-VAR(I)%ALFA(2))*
c     +            VOL(I)/DTM * PRO(I)%RO(2)
c         DH(I)  = VAR(I)%DH(IPHASE) 
c     +          - APSAT*(PRO(I)%DTEMP(2) - PRO(I)%TSAT)
c     +          * PRO(I)%CP(2)


c         VAR(I)%DE(IPHASE) =  VAR(I)%DE(IPHASE) ! kokeile ja merkin valinta
c     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
c     +          * PRO(I)%CP(IPHASE)



      ENDDO

      ELSE IF(IPHASE == 2) THEN
           
         ALIMIT = .2E-5 ! Edellinen .1E-4
         QLIMIT = 100.E-4
      DO I = 1,NTOT
         DTM    = DTL(I)*.1 + EPS10
         APSAT  = MAX(0.,ALIMIT-VAR(I)%ALFA(IPHASE))*
     +            VOL(I)/DTM * PRO(I)%RO(IPHASE)
         DH(I)  = VAR(I)%DH(IPHASE) 
     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
     +          * PRO(I)%CP(IPHASE)
         AP(I)  = AP(I) + APSAT
      ENDDO

      ENDIF

      ENDIF ! .NOT.MULPHL .ELSEIF. MULPHL

C ... Solve using an algebraic multigrid

c      APU(1:NTOT) = VOL(1:NTOT)*DP(1:NTOT)/DTL(1:NTOT)*ALFAP
c      CALL ADJVIN(DH,APU,IMAX,JMAX,KMAX,IN,JN,KN)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE NORMALIZATION'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

C ... This will limit the temeprature change to one degree
      IF(MULPHL) THEN
      AP(1:NTOT) = AP(1:NTOT) + ABS(DH(1:NTOT)/PRO(1:NTOT)%CP(IPHASE))
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT) = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.
  
C ... oli 12
c       MGRID = 1
      CALL SOLAMG(DTEMP,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,5,0,NGL,ICYCLE) ! Oli MGRID

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,DHC,      DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(432+NGL,*)'                         '
         WRITE(432+NGL,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(432+NGL,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(432+NGL,DHC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,
     +        NGL,M)
      ENDIF

      IF(.NOT.MULPHL) THEN
        DTEMP(1:NTOT) = DTEMP(1:NTOT)/CP(1:NTOT)
        DRO(1:NTOT)   = DTEMP(1:NTOT)       ! Back-substitution (mersu)
      ELSE IF(MULPHL) THEN
        DTEMP(1:NTOT) = DTEMP(1:NTOT)/PRO(1:NTOT)%CP(IPHASE)
        VAR(1:NTOT)%DTEMP(IPHASE) = DTEMP(1:NTOT)
      ENDIF

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
      WRITE(432+NGL,*)'                         '
      WRITE(432+NGL,*)' RESULTS AFTER DTEMP/Cp'
      CALL PRINYS(432+NGL,TIC,     DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      PTUL(1:NTOT) =  PRC(1:NTOT)%F2R
      CALL PRINYS(432+NGL,F2RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
      CALL PRINYS(432+NGL,DROC,   DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(432+NGL,DEC,     DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(432+NGL,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      IF(MULPHL) THEN
         PTUL(1:NTOT) =  PRO(1:NTOT)%TSAT
         CALL PRINYS(432+NGL,TEMPC,PTUL,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      ENDIF
      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER DTEMP/Cp'
       IF(.NOT.MULPHL) THEN
         CALL PRINYS(3,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
       ELSE IF(MULPHL) THEN
         PTUL(1:NTOT) =  VAR(1:NTOT)%DTEMP(IPHASE)
         CALL PRINYS(3,CHAR_VAR(4,IPHASE), PTUL,IT,IL,0,IK,
     +   IMAX,JMAX,KMAX,NGL,M)
       ENDIF
      ENDIF ! ICYCLE == IPRINT

      ENDDO ! IPHASES

C **********************************************************************
C ... Solve the mass fractions
C **********************************************************************

      IF(MULPHL) THEN ! Cavitation model

      DIFF(1:NTOT) = 0.

      CFA = 1.*CFLL ! Increase CFL
      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFF)
 
cc      DRO1(1:NTOT) = 0.
c      DO IPHASE = 2,2!!NPHASE
cc      DRO1(1:NTOT) = DRO1(1:NTOT) - (1.-VAR(1:NTOT)%X(IPHASE))*
cc      DRO1(1:NTOT) = DRO1(1:NTOT) -  (1.-VAR(1:NTOT)%ALFA(IPHASE))*
cc     +                               VAR(1:NTOT)%DX(IPHASE)
C ... Yllä potaskaa

      DRO1(1:NTOT) = -VAR(1:NTOT)%DX(2)*VAR(1:NTOT)%X(1)
     +               +VAR(1:NTOT)%DX(1)*VAR(1:NTOT)%X(2)
c      ENDDO

c      CALL EVAP_CORR(PRO,VAR,BLKS,U,V,W,P,PDIFF,VIS,EPS2,FRSDEN,FRSVEL,
c     + IMAX,JMAX,KMAX,IN,JN,KN,DT,DTL,NGL,M,PRT)

      IF(ICYCLE == IPRINT .AND. MULPHL) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE AAMAMA'
         ICHARP = BLKS(NGL)%ICHAR(2)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),DRO1,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      CALL AAMAMA(AP,AN,AS,AE,AW,AT,AD,VAR,PRC,D1,D2,D3,RO,DRO,VOL,DIFF,
     +  A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,4,DRO1,RKSI,NCHIM,
     + BLKS(NGL)%VOIDTT,DT,TIMEL,BLKS(NGL)%SOLUTION_TYPE,NPHASE,NGL,VIS,
     + FREDIF,PRO)

      IF(ICYCLE == IPRINT .AND. MULPHL) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER AAMAMA'
         CALL PRINYS(3,APC, AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(2)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),DRO1,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

       IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(800+NGL,*) 'ICYCLE = ', ICYCLE,' NGL = ',NGL
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,800+NGL)
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DRO1,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,7,0)
     + IN,JN,KN,4,1,5,MGRID,7,0,NGL,ICYCLE)
c     + IN,JN,KN,100,1,5,1,ICONV3,0)

      VAR(1:NTOT)%DX(1) = DTEMP(1:NTOT) !*
c     +                     PRO(1:NTOT)%RO(1)/RO(1:NTOT)
      VAR(1:NTOT)%DX(2) =-DTEMP(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,TIC,  DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DX(IPHASE)                    
         CALL PRINYS(3,APC, AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
      ENDIF

      ENDIF ! MULPHL

C **********************************************************************
C ... Solve the turbulence equations
C **********************************************************************

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Not an algebraic model

      IF(.NOT.BLKS(NGL)%COMPCORR) THEN
      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFFT)
C ... This does not change much, about CFL=1.5 is possible for Roe
      ELSEIF(BLKS(NGL)%COMPCORR) THEN
      PVISC = 1.8 ! Or does it?
      CALL TIME3(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDM1,IDM2,IDM3,M,TIMEL,
     + PVISC,PTUL,1,FRSPRE,CDIFFT)
       ENDIF
       CALL EXTEND(DTL,IMAX,JMAX,KMAX,IN,JN,KN)

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID'
       CALL PRINYS(3,DKC,   DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
       CALL PRINYS(3,DEPSC,DEPS, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      DTUR(1:NTOT) = 0.

      IF(ITURB < 8) THEN ! Solve k-equation for 2-eq. models

      CALL TURBCO(ITURB,JRDIS,JRPRE,C1,C2,C3,C21,CMU,CTA,PSIGK,
     +     PSIGE,AA1,ETA0)

      CALL DIFFUK(DIFF,VIST,VIS,PSIGK,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +     EPS2)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DRK,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ITURB <= 5) THEN
         ICASE = 1 ! k-epsilon model
      ELSE IF(ITURB < 9) THEN
         ICASE = 3 ! k-omega model
      ENDIF

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     +    M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,ICASE,TRM)
C ... Solve using an algebraic multigrid
      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DRK(1:NTOT) = DRK(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DRK,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     +    IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)
      ENDIF ! ITURB < 9

      DRK(1:NTOT) = DTUR(1:NTOT)

      IF(ITURB < 8) THEN      ! Solve K-EPSILON or K-OMEGA
         PSIGG = PSIGE
         CALL DIFFUK(DIFF,VIST,VIS,PSIGG,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +   EPS2)
      ELSE IF(ITURB == 9) THEN   ! Spalart-Allmaras
         PSIGG = SRNU
         CALL DIFFUS(DIFF,REPS,VIS,PSIGG,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)
      ENDIF

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DEPS,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ITURB <= 5) THEN
         ICASE = 2 ! k-epsilon model
      ELSE IF(ITURB < 9) THEN
         ICASE = 4 ! k-omega model
      ELSE         ! Spalart-Allmaras
         ICASE = 6
      ENDIF

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,ICASE,TRM)
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(970+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(970+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DEPS,
     +   IMAX,JMAX,KMAX,IN,JN,KN,970+NGL+(M-1)*1000)
      ENDIF

C ... Solve using an algebraic multigrid
      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DEPS(1:NTOT) = DEPS(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DEPS,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,2,1,5,MGRID,9,0,NGL,ICYCLE)
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)
       DEPS(1:NTOT) = DTUR(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
      CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(3,DEPSC,DEPS, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      ENDIF ! ITURB >= 3      
             
C **********************************************************************
C ... Solve the transition model equations
C **********************************************************************

      IF(TRANSL) THEN ! Intermittency variables

C ... Solve G
      DH(1:NTOT)   = TRM(1:NTOT)%DG

      CALL DIFFTR(DIFF,VIS,EPS2,TRM,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,7,TRM)

      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT)  = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

C ... Solve using an algebraic multigrid
      DTUR(1:NTOT) = 0.
      CALL SOLAMG(DTUR,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     +    IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)

      TRM(1:NTOT)%DG = DTUR(1:NTOT)

C ... Solve RET
      DH(1:NTOT)   = TRM(1:NTOT)%DRET

      CALL DIFFTR(DIFF,VIS,EPS2,TRM,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,2)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,8,TRM)

      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT)  = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

C ... Solve using an algebraic multigrid
      DTUR(1:NTOT) = 0.
      CALL SOLAMG(DTUR,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     +    IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)

      TRM(1:NTOT)%DRET = DTUR(1:NTOT)

      ENDIF ! TRANSL

C **********************************************************************
C ... Solve the scalar equations. Reynolds stresses are not included to scalars
C **********************************************************************

      DO NS = 1,KSCAL

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID'
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
      ENDIF

c      DTUR(1:NTOT) = 0. ! Vector for solution
      DIFF(1:NTOT) = 0.

      CALL DIFFSC(DIFF,VIS,EPS2,PSIGSC,PSIGS2,IMAX,JMAX,KMAX,
     +     NTOT,IN,JN,KN,NS)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DTUR,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURSC(AP,SFI(1,NS),FI(1,NS),VOL,DTL,NTOT,NS,EPS6)

C ... Solve using an algebraic multigrid
      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DTEMP(1:NTOT) = DFI(1:NTOT,NS)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DTEMP,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)

      DFI(1:NTOT,NS) = DTUR(1:NTOT) * RO(1:NTOT) ! Definition of scalars

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      ENDDO ! KSCAL >= 0
             
C **********************************************************************
C ... Finally, store the pressure correction to 'DE'
C **********************************************************************

c      DE(1:NTOT)    = MIN(ABS(ALFAP*DP(1:NTOT)),ABS(0.05*P(1:NTOT)))
c      DE(1:NTOT)    = SIGN(DE(1:NTOT),DP(1:NTOT)) 

      DO I = 1,NTOT
         ALFAPP = ALFAP
         IF(RKSI(I) > 1.) ALFAPP = 1.
         DE(I)  = MIN(ABS(ALFAPP*DP(I)),ABS(0.05*P(I)))
c         DE(I)  = MIN(ABS(ALFAPP*DP(I)),ABS(0.05*FRSPRE)) ! Not robust
         DE(I)  = SIGN(DE(I),DP(I))
         IF(M > 1) DE(I) = 0 ! Occasionally stabilizing
         PRC(I)%DP = DE(I)
      ENDDO

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DRK(1:NTOT)   = DRK(1:NTOT) *RO(1:NTOT)
      DEPS(1:NTOT)  = DEPS(1:NTOT)*RO(1:NTOT)
      ELSE
      DRK(1:NTOT)   = 0.
      DEPS(1:NTOT)  = 0.
      ENDIF

      IF(ITURB == 9) DRK(1:NTOT) = 0.

      RETURN
      END SUBROUTINE IMPSEG
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE PREVEL1(NGL,M,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,
     + DW,DP,DE,RO,U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,
     + VIST,ITURB,JRDIS,JRPRE,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,
     + XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,ITERAC,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,
     + FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,KOVER,
     + ALFAP,MULPHL,BLKS,
     + NPATCH,ICON,IHF,XCP,YCP,ZCP,XCO,YCO,ZCO,FPRINTL,FRSPRE,
     + MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,PDIFF,RKSI,PDFOR,APU,TIMEL,
     + DH,IFLUX) 

      USE TYPE_ARRAYS
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS, SRNU, EPS20
      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IHF(*)

      INTEGER :: M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,I,II,ITERAD,MCYCAM,ITERAC,MGRIDA,ICONV3,IPOST,
     +     IDM1,IDM2,IDM3,NGL,JRDIS,JRPRE,KOVER,ICASE,NS,NSCAL,MAXSB,
     +     MGRID,IPHASE,NPATCH,NPHASE,ICHARP,III,ITERMAX,J,K,JJ,KSTATE,
     +     IFLUX

      REAL RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),DP(*),DRO1(*),
     +     A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),A2XA(*),A2YA(*),A1ZA(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +     D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),
     +     RESP(*),DIFF(*),DTEMP(*),DRO(*),DM(*),DN(*),DW(*),
     +     U(*),V(*),W(*),C(*),P(*),
     +     UROT(*),VROT(*),WROT(*),SRK(*),SEPS(*),PTUR(*),
     +     XFC(*),YFC(*),ZFC(*),PDIFF(*),RKSI(*),
     +     PDFOR(*),APU(*),DH(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*),
     +        XCP(*),YCP(*),ZCP(*)

      REAL :: DT,FRSDEN,FRSVEL,ALFAP,ALFAU,RKMAX,ALFAPP,ARTSSP,
     +        FRSPRE,RLAM,RELAX

      REAL, ALLOCATABLE :: APP(:),APN(:),APS(:),APE(:),APW(:),APT(:),
     +      APD(:),APV(:),DPDX(:),DPDY(:),DPDZ(:),DVP(:),
     +      BETAI(:),BETAJ(:),BETAK(:),BETAIL(:),BETAJL(:),BETAKL(:),
     +      DROB(:),APUU(:)

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)

C ... Pressure correction by updating the pressure field (IPRESC = 1,3)
C ... Perform firstly the pressure solution and then correct momentum resuduals

      LOGICAL :: MOMCOL,MULPHL,FPRINTL,AMGDIVG,TIMEL
         
      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KSTR    = ISTR*JSTR

      MGRID   = MAX(1,MGRIDA+1-M)
      NPHASE  = BLKS(NGL)%NPHASE

      ALLOCATE(APP(NTOT),APE(NTOT),APW(NTOT),APN(NTOT),APS(NTOT),
     +         APT(NTOT),APD(NTOT),DPDX(NTOT),DPDY(NTOT),
     +         DPDZ(NTOT),DVP(NTOT),DROB(NTOT),APUU(NTOT))  ! ,APV(NTOT))      

      DPDX = 0. ; DPDY = 0. ; DPDZ = 0. ; DVP = 0.  

C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG
c     ICONV3  = 0 ! No convergence monitoring in AMG

      IPOST   = 0 ! No post smoothing

C ... Form the matrix for the implicit stage of the momentum equations
          
c      IF(MOMCOL) THEN ! Extend the residuals to ghost cells (no effect)
c         CALL EXTEND(DM,IMAX,JMAX,KMAX,IN,JN,KN)
c         CALL EXTEND(DN,IMAX,JMAX,KMAX,IN,JN,KN)
c         CALL EXTEND(DW,IMAX,JMAX,KMAX,IN,JN,KN)
c      ENDIF ! MOMCOL

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DM,DN,DW,
     +    MOMCOL,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      DIFF(1:NTOT) = AP(1:NTOT) ! Store AP for velocity correction

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (MOMENTUM)'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(950+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(950+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,
     +   IMAX,JMAX,KMAX,IN,JN,KN,950+NGL+(M-1)*1000)
      ENDIF


C ... Solve the Poisson equation for pressure (store AP to DIFF in SIMPLEC)

      ALFAPP = 1. ! Diagonal division?! Does not work

       IF(BLKS(NGL)%PREVEL == 1) THEN ! Do not apply physical DT here.
!           APU(1:NTOT) = AP(1:NTOT) - 1.5*VOL(1:NTOT)*RO(1:NTOT)/DT
       ELSE IF(BLKS(NGL)%PREVEL == 3) THEN
           APU(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/(DTL(1:NTOT)+EPS) ! simplec yritys
cc           APV(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/DTL(1:NTOT)  ! 2. simplec yritys
cc           CALL ADJVIN(APU,APV,IMAX,JMAX,KMAX,IN,JN,KN)      ! 2. simplec yritys
      ENDIF 

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN
     
      CALL AAMATP(APP,APN,APS,APE,APW,APT,APD,APU,PRC,PRO,VAR,RO,DRO,
     + VOL,P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,
     + NGL,XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)
c      CALL AAMATP_DIST(AP,AN,AS,AE,AW,AT,AD,DIFF,PRC,PRO,VAR,RO,DRO,VOL,
c     + P,C,A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,XC,YC,
c     + ZC,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,NGL,D1,D2,D3)

      ELSEIF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

c      APU(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/DTL(1:NTOT) ! simplec yritys
c      CALL AAMATP_CA(APP,APN,APS,APE,APW,APT,APD,APU,PRC,PRO,VAR,RO,DRO,
c     + VOL,P,PDIFF,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
c     + ALFAPP,M,NGL,TIMEL,DT)
      CALL AAMATP_CA(APP,APN,APS,APE,APW,APT,APD,APU,PRC,PRO,VAR,RO,DRO,
     + VOL,P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,
     + NGL,XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)


      ELSE

         STOP'No such flow model'

      ENDIF ! SOLUTION TYPE
         
c      IF(M == 1) CALL AAMATP_WALLS(AP,AN,AS,AE,AW,AT,AD,DRO,DRO1,
       CALL AAMATP_WALLS(APP,APN,APS,APE,APW,APT,APD,DRO,DRO1,NPATCH,
     +    ICON,IHF,XC,YC,ZC,XCO,YCO,ZCO,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     +    A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,M,NGL)

c      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS)
c      AP(1:NTOT) = 1.

c      DRO1(1:NTOT) = 0.
 
      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN
         DRO1(1:NTOT) = DRO(1:NTOT) ! /RO(1:NTOT) ! Muutettu AAMATP
      
c      ELSE
c      DO IPHASE = 1,NPHASE
c         DRO1(1:NTOT) = VAR(1:NTOT)%DX(IPHASE)/PRO(1:NTOT)%RO(IPHASE)
ccc     +   *VAR(1:NTOT)%X(IPHASE)
c     +             + DRO1(1:NTOT)
c      ENDDO
      ENDIF

      IF(NCHIM > 0) THEN ! Fortified algorithm
         DRO1(1:NTOT) = DRO1(1:NTOT) - 
     +                  RKSI(1:NTOT)*(PDFOR(1:NTOT)-PDIFF(1:NTOT))
      ENDIF

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(900+NGL,*) 'ICYCLE = ', ICYCLE,' NGL = ',NGL
         CALL IPRINT2(APP,APN,APS,APE,APW,APT,APD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,900+NGL)
      ENDIF

      APN(1:NTOT) = APN(1:NTOT)/(APP(1:NTOT)+EPS20)
      APS(1:NTOT) = APS(1:NTOT)/(APP(1:NTOT)+EPS20)
      APW(1:NTOT) = APW(1:NTOT)/(APP(1:NTOT)+EPS20)
      APE(1:NTOT) = APE(1:NTOT)/(APP(1:NTOT)+EPS20)
      APT(1:NTOT) = APT(1:NTOT)/(APP(1:NTOT)+EPS20)
      APD(1:NTOT) = APD(1:NTOT)/(APP(1:NTOT)+EPS20)
      DRO1(1:NTOT) = DRO1(1:NTOT)/(APP(1:NTOT)+EPS20)
      APP(1:NTOT) = 1.
      DP(1:NTOT)  = 0.

      IF(MGRID == 1) THEN
         ITERMAX = 200
      ELSE IF(MGRID == 2) THEN
         ITERMAX = 30
      ELSE
         ITERMAX = 20
      ENDIF

C ... Increase AMG cycles with the skewness correction

      IF(M == 1 .AND. BLKS(NGL)%SKEWCORR) THEN

      RELAX = MIN(BLKS(NGL)%SKEWIMIN,BLKS(NGL)%SKEWJMIN,
     +            BLKS(NGL)%SKEWKMIN)
      RELAX = MIN(MAX(RELAX,0.06),BLKS(NGL)%ALFASKEW)

      IF(RELAX < 0.1) THEN
         ITERMAX = MAX(100,ITERMAX)
      ELSE IF(RELAX >= 0.1 .AND. RELAX < 0.2) THEN
         ITERMAX = MAX(60,ITERMAX)
      ELSE IF(RELAX >= 0.2 .AND. RELAX < 0.3) THEN
         ITERMAX = MAX(30,ITERMAX)
      ENDIF

      CALL SOLAMG_FULL(DP,DRO1,APP,APN,APS,APE,APW,APT,APD,IMAX,JMAX,
     + KMAX,IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,ICON,NPATCH,DPDX,DPDY,
     + DPDZ,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,
     + XFC,YFC,ZFC,SDI,BLKS,NGL,AMGDIVG)

      ELSE

      CALL SOLAMG(DP,DRO1,APP,APN,APS,APE,APW,APT,APD,IMAX,JMAX,KMAX,
     + IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,NGL,ICYCLE)

      ENDIF
 
      CALL MIRP(DP,NPATCH,ICON,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER POISSON EQUATION'
         CALL PRINYS(3,DROC,  DRO1, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPC,   DP,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) CALL IPRINT1(DP,IMAX,JMAX,KMAX,IN,JN,KN,532+NGL)
      ENDIF


C ... Solve the momentum equations using an algebraic multigrid

      ALFAU = BLKS(NGL)%ALFAU

C ... Uudestaan (mersu). DIFF contains AP for momentum equations

      CALL GRADP(DPDX,DPDY,DPDZ,DP,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1)


      IF(ICYCLE == IPRINT) THEN
         CALL PRINYS(3,DPDXC, DPDX,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDYC, DPDY,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDZC, DPDZ,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

c      CALL CORRVY(DIFF,DM,DN,DW,DP,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
c     + A2ZA,A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
c     + XFC,YFC,ZFC,DPDX,DPDY,DPDZ,VAR,BLKS(NGL)%SOLUTION_TYPE)

      CALL CORRESMY(DIFF,DM,DN,DW,DP,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,VOL,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + DPDX,DPDY,DPDZ,DH,DVP,1)

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DM(1:NTOT) = DM(1:NTOT)/(AP(1:NTOT)+EPS)
      DN(1:NTOT) = DN(1:NTOT)/(AP(1:NTOT)+EPS)
      DW(1:NTOT) = DW(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.

      CALL SOLAMG(DTEMP,DM,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX, 
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,2,0)
     + IN,JN,KN,2,1,5,MGRID,2,0,NGL,ICYCLE)
      DM(1:NTOT)    = DTEMP(1:NTOT)
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DN,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,3,0)
     + IN,JN,KN,2,1,5,MGRID,3,0,NGL,ICYCLE)
      DN(1:NTOT)    = DTEMP(1:NTOT)
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DW,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,3,0)
     + IN,JN,KN,2,1,5,MGRID,3,0,NGL,ICYCLE)
      DW(1:NTOT)    = DTEMP(1:NTOT)
      DTEMP(1:NTOT) = 0.

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID IN PREVEL1'
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DUC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DVC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      RETURN
      END SUBROUTINE PREVEL1
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE PREVEL2(NGL,M,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,
     + DW,DP,DE,RO,U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,
     + VIST,ITURB,JRDIS,JRPRE,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,
     + XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,ITERAC,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,
     + FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,KOVER,
     + ALFAP,MULPHL,BLKS,
     + NPATCH,ICON,IHF,XCP,YCP,ZCP,XCO,YCO,ZCO,FPRINTL,FRSPRE,
     + MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,PDIFF,RKSI,PDFOR,APU,TIMEL,
     + DH,IFLUX) 

      USE TYPE_ARRAYS
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS, SRNU
      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IHF(*)

      INTEGER :: M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,I,II,ITERAD,MCYCAM,ITERAC,MGRIDA,ICONV3,IPOST,ICER,
     +     IDM1,IDM2,IDM3,NGL,JRDIS,JRPRE,KOVER,ICASE,NS,NSCAL,MAXSB,
     +     MGRID,IPHASE,NPATCH,NPHASE,ICHARP,III,ITERMAX,J,K,JJ,ICERO,
     +     IFLUX

      INTEGER :: KAPUALA

      REAL :: RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),DP(*),DRO1(*),
     +     A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),A2XA(*),A2YA(*),A1ZA(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +     D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),
     +     RESP(*),DIFF(*),DTEMP(*),DRO(*),DM(*),DN(*),DW(*),
     +     U(*),V(*),W(*),C(*),P(*),
     +     UROT(*),VROT(*),WROT(*),SRK(*),SEPS(*),PTUR(*),
     +     XFC(*),YFC(*),ZFC(*),PDIFF(*),RKSI(*),
     +     PDFOR(*),APU(*),DH(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*),
     +        XCP(*),YCP(*),ZCP(*)

      REAL :: DT,FRSDEN,FRSVEL,ALFAP,ALFAU,RKMAX,ALFAPP,ARTSSP,
     +        FRSPRE,RLAM,RELAX,CCC

      REAL, ALLOCATABLE :: APP(:),APN(:),APS(:),APE(:),APW(:),APT(:),
     +      APD(:),APV(:),DPDX(:),DPDY(:),DPDZ(:),DPA(:),DROA(:),
     +      BETAI(:),BETAJ(:),BETAK(:),BETAIL(:),BETAJL(:),BETAKL(:),
     +      DVP(:),DROB(:),APUU(:)

      REAL, ALLOCATABLE :: DPADP(:),DRO1DP(:),APDP(:),ANDP(:),
     +                     ASDP(:),AEDP(:),AWDP(:),ATDP(:),ADDP(:)

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)

C ... Solve the velocities, then pressures and finally correct velocities
C ... A traditional SIMPLE-type approach (IPRESC = 2,4)

      LOGICAL :: MOMCOL,MULPHL,FPRINTL,AMGDIVG,TIMEL
         
      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KSTR    = ISTR*JSTR

      MGRID   = MAX(1,MGRIDA+1-M)
      NPHASE  = BLKS(NGL)%NPHASE

c      ICER    = 0  ! Counter for the internal pressure correction steps

      ALLOCATE(APP(NTOT),DROB(NTOT),APUU(NTOT)) 
!,APE(NTOT),APW(NTOT),APN(NTOT),APS(NTOT),
c     +         APT(NTOT),APD(NTOT),APU(NTOT),APV(NTOT))   

      ALLOCATE(DPDX(NTOT),DPDY(NTOT),DPDZ(NTOT),DROA(NTOT),DPA(NTOT),
     +         DVP(NTOT))  

      DPDX = 0. ; DPDY = 0. ; DPDZ = 0. ; DVP = 0. ; APUU = 0.

c      ALLOCATE(DPADP(NTOT),DRO1DP(NTOT),APDP(NTOT),ANDP(NTOT),
c     +         ASDP(NTOT),AEDP(NTOT),AWDP(NTOT),ATDP(NTOT),ADDP(NTOT))
    

C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG
c     ICONV3  = 0 ! No convergence monitoring in AMG

      IPOST   = 0 ! No post smoothing


C ... Form the matrix for the implicit stage of the momentum equations

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DM,DN,DW,
     +    MOMCOL,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (MOMENTUM)'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(950+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(950+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,
     +   IMAX,JMAX,KMAX,IN,JN,KN,950+NGL+(M-1)*1000)
      ENDIF

C ... Solve the momentum equations using an algebraic multigrid
c       ii = 4-m

c      DIFF(1:NTOT) = APU(1:NTOT) ! Store AP in SIMPLE in AAMATF
            
      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DM(1:NTOT) = DM(1:NTOT)/(AP(1:NTOT)+EPS)
      DN(1:NTOT) = DN(1:NTOT)/(AP(1:NTOT)+EPS)
      DW(1:NTOT) = DW(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.

      CALL SOLAMG(DTEMP,DM,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX, 
     + IN,JN,KN,4,1,5,MGRID,2,0,NGL,ICYCLE)

      DM(1:NTOT) = DTEMP(1:NTOT) ;   DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DN,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,3,0,NGL,ICYCLE)

      DN(1:NTOT)    = DTEMP(1:NTOT) ;   DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DW,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,4,0,NGL,ICYCLE)

      DW(1:NTOT)    = DTEMP(1:NTOT) ;   DTEMP(1:NTOT) = 0.

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DUC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DVC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

C ... Solve the Poisson equation for pressure (store AP to DIFF in SIMPLEC)

c      DIFF(1:NTOT) = 1. !AP(1:NTOT)
c      ALFAU  = 1.
      ALFAPP = 1. ! Diagonal division?! Does not work

       IF(BLKS(NGL)%PREVEL == 2 .AND. TIMEL) THEN      ! SIMPLE
c       DIFF(1:NTOT) = DIFF(1:NTOT) - 1.5*VOL(1:NTOT)*RO(1:NTOT)/DT 
! Use AP from the momentum equation (stored in DIFF). Remove DT-term !
       ELSE IF(BLKS(NGL)%PREVEL == 4) THEN ! SIMPLEC
           APUU(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/(DTL(1:NTOT)+EPS)
c           APUU(1:NTOT) = VOL(1:NTOT)*(PRO(1:NTOT)%RO(1)+
c     2      PRO(1:NTOT)%RO(2))      /(DTL(1:NTOT)+EPS)
c      diff(1:ntot)=diff(1:ntot)-vol(1:ntot)*ro(1:ntot)/(dtl(1:ntot)+eps)
c      diff(1:ntot) = vol(1:ntot)*ro(1:ntot)/(dtl(1:ntot)+eps)
c      diff(1:ntot)= .02 !  diff(1:ntot)   
c           APU(1:NTOT) = VOL(1:NTOT)/(DTL(1:NTOT)+EPS) ! simplec yritys
cc           APV(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/DTL(1:NTOT)  ! 2. simplec yritys
cc           CALL ADJVIN(APU,APV,IMAX,JMAX,KMAX,IN,JN,KN)      ! 2. simplec yritys
      ENDIF 


C ... Build up the pressure correction equation

c      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS)
c      AP(1:NTOT) = 1.

c      DRO1(1:NTOT) = 0.

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN
         DRO1(1:NTOT) = DRO(1:NTOT) ! /RO(1:NTOT) ! Muutettu AAMATP! /RO(1:NTOT) ! Muutettu AAMATP
c      ELSE
c      DO IPHASE = 1,NPHASE
c         DRO1(1:NTOT) = VAR(1:NTOT)%DX(IPHASE)/PRO(1:NTOT)%RO(IPHASE)
ccc     +   *VAR(1:NTOT)%X(IPHASE)
c     +             + DRO1(1:NTOT)
c      ENDDO
      ENDIF ! FLUID

      IF(NCHIM > 0) THEN ! Fortified algorithm
         DRO1(1:NTOT) = DRO1(1:NTOT) - 
     +                  RKSI(1:NTOT)*(PDFOR(1:NTOT)-PDIFF(1:NTOT))
      ENDIF

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(2950+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(2950+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(APUU,AN,AS,AE,AW,AT,AD,DRO,
     +   IMAX,JMAX,KMAX,IN,JN,KN,2950+NGL+(M-1)*1000)
      ENDIF

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN !
      CALL AAMATP(AP,AN,AS,AE,AW,AT,AD,APUU,PRC,PRO,VAR,RO,DRO,VOL,
     + P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,NGL,
     + XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)
c      CALL AAMATP_DIST(AP,AN,AS,AE,AW,AT,AD,APUU,PRC,PRO,VAR,RO,DRO,VOL,
c     + P,C,A1,A2,A3,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,XC,YC,
c     + ZC,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,NGL,D1,D2,D3)

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT' .OR.
     +        BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN ! Temporary

      CALL AAMATP_CA(AP,AN,AS,AE,AW,AT,AD,APUU,PRC,PRO,VAR,RO,DRO,VOL,
     + P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,NGL,
     + XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)

      ELSE
         STOP'No such flow model'
      ENDIF ! SOLUTION TYPE
         
c      IF(M == 1) CALL AAMATP_WALLS(AP,AN,AS,AE,AW,AT,AD,DRO,DRO1,
       CALL AAMATP_WALLS(AP,AN,AS,AE,AW,AT,AD,DRO,DRO1,NPATCH,ICON,IHF,
     +    XC,YC,ZC,XCO,YCO,ZCO,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,
     +    A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,M,NGL)


      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(900+NGL,*) 'ICYCLE = ', ICYCLE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,900+NGL)
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      APP(1:NTOT)= AP(1:NTOT)

C uusi päivityskokeilu
c      CALL CORRM(DRO1,DM,DN,DW,RO,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
c     + A2ZA,A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,
c     + PRC,VOL,NGL)

      DROA(1:NTOT) = DRO1(1:NTOT)
      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS)


      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      KSTR = ISTR*JSTR

c      DO K = -1,KMAX+2
c      JJ   = (KN+K-1)*KSTR
c      DO J = -1,JMAX+2
c      II   = (JN+J-1)*ISTR + JJ + IN
c      DO I = -1,IMAX+2
c       ccc= min(ccc,app(i+ii))
c      DRO1(I+II) = DRO1(I+II)/(AP(I+II)+EPS)
c      ENDDO; ENDDO; ENDDO

      AP(1:NTOT) =MAX(1.,-(AN(1:NTOT)+AS(1:NTOT)+AE(1:NTOT)+AW(1:NTOT)
     +                   +AT(1:NTOT)+AD(1:NTOT)))   ! cc

      IF(MGRID == 1) THEN
         ITERMAX = 200
      ELSE IF(MGRID == 2) THEN
         ITERMAX = 30
      ELSE
         ITERMAX = 20
      ENDIF

      IF(M == 1) THEN
         ICERO = BLKS(NGL)%ICERO
      ELSE
         ICERO = 1
      ENDIF

      ICER    = 0
      ITERMAX = ITERMAX*BLKS(NGL)%ITERMP ! Kokeilu
C ... Increase AMG cycles with the skewness correction

      DP(1:NTOT)  = 0.
      DPA(1:NTOT) = 0.
         
      DO WHILE (ICER < ICERO)

      ICER = ICER + 1

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(1900+NGL*ICER,*) 'ICYCLE = ', ICYCLE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,1900+NGL*ICER)
      ENDIF

      IF(M == 1 .AND. BLKS(NGL)%SKEWCORR) THEN ! Skewness correction

      RELAX = MIN(BLKS(NGL)%SKEWIMIN,BLKS(NGL)%SKEWJMIN,
     +            BLKS(NGL)%SKEWKMIN)
      RELAX = MIN(MAX(RELAX,0.06),BLKS(NGL)%ALFASKEW)

      IF(RELAX < 0.1) THEN
         ITERMAX = MAX(100,ITERMAX)
      ELSE IF(RELAX >= 0.1 .AND. RELAX < 0.2) THEN
         ITERMAX = MAX(60,ITERMAX)
      ELSE IF(RELAX >= 0.2 .AND. RELAX <  0.3) THEN
         ITERMAX = MAX(30,ITERMAX)
      ENDIF

c      DPA(1:NTOT) = 0.
 
      CALL SOLAMG_FULL(DPA,DRO1,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,ICON,NPATCH,DPDX,DPDY,DPDZ,
     + VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,
     + XFC,YFC,ZFC,SDI,BLKS,ngl,AMGDIVG)

      ELSE

c      DPADP(1:NTOT)  = 0.
c      DRO1DP(1:NTOT) = DRO1(1:NTOT)
c      APDP(1:NTOT)   = AP(1:NTOT)
c      ANDP(1:NTOT)   = AN(1:NTOT)
c      ASDP(1:NTOT)   = AS(1:NTOT)
c      AEDP(1:NTOT)   = AE(1:NTOT)
c      AWDP(1:NTOT)   = AW(1:NTOT)
c      ATDP(1:NTOT)   = AT(1:NTOT)
c      ADDP(1:NTOT)   = AD(1:NTOT)

c      ITERMAX = 100 ; ICONV3 = 0
      CALL SOLAMG(DPA,DRO1,AP,AN,AS,AE,AW,AT,AD,
     + IMAX,JMAX,KMAX,IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,NGL,ICYCLE)

      ENDIF

      CALL MIRP(DPA,NPATCH,ICON,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER POISSON EQUATION'
         WRITE(3,*)' DRIVING DRO, PRESSURE CHANGE AT ICERO =',ICER
         CALL PRINYS(3,APC,   APP,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC,  DRO1, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPC,   DPA,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) CALL IPRINT1(DPA,IMAX,JMAX,KMAX,IN,JN,KN,532+NGL)
c        CALL PRINYS(532+NGL,DNC,DN,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
c        CALL PRINYS(532+NGL,DWC,DW,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

C ... Correct the velocities

      ALFAU = BLKS(NGL)%ALFAU

      CALL GRADP(DPDX,DPDY,DPDZ,DPA,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,0)

      IF(ICYCLE == IPRINT) THEN
         CALL PRINYS(3,DPDXC, DPDX,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDYC, DPDY,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDZC, DPDZ,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      IF(M == 1) THEN 
      CALL CORRVY(APUU,DM,DN,DW,DPA,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,DPDX,DPDY,DPDZ,VOL,RKSI,VAR,BLKS(NGL)%SOLUTION_TYPE)
      ELSE
      CALL CORRV(APUU,DM,DN,DW,DPA,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,RKSI,VAR,BLKS(NGL)%SOLUTION_TYPE)
      ENDIF

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' DENS. RES. BEFORE CORRECTION IN PREVEL2'
         CALL PRINYS(3,DROC, DROA, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      IF(ICERO > 1) THEN ! New mass balance
      CALL CORRM2(DROA,DM,DN,DW,RO,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,
     + A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,DPA,APUU,
     + PRC,DROB,VOL,NGL,ICERO)
      ENDIF

      DRO1(1:NTOT) = DROA(1:NTOT)/(APP(1:NTOT)+EPS)
      DP(1:NTOT)   = DP(1:NTOT) + DPA(1:NTOT) !*ALFAP

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' FINAL STATUS AT ICERO =',ICER
         WRITE(3,*)' RESULTS AFTER DENS. CORRECTION IN PREVEL2'
         CALL PRINYS(3,DROC,DROA, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPC,   DP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DUC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DVC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,UC,     U, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,VC,     V, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,WC,     W, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
c      WRITE(3,*)
c      WRITE(3,*) 'ICERO CYCLE =',icero
c      if(m == 1 .and. icero <= 20) go to 1111      

      ENDDO ! ICER < ICERO
        
      RETURN
      END SUBROUTINE PREVEL2
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE PREVELMF(NGL,M,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,
     + DW,DP,DE,RO,UG,VG,WG,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,
     + VIST,ITURB,JRDIS,JRPRE,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,
     + XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,ITERAC,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,
     + FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,KOVER,
     + ALFAP,MULPHL,BLKS,
     + NPATCH,ICON,IHF,XCP,YCP,ZCP,XCO,YCO,ZCO,FPRINTL,FRSPRE,
     + MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,PDIFF,RKSI,PDFOR,APU,TIMEL,
     + DH,IFLUX) 

      USE TYPE_ARRAYS
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS, SRNU, EPS20
      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IHF(*)

      INTEGER :: M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,I,II,ITERAD,MCYCAM,ITERAC,MGRIDA,ICONV3,IPOST,ICER,
     +     IDM1,IDM2,IDM3,NGL,JRDIS,JRPRE,KOVER,ICASE,NS,NSCAL,MAXSB,
     +     MGRID,IPHASE,NPATCH,NPHASE,ICHARP,III,ITERMAX,J,K,JJ,ICERO,
     +     IFLUX,IPH,ICHA

      INTEGER :: KAPUALA

      REAL :: RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),DP(*),DRO1(*),
     +     A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),A2XA(*),A2YA(*),A1ZA(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +     D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),
     +     RESP(*),DIFF(*),DTEMP(*),DRO(*),DM(*),DN(*),DW(*),
     +     UG(*),VG(*),WG(*),C(*),P(*),
     +     UROT(*),VROT(*),WROT(*),SRK(*),SEPS(*),PTUR(*),
     +     XFC(*),YFC(*),ZFC(*),PDIFF(*),RKSI(*),
     +     PDFOR(*),APU(*),DH(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*),
     +        XCP(*),YCP(*),ZCP(*)

      REAL :: DT,FRSDEN,FRSVEL,ALFAP,ALFAU,RKMAX,ALFAPP,ARTSSP,
     +        FRSPRE,RLAM,RELAX,CCC,REFVEL

      REAL, ALLOCATABLE :: APP(:),DVP(:),DROB(:),APUU(:),
     +      DPDX(:),DPDY(:),DPDZ(:),DPA(:),DROA(:),
     +      U(:),V(:),W(:),PTUL(:)

      REAL, ALLOCATABLE :: APL(:),AEL(:),AWL(:),ANL(:),ASL(:),ATL(:),
     +      ADL(:),APG(:),AEG(:),AWG(:),ANG(:),ASG(:),ATG(:),ADG(:),
     +      APUULG(:,:),DMLG(:,:),DNLG(:,:),DWLG(:,:),
     +      DML(:),DNL(:),DWL(:),AP1X(:,:),AP1Y(:,:),AP1Z(:,:),
     +      AP2X(:,:),AP2Y(:,:),AP2Z(:,:),AP3X(:,:),AP3Y(:,:),AP3Z(:,:),
     +      ALF1LG(:,:),ALF2LG(:,:),ALF3LG(:,:),
     +      AP1(:,:),AP2(:,:),AP3(:,:)

      REAL, ALLOCATABLE :: DPADP(:),DRO1DP(:),APDP(:),ANDP(:),
     +                     ASDP(:),AEDP(:),AWDP(:),ATDP(:),ADDP(:)

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)

C ... Solve the velocities, then pressures and finally correct velocities
C ... A traditional SIMPLE-type approach (IPRESC = 2,4)

      LOGICAL :: MOMCOL,MULPHL,FPRINTL,AMGDIVG,TIMEL
         
      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KSTR    = ISTR*JSTR

      MGRID   = MAX(1,MGRIDA+1-M)
      NPHASE  = BLKS(NGL)%NPHASE
      REFVEL  = BLKS(NGL)%REFVEL

c      ICER    = 0  ! Counter for the internal pressure correction steps

      ALLOCATE(APP(NTOT),DROB(NTOT),APUU(NTOT),PTUL(NTOT))

      ALLOCATE(DPDX(NTOT),DPDY(NTOT),DPDZ(NTOT),DROA(NTOT),DPA(NTOT),
     +         DVP(NTOT))  

      ALLOCATE(U(NTOT),V(NTOT),W(NTOT))
      ALLOCATE(DML(NTOT),DNL(NTOT),DWL(NTOT))

      ALLOCATE(APL(NTOT),AEL(NTOT),AWL(NTOT),ANL(NTOT),ASL(NTOT),
     +         ATL(NTOT),ADL(NTOT),APUULG(NTOT,NPHASE)) !,APUUL(NTOT))

      ALLOCATE(APG(NTOT),AEG(NTOT),AWG(NTOT),ANG(NTOT),ASG(NTOT),
     +         ATG(NTOT),ADG(NTOT))!,APUUG(NTOT))

      ALLOCATE(DMLG(NTOT,NPHASE),DNLG(NTOT,NPHASE),DWLG(NTOT,NPHASE),
     +         AP1X(NTOT,NPHASE),AP2X(NTOT,NPHASE),AP3X(NTOT,NPHASE),
     +         AP1Y(NTOT,NPHASE),AP2Y(NTOT,NPHASE),AP3Y(NTOT,NPHASE),
     +         AP1Z(NTOT,NPHASE),AP2Z(NTOT,NPHASE),AP3Z(NTOT,NPHASE),
     +         AP1(NTOT,NPHASE), AP2(NTOT,NPHASE), AP3(NTOT,NPHASE),
     +      ALF1LG(NTOT,NPHASE),ALF2LG(NTOT,NPHASE),ALF3LG(NTOT,NPHASE))

      DPDX = 0. ; DPDY = 0. ; DPDZ = 0. ; DVP = 0. ; APUU = 0.


      IF(ICYCLE == IPRINT) THEN
         DO IPHASE = 1,2 !BLKS(NGL)%NPHASE
         WRITE(3,*)'                         '
         WRITE(3,*)' AT START OF PREVELMF. PHASE =',IPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DMC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)

         IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(11990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'DM, ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(11990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(PTUL,IMAX,JMAX,KMAX,IN,JN,KN,
     +     11990+950+NGL+(M-1)*1000+10*(IPHASE-1)) ! Tässä ok ja kaasu = 0
         WRITE(10990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'DM, ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(10990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(DM,IMAX,JMAX,KMAX,IN,JN,KN,
     +     10990+950+NGL+(M-1)*1000+10*(IPHASE-1)) !
      ENDIF

         ENDDO
       ENDIF

      DO IPHASE = 1,2  !NPHASE

      U(1:NTOT) = VAR(1:NTOT)%U(IPHASE)
      V(1:NTOT) = VAR(1:NTOT)%V(IPHASE)
      W(1:NTOT) = VAR(1:NTOT)%W(IPHASE)

      DML(1:NTOT) = DM(1:NTOT)
      DNL(1:NTOT) = DN(1:NTOT)
      DWL(1:NTOT) = DW(1:NTOT)


C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG
c     ICONV3  = 0 ! No convergence monitoring in AMG

      IPOST   = 0 ! No post smoothing

C ... Form the matrix for the implicit stage of the momentum equations

       CALL DIFMFM(DIFF,VIS,EPS2,PRO,VAR,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

C ... AP diagonal, APU without true-time derivative, APUU without t.d and RKSI
   
      CALL AAMATMF(AP,AN,AS,AE,AW,AT,AD,PRC,PRO,VAR,
     +    D1,D2,D3,RO,DRO,VOL,DIFF,A1,A2,A3,DTL,C,
     +    IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1,IPHASE,FRSDEN,UG,VG,WG, ! Global now
     +    A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DM,DN,DW,
     +    MOMCOL,RKSI,NCHIM,DT,TIMEL,APU,APUU,REFVEL)

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         IPH = 1
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'ICYCLE = ',ICYCLE, 'AAMATF jälkeen: IPHASE = ',IPHASE
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,IMAX,JMAX,KMAX,
     +   IN,JN,KN,2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1))
      ENDIF

C ... Jacobian matrices without the interfacial friction are stored

      IF(IPHASE == 1) THEN
         APL(1:NTOT) = AP(1:NTOT) ;  AEL(1:NTOT)   = AE(1:NTOT)
         AWL(1:NTOT) = AW(1:NTOT) ;  ANL(1:NTOT)   = AN(1:NTOT)
         ASL(1:NTOT) = AS(1:NTOT) ;  ATL(1:NTOT)   = AT(1:NTOT)
         ADL(1:NTOT) = AD(1:NTOT) ;  APUULG(1:NTOT,IPHASE)= APUU(1:NTOT)
      ELSEIF(IPHASE == 2) THEN
         APG(1:NTOT) = AP(1:NTOT) ;  AEG(1:NTOT)   = AE(1:NTOT)
         AWG(1:NTOT) = AW(1:NTOT) ;  ANG(1:NTOT)   = AN(1:NTOT)
         ASG(1:NTOT) = AS(1:NTOT) ;  ATG(1:NTOT)   = AT(1:NTOT)
         ADG(1:NTOT) = AD(1:NTOT) ;  APUULG(1:NTOT,IPHASE)= APUU(1:NTOT)
      ELSE
         STOP ' No more than two phases at the moment. Exiting...'
      ENDIF

      ENDDO ! IPHASE

      DO IPHASE = 1,NPHASE

      CALL IPSAP(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,ADL,
     +    APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DM,NTOT,IMAX,
     +    JMAX,KMAX,NGL,IPHASE,AP1,AP2,AP3,ALF1LG,ALF2LG,ALF3LG,VOL,DTL,
     +    RO)

      ENDDO ! IPHASE

C ... Solve the x-momentum equation using an algebraic multigrid

      DO IPHASE = 1,NPHASE

      CALL IPSA(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,ADL,
     +    APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DM,NTOT,IMAX,
     +    JMAX,KMAX,1,NGL,IPHASE,AP1X,AP1Y,AP1Z,ALF1LG,VOL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         IF(IPHASE == 1) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (UL-MOMENTUM)'
         CALL PRINYS(3,APC,  APL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ELSEIF(IPHASE == 2) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (UG-MOMENTUM)'
         CALL PRINYS(3,APC,  APG, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF

         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         WRITE(3,*)
         WRITE(3,*) ' WITHOUT INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DU(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' WITH INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),DM, IT,IL,0,IK,
     +    IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         IPH = 1
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'ICYCLE = ',ICYCLE, 'IPSAN jälkeen: IPHASE = ',IPHASE
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,IMAX,JMAX,KMAX,
     +   IN,JN,KN,2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1))
         WRITE(11990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'DM, ICYCLE = ',ICYCLE, 'IPSAn jälkeen: IPHASE = ',IPHASE
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         WRITE(11990+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(PTUL,IMAX,JMAX,KMAX,IN,JN,KN,
     +     11990+950+NGL+(M-1)*1000+10*(IPHASE-1))
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS20)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS20)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS20)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS20)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS20)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS20)
      DM(1:NTOT) = DM(1:NTOT)/(AP(1:NTOT)+EPS20) ! DM, DN or DW
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DM,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX, 
     + IN,JN,KN,4,1,5,MGRID,2,0,NGL,ICYCLE)

      DMLG(1:NTOT,IPHASE) = DTEMP(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  DMLG(1:NTOT,IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) THEN
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'AFTER AMG: DM, ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(DTEMP,IMAX,JMAX,KMAX,IN,JN,KN,
     +     1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1))
         ENDIF
      ENDIF

C ... Solve the momentum y-equation using an algebraic multigrid

      CALL IPSA(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,ADL,
     +    APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DN,NTOT,IMAX,
     +    JMAX,KMAX,2,NGL,IPHASE,AP2X,AP2Y,AP2Z,ALF2LG,VOL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         IF(IPHASE == 1) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (VL-MOMENTUM)'
         CALL PRINYS(3,APC,  APL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ELSEIF(IPHASE == 2) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (VG-MOMENTUM)'
         CALL PRINYS(3,APC,  APG, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF

         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         WRITE(3,*)
         WRITE(3,*) ' WITHOUT INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DV(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' WITH INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DN(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),DN, IT,IL,0,IK,
     +    IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         IPH  = 2
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*) 
     +     'ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,IMAX,JMAX,KMAX,
     +IN,JN,KN,2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1))
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS20)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS20)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS20)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS20)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS20)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS20)
      DN(1:NTOT) = DN(1:NTOT)/(AP(1:NTOT)+EPS20) ! DM, DN or DW
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DN,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,3,0,NGL,ICYCLE)

      DNLG(1:NTOT,IPHASE) = DTEMP(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         PTUL(1:NTOT) =  DNLG(1:NTOT,IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) THEN
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'AFTER AMG: DN, ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(DTEMP,IMAX,JMAX,KMAX,IN,JN,KN,
     +     1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1))
         ENDIF
      ENDIF

C ... Solve the z-momentum equation using an algebraic multigrid

      CALL IPSA(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,ADL,
     +    APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DW,NTOT,IMAX,
     +    JMAX,KMAX,3,NGL,IPHASE,AP3X,AP3Y,AP3Z,ALF3LG,VOL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         IF(IPHASE == 1) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (WL-MOMENTUM)'
         CALL PRINYS(3,APC,  APL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ELSEIF(IPHASE == 2) THEN
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (WG-MOMENTUM)'
         CALL PRINYS(3,APC,  APG, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF

         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         WRITE(3,*)
         WRITE(3,*) ' WITHOUT INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DWW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' WITH INTERFACIAL FRICTION'
         PTUL(1:NTOT) =  VAR(1:NTOT)%DW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         WRITE(3,*)
         WRITE(3,*) ' COMBINED AFTER IPSA'
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),DW, IT,IL,0,IK,
     +    IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         IPH = 3
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DM,IMAX,JMAX,KMAX,
     +   IN,JN,KN,2000+1000*IPH+950+NGL+(M-1)*1000+10*(IPHASE-1))
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS20)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS20)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS20)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS20)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS20)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS20)
      DW(1:NTOT) = DW(1:NTOT)/(AP(1:NTOT)+EPS20) ! DM, DN or DW
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DW,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,4,0,NGL,ICYCLE)
       
      DWLG(1:NTOT,IPHASE) = DTEMP(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         PTUL(1:NTOT) =  DWLG(1:NTOT,IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) THEN
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
     +     'AFTER AMG: DW, ICYCLE = ',ICYCLE, 'IPHASE = ',IPHASE
         WRITE(1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1),*)
         CALL IPRINT1(DTEMP,IMAX,JMAX,KMAX,IN,JN,KN,
     +     1970+1000+950+NGL+(M-1)*1000+10*(IPHASE-1))
         ENDIF
      ENDIF

      ENDDO ! IPHASE

      DO IPHASE = 1,NPHASE
      VAR(1:NTOT)%DM(IPHASE) = DMLG(1:NTOT,IPHASE)
      VAR(1:NTOT)%DN(IPHASE) = DNLG(1:NTOT,IPHASE)
      VAR(1:NTOT)%DW(IPHASE) = DWLG(1:NTOT,IPHASE)
      ENDDO

C ... Solve the Poisson equation for pressure (store AP to DIFF in SIMPLEC)

      ALFAPP = 1. ! Obsolate. Diagonal division does not work

       IF(BLKS(NGL)%PREVEL == 2 .AND. TIMEL) THEN      ! SIMPLE
c       DIFF(1:NTOT) = DIFF(1:NTOT) - 1.5*VOL(1:NTOT)*RO(1:NTOT)/DT 
! Use AP from the momentum equation (stored in DIFF). Remove DT-term !
       ELSE IF(BLKS(NGL)%PREVEL == 4) THEN ! SIMPLEC
           DIFF(1:NTOT) = VOL(1:NTOT)*RO(1:NTOT)/(DTL(1:NTOT)+EPS)
      ENDIF 

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN
         DRO1(1:NTOT) = DRO(1:NTOT) ! /RO(1:NTOT) ! Muutettu AAMATP
!                                     /RO(1:NTOT) ! Muutettu AAMATP
      ENDIF ! FLUID

      IF(NCHIM > 0) THEN ! Fortified algorithm
         DRO1(1:NTOT) = DRO1(1:NTOT) - 
     +                  RKSI(1:NTOT)*(PDFOR(1:NTOT)-PDIFF(1:NTOT))
      ENDIF

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(2950+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(2950+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(APUU,AN,AS,AE,AW,AT,AD,DRO,
     +   IMAX,JMAX,KMAX,IN,JN,KN,2950+NGL+(M-1)*1000)
      ENDIF

      IF(BLKS(NGL)%SOLUTION_TYPE == 'FLUID') THEN   ! Obsolate here

      CALL AAMATP(AP,AN,AS,AE,AW,AT,AD,APUU,PRC,PRO,VAR,RO,DRO,VOL,
     + P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAPP,M,NGL,
     + XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)

      ELSE IF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT' .OR.
     +        BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN ! Default

      CALL AAMATP_MF(AP,AN,AS,AE,AW,AT,AD,APUU,PRC,PRO,VAR,RO,DRO,VOL,
     + P,PDIFF,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     + ALFAPP,M,NGL,TIMEL,DT,APUULG,NPHASE,AP1X,AP1Y,AP1Z,AP2X,AP2Y,
     + AP2Z,AP3X,AP3Y,AP3Z,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,
     + ALF1LG,ALF2LG,ALF3LG,AP1,AP2,AP3)

      ELSE
         STOP'No such flow model'
      ENDIF ! SOLUTION TYPE
         
       CALL AAMATP_WALLS(AP,AN,AS,AE,AW,AT,AD,DRO,DRO1,NPATCH,ICON,IHF,
     +    XC,YC,ZC,XCO,YCO,ZCO,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,
     +    A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,M,NGL)

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(900+NGL,*) 'ICYCLE = ', ICYCLE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,900+NGL)
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS20)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS20)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS20)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS20)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS20)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS20)
      APP(1:NTOT)= AP(1:NTOT)

      DROA(1:NTOT) = DRO1(1:NTOT)
      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS20)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      KSTR = ISTR*JSTR

      AP(1:NTOT) =MAX(1.,-(AN(1:NTOT)+AS(1:NTOT)+AE(1:NTOT)+AW(1:NTOT)
     +                   +AT(1:NTOT)+AD(1:NTOT)))   ! cc

      IF(MGRID == 1) THEN
         ITERMAX = 200
      ELSE IF(MGRID == 2) THEN
         ITERMAX = 30
      ELSE
         ITERMAX = 20
      ENDIF

      IF(M == 1) THEN
         ICERO = BLKS(NGL)%ICERO
      ELSE
         ICERO = 1
      ENDIF

      ICER    = 0
      ITERMAX = ITERMAX*BLKS(NGL)%ITERMP ! Kokeilu
C ... Increase AMG cycles with the skewness correction

      DP(1:NTOT)  = 0.
      DPA(1:NTOT) = 0.
         
      DO WHILE (ICER < ICERO)

      ICER = ICER + 1

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(1900+NGL*ICER,*) 'ICYCLE = ', ICYCLE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,1900+NGL*ICER)
      ENDIF

      IF(M == 1 .AND. BLKS(NGL)%SKEWCORR) THEN ! Skewness correction

      RELAX = MIN(BLKS(NGL)%SKEWIMIN,BLKS(NGL)%SKEWJMIN,
     +            BLKS(NGL)%SKEWKMIN)
      RELAX = MIN(MAX(RELAX,0.06),BLKS(NGL)%ALFASKEW)

      IF(RELAX < 0.1) THEN
         ITERMAX = MAX(100,ITERMAX)
      ELSE IF(RELAX >= 0.1 .AND. RELAX < 0.2) THEN
         ITERMAX = MAX(60,ITERMAX)
      ELSE IF(RELAX >= 0.2 .AND. RELAX <  0.3) THEN
         ITERMAX = MAX(30,ITERMAX)
      ENDIF

c      DPA(1:NTOT) = 0.
 
      CALL SOLAMG_FULL(DPA,DRO1,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,ICON,NPATCH,DPDX,DPDY,DPDZ,
     + VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,
     + XFC,YFC,ZFC,SDI,BLKS,ngl,AMGDIVG)

      ELSE

c      DPADP(1:NTOT)  = 0.
c      DRO1DP(1:NTOT) = DRO1(1:NTOT)
c      APDP(1:NTOT)   = AP(1:NTOT)
c      ANDP(1:NTOT)   = AN(1:NTOT)
c      ASDP(1:NTOT)   = AS(1:NTOT)
c      AEDP(1:NTOT)   = AE(1:NTOT)
c      AWDP(1:NTOT)   = AW(1:NTOT)
c      ATDP(1:NTOT)   = AT(1:NTOT)
c      ADDP(1:NTOT)   = AD(1:NTOT)

c      ITERMAX = 100 ; ICONV3 = 0
      CALL SOLAMG(DPA,DRO1,AP,AN,AS,AE,AW,AT,AD,
     + IMAX,JMAX,KMAX,IN,JN,KN,ITERMAX,1,5,MGRID,ICONV3,0,NGL,ICYCLE)
      ENDIF

      CALL MIRP(DPA,NPATCH,ICON,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER POISSON EQUATION'
         WRITE(3,*)' DRIVING DRO, PRESSURE CHANGE AT ICERO =',ICER
         CALL PRINYS(3,APC,   APP,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC,  DRO1, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPC,   DPA,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         DO IPHASE = 1,2 !BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
         IF(FPRINTL) CALL IPRINT1(DPA,IMAX,JMAX,KMAX,IN,JN,KN,532+NGL)
c        CALL PRINYS(532+NGL,DNC,DN,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
c        CALL PRINYS(532+NGL,DWC,DW,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

C ... Correct the velocities

      ALFAU = BLKS(NGL)%ALFAU
   
      CALL GRADP(DPDX,DPDY,DPDZ,DPA,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,0)

      DO IPHASE = 1,NPHASE

c      CALL GRADAP(DPDX,DPDY,DPDZ,DPA,VOL,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
c     + A2ZA,A3,A3XA,A3YA,A3ZA,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
c     + VAR,NPHASE,ALF1LG,ALF2LG,ALF3LG,NGL,IPHASE)

      IF(ICYCLE == IPRINT) THEN
         CALL PRINYS(3,DPDXC, DPDX,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDYC, DPDY,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPDZC, DPDZ,   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      IF(M == 1) THEN 

c      CALL CORRVYMFA(IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,
c     + DPDX,DPDY,DPDZ,VOL,RKSI,VAR,AP1,AP2,AP3,NPHASE,NGL,IPHASE)
      CALL CORRVYMF(IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,
     + DPDX,DPDY,DPDZ,VOL,RKSI,VAR,AP1,AP2,AP3,
     + NPHASE,ALF1LG,ALF2LG,ALF3LG,NGL,IPHASE)

      ELSE ! Does not worka yet
      CALL CORRVMF(APUU,DM,DN,DW,DPA,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,
     + A2ZA,A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,RKSI,VAR,BLKS(NGL)%SOLUTION_TYPE,APUULG)
      ENDIF

      ENDDO

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' DENS. RES. BEFORE CORRECTION IN PREVELMF'
         CALL PRINYS(3,DROC, DROA, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(FPRINTL) CALL IPRINT1(DM,IMAX,JMAX,KMAX,IN,JN,KN,632+NGL)
         IF(FPRINTL) CALL IPRINT1(DN,IMAX,JMAX,KMAX,IN,JN,KN,642+NGL)
         IF(FPRINTL) CALL IPRINT1(DW,IMAX,JMAX,KMAX,IN,JN,KN,648+NGL)
      ENDIF

      IF(ICERO > 1) THEN ! New mass balance
      CALL CORRM2(DROA,DM,DN,DW,RO,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,
     + A3,A3XA,A3YA,A3ZA,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,DPA,APUU,
     + PRC,DROB,VOL,NGL,ICERO)
      ENDIF

      DRO1(1:NTOT) = DROA(1:NTOT)/(APP(1:NTOT)+EPS)
      DP(1:NTOT)   = DP(1:NTOT) + DPA(1:NTOT) !*ALFAP

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' FINAL STATUS AT ICERO =',ICER
         WRITE(3,*)' RESULTS AFTER DENS. CORRECTION IN PREVELMF'
         CALL PRINYS(3,DROC,DROA, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DPC,   DP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DN(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         U(1:NTOT) = VAR(1:NTOT)%U(IPHASE)
         V(1:NTOT) = VAR(1:NTOT)%V(IPHASE)
         W(1:NTOT) = VAR(1:NTOT)%W(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),U,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +   NGL,M)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),V,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +   NGL,M)
         CALL PRINYS(3,CHAR_VAR(35,ICHARP),W,IT,IL,0,IK,IMAX,JMAX,KMAX,
     +   NGL,M)
         ENDDO
      ENDIF

      DM(1:NTOT) = VAR(1:NTOT)%DM(1) ! Antakee mersu
      DN(1:NTOT) = VAR(1:NTOT)%DN(1)
      DW(1:NTOT) = VAR(1:NTOT)%DW(1)

c      VAR(1:NTOT)%DM(2) = VAR(1:NTOT)%DM(1) ! Konvergenssikokeilu
c      VAR(1:NTOT)%DN(2) = VAR(1:NTOT)%DN(1)
c      VAR(1:NTOT)%DW(2) = VAR(1:NTOT)%DW(1)

      ENDDO ! ICER < ICERO
        
      RETURN
      END SUBROUTINE PREVELMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE COPRMF(DH,DRO,DM,DN,DW,DK,DEPS,RO,E,P,U,V,W,RK,REPS,
     + VOL,PDIFF,DRDP,DRDH,IMAX,JMAX,KMAX,IN,JN,KN,PRO,VAR,NPHASE,
     + MULPHL,DRO1,DFI,FI,KSCAL,MAXSB,TRANSL,TRM,KSTATE,SOLUTION_TYPE)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,KSTR,I,J,K,IL,II,JJ,
     +     NPHASE,IPH,KSCAL,ISC,MAXSB,KSTATE

      REAL DH(*),DRO(*),DM(*),DN(*),DW(*),RO(*),E(*),VOL(*),HTOT,
     +     P(*),U(*),V(*),W(*),DK(*),RK(*),DEPS(*),REPS(*),PDIFF(*),
     +     DRDH(*),DRDP(*),DKIN,apu,DRO1(*),DFI(MAXSB,*),FI(MAXSB,*)
      CHARACTER(*) :: SOLUTION_TYPE
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(INTERMITTENCY) TRM(*)
      LOGICAL MULPHL,TRANSL


      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX
       DM(I+II)  = -(DM(I+II)  -  U(I+II) *DRO(I+II))
       DN(I+II)  = -(DN(I+II)  -  V(I+II) *DRO(I+II))
       DW(I+II)  = -(DW(I+II)  -  W(I+II) *DRO(I+II))

       IF(SOLUTION_TYPE == 'MULTI') THEN
          DO IPH = 1,NPHASE
          APU = VAR(I+II)%DM(IPH)
          VAR(I+II)%DM(IPH) = -(VAR(I+II)%DM(IPH) - 
     +                          VAR(I+II)%U(IPH)*VAR(I+II)%DX(IPH))
          VAR(I+II)%DN(IPH) = -(VAR(I+II)%DN(IPH) - 
     +                          VAR(I+II)%V(IPH)*VAR(I+II)%DX(IPH))
          VAR(I+II)%DW(IPH) = -(VAR(I+II)%DW(IPH) - 
     +                          VAR(I+II)%W(IPH)*VAR(I+II)%DX(IPH))
          VAR(I+II)%DU(IPH) = -(VAR(I+II)%DU(IPH) - 
     +                          VAR(I+II)%U(IPH)*VAR(I+II)%DX(IPH))
          VAR(I+II)%DV(IPH) = -(VAR(I+II)%DV(IPH) - 
     +                          VAR(I+II)%V(IPH)*VAR(I+II)%DX(IPH))
          VAR(I+II)%DWW(IPH) = -(VAR(I+II)%DWW(IPH) - 
     +                          VAR(I+II)%W(IPH)*VAR(I+II)%DX(IPH))
          ENDDO
       ENDIF

       DK(I+II)  = -(DK(I+II)  -  RK(I+II)*DRO(I+II)/RO(I+II))
       DEPS(I+II)= -(DEPS(I+II)-REPS(I+II)*DRO(I+II)/RO(I+II))

       IF(TRANSL) THEN ! Intermittency variables
          TRM(I+II)%DG   = -(TRM(I+II)%DG   - TRM(I+II)%G   *DRO(I+II))
          TRM(I+II)%DRET = -(TRM(I+II)%DRET - TRM(I+II)%RET *DRO(I+II))
       ENDIF

       DO ISC = 1,KSCAL
         DFI(I+II,ISC)=-(DFI(I+II,ISC)-FI(I+II,ISC)*DRO(I+II)/RO(I+II))
       ENDDO

         HTOT     = (E(I+II) + P(I+II))/RO(I+II)
         DKIN      = DM(I+II)*U(I+II)+DN(I+II)*V(I+II)+DW(I+II)*W(I+II)
     +           + DK(I+II)
         DH(I+II) = -(DH(I+II)-DRO(I+II)*HTOT) +  DKIN
      DO IPH = 1,NPHASE
          HTOT    = VAR(I+II)%ETOT(IPH) + P(I+II)/PRO(I+II)%RO(IPH)
c          DKIN    =(DM(I+II)*U(I+II)+DN(I+II)*V(I+II)+DW(I+II)*W(I+II)
c     +            + DK(I+II))*VAR(I+II)%X(IPH)  
          DKIN    =(VAR(I+II)%DM(IPH)*VAR(I+II)%U(IPH) +
     +              VAR(I+II)%DN(IPH)*VAR(I+II)%V(IPH)+
     +              VAR(I+II)%DW(IPH)*VAR(I+II)%W(IPH)
     +            + DK(I+II))*VAR(I+II)%X(IPH)  

c         VAR(I+II)%DH(IPH) = -VAR(I+II)%DH(IPH) ! Do not test this
         VAR(I+II)%DH(IPH) = -(VAR(I+II)%DH(IPH) -
     +                         VAR(I+II)%DX(IPH) * HTOT) +  DKIN

          VAR(I+II)%DX(IPH) =-(VAR(I+II)%DX(IPH)-VAR(I+II)%X(IPH)
     +                   * DRO(I+II))
         ENDDO ! NPHASE

      IF(KSTATE /= 10) THEN
         DRO(I+II) = -DRO(I+II) - 1./RO(I+II)*DRDH(I+II)*DH(I+II)
      ELSE
       DRO(I+II) = -DRO(I+II) 
      ENDIF

      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE COPRMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE IMPSMF(NGL,M,CFM,SMAX,TEMP,A1,A1XA,A1YA,A1ZA,A2,
     + A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DM,DN,
     + DW,DE,DRK,DEPS,RO,U,V,W,E,C,RK,REPS,CH,CP,PRO,VAR,VOL,VIS,EPS2,
     + VIST,ITURB,JRDIS,JRPRE,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     + KBOT,KTOP,ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,PR,PRT,
     + XC,YC,ZC,TIMEL,MCYCAM,MGRIDA,ITERAD,ITERAC,NCHIM,AP,AN,AS,AE,
     + AW,AT,AD,RESP,DIFF,DH,DTEMP,PRC,LUSGS,P,OMEGA,OMEX,OMEY,OMEZ,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,
     + CFL,FRSDEN,FRSVEL,RKMAX,CSIMPS,SRK,SEPS,PTUR,KOVER,
     + PDIFF,DRDP,DRDH,ALFAP,FI,DFI,SFI,MAXSB,NSCAL,MULPHL,BLKS,
     + NPATCH,ICON,IHF,XCP,YCP,ZCP,XCO,YCO,ZCO,FPRINTL,FRSPRE,
     + XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,RKSI,PDFOR,PSIGSC,PSIGS2,
     + TRANSL,TRM,IFLUX,CDIFF,CDIFFT,FREDIF) 

      USE TYPE_ARRAYS
      USE CHARACTERS
      USE CONSTANTS, ONLY : EPS, SRNU, EPS10, EPS6, EPS20
      USE NS3CO, ONLY : IN, JN, KN, IC9, TWO_FLUIDL

      IMPLICIT NONE

      INTEGER :: ICON(IC9,*),IHF(*)

      INTEGER :: M,ITURB,IMAX,JMAX,KMAX,IBOT,JBOT,KBOT,ITOP,JTOP,KTOP,
     +     ICYCLE,IPRINT,IT,IL,IK,IDI1,IDI2,IDI3,NTOT,NCHIM,ISTR,JSTR,
     +     KSTR,I,II,ITERAD,MCYCAM,ITERAC,MGRIDA,ICONV3,IPOST,LUSGS,
     +     IDM1,IDM2,IDM3,NGL,JRDIS,JRPRE,KOVER,ICASE,NS,NSCAL,MAXSB,
     +     MGRID,IPHASE,NPATCH,NPHASE,ICHARP,III,ITERMAX,IEVAP,KSCAL,
     +     KBEGIN,lll,mmm,nnn,IFLUX,KSTATE,JJ,KK,J,K,L

      REAL :: TEMP(*),RO(*),CP(*),CH(*),DE(*),VOL(*),DTL(*),
     +     A1(*),A1XA(*),A1YA(*),A2ZA(*),A2(*),A2XA(*),A2YA(*),A1ZA(*),
     +     A3(*),A3XA(*),A3YA(*),A3ZA(*),VIS(*),EPS2(*),VIST(*),
     +     D1(*),D2(*),D3(*),AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),
     +     RESP(*),DIFF(*),DH(*),DTEMP(*),DRO(*),DM(*),DN(*),DW(*),
     +     DRK(*),DEPS(*),U(*),V(*),W(*),E(*),C(*),RK(*),REPS(*),
     +     P(*),UROT(*),VROT(*),WROT(*),SRK(*),SEPS(*),PTUR(*),
     +     FI(MAXSB,MAX(NSCAL,1)),
     ?     DFI(MAXSB,MAX(NSCAL,1)),SFI(MAXSB,MAX(NSCAL,1)),PDIFF(*),
     +     DRDP(*),DRDH(*),XFC(*),YFC(*),ZFC(*),
     +     RKSI(*),PDFOR(*),PSIGSC(*),PSIGS2(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*),
     +        XCP(*),YCP(*),ZCP(*)

      REAL ::DT,CFM,PR,PRT,SMAX,PVISC,CFL,FRSDEN,FRSVEL,CFLL,
     +     CMU,CTA,PSIGK,PSIGE,AA1,ETA0,RKMAX,DTURB,CSIMPS,ALFAP,ALFAU,
     +     ALFAPP,ARTSSP,OMEGA,OMEX,OMEY,OMEZ,PSIGG,FRSPRE,DTM,
     +     APSAT,ALIMIT,APU3,QLIMIT,C1,C2,C21,C3,CDIFF,CDIFFT,ALFA,
     +     FREDIF

      REAL, ALLOCATABLE :: DP(:),DRO1(:),DTUR(:) ! Mersu, all aux arrays spent
      REAL, ALLOCATABLE :: APP(:),APN(:),APS(:),APE(:),APW(:),APT(:),
     +      APD(:),APU(:),APV(:),DRHOL(:),DRHOG(:),PTUL(:),APUU(:)

      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)
      TYPE(INTERMITTENCY) TRM(*)

C ... Solution using a segragated approach

      LOGICAL TIMEL,MULPHL,FPRINTL,MOMCOL,AMGDIVG,TRANSL
         
      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KSTR    = ISTR*JSTR
      MGRID   = MAX(1,MGRIDA+1-M)
      NPHASE  = BLKS(NGL)%NPHASE

      ALLOCATE (DP(NTOT),DRO1(NTOT),DTUR(NTOT),DRHOL(NTOT),DRHOG(NTOT))   

      ALLOCATE(APP(NTOT),APE(NTOT),APW(NTOT),APN(NTOT),APS(NTOT),
     +         APT(NTOT),APD(NTOT),APU(NTOT),APV(NTOT),APUU(NTOT))      

      ALLOCATE(PTUL(NTOT))     
           
C ... These values can be modified here:

c     MCYCLE  = 2 ! Number of multigrid cycles
c     ITERM   = 1 ! LGS-sweeps on the first level
c     ITERH   = 5 ! LGS- sweeps on the coarse levels
c     MGRID       ! Maximum number of multigrid levels is used in AMG
c     ICONV3  = 0 ! No convergence monitoring in AMG

      IPOST   = 0 ! No post smoothing
      KSCAL   = NSCAL ! No Reynolds stresses in this system
      KBEGIN  = 1
       
C ... Define the local time step
          
c      CALL TIMEHS(DTL,CH,RO,CP,D1,D2,D3,.01,IMAX,JMAX,KMAX,IN,JN,KN)

      PVISC = 1.8E0*6.**(M-1) ! TESTED 1.5.2008, 6. based on Wigley hull
 
      CFLL   = cfm !  CFL
      ARTSSP = BLKS(NGL)%ARTSSP
      KSTATE = BLKS(NGL)%ISTATE(1)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      KSTR = ISTR*JSTR

      DO I = 1,NTOT
      U(I) = MAX(VAR(I)%U(1),VAR(I)%U(2))
      V(I) = MAX(VAR(I)%V(1),VAR(I)%V(2))
      W(I) = MAX(VAR(I)%W(1),VAR(I)%W(2))
      ENDDO
       
      IF(.NOT.BLKS(NGL)%COMPCORR) THEN
      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFF)
C ... This does not change much, about CFL=1.5 is possible for Roe
      ELSEIF(BLKS(NGL)%COMPCORR) THEN
      PVISC = 1.8 ! Or does it?
      CALL TIME3(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDM1,IDM2,IDM3,M,TIMEL,
     + PVISC,PTUL,1,FRSPRE,CDIFF)
       ENDIF
       CALL EXTEND(DTL,IMAX,JMAX,KMAX,IN,JN,KN)

      U(1:NTOT) = 0.; V(1:NTOT) = 0.; W(1:NTOT) = 0.

C ... To the primitive variables
         
      DO I = 1,NTOT
      DRO(I)   = DRO(I)*VOL(I)
      DM(I)    = DM(I) *VOL(I)
      DN(I)    = DN(I) *VOL(I)
      DW(I)    = DW(I) *VOL(I)
      DH(I)    = DE(I) *VOL(I)
      DTEMP(I) = 0.
      DP(I)    = 0.
      DTUR(I)  = 0.
      ENDDO

      DO I = 1,NTOT
      PRC(I)%DRO   = -PRC(I)%DRO*VOL(I) ! The other DRO is in COPRPC
      ENDDO

      IF(TWO_FLUIDL) THEN

      DO IPHASE = 1,BLKS(NGL)%NPHASE

c      DO K = 1,KMAX
c      JJ      = (KN+K-1)*KSTR
c      DO J = 1,JMAX
c      II      = (JN+J-1)*ISTR + JJ + IN
c      DO I = 1,IMAX
c         L = I + II
        do l= 1,ntot
c        CALL IJKPAI(l,IMAX,JMAX,KMAX,II,JJ,KK)
         ALFA = VAR(L)%ALFA(IPHASE)
         VAR(L)%DM(IPHASE)  = VAR(L)%DM(IPHASE)  - ALFA*PRC(L)%DPDX
         VAR(L)%DN(IPHASE)  = VAR(L)%DN(IPHASE)  - ALFA*PRC(L)%DPDY
         VAR(L)%DW(IPHASE)  = VAR(L)%DW(IPHASE)  - ALFA*PRC(L)%DPDZ
         VAR(L)%DU(IPHASE)  = VAR(L)%DU(IPHASE)  - ALFA*PRC(L)%DPDX
         VAR(L)%DV(IPHASE)  = VAR(L)%DV(IPHASE)  - ALFA*PRC(L)%DPDY
         VAR(L)%DWW(IPHASE) = VAR(L)%DWW(IPHASE) - ALFA*PRC(L)%DPDZ
c         VAR(L)%DE(IPHASE)  = VAR(L)%DE(IPHASE)  - ALFA*(
c     +                        VAR(L)%U(IPHASE)*PRC(L)%DPDX
c     +                      + VAR(L)%V(IPHASE)*PRC(L)%DPDY
c     +                      + VAR(L)%W(IPHASE)*PRC(L)%DPDZ)
c      write(8000+ngl*10+kk,*)ii,jj,PRC(L)%DPDX,PRC(L)%DPDY,
c     1 PRC(L)%DPDZ
      ENDDO ; ENDDO!; ENDDO; ENDDO

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ADDING PRESSURE GRADIENT'
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
      ENDIF

      ENDIF ! TWO_FLUIDL

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO I = 1,NTOT
      DRK(I)   = DRK(I) *VOL(I)
      DEPS(I)  = DEPS(I)*VOL(I)
      ENDDO
      ENDIF

      IF(MULPHL) THEN
         DO I = 1,NTOT
         VAR(I)%DX(1) = VAR(I)%DX(1) * VOL(I)
         VAR(I)%DX(2) = VAR(I)%DX(2) * VOL(I)
         DRHOL(I)     = VAR(I)%DX(1)
         DRHOG(I)     = VAR(I)%DX(2)
         DRO1(I) = -VAR(I)%DX(1)/PRO(I)%RO(1)- VAR(I)%DX(2)/PRO(I)%RO(2)
c     +   - DRO(I)/RO(I) !)*.5  ! April. Yllä tärkeää
         VAR(I)%DH(1) = VAR(I)%DE(1) * VOL(I)
         VAR(I)%DH(2) = VAR(I)%DE(2) * VOL(I)
         VAR(I)%DM(1) = VAR(I)%DM(1) * VOL(I)
         VAR(I)%DM(2) = VAR(I)%DM(2) * VOL(I)
         VAR(I)%DN(1) = VAR(I)%DN(1) * VOL(I)
         VAR(I)%DN(2) = VAR(I)%DN(2) * VOL(I)
         VAR(I)%DW(1) = VAR(I)%DW(1) * VOL(I)
         VAR(I)%DW(2) = VAR(I)%DW(2) * VOL(I)
         VAR(I)%DU(1) = VAR(I)%DU(1) * VOL(I)
         VAR(I)%DU(2) = VAR(I)%DU(2) * VOL(I)
         VAR(I)%DV(1) = VAR(I)%DV(1) * VOL(I)
         VAR(I)%DV(2) = VAR(I)%DV(2) * VOL(I)
         VAR(I)%DWW(1) = VAR(I)%DWW(1) * VOL(I)
         VAR(I)%DWW(2) = VAR(I)%DWW(2) * VOL(I)
         ENDDO
      ENDIF ! MULPHL
       
      IF(KSCAL > 0) THEN
         DO NS = 1,NSCAL
         DO I  = 1,NTOT
            DFI(I,NS) = DFI(I,NS) * VOL(I)
         ENDDO ; ENDDO
      ENDIF ! KSCAL

      IF(TRANSL) THEN ! Intermittency variables
         DO I = 1,NTOT
            TRM(I)%DG   = TRM(I)%DG   * VOL(I)
            TRM(I)%DRET = TRM(I)%DRET * VOL(I)
         ENDDO
      ENDIF

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE CONSERVATIVE TO PRIMITIVE (COPRMF)'
         CALL PRINYS(3,TIC,  DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC, DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DX(IPHASE)
         CALL PRINYS(3,CHAR_VAR(13,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DH(IPHASE)
         CALL PRINYS(3,CHAR_VAR(14,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
	   ENDIF ! MULPHL
         DO IPHASE = 1,2 !BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%V(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%W(IPHASE)

         CALL PRINYS(3,CHAR_VAR(35,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO

         DO IPHASE = 1,1 !BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),U(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(34,ICHARP),V(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,CHAR_VAR(35,ICHARP),W(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO

         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE CONSERVATIVE TO PRIMITIVE (COPRMF)'
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DMC,   DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DN(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DNC,   DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,   DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
         CALL PRINYS(3,DEC,   DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
C ... FOR SCALAR EQ PPR 24.2
         DO NS = 1,NSCAL
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
         ENDDO
      ENDIF

C-----SOURCES-----------------------------------------------------------

      IF(OMEGA /= 0. .AND. .NOT.TIMEL) THEN
         CALL ROTDIAMF(VAR,DTL,OMEGA,OMEX,OMEY,OMEZ,NTOT,NPHASE)
         IF(ICYCLE == IPRINT) THEN
            WRITE(3,*)'                         '
            WRITE(3,*)' RESULTS AFTER ROTATIONAL CORRECTION'
            CALL PRINYS(3,DMC,  DM,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DNC,  DN,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DWC,  DW,  IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
      ENDIF

C***********************************************************************
C ... From the conservative to the primitive variables

      CALL COPRMF(DH,DRO,DM,DN,DW,DRK,DEPS,RO,E,P,U,V,W,RK,
     +  REPS,VOL,PDIFF,DRDP,DRDH,IMAX,JMAX,KMAX,IN,JN,KN,PRO,VAR,
     +  NPHASE,MULPHL,DRO1,DFI(1,KBEGIN),FI(1,KBEGIN),KSCAL,MAXSB,
     +  TRANSL,TRM,KSTATE,BLKS(NGL)%SOLUTION_TYPE)

c      IF (KSCAL >= 1) THEN
c          CALL COPRFI(DFI(1,KBEGIN),FI(1,KBEGIN),RO,DRO,NTOT,KSCAL,
c     +    MAXSB)
c      ENDIF

C ... Print out the mass fluxes of multi-phase flow

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID (AFTER COPRMF)'
         WRITE(3,*)'                         '
         CALL PRINYS(3,ROC,  RO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF(MULPHL) THEN
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = PRC(1:NTOT)%F1R
         CALL PRINYS(3,F1RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F1R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(5,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRC(1:NTOT)%F2R
         CALL PRINYS(3,F2RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F2R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(6,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRC(1:NTOT)%F3R
         CALL PRINYS(3,F3RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%F3R(IPHASE)
         CALL PRINYS(3,CHAR_VAR(7,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DROC,DRO,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DX(IPHASE)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DH(IPHASE)
         CALL PRINYS(3,CHAR_VAR(14,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRO(1:NTOT)%CP(IPHASE)
         CALL PRINYS(3,CHAR_PH(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO

        DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DMC,  DM, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DN(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DNC,  DN, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL(1),
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DWC,  DW, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
         WRITE(3,*)'                         '
         WRITE(3,*)' IMPSMF: End of multiphase variables'
         WRITE(3,*)'                         '
         ENDIF ! MULPHL
         CALL PRINYS(3,DHC,  DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         IF (ITURB >= 3 .AND. ITURB /= 8) THEN
            CALL PRINYS(3,DKC,  DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
            CALL PRINYS(3,DEPSC,DEPS,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDIF
      ENDIF

C **********************************************************************
C ... Solve the pressure-velocity coupling
C **********************************************************************

C ... Calculate diffusion coefficient for the momentum equations

c      CALL DIFFUM(DIFF,VIST,EPS2,VIS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)
      ICASE = BLKS(NGL)%PREVEL
 
      SELECT CASE(ICASE)

      CASE(1,3) ! Pressures first, a SIMPLEC type solution is possible

      STOP ' The pre-vel coupling option does not exist, use 2'

      MOMCOL = .true. !.FALSE. ! The effect is not clear
       
      CALL PREVEL1(NGL,M,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,
     + A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,DW,DP,DE,RO,
     + U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,VIST,ITURB,JRDIS,JRPRE,
     + IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,ICYCLE,IPRINT,
     + IT,IL,IK,IDI1,IDI2,IDI3,NTOT,XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,
     + ITERAC,NCHIM,AP,AN,AS,AE,AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,
     + KOVER,ALFAP,MULPHL,BLKS,NPATCH,ICON,IHF,XCP,YCP,ZCP,
     + XCO,YCO,ZCO,FPRINTL,FRSPRE,MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,
     + PDIFF,RKSI,PDFOR,APU,TIMEL,DH,IFLUX) 

      CASE(2,4) ! A traditional SIMPLE with velocity correction

      MOMCOL = .false.  ! true Decreases stability?

      CALL PREVELMF(NGL,M,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,
     + A3XA,A3YA,A3ZA,D1,D2,D3,DTL,DT,DRO,DRO1,DM,DN,DW,DP,DE,RO,
     + U,V,W,P,C,CH,CP,PRO,VAR,VOL,VIS,EPS2,VIST,ITURB,JRDIS,JRPRE,
     + IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,ICYCLE,IPRINT,
     + IT,IL,IK,IDI1,IDI2,IDI3,NTOT,XC,YC,ZC,MCYCAM,MGRIDA,ITERAD,
     + ITERAC,NCHIM,AP,AN,AS,AE,AW,AT,AD,RESP,DIFF,DTEMP,PRC,
     + UROT,VROT,WROT,IDM1,IDM2,IDM3,FRSDEN,FRSVEL,RKMAX,SRK,SEPS,PTUR,
     + KOVER,ALFAP,MULPHL,BLKS,NPATCH,ICON,IHF,XCP,YCP,ZCP,
     + XCO,YCO,ZCO,FPRINTL,FRSPRE,MOMCOL,XFC,YFC,ZFC,SDI,AMGDIVG,ICONV3,
     + PDIFF,RKSI,PDFOR,APU,TIMEL,DH,IFLUX) 

      END SELECT

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*) ' After PREVELMF'
         CALL PRINYS(3,DEC,   DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,APC, DIFF, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         DO IPHASE = 1,BLKS(NGL)%NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = VAR(1:NTOT)%DM(IPHASE)
         CALL PRINYS(3,CHAR_VAR(30,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DN(IPHASE)
         CALL PRINYS(3,CHAR_VAR(31,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%DW(IPHASE)
         CALL PRINYS(3,CHAR_VAR(32,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) =  VAR(1:NTOT)%U(IPHASE)
         CALL PRINYS(3,CHAR_VAR(33,ICHARP),PTUL,
     +    IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO

      ENDIF

C **********************************************************************
C ... Solve the energy equation(s)
C **********************************************************************

      DO IPHASE = 1,BLKS(NGL)%NPHASE ! Korjaa

      IF(.NOT.MULPHL) THEN

      CALL DIFFUE(DIFF,VIST,CH,CP,PRT,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +     EPS2,VIS)
c         CALL PRINYS(3,DHC,      DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,2,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(956+NGL,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DH,
     +   IMAX,JMAX,KMAX,IN,JN,KN,956+NGL)
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER AAMATF'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
          
      ELSE IF(MULPHL) THEN

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' FOR PHASE :',IPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         PTUL(1:NTOT) = PRO(1:NTOT)%CP(IPHASE)
         CALL PRINYS(3,CHAR_PH(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = PRO(1:NTOT)%CH(IPHASE)
         CALL PRINYS(3,CHAR_PH( 6,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) = VAR(1:NTOT)%ETOT(IPHASE) +
     +   P(1:NTOT)/PRO(1:NTOT)%RO(IPHASE)
         CALL PRINYS(3,CHAR_PH(16,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
       
      CALL DIFMFE(DIFF,VIS,EPS2,PRO,VAR,PRT,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

      IEVAP = BLKS(NGL)%EVAP_TYPE

      CALL AAMAMF(AP,AN,AS,AE,AW,AT,AD,PRO,VAR,D1,D2,D3,RO,DRO,VOL,
     +     DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,M,IN,JN,KN,2,IPHASE,
     +     IEVAP,RKSI,NCHIM,ICYCLE,u,v,w,pdiff,DRHOL,DRHOG,APU,TIMEL,DT)
        
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
        WRITE(960+NGL+(M-1)*1000,*)'               PHASE :',IPHASE
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DH,
     +   IMAX,JMAX,KMAX,IN,JN,KN,956+NGL)
      ENDIF

      IF(IPHASE == 1) THEN
         ALIMIT = 10.E-4 ! Edellinen 1.E-4
         QLIMIT = .2E-5

      DO I = 1,NTOT
         DTM    = DTL(I)*.1 + EPS10
         APSAT  = MAX(0.,ALIMIT-VAR(I)%ALFA(IPHASE))*
     +            VOL(I)/DTM * PRO(I)%RO(IPHASE)
         DH(I)  = VAR(I)%DH(IPHASE) 
     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
     +          * PRO(I)%CP(IPHASE)
         AP(I)  = AP(I) + APSAT

c         APSAT  = MAX(0.,QLIMIT-VAR(I)%ALFA(2))*
c     +            VOL(I)/DTM * PRO(I)%RO(2)
c         DH(I)  = VAR(I)%DH(IPHASE) 
c     +          - APSAT*(PRO(I)%DTEMP(2) - PRO(I)%TSAT)
c     +          * PRO(I)%CP(2)


c         VAR(I)%DE(IPHASE) =  VAR(I)%DE(IPHASE) ! kokeile ja merkin valinta
c     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
c     +          * PRO(I)%CP(IPHASE)



      ENDDO

      ELSE IF(IPHASE == 2) THEN
           
         ALIMIT = .2E-5 ! Edellinen .1E-4
         QLIMIT = 100.E-4
      DO I = 1,NTOT
         DTM    = DTL(I)*.1 + EPS10
         APSAT  = MAX(0.,ALIMIT-VAR(I)%ALFA(IPHASE))*
     +            VOL(I)/DTM * PRO(I)%RO(IPHASE)
         DH(I)  = VAR(I)%DH(IPHASE) 
     +          + APSAT*(PRO(I)%DTEMP(IPHASE) - PRO(I)%TSAT)
     +          * PRO(I)%CP(IPHASE)
         AP(I)  = AP(I) + APSAT
      ENDDO

      ENDIF

      ENDIF ! .NOT.MULPHL .ELSEIF. MULPHL

C ... Solve using an algebraic multigrid

c      APU(1:NTOT) = VOL(1:NTOT)*DP(1:NTOT)/DTL(1:NTOT)*ALFAP
c      CALL ADJVIN(DH,APU,IMAX,JMAX,KMAX,IN,JN,KN)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE NORMALIZATION'
         CALL PRINYS(3,APC,   AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

C ... This will limit the temeprature change to one degree
      IF(MULPHL) THEN
      AP(1:NTOT) = AP(1:NTOT) + ABS(DH(1:NTOT)/PRO(1:NTOT)%CP(IPHASE))
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT) = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

C ... oli 12
c       MGRID = 1
      CALL SOLAMG(DTEMP,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,4,1,5,MGRID,5,0,NGL,ICYCLE) ! Oli MGRID

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(3,DHC,      DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(432+NGL,*)'                         '
         WRITE(432+NGL,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         CALL PRINYS(432+NGL,DHC,   DH, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(432+NGL,DHC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,
     +        NGL,M)
      ENDIF

      IF(.NOT.MULPHL) THEN
        DTEMP(1:NTOT) = DTEMP(1:NTOT)/CP(1:NTOT)
        DRO(1:NTOT)   = DTEMP(1:NTOT)       ! Back-substitution (mersu)
      ELSE IF(MULPHL) THEN
        DTEMP(1:NTOT) = DTEMP(1:NTOT)/PRO(1:NTOT)%CP(IPHASE)
        VAR(1:NTOT)%DTEMP(IPHASE) = DTEMP(1:NTOT)
      ENDIF
 
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
      WRITE(432+NGL,*)'                         '
      WRITE(432+NGL,*)' RESULTS AFTER DTEMP/Cp'
      CALL PRINYS(432+NGL,TIC,     DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      PTUL(1:NTOT) =  PRC(1:NTOT)%F2R
      CALL PRINYS(432+NGL,F2RC,PTUL, IT,IL,0,IK,IMAX,JMAX,
     +   KMAX,NGL,M)
      CALL PRINYS(432+NGL,DROC,   DRO, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(432+NGL,DEC,     DE, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(432+NGL,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      IF(MULPHL) THEN
         PTUL(1:NTOT) =  PRO(1:NTOT)%TSAT
         CALL PRINYS(432+NGL,TEMPC,PTUL,IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF
      ENDIF
      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER DTEMP/Cp'
       IF(.NOT.MULPHL) THEN
         CALL PRINYS(3,DTEMPC,DTEMP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
       ELSE IF(MULPHL) THEN
         PTUL(1:NTOT) =  VAR(1:NTOT)%DTEMP(IPHASE)
         CALL PRINYS(3,CHAR_VAR(4,IPHASE), PTUL,IT,IL,0,IK,
     +   IMAX,JMAX,KMAX,NGL,M)
       ENDIF
      ENDIF ! ICYCLE == IPRINT

      ENDDO ! IPHASES
       
C **********************************************************************
C ... Solve the mass fractions
C **********************************************************************

      IF(MULPHL) THEN ! Obsolate announcement: Cavitation model

      DIFF(1:NTOT) = 0.

      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFF)

cc      DRO1(1:NTOT) = 0.
c      DO IPHASE = 2,2!!NPHASE
cc      DRO1(1:NTOT) = DRO1(1:NTOT) - (1.-VAR(1:NTOT)%X(IPHASE))*
cc      DRO1(1:NTOT) = DRO1(1:NTOT) -  (1.-VAR(1:NTOT)%ALFA(IPHASE))*
cc     +                               VAR(1:NTOT)%DX(IPHASE)
C ... Yllä potaskaa

      DRO1(1:NTOT) = -VAR(1:NTOT)%DX(2)*VAR(1:NTOT)%X(1)
     +               +VAR(1:NTOT)%DX(1)*VAR(1:NTOT)%X(2)
c     +              -VAR(1:NTOT)%DX(2) ! Ei muuta mitään
c      ENDDO

c      CALL EVAP_CORR(PRO,VAR,BLKS,U,V,W,P,PDIFF,VIS,EPS2,FRSDEN,FRSVEL,
c     + IMAX,JMAX,KMAX,IN,JN,KN,DT,DTL,NGL,M,PRT)

      CALL AAMAMA(AP,AN,AS,AE,AW,AT,AD,VAR,PRC,D1,D2,D3,RO,DRO,VOL,DIFF,
     +  A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,4,DRO1,RKSI,NCHIM,
     + BLKS(NGL)%VOIDTT,DT,TIMEL,BLKS(NGL)%SOLUTION_TYPE,NPHASE,NGL,VIS,
     + FREDIF,PRO)

       IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(800+NGL,*) 'ICYCLE = ', ICYCLE,' NGL = ',NGL
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO1,
     +   IMAX,JMAX,KMAX,IN,JN,KN,800+NGL)
      ENDIF

      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DRO1(1:NTOT) = DRO1(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.
      DTEMP(1:NTOT) = 0.

      CALL SOLAMG(DTEMP,DRO1,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,7,0)
     + IN,JN,KN,4,1,5,MGRID,7,0,NGL,ICYCLE)
c     + IN,JN,KN,100,1,5,1,ICONV3,0)

      VAR(1:NTOT)%DX(1) = DTEMP(1:NTOT) !*
c     +                     PRO(1:NTOT)%RO(1)/RO(1:NTOT)
      VAR(1:NTOT)%DX(2) =-DTEMP(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
         DO IPHASE = 1,NPHASE
         ICHARP = BLKS(NGL)%ICHAR(IPHASE)
         CALL PRINYS(3,TIC,  DTL, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         PTUL(1:NTOT) =  VAR(1:NTOT)%DX(IPHASE)                    
         CALL PRINYS(3,APC, AP, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         CALL PRINYS(3,CHAR_VAR(11,ICHARP),PTUL,
     +   IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
         ENDDO
      ENDIF

      ENDIF ! MULPHL

C **********************************************************************
C ... Solve the turbulence equations
C **********************************************************************

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Not an algebraic model

      IF(.NOT.BLKS(NGL)%COMPCORR) THEN
      CALL TIMEIN(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFLL,VIS,VIST,
     + UROT,VROT,WROT,FRSDEN,FRSVEL,ARTSSP,PR,PRT,IDM1,IDM2,IDM3,M,
     + TIMEL,PVISC,DTUR,1,CDIFFT)
C ... This does not change much, about CFL=1.5 is possible for Roe
      ELSEIF(BLKS(NGL)%COMPCORR) THEN
      PVISC = 1.8 ! Or does it?
      CALL TIME3(A1,A2,A3,D1,D2,D3,VOL,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,
     + A3XA,A3YA,A3ZA,DTL,DT,RO,U,V,W,C,P,IMAX,JMAX,KMAX,CFL,VIS,
     + VIST,UROT,VROT,WROT,FRSDEN,FRSVEL,PR,PRT,IDM1,IDM2,IDM3,M,TIMEL,
     + PVISC,PTUL,1,FRSPRE,CDIFFT)
       ENDIF
       CALL EXTEND(DTL,IMAX,JMAX,KMAX,IN,JN,KN)

      IF(ICYCLE == IPRINT) THEN
         WRITE(3,*)'                         '
         WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID'
       CALL PRINYS(3,DKC,   DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
       CALL PRINYS(3,DEPSC,DEPS, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      DTUR(1:NTOT) = 0.

      IF(ITURB < 8) THEN ! Solve k-equation for 2-eq. models

      CALL TURBCO(ITURB,JRDIS,JRPRE,C1,C2,C3,C21,CMU,CTA,PSIGK,
     +     PSIGE,AA1,ETA0)

      CALL DIFFUK(DIFF,VIST,VIS,PSIGK,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +     EPS2)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DRK,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ITURB <= 5) THEN
         ICASE = 1 ! k-epsilon model
      ELSE IF(ITURB < 9) THEN
         ICASE = 3 ! k-omega model
      ENDIF

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,ICASE,TRM)
C ... Solve using an algebraic multigrid
      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DRK(1:NTOT) = DRK(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DRK,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     +    IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)
      ENDIF ! ITURB < 9

      DRK(1:NTOT) = DTUR(1:NTOT)

      IF(ITURB < 8) THEN      ! Solve K-EPSILON or K-OMEGA
         PSIGG = PSIGE
         CALL DIFFUK(DIFF,VIST,VIS,PSIGG,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +    EPS2)
      ELSE IF(ITURB == 9) THEN   ! Spalart-Allmaras
         PSIGG = SRNU
         CALL DIFFUS(DIFF,REPS,VIS,PSIGG,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)
      ENDIF

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DEPS,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      IF(ITURB <= 5) THEN
         ICASE = 2 ! k-epsilon model
      ELSE IF(ITURB < 9) THEN
         ICASE = 4 ! k-omega model
      ELSE         ! Spalart-Allmaras
         ICASE = 6
      ENDIF

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,ICASE,TRM)
      IF(ICYCLE == IPRINT .AND. FPRINTL) THEN
         WRITE(970+NGL+(M-1)*1000,*) 'ICYCLE = ', ICYCLE
         WRITE(970+NGL+(M-1)*1000,*)'                         '
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DEPS,
     +   IMAX,JMAX,KMAX,IN,JN,KN,970+NGL+(M-1)*1000)
      ENDIF

C ... Solve using an algebraic multigrid
      AN(1:NTOT) = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT) = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT) = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT) = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT) = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT) = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DEPS(1:NTOT) = DEPS(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT) = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DEPS,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,2,1,5,MGRID,9,0,NGL,ICYCLE)
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)
       DEPS(1:NTOT) = DTUR(1:NTOT)

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
      CALL PRINYS(3,DKC,   DRK, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      CALL PRINYS(3,DEPSC,DEPS, IT,IL,0,IK,IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      ENDIF ! ITURB >= 3      

C **********************************************************************
C ... Solve the transition model equations
C **********************************************************************

      IF(TRANSL) THEN ! Intermittency variables

C ... Solve G
      DH(1:NTOT)   = TRM(1:NTOT)%DG

      CALL DIFFTR(DIFF,VIS,EPS2,TRM,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,1)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,7,TRM)

      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT)  = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

C ... Solve using an algebraic multigrid
      DTUR(1:NTOT) = 0.
      CALL SOLAMG(DTUR,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     + IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)

      TRM(1:NTOT)%DG = DTUR(1:NTOT)

C ... Solve RET
      DH(1:NTOT)   = TRM(1:NTOT)%DRET

      CALL DIFFTR(DIFF,VIS,EPS2,TRM,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,2)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DH,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,C1,C2,KOVER,8,TRM)

      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DH(1:NTOT)  = DH(1:NTOT)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

C ... Solve using an algebraic multigrid
      DTUR(1:NTOT) = 0.
      CALL SOLAMG(DTUR,DH,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)
c     +    IN,JN,KN,MCYCAM,ITERAD,ITERAC,MGRID,ICONV3,0)

      TRM(1:NTOT)%DRET = DTUR(1:NTOT)

      ENDIF ! TRANSL

C **********************************************************************
C ... Solve the scalar equations. Reynolds stresses are not included to scalars
C **********************************************************************

      DO NS = 1,KSCAL

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS BEFORE ALGEBRAIC MULTIGRID'
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
      ENDIF

c      DTUR(1:NTOT) = 0. ! Vector for solution
      DIFF(1:NTOT) = 0.

      CALL DIFFSC(DIFF,VIS,EPS2,PSIGSC,PSIGS2,IMAX,JMAX,KMAX,
     +     NTOT,IN,JN,KN,NS)

      CALL AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,VOL,DIFF,
     +    A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,3,FRSDEN,
     +    U,V,W,A1XA,A1YA,A1ZA,A2XA,A2YA,A2ZA,A3XA,A3YA,A3ZA,DTUR,DN,DW,
     +    .FALSE.,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      CALL SOURSC(AP,SFI(1,NS),FI(1,NS),VOL,DTL,NTOT,NS,EPS6)

C ... Solve using an algebraic multigrid
      AN(1:NTOT)  = AN(1:NTOT)/(AP(1:NTOT)+EPS)
      AS(1:NTOT)  = AS(1:NTOT)/(AP(1:NTOT)+EPS)
      AW(1:NTOT)  = AW(1:NTOT)/(AP(1:NTOT)+EPS)
      AE(1:NTOT)  = AE(1:NTOT)/(AP(1:NTOT)+EPS)
      AT(1:NTOT)  = AT(1:NTOT)/(AP(1:NTOT)+EPS)
      AD(1:NTOT)  = AD(1:NTOT)/(AP(1:NTOT)+EPS)
      DTEMP(1:NTOT) = DFI(1:NTOT,NS)/(AP(1:NTOT)+EPS)
      AP(1:NTOT)  = 1.

      DTUR(1:NTOT) = 0.

      CALL SOLAMG(DTUR,DTEMP,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     +    IN,JN,KN,2,1,5,MGRID,8,0,NGL,ICYCLE)

      DFI(1:NTOT,NS) = DTUR(1:NTOT) * RO(1:NTOT) ! Definition of scalars

      IF(ICYCLE == IPRINT) THEN
       WRITE(3,*)'                         '
       WRITE(3,*)' RESULTS AFTER ALGEBRAIC MULTIGRID'
            CALL PRINYS(3,DFIC(NS),DFI(1,NS),IT,IL,0,IK,
     +           IMAX,JMAX,KMAX,NGL,M)
      ENDIF

      ENDDO ! KSCAL >= 0

C **********************************************************************
C ... Finally, store the pressure correction to 'DE'
C **********************************************************************

c      DE(1:NTOT)    = MIN(ABS(ALFAP*DP(1:NTOT)),ABS(0.05*P(1:NTOT)))
c      DE(1:NTOT)    = SIGN(DE(1:NTOT),DP(1:NTOT)) 

      DO I = 1,NTOT
         ALFAPP = ALFAP
         IF(RKSI(I) > 1.) ALFAPP = 1.
         DE(I)  = MIN(ABS(ALFAPP*DP(I)),ABS(0.05*P(I)))
         DE(I)  = SIGN(DE(I),DP(I))
         IF(M > 1) DE(I) = 0 ! Occasionally stabilizing
         PRC(I)%DP = DE(I)
      ENDDO

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DRK(1:NTOT)   = DRK(1:NTOT) *RO(1:NTOT)
      DEPS(1:NTOT)  = DEPS(1:NTOT)*RO(1:NTOT)
      ELSE
      DRK(1:NTOT)   = 0.
      DEPS(1:NTOT)  = 0.
      ENDIF

      IF(ITURB == 9) DRK(1:NTOT) = 0.
      RETURN
      END SUBROUTINE IMPSMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE IPSAP(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,
     +   ADL,APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DM,NTOT,IMAX,
     +   JMAX,KMAX,NGL,IPHASE,AP1,AP2,AP3,ALF1LG,ALF2LG,ALF3LG,VOL,DTL,
     +   RO)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND, IN, JN, KN, NPHASE
      USE CONSTANTS, ONLY : EPS10, EPS20

      IMPLICIT NONE

      INTEGER :: NTOT,N,IMAX,JMAX,KMAX,II,JJ,KK,IDIR,I,J,K,IL,ISTR,JSTR,
     +           NGL,IPHASE,IP1,IP2
      REAL :: DMNW,RKIF,RKIFD,RKIFE
      REAL :: AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),APL(*),ANL(*),
     +    ASL(*),AEL(*),AWL(*),ATL(*),ADL(*),APG(*),ANG(*),ASG(*),
     +    AEG(*),AWG(*),ATG(*),ADG(*),APUU(*),APUULG(NTOT,NPHASE),
     +    VOL(*),DTL(*),RO(*),
     +    DM(*),AP1(NTOT,NPHASE),AP2(NTOT,NPHASE),AP3(NTOT,NPHASE),
     +    ALF1LG(NTOT,NPHASE),ALF2LG(NTOT,NPHASE),ALF3LG(NTOT,NPHASE)
      TYPE(MPHASE_VARIABLES) VAR(*)

C ... Manipulate Jacobian matrices in interphase slip algorithm (IPSA)
C ... according to linearized interfacial friction term (xyz directions)
C ... for the pressure correction equation (diagonal terms only)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      IF(IPHASE == 1) THEN
         IP1 = 1
         IP2 = 2
      ELSEIF(IPHASE == 2) THEN
         IP1 = 2
         IP2 = 1
      ELSE
         STOP' No such phase in IPSAP'
      ENDIF

c      DO N = 1,NTOT
c         CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c         write(6644++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE
c         write(6647++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE

      DO IDIR = 1,3

      SELECT CASE(IDIR)

      CASE(1) ! X-direction

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1

         RKIF = VOL(I+II)*(2.*VAR(I+II)%RKIF(IDIR)) ! 1.E12
!     +        + MAX(VAR(I+II)%EVAPR(IPHASE),0.))
         RKIFD = RKIF !!+ VOL(I+II)*MAX(VAR(I+II)%EVAPR(2),0.)
!         RKIF = 0. !1.E12
        RKIFE = APUULG(I+II,IP2)/MAX(RKIFD,EPS20)

         AP1(I+II,IPHASE) =APUULG(I+II,IP1)
     +                + APUULG(I+II,IP2)/(RKIFE + 1.)
         ALF1LG(I+II,IPHASE)= VAR(I+II)%ALFA(IP1) + 
     +               VAR(I+II)%ALFA(IP2)/(RKIFE + 1.)
c      IF(IPHASE == 1) THEN
c      write(11000+idir,*) i+ii,AP1(I+II,IPHASE)
c      write(14000+IDIR,*) i+ii,ALF1LG(I+II,IPHASE)
c      ELSEIF(IPHASE == 2) THEN
c      write(15000+idir,*) i+ii,AP1(I+II,IPHASE)
c      write(18000+IDIR,*) i+ii,ALF1LG(I+II,IPHASE)
c      ENDIF
      ENDDO; ENDDO; ENDDO

      CASE(2) ! Y-direction

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1

         RKIF = VOL(I+II)*(2.*VAR(I+II)%RKIF(IDIR)) ! 1.E12
!     +        + MAX(VAR(I+II)%EVAPR(IPHASE),0.))
         RKIFD = RKIF !!+ VOL(I+II)*MAX(VAR(I+II)%EVAPR(2),0.)
        RKIFE = APUULG(I+II,IP2)/MAX(RKIFD,EPS20)

 !        RKIF = 0. !1.E12
         AP2(I+II,IPHASE) =APUULG(I+II,IP1)
     +                + APUULG(I+II,IP2)/(RKIFE + 1.)
         ALF2LG(I+II,IPHASE)= VAR(I+II)%ALFA(IP1) + 
     +               VAR(I+II)%ALFA(IP2)/(RKIFE + 1.)
c      IF(IPHASE == 1) THEN
c      write(12000+idir,*) i+ii,AP2(I+II,IPHASE)
c      write(14000+IDIR,*) i+ii,ALF2LG(I+II,IPHASE)
c      ELSEIF(IPHASE == 2) THEN
c      write(16000+idir,*) i+ii,AP2(I+II,IPHASE)
c      write(18000+IDIR,*) i+ii,ALF2LG(I+II,IPHASE)
c      ENDIF

      ENDDO; ENDDO; ENDDO

      CASE(3) ! Z-direction      

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1

         RKIF = VOL(I+II)*(2.*VAR(I+II)%RKIF(IDIR)) ! 1.E12
!     +        + MAX(VAR(I+II)%EVAPR(IPHASE),0.))
         RKIFD = RKIF !!+ VOL(I+II)*MAX(VAR(I+II)%EVAPR(2),0.)
!         RKIF = 0. !1.E12
       RKIFE = APUULG(I+II,IP2)/MAX(RKIFD,EPS20)

         AP3(I+II,IPHASE) =APUULG(I+II,IP1)
     +                + APUULG(I+II,IP2)/(RKIFE + 1.)
         ALF3LG(I+II,IPHASE)= VAR(I+II)%ALFA(IP1) + 
     +               VAR(I+II)%ALFA(IP2)/(RKIFE + 1.)
c      IF(IPHASE == 1) THEN
c      write(13000+idir,*) i+ii,AP3(I+II,IPHASE)
c      write(14000+IDIR,*) i+ii,ALF3LG(I+II,IPHASE)
c      ELSEIF(IPHASE == 2) THEN
c      write(17000+idir,*) i+ii,AP3(I+II,IPHASE)
c      write(18000+IDIR,*) i+ii,ALF3LG(I+II,IPHASE)
c      ENDIF

      ENDDO; ENDDO; ENDDO

      END SELECT

      ENDDO

      RETURN
      END SUBROUTINE IPSAP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE IPSA2(VAR,PRO,APMF,NTOT,IMAX,JMAX,KMAX,NGL,IPHASE,
     +     AP1X,AP2X,AP3X,ALF1X,ALF2X,ALF3X,VOL,DTL,RO)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND, IN, JN, KN, NPHASE
      USE CONSTANTS, ONLY : EPS6

      IMPLICIT NONE

      INTEGER :: NTOT,N,IMAX,JMAX,KMAX,II,JJ,KK,IDIR,I,J,K,IL,ISTR,JSTR,
     +           NGL,IPHASE
      REAL :: DMNW,RKIF,RKIFG,RKIFL
      REAL :: VOL(*),APMF(NTOT,IPHASE),RO(*),DTL(*),
     +    AP1X(NTOT,NPHASE),AP2X(NTOT,NPHASE),AP3X(NTOT,NPHASE),
     +    ALF1X(NTOT,IPHASE),ALF2X(NTOT,IPHASE),ALF3X(NTOT,IPHASE)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)

C ... Manipulate Jacobian matrices in interphase slip algorithm (IPSA)
C ... according to linearized interfacial friction term (xyz directions)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      DO N = 1,NTOT
c         CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c         write(6644++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE
c         write(6647++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      RKIFG = VOL(I+II)*(2.*VAR(I+II)%RKIF(1))/
     2        (PRO(I+II)%RO(2)*MAX(VAR(I+II)%ALFA(2),EPS6))
      RKIFL = VOL(I+II)*(2.*VAR(I+II)%RKIF(1))/
     2        (PRO(I+II)%RO(1)*MAX(VAR(I+II)%ALFA(1),EPS6))
         AP1X(I+II,1) = APMF(I+II,1) +
     +                      RKIFL/(1+RKIFG/APMF(I+II,2))
         AP2X(I+II,1) = AP1X(I+II,1)
         AP3X(I+II,1) = AP1X(I+II,1) ! Mersu

         ALF1X(I+II,1)= VAR(I+II)%ALFA(1) + 
     +                       VAR(I+II)%ALFA(2)/(1+APMF(I+II,2)/RKIFG)
         ALF2X(I+II,1)= ALF1X(I+II,1)
         ALF3X(I+II,1)= ALF1X(I+II,1)

         AP1X(I+II,2) = APMF(I+II,2) +
     +                      RKIFG/(1+RKIFL/APMF(I+II,1))
         AP2X(I+II,2) = AP1X(I+II,2)
         AP3X(I+II,2) = AP1X(I+II,2) ! Mersu

         ALF1X(I+II,2)= VAR(I+II)%ALFA(2) + 
     +                       VAR(I+II)%ALFA(1)/(1+APMF(I+II,1)/RKIFL)
         ALF2X(I+II,2)= ALF1X(I+II,2)
         ALF3X(I+II,2)= ALF1X(I+II,2)

      ENDDO; ENDDO; ENDDO

      RETURN
      END SUBROUTINE IPSA2
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE IPSA(VAR,AP,AN,AS,AE,AW,AT,AD,APL,ANL,ASL,AEL,AWL,ATL,
     +    ADL,APG,ANG,ASG,AEG,AWG,ATG,ADG,APUU,APUULG,DM,NTOT,IMAX,
     +    JMAX,KMAX,IDIR,NGL,IPHASE,APX,APY,APZ,ALFLG,VOL)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND, IN, JN, KN, NPHASE
      USE CONSTANTS, ONLY : EPS10, EPS20

      IMPLICIT NONE

      INTEGER :: NTOT,N,IMAX,JMAX,KMAX,II,JJ,KK,IDIR,I,J,K,IL,ISTR,JSTR,
     +           NGL,IPHASE
      REAL :: DMNW,RKIF,RKIFD,RKIFE
      REAL :: AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),APL(*),ANL(*),
     +    ASL(*),AEL(*),AWL(*),ATL(*),ADL(*),APG(*),ANG(*),ASG(*),
     +    AEG(*),AWG(*),ATG(*),ADG(*),APUU(*),APUULG(NTOT,NPHASE),
     +    VOL(*),
     +    DM(*),APX(NTOT,NPHASE),APY(NTOT,NPHASE),APZ(NTOT,NPHASE),
     +    ALFLG(NTOT,IPHASE)
      TYPE(MPHASE_VARIABLES) VAR(*)

C ... Manipulate Jacobian matrices in interphase slip algorithm (IPSA)
C ... according to linearized interfacial friction term (xyz directions)
C ... for the momentum equations

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      DO N = 1,NTOT
c         CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c         write(6644++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE
c         write(6647++NGL+(IPHASE-1)*10,*) 'IPHASE=',IPHASE

      IF(IPHASE == 1) THEN

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
c        CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c      write(666,*) i+ii,i,j,k,var(i+ii)%rkif(idir)
c         RKIF = MAX(1.E12,VAR(I+II)%RKIF(IDIR)) ! 1.E12
         RKIF = VOL(I+II)*(2.*VAR(I+II)%RKIF(IDIR)) ! 1.E12
!     +        + MAX(VAR(I+II)%EVAPR(IPHASE),0.))
         RKIFD = RKIF !!+ VOL(I+II)*MAX(VAR(I+II)%EVAPR(2),0.)
        RKIFE = APG(I+II)/MAX(RKIFD,EPS20)

!         RKIF  = 0.
c         RKIF = 1.E12
c       write(666+idir,*)i,j,k,apl(i+ii),apg(i+ii),
c     + APG(I+II)*RKIF/(APG(I+II) + RKIF),
c     + APL(I+II) + APG(I+II)*RKIF/(APG(I+II) + RKIF)
         AP(I+II)   = APL(I+II) + APG(I+II)/(RKIFE + 1.)
         AN(I+II)   = ANL(I+II) + ANG(I+II)/(RKIFE + 1.)
         AS(I+II)   = ASL(I+II) + ASG(I+II)/(RKIFE + 1.)
         AE(I+II)   = AEL(I+II) + AEG(I+II)/(RKIFE + 1.)
         AW(I+II)   = AWL(I+II) + AWG(I+II)/(RKIFE + 1.)
         AT(I+II)   = ATL(I+II) + ATG(I+II)/(RKIFE + 1.)
         AD(I+II)   = ADL(I+II) + ADG(I+II)/(RKIFE + 1.)
         APUU(I+II) = APUULG(I+II,1) + APUULG(I+II,2)
C ... APL/APG and APUULG are currently the same

c      write(1000+idir,*) i+ii,APX(I+II,IPHASE)
c      write(2000+idir,*) i+ii,APY(I+II,IPHASE)
c      write(3000+idir,*) i+ii,APZ(I+II,IPHASE)
c      write(4000+IDIR,*) i+ii,ALFLG(I+II,IPHASE)

c         write(6668+idir,*) i,j,APX(I+II,1),APY(I+II,1),APZ(I+II,1),
c     +   ALFLG(I+II,IPHASE)

c         if(VAR(I+II)%alfa(iphase) > .1) THEN
c        write(6050+idir,*) I,J,ALFLG(I+II,IPHASE),VAR(I+II)%alfa(2),
c     +   VAR(I+II)%ALFA(1)*RKIFD/(APL(I+II) + RKIF)
c         write(6060+idir,*) I,J,APX(I+II,IPHASE),
c     +   APY(I+II,IPHASE),APZ(I+II,IPHASE)
c         endif

c         write(6049,*) I,J,VAR(I+II)%ALFA(2)
c         write(6050+iphase,*) i,j,alflg(i+ii,iphase)
c          write(666,*) rkif,idir
         IF(IDIR  == 1) THEN ! Does not help
         DMNW    = VAR(I+II)%DM(1) + 
     +             VAR(I+II)%DM(2)/(RKIFE + 1.)  !-
!     +             APG(I+II)*(VAR(I+II)%U(1) - VAR(I+II)%U(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DM(1) 
c         write(6644+NGL+(IPHASE-1)*10,*) i,j,k,DMNW
c         write(6647+NGL+(IPHASE-1)*10,*) i,j,k,DM(I+II)
         ELSEIF(IDIR == 2) THEN
         DMNW    = VAR(I+II)%DN(1) + 
     +             VAR(I+II)%DN(2)/(RKIFE + 1.)  !-
!     +             APG(I+II)*(VAR(I+II)%V(1) - VAR(I+II)%V(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DN(1) 
c         DM(I+II) = DMNW
         ELSEIF(IDIR  == 3) THEN
         DMNW    = VAR(I+II)%DW(1) + 
     +             VAR(I+II)%DW(2)/(RKIFE + 1.)  !-
!     +             APG(I+II)*(VAR(I+II)%W(1) - VAR(I+II)%W(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DW(1) 
c         DM(I+II) = DMNW
         ENDIF
         DM(I+II) = DMNW

      ENDDO; ENDDO; ENDDO
c         write(6644++NGL+(IPHASE-1)*10,*) '----------------------'
c         write(6647++NGL+(IPHASE-1)*10,*) '----------------------'

      ELSEIF(IPHASE == 2) THEN

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
c        CALL IJKPAI(N,IMAX,JMAX,KMAX,II,JJ,KK)
c      write(667,*) i+ii,i,j,k,var(i+ii)%rkif(idir)

c         RKIF = MAX(1.E12,VAR(I+II)%RKIF(IDIR)) ! 1.E12
c         RKIF = VOL(I+II)*VAR(I+II)%RKIF(IDIR) ! 1.E12

         RKIF  = VOL(I+II)*(2.*VAR(I+II)%RKIF(IDIR)) ! 1.E12
!     +        + MAX(VAR(I+II)%EVAPR(IPHASE),0.))
         RKIFD = RKIF !!+ VOL(I+II)*MAX(VAR(I+II)%EVAPR(1),0.)
         RKIFE = APL(I+II)/MAX(RKIFD,EPS20)

!         RKIF  = 0.
c         RKIF = 1.E12
c         write(676+idir,*)i,j,k,apl(i+ii),apg(i+ii),
c     + APL(I+II)*RKIF/(APL(I+II) + RKIF),
c     + APL(I+II)*RKIF/(APL(I+II) + RKIF) + APG(I+II)
c         write(696+idir,*)i,j,k,anl(i+ii),ang(i+ii),
c     + ANL(I+II)*RKIF/(ANL(I+II) + RKIF),
c     + ANL(I+II)*RKIF/(ANL(I+II) + RKIF) + ANG(I+II)

         AP(I+II)   = APL(I+II)/(RKIFE + 1.) + APG(I+II)
         AN(I+II)   = ANL(I+II)/(RKIFE + 1.) + ANG(I+II)
         AS(I+II)   = ASL(I+II)/(RKIFE + 1.) + ASG(I+II)
         AE(I+II)   = AEL(I+II)/(RKIFE + 1.) + AEG(I+II)
         AW(I+II)   = AWL(I+II)/(RKIFE + 1.) + AWG(I+II)
         AT(I+II)   = ATL(I+II)/(RKIFE + 1.) + ATG(I+II)
         AD(I+II)   = ADL(I+II)/(RKIFE + 1.) + ADG(I+II)
         APUU(I+II) = APUULG(I+II,1) + APUULG(I+II,2)
c          write(667,*) rkif,idir

c         write(6050+iphase,*) i,j,alflg(i+ii,iphase)

c      write(5000+idir,*) i+ii,APX(I+II,IPHASE)
c      write(6000+idir,*) i+ii,APY(I+II,IPHASE)
c      write(7000+idir,*) i+ii,APZ(I+II,IPHASE)
c      write(8000+IDIR,*) i+ii,ALFLG(I+II,IPHASE)

         IF(IDIR  == 1) THEN
         DMNW    = VAR(I+II)%DM(1)/(RKIFE + 1.) + 
     +             VAR(I+II)%DM(2) !+
!     +             APL(I+II)*(VAR(I+II)%U(1) - VAR(I+II)%U(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DM(2) 

         ELSEIF(IDIR == 2) THEN
         DMNW    = VAR(I+II)%DN(1)/(RKIFE + 1.) + 
     +             VAR(I+II)%DN(2) !+
!     +             APL(I+II)*(VAR(I+II)%V(1) - VAR(I+II)%V(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DN(2) 

         ELSEIF(IDIR  == 3) THEN
         DMNW    = VAR(I+II)%DW(1)/(RKIFE + 1.) + 
     +             VAR(I+II)%DW(2) !+
!     +             APL(I+II)*(VAR(I+II)%W(1) - VAR(I+II)%W(2)) /
!     +             (RKIFE + 1.)!*.5
c         DMNW = VAR(I+II)%DW(2) 

         ENDIF
         DM(I+II) = DMNW

      ENDDO; ENDDO; ENDDO


      ELSE
         STOP 'No such phase in IPSA. Exiting...'
      ENDIF

      RETURN
      END SUBROUTINE IPSA
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATMF(AP,AN,AS,AE,AW,AT,AD,PRC,PRO,VAR,D1,D2,D3,RO,
     + DRO,VOL,DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,IPH,
     + FRSDEN,U,V,W,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,DU,DV,DW,
     + MOMCOL,RKSI,NCHIM,DT,TIMEL,APU,APUU,REFVEL)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND
      USE CONSTANTS, ONLY : EPS6, EPS10, EPS20

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,NCHIM,IPH
      REAL :: DDISTX,DDISTY,DDISTZ,RLAM,DTT,APURO,RMASS,FRSDEN,ABSVEL,
     + FORCE,ROIJ,RUIJ,RVIJ,RWIJ,EPS,DROII,DT,AAVA,ALFA1,ALFA2,BBB,
     + REFVEL
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),D1(*),D2(*),D3(*),RO(*),VOL(*),DIFF(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),C(*),U(*),V(*),W(*),A1X(*),A1Y(*),
     + A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DU(*),DV(*),
     + DW(*),RKSI(*),APU(*),APUU(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
      LOGICAL :: MOMCOL, TIMEL
    
C
C ... DISCRETIZE MOMENTUM EQUATIONS FOR A TWO-FLUID MODEL
 
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      IF(IEQ == 1) THEN
         DTT = 1.! Momentum
      ELSE IF(IEQ == 2) THEN
         DTT = 1. ! Energy
      ELSE IF(IEQ == 3) THEN
         DTT = 1. ! Turbulence
      ELSE
         WRITE(*,*) 'No such equation type in AAMATMF. Exiting...'
         STOP
      ENDIF
     
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
      EPS  = 1.E-20

C ... K-direction

      DO K = 0,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ
      AD(I+II+IL) = AT(I+II)
      AP(I+II)= AP(I+II) - AT(I+II)
      AP(I+II+IL)= AP(I+II+IL) - AT(I+II)
      ENDDO ; ENDDO ; ENDDO

C ... J-direction

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY
      AS(I+II+ISTR) = AN(I+II)
      AP(I+II)= AP(I+II) - AN(I+II)
      AP(I+II+ISTR)= AP(I+II+ISTR) - AN(I+II)

      ENDDO ; ENDDO ; ENDDO

C ... I-direction

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX

      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX
      AW(I+II+1) = AE(I+II)
c      IF(IEQ == 1) AE(I+II) = AE(I+II)-A1(I+II+1)*RO(I+II)*C(I+II)
      AP(I+II)= AP(I+II) - AE(I+II)
      AP(I+II+1)= AP(I+II+1) - AE(I+II)

      ENDDO ; ENDDO ; ENDDO


      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      AN(I+II)= AN(I+II) - MAX(0.,-VAR(I+II+ISTR)%F2R(IPH))
      AS(I+II)= AS(I+II) - MAX(0.,VAR(I+II)%F2R(IPH))
      AE(I+II)= AE(I+II) - MAX(0.,-VAR(I+II+1)%F1R(IPH))
      AW(I+II)= AW(I+II) - MAX(0.,VAR(I+II)%F1R(IPH))
      AT(I+II)= AT(I+II) - MAX(0.,-VAR(I+II+IL)%F3R(IPH))
      AD(I+II)= AD(I+II) - MAX(0.,VAR(I+II)%F3R(IPH))
      AP(I+II)= AP(I+II) 
     +  + MAX(0.,VAR(I+II+ISTR)%F2R(IPH))
     +  + MAX(0.,-VAR(I+II)%F2R(IPH))
     +  + MAX(0.,VAR(I+II+1)%F1R(IPH))
     +  + MAX(0.,-VAR(I+II)%F1R(IPH))
     +  + MAX(0.,VAR(I+II+IL)%F3R(IPH))
     +  + MAX(0.,-VAR(I+II)%F3R(IPH))


c      AN(I+II)= AN(I+II) - MAX(0.,-PRC(I+II+ISTR)%F2R)
c      AS(I+II)= AS(I+II) - MAX(0.,PRC(I+II)%F2R)
c      AE(I+II)= AE(I+II) - MAX(0.,-PRC(I+II+1)%F1R)
c      AW(I+II)= AW(I+II) - MAX(0.,PRC(I+II)%F1R)
c      AT(I+II)= AT(I+II) - MAX(0.,-PRC(I+II+IL)%F3R)
c      AD(I+II)= AD(I+II) - MAX(0.,PRC(I+II)%F3R)

c      IF(IEQ == 1) THEN
c      FORCE   = VOL(I+II)*ABS((RO(I+II)-FRSDEN)*(GX/(U(I+II)+1.E-5) + 
c     +          GY/(V(I+II)+1.E-5) + GZ/(W(I+II)+1.E-5))) *.05
c      ELSE
c      FORCE = 0.
c      ENDIF
c      ABSVEL = SQRT(VAR(I+II)%U(IPH)**2 + VAR(I+II)%V(IPH)**2 +
c     +              VAR(I+II)%W(IPH)**2 + EPS6)
      FORCE =PRO(I+II)%FITOT(1)**2+PRO(I+II)%FITOT(2)**2+
     +       PRO(I+II)%FITOT(3)**2 
      FORCE = 100.*VOL(I+II)*SQRT(FORCE)/REFVEL !/ABSVEL ! 265

c       ALFA1 = MAX(VAR(N)%ALFA(1),EPS6)
c       ALFA2 = MAX(VAR(N)%ALFA(2),EPS6)
c       BBB = (PRO(N)%RO(1)+PRO(N)%RO(2))*ALFA1*ALFA2


      APURO   = -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II)
     +          -AD(I+II)
      AP(I+II)= APURO +  ! Marginally better than MAX
     +  VOL(I+II)*RO(I+II)/DTL(I+II)*.5! Mersu pelaa
c     +  VOL(I+II)*BBB/DTL(I+II) ! Mersu pelaa
!     +  +VOL(I+II)*VAR(I+II)%RKIF(IPH)/.001
c     +    0.1*VAR(I+II)%ALFA(1)*PRO(I+II)%RO(1)*VOL(I+II)/DTL(I+II) *DTT
c     +  + VAR(I+II)%ALFA(IPH)*(PRO(I+II)%RO(1)+PRO(I+II)%RO(2))
c     2   *VOL(I+II)/DTL(I+II) *DTT 
     +  + VAR(I+II)%ALFA(IPH)*PRO(I+II)%RO(IPH)*VOL(I+II)/DTL(I+II)*.5!DTT 
c     +  + PRO(I+II)%RO(1)*VOL(I+II)/DTL(I+II)*DTT  ! Ei potkurilla
     +  + EPS10 !+ VOL(I+II)*ABS(VAR(I+II)%EVAPR(IPH))!
     +  + FORCE ! Test if not participates to pressure equation
c      write(999,*)i,j,force,ap(i+II),force/ap(i+II)
c      IF(IPH == 2) AP(I+II) = AP(I+II) + ABS(DRO(I+II)) +
c      write(6666+iph,*) i,j,apuro,apuro+
c    + VAR(I+II)%ALFA(IPH)*PRO(I+II)%RO(IPH)*VOL(I+II)/DTL(I+II)
c     + ap(i+ii),ABS(VAR(I+II)%DX(IPH))
!      AP(I+II) = AP(I+II) + ABS(VAR(I+II)%DX(IPH))
cc      IF(IPH == 2) AP(I+II) = AP(I+II) +MAX(
cc     + VAR(I+II)%X(IPH),1.)*RO(I+II)*VOL(I+II)/DTL(I+II) 
!     + - VOL(I+II)*MIN(VAR(I+II)%EVAPR(IPH),0.)
!!!      APUU(I+II) = AP(I+II)
      IF(IPH == 2) THEN
      IF(VAR(I+II)%ALFA(IPH) <= .5) THEN
         AAVA    = MAX(VAR(I+II)%ALFA(IPH),.1)
      ELSE
         AAVA    = MAX(VAR(I+II)%ALFA(1),  .1)
      ENDIF
cc      AP(I+II) = AP(I+II) + ! Oli AP AP
cc     + 0.1*PRO(I+II)%RO(1)*VOL(I+II)/DTL(I+II)
c     + 0.1*PRO(I+II)%RO(1)*VOL(I+II)/DTL(I+II)
      ENDIF

      APUU(I+II) = AP(I+II) 
      AP(I+II)   = AP(I+II) + FORCE ! Test if not participates pressure
c      if(iph == 2) APUU(I+II)= AP(I+II) -VOL(I+II)*RO(I+II)/DTL(I+II)

      IF(NCHIM > 0) THEN
         AP(I+II)= AP(I+II) + VOL(I+II)*RKSI(I+II)*RO(I+II)
      ENDIF
         APU(I+II) = AP(I+II)
      IF(TIMEL) THEN
         AP(I+II)= AP(I+II) + 1.5*VOL(I+II)*VAR(I+II)%ARO(IPH)/DT
      ENDIF
c      IF(IEQ == 1) PRC(I+II)%AP = APURO ! Store for the Rhie-Chow term
      
1000  CONTINUE

C ... Boundary cells (Von Neumann condition can be activated)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
c      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ-MAX(0.,-PRC(I+II+IL)%F3R)
      AP(I+II)= AP(I+II+IL)!-AT(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+IL)
      APUU(I+II) = AP(I+II)
      AT(I+II)= 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+IL)%AP ! Store AP
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-IL))
c      AD(I+II)= -RLAM*A3(I+II)/DDISTZ-MAX(0.,PRC(I+II)%F3R)
      AP(I+II)= AP(I+II-IL)!-AD(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-IL)
      APUU(I+II) = AP(I+II)
      AD(I+II)=  0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-IL)%AP ! Store AP
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
c      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY-MAX(0.,-PRC(I+II+ISTR)%F2R)
      AP(I+II)= AP(I+II+ISTR)!-AN(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+ISTR)
      APUU(I+II) = AP(I+II)
      AN(I+II) = 0. ! Dirichlet condition
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+ISTR)%AP ! Store AP
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-ISTR))
c      AS(I+II)= -RLAM*A2(I+II)/DDISTY-MAX(0.,PRC(I+II)%F2R)
      AP(I+II)= AP(I+II-ISTR)!-AS(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-ISTR)
      APUU(I+II) = AP(I+II)
      AS(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-ISTR)%AP ! Store AP
2300  CONTINUE
      DO 2400 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
c      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX-MAX(0.,-PRC(I+II+1)%F1R)
      AP(I+II)= AP(I+II+1) !-AE(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+1)
      APUU(I+II) = AP(I+II)
      AE(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+1)%AP ! Store AP
2400  CONTINUE
      DO 2500 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-1))
c      AW(I+II)= -RLAM*A1(I+II)/DDISTX-MAX(0.,PRC(I+II)%F1R)
      AP(I+II)= AP(I+II-1) !-AW(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-1)
      APUU(I+II) = AP(I+II)
      AW(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
c      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-1)%AP ! Store AP
2500  CONTINUE

c      APUU(1:NTOT) = APU(1:NTOT)

c      DO N = 1,NTOT
c      AP(N)= AP(N)/RO(N)
c      AE(N)= AE(N)/RO(N)
c      AW(N)= AW(N)/RO(N)
c      AN(N)= AN(N)/RO(N)
c      AS(N)= AS(N)/RO(N)
c      AT(N)= AT(N)/RO(N)
c      AD(N)= AD(N)/RO(N)
c      ENDDO
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   

      IF(MOMCOL) THEN ! Subtract the momentum residual

      STOP 'Not possible with the two-fluid model'

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX
         ROIJ = .5*(RO(I+II)+RO(I+II+1))
         DROII= 0.
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1X(I+II+1)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1Y(I+II+1)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1Z(I+II+1)
         DROII= DROII -.5*A1(I+II+1)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-1))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1Z(I+II)
         DROII= DROII +.5*A1(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         ROIJ = .5*(RO(I+II)+RO(I+II+ISTR))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2X(I+II+ISTR)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2Y(I+II+ISTR)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2Z(I+II+ISTR)
         DROII= DROII -.5*A2(I+II+ISTR)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-ISTR))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2Z(I+II)
         DROII= DROII +.5*A2(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         ROIJ = .5*(RO(I+II)+RO(I+II+IL))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3X(I+II+IL)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3Y(I+II+IL)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3Z(I+II+IL)
         DROII= DROII -.5*A3(I+II+IL)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-IL))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3Z(I+II)
         DROII= DROII +.5*A3(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         DRO(I+II) = DRO(I+II) + 1.*DROII ! Merkki (ei ollut) testattu
c         write(666,*) i+ii, DRO(I+II),DROII
      ENDDO
      ENDDO
      ENDDO
      ENDIF ! MOMCOL

      END SUBROUTINE AAMATMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATF(AP,AN,AS,AE,AW,AT,AD,PRC,VAR,D1,D2,D3,RO,DRO,
     + VOL,DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,
     + FRSDEN,U,V,W,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,DU,DV,DW,
     + MOMCOL,RKSI,NCHIM,DT,TIMEL,APU,APUU,NGL)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : GX, GY, GZ, GROUND

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,NCHIM,NGL
      REAL :: DDISTX,DDISTY,DDISTZ,RLAM,DTT,APURO,RMASS,FRSDEN,
     + FORCE,ROIJ,RUIJ,RVIJ,RWIJ,EPS,DROII,DT
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),D1(*),D2(*),D3(*),RO(*),VOL(*),DIFF(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),C(*),U(*),V(*),W(*),A1X(*),A1Y(*),
     + A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DU(*),DV(*),
     + DW(*),RKSI(*),APU(*),APUU(*)
      TYPE(PRE_COR) PRC(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      LOGICAL :: MOMCOL, TIMEL
    
C
C ... DISCRETIZE CONVECTION-DIFFUSION EQUATION
 
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      IF(IEQ == 1) THEN
         DTT = 1.! Momentum
      ELSE IF(IEQ == 2) THEN
         DTT = 1. ! Energy
      ELSE IF(IEQ == 3) THEN
         DTT = 1. ! Turbulence
      ELSE
         WRITE(*,*) 'No such equation type in AAMATF. Exiting...'
         STOP
      ENDIF
     
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
      EPS  = 1.E-20

C ... K-direction

      DO K = 0,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ
      AD(I+II+IL) = AT(I+II)

      ENDDO ; ENDDO ; ENDDO

C ... J-direction

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY
      AS(I+II+ISTR) = AN(I+II)

      ENDDO ; ENDDO ; ENDDO

C ... I-direction

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX

      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX
      AW(I+II+1) = AE(I+II)
c      IF(IEQ == 1) AE(I+II) = AE(I+II)-A1(I+II+1)*RO(I+II)*C(I+II)

      ENDDO ; ENDDO ; ENDDO


      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      AN(I+II)= AN(I+II) - MAX(0.,-PRC(I+II+ISTR)%F2R)
c      IF(IEQ == 1) AN(I+II) = AN(I+II)-A2(I+II+ISTR)*RO(I+II)*C(I+II)
      AS(I+II)= AS(I+II) - MAX(0.,PRC(I+II)%F2R)
c      IF(IEQ == 1) AS(I+II) = AS(I+II)-A2(I+II)*RO(I+II)*C(I+II)

      AE(I+II)= AE(I+II) - MAX(0.,-PRC(I+II+1)%F1R)
c      IF(IEQ == 1) AE(I+II) = AE(I+II)-A1(I+II+1)*RO(I+II)*C(I+II)

      AW(I+II)= AW(I+II) - MAX(0.,PRC(I+II)%F1R)
c      IF(IEQ == 1)AW(I+II) = AW(I+II)-A1(I+II)*RO(I+II)*C(I+II)

      AT(I+II)= AT(I+II) - MAX(0.,-PRC(I+II+IL)%F3R)
c      IF(IEQ == 1) AT(I+II) = AT(I+II)-A3(I+II+IL)*RO(I+II)*C(I+II)

      AD(I+II)= AD(I+II) - MAX(0.,PRC(I+II)%F3R)
c      IF(IEQ == 1) AD(I+II) = AD(I+II)-A3(I+II)*RO(I+II)*C(I+II)

c      IF(IEQ == 1) THEN
c      FORCE   = VOL(I+II)*ABS((RO(I+II)-FRSDEN)*(GX/(U(I+II)+1.E-5) + 
c     +          GY/(V(I+II)+1.E-5) + GZ/(W(I+II)+1.E-5))) *.05
c      ELSE
      FORCE = 0.
c      ENDIF

      APURO   = -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II)
     +          -AD(I+II)
      AP(I+II)= APURO + RO(I+II)*VOL(I+II)/DTL(I+II) *DTT
     +          +ABS(DRO(I+II)) + FORCE
      APUU(I+II) = AP(I+II)
      IF(NCHIM > 0) THEN
         AP(I+II)= AP(I+II) + VOL(I+II)*RKSI(I+II)*RO(I+II)
      ENDIF
         APU(I+II) = AP(I+II)
      IF(TIMEL) THEN
         AP(I+II)= AP(I+II) + 1.5*VOL(I+II)*RO(I+II)/DT
      ENDIF
      IF(IEQ == 1) PRC(I+II)%AP = APURO ! Store for the Rhie-Chow term
      
c     AP(I+II) = 2.*AP(I+II) ! Useless underrelaxation?
1000  CONTINUE

C ... Boundary cells (Von Neumann condition can be activated)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
c      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ-MAX(0.,-PRC(I+II+IL)%F3R)
      AP(I+II)= AP(I+II+IL)!-AT(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+IL)
      APUU(I+II) = AP(I+II)
      AT(I+II)= 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+IL)%AP ! Store AP
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-IL))
c      AD(I+II)= -RLAM*A3(I+II)/DDISTZ-MAX(0.,PRC(I+II)%F3R)
      AP(I+II)= AP(I+II-IL)!-AD(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-IL)
      APUU(I+II) = AP(I+II)
      AD(I+II)=  0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-IL)%AP ! Store AP
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
c      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY-MAX(0.,-PRC(I+II+ISTR)%F2R)
      AP(I+II)= AP(I+II+ISTR)!-AN(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+ISTR)
      APUU(I+II) = AP(I+II)
      AN(I+II) = 0. ! Dirichlet condition
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+ISTR)%AP ! Store AP
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-ISTR))
c      AS(I+II)= -RLAM*A2(I+II)/DDISTY-MAX(0.,PRC(I+II)%F2R)
      AP(I+II)= AP(I+II-ISTR)!-AS(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-ISTR)
      APUU(I+II) = AP(I+II)
      AS(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-ISTR)%AP ! Store AP
2300  CONTINUE
      DO 2400 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
c      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX-MAX(0.,-PRC(I+II+1)%F1R)
      AP(I+II)= AP(I+II+1) !-AE(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II+1)
      APUU(I+II) = AP(I+II)
      AE(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II+1)%AP ! Store AP
2400  CONTINUE
      DO 2500 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-1))
c      AW(I+II)= -RLAM*A1(I+II)/DDISTX-MAX(0.,PRC(I+II)%F1R)
      AP(I+II)= AP(I+II-1) !-AW(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      APU(I+II) = APU(I+II-1)
      APUU(I+II) = AP(I+II)
      AW(I+II) = 0.
      DU(I+II) = 0. ; DV(I+II) = 0. ; DW(I+II) = 0.
      IF(IEQ == 1) PRC(I+II)%AP = PRC(I+II-1)%AP ! Store AP
2500  CONTINUE

c      APUU(1:NTOT) = APU(1:NTOT)

c      DO N = 1,NTOT
c      AP(N)= AP(N)/RO(N)
c      AE(N)= AE(N)/RO(N)
c      AW(N)= AW(N)/RO(N)
c      AN(N)= AN(N)/RO(N)
c      AS(N)= AS(N)/RO(N)
c      AT(N)= AT(N)/RO(N)
c      AD(N)= AD(N)/RO(N)
c      ENDDO
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   

      IF(MOMCOL) THEN ! Subtract the momentum residual

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX
         ROIJ = .5*(RO(I+II)+RO(I+II+1))
         DROII= 0.
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1X(I+II+1)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1Y(I+II+1)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+1)/(AP(I+II+1)
     +            +EPS))*A1Z(I+II+1)
         DROII= DROII -.5*A1(I+II+1)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-1))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-1)/(AP(I+II-1)
     +            +EPS))*A1Z(I+II)
         DROII= DROII +.5*A1(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         ROIJ = .5*(RO(I+II)+RO(I+II+ISTR))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2X(I+II+ISTR)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2Y(I+II+ISTR)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+ISTR)/(AP(I+II+ISTR)
     +            +EPS))*A2Z(I+II+ISTR)
         DROII= DROII -.5*A2(I+II+ISTR)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-ISTR))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-ISTR)/(AP(I+II-ISTR)
     +            +EPS))*A2Z(I+II)
         DROII= DROII +.5*A2(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         ROIJ = .5*(RO(I+II)+RO(I+II+IL))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3X(I+II+IL)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3Y(I+II+IL)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II+IL)/(AP(I+II+IL)
     +            +EPS))*A3Z(I+II+IL)
         DROII= DROII -.5*A3(I+II+IL)*ROIJ*(RUIJ+RVIJ+RWIJ)
         ROIJ = .5*(RO(I+II)+RO(I+II-IL))
         RUIJ = (DU(I+II)/(AP(I+II)+EPS)+DU(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3X(I+II)
         RVIJ = (DV(I+II)/(AP(I+II)+EPS)+DV(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3Y(I+II)
         RWIJ = (DW(I+II)/(AP(I+II)+EPS)+DW(I+II-IL)/(AP(I+II-IL)
     +            +EPS))*A3Z(I+II)
         DROII= DROII +.5*A3(I+II)*ROIJ*(RUIJ+RVIJ+RWIJ)

         DRO(I+II) = DRO(I+II) + 1.*DROII ! Merkki (ei ollut) testattu
c         write(666,*) i+ii, DRO(I+II),DROII
      ENDDO
      ENDDO
      ENDDO
      ENDIF ! MOMCOL

      END SUBROUTINE AAMATF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMAMF(AP,AN,AS,AE,AW,AT,AD,PRO,VAR,D1,D2,D3,RO,DRO,
     + VOL,DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,M,IN,JN,KN,IEQ,IPH,
     + IEVAP,RKSI,NCHIM,ICYCLE,u,v,w,pdiff,DRHOL,DRHOG,APU,TIMEL,DT)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : EPS10

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,IPH,IEVAP,NCHIM
      REAL DDISTX,DDISTY,DDISTZ,RLAM,DTT,DQIFDT,DHMAX,ALMIN,ALMAX,DIAGH,
     + AAVE,apu1,apu2,HSEFF,EQEFF,DT,APURO
      real u(*),v(*),W(*),pdiff(*)
      INTEGER N,M,ISTR,JSTR,IL,I,J,K,II,JJ,ICYCLE
      REAL AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),D1(*),D2(*),D3(*),RO(*),VOL(*),DIFF(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),C(*),RKSI(*),DRHOL(*),DRHOG(*),APU(*)
      LOGICAL TIMEL
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
C
C ... DISCRETIZE CONVECTION-DIFFUSION EQUATION FOR HOMOGENEOUS
C ... MULTI-PHASE FLOW (CAVIT)
C
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      IF(IEQ == 1) THEN
         DTT = 1.! Momentum
      ELSE IF(IEQ == 2) THEN
         DTT = 1. ! Energy
      ELSE IF(IEQ == 3) THEN
         DTT = 1. ! Turbulence
      ELSE
         WRITE(*,*) 'No such equation type in AAMAMF. Exiting...'
         STOP
      ENDIF

c      IF(IPH == 2) DTT = 10. ! Stabiilimpi? (Ei vaikuta, uudellakaan testillä)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      ALMIN   = 1.E-5
      ALMAX   = .99999
      DQIFDT  = 0.
      AP(I+II)= 0.
      apu1 = 0. ; apu2 = 0. ; APURO = 0.

      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY-MAX(0.,-VAR(I+II+ISTR)%
     +                                                       F2R(IPH))
      AP(I+II)=  RLAM*A2(I+II+ISTR)/DDISTY+MAX(0.,VAR(I+II+ISTR)%
     +                                                       F2R(IPH))

      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-ISTR))
      AS(I+II)= -RLAM*A2(I+II)/DDISTY-MAX(0.,VAR(I+II)%F2R(IPH))
      AP(I+II)= AP(I+II) +
     + RLAM*A2(I+II)/DDISTY + MAX(0.,-VAR(I+II)%F2R(IPH))

      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX-MAX(0.,-VAR(I+II+1)%F1R(IPH))
      AP(I+II)= AP(I+II) +
     + RLAM*A1(I+II+1)/DDISTX+MAX(0.,VAR(I+II+1)%F1R(IPH))


      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-1))
      AW(I+II)= -RLAM*A1(I+II)/DDISTX-MAX(0.,VAR(I+II)%F1R(IPH))
      AP(I+II)= AP(I+II) +
     + RLAM*A1(I+II)/DDISTX+MAX(0.,-VAR(I+II)%F1R(IPH))

      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ-MAX(0.,-VAR(I+II+IL)%F3R(IPH))
      AP(I+II)= AP(I+II) +
     + RLAM*A3(I+II+IL)/DDISTZ + MAX(0.,VAR(I+II+IL)%F3R(IPH))

      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-IL))
      AD(I+II)= -RLAM*A3(I+II)/DDISTZ-MAX(0.,VAR(I+II)%F3R(IPH))
      AP(I+II)= AP(I+II) +
     + RLAM*A3(I+II)/DDISTZ + MAX(0.,-VAR(I+II)%F3R(IPH))


      IF(IPH == 1) THEN
      DHMAX = PRO(I+II)%CP(IPH)*ABS(PRO(I+II)%TSAT-PRO(I+II)%DTEMP(IPH))
     +  *.1 + 0.0001   ! 10. ! .5 ja 0.0001
      ELSE IF(IPH == 2) THEN
      DHMAX = PRO(I+II)%CP(IPH)*ABS(PRO(I+II)%TSAT-PRO(I+II)%DTEMP(IPH))
     +  *.1 + 0.0001  ! 10. ! .5 ja 0.0001
      ENDIF

      ALMIN   = MAX(VAR(I+II)%ALFA(IPH),ALMIN)
      ALMAX   = 1. - ALMIN !MIN(VAR(I+II)%ALFA(IPH),ALMAX)
cc      DQIFDT  = VOL(I+II)*(VAR(I+II)%ALFA(IPH)*PRO(I+II)%RO(IPH) /

      SELECT CASE(IEVAP)

      CASE(0) ! No phase change
         DQIFDT = 0.
      CASE(1,6) ! Cavitation model of Merkle (default)
c         DQIFDT =VOL(I+II)*(VAR(N)%EVAPR(IPH)+5.*ABS(VAR(N)%EVAPR(IPH)))
c         IF(IPH == 2) DQIFDT = DQIFDT + VOL(I+II)*VAR(N)%ALFA(2)*1.E6
c         DIAGH  = ABS(VAR(I+II)%EVAPR(IPH)) ! /DHMAX ! Mistä??
         DIAGH  = MAX(0.,VAR(I+II)%EVAPR(IPH))
         DQIFDT = VOL(I+II)*(VAR(I+II)%EQL(IPH)/PRO(I+II)%CP(IPH)+DIAGH)

      CASE(2) ! Cavitation model of Merkle with modified condensation
         DQIFDT = 0.
      CASE(3) ! Houdayer
         DQIFDT  = VOL(I+II)*(ALMIN*ALMAX*PRO(I+II)%RO(IPH) /
     +            PRO(I+II)%TAUF(IPH))
c020206     +        + ((ABS(VAR(I+II)%EVAPR(IPH))))/DHMAX)*2.  !+
c     +          ABS(PRO(I+II)%QIF(IPH)/PRO(I+II)%HSAT(IPH)))/
c     +          (1.*PRO(I+II)%CP(IPH))*PRO(I+II)%HSAT(IPH)))
      CASE(4) ! 2012 evaporation model
         DQIFDT  = VOL(I+II)*VAR(I+II)%EQL(IPH)
c     +             +  ABS(VAR(N)%EVAPR(IPH))/DHMAX)
      CASE(5,7,8) ! Merkle with temperatures
         APU1 = 1.
         IF(IPH == 2) 
     +  APU1 = 1. !+
c     + VOL(I+II)*VAR(I+II)%EQL(IPH)/(PRO(I+II)%RO(2)*(1.E-6+APU(I+II)))
c     + VOL(I+II)*VAR(I+II)%EQL(1)*PRO(I+II)%HSAT(1)*PRO(I+II)%HSAT(2)/
c     + PRO(I+II)%CP(1)/(1.E-6+APU(I+II))
         DIAGH  = ABS(VAR(I+II)%EVAPR(IPH))/DHMAX*PRO(I+II)%HSAT(IPH) ! Mistä??
c         DIAGH  = MAX(0.,VAR(I+II)%EVAPR(IPH))
         DQIFDT = VOL(I+II)*(VAR(I+II)%EQL(IPH)/PRO(I+II)%CP(IPH)*APU1
     +  + DIAGH)
      END SELECT

      AP(I+II) = AP(I+II) + DQIFDT
     +          +(1.E-4+VAR(I+II)%ALFA(IPH))*PRO(I+II)%RO(IPH)*VOL(I+II)
     +          / DTL(I+II) *DTT! *.5
!     +     +  VOL(I+II)*RO(I+II)/DTL(I+II)*.5! Mersu pelaako. Ei CAVIT


      IF(IPH == 1) THEN
         AP(I+II) = AP(I+II) + ABS(DRHOL(I+II))
      ELSE IF(IPH == 2) THEN
         AP(I+II) = AP(I+II) + ABS(DRHOG(I+II))
      ENDIF 
                  apu(i+ii) = ap(i+ii) ! talteen
       IF(NCHIM > 0) THEN
          AP(I+II) = AP(I+II) + VOL(I+II)*RKSI(I+II)
       ENDIF
       IF(TIMEL) THEN
          AP(I+II) = AP(I+II) + 1.5*VOL(I+II)*VAR(I+II)%ARO(IPH)/DT
       ENDIF

1000  CONTINUE

C ... Boundary cells (Von Neumann condition can be activated)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+IL))
c      AT(I+II)= -RLAM*A3(I+II+IL)/DDISTZ-MAX(0.,-VAR(I+II+IL)%F3R(IPH))
      AP(I+II)= AP(I+II+IL)!-AT(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AT(I+II) = 0.
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-IL))
c      AD(I+II)= -RLAM*A3(I+II)/DDISTZ-MAX(0.,VAR(I+II)%F3R(IPH))
      AP(I+II)= AP(I+II-IL)!-AD(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AD(I+II) = 0.
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+ISTR))
c      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY-MAX(0.,-VAR(I+II+ISTR)%F2R(IPH))
      AP(I+II)= AP(I+II+ISTR)!-AN(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AN(I+II) = 0. ! Dirichlet condition
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-ISTR))
c      AS(I+II)= -RLAM*A2(I+II)/DDISTY-MAX(0.,VAR(I+II)%F2R(IPH))
      AP(I+II)= AP(I+II-ISTR)!-AS(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AS(I+II) = 0.
2300  CONTINUE
      DO 2400 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II+1))
c      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX-MAX(0.,-VAR(I+II+1)%F1R(IPH))
      AP(I+II)= AP(I+II+1) !-AE(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0.
2400  CONTINUE
      DO 2500 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(DIFF(I+II)+DIFF(I+II-1))
c      AW(I+II)= -RLAM*A1(I+II)/DDISTX-MAX(0.,VAR(I+II)%F1R(IPH))
      AP(I+II)= AP(I+II-1) !-AW(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.
2500  CONTINUE
c      DO N = 1,NTOT
c      AP(N)= AP(N)/RO(N)
c      AE(N)= AE(N)/RO(N)
c      AW(N)= AW(N)/RO(N)
c      AN(N)= AN(N)/RO(N)
c      AS(N)= AS(N)/RO(N)
c      AT(N)= AT(N)/RO(N)
c      AD(N)= AD(N)/RO(N)
c      ENDDO
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   
      END SUBROUTINE AAMAMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMAMA(AP,AN,AS,AE,AW,AT,AD,VAR,PRC,D1,D2,D3,RO,DRO,
     + VOL,DIFF,A1,A2,A3,DTL,C,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,
     + DRO1,RKSI,NCHIM,VOIDTT,DT,TIMEL,MULPHC,NPHASE,NGL,VIS,FREDIF,PRO)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IEQ,IPHASE,NCHIM,NPHASE,NGL
      REAL DDISTX,DDISTY,DDISTZ,RLAM,DTT,QLMAX,AADIA,QLMIN,VOIDTT,DT,
     + FREDIF,DISTJ,DISTJP,VISJ,VISJP,APU,APU1,APU2,ROPJ,ROPJP
      INTEGER N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),D1(*),D2(*),D3(*),RO(*),VOL(*),DIFF(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),C(*),DRO1(*),RKSI(*),VIS(*)
      LOGICAL TIMEL
      CHARACTER(10)          :: MULPHC

      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(PRE_COR) PRC(*)
C
C ... DISCRETIZE THE QUALITY SOLUTION FOR HOMOGENEOUS MULTIPHASE FLOW (CAVIT)

      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      IF(IEQ == 1) THEN
         DTT = 1.! Momentum
      ELSE IF(IEQ == 2) THEN
         DTT = 1. ! Energy
      ELSE IF(IEQ == 3) THEN
         DTT = 1. ! Turbulence
      ELSE IF(IEQ == 4) THEN
         DTT = VOIDTT ! Void
      ELSE
         WRITE(*,*) 'No such equation type in AAMAMA. Exiting...'
         STOP
      ENDIF

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      IF(MULPHC == 'CAVIT') THEN

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... Mersu, laskee kahdesti
c      if(abs(dro1(i+ii))  <= 1.e-9) dro1(i+ii) = 0.
      AN(I+II)= -MAX(0.,-PRC(I+II+ISTR)%F2R)
      AS(I+II)= -MAX(0., PRC(I+II)%F2R)
      AE(I+II)= -MAX(0.,-PRC(I+II+1)%F1R)
      AW(I+II)= -MAX(0., PRC(I+II)%F1R)
      AT(I+II)= -MAX(0.,-PRC(I+II+IL)%F3R)
      AD(I+II)= -MAX(0., PRC(I+II)%F3R)
c      AN(I+II)= -MAX(0.,-VAR(I+II+ISTR)%F2A(IPHASE))
c      AS(I+II)= -MAX(0., VAR(I+II)%F2A(IPHASE))
c      AE(I+II)= -MAX(0.,-VAR(I+II+1)%F1A(IPHASE))
c      AW(I+II)= -MAX(0., VAR(I+II)%F1A(IPHASE))
c      AT(I+II)= -MAX(0.,-VAR(I+II+IL)%F3A(IPHASE))
c      AD(I+II)= -MAX(0., VAR(I+II)%F3A(IPHASE))

      ENDDO; ENDDO; ENDDO

      ELSEIF(MULPHC == 'MULTI') THEN

      DO IPHASE = 1,NPHASE

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... Mersu, laskee kahdesti

      AN(I+II)= -MAX(0.,-VAR(I+II+ISTR)%F2X(IPHASE)) + AN(I+II)
      AS(I+II)= -MAX(0., VAR(I+II)%F2X(IPHASE))      + AS(I+II)
      AE(I+II)= -MAX(0.,-VAR(I+II+1)%F1X(IPHASE))    + AE(I+II)
      AW(I+II)= -MAX(0., VAR(I+II)%F1X(IPHASE))      + AW(I+II)
      AT(I+II)= -MAX(0.,-VAR(I+II+IL)%F3X(IPHASE))   + AT(I+II)
      AD(I+II)= -MAX(0., VAR(I+II)%F3X(IPHASE))      + AD(I+II)

      ENDDO; ENDDO; ENDDO
      ENDDO ! IPHASE

      ELSE

         STOP 'No such choice in AAMAMA. Exiting...'

      ENDIF ! MULPHC

      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

      QLMIN   = 1.E-6  !1.E-12
c      IF(VAR(I+II)%EVAPR(1) > 0.) QLMIN = 1.E-12
      QLMAX   = MIN(VAR(I+II)%X(1),VAR(I+II)%X(2))*.05 !+ 1.E-6
      IF(VAR(I+II)%X(1) < .5 .AND. VAR(I+II)%EVAPR(2) > 0.) 
     +QLMAX   = QLMAX*16.*VAR(I+II)%X(1)**4      
      AADIA   = VOL(I+II)*ABS(VAR(I+II)%EVAPR(1))/(QLMAX+QLMIN) ! Oli 1.E-6
      AP(I+II)= - AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II)
     +          - AD(I+II) !+ ABS(DRO1(I+II))
     +          + RO(I+II)*VOL(I+II)/DTL(I+II) *DTT
     +          + AADIA
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + VOL(I+II)*RKSI(I+II)
      ENDIF
      IF(TIMEL) THEN
         AP(I+II) = AP(I+II) + 1.5*VOL(I+II)*RO(I+II)/DT
      ENDIF
c      AP(I+II) = MAX(AP(I+II),1.E-10)
c      AP(I+II) = 5.*AP(I+II)
1000  CONTINUE

C ... Boundary cells (Von Neumann condition can be activated)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
      AP(I+II)= AP(I+II+IL)!-AT(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AT(I+II) = 0.
      DRO1(I+II) = 0.
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
      AP(I+II)= AP(I+II-IL)!-AD(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AD(I+II) = 0.
      DRO1(I+II) = 0.
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
      AP(I+II)= AP(I+II+ISTR)!-AN(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AN(I+II) = 0. ! Dirichlet condition
      DRO1(I+II) = 0.
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
      AP(I+II)= AP(I+II-ISTR)!-AS(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AS(I+II) = 0.
      DRO1(I+II) = 0.
2300  CONTINUE
      DO 2400 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      AP(I+II)= AP(I+II+1) !-AE(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0. !-AP(I+II) !0.
      DRO1(I+II) = 0.
2400  CONTINUE
      DO 2500 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
      AP(I+II)= AP(I+II-1) !-AW(I+II) !+ RO(I+II)*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.! -AP(I+II)  !0.
      DRO1(I+II) = 0.
2500  CONTINUE
c      DO N = 1,NTOT
c      AP(N)= AP(N)/RO(N)
c      AE(N)= AE(N)/RO(N)
c      AW(N)= AW(N)/RO(N)
c      AN(N)= AN(N)/RO(N)
c      AS(N)= AS(N)/RO(N)
c      AT(N)= AT(N)/RO(N)
c      AD(N)= AD(N)/RO(N)
c      ENDDO
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   
      END SUBROUTINE AAMAMA
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATP_DIST(AP,AN,AS,AE,AW,AT,AD,APU,PRC,PRO,VAR,RO,
     + DRO,VOL,P,C,A1,A2,A3,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     + XC,YC,ZC,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAP,M,NGL,
     + D1,D2,D3)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      INTEGER :: N,M,ISTR,JSTR,IL,I,J,K,II,JJ,NGL

      REAL    :: ALFAP,ALFAPT,APJ,ROJ,CCC,YPC1,YPC2,GAMMA,DPMAX,EE,SLEN,
     +           SXI,SYI,SZI,ALFA

      REAL    :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),APU(*),RO(*),VOL(*),C(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),P(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),D1(*),D2(*),D3(*)

      REAL :: XC(*),YC(*),ZC(*)

      TYPE(PRE_COR) PRC(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(BLOCKS) BLKS(*)

C
C ... DISCRETIZE POISSON EQUATION

      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      ALFAPT = .000001*M  !Vedellä nolla pois
     	ALFAPT = .01*M

      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti
       EE = 0.


      SXI    = XC(I+II+ISTR) - XC(I+II)
      SYI    = YC(I+II+ISTR) - YC(I+II)
      SZI    = ZC(I+II+ISTR) - ZC(I+II)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN
      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A2X(I+II+ISTR)+SYI*A2Y(I+II+ISTR)+SZI*A2Z(I+II+ISTR)
c      if(alfa <= 0.)       write(777,*) i+ii,alfa,slen,'north'
c      APJ     = .5*(VOL(I+II)/APU(I+II) + VOL(I+II+ISTR)/APU(I+II+ISTR))
      APJ     = .5*(1./APU(I+II) + 1./APU(I+II+ISTR))
      ROJ     = .5*(RO(I+II)  + RO(I+II+ISTR))
      SLEN=.5*(D2(I+II)+D2(I+II+ISTR))
      AN(I+II)= -A2(I+II+ISTR)**2*ROJ*APJ !/SLEN
c     + *.5*(VOL(I+II)+VOL(I+II+ISTR))

      SXI    = XC(I+II) - XC(I+II-ISTR)
      SYI    = YC(I+II) - YC(I+II-ISTR)
      SZI    = ZC(I+II) - ZC(I+II-ISTR)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN

      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A2X(I+II)+SYI*A2Y(I+II)+SZI*A2Z(I+II)
c      if(alfa <= 0.)        write(777,*) i+ii,alfa,slen,'south'

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      ROJ     = .5*(RO(I+II)  + RO(I+II-ISTR))
      SLEN=.5*(D2(I+II)+D2(I+II-ISTR))

      AS(I+II)= -A2(I+II)**2*ROJ*APJ!/SLEN
c     + *.5*(VOL(I+II)+VOL(I+II-ISTR))


      SXI    = XC(I+II+1) - XC(I+II)
      SYI    = YC(I+II+1) - YC(I+II)
      SZI    = ZC(I+II+1) - ZC(I+II)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN
      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A1X(I+II+1)+SYI*A1Y(I+II+1)+SZI*A1Z(I+II+1)
c      if(alfa <= 0.)        write(777,*) i+ii,alfa,slen,'east'

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II+1))
      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
c      SLEN=.5*(D1(I+II)+D1(I+II+1))

      AE(I+II)= -A1(I+II+1)**2*ROJ*APJ !/SLEN
!     + *.5*(VOL(I+II)+VOL(I+II+1))
!       if(i == 1) AE(I+II) = 0.


      SXI    = XC(I+II) - XC(I+II-1)
      SYI    = YC(I+II) - YC(I+II-1)
      SZI    = ZC(I+II) - ZC(I+II-1)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN
      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A1X(I+II)+SYI*A1Y(I+II)+SZI*A1Z(I+II)
c      if(alfa <= 0.)       write(777,*) i+ii,alfa,slen,'west'

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
      SLEN=.5*(D1(I+II)+D1(I+II-1))

      AW(I+II)= -A1(I+II)**2*ROJ*APJ!/SLEN
!     + *.5*(VOL(I+II)+VOL(I+II-1))


      SXI    = XC(I+II+IL) - XC(I+II)
      SYI    = YC(I+II+IL) - YC(I+II)
      SZI    = ZC(I+II+IL) - ZC(I+II)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN
      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A3X(I+II+IL)+SYI*A3Y(I+II+IL)+SZI*A3Z(I+II+IL)
c      if(alfa <= 0.)        write(777,*) i+ii,alfa,slen,'top'

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II+IL))
      ROJ     = .5*(RO(I+II)  + RO(I+II+IL))
      SLEN=.5*(D3(I+II)+D3(I+II+IL))

      AT(I+II)= -A3(I+II+IL)**2*ROJ*APJ!/SLEN
!     + *.5*(VOL(I+II)+VOL(I+II+IL))


      SXI    = XC(I+II) - XC(I+II-IL)
      SYI    = YC(I+II) - YC(I+II-IL)
      SZI    = ZC(I+II) - ZC(I+II-IL)
      SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
      SXI    = SXI/SLEN
      SYI    = SYI/SLEN
      SZI    = SZI/SLEN
      ALFA   = SXI*A3X(I+II)+SYI*A3Y(I+II)+SZI*A3Z(I+II)
c      if(alfa <= 0.)        write(777,*) i+ii,alfa,slen,'down'

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
      SLEN=.5*(D3(I+II)+D3(I+II-IL))

      AD(I+II)= -A3(I+II)**2*ROJ*APJ!/SLEN
!     + *.5*(VOL(I+II)+VOL(I+II-IL))


c ... test
c          stop 'aamatp'
      CCC     = C(I+II)
      GAMMA   = 0.
      DPMAX   = 1000.
      ee      = 0.

      IF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN   !IF(NPHASES > 0) THEN
      DPMAX   = 1.*ABS(P(I+II)-PRO(I+II)%PSAT) + 100.
      GAMMA   = VOL(I+II)*ABS(VAR(I+II)%EVAPR(1))*(-1./PRO(I+II)%RO(1)+
     +        1./PRO(I+II)%RO(2))/DPMAX
cc       GAMMA  = VOL(I+II)*(1./PRO(I+II)%DHSDP(1)+1./PRO(I+II)%DHSDP(2))/
cc     +  (PRO(I+II)%HSAT(2) - PRO(I+II)%HSAT(1))
       IF(ABS(VAR(I+II)%EVAPR(1)) < 1.E-6) Gamma = 0.
      ENDIF
            
      EE = EE + VOL(I+II)/DTL(I+II) 
      AP(I+II)= -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II)
     + - AD(I+II) + 1./(CCC**2)*EE*ALFAPT
     + + GAMMA    ! Vesi
cc     + -AD(I+II) + 1./(450.**2*RO(I+II))*VOL(I+II)/DTL(I+II)*.0001 ! höyry
c     + + 1./450.**2*VOL(I+II)/1.
cc      AP(I+II)= AP(I+II) + ABS(DRO(I+II))/1.E5
      AP(I+II)= AP(I+II)/ALFAP

1000  CONTINUE
         
C ... Boundary cells (Von Neumann condition for the implicit stage)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+IL))
      APJ     = 1./APU(I+II+IL) ! APU(I+II) N.A. everywhere along the boundary
      ROJ     = .5*(RO(I+II)  + RO(I+II+IL))
c      AT(I+II)= AD(I+II+IL)
      AT(I+II)= -A3(I+II+IL)**2*ROJ*APJ
      AP(I+II)= -AT(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AT(I+II) = 0.
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      APJ     = 1./APU(I+II-IL)
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
c      AD(I+II)= AT(I+II-IL)
      AD(I+II)= -A3(I+II)**2*ROJ*APJ
      AP(I+II)= -AD(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AD(I+II) = 0.
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+ISTR))
      APJ     = 1./APU(I+II+ISTR)
      ROJ     = .5*(RO(I+II)  + RO(I+II+ISTR))
c      AN(I+II)= AS(I+II+ISTR)
      AN(I+II)= -A2(I+II+ISTR)**2*ROJ*APJ
      AP(I+II)= -AN(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AN(I+II) = 0. ! Dirichlet condition
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      APJ     = 1./APU(I+II-ISTR)
      ROJ     = .5*(RO(I+II)  + APU(I+II-ISTR))
c      AS(I+II)= AN(I+II-ISTR)
      AS(I+II)= -A2(I+II)**2*ROJ*APJ
      AP(I+II)= -AS(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AS(I+II) = 0.
2300  CONTINUE
      DO 2400 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+1))
      APJ     = 1./APU(I+II+1)
      ROJ     = .5*(RO(I+II)  + APU(I+II+1))
c      AE(I+II)= AW(I+II+1)
      AE(I+II)= -A1(I+II+1)**2*ROJ*APJ
      AP(I+II)= -AE(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0.
2400  CONTINUE
      DO 2500 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      APJ     = 1./APU(I+II-1)
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
c      AW(I+II)= AE(I+II-1)
      AW(I+II)= -A1(I+II)**2*ROJ*APJ
      AP(I+II)= -AW(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.
2500  CONTINUE
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   
      END SUBROUTINE AAMATP_DIST
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATP(AP,AN,AS,AE,AW,AT,AD,APU,PRC,PRO,VAR,RO,DRO,VOL,
     + P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAP,M,NGL,
     + XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NCHIM
      INTEGER :: N,M,ISTR,JSTR,IL,I,J,K,II,JJ,NGL

      REAL :: ALFAP,ALFAPT,APJ,ROJ,CCC,YPC1,YPC2,GAMMA,DPMAX,EE,
     +        GAMMA2,DT

      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),APU(*),RO(*),VOL(*),C(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),P(*),XFC(*),YFC(*),ZFC(*),PDIFF(*),
     + RKSI(*)

      TYPE(PRE_COR) PRC(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)

      
      LOGICAL :: COMPCORR,TIMEL

C
C ... DISCRETIZE POISSON EQUATION
C

      COMPCORR = BLKS(NGL)%COMPCORR
 
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      ALFAPT = .000001*M  !Vedellä nolla pois
      ALFAPT = .01*M ! Oli 0.01*M

C ... Symmetric part of the Poisson equation

      DO K = 1,KMAX + 1
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
      AD(I+II)= -A3(I+II)**2*ROJ*APJ
      AT(I+II-IL) = AD(I+II) ! uusi

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      ROJ     = .5*(RO(I+II)  + RO(I+II-ISTR))
      AS(I+II)= -A2(I+II)**2*ROJ*APJ
      AN(I+II-ISTR) = AS(I+II) ! uusi

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX + 1

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
      AW(I+II)= -A1(I+II)**2*ROJ*APJ  
      AE(I+II-1) = AW(I+II)

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

C ... Diagonal element is added together

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      AP(I+II) = -AN(I+II)-AS(I+II)-AE(I+II)-AW(I+II)-AT(I+II)-AD(I+II)
c      write(7667,*) AP(I+II),AT(I+II),AD(I+II),AN(I+II),AS(I+II),
c     +              AE(I+II),AW(I+II)

      ENDDO ; ENDDO ; ENDDO

C ... Compressibility correction and mass-transfer relaxation

      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      EE = 0.

      IF(COMPCORR) THEN

      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
c      ROJ = RO(I+II)
c      IF(PRC(I+II+ISTR)%F2R < 0.) ROJ = RO(I+II+ISTR)
      EE      = EE + MAX(0.,PRC(I+II+ISTR)%F2R)/ROJ
      AN(I+II)= AN(I+II) - MAX(0.,-PRC(I+II+ISTR)%F2R)/ROJ/C(I+II)**2
c      AP(I+II) = AP(I+II) + A2(I+II+ISTR)**2*ROJ*APJ

      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
c      ROJ = RO(I+II)
c      IF(PRC(I+II+ISTR)%F2R > 0.) ROJ = RO(I+II-ISTR)
      EE      = EE + MAX(0.,-PRC(I+II)%F2R)/ROJ
      AS(I+II)= AS(I+II) - MAX(0.,PRC(I+II)%F2R)/ROJ/C(I+II)**2

      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
c      ROJ = RO(I+II)
c      IF(PRC(I+II+1)%F1R < 0.) ROJ = RO(I+II+1)
      EE      = EE + MAX(0.,PRC(I+II+1)%F1R)/ROJ
      AE(I+II)= AE(I+II) - MAX(0.,-PRC(I+II+1)%F1R)/ROJ/C(I+II)**2

      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
c      ROJ = RO(I+II)
c      IF(PRC(I+II)%F1R > 0.) ROJ = RO(I+II-1)
      EE      = EE + MAX(0.,-PRC(I+II)%F1R)/ROJ
      AW(I+II)= AW(I+II) - MAX(0.,PRC(I+II)%F1R)/ROJ/C(I+II)**2
c      AP(I+II) = AP(I+II) + A1(I+II)**2*ROJ*APJ

      ROJ     = .5*(RO(I+II)  + RO(I+II+IL))
c      ROJ = RO(I+II)
c      IF(PRC(I+II+IL)%F3R < 0.) ROJ = RO(I+II+IL)
      EE      = EE + MAX(0.,PRC(I+II+IL)%F3R)/ROJ
      AT(I+II)= AT(I+II)- MAX(0.,-PRC(I+II+IL)%F3R)/ROJ/C(I+II)**2

      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
c      ROJ = RO(I+II)
c      IF(PRC(I+II)%F3R > 0.) ROJ = RO(I+II-IL)
      EE      = EE + MAX(0.,-PRC(I+II)%F3R)/ROJ
      AD(I+II)= AD(I+II) - MAX(0.,PRC(I+II)%F3R)/ROJ/C(I+II)**2

      ENDIF ! COMPCORR

      CCC     = C(I+II)
      GAMMA   = 0.

      IF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN   !IF(NPHASES > 0) THEN
c      DPMAX   = 1.*ABS(P(I+II)-PRO(I+II)%PSAT) + 100.
      DPMAX   = ABS(PRO(I+II)%DPSDT*(PRO(I+II)%TEMP(1)-PRO(I+II)%TSAT))
     +          *1.+ .1 !1.E-2

      GAMMA   = VOL(I+II)*ABS(VAR(I+II)%EVAPR(1))*(-1./PRO(I+II)%RO(1)+
     +        1./PRO(I+II)%RO(2))/DPMAX
cc       GAMMA  = VOL(I+II)*(1./PRO(I+II)%DHSDP(1)+1./PRO(I+II)%DHSDP(2))/
cc     +  (PRO(I+II)%HSAT(2) - PRO(I+II)%HSAT(1))
       IF(ABS(VAR(I+II)%EVAPR(1)) < 1.E-6) Gamma = 0.

      ENDIF ! CAVIT

      IF(TIMEL) THEN    ! Include time derivative
         AP(I+II) = AP(I+II) + VOL(I+II)/(DT*CCC**2)
      ENDIF
            
      AP(I+II) = AP(I+II)+1./(CCC**2)*(EE + VOL(I+II)/DTL(I+II)*ALFAPT)
     + + GAMMA !+ VOL(I+II)/(DT*340**2) ! In the incompressible solution?

C ... Check the diagonal dominance

      AP(I+II) = MAX(AP(I+II),-AN(I+II)-AS(I+II)-AE(I+II)-AW(I+II)
     + -AT(I+II)-AD(I+II))

      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF

1000  CONTINUE

************************************************************************

C ... Boundary cells (Von Neumann condition for the implicit stage)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+IL))
      APJ     = 1./APU(I+II+IL) ! APU(I+II) N.A. everywhere along the boundary
      ROJ     = .5*(RO(I+II)  + RO(I+II+IL))
      AT(I+II)= -A3(I+II+IL)**2*ROJ*APJ
c      AT(I+II)= -A3(I+II+IL)**2*APJ
      AP(I+II)= -AT(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AT(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      APJ     = 1./APU(I+II-IL)
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
      AD(I+II)= -A3(I+II)**2*ROJ*APJ
c      AD(I+II)= -A3(I+II)**2*APJ
      AP(I+II)= -AD(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AD(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+ISTR))
      APJ     = 1./APU(I+II+ISTR)
      ROJ     = .5*(RO(I+II)  + RO(I+II+ISTR))
      AN(I+II)= -A2(I+II+ISTR)**2*ROJ*APJ
c      AN(I+II)= -A2(I+II+ISTR)**2*APJ
      AP(I+II)= -AN(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AN(I+II) = 0. ! Dirichlet condition
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      APJ     = 1./APU(I+II-ISTR)
      ROJ     = .5*(RO(I+II)  + RO(I+II-ISTR))
      AS(I+II)= -A2(I+II)**2*ROJ*APJ
c      AS(I+II)= -A2(I+II)**2*APJ
      AP(I+II)= -AS(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AS(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2300  CONTINUE
      DO 2400 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+1))
      APJ     = 1./APU(I+II+1)
      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
      AE(I+II)= -A1(I+II+1)**2*ROJ*APJ
c      AE(I+II)= -A1(I+II+1)**2*APJ
      AP(I+II)= -AE(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2400  CONTINUE
      DO 2500 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      APJ     = 1./APU(I+II-1)
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
      AW(I+II)= -A1(I+II)**2*ROJ*APJ
c      AW(I+II)= -A1(I+II)**2*APJ
      AP(I+II)= -AW(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2500  CONTINUE
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   
      END SUBROUTINE AAMATP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATP_CA(AP,AN,AS,AE,AW,AT,AD,APU,PRC,PRO,VAR,RO,DRO,
     + VOL,P,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAP,M,
     + NGL,XFC,YFC,ZFC,PDIFF,RKSI,NCHIM,SDI,TIMEL,DT)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NCHIM
      INTEGER :: N,M,ISTR,JSTR,IL,I,J,K,II,JJ,NGL,IEVAP

      REAL :: ALFAP,ALFAPT,APJ,ROJ,CCC,YPC1,YPC2,GAMMA,DPMAX,EE,GAMMA2,
     + DPLIQ,DPGAS,HFG,DT
      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),APU(*),RO(*),VOL(*),C(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),P(*),XFC(*),YFC(*),ZFC(*),PDIFF(*),
     + RKSI(*)

      TYPE(PRE_COR) PRC(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(BLOCKS) BLKS(*)
      TYPE(DIFF_COR) SDI(*)

      LOGICAL :: COMPCORR, TIMEL

C
C ... DISCRETIZE POISSON EQUATION
C

      COMPCORR = BLKS(NGL)%COMPCORR
      IEVAP = BLKS(NGL)%EVAP_TYPE
 
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      ALFAPT = .000001*M  !Vedellä nolla pois
      ALFAPT = .01*M ! Oli 0.01*M

C ... Symmetric part of the Poisson equation

      DO K = 1,KMAX + 1
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
      AD(I+II)= -A3(I+II)**2*APJ
      AT(I+II-IL) = AD(I+II) ! uusi

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      ROJ     = .5*(RO(I+II)  + RO(I+II-ISTR))
      AS(I+II)= -A2(I+II)**2*APJ
      AN(I+II-ISTR) = AS(I+II) ! uusi

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX + 1

      APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
      AW(I+II)= -A1(I+II)**2*APJ  
      AE(I+II-1) = AW(I+II)

      ENDDO ; ENDDO ; ENDDO ! Mersu, tonttu hävisi

C ... Diagonal element is added together

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      AP(I+II) = -AN(I+II)-AS(I+II)-AE(I+II)-AW(I+II)-AT(I+II)-AD(I+II)
c      write(7667,*) AP(I+II),AT(I+II),AD(I+II),AN(I+II),AS(I+II),
c     +              AE(I+II),AW(I+II)

      ENDDO ; ENDDO ; ENDDO

C ... Compressibility correction and mass-transfer relaxation

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... Mersu, laskee kahdesti

      EE = 0.

      CCC     = C(I+II)
      GAMMA   = 0.

      IF(BLKS(NGL)%SOLUTION_TYPE == 'CAVIT') THEN

      SELECT CASE(IEVAP)

      CASE(0) ! No phase change

         GAMMA = 0.

      CASE(1,2,6) ! Standard cavitation

      DPMAX   = 1.*ABS(P(I+II)-PRO(I+II)%PSAT) + 1.E-3 !1.
c      DPMAX   = 1.*ABS(PDIFF(I+II)-PRO(I+II)%DPSAT) + 1.E-2 !1.
      GAMMA   = VOL(I+II)*ABS(VAR(I+II)%EVAPR(1))*(-1./PRO(I+II)%RO(1)+
     +        1./PRO(I+II)%RO(2))/DPMAX

      CASE(3,4,5,7) ! Temperature-based

c      DPMAX   = ABS(PRO(I+II)%DPSDT*(PRO(I+II)%TEMP(1)-PRO(I+II)%TSAT))
c     +          *1.+ .1 !1.E-2
      DPLIQ = ABS(PRO(I+II)%DPSDT*(PRO(I+II)%DTEMP(1)-PRO(I+II)%TSAT))
     +          *2.+ 1. !10.

      DPGAS =  ABS(PRO(I+II)%DPSDT*(PRO(I+II)%DTEMP(2)-PRO(I+II)%TSAT))
     +          *2. + 1. !10.  !.0001 

      GAMMA   = VOL(I+II)*(ABS(VAR(I+II)%EVAPr(1))/DPLIQ +  ! edellinen
     +  ABS(VAR(I+II)%EVAPr(2))/DPGAS) * 
     +  (-1./PRO(I+II)%RO(1) + 1./PRO(I+II)%RO(2)) ! ei ole kokeiltu

      END SELECT

      ELSEIF(BLKS(NGL)%SOLUTION_TYPE == 'MULTI') THEN
         STOP'Wrong solution option in AAMATP_CA'
      ENDIF
            
      EE = EE + VOL(I+II)/(RO(I+II)*DTL(I+II))
      AP(I+II)= -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II) -
     + AD(I+II) + 1./(CCC**2)*EE*ALFAPT !+ VOL(I+II)/(RO(I+II)*CCC**2*DT)
     + + GAMMA    ! Vesi                  This was added twice??

      IF(TIMEL) AP(I+II) = AP(I+II) + VOL(I+II)/(RO(I+II)*CCC**2*DT)

      AP(I+II)= AP(I+II)/ALFAP

c      ENDDO; ENDDO; ENDDO
            
c      AP(I+II) = AP(I+II)+1./(CCC**2)*(EE + VOL(I+II)/DTL(I+II)*ALFAPT)
c     + + GAMMA !+ VOL(I+II)/(DT*340**2) ! In the incompressible solution?

C ... Check the diagonal dominance

      AP(I+II) = MAX(AP(I+II),-AN(I+II)-AS(I+II)-AE(I+II)-AW(I+II)
     + -AT(I+II)-AD(I+II))

      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF

      ENDDO; ENDDO; ENDDO

1000  CONTINUE

************************************************************************

C ... Boundary cells (Von Neumann condition for the implicit stage)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+IL))
      APJ     = 1./APU(I+II+IL) ! APU(I+II) N.A. everywhere along the boundary
      ROJ     = .5*(RO(I+II)  + RO(I+II+IL))
      AT(I+II)= -A3(I+II+IL)**2*APJ
c      AT(I+II)= -A3(I+II+IL)**2*APJ
      AP(I+II)= -AT(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AT(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      APJ     = 1./APU(I+II-IL)
      ROJ     = .5*(RO(I+II)  + RO(I+II-IL))
      AD(I+II)= -A3(I+II)**2*APJ
c      AD(I+II)= -A3(I+II)**2*APJ
      AP(I+II)= -AD(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AD(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+ISTR))
      APJ     = 1./APU(I+II+ISTR)
      ROJ     = .5*(RO(I+II)  + RO(I+II+ISTR))
      AN(I+II)= -A2(I+II+ISTR)**2*APJ
c      AN(I+II)= -A2(I+II+ISTR)**2*APJ
      AP(I+II)= -AN(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AN(I+II) = 0. ! Dirichlet condition
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      APJ     = 1./APU(I+II-ISTR)
      ROJ     = .5*(RO(I+II)  + RO(I+II-ISTR))
      AS(I+II)= -A2(I+II)**2*APJ
c      AS(I+II)= -A2(I+II)**2*APJ
      AP(I+II)= -AS(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AS(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2300  CONTINUE
      DO 2400 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+1))
      APJ     = 1./APU(I+II+1)
      ROJ     = .5*(RO(I+II)  + RO(I+II+1))
      AE(I+II)= -A1(I+II+1)**2*APJ
c      AE(I+II)= -A1(I+II+1)**2*APJ
      AP(I+II)= -AE(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2400  CONTINUE
      DO 2500 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      APJ     = 1./APU(I+II-1)
      ROJ     = .5*(RO(I+II)  + RO(I+II-1))
      AW(I+II)= -A1(I+II)**2*APJ
c      AW(I+II)= -A1(I+II)**2*APJ
      AP(I+II)= -AW(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.
      IF(NCHIM > 0) THEN
         AP(I+II) = AP(I+II) + RKSI(I+II)
      ENDIF
2500  CONTINUE
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
      END SUBROUTINE AAMATP_CA
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
       SUBROUTINE AAMATP_MF(AP,AN,AS,AE,AW,AT,AD,APU,PRC,PRO,VAR,RO,DRO,
     + VOL,P,PDIFF,C,A1,A2,A3,DTL,BLKS,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     + ALFAP,M,NGL,TIMEL,DT,APUULG,NPHASE,AP1X,AP1Y,AP1Z,AP2X,AP2Y,AP2Z,
     + AP3X,AP3Y,AP3Z,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,ALF1LG,
     + ALF2LG,ALF3LG,AP1,AP2,AP3)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : EPS6

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NPHASE,IPHASE
      INTEGER N,M,ISTR,JSTR,IL,I,J,K,II,JJ,NGL,IEVAP

      REAL ALFAP,ALFAPT,APJ,ROJ,CCC,YPC1,YPC2,GAMMA,DPMAX,EE,GAMMA2,
     + DPLIQ,DPGAS,HFG,DT,AAVE
      REAL NXI,NYI,NZI,NXIP1,NYIP1,NZIP1 ! Antakee mersu
      REAL AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),APU(*),RO(*),VOL(*),C(*),DTL(*),
     + A1(*),A2(*),A3(*),DRO(*),P(*),PDIFF(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)
      REAL APUULG(NTOT,NPHASE),AP1X(NTOT,NPHASE),AP1Y(NTOT,NPHASE),
     + AP1Z(NTOT,NPHASE),AP2X(NTOT,NPHASE),AP2Y(NTOT,NPHASE),
     + AP2Z(NTOT,NPHASE),AP3X(NTOT,NPHASE),AP3Y(NTOT,NPHASE),
     + AP3Z(NTOT,NPHASE),ALF1LG(NTOT,NPHASE),ALF2LG(NTOT,NPHASE),
     + ALF3LG(NTOT,NPHASE),AP1(NTOT,NPHASE),AP2(NTOT,NPHASE),
     + AP3(NTOT,NPHASE)
  
      TYPE(PRE_COR) PRC(*)
      TYPE(PROPERTIES) PRO(*)
      TYPE(MPHASE_VARIABLES) VAR(*)
      TYPE(BLOCKS) BLKS(*)
      LOGICAL TIMEL

C
C ... DISCRETIZE POISSON EQUATION. No Chimera approach yet.
C

      IEVAP = BLKS(NGL)%EVAP_TYPE
     
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AW(N) = 0.
      AE(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) ! In this solution system 
      ENDDO ! Antakee mersu

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

c      ALFAPT = .000001*M  !Vedellä nolla pois
     	ALFAPT = 0.1*M ! oli 0.01

      DO IPHASE = 1,NPHASE

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX
c        write(6990,*) i,j,apu(i+ii)
c        write(6991,*) i,j,apuulg(i+ii,1)
c        write(6992,*) i,j,apuulg(i+ii,2)
C ... Mersu, laskee kahdesti
       EE = 0.

C ap2
      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II+ISTR,IPHASE)/AP1(I+II+ISTR,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP2X(I+II,IPHASE) +
c     +         VAR(I+II+ISTR)%ALFA(IPHASE)/AP2X(I+II+ISTR,IPHASE))
c     +        +  NYIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP2Y(I+II,IPHASE) +
c     +         VAR(I+II+ISTR)%ALFA(IPHASE)/AP2Y(I+II+ISTR,IPHASE))
c     +        +  NZIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP2Z(I+II,IPHASE) +
c     +         VAR(I+II+ISTR)%ALFA(IPHASE)/AP2Z(I+II+ISTR,IPHASE)))
c      AAVE  = .5*(VAR(I+II)%ALFA(IPHASE) + VAR(I+II+ISTR)%ALFA(IPHASE))
      AN(I+II)=-A2(I+II+ISTR)*VAR(I+II+ISTR)%F2A(IPHASE)*APJ + AN(I+II)

      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II-ISTR,IPHASE)/AP1(I+II-ISTR,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXI**2*(VAR(I+II)%ALFA(IPHASE)/AP2X(I+II,IPHASE) +
c     +         VAR(I+II-ISTR)%ALFA(IPHASE)/AP2X(I+II-ISTR,IPHASE))
c     +        +  NYI**2*(VAR(I+II)%ALFA(IPHASE)/AP2Y(I+II,IPHASE) +
c     +         VAR(I+II-ISTR)%ALFA(IPHASE)/AP2Y(I+II-ISTR,IPHASE))
c     +        +  NZI**2*(VAR(I+II)%ALFA(IPHASE)/AP2Z(I+II,IPHASE) +
c     +         VAR(I+II-ISTR)%ALFA(IPHASE)/AP2Z(I+II-ISTR,IPHASE)))
      AS(I+II)= -A2(I+II)*VAR(I+II)%F2A(IPHASE)*APJ + AS(I+II)

c ap1
      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II+1,IPHASE)/AP1(I+II+1,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP1X(I+II,IPHASE) +
c     +         VAR(I+II+1)%ALFA(IPHASE)/AP2X(I+II+1,IPHASE))
c     +        +  NYIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP1Y(I+II,IPHASE) +
c     +         VAR(I+II+1)%ALFA(IPHASE)/AP2Y(I+II+1,IPHASE))
c     +        +  NZIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP1Z(I+II,IPHASE) +
c     +         VAR(I+II+1)%ALFA(IPHASE)/AP2Z(I+II+1,IPHASE)))
c      AAVE  = .5*(VAR(I+II)%ALFA(IPHASE) + VAR(I+II+1)%ALFA(IPHASE))
      AE(I+II)= -A1(I+II+1)*VAR(I+II+1)%F1A(IPHASE)*APJ + AE(I+II)
c       write(7060+ngl,*) i,j,VAR(I+II+1)%F1A(IPHASE),
c     2 A1(I+II+1),APJ,VAR(I+II)%ALFA(IPHASE)/APUULG(I+II,IPHASE),
c     3 VAR(I+II+1)%ALFA(IPHASE)/APUULG(I+II+1,IPHASE),AE(I+II),apj

      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II-1,IPHASE)/AP1(I+II-1,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXI**2*(VAR(I+II)%ALFA(IPHASE)/AP1X(I+II,IPHASE) +
c     +         VAR(I+II-1)%ALFA(IPHASE)/AP2X(I+II-1,IPHASE))
c     +        +  NYI**2*(VAR(I+II)%ALFA(IPHASE)/AP1Y(I+II,IPHASE) +
c     +         VAR(I+II-1)%ALFA(IPHASE)/AP2Y(I+II-1,IPHASE))
c     +        +  NZI**2*(VAR(I+II)%ALFA(IPHASE)/AP1Z(I+II,IPHASE) +
c     +         VAR(I+II-1)%ALFA(IPHASE)/AP2Z(I+II-1,IPHASE)))
      AW(I+II)= -A1(I+II)*VAR(I+II)%F1A(IPHASE)*APJ + AW(I+II)
c       write(7070+ngl,*) i,j,VAR(I+II)%F1A(IPHASE),
c     2 A1(I+II),APJ,VAR(I+II)%ALFA(IPHASE)/APUULG(I+II,IPHASE),
c     3 VAR(I+II-1)%ALFA(IPHASE)/APUULG(I+II-1,IPHASE),AW(I+II),apj

c ap3
      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II+IL,IPHASE)/AP1(I+II+IL,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3X(I+II,IPHASE) +
c     +         VAR(I+II+IL)%ALFA(IPHASE)/AP2X(I+II+IL,IPHASE))
c     +        +  NYIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3Y(I+II,IPHASE) +
c     +         VAR(I+II+IL)%ALFA(IPHASE)/AP2Y(I+II+IL,IPHASE))
c     +        +  NZIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3Z(I+II,IPHASE) +
c     +         VAR(I+II+IL)%ALFA(IPHASE)/AP2Z(I+II+IL,IPHASE)))
      AT(I+II)= -A3(I+II+IL)*VAR(I+II+IL)%F3A(IPHASE)*APJ + AT(I+II)

      APJ   =.5*(ALF1LG(I+II,IPHASE)/AP1(I+II,IPHASE) +
     +         ALF1LG(I+II-IL,IPHASE)/AP1(I+II-IL,IPHASE))
     +         /PRO(I+II)%RO(IPHASE)

c      APJ   =.5*(NXIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3X(I+II,IPHASE) +
c     +         VAR(I+II-IL)%ALFA(IPHASE)/AP2X(I+II-IL,IPHASE))
c     +        +  NYIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3Y(I+II,IPHASE) +
c     +         VAR(I+II-IL)%ALFA(IPHASE)/AP2Y(I+II-IL,IPHASE))
c     +        +  NZIP1**2*(VAR(I+II)%ALFA(IPHASE)/AP3Z(I+II,IPHASE) +
c     +         VAR(I+II-IL)%ALFA(IPHASE)/AP2Z(I+II-IL,IPHASE)))
      AD(I+II)= -A3(I+II)*VAR(I+II)%F3A(IPHASE)*APJ + AD(I+II)
c       write(7000+ngl,*) i,j,ae(i+ii),aw(i+ii),an(i+ii),as(i+ii),
c     2 at(i+ii),ad(i+ii)

      ENDDO; ENDDO; ENDDO ! Voi mersu

      ENDDO ! IPHASE

c ... test
c          stop 'aamatp'
      CCC     = C(I+II)
      GAMMA   = 0.
      DPMAX   = 1000.
      ee      = 0.


      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

      CCC     = C(I+II)
      ee      = 0.


      SELECT CASE(IEVAP)

      CASE(0) ! No phase change

         GAMMA = 0.

      CASE(1,2,6) ! Standard cavitation

      DPMAX   = 1.*ABS(P(I+II)-PRO(I+II)%PSAT) + 1.E-3 !1.
c      DPMAX   = 1.*ABS(PDIFF(I+II)-PRO(I+II)%DPSAT) + 1.E-2 !1.
      GAMMA   = VOL(I+II)*ABS(VAR(I+II)%EVAPR(1))*(-1./PRO(I+II)%RO(1)+
     +        1./PRO(I+II)%RO(2))/DPMAX

      CASE(3,4,5,7) ! Temperature-based

c      DPMAX   = ABS(PRO(I+II)%DPSDT*(PRO(I+II)%TEMP(1)-PRO(I+II)%TSAT))
c     +          *1.+ .1 !1.E-2
      DPLIQ = ABS(PRO(I+II)%DPSDT*(PRO(I+II)%DTEMP(1)-PRO(I+II)%TSAT))
     +          *2.+ 1. !10.

      DPGAS =  ABS(PRO(I+II)%DPSDT*(PRO(I+II)%DTEMP(2)-PRO(I+II)%TSAT))
     +          *2. + 1. !10.  !.0001 

      GAMMA   = VOL(I+II)*(ABS(VAR(I+II)%EVAPr(1))/DPLIQ +  ! edellinen
     +  ABS(VAR(I+II)%EVAPr(2))/DPGAS) * 
     +  (-1./PRO(I+II)%RO(1) + 1./PRO(I+II)%RO(2)) ! ei ole kokeiltu

      END SELECT

      EE = EE + VOL(I+II)/(RO(I+II)*DTL(I+II))
      AP(I+II)= -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II) -
     + AD(I+II) + 1./(CCC**2)*EE*ALFAPT !+ VOL(I+II)/(RO(I+II)*CCC**2*DT)
     + + GAMMA    ! Vesi                  This was added twice??

      IF(TIMEL) AP(I+II) = AP(I+II) + VOL(I+II)/(RO(I+II)*CCC**2*DT)
c       write(7010+ngl,*) i,j,ae(i+ii),aw(i+ii),an(i+ii),as(i+ii),
c     2 at(i+ii),ad(i+ii)

      AP(I+II)= AP(I+II)/ALFAP

1000  CONTINUE
         
C ... Boundary cells (Von Neumann condition for the implicit stage)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+IL))
      APJ     = 1./APU(I+II+IL) ! APU(I+II) N.A. everywhere along the boundary
c      AT(I+II)= -A3(I+II+IL)**2*ROJ*APJ
      AT(I+II)= -A3(I+II+IL)**2*APJ
      AP(I+II)= AP(I+II+IL) !-AT(I+II)!+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AT(I+II) = 0.
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-IL))
      APJ     = 1./APU(I+II-IL)
c      AD(I+II)= -A3(I+II)**2*ROJ*APJ
      AD(I+II)= -A3(I+II)**2*APJ
      AP(I+II)= AP(I+II-IL) !-AD(I+II)!+ C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AD(I+II) = 0.
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+ISTR))
      APJ     = 1./APU(I+II+ISTR)
c      AN(I+II)= -A2(I+II+ISTR)**2*ROJ*APJ
      AN(I+II)= -A2(I+II+ISTR)**2*APJ
      AP(I+II)= AP(I+II+ISTR)!-AN(I+II)!+C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AN(I+II) = 0. ! Dirichlet condition
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-ISTR))
      APJ     = 1./APU(I+II-ISTR)
c      AS(I+II)= -A2(I+II)**2*ROJ*APJ
      AS(I+II)= -A2(I+II)**2*APJ
      AP(I+II)= AP(I+II-ISTR)!-AS(I+II)!+C(I+II)**2*VOL(I+II)/DTL(I+II)*.00001
      AS(I+II) = 0.
2300  CONTINUE
      DO 2400 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2400 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II+1))
      APJ     = 1./APU(I+II+1)
c      AE(I+II)= -A1(I+II+1)**2*ROJ*APJ
      AE(I+II)= -A1(I+II+1)**2*APJ
      AP(I+II)= AP(I+II+1)!-AE(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AE(I+II) = 0.
2400  CONTINUE
      DO 2500 K = 0,KMAX+1!1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2500 J = 0,JMAX+1!1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
c     APJ     = .5*(1./APU(I+II) + 1./APU(I+II-1))
      APJ     = 1./APU(I+II-1)
c      AW(I+II)= -A1(I+II)**2*ROJ*APJ
      AW(I+II)= -A1(I+II)**2*APJ
      AP(I+II)= AP(I+II-1)!-AW(I+II) !+ C(I+II)**2*VOL(I+II)/DTL(I+II)
      AW(I+II) = 0.
2500  CONTINUE
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
      END SUBROUTINE AAMATP_MF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C       
      SUBROUTINE AAMATP_WALLS(AP,AN,AS,AE,AW,AT,AD,DRO,DRO1,NPATCH,ICON,
     + IHF,XC,YC,ZC,XCO,YCO,ZCO,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     + IMAX,JMAX,KMAX,NTOT,INRE,JNRE,KNRE,M,NGL)

      USE NS3CO, ONLY : LN, IC9
      
      IMPLICIT NONE

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     + NTOT,ISTRID,JSTRID,KSTRID,
     + IL,IP,IBC,IDIR,IFACE,KX1,KX2,KY1,KY2,ISTR,JSTR,KSTR,I,J,K,ID,
     + LSTIRD,NFL,KA,KBIND,IJ,L,LSTRID,II,II1,M,IPL,NR,NP,L1,L2,L3,L4,
     + NGL

      INTEGER :: ICON(IC9,*), IHF(*)

      REAL :: ALFA,SXI,SYI,SZI,SLEN,XCPNP,YCPNP,ZCPNP

      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),
     + AD(*),DRO(*),DRO1(*),
     + A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),A2Z(*),
     + A3X(*),A3Y(*),A3Z(*)

      REAL :: XC(*),YC(*),ZC(*),XCO(*),YCO(*),ZCO(*)

      LOGICAL :: MIRL

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      IF(M == 1) THEN
         MIRL = .TRUE. 
      ELSE
         MIRL = .FALSE. ! Avoid mirror at lower levels
      ENDIF
      
C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH
C ... A patch type and indicator for solids
         IBC     = ICON(1,IP)
           
         IF(IBC >= 7.AND.IBC <= 10  .OR.  IBC == 15 !) then ! FRE away 8.11.2010
     +    .OR. IBC == 4) THEN ! A solid type patch, !MIR away 11.11.2008
                                                    !INL away 21.10.2011
         IPL     = ICON(2,IP) ! Proces local patch number
c        IPG     = ICON(25,IP)! Global patch number
         IDIR    = 1
         IFACE   = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1

C ... Patch indeces

         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
             
C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K  = IMAX + 1
C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K  = JMAX +1  
C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K  = KMAX +1
         IF(M > 1 .AND. KMAX == 1) MIRL = .TRUE. ! For 2D cases
      ENDIF

      IF(IFACE == 2 .OR. IFACE == 5) THEN
         KY1  = ICON(4,IP)
         KY2  = ICON(5,IP)
         KX1  = ICON(6,IP)
         KX2  = ICON(7,IP)
      ENDIF

      LSTRID  = KSTR*IDIR
      KA      = (KN+K-1)*KSTR
      KBIND   = 0
      IF(IFACE >= 4) KBIND = -KSTR

C ... Calculate coefficients for the Poisson equation on the walls 
C ... Von Neumann condition is applied

      IF(IBC /=7 .AND. IBC /=4 .AND. MIRL) THEN ! Solid type patch

      DO J = KY1,KY2
          IJ      = (JN+J-1)*JSTR + KA
          DO I = KX1,KX2
            II1   = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II    = II1 - LSTRID           ! UNDER THE SURFACE
            IF(IFACE == 1) THEN
            AE(II) = -AP(II)
            ELSE IF(IFACE == 2) THEN
            AN(II) = -AP(II)
            ELSE IF(IFACE == 3) THEN
            AT(II) = -AP(II)
            ELSE IF(IFACE == 4) THEN
            AW(II) = -AP(II)
            ELSE IF(IFACE == 5) THEN
            AS(II) = -AP(II)
            ELSE IF(IFACE == 6) THEN
            AD(II) = -AP(II)
            ENDIF
            DRO(II) = 0.
            DRO1(II)= 0.
          ENDDO
      ENDDO

      ELSE IF(IBC == 7.AND. MIRL) THEN ! Singularity (with mirror?)

      DO J = KY1,KY2
          IJ      = (JN+J-1)*JSTR + KA 
          DO I = KX1,KX2
            II1   = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II    = II1 - LSTRID           ! UNDER THE SURFACE
            IF(IFACE == 1) THEN
            AP(II) = AP(II1)
            AE(II) = -AP(II)
            ELSE IF(IFACE == 2) THEN
            AP(II) = AP(II1)
            AN(II) = -AP(II)
            ELSE IF(IFACE == 3) THEN
            AP(II) = AP(II1)
            AT(II) = -AP(II)
            ELSE IF(IFACE == 4) THEN
            AP(II) = AP(II1)
            AW(II) = -AP(II)
            ELSE IF(IFACE == 5) THEN
            AP(II) = AP(II1)
            AS(II) = -AP(II)
            ELSE IF(IFACE == 6) THEN
            AP(II) = AP(II1)
            AD(II) = -AP(II)
            ENDIF
            DRO(II) = 0.
            DRO1(II)= 0.
          ENDDO
      ENDDO

      ELSE IF(IBC == 4 .AND. MIRL) THEN         ! Mirror added 9.4.2010

      DO J = KY1,KY2
          IJ      = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
          DO I = KX1,KX2
            L1      = 1 + (IN+I-1)*ISTR + IJ  ! CELL INDEX UP (STRIDE=ISTR)
            L2      = 1 + (IN+I  )*ISTR + IJ
            L3      = 1 + (IN+I-1)*ISTR + IJ + JSTR
            L4      = 1 + (IN+I  )*ISTR + IJ + JSTR

            II1   = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II    = II1 - LSTRID          ! UNDER THE SURFACE
            NP  = I + NR                     ! Patch index with LN ghost cells
            ZCPNP = .25*(ZCO(L1) + ZCO(L2) + ZCO(L3) + ZCO(L4)) 
            YCPNP = .25*(YCO(L1) + YCO(L2) + YCO(L3) + YCO(L4))
            XCPNP = .25*(XCO(L1) + XCO(L2) + XCO(L3) + XCO(L4))

            SXI    = XCPNP - XC(II1)
            SYI    = YCPNP - YC(II1)
            SZI    = ZCPNP - ZC(II1)
            SLEN   = SQRT(SXI**2 + SYI**2 + SZI**2) + 1.E-21
            SXI    = SXI/SLEN
            SYI    = SYI/SLEN
            SZI    = SZI/SLEN

            IF(IFACE == 1) THEN
c            AP(II1) = AP(II1) + AW(II1)  ! Wigley hull diverges
c            AW(II1)= 0.
            ALFA   = SXI*A1X(II1) + SYI*A1Y(II1) + SZI*A1Z(II1)
            AP(II) = AP(II1)
            AE(II) = -AP(II) * ABS(ALFA)
            ELSE IF(IFACE == 2) THEN
            ALFA   = SXI*A2X(II1) + SYI*A2Y(II1) + SZI*A2Z(II1)
c            AS(II1)= 0.
            AP(II) = AP(II1)
            AN(II) = -AP(II) * ABS(ALFA)
            ELSE IF(IFACE == 3) THEN
c           AD(II1)= 0.
            ALFA   = SXI*A3X(II1) + SYI*A3Y(II1) + SZI*A3Z(II1)
            AP(II) = AP(II1)
            AT(II) = -AP(II) * ABS(ALFA)
            ELSE IF(IFACE == 4) THEN
c            AE(II1)= 0.
            ALFA   = SXI*A1X(II) + SYI*A1Y(II) + SZI*A1Z(II)
            AP(II) = AP(II1)
            AW(II) = -AP(II) * ABS(ALFA)
            ELSE IF(IFACE == 5) THEN
c            AN(II1)= 0.
            ALFA   = SXI*A2X(II) + SYI*A2Y(II) + SZI*A2Z(II)
            AP(II) = AP(II1)
            AS(II) = -AP(II) * ABS(ALFA)
            ELSE IF(IFACE == 6) THEN
            ALFA   = SXI*A3X(II) + SYI*A3Y(II) + SZI*A3Z(II)
cc            AP(II1)= AP(II1) + AT(II1)
cc            AT(II1)= 0.
            AP(II) = AP(II1)
            AD(II) = -AP(II) * ABS(ALFA)
            ENDIF
            DRO(II) = 0.
            DRO1(II)= 0.
          ENDDO
      ENDDO
      ENDIF

      ENDIF
8000  CONTINUE

      RETURN
      END SUBROUTINE AAMATP_WALLS
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE MIRP(PDIF,NPATCH,ICON,IMAX,JMAX,KMAX,NTOT,
     +                INRE,JNRE,KNRE,NGL)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: ICON(IC9,*)

      INTEGER :: NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     + NTOT,ISTRID,JSTRID,KSTRID,NGL,
     + IL,IP,IBC,IDIR,IFACE,KX1,KX2,KY1,KY2,ISTR,JSTR,KSTR,I,J,K,ID,
     + LSTIRD,NFL,KA,KBIND,IJ,L,LSTRID,II,II1

      REAL :: PDIF(NTOT)

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID
      
C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH
C ... A patch type and indicator for solids
         IBC     = ICON(1,IP)
           
        IF(IBC >= 7.AND.IBC <= 10  .OR.  IBC == 13 .OR. IBC == 15 !) testing
c        IF(IBC >= 7.AND.IBC <= 10 .OR. IBC == 15 !) ! Without a free surface???
c     &  .OR. IBC == 4 .OR. IBC == 3 .OR. IBC == 5) THEN ! Mirror the pressure. Inlet???(6.7.11)
C ... Mirroring the outlet is extremely harmful, inlet (IBC=3) must be studied
     &  .OR. IBC == 4 .OR. IBC == 3) THEN ! Mirror the pressure.

         IDIR    = 1
         IFACE   = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1
        
C ... Patch indeces

         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
             
C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K  = IMAX + 1
C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K  = JMAX +1  
C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K  = KMAX +1 
      ENDIF

      IF(IFACE ==  2 .OR. IFACE == 5) THEN
         KY1  = ICON(4,IP)
         KY2  = ICON(5,IP)
         KX1  = ICON(6,IP)
         KX2  = ICON(7,IP)
      ENDIF

      LSTRID  = KSTR*IDIR
      KA      = (KN+K-1)*KSTR
      KBIND   = 0
      IF(IFACE >= 4) KBIND = - KSTR

C ... Von Neumann condition is applied

      DO J = KY1,KY2
          IJ      = (JN+J-1)*JSTR + KA 
          DO I = KX1,KX2
            II1   = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II    = II1 - LSTRID           ! UNDER THE SURFACE
            PDIF(II) = PDIF(II1)
          ENDDO
      ENDDO

      ENDIF
8000  CONTINUE

      RETURN
      END SUBROUTINE MIRP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESCOR(Q,QC,QM,DP,DPDX,DPDY,DPDZ,AP,AN,AS,AE,AW,AT,AD,
     + XFC,YFC,ZFC,SDI,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,BLKS,IMAX,
     + JMAX,KMAX,IN,JN,KN,NGL,MULT,AMGDIVL,PRINTOUTL)

C ... Correct the pressure gradient 'dot' normal vector (laskee kahdesti)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,I,J,K,
     +           ISTR,JSTR,IL,II,JJ,ngl,mult

      REAL :: SNA1X,SNA1Y,SNA1Z,SNA2X,SNA2Y,SNA2Z,SNA3X,SNA3Y,SNA3Z,
     + DPDE,DPDW,DPDN,DPDS,DPDT,DPDD,ALFAS,QAPU,QAPU2,RELAX

      REAL :: Q(*),QC(*),QM(*),DPDX(*),DPDY(*),DPDZ(*),
     + AP(*),AN(*),AS(*),
     + AE(*),AW(*),AT(*),AD(*),XFC(*),YFC(*),ZFC(*),DP(*),A1X(*),A1Y(*),
     + A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      TYPE(DIFF_COR) :: SDI(*)
      TYPE(BLOCKS)   :: BLKS(*)

      LOGICAL :: AMGDIVL, PRINTOUTL

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      ALFAS = 1.  !1.0

      RELAX = MIN(BLKS(NGL)%SKEWIMIN,BLKS(NGL)%SKEWJMIN,
     +            BLKS(NGL)%SKEWKMIN)
      RELAX = MIN(MAX(RELAX,0.06),BLKS(NGL)%ALFASKEW)

      QAPU2 = 0.

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      SNA1X = A1X(I+II) - XFC(I+II)*SDI(I+II)%S1X
      SNA1Y = A1Y(I+II) - XFC(I+II)*SDI(I+II)%S1Y
      SNA1Z = A1Z(I+II) - XFC(I+II)*SDI(I+II)%S1Z

      SNA2X = A2X(I+II) - YFC(I+II)*SDI(I+II)%S2X
      SNA2Y = A2Y(I+II) - YFC(I+II)*SDI(I+II)%S2Y
      SNA2Z = A2Z(I+II) - YFC(I+II)*SDI(I+II)%S2Z

      SNA3X = A3X(I+II) - ZFC(I+II)*SDI(I+II)%S3X
      SNA3Y = A3Y(I+II) - ZFC(I+II)*SDI(I+II)%S3Y
      SNA3Z = A3Z(I+II) - ZFC(I+II)*SDI(I+II)%S3Z

      DPDW  = .5*(SNA1X*(DPDX(I+II-1)+DPDX(I+II)) +
     +            SNA1Y*(DPDY(I+II-1)+DPDY(I+II)) +
     +            SNA1Z*(DPDZ(I+II-1)+DPDZ(I+II)))
      DPDS  = .5*(SNA2X*(DPDX(I+II-ISTR)+DPDX(I+II)) +
     +            SNA2Y*(DPDY(I+II-ISTR)+DPDY(I+II)) +
     +            SNA2Z*(DPDZ(I+II-ISTR)+DPDZ(I+II)))
      DPDD  = .5*(SNA3X*(DPDX(I+II-IL)+DPDX(I+II)) +
     +            SNA3Y*(DPDY(I+II-IL)+DPDY(I+II)) +
     +            SNA3Z*(DPDZ(I+II-IL)+DPDZ(I+II)))

      if(i == 1) dpdw = 0.
      if(j == 1) dpds = 0.
      if(k == 1) dpdd = 0.
C ... uudestaan, korjattava
      SNA1X = A1X(I+II+1) - XFC(I+II+1)*SDI(I+II+1)%S1X
      SNA1Y = A1Y(I+II+1) - XFC(I+II+1)*SDI(I+II+1)%S1Y
      SNA1Z = A1Z(I+II+1) - XFC(I+II+1)*SDI(I+II+1)%S1Z

      SNA2X = A2X(I+II+ISTR) - YFC(I+II+ISTR)*SDI(I+II+ISTR)%S2X
      SNA2Y = A2Y(I+II+ISTR) - YFC(I+II+ISTR)*SDI(I+II+ISTR)%S2Y
      SNA2Z = A2Z(I+II+ISTR) - YFC(I+II+ISTR)*SDI(I+II+ISTR)%S2Z


      SNA3X = A3X(I+II+IL) - ZFC(I+II+IL)*SDI(I+II+IL)%S3X
      SNA3Y = A3Y(I+II+IL) - ZFC(I+II+IL)*SDI(I+II+IL)%S3Y
      SNA3Z = A3Z(I+II+IL) - ZFC(I+II+IL)*SDI(I+II+IL)%S3Z

      DPDE  = .5*(SNA1X*(DPDX(I+II+1)+DPDX(I+II)) +
     +            SNA1Y*(DPDY(I+II+1)+DPDY(I+II)) +
     +            SNA1Z*(DPDZ(I+II+1)+DPDZ(I+II)))
      DPDN  = .5*(SNA2X*(DPDX(I+II+ISTR)+DPDX(I+II)) +
     +            SNA2Y*(DPDY(I+II+ISTR)+DPDY(I+II)) +
     +            SNA2Z*(DPDZ(I+II+ISTR)+DPDZ(I+II)))
      DPDT  = .5*(SNA3X*(DPDX(I+II+IL)+DPDX(I+II)) +
     +            SNA3Y*(DPDY(I+II+IL)+DPDY(I+II)) +
     +            SNA3Z*(DPDZ(I+II+IL)+DPDZ(I+II)))

      if(i == imax) dpde = 0.
      if(j == jmax) dpdn = 0.
      if(k == kmax) dpdt = 0.

      DPDE = ALFAS*DPDE
      DPDW = ALFAS*DPDW
      DPDN = ALFAS*DPDN
      DPDS = ALFAS*DPDS
      DPDT = ALFAS*DPDT
      DPDD = ALFAS*DPDD
 
      DPDE = ALFAS*DPDE
      DPDW = ALFAS*DPDW
c      if(ngl <= 0)then
c      dpdn= 0.
c      dpds=0.
c      dpdt=0.
c      dpdd=0.
c       else
c      DPDN = ALFAS*DPDN
c      DPDS = ALFAS*DPDS
c      DPDT = ALFAS*DPDT
c      DPDD = ALFAS*DPDD
c        endif

      QAPU = ABS(QC(I+II))

      QC(I+II) = (1.-RELAX)*QC(I+II) + RELAX*(               ! Malliskaalassa .3/.7
     + + AE(I+II)*(-(1.-XFC(I+II+1))   *(DP(I+II+1)- DP(I+II))
     + + SDI(I+II+1)%S1 *DPDE)
     + + AW(I+II)*(-(1.-XFC(I+II))     *(DP(I+II-1)- DP(I+II))
     + - SDI(I+II)%S1   *DPDW)
     + + AN(I+II)*(-(1.-YFC(I+II+ISTR))*(DP(I+II+ISTR)- DP(I+II))
     + + SDI(I+II+ISTR)%S2*DPDN)
     + + AS(I+II)*(-(1.-YFC(I+II))     *(DP(I+II-ISTR)- DP(I+II))
     + - SDI(I+II)%S2   *DPDS)
     + + AT(I+II)*(-(1.-ZFC(I+II+IL))  *(DP(I+II+IL)- DP(I+II))
     + + SDI(I+II+IL)%S3*DPDT)
     + + AD(I+II)*(-(1.-ZFC(I+II))     *(DP(I+II-IL)- DP(I+II))
     + - SDI(I+II)%S3   *DPDD))

      QAPU2   = QAPU2 + QC(I+II)

      Q(I+II) = QM(I+II) + QC(I+II)
     
      ENDDO ; ENDDO ; ENDDO

      IF(ABS(QAPU2) > 1.E10) THEN
         AMGDIVL = .TRUE.
      IF(PRINTOUTL) THEN
         WRITE(13,'(A,/A,I4,2A,F5.3,/A)')
     +   '  AMG for pressure diverges:',
     +   '  For block no.',NGL,' a skewness ',
     +   'correction in RESCOR diverges with ALFASKEW = ',RELAX,
     +   '  Give a smaller value for this.'
         PRINTOUTL = .FALSE.
      ENDIF
      ENDIF ! QAPU2

      RETURN 

      END SUBROUTINE RESCOR
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRM(DRO,DU,DV,DW,RO,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,PRC,
     + VOL,NGL)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL :: ALFAU,RON,ROS,ROE,ROW,ROT,ROD,
     +        DUN,DVN,DWN,DUS,DVS,DWS,DUE,DVE,DWE,DUW,DVW,DWW,
     +        DUT,DVT,DWT,DUD,DVD,DWD
      REAL :: DRO(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),RO(*),
     + A1(*),A2(*),A3(*),VOL(*)
      TYPE(PRE_COR) :: PRC(*)

C
C ... Perform the velocity correction
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      RON     = .5*(RO(I+II) + RO(I+II+ISTR))
      ROS     = .5*(RO(I+II) + RO(I+II-ISTR))
      ROE     = .5*(RO(I+II) + RO(I+II+1))
      ROW     = .5*(RO(I+II) + RO(I+II-1))
      ROT     = .5*(RO(I+II) + RO(I+II+IL))
      ROD     = .5*(RO(I+II) + RO(I+II-IL))
c      IF(J == 1 .OR. I == 1 .OR.
c     +   J == JMAX .OR. I == IMAX) THEN
c      APN     = APU(I+II)
c      APS     = APU(I+II)
c      APE     = APU(I+II)
c      APW     = APU(I+II)
c      APT     = APU(I+II)
c      APD     = APU(I+II) 
c      ENDIF
      DRO(I+II)=  DRO(I+II) + 1.*(
     + - A2(I+II+ISTR)*A2X(I+II+ISTR)*RON*.5*(DU(I+II+ISTR)+ DU(I+II))
     + + A2(I+II)     *A2X(I+II)*ROS*     .5*(DU(I+II)+ DU(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1)*ROE*   .5*(DU(I+II+1)+ DU(I+II))
     + + A1(I+II)     *A1X(I+II)*ROW*     .5*(DU(I+II)+ DU(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL)*ROT*  .5*(DU(I+II+IL)+ DU(I+II))
     + + A3(I+II)     *A3X(I+II)*ROD*     .5*(DU(I+II)+ DU(I+II-IL))
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR)*RON*.5*(DV(I+II+ISTR)+ DV(I+II))
     + + A2(I+II)     *A2Y(I+II)*ROS*     .5*(DV(I+II)+ DV(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1)*ROE*   .5*(DV(I+II+1)+ DV(I+II))
     + + A1(I+II)     *A1Y(I+II)*ROW*     .5*(DV(I+II)+ DV(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL)*ROT*  .5*(DV(I+II+IL)+ DV(I+II))
     + + A3(I+II)     *A3Y(I+II)*ROD*     .5*(DV(I+II)+ DV(I+II-IL))
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR)*RON*.5*(DW(I+II+ISTR)+ DW(I+II))
     + + A2(I+II)     *A2Z(I+II)*ROS*     .5*(DW(I+II)+ DW(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1)*ROE*   .5*(DW(I+II+1)+ DW(I+II))
     + + A1(I+II)     *A1Z(I+II)*ROW*     .5*(DW(I+II)+ DW(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL)*ROT*  .5*(DW(I+II+IL)+ DW(I+II))
     + + A3(I+II)     *A3Z(I+II)*ROD*     .5*(DW(I+II)+ DW(I+II-IL)))

1000  CONTINUE

c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(DIFF,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,DRO,IMAX,JMAX,KMAX,IN,JN,KN,630)
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   

      RETURN
      END SUBROUTINE CORRM
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRM2(DRO,DU,DV,DW,RO,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,DP,APU,PRC,
     + DROB,VOL,NGL,ICERO)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL,ICERO
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL :: ALFAU,RON,ROS,ROE,ROW,ROT,ROD
      REAL :: DPA,DUN,DVN,DWN,DUS,DVS,DWS,DUE,DVE,DWE,DUW,DVW,DWW,
     +     DUT,DVT,DWT,DUD,DVD,DWD,APN,APS,APE,APW,APT,APD,DROA
      REAL :: DRO(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),RO(*),
     + A1(*),A2(*),A3(*),DP(*),APU(*),DROB(*),VOL(*)
      TYPE(PRE_COR) :: PRC(*)
C
C ... Perform the velocity correction for the face mass fluxes
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      IF(ICERO == 1) THEN
         PRC(1:NTOT)%FIR  = 0.
         PRC(1:NTOT)%FJR  = 0.
         PRC(1:NTOT)%FKR  = 0.
      ELSEIF(ICERO >= 2) THEN
         PRC(1:NTOT)%F1R  = PRC(1:NTOT)%FIR
         PRC(1:NTOT)%F2R  = PRC(1:NTOT)%FJR
         PRC(1:NTOT)%F3R  = PRC(1:NTOT)%FKR
      ENDIF

c      CALL IPRINT1(DRO,IMAX,JMAX,KMAX,IN,JN,KN,632)   

C **********************************************************************        
      K    = 0
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... Mersu, ei laske kahdesti

      ROT     = .5*(RO(I+II) + RO(I+II+IL))
      APT     = 1./APU(I+II) + 1./APU(I+II+IL)

      DPA     =   A3(I+II+IL) * (DP(I+II+IL)- DP(I+II))
      DUT     = -.5* ALFAU * A3X(I+II+IL)* APT * DPA
      DVT     = -.5* ALFAU * A3Y(I+II+IL)* APT * DPA
      DWT     = -.5* ALFAU * A3Z(I+II+IL)* APT * DPA
      DUT     = A3(I+II+IL)* A3X(I+II+IL)* ROT * DUT
      DVT     = A3(I+II+IL)* A3Y(I+II+IL)* ROT * DVT
      DWT     = A3(I+II+IL)* A3Z(I+II+IL)* ROT * DWT

      PRC(I+II+IL)%FKR   = PRC(I+II+IL)%F3R   + DUT + DVT + DWT
      DRO(I+II+IL)=  DRO(I+II+IL) + (- DUT - DVT - DWT)

      ENDDO; ENDDO

      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... Mersu, ei laske kahdesti

      ROT     = .5*(RO(I+II) + RO(I+II+IL))
      APT     = 1./APU(I+II) + 1./APU(I+II+IL)

      DPA     =   A3(I+II+IL) * (DP(I+II+IL)- DP(I+II))
      DUT     = -.5* ALFAU * A3X(I+II+IL)* APT * DPA
      DVT     = -.5* ALFAU * A3Y(I+II+IL)* APT * DPA
      DWT     = -.5* ALFAU * A3Z(I+II+IL)* APT * DPA
      DUT     = A3(I+II+IL)* A3X(I+II+IL)* ROT * DUT
      DVT     = A3(I+II+IL)* A3Y(I+II+IL)* ROT * DVT
      DWT     = A3(I+II+IL)* A3Z(I+II+IL)* ROT * DWT

      PRC(I+II+IL)%FKR   = PRC(I+II+IL)%F3R   + DUT + DVT + DWT
      DRO(I+II)=  DRO(I+II) - (- DUT - DVT - DWT)
      DRO(I+II+IL)=  DRO(I+II+IL) + (- DUT - DVT - DWT)
      PRC(I+II)%DRO = dro(i+ii)/VOL(I+II)   !(droa+1.e-12)
c       write(8000,*)i,j,droa,dro(i+ii),dro(i+ii)/(droa+1.e-12),dp(i+ii)
      ENDDO; ENDDO; ENDDO
        
C ***********************************************************************

      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTR + JJ + IN
      I    = 0

C ... I-direction, first row

      ROE     = .5*(RO(I+II) + RO(I+II+1))
      APE     = 1./APU(I+II) + 1./APU(I+II+1)

      DPA     =   A1(I+II+1) * (DP(I+II+1)- DP(I+II))
      DUE     = -.5* ALFAU * A1X(I+II+1)* APE * DPA
      DVE     = -.5* ALFAU * A1Y(I+II+1)* APE * DPA
      DWE     = -.5* ALFAU * A1Z(I+II+1)* APE * DPA
      DUE     = A1(I+II+1) * A1X(I+II+1)* ROE * DUE
      DVE     = A1(I+II+1) * A1Y(I+II+1)* ROE * DVE
      DWE     = A1(I+II+1) * A1Z(I+II+1)* ROE * DWE

      PRC(I+II+1)%FIR    = PRC(I+II+1)%F1R    + DUE + DVE + DWE
c      DROA = DROB(I+II)
      DRO(I+II+1)=  DRO(I+II+1) + (- DUE - DVE  - DWE)
      ENDDO; ENDDO
        
      DO K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      ROE     = .5*(RO(I+II) + RO(I+II+1))
      APE     = 1./APU(I+II) + 1./APU(I+II+1)

      DPA     =   A1(I+II+1) * (DP(I+II+1)- DP(I+II))
      DUE     = -.5* ALFAU * A1X(I+II+1)* APE * DPA
      DVE     = -.5* ALFAU * A1Y(I+II+1)* APE * DPA
      DWE     = -.5* ALFAU * A1Z(I+II+1)* APE * DPA
      DUE     = A1(I+II+1) * A1X(I+II+1)* ROE * DUE
      DVE     = A1(I+II+1) * A1Y(I+II+1)* ROE * DVE
      DWE     = A1(I+II+1) * A1Z(I+II+1)* ROE * DWE

      PRC(I+II+1)%FIR    = PRC(I+II+1)%F1R    + DUE + DVE + DWE
c      DROA = DROB(I+II)
c      IF(I == IMAX) PRC(I+II+1)%FIR    = 0. ! Test a recalculation
      DRO(I+II)=  DRO(I+II) - (- DUE - DVE - DWE)
      DRO(I+II+1)=  DRO(I+II+1) + (- DUE - DVE - DWE)
      PRC(I+II)%DRO = dro(i+ii)/VOL(I+II)   !(droa+1.e-12)
      ENDDO; ENDDO; ENDDO

C ***********************************************************************
        
      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      J    = 0
      II   = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

C ... J-direction, first row

      RON     = .5*(RO(I+II) + RO(I+II+ISTR))
      APN     = 1./APU(I+II) + 1./APU(I+II+ISTR)

      DPA     =   A2(I+II+ISTR) * (DP(I+II+ISTR)- DP(I+II))
      DUN     = -.5* ALFAU * A2X(I+II+ISTR)* APN * DPA
      DVN     = -.5* ALFAU * A2Y(I+II+ISTR)* APN * DPA
      DWN     = -.5* ALFAU * A2Z(I+II+ISTR)* APN * DPA
      DUN     =   A2(I+II+ISTR)*A2X(I+II+ISTR)*RON*DUN
      DVN     =   A2(I+II+ISTR)*A2Y(I+II+ISTR)*RON*DVN
      DWN     =   A2(I+II+ISTR)*A2Z(I+II+ISTR)*RON*DWN
      PRC(I+II+ISTR)%FJR = PRC(I+II+ISTR)%F2R + DUN + DVN + DWN
C      DROA = DROB(I+II)
      DRO(I+II+ISTR)=  DRO(I+II+ISTR) + (- DUN - DVN - DWN)
      ENDDO; ENDDO
        
      DO K = 1,KMAX
      JJ   = (KN+K-1)*IL
      DO J = 1,JMAX
      II   = (JN+J-1)*ISTR + JJ + IN
      DO I = 1,IMAX

      RON     = .5*(RO(I+II) + RO(I+II+ISTR))
      APN     = 1./APU(I+II) + 1./APU(I+II+ISTR)

      DPA     =   A2(I+II+ISTR) * (DP(I+II+ISTR)- DP(I+II))
      DUN     = -.5* ALFAU * A2X(I+II+ISTR)* APN * DPA
      DVN     = -.5* ALFAU * A2Y(I+II+ISTR)* APN * DPA
      DWN     = -.5* ALFAU * A2Z(I+II+ISTR)* APN * DPA
      DUN     =   A2(I+II+ISTR)*A2X(I+II+ISTR)*RON*DUN
      DVN     =   A2(I+II+ISTR)*A2Y(I+II+ISTR)*RON*DVN
      DWN     =   A2(I+II+ISTR)*A2Z(I+II+ISTR)*RON*DWN
      PRC(I+II+ISTR)%FJR = PRC(I+II+ISTR)%F2R + DUN + DVN + DWN
C      DROA = DROB(I+II)
      DRO(I+II)=  DRO(I+II) - (- DUN - DVN - DWN)
      DRO(I+II+ISTR)=  DRO(I+II+ISTR) + (- DUN - DVN - DWN)
      PRC(I+II)%DRO = dro(i+ii)/VOL(I+II)   !(droa+1.e-12)
      ENDDO; ENDDO; ENDDO

      RETURN
      END SUBROUTINE CORRM2
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRESMY(APU,DU,DV,DW,DP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,VOL,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + DPDX,DPDY,DPDZ,DH,DVP,IPT)

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IPT
      REAL ALFAU,APN,APS,APE,APW,APT,APD
      INTEGER N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL APU(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DP(*),DH(*),
     + A1(*),A2(*),A3(*),RO(*),DPDX(*),DPDY(*),DPDZ(*),VOL(*),DVP(*)
C
C ... Perform the correction of the velocity residual
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti (eikä laske)

      DU(I+II)=  DU(I+II) + ALFAU*DPDX(I+II)*VOL(I+II)
      DV(I+II)=  DV(I+II) + ALFAU*DPDY(I+II)*VOL(I+II)
      DW(I+II)=  DW(I+II) + ALFAU*DPDZ(I+II)*VOL(I+II)
      IF(IPT == 1) THEN
      DH(I+II)=  DH(I+II) + ALFAU*DVP(I+II)
      ENDIF
        
1000  CONTINUE

      RETURN
      END SUBROUTINE CORRESMY
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE GRADP(DPDX,DPDY,DPDZ,DP,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IPT)

      USE TYPE_ARRAYS

      IMPLICIT NONE


      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,IPT
      REAL VOLP,APN,APS,APE,APW,APT,APD,ROIP,ROIM,ROJP,ROJM,ROKP,ROKM
      INTEGER N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL DP(*),DPDX(*),DPDY(*),DPDZ(*),VOL(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),A1(*),A2(*),A3(*),
     + RO(*),DVP(*)
      TYPE(PRE_COR) PRC(*)

C
C ... Calculate the pressure gradient for internal cells
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      VOLP = 1./VOL(I+II)

      DPDX(I+II) =  - (
     + - A2(I+II+ISTR)*A2X(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2X(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1X(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3X(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))*VOLP
      DPDY(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Y(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Y(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Y(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))*VOLP
      DPDZ(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Z(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Z(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Z(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))*VOLP
        
      IF(IPT == 1) THEN ! Correct the energy residual

      ROIP = RO(I+II+1)+ RO(I+II)    ;  ROIM = RO(I+II-1)+ RO(I+II)
      ROJP = RO(I+II+ISTR)+ RO(I+II) ;  ROJM = RO(I+II-ISTR)+ RO(I+II)
      ROKP = RO(I+II+IL)+ RO(I+II)   ;  ROKM = RO(I+II-IL)+ RO(I+II)

      DVP(I+II) =  - ( 
     + - PRC(I+II+ISTR)%F2R *(DP(I+II+ISTR)+ DP(I+II))/ROJP
     + + PRC(I+II)%F2R * (DP(I+II)+ DP(I+II-ISTR))/ROJM
     + - PRC(I+II+1)%F1R*(DP(I+II+1)+ DP(I+II))/ROIP
     + + PRC(I+II)%F1R *(DP(I+II)+ DP(I+II-1))/ROIM
     + - PRC(I+II+IL)%F3R*(DP(I+II+IL)+ DP(I+II))/ROKP
     + + PRC(I+II)%F3R * (DP(I+II)+ DP(I+II-IL))*ROKM)
      ENDIF
        
1000  CONTINUE

      RETURN
      END SUBROUTINE GRADP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE GRADAP(DPDX,DPDY,DPDZ,DP,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,DVP,RO,PRC,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,VAR,
     + NPHASE,ALF1LG,ALF2LG,ALF3LG,NGL,IPH)

      USE TYPE_ARRAYS

      IMPLICIT NONE


      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NPHASE,NGL,IPH
      REAL VOLP,APN,APS,APE,APW,APT,APD,ROIP,ROIM,ROJP,ROJM,ROKP,ROKM
      INTEGER N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL DP(*),DPDX(*),DPDY(*),DPDZ(*),VOL(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),A1(*),A2(*),A3(*),
     + RO(*),DVP(*),
     + ALF1LG(NTOT,NPHASE),ALF2LG(NTOT,NPHASE),ALF3LG(NTOT,NPHASE)

      TYPE(PRE_COR) PRC(*)
      TYPE(MPHASE_VARIABLES) VAR(*)

C
C ... Calculate the pressure gradient for internal cells
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      VOLP = 1./VOL(I+II)

      DPDX(I+II) =  - (
     + - A2(I+II+ISTR)*A2X(I+II+ISTR) *.5*
     + (DP(I+II+ISTR)*ALF2LG(I+II+ISTR,IPH)+ DP(I+II)*ALF2LG(I+II,IPH))
     + + A2(I+II)     *A2X(I+II) *     .5*
     + (DP(I+II)*ALF2LG(I+II,IPH)+ DP(I+II-ISTR)*ALF2LG(I+II-ISTR,IPH))
     + - A1(I+II+1)   *A1X(I+II+1) *   .5*
     +  (DP(I+II+1)*ALF1LG(I+II+1,IPH)+ DP(I+II)*ALF1LG(I+II,IPH))
     + + A1(I+II)     *A1X(I+II) *     .5*
     +  (DP(I+II)*ALF1LG(I+II,IPH)+ DP(I+II-1)*ALF1LG(I+II-1,IPH))
     + - A3(I+II+IL)  *A3X(I+II+IL) *  .5*
     +  (DP(I+II+IL)*ALF3LG(I+II+IL,IPH)+ DP(I+II)*ALF3LG(I+II,IPH))
     + + A3(I+II)     *A3X(I+II) *     .5*
     +  (DP(I+II)*ALF3LG(I+II,IPH)+ DP(I+II-IL))*ALF3LG(I+II-IL,IPH))
     + *VOLP
      DPDY(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR) *.5*
     + (DP(I+II+ISTR)*ALF2LG(I+II+ISTR,IPH)+ DP(I+II)*ALF2LG(I+II,IPH))
     + + A2(I+II)     *A2Y(I+II) *     .5*
     + (DP(I+II)*ALF2LG(I+II,IPH)+ DP(I+II-ISTR)*ALF2LG(I+II-ISTR,IPH))
     + - A1(I+II+1)   *A1Y(I+II+1) *   .5*
     +  (DP(I+II+1)*ALF1LG(I+II+1,IPH)+ DP(I+II)*ALF1LG(I+II,IPH))
     + + A1(I+II)     *A1Y(I+II) *     .5*
     +  (DP(I+II)*ALF1LG(I+II,IPH)+ DP(I+II-1)*ALF1LG(I+II-1,IPH))
     + - A3(I+II+IL)  *A3Y(I+II+IL) *  .5*
     +  (DP(I+II+IL)*ALF3LG(I+II+IL,IPH)+ DP(I+II)*ALF3LG(I+II,IPH))
     + + A3(I+II)     *A3Y(I+II) *     .5*
     +  (DP(I+II)*ALF3LG(I+II,IPH)+ DP(I+II-IL))*ALF3LG(I+II-IL,IPH))
     + *VOLP
      DPDZ(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR) *.5*
     + (DP(I+II+ISTR)*ALF2LG(I+II+ISTR,IPH)+ DP(I+II)*ALF2LG(I+II,IPH))
     + + A2(I+II)     *A2Z(I+II) *     .5*
     + (DP(I+II)*ALF2LG(I+II,IPH)+ DP(I+II-ISTR)*ALF2LG(I+II-ISTR,IPH))
     + - A1(I+II+1)   *A1Z(I+II+1) *   .5*
     +  (DP(I+II+1)*ALF1LG(I+II+1,IPH)+ DP(I+II)*ALF1LG(I+II,IPH))
     + + A1(I+II)     *A1Z(I+II) *     .5*
     +  (DP(I+II)*ALF1LG(I+II,IPH)+ DP(I+II-1)*ALF1LG(I+II-1,IPH))
     + - A3(I+II+IL)  *A3Z(I+II+IL) *  .5*
     +  (DP(I+II+IL)*ALF3LG(I+II+IL,IPH)+ DP(I+II)*ALF3LG(I+II,IPH))
     + + A3(I+II)     *A3Z(I+II) *     .5*
     +  (DP(I+II)*ALF3LG(I+II,IPH)+ DP(I+II-IL))*ALF3LG(I+II-IL,IPH))
     + *VOLP
        
1000  CONTINUE

      RETURN
      END SUBROUTINE GRADAP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE GRADP_tul(DPDX,DPDY,DPDZ,DP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL :: APN,APS,APE,APW,APT,APD
      REAL :: DP(*),DPDX(*),DPDY(*),DPDZ(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),A1(*),A2(*),A3(*)
C
C ... Perform the correction of the velocity residual
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      DPDX(I+II) =  - (
     + - A2(I+II+ISTR)*A2X(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2X(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1X(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3X(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))
      DPDY(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Y(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Y(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Y(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))
      DPDZ(I+II) =  - ( 
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR) *.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Z(I+II) *     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1) *   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Z(I+II) *     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL) *  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Z(I+II) *     .5*(DP(I+II)+ DP(I+II-IL)))
c        if(i == 2 .and. j == 2 .and. k == 1) then
c            write(6666,*) i+ii,dpdx(i+ii),dpdy(i+ii),dpdz(i+ii)
c            write(6667,*)dp(i+ii),dp(i+ii-istr),dp(i+ii+istr),
c     +      dp(i+ii-1),dp(i+ii+1),dp(i+ii-il),dp(i+ii-il)
c        endif
1000  CONTINUE

      RETURN
      END SUBROUTINE GRADP_tul
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C            
      SUBROUTINE CORRV(APU,DU,DV,DW,DP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,RKSI,VAR,MULPHC)

      USE TYPE_ARRAYS
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ,IPHASE
      REAL :: ALFAU,APN,APS,APE,APW,APT,APD,ARKSI
      REAL :: APU(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DP(*),
     + A1(*),A2(*),A3(*),RO(*),XFC(*),YFC(*),ZFC(*),RKSI(*)
      CHARACTER (LEN=10) :: MULPHC
       TYPE(MPHASE_VARIABLES) VAR(*)

C
C ... Perform the velocity correction by calculating the pressure gradient
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

c      APN     = .5*(APU(I+II) + APU(I+II+ISTR))
c      APS     = .5*(APU(I+II) + APU(I+II-ISTR))
c      APE     = .5*(APU(I+II) + APU(I+II+1))
c      APW     = .5*(APU(I+II) + APU(I+II-1))
c      APT     = .5*(APU(I+II) + APU(I+II+IL))
c      APD     = .5*(APU(I+II) + APU(I+II-IL))

      APN     = APU(I+II)
      APS     = APU(I+II)
      APE     = APU(I+II)
      APW     = APU(I+II)
      APT     = APU(I+II)
      APD     = APU(I+II) ! Mersu

      ARKSI   = 1./(1.+RKSI(I+II)) ! This should make the correction zero

      DU(I+II)=  DU(I+II) + ALFAU*ARKSI*(
     + - A2(I+II+ISTR)*A2X(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2X(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1X(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3X(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      DV(I+II)=  DV(I+II) + ALFAU*ARKSI*( 
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Y(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Y(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL)/APT*    .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Y(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      DW(I+II)=  DW(I+II) + ALFAU*ARKSI*( 
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Z(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Z(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Z(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
        
1000  CONTINUE

c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   

      RETURN
      END SUBROUTINE CORRV
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRVMF(APU,DU,DV,DW,DP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,RKSI,VAR,MULPHC,APUULG)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ,IPHASE
      REAL :: ALFAU,APN,APS,APE,APW,APT,APD,ARKSI
      REAL :: APU(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DP(*),
     + A1(*),A2(*),A3(*),RO(*),XFC(*),YFC(*),ZFC(*),RKSI(*),
     + APUULG(NTOT,NPHASES)

      CHARACTER (LEN=10) :: MULPHC

      TYPE(MPHASE_VARIABLES) VAR(*)

C
C ... Perform the velocity correction by calculating the pressure gradient
C

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

c      APN     = .5*(APU(I+II) + APU(I+II+ISTR))
c      APS     = .5*(APU(I+II) + APU(I+II-ISTR))
c      APE     = .5*(APU(I+II) + APU(I+II+1))
c      APW     = .5*(APU(I+II) + APU(I+II-1))
c      APT     = .5*(APU(I+II) + APU(I+II+IL))
c      APD     = .5*(APU(I+II) + APU(I+II-IL))

      APN     = APU(I+II)
      APS     = APU(I+II)
      APE     = APU(I+II)
      APW     = APU(I+II)
      APT     = APU(I+II)
      APD     = APU(I+II) ! Mersu

      ARKSI   = 1./(1.+RKSI(I+II)) ! This should make the correction zero

      DU(I+II)=  DU(I+II) + ALFAU*ARKSI*(
     + - A2(I+II+ISTR)*A2X(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2X(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1X(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3X(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      DV(I+II)=  DV(I+II) + ALFAU*ARKSI*( 
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Y(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Y(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL)/APT*    .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Y(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      DW(I+II)=  DW(I+II) + ALFAU*ARKSI*( 
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Z(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Z(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Z(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
        
      IF(MULPHC == 'MULTI') THEN
      DO IPHASE = 1,NPHASES
      APN     = APUULG(I+II,IPHASE)
      APS     = APUULG(I+II,IPHASE)
      APE     = APUULG(I+II,IPHASE)
      APW     = APUULG(I+II,IPHASE)
      APT     = APUULG(I+II,IPHASE)
      APD     = APUULG(I+II,IPHASE)

      VAR(I+II)%DM(IPHASE) =  VAR(I+II)%DM(IPHASE) + ALFAU*ARKSI*(
     + - A2(I+II+ISTR)*A2X(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2X(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1X(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1X(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3X(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3X(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      VAR(I+II)%DN(IPHASE) =  VAR(I+II)%DN(IPHASE) + ALFAU*ARKSI*(
     + - A2(I+II+ISTR)*A2Y(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Y(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Y(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Y(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Y(I+II+IL)/APT*    .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Y(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      VAR(I+II)%DW(IPHASE) =  VAR(I+II)%DW(IPHASE) + ALFAU*ARKSI*(
     + - A2(I+II+ISTR)*A2Z(I+II+ISTR)/APN*.5*(DP(I+II+ISTR)+ DP(I+II))
     + + A2(I+II)     *A2Z(I+II)/APS*     .5*(DP(I+II)+ DP(I+II-ISTR))
     + - A1(I+II+1)   *A1Z(I+II+1)/APE*   .5*(DP(I+II+1)+ DP(I+II))
     + + A1(I+II)     *A1Z(I+II)/APW*     .5*(DP(I+II)+ DP(I+II-1))
     + - A3(I+II+IL)  *A3Z(I+II+IL)/APT*  .5*(DP(I+II+IL)+ DP(I+II))
     + + A3(I+II)     *A3Z(I+II)/APD*     .5*(DP(I+II)+ DP(I+II-IL)))
      ENDDO ! IPHASE
      ENDIF ! MULPHC        
1000  CONTINUE

c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   

      RETURN
      END SUBROUTINE CORRVMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRVY(APU,DU,DV,DW,DP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,RO,
     + XFC,YFC,ZFC,DPDX,DPDY,DPDZ,VOL,RKSI,VAR,MULPHC)

      USE TYPE_ARRAYS
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ,IPHASE
      REAL :: ALFAU,APN,APS,APE,APW,APT,APD,ARKSI,DUS,DVS,DWS,
     + APX,APY,APZ
      REAL :: APU(*),DU(*),DV(*),DW(*),A1X(*),A1Y(*),A1Z(*),
     + A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),DP(*),
     + A1(*),A2(*),A3(*),RO(*),XFC(*),YFC(*),ZFC(*),
     + DPDX(*),DPDY(*),DPDZ(*),VOL(*),RKSI(*)
       TYPE(MPHASE_VARIABLES) VAR(*)
       CHARACTER (LEN=10) :: MULPHC
C
C ... Perform the velocity correction using known pressure gradient
C
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
        
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

      APN     = APU(I+II)
      APS     = APU(I+II)
      APE     = APU(I+II)
      APW     = APU(I+II)
      APT     = APU(I+II)
      APD     = APU(I+II) ! Mersu
      ARKSI   = 1./(1.+RKSI(I+II)) ! This should make the correction zero

      DU(I+II)=  DU(I+II) - ALFAU*ARKSI*DPDX(I+II)/APN*VOL(I+II)
      DV(I+II)=  DV(I+II) - ALFAU*ARKSI*DPDY(I+II)/APN*VOL(I+II)
      DW(I+II)=  DW(I+II) - ALFAU*ARKSI*DPDZ(I+II)/APN*VOL(I+II)
        
1000  CONTINUE

      RETURN
      END SUBROUTINE CORRVY
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRVYMF(IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,
     + DPDX,DPDY,DPDZ,VOL,RKSI,VAR,AP1,AP2,AP3,
     + NPHASE,ALF1LG,ALF2LG,ALF3LG,NGL,IPHASE)

      USE TYPE_ARRAYS
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ,IPHASE,NPHASE
      REAL :: ALFAU,APN,APS,APE,APW,APT,APD,ARKSI,DUS,DVS,DWS,
     + APX,APY,APZ
      REAL :: DPDX(*),DPDY(*),DPDZ(*),VOL(*),RKSI(*)
      REAL:: AP1(NTOT,NPHASE),AP2(NTOT,NPHASE),AP3(NTOT,NPHASE),
     +       ALF1LG(NTOT,NPHASE),ALF2LG(NTOT,NPHASE),ALF3LG(NTOT,NPHASE)
      TYPE(MPHASE_VARIABLES) VAR(*)
C
C ... Velocity correction using known pressure gradient for a two-fluid
C
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
      
      DO 2000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 1,IMAX

c      APN     = APU(I+II) !ULG(I+II,IPHASE)

      APX     = AP1(I+II,IPHASE)
      APY     = AP2(I+II,IPHASE)
      APZ     = AP3(I+II,IPHASE)

      ARKSI   = 1./(1.+RKSI(I+II)) ! This should make the correction zero
c      if(ngl  == 2) then
c       write(7002,*) i,j,k,VAR(I+II)%DM(IPHASE),
c     + ALFLG(I+II,IPHASE)*ALFAU*ARKSI*DPDX(I+II)/APX*VOL(I+II)
c      endif
      VAR(I+II)%DM(IPHASE) =  VAR(I+II)%DM(IPHASE) 
     +        - ALF1LG(I+II,IPHASE)*ALFAU*ARKSI*DPDX(I+II)/APX*VOL(I+II)
      VAR(I+II)%DN(IPHASE) =  VAR(I+II)%DN(IPHASE) 
     +        - ALF2LG(I+II,IPHASE)*ALFAU*ARKSI*DPDY(I+II)/APY*VOL(I+II)
      VAR(I+II)%DW(IPHASE) =  VAR(I+II)%DW(IPHASE) 
     +        - ALF3LG(I+II,IPHASE)*ALFAU*ARKSI*DPDZ(I+II)/APZ*VOL(I+II)
        
2000  CONTINUE

      RETURN
      END SUBROUTINE CORRVYMF
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRVYMFA(IMAX,JMAX,KMAX,NTOT,IN,JN,KN,ALFAU,
     + DPDX,DPDY,DPDZ,VOL,RKSI,VAR,AP1,AP2,AP3,
     + NPHASE,NGL,IPHASE)

      USE TYPE_ARRAYS
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,NGL
      INTEGER :: N,ISTR,JSTR,IL,I,J,K,II,JJ,IPHASE,NPHASE
      REAL :: ALFAU,ARKSI,APX,APY,APZ
      REAL ::  DPDX(*),DPDY(*),DPDZ(*),VOL(*),RKSI(*)
      REAL:: AP1(NTOT,NPHASE),AP2(NTOT,NPHASE),AP3(NTOT,NPHASE)
      TYPE(MPHASE_VARIABLES) VAR(*)
C
C ... Velocity correction using void-pressure gradient for a two-fluid
C
      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
      
      DO 2000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 1,IMAX

c      APN     = APU(I+II) !ULG(I+II,IPHASE)

      APX     = AP1(I+II,IPHASE)
      APY     = AP2(I+II,IPHASE)
      APZ     = AP3(I+II,IPHASE)

      ARKSI   = 1./(1.+RKSI(I+II)) ! This should make the correction zero

      VAR(I+II)%DM(IPHASE) =  VAR(I+II)%DM(IPHASE) 
     +        - ALFAU*ARKSI*DPDX(I+II)/APX*VOL(I+II)
      VAR(I+II)%DN(IPHASE) =  VAR(I+II)%DN(IPHASE) 
     +        - ALFAU*ARKSI*DPDY(I+II)/APY*VOL(I+II)
      VAR(I+II)%DW(IPHASE) =  VAR(I+II)%DW(IPHASE) 
     +        - ALFAU*ARKSI*DPDZ(I+II)/APZ*VOL(I+II)
c        if(ngl == 2 .and.i == 18) write(9002,*)
c     + j,iphase,VAR(I+II)%DM(IPHASE),ALF1LG(I+II,IPHASE),DPDX(I+II),
c     + APX,alfau
        
2000  CONTINUE

      RETURN
      END SUBROUTINE CORRVYMFA
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOURAP(AP,RK,REPS,E,SRK,SEPS,PTUR,EPS2,VOL,DTL,NTOT,
     + M,RO,VIS,CSIMPS,RKMAX,DTURB,CE1,CE2,KOVER,ICASE,TRM)	

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ICASE,NTOT,M,N,KOVER
      REAL :: RKN,ELIM1,ELIM2,SOURCE,DTURB,CSIMPS,CE1,CE2,RKMAX,RKLIM,
     + ELIM,EPS,A1KLEB,BSTAR,CKAPPA,BETA1,SIGRK1,SIGOM1,BETA2,SIGRK2,
     + SIGOM2,CKMAX,EMEANF,DRKMAX,PROL,ATPROD,COMAX,GAMKOM,DOMMAX,
     + PROL2,RNUTRB,EPSLIM,REPSA,SIGPHI,ZETA2
      REAL :: AP(*),RK(*),REPS(*),SRK(*),SEPS(*),VOL(*),DTL(*),
     + PTUR(*),RO(*),VIS(*),E(*),EPS2(*)

      TYPE(INTERMITTENCY) TRM(*)

      EPS     = 1.E-20
      RKLIM   = 0.
      EPSLIM  = 1.E-10

      SELECT CASE(ICASE)

      CASE(1)  ! K-epsilon model, k-equation

      DO N = 1,NTOT
         RKN     = 1./(RK(N)+EPS)
         ELIM1   = 1./(ABS(0.1*E(N)-RK(N))+EPS)
         ELIM2   = 1./(ABS(RKLIM-RK(N))+EPS)
         ELIM    = MAX(ELIM1,ELIM2)
         SOURCE  = CSIMPS*VOL(N)*ABS(PTUR(N))*MAX(RKN,ELIM)
         REPSA   = ABS(REPS(N))
c         AP(N)   = AP(N) + 2.*VOL(N)*REPS(N)*RKN + 2.*SOURCE*RO(N)
         AP(N)   = AP(N) + (2.*VOL(N)*REPSA*RKN + SOURCE)*RO(N) ! Added 2 29.10.2010
      ENDDO

      CASE(2)  ! K-epsilon model, epsilon-equation

      DO N = 1,NTOT
         RKN     = 1./(RK(N)+EPS)
         ELIM1   = 1./(ABS(0.1*E(N)-RK(N))+EPS)
         ELIM2   = 1./(ABS(RKLIM-RK(N))+EPS)
         ELIM    = MAX(ELIM1,ELIM2)
         SOURCE  = CSIMPS*VOL(N)*ABS(PTUR(N))*MAX(RKN,ELIM)
         REPSA   = ABS(REPS(N))
         AP(N)   = AP(N) + (2.*CE2*VOL(N)*REPSA*RKN+CE1*SOURCE)*RO(N)
      ENDDO

      CASE(3)  ! K-omega-SST model, k-equation (testaamatta)

C     Model coefficients:
      CALL COEFKO(6,A1KLEB,BSTAR,CKAPPA,
     &     BETA1,SIGRK1,SIGOM1,BETA2,SIGRK2,SIGOM2,
     &     SIGPHI,ZETA2,KOVER)
      CKMAX  = 1.0/(2.0*CSIMPS)

      DO N = 1,NTOT
         PROL   = 1./RO(N)
         EMEANF = E(N) - RK(N)
         DRKMAX = MAX(CKMAX*MIN(RK(N),ABS(0.1*EMEANF-0.9*RK(N))),
     &           EPSLIM)*PROL
         ATPROD = ABS(PTUR(N))
         AP(N)  = AP(N) + VOL(N)*(ATPROD/DRKMAX+BSTAR*ABS(REPS(N)))
      ENDDO

      CASE(4)  ! K-omega-SST model, omega-equation (testaamatta)

C     Model coefficients:
      CALL COEFKO(6,A1KLEB,BSTAR,CKAPPA,
     &     BETA1,SIGRK1,SIGOM1,BETA2,SIGRK2,SIGOM2,
     &     SIGPHI,ZETA2,KOVER)
      GAMKOM = BETA1/BSTAR - SIGOM1*CKAPPA**2/(SQRT(BSTAR))
      COMAX  = 1.0/(2.0*CSIMPS)

      DO N = 1,NTOT
         PROL   = 1./RO(N)
         PROL2  = 2.0*PROL
         DOMMAX = MAX(COMAX*REPS(N)*PROL,EPSLIM)
         ATPROD = ABS(PTUR(N))
         RNUTRB = MAX((EPS2(N)-1.0),EPSLIM)*VIS(N)*PROL
         AP(N)  = AP(N) + VOL(N)*(ABS(GAMKOM*ATPROD/(RNUTRB*DOMMAX)) 
     &          + 2.*BETA1*ABS(REPS(N)) )   ! F4 = 1 forced
c    &          + 2.*BETA1*ABS(REPS(N)) + ABS(SEPS(N)/DOMMAX))  ! F4 = 1 forced, kokeiltu
      ENDDO

      CASE(5)  ! Spalart-Allmaras model, unused

      CASE(6)  ! Spalart-Allmaras model

      DO N = 1,NTOT
         RKN     = 1./(REPS(N)+EPS)
         ELIM2   = 1./(ABS(EPSLIM-REPS(N))+EPS)*5.
c         AP(N)   = AP(N) + VOL(N)*ABS(SEPS(N))*ELIM2
         AP(N)   = AP(N) + 
     &   VOL(N)*(ABS(PTUR(N))*ELIM2 + 2.*ABS(SRK(N)/(REPS(N)+EPS)))
      ENDDO

C     Intermittency variables

      CASE(7)  !G

      DO N = 1,NTOT
         AP(N)   = AP(N) + VOL(N)*(
     &             ABS( TRM(N)%SG/(2.*0.1+EPSLIM) )
     &           + ABS( TRM(N)%LG ) )
      ENDDO

      CASE(8)  !RET

      DO N = 1,NTOT
         AP(N)   = AP(N) + VOL(N)*(
     &             ABS( TRM(N)%SRET/(2000.*0.1+EPSLIM) )
     &           + ABS( TRM(N)%LRET ) )
      ENDDO


      END SELECT
      END SUBROUTINE SOURAP
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOURSC(AP,SFI,FI,VOL,DTL,NTOT,NS,SCLIM)

C ... Source term linearization for the scalar equations

      IMPLICIT NONE

      INTEGER :: NS,NTOT,N
      REAL :: SCLIM
      REAL :: AP(*),SFI(*),FI(*),VOL(*),DTL(*)

      DO N = 1,NTOT
         AP(N)   = AP(N) + VOL(N)*ABS(SFI(N))/(FI(N) + SCLIM)
      ENDDO

      RETURN
      END SUBROUTINE SOURSC
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOURTR(AP,TRM,VOL,DTL,NTOT,TRLIM,ICASE)

C ... Source term linearization for the transition equations

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ICASE,NTOT,N
      REAL :: TRLIM
      REAL :: AP(*),VOL(*),DTL(*)

      TYPE(INTERMITTENCY) :: TRM(*)

      SELECT CASE(ICASE)

      CASE(1)  ! G-equation

      DO N = 1,NTOT
         AP(N)   = AP(N) + VOL(N)*ABS(TRM(N)%SG)/(TRM(N)%G + TRLIM)
      ENDDO

      CASE(2)  ! RET-equation

      DO N = 1,NTOT
         AP(N)   = AP(N) + VOL(N)*ABS(TRM(N)%SRET)/(TRM(N)%RET + TRLIM)
      ENDDO

      END SELECT

      END SUBROUTINE SOURTR
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFUM(DIFF,VIST,EPS2,VIS,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN)
C ... Diffusion multipliers for the momentum equation

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL
      REAL :: DIFF(*),VIST(*),VIS(*),EPS2(*)
      REAL :: PRT, APU

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = -1,KMAX+2 ! Extended just to be sure
      JJ      = (KN+K-1)*IL
      DO J = -1,JMAX+2
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = -1,IMAX+2
      DIFF(I+II) = VIS(I+II)*EPS2(I+II)  !  + VIST(I+II)
c      APU  = (VIS(I+II)*EPS2(I+II) - VIST(I+II))/VIS(I+II)
C      IF(ABS(APU) > 1.E-3) 
c     +write(9600+k,*)i,j,vist(i+ii),(eps2(i+ii)-1.)*vis(i+ii),eps2(i+II)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIFFUM
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFUE(DIFF,VIST,CH,CP,PRT,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,EPS2,VIS)
C ... Diffusion multipliers for the energy equation

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL
      REAL :: DIFF(*),VIST(*),CH(*),CP(*),EPS2(*),VIS(*)
      REAL :: PRT

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = CH(I+II)/CP(I+II) + (EPS2(I+II)-1.)*VIS(I+II)/PRT  
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFFUE
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFUK(DIFF,VIST,VIS,SC,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,EPS2)

C ... Diffusion multipliers for the turbulence equations

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL
      REAL :: DIFF(*),VIST(*),VIS(*),EPS2(*)
      REAL :: SC

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = VIS(I+II) + (EPS2(I+II)-1.)*VIS(I+II) *SC ! Corr. 2.11.10
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFFUK
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFUS(DIFF,REPS,VIS,SC,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN)

C ... Diffusion multipliers for the Spalart-Allmaras model

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL
      REAL :: DIFF(*),REPS(*),VIS(*)
      REAL :: SC

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = (VIS(I+II) + REPS(I+II))/SC
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFFUS
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFMFM(DIFF,VIS,EPS2,PRO,VAR,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

C ... Diffusion multipliers for the two-fluid momentum equations

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL,
     + IPHASE
      REAL :: DIFF(*),VIS(*),EPS2(*)
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = -1,KMAX+2
      JJ      = (KN+K-1)*IL
      DO J = -1,JMAX+2
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = -1,IMAX+2
      DIFF(I+II) = PRO(I+II)%VIS(IPHASE)*EPS2(I+II) *
     +            VAR(I+II)%ALFA(IPHASE) 
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFMFM
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFMFE(DIFF,VIS,EPS2,PRO,VAR,PRT,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

C ... Diffusion multipliers for the multiphase energy equations

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL,
     + IPHASE
      REAL :: DIFF(*),VIS(*),EPS2(*)
      REAL :: PRT
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = (PRO(I+II)%CH(IPHASE)/PRO(I+II)%CP(IPHASE) + 
     +           (EPS2(I+II)-1.)*PRO(I+II)%VIS(IPHASE)/PRT)   * 
     +            VAR(I+II)%ALFA(IPHASE)
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFMFE
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFMFK(DIFF,VIS,EPS2,PRO,VAR,SC,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,IPHASE)

C ... Diffusion multipliers for the two-fluid turbulence equations (unused)

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL,
     + IPHASE
      REAL :: DIFF(*),VIS(*),EPS2(*)
      REAL :: SC
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = PRO(I+II)%VIS(IPHASE) + 
     +             (EPS2(I+II)-1.)*PRO(I+II)%VIS(IPHASE)/SC
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFMFK
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFSC(DIFF,VIS,EPS2,PSIGSC,PSIGS2,IMAX,JMAX,KMAX,NTOT,
     + IN,JN,KN,NS)

C ... Diffusion multipliers for the scalar equations

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL,
     + NS
      REAL :: DIFF(*),VIS(*),EPS2(*),PSIGSC(*),PSIGS2(*)
      REAL :: PSIGE1,PSIGE2

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

      PSIGE1  = 1./PSIGSC(NS) ! Turbulent
      PSIGE2  = 1./PSIGS2(NS) ! Laminar Schmidt number

      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1

      DIFF(I+II) = VIS(I+II)*PSIGE2 + (EPS2(I+II)-1.)*VIS(I+II)*PSIGE1
      ENDDO ; ENDDO ; ENDDO

      RETURN
      END SUBROUTINE DIFFSC
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIFFTR(DIFF,VIS,EPS2,TRM,IMAX,JMAX,KMAX,NTOT,IN,JN,KN,
     +     ICASE)

C ... Diffusion multipliers for the transition equations

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,NTOT,IN,JN,KN,I,J,K,II,JJ,ISTR,JSTR,IL,
     + ICASE
      REAL :: DIFF(*),VIS(*),EPS2(*)
      REAL :: CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1

      TYPE(INTERMITTENCY) :: TRM(*)

      CALL  TRANVAR(CE1,CE2,CA1,CA2,SIGF,COT,SIGOT,S1)

      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR

C ... For G
      IF(ICASE == 1) THEN 
      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) =   VIS(I+II) + (EPS2(I+II)-1.)*VIS(I+II) / SIGF
      ENDDO ; ENDDO ; ENDDO

C ... For RET
      ELSE IF(ICASE == 2) THEN
      DO K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      DO J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO I = 0,IMAX+1
      DIFF(I+II) = ( VIS(I+II) + (EPS2(I+II)-1.)*VIS(I+II) ) * SIGOT
      ENDDO ; ENDDO ; ENDDO
	
      ELSE
         WRITE(*,*) 'No such case number in intermittency variables'
         STOP
      ENDIF

      RETURN
      END SUBROUTINE DIFFTR
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      
      SUBROUTINE EVAP_CORR(PRO,VAR,BLKS,U,V,W,P,PDIFF,VIS,EPS2,FRSDEN,
     + FRSVEL,IMAX,JMAX,KMAX,IN,JN,KN,DT,DTL,NGL,M,PRT,VOL,DH,AP,ALIMIT,
     + IPHASE)

      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : PII,EPS20

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,I,K,N,ISTRID,JSTRID,IL,IA,
     + IPHASE,NGL,II,JJ,J,ISTR,JSTR,KSTR,IJ,J1,IEVAP,mmm,nnn,lll,
     + M
      REAL :: U(*),V(*),W(*),P(*),DTL(*),PDIFF(*),VIS(*),EPS2(*),DH(*),
     + AP(*),VOL(*)
      REAL :: FRSDEN,VEL2,DT,EVAPH,EVAPR,APU,TEMDIF,QILMIN,QIGMAX,
     + HFG,PPSAT,RKF,RKB,FRSVEL,CHLREF,EVAPOLD,REFVEL,VELKIN,YPLS,
     + RLAM,B1E,ALMAX,ALMIN,PRT,EVAPMAX,DTM,ALIMIT,APSAT
      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(BLOCKS)           :: BLKS(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      IL      = ISTRID*JSTRID
        
      ISTR    = 1
      JSTR    = ISTRID
      KSTR    = IL
      IEVAP   = BLKS(NGL)%EVAP_TYPE
      REFVEL  = BLKS(NGL)%REFCAV
      CHLREF  = BLKS(NGL)%CHLCAV

      DO K = 1,KMAX
      IA = (K+KN-1)*KSTR + JN*JSTR + (IN-1)*ISTR + 1
      DO IJ = 1,IMAX*JMAX
         J1 = (IJ-1)/IMAX           ! J1 = J-1
         I  = IJ - J1*IMAX
         N  = IA + J1*JSTR + I*ISTR ! (I-1)*ISTR
         DTM    = DTL(N)*1.E-4 + EPS20
         APSAT  = MAX(0.,ALIMIT-VAR(N)%ALFA(IPHASE))/(ALIMIT+EPS20)*
     +            VOL(N)/DTM * PRO(N)%RO(IPHASE)
         VAR(N)%ADDE(IPHASE) = APSAT*(PRO(N)%DTEMP(IPHASE)-PRO(N)%TSAT)
     +          * PRO(N)%CP(IPHASE)
         DH(N)  = VAR(N)%DH(IPHASE)  + VAR(N)%ADDE(IPHASE)
c     +          + APSAT*(PRO(N)%DTEMP(IPHASE) - PRO(N)%TSAT)
c     +          * PRO(N)%CP(IPHASE)
         AP(N)  = AP(N) + APSAT
c      write(7001,*) iphase,i,j1+1,VAR(N)%DH(IPHASE),
c     + VAR(N)%ADDE(IPHASE)
      ENDDO ; ENDDO

      RETURN
      END SUBROUTINE EVAP_CORR
C
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------     
C ----------------------------------------------------------------------
C      

