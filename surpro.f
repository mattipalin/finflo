C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SURPRO(UWALL,VWALL,WWALL,XCP,YCP,ZCP,
     +     XCO,YCO,ZCO,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     +     A3X,A3Y,A3Z,A1,A2,A3,OMEGA,OMEX,OMEY,OMEZ,
     +     CENX,CENY,CENZ,IMAX,JMAX,KMAX,
     +     IDI1,IDI2,IDI3,INRE,JNRE,KNRE,NBL,M,NPATCH,
     +     ICON,RCON,IHF,NSPTOT,NSOLPA,NORMAL,
     +     IDIMS,JDIMS,KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,
     +     ISTRS,JSTRS,KSTRS,NSPB,NSPP,XC,YC,ZC,
     +     NSPG,APATCH,IREPEA,NREPEA,XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,DT,
     +     TIMEL,IGRID,POROS,RMLOSS,IDIMSG,JDIMSG,BOUNDN,IACTNO,
     +     SURFA,SURFT,IALL,IPRO,DTOLD,GVEX,GVEY,GVEZ)

      USE FLIGHT, ONLY    : IFA
      USE NS3CO, ONLY     : IC9, LN
      USE CONSTANTS, ONLY : EPS10,EPS6

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +     IDI1,IDI2,IDI3,NBL,M,NPATCH,ISTRID,JSTRID,KSTRID,
     +     IL,IP,IBC,IPL,IDIR,IFACE,ISTR,JSTR,KSTR,I,J,K,ID,NFL,
     +     KX1,KX2,KY1,KY2,LSTRID,KA,II,IJ,L,L2,L3,L4,JST2,
     +     IMOV,IAP,NR,NP,NSPTOT,NSOLPA,NROT,NMOV1,NMOV2,NMOV3,ISP,
     +     ITYPE,NORM,II1,KBIND,NTOT,IPG,NMOV4,IGRID,
     +     IDIMP,LENGTH,IACTU,IACTY,IACTNO,IALL,IPRO
      INTEGER :: ICON(IC9,*),IHF(*),IREPEA(*),NREPEA(*)
      INTEGER :: NORMAL(*),IDIMS(*),JDIMS(*),KX1S(*),KX2S(*),KY1S(*),
     +     KY2S(*),KZ1S(*),KZ2S(*),ISTRS(*),JSTRS(*),KSTRS(*),NSPB(*),
     +     NSPP(*),NSPG(*),IDIMSG(*),JDIMSG(*)

      INTEGER :: III
      
      REAL :: YCO(*), ZCO(*), XCO(*), XC(*), YC(*), ZC(*)
      REAL :: XLE2(*), YLE2(*), ZLE2(*)
      REAL :: XLE3(*), YLE3(*), ZLE3(*)
      REAL :: XCP(*),YCP(*),ZCP(*)
      REAL :: XP1, YP1, ZP1, XP0, YP0, ZP0

      REAL :: RCON(IC9,*),UWALL(*),VWALL(*),
     +     WWALL(*),
     +     A1X(*),A1Y(*),A1Z(*),A2X(*),A2Y(*),
     +     A2Z(*),A3X(*),A3Y(*),A3Z(*),A1(*),A2(*),A3(*),APATCH(*),
     +     POROS(*),RMLOSS(*),SURFA(*),SURFT(*)
      REAL :: OMEGA,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,UBIU,UBIV,UBIW,VELAVE,
     +     XAVE,YAVE,ZAVE,XVEL,YVEL,ZVEL,VEL,XN1,XN2,XN3,XA1,XA2,XA3,
     +     XNORM,XNOR,R1,R2,R3,B1,B2,B3,ARGSQR,DT,XCPLE2,YCPLE2,ZCPLE2,
     +     XCPLE3,YCPLE3,ZCPLE3,DTOLD,GVEX,GVEY,GVEZ

      REAL, ALLOCATABLE :: AAX(:,:), AAA(:,:)

      LOGICAL :: HAVRES,FIRSTL,TIMEL,THERE  

      CHARACTER(LEN=80)  :: BOUNDN(*)
      CHARACTER(LEN=180) :: LINE1
 

C ... Calculate surface velocities for the solid walls. The patches
C ... are extended by one row of cells for postprocessing
C ... (Moving solids are not tested yet 20.2.2003...)
      
      KN = KNRE
      JN = JNRE
      IN = INRE

      KSTRID = KMAX + 2*KN
      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      NTOT   = ISTRID*JSTRID*KSTRID
      IL     = ISTRID*JSTRID
      ISP    = NSPTOT
      NSOLPA = 0
      NROT   = 0
      NMOV1  = 0
      NMOV2  = 0
      NMOV3  = 0
      NMOV4  = 0
      HAVRES = .TRUE.
      FIRSTL = .TRUE.
      LENGTH = 0

C ... Antakee mersu

      ALLOCATE(AAX(NTOT,9))
      ALLOCATE(AAA(NTOT,3))
      DO L = 1,NTOT
         AAX(L,1) = A1X(L)
         AAX(L,2) = A1Y(L)
         AAX(L,3) = A1Z(L)
         AAX(L,4) = A2X(L)
         AAX(L,5) = A2Y(L)
         AAX(L,6) = A2Z(L)
         AAX(L,7) = A3X(L)
         AAX(L,8) = A3Y(L)
         AAX(L,9) = A3Z(L)
         AAA(L,1) = A1(L)
         AAA(L,2) = A2(L)
         AAA(L,3) = A3(L)
      ENDDO
         
C ... Obsolate block above, AAX could be replaced by a patch-sized array

C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH
C ... A patch type and indicator for solids
         IBC = ICON(1,IP)
         IF(IBC >= 8 .AND. IBC <= 10 .OR. IBC == 13 .OR.
     &      IBC >= 3 .AND. IBC <=  5 .OR. IBC == 15 .OR.
     &      IBC == 16.OR.  IBC == 17 .OR.
     &      IBC == 1 .AND. ICON(20,IP) ==  1        .OR. ! Data is stored
     &      IBC == 1 .AND. ICON(20,IP) ==  2        .OR. ! Data is stored?
     &      IBC == 1 .AND. ICON(20,IP) ==  3        .OR. ! Data is stored?
     &      IBC == 1 .AND. ICON(20,IP) ==  7        .OR. ! Data is stored
     &      IBC == 1 .AND. ICON(20,IP) ==  8)       THEN ! Data is stored

         IPL = ICON(2,IP) ! Proces local patch number
         IPG = ICON(25,IP)! Global patch number

         IF(IACTNO /= 0) THEN
            IACTU   = ICON(26,IP)! Actuator disc number. IACTU = IACTNO
            IACTY   = IFA(IACTNO)
         ELSE
            IACTU = 0
         ENDIF

         IDIR    = 1
         IFACE   = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1
         IMOV    = ICON(9,IP)

         IF(IBC == 8) ITYPE = 1
         IF(IBC == 9) NROT = NROT + 1
         IF(IBC == 9) ITYPE = 2
         IF(IBC == 10 .AND. IGRID == 31) NROT  = NROT + 1
         IF(IBC == 10 .AND. IMOV  ==  1) NMOV1 = NMOV1 + 1
         IF(IBC == 10 .AND. IMOV  ==  1) ITYPE = 3
         IF(IBC == 10 .AND. IMOV  ==  2) NMOV2 = NMOV2 + 1
         IF(IBC == 10 .AND. IMOV  ==  2) ITYPE = 4
         IF(IBC == 10 .AND. IMOV  ==  3) NMOV3 = NMOV3 + 1
         IF(IBC == 10 .AND. IMOV  ==  3) ITYPE = 5
         IF(IBC == 10 .AND. IMOV  ==  4) NMOV4 = NMOV4 + 1
         IF(IBC == 10 .AND. IMOV  ==  4) ITYPE = 6
         IF(IBC == 15) ITYPE = 10
         IF(IBC == 13) ITYPE = 11
         
C ... A global solid patch number in this process including the other blocks
        
         IF(M == 1) THEN
            NSOLPA = NSOLPA + 1
            ISP    = ISP + 1
         ENDIF

         XVEL = RCON(1,IP)
         YVEL = RCON(2,IP)
         ZVEL = RCON(3,IP)
         VEL  = RCON(4,IP)

         ARGSQR = XVEL**2+YVEL**2+ZVEL**2
         IF(ARGSQR > 0.0) THEN
             XNOR = SQRT(ARGSQR)
         ELSE
            XNOR = 1.E-10
         ENDIF 

         R1 = XVEL/XNOR
         R2 = YVEL/XNOR
         R3 = ZVEL/XNOR
             
C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         NORM = 1
         IF(IFACE == 4) K    = IMAX + 1
         IF(IFACE == 4) NORM = -1
         ID   = IDI1
         JST2 = JSTRID
         IAP  = 1

C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         NORM = 1
         IF(IFACE == 5) K    = JMAX + 1
         IF(IFACE == 5) NORM = -1
         ID   = IDI2
         JST2 = KSTRID
         IAP  = 2

C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         NORM = 1
         IF(IFACE == 6) K    = KMAX + 1
         IF(IFACE == 6) NORM = -1
         ID   = IDI3
         JST2 = ISTRID
         IAP  = 3
      ENDIF

C ... Store patch data to produce boundary files

         KX1 = ICON(4,IP)
         KX2 = ICON(5,IP)
         KY1 = ICON(6,IP)
         KY2 = ICON(7,IP)

      IF(M == 1 .AND. IALL == 0) THEN

         ISTRS(ISP) = ISTR
         JSTRS(ISP) = JSTR
         KSTRS(ISP) = KSTR
         NORMAL(ISP)= NORM

         KX1S(ISP)  = KX1
         KX2S(ISP)  = KX2
         KY1S(ISP)  = KY1
         KY2S(ISP)  = KY2
         KZ1S(ISP)  = K
         KZ2S(ISP)  = K
         IDIMS(ISP) = KX2 - KX1 + 1
         JDIMS(ISP) = KY2 - KY1 + 1

         IDIMSG(IP) = IDIMS(ISP)  ! Antakee mersu
         JDIMSG(IP) = JDIMS(ISP)
         NSPB(ISP)  = ICON(23,IP)
         NSPP(ISP)  = IPL
         NSPG(ISP)  = IPG

      ENDIF

      
      IDIMP = IDIMSG(IP)

      IF(IFACE == 2 .OR. IFACE == 5) THEN  ! Change the order
         KX1 = ICON(6,IP)
         KX2 = ICON(7,IP)
         KY1 = ICON(4,IP)
         KY2 = ICON(5,IP)
      ENDIF

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+K-1)*KSTR
      KBIND   = 0
      IF(IFACE >= 4) KBIND = -KSTR

C ... Calculate the centerpoint coordinates for all solid patches

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA 
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL)
         DO I = KX1,KX2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            L2 = 1 + (IN+I  )*ISTR + IJ
            L3 = 1 + (IN+I-1)*ISTR + IJ + JSTR
            L4 = 1 + (IN+I  )*ISTR + IJ + JSTR
            NP = I + NR                      ! Patch index with LN ghost cells

            ZCP(NP) = .25*(ZCO(L) + ZCO(L2) + ZCO(L3) + ZCO(L4)) 
            YCP(NP) = .25*(YCO(L) + YCO(L2) + YCO(L3) + YCO(L4))
            XCP(NP) = .25*(XCO(L) + XCO(L2) + XCO(L3) + XCO(L4))

          ENDDO
      ENDDO
      
C ... Reflect the cell centerpoint coordinates under the patches
      
      IF(M == 1) THEN   ! Centerpoints are on the first level

      APATCH(IPL) = 0.  ! Antakee mersu, onko nyt oikein?
      
      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL)
         DO I = KX1,KX2
C ...       Cell II index takes into account the top surface via KBIND
            II  = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II1 = II - LSTRID                ! UNDER THE SURFACE
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells

            XP1 = XC(II) - XCP(NP)
            YP1 = YC(II) - YCP(NP)
            ZP1 = ZC(II) - ZCP(NP)

C ... VEC{N} = (This is single precision)
            XN1 = AAX(L,1+(IAP-1)*3)
            XN2 = AAX(L,2+(IAP-1)*3)
            XN3 = AAX(L,3+(IAP-1)*3)
            XP0 = (1.-2.*XN1**2)*XP1-2.*XN1*XN2*YP1-2.*XN1*XN3*ZP1
            YP0 =-2.*XN1*XN2*XP1+(1.-2.*XN2**2)*YP1-2.*XN2*XN3*ZP1
            ZP0 =-2.*XN1*XN3*XP1-2.*XN2*XN3*YP1+(1.-2.*XN3**2)*ZP1
            XP0 = XP0 + XCP(NP)
            YP0 = YP0 + YCP(NP)
            ZP0 = ZP0 + ZCP(NP)
            XC(II1) = XP0
            YC(II1) = YP0
            ZC(II1) = ZP0          
            APATCH(IPL)= APATCH(IPL) + AAA(L,IAP)
         ENDDO
      ENDDO
      ENDIF

C*******************************************************************************

C ... Store the K-values and porosity for a pressure loss surface

      IF(IBC == 1 .AND. ICON(20,IP) == 7 .AND. IALL == 0) THEN  ! Pressure loss

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL)
         DO I = KX1,KX2
C           Cell II index takes into account the top surface via KBIND
            II  = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II1 = II - LSTRID                ! UNDER THE SURFACE
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            RMLOSS(NP) = RCON(1,IP)
            POROS(NP)  = RCON(2,IP)
c            write(671,*) I,J,NP,NR,IHF(IPL)
         ENDDO
      ENDDO

      ENDIF
      
C ******************************************************************************

      IF(IBC == 1 .AND. ICON(20,IP) == 8 .AND. IALL == 0
     +       .AND. IACTY >= 7 .AND. IACTY <= 10 ) THEN  ! An actuator disc

      DO J = KY1,KY2
          IJ = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL)
c          write(666,*) nr,nr-ihf(ipl),kx1+1,kx2-1,ky1+1,ky2-1, 'voim'
          DO I = KX1,KX2
C           Cell II index takes into account the top surface via KBIND
            II  = 1 + (IN+I-1)*ISTR + IJ + KBIND
            II1 = II - LSTRID                ! UNDER THE SURFACE
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            L2 = 1 + (IN+I  )*ISTR + IJ
            L3 = 1 + (IN+I-1)*ISTR + IJ + JSTR
            L4 = 1 + (IN+I  )*ISTR + IJ + JSTR
            NP = I + NR                      ! Patch index with LN ghost cells
            RMLOSS(NP) = RCON(1,IP)
            POROS(NP)  = RCON(2,IP)

            ZCP(NP) = .25*(ZCO(L) + ZCO(L2) + ZCO(L3) + ZCO(L4)) 
            YCP(NP) = .25*(YCO(L) + YCO(L2) + YCO(L3) + YCO(L4))
            XCP(NP) = .25*(XCO(L) + XCO(L2) + XCO(L3) + XCO(L4))
*            write(*,*) I,J,NP,NR,IHF(IPL),xcp(np),ycp(np),zcp(np)
          ENDDO
      ENDDO

C ... Find the distribution file for the actuator disc

      INQUIRE(FILE=BOUNDN(IPL),EXIST=THERE)
      IF(.NOT. THERE) THEN
         WRITE(*,*) 'Cannot find file:'
         WRITE(*,'(A50)') BOUNDN(IPL) ! Add IFA test here if an even distribution
         WRITE(*,*) NBL !,IBWALL
         WRITE(*,*) ' Exiting in BOUNDP...'
         STOP
      ENDIF

      OPEN(51,FILE=BOUNDN(IPL),STATUS='OLD',FORM='FORMATTED')
      DO III = 1,100000 ! Mersu
      READ(51,'(1A180)',END=2021,ERR=2023) LINE1
      LENGTH = LENGTH + 1
      ENDDO
2023  CONTINUE
      WRITE(*,*) 'Actuator disc went into a forest in SURPRO'
      WRITE(*,*) LINE1
      STOP
2021  CONTINUE
      REWIND(51) 

C ... Extract the distribution to the work arrays on each grid level (mersu)

      IF(IACTY >= 7 .AND. IACTY <= 10) THEN ! Read the distribution

         CALL READAD(LENGTH,IACTU,IACTY,XCP,YCP,ZCP,KA,KBIND,KX1,KY1,
     +        KX2,KY2,ISTR,JSTR,IHF(IPL),SURFA,SURFT,IPRO)

      ENDIF ! IACTY >= 7

      ENDIF ! IBC == 1 .AND. ICON(20,IP) == 8 (Actuator disc)

C ******************************************************************************

      IF(IBC >= 8 .AND. IBC <= 10 .OR. IBC == 15) THEN  ! A solid wall patch

C ... Nonvectorizable loop for the solid wall velocities
       
      VELAVE = 0.
      
*      DO 7000 J = KY1-1,KY2+1
      DO 7000 J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
*         DO I = KX1-1,KX2+1
         DO I = KX1,KX2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells

          IF(IBC == 9 .OR. (IBC == 10 .AND. IGRID >= 31
     +           .AND. IGRID <= 39)) THEN     ! A rotating solid

            ZAVE = ZCP(NP) - CENZ 
            YAVE = YCP(NP) - CENY
            XAVE = XCP(NP) - CENX
            UBIU = OMEGA*(OMEY*ZAVE - OMEZ*YAVE)
            UBIV = OMEGA*(OMEZ*XAVE - OMEX*ZAVE)
            UBIW = OMEGA*(OMEX*YAVE - OMEY*XAVE)
            UWALL(NP) = UBIU
            VWALL(NP) = UBIV
            WWALL(NP) = UBIW
            IF(J >= KY1 .AND. J <= KY2 .AND. I >= KX1 .AND. I <= KX2)
     +      VELAVE = VELAVE + SQRT(UBIU**2 + UBIV**2 + UBIW**2)

          ELSE IF(IBC == 10 .AND. IMOV == 1) THEN  ! A moving solid (1)

            UWALL(NP) = XVEL
            VWALL(NP) = YVEL
            WWALL(NP) = ZVEL
          
          ELSE IF(IBC == 10 .AND. IMOV == 2) THEN  ! A moving solid (2)

C ... VEC{N} =
            XN1       = AAX(L,1+(IAP-1)*3)
            XN2       = AAX(L,2+(IAP-1)*3)
            XN3       = AAX(L,3+(IAP-1)*3)
C ... VEC{R}xVEC{N}=
            XA1       = R2*XN3 - R3*XN2
            XA2       = R3*XN1 - R1*XN3
            XA3       = R1*XN2 - R2*XN1
            XNORM     = SQRT(XA1**2 + XA2**2 + XA3**2)
            IF(XNORM <= 1.E-6) THEN
               WRITE(*,*)'Unknown direction for a moving solid in',
     +         ' patch',IP
               WRITE(*,*)'A1,A2,A3:',XA1,XA2,XA3
               WRITE(*,*)'R1,R2,R3:',R1,R2,R3
               WRITE(*,*)'in place (i,j,k)',I,J,K
               WRITE(*,*)'Exiting ...'
               STOP
            ENDIF
            XA1       = XA1/XNORM
            XA2       = XA2/XNORM
            XA3       = XA3/XNORM
C ... VEC{V} = NORM{VEC{N}x(VEC{R}xVEC{N})}* Stymied. What is this?
            B1        = R2*XA3 - R3*XA2
            B2        = R3*XA1 - R1*XA3
            B3        = R1*XA2 - R2*XA1
            XNORM     = SQRT(B1**2 + B2**2 + B3**2)
            IF(ABS(XNORM-1.) >= 1.E-5) THEN
               WRITE(*,*) 'Something wrong in moving solid',
     +              ' in place (i,j,k)',I,J,K
               WRITE(*,*) 'Exiting ...'
               STOP
            ENDIF
            UWALL(NP)  = XA1*VEL
            VWALL(NP)  = XA2*VEL
            WWALL(NP)  = XA3*VEL
        ELSE IF(IBC == 10 .AND. IMOV == 3) THEN ! A moving solid (3)

C ... DIRECTION IS CALCULATED AS: 
C ... VEC{V} = NORM{(VEC{R}xVEC{X}}*VEL

            ZAVE      = ZCP(NP) - CENZ 
            YAVE      = YCP(NP) - CENY
            XAVE      = XCP(NP) - CENX
            UWALL(NP) = VEL*(R2*ZAVE - R3*YAVE)
            VWALL(NP) = VEL*(R3*XAVE - R1*ZAVE)
            WWALL(NP) = VEL*(R1*YAVE - R2*XAVE)
c            UWALL(NP) = VEL*(R1*YAVE - R2*XAVE)
c            VWALL(NP) = VEL*(R2*ZAVE - R3*YAVE)
c            WWALL(NP) = VEL*(R3*XAVE - R1*ZAVE)
            IF(J >= KY1 .AND. J <= KY2 .AND. I >= KX1 .AND. I <= KX2)
     +      VELAVE    = VELAVE + SQRT(UWALL(NP)**2+VWALL(NP)**2+
     +                  WWALL(NP)**2)

        ELSE IF(IBC == 10 .AND. IMOV == 4) THEN  ! A moving solid (4)

C ... A moving solid (1)
           
           IF (.NOT. TIMEL) THEN
              
              UWALL(NP) = GVEX
              VWALL(NP) = GVEY
              WWALL(NP) = GVEZ
              IF(J >= KY1 .AND. J <= KY2 .AND. I >= KX1 .AND. I <= KX2)
     +             VELAVE    = VELAVE + SQRT(UWALL(NP)**2+VWALL(NP)**2+
     +             WWALL(NP)**2)

           ENDIF
           
C ... Direction is calculated from the swept volume

           IF (TIMEL) THEN

            L2 = 1 + (IN+I  )*ISTR + IJ
            L3 = 1 + (IN+I-1)*ISTR + IJ + JSTR
            L4 = 1 + (IN+I  )*ISTR + IJ + JSTR

            ZCPLE2 = .25*(ZLE2(L)+ZLE2(L2)+ZLE2(L3)+ZLE2(L4)) 
            YCPLE2 = .25*(YLE2(L)+YLE2(L2)+YLE2(L3)+YLE2(L4))
            XCPLE2 = .25*(XLE2(L)+XLE2(L2)+XLE2(L3)+XLE2(L4))
            ZCPLE3 = .25*(ZLE3(L)+ZLE3(L2)+ZLE3(L3)+ZLE3(L4)) 
            YCPLE3 = .25*(YLE3(L)+YLE3(L2)+YLE3(L3)+YLE3(L4))
            XCPLE3 = .25*(XLE3(L)+XLE3(L2)+XLE3(L3)+XLE3(L4))

c           UWALL(NP) = (1.5*XCP(NP)-2.*XCPLE2+.5*XCPLE3)/DT
c           VWALL(NP) = (1.5*YCP(NP)-2.*YCPLE2+.5*YCPLE3)/DT
c           WWALL(NP) = (1.5*ZCP(NP)-2.*ZCPLE2+.5*ZCPLE3)/DT
*     write(*,*) UWALL(NP),VWALL(NP),WWALL(NP)

            IF(ABS(XCPLE2 - XCPLE3) <= EPS6) THEN
               UWALL(NP) = (XCP(NP) - XCPLE2)/DT
            ELSE
               UWALL(NP) = 1.5*(XCP(NP) - XCPLE2)/DT -
     +                   .5*(XCPLE2  - XCPLE3)/DTOLD
            ENDIF
            IF(ABS(YCPLE2 - YCPLE3) <= EPS6) THEN
               VWALL(NP) = (YCP(NP) - YCPLE2)/DT 
            ELSE
               VWALL(NP) = 1.5*(YCP(NP) - YCPLE2)/DT -
     +                   .5*(YCPLE2  - YCPLE3)/DTOLD
            ENDIF
            IF(ABS(ZCPLE2 - ZCPLE3) <= EPS6) THEN
               WWALL(NP) = (ZCP(NP) - ZCPLE2)/DT 
            ELSE
               WWALL(NP) = 1.5*(ZCP(NP) - ZCPLE2)/DT -
     +                   .5*(ZCPLE2  - ZCPLE3)/DTOLD
            ENDIF

            IF(J >= KY1 .AND. J <= KY2 .AND. I >= KX1 .AND. I <= KX2)
     +      VELAVE    = VELAVE + SQRT(UWALL(NP)**2+VWALL(NP)**2+
     +                  WWALL(NP)**2)
            
            ENDIF ! IF (TIMEL) THEN

        ENDIF  ! IBC == 9 (Solid wall)

      ENDDO
7000  CONTINUE  ! End of non-vectorizable patch loop

C *****************************************************************

C ... Print to the MEMORY file 
  
      IF(IBC == 9 .OR. (IBC == 10 .AND. IGRID == 31)) THEN  ! A rotating solid

      VELAVE = VELAVE/REAL((KY2-KY1+1)*(KX2-KX1+1))
      IF(FIRSTL) THEN
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9201) NBL,M
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9202)
         FIRSTL = .FALSE.
      ENDIF
      IF(IREPEA(3) <= NREPEA(3))
     +WRITE(45,*)'Average velocity of the rotating patch'
     +       ,IP,' is',VELAVE
c      WRITE(45,*)

      ELSE IF(IBC == 10 .AND. IMOV == 1) THEN  ! A moving solid (1)

      IF(FIRSTL) THEN
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9201) NBL,M
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9202)
         FIRSTL = .FALSE.
      ENDIF
c9201  FORMAT(' Block = ',I3,' MG = ',I1)
9201  FORMAT(' Block = ',I6,' MG = ',I1)
9202  FORMAT(' ------------------')
      IF(IREPEA(3) <= NREPEA(3))
     +WRITE(45,101) IP,XVEL,YVEL,ZVEL,SQRT(XVEL**2+YVEL**2+ZVEL**2)
 101  FORMAT(/'Patch' ,I4,' is moving solid with velocity:',3F10.5,'.'/
     +     'Total velocity is ',F10.5/)
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,*)

      ELSE IF(IBC == 10 .AND. IMOV == 2) THEN ! A moving solid (2)

      IF(FIRSTL) THEN
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9201) NBL,M
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9202)
         FIRSTL = .FALSE.
      ENDIF
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,111) IP,R1,R2,R3,VEL
 111  FORMAT(/'Patch' ,I4,' is a tangentially moving solid.',
     +     ' Direction of solid is:',3F10.5,'.'/
     +     'Total velocity is ',F10.5/)
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,*)

      ELSE IF(IBC == 10 .AND. IMOV == 3) THEN  ! A moving solid (3)

      IF(FIRSTL) THEN
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9201) NBL,M
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9202)
         FIRSTL = .FALSE.
      ENDIF
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,113) IP,R1,R2,R3,VEL
 113  FORMAT(/'Patch' ,I4,' is a rotationally moving solid.'/
     +     ' Rotational axis is:',3F8.3,'.'/
     +     'Angular velocity is ',F10.5,' 1/rad')
      VELAVE = VELAVE/REAL((KY2-KY1+1)*(KX2-KX1+1))
      IF(IREPEA(3) <= NREPEA(3))
     +WRITE(45,*)'Average velocity of the moving patch',IP,' is',VELAVE
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,*)

      ELSE IF(IBC == 10 .AND. IMOV == 4) THEN  ! A moving solid (4)

      IF(FIRSTL) THEN
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9201) NBL,M
         IF(IREPEA(3) <= NREPEA(3)) WRITE(45,9202)
         FIRSTL = .FALSE.
      ENDIF
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,114) IP
 114  FORMAT(/'Patch' ,I4,' is a freely moving solid.'/)
      VELAVE = VELAVE/REAL((KY2-KY1+1)*(KX2-KX1+1))
      IF(IREPEA(3) <= NREPEA(3))
     +WRITE(45,*)'Average velocity of the moving patch',IP,' is',VELAVE
      IF(IREPEA(3) <= NREPEA(3)) WRITE(45,*)

      ELSE IF(IBC /= 8 .AND. IBC /=15) THEN

         WRITE(*,*) 'Error in patch',IP,' Type',IBC
         WRITE(*,*) 'Non existing moving solid type:',IMOV
         WRITE(*,*) 'Exiting ...'
         STOP

      ENDIF     ! IBC == 9 (Print to the memory file)

C *****************************************************************
      ENDIF     ! IBC >= 8 .AND. IBC <= 10
      ENDIF     ! IBC >= 8 .AND. IBC <= 10 .OR. IBC == 3 .OR. IBC == 5 .OR.
                ! IBC == 15 .OR. IBC == 16 .OR. IBC == 17

8000  CONTINUE

      IF(NSOLPA > 0 .AND. IREPEA(3) <= NREPEA(3)) THEN  ! Report solid patches
         WRITE(45,*)
         WRITE(45,9100) NBL,M,NSOLPA
         WRITE(45,9106)
         WRITE(45,9101) NSOLPA - NROT - NMOV1 - NMOV2 - NMOV3
         WRITE(45,9102) NROT
         WRITE(45,9103) NMOV1
         WRITE(45,9104) NMOV2
         WRITE(45,9105) NMOV3
         WRITE(45,*)
c9100     FORMAT(1X,'Block = ',I3,' MG = ',I1,' found',I4,
c     +   ' solid patches')
9100     FORMAT(1X,'Block = ',I6,' MG = ',I1,' found',I7,
     +   ' solid patches')
9106     FORMAT(' ------------------------------------------')
9101     FORMAT(' Stationary = ',I4)
9102     FORMAT(' Rotating   = ',I4)
9103     FORMAT(' Moving(1)  = ',I4)
9104     FORMAT(' Moving(2)  = ',I4)
9105     FORMAT(' Moving(3)  = ',I4)
      ENDIF ! NSOLPA > 0

      DEALLOCATE(AAX,AAA)
       
      RETURN
      END SUBROUTINE SURPRO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SURHTP(UWALL,TWALL,HFLUX,WMFLUX,POROS,
     +     WHSTAG,WTEMP,RSDIRX,RSDIRY,RSDIRZ,RBK,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,
     +     INRE,JNRE,KNRE,NBL,M,NPATCH,ICON,RCON,
     +     IHF,NSOLPA,ISTRS,JSTRS,KSTRS,NSPB,ISTATE,MAXB,
     +     IREPEA,NREPEA,WTRAN,ngl)

      USE NS3CO, ONLY : LN, IC9

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +     IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     +     IDI1,IDI2,IDI3,NBL,M,NPATCH,ISTATE,MAXB,ISTRID,JSTRID,KSTRID,
     +     IL,IP,IBC,IPL,IHEAT,IDIR,IFACE,ISTR,JSTR,KSTR,I,J,K,ID,NFL,
     +     KX1,KX2,KY1,KY2,LSTRID,KA,IJ,L,JST2,
     +     IAP,NR,NP,NSOLPA,KBIND,IF1,IPG,ngl
      INTEGER :: ICON(IC9,*),IHF(*),IREPEA(*),NREPEA(*)
      INTEGER :: ISTRS(*),JSTRS(*),KSTRS(*),NSPB(*)

      REAL :: RCON(IC9,*),UWALL(*),TWALL(*),
     +     HFLUX(*),WMFLUX(*),POROS(*),A1X(*),A1Y(*),
     +     A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +     WHSTAG(*),WTEMP(*),RSDIRX(*),RSDIRY(*),RSDIRZ(*),RBK(*),
     +     WTRAN(*)

      REAL :: XN1,XN2,XN3,RSLEN

      REAL, ALLOCATABLE :: AAX(:,:)

      LOGICAL :: THERE

C ... Definitions to be moved to file reading subroutine

      INTEGER :: IPGLO,IMAX1,JMAX1,NTOT
          
      CHARACTER(LEN=3) :: FILEE
     
C ... Calculate heat transfer properties for the solid walls. The patches
C ... are extended by one row of cells for postprocessing

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
      IL      = ISTRID*JSTRID

      IF(M == 1 .AND. NSOLPA > 0 .AND. IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*)
         WRITE(45,9200) NBL,M
         WRITE(45,9201)
      ENDIF
c 9200 FORMAT(' Block = ',I3,' MG = ',I1)
 9200 FORMAT(' Block = ',I6,' MG = ',I1)
 9201 FORMAT(' ------------------')

C ... Antakee mersu

      ALLOCATE(AAX(NTOT,9))
      DO L = 1,NTOT
         AAX(L,1) = A1X(L)
         AAX(L,2) = A1Y(L)
         AAX(L,3) = A1Z(L)
         AAX(L,4) = A2X(L)
         AAX(L,5) = A2Y(L)
         AAX(L,6) = A2Z(L)
         AAX(L,7) = A3X(L)
         AAX(L,8) = A3Y(L)
         AAX(L,9) = A3Z(L)
      ENDDO

C ... Obsolate block above, AAX could be replaced by a patch-sized array

C ... LOOP OVER THE PATCHES OF THE BLOCK 

C ... heat flux calculations (The same as the previous subroutine SURVAR)

      DO 8000 IP = 1,NPATCH

      IBC   = ICON(1,IP)
      IPL   = ICON(2,IP) ! proces local patch number
      IHEAT = ICON(20,IP)
      IPG   = ICON(25,IP)! Global patch number
      IDIR   = 1
       
      IF(IBC >= 8 .AND. IBC <= 10 .OR. IBC == 15) THEN

      IFACE   = ICON(3,IP)
      IF(IFACE >= 4) IDIR = -1

C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = IBOT
         ID   = IDI1
         JST2 = JSTRID
         IAP  = 1
         IF(IFACE == 4) K   = ITOP + 1

C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = JBOT
         ID   = IDI2
         JST2 = KSTRID
         IAP  = 2
         IF(IFACE == 5) K   = JTOP + 1

C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K = KBOT
         ID   = IDI3
         JST2 = ISTRID
         IAP  = 3
         IF(IFACE == 6) K   = KTOP + 1
      ENDIF

      KX1 = ICON(4,IP)
      KX2 = ICON(5,IP)
      KY1 = ICON(6,IP)
      KY2 = ICON(7,IP)

      IF(ICON(3,IP) == 2 .OR. ICON(3,IP) == 5) THEN
         KY1 = ICON(4,IP)
         KY2 = ICON(5,IP)
         KX1 = ICON(6,IP)
         KX2 = ICON(7,IP)
      ENDIF

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA      = (KN+K-1)*KSTR
      KBIND   = 0
      IF(IFACE >= 4) KBIND = - KSTR

      IPGLO = ICON(25,IP)                     ! global patch number
      IMAX1 = KX2-KX1+1
      JMAX1 = KY2-KY1+1
      IF1   = IHF(IPL)

C ... Give surface roughness

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            RBK(NP) = RCON(13,IP)
         ENDDO
      ENDDO

C ... Type of boundary layer

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            WTRAN(NP) = RCON(14,IP)
         ENDDO                     
      ENDDO

C ... Specify injected mass flux

      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            WMFLUX(NP) = RCON(7,IP)
         ENDDO
      ENDDO

C ... Specify wall porosity

      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            POROS(NP) = RCON(8,IP)
         ENDDO
       ENDDO

C ... Specify a possible injection vector

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX UP (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
C ... VEC{N} =
            XN1 = AAX(L,1+(IAP-1)*3)
            XN2 = AAX(L,2+(IAP-1)*3)
            XN3 = AAX(L,3+(IAP-1)*3)
            IF(RCON(10,IP) == 0. .AND. RCON(11,IP) == 0. .AND.
     +      RCON(12,IP) == 0.) THEN ! Normal direction
               RSDIRX(NP) = XN1
               RSDIRY(NP) = XN2
               RSDIRZ(NP) = XN3
            ELSE                    ! Tangential direction
               RSDIRX(NP) = XN2*RCON(12,IP) - XN3*RCON(11,IP)
               RSDIRY(NP) = XN3*RCON(10,IP) - XN1*RCON(12,IP)
               RSDIRZ(NP) = XN1*RCON(11,IP) - XN2*RCON(10,IP)
               RSLEN      = SQRT(RSDIRX(NP)**2 + RSDIRY(NP)**2
     +                    +      RSDIRZ(NP)**2) + 1.E-20
               RSDIRX(NP) = RSDIRX(NP)/RSLEN
               RSDIRY(NP) = RSDIRY(NP)/RSLEN
               RSDIRZ(NP) = RSDIRZ(NP)/RSLEN
            ENDIF
         ENDDO
      ENDDO

C ... Specify the injection temperature ? 

      DO J = KY1,KY2
      NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
c            WTEMP(NP) = RCON(14,IP)         ! Unused
         ENDDO
      ENDDO

C ... Specify the injection stagnation enthalpy

      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            WHSTAG(NP) = RCON(9,IP)
         ENDDO
      ENDDO


C ***********************************************************************

C ... Specify the heat transfer type

      IF(IHEAT == 0) THEN ! An adiabatic surface is assumed

         IF(IREPEA(3) <= NREPEA(3))
     &   WRITE(45,*)'Adiabatic surface is applied for patch =',IP 

      ELSE IF(IHEAT == 1) THEN ! Surface temperature is given

      IF(IREPEA(3) <= NREPEA(3))
     &WRITE(45,*)'For patch =',IP,'the surface temperature is =',
     &RCON(5,IP) 
      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            TWALL(NP) = RCON(5,IP)
         ENDDO
      ENDDO

      ELSE IF(IHEAT == 2) THEN ! Heat flux is given

      IF(IREPEA(3) <= NREPEA(3))
     &WRITE(45,*)'For patch =',IP,'the surface heat flux is =',
     &RCON(6,IP) 
      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1,KX2
            NP = I + NR                      ! Patch index with LN ghost cells
            HFLUX(NP) = RCON(6,IP)
         ENDDO
      ENDDO

      ELSE IF(IHEAT == 3) THEN ! Stagnation temperature assumed 

      IF(IREPEA(3) <= NREPEA(3))
     &WRITE(45,*)'Stagnation temperature assumed for patch =',IP 

      ELSE IF(IHEAT == 4) THEN ! Coupling between solid and fluid 

      IF(IREPEA(3) <= NREPEA(3))
     &WRITE(45,*)'Solid-fluid coupling assumed for patch =',IP 

C ... The following options are not tested!!! 

      ELSE IF(IHEAT == 11) THEN  ! Read the wall temperature from a file

       CALL NUMCH3(FILEE,IPGLO)
       OPEN(51,FILE='THERMO'//FILEE,STATUS='OLD',FORM='FORMATTED')
       CALL READHT(TWALL(IF1),IMAX,JMAX,IMAX1,JMAX1)      
       CLOSE(51)

      ELSE IF(IHEAT == 12) THEN ! Read the wall heat flux from a file

       CALL NUMCH3(FILEE,IPGLO)
       OPEN(51,FILE='THERMO'//FILEE,STATUS='OLD',FORM='FORMATTED')
       CALL READHT(HFLUX(IF1),IMAX,JMAX,IMAX1,JMAX1)      
       CLOSE(51)

      ELSE IF(IHEAT == 7) THEN ! Read the porosities from a file


      ELSE IF(IHEAT > 0) THEN  ! Unknown boundary type 
                    
         WRITE(*,*) ' No such solid tempperature boundary condition'
         WRITE(*,*) ' Exiting ...'
         STOP
      ENDIF                     !(IHEAT == 1)

C ***********************************************************************

C ... Check whether a transition file exists for this patch

      CALL NUMCH3(FILEE,IPGLO)
      INQUIRE(FILE='WTRAN'//FILEE,EXIST=THERE)

      IF(THERE) THEN
      WRITE(45,9202) IP,RCON(13,IP)
 9202 FORMAT((' Surface roughness of patch ='),
     &     1X,I6,(' is read from a file instead of RBK ='),E11.4)
         OPEN(51,FILE='RBK'//FILEE,STATUS='OLD',FORM='FORMATTED')
         CALL READHT(RBK(IF1),IMAX,JMAX,IMAX1,JMAX1)      
         CLOSE(51)
      ENDIF ! THERE

C ... Check whether a surface-roughness file exists for this patch

      CALL NUMCH3(FILEE,IPGLO)
      INQUIRE(FILE='RBK'//FILEE,EXIST=THERE)

      IF(THERE) THEN
      WRITE(45,9302) IP,RCON(14,IP)
 9302 FORMAT((' Flow type of patch (0=laminar, 1=fully turbulent) ='),
     &     1X,I6,(' is '),E11.4)
         OPEN(51,FILE='WTRAN'//FILEE,STATUS='OLD',FORM='FORMATTED')
         CALL READHT(WTRAN(IF1),IMAX,JMAX,IMAX1,JMAX1)      
         CLOSE(51)
      ENDIF ! THERE

C ... Extend patches into the ghost cells

      CALL EXTPAT(WTRAN(IF1),KX1,KX2,KY1,KY2)
      CALL EXTPAT(  RBK(IF1),KX1,KX2,KY1,KY2)
      CALL EXTPAT(TWALL(IF1),KX1,KX2,KY1,KY2)
      CALL EXTPAT(WTEMP(IF1),KX1,KX2,KY1,KY2)

C *****************************************************************

C ... Print to the MEMORY file

      IF(M == 1 .AND. IREPEA(3) <= NREPEA(3)) THEN  !   ! Printing for level 1

      WRITE(45,9203) IP,RCON(13,IP)
9203  FORMAT((' Surface roughness of patch'),6X,I6,(' is '),E10.3)

      IF(RCON(7,IP) /= 0.) THEN
         WRITE(45,9207) IP,RCON(7,IP)
9207     FORMAT((' Injected mass flux'),15X,I6,(' is '),E11.4,
     +       (' kg/m2s'))
         WRITE(45,9208) IP,RCON(7,IP)
9208     FORMAT((' Wall porosity'),20X,I6,(' is '),E11.4)
      ENDIF

      IF(IHEAT == 1) THEN                      ! Specified temperature
         WRITE(45,9204) IP,RCON(5,IP)
9204     FORMAT(('Surface temperature of patch'),5X,I6,(' is '),
     +   F8.3,(' K'))
      ELSE IF(IHEAT == 2) THEN                 ! Specified heat flux
         WRITE(45,9205) IP,RCON(6,IP)
9205     FORMAT(('Surface heat flux of patch'),7X,I6,(' is '),
     +   E11.4,(' W/m2'))
      ELSE 
         WRITE(45,9206)IP
9206     FORMAT((' No heat transfer data for patch'),1X,I6)
      ENDIF     ! IHEAT == 1 (Print to the memory file)

      WRITE(45,*)
      ENDIF     ! M == 1

C *****************************************************************
      ENDIF     ! IBC >= 8 .AND. IBC <= 10
8000  CONTINUE

      DEALLOCATE(AAX)
      RETURN
      END SUBROUTINE SURHTP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE SURTRA(VTRAN,WTRAN,NPATCH,ICON,IHF,IMAX,JMAX,KMAX,
     +     INRE,JNRE,KNRE)

      USE NS3CO, ONLY : LN, IC9
      
      IMPLICIT NONE

      INTEGER :: IP,NPATCH,IBC,IPL,IFACE,ISTR,JSTR,KSTR,ISTRID,JSTRID,
     +     KSTRID,IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,I,J,K,IL,KBIND,
     +     KX1,KX2,KY1,KY2,IF1,KA,II,IJ,NP,NR,IPG

      INTEGER :: ICON(IC9,*), IHF(*)

      REAL :: VTRAN(*), WTRAN(*)

C ... LOOP OVER THE PATCHES OF THE BLOCK 

C ... Update transition information inside a block

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      DO 8000 IP = 1,NPATCH

      IBC = ICON(1,IP)
      IPL = ICON(2,IP) ! proces local patch number
      IPG = ICON(25,IP)! Global patch number

      IF(IBC >= 8 .AND. IBC <= 10 .OR. IBC == 15) THEN

      IFACE   = ICON(3,IP)

C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         KBIND= IMAX

C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         KBIND= JMAX

C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         KBIND= KMAX
      ENDIF

      KX1 = ICON(4,IP)
      KX2 = ICON(5,IP)
      KY1 = ICON(6,IP)
      KY2 = ICON(7,IP)

      IF(ICON(3,IP) == 2 .OR. ICON(3,IP) == 5) THEN
         KY1 = ICON(4,IP)
         KY2 = ICON(5,IP)
         KX1 = ICON(6,IP)
         KX2 = ICON(7,IP)
      ENDIF

      IF1   = IHF(IPL)

C ... Update transition data

      DO K = 0,KBIND+1

      KA = (KN+K-1)*KSTR

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL)
         DO I = KX1,KX2
            II = 1 + (IN+I-1)*ISTR + IJ
            NP = I + NR                      ! Patch index with LN ghost cells
            VTRAN(II) = WTRAN(NP)
         ENDDO
      ENDDO

      ENDDO

      ENDIF    !  IBC >= 8.AND.IBC <= 10 .OR. IBC == 15

8000  CONTINUE

      RETURN
      END SUBROUTINE SURTRA
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE BOUNDP(UWALL,VWALL,WWALL,TWALL,HFLUX,WMFLUX,POROS,
     +     WHSTAG,WTEMP,RBK,XCP,YCP,ZCP,
     +     XCO,YCO,ZCO,VOL,A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     +     A3X,A3Y,A3Z,A1,A2,A3,UROT,VROT,WROT,
     +     IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,
     +     KTOP,IDI1,IDI2,IDI3,INRE,JNRE,KNRE,NGL,NBL,M,NPATCH,
     +     ICON,RCON,IHF,NSPTOT,NSOLPA,NORMAL,IDIMS,JDIMS,
     +     KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,ISTRS,JSTRS,KSTRS,
     +     NSPB,NSPP,ISTATE,JSTATE,NBLOCG,
     +     MAXB,XC,YC,ZC,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNMF,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,FRSVEL,ALPHA,BETA,FRSDEN,FRSSIE,FRSTEM,
     +     FRSPRE,FRSRK,FRSEPS,ITURB,INITCC,BOUNDF,BOUNPD,BOUNA1,BOUNA2,
     +     BOUNFI,BOUNBI,LEVEL,NSCAL,MAXEB,MAXSB,IBF,ISTRES,RGAS,GAMMA,
     +     E0REF,T0REF,MULPHL,NPHASE,FRSALFA,APATCH,IREPEA,NREPEA,
     +     BOUNG,BOUNRET,TRANSL,FRSTUR,FRSVIS,IUPPTL,
     +     BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2,FRADEN)

      USE NS3CO, ONLY : LN, IC9, GX, GY, GZ, GROUND, KSCAL

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,INRE,JNRE,KNRE,
     +     IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     +     IDI1,IDI2,IDI3,NBL,M,NPATCH,ISTATE,MAXB,ISTRID,JSTRID,KSTRID,
     +     IL,IP,IBC,IPL,IHEAT,IDIR,IFACE,ISTR,JSTR,KSTR,I,J,K,ID,NFL,
     +     KX1,KX2,KY1,KY2,LSTRID,KA,IJ,L,JST2,
     +     IAP,NR,NP,NSPTOT,NSOLPA,
     +     KBIND,IF1,INITCC,IMAX0,JMAX0,LEVEL,NSCAL,
     +     IFILE,MAXEB,MAXSB,IBF,NBLOCG,NPHASE,NGL,IPHASE,ICHAN,IUPPTL
      INTEGER :: ICON(IC9,*),IHF(*),IREPEA(*),NREPEA(*)
      INTEGER :: NORMAL(*),IDIMS(*),JDIMS(*),KX1S(*),KX2S(*),KY1S(*),
     +     KY2S(*),KZ1S(*),KZ2S(*),ISTRS(*),JSTRS(*),KSTRS(*),NSPB(*),
     +     NSPP(*),JSTATE(NBLOCG,*)

      REAL :: YCO(*),ZCO(*),XCO(*),XC(*),YC(*),ZC(*),
     +                    XCP(*),YCP(*),ZCP(*)

      REAL :: RCON(IC9,*),VOL(*),UWALL(*),VWALL(*),WWALL(*),UROT(*),
     +     VROT(*),WROT(*),TWALL(*),HFLUX(*),
     +     WMFLUX(*),POROS(*),A1X(*),A1Y(*),
     +     A1Z(*),A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),
     +     WHSTAG(*),WTEMP(*),RBK(*),
     +     BOUNR(*),BOUNU(*),BOUNV(*),BOUNW(*),BOUNE(*),BOUNMF(*),
     +     BOUNT(*),BOUNP(*),BOUNRK(*),BOUNEP(*),BOUNPD(*),BOUNA1(*),
     +     BOUNFI(IBF,*),BOUNBI(IBF,*),BOUNA2(*),A1(*),A2(*),A3(*),
     +     APATCH(*),BOUNG(*),BOUNRET(*),BOUNU1(*),BOUNU2(*),BOUNV1(*),
     +     BOUNV2(*),BOUNW1(*),BOUNW2(*),FRADEN(*)
      REAL :: 
     +     XN1,XN2,XN3,
     +     U10,V10,W10,UU10,VV10,WW10,FRSVEL,ALPHA,BETA,FRSDEN,FRSSIE,
     +     FRSTEM,FRSPRE,FRSRK,FRSEPS,RGAS,GAMMA,E0REF,T0REF,E1(3),
     +     UK,VK,WK,FRSALFA(3),RO1(3),RETIN,FRSTUR,RMUL1,RMUL2,
     +     TEMP,P,FRSVIS

      REAL, ALLOCATABLE :: AAX(:,:),UUROT(:,:),AAA(:,:), APU(:)

      REAL :: BIJ(5), FISET(10)

C ... Definitions to be moved to file reading subroutine

      INTEGER :: IPGLO,IMAX1,JMAX1,
     +           NTOT,ITURB,IBTYPE,ISTRES,IPG

      INTEGER :: INPTYP
      REAL    :: XMASS, HTOT, PTOT, TURBLE, RMUINI, LRATUP, TG, TRET
          
      CHARACTER(LEN=80) :: BOUNDF(*)

      LOGICAL :: THERE, MULPHL, TRANSL
     
C ... Calculate boundary conditions for the inlets and outlets. The patches
C ... are extended by one row of cells for postprocessing

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID  = KMAX + 2*KN
      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      NTOT    = ISTRID*JSTRID*KSTRID
      IL      = ISTRID*JSTRID

      IF(M == 1 .AND. NSOLPA > 0 .AND. IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*)
         WRITE(45,9201) NBL,M
         WRITE(45,9202)
      ENDIF
c9201  FORMAT(' Block = ',I3,' MG = ',I1)
9201  FORMAT(' Block = ',I6,' MG = ',I1)
9202  FORMAT(' ------------------')

C ... Antakee 2 mersu
      ALLOCATE(AAX(NTOT,9))
      ALLOCATE(UUROT(NTOT,3))
      ALLOCATE(AAA(NTOT,3))
      ALLOCATE(APU(NTOT))

      DO L = 1,NTOT
         AAX(L,1) = A1X(L)
         AAX(L,2) = A1Y(L)
         AAX(L,3) = A1Z(L)
         AAX(L,4) = A2X(L)
         AAX(L,5) = A2Y(L)
         AAX(L,6) = A2Z(L)
         AAX(L,7) = A3X(L)
         AAX(L,8) = A3Y(L)
         AAX(L,9) = A3Z(L)
         AAA(L,1) = A1(L)
         AAA(L,2) = A2(L)
         AAA(L,3) = A3(L)
         UUROT(L,1)= UROT(L)
         UUROT(L,2)= VROT(L)
         UUROT(L,3)= WROT(L)
      ENDDO

C ... Obsolate block above, AAX could be replaced by a patch-sized array

C ... LOOP OVER THE PATCHES OF THE BLOCK 

C ... Onlet/outlet is read in

      DO 8000 IP = 1,NPATCH

      IBC    = ICON(1,IP)
      IPL    = ICON(2,IP)   ! proces local patch number
      IPG    = ICON(25,IP)  ! Process global block number
      IBTYPE = ICON(27,IP)  ! Boundary condition type for inlets/outlets
      IHEAT  = ICON(20,IP)  ! Heat transfer type
      IDIR   = 1

      IF(IBC == 3 .OR. IBC == 5) THEN  ! An inlet or an outlet

      IFACE   = ICON(3,IP)
      IF(IFACE >= 4) IDIR = -1

C ... XI-DIRECTION
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         ID   = IDI1
         JST2 = JSTRID
         IAP  = 1
         IF(IFACE == 4) K = IMAX + 1

C ... ETA-DIRECTION
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         ID   = IDI2
         JST2 = KSTRID
         IAP  = 2
         IF(IFACE == 5) K = JMAX + 1

C ... ZETA DIRECTION
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         ID   = IDI3
         JST2 = ISTRID
         IAP  = 3
         IF(IFACE == 6) K = KMAX + 1
      ENDIF

      KX1 = ICON(4,IP)
      KX2 = ICON(5,IP)
      KY1 = ICON(6,IP)
      KY2 = ICON(7,IP)

      IF(ICON(3,IP) == 2 .OR. ICON(3,IP) == 5) THEN
         KY1 = ICON(4,IP)
         KY2 = ICON(5,IP)
         KX1 = ICON(6,IP)
         KX2 = ICON(7,IP)
      ENDIF

      LSTRID  = KSTR*IDIR
      IF(IDIR ==  1) NFL = 0
      IF(IDIR == -1) NFL = -LSTRID
      KA    = (KN+K-1)*KSTR
      KBIND = 0
      IF(IFACE >= 4) KBIND = - KSTR

      IPGLO = ICON(25,IP)                     ! global patch number
      IMAX1 = KX2-KX1+1
      JMAX1 = KY2-KY1+1
      IMAX0 = 2**(LEVEL-1)*2**(M-1)*IMAX1     ! size on the densest level
      JMAX0 = 2**(LEVEL-1)*2**(M-1)*JMAX1
      IF1   = IHF(IPL)


      IF(IBTYPE > 20 .AND. IBTYPE <= 30) THEN ! Read INL/OUT parameters

      INQUIRE(FILE=BOUNDF(IPL),EXIST=THERE)
      IF(.NOT. THERE) THEN
         WRITE(*,*) 'Cannot find file:'
         WRITE(*,'(A50)') BOUNDF(IPL)
         WRITE(*,*) NBL
         WRITE(*,*) ' Exiting in BOUNDP...'
         STOP
      ENDIF

      IFILE = 0  ! Primitive variables
      ICHAN = 51

      OPEN(ICHAN,FILE=BOUNDF(IPL),STATUS='OLD',FORM='FORMATTED')

C ... Read INL/OUT flow condition parameters (Two-fluid is missing)

      CALL READIOP(ICHAN,XMASS,HTOT,P,PTOT,TEMP,BIJ,
     &             TURBLE,RMUINI,RMUL1,RMUL2,TG,TRET,
     &             FISET,INPTYP,TRANSL,KSCAL,LRATUP)

C ... Calculate the boundary values from user given parameters

      IF(IFACE == 1 .OR. IFACE == 4) THEN

         CALL CALCIOP(BOUNMF(IF1),BOUNR(IF1),BOUNU(IF1),BOUNV(IF1),
     &        BOUNW(IF1),BOUNE(IF1),BOUNT(IF1),BOUNP(IF1),BOUNRK(IF1),
     &        BOUNEP(IF1),BOUNFI(IF1,1),BOUNBI(IF1,1),
     &        BOUNA1(IF1),BOUNA2(IF1),BOUNG(IF1),BOUNRET(IF1),
     &        MAXEB,MAXSB,IF1,KSCAL,ITURB,TRANSL,MULPHL,
     &        XCO,YCO,ZCO,UROT,KX1,KX2,KY1,KY2,
     &        A1,A1X,A1Y,A1Z,IMAX,JMAX,KMAX,NGL,IFACE,
     &        XMASS,HTOT,PTOT,P,TEMP,TURBLE,RMUINI,RMUL1,RMUL2,
     &        TG,TRET,LRATUP,IFILE,INPTYP,APU,FRSPRE,FRSTEM,FRSVEL,
     &        FRSDEN,FRSSIE,FRSVIS,IUPPTL)

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

         CALL CALCIOP(BOUNMF(IF1),BOUNR(IF1),BOUNU(IF1),BOUNV(IF1),
     &        BOUNW(IF1),BOUNE(IF1),BOUNT(IF1),BOUNP(IF1),BOUNRK(IF1),
     &        BOUNEP(IF1),BOUNFI(IF1,1),BOUNBI(IF1,1),
     &        BOUNA1(IF1),BOUNA2(IF1),BOUNG(IF1),BOUNRET(IF1),
     &        MAXEB,MAXSB,IF1,KSCAL,ITURB,TRANSL,MULPHL,
     &        XCO,YCO,ZCO,VROT,KX1,KX2,KY1,KY2,
     &        A2,A2X,A2Y,A2Z,IMAX,JMAX,KMAX,NGL,IFACE,
     &        XMASS,HTOT,PTOT,P,TEMP,TURBLE,RMUINI,RMUL1,RMUL2,
     &        TG,TRET,LRATUP,IFILE,INPTYP,APU,FRSPRE,FRSTEM,FRSVEL,
     &        FRSDEN,FRSSIE,FRSVIS,IUPPTL)

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         CALL CALCIOP(BOUNMF(IF1),BOUNR(IF1),BOUNU(IF1),BOUNV(IF1),
     &        BOUNW(IF1),BOUNE(IF1),BOUNT(IF1),BOUNP(IF1),BOUNRK(IF1),
     &        BOUNEP(IF1),BOUNFI(IF1,1),BOUNBI(IF1,1),
     &        BOUNA1(IF1),BOUNA2(IF1),BOUNG(IF1),BOUNRET(IF1),
     &        MAXEB,MAXSB,IF1,KSCAL,ITURB,TRANSL,MULPHL,
     &        XCO,YCO,ZCO,WROT,KX1,KX2,KY1,KY2,
     &        A3,A3X,A3Y,A3Z,IMAX,JMAX,KMAX,NGL,IFACE,
     &        XMASS,HTOT,PTOT,P,TEMP,TURBLE,RMUINI,RMUL1,RMUL2,
     &        TG,TRET,LRATUP,IFILE,INPTYP,APU,FRSPRE,FRSTEM,FRSVEL,
     &        FRSDEN,FRSSIE,FRSVIS,IUPPTL)

      ENDIF

      CLOSE(ICHAN)
        
      ELSEIF(IBTYPE > 10 .AND. IBTYPE <= 20) THEN 

C ... Calculate the boundary values from free-stream condition
C ... Two-fluid velocities are put to the same values (has to be changed)

      U10 = FRSVEL*COS(ALPHA)*COS(BETA)
      V10 = FRSVEL*SIN(ALPHA)
      W10 = FRSVEL*COS(ALPHA)*SIN(BETA)

      IF(INITCC /= 2) THEN
        UU10  = U10
        VV10  = V10
        WW10  = W10
      ELSE IF(INITCC == 2) THEN ! A pull-up case
        UU10  = 0. 
        VV10  = 0.
        WW10  = 0.
      ENDIF

      IF(M == 1) APATCH(IPL) = 0.

       DO J = KY1,KY2
          IJ = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
          DO I = KX1,KX2
             L  = 1 + (IN+I-1)*ISTR + IJ     ! CELL INDEX UP (STRIDE=ISTR)
             NP = I + NR                     ! Patch index with LN ghost cells
             XN1 = AAX(L,1+(IAP-1)*3)
             XN2 = AAX(L,2+(IAP-1)*3)
             XN3 = AAX(L,3+(IAP-1)*3)
C??             IF(M == 1) APATCH(IPG)= APATCH(IPL) + AAA(L,IAP)
             IF(M == 1) APATCH(IPL)= APATCH(IPL) + AAA(L,IAP)
             BOUNR(NP)  = FRSDEN
             BOUNU(NP)  = U10
             BOUNV(NP)  = V10
             BOUNW(NP)  = W10
             BOUNMF(NP) = (XN1*BOUNU(NP) + XN2*BOUNV(NP) + XN3*BOUNW(NP)
     &                    -UUROT(L,IAP))*BOUNR(NP)
             BOUNE(NP)  = FRSDEN*(FRSSIE +
     &       .5*(UU10**2 + VV10**2 + WW10**2))
             BOUNT(NP)  = FRSTEM
             BOUNP(NP)  = FRSPRE
     &              + FRSDEN*(GX*(XC(NP) - GROUND)
     &              +         GY*(YC(NP) - GROUND)
     &              +         GZ*(ZC(NP) - GROUND))
             BOUNPD(NP) = 0.

C ... Void fractions for the cavitation and two-fluid models
             BOUNA2(NP) = FRSALFA(2)
             BOUNA1(NP) = 1.- FRSALFA(2)

c             BOUNPD(NP) = BOUNP(NP) - FRSPRE ! Haasteellinen idea

               IF(.NOT. MULPHL) THEN
                UK = BOUNU(NP)
                VK = BOUNV(NP)
                WK = BOUNW(NP)
               CALL ROFPT(BOUNT(NP),BOUNP(NP),BOUNR(NP),1,ISTATE,RGAS,
     +          GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
                CALL EFPT(E1(1),BOUNP(NP),BOUNT(NP),1,ISTATE,RGAS,GAMMA,
     +          FRSPRE,FRSSIE,E0REF,T0REF)
                IF(ITURB >= 3 .AND. ITURB /= 8) 
     +          E1(1) = E1(1) + BOUNRK(NP)/BOUNR(NP)
                BOUNE(NP)  = BOUNR(NP)*(E1(1)+.5*(UK**2+VK**2+WK**2))
               ELSE IF(MULPHL) THEN
                DO IPHASE = 1,NPHASE
                ISTATE    = JSTATE(NGL,IPHASE)
               CALL ROFPT(BOUNT(NP),BOUNP(NP),RO1(IPHASE),1,ISTATE,RGAS,
     +          GAMMA,FRADEN(IPHASE),FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
                CALL EFPT(E1(IPHASE),BOUNP(NP),BOUNT(NP),1,ISTATE,RGAS,
     +          GAMMA,FRSPRE,FRSSIE,E0REF,T0REF)
                ENDDO
                BOUNR(NP) = BOUNA1(NP)*RO1(1) + BOUNA2(NP)*RO1(2)
                BOUNE(NP) = BOUNA1(NP)*RO1(1)*E1(1) + BOUNA2(NP)*
     +          RO1(2)*E1(2) + .5*BOUNR(NP)*(UK**2+VK**2+WK**2)
                IF(ITURB >= 3 .AND. ITURB /= 8) 
     +          BOUNE(NP) = BOUNE(NP) + BOUNRK(NP)
               ENDIF ! .NOT. MULPHL

C ... Two-fluid velocities
             BOUNU1(NP) = BOUNU(NP);  BOUNU2(NP) = BOUNU(NP)
             BOUNV1(NP) = BOUNV(NP);  BOUNV2(NP) = BOUNV(NP)
             BOUNW1(NP) = BOUNW(NP);  BOUNW2(NP) = BOUNW(NP)

             IF(ITURB >= 3 .AND. ITURB /= 8) THEN ! Turbulence variables
             BOUNRK(NP) = FRSRK
             BOUNEP(NP) = FRSEPS
             ENDIF
             IF(TRANSL) THEN ! Intermittency variables
             CALL INLETRET(FRSTUR,RETIN)
             BOUNG(NP)   = 1.
             BOUNRET(NP) = RETIN
             ENDIF

             IF(ISTRES > 0) THEN
                WRITE(*,*) 'Cannot establish b_ij values'
                WRITE(*,*) 'Use boundary files to give ditributions'
                WRITE(*,*) 'Exiting in BOUNDP...'
                STOP
             ENDIF
          ENDDO
       ENDDO

      ELSE IF(IBTYPE <= 10) THEN

      INQUIRE(FILE=BOUNDF(IPL),EXIST=THERE)
      IF(.NOT. THERE) THEN
         WRITE(*,*) 'Cannot find file:'
         WRITE(*,'(A50)') BOUNDF(IPL)
         WRITE(*,*) NBL !,IBWALL
         WRITE(*,*) ' Exiting in BOUNDP...'
         STOP
      ENDIF

C ... Read the boundary values from a file (Two-fluid still missing)

      OPEN(51,FILE=BOUNDF(IPL),STATUS='OLD',FORM='FORMATTED')

      CALL READIO(BOUNMF(IF1),BOUNR(IF1),BOUNU(IF1),BOUNV(IF1),
     +   BOUNW(IF1),BOUNE(IF1),BOUNT(IF1),BOUNP(IF1),BOUNRK(IF1),
     +   BOUNEP(IF1),BOUNFI(IF1,1),BOUNBI(IF1,1),ITURB,BOUNDF(IPL),
     +   BOUNA1(IF1),BOUNA2(IF1),BOUNG(IF1),BOUNRET(IF1),
     +   IMAX0,JMAX0,IMAX1,JMAX1,IFILE,NSCAL,
     +   MAXEB,MAXSB,ISTRES,IN,JN,KN,MULPHL,TRANSL,
     +   ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSVEL,
     +   FRSTEM,IUPPTL)      
      CLOSE(51)

       DO J = KY1,KY2
          IJ = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
          DO I = KX1,KX2
             L  = 1 + (IN+I-1)*ISTR + IJ     ! CELL INDEX UP (STRIDE=ISTR)
             NP = I + NR                     ! Patch index with LN ghost cells
             IF(IFILE == 0) THEN             ! Primitive variables
                UK = BOUNU(NP)
                VK = BOUNV(NP)
                WK = BOUNW(NP)
                BOUNPD(NP) = BOUNP(NP) - FRSPRE
     &              - FRSDEN*(GX*(XC(NP) - GROUND)
     &              +         GY*(YC(NP) - GROUND)
     &              +         GZ*(ZC(NP) - GROUND))
               IF(.NOT. MULPHL) THEN
                CALL ROFPT(BOUNT(NP),BOUNP(NP),BOUNR(NP),1,ISTATE,RGAS,
     +          GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
                CALL EFPT(E1(1),BOUNP(NP),BOUNT(NP),1,ISTATE,RGAS,GAMMA,
     +          FRSPRE,FRSSIE,E0REF,T0REF)
                IF(ITURB >= 3 .AND. ITURB /= 8) 
     +          E1(1) = E1(1) + BOUNRK(NP)/BOUNR(NP)
                BOUNE(NP)  = BOUNR(NP)*(E1(1)+.5*(UK**2+VK**2+WK**2))
               ELSE IF(MULPHL) THEN
                DO IPHASE = 1,NPHASE
                ISTATE    = JSTATE(NGL,IPHASE)
               CALL ROFPT(BOUNT(NP),BOUNP(NP),RO1(IPHASE),1,ISTATE,RGAS,
     +          GAMMA,FRADEN(IPHASE),FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
                CALL EFPT(E1(IPHASE),BOUNP(NP),BOUNT(NP),1,ISTATE,RGAS,
     +          GAMMA,FRSPRE,FRSSIE,E0REF,T0REF)
                ENDDO
                BOUNR(NP) = BOUNA1(NP)*RO1(1) + BOUNA2(NP)*RO1(2)
                BOUNE(NP) = BOUNA1(NP)*RO1(1)*E1(1) + BOUNA2(NP)*
     +          RO1(2)*E1(2) + .5*BOUNR(NP)*(UK**2+VK**2+WK**2)
                IF(ITURB >= 3 .AND. ITURB /= 8) 
     +          BOUNE(NP) = BOUNE(NP) + BOUNRK(NP)
               ENDIF ! .NOT. MULPHL
             ELSE ! Conservative variables were given (ideal gas law)
                IF(ISTATE /= 1) THEN
                   WRITE(*,*) ' ISTATE /= 1 requires currently',
     +             ' IFILE = 0. Exiting in BOUNDP...'
                   STOP
                ENDIF
                UK = BOUNU(NP)/BOUNR(NP)
                VK = BOUNV(NP)/BOUNR(NP)
                WK = BOUNW(NP)/BOUNR(NP)
                BOUNU(NP) = UK
                BOUNV(NP) = VK
                BOUNW(NP) = WK
                E1(1)     = BOUNE(NP)/BOUNR(NP) -.5*(UK**2+VK**2+WK**2)
                IF(ITURB >= 3 .AND. ITURB /= 8)
     +          E1(1) = E1(1) - BOUNRK(NP)/BOUNR(NP)
                BOUNPD(NP) = (GAMMA-1.)*BOUNR(NP)*E1(1) - FRSPRE
                BOUNP(NP)  = (GAMMA-1.)*BOUNR(NP)*E1(1)
                CALL TFROE(BOUNR(NP),E1(1),BOUNT(NP),1,ISTATE,RGAS,
     +          GAMMA,FRSPRE,E0REF,T0REF)
                IF(ISTATE >= 6) THEN ! viritys maximus ei pellaa enää..
                   BOUNPD(NP) = 0.
                   BOUNP(NP)  = BOUNPD(NP) + FRSPRE
                   CALL TFROE(BOUNR(NP),E1(1),BOUNT(NP),1,ISTATE,
     +             RGAS,GAMMA,FRSPRE,E0REF,T0REF)                
                ENDIF
             ENDIF ! IFILE == 0
             XN1        = AAX(L,1+(IAP-1)*3)
             XN2        = AAX(L,2+(IAP-1)*3)
             XN3        = AAX(L,3+(IAP-1)*3)
             BOUNMF(NP) = (XN1*BOUNU(NP) + XN2*BOUNV(NP) + XN3*BOUNW(NP)
     +                    -UUROT(L,IAP))*BOUNR(NP)
          ENDDO
       ENDDO

      ELSE

         WRITE(*,*) 'Strange boundary type in BOUNDP'
         WRITE(*,*) IBTYPE
         WRITE(*,*) ' Exiting in BOUNDP...'

      ENDIF  ! IBTYPE > 10

      ENDIF  ! IBC == 3 .OR. IBC == 5

 8000 CONTINUE

      DEALLOCATE(AAX,AAA,UUROT,APU)

      RETURN
      END SUBROUTINE BOUNDP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE READHT(HFLUX,IMAX,JMAX,IMAX1,JMAX1)

      USE NS3CO, ONLY : IN, JN, KN
      
      IMPLICIT NONE

      INTEGER :: IMAX1,JMAX1,IMAX0,JMAx0,NDIVI,NDIVJ,NTWO,NONE1I,
     +     NONE1J,N,IG,I2,J2,IA,I,J,IMAX,JMAX
      
      REAL :: HFLUX(*), APU1

      CHARACTER(LEN=3) :: FILEE

C ... Read the boundary surface distribution from a file (not tested)

      READ(51,*) IMAX0,JMAX0
        
      IF ((MOD(IMAX0,IMAX1) /= 0) .OR. (MOD(JMAX0,JMAX1) /= 0)) THEN
         WRITE(*,*) 'In inlet file:',IMAX0,JMAX0
         WRITE(*,*) 'In INPUT file:',IMAX ,JMAX
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',
     +        'THERMO'//FILEE
         WRITE(*,*) 'EXITING...'
         STOP
      ENDIF

      IF (IMAX0/IMAX1 /= JMAX0/JMAX1 .AND. IMAX0 /= 1 .AND.
     +      JMAX0 /= 1) THEN
         WRITE(*,*) 'In inlet file:',imax0,jmax0
         WRITE(*,*) 'In INPUT file:',imax ,jmax
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',
     +        'THERMO'//FILEE
         WRITE(*,*) 'EXITING...'
         STOP
      ENDIF
 
      DO N = 1,IMAX1*JMAX1  ! initialized
         HFLUX(N)  = 0.
      ENDDO
      
      NDIVI  = INT(IMAX0/IMAX1)
      NDIVJ  = INT(JMAX0/JMAX1)
      NTWO   = NDIVI*NDIVJ
      NONE1I = NDIVI-1
      NONE1J = NDIVJ-1
 
C ... Read the heat flux in

      DO J = 1,JMAX1
c         IG = (J-1)*IMAX1 ! Old without the ghhost cells
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  READ(51,*) APU1
                  HFLUX(IA)  = HFLUX(IA) + APU1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      DO N  = 1,IMAX1*JMAX1
         HFLUX(N)  = HFLUX(N)/NTWO
      ENDDO

      RETURN
      END SUBROUTINE READHT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE READIO(BOUNMF,BOUNR,BOUNU,BOUNV,BOUNW,BOUNE,BOUNT,
     +     BOUNP,BOUNRK,BOUNEP,BOUNFI,BOUNBI,ITURB,BOUNDF,BOUNA1,
     +     BOUNA2,BOUNG,BOUNRET,IMAX,JMAX,IMAX1,JMAX1,IFILE,NSCAL,MAXEB,
     +     MAXSB,ISTRES,IN,JN,KN,MULPHL,TRANSL,
     +     ISTATE,RGAS,GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSVEL,FRSTEM,
     +     IUPPTL)  

      USE CONSTANTS, ONLY : EPS
     
      IMPLICIT NONE

      INTEGER :: IMAX1,JMAX1,IMAX0,JMAX0,NDIVI,NDIVJ,NTWO,NONE1I,
     +     NONE1J,N,IG,I2,J2,IA,I,J,IMAX,JMAX,IFILE,IFILE1,NSCAL1,
     +     NSCAL,MAXEB,MAXSB,ITURB,NS,ISTRES,IN,JN,KN,ISTATE,IUPPTL
      
      REAL :: BOUNMF(*),BOUNR(*),BOUNU(*),BOUNV(*),BOUNW(*),BOUNE(*),
     +     BOUNT(*),BOUNP(*),BOUNRK(*),BOUNEP(*),BOUNA1(*),
     +     BOUNA2(*),BOUNFI(MAXSB,*),BOUNBI(MAXEB,*),BOUNG(*),BOUNRET(*)

      REAL :: APU1, APU2, APU3, APU4, APU5, ROIN, RGAS, GAMMA,
     +     FRSDEN, FRSPRE, E0REF, T0REF, TU, FRSVEL, FRSTEM

      LOGICAL :: MULPHL, TRANSL

      CHARACTER(LEN=3)  :: FILEE
      CHARACTER(LEN=80) :: BOUNDF, LINE1

C ... Read the inlet/outlet boundary surface distribution from a file

      IFILE1 = 1 ! Old default value

      READ(51,'(1A80)') LINE1

         READ(LINE1,*,END=2013,ERR=2011) NSCAL1,IFILE 
         GO TO 2012
2011     WRITE(*,*) 'The following card is broken or erronous:'
         WRITE(*,*) LINE1
         STOP
2013     IFILE = IFILE1 ! With old files, the old default value
2012  CONTINUE     ! This line succesfully read


      IF (NSCAL1 < NSCAL) THEN
         WRITE(*,*) 'WRONG NUMBER OF SCALARS AT BOUNDARIES AT FILE',
     +        BOUNDF
         WRITE(*,*) 'EXITING...'
*         CALL EXIT
         STOP
      ENDIF

      READ(51,*) IMAX0,JMAX0
        
      IF ((MOD(IMAX0,IMAX1) /= 0) .OR. (MOD(JMAX0,JMAX1) /= 0)) THEN
         WRITE(*,*) 'In inlet file:',IMAX0,JMAX0
         WRITE(*,*) 'In INPUT file:',IMAX ,JMAX
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',
     +        'THERMO'//FILEE
         WRITE(*,*) 'EXITING...'
         STOP
      ENDIF

      IF (IMAX0/IMAX1 /= JMAX0/JMAX1 .AND. IMAX0 /= 1 .AND.
     +      JMAX0 /= 1) THEN
         WRITE(*,*) 'In inlet file:',imax0,jmax0
         WRITE(*,*) 'In INPUT file:',imax ,jmax
         WRITE(*,*) 'DIMENSIONS AND LEVEL DO NOT FIT IN FILE ',
     +        'THERMO'//FILEE
         WRITE(*,*) 'EXITING...'
         STOP
      ENDIF

C ... Initialization of the main variable

      DO 100 N = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
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
 100  CONTINUE

C ... Initialization of the turbulence variables

      IF (ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 200 N = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNRK(N) = 0.
         BOUNEP(N) = 0.
 200  CONTINUE
      ENDIF

C ... Initialization of the Intermittency variables

      IF (TRANSL) THEN
      DO 250 N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNG(N)  = 0.
         BOUNRET(N)= 0.
 250  CONTINUE
      ENDIF

C ... Initialization of scalars

      DO 300 NS = 1,NSCAL
      DO 300 N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNFI(N,NS) = 0.
 300  CONTINUE

C ... Initialization of EARSM anisotropic components

      IF(ISTRES > 0) THEN
      DO 350 NS = 1,6
      DO 350 N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNBI(N,NS) = 0.
 350  CONTINUE
      ENDIF

      NDIVI = INT(IMAX0/IMAX1)
      NDIVJ = INT(JMAX0/JMAX1)
      NTWO    = NDIVI*NDIVJ
      NONE1I  = NDIVI-1
      NONE1J  = NDIVJ-1
      
C ... Read primary variables in

      DO J = 1,JMAX1
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  READ(51,*) APU1,APU2,APU3,APU4,APU5
                  IF(IFILE == 0) THEN ! Primitive variables were given
                     BOUNP(IA)  = BOUNP(IA) + APU1
                     BOUNU(IA)  = BOUNU(IA) + APU2
                     BOUNV(IA)  = BOUNV(IA) + APU3
                     BOUNW(IA)  = BOUNW(IA) + APU4
                     BOUNT(IA)  = BOUNT(IA) + APU5
                     CALL ROFPT(APU5,APU1,ROIN,1,ISTATE,RGAS,
     +               GAMMA,FRSDEN,FRSPRE,E0REF,T0REF,FRSTEM,IUPPTL)
                     BOUNR(IA)  = BOUNR(IA) + ROIN
                  ELSE ! Conservative variables are given
                     BOUNR(IA)  = BOUNR(IA) + APU1
                     BOUNU(IA)  = BOUNU(IA) + APU2
                     BOUNV(IA)  = BOUNV(IA) + APU3
                     BOUNW(IA)  = BOUNW(IA) + APU4
                     BOUNE(IA)  = BOUNE(IA) + APU5
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNR(N)  = BOUNR(N)/NTWO
         BOUNU(N)  = BOUNU(N)/NTWO
         BOUNV(N)  = BOUNV(N)/NTWO
         BOUNW(N)  = BOUNW(N)/NTWO
         BOUNE(N)  = BOUNE(N)/NTWO
         BOUNP(N)  = BOUNP(N)/NTWO
         BOUNT(N)  = BOUNT(N)/NTWO
      ENDDO

C ... Read turbulence varibales in

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO J = 1,JMAX1
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  IF(TRANSL) THEN ! Include intermittency variables
                  READ(51,*) APU1,APU2    !,APU3,APU4
                  BOUNRK(IA)  = BOUNRK(IA) + APU1
                  BOUNEP(IA)  = BOUNEP(IA) + APU2
                  BOUNG(IA)   = 1.    !BOUNG(IA)  + APU3 ! To NEWBOUND?
C                  BOUNRET(IA) = BOUNRET(IA)+ APU4
                  ELSE ! Only two-equation model
                  READ(51,*) APU1,APU2
                  BOUNRK(IA)  = BOUNRK(IA) + APU1
                  BOUNEP(IA)  = BOUNEP(IA) + APU2
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNRK(N)  = BOUNRK(N)/NTWO
         BOUNEP(N)  = BOUNEP(N)/NTWO
         IF(TRANSL) THEN
         BOUNG(N)   = BOUNG(N)/NTWO
         TU         = SQRT(2./3.*BOUNRK(N)/(BOUNR(N)+EPS))/FRSVEL
         CALL INLETRET(TU,APU4)
         BOUNRET(N) = APU4
         ENDIF
      ENDDO

      IF (ISTRES > 0) THEN ! b_ij distributions
      DO J = 1,JMAX1
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  READ(51,*) APU1,APU2,APU3,APU4,APU5
                  BOUNBI(IA,1) = BOUNBI(IA,1) + APU1
                  BOUNBI(IA,2) = BOUNBI(IA,2) + APU2
                  BOUNBI(IA,3) = BOUNBI(IA,3) + APU3
                  BOUNBI(IA,4) = BOUNBI(IA,4) + APU4
                  BOUNBI(IA,5) = BOUNBI(IA,5) + APU5
                  BOUNBI(IA,6) =-BOUNBI(IA,1) - BOUNBI(IA,4)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO NS = 1,6
      DO N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNBI(N,NS) = BOUNBI(N,NS)/NTWO
      ENDDO
      ENDDO
      ENDIF ! ISTRES > 0
      ENDIF ! ITURB >= 3

C ... read the multiphase varibales in

      IF(MULPHL) THEN
      DO J = 1,JMAX1
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  READ(51,*) APU1,APU2 ! Maybe something else
                  BOUNA1(IA)  = BOUNA1(IA) + APU1
                  BOUNA2(IA)  = BOUNA2(IA) + APU2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNA1(N)  = BOUNA1(N)/NTWO
         BOUNA2(N)  = BOUNA2(N)/NTWO
      ENDDO
      ENDIF ! MULPHL

C ... Read scalar variables in

      DO NS = 1,NSCAL
      DO J = 1,JMAX1
         IG = (J+JN-1)*(IMAX1+2*IN) + IN
         DO J2= 0,NONE1J
            DO I = 1,IMAX1
               IA = IG + I
               DO I2= 0,NONE1I
                  READ(51,*) APU1
                  BOUNFI(IA,NS) = BOUNFI(IA,NS) + APU1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      DO NS = 1,NSCAL
      DO N  = 1,(IMAX1+2*IN)*(JMAX1+2*JN)
         BOUNFI(N,NS) = BOUNFI(N,NS)/NTWO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE READIO
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE READIOP(ICHAN,XMASS,HTOT,P,PTOT,TEMP,BIJ,
     &                   TURBLE,RMUINI,RMUL1,RMUL2,TG,TRET,
     &                   FISET,INPTYP,TRANSL,KSCAL,LRATUP)

C ... Subroutine reads INL/OUT flow parameters instead of all boundary
C ... values (see subroutine READIO).

      IMPLICIT NONE

      INTEGER :: ICHAN, INPTYP, I, KSCAL
      REAL :: APU1, APU2, APU3, TURBLE, RMUINI, RMUL1, RMUL2, LRATUP
      REAL :: BIJ(*), FISET(*)
      REAL :: XMASS, HTOT, P, PTOT, TEMP, TG, TRET

      LOGICAL :: TRANSL

      CHARACTER(80) :: LINE

      XMASS = 0.0; HTOT = 0.0; P = 0.0; PTOT = 0.0; TEMP = 0.0
      

C ... Read basic variables 
      
 10   READ(UNIT=ICHAN,FMT='(A)') LINE
      IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 10
      BACKSPACE ICHAN
      READ(ICHAN,*) INPTYP, APU1, APU2, APU3

      SELECT CASE(INPTYP)
      CASE(1)
         XMASS = APU1
         HTOT  = APU2
         P     = APU3
      CASE(2)
         XMASS = APU1
         HTOT  = APU2
         PTOT  = APU3
      CASE(3)
         XMASS = APU1
         HTOT  = APU2
         TEMP  = APU3
      CASE(4)
         XMASS = APU1
         P     = APU2
         TEMP  = APU3
      END SELECT
      

C ... Read turbules variables 

      IF(TRANSL) THEN
 20      READ(UNIT=ICHAN,FMT='(A)') LINE
         IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 20
         BACKSPACE ICHAN
         READ(ICHAN,*) TURBLE, RMUINI, TG, TRET
      ELSE
 21      READ(UNIT=ICHAN,FMT='(A)') LINE
         IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 21
         BACKSPACE ICHAN
         READ(ICHAN,*) TURBLE, RMUINI
      ENDIF


C ... Read optional Reynolds stesses
      
 30   READ(UNIT=ICHAN,FMT='(A)') LINE
      IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 30
      BACKSPACE ICHAN
      READ(ICHAN,*) (BIJ(I),I=1,5)


C .. Read optional multiphase variables
      
 40   READ(UNIT=ICHAN,FMT='(A)') LINE
      IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 40
      BACKSPACE ICHAN
      READ(ICHAN,*) RMUL1, RMUL2


C ... Read a representative relative length of the tube upstream of
C ... the boundary surface to enable a reasonable flow distribution
C ... to be generated as an output option. L/A = LRATUP < 0.1 means
C ... bulk flow.      

 50   READ(UNIT=ICHAN,FMT='(A)') LINE
      IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 50
      BACKSPACE ICHAN
      READ(ICHAN,*) LRATUP


C ... Read optional scalars
      
      IF(KSCAL > 0) THEN
 60      READ(UNIT=ICHAN,FMT='(A)') LINE
         IF(LINE(1:1) == '#' .OR. LEN(TRIM(LINE)) == 0) GOTO 60
         BACKSPACE ICHAN
         DO I=1,KSCAL
            READ(ICHAN,*) APU1
            FISET(I) = APU1
         ENDDO
      ENDIF
        
      RETURN
      END SUBROUTINE READIOP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE READAD(LENGTH,IACTU,IACTY,XCP,YCP,ZCP,KA,KBIND,
     +      KX1,KY1,KX2,KY2,ISTR,JSTR,IHFIPL,SURFA,SURFT,IPRO)

      USE CONSTANTS, ONLY : PII

      USE NS3CO, ONLY : IN, JN, KN, LN

      USE FLIGHT, ONLY : THRUST,TORQUE,ACTUA,OSKU,CONEA,
     +   ROTA1,XCG,YCG,ZCG,ROTB1,RIN,ROUT

      IMPLICIT NONE

      INTEGER :: III,LENGTH,IACTU,IACTY,NINTERP,KA,KBIND,TQINT,
     +           KX1,KY1,KX2,KY2,IHFIPL,I,J,IJ,II,NR,NP,ISTR,JSTR,IPRO

      REAL :: RR1,RR2,DISTT1,DISTQ1,DISTT2,DISTQ2,DIST,PLOSS2,RR,
     +        RVECLEN,R2,XAC,YAC,ZAC,A0,A1,B1,XLOC,YLOC,ZLOC,RX,RY,RZ,
     +        CLDIST,THEANG1,CLKOE,PLOSS2M,R1

      REAL :: XFIN, YFIN, ZFIN

      REAL :: XCP(*),YCP(*),ZCP(*)

      REAL :: SURFA(*),SURFT(*)

      CHARACTER (LEN=180) :: LINE1, HEADER

      LOGICAL :: THERE

      REAL, ALLOCATABLE ::ACTX(:),ACTY(:),ACTZ(:),ACTFA(:),ACTFT(:),
     +       rR3(:),THEANG(:),DISTT3(:),DISTQ3(:)

C ... Actuatuator force distribution is read in

C ... Allocate temporary arrays for read in the distributions

      ALLOCATE(ACTX(LENGTH),ACTY(LENGTH),ACTZ(LENGTH),ACTFA(LENGTH),
     + ACTFT(LENGTH),THEANG(LENGTH),DISTT3(LENGTH),DISTQ3(LENGTH),
     + rR3(LENGTH))

      R2 = ROUT(IACTU)
      R1 = RIN(IACTU)
      A0 = CONEA(IACTU)
      A1 = ROTA1(IACTU)
      B1 = ROTB1(IACTU) 


      SELECT CASE(IACTY)

         CASE(0:3) ! No force, constant or Bramwell distribution

         CASE(6)   ! Miklos' diploma thesis

C****************************************************************************

C ... Check if "T_Q_DIST.INPUT" file exsists and open it to channel 520+NGL
          INQUIRE(FILE='T_Q_DIST.INPUT',EXIST=THERE)
          IF(THERE) THEN          
             OPEN(51, FILE='T_Q_DIST.INPUT',
     &            STATUS='UNKNOWN', FORM='FORMATTED')
          ELSE 
             WRITE(*,*)' '
             WRITE(*,*)' Warning: T_Q_DIST.INPUT',
     &            ' was not found'
             WRITE(*,*)' in ACTUATOR_FORCES Subroutine.' 
             WRITE(*,*)' '
             WRITE(*,*)' Exiting...'
             WRITE(*,*)' '
             STOP
          ENDIF
C ... Read data in
          READ(51,'(A5)')HEADER
          READ(51,*)NINTERP
          READ(51,'(A5)')HEADER
          READ(51,*)rR1,DISTT1,DISTQ1
          READ(51,*)rR2,DISTT2,DISTQ2

C ... Apply a closest radius search

      DO J = KY1,KY2
          IJ = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHFIPL
        DO I = KX1,KX2
C         Cell II index takes into account the top surface via KBIND
          II = 1 + (IN+I-1)*ISTR + IJ + KBIND
          NP = I + NR                        ! Patch index with LN ghost cells

C ... Set rR
          RX       = XCP(NP)-XCG(IACTU)
          RY       = YCP(NP)-YCG(IACTU)
          RZ       = ZCP(NP)-ZCG(IACTU)
          RVECLEN  = SQRT(RX**2+RY**2+RZ**2)
          rR=RVECLEN/R2

C ... Search correct location
          DO TQINT=1,NINTERP-2
             IF (rR == rR1)THEN
                DIST= DISTT1
             ELSEIF(rR == rR2)THEN
                DIST= DISTT2
             ELSEIF(rR > rR1.AND.rR < rR2)THEN
                DIST= DISTT1+(DISTT2-DISTT1)/(rR2-rR1)*(rR-rR1)
                GO TO 115
             ELSE
                rR1=Rr2
                DISTT1=DISTT2
                READ(51,*)rR2,DISTT2,DISTQ2 
             ENDIF
          ENDDO
C... Update PLOSS2 and PLOSS2M
 115   PLOSS2    = DIST*THRUST(IACTU)/ACTUA(IACTU)
       PLOSS2M   = DIST*1.5*TORQUE(IACTU)/(PII*(R2**3-R1**3))
       SURFA(NP) = PLOSS2
       SURFT(NP) = PLOSS2M
      ENDDO; ENDDO

      CLOSE (51)


      CASE(7,8)   ! Uneven distribution developed from (7)

C****************************************************************************
C ... Reads distribution data in
C ... Header lines

      READ(51,'(A5)')HEADER
      READ(51,'(A5)')HEADER
      READ(51,'(A5)')HEADER
      READ(51,'(A5)')HEADER

C ... Read data for thrust and torque

      DO III = 1,LENGTH
          READ(51,*,ERR= 118,END=119) rR1,THEANG1,DISTT1,DISTQ1
          rR3(III) = rR1
          THEANG(III) = THEANG1
          DISTT3(III) = DISTT1
          DISTQ3(III) = DISTQ1
      ENDDO

118   WRITE(*,*) ' End of the actuator-disc data detected in READAD'
      STOP
119   WRITE(4,"(2X,A,I4)") 'Acturator data succesfully read. Length =',
     +      LENGTH
      CLOSE (51)

C ... Apply a closest point search

      DO J = KY1,KY2
          IJ = (JN+J-1)*JSTR + KA
          NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHFIPL
        DO I = KX1,KX2
C         Cell II index takes into account the top surface via KBIND
          II = 1 + (IN+I-1)*ISTR + IJ + KBIND
          NP = I + NR                        ! Patch index with LN ghost cells
          CLDIST = 1.E6

        DO III = 1,LENGTH-4 ! Headers are substracted from the loop
          Xloc=0.0
          Yloc=SIN(THEANG(III)*(PII/180.))*rR3(III)*R2
          Zloc=COS(THEANG(III)*(PII/180.))*rR3(III)*R2
C ... Changes cooridate system
          CALL EULFIN(0.0,A1,B1,Yloc,-Zloc,Xloc,Xfin,Yfin,Zfin,1)
C ... Search the closest point and stores data (+normalization)
          Xac=Xfin+XCG(IACTU)
          Yac=Yfin+YCG(IACTU)
          Zac=Zfin+ZCG(IACTU)
          CLKOE = (Xac-XCP(NP))**2+(Yac-YCP(NP))**2+(Zac-ZCP(NP))**2 
          IF(CLKOE <= CLDIST)  THEN
             DISTT2=DISTT3(III)*OSKU(IACTU)%DAMPT
             DISTQ2=DISTQ3(III)*OSKU(IACTU)%DAMPN
             CLDIST=CLKOE
          ELSE 
             CONTINUE
          ENDIF

        ENDDO ! DO III

C... Update PLOSS2-arrays

       PLOSS2    = DISTT2*THRUST(IACTU)/ACTUA(IACTU)
       PLOSS2M   = DISTQ2*1.5*TORQUE(IACTU)/(PII*(R2**3-R1**3))
       SURFA(NP) = PLOSS2
       SURFT(NP) = PLOSS2M

      ENDDO; ENDDO

      CASE(9,10)   ! A general interpolation method

C****************************************************************************

C ... Read the distribution in a xyz-form

      READ(51,*,ERR=218,END=219) LINE1 ! Skip the heading
      DO III = 1,LENGTH-1
      READ(51,*,ERR=218,END=219) ACTX(III),ACTY(III),ACTZ(III),
     &   ACTFA(III),ACTFT(III)
      ENDDO
218   WRITE(*,*) ' End of the actuator-disc data detected in READAD'
      STOP
219   WRITE(4,"(2X,A,I4)") 'Acturator data succesfully read. Length =',
     +      LENGTH
      CLOSE(51) ! The same file is read for each grid level !!

C ... Apply a general interpolation method
      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHFIPL
         DO I = KX1,KX2
C ...       Cell II index takes into account the top surface via KBIND
            II = 1 + (IN+I-1)*ISTR + IJ + KBIND
            NP = I + NR                      ! Patch index with LN ghost cells
            CLDIST = 1.E6

C ... Update PLOSS2-arrays

            PLOSS2    = 0.  ! Saab
            PLOSS2M   = 0.  ! Interpolation to these using ACT*-arrays
            SURFA(NP) = PLOSS2
            SURFT(NP) = PLOSS2M

         ENDDO
      ENDDO

      DEALLOCATE(ACTX,ACTY,ACTZ,ACTFA,ACTFT,THEANG,DISTT3,DISTQ3,rR3)

      END SELECT

      RETURN
      END SUBROUTINE READAD
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE EXTPAT(QWALL,KX1,KX2,KY1,KY2)

      USE NS3CO, ONLY : LN
      
      IMPLICIT NONE

      INTEGER :: JSTRID, KX1, KX2, KY1, KY2, I, J, NR, NP

      REAL :: QWALL(*)

C ... Solid patches are extended into a first row of ghost cells
C ... Start with I = KX1 surface

      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
         I  = KX1
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2) THEN
            QWALL(NP-1) = 2.*QWALL(NP) - QWALL(NP+1)
         ELSE
            QWALL(NP-1) = QWALL(NP)
         ENDIF
      ENDDO

C ... Continue with I = KX2 surface

      DO J = KY1,KY2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
         I  = KX2
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2) THEN
            QWALL(NP+1) = 2.*QWALL(NP) - QWALL(NP-1)
         ELSE
            QWALL(NP+1) = QWALL(NP)
         ENDIF
      ENDDO

C ... Then continue with J = KY1 surface

      J      = KY1
      JSTRID = KX2-KX1+1+2*LN
      NR     = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      DO I  = KX1,KX2
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KY1 /= KY2) THEN
            QWALL(NP-JSTRID) = 2.*QWALL(NP) - QWALL(NP+JSTRID)
         ELSE
            QWALL(NP-JSTRID) = QWALL(NP)
         ENDIF
      ENDDO

C ... Then J = KY2 surface

      J      = KY2
      JSTRID = KX2-KX1+1+2*LN
      NR     = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      DO I = KX1,KX2
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KY1 /= KY2) THEN
            QWALL(NP+JSTRID) = 2.*QWALL(NP) - QWALL(NP-JSTRID)
         ELSE
            QWALL(NP+JSTRID) = QWALL(NP)
         ENDIF
      ENDDO

C ... Then (I,J) = KX1-1,KY1-1

      J = KY1-1
      JSTRID = KX2-KX1+1+2*LN
      NR = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      I  = KX1-1
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2 .AND. KY1 /= KY2) THEN
            QWALL(NP) = 2.*QWALL(NP+JSTRID+1) - QWALL(NP+2*JSTRID+2)
         ELSE
            QWALL(NP) = QWALL(NP+JSTRID+1)
         ENDIF 

C ... Then (I,J) = KX2+1,KY1-1

      J = KY1-1
      JSTRID = KX2-KX1+1+2*LN
      NR = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      I  = KX2+1
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2 .AND. KY1 /= KY2) THEN
            QWALL(NP) = 2.*QWALL(NP+JSTRID-1) - QWALL(NP+2*JSTRID-2)
         ELSE
            QWALL(NP) = QWALL(NP+JSTRID-1)
         ENDIF 

C ... Then (I,J) = KX2+1,KY2+1

      J = KY2+1
      JSTRID = KX2-KX1+1+2*LN
      NR = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      I  = KX2+1
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2 .AND. KY1 /= KY2) THEN
            QWALL(NP) = 2.*QWALL(NP-JSTRID-1) - QWALL(NP-2*JSTRID-2)
         ELSE
            QWALL(NP) = QWALL(NP-JSTRID-1)
         ENDIF 

C ... Then (I,J) = KX1-1,KY2+1

      J = KY2+1
      JSTRID = KX2-KX1+1+2*LN
      NR = (LN+J-KY1)*JSTRID - KX1 + LN + 1
      I  = KX1-1
         NP = I + NR                         ! Patch index with LN ghost cells
         IF(KX1 /= KX2 .AND. KY1 /= KY2) THEN
            QWALL(NP) = 2.*QWALL(NP-JSTRID+1) - QWALL(NP-2*JSTRID+2)
         ELSE
            QWALL(NP) = QWALL(NP-JSTRID+1)
         ENDIF 

      RETURN
      END SUBROUTINE EXTPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE REDPAT(QWALL,ZZZ,KX1,KX2,KY1,KY2)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: KX1,KX2,KY1,KY2,I,J,JPSTR,NP,NR,NH,NN

      REAL :: QWALL(*), ZZZ(*)

C ... Solid patches are changed into a cell-vertex format

      DO J = KY1-2,KY2+2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
         DO I = KX1-2,KX2+2
            NP = I + NR                      ! Patch index with LN ghost cells
            ZZZ(NP) = QWALL(NP)
         ENDDO
      ENDDO

      DO J = KY1,KY2+1                       ! +1 Feb 17, 2015 ESa, see VERPAT

         NN    = (J-KY1)*(KX2-KX1+1+1) - KX1 + 1
         JPSTR = KX2-KX1+1+2*LN
         NR    = (LN+J-KY1)*JPSTR - KX1 + LN + 1

         DO I = KX1,KX2+1                    ! +1 Feb 17, 2015 ESa, see VERPAT

            NH = I + NN                      ! New patch index
            NP = I + NR                      ! Patch index with LN ghost cells
            QWALL(NH)= .25*(ZZZ(NP)      + ZZZ(NP-1)+
     &           ZZZ(NP-JPSTR)+ ZZZ(NP-JPSTR-1))

         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE REDPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE VERPAT(QWALL,ZZZ,KX1,KX2,KY1,KY2)

      USE NS3CO, ONLY : LN

      IMPLICIT NONE

      INTEGER :: KX1, KX2, KY1, KY2, I, J, JPSTR, NP, NR

      REAL :: QWALL(*), ZZZ(*)

C ... Solid patches are changed into a cell-vertex format,
C ... Output QWALL
      
      DO J = KY1,KY2+1
         JPSTR = KX2-KX1+1+2*LN
         NR = (LN+J-KY1)*JPSTR - KX1 + LN + 1
         DO I = KX1,KX2+1
            NP = I + NR                      ! Patch index with LN ghost cells
            QWALL(NP)= .25*(ZZZ(NP) + ZZZ(NP-1)+
     &                      ZZZ(NP-JPSTR) + ZZZ(NP-JPSTR-1))
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE VERPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE YPLUS1(YPLUS,OHMI,D,RO,VIS,
     + ISTR,JSTR,KSTR,IN,JN,KN,KX1,KX2,KY1,KY2,KZ1,KZ2,
     + IP,IDIR,UTAU)

      USE CONSTANTS, ONLY : EPS      
      USE NS3CO, ONLY : LN, WALLFUNL
      
      IMPLICIT NONE

      INTEGER :: ISTR,JSTR,KSTR,IN,JN,KN,I,J,IJ,KA,
     +           KX1,KX2,KY1,KY2,KZ1,KZ2,NP,NR,L,IP,IDIR

      REAL :: OHMI(*),YPLUS(*),D(*),RO(*),VIS(*),UTAU(*),
     +        DIST,OHMIL

C ... Calculates the y+ values in the first row of cells

      IF(KZ1 /= KZ2) THEN
         WRITE(*,*) ' KZ1 /= KZ2 in YPLUS1. Odd, returning..'
         RETURN
      ENDIF

      KA = (KN+KZ1-1)*KSTR

C ... Calculate the y+ values in the first cell for all solid patches

      DO J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
         DO I = KX1,KX2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! =II,CELL INDEX(STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            DIST  = .5*D(L)                  ! Accurate enough, notice .5
            OHMIL     = MAX(.5*(OHMI(L)+OHMI(L-KSTR)),0.)
            IF(.NOT.WALLFUNL) THEN
            YPLUS(NP) = DIST*SQRT(RO(L)*OHMIL/(VIS(L)+EPS))
            ELSEIF(WALLFUNL) THEN
            YPLUS(NP) = DIST*UTAU(NP)*RO(L)/(VIS(L)+EPS)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE YPLUS1
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE WRIPAT(QWALL,KX1,KX2,KY1,KY2,M,IN,JN,KN,frspre)

      USE NS3CO, ONLY : LN
      
      IMPLICIT NONE

      INTEGER :: KX1, KX2, KY1, KY2, M, IN, JN, KN, I, J, NR

      REAL :: QWALL(*), frspre

C ... Solid patches are printed into an ascii file for checking

      DO J = KY1-2,KY2+2
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + 1
                                             ! Patch index with LN ghost cells
         WRITE(98,*) (QWALL(I)/frspre,I,I=kx1-2+nr,kx2+2+nr)
      ENDDO

      RETURN
      END SUBROUTINE WRIPAT
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE PATVOL(WAVEH,DPU1,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     &     NBL,M,NPATCH,ICON,IHF,MAXB,IREPEA,NREPEA,SURCHA,BACK)

      USE NS3CO, ONLY : LN, IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,IL,NBL,M,NPATCH,MAXB,
     &           IN,JN,KN,INRE,JNRE,KNRE,
     &           IP,IBC,IPL,IPG,IDIR,IFACE,KX1,KX2,KY1,KY2,ISTR,JSTR,
     &           KSTR,I,J,K,L,KA,IJ,NR,NP
      INTEGER :: ICON(IC9,*), IHF(*), IREPEA(*), NREPEA(*)

      REAL :: WAVEH(*), DPU1(*)
      LOGICAL :: SOLID, INLOUT, SURFACE, BACK
      CHARACTER(LEN=3) :: SURCHA

C ... Put the patch values into a 3D block under the surfaces for connection

      IN = INRE
      JN = JNRE
      KN = KNRE

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID      
         
C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH

C ... A patch type and indicator for solids

         IBC = ICON(1,IP)

         SOLID   = IBC >=  8 .AND. IBC <= 10 .OR. IBC == 15
         INLOUT  = IBC >=  3 .AND. IBC <=  5  ! I/O or symmetry
         SURFACE = IBC == 13

         IF(SURFACE .AND. SURCHA == 'FRE' .OR.
     &      SOLID   .AND. SURCHA == 'SOL' .OR.
     &      INLOUT  .AND. SURCHA == 'MIR' .OR.
     &      INLOUT  .AND. SURCHA == 'INL' .OR.
     &      INLOUT  .AND. SURCHA == 'OUT') THEN ! Data is stored

         IPL   = ICON(2,IP) ! Proces local patch number
         IPG   = ICON(25,IP)! Global patch number
         IDIR  = 1
         IFACE = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1
         
C ... Patch indeces

         KX1 = ICON(4,IP)
         KX2 = ICON(5,IP)
         KY1 = ICON(6,IP)
         KY2 = ICON(7,IP)
             
         IF(IFACE == 2 .OR. IFACE == 5) THEN
           KY1 = ICON(4,IP)
           KY2 = ICON(5,IP)
           KX1 = ICON(6,IP)
           KX2 = ICON(7,IP)
         ENDIF

C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K = IMAX 
C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K = JMAX 
C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K = KMAX 
      ENDIF

      KA = (KN+K-1)*KSTR

C ... Loop over this patch

      IF(BACK) THEN         

         DO J = KY1-2,KY2+2
            IJ = (JN+J-1)*JSTR + KA
            NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
            DO I = KX1-2,KX2+2
               L  = 1 + (IN+I-1)*ISTR + IJ   ! CELL INDEX (STRIDE=ISTR)
               NP = I + NR                   ! Patch index with LN ghost cells
               WAVEH(NP) = DPU1(L)
            ENDDO 
         ENDDO

      ELSE IF(.NOT.BACK) THEN

      DO J = KY1-1,KY2+1
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1-1,KX2+1
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            DPU1(L) = WAVEH(NP)
         ENDDO 
      ENDDO

      ENDIF

      ENDIF ! SURFACE .AND. SURCHA == 'FRE'

 8000 CONTINUE

      RETURN
      END SUBROUTINE PATVOL
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE PATVOL_SP(WAVEH,DPU1,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     &     NBL,M,NPATCH,ICON,IHF,MAXB,IREPEA,NREPEA,SURCHA,BACK)

      USE NS3CO, ONLY : LN, IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,IL,NBL,M,NPATCH,MAXB,
     &           IN,JN,KN,INRE,JNRE,KNRE,
     &           IP,IBC,IPL,IPG,IDIR,IFACE,KX1,KX2,KY1,KY2,ISTR,JSTR,
     &           KSTR,I,J,K,L,KA,IJ,NR,NP
      INTEGER :: ICON(IC9,*), IHF(*), IREPEA(*), NREPEA(*)

      REAL :: WAVEH(*), DPU1(*)
      LOGICAL :: SOLID, INLOUT, SURFACE, BACK
      CHARACTER(LEN=3) :: SURCHA

C ... Put the patch values into a 3D block under the surfaces for connection

      IN = INRE
      JN = JNRE
      KN = KNRE

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      IL     = ISTRID*JSTRID      
         
C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH

C ... A patch type and indicator for solids

         IBC     = ICON(1,IP)

         SOLID   = IBC >=  8 .AND. IBC <= 10 .OR. IBC ==15
         INLOUT  = IBC >=  3 .AND. IBC <=  5  ! I/O or symmetry
         SURFACE = IBC == 13
           
         IF(SURFACE .AND. SURCHA == 'FRE' .OR.
     &      SOLID   .AND. SURCHA == 'SOL' .OR.
     &      INLOUT  .AND. SURCHA == 'MIR' .OR.
     &      INLOUT  .AND. SURCHA == 'INL' .OR.
     &      INLOUT  .AND. SURCHA == 'OUT') THEN ! Data is stored

         IPL    = ICON(2,IP) ! Proces local patch number
         IPG    = ICON(25,IP)! Global patch number
         IDIR   = 1
         IFACE  = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1
         
C ... Patch indeces

         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
             
         IF(IFACE == 2 .OR. IFACE == 5) THEN
           KY1  = ICON(4,IP)
           KY2  = ICON(5,IP)
           KX1  = ICON(6,IP)
           KX2  = ICON(7,IP)
         ENDIF

C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K = IMAX 

C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K = JMAX 

C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K = KMAX 
      ENDIF

      KA = (KN+K-1)*KSTR

C ... Loop over this patch

      IF(BACK) THEN       
      DO J = KY1-2,KY2+2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1-2,KX2+2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            WAVEH(NP) = DPU1(L)
         ENDDO
      ENDDO

      ELSE IF(.NOT.BACK) THEN ! Should take the actual patch size??

      DO J = KY1-1,KY2+1
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1-1,KX2+1
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            DPU1(L) = WAVEH(NP)
         ENDDO
      ENDDO

      ENDIF

      ENDIF ! SURFACE .AND. SURCHA == 'FRE'

 8000 CONTINUE

      RETURN
      END SUBROUTINE PATVOL_SP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
      SUBROUTINE PATVSP(WAVEH,APU1,APU2,IMAX,JMAX,KMAX,INRE,JNRE,KNRE,
     &     NBL,M,NPATCH,ICON,IHF,SURCHA,BACK)

      USE NS3CO, ONLY : LN, IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,ISTRID,JSTRID,IL,NBL,M,NPATCH,
     &           IN,JN,KN,INRE,JNRE,KNRE,
     &           IP,IBC,IPL,IPG,IDIR,IFACE,KX1,KX2,KY1,KY2,ISTR,JSTR,
     &           KSTR,I,J,K,L,KA,IJ,NR,NP,KALA,KYLA
      INTEGER :: ICON(IC9,*), IHF(*)

      REAL :: WAVEH(*)

      REAL :: APU1(*), APU2(*)

      LOGICAL :: SOLID, INLOUT, SURFACE, BACK
      CHARACTER(LEN=3) :: SURCHA

C ... Put the patch values into a 3D block under the surfaces for connection

      IN = INRE
      JN = JNRE
      KN = KNRE

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID      
         
C ... Loop over patches of this block

      DO 8000 IP = 1,NPATCH
C ... A patch type and indicator for solids
         IBC     = ICON(1,IP)

         SOLID   = IBC >= 8 .AND. IBC <= 10 .OR. IBC == 15
         INLOUT  = IBC >= 3 .AND. IBC <= 5  ! I/O
         SURFACE = IBC == 13
           
         IF(SURFACE .AND. SURCHA == 'FRE' .OR.
     &      SOLID   .AND. SURCHA == 'SOL' .OR.
     &      INLOUT  .AND. SURCHA == 'MIR' .OR.
     &      INLOUT  .AND. SURCHA == 'INL' .OR.
     &      INLOUT  .AND. SURCHA == 'OUT') THEN ! Data is stored

         IPL   = ICON(2,IP)   ! Proces local patch number
         IPG   = ICON(25,IP)  ! Global patch number
         IDIR  = 1
         IFACE = ICON(3,IP)
         IF(IFACE >= 4) IDIR = -1
         
C ... Patch indeces

         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
             
         IF(IFACE == 2 .OR. IFACE == 5) THEN
           KY1  = ICON(4,IP)
           KY2  = ICON(5,IP)
           KX1  = ICON(6,IP)
           KX2  = ICON(7,IP)
         ENDIF

C ... Xi-direction
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IN   = JNRE
         JN   = KNRE
         KN   = INRE
         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K = IMAX 
         KALA = 1
         KYLA = IMAX

C ... Eta-direction
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IN   = KNRE
         JN   = INRE
         KN   = JNRE
         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K = JMAX
         KALA = 1
         KYLA = JMAX

C ... Zeta direction
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IN   = INRE
         JN   = JNRE
         KN   = KNRE
         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K = KMAX
         KALA = 1
         KYLA = KMAX 
      ENDIF

C ... Loop over this patch

      IF(BACK) THEN                       ! Iclude ghostcells

      KA = (KN+K-1)*KSTR

      DO J = KY1-2,KY2+2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1-2,KX2+2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            APU2(NP) = APU1(L)               ! 3D real to 2d real
         ENDDO
      ENDDO

      ELSE IF(.NOT.BACK) THEN 

      DO K = KALA-KN,KYLA+KN

      KA = (KN+K-1)*KSTR

      DO J = KY1-2,KY2+2
         IJ = (JN+J-1)*JSTR + KA
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO I = KX1-2,KX2+2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! CELL INDEX (STRIDE=ISTR)
            NP = I + NR                      ! Patch index with LN ghost cells
            APU1(L) = WAVEH(NP)      
         ENDDO
      ENDDO

      ENDDO
      ENDIF ! BACK

      ENDIF ! SURFACE .AND. SURCHA == 'FRE'

 8000 CONTINUE

      RETURN
      END SUBROUTINE PATVSP
C
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C ----------------------------------------------------------------------      
C
