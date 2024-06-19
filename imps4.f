C *******************************************************************
C     ROUTINES FOR ALGEBRAIC MODELS                                 *
C *******************************************************************

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SYI(W1,W2,W3,W4,W5,X1,X2,X3,X4,X5,C,RO,
     2     A2X,A2Y,A2Z,A21,A22,A23,A31,A32,A33)
C
C ... MULTIPLICATION OF Y BY SY (INDEPENDENT OF THE STATE EQUATION)
C
      Z1I     = X1
      Z5I     = X5
C
C ... FROM THE PRIMITIVE (RO,U,V,W,P) TO THE CONTRAVARIANT VARIABLES
C
      A11 = A2X
      A12 = A2Y
      A13 = A2Z

      CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)

      Z2I      = A2X*X2 + A2Y*X3 + A2Z*X4
      Z3I      = A21*X2 + A22*X3 + A23*X4
      Z4I      = A31*X2 + A32*X3 + A33*X4
C
C ... FROM THE CONTRAVARIANT TO THE CHARACTERISTIC VARIABLES
C
      ROC  = RO*C
      W1   = Z1I - Z5I/RO
      W2   = ROC*Z2I + Z5I
      W3   = Z3I
      W4   = Z4I
      W5   = -ROC*Z2I + Z5I

      RETURN
      END SUBROUTINE SYI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SYM1I(X1,X2,X3,X4,X5,W1,W2,W3,W4,W5,C,RO,
     2     A2X,A2Y,A2Z,A21,A22,A23,A31,A32,A33)

C
C ... MULTIPLICATION OF Y BY SYM1 (INDEPENDENT OF THE STATE EQUATION)
C
C ... FROM THE CHARACTERISTIC TO THE CONTRAVARIANT VARIABLES
C
      YPRO    = 1./RO
      PC2     = .5*YPRO
      PROC    = PC2/C
      X1      = W1 + PC2*(W2 + W5)
      X2I     = PROC*(W2 - W5)
      X3I     = W3
      X4I     = W4
      X5   = .5*(W2 + W5)
C
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES
C
      X2   = A2X*X2I + A21*X3I + A31*X4I
      X3   = A2Y*X2I + A22*X3I + A32*X4I
      X4   = A2Z*X2I + A23*X3I + A33*X4I

      RETURN
      END SUBROUTINE SYM1I
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,
     3 IDI1,IDI2,IDI3,VIS,VIST,PR,PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,
     4 ZZZ,MAW,MAXW)

      INTEGER, PARAMETER :: IROW = 1000

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),RKSI(*)

      REAL, POINTER ::
     +     CDRO(:),CDM(:), CDN(:), CDW(:), CDE(:),
     +     XI1(:), XI2(:), XI5(:),SCINP(:),
     +     XJ1(:), XJ2(:), XJ5(:),SCINV(:), 
     +     XK1(:), XK2(:), XK5(:),DT(:)
      REAL, TARGET ::  ZZZ(MAXW)

      REAL :: YDROI(IROW),YDMI(IROW),YDNI(IROW),YDWI(IROW),YDEI(IROW),
     6        BDRO(IROW), BDM(IROW), BDN(IROW), BDW(IROW), BDE(IROW)
          
      XK1   => ZZZ( 0*MAW+1: 1*MAW); XK2  => ZZZ( 1*MAW+1: 2*MAW)
      XJ1   => ZZZ( 2*MAW+1: 3*MAW); XJ2  => ZZZ( 3*MAW+1: 4*MAW)
      XK5   => ZZZ( 4*MAW+1: 5*MAW);CDRO  => ZZZ( 5*MAW+1: 6*MAW)
      CDM   => ZZZ( 6*MAW+1: 7*MAW); CDN  => ZZZ( 7*MAW+1: 8*MAW)
      CDE   => ZZZ( 8*MAW+1: 9*MAW); CDW  => ZZZ( 9*MAW+1:10*MAW) 
       DT   => ZZZ(10*MAW+1:11*MAW); XI2  => ZZZ(11*MAW+1:12*MAW)
      XI1   => ZZZ(12*MAW+1:13*MAW);XJ5   => ZZZ(13*MAW+1:14*MAW)
      XI5   => ZZZ(14*MAW+1:15*MAW);SCINV => ZZZ(15*MAW+1:16*MAW)
      SCINP => ZZZ(16*MAW+1:17*MAW)
     
C
C ... CORRECTOR STEP OF THE LU-SGS METHOD 
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = JSTRID*ISTRID

      CDIFF   = 2.
C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
      CDRO(I) = 0.
      CDM(I)  = 0.
      CDN(I)  = 0.
      CDW(I)  = 0.
900   CDE(I)  = 0.

      DO 1000 K = 1,KMAX
      KK      = (KN+K-1)*IL
      FFK     = 2.
C     IF(K == 1) FFK = 4.        ! VAIKUTUS TESTAAMATTA
      IF(K > IDI3) FFK = 0.
      CDIFFK  = FFK*CDIFF

      DO 1100 J = 1,JMAX
      NN      = KK + (JN+J-1)*ISTRID + IN
      FFJ     = 2.
C     IF(J == 1) FFJ = 4.
      IF(J > IDI2) FFJ = 0.
      CDIFFJ  = CDIFF*FFJ

      DO 1100 I = 1,IMAX
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      FFI     = 2.
C     IF(I == 1) FFI = 4.
      IF(I > IDI1) FFI = 0.
      CDIFFI  = FFI*CDIFF

      UC      = UROT(N+1)
      RL1     = A1X(N+1)*U(N) + A1Y(N+1)*V(N) + A1Z(N+1)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+1))*(RO(N)+RO(N+1))
      T       = CDIFFI*A1(N)/RNUM*(VIST(N)+VIS(N)+VIST(N+1)+VIS(N+1))

      RIMAX   = 1.00*(ABS(RL1) + C(N))
      XI1(L)  = (RL1 + RIMAX)*.5 + T
      XI2(L)  = (RL2 + RIMAX)*.5 + T
      XI5(L)  = (RL5 + RIMAX)*.5 + T
      XIA     = RIMAX + T

      UC      = VROT(N+ISTRID)
      RL1     = A2X(N+ISTRID)*U(N)+A2Y(N+ISTRID)*V(N)+A2Z(N+ISTRID)*W(N)
     +        - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+ISTRID))*(RO(N)+RO(N+ISTRID))
      T       = CDIFFJ*A2(N)/RNUM*
     +         (VIST(N)+VIS(N)+VIST(N+ISTRID)+VIS(N+ISTRID))

      RJMAX   = 1.00*(ABS(RL1) + C(N))
      XJ1(L)  = (RL1 + RJMAX)*.5 + T
      XJ2(L)  = (RL2 + RJMAX)*.5 + T
      XJ5(L)  = (RL5 + RJMAX)*.5 + T
      XJA     = RJMAX + T

      UC      = WROT(N+IL)
      RL1     = A3X(N+IL)*U(N) + A3Y(N+IL)*V(N) + A3Z(N+IL)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
      T       = CDIFFK*A3(N)/RNUM*(VIST(N)+VIS(N)+VIST(N+IL)+VIS(N+IL))

      RKMAX   = 1.00*(ABS(RL1) + C(N))
      XK1(L)  = (RL1 + RKMAX)*.5 + T
      XK2(L)  = (RL2 + RKMAX)*.5 + T
      XK5(L)  = (RL5 + RKMAX)*.5 + T
      XKA     = RKMAX + T

      DT1     = THETA*DTL(N)
      DT(L)   = DT1
      CDRO(L) = DT1*A3(N)*CDRO(L) + VOL(N)*DRO(N)
      CDM(L)  = DT1*A3(N)*CDM(L)  + VOL(N)*DM(N)
      CDN(L)  = DT1*A3(N)*CDN(L)  + VOL(N)*DN(N)
      CDW(L)  = DT1*A3(N)*CDW(L)  + VOL(N)*DW(N)
      CDE(L)  = DT1*A3(N)*CDE(L)  + VOL(N)*DE(N)

      SCINV(L)= 1./(VOL(N)+DT1*(XIA*A1(N+1)+XJA*A2(N+ISTRID)+
     +            XKA*A3(N+IL)) + VOL(N)*RKSI(N))
      SCINP(L)= 1./(VOL(N)+DT1*(RIMAX*A1(N)+RJMAX*A2(N)+RKMAX*A3(N))
     +              + VOL(N)*RKSI(N))
1100  CONTINUE

C ... STARTING FLUXES INSIDE THE SLAB

      DO 1200  J = 1,IMAX+1
      YDROI(J) = 0.
      YDMI(J)  = 0.
      YDNI(J)  = 0.
      YDWI(J)  = 0.
1200  YDEI(J)  = 0.
      
      DO 1300  I = 1,IMAX+1
      BDRO(I) = 0.
      BDM(I)  = 0.
      BDN(I)  = 0.
      BDW(I)  = 0.
1300  BDE(I)  = 0.
      
C ... SLABWISE SWEEP. SHORTEST LOOP IN THIS CODE. SUBSTITUTION
C ... OVER m = i+j+k IS LEFT FOR A LATER EXCERCISE

      DO 1500 M = 2,IMAX+JMAX

      IALKU   = MAX0(1,M-JMAX)
      ILOPPU  = MIN0(IMAX,M-1)
      JALKU   = M - IALKU
      JLOPPU  = M - ILOPPU
      DO 1400 I = IALKU,ILOPPU
      J       = M - I
      NN      = KK + (JN+J-1)*ISTRID + IN
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      DT1     = DT(L)

      CDRO(L) = CDRO(L)+ DT1*(A1(N)*YDROI(I)+ A2(N)*BDRO(I))
      CDM(L)  = CDM(L) + DT1*(A1(N)*YDMI(I) + A2(N)*BDM(I))
      CDN(L)  = CDN(L) + DT1*(A1(N)*YDNI(I) + A2(N)*BDN(I))
      CDW(L)  = CDW(L) + DT1*(A1(N)*YDWI(I) + A2(N)*BDW(I))
      CDE(L)  = CDE(L) + DT1*(A1(N)*YDEI(I) + A2(N)*BDE(I))

      CDRO(L) = CDRO(L)*SCINV(L)
      CDM(L)  = CDM(L) *SCINV(L)
      CDN(L)  = CDN(L) *SCINV(L)
      CDW(L)  = CDW(L) *SCINV(L)
      CDE(L)  = CDE(L) *SCINP(L)

      DRO(N ) = CDRO(L)
      DM(N )  = CDM(L)
      DN(N )  = CDN(L)
      DW(N )  = CDW(L)
      DE(N )  = CDE(L)


      CALL SYI(W1,W2,W3,W4,W5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A1X(N+1),A1Y(N+1),A1Z(N+1),
     +     A21,A22,A23,A31,A32,A33)

      CALL SYI(X1,X2,X3,X4,X5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A2X(N+ISTRID),A2Y(N+ISTRID),A2Z(N+ISTRID),
     +     B21,B22,B23,B31,B32,B33)

      W1   = XI1(L)*W1
      W2   = XI2(L)*W2
      W3   = XI1(L)*W3
      W4   = XI1(L)*W4
      W5   = XI5(L)*W5
      
      X1   = XJ1(L)*X1
      X2   = XJ2(L)*X2
      X3   = XJ1(L)*X3
      X4   = XJ1(L)*X4
      X5   = XJ5(L)*X5

      CALL SYM1I(YDROX,YDMX,YDNX,YDWX,YDEX,W1,W2,W3,W4,W5,
     +     C(N),RO(N),A1X(N+1),A1Y(N+1),A1Z(N+1),
     +     A21,A22,A23,A31,A32,A33)

      YDROI(I) = YDROX ! this must be done this way to INLINE pointers
      YDMI(I)  = YDMX  ! in SGI f90 Version 7.2.1.1m
      YDNI(I)  = YDNX
      YDWI(I)  = YDWX
      YDEI(I)  = YDEX


      CALL SYM1I(BDROX,BDMX,BDNX,BDWX,BDEX,X1,X2,X3,X4,X5,
     +     C(N),RO(N),A2X(N+ISTRID),A2Y(N+ISTRID),A2Z(N+ISTRID),
     +     B21,B22,B23,B31,B32,B33)

      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX

1400  CONTINUE
      DO 1460 I = ILOPPU+1,MAX0(IALKU+1,2),-1
      YDROI(I) = YDROI(I-1)
      YDMI(I)  = YDMI(I-1) 
      YDNI(I)  = YDNI(I-1) 
      YDWI(I)  = YDWI(I-1) 
      YDEI(I)  = YDEI(I-1) 
 1460 CONTINUE

      YDROI(IALKU) = 0.
      YDMI(IALKU)  = 0.
      YDNI(IALKU)  = 0.
      YDWI(IALKU)  = 0.
      YDEI(IALKU)  = 0.
1500  CONTINUE

C ... END OF SLABWISE SWEEP. ESTABLISH FLUX CORRECTION FOR THE NEXT SLAB

      DO 1600 I = 1,JMAX*ISTRID
      N       = KK + I + JN*ISTRID
      CALL SYI(W1,W2,W3,W4,W5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A3X(N+IL),A3Y(N+IL),A3Z(N+IL),A21,A22,A23,
     +     A31,A32,A33)


      W1  = XK1(I)*W1
      W2  = XK2(I)*W2
      W3  = XK1(I)*W3
      W4  = XK1(I)*W4
      W5  = XK5(I)*W5

      CALL SYM1I(CDROX,CDMX,CDNX,CDWX,CDEX,W1,W2,W3,W4,W5,
     +     C(N),RO(N),A3X(N+IL),A3Y(N+IL),A3Z(N+IL),A21,A22,A23,
     +     A31,A32,A33)

      CDRO(I) = CDROX ! this must be done this way to INLINE pointers
      CDM(I)  = CDMX  ! in SGI f90 Version 7.2.1.1m
      CDN(I)  = CDNX
      CDW(I)  = CDWX
      CDE(I)  = CDEX

 1600 CONTINUE

1000  CONTINUE
      RETURN
      END SUBROUTINE CORRGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PREDGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,
     3 IDI1,IDI2,IDI3,VIS,VIST,PR,PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,
     4 ZZZ,MAW,MAXW)
      PARAMETER(IROW=1000)
      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),RKSI(*)

      REAL, POINTER ::
     +     CDRO(:),CDM(:), CDN(:), CDW(:), CDE(:),
     +     XI1(:), XI2(:), XI5(:),SCINP(:),
     +     XJ1(:), XJ2(:), XJ5(:),SCINV(:), 
     +     XK1(:), XK2(:), XK5(:),   DT(:)
      REAL, TARGET ::  ZZZ(MAXW)

      REAL BDRO(IROW), BDM(IROW), BDN(IROW), BDW(IROW), BDE(IROW),
     2     YDROI(IROW),YDMI(IROW),YDNI(IROW),YDWI(IROW),YDEI(IROW)
          
      XK1   => ZZZ( 0*MAW+1: 1*MAW); XK2  => ZZZ( 1*MAW+1: 2*MAW)
      XJ1   => ZZZ( 2*MAW+1: 3*MAW); XJ2  => ZZZ( 3*MAW+1: 4*MAW)
      XK5   => ZZZ( 4*MAW+1: 5*MAW);CDRO  => ZZZ( 5*MAW+1: 6*MAW)
      CDM   => ZZZ( 6*MAW+1: 7*MAW); CDN  => ZZZ( 7*MAW+1: 8*MAW)
      CDE   => ZZZ( 8*MAW+1: 9*MAW); XJ5  => ZZZ( 9*MAW+1:10*MAW)
       DT   => ZZZ(10*MAW+1:11*MAW); CDW  => ZZZ(11*MAW+1:12*MAW)
      XI1   => ZZZ(12*MAW+1:13*MAW); XI2  => ZZZ(13*MAW+1:14*MAW)
      XI5   => ZZZ(14*MAW+1:15*MAW);SCINP => ZZZ(15*MAW+1:16*MAW)
      SCINV => ZZZ(16*MAW+1:17*MAW)

     
C
C ... PREDICTOR STEP OF THE LU-SGS METHOD 
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = JSTRID*ISTRID
      CDIFF   = 2.

C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
      CDRO(I) = 0.
      CDM(I)  = 0.
      CDN(I)  = 0.
      CDW(I)  = 0.
900   CDE(I)  = 0.

      KLOP    = 1
      DO 1000 K = KMAX,KLOP,-1
      KK      = (KN+K-1)*IL
      FFK     = 2.
      IF(K == 1) FFK = 4.
      IF(K > IDI3) FFK = 0.
      CDIFFK  = FFK*CDIFF

      DO 1100 J = 1,JMAX
      NN      = KK + (JN+J-1)*ISTRID + IN
      FFJ     = 2.
      IF(J == 1) FFJ = 4.
      IF(J > IDI2) FFJ = 0.
      CDIFFJ  = CDIFF*FFJ

      DO 1100 I = 1,IMAX
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      FFI     = 2.
      IF(I == 1) FFI = 4.
      IF(I > IDI1) FFI = 0.
      CDIFFI  = FFI*CDIFF

      UC      = UROT(N)
      RL1     = A1X(N)*U(N) + A1Y(N)*V(N) + A1Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-1))*(RO(N)+RO(N-1))
      T       = CDIFFI*A1(N)/RNUM*(VIST(N)+VIS(N)+VIST(N-1)+VIS(N-1))

      RIMAX   = 1.00*(ABS(RL1) + C(N))
      XI1(L)  = (RL1 - RIMAX)*.5 - T
      XI2(L)  = (RL2 - RIMAX)*.5 - T
      XI5(L)  = (RL5 - RIMAX)*.5 - T
      XIA     = RIMAX + T

      UC      = VROT(N)
      RL1     = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-ISTRID))*(RO(N)+RO(N-ISTRID))
      T       = CDIFFJ*A2(N)/RNUM*
     +         (VIST(N)+VIS(N)+VIST(N-ISTRID)+VIS(N-ISTRID))

      RJMAX   = 1.00*(ABS(RL1) + C(N))
      XJ1(L)  = (RL1 - RJMAX)*.5 - T
      XJ2(L)  = (RL2 - RJMAX)*.5 - T
      XJ5(L)  = (RL5 - RJMAX)*.5 - T
      XJA     = RJMAX + T

      UC      = WROT(N)
      RL1     = A3X(N)*U(N) + A3Y(N)*V(N) + A3Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-IL))*(RO(N)+RO(N-IL))
      T       = CDIFFK*A3(N)/RNUM*(VIST(N)+VIS(N)+VIST(N-IL)+VIS(N-IL))

      RKMAX   = 1.00*(ABS(RL1) + C(N))
      XK1(L)  = (RL1 - RKMAX)*.5 - T
      XK2(L)  = (RL2 - RKMAX)*.5 - T
      XK5(L)  = (RL5 - RKMAX)*.5 - T
      XKA     = RKMAX + T

      DT1     = THETA*DTL(N)
      DT(L)   = DT1
      CDRO(L) = -DT1*A3(N+IL)*CDRO(L) + VOL(N)*DRO(N)
      CDM(L)  = -DT1*A3(N+IL)*CDM(L)  + VOL(N)*DM(N)
      CDN(L)  = -DT1*A3(N+IL)*CDN(L)  + VOL(N)*DN(N)
      CDW(L)  = -DT1*A3(N+IL)*CDW(L)  + VOL(N)*DW(N)
      CDE(L)  = -DT1*A3(N+IL)*CDE(L)  + VOL(N)*DE(N)
      SCINV(L)= 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N))
     +              + VOL(N)*RKSI(N))
      SCINP(L)= 1./(VOL(N)+DT1*(RIMAX*A1(N)+RJMAX*A2(N)+RKMAX*A3(N))
     +              + VOL(N)*RKSI(N))
1100  CONTINUE

C ... STARTING FLUXES INSIDE THE SLAB

      DO 1200  J = 1,IMAX+1
      YDROI(J) = 0.
      YDMI(J)  = 0.
      YDNI(J)  = 0.
      YDWI(J)  = 0.
1200  YDEI(J)  = 0.
      
      DO 1300  I = 1,IMAX+1
      BDRO(I) = 0.
      BDM(I)  = 0.
      BDN(I)  = 0.
      BDW(I)  = 0.
1300  BDE(I)  = 0.
      
C ... SLABWISE SWEEP. SHORTEST LOOP IN THIS CODE. SUBSTITUTION
C ... OVER m = i+j+k IS LEFT FOR A LATER EXCERCISE

      DO 1500 M = IMAX+JMAX,2,-1

      IALKU   = MAX0(1,M-JMAX)
      ILOPPU  = MIN0(IMAX,M-1)
      JALKU   = M - IALKU
      JLOPPU  = M - ILOPPU
      DO 1400 I = IALKU,ILOPPU
      J       = M - I
      NN      = KK + (JN+J-1)*ISTRID + IN
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      DT1     = DT(L)

      CDRO(L) = CDRO(L)- DT1*(A1(N+1)*YDROI(I)+A2(N+ISTRID)*BDRO(I))
      CDM(L)  = CDM(L) - DT1*(A1(N+1)*YDMI (I) +A2(N+ISTRID)*BDM(I))
      CDN(L)  = CDN(L) - DT1*(A1(N+1)*YDNI (I) +A2(N+ISTRID)*BDN(I))
      CDW(L)  = CDW(L) - DT1*(A1(N+1)*YDWI (I) +A2(N+ISTRID)*BDW(I))
      CDE(L)  = CDE(L) - DT1*(A1(N+1)*YDEI (I) +A2(N+ISTRID)*BDE(I))

      CDRO(L) = CDRO(L)*SCINV(L)
      CDM(L)  = CDM(L) *SCINV(L)
      CDN(L)  = CDN(L) *SCINV(L)
      CDW(L)  = CDW(L) *SCINV(L)
      CDE(L)  = CDE(L) *SCINP(L)

      DRO(N ) = CDRO(L)
      DM(N )  = CDM(L)
      DN(N )  = CDN(L)
      DW(N )  = CDW(L)
      DE(N )  = CDE(L)

      CALL SYI(W1,W2,W3,W4,W5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A1X(N),A1Y(N),A1Z(N),A21,A22,A23,A31,A32,A33)

      CALL SYI(X1,X2,X3,X4,X5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A2X(N),A2Y(N),A2Z(N),B21,B22,B23,B31,B32,B33)

      W1   = XI1(L)*W1
      W2   = XI2(L)*W2
      W3   = XI1(L)*W3
      W4   = XI1(L)*W4
      W5   = XI5(L)*W5
      
      X1   = XJ1(L)*X1
      X2   = XJ2(L)*X2
      X3   = XJ1(L)*X3
      X4   = XJ1(L)*X4
      X5   = XJ5(L)*X5

      CALL SYM1I(YDROX,YDMX,YDNX,YDWX,YDEX,W1,W2,W3,W4,W5,
     +     C(N),RO(N),A1X(N),A1Y(N),A1Z(N),A21,A22,A23,A31,A32,A33)

      YDROI(I) = YDROX   ! this must be done this way to INLINE pointers
      YDMI(I)  = YDMX    ! in SGI f90 Version 7.2.1.1m
      YDNI(I)  = YDNX
      YDWI(I)  = YDWX
      YDEI(I)  = YDEX


      CALL SYM1I(BDROX,BDMX,BDNX,BDWX,BDEX,X1,X2,X3,X4,X5,
     +     C(N),RO(N),A2X(N),A2Y(N),A2Z(N),B21,B22,B23,B31,B32,B33)

      BDRO(I) = BDROX   ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX    ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
1400  CONTINUE

      DO 1460 I = IALKU,ILOPPU-1
      YDROI(I) = YDROI(I+1)
      YDMI(I)  = YDMI(I+1) 
      YDNI(I)  = YDNI(I+1) 
      YDWI(I)  = YDWI(I+1) 
      YDEI(I)  = YDEI(I+1) 
 1460 CONTINUE

      YDROI(ILOPPU) = 0.
      YDMI(ILOPPU)  = 0.
      YDNI(ILOPPU)  = 0.
      YDWI(ILOPPU)  = 0.
      YDEI(ILOPPU)  = 0.
1500  CONTINUE

C ... END OF SLABWISE SWEEP. ESTABLISH FLUX CORRECTION FOR THE NEXT SLAB

      DO 1600 I = 1,JMAX*ISTRID
      N       = KK + I + JN*ISTRID
      CALL SYI(W1,W2,W3,W4,W5,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     C(N),RO(N),A3X(N),A3Y(N),A3Z(N),A21,A22,A23,A31,A32,A33)

      W1  = XK1(I)*W1
      W2  = XK2(I)*W2
      W3  = XK1(I)*W3
      W4  = XK1(I)*W4
      W5  = XK5(I)*W5

      CALL SYM1I(CDROX,CDMX,CDNX,CDWX,CDEX,W1,W2,W3,W4,W5,
     +     C(N),RO(N),A3X(N),A3Y(N),A3Z(N),A21,A22,A23,A31,A32,A33)

      CDRO(I) = CDROX   ! this must be done this way to INLINE pointers
      CDM(I)  = CDMX    ! in SGI f90 Version 7.2.1.1m
      CDN(I)  = CDNX
      CDW(I)  = CDWX
      CDE(I)  = CDEX

 1600 CONTINUE

1000  CONTINUE
      RETURN
      END SUBROUTINE PREDGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIAGGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,
     3 IDI1,IDI2,IDI3,VIS,VIST,PR,PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,
     4 ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),RKSI(*)
C
C ... MULTPLICATION WITH A DIAGONAL MATRIX IN LU-SGS
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN

      IL      = JSTRID*ISTRID
      NALKU   = KN*IL + JN*ISTRID + 1
      NLOPP   = (KMAX+KN-1)*IL + (JN+JMAX)*ISTRID

      DO 1200 N = NALKU,NLOPP

      UC      = UROT(N+1)
      UCONN   = A1X(N+1)*U(N) + A1Y(N+1)*V(N) + A1Z(N+1)*W(N) - UC
      R1MAXN  = 1.00*(ABS(UCONN) + C(N))

      UC      = UROT(N)
      UCONP   = A1X(N)*U(N) + A1Y(N)*V(N) + A1Z(N)*W(N) - UC
      R1MAXP  = 1.00*(ABS(UCONP) + C(N))

      UC      = VROT(N+ISTRID)
      UCONN   = A2X(N+ISTRID)*U(N)+A2Y(N+ISTRID)*V(N)+A2Z(N+ISTRID)*W(N)
     +        - UC
      R2MAXN  = 1.00*(ABS(UCONN) + C(N))

      UC      = VROT(N)
      UCONP   = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      R2MAXP  = 1.00*(ABS(UCONP) + C(N))

      UC      = WROT(N+IL)
      UCONN   = A3X(N+IL)*U(N) + A3Y(N+IL)*V(N) + A3Z(N+IL)*W(N) - UC
      R3MAXN  = 1.00*(ABS(UCONN) + C(N))

      UC      = WROT(N)
      UCONP   = A3X(N)*U(N) + A3Y(N)*V(N) + A3Z(N)*W(N) - UC
      R3MAXP  = 1.00*(ABS(UCONP) + C(N))

      DT1     = THETA*DTL(N)/VOL(N)
      SCINV   = .5*(A1(N)*R1MAXP + A1(N+1)     *R1MAXN + 
     +              A2(N)*R2MAXP + A2(N+ISTRID)*R2MAXN +
     +              A3(N)*R3MAXP + A3(N+IL)    *R3MAXN)

      DRO(N)  = DRO(N)*(1.+ DT1*SCINV + RKSI(N))
      DM(N)   = DM(N)* (1.+ DT1*SCINV + RKSI(N))
      DN(N)   = DN(N)* (1.+ DT1*SCINV + RKSI(N))
      DW(N)   = DW(N)* (1.+ DT1*SCINV + RKSI(N))
      DE(N)   = DE(N)* (1.+ DT1*SCINV + RKSI(N))
1200  CONTINUE

1000  CONTINUE
      RETURN
      END SUBROUTINE DIAGGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
     
      SUBROUTINE CORR(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 PR,PRT,IN,JN,KN,UROT,ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*)
      REAL, POINTER :: BDRO(:),BDM(:),BDN(:),BDW(:),BDE(:)
      REAL, TARGET ::  ZZZ(MAXW)


      BDRO => ZZZ( 0*MAW+1: 1*MAW);BDM => ZZZ( 1*MAW+1: 2*MAW)
      BDN  => ZZZ( 2*MAW+1: 3*MAW);BDW => ZZZ( 3*MAW+1: 4*MAW)
      BDE  => ZZZ( 4*MAW+1: 5*MAW)

C
C ... CORRECTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C
      ISTRID = IMAX + 2*IN
      ITOT   = JMAX*ISTRID
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

C ... STARTING FLUX
      DO I = 1,ITOT
         BDRO(I) = 0.
         BDM(I)  = 0.
         BDN(I)  = 0.
         BDW(I)  = 0.
         BDE(I)  = 0.
      ENDDO

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 K = 1,KMAX

      IF(K <= IDI) THEN
         CDIFF = 2.
      ELSE
         CDIFF = 0.
      ENDIF

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1400 I = IN+1,ITOT,ISTEP

      N       = JJ + I
      DT1     = THETA*DTL(N)
CC    DT1     = DTG
      BDROX = DT1*A2(N)*BDRO(I) + VOL(N)*DRO(N)
      BDMX  = DT1*A2(N)*BDM(I)  + VOL(N)*DM(N)
      BDNX  = DT1*A2(N)*BDN(I)  + VOL(N)*DN(N)
      BDWX  = DT1*A2(N)*BDW(I)  + VOL(N)*DW(N)
      BDEX  = DT1*A2(N)*BDE(I)  + VOL(N)*DE(N)
      DTX   = DT1

      CALL SYI(YDROX,YDMX,YDNX,YDWX,YDEX,BDROX,BDMX,BDNX,BDWX,BDEX,
     2     C(N),RO(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21X,A22X,A23X,
     3     A31X,A32X,A33X)

      UC      = UROT(N+IL)
      RL1     = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
      RL2     = RL1  + C(N)
      RL3     = RL1
      RL4     = RL1
      RL5     = RL1  - C(N)
      RNUM    = .5*(VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
c         IF(J == kmax) THEN
c            RNUM = .5*(VOL(N)+1.)*(RO(N)+RO(N+IL))
c         endif
c ... saattaa vaikuttaa convergenssiin.
      T    = CDIFF*A2(N+IL)/RNUM*(VIST(N)+VIS(N)+VIST(N+IL)+VIS(N+IL))

      X1X   = MAX(RL1,0.) + T
      X2X   = MAX(RL2,0.) + T
      X5X   = MAX(RL5,0.) + T
      X1AX  = ABS(RL1) + T
      X2AX  = ABS(RL2) + T
      X5AX  = ABS(RL5) + T

      YDROX = YDROX/(VOL(N) + DTX*X1AX*A2(N+IL))
      YDMX  = YDMX/(VOL(N)  + DTX*X2AX*A2(N+IL))
      YDNX  = YDNX/(VOL(N)  + DTX*X1AX*A2(N+IL))
      YDWX  = YDWX/(VOL(N)  + DTX*X1AX*A2(N+IL))
      YDEX  = YDEX/(VOL(N)  + DTX*X5AX*A2(N+IL))

      CALL SYM1I(DRO(N),DM(N),DN(N),DW(N),DE(N),YDROX,YDMX,YDNX,YDWX,
     +     YDEX,C(N),RO(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21X,A22X,A23X,
     +     A31X,A32X,A33X)

      YDROX = X1X*YDROX
      YDMX  = X2X*YDMX
      YDNX  = X1X*YDNX
      YDWX  = X1X*YDWX
      YDEX  = X5X*YDEX

      CALL SYM1I(BDROX,BDMX,BDNX,BDWX,BDEX,YDROX,YDMX,YDNX,YDWX,YDEX,
     +     C(N),RO(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21X,A22X,A23X,
     +     A31X,A32X,A33X)

      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
 1400 CONTINUE
1000  CONTINUE
      RETURN
      END SUBROUTINE CORR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PRED(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 PR,PRT,IN,JN,KN,UROT,ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*)

      REAL, POINTER :: BDRO(:),BDM(:),BDN(:),BDW(:),BDE(:)
      REAL, TARGET ::  ZZZ(MAXW)

      BDRO => ZZZ( 0*MAW+1: 1*MAW);BDM => ZZZ( 1*MAW+1: 2*MAW)
      BDN  => ZZZ( 2*MAW+1: 3*MAW);BDW => ZZZ( 3*MAW+1: 4*MAW)
      BDE  => ZZZ( 4*MAW+1: 5*MAW)

C
C ... PREDICTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C
      ISTRID = IMAX + 2*IN
      ITOT   = JMAX*ISTRID
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = JSTRID*ISTRID

C ... STARTING FLUX
      DO I = 1,ITOT
         BDRO(I) = 0.
         BDM(I)  = 0.
         BDN(I)  = 0.
         BDW(I)  = 0.
         BDE(I)  = 0.
      ENDDO

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      KLOP = 1

      DO 1000 K = KMAX,KLOP,-1

      FF      = .5
      IF(K == 1) FF = 0.25
      IF(K <= IDI) THEN
         CDIFF= 2.
      ELSE
         CDIFF= 0.
      ENDIF

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1400 I = IN+1,ITOT,ISTEP

      N       = JJ + I
      DT1     = THETA*DTL(N)
CC    DT1     = DTG

      BDROX = DT1*A2(N+IL)*BDRO(I) + VOL(N)*DRO(N)
      BDMX  = DT1*A2(N+IL)*BDM(I)  + VOL(N)*DM(N)
      BDNX  = DT1*A2(N+IL)*BDN(I)  + VOL(N)*DN(N)
      BDWX  = DT1*A2(N+IL)*BDW(I)  + VOL(N)*DW(N)
      BDEX  = DT1*A2(N+IL)*BDE(I)  + VOL(N)*DE(N)
      DTX   = DT1

      CALL SYI(YDROX,YDMX,YDNX,YDWX,YDEX,BDROX,BDMX,BDNX,BDWX,BDEX,
     2     C(N),RO(N),A2X(N),A2Y(N),A2Z(N),A21X,A22X,
     3     A23X,A31X,A32X,A33X)

      UC      = UROT(N)
      RL1     = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL3     = RL1
      RL4     = RL1
      RL5     = RL1  - C(N)
      RNUM = FF*(VOL(N)+VOL(N-IL))*(RO(N)+RO(N-IL))
c         IF(J == 1) THEN
c            RNUM = FF*(VOL(N)+1.)*(RO(N)+RO(N-IL))
c         endif
      T    = CDIFF*A2(N)/RNUM*(VIST(N)+VIS(N)+VIST(N-IL)+VIS(N-IL))

      X1X   = MAX(-RL1,0.) + T
      X2X   = MAX(-RL2,0.) + T
      X3X   = MAX(-RL3,0.) + T
      X4X   = MAX(-RL4,0.) + T
      X5X   = MAX(-RL5,0.) + T
      X1AX  = ABS(RL1) + T
      X2AX  = ABS(RL2) + T
      X5AX  = ABS(RL5) + T

      YDROX = YDROX/(VOL(N) + DTX*X1AX*A2(N))
      YDMX  = YDMX/(VOL(N)  + DTX*X2AX*A2(N))
      YDNX  = YDNX/(VOL(N)  + DTX*X1AX*A2(N))
      YDWX  = YDWX/(VOL(N)  + DTX*X1AX*A2(N))
      YDEX  = YDEX/(VOL(N)  + DTX*X5AX*A2(N))

      CALL SYM1I(DRO(N),DM(N),DN(N),DW(N),DE(N),YDROX,YDMX,YDNX,YDWX,
     +     YDEX,C(N),RO(N),A2X(N),A2Y(N),A2Z(N),A21X,A22X,A23X,
     +     A31X,A32X,A33X)

      YDROX = X1X*YDROX
      YDMX  = X2X*YDMX
      YDNX  = X3X*YDNX
      YDWX  = X4X*YDWX
      YDEX  = X5X*YDEX

      CALL SYM1I(BDROX,BDMX,BDNX,BDWX,BDEX,YDROX,YDMX,YDNX,YDWX,YDEX,
     +     C(N),RO(N),A2X(N),A2Y(N),A2Z(N),A21X,A22X,A23X,
     +     A31X,A32X,A33X)

      BDRO(I) = BDROX
      BDM(I)  = BDMX
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX

 1400 CONTINUE
1000  CONTINUE
      RETURN
      END SUBROUTINE PRED
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIAG(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P, IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 PR,PRT,IN,JN,KN,UROT)

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*)
C
C ... MULTPLICATION WITH DIAGONAL MATRIX IN DD-ADI
C
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID
      CDIFF  = 2.

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 K = 1,KMAX

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1200 I = IN+1,JMAX*ISTRID,ISTEP

      N   = JJ + I
      DTX = THETA*DTL(N)

      CALL SYI(YDROX,YDMX,YDNX,YDWX,YDEX,DRO(N),DM(N),DN(N),DW(N),
     +     DE(N),C(N),RO(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21X,A22X,
     +     A23X,A31X,A32X,A33X)

      UC      = UROT(N+IL)
      UCON    = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
      RL2N    = MAX(-(UCON+C(N)), 0.)
      RL5N    = MAX(-(UCON-C(N)), 0.)
      RL1N    = MAX(- UCON,0.)

      UC      = UROT(N)
      UCON    = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      RL2P    = MAX(UCON+C(N), 0.)
      RL5P    = MAX(UCON-C(N), 0.)
      RL1P    = MAX(UCON,0.)

      DT1     = DTX/VOL(N)
      YDROX = YDROX*(1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDMX  = YDMX* (1.+ DT1*(A2(N)*RL2P + A2(N+IL)*RL2N))
      YDNX  = YDNX* (1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDWX  = YDWX* (1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDEX  = YDEX* (1.+ DT1*(A2(N)*RL5P + A2(N+IL)*RL5N))

      CALL SYM1I(DRO(N),DM(N),DN(N),DW(N),DE(N),YDROX,YDMX,YDNX,YDWX,
     +     YDEX,C(N),RO(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21X,A22X,A23X,
     +     A31X,A32X,A33X)

1200  CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE DIAG

C *********************************************************************
C     K-EPSILON ROUTINES                                              *
C *********************************************************************

C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SYIKE(W1,W2,W3,W4,W5,W6,W7,X1,X2,X3,X4,X5,X6,X7,
     2     RK,C,RO,P,DRDP,DRDH,A2X,A2Y,A2Z,A21,A22,A23,
     3     A31,A32,A33)
C
C ... MULTIPLICATION OF Y BY SY (INDEPENDENT OF THE STATE EQUATION)
C
      Z1    = X1
      Z5    = X5
      Z6    = X6
      Z7    = X7
C
C ... FROM THE PRIMITIVE (RO,U,V,W,P) TO THE CONTRAVARIANT VARIABLES
C

      A11 = A2X
      A12 = A2Y
      A13 = A2Z

      CALL SURFVECTS(A11,A12,A13,A21,A22,A23,A31,A32,A33)


      Z2    = A2X*X2 + A2Y*X3 + A2Z*X4
      Z3    = A21*X2 + A22*X3 + A23*X4
      Z4    = A31*X2 + A32*X3 + A33*X4
C
C ... FROM THE CONTRAVARIANT TO THE CHARACTERISTIC VARIABLES
C
      DRODH = DRDH   
      DRODP = DRDP

      ROC   = RO*C
      PC2   = 1./RO
      YPRO  = PC2
      DRH   = .667*RK*YPRO**2*DRODH
      DRP   = 1.+.667*RK*YPRO*DRODP
      DPRES = Z5     ! hahahahahahahaa !

      W1   = Z1 - PC2*DPRES
      W2   = ROC*Z2 + DPRES
      W3   = Z3
      W4   = Z4
      W6   = Z6
      W7   = Z7
      W5   =-ROC*Z2 + DPRES

      RETURN
      END SUBROUTINE SYIKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SYM1KE(X1,X2,X3,X4,X5,X6,X7,W1,W2,W3,W4,W5,W6,W7,
     2     RK,C,RO,P,DRDP,DRDH,A2X,A2Y,A2Z,A21,A22,A23,A31,A32,A33)
C
C ... MULTIPLICATION OF Y BY SYM1 (INDEPENDENT OF THE STATE EQUATION)
C
C
C ... FROM THE CHARACTERISTIC TO THE CONTRAVARIANT VARIABLES
C
      DRODH = DRDH   
      DRODP = DRDP

      YPRO  = 1./RO
      PC2   = .5*YPRO
      PROC  = PC2/C
      DRH   = .667*RK*YPRO**2*DRODH
      YDRP  = 1./(1.+.667*RK*YPRO*DRODP)

      DH    = W1 + PC2*(W2 + W5)
      RKI   = RK*YPRO
      X1    = DH
      X2I   = PROC*(W2 - W5)
      X3I   = W3
      X4I   = W4

      X5    = .5*(W2+W5)
      X5    = YDRP*(X5 - DRH*RO*DH -.667*RO*W6)

      X6    = W6
      X7    = W7
C
C ... FROM THE CONTRAVARIANT TO THE PRIMITIVE VARIABLES
C

      X2    = A2X*X2I + A21*X3I + A31*X4I
      X3    = A2Y*X2I + A22*X3I + A32*X4I
      X4    = A2Z*X2I + A23*X3I + A33*X4I

      RETURN
      END SUBROUTINE SYM1KE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE COKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P,DRDP,DRDH, IMAX,JMAX,KMAX,
     3 THETA,DTG,IDI1,IDI2,IDI3,VIS,VIST,DRK,DEPS,RK,REPS,SRK,
     4 SEPS,PR,PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,ZZZ,MAW,MAXW,
     5 TRANSL,DGI,DRETI)
      PARAMETER(IROW=1000)
      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),
     5 SEPS(*),DRDP(*),DRDH(*),RKSI(*),DGI(*),DRETI(*)

      REAL, POINTER ::
     +     CDRO(:),CDM(:), CDN(:), CDW(:), CDE(:),CDK(:),CDEP(:), 
     +     XI1(:), XI2(:), XI5(:),SCINP(:),
     +     XJ1(:), XJ2(:), XJ5(:),SCINV(:), 
     +     XK1(:), XK2(:), XK5(:),DT(:),
     +     CDG(:),CDRET(:)
      REAL, TARGET ::  ZZZ(MAXW)

      REAL YDROI(IROW),YDMI(IROW),YDNI(IROW),YDWI(IROW),YDEI(IROW),
     6      BDRO(IROW), BDM(IROW), BDN(IROW), BDW(IROW), BDE(IROW),
     7      YDKI(IROW), BDK(IROW),YDEPI(IROW),
     9      BDEP(IROW), 
     +      YDGI(IROW), BDG(IROW),YDRET(IROW),BDRET(IROW)

      LOGICAL :: TRANSL
          
      XK1   => ZZZ( 0*MAW+1: 1*MAW); XK2  => ZZZ( 1*MAW+1: 2*MAW)
      XJ1   => ZZZ( 2*MAW+1: 3*MAW); XJ2  => ZZZ( 3*MAW+1: 4*MAW)
      XK5   => ZZZ( 4*MAW+1: 5*MAW);CDRO  => ZZZ( 5*MAW+1: 6*MAW)
      CDM   => ZZZ( 6*MAW+1: 7*MAW); CDN  => ZZZ( 7*MAW+1: 8*MAW)
      CDE   => ZZZ( 8*MAW+1: 9*MAW); CDW  => ZZZ( 9*MAW+1:10*MAW) 
       DT   => ZZZ(10*MAW+1:11*MAW); XI2  => ZZZ(11*MAW+1:12*MAW)
      XI1   => ZZZ(12*MAW+1:13*MAW);XJ5   => ZZZ(13*MAW+1:14*MAW)
      XI5   => ZZZ(14*MAW+1:15*MAW);SCINV => ZZZ(15*MAW+1:16*MAW)
      SCINP => ZZZ(16*MAW+1:17*MAW); CDEP => ZZZ(17*MAW+1:18*MAW)
      CDK   => ZZZ(18*MAW+1:19*MAW);CDG   => ZZZ(19*MAW+1:20*MAW)
      CDRET => ZZZ(20*MAW+1:21*MAW)

C
C ... CORRECTOR STEP OF THE LU-SGS METHOD 
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = JSTRID*ISTRID

      CDIFF   = 2.
C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
      CDRO(I) = 0.
      CDM(I)  = 0.
      CDN(I)  = 0.
      CDW(I)  = 0.
      CDE(I)  = 0.
      CDK(I)  = 0.; CDG(I) = 0.;  CDRET(I) = 0.
900   CDEP(I) = 0.

      DO 1000 K = 1,KMAX
      KK      = (KN+K-1)*IL
      FFK     = 2.
C     IF(K == 1) FFK = 4.        ! VAIKUTUS TESTAAMATTA
      IF(K > IDI3) FFK = 0.
      CDIFFK  = FFK*CDIFF

      DO 1100 J = 1,JMAX
      NN      = KK + (JN+J-1)*ISTRID + IN
      FFJ     = 2.
C     IF(J == 1) FFJ = 4.
      IF(J > IDI2) FFJ = 0.
      CDIFFJ  = CDIFF*FFJ

      DO 1100 I = 1,IMAX
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      FFI     = 2.
C     IF(I == 1) FFI = 4.
      IF(I > IDI1) FFI = 0.
      CDIFFI  = FFI*CDIFF

      UC      = UROT(N+1)
      RL1     = A1X(N+1)*U(N) + A1Y(N+1)*V(N) + A1Z(N+1)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+1))*(RO(N)+RO(N+1))
      T       = CDIFFI*A1(N)/RNUM*(VIST(N)+VIS(N)+VIST(N+1)+VIS(N+1))

      RIMAX   = 1.00*(ABS(RL1) + C(N))
      XI1(L)  = (RL1 + RIMAX)*.5 + T
      XI2(L)  = (RL2 + RIMAX)*.5 + T
      XI5(L)  = (RL5 + RIMAX)*.5 + T
      XIA     = RIMAX + T

      UC      = VROT(N+ISTRID)
      RL1     = A2X(N+ISTRID)*U(N)+A2Y(N+ISTRID)*V(N)+A2Z(N+ISTRID)*W(N)
     +        - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+ISTRID))*(RO(N)+RO(N+ISTRID))
      T       = CDIFFJ*A2(N)/RNUM*
     +         (VIST(N)+VIS(N)+VIST(N+ISTRID)+VIS(N+ISTRID))
C      APUL    = A2(N)*
C     +         (VIST(N)+VIS(N)+VIST(N+ISTRID)+VIS(N+ISTRID))

      RJMAX   = 1.00*(ABS(RL1) + C(N))
      XJ1(L)  = (RL1 + RJMAX)*.5 + T
      XJ2(L)  = (RL2 + RJMAX)*.5 + T
      XJ5(L)  = (RL5 + RJMAX)*.5 + T
      XJA     = RJMAX + T

      UC      = WROT(N+IL)
      RL1     = A3X(N+IL)*U(N) + A3Y(N+IL)*V(N) + A3Z(N+IL)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
      T       = CDIFFK*A3(N)/RNUM*(VIST(N)+VIS(N)+VIST(N+IL)+VIS(N+IL))

      RKMAX   = 1.00*(ABS(RL1) + C(N))
      XK1(L)  = (RL1 + RKMAX)*.5 + T
      XK2(L)  = (RL2 + RKMAX)*.5 + T
      XK5(L)  = (RL5 + RKMAX)*.5 + T
      XKA     = RKMAX + T

      DT1     = THETA*DTL(N)
      DT(L)   = DT1
      CDRO(L) = DT1*A3(N)*CDRO(L) + VOL(N)*DRO(N)
      CDM(L)  = DT1*A3(N)*CDM(L)  + VOL(N)*DM(N)
      CDN(L)  = DT1*A3(N)*CDN(L)  + VOL(N)*DN(N)
      CDW(L)  = DT1*A3(N)*CDW(L)  + VOL(N)*DW(N)
      CDE(L)  = DT1*A3(N)*CDE(L)  + VOL(N)*DE(N)
      CDK(L)  = DT1*A3(N)*CDK(L)  + VOL(N)*DRK(N)
      CDEP(L) = DT1*A3(N)*CDEP(L) + VOL(N)*DEPS(N)

      IF(TRANSL) THEN
         CDG(L)   = DT1*A3(N)*CDG(L)   + VOL(N)*DGI(N)
         CDRET(L) = DT1*A3(N)*CDRET(L) + VOL(N)*DRETI(N)
      ENDIF

c       SRKA    = ABS(SRK(N))/(RK(N)+1.E-8)*.5
c       SEPSA   = ABS(SEPS(N))/(REPS(N)+1.E-8)*.5
      SCINV(L)   = 1./(VOL(N)+DT1*(XIA*A1(N+1)+XJA*A2(N+ISTRID)+
     +     XKA*A3(N+IL)) + VOL(N)*RKSI(N))
      SCINP(L)= 1./(VOL(N)+DT1*(RIMAX*A1(N)+RJMAX*A2(N)+RKMAX*A3(N))
     +     + VOL(N)*RKSI(N))
c      SCINVK(L)  = 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N)
c     2            +  SRKA))
c      SCINVE(L)  = 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N)
c     2            +  SEPSA))

1100  CONTINUE

C ... STARTING FLUXES INSIDE THE SLAB

      DO 1200  J = 1,JMAX+1
      YDROI(J) = 0.
      YDMI(J)  = 0.
      YDNI(J)  = 0.
      YDWI(J)  = 0.
      YDEI(J)  = 0.
      YDKI(J)  = 0.
      YDGI(J)  = 0.
      YDRET(J) = 0.
1200  YDEPI(J) = 0.
      
      DO 1300  I = 1,IMAX+1
      BDRO(I) = 0.
      BDM(I)  = 0.
      BDN(I)  = 0.
      BDW(I)  = 0.
      BDE(I)  = 0.
      BDK(I)  = 0.
      BDG(I)  = 0.
      BDRET(I)= 0.
1300  BDEP(I) = 0.
      
C ... SLABWISE SWEEP. SHORTEST LOOP IN THIS CODE. SUBSTITUTION
C ... OVER m = i+j+k IS LEFT FOR A LATER EXCERCISE

      DO 1500 M = 2,IMAX+JMAX

      IALKU   = MAX0(1,M-JMAX)
      ILOPPU  = MIN0(IMAX,M-1)
      JALKU   = M - IALKU
      JLOPPU  = M - ILOPPU
      DO 1400 I = IALKU,ILOPPU
      J       = M - I
      NN      = KK + (JN+J-1)*ISTRID + IN
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      DT1     = DT(L)

      CDRO(L) = CDRO(L) + DT1*(A1(N)*YDROI(I) + A2(N)*BDRO(I))
      CDM(L)  = CDM(L)  + DT1*(A1(N)*YDMI(I)  + A2(N)*BDM(I))
      CDN(L)  = CDN(L)  + DT1*(A1(N)*YDNI(I)  + A2(N)*BDN(I))
      CDW(L)  = CDW(L)  + DT1*(A1(N)*YDWI(I)  + A2(N)*BDW(I))
      CDE(L)  = CDE(L)  + DT1*(A1(N)*YDEI(I)  + A2(N)*BDE(I))
      CDK(L)  = CDK(L)  + DT1*(A1(N)*YDKI(I)  + A2(N)*BDK(I))
      CDEP(L) = CDEP(L) + DT1*(A1(N)*YDEPI(I) + A2(N)*BDEP(I))

      IF(TRANSL) THEN
         CDG(L)   = CDG(L)   + DT1*(A1(N)*YDGI(I)  + A2(N)*BDG(I))
         CDRET(L) = CDRET(L) + DT1*(A1(N)*YDRET(I) + A2(N)*BDRET(I))
      ENDIF

      CDRO(L) = CDRO(L)*SCINV(L)
      CDM(L)  = CDM(L) *SCINV(L)
      CDN(L)  = CDN(L) *SCINV(L)
      CDW(L)  = CDW(L) *SCINV(L)
      CDE(L)  = CDE(L) *SCINP(L)
      CDK(L)  = CDK(L) *SCINV(L) !K
      CDEP(L) = CDEP(L)*SCINV(L) !E

      IF(TRANSL) THEN
         CDG(L)  = CDG(L)  *SCINV(L) !G
         CDRET(L)= CDRET(L)*SCINV(L) !Re_t
      ENDIF

      DRO(N ) = CDRO(L)
      DM(N )  = CDM(L)
      DN(N )  = CDN(L)
      DW(N )  = CDW(L)
      DE(N )  = CDE(L)
      DRK(N ) = CDK(L)
      DEPS(N )= CDEP(L)

      IF(TRANSL) THEN
         DGI(N ) = CDG(L)
         DRETI(N)= CDRET(L)
         W8      = DGI(N); W9 = DRETI(N) ! Does not anything, mersu
         X8      = DGI(N); X9 = DRETI(N) ! SYIKE useless here !?
      ENDIF

      CALL SYIKE(W1,W2,W3,W4,W5,W6,W7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A1X(N+1),A1Y(N+1),A1Z(N+1),A21,A22,A23,A31,A32,A33)

      CALL SYIKE(X1,X2,X3,X4,X5,X6,X7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A2X(N+ISTRID),A2Y(N+ISTRID),A2Z(N+ISTRID),
     +     B21,B22,B23,B31,B32,B33)


      W1 = XI1(L)*W1
      W2 = XI2(L)*W2
      W3 = XI1(L)*W3
      W4 = XI1(L)*W4
      W5 = XI5(L)*W5
      W6 = XI1(L)*W6
      W7 = XI1(L)*W7
      IF(TRANSL) THEN		! Matti Palin 25.09.2017
         W8 = XI1(L)*W8
         W9 = XI1(L)*W9
      ENDIF

      X1 = XJ1(L)*X1
      X2 = XJ2(L)*X2
      X3 = XJ1(L)*X3
      X4 = XJ1(L)*X4
      X5 = XJ5(L)*X5
      X6 = XJ1(L)*X6
      X7 = XJ1(L)*X7
      IF(TRANSL) THEN		! Esa 2.2.2018
         X8 = XJ1(L)*X8
         X9 = XJ1(L)*X9
      ENDIF

      CALL SYM1KE(YDROX,YDMX,YDNX,YDWX,YDEX,YDKX,YDEPX,
     +     W1,W2,W3,W4,W5,W6,W7,RK(N),C(N),RO(N),P(N),DRDP(N),
     +     DRDH(N),A1X(N+1),A1Y(N+1),A1Z(N+1),
     +     A21,A22,A23,A31,A32,A33)

      YDROI(I) = YDROX ! this must be done this way to INLINE pointers
      YDMI(I)  = YDMX  ! in SGI f90 Version 7.2.1.1m
      YDNI(I)  = YDNX
      YDWI(I)  = YDWX
      YDEI(I)  = YDEX
      YDKI(I)  = YDKX
      YDEPI(I) = YDEPX

      IF(TRANSL) THEN
         YDGI(I)  = W8
         YDRET(I) = W9
      ENDIF

      CALL SYM1KE(BDROX,BDMX,BDNX,BDWX,BDEX,BDKX,BDEPX,
     +     X1,X2,X3,X4,X5,X6,X7,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A2X(N+ISTRID),A2Y(N+ISTRID),A2Z(N+ISTRID),
     +     B21,B22,B23,B31,B32,B33)

      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
      BDK(I)  = BDKX
      BDEP(I) = BDEPX

      IF(TRANSL) THEN
         BDG(I)   = X8
         BDRET(I) = X9
      ENDIF

1400  CONTINUE

      DO 1460 I = ILOPPU+1,MAX0(IALKU+1,2),-1
      YDROI(I) = YDROI(I-1)
      YDMI(I)  = YDMI(I-1)
      YDNI(I)  = YDNI(I-1)
      YDWI(I)  = YDWI(I-1)
      YDEI(I)  = YDEI(I-1)
      YDKI(I)  = YDKI(I-1)
      YDEPI(I) = YDEPI(I-1)
 1460 CONTINUE

      IF(TRANSL) THEN
         DO I = ILOPPU+1,MAX0(IALKU+1,2),-1
            YDGI(I)  = YDGI(I-1)
            YDRET(I) = YDRET(I-1)
         ENDDO
      ENDIF

      YDROI(IALKU) = 0.
      YDMI(IALKU)  = 0.
      YDNI(IALKU)  = 0.
      YDWI(IALKU)  = 0.
      YDEI(IALKU)  = 0.
      YDKI(IALKU)  = 0.
      YDEPI(IALKU) = 0.
      IF(TRANSL) THEN
         YDGI(IALKU)  = 0.
         YDRET(IALKU) = 0.
      ENDIF
1500  CONTINUE

C ... END OF SLABWISE SWEEP. ESTABLISH FLUX CORRECTION FOR THE NEXT SLAB


      DO 1600 I = 1,JMAX*ISTRID
      N       = KK + I + JN*ISTRID
      CALL SYIKE(W1,W2,W3,W4,W5,W6,W7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A3X(N+IL),A3Y(N+IL),A3Z(N+IL),A21,A22,A23,A31,A32,A33)

      W1 = XK1(I)*W1
      W2 = XK2(I)*W2
      W3 = XK1(I)*W3
      W4 = XK1(I)*W4
      W5 = XK5(I)*W5
      W6 = XK1(I)*W6
      W7 = XK1(I)*W7
      IF(TRANSL) THEN
         W8  = XK1(I)*DGI(N)
         W9  = XK1(I)*DRETI(N)
      ENDIF

      CALL SYM1KE(CDROX,CDMX,CDNX,CDWX,CDEX,CDKX,CDEPX,
     +     W1,W2,W3,W4,W5,W6,W7,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A3X(N+IL),A3Y(N+IL),A3Z(N+IL),A21,A22,A23,A31,A32,A33)

      CDRO(I) = CDROX ! this must be done this way to INLINE pointers
      CDM(I)  = CDMX  ! in SGI f90 Version 7.2.1.1m
      CDN(I)  = CDNX
      CDW(I)  = CDWX
      CDE(I)  = CDEX
      CDK(I)  = CDKX
      CDEP(I) = CDEPX

      IF(TRANSL) THEN
         CDG(I)   = W8
         CDRET(I) = W9
      ENDIF

 1600 CONTINUE

1000  CONTINUE

      RETURN
      END SUBROUTINE COKEGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C            
      SUBROUTINE PRKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P,DRDP,DRDH, IMAX,JMAX,KMAX,
     3 THETA,DTG,IDI1,IDI2,IDI3,VIS,VIST,DRK,DEPS,RK,REPS,SRK,
     4 SEPS,PR,PRT,IN,JN,KN,UROT,VROT,WROT,RKSI,ZZZ,MAW,MAXW,
     5 TRANSL,DGI,DRETI)
      PARAMETER(IROW=1000)
      REAL X8,X9,W8,W9
      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),
     5 SEPS(*),DRDP(*),DRDH(*),RKSI(*),DGI(*),DRETI(*)

      REAL, POINTER ::
     +     CDRO(:),CDM(:), CDN(:), CDW(:), CDE(:),CDK(:),CDEP(:),
     +     XI1(:), XI2(:), XI5(:),SCINP(:),
     +     XJ1(:), XJ2(:), XJ5(:),SCINV(:), 
     +     XK1(:), XK2(:), XK5(:),   DT(:),
     +     CDG(:),CDRET(:)
      REAL, TARGET ::  ZZZ(MAXW)

      REAL YDROI(IROW),YDMI(IROW),YDNI(IROW),YDWI(IROW),YDEI(IROW),
     4      BDRO(IROW), BDM(IROW), BDN(IROW), BDW(IROW), BDE(IROW),
     7      YDKI(IROW), BDK(IROW),YDEPI(IROW),BDEP(IROW),
     9      YDGI(IROW), BDG(IROW),YDRET(IROW),BDRET(IROW)

      LOGICAL :: TRANSL
          
      XK1   => ZZZ( 0*MAW+1: 1*MAW); XK2  => ZZZ( 1*MAW+1: 2*MAW)
      XJ1   => ZZZ( 2*MAW+1: 3*MAW); XJ2  => ZZZ( 3*MAW+1: 4*MAW)
      XK5   => ZZZ( 4*MAW+1: 5*MAW);CDRO  => ZZZ( 5*MAW+1: 6*MAW)
      CDM   => ZZZ( 6*MAW+1: 7*MAW); CDN  => ZZZ( 7*MAW+1: 8*MAW)
      CDE   => ZZZ( 8*MAW+1: 9*MAW); XJ5  => ZZZ( 9*MAW+1:10*MAW)
       DT   => ZZZ(10*MAW+1:11*MAW); CDW  => ZZZ(11*MAW+1:12*MAW)
      XI1   => ZZZ(12*MAW+1:13*MAW); XI2  => ZZZ(13*MAW+1:14*MAW)
      XI5   => ZZZ(14*MAW+1:15*MAW);SCINP => ZZZ(15*MAW+1:16*MAW)
      SCINV => ZZZ(16*MAW+1:17*MAW);CDK   => ZZZ(17*MAW+1:18*MAW)
      CDEP  => ZZZ(18*MAW+1:19*MAW);CDG   => ZZZ(19*MAW+1:20*MAW)
      CDRET => ZZZ(20*MAW+1:21*MAW)

C
C ... PREDICTOR STEP OF THE LU-SGS METHOD 
C

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = JSTRID*ISTRID
      CDIFF   = 2.

C ... STARTING FLUX
      DO 900  I = 1,JMAX*ISTRID
      CDRO(I) = 0.
      CDM(I)  = 0.
      CDN(I)  = 0.
      CDW(I)  = 0.
      CDE(I)  = 0.
      CDK(I)  = 0.; CDG(I) = 0.;  CDRET(I) = 0.
900   CDEP(I) = 0.

      KLOP    = 1
      DO 1000 K = KMAX,KLOP,-1
      KK      = (KN+K-1)*IL
      FFK     = 2.
      IF(K == 1) FFK = 4.
      IF(K > IDI3) FFK = 0.
      CDIFFK  = FFK*CDIFF

      DO 1100 J = 1,JMAX
      NN      = KK + (JN+J-1)*ISTRID + IN
      FFJ     = 2.
      IF(J == 1) FFJ = 4.
      IF(J > IDI2) FFJ = 0.
      CDIFFJ  = CDIFF*FFJ

      DO 1100 I = 1,IMAX
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      FFI     = 2.
      IF(I == 1) FFI = 4.
      IF(I > IDI1) FFI = 0.
      CDIFFI  = FFI*CDIFF

      UC      = UROT(N)
      RL1     = A1X(N)*U(N) + A1Y(N)*V(N) + A1Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-1))*(RO(N)+RO(N-1))
      T       = CDIFFI*A1(N)/RNUM*(VIST(N)+VIS(N)+VIST(N-1)+VIS(N-1))

      RIMAX   = 1.00*(ABS(RL1) + C(N))
      XI1(L)  = (RL1 - RIMAX)*.5 - T
      XI2(L)  = (RL2 - RIMAX)*.5 - T
      XI5(L)  = (RL5 - RIMAX)*.5 - T
      XIA     = RIMAX + T

      UC      = VROT(N)
      RL1     = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-ISTRID))*(RO(N)+RO(N-ISTRID))
      T       = CDIFFJ*A2(N)/RNUM*
     +         (VIST(N)+VIS(N)+VIST(N-ISTRID)+VIS(N-ISTRID))

      RJMAX   = 1.00*(ABS(RL1) + C(N))
      XJ1(L)  = (RL1 - RJMAX)*.5 - T
      XJ2(L)  = (RL2 - RJMAX)*.5 - T
      XJ5(L)  = (RL5 - RJMAX)*.5 - T
      XJA     = RJMAX + T

      UC      = WROT(N)
      RL1     = A3X(N)*U(N) + A3Y(N)*V(N) + A3Z(N)*W(N) - UC
      RL2     = RL1  + C(N)
      RL5     = RL1  - C(N)
      RNUM    = (VOL(N)+VOL(N-IL))*(RO(N)+RO(N-IL))
      T       = CDIFFK*A3(N)/RNUM*(VIST(N)+VIS(N)+VIST(N-IL)+VIS(N-IL))

      RKMAX   = 1.00*(ABS(RL1) + C(N))
      XK1(L)  = (RL1 - RKMAX)*.5 - T
      XK2(L)  = (RL2 - RKMAX)*.5 - T
      XK5(L)  = (RL5 - RKMAX)*.5 - T
      XKA     = RKMAX + T

      DT1     = THETA*DTL(N)
      DT(L)   = DT1
      CDRO(L) = -DT1*A3(N+IL)*CDRO(L) + VOL(N)*DRO(N)
      CDM(L)  = -DT1*A3(N+IL)*CDM(L)  + VOL(N)*DM(N)
      CDN(L)  = -DT1*A3(N+IL)*CDN(L)  + VOL(N)*DN(N)
      CDW(L)  = -DT1*A3(N+IL)*CDW(L)  + VOL(N)*DW(N)
      CDE(L)  = -DT1*A3(N+IL)*CDE(L)  + VOL(N)*DE(N)
      CDK(L)  = -DT1*A3(N+IL)*CDK(L)  + VOL(N)*DRK(N)
      CDEP(L) = -DT1*A3(N+IL)*CDEP(L) + VOL(N)*DEPS(N)

      IF(TRANSL) THEN
         CDG(L)  = -DT1*A3(N+IL)*CDG(L)  + VOL(N)*DGI(N)
         CDRET(L)= -DT1*A3(N+IL)*CDRET(L) + VOL(N)*DRETI(N)
      ENDIF

c       SRKA    = ABS(SRK(N))/(RK(N)+1.E-8)*.5
c       SEPSA   = ABS(SEPS(N))/(REPS(N)+1.E-8)*.5
      SCINV(L)= 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N))
     +     + VOL(N)*RKSI(N))
      SCINP(L)= 1./(VOL(N)+DT1*(RIMAX*A1(N)+RJMAX*A2(N)+RKMAX*A3(N))
     +     + VOL(N)*RKSI(N))

c      SCINVK  = 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N) + SRKA))
c      SCINVE  = 1./(VOL(N)+DT1*(XIA*A1(N)+XJA*A2(N)+XKA*A3(N) + SEPSA))
1100  CONTINUE

C ... STARTING FLUXES INSIDE THE SLAB

      DO 1200  J = 1,IMAX+1
      YDROI(J) = 0.
      YDMI(J)  = 0.
      YDNI(J)  = 0.
      YDWI(J)  = 0.
      YDEI(J)  = 0.
      YDKI(J)  = 0.
      YDGI(J)  = 0.
      YDRET(J) = 0.
1200  YDEPI(J) = 0.
      
      DO 1300  I = 1,IMAX+1
      BDRO(I) = 0.
      BDM(I)  = 0.
      BDN(I)  = 0.
      BDW(I)  = 0.
      BDE(I)  = 0.
      BDK(I)  = 0.
      BDG(I)  = 0.
      BDRET(I)= 0.
1300  BDEP(I) = 0.
      
C ... SLABWISE SWEEP. SHORTEST LOOP IN THIS CODE. SUBSTITUTION
C ... OVER m = i+j+k IS LEFT FOR A LATER EXCERCISE

      DO 1500 M = IMAX+JMAX,2,-1

      IALKU   = MAX0(1,M-JMAX)
      ILOPPU  = MIN0(IMAX,M-1)
      JALKU   = M - IALKU
      JLOPPU  = M - ILOPPU
      DO 1400 I = IALKU,ILOPPU
      J       = M - I
      NN      = KK + (JN+J-1)*ISTRID + IN
      N       = NN + I                     ! GLOBAL INDEX
      L       = (J-1)*ISTRID + IN + I      ! LOCAL SLABWISE
      DT1     = DT(L)

      CDRO(L) = CDRO(L)- DT1*(A1(N+1)*YDROI(I)+A2(N+ISTRID)*BDRO(I))
      CDM(L)  = CDM(L) - DT1*(A1(N+1)*YDMI(I) +A2(N+ISTRID)*BDM(I))
      CDN(L)  = CDN(L) - DT1*(A1(N+1)*YDNI(I) +A2(N+ISTRID)*BDN(I))
      CDW(L)  = CDW(L) - DT1*(A1(N+1)*YDWI(I) +A2(N+ISTRID)*BDW(I))
      CDE(L)  = CDE(L) - DT1*(A1(N+1)*YDEI(I) +A2(N+ISTRID)*BDE(I))
      CDK(L)  = CDK(L) - DT1*(A1(N+1)*YDKI(I) +A2(N+ISTRID)*BDK(I))
      CDEP(L) = CDEP(L)- DT1*(A1(N+1)*YDEPI(I)+A2(N+ISTRID)*BDEP(I))
      IF(TRANSL) THEN
      CDG(L)  = CDG(L)  - DT1*(A1(N+1)*YDGI(I) +A2(N+ISTRID)*BDG(I))
      CDRET(L)= CDRET(L)- DT1*(A1(N+1)*YDRET(I)+A2(N+ISTRID)*BDRET(I))
      ENDIF

      CDRO(L) = CDRO(L)*SCINV(L)
      CDM(L)  = CDM(L) *SCINV(L)
      CDN(L)  = CDN(L) *SCINV(L)
      CDW(L)  = CDW(L) *SCINV(L)
      CDE(L)  = CDE(L) *SCINP(L)
      CDK(L)  = CDK(L) *SCINV(L) !K
      CDEP(L) = CDEP(L)*SCINV(L) !E
      IF(TRANSL) THEN
         CDG(L)  = CDG(L)  *SCINV(L) !G
         CDRET(L)= CDRET(L)*SCINV(L) !Re_t
      ENDIF
c      call ijkpai(n,imax,jmax,kmax,mmm,nnn,lll)
c      write(401,*) mmm,nnn,lll,n,l
c      write(301,*)i,j,k,imax,jmax,kmax

      DRO(N ) = CDRO(L)
      DM(N )  = CDM(L)
      DN(N )  = CDN(L)
      DW(N )  = CDW(L)
      DE(N )  = CDE(L)
      DRK(N ) = CDK(L)
      DEPS(N )= CDEP(L)
      IF(TRANSL) THEN
         DGI(N ) = CDG(L)
         DRETI(N)= CDRET(L)
         W8      = DGI(N); W9 = DRETI(N) ! Does not anything, mersu
         X8      = DGI(N); X9 = DRETI(N) ! SYIKE useless here !?
      ENDIF

      CALL SYIKE(W1,W2,W3,W4,W5,W6,W7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A1X(N),A1Y(N),A1Z(N),A21,A22,A23,A31,A32,A33)

      CALL SYIKE(X1,X2,X3,X4,X5,X6,X7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A2X(N),A2Y(N),A2Z(N),
     +     B21,B22,B23,B31,B32,B33)


      W1 = XI1(L)*W1
      W2 = XI2(L)*W2
      W3 = XI1(L)*W3
      W4 = XI1(L)*W4
      W5 = XI5(L)*W5
      W6 = XI1(L)*W6
      W7 = XI1(L)*W7
      IF(TRANSL) THEN		! Matti Palin 25.09.2017
         W8 = XI1(L)*W8
         W9 = XI1(L)*W9
      ENDIF

      X1 = XJ1(L)*X1
      X2 = XJ2(L)*X2
      X3 = XJ1(L)*X3
      X4 = XJ1(L)*X4
      X5 = XJ5(L)*X5
      X6 = XJ1(L)*X6
      X7 = XJ1(L)*X7
      IF(TRANSL) THEN		! Esa 2.02.2018
         X8 = XJ1(L)*X8
         X9 = XJ1(L)*X9
      ENDIF

      CALL SYM1KE(YDROX,YDMX,YDNX,YDWX,YDEX,YDKX,YDEPX,
     +     W1,W2,W3,W4,W5,W6,W7,RK(N),C(N),RO(N),P(N),DRDP(N),
     +     DRDH(N),A1X(N),A1Y(N),A1Z(N),
     +     A21,A22,A23,A31,A32,A33)

      YDROI(I) = YDROX   ! this must be done this way to INLINE pointers
      YDMI(I)  = YDMX    ! in SGI f90 Version 7.2.1.1m
      YDNI(I)  = YDNX
      YDWI(I)  = YDWX
      YDEI(I)  = YDEX
      YDKI(I)  = YDKX
      YDEPI(I) = YDEPX
      IF(TRANSL) THEN
         YDGI(I) = W8
         YDRET(I)= W9
      ENDIF

      CALL SYM1KE(BDROX,BDMX,BDNX,BDWX,BDEX,BDKX,BDEPX,
     +     X1,X2,X3,X4,X5,X6,X7,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A2X(N),A2Y(N),A2Z(N),B21,B22,B23,B31,B32,B33)

      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
      BDK(I)  = BDKX
      BDEP(I) = BDEPX
      IF(TRANSL) THEN
         BDG(I)  = X8
         BDRET(I)= X9
      ENDIF


1400  CONTINUE

      DO 1460 I = MAX0(IALKU-1,1),ILOPPU-1
      YDROI(I) = YDROI(I+1)
      YDMI(I)  = YDMI(I+1)
      YDNI(I)  = YDNI(I+1)
      YDWI(I)  = YDWI(I+1)
      YDEI(I)  = YDEI(I+1)
      YDKI(I)  = YDKI(I+1)
      YDEPI(I) = YDEPI(I+1)
 1460 CONTINUE

      IF(TRANSL) THEN
         DO I = MAX0(IALKU-1,1),ILOPPU-1
            YDGI(I)  = YDGI(I+1)
            YDRET(I) = YDRET(I+1)
         ENDDO
      ENDIF

      YDROI(ILOPPU) = 0.
      YDMI(ILOPPU)  = 0.
      YDNI(ILOPPU)  = 0.
      YDWI(ILOPPU)  = 0.
      YDEI(ILOPPU)  = 0.
      YDKI(ILOPPU)  = 0.
      YDEPI(ILOPPU) = 0.
      IF(TRANSL) THEN
         YDGI(ILOPPU)  = 0.
         YDRET(ILOPPU) = 0.
      ENDIF


1500  CONTINUE

C ... END OF SLABWISE SWEEP. ESTABLISH FLUX CORRECTION FOR THE NEXT SLAB

      DO 1600 I = 1,JMAX*ISTRID
      N       = KK + I + JN*ISTRID
      CALL SYIKE(W1,W2,W3,W4,W5,W6,W7,DRO(N),DM(N),DN(N),DW(N),DE(N),
     +     DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A3X(N),A3Y(N),A3Z(N),A21,A22,A23,A31,A32,A33)

      W1 = XK1(I)*W1
      W2 = XK2(I)*W2
      W3 = XK1(I)*W3
      W4 = XK1(I)*W4
      W5 = XK5(I)*W5
      W6 = XK1(I)*W6
      W7 = XK1(I)*W7
      IF(TRANSL) THEN
         W8 = XK1(I)*DGI(N)
         W9 = XK1(I)*DRETI(N)
      ENDIF

      CALL SYM1KE(CDROX,CDMX,CDNX,CDWX,CDEX,CDKX,CDEPX,
     +     W1,W2,W3,W4,W5,W6,W7,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     2     A3X(N),A3Y(N),A3Z(N),A21,A22,A23,A31,A32,A33)

      CDRO(I) = CDROX ! this must be done this way to INLINE pointers
      CDM(I)  = CDMX  ! in SGI f90 Version 7.2.1.1m
      CDN(I)  = CDNX
      CDW(I)  = CDWX
      CDE(I)  = CDEX
      CDK(I)  = CDKX
      CDEP(I) = CDEPX

      IF(TRANSL) THEN
         CDG(I)   = W8
         CDRET(I) = W9
      ENDIF

 1600 CONTINUE

1000  CONTINUE
      RETURN
      END SUBROUTINE PRKEGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIKEGS(DRO,DM,DN,DW,DE,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     2 A2Z,A3,A3X,A3Y,A3Z,DTL,U,V,W,C,RO,P,DRDP,DRDH, IMAX,JMAX,KMAX,
     3 THETA,VIS,VIST,DRK,DEPS,RK,REPS,SRK,SEPS,
     4 IN,JN,KN,UROT,VROT,WROT,RKSI,TRANSL,DGI,DRETI)

      IMPLICIT NONE

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),UROT(*),P(*),A1(*),A1X(*),A1Y(*),A1Z(*),A3(*),A3X(*),A3Y(*),
     4 A3Z(*),VROT(*),WROT(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),
     5 SEPS(*),DRDP(*),DRDH(*),RKSI(*),DGI(*),DRETI(*)

      REAL :: UC, UCONN, UCONP, DT1, SCINV, THETA, 
     &        R1MAXN, R1MAXP, R2MAXN, R2MAXP, R3MAXN, R3MAXP
 
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN,
     &           ISTRID, JSTRID, KSTRID, N, IL, I, K, KK

      LOGICAL :: TRANSL
C
C ... MULTIPLICATION WITH A DIAGONAL MATRIX IN LU-SGS
C
      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = JSTRID*ISTRID

      DO 1000 K = 1,KMAX
      KK      = (K+KN-1)*IL + JN*ISTRID

c      DO 1200 J = 1,JMAX
c      JJ   = (J-1)*ISTRID + IN + KK
c      DO 1200 I = 1,IMAX
c      N       = JJ + I
      DO 1200 I= 1,JMAX*ISTRID
      N       = KK + I
*      UC      = WROT(N+1)
      UC      = UROT(N+1)
      UCONN   = A1X(N+1)*U(N) + A1Y(N+1)*V(N) + A1Z(N+1)*W(N) - UC
      R1MAXN  = 1.00*(ABS(UCONN) + C(N))

      UC      = UROT(N)
      UCONP   = A1X(N)*U(N) + A1Y(N)*V(N) + A1Z(N)*W(N) - UC
      R1MAXP  = 1.00*(ABS(UCONP) + C(N))

      UC      = VROT(N+ISTRID)
      UCONN   = A2X(N+ISTRID)*U(N)+A2Y(N+ISTRID)*V(N)+A2Z(N+ISTRID)*W(N)
     +        - UC
      R2MAXN  = 1.00*(ABS(UCONN) + C(N))

      UC      = VROT(N)
*      UCONP   = A3X(N)*U(N) + A3Y(N)*V(N) + A3Z(N)*W(N) - UC
      UCONP   = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      R2MAXP  = 1.00*(ABS(UCONP) + C(N))

      UC      = WROT(N+IL)
      UCONN   = A3X(N+IL)*U(N) + A3Y(N+IL)*V(N) + A3Z(N+IL)*W(N) - UC
      R3MAXN  = 1.00*(ABS(UCONN) + C(N))

*      UC      = UROT(N)
      UC      = WROT(N)
      UCONP   = A3X(N)*U(N) + A3Y(N)*V(N) + A3Z(N)*W(N) - UC
      R3MAXP  = 1.00*(ABS(UCONP) + C(N))

c       SRKA    = ABS(SRK(N))/(RK(N)+1.E-8)*.5
c       SEPSA   = ABS(SEPS(N))/(REPS(N)+1.E-8)*.5

      DT1     = THETA*DTL(N)/VOL(N)
      SCINV   = .5*(A1(N)*R1MAXP + A1(N+1)     *R1MAXN + 
     +              A2(N)*R2MAXP + A2(N+ISTRID)*R2MAXN +
     +              A3(N)*R3MAXP + A3(N+IL)    *R3MAXN)

      DRO(N)  = DRO(N)* (1.+ DT1*SCINV + RKSI(N))
      DM(N)   = DM(N)*  (1.+ DT1*SCINV + RKSI(N))
      DN(N)   = DN(N)*  (1.+ DT1*SCINV + RKSI(N))
      DW(N)   = DW(N)*  (1.+ DT1*SCINV + RKSI(N))
      DE(N)   = DE(N)*  (1.+ DT1*SCINV + RKSI(N))
      DRK(N)  = DRK(N)* (1.+ DT1*SCINV + RKSI(N))
      DEPS(N) = DEPS(N)*(1.+ DT1*SCINV + RKSI(N))

      IF(TRANSL) THEN ! Ineffective, requires a separate loop
         DGI(N)   = DGI(N)   * (1.+ DT1*SCINV + RKSI(N))
         DRETI(N) = DRETI(N) * (1.+ DT1*SCINV + RKSI(N))
      ENDIF

1200  CONTINUE

1000  CONTINUE

      RETURN

      END SUBROUTINE DIKEGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE CORRKE(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,IN,JN,KN,UROT,DTURB,
     4 ZZZ,MAW,MAXW)

      REAL DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),A2(*),
     2 VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),P(*),
     4 DRDP(*),DRDH(*)
      INTEGER KCP(*)
      REAL, POINTER :: BDRO(:),BDM(:),BDN(:),BDW(:),BDE(:),BDK(:),
     +     BDEP(:)
      REAL, TARGET ::  ZZZ(MAXW)

      BDRO => ZZZ( 0*MAW+1: 1*MAW);BDM => ZZZ( 1*MAW+1: 2*MAW)
      BDN  => ZZZ( 2*MAW+1: 3*MAW);BDW => ZZZ( 3*MAW+1: 4*MAW)
      BDE  => ZZZ( 4*MAW+1: 5*MAW);BDK => ZZZ( 5*MAW+1: 6*MAW)
      BDEP => ZZZ( 6*MAW+1: 7*MAW)

C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... CORRECTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

      CDIFF   = 2.
      C1      = 0.5
      EPS     = 1.E-6

C ... STARTING FLUX ZERO
      DO I = 1,JMAX*ISTRID
         BDRO(I) = 0.
         BDM(I)  = 0.
         BDN(I)  = 0.
         BDW(I)  = 0.
         BDE(I)  = 0.
         BDK(I)  = 0.
         BDEP(I) = 0.
      ENDDO

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1
     
      DO 1000 K = 1,KMAX

      IF(K <= IDI) THEN
         CDIFF = 2.
      ELSE
         CDIFF = 0.
      ENDIF

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1400 I = IN+1,JMAX*ISTRID,ISTEP

         N      = JJ + I
         DT1    = THETA*DTL(N)
CC       DT1    = DTG
         BDROI  = DT1*A2(N)*BDRO(I) + VOL(N)*DRO(N)
         BDMI   = DT1*A2(N)*BDM(I)  + VOL(N)*DM(N)
         BDNI   = DT1*A2(N)*BDN(I)  + VOL(N)*DN(N)
         BDWI   = DT1*A2(N)*BDW(I)  + VOL(N)*DW(N)
         BDEI   = DT1*A2(N)*BDE(I)  + VOL(N)*DE(N)
         BDKI   = DTURB*DT1*A2(N)*BDK(I)  + VOL(N)*DRK(N)
         BDEPI  = DTURB*DT1*A2(N)*BDEP(I) + VOL(N)*DEPS(N)
         DT     = DT1

      CALL SYIKE(YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,BDROI,BDMI,BDNI,BDWI,
     +        BDEI,BDKI,BDEPI,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +        A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

         UC     = UROT(N+IL)
         RL1    = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
         RL2    = RL1  + C(N)
         RL5    = RL1  - C(N)
         RNUM   = .5*(VOL(N)+VOL(N+IL))*(RO(N)+RO(N+IL))
c         if(rnum  <= 0.) write(6,*)
c         if(n >= 927.and. n <= 934 .and. il == 60) write(76,*)
c     2   vol(n),vol(n-il),ro(n),ro(n+il),n,il,imax,jmax,kmax
c          call ijkpai(n,imax,jmax,kmax,mmm,nnn,lll)  
c           write(76,*) mmm,nnn,lll       
c         IF(J == kmax) THEN
c            RNUM = .5*(VOL(N)+1.)*(RO(N)+RO(N+IL))
c         endif
c ... some times calculation is more stable with this???? ppr 4.10.95
         T       = CDIFF*A2(N+IL)/RNUM*
     2        (VIST(N)+VIS(N) + VIST(N+IL)+VIS(N+IL))
       
         X1    = MAX(RL1,0.) + T
         X2    = MAX(RL2,0.) + T
         X5    = MAX(RL5,0.) + T
         X1A   = ABS(RL1) + T
         X2A   = ABS(RL2) + T
         X5A   = ABS(RL5) + T
C        SRKA  = ABS(SRK(N))/(RK(N)+EPS)*C1*0.
C        SEPSA = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
         SRKA  = 0.
         SEPSA = 0.
         VREL = 1.E-21 ! For very small volumes
         YDRO = YDRO/(VOL(N) + VREL + DT*X1A*A2(N+IL))
         YDM  = YDM/(VOL(N)  + VREL + DT*X2A*A2(N+IL))
         YDN  = YDN/(VOL(N)  + VREL + DT*X1A*A2(N+IL))
         YDW  = YDW/(VOL(N)  + VREL + DT*X1A*A2(N+IL))
         YDE  = YDE/(VOL(N)  + VREL + DT*X5A*A2(N+IL))
         YDK  = YDK/(VOL(N)  + VREL + DTURB*DT*(X1A*A2(N+IL)+SRKA))
         YDEP = YDEP/(VOL(N) + VREL + DTURB*DT*(X1A*A2(N+IL)+SEPSA))

      CALL SYM1KE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),YDRO,
     +        YDM,YDN,YDW,YDE,YDK,YDEP,RK(N),C(N),RO(N),P(N),DRDP(N),
     +        DRDH(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,
     +        A31,A32,A33)

         
         YDRO = X1*YDRO
         YDM  = X2*YDM
         YDN  = X1*YDN
         YDW  = X1*YDW
         YDE  = X5*YDE
         YDK  = X1*YDK
         YDEP = X1*YDEP

      CALL SYM1KE(BDROX,BDMX,BDNX,BDWX,BDEX,BDKX,BDEPX,
     +     YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,RK(N),C(N),RO(N),P(N),DRDP(N),
     +     DRDH(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,
     +     A31,A32,A33)
      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
      BDK(I)  = BDKX
      BDEP(I) = BDEPX
 1400 CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE CORRKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PREDKE(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,IN,JN,KN,UROT,DTURB,
     4 ZZZ,MAW,MAXW)

      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),P(*),
     4 DRDP(*),DRDH(*)
      INTEGER :: KCP(*)
      REAL, POINTER :: BDRO(:),BDM(:),BDN(:),BDW(:),BDE(:),BDK(:),
     +     BDEP(:)
      REAL, TARGET ::  ZZZ(MAXW)

      BDRO => ZZZ( 0*MAW+1: 1*MAW);BDM => ZZZ( 1*MAW+1: 2*MAW)
      BDN  => ZZZ( 2*MAW+1: 3*MAW);BDW => ZZZ( 3*MAW+1: 4*MAW)
      BDE  => ZZZ( 4*MAW+1: 5*MAW);BDK => ZZZ( 5*MAW+1: 6*MAW)
      BDEP => ZZZ( 6*MAW+1: 7*MAW)

C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... PREDICTOR STEP OF MACCORMACK METHOD (Y-DIRECTION)
C

      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ISTRID = IMAX + 2*IN
      IL     = JSTRID*ISTRID

      CDIFF = 2.
      C1    = 0.5
      EPS   = 1.E-6

C ... STARTING FLUX

C ... Zero boundary condition
      DO 900  I = 1,JMAX*ISTRID
         BDRO(I) = 0.
         BDM(I)  = 0.
         BDN(I)  = 0.
         BDW(I)  = 0.
         BDE(I)  = 0.
         BDK(I)  = 0.
         BDEP(I) = 0.
 900  CONTINUE

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      KLOP = 1

      DO 1000 K = KMAX,KLOP,-1

      JJ = (K+KN-1)*IL + JN*ISTRID
      FF = .5

      IF(K == 1) FF = .25
c      IF(J == 1) FF = 10000.25
c     mita helvettia. Laskut stabiilimpija kun ff = suuri.
C ... NACA laskut eivat pysy pystyssa 1.tasolla jos F1 on mukana.
C ... PPR 30.10.96
c      F1 = 0.
c      IF(J == 1) F1 = 1.

      IF(K <= IDI) THEN
         CDIFF = 2.
      ELSE
         CDIFF = 0.
      ENDIF

      DO 1400 I = IN+1,JMAX*ISTRID,ISTEP

         N       = JJ + I
         DT1     = THETA*DTL(N)
CC       DT1     = DTG
         BDROI = DT1*A2(N+IL)*BDRO(I) + VOL(N)*DRO(N)
         BDMI  = DT1*A2(N+IL)*BDM(I)  + VOL(N)*DM(N)
         BDNI  = DT1*A2(N+IL)*BDN(I)  + VOL(N)*DN(N)
         BDWI  = DT1*A2(N+IL)*BDW(I)  + VOL(N)*DW(N)
         BDEI  = DT1*A2(N+IL)*BDE(I)  + VOL(N)*DE(N)
         BDKI  = DTURB*DT1*A2(N+IL)*BDK(I)  + VOL(N)*DRK(N)
         BDEPI = DTURB*DT1*A2(N+IL)*BDEP(I) + VOL(N)*DEPS(N)
         DT    = DT1

      CALL SYIKE(YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,BDROI,BDMI,BDNI,BDWI,
     +        BDEI,BDKI,BDEPI,RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +        A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

         UC    = UROT(N)
         RL1   = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
         RL2   = RL1  + C(N)
         RL5   = RL1  - C(N)
         RNUM  = FF*(VOL(N)+VOL(N-IL))*(RO(N)+RO(N-IL))
c ... tällä testillä oli jokin merkitys!!!
c         if(rnum  == 0) write(6,*)vol(n),vol(n-il),n,il,imax,jmax,kmax
C         RNUM  = FF*(VOL(N)+VOL(N-IL)+F1)*(RO(N)+RO(N-IL))
         T     = CDIFF*A2(N)/RNUM*
     2        (VIST(N)+VIS(N) + VIST(N-IL)+VIS(N-IL))
         
         X1    = MAX(-RL1,0.) + T
         X2    = MAX(-RL2,0.) + T
         X5    = MAX(-RL5,0.) + T
         X1A   = ABS(RL1) + T
         X2A   = ABS(RL2) + T
         X5A   = ABS(RL5) + T

C        SRKA  = ABS(SRK(N))/(RK(N)+EPS)*C1*0.    ! huomaa naita ei sortata
C        SEPSA = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0. ! huomaa naita ei sortata
         SRKA  = 0.
         SEPSA = 0.
         VREL  = 1.E-21 ! For very small volumes
         YDRO  = YDRO/(VOL(N) + VREL + DT*X1A*A2(N))
         YDM   = YDM/(VOL(N)  + VREL + DT*X2A*A2(N))
         YDN   = YDN/(VOL(N)  + VREL + DT*X1A*A2(N))
         YDW   = YDW/(VOL(N)  + VREL + DT*X1A*A2(N))
         YDE   = YDE/(VOL(N)  + VREL + DT*X5A*A2(N))
         YDK   = YDK/(VOL(N)  + VREL + DTURB*DT*(X1A*A2(N)+SRKA))
         YDEP  = YDEP/(VOL(N) + VREL + DTURB*DT*(X1A*A2(N)+SEPSA))

      CALL SYM1KE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),YDRO,
     +        YDM,YDN,YDW,YDE,YDK,YDEP,RK(N),C(N),RO(N),P(N),DRDP(N),
     +        DRDH(N),A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

         YDRO = X1*YDRO
         YDM  = X2*YDM
         YDN  = X1*YDN
         YDW  = X1*YDW
         YDE  = X5*YDE
         YDK  = X1*YDK
         YDEP = X1*YDEP

      CALL SYM1KE(BDROX,BDMX,BDNX,BDWX,BDEX,BDKX,BDEPX,
     +     YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,RK(N),C(N),RO(N),P(N),DRDP(N),
     +     DRDH(N),A2X(N),A2Y(N),A2Z(N),A21,A22,A23,A31,A32,A33)

      BDRO(I) = BDROX ! this must be done this way to INLINE pointers
      BDM(I)  = BDMX  ! in SGI f90 Version 7.2.1.1m
      BDN(I)  = BDNX
      BDW(I)  = BDWX
      BDE(I)  = BDEX
      BDK(I)  = BDKX
      BDEP(I) = BDEPX
 1400 CONTINUE
1000  CONTINUE
      RETURN
      END SUBROUTINE PREDKE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE DIAGKE(DRO,DM,DN,DW,DE,A2,VOL,A2X,A2Y,A2Z,DTL,U,V,W,
     2 C,RO,P,DRDP,DRDH,IMAX,JMAX,KMAX,THETA,DTG,IDI,VIS,VIST,
     3 DRK,DEPS,RK,REPS,SRK,SEPS,PR,PRT,KCP,IN,JN,KN,UROT,DTURB,
     4 ZZZ,MAW,MAXW)

      INTEGER :: KCP(*)
      REAL :: DRO(*),DM(*),DN(*),DW(*),DE(*),VIST(*),VIS(*),DTL(*),
     2 A2(*),VOL(*),U(*),V(*),C(*),RO(*),A2X(*),A2Y(*),A2Z(*),
     3 W(*),DRK(*),DEPS(*),RK(*),REPS(*),SRK(*),SEPS(*),UROT(*),P(*),
     4 DRDP(*),DRDH(*)

C
C ... Y IS CHARACTERISTISC VARIABLE AND W IS PRIMITIVE
C
C ... MULTPLICATION WITH DIAGONAL MATRIX IN DD-ADI
C
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      ISTRID  = IMAX + 2*IN
      IL      = JSTRID*ISTRID
      CDIFF   = 2.
      C1      = 0.5
      EPS     = 1.E-6

      ISTEP = 1
      IF(IMAX == 1) ISTEP = 2*IN+1

      DO 1000 K = 1,KMAX

      JJ = (K+KN-1)*IL + JN*ISTRID

      DO 1100 I = IN+1,JMAX*ISTRID,ISTEP

      N  = JJ + I
      DT = THETA*DTL(N)

      CALL SYIKE(YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,DRO(N),DM(N),DN(N),DW(N),
     +     DE(N),DRK(N),DEPS(N),RK(N),C(N),RO(N),P(N),DRDP(N),DRDH(N),
     +     A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,A31,A32,A33)

      UC    = UROT(N+IL)
      UCON  = A2X(N+IL)*U(N) + A2Y(N+IL)*V(N) + A2Z(N+IL)*W(N) - UC
      RL2N  = MAX(-(UCON+C(N)), 0.)
      RL5N  = MAX(-(UCON-C(N)), 0.)
      RL1N  = MAX(- UCON,0.)

      UC    = UROT(N)
      UCON  = A2X(N)*U(N) + A2Y(N)*V(N) + A2Z(N)*W(N) - UC
      RL2P  = MAX(UCON+C(N), 0.)
      RL5P  = MAX(UCON-C(N), 0.)
      RL1P  = MAX(UCON,0.)

      DT1   = DT/(VOL(N)+1.E-21) ! For very small volumes
C     SRKA  = ABS(SRK(N))/(RK(N)+EPS)*C1*0.
C     SEPSA = ABS(SEPS(N))/(REPS(N)+EPS)*C1*0.
      SRKA  = 0.
      SEPSA = 0.
      YDRO  = YDRO*(1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDM   = YDM* (1.+ DT1*(A2(N)*RL2P + A2(N+IL)*RL2N))
      YDN   = YDN* (1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDW   = YDW* (1.+ DT1*(A2(N)*RL1P + A2(N+IL)*RL1N))
      YDE   = YDE* (1.+ DT1*(A2(N)*RL5P + A2(N+IL)*RL5N))
      YDK   = YDK* (1.+ DTURB*DT1*(A2(N)*RL1P + A2(N+IL)*RL1N+SRKA))
      YDEP  = YDEP*(1.+ DTURB*DT1*(A2(N)*RL1P + A2(N+IL)*RL1N+SEPSA))

      CALL SYM1KE(DRO(N),DM(N),DN(N),DW(N),DE(N),DRK(N),DEPS(N),
     +     YDRO,YDM,YDN,YDW,YDE,YDK,YDEP,RK(N),C(N),RO(N),P(N),
     2     DRDP(N),DRDH(N),A2X(N+IL),A2Y(N+IL),A2Z(N+IL),A21,A22,A23,
     3     A31,A32,A33)

1100  CONTINUE
1000  CONTINUE

      RETURN
      END SUBROUTINE DIAGKE

C *******************************************************************
C     ROUTINES FOR THE SOLID BLOCKS                                 *
C *******************************************************************

C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE AAMATG(AP,AN,AS,AE,AW,AT,AD,D1,D2,D3,RO,CP,VOL,CH,
     + A1,A2,A3,DTL,DE,IMAX,JMAX,KMAX,NTOT,IN,JN,KN)

      IMPLICIT NONE

      INTEGER IMAX,JMAX,KMAX,NTOT,IN,JN,KN
      REAL DDISTX,DDISTY,DDISTZ,RLAM
      INTEGER N,ISTR,JSTR,IL,I,J,K,II,JJ
      REAL AP(NTOT),AN(NTOT),AS(NTOT),AW(NTOT),AE(NTOT),AT(NTOT),
     + AD(NTOT),D1(*),D2(*),D3(*),RO(*),CP(*),VOL(*),CH(*),DTL(*),
     + A1(*),A2(*),A3(*),DE(*)
C
C ... DISCRETIZE POISSON EQUATION
C
      DO N = 1,NTOT
      AP(N) = 0.
      AN(N) = 0.
      AS(N) = 0.
      AE(N) = 0.
      AW(N) = 0.
      AT(N) = 0.
      AD(N) = 0.
c      DE(N) = -DE(N)*VOL(N) !- VOL(N)! In this solution system 
      ENDDO ! Antakee mersu


      ISTR = IMAX + 2*IN
      JSTR = JMAX + 2*JN
      IL   = ISTR*JSTR
 
      DO 1000 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX

C ... Mersu, laskee kahdesti

      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(CH(I+II)+CH(I+II+ISTR))
      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY

      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(CH(I+II)+CH(I+II-ISTR))
      AS(I+II)= -RLAM*A2(I+II)/DDISTY

      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(CH(I+II)+CH(I+II+1))
      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX 

      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(CH(I+II)+CH(I+II-1))
      AW(I+II)= -RLAM*A1(I+II)/DDISTX 

      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(CH(I+II)+CH(I+II+IL))
      AT(I+II)= -RLAM*A1(I+II+1)/DDISTZ

      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(CH(I+II)+CH(I+II-IL))
      AD(I+II)= -RLAM*A1(I+II)/DDISTZ

      AP(I+II)= -AN(I+II) - AS(I+II) - AE(I+II) - AW(I+II) - AT(I+II)
     +          -AD(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
1000  CONTINUE

C ... Boundary cells (Von Neumann condition for the implicit stage)

      K       = 0
      JJ      = (KN+K-1)*IL
      DO 2000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2000 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II+IL))
      RLAM    = .5*(CH(I+II)+CH(I+II+IL))
      AT(I+II)= -RLAM*A3(I+II+1)/DDISTZ
      AP(I+II)= -AT(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2000  CONTINUE
      K       = KMAX+1
      JJ      = (KN+K-1)*IL
      DO 2100 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2100 I = 0,IMAX+1
      DDISTZ  = .5*(D3(I+II)+D3(I+II-IL))
      RLAM    = .5*(CH(I+II)+CH(I+II-IL))
      AD(I+II)= -RLAM*A3(I+II)/DDISTZ
      AP(I+II)= -AD(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2100  CONTINUE
      DO 2200 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2200 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II+ISTR))
      RLAM    = .5*(CH(I+II)+CH(I+II+ISTR))
      AN(I+II)= -RLAM*A2(I+II+ISTR)/DDISTY
      AP(I+II)= -AN(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2200  CONTINUE
      DO 2300 K = 0,KMAX+1
      JJ      = (KN+K-1)*IL
      J       = JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 2300 I = 0,IMAX+1
      DDISTY  = .5*(D2(I+II)+D2(I+II-ISTR))
      RLAM    = .5*(CH(I+II)+CH(I+II-ISTR))
      AS(I+II)= -RLAM*A2(I+II)/DDISTY
      AP(I+II)= -AS(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2300  CONTINUE
      DO 2400 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2400 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      DDISTX  = .5*(D1(I+II)+D1(I+II+1))
      RLAM    = .5*(CH(I+II)+CH(I+II+1))
      AE(I+II)= -RLAM*A1(I+II+1)/DDISTX
      AP(I+II)= -AE(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2400  CONTINUE
      DO 2500 K = 1,KMAX
      JJ      = (KN+K-1)*IL
      DO 2500 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX+1
      DDISTX  = .5*(D1(I+II)+D1(I+II-1))
      RLAM    = .5*(CH(I+II)+CH(I+II-1))
      AW(I+II)= -RLAM*A1(I+II)/DDISTX
      AP(I+II)= -AW(I+II) + RO(I+II)*CP(I+II)*VOL(I+II)/DTL(I+II)
2500  CONTINUE
c      CALL IPRINT1(AP,IMAX,JMAX,KMAX,IN,JN,KN,631)   
c      CALL IPRINT1(AN,IMAX,JMAX,KMAX,IN,JN,KN,632)   
c      CALL IPRINT1(AS,IMAX,JMAX,KMAX,IN,JN,KN,633)   
c      CALL IPRINT1(AE,IMAX,JMAX,KMAX,IN,JN,KN,634)   
c      CALL IPRINT1(AW,IMAX,JMAX,KMAX,IN,JN,KN,635)   
c      CALL IPRINT1(AT,IMAX,JMAX,KMAX,IN,JN,KN,636)   
c      CALL IPRINT1(AD,IMAX,JMAX,KMAX,IN,JN,KN,637)   
c      CALL IPRINT1(CH,IMAX,JMAX,KMAX,IN,JN,KN,731)   
c      CALL IPRINT1(VOL,IMAX,JMAX,KMAX,IN,JN,KN,732)   
c      CALL IPRINT1(D1,IMAX,JMAX,KMAX,IN,JN,KN,733)   
c      CALL IPRINT1(D2,IMAX,JMAX,KMAX,IN,JN,KN,734)   
c      CALL IPRINT1(D3,IMAX,JMAX,KMAX,IN,JN,KN,735)   
c      CALL IPRINT1(A1,IMAX,JMAX,KMAX,IN,JN,KN,736)   
c      CALL IPRINT1(A2,IMAX,JMAX,KMAX,IN,JN,KN,737)   
c      CALL IPRINT1(A3,IMAX,JMAX,KMAX,IN,JN,KN,738)   
      END SUBROUTINE AAMATG


C *******************************************************************
C     ROUTINES FOR AUXILIARY PRINTING (UNUSED, FOR TESTING)         *
C *******************************************************************

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE IPRINT0(DTEMPP,IMAX,JMAX,KMAX,IN,JN,KN,LAITE)

      REAL :: DTEMPP(*)

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR

      DO 8000 K = 1,KMAX        !0,KMAX+1
      WRITE(LAITE,*) 'K = ',K
      WRITE(LAITE,*) '_______'
      JJ      = (KN+K-1)*IL
      DO 8000 J = 0,JMAX+1
      WRITE(LAITE,*) 'J = ',J
      II      = (JN+J-1)*ISTR + JJ + IN
      WRITE(LAITE,9000) (REAL(DTEMPP(I+II),4),I=0,IMAX+1)
9000  FORMAT(/12E11.4)
8000  CONTINUE

      RETURN
      END SUBROUTINE IPRINT0
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE IPRINT1(DTEMPP,IMAX,JMAX,KMAX,IN,JN,KN,LAITE)

      REAL :: DTEMPP(*)

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR
      DO 8000 K = 1,KMAX  !0,KMAX+1
      WRITE(LAITE,*) 'K = ',K
      WRITE(LAITE,*) '_______'
      JJ      = (KN+K-1)*IL
      DO 8000 J = 0,JMAX+1
      WRITE(LAITE,*) 'J = ',J
      II      = (JN+J-1)*ISTR + JJ + IN
      WRITE(LAITE,9000) (DTEMPP(I+II),I=0,IMAX+1)
9000  FORMAT(/12E11.4)
8000  CONTINUE

      RETURN
      END SUBROUTINE IPRINT1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C            
      SUBROUTINE IPRINT2(AP,AN,AS,AE,AW,AT,AD,Q,
     + IMAX,JMAX,KMAX,IN,JN,KN,LAITE)

      REAL :: AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),Q(*)

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR
      LAITE2 = LAITE ! Antakee mersu, muuten sekoamisvaara
        
      DO 8000 K = 0,KMAX+1
      WRITE(LAITE2,*) 'K = ',K
      WRITE(LAITE2,*) '_______'
      JJ      = (KN+K-1)*IL
      DO 8000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
       
      WRITE(LAITE2,*)
c       if(ilpo == 901) write(200,*) laite
      WRITE(LAITE2,*) 'AP, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AP(I+II),I=0,IMAX+1)
c       if(ilpo == 901) write(200,*) laite,'ap tehty'
      WRITE(LAITE2,*) 'AN, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AN(I+II),I=0,IMAX+1)
c       if(ilpo == 901) write(200,*) laite,'an tehty'
      WRITE(LAITE2,*) 'AS, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AS(I+II),I=0,IMAX+1)
c       if(ilpo == 901) write(200,*) laite,'as tehty'
      WRITE(LAITE2,*) 'AE, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AE(I+II),I=0,IMAX+1)
      WRITE(LAITE2,*) 'AW, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AW(I+II),I=0,IMAX+1)
      WRITE(LAITE2,*) 'AT, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AT(I+II),I=0,IMAX+1)
      WRITE(LAITE2,*) 'AD, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (AD(I+II),I=0,IMAX+1)
      WRITE(LAITE2,*) ' Q, J=',J
      WRITE(LAITE2,*) '--'
      WRITE(LAITE2,9000) (Q(I+II),I=0,IMAX+1)
9000  FORMAT(/12E11.4)
8000  CONTINUE

      RETURN
      END SUBROUTINE IPRINT2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C            
      SUBROUTINE IPRINT3(AP,AN,AS,AE,AW,AT,AD,Q,
     + IMAX,JMAX,KMAX,IN,JN,KN,LAITE)

      REAL :: AP(*),AN(*),AS(*),AE(*),AW(*),AT(*),AD(*),Q(*)

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR

      DO 8000 K = -1,KMAX+2
      WRITE(LAITE,*) 'K = ',K
      WRITE(LAITE,*) '_______'
      JJ      = (KN+K-1)*IL
      DO 8000 J = -1,JMAX+2
      II      = (JN+J-1)*ISTR + JJ + IN
      WRITE(LAITE,*)
      WRITE(LAITE,*) 'AP, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AP(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AN, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AN(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AS, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AS(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AE, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AE(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AW, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AW(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AT, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AT(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) 'AD, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (AD(I+II),I=-1,IMAX+2)
      WRITE(LAITE,*) ' Q, J=',J
      WRITE(LAITE,*) '--'
      WRITE(LAITE,9000) (Q(I+II),I=-1,IMAX+2)
9000  FORMAT(/12E11.4)
8000  CONTINUE

      RETURN
      END SUBROUTINE IPRINT3


C *******************************************************************
C     ROUTINES FOR THE SOLUTION OF LINEAR EQUATION SETS             *
C *******************************************************************

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOLAMG(TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAXX,JMAXX,
     + KMAXX,IN,JN,KN,MCYCLE,ITERM,ITERH,MGRID,ICONV3,ipostt,NGL,
     + ICYCLE)

C ... Multigrid solver

      USE INTEGERS, ONLY : IREPEA,NREPEA

      IMPLICIT NONE

      INTEGER :: IPOSTT,IPOST,ICONV3,MGRID,ITERH,ITERM,MCYCLE,IN,JN,KN,
     + IMAXX,JMAXX,KMAXX,M,INDD,ISO,LINE,KOKO,IPI,IPL,ISTR,JSTR,IL,I,
     + MULT,IG1,IH1,IH2,KOK,KOP,IG2,N,NGL,ICYCLE,ITERHH
      INTEGER :: IMAX(10),JMAX(10),KMAX(10),IG(10),IH(10)
      REAL :: CONV2,CONV
      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),AD(*),TEMP(*),Q(*)
      REAL, ALLOCATABLE :: APP(:),APN(:),APS(:),APE(:),APW(:),
     +      APT(:),APD(:),RESP(:),DTEMPP(:),RES(:),TEMPP(:),
     +      AA(:),BB(:),CC(:),DD(:),TEMPOLD(:)

      indd = ipostt
      IPOST   = 0 ! No post smoothing in this version (IPOST = 0)

C ... Problem dimensions and array pointers

      IMAX    = 0
      JMAX    = 0
      KMAX    = 0
      IG      = 0
      IH      = 0
      IH(1)   = 1
      IH(2)   = 1
      IG(1)   = 1
      ISO     = (IMAXX + 2*IN)*(JMAXX + 2*JN)*(KMAXX + 2*KN)
      IMAX(1) = IMAXX
      JMAX(1) = JMAXX
      KMAX(1) = KMAXX
      LINE    = AMAX0(IMAXX+2*IN,JMAXX+2*JN,KMAXX+2*KN) 

      DO 1000 M = 2,MGRID

         IMAX(M) = (IMAX(M-1)+1)/2
         JMAX(M) = (JMAX(M-1)+1)/2
         KMAX(M) = (KMAX(M-1)+1)/2
         KMAX(M) = MAX(1,KMAX(M)) ! Minimum number of slabs
         KOKO    = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
         IG(M)   = IG(M-1) + KOKO + 10
         IF(M > 2) IH(M) = IH(M-1) + KOKO + 10

1000  CONTINUE
      
      KOKO = (IMAX(MGRID)+2*IN)*(JMAX(MGRID)+2*JN)*(KMAX(MGRID)+2*KN)
      IPI  = IG(MGRID) + KOKO
      IPL  = IPI - ISO - 10

C ... Problem size determined. Start to allocate

      ALLOCATE (RES(IPI),TEMPP(IPI),DTEMPP(IPI),TEMPOLD(IPI),
     + APP(IPL),APN(IPL),APS(IPL),APE(IPL),APW(IPL),RESP(IPL),
     + APT(IPL),APD(IPL),AA(LINE),BB(LINE),CC(LINE),DD(LINE))

      ISTR    = IMAX(1) + 2*IN
      JSTR    = JMAX(1) + 2*JN
      IL      = ISTR*JSTR

      DTEMPP  = 0.
      RESP    = 0.
      RES     = 0.
      MULT    = 0
      CONV2   = 1.E9 ; CONV = 1.E9

      DO I = 1,ISO ! First level size
      TEMPP(I)   = TEMP(I)
      TEMPOLD(I) = TEMP(I)
      ENDDO

c      CALL SECAMG(RTIME)
c      RTIME1  = RTIME

C
C ... Iteration cycle begins

C ... Make an exception if a maximum change is greater than 100 (Pa)?

      DO WHILE(MULT < MCYCLE .OR. CONV > 100.)

      MULT    = MULT + 1
      IG1     = IG(1) ! Useless, starting addrres = 1
      RES     = 0.
      CALL RESAMG(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ICONV3,mult,11,NGL,CONV2)
      DTEMPP = 0.
c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,960)
c       if(iconv3 == 3) then
c      CALL LGSAMG_tul(DTEMPP,AP,AE,AW,AS,AN,AT,AD,RES,AA,BB,CC,DD,
c     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ITERM,0,2,indd)
c      else
      CALL LGSAMG(DTEMPP,AP,AE,AW,AS,AN,AT,AD,RES,AA,BB,CC,DD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ITERM,0,2,indd)
c      endif
      CALL SUBAMG(TEMPP,DTEMPP,ISO)

c       write(41+kmax(1),*) 'tihein taso kirros ',mult
c      call iprint1(tempp,imax(1),jmax(1),kmax(1),in,jn,kn,41+kmax(1))

C ... Perform multigrid?
 
      IF(MGRID > 1) THEN ! Start multigrid
      IF(MULT == 1) THEN
C ... Generate the matrix recursively on the coarse levels
      DO 4000 M = 2,MGRID
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOK     = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
      KOP     = (IMAX(M)  +2*IN)*(JMAX(M)  +2*JN)*(KMAX(M)  +2*KN)
      IF(M == 2) THEN ! Use the first level arrays
      CALL AASMA(AP,AN,AS,AE,AW,AT,AD,APP(IH2),APN(IH2),APS(IH2),
     + APE(IH2),APW(IH2),APT(IH2),APD(IH2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN,KOK,KOP,1,MULT,mcycle,M)
      ELSE ! In this way we will save some space 
      CALL AASMA(APP(IH1),APN(IH1),APS(IH1),APE(IH1),APW(IH1),APT(IH1),
     + APD(IH1),APP(IH2),APN(IH2),APS(IH2),APE(IH2),APW(IH2),APT(IH2),
     + APD(IH2),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),IN,JN,KN,KOK,KOP,1,
     + MULT,mcycle,M)
      ENDIF ! M == 2
4000  CONTINUE
      ENDIF ! MULT == 1

C ... V-cycle down

      DO 5000 M = 2,MGRID
      ITERHH  = ITERH
      IF(M >= 3) ITERHH = ITERH**2
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOP     = (IMAX(M)+2*IN)*(JMAX(M)+2*JN)*(KMAX(M)+2*KN)

C ... Recalculate the residual on the previous level.
 
      IF(M == 2) THEN ! First level residual
      RES     = 0.
      CALL RESAMG(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,12,NGL,CONV2)
c      IF(MULT == 1 .and. iconv3 == 0) THEN
c         WRITE(931,*) 'Ed.taso uusi, ICYCLE = ', ICYCLE, 'LEVEL = ',M-1
c         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,RES,
c     +   IMAX(M-1),JMAX(M-1),KMAX(M-1),IN,JN,KN,931)
c      ENDIF

c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,961)
      ELSE ! In this way we will save some space
      CALL RESAMG(RES(IG1),TEMPP(IG1),RESP(IH1),APP(IH1),APN(IH1),
     + APS(IH1),APE(IH1),APW(IH1),APT(IH1),APD(IH1),IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,13,NGL,CONV2)

      ENDIF ! M == 2

c       write(51+kmax(1),*) 'residuaali tasolla ',m-1
c      call iprint1(res(ig1),imax(m-1),jmax(m-1),kmax(m-1),in,jn,kn,
c     + 51+kmax(1))

      CALL ZERO(RESP(IH2),KOP)
      CALL RESTRI(RESP(IH2),RES(IG1),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),
     + IN,JN,KN,MULT,mcycle,APP(IH2))
c      IF(MULT == 1 .and. iconv3 == 0) THEN
c         WRITE(931,*) 'ICYCLE = ', ICYCLE, 'LEVEL = ',M
c         CALL IPRINT2(APP(IH2),APN(IH2),APS(IH2),APE(IH2),APW(IH2),
c     +   APT(IH2),APD(IH2),RESP(IH2),
c     +   IMAX(M),JMAX(M),KMAX(M),IN,JN,KN,931)
c      ENDIF

c       write(51+kmax(1),*) 'restriktioitu residuaali tasolla ',m-1
c      call iprint1(resp(ih2),imax(m),jmax(m),kmax(m),in,jn,kn,
c     + 51+kmax(1))
      CALL ZERO(TEMPP(IG2),KOP)
c      IF(iconv3 == 0) Then
c      write(8002,*) 'MULT = ',MULT
c      CALL LGSAMG_tul(TEMPP(IG2),APP(IH2),APE(IH2),APW(IH2),APS(IH2),
c     + APN(IH2),APT(IH2),APD(IH2),RESP(IH2),AA,BB,CC,DD,IMAX(M),
c     + JMAX(M),KMAX(M),IN,JN,KN,ITERHH,0,2,indd)
c       ELSE
      CALL LGSAMG(TEMPP(IG2),APP(IH2),APE(IH2),APW(IH2),APS(IH2),
     + APN(IH2),APT(IH2),APD(IH2),RESP(IH2),AA,BB,CC,DD,IMAX(M),
     + JMAX(M),KMAX(M),IN,JN,KN,ITERHH,0,2,indd)
c       ENDIF
c       write(41+kmax(1),*) 'korjaukset tasolla ',m
c      call iprint1(tempp(ig2),imax(m),jmax(m),kmax(m),in,jn,kn,
c     + 41+kmax(1))
5000  CONTINUE

C ... V-cycle up

      DO 6000 M = MGRID,2,-1
      ITERHH  = ITERH
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      CALL PROLON(TEMPP(IG1),TEMPP(IG2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN)
c ... Jälkitasoitus ei yleensä kannata, voi aktivoida halutessa

      IF(IPOST > 0) THEN
      IF(M >= 3) THEN
      CALL LGSAMG(TEMPP(IG1),APP(IH1),APE(IH1),APW(IH1),APS(IH1),
     + APN(IH1),APT(IH1),APD(IH1),RESP(IH1),AA,BB,CC,DD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,ITERHH,0,2,indd)
c      ELSE ! M = 2 (Smoothing on the first level prohibited)
c      CALL LGSAMG(TEMPP(IG1),AP(IG1),AE(IG1),AW(IG1),AS(IG1),AN(IG1),
c     + APT(IG1),APD(IG1),Q(IG1),AA,BB,CC,DD,IMAX(M-1),JMAX(M-1),
c     +KMAX(M-1),IN,JN,KN,ITERH,0,2)
      ENDIF ! M <= 3
      ENDIF ! IPOST > 0
6000  CONTINUE
      ENDIF ! MGRID > 1 (End of multigrid cycling)

c         CONV3 = SQRT(CONV3) / ((IMAX(1)+2)*(JMAX(1)+2))
c         WRITE(189,*) CONV3
c         WRITE(6,*) MULT,CONV3
C7000  CONTINUE ! End of iteration loop

      CONV = 0.

      DO N = 1,ISO ! Find the maximum change
      CONV = MAX(CONV,ABS(TEMPOLD(N)-TEMPP(N)))
      TEMPOLD(N)= TEMPP(N)
      ENDDO

c      if(icycle == 308)
c     + write(666,*)icycle,mult,5*mcycle,conv,tempp(n) ! testaa tämä. Jos conv jää
C      nakuttamaan (+/-) voidaan yrittää alirelaksoida muutosta (mersu)

      IF(MULT >=  5*MCYCLE) THEN ! 5-10 (ehdoton maksimi ajankulutukselta)
         WRITE(13,'(A,I3,2A,I6,A,I4,A,E10.2)') 
     +   '  Hermit:',MULT,' iterations made in SOLAMG.'
     +   ,' ICYCLE =',ICYCLE,', BLOCK =',NGL,', CONV =',CONV
         GO TO 7490 ! Antakee mersu, DO WHILE -rakenne ei taivu katkaisuun
      ENDIF

      ENDDO ! DO WHILE
 
7490  CONTINUE ! How to break the loop in a better way?
 
      DO 7500 N = 1,ISO ! Put the solution into the first level array
      TEMP(N)= TEMPP(N)
7500  CONTINUE

C ... Limit the maximum number of iterations to a five-fold value

      IF(MULT > MCYCLE) IREPEA(10) = IREPEA(10) + 1

      IF(IREPEA(10) >= NREPEA(10)) THEN
         WRITE(14,*)
         WRITE(14,'(" From SOLAMG: I have exeeded the given number of ",
     +   "AMG CYCLES",1X,I3," times.")') NREPEA(10)
         WRITE(14,'(" In block: ",I3," at cycle ",I5)') NGL,ICYCLE
         IREPEA(10) = 0
      ENDIF

c      RTIME   = 100.
c      CALL SECAMG(RTIME)
c      RTIME   = RTIME - RTIME1
c      RTIME1  = RTIME/MCYCLE*1000.
c      WRITE(6,*)
c      WRITE(6,9100) RTIME,RTIME1
c9100  FORMAT('KOKONAISAIKA = ',E10.5,'  KIERROS =',F7.2,' ms')
      DEALLOCATE (APP,APN,APS,APE,APW,APT,APD,RES,RESP,TEMPP,DTEMPP,
     + AA,BB,CC,DD,TEMPOLD)
c     CALL EXIT ! Joku vanha lopetuspaikka

      RETURN
      END SUBROUTINE SOLAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOLAMG_FULL(TEMP,QM,AP,AN,AS,AE,AW,AT,AD,IMAXX,JMAXX,
     + KMAXX,IN,JN,KN,MCYCLE,ITERM,ITERH,MGRID,ICONV3,ipostt,ICON,
     + NPATCH,DPDX,DPDY,DPDZ,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,A2Z,A3,A3X,
     + A3Y,A3Z,XFC,YFC,ZFC,SDI,BLKS,NGL,AMGDIVG)

C ... Multigrid solver for the pressure correction equation with
C ... skewness correction

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER IPOSTT,IPOST,ICONV3,MGRID,ITERH,ITERM,MCYCLE,IN,JN,KN,
     + IMAXX,JMAXX,KMAXX,M,INDD,ISO,LINE,KOKO,IPI,IPL,ISTR,JSTR,IL,I,
     + MULT,IG1,IH1,IH2,KOK,KOP,IG2,N,NPATCH,ngl
      INTEGER ICON(IC9,*)
      INTEGER IMAX(10),JMAX(10),KMAX(10),IG(10),IH(10)
      REAL AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),AD(*),TEMP(*),QM(*)
      REAL A1(*),A1X(*),A1Y(*),A1Z(*),A2(*),A2X(*),A2Y(*),A2Z(*),
     + A3(*),A3X(*),A3Y(*),A3Z(*),DPDX(*),DPDY(*),DPDZ(*),XFC(*),
     + YFC(*),ZFC(*),VOL(*)
      REAL, ALLOCATABLE::APP(:),APN(:),APS(:),APE(:),APW(:),
     + APT(:),APD(:),RESP(:),DTEMPP(:),RES(:),TEMPP(:),
     + AA(:),BB(:),CC(:),DD(:),QC(:),Q(:),APRINT(:)

      TYPE(DIFF_COR) SDI(*)
      TYPE(BLOCKS) BLKS

      LOGICAL AMGDIVG,PRINTOUTL

      indd = ipostt
      IPOST   = 0 ! No post smoothing in this version (IPOST = 0)

C ... Problem dimensions and array pointers

      IMAX    = 0
      JMAX    = 0
      KMAX    = 0
      IG      = 0
      IH      = 0
      IH(1)   = 1
      IH(2)   = 1
      IG(1)   = 1
      ISO     = (IMAXX + 2*IN)*(JMAXX + 2*JN)*(KMAXX + 2*KN)
      IMAX(1) = IMAXX
      JMAX(1) = JMAXX
      KMAX(1) = KMAXX
      LINE    = MAX(IMAXX+2*IN,JMAXX+2*JN,KMAXX+2*KN) 
         
      DO 1000 M = 2,MGRID

         IMAX(M) = (IMAX(M-1)+1)/2
         JMAX(M) = (JMAX(M-1)+1)/2
         KMAX(M) = (KMAX(M-1)+1)/2
         KMAX(M) = MAX(1,KMAX(M)) ! Minimum number of slabs
         KOKO    = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
         IG(M)   = IG(M-1) + KOKO + 10
         IF(M > 2) IH(M) = IH(M-1) + KOKO + 10

1000  CONTINUE
      
      KOKO = (IMAX(MGRID)+2*IN)*(JMAX(MGRID)+2*JN)*(KMAX(MGRID)+2*KN)
      IPI  = IG(MGRID) + KOKO
      IPL  = IPI - ISO - 10

C ... Problem size determined. Start to allocate

      ALLOCATE (RES(IPI),TEMPP(IPI),DTEMPP(IPI),QC(IPI),Q(IPI),
     + APP(IPL),APN(IPL),APS(IPL),APE(IPL),APW(IPL),RESP(IPL),
     + APT(IPL),APD(IPL),AA(LINE),BB(LINE),CC(LINE),DD(LINE),
     + APRINT(IPI))

      ISTR    = IMAX(1) + 2*IN
      JSTR    = JMAX(1) + 2*JN
      IL      = ISTR*JSTR

      DTEMPP  = 0.
      RESP    = 0.
      RES     = 0.
      QC      = 0.
      PRINTOUTL = .TRUE.

      DO I = 1,ISO ! First level size
      TEMPP(I)= TEMP(I)
      Q(I)    = QM(I)
      APRINT(I) =AE(I)+AW(I)+AN(I)+AS(I)+AT(I)+AD(I)
      ENDDO

c      CALL SECAMG(RTIME)
c      RTIME1  = RTIME

C
C ... Iteration cycle begins

      DO 7000 MULT = 1,MCYCLE
      IG1     = IG(1) ! Useless, starting addrres = 1
      RES     = 0.
      CALL RESAMGF(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ICONV3,mult,21,NGL)
      DTEMPP = 0.
c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,960)
      CALL LGSAMG(DTEMPP,AP,AE,AW,AS,AN,AT,AD,RES,AA,BB,CC,DD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ITERM,0,2,indd)
      CALL SUBAMG(TEMPP,DTEMPP,ISO)

C ... Calculate a gradient of pressure correction and extend on boundaries

      CALL GRADP(DPDX,DPDY,DPDZ,TEMPP,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX(1),JMAX(1),KMAX(1),IPI,IN,JN,KN)
      CALL EXTEND(DPDX,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN)
      CALL EXTEND(DPDY,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN)
      CALL EXTEND(DPDZ,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN)
      
      CALL RESCOR(Q,QC,QM,TEMPP,DPDX,DPDY,DPDZ,AP,AN,AS,AE,AW,AT,AD,
     + XFC,YFC,ZFC,SDI,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,BLKS,IMAX(1),
     + JMAX(1),KMAX(1),IN,JN,KN,ngl,mult,AMGDIVG,PRINTOUTL)

C ... Perform multigrid?
 
      IF(MGRID > 1) THEN ! Start multigrid
      IF(MULT == 1) THEN
C ... Generate the matrix recursively on the coarse levels
      DO 4000 M = 2,MGRID
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOK     = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
      KOP     = (IMAX(M)  +2*IN)*(JMAX(M)  +2*JN)*(KMAX(M)  +2*KN)
      IF(M == 2) THEN ! Use the first level arrays
      CALL AASMA(AP,AN,AS,AE,AW,AT,AD,APP(IH2),APN(IH2),APS(IH2),
     + APE(IH2),APW(IH2),APT(IH2),APD(IH2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN,KOK,KOP,1,MULT,mcycle,M)
      ELSE ! In this way we will save some space 
      CALL AASMA(APP(IH1),APN(IH1),APS(IH1),APE(IH1),APW(IH1),APT(IH1),
     + APD(IH1),APP(IH2),APN(IH2),APS(IH2),APE(IH2),APW(IH2),APT(IH2),
     + APD(IH2),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),IN,JN,KN,KOK,KOP,1,
     + MULT,mcycle,M)
      ENDIF ! M == 2
4000  CONTINUE
      ENDIF ! MULT == 1

C ... V-cycle down

      DO 5000 M = 2,MGRID
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOP     = (IMAX(M)+2*IN)*(JMAX(M)+2*JN)*(KMAX(M)+2*KN)

C ... Recalculate the residual on the previous level.

      IF(M == 2) THEN ! First level residual
      RES     = 0.
      CALL RESAMGF(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,22,NGL)
      ELSE ! In this way we will save some space
      CALL RESAMGF(RES(IG1),TEMPP(IG1),RESP(IH1),APP(IH1),APN(IH1),
     + APS(IH1),APE(IH1),APW(IH1),APT(IH1),APD(IH1),IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,23,NGL)
      ENDIF ! M == 2

      CALL ZERO(RESP(IH2),KOP)
      CALL RESTRI(RESP(IH2),RES(IG1),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),
     + IN,JN,KN,MULT,mcycle,APP(IH2))
      CALL ZERO(TEMPP(IG2),KOP)
      CALL LGSAMG(TEMPP(IG2),APP(IH2),APE(IH2),APW(IH2),APS(IH2),
     + APN(IH2),APT(IH2),APD(IH2),RESP(IH2),AA,BB,CC,DD,IMAX(M),
     + JMAX(M),KMAX(M),IN,JN,KN,ITERH,0,2,indd)
5000  CONTINUE

C ... V-cycle up

      DO 6000 M = MGRID,2,-1
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      CALL PROLON(TEMPP(IG1),TEMPP(IG2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN)
c ... Jälkitasoitus ei yleensä kannata, voi aktivoida halutessa
      IF(IPOST > 0) THEN
      IF(M >= 3) THEN
      CALL LGSAMG(TEMPP(IG1),APP(IH1),APE(IH1),APW(IH1),APS(IH1),
     + APN(IH1),APT(IH1),APD(IH1),RESP(IH1),AA,BB,CC,DD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,ITERH,0,2,indd)
c      ELSE ! M = 2 (Smoothing on the first level prohibited)
c      CALL LGSAMG(TEMPP(IG1),AP(IG1),AE(IG1),AW(IG1),AS(IG1),AN(IG1),
c     + APT(IG1),APD(IG1),Q(IG1),AA,BB,CC,DD,IMAX(M-1),JMAX(M-1),
c     +KMAX(M-1),IN,JN,KN,ITERH,0,2)
      ENDIF ! M <= 3
      ENDIF ! IPOST > 0
6000  CONTINUE
      ENDIF ! MGRID > 1 (End of multigrid cycling)

c         CONV3 = SQRT(CONV3) / ((IMAX(1)+2)*(JMAX(1)+2))
c         WRITE(189,*) CONV3
c         WRITE(6,*) MULT,CONV3
7000  CONTINUE ! End of iteration loop
 
      DO 7500 N = 1,ISO ! Put the solution into the first level array
      TEMP(N)= TEMPP(N)
7500  CONTINUE

C ... Calculate the final gradient
c       write(6666,*) 'voihan vuohi'; write(6667,*)'voihan vuohi'
      CALL MIRP(TEMPP,NPATCH,ICON,IMAX,JMAX,KMAX,IPI,IN,JN,KN)
c      CALL GRADP_tul(DPDX,DPDY,DPDZ,TEMPP,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
      CALL GRADP(DPDX,DPDY,DPDZ,TEMPP,VOL,A1,A1X,A1Y,A1Z,A2,A2X,A2Y,
     + A2Z,A3,A3X,A3Y,A3Z,IMAX(1),JMAX(1),KMAX(1),IPI,IN,JN,KN)

c      RTIME   = 100.
c      CALL SECAMG(RTIME)
c      RTIME   = RTIME - RTIME1
c      RTIME1  = RTIME/MCYCLE*1000.
c      WRITE(6,*)
c      WRITE(6,9100) RTIME,RTIME1
c9100  FORMAT('KOKONAISAIKA = ',E10.5,'  KIERROS =',F7.2,' ms')
      DEALLOCATE (APP,APN,APS,APE,APW,APT,APD,RES,RESP,TEMPP,DTEMPP,
     + AA,BB,CC,DD)
c     CALL EXIT ! Joku vanha lopetuspaikka
      END SUBROUTINE SOLAMG_FULL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESAMG(RES,TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,ICONV3,mult,ICONV,NGL,CONV)

      REAL :: RES(*),TEMP(*),AP(*),AN(*),AS(*),AE(*),AW(*),AD(*),AT(*),
     + Q(*)

C ... Calculate the residual of the equation to be solved

      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KIND    = 1
c      KIND    = 2
c      IF(KMAX == 1) KIND = 1
      IL      = ISTR*JSTR
C ... Internal points
      DO 1000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = (AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+
     + AS(IJ)*TEMP(IJ-ISTR)+AE(IJ)*TEMP(IJ+1)+AW(IJ)*TEMP(IJ-1)+
     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ))
c      if(i == 10 .and.j == 1 .and. k == 1) 
c     + write(808,*) ij,res(ij),temp(ij),q(ij),an(ij),as(ij)
1000  CONTINUE
C ... Sides
      DO 1100 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1100 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AE(IJ)*TEMP(IJ+1)+Q(IJ)
1100  CONTINUE
      DO 1200 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1200 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX + 1
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AW(IJ)*TEMP(IJ-1)+Q(IJ)
1200  CONTINUE
      DO 1300 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1300 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+Q(IJ)
1300  CONTINUE
      DO 1400 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1400 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AS(IJ)*TEMP(IJ-ISTR)+Q(IJ)
1400  CONTINUE
      K       = 0
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1500 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1500 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AT(IJ)*TEMP(IJ+IL)+Q(IJ)
1500  CONTINUE
      K       = KMAX + 1
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1600 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1600 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AD(IJ)*TEMP(IJ-IL)+Q(IJ)
c      if(mult ==  0) write(6,*) res(ij),ij,ap(ij),ad(ij),temp(ij),
c     + temp(ij-il),q(ij)
1600  CONTINUE

C ... Convergence

      IF(ICONV3 > 0) THEN
      CONV    = 0.
      CONV2   = 0.
c      if(iconv == 1 .and. mult == 50) then
c      WRITE(71,*)
c      WRITE(71,*) 'RESIDUAL'
c      CALL IPRINT1(RES,IMAX,JMAX,KMAX,IN,JN,KN,71)
c      ENDIF 
      DO 8000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 8000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 8000 I = 0,IMAX+1
      CONV    = CONV + RES(I+II)**2
      CONV2   = MAX(CONV2,ABS(RES(I+II))) !Maximum value
c      CONV2   = CONV2 + RES(I+II) ! Mielenkiintoinen ilmiö
8000  CONTINUE
9000  FORMAT(/12F8.3)
      IF(CONV > 1.E-30) THEN
      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2)*(KMAX+2))
      ELSE
      CONV    = 0.
      ENDIF
c      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2))
      IF(ICONV3 == 1) WRITE(188+NGL,*) NGL,ICONV,CONV,CONV2
c      IF(ICONV3 == 2) WRITE(189,*) CONV,CONV2
      ENDIF
      END SUBROUTINE RESAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESAMGF(RES,TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,KMAX,
     + IN,JN,KN,ICONV3,mult,ICONV,NGL)
      REAL RES(*),TEMP(*),AP(*),AN(*),AS(*),AE(*),AW(*),AD(*),AT(*),
     + Q(*)

C ... For SOLAMG_FULL
C ... Calculate the residual of the equation to be solved

      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KIND    = 1
c      KIND    = 2
c      IF(KMAX == 1) KIND = 1
      IL      = ISTR*JSTR
C ... Internal points
      DO 1000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = (AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+
     + AS(IJ)*TEMP(IJ-ISTR)+AE(IJ)*TEMP(IJ+1)+AW(IJ)*TEMP(IJ-1)+
     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ) )
c     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ) +.005*TEMP(IJ))
c      if(i == 10 .and.j == 1 .and. k == 1) 
c     + write(808,*) ij,res(ij),temp(ij),q(ij),an(ij),as(ij)
1000  CONTINUE
C ... Sides
      DO 1100 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1100 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AE(IJ)*TEMP(IJ+1)+Q(IJ)
1100  CONTINUE
      DO 1200 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1200 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX + 1
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AW(IJ)*TEMP(IJ-1)+Q(IJ)
1200  CONTINUE
      DO 1300 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1300 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+Q(IJ)
1300  CONTINUE
      DO 1400 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1400 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AS(IJ)*TEMP(IJ-ISTR)+Q(IJ)
1400  CONTINUE
      K       = 0
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1500 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1500 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AT(IJ)*TEMP(IJ+IL)+Q(IJ)
1500  CONTINUE
      K       = KMAX + 1
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1600 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1600 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AD(IJ)*TEMP(IJ-IL)+Q(IJ)
c      if(mult ==  0) write(6,*) res(ij),ij,ap(ij),ad(ij),temp(ij),
c     + temp(ij-il),q(ij)
1600  CONTINUE

C ... Convergence

      IF(ICONV3 > 0) THEN
      CONV    = 0.
      CONV2   = 0.
c      if(iconv == 1 .and. mult == 50) then
c      WRITE(71,*)
c      WRITE(71,*) 'RESIDUAL'
c      CALL IPRINT1(RES,IMAX,JMAX,KMAX,IN,JN,KN,71)
c      ENDIF 
      DO 8000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 8000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 8000 I = 0,IMAX+1
      CONV    = CONV + RES(I+II)**2
      CONV2   = MAX(CONV2,ABS(RES(I+II))) !Maximum value
c      CONV2   = CONV2 + RES(I+II) ! Mielenkiintoinen ilmiö
8000  CONTINUE
9000  FORMAT(/12F8.3)
      IF(CONV > 1.E-30) THEN
      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2)*(KMAX+2))
      ELSE
      CONV    = 0.
      ENDIF
      IF(ICONV3 == 1) WRITE(188+NGL,*) NGL,ICONV,CONV,CONV2
c      IF(ICONV3 == 2) WRITE(189,*) CONV,CONV2
      ENDIF
      END SUBROUTINE RESAMGF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOLAMG_p(TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAXX,JMAXX,
     + KMAXX,IN,JN,KN,MCYCLE,ITERM,ITERH,MGRID,ICONV3,ipostt,ngl)

C ... Multigrid solver

      IMPLICIT NONE

      INTEGER :: IPOSTT,IPOST,ICONV3,MGRID,ITERH,ITERM,MCYCLE,IN,JN,KN,
     + IMAXX,JMAXX,KMAXX,M,INDD,ISO,LINE,KOKO,IPI,IPL,ISTR,JSTR,IL,I,
     + MULT,IG1,IH1,IH2,KOK,KOP,IG2,N,ngl
      INTEGER :: IMAX(10),JMAX(10),KMAX(10),IG(10),IH(10)
      REAL    :: CONV2
      REAL    :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),AD(*),TEMP(*),Q(*)
      REAL, ALLOCATABLE::APP(:),APN(:),APS(:),APE(:),APW(:),
     +      APT(:),APD(:),RESP(:),DTEMPP(:),RES(:),TEMPP(:),
     +      AA(:),BB(:),CC(:),DD(:)

      indd = ipostt
      IPOST   = 0 ! No post smoothing in this version (IPOST = 0)

C ... Problem dimensions and array pointers

      IMAX    = 0
      JMAX    = 0
      KMAX    = 0
      IG      = 0
      IH      = 0
      IH(1)   = 1
      IH(2)   = 1
      IG(1)   = 1
      ISO     = (IMAXX + 2*IN)*(JMAXX + 2*JN)*(KMAXX + 2*KN)
      IMAX(1) = IMAXX
      JMAX(1) = JMAXX
      KMAX(1) = KMAXX
      LINE    = AMAX0(IMAXX+2*IN,JMAXX+2*JN,KMAXX+2*KN) 

      DO 1000 M = 2,MGRID
      IMAX(M) = IMAX(M-1)/2
      JMAX(M) = JMAX(M-1)/2
      KMAX(M) = KMAX(M-1)/2
      KMAX(M) = AMAX0(1,KMAX(M)) ! Minimum number of slabs
      KOKO    = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
      IG(M)   = IG(M-1) + KOKO + 10
      IF(M > 2) IH(M) = IH(M-1) + KOKO + 10
1000  CONTINUE
      
      KOKO    = (IMAX(MGRID)+2*IN)*(JMAX(MGRID)+2*JN)*
     + (KMAX(MGRID)+2*KN)
      IPI     = IG(MGRID) + KOKO
      IPL     = IPI - ISO - 10

C ... Problem size determined. Start to allocate
c         if(mgrid == 3) then
c         write(645,*) imaxx,jmaxx,kmaxx
c         write(645,*) imax(1),imax(2),imax(3),jmax(1),jmax(2),jmax(3),
c     + kmax(1),kmax(2),kmax(3)
c         write(645,*) ig(1),ig(2),ig(3),ih(1),ih(2),ih(3),ipi,ipl
c         endif
      ALLOCATE (RES(IPI),TEMPP(IPI),DTEMPP(IPI),
     + APP(IPL),APN(IPL),APS(IPL),APE(IPL),APW(IPL),RESP(IPL),
     + APT(IPL),APD(IPL),AA(LINE),BB(LINE),CC(LINE),DD(LINE))

      ISTR    = IMAX(1) + 2*IN
      JSTR    = JMAX(1) + 2*JN
      IL      = ISTR*JSTR

      DTEMPP  = 0.
      RESP    = 0.
      RES     = 0.

      DO I = 1,ISO ! First level size
      TEMPP(I)= TEMP(I)
      ENDDO

c      CALL SECAMG(RTIME)
c      RTIME1  = RTIME

C
C ... Iteration cycle begins

      DO 7000 MULT = 1,MCYCLE
      IG1     = IG(1) ! Useless, starting addrres = 1
      RES     = 0.
      CALL RESAMG(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,1,mult,ICONV3,NGL,CONV2)
      DTEMPP = 0.
c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,960)
      CALL LGSAMG(DTEMPP,AP,AE,AW,AS,AN,AT,AD,RES,AA,BB,CC,DD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ITERM,0,2,indd)
      CALL SUBAMG(TEMPP,DTEMPP,ISO)
c       write(41+kmax(1),*) 'tihein taso kirros ',mult
c      call iprint1(tempp,imax(1),jmax(1),kmax(1),in,jn,kn,41+kmax(1))

C ... Perform multigrid?
 
      IF(MGRID > 1) THEN ! Start multigrid
      IF(MULT == 1) THEN
C ... Generate the matrix recursively on the coarse levels
      DO 4000 M = 2,MGRID
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOK     = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
      KOP     = (IMAX(M)  +2*IN)*(JMAX(M)  +2*JN)*(KMAX(M)  +2*KN)
      IF(M == 2) THEN ! Use the first level arrays
      CALL AASMA(AP,AN,AS,AE,AW,AT,AD,APP(IH2),APN(IH2),APS(IH2),
     + APE(IH2),APW(IH2),APT(IH2),APD(IH2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN,KOK,KOP,1,MULT,mcycle,M)
         CALL IPRINT2(AP,AN,AS,AE,AW,AT,AD,q,
     +   IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,600)
         write(600,*)'==============================================='
         CALL IPRINT2(APP(IH2),APN(IH2),APS(IH2),
     +   APE(IH2),APW(IH2),APT(IH2),APD(IH2),q,
     +   IMAX(2),JMAX(2),KMAX(2),IN,JN,KN,601)
         write(601,*)'==============================================='

      ELSE ! In this way we will save some space 
      CALL AASMA(APP(IH1),APN(IH1),APS(IH1),APE(IH1),APW(IH1),APT(IH1),
     + APD(IH1),APP(IH2),APN(IH2),APS(IH2),APE(IH2),APW(IH2),APT(IH2),
     + APD(IH2),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),IN,JN,KN,KOK,KOP,1,
     + MULT,mcycle,M)
         CALL IPRINT2(APP(IH2),APN(IH2),APS(IH2),
     +   APE(IH2),APW(IH2),APT(IH2),APD(IH2),q,
     +   IMAX(3),JMAX(3),KMAX(3),IN,JN,KN,602)
         write(602,*)'==============================================='
      ENDIF ! M == 2
4000  CONTINUE
      ENDIF ! MULT == 1

C ... V-cycle down

      DO 5000 M = 2,MGRID
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOP     = (IMAX(M)+2*IN)*(JMAX(M)+2*JN)*(KMAX(M)+2*KN)

C ... Recalculate the residual on the previous level.

      IF(M == 2) THEN ! First level residual
      RES     = 0.
      CALL RESAMG(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,ICONV3,NGL,CONV2)
      IF(NGL == 3)
     +CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,963)
      ELSE ! In this way we will save some space
      CALL RESAMG(RES(IG1),TEMPP(IG1),RESP(IH1),APP(IH1),APN(IH1),
     + APS(IH1),APE(IH1),APW(IH1),APT(IH1),APD(IH1),IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,ICONV3,NGL,CONV2)
      IF(NGL == 3) THEN
      write(963,*) 'residuaali tasolla ',m-1,'mult=',mult
      CALL IPRINT1(RES(IG1),IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,963)
      endif
      ENDIF ! M == 2

      CALL ZERO(RESP(IH2),KOP)
      CALL RESTRI(RESP(IH2),RES(IG1),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),
     + IN,JN,KN,MULT,mcycle,APP(IH2))
       IF(NGL == 3) THEN
       write(151+kmax(1),*) 'restriktioitu residuaali tasolla ',m,
     + 'mult=',mult
      call iprint1(resp(ih2),imax(m),jmax(m),kmax(m),in,jn,kn,
     + 151+kmax(1))
       ENDIF
      CALL ZERO(TEMPP(IG2),KOP)
      CALL LGSAMG(TEMPP(IG2),APP(IH2),APE(IH2),APW(IH2),APS(IH2),
     + APN(IH2),APT(IH2),APD(IH2),RESP(IH2),AA,BB,CC,DD,IMAX(M),
     + JMAX(M),KMAX(M),IN,JN,KN,ITERH,0,2,indd)
c       write(41+kmax(1),*) 'korjaukset tasolla ',m
c      call iprint1(tempp(ig2),imax(m),jmax(m),kmax(m),in,jn,kn,
c     + 41+kmax(1))
5000  CONTINUE

C ... V-cycle up

      DO 6000 M = MGRID,2,-1
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      CALL PROLON(TEMPP(IG1),TEMPP(IG2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN)
c ... Jälkitasoitus ei yleensä kannata, voi aktivoida halutessa
      IF(IPOST > 0) THEN
      IF(M >= 3) THEN
      CALL LGSAMG(TEMPP(IG1),APP(IH1),APE(IH1),APW(IH1),APS(IH1),
     + APN(IH1),APT(IH1),APD(IH1),RESP(IH1),AA,BB,CC,DD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,ITERH,0,2,indd)
c      ELSE ! M = 2 (Smoothing on the first level prohibited)
c      CALL LGSAMG(TEMPP(IG1),AP(IG1),AE(IG1),AW(IG1),AS(IG1),AN(IG1),
c     + APT(IG1),APD(IG1),Q(IG1),AA,BB,CC,DD,IMAX(M-1),JMAX(M-1),
c     +KMAX(M-1),IN,JN,KN,ITERH,0,2)
      ENDIF ! M <= 3
      ENDIF ! IPOST > 0
6000  CONTINUE
      ENDIF ! MGRID > 1 (End of multigrid cycling)

c         CONV3 = SQRT(CONV3) / ((IMAX(1)+2)*(JMAX(1)+2))
c         WRITE(189,*) CONV3
c         WRITE(6,*) MULT,CONV3
7000  CONTINUE ! End of iteration loop
 
      DO 7500 N = 1,ISO ! Put the solution into the first level array
      TEMP(N)= TEMPP(N)
7500  CONTINUE

c      RTIME   = 100.
c      CALL SECAMG(RTIME)
c      RTIME   = RTIME - RTIME1
c      RTIME1  = RTIME/MCYCLE*1000.
c      WRITE(6,*)
c      WRITE(6,9100) RTIME,RTIME1
c9100  FORMAT('KOKONAISAIKA = ',E10.5,'  KIERROS =',F7.2,' ms')
      DEALLOCATE (APP,APN,APS,APE,APW,APT,APD,RES,RESP,TEMPP,DTEMPP,
     + AA,BB,CC,DD)
c     CALL EXIT ! Joku vanha lopetuspaikka
      END SUBROUTINE SOLAMG_p
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE AASMA(AP,AN,AS,AE,AW,AT,AD,APP,APN,APS,APE,APW,APT,
     + APD,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN,ISO,IPI,IFILE,MCYCLE,iterm,M)

      IMPLICIT NONE

      INTEGER IMAX2,JMAX2,KMAX2,IN,JN,KN,ISO,IPI,IFILE,MCYCLE,
     + iterm,IMAX,JMAX,KMAX,ISTR,JSTR,IMAX22,JMAX22,IMJM,IMCJMC,KIND,
     + IMJM2,I,J,K,KKC,KKF,IIC,IIF,IP,IP1,IPN,IPN1,IPK,IP1K,IPN1K,IPNK,
     + I2,M,IP1N

      REAL AP(ISO),AN(ISO),AS(ISO),AW(ISO),AE(ISO),AT(ISO),AD(ISO),
     + APP(IPI),APN(IPI),APS(IPI),APE(IPI),APW(IPI),APT(IPI),APD(IPI)
      REAL SUMA,DIAGCC


      APP     = 0.
      APN     = 0.
      APS     = 0.
      APE     = 0.
      APW     = 0.
      APT     = 0.
      APD     = 0.

      IF(M == 2) DIAGCC = 1.001
      IF(M == 3) DIAGCC = 1.01
      IF(M >= 4) DIAGCC = 1.1


      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2

      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTR*JSTR
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IMJM2   = IMJM
      IF(KMAX == 1) THEN
        KIND = 1
        IMJM2= 0
      ENDIF

C ... Internal points (Last rows for odd numbers)

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5000 I = 1,IMAX2
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IPN1    = IPN  + 1
      IPK     = IP   + IMJM2                       
      IP1K    = IP1  + IMJM2                       
      IPNK    = IPN  + IMJM2                       
      IPN1K   = IPN1 + IMJM2                      
      I2      = IIC  + I

      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IPN1) +
     +          AP(IPK)+ AP(IP1K)+ AP(IPNK)+ AP(IPN1K)+
     +          AN(IP) + AN(IP1) + AN(IPK) + AN(IP1K) +
     +          AS(IPN)+ AS(IPN1)+ AS(IPNK)+ AS(IPN1K)+
     +          AE(IP) + AE(IPN) + AE(IPK) + AE(IPNK) +
     +          AW(IP1)+ AW(IPN1)+ AW(IP1K)+ AW(IPN1K)+
     +          AT(IP) + AT(IP1) + AT(IPN) + AT(IPN1) +
     +          AD(IPK)+ AD(IP1K)+ AD(IPNK)+ AD(IPN1K)
      APW(I2) = AW(IP) + AW(IPN) + AW(IPK) + AW(IPNK)
      APE(I2) = AE(IP1)+ AE(IPN1)+ AE(IP1K)+ AE(IPN1K)
      APN(I2) = AN(IPN)+ AN(IPN1)+ AN(IPNK)+ AN(IPN1K)
      APS(I2) = AS(IP) + AS(IP1) + AS(IPK) + AS(IP1K)
      APT(I2) = AT(IPK)+ AT(IP1K)+ AT(IPNK)+ AT(IPN1K)
      APD(I2) = AD(IP) + AD(IP1) + AD(IPN) + AD(IPN1)
      SUMA    =-APW(I2)-APE(I2)-APN(I2)-APS(I2)-APT(I2)-APD(I2)
c     if(app(i2) <=  1.0001*suma) write(443,*) i,j,k,app(i2),suma
      APP(I2) = MAX(APP(I2),SUMA) ! Should be useless
      APP(I2) = DIAGCC*APP(I2)
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IP1) + AP(IPK) + AP(IP1K)
      APN(I2) = AN(IP) + AN(IP1) + AN(IPK) + AN(IP1K)
      APP(I2) = MAX(-APN(I2)*1.00001,APP(I2))
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IP1) + AP(IPK) + AP(IP1K)
      APS(I2) = AS(IP) + AS(IP1) + AS(IPK) + AS(IP1K)
      APP(I2) = MAX(-APS(I2)*1.00001,APP(I2))
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IPN) + AP(IPK) + AP(IPNK)
      APE(I2) = AE(IP) + AE(IPN) + AE(IPK) + AE(IPNK)
      APP(I2) = MAX(-APE(I2)*1.00001,APP(I2))
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IPN) + AP(IPK) + AP(IPNK)
      APW(I2) = AW(IP) + AW(IPN) + AW(IPK) + AW(IPNK)
      APP(I2) = MAX(-APW(I2)*1.00001,APP(I2))
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM + IMJM
      IF(KIND == 1) KKF = (KN+KIND*K-KIND)*IMJM  !  Eturivi oli + IMJM. Antakee vanha mersu!
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IP1N)
      APT(I2) = AT(IP) + AT(IP1) + AT(IPN) + AT(IP1N)
      APP(I2) = MAX(-APT(I2)*1.00001,APP(I2))
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  ! Takarivi  
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IP1N)
      APD(I2) = AD(IP) + AD(IP1) + AD(IPN) + AD(IP1N)
      APP(I2) = MAX(-APD(I2)*1.00001,APP(I2))
c     IF(APP(I2) < APD(I2)) write(444,*)'apd',app(i2),apd(i2)
5600  CONTINUE
6000  CONTINUE
C ... Corners
      DO 6100 K = 1,KMAX2
      I       = 0
      J       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2+IMAX22))
6100  CONTINUE
      DO 6200 K = 1,KMAX2
      I       = IMAX2 + 1
      J       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2+IMAX22))
6200  CONTINUE
      DO 6300 K = 1,KMAX2
      I       = 0
      J       = JMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2-IMAX22))
6300  CONTINUE
      DO 6400 K = 1,KMAX2
      I       = IMAX2 + 1
      J       = JMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2-IMAX22))
6400  CONTINUE
      DO 7100 J = 1,JMAX2
      I       = 0
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2+IMCJMC))
7100  CONTINUE
      DO 7200 J = 1,JMAX2
      I       = IMAX2 + 1
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2+IMCJMC))
7200  CONTINUE
      DO 7300 J = 1,JMAX2
      I       = 0
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2-IMCJMC))
7300  CONTINUE 
      DO 7400 J = 1,JMAX2
      I       = IMAX2 + 1
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2-IMCJMC))
7400  CONTINUE
      DO 8100 I = 0,IMAX2+1
      J       = 0
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+IMCJMC) + APP(I2+IMAX22))
8100  CONTINUE
      DO 8200 I = 0,IMAX2+1
      J       = JMAX2 + 1
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+IMCJMC) + APP(I2-IMAX22))
8200  CONTINUE
      DO 8300 I = 0,IMAX2+1
      J       = 0
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-IMCJMC) + APP(I2+IMAX22))
8300  CONTINUE 
      DO 8400 I = 0,IMAX2+1
      J       = JMAX2 + 1
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-IMCJMC) + APP(I2-IMAX22))
8400  CONTINUE

C ... Put them together and print

4400  CONTINUE
9120  CONTINUE

c      IF(MCYCLE == 1) THEN ! This can be activated
c      CALL IPRINT1(APP,IMAX2,JMAX2,KMAX2,IN,JN,KN,931)   
c      CALL IPRINT1(APN,IMAX2,JMAX2,KMAX2,IN,JN,KN,932)   
c      CALL IPRINT1(APS,IMAX2,JMAX2,KMAX2,IN,JN,KN,933)   
c      CALL IPRINT1(APE,IMAX2,JMAX2,KMAX2,IN,JN,KN,934)   
c      CALL IPRINT1(APW,IMAX2,JMAX2,KMAX2,IN,JN,KN,935)   
c      CALL IPRINT1(APT,IMAX2,JMAX2,KMAX2,IN,JN,KN,936)   
c      CALL IPRINT1(APD,IMAX2,JMAX2,KMAX2,IN,JN,KN,937)   
c      ENDIF ! MCYCLE == 1
9000  FORMAT(/16F7.3)
      END SUBROUTINE AASMA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESTRI(RESP,RES,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN,
     + MCYCLE,ITERM,APP)
      REAL RESP(*),RES(*),APP(*)

C ... Restrict the residual (Last rows for odd numbers)

      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2

      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTR*JSTR
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IMJM2   = IMJM
      IF(KMAX == 1) THEN
        KIND = 1
        IMJM2= 0
      ENDIF

C ... Internal points

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5000 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      IPK     = IP    + IMJM2                       
      IP1K    = IP1   + IMJM2                       
      IPNK    = IPN   + IMJM2                       
      IPN1K   = IPN1  + IMJM2                      
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IPN1) + 
     +          RES(IPK)+ RES(IP1K)+ RES(IPNK)+ RES(IPN1K)
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPK) + RES(IP1K)
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPK) + RES(IP1K)
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IPN) + RES(IPK) + RES(IPNK)
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IPN) + RES(IPK) + RES(IPNK)
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM + (KIND-1)*IMJM !+ IMJM ! Eturivi  
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IP1N)
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IP1N)
5600  CONTINUE
C ... Corners are put to zero in multi
      END SUBROUTINE RESTRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PROLON(TEMP,TEMP2,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN)
      REAL TEMP(*),TEMP2(*)

C ... Prolongation

      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2
      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMJM    = ISTR*JSTR

      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IF(KMAX == 1) KIND = 1

C ... Internal points (Last row for odd numbers)

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5010 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
5010  CONTINUE
      IF(KMAX > 1) THEN
      DO 5020 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      IPK     = IP    + IMJM                       
      IP1K    = IP1   + IMJM                       
      IPNK    = IPN   + IMJM                       
      IPN1K   = IPN1  + IMJM                      
      I2      = IIC + I
      TEMP(IPK)  = TEMP(IPK)   + TEMP2(I2)
      TEMP(IP1K) = TEMP(IP1K)  + TEMP2(I2)
      TEMP(IPNK) = TEMP(IPNK)  + TEMP2(I2)
      TEMP(IPN1K)= TEMP(IPN1K) + TEMP2(I2)
5020  CONTINUE
      ENDIF ! KMAX2 > 1
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM
      IP1K    = IP1+ IMJM
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IP1K)= TEMP(IP1K) + TEMP2(I2)
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM
      IP1K    = IP1+ IMJM
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IP1K)= TEMP(IP1K) + TEMP2(I2)
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM
      IPNK    = IPN+ IMJM
      I2      = IIC + I
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IPNK)= TEMP(IPNK) + TEMP2(I2)
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM
      IPNK    = IPN+ IMJM
      I2      = IIC + I
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IPNK)= TEMP(IPNK) + TEMP2(I2)
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  + (KIND-1)*IMJM !+ IMJM ! Eturivi 
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      TEMP(IP) = TEMP(IP)  + TEMP2(I2)
      TEMP(IP1)= TEMP(IP1) + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IP1N)= TEMP(IP1N) + TEMP2(I2)
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      TEMP(IP) = TEMP(IP)  + TEMP2(I2)
      TEMP(IP1)= TEMP(IP1) + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IP1N)= TEMP(IP1N) + TEMP2(I2)
5600  CONTINUE
      END SUBROUTINE PROLON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SUBAMG(TEMP,DTEMP,KOK)

      REAL :: TEMP(KOK), DTEMP(KOK)

      TEMP  = TEMP + DTEMP

      RETURN
      END SUBROUTINE SUBAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE ZERO(TEMP,KOK)

      REAL :: TEMP(KOK)

      TEMP = 0.

      RETURN
      END SUBROUTINE ZERO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      

C
C ... 3D routines for the solution of linear equations

      SUBROUTINE SYAMG(IL,IU,BB,DD,AA,CC)

      REAL :: AA(*),BB(*),CC(*),DD(*)
C ...
C ... SUBROUTINE SY SOLVES TRIDIAGONAL SYSTEM BY ELIMINATION
C ... IL = SUBSCRIPT OF FIRST EQUATION
C ... IU = SUBSCRIPT OF LAST EQUATION
C ... BB = COEFFICIENT BEHIND DIAGONAL
C ... DD = COEFFICIENT ON DIAGONAL
C ... COEFFICIENT AHEAD OF DIAGONAL
C ... ELEMENT OF CONSTANT VECTOR
C ...
C ... ESTABLISH UPPER TRIANGULAR MATRIX
C ...
      LP = IL+1
      DO 10 I = LP,IU
      R = BB(I)/DD(I-1)
      DD(I) = DD(I)-R*AA(I-1)
   10 CC(I) = CC(I)-R*CC(I-1)
C ...
C ... BACK SUBSTITUTION
C ...
      CC(IU) = CC(IU)/DD(IU)
      DO 20 I = LP,IU
      J = IU-I+IL
   20 CC(J) = (CC(J)-AA(J)*CC(J+1))/DD(J)
C ...
C ... SOLUTION STORED IN CC
C ...
      RETURN
      END SUBROUTINE SYAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE LGSAMG(DPCORR,AP,AE,AW,AS,AN,AT,AD,RMASS,AA,BB,CC,DD,
     +IMAX,JMAX,KMAX,IN,JN,KN,ITERM,ICON,KIERR,INDD)
      REAL DPCORR(*),AP(*),AE(*),AW(*),AS(*),AN(*),AT(*),AD(*),RMASS(*)
     + ,DD(*),AA(*),BB(*),CC(*)
      REAL, ALLOCATABLE :: ED(:) ! Voi aktivoida

C ... LINE GAUSS-SEIDEL ITERATION

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR
      NTOT   = ISTR*JSTR*(KMAX+2*KN)
      ALLOCATE (ED(NTOT))
c      DO 8000 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(71,9001) (DPCORR(I+II),I=1,IMAX+2)
9001  FORMAT(/16F7.3)
c8000  CONTINUE

      ED     = 0.

      DO 9000 ITER = 1,ITERM
      DPCORM = 0.
      REL    = 0.0

      IF(INDD == 13) go to 3999
      DO 2000 K = 1,KMAX!+1 ! 2D
      JJ      = (KN+K-1)*IL
      DO 2000 J = 1,JMAX!+1 ! Sweeppaa jk-suunnassa
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 0,IMAX+1 ! Aidosti yli laitojen 
      M       = I + 1
      IJ      = I + II
      IF(J <= JMAX) THEN
         DPCP = DPCORR(IJ+ISTR)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(J > 0) THEN
         DPCM = DPCORR(IJ-ISTR)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(K <= KMAX) THEN
         DPCT = DPCORR(IJ+IL)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(K > 0) THEN
         DPCD = DPCORR(IJ-IL)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(M)   = -RMASS(IJ)-AN(IJ)*DPCP-AS(IJ)*DPCM
     +                    -AT(IJ)*DPCT-AD(IJ)*DPCD
      AA(M)   = AE(IJ)
      DD(M)   = AP(IJ)  
      BB(M)   = AW(IJ)
1000  CONTINUE
c      CALL SYAMG(2,IMAX+3,BB,DD,AA,CC)
      CALL SYAMG(1,IMAX+2,BB,DD,AA,CC)
      DO 1500 I = 0,IMAX+1 ! Aidosti yli laitojen 
      M       = I + 1
      IJ      = I + II
      DPCORR(IJ) = CC(M)
1500  CONTINUE
2000  CONTINUE
c      DO 8001 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(72,*) ((DPCORR(I+II),i,j),I=1,IMAX+2)
9002  FORMAT(/16F7.3,2I2)
c8001  CONTINUE

3999  CONTINUE
      IF(INDD == 13) go to 5999
      DO 4000 K = 1,KMAX!+1
      JJ      = (KN+K-1)*IL
      DO 4000 I = 1,IMAX!+1 ! Toinen kierros sweeppaa ik-suunnassa
      M       = I + 1 
      DO 3000 J = 0,JMAX+1 ! Aidosti yli laitojen 
      N       = J + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      IF(I <= IMAX) THEN
         DPCP = DPCORR(IJ+1)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(I > 0) THEN
         DPCM = DPCORR(IJ-1)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(K <= KMAX) THEN
         DPCT = DPCORR(IJ+IL)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(K > 0) THEN
         DPCD = DPCORR(IJ-IL)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(N)   = -RMASS(IJ)-AE(IJ)*DPCP-AW(IJ)*DPCM
     +                    -AT(IJ)*DPCT-AD(IJ)*DPCD
      AA(N)   = AN(IJ)
      DD(N)   = AP(IJ)  
      BB(N)   = AS(IJ)
c      write(888,*) n,aa(n),dd(n),BB(N),cc(n)
3000  CONTINUE
c      write(888,*) '------------'
      CALL SYAMG(1,JMAX+2,BB,DD,AA,CC)
      DO 3500 J = 0,JMAX+1 ! Aidosti yli laitojen 
      N       = J + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      DPCORR(IJ) = CC(N)
cc      write(888,*) n,cc(n)
3500  CONTINUE
4000  CONTINUE
5999  CONTINUE
      
      IF(KIERR > 1) THEN
      DO 6000 J = 1,JMAX!+1 
      DO 6000 I = 1,IMAX!+1 ! Kolmas kierros sweeppaa ij-suunnassa
      DO 5000 K = 0,KMAX+1 ! Aidosti yli laitojen
      L       = K + 1
      JJ      = (KN+K-1)*IL
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      IF(I <= IMAX) THEN
         DPCP = DPCORR(IJ+1)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(I > 0) THEN
         DPCM = DPCORR(IJ-1)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(J <= JMAX) THEN
         DPCT = DPCORR(IJ+ISTR)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(J > 0) THEN
         DPCD = DPCORR(IJ-ISTR)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(L)   = -RMASS(IJ)-AN(IJ)*DPCT-AS(IJ)*DPCD
     +                    -AE(IJ)*DPCP-AW(IJ)*DPCM
      AA(L)   = AT(IJ)
      DD(L)   = AP(IJ)  
      BB(L)   = AD(IJ)

5000  CONTINUE
      CALL SYAMG(1,KMAX+2,BB,DD,AA,CC)
      DO 5500 K = 0,KMAX+1 ! Aidosti yli laitojen 
      L       = K + 1
      JJ      = (KN+K-1)*IL
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      DPCORR(IJ) = CC(L)
5500  CONTINUE
6000  CONTINUE
      ENDIF ! KMAX > 1


c      DO 4005 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(1073,*) (DPCORR(I+II),I=1,IMAX+2)
c4005  CONTINUE 
c      DO 7000 J = 0,JMAX+1  ! Aidosti yli laitojen 
c      N       = J + 1 
c      II      = (JN+J-1)*ISTR + IN
c      DO 7000 I = 0,IMAX+1
c      M       = I + 1
cc      IJ      = I + II
c        DPCORM    = DPCORM + ABS(ED(IJ)-DPCORR(IJ))
c         ED(IJ)    = DPCORR(IJ)
c7000  CONTINUE
c      DPCORM  = DPCORM / ((IMAX+2)*(JMAX+2))
c      write(1088,*) dpcorm!Convergence monitor can be activated

9000  CONTINUE ! END OF ITERATION
c         write(1088,*) '-------------------------'
      END SUBROUTINE LGSAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE LGSRSM(DPCORR,D,BM,BN,BW,AM,AN,AW,RMASS,IMAX,JMAX,KMAX,
     + IN,JN,KN,ISTRID,JSTRID,ITERM)

      REAL DPCORR(*),D(*),AM(*),AN(*),AW(*),BM(*),BN(*),BW(*),RMASS(*),
     +     AA(300),BB(300),CC(300),DD(300),ED(300*300)

C ... RMASS   = RIGHT hand side
C ... DPCORR  = RESULT vector
C ... D       = DIAGONAL vector
C ... AM,BM   = i-direction
C ... AN,BN   = j-direction
C ... AW,BW   = k-direction

      IL      = ISTRID*JSTRID
      ISTR   = 1
      JSTR   = IMAX + 2*IN
      KSTR   = IL
      IF(KMAX == 1) KSTR  = 0

      DO K = 0,KMAX+1
         KA      = (KN+K-1)*IL
         DO J = 0,JMAX+1
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO I = 0,IMAX+1
               N     = JJ + I
c               DPCORR(N) = 0.
               ED(N) = DPCORR(N)      
            ENDDO
         ENDDO
      ENDDO

C ... LINE GAUSS-SEIDEL ITERATION

      IL      = ISTRID*JSTRID
      RELAX = 1.4              !over relaxation
      EPS   = 1.E-8
      DO 3000 ITER = 1,ITERM
      dpcorm = 0.
C ... i-direction
      DO K = 1,KMAX
         KA      = (KN+K-1)*IL
         DO J = 1,JMAX
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO I = 0,IMAX+1
               N     = JJ + I
               M     = I  + 2
               CC(M) = RMASS(N) - 
     +              BN(N)*DPCORR(N-JSTR) - AN(N)*DPCORR(N+JSTR) -
     +              BW(N)*DPCORR(N-KSTR) - AW(N)*DPCORR(N+KSTR)
               AA(M) = AM(N)
               DD(M) = D(N)
               BB(M) = BM(N)
            ENDDO
            CALL SY(2,IMAX+3,BB,DD,AA,CC)
            DO I = 0,IMAX+1
               N     = JJ + I
               M     = I  + 2
               DPCORR(N) = DPCORR(N) + RELAX*(CC(M) - DPCORR(N))
            ENDDO
         ENDDO
      ENDDO
 112  format(2I4,8E12.4)

C ... j-direction
      DO K = 1,KMAX
         KA      = (KN+K-1)*IL
         DO I = 1,IMAX
            DO J = 0,JMAX+1
               JJ    = (JN+J-1)*ISTRID + IN + KA
               N     = JJ + I
               M     = J  + 2
               CC(M) = RMASS(N) - 
     +              BM(N)*DPCORR(N-ISTR) - AM(N)*DPCORR(N+ISTR) -
     +              BW(N)*DPCORR(N-KSTR) - AW(N)*DPCORR(N+KSTR)
               AA(M) = AN(N)
               DD(M) = D(N)
               BB(M) = BN(N)
c               if(i == 40 .and. iter == 100) 
c     +              write(*,111) j,BB(M),DD(M),AA(M),CC(M)
 111  format(1I4,8E12.4)
            ENDDO
            CALL SY(2,JMAX+3,BB,DD,AA,CC)
            DO J = 0,JMAX+1
               JJ    = (JN+J-1)*ISTRID + IN + KA
               N     = JJ + I
               M     = J  + 2
               DPCORR(N) = DPCORR(N) + RELAX*(CC(M) - DPCORR(N))
            ENDDO
         ENDDO
      ENDDO

      IF(KMAX /= 1) THEN
C ... k-direction
      DO J = 1,JMAX
         DO I = 1,IMAX
            DO K = 0,KMAX+1
               KA      = (KN+K-1)*IL
               JJ      = (JN+J-1)*ISTRID + IN + KA
               N     = JJ + I
               M     = K  + 2
               CC(M) = RMASS(N) - 
     +              BN(N)*DPCORR(N-JSTR) - AN(N)*DPCORR(N+JSTR) -
     +              BM(N)*DPCORR(N-ISTR) - AM(N)*DPCORR(N+ISTR)
               AA(M) = AW(N)
               DD(M) = D(N)
               BB(M) = BW(N)
            ENDDO
            CALL SY(2,KMAX+3,BB,DD,AA,CC)
            DO K = 0,KMAX+1
               KA    = (KN+K-1)*IL
               JJ    = (JN+J-1)*ISTRID + IN + KA
               N     = JJ + I
               M     = K  + 2
               DPCORR(N) = DPCORR(N) + RELAX*(CC(M) - DPCORR(N))
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      DO K = 1,KMAX
         KA      = (KN+K-1)*IL
         DO J = 1,JMAX
            JJ      = (JN+J-1)*ISTRID + IN + KA
            DO I = 1,IMAX
               N     = JJ + I
               DPCORM = DPCORM + (ED(N)-DPCORR(N))**2/(ED(N)**2+EPS)
               ED(N)  = DPCORR(N)
            ENDDO
         ENDDO
      ENDDO
      DPCORM = SQRT(dpcorm/(IMAX*JMAX*KMAX))
c      if(imax >= 40) write(*,*) iter,dpcorm
c      stop
      IF(DPCORM <= 1.e-6) GOTO 3003

3000  CONTINUE ! END OF ITERATION
 3003 CONTINUE
c     if(imax >= 40 .and. iter /= 1) write(88,*) imax,iter,DPCORM
      RETURN
      END SUBROUTINE LGSRSM
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE LGS(DPCORR,AP,AE,AW,AS,AN,RMASS,IMAX,JMAX,MAX,ITERM)

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: DPCORR(MAX,*), AP(MAX,*), AE(MAX,*), AW(MAX,*), AS(MAX,*),
     + AN(MAX,*), RMASS(MAX,*), ED(300,300), AA(300), BB(300), CC(300),
     + DD(300)


C ... LINE GAUSS-SEIDEL ITERATION

      DO J = 1,JMAX+4
         DO I = 1,IMAX+4
c            DPCORR(I,J) = 0.
            ED(I,J)     = DPCORR(I,J)
         ENDDO
      ENDDO

      DO 9000 ITER = 1,ITERM
      DPCORM = 0.
      DPTOT  = 0.
      REL    = 0.0

      DO 2000 N = 1,JMAX
      J       = N + 2 
      DO 1000 M = 0,IMAX+1 ! Aidosti yli laitojen 
      I       = M + 2
      CC(I)   = -RMASS(I,J)-AN(I,J)*DPCORR(I,J+1)-AS(I,J)*DPCORR(I,J-1)
      AA(I)   = AE(I,J)
      DD(I)   = AP(I,J)
      BB(I)   = AW(I,J)
1000  CONTINUE
      CALL SY(2,IMAX+3,BB,DD,AA,CC)
      DO 1500 M = 0,IMAX+1 ! Aidosti yli laitojen 
      I       = M + 2
      DPCORR(I,J) = CC(I)
1500  CONTINUE
2000  CONTINUE

      DO 4000 M = 1,IMAX
      I       = M + 2 
      DO 3000 N = 0,JMAX+1 ! Aidosti yli laitojen 
      J       = N + 2
      CC(J)   = -RMASS(I,J)-AE(I,J)*DPCORR(I+1,J)-AW(I,J)*DPCORR(I-1,J)
      AA(J)   = AN(I,J)
      DD(J)   = AP(I,J)
      BB(J)   = AS(I,J)
3000  CONTINUE
      CALL SY(2,JMAX+3,BB,DD,AA,CC)
      DO 3500 N = 0,JMAX+1 ! Aidosti yli laitojen 
      J       = N + 2
      DPCORR(I,J) = CC(J)
3500  CONTINUE
4000  CONTINUE

      DO 5000 N = 1,JMAX 
      J       = N + 2 
      DO 5000 M = 1,IMAX
      I       = M + 2
         DPCORM     = DPCORM + ABS(ED(I,J)-DPCORR(I,J))
         DPTOT      = DPTOT  + ABS(DPCORR(I,J))
         ED(I,J)    = DPCORR(I,J)
5000  CONTINUE

9000  CONTINUE ! END OF ITERATION


 333  CONTINUE

C     write(88,*) iter,iterm,dpcorm,dpcorm/(DPTOT+1.E-5) ! Convergence monitor
      CLOSE(89)
      RETURN
      END SUBROUTINE LGS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SY(IL,IU,BB,DD,AA,CC)

      REAL :: AA(*),BB(*),CC(*),DD(*)
C ...
C ... SUBROUTINE SY SOLVES TRIDIAGONAL SYSTEM BY ELIMINATION
C ... IL = SUBSCRIPT OF FIRST EQUATION
C ... IU = SUBSCRIPT OF LAST EQUATION
C ... BB = COEFFICIENT BEHIND DIAGONAL
C ... DD = COEFFICIENT ON DIAGONAL
C ... AA = COEFFICIENT AHEAD OF DIAGONAL
C ... CC = ELEMENT OF CONSTANT VECTOR
C ... OR
C ... DD(1)  AA(1)  0     0 ...  0
C ... BB(2)  DD(2)  AA(2) 0 ...  0
C ...   .     .      .    .      .
C ...   .     .      .    .      .
C ...  0      0     ...   BB(N) DD(N)     

C ...
C ... ESTABLISH UPPER TRIANGULAR MATRIX
C ...
      LP = IL+1
      DO 10 I = LP,IU
      R = BB(I)/DD(I-1)
      DD(I) = DD(I)-R*AA(I-1)
   10 CC(I) = CC(I)-R*CC(I-1)
C ...
C ... BACK SUBSTITUTION
C ...
      CC(IU) = CC(IU)/DD(IU)
      DO 20 I = LP,IU
      J = IU-I+IL
   20 CC(J) = (CC(J)-AA(J)*CC(J+1))/DD(J)
C ...
C ... SOLUTION STORED IN CC
C ...
      RETURN
      END SUBROUTINE SY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SECAMG(RTIME)
C
C ... Monitor CPU-time on the IRIS
C
      REAL TARRAY(2)
      APU = RTIME
*      RTIME = DTIME(TARRAY)
      RTIME = DPTIME(TARRAY)
      IF(APU < 1.E-10) RTIME = 0.
      END SUBROUTINE SECAMG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PRECOR(NBL,M,RO,RM,RN,RW,E,P,PDIFF,U,V,W,C,DU,DV,DW,
     1 DP,DRODP,A1,A1XA,A1YA,A1ZA,A2,A2XA,A2YA,A2ZA,A3,A3XA,A3YA,A3ZA,
     2 VOL,VIS,VIST,DTL,ABU,ABV,ABW,DPCORR,RMASS,A1PP,A1RM,A1RN,A1RW,
     3 ICP,JCP,KCP,IMAX,JMAX,KMAX,ICYCLE,INTERI,INTERJ,INTERK,JBOT,
     4 JTOP,IDI1,IDI2,IDI3,IT,IL,IK,IPRINT,ICON,NPATCH,IBOT,ITOP,KBOT,
     5 KTOP,IFLUX,UROT,VROT,WROT,ITURB,FRSDEN,FRSPRE,ITERM,TOLER,
     6 ZZZ,MAXW)

      USE CHARACTERS
      USE NS3CO, ONLY : IN, JN, KN, IC9
      
cc      INCLUDE 'NS3CH.C'

      DIMENSION RO(*),RM(*),RN(*),RW(*),E(*),P(*),U(*),V(*),W(*),C(*),
     2 A1(*),A1XA(*),A1YA(*),A1ZA(*),A2(*),A2XA(*),A2YA(*),A2ZA(*),
     3 A3(*),A3XA(*),A3YA(*),A3ZA(*),VOL(*),PDIFF(*),DPCORR(*),
     4 DRODP(*),ABU(*),ABV(*),ABW(*),A1PP(*),A1RM(*),A1RN(*),A1RW(*),
     5 ICP(*),JCP(*),KCP(*),UROT(*),VROT(*),WROT(*),ICON(IC9,*),
     6 DU(*),DV(*),DW(*),DP(*),VIS(*),VIST(*),DTL(*),RMASS(*)
       REAL AA(100),BB(100),CCC(100),DD(100),ZZZ(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN

      NTOT    = ISTRID*JSTRID*KSTRID
c       write(6,*) imax,jmax,kmax,icycle,interi,interj,interk,idi1,idi2,
c     + idi3,it,il,ik,iprint,npatch,iflux
C ... ADJUST CELLS WITH FRICTION TO THE GRID LEVEL USED

      IDM1    = (IDI1-1)/2**(M-1) + 1
      IDM2    = (IDI2-1)/2**(M-1) + 1
      IDM3    = (IDI3-1)/2**(M-1) + 1

C **********************************************************
C                                                          *
C ... CALCULATION OF PRESSURE CORRECTION                   *
C                                                          *
C **********************************************************

      WRITE(*,*) 'This pressure correction method does not work'
      WRITE(*,*) ' any more. Exiting...'
      STOP 'PRECOR'
      END SUBROUTINE PRECOR

C *******************************************************************
C     ROUTINES FOR THE SOLUTION OF LINEAR EQUATION SETS IN d_p      *
C *******************************************************************

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SOLAMG_DP(TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAXX,JMAXX,
     + KMAXX,IN,JN,KN,MCYCLE,ITERM,ITERH,MGRID,ICONV3,ipostt,NGL,ICYCLE)

C ... Multigrid solver

      IMPLICIT NONE

      INTEGER :: IPOSTT,IPOST,ICONV3,MGRID,ITERH,ITERM,MCYCLE,IN,JN,KN,
     + IMAXX,JMAXX,KMAXX,M,INDD,ISO,LINE,KOKO,IPI,IPL,ISTR,JSTR,IL,I,
     + MULT,IG1,IH1,IH2,KOK,KOP,IG2,N,NGL,ICYCLE

      INTEGER :: IMAX(10),JMAX(10),KMAX(10),IG(10),IH(10)

      REAL :: AP(*),AN(*),AS(*),AW(*),AE(*),AT(*),AD(*),TEMP(*),Q(*)
      REAL, ALLOCATABLE::APP(:),APN(:),APS(:),APE(:),APW(:),
     +      APT(:),APD(:),RESP(:),DTEMPP(:),RES(:),TEMPP(:),
     +      AA(:),BB(:),CC(:),DD(:)

      indd = ipostt
      IPOST   = 0 ! No post smoothing in this version (IPOST = 0)

C ... Problem dimensions and array pointers

      IMAX    = 0
      JMAX    = 0
      KMAX    = 0
      IG      = 0
      IH      = 0
      IH(1)   = 1
      IH(2)   = 1
      IG(1)   = 1
      ISO     = (IMAXX + 2*IN)*(JMAXX + 2*JN)*(KMAXX + 2*KN)
      IMAX(1) = IMAXX
      JMAX(1) = JMAXX
      KMAX(1) = KMAXX
      LINE    = AMAX0(IMAXX+2*IN,JMAXX+2*JN,KMAXX+2*KN) 

      DO 1000 M = 2,MGRID

         IMAX(M) = (IMAX(M-1)+1)/2
         JMAX(M) = (JMAX(M-1)+1)/2
         KMAX(M) = (KMAX(M-1)+1)/2
         KMAX(M) = MAX(1,KMAX(M)) ! Minimum number of slabs
         KOKO    = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
         IG(M)   = IG(M-1) + KOKO + 10
         IF(M > 2) IH(M) = IH(M-1) + KOKO + 10

1000  CONTINUE
      
      KOKO = (IMAX(MGRID)+2*IN)*(JMAX(MGRID)+2*JN)*(KMAX(MGRID)+2*KN)
      IPI  = IG(MGRID) + KOKO
      IPL  = IPI - ISO - 10

C ... Problem size determined. Start to allocate

      ALLOCATE (RES(IPI),TEMPP(IPI),DTEMPP(IPI),
     + APP(IPL),APN(IPL),APS(IPL),APE(IPL),APW(IPL),RESP(IPL),
     + APT(IPL),APD(IPL),AA(LINE),BB(LINE),CC(LINE),DD(LINE))

      ISTR    = IMAX(1) + 2*IN
      JSTR    = JMAX(1) + 2*JN
      IL      = ISTR*JSTR

      DTEMPP  = 0.
      RESP    = 0.
      RES     = 0.

      DO I = 1,ISO ! First level size
      TEMPP(I)= TEMP(I)
      ENDDO

c      CALL SECAMG(RTIME)
c      RTIME1  = RTIME
C
C ... Iteration cycle begins

      DO 7000 MULT = 1,MCYCLE
      IG1     = IG(1) ! Useless, starting addrres = 1
      RES     = 0.
      CALL RESAMG_DP(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ICONV3,mult,11,NGL)

      DTEMPP = 0.
c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,960)
      CALL LGSAMG_DP(DTEMPP,AP,AE,AW,AS,AN,AT,AD,RES,AA,BB,CC,DD,
     + IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,ITERM,0,2,indd)
      CALL SUBAMG_DP(TEMPP,DTEMPP,ISO)
c       write(41+kmax(1),*) 'tihein taso kirros ',mult
c      call iprint1(tempp,imax(1),jmax(1),kmax(1),in,jn,kn,41+kmax(1))

C ... Perform multigrid?
 
      IF(MGRID > 1) THEN ! Start multigrid
      IF(MULT == 1) THEN
C ... Generate the matrix recursively on the coarse levels
      DO 4000 M = 2,MGRID
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOK     = (IMAX(M-1)+2*IN)*(JMAX(M-1)+2*JN)*(KMAX(M-1)+2*KN)
      KOP     = (IMAX(M)  +2*IN)*(JMAX(M)  +2*JN)*(KMAX(M)  +2*KN)
      IF(M == 2) THEN ! Use the first level arrays
      CALL AASMA_DP(AP,AN,AS,AE,AW,AT,AD,APP(IH2),APN(IH2),APS(IH2),
     + APE(IH2),APW(IH2),APT(IH2),APD(IH2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN,KOK,KOP,1,MULT)
      ELSE ! In this way we will save some space 
      CALL AASMA_DP(APP(IH1),APN(IH1),APS(IH1),APE(IH1),APW(IH1),
     + APT(IH1),APD(IH1),APP(IH2),APN(IH2),APS(IH2),APE(IH2),APW(IH2),
     + APT(IH2),APD(IH2),IMAX(M),JMAX(M),KMAX(M),KMAX(M-1),IN,JN,KN,KOK,
     + KOP,1,MULT)
      ENDIF ! M == 2
4000  CONTINUE
      ENDIF ! MULT == 1

C ... V-cycle down

      DO 5000 M = 2,MGRID
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      IH2     = IH(M)
      KOP     = (IMAX(M)+2*IN)*(JMAX(M)+2*JN)*(KMAX(M)+2*KN)

C ... Recalculate the residual on the previous level.
 
      IF(M == 2) THEN ! First level residual
      RES     = 0.
      CALL RESAMG_DP(RES,TEMPP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,12,NGL)
c      CALL IPRINT1(RES,IMAX(1),JMAX(1),KMAX(1),IN,JN,KN,961)
      ELSE ! In this way we will save some space
      CALL RESAMG_DP(RES(IG1),TEMPP(IG1),RESP(IH1),APP(IH1),APN(IH1),
     + APS(IH1),APE(IH1),APW(IH1),APT(IH1),APD(IH1),IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,0,mult,13,NGL)
      ENDIF ! M == 2
c       write(51+kmax(1),*) 'residuaali tasolla ',m-1
c      call iprint1(res(ig1),imax(m-1),jmax(m-1),kmax(m-1),in,jn,kn,
c     + 51+kmax(1))

      CALL ZERO_DP(RESP(IH2),KOP)
      CALL RESTRI_DP(RESP(IH2),RES(IG1),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN)
c       write(51+kmax(1),*) 'RESTRI_DPktioitu residuaali tasolla ',m-1
c      call iprint1(resp(ih2),imax(m),jmax(m),kmax(m),in,jn,kn,
c     + 51+kmax(1))
      CALL ZERO_DP(TEMPP(IG2),KOP)
      CALL LGSAMG_DP(TEMPP(IG2),APP(IH2),APE(IH2),APW(IH2),APS(IH2),
     + APN(IH2),APT(IH2),APD(IH2),RESP(IH2),AA,BB,CC,DD,IMAX(M),
     + JMAX(M),KMAX(M),IN,JN,KN,ITERH,0,2,indd)
c       write(41+kmax(1),*) 'korjaukset tasolla ',m
c      call iprint1(tempp(ig2),imax(m),jmax(m),kmax(m),in,jn,kn,
c     + 41+kmax(1))
5000  CONTINUE

C ... V-cycle up

      DO 6000 M = MGRID,2,-1
      IG1     = IG(M-1)
      IG2     = IG(M)
      IH1     = IH(M-1)
      CALL PROLON_DP(TEMPP(IG1),TEMPP(IG2),IMAX(M),JMAX(M),KMAX(M),
     + KMAX(M-1),IN,JN,KN)
c ... Jälkitasoitus ei yleensä kannata, voi aktivoida halutessa
      IF(IPOST > 0) THEN
      IF(M >= 3) THEN
      CALL LGSAMG_DP(TEMPP(IG1),APP(IH1),APE(IH1),APW(IH1),APS(IH1),
     + APN(IH1),APT(IH1),APD(IH1),RESP(IH1),AA,BB,CC,DD,IMAX(M-1),
     + JMAX(M-1),KMAX(M-1),IN,JN,KN,ITERH,0,2,indd)
c      ELSE ! M = 2 (Smoothing on the first level prohibited)
c      CALL LGSAMG_DP(TEMPP(IG1),AP(IG1),AE(IG1),AW(IG1),AS(IG1),AN(IG1),
c     + APT(IG1),APD(IG1),Q(IG1),AA,BB,CC,DD,IMAX(M-1),JMAX(M-1),
c     +KMAX(M-1),IN,JN,KN,ITERH,0,2)
      ENDIF ! M <= 3
      ENDIF ! IPOST > 0
6000  CONTINUE
      ENDIF ! MGRID > 1 (End of multigrid cycling)

c         CONV3 = SQRT(CONV3) / ((IMAX(1)+2)*(JMAX(1)+2))
c         WRITE(189,*) CONV3
c         WRITE(6,*) MULT,CONV3
7000  CONTINUE ! End of iteration loop
 
      DO 7500 N = 1,ISO ! Put the solution into the first level array
      TEMP(N)= TEMPP(N)
7500  CONTINUE

c      RTIME   = 100.
c      CALL SECAMG(RTIME)
c      RTIME   = RTIME - RTIME1
c      RTIME1  = RTIME/MCYCLE*1000.
c      WRITE(6,*)
c      WRITE(6,9100) RTIME,RTIME1
c9100  FORMAT('KOKONAISAIKA = ',E10.5,'  KIERROS =',F7.2,' ms')
      DEALLOCATE (APP,APN,APS,APE,APW,APT,APD,RES,RESP,TEMPP,DTEMPP,
     + AA,BB,CC,DD)
c     CALL EXIT ! Joku vanha lopetuspaikka
      END SUBROUTINE SOLAMG_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESAMG_DP(RES,TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,
     + KMAX,IN,JN,KN,ICONV3,mult,ICONV,NGL)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ICONV3,MULT,ICONV,NGL,ISTR,
     + JSTR,KIND,I,J,K,II,JJ,IJ,IL

      REAL :: CONV,CONV2

      REAL :: RES(*),TEMP(*),AP(*),AN(*),AS(*),AE(*),AW(*),
     +        AD(*),AT(*),Q(*)

C ... Calculate the residual of the equation to be solved

      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KIND    = 1
c      KIND    = 2
c      IF(KMAX == 1) KIND = 1
      IL      = ISTR*JSTR
C ... Internal points
      DO 1000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = (AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+
     + AS(IJ)*TEMP(IJ-ISTR)+AE(IJ)*TEMP(IJ+1)+AW(IJ)*TEMP(IJ-1)+
     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ))
c      if(i == 10 .and.j == 1 .and. k == 1) 
c      write(808,*) ij,res(ij),temp(ij),q(ij),an(ij),as(ij)
1000  CONTINUE
C ... Sides
      DO 1100 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1100 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AE(IJ)*TEMP(IJ+1)+Q(IJ)
1100  CONTINUE
      DO 1200 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1200 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX + 1
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AW(IJ)*TEMP(IJ-1)+Q(IJ)
1200  CONTINUE
      DO 1300 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1300 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+Q(IJ)
1300  CONTINUE
      DO 1400 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1400 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AS(IJ)*TEMP(IJ-ISTR)+Q(IJ)
1400  CONTINUE
      K       = 0
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1500 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1500 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AT(IJ)*TEMP(IJ+IL)+Q(IJ)
1500  CONTINUE
      K       = KMAX + 1
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1600 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1600 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AD(IJ)*TEMP(IJ-IL)+Q(IJ)
c      if(mult ==  0) write(6,*) res(ij),ij,ap(ij),ad(ij),temp(ij),
c     + temp(ij-il),q(ij)
1600  CONTINUE

C ... Convergence
      CONV    = 0.
      CONV2   = 0.
c      if(iconv == 1 .and. mult == 50) then
c      WRITE(71,*)
c      WRITE(71,*) 'RESIDUAL'
c      CALL IPRINT1(RES,IMAX,JMAX,KMAX,IN,JN,KN,71)
c      ENDIF 
      DO 8000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 8000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 8000 I = 0,IMAX+1
      CONV    = CONV + RES(I+II)**2
      CONV2   = MAX(CONV2,ABS(RES(I+II))) !Maximum value
c      CONV2   = CONV2 + RES(I+II) ! Mielenkiintoinen ilmiö
8000  CONTINUE
9000  FORMAT(/12F8.3)
      IF(CONV > 1.E-30) THEN
      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2)*(KMAX+2))
      ELSE
      CONV    = 0.
      ENDIF
c      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2))
      IF(ICONV3 == 1) WRITE(188,*) NGL,ICONV,CONV,CONV2
c      IF(ICONV == 2) WRITE(189,*) CONV,CONV2
      END SUBROUTINE RESAMG_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESAMG_DPF(RES,TEMP,Q,AP,AN,AS,AE,AW,AT,AD,IMAX,JMAX,
     + KMAX,IN,JN,KN,ICONV3,mult,ICONV,NGL)
      REAL RES(*),TEMP(*),AP(*),AN(*),AS(*),AE(*),AW(*),AD(*),AT(*),
     + Q(*)

C ... Calculate the residual of the equation to be solved

      ISTR    = IMAX + 2*IN
      JSTR    = JMAX + 2*JN
      KIND    = 1
c      KIND    = 2
c      IF(KMAX == 1) KIND = 1
      IL      = ISTR*JSTR
C ... Internal points
      DO 1000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1000 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = (AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+
     + AS(IJ)*TEMP(IJ-ISTR)+AE(IJ)*TEMP(IJ+1)+AW(IJ)*TEMP(IJ-1)+
     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ) )
c     + AT(IJ)*TEMP(IJ+IL)+AD(IJ)*TEMP(IJ-IL)+Q(IJ) +.005*TEMP(IJ))
c      if(i == 10 .and.j == 1 .and. k == 1) 
c     + write(808,*) ij,res(ij),temp(ij),q(ij),an(ij),as(ij)
1000  CONTINUE
C ... Sides
      DO 1100 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1100 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = 0
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AE(IJ)*TEMP(IJ+1)+Q(IJ)
1100  CONTINUE
      DO 1200 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1200 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      I       = IMAX + 1
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AW(IJ)*TEMP(IJ-1)+Q(IJ)
1200  CONTINUE
      DO 1300 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = 0
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1300 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AN(IJ)*TEMP(IJ+ISTR)+Q(IJ)
1300  CONTINUE
      DO 1400 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      J       = JMAX + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1400 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AS(IJ)*TEMP(IJ-ISTR)+Q(IJ)
1400  CONTINUE
      K       = 0
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1500 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1500 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AT(IJ)*TEMP(IJ+IL)+Q(IJ)
1500  CONTINUE
      K       = KMAX + 1
      JJ      = (KN+KIND*K-KIND)*IL
      DO 1600 J = 1,JMAX
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1600 I = 1,IMAX
      IJ      = I + II
      RES(IJ) = AP(IJ)*TEMP(IJ)+AD(IJ)*TEMP(IJ-IL)+Q(IJ)
c      if(mult ==  0) write(6,*) res(ij),ij,ap(ij),ad(ij),temp(ij),
c     + temp(ij-il),q(ij)
1600  CONTINUE

C ... Convergence
      CONV    = 0.
      CONV2   = 0.
c      if(iconv == 1 .and. mult == 50) then
c      WRITE(71,*)
c      WRITE(71,*) 'RESIDUAL'
c      CALL IPRINT1(RES,IMAX,JMAX,KMAX,IN,JN,KN,71)
c      ENDIF 
      DO 8000 K = 1,KMAX
      JJ      = (KN+KIND*K-KIND)*IL
      DO 8000 J = 0,JMAX+1
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 8000 I = 0,IMAX+1
      CONV    = CONV + RES(I+II)**2
      CONV2   = MAX(CONV2,ABS(RES(I+II))) !Maximum value
c      CONV2   = CONV2 + RES(I+II) ! Mielenkiintoinen ilmiö
8000  CONTINUE
9000  FORMAT(/12F8.3)
      IF(CONV > 1.E-30) THEN
      CONV    = SQRT(CONV)/((IMAX+2)*(JMAX+2)*(KMAX+2))
      ELSE
      CONV    = 0.
      ENDIF
      IF(ICONV3 == 1) WRITE(188,*) NGL,ICONV,CONV,CONV2
c      IF(ICONV == 2) WRITE(189,*) CONV,CONV2
      END SUBROUTINE RESAMG_DPF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE AASMA_DP(AP,AN,AS,AE,AW,AT,AD,APP,APN,APS,APE,APW,APT,
     + APD,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN,ISO,IPI,IFILE,MCYCLE)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,I,J,K,IMAX2,JMAX2,KMAX2,IN,JN,KN,ISO,
     + IPI,IFILE,MCYCLE,ISTR,JSTR,IMAX22,JMAX22,IMJM,IMCJMC,KIND,IMJM2,
     + KKC,KKF,IIC,IIF,IP,IP1,IPN,IPN1,IPK,I2,IP1N,IP1K,IPNK,IPN1K

      REAL :: AP(ISO),AN(ISO),AS(ISO),AW(ISO),AE(ISO),AT(ISO),
     + AD(ISO),APP(IPI),APN(IPI),APS(IPI),APE(IPI),APW(IPI),APT(IPI),
     + APD(IPI)
      REAL :: SUMA

      APP     = 0.
      APN     = 0.
      APS     = 0.
      APE     = 0.
      APW     = 0.
      APT     = 0.
      APD     = 0.

      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2

      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTR*JSTR
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IMJM2   = IMJM
      IF(KMAX == 1) THEN
        KIND = 1
        IMJM2= 0
      ENDIF

C ... Internal points (Last rows for odd numbers)

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5000 I = 1,IMAX2
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IPN1    = IPN  + 1
      IPK     = IP   + IMJM2                       
      IP1K    = IP1  + IMJM2                       
      IPNK    = IPN  + IMJM2                       
      IPN1K   = IPN1 + IMJM2                      
      I2      = IIC  + I
      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IPN1) +
     +          AP(IPK)+ AP(IP1K)+ AP(IPNK)+ AP(IPN1K)+
     +          AN(IP) + AN(IP1) + AN(IPK) + AN(IP1K) +
     +          AS(IPN)+ AS(IPN1)+ AS(IPNK)+ AS(IPN1K)+
     +          AE(IP) + AE(IPN) + AE(IPK) + AE(IPNK) +
     +          AW(IP1)+ AW(IPN1)+ AW(IP1K)+ AW(IPN1K)+
     +          AT(IP) + AT(IP1) + AT(IPN) + AT(IPN1) +
     +          AD(IPK)+ AD(IP1K)+ AD(IPNK)+ AD(IPN1K)
      APW(I2) = AW(IP) + AW(IPN) + AW(IPK) + AW(IPNK)
      APE(I2) = AE(IP1)+ AE(IPN1)+ AE(IP1K)+ AE(IPN1K)
      APN(I2) = AN(IPN)+ AN(IPN1)+ AN(IPNK)+ AN(IPN1K)
      APS(I2) = AS(IP) + AS(IP1) + AS(IPK) + AS(IP1K)
      APT(I2) = AT(IPK)+ AT(IP1K)+ AT(IPNK)+ AT(IPN1K)
      APD(I2) = AD(IP) + AD(IP1) + AD(IPN) + AD(IPN1)
      SUMA    =-APW(I2)-APE(I2)-APN(I2)-APS(I2)-APT(I2)-APD(I2)
c     if(app(i2) <=  1.0001*suma) write(443,*) i,j,k,app(i2),suma
      APP(I2) = MAX(APP(I2),SUMA) ! Should be useless
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IP1) + AP(IPK) + AP(IP1K)
      APN(I2) = AN(IP) + AN(IP1) + AN(IPK) + AN(IP1K)
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IP1) + AP(IPK) + AP(IP1K)
      APS(I2) = AS(IP) + AS(IP1) + AS(IPK) + AS(IP1K)
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IPN) + AP(IPK) + AP(IPNK)
      APE(I2) = AE(IP) + AE(IPN) + AE(IPK) + AE(IPNK)
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      APP(I2) = AP(IP) + AP(IPN) + AP(IPK) + AP(IPNK)
      APW(I2) = AW(IP) + AW(IPN) + AW(IPK) + AW(IPNK)
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM !+ IMJM !  Eturivi 
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IP1N)
      APT(I2) = AT(IP) + AT(IP1) + AT(IPN) + AT(IP1N)
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  ! Takarivi  
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      APP(I2) = AP(IP) + AP(IP1) + AP(IPN) + AP(IP1N)
      APD(I2) = AD(IP) + AD(IP1) + AD(IPN) + AD(IP1N)
c     IF(APP(I2) < APD(I2)) write(444,*)'apd',app(i2),apd(i2)
5600  CONTINUE
6000  CONTINUE
C ... Corners
      DO 6100 K = 1,KMAX2
      I       = 0
      J       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2+IMAX22))
6100  CONTINUE
      DO 6200 K = 1,KMAX2
      I       = IMAX2 + 1
      J       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2+IMAX22))
6200  CONTINUE
      DO 6300 K = 1,KMAX2
      I       = 0
      J       = JMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2-IMAX22))
6300  CONTINUE
      DO 6400 K = 1,KMAX2
      I       = IMAX2 + 1
      J       = JMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2-IMAX22))
6400  CONTINUE
      DO 7100 J = 1,JMAX2
      I       = 0
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2+IMCJMC))
7100  CONTINUE
      DO 7200 J = 1,JMAX2
      I       = IMAX2 + 1
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2+IMCJMC))
7200  CONTINUE
      DO 7300 J = 1,JMAX2
      I       = 0
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+1) + APP(I2-IMCJMC))
7300  CONTINUE 
      DO 7400 J = 1,JMAX2
      I       = IMAX2 + 1
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-1) + APP(I2-IMCJMC))
7400  CONTINUE
      DO 8100 I = 0,IMAX2+1
      J       = 0
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+IMCJMC) + APP(I2+IMAX22))
8100  CONTINUE
      DO 8200 I = 0,IMAX2+1
      J       = JMAX2 + 1
      K       = 0
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2+IMCJMC) + APP(I2-IMAX22))
8200  CONTINUE
      DO 8300 I = 0,IMAX2+1
      J       = 0
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-IMCJMC) + APP(I2+IMAX22))
8300  CONTINUE 
      DO 8400 I = 0,IMAX2+1
      J       = JMAX2 + 1
      K       = KMAX2 + 1
      I2      = (KN+K-1)*IMCJMC + (JN+J-1)*IMAX22  + IN + I
      APP(I2) = .5*(APP(I2-IMCJMC) + APP(I2-IMAX22))
8400  CONTINUE

C ... Put them together and print

4400  CONTINUE
9120  CONTINUE

c      IF(MCYCLE == 1) THEN ! This can be activated
c      CALL IPRINT1(APP,IMAX2,JMAX2,KMAX2,IN,JN,KN,931)   
c      CALL IPRINT1(APN,IMAX2,JMAX2,KMAX2,IN,JN,KN,932)   
c      CALL IPRINT1(APS,IMAX2,JMAX2,KMAX2,IN,JN,KN,933)   
c      CALL IPRINT1(APE,IMAX2,JMAX2,KMAX2,IN,JN,KN,934)   
c      CALL IPRINT1(APW,IMAX2,JMAX2,KMAX2,IN,JN,KN,935)   
c      CALL IPRINT1(APT,IMAX2,JMAX2,KMAX2,IN,JN,KN,936)   
c      CALL IPRINT1(APD,IMAX2,JMAX2,KMAX2,IN,JN,KN,937)   
c      ENDIF ! MCYCLE == 1
9000  FORMAT(/16F7.3)
      END SUBROUTINE AASMA_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE RESTRI_DP(RESP,RES,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,IMAX2,JMAX2,ISTR,JSTR,IMAX22,JMAX22,IMJM,
     + IMCJMC,KIND,IMJM2,I,J,K,KKC,KKF,IIC,IIF,IP,IP1,IPN,IPN1,IPK,IP1K,
     + IPNK,IPN1K,I2,IN,JN,KN,KMAX2,IP1N,KMAX

      REAL :: RESP(*),RES(*)

C ... RESTRI_DPct the residual (Last rows for odd numbers)

      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2

      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMJM    = ISTR*JSTR
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IMJM2   = IMJM
      IF(KMAX == 1) THEN
        KIND = 1
        IMJM2= 0
      ENDIF

C ... Internal points

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5000 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      IPK     = IP    + IMJM2                       
      IP1K    = IP1   + IMJM2                       
      IPNK    = IPN   + IMJM2                       
      IPN1K   = IPN1  + IMJM2                      
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IPN1) + 
     +          RES(IPK)+ RES(IP1K)+ RES(IPNK)+ RES(IPN1K)
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPK) + RES(IP1K)
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM2
      IP1K    = IP1+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPK) + RES(IP1K)
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IPN) + RES(IPK) + RES(IPNK)
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM2
      IPNK    = IPN+ IMJM2
      I2      = IIC + I
      RESP(I2)= RES(IP) + RES(IPN) + RES(IPK) + RES(IPNK)
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM + (KIND-1)*IMJM !+ IMJM ! Eturivi  
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IP1N)
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      RESP(I2)= RES(IP) + RES(IP1) + RES(IPN) + RES(IP1N)
5600  CONTINUE
C ... Corners are put to ZERO_DP in multi
      END SUBROUTINE RESTRI_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE PROLON_DP(TEMP,TEMP2,IMAX2,JMAX2,KMAX2,KMAX,IN,JN,KN)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IMAX2,JMAX2,KMAX2,IMAX22,JMAX22,KMAX22,
     + IMCJMC,KIND,IN,JN,KN,ISTR,JSTR,IMJM,I,J,K,KKC,KKF,IIC,IIF,IP,
     + IP1,IPN,IPN1,IP1K,IPNK,IPN1K,IP1N,I2,IPK

      REAL :: TEMP(*),TEMP2(*)

C ... PROLON_DPgation

      IMAX    = 2*IMAX2
      JMAX    = 2*JMAX2
      ISTR    = 2*IMAX2 + 2*IN
      JSTR    = 2*JMAX2 + 2*JN
      IMJM    = ISTR*JSTR

      IMAX22  = IMAX2 + 2*IN
      JMAX22  = JMAX2 + 2*JN
      IMCJMC  = IMAX22*JMAX22
      KIND    = 2
      IF(KMAX == 1) KIND = 1

C ... Internal points (Last row for odd numbers)

      DO 5000 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5000 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5010 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
5010  CONTINUE
      IF(KMAX > 1) THEN
      DO 5020 I = 1,IMAX2
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPN     = IP + ISTR
      IPN1    = IPN + 1
      IPK     = IP    + IMJM                       
      IP1K    = IP1   + IMJM                       
      IPNK    = IPN   + IMJM                       
      IPN1K   = IPN1  + IMJM                      
      I2      = IIC + I
      TEMP(IPK)  = TEMP(IPK)   + TEMP2(I2)
      TEMP(IP1K) = TEMP(IP1K)  + TEMP2(I2)
      TEMP(IPNK) = TEMP(IPNK)  + TEMP2(I2)
      TEMP(IPN1K)= TEMP(IPN1K) + TEMP2(I2)
5020  CONTINUE
      ENDIF ! KMAX2 > 1
5000  CONTINUE
C ... Surface J = 0
      DO 5100 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J       = 0
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5100 I = 1,IMAX2
      IP      = 2*I + IIF - 1 + ISTR
      IP1     = IP + 1
      IPK     = IP + IMJM
      IP1K    = IP1+ IMJM
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IP1K)= TEMP(IP1K) + TEMP2(I2)
5100  CONTINUE
C ... Surface J = JMAX+1
      DO 5200 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      J = JMAX2 + 1
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5200 I = 1,IMAX2 
      IP      = 2*I + IIF - 1
      IP1     = IP + 1
      IPK     = IP + IMJM
      IP1K    = IP1+ IMJM
      I2      = IIC + I
      TEMP(IP)  = TEMP(IP)   + TEMP2(I2)
      TEMP(IP1) = TEMP(IP1)  + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IP1K)= TEMP(IP1K) + TEMP2(I2)
5200  CONTINUE
C ... Surface I = 0
      DO 5300 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = 0
      DO 5300 J = 1,JMAX2 ! eturivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF
      IPN     = IP + ISTR
      IPK     = IP + IMJM
      IPNK    = IPN+ IMJM
      I2      = IIC + I
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IPNK)= TEMP(IPNK) + TEMP2(I2)
5300  CONTINUE
C ... Surface I = IMAX + 1
      DO 5400 K = 1,KMAX2
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  
      I       = IMAX2 + 1
      DO 5400 J = 1,JMAX2 ! takarivi
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      IP      = 2*I + IIF - 1
      IPN     = IP + ISTR
      IPK     = IP + IMJM
      IPNK    = IPN+ IMJM
      I2      = IIC + I
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IPN1)= TEMP(IPN1) + TEMP2(I2)
      TEMP(IPK) = TEMP(IPK)  + TEMP2(I2)
      TEMP(IPNK)= TEMP(IPNK) + TEMP2(I2)
5400  CONTINUE
C ... Surface K = 0
      K       = 0
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM  + (KIND-1)*IMJM !+ IMJM ! Eturivi 
      DO 5500 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5500 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      TEMP(IP) = TEMP(IP)  + TEMP2(I2)
      TEMP(IP1)= TEMP(IP1) + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IP1N)= TEMP(IP1N) + TEMP2(I2)
5500  CONTINUE
C ... Surface K = KMAX + 1
      K       = KMAX2 + 1
      KKC     = (KN+K-1)*IMCJMC
      KKF     = (KN+KIND*K-KIND)*IMJM 
      DO 5600 J = 1,JMAX2
      IIC     = (JN+J-1)*IMAX22 + IN + KKC
      IIF     = (JN+2*J-2)*ISTR + IN + KKF
      DO 5600 I = 1,IMAX2
      I2      = I    + IIC
      IP      = 2*I  + IIF - 1
      IP1     = IP   + 1
      IPN     = IP   + ISTR
      IP1N    = IP1  + ISTR
      TEMP(IP) = TEMP(IP)  + TEMP2(I2)
      TEMP(IP1)= TEMP(IP1) + TEMP2(I2)
      TEMP(IPN) = TEMP(IPN)  + TEMP2(I2)
      TEMP(IP1N)= TEMP(IP1N) + TEMP2(I2)
5600  CONTINUE

      RETURN
      END SUBROUTINE PROLON_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE SUBAMG_DP(TEMP,DTEMP,KOK)

      REAL :: TEMP(KOK),DTEMP(KOK)
      TEMP  = TEMP + DTEMP

      RETURN
      END SUBROUTINE SUBAMG_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE ZERO_DP(TEMP,KOK)

      REAL :: TEMP(KOK)
      TEMP = 0.

      RETURN
      END SUBROUTINE ZERO_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      

C
C ... 3D routines for the solution of linear equations

      SUBROUTINE SYAMG_DP(IL,IU,BB,DD,AA,CC)

      IMPLICIT NONE

      INTEGER :: IL,IU,LP,I,J
      REAL :: AA(*),BB(*),CC(*),DD(*),R
C ...
C ... SUBROUTINE SY SOLVES TRIDIAGONAL SYSTEM BY ELIMINATION
C ... IL = SUBSCRIPT OF FIRST EQUATION
C ... IU = SUBSCRIPT OF LAST EQUATION
C ... BB = COEFFICIENT BEHIND DIAGONAL
C ... DD = COEFFICIENT ON DIAGONAL
C ... COEFFICIENT AHEAD OF DIAGONAL
C ... ELEMENT OF CONSTANT VECTOR
C ...
C ... ESTABLISH UPPER TRIANGULAR MATRIX
C ...
      LP = IL+1
      DO 10 I = LP,IU
      R = BB(I)/DD(I-1)
      DD(I) = DD(I)-R*AA(I-1)
   10 CC(I) = CC(I)-R*CC(I-1)
C ...
C ... BACK SUBSTITUTION
C ...
      CC(IU) = CC(IU)/DD(IU)
      DO 20 I = LP,IU
      J = IU-I+IL
   20 CC(J) = (CC(J)-AA(J)*CC(J+1))/DD(J)
C ...
C ... SOLUTION STORED IN CC
C ...
      RETURN
      END SUBROUTINE SYAMG_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      
      SUBROUTINE LGSAMG_DP(DPCORR,AP,AE,AW,AS,AN,AT,AD,RMASS,AA,BB,CC,
     +DD,IMAX,JMAX,KMAX,IN,JN,KN,ITERM,ICON,KIERR,INDD)

      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IN,JN,KN,ISTR,JSTR,IL,NTOT,ITER,
     +           I,J,K,L,N,ICON,KIERR,INDD,M,II,JJ,IJ,ITERM

      REAL :: DPCORR(*),AP(*),AE(*),AW(*),AS(*),AN(*),AT(*),
     + AD(*),RMASS(*),DD(*),AA(*),BB(*),CC(*),DPCORM,REL,DPCP,DPCM,
     + DPCT,DPCD
      REAL , ALLOCATABLE :: ED(:)  ! Voi aktivoida

C ... LINE GAUSS-SEIDEL ITERATION

      ISTR   = IMAX + 2*IN
      JSTR   = JMAX + 2*JN
      IL     = ISTR*JSTR
      NTOT   = ISTR*JSTR*(KMAX+2*KN)
      ALLOCATE (ED(NTOT))
c      DO 8000 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(71,9001) (DPCORR(I+II),I=1,IMAX+2)
9001  FORMAT(/16F7.3)
c8000  CONTINUE

      ED     = 0.

      DO 9000 ITER = 1,ITERM
      DPCORM = 0.
      REL    = 0.0

      IF(INDD == 13) go to 3999
      DO 2000 K = 1,KMAX!+1 ! 2D
      JJ      = (KN+K-1)*IL
      DO 2000 J = 1,JMAX!+1 ! Sweeppaa jk-suunnassa
      II      = (JN+J-1)*ISTR + JJ + IN
      DO 1000 I = 0,IMAX+1 ! Aidosti yli laitojen 
      M       = I + 1
      IJ      = I + II
      IF(J <= JMAX) THEN
         DPCP = DPCORR(IJ+ISTR)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(J > 0) THEN
         DPCM = DPCORR(IJ-ISTR)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(K <= KMAX) THEN
         DPCT = DPCORR(IJ+IL)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(K > 0) THEN
         DPCD = DPCORR(IJ-IL)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(M)   = -RMASS(IJ)-AN(IJ)*DPCP-AS(IJ)*DPCM
     +                    -AT(IJ)*DPCT-AD(IJ)*DPCD
      AA(M)   = AE(IJ)
      DD(M)   = AP(IJ)  
      BB(M)   = AW(IJ)
1000  CONTINUE
c      CALL SYAMG(2,IMAX+3,BB,DD,AA,CC)
      CALL SYAMG_DP(1,IMAX+2,BB,DD,AA,CC)
      DO 1500 I = 0,IMAX+1 ! Aidosti yli laitojen 
      M       = I + 1
      IJ      = I + II
      DPCORR(IJ) = CC(M)
1500  CONTINUE
2000  CONTINUE
c      DO 8001 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(72,*) ((DPCORR(I+II),i,j),I=1,IMAX+2)
9002  FORMAT(/16F7.3,2I2)
c8001  CONTINUE

3999  CONTINUE
      IF(INDD == 13) go to 5999
      DO 4000 K = 1,KMAX!+1
      JJ      = (KN+K-1)*IL
      DO 4000 I = 1,IMAX!+1 ! Toinen kierros sweeppaa ik-suunnassa
      M       = I + 1 
      DO 3000 J = 0,JMAX+1 ! Aidosti yli laitojen 
      N       = J + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      IF(I <= IMAX) THEN
         DPCP = DPCORR(IJ+1)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(I > 0) THEN
         DPCM = DPCORR(IJ-1)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(K <= KMAX) THEN
         DPCT = DPCORR(IJ+IL)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(K > 0) THEN
         DPCD = DPCORR(IJ-IL)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(N)   = -RMASS(IJ)-AE(IJ)*DPCP-AW(IJ)*DPCM
     +                    -AT(IJ)*DPCT-AD(IJ)*DPCD
      AA(N)   = AN(IJ)
      DD(N)   = AP(IJ)  
      BB(N)   = AS(IJ)
c      write(888,*) n,aa(n),dd(n),BB(N),cc(n)
3000  CONTINUE
c      write(888,*) '------------'
      CALL SYAMG_DP(1,JMAX+2,BB,DD,AA,CC)
      DO 3500 J = 0,JMAX+1 ! Aidosti yli laitojen 
      N       = J + 1
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      DPCORR(IJ) = CC(N)
cc      write(888,*) n,cc(n)
3500  CONTINUE
4000  CONTINUE
5999  CONTINUE
      
      IF(KIERR > 1) THEN
      DO 6000 J = 1,JMAX!+1 
      DO 6000 I = 1,IMAX!+1 ! Kolmas kierros sweeppaa ij-suunnassa
      DO 5000 K = 0,KMAX+1 ! Aidosti yli laitojen
      L       = K + 1
      JJ      = (KN+K-1)*IL
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      IF(I <= IMAX) THEN
         DPCP = DPCORR(IJ+1)
         ELSE
         DPCP  = 0.
      ENDIF
      IF(I > 0) THEN
         DPCM = DPCORR(IJ-1)
         ELSE
         DPCM  = 0.
      ENDIF
      IF(J <= JMAX) THEN
         DPCT = DPCORR(IJ+ISTR)
         ELSE
         DPCT  = 0.
      ENDIF
      IF(J > 0) THEN
         DPCD = DPCORR(IJ-ISTR)
         ELSE
         DPCD  = 0.
      ENDIF
      CC(L)   = -RMASS(IJ)-AN(IJ)*DPCT-AS(IJ)*DPCD
     +                    -AE(IJ)*DPCP-AW(IJ)*DPCM
      AA(L)   = AT(IJ)
      DD(L)   = AP(IJ)  
      BB(L)   = AD(IJ)

5000  CONTINUE
      CALL SYAMG_DP(1,KMAX+2,BB,DD,AA,CC)
      DO 5500 K = 0,KMAX+1 ! Aidosti yli laitojen 
      L       = K + 1
      JJ      = (KN+K-1)*IL
      II      = (JN+J-1)*ISTR + JJ + IN
      IJ      = I + II
      DPCORR(IJ) = CC(L)
5500  CONTINUE
6000  CONTINUE
      ENDIF ! KMAX > 1


c      DO 4005 J = 1,JMAX+2
c      II      = (J-1)*ISTR
c      WRITE(1073,9001) (DPCORR(I+II),I=1,IMAX+2)
c4005  CONTINUE 
c      DO 7000 J = 0,JMAX+1  ! Aidosti yli laitojen 
c      N       = J + 1 
c      II      = (JN+J-1)*ISTR + IN
c      DO 7000 I = 0,IMAX+1
c      M       = I + 1
cc      IJ      = I + II
c        DPCORM    = DPCORM + ABS(ED(IJ)-DPCORR(IJ))
c         ED(IJ)    = DPCORR(IJ)
c7000  CONTINUE
c      DPCORM  = DPCORM / ((IMAX+2)*(JMAX+2))
c      write(1088,*) dpcorm!Convergence monitor can be activated

9000  CONTINUE ! END OF ITERATION

      RETURN
      END SUBROUTINE LGSAMG_DP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      

