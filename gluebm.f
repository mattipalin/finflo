C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REDUCE(RO,ZZ,IMAX,JMAX,KMAX)

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      INTEGER :: IMAXP1, JMAXP1, KMAXP1, IMAXP2, JMAXP2, KMAXP2
      INTEGER :: ISTRID, JSTRID, KSTRID, IMAX, JMAX, KMAX
      INTEGER :: I, J, K, N, IJK, LMN, J1, J2, K1, K2, IJ, KK

      REAL :: RO(*), ZZ(*)

C ... TRANSFER OF DATA FOR THE VISUALIZATION (TO CELL VERTEX SYSTEM)

      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1
      IMAXP2 = IMAX + 2
      JMAXP2 = JMAX + 2
      KMAXP2 = KMAX + 2
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = IMAXP2*JMAXP2

      N = 0

      DO 3620 K = 0,KMAX + 1
      KK      = (KN+K-1)*ISTRID*JSTRID
      DO 3620 J = 0,JMAX + 1
      DO 3620 I = 0,IMAX + 1
      IJ      = (JN+J-1)*ISTRID + I + IN + KK
      N       = N + 1
3620  ZZ(N)   = RO(IJ)

      DO 111  K = 1,KMAXP1
      K1      = (K-1)*IMAXP1*JMAXP1
      K2      = (K-1)*KSTRID
      DO 111  J = 1,JMAXP1
      J1      = K1 + (J-1)*IMAXP1
      J2      = K2 + (J-1)*IMAXP2
      DO 111  I = 1,IMAXP1
      IJK     = J1 + I
      LMN     = J2 + I
      RO(IJK) = .125*(ZZ(LMN) + ZZ(LMN+1) + ZZ(LMN+IMAXP2)+
     &            ZZ(LMN+IMAXP2+1) + ZZ(LMN+KSTRID) + ZZ(LMN+1+KSTRID)
     &           +ZZ(LMN+IMAXP2+KSTRID) + ZZ(LMN+IMAXP2+1+KSTRID))
 111  CONTINUE

      RETURN
      END SUBROUTINE REDUCE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REDUCE_SRF(UROT,VROT,WROT,XCO,YCO,ZCO,XLE2,YLE2,ZLE2,
     & XLE3,YLE3,ZLE3,GVEX,GVEY,GVEZ,OMEGA,OMEX,OMEY,OMEZ,CENX,CENY,
     & CENZ,DT,DTOLD,IMAX,JMAX,KMAX,IN,JN,KN,TIMEL,IGRID,NGL)

      IMPLICIT NONE

      INTEGER :: IMAXP1, JMAXP1, KMAXP1, NTOT, IGRID, ICASE, NGL
      INTEGER :: ISTRID, JSTRID, KSTRID, IMAX, JMAX, KMAX, IN, JN, KN
      INTEGER :: I, J, K, N, IJK, LMN, J1, J2, K1, K2, IJ, KK
      REAL    :: GVEX, GVEY, GVEZ, OMEGA, OMEX, OMEY, OMEZ, CENX, CENY,
     &           CENZ, XAVE, YAVE, ZAVE, DT, DTOLD

      REAL    :: UROT(*), VROT(*), WROT(*)

      REAL :: XCO(*), YCO(*), ZCO(*), XLE2(*), YLE2(*),ZLE2(*), 
     &        XLE3(*), YLE3(*), ZLE3(*)

      LOGICAL :: TIMEL

C ... Calculate grid-point velocities depending on the grid motion

      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1
      NTOT   = IMAXP1*JMAXP1*KMAXP1
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = ISTRID*JSTRID

      ICASE  = 0
      IF(IGRID >= 1  .AND. IGRID <= 10) ICASE = 2
      IF(IGRID == 21 .OR.  IGRID == 22) ICASE = 2
      IF(IGRID >= 11 .AND. IGRID <= 20 .AND. TIMEL) ICASE = 3
      IF(SQRT(GVEX**2 + GVEY**2 + GVEZ**2) >= 1.E-5 .AND. .NOT.TIMEL)
     &   ICASE = 1
************************************************************************
      ICASE = 3  ! Pulttaus
************************************************************************
C ... It seems that the last alternative dominates (antakee mersu)

      SELECT CASE(ICASE)

      CASE(0)  ! Stationary grid

c         DO K = 1,KMAXP1
c            K1 = (K-1)*IMAXP1*JMAXP1
c            DO J = 1,JMAXP1
c               J1 = K1 + (J-1)*IMAXP1
c               DO I = 1,IMAXP1
c                  IJK = J1 + I
c                  UROT(IJK) = 0.
c                  VROT(IJK) = 0.
c                  WROT(IJK) = 0.
c               ENDDO
c            ENDDO
c         ENDDO

         UROT(1:NTOT) = 0.
         VROT(1:NTOT) = 0.
         WROT(1:NTOT) = 0.
         
         
      CASE(1)  ! Linear motion

c         DO K = 1,KMAXP1
c            K1 = (K-1)*IMAXP1*JMAXP1
c            K2 = (K+KN-1)*KSTRID
c            DO J = 1,JMAXP1
c               J1 = K1 + (J-1)*IMAXP1
c               J2 = K2 + (J+JN-1)*ISTRID + IN
c               DO I = 1,IMAXP1
c                  IJK = J1 + I
c                  LMN = J2 + I
c                  UROT(IJK) = GVEX
c                  VROT(IJK) = GVEY
c                  WROT(IJK) = GVEZ
c               ENDDO
c            ENDDO
c         ENDDO

         UROT(1:NTOT) = GVEX
         VROT(1:NTOT) = GVEY
         WROT(1:NTOT) = GVEZ

         
      CASE(2)  ! Rotating grid

         DO K = 1,KMAXP1
            K1 = (K-1)*IMAXP1*JMAXP1
            K2 = (K+KN-1)*KSTRID
            DO J = 1,JMAXP1
               J1 = K1 + (J-1)*IMAXP1
               J2 = K2 + (J+JN-1)*ISTRID + IN
               DO I = 1,IMAXP1
                  IJK = J1 + I
                  LMN = J2 + I
                  ZAVE = ZCO(LMN) - CENZ 
                  YAVE = YCO(LMN) - CENY
                  XAVE = XCO(LMN) - CENX
                  UROT(IJK) = OMEGA*(OMEY*ZAVE - OMEZ*YAVE)
                  VROT(IJK) = OMEGA*(OMEZ*XAVE - OMEX*ZAVE)
                  WROT(IJK) = OMEGA*(OMEX*YAVE - OMEY*XAVE)
               ENDDO
            ENDDO
         ENDDO


      CASE(3)  ! General grid movement. Needs two time levels

         IJK = 0
************************************************************************
*         if(NGL == 6) write(777,*) 'Block 6, I = 4, K = 28'
************************************************************************
         DO K = 1,KMAXP1
            DO J = 1,JMAXP1
               DO I = 1,IMAXP1
                  IJK = IJK + 1
                  LMN = I+IN + (J+JN-1)*ISTRID + (K+KN-1)*KSTRID
                  UROT(IJK) = 1.5*(XCO(LMN)  - XLE2(LMN))/DT -
     +                         .5*(XLE2(LMN) - XLE3(LMN))/DTOLD
                  VROT(IJK) = 1.5*(YCO(LMN)  - YLE2(LMN))/DT -
     +                         .5*(YLE2(LMN) - YLE3(LMN))/DTOLD
                  WROT(IJK) = 1.5*(ZCO(LMN)  - ZLE2(LMN))/DT -
     +                         .5*(ZLE2(LMN) - ZLE3(LMN))/DTOLD
************************************************************************
*                  if(NGL == 6 .AND. I == 4 .AND. K == 28) THEN
*                     write(777,*) J,REAL(YCO(LMN),4),REAL(YLE2(LMN),4),
*     +                              REAL(YLE3(LMN),4),                
*     +                              REAL(VROT(IJK),4)                
*
*                  endif
************************************************************************
               ENDDO
            ENDDO
         ENDDO

         
      END SELECT

      RETURN
      END SUBROUTINE REDUCE_SRF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REDUCE_SRFO(RO,ZZ,IMAX,JMAX,KMAX,IN,JN,KN,IDIR)

      IMPLICIT NONE

      INTEGER :: IMAXP1, JMAXP1, KMAXP1, IMAXP2, JMAXP2, KMAXP2, IDIR
      INTEGER :: ISTRID, JSTRID, KSTRID, IMAX, JMAX, KMAX, IN, JN, KN
      INTEGER :: I, J, K, N, IJK, LMN, J1, J2, K1, K2, IJ, KK

      REAL :: RO(*), ZZ(*)

C ... TRANSFER OF 2D DATA FOR THE VISUALIZATION (TO CELL VERTEX SYSTEM)

      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1
      IMAXP2 = IMAX + 2
      JMAXP2 = JMAX + 2
      KMAXP2 = KMAX + 2
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = IMAXP2*JMAXP2

      N = 0

      DO 3620 K = 0,KMAX + 1
      KK      = (KN+K-1)*ISTRID*JSTRID
      DO 3620 J = 0,JMAX + 1
      DO 3620 I = 0,IMAX + 1
      IJ      = (JN+J-1)*ISTRID + I + IN + KK
3620  CONTINUE

      SELECT CASE(IDIR)

      CASE(1) ! I-direction

      DO 111  K = 1,KMAXP1
      K1      = (K-1)*IMAXP1*JMAXP1
      K2      = (K-1)*KSTRID
      DO 111  J = 1,JMAXP1
      J1      = K1 + (J-1)*IMAXP1
      J2      = K2 + (J-1)*IMAXP2
      DO 111  I = 1,IMAXP1
      IJK     = J1 + I
      LMN     = J2 + I
      RO(IJK) = .25*(ZZ(LMN)  + ZZ(LMN+IMAXP2)+
     &             ZZ(LMN+KSTRID) + ZZ(LMN+IMAXP2+KSTRID))
 111  CONTINUE

      CASE(2) ! J-direction

      DO 222  K = 1,KMAXP1
      K1      = (K-1)*IMAXP1*JMAXP1
      K2      = (K-1)*KSTRID
      DO 222  J = 1,JMAXP1
      J1      = K1 + (J-1)*IMAXP1
      J2      = K2 + (J-1)*IMAXP2
      DO 222  I = 1,IMAXP1
      IJK     = J1 + I
      LMN     = J2 + I
      RO(IJK) = .25*(ZZ(LMN) + ZZ(LMN+1)+
     &             ZZ(LMN+KSTRID) + ZZ(LMN+1+KSTRID))
 222  CONTINUE

      CASE(3) ! K-direction

      DO 333  K = 1,KMAXP1
      K1      = (K-1)*IMAXP1*JMAXP1
      K2      = (K-1)*KSTRID
      DO 333  J = 1,JMAXP1
      J1      = K1 + (J-1)*IMAXP1
      J2      = K2 + (J-1)*IMAXP2
      DO 333  I = 1,IMAXP1
      IJK     = J1 + I
      LMN     = J2 + I
      RO(IJK) = .25*(ZZ(LMN) + ZZ(LMN+1) + ZZ(LMN+IMAXP2)+
     &            ZZ(LMN+IMAXP2+1))
 333  CONTINUE

      END SELECT

      RETURN
      END SUBROUTINE REDUCE_SRFO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRITER(LAITE,IFORM,RO,APUP2,IMAX,JMAX,KMAX,IN,JN,KN)

      REAL :: RO(*), APUP2(*)

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      N       = 0
      DO 3620 K = 0,KMAX + 1
      KK      = (KN+K-1)*ISTRID*JSTRID
      DO 3620 J = 0,JMAX + 1
      DO 3620 I = 0,IMAX + 1
      IJ      = (JN+J-1)*ISTRID + I + IN + KK
      N       = N + 1
3620  APUP2(N)= RO(IJ)
      IF(IFORM == 0) THEN
           WRITE(LAITE)   (APUP2(N),N=1,(IMAX+2)*(JMAX+2)*(KMAX+2))
      ELSE
           WRITE(LAITE,*) (APUP2(N),N=1,(IMAX+2)*(JMAX+2)*(KMAX+2))
      ENDIF

      RETURN
      END SUBROUTINE WRITER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C***********************************************************************
C ... ROUTINES FOR TRANSFORMING DATA TO GHOST CELLS AT OUTPUT STAGE
C***********************************************************************

      SUBROUTINE VIEWRF(OMX,OMY,OMZ,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,
     2 KBOT,KTOP,INRE,JNRE,KNRE,ICON,NPATCH)

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: OMX(*), OMY(*), OMZ(*)
      INTEGER   :: ICON(IC9,*)
      INTEGER   :: IMAX, JMAX, KMAX, IN, JN, KN, INRE, JNRE, KNRE
C
C ... REFLECT THE BOUNDARY VALUES BELOW THE SOLID SURFACE (FOR FLOVIS)
C
C ... MODIFIED VIEWRF. DO REFLECTION ONLY FOR VORTICITY COMPONENTS

      IN = INRE
      JN = JNRE
      KN = KNRE

      JSTRID  = JMAX + 2*JN
      ISTRID  = IMAX + 2*IN
      IL      = ISTRID*JSTRID

      DO 7000 IP = 1,NPATCH
      IBC     = ICON(1,IP)
      IF(IBC >= 7 .AND. IBC <= 10) THEN
         IFACE   = ICON(3,IP)
         IXLO    = ICON(4,IP)
         IXUP    = ICON(5,IP)
         IYLO    = ICON(6,IP)
         IYUP    = ICON(7,IP)
         IDIR    = 1
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
            I1   = IXLO
            I2   = IXUP
            J1   = IYLO
            J2   = IYUP
C ... EXTEND PATCHES 
            IF(I2 == JMAX) I2 = JMAX + 1
            IF(J2 == KMAX) J2 = KMAX + 1
         ENDIF
         IF(IFACE == 4) THEN
            K   = ITOP
            IF(IBC == 7) K = IMAX
         ENDIF
C ... ETA-DIRECTION
         IF(IFACE == 2 .OR. IFACE == 5) THEN
            ISTR = IL
            JSTR = 1
            KSTR = ISTRID
            IN   = KNRE
            JN   = INRE
            KN   = JNRE
            K    = JBOT
            I1   = IYLO
            I2   = IYUP
            J1   = IXLO
            J2   = IXUP
C ... EXTEND PATCHES 
            IF(I2 == KMAX) I2 = KMAX + 1
            IF(J2 == IMAX) J2 = IMAX + 1
         ENDIF
         IF(IFACE == 5) THEN
            K   = JTOP
            IF(IBC == 7) K = JMAX
         ENDIF
C ... ZETA DIRECTION
         IF(IFACE == 3 .OR. IFACE == 6) THEN
            IN   = INRE
            JN   = JNRE
            KN   = KNRE
            ISTR = 1
            JSTR = ISTRID
            KSTR = IL
            K = KBOT
            I1   = IXLO
            I2   = IXUP
            J1   = IYLO
            J2   = IYUP
C ... EXTEND PATCHES 
            IF(I2 == IMAX) I2 = IMAX + 1
            IF(J2 == JMAX) J2 = JMAX + 1
         ENDIF
         IF(IFACE == 6) THEN
            K   = KTOP
            IF(IBC == 7) K = KMAX
         ENDIF

         IF(IBC == 7.AND.IFACE <= 3) K = 1
C ... EXTEND PATCHES 
         IF(I1 == 1) I1 = 0
         IF(J1 == 1) J1 = 0

C **********************************************************************
         DO 1000 J = J1,J2
         KK      = (JN+J-1)*JSTR + (KN+K-1)*KSTR
         DO 1000 I = I1,I2
            L        = 1 + (IN+I-1)*ISTR + KK
            OMX(L-KSTR*IDIR) = 2.*OMX(L) - OMX(L+KSTR*IDIR)
            OMY(L-KSTR*IDIR) = 2.*OMY(L) - OMY(L+KSTR*IDIR)
            OMZ(L-KSTR*IDIR) = 2.*OMZ(L) - OMZ(L+KSTR*IDIR)
 1000    CONTINUE
C **********************************************************************

      ENDIF
7000  CONTINUE

      RETURN
      END SUBROUTINE VIEWRF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE VIEWSC(FI,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     2 IDI1,IDI2,IDI3,INRE,JNRE,KNRE,ICON,NPATCH,ITURB,MAXSB,NSCAL,IWAY)

      USE NS3CO, ONLY : IC9
      
      INTEGER   :: IMAX, JMAX, KMAX, IN, JN, KN, INRE, JNRE, KNRE
      DIMENSION :: FI(MAXSB,MAX(1,NSCAL))
      INTEGER   :: ICON(IC9,*)
      LOGICAL   :: REYN
C
C ... REFLECT THE SCALAR VALUES BELOW THE SOLID SURFACE (FOR FLOVIS)
C ... PPR 16.12.94
C
      IN = INRE
      JN = JNRE
      KN = KNRE

      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      IL     = ISTRID*JSTRID

      DO 313 NS = 1,NSCAL

         IF (ITURB >= 10) THEN
            REYN = .TRUE.
         ELSE
            REYN = .FALSE.
         ENDIF

         DO 7000 IP = 1,NPATCH

            IBC = ICON(1,IP)

            IF(IBC >= 7 .AND. IBC <= 10) THEN

               IFACE   = ICON(3,IP)
               IXLO    = ICON(4,IP)
               IXUP    = ICON(5,IP)
               IYLO    = ICON(6,IP)
               IYUP    = ICON(7,IP)
               IDIR    = 1

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
                  I1   = IXLO
                  I2   = IXUP
                  J1   = IYLO
                  J2   = IYUP
                  ID   = IDI1

C ... EXTEND PATCHES 

                  IF(I2 == JMAX) I2 = JMAX + 1

                  IF(J2 == KMAX) J2 = KMAX + 1

               ENDIF

               IF(IFACE == 4) THEN
                  K = ITOP
                  IF(IBC == 7) K = IMAX
               ENDIF

C ... ETA-DIRECTION

               IF(IFACE == 2 .OR. IFACE == 5) THEN

                  IN   = KNRE
                  JN   = INRE
                  KN   = JNRE
                  ISTR = IL
                  JSTR = 1
                  KSTR = ISTRID
                  K    = JBOT
                  I1   = IYLO
                  I2   = IYUP
                  J1   = IXLO
                  J2   = IXUP
                  ID   = IDI2

C ... EXTEND PATCHES 

                  IF(I2 == KMAX) I2 = KMAX + 1
                  IF(J2 == IMAX) J2 = IMAX + 1

               ENDIF

               IF(IFACE == 5) THEN
                  K = JTOP
                  IF(IBC == 7) K = JMAX
               ENDIF

C ... ZETA DIRECTION

               IF(IFACE == 3 .OR. IFACE == 6) THEN 

                  IN   = INRE
                  JN   = JNRE
                  KN   = KNRE
                  ISTR = 1
                  JSTR = ISTRID
                  KSTR = IL
                  K    = KBOT
                  I1   = IXLO
                  I2   = IXUP
                  J1   = IYLO
                  J2   = IYUP
                  ID   = IDI3

C ... EXTEND PATCHES 

                  IF(I2 == IMAX) I2 = IMAX + 1
                  IF(J2 == JMAX) J2 = JMAX + 1

               ENDIF

               IF(IFACE == 6) THEN
                  K = KTOP
                  IF(IBC == 7) K = KMAX
               ENDIF

               IF(IBC == 7 .AND. IFACE <= 3) K = 1

C ... EXTEND PATCHES 

               IF(I1 == 1) I1 = 0
               IF(J1 == 1) J1 = 0

               IWAY2 = IWAY          

               IF(IWAY2 == 0) THEN

                  IF(((ID /= 0 .AND. IFACE <= 3) .OR.
     &                 (ID > K .AND. IFACE >= 4))
     &                         .AND. IBC /= 7 .AND. REYN) THEN
                     IWAY2 = 1
                  ELSE
                     IWAY2 = 2
                  ENDIF

               ENDIF 

               IF(IWAY2 == 1) THEN  ! Zero at the wall

               DO J = J1,J2
                  KK = (JN+J-1)*JSTR + (KN+K-1)*KSTR
                  DO I = I1,I2
                     L = 1 + (IN+I-1)*ISTR + KK
                     FI(L-KSTR*IDIR,NS)  = -FI(L,NS)
                  ENDDO
               ENDDO
         
            ELSE IF(IWAY2 == 2) THEN  ! Extrapolate from the domain

               DO J = J1,J2
                  KK = (JN+J-1)*JSTR + (KN+K-1)*KSTR
                  DO I = I1,I2
                     L = 1 + (IN+I-1)*ISTR + KK
                     FI(L-KSTR*IDIR,NS) = 2.*FI(L,NS)-FI(L+KSTR*IDIR,NS)
                  ENDDO
               ENDDO
            
            ELSE IF(IWAY2 == 3) THEN  ! The same as in the first cell

               DO J = J1,J2
                  KK = (JN+J-1)*JSTR + (KN+K-1)*KSTR
                  DO I = I1,I2
                     L = 1 + (IN+I-1)*ISTR + KK
                     FI(L-KSTR*IDIR,NS)  = FI(L,NS)
                  ENDDO
               ENDDO
      
            ELSE

               WRITE(*,*) 'Something wrong in VIEWSC, returning...'

            ENDIF

         ENDIF

 7000 CONTINUE

 313  CONTINUE

      RETURN
      END SUBROUTINE VIEWSC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE VIEWTE(TEMP,TWALL,IHF,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,
     2 JTOP,KBOT,KTOP,IDI1,IDI2,IDI3,INRE,JNRE,KNRE,
     3     ICON,NPATCH,ITURB,MAXSB,NSCAL,IWAY)

      USE NS3CO, ONLY : LN, IC9

      IMPLICIT NONE

      REAL    :: TEMP(*),TWALL(*)

      INTEGER :: ICON(IC9,*), IHF(*)

      INTEGER :: IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     &           IDI1,IDI2,IDI3,IN,JN,KN,INRE, JNRE, KNRE,
     &           NPATCH,ITURB,MAXSB,NSCAL,IWAY,
     &           ISTRID,JSTRID,KSTRID,IL,IP,IPL,IBC,KX1,KX2,KY1,KY2,
     &           IFACE,IDIR,ISTR,JSTR,KSTR,I,J,K,KA,NR,L,NP,LB,IJ

C ... REFLECT THE SCALAR VALUE BELOW THE SOLID SURFACE USING WALL VALUE
C ... TSii 16.11.04

      IN = INRE
      JN = JNRE
      KN = KNRE

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID

C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH

      IBC     = ICON(1,IP)

      IF(IBC. EQ. 15) THEN

      IPL   = ICON(2,IP) ! proces local patch number
      IFACE = ICON(3,IP)
C ... So-called standard wall-indeces
         KX1  = ICON(4,IP)
         KX2  = ICON(5,IP)
         KY1  = ICON(6,IP)
         KY2  = ICON(7,IP)
      IF(IFACE ==  2 .OR. IFACE == 5) THEN
         KY1  = ICON(4,IP)
         KY2  = ICON(5,IP)
         KX1  = ICON(6,IP)
         KX2  = ICON(7,IP)
      ENDIF

      IDIR    = 1
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
         IF(IFACE == 4) K   = ITOP

C ... ETA-DIRECTION

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

         IN   = KNRE
         JN   = INRE
         KN   = JNRE

         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = JBOT
         IF(IFACE == 5) K   = JTOP

C ... ZETA DIRECTION

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         IN   = INRE
         JN   = JNRE
         KN   = KNRE

         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = KBOT
         IF(IFACE == 6) K   = KTOP
      ENDIF
         
C ... Use the wall value to get the correct reflected value
       
      KA = (KN+K-1)*KSTR

         DO 500 J = KY1,KY2
         IJ = (JN+J-1)*JSTR + KA 
         NR = (LN+J-KY1)*(KX2-KX1+1+2*LN) - KX1 + LN + IHF(IPL) 
         DO 500 I = KX1,KX2
            L  = 1 + (IN+I-1)*ISTR + IJ      ! Cell index
            NP = I + NR                      ! Patch index with LN ghost cells
            LB = L - IDIR*KSTR
            TEMP(LB)= 2.*TWALL(NP) - TEMP(L)
 500    CONTINUE

      ENDIF                     ! IBC == 15
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE VIEWTE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRDTE(TEMP,INTO,IHF,IMAX,JMAX,KMAX,IBOT,
     &                  ITOP,JBOT,JTOP,KBOT,KTOP,M,IDI1,IDI2,IDI3,
     &                  INRE,JNRE,KNRE,ICON,NPATCH,
     &                  ITURB,MAXSB,NSCAL,IWAY)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: ICON(IC9,*), IHF(*)

      INTEGER :: IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     &           IDI1,IDI2,IDI3,IN,JN,KN,INRE, JNRE, KNRE, M,
     &           NPATCH,ITURB,MAXSB,NSCAL,IWAY,
     &           ISTRID,JSTRID,KSTRID,IL,IP,IPL,KX1,KX2,KY1,KY2,IDIR,
     &           IBC,IFACE,ISTR,JSTR,KSTR,I,J,K,KA,NR,L,NP,LB,IJ,INTO,
     &           IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2

      REAL :: TEMP(*)

C ... REFLECT THE SCALAR VALUE BELOW THE 'INTO' SURFACE
C ... TSii 5.9.06

      IN = INRE
      JN = JNRE
      KN = KNRE

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      IL     = ISTRID*JSTRID
      
C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
       
      IBC = ICON(1,IP)

      IF(IBC == INTO) THEN

      CALL PATCHE(ICON,IP,M,IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2)

      IPL   = ICON(2,IP) ! proces local patch number
      IFACE = ICON(3,IP)

C ... So-called standard wall-indeces
         KX1 = IA1
         KX2 = IM1
         KY1 = JA1
         KY2 = JM1

      IF(IFACE == 2 .OR. IFACE == 5) THEN
         KY1 = IA1
         KY2 = IM1
         KX1 = JA1
         KX2 = JM1
      ENDIF

      IDIR = 1
      IF(IFACE >= 4) IDIR = -1

C ... XI-DIRECTION

      IF(IFACE == 1 .OR. IFACE == 4) THEN

         IN = JNRE
         JN = KNRE
         KN = INRE

         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K = IMAX

C ... ETA-DIRECTION

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

         IN = KNRE
         JN = INRE
         KN = JNRE

         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K = JMAX

C ... ZETA DIRECTION

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         IN = INRE
         JN = JNRE
         KN = KNRE

         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K = KMAX
      ENDIF
         
C ... Use the wall value to get the correct reflected value
       
      KA      = (KN+K-1)*KSTR

         DO 500 J = KY1,KY2
         IJ       = (JN+J-1)*JSTR + KA 
         DO 500 I = KX1,KX2
            L     = 1 + (IN+I-1)*ISTR + IJ  ! Cell index
            LB    = L - IDIR*KSTR
            TEMP(LB) = TEMP(L)
            TEMP(LB-IDIR*KSTR) = TEMP(L+IDIR*KSTR)
 500    CONTINUE

      ENDIF                     ! IBC == INTO
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE MIRDTE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRDTE_SP(TEMP,INTO,IHF,IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,
     2 JTOP,KBOT,KTOP,M,IDI1,IDI2,IDI3,INRE,JNRE,KNRE,
     3 ICON,NPATCH,ITURB,MAXSB,NSCAL,IWAY)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX,JMAX,KMAX,IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     2           IDI1,IDI2,IDI3,IN,JN,KN,INRE,JNRE,KNRE,M,
     3           NPATCH,ITURB,MAXSB,NSCAL,IWAY,
     3           ISTRID,JSTRID,KSTRID,IL,IP,IPL,IBC,IFACE,
     4           KX1,KX2,KY1,KY2,IDIR,ISTR,
     5           JSTR,KSTR,I,J,K,KA,NR,L,NP,LB,IJ,INTO,
     6           IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2
      
      REAL :: TEMP(*)
      INTEGER :: ICON(IC9,*), IHF(*)

C ... REFLECT THE SCALAR VALUE BELOW THE 'INTO' SURFACE
C ... TSii 5.9.06

      IN = INRE
      JN = JNRE
      KN = KNRE

      KSTRID = KMAX + 2*KN
      JSTRID = JMAX + 2*JN
      ISTRID = IMAX + 2*IN
      IL     = ISTRID*JSTRID
      
C ... LOOP OVER THE PATCHES OF THE BLOCK 

      DO 7000 IP = 1,NPATCH
       
      IBC = ICON(1,IP)

      IF(IBC == INTO) THEN

      CALL PATCHE(ICON,IP,M,IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2)

      IPL   = ICON(2,IP) ! proces local patch number
      IFACE = ICON(3,IP)

C ... So-called standard wall-indeces
         KX1 = IA1
         KX2 = IM1
         KY1 = JA1
         KY2 = JM1

      IF(IFACE == 2 .OR. IFACE == 5) THEN
         KY1 = IA1
         KY2 = IM1
         KX1 = JA1
         KX2 = JM1
      ENDIF
      
      IDIR = 1
      IF(IFACE >= 4) IDIR = -1

C ... XI-DIRECTION

      IF(IFACE == 1 .OR. IFACE == 4) THEN

         IN = JNRE
         JN = KNRE
         KN = INRE

         ISTR = ISTRID
         JSTR = IL
         KSTR = 1
         K    = 1
         IF(IFACE == 4) K = IMAX

C ... ETA-DIRECTION

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN

         IN = KNRE
         JN = INRE
         KN = JNRE

         ISTR = IL
         JSTR = 1
         KSTR = ISTRID
         K    = 1
         IF(IFACE == 5) K = JMAX

C ... ZETA DIRECTION

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN

         IN = INRE
         JN = JNRE
         KN = KNRE

         ISTR = 1
         JSTR = ISTRID
         KSTR = IL
         K    = 1
         IF(IFACE == 6) K = KMAX
      ENDIF
         
C ... Use the wall value to get the correct reflected value
       
      KA = (KN+K-1)*KSTR

         DO 500 J = KY1,KY2
         IJ       = (JN+J-1)*JSTR + KA 
         DO 500 I = KX1,KX2
            L     = 1 + (IN+I-1)*ISTR + IJ  ! Cell index
            LB    = L - IDIR*KSTR
            TEMP(LB) = TEMP(L)
            TEMP(LB-IDIR*KSTR) = TEMP(L+IDIR*KSTR)
 500    CONTINUE

      ENDIF                     ! IBC == INTO
7000  CONTINUE ! LOOP OVER THE PATCHES

      RETURN
      END SUBROUTINE MIRDTE_SP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C***********************************************************************
C ... NEW ROUTINES FOR GRID PROPERTIES
C***********************************************************************

      SUBROUTINE GRIDEX(XCO,YCO,ZCO,IMAX,JMAX,KMAX,ICON,NPATCH,XCOL,
     &                  IN,JN,KN)

      USE NS3CO, ONLY : IC9
      
C ... THIS SUB. EXTRAPOLATES CORNER COORDINATES OVER BOUNDARIES.

      REAL    :: XCO(*), YCO(*), ZCO(*)
      REAL    :: XCOL(3,4,*)
      INTEGER :: ICON(IC9,*)
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN

      ISTRID  = IMAX + 2*IN
      JSTRID  = JMAX + 2*JN
      KSTRID  = KMAX + 2*KN
      IL      = ISTRID*JSTRID

      IMAXP   = IMAX + 1
      JMAXP   = JMAX + 1
      KMAXP   = KMAX + 1
      ISTR    = 1
      JSTR    = ISTRID
      KSTR    = IL

C ... REFLECT BLOCK FACES AND EDGES
C ... ZETA-DIRECTION 
      CALL REFXYZ(XCO,YCO,ZCO,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,IN,JN,KN)
C ... XI-DIRECTION
      CALL REFXYZ(XCO,YCO,ZCO,JMAX,KMAX,IMAX,JSTR,KSTR,ISTR,JN,KN,IN)
C ... ETA-DIRECTION
      CALL REFXYZ(XCO,YCO,ZCO,KMAX,IMAX,JMAX,KSTR,ISTR,JSTR,KN,IN,JN)

C ... CORNER POINTS
C ... G ---- GHOST POINT     1 --- WALL POINT    2 ---- SECOND POINT IN DOMAIN
C ... LLL --- LOWER-LOWER-LOWER     LUU -- LOWER-UPPER-UPPER etc.
C ...               ... g g g G   in figure: g - normal ghost point
C ...               ... x x 1 g              x - normal domain point
C ...               ... x 2 x g              G - calculated ghost point
C ...               ... x x x g              1 - first point in domain
C ...               ... . . . .              2 - second point in domain
C ...               ... . . . . 
C ... NOTE THAT G, 1 AND 2 ARE IN DIFFERENT LEVELS
      ILG     = (0      -1+IN)*ISTR
      IL1     = (1      -1+IN)*ISTR
      IL2     = (2      -1+IN)*ISTR

      IUG     = (IMAXP+1-1+IN)*ISTR
      IU1     = (IMAXP  -1+IN)*ISTR
      IU2     = (IMAXP-1-1+IN)*ISTR

      JLG     = (0      -1+JN)*JSTR
      JL1     = (1      -1+JN)*JSTR
      JL2     = (2      -1+JN)*JSTR

      JUG     = (JMAXP+1-1+JN)*JSTR
      JU1     = (JMAXP  -1+JN)*JSTR
      JU2     = (JMAXP-1-1+JN)*JSTR

      KLG     = (0      -1+KN)*KSTR
      KL1     = (1      -1+KN)*KSTR
      KL2     = (2      -1+KN)*KSTR

      KUG     = (KMAXP+1-1+KN)*KSTR
      KU1     = (KMAXP  -1+KN)*KSTR
      KU2     = (KMAXP-1-1+KN)*KSTR


      LLLLG   = 1 + ILG + JLG + KLG
      LLLL1   = 1 + IL1 + JL1 + KL1
      LLLL2   = 1 + IL2 + JL2 + KL2

      LLLUG   = 1 + ILG + JLG + KUG
      LLLU1   = 1 + IL1 + JL1 + KU1
      LLLU2   = 1 + IL2 + JL2 + KU2

      LLULG   = 1 + ILG + JUG + KLG
      LLUL1   = 1 + IL1 + JU1 + KL1
      LLUL2   = 1 + IL2 + JU2 + KL2

      LLUUG   = 1 + ILG + JUG + KUG
      LLUU1   = 1 + IL1 + JU1 + KU1
      LLUU2   = 1 + IL2 + JU2 + KU2

      LULLG   = 1 + IUG + JLG + KLG
      LULL1   = 1 + IU1 + JL1 + KL1
      LULL2   = 1 + IU2 + JL2 + KL2

      LULUG   = 1 + IUG + JLG + KUG
      LULU1   = 1 + IU1 + JL1 + KU1
      LULU2   = 1 + IU2 + JL2 + KU2

      LUULG   = 1 + IUG + JUG + KLG
      LUUL1   = 1 + IU1 + JU1 + KL1
      LUUL2   = 1 + IU2 + JU2 + KL2

      LUUUG   = 1 + IUG + JUG + KUG
      LUUU1   = 1 + IU1 + JU1 + KU1
      LUUU2   = 1 + IU2 + JU2 + KU2

      XCO(LLLLG) = 2.*XCO(LLLL1) - XCO(LLLL2)
      YCO(LLLLG) = 2.*YCO(LLLL1) - YCO(LLLL2)
      ZCO(LLLLG) = 2.*ZCO(LLLL1) - ZCO(LLLL2)

      XCO(LLLUG) = 2.*XCO(LLLU1) - XCO(LLLU2)
      YCO(LLLUG) = 2.*YCO(LLLU1) - YCO(LLLU2)
      ZCO(LLLUG) = 2.*ZCO(LLLU1) - ZCO(LLLU2)

      XCO(LLULG) = 2.*XCO(LLUL1) - XCO(LLUL2)
      YCO(LLULG) = 2.*YCO(LLUL1) - YCO(LLUL2)
      ZCO(LLULG) = 2.*ZCO(LLUL1) - ZCO(LLUL2)

      XCO(LLUUG) = 2.*XCO(LLUU1) - XCO(LLUU2)
      YCO(LLUUG) = 2.*YCO(LLUU1) - YCO(LLUU2)
      ZCO(LLUUG) = 2.*ZCO(LLUU1) - ZCO(LLUU2)

      XCO(LULLG) = 2.*XCO(LULL1) - XCO(LULL2)
      YCO(LULLG) = 2.*YCO(LULL1) - YCO(LULL2)
      ZCO(LULLG) = 2.*ZCO(LULL1) - ZCO(LULL2)

      XCO(LULUG) = 2.*XCO(LULU1) - XCO(LULU2)
      YCO(LULUG) = 2.*YCO(LULU1) - YCO(LULU2)
      ZCO(LULUG) = 2.*ZCO(LULU1) - ZCO(LULU2)

      XCO(LUULG) = 2.*XCO(LUUL1) - XCO(LUUL2)
      YCO(LUULG) = 2.*YCO(LUUL1) - YCO(LUUL2)
      ZCO(LUULG) = 2.*ZCO(LUUL1) - ZCO(LUUL2)

      XCO(LUUUG) = 2.*XCO(LUUU1) - XCO(LUUU2)
      YCO(LUUUG) = 2.*YCO(LUUU1) - YCO(LUUU2)
      ZCO(LUUUG) = 2.*ZCO(LUUU1) - ZCO(LUUU2)

C ... CALCULATE PATCH CORNER FOR ORIENTATION

      DO 175 L = 1,NPATCH
           CALL CORNER(XCO,YCO,ZCO,XCOL(1,1,L),
     &     ICON(3,L),ICON(4,L),ICON(5,L),ICON(6,L),ICON(7,L),
     &     IMAX,JMAX,KMAX,IN,JN,KN,L)
  175 CONTINUE

      RETURN
      END SUBROUTINE GRIDEX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE REFXYZ(XCO,YCO,ZCO,IMAX,JMAX,KMAX,ISTR,JSTR,KSTR,
     +     IN,JN,KN)

      REAL :: XCO(*), YCO(*), ZCO(*)

C ... REFLECT BLOCK FACES AND EDGES

C     ****************** REFLEC BLOCK FACES *************************
C ... LBW ---- LOWER BOUNDARY        LB2 ------ SECOND POINT IN DOMAIN(BOT)
C ... LBG ---- LOWER GHOST POINT
C ... LTW ---- UPPER BOUNDARY        LT2 ------ SECOND POINT IN DOMAIN(TOP)
C ... LTG ---- UPPER GHOST POINT

      IMAXP = IMAX + 1
      JMAXP = JMAX + 1
      KMAXP = KMAX + 1

      DO 100 I = 1,IMAXP
      DO 100 J = 1,JMAXP
         LL  = 1 + (I-1+IN)*ISTR + (J-1+JN)*JSTR
         LBG = (0-1+KN)*KSTR + LL
         LBW = (1-1+KN)*KSTR + LL
         LB2 = (2-1+KN)*KSTR + LL
         
         LTG = (KMAXP+1-1+KN)*KSTR + LL
         LTW = (KMAXP  -1+KN)*KSTR + LL
         LT2 = (KMAXP-1-1+KN)*KSTR + LL

         XCO(LBG) = 2.*XCO(LBW) - XCO(LB2)
         YCO(LBG) = 2.*YCO(LBW) - YCO(LB2)
         ZCO(LBG) = 2.*ZCO(LBW) - ZCO(LB2)
         
         XCO(LTG) = 2.*XCO(LTW) - XCO(LT2)
         YCO(LTG) = 2.*YCO(LTW) - YCO(LT2)
         ZCO(LTG) = 2.*ZCO(LTW) - ZCO(LT2)
 100  CONTINUE

C     ****************** REFLEC BLOCK EDGES *************************
C ... G ---- GHOST POINT     1 --- WALL POINT    2 ---- SECOND POINT IN DOMAIN
C ... LL --- LOWER-LOWER     LU -- LOWER-UPPER etc.
C ...               ... g g g G   in figure: g - normal ghost point
C ...               ... x x 1 g              x - normal domain point
C ...               ... x 2 x g              G - calculated ghost point
C ...               ... x x x g              1 - first point in domain
C ...               ... . . . .              2 - second point in domain
C ...               ... . . . . 

      DO 1000 K = 1,KMAXP
         KK     = (K-1+KN)*KSTR
         LLLG   = 1+(0      -1+IN)*ISTR+(0      -1+JN)*JSTR + KK
         LLL1   = 1+(1      -1+IN)*ISTR+(1      -1+JN)*JSTR + KK
         LLL2   = 1+(2      -1+IN)*ISTR+(2      -1+JN)*JSTR + KK

         LLUG   = 1+(0      -1+IN)*ISTR+(JMAXP+1-1+JN)*JSTR + KK
         LLU1   = 1+(1      -1+IN)*ISTR+(JMAXP  -1+JN)*JSTR + KK
         LLU2   = 1+(2      -1+IN)*ISTR+(JMAXP-1-1+JN)*JSTR + KK

         LULG   = 1+(IMAXP+1-1+IN)*ISTR+(0      -1+JN)*JSTR + KK
         LUL1   = 1+(IMAXP  -1+IN)*ISTR+(1      -1+JN)*JSTR + KK
         LUL2   = 1+(IMAXP-1-1+IN)*ISTR+(2      -1+JN)*JSTR + KK

         LUUG   = 1+(IMAXP+1-1+IN)*ISTR+(JMAXP+1-1+JN)*JSTR + KK
         LUU1   = 1+(IMAXP  -1+IN)*ISTR+(JMAXP  -1+JN)*JSTR + KK
         LUU2   = 1+(IMAXP-1-1+IN)*ISTR+(JMAXP-1-1+JN)*JSTR + KK

         XCO(LLLG) = 2.*XCO(LLL1) - XCO(LLL2)
         YCO(LLLG) = 2.*YCO(LLL1) - YCO(LLL2)
         ZCO(LLLG) = 2.*ZCO(LLL1) - ZCO(LLL2)

         XCO(LLUG) = 2.*XCO(LLU1) - XCO(LLU2)
         YCO(LLUG) = 2.*YCO(LLU1) - YCO(LLU2)
         ZCO(LLUG) = 2.*ZCO(LLU1) - ZCO(LLU2)

         XCO(LULG) = 2.*XCO(LUL1) - XCO(LUL2)
         YCO(LULG) = 2.*YCO(LUL1) - YCO(LUL2)
         ZCO(LULG) = 2.*ZCO(LUL1) - ZCO(LUL2)

         XCO(LUUG) = 2.*XCO(LUU1) - XCO(LUU2)
         YCO(LUUG) = 2.*YCO(LUU1) - YCO(LUU2)
         ZCO(LUUG) = 2.*ZCO(LUU1) - ZCO(LUU2)

 1000 CONTINUE

      RETURN
      END SUBROUTINE REFXYZ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GRID1B(N,A1X,A1Y,A1Z,A1,A2X,A2Y,A2Z,A2,
     &     A3X,A3Y,A3Z,A3,V,XK,YK,ZK,UROT,VROT,WROT,XCO,YCO,ZCO,
c    &     XFIC,XFJC,XFKC,YFIC,YFJC,YFKC,ZFIC,ZFJC,ZFKC, ! face centers
     &     IMAX,JMAX,KMAX,IG,MGRID,LFINGD,ICON,NPATCH,GRILEN,
     &     VOLN,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,XCOR,NEGVOL,PRN,
     &     IREPEA,NREPEA,XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,IGRID,IGRID2,
     &     DT,DTOLD,TIMEL,XFC,YFC,ZFC)

      USE NS3CO, ONLY : IN, JN, KN, IC9

      IMPLICIT NONE

      INTEGER :: I1,I2,J1,J2,IFACE,IWALL,IGRID,NTOT1,MNDSTP,
     &     INDSTP,JNDSTP,KNDSTP,IMAXP1,JMAXP1,KMAXP1,ISTRID,JSTRID,
     &     KSTRID,ILE,ISTEP,JSTEP,KSTEP,ISTR,JSTR,KSTR,I,J,K,IT1,
     &     IT0,LEVE,MSTP,KMSTP,NTOT,IT,NTOT2,L,IMH1,JMH1,KMH1,ISH,
     &     JSH,KSH,IMT1,JMT1,KMT1,IST,JST,KST,IG1,IG2,LEVEL,N,
     &     IMAX,JMAX,KMAX,MGRID,LFINGD,NPATCH,IGRID2

      INTEGER :: IG(*), ICON(IC9,*), IREPEA(*), NREPEA(*)

      REAL :: GRILEN,VOLN,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,DT,DTOLD

      REAL :: XCO(*), YCO(*), ZCO(*), XK(*), YK(*), ZK(*),
     &        XLE2(*), YLE2(*), ZLE2(*), XLE3(*), YLE3(*), ZLE3(*)

      REAL :: UROT(*),VROT(*),WROT(*),V(*),
     &        A1X(*),A1Y(*),A1Z(*),A1(*),A2X(*),A2Y(*),A2Z(*),A2(*),
     &        A3X(*),A3Y(*),A3Z(*),A3(*),XCOR(6),XFC(*),YFC(*),ZFC(*)

      REAL, ALLOCATABLE :: VOL1(:),VOL2(:)

      LOGICAL :: NEGVOL, TIMEL
      CHARACTER(LEN=3) :: PRN

      MNDSTP = 2**(LFINGD-1)
      INDSTP = MNDSTP
      JNDSTP = MNDSTP
      KNDSTP = MNDSTP

      IF(IMAX == 2) INDSTP = 1
      IF(JMAX == 2) JNDSTP = 1
      IF(KMAX <= 2) KNDSTP = 1

      IMAXP1 = IMAX + 1
      JMAXP1 = JMAX + 1
      KMAXP1 = KMAX + 1

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN

      ILE    = ISTRID*JSTRID
      NTOT1  = ISTRID*JSTRID*KSTRID

      ISTEP  = 1
      JSTEP  = IMAXP1
      KSTEP  = IMAXP1*JMAXP1

      ISTR   = 1
      JSTR   = IMAXP1-1+2*IN
      KSTR   = (IMAXP1-1+2*IN)*(JMAXP1-1+2*JN)

C ... Build the smallest cartesian cube around a block.
C ... Could be accelerated by scaning only block faces, ESa 3.2.98.

      I = 1
      J = 1
      K = 1
      IT1 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*ILE
      XCOR(1) = XCO(IT1)
      XCOR(2) = XCO(IT1)
      XCOR(3) = YCO(IT1)
      XCOR(4) = YCO(IT1)
      XCOR(5) = ZCO(IT1)
      XCOR(6) = ZCO(IT1)
      DO K = 1,KMAXP1
      DO J = 1,JMAXP1
      DO I = 1,IMAXP1
         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*ILE
         XCOR(1) = MIN(XCOR(1),XCO(IT0))
         XCOR(2) = MAX(XCOR(2),XCO(IT0))
         XCOR(3) = MIN(XCOR(3),YCO(IT0))
         XCOR(4) = MAX(XCOR(4),YCO(IT0))
         XCOR(5) = MIN(XCOR(5),ZCO(IT0))
         XCOR(6) = MAX(XCOR(6),ZCO(IT0))
      ENDDO
      ENDDO
      ENDDO

      DO 110 LEVEL = 1,MGRID
         MSTP  = 2**(LEVEL-1)
         KMSTP = MSTP                                    !PPR
         IF ((KMAXP1-1)/KMSTP <= 1) KMSTP =  (KMAXP1-1)!PPR
C ... Tämä toimii väärin. Merkitys??
C         IF (KMAX <= 2) KMSTP = 1                        !PPR
         NTOT = ((IMAXP1-1)/MSTP+2*IN)*
     &          ((JMAXP1-1)/MSTP+2*JN)*
     &          ((KMAXP1-1)/KMSTP+2*KN)
         IT = IG(LEVEL)-IG(1)
         DO 100  I = 1,NTOT
            A1X(I+IT) = 1.
            A1Y(I+IT) = 0.
            A1Z(I+IT) = 0.
            A2X(I+IT) = 0.
            A2Y(I+IT) = 1.
            A2Z(I+IT) = 0.
            A3X(I+IT) = 0.
            A3Y(I+IT) = 0.
            A3Z(I+IT) = 1.
 100     CONTINUE
 110  CONTINUE

C ... Compute The Properties Of Cell Faces Related To Each Index
C ... Direction (1,2,3) At The Finest Grid Level.

      CALL AREVC1(XCO,YCO,ZCO,IMAXP1+1,JMAXP1,KMAXP1,A1X,A1Y,A1Z,
     +     A1,UROT,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,1,XK,YK,ZK,V,
     +     ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,ILE)
      CALL AREVC1(XCO,YCO,ZCO,IMAXP1,JMAXP1+1,KMAXP1,A2X,A2Y,A2Z,
     +     A2,VROT,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,2,XK,YK,ZK,V,
     +     JSTR,KSTR,ISTR,IN,JN,KN,ISTRID,ILE)
      CALL AREVC1(XCO,YCO,ZCO,IMAXP1,JMAXP1,KMAXP1+1,A3X,A3Y,A3Z,
     +     A3,WROT,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,3,XK,YK,ZK,V,
     +     KSTR,ISTR,JSTR,IN,JN,KN,ISTRID,ILE)

C ... Compute the cell face velocities from the swept volumes

      IF (TIMEL) THEN

         IF(IGRID >= 11 .AND. IGRID <= 20 .OR.
     +      IGRID >= 23 .AND. IGRID <= 39 .OR. 
     +      IGRID2 /= 0)                  THEN ! A moving object

            ALLOCATE (VOL1(NTOT1),VOL2(NTOT1))

            CALL AREVC2(XCO,YCO,ZCO,XLE2,YLE2,ZLE2,IMAXP1+1,JMAXP1,
     +           KMAXP1,1,VOL1,ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,ILE)
            CALL AREVC2(XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,IMAXP1+1,JMAXP1,
     +           KMAXP1,1,VOL2,ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,ILE)
            CALL FACEV(UROT,VOL1,VOL2,A1,DT,DTOLD,IMAXP1+1,JMAXP1,
     +           KMAXP1,1,ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,ILE)
                      
            CALL AREVC2(XCO,YCO,ZCO,XLE2,YLE2,ZLE2,IMAXP1,JMAXP1+1,
     +           KMAXP1,2,VOL1,JSTR,KSTR,ISTR,IN,JN,KN,ISTRID,ILE)
            CALL AREVC2(XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,IMAXP1,JMAXP1+1,
     +           KMAXP1,2,VOL2,JSTR,KSTR,ISTR,IN,JN,KN,ISTRID,ILE)
            CALL FACEV(VROT,VOL1,VOL2,A2,DT,DTOLD,IMAXP1,JMAXP1+1,
     +           KMAXP1,2,JSTR,KSTR,ISTR,IN,JN,KN,ISTRID,ILE)
                  
            CALL AREVC2(XCO,YCO,ZCO,XLE2,YLE2,ZLE2,IMAXP1,JMAXP1,
     +           KMAXP1+1,3,VOL1,KSTR,ISTR,JSTR,IN,JN,KN,ISTRID,ILE)
            CALL AREVC2(XLE2,YLE2,ZLE2,XLE3,YLE3,ZLE3,IMAXP1,JMAXP1,
     +           KMAXP1+1,3,VOL2,KSTR,ISTR,JSTR,IN,JN,KN,ISTRID,ILE)
            CALL FACEV(WROT,VOL1,VOL2,A3,DT,DTOLD,IMAXP1,JMAXP1,
     +           KMAXP1+1,3,KSTR,ISTR,JSTR,IN,JN,KN,ISTRID,ILE)
            
            DEALLOCATE (VOL1,VOL2)

         ENDIF ! IF (IGRID >= 11 .AND.

      ENDIF  ! IF (TIMEL) THEN

C     Clean the singularity boundaries on the finest grid level

      DO 200 L = 1,NPATCH
         IF(ICON(1,L) == 7) THEN
            IFACE = ICON(3,L)
            I1    = ICON(4,L)
            I2    = ICON(5,L)
            J1    = ICON(6,L)
            J2    = ICON(7,L)
            IF (IFACE >= 1 .AND. IFACE <= 3) THEN
               IWALL = 1
            ELSE
               IF (IFACE == 4) IWALL=IMAXP1
               IF (IFACE == 5) IWALL=JMAXP1
               IF (IFACE == 6) IWALL=KMAXP1
            ENDIF

            IF (IFACE == 1 .OR. IFACE == 4) THEN
               CALL AXISW(A1X,A1Y,A1Z,A1,UROT,IWALL,IFACE,
     &              I1,I2,J1,J2,ISTR,JSTR,KSTR,IN,JN,KN,0)
            ENDIF
            IF (IFACE == 2 .OR. IFACE == 5) THEN
               CALL AXISW(A2X,A2Y,A2Z,A2,VROT,IWALL,IFACE,
     &              I1,I2,J1,J2,JSTR,ISTR,KSTR,JN,IN,KN,0)
            ENDIF
            IF (IFACE == 3 .OR. IFACE == 6) THEN
               CALL AXISW(A3X,A3Y,A3Z,A3,WROT,IWALL,IFACE,
     &              I1,I2,J1,J2,KSTR,ISTR,JSTR,KN,IN,JN,0)
            ENDIF
         ENDIF
 200  CONTINUE
C     End of Patch-loop

C ... Multigrid Loop Begins.
      IF (MGRID >= 2) THEN
         DO 350 LEVEL = 2,MGRID
            MSTP  = 2**(LEVEL-1)

C     Max+1 On The Coarse Grid Level
            IMH1 = (IMAXP1-1)/MSTP + 1
            JMH1 = (JMAXP1-1)/MSTP + 1
            KMH1 = (KMAXP1-1)/MSTP + 1
            IF ((KMAXP1-1)/MSTP <= 1) KMH1 = 2             !PPR

C     Increments On The Coarse Grid Level
            ISH  = 1
            JSH  = IMH1-1+2*IN
            KSH  = (IMH1-1+2*IN)*(JMH1-1+2*JN)

C     Max+1 On The Dense Grid Level
            IMT1 = (IMAXP1-1)/(MSTP/2) + 1
            JMT1 = (JMAXP1-1)/(MSTP/2) + 1
            KMT1 = (KMAXP1-1)/(MSTP/2) + 1
            IF ((KMAXP1-1)/MSTP <= 1) KMT1 = 3  !PPR
C     Increments On The Dense Grid Level
            IST  = 1
            JST  = IMT1-1+2*IN
            KST  = (IMT1-1+2*IN)*(JMT1-1+2*JN)

C ... Compute The Properties At The Coarse Grid Levels By Weighted
C ... Averaging And Lumping.

            IG1 = IG(LEVEL-1)-IG(1)
            IG2 = IG(LEVEL)-IG(1)
            CALL LUMPA1(A1X,A1Y,A1Z,A1,UROT,IMH1,JMH1,KMH1,ISH,JSH,KSH,
     &           IMT1,JMT1,KMT1,IST,JST,KST,IN,JN,KN,1,IG1,IG2)
            CALL LUMPA1(A2X,A2Y,A2Z,A2,VROT,JMH1,KMH1,IMH1,JSH,KSH,ISH,
     &           JMT1,KMT1,IMT1,JST,KST,IST,JN,KN,IN,2,IG1,IG2)
            CALL LUMPA1(A3X,A3Y,A3Z,A3,WROT,KMH1,IMH1,JMH1,KSH,ISH,JSH,
     &           KMT1,IMT1,JMT1,KST,IST,JST,KN,IN,JN,3,IG1,IG2)

C     Clean the singularity boundaries on the coarser grid levels
         
            DO 300 L = 1,NPATCH
               IF(ICON(1,L) == 7) THEN
                  IFACE = ICON(3,L)
                  I1    = (ICON(4,L)-1)/MSTP + 1
                  I2    = ICON(5,L)/MSTP
                  J1    = (ICON(6,L)-1)/MSTP + 1
                  J2    = ICON(7,L)/MSTP
                  IF (IFACE >= 1 .AND. IFACE <= 3) THEN
                     IWALL = 1
                  ELSE
                     IF (IFACE == 4) IWALL = IMH1
                     IF (IFACE == 5) IWALL = JMH1
                     IF (IFACE == 6) IWALL = KMH1
                  ENDIF

                  IF (IFACE == 1 .OR. IFACE == 4) THEN
                     CALL AXISW(A1X,A1Y,A1Z,A1,UROT,IWALL,IFACE,
     &                    I1,I2,J1,J2,ISH,JSH,KSH,IN,JN,KN,IG2)
                  ENDIF
                  IF (IFACE == 2 .OR. IFACE == 5) THEN
                     CALL AXISW(A2X,A2Y,A2Z,A2,VROT,IWALL,IFACE,
     &                    I1,I2,J1,J2,JSH,ISH,KSH,JN,IN,KN,IG2)
                  ENDIF
                  IF (IFACE == 3 .OR. IFACE == 6) THEN
                     CALL AXISW(A3X,A3Y,A3Z,A3,WROT,IWALL,IFACE,
     &                    I1,I2,J1,J2,KSH,ISH,JSH,KN,IN,JN,IG2)
                  ENDIF
               ENDIF
 300        CONTINUE
C     End of Patch-loop

 350     CONTINUE
      ENDIF


C ... once upon time this was needed. Nowadays the meening of this has
C ... been forget. PPR 1.9.99
C      IT = IMAXP1+IN+(JMAXP1-1+JN)*(IMAXP1-1+2*IN)+
C     &     (KMAXP1-1+IN)*(IMAXP1-1+2*JN)*(JMAXP1-1+2*JN)
C      V(IT)=1.0
C      XK(IT)=1.0
C      YK(IT)=1.0
C      ZK(IT)=1.0

C ... Check The Validity Of The Grid At The Finest Grid Level.
      CALL CHECKV(V,IMAX,JMAX,KMAX,N,VOLN,NEGVOL,PRN,IREPEA,NREPEA)

C ... Cell Volume Multigrid Loop

      IF(IREPEA(3) <= NREPEA(3)) THEN
      WRITE(45,*)
      WRITE(45,*) 'In GRID1B:'
      WRITE(45,*) '----------'
      ENDIF

      DO 400 LEVEL = 1,MGRID
         MSTP = 2**(LEVEL-1)
         KMSTP = MSTP                                    !PPR
         IF ((KMAXP1-1)/MSTP  <= 1) KMSTP = (KMAXP1-1) !PPR
C ... Tämä toimii väärin. Merkitys??
C         IF (KMAX == 2) KMSTP = 1                        !PPR
         IF(IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,9877) (IMAXP1-1)/MSTP,(JMAXP1-1)/MSTP,
     &        (KMAXP1-1)/KMSTP,LEVEL
         ENDIF

C ... Compute The Volumes And Cell Centerpoint Coordinates On The
C ... Coarse Grid Levels.

         IF(LEVEL >= 2) THEN
            IG1 = IG(LEVEL-1)-IG(1)
            IG2 = IG(LEVEL)-IG(1)
            CALL LUMPVS(V,XK,YK,ZK,IMAXP1,JMAXP1,KMAXP1,MSTP,IG1,IG2)
         ENDIF
 400  CONTINUE

 9877 FORMAT('IMAX = ',I4,' JMAX = ',I4,' KMAX = ',I4,' LEVEL = ',I4)

      RETURN
      END SUBROUTINE GRID1B
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AREVC1(X,Y,Z,IMAXP1,JMAXP1,KMAXP1,AX,AY,
     +     AZ,A,ADR,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,IFACE,XK,YK,ZK,V,
     +     ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,IL)

C ... Cell Face Geometrical Properties, Cell Volumes And Cell
C ... Centerpoint Coordinates Are Defined.

      USE NS3CO, ONLY : NEGVL, NEGV

      REAL :: X(*), Y(*), Z(*), XK(*), YK(*), ZK(*)
      REAL :: AX(*), AY(*), AZ(*), A(*), ADR(*), V(*)
      REAL :: VDIV(3), VSUM(3)

      REAL :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,YA,ZA,XA,
     +     AXX,AYY,AZZ,AA,RDIAGX,RDIAGY,RDIAGZ,RD2AGX,RD2AGY,RD2AGZ,
     +     AXX2,AYY2,AZZ2,AA2,X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,
     +     XCE,YCE,ZCE,OMX,OMY,OMZ
             
      VDIV(1) = 1.0
      VDIV(2) = 1.0
      VDIV(3) = 3.0
      VSUM(1) = 0.0
      VSUM(2) = 1.0
      VSUM(3) = 1.0

C ... ROTATIONAL AXIS AND THE PLACE

      OMX = OMEX
      OMY = OMEY
      OMZ = OMEZ

C ... Cell Face A1 (I-Direction):

      IMP1 = IMAXP1
      JMP1 = JMAXP1
      KMP1 = KMAXP1

      IF(IFACE == 1) IMP1 = IMAXP1 - 1
      IF(IFACE == 2) JMP1 = JMAXP1 - 1
      IF(IFACE == 3) KMP1 = KMAXP1 - 1

      DO 100 K = 0,KMP1
      DO 110 J = 0,JMP1
      DO 120 I = 0,IMP1

         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*IL

         IT1 = IT0
         IT2 = IT0               + KSTR
         IT3 = IT0        + JSTR + KSTR
         IT4 = IT0        + JSTR

         IT5 = IT0 + ISTR
         IT6 = IT0 + ISTR        + KSTR
         IT7 = IT0 + ISTR + JSTR + KSTR
         IT8 = IT0 + ISTR + JSTR

C ... Corner Points Coordinates:
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

C ... Diagonal Components:
         RDIAGX = XCE-.5*(X1 + X2)
         RDIAGY = YCE-.5*(Y1 + Y2)
         RDIAGZ = ZCE-.5*(Z1 + Z2)
         RD2AGX = XCE-.5*(X7 + X8)
         RD2AGY = YCE-.5*(Y7 + Y8)
         RD2AGZ = ZCE-.5*(Z7 + Z8)
C ... Cell Face Area Components:
         AXX = 0.5*((Z3-Z1)*(Y4-Y2)-(Y3-Y1)*(Z4-Z2))
         AYY = 0.5*((X3-X1)*(Z4-Z2)-(Z3-Z1)*(X4-X2))
         AZZ = 0.5*((Y3-Y1)*(X4-X2)-(X3-X1)*(Y4-Y2))
         AA  = SQRT(AXX**2+AYY**2+AZZ**2)

         AXX2 = -0.5*((Z7-Z5)*(Y8-Y6)-(Y7-Y5)*(Z8-Z6))
         AYY2 = -0.5*((X7-X5)*(Z8-Z6)-(Z7-Z5)*(X8-X6))
         AZZ2 = -0.5*((Y7-Y5)*(X8-X6)-(X7-X5)*(Y8-Y6))
         AA2  = SQRT(AXX2**2+AYY2**2+AZZ2**2)
C ... Normal Components:
         AX(IT0) = AXX/(AA+1.0E-20)
         AY(IT0) = AYY/(AA+1.0E-20)
         AZ(IT0) = AZZ/(AA+1.0E-20)
         AX(IT5) = -AXX2/(AA2+1.0E-20)
         AY(IT5) = -AYY2/(AA2+1.0E-20)
         AZ(IT5) = -AZZ2/(AA2+1.0E-20)
C ... Total Area:
         A(IT0) = AA
         A(IT5) = AA2

         XK(IT0) = XCE
         YK(IT0) = YCE
         ZK(IT0) = ZCE
         
C ... Volume Increment:
         V(IT0) = (VSUM(IFACE)*V(IT0) + 
     &        RDIAGX*AXX + RDIAGY*AYY +RDIAGZ*AZZ + 
     &        RD2AGX*AXX2+ RD2AGY*AYY2+RD2AGZ*AZZ2)/VDIV(IFACE)

ccccc      IF(V(IT0) < 1.E-11 .AND. V(IT0) > -1.E-9) V(IT0) = 1.E-11 
         IF (NEGVL) THEN
            IF(V(IT0) < 1.E-11 .AND. V(IT0) >= NEGV) V(IT0) = 1.E-11
         ENDIF

 120  CONTINUE
 110  CONTINUE
 100  CONTINUE

      DO K = 0,KMAXP1
      DO J = 0,JMAXP1
      DO I = 0,IMAXP1

         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*IL

         IF(A(IT0) < 1.E-20) THEN
            AX(IT0) = 1.
            AY(IT0) = 0.
            AZ(IT0) = 0.
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      DO K = 0,KMAXP1
      DO J = 0,JMAXP1
      DO I = 0,IMAXP1

         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*IL

         IT1 = IT0
         IT2 = IT0               + KSTR
         IT3 = IT0        + JSTR + KSTR
         IT4 = IT0        + JSTR

         IT5 = IT0 + ISTR

C ... Corner Points Coordinates:

         X1 = X(IT1)-CENX
         Y1 = Y(IT1)-CENY
         Z1 = Z(IT1)-CENZ
         X2 = X(IT2)-CENX
         Y2 = Y(IT2)-CENY
         Z2 = Z(IT2)-CENZ
         X3 = X(IT3)-CENX
         Y3 = Y(IT3)-CENY
         Z3 = Z(IT3)-CENZ
         X4 = X(IT4)-CENX
         Y4 = Y(IT4)-CENY
         Z4 = Z(IT4)-CENZ

         AA = A(IT0)

C ... Static Moment Around The X-Axis (For Rotating Coordinates):
         YA = ((Y1+Y2)*(Y1*X2-Y2*X1)+(Y2+Y3)*(Y2*X3-Y3*X2)
     2        +(Y3+Y4)*(Y3*X4-Y4*X3)+(Y4+Y1)*(Y4*X1-Y1*X4))
         ZA = ((Z1+Z2)*(Z2*X1-Z1*X2)+(Z2+Z3)*(Z3*X2-Z2*X3)
     2        +(Z3+Z4)*(Z4*X3-Z3*X4)+(Z4+Z1)*(Z1*X4-Z4*X1))
         AXX= (YA-ZA)

C ... Static Moment Around The Y-Axis (For Rotating Coordinates):
         ZA = ((Z1+Z2)*(Z1*Y2-Z2*Y1)+(Z2+Z3)*(Z2*Y3-Z3*Y2)
     2        +(Z3+Z4)*(Z3*Y4-Z4*Y3)+(Z4+Z1)*(Z4*Y1-Z1*Y4))
         XA = ((X1+X2)*(X2*Y1-X1*Y2)+(X2+X3)*(X3*Y2-X2*Y3)
     2        +(X3+X4)*(X4*Y3-X3*Y4)+(X4+X1)*(X1*Y4-X4*Y1))
         AYY= (ZA-XA)

C ... Static Moment Around The Z-Axis (For Rotating Coordinates):
         XA = ((X1+X2)*(X1*Z2-X2*Z1)+(X2+X3)*(X2*Z3-X3*Z2)
     2        +(X3+X4)*(X3*Z4-X4*Z3)+(X4+X1)*(X4*Z1-X1*Z4))
         YA = ((Y1+Y2)*(Y2*Z1-Y1*Z2)+(Y2+Y3)*(Y3*Z2-Y2*Z3)
     2        +(Y3+Y4)*(Y4*Z3-Y3*Z4)+(Y4+Y1)*(Y1*Z4-Y4*Z1))
         AZZ= (XA-YA)

         ADR(IT0) = (OMX*AXX + OMY*AYY + OMZ*AZZ)/
     +        (6.*AA+1.0E-20)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE AREVC1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AREVC2(X,Y,Z,XLE2,YLE2,ZLE2,IMAXP1,JMAXP1,KMAXP1,
     +     IFACE,V,ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,IL)

C ... Cell Face velocity is based on the volume swept byt the
C ... surface

      IMPLICIT NONE

      INTEGER :: IMAXP1,JMAXP1,KMAXP1,IFACE,ISTR,JSTR,KSTR,IN,JN,KN,
     +           ISTRID,IL,IMP1,JMP1,KMP1,I,J,K,IT0,IT1,IT2,IT3,IT4

      REAL :: X(*), Y(*), Z(*), XLE2(*), YLE2(*), ZLE2(*)

      REAL :: V(*)

      REAL :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                    X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,VIT0

C ... Cell Face A1 (I-Direction):

      IMP1 = IMAXP1
      JMP1 = JMAXP1
      KMP1 = KMAXP1

      IF(IFACE == 1) IMP1 = IMAXP1 - 1
      IF(IFACE == 2) JMP1 = JMAXP1 - 1
      IF(IFACE == 3) KMP1 = KMAXP1 - 1

      DO 100 K = 0,KMP1
      DO 110 J = 0,JMP1
      DO 120 I = 0,IMP1

         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*IL

         IT1 = IT0
         IT2 = IT0               + KSTR
         IT3 = IT0        + JSTR + KSTR
         IT4 = IT0        + JSTR

C ... Corner Points Coordinates:
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

         X5 = XLE2(IT1)
         Y5 = YLE2(IT1)
         Z5 = ZLE2(IT1)
         X6 = XLE2(IT2)
         Y6 = YLE2(IT2)
         Z6 = ZLE2(IT2)
         X7 = XLE2(IT3)
         Y7 = YLE2(IT3)
         Z7 = ZLE2(IT3)
         X8 = XLE2(IT4)
         Y8 = YLE2(IT4)
         Z8 = ZLE2(IT4)

      CALL VOLINC(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &            X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,VIT0,1)

      CALL VOLINC(X1,Y1,Z1,X5,Y5,Z5,X6,Y6,Z6,X2,Y2,Z2,
     &            X4,Y4,Z4,X8,Y8,Z8,X7,Y7,Z7,X3,Y3,Z3,VIT0,2)

      CALL VOLINC(X5,Y5,Z5,X1,Y1,Z1,X4,Y4,Z4,X8,Y8,Z8,
     &            X6,Y6,Z6,X2,Y2,Z2,X3,Y3,Z3,X7,Y7,Z7,VIT0,3)

      V(IT0) = VIT0

 120  CONTINUE
 110  CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE AREVC2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE VOLINC(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &                  X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,VITO,IFACE)

C ... Volume increment.

      IMPLICIT NONE

      INTEGER :: IFACE

      REAL :: VITO, VDIV(3),VSUM(3)

      REAL :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +     AXX,AYY,AZZ,AA,RDIAGX,RDIAGY,RDIAGZ,RD2AGX,RD2AGY,RD2AGZ,
     +     AXX2,AYY2,AZZ2,AA2,X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,
     +     XCE,YCE,ZCE

      VDIV(1) = 1.0
      VDIV(2) = 1.0
      VDIV(3) = 3.0
      VSUM(1) = 0.0
      VSUM(2) = 1.0
      VSUM(3) = 1.0

      XCE = .125*(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8)
      YCE = .125*(Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8)
      ZCE = .125*(Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8)

C ... Diagonal Components:

      RDIAGX = XCE-.5*(X1 + X2)
      RDIAGY = YCE-.5*(Y1 + Y2)
      RDIAGZ = ZCE-.5*(Z1 + Z2)
      RD2AGX = XCE-.5*(X7 + X8)
      RD2AGY = YCE-.5*(Y7 + Y8)
      RD2AGZ = ZCE-.5*(Z7 + Z8)

C ... Cell Face Area Components:

      AXX = 0.5*((Z3-Z1)*(Y4-Y2)-(Y3-Y1)*(Z4-Z2))
      AYY = 0.5*((X3-X1)*(Z4-Z2)-(Z3-Z1)*(X4-X2))
      AZZ = 0.5*((Y3-Y1)*(X4-X2)-(X3-X1)*(Y4-Y2))

      AXX2 = -0.5*((Z7-Z5)*(Y8-Y6)-(Y7-Y5)*(Z8-Z6))
      AYY2 = -0.5*((X7-X5)*(Z8-Z6)-(Z7-Z5)*(X8-X6))
      AZZ2 = -0.5*((Y7-Y5)*(X8-X6)-(X7-X5)*(Y8-Y6))
         
C ... Volume Increment:
      VITO = (VSUM(IFACE)*VITO + 
     &        RDIAGX*AXX + RDIAGY*AYY +RDIAGZ*AZZ + 
     &        RD2AGX*AXX2+ RD2AGY*AYY2+RD2AGZ*AZZ2)/VDIV(IFACE)
           
      RETURN
      END SUBROUTINE VOLINC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FACEV(UROT,VOL1,VOL2,A1,DT,DTOLD,IMAXP1,JMAXP1,KMAXP1,
     +     IFACE,ISTR,JSTR,KSTR,IN,JN,KN,ISTRID,IL)

C ... Cell Face velocities form the swept volumes

      IMPLICIT NONE

      REAL :: UROT(*),VOL1(*),VOL2(*),A1(*)
      REAL :: DT,DTOLD
      INTEGER :: IMAXP1,JMAXP1,KMAXP1,IMP1,JMP1,KMP1,I,J,K,IT0,IT5,ISTR,
     +           IFACE,ISTRID,JSTR,KSTR,IN,JN,KN,IL

C ... Cell Face A1 (I-Direction):

      IMP1 = IMAXP1
      JMP1 = JMAXP1
      KMP1 = KMAXP1
      IF(IFACE == 1) IMP1 = IMAXP1 - 1
      IF(IFACE == 2) JMP1 = JMAXP1 - 1
      IF(IFACE == 3) KMP1 = KMAXP1 - 1

      DO 100 K = 0,KMP1
      DO 110 J = 0,JMP1
      DO 120 I = 0,IMP1
         IT0 = 1+(I-1  +IN) +(J-1  +JN)*ISTRID+(K-1  +KN)*IL
         IT5 = IT0 + ISTR
              
         IF(ABS(VOL2(IT0)/(A1(IT0)+1.E-20)) <= 1.E-6) THEN
         UROT(IT0) =-VOL1(IT0) / (DT * (A1(IT0)+1.E-20))
            ELSE
         UROT(IT0) =-(1.5*VOL1(IT0)/DT -.5*VOL2(IT0)/DTOLD) /
     +               ((A1(IT0)+1.E-20))
         ENDIF
         
c       if(k == 3 .and. J == 5 .and. i == 10)
c      if(iface == 1)
c     +  write(666,*)i,j,k, urot(it0),vol1(it0),vol2(it0),a1(it0)
c        UROT(IT5) = (1.5*VOL1(IT5) -.5*VOL2(IT5))/((A1(IT5)+1.E-20)*DT)! Mersu
            
120   CONTINUE
110   CONTINUE
100   CONTINUE

      RETURN
      END SUBROUTINE FACEV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AXISW(AX,AY,AZ,A,ADR,IWALL,IFACE,J1,J2,K1,K2,
     &                 ISTR,JSTR,KSTR,IN,JN,KN,IGL)

C ... Reset the grid properties on singularity patches.

      REAL    :: AX(*), AY(*), AZ(*), A(*), ADR(*)
      INTEGER :: IWALL,IFACE,J1,J2,K1,K2,I0,I1,I2,IT0,IT1,IT2
      INTEGER :: IN, JN, KN

C ... Cell Face A1 (I-Direction):

      IF (IFACE == 1 .OR. IFACE == 2 .OR. IFACE == 3) THEN
         I0 = IWALL
         I1 = IWALL+1
         I2 = IWALL+2
      ELSE
         I0 = IWALL
         I1 = IWALL-1
         I2 = IWALL-2
      ENDIF

      DO 100 K = K1,K2
         DO 110 J = J1,J2
            IT0 = IGL+1+(I0-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
            IT1 = IGL+1+(I1-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR
            IT2 = IGL+1+(I2-1+IN)*ISTR+(J-1+JN)*JSTR+(K-1+KN)*KSTR

            A(IT0)   = 0.0
            ADR(IT0) = 0.0
            AXX      = 2.0*AX(IT1)-1.0*AX(IT2)
            AYY      = 2.0*AY(IT1)-1.0*AY(IT2)
            AZZ      = 2.0*AZ(IT1)-1.0*AZ(IT2)
            AA       = SQRT(AXX**2+AYY**2+AZZ**2)
            AX(IT0)  = AXX/(AA+1.0E-20)
            AY(IT0)  = AYY/(AA+1.0E-20)
            AZ(IT0)  = AZZ/(AA+1.0E-20)
 110     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE AXISW
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CHECKV(V,IMAX,JMAX,KMAX,N,VOLN,NEGVOL,PRN,
     +                  IREPEA,NREPEA)

C ... The Validity Of The Grid Is Checked By Detecting Possible
C ... Negative Volumes. The User Is Warned In Case Of An Invalid Grid.

      USE NS3CO, ONLY : IN, JN, KN
      USE INTEGERS, ONLY : IPRO
      USE MAIN_ARRAYS, ONLY : NPROCE
    
      IMPLICIT NONE

      REAL :: V(*), VOLN
      LOGICAL :: NEGVOL, NEGLOC
      CHARACTER (LEN=3) :: PRN
      INTEGER :: IREPEA(*), NREPEA(*)
      INTEGER :: IMAX, JMAX, KMAX, I, J, K, N, IT
      INTEGER :: ISTRID, JSTRID, LEN, IRR, JRR, KRR
      INTEGER :: LERROR, IUERR, ILERR, JUERR, JLERR, KUERR, KLERR

      NEGLOC = .FALSE.
      LERROR = 0
      IUERR  = 0
      ILERR  = IMAX + 1
      JUERR  = 0
      JLERR  = JMAX + 1
      KUERR  = 0
      KLERR  = KMAX + 1

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      LEN = 0

      DO K = 1,KMAX
         DO J = 1,JMAX
            DO I = 1,IMAX
               IT = I+IN + (J-1+JN)*ISTRID + (K-1+KN)*ISTRID*JSTRID
               IF(V(IT) <= 0.) THEN
                  LERROR = 1
                  IF(I > IUERR) IUERR = I
                  IF(I < ILERR) ILERR = I
                  IF(J > JUERR) JUERR = J
                  IF(J < JLERR) JLERR = J
                  IF(K > KUERR) KUERR = K
                  IF(K < KLERR) KLERR = K
                  IRR  = I
                  JRR  = J
                  KRR  = K
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(LERROR /= 0 .AND. IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*)
         WRITE(45,500) ILERR, IUERR, JLERR, JUERR, KLERR, KUERR
         WRITE(45,600) N, IRR, JRR, KRR
      ENDIF

      LERROR = 0
      IUERR  = 0
      ILERR  = IMAX + 1 
      JUERR  = 0
      JLERR  = JMAX + 1
      KUERR  = 0
      KLERR  = KMAX + 1

      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      VOLN  = 0.

      DO K = 1,KMAX
         DO J = 1,JMAX
            DO I = 1,IMAX
               IT = I+IN + (J-1+JN)*ISTRID + (K-1+KN)*ISTRID*JSTRID
               VOLN = VOLN + V(IT)
               IF(V(IT) <= 0.) THEN
                  LERROR = 1
                  IF(I > IUERR) IUERR = I
                  IF(I < ILERR) ILERR = I
                  IF(J > JUERR) JUERR = J
                  IF(J < JLERR) JLERR = J
                  IF(K > KUERR) KUERR = K
                  IF(K < KLERR) KLERR = K
                  IRR  = I
                  JRR  = J
                  KRR  = K
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      IF(LERROR /= 0 .AND. IREPEA(3) <= NREPEA(3)) THEN
         WRITE(45,*) 
         WRITE(45,501) ILERR, IUERR, JLERR, JUERR, KLERR, KUERR
         WRITE(45,600) N, IRR, JRR, KRR
      ENDIF

 500  FORMAT(1X,'Negative volumes:',1X,'I between ',I3,' to ',I3,
     &     ', J  ',I3,' to ',I3,1X,'and K  ',I3,' to ',I3,'.')
 600  FORMAT(1X,'Last error was found at block ',I3,' in point (i,j,k) '
     &     ,'(',I4,',',I4,',',I4,')')
 501  FORMAT(1X,'Negative cell volumes:',1X,'I between ',I3,' to ',I3,
     &     ', J  ',I3,' to ',I3,1X,'and K  ',I3,' to ',I3,'.')
 601  FORMAT(1X,'Last error was found at block ',I3,' in point (i,j,k) '
     &     ,'(',I4,',',I4,',',I4,')')

      DO K = 1,KMAX
         DO J = 1,JMAX
            DO I = 1,IMAX
               IT = I+IN + (J-1+JN)*ISTRID + (K-1+KN)*ISTRID*JSTRID
               IF(V(IT) <= 0.) THEN
                  NEGVOL = .TRUE.
                  NEGLOC = .TRUE.
               ENDIF
            ENDDO
         ENDDO
      ENDDO
         
      IF(NEGLOC) THEN

         OPEN(47,FILE= 'NEGVOL'//PRN,FORM='FORMATTED',STATUS='UNKNOWN')
         WRITE(47,*)  
*         WRITE(47,*) 'BLOCK =',N
         WRITE(47,*) 'BLOCK =',NPROCE(1+N,IPRO)  ! Global block number
         WRITE(47,*) '*******************'
         WRITE(47,*)    
      
         DO K = 1,KMAX
            DO J = 1,JMAX
               DO I = 1,IMAX
                  IT = I+IN + (J-1+JN)*ISTRID + (K-1+KN)*ISTRID*JSTRID
                  IF(V(IT) <= 0. .AND. LEN < 10000) THEN
                     WRITE(47,*) I, J, K, V(IT)
                     LEN = LEN + 1
                  ENDIF 
               ENDDO
            ENDDO
         ENDDO

         IF(LEN >= 10000) WRITE(47,*) '...continuing'

      ENDIF

      RETURN
      END SUBROUTINE CHECKV
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LUMPA1(AX,AY,AZ,A,ADR,IMH1,JMH1,KMH1,ISH,JSH,KSH,
     &                  IMT1,JMT1,KMT1,IST,JST,KST,
     &                  IN,JN,KN,IFACE,IG1,IG2)

C ... The Cell Face Areas And Their Components Are Computed By
C ... Lumping And Weighted Averaging Of The Results At The Previous
C ... Grid Level.

      REAL :: AX(*), AY(*), AZ(*), A(*), ADR(*)
      REAL :: AXLUMP, AYLUMP, AZLUMP, AXYZL
      INTEGER :: IN, JN, KN


C ... Cell Face A1 (I-Direction):

      JNC = 1                 !PPR
      KNC = 1                 !PPR
      INC = 1                 !PPR
      IF (JMH1 == 2) JNC = 0  !PPR
      IF (KMH1 == 2) KNC = 0  !PPR
      IF (IMH1 == 2) INC = 0  !PPR
      DO 100 K = 1,KMT1-2,2
         DO 110 J = 1,JMT1-2,2
            DO 120 I = 1,IMT1,2
               IT0 = IG2+1+((I-1)/2+IN)*ISH+((J-1)/2+JN)*JSH+
     &              ((K-1)/2+KN)*KSH
C               II = (I-1+IN)*IST
               II = (I-1)*IST*INC+IN*IST
               IT1=IG1+1+II+(J-1    +JN)*JST+(K-1    +KN)*KST
               IT2=IG1+1+II+(J-1+JNC+JN)*JST+(K-1    +KN)*KST
               IT3=IG1+1+II+(J-1    +JN)*JST+(K-1+KNC+KN)*KST
               IT4=IG1+1+II+(J-1+JNC+JN)*JST+(K-1+KNC+KN)*KST
               AXLUMP = A(IT1)*AX(IT1)+A(IT2)*AX(IT2)+
     &                  A(IT3)*AX(IT3)+A(IT4)*AX(IT4)
               AYLUMP = A(IT1)*AY(IT1)+A(IT2)*AY(IT2)+
     &                  A(IT3)*AY(IT3)+A(IT4)*AY(IT4)
               AZLUMP = A(IT1)*AZ(IT1)+A(IT2)*AZ(IT2)+
     &                  A(IT3)*AZ(IT3)+A(IT4)*AZ(IT4)
               ALUMP  = A(IT1)+A(IT2)+A(IT3)+A(IT4) + 1.E-20
c               ADRLMP = ADR(IT1)+ADR(IT2)+ADR(IT3)+ADR(IT4)
               ADRLMP =(ADR(IT1)*A(IT1)+ADR(IT2)*A(IT2)+ADR(IT3)*A(IT3)
     &                 +ADR(IT4)*A(IT4))/ALUMP
               AXYZL  = SQRT(AXLUMP**2+AYLUMP**2+AZLUMP**2)
               AX(IT0) = AXLUMP/(AXYZL+1.E-20)
               AY(IT0) = AYLUMP/(AXYZL+1.E-20)
               AZ(IT0) = AZLUMP/(AXYZL+1.E-20)

C ... Modify Ax To Suit The Flow Solver If AXYZL Is Zero.
               IF(AXYZL < 1.E-20) AX(IT0) = 1.
               A(IT0)   = ALUMP
               ADR(IT0) = ADRLMP
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE LUMPA1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LUMPVS(V,XK,YK,ZK,IMAXP1,JMAXP1,KMAXP1,MSTP,IG1,IG2)

C ... The Volumes At Coarse Grid Levels Are Computed By Lumping
C ... The Values At The Previous Grid Level.

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: V(*)
      REAL :: XK(*), YK(*), ZK(*)

      IMH1 = (IMAXP1-1)/MSTP + 1
      JMH1 = (JMAXP1-1)/MSTP + 1
      KMH1 = (KMAXP1-1)/MSTP + 1

      ISH  = IMH1-1+2*IN
      JSH  = JMH1-1+2*JN
      KSH  = KMH1-1+2*KN

      IMT1 = (IMAXP1-1)/(MSTP/2) + 1
      JMT1 = (JMAXP1-1)/(MSTP/2) + 1
      KMT1 = (KMAXP1-1)/(MSTP/2) + 1
      IF (KMT1 <= 3) KMT1 = 3     !PPR

      IST  = IMT1-1+2*IN
      JST  = JMT1-1+2*JN
      KST  = KMT1-1+2*KN

      INC = 1
      KNC = 1
      IF (KMT1 <=  3) KNC = 0       !PPR (ei tää oo kylläkään käytössä)
      DO 100 K = 1,KMT1-2,2
         DO 110 J = 1,JMT1-2,2
ccc      DO 120 I = 1,IMAXP1-2,2 ! Tämä oli tässä yli 5 vuotta!
         DO 120 I = 1,IMT1-2,2
               IT0 = IG2+1+((I-1)/2+IN)+((J-1)/2+JN)*ISH+
     &              ((K-1)/2+KN)*ISH*JSH
               IT1 = IG1+I+IN+    (J-1    +JN)*IST+(K-1    +KN)*IST*JST
               IT2 = IG1+I+IN+    (J-1+INC+JN)*IST+(K-1    +KN)*IST*JST
               IT3 = IG1+I+IN+    (J-1    +JN)*IST+(K-1+KNC+KN)*IST*JST
               IT4 = IG1+I+IN+    (J-1+INC+JN)*IST+(K-1+KNC+KN)*IST*JST
               IT5 = IG1+I+IN+INC+(J-1    +JN)*IST+(K-1    +KN)*IST*JST
               IT6 = IG1+I+IN+INC+(J-1+INC+JN)*IST+(K-1    +KN)*IST*JST
               IT7 = IG1+I+IN+INC+(J-1    +JN)*IST+(K-1+KNC+KN)*IST*JST
               IT8 = IG1+I+IN+INC+(J-1+INC+JN)*IST+(K-1+KNC+KN)*IST*JST

               VLUMP = V(IT1)+V(IT2)+V(IT3)+V(IT4)+V(IT5)+V(IT6)+V(IT7)+
     &              V(IT8)
               XKLMP = V(IT1)*XK(IT1)+V(IT2)*XK(IT2)+V(IT3)*XK(IT3)+
     &              V(IT4)*XK(IT4)+V(IT5)*XK(IT5)+V(IT6)*XK(IT6)+
     &              V(IT7)*XK(IT7)+V(IT8)*XK(IT8)
               YKLMP = V(IT1)*YK(IT1)+V(IT2)*YK(IT2)+V(IT3)*YK(IT3)+
     &              V(IT4)*YK(IT4)+V(IT5)*YK(IT5)+V(IT6)*YK(IT6)+
     &              V(IT7)*YK(IT7)+V(IT8)*YK(IT8)
               ZKLMP = V(IT1)*ZK(IT1)+V(IT2)*ZK(IT2)+V(IT3)*ZK(IT3)+
     &              V(IT4)*ZK(IT4)+V(IT5)*ZK(IT5)+V(IT6)*ZK(IT6)+
     &              V(IT7)*ZK(IT7)+V(IT8)*ZK(IT8)
               V(IT0)  = VLUMP
               XK(IT0) = XKLMP/VLUMP
               YK(IT0) = YKLMP/VLUMP
               ZK(IT0) = ZKLMP/VLUMP
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      RETURN
      END SUBROUTINE LUMPVS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE D1D2D3(NBL,IMAX,JMAX,KMAX,MGRID,IG,A1,A2,A3,V,
     +                  D1,D2,D3,MGM)

      USE NS3CO, ONLY : IN, JN, KN

      INTEGER :: IMAX(MGM,*),JMAX(MGM,*),KMAX(MGM,*),MGRID(*),IG(MGM,*)
      REAL    :: A1(*), A2(*), A3(*), V(*), D1(*), D2(*), D3(*)

      EPS     = 1.E-20

      DO 8000 MG = 1,MGRID(NBL)
         ISTRID  = IMAX(MG,NBL) + 2*IN
         JSTRID  = JMAX(MG,NBL) + 2*JN
         KSTRID  = KMAX(MG,NBL) + 2*KN
         IL      = ISTRID*JSTRID
         NTOT1   = ISTRID*JSTRID*KSTRID
         IT2     = IG(MG,NBL)-1
C     IT1     = IG(MG,NBL)-IG(1,NBL)

      DO 6000 N = 1,NTOT1
         D3(IT2+N) = 1.E4
         D2(IT2+N) = 1.E4
         D1(IT2+N) = 1.E4
 6000 CONTINUE
C ... LISASIN MINIT KOSKA HAAMUKOPEISSA YLIPITKIA KOPPEJA 2D tapaus

      DO 7000 K = -1,KMAX(MG,NBL)+1
      KA      = (KN+K-1)*IL
      DO 7000 J = -1,JMAX(MG,NBL)+1
      JJ      = (JN+J-1)*ISTRID + IN + KA
      DO 7000 I = -1,IMAX(MG,NBL)+1
      KK      = JJ + I + IT2
         D3(KK) = MIN(
     +     1.E4,(V(KK)+EPS)/(A3(KK)+A3(KK+    IL)+EPS)*2.)
         D2(KK) = MIN(
     +     1.E4,(V(KK)+EPS)/(A2(KK)+A2(KK+ISTRID)+EPS)*2.)
         D1(KK) = MIN(
     +     1.E4,(V(KK)+EPS)/(A1(KK)+A1(KK+     1)+EPS)*2.)
7000  CONTINUE
8000  CONTINUE

      RETURN
      END SUBROUTINE D1D2D3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LUMXCO(X2,Y2,Z2,X1,Y1,Z1,IMAX,JMAX,KMAX,KMAX0)

C ... LUMP CORNER POINTS TO THE LOWER GRID LEVELS

      USE NS3CO, ONLY : IN, JN, KN

      REAL :: X2(*), Y2(*), Z2(*), X1(*), Y1(*), Z1(*)

      ISTR2 = IMAX + 2*IN
      JSTR2 = JMAX + 2*JN
      IL2   = ISTR2*JSTR2
      ISTR1 = 2*IMAX + 2*IN
      JSTR1 = 2*JMAX + 2*JN
      IL1   = ISTR1*JSTR1

      KSTR = 2
      IF(KMAX0 == 1) KSTR = 1

      DO K2 = 1,KMAX+1
      DO J2 = 1,JMAX+1
      DO I2 = 1,IMAX+1
         K1 = (K2-1)*KSTR + 1
         J1 = (J2-1)*2    + 1
         I1 = (I2-1)*2    + 1
         IA2 = 1 + (I2-1+IN) + (J2-1+JN)*ISTR2 + (K2-1+KN)*IL2
         IA1 = 1 + (I1-1+IN) + (J1-1+JN)*ISTR1 + (K1-1+KN)*IL1
         X2(IA2) = X1(IA1)
         Y2(IA2) = Y1(IA1)
         Z2(IA2) = Z1(IA1)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE LUMXCO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE COVEC(CNX,CNY,CNZ,X,Y,Z,NI,NJ,NK,ISTRID,JSTRID,KSTRID,
     +                 NISTR,NJSTR,NKSTR)

      USE NS3CO, ONLY : IN, JN, KN

      IMPLICIT NONE

      REAL :: CNX(*), CNY(*), CNZ(*), RLEN
      REAL :: X(*), Y(*), Z(*)
      INTEGER :: NI,NJ,NK,ISTRID,JSTRID,KSTRID,I,J,K,II,IJ,IM,
     +           IJ2,NISTR,NJSTR,NKSTR

C ... Calculate covariant unit vectors

      DO K = 1,NK
      DO J = 1,NJ

         IJ  = (J+JN-1)*JSTRID + (K+KN-1)*KSTRID
         IJ2 = (J-1)*NJSTR + (K-1)*NKSTR

         DO I = 1,NI

            II  = 1 + (I+IN-1)*ISTRID + IJ
            IM  = 1 + (I-1)*NISTR + IJ2

C ... In i-direction

            IF(I == 1) THEN
               CNX(IM) = X(II+ISTRID)-X(II)
               CNY(IM) = Y(II+ISTRID)-Y(II)
               CNZ(IM) = Z(II+ISTRID)-Z(II)
            ELSE IF(I == NI) THEN
               CNX(IM) = X(II)-X(II-ISTRID)
               CNY(IM) = Y(II)-Y(II-ISTRID)
               CNZ(IM) = Z(II)-Z(II-ISTRID)
            ELSE
               CNX(IM) = .5*(X(II+ISTRID)-X(II-ISTRID))
               CNY(IM) = .5*(Y(II+ISTRID)-Y(II-ISTRID))
               CNZ(IM) = .5*(Z(II+ISTRID)-Z(II-ISTRID))
            ENDIF

            RLEN = SQRT(CNX(IM)**2 + CNY(IM)**2 + CNZ(IM)**2) + 1.E-20

            CNX(IM) = CNX(IM)/RLEN
            CNY(IM) = CNY(IM)/RLEN
            CNZ(IM) = CNZ(IM)/RLEN

         ENDDO

      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE COVEC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C **********************************************************************
C                                                                      *
C     NEW CONNECTION ROUTINES                                          *
C                                                                      *
C **********************************************************************

      SUBROUTINE CORNER(XCO,YCO,ZCO,XCOL,IFACE,IXLQ,IXUQ,IYLQ,IYUQ,
     +                  IMAX,JMAX,KMAX,IN,JN,KN,L)

      REAL :: XCO(*), YCO(*), ZCO(*)
      REAL :: XCOL(3,4)
C
C ... FIND THE CORNER POINTS
C
      IXLO   = IXLQ
      IYLO   = IYLQ
      IXUP   = IXUQ + 1
      IYUP   = IYUQ + 1
      ISTRID = IMAX + 2*IN
      JSTRID = JMAX + 2*JN
      KSTRID = KMAX + 2*KN
      ILE    = ISTRID*JSTRID
      IF(IFACE == 1)      THEN
         IC1 = 1+IN + (IXLO-1+JN)*ISTRID + (IYLO-1+KN)*ILE
         IC2 = 1+IN + (IXLO-1+JN)*ISTRID + (IYUP-1+KN)*ILE
         IC3 = 1+IN + (IXUP-1+JN)*ISTRID + (IYUP-1+KN)*ILE
         IC4 = 1+IN + (IXUP-1+JN)*ISTRID + (IYLO-1+KN)*ILE
      ELSE IF(IFACE == 2) THEN
         IC1 = IXLO+IN + (1-1+JN)*ISTRID +(IYLO-1+KN)*ILE
         IC2 = IXUP+IN + (1-1+JN)*ISTRID +(IYLO-1+KN)*ILE
         IC3 = IXUP+IN + (1-1+JN)*ISTRID +(IYUP-1+KN)*ILE
         IC4 = IXLO+IN + (1-1+JN)*ISTRID +(IYUP-1+KN)*ILE
      ELSE IF(IFACE == 3) THEN
         IC1 = IXLO+IN + (IYLO-1+JN)*ISTRID + (1-1+KN)*ILE
         IC2 = IXLO+IN + (IYUP-1+JN)*ISTRID + (1-1+KN)*ILE
         IC3 = IXUP+IN + (IYUP-1+JN)*ISTRID + (1-1+KN)*ILE
         IC4 = IXUP+IN + (IYLO-1+JN)*ISTRID + (1-1+KN)*ILE
      ELSE IF(IFACE == 4) THEN
         IC1 = IMAX+1+IN + (IXLO-1+JN)*ISTRID + (IYLO-1+KN)*ILE
         IC2 = IMAX+1+IN + (IXUP-1+JN)*ISTRID + (IYLO-1+KN)*ILE
         IC3 = IMAX+1+IN + (IXUP-1+JN)*ISTRID + (IYUP-1+KN)*ILE
         IC4 = IMAX+1+IN + (IXLO-1+JN)*ISTRID + (IYUP-1+KN)*ILE
      ELSE IF(IFACE == 5) THEN
         IC1 = IXLO+IN + (JMAX+1-1+JN)*ISTRID + (IYLO-1+KN)*ILE
         IC2 = IXLO+IN + (JMAX+1-1+JN)*ISTRID + (IYUP-1+KN)*ILE
         IC3 = IXUP+IN + (JMAX+1-1+JN)*ISTRID + (IYUP-1+KN)*ILE
         IC4 = IXUP+IN + (JMAX+1-1+JN)*ISTRID + (IYLO-1+KN)*ILE
      ELSE IF(IFACE == 6) THEN
         IC1 = IXLO+IN + (IYLO-1+JN)*ISTRID + (KMAX+1-1+KN)*ILE
         IC2 = IXUP+IN + (IYLO-1+JN)*ISTRID + (KMAX+1-1+KN)*ILE
         IC3 = IXUP+IN + (IYUP-1+JN)*ISTRID + (KMAX+1-1+KN)*ILE
         IC4 = IXLO+IN + (IYUP-1+JN)*ISTRID + (KMAX+1-1+KN)*ILE
      ELSE
         WRITE(*,1010) IFACE,L,IMAX,JMAX,KMAX
 1010    FORMAT('CUMBERSOME CORNER DATA. IFACE =',I3,'PATCH=',I3,
     +        ' (IMAX,JMAX,KMAX) = (',3I3,')')
         CALL EXIT
      ENDIF

      XCOL(1,1) = XCO(IC1)
      XCOL(2,1) = YCO(IC1)
      XCOL(3,1) = ZCO(IC1)
      XCOL(1,2) = XCO(IC2)
      XCOL(2,2) = YCO(IC2)
      XCOL(3,2) = ZCO(IC2)
      XCOL(1,3) = XCO(IC3)
      XCOL(2,3) = YCO(IC3)
      XCOL(3,3) = ZCO(IC3)
      XCOL(1,4) = XCO(IC4)
      XCOL(2,4) = YCO(IC4)
      XCOL(3,4) = ZCO(IC4)

      RETURN
      END SUBROUTINE CORNER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE READBC(IUNIT,NBLOCK,IX,IY,IZ,RCON,NBCS,BCFILE,
     &     TBC,IBC,IPRO,MGRID,INPAR,GROUP,TIMEL,ITIMES,LEVEL,
     &     FRESUL,COORL,IREPEA,NREPEA,BNFILE)
C ...
C     This subroutine reads the boundary condition file.
C
C     Valid boundary condition types are:
C
C     DUM : Dummy           (=  0)
C     CON : Connectivity    (=  1)
C     PER : Periodic        (=  1) +icon(20) = 1
C     EXT : External        (=  2)
C     INL : Inlet           (=  3)
C     MIR : Mirror          (=  4)
C     OUT : Outlet          (=  5)
C     CYC : Cyclic          (=  6)
C     SNG : Singularity     (=  7)
C     SOL : Solid           (=  8)
C     ROT : Rotating solid  (=  9)
C     MOV : Moving solid    (= 10)
C     SLD : Sliding patch   (= 11)
C     CHI : Chimera         (= 12)
C     FRE : Free surface    (= 13)
C     MIX : Mixing plane    (= 14)
C     HTS : Coupling with a solid block  (= 15)
C     SLI : Slip boundary   (=16)
C     INP : Cycling inlet   (=17)
C     CML : Connection with a pressure loss (= 18)
C     ACT : Actuator disc (= 19)
C
C
C     Input parameters:
C
C     IUNIT    : I/O unit number of the boundary condition file
C     NBLOCK   : Total number of grid blocks
C     IX,IY,IZ : Block dimensions
C
C
C     Output parameters:
C
C     IBC      : Boundary condition integer vector (size = NBCS*IC9)
C     RCON     : Boundary condition real vector    (size = NBCS*IC9)
C     NBCS     : Total number of BC pathes
C     BCFILE   : Force group files (voi mersu)
C     BNFILE   : Boundary condition files (in case of INL or OUT)
C ...

      USE NS3CO, ONLY : IC9
      USE MAIN_ARRAYS, ONLY : IGRID

      IMPLICIT NONE

      CHARACTER(LEN=80) :: LINE1, LINE2, LINE3
      CHARACTER(LEN=80) :: BCFILE(*), BNFILE(*)
      CHARACTER(LEN=80) :: ENSBND_FILE, FVBND_FILE
      CHARACTER(LEN=3)  :: CBC, TBC(*)

      INTEGER :: IFOUND,ITYPE,IUNIT,NBLOCK,IPRO,INPAR,ITIMES,
     &           LEVEL,I,J,K,IT,NBLKS,IB,IDIM,JDIM,KDIM,IFACE,NOFBCS,
     &           IFSIZE,ISTART,JSTART,M,LFSIZE,ILEN,IR,IPAREA,JPAREA,
     &           ITY,NGPN,NP,N,L,IP,IBEDEL,NPATCH,II,JJ,III

      INTEGER :: IPTRFG(52)                    ! Force group pointer
      INTEGER :: IBC(IC9,*), IHPEU(80), NBCS   ! BC Patches temporarily
      INTEGER :: IX(NBLOCK),IY(NBLOCK),IZ(NBLOCK),MGRID(INPAR,*),
     &           IREPEA(*),NREPEA(*)

      LOGICAL :: SPACES,GROUP,TIMEL,FRESUL,COORL,INLOUT,SOLIDS

      REAL    :: RCON(IC9,*),WTEMP,WFLUX,WINJM,WPORO,WINJT,WINJRX,
     &           WINJRY,WINJRZ,WRBK,WTRAN


C ... Initialize arrays

      DO I=1,NBCS

         BCFILE(I)(1:80) = ' '
         BNFILE(I)(1:80) = ' '

         DO J=1,IC9
            IBC(J,I) = 0
            RCON(J,I)= 0.
         ENDDO

      ENDDO


      NBCS = 0

                                      
      GROUP = .FALSE.
      IT = 45                              ! Output unit

      REWIND IUNIT                         ! Input unit
       
      READ(IUNIT,'(1A80)') LINE1           ! Data description
      READ(IUNIT,*) NBLKS                  ! Total number of blocks

      IF(NBLKS /= NBLOCK) THEN
         WRITE(*,313)  NBLOCK,NBLKS
         WRITE(IT,313) NBLOCK,NBLKS
         STOP
      ENDIF

 313  FORMAT(/'READBC: Number of grid blocks in the grid file (',
     +     I2,') and '/
     +        '        in the BC file (',I2,
     +     ') do not match ! Aborting ...'//)
      
C ... Start the block loop

      DO 13 IB=1,NBLKS                     

        READ(IUNIT,'(1A80)') LINE1         ! Comment line only
        READ(IUNIT,*) IDIM, JDIM, KDIM     ! Block size (in cells)

        IF(IDIM+1 /= IX(IB) .OR.
     &     JDIM+1 /= IY(IB) .OR.
     &     KDIM+1 /= IZ(IB))
     &  THEN
          WRITE(*,*)
          WRITE(*,1101) IB
          WRITE(*,*)
          WRITE(IT,*)
          WRITE(IT,1101) IB
          WRITE(IT,*)
          IDIM = IX(IB)-1
          JDIM = IY(IB)-1
          KDIM = IZ(IB)-1
        ENDIF

 1101 FORMAT(' READBC: ERROR in the block dimensions in the BC file: ',
     &      'block ',I3/'         Using dimensions from the grid file.')

C ... Start the face loop inside the block loop

      DO 12 I=1,6                 
      READ(IUNIT,*) IFACE,NOFBCS       ! Face number, NOF BC patches
      IFSIZE     = 0
      DO 11 J=1,NOFBCS                 ! Start BC patch loop

            READ(IUNIT,'(1A80)') LINE1

C ... Skip extra spaces in the begining of an input line

            ISTART = 0
            DO 10 M=1,80
              ISTART = ISTART + 1
              IF(LINE1(M:M) /= ' ') GOTO 19
 10         CONTINUE
 19         CONTINUE

            LINE2(1:80) = ' '
            LINE2(1:78-ISTART) = LINE1(ISTART+3:80)
            CBC                = LINE1(ISTART:ISTART+3)
            IF(IFACE == 1 .OR. IFACE == 4) THEN
              LFSIZE = JDIM * KDIM
            ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
              LFSIZE = IDIM * KDIM
            ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
              LFSIZE = IDIM * JDIM
            ELSE
              WRITE(*,*)
              WRITE(*,1001) IFACE,IB
              WRITE(*,*)
              WRITE(IT,*)
              WRITE(IT,1001) IFACE,IB
              WRITE(IT,*)
              STOP
 1001         FORMAT(' READBC: Illegal face number ',I3,' found in ',
     &               'block ',I3,'.')
            ENDIF

            NBCS = NBCS + 1
            TBC(NBCS) = CBC

C ... Pick up the force group identifiers from the input line

            IF(  CBC == 'SOL' .OR. CBC == 'ROT' .OR. CBC == 'MOV' 
     &      .OR. CBC == 'HTS' .OR. CBC == 'SLI' .OR. CBC == 'CML'
     &      .OR. CBC == 'INL' .OR. CBC == 'OUT' .OR. CBC == 'ACT')
     &      CALL FGIDS(LINE2,BCFILE,NBCS,GROUP,TIMEL,ITIMES)

            IF    (CBC == 'DUM') THEN
              IBC(1,NBCS) =  0
            ELSEIF(CBC == 'CON') THEN
              IBC(1,NBCS) =  1
            ELSEIF(CBC == 'PER') THEN
              IBC(1,NBCS) =  1
              IBC(20,NBCS)=  1
              CBC         = 'CON'
            ELSEIF(CBC == 'EXT') THEN
              IBC(1,NBCS) =  2
            ELSEIF(CBC == 'INL') THEN
              IBC(1,NBCS) =  3
            ELSEIF(CBC == 'MIR') THEN
              IBC(1,NBCS) =  4
            ELSEIF(CBC == 'OUT') THEN
              IBC(1,NBCS) =  5
            ELSEIF(CBC == 'CYC') THEN
              IBC(1,NBCS) =  6
            ELSEIF(CBC == 'SNG') THEN
              IBC(1,NBCS) =  7
            ELSEIF(CBC == 'SOL') THEN
              IBC(1,NBCS) =  8
            ELSEIF(CBC == 'ROT') THEN
              IBC(1,NBCS) =  9
            ELSEIF(CBC == 'MOV') THEN
              IBC(1,NBCS) = 10
            ELSEIF(CBC == 'SLD') THEN
              IBC(1,NBCS) = 11
            ELSEIF(CBC == 'CHI') THEN
              IBC(1,NBCS) = 12
            ELSEIF(CBC == 'FRE') THEN
              IBC(1,NBCS) = 13
              IF(.NOT.FRESUL) THEN
                 FRESUL = .TRUE.
                 WRITE(45,*) '  A free-surface found in',
     &           ' READBC, FRESUL and COORL were changed to TRUE'
                 WRITE(*,*)  '  A free-surface found in',
     &           ' READBC, FRESUL and COORL were changed to TRUE'
              ENDIF
            ELSEIF(CBC == 'MIX') THEN
              IBC(1,NBCS) = 14
            ELSEIF(CBC == 'HTS') THEN
              IBC(1,NBCS) = 15
              IBC(20,NBCS)=  4  ! Coupling with a heat transfer surface
            ELSEIF(CBC == 'SLI') THEN
              IBC(1,NBCS) = 16
            ELSEIF(CBC == 'INP') THEN
              IBC(1,NBCS) = 17
            ELSEIF(CBC == 'CML') THEN
              IBC(1,NBCS) =  1
              IBC(20,NBCS)=  7
              CBC         = 'CON'
            ELSEIF(CBC == 'ACT') THEN
              IBC(1,NBCS) =  1
              IBC(20,NBCS)=  8
              CBC         = 'CON'
            ELSE
              WRITE(*,1111)  CBC,NBCS
              WRITE(IT,1111) CBC,NBCS
 1111         FORMAT(/'READBC: Unknown boundary condition ',A3,
     +             'in patch',I5)
              STOP
            ENDIF

            IBC(2,NBCS) = IB
            IBC(3,NBCS) = IFACE
            IBC(26,NBCS) = IGRID(IB,3) ! Number of object is stored here

            IF(IFACE == 4) THEN
               IBC(21,NBCS) = IDIM
            ELSEIF(IFACE == 5) THEN
               IBC(21,NBCS) = JDIM
            ELSEIF(IFACE == 6) THEN
               IBC(21,NBCS) = KDIM
            ELSE
               IBC(21,NBCS) = 1
            ENDIF
               
            IF(CBC == 'CON' .OR. CBC == 'CYC' .OR. CBC == 'HTS' .OR.
     +         CBC == 'SLD' .OR. CBC == 'MIX') THEN ! A connective patch

               CALL INLINE(LINE2,8,IBC(4,NBCS),IFOUND)

*               IF (IFOUND /= 7 .AND. IFOUND /= 8) THEN
*                 WRITE(*,*) 'READBC: CALL INLINE(7)'
*                 WRITE(*,*) 'IN PATCH',NBCS,'IFOUND=',IFOUND
*                 WRITE(*,*) 'The following card is broken or erroneous:'
*                 WRITE(*,*) LINE1
*                 STOP
*               ENDIF
               IF(IFOUND < 7) THEN  ! Force an automatic connection
                  IBC(8,NBCS) = 0
                  IFOUND = 7
               ENDIF
               IBC(28,NBCS) = IBC(8,NBCS)  ! Global connective block #

               IF(CBC == 'CON' .AND. IBC(20,NBCS) == 7) THEN 
                  READ(IUNIT,'(1A80)') LINE1 ! Read a pressure loss

               READ(LINE1,*,END=2021,ERR=2023) 
     &         (RCON(IR,NBCS),IR=1,2)  ! K-value and porosity
               GO TO 2025
2021           READ(LINE1,*,END=2022,ERR=2023)
     & 	       (RCON(IR,NBCS),IR=1,1)  ! Pressure loss 
                RCON(1,2) = 1.   ! Porosity = 1.
                NREPEA(7) = 5
                IREPEA(7) = IREPEA(7) + 1

               IF(IREPEA(7) <= NREPEA(7)) THEN
               WRITE(13,'(2A,I4)')' Porosity is missing from BC-file.',
     &                ' POROS = 1. is assumed for patch',NBCS

               ELSEIF(IREPEA(7) == NREPEA(7)+1) THEN 
                  WRITE(13,'(A,I2,A)')'  THIS HAS BEEN REPEATED',
     &            IREPEA(7)-1,' TIMES. SILENCE'
               ENDIF

               GO TO 2025	 

 2022          WRITE(*,*) 'Wrong number of parameter in patch',NBCS
               STOP "READBC: After 'CML' not enough data"
2023           WRITE(*,*) 'The following card is broken or erroneous:'
               WRITE(*,*) LINE1
               STOP 'READBC'
2025           CONTINUE               ! This line succesfully read
               ENDIF ! CBC == CON .AND. IBC(20,NBCS) == 7

               IF(CBC == 'CON' .AND. IBC(20,NBCS) == 8) THEN
               READ(IUNIT,'(1A80)') LINE1  ! Read actuator disc data
c               READ(LINE1,*,END=2031,ERR=2033) 
c     & 	       (RCON(IR,NBCS),IR=1,2)      ! K-value and porosity
c               GO TO 2035
c2031           READ(LINE1,*,END=2032,ERR=2033)
c     & 	       (RCON(IR,NBCS),IR=1,1)      ! Pressure loss 
c                RCON(1,2) = 1.             ! Porosity = 1.
c                NREPEA(9) = 5
c                IREPEA(9) = IREPEA(9) + 1
c               IF(IREPEA(9) <= NREPEA(9)) THEN
c               WRITE(13,'(2A,I4)')' Porosity is missing from BC-file.',
c     &                ' POROS = 1. is assumed for patch',NBCS
c               ELSEIF(IREPEA(9) == NREPEA(9)+1) THEN 
c               WRITE(13,'(A,I2,A)')'  THIS HAS BEEN REPEATED',
c     &                IREPEA(9)-1,' TIMES. SILENCE'
c               ENDIF
c               GO TO 2035	 
c2032           WRITE(*,*) 'Wrong number of parameter in patch',NBCS
c               STOP "READBC: After 'ACT' not enough data"
c2033           WRITE(*,*) 'The following card is broken or erroneous:'c
c	       WRITE(*,*) LINE1
c               STOP 'READBC'
c2035           CONTINUE               ! This line succesfully read


C ... It is assumed that there is only the file name. No diagnostics,
C ... if the card is missing for the actuator disc!! 

                RCON(1,NBCS) = 1 ! K-value
                RCON(2,NBCS) = 1.! Porosity, activate, if needed

                JSTART=0
                ILEN  = 78
                DO 2030 M=1,ILEN
                  JSTART = JSTART + 1
                  IF(LINE1(M:M) /= ' ') GOTO 2031
2030             CONTINUE
2031             CONTINUE
            
                ILEN=81
                DO 2032 M=80,1,-1
                  ILEN = ILEN - 1
c                  WRITE(678,*) M,LINE1(M:M)
                  IF(LINE1(M:M) /= ' ') GOTO 2033
2032             CONTINUE
2033             CONTINUE

         BNFILE(NBCS)(1:ILEN) = LINE1(JSTART:JSTART+ILEN-1)

               ENDIF ! CBC == CON .AND. IBC(20,NBCS) == 8


C ... DIRECTIVES FOR MULTIPROCESSING AND ACTUATOR DISC FILE

              IBC(17,NBCS) = IPRO
              IBC(18,NBCS) = 1
              IBC(14,NBCS) = J

              IF (IFOUND == 8 .AND.TBC(NBCS) /= 'ACT') THEN ! What is this?
                 IBC(18,NBCS) = IBC(11,NBCS) 
              ENDIF
            ENDIF ! A connective patch

C ... Read Inlet or Outlet file name

            IF(CBC == 'INL' .OR. CBC == 'OUT') THEN
              IF(NOFBCS > 1) THEN
                WRITE(*,*)
                WRITE(*,*)  'READBC: INLet or OUTlet does not cover ',
     &                      'whole block face! Aborting ...'
                WRITE(*,*)
                WRITE(IT,*)
                WRITE(IT,*) 'READBC: INLet or OUTlet does not cover ',
     &                      'whole block face! Aborting ...'
                WRITE(IT,*)
                STOP
              ELSE
                JSTART=0
                ILEN  = 78-ISTART
                DO 70 M=1,ILEN
                  JSTART = JSTART + 1
                  IF(LINE2(M:M) /= ' ') GOTO 71
 70             CONTINUE
 71             CONTINUE

               DO II = JSTART,ILEN
               JJ    = II + 1 - JSTART
               LINE3(JJ:JJ) = LINE2(II:II)
               ENDDO

               CALL INLINE(LINE3,1,IHPEU,IFOUND)
c               CALL INLINE(LINE2(JSTART:ILEN),1,IHPEU,IFOUND)

                IF (IFOUND /= 1) THEN
                   IBC(8,NBCS)=1
                ELSE
c ... varo tätä
                   IBC(8,NBCS)=IHPEU(1)
                   SPACES = .FALSE.

C ... Store the inlet or outlet type

                   IF(CBC == 'INL' .OR. CBC == 'OUT') THEN
                      IBC(27,NBCS) = IHPEU(1)
                   ENDIF

                   DO 76 M=JSTART+1,ILEN
                      IF (LINE2(M:M) == ' ') SPACES = .TRUE.
                      IF (LINE2(M:M) /= ' ' .AND. SPACES) GOTO 77
 76                CONTINUE
                   M = M-1
 77                CONTINUE
                   JSTART = M
                ENDIF

                BCFILE(NBCS)(1:ILEN-JSTART-1) = LINE2(JSTART:ILEN)

             ENDIF
            ENDIF

C ... Move IMOV to second line soon.

            IF(CBC == 'MOV') THEN
c              CALL RNLINE(LINE2,9,RCON(1,NBCS),IFOUND)
               CALL INLINE(LINE2,5,IBC(4,NBCS),IFOUND)
               IF (IFOUND /= 5) THEN ! Drop 5 to 4
                  WRITE(*,*) LINE2
                  WRITE(*,*) 'READBC: CALL INLINE(6)'
                  WRITE(*,*) 'IN PATCH',NBCS,'IFOUND=',IFOUND
                  STOP
               ENDIF

               IBC(9,NBCS)  = IBC(8,NBCS) ! IMOV stored to place 9
               IBC(8,NBCS)  = 0           ! Obsolate

c               IBC(4,NBCS) = NINT(RCON(1,NBCS))
c               IBC(5,NBCS) = NINT(RCON(2,NBCS))
c               IBC(6,NBCS) = NINT(RCON(3,NBCS))
c               IBC(7,NBCS) = NINT(RCON(4,NBCS))
c               IBC(9,NBCS) = NINT(RCON(5,NBCS))
c               RCON(1,NBCS)=      RCON(6,NBCS)
c               RCON(2,NBCS)=      RCON(7,NBCS)
c               RCON(3,NBCS)=      RCON(8,NBCS)
c               RCON(4,NBCS)=      RCON(9,NBCS)
c               RCON(5,NBCS)= 0.
c               RCON(6,NBCS)= 0.
c               RCON(7,NBCS)= 0.
c               RCON(8,NBCS)= 0.
c               RCON(9,NBCS)= 0.
            ENDIF

C ... If the whole face is covered by the same BC type, we can overwrite
C ... the data given by the user thus making sure that the indices are right.

            IF(NOFBCS == 1) THEN
              IF(IFACE == 1 .OR. IFACE == 4) THEN
                IBC(4,NBCS) = 1
                IBC(5,NBCS) = JDIM
                IBC(6,NBCS) = 1
                IBC(7,NBCS) = KDIM
              ENDIF
              IF(IFACE == 2 .OR. IFACE == 5) THEN
                IBC(4,NBCS) = 1
                IBC(5,NBCS) = IDIM
                IBC(6,NBCS) = 1
                IBC(7,NBCS) = KDIM
              ENDIF
              IF(IFACE == 3 .OR. IFACE == 6) THEN
                IBC(4,NBCS) = 1
                IBC(5,NBCS) = IDIM
                IBC(6,NBCS) = 1
                IBC(7,NBCS) = JDIM
              ENDIF
            ENDIF

      IF((CBC /= 'CON' .AND. CBC /= 'CYC' .AND. CBC /= 'SLD'
     & .AND. CBC /= 'MIX' .AND. CBC /= 'MOV' .AND. CBC /= 'HTS')
     & .AND. NOFBCS > 1) THEN
         CALL INLINE(LINE2,4,IBC(4,NBCS),IFOUND)
         IF (IFOUND /= 4) THEN
            WRITE(*,*) 'Wrong number of parameter in patch',NBCS
            WRITE(*,*) 'The following card is broken or erroneous:'
            WRITE(*,*) LINE1
            STOP 'READBC: CALL INLINE(4)'
         ENDIF ! IFOUND
      ENDIF

      IFSIZE = IFSIZE + (IBC(5,NBCS) - IBC(4,NBCS) + 1)
     &                * (IBC(7,NBCS) - IBC(6,NBCS) + 1)

C ************************************************************************
C ... Read in the extra line for solid boundaries (ancient THERMO.BC info)
C ************************************************************************

      IF(IBC(1,NBCS) /= 1 .AND. IBC(1,NBCS) /= 15) 
     &IBC(20,NBCS)  = 0                     ! Nämät olivat ennen RRCONissa?
      RCON( 5,NBCS) = 0.
      RCON( 6,NBCS) = 0.

C ... Extra line for solids

      IF(CBC == 'SOL' .OR. CBC == 'ROT' .OR. CBC == 'MOV') THEN

         READ(IUNIT,'(1A80)') LINE1 
         READ(LINE1,*,END=2011,ERR=2013) 
     &   IBC(20,NBCS),(RCON(IR,NBCS),IR=1,14) ! Add IMOV here on the 2nd line
         GO TO 2015

2011     READ(LINE1,*,END=2012,ERR=2013)
     &   IBC(20,NBCS),(RCON(IR,NBCS),IR=1,13) 
         RCON(14,NBCS)=1        ! WTRAN=1, fully turbulent flow   
         NREPEA(5) = 5
         IREPEA(5) = IREPEA(5) + 1
         IF(IREPEA(5) <= NREPEA(5)) THEN
            WRITE(13,'(2A,I4)')'  WTRAN is missing from the BC-file.',
     &          ' Turbulence production is not blocked for patch',NBCS
         ELSEIF(IREPEA(5) == NREPEA(5)+1) THEN 
            WRITE(13,'(A,I2,A)')'  THIS HAS BEEN REPEATED',IREPEA(5)-1,
     &                ' TIMES. SILENCE'
         ENDIF

         GO TO 2015	 

2012     WRITE(*,*) 'Wrong number of parameter in patch',NBCS
         STOP "READBC: After 'SOL' or 'ROT' or 'MOV' not enough data"

2013     WRITE(*,*) 'The following card is broken or erroneous:'
         WRITE(*,*) LINE1
         STOP

2015     CONTINUE               ! This line succesfully read without RNLINE

      ENDIF        ! Solid surface ('SOL', 'ROT', 'MOV')

C ... Extra line for INL/OUT if ITYPE <= 10 or 20 < ITYPE <= 30

      ITYPE = IBC(27,NBCS)

      IF((CBC == 'INL' .OR. CBC == 'OUT') .AND. 
     &  (ITYPE <= 10 .OR. ITYPE > 20 .AND. ITYPE <=30)) THEN

         READ(IUNIT,'(1A80)') LINE1 

C ... It is assumed that there is only the file name. No diagnostics,
C ... if the card is missing!!

                JSTART=0
                ILEN  = 78-ISTART
                DO 2070 M=1,ILEN
                  JSTART = JSTART + 1
                  IF(LINE1(M:M) /= ' ') GOTO 2071
2070             CONTINUE
2071             CONTINUE
            
                ILEN=81
                DO 2072 M=80,1,-1
                  ILEN = ILEN - 1
c                  WRITE(678,*) M,LINE1(M:M)
                  IF(LINE1(M:M) /= ' ') GOTO 2073
2072             CONTINUE
2073             CONTINUE

         BNFILE(NBCS)(1:ILEN) = LINE1(JSTART:JSTART+ILEN-1)

      ENDIF        ! INLET/OUTLET surface ('INL', 'OUT')

 11   CONTINUE     ! End BC Patch loop

**********************************************************************

          IF(LFSIZE /= IFSIZE) THEN
            WRITE(*,*)
            WRITE(*,1002) IB,IFACE
            WRITE(*,*)
            WRITE(IT,*)
            WRITE(IT,1002) IB,IFACE
            WRITE(IT,*)
            STOP
 1002       FORMAT(' READBC: Error in the BC patch dimensions in ',
     &               'block ',I3,', face ',I2,'.')
          ENDIF
         

 12     CONTINUE       ! End of face loop
 13   CONTINUE         ! End of block loop


      CALL AUTOCON(IBC,NBCS) 


C ... Complete the connectivity array, i.e., replace the local patch
C ... numbers with the global batch numbers and the patch dimensions

      DO 32 I = 1,NBCS
        DO 31 J = 1,NBCS
          IF((IBC( 1,J) ==  1 .OR. IBC( 1,J) == 6   .OR. 
     &        IBC( 1,J) == 11 .OR. IBC( 1,J)  == 14   .OR.
     &        IBC( 1,J) == 15)  .AND.
     &        IBC( 2,J)  == IBC( 8,I) .AND.
     &        IBC( 3,J)  == IBC( 9,I) .AND.
     &   IABS(IBC(14,J)) == IBC(10,I)) THEN
              IBC(10,I) =  IBC( 4,J)
              IBC(11,I) =  IBC( 5,J)
              IBC(12,I) =  IBC( 6,J)
              IBC(13,I) =  IBC( 7,J)
              IPAREA = (IBC(5,I)-IBC(4,I)-1)*(IBC(7,I)-IBC(6,I)-1)
              JPAREA = (IBC(5,J)-IBC(4,J)-1)*(IBC(7,J)-IBC(6,J)-1)
              IF(IPAREA /= JPAREA) THEN
                 IF(IBC(1,J) == 1) THEN
                    IBC(22,J) = 1 ! nonmatching connection
                 ELSE
                WRITE(*,*)
                IF(IBC(1,I) == 1) WRITE(*,1003) IBC(2,I), IBC(3,I)
                IF(IBC(1,I) == 6) WRITE(*,1033) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 11) WRITE(*,1333) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 14) WRITE(*,1334) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 15) WRITE(*,1335) IBC(2,I), IBC(3,I)
                WRITE(*,*) '        Patch areas do not match!'
                WRITE(*,*)
                WRITE(IT,*)
                IF(IBC(1,I) == 1) WRITE(IT,1003) IBC(2,I), IBC(3,I)
                IF(IBC(1,I) == 6) WRITE(IT,1033) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 11) WRITE(IT,1333) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 14) WRITE(IT,1334) IBC(2,I), IBC(3,I)
                IF(IBC(1,I)  == 15) WRITE(IT,1335) IBC(2,I), IBC(3,I)
                WRITE(IT,*) '        Patch areas do not match!'
                WRITE(IT,*)
                STOP
              ENDIF
               ENDIF
              IBC(14,I) = -IBC(14,I)
              IBC(15,I) =  J
              GOTO 32
          ENDIF
 31    CONTINUE
 32   CONTINUE

C ... Complete the mirror data

      DO 55 I = 1,NBCS
         ITY  = IBC(1,I)
         IF(ITY == 4 .OR. ITY == 12 .OR. ITY == 13) THEN
            IBC( 8,I) = IBC( 2,I)
            IBC( 9,I) = IBC( 3,I)
            IBC(10,I) = IBC( 4,I)
            IBC(11,I) = IBC( 5,I)
            IBC(12,I) = IBC( 6,I)
            IBC(13,I) = IBC( 7,I)
            IBC(14,I) = 0
            IBC(15,I) = I
            IBC(18,I) = 1
         ENDIF
C ... SAVE THE MGRID LEVELS OF THE CONNECTIVE BLOCK
         IF(ITY ==  1 .OR. ITY ==  4 .OR. ITY ==  6 .OR. 
     +      ITY == 11 .OR. ITY == 13 .OR. ITY == 14 .OR.
     +      ITY == 15) THEN
            IBC(19,I) = MGRID(1,IBC(8,I))
         ENDIF

 55   CONTINUE

C ... Check that all connectivities were found

      DO 14 I = 1,NBCS
        IF(IBC(14,I) > 0) THEN
          WRITE(*,*)
          IF(IBC(1,I) == 1) WRITE(*,1003) IBC(2,I), IBC(3,I)
          IF(IBC(1,I) == 6) WRITE(*,1033) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 11) WRITE(*,1333) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 14) WRITE(*,1334) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 15) WRITE(*,1335) IBC(2,I), IBC(3,I)
          WRITE(*,*)
          WRITE(IT,*)
          IF(IBC(1,I) == 1) WRITE(IT,1003) IBC(2,I), IBC(3,I)
          IF(IBC(1,I) == 6) WRITE(IT,1033) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 11) WRITE(IT,1333) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 14) WRITE(IT,1334) IBC(2,I), IBC(3,I)
          IF(IBC(1,I)  == 15) WRITE(IT,1335) IBC(2,I), IBC(3,I)
          WRITE(IT,*)
          STOP'vika 2'
 1003     FORMAT(' READBC: Error in the block connectivity data: ',
     &                 'block ',I3,', face ',I2,'.')
 1033     FORMAT(' READBC: Error in the cyclic connectivities: ',
     &                 'block ',I3,', face ',I2,'.')
 1333     FORMAT(' READBC: Error in the sliding connectivities: ',
     &                 'block ',I3,', face ',I2,'.')
 1334     FORMAT(' READBC: Error in the mixing connectivities: ',
     &                 'block ',I3,', face ',I2,'.')
 1335     FORMAT(' READBC: Error in the HTS connectivities: ',
     &                 'block ',I3,', face ',I2,'.')
        ENDIF
 14   CONTINUE

C ... CHECK THAT ALL CONNECTIONS MATCH EACH OTHER

      DO 140 I = 1,NBCS
        IF(IBC(1,I) == 1 .OR. IBC(1,I) == 6 .OR. IBC(1,I) == 11 .OR.
     &     IBC(1,I) == 14) THEN
           NGPN = IBC(15,I)
           IF(IBC(8,I) /= IBC(2,NGPN)) THEN
              WRITE(*,1203) I,IBC(2,I),NGPN,IBC(8,I),
     &             IBC(15,NGPN),IBC(8,NGPN)
              WRITE(IT,1203) I,IBC(2,I),NGPN,IBC(8,I),
     &             IBC(15,NGPN),IBC(8,NGPN)
              STOP
           ENDIF
           IF(IBC(9,I) /= IBC(3,NGPN)) THEN
              WRITE(*,1204) IBC(3,I),IBC(2,I),IBC(3,NGPN),IBC(8,I),
     &             IBC(9,NGPN),IBC(8,NGPN)
              WRITE(IT,1204) IBC(3,I),IBC(2,I),IBC(3,NGPN),IBC(8,I),
     &             IBC(9,NGPN),IBC(8,NGPN)
              STOP
           ENDIF
            IF(I /= IBC(15,NGPN)) THEN
              WRITE(*,1205) I,IBC(2,I),NGPN,IBC(8,I),
     &             IBC(15,NGPN),IBC(8,NGPN)
              WRITE(IT,1205) I,IBC(2,I),NGPN,IBC(8,I),
     &             IBC(15,NGPN),IBC(8,NGPN)
              STOP
           ENDIF
           IF(IBC(1,I) /= IBC(1,NGPN)) THEN
              WRITE(*,1206) I,IBC(2,I),IBC(1,I),NGPN,IBC(8,I),
     &             IBC(1,NGPN)
              WRITE(IT,1206) I,IBC(2,I),IBC(1,I),NGPN,IBC(8,I),
     &             IBC(1,NGPN)
              STOP
           ENDIF
             
 1203      FORMAT(' READBC: Error in the block connectivity data: '/
     &          'patch ',I3,', in block ',I2,' is connected in'/
     &          'patch ',I3,', in block ',I2,', which is connected to'/
     &          'patch ',I3,', in block ',I2,'.'/'Exiting ...')
 1204      FORMAT(' READBC: Error in the face connectivity data: '/
     &          'face  ',I3,', in block ',I2,' is connected in'/
     &          'face  ',I3,', in block ',I2,', which is connected to'/
     &          'face  ',I3,', in block ',I2,'.'/'Exiting ...')
 1205      FORMAT(' READBC: Error in the patch connectivity data: '/
     &          'patch ',I3,', in block ',I2,' is connected in'/
     &          'patch ',I3,', in block ',I2,', which is connected to'/
     &          'patch ',I3,', in block ',I2,'.'/'Exiting ...')
 1206      FORMAT(' READBC: Error in the patch connectivity data: '/
     &          'patch ',I3,', in block ',I2,' which is type ',I2,
     &          'is connected in'/
     &          'patch ',I3,', in block ',I2,' which is type ',I2,
     &          '.'/'Exiting ...')
 1233      FORMAT(' READBC: Error in the cyclic connectivities: ',
     &          'block ',I3,', face ',I2,'.')
        ENDIF
 140  CONTINUE


C ... For clarity, clean the orientation patch numbers and save
C ... global patch numbers

      DO 35 I=1,NBCS
        IBC(14,I) = 0
        IBC(16,I) = IBC(15,I)
        IBC(24,I) = IBC(2,I)   !Global block number
        IBC(25,I) = I          !Global patch number
 35   CONTINUE
C ... find the block local patch number
       NP         = 0
       DO 2100 N  = 1,NBLOCK
          L = 0
 1110     L = L + 1
          NP         = NP + 1
          IF(NP > NBCS) THEN
             WRITE(*,*)'Error in finding local pathc number'
             WRITE(*,*)'Exiting ...'
             STOP
          ENDIF
          IF(IBC(16,NP) > 0) THEN
             IP     = IBC(16,NP)
             IBC(15,IP) = L
          ENDIF
          IF(N ==  IBC(2,NP+1))GOTO 1110
 2100  CONTINUE

C ... Write out the boundary condition patches for checking purposes

      WRITE(IT,1102)
 1102 FORMAT(//4X,'BOUNDARY CONDITION PATCHES ON THE FINEST ',
     &            'GRID LEVEL:',
     &       //4X,'CON = Connectivity    EXT = External     ',
     &            'INL = Inlet',
     &        /4X,'MIR = Mirror          OUT = Outlet       ',
     &            'CYC = Cyclic',
     &        /4X,'SNG = Singularity     SOL = Solid        ',
     &            'ROT = Rotating solid',
     &        /4X,'MOV = Moving solid    SLD = Sliding      ',
     &            'CHI = Chimera',
     &        /4X,'FRE = Free surface    MIX = Mixing       ',
     &            'HTS = Coupling with solid',
     &        /4X,'SLI = Slip surface    INP = Circulating inlet',
     &        /4X,'CML = Minor loss      ACT = Actuator disc')
c 1103 FORMAT(//4X,'BLOCK NO. ',I4,'                     CONNECTIVITY'//,
c     &         4X,'BC   BLK FACE  XLO  XUP  YLO  YUP BLK  FACE ',
c     &'XLO  XUP  YLO  YUP EMP GPN GPN PRO  TY MGS PER'/)
 1103 FORMAT(//4X,'BLOCK NO.',I6,'                    CONNECTIVITY'//,
     &         4X,'BC   BLK FACE XLO  XUP  YLO  YUP    BLK  FACE ',
     &'XLO  XUP  YLO  YUP EMP GPN   GPN PRO  TY MGS PER'/)

      IBEDEL = 0
      DO 17 I=1,NBCS
        IF(IBC(2,I) > IBEDEL) THEN
          WRITE(IT,1103) IBC(2,I)
          IBEDEL = IBC(2,I)
        ENDIF
        WRITE(IT,111) TBC(I),(IBC(J,I),J=2,20)
c 111    FORMAT(4X,1A3,I4,1X,I4,4I5,2X,I3,2X,I3,4I5,I4,I4,I4,4I4)
 111    FORMAT(3X,1A3,I6,I4,4I5,1X,I6,I5,4I5,I4,I4,I6,4I4)
 17   CONTINUE
      WRITE(IT,*)
      WRITE(IT,*)

C ... Do intructions for a super connection. This means that
C ... every connection has a pair that will be threated same order
C ... in different processors. Here must put boundary conditions that
C ... needed to be connected. Now CON, CYC, SLD and MIX are connections like
C ... that.


C ... INSIDE OF RCON VECTOR
C ... 1,2,3,4 = SPEFICATION OF THE MOVING SOLID
C ... 5 = WALL TEMPERATURE
C ... 6 = WALL HEAT FLUX
C ... 7 = WALL INJECTION OR SUCTION (NEGATIVE)
C ... 8 = WALL POROSITY
C ... 9 = THE STAGNATION ENTHALPY OF THE INJECTED MASS FLOW (WITH INJECTION)
C ... 10...12 = THE R_suc vector (WITH INJECTION)
C ... 13 = SURFACE ROUGHNESS
C ... 14 = TYPE OF THE BOUNDARY LAYER 0= FULLY LAMINAR 1=FULLY TURBULENT 
ccc      CALL RRCON(IBC,RCON,NBCS,BCFILE)

      IBEDEL = 0
      WRITE(IT,1142)
 1142 FORMAT(//4X,'SOLID-WALL RCON ARRAY FOR THE FINEST GRID LEVEL:')
 1143 FORMAT(/4X,'BLOCK NO. ',I4//,
     &         4X,'PT ITYP    R1        R2        R3        VEL    ',
     &'  WALLTE WALL HEAT      SUC/INJ    POROSITY INJ.TEMP.  RSUCX',
     &'    RSUCY      RSUCZ   RBK       WTRAN')
 1144 FORMAT(4X,130('='))

      DO I=1,NBCS
        IF(IBC(2,I) > IBEDEL .AND.
     &	IBC(1,I) >= 8 .AND. IBC(1,I) <= 10) THEN
          WRITE(IT,*)
          WRITE(IT,1143) IBC(2,I)
          WRITE(IT,1144)
c	         WRITE(IT,*)
          IBEDEL = IBC(2,I)
        ENDIF
        IF(IBC(1,I) >= 8 .AND. IBC(1,I) <= 10) THEN
           WRITE(IT,141) I,IBC(9,I),(RCON(J,I),J=1,14)
        ENDIF
c 141    FORMAT(2X,I4,1X,I2,5(1X,F8.2),1X,E10.5,6(1X,F8.2),(1X,F8.6),
c     &        (2X,F4.2))
 141    FORMAT(I6,1X,I2,5(1X,F9.2),1X,E11.4,6(1X,F9.2),(1X,F11.4),
     &        (2X,F4.2))
      ENDDO
      WRITE(IT,*)

C ... Minor losses

      IBEDEL = 0
      WRITE(IT,1242)
 1242 FORMAT(//4X,'MINOR-LOSS RCON ARRAY FOR THE FINEST GRID LEVEL:')
 1243 FORMAT(/4X,'BLOCK NO. ',I4//,
     &         4X,'PT ITYP    RK     POROSITY')
 1244 FORMAT(4X,26('='))

      DO I=1,NBCS
        IF(IBC(2,I) > IBEDEL .AND.
     &	IBC(1,I) == 1 .AND. IBC(20,I) == 7) THEN
          WRITE(IT,*)
          WRITE(IT,1243) IBC(2,I)
          WRITE(IT,1244)
c	         WRITE(IT,*)
          IBEDEL = IBC(2,I)
        ENDIF
        IF(IBC(1,I) == 1 .AND. IBC(20,I) == 7) THEN
           WRITE(IT,142) I,IBC(20,I),(RCON(J,I),J=1,2)
        ENDIF
 142    FORMAT(2X,I4,1X,I2,2(1X,F8.2))
      ENDDO
      WRITE(IT,*)

C ... Heat transfer

      DO I = 1,NBCS

       N     = IBC(2,I)
       IFACE = IBC(3,I)
       NPATCH= I
       ITYPE = IBC(20,I)
       WTEMP = RCON(5,I)
       WFLUX = RCON(6,I)
       WINJM = RCON(7,I)
       WPORO = RCON(8,I)
       WINJT = RCON(9,I)
       WINJRX= RCON(10,I)
       WINJRY= RCON(11,I)
       WINJRZ= RCON(12,I)
       WRBK  = RCON(13,I)
       WTRAN = RCON(14,I)

       IF(ITYPE /= 0) THEN
         IF(ITYPE == 1) THEN
           WRITE(IT,'(3I4,A48,F9.2,A1)') N,IFACE,NPATCH,
     +    ' Block, face, patch has temperature specified as ',WTEMP,'K'   
         ELSEIF(ITYPE == 2) THEN
           WRITE(IT,'(3I4,A47,F9.2,A5)') N,IFACE,NPATCH,
     +   ' Block, face, patch has heat flux specified as ',WFLUX,'W/m^2'
         ELSEIF(ITYPE == 3) THEN
           WRITE(IT,'(3I4,A46,F9.2,A5)') N,IFACE,NPATCH,
     +   ' Block, face, patch have free-stream stagnation temperature '
         ELSEIF(ITYPE == 4) THEN
           WRITE(IT,'(3I4,A46,F9.2,A5)') N,IFACE,NPATCH,
     +   ' Block, face, patch are coupled with solid block '
         ELSEIF(ITYPE == 7) THEN
           WRITE(IT,'(3I4,A52)')         N,IFACE,NPATCH,
     +   ' Block, face, a connecting patch has a pressure loss '
         ELSEIF(ITYPE == 8) THEN
           WRITE(IT,'(3I4,A56,I2)')      N,IFACE,NPATCH,
     +   ' Block, face, a connecting patch has an actuator disc no.',
     +     IBC(26,I)
         ELSEIF(ITYPE == 11) THEN
           WRITE(IT,'(3I4,A45)')         N,IFACE,NPATCH,
     +   ' Block, face, patch has temperature specified in file '
         ELSEIF(ITYPE == 12) THEN
           WRITE(IT,'(3I4,A43)')         N,IFACE,NPATCH,
     +   ' Block, face, patch has heat flux specified in file '
         ELSEIF(ITYPE == 13 .OR. ITYPE == 14) THEN
           WRITE(IT,'(3I4,A42)')         N,IFACE,NPATCH,
     +   ' Block, face, patch is special. Look at SUB SURVEL '
         ELSE
           WRITE(*,*) 'No such heat wall type in Block, face, patch ',
     +     N,IFACE,NPATCH,ITYPE
           WRITE(*,*) 'Error in input file THERMO.BC'
           WRITE(*,*) 'Exiting ...'
           STOP
         ENDIF ! ITYPE == 1
       ENDIF   ! ITYPE == 0

       IF(WINJM /= 0.) THEN
            WRITE(IT,*) 
            WRITE(IT,'(3I4,A43,F9.2,A7)') N,IFACE,NPATCH,
     +     ' Block, face, patch has injection mass flux',WINJM,'kg/m^2s'   
            WRITE(IT,*) 
     +	   '         The corresponding porosity is',WPORO*100.,'%'   
            WRITE(IT,*)
     +     '         The stagnation entahlpy of the injection flow is ',
     +      '(instead of temperature)',WINJT,'K'  
c            WRITE(IT,*)
c     +     '         The temperature is not used if WINJM < 0 (suction)'   
            WRITE(IT,*) 
     +	   '         The R_suc vector is (',WINJRX,WINJRY,WINJRZ,')'   
       ENDIF      ! WINJM /= 0

      ENDDO ! NBCS into the MEMORY file
      WRITE(IT,*)      

C ... Write boundary information for Ensight and Fieldview

      ENSBND_FILE = 'XYZ.BIN.ensbnd'
      FVBND_FILE  = 'XYZ.BIN.fvbnd'

      CALL WRITE_BNDFILE(NBLOCK,IBC,NBCS,BCFILE,IPTRFG,LEVEL,
     &                   ENSBND_FILE,FVBND_FILE)

      DO I = 1,NBCS
         IBC(21,I) = 0
      ENDDO

      RETURN
      END SUBROUTINE READBC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RRCON(ICON,RCON,NBCS,BCFILE)

C ... THIS SUBROUTINE READ IN DATA FOR PATHC REAL ARRAY RCON

      USE NS3CO, ONLY : IC9
      
      DIMENSION :: ICON(IC9,*), RCON(IC9,*)
      CHARACTER(LEN=80) :: BCFILE(*)
      LOGICAL :: THERE

      DO 5 N =1,NBCS
         IF(ICON(1,N) /= 1) ICON(20,N) = 0
         RCON( 5,N) = 0.
         RCON( 6,N) = 0.
 5    CONTINUE

      INQUIRE(FILE='THERMO.BC',EXIST=THERE)
      WRITE(45,*)
      WRITE(45,*)

      IF(THERE) THEN
         WRITE(45,*) 'READING WALL THERMODYNAMICS FROM FILE THERMO.BC'
         OPEN(65,FILE='THERMO.BC',STATUS='OLD',FORM='FORMATTED')
 10      CONTINUE
         READ(65,*,END=100,ERR=12) NBL1,NPTC1,ITYPE,WTEMP,WFLUX
         NPTC2 = 0
         DO I = 1,NBCS
            NBL2 = ICON(2,I)
            IF(NBL2 == NBL1) THEN
               NPTC2 = NPTC2 + 1
               IF(NPTC2 == NPTC1) THEN
                  NPATCH = I
                  GOTO 11
               ENDIF
            ENDIF
         ENDDO
         WRITE(*,*) 'Does not find patch',NPTC1,' in block',NBL1
         WRITE(*,*) 'Exiting ...'
         STOP

 12   WRITE(*,*) 'Format of THERMO.BC is changed. Give file in format:'
         WRITE(*,*) 'Block_number block_local_patch_number ITYPE....'
         WRITE(*,*) 'Exiting ...'
         STOP

 11      CONTINUE
         IF(NPATCH > NBCS) THEN
            WRITE(*,*) 'There is not a pathc number of',NPATCH
            WRITE(*,*) 'Error in input file THERMO.BC'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
         IF(ICON(1,NPATCH) < 7 .OR. ICON(1,NPATCH) > 10) THEN
            WRITE(*,*) 
     +       'You give heat values in-correct wall at file THERMO.BC'
            WRITE(*,*) 'BC-TYPE',ICON(1,NPATCH),' cannot have heat flux'
            WRITE(*,*) 'Error in input file THERMO.BC'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
         ICON(20,NPATCH) = ITYPE
         RCON(5,NPATCH)  = WTEMP
         RCON(6,NPATCH)  = WFLUX
         N               = ICON(2,NPATCH)
         IFACE           = ICON(3,NPATCH)
         IF(ITYPE == 1) THEN
            WRITE(45,*) N,IFACE,NPATCH,
     +     ' Block, face, patch has temperature specified as',WTEMP,'C'   
         ELSEIF(ITYPE == 2) THEN
            WRITE(45,*) N,IFACE,NPATCH,
     + ' Block, face, patch has heat flux specified   as',WTEMP,'W/m^2'
         ELSEIF(ITYPE == 11) THEN
            WRITE(45,*) N,IFACE,NPATCH,
     +         ' Block, face, patch has temperature specified in file'
         ELSEIF(ITYPE == 12) THEN
            WRITE(45,*) N,IFACE,NPATCH,
     +      ' Block, face, patch has heat flux specified in file'
         ELSEIF(ITYPE == 13 .OR. ITYPE == 14) THEN
            WRITE(45,*) N,IFACE,NPATCH,
     +           ' Block, face, patch is special. look SUB SURVEL'
         ELSE
            WRITE(*,*) 'No such heat wall type in Block, face, pathc',
     +           N,IFACE,NPATCH,ITYPE
            WRITE(*,*) 'Error in input file THERMO.BC'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
         GOTO 10
 100     CLOSE(65)
      ENDIF  ! (THERE)
C ... tama on tehty ydinvoimalan polttoainesauvojen lammitykseen. PR 22.7.98
c      IAPU = 0
c      DO NP = 1,NBCS
c         ITYP = ICON(1,NP)
c         IF(ITYP >= 8 .AND. ITYP <= 10) THEN
c            DO L=1,20
c               IF(BCFILE(NP)(L:L) /= ' ') THEN
c                  ICON(20,NP) = 14
c                  RCON(6,NP)  = L + .5
c                  WRITE(45,*) NP,' is special. look SUB SURVEL'
c                  IAPU = IAPU + 1
c                  IF(L == 20) ICON(20,NP) = 15   ! T kirjain = reuna
c                  GOTO 110
c               ENDIF
c            ENDDO
c            BCFILE(NP)(25:25) = 'Y'
c 110        CONTINUE
c         ENDIF
c      ENDDO
c      IF(IAPU == 0 .AND. .NOT. THERE) THEN
      IF(THERE) THEN
         WRITE(45,*) 'No THERMO.BC-file or force group heats.'
         WRITE(45,*) 'All walls are adiabadic.'
c      ELSEIF(IAPU /= 0) THEN
c         WRITE(45,*) IAPU,' heated walls'
      ENDIF

      WRITE(45,*)
      WRITE(45,*)

      RETURN
      END
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      SUBROUTINE FGIDS(LINE2,BCFILE,NBCS,GROUP,TIMEL,ITIMES)

C ... Pick the force group identifiers from the input line and save
C ... them in the BCFILE. ASCII character set assumed.
C ... also open the channels for timeaccurate calculation

      CHARACTER(LEN=80) :: LINE2, BCFILE(*)
      LOGICAL :: GROUP, TIMEL, AVATTU, OPENED

      DO I=1,80

        IF(ICHAR(LINE2(I:I)) >= 65 .AND. ICHAR(LINE2(I:I)) <=  90) THEN
          L = ICHAR(LINE2(I:I))-64
          BCFILE(NBCS)(L:L) = LINE2(I:I)
          LINE2(I:I) = ' '
          GROUP = .TRUE.
          IF(TIMEL) THEN
             INQUIRE(UNIT=L+102,OPENED=AVATTU)
             IF(.NOT. AVATTU) THEN
                OPEN(L+102,FILE='FORCES.'//CHAR(L+64))
                IF(ITIMES == 0) THEN
                   WRITE(L+102,111)
                   WRITE(L+102,112)
                ELSE
                   READ(L+102,*)
                   READ(L+102,*)
		   ! CALL EOF(L+102,ITIMES+1)
                   CALL WTCFOR(L+102,ITIMES)
                ENDIF
             ENDIF
          ENDIF
        ENDIF

        IF(ICHAR(LINE2(I:I)) >= 97 .AND. ICHAR(LINE2(I:I)) <= 122) THEN
          L = ICHAR(LINE2(I:I))-70
          BCFILE(NBCS)(L:L) = LINE2(I:I)
          LINE2(I:I) = ' '
          GROUP = .TRUE.
          IF(TIMEL) THEN
             INQUIRE(UNIT=L+102,OPENED=AVATTU)
             IF(.NOT. AVATTU) THEN
                OPEN(L+102,FILE='FORCES.'//CHAR(L+70))
                IF(ITIMES == 0) THEN
                   WRITE(L+102,111)
                   WRITE(L+102,112)
                ELSE
                   READ(L+102,*)
                   READ(L+102,*)
                   CALL EOF(L+102,ITIMES+1)
                ENDIF
             ENDIF
          ENDIF
        ENDIF

      ENDDO

C ... Convergence history of force groups
      IF(GROUP) THEN  
         OPEN(155,FILE='IBDIAG_FGS',STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(156,FILE='TBDIAG_FGS',STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF

 111  FORMAT('#  Time (s)  ICYCLE  CL(N)       CD(N)       CS(N)',7X,
     &    'CX(N)       CY(N)       CZ(N)       DX(N)       DY(N)',7X,
     &    'DZ(N)       CMX(Nm)     CMY(Nm)     CMZ(Nm)',6X,
     &    'QT(W)       QW(W)       QH(W)     QT/A(W/M2)  QW/A(W/M2)',3X,
     &    'QH/A(W/m2)  TOM(W)      AREA(M2)')
 112  FORMAT('#  ',260('='))

      RETURN
      END  SUBROUTINE FGIDS           
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C

C     ------------------------------------------------------------------
C     The following three subroutines read INTEGER and REAL variables
C     from a CHARACTER*80 string to an INTEGER(*) or REAL(*) array.
C
C     INPUT:
C        LINE  (Character string to be parsed)
C        NV    (Integer, number of variables searched after)
C     OUTPUT
C        IVARS (Integer array, target array for the reading)
C        RVARS (Real    array, target array for the reading)
C        NVF   (Integer, number of variables found)
C
C     Petri Kaurinkoski, 2.2.1995
C     ------------------------------------------------------------------
C     ICON(IC,NP) - array deals with the data transport between patches
C     NP is a patch number and IC is various information:
C     IC    Type of information
C     1     Wall type (1 = Connectivy, 2 = External, 3 = Inlet, 
C                      4 = Mirror, 5 = Outlet, 6=Cyclic, 
C                      7=Singularity, 8 = Solid, 9 = Rotating solid,
C                      10 = Moving solid, 11 = Sliding patch,
C                      12 = Chimera, 13 = Free surface, 14 = Mixing surface,
C                      15 = Heat transfer surface)
C     2     Block number
C     3     Face number
C     4,5,6 & 7   xlow, xup, ylow & yup
C     8     Connective block number
C     9     Connective face number
C     10,11,12 & 13 Connective  xlow, xup, ylow & yup
C     14    Orientation
C     15    Connective patch number (global)
C     16    Connective prosess number
C     17    Method of tranport data between prosesses
C            (1=ZZZ-array, 2-testing)
C
C     PPR 11.10.1995
C     -----------------------------------------------------------------


      SUBROUTINE INLINE(LINE,NV,IVARS,NVF)

C ... READ IN INTEGER VARIABLES FROM THE LINE

      CHARACTER(LEN=80) :: LINE
      REAL :: REALIN
      INTEGER :: IVARS(*)
      LOGICAL :: ACCEPT

      LAST = 80
      IND1 = 1
      I    = 1
      NVF  = 0

 1    CONTINUE

C ... DEAL WITH THE REAL VALUES

      CALL UNWRAP(LINE,REALIN,IND1,LAST,LEN,ACCEPT)
      IF (ACCEPT) THEN
         NVF      = I
         IVARS(I) = NINT(REALIN)
         I        = I+1
         IND1     = IND1 + LEN
         IF (I <= NV) GOTO 1
      ENDIF

      RETURN
      END SUBROUTINE INLINE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RNLINE(LINE,NV,RVARS,NVF)

C ... READ IN REAL VARIABLES FROM THE LINE

      CHARACTER(LEN=80) :: LINE
      REAL :: RVARS(*), REALIN
      LOGICAL :: ACCEPT

      LAST = 80
      IND1 = 1
      I    = 1
      NVF  = 0

 1    CONTINUE

C ... DEAL WITH THE REAL VALUES

      CALL UNWRAP(LINE,REALIN,IND1,LAST,LEN,ACCEPT)
      IF (ACCEPT) THEN
         NVF      = I
         RVARS(I) = REALIN
         I        = I+1
         IND1     = IND1 + LEN
         IF (I <= NV) GOTO 1
      ENDIF

      RETURN
      END SUBROUTINE RNLINE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE UNWRAP(LINE,FUN,IFIRST,ILAST,LEN,ACCEPT)

C ... Find numbers in an input line (ASCII data set assumed)

      INTEGER :: IFUN, I, DES, DDES, NEG
      REAL :: FUN
      CHARACTER(LEN=80) :: LINE
      LOGICAL :: ACCEPT

      IFUN   = 0
      DES    = 0
      DDES   = 0
      NEG    = 1
      I      = IFIRST
      ACCEPT =.FALSE.
 1    CONTINUE
C      IF (LINE(I:I) == ' ') THEN
      IF ((ICHAR(LINE(I:I)) < 48 .OR. ICHAR(LINE(I:I)) > 57)
     & .AND. LINE(I:I) /= '.' .AND. LINE(I:I) /= '-') THEN
         I=I+1
         IF (I > ILAST) GOTO 10
         GOTO 1
      ENDIF

 2    CONTINUE
      IF (((ICHAR(LINE(I:I)) < 48 .OR. ICHAR(LINE(I:I)) > 57)
     & .AND. LINE(I:I) /= '.' .AND. LINE(I:I) /= '-').OR. I > ILAST)
     & THEN
         GOTO 10
       ELSE
         IF (LINE(I:I) == '.') THEN
            DDES=1
         ELSE
            IF (LINE(I:I) == '-') THEN
               NEG=-1
            ELSE
               IFUN   = IFUN*10 + (ICHAR(LINE(I:I)) - 48)
               DES    = DES+DDES
               ACCEPT =.TRUE.
            ENDIF
         ENDIF
         I=I+1
         GOTO 2
      ENDIF
 10   FUN=REAL(NEG*IFUN)/REAL(10**DES)

      LEN=I-IFIRST

      RETURN
      END SUBROUTINE UNWRAP
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRROR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,
     +     PRO,VAR,TRM,FI,KSCAL,ITURB,MULPHL,TRANSL,MAXSB,NX,NY,NZ,
     +     I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IN,JN,KN,IWALL,XVEL,ISCL,NSCAL,
     +     IWAVEB,FRESUL,IBTYPE,MULPHC,PRC,IPRESC)

C ... This subroutine rotates the vector components first into a
C ... coordinate system that has a coordinate direction normal to
C ... the boundary face, next the signs of the components perpendicular
C ... to the surface are changed (hence imposing no flow through the
C ... face) and last the components are rotated back to the original
C ... coordinate system. All this is accomplished by one matrix
C ... multiplication.

      USE TYPE_ARRAYS

      IMPLICIT NONE
 
      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF
      INTEGER :: IT,IT1,IT2,JL,IL,KSCAL,MAXSB,ISCL,NSCAL,KT1,KT2,
     +           I,J,IC1,IC2,IPHASE,NS,ITURB,IBTYPE,KB,KB1,KB2,IPRESC

      INTEGER :: IN, JN, KN, imax,jmax,kmax,ii,jj,kk

      INTEGER :: IWAVEB(*)

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),TTS(*),FI(MAXSB,*),NX(*), NY(*), NZ(*)
      REAL :: UMIR1,UMIR2,VMIR1,VMIR2,WMIR1,WMIR2,XVEL
      REAL :: A11,A12,A13,A21,A22,A23,A31,A32,A33

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(PRE_COR)          :: PRC(*)

      LOGICAL :: MULPHL, FRESUL, TRANSL

      CHARACTER(*) :: MULPHC

      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         STOP 'Illegal wall specification for subroutine mirror !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR  ! WALL VALUE FOR GEOMETRIC
      KB1 = (KW +   KDIR - 1 + KN)*KSTR  ! GHOST CELL 1            
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR  ! GHOST CELL 2            
      KT1 = (KW          - 1 + KN)*KSTR  ! COMPUTATIONAL DOMAIN 1  
      KT2 = (KW -   KDIR - 1 + KN)*KSTR  ! COMPUTATIONAL DOMAIN 2  

      IF(ISCL == 0) THEN
      DO 100 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB   ! WALL VALUE FOR GEOMETRIC
            IT1 = 1 + IL + KB1  ! GHOST CELL 1
            IT2 = 1 + IL + KB2  ! GHOST CELL 2
            IC1 = 1 + IL + KT1  ! COMPUTATIONAL DOMAIN 1
            IC2 = 1 + IL + KT2  ! COMPUTATIONAL DOMAIN 2
            A11 = XVEL*(-2.*NX(IT)*NX(IT) + 1.0)
            A12 = XVEL*(-2.*NX(IT)*NY(IT))
            A13 = XVEL*(-2.*NX(IT)*NZ(IT))
            A21 = XVEL*(-2.*NX(IT)*NY(IT))
            A22 = XVEL*(-2.*NY(IT)*NY(IT) + 1.0)
            A23 = XVEL*(-2.*NY(IT)*NZ(IT))
            A31 = XVEL*(-2.*NX(IT)*NZ(IT))
            A32 = XVEL*(-2.*NY(IT)*NZ(IT))
            A33 = XVEL*(-2.*NZ(IT)*NZ(IT) + 1.0)

            UMIR1 = U(IC1)
            VMIR1 = V(IC1)
            WMIR1 = W(IC1)

            UMIR2 = U(IC2)
            VMIR2 = V(IC2)
            WMIR2 = W(IC2)

            U(IT1) = A11*UMIR1 + A12*VMIR1 + A13*WMIR1
            V(IT1) = A21*UMIR1 + A22*VMIR1 + A23*WMIR1
            W(IT1) = A31*UMIR1 + A32*VMIR1 + A33*WMIR1

            U(IT2) = A11*UMIR2 + A12*VMIR2 + A13*WMIR2
            V(IT2) = A21*UMIR2 + A22*VMIR2 + A23*WMIR2
            W(IT2) = A31*UMIR2 + A32*VMIR2 + A33*WMIR2

            IF(IPRESC == 1) THEN ! Do always!

            UMIR1 = PRC(IC1)%DPDX
            VMIR1 = PRC(IC1)%DPDY
            WMIR1 = PRC(IC1)%DPDZ

            UMIR2 = PRC(IC2)%DPDX
            VMIR2 = PRC(IC2)%DPDY
            WMIR2 = PRC(IC2)%DPDZ
            PRC(IT1)%DPDX = A11*UMIR1 + A12*VMIR1 + A13*WMIR1
            PRC(IT1)%DPDY = A21*UMIR1 + A22*VMIR1 + A23*WMIR1
            PRC(IT1)%DPDZ = A31*UMIR1 + A32*VMIR1 + A33*WMIR1

            PRC(IT2)%DPDX = A11*UMIR2 + A12*VMIR2 + A13*WMIR2
            PRC(IT2)%DPDY = A21*UMIR2 + A22*VMIR2 + A23*WMIR2
            PRC(IT2)%DPDZ = A31*UMIR2 + A32*VMIR2 + A33*WMIR2
            ENDIF

            IF(MULPHC == 'MULTI' .AND. XVEL > 0.) THEN ! Do not always!

            DO IPHASE = 1,NPHASES

            UMIR1 = VAR(IC1)%U(IPHASE)
            VMIR1 = VAR(IC1)%V(IPHASE)
            WMIR1 = VAR(IC1)%W(IPHASE)

            UMIR2 = VAR(IC2)%U(IPHASE)
            VMIR2 = VAR(IC2)%V(IPHASE)
            WMIR2 = VAR(IC2)%W(IPHASE)
c      call ijkpai(it1,imax,jmax,kmax,ii,jj,kk)
            VAR(IT1)%U(IPHASE) = A11*UMIR1 + A12*VMIR1 + A13*WMIR1
            VAR(IT1)%V(IPHASE) = A21*UMIR1 + A22*VMIR1 + A23*WMIR1
            VAR(IT1)%W(IPHASE) = A31*UMIR1 + A32*VMIR1 + A33*WMIR1

            VAR(IT2)%U(IPHASE) = A11*UMIR2 + A12*VMIR2 + A13*WMIR2
            VAR(IT2)%V(IPHASE) = A21*UMIR2 + A22*VMIR2 + A23*WMIR2
            VAR(IT2)%W(IPHASE) = A31*UMIR2 + A32*VMIR2 + A33*WMIR2
        
            ENDDO
            ENDIF ! MULTI

            RO   (IT1) =    RO(IC1)
            E    (IT1) =     E(IC1)
            P    (IT1) =     P(IC1)
            PDIFF(IT1) = PDIFF(IC1)
            EPS2 (IT1) =  EPS2(IC1)
            VIST (IT1) =  VIST(IC1)

            RO   (IT2) =    RO(IC2)
            E    (IT2) =     E(IC2)
            P    (IT2) =     P(IC2)
            PDIFF(IT2) = PDIFF(IC2)
            EPS2 (IT2) =  EPS2(IC2)
            VIST (IT2) =  VIST(IC2)

            IF(IBTYPE == 16) THEN ! Slip boundary

            A11 = XVEL*(-2.*NX(IT)*NX(IT) + 1.0)
            A12 = XVEL*(-2.*NX(IT)*NY(IT))
            A13 = XVEL*(-2.*NX(IT)*NZ(IT))
            A21 = XVEL*(-2.*NX(IT)*NY(IT))
            A22 = XVEL*(-2.*NY(IT)*NY(IT) + 1.0)
            A23 = XVEL*(-2.*NY(IT)*NZ(IT))
            A31 = XVEL*(-2.*NX(IT)*NZ(IT))
            A32 = XVEL*(-2.*NY(IT)*NZ(IT))
            A33 = XVEL*(-2.*NZ(IT)*NZ(IT) + 1.0)

            UMIR1 = U(IC1) - U(IC2)
            VMIR1 = V(IC1)
            WMIR1 = W(IC1)

            U(IT1) = U(IT1) + A11*UMIR1 + A12*VMIR1 + A13*WMIR1
            V(IT1) = A21*UMIR1 + A22*VMIR1 + A23*WMIR1
            W(IT1) = A31*UMIR1 + A32*VMIR1 + A33*WMIR1

            UMIR2 = U(IC2)
            VMIR2 = V(IC2)
            WMIR2 = W(IC2)

            U(IT2) = A11*UMIR2 + A12*VMIR2 + A13*WMIR2
            V(IT2) = A21*UMIR2 + A22*VMIR2 + A23*WMIR2
            W(IT2) = A31*UMIR2 + A32*VMIR2 + A33*WMIR2

            RO   (IT1) =    2.*RO(IC1) - RO(IC2)
            E    (IT1) =     E(IC1)
            P    (IT1) =     P(IC1)
            PDIFF(IT1) = PDIFF(IC1)
            EPS2 (IT1) =  EPS2(IC1)
            VIST (IT1) =  VIST(IC1)

            RO   (IT2) =    RO(IC2)
            E    (IT2) =     E(IC2)
            P    (IT2) =     P(IC2)
            PDIFF(IT2) = PDIFF(IC2)
            EPS2 (IT2) =  EPS2(IC2)
            VIST (IT2) =  VIST(IC2)

            ENDIF
            IF(FRESUL) IWAVEB(IT1) = 4

 200     CONTINUE
 100  CONTINUE

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 400 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 300 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            RK  (IT1) =    RK(IC1)
            REPS(IT1) =  REPS(IC1)
            RK  (IT2) =    RK(IC2)
            REPS(IT2) =  REPS(IC2)
            TTS (IT2) =  TTS(IC2)
            TTS(IT2)  =  TTS(IC2)
 300     CONTINUE
 400  CONTINUE
      ENDIF
      
      IF(MULPHL) THEN
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            DO IPHASE = 1,NPHASES
            PRO(IT1)%RO(IPHASE)   =    PRO(IC1)%RO(IPHASE)
            PRO(IT2)%RO(IPHASE)   =    PRO(IC2)%RO(IPHASE)
            PRO(IT1)%TEMP(IPHASE) =    PRO(IC1)%TEMP(IPHASE)
            PRO(IT2)%TEMP(IPHASE) =    PRO(IC2)%TEMP(IPHASE)
            PRO(IT1)%DTEMP(IPHASE)=    PRO(IC1)%DTEMP(IPHASE)
            PRO(IT2)%DTEMP(IPHASE)=    PRO(IC2)%DTEMP(IPHASE)
            PRO(IT1)%E(IPHASE)    =    PRO(IC1)%E(IPHASE)
            PRO(IT2)%E(IPHASE)    =    PRO(IC2)%E(IPHASE)
            PRO(IT1)%VIS(IPHASE)  =    PRO(IC1)%VIS(IPHASE)
            PRO(IT2)%VIS(IPHASE)  =    PRO(IC2)%VIS(IPHASE)
            PRO(IT1)%CH(IPHASE)   =    PRO(IC1)%CH(IPHASE)
            PRO(IT2)%CH(IPHASE)   =    PRO(IC2)%CH(IPHASE)

            PRO(IT1)%RO(IPHASE)   =    PRO(IC1)%RO(IPHASE)
            PRO(IT2)%RO(IPHASE)   =    PRO(IC2)%RO(IPHASE)
            PRO(IT1)%TEMP(IPHASE) =    PRO(IC1)%TEMP(IPHASE)
            PRO(IT2)%TEMP(IPHASE) =    PRO(IC2)%TEMP(IPHASE)
            PRO(IT1)%DTEMP(IPHASE)=    PRO(IC1)%DTEMP(IPHASE)
            PRO(IT2)%DTEMP(IPHASE)=    PRO(IC2)%DTEMP(IPHASE)
            PRO(IT1)%E(IPHASE)    =    PRO(IC1)%E(IPHASE)
            PRO(IT2)%E(IPHASE)    =    PRO(IC2)%E(IPHASE)
            PRO(IT1)%VIS(IPHASE)  =    PRO(IC1)%VIS(IPHASE)
            PRO(IT2)%VIS(IPHASE)  =    PRO(IC2)%VIS(IPHASE)
            PRO(IT1)%CH(IPHASE)   =    PRO(IC1)%CH(IPHASE)
            PRO(IT2)%CH(IPHASE)   =    PRO(IC2)%CH(IPHASE)

            VAR(IT1)%ALFA(IPHASE) =    VAR(IC1)%ALFA(IPHASE)
            VAR(IT2)%ALFA(IPHASE) =    VAR(IC2)%ALFA(IPHASE)
            VAR(IT1)%X(IPHASE)    =    VAR(IC1)%X(IPHASE)
            VAR(IT2)%X(IPHASE)    =    VAR(IC2)%X(IPHASE)
            VAR(IT1)%EVAPR(IPHASE)=    VAR(IC1)%EVAPR(IPHASE)
            VAR(IT2)%EVAPR(IPHASE)=    VAR(IC2)%EVAPR(IPHASE) ! Mersu
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      IF(TRANSL) THEN ! Intermittency variables
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            TRM(IT1) = TRM(IC1)
            TRM(IT2) = TRM(IC2)
         ENDDO 
      ENDDO
      ENDIF
      
      IF(KSCAL >= 1) THEN
      DO NS = 1,KSCAL
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            FI(IT1,NS) =  FI(IC1,NS)
            FI(IT2,NS) =  FI(IC2,NS)
         ENDDO
      ENDDO
      ENDDO
      ENDIF
      ELSE                      !(ISCL == 0)

      IF(NSCAL >= 1) THEN
      DO NS = 1,NSCAL
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            FI(IT1,NS) =  FI(IC1,NS)
            FI(IT2,NS) =  FI(IC2,NS)
         ENDDO
      ENDDO
      ENDDO
      ENDIF

      ENDIF                     !(ISCL == 0)

      RETURN
      END SUBROUTINE MIRROR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRREY(UU,UV,UW,VV,VW,WW,FUU,FUV,FUW,FVV,FVW,FWW,
     &     NX,NY,NZ,I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IN,JN,KN,IWALL,MAXSB)

C ... This subroutine rotates the stress tensor first into a
C ... coordinate system that has a coordinate direction normal to
C ... the boundary face, next the signs of the cross components perpendicular
C ... to the surface are put to zero and cross component along to the surface
C ... is reflected and last the components are rotated back to the original
C ... coordinate system. PPR 21.12.1994
C ... MAXSB changed to NTOT in the subroutine call. TSii 8.8.2005

      REAL :: UU(*), VV(*), WW(*)
      REAL :: UV(*), UW(*), VW(*)

      REAL :: FUU(*), FVV(*), FWW(*)
      REAL :: FUV(*), FUW(*), FVW(*)

      REAL :: NX(*), NY(*), NZ(*)
      
      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF
      INTEGER :: IT,IT1,IT2,JL,IL,KB,KB1,KB2

      DO I = 1,MAXSB
         FUU(I) = UU(I)
         FUV(I) = UV(I)
         FUW(I) = UW(I)
         FVV(I) = VV(I)
         FVW(I) = VW(I)
         FWW(I) = WW(I)
      ENDDO

      DO I = 1,MAXSB
         UU(I) = 0.
         UV(I) = 0.
         WW(I) = 0.
         UW(I) = 0.
         VV(I) = 0.
         VW(I) = 0.
      ENDDO
      
      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine mirrey !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR
      KT1 = (KW          - 1 + KN)*KSTR
      KT2 = (KW -   KDIR - 1 + KN)*KSTR

      DO 100 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            U11 =  FUU(IC1)
            U12 = -FUV(IC1)
            U13 = -FUW(IC1)
            U22 =  FVV(IC1)
            U23 = -FVW(IC1)
            U33 =  FWW(IC1)
            U21 =  U12
            U31 =  U13
            U32 =  U23

            V11 =  FUU(IC2)
            V12 = -FUV(IC2)
            V13 = -FUW(IC2)
            V22 =  FVV(IC2)
            V23 = -FVW(IC2)
            V33 =  FWW(IC2)
            V21 =  V12
            V31 =  V13
            V32 =  V23

            A11 = -2.*NX(IT)*NX(IT) + 1.
            A12 = -2.*NX(IT)*NY(IT)
            A13 = -2.*NX(IT)*NZ(IT)
            A21 = -2.*NX(IT)*NY(IT)
            A22 = -2.*NY(IT)*NY(IT) + 1.
            A23 = -2.*NY(IT)*NZ(IT)
            A31 = -2.*NX(IT)*NZ(IT)
            A32 = -2.*NY(IT)*NZ(IT)
            A33 = -2.*NZ(IT)*NZ(IT) + 1.

            B11 = A11*U11 + A12*U21 + A13*U31
            B12 = A11*U12 + A12*U22 + A13*U32
            B13 = A11*U13 + A12*U23 + A13*U33
            B21 = A21*U11 + A22*U21 + A23*U31
            B22 = A21*U12 + A22*U22 + A23*U32
            B23 = A21*U13 + A22*U23 + A23*U33
            B31 = A31*U11 + A32*U21 + A33*U31
            B32 = A31*U12 + A32*U22 + A33*U32
            B33 = A31*U13 + A32*U23 + A33*U33
                  
            F11 = A11*V11 + A12*V21 + A13*V31
            F12 = A11*V12 + A12*V22 + A13*V32
            F13 = A11*V13 + A12*V23 + A13*V33
            F21 = A21*V11 + A22*V21 + A23*V31
            F22 = A21*V12 + A22*V22 + A23*V32
            F23 = A21*V13 + A22*V23 + A23*V33
            F31 = A31*V11 + A32*V21 + A33*V31
            F32 = A31*V12 + A32*V22 + A33*V32
            F33 = A31*V13 + A32*V23 + A33*V33

            C11 = B11*A11 + B12*A21 + B13*A31
            C12 = B11*A12 + B12*A22 + B13*A32
            C13 = B11*A13 + B12*A23 + B13*A33
            C21 = B21*A11 + B22*A21 + B23*A31
            C22 = B21*A12 + B22*A22 + B23*A32
            C23 = B21*A13 + B22*A23 + B23*A33
            C31 = B31*A11 + B32*A21 + B33*A31
            C32 = B31*A12 + B32*A22 + B33*A32
            C33 = B31*A13 + B32*A23 + B33*A33
                                             
            G11 = F11*A11 + F12*A21 + F13*A31
            G12 = F11*A12 + F12*A22 + F13*A32
            G13 = F11*A13 + F12*A23 + F13*A33
            G21 = F21*A11 + F22*A21 + F23*A31
            G22 = F21*A12 + F22*A22 + F23*A32
            G23 = F21*A13 + F22*A23 + F23*A33
            G31 = F31*A11 + F32*A21 + F33*A31
            G32 = F31*A12 + F32*A22 + F33*A32
            G33 = F31*A13 + F32*A23 + F33*A33

 111        format(1x,a18,1X,3E15.7,1x,2i4)
            FUU(IT1) =  C11
            FUV(IT1) = -C12
            FUW(IT1) = -C13
            FVV(IT1) =  C22
            FVW(IT1) = -C23
            FWW(IT1) =  C33

            FUU(IT2) =  G11
            FUV(IT2) = -G12
            FUW(IT2) = -G13
            FVV(IT2) =  G22
            FVW(IT2) = -G23
            FWW(IT2) =  G33
 200     CONTINUE
 100  CONTINUE

      DO I = 1,MAXSB
          UU(I) = FUU(I) 
          UV(I) = FUV(I)  
          UW(I) = FUW(I)  
          VV(I) = FVV(I)  
          VW(I) = FVW(I)  
          WW(I) = FWW(I)  
      ENDDO
      
      RETURN
      END SUBROUTINE MIRREY
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRGRI(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A1,NX,NY,NZ,
     &             I1,I2,J1,J2,KW,ISTRID,JSTRID,KSTRID,IN,JN,KN,IWALL)

C ... This subroutine moves the grid center points components 
C ... first into a coordinate system that has a zero in the boundary and 
C ... coordinate direction normal to the boundary face, 
C ... next the signs of the components perpendicular
C ... to the surface are changed 
C ... and then the components are rotated back to the original
C ... coordinate system. All this is accomplished by one matrix
C ... multiplication. After that center points of the ghost shells can 
C ... be obtained.

      REAL :: VOL(*), D1(*), D2(*), D3(*), DISTW(*), BLANK(*),
     +        A1(*),NX(*),NY(*),NZ(*)
      REAL :: UMIR1, UMIR2, VMIR1, VMIR2, WMIR1, WMIR2

      REAL :: XC(*), YC(*), ZC(*), 
     +        XMIR1,YMIR1,ZMIR1,XMIR2,YMIR2,ZMIR2

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF
      INTEGER :: KB,KB1,KB2,IT,IT1,IT2,JL,IL

      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine mirgri !'
      ENDIF

      KB   = (KW +   KOFF - 1 + KN)*KSTRID  ! WALL VALUE FOR GEOMETRIC
      KB1  = (KW +   KDIR - 1 + KN)*KSTRID  ! GHOST CELL 1            
      KB2  = (KW + 2*KDIR - 1 + KN)*KSTRID  ! GHOST CELL 2            
      KT1  = (KW          - 1 + KN)*KSTRID  ! COMPUTATIONAL DOMAIN 1  
      KT2  = (KW -   KDIR - 1 + KN)*KSTRID  ! COMPUTATIONAL DOMAIN 2  
      DO 100 J = J1,J2
         JL    = (J-1+JN)*JSTRID
         DO 200 I = I1,I2
            IL    = JL  + (I-1+IN)*ISTRID
            IT    = 1 + IL + KB    ! WALL VALUE FOR GEOMETRIC
            IT1   = 1 + IL + KB1   ! GHOST CELL 1            
            IT2   = 1 + IL + KB2   ! GHOST CELL 2            
            IC1   = 1 + IL + KT1   ! COMPUTATIONAL DOMAIN 1  
            IC2   = 1 + IL + KT2   ! COMPUTATIONAL DOMAIN 2  

            A11   = -2.*NX(IT)*NX(IT) + 1.0
            A12   = -2.*NX(IT)*NY(IT)
            A13   = -2.*NX(IT)*NZ(IT)
            A21   = -2.*NX(IT)*NY(IT)
            A22   = -2.*NY(IT)*NY(IT) + 1.0
            A23   = -2.*NY(IT)*NZ(IT)
            A31   = -2.*NX(IT)*NZ(IT)
            A32   = -2.*NY(IT)*NZ(IT)
            A33   = -2.*NZ(IT)*NZ(IT) + 1.0

            XAUX   = .5
            XC0   =  XC(IC1) + XAUX*(XC(IC1)-XC(IC2))
            YC0   =  YC(IC1) + XAUX*(YC(IC1)-YC(IC2))
            ZC0   =  ZC(IC1) + XAUX*(ZC(IC1)-ZC(IC2))
C ... TO THE ZERO CO-ORDINATE
            XMIR1 = XC(IC1) - XC0
            YMIR1 = YC(IC1) - YC0
            ZMIR1 = ZC(IC1) - ZC0

            XMIR2 = XC(IC2) - XC0
            YMIR2 = YC(IC2) - YC0
            ZMIR2 = ZC(IC2) - ZC0

            XC(IT1) = A11*XMIR1 + A12*YMIR1 + A13*ZMIR1 + XC0
            YC(IT1) = A21*XMIR1 + A22*YMIR1 + A23*ZMIR1 + YC0
            ZC(IT1) = A31*XMIR1 + A32*YMIR1 + A33*ZMIR1 + ZC0

            XC(IT2) = A11*XMIR2 + A12*YMIR2 + A13*ZMIR2 + XC0
            YC(IT2) = A21*XMIR2 + A22*YMIR2 + A23*ZMIR2 + YC0
            ZC(IT2) = A31*XMIR2 + A32*YMIR2 + A33*ZMIR2 + ZC0

            VOL(IT1) = VOL(IC1)
            D1(IT1)  = D1(IC1)
            D2(IT1)  = D2(IC1)
            D3(IT1)  = D3(IC1)
            DISTW(IT1)  = DISTW(IC1)
            BLANK(IT1)  = BLANK(IC1)

            VOL(IT2) = VOL(IC2)
            D1(IT2)  = D1(IC2)
            D2(IT2)  = D2(IC2)
            D3(IT2)  = D3(IC2)
            BLANK(IT2)  = BLANK(IC2)

 200     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE MIRGRI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRGRC(XCO,YCO,ZCO,I1,I2,J1,J2,KW,
     &                  ISTRID,JSTRID,KSTRID,IN,JN,KN,IWALL)

C ... This subroutine moves the grid corner points components 
C ... first into a coordinate system that has a zero in the boundary and 
C ... coordinate direction normal to the boundary face, 
C ... next the signs of the components perpendicular
C ... to the surface are changed 
C ... and then the components are rotated back to the original
C ... coordinate system. All this is accomplished by one matrix
C ... multiplication. After that the corner points of the ghost shells can 
C ... be obtained.

      REAL :: XCO(*), YCO(*), ZCO(*)
      REAL :: NX,NY,NZ,X1,X2,Y1,Y2,Z1,Z2,PALPO,
     +        A11,A12,A13,A21,A22,A23,A31,A32,A33,XMIR1,YMIR1,ZMIR1,
     +        XCIT,YCIT,ZCIT,EPS

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF
      INTEGER :: KB,KB1,KB2,IT,IT1,IT2,JL,IL

      EPS  = 1.0E-30

      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine mirgrc !'
      ENDIF

      KCC  = (KW +   KOFF - 1 + KN - KDIR)*KSTRID
      KC   = (KW           -1 + KN)*KSTRID
      KB   = (KW +   KOFF - 1 + KN)*KSTRID
      KB1  = (KW +   KDIR + KOFF - 1 + KN)*KSTRID
      DO 100 J = J1,J2
         JL    = (J-1+JN)*JSTRID
         DO 200 I = I1,I2

            IL    = JL  + (I-1+IN)*ISTRID
            IC    = 1 + IL  + KC
            ICC   = 1 + IL  + KCC
            IT    = 1 + IL  + KB
            IT1   = 1 + IL  + KB1
            II1   = IT + ISTRID
            II2   = IT - ISTRID
            II3   = IT + JSTRID
            II4   = IT - JSTRID

            X1    = XCO(II1) - XCO(II2)
            Y1    = YCO(II1) - YCO(II2)
            Z1    = ZCO(II1) - ZCO(II2)

            X2    = XCO(II3) - XCO(II4)
            Y2    = YCO(II3) - YCO(II4)
            Z2    = ZCO(II3) - ZCO(II4)

            NX    = Y1*Z2 - Z1*Y2
            NY    = Z1*X2 - X1*Z2
            NZ    = X1*Y2 - Y1*X2
            PALPO = SQRT(NX**2 + NY**2 + NZ**2 + EPS)
            NX    = NX/PALPO
            NY    = NY/PALPO
            NZ    = NZ/PALPO

c            call ijkpai(it1,48,32,1,mmm,lll,nnn)
c            write(77,111) i,j,mmm,lll,nnn,x1,y1,z1,palpo,x2,y2,z2
c 111        format(5i4,7E10.2)

            A11   = -2.*NX*NX + 1.0
            A12   = -2.*NX*NY
            A13   = -2.*NX*NZ
            A21   = -2.*NX*NY
            A22   = -2.*NY*NY + 1.0
            A23   = -2.*NY*NZ
            A31   = -2.*NX*NZ
            A32   = -2.*NY*NZ
            A33   = -2.*NZ*NZ + 1.0

C ... TO THE ZERO CO-ORDINATE
            XMIR1 = XCO(IT) - XCO(ICC)
            YMIR1 = YCO(IT) - YCO(ICC)
            ZMIR1 = ZCO(IT) - ZCO(ICC)

            XCIT = A11*XMIR1 + A12*YMIR1 + A13*ZMIR1
            YCIT = A21*XMIR1 + A22*YMIR1 + A23*ZMIR1
            ZCIT = A31*XMIR1 + A32*YMIR1 + A33*ZMIR1

C ... TO NORMAL COORDINATE
            XCO(IT1) = XCO(IT) - XCIT
            YCO(IT1) = YCO(IT) - YCIT
            ZCO(IT1) = ZCO(IT) - ZCIT
 200     CONTINUE
 100  CONTINUE

 900  continue

      RETURN
      END SUBROUTINE MIRGRC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,BIJ,
     +     KSCAL,ITURB,MAXEB,MAXSB,NX,NY,NZ,WAVEH,IPL,IHF,I1,I2,J1,J2,
     +     KW,ISTR,JSTR,KSTR,IN,JN,KN,IWALL,XVEL,FRSPRE,FRSDEN,
     +     ISTRES,POLD,DSURLE,IWAVEB,NGL)

C ... FREE SURFACE

C ... This subroutine rotates the vector components first into a
C ... coordinate system that has a coordinate direction normal to
C ... the boundary face, next the signs of the components perpendicular
C ... to the surface are changed (hence imposing no flow through the
C ... face) and last the components are rotated back to the original
C ... coordinate system. All this is accomplished by one matrix
C ... multiplication.

      USE NS3CO, ONLY : LN, GVEX, GVEY, GVEZ, G0
      USE MAIN_ARRAYS, ONLY : UROTCP, VROTCP, WROTCP
      
      IMPLICIT NONE

      INTEGER :: IN, JN, KN

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF,KB,
     +           KB1,KB2,IT,IT1,IT2,JL,IL,MAXEB,KT3,KSCAL,ITURB,MAXSB,
     +           IPL,ISTRES,I,J,IC1,IC2,IC3,KT1,KT2,NN,NP,NPP,NS,NGL

      INTEGER :: IHF(*),IWAVEB(*)

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),NX(*), NY(*), NZ(*),
     +        BIJ(MAXEB,*),POLD(*),DSURLE(*)
      REAL :: PDPAT((I2-I1+1+2*LN)*(J2-J1+1+2*LN))
      REAL :: DELTAP,UMIR1,UMIR2,VMIR1,VMIR2,WMIR1,WMIR2
      REAL :: FRSPRE,FRSDEN,XVEL,A11,A12,A13,A21,A22,A23,A31,A32,A33

      REAL :: WAVEH(*)


      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine FRESUR !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR
      KT1 = (KW          - 1 + KN)*KSTR
      KT2 = (KW -   KDIR - 1 + KN)*KSTR


      DO 100 J=J1,J2
*         NN = (JN+J-J1)*(I2-I1+1+2*IN) - I1 + IN + IHF(IPL)
         NN = (LN+J-J1)*(I2-I1+1+2*LN) - I1 + LN + IHF(IPL)
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2
            NP  = NN + I                     ! Patch index with ghost cells
*            NP  = NN + I                    ! Patch index with LN ghost cells

            NPP = NN + I - IHF(IPL)          ! Local tempperature index


            A11 = XVEL*(-2.*NX(IT)*NX(IT) + 1.0)
            A12 = XVEL*(-2.*NX(IT)*NY(IT))
            A13 = XVEL*(-2.*NX(IT)*NZ(IT))
            A21 = XVEL*(-2.*NX(IT)*NY(IT))
            A22 = XVEL*(-2.*NY(IT)*NY(IT) + 1.0)
            A23 = XVEL*(-2.*NY(IT)*NZ(IT))
            A31 = XVEL*(-2.*NX(IT)*NZ(IT))
            A32 = XVEL*(-2.*NY(IT)*NZ(IT))
            A33 = XVEL*(-2.*NZ(IT)*NZ(IT) + 1.0)

            UMIR1 = U(IC1) - GVEX - UROTCP(IC1)
            VMIR1 = V(IC1) - GVEY - VROTCP(IC1)
            WMIR1 = W(IC1) - GVEZ - WROTCP(IC1)

            UMIR2 = U(IC2) - GVEX - UROTCP(IC2)
            VMIR2 = V(IC2) - GVEY - VROTCP(IC2)
            WMIR2 = W(IC2) - GVEZ - WROTCP(IC2)
  
            U(IT1) = A11*UMIR1 + A12*VMIR1 + A13*WMIR1 +GVEX+UROTCP(IT1)
            V(IT1) = A21*UMIR1 + A22*VMIR1 + A23*WMIR1 +GVEY+VROTCP(IT1)
            W(IT1) = A31*UMIR1 + A32*VMIR1 + A33*WMIR1 +GVEZ+WROTCP(IT1)

            U(IT2) = A11*UMIR2 + A12*VMIR2 + A13*WMIR2 +GVEX+UROTCP(IT2)
            V(IT2) = A21*UMIR2 + A22*VMIR2 + A23*WMIR2 +GVEY+VROTCP(IT2)
            W(IT2) = A31*UMIR2 + A32*VMIR2 + A33*WMIR2 +GVEZ+WROTCP(IT2)

c           apu1 = nx(it)*u(it1)+ny(it)*v(it1)+nz(it)*w(it1)
c           apu2 = nx(it)*u(ic1)+ny(it)*v(ic1)+nz(it)*w(ic1)
C ... testi taysi pelaus
c            U(IT1) = UMIR1
c            V(IT1) = VMIR1
c            W(IT1) = WMIR1

c            U(IT2) = UMIR2
c            V(IT2) = VMIR2 
c            W(IT2) = WMIR2

            RO(IT1)    =    RO(IC1)
            E (IT1)    =     E(IC1)
            EPS2(IT1)  =  EPS2(IC1)
            VIST(IT1)  =  VIST(IC1)

            RO   (IT2) =    RO(IC2)
            E    (IT2) =     E(IC2)
            EPS2 (IT2) =  EPS2(IC2)
            VIST (IT2) =  VIST(IC2)

C ... THESE ARE SPECIAL TREATMENT FOR FREE SURFACE
            P(IT1)     = FRSPRE
            P(IT2)     = FRSPRE

c            DELTAP = FRSDEN*G0*DSURLE(NP) !WAVEH(NP)
            DELTAP     = FRSDEN*G0*WAVEH(NP)
            PDIFF(IT1) = 2.*DELTAP - PDIFF(IC1)!+1./6.*PDIFF(IC2)
c     +                 +  PDIFF(IT1))*.5 ! Tryied to under-relax
            PDPAT(NPP) = PDIFF(IC1)-POLD(IC1)
c            PDIFF(IT1) = 2.*DELTAP -1.*(PDIFF(IC1)-POLD(IC1))
c     +      +PDIFF(IT1)
       
c            PDIFF(IT1) = .5*(2.*DELTAP - PDIFF(IC1)+PDIFF(IT1))
c            pdiff(it1) = 2.*pdiff(ic1) - pdiff(ic2)
c            u(it1) = 2.*u(ic1) - u(ic2)
c            v(it1) = 2.*v(ic1) - v(ic2)
c            w(it1) = 2.*w(ic1) - w(ic2)
! taysi oli liikaa??
            PDIFF(IT2) = 2.*PDIFF(IT1)-PDIFF(IC1) ! + PDIFF(IT2))*.5
c             PDIFF(IT2) = 3.*PDIFF(IT1)-3.*PDIFF(IC1)+PDIFF(IC2)
c            u(IT2) = 2.*u(IT1)-u(IC1)
c            v(IT2) = 2.*v(IT1)-v(IC1)
c            w(IT2) = 2.*w(IT1)-w(IC1)
*           PDIFF(IT2) = 4.*DELTAP - PDIFF(IC2) ! testaa
*         write(571+ngl,*) i,j,ic1,it1,pdiff(it1),pdiff(ic1),
*     +  pdiff(it2),real(waveh(np))
c       write(568,*) i,j,np,ic1,deltap,waveh(np),pdiff(ic1),
c     + apu1,apu2
       
             IWAVEB(IT1) = 13

 200     CONTINUE
 100  CONTINUE

c      DO J=J1,J2
c         NN = (JN+J-J1)*(I2-I1+1+2*IN) - I1 + IN + IHF(IPL)
c         JL = (J-1+JN)*JSTR
c         DO I=I1,I2
c            IL  = JL + (I-1+IN)*ISTR
c            IT  = 1 + IL + KB
c            IT1 = 1 + IL + KB1
c            IT2 = 1 + IL + KB2
c            IC1 = 1 + IL + KT1
c            IC2 = 1 + IL + KT2
c            NP  = NN + I + 1          ! Patch tempperature index
c            NPP = NN + I - IHF(IPL)! Local tempperature indexccccc
c
c            PDIFF(IT1) = PDIFF(IT1) + 1.*    !FRSDEN*G0*
c     +      (4.*PDPAT(NPP)-PDPAT(NPP-1)-PDPAT(NPP+1)
c     +       -PDPAT(NPP-(I2-I1+1+2*IN))-PDPAT(NPP+(I2-I1+1+2*IN)))
c     +      (4.*WAVEH(NP)-WAVEH(NP-1)-WAVEH(NP+1)
c     +       -WAVEH(NP-(I2-I1+1+2*IN))-WAVEH(NP+(I2-I1+1+2*IN)))
c            PDIFF(IT2) = 2.*PDIFF(IT1)-PDIFF(IC1)
c         ENDDO
c      ENDDO

      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      DO 400 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 300 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            RK  (IT1) =    RK(IC1)
            REPS(IT1) =  REPS(IC1)
            RK  (IT2) =    RK(IC2)
            REPS(IT2) =  REPS(IC2)
 300     CONTINUE
 400  CONTINUE
      ENDIF

      IF(KSCAL >= 1) THEN
      DO NS = 1,KSCAL
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            FI(IT1,NS) =  FI(IC1,NS)
            FI(IT2,NS) =  FI(IC2,NS)
         ENDDO
      ENDDO
      ENDDO
      ENDIF

      IF(ISTRES >= 1) THEN
      DO NS = 1,6
      DO J=J1,J2
         JL = (J-1+JN)*JSTR
         DO I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2
            IC1 = 1 + IL + KT1
            IC2 = 1 + IL + KT2

            BIJ(IT1,NS) =  BIJ(IC1,NS)
            BIJ(IT2,NS) =  BIJ(IC2,NS)
         ENDDO
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FRESUR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERMUD(IND,R,Q,IMAX,JMAX,IN,JN,KSLBS)

C ... Driver for the Permutation Subroutines

      REAL :: R(*), Q(*)
      INTEGER :: IMAX,JMAX,IN,JN,I,J,K,IR,IQ,KSLBS,IND

      IF(MOD(IND,4) == 0) RETURN
      IF(MOD(IND,4) == 1) CALL PERD1(R,Q,JMAX,IMAX,IN,JN,KSLBS)
      IF(MOD(IND,4) == 2) CALL PERD2(R,Q,IMAX,JMAX,IN,JN,KSLBS)
      IF(MOD(IND,4) == 3) CALL PERD3(R,Q,JMAX,IMAX,IN,JN,KSLBS)

      RETURN
      END SUBROUTINE PERMUD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERMUT(IND,R,Q,IMAX,JMAX,IN,JN,KSLBS)

C ... Driver for the Permutation Subroutines

      REAL    :: R(*), Q(*)
      INTEGER :: IMAX,JMAX,IN,JN,I,J,K,IR,IQ,KSLBS,IND

      IF (MOD(IND,4) == 0) RETURN
      IF (MOD(IND,4) == 1) CALL PERM1(R,Q,JMAX,IMAX,IN,JN,KSLBS)
      IF (MOD(IND,4) == 2) CALL PERM2(R,Q,IMAX,JMAX,IN,JN,KSLBS)
      IF (MOD(IND,4) == 3) CALL PERM3(R,Q,JMAX,IMAX,IN,JN,KSLBS)

      RETURN
      END SUBROUTINE PERMUT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERM1(R,Q,IMAX,JMAX,IN,JN,KSLBS)

C     Permuts data       1234           482
C     in the following   5678    -->    371
C     way:               9012           260
C                                       159
C                       Im x Jm       Jm x Im
C     I' = J
C     J' = IMAX+2*IN-I+1
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
               IQ = (ISTR-I)*JSTR+J+(K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE

      RETURN
      END SUBROUTINE PERM1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERM2(R,Q,IMAX,JMAX,IN,JN,KSLBS)

C     Permuts data       1234          2109
C     in the following   5678    -->   8765
C     way:               9012          4321
C                       Im x Jm       Im x Jm
C     I' = IMAX+2*IN-I+1
C     J' = JMAX+2*JN-J+1
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
               IQ = 1 + (JSTR-J)*ISTR+(ISTR-I)+(K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE

      RETURN
      END SUBROUTINE PERM2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERM3(R,Q,IMAX,JMAX,IN,JN,KSLBS)

C     Permuts data       1234           951
C     in the following   5678    -->    062
C     way:               9012           173
C                                       284
C                       Im x Jm       Jm x Im
C     I' = JMAX+2*JN-J+1
C     J' = I
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
C               IQ = (ISTR-I)*JSTR+J+(K-1)*ISTR*JSTR
               IQ = 1 + (JSTR-J) + (I-1)*JSTR + (K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE

      RETURN
      END SUBROUTINE PERM3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
      SUBROUTINE PERD1(R,Q,IMAX,JMAX,IN,JN,KSLBS)
C     ****************
C     Permuts data       1234           482
C     in the following   5678    -->    371
C     way:               9012           260
C                                       159
C                       Im x Jm       Jm x Im
C     I' = J
C     J' = IMAX+2*IN-I+1
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
               IQ = (ISTR-I)*JSTR+J+(K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE

      RETURN
      END SUBROUTINE PERD1
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERD2(R,Q,IMAX,JMAX,IN,JN,KSLBS)

C     Permuts data       1234          2109
C     in the following   5678    -->   8765
C     way:               9012          4321
C                       Im x Jm       Im x Jm
C     I' = IMAX+2*IN-I+1
C     J' = JMAX+2*JN-J+1
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
               IQ = 1 + (JSTR-J)*ISTR+(ISTR-I)+(K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE
      RETURN
      END SUBROUTINE PERD2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PERD3(R,Q,IMAX,JMAX,IN,JN,KSLBS)

C     Permuts data       1234           951
C     in the following   5678    -->    062
C     way:               9012           173
C                                       284
C                       Im x Jm       Jm x Im
C     I' = JMAX+2*JN-J+1
C     J' = I
C
      REAL :: R(*), Q(*)
      INTEGER :: IMAX, JMAX, IN, JN, I, J, K, IR, IQ, KSLBS

      ISTR = IMAX+2*IN
      JSTR = JMAX+2*JN
      DO 100 K=1,KSLBS
         DO 200 J=1,JMAX+2*JN
            DO 300 I=1,IMAX+2*IN
C               IQ = (ISTR-I)*JSTR+J+(K-1)*ISTR*JSTR
               IQ = 1 + (JSTR-J) + (I-1)*JSTR + (K-1)*ISTR*JSTR
               IR = 1 + (I-1) + (J-1)*ISTR + (K-1)*ISTR*JSTR
               Q(IQ) = R(IR)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 400 I=1,ISTR*JSTR*KSLBS
         R(I) = Q(I)
 400  CONTINUE

      RETURN
      END SUBROUTINE PERD3
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPORT(R,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN)
C     *****************
C     Indices on various walls are as follows
C     Wall No 1 -->    K,J,-I
C     Wall No 2 -->    I,K,-J
C     Wall No 3 -->    J,I,-K
C     Wall No 4 -->    J,K,I
C     Wall No 5 -->    K,I,J
C     Wall No 6 -->    I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL    :: R(*), Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR

      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF
      DO 100 K=K1,K2,KSTP
         DO 200 J=J1,J2
            DO 300 I=I1,I2
               IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IT = (K-K1)*KSTP*(J2-J1+1)*(I2-I1+1)+
     &              (J-J1)*(I2-I1+1)+(I-I1)+1
               Q(IT) = R(IL)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE EXPORT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPORD(R,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN)

C     Indices on various walls are as follows
C     Wall No 1 -->    K,J,-I
C     Wall No 2 -->    I,K,-J
C     Wall No 3 -->    J,I,-K
C     Wall No 4 -->    J,K,I
C     Wall No 5 -->    K,I,J
C     Wall No 6 -->    I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL :: R(*), Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR

      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF
      DO 100 K=K1,K2,KSTP
         DO 200 J=J1,J2
            DO 300 I=I1,I2
               IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IT = (K-K1)*KSTP*(J2-J1+1)*(I2-I1+1)+
     &              (J-J1)*(I2-I1+1)+(I-I1)+1
               Q(IT) = R(IL)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE EXPORD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPORT(R,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN)

C     K1 = K-value at the First Slab Outside the Block
C     K2 = K-value at the Second Slab Outside the Block
C     Note the inversion of the IJ-Base to account for the transposed
C     indices in the incoming data.
C     Indices on various walls are as follows
C     Wall No 1 -->    K,J,-I
C     Wall No 2 -->    I,K,-J
C     Wall No 3 -->    J,I,-K
C     Wall No 4 -->    J,K,I
C     Wall No 5 -->    K,I,J
C     Wall No 6 -->    I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL    :: R(*), Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR

      IF (K1 > K2) THEN
         KSTP = -1
      ELSE
         KSTP = 1
      ENDIF
      DO 100 K=K1,K2,KSTP
         DO 200 J=J1,J2
            DO 300 I=I1,I2
               IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IT = (K-K1)*KSTP*(J2-J1+1)*(I2-I1+1)+
     &              (I-I1)*(J2-J1+1)+(J-J1) + 1
               R(IL) = Q(IT)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END SUBROUTINE IMPORT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPORD(R,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN)

C     K1 = K-value at the First Slab Outside the Block
C     K2 = K-value at the Second Slab Outside the Block
C     Note the inversion of the IJ-Base to account for the transposed
C     indices in the incoming data.
C     Indices on various walls are as follows
C     Wall No 1 --> K,J,-I
C     Wall No 2 --> I,K,-J
C     Wall No 3 --> J,I,-K
C     Wall No 4 --> J,K,I
C     Wall No 5 --> K,I,J
C     Wall No 6 --> I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL :: R(*), Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR

      IF (K1 > K2) THEN
         KSTP = -1
      ELSE
         KSTP = 1
      ENDIF
      DO 100 K=K1,K2,KSTP
         DO 200 J=J1,J2
            DO 300 I=I1,I2
               IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IT = (K-K1)*KSTP*(J2-J1+1)*(I2-I1+1)+
     &              (I-I1)*(J2-J1+1)+(J-J1) + 1
               R(IL) = Q(IT)
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE IMPORD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE APPRAI(XC,YC,ZC,ZZZ,APU,APU2,JLOC,JTRA,APP,
     +     ICON,JAPP,IWAPP,NPATCH,IW,IMAX,JMAX,
     +     KMAX,ITAG,IN,JN,KN,NPPV,IPRO,N,M,MAXB)

      USE NS3CO, ONLY : IC9
      
      REAL :: APP(*), APU(*)
      REAL :: XC(*), YC(*), ZC(*), APU2(*), ZZZ(*)
      INTEGER :: IW(NPPV,*), ITAG(*), JAPP(*), IWAPP(NPPV,*)
      INTEGER :: ICON(IC9,*), JLOC(*), JTRA(*)
      INTEGER :: L, M, N
      INTEGER :: IA1, IM1, JA1, JM1, IA2, IM2, JA2, JM2
      LOGICAL :: CONN

C ... Appraise the fluxes over block connection in case of non-matching
C ... connection.  

C ... JLOC  is local grid point number
C ... JTRA  is connective grid point number
C ... APP   is weighting of the JTRA to the JLOC
C ... JAPP  is number of connectivities
C ... IWAPP is starting number for this patch in arrays JLOC,JTRA & APP

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR


      DO 1000 L = 1,NPATCH

      IF(ITAG(L) == 1) THEN

      ITYPE   = ICON(1,L)  ! Boundary condition type

      CONN    = (ITYPE == 1 .AND. ICON(22,L) /= 0) 
      
      IF(CONN) THEN

      IGP  = ICON(2,L)   ! Process local patch (or block) number
      LV   = ICON(8,L)   ! number of connective block (process local)
      LP   = ICON(15,L)  ! connective patch n. (block l.)
      JP   = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO = ICON(17,L)  ! number of connectivity process
      IOTY = ICON(18,L)  ! connectivity method (1=normal,3=MPI)
      NBG  = ICON(24,L)

      IW1  = IW(LP,LV)   ! IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2  = IW(L,N)     ! IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI

      IWA  = IWAPP(L,N)  ! Starting address of the transfer matrises

      CALL PATCHE(ICON,L,M,IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2)  ! Extend

      IA  = IA1
      IM  = IM1 + 1
      JA  = JA1
      JM  = JM1 + 1
      IMA = IM2 - IA2 + 2
      JMA = JM2 - JA2 + 2

      KA     = 1
      KM     = 1
      IOTOT  = IABS(IMA*JMA)
      IOSIZE = MAX0(IOTOT,(IM - IA + 1)*(JM - JA + 1))

      IF(IOTY <= 1) CALL SETVDD(APU2,ZZZ(IW1),2*3*IOTOT)
      IF(IOTY >= 2) CALL SETVDD(APU2,ZZZ(IW2),2*3*IOTOT)

C ... RECEIVE FROM NODE ICON(8,L) DATA TO ZZZ-ARRAY WITH DISTR.MEM.

      MAX2  = MAXB/2
      IN1 = 1
      IF(ICON(3,L) == 1) THEN                          ! WALL NO. 1
         CALL APPRAX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,KN,JN,IN)

      ELSEIF(ICON(3,L) == 2) THEN                      ! WALL NO. 2
         CALL APPRAX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,IN,KN,JN)
         
      ELSEIF(ICON(3,L) == 3) THEN                      ! WALL NO. 3
         CALL APPRAX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,JN,IN,KN)

      ELSEIF(ICON(3,L) == 4) THEN                      ! WALL NO. 4
         KA = IMAX + 1
         KM = IMAX + 1 
         CALL APPRAX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,JN,KN,IN)

      ELSEIF(ICON(3,L) == 5) THEN                      ! WALL NO. 5
         KA = JMAX + 1
         KM = JMAX + 1
         CALL APPRAX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,KN,IN,JN)

      ELSEIF(ICON(3,L) == 6) THEN                      ! WALL NO. 6
         KA = KMAX + 1
         KM = KMAX + 1
         CALL APPRAX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,
     +        JLOC(IWA),JTRA(IWA),APP(IWA),JAPP(L),APU,APU(MAX2),MAX2,
     +               ICON( 9,L),IA2,IM2,JA2,JM2,ICON(14,L),
     +               IOTOT,IOSIZE,NBG,L,NBG,IWA,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such wall. Exiting ...', ICON(3,L),L
         STOP
      ENDIF

      DO I = L+1,NPPV
         IWAPP(I,N) = IWAPP(L,N) + JAPP(L)
      ENDDO

      IWAPP(1,N+1) = IWAPP(L,N) + JAPP(L)

      IF(IWAPP(1,N+1) >= MAXB) WRITE(*,333) MAXB,IWAPP(1,N+1)
 333  FORMAT(' APPRAI: there is problem with size of MAXB (=',I7,').'/
     +       '         IWAPPMAX is too large',I7)
      ENDIF                                   ! END OF CONNECTIONS

      ENDIF

1000  CONTINUE

      RETURN
      END SUBROUTINE APPRAI
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPOXX(XC,YC,ZC,XX,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,
     +                  IN,JN,KN,D)

C     K1 = K-value at the First Slab Outside the Block
C     K2 = K-value at the Second Slab Outside the Block
C     Note the inversion of the IJ-Base to account for the transposed
C     indices in the incoming data.
C     Indices on various walls are as follows
C     Wall No 1 -->    K,J,-I
C     Wall No 2 -->    I,K,-J
C     Wall No 3 -->    J,I,-K
C     Wall No 4 -->    J,K,I
C     Wall No 5 -->    K,I,J
C     Wall No 6 -->    I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL    :: XC(*), YC(*), ZC(*), XX(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR
      REAL    :: DIST1, DIST

      IOTOT   = ABS((I2-I1+1)*(J2-J1+1))
          
C ... Calculate reference lenght

      KAVE = (K1 + K2)/2
      JAVE = (J1 + J2)/2
      IAVE = (I1 + I2)/2
      ICT  = (KAVE-1+KN)*KSTR+(JAVE-1+JN)*JSTR+(IAVE-1+IN)*ISTR+1

      I1J1 = (K1-1+KN)*KSTR+(J1-1+JN)*JSTR+(I1-1+IN)*ISTR+1
      I2J1 = (K1-1+KN)*KSTR+(J1-1+JN)*JSTR+(I2-1+IN)*ISTR+1
      I1J2 = (K1-1+KN)*KSTR+(J2-1+JN)*JSTR+(I1-1+IN)*ISTR+1
      I2J2 = (K1-1+KN)*KSTR+(J2-1+JN)*JSTR+(I2-1+IN)*ISTR+1

      DIST11 = SQRT((XC(ICT) - XC(I1J1))**2 + 
     +     (YC(ICT) - YC(I1J1))**2 + (ZC(ICT) - ZC(I1J1))**2)
      DIST12 = SQRT((XC(ICT) - XC(I1J2))**2 + 
     +     (YC(ICT) - YC(I1J2))**2 + (ZC(ICT) - ZC(I1J2))**2)
      DIST21 = SQRT((XC(ICT) - XC(I2J1))**2 + 
     +     (YC(ICT) - YC(I2J1))**2 + (ZC(ICT) - ZC(I2J1))**2)
      DIST22 = SQRT((XC(ICT) - XC(I2J2))**2 + 
     +     (YC(ICT) - YC(I2J2))**2 + (ZC(ICT) - ZC(I2J2))**2)

      REFLEN = DIST11 + DIST12 + DIST21 + DIST22

      DIST = 0.
      K=K1

      DO 200 J=J1,J2
         DO 300 I=I1,I2
            IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IT = (I-I1)*(J2-J1+1)+(J-J1) + 1
            DIST1 = (XC(IL) - XX(IT))**2 + (YC(IL) - XX(IT+IOTOT))**2 + 
     +           (ZC(IL) - XX(IT+2*IOTOT))**2 
            DIST = DIST + DIST1
 300     CONTINUE
 200  CONTINUE

      D = SQRT(DIST/IOTOT)/REFLEN

      RETURN
      END SUBROUTINE IMPOXX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE APPRAX(XC,YC,ZC,XX,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,
     +     JLOC,JTRA,APP,JAPP,APU,RIAPU,MAXB,
     +    IFACEB,IB1,IB2,JB1,JB2,IORI,IOTOT,IOSIZE,N,L,NBG,IWA,IN,JN,KN)

C     K1 = K-value at the First Slab Outside the Block
C     K2 = K-value at the Second Slab Outside the Block
C     Note the inversion of the IJ-Base to account for the transposed
C     indices in the incoming data.
C     Indices on various walls are as follows
C     Wall No 1 -->    K,J,-I
C     Wall No 2 -->    I,K,-J
C     Wall No 3 -->    J,I,-K
C     Wall No 4 -->    J,K,I
C     Wall No 5 -->    K,I,J
C     Wall No 6 -->    I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      REAL :: APP(*), APU(IOSIZE,7), RESIDU, RIAPU(*) ! RIAPU is not used
      REAL :: XC(*),YC(*),ZC(*),XX(IOTOT,3)
      REAL :: A1(3,5), B1(3,5)
      REAL :: RX, RY, RZ, SX, SY, SZ, TX, TY, TZ
      INTEGER :: JLOC(*), JTRA(*), IAPU(1,1) !IAPU(MAXB/3,*) IAPU is not active
      INTEGER :: I1, I2, J1, J2, K1, K2, ISTR, JSTR, KSTR
      INTEGER :: I1B, I2B, J1B, J2B, ISTRB, JSTRB, KSTRB
      LOGICAL :: THERE

C ... APU(1,1)  Cell face areas on patch A
C ... APU(1,2)  Cell face areas on patch B
C ... APU(1,3)  Cell face areas on patch A for checking residual
C ... APU(1,4)  Cell face areas on patch B for checking residual
C ... APU(1,5)  for sortting (not used any more)
C ... IAPU(1,1)  for sortting
C ... IAPU(1,2)  for sortting
C ... IAPU(1,3)  for sortting

      IF(4*IOSIZE > MAXB) THEN
         WRITE(*,*) 'APPRAX: CONTACT PATRIK. EXITING ...',IOTOT,MAXB
         STOP
      ENDIF 


      IF(MOD(IFACEB,2) == 0) THEN
        I1B = 1
        I2B = IB2-IB1+2
        J1B = 1        
        J2B = JB2-JB1+2
      ELSE
        I1B = 1
        I2B = JB2-JB1+2
        J1B = 1        
        J2B = IB2-IB1+2
      ENDIF
      ISTRB = 1
      JSTRB = I2B
      KSTRB = 1

      RESIDU = 1.0E-7

C ... Check the directions of the face normals

C ... Face A normal

      I   = I1
      J   = J1
      K   = K1
      IL  = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
      IL1 = IL
      IL2 = IL+JSTR
      IL3 = IL+JSTR+ISTR
      IL4 = IL+ISTR
      P1X = XC(IL1)
      P1Y = YC(IL1)
      P1Z = ZC(IL1)
      P2X = XC(IL2)
      P2Y = YC(IL2)
      P2Z = ZC(IL2)
      P3X = XC(IL3)
      P3Y = YC(IL3)
      P3Z = ZC(IL3)
      P4X = XC(IL4)
      P4Y = YC(IL4)
      P4Z = ZC(IL4)

      AX  = P3X - P1X
      AY  = P3Y - P1Y
      AZ  = P3Z - P1Z

      BX  = P2X - P4X
      BY  = P2Y - P4Y
      BZ  = P2Z - P4Z

      EX  = AY*BZ - AZ*BY
      EY  = AZ*BX - AX*BZ
      EZ  = AX*BY - AY*BX

C ... Face B normal 

      IF(IORI == 0) THEN
        IR = 0
        JR = 0
      ELSEIF(IORI == 1) THEN
        IF(IFACEB == 1 .OR. IFACEB == 3 .OR. IFACEB == 5) THEN
          IR = 0
          JR = 1
        ELSE
          IR = 1
          JR = 0
        ENDIF
      ELSEIF(IORI == 2) THEN
        IR = 1
        JR = 1
      ELSEIF(IORI == 3) THEN
        IF(IFACEB == 1 .OR. IFACEB == 3 .OR. IFACEB == 5) THEN
          IR = 1
          JR = 0
        ELSE
          IR = 0
          JR = 1
        ENDIF
      ENDIF

      IC  = 1 + IR*(I2B-2)
      JC  = 1 + JR*(J2B-2)

      IT  = (JC-1)*JSTRB+IC
      IT1 = IT
      IT2 = IT+JSTRB
      IT3 = IT+JSTRB+ISTRB
      IT4 = IT+ISTRB

      P1X = XX(IT1,1)
      P1Y = XX(IT1,2)
      P1Z = XX(IT1,3)
      P2X = XX(IT2,1)
      P2Y = XX(IT2,2)
      P2Z = XX(IT2,3)
      P3X = XX(IT3,1)
      P3Y = XX(IT3,2)
      P3Z = XX(IT3,3)
      P4X = XX(IT4,1)
      P4Y = XX(IT4,2)
      P4Z = XX(IT4,3)

      AX  = P3X - P1X
      AY  = P3Y - P1Y
      AZ  = P3Z - P1Z

      BX  = P2X - P4X
      BY  = P2Y - P4Y
      BZ  = P2Z - P4Z

      FX  = AY*BZ - AZ*BY
      FY  = AZ*BX - AX*BZ
      FZ  = AX*BY - AY*BX

      DIR = SIGN(1.0,EX*FX + EY*FY + EZ*FZ)

C ... Compute the cell face areas on both patches.

C ... Cell face areas: patch A.

      K=K1
      ART = 0.0
      DO 151 J=J1,J2-1
         DO 152 I=I1,I2-1
            IK  = (J-J1)*(I2-I1)+I-I1+1
            IL0 = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IL1 = IL0
            IL2 = IL0+JSTR
            IL3 = IL0+JSTR+ISTR
            IL4 = IL0+ISTR
            A1(1,1) = XC(IL1)
            A1(2,1) = YC(IL1)
            A1(3,1) = ZC(IL1)
            A1(1,2) = XC(IL2)
            A1(2,2) = YC(IL2)
            A1(3,2) = ZC(IL2)
            A1(1,3) = XC(IL3)
            A1(2,3) = YC(IL3)
            A1(3,3) = ZC(IL3)
            A1(1,4) = XC(IL4)
            A1(2,4) = YC(IL4)
            A1(3,4) = ZC(IL4)

            RX = A1(1,3) - A1(1,1)
            RY = A1(2,3) - A1(2,1)
            RZ = A1(3,3) - A1(3,1)

            SX = A1(1,2) - A1(1,4)
            SY = A1(2,2) - A1(2,4)
            SZ = A1(3,2) - A1(3,4)

            TX = RY*SZ   - RZ*SY
            TY = RZ*SX   - RX*SZ
            TZ = RX*SY   - RY*SX

            APU(IK,1) = 0.5*SQRT(TX**2 + TY**2 + TZ**2)
            APU(IK,3) = APU(IK,1)
            ART       = ART + APU(IK,1)  ! Total area of patch A 
 152     CONTINUE
 151  CONTINUE

c      write(*,*) 'Patch A size = ',ART

C ... Cell face areas: patch B.

      BRT = 0.0
      DO 153 JB=J1B,J2B-1
         DO 154 IB=I1B,I2B-1
            JK  = (JB-1)*(JSTRB-1)+IB
            IT  = (JB-1)*JSTRB+IB
            IT1 = IT
            IT2 = IT+JSTRB
            IT3 = IT+JSTRB+ISTRB
            IT4 = IT+ISTRB
            B1(1,1) = XX(IT1,1)
            B1(2,1) = XX(IT1,2)
            B1(3,1) = XX(IT1,3)
            B1(1,2) = XX(IT2,1)
            B1(2,2) = XX(IT2,2)
            B1(3,2) = XX(IT2,3)
            B1(1,3) = XX(IT3,1)
            B1(2,3) = XX(IT3,2)
            B1(3,3) = XX(IT3,3)
            B1(1,4) = XX(IT4,1)
            B1(2,4) = XX(IT4,2)
            B1(3,4) = XX(IT4,3)

            RX = B1(1,3) - B1(1,1)
            RY = B1(2,3) - B1(2,1)
            RZ = B1(3,3) - B1(3,1)

            SX = B1(1,2) - B1(1,4)
            SY = B1(2,2) - B1(2,4)
            SZ = B1(3,2) - B1(3,4)

            TX = RY*SZ   - RZ*SY
            TY = RZ*SX   - RX*SZ
            TZ = RX*SY   - RY*SX

            APU(JK,2) = 0.5*SQRT(TX**2 + TY**2 + TZ**2)
            APU(JK,4) = APU(JK,2) 
            BRT       = BRT + APU(JK,2)  ! Total area of patch B 
 154     CONTINUE
 153  CONTINUE


C ... Compute the overlaping areas.

      IF(IR == 0) THEN
         IA =  I1B
         IE =  I2B-1
         IS =  1
      ELSE                   ! Reverse (local) I-loop
         IA =  I2B-1
         IE =  I1B
         IS = -1
      ENDIF
      
      IF(JR == 0) THEN
         JA =  J1B
         JE =  J2B-1
         JS =  1
      ELSE                   ! Reverse (local) J-loop
         JA =  J2B-1
         JE =  J1B
         JS = -1
      ENDIF

      JAPP = 0

      K=K1
      K2=K1
      IF(K == 1) K2 = 0
      DO 200 J=J1+1,J2-2
         DO 300 I=I1+1,I2-2
            IK  = (J-J1)*(I2-I1)+I-I1+1
            IL  = (K2-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IL0 = (K -1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IL1 = IL0
            IL2 = IL0+JSTR
            IL3 = IL0+JSTR+ISTR
            IL4 = IL0+ISTR
            A1(1,1) = XC(IL1)
            A1(2,1) = YC(IL1)
            A1(3,1) = ZC(IL1)
            A1(1,2) = XC(IL2)
            A1(2,2) = YC(IL2)
            A1(3,2) = ZC(IL2)
            A1(1,3) = XC(IL3)
            A1(2,3) = YC(IL3)
            A1(3,3) = ZC(IL3)
            A1(1,4) = XC(IL4)
            A1(2,4) = YC(IL4)
            A1(3,4) = ZC(IL4)

C ... first sweep is done for inner surfaces

            ISUM  = 0
            DO 120 JB=JA+JS,JE-JS,JS
               DO 130 IB=IA+IS,IE-IS,IS
                  JK  = (JB-1)*(JSTRB-1)+IB
                  IT  = (JB-1)*JSTRB+IB
                  IT1 = IT
                  IT2 = IT+JSTRB
                  IT3 = IT+JSTRB+ISTRB
                  IT4 = IT+ISTRB
                  IF(APU(JK,4) > RESIDU*APU(JK,2)) THEN
                     B1(1,1) = XX(IT1,1)
                     B1(2,1) = XX(IT1,2)
                     B1(3,1) = XX(IT1,3)
                     B1(1,2) = XX(IT2,1)
                     B1(2,2) = XX(IT2,2)
                     B1(3,2) = XX(IT2,3)
                     B1(1,3) = XX(IT3,1)
                     B1(2,3) = XX(IT3,2)
                     B1(3,3) = XX(IT3,3)
                     B1(1,4) = XX(IT4,1)
                     B1(2,4) = XX(IT4,2)
                     B1(3,4) = XX(IT4,3)

                     CALL OVERLA(A1,B1,APT,DIR)

                     IF(APT > RESIDU*APU(IK,1)) THEN
                        APU(IK,3) = APU(IK,3) - APT
                        APU(JK,4) = APU(JK,4) - APT
                        APT       = APT/APU(IK,1)
                        JAPP       = JAPP + 1
                        ISUM       = ISUM + 1
                        JLOC(JAPP) = IL   ! block local
                        JTRA(JAPP) = JK
                        APP(JAPP)  = APT

                     ENDIF
                        IF(APU(IK,3) < RESIDU*APU(IK,1)) GOTO 29
                     ENDIF
 130           CONTINUE
 120        CONTINUE
 29         CONTINUE

 300     CONTINUE
 200  CONTINUE

*C ... second is similar as first but it is done for whole patch
C ... second is done only for the ghost cells. (ESa 9.11.1998)

      DO 1200 J=J1,J2-1
         DO 1300 I=I1,I2-1
            IK  = (J-J1)*(I2-I1)+I-I1+1
            IL  = (K2-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IL0 = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IL1 = IL0
            IL2 = IL0+JSTR
            IL3 = IL0+JSTR+ISTR
            IL4 = IL0+ISTR
            IF(APU(IK,3) > RESIDU*APU(IK,1)) THEN
            A1(1,1) = XC(IL1)
            A1(2,1) = YC(IL1)
            A1(3,1) = ZC(IL1)
            A1(1,2) = XC(IL2)
            A1(2,2) = YC(IL2)
            A1(3,2) = ZC(IL2)
            A1(1,3) = XC(IL3)
            A1(2,3) = YC(IL3)
            A1(3,3) = ZC(IL3)
            A1(1,4) = XC(IL4)
            A1(2,4) = YC(IL4)
            A1(3,4) = ZC(IL4)

            ISUM  = 0
            DO 1120 JB=JA,JE,JS
               DO 1130 IB=IA,IE,IS
                  JK  = (JB-1)*(JSTRB-1)+IB
                  IT  = (JB-1)*JSTRB+IB
                  IT1 = IT
                  IT2 = IT+JSTRB
                  IT3 = IT+JSTRB+ISTRB
                  IT4 = IT+ISTRB
*                  IF(APU(JK,4) > RESIDU*APU(JK,2)) THEN
                     B1(1,1) = XX(IT1,1)
                     B1(2,1) = XX(IT1,2)
                     B1(3,1) = XX(IT1,3)
                     B1(1,2) = XX(IT2,1)
                     B1(2,2) = XX(IT2,2)
                     B1(3,2) = XX(IT2,3)
                     B1(1,3) = XX(IT3,1)
                     B1(2,3) = XX(IT3,2)
                     B1(3,3) = XX(IT3,3)
                     B1(1,4) = XX(IT4,1)
                     B1(2,4) = XX(IT4,2)
                     B1(3,4) = XX(IT4,3)
                        
                     CALL OVERLA(A1,B1,APT,DIR)

                     IF(APT > RESIDU*APU(IK,1)) THEN
                        APU(IK,3) = APU(IK,3) - APT
                        APU(JK,4) = APU(JK,4) - APT
                        APT       = APT/APU(IK,1)
                        JAPP       = JAPP + 1
                        ISUM       = ISUM + 1
                        JLOC(JAPP) = IL   ! block local
                        JTRA(JAPP) = JK
                        APP(JAPP)  = APT

                     ENDIF
                        IF(APU(IK,3) < RESIDU*APU(IK,1)) GOTO 129
*                  ENDIF
 1130           CONTINUE
 1120        CONTINUE
 129         CONTINUE

          ENDIF
 1300     CONTINUE
 1200  CONTINUE


C ... Compute the area residuals for checking purposes.

      ARR = 0.0
      DO J=J1,J2-1
        DO I=I1,I2-1
          IK  = (J-J1)*(I2-I1)+I-I1+1
          ARR = ARR + APU(IK,3)
c          if(APU(IK,3) < 0.0) write(*,*) APU(IK,3)
        ENDDO
      ENDDO

      BRR = 0.0
      DO JB=J1B,J2B-1
        DO IB=I1B,I2B-1
          JK  = (JB-1)*(I2B-1)+IB
          BRR = BRR + APU(JK,4)
c          if(APU(JK,4) < 0.0) write(*,*) APU(JK,4)
        ENDDO
      ENDDO

C ... put arrays in order of JLOC

 110  THERE = .FALSE. ! kupla sorttaus
      DO J = 1,JAPP-1
         IF(JLOC(J) > JLOC(J+1)) THEN
            THERE     = .TRUE.
            JLOC2     = JLOC(J)
            APP2      = APP(J) 
            JTRA2     = JTRA(J)
            JLOC(J)   = JLOC(J+1)
            APP(J)    = APP(J+1) 
            JTRA(J)   = JTRA(J+1)
            JLOC(J+1) = JLOC2
            APP(J+1)  = APP2
            JTRA(J+1) = JTRA2
         ENDIF
      ENDDO
      IF(THERE) GOTO 110

      
 111  THERE = .FALSE. ! kupla sorttaus
      DO J = 1,JAPP-1
         IF(JLOC(J) == JLOC(J+1) .AND. JTRA(J) > JTRA(J+1)) THEN
            THERE     = .TRUE.
            JLOC2     = JLOC(J)
            APP2      = APP(J) 
            JTRA2     = JTRA(J)
            JLOC(J)   = JLOC(J+1)
            APP(J)    = APP(J+1) 
            JTRA(J)   = JTRA(J+1)
            JLOC(J+1) = JLOC2
            APP(J+1)  = APP2
            JTRA(J+1) = JTRA2
         ENDIF
      ENDDO
      IF(THERE) GOTO 111

C ... remove double data
      IRES = 0
      DO J = 1,JAPP-1
         IF(JLOC(J) == JLOC(J+1) .AND. JTRA(J) == JTRA(J+1)) THEN
            IRES = IRES + 1
            APP(J)  = .5*(APP(J) + APP(J+1))
            DO I = J+1,JAPP-IRES
               JLOC(I) = JLOC(I+1)
               JTRA(I) = JTRA(I+1)
               APP(I)  = APP(I+1)
            ENDDO
         ENDIF
         IF(J+IRES >= JAPP-1) GOTO 145
      ENDDO
 145  JAPP = JAPP - IRES
               
C ... make normalization. The reasoning for this step is bit vague.
      II = 1

 113  ISTART = JLOC(II)
      XSUM   = 0.
      ISUM   = 0
      DO III = II,JAPP
         IF(JLOC(III) /= ISTART) GOTO 112
         ISUM = ISUM + 1
         XSUM = XSUM + APP(III)
      ENDDO
 112  DO III2 = II,MIN0(JAPP,ISUM+II-1)
         APP(III2) = APP(III2)/XSUM
      ENDDO
      II = II + ISUM
      IF(II <=  JAPP) GOTO 113

      INOHIT = 0
      DO J=J1,J2-1
         DO I=I1,I2-1
            IL  = (K2-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
            IHIT = 0
            DO JJ = 1,JAPP
               IF(JLOC(JJ) == IL) IHIT = IHIT + 1
            ENDDO
            IF(IHIT == 0) INOHIT = INOHIT + 1
         ENDDO
      ENDDO

C ... Check the area residuals.

      APR = 100.0 - ARR/ART*100.0
      BPR = 100.0 - BRR/BRT*100.0

      IF(IWA == 1) WRITE(*,312)
      WRITE(*,313) N,L,(J2-J1)*(I2-I1),JAPP,INOHIT,APR,BPR,DIR

 312  FORMAT(' NON-MATCHING BOUNDARY CONDITIONS ARE FOUND!!!!'/
     +     '    N  Patch    size  n.o.con nohit t.area    g.area',
     +                                     '   Direction')
 313  FORMAT(2(1X,I4),2(1X,I8),1X,I4,4F10.4)

      IF(INOHIT /= 0) THEN
         WRITE(*,314) 
         WRITE(13,314) 
         WRITE(*,312)
         WRITE(13,313) NBG,L,(J2-J1)*(I2-I1),JAPP,INOHIT,APR,BPR,DIR
 314     FORMAT(/' Warning APPRAX: all boundary areas have not ',
     +        'connection'/)

         WRITE(13,316) 
 316     FORMAT(' Corners of the nohit patch:'/
     +        '    I    J     X',13X,'Y',13X,'Z')
         CALL WRERRO(13,I1,J1,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)
         CALL WRERRO(13,I2,J1,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)
         CALL WRERRO(13,I2,J2,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)
         CALL WRERRO(13,I1,J2,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)
         WRITE(13,*)
         WRITE(13,*) ' Nohit locations:'
         DO J=J1,J2-1
            DO I=I1,I2-1
               IL  = (K2-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IHIT = 0
               DO JJ = 1,JAPP
                  IF(JLOC(JJ) == IL) IHIT = IHIT + 1
               ENDDO
               IF(IHIT == 0) THEN
                 CALL WRERRO(13,I,J,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE APPRAX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WRERRO(IOOUT,I,J,K2,ISTR,JSTR,KSTR,IN,JN,KN,XC,YC,ZC)

      REAL :: XC(*), YC(*), ZC(*)

      IL  = (K2-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
      
      WRITE(IOOUT,315) I,J,REAL(XC(IL),4),REAL(YC(IL),4),REAL(ZC(IL),4)
 315  FORMAT(2(1X,I4),3(1X,3E14.5))

      RETURN
      END SUBROUTINE WRERRO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPVOL(VOL,ZZZ,APU,ICON,NPATCH,IW,IMAX,JMAX,KMAX,ITAG,
     + IN,JN,KN,NPPV,IEXT,M,MPE,IPRO,N,L)

      USE MPI
      USE NS3CO, ONLY : IC9
      
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: VOL(*), ZZZ(*), APU(*)
      INTEGER :: ICON(IC9,*), ITAG(NPPV,*)
      LOGICAL :: CONN, CORN
      INTEGER :: IFACE
C
C ... EXPORT ALL THE BOUNDARIES OF THE BLOCK
C
      ISTR =  1
      JSTR =  IMAX+2*IN
      KSTR = (JMAX+2*JN)*JSTR

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = ((MPE ==  0 .AND. (ITYPE ==  1 .OR. ITYPE == 6  
     +   .OR. ITYPE == 11 .OR.   ITYPE == 14 .OR. ITYPE == 15)) 
     +   .OR.  (MPE ==  1 .AND. (ITYPE ==  1 .OR. ITYPE == 11 
     +   .OR. ITYPE == 14 .OR.   ITYPE == 15) .AND. ICON(20,L) == 0)
     +   .OR.  (MPE ==  2 .AND.  ITYPE == 6))
      
      CORN    = (.NOT. (ITYPE == 6 .AND. IEXT == 1))
      IF(CORN .AND. CONN) THEN

C ... CHECK THAT RECEIVING PATCH HAS AS HIGH MULTIGRID LEVEL
C ... (vois ehka laittaa impvolliinki. En pistanyt koska ei 
C ... testitapausta) PPR 24.11.95

      IF(M <= ICON(19,L)) THEN 
      IW1   = IW
      IFACE = ICON(3,L)   
      IVOL  = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)
      NBGL  = ICON(24,L)  ! global block number

C ... Extended connection

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + IEXT
      JM = JM + IEXT

      KA = 2
      KM = 1 + IEXT
      
      IOTOT = ABS((IM-IA+1)*(JM-JA+1)*(ABS(KA-KM)+1))
      IOCON = (IM2-IA2+1)*(JM2-JA2+1)*2

      IF(IVOL > 0 .AND. IOTY == 1) THEN
         ITAG(LP,IVOL) = 1
      ENDIF

      IF(IFACE == 1) THEN
         CALL EXPORT(VOL,APU,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN
         CALL EXPORT(VOL,APU,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN
         CALL EXPORT(VOL,APU,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN
         KA = IMAX - 1 + IEXT
         KM = IMAX
         CALL EXPORT(VOL,APU,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN
         KA = JMAX - 1 + IEXT
         KM = JMAX
         CALL EXPORT(VOL,APU,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN
         KA = KMAX - 1 + IEXT
         KM = KMAX 
         CALL EXPORT(VOL,APU,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such face. Exiting ...',ICON(3,L)
         STOP
      ENDIF

      IF(IOTY <= 1) CALL SETV12(ZZZ(IW1),APU,IOTOT)

      IF(IOTY >= 2) THEN
      IGP = ICON(2,L)   ! Process local patch (or block) number
      IGP = ICON(2,L)
        IF(IPRO <  JPRO) CALL MPI_RECV(ZZZ(IW1),IOCON,MPI_REAL8,
     &                   JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
        CALL MPI_SSEND(APU,IOTOT,MPI_REAL8,JPRO-1,JP,
     &       MPI_COMM_WORLD,IERR)
        IF(IPRO >= JPRO) CALL MPI_RECV(ZZZ(IW1),IOCON,MPI_REAL8,
     &                   JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
      ENDIF

      ENDIF                     ! Multigrid level check
      ENDIF                     ! Connections
1000  CONTINUE

      RETURN
      END SUBROUTINE EXPVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPVOL(VOL,ZZZ,APU,APU2,ICON,NPATCH,IW,IMAX,JMAX,KMAX,
     + ITAG,IN,JN,KN,NPPV,IEXT,IPRO,MPE,N,M,JLOC,JTRA,APP,IWAPP,JAPP)

      USE MPI
      USE NS3CO, ONLY : IC9
      
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: VOL(*), ZZZ(*), APU(*), APU2(*), APP(*)
      INTEGER :: JLOC(*), JTRA(*), IWAPP(*), JAPP(*)
      INTEGER :: IW(NPPV,*), ICON(IC9,*), ITAG(*)
      LOGICAL :: CONN, CORN
      INTEGER :: IFACE
C
C ... IMPORT DATA FOR ALL THE BOUNDARIES (= NPATCH) OF THE BLOCK
C
C ... IMPORT FROM BLOCK (= PROCESSOR) LV, (LOCAL) PATCH NO. LP
C
      ISTR =  1
      JSTR =  IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
      IF(ITAG(L) == 1) THEN

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = ((MPE ==  0 .AND. (ITYPE ==  1  .OR. ITYPE ==  6 
     +   .OR. ITYPE == 11 .OR.   ITYPE == 14  .OR. ITYPE == 15)) 
     +   .OR.  (MPE ==  1 .AND. (ITYPE ==  1  .OR. ITYPE == 11
     +   .OR. ITYPE == 14 .OR.   ITYPE == 15) .AND. ICON(20,L) == 0)
     +   .OR.  (MPE ==  2 .AND.  ITYPE ==  6))
      
C ... FOR CYCLIC PATCHES USE ONLY EXTRAPOLATION FOR CORNER POINTS
C ... OTHERWISE ALGORITHM IS COMPLEX AND IN GENERAL CASE ADVANTAGES
C ... ARE MINOR. PPR 28.9.1995
      CORN    = (.NOT. (ITYPE == 6 .AND. IEXT == 1))
      IF(CORN .AND. CONN) THEN
C ------------------------------------------------------------------

      IGP   = ICON(2,L)   ! Process local patch (or block) number
      IFACE = ICON(3,L)   ! patch face
      LV    = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)

      IW1   = IW(LP,LV)   ! IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2   = IW(L,N)     ! IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI
      NBLG  = ICON(24,L)  ! global block number

C ... Extended connection

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + IEXT
      JM = JM + IEXT
      
      IMA = IM - IA + 1
      JMA = JM - JA + 1
      
      KA = -1 + IEXT
      KM = 0
      IOTOT = IABS((IM2-IA2+1+IEXT)*(JM2-JA2+1+IEXT)*(KM-KA+1))
      
      IW1 = IW(LP,LV)     !IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2 = IW(L,N)       !IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI

      IF(IOTY <= 1) CALL SETV12(APU2,ZZZ(IW1),IOTOT)
      IF(IOTY >= 2) CALL SETV12(APU2,ZZZ(IW2),IOTOT)
C ... next line can be used if sendrecv command is used in expvol
C      IF(IOTY >= 2) CALL SETV12(APU2,ZZZ(IW2),IOTOT)
C ... RECEIVE FROM NODE ICON(8,L) DATA TO ZZZ-ARRAY WITH DISTR.MEM.

      IF(ICON(22,L) /= 0) THEN  ! NONMATCHING BC
         IF(IEXT /= 1 .AND. M == 1) THEN 
         KSTR2 = KSTR
         IF(IFACE == 1 .OR. IFACE == 4) KSTR2 = ISTR
         IF(IFACE == 2 .OR. IFACE == 5) KSTR2 = JSTR
         KSF = -KSTR2
         IF(IFACE >= 4) KSF = KSTR2
         KTRA = IOTOT/2
         DO I=1,JAPP(L)         ! but zeros
            I2 = IWAPP(L) + I - 1
            VOL(JLOC(I2)    ) = 0.
            VOL(JLOC(I2)+KSF) = 0.
         ENDDO
         DO I=1,JAPP(L)
            I2 = IWAPP(L) + I - 1
            KLOC = JLOC(I2)
            VOL(KLOC    ) = APP(I2)*APU2(JTRA(I2)+KTRA) + VOL(KLOC    )
            VOL(KLOC+KSF) = APP(I2)*APU2(JTRA(I2)     ) + VOL(KLOC+KSF)
         ENDDO
         ENDIF
      ELSE

C ------------------------------------------------------------------
      IN1 = IN
      IF(IEXT == 1) IN1 = 1
      IF(IFACE == 1) THEN                                ! WALL NO. 1
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN                            ! WALL NO. 2
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN                            ! WALL NO. 3
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN                            ! WALL NO. 4
         KA = IMAX + 2
         KM = IMAX + 1 + IEXT
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN                            ! WALL NO. 5
         KA = JMAX + 2
         KM = JMAX + 1 + IEXT
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN                            ! WALL NO. 6
         KA = KMAX + 2
         KM = KMAX + 1 + IEXT
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such wall. Exiting ...', ICON(3,L),L
         STOP
      ENDIF
      ENDIF ! ICON(22,L) /= 0

C ------------------------------------------------------------------
      ENDIF
      ENDIF

 1000 CONTINUE

      RETURN
      END SUBROUTINE IMPVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPVOD(VOL,ZZZ,APU,ICON,NPATCH,IW,IMAX,JMAX,KMAX,ITAG,
     + IN,JN,KN,NPPV,M,IPRO,N,L,IEXT)

      USE MPI
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: IEXT,NPATCH,IW,IMAX,JMAX,KMAX,IN,JN,KN,NPPV,M,IPRO,N,L,
     +     ISTR,JSTR,KSTR,ITYPE,IW1,IVOL,LP,JP,JPRO,IOTY,
     +     IA,IM,JA,JM,KA,KM,IA2,IM2,JA2,JM2,IOTOT,IOCON,IGP,IERR,IFACE
      INTEGER :: STATUS(MPI_STATUS_SIZE)

      REAL :: VOL(*), APU(*), ZZZ(*)

      INTEGER :: ICON(IC9,*), ITAG(NPPV,*)

      LOGICAL :: CONN

C ... EXPORT ALL THE BOUNDARIES OF THE BLOCK

      ISTR =  1
      JSTR =  IMAX+2*IN
      KSTR = (JMAX+2*JN)*JSTR

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = ((ITYPE ==  1 .OR. ITYPE == 11  .OR.
     &          ITYPE == 14 .OR. ITYPE == 15) .AND.ICON(20,L) == 0)
      
      IF(CONN)THEN

C ... CHECK THAT RECEIVING PATCH HAS AS HIGH MULTIGRID LEVEL
C ... (vois ehka laittaa impvolliinki. En pistanyt koska ei 
C ... testitapausta) PPR 24.11.95

      IF(M <= ICON(19,L)) THEN 

      IW1   = IW
      IFACE = ICON(3,l)   
      IVOL  = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)

C ... Extended connection 

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + IEXT
      JM = JM + IEXT

      KA = 2
      KM = 1 + IEXT
      
      IOTOT = ABS((IM-IA+1)*(JM-JA+1)*(ABS(KA-KM)+1))
      IOCON = (IM2-IA2+1+IEXT)*(JM2-JA2+1+IEXT)*2

      IF(IVOL > 0 .AND. IOTY == 1) THEN
         ITAG(LP,IVOL) = 1
      ENDIF

      IF(IFACE == 1) THEN
         CALL EXPORD(VOL,APU,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN
         CALL EXPORD(VOL,APU,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN
         CALL EXPORD(VOL,APU,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN
         KA = IMAX - 1 + IEXT
         KM = IMAX
         CALL EXPORD(VOL,APU,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN
         KA = JMAX - 1 + IEXT
         KM = JMAX
         CALL EXPORD(VOL,APU,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN
         KA = KMAX - 1 + IEXT
         KM = KMAX 
         CALL EXPORD(VOL,APU,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such face. Exiting ...',ICON(3,L)
         STOP
      ENDIF

      IF(IOTY <= 1) CALL SETVDD(ZZZ(IW1),APU,2*IOTOT)

      IF(IOTY >= 2) THEN
      IGP = ICON(2,L)  ! Process local patch (or block) number
      IGP = ICON(2,L)
        IF(IPRO < JPRO) CALL MPI_RECV(ZZZ(IW1),IOCON,
     +     MPI_REAL8,JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
        CALL MPI_SSEND(APU,IOTOT,MPI_REAL8,
     +     JPRO-1,JP,MPI_COMM_WORLD,IERR)
        IF(IPRO >= JPRO) CALL MPI_RECV(ZZZ(IW1),IOCON,
     +     MPI_REAL8,JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
      ENDIF                     ! IOTY >= 2

      ENDIF                     ! Multigrid level check
      ENDIF                     ! CONN

      RETURN
      END SUBROUTINE EXPVOD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPVOD(VOL,ZZZ,APU,APU2,ICON,NPATCH,IW,IMAX,JMAX,KMAX,
     + ITAG,IN,JN,KN,NPPV,IPRO,N,M,IEXT)

      USE MPI
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: IEXT,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,NPPV,IPRO,N,M,ISTR,
     +   JSTR,KSTR,L,ITYPE,IGP,IFACE,LV,LP,JP,JPRO,IOTY,IW1,IW2,
     +   IA,IM,JA,JM,KA,KM,IA2,IM2,JA2,JM2,IMA,JMA,IOTOT,IN1
      INTEGER :: STATUS(MPI_STATUS_SIZE)

      REAL :: VOL(*), APU(*), APU2(*), ZZZ(*)

      INTEGER :: IW(NPPV,*), ICON(IC9,*), ITAG(*)

      LOGICAL :: CONN, CORN
C
C ... IMPORT DATA FOR ALL THE BOUNDARIES (= NPATCH) OF THE BLOCK
C
C ... IMPORT FROM BLOCK (= PROCESSOR) LV, (LOCAL) PATCH NO. LP
C
      ISTR=1
      JSTR=IMAX +2*IN
      KSTR=(JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
      IF(ITAG(L) == 1) THEN

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = ((ITYPE ==  1  .OR. ITYPE == 11 .OR. ITYPE == 14
     &     .OR. ITYPE == 15) .AND. ICON(20,L) == 0)
      
      IF(CONN) THEN
C ------------------------------------------------------------------

      IGP   = ICON(2,L)   ! Process local patch (or block) number
      IFACE = ICON(3,L)   ! patch face
      LV    = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)

      IW1   = IW(LP,LV)   ! IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2   = IW(L,N)     ! IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI

C ... Extended connection 

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + IEXT
      JM = JM + IEXT
      
      IMA = IM - IA + 1
      JMA = JM - JA + 1
      
      KA = -1 + IEXT
      KM = 0
      IOTOT = IABS((IM2-IA2+1+IEXT)*(JM2-JA2+1+IEXT)*(KM-KA+1))
      
      IW1 = IW(LP,LV)     !IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2 = IW(L,N)       !IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI

      IF(IOTY <= 1) CALL SETVDD(APU2,ZZZ(IW1),IOTOT)
      IF(IOTY >= 2) CALL SETVDD(APU2,ZZZ(IW2),IOTOT)

      IF(ICON(22,L) == 0) THEN  ! MATCHING BC

C ------------------------------------------------------------------
      IN1 = IN
      IF(IEXT == 1) IN1 = 1
      IF(IFACE == 1) THEN                          ! WALL NO. 1
         CALL PERMUD(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN                      ! WALL NO. 2
         CALL PERMUD(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN                      ! WALL NO. 3
         CALL PERMUD(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN                      ! WALL NO. 4
         KA = IMAX + 2
         KM = IMAX + 1 + IEXT
         CALL PERMUD(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN                      ! WALL NO. 5
         KA = JMAX + 2
         KM = JMAX + 1 + IEXT
         CALL PERMUD(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN                      ! WALL NO. 6
         KA = KMAX + 2
         KM = KMAX + 1 + IEXT
         CALL PERMUD(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORD(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such wall. Exiting ...', ICON(3,L),L
         STOP
      ENDIF
      ENDIF ! ICON(22,L) /= 0

C ------------------------------------------------------------------
C      ENDIF                                   ! CYCLIC CORNERS
      ENDIF                                   ! END OF CONNECTIONS
      ENDIF

 1000 CONTINUE

      RETURN
      END SUBROUTINE IMPVOD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPVXX(XC,YC,ZC,ZZZ,APU,ICON,NPATCH,IW,IMAX,JMAX,KMAX,
     + ITAG,IN,JN,KN,NPPV,M,IPRO,N,L)

      USE MPI
      USE NS3CO, ONLY : IC9

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: XC(*), YC(*), ZC(*), APU(*), ZZZ(*)
      INTEGER :: ICON(IC9,*), ITAG(NPPV,*)
      LOGICAL :: CONN
      INTEGER :: IFACE

C ... EXPORT ALL THE BOUNDARIES OF THE BLOCK
   
      ISTR =  1
      JSTR =  IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = (ITYPE == 1)
      
      IF(CONN) THEN

C ... CHECK THAT RECEIVING PATCH HAS AS HIGH MULTIGRID LEVEL
C ... (vois ehka laittaa impvolliinki. En pistanyt koska ei 
C ... testitapausta) PPR 24.11.95

      IF(M <= ICON(19,L)) THEN 

      IW1   = IW
      IFACE = ICON(3,L)   
      IVOL  = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)

C ... Extended connection 

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + 1
      JM = JM + 1

      KA = 1
      KM = 1

      IOTOT = ABS((IM-IA+1)*(JM-JA+1)*(ABS(KA-KM)+1))
      IOTO2 = 2*IOTOT
      IOCON = (IM2-IA2+2)*(JM2-JA2+2)*1

      IF(IVOL > 0 .AND. IOTY == 1) THEN
         ITAG(LP,IVOL) = 1
      ENDIF

      IF(IFACE == 1) THEN
         CALL EXPORD(XC,APU(1      ),JA,JM,IA,IM,KA,KM,
     +                               KSTR,JSTR,ISTR,KN,JN,IN)
         CALL EXPORD(YC,APU(1+IOTOT),JA,JM,IA,IM,KA,KM,
     +                               KSTR,JSTR,ISTR,KN,JN,IN)
         CALL EXPORD(ZC,APU(1+IOTO2),JA,JM,IA,IM,KA,KM,
     +                               KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN
         CALL EXPORD(XC,APU(1      ),IA,IM,JA,JM,KA,KM,
     +                               ISTR,KSTR,JSTR,IN,KN,JN)
         CALL EXPORD(YC,APU(1+IOTOT),IA,IM,JA,JM,KA,KM,
     +                               ISTR,KSTR,JSTR,IN,KN,JN)
         CALL EXPORD(ZC,APU(1+IOTO2),IA,IM,JA,JM,KA,KM,
     +                               ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN
         CALL EXPORD(XC,APU(1      ),JA,JM,IA,IM,KA,KM,
     +                               JSTR,ISTR,KSTR,JN,IN,KN)
         CALL EXPORD(YC,APU(1+IOTOT),JA,JM,IA,IM,KA,KM,
     +                               JSTR,ISTR,KSTR,JN,IN,KN)
         CALL EXPORD(ZC,APU(1+IOTO2),JA,JM,IA,IM,KA,KM,
     +                               JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN
         KA = IMAX + 1
         KM = IMAX + 1
         CALL EXPORD(XC,APU(1      ),IA,IM,JA,JM,KA,KM,
     +                               JSTR,KSTR,ISTR,JN,KN,IN)
         CALL EXPORD(YC,APU(1+IOTOT),IA,IM,JA,JM,KA,KM,
     +                               JSTR,KSTR,ISTR,JN,KN,IN)
         CALL EXPORD(ZC,APU(1+IOTO2),IA,IM,JA,JM,KA,KM,
     +                               JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN
         KA = JMAX + 1
         KM = JMAX + 1
         CALL EXPORD(XC,APU(1      ),JA,JM,IA,IM,KA,KM,
     +                               KSTR,ISTR,JSTR,KN,IN,JN)
         CALL EXPORD(YC,APU(1+IOTOT),JA,JM,IA,IM,KA,KM,
     +                               KSTR,ISTR,JSTR,KN,IN,JN)
         CALL EXPORD(ZC,APU(1+IOTO2),JA,JM,IA,IM,KA,KM,
     +                               KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN
         KA = KMAX + 1
         KM = KMAX + 1
         CALL EXPORD(XC,APU(1      ),IA,IM,JA,JM,KA,KM,
     +                               ISTR,JSTR,KSTR,IN,JN,KN)
         CALL EXPORD(YC,APU(1+IOTOT),IA,IM,JA,JM,KA,KM,
     +                               ISTR,JSTR,KSTR,IN,JN,KN)
         CALL EXPORD(ZC,APU(1+IOTO2),IA,IM,JA,JM,KA,KM,
     +                               ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such face. Exiting ...',ICON(3,L)
         STOP
      ENDIF

      IGP = ICON(2,L)  ! Process local patch (or block) number

      IF(IOTY <= 1) CALL SETVDD(ZZZ(IW1),APU,6*IOTOT)

      IF(IOTY >= 2) THEN
      IGP = ICON(2,L)
        IF(IPRO < JPRO) CALL MPI_RECV(ZZZ(IW1),6*IOCON,
     +     MPI_REAL8,JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
        CALL MPI_SSEND(APU,6*IOTOT,MPI_REAL8,
     +     JPRO-1,JP,MPI_COMM_WORLD,IERR)
        IF(IPRO >= JPRO) CALL MPI_RECV(ZZZ(IW1),6*IOCON,
     +     MPI_REAL8,JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
      ENDIF

      ENDIF ! MULTIGRID LEVEL CHECK
      ENDIF ! 
1000  CONTINUE

      RETURN
      END SUBROUTINE EXPVXX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPVXX(XC,YC,ZC,ZZZ,APU,APU2,ICON,NPATCH,IW,IMAX,JMAX,
     +     KMAX,ITAG,IN,JN,KN,NPPV,IPRO,M,N,NOMATC)

      USE NS3CO, ONLY : IC9
      
      REAL    :: XC(*), YC(*), ZC(*), APU(*), APU2(*), ZZZ(*)
      INTEGER :: IW(NPPV,*), ICON(IC9,*), ITAG(*)
      LOGICAL :: CONN, NOMATC
      INTEGER :: IFACE

C ... IMPORT DATA FOR ALL THE BOUNDARIES (= NPATCH) OF THE BLOCK

C ... IMPORT FROM BLOCK (= PROCESSOR) LV, (LOCAL) PATCH NO. LP

      ISTR =  1
      JSTR =  IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH

      IF(ITAG(L) == 1) THEN

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = (ITYPE == 1)
      
      IF(CONN) THEN
C ------------------------------------------------------------------

      IFACE = ICON(3,L)
      LV    = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)
      IORI  = ICON(14,L)  ! orientation

      IW1   = IW(LP,LV)   ! IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2   = IW(L,N)     ! IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI

C ... Extended connection 

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IM = IM + 1
      JM = JM + 1
      
      IMA = IM - IA + 1
      JMA = JM - JA + 1 
      
      KA = 1
      KM = 1

      IOTOT = IABS((IM2-IA2+2)*(JM2-JA2+2))
      IOTO2 = 2*IOTOT

      IF(IOTY <= 1) CALL SETVDD(APU2,ZZZ(IW1),6*IOTOT)
      IF(IOTY >= 2) CALL SETVDD(APU2,ZZZ(IW2),6*IOTOT)


C ... RECEIVE FROM NODE ICON(8,L) DATA TO ZZZ-ARRAY WITH DISTR.MEM.
C ------------------------------------------------------------------
      IF(ICON(22,L) == 0) THEN
      IN1 = 1
      IF(IFACE == 1) THEN                          ! WALL NO. 1
         CALL PERMUD(ICON(14,L),APU2         ,APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,IMA,JMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,
     +               KN,JN,IN,DD)

      ELSEIF(IFACE == 2) THEN                      ! WALL NO. 2
         CALL PERMUD(ICON(14,L),APU2         ,APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,JMA,IMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,
     +               IN,KN,JN,DD)

      ELSEIF(IFACE == 3) THEN                      ! WALL NO. 3
         CALL PERMUD(ICON(14,L),APU2         ,APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,IMA,JMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,
     +               JN,IN,KN,DD)

      ELSEIF(IFACE == 4) THEN                      ! WALL NO. 4
         KA = IMAX + 1
         KM = IMAX + 1 
         CALL PERMUD(ICON(14,L),APU2         ,APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,JMA,IMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,
     +               JN,KN,IN,DD)

      ELSEIF(IFACE == 5) THEN                      ! WALL NO. 5
         KA = JMAX + 1
         KM = JMAX + 1
         CALL PERMUD(ICON(14,L),APU2         ,APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,IMA,JMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,IMA,JMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,
     +               KN,IN,JN,DD)

      ELSEIF(IFACE == 6) THEN                      ! WALL NO. 6
         KA = KMAX + 1
         KM = KMAX + 1
         CALL PERMUD(ICON(14,L),APU2         ,APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTOT),APU,JMA,IMA,0,0,IN1)
         CALL PERMUD(ICON(14,L),APU2(1+IOTO2),APU,JMA,IMA,0,0,IN1)
         CALL IMPOXX(XC,YC,ZC,APU2,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,
     +               IN,JN,KN,DD)
      ELSE
         WRITE(*,*) 'No such wall. Exiting ...', ICON(3,L),L
         STOP
      ENDIF

C     ICON(22,L) = 1  ! warning all connections are made by a new method

      IF(ICON(20,L) /= 0) DD = 0.
      IF(DD > 1.5E-4) ICON(22,L) = 1
      IF(DD > 1.5E-4) write(*,*) 'block',ICON(24,L),' patch',l,
     +     ' dist',dd
      IF(ICON(22,L) == 1) NOMATC = .TRUE.
      ENDIF                     ! (ICON(22,L) == 0)

      ENDIF
      ENDIF                     ! END OF CONNECTIONS

 1000 CONTINUE

      RETURN
      END SUBROUTINE IMPVXX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPTOT(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,
     +    ITURB,ISTRES,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,MAXSS,NSCAL,
     +    ZZZ,APU,ICON,NPATCH,IW,IMAX,JMAX,KMAX,ITAG,IN,JN,KN,NPPV,M,
     +    IPRO,N,L,MULPHL,TRANSL,PEUDO1,MULPHC,TWO_FLUIDL)

      USE MPI
      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: ITURB,ISTRES,MAXEB,MAXSB,MAXSS,NSCAL,NPATCH,IW,
     +           IMAX,JMAX,KMAX,IN,JN,KN,NPPV,M,IPRO,N,L,IFACE
      INTEGER :: NTOT,ILENG,ILNG,ISTR,JSTR,KSTR,ITYPE,IW1,IVOL,LP,JP,
     +           JPRO,IOTY,NBGL,IA,IM,JA,JM,KA,KM,IA2,IM2,JA2,JM2,
     +           IOTOT,IOCON,IGP,IERR
      LOGICAL :: MULPHL, TRANSL, TWO_FLUIDL
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: RO(*),RM(*),RN(*),RW(*),E(*),EPS2(*),VIST(*),P(*),
     +           PDIFF(*),RK(*),REPS(*),FI(MAXSB,*),S11(MAXSS,*),
     +           ZZZ(*),APU(*),CH(*),BIJ(MAXEB,*),PEUDO1(*)
      INTEGER :: ICON(IC9,*), ITAG(NPPV,*)
      CHARACTER(*) :: MULPHC

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

C
C ... EXPORT ALL THE BOUNDARIES OF THE BLOCK (IOTY = 1, 6, 11, 14 or 15)
C
C ... CHECK OUT LENGHT OF THE IO PATHC

      CALL IWLENG(ILENG,ITURB,NSCAL,ISTRES,MULPHL,NPHASES,TRANSL,
     +   TWO_FLUIDL)

      ILNG = ILENG

      ISTR =  1
      JSTR =  IMAX+2*IN
      KSTR = (JMAX+2*JN)*JSTR
      NTOT = (IMAX+2*IN)*(JMAX+2*JN)*(KMAX+2*KN)

      ITYPE = ICON(1,L)  ! Boundary condition type
C ... PERIODIC VELOCITIES ARE NOT CONNECTED
      IF(ITYPE == 6) ILNG = ILENG - 3

C ... CHECK THAT RECEIVING PATCH HAS AS HIGH MULTIGRID LEVEL
C ... (vois ehka laittaa impvolliinki. En pistanyt koska ei 
C ... testitapausta) PPR 24.11.95

      IF(M <= ICON(19,L)) THEN 

      IW1   = IW
         
      IFACE = ICON(3,L)
      IVOL  = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)
      NBGL  = ICON(24,L)  ! global block number

C ... Extended connection

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)      
      
      KA = 2
      KM = 1
      
      IOTOT = ABS((IM-IA+1)*(JM-JA+1)*(KA-KM+1))
      IOCON = (IM2-IA2+1)*(JM2-JA2+1)*2

      IF(IVOL > 0 .AND. IOTY == 1) THEN
         ITAG(LP,IVOL) = 1
      ENDIF

      IF(IFACE == 1) THEN
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,JA,JM,IA,IM,KA,KM,
     +        KSTR,JSTR,ISTR,KN,JN,IN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)

      ELSEIF(IFACE == 2) THEN
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,IA,IM,JA,JM,KA,KM,
     +        ISTR,KSTR,JSTR,IN,KN,JN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)

      ELSEIF(IFACE == 3) THEN
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,JA,JM,IA,IM,KA,KM,
     +        JSTR,ISTR,KSTR,JN,IN,KN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)

      ELSEIF(IFACE == 4) THEN
         KA = IMAX - 1
         KM = IMAX
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,IA,IM,JA,JM,KA,KM,
     +        JSTR,KSTR,ISTR,JN,KN,IN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)

      ELSEIF(IFACE == 5) THEN
         KA = JMAX - 1
         KM = JMAX
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,JA,JM,IA,IM,KA,KM,
     +        KSTR,ISTR,JSTR,KN,IN,JN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)


      ELSEIF(IFACE == 6) THEN
         KA = KMAX - 1
         KM = KMAX 
         CALL EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,ITURB,
     +        ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +        MAXSS,NSCAL,IOTOT,ITYPE,APU,IA,IM,JA,JM,KA,KM,
     +        ISTR,JSTR,KSTR,IN,JN,KN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)
      ELSE
         WRITE(*,*) 'No such face. Exiting ...',ICON(3,L)
         STOP
      ENDIF

      IF(IOTY <= 1) CALL SETV12(ZZZ(IW1),APU,ILNG*IOTOT)

      IF(IOTY >= 2) THEN
      IGP = ICON(2,L)
        IF(IPRO < JPRO) CALL MPI_RECV(ZZZ(IW1),ILNG*IOCON,MPI_REAL8,
     +     JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
        CALL MPI_SSEND(APU,ILNG*IOTOT,MPI_REAL8,
     +     JPRO-1,JP,MPI_COMM_WORLD,IERR)
        IF(IPRO >= JPRO) CALL MPI_RECV(ZZZ(IW1),ILNG*IOCON,MPI_REAL8,
     +     JPRO-1,IGP,MPI_COMM_WORLD,STATUS,IERR)
      ENDIF

      ENDIF

      RETURN
      END SUBROUTINE EXPTOT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IMPVO2(VOL,ZZZ,APU,APU2,ICON,NPATCH,IW,IMAX,JMAX,KMAX,
     + ITAG,IN,JN,KN,NPPV,IPRO,N,M,JLOC,JTRA,APP,IWAPP,JAPP,NCAS2)

      USE MPI
      USE NS3CO, ONLY : IC9

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL    :: VOL(*), ZZZ(*), APU(*), APU2(*), APP(*)
      INTEGER :: JLOC(*), JTRA(*), IWAPP(*), JAPP(*)
      INTEGER :: IW(NPPV,*), ICON(IC9,*), ITAG(*)
      LOGICAL :: CONN, CORN
      INTEGER :: IFACE
C
C ... IMPORT DATA FOR ALL THE BOUNDARIES (= NPATCH) OF THE BLOCK
C
C ... IMPORT FROM BLOCK (= PROCESSOR) LV, (LOCAL) PATCH NO. LP
C
      ISTR =  1
      JSTR =  IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH

      NCASE = NCAS2
      IF(ITAG(L) == 1) THEN

      ITYPE = ICON(1,L)  ! Boundary condition type
      CONN  = (ITYPE == 1  .OR. ITYPE == 11 .OR. ITYPE == 14
     +      .OR. ITYPE == 15)  .OR.  (ITYPE == 6 .AND. 
     +     (NCAS2 /= 1 .AND. NCAS2 /= 2 .AND. NCAS2 /= 3))
      IF(CONN) THEN
C ------------------------------------------------------------------
      IF(ITYPE == 6 .AND. NCASE /= 0) NCASE = NCAS2 - 3 

      IGP   = ICON(2,L)   ! Process local patch (or block) number
      IFACE = ICON(3,L)   ! patch face
      LV    = ICON(8,L)   ! number of connective block (process local)
      LP    = ICON(15,L)  ! connective patch n. (block l.)
      JP    = ICON(16,L)  ! connective patch n. (proces l.)
      JPRO  = ICON(17,L)  ! number of connectivity process
      IOTY  = ICON(18,L)  ! connectivity method (1=normal,3=MPI)
      NBLG  = ICON(24,L)  ! global block number

C ... Entended connction

      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)
      
      IMA = IM - IA + 1
      JMA = JM - JA + 1
      
      KA = -1
      KM = 0
      IOTOT = IABS((IM2-IA2+1)*(JM2-JA2+1)*(KM-KA+1))
      
      IAL = NCASE*IOTOT
      IW1 = IW(LP,LV) + IAL ! IW1 IS THE ADDRESS OF THE DESIRED PATCH
      IW2 = IW(L,N)   + IAL ! IW2 IS THE ADDRESS OF THE DESIRED PATCH IN MPI
      IF(IOTY <= 1) CALL SETV12(APU2,ZZZ(IW1),IOTOT)
      IF(IOTY >= 2) CALL SETV12(APU2,ZZZ(IW2),IOTOT)
C ... next line can be used if sendrecv command is used in expvol
C      IF(IOTY >= 2) CALL SETV12(APU2,ZZZ(IW2),IOTOT)
C ... RECEIVE FROM NODE ICON(8,L) DATA TO ZZZ-ARRAY WITH DISTR.MEM.

      IF(ICON(22,L) /= 0) THEN  ! NONMATCHING BC
         IF(M == 1) THEN 
         KSTR2 = KSTR
         IF(IFACE == 1 .OR. IFACE == 4) KSTR2 = ISTR
         IF(IFACE == 2 .OR. IFACE == 5) KSTR2 = JSTR
         KSF = -KSTR2
         IF(IFACE >= 4) KSF = KSTR2
         KTRA = IOTOT/2
         DO I=1,JAPP(L)         ! but zeros
            I2 = IWAPP(L) + I - 1
            VOL(JLOC(I2)    ) = 0.
            VOL(JLOC(I2)+KSF) = 0.
         ENDDO
         DO I=1,JAPP(L)
            I2 = IWAPP(L) + I - 1
            KLOC = JLOC(I2)
            VOL(KLOC    ) = APP(I2)*APU2(JTRA(I2)+KTRA) + VOL(KLOC    )
            VOL(KLOC+KSF) = APP(I2)*APU2(JTRA(I2)     ) + VOL(KLOC+KSF)
         ENDDO
         ENDIF
      ELSE

C ------------------------------------------------------------------
      IN1 = IN
      IF(IFACE == 1) THEN                          ! WALL NO. 1
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,JSTR,ISTR,KN,JN,IN)

      ELSEIF(IFACE == 2) THEN                      ! WALL NO. 2
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,KSTR,JSTR,IN,KN,JN)

      ELSEIF(IFACE == 3) THEN                      ! WALL NO. 3
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,JSTR,ISTR,KSTR,JN,IN,KN)

      ELSEIF(IFACE == 4) THEN                      ! WALL NO. 4
         KA = IMAX + 2
         KM = IMAX + 1
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,JSTR,KSTR,ISTR,JN,KN,IN)

      ELSEIF(IFACE == 5) THEN                      ! WALL NO. 5
         KA = JMAX + 2
         KM = JMAX + 1
         CALL PERMUT(ICON(14,L),APU2,APU,IMA,JMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,JA,JM,IA,IM,KA,KM,KSTR,ISTR,JSTR,KN,IN,JN)

      ELSEIF(IFACE == 6) THEN                      ! WALL NO. 6
         KA = KMAX + 2
         KM = KMAX + 1
         CALL PERMUT(ICON(14,L),APU2,APU,JMA,IMA,0,0,IN1)
         CALL IMPORT(VOL,APU2,IA,IM,JA,JM,KA,KM,ISTR,JSTR,KSTR,IN,JN,KN)
      ELSE
         WRITE(*,*) 'No such wall. Exiting ...', ICON(3,L),L
         STOP
      ENDIF

      ENDIF                                   ! end nonmatching
C ------------------------------------------------------------------
C      ENDIF                                   ! CYCLIC CORNERS
      ENDIF                                   ! END OF CONNECTIONS
      ENDIF

 1000 CONTINUE

      RETURN
      END SUBROUTINE IMPVO2
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE EXPTOR(RO,RM,RN,RW,E,EPS2,VIST,CH,P,PDIFF,RK,REPS,
     +   ITURB,ISTRES,MULPHL,TRANSL,FI,S11,BIJ,PRO,VAR,TRM,MAXEB,MAXSB,
     +   MAXSS,NSCAL,IOTOT,ITYPE,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,
     +   IN,JN,KN,NTOT,PEUDO1,MULPHC,TWO_FLUIDL)

C     Indices on various walls are as follows
C     Wall No 1 --> K,J,-I
C     Wall No 2 --> I,K,-J
C     Wall No 3 --> J,I,-K
C     Wall No 4 --> K,I,J
C     Wall No 5 --> J,K,I
C     Wall No 6 --> I,J,K
C
C     IT = Pointer for the Temporary Array
C     IL = Pointer for the Local Array

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: ITURB,ISTRES,MAXEB,MAXSB,MAXSS,NSCAL,IOTOT,ITYPE,I1,I2,
     +           J1,J2,K1,K2,ISTR,JSTR,KSTR,NTOT,IPHASE
      INTEGER :: KSTP, IOT, NEX, NS, IN, JN, KN
      LOGICAL :: MULPHL, TRANSL, TWO_FLUIDL
      REAL    :: RO(*),RM(*),RN(*),RW(*),E(*),EPS2(*),VIST(*),P(*),
     +           PDIFF(*),RK(*),REPS(*),FI(MAXSB,*),S11(MAXSS,*),Q(*),
     +           BIJ(MAXEB,*),CH(*)
      REAL    :: PEUDO1(*)

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)

      CHARACTER(*)           :: MULPHC

        
      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF

      IOT = IOTOT

C ... IOTY = 1  <=> Connectivity
C ... IOTY = 6  <=> Cyclic
C ... IOTY = 11 <=> Sliding surface
C ... IOTY = 14 <=> Mixing surface
C ... IOTY = 15 <=> Connection with a solid block

      IF(ITYPE == 1 .OR. ITYPE == 11 .OR. ITYPE == 14) THEN
      NEX = 10
      CALL RTOQ(RO,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,0)
      CALL RTOQ(RM,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,1)
      CALL RTOQ(RN,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,2)
      CALL RTOQ(RW,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,3)
      CALL RTOQ(E,    Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,4)
      CALL RTOQ(EPS2 ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,5)
      CALL RTOQ(VIST ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,6)
      CALL RTOQ(CH   ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,7)
      CALL RTOQ(P    ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,8)
      CALL RTOQ(PDIFF,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,9)
      ELSE IF(ITYPE == 6) THEN
      NEX = 7
      CALL RTOQ(RO,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,0)
      CALL RTOQ(E,    Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,1)
      CALL RTOQ(EPS2 ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,2)
      CALL RTOQ(VIST ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,3)
      CALL RTOQ(CH   ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,4)
      CALL RTOQ(P    ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,5)
      CALL RTOQ(PDIFF,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,6)
      ELSE IF(ITYPE == 15) THEN
      NEX = 10
      CALL RTOQ(RO,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,0)
      CALL RTOQ(RM,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,1)
      CALL RTOQ(RN,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,2)
      CALL RTOQ(RW,   Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,3)
      CALL RTOQ(E,    Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,4)
      CALL RTOQ(EPS2 ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,5)
      CALL RTOQ(VIST ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,6)
      CALL RTOQ(CH   ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,7)
      CALL RTOQ(P    ,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,8)
      CALL RTOQ(PDIFF,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,9)
      ENDIF

C ... K-E MODEL
      IF(ITURB >= 3 .AND. ITURB /= 8) THEN
      CALL RTOQ(RK,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX)
      CALL RTOQ(REPS,Q,I1,I2,J1,J2,K1,K2,
     &          ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+1)
      NEX = NEX + 2
      ENDIF

C ... SCALARS
      IF(NSCAL >= 1) THEN
      DO 720 NS = 1,NSCAL
         CALL RTOQ(FI(1,NS),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX)
      NEX = NEX + 1
 720  CONTINUE
      ENDIF

C ... ARSM
      IF(ITURB >= 10 .AND. ITURB <= 19 .OR. ISTRES > 0) THEN
        CALL RTOQ(S11(1,1),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX)
        CALL RTOQ(S11(1,2),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+1)
        CALL RTOQ(S11(1,3),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+2)
        CALL RTOQ(S11(1,4),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+3)
        CALL RTOQ(S11(1,5),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+4)
        CALL RTOQ(S11(1,6),Q,I1,I2,J1,J2,K1,K2,
     &            ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+5)
         NEX = NEX + 6
      ENDIF

C ... EARSM
      IF(ISTRES > 0) THEN
         CALL RTOQ(BIJ(1,1),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX)
         CALL RTOQ(BIJ(1,2),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+1)
         CALL RTOQ(BIJ(1,3),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+2)
         CALL RTOQ(BIJ(1,4),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+3)
         CALL RTOQ(BIJ(1,5),Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX+4)
         NEX = NEX + 5
      ENDIF
         
C ... Multiphase variables (Current ILENG = 3-6*NPHASES)
      IF(MULPHL) THEN

      DO IPHASE = 1,NPHASES
         CALL TTOQ(PRO,VAR,Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,IPHASE,NEX)
         NEX = NEX + 3
         IF(TWO_FLUIDL) THEN
           CALL TTUQ(PRO,VAR,Q,I1,I2,J1,J2,K1,K2,
     &             ISTR,JSTR,KSTR,IN,JN,KN,IOT,IPHASE,NEX)
           NEX = NEX + 3
         ENDIF
         ENDDO
      ENDIF ! MULPHL

C ... Intermittency variables
      IF(TRANSL) THEN
      CALL TROQ(TRM,Q,I1,I2,J1,J2,K1,K2,ISTR,JSTR,KSTR,IN,JN,KN,IOT,NEX)
      NEX = NEX + 2
      ENDIF ! TRANSL

      RETURN
      END SUBROUTINE EXPTOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE RTOQ(R,Q,I1,I2,J1,J2,K1,K2,
     &                ISTR,JSTR,KSTR,IN,JN,KN,IOT,N)

      IMPLICIT NONE

      REAL    :: R(*),Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR,IOT,N
      INTEGER :: I,J,K,IL,IT,I2I1,J2J1,ILK,ITK,ILJ,ITJ,IN,JN,KN
      
      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF

      I2I1 = (I2-I1+1)
      J2J1 = (J2-J1+1)

      DO K=K1,K2,KSTP
         ILK = (K-1+KN)*KSTR+(-1+JN)*JSTR + (-1+IN)*ISTR+1
         ITK = (K-K1)*KSTP*J2J1*I2I1 - J1*I2I1 - I1 + 1
         DO J=J1,J2
            ILJ = ILK + J*JSTR
            ITJ = ITK + J*I2I1
            DO I=I1,I2
               IL = ILJ + I*ISTR
               IT = ITJ + I
               Q(IT+N*IOT)   =  R(IL)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RTOQ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TTOQ(T1,T2,Q,I1,I2,J1,J2,K1,K2,
     &                ISTR,JSTR,KSTR,IN,JN,KN,IOT,IPHASE,N)

C ... Put the phase variales to vector Q

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR,IOT,IPHASE,N
      INTEGER :: I,J,K,IL,IT,I2I1,J2J1,ILK,ITK,ILJ,ITJ,IN,JN,KN
      TYPE(PROPERTIES) T1(*)
      TYPE(MPHASE_VARIABLES) T2(*)
      
      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF

      I2I1 = (I2-I1+1)
      J2J1 = (J2-J1+1)
       
      DO K=K1,K2,KSTP
         ILK = (K-1+KN)*KSTR+(-1+JN)*JSTR + (-1+IN)*ISTR+1
         ITK = (K-K1)*KSTP*J2J1*I2I1 - J1*I2I1 - I1 + 1
         
         DO J=J1,J2
            
            ILJ = ILK + J*JSTR
            ITJ = ITK + J*I2I1
            DO I=I1,I2
               IL = ILJ + I*ISTR
               IT = ITJ + I
               Q(IT+N*IOT)     =  T1(IL)%TEMP(IPHASE)
               Q(IT+(N+1)*IOT) =  T2(IL)%ALFA(IPHASE)
               Q(IT+(N+2)*IOT) =  T2(IL)%X(IPHASE)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TTOQ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TTUQ(T1,T2,Q,I1,I2,J1,J2,K1,K2,
     &                ISTR,JSTR,KSTR,IN,JN,KN,IOT,IPHASE,N)

C ... Put the phase velocities to vector Q

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR,IOT,IPHASE,N
      INTEGER :: I,J,K,IL,IT,I2I1,J2J1,ILK,ITK,ILJ,ITJ,IN,JN,KN
      TYPE(PROPERTIES) T1(*)
      TYPE(MPHASE_VARIABLES) T2(*)
      
      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF

      I2I1 = (I2-I1+1)
      J2J1 = (J2-J1+1)

      DO K=K1,K2,KSTP
         ILK = (K-1+KN)*KSTR+(-1+JN)*JSTR + (-1+IN)*ISTR+1
         ITK = (K-K1)*KSTP*J2J1*I2I1 - J1*I2I1 - I1 + 1
         
         DO J=J1,J2
            
            ILJ = ILK + J*JSTR
            ITJ = ITK + J*I2I1
            DO I=I1,I2
               IL = ILJ + I*ISTR
               IT = ITJ + I
               Q(IT+N*IOT)     =  T2(IL)%U(IPHASE)
               Q(IT+(N+1)*IOT) =  T2(IL)%V(IPHASE)
               Q(IT+(N+2)*IOT) =  T2(IL)%W(IPHASE)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TTUQ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE TROQ(TRM,Q,I1,I2,J1,J2,K1,K2,
     &                ISTR,JSTR,KSTR,IN,JN,KN,IOT,N)

C ... Put the intermittency variales to vector Q

      USE TYPE_ARRAYS

      IMPLICIT NONE

      REAL    :: Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR,IOT,IPHASE,N
      INTEGER :: I,J,K,IL,IT,I2I1,J2J1,ILK,ITK,ILJ,ITJ,IN,JN,KN
      TYPE(INTERMITTENCY) TRM(*)
      
      IF (K2 > K1) THEN
         KSTP = 1
      ELSE
         KSTP = -1
      ENDIF

      I2I1 = (I2-I1+1)
      J2J1 = (J2-J1+1)
       
      DO K=K1,K2,KSTP
         ILK = (K-1+KN)*KSTR+(-1+JN)*JSTR + (-1+IN)*ISTR+1
         ITK = (K-K1)*KSTP*J2J1*I2I1 - J1*I2I1 - I1 + 1
         
         DO J=J1,J2
            
            ILJ = ILK + J*JSTR
            ITJ = ITK + J*I2I1
            DO I=I1,I2
               IL = ILJ + I*ISTR
               IT = ITJ + I
               Q(IT+N*IOT)     =  TRM(IL)%G
               Q(IT+(N+1)*IOT) =  TRM(IL)%RET
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TROQ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE QTOR(R,Q,I1,I2,J1,J2,K1,K2,
     &                ISTR,JSTR,KSTR,IN,JN,KN,IOT,N)

      REAL    :: R(*), Q(*)
      INTEGER :: I1,I2,J1,J2,K1,K2,KSTP,ISTR,JSTR,KSTR,IN,JN,KN
      
      IF (K1 > K2) THEN
         KSTP = -1
      ELSE
         KSTP = 1
      ENDIF

      DO K=K1,K2,KSTP
         DO J=J1,J2
            DO I=I1,I2
               IL = (K-1+KN)*KSTR+(J-1+JN)*JSTR+(I-1+IN)*ISTR+1
               IT = (K-K1)*KSTP*(J2-J1+1)*(I2-I1+1)+
     &              (I-I1)*(J2-J1+1)+(J-J1) + 1
               R(IL)   =  Q(IT+N*IOT)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE QTOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE IWLENG(ILENG,ITURB,NSCAL,ISTRES,MULPHL,NPHASES,TRANSL,
     +   TWO_FLUIDL)

C ... CHECK OUT LENGHT OF THE IO PATHC

      IMPLICIT NONE

      INTEGER :: ILENG,ITURB,NSCAL,ISTRES,IKE,IAL,IEARSM,IMPHASE,
     +           NPHASES,ITRANSL,ITWOFLUID
      LOGICAL :: MULPHL, TRANSL, TWO_FLUIDL

      IKE = 0
      IAL = 0
      IEARSM = 0
      IMPHASE= 0
      ITRANSL= 0
      ITWOFLUID = 0	! Matti Palin 25.09.2017

      IF (ITURB >= 3 .AND. ITURB /= 8) IKE = 2
      IF (ITURB >= 10.AND.ITURB <= 19 ) IAL = 6
      IF(ISTRES >= 1) IEARSM = 12
      IF(MULPHL) IMPHASE = 3
      IF(TWO_FLUIDL) ITWOFLUID = 3
      IF(TRANSL) ITRANSL = 2
      ILENG = 10 + IKE + NSCAL + IAL + IEARSM + IMPHASE*NPHASES +ITRANSL
     +      + ITWOFLUID*NPHASES

      RETURN
      END SUBROUTINE IWLENG
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRVOL(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,FI,
     +   PRO,VAR,TRM,MAXSB,NSCAL,ITURB,MULPHL,TRANSL,A1X,A1Y,A1Z,A2X,
     +   A2Y,A2Z,A3X,A3Y,A3Z,ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,XVEL,
     +   ISCL,IWAVEB,FRESUL,M,MULPHC,PRC,IPRESC)

      USE TYPE_ARRAYS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: MAXSB,NSCAL,ITURB,NPATCH,IMAX,JMAX,KMAX,IPRESC,
     +           IN,JN,KN,KBEGIN,KSCAL,ISTR,JSTR,KSTR,L,IA,IM,JA,JM,
     +           IA2,IM2,JA2,JM2,KM,IFACE,ISCL,IBTYPE,M

      INTEGER :: ICON(IC9,*), IWAVEB(*)

      REAL :: XVEL

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),FI(MAXSB,*),A1X(*),A1Y(*),A1Z(*),TTS(*),
     +        A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)

      LOGICAL :: MULPHL, FRESUL, TRANSL

      TYPE(PROPERTIES)       :: PRO(*)
      TYPE(MPHASE_VARIABLES) :: VAR(*)
      TYPE(INTERMITTENCY)    :: TRM(*)
      TYPE(PRE_COR)          :: PRC(*)

      CHARACTER(*) :: MULPHC

C
C ... MIRROR THE VECTORS ON THE REQUIRED  BOUNDARIES OF THE BLOCK
C
C ... XVEL is 1 for velocities and -1 for vorticities

C ... IF REYNOLDS STRESSES THERE IS NO NEED TO MIRROR THEM
      KBEGIN = 1
      KSCAL = NSCAL
      IF(NSCAL /= 0) THEN
         IF(ITURB >= 20) THEN
            KBEGIN = 7
            KSCAL = NSCAL - 6
         ENDIF
      ENDIF
      IF(ISCL == 1) KBEGIN = 1

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 4 .OR. ICON(1,L) == 7 .OR. ICON(1,L) == 16) THEN
C ------------------------------------------------------------------

! Singularity added 14.12.2010

      IA     = ICON(4,L)
      IM     = ICON(5,L)
      JA     = ICON(6,L)
      JM     = ICON(7,L)
      KM     = 1
      IFACE  = ICON(3,L)
      IBTYPE = ICON(1,L)

C ... Extended MIR
      
      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)
      
      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IF(IFACE == 4) KM = IMAX
         CALL MIRROR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,PRO,
     +        VAR,TRM,FI(1,KBEGIN),KSCAL,ITURB,MULPHL,TRANSL,MAXSB,
     +        A1X,A1Y,A1Z,IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,JN,KN,IN,
     +        IFACE,XVEL,ISCL,NSCAL,IWAVEB,FRESUL,IBTYPE,MULPHC,
     +        PRC,IPRESC)

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IF(IFACE == 5) KM = JMAX
         CALL MIRROR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,PRO,
     +        VAR,TRM,FI(1,KBEGIN),KSCAL,ITURB,MULPHL,TRANSL,MAXSB,
     +        A2X,A2Y,A2Z,JA,JM,IA,IM,KM,KSTR,ISTR,JSTR,KN,IN,JN,
     +        IFACE,XVEL,ISCL,NSCAL,IWAVEB,FRESUL,IBTYPE,MULPHC,
     +        PRC,IPRESC)

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IF(IFACE == 6) KM = KMAX
         CALL MIRROR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,TTS,PRO,
     +        VAR,TRM,FI(1,KBEGIN),KSCAL,ITURB,MULPHL,TRANSL,MAXSB,
     +        A3X,A3Y,A3Z,IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IN,JN,KN,
     +        IFACE,XVEL,ISCL,NSCAL,IWAVEB,FRESUL,IBTYPE,MULPHC,
     +        PRC,IPRESC)
      ENDIF

C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS
1000  CONTINUE

      RETURN
      END SUBROUTINE MIRVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRREL(UU,UV,UW,VV,VW,WW,FUU,FUV,FUW,FVV,FVW,FWW,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,MAXSB,M)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE
      
      INTEGER :: IA, IM, JA, JM, IA2, IM2, JA2, JM2, KM
      INTEGER :: L, M, NPATCH, NTOT, IFACE, MAXSB
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, ISTR, JSTR, KSTR
      INTEGER :: ICON(IC9,*)
      
      REAL :: UU(*), VV(*), WW(*)
      REAL :: UV(*), UW(*), VW(*)

      REAL :: FUU(*), FVV(*), FWW(*)
      REAL :: FUV(*), FUW(*), FVW(*)

      REAL :: A1X(*), A1Y(*), A1Z(*)
      REAL :: A2X(*), A2Y(*), A2Z(*)
      REAL :: A3X(*), A3Y(*), A3Z(*)

C
C ... MIRROR THE REYNOLDS STRESS TENSOR AND ALSO FOR FREE SURFACE
C
      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR
      NTOT = (IMAX+2*IN)*(JMAX+2*JN)*(KMAX+2*KN)
       
      DO 1000 L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 4 .OR. ICON(1,L) == 13) THEN
C ------------------------------------------------------------------
      IA = ICON(4,L)
      IM = ICON(5,L)
      JA = ICON(6,L)
      JM = ICON(7,L)
      KM = 1

      IFACE = ICON(3,L)

C ... Extended MIR      
      
      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)

      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IF(IFACE == 4) KM = IMAX
         CALL MIRREY(UU,UV,UW,VV,VW,WW,FUU,FUV,FUW,FVV,FVW,FWW,
     +               A1X,A1Y,A1Z,IA,IM,JA,JM,KM,
     +               JSTR,KSTR,ISTR,JN,KN,IN,IFACE,NTOT)
      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IF(IFACE == 5) KM = JMAX
         CALL MIRREY(UU,UV,UW,VV,VW,WW,FUU,FUV,FUW,FVV,FVW,FWW,
     +               A2X,A2Y,A2Z,IA,IM,JA,JM,KM,
     +               ISTR,KSTR,JSTR,IN,KN,JN,IFACE,NTOT)
      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IF(IFACE == 6) KM = KMAX
         CALL MIRREY(UU,UV,UW,VV,VW,WW,FUU,FUV,FUW,FVV,FVW,FWW,
     +               A3X,A3Y,A3Z,IA,IM,JA,JM,KM,
     +               ISTR,JSTR,KSTR,IN,JN,KN,IFACE,NTOT)
      ENDIF

C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS
1000  CONTINUE

      RETURN
      END SUBROUTINE MIRREL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRCOR(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A1,A2,A3,
     +     A1X,A1Y,A1Z,A2X,A2Y,A2Z,
     + A3X,A3Y,A3Z,ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,M)

      USE NS3CO, ONLY : IC9
      
      INTEGER :: ICON(IC9,*)
      REAL    :: VOL(*),D1(*),D2(*),D3(*),DISTW(*),BLANK(*),
     +           A1(*),A2(*),A3(*),A1X(*),A1Y(*),A1Z(*),
     +           A2X(*),A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*)
      REAL    :: XC(*),YC(*),ZC(*)
      INTEGER :: NPATCH, IMAX, JMAX, KMAX, IN, JN, KN, M
      INTEGER :: ISTR, JSTR, KSTR
C
C ... MIRROR THE CENTER POINTS OF THE MESH ON THE REQUIRED  BOUNDARIES 
C ... OF THE BLOCK
C

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH

      IF(ICON(1,L) == 4) THEN

      IA = ICON(4,L)
      IM = ICON(5,L)
      JA = ICON(6,L)
      JM = ICON(7,L)
      KM = 1

C ... Extended MIR
      
      CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)
      
      IFACE = ICON(3,L)

      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IF(IFACE == 4) KM = IMAX
         CALL MIRGRI(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A1,A1X,A1Y,A1Z,
     &        IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,JN,KN,IN,IFACE)
      

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IF(IFACE == 5) KM = JMAX
         CALL MIRGRI(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A2,A2X,A2Y,A2Z,
     &        IA,IM,JA,JM,KM,ISTR,KSTR,JSTR,IN,KN,JN,IFACE)
      

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IF(IFACE == 6) KM = KMAX
         CALL MIRGRI(XC,YC,ZC,VOL,D1,D2,D3,DISTW,BLANK,A3,A3X,A3Y,A3Z,
     &        IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IN,JN,KN,IFACE)
      ENDIF

C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS
1000  CONTINUE

      RETURN
      END SUBROUTINE MIRCOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE MIRCCC(XCO,YCO,ZCO,ICON,NPATCH,IW,
     +     IMAX,JMAX,KMAX,IN,JN,KN,M)

      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN, L, M, KA, KM, NPATCH
      INTEGER :: ISTR, JSTR, KSTR, IFACE
      INTEGER :: IA, IM, JA, JM, IA1, IM1, JA1, JM1, IW1
      INTEGER :: IW(*), ICON(IC9,*)
      REAL    :: XCO(*), YCO(*), ZCO(*)
      integer :: IERI, IERJ, IERK, IEB4, IEB6
C
C ... MIRROR THE CORNER POINTS OF THE MESH ON THE REQUIRED  BOUNDARIES 
C ... OF THE BLOCK. ICON must not have patch extensions.
C

      ISTR =  1
      JSTR =  IMAX+2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH

         IW1 = IW(L)

         IF(ICON(1,L) == 4) THEN

            IFACE = ICON(3,L)
            IA1   = ICON(4,L) 
            IM1   = ICON(5,L) + 1  ! Number of grid points 
            JA1   = ICON(6,L) 
            JM1   = ICON(7,L) + 1  ! Number of grid points

            IF(IFACE == 1) THEN
               JA = IA1
               JM = IM1
               KA = JA1
               KM = JM1
               IM = 1
               CALL MIRGRC(XCO,YCO,ZCO,KA,KM,JA,JM,IM,
     &              KSTR,JSTR,ISTR,KN,JN,IN,1)
            ENDIF

            IF(IFACE == 2) THEN
               IA = IA1
               IM = IM1
               KA = JA1
               KM = JM1
               JM = 1
               CALL MIRGRC(XCO,YCO,ZCO,IA,IM,KA,KM,JM,
     &              ISTR,KSTR,JSTR,IN,KN,JN,2)
            ENDIF

            IF(IFACE == 3) THEN
               IA = IA1
               IM = IM1
               JA = JA1
               JM = JM1
               KM = 1
               CALL MIRGRC(XCO,YCO,ZCO,JA,JM,IA,IM,KM,
     &              JSTR,ISTR,KSTR,JN,IN,KN,3)
            ENDIF
            
            IF(IFACE == 4) THEN
               JA = IA1
               JM = IM1
               KA = JA1
               KM = JM1
               IM = IMAX
               CALL MIRGRC(XCO,YCO,ZCO,JA,JM,KA,KM,IM,
     &              JSTR,KSTR,ISTR,JN,KN,IN,4)
            ENDIF

            IF(IFACE == 5) THEN
               IA = IA1
               IM = IM1
               KA = JA1
               KM = JM1
               JM = JMAX
               CALL MIRGRC(XCO,YCO,ZCO,KA,KM,IA,IM,JM,
     &              KSTR,ISTR,JSTR,KN,IN,JN,5)
            ENDIF

            IF(IFACE == 6) THEN
               IA = IA1
               IM = IM1
               JA = JA1
               JM = JM1
               KM = KMAX
               CALL MIRGRC(XCO,YCO,ZCO,IA,IM,JA,JM,KM,
     &               ISTR,JSTR,KSTR,IN,JN,KN,6)
            ENDIF

         ENDIF

 1000 CONTINUE

      RETURN
      END  SUBROUTINE MIRCCC
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYCVOL(U,V,W,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     + AAA,ICON,NPATCH,IW,IMAX,JMAX,KMAX,IN,JN,KN,M,N,IT,IL,IK)

      USE CHARACTERS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE
      
      INTEGER :: IA, IM, JA, JM, KA, KM, ISTR, JSTR, KSTR
      INTEGER :: IA2, IM2, JA2, JM2, KA2, KM2
      INTEGER :: NPATCH, IMAX, JMAX, KMAX, IN, JN, KN, M, N, L
      INTEGER :: IT, IL, IK, IW1
      INTEGER :: IW(*), ICON(IC9,*)
      REAL    :: U(*), V(*), W(*), AAA(9,*),
     &           A1X(*), A1Y(*), A1Z(*),
     &           A2X(*), A2Y(*), A2Z(*),
     &           A3X(*), A3Y(*), A3Z(*)

C ... CYCLIC VECTORS ON THE REQUIRED BOUNDARIES OF THE BLOCK

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
         IW1     = IW(L)

         IF(ICON(1,L) == 6) THEN

            IF(ICON(3,L) == 1) THEN
               JA = ICON(4,L)
               JM = ICON(5,L)
               KA = ICON(6,L)
               KM = ICON(7,L)
               IM = 1
               CALL PATCHE(ICON,L,M,JA,JM,KA,KM,JA2,JM2,KA2,KM2)
               CALL CYCROT(U,V,W,A1X,A1Y,A1Z,KA2,KM2,JA2,JM2,IM,
     &              KSTR,JSTR,ISTR,KN,JN,IN,AAA(1,L),1)
            ENDIF

            IF(ICON(3,L) == 2) THEN
               IA = ICON(4,L)
               IM = ICON(5,L)
               KA = ICON(6,L)
               KM = ICON(7,L)
               JM = 1
               CALL PATCHE(ICON,L,M,IA,IM,KA,KM,IA2,IM2,KA2,KM2)
               CALL CYCROT(U,V,W,A2X,A2Y,A2Z,IA2,IM2,KA2,KM2,JM,
     &              ISTR,KSTR,JSTR,IN,KN,JN,AAA(1,L),2)
            ENDIF

            IF(ICON(3,L) == 3) THEN
               IA = ICON(4,L)
               IM = ICON(5,L)
               JA = ICON(6,L)
               JM = ICON(7,L)
               KM = 1
               CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)
               CALL CYCROT(U,V,W,A3X,A3Y,A3Z,JA2,JM2,IA2,IM2,KM,
     &              JSTR,ISTR,KSTR,JN,IN,KN,AAA(1,L),3)
            ENDIF

            IF(ICON(3,L) == 4) THEN
               JA = ICON(4,L)
               JM = ICON(5,L)
               KA = ICON(6,L)
               KM = ICON(7,L)
               IM = IMAX
               CALL PATCHE(ICON,L,M,JA,JM,KA,KM,JA2,JM2,KA2,KM2)
               CALL CYCROT(U,V,W,A1X,A1Y,A1Z,JA2,JM2,KA2,KM2,IM,
     &              JSTR,KSTR,ISTR,JN,KN,IN,AAA(1,L),4)
            ENDIF

            IF(ICON(3,L) == 5) THEN
               IA = ICON(4,L)
               IM = ICON(5,L)
               KA = ICON(6,L)
               KM = ICON(7,L)
               JM = JMAX
               CALL PATCHE(ICON,L,M,IA,IM,KA,KM,IA2,IM2,KA2,KM2)
               CALL CYCROT(U,V,W,A2X,A2Y,A2Z,KA2,KM2,IA2,IM2,JM,
     &              KSTR,ISTR,JSTR,KN,IN,JN,AAA(1,L),5)
            ENDIF

            IF(ICON(3,L) == 6) THEN
               IA = ICON(4,L)
               IM = ICON(5,L)
               JA = ICON(6,L)
               JM = ICON(7,L)
               KM = KMAX
               CALL PATCHE(ICON,L,M,IA,IM,JA,JM,IA2,IM2,JA2,JM2)
               CALL CYCROT(U,V,W,A3X,A3Y,A3Z,IA2,IM2,JA2,JM2,KM,
     &              ISTR,JSTR,KSTR,IN,JN,KN,AAA(1,L),6)
            ENDIF

         ENDIF

1000  CONTINUE

      RETURN
      END SUBROUTINE CYCVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE FREVOL(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI,BIJ,
     +     MAXEB,MAXSB,NSCAL,ITURB,A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z,
     +     XC,YC,ZC,FRSDEN,WAVEH,IHF,M,
     +     ICON,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,XVEL,FRSPRE,ISTRES,POLD,
     +     DSURLE,IWAVEB,NGL)

      USE CHARACTERS
      USE NS3CO, ONLY : IC9

      IMPLICIT NONE

      INTEGER :: MAXEB,MAXSB,NSCAL,ITURB,NPATCH,IMAX,JMAX,KMAX,IN,JN,KN,
     +           KBEGIN,KSCAL,ISTRES,ISTR,JSTR,KSTR,L,IA,IM,JA,JM,KM,
     +           IFACE,NGL,IPL,M

      INTEGER :: ICON(IC9,*), IHF(*), IWAVEB(*)

      REAL :: FRSDEN, FRSPRE, XVEL

      REAL :: RO(*),U(*),V(*),W(*),E(*),EPS2(*),VIST(*),P(*),PDIFF(*),
     +        RK(*),REPS(*),BIJ(MAXEB,*),A1X(*),A1Y(*),A1Z(*),A2X(*),
     +        A2Y(*),A2Z(*),A3X(*),A3Y(*),A3Z(*),POLD(*),DSURLE(*),
     +        FI(MAXSB,*)

      REAL :: WAVEH(*), XC(*), YC(*), ZC(*)

C ... XVEL is 1 for velocities and -1 for vorticities

C ... IF REYNOLDS STRESSES THERE IS NO NEED TO free surface THEM
      KBEGIN = 1
      KSCAL = NSCAL
      IF(NSCAL /= 0) THEN
         IF(ITURB >= 20) THEN
            KBEGIN = 7
            KSCAL = NSCAL - 6
         ENDIF
      ENDIF

      ISTR = 1
      JSTR = IMAX +2*IN
      KSTR = (JMAX+2*JN)*JSTR

      DO 1000 L = 1,NPATCH
C ------------------------------------------------------------------
      IF(ICON(1,L) == 13) THEN
C ------------------------------------------------------------------

      IA    = ICON(4,L)
      IM    = ICON(5,L)
      JA    = ICON(6,L)
      JM    = ICON(7,L)
      KM    = 1
      IFACE = ICON(3,L)
      IPL   = ICON(2,L)  ! Proces local patch number

      IF(IFACE == 1 .OR. IFACE == 4) THEN
         IF(IFACE == 4) KM = IMAX
         CALL FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        BIJ,KSCAL,ITURB,MAXEB,MAXSB,A1X,A1Y,A1Z,WAVEH,IPL,IHF,
     +        IA,IM,JA,JM,KM,JSTR,KSTR,ISTR,JN,KN,IN,IFACE,XVEL,
     +        FRSPRE,FRSDEN,ISTRES,POLD,DSURLE,IWAVEB,NGL)

      ELSEIF(IFACE == 2 .OR. IFACE == 5) THEN
         IF(IFACE == 5) KM = JMAX
         CALL FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        BIJ,KSCAL,ITURB,MAXEB,MAXSB,A2X,A2Y,A2Z,WAVEH,IPL,IHF,
*     +        JA,JM,IA,IM,KM,KSTR,ISTR,JSTR,IN,KN,JN,IFACE,XVEL,
     +        JA,JM,IA,IM,KM,KSTR,ISTR,JSTR,KN,IN,JN,IFACE,XVEL,
     +        FRSPRE,FRSDEN,ISTRES,POLD,DSURLE,IWAVEB,NGL)

      ELSEIF(IFACE == 3 .OR. IFACE == 6) THEN
         IF(IFACE == 6) KM = KMAX
         CALL FRESUR(RO,U,V,W,E,EPS2,VIST,P,PDIFF,RK,REPS,FI(1,KBEGIN),
     +        BIJ,KSCAL,ITURB,MAXEB,MAXSB,A3X,A3Y,A3Z,WAVEH,IPL,IHF,
     +        IA,IM,JA,JM,KM,ISTR,JSTR,KSTR,IN,JN,KN,IFACE,XVEL,
     +        FRSPRE,FRSDEN,ISTRES,POLD,DSURLE,IWAVEB,NGL)
      ENDIF


C ------------------------------------------------------------------
      ENDIF                                   ! END OF CONNECTIONS
1000  CONTINUE

      RETURN
      END SUBROUTINE FREVOL
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CYCROT(U,V,W,NX,NY,NZ,I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,
     &                  IN,JN,KN,AAA,IWALL)

      REAL :: U(*), V(*), W(*), NX(*), NY(*), NZ(*), AAA(*)
      REAL :: UCYC1, UCYC2, VCYC1, VCYC2, WCYC1, WCYC2

      INTEGER :: I1,I2,J1,J2,KW,ISTR,JSTR,KSTR,IWALL,KDIR,KOFF
      INTEGER :: IT,IT1,IT2,JL,IL,KB,KB1,KB2

      IF ((IWALL >= 1) .AND. (IWALL <= 3)) THEN
         KDIR = -1
         KOFF = 0
      ELSEIF ((IWALL >= 4) .AND. (IWALL <= 6)) THEN
         KDIR = 1
         KOFF = 1
      ELSE
         WRITE(*,*) IWALL
         STOP 'Illegal wall specification for subroutine CYCROT !'
      ENDIF

      KB  = (KW +   KOFF - 1 + KN)*KSTR
      KB1 = (KW +   KDIR - 1 + KN)*KSTR
      KB2 = (KW + 2*KDIR - 1 + KN)*KSTR

      DO 100 J=J1,J2
         JL = (J-1+JN)*JSTR
         DO 200 I=I1,I2
            IL  = JL + (I-1+IN)*ISTR
            IT  = 1 + IL + KB
            IT1 = 1 + IL + KB1
            IT2 = 1 + IL + KB2

            A11 = AAA(1)
            A12 = AAA(2)
            A13 = AAA(3)
            A21 = AAA(4)
            A22 = AAA(5)
            A23 = AAA(6)
            A31 = AAA(7)
            A32 = AAA(8)
            A33 = AAA(9)

            UCYC1 = U(IT1)
            VCYC1 = V(IT1)
            WCYC1 = W(IT1)

            UCYC2 = U(IT2)
            VCYC2 = V(IT2)
            WCYC2 = W(IT2)

            U(IT1) = A11*UCYC1 + A12*VCYC1 + A13*WCYC1
            V(IT1) = A21*UCYC1 + A22*VCYC1 + A23*WCYC1
            W(IT1) = A31*UCYC1 + A32*VCYC1 + A33*WCYC1

            U(IT2) = A11*UCYC2 + A12*VCYC2 + A13*WCYC2
            V(IT2) = A21*UCYC2 + A22*VCYC2 + A23*WCYC2
            W(IT2) = A31*UCYC2 + A32*VCYC2 + A33*WCYC2

 200     CONTINUE
 100  CONTINUE

      RETURN
      END SUBROUTINE CYCROT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PATCHO(ZZZ,APU,XCOL,XCOG,AAA,NBCS,ICON,NPATCH,ROTANG,
     +     NPROCE,IPRO,MBPRO,CHLREF)

C ... Patch orientations in case of connectivity boundary condition
C ... and transformation matrices in case of cyclic boundary
C ... condition.

      USE MAIN_ARRAYS, ONLY : OMEGA, OMEGAX, OMEGAY, OMEGAZ,
     &                        CENAX, CENAY, CENAZ
      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER :: IP,IPRO,NBG,ITYPE,JP,II,JJ,NBGCON,IORI,NBCS,MBPRO
      INTEGER :: ICON(IC9,*), NPATCH(*)
      REAL    :: CHLREF,ROTA1,ROTA2,SLIDE,X1JP,Y1JP,Z1JP,X2JP,Y2JP,Z2JP,
     &           X3JP,Y3JP,Z3JP,X4JP,Y4JP,Z4JP,
     &           A11,A12,A13,A21,A22,A23,A31,A32,A33
      REAL    :: AAA(9,*), ZZZ(*), APU(*)
      REAL    :: XCOL(3,4, NBCS), XCOG(3,4,NBCS), XG(3), XL(3)
      REAL    :: ROTANG(*)
      INTEGER :: NPROCE(MBPRO+1,*)

      REAL, DIMENSION(3) :: RA 
      REAL :: RAL

      CALL CONCOR(ZZZ,APU,XCOL,XCOG,ICON,IPRO,IP,NPATCH,NBCS)

      DO 10 IP=1,NBCS  ! Loop over all patches

         NBG   = ICON(24,IP)                      ! global block number
         ITYPE = ICON(1,IP)
         JP    = ICON(16,IP)
C ... for extend patches the orientation for mirror or free surface 
C ... should be 3. 11.10.96 PPR (tai sit Esa).
         IF(ITYPE == 4 .OR. ITYPE == 12 .OR. ITYPE == 13) ICON(14,IP)=3

         IF (ITYPE ==  1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &       ITYPE == 14 .OR. ITYPE == 15) THEN

         IF(ITYPE == 1 .AND. ICON(20,IP) == 1) THEN ! periodic

            XG(1) = XCOG(1,1,IP)+XCOG(1,2,IP)+XCOG(1,3,IP)+XCOG(1,4,IP)
            XG(2) = XCOG(2,1,IP)+XCOG(2,2,IP)+XCOG(2,3,IP)+XCOG(2,4,IP)
            XG(3) = XCOG(3,1,IP)+XCOG(3,2,IP)+XCOG(3,3,IP)+XCOG(3,4,IP)
            
            XL(1) = XCOL(1,1,IP)+XCOL(1,2,IP)+XCOL(1,3,IP)+XCOL(1,4,IP)
            XL(2) = XCOL(2,1,IP)+XCOL(2,2,IP)+XCOL(2,3,IP)+XCOL(2,4,IP)
            XL(3) = XCOL(3,1,IP)+XCOL(3,2,IP)+XCOL(3,3,IP)+XCOL(3,4,IP)

            DO II = 1,3
            DO JJ = 1,4
               XCOG(II,JJ,IP) = XCOG(II,JJ,IP) + .25*(XL(II)-XG(II))
            ENDDO
            ENDDO
         ENDIF ! periodic

        IF(ITYPE == 1 .OR. ITYPE == 15) THEN   ! Connectivity BC
          NBGCON= NPROCE(1+ICON(8,IP),ICON(17,IP)) ! g. connectivy block n.
          CALL ORIENT(XCOL(1,1,IP),XCOG(1,1,IP),IP,IORI,ICON(22,IP),
     +    NBG,NBGCON,CHLREF)
          ICON(14,IP) = IORI
        ENDIF

        IF(ITYPE == 6 .OR. ITYPE == 11 .OR.
     +     ITYPE == 14) THEN   ! Cyclic BC or SLD or MIX

C ... Solve the cyclicity matrix or the SLD- or MIX-matrix

        IF(ITYPE == 6) THEN
          CALL AAAMAT(XCOL(1,1,IP),XCOG(1,1,IP),AAA,NBCS,IP,IORI)
        ELSEIF(ITYPE == 11 .OR. ITYPE == 14) THEN

           RA(1) = OMEGAX(NBG) ! Rotation axis
           RA(2) = OMEGAY(NBG)
           RA(3) = OMEGAZ(NBG)
           RAL = SQRT(RA(1)**2 + RA(2)**2 + RA(3)**2)

           NBGCON = NPROCE(1+ICON(8,IP),ICON(17,IP)) ! g. connectivy block n.

           IF(OMEGA(NBG) == 0.0) THEN
              RA(1) = OMEGAX(NBGCON) 
              RA(2) = OMEGAY(NBGCON)
              RA(3) = OMEGAZ(NBGCON)
              RAL = SQRT(RA(1)**2 + RA(2)**2 + RA(3)**2)
           ENDIF

           RA = RA/RAL
            
           ROTA1  = ROTANG(NBG) ! ANGLE OF THIS BLOCK
           ROTA2  = ROTANG(NBGCON) ! ANGLE OF CONNECTIVE BLOCK
           SLIDE  = ROTA2-ROTA1 ! SLIDE ANGLE
           CALL SLDMAT(AAA,SLIDE,RA,IP)
        ENDIF

C ... Solve the path orientation in cyclic or sliding case. Temporarily
C ... rotate the patch. 

          X1JP = XCOG(1,1,IP)
          Y1JP = XCOG(2,1,IP)
          Z1JP = XCOG(3,1,IP)
          X2JP = XCOG(1,2,IP)
          Y2JP = XCOG(2,2,IP)
          Z2JP = XCOG(3,2,IP)
          X3JP = XCOG(1,3,IP)
          Y3JP = XCOG(2,3,IP)
          Z3JP = XCOG(3,3,IP)
          X4JP = XCOG(1,4,IP)
          Y4JP = XCOG(2,4,IP)
          Z4JP = XCOG(3,4,IP)

          A11  = AAA(1,IP)
          A12  = AAA(2,IP)
          A13  = AAA(3,IP)
          A21  = AAA(4,IP)
          A22  = AAA(5,IP)
          A23  = AAA(6,IP)
          A31  = AAA(7,IP)
          A32  = AAA(8,IP)
          A33  = AAA(9,IP)

          XCOG(1,1,IP) = A11*X1JP + A12*Y1JP + A13*Z1JP
          XCOG(2,1,IP) = A21*X1JP + A22*Y1JP + A23*Z1JP
          XCOG(3,1,IP) = A31*X1JP + A32*Y1JP + A33*Z1JP
          XCOG(1,2,IP) = A11*X2JP + A12*Y2JP + A13*Z2JP
          XCOG(2,2,IP) = A21*X2JP + A22*Y2JP + A23*Z2JP
          XCOG(3,2,IP) = A31*X2JP + A32*Y2JP + A33*Z2JP
          XCOG(1,3,IP) = A11*X3JP + A12*Y3JP + A13*Z3JP
          XCOG(2,3,IP) = A21*X3JP + A22*Y3JP + A23*Z3JP
          XCOG(3,3,IP) = A31*X3JP + A32*Y3JP + A33*Z3JP
          XCOG(1,4,IP) = A11*X4JP + A12*Y4JP + A13*Z4JP
          XCOG(2,4,IP) = A21*X4JP + A22*Y4JP + A23*Z4JP
          XCOG(3,4,IP) = A31*X4JP + A32*Y4JP + A33*Z4JP

          CALL ORIENT(XCOL(1,1,IP),XCOG(1,1,IP),IP,IORI,ICON(22,IP),
     +    NBG,NBGCON,CHLREF)

          ICON(14,IP) = IORI

          XCOG(1,1,IP) = X1JP
          XCOG(2,1,IP) = Y1JP
          XCOG(3,1,IP) = Z1JP
          XCOG(1,2,IP) = X2JP
          XCOG(2,2,IP) = Y2JP
          XCOG(3,2,IP) = Z2JP
          XCOG(1,3,IP) = X3JP
          XCOG(2,3,IP) = Y3JP
          XCOG(3,3,IP) = Z3JP
          XCOG(1,4,IP) = X4JP
          XCOG(2,4,IP) = Y4JP
          XCOG(3,4,IP) = Z4JP

        ENDIF
        ENDIF !(ITYPE == 1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR. ITYPE == 14 
!                                                                ITYPE == 15)

 10   CONTINUE

      RETURN
      END SUBROUTINE PATCHO
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AAAMAT(XCOL,XCOG,AAA,NBCS,IP,IORI)

C ... This subroutine solves the patch transformation matrix in case
C ... of cyclic boundary condition.

      REAL :: AAA(9,*)
      REAL :: XCOL(3,4), XCOG(3,4)

C ... Surface normal for local patch

      PX  = XCOL(1,3) - XCOL(1,1)
      PY  = XCOL(2,3) - XCOL(2,1)
      PZ  = XCOL(3,3) - XCOL(3,1)
      QX  = XCOL(1,4) - XCOL(1,2)
      QY  = XCOL(2,4) - XCOL(2,2)
      QZ  = XCOL(3,4) - XCOL(3,2)
      RX1 = PY*QZ - PZ*QY
      RY1 = PZ*QX - PX*QZ
      RZ1 = PX*QY - PY*QX
      DR1 = SQRT(RX1*RX1+RY1*RY1+RZ1*RZ1)
      RX1 = RX1/DR1
      RY1 = RY1/DR1
      RZ1 = RZ1/DR1

C ... Centre of the patch

      XI5 = (XCOL(1,1)+XCOL(1,2)+XCOL(1,3)+XCOL(1,4))/4.0
      YI5 = (XCOL(2,1)+XCOL(2,2)+XCOL(2,3)+XCOL(2,4))/4.0
      ZI5 = (XCOL(3,1)+XCOL(3,2)+XCOL(3,3)+XCOL(3,4))/4.0

C ... Surface normal for connectivity patch

      PX  = XCOG(1,3) - XCOG(1,1)
      PY  = XCOG(2,3) - XCOG(2,1)
      PZ  = XCOG(3,3) - XCOG(3,1)
      QX  = XCOG(1,2) - XCOG(1,4)
      QY  = XCOG(2,2) - XCOG(2,4)
      QZ  = XCOG(3,2) - XCOG(3,4)
      RX2 = PY*QZ - PZ*QY
      RY2 = PZ*QX - PX*QZ
      RZ2 = PX*QY - PY*QX
      DR2 = SQRT(RX2*RX2+RY2*RY2+RZ2*RZ2)
      RX2 = RX2/DR2
      RY2 = RY2/DR2
      RZ2 = RZ2/DR2

C ... Centre of the patch

      XJ5 = (XCOG(1,1)+XCOG(1,2)+XCOG(1,3)+XCOG(1,4))/4.0
      YJ5 = (XCOG(2,1)+XCOG(2,2)+XCOG(2,3)+XCOG(2,4))/4.0
      ZJ5 = (XCOG(3,1)+XCOG(3,2)+XCOG(3,3)+XCOG(3,4))/4.0
      

      CX1 = XI5 + RX1
      CY1 = YI5 + RY1
      CZ1 = ZI5 + RZ1

      CX2 = XJ5 + RX2
      CY2 = YJ5 + RY2
      CZ2 = ZJ5 + RZ2

      DX1 = XJ5 - XI5
      DY1 = YJ5 - YI5    ! Vector between patch centre points
      DZ1 = ZJ5 - ZI5
      DD1 = SQRT(DX1*DX1 + DY1*DY1 + DZ1*DZ1)

      DX2 = CX2 - CX1
      DY2 = CY2 - CY1    ! Vector between patch normal heads
      DZ2 = CZ2 - CZ1

      DOT = ABS(RX1*RX2 + RY1*RY2 + RZ1*RZ2 - 1.)

      IF(DD1 < 1.E-6 .OR. DOT < 1.E-6) THEN  ! Cyclicity angle = 0
 
          A = 1.0
          B = 0.0
          C = 0.0
          THETA = 0.0

      ELSE

C ... Rotation axis vector

        A = DY1*DZ2 - DZ1*DY2
        B = DZ1*DX2 - DX1*DZ2
        C = DX1*DY2 - DY1*DX2
        D = SQRT(A*A+B*B+C*C)

        IF(D < 1E-06) THEN
          WRITE(*,*) 'AAAMAT: Cannot solve the cyclicity matrix ',
     &               'of patch',IP
          WRITE(*,*) '        Aborting ...'
          WRITE(45,*) 'AAAMAT: Cannot solve the cyclicity matrix ',
     &                'of patch',IP
          WRITE(45,*) '        Aborting ...'
          STOP
        ENDIF

        EX1 = B*RZ1 - C*RY1
        EY1 = C*RX1 - A*RZ1
        EZ1 = A*RY1 - B*RX1
        DE1 = SQRT(EX1*EX1+EY1*EY1+EZ1*EZ1)

        EX2 = B*RZ2 - C*RY2
        EY2 = C*RX2 - A*RZ2
        EZ2 = A*RY2 - B*RX2
        DE2 = SQRT(EX2*EX2+EY2*EY2+EZ2*EZ2)

        IF(DE1 < 1E-06 .OR. DE2 < 1E-06) THEN
          WRITE(*,*) 'AAAMAT: Cannot solve the cyclicity matrix ',
     &               'of patches',IP
          WRITE(*,*) '        Aborting ...'
          WRITE(45,*) 'AAAMAT: Cannot solve the cyclicity matrix ',
     &                'of patches',IP
          WRITE(45,*) '        Aborting ...'
          STOP
        ENDIF

        EX1 = EX1/DE1
        EY1 = EY1/DE1
        EZ1 = EZ1/DE1

        EX2 = EX2/DE2
        EY2 = EY2/DE2
        EZ2 = EZ2/DE2

        AA = EY1*EZ2 - EZ1*EY2
        BB = EZ1*EX2 - EX1*EZ2
        CC = EX1*EY2 - EY1*EX2
        DD = SQRT(AA*AA+BB*BB+CC*CC)

        IF(DD < 1E-06) THEN
          A = 1.0
          B = 0.0
          C = 0.0
          D = 1.0
        ELSE
          A = AA
          B = BB
          C = CC
          D = DD
        ENDIF

        A = A/D
        B = B/D
        C = C/D

C ... Rotation angle

        DOTVE = EX1*EX2+EY1*EY2+EZ1*EZ2
        ANGLE = SIGN(MIN(ABS(DOTVE),1.0),DOTVE)
        THETA = ACOS(ANGLE)

      ENDIF
         
C ... Cyclicity matrix

      COSTHE = COS(THETA)
      SINTHE = SIN(THETA)

      AAA(1,IP) = A*A*(1.-COSTHE) +   COSTHE
      AAA(2,IP) = B*A*(1.-COSTHE) + C*SINTHE
      AAA(3,IP) = C*A*(1.-COSTHE) - B*SINTHE
      AAA(4,IP) = A*B*(1.-COSTHE) - C*SINTHE
      AAA(5,IP) = B*B*(1.-COSTHE) +   COSTHE
      AAA(6,IP) = C*B*(1.-COSTHE) + A*SINTHE
      AAA(7,IP) = A*C*(1.-COSTHE) + B*SINTHE
      AAA(8,IP) = B*C*(1.-COSTHE) - A*SINTHE
      AAA(9,IP) = C*C*(1.-COSTHE) +   COSTHE

      RETURN
      END SUBROUTINE AAAMAT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SLDMAT(AAA,SLIDE,RA,IP)

C ... Rotation matrix for sliding patches.

      IMPLICIT NONE

      REAL, DIMENSION(3) :: RA

      REAL :: AAA(9,*)
      REAL :: SLIDE, A, B, C, COSTHE, SINTHE

      INTEGER :: IP

      A = RA(1)   ! Rotation axis
      B = RA(2)
      C = RA(3)

      COSTHE = COS(SLIDE)
      SINTHE = SIN(SLIDE)

      AAA(1,IP) = A*A*(1.-COSTHE) +   COSTHE
      AAA(2,IP) = B*A*(1.-COSTHE) + C*SINTHE
      AAA(3,IP) = C*A*(1.-COSTHE) - B*SINTHE
      AAA(4,IP) = A*B*(1.-COSTHE) - C*SINTHE
      AAA(5,IP) = B*B*(1.-COSTHE) +   COSTHE
      AAA(6,IP) = C*B*(1.-COSTHE) + A*SINTHE
      AAA(7,IP) = A*C*(1.-COSTHE) + B*SINTHE
      AAA(8,IP) = B*C*(1.-COSTHE) - A*SINTHE
      AAA(9,IP) = C*C*(1.-COSTHE) +   COSTHE

      RETURN
      END SUBROUTINE SLDMAT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ORIENT(XCOL,XCOG,IP,IORI,NONMAT,NBG,NBGCON,CHLREF)

C ... Solve the orientation of patches in patch connectivity. Solution
C ... is based on the directions of the patch diagonals.

C     4 ----- 3     2 ----- 3   3 ----- 4   4 ----- 1   1 ----- 2
C     |       |     |       |   |       |   |       |   |       |
C     |   A   |     |   B   |   |   B   |   |   B   |   |   B   |
C     |       |     |       |   |       |   |       |   |       |
C     1 ----- 2     1 ----- 4   2 ----- 1   3 ----- 2   4 ----- 3
C
C                   IORI = 0    IORI = 1    IORI = 2    IORI = 3

      REAL :: XCOL(3,4), XCOG(3,4)

      A13X = XCOL(1,3) - XCOL(1,1)
      A13Y = XCOL(2,3) - XCOL(2,1)
      A13Z = XCOL(3,3) - XCOL(3,1)

      A24X = XCOL(1,4) - XCOL(1,2)
      A24Y = XCOL(2,4) - XCOL(2,2)
      A24Z = XCOL(3,4) - XCOL(3,2)

      B13X = XCOG(1,3) - XCOG(1,1)
      B13Y = XCOG(2,3) - XCOG(2,1)
      B13Z = XCOG(3,3) - XCOG(3,1)

      B24X = XCOG(1,4) - XCOG(1,2)
      B24Y = XCOG(2,4) - XCOG(2,2)
      B24Z = XCOG(3,4) - XCOG(3,2)

      B31X = XCOG(1,1) - XCOG(1,3)
      B31Y = XCOG(2,1) - XCOG(2,3)
      B31Z = XCOG(3,1) - XCOG(3,3)

      B42X = XCOG(1,2) - XCOG(1,4)
      B42Y = XCOG(2,2) - XCOG(2,4)
      B42Z = XCOG(3,2) - XCOG(3,4)

      A13B13 = SQRT((A13X-B13X)**2 + (A13Y-B13Y)**2 + (A13Z-B13Z)**2)
      A24B42 = SQRT((A24X-B42X)**2 + (A24Y-B42Y)**2 + (A24Z-B42Z)**2)
      S0     = A13B13 + A24B42

      A13B24 = SQRT((A13X-B24X)**2 + (A13Y-B24Y)**2 + (A13Z-B24Z)**2)
      A24B13 = SQRT((A24X-B13X)**2 + (A24Y-B13Y)**2 + (A24Z-B13Z)**2)
      S1     = A13B24 + A24B13

      A13B31 = SQRT((A13X-B31X)**2 + (A13Y-B31Y)**2 + (A13Z-B31Z)**2)
      A24B24 = SQRT((A24X-B24X)**2 + (A24Y-B24Y)**2 + (A24Z-B24Z)**2)
      S2     = A13B31 + A24B24

      A13B42 = SQRT((A13X-B42X)**2 + (A13Y-B42Y)**2 + (A13Z-B42Z)**2)
      A24B31 = SQRT((A24X-B31X)**2 + (A24Y-B31Y)**2 + (A24Z-B31Z)**2)
      S3     = A13B42 + A24B31


      IORI = -1
      SS = 1.E-3*CHLREF   ! New value for testing
      IF(S0 < SS) THEN
        IORI = 0
        SS   = S0
      ENDIF
      IF(S1 < SS) THEN
        IORI = 1
        SS   = S1
      ENDIF
      IF(S2 < SS) THEN
        IORI = 2
        SS   = S2
      ENDIF
      IF(S3 < SS) THEN
        IORI = 3
        SS   = S3
      ENDIF
      write(45,*) 'S0,S1,S2,S3 =',
     +             REAL(S0,4),REAL(S1,4),REAL(S2,4),REAL(S3,4)

c     print *
c      print *,A13B13,A24B42,S0
c      print *,A13B24,A24B13,S1
c      print *,A13B31,A24B24,S2
c      print *,A13B42,A24B31,S3
c      print *,IP,IORI

      IF(IORI < 0 .AND. NONMAT /= 1) THEN
        WRITE(*,*)
        WRITE(*,*)  'ORIENT: Cannot solve the patch orientation ????'
        WRITE(*,100) IP
        WRITE(45,*)
        WRITE(45,*) 'ORIENT: Cannot solve the patch orientation ????'
        DO J = 1,4
           WRITE(45,*) 'Corner',J,' Iori=',IORI
           WRITE(45,*) 'Local x,y,z',(XCOL(I,J),I=1,3)
           WRITE(45,*) 'Conne x,y,z',(XCOG(I,J),I=1,3)
           WRITE(45,*)
        ENDDO
        WRITE(45,101) NBG,NBGCON
        WRITE(45,100) IP
        STOP
      ELSEIF(IORI < 0 .AND. NONMAT == 1) THEN
        IORI = 0  ! Non-matching connection
      ENDIF

100   FORMAT('          Connectivity patch',I4,/
     +       '          Aborting...')
101   FORMAT('          Between the blocks',I3,' and ',I3)

      RETURN
      END SUBROUTINE ORIENT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CONCOR(ZZZ,APU,XCOL,XCOG,ICON,IPRO,IP,NPATCH,NBCS)

C ... THIS SUBROUTINE SOLVES CORNER POINTS OF TWO CONNECTIONS

      USE NS3CO, ONLY : IC9
      
      REAL :: XCOL(3,4,NBCS), XCOG(3,4,NBCS), ZZZ(*), APU(*)
 
      INTEGER :: ICON(IC9,*), NPATCH(*)

      DO 10 IP=1,NBCS  ! Loop over all patches
         ITYPE = ICON(1,IP)
         JP    = ICON(16,IP)
         IW1   = (IP-1)*12 + 1
         IOTY  = ICON(18,IP)
         JPRO  = ICON(17,IP)
         IF (ITYPE == 1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &       ITYPE == 14 .OR. ITYPE == 15) THEN
 
C ... CONNECTIVE CORNERS
            ZZZ(IW1   ) = XCOL(1,1,IP)
            ZZZ(IW1+ 1) = XCOL(2,1,IP)
            ZZZ(IW1+ 2) = XCOL(3,1,IP)
            ZZZ(IW1+ 3) = XCOL(1,2,IP)
            ZZZ(IW1+ 4) = XCOL(2,2,IP)
            ZZZ(IW1+ 5) = XCOL(3,2,IP)
            ZZZ(IW1+ 6) = XCOL(1,3,IP)
            ZZZ(IW1+ 7) = XCOL(2,3,IP)
            ZZZ(IW1+ 8) = XCOL(3,3,IP)
            ZZZ(IW1+ 9) = XCOL(1,4,IP)
            ZZZ(IW1+10) = XCOL(2,4,IP)
            ZZZ(IW1+11) = XCOL(3,4,IP)
            IF(IOTY >= 2) CALL SEFILE(ZZZ(IW1),JPRO,JP,12,IOTY)
         ENDIF
 10   CONTINUE
      DO 20 IP=1,NBCS  ! Loop over all patches
         ITYPE = ICON(1,IP)
         JP    = ICON(16,IP)
         IW1   = (JP-1)*12 + 1
         IOTY  = ICON(18,IP)
         JPRO  = ICON(17,IP)
         IF (ITYPE == 1 .OR. ITYPE == 6 .OR. ITYPE == 11 .OR.
     &       ITYPE == 14 .OR. ITYPE == 15) THEN

            IF(IOTY == 1) CALL SETV12(APU,ZZZ(IW1),12)
            IF(IOTY >= 2) CALL REFILE(APU,IPRO,JPRO,IP,12,IOTY)
            XCOG(1,1,IP) = APU( 1)
            XCOG(2,1,IP) = APU( 2)
            XCOG(3,1,IP) = APU( 3)
            XCOG(1,2,IP) = APU( 4)
            XCOG(2,2,IP) = APU( 5)
            XCOG(3,2,IP) = APU( 6)
            XCOG(1,3,IP) = APU( 7)
            XCOG(2,3,IP) = APU( 8)
            XCOG(3,3,IP) = APU( 9)
            XCOG(1,4,IP) = APU(10)
            XCOG(2,4,IP) = APU(11)
            XCOG(3,4,IP) = APU(12)
            IF(IOTY >= 4) STOP 'EI OLE TYYPPIA'
         ENDIF
 20   CONTINUE

      RETURN 
      END SUBROUTINE CONCOR
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE OVERLA(D,E,AP,DIR)      

C ... This subroutine computes the common area of two overlaping
C ... polygons. Polygons have four corners.

      IMPLICIT REAL (A-H,O-Z)

      REAL :: A(3,5),B(3,5),C(4,25),D(3,4),E(3,4),
     +        F(3,3),G(4),XCD,YCD,ZCD,D1,D2,D3,D4

      REAL :: AP, DIR, XCA, YCA, ZCA, XCB, YCB, ZCB, RA, RB

      LOGICAL :: SORTED
     
      EPS = 1.0E-10


C ... Check first if the circles drawn around the polygons overlap
C ... each other or not. If not, we know that the polygons cannot
C ... have common area and we can just return to the calling routine.

      XCD = (D(1,1) + D(1,2) + D(1,3) + D(1,4))/4.0
      YCD = (D(2,1) + D(2,2) + D(2,3) + D(2,4))/4.0
      ZCD = (D(3,1) + D(3,2) + D(3,3) + D(3,4))/4.0

      XCA = XCD
      YCA = YCD
      ZCA = ZCD

      D1 = (D(1,1)-XCD)**2 + (D(2,1)-YCD)**2 + (D(3,1)-ZCD)**2 
      D2 = (D(1,2)-XCD)**2 + (D(2,2)-YCD)**2 + (D(3,2)-ZCD)**2 
      D3 = (D(1,3)-XCD)**2 + (D(2,3)-YCD)**2 + (D(3,3)-ZCD)**2
      D4 = (D(1,4)-XCD)**2 + (D(2,4)-YCD)**2 + (D(3,4)-ZCD)**2

      RA  = SQRT(MAX(D1,D2,D3,D4))

      XCD = (E(1,1) + E(1,2) + E(1,3) + E(1,4))/4.0
      YCD = (E(2,1) + E(2,2) + E(2,3) + E(2,4))/4.0
      ZCD = (E(3,1) + E(3,2) + E(3,3) + E(3,4))/4.0

      XCB = XCD
      YCB = YCD
      ZCB = ZCD

      D1 = (E(1,1)-XCD)**2 + (E(2,1)-YCD)**2 + (E(3,1)-ZCD)**2 
      D2 = (E(1,2)-XCD)**2 + (E(2,2)-YCD)**2 + (E(3,2)-ZCD)**2 
      D3 = (E(1,3)-XCD)**2 + (E(2,3)-YCD)**2 + (E(3,3)-ZCD)**2
      D4 = (E(1,4)-XCD)**2 + (E(2,4)-YCD)**2 + (E(3,4)-ZCD)**2

      RB  = SQRT(MAX(D1,D2,D3,D4))
 

      IF(RA+RB < SQRT((XCA-XCB)**2+(YCA-YCB)**2+(ZCA-ZCB)**2)) THEN
        AP = 0.0
        RETURN
      ENDIF


C ... If the polygons do not lie exactly on the same plane, 
C ... they must be first projected onto the mean plane
C ... between the polygons.


C ... Normal of polygon D  

      CX = D(1,3) - D(1,1)
      CY = D(2,3) - D(2,1)
      CZ = D(3,3) - D(3,1)

      DX = D(1,4) - D(1,2)
      DY = D(2,4) - D(2,2)
      DZ = D(3,4) - D(3,2)

      EX = CY*DZ - CZ*DY
      EY = CZ*DX - CX*DZ
      EZ = CX*DY - CY*DX
      EE = SQRT(EX**2 + EY**2 + EZ**2)
      EX = EX/EE
      EY = EY/EE
      EZ = EZ/EE


C ... Normal of polygon E  

      CX = E(1,3) - E(1,1)
      CY = E(2,3) - E(2,1)
      CZ = E(3,3) - E(3,1)

      DX = E(1,4) - E(1,2)
      DY = E(2,4) - E(2,2)
      DZ = E(3,4) - E(3,2)

      FX = CY*DZ - CZ*DY
      FY = CZ*DX - CX*DZ
      FZ = CX*DY - CY*DX
      FF = SQRT(FX**2 + FY**2 + FZ**2)
      FX = FX/FF
      FY = FY/FF
      FZ = FZ/FF

      IF(EX*DIR*FX+EY*DIR*FY+EZ*DIR*FZ < 0.0E+00) THEN
        AP = 0.0
        RETURN
      ENDIF


C ... Sum of the normals (mean normal)

      GX = EX + DIR*FX 
      GY = EY + DIR*FY 
      GZ = EZ + DIR*FZ
      GG = SQRT(GX**2 + GY**2 + GZ**2) 
      GX = GX/GG
      GY = GY/GG
      GZ = GZ/GG


C ... Centers of the polygons

      CXD = (D(1,1) + D(1,2) + D(1,3) + D(1,4))/4.0
      CYD = (D(2,1) + D(2,2) + D(2,3) + D(2,4))/4.0
      CZD = (D(3,1) + D(3,2) + D(3,3) + D(3,4))/4.0

      CXE = (E(1,1) + E(1,2) + E(1,3) + E(1,4))/4.0
      CYE = (E(2,1) + E(2,2) + E(2,3) + E(2,4))/4.0
      CZE = (E(3,1) + E(3,2) + E(3,3) + E(3,4))/4.0


C ... A point between the polygons (point on the mean plane)

      P1X = (CXD + CXE)/2.0
      P1Y = (CYD + CYE)/2.0
      P1Z = (CZD + CZE)/2.0


C ... An additional point on the plane

      P2X = P1X + CY*GZ - CZ*GY
      P2Y = P1Y + CZ*GX - CX*GZ
      P2Z = P1Z + CX*GY - CY*GX


C ... Project the corner points of the original polygons onto the
C ... mean plane.

      DO K = 1,4
         S = GX*(P2X-D(1,K))+GY*(P2Y-D(2,K))+GZ*(P2Z-D(3,K))
         A(1,K) = D(1,K) + S*GX
         A(2,K) = D(2,K) + S*GY
         A(3,K) = D(3,K) + S*GZ
         S = GX*(P2X-E(1,K))+GY*(P2Y-E(2,K))+GZ*(P2Z-E(3,K))
         B(1,K) = E(1,K) + S*GX
         B(2,K) = E(2,K) + S*GY
         B(3,K) = E(3,K) + S*GZ
      ENDDO


C ... Close the polygons

      A(1,5) = A(1,1)
      A(2,5) = A(2,1)
      A(3,5) = A(3,1)

      B(1,5) = B(1,1)
      B(2,5) = B(2,1)
      B(3,5) = B(3,1)


      NPOINT = 0

C ... Check first if the corner(s) of polygon A lie inside
C ... polygon B and vice versa.

      DO N=1,4
         DO K=1,4 
            L = K + 1 
            AX = B(1,L) - B(1,K)
            AY = B(2,L) - B(2,K)
            AZ = B(3,L) - B(3,K)
            BX = A(1,N) - B(1,K)
            BY = A(2,N) - B(2,K)
            BZ = A(3,N) - B(3,K)
            SX = AY*BZ  - AZ*BY
            SY = AZ*BX  - AX*BZ
            SZ = AX*BY  - AY*BX
            G(K) = GX*SX + GY*SY + GZ*SZ
         ENDDO
         GMIN = MIN(G(1),G(2),G(3),G(4))
         GMAX = MAX(G(1),G(2),G(3),G(4))
         IF(GMIN*GMAX >= 0.0) THEN
            NPOINT = NPOINT + 1
            C(1,NPOINT) = A(1,N)
            C(2,NPOINT) = A(2,N) 
            C(3,NPOINT) = A(3,N)
         ENDIF 
      ENDDO

      DO N=1,4
         DO K=1,4 
            L = K + 1 
            AX = A(1,L) - A(1,K)
            AY = A(2,L) - A(2,K)
            AZ = A(3,L) - A(3,K)
            BX = B(1,N) - A(1,K)
            BY = B(2,N) - A(2,K)
            BZ = B(3,N) - A(3,K)
            SX = AY*BZ  - AZ*BY
            SY = AZ*BX  - AX*BZ
            SZ = AX*BY  - AY*BX
            G(K) = GX*SX + GY*SY + GZ*SZ
         ENDDO
         GMIN = MIN(G(1),G(2),G(3),G(4))
         GMAX = MAX(G(1),G(2),G(3),G(4))
         IF(GMIN*GMAX >= 0.0) THEN
          NPOINT = NPOINT + 1
            C(1,NPOINT) = B(1,N)
            C(2,NPOINT) = B(2,N) 
            C(3,NPOINT) = B(3,N)
        ENDIF 
      ENDDO


C ... Search next locations where the polygon edges intersect each other.

      DO N=1,4
        M = N + 1
         DO L=1,4
          K = L + 1

            AX = A(1,L) - A(1,K)                
            AY = A(2,L) - A(2,K)                
            AZ = A(3,L) - A(3,K)                

            BX = B(1,M) - B(1,N)                
            BY = B(2,M) - B(2,N)                
            BZ = B(3,M) - B(3,N)                

            DX = B(1,M) - A(1,K)                
            DY = B(2,M) - A(2,K)
            DZ = B(3,M) - A(3,K)

            UX = DY*BZ  - DZ*BY
            UY = DZ*BX  - DX*BZ
            UZ = DX*BY  - DY*BX
            UL = SQRT(UX**2 + UY**2 + UZ**2)

            VX = AY*BZ  - AZ*BY
            VY = AZ*BX  - AX*BZ
            VZ = AX*BY  - AY*BX
            VL = SQRT(VX**2 + VY**2 + VZ**2)

            PX = B(1,N) - A(1,K) 
            PY = B(2,N) - A(2,K) 
            PZ = B(3,N) - A(3,K) 

            QX = A(1,L) - B(1,N) 
            QY = A(2,L) - B(2,N) 
            QZ = A(3,L) - B(3,N) 

            RX = B(1,M) - A(1,L) 
            RY = B(2,M) - A(2,L) 
            RZ = B(3,M) - A(3,L) 

            SX = A(1,K) - B(1,M) 
            SY = A(2,K) - B(2,M) 
            SZ = A(3,K) - B(3,M) 

            PRQX = PY*QZ - PZ*QY
            PRQY = PZ*QX - PX*QZ
            PRQZ = PX*QY - PY*QX

            QRRX = QY*RZ - QZ*RY
            QRRY = QZ*RX - QX*RZ
            QRRZ = QX*RY - QY*RX

            RRSX = RY*SZ - RZ*SY
            RRSY = RZ*SX - RX*SZ
            RRSZ = RX*SY - RY*SX

            SRPX = SY*PZ - SZ*PY
            SRPY = SZ*PX - SX*PZ
            SRPZ = SX*PY - SY*PX

            G(1) = GX*PRQX + GY*PRQY + GZ*PRQZ           
            G(2) = GX*QRRX + GY*QRRY + GZ*QRRZ           
            G(3) = GX*RRSX + GY*RRSY + GZ*RRSZ           
            G(4) = GX*SRPX + GY*SRPY + GZ*SRPZ           

            GMIN = MIN(G(1),G(2),G(3),G(4))
            GMAX = MAX(G(1),G(2),G(3),G(4))

            IF(GMIN*GMAX > 0.0) THEN
               IF(VL /= 0.0) THEN
                  S = UL/VL
                  IF(S > 0.0 .AND. S < 1.0) THEN
            NPOINT = NPOINT + 1
                     C(1,NPOINT) = A(1,K) + S*AX 
                     C(2,NPOINT) = A(2,K) + S*AY 
                     C(3,NPOINT) = A(3,K) + S*AZ 
                  ENDIF
               ENDIF
            ENDIF                
         ENDDO
      ENDDO


C ... Search the points which were found twice or even more often

      DO IP=1,NPOINT-1
         DO JP=IP+1,NPOINT
            DD = SQRT((C(1,JP)-C(1,IP))**2 + (C(2,JP)-C(2,IP))**2
     &                                      + (C(3,JP)-C(3,IP))**2)
            IF(DD < EPS) THEN
               C(1,JP) = 99999.0 
               C(2,JP) = 99999.0
               C(3,JP) = 99999.0
          ENDIF 
        ENDDO
      ENDDO

      KP = 0
      DO IP=1,NPOINT
         IF(C(1,IP) < 99998.0) THEN 
            KP = KP + 1
            C(1,KP) = C(1,IP)
            C(2,KP) = C(2,IP)
            C(3,KP) = C(3,IP)
         ENDIF
      ENDDO

      NPOINT = KP

      IF(NPOINT < 3) THEN
         AP = 0.0
         RETURN
      ENDIF


C ... Compute the center of the overlaping area.

      XCP = 0.0
      YCP = 0.0
      ZCP = 0.0
      DO NP = 1,NPOINT
         XCP = XCP + C(1,NP)
         YCP = YCP + C(2,NP)
         ZCP = ZCP + C(3,NP)
      ENDDO
      XCP = XCP/NPOINT
      YCP = YCP/NPOINT 
      ZCP = ZCP/NPOINT


C ... Reorder the polygon corner points

      CX = C(1,1) - XCP  ! Vector from center to point 1 (vector a)
      CY = C(2,1) - YCP
      CZ = C(3,1) - ZCP
      CC = SQRT(CX**2 + CY**2 + CZ**2)
      C(4,1) = 0.0

      DO I = 2,NPOINT
         DX = C(1,I) - XCP ! Vector from center to point i (vector b)
         DY = C(2,I) - YCP
         DZ = C(3,I) - ZCP
         DD = SQRT(DX**2 + DY**2 + DZ**2)
         HX = CY*DZ-CZ*DY
         HY = CZ*DX-CX*DZ
         HZ = CX*DY-CY*DX
         HH = GX*HX + GY*HY + GZ*HZ
         AA = MAX((CX*DX+CY*DY+CZ*DZ)/(CC*DD), -1.0)
         AA = MIN(AA, 1.0)
         C(4,I) = SIGN(ACOS(AA),HH)  ! Angle between vectors a and b
      ENDDO


C ... Bubble sort (use angle as the sorting key)

      I = 1
 10   CONTINUE
         I = I + 1
         SORTED = .TRUE.
         DO J=NPOINT,I,-1
            IF(C(4,J-1) > C(4,J)) THEN
               SORTED   = .FALSE.
               XAPU     = C(1,J-1)
               YAPU     = C(2,J-1)
               ZAPU     = C(3,J-1)
               AAPU     = C(4,J-1)
               C(1,J-1) = C(1,J)
               C(2,J-1) = C(2,J)
               C(3,J-1) = C(3,J)
               C(4,J-1) = C(4,J)
               C(1,J)   = XAPU 
               C(2,J)   = YAPU 
               C(3,J)   = ZAPU 
               C(4,J)   = AAPU 
            ENDIF
         ENDDO
         IF(SORTED) GOTO 20
      GOTO 10
 20   CONTINUE


C ... Compute the area of the polygon (overlaping area)

      AP = 0.0
      F(1,1) = C(1,1)
      F(2,1) = C(2,1)
      F(3,1) = C(3,1)
      DO I = 2,NPOINT-1
         F(1,2) = C(1,I)
         F(2,2) = C(2,I)
         F(3,2) = C(3,I)
         F(1,3) = C(1,I+1)
         F(2,3) = C(2,I+1)
         F(3,3) = C(3,I+1)

         AX = F(1,2) - F(1,1)
         AY = F(2,2) - F(2,1)
         AZ = F(3,2) - F(3,1)

         BX = F(1,3) - F(1,1)
         BY = F(2,3) - F(2,1)
         BZ = F(3,3) - F(3,1)
      
         AP = AP + 0.5*SQRT((AY*BZ-AZ*BY)**2
     &                    + (AZ*BX-AX*BZ)**2
     &                    + (AX*BY-AY*BX)**2)

      ENDDO


C ... Project the overlaping area back to the original surface.

      AP = AP/(EX*GX + EY*GY + EZ*GZ)

      RETURN
      END SUBROUTINE OVERLA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ... WALL DISTANCE SUBROUTINES....

      SUBROUTINE DISTANCES(XC,YC,ZC,DISTW,XW,YW,ZW,IMAX,JMAX,KMAX,
     &     NTOT,IG,MGRID,IN,JN,KN,ICON,ICS,NPATCH,NB,MGM,NBLOCK,
     &     LEVEL,IGRI,myid,IMAXG,JMAXG,KMAXG,NLOCAL,NBLOCG,NPNUM,
     &     ntasks,LOCDIS,JET,IREPEA,NREPEA,WDISTS_FORCED)

C ... Subroutine to calculate the distance to the nearest wall.
C ... PK 14.2.1997

      USE MPI
      
      IMPLICIT NONE

      INTEGER :: NB,MGM,ICS,NBLOCK,IN,JN,KN,IGRI

      INTEGER :: STATUS(MPI_STATUS_SIZE)
      
      REAL :: XC(*), YC(*), ZC(*)

      REAL :: DISTW(*), XW(*), YW(*), ZW(*), DISTINI
      INTEGER :: LOCDIS(*), JET(*), IREPEA(*), NREPEA(*)
      REAL, ALLOCATABLE,DIMENSION(:) :: B
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCDIB

      INTEGER :: MGRID(*), NPATCH(*), ICON(ICS,*)
      INTEGER :: IMAX(MGM,NB), JMAX(MGM,NB), KMAX(MGM,NB), NTOT(MGM,NB),
     &           IG(MGM,NB)
      INTEGER :: N,M,NBL,NP,NPBASE,I,ISTR,JSTR,KSTR,MG,ITYPE,IG1,IG2,
     &           LEVEL,ISPEC,IDSRC,ISTAG,NTOTMX,NT,ERCODE
      INTEGER :: IMAXG(*),JMAXG(*),KMAXG(*),NLOCAL(*),myid,ntasks,ierr, 
     &           NBLOCG,NTOTG,NPNUM(*),NG

      CHARACTER(LEN=120) :: RIVI,RIVI2
      CHARACTER(LEN=1)   :: EXT
      CHARACTER(LEN=8)   :: FNAM
      LOGICAL            :: THERE, WDISTS_FORCED
                     
      ERCODE = 0
      MG     = 1
      NPBASE = 1
      NTOTMX = 1
c      DISTINI= HUGE(DISTW)*0.0001
      DISTINI= 1.E+8

C     Initialize wall distance array

      DO NBL=1,NBLOCK
            IG1 = IG(MG,NBL)
            IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
            DISTW(IG1:IG2) = DISTINI
            LOCDIS(IG1:IG2)= 0
                  
C     Calculate amount of extra memory.
C     Note: Stuff of pid=0 doesn't need it
      END DO

      DO NBL = 1,NBLOCG
         NTOTG = (IMAXG(NBL)+2*IN)*(JMAXG(NBL)+2*JN)*(KMAXG(NBL)+2*KN)
         NTOTMX = MAX0(NTOTMX,NTOTG)
      ENDDO

C     Allocate memory for messages

      IF (ntasks > 1) THEN
         ERCODE = 0
         IF (myid == 0) THEN
            ALLOCATE(B(NTOTMX), LOCDIB(NTOTMX), STAT=ERCODE)
         END IF
         CALL MPI_BCAST(ERCODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF (ERCODE /= 0) THEN
            WRITE(*,*) 'DISTANCES: Not enough memory. Aborting.'
            CALL MPI_FINALIZE(ierr)
         ELSE
            IF(IREPEA(2) < NREPEA(2)) WRITE(45,9050)
         END IF
      END IF

9050  FORMAT('In DISTANCES: Allocated memory for messages')

C     Check for a file DISTW_L with pre-calulated values

      EXT  = CHAR(48+LEVEL)
      FNAM = 'DISTW_L'//EXT

      IF (myid == 0) INQUIRE(FILE=FNAM,EXIST=THERE)
      IF (ntasks > 1) THEN
         CALL MPI_BCAST(THERE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      END IF

      IF (myid == 0) THEN
         OPEN(IGRI,FILE=FNAM,STATUS='UNKNOWN',FORM='UNFORMATTED')
      END IF

C     Read the values from DISTW_L* file (Fine grid only)
                 
      IF (IREPEA(2) <= NREPEA(2)) WRITE(*,*)

      IF (THERE .AND. .NOT.WDISTS_FORCED) THEN

         IF (myid == 0) THEN
            RIVI="('  DISTANCES: Reading wall distances from file ',A8)"
            RIVI2="('  MESSAGE: Reading wall distances from file ',
     &      A8,' - repeated more than',I3,' times. Silence')"
            IF(IREPEA(2) < NREPEA(2)) THEN
               WRITE(*,RIVI) FNAM
            ELSE
               IF(IREPEA(2) == NREPEA(2))WRITE(*,RIVI2) FNAM,IREPEA(2)
            ENDIF
         END IF
         DO NG=1,NBLOCG
            IDSRC = NPNUM(NG) - 1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
            IF (myid == 0) THEN
               IF (myid == IDSRC) THEN
                  IG1   = IG(MG,N)-IG(MG,1)
                  READ(IGRI) (DISTW(IG1+I),I=1,NT)
                  READ(IGRI) (LOCDIS(IG1+I),I=1,NT)
               ELSE
                  READ(IGRI) (B(I),I=1,NT)
                  READ(IGRI) (LOCDIB(I),I=1,NT)
                  CALL MPI_SEND(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)
                  CALL MPI_SEND(LOCDIB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)
               END IF
            ELSEIF (myid == IDsrc) THEN
               IG1   = IG(MG,N)-IG(MG,1)
               CALL MPI_RECV(DISTW(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
               CALL MPI_RECV(LOCDIS(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
            END IF
         END DO
      ELSE                      ! We really have to calculate them.
                                ! Better luck next time.
         IF (myid == 0) THEN
            RIVI="('  DISTANCES: Calculating wall distances.')"
            WRITE(*,RIVI)
            RIVI ="('Block ',I5)"
         END IF
         DO NG=1,NBLOCG
            NBL = NLOCAL(NG)
            IF (myid == 0) WRITE(*,RIVI) NG
            IF (NPNUM(NG)-1 == myid) THEN
               ISTR = 1
               JSTR = IMAX(MG,NBL) + 2*IN
               KSTR = (JMAX(MG,NBL) + 2*JN)*JSTR
            ENDIF
C            DO NP = NPBASE+1 , NPBASE+NPATCH(NBL)
C               ITYPE = ICON(1,NP)
C               ISPEC = ICON(16,NP)
            CALL ONE_BLOCK(XC,YC,ZC,DISTW,XW,YW,ZW,IMAX,JMAX,
     &           KMAX,NTOT,IG,MGRID,IN,JN,KN,ICON(1,NPBASE),ICS,NPATCH,
     &           NB,MGM,NBLOCK,NBL,MG,ISTR,JSTR,KSTR,NPNUM,NG,
     &           myid,ntasks,LOCDIS,JET)
C         END DO
            IF (NPNUM(NG)-1 == myid) NPBASE = NPBASE + NPATCH(NBL)
         END DO


         DO NBL=1,NBLOCK
            IG1 = IG(MG,NBL)
            IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
            DISTW(IG1:IG2) = SQRT(DISTW(IG1:IG2))
         END DO


C     Write the values from DISTW_L* file (Fine grid only)

         IF (myid == 0) THEN
            RIVI="('  DISTANCES: Writing wall distances to file ',A8)"
            WRITE(*,RIVI) FNAM
            WRITE(*,*)
         END IF

         DO NG=1,NBLOCG
            IDSRC = NPNUM(NG)-1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
            IF (myid == IDSRC) IG1   = IG(MG,N)-IG(MG,1)
            IF (myid == 0) THEN
               IF (myid == IDSRC) THEN
                  WRITE(IGRI) (DISTW(IG1+I),I=1,NT)
                  WRITE(IGRI) (LOCDIS(IG1+I),I=1,NT)
               ELSE
                  CALL MPI_RECV(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  CALL MPI_RECV(LOCDIB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  WRITE(IGRI) (B(I),I=1,NT)
                  WRITE(IGRI) (LOCDIB(I),I=1,NT)
               END IF
            ELSEIF (myid == IDsrc) THEN
               CALL MPI_SEND(DISTW(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,ierr)
               CALL MPI_SEND(LOCDIS(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,ierr)
            END IF
         END DO
      END IF

      DO N=1,NBLOCK
         IF (MGRID(N) > 1) THEN
            DO M=2,MGRID(N)
               IG1 = IG(M-1,N)
               IG2 = IG(M,N)
               CALL LUMP_DIST(DISTW(IG1),DISTW(IG2),
     &              IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),IN,JN,KN)
            END DO
         END IF
      END DO

      IF (myid == 0) THEN
         CLOSE(IGRI)
      END IF

      IF(NTASKS > 1 .AND. MYID == 0) THEN
         DEALLOCATE (B,LOCDIB)
         WRITE(45,9070)
      ENDIF
9070  FORMAT('In DISTANCES: Deallocated memory for B')

      IREPEA(2) = IREPEA(2) + 1
      RETURN
      END SUBROUTINE DISTANCES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE ONE_BLOCK(XC,YC,ZC,DISTW,XW,YW,ZW,IMAX,JMAX,KMAX,
     &     NTOT,IG,MGRID,IN,JN,KN,ICON,ICS,NPATCH,NB,MGM,NBLOCK,
     &     NBL,MG,ISTR,JSTR,KSTR,NPNUM,NG,myid,ntasks,LOCDIS,JET)

      USE MPI

      IMPLICIT NONE

      REAL    :: XC(*),YC(*),ZC(*)
      REAL    :: DISTW(*),XW(*),YW(*),ZW(*)
      INTEGER :: LOCDIS(*),JET(*)
      INTEGER :: NB,MGM,ICS,NBLOCK,IN,JN,KN,NBL,NP,MG,ISTR,JSTR,KSTR
      INTEGER :: MGRID(*),NPATCH(*),ICON(ICS,*)
      INTEGER :: IMAX(MGM,NB),JMAX(MGM,NB),KMAX(MGM,NB),NTOT(MGM,NB),
     &     IG(MGM,NB),NPNUM(*),ii

      INTEGER :: IMI,IMA,JMI,JMA,KMI,KMA,IST,JST,KST,IG1,IG2,ierr,
     &     IN1,JN1,KN1,KMX,IJSLB,ITYPE,IWALL,N,NT,M,myid,ntasks,NG,
     &     imast,I,NPNG

      IJSLB = 2
      IF (NPNUM(NG)-1 == myid) THEN
         IJSLB = 1
         DO NP = 1,NPATCH(NBL)
            ITYPE = ICON(1,NP)
            NPNG  = ICON(24,NP)
            IF (ITYPE >= 8 .AND. ITYPE <= 10) THEN
            CALL DW_PARA(ICON,ICS,IWALL,IMAX(MG,NBL),JMAX(MG,NBL),
     &           KMAX(MG,NBL),IMI,IMA,JMI,JMA,KMI,KMA,IST,JST,KST,
     &           ISTR,JSTR,KSTR,IN1,JN1,KN1,IN,JN,KN,NP)
            IG1 = IG(MG,NBL)
            CALL PICK_XYZ(IG1,XC,YC,ZC,XW(IJSLB),YW(IJSLB),ZW(IJSLB),
     &           IMI,IMA,JMI,JMA,KMI,KMA,IST,JST,KST,IN1,JN1,KN1,
     &           JET(IJSLB),NPNG)
            IJSLB = IJSLB + (IMA-IMI+1)*(JMA-JMI+1)
            ENDIF
         ENDDO
      END IF
      IJSLB = IJSLB -1

**        write(87,*) 'ijslb=',ijslb
**        do i=1,ijslb
**        write(87,*) xw(i),Yw(i),zw(i),jet(i)
**        enddo

      IF (ntasks > 1) THEN
         imast = NPNUM(NG)-1
      CALL MPI_BCAST(IJSLB,1,MPI_INTEGER,imast,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XW,IJSLB,MPI_REAL8,imast,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(YW,IJSLB,MPI_REAL8,imast,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ZW,IJSLB,MPI_REAL8,imast,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(JET,IJSLB,MPI_INTEGER,imast,MPI_COMM_WORLD,ierr)
      END IF

      DO N=1,NBLOCK
         IG1 = IG(MG,N)
         NT  = NTOT(MG,N)
**          write(6,*) 'alku=',ig1
**          write(89,*) 'alku=',ig1
         CALL WALL_DIST(XC(IG1),YC(IG1),ZC(IG1),DISTW(IG1),XW,YW,ZW,
     &        LOCDIS(IG1),IJSLB,NT,N,JET)
**         DO I = 1,NT
**         write(89,*) ig1-1+i,locdis(ig1-1+i),n,distw(ig1-1+i)
c         IF(LOCDIS(I) == N) THEN
c            LOCDIS(I) = 1
c         ELSE
c            LOCDIS(I) = 0
c         ENDIF
**         ENDDO ! Voi mersu
      END DO

      RETURN
      END SUBROUTINE ONE_BLOCK
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DW_PARA(ICON,ICS,IWALL,IMAX,JMAX,KMAX,IMI,IMA,
     &     JMI,JMA,KMI,KMA,IST,JST,KST,ISTR,JSTR,KSTR,
     &     IN1,JN1,KN1,IN,JN,KN,NP)

      IMPLICIT NONE
      INTEGER :: IMI,IMA,JMI,JMA,KMI,KMA,IST,JST,KST,IN1,JN1,KN1,IWALL
      INTEGER :: ISTR,JSTR,KSTR,IN,JN,KN,NP,ICS,IMAX,JMAX,KMAX
      INTEGER :: ICON(ICS,*)
      INTEGER :: KMX

      IWALL = ICON(3,NP)
      IMI   = ICON(4,NP)
      IMA   = ICON(5,NP)
      JMI   = ICON(6,NP)
      JMA   = ICON(7,NP)
      SELECT CASE (IWALL)
      CASE (1,4)
         KMX = IMAX
         IST = JSTR
         JST = KSTR
         KST = ISTR
         IN1 = JN
         JN1 = KN
         KN1 = IN
      CASE (2,5)
         KMX = JMAX
         IST = ISTR
         JST = KSTR
         KST = JSTR
         IN1 = IN
         JN1 = KN
         KN1 = JN
      CASE (3,6)
         KMX = KMAX
         IST = ISTR
         JST = JSTR
         KST = KSTR
         IN1 = IN
         JN1 = JN
         KN1 = KN
      END SELECT
      IF (IWALL >= 1 .AND. IWALL <= 3) THEN
         KMI = 0
         KMA = 1
      ELSE IF (IWALL >= 4 .AND. IWALL <= 6) THEN
         KMI = KMX
         KMA = KMX+1
      END IF

      RETURN
      END SUBROUTINE DW_PARA
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PICK_XYZ(IG1,XC,YC,ZC,XW,YW,ZW,
     &     IMI,IMA,JMI,JMA,KMI,KMA,ISTR,JSTR,KSTR,IN,JN,KN,JET,NBL)

      IMPLICIT NONE

      REAL    :: XC(*), YC(*), ZC(*)
      REAL    :: XW(*), YW(*), ZW(*)
      INTEGER :: IMI,IMA,JMI,JMA,KMI,KMA,ISTR,JSTR,KSTR,IN,JN,KN,IG1
      INTEGER :: ITL,ITG,ITL1,ITG1,I,J,NBL,JET(*)

      DO J = JMI,JMA
         ITL1 = (J-JMI)*(IMA-IMI+1)
         ITG1 = (KN-1)*KSTR + (J-1+JN)*JSTR + (IN-1)*ISTR + IG1
         DO I = IMI,IMA
            ITL = ITL1 + I - IMI + 1 
            ITG = ITG1 + I*ISTR

            XW(ITL) = 0.5*(XC(ITG+KMI*KSTR) + XC(ITG+KMA*KSTR))
            YW(ITL) = 0.5*(YC(ITG+KMI*KSTR) + YC(ITG+KMA*KSTR))
            ZW(ITL) = 0.5*(ZC(ITG+KMI*KSTR) + ZC(ITG+KMA*KSTR))
            JET(ITL)= NBL
         END DO
      END DO

      RETURN
      END SUBROUTINE PICK_XYZ
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE WALL_DIST(XC,YC,ZC,DISTW,XW,YW,ZW,LOCDIS,IJSLB,
     &                     NTOT,N,JET)

      IMPLICIT NONE

      REAL :: XC(*), YC(*), ZC(*)
      REAL :: XW(*), YW(*), ZW(*), DISTW(*)
      REAL :: R
      INTEGER :: IJSLB, NTOT, IJ, I, N, LOCDIS(*), JET(*)

      DO IJ=1,IJSLB
         DO I=1,NTOT
            R =  (XC(I)-XW(IJ))**2 +
     &           (YC(I)-YW(IJ))**2 +
     &           (ZC(I)-ZW(IJ))**2
             IF(R < DISTW(I)) THEN
                DISTW(I) = R
                LOCDIS(I)= JET(IJ)
             ENDIF
         END DO
      END DO

      RETURN
      END SUBROUTINE WALL_DIST
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C 
      SUBROUTINE LUMP_DIST(DISTW1,DISTW2,IMAX,JMAX,KMAX,IN,JN,KN)

C     The Wall-Distances at coarse grid levels are computed by lumping
C     the values at the previous grid level.

      REAL    :: DISTW1(*), DISTW2(*)
      REAL    :: RA
      INTEGER :: IMAX, JMAX, KMAX, IN, JN, KN

      IMH1 = (IMAX-1)/2 + 2
      JMH1 = (JMAX-1)/2 + 2
      KMH1 = (KMAX-1)/2 + 2

      ISH  = IMH1-1 + 2*IN
      JSH  = JMH1-1 + 2*JN
      KSH  = KMH1-1 + 2*KN

      IMT1 = IMAX+1
      JMT1 = JMAX+1
      KMT1 = KMAX+1

      IST  = IMAX + 2*IN
      JST  = JMAX + 2*JN
      KST  = KMAX + 2*KN

      INC = 1

C      DO K = 1,KMT1-2,2
C         DO J = 1,JMT1-2,2
C            DO I = 1,IMT1-2,2
C               IT0 = 1+((I-1)/2+IN)+((J-1)/2+JN)*ISH+
C     &              ((K-1)/2+KN)*ISH*JSH
C               IT1 = I+IN+    (J-1    +JN)*IST+(K-1    +KN)*IST*JST
C               IT2 = I+IN+    (J-1+INC+JN)*IST+(K-1    +KN)*IST*JST
C               IT3 = I+IN+    (J-1    +JN)*IST+(K-1+INC+KN)*IST*JST
C               IT4 = I+IN+    (J-1+INC+JN)*IST+(K-1+INC+KN)*IST*JST
C               IT5 = I+IN+INC+(J-1    +JN)*IST+(K-1    +KN)*IST*JST
C               IT6 = I+IN+INC+(J-1+INC+JN)*IST+(K-1    +KN)*IST*JST
C               IT7 = I+IN+INC+(J-1    +JN)*IST+(K-1+INC+KN)*IST*JST
C               IT8 = I+IN+INC+(J-1+INC+JN)*IST+(K-1+INC+KN)*IST*JST

      DO K = 1 , MAX( 1 , KMT1-2 ) , 2
         K1 = (K-1+KN)               *IST*JST
         K2 = (MIN(K+1, KMT1-1)-1+KN)*IST*JST
         DO J = 1 , MAX(1,JMT1-2) , 2
            J1 = (J-1+JN)               *IST
            J2 = (MIN(J+1, JMT1-1)-1+JN)*IST
            DO I = 1, MAX( 1 ,IMT1-2) , 2
               I1 = 1 + (I-1+IN)
               I2 = 1 + (MIN(I+1,IMT1-1)-1+IN)

               IH0 = 1+((I-1)/2+IN)+((J-1)/2+JN)*ISH+
     &              ((K-1)/2+KN)*ISH*JSH
               IT1 = I1 + J1 + K1
               IT2 = I1 + J2 + K1
               IT3 = I1 + J1 + K2
               IT4 = I1 + J2 + K2
               IT5 = I2 + J1 + K1
               IT6 = I2 + J2 + K1
               IT7 = I2 + J1 + K2
               IT8 = I2 + J2 + K2

               RA = DISTW1(IT1)+DISTW1(IT2)+DISTW1(IT3)+DISTW1(IT4)+
     &              DISTW1(IT5)+DISTW1(IT6)+DISTW1(IT7)+DISTW1(IT8)
               DISTW2(IH0) = 0.125*RA

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE LUMP_DIST
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ... New subroutines for calculating grid point distances from
C ... closest skin points. This information can be utilized in remeshing
C ... algorithms in case of deforming grids. ESa May, 2023.

      SUBROUTINE GP_DISTANCES(IGRI,WDISTS_FORCED)

C ... This subroutine computes true wall distances.

      USE MPI

      USE NS3CO, ONLY       : NBLOCK, LEVEL, NBLOCG, 
     &                        PARALLEL, IN, JN, KN

      USE INTEGERS, ONLY    : NB, MGM, IPRO, NPRO, IREPEA, NREPEA 

      USE MAIN_ARRAYS, ONLY : DISTP, IDP, IMAX, JMAX, KMAX, NTOT,
     &                        IG, MGRID, IMAXG, JMAXG, KMAXG,
     &                        NLOCAL, NPNUM
      
      IMPLICIT NONE

      INTEGER STATUS(MPI_STATUS_SIZE)

      INTEGER :: myid, ierr, IDsrc, IStag

      INTEGER :: NBL, NTOTG, NTOTMX, ERCODE, NG, MG, N, M, II, JJ,
     &           IG1, IG2, IGRI, NT, I, J, K, ISTR, JSTR, KSTR, NTMLOC

      INTEGER(KIND=8) :: NTSUM
      
      LOGICAL :: THERE, MASTER, WDISTS_FORCED
 
      REAL    :: DISTINI

      REAL,    ALLOCATABLE, DIMENSION(:) :: B
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDPB, BAPU

      CHARACTER(LEN=120) :: RIVI, RIVI2
      CHARACTER(LEN=1)   :: EXT
      CHARACTER(LEN=8)   :: FNAM


      myid     = IPRO - 1
      MASTER   = IPRO == 1
                           
      ERCODE   = 0
      MG       = 1
      NTOTMX   = 1
      DISTINI  = 1.E+8

C ... Initialize wall distance array

      DO NBL=1,NBLOCK
         IG1 = IG(MG,NBL)
         IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
         DISTP(IG1:IG2) = DISTINI
         IDP(IG1:IG2) = 0
      ENDDO

      DO NBL = 1,NBLOCG
         NTOTG = (IMAXG(NBL)+2*IN)*(JMAXG(NBL)+2*JN)*(KMAXG(NBL)+2*KN)
         NTOTMX = MAX(NTOTMX,NTOTG)
      ENDDO

      IF(MASTER) ALLOCATE(BAPU(NBLOCG), STAT=ERCODE)

C ... Allocate memory for messages

      IF(PARALLEL) THEN
         ERCODE = 0
         IF(MASTER) THEN
            ALLOCATE(B(NTOTMX), IDPB(NTOTMX), STAT=ERCODE)
         ENDIF
         CALL MPI_BCAST(ERCODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(ERCODE /= 0) THEN
            WRITE(*,*) 'GP_DISTANCES: Not enough memory. Aborting.'
            CALL MPI_FINALIZE(ierr)
         ELSE
            IF(IREPEA(2) < NREPEA(2)) WRITE(45,9050)
         ENDIF
      ENDIF

9050  FORMAT('In GP_DISTANCES: Allocated memory for messages')

C ... Check for a file DISTG_L with pre-calulated values

      EXT  = CHAR(48+LEVEL)
      FNAM = 'DISTG_L'//EXT

      IF(MASTER) INQUIRE(FILE=FNAM,EXIST=THERE)
      IF(PARALLEL) THEN
         CALL MPI_BCAST(THERE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ENDIF

      IF(MASTER) THEN
         OPEN(IGRI,FILE=FNAM,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF

C ... Read the values from DISTG_L* file (Fine grid only)
                 
      IF(IREPEA(2) <= NREPEA(2) .AND. MASTER) WRITE(*,*)

      IF(THERE .AND. .NOT.WDISTS_FORCED) THEN

         IF(MASTER) THEN
            RIVI="('  GP_DISTANCES: Reading grid point distances ',
     &             'from file ',A8,/)"
            RIVI2="('  MESSAGE: Reading grid point distances from file ',
     &      A8,' - repeated more than',I3,' times. Silence')"
            IF(IREPEA(2) < NREPEA(2)) THEN
               WRITE(*,RIVI) FNAM
            ELSE
               IF(IREPEA(2) == NREPEA(2))WRITE(*,RIVI2) FNAM,IREPEA(2)
            ENDIF
         ENDIF

         DO NG=1,NBLOCG

            IDSRC = NPNUM(NG) - 1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)

            IF(MASTER) THEN

               IF(myid == IDsrc) THEN

                  IG1 = IG(MG,N)-IG(MG,1)
                  READ(IGRI) (DISTP(IG1+I),I=1,NT)
                  READ(IGRI) (IDP(IG1+I),I=1,NT)

               ELSE

                  READ(IGRI) (B(I),I=1,NT)
                  READ(IGRI) (IDPB(I),I=1,NT)
                  CALL MPI_SEND(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)
                  CALL MPI_SEND(IDPB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)

               ENDIF

            ELSEIF(myid == IDsrc) THEN

               IG1 = IG(MG,N)-IG(MG,1)
               CALL MPI_RECV(DISTP(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
               CALL MPI_RECV(IDP(IG1+1),NT,MPI_INTEGER,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)

            ENDIF

         ENDDO

      ELSE      ! Wall distances must be calculated.


C ... Total number of distances to be calculated in this process.

         NTSUM=0

         DO N=1,NBLOCK
            NTSUM = NTSUM + NTOT(MG,N)         
         ENDDO


         IF(MASTER) THEN
            RIVI="('  GP_DISTANCES: Calculating grid point distances.')"
            WRITE(*,RIVI)
         ENDIF

         CALL DOBLOCKS(MG,NTSUM)

         DO NBL=1,NBLOCK
            IG1 = IG(MG,NBL)
            IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
         ENDDO


C ... Write the values from DISTG_L* file (Fine grid only)

         RIVI="('  GP_DISTANCES: Writing distances to file ',A8)"

         IF(MASTER) THEN
            WRITE(*,RIVI) FNAM
            WRITE(*,*)
         ENDIF

         DO NG=1,NBLOCG

            IDSRC = NPNUM(NG)-1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)

            IF(myid == IDSRC) IG1   = IG(MG,N)-IG(MG,1)

            IF(MASTER) THEN

               IF(myid == IDsrc) THEN

                  WRITE(IGRI) (DISTP(IG1+I),I=1,NT)
                  WRITE(IGRI) (IDP(IG1+I),I=1,NT)

               ELSE

                  CALL MPI_RECV(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  CALL MPI_RECV(IDPB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  WRITE(IGRI) (B(I),I=1,NT)
                  WRITE(IGRI) (IDPB(I),I=1,NT)

               ENDIF

            ELSEIF (myid == IDsrc) THEN

               CALL MPI_SEND(DISTP(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,ierr)
               CALL MPI_SEND(IDP(IG1+1),NT,MPI_INTEGER,0,IStag,
     &              MPI_COMM_WORLD,ierr)

            ENDIF

         ENDDO

      ENDIF

*      DO N=1,NBLOCK
*         IF (MGRID(N) > 1) THEN
*            DO M=2,MGRID(N)
*               IG1 = IG(M-1,N)
*               IG2 = IG(M,N)
*               CALL LUMP_DIST(DISTP(IG1),DISTP(IG2),
*     &              IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),IN,JN,KN)
*            END DO
*         END IF
*      END DO

      IF (MASTER) THEN
         DEALLOCATE(BAPU)
         CLOSE(IGRI)
      END IF

      IF(PARALLEL .AND. MASTER) THEN
         DEALLOCATE (B,IDPB)
         WRITE(45,9070)
      ENDIF

9070  FORMAT('In GP_DISTANCES: Deallocated memory for B')

      IREPEA(2) = IREPEA(2) + 1

      RETURN
      END SUBROUTINE GP_DISTANCES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE DOBLOCKS(MG,NTSUM)

C ... Calculate grid point distances within blocks one by one. First group
C ... the skin elements according to their locations in space in order to
C ... accelerate the minimum distance search. 

      USE MPI

      USE MAIN_ARRAYS, ONLY : XORI, YORI, ZORI, DISTP, IG, NTOT, IDP,
     &                        SKINX, SKINY, SKINZ,      
     &                        IMAX, JMAX, KMAX, ICON
      USE NS3CO, ONLY       : NSKIN, NBLOCK, LEVEL, PARALLEL
      USE INTEGERS, ONLY    : IPRO
      USE CONSTANTS, ONLY   : EPS6
      USE MPCCIVARS

      IMPLICIT NONE

C ... Sorting cells must be cubic, i.e., DX = DY = DZ

      INTEGER :: NBOXMAX

      INTEGER :: IG1, NT, M, N, MG, IERR, NOCCUPIED1 
      INTEGER :: NCB, IERRCODE
      INTEGER :: I, J, K, L, IT, JT, IJ, IC, NR, NBOXTOT, NBOXSCL
      INTEGER :: NXBOX1, NYBOX1, NZBOX1
      INTEGER :: NXBOX2, NYBOX2, NZBOX2
      INTEGER :: NXBOX3, NYBOX3, NZBOX3
      INTEGER(KIND=8) :: NTSUMI, NTSUM, KPROS

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OCCUPIED1, CLOSEST1
      INTEGER, ALLOCATABLE :: IBOXES1(:,:,:,:)
      INTEGER, ALLOCATABLE :: IBOXES2(:,:,:,:)
      INTEGER, ALLOCATABLE :: IBOXES3(:,:,:,:)

      INTEGER, ALLOCATABLE :: KSC1(:)  ! Skin element pointer vector
      INTEGER, ALLOCATABLE :: KSC2(:)  ! Skin element pointer vector
      INTEGER, ALLOCATABLE :: KSC3(:)  ! Skin element pointer vector

      REAL :: XMIN,  XMAX,  YMIN,  YMAX,  ZMIN,  ZMAX, TMAX
      REAL :: XMIN1, XMAX1, YMIN1, YMAX1, ZMIN1, ZMAX1
      REAL :: DX, DY, DZ, DXR, DYR, DZR, DS
 
      LOGICAL :: MASTER, TWODIM

      MASTER = IPRO == 1
      NTSUMI = 0
      KPROS  = 0
      TWODIM = KMAX(1,1) == 1


C ... Sorting box base resolution.

      NBOXMAX = 32

*      IF(NSKIN < 64000) NBOXMAX = 24
*      IF(NSKIN < 16000) NBOXMAX = 16

      IF(fnelems < 64000) NBOXMAX = 24
      IF(fnelems < 16000) NBOXMAX = 16

C ... Check the limits of the wetted boundary.

      xmin1 = +1000000.0
      xmax1 = -1000000.0
      ymin1 = +1000000.0
      ymax1 = -1000000.0
      zmin1 = +1000000.0
      zmax1 = -1000000.0
      
      do ij = 1,fnelems
         do ic = 1,4

            nr = felemnodes(ij,ic)  ! Wetted element node number      

            XMIN1 = MIN(fnodecoor_orig(nr,1),xmin1)
            XMAX1 = MAX(fnodecoor_orig(nr,1),xmax1)
            YMIN1 = MIN(fnodecoor_orig(nr,2),ymin1)
            YMAX1 = MAX(fnodecoor_orig(nr,2),ymax1)
            ZMIN1 = MIN(fnodecoor_orig(nr,3),zmin1)
            ZMAX1 = MAX(fnodecoor_orig(nr,3),zmax1)

       enddo
      enddo
      
C ... Check the limits of the whole skin boundary.

*      XMIN1 = MINVAL(SKINX)
*      XMAX1 = MAXVAL(SKINX)
*      YMIN1 = MINVAL(SKINY) 
*      YMAX1 = MAXVAL(SKINY)
*      ZMIN1 = MINVAL(SKINZ) 
*      ZMAX1 = MAXVAL(SKINZ)
      
C ... Sort the wetted elements according to their locations in space in 
C ... order to accelerate the search process.

      TMAX = MAX(XMAX1-XMIN1,YMAX1-YMIN1,ZMAX1-ZMIN1)

      DX = TMAX/NBOXMAX
      DY = DX
      DZ = DX

      XMIN = 0.5*(XMIN1+XMAX1) - (INT((XMAX1-XMIN1)/DX)+2)/2*DX - EPS6
      XMAX = 0.5*(XMIN1+XMAX1) + (INT((XMAX1-XMIN1)/DX)+2)/2*DX - EPS6
      YMIN = 0.5*(YMIN1+YMAX1) - (INT((YMAX1-YMIN1)/DY)+2)/2*DY - EPS6
      YMAX = 0.5*(YMIN1+YMAX1) + (INT((YMAX1-YMIN1)/DY)+2)/2*DY - EPS6
      ZMIN = 0.5*(ZMIN1+ZMAX1) - (INT((ZMAX1-ZMIN1)/DZ)+2)/2*DZ - EPS6
      ZMAX = 0.5*(ZMIN1+ZMAX1) + (INT((ZMAX1-ZMIN1)/DZ)+2)/2*DZ - EPS6

      NXBOX1 = NINT((XMAX-XMIN)/DX)
      NYBOX1 = NINT((YMAX-YMIN)/DY)
      NZBOX1 = NINT((ZMAX-ZMIN)/DZ)

      NBOXTOT = NXBOX1*NYBOX1*NZBOX1

      NBOXSCL = MAX((NBOXMAX**3/NBOXTOT)**(1.0/3.0),1.0)

      IF(TWODIM) THEN
         NBOXTOT = NXBOX1*NYBOX1
         NBOXSCL = MAX((NBOXMAX**2/NBOXTOT)**(1.0/2.0),1.0)
      ENDIF

      NXBOX1 = NBOXSCL * NXBOX1
      NYBOX1 = NBOXSCL * NYBOX1
      NZBOX1 = NBOXSCL * NZBOX1

      IF(TWODIM) THEN 
         NZBOX1 = 1
      ENDIF
      
      NXBOX2 = 4*NXBOX1
      NYBOX2 = 4*NYBOX1
      NZBOX2 = 4*NZBOX1

      NXBOX3 = 16*NXBOX1
      NYBOX3 = 16*NYBOX1
      NZBOX3 = 16*NZBOX1

C ... Sort the wetted elements differently for each block based on
C ... different overlapping grid priorities. This increases computing       
C ... time but needs to be done only once in the beginning of the 
C ... simulation.
      
      DO N=1,NBLOCK

         IG1 = IG(MG,N)
         NT  = NTOT(MG,N)

         ALLOCATE(IBOXES1(NXBOX1,NYBOX1,NZBOX1,2),STAT=IERRCODE)
         ALLOCATE(IBOXES2(NXBOX2,NYBOX2,NZBOX2,2),STAT=IERRCODE)
         ALLOCATE(IBOXES3(NXBOX3,NYBOX3,NZBOX3,2),STAT=IERRCODE)

         DX  = (XMAX-XMIN)/(NXBOX1)
         DY  = (YMAX-YMIN)/(NYBOX1)
         DZ  = (ZMAX-ZMIN)/(NZBOX1)

         DXR = 1.0/DX
         DYR = 1.0/DY
         DZR = 1.0/DZ

***** Coarse skin element sorting **************************************

         ALLOCATE(KSC1(1))

         CALL BOXEL_WETTED(N,KSC1,NCB,IBOXES1,NXBOX1,NYBOX1,NZBOX1,
     &                     XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)

         DEALLOCATE(KSC1)
         ALLOCATE(KSC1(NCB),STAT=IERRCODE)

         IF(IERRCODE > 0) THEN                                          
            WRITE(*,*)'DOBLOCKS:  Not enough memory, aborting ...'             
            STOP                                                           
         ENDIF                                                    

         CALL BOXEL_WETTED(N,KSC1,NCB,IBOXES1,NXBOX1,NYBOX1,NZBOX1,
     &                     XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)

***** 4 x finer element sorting ****************************************

         DS = 4.0

         ALLOCATE(KSC2(1))

         CALL BOXEL_WETTED(N,KSC2,NCB,IBOXES2,NXBOX2,NYBOX2,NZBOX2,
     &                     XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,1)

         DEALLOCATE(KSC2)
         ALLOCATE(KSC2(NCB),STAT=IERRCODE)

         IF(IERRCODE > 0) THEN                                          
            WRITE(*,*)'DOBLOCKS:  Not enough memory, aborting ...'             
            STOP                                                           
         ENDIF                                                    

         CALL BOXEL_WETTED(N,KSC2,NCB,IBOXES2,NXBOX2,NYBOX2,NZBOX2,
     &                     XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,2)

***** 16 x finer element sorting ***************************************

         DS = 16.0

         ALLOCATE(KSC3(1))

         CALL BOXEL_WETTED(N,KSC3,NCB,IBOXES3,NXBOX3,NYBOX3,NZBOX3,
     &                     XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,1)

         DEALLOCATE(KSC3)
         ALLOCATE(KSC3(NCB),STAT=IERRCODE)

         IF(IERRCODE > 0) THEN                                          
            WRITE(*,*)'DOBLOCKS:  Not enough memory, aborting ...'             
            STOP                                                           
         ENDIF                                                    

         CALL BOXEL_WETTED(N,KSC3,NCB,IBOXES3,NXBOX3,NYBOX3,NZBOX3,
     &                     XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,2)

      
C ... Now the skin elements are grouped according to their locations in space. 
C ... These groups are utilized in the closest skin element search.

C ... Occupied voxels in the coarse sorting grid. 

         NOCCUPIED1 = 0
 
         DO K = 1,NZBOX1
            DO J = 1,NYBOX1
               DO I = 1,NXBOX1
               
                  IF(IBOXES1(I,J,K,1) == 0) CYCLE
               
                  NOCCUPIED1 = NOCCUPIED1 + 1
               
               ENDDO
            ENDDO
         ENDDO

         ALLOCATE(OCCUPIED1(3,NOCCUPIED1),STAT=IERRCODE)
         ALLOCATE(CLOSEST1(3,NOCCUPIED1),STAT=IERRCODE)
      
         M = 0
      
         DO K = 1,NZBOX1
            DO J = 1,NYBOX1
               DO I = 1,NXBOX1
               
                  IF(IBOXES1(I,J,K,1) > 0) THEN
                     M = M + 1
                     OCCUPIED1(1,M) = I
                     OCCUPIED1(2,M) = J
                     OCCUPIED1(3,M) = K
                  ENDIF
               
               ENDDO
            ENDDO
         ENDDO

C ... Block loop runs here

*      DO N=1,NBLOCK
*
*         IG1 = IG(MG,N)
*         NT  = NTOT(MG,N)

         CALL GPDISTS(N,XORI(IG1),YORI(IG1),ZORI(IG1),NT,DISTP(IG1),
     &        IDP(IG1),NTSUMI,NTSUM,KPROS,TWODIM,
     &        IBOXES1,NXBOX1,NYBOX1,NZBOX1,KSC1,
     &        IBOXES2,NXBOX2,NYBOX2,NZBOX2,KSC2,
     &        IBOXES3,NXBOX3,NYBOX3,NZBOX3,KSC3,fnelems,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,DXR,DYR,DZR,
     &        XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX)

*      ENDDO

         DEALLOCATE(OCCUPIED1,CLOSEST1)
         DEALLOCATE(KSC1,IBOXES1,KSC2,IBOXES2,KSC3,IBOXES3)

      ENDDO

      RETURN      
      END SUBROUTINE DOBLOCKS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE GPDISTS(N,XORI,YORI,ZORI,NT,DISTP,
     &        IDP,NTSUMI,NTSUM,KPROS,TWODIM,
     &        IBOXES1,NXBOX1,NYBOX1,NZBOX1,KSC1,
     &        IBOXES2,NXBOX2,NYBOX2,NZBOX2,KSC2,
     &        IBOXES3,NXBOX3,NYBOX3,NZBOX3,KSC3,NSKIN,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,DXR,DYR,DZR,
     &        XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX)


C ... Calculate grid point distances within one block.

      USE MPCCIVARS, ONLY   : fnelems, felemnodes, fnodecoor_orig,
     &                        fnodechimt

      USE MAIN_ARRAYS, ONLY : SKINX, SKINY, SKINZ, ITSKIN, 
     &                        NPROCE, NCHIMT, IMAX, JMAX, KMAX    
      USE INTEGERS, ONLY    : IPRO 
      USE NS3CO, ONLY       : PARALLEL, TIMEL
      USE MPI

      IMPLICIT NONE

      INTEGER :: NT, NSKIN, I, J, K, L, N, NG, IJ, CHT
      INTEGER :: NOCCUPIED, NCLOSEST
      INTEGER :: NOCCUPIED1, NCLOSEST1
      INTEGER :: II, JJ, KK, IP, JP, KP
      INTEGER :: IERRCODE, ISIDE, ISIDEL, ITYPEL
      INTEGER :: LBS, LBE, IC, LC 
      INTEGER :: IL, IU, JL, JU, KL, KU

      INTEGER :: NXBOX1,  NYBOX1,  NZBOX1
      INTEGER :: NXBOX2,  NYBOX2,  NZBOX2
      INTEGER :: NXBOX3,  NYBOX3,  NZBOX3
      INTEGER(KIND=8) :: NTSUMI, NTSUM, JPROS, KPROS
 
      INTEGER, DIMENSION(NXBOX1,NYBOX1,NZBOX1,2) :: IBOXES1
      INTEGER, DIMENSION(NXBOX2,NYBOX2,NZBOX2,2) :: IBOXES2
      INTEGER, DIMENSION(NXBOX3,NYBOX3,NZBOX3,2) :: IBOXES3
      INTEGER, DIMENSION(3,NOCCUPIED1) :: OCCUPIED1, CLOSEST1
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OCCUPIED, CLOSEST
      INTEGER, DIMENSION(*) :: KSC1, KSC2, KSC3, IDP

      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
      REAL :: DXR, DYR, DZR, DS
      REAL :: DISTP(*)

      REAL, DIMENSION(3) :: GP, SKINP
      REAL, DIMENSION(*) :: XORI, YORI, ZORI
      REAL :: DD, DISTWL

      LOGICAL :: MASTER, OVERLA, TWODIM

      LOGICAL, ALLOCATABLE, DIMENSION(:) :: DONE

      
      MASTER = IPRO == 1

      NG = NPROCE(N+1,IPRO)  ! Global block number

      ALLOCATE(DONE(NSKIN),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'COMPDISTS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    
      
C ... Grid point loop starts here

      DO IP=1,NT
        
         JPROS = 100-100*(NTSUM-NTSUMI)/NTSUM

         IF(JPROS >= KPROS) THEN
            IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERRCODE)
            IF(KPROS == 0 .AND. MASTER) write(*,*)
            IF(MASTER) WRITE(*,'(I5," % completed")') KPROS
            IF(KPROS == 100 .AND. MASTER) write(*,*)
            KPROS = KPROS + 5
         ENDIF

         NTSUMI = NTSUMI + 1
         
         DONE(1:NSKIN) = .FALSE.

         DISTWL = DISTP(IP)

         GP(1) = XORI(IP); GP(2) = YORI(IP); GP(3) = ZORI(IP)
         
C ... Closest voxels search on the start level.

         CALL CLOSEST_VOXELS(GP,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &        NXBOX1,NYBOX1,NZBOX1,IBOXES1,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,NCLOSEST1)

C ... Closest voxels search on the 4 x finer level.

         DS = 4.0

         ALLOCATE(OCCUPIED(3,64*NOCCUPIED1),STAT=IERRCODE)

         NOCCUPIED = 0
 
         DO L=1,NCLOSEST1

            I = CLOSEST1(1,L)
            J = CLOSEST1(2,L)
            K = CLOSEST1(3,L)

            IU = 4*I
            IL = IU-3
            JU = 4*J
            JL = JU-3
            KU = 4*K
            KL = KU-3

            DO KK=KL,KU
               DO JJ=JL,JU
                  DO II=IL,IU
               
                     IF(IBOXES2(II,JJ,KK,1) == 0) CYCLE

                     NOCCUPIED = NOCCUPIED + 1

                     OCCUPIED(1,NOCCUPIED) = II
                     OCCUPIED(2,NOCCUPIED) = JJ
                     OCCUPIED(3,NOCCUPIED) = KK

                  ENDDO
               ENDDO
            ENDDO

         ENDDO

         ALLOCATE(CLOSEST(3,NOCCUPIED),STAT=IERRCODE)

         CALL CLOSEST_VOXELS(GP,XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,
     &        NXBOX2,NYBOX2,NZBOX2,IBOXES2,
     &        OCCUPIED,NOCCUPIED,CLOSEST,NCLOSEST)

         DEALLOCATE(OCCUPIED)


C ... Closest voxels search on the 16 x finer level.

         DS = 16.0

         ALLOCATE(OCCUPIED(3,64*NOCCUPIED),STAT=IERRCODE)

         NOCCUPIED = 0
 
         DO L=1,NCLOSEST

            I = CLOSEST(1,L)
            J = CLOSEST(2,L)
            K = CLOSEST(3,L)

            IU = 4*I
            IL = IU-3
            JU = 4*J
            JL = JU-3
            KU = 4*K
            KL = KU-3

            DO KK=KL,KU
               DO JJ=JL,JU
                  DO II=IL,IU
               
                     IF(IBOXES3(II,JJ,KK,1) == 0) CYCLE

                     NOCCUPIED = NOCCUPIED + 1

                     OCCUPIED(1,NOCCUPIED) = II
                     OCCUPIED(2,NOCCUPIED) = JJ
                     OCCUPIED(3,NOCCUPIED) = KK

                  ENDDO
               ENDDO
            ENDDO

         ENDDO

         DEALLOCATE(CLOSEST)         
         ALLOCATE(CLOSEST(3,NOCCUPIED),STAT=IERRCODE)

         CALL CLOSEST_VOXELS(GP,XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,
     &        NXBOX3,NYBOX3,NZBOX3,IBOXES3,
     &        OCCUPIED,NOCCUPIED,CLOSEST,NCLOSEST)

         DEALLOCATE(OCCUPIED)

         
C ... Scan through the skin elements in the closest sorting cells.

         DO J=1,NCLOSEST
            
            II = CLOSEST(1,J)
            JJ = CLOSEST(2,J)
            KK = CLOSEST(3,J)

            LBS = IBOXES3(II,JJ,KK,2)            ! Start
            LBE = LBS + IBOXES3(II,JJ,KK,1) - 1  ! End

            DO LC = LBS,LBE
          
               IJ = KSC3(LC)                     ! Wetted element number

               IF(DONE(IJ)) CYCLE
               
               DO IC=1,4  ! Wetted element corner point loop

                  SKINP(1) = fnodecoor_orig(felemnodes(IJ,IC),1)
                  SKINP(2) = fnodecoor_orig(felemnodes(IJ,IC),2) 
                  SKINP(3) = fnodecoor_orig(felemnodes(IJ,IC),3)
                  CHT      = fnodechimt(felemnodes(IJ,IC))
                  
                  DD = SQRT(SUM((GP-SKINP)**2))
              
                  IF(DD < DISTWL) THEN
                     DISTWL    = DD
                     DISTP(IP) = DISTWL
                     IDP(IP)   = felemnodes(IJ,IC)
                  ENDIF 

               ENDDO  ! IC=1,4 Wetted element corner point loop

               DONE(IJ) = .TRUE.

            ENDDO

         ENDDO             

         DEALLOCATE(CLOSEST)

      END DO  ! IP=1,NT grid point loop

      DEALLOCATE(DONE)

      RETURN
      END SUBROUTINE GPDISTS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

C ... New wall distance subroutines (TRUE DISTANCES) ESa Nov, 2008

      SUBROUTINE TRUE_DISTANCES(IGRI,WDISTS_FORCED)

C ... This subroutine computes true wall distances for cell center points.

      USE MPI

      USE NS3CO, ONLY       : NBLOCK, LEVEL, NBLOCG, MPCCIL, MODALFSIL, 
     &                        NSKIN, TIMEL, PARALLEL,
     &                        IN, JN, KN

      USE INTEGERS, ONLY    : NB, MGM, IPRO, NPRO, IREPEA, NREPEA 

      USE MAIN_ARRAYS, ONLY : DISTW, IMAX, JMAX, KMAX, NTOT,
     &                        IG, MGRID, IMAXG, JMAXG, KMAXG,
     &                        NLOCAL, NPNUM, LOCDIS, BLANK, NEARBLOCK    
      
      IMPLICIT NONE

      INTEGER STATUS(MPI_STATUS_SIZE)

      INTEGER :: myid, ierr, IDsrc, IStag

      INTEGER :: NBL, NTOTG, NTOTMX, ERCODE, NG, MG, N, M, II, JJ,
     &           IG1, IG2, IGRI, NT, I, J, K, ISTR, JSTR, KSTR, NTMLOC

      INTEGER(KIND=8) :: NTSUM, NTSUMG
      
      LOGICAL :: THERE, MASTER, ITSME, WDISTS_FORCED

      REAL    :: DISTINI

      REAL,    ALLOCATABLE, DIMENSION(:) :: B, BLANKB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCDIB, BAPU

      INTEGER(KIND=8), DIMENSION(2) :: NTMAX_LOCAL, NTMAX_GLOBAL

      CHARACTER(LEN=120) :: RIVI, RIVI2
      CHARACTER(LEN=1)   :: EXT
      CHARACTER(LEN=8)   :: FNAM


      myid     = IPRO - 1
      MASTER   = IPRO == 1
                           
      ERCODE   = 0
      MG       = 1
      NTOTMX   = 1
      DISTINI  = 1.E+8

C ... Initialize wall distance array

      DO NBL=1,NBLOCK
         IG1 = IG(MG,NBL)
         IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
         DISTW(IG1:IG2) = DISTINI
         LOCDIS(IG1:IG2)= 0
      END DO

      DO NBL = 1,NBLOCG
         NTOTG = (IMAXG(NBL)+2*IN)*(JMAXG(NBL)+2*JN)*(KMAXG(NBL)+2*KN)
         NTOTMX = MAX(NTOTMX,NTOTG)
      ENDDO

      IF(MASTER) ALLOCATE(BAPU(NBLOCG), STAT=ERCODE)

C ... Allocate memory for messages

      IF (PARALLEL) THEN
         ERCODE = 0
         IF (MASTER) THEN
            ALLOCATE(B(NTOTMX), LOCDIB(NTOTMX),
     &               BLANKB(NTOTMX), STAT=ERCODE)
         END IF
         CALL MPI_BCAST(ERCODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF (ERCODE /= 0) THEN
            WRITE(*,*) 'DISTANCES: Not enough memory. Aborting.'
            CALL MPI_FINALIZE(ierr)
         ELSE
            IF(IREPEA(2) < NREPEA(2)) WRITE(45,9050)
         END IF
      END IF

9050  FORMAT('In TRUE_DISTANCES: Allocated memory for messages')

C ... Check for a file DISTW_L with pre-calulated values

      EXT  = CHAR(48+LEVEL)
      FNAM = 'DISTW_L'//EXT

      IF (MASTER) INQUIRE(FILE=FNAM,EXIST=THERE)
      IF (PARALLEL) THEN
         CALL MPI_BCAST(THERE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      END IF

      IF (MASTER) THEN
         OPEN(IGRI,FILE=FNAM,STATUS='UNKNOWN',FORM='UNFORMATTED')
      END IF

C ... Read the values from DISTW_L* file (Fine grid only)
                 
      IF (IREPEA(2) <= NREPEA(2) .AND. MASTER) WRITE(*,*)

      IF (THERE .AND. .NOT.WDISTS_FORCED) THEN

         IF (MASTER) THEN
            RIVI="('  TRUE_DISTANCES: Reading wall distances ',
     &             'from file ',A8)"
            RIVI2="('  MESSAGE: Reading wall distances from file ',
     &      A8,' - repeated more than',I3,' times. Silence')"
            IF(IREPEA(2) < NREPEA(2)) THEN
               WRITE(*,RIVI) FNAM
            ELSE
               IF(IREPEA(2) == NREPEA(2))WRITE(*,RIVI2) FNAM,IREPEA(2)
            ENDIF
         END IF
         DO NG=1,NBLOCG
            IDSRC = NPNUM(NG) - 1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
            IF (MASTER) THEN
               IF (myid == IDsrc) THEN
                  IG1   = IG(MG,N)-IG(MG,1)
                  READ(IGRI) (DISTW(IG1+I),I=1,NT)
                  READ(IGRI) (LOCDIS(IG1+I),I=1,NT)
                  READ(IGRI) (BLANK(IG1+I),I=1,NT)

                  DO I=1,NT
                     NEARBLOCK((NG-1)*NBLOCG + LOCDIS(IG1+I)) = 1
                  ENDDO

               ELSE
                  READ(IGRI) (B(I),I=1,NT)
                  READ(IGRI) (LOCDIB(I),I=1,NT)
                  READ(IGRI) (BLANKB(I),I=1,NT)
                  CALL MPI_SEND(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)
                  CALL MPI_SEND(LOCDIB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)
                  CALL MPI_SEND(BLANKB,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,ierr)

                  DO I=1,NT
                     NEARBLOCK((NG-1)*NBLOCG + LOCDIB(I)) = 1
                  ENDDO

               END IF
            ELSEIF (myid == IDsrc) THEN
               IG1   = IG(MG,N)-IG(MG,1)
               CALL MPI_RECV(DISTW(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
               CALL MPI_RECV(LOCDIS(IG1+1),NT,MPI_INTEGER,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
               CALL MPI_RECV(BLANK(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,STATUS,ierr)
            END IF
         END DO


C ... A clumsy way to identify potential Chimera block pairs based on 
C ... common closest wall.

         IF(MASTER) THEN
            DO I=1,NBLOCG
               K = 0
               DO J=1,NBLOCG
                  IF(NEARBLOCK((J-1)*NBLOCG+I) == 1) THEN
                     K = K + 1
                     BAPU(K) = J
                  ENDIF
               ENDDO
               DO II = 1,K
                  DO JJ = II,K
                     IF(NEARBLOCK((BAPU(II)-1)*NBLOCG + BAPU(JJ)) == 0)
     &                  NEARBLOCK((BAPU(II)-1)*NBLOCG + BAPU(JJ)) = 2
                     IF(NEARBLOCK((BAPU(JJ)-1)*NBLOCG + BAPU(II)) == 0)
     &                  NEARBLOCK((BAPU(JJ)-1)*NBLOCG + BAPU(II)) = 2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF (PARALLEL) THEN
            CALL MPI_BCAST(NEARBLOCK,NBLOCG*NBLOCG,MPI_INTEGER,0,
     &                     MPI_COMM_WORLD,IERR)
         END IF

      ELSE      ! Wall distances must be calculated.


C ... Total number of distances to be calculated in this process.

         NTSUM=0
         DO N=1,NBLOCK
            NTSUM = NTSUM + NTOT(MG,N)         
         ENDDO

C ... Which proces has the maximum task to do? Actually, this information
C ... is not needed any more. This is just a reminder how to do it, if it's
C ... needed one day.

         IF(PARALLEL) THEN
            NTMAX_LOCAL(1) = NTSUM
            NTMAX_LOCAL(2) = IPRO-1
            CALL MPI_REDUCE(NTMAX_LOCAL,NTMAX_GLOBAL,1,MPI_2INTEGER,
     &                      MPI_MAXLOC,0,MPI_COMM_WORLD,IERR)
            IF(MASTER) THEN
               NTSUMG = NTMAX_GLOBAL(1) ! max value across all processes
               NTMLOC = NTMAX_GLOBAL(2) ! process owning max value
            ENDIF
            CALL MPI_BCAST(NTMLOC,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(NTSUMG,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
            ITSME = NTMLOC == IPRO-1 
         ELSE
            NTSUMG = NTSUM
            ITSME  = .TRUE.
         ENDIF


         IF (MASTER) THEN
         RIVI="('  TRUE_DISTANCES: Calculating true wall distances.')"
         WRITE(*,RIVI)
         END IF

         CALL LOOPBLOCKS(MG,ITSME,NTSUM,NTSUMG)

         DO NBL=1,NBLOCK
            IG1 = IG(MG,NBL)
            IG2 = IG(MG,NBL) + NTOT(MG,NBL) - 1
         END DO


C ... Write the values from DISTW_L* file (Fine grid only)

         IF(.NOT.MPCCIL .AND. .NOT.MODALFSIL) THEN
         
         RIVI="('  TRUE_DISTANCES: Writing wall distances to file ',A8)"
         IF (MASTER) THEN
            WRITE(*,RIVI) FNAM
            WRITE(*,*)
         END IF

         DO NG=1,NBLOCG
            IDSRC = NPNUM(NG)-1
            ISTAG = 500 + NG
            N     = NLOCAL(NG)
            NT    = (IMAXG(NG)+2*IN)*(JMAXG(NG)+2*JN)*(KMAXG(NG)+2*KN)
            IF (myid == IDSRC) IG1   = IG(MG,N)-IG(MG,1)
            IF (MASTER) THEN
               IF (myid == IDsrc) THEN
                  WRITE(IGRI) (DISTW(IG1+I),I=1,NT)
                  WRITE(IGRI) (LOCDIS(IG1+I),I=1,NT)
                  WRITE(IGRI) (BLANK(IG1+I),I=1,NT)

                  DO I=1,NT
                     NEARBLOCK((NG-1)*NBLOCG + LOCDIS(IG1+I)) = 1
                  ENDDO

               ELSE
                  CALL MPI_RECV(B,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  CALL MPI_RECV(LOCDIB,NT,MPI_INTEGER,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  CALL MPI_RECV(BLANKB,NT,MPI_REAL8,IDsrc,IStag,
     &                 MPI_COMM_WORLD,STATUS,ierr)
                  WRITE(IGRI) (B(I),I=1,NT)
                  WRITE(IGRI) (LOCDIB(I),I=1,NT)
                  WRITE(IGRI) (BLANKB(I),I=1,NT)

                  DO I=1,NT
                     NEARBLOCK((NG-1)*NBLOCG + LOCDIB(I)) = 1
                  ENDDO

               END IF
            ELSEIF (myid == IDsrc) THEN
               CALL MPI_SEND(DISTW(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,ierr)
               CALL MPI_SEND(LOCDIS(IG1+1),NT,MPI_INTEGER,0,IStag,
     &              MPI_COMM_WORLD,ierr)
               CALL MPI_SEND(BLANK(IG1+1),NT,MPI_REAL8,0,IStag,
     &              MPI_COMM_WORLD,ierr)
            END IF
         END DO

         ENDIF
         
C ... A clumsy way to identify potential Chimera block pairs based on 
C ... common closest wall.

         IF(MASTER) THEN
            DO I=1,NBLOCG
               K = 0
               DO J=1,NBLOCG
                  IF(NEARBLOCK((J-1)*NBLOCG+I) == 1) THEN
                     K = K + 1
                     BAPU(K) = J
                  ENDIF
               ENDDO
               DO II = 1,K
                  DO JJ = II,K
                     IF(NEARBLOCK((BAPU(II)-1)*NBLOCG + BAPU(JJ)) == 0)
     &                  NEARBLOCK((BAPU(II)-1)*NBLOCG + BAPU(JJ)) = 2
                     IF(NEARBLOCK((BAPU(JJ)-1)*NBLOCG + BAPU(II)) == 0)
     &                  NEARBLOCK((BAPU(JJ)-1)*NBLOCG + BAPU(II)) = 2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF (PARALLEL) THEN
            CALL MPI_BCAST(NEARBLOCK,NBLOCG*NBLOCG,MPI_INTEGER,0,
     &                     MPI_COMM_WORLD,IERR)
         END IF

      END IF

      DO N=1,NBLOCK
         IF (MGRID(N) > 1) THEN
            DO M=2,MGRID(N)
               IG1 = IG(M-1,N)
               IG2 = IG(M,N)
               CALL LUMP_DIST(DISTW(IG1),DISTW(IG2),
     &              IMAX(M-1,N),JMAX(M-1,N),KMAX(M-1,N),IN,JN,KN)
            END DO
         END IF
      END DO

      IF (MASTER) THEN
         DEALLOCATE(BAPU)
         CLOSE(IGRI)
      END IF

      IF(PARALLEL .AND. MASTER) THEN
         DEALLOCATE (B,LOCDIB)
         WRITE(45,9070)
      ENDIF

9070  FORMAT('In TRUE_DISTANCES: Deallocated memory for B')

      IREPEA(2) = IREPEA(2) + 1

      RETURN
      END SUBROUTINE TRUE_DISTANCES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE LOOPBLOCKS(MG,ITSME,NTSUM,NTSUMG)

C ... Calculate true wall distances within blocks one by one. First group
C ... the skin elements according to their locations in space in order to
C ... accelerate the minimum distance search. 

      USE MPI

      USE MAIN_ARRAYS, ONLY : XC, YC, ZC, DISTW, IG, NTOT, LOCDIS,
     &                        SKINX, SKINY, SKINZ, IBSKIN,     
     &                        IMAX, JMAX, KMAX, BLANK, ICON
      USE NS3CO, ONLY       : NSKIN, NBLOCK, LEVEL, PARALLEL
      USE INTEGERS, ONLY    : IPRO
      USE CONSTANTS, ONLY   : EPS6

      IMPLICIT NONE

C ... Sorting cells must be cubic, i.e., DX = DY = DZ

      INTEGER :: NBOXMAX

      INTEGER :: IG1, NT, M, N, MG, IERR, NOCCUPIED1 
      INTEGER :: NCB, IERRCODE
      INTEGER :: I, J, K, L, IT, JT, NBOXTOT, NBOXSCL
      INTEGER :: NXBOX1, NYBOX1, NZBOX1
      INTEGER :: NXBOX2, NYBOX2, NZBOX2
      INTEGER :: NXBOX3, NYBOX3, NZBOX3
      INTEGER(KIND=8) :: NTSUMI, NTSUM, NTSUMG, KPROS

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OCCUPIED1, CLOSEST1
      INTEGER, ALLOCATABLE :: IBOXES1(:,:,:,:)
      INTEGER, ALLOCATABLE :: IBOXES2(:,:,:,:)
      INTEGER, ALLOCATABLE :: IBOXES3(:,:,:,:)

      INTEGER, ALLOCATABLE :: KSC1(:)  ! Skin element pointer vector
      INTEGER, ALLOCATABLE :: KSC2(:)  ! Skin element pointer vector
      INTEGER, ALLOCATABLE :: KSC3(:)  ! Skin element pointer vector

      REAL, ALLOCATABLE :: SKINXCP(:),SKINYCP(:),SKINZCP(:) 
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TRNX, TRNY, TRNZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TRCX, TRCY, TRCZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TRR

      REAL, DIMENSION(3) :: S1, S2, S3, S4, TC 

      REAL :: XMIN,  XMAX,  YMIN,  YMAX,  ZMIN,  ZMAX, TMAX
      REAL :: XMIN1, XMAX1, YMIN1, YMAX1, ZMIN1, ZMAX1
      REAL :: DX, DY, DZ, DXR, DYR, DZR, DS
 
      LOGICAL :: MASTER, ITSME, TWODIM

      MASTER = IPRO == 1
      NTSUMI = 0
      KPROS  = 0
      TWODIM = KMAX(1,1) == 1


C ... Sorting box base resolution.

      NBOXMAX = 32
      IF(NSKIN < 64000) NBOXMAX = 24
      IF(NSKIN < 16000) NBOXMAX = 16
      
C ... Check the limits of the whole skin boundary.

      XMIN1 = MINVAL(SKINX)
      XMAX1 = MAXVAL(SKINX)
      YMIN1 = MINVAL(SKINY) 
      YMAX1 = MAXVAL(SKINY)
      ZMIN1 = MINVAL(SKINZ) 
      ZMAX1 = MAXVAL(SKINZ)

C ... Sort the skin elements according to their locations in space in 
C ... order to accelerate the search process.

      TMAX = MAX(XMAX1-XMIN1,YMAX1-YMIN1,ZMAX1-ZMIN1)

      DX = TMAX/NBOXMAX
      DY = DX
      DZ = DX

      XMIN = 0.5*(XMIN1+XMAX1) - (INT((XMAX1-XMIN1)/DX)+2)/2*DX - EPS6
      XMAX = 0.5*(XMIN1+XMAX1) + (INT((XMAX1-XMIN1)/DX)+2)/2*DX - EPS6
      YMIN = 0.5*(YMIN1+YMAX1) - (INT((YMAX1-YMIN1)/DY)+2)/2*DY - EPS6
      YMAX = 0.5*(YMIN1+YMAX1) + (INT((YMAX1-YMIN1)/DY)+2)/2*DY - EPS6
      ZMIN = 0.5*(ZMIN1+ZMAX1) - (INT((ZMAX1-ZMIN1)/DZ)+2)/2*DZ - EPS6
      ZMAX = 0.5*(ZMIN1+ZMAX1) + (INT((ZMAX1-ZMIN1)/DZ)+2)/2*DZ - EPS6

      NXBOX1 = NINT((XMAX-XMIN)/DX)
      NYBOX1 = NINT((YMAX-YMIN)/DY)
      NZBOX1 = NINT((ZMAX-ZMIN)/DZ)

      NBOXTOT = NXBOX1*NYBOX1*NZBOX1
      NBOXSCL = MAX((NBOXMAX**3/NBOXTOT)**(1.0/3.0),1.0)

      IF(TWODIM) THEN
         NBOXTOT = NXBOX1*NYBOX1
         NBOXSCL = MAX((NBOXMAX**2/NBOXTOT)**(1.0/2.0),1.0)
      ENDIF

      NXBOX1 = NBOXSCL * NXBOX1
      NYBOX1 = NBOXSCL * NYBOX1
      NZBOX1 = NBOXSCL * NZBOX1

      IF(TWODIM) THEN 
         NZBOX1 = 1
      ENDIF
      
      NXBOX2 = 4*NXBOX1
      NYBOX2 = 4*NYBOX1
      NZBOX2 = 4*NZBOX1

      NXBOX3 = 16*NXBOX1
      NYBOX3 = 16*NYBOX1
      NZBOX3 = 16*NZBOX1

      ALLOCATE(IBOXES1(NXBOX1,NYBOX1,NZBOX1,2),STAT=IERRCODE)
      ALLOCATE(IBOXES2(NXBOX2,NYBOX2,NZBOX2,2),STAT=IERRCODE)
      ALLOCATE(IBOXES3(NXBOX3,NYBOX3,NZBOX3,2),STAT=IERRCODE)

      DX  = (XMAX-XMIN)/(NXBOX1)
      DY  = (YMAX-YMIN)/(NYBOX1)
      DZ  = (ZMAX-ZMIN)/(NZBOX1)

      DXR = 1.0/DX
      DYR = 1.0/DY
      DZR = 1.0/DZ

***** Coarse skin element sorting **************************************

      ALLOCATE(KSC1(1))

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC1,NCB,
     &                 IBOXES1,NXBOX1,NYBOX1,NZBOX1,
     &                 XMIN,YMIN,ZMIN,DXR,DYR,DZR,1)

      DEALLOCATE(KSC1)
      ALLOCATE(KSC1(NCB),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'LOOPBLOCKS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC1,NCB,
     &                 IBOXES1,NXBOX1,NYBOX1,NZBOX1,
     &                 XMIN,YMIN,ZMIN,DXR,DYR,DZR,2)

***** 4 x finer element sorting ****************************************

      DS = 4.0

      ALLOCATE(KSC2(1))

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC2,NCB,
     &                 IBOXES2,NXBOX2,NYBOX2,NZBOX2,
     &                 XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,1)

      DEALLOCATE(KSC2)
      ALLOCATE(KSC2(NCB),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'LOOPBLOCKS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC2,NCB,
     &                 IBOXES2,NXBOX2,NYBOX2,NZBOX2,
     &                 XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,2)

***** 16 x finer element sorting ***************************************

      DS = 16.0

      ALLOCATE(KSC3(1))

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC3,NCB,
     &                 IBOXES3,NXBOX3,NYBOX3,NZBOX3,
     &                 XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,1)

      DEALLOCATE(KSC3)
      ALLOCATE(KSC3(NCB),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'LOOPBLOCKS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      CALL BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC3,NCB,
     &                 IBOXES3,NXBOX3,NYBOX3,NZBOX3,
     &                 XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,2)

************************************************************************
      
C ... Now the skin elements are grouped according to their locations in space. 
C ... These groups are utilized in the closest skin element search.


C ... Skin element center points.

      ALLOCATE(SKINXCP(NSKIN),STAT=IERRCODE)
      ALLOCATE(SKINYCP(NSKIN),STAT=IERRCODE)
      ALLOCATE(SKINZCP(NSKIN),STAT=IERRCODE)

C ... Triangle normals (trapezoids divided into four triangles).
      
      ALLOCATE(TRNX(4,NSKIN),STAT=IERRCODE)
      ALLOCATE(TRNY(4,NSKIN),STAT=IERRCODE)
      ALLOCATE(TRNZ(4,NSKIN),STAT=IERRCODE)

C ... Triangle center points.
      
      ALLOCATE(TRCX(5,NSKIN),STAT=IERRCODE)
      ALLOCATE(TRCY(5,NSKIN),STAT=IERRCODE)
      ALLOCATE(TRCZ(5,NSKIN),STAT=IERRCODE)

C ... Triangle radiuses.
      
      ALLOCATE(TRR(5,NSKIN),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'CALCDISTS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    

      
      DO L = 1,NSKIN

C ... Element center points.
         
      SKINXCP(L) = 0.25*(SKINX(1,L)+SKINX(2,L)+SKINX(3,L)+SKINX(4,L))
      SKINYCP(L) = 0.25*(SKINY(1,L)+SKINY(2,L)+SKINY(3,L)+SKINY(4,L))
      SKINZCP(L) = 0.25*(SKINZ(1,L)+SKINZ(2,L)+SKINZ(3,L)+SKINZ(4,L))

      
C ... Surface element center points and radiuses
         
      S1(1) = SKINX(1,L)
      S1(2) = SKINY(1,L)
      S1(3) = SKINZ(1,L)
         
      S2(1) = SKINX(2,L)
      S2(2) = SKINY(2,L)
      S2(3) = SKINZ(2,L)
         
      S3(1) = SKINX(3,L)
      S3(2) = SKINY(3,L)
      S3(3) = SKINZ(3,L)
         
      S4(1) = SKINX(4,L)
      S4(2) = SKINY(4,L)
      S4(3) = SKINZ(4,L)
         
      TC = (S1 + S2 + S3 + S4) / 4.0  
         
      TRCX(5,L) = TC(1) 
      TRCY(5,L) = TC(2)  
      TRCZ(5,L) = TC(3)

      TRR(5,L) = MAX(SUM((S1-TC)**2),
     &               SUM((S2-TC)**2),
     &               SUM((S3-TC)**2),
     &               SUM((S4-TC)**2))
      TRR(5,L) = SQRT(TRR(5,L))


C ... Surface elements splitted into four triangles,
C ... normals, center points and radiuses      

      DO IT = 1,4

      JT = IT + 1 ; IF(JT > 4) JT = 1 
         
      S1(1) = SKINX(IT,L)
      S1(2) = SKINY(IT,L)
      S1(3) = SKINZ(IT,L)

      S2(1) = SKINX(JT,L)
      S2(2) = SKINY(JT,L)
      S2(3) = SKINZ(JT,L)

      S3(1) = SKINXCP(L)
      S3(2) = SKINYCP(L)
      S3(3) = SKINZCP(L)
         
      TRNX(IT,L)=(S2(2)-S1(2))*(S3(3)-S1(3))-(S2(3)-S1(3))*(S3(2)-S1(2))   
      TRNY(IT,L)=(S2(3)-S1(3))*(S3(1)-S1(1))-(S2(1)-S1(1))*(S3(3)-S1(3))  
      TRNZ(IT,L)=(S2(1)-S1(1))*(S3(2)-S1(2))-(S2(2)-S1(2))*(S3(1)-S1(1)) 
         
      TC = (S1 + S2 + S3) / 3.0  
         
      TRCX(IT,L) = TC(1) 
      TRCY(IT,L) = TC(2)  
      TRCZ(IT,L) = TC(3)

      TRR(IT,L) = MAX(SUM((S1-TC)**2),SUM((S2-TC)**2),SUM((S3-TC)**2))
      TRR(IT,L) = SQRT(TRR(IT,L))

      ENDDO
      
      ENDDO 

C ... Occupied voxels in the coarse sorting grid. 

      NOCCUPIED1 = 0
 
      DO K = 1,NZBOX1
         DO J = 1,NYBOX1
            DO I = 1,NXBOX1
               
               IF(IBOXES1(I,J,K,1) == 0) CYCLE
               
               NOCCUPIED1 = NOCCUPIED1 + 1
               
            ENDDO
         ENDDO
      ENDDO

      ALLOCATE(OCCUPIED1(3,NOCCUPIED1),STAT=IERRCODE)
      ALLOCATE(CLOSEST1(3,NOCCUPIED1),STAT=IERRCODE)
      
      M = 0
      
      DO K = 1,NZBOX1
         DO J = 1,NYBOX1
            DO I = 1,NXBOX1
               
               IF(IBOXES1(I,J,K,1) > 0) THEN
                  M = M + 1
                  OCCUPIED1(1,M) = I
                  OCCUPIED1(2,M) = J
                  OCCUPIED1(3,M) = K
               ENDIF
               
            ENDDO
         ENDDO
      ENDDO

C ... Block loop runs here

      DO N=1,NBLOCK

         IG1 = IG(MG,N)
         NT  = NTOT(MG,N)
         CALL CALCDISTS(N,XC(IG1),YC(IG1),ZC(IG1),NT,DISTW(IG1),
     &        LOCDIS(IG1),NTSUMI,NTSUM,NTSUMG,KPROS,TWODIM,
     &        IBOXES1,NXBOX1,NYBOX1,NZBOX1,KSC1,
     &        IBOXES2,NXBOX2,NYBOX2,NZBOX2,KSC2,
     &        IBOXES3,NXBOX3,NYBOX3,NZBOX3,KSC3,
     &        NSKIN,SKINXCP,SKINYCP,SKINZCP,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,
     &        TRNX,TRNY,TRNZ,TRCX,TRCY,TRCZ,TRR,DXR,DYR,DZR,
     &        XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,BLANK(IG1))

      ENDDO

      DEALLOCATE(OCCUPIED1,CLOSEST1)
      DEALLOCATE(KSC1,IBOXES1,KSC2,IBOXES2,KSC3,IBOXES3)
      DEALLOCATE(SKINXCP,SKINYCP,SKINZCP)
      DEALLOCATE(TRNX,TRNY,TRNZ,TRCX,TRCY,TRCZ,TRR)

      RETURN      
      END SUBROUTINE LOOPBLOCKS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CALCDISTS(N,XC,YC,ZC,NT,DISTW,
     &        LOCDIS,NTSUMI,NTSUM,NTSUMG,KPROS,TWODIM,
     &        IBOXES1,NXBOX1,NYBOX1,NZBOX1,KSC1,
     &        IBOXES2,NXBOX2,NYBOX2,NZBOX2,KSC2,
     &        IBOXES3,NXBOX3,NYBOX3,NZBOX3,KSC3,
     &        NSKIN,SKINXCP,SKINYCP,SKINZCP,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,
     &        TRNX,TRNY,TRNZ,TRCX,TRCY,TRCZ,TRR,DXR,DYR,DZR,
     &        XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,BLANK)


C ... Calculate true wall distances within one block.

      USE MAIN_ARRAYS, ONLY : SKINX, SKINY, SKINZ, IBSKIN, ITSKIN, 
     &                        NPROCE, NCHIMT, IMAX, JMAX, KMAX    
      USE INTEGERS, ONLY    : IPRO 
      USE NS3CO, ONLY       : PARALLEL, TIMEL
      USE MPI

      IMPLICIT NONE

      INTEGER :: NT, NSKIN, I, J, K, L, N, NG, IJ
      INTEGER :: NOCCUPIED, NCLOSEST
      INTEGER :: NOCCUPIED1, NCLOSEST1
      INTEGER :: II, JJ, KK, IP, JP, KP
      INTEGER :: IERRCODE, ISIDE, ISIDEL, ITYPEL
      INTEGER :: LBS, LBE, IC, LC 
      INTEGER :: IL, IU, JL, JU, KL, KU

      INTEGER :: NXBOX1,  NYBOX1,  NZBOX1
      INTEGER :: NXBOX2,  NYBOX2,  NZBOX2
      INTEGER :: NXBOX3,  NYBOX3,  NZBOX3
      INTEGER(KIND=8) :: NTSUMI, NTSUM, NTSUMG, JPROS, KPROS
 
      INTEGER, DIMENSION(NXBOX1,NYBOX1,NZBOX1,2) :: IBOXES1
      INTEGER, DIMENSION(NXBOX2,NYBOX2,NZBOX2,2) :: IBOXES2
      INTEGER, DIMENSION(NXBOX3,NYBOX3,NZBOX3,2) :: IBOXES3
      INTEGER, DIMENSION(3,NOCCUPIED1) :: OCCUPIED1, CLOSEST1
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OCCUPIED, CLOSEST
      INTEGER, DIMENSION(*) :: KSC1, KSC2, KSC3, LOCDIS

      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
      REAL :: DXR, DYR, DZR, DS
      REAL :: DISTW(*), BLANK(*)

      REAL, DIMENSION(NSKIN)   :: SKINXCP, SKINYCP, SKINZCP
      REAL, DIMENSION(4,NSKIN) :: TRNX, TRNY, TRNZ 
      REAL, DIMENSION(5,NSKIN) :: TRCX, TRCY, TRCZ, TRR 
      REAL, DIMENSION(3) :: C1, C2, C3, C4, C5, PP, QQ
      REAL, DIMENSION(3) :: T1N, T2N, T3N, T4N
      REAL, DIMENSION(3) :: T1C, T2C, T3C, T4C, T5C
      REAL, DIMENSION(4) :: TR
      REAL, DIMENSION(*) :: XC, YC, ZC
      REAL :: DD, DISTWL

      LOGICAL :: MASTER, OVERLA, TWODIM

      LOGICAL, ALLOCATABLE, DIMENSION(:) :: DONE

      
      MASTER = IPRO == 1

      NG = NPROCE(N+1,IPRO)  ! Global block number

      ALLOCATE(DONE(NSKIN),STAT=IERRCODE)

      IF(IERRCODE > 0) THEN                                          
         WRITE(*,*)'CALCDISTS:  Not enough memory, aborting ...'             
         STOP                                                           
      ENDIF                                                    
      
C ... Cell center point loop starts here

      DO IC=1,NT
        
         JPROS = 100-100*(NTSUM-NTSUMI)/NTSUM
         IF(JPROS >= KPROS) THEN
            IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERRCODE)
            IF(KPROS == 0 .AND. MASTER) write(*,*)
            IF(MASTER) WRITE(*,'(I5," % completed")') KPROS
            IF(KPROS == 100 .AND. MASTER) write(*,*)
            KPROS = KPROS + 5
         ENDIF

         NTSUMI = NTSUMI + 1
         
         DONE(1:NSKIN) = .FALSE.

         DISTWL = DISTW(IC)
         BLANK(IC) = 1.0

         PP(1) = XC(IC); PP(2) = YC(IC); PP(3) = ZC(IC)


C ... Closest voxels search on the start level.

         CALL CLOSEST_VOXELS(PP,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &        NXBOX1,NYBOX1,NZBOX1,IBOXES1,
     &        OCCUPIED1,NOCCUPIED1,CLOSEST1,NCLOSEST1)


C ... Closest voxels search on the 4 x finer level.

         DS = 4.0

         ALLOCATE(OCCUPIED(3,64*NOCCUPIED1),STAT=IERRCODE)

         NOCCUPIED = 0
 
         DO L=1,NCLOSEST1

            I = CLOSEST1(1,L)
            J = CLOSEST1(2,L)
            K = CLOSEST1(3,L)

            IU = 4*I
            IL = IU-3
            JU = 4*J
            JL = JU-3
            KU = 4*K
            KL = KU-3

            DO KK=KL,KU
               DO JJ=JL,JU
                  DO II=IL,IU
               
                     IF(IBOXES2(II,JJ,KK,1) == 0) CYCLE

                     NOCCUPIED = NOCCUPIED + 1

                     OCCUPIED(1,NOCCUPIED) = II
                     OCCUPIED(2,NOCCUPIED) = JJ
                     OCCUPIED(3,NOCCUPIED) = KK

                  ENDDO
               ENDDO
            ENDDO

         ENDDO

         ALLOCATE(CLOSEST(3,NOCCUPIED),STAT=IERRCODE)

         CALL CLOSEST_VOXELS(PP,XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,
     &        NXBOX2,NYBOX2,NZBOX2,IBOXES2,
     &        OCCUPIED,NOCCUPIED,CLOSEST,NCLOSEST)

         DEALLOCATE(OCCUPIED)


C ... Closest voxels search on the 16 x finer level.

         DS = 16.0

         ALLOCATE(OCCUPIED(3,64*NOCCUPIED),STAT=IERRCODE)

         NOCCUPIED = 0
 
         DO L=1,NCLOSEST

            I = CLOSEST(1,L)
            J = CLOSEST(2,L)
            K = CLOSEST(3,L)

            IU = 4*I
            IL = IU-3
            JU = 4*J
            JL = JU-3
            KU = 4*K
            KL = KU-3

            DO KK=KL,KU
               DO JJ=JL,JU
                  DO II=IL,IU
               
                     IF(IBOXES3(II,JJ,KK,1) == 0) CYCLE

                     NOCCUPIED = NOCCUPIED + 1

                     OCCUPIED(1,NOCCUPIED) = II
                     OCCUPIED(2,NOCCUPIED) = JJ
                     OCCUPIED(3,NOCCUPIED) = KK

                  ENDDO
               ENDDO
            ENDDO

         ENDDO

         DEALLOCATE(CLOSEST)         
         ALLOCATE(CLOSEST(3,NOCCUPIED),STAT=IERRCODE)

         CALL CLOSEST_VOXELS(PP,XMIN,YMIN,ZMIN,DS*DXR,DS*DYR,DS*DZR,
     &        NXBOX3,NYBOX3,NZBOX3,IBOXES3,
     &        OCCUPIED,NOCCUPIED,CLOSEST,NCLOSEST)

         DEALLOCATE(OCCUPIED)

         
C ... Scan through the skin elements in the closest sorting cells.

         DO J=1,NCLOSEST
            
           II = CLOSEST(1,J)
           JJ = CLOSEST(2,J)
           KK = CLOSEST(3,J)

           LBS = IBOXES3(II,JJ,KK,2)             ! Start
           LBE = LBS + IBOXES3(II,JJ,KK,1) - 1   ! End

           DO LC = LBS,LBE
          
           IJ = KSC3(LC)                         ! Skin element number

           IF(DONE(IJ)) CYCLE

C ... A sphere around the surface element lies outside the minimum
C ... distance so far?
           
           T5C(1) = TRCX(5,IJ); T5C(2) = TRCY(5,IJ); T5C(3) = TRCZ(5,IJ)
           IF(SQRT(SUM((PP-T5C)**2)) - TRR(5,IJ) > DISTWL) CYCLE      
           
           C1(1) = SKINX(1,IJ); C1(2) = SKINY(1,IJ); C1(3) = SKINZ(1,IJ) 
           C2(1) = SKINX(2,IJ); C2(2) = SKINY(2,IJ); C2(3) = SKINZ(2,IJ) 
           C3(1) = SKINX(3,IJ); C3(2) = SKINY(3,IJ); C3(3) = SKINZ(3,IJ) 
           C4(1) = SKINX(4,IJ); C4(2) = SKINY(4,IJ); C4(3) = SKINZ(4,IJ) 
           C5(1) = SKINXCP(IJ); C5(2) = SKINYCP(IJ); C5(3) = SKINZCP(IJ)

           T1N(1) = TRNX(1,IJ); T1N(2) = TRNY(1,IJ); T1N(3) = TRNZ(1,IJ)
           T2N(1) = TRNX(2,IJ); T2N(2) = TRNY(2,IJ); T2N(3) = TRNZ(2,IJ)
           T3N(1) = TRNX(3,IJ); T3N(2) = TRNY(3,IJ); T3N(3) = TRNZ(3,IJ)
           T4N(1) = TRNX(4,IJ); T4N(2) = TRNY(4,IJ); T4N(3) = TRNZ(4,IJ)

           T1C(1) = TRCX(1,IJ); T1C(2) = TRCY(1,IJ); T1C(3) = TRCZ(1,IJ)
           T2C(1) = TRCX(2,IJ); T2C(2) = TRCY(2,IJ); T2C(3) = TRCZ(2,IJ)
           T3C(1) = TRCX(3,IJ); T3C(2) = TRCY(3,IJ); T3C(3) = TRCZ(3,IJ)
           T4C(1) = TRCX(4,IJ); T4C(2) = TRCY(4,IJ); T4C(3) = TRCZ(4,IJ)

           TR(1) = TRR(1,IJ)
           TR(2) = TRR(2,IJ)
           TR(3) = TRR(3,IJ)
           TR(4) = TRR(4,IJ)

C ... Calculate the true distance.

           CALL PANELD(PP,C1,C2,C3,C4,C5,
     &                 T1N,T2N,T3N,T4N,T1C,T2C,T3C,T4C,TR,
     &                 QQ,DD,DISTWL,IJ,ISIDE)
          
C ...      DD = SUM((PP-C5)**2)  ! 'Wall distance' to the element center.
           
           DONE(IJ) = .TRUE.

           IF(DD < DISTWL .OR. DD == DISTWL .AND. ISIDE == 1) THEN
*           IF(DD < DISTWL) THEN
              DISTWL     = DD
              DISTW(IC)  = DISTWL
              LOCDIS(IC) = IBSKIN(IJ)
              ISIDEL     = ISIDE
              ITYPEL     = ITSKIN(IJ)
              OVERLA     = NCHIMT(NG) /= NCHIMT(LOCDIS(IC))
           ENDIF 

           ENDDO

         ENDDO             

         IF(ISIDEL == -1.AND.ITYPEL /= 19 .AND. OVERLA) BLANK(IC) = -8.0

         DEALLOCATE(CLOSEST)
         
      END DO  ! IC=1,NT cell center point loop

      DEALLOCATE(DONE)

      RETURN
      END SUBROUTINE CALCDISTS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CLOSEST_VOXELS(PP,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &           NXBOX,NYBOX,NZBOX,IBOXES,
     &           OCCUPIED,NOCCUPIED,CLOSEST,NCLOSEST)

      IMPLICIT NONE
      
      REAL, DIMENSION(3) :: PP

      REAL :: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, DXR, DYR, DZR

      INTEGER :: I1F, J1F, K1F
      INTEGER :: I3F, I3L, J3F, J3L, K3F, K3L
      INTEGER :: NXBOX,  NYBOX,  NZBOX
      INTEGER :: NOCCUPIED, NCLOSEST, ICLOSEST
      INTEGER :: ISAVE, JSAVE, KSAVE, IERRCODE

      INTEGER(KIND=8) :: I, J, K, L, II, JJ, KK, IP, JP, KP
      INTEGER(KIND=8) :: IDIS, JDIS, KDIS, IDISP1, JDISP1, KDISP1
      INTEGER(KIND=8) :: INTDISTMIN, INTDISTMINP1 
      INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: DOCCUPIED

      INTEGER, DIMENSION(NXBOX,NYBOX,NZBOX,2) :: IBOXES
      INTEGER, DIMENSION(3,NOCCUPIED), INTENT(IN)  :: OCCUPIED
      INTEGER, DIMENSION(3,NOCCUPIED), INTENT(OUT) :: CLOSEST

      LOGICAL :: ITISIN, NOSEARCH

      ALLOCATE(DOCCUPIED(NOCCUPIED),STAT=IERRCODE)

      IP = INT((PP(1)-XMIN)*DXR)
      JP = INT((PP(2)-YMIN)*DYR)
      KP = INT((PP(3)-ZMIN)*DZR)

      IF(PP(1) > XMIN) IP = IP + 1
      IF(PP(2) > YMIN) JP = JP + 1
      IF(PP(3) > ZMIN) KP = KP + 1
         

      ITISIN = 1 <= IP .AND. IP <= NXBOX
     &   .AND. 1 <= JP .AND. JP <= NYBOX
     &   .AND. 1 <= KP .AND. KP <= NZBOX

      NOSEARCH = .FALSE.
      IF(ITISIN) THEN
         I1F = IP
         J1F = JP
         K1F = KP
         NOSEARCH = IBOXES(I1F,J1F,K1F,1) > 0
      ENDIF

         
C ... Find the closest occupied sorting cells.          

C ... If we are lucky, the cell center point lies in an occupied sorting 
C ... cell (a cell that contains surface elements). The closest occupied 
C ... sorting cells are known without searching.

      IF(NOSEARCH) THEN

         I3F = MAX(I1F-2,1)
         J3F = MAX(J1F-2,1)
         K3F = MAX(K1F-2,1)
         I3L = MIN(I1F+2,NXBOX)
         J3L = MIN(J1F+2,NYBOX)
         K3L = MIN(K1F+2,NZBOX)  

C ... Put the closest sorting cell first.
            
         NCLOSEST = 1

         CLOSEST(1,NCLOSEST) = I1F
         CLOSEST(2,NCLOSEST) = J1F
         CLOSEST(3,NCLOSEST) = K1F
         
         DO KK = K3F,K3L
            DO JJ = J3F,J3L
               DO II = I3F,I3L
                    
                  IF(IBOXES(II,JJ,KK,1) == 0) CYCLE

                  NCLOSEST = NCLOSEST + 1

                  CLOSEST(1,NCLOSEST) = II
                  CLOSEST(2,NCLOSEST) = JJ
                  CLOSEST(3,NCLOSEST) = KK

               ENDDO
            ENDDO
         ENDDO

      ELSE

C ... The closest occupied sorting cells must be searched. At the first 
C ... loop determine only the (minimum + sqrt(3.0))**2 distance. 

         INTDISTMINP1 = HUGE(INTDISTMINP1)
           
         DO L=1,NOCCUPIED

            I = OCCUPIED(1,L); IDIS=(IP-I)**2; IDISP1=(ABS(IP-I)+1)**2
            J = OCCUPIED(2,L); JDIS=(JP-J)**2; JDISP1=(ABS(JP-J)+1)**2
            K = OCCUPIED(3,L); KDIS=(KP-K)**2; KDISP1=(ABS(KP-K)+1)**2

            DOCCUPIED(L) = IDIS + JDIS + KDIS
        
            INTDISTMINP1 = MIN(INTDISTMINP1,IDISP1+JDISP1+KDISP1)

         ENDDO
                      
C ... Now the minimum distance to the closest occupied sorting cells is known. 
C ... All further away cells can be skipped. Record the closest cells.
            
         INTDISTMIN = INTDISTMINP1

         ICLOSEST = 1
         NCLOSEST = 0

         DO L=1,NOCCUPIED
 
            IF(DOCCUPIED(L) <= INTDISTMINP1) THEN

               NCLOSEST = NCLOSEST + 1
               
               CLOSEST(1,NCLOSEST) = OCCUPIED(1,L)
               CLOSEST(2,NCLOSEST) = OCCUPIED(2,L)
               CLOSEST(3,NCLOSEST) = OCCUPIED(3,L)

               IF(DOCCUPIED(L) < INTDISTMIN) THEN
                  INTDISTMIN = DOCCUPIED(L)
                  ICLOSEST = NCLOSEST
               ENDIF

            ENDIF

         ENDDO
         
C ... Swap the closest occupied sorting cell to the top.            

         IF(ICLOSEST /= 1) THEN
            ISAVE = CLOSEST(1,ICLOSEST)
            JSAVE = CLOSEST(2,ICLOSEST)
            KSAVE = CLOSEST(3,ICLOSEST)
            CLOSEST(1,ICLOSEST) = CLOSEST(1,1)
            CLOSEST(2,ICLOSEST) = CLOSEST(2,1)
            CLOSEST(3,ICLOSEST) = CLOSEST(3,1)
            CLOSEST(1,1) = ISAVE
            CLOSEST(2,1) = JSAVE
            CLOSEST(3,1) = KSAVE
         ENDIF

      ENDIF

      DEALLOCATE(DOCCUPIED)
      
      RETURN
      END SUBROUTINE CLOSEST_VOXELS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PANELD(PP,C1,C2,C3,C4,C5,
     &                  T1N,T2N,T3N,T4N,T1C,T2C,T3C,T4C,TR,
     &                  Q,DD,DISTWL,IJ,ISIDE)

C ... Minimum distance (square) between a point and a four corner
C ... surface element. Split the four corner element first into 
C ... four triangles.

      USE MAIN_ARRAYS, ONLY : SKINXN, SKINYN, SKINZN

      IMPLICIT NONE

      REAL, DIMENSION(3) :: C1, C2, C3, C4, C5, PP, Q, H, N
      REAL, DIMENSION(3) :: T1N, T2N, T3N, T4N
      REAL, DIMENSION(3) :: T1C, T2C, T3C, T4C
      REAL, DIMENSION(4) :: TR
      REAL :: R1, R2, R3, R4, RINI, DD
      REAL :: DISTWL

      INTEGER :: ISIDE, IJ, ISTYPE, JSTYPE

      ISTYPE = 10

*      RINI = 1.0D+10 
      RINI = HUGE(RINI)

      DD = RINI; R1 = RINI; R2 = RINI; R3 = RINI; R4 = RINI

      IF(SQRT(SUM((PP-T1C)**2)) - TR(1) < DISTWL) THEN      
      CALL POINT_TRIANGLE(PP,C1,C2,C5,T1N,T1C,TR(1),H,R1,1,JSTYPE)
         IF(R1 < DD) THEN
            DD = R1
            Q  = H
            ISTYPE = JSTYPE
         ENDIF
      ENDIF

      IF(SQRT(SUM((PP-T2C)**2)) - TR(2) < DISTWL) THEN      
      CALL POINT_TRIANGLE(PP,C2,C3,C5,T2N,T2C,TR(2),H,R2,2,JSTYPE)
         IF(R2 < DD) THEN
            DD = R2
            Q  = H
            ISTYPE = JSTYPE
         ENDIF
      ENDIF

      IF(SQRT(SUM((PP-T3C)**2)) - TR(3) < DISTWL) THEN      
      CALL POINT_TRIANGLE(PP,C3,C4,C5,T3N,T3C,TR(3),H,R3,3,JSTYPE)
         IF(R3 < DD) THEN
            DD = R3
            Q  = H
            ISTYPE = JSTYPE
         ENDIF
      ENDIF

      IF(SQRT(SUM((PP-T4C)**2)) - TR(4) < DISTWL) THEN      
      CALL POINT_TRIANGLE(PP,C4,C1,C5,T4N,T4C,TR(4),H,R4,4,JSTYPE)
         IF(R4 < DD) THEN
            DD = R4
            Q  = H
            ISTYPE = JSTYPE
         ENDIF
      ENDIF

      DD = SQRT(DD)
         
C ... Does the point PP lie inside solid body?

      IF(ISTYPE < 10) THEN

         N(1) = SKINXN(ISTYPE,IJ)
         N(2) = SKINYN(ISTYPE,IJ)
         N(3) = SKINZN(ISTYPE,IJ)

         IF(DOT_PRODUCT(N,PP-Q) < 0.0) THEN
            ISIDE = -1
         ELSE
            ISIDE =  1
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE PANELD
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE POINT_TRIANGLE(PP,D1,D2,D3,TN,TC,TR,
     &                          Q,RR,ITRI,JSTYPE)      

      IMPLICIT NONE

      INTEGER :: ITRI, ISIDE, ICODE, JSTYPE 

      REAL, DIMENSION(3) :: PP, D1, D2, D3, Q, H, TN, TC
      REAL :: R1, R2, R3, R4, RR, DD, TR

      DD = RR; R1 = RR; R2 = RR; R3 = RR; R4 = RR

      
C ... Normal distance (square) to the triangle surface.

      CALL POINT_PLANE(PP,D1,D2,D3,TN,TC,TR,DD,H,ICODE)

      IF(ICODE == 1) THEN
         JSTYPE = 1
         Q      = H
         RR     = DD
         RETURN
      ENDIF

C ... Distances (squares) to the triangle border segments.

      CALL POINT_SEGMENT(PP,D1,D2,H,R1)
      IF(R1 < RR) THEN
         JSTYPE = ITRI + 5
         Q      = H
         RR     = R1
      ENDIF
      CALL POINT_SEGMENT(PP,D2,D3,H,R2)
      IF(R2 < RR) THEN
         JSTYPE = 1
         Q      = H
         RR     = R2
      ENDIF

C ... Distances (squares) to the triangle corner points. 
      
      CALL POINT_POINT(PP,D1,R3)
      IF(R3 < RR) THEN
         JSTYPE = ITRI + 1
         Q      = D1
         RR     = R3
      ENDIF

      IF(ITRI == 1) THEN
         CALL POINT_POINT(PP,D3,R4)
         IF(R4 < RR) THEN
            JSTYPE = 1
            Q      = D3
            RR     = R4
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE POINT_TRIANGLE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE POINT_PLANE(P,S1,S2,S3,TN,TC,TR,D,Q,ICODE)

C ... This subroutine calculates the intersection point Q of a line 
C     segment defined by normal N and point P and a plane defined by 
C     three points S1, S2 and S3. Subroutine also checks whether the
C     intersection point lies within the triangle or not.
C
C       ICODE = 1:  The line and the plane have an intersection point
C                   which lies within the triangle determined by the
C                   points S1, S2, and S3.
C
C       ICODE = 2:  The line and the plane have an intersection point but
C                   it does not lie within the triangle determined by
C                   the points S1, S2, and S3.
C ...               

      IMPLICIT NONE

      REAL, PARAMETER    :: EPS=1.0E-25

      REAL, DIMENSION(3) :: P, Q, TN, TC, S1, S2, S3,
     &                      S12, S23, S31, S1Q, S2Q, S3Q

      REAL :: A, B, C, D, DD, TR
      REAL :: D1, D2, D3, D4

      INTEGER :: ICODE

C ... Check the validity of the triangle geometry

       IF(SUM(TN**2) < EPS) RETURN


C ... Intersection point of a line and a plane

      Q = P + DOT_PRODUCT(TN,S1-P)/DOT_PRODUCT(TN,TN)*TN      

      IF(SQRT(SUM((Q-TC)**2)) > TR) THEN
         ICODE = 2
         RETURN
      ENDIF


C ... Check if the line intersects the triangle

      S1Q = S1 - Q
      S2Q = S2 - Q
      S3Q = S3 - Q

      CALL CROSS_PRODUCT(S1Q,S2Q,S12)
      CALL CROSS_PRODUCT(S2Q,S3Q,S23)

      A = DOT_PRODUCT(TN,S12) 
      B = DOT_PRODUCT(TN,S23) 

      IF(A*B <= 0.0) THEN
         ICODE = 2
         RETURN
      ENDIF

      CALL CROSS_PRODUCT(S3Q,S1Q,S31)

      C = DOT_PRODUCT(TN,S31) 

      IF(A*C <= 0.0) THEN
         ICODE = 2
      ELSE
         ICODE = 1
      ENDIF
         
      D = SUM((P-Q)**2)  ! Perpendicular distance

      RETURN
      END SUBROUTINE POINT_PLANE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE POINT_SEGMENT(P,S1,S2,H,D)

C ... Distance (square) between a point and a line segment. Distance 
C ... D is calculated only if the perpendicular distance from the point 
C ... to the line lies between S1 and S2.

      IMPLICIT NONE

      REAL, PARAMETER    :: EPS=1.0E-25
      REAL, DIMENSION(3) :: P, S1, S2, C, H, S2S1, S1P, PS1
      REAL :: T, D, D21

      LOGICAL :: OK
      
      S2S1 = S2 - S1
      S1P  = S1 - P
      PS1  = P  - S1
      D21  = SUM(S2S1**2)
      
      IF(D21 < EPS) RETURN

      T = -DOT_PRODUCT(S1P,S2S1) / D21

      OK = T > 0.0 .AND. T < 1.0

      IF(OK) THEN
 
         CALL CROSS_PRODUCT(S2S1,PS1,C)
                      
         D = SUM(C**2) / D21
 
         H = S1 + T*S2S1

      ENDIF

      
      RETURN
      END SUBROUTINE POINT_SEGMENT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE POINT_POINT(P,S1,D)

C ... Distance (square) between two points.

      IMPLICIT NONE

      REAL, DIMENSION(3) :: P, S1
      REAL :: D

      D = SUM((P-S1)**2)

      RETURN
      END SUBROUTINE POINT_POINT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOXEL_EXACT(NSKIN,SKINX,SKINY,SKINZ,KSC,NCB,
     &                       IBOXES,NXBOX,NYBOX,NZBOX,
     &                       XMIN,YMIN,ZMIN,DXR,DYR,DZR,MT)

C ... This subroutine shares the skin elements into groups in order to 
C ... accelerate the search process. If MT=1, only the space
C ... requirement is calculated. In case MT=2, the elements are
C ... actually shared. This special version marks only the grouping cells
C ... that are actually hit by the element edges.
      

      IMPLICIT NONE

      INTEGER :: NSKIN, NXBOX, NYBOX, NZBOX
      INTEGER, DIMENSION(NXBOX,NYBOX,NZBOX,2) :: IBOXES
      INTEGER, DIMENSION(*) :: KSC
      INTEGER :: I1, J1, K1, MT, N, ISKIN, NCB, ICP, INODE, ICODE
      INTEGER :: IB1, IB2, IB3, IB4
      INTEGER :: JB1, JB2, JB3, JB4
      INTEGER :: KB1, KB2, KB3, KB4
      INTEGER :: IBS, IBE, JBS, JBE, KBS, KBE

      REAL, DIMENSION(4,NSKIN) :: SKINX, SKINY, SKINZ
      REAL :: X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3, X4,Y4,Z4
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR

      IF(MT == 1) THEN     ! Initialize the box pointer array
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,1) = 0
                  IBOXES(I1,J1,K1,2) = 0
               ENDDO
            ENDDO
         ENDDO

         NCB = 0  ! Total number of elements in groups. Note that one 
                  ! element can belong to several groups.
      ENDIF


      IF(MT == 2) THEN     ! Starting pointer for each group ('box')
         N = 1
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,2) = N
                  N = N + IBOXES(I1,J1,K1,1)
                  IBOXES(I1,J1,K1,1) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO ISKIN=1,NSKIN

         X1 = SKINX(1,ISKIN)
         Y1 = SKINY(1,ISKIN)
         Z1 = SKINZ(1,ISKIN)

         X2 = SKINX(2,ISKIN)
         Y2 = SKINY(2,ISKIN)
         Z2 = SKINZ(2,ISKIN)

         X3 = SKINX(3,ISKIN)
         Y3 = SKINY(3,ISKIN)
         Z3 = SKINZ(3,ISKIN)

         X4 = SKINX(4,ISKIN)
         Y4 = SKINY(4,ISKIN)
         Z4 = SKINZ(4,ISKIN)
         
         IB1 = INT((X1-XMIN)*DXR) + 1
         IB2 = INT((X2-XMIN)*DXR) + 1
         IB3 = INT((X3-XMIN)*DXR) + 1
         IB4 = INT((X4-XMIN)*DXR) + 1
         
         JB1 = INT((Y1-YMIN)*DYR) + 1
         JB2 = INT((Y2-YMIN)*DYR) + 1
         JB3 = INT((Y3-YMIN)*DYR) + 1
         JB4 = INT((Y4-YMIN)*DYR) + 1
         
         KB1 = INT((Z1-ZMIN)*DZR) + 1
         KB2 = INT((Z2-ZMIN)*DZR) + 1
         KB3 = INT((Z3-ZMIN)*DZR) + 1
         KB4 = INT((Z4-ZMIN)*DZR) + 1
         
         IBS = MIN(IB1,IB2,IB3,IB4)
         IBE = MAX(IB1,IB2,IB3,IB4)
         JBS = MIN(JB1,JB2,JB3,JB4)
         JBE = MAX(JB1,JB2,JB3,JB4)
         KBS = MIN(KB1,KB2,KB3,KB4)
         KBE = MAX(KB1,KB2,KB3,KB4)
         IBS = MAX(IBS,1)
         IBE = MIN(IBE,NXBOX)
         JBS = MAX(JBS,1)
         JBE = MIN(JBE,NYBOX)
         KBS = MAX(KBS,1)
         KBE = MIN(KBE,NZBOX)
            
         DO I1=IBS,IBE
            DO J1=JBS,JBE
               DO K1=KBS,KBE

                  IF(IBS == IBE .AND. JBS == JBE .AND. KBS == KBE) THEN
                     ICODE = 0
                  ELSE
                  CALL CHECK_VOXEL1(I1,J1,K1,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ICODE)
                  IF(ICODE /= 0) 
     &            CALL CHECK_VOXEL2(I1,J1,K1,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &                 X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ICODE)
                  ENDIF
                  
                  IF(MT == 1 .AND. ICODE == 0) THEN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                     NCB = NCB + 1
                  ELSEIF(MT == 2 .AND. ICODE == 0) THEN
                     ICP = IBOXES(I1,J1,K1,1)+IBOXES(I1,J1,K1,2)
                     KSC(ICP) = ISKIN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE BOXEL_EXACT
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE BOXEL_WETTED(IBL,KSC,NCB,IBOXES,NXBOX,NYBOX,NZBOX,
     &                        XMIN,YMIN,ZMIN,DXR,DYR,DZR,MT)

C ... This subroutine shares the wetted elements into groups in order to 
C ... accelerate the search process. If MT=1, only the space
C ... requirement is calculated. In case MT=2, the elements are
C ... actually shared. This special version marks only the grouping cells
C ... that are actually hit by the element edges.
      
      USE MPCCIVARS, ONLY : fnelems, felemnodes, fnodecoor_orig,
     &                      fnodechimt

      USE MAIN_ARRAYS, ONLY : NPROCE, NCHIMT

      USE INTEGERS, ONLY : IPRO
      
      IMPLICIT NONE

      INTEGER :: NXBOX, NYBOX, NZBOX
      INTEGER, DIMENSION(NXBOX,NYBOX,NZBOX,2) :: IBOXES
      INTEGER, DIMENSION(*) :: KSC
      INTEGER :: I1, J1, K1, MT, N, ISKIN, NCB, ICP, INODE, ICODE
      INTEGER :: IB1, IB2, IB3, IB4
      INTEGER :: JB1, JB2, JB3, JB4
      INTEGER :: KB1, KB2, KB3, KB4
      INTEGER :: IBS, IBE, JBS, JBE, KBS, KBE, IBL, IGL 

      REAL :: X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3, X4,Y4,Z4
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR

      IGL = NPROCE(1+IBL,IPRO)
      
      IF(MT == 1) THEN     ! Initialize the box pointer array
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,1) = 0
                  IBOXES(I1,J1,K1,2) = 0
               ENDDO
            ENDDO
         ENDDO

         NCB = 0  ! Total number of elements in groups. Note that one 
                  ! element can belong to several groups.
      ENDIF


      IF(MT == 2) THEN     ! Starting pointer for each group ('box')
         N = 1
         DO I1 = 1,NXBOX
            DO J1 = 1,NYBOX
               DO K1 = 1,NZBOX
                  IBOXES(I1,J1,K1,2) = N
                  N = N + IBOXES(I1,J1,K1,1)
                  IBOXES(I1,J1,K1,1) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO ISKIN=1,fnelems

************************************************************************
         if(fnodechimt(ISKIN) /= NCHIMT(IGL)) CYCLE         
************************************************************************
         
         X1 = fnodecoor_orig(felemnodes(ISKIN,1),1)
         Y1 = fnodecoor_orig(felemnodes(ISKIN,1),2) 
         Z1 = fnodecoor_orig(felemnodes(ISKIN,1),3)

         X2 = fnodecoor_orig(felemnodes(ISKIN,2),1)
         Y2 = fnodecoor_orig(felemnodes(ISKIN,2),2)
         Z2 = fnodecoor_orig(felemnodes(ISKIN,2),3)

         X3 = fnodecoor_orig(felemnodes(ISKIN,3),1)
         Y3 = fnodecoor_orig(felemnodes(ISKIN,3),2)
         Z3 = fnodecoor_orig(felemnodes(ISKIN,3),3)

         X4 = fnodecoor_orig(felemnodes(ISKIN,4),1)
         Y4 = fnodecoor_orig(felemnodes(ISKIN,4),2)
         Z4 = fnodecoor_orig(felemnodes(ISKIN,4),3)

         IB1 = INT((X1-XMIN)*DXR) + 1
         IB2 = INT((X2-XMIN)*DXR) + 1
         IB3 = INT((X3-XMIN)*DXR) + 1
         IB4 = INT((X4-XMIN)*DXR) + 1
         
         JB1 = INT((Y1-YMIN)*DYR) + 1
         JB2 = INT((Y2-YMIN)*DYR) + 1
         JB3 = INT((Y3-YMIN)*DYR) + 1
         JB4 = INT((Y4-YMIN)*DYR) + 1
         
         KB1 = INT((Z1-ZMIN)*DZR) + 1
         KB2 = INT((Z2-ZMIN)*DZR) + 1
         KB3 = INT((Z3-ZMIN)*DZR) + 1
         KB4 = INT((Z4-ZMIN)*DZR) + 1
         
         IBS = MIN(IB1,IB2,IB3,IB4)
         IBE = MAX(IB1,IB2,IB3,IB4)
         JBS = MIN(JB1,JB2,JB3,JB4)
         JBE = MAX(JB1,JB2,JB3,JB4)
         KBS = MIN(KB1,KB2,KB3,KB4)
         KBE = MAX(KB1,KB2,KB3,KB4)
         IBS = MAX(IBS,1)
         IBE = MIN(IBE,NXBOX)
         JBS = MAX(JBS,1)
         JBE = MIN(JBE,NYBOX)
         KBS = MAX(KBS,1)
         KBE = MIN(KBE,NZBOX)
            
         DO I1=IBS,IBE
            DO J1=JBS,JBE
               DO K1=KBS,KBE

                  IF(IBS == IBE .AND. JBS == JBE .AND. KBS == KBE) THEN
                     ICODE = 0
                  ELSE
                  CALL CHECK_VOXEL1(I1,J1,K1,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ICODE)
                  IF(ICODE /= 0) 
     &            CALL CHECK_VOXEL2(I1,J1,K1,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &                 X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,ICODE)
                  ENDIF
                  
                  IF(MT == 1 .AND. ICODE == 0) THEN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                     NCB = NCB + 1
                  ELSEIF(MT == 2 .AND. ICODE == 0) THEN
                     ICP = IBOXES(I1,J1,K1,1)+IBOXES(I1,J1,K1,2)
                     KSC(ICP) = ISKIN
                     IBOXES(I1,J1,K1,1) = IBOXES(I1,J1,K1,1) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE BOXEL_WETTED
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CHECK_VOXEL1(I,J,K,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &     XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,ICODE)

C ... Check if the surface element edges hit the voxel faces (ICODE == 0).  

      IMPLICIT NONE

      INTEGER :: I, J, K, ICODE, IFACE, IEDGE
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR, DX, DY, DZ
      REAL, PARAMETER :: EPS=1.0E-12
      REAL :: XX1,YY1,ZZ1, XX2,YY2,ZZ2,XX3,YY3,ZZ3, XX4,YY4,ZZ4
      REAL :: XC1,YC1,ZC1, XC8,YC8,ZC8
      REAL :: X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3
      REAL :: XN,YN,ZN,XI,YI,ZI,DD,CE,CD,DP 

C ... Voxel edge lengths.
      
      DX = 1.0/DXR
      DY = 1.0/DYR
      DZ = 1.0/DZR      

C ... Voxel (two opposite) corners.

      XC1 = XMIN + (I-1)*DX
      YC1 = YMIN + (J-1)*DY
      ZC1 = ZMIN + (K-1)*DZ

      XC8 = XMIN + (I-0)*DX
      YC8 = YMIN + (J-0)*DY
      ZC8 = ZMIN + (K-0)*DZ
     
      ICODE = -1
      
      DO IEDGE = 1,4  ! Trapezoid edges.

         SELECT CASE(IEDGE)
            CASE(1)
               X1 = XX1
               Y1 = YY1
               Z1 = ZZ1
               X2 = XX2
               Y2 = YY2
               Z2 = ZZ2
            CASE(2)
               X1 = XX2
               Y1 = YY2
               Z1 = ZZ2
               X2 = XX3
               Y2 = YY3
               Z2 = ZZ3
            CASE(3)
               X1 = XX3
               Y1 = YY3
               Z1 = ZZ3
               X2 = XX4
               Y2 = YY4
               Z2 = ZZ4
            CASE(4)
               X1 = XX4
               Y1 = YY4
               Z1 = ZZ4
               X2 = XX1
               Y2 = YY1
               Z2 = ZZ1
         END SELECT

C ... Edge length.
            
         DD = SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1))
         IF(DD < EPS) CYCLE 

         DO IFACE = 1,6

C ... One point on the face and the face normal are needed.

            SELECT CASE(IFACE)
               CASE(1)
                  X3 = XC1
                  Y3 = YC1
                  Z3 = ZC1
                  XN = 1.0
                  YN = 0.0
                  ZN = 0.0
               CASE(2)
                  X3 = XC1
                  Y3 = YC1
                  Z3 = ZC1
                  XN = 0.0
                  YN = 1.0
                  ZN = 0.0
               CASE(3)
                  X3 = XC1
                  Y3 = YC1
                  Z3 = ZC1
                  XN = 0.0
                  YN = 0.0
                  ZN = 1.0
               CASE(4)
                  X3 = XC8
                  Y3 = YC8
                  Z3 = ZC8
                  XN = 1.0
                  YN = 0.0
                  ZN = 0.0
               CASE(5)
                  X3 = XC8
                  Y3 = YC8
                  Z3 = ZC8
                  XN = 0.0
                  YN = 1.0
                  ZN = 0.0
               CASE(6)
                  X3 = XC8
                  Y3 = YC8
                  Z3 = ZC8
                  XN = 0.0
                  YN = 0.0
                  ZN = 1.0
            END SELECT

            CE = XN*(X3-X1)+YN*(Y3-Y1)+ZN*(Z3-Z1)
            CD = XN*(X2-X1)+YN*(Y2-Y1)+ZN*(Z2-Z1)

            CD = SIGN(MAX(ABS(CD),1.E-16),CD)

C ... Intersection point of the trapezoid edge and the voxel face.

            XI = X1 + CE/CD*(X2-X1)
            YI = Y1 + CE/CD*(Y2-Y1)
            ZI = Z1 + CE/CD*(Z2-Z1)

            IF(IFACE == 1 .OR. IFACE == 4) THEN
               IF(YI < YC1 .OR. YI > YC8) CYCLE
               IF(ZI < ZC1 .OR. ZI > ZC8) CYCLE
            ENDIF
            IF(IFACE == 2 .OR. IFACE == 5) THEN
               IF(XI < XC1 .OR. XI > XC8) CYCLE
               IF(ZI < ZC1 .OR. ZI > ZC8) CYCLE
            ENDIF
            IF(IFACE == 3 .OR. IFACE == 6) THEN
               IF(XI < XC1 .OR. XI > XC8) CYCLE
               IF(YI < YC1 .OR. YI > YC8) CYCLE
            ENDIF

C ... Do the edge end points lie on the same side of the voxel face?
                     
            DP = (X1-XI)*(X2-XI)
     &         + (Y1-YI)*(Y2-YI)                     
     &         + (Z1-ZI)*(Z2-ZI)                     

            IF(DP > 0.0) CYCLE

            ICODE = 0
            GOTO 100
                  
         ENDDO

      ENDDO

 100  CONTINUE
      
      RETURN
      END SUBROUTINE CHECK_VOXEL1
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CHECK_VOXEL2(I,J,K,XMIN,YMIN,ZMIN,DXR,DYR,DZR,
     &     XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,ICODE)

C ... Check if the voxel edges (twelve) hit the trapezoid divided
C ... into four triangles.       

      IMPLICIT NONE

      INTEGER :: I, J, K, ICODE, ITRI, IEDGE
      REAL :: XMIN, YMIN, ZMIN, DXR, DYR, DZR, DX, DY, DZ
      REAL, PARAMETER :: TL=1.0E-3
      REAL :: XX1,YY1,ZZ1, XX2,YY2,ZZ2,
     &        XX3,YY3,ZZ3, XX4,YY4,ZZ4
      REAL :: XC1,YC1,ZC1, XC2,YC2,ZC2,
     &        XC3,YC3,ZC3, XC4,YC4,ZC4,
     &        XC5,YC5,ZC5, XC6,YC6,ZC6,
     &        XC7,YC7,ZC7, XC8,YC8,ZC8
      REAL ::
     &     XI,YI,ZI, X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3, X4,Y4,Z4, X5,Y5,Z5

C ... Voxel edge lengths.
      
      DX = 1.0/DXR
      DY = 1.0/DYR
      DZ = 1.0/DZR      

C ... Trapezoid centre point.

      X5 = (XX1 + XX2 + XX3 + XX4)/4.0
      Y5 = (YY1 + YY2 + YY3 + YY4)/4.0
      Z5 = (ZZ1 + ZZ2 + ZZ3 + ZZ4)/4.0

C ... Voxel corners.

      XC1 = XMIN + (I-1)*DX
      XC2 = XMIN + (I-0)*DX
      XC3 = XC1
      XC4 = XC2
      XC5 = XC1
      XC6 = XC2
      XC7 = XC1
      XC8 = XC2
 
      YC1 = YMIN + (J-1)*DY
      YC2 = YC1
      YC3 = YMIN + (J-0)*DY
      YC4 = YC3
      YC5 = YC1
      YC6 = YC2
      YC7 = YC3
      YC8 = YC4
 
      ZC1 = ZMIN + (K-1)*DZ
      ZC2 = ZC1
      ZC3 = ZC1
      ZC4 = ZC1
      ZC5 = ZMIN + (K-0)*DZ
      ZC6 = ZC5
      ZC7 = ZC5
      ZC8 = ZC5
     
      ICODE = -1
      
      DO ITRI = 1,4  ! Triangle corners.

         SELECT CASE(ITRI)
            CASE(1)
               X3 = XX1
               Y3 = YY1
               Z3 = ZZ1
               X4 = XX2
               Y4 = YY2
               Z4 = ZZ2
            CASE(2)
               X3 = XX2
               Y3 = YY2
               Z3 = ZZ2
               X4 = XX3
               Y4 = YY3
               Z4 = ZZ3
            CASE(3)
               X3 = XX3
               Y3 = YY3
               Z3 = ZZ3
               X4 = XX4
               Y4 = YY4
               Z4 = ZZ4
            CASE(4)
               X3 = XX4
               Y3 = YY4
               Z3 = ZZ4
               X4 = XX1
               Y4 = YY1
               Z4 = ZZ1
            END SELECT

            DO IEDGE = 1,12

               SELECT CASE(IEDGE)
                  CASE(1)
                     X1 = XC1
                     Y1 = YC1
                     Z1 = ZC1
                     X2 = XC2
                     Y2 = YC2
                     Z2 = ZC2
                  CASE(2)
                     X1 = XC3
                     Y1 = YC3
                     Z1 = ZC3
                     X2 = XC4
                     Y2 = YC4
                     Z2 = ZC4
                  CASE(3)
                     X1 = XC5
                     Y1 = YC5
                     Z1 = ZC5
                     X2 = XC6
                     Y2 = YC6
                     Z2 = ZC6
                  CASE(4)
                     X1 = XC7
                     Y1 = YC7
                     Z1 = ZC7
                     X2 = XC8
                     Y2 = YC8
                     Z2 = ZC8
                  CASE(5)
                     X1 = XC1
                     Y1 = YC1
                     Z1 = ZC1
                     X2 = XC3
                     Y2 = YC3
                     Z2 = ZC3
                  CASE(6)
                     X1 = XC2
                     Y1 = YC2
                     Z1 = ZC2
                     X2 = XC4
                     Y2 = YC4
                     Z2 = ZC4
                  CASE(7)
                     X1 = XC5
                     Y1 = YC5
                     Z1 = ZC5
                     X2 = XC7
                     Y2 = YC7
                     Z2 = ZC7
                  CASE(8)
                     X1 = XC6
                     Y1 = YC6
                     Z1 = ZC6
                     X2 = XC8
                     Y2 = YC8
                     Z2 = ZC8
                  CASE(9)
                     X1 = XC1
                     Y1 = YC1
                     Z1 = ZC1
                     X2 = XC5
                     Y2 = YC5
                     Z2 = ZC5
                  CASE(10)
                     X1 = XC2
                     Y1 = YC2
                     Z1 = ZC2
                     X2 = XC6
                     Y2 = YC6
                     Z2 = ZC6
                  CASE(11)
                     X1 = XC3
                     Y1 = YC3
                     Z1 = ZC3
                     X2 = XC7
                     Y2 = YC7
                     Z2 = ZC7
                  CASE(12)
                     X1 = XC4
                     Y1 = YC4
                     Z1 = ZC4
                     X2 = XC8
                     Y2 = YC8
                     Z2 = ZC8
               END SELECT

                  CALL IPOINT
     &           (XI,YI,ZI,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     &            X4,Y4,Z4,X5,Y5,Z5,TL,ICODE)
                  IF(ICODE == 0) GOTO 100
                  
            ENDDO

         ENDDO

 100     CONTINUE
      
      RETURN
      END SUBROUTINE CHECK_VOXEL2
C     
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE AUTOCON(KCON,KNBCS) 

C ... Find the connectivity boundary condition patches automatically.

      USE MAIN_ARRAYS, ONLY : CENAX, CENAY, CENAZ, NCHIMT
      USE NS3CO, ONLY : IC9, AUTOCONL
      
      IMPLICIT NONE

      INTEGER :: NB, IB, JB, KNBCS, NNODE, IL, NI, NJ, NK, IERRCODE
      INTEGER :: I, J, K, KK, KP, LP, IFACE, IEKA
      INTEGER, DIMENSION(IC9,KNBCS) :: KCON 

      REAL :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, CC
      REAL, DIMENSION(3) :: A, B


      REAL, ALLOCATABLE :: X(:), Y(:), Z(:)
      REAL, ALLOCATABLE :: XC1(:), YC1(:), ZC1(:), AREAP(:)
      INTEGER, ALLOCATABLE :: IDMN(:),  JDMN(:),  KDMN(:) 
      INTEGER, ALLOCATABLE :: NBCSF(:,:)
      LOGICAL :: SAMESIZE, SAMENCHIMT, FOUND, CGNEAR, CONTYPE


C ... Read the grid file. Master process reads temporarily the whole grid. 

      READ(41) NB

      ALLOCATE(IDMN(NB),JDMN(NB),KDMN(NB),STAT=IERRCODE)
      ALLOCATE(NBCSF(6,NB),STAT=IERRCODE)

      READ(41) (IDMN(IB),JDMN(IB),KDMN(IB), IB=1,NB)

      NNODE = 0

      DO IB=1,NB
         NNODE = NNODE + IDMN(IB)*JDMN(IB)*KDMN(IB)
      ENDDO

      ALLOCATE(X(NNODE),Y(NNODE),Z(NNODE),STAT=IERRCODE)

      IL = 0

      DO IB=1,NB    
         NI = IDMN(IB)
         NJ = JDMN(IB)
         NK = KDMN(IB)
         READ(41)   
     &   (((X(IL+NI*NJ*(K-1)+NI*(J-1)+I),I=1,NI),J=1,NJ),K=1,NK),
     &   (((Y(IL+NI*NJ*(K-1)+NI*(J-1)+I),I=1,NI),J=1,NJ),K=1,NK),
     &   (((Z(IL+NI*NJ*(K-1)+NI*(J-1)+I),I=1,NI),J=1,NJ),K=1,NK)
         IL = IL + IDMN(IB)*JDMN(IB)*KDMN(IB)
      ENDDO

      REWIND 41

      ALLOCATE(XC1(KNBCS),YC1(KNBCS),ZC1(KNBCS),STAT=IERRCODE)
      ALLOCATE(AREAP(KNBCS),STAT=IERRCODE)


C ... Compute the number of BC patches on each block face

      DO IB=1,NB
         DO IFACE=1,6
            K = 0
            DO KK=1,KNBCS
               IF(KCON(2,KK) == IB .AND. KCON(3,KK) == IFACE) THEN 
                  K = K + 1
               ENDIF  
            ENDDO
            NBCSF(IFACE,IB) = K
         ENDDO
      ENDDO


C ... Find the patch center point coordinates

      DO KP=1,KNBCS

         IB = KCON(2,KP)
         IEKA = 1

         DO JB=1,IB-1
            IEKA = IEKA + IDMN(JB)*JDMN(JB)*KDMN(JB)
         ENDDO

         CALL CENTER(X(IEKA),Y(IEKA),Z(IEKA),
     &               XC1(KP),YC1(KP),ZC1(KP),AREAP(KP),
     &   KCON(3,KP),KCON(4,KP),KCON(5,KP),KCON(6,KP),KCON(7,KP),
     &   IDMN(IB),JDMN(IB),KDMN(IB))

      ENDDO


C ... Update the connectivities other than PER

      DO KP=1,KNBCS

         CONTYPE = (KCON(1,KP) ==  1  .AND. KCON(20,KP) /= 1) .OR. ! not PER
     &              KCON(1,KP) ==  6  .OR.  KCON(1,KP) == 11  .OR.
     &              KCON(1,KP) == 14  .OR.  KCON(1,KP) == 15

         IF(.NOT.AUTOCONL) CONTYPE = CONTYPE .AND. KCON(8,KP) == 0

         IF(CONTYPE) THEN
         
            FOUND = .FALSE.

            DO LP=1,KNBCS

              IF(LP == KP) CYCLE

              IF((KCON(1,LP) ==  1  .AND. KCON(20,LP) /= 1) .OR. 
     &            KCON(1,LP) ==  6  .OR.
     &            KCON(1,LP) == 11  .OR.
     &            KCON(1,LP) == 14  .OR. 
     &            KCON(1,LP) == 15) THEN

                  X1 = XC1(KP)
                  Y1 = YC1(KP)
                  Z1 = ZC1(KP)

                  X2 = XC1(LP)
                  Y2 = YC1(LP)
                  Z2 = ZC1(LP)

                  IF(KCON(1,KP) == 6 .AND. KCON(1,LP) == 6) THEN  ! CYC
 
                     X3 = 0.5*(X1+X2)
                     Y3 = 0.5*(Y1+Y2)
                     Z3 = 0.5*(Z1+Z2)

                     A(1) = X2 - X1
                     A(2) = Y2 - Y1
                     A(3) = Z2 - Z1

                     B(1) = X3 - CENAX(KCON(2,LP))
                     B(2) = Y3 - CENAY(KCON(2,LP))
                     B(3) = Z3 - CENAZ(KCON(2,LP))

                     CC = ABS(DOT_PRODUCT(A,B))

                  ELSE

                     CC = SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)

                  ENDIF

C ... The distance between CON patch CGs must be less than 0.5 % of the patch
C ... edge length and the patch areas must not deviate more than 1 %. 

*                  CGNEAR   = CC < 0.01*SQRT(AREAP(KP)) 
*                  SAMESIZE = (ABS((AREAP(KP)-AREAP(LP))/AREAP(KP))<0.02)
                 CGNEAR   = CC < 0.005*SQRT(AREAP(KP)) 
                 SAMESIZE = (ABS((AREAP(KP)-AREAP(LP))/AREAP(KP))<0.01)
               SAMENCHIMT = NCHIMT(KCON(2,LP)) == NCHIMT(KCON(2,KP))

                  IF(CGNEAR .AND. SAMESIZE .AND. SAMENCHIMT) THEN

                     K = 0
                     DO IB=1,KCON(2,LP)-1
                        DO IFACE=1,6
                           K = K + NBCSF(IFACE,IB)
                        ENDDO
                     ENDDO

                     DO WHILE (KCON(3,K+1) /= KCON(3,LP))
                        K = K + 1
                     ENDDO

                     FOUND = .TRUE.
                     KCON( 8,KP) = KCON(2,LP)
                     KCON( 9,KP) = KCON(3,LP)
                     KCON(10,KP) = LP - K
                     KCON(28,KP) = KCON( 8,KP)  ! Global connective block # 

                  ENDIF
             
               ENDIF

            ENDDO

            IF(.NOT.FOUND) THEN
               WRITE(*,111) KCON(2,KP), KCON(3,KP)
               STOP
            ENDIF

         ENDIF

      ENDDO


      DEALLOCATE(IDMN,JDMN,KDMN)
      DEALLOCATE(X,Y,Z)
      DEALLOCATE(XC1,YC1,ZC1,AREAP)
      DEALLOCATE(NBCSF)

 111  FORMAT('AUTOCON failed: Block ',I4,',    Face ',I1)

      RETURN
      END SUBROUTINE AUTOCON
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE CENTER(XCO,YCO,ZCO,XC1,YC1,ZC1,ARE,
     &                  IFACE,IXLO,IXUQ,IYLO,IYUQ,IMAX,JMAX,KMAX)

C ... Calculate the center of gravity for each BC patch. 

      IMPLICIT NONE


      REAL, DIMENSION(*) :: XCO, YCO, ZCO

      REAL :: XC1, YC1, ZC1, XCP, YCP, ZCP
      REAL :: AX, AY, AZ, BX, BY, BZ, AR, ARE

      INTEGER :: IFACE, IXLO, IXUP, IYLO, IYUP, IXUQ, IYUQ
      INTEGER :: IMAX, JMAX, KMAX, ISTRID, JSTRID, KSTRID, ILE
      INTEGER :: IDIM, JDIM, KDIM, IC1, IC2, I, J, IST, JST


      IXUP    = IXUQ + 1
      IYUP    = IYUQ + 1
      ISTRID  = IMAX
      JSTRID  = JMAX
      KSTRID  = KMAX
      ILE     = ISTRID*JSTRID
      IDIM    = IXUP - IXLO
      JDIM    = IYUP - IYLO

      IF(IFACE == 1)      THEN
         IC1 =    1 + (IXLO-1)*ISTRID + (IYLO-1)*ILE
         IC2 =    1 + (IXLO-1)*ISTRID + (IYUP-1)*ILE
         IST = ISTRID
         JST = ILE
      ELSE IF(IFACE == 2) THEN
         IC1 = IXLO                   + (IYLO-1)*ILE
         IC2 = IXLO                   + (IYUP-1)*ILE
         IST = 1
         JST = ILE
      ELSE IF(IFACE == 3) THEN
         IC1 = IXLO + (IYLO-1)*ISTRID
         IC2 = IXLO + (IYUP-1)*ISTRID
         IST = 1
         JST = ISTRID
      ELSE IF(IFACE == 4) THEN
         IC1 = IMAX + (IXLO-1)*ISTRID + (IYLO-1)*ILE
         IC2 = IMAX + (IXLO-1)*ISTRID + (IYUP-1)*ILE
         IST = ISTRID
         JST = ILE
      ELSE IF(IFACE == 5) THEN
         IC1 = IXLO + (JMAX-1)*ISTRID + (IYLO-1)*ILE
         IC2 = IXLO + (JMAX-1)*ISTRID + (IYUP-1)*ILE
         IST = 1
         JST = ILE
      ELSE IF(IFACE == 6) THEN
         IC1 = IXLO + (IYLO-1)*ISTRID + (KMAX-1)*ILE
         IC2 = IXLO + (IYUP-1)*ISTRID + (KMAX-1)*ILE
         IST = 1
         JST = ISTRID
      ENDIF

      XC1 = 0.0
      YC1 = 0.0
      ZC1 = 0.0
      ARE = 0.0

      DO J=IC1,IC2-JST,JST
        DO I=J,J+(IDIM+1)*IST-1-IST,IST

          AX = XCO(I+JST+IST) - XCO(I)
          AY = YCO(I+JST+IST) - YCO(I) 
          AZ = ZCO(I+JST+IST) - ZCO(I) 

          BX = XCO(I+IST) - XCO(I+JST)
          BY = YCO(I+IST) - YCO(I+JST)
          BZ = ZCO(I+IST) - ZCO(I+JST)


          AR = 0.5*SQRT((AY*BZ-AZ*BY)**2
     &                + (AZ*BX-AX*BZ)**2
     &                + (AX*BY-AY*BX)**2)

          ARE = ARE + AR

          XCP = 0.25*(XCO(I) + XCO(I+JST) 
     &              + XCO(I+JST+IST) + XCO(I+IST))
          YCP = 0.25*(YCO(I) + YCO(I+JST) 
     &              + YCO(I+JST+IST) + YCO(I+IST))
          ZCP = 0.25*(ZCO(I) + ZCO(I+JST) 
     &              + ZCO(I+JST+IST) + ZCO(I+IST))
         
          XC1 = XC1 + XCP*AR
          YC1 = YC1 + YCP*AR 
          ZC1 = ZC1 + ZCP*AR 

        ENDDO
      ENDDO

      XC1 = XC1/ARE
      YC1 = YC1/ARE
      ZC1 = ZC1/ARE

      RETURN
      END SUBROUTINE CENTER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE PATCHE(ICON,L,M,IA1,IM1,JA1,JM1,IA2,IM2,JA2,JM2)

C ... Patches are extended outside block boundaries for connections.
C ... This subroutine handles only one patch and its counterpart (CON)
C ... at a time.
      
      USE INTEGERS, ONLY    : IPRO, MBPRO, NPRO
      USE MAIN_ARRAYS, ONLY : NPROCE, IMAXG, JMAXG, KMAXG
      USE NS3CO, ONLY : IC9
      
      IMPLICIT NONE

      INTEGER, DIMENSION(IC9,*) :: ICON

      INTEGER :: M, NDIV, L, II
      INTEGER :: ITY, IFACE1, IFACE2, IBLO1, IMAX, JMAX, KMAX
      INTEGER :: IORI, NSI1, NSI2, IBLO2, IMAX2, JMAX2, KMAX2
      INTEGER :: IXUP1, IYUP1, IXUP2, IYUP2, IXR, IYR
      INTEGER :: IXL, IXU, IXLL, IXUL, IXLE, IXUE
      INTEGER :: IA1, IM1, JA1, JM1, IA2, IM2, JA2, JM2  
      INTEGER :: IYL, IYU, IYLL, IYUL, IYLE, IYUE  
      LOGICAL :: PROB


C ... CON=1, MIR=4, CYC=6, SNG=7, SLD=11, CHI=12, FRE=13, MIX=14, HTS=15

      NDIV = 2**(M-1)  
      ITY  = ICON(1,L)

      IF(ITY ==  1 .OR. ITY ==  4 .OR. ITY ==  6 .OR. ITY ==  7 .OR. 
     +   ITY == 11 .OR. ITY == 12 .OR. ITY == 13 .OR. ITY == 14 .OR.
     +   ITY == 15) THEN

      IFACE1 = ICON(3,L)
      IFACE2 = ICON(9,L)

      IBLO1 = ICON(24,L)
      IMAX  = IMAXG(IBLO1)/NDIV
      JMAX  = JMAXG(IBLO1)/NDIV
      KMAX  = MAX0(1,KMAXG(IBLO1)/NDIV)

      IXUP1 = 0
      IYUP1 = 0
      IXUP2 = 0
      IYUP2 = 0      

C ... Connectivity block is needed (CON, CYC, SLD, MIX and HTS)

      IF(ITY ==  1 .OR. ITY == 6 .OR. ITY == 11 .OR.ITY == 14 .OR.
     +   ITY == 15) THEN
         IBLO2 = NPROCE(1+ICON(8,L),ICON(17,L))
         IMAX2 = IMAXG(IBLO2)/NDIV
         JMAX2 = JMAXG(IBLO2)/NDIV
         KMAX2 = MAX0(1,KMAXG(IBLO2)/NDIV)
      ELSE
         IBLO2 = IBLO1
         IMAX2 = IMAX
         JMAX2 = JMAX
         KMAX2 = KMAX
      ENDIF


      IORI = ICON(14,L)

       IF(IORI == 0) THEN
        IXR = 0
        IYR = 0
      ELSEIF(IORI == 1) THEN
        IF(IFACE1 == 1 .OR. IFACE1 == 3 .OR. IFACE1 == 5) THEN
          IXR = 0
          IYR = 1
        ELSE
          IXR = 1
          IYR = 0
        ENDIF
      ELSEIF(IORI == 2) THEN
        IXR = 1
        IYR = 1
      ELSEIF(IORI == 3) THEN
        IF(IFACE1 == 1 .OR. IFACE1 == 3 .OR. IFACE1 == 5) THEN
          IXR = 1
          IYR = 0
        ELSE
          IXR = 0
          IYR = 1
        ENDIF
      ENDIF

      
      IF(IFACE1 == 1 .OR. IFACE1 == 4) THEN
         IXUP1 = JMAX
         IYUP1 = KMAX
      ENDIF
      IF(IFACE2 == 1 .OR. IFACE2 == 4) THEN
         IXUP2 = JMAX2
         IYUP2 = KMAX2
      ENDIF

      IF(IFACE1 == 2 .OR. IFACE1 == 5) THEN
         IXUP1 = IMAX
         IYUP1 = KMAX
      ENDIF
      IF(IFACE2 == 2 .OR. IFACE2 == 5) THEN
         IXUP2 = IMAX2
         IYUP2 = KMAX2
      ENDIF

      IF(IFACE1 == 3 .OR. IFACE1 == 6) THEN
         IXUP1 = IMAX
         IYUP1 = JMAX
      ENDIF
      IF(IFACE2 == 3 .OR. IFACE2 == 6) THEN
         IXUP2 = IMAX2
         IYUP2 = JMAX2
      ENDIF


      IF(MOD(IFACE1+IFACE2+IORI,2) == 0) THEN

        IF(IXR == 0) THEN
          IXL  =  12
          IXU  =  13
          IXLL =  1
          IXUL =  IYUP2
          IXLE = -1
          IXUE = +1
        ELSE
          IXL  =  13
          IXU  =  12
          IXLL =  IYUP2
          IXUL =  1
          IXLE = +1
          IXUE = -1
        ENDIF
        IF(IYR == 0) THEN
          IYL  =  10
          IYU  =  11
          IYLL =  1
          IYUL =  IXUP2
          IYLE = -1
          IYUE = +1
        ELSE
          IYL  =  11
          IYU  =  10
          IYLL =  IXUP2
          IYUL =  1
          IYLE = +1
          IYUE = -1
        ENDIF
      ELSE
        IF(IXR == 0) THEN
          IXL  =  10
          IXU  =  11
          IXLL =  1
          IXUL =  IXUP2
          IXLE = -1
          IXUE = +1
        ELSE
          IXL  =  11
          IXU  =  10
          IXLL =  IXUP2
          IXUL =  1
          IXLE = +1
          IXUE = -1
        ENDIF
        IF(IYR == 0) THEN
          IYL  =  12
          IYU  =  13
          IYLL =  1
          IYUL =  IYUP2
          IYLE = -1
          IYUE = +1
        ELSE
          IYL  =  13
          IYU  =  12
          IYLL =  IYUP2
          IYUL =  1
          IYLE = +1
          IYUE = -1
        ENDIF
      ENDIF

      IF(IBLO2 /= 0 .AND. ITY ==  1 .OR. ITY == 4  .OR. ITY ==  6 .OR.
     +     ITY == 11 .OR. ITY == 12 .OR. ITY == 14 .OR. ITY == 15) THEN
C ... Used to be without HTS (ITY == 15) because of FIELDVIEW.
C ... Changed Feb 15, 2015 by Esa.         

C ... Patch borders without extension.

        IA1 = ICON(  4,L)
        IA2 = ICON(IXL,L)
        IM1 = ICON(  5,L)
        IM2 = ICON(IXU,L)
        JA1 = ICON(  6,L)
        JA2 = ICON(IYL,L)
        JM1 = ICON(  7,L)
        JM2 = ICON(IYU,L)
           
C ... Extend patch and its counterpart (CON).
         
        IF(ICON(4,L) == 1     .AND. ICON(IXL,L) == IXLL) THEN
           IA1 = ICON(  4,L) - 1
           IA2 = ICON(IXL,L) + IXLE
        ENDIF

        IF(ICON(5,L) == IXUP1 .AND. ICON(IXU,L) == IXUL) THEN
           IM1 = ICON(  5,L) + 1
           IM2 = ICON(IXU,L) + IXUE
        ENDIF

        IF(ICON(6,L) == 1     .AND. ICON(IYL,L) == IYLL) THEN
           JA1 = ICON(  6,L) - 1
           JA2 = ICON(IYL,L) + IYLE
        ENDIF

        IF(ICON(7,L) == IYUP1 .AND. ICON(IYU,L) == IYUL) THEN
           JM1 = ICON(  7,L) + 1
           JM2 = ICON(IYU,L) + IYUE
        ENDIF

      ENDIF

      IF(IA1 > IM1) CALL SWAP_INDEXES(IA1,IM1)
      IF(JA1 > JM1) CALL SWAP_INDEXES(JA1,JM1)
      IF(IA2 > IM2) CALL SWAP_INDEXES(IA2,IM2)
      IF(JA2 > JM2) CALL SWAP_INDEXES(JA2,JM2)

      IF(ITY == 7) THEN                                  ! Singularity
         IF(ICON(4,L) ==     1) IA1 = ICON( 4,L) - 1
         IF(ICON(5,L) == IXUP1) IM1 = ICON( 5,L) + 1
         IF(ICON(6,L) ==     1) JA1 = ICON( 6,L) - 1
         IF(ICON(7,L) == IYUP1) JM1 = ICON( 7,L) + 1
      ENDIF
c      write(*,*) M,IA1,IM1,JA1,JM1
      ENDIF
      
      RETURN
      END SUBROUTINE PATCHE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      SUBROUTINE SWAP_INDEXES(IA,IM)

      IMPLICIT NONE

      INTEGER :: IA, IM, IS

      IS = IM
      IM = IA
      IA = IS

      RETURN
      END SUBROUTINE SWAP_INDEXES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C


