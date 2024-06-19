C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
      SUBROUTINE SECOND(RTIME)
C
C ... MONITOR CPU-TIME ON THE IRIS
C
      
      REAL, DIMENSION (2) :: TARRAY
      REAL :: APU, RTIME
      
      APU = RTIME
      TARRAY(1) = 313.
      RTIME = DPTIME(TARRAY)

      IF(APU < 1.E-10) RTIME = 0.

      RETURN
      END SUBROUTINE SECOND
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C           
! DTIME in standard Fortran.
      Real Function dptime(time)
         Real time(2)
         Double Precision,Save :: last_time = 0
         Double Precision this_time
         Intrinsic Cpu_Time
         Call Cpu_Time(this_time)
         time(1) = this_time - last_time
         time(2) = 0
         dtime = time(1)
         last_time = this_time
      End Function dptime
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C           
      REAL FUNCTION CDELTA(I,J)
C     
C ... CRONECKER'S DELTA
C

      IMPLICIT NONE

      INTEGER :: I, J

      IF(I == J) THEN
         CDELTA = 1.
      ELSE
         CDELTA = 0.
      ENDIF

      RETURN
      END  FUNCTION CDELTA
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C           
      SUBROUTINE WALL_CLOCK_TIME(VALUES,IMONTH,IPAIVA,IHOUR,IMINUT,
     + ISECON)

C ... Print out the wall clock time on the LOG-file

      IMPLICIT NONE

      INTEGER VALUES(*),IMONTH,IPAIVA,IHOUR,IMINUT,ISECON,IWALL(5),
     + MUISTI,MONTH(12),ILPO

      MONTH(1) = 31
      MONTH(2) = 28
      IF(VALUES(1) == 2008 .OR. VALUES(1) == 2012 .OR.
     +   VALUES(1) == 2016 .OR. VALUES(1) == 2020 .OR.
     +   VALUES(1) == 2024 .OR. VALUES(1) == 2028 .OR.
     +   VALUES(1) == 2032 .OR. VALUES(1) == 2036) MONTH(2) = 29   ! Mersu
      MONTH(3) = 31
      MONTH(4) = 30
      MONTH(5) = 31
      MONTH(6) = 30
      MONTH(7) = 31
      MONTH(8) = 31
      MONTH(9) = 30
      MONTH(10)= 31
      MONTH(11)= 30
      MONTH(12)= 31

C ... Seconds

      IWALL(1) = VALUES(7) - ISECON
      IF(IWALL(1) >= 0.) THEN
         MUISTI = 0
      ELSE
         MUISTI   = -1
         IWALL(1) = IWALL(1) + 60
      ENDIF

C ... Minutes

      IWALL(2) = VALUES(6) - IMINUT + MUISTI
      IF(IWALL(2) >= 0.) THEN
         MUISTI = 0
      ELSE
         MUISTI   = -1
         IWALL(2) = IWALL(2) + 60
      ENDIF

C ... Hours

      IWALL(3) = VALUES(5) - IHOUR + MUISTI
      IF(IWALL(3) >= 0.) THEN
         MUISTI = 0
      ELSE
         MUISTI   = -1
         IWALL(3) = IWALL(3) + 24
      ENDIF

C ... Days

      IWALL(4) = VALUES(3) - IPAIVA + MUISTI

      IF(IMONTH /= VALUES(2)) THEN
         ILPO     = IMONTH 
         IWALL(4) = IWALL(4) + MONTH(ILPO)
         ILPO     = ILPO + 1
         IF(ILPO > 12) ILPO = 1
         IF(VALUES(2) /= ILPO) THEN
           IWALL(4) = IWALL(4) + MONTH(ILPO)
           ILPO     = ILPO + 1
           IF(ILPO > 12) ILPO = 1
           IF(VALUES(2) /= ILPO) THEN
             IWALL(4) = IWALL(4) + MONTH(ILPO)
             ILPO     = ILPO + 1
             IF(ILPO > 12) ILPO = 1
             IF(VALUES(2) /= ILPO) THEN
               IWALL(4) = IWALL(4) + MONTH(ILPO)
               ILPO     = ILPO + 1
               IF(ILPO > 12) ILPO = 1
               IF(VALUES(2) /= ILPO) THEN
                 IWALL(4) = IWALL(4) + MONTH(ILPO)
               ELSE
                 GO TO 1000
               ENDIF
             ELSE
               GO TO 1000
             ENDIF
           ELSE
             GO TO 1000
           ENDIF
         ELSE
           GO TO 1000
         ENDIF         
      ELSE
         CONTINUE
      ENDIF

1000  CONTINUE ! Mersu

      IF(IWALL(4) >= 100) THEN
      WRITE(13,'(A,I3,A,I2,A,I2,A,I2,A)') '  Wall-clock time is ',
     + IWALL(4),' days ',IWALL(3),' hours ',IWALL(2),' minutes ',
     + IWALL(1),' seconds'

      ELSE IF(IWALL(4) > 0 .AND. IWALL(4) < 100) THEN
      WRITE(13,'(A,I2,A,I2,A,I2,A,I2,A)') '  Wall-clock time is ',
     + IWALL(4),' days ',IWALL(3),' hours ',IWALL(2),' minutes and ',
     + IWALL(1),' seconds'

      ELSE IF(IWALL(3) > 0) THEN
      WRITE(13,'(A,I2,A,I2,A,I2,A)') '  Wall-clock time is ',
     + IWALL(3),' hours ',IWALL(2),' minutes and ',IWALL(1),' seconds'

      ELSE IF(IWALL(2) > 0) THEN
      WRITE(13,'(A,I2,A,I2,A)') '  Wall-clock time is ',
     + IWALL(2),' minutes and ',IWALL(1),' seconds'

      ELSE
      WRITE(13,'(A,I2,A)') '  Wall-clock time is ',
     + IWALL(1),' seconds'
      ENDIF

      RETURN
      END SUBROUTINE WALL_CLOCK_TIME
C     
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C     
