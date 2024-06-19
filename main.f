      PROGRAM FINFLO
C     **************
C     +----------------------------------------------------------------+
C     |                                                                |
C     |    A LU-factorized finite-volume based 3D-Navier-Stokes code   |
C     |    applying multi-block grids, multigrid cycling and           |
C     |    local time-stepping.                                        |
C     |                                                                |
C     |    Original code developed by Dr. Timo Siikonen                |
C     |    at Helsinki University of Technology,                       |
C     |    Laboratory of Aerodynamics.                                 |
C     |                                                                |
C     |    Reynolds Stress model and scalar solver developed by        |
C     |    Patrik Rautaheimo.                                          |
C     |                                                                |
C     |    Patch-oriented connectivity and debugging by Esa Salminen.  |
C     |                                                                |
C     |    Dynamic memory allocation, new boundary routines, dynamic   |
C     |    grid generation and parallelization by Petri Kaurinkoski.   |
C     |                                                                |
C     |    Copyright (C) 1988-1995 AeRoLa     Version 1.0              |
C     |                                                                |
C     |    MPI-standart and new turbulence models (Exp. RSM by         |
C     |    Speziale and Cross-diffusion model) included by             |
C     |    Patrik Rautaheimo.                                          |
C     |                                                                |
C     |    Chimera and sliding mesh techniques included by Esa         |
C     |    Salminen and Timo Siikonen.                                 |
C     |                                                                |
C     |    Pseudo compressibility and time integration included by     |
C     |    Timo Siikonen.                                              |
C     |                                                                |
C     |    Version 2.0  January, 1996                                  |
C     |    Copyright (C) 1995-1996 Thermo/HUT                          |
C     |                                                                |
C     |    Version 2.1  January, 1996                                  |
C     |    Copyright (C) 1995-1996 Thermo/HUT                          |
C     |    Flexible solid walls.                                       |
C     |                                                                |
C     |    Version 3.0 November, 1997 FINFLO2000                       |
C     |    Copyright (C) 1995-2001 Thermo/HUT                          |
C     |    Non-matching boundaries by Esa Salminen and Patrik          |
C     |    Rautaheimo 20.11.97                                         |
C     |                                                                |
C     |    Grid in DOUBLE PRECISION                                    |
C     |    Rautaheimo 28.1.98                                          |
C     |                                                                |
C     |    Block and BC data is allocated dynamically. Also main.f     |
C     |    and ns3c.f are using f90-features.                          |
C     |    Rautaheimo 10.8.98                                          |
C     |                                                                |
C     |    Working array ZZZ is allocated dynamically.                 |
C     |    Rautaheimo 3.11.98                                          |
C     |                                                                |
C     |    Copyright (C) 2001- 2016 Finflo Ltd.                        |
C     |                                                                |
C     |    Two-phase model added   1.10.2003                           |
C     |                                                                |
C     |    Version 7.3 (2.2.2006)                                      |
C     |                                                                |
C     |    Version 8.0 (1.4.2006) includes multi-phase and             |
C     |    free-surface models                                         |
C     |                                                                |
C     |    Version 8.1 (19.6.2008) includes simple cavitation model    |
C     |    and trajectory calculation (6DOF)                           |
C     |                                                                |
C     |    Finflo-8.2.0 started 11.3.2009                              |
C     |                                                                |
C     +----------------------------------------------------------------+
C
C     Memory is allocated dynamically for the following variables:
C     ------------------------------------------------------------
C     SIZE MAXB:
C        A1,A2,A3,VOL,
C        A1XA,A1YA,A1ZA,
C        A2XA,A2YA,A2ZA,
C        A3XA,A3YA,A3ZA,
C        D1,D2,D3,XFC,YFC,ZFC,
C        UROT,VROT,WROT,DISTW,DISTP,LOCDIS,IDP,RLD,RLDIST
C        UROTCP,VROTCP,WROTCP
C           SUBTOTAL 31  VARIABLES
C
C        RO,RM,RN,RW,E,
C        DRO,DM,DN,DW,DE,
C        TOLD,UOLD,VOLD,WOLD,POLD,
C        OHMI,TIJ,U,V,W,P,PDIFF,EPS2,VIST,C,TEMP,
C        F1R,F1RM,F1RN,F1RW,F1E,VIS,DTL,CP,CH,DRDP,DRDH,
C        HAT1,HAT2,HAT3,HAT4,F1RK,F1EPS,RKSI,BLANK,ROLD 
C        RUVAV,RUWAV,RVWAV,                       (cross terms of Reynolds
C        JET                                       stresses)
C           SUBTOTAL 50  VARIABLES
C
C        STRAIN,SHEAR(6),VORT(3)
C           SUBTOTAl 10 VARIABLES
C
C           ---   TOTAL 83  VARIABLES OF SIZE MAXB  ---  (Minimum required)
C
C        RK,REPS,DDEPS,DRK,DEPS,RKOLD,EPSOLD      (k-epsilon)
C        SRK,SEPS,PTUR,FUN1,TTS,VTRAN,RNUT,
C        QSAS,VELLAP
C           SUBTOTAL 16 VARIABLES
C     
C        ALFA,X,DX,XOLD,DE,DTEMP,TOLD,EVAPR,      (Two-phase model)
C        F1R,F2R,F3R,FRO,FE,FA,F1A,F2A,F3A
C           SUBTOTAL 17*NPHASES VARIABLES
C     
C        TEMP,RO,E,DRODP,DRODH,VIS,CH,C,QIF,      (Properties for the two-phase
C        HSAT,DHSDP,DHSDH,CP,TAUF,BUBBLES,alpo             model)
C        TSAT,PSAT
C           SUBTOTAL 15*NPHASES + 2 VARIABLES
C
C        F1R,F2R,F3R                              (Pressure correction)
C           SUBTOTAL 3 VARIABLES
C
C        S11(6)                                   (shear stresses, optional)
C           SUBTOTAL 6 VARIABLES                  (STRESL = .true. or .false.)
C
C        BIJ(6),WIR(3)                            (EARSM, optional)
C           SUBTOTAL 9 VARIABLES
C
C        PROD(6),SPI(6),DIF(6),DIS(6),VVIS(6)     (RSM)
C           SUBTOTAL 30 VARIABLES
C
C        FWLL                                     (RSM)
C           SUBTOTAL 1 VARIABLES
C
C        FI(NSCAL),DFI(NSCAL),FIOLD(NSCAL),
C        SFI(NSCAL),F1FI(NSCAL)                   (scalar equations)
C           SUBTOTAL NSCAL*5 VARIABLES
C
C        ROLE2,RMLE2,RNLE2,RWLE2,ELE2             (time accurate)
C        ROLE3,RMLE3,RNLE3,RWLE3,ELE3
C        ROAV1,RMAV1,RNAV1,RWAV1,EAV1,PAV1,TAV1,  (time averaged)
C        ROAV2,RMAV2,RNAV2,RWAV2,EAV2             (average of squared)
C        (RUVAV,RUWAV,RVWAV)                      (cross terms above)
C           SUBTOTAL 22 VARIABLES
C
C        RKLE2,EPSLE2
C        RKLE3,EPSLE3
C        RKAV1,EPSAV1
C        RKAV2,EPSAV2
C           SUBTOTAL 8 VARIABLES
C
C        FILE2(NSCAL)
C        FILE3(NSCAL)
C        FIAV1(NSCAL)
C        FIAV2(NSCAL)
C           SUBTOTAL 4*NSCAL VARIABLES
C
C        ROFOR,RMFOR,RNFOR,RWFOR,EFOR,PDFOR       (Chimera)
C        WGH1,WGH2,WGH3,WGH4
C           SUBTOTAL 10
C
C        RKFOR,REFOR                              (Chimera+k-e)
C           SUBTOTAL 2
C
C        FIFOR(NSCAL)                             (Chimera+scalar)
C           SUBTOTAL NSCAL
C
C        II1,II2,II3,II4,INTF                     (Chimera interpolation) 
C           SUBTOTAL 5
C
C        XXI,YXI,ZXI,XETA,YETA,WFS,XI,ETA,
C        ZETA,XHULL,YHULL,ZHULL,XVL,YVL,ZVL,
C        ZZTOP,UHULL,VHULL,WHULL,XSP,YSP,ZSP
C        F1H,SUROLD
C           SUBTOTAL 22 + 2 VARIABLES
C
C        IHULL,JHULL, IWAVEB
C           SUBTOTAL 3 VARIABLES
C    
C
C           min   k-e   shear   EARSM RSM   scalars   time k-e  scalar  chimera  
C --- TOTAL (73 + 16 + 10 + 6  +  9  + 30+1 + NSCAL*5 + 22 + 8 + 4*NSCAL + 10 +
C           2 +  NSCAL +         5   + 17*NPHASES + 15*NPHASES+2 +3 + 27*MAXFS)
C        chim+k-e chim+scalar chim+int multiphase properties precor  fresur
C               VARIABLES OF SIZE MAXB ---    (Maximum required)
C
C
C     Memory is allocated dynamically for the following variables:
C     ---------------------------------------------------------------------
C        XCO,YCO,ZCO,XC,YC,ZC   
C           SUBTOTAL 6  VARIABLES (double precission)
C
C        XLE2,YLE2,ZLE2                       (moving grid)
C        XLE3,YLE3,ZLE3
C        XGRI,YGRI,ZGRI
C        XORI,YORI,ZORI
C           SUBTOTAL 12 VARIABLES  (double precission)
C
C     ---------------------------------------------------------------------
C
C
C     SIZE MAX11:   (greatest block with multi grid)
C        W12(3),SIJ(6),GRADT(3),GRADK(3),GREPS(3)
C           SUBTOTAL  18 VARIABLES
C
C
C     SIZE MAX2:
C        ROP2H,RMP2H,RNP2H,RWP2H,EP2H,SURH2
C           SUBTOTAL  6 VARIABLES
C
C        RKP2H, EPSP2H                      (k-epsilon)
C           SUBTOTAL  2 VARIABLES
C
C        VAR%ROP2H, VAR%EP2H         (Multiphase variables)
C           SUBTOTAL 2'*NPHASES VARIABLES
C        FIP2H(NSCAL)                (scalar equations)
C           SUBTOTAL  NSCAL VARIABLES
C
C        --- TOTAL  7+NSCAL+2*NPHASES OF SIZE MAX2 ---
C
C     SIZE IB:
C        (TA1BI,TA1BJ,TA1BK,TA1TI,TA1TJ,TA1TK,
C        TA2BI,TA2BJ,TA2BK,TA2TI,TA2TJ,TA2TK,
C        QA1BI,QA1BJ,QA1BK,QA1TI,QA1TJ,QA1TK,
C        QA2BI,QA2BJ,QA2BK,QA2TI,QA2TJ,QA2TK),
C        ICP,JCP,KCP,IJMASK,
C        UBI,VBI,WBI,UBJ,VBJ,WBJ,UBK,VBK,WBK,
C        UTI,VTI,WTI,UTJ,VTJ,WTJ,UTK,VTK,WTK
C
C           ---   TOTAL 22 VARIABLES OF SIZE IB    ---
C
C     SIZE IBF:                (this size is for solid+FRE pathches)
C        HFLUX,TWALL,CPWALL,TAUW1,TAUW2,QWALL,QWFRIC,
C        UUALL,VWALL,WWALL,SURFX,SURFY,SURFZ, ! (UUALL is 'UWALL' outside main)
C        WMFLUX,POROS,RMLOSS,WHSTAG,WTEMP,RSDIRX,RSDIRY,RSDIRZ,
C        RBK,XCP,YCP,ZCP,TAUWX,TAUWY,TAUWZ,SURLE,DSURLE,WH
C
C        BOUNMF,NOUNU,BOUNV,BOUNW,BOUNT,BOUNP,BOUNR,BOUNE,
C        BOUNRK,BOUNEP,BOUNPD,BOUNA1,BOUNA2,BOUNFI(NSCAL),BOUNBI(6),
C        BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2
C
C           ---   TOTAL 31 + 25 + 6 + NSCAL VARIABLES OF SIZE IBF    ---
C
C        WAVES,WAVEH
C           ---   TOTAL 2 VARIABLES (double precision)
C
C     Maximum memory used by the main arrays is
C     ((190+NSCAL*10+32*NPHASES)*MAXB + 18*MAX11 + (8+NSCAL+2*NPHASES)*MAX2 + 
C     22*IB + (62+NSCAL)*IBF) + 27*MAXFS   words
C
C     Maximum memory used by the double precision arrays is
C     (15*MAXB)*2 + (2*IBF)*2 words
C
C
C     ------------------------------------------------------------------
C
C     The number of scalar equations is NSCAL=MAX0(0,IBREY*6+MSCAL),
C     where IBREY is 1 if the RSM is applied, otherwise 0. MSCAL is the
C     number of additional scalar equations.
C
C     ------------------------------------------------------------------
C
C        VARIABLE MENU IN FINFLO
C        =======================
C
C     WHOLE MESH SIZE ARRAYS (MAXB)
C
C     GEOMETRY VARIABLES:
C          A1,A2,A3       - areas of the face surface of the cells
C          VOL            - volume of the cells
C          A1XA,A1YA,A1ZA - unit vectors of surface 1 of the cells
C          A2XA,A2YA,A2ZA - unit vectors of surface 2 of the cells
C          A3XA,A3YA,A3ZA - unit vectors of surface 3 of the cells
C          D1,D2,D3       - height of the cells
C          UROT,VROT,WROT - static moment of the cell face
C          UROTC,VROTC,WROTC - cell center-point grid velocities
C          XC,YC,ZC       - center points of the grid (DBLE)
C          DISTW          - global distance from the wall
C          LOCDIS         - indicator for the distances
C          RLD, RLDIST    - differential equation based wall distance
C          XCO,YCO,ZCO    - corner points of the grid (DBLE)
C          XGRI,YGRI,ZGRI - stationary grid (DBLE, for moving grid)
C          XORI,YORI,ZORI - original grid (DBLE, for moving grid)
C          XFC,YFC,ZFC    - unused
C     
C
C     BASIC FLOW VARIABLES:
C          RO             - density
C          RM,RN,RW       - momentum components
C          E              - total energy
C          DRO            - density residual (some times temperature)
C          DM,DN,DW       - residual of the momentum components 
C                           (some times velocity)
C          DE             - residual of the total energy(s. t. pressure)
C          TOLD           - old temperature
C          UOLD,VOLD,WOLD 
C                         - old velocity components
C          POLD           - old static pressure
C          ROLD           - old density
C          OHMI           - absolute value of the vorticity
C          TIJ            - Lighthill's tensor gradient
C          U,V,W          - velocity components
C          P              - pressure
C          PDIFF          - pressure difference (P-PREF)
C          EPS2           - turbulent over molecylar viscosity plus 1
C          VIST           - dimensional turbulent eddy viscosity
C          C              - speed of sound
C          TEMP           - tempperature
C          F1R            - density flux
C          F1RM,F1RN,F1RW - momentum fluxes
C          F1E            - total energy flux
C          F1H            - wall-height flux
C          VIS            - molecylar viscosity
C          DTL            - local time step
C          CH             - thermal conductivity (=k,W/(mK))
C          DRDP,DRDH      - thermodynamical derivatives
C          SUROLD         - old surface level (for multigrid, unused)
C
C          HAT1,HAT2,HAT3,HAT4
C                         - auxliarly arrays for scalar fluxes
C          RKSI           - limiter function for chimera treatment
C          BLANK          - blanking vector in PLOT3D format
C
C     K-EPSILON VARIABLES:
C          RK             - kinectic energy of turbulence
C          REPS           - dissipation of RK
C          DDEPS          - difference between diss. tilde and diss.
C          DRK,DEPS       - residuals of    the RK and REPS
C          RKOLD,EPSOLD   - old values of   the RK and REPS
C          SRK,SEPS       - source terms of the RK and REPS
C          F1RK,F1EPS     - fluxes of       the RK and REPS
C          PTUR           - production of the turbulent energy
C          FUN1           - mixing function for k-omega equations
C          TTS            - time scale of turbulence
C          VTRAN          - Multiplier for transition
C          RNUT           - Variable in Spalart-Allmaras model
C          QSAS           - SAS source term
C          VELLAP         - velocity Laplacian
C
C     TWO-PHASE MODEL VARIABLES (NPHASES):
C          ALFA           - volume fraction
C          X              - quality
C          DX             - quality residual
C          XOLD           - old residual
C          DE             - energy residual
C          DTEMP          - change in temperature
C          TOLD           - old temperature
C          EVAPR          - mass transfer
C          F1R,F2R,F3R    - mass fluxes
C          FRO,FE         - auxiliary mass flux, energy flux
C          FA,F1A,F2A,F3A - auxiliary volume flux ,volume fluxes
C
C     TWO-PHASE PROPERTIES (NPHASES):
C          TEMP           - temperature
C          RO             - density
C          E              - specific internal energy
C          DRODP          - derivative with respect to pressure
C          DRODH          - derivative with respect to enthalpy
C          VIS            - viscosity
C          CH             - heat conductivity
C          C              - sound speed
C          QIF            - heat transfer rate from the interface to fluid
C          HSAT           - saturation enthalpy
C          DHSDP          - derivative of saturation enthalpy
C          DHSDH          - dummy 
C          CP             - specific heat
C          TAUF           - evaporation time constant
C          N              - number of bubbles or droplets / volume
C          alpo           - testing
C
C     TWO-PHASE PROPERTIES:     
C          TSAT           - saturation temperature
C          PSAT           - saturation pressure
C
C     PRESSURE CORRECTION VARIABLES:
C          F1R,F2R,F3R    - mass fluxes
C
C     OPTINAL VARIABLES:
C          S11(6)         - shear stresses, optional
C
C     EARSM MODEL VARIABLES:
C          BIJ(6)         - anisotropy components b_ij
C          WIR(3)         - for rotational correction
C
C     RSM MODEL VARIABLES:
C          PROD(6)        - prodution of the RS (Reynolds stresses)
C          SPI(6)         - pressure redistribution of the RS
C          DIF(6),VVIS(6) - turbulent and molecylar diffusion of the RS
C          DIS(6)         - dissipation of the RS
C          FWLL           - wall function
C
C     RSM MODEL AND SCALAR VARIABLES
C          FI(NSCAL)      - conservative variable
C          DFI(NSCAL)     - residual  of the FI
C          FIOLD(NSCAL)   - old value of the FI
C          SFI(NSCAL)     - source    of the FI
C          F1FI(NSCAL)    - flux      of the FI
C     If RSM model is applied scalar numer 1,2,3,4,5,6 are respectly
C     uu,uv,uw,vv,vw,ww
C
C     TIME LEVEL VARIABLES (TIMEL == TRUE)
C           ROLE2,RMLE2,..,
C           ROLE3,RMLE3,..,
C           FILE2(NSCAL)  - time level 2 variables
C           FILE3(NSCAL)  - time level 3 variables
C           XLE2,YLE2,ZLE2- corner points of the old grid (DBLE)
C           XLE3,YLE3,ZLE3- corner points of the old grid (DBLE)
C           ROAV1         - time averaged for period TAVER
C           ROAV2         - time averaged for period TAVER for square values
C           RUVAV...      - cross terms of the Reynolds stresses
C
C     CHIMERA VARIABLES (CHIMER == TRUE)
C           ROFOR,RMFOR,RNFOR,RWFOR,EFOR
C                         - fortifite function with chimera
C           RKFOR,REFOR   - fortifite function with chimera + k-e
C           WGH1,WGH2...  -
C           II1,II2...    - 
C           FIFOR(NSCAL)  - fortifite function with chimera + scalar
C
C     SIZE MAX11 VARIABLES
C
C     STRAIN RATES AND VORTICITY RATES
C           W12(3)        - vorticity tensor
C           SIJ(6)        - shear rate tensor
C           GRADT(3)      - temperature gradient
C           GRADK(3)      - gradient of RK
C           GREPS(3)      - gradient of epsilon-like variable
C           STRAIN        - shear strain
C           SHEAR(6)      - components of shear strain
C           VORT(3)       - components of vorticity tensor (2xvector)!
C
C
C     max2 and ib size variables I or someone else can make later
C     ppr 5.5.1995
C
C     Later (19.7.2006):
C
C     SIZE IBF VARIABLES
C           HFLUX         - wall heat flux given in a boundary file
C           QWALL         - calculated wall heat flux
C           QWFRIC        - heat flux produced by friction
C           CPWALL        - pressure coefficient
C           UWALL,VWALL.. - surface velocity
C           TWALL         - surface temperature
C           TAUW1,TAUW2   - shear stresses along the grid lines
C           SURFX,SURFY.. - surface force
C           XCP,YCP,ZCP   - centrepoint coordinates
C           WMFLUX        - calculated mass flux
C           POROS         - porosity
C           RMLOSS        - K-value for a pressure loss
C           WHSTAG        - stagnation enthalpy with injection
C           WTEMP         - unused wall temperature
C           RSDIRX,..     - dircetion of the wall injection
C           RBK           - equivalent surface roughness
C           TAUWX,TAUWY.. - shear stresses on the wall
C           SURLE,DSURLE  - free-surface related variables
C           UTAUM         - Friction velocity of the first cell
C
C           WAVES,WAVEH   - wave height (nodes and faces, DBLE)
C
C           BOUNMF        - mass flux based on a boundary file
C           BOUNU,BOUNV.. - velocities -''-
C           BOUNP         - surface pressure given in a boundary file
C           BOUNT         - surface temperature given in a boundary file
C           BOUNR         - surface density
C           BOUNE         - total energy on the surface
C           BOUNRK,BOUNEP - turbulence variables
C           BOUNPD        - pressure diffeence on the surface
C           BOUNA1,BOUNA2 - volume fractions on the surface
C           BOUNFI(NSCAL) - scalar variables
C           BOUNBI(6)     - Reynolds stresses in EARSM
C
C     SIZE IB VARIABLES
C           UBI,VBI,VBK.. - surface velocity (obsolate, in algebraic turb.)
C           ICP,JCP,KCP   - indeces for algebraic models
C           IJMASK        - mask for solids (unused, obsolate)
C
C ***************************************
C *** EXPLANATION SOME OF THE INDEXES ***  
C ***************************************
C     
C     Chimera valiables:
C        NCHIMT(NBGG) - level of chimera blocks. 0-ordinary 1 or high chim
C        NCHOR(NBGG)  - order of blocks in chimera. Smaller NCHIMT first.
C                       Size NBLOCK. Ordinary blocks are first
C
C     Parallelization varibles
C        NLOCAL(NBGG) - index is global block number and value is local block
C                     number. Value is -1 if block is not in the proces
C        NPNUM(NBGG)  - index is global block number and value is process 
C                     number.
C        TAVER        - time of the averaged is done for time accurate 
C                       simulation
C        NAVER        - number of time steps for averaging
C
C ************************************************
C *** EXPLANATION SOME OF THE LOGICAL VARIABLE ***  
C ************************************************
C     STARTL - write RSTART       STRESL - is strain rate tensor (S11) in use
C     FULLNS - exact viscous term SOURL  - source terms in RSM
C     STATEL - default ideal gas  TIMEL  - time accurate calculation
C     COORL  - moving coordinates PARALLEL  - MPI in use
C     PRESL  - obsolate, old??    REANEW - if EARSM is started from 2.eq
C     GRAVIL - gravitation        CONVL  - follow convergence with MPI
C     SPLIT  - DIVP3D.INFO pres.  
C     CHIMEL - chimera in use     GROUP  - force groups in use
C     PERCHL - periodic bc        MULPHL - two-phase flow enabled
C     FRESUR - free-surface       CAVITL - cavitation model
C     TRUE_DISTL - true wall distances (default)
C     FPRINTL- blocks unnecessary temporary printing
C     INOUTL - contains inlets or outlets
C     TURDESL- DES, currently Spalart-Allmaras or SST K-omega
C     INCHIML- Chimera interpolation first-order in the boundary
C     REFLECL- Instead of wall velocities extrapolation from the domain
C     ENTROPY_FIX - Make the entropy fix in FLUXKE (only for Roe's method)
C     TRANSL - Activate transition model
C ************************************************
C *** EXPLANATION SOME OF THE INTEGER VARIABLE ***  
C ************************************************
C     IROTCO - rotational correction 2.eq
C *****************************
C *** EXPLANATION WORK ARRAY  *
C *****************************
C
C     ZZZ  - is a work array. It is size of the 46*biggest slab (FLUMPH).
C     it also has some other disdrictions. It is dynamically allocated
C     and the size is MAXW. PPR 3.11.98
C
C ********************************
C *** END OF THE VARIABLE MENU ***
C ********************************
C
C     MONIHILARESERVI     = 58 (= MBSTP belove)  This is a piece of history.
C     MONILOHKORESERVI    = 10 (= MGSTP belove)  This is a piece of history.
C     MONIHILATASOMAKSIMI = 6  (= MGM in NS3CO)
C
C     MREQ is the memory requirement in words. On SGI computers the word
C     length is 4 bytes. Therefore the memory requested is MREQ*4 bytes.
C

      USE MPI
      USE CHARACTERS
      USE INTEGERS
      USE MAIN_ARRAYS
      USE TYPE_ARRAYS
      USE CONSTANTS, ONLY : PII, DEG2RAD
      USE NS3CO
      USE GM
      USE FLIGHT
      USE BLADE_VARIABLES, ONLY : TDBL,DTB
*      USE BLADE_VARIABLES, ONLY : RCG,RITH,RIZE,RIBE,RHINGE

      IMPLICIT NONE

      INTEGER, PARAMETER :: MGSTP = 0, MBSTP = 0 ! 100 ! Nice piece of history.

C         PARAMETER (MBITS = 4,NPTIM = 2)     ! SGI
CRAY      PARAMETER (MBITS = 8,NPTIM = 1)     ! CRAY - (sniff, sniff)
C         NPTIM = 2 for GRIDS with Double precission
C         NPTIM = 1 for GRIDS with READ. Option -dp

      INTEGER :: IBCOM(NBGGG), IBCOM2(NBGGG), IBSLB(NBGGG)
      INTEGER :: NPNTS

      INTEGER :: MAX2, MAXTB, MAXT2, MAXRB, NMAXSB,
     2           NMAXS2,NMAXTS,MAXSCH,MAXP2

      INTEGER(KIND=8) :: 
     & IM1,  IE1,  IM11,  IM2,   IM3,  IM4,  IM5,   IM4B, ! IM4C, IM4D, 
     & IM4E, IM4F,  IM4G,  IM4H, IM4I, IM6,   IM7,
     & IM60,IM6A,IM6B,IM6C,IM10,IMD1,IMD1A,IIB,
     & IIB2,IIB3,IPGRI,MREQB,IM12,IM13,IM14,IIB22,
     & IIB4,IIB5,IIB6

      CHARACTER(LEN=80) :: RIVI, FILETYPE, GRIDFI, BCFILE, OUTPUT
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      CHARACTER(LEN=4)  :: RUN
      CHARACTER(LEN=3)  :: CN
      INTEGER :: VALUES(8), IPAIVA, IMONTH, IHOUR, IMINUT, ISECON
      LOGICAL :: LRIVI, NAMEL, MASTER
      INTEGER :: IOLD1, IERR, IDIS0, MAXWNE, MAWNE, NBMAX, NPPN,
     &           I, J, I0, MAXCO, NMAXSC, IERRCODE, NP, IERR3,
     &           ERRORCODE, RC
      REAL    :: OMEGAREF

      NAMELIST /WORKS/ IOLD1, GRIDFI, BCFILE, OUTPUT, OUTVAR

      DATA MAX2,MAXTB,MAXT2,MAXRB,NMAXSB/
     &     0,0,0,0,0/

      MAXB  = 0
      MAXSB = 0
      MAXS2 = 0
      MAXCH = 0   
*      MBITS = 4
      MBITS = 8  ! Double precision version
      NPTIM = 2
        
C ... Check the input file type from the command line

      I0    = IARGC()
      NAMEL = .TRUE.
        
      IF (I0 > 0) THEN
         DO I=1,I0
            CALL GETARG(I,FILETYPE)
            IF(FILETYPE(1:2) == '-n') NAMEL = .TRUE.
*            IF(FILETYPE(1:2) == '-l') NAMEL = .FALSE.
            IF(FILETYPE(1:1) /= '-')  GOTO 10
         ENDDO
      ENDIF
 10   CONTINUE

C ... Find the number of processes (NPRO)
           
      CALL MPRO1(IPRO,PRN,NPRO,IOTY)
         
      MASTER = IPRO == 1

      NPNTS = 1
       
      CALL OPENER(PRN)
        
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

      ISECON  = VALUES(7)
      IMINUT  = VALUES(6)
      IHOUR   = VALUES(5)
      IPAIVA  = VALUES(3)
      IMONTH  = VALUES(2)

      WRITE(13,*)
      IF(VALUES(2) >= 10) THEN
         WRITE(13,'(A,I2,A1,I2,A1,I4)') '  Today is ',
     +   VALUES(3),'.',VALUES(2),'.',VALUES(1)
      ELSE
         WRITE(13,'(A,I2,A1,I1,A1,I4)') '  Today is ',
     +   VALUES(3),'.',VALUES(2),'.',VALUES(1)
      ENDIF
      WRITE(13,'(A,I2,A1,I2,A1,I2,A)') '  It is ',
     + VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'

      IF(PARALLEL)  THEN
         CALL MPICHE(IPRO)
         IF(.NOT.MASTER) THEN
            WRITE(13,'(A,I3)') '  This is process no.',IPRO
         ELSE
            WRITE(13,'(A,I3)') '  This is the main process'
         ENDIF
         WRITE(13,*)
      ENDIF

C ... Restart check

      IF(MASTER) THEN
*         IF(NAMEL) THEN ! Namelist input
            READ(2,NML=WORKS)
*         ELSE ! Traditional input file
*            READ(2,9945) NAME
*            IF(NAME(1:14) == 'FINFLO VERSION') READ(2,*) ! if to be removed
*            READ(2,*) IOLD1
*         ENDIF
         REWIND 2
      ENDIF
*9945  FORMAT(A80)
       
      IF(PARALLEL) CALL MPI_BCAST(IOLD1,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     +     IERR)

      IF(IOLD1 > 0) THEN ! READ RESTART FILE
         IF(MASTER) THEN
            OPEN(20,FILE='RSTART',STATUS='UNKNOWN',FORM='UNFORMATTED')

            READ(20)

     &         TDBL,    DTB,      THETA,   H0,      DROMAX,   
     &         CD,      CL,      CS,      
     &         CX,      CY,      CZ,      CMX,     CMY,     CMZ,   
     &         RMUFRS,  FRSMUT,  TOM,    
     &         FRSVIS,  FRSSIE,  T0,      ANGLE,   
     &         VOLTOT,  QMFIN,   QMEIN,   QMFOUT,  QMEOUT,  RKMAX,  
     &         TAVER,   NBLOCG,  ROTANG(1:NBLOCG),        

     &         ICYCLE,  IPRINT,  IXERR, 
     &         INERR,   IMACH,   NBLOCK,  ICYOLD, 
     &         IUPDAT,  NSCAL,   JSCAL,   IBOUT,   ITIMES,  ICYTOT,
     &         NAVER

         ENDIF
 
         IF(PARALLEL) THEN

            CALL MPI_BCAST(TDBL   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(DTB    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(THETA  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(H0     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(DROMAX ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CD     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CL     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CS     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CX     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CY     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CZ     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CMX    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CMY    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(CMZ    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(RMUFRS ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FRSMUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(TOM    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FRSVIS ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(FRSSIE ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(T0     ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ANGLE  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(VOLTOT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(QMFIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(QMEIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(QMFOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(QMEOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(RKMAX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(TAVER  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(NBLOCG ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ROTANG ,NBLOCG,
     &                               MPI_REAL8,0,MPI_COMM_WORLD,IERR)

            CALL MPI_BCAST(ICYCLE ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(IPRINT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(IXERR  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(INERR  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(IMACH  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(NBLOCK ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ICYOLD ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(IUPDAT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(NSCAL  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(JSCAL  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(IBOUT  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ITIMES ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(ICYTOT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(NAVER  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

         ENDIF
      ENDIF
                   
      IOLD  = IOLD1
      DTOLD = DTB
      MAXW  = 0
      IDIS0 = IDIS  ! old length scale equation
      
      CALL COMSI(MGSTP,MBSTP,IB,MAXB,MAX11,MAX2,MAXWNE,MAWNE,NPNTS,
     +     IBCOM,IBCOM2,IBSLB,NBMAX,IBF,OMEGAREF,NAMEL,NPPN,IOLD1,
     +     GRIDFI,BCFILE,OUTPUT)

*      MAXW = MAX0(MAXW,MAXWNE,NPPN*25)
      MAXW = MAX0(MAXW,MAXWNE,NPPN*IC9)
      MAW  = MAWNE
      
 2201 FORMAT(/9X,'MBPRO (=',I4,') is too small.'/
     &     9X,'Should be at least same as number of blocks (=',I4,')')

      WRITE(45,*)
      DO 100 I=1,NBLOCK
         WRITE(45,1000) 'IG = ',(IG(J,I),J=1,MGRID(I))
         WRITE(45,1000) 'JF = ',(JF(J,I),J=1,MGRID(I))
         WRITE(45,1000) 'IH = ',(IH(J,I),J=1,MGRID(I))
 100  CONTINUE
 1000 FORMAT(9X,A5,5I10)
       
      IF (ITURB >= 3) THEN
         MAXTB = MAXB
         MAXT2 = MAX2
      ELSE
         MAXTB = 1
         MAXT2 = 1
      ENDIF

      IF (ITURB >= 20) THEN
         MAXRB = MAXB
      ELSE
         MAXRB = 1
      ENDIF

      IF(ISTRES > 0) THEN
         MAXEB = MAXB
      ELSE
         MAXEB = 1
      ENDIF

      IF(IPRESC == 1 .OR. FRESUL) THEN
         MAXPC = MAXB
      ELSE
         MAXPC = 1
      ENDIF
          
      IF (MULPHL) THEN
         MAXMP = MAXB
         MAXP2 = MAX2
      ELSE
         MAXMP = 1
         MAXP2 = 1
      ENDIF

      IF (NSCAL > 0) THEN
         MAXSB = MAXB
         MAXS2 = MAX2
      ELSE
         MAXSB = 1
         MAXS2 = 1
      ENDIF

      NMAXSB   = MAX0(MAX0(0,NSCAL)*MAXSB,1)
      NMAXS2   = MAX0(MAXS2*NSCAL  ,1)

      IF (TIMEL) THEN
         MAXTI = MAXB
         MAXTS = MAXSB
         NMAXTS= NMAXSB
      ELSE
         MAXTI = 1
         MAXTS = 1
         NMAXTS= 1
      ENDIF

*      IF (COORL) THEN
*         MAXCO = MAXB
*      ELSE
*         MAXCO = 1
*      ENDIF

      MAXCO = MAXB
      
      IF (STRESL) THEN
         MAXSS = MAXB
      ELSE
         MAXSS = 1
      ENDIF

      IF(NCHIM > 0) THEN  ! CHIMERA IS IN USE
         MAXCH = MAXB
         MAXSCH= MAXCH
      ELSE
         MAXCH = 1
         MAXSCH= 1
      ENDIF

      IF(FRESUL) THEN
         MAXFS = MAXB
      ELSE
         MAXFS = 1
      ENDIF

      IF(TRANSL) THEN
         MAXTR = MAXB
      ELSE
         MAXTR  = 1
      ENDIF

      NMAXSC = MAX0(MAX0(0,NSCAL)*MAXSCH,1)

C MEMORY FOR XHULL,YHULL,ZHULL
      IF(IFSBC == 1) CALL MEMHUL(MAXHB,LEVEL,MGRID)

C ... Estimate the size of the work space

      IP  = 1    + 29*MAXB                  ! GEOMETRY
      IM1 = IP   + 50*MAXB + 16*MAXB        ! BASIC FLOW + K-E
      IE1 = IM1  +  9*MAXEB                 ! EARSM
      IM2 = IE1  +  5*NMAXSB                ! REYNOLDS STRESSES OR SCALARS
      IM3 = IM2  +  5*6*MAXRB + MAXRB       ! SOURCE TERMS IN RSM
      IM4 = IM3  +  6*MAXSS + 10*MAXB       ! SHEAR STRESSES AND VORTICITY
      IM5 = IM4  +  23*MAXTI+8*MAXTI+4*NMAXTS ! Time int. block

      IM4B= IM4  + 31*MAXTI                 ! Start of scalars
c     IM4C= IM4  + 7*MAXTI+2*MAXTI+4*NMAXTS ! Start of third time level
c     IM4D= IM4C + 5*MAXTI                  ! Start of turbulence quantities

      IM4E= IM4  + 2*(7*MAXTI)              ! ROAVE1 starting address
      IM4F= IM4E + 7*MAXTI                  ! RKAVE1 starting address
      IM4G= IM4E + 5*MAXTI+2*MAXTI          ! ROAVE2 starting address
      IM4H= IM4G + 5*MAXTI                  ! RKAVE2 starting address
      IM4I= IM1  - 5*MAXB                   ! RUVAV starting address

      IM6 = IM5                             ! old XLE2,.. now dummy
      IM7 = IM6  + 8*MAXCH+4*MAXCH+NMAXSC +
     +      5*MAXCH                         ! ROFOR blocksize
      IM60= IM6                             ! ROFOR
      IM6A= IM60 + 5*MAXCH                  ! RKFOR 
      IM6B= IM6A + 7*MAXCH                  ! FIFOR
      IM6C= IM6B + NMAXSC 
      IM10= IM7  + 18*MAX11                 ! GREPS 
      IM11= IM10 + (24+3)*MAXFS             ! Hull and Free surface
      IM12= IM11 + 17*MAXMP*NPHASES         ! Multiphase varibales
      IM13= IM12 + 16*MAXMP*NPHASES + 2*MAXMP ! Multiphase properties
      IM14= IM13 + 3*MAXPC                  ! Pressure correction variables
 
C ... double precision variables

      IMD1  = 1 + 6*MAXB +  12*MAXCO 
      IMD1A = 1 + 6*MAXB

      IIB   = IM14 +  (5+1)*MAX2 + 2*MAX2   ! add 2. level variables
      IIB2  = IIB  +  NMAXS2                ! add 2. level scalars
      IIB22 = IIB2 +  2*NPHASES*MAXP2       ! add 2. level multiphase variables
     
      IIB3  = IIB22+  22*IB                 ! add size IB
      IIB4  = IIB3 + MAXB                   ! add size JET (mersu)
      IIB5  = IIB4 + 
     +       (47+25+MAX(NSCAL,1))*IBF       ! add size IBF
      IIB6  = IIB5 + IBF                    ! add size WH(IBF)

      MREQ  = IIB6 - 1 + 
     +       6*NPTIM*MAXB+12*NPTIM*MAXCO    ! add size XCO and XLE2

      IF (NPNTS < IXMA*IYMA*IZMA) THEN
         WRITE(45,*) ' IXMA*IYMA*IZMA is larger than NPNTS',
     +   IXMA*IYMA*IZMA,NPNTS
         WRITE(45,*) ' NPNTS = IXMA*IYMA*IZMA'
         NPNTS = IXMA*IYMA*IZMA
      ENDIF

      IPGRI = 2*(IP/2) + 1

      IF (MREQ <= (NPNTS*3*NPTIM+IPGRI-1)) THEN
         WRITE(45,*)'NPNTS memory array larger (NPNTS*3*NPTIM+IPGRI-1)',
     +        ',(MREQ)',
     +        (NPNTS*3*NPTIM+IPGRI-1),MREQ
      ENDIF

      MREQ  = MAX(MREQ,(NPNTS*3*NPTIM+IPGRI-1))
      MREQB = (MREQ+IMD1*NPTIM)*MBITS

      CALL PRIME(MREQ+IMD1*NPTIM,IPRO,NPRO,MBITS)

      IERRCODE = 3 ! for ftnchek
      IERR3    = 3 ! for ftnchek

C ... Start to allocate the main arrays

      CALL ALLOCA(MAXCO,MAXRB,NMAXTS,MAX2,MAXP2,IERRCODE)

C ... Put starting values to zero. (A historical note)?
        
      CALL ZEROMA

      IF (IERRCODE > 0) THEN
         WRITE( *,1022)  IPRO,REAL(MREQB)/1048576.
         WRITE(45,1022)  IPRO,REAL(MREQB)/1048576.
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

      IF(PARALLEL .AND. MASTER) THEN
         WRITE( *,1030)
         WRITE(45,1030)
      ELSE IF(.NOT.PARALLEL) THEN
         WRITE( *,1031)
         WRITE(45,1031)
      ENDIF
      
      WRITE(45,1040) IB,MAXB,MAX2,MAX11,MAXTB,MAXT2,MAXRB,MAXEB,
     &     MAXSB,MAXSS,MAXS2,NMAXSB,MAXPC,MAXMP,MAXP2,MAXTI,MAXTS,
     &     MAXTR,NMAXTS,MAXCH,MAXFS,
     &     NPNTS,IP,IM1,IM2,IM3,IM4,IM5,IM6,IM7,IM10,IM11,IM12,
     &     IM13,IM14,
     &     IIB,IIB2,IIB3,IIB4,IIB5,IIB6,1,IP-1,24,    ! VOL,A1
     &     IP,47*MAXB,47,IP+47*MAXB,12*MAXB,12,       ! RK
     &     IM1,9*MAXEB,9,IE1,5*NMAXSB,5*MAX(NSCAL,1),IM2,31*MAXRB,31,
     &     IM3,6*MAXSS,6,IM3+6*MAXSS,10*MAXB,10,      ! S11,VORT
     &     IM4,33*MAXTI,33,
     &     IM4B,4*NMAXTS,4*MAX(NSCAL,1),
     &     IM6,5*MAXCH,5,IM6A,2*MAXCH,2,              ! ROFOR, RKFOR
     &     IM6B,NMAXSC,MAX(NSCAL,1),                  ! FIFOR
     &     IM6C,5*MAXCH,5,                            ! II1
     &     IM7,18*MAX11,18,IM10,24*MAXFS,27,          ! XXI
     &     IM11,17*MAXMP*NPHASES,17*NPHASES,          ! VAR
     &     IM12,15*MAXMP*NPHASES+2*MAXMP,15*NPHASES+2,! PRO
     &     IM13,3*MAXPC,3,IM14,5*MAX2,5,              ! PRC,2.level
     &     IM14+5*MAX2,2*MAX2,2,IIB,NMAXS2,MAX(NSCAL,1),
     &     IIB2,2*NPHASES*MAXP2,2*NPHASES,            ! P2H
     &     IIB22,18*IB,18,                            ! UBI (mersu)
     &     IIB22+18*IB,4*IB,4,                        ! ICP
     &     IIB3,MAXB,1,                               ! JET
     &     IIB4,47*IBF,47,                            ! HFLUX
     &     IIB4+47*IBF,(25+MAX(NSCAL,1))*IBF,25+MAX(NSCAL,1),! BOUNMF
     &     IIB5,IBF,1,                                ! WH
     &     IIB6,6*NPTIM*MAXB,6,                       ! XCO
     &     IIB6+6*NPTIM*MAXB,12*NPTIM*MAXCO,12        ! XLE2
 1040 FORMAT(9X,'IB     =',I10,
     &     /9X,'MAXB   =',I10,'  MAX2   =',I10,'  MAX11  =',I10,
     &     /9X,'MAXTB  =',I10,'  MAXT2  =',I10,'  MAXRB  =',I10,
     &     /9X,'MAXEB  =',I10,'  MAXSB  =',I10,'  MAXSS  =',I10,
     &     /9X,'MAXS2  =',I10,'  NMAXSB =',I10,'  MAXPC  =',I10,
     &     /9X,'MAXMP  =',I10,'  MAXP2  =',I10,'  MAXTI  =',I10,
     &     /9X,'MAXTS  =',I10,'  MAXTR  =',I10,'  NMAXTS =',I10,
     &     /9X,'MAXCH  =',I10,'  MAXFS  =',I10,
     &     /9X,'NPNTS  =',I10,/
     &     /9X,'IP     =',I10,
     &     /9X,'IM1    =',I10,'  IM2    =',I10,'  IM3    =',I10,
     &     /9X,'IM4    =',I10,'  IM5    =',I10,'  IM6    =',I10,
     &     /9X,'IM7    =',I10,'  IM10   =',I10,'  IM11   =',I10,
     &     /9X,'IM12   =',I10,'  IM13   =',I10,'  IM14   =',I10,
     &     /9X,'IIB    =',I10,'  IIB2   =',I10,'  IIB3   =',I10,
     &     /9X,'IIB4   =',I10,'  IIB5   =',I10,'  IIB6   =',I10,//
     &     /9X,'CURRENT MEMORY ADDRESS MAPPING',
     &     /9X,'==============================',/
     &     /9X,'Variable  Base Address A()    Size (words)    Variables',
     &     /9X,60('='),
     &     /9X,'VOL        ',2I15,I14,
     &     /9X,'RO         ',2I15,I14,
     &     /9X,'RK         ',2I15,I14,
     &     /9X,'BIJ        ',2I15,I14,
     &     /9X,'FI         ',2I15,I14,
     &     /9X,'PROD       ',2I15,I14,
     &     /9X,'S11        ',2I15,I14
     &     /9X,'VORT       ',2I15,I14
     &     /9X,'ROLE2      ',2I15,I14,
     &     /9X,'FILE2      ',2I15,I14,
     &     /9X,'ROFOR      ',2I15,I14,
     &     /9X,'RKFOR      ',2I15,I14,
     &     /9X,'FIFOR      ',2I15,I14,
     &     /9X,'II1        ',2I15,I14,
     &     /9X,'W12        ',2I15,I14,
     &     /9X,'XXI        ',2I15,I14,
     &     /9X,'VAR        ',2I15,I14,
     &     /9X,'PRO        ',2I15,I14,
     &     /9X,'PRC        ',2I15,I14,
     &     /9X,'ROP2H      ',2I15,I14,
     &     /9X,'RKP2H      ',2I15,I14,
     &     /9X,'FIP2H      ',2I15,I14,
     &     /9X,'P2H        ',2I15,I14,
     &     /9X,'UBI        ',2I15,I14,
     &     /9X,'ICP        ',2I15,I14,
     &     /9X,'JET        ',2I15,I14,
     &     /9X,'HFLUX      ',2I15,I14,
     &     /9X,'BOUNMF     ',2I15,I14,
     &     /9X,'WH         ',2I15,I14,
     &     /9X,60('='),
     &     /9X,'Variable Base Address AD()  Size (words)     Variables',
     &     /9X,60('='),
     &     /9X,'XCO        ',2I15,I14,
     &     /9X,'XLE2       ',2I15,I14//)

C     -----  Replacement of the Subroutine ZERO, Part 2 -----
      DO I = 1,MAXB
        VOL(I) = 1.0 ! This prevents unrealistic geometry in ghost-cells
      ENDDO
C     -----  Replacement #2 Ends Here -----------------------
         
*******************************************************************************
C ... Execution in NS3D begins
*******************************************************************************
  
      CALL NS3D
           
*******************************************************************************
*******************************************************************************

C ... Simulation is completed, start to finalize

      IF(MASTER) THEN
      OPEN(UNIT=1,FILE='RUN',STATUS='UNKNOWN')
      WRITE(1,'(A5)') 'norun'
      WRITE(1,'(I7,3X,A5)') 0,'ICMAX'
      CLOSE(1)
      ENDIF

      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

      WRITE(13,*)
      IF(.NOT. TIMEL) THEN
      WRITE(13,'(A,I7)') '  The current cycle is ',ICYCLE
      ELSE IF(TIMEL) THEN
      WRITE(13,'(A,F6.3)') '  The current time is T =',T
      ENDIF

      IF(VALUES(3) == IPAIVA) THEN 
      WRITE(13,'(A,I2,A1,I2,A1,I2,A)') '  Simulation is finished at ',
     + VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ELSE IF (VALUES(2) >= 10) THEN
      WRITE(13,'(A,I2,A1,I2,A1,I4,A,I2,A1,I2,A1,I2,A)') 
     + '  Simulation is finished ',
     + VALUES(3),'.',VALUES(2),'.',VALUES(1),' at ',
     + VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ELSE IF (VALUES(2) < 10) THEN
      WRITE(13,'(A,I2,A1,I1,A1,I4,A,I2,A1,I2,A1,I2,A)') 
     + '  Simulation is finished ',
     + VALUES(3),'.',VALUES(2),'.',VALUES(1),' at ',
     + VALUES(5),'.',VALUES(6),'.',VALUES(7),' oclock'
      ENDIF

      CALL WALL_CLOCK_TIME(VALUES,IMONTH,IPAIVA,IHOUR,IMINUT,ISECON)
         
      IF(VALUES(3)-IPAIVA  <= 2) THEN
         WRITE(13,'(A)') '  Hermit : I will be resting now... '
      ELSE
         WRITE(13,'(A)') '  Hermit : A long run, I need rest now... '
      ENDIF

      IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      WRITE(13,*) 'END.'
      CLOSE(13)

      IF(PARALLEL) CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      IF(PARALLEL .AND. MASTER) THEN ! collect the ERROR.LOG
         OPEN(13,FILE='RUN.LOG',STATUS='UNKNOWN',FORM='FORMATTED')
         DO NP = 1,NPRO
            LRIVI = .TRUE.
            CALL NUMCH3(CN,NP)
            OPEN(14,FILE='RUN.LOG'//CN,STATUS='OLD',FORM='FORMATTED')
c 113        READ(14,'(A80)') RIVI
 113        READ(14,'(A80)',END=123) RIVI
            IF(RIVI(1:2) ==  '**') THEN
               IF(LRIVI) WRITE(13,*) 'Errors in proces: ',NP
               LRIVI = .FALSE.
               WRITE(13,'(A80)') RIVI
               GOTO 113
            ENDIF
 123        continue
            CLOSE(14)
         ENDDO
         WRITE(13,*) 'END.'
         CLOSE(13)
      ENDIF

      IF(MASTER) THEN

         OPEN(13,FILE='RUN.LOG',STATUS='OLD',FORM='FORMATTED')
         READ(13,*,ERR=111,END=111) RUN

         IF(RUN(1:2) == '**') GOTO 111

         GOTO 112

 111     CONTINUE

         WRITE(*,*) run(1:2)
         WRITE(*,*) 'Some errors have been found. Look more details'
         WRITE(*,*) 'from file RUN.LOG'
         
 112     CONTINUE

         CLOSE(13)

      ENDIF

      CALL DEALLO
      CALL DEALNS

      IF(PARALLEL) CALL MPI_FINALIZE(RC)

      WRITE(45,1050)
      IF(MASTER) WRITE(*,1050)
c 1022 FORMAT(/'  MAIN :  Not enough memory in process ',I3,
c     +     '. Desired ',F6.2,'MB. aborting...'/)
c 1030 FORMAT(/'  MAIN :  In master the memory is successfully allocated'
c     + ' for the work arrays!')
c 1031 FORMAT(/'  MAIN :  The memory is successfully allocated '
c     + 'for the work arrays!')
c 1050 FORMAT(/'  MAIN :  The allocated memory has been released!'/)
 1022 FORMAT(/'  MAIN :    Not enough memory in process ',I3,
     +     '. Desired ',F6.2,'MB. aborting...'/)
 1030 FORMAT(/'  MAIN :    Memory successfully allocated'
     + ' for the work arrays in master!')
 1031 FORMAT(/'  MAIN :    Memory successfully allocated '
     + 'for the work arrays!')
 1050 FORMAT(/'  MAIN :    The allocated memory has been released!'/)
c        write(6,*) 'voihan vuohi' ! (Tarvitset tätä vielä joskus)!

      CLOSE(4)
      CLOSE(45)
      CLOSE(13)

      STOP
      END PROGRAM FINFLO

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE COMSI(MGSTP,MBSTP,IB,MAXB,MAX11,MAX2,MAXWNE,MAWNE,
     +                 NPNTS,IBCOM,IBCOM2,IBSLB,NBMAX,IBF,OMEGAREF,
     +                 NAMEL,NPPN,IOLD1,GRIDFI,BCFILE,OUTPUT)

      USE INTEGERS, ONLY : NBGG, NB, MGM, NPPV, MBPRO, NPRO

      IMPLICIT NONE

      INTEGER :: IBCOM(*), IBCOM2(*), IBSLB(*)
      INTEGER :: MGSTP, MBSTP, IB, MAXB, MAX11, MAX2, MAXWNE, MAWNE,
     &           NPNTS, NBMAX, IBF, MB, MS, NPPN, IOLD1
      REAL    :: OMEGAREF
      LOGICAL :: NAMEL
      CHARACTER(LEN=80) :: NAME, GRIDFI, BCFILE, OUTPUT

      CALL INPUT(NPPN,NAME,NAMEL,IOLD1,GRIDFI,BCFILE,OUTPUT,OMEGAREF)
 
      CALL INPUT2(NPPN,OMEGAREF,NAME)
 
      CALL COMPUT(IBCOM,IBCOM2,IBSLB,MB,MS,MGSTP,MBSTP,NPNTS,IBF,
     &            MGM)
 
      CALL FINAL(IBCOM,IBCOM2,IBSLB,MB,MS,MAXB,MAX11,MAX2,IB,MAXWNE,
     &     MAWNE,NBMAX)

      RETURN
      END SUBROUTINE COMSI

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE INPUT(NPPN,NAME,NAMEL,IOLD1,GRIDFI,BCFILE,
     &                 OUTPUT,OMEGAREF)

      USE NS3CO
      USE CONSTANTS, ONLY : DEG2RAD, PII, EPS6
      USE BLADE_VARIABLES, ONLY : SHAFT, RCG, RITH, RIZE, RIBE,DTB,
     &                            THETACOL,THCYCLON,THCYCLAT,ETIP,CBLADE
      USE FORCE_GROUPS
      USE INTEGERS, ONLY : IXMA, IYMA, IZMA, IPRO, NBCS, IREPEA, NREPEA,
     &                     MCYCAM, ITERAD, ITERAC, NUMCOL, NMOV, 
     &                     NBGG, NB, MGM, NPPV, MBPRO, NPRO
      USE MPI
      USE MAIN_ARRAYS, ONLY : JSTATE, IMAXG, JMAXG, KMAXG,  
     &                        IDOB,  JDOB,  KDOB, NSSB, ISSB, 
     &                        NCHIMT, NCHOR, NLOCAL, NPNUM, OMEGA,
     &                        OMEGAX, OMEGAY, OMEGAZ, 
     &                        CENAX, CENAY, CENAZ, AMPL, OSCLS,
     &                        BXMIN, BXMAX, BYMIN, BYMAX, BZMIN, BZMAX, 
     &                        IGRID, MGRIDA, BLKS, ICNH,
     &                        IBTGR, ITPGR, JBTGR, JTPGR, KBTGR, KTPGR,
     &                        LHULL, MHULL, MMHUL, ICGH, ICOGH, ICONH,
     &                        IBOTGR, ITOPGR, 
     &                        JBOTGR, JTOPGR, 
     &                        KBOTGR, KTOPGR,
     &                        NPROCE, ICON, RCON, 
     &                        IMAX, JMAX, KMAX,
     &                        INTERI, INTERJ, INTERK, IDER, LAMIN, 
     &                        IT, IL, IK, INITC, IDI1, IDI2, IDI3,
     &                        MOV, MGRID, IROTVE, 
     &                        MIB, MJB, MKB, MIT, MJT, MKT,
     &                        IMAXM, JMAXM, KMAXM, IMINM, JMINM, KMINM,
     &                        MOVPO, RLOLIM, UPPLIM, 
     &                        PSIGSC, PSIGS2, IUPPT, NEARBLOCK

      USE FLIGHT, ONLY : NGRIFL, OSKU, TRMODE,
     &                   XCGI, YCGI, ZCGI, XCG, YCG, ZCG,
     &                   PSIR, THETAR, PHIR, PSIRI, THETARI, PHIRI, 
     &                   PSIM, XCGIS, YCGIS, ZCGIS, DAMPC1, DAMPC2,
     &                   DRAUGHTI ,DRAUGHT, SINKI, SINK, TRIMAI, TRIMA,
     &                   RIN, ROUT, THRUST, TORQUE, ADV, 
     &                   ROTA1, ROTB1, ROTA1I, ROTB1I, CONEA, CONEAI,
     &                   ACTUA, VTIP, CDBLADE, CFBLADE,SIGMA, CFTI,
     &                   IFA, IFT, IFR, NSERIES, IGRSLAVE, 
     &                   ROTA1S, ROTB1S, MVSHIP, SHIP, SHIPWAVE,
     &                   ACTDISK, HROTOR, HMVBLADE, FLYOBJ, SHIPPROP,
     &                   RTMSP,FDSP,UTSP,FXSP,FXTSP,QFACT,NBLADE,RGML,
     &                   XFAKEP,YFAKEP,ZFAKEP,ROTB1FAKEP,IFAKEP
   
      USE CHARACTERS, ONLY : BOUNDF, BOUNDN

      USE MPCCIVARS, ONLY : JOINTS_FILE, MODE_SHAPES_FILE

      USE GM


      IMPLICIT NONE

      INTEGER, PARAMETER :: INPAR=67, IRNPAR=60,  LRNPAR=30, INSCAL=36,
     +                      IRSCAL=76, ICNPAR=6, IFNPAR=108, INFNPAR=6

      INTEGER :: NPPN, IPHASE
      INTEGER :: IX(NBGGG), IY(NBGGG), IZ(NBGGG), ISCALA(INSCAL)
      INTEGER :: STATUS(MPI_STATUS_SIZE)

      CHARACTER(LEN=80) :: NAME, GRIDFI, BCFILE, OUTPUT
      CHARACTER(LEN=80) :: LINE1,RIVI
      CHARACTER(LEN=3)  :: LINE3

      REAL ::  RSCALA(IRSCAL)
      REAL, ALLOCATABLE :: RBC(:,:),RBC2(:,:),RNPARA(:,:),FNPARA(:,:)
      INTEGER, ALLOCATABLE :: IBC(:,:),IBC2(:,:)

      CHARACTER(LEN=80), ALLOCATABLE :: BCAPU(:),BCAPU2(:)
      CHARACTER(LEN=80), ALLOCATABLE :: BNAPU(:),BNAPU2(:)
      CHARACTER(LEN=80), ALLOCATABLE :: CNPARA(:,:)
      CHARACTER(LEN=3),  ALLOCATABLE :: TBC(:)
      CHARACTER(LEN=1),  ALLOCATABLE :: CFPARA(:)

      INTEGER,      ALLOCATABLE :: INPARA(:,:),INFPAR(:,:)

      LOGICAL :: THERE, LVERSIO, IFOUND, MASTER
      LOGICAL :: LNPARA(LRNPAR), NAMEL, TRANSIT

      INTEGER :: ITURBO, N, NBGRID, ICHA, NDIV, MGRIDN,
     &           IMAXA, JMAXA, KMAXA, NTMAX, NBLKS, I, J, K, 
     &           NPPVS, IFACE, 
     &           NOFBCS, LAMAPU, I0, NBL, IBT, NGL,
     &           IERR1, IERR2, IERR3, IERR4,  
     &           ICYCTI, ICYTIT, MAXNSC, NS, IWR77, NMOVA, III,
     &           ERRORCODE, IERR, IB, IP, NP, IOTY, IOLD1, ITYPE

      REAL :: EPS, FRPR, FRSCP, FRSCH, FRSMUD, RKII, REPSI, RMULIM,
     &        RUUPRE, RKIIB, FRSMUTB, REPSIB, RMUINIB, RMULIMB, RKHI,
     &        PSEUCOB, OMEGAREF

C ... These are for a proper use of water and steam state routines

      REAL, DIMENSION(NPHASES) :: FRAVIS, FRAPRE, FRADEN, FRASIE, 
     &                            FRATEM, FRACP, FRACH, FRDRDP, FRDRDH,
     &                            FRASSP, FRARE, FRAMA

      NAMELIST /WORKS/ IOLD1, GRIDFI, BCFILE, OUTPUT, OUTVAR

      MASTER = IPRO == 1

C ... INPUT FOR THE THREE-DIMENSIONAL FLOW SOLVER

      EPS    = 1.E-6
      IXMA   = 0
      IYMA   = 0
      IZMA   = 0
      ITURBO = ITURB
      WOODYL = .FALSE. ! Hoida tämä inputtiin. En taida hoitaa..
      IFOUND = .FALSE.
         
C ********* READ INPUT. IN MPI ONLY THE MASTER PROCESS READS THIS *****

      IF(MASTER) THEN

C ... Set output format defaults

      OUTPUT='Plot3D'  ! Default output format is Plot3D (aka Fieldview)
      PLOT3D = .TRUE.
      ENSIGHT = .FALSE.

      OUTVAR = 0
      OUTVAR(1:73) = 1  ! Output of all currently available FUN variables
                        ! at the end of the simulation and during the time
                        ! accurate loop (if MOVIE files are requested) 
                        ! is the default 

C *** READ GRID FILE to get number of blocks and allocate memory for INPUT

      IF(NAMEL) READ(2,NML=WORKS)

      REWIND 2

      IF(OUTPUT == 'Plot3D' .OR. OUTPUT == 'PLOT3D' .OR. 
     &   OUTPUT == 'plot3d') THEN
         PLOT3D  = .TRUE.
         ENSIGHT = .FALSE.
      ENDIF

      IF(OUTPUT == 'EnSight' .OR. OUTPUT == 'ENSIGHT' .OR. 
     &   OUTPUT == 'ensight') THEN
         PLOT3D  = .FALSE.
         ENSIGHT = .TRUE.
      ENDIF

      IF(OUTPUT == 'Both' .OR. OUTPUT == 'BOTH' .OR. 
     &   OUTPUT == 'both') THEN
         PLOT3D  = .TRUE.
         ENSIGHT = .TRUE.
      ENDIF

       
      OPEN(41,FILE=GRIDFI,STATUS='OLD',FORM='UNFORMATTED')
          
      READ (41) NBLOCK
      WRITE(13,'(A,A70)')'  Opening ',GRIDFI
      WRITE(13,'(A,I4)') '  Number of grid blocks = ',NBLOCK
      WRITE(45,'(A,A70)')'  Opening ',GRIDFI
      WRITE(45,'(A,I4)') '  Number of grid blocks = ',NBLOCK
      WRITE(45,*)
      READ (41) (IX(N),IY(N),IZ(N),N=1,NBLOCK)
          
      IF(NBLOCK > NBGGG) THEN
         WRITE(*,*) 'Must increase parameter NBGGG in character.f'
         WRITE(*,*) 'at least value of ',NBLOCK
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      NBGG    = NBLOCK + 1  ! MAXIMUM NUMBER OF BLOCKS IN ALL PROCESSES
      NBGRID  = NBLOCK      ! something to do with BC file. Can remember? 
      NBLOCG  = NBLOCK      ! global number of blocks


C ... Number of processes must not exceed the total number of grid blocks 
C ... in parallel runs. 

      IF(PARALLEL .AND. NPRO > NBLOCG) THEN
         WRITE(*,*)
         WRITE(*,*) 
     &        'JOB ABORTED: An automatic parallelization is not ',
     &        'possible due to the required number of procceses ', 
     &        'exceeding the number of grid blocks!'
         WRITE(*,*) 
         CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

C ... Allocation for process no. 1

      ALLOCATE(INPARA(INPAR ,NBLOCG))
      ALLOCATE(RNPARA(IRNPAR,NBLOCG))
      ALLOCATE(FNPARA(IFNPAR,19))
      ALLOCATE(INFPAR(INFNPAR,19))
      ALLOCATE(CNPARA(ICNPAR ,NBLOCG))
      ALLOCATE(CFPARA(19))
      ALLOCATE(JSTATE(NBLOCG,3))
      ALLOCATE(IMAXG(NBLOCG))
      ALLOCATE(JMAXG(NBLOCG))
      ALLOCATE(KMAXG(NBLOCG))
      ALLOCATE(IDOB(NBLOCG))
      ALLOCATE(JDOB(NBLOCG))
      ALLOCATE(KDOB(NBLOCG))
      ALLOCATE(NSSB(NBLOCG))
      ALLOCATE(ISSB(7,NBLOCG))
      ALLOCATE(NCHIMT(NBLOCG))
      ALLOCATE(NCHOR(NBLOCG))
      ALLOCATE(NLOCAL(NBLOCG))
      ALLOCATE(NPNUM(NBLOCG))
      ALLOCATE(OMEGA(NBLOCG))
      ALLOCATE(OMEGAX(NBLOCG))
      ALLOCATE(OMEGAY(NBLOCG))
      ALLOCATE(OMEGAZ(NBLOCG))
      ALLOCATE(CENAX(NBLOCG))
      ALLOCATE(CENAY(NBLOCG))
      ALLOCATE(CENAZ(NBLOCG))
      ALLOCATE(AMPL(NBLOCG))
      ALLOCATE(OSCLS(NBLOCG))
      ALLOCATE(BXMIN(NBLOCG))
      ALLOCATE(BXMAX(NBLOCG))
      ALLOCATE(BYMIN(NBLOCG))
      ALLOCATE(BYMAX(NBLOCG))
      ALLOCATE(BZMIN(NBLOCG))
      ALLOCATE(BZMAX(NBLOCG))
      ALLOCATE(IGRID(NBLOCG,5))
      ALLOCATE(MGRIDA(NBLOCG))
      ALLOCATE(BLKS(NBGG))
      ALLOCATE(NEARBLOCK(NBLOCG*NBLOCG))

      ALLOCATE(ICNH(NBLOCG))
      ALLOCATE(IBTGR(NBLOCG))
      ALLOCATE(ITPGR(NBLOCG))
      ALLOCATE(JBTGR(NBLOCG))
      ALLOCATE(JTPGR(NBLOCG))
      ALLOCATE(KBTGR(NBLOCG))
      ALLOCATE(KTPGR(NBLOCG))
      ALLOCATE(LHULL(NBGG))
      ALLOCATE(MMHUL(NBGG))
      ALLOCATE(MHULL(NBGG))
      ALLOCATE(ICGH(NBGG))
      ALLOCATE(ICOGH(NBGG))
      ALLOCATE(ICONH(NBGG))
      ALLOCATE(IBOTGR(NBGG))
      ALLOCATE(ITOPGR(NBGG))
      ALLOCATE(JBOTGR(NBGG))
      ALLOCATE(JTOPGR(NBGG))
      ALLOCATE(KBOTGR(NBGG))
      ALLOCATE(KTOPGR(NBGG))
      ALLOCATE(IUPPT(NBGG))

C ... Set a default value for the solution type
   
      INPARA(:,:) = 0
      RNPARA(:,:) = 0.
      JSTATE(:,:) = 0
      MGRIDA(:)   = 0

      BLKS(1:NBLOCG)%SOLUTION_TYPE = 'FLUID'
      BLKS(1:NBLOCG)%MATERIAL(1)   = 'DUMMY'
      BLKS(1:NBLOCG)%MATERIAL(2)   = 'DUMMY'
C     BLKS(1:NBLOCG)%MATERIAL(3)   = 'DUMMY'
      BLKS(1:NBLOCG)%NPHASE = 1
      BLKS(1:NBLOCG)%IDIFF  = 1
      NPHASE = 1
      TWO_FLUIDL = .FALSE.

C ... Memory is allocated for global arrays

C ... READ INPUT in subprogram RINP or from a namelist

      ICHA = 4
      LVERSIO = .FALSE.


      IF(NAMEL) THEN

         CALL SET_RINP(INPARA,INPAR,RNPARA,IRNPAR,LNPARA,LRNPAR,
     +    ISCALA,INSCAL,RSCALA,IRSCAL,CNPARA,ICNPAR,ISTATE,IFLUX,
     +    IX,IY,IZ,NAME,NPHASES,FNPARA,CFPARA,IFNPAR,NGRIFL,INFPAR,
     +    INFNPAR)

         CALL BACK_RINP(ISCALA,INSCAL,RSCALA,IRSCAL,LNPARA,LRNPAR)
           
         CALL BACK_FLIGHT(FNPARA,CFPARA,IFNPAR,INPARA,
     +                    INPAR,INFPAR,INFNPAR)

      ENDIF


      WRITE(45,*) ' Grid dimensions at the finest grid level:'
      WRITE(45,*) ' Block   IMAX      JMAX       KMAX   =  SIGMA'
      WRITE(45,*) ' ============================================'

      WRITE(45,*) ' Assumed calculation level:',LEVEL
      IF(LEVEL > 5) WRITE(45,*) 'Assumed grid level is high?!'
      IF(LEVEL > 9) THEN
         WRITE(45,*) ' Assumed grid level is too high'
         WRITE(13,*) ' Assumed grid level is too high. I stopped'
         WRITE(*,*)  ' Assumed grid level is too high. I stopped'
         STOP
      ENDIF
      IF(LEVEL == 0) THEN
         WRITE(45,*) ' Assumed grid level is zero. I stopped.'
         WRITE(13,*) ' Assumed grid level is zero, check input data.'
         WRITE(*,*)  ' Assumed grid level is zero, check input data.'
         STOP
      ENDIF

      MTOT = 0
      MGM  = 0
      OMEGAREF = 0.
      
      DO 1000 N = 1,NBLOCG ! Global N

      IF(NAMEL) THEN ! Global arrays, local ones in MPITOU
        
         JSTATE(N,1)    = INPARA(25,N)
         JSTATE(N,2)    = INPARA(26,N)
         JSTATE(N,3)    = INPARA(27,N)
 ! INPARA(28,N) is free now, IGRID is changed to positions between 44...48
         NCHIMT(N)      = INPARA(29,N)
         BLKS(N)%NPHASE = INPARA(30,N)
         BLKS(N)%IPRESC = INPARA(31,N)
         BLKS(N)%LUSGS  = INPARA(32,N)
         BLKS(N)%IFLUX  = INPARA(33,N)
         ICNH(N)        = INPARA(34,N)
         IBTGR(N)       = INPARA(35,N)
         ITPGR(N)       = INPARA(36,N)
         JBTGR(N)       = INPARA(37,N)
         JTPGR(N)       = INPARA(38,N)
         KBTGR(N)       = INPARA(39,N)
         KTPGR(N)       = INPARA(40,N)
         BLKS(N)%EVAP_TYPE = INPARA(41,N)
         BLKS(N)%IDIFF  = INPARA(42,N)
         BLKS(N)%MPCASE = INPARA(43,N)
         IGRID(N,1)     = INPARA(44,N)
         IGRID(N,2)     = INPARA(45,N)
         IGRID(N,3)     = INPARA(46,N)
         IGRID(N,4)     = INPARA(47,N)
         IGRID(N,5)     = INPARA(48,N)
         BLKS(N)%PREVEL = INPARA(49,N)
         BLKS(N)%INLRC  = INPARA(50,N)
         BLKS(N)%OUTRC  = INPARA(51,N)
         BLKS(N)%ICERO  = INPARA(52,N)
         BLKS(N)%ITERMP = INPARA(53,N)
         BLKS(N)%IUPTEM(1) = INPARA(54,N)
         BLKS(N)%IUPTEM(2) = INPARA(55,N)
C         BLKS(N)%IUPTEM(3) = INPARA(56,N)
         IUPPT(N)       = INPARA(57,N)
         BLKS(N)%INTERI = INPARA(58,N)
         BLKS(N)%INTERJ = INPARA(59,N)
         BLKS(N)%INTERK = INPARA(60,N)
         BLKS(N)%TURMUL = INPARA(61,N)
         BLKS(N)%DRAGMUL= INPARA(62,N)
         BLKS(N)%LIFTMUL= INPARA(63,N)
         BLKS(N)%DISPERMUL = INPARA(64,N)
         BLKS(N)%VMASS  = INPARA(65,N)
         BLKS(N)%WFORCEMUL = INPARA(66,N)
         OSCLS(N)       = INPARA(67,N) 

         ROTAT(N)       = RNPARA(1,N)
         OMEGA(N)       = RNPARA(2,N)
         OMEGAX(N)      = RNPARA(3,N)
         OMEGAY(N)      = RNPARA(4,N)
         OMEGAZ(N)      = RNPARA(5,N)
         CENAX(N)       = RNPARA(6,N)
         CENAY(N)       = RNPARA(7,N)
         CENAZ(N)       = RNPARA(8,N)
         AMPL(N)        = RNPARA(60,N)
         BLKS(N)%SOLUTION_TYPE = CNPARA(1,N)
         CALL YESNO(CNPARA(2,N),BLKS(N)%CONVL)
         CALL YESNO(CNPARA(3,N),BLKS(N)%COMPCORR)
         CALL YESNO(CNPARA(4,N),BLKS(N)%KATOL)
         CALL YESNO(CNPARA(5,N),BLKS(N)%ZEROV)
         CALL YESNO(CNPARA(6,N),BLKS(N)%FLUXCORR)
         BLKS(N)%FRSDEN = RNPARA(9,N)
         BLKS(N)%FRSPRE = RNPARA(10,N) ! Number is misleading
         BLKS(N)%FRSTEM = RNPARA(26,N)
         BLKS(N)%DFRSTEM= RNPARA(26,N)
         BLKS(N)%FRSVEL = RNPARA(11,N)
         BLKS(N)%ARTSSP = RNPARA(12,N)
         BLKS(N)%SOLTEM = RNPARA(13,N)
         BLKS(N)%FRSTUR = RNPARA(14,N)
         BLKS(N)%FRSMUT = RNPARA(15,N)
         BLKS(N)%TURLIM = RNPARA(16,N)
         BLKS(N)%RE     = RNPARA(17,N)
         BLKS(N)%PR     = RNPARA(18,N)
         BLKS(N)%PRT    = RNPARA(19,N)
         BLKS(N)%CHLREF = RNPARA(20,N)
         BLKS(N)%RMACH  = RNPARA(21,N)
         BLKS(N)%RMUINI = RNPARA(22,N)
         BLKS(N)%AREF   = RNPARA(23,N)
         BLKS(N)%REFPRE = RNPARA(24,N)
         BLKS(N)%DIFPRE = RNPARA(25,N)
         BLKS(N)%TEMINI = RNPARA(27,N)
         BLKS(N)%RMULTV = RNPARA(28,N)
         BLKS(N)%REFVEL = RNPARA(37,N)
         BLKS(N)%ALFAP  = RNPARA(38,N)
         BLKS(N)%ALFAU  = RNPARA(39,N)
         BLKS(N)%CHLCAV = RNPARA(40,N)
         BLKS(N)%REFCAV = RNPARA(41,N)
         BLKS(N)%ALFASKEW = RNPARA(42,N)
         BLKS(N)%GVEX   = RNPARA(49,N)
         BLKS(N)%GVEY   = RNPARA(50,N)
         BLKS(N)%GVEZ   = RNPARA(51,N) 
         BLKS(N)%FREDIFI= RNPARA(52,N)
         BLKS(N)%FREDIFJ= RNPARA(53,N)
         BLKS(N)%FREDIFK= RNPARA(54,N) 
         BLKS(N)%RINTF  = RNPARA(55,N)
         BLKS(N)%SKEWCORR = .FALSE.
         IF(BLKS(N)%ALFASKEW > 0.) BLKS(N)%SKEWCORR = .TRUE.

         DO IPHASE = 1,NPHASES
         BLKS(N)%FRSALFA(IPHASE) = RNPARA(27+IPHASE,N) ! 28...30
         ENDDO
         BLKS(N)%TLOLIM = RNPARA(31,N)
         BLKS(N)%TUPLIM = RNPARA(32,N)
         BLKS(N)%CAVNO  = RNPARA(33,N)
         DO IPHASE = 1,NPHASES
         BLKS(N)%TAUF(IPHASE) = RNPARA(33+IPHASE,N)    ! 33...36
         BLKS(N)%N(IPHASE)    = RNPARA(42+IPHASE,N)    ! 43...45
         BLKS(N)%RK(IPHASE)   = RNPARA(55+IPHASE,N) ! 56...57
         BLKS(N)%FRADEN(IPHASE) = RNPARA(57+IPHASE,N) ! 58...59
         ENDDO
         BLKS(N)%FRSG   = RNPARA(46,N)
         BLKS(N)%FRSRET = RNPARA(47,N)
         BLKS(N)%VOIDTT = RNPARA(48,N)
      ELSE ! Defaults using global values in the old input file
         BLKS(N)%SOLUTION_TYPE = 'FLUID'
         BLKS(N)%CONVL  = CONVL
         BLKS(N)%NPHASE = NPHASE
         BLKS(N)%IPRESC = IPRESC
         BLKS(N)%LUSGS  = LUSGS
         BLKS(N)%IFLUX  = IFLUX
         BLKS(N)%FRSDEN = FRSDEN
         BLKS(N)%FRSPRE = FRSPRE
         BLKS(N)%FRSTEM = FRSTEM ! Number is misleading
         BLKS(N)%DFRSTEM= FRSTEM
         BLKS(N)%FRSVEL = FRSVEL
         BLKS(N)%ARTSSP = ARTSSP
         BLKS(N)%SOLTEM = SOLTEM
         BLKS(N)%FRSTUR = RKLIM
         BLKS(N)%FRSMUT = EPSLIM ! Mersu
         BLKS(N)%TURLIM = TURLIM
         BLKS(N)%RE     = RE
         BLKS(N)%PR     = PR
         BLKS(N)%PRT    = PRT
         BLKS(N)%CHLREF = CHLREF
         BLKS(N)%RMACH  = RMACH
         BLKS(N)%RMUINI = RMUINI
         BLKS(N)%AREF   = AREF
         BLKS(N)%REFPRE = REFPRE
         BLKS(N)%REFVEL = REFVEL
         BLKS(N)%DIFPRE = DIFPRE
         BLKS(N)%TEMINI = TEMINI
         BLKS(N)%RMULTV = RMULTV
         BLKS(N)%ALFAP  = ALFAP
         BLKS(N)%ALFAU  = ALFAU
         BLKS(N)%CHLCAV = CHLREF
         BLKS(N)%REFCAV = REFVEL
         BLKS(N)%GVEX   = GVEX
         BLKS(N)%GVEY   = GVEY
         BLKS(N)%GVEZ   = GVEZ         
c         BLKS(N)%FRSALFA(1) = 1.
         DO IPHASE = 1,NPHASES
            BLKS(N)%FRSALFA(IPHASE) = 1.
            BLKS(N)%FRADEN(IPHASE) = FRSDEN
         ENDDO
C ... Blockwise limits only from the namelist input. Defaults at line 2399
C        BLKS(N)%TLOLIM = TLOLIM
C        BLKS(N)%TUPLIM = TUPLIM
      ENDIF ! NAMEL
        
      NPHASE = MAX(NPHASE,BLKS(N)%NPHASE)     ! Find the number of phases

      IF(.NOT.MULPHL .AND. BLKS(N)%SOLUTION_TYPE == 'MULTI')THEN
         MULPHL = .TRUE.
         WRITE(4,*) 'SOLUTION TYPE equal to "MULTI" in BLOCK',N
         WRITE(4,*) 'MULPHL was changed to "TRUE"' 
         WRITE(*,*) 'SOLUTION TYPE equal to "MULTI" in BLOCK',N
         WRITE(*,*) 'MULPHL was changed to "TRUE"' 
      ENDIF
 
      IF(.NOT.TUR_MULPHL .AND. BLKS(N)%TURMUL /= 0) THEN
         TUR_MULPHL = .TRUE.
         WRITE(4,*) 'TURMUL not equal to zero in BLOCK',N
         WRITE(4,*) 'TUR_MULPHL was changed to "TRUE"' 
         WRITE(*,*) 'TURMUL not equal to zero in BLOCK',N
         WRITE(*,*) 'TUR_MULPHL was changed to "TRUE"' 
      ENDIF 
 
      IF(.NOT.TWO_FLUIDL .AND. BLKS(N)%SOLUTION_TYPE == 'MULTI')THEN
         TWO_FLUIDL = .TRUE.
      ENDIF 

      IF(.NOT.MULPHL .AND. BLKS(N)%SOLUTION_TYPE == 'CAVIT')THEN
         MULPHL = .TRUE.
         WRITE(4,*) ' SOLUTION TYPE equal to "CAVIT" in BLOCK',N
         WRITE(4,*) ' MULPHL was changed to "TRUE"' 
         WRITE(*,*) ' SOLUTION TYPE equal to "CAVIT" in BLOCK',N
         WRITE(*,*) ' MULPHL was changed to "TRUE"' 
      ENDIF

      IF(BLKS(N)%SOLUTION_TYPE == 'FLUID' .AND. BLKS(N)%NPHASE > 1 .OR.
     &   BLKS(N)%SOLUTION_TYPE == 'SOLID' .AND. BLKS(N)%NPHASE > 1)THEN
         WRITE(4,*) ' SOLUTION TYPE differs from NPHASE in BLOCK',N
         WRITE(4,*) ' Number of phases and type should match. Exiting..' 
         WRITE(*,*) ' SOLUTION TYPE differs from NPHASE in BLOCK',N
         WRITE(*,*) ' Number of phases and type should match. Exiting..'
         WRITE(13,*)' SOLUTION TYPE differs from NPHASE in BLOCK',N
         WRITE(13,*)' Number of phases and type should match. Exiting..'
         STOP
      ENDIF

      IF(.NOT.FULLNS .AND. BLKS(N)%IDIFF == 3)  THEN
         WRITE(4,*) 'Thin shear layer approximation selected in block',N
         WRITE(4,*) 'IDIFF was changed from 3 to 1' 
         WRITE(*,*) 'Thin shear layer approximation selected in block',N
         WRITE(*,*) 'IDIFF was changed from 3 to 1' 
      ENDIF 
      IF(.NOT.FULLNS .AND. BLKS(N)%IDIFF == 6)  THEN
         WRITE(4,*) 'Thin shear-layer approximation selected in block',N
         WRITE(4,*) 'IDIFF was changed from 6 to 5' 
         WRITE(*,*) 'Thin shear layer approximation selected in block',N
         WRITE(*,*) 'IDIFF was changed from 6 to 5' 
      ENDIF 
      IF(    BLKS(N)%IPRESC == 1 .AND. BLKS(N)%IFLUX <= 3 .AND.
     & .NOT. BLKS(N)%COMPCORR)  THEN
         WRITE(13,'(2A,I3,A/,2A)')'  Warning : A compressible flux-type'
     &  ,' is selected for block',N,', but it is','  not taken into'
     &  ,' account in pressure correction (use %COMPCORR if needed).'
      ENDIF 


      IF(OMEGA(N) /= 0.) OMEGAREF = OMEGA(N)  ! To store a reference value

      INPARA(18,N) = INPARA(18,N) - INPARA(1,N) !MIT(N)   = MIT(N) - IMAX(1,N)
      INPARA(20,N) = INPARA(20,N) - INPARA(2,N) !MJT(N)   = MJT(N) - JMAX(1,N)
      INPARA(22,N) = INPARA(22,N) - INPARA(3,N) !MKT(N)   = MKT(N) - KMAX(1,N)
      NDIV         = 2**(LEVEL-1)
      IMAXG(N)     = INPARA(1,N)/NDIV
      JMAXG(N)     = INPARA(2,N)/NDIV
      KMAXG(N)     = MAX0(1,INPARA(3,N)/NDIV)

      INPARA(18,N) = AMIN0(0,INPARA(18,N))
      INPARA(20,N) = AMIN0(0,INPARA(20,N))
      INPARA(22,N) = AMIN0(0,INPARA(22,N))
      INPARA(18,N) = AMAX0(-IMAXG(N),INPARA(18,N))
      INPARA(20,N) = AMAX0(-JMAXG(N),INPARA(20,N))
      INPARA(22,N) = AMAX0(-KMAXG(N),INPARA(22,N))
      INPARA(17,N) = AMAX0(1,INPARA(17,N))
      INPARA(19,N) = AMAX0(1,INPARA(19,N))
      INPARA(21,N) = AMAX0(1,INPARA(21,N))
      INPARA(17,N) = AMIN0(IMAXG(N),INPARA(17,N))
      INPARA(19,N) = AMIN0(JMAXG(N),INPARA(19,N))
      INPARA(21,N) = AMIN0(KMAXG(N),INPARA(21,N))

c      IF(IGRID(N) == 0) IGRID(N) = 5  ! Default value, if not given

      IF(IGRID(N,2) == 9) THEN ! A traditional FINFLO-SHIP treatment
         IBTGR(N)     = IBTGR(N)/NDIV               
         JBTGR(N)     = JBTGR(N)/NDIV               
         KBTGR(N)     = KBTGR(N)/NDIV               
         ITPGR(N)     = ITPGR(N)/NDIV               
         JTPGR(N)     = JTPGR(N)/NDIV               
         KTPGR(N)     = KTPGR(N)/NDIV
      ENDIF

C ... Take information from another block marked with a negative number
          
      MGRIDN       = INPARA(16,N)
      IF(MGRIDN <= -N) STOP 'Negative MGRID'
      IF(MGRIDN <= -1) INPARA(16,N) = INPARA(16,-INPARA(16,N))
C     IF(MGRIDN <= -1) INPARA(16,N) = MGRID(-MGRID(N))

      MGM  = MAX0(MGM,INPARA(16,N)+1) ! maximum multigrid
      WRITE(45,'(2I7,2I9,I10)') N,(INPARA(I,N),I=1,3),
     +     INPARA(1,N)*INPARA(2,N)*INPARA(3,N)
      MTOT = MTOT + INPARA(1,N)*INPARA(2,N)*INPARA(3,N)
      IF(N == NBLOCK) WRITE(45,*) 
     +     ' ====== total in comp. blocks ======',MTOT 
      IF(INPARA(12,N) < INPARA(1,N)+1) 
     + WRITE(*,*) 'Inviscid calculation in I-dir in block and index',N,
     + INPARA(12,N)
      IF(INPARA(13,N) < INPARA(2,N)+1) 
     + WRITE(*,*) 'Inviscid calculation in J-dir in block and index',N,
     + INPARA(13,N)
      IF(INPARA(14,N) < INPARA(3,N)+1) 
     + WRITE(*,*) 'Inviscid calculation in K-dir in block and index',N,
     + INPARA(14,N)

      IF(IGRID(N,1) >= 40 .AND. IGRID(N,1) < 50) THEN
         WRITE(13,*) ' Actuator disk found in block',N
      ENDIF

C ... Set the material property table
      
      DO I = 1,BLKS(N)%NPHASE

      BLKS(N)%ISTATE(I)   = JSTATE(N,I)

      SELECT CASE(JSTATE(N,I))
      CASE(1)
      BLKS(N)%MATERIAL(I) = 'PERF. GAS'
      BLKS(N)%ICHAR(I)    = 2
      CASE(2)
      BLKS(N)%MATERIAL(I) = 'EQUIL. AIR'
      BLKS(N)%ICHAR(I)    = 2
      CASE(3)
      BLKS(N)%MATERIAL(I) = 'PGAS (CP)'
      BLKS(N)%ICHAR(I)    = 2
      CASE(4)
      BLKS(N)%MATERIAL(I) = 'MOD. PGAS'
      BLKS(N)%ICHAR(I)    = 2
      CASE(5)
      BLKS(N)%MATERIAL(I) = 'AIR'
      BLKS(N)%ICHAR(I)    = 2
      CASE(6)
      BLKS(N)%MATERIAL(I) = 'WATER'
      BLKS(N)%ICHAR(I)    = 1
      CASE(7)
      BLKS(N)%MATERIAL(I) = 'PWATER'
      BLKS(N)%ICHAR(I)    = 1
      CASE(8)
      BLKS(N)%MATERIAL(I) = 'WATER'
      BLKS(N)%ICHAR(I)    = 1
      CASE(9)
      BLKS(N)%MATERIAL(I) = 'STEAM'
      BLKS(N)%ICHAR(I)    = 2
      CASE(10)
      BLKS(N)%MATERIAL(I) = 'INCOMPRESSIBLE WATER'
      BLKS(N)%ICHAR(I)    = 1
      CASE(11)
      BLKS(N)%MATERIAL(I) = 'STEAM'
      BLKS(N)%ICHAR(I)    = 1
      CASE(12)
      BLKS(N)%MATERIAL(I) = 'WATER'
      BLKS(N)%ICHAR(I)    = 1
      CASE(13)
      BLKS(N)%MATERIAL(I) = 'STEAM'
      BLKS(N)%ICHAR(I)    = 2
      CASE(101)
      BLKS(N)%MATERIAL(I) = 'CAST IRON'
      BLKS(N)%ICHAR(I)    = 3
      CASE(201)
      BLKS(N)%MATERIAL(I) = 'WALL DIST'
      BLKS(N)%ICHAR(I)    = 3
      END SELECT

      ENDDO     
 1000 CONTINUE

      WRITE(45,*) ' ==================================================='
      WRITE(45,*) '                                    ',MTOT   

C *********************** INPUT FILE READ ******************************

C ... make order of chimera blocks
        
      CALL SUPCHI(NCHIMT,NBLOCK,NCHIM,NCHOR)
      CHIMEL = NCHIM > 0

      IF(CHIMEL) THEN
         IF(INCHIML) THEN
         WRITE(13,'(1X,2A,L1)') ' First-order interpolation is applied',
     +    ' in the Chimera frontier. INCHIML: ',INCHIML
         ELSE
         WRITE(13,'(1X,2A,L1)')' Second-order interpolation is applied',
     +    ' in the Chimera frontier. INCHIML: ',INCHIML
         ENDIF
      ELSE
         IF(INCHIML) THEN
            INCHIML = .FALSE.
            WRITE(13,'(1X,2A)')' Chimera is not used, INCHIML is ',
     +      'put to FALSE'
         ENDIF
      ENDIF
       
      IF(REFLECL) THEN
         WRITE(13,'(1X,2A,L1)') ' Surface velocities are reflected',
     +    ' from the domain. REFLECL: ',REFLECL
      ENDIF
       
      IF(ENTROPY_FIX .AND. IFLUX == 1) THEN
         WRITE(13,'(1X,2A,L1)') ' Entropy fix is made in FLUXKE',
     +    ' (Roe). ENTROPY_FIX: ',ENTROPY_FIX
      ELSEIF(.NOT.ENTROPY_FIX .AND. IFLUX == 1) THEN
         WRITE(13,'(1X,2A,L1)') ' Entropy fix is not made in FLUXKE,',
     +    ' in a case of Roe flux. ENTROPY_FIX: ',ENTROPY_FIX
      ENDIF

      IF(TUR_FRS_SOURCE) THEN
         WRITE(13,'(1X,2A,L1)') ' Turbulence sources are corrected',
     +   ' in order to preserve the free-stream values: ',TUR_FRS_SOURCE
      ELSEIF(.NOT.TUR_FRS_SOURCE) THEN
         WRITE(13,'(1X,2A,L1)') ' Turbulence sources are not corrected',
     +   ' in order to preserve the free-stream values: ',TUR_FRS_SOURCE
      ENDIF
      
      IF(TUR_MULPHL) THEN
         WRITE(13,'(1X,2A,L1)') ' Eddy viscosity is modified',
     +   ' because of the bubble motion: ',               TUR_MULPHL
      ELSEIF(.NOT.TUR_MULPHL .AND. MULPHL) THEN
         WRITE(13,'(1X,2A,L1)') ' Eddy viscosity is not modified',
     +   ' in order to take into account the bubbles: ',  TUR_MULPHL
      ENDIF
      
c ... The first processor reads . Ok to this point

      INQUIRE(FILE=GRIDFI,EXIST=THERE)
      IF (.NOT.THERE) THEN
         WRITE(*,470) GRIDFI
         IF(MASTER) WRITE(4,470) GRIDFI
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF
 470  FORMAT(9X,'***   ERROR   ***'/
     &     9X,'CANNOT OPEN GRID FILE ',A80)

      INQUIRE(FILE=BCFILE,EXIST=THERE)
      IF (.NOT.THERE) THEN
         WRITE(*,475) BCFILE
         IF(MASTER) WRITE(4,475) BCFILE
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF
 475  FORMAT(9X,'***   ERROR   ***'/
     &     9X,'CANNOT OPEN CONNECTIVITY FILE ',A80)

C *** CHECK GRID FILE ***************************************************

      DO 476 N = 1,NBLOCK
         IMAXA = INPARA(1,N)
         JMAXA = INPARA(2,N)
         KMAXA = INPARA(3,N)
         IF(IX(N) /= IMAXA+1 .OR. IY(N) /= JMAXA+1 .OR.
     +        IZ(N) /= KMAXA+1) THEN
            WRITE(*,*) 'Grid dimensions do not match in block',N
            WRITE(*,*) 'IMAX,JMAX,KMAX',IMAXA,JMAXA,KMAXA
            WRITE(*,*) ' IX,   IY,  IZ',IX(N),IY(N),IZ(N)
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
C ... SET UP COORL (MOVING OR DEFORMING GRID)
C ... MOVING TEST
         IF(TIMEL) THEN
            IF(IGRID(N,1) > 0) THEN
               COORL = .TRUE.
               IF (SQRT(GVEX**2 + GVEY**2 + GVEZ**2) >= 1.E-5) THEN
                  WRITE(*,*) 'GVEX,GVEY,GVEZ must be zeros' 
                  WRITE(*,*) 'if IGRID(N,1) is greater than zero and'
                  WRITE(*,*) 'TIMEL is TRUE, exiting...'
                  STOP
               ENDIF
                  
            ENDIF
         ELSEIF(.NOT. TIMEL) THEN
            IF(IGRID(N,1) > 1 .AND. IGRID(N,1) <= 49) THEN 
               COORL = .TRUE.                             
            ENDIF
         ENDIF
C ... DEFORMATION TEST
*         IF(IGRID(N,2) /= 0) THEN
         IF(IGRID(N,2) /= 0 .OR. MPCCIL .OR. MODALFSIL) THEN
            COORL = .TRUE.
         ENDIF
         IF((MPCCIL .OR. MODALFSIL) .AND. IGRID(N,2) == 0) THEN
            IGRID(N,2) = 1
	 ENDIF
C         IF(IGRID(N,1) > 1 .AND. IGRID(N,1) <= 39) COORL = .TRUE.
C         IF(IGRID(N,1) /= 0 .AND. TIMEL)            COORL = .TRUE.
C         IF(IGRID(N,2) /= 0)                        COORL = .TRUE.
 476  CONTINUE

C *** FIND MAXIMUN OF IX,IY,IZ

      NTMAX = 0
      DO 480 N = 1,NBLOCK
         IF((IX(N)*IY(N)*IZ(N)) > NTMAX) THEN
            IXMA = IX(N)
            IYMA = IY(N)
            IZMA = IZ(N)
            NTMAX = IX(N)*IY(N)*IZ(N)
         ENDIF
 480  CONTINUE

      REWIND(41)

C *** READ IN BOUNDARY FILE ********************************************

      OPEN(49,FILE=BCFILE,STATUS='OLD',FORM='FORMATTED')

C ... calculate total number of patches

      NBCS = 0
      NPPV = 0
      READ(49,'(1A80)') LINE1   ! Data description
      READ(49,*) NBLKS          ! Total number of blocks
      DO 13 IB=1,NBLKS          ! Start block loop
         READ(49,'(1A80)') LINE1
         READ(49,*) I, J, K     ! Block size (in cells)
         NPPVS = 0
         DO 12 I=1,6            ! Start face loop
            READ(49,*) IFACE,NOFBCS ! Face number, NOF BC patches
            DO 11 J=1,NOFBCS    ! Start BC patch loop
               READ(49,'(1A80)') LINE1
	       READ(LINE1,*) LINE3
	       IF(LINE3 == 'SOL' .OR. LINE3 == 'ROT' .OR. LINE3 ==  ! BCP
     &         'MOV' .OR. LINE3 == 'CML' .OR. LINE3 == 'ACT') 
     &         READ(49,'(1A80)') LINE1
               IF(LINE3 == 'INL' .OR. LINE3 == 'OUT') THEN
                   READ(LINE1,*) LINE3,ITYPE
                   IF(ITYPE <= 10 .OR. ITYPE > 20 .AND. ITYPE <= 30) 
     &                READ(49,'(1A80)') LINE1
               ENDIF
               NBCS  = NBCS + 1
               NPPVS = NPPVS + 1
 11         CONTINUE            ! End BC Patch loop
 12      CONTINUE               ! End face loop
         NPPV = MAX0(NPPV,NPPVS+1)
 13   CONTINUE                  ! End block loop
      
      ALLOCATE(RBC(IC9,NBCS+1),STAT=IERR1)
      ALLOCATE(IBC(IC9,NBCS+1),STAT=IERR2)
      ALLOCATE(BCAPU (NBCS+1),STAT=IERR3)
      ALLOCATE(BNAPU (NBCS+1),STAT=IERR3)
      ALLOCATE(TBC   (NBCS+1),STAT=IERR4)
      IF(IERR1 > 0 .OR. IERR2 > 0 .OR. IERR3 > 0 .OR. IERR4 > 0)THEN
         WRITE(*,*) ' Not enough memory to read BC file. Patches',NBCS
         WRITE(*,*) ' Exiting...'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF
      ALLOCATE(RBC2(IC9,NBCS+1),STAT=IERR1)
      ALLOCATE(IBC2(IC9,NBCS+1),STAT=IERR2)
      ALLOCATE(BCAPU2(NBCS+1),STAT=IERR3)
      ALLOCATE(BNAPU2(NBCS+1),STAT=IERR3)
      IF(IERR1 > 0 .OR. IERR2 > 0 .OR. IERR3 > 0)THEN
         WRITE(*,*) ' Not enough memory to read BC 2 file. Patches',NBCS
         WRITE(*,*) ' Exiting...'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

      CALL READBC(49,NBGRID,IX,IY,IZ,RBC,NBCS,BCAPU,TBC,IBC,
     +     IPRO,INPARA(16,1),INPAR,GROUP,TIMEL,ITIMES,LEVEL,
     +     FRESUL,COORL,IREPEA,NREPEA,BNAPU)

      IF(FRESUL) THEN

         DO N = 1,NBLOCK ! IGRID must be given
            IF(IGRID(N,2) /= 0) IFOUND = .TRUE.
         ENDDO

      IF(IFSBC == 0) THEN
         WRITE(*,*) ' A free surface was found, but IFSBC was not given'
         WRITE(*,*) ' Exiting...'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF
         
       IF(MOD(JFIRST,NDIV) /= 0.0) WRITE(*,*)
     *'WARNING!! JFIRST IS NOT SUFFICIENTLY EVEN FOR THIS GRID LEVEL' 
      JFIRST       = JFIRST/NDIV                 
       
      IF(.NOT.IFOUND) THEN
         WRITE(*,*) 'Warning: A free surface was found, but',
     &   ' IGRID(N,2)=5..9 was not given for any block.'// ' A surface',
     &   ' will not be deformed.'
         WRITE(*,*) 'This time continuing ...'
c        STOP
      ENDIF

      ENDIF ! FRESUL
        

C ... check if solid walls are included but LAMIN is not correct

      DO I = 1,NBCS
         IF (IBC(1,I) >= 8 .AND. IBC(1,I) <= 10) THEN
            N     = IBC(24,I)
            IF(N <= NBLOCK) THEN
            LAMAPU = INPARA(7,N)
            IF(LAMAPU < 100 .AND. (IBC(3,I) == 1 .OR. IBC(3,I) == 4))
     +           THEN
               WRITE(13,4444) N,I,IBC(1,I),IBC(3,I),INPARA(7,N)
            ENDIF
            IF(LAMAPU >= 100) LAMAPU = LAMAPU - 100
            IF(LAMAPU < 10 .AND. (IBC(3,I) == 2 .OR. IBC(3,I) == 5))
     +           THEN
               WRITE(13,4444) N,I,IBC(1,I),IBC(3,I),INPARA(7,N)
            ENDIF
            IF(LAMAPU >= 10) LAMAPU = LAMAPU - 10
            IF(LAMAPU < 1 .AND. (IBC(3,I) == 3 .OR. IBC(3,I) == 6))
     +           THEN
               WRITE(13,4444) N,I,IBC(1,I),IBC(3,I),INPARA(7,N)
            ENDIF
            ENDIF
         ENDIF
      ENDDO

 4444 FORMAT(//' WARNING:'/
     + 'There is solid, but LAMIN is not described in that direction.'/
     + ' Block, patch number, bc-number, iface and lamin',5I5//)

C ... multiproces directives ============================

      MBPRO = NBLOCG            ! NBLOCG is bit too large but...
      ALLOCATE(NPROCE(MBPRO+1,NPRO))

      IF(PARALLEL)THEN
        
C     Pick the file name from command-line, if one is available

         IOTY = 3
         I0 = IARGC()

         IF (I0 > 0) THEN
            DO I=1,I0
               CALL GETARG(I,LINE1)
               IF (LINE1(1:2) == '-a') THEN ! automatic parallel run
                  IOTY = 4
                  GOTO 111
               ENDIF
            ENDDO
 111        CONTINUE
         ENDIF
            
         IF(IOTY == 4) THEN
            IOTY = 3
            CALL MKPROC(NPROCE,MBPRO,NPRO,IBC,NBCS,NBLOCK,IX,IY,IZ,
     &                  INPARA,INPAR,LEVEL)
         ELSE
            CALL MPRO2(NPROCE,MBPRO)
         ENDIF ! IOTY == 4

      ELSE ! A single-processor run

         NPROCE(1,1) = NBLOCK
         DO N = 1,MBPRO
            NPROCE(N+1,1) = N
         ENDDO
      ENDIF ! Parallel

      IF(MASTER .AND. PARALLEL) THEN
         WRITE(45,*) 'PROCES file read:'
         WRITE(45,*) NPRO,IOTY
         DO N = 1,NPRO
            WRITE(45,*)  NPROCE(1,N)
            WRITE(45,*) (NPROCE(J,N),J=2,NPROCE(1,N)+1)
         ENDDO
      ENDIF

C ... Multi-process directives done

C ... Set order of connections

      CALL CONORD(IBC,NBCS,NPROCE,MBPRO,NPRO,IOTY,ITURB,ISTRES,NSCAL,
     +     INPARA(1,1),INPARA(2,1),INPARA(3,1),INPAR,INPARA(16,1),LEVEL,
     +     MULPHL,NPHASES,TRANSL,TWO_FLUIDL)

      CALL PRIICO(IBC,NBCS)

      IF(PARALLEL) THEN
      DO NP = 2,NPRO
         CALL MPI_SEND(NPROCE(1,NP),1,MPI_INTEGER,NP-1,50,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(NBLOCG,1,MPI_INTEGER,NP-1,51,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(MGM,1,MPI_INTEGER,NP-1,52,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(NPPV,1,MPI_INTEGER,NP-1,53,
     +        MPI_COMM_WORLD,IERR)
      ENDDO

      CALL MPITBC(IBC,RBC,NBCS,IBC2,RBC2,NPROCE,MBPRO,NPRO,
     +            NBLOCK,BCAPU,IOTY,BCAPU2,BNAPU,BNAPU2)

      CALL MPITIN(NBLOCK,INPARA,INPAR,NPROCE,MBPRO,NPRO)


      ENDIF ! PARALLEL    

C *********************** WRITE BEGINING OF CALCULATIONS ***************
       
      CALL HEADING(6)
      CALL HEADING(13)

      WRITE( 6,9946) NAME(1:60)
      WRITE(44,9946) NAME(1:60)
      WRITE(45,9946) NAME(1:60)

      IF(IOLD == 0 ) THEN
         IF(MASTER)  WRITE(4,9946)  NAME(1:60)
         WRITE(13,9950) NAME(1:60)
         WRITE(13,*)
      ELSE IF(IOLD > 0) THEN
         IF(MASTER)  WRITE(4,9946)  NAME(1:60)
         WRITE(13,9951) NAME(1:60)
         WRITE(13,*)
      ELSE IF(IOLD < 0) THEN
         IF(MASTER)  WRITE(4,9946)  NAME(1:60)
         WRITE(13,9952) NAME(1:60)
         WRITE(4,9953)
         WRITE(13,*)
      ENDIF

9946  FORMAT('  SIMULATION OF THE ',A60)
9950  FORMAT('  I am starting a simulation of ',A60)
9951  FORMAT('  I am continuing the simulation of ',A60)
9952  FORMAT('  I am starting a simulation of ',A60)
9953  FORMAT('  Initializing from the previous level')
      WRITE(13,'(A,I4)') '  Grid level = ',LEVEL
      WRITE(45,'(A,I4)') '  Grid level = ',LEVEL

      TRANSIT = IFLUX == 1 .OR. IFLUX == 4 .OR. IFLUX == 8 ! Check transition

      IF(IFLUX == 0) THEN ! VAN LEER
         WRITE( 4,9150)
         WRITE( 6,9150)
         WRITE(45,9150)
         WRITE(13,9150)
      ELSEIF(IFLUX == 1) THEN ! ROE
         WRITE( 4,9151)
         WRITE( 6,9151)
         WRITE(45,9151)
         WRITE(13,9151)
      ELSE IF(IFLUX == 2) THEN !  WAVE/PARTICLE SPLITTING SCHEME
         WRITE( 4,9152)
         WRITE( 6,9152)
         WRITE(45,9152)
         WRITE(13,9152)
      ELSE IF(IFLUX == 3) THEN 
C ... WAVE/PARTICLE SPLITTING SCHEME USING AUSM IMPLEMENTATION
         WRITE( 4,9153)
         WRITE( 6,9153)
         WRITE(45,9153)
         WRITE(13,9153)
      ELSE IF(IFLUX == 4) THEN 
C ... APPROXIMATIVE INCOMPRESSIBLE FLUX. PSEUDOCOMPRESSIBILITY
         WRITE( 4,9154)
         WRITE( 6,9154)
         WRITE(45,9154)
         WRITE(13,9154)
      ELSE IF(IFLUX == 5) THEN ! APPROXIMATIVE Roe for RSM
         WRITE( 4,9155)
         WRITE( 6,9155)
         WRITE(45,9155)
         WRITE(13,9155)
         IFLUX = 5
         IF(ITURB < 20) THEN
            WRITE(*,*) 'You have not activated RSM. You cannot have '
     +           ,'flux-splitting ROERSM'
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF
      ELSE IF(IFLUX == 6) THEN 
C ... MULTI-PHASE FLOW
         WRITE( 4,9156)
         WRITE( 6,9156)
         WRITE(45,9156)
         WRITE(13,9156)
      ELSEIF(IFLUX == 7) THEN ! Kurganov-Tadmor
         WRITE( 4,9157)
         WRITE( 6,9157)
         WRITE(45,9157)
         WRITE(13,9157)
      ELSEIF(IFLUX == 8) THEN ! H-CUSP
         WRITE( 4,9158)
         WRITE( 6,9158)
         WRITE(45,9158)
         WRITE(13,9158)
      ENDIF
      IF(.NOT.TRANSIT .AND. TRANSL) THEN
         WRITE( 4,9159)
         WRITE( 6,9159)
         WRITE(45,9159)
         WRITE(13,9159)
            WRITE(*,*) 'Exiting ...'
            STOP
         ENDIF

 9150 FORMAT('  USING FLUX-VECTOR SPLITTING OF VAN LEER FOR INVISCID ',
     +     'FLUXES')
 9151 FORMAT('  USING FLUX-DIFFERENCE SPLITTING OF ROE FOR INVISCID ',
     +     'FLUXES')
 9152 FORMAT('  USING WAVE/PARTICLE SPLITTING SCHEME FOR INVISCID ',
     +     'FLUXES')
 9153 FORMAT('  USING WPS/AUSM SCHEME FOR INVISCID FLUXES')
 9154 FORMAT('  USING MUSCL-BASED SCHEME FOR INVISCID FLUXES ')
 9155 FORMAT('  USING FLUX-DIFFERENCE SPLITTING OF ROE FOR INVISCID ',
     +     'FLUXES (RSM)')
 9156 FORMAT('  USING MULTI-PHASE FLUXES FOR INVISCID ',
     +     'FLUXES')
 9157 FORMAT('  USING RUSANOV (KURGANOV-TADMOR) SCHEME FOR INVISCID ',
     +     'FLUXES')
 9158 FORMAT('  USING H-CUSP SCHEME FOR INVISCID ',
     +     'FLUXES')
 9159 FORMAT('  TRANSITION MODEL WORKS CURRENTLY ONLY WITH ROE, H-CUSP',
     +     ' AND',/,'  INCOMPRESSIBLE FLUX SCHEMES')
      IF(IPRO == 1)WRITE(4,*)
      WRITE(45,*)
      WRITE(13,*)
      IF (TIMEL) THEN
         IF(MASTER) WRITE(4, 7107)
         WRITE(45,7107)
         WRITE(13,7107)
      ELSE
         IF(MASTER) WRITE(4, 7117)
         WRITE(45,7117)
         WRITE(13,7117)
      ENDIF ! TIMEL
 7107 FORMAT('  Time-accurate, saving old timelevels')
 7117 FORMAT('  Pseudo time integration')
      IF (COORL) THEN
         IF(MASTER) WRITE(4, 7127)
         WRITE(45,7127)
         WRITE(13,7127)
      ELSE
         IF(MASTER) WRITE(4, 7137)
         WRITE(45,7137)
         WRITE(13,7137)
      ENDIF ! COORL
      IF (MULPHL .AND. BLKS(1)%SOLUTION_TYPE == 'CAVIT') THEN
         IF(MASTER) WRITE(6, 7128)
c         WRITE(45,7128)
c         WRITE(13,7128)
      ELSEIF (MULPHL .AND. BLKS(1)%SOLUTION_TYPE == 'MULTI') THEN
         IF(MASTER) WRITE(6, 7138)
c         WRITE(45,7138)
c         WRITE(13,7138)
      ENDIF ! CAVIT

 7127 FORMAT('  Calculation of changing grid')
 7137 FORMAT('  Grid is unchanged during calculation')
 7128 FORMAT('  TWO-PHASE FLOW WITH EQUAL VELOCITIES IS ASSUMED')
 7138 FORMAT('  TWO-PHASE FLOW WITH NON-EQUAL VELOCITIES IS ASSUMED')
      
      IF(IPRESC == 2) THEN
         IF(MASTER) WRITE(4, 7147)
         WRITE(45,7147)
         WRITE(13,7147)
      ELSE IF(IPRESC == 1) THEN          
         IF(MASTER) WRITE(4, 7157)
         WRITE(45,7157)
         WRITE(13,7157)
      ENDIF ! IPRESC == 2
 7147 FORMAT('  Pseudocompressibility approach is utilized')
 7157 FORMAT('  Pressure correction method')

      IF (MULPHL .AND. .NOT.TWO_FLUIDL) THEN
         IF(MASTER) WRITE(4, 7166)
         WRITE(45,7166)
         WRITE(13,7166)
      ELSEIF (MULPHL .AND. TWO_FLUIDL) THEN
         IF(MASTER) WRITE(4, 7266)
         WRITE(45,7266)
         WRITE(13,7266)
      ELSE
         IF(MASTER) WRITE(4, 7177)
         WRITE(45,7177)
         WRITE(13,7177)
      ENDIF ! MULPHL
 7166 FORMAT('  Two-phase flow is enabled')
 7266 FORMAT('  Two-phase flow is enabled with unequal velocities')
 7177 FORMAT('  A single-phase flow is assumed')
        
      IF (XXTRAL) THEN
         IF(MASTER) WRITE(4, 7186)
         WRITE(45,7186)
         WRITE(13,7186)
      ELSE
         IF(MASTER) WRITE(4, 7187)
         WRITE(45,7187)
         WRITE(13,7187)
      ENDIF ! XXTRAL
 7186 FORMAT('  A distance-based MUSCL interpolation is assumed')
 7187 FORMAT('  A distance-based MUSCL interpolation is ignored')
   
      IF (FULLNS) THEN
         IF(MASTER) WRITE(4, 7188)
         WRITE(45,7188)
         WRITE(13,7188)
      ELSE
         IF(MASTER) WRITE(4, 7189)
         WRITE(45,7189)
         WRITE(13,7189)
      ENDIF ! FULNS
 7188 FORMAT('  A full friction term is applied according to IDIFF')
 7189 FORMAT('  A thin-layer aprroximation is used for friction')
   
      IF (FRESUL) THEN
         IPRIFS  = 0
         IF(MASTER) WRITE(4, 7196)
         WRITE(45,7196)
         WRITE(13,7196)
         IF(IFSBC == 0) THEN
           WRITE( 4,7198)
           WRITE(13,7198)
           WRITE(45,7198)
           WRITE(*,7198)
           STOP
         ELSE IF(IFSBC == 1) THEN
c          WRITE( 4,7199) IFSBC,JFIRST
           WRITE(45,7199) IFSBC,JFIRST
         ELSE
c          WRITE( 4,7200) IFSBC
           WRITE(45,7200) IFSBC
         ENDIF
      ELSE
         IF(MASTER) WRITE(4, 7197)
         WRITE(13,7197)
         WRITE(45,7197)
      ENDIF ! FRESUL
 7196 FORMAT('  A free-surface model is enabled')
 7197 FORMAT('  A free-surface model is disabled')
 7198 FORMAT('  No free-surface option IFSBC specified. Exiting...')
 7199 FORMAT('  IFSBC = ',I2,'  JFIRST = ',I2)
 7200 FORMAT('  IFSBC = ',I2)
C ***************** SET ARTIFICIAL COMPRESSIBILITY *******************

      PSEUCO  = 10.

C ***************** SET THE EQUATION OF STATE  ***********************

      IF(ISTATE == 1) THEN
         IF(MASTER) WRITE(4,7167)
      ELSE IF(ISTATE == 2) THEN
         IF(MASTER) WRITE(4,7168)
      ELSE IF(ISTATE == 3) THEN
         IF(MASTER) WRITE(4,7167)
      ELSE IF(ISTATE == 4) THEN
         IF(MASTER) WRITE(4,7170)
      ELSE IF(ISTATE == 5) THEN
         WRITE(*,*)' WRONG EQUATION OF STATE'
         WRITE(*,*) 'ISTATE = 5 not in use anymore.'
         WRITE(*,*) 'Use ISTATE = 1.  Exiting ...'
*        CALL EXIT
         STOP
      ELSE IF(ISTATE == 6) THEN
         IF(MASTER) WRITE(4,7169)
         IF(IPRESC == 0) THEN
          WRITE(*,*)' Do not possible work without pseudo compressible'
         ENDIF
      ELSE IF(ISTATE == 7) THEN
         IF(MASTER) WRITE(4,7173)
         IF(IPRESC == 0) THEN
          WRITE(*,*)' Do not possible work without pseudo compressible'
         ENDIF
      ELSE IF(ISTATE == 8) THEN
         IF(MASTER) WRITE(4,7174)
         IF(IPRESC == 0) THEN
          WRITE(*,*)' Do not possible work without pseudo compressible'
         ENDIF
      ELSE IF(ISTATE == 9) THEN
         IF(MASTER) WRITE(4,7175)
      ELSE IF(ISTATE == 10) THEN
         IF(MASTER) WRITE(4,7178)
      ELSE IF(ISTATE == 11) THEN
         IF(MASTER) WRITE(4,7178) 
      ELSE IF(ISTATE == 12) THEN
         IF(MASTER) WRITE(4,7179)
      ELSE IF(ISTATE == 13) THEN
         IF(MASTER) WRITE(4,7180)

      ELSE IF(ISTATE == 0) THEN
         WRITE(4,7176) 
         DO IPHASE = 1,NPHASE
         WRITE(4,*) '   N   NGL JSTATE'  
         WRITE(4,*) ' ================================================',
     2              '=================='
         CALL STATE_REPORT(NBLOCK,NPROCE,NBLOCG,NPRO,JSTATE,IPRO,
     2   MULPHL,IPHASE,IPRESC)    ! In outp3.f
         WRITE(4,*) ' ================================================',
     2              '=================='
         WRITE(4,*)
         ENDDO
      ELSE
         WRITE(4,*) ' No such equation of state'
         WRITE(*,*) ' No such equation of state. Exiting ...'
      ENDIF                     ! ISTATE == 1

 7167 FORMAT('  PERFECT GAS AS A FUNCTION OF RHO AND E')
 7168 FORMAT('  CHEMICALLY REACTING EQUILIBRIUM AIR')
 7169 FORMAT('  WATER AS A FUNCTION OF P AND T')
 7170 FORMAT('  SEMI PERFECT GAS')
 7171 FORMAT('  PERFECT GAS A FUNCTION OF P AND E')
 7172 FORMAT('  EQUATION OF STATE DOES NOT EXIST'//'  Exiting ...')
 7173 FORMAT('  WATER AS A FUNCTION OF P AND T. MODIFIED VISCOSITY!')
 7174 FORMAT('  WATER AS A FUNCTION OF P AND T. RANGE 0.5 - 1.5 BAR')
 7175 FORMAT('  STEAM AS A FUNCTION OF P AND T. RANGE 0.5 - 1.5 BAR')
 7178 FORMAT('  INCOMPRESSIBLE FLOW BASED ON THE FREE-STREAM VALUES')
 7179 FORMAT('  CONSTANT PROPERTIES FOR R-12 LIQUID')
 7180 FORMAT('  CONSTANT PROPERTIES FOR R-12 GAS')

 7176 FORMAT(/'  EQUATION OF STATE IS DEFINED BLOCKWISE'/)
      
C ***************** SET THE TURBULENCE MODEL   ***********************

      SELECT CASE(ITURB)
      CASE(:0)
         RIVI='ASSUMING COMPLETELY LAMINAR FLOW'
      CASE(1)
         RIVI='EMPLOYING ALGEBRAIC BALDWIN-LOMAX TURBULENCE MODELS'
      CASE(2)
         RIVI='EMPLOYING ALGEBRAIC CEBECI-SMITH TURBULENCE MODELS'
      CASE(3)
         RIVI='EMPLOYING THE K-EPSILON TURBULENCE MODEL'
      CASE(4)
         RIVI='THIS TURBULENCE MODEL IS NOT USED (ITURB=4)'
      CASE(5)
         RIVI='EMPLOYING THE RNG-K-EPSILON TURBULENCE MODEL '//'
     +        (warning: not tested)'
      CASE(6)
       IF(ISTRES == 0 .AND. .NOT.TURDESL) THEN
         RIVI='EMPLOYING THE SST K-OMEGA TURBULENCE MODEL'
         IF ((KOVER >= 1 .AND. KOVER <= 4) .OR. KOVER >= 6) THEN
            RIVI = 'EMPLOYING THE K-OMEGA-RCSST TURBULENCE MODEL'
         ELSE IF (TRANSL .AND. .NOT.LONG_TRANSL) THEN
            RIVI = 'EMPLOYING THE SST K-OMEGA MODEL WITH TRANSITION'
         ELSE IF (TRANSL .AND. LONG_TRANSL) THEN
            RIVI = 'EMPLOYING THE SST MODEL AND STARTING TRANSITION'
         ELSE IF(KOVER == 5) THEN
            RIVI='EMPLOYING THE K-OMEGA-BSL TURBULENCE MODEL'
         END IF
       ELSE IF(ISTRES == 1 .AND. .NOT.TURDESL) THEN
         RIVI='EMPLOYING THE SST K-OMEGA TURBULENCE MODEL WITH EARSM'
         IF ((KOVER >= 1 .AND. KOVER <= 4) .OR. KOVER >= 6) THEN
            RIVI = 'EMPLOYING K-OMEGA-RCSST TURBULENCE MODEL WITH EARSM'
         ELSE IF(KOVER == 5) THEN
            RIVI='EMPLOYING THE K-OMEGA-BSL TURBULENCE MODEL WITH EARSM'
         END IF
       ELSE IF(ISTRES == 0 .AND. TURDESL) THEN
         RIVI='EMPLOYING DES BASED ON THE SST K-OMEGA TURBULENCE MODEL'
         IF (KOVER == 2) THEN
            RIVI = 'EMPLOYING DDES BASED ON THE SST K-OMEGA MODEL'
         ELSE IF(KOVER == 3) THEN
            RIVI='EMPLOYING SAS BASED ON THE SST K-OMEGA MODEL'
         ELSE IF(KOVER == 5) THEN
            RIVI='EMPLOYING DES BASED ON THE K-OMEGA-BSL MODEL'
         ELSE IF(KOVER == 6) THEN
            RIVI='EMPLOYING DDES BASED ON THE K-OMEGA-BSL MODEL'
         ELSE IF(KOVER == 8) THEN
            RIVI='EMPLOYING SAS BASED ON THE K-OMEGA-BSL MODEL'
         ELSE IF(KOVER /= 0) THEN
            RIVI='WRONG ITURB VERSION, SORRY ABOUT THAT CHIEF'
         ENDIF
       ELSE IF(ISTRES /= 0 .AND. TURDESL) THEN
         WRITE(*,*) ' Dubious ISTRES option with DES. Warning...'
         WRITE(13,*)' Dubious ISTRES option found with DES. Warning...'
c         STOP
       ENDIF
      CASE(7)
         RIVI='EMPLOYING THE K-EPSILON MODEL WITH CROSSDIFFUSION'
      CASE(8)
         RIVI='EMPLOYING THE SMAGORINSKY SGS-model for LES'
      CASE(9)
         RIVI='EMPLOYING THE SPALART-ALLMARAS TURBULENCE MODEL'
         IF ((KOVER >= 1 .AND. KOVER <= 4) .AND. TURDESL) THEN
            RIVI = 'EMPLOYING DES (S-A) WITH DDES OPTION (n.a.)'
         ELSE IF(KOVER == 5) THEN
            RIVI='UNUSED'
         END IF
      CASE(10)
         RIVI='EMPLOYING THE K-EPSILON MODEL WITH AERMS'
      CASE(11)
         RIVI='EMPLOYING THE K-EPSILON MODEL (of Chien) WITH AERMS'
      CASE(12:19)
         RIVI='INVALID TURBULENCE MODEL SELECTED'
      CASE(20:23)
         RIVI='EMPLOYING THE UNCOUPLED REYNOLDS STRESS MODEL'
      CASE(24:)
         IF(JRPRE == 1) THEN
            RIVI='EMPLOYING SHIMAS REYNOLDS STRESS MODEL'
         ELSEIF(JRPRE == 2) THEN
            RIVI='EMPLOYING SSG/SHIMAS REYNOLDS STRESS MODEL'
         ENDIF
      END SELECT

      IF(TRANSL .AND. ITURB /= 6) THEN	!Intermittency variables
         WRITE(*,*) ' Only ITURB=6 allowed with TRANSC. Exiting...'
         WRITE(13,*)' ONLY ITURB=6 allowed with TRANSC. Exiting...'
         STOP
      ENDIF

      IF (MASTER) THEN       ! myid=0 Can only write to the stdout
         WRITE(*,'(2X,A)') TRIM(RIVI)
         IF(TURDESL) WRITE(*,'(2X,A)') 'with detached eddy simulation'
      END IF
      WRITE(4,'(2X,A)') TRIM(RIVI)
      WRITE(13,'(2X,A)') TRIM(RIVI)
      IF (TRANSL .AND. LONG_TRANSL) WRITE(13,'(2X,2A)') 'Warning:',
     +               ' Gamma and Re_t are initialized'
      IF(TURDESL) WRITE(4,'(2X,A)')'with detached eddy simulation (DES)'

      IF(TURCOR .AND. ITURB >= 3) THEN
          WRITE(4,*) ' TURBULENCE VARIABLE CORRECTIONS ON THE WALLS',
     +               ' ENABLED'
      ELSE
          WRITE(4,*) ' TURBULENCE VARIABLE CORRECTIONS ON THE WALLS',
     +               ' DISABLED'
      ENDIF

      IF(ITURB /= 0) THEN

      IF(INTERTU /= 0) THEN
          WRITE(4,*) ' TURBULENCE VARIABLES INTERPOLATED ACCORDING',
     +               ' TO INTER=',INTERTU
      ELSE
          WRITE(4,*) ' TURBULENCE VARIABLES INTERPOLATED AS INDICATED',
     +               ' BY MAIN VARIABLES'
      ENDIF
      ENDIF ! ITRUB /= 0

      IF(IOLD == 0) THEN
        WRITE(*,9947) LEVEL
      ELSE IF(IOLD > 0 .AND. LONG_STARTL) THEN
        WRITE(*,9948) LEVEL
      ELSE IF(IOLD > 0 .AND. .NOT.LONG_STARTL) THEN
        WRITE(*,9941) LEVEL
        WRITE(13,9942)
      ELSE
        WRITE(*,9949) LEVEL
      ENDIF ! IOLD == 0

9947  FORMAT('  ON LEVEL',I2,' BEGINS')
9948  FORMAT('  ON LEVEL',I2,' CONTINUES')
9949  FORMAT('  ON LEVEL',I2,' BEGINS. INITIALIZATION PERFORMED FROM ',
     +     'THE PREVIOUS GRID LEVEL')
9941  FORMAT('  ON LEVEL',I2,' CONTINUES. A SHORT RESTART ',
     +     'IS ASSUMED (MASS FLOWS ARE MISSED)')
9942  FORMAT('  A short restart list is assumed. Mass flows are not ',
     +     'read in. Are you ','sure about this? :-)')

      WRITE(*,*)'*****************************************************',
     +     '********************'

      ENDIF  ! (MASTER)

C ******** RECEIVE THE INPUT AND BC DATA FROM PROCES 0 ***************

      IF(PARALLEL) THEN 
         IF(.NOT.MASTER) THEN
         CALL MPI_RECV(NBLOCK,1,MPI_INTEGER,0,50,MPI_COMM_WORLD,STATUS,
     +     IERR)
         CALL MPI_RECV(NBLOCG,1,MPI_INTEGER,0,51,MPI_COMM_WORLD,STATUS,
     +     IERR)
         CALL MPI_RECV(MGM,1,MPI_INTEGER,0,52,MPI_COMM_WORLD,STATUS,
     +     IERR)
         CALL MPI_RECV(NPPV,1,MPI_INTEGER,0,53,MPI_COMM_WORLD,STATUS,
     +     IERR)
         CALL MPI_RECV(NBCS,1,MPI_INTEGER,0,100+IPRO,MPI_COMM_WORLD,
     +        STATUS,IERR)
         ELSE
            NBLOCK = NPROCE(1,1)
         ENDIF

      ENDIF

C ... Allocation for other processes (2... NBGG)
        
*      CALL ALLONS(NBLOCK,NBLOCG,MGM,NBCS,IPRO,NPPV,NPPN,NB,MBPRO,NPRO,
*     +     NBGG)
      CALL ALLONS(NPPN)

      IF(.NOT.MASTER) THEN
         CALL MPI_RECV(ICON,IC9*NBCS,MPI_INTEGER,0,300+IPRO,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(RCON,IC9*NBCS,MPI_REAL8,0,1300+IPRO,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(BOUNDF(1),80*NBCS,MPI_CHARACTER,0,500+IPRO,
     +        MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(BOUNDN(1),80*NBCS,MPI_CHARACTER,0,1500+IPRO,
     +        MPI_COMM_WORLD,STATUS,IERR)
      ELSE

C ... Save the patch data into the ICON vector

         DO 25 IP = 1,NBCS
            DO I=1,IC9
               ICON(I+(IP-1)*IC9)  = IBC(I,IP)
               RCON(I+(IP-1)*IC9)  = RBC(I,IP)
            ENDDO
            BOUNDF(IP)(1:80)  = ' '
            BOUNDN(IP)(1:80)  = ' '
            BOUNDF(IP) = BCAPU(IP)
            BOUNDN(IP) = BNAPU(IP)
 25      CONTINUE
         DEALLOCATE(RBC,IBC,BCAPU,TBC,RBC2,IBC2,BCAPU2,BNAPU,BNAPU2)
      ENDIF

      IF(PARALLEL) THEN 
      
         IF(.NOT.MASTER) ALLOCATE(INPARA(1,1))
         CALL MPITOU(NBLOCK,IPRO,NPROCE,MBPRO,NPRO,IMAX,JMAX,KMAX,
     +        INTERI,INTERJ,INTERK,IDER,LAMIN,INITC,IT,IL,IK,IDI1,IDI2,
     +        IDI3,MOV,MGRID,MIB,MIT,MJB,MJT,MKB,MKT,IROTVE,MGM,
     *        INPAR,INPARA)
	     
         CALL MPI_BCAST(NPROCE,(MBPRO+1)*NPRO,MPI_INTEGER,0,
     +        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NAME,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(OUTVAR,NOVARMAX,MPI_INTEGER,0,
     +                  MPI_COMM_WORLD,IERR)
	       
         CALL MPI_BCAST(CFL   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CFLL  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)  
         CALL MPI_BCAST(DROLIM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TMAX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR) 
         CALL MPI_BCAST(TROT  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR) 
         CALL MPI_BCAST(DTB   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RMACH ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ALPHA ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BETA  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RE    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PR    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PRT   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RGAS  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TSU0  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(VISU0 ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(EXPSU ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GAMMA ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSTEM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSVEL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSDEN,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(EPSLIM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RKLIM ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TURLIM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TURBLE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RMUINI,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CMGK  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CMGEPS,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CSIMPS,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CC1   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CC2   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(AREF  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CHLREF,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GRILEN,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XMOM  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YMOM  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZMOM  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(REFPRE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DIFPRE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSEUCO,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GROUND,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GX    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GY    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GZ    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(G0    ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ALTITUDE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CAVLEV,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSPRE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TOLER ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TEMINI,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SIEINI,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TLOLIM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TUPLIM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRSALF,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SMAX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CDIFF ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CDIFFT,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(IOLD  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(LEVEL ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITURB ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KSCAL ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFLUX ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ICMAX ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KP    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MPRINT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ISTATE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JRDIF ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JRDIS ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JRPRE ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JRIMP ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRESC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IDRXX ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(LUSGS ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IEPSMA,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITERMA,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NOB   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NBLOCG,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NCHIM ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KOVER ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IDIS  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NPHASE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MCYCLE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITERHA,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FEMXYZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RMESHT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(IROTCO ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KRRINT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NSPTOT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ISTRES ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NSKIN  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INTERTU,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ICONV3 ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
	      
         CALL MPI_BCAST(STARTL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(STRESL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FULLNS     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SOURL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(STATEL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TIMEL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(COORL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PARALLEL   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PRESL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(REANEW     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GRAVIL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CONVL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SPLIT      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CHIMEL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GROUP      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PERCHL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MULPHL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TURCOR     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XXTRAL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FRESUL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WOODYL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(LDISTL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CAVITL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TRUE_DISTL ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FPRINTL    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INOUTL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TURDESL    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INCHIML    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(REFLECL    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ENTROPY_FIX,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(SURFACE_FORCES,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TUR_FRS_SOURCE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TRANSL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TUR_MULPHL ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NEGVL      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MPCCIL     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MODALFSIL  ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(LONG_STARTL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(LONG_TRANSL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(AVER_STARTL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TWO_FLUIDL ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WALLFUNL   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(JFIRST,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NFSD  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ICFST ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INWH  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFSBC ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INDXTS,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IGLOBL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NGLOBL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRIFS,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(FROUDE ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DTWMAX ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR) 
         CALL MPI_BCAST(DWMV   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DWMAX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FLODWH ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WHMAX  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WHMIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FREDIF ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SUMDWH ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CFLFRE ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DTWMIN ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(AGAMMA ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ABANK  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ADVANJ ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(REFVEL ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTORW ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WHEIGHT,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XBULB  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YBULB  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ABULB  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QGFIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QGEIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QGFOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QGEOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QLFIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QLEIN  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QLFOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QLEOUT ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NEGV   ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(SOLTEM,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ARTSSP,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RMULTV,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ALFAP ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ALFAU ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RJK2  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RJK4  ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(MCYCAM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITERAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITERAC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NUMCOL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SMAX,  1,MPI_REAL8,   0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ALFAP, 1,MPI_REAL8,   0,MPI_COMM_WORLD,IERR)
       
         CALL MPI_BCAST(GVEX,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GVEY,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GVEZ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)         
         
c         CALL MPI_BCAST(STARTL,16,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
C ...                          ** huom obs !
         CALL MPI_BCAST(IMAXG, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JMAXG, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KMAXG, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(IUPPT, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         
         IF(FRESUL) THEN
         CALL MPI_BCAST(IBTGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JBTGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KBTGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ITPGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JTPGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KTPGR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ICNH,  NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GML,   NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         ENDIF
	     
         CALL MPI_BCAST(OMEGA, NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OMEGAX,NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OMEGAY,NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OMEGAZ,NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
	        
         CALL MPI_BCAST(CENAX, NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CENAY, NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CENAZ, NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(AMPL,  NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSCLS, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTAT ,NBLOCG,MPI_REAL8,  0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(BLKS%SOLUTION_TYPE,10*NBLOCG,MPI_CHARACTER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%MATERIAL(1),  10*NBLOCG,MPI_CHARACTER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%MATERIAL(2),  10*NBLOCG,MPI_CHARACTER,0,
     +   MPI_COMM_WORLD,IERR)
        IF(NPHASES == 3) THEN
        CALL MPI_BCAST(BLKS%MATERIAL(NPHASES),10*NBLOCG,MPI_CHARACTER,0,
     +   MPI_COMM_WORLD,IERR)
        ENDIF

         CALL MPI_BCAST(BLKS%NPHASE,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         
         DO IPHASE = 1,NPHASE
            CALL MPI_BCAST(BLKS%ISTATE(IPHASE),NBLOCG,MPI_INTEGER,0,
     +      MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(BLKS%IUPTEM(IPHASE),NBLOCG,MPI_INTEGER,0,
     +      MPI_COMM_WORLD,IERR)
         ENDDO
          
         CALL MPI_BCAST(BLKS%IPRESC, NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%LUSGS,  NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%IFLUX,  NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%CONVL,NBLOCG,MPI_LOGICAL,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSDEN, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSPRE, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSTEM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%TEMINI, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSVEL, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%ARTSSP, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%RMULTV, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%SOLTEM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSTUR, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSMUT, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%TURLIM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%RE,     NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%PR,     NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%PRT,    NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%CHLREF, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%RMACH,  NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%RMUINI, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%AREF,   NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%REFPRE, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%REFVEL, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%DIFPRE, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSSSP, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%FRSLEN, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%RKLIM,  NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%EPSLIM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%TUPLIM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%TLOLIM, NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%ALFAP,  NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%ALFAU,  NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%ALFASKEW,NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%GVEX,   NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%GVEY,   NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(BLKS%GVEZ,   NBLOCG,MPI_REAL8,0,MPI_COMM_WORLD,
     +   IERR)         
         
         CALL MPI_BCAST(BLKS%SKEWCORR, NBLOCG,MPI_LOGICAL,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%COMPCORR, NBLOCG,MPI_LOGICAL,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%KATOL,    NBLOCG,MPI_LOGICAL,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%ZEROV,    NBLOCG,MPI_LOGICAL,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FLUXCORR, NBLOCG,MPI_LOGICAL,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%EVAP_TYPE,NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%CAVNO,    NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%IDIFF,    NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%MPCASE,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%PREVEL,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%INLRC,    NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%OUTRC,    NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%ICERO,    NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%ITERMP,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FRSG,     NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FRSRET,   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%CHLCAV,   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%REFCAV,   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%VOIDTT,   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FREDIFI,  NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FREDIFJ,  NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FREDIFK,  NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)         
         CALL MPI_BCAST(BLKS%INTERI,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%INTERJ,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%INTERK,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%TURMUL,   NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%DRAGMUL,  NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%LIFTMUL,  NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%DISPERMUL,NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%VMASS,    NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%WFORCEMUL,NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)

         
         DO IPHASE = 1,NPHASES
         CALL MPI_BCAST(BLKS%FRSALFA(IPHASE),NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FRSX(IPHASE),   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%TAUF(IPHASE),   NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%N(IPHASE),      NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%RK(IPHASE),     NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%ICHAR(IPHASE),  NBLOCG,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(BLKS%FRADEN(IPHASE), NBLOCG,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         ENDDO
         
         CALL MPI_BCAST(NGRIFL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         IF(NGRIFL > 0) THEN
         CALL MPI_BCAST(TRMODE,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%RMASS,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%RCG,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%AREF,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%CHLREF,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IX,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IY,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IZ,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IXY,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IXZ,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%IYZ,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%TDEL,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%TSTEER,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%VX,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%VY,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%VZ,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%ROTX,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%ROTY,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%ROTZ,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%DAMPN,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(OSKU%DAMPT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCGI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YCGI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCGI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCG,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YCG,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCG,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIR,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETAR,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIR,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIRI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETARI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PHIRI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(PSIM,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XCGIS,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YCGIS,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZCGIS,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DAMPC1,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DAMPC2,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DRAUGHTI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DRAUGHT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SINKI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SINK,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TRIMAI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TRIMA,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RIN,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROUT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THRUST,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(TORQUE,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ADV,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTA1,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTB1,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTA1I,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTB1I,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CONEA,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CONEAI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ACTUA,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(VTIP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CDBLADE,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CFBLADE,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SIGMA,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CFTI,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFA,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFT,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFR,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NSERIES,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IGRSLAVE,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTA1S,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ROTB1S,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RTMSP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FDSP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(UTSP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FXSP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FXTSP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(QFACT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RGML,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XFAKEP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YFAKEP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZFAKEP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)         
         CALL MPI_BCAST(ROTB1FAKEP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NBLADE,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IFAKEP,NGRIFL+10,MPI_INTEGER,0,
     +   MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(SHAFT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THETACOL,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THCYCLON,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(THCYCLAT,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ETIP,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CBLADE,NGRIFL+10,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(OSKU%CHAR,NGRIFL+10,MPI_CHARACTER,0,
     +   MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(MVSHIP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SHIP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SHIPWAVE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ACTDISK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(HROTOR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(HMVBLADE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FLYOBJ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(SHIPPROP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

         CALL MPI_BCAST(RCG,4,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
c         CALL MPI_BCAST(RM,4,MPI_REAL8,0,
c     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RITH,4,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RIZE,4,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(RIBE,4,MPI_REAL8,0,
     +   MPI_COMM_WORLD,IERR)
         ENDIF

         CALL MPI_BCAST(IGRID(:,1),NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(IGRID(:,2),NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(IGRID(:,3),NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(IGRID(:,4),NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(IGRID(:,5),NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
c         CALL MPI_BCAST(JSTATE,3*NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
c     +   IERR)
         CALL MPI_BCAST(JSTATE,3*NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,
     +   IERR)
         CALL MPI_BCAST(NCHIMT,NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NCHOR, NBLOCG,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IDOB,NOB,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(JDOB,NOB,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KDOB,NOB,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NSSB,NOB,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ISSB,NBGG*7,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
           
C ... Force group info
         CALL MPI_BCAST(FORCE_GROUP_SHORT_NAME,52*1,MPI_CHARACTER,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(FORCE_GROUP_FULL_NAME,52*24,MPI_CHARACTER,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XMOMR,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YMOMR,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZMOMR,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(XMOMA,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YMOMA,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(ZMOMA,52,MPI_REAL8,0,
     +                  MPI_COMM_WORLD,IERR)
           
C         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR) !just in case ...
      ELSE
         CALL INPUTO(NBLOCK,IMAX,JMAX,KMAX,INTERI,INTERJ,
     +        INTERK,IDER,LAMIN,INITC,IT,IL,IK,IDI1,IDI2,IDI3,
     +        MOV,MGRID,MIB,MIT,MJB,MJT,MKB,MKT,
     +        IROTVE,MGM,INPARA,INPAR)
      ENDIF ! PARALLEL

      IF(MASTER) DEALLOCATE(INPARA,RNPARA,CNPARA)
C ******** END OF RECEIVED BLOCK FROM PROCESS 0 *************************
C ... Calculate NMOV size (connected to sub INPGRI)
      NMOV = 0
      DO N = 1,NBLOCK
         IF(MOV(N) == 1) THEN
            NMOV = NMOV + 1
         ELSEIF(MOV(N) == 2) THEN
            DO IP = 1,NBCS
               NBL    = ICON(2+(IP-1)*IC9) ! BLOCK NUMBER
               IF(N == NBL) THEN
               IBT = ICON(1+(IP-1)*IC9) ! BC TYPE
               IF((IBT >= 8 .AND. IBT <= 10).OR.IBT == 4) NMOV=NMOV+1
               IF(IBT == 13) NMOV=NMOV+1
               ENDIF
            ENDDO
         ELSEIF(MOV(N) == 3) THEN
            NMOV = NMOV + 1
         ENDIF
      ENDDO ! N = 1,NBLOCK

C ... The minimum length of these arrays should be NBLOCK (TSii 16.1.2004)
C ... In a case of a movie file NMOV (TSii 3.4.2008)

      NMOVA = MAX(NMOV,NBLOCK+1)

      ALLOCATE(IMAXM(NMOVA),JMAXM(NMOVA),KMAXM(NMOVA),
     +  IMINM(NMOVA),JMINM(NMOVA),KMINM(NMOVA),MOVPO(NMOVA))

      IF(NMOV == 0) THEN
         WRITE(45,'(/A,I6/)')'No movies, MOVPO etc. allocated to',NMOVA
      ELSE IF(NMOV > 0) THEN
         WRITE(45,'(/A,I6/)')'Movies specified, MOVPO etc. allocated to'
     +   ,NMOVA
         WRITE(13,'(2A,I6/)')'  I will make movies to MOVIE directory,'
     +   ,' NMOV= ',NMOVA
         WRITE(4,'(/A,A,/A,I6/)')'  MOVIE files are made to MOVIE ',
     +  'directory','  Blockwise information in MEMORY-file, NMOV = ',
     +   NMOVA
         WRITE(45,'(A)') ' Blockwise movie information:'
         WRITE(45,9630)
9630     FORMAT(1X,28('*')/)
9631     FORMAT(2X,13('='))
         WRITE(45,'(A)') '  Block  MOV(N)'
         WRITE(45,9631)
         DO N = 1,NBLOCK
         WRITE(45,'(1X,I4,6X,I1)') N,MOV(N)
         ENDDO
         WRITE(45,9631)
      ELSE
         WRITE(*,*) 'Negative NMOV found. Exiting...'
         STOP
      ENDIF

C ********************************************************************

      CALL NLOCAS(NLOCAL,NPNUM,NBLOCG,NBLOCK,NPROCE,MBPRO,IPRO,NPRO)

C ... REANEW CONTROLS IF ASM IS STARTED (FALSE) OR JUST CONTINUE (TRUE)

      REANEW = .FALSE.
      IF((ITURB  >= 10 .AND.ITURB <= 19) .AND. 
     +   (ITURBO < 10 .OR. ITURB > 19)) THEN
         REANEW = .TRUE.
         IF(MASTER) WRITE(4,*) ' ASM MODEL IS STARTED'
      ENDIF
      IF(IOLD <= 0) THEN
         DO 3046 N = 1,NBLOCG
C            IF(OMEGA(N) /= 0.) ROTANG(N) = ROTAT(N)              
            IF(IGRID(N,1) == 1) ROTANG(N) = ROTAT(N)
 3046    CONTINUE
      ENDIF
       
      IOLD1    = IOLD

 9104 FORMAT('  NUMBER OF GRID LEVELS WAS CHANGED FROM',I3,1X,'TO',I3,
     +1X,'IN BLOCK',I3)

      IF(.NOT. STARTL .AND. MASTER) THEN
         WRITE(*,*)'  RSTART-file will not be written'
      ENDIF
      
C ***************** SET UP IOS ****************************************

      ICYCTI = ICYCLE
      ICYTIT = ICYTOT

      IF(TIMEL) THEN
        ICYCTI = ITIMES
        ICYTIT = ICYTOT + 3*ITIMES
        ICYCLE = 0
      ENDIF

      IF(IOLD1 > 0) THEN

         IF(.NOT.TIMEL .AND. MASTER) THEN
            CALL WTC(9,ICYTOT)
            CALL WTCFGS(155,ICYTOT)
         ENDIF

         CALL LOCATE_LINE(7,ICYTIT)

         IF(MASTER) THEN
            CALL WTC(59,ICYCTI)
            CALL WTCFGS(156,ICYTOT)
         ENDIF

         IF(PARALLEL .AND. MASTER) CALL LOCATE_LINE(11,ICYTIT)

         ICYOLD = ICYTOT

      ENDIF

C ****************** DO SOME PRE CALCULATIONS *************************

      ALPHA   = ALPHA *DEG2RAD
      BETA    = BETA  *DEG2RAD
      AGAMMA  = AGAMMA*DEG2RAD
      ABANK   = ABANK *DEG2RAD

C ... PREDETERMINED VALUE FOR THE IMPLICIT STAGE

      THETA   = 1.5

C *********************************************************************

      IF(ITURB >= 10.AND.ITURB <= 19 .OR. ISTRES > 0) THEN ! Muuta testit näin
         STRESL = .TRUE.
         IF(MASTER) 
     +   WRITE(4,*)' Strain components are needed if ASM is used'
      ENDIF

C ... IDIS = 0  no differential lenght scale
C ... IDIS = 1  epsilon lenght scale
C ... IDIS = 2  omega lenght scale

      IDIS = 0
      IF(ITURB >= 3) IDIS = 1
      IF(ITURB == 6) IDIS = 2

      NSCAL     = KSCAL                    ! NO. OF SCALARS
c      IF(ITURB == 9)  NSCAL  = KSCAL + 1 ! Spalart-Allmaras + scalars
      IF(ITURB >= 21) NSCAL  = KSCAL + 6 ! REYNOLDS STRESSES + SCALARS
      IF(ITURB <= 24) KSCAL  = NSCAL     ! FOR UNCOUPLED REYNOLDS STRESSES
      IF(IFLUX == 5)  KSCAL  = NSCAL - 6 ! FULLY COUPLED
      JSCAL     = NSCAL                    ! NO. OF UPPDATED SCALARS
      IF(ITURB >= 21 .AND. ITURB <= 22) JSCAL  = NSCAL - 6
      MAXNSC = NSCAL + 1

      ALLOCATE(RLOLIM(MAXNSC),UPPLIM(MAXNSC),PSIGSC(MAXNSC),
     +         PSIGS2(MAXNSC))

      IF (NSCAL > 0) THEN
       IF(MASTER) THEN
       DO 3150 NS = 1,NSCAL
          READ(43,*) RLOLIM(NS),UPPLIM(NS),PSIGSC(NS),PSIGS2(NS)
 3150  CONTINUE
       CLOSE(43)
       ENDIF
      
       IF(PARALLEL) THEN
       CALL MPI_BCAST(RLOLIM(1),NSCAL,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(UPPLIM(1),NSCAL,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(PSIGSC(1),NSCAL,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(PSIGS2(1),NSCAL,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
       ENDIF

      ENDIF ! NSCAL > 0

      IF (ITURB >= 21) THEN
         IF(JRDIF == 0 .OR.JRDIS == 0 .OR.JRPRE == 0 )THEN
            WRITE(*,*) 'RSM is not fully defined.'
            WRITE(*,*) 'Aborting....'
            STOP
         ENDIF
         IF (JRIMP >= 2 .OR. JRIMP < 0) THEN
            WRITE(*,*) 'RSM time step is not defined.'
            WRITE(*,*) 'Aborting....'
            STOP
         ENDIF
      ENDIF ! ITURB >= 21

      IF (ITURB <= 23) JRIMP = 0
      IF (ITURB < 20) THEN
         IF(JRDIF /= 0 .OR.JRDIS /= 0 .OR.JRPRE /= 0.OR.JRIMP /= 0)THEN
            WRITE(*,*) 'You cannot have any diffusion of Reynolds'//
     +           ' stesses'
            WRITE(*,*) 'Putting zeros....'
            JRDIF = 0
            JRDIS = 0
            JRPRE = 0
         ENDIF
      ENDIF ! ITURB < 20

      IF(ITURB >= 21 .AND. IROTCO >= 1) THEN
         WRITE(13,*) 
         WRITE(13,*)'Warning IROTCO is',IROTCO,' and RSM is used'
         WRITE(13,*)'In general, RSM do not need rotational corrections'
         WRITE(13,*) 
      ENDIF
      IF(ITURB >= 10 .AND. FULLNS) THEN
         WRITE(*,*)  'Warning, no exact diffusion for Reynolds stresses'
         WRITE(13,*) 'Warning, no exact diffusion for Reynolds stresses'
      ENDIF

C *********************************************************************
C ... Check convergence history treatment is leagal

      IF(IDRXX == 1 .AND. ISTATE == 10) THEN
         WRITE(*,*) 'With an incompressible simulation you cannot'//
     +        ' ask maximum change of density'
         WRITE(*,*) 'in convergence history.'
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

      IF(IDRXX >= 7 .AND. IDRXX <= 8 .AND. ITURB <= 2) THEN
         WRITE(*,*) 'With algebraic turbulence model you cannot'//
     +        ' ask maximum change of RK or REPS'
         WRITE(*,*) 'in convergence history.'
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF
      IF(IDRXX-8 > NSCAL) THEN
         WRITE(*,*) 'You have only',NSCAL,' scalars and you asked'
         WRITE(*,*) IDRXX-8,'th scalar in convergence history.'
         WRITE(*,*) 'Exiting ...'
         STOP
      ENDIF

C **********************************************************************

      IF(MASTER) THEN
         OPEN(UNIT=1,FILE='RUN',STATUS='UNKNOWN')
         WRITE(1,'(A3)') 'run'
         WRITE(1,'(I7,3X,A5)') ICMAX,'ICMAX'
         CLOSE(1)
      ENDIF

C ***********************************************************************
C
C ... SPECIFY REFERENCE AND FREE-STREAM QUANTITIES
C ... istate=6 (=vesi) lisatty 12.1.96 KPe
C

C ... AND FOR THE EQUATION OF STATE (PERFECT GAS)

      IF(STATEL .AND. ISTATE <= 5) THEN
        IF(MASTER .AND. ISTATE == 0) WRITE(4,9126) ! Unused?
        RGAS  = 287.05287
        GAMMA = 1.4
        VISU0 = 1.458E-6
        EXPSU = 1.5
        TSU0  = 110.4
        T0REF = 0.
        E0REF = 0.
        IF(ISTATE == 0) ISTATE = 1
        IF(MASTER) WRITE(4,9127)
        WRITE(4,*)'****************************************************'
        WRITE(4,*) ' GAS CONSTANT =',REAL(RGAS,4), ' K'
        WRITE(4,*) ' GAMMA        =',REAL(GAMMA,4)
        WRITE(4,*) ' MU0          =',REAL(VISU0,4),' K'
        WRITE(4,*) ' EXP          =',REAL(EXPSU,4)
        WRITE(4,*) ' T0           =',REAL(TSU0,4), ' K'
        WRITE(4,*) ' E0REF        =',REAL(E0REF,4), ' (m/s)**2'
        WRITE(4,*) ' T0REF        =',REAL(T0REF,4), ' K'
      ELSE IF(ISTATE <= 5) THEN
        IF(MASTER) THEN
        WRITE(4,9128)
        WRITE(4,*)'****************************************************'
        WRITE(4,*) ' GAS CONSTANT =',REAL(RGAS,4), ' K'
        WRITE(4,*) ' GAMMA        =',REAL(GAMMA,4)
        WRITE(4,*) ' MU0          =',REAL(VISU0,4),' K'
        WRITE(4,*) ' EXP          =',REAL(EXPSU,4)
        WRITE(4,*) ' T0           =',REAL(TSU0,4), ' K'
        WRITE(4,*) ' E0REF        =',REAL(E0REF,4), ' (m/s)**2'
        WRITE(4,*) ' T0REF        =',REAL(T0REF,4), ' K'
      ENDIF
      ENDIF ! STATEL
      
9126  FORMAT('  NO SELECTION FOR THE FREE-STREAM')
9127  FORMAT('  STANDARD PERFECT GAS MODEL IS ASSUMED FOR FREE-',
     + 'STREAM')
9128  FORMAT('  NON-STANDARD PERFECT GAS MODEL IS APPLIED FOR FREE-',
     + 'STREAM')
     
      WRITE(4,'(/,2X,A)') 'Default initialization of the free-stream:'
      WRITE(4,'(2X,A)')   '========================================='
      WRITE(4,'(2X,A,I2)') 'Assuming ISTATE =',ISTATE

      IF (ISTATE == 6) THEN
        CALL BCINI2(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,SIEINI,
     +  FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,VISU0,EXPSU,TSU0,
     +  E0REF,T0REF,RMULTV,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI2 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI2'
      ELSEIF (ISTATE == 7) THEN
        CALL BCINI3(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,SIEINI,
     +  FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,VISU0,EXPSU,TSU0,
     +  E0REF,T0REF,RMULTV,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI3 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI3'
      ELSEIF (ISTATE == 8) THEN
        FRAPRE(1) = FRSPRE
        FRATEM(1) = FRSTEM
        CALL WATER2(FRAVIS,FRAPRE,FRADEN,FRASIE,FRATEM,FRACP,FRACH,
     +  FRDRDP,FRDRDH,FRPR,1,1,1,1,0,0,0,0)
        FRSVIS = FRAVIS(1)
        FRSDEN = FRADEN(1)
        FRSSIE = FRASIE(1)
        FRSCP  = FRACP(1)
        FRSCH  = FRACH(1)
        FRSSSP = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRSDEN)
        CALL BCINI4(FRSVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,SIEINI,FRSSSP,
     +  FRSVIS,T0,RE,RMACH,CHLREF,FRSCP,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI4 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI4 (ISTATE=8'
      ELSEIF (ISTATE == 9) THEN 
        FRAPRE(1) = FRSPRE
        FRATEM(1) = FRSTEM
        CALL STEAM3(FRAVIS,FRAPRE,FRADEN,FRASIE,FRATEM,FRACP,FRACH,
     +  FRDRDP,FRDRDH,FRPR,1,1,1,1,0,0,0,0)
        FRSVIS = FRAVIS(1)
        FRSDEN = FRADEN(1)
        FRSSIE = FRASIE(1)
        FRSCP  = FRACP(1)
        FRSCH  = FRACH(1)
        FRSSSP = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRSDEN)
        CALL BCINI4(FRSVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,SIEINI,FRSSSP,
     +  FRSVIS,T0,RE,RMACH,CHLREF,FRSCP,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI4 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI4 (ISTATE=9)'
      ELSEIF (ISTATE == 10) THEN 
        CALL BCINI5(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,SIEINI,
     +  FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,VISU0,EXPSU,TSU0,
     +  E0REF,T0REF,RMULTV,PSEUCO,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI5 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI5 (ISTATE=10'
      ELSEIF (ISTATE == 11) THEN 
        CALL BCINI6(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,SIEINI,
     +  FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,VISU0,EXPSU,TSU0,
     +  E0REF,T0REF,RMULTV,PSEUCO,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI6 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI6 (ISTATE=11)'
      ELSEIF (ISTATE == 12) THEN
        FRAPRE(1) = FRSPRE
        FRATEM(1) = FRSTEM
        CALL R12LIQ(FRAVIS,FRAPRE,FRADEN,FRASIE,FRATEM,FRACP,FRACH,
     +  FRDRDP,FRDRDH,FRPR,1,1,1,1,0,0,0,0)
        FRSVIS = FRAVIS(1)
        FRSDEN = FRADEN(1)
        FRSSIE = FRASIE(1)
        FRSCP  = FRACP(1)
        FRSCH  = FRACH(1)
        FRSSSP = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRSDEN)
        CALL BCINI4(FRSVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,SIEINI,FRSSSP,
     +  FRSVIS,T0,RE,RMACH,CHLREF,FRSCP,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI4 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI4 (ISTATE=12)'
      ELSEIF (ISTATE == 13) THEN 
        FRAPRE(1) = FRSPRE
        FRATEM(1) = FRSTEM
        CALL R12GAS(FRAVIS,FRAPRE,FRADEN,FRASIE,FRATEM,FRACP,FRACH,
     +  FRDRDP,FRDRDH,FRPR,1,1,1,1,0,0,0,0)
        FRSVIS = FRAVIS(1)
        FRSDEN = FRADEN(1)
        FRSSIE = FRASIE(1)
        FRSCP  = FRACP(1)
        FRSCH  = FRACH(1)
        FRSSSP = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRSDEN)
        CALL BCINI4(FRSVEL,FRSPRE,FRSDEN,FRSTEM,FRSSIE,SIEINI,FRSSSP,
     +  FRSVIS,T0,RE,RMACH,CHLREF,FRSCP,REFVEL)
      WRITE(13,'(A)')'  Initialization was made in BCINI4 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI4 (ISTATE=13)'
      ELSE
        CALL BCINI1(FRSVEL,FRSPRE,FRSDEN,FRSTEM,TEMINI,FRSSIE,SIEINI,
     +  FRSSSP,FRSVIS,T0,RE,RMACH,CHLREF,RGAS,GAMMA,VISU0,EXPSU,TSU0,
     +  E0REF,T0REF,RMULTV,REFVEL,GROUND,ALTITUDE)
      WRITE(13,'(A)')'  Initialization was made in BCINI1 (see FORCES)'
      WRITE(4,'(A)')'  Initialization was made in BCINI1'

      ENDIF ! ISTATE == 6
      
      IF(ABS(REFVEL) < EPS) REFVEL = FRSVEL

      FRSMUD = 1.0E-3           ! Default value for free-stream mut/mu
      FROUDE = FRSVEL/SQRT(G0*CHLREF)
      IF (REFVEL /= 0) FROUDE = REFVEL/SQRT(G0*CHLREF)
      IF(OMEGAREF /= 0.) THEN
         ADVANJ = 1./ABS(OMEGAREF/(2.*PII)*CHLREF/FRSVEL) ! (Mersu)
         IF (REFVEL /= 0) ADVANJ = 1./ABS(OMEGAREF/(2.*PII)*
     +        CHLREF/REFVEL)    ! (Mersu)
      ELSE
         ADVANJ = 0.
      ENDIF
       
      IF(ITURB >= 3 .AND. ITURB /= 9) THEN
         RKII   = 1.5*FRSDEN*(RKLIM*REFVEL)**2
         EPSLIM = MAX(EPSLIM,EPS6)
         FRSMUT = EPSLIM*FRSVIS                  ! IS READ IN AS EPSLIM
         REPSI  = 0.09*RKII**2/FRSMUT  
         RMUINI = RMUINI*FRSVIS
         IF(RKLIM < 1.E-15) THEN
            FRSMUT = FRSMUD*FRSVIS
            RMULIM = FRSMUT 
            RKII   = FRSMUT*1.0*REFVEL/CHLREF
            REPSI  = 0.09*RKII**2/RMULIM
         ENDIF
C      IF(IDIS == 2) REPSI = .001*FRSDEN*RKII/FRSMUT
         IF(IDIS == 2) THEN
            REPSI  = FRSDEN*RKII/FRSMUT  !        ! AHe 28.2.2001
         ENDIF ! IDIS  == 2
         FRSLEN    = (FRSMUT/FRSVIS)*(FRSVIS/(SQRT(FRSDEN)*SQRT(RKII)))
      ELSE IF(ITURB == 9) THEN ! Spalart-Allmaras
         RKII   = 1.E-10
         RMUINI = RMUINI*FRSVIS
         FRSMUT = EPSLIM*FRSVIS
         RMULIM = 1.E-10
         RKHI   = 1.
         DO III = 1,20
            RKHI= SQRT(SQRT(FRSMUT/FRSVIS*(RKHI**3+7.1**3)))
         ENDDO
         REPSI  = 1.E-12  ! RKHI*FRSVIS ! Free stream is a minimum? 
      ELSE ! An algebraic model
         RKII   = 1.E-10
         REPSI  = 1.E-10
         RMUINI = 1.E-10
         FRSMUT = 1.E-10
         RMULIM = 1.E-10
      ENDIF ! ITURB >= 3
C
C ... Reference pressures and limitations

c      IF (REFPRE == 0. .AND. OMEGAREF /=  0.) THEN ! Rotating flow case
c         REFPRE = .5*FRSDEN*(OMEGAREF/(2.*PII)*CHLREF)**2
c         RE     = FRSDEN*OMEGAREF/(2.*PII)*CHLREF/FRSVIS
c      ELSE IF (REFPRE == 0.) THEN

c      IF (REFVEL == 0.) THEN
c         REFVEL = FRSVEL
c      ENDIF

         REFPRE = .5*FRSDEN*REFVEL**2

c      ENDIF

      IF (DIFPRE < 0.) DIFPRE = FRSPRE

      RKLIM     = RKII
      EPSLIM    = REPSI
      TLOLIM    = 0.05*FRSTEM
      TUPLIM    = 6.*FRSTEM

      IF(ARTSSP <= 0.) THEN
         ARTSSP = PSEUCO*REFVEL
      ELSE IF(REFVEL > 0.) THEN
         PSEUCO = ARTSSP/REFVEL
      ELSE
         WRITE(*,*) ' REFVEL not defined for PSEUCO, Exiting..'
         WRITE(13,'(2X,2A,/,A)')'REFVEL is not defined, but ARTSSP is',
     +   ' given for PSEUCO.',' Check the code, now exiting..'
      ENDIF

      IF(FRESUL) FREDIF = FREDIF*0.01*(FRSVEL+REFVEL)*CHLREF/SQRT(RE)

C ... Blockwise values for the free-stream variables and limitations

C **********************************************************************

      WRITE(4,'(/,2X,A)') 'Applied blockwise initialization for the free
     +-stream:'
      WRITE(4,'(2X,A)')   '=============================================
     +======='
  
      DO N = 1,NBLOCK

      NGL    = NPROCE(N+1,IPRO)

      DO IPHASE = 1,NPHASE ! NPHASES !BLKS(NGL)%NPHASE

        FRAPRE(IPHASE) = BLKS(NGL)%FRSPRE
        FRATEM(IPHASE) = BLKS(NGL)%FRSTEM
        FRAMA(IPHASE)  = BLKS(NGL)%RMACH
        FRARE(IPHASE)  = BLKS(NGL)%RE
        FRADEN(IPHASE) = BLKS(NGL)%FRADEN(IPHASE) !BLKS(NGL)%FRSDEN
cc        IF(JSTATE(NGL,IPHASE) == 3 .AND. FRAPRE(IPHASE) /= 0. .AND.
cc     +  FRATEM(IPHASE) /= 0.) FRADEN(IPHASE) = 0.

        IF(NPHASE > 1) THEN
        IF (ABS(BLKS(NGL)%FRADEN(1)*BLKS(NGL)%FRSALFA(1)+
     +       BLKS(NGL)%FRADEN(NPHASE)*BLKS(NGL)%FRSALFA(NPHASE)-
     +       BLKS(NGL)%FRSDEN)/BLKS(NGL)%FRSDEN >= 0.05) THEN
           WRITE(*,*)'!!!!!!!!!!!! WARNING WARNING WARNING !!!!!!!!!!!!'
           WRITE(*,*)'!!!!!!!!!!!! WARNING WARNING WARNING !!!!!!!!!!!!'
           WRITE(*,*)'Check FRSDENB and FRADENB values in block',NGL
           WRITE(*,*)' '
        ENDIF
        ENDIF
        
      IF(IPHASE == 1 .AND. BLKS(NGL)%NPHASE == 1) THEN
c         WRITE(4,'(2X,A,I3)') 'In (global) block',NGL
         WRITE(4,'(2X,A,I6)') 'In (global) block',NGL
      ELSE IF(IPHASE == 1 .AND. BLKS(NGL)%NPHASE > 1) THEN
c         WRITE(4,'(2X,A,I3,A,I3,A,I2)') 
         WRITE(4,'(2X,A,I6,A,I3,A,I2)') 
     +  'In (global) block',NGL,', phase',IPHASE,', JSTATE =',
     +   JSTATE(NGL,IPHASE)
      ELSE IF(BLKS(NGL)%NPHASE > 1) THEN
c         WRITE(4,'(2X,A,I3,A,I3,A,I2)') 
         WRITE(4,'(2X,A,I6,A,I3,A,I2)') 
     +  'In (global) block',NGL,', phase',IPHASE,', JSTATE =',
     +   JSTATE(NGL,IPHASE)
      ENDIF

      IF (JSTATE(NGL,IPHASE) == 6) THEN ! Water (old version)
        CALL BCINI2(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,FRASIE(IPHASE),
     +  FRASIE(IPHASE),FRASSP(IPHASE),FRAVIS(IPHASE),
     +  BLKS(NGL)%T0,FRARE(IPHASE),FRAMA(IPHASE),BLKS(NGL)%CHLREF,
     +  RGAS,GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,BLKS(NGL)%RMULTV,REFVEL)
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)
      ELSEIF (JSTATE(NGL,IPHASE) == 7) THEN ! Water, modified viscosity
        CALL BCINI3(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,FRASIE(IPHASE),
     +  FRASIE(IPHASE),FRASSP(IPHASE),FRAVIS(IPHASE),
     +  BLKS(NGL)%T0,FRARE(IPHASE),FRAMA(IPHASE),BLKS(NGL)%CHLREF,
     +  RGAS,GAMMA,VISU0,EXPSU,TSU0,E0REF,T0REF,BLKS(NGL)%RMULTV,REFVEL)
      ELSEIF (JSTATE(NGL,IPHASE) == 8) THEN ! Water as f(p,T)
        CALL WATER2(FRAVIS(IPHASE),FRAPRE,FRADEN(IPHASE),FRASIE(IPHASE),
     +  FRATEM,FRACP(IPHASE),FRACH(IPHASE),FRDRDP,FRDRDH,FRPR,1,1,1,1,
     +  0,0,0,0)
        FRASSP(IPHASE) = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRADEN(IPHASE))
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)
        CALL BCINI4(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,FRASIE(IPHASE),BLKS(NGL)%SIEINI,
     +  FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,FRARE(IPHASE),
     +  FRAMA(IPHASE),BLKS(NGL)%CHLREF,FRACP(IPHASE),REFVEL)
      ELSEIF (JSTATE(NGL,IPHASE) == 9) THEN ! Steam as f(p,T)
        CALL STEAM3(FRAVIS(IPHASE),FRAPRE,FRADEN(IPHASE),FRASIE(IPHASE),
     +  FRATEM,FRACP(IPHASE),FRACH(IPHASE),FRDRDP,FRDRDH,FRPR,1,1,1,1,
     +  0,0,0,0)
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)
        FRASSP(IPHASE) = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRADEN(IPHASE))
        CALL BCINI4(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,FRASIE(IPHASE),BLKS(NGL)%SIEINI,
     +  FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,FRARE(IPHASE),
     +  FRAMA(IPHASE),BLKS(NGL)%CHLREF,FRACP(IPHASE),REFVEL)
      ELSEIF(JSTATE(NGL,IPHASE) == 1 .OR. JSTATE(NGL,IPHASE) == 3) THEN ! Perf. gas
        CALL BCINI1(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,FRASIE(IPHASE),
     +  FRASIE(IPHASE),FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,
     +  FRARE(IPHASE),FRAMA(IPHASE),BLKS(NGL)%CHLREF,RGAS,GAMMA,VISU0,
     +  EXPSU,TSU0,E0REF,T0REF,RMULTV,REFVEL,GROUND,ALTITUDE)
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)
      ELSEIF(JSTATE(NGL,IPHASE) == 10) THEN ! Totally incompressible (p unlimited)
        IF(BLKS(NGL)%IPRESC /= 1) THEN ! Only pressure correction (10)
c        WRITE(4,'(2X,A,I3,2A,/,2A)') 'In (global) block',NGL,' a fully',
        WRITE(4,'(2X,A,I6,2A,/,2A)') 'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution','  type must be',
     +  ' IPRESC=1 (pressure correction method)'
c        WRITE(*,'(2X,A,I3,2A,/,A)') 'In (global) block',NGL,' a fully',
        WRITE(*,'(2X,A,I6,2A,/,A)') 'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution type',
     +  '  must be IPRESC=1. Exiting..'
        WRITE(13,'(2X,A)')'must be IPRESC=1. Exiting..'
c        WRITE(13,'(2X,A,I3,2A,/,A)')'In (global) block',NGL,' a fully',
        WRITE(13,'(2X,A,I6,2A,/,A)')'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution',
     +  '  type must be IPRESC=1. Exiting..'
        WRITE(13,'(2X,A)') 'Sorry about that chief.'
        STOP
        ENDIF
        CALL BCINI5(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,FRASIE(IPHASE),
     +  FRASIE(IPHASE),FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,
     +  FRARE(IPHASE),FRAMA(IPHASE),BLKS(NGL)%CHLREF,RGAS,GAMMA,VISU0,
     +  EXPSU,TSU0,E0REF,T0REF,RMULTV,PSEUCO,REFVEL)
      ELSEIF(JSTATE(NGL,IPHASE) == 11) THEN ! Sea water at constant density
        IF(BLKS(NGL)%IPRESC /= 1) THEN ! Only pressure correction (10)
c        WRITE(4,'(2X,A,I3,2A,/,2A)') 'In (global) block',NGL,' a fully',
        WRITE(4,'(2X,A,I6,2A,/,2A)') 'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution','  type must be',
     +  ' IPRESC=1 (pressure correction method)'
c        WRITE(*,'(2X,A,I3,2A,/,A)') 'In (global) block',NGL,' a fully',
        WRITE(*,'(2X,A,I6,2A,/,A)') 'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution type',
     +  '  must be IPRESC=1. Exiting..'
        WRITE(13,'(2X,A)')'must be IPRESC=1. Exiting..'
c        WRITE(13,'(2X,A,I3,2A,/,A)')'In (global) block',NGL,' a fully',
        WRITE(13,'(2X,A,I6,2A,/,A)')'In (global) block',NGL,' a fully',
     +  ' incompressible flow was assumed, solution',
     +  '  type must be IPRESC=1. Exiting..'
        WRITE(13,'(2X,A)') 'Sorry about that chief.'
        STOP
        ENDIF
        CALL BCINI6(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,BLKS(NGL)%TEMINI,FRASIE(IPHASE),
     +  FRASIE(IPHASE),FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,
     +  FRARE(IPHASE),FRAMA(IPHASE),BLKS(NGL)%CHLREF,RGAS,GAMMA,VISU0,
     +  EXPSU,TSU0,E0REF,T0REF,RMULTV,PSEUCO,REFVEL)
      ELSEIF (JSTATE(NGL,IPHASE) == 12) THEN ! Liquid freon R12
        CALL R12LIQ(FRAVIS(IPHASE),FRAPRE,FRADEN(IPHASE),FRASIE(IPHASE),
     +  FRATEM,FRACP(IPHASE),FRACH(IPHASE),FRDRDP,FRDRDH,FRPR,1,1,1,1,
     +  0,0,0,0)
        FRASSP(IPHASE) = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRADEN(IPHASE))
        CALL BCINI4(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,FRASIE(IPHASE),BLKS(NGL)%SIEINI,
     +  FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,FRARE(IPHASE),
     +  FRAMA(IPHASE),BLKS(NGL)%CHLREF,FRACP(IPHASE),REFVEL)
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)
      ELSEIF (JSTATE(NGL,IPHASE) == 13) THEN ! Gaseous freon R12
        CALL R12GAS(FRAVIS(IPHASE),FRAPRE,FRADEN(IPHASE),FRASIE(IPHASE),
     +  FRATEM,FRACP(IPHASE),FRACH(IPHASE),FRDRDP,FRDRDH,FRPR,1,1,1,1,
     +  0,0,0,0)
        FRASSP(IPHASE) = 1./SQRT(FRDRDP(1)+FRDRDH(1)/FRADEN(IPHASE))
        CALL BCINI4(BLKS(NGL)%FRSVEL,BLKS(NGL)%FRSPRE,FRADEN(IPHASE),
     +  BLKS(NGL)%FRSTEM,FRASIE(IPHASE),BLKS(NGL)%SIEINI,
     +  FRASSP(IPHASE),FRAVIS(IPHASE),BLKS(NGL)%T0,FRARE(IPHASE),
     +  FRAMA(IPHASE),BLKS(NGL)%CHLREF,FRACP(IPHASE),REFVEL)
        BLKS(NGL)%FRADEN(IPHASE) = FRADEN(IPHASE)

      ELSEIF(JSTATE(NGL,IPHASE) == 101) THEN ! Cast iron
        CALL CAST_IRON(FRADEN(IPHASE),FRASIE(IPHASE),FRATEM,
     +  FRACP(IPHASE),FRACH(IPHASE),1,1,1,1,0,0,0,0)
        BLKS(NGL)%T0     = 0.
        FRASSP(IPHASE)   = 7500.
        WRITE(4,*)' SOLID BLOCK PROPERTIES BASED ON A GIVEN TEMPERATURE'
      ENDIF ! ISTATE == 6

      BLKS(NGL)%FRSX(IPHASE)= BLKS(NGL)%FRSALFA(IPHASE)*FRADEN(IPHASE)

      ENDDO ! IPHASE

      BLKS(NGL)%FRSVIS   = 0.
      BLKS(NGL)%FRSDEN   = 0. 
      BLKS(NGL)%FRSSIE   = 0.
      BLKS(NGL)%FRSCP    = 0.
      BLKS(NGL)%FRSCH    = 0.
      BLKS(NGL)%FRSSSP   = 0.
      BLKS(NGL)%SIEINI   = 0.

      IF(.NOT. NAMEL) THEN
         BLKS(NGL)%TLOLIM = TLOLIM
         BLKS(NGL)%TUPLIM = TUPLIM
      ENDIF

      DO IPHASE = 1,BLKS(NGL)%NPHASE
          
        BLKS(NGL)%FRSVIS = BLKS(NGL)%FRSVIS + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRAVIS(IPHASE)
        BLKS(NGL)%FRSDEN = BLKS(NGL)%FRSDEN + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRADEN(IPHASE)
        BLKS(NGL)%FRSSIE = BLKS(NGL)%FRSSIE + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRASIE(IPHASE)
        BLKS(NGL)%SIEINI = BLKS(NGL)%SIEINI + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRASIE(IPHASE)
        BLKS(NGL)%FRSCP  = BLKS(NGL)%FRSCP  + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRACP(IPHASE)
        BLKS(NGL)%FRSCH  = BLKS(NGL)%FRSCH  + 
     +                     BLKS(NGL)%FRSALFA(IPHASE)*FRACH(IPHASE)
        BLKS(NGL)%FRSSSP = BLKS(NGL)%FRSSSP + BLKS(NGL)%FRSALFA(IPHASE)/
     +                    (FRADEN(IPHASE)*FRASSP(IPHASE)**2)
      ENDDO
            
      DO IPHASE = 1,BLKS(NGL)%NPHASE
      BLKS(NGL)%FRSX(IPHASE) = BLKS(NGL)%FRSX(IPHASE)/BLKS(NGL)%FRSDEN 
      ENDDO

      IF(ABS(BLKS(NGL)%REFVEL) < EPS) THEN ! Added 3.12.2013
         BLKS(NGL)%REFVEL = REFVEL
      ENDIF
      IF(ABS(BLKS(NGL)%REFVEL) < EPS) THEN
         BLKS(NGL)%REFVEL = BLKS(NGL)%FRSVEL
      ENDIF

      IF(ABS(BLKS(NGL)%REFCAV) < EPS) THEN
         BLKS(NGL)%REFCAV = SQRT(BLKS(NGL)%FRSVEL**2 +
     +                      (0.7*.5*BLKS(NGL)%CHLREF*OMEGA(NGL))**2)
      ENDIF

      IF(ABS(BLKS(NGL)%CHLCAV) < EPS) THEN
         BLKS(NGL)%CHLCAV = BLKS(NGL)%CHLREF
      ENDIF


      BLKS(NGL)%FRSSSP   = BLKS(NGL)%FRSDEN * BLKS(NGL)%FRSSSP
      BLKS(NGL)%FRSSSP   = 1./SQRT(BLKS(NGL)%FRSSSP)
      BLKS(NGL)%RMACH    = BLKS(NGL)%REFVEL/BLKS(NGL)%FRSSSP
      BLKS(NGL)%RE       = BLKS(NGL)%FRSDEN*BLKS(NGL)%REFVEL*
     +                     BLKS(NGL)%CHLREF/BLKS(NGL)%FRSVIS

      FRSMUD = 1.0E-3           ! Default value for free-stream mut/mu

      IF(JSTATE(NGL,1) <= 99) THEN ! Fluid block

      IF(ITURB >= 3 .AND. ITURB /= 9) THEN
         IF(BLKS(NGL)%FRSTUR >= 1.E-15) THEN
         RKIIB   = 1.5*BLKS(NGL)%FRSDEN*(BLKS(NGL)%FRSTUR * 
     +              BLKS(NGL)%REFVEL)**2
         FRSMUTB = BLKS(NGL)%FRSMUT * BLKS(NGL)%FRSVIS
         REPSIB  = 0.09*RKIIB**2/FRSMUTB  
         RMUINIB = BLKS(NGL)%RMUINI * BLKS(NGL)%FRSVIS
         ELSEIF(BLKS(NGL)%FRSTUR < 1.E-15) THEN
            FRSMUTB = FRSMUD*BLKS(NGL)%FRSVIS
            RMULIMB = FRSMUTB 
            RKIIB   = FRSMUTB*1.0*BLKS(NGL)%REFVEL/BLKS(NGL)%CHLREF
            REPSIB  = 0.09*RKIIB**2/RMULIMB
            BLKS(NGL)%FRSTUR = SQRT(2./3.*RKIIB/BLKS(NGL)%FRSDEN)/REFVEL
         ENDIF
         IF(IDIS == 2) THEN
            REPSIB  = BLKS(NGL)%FRSDEN*RKIIB/FRSMUTB  
         ENDIF ! IDIS  == 2
      ELSE IF(ITURB == 9) THEN ! Spalart-Allmaras
         RKIIB   = 1.E-10
         RMUINIB = BLKS(NGL)%RMUINI * BLKS(NGL)%FRSVIS
         FRSMUTB = BLKS(NGL)%FRSMUT * BLKS(NGL)%FRSVIS
         RKHI    = 1.
         DO III  = 1,20
            RKHI = SQRT(SQRT(FRSMUT/FRSVIS*(RKHI**3+7.1**3)))
         ENDDO
         REPSIB  = 1.E-12  !RKHI*BLKS(NGL)%FRSVIS ! Free stream is a minimum? 

C     ELSE ! An algebraic model previously in global block
      ENDIF ! ITURB >= 3
       
      BLKS(NGL)%FRSLEN = (FRSMUTB/BLKS(NGL)%FRSVIS) *
     +  (BLKS(NGL)%FRSVIS / (SQRT(BLKS(NGL)%FRSDEN) * SQRT(RKIIB)))
      BLKS(NGL)%RMUINI = RMUINIB
      BLKS(NGL)%FRSMUT = FRSMUTB
      BLKS(NGL)%RKLIM  = RKIIB
      BLKS(NGL)%EPSLIM = REPSIB
      BLKS(NGL)%TLOLIM = MAX(0.05*BLKS(NGL)%FRSTEM,BLKS(NGL)%TLOLIM)
      
      ELSE ! Solid block

      BLKS(NGL)%FRSLEN = 0.
      BLKS(NGL)%RMUINI = 0.
      BLKS(NGL)%FRSMUT = 0.
      BLKS(NGL)%RKLIM  = 0.
      BLKS(NGL)%EPSLIM = 0.
      BLKS(NGL)%FRSVIS = 0.
      BLKS(NGL)%IPRESC = 3     ! For block report
      ENDIF
      
c      IF(JSTATE(NGL,1) == 6 .OR. JSTATE(NGL,2) == 6 .OR.
c     +   JSTATE(NGL,3) == 6)     BLKS(NGL)%TLOLIM = 274. ! Does not freeze
c      IF(JSTATE(NGL,1) == 7 .OR. JSTATE(NGL,2) == 7 .OR.
c     +   JSTATE(NGL,3) == 7)     BLKS(NGL)%TLOLIM = 274. ! Does not freeze
c      IF(JSTATE(NGL,1) == 8 .OR. JSTATE(NGL,2) == 8 .OR.
c     +   JSTATE(NGL,3) == 8)     BLKS(NGL)%TLOLIM = 274. 
      BLKS(NGL)%TUPLIM = MIN(16.*BLKS(NGL)%FRSTEM,BLKS(NGL)%TUPLIM) ! Used to be 6!

C ... Reference pressures and limitations

c      IF(BLKS(NGL)%REFPRE == 0. .AND. OMEGA(NGL) /= 0.)  THEN
c         BLKS(NGL)%REFPRE = .5*BLKS(NGL)%FRSDEN*(OMEGA(NGL)/(2.*PII)*
c     +                       CHLREF)**2
c         BLKS(NGL)%RE     = BLKS(NGL)%FRSDEN*OMEGA(NGL)/(2.*PII)*
c     +                      BLKS(NGL)%CHLREF/BLKS(NGL)%FRSVIS
c      ELSE IF(BLKS(NGL)%REFPRE == 0.) THEN

      IF (BLKS(NGL)%REFVEL == 0.) THEN
         BLKS(NGL)%REFVEL = REFVEL  ! BLKS(NGL)%FRSVEL Global REFVEL?
      ENDIF
         BLKS(NGL)%REFPRE = .5*BLKS(NGL)%FRSDEN * BLKS(NGL)%REFVEL**2

c      ENDIF
      IF(BLKS(NGL)%DIFPRE < 0.) BLKS(NGL)%DIFPRE = BLKS(NGL)%FRSPRE

      PSEUCOB = 0.

      IF(BLKS(NGL)%ARTSSP <= 0.) THEN
         BLKS(NGL)%ARTSSP = BLKS(NGL)%REFVEL * PSEUCO
      ELSE IF(BLKS(NGL)%REFVEL > 0.) THEN
         PSEUCOB = MAX(PSEUCOB,BLKS(NGL)%ARTSSP/BLKS(NGL)%REFVEL)
      ELSE
         WRITE(*,*)' Blockwise REFVEL not defined for PSEUCO, Exiting..'
         WRITE(13,'(2X,2A,I3,/,A)')'REFVEL is not defined, but ARTSSP',
     +   ' is given for PSEUCO for block ',N,' Exiting..'
      ENDIF

      PSEUCO = MAX(PSEUCO,PSEUCOB)

      ENDDO ! N = 1,NBLOCK Mersu

      IF(NBLOCK >= 5)WRITE(4,'(A)')'  (Sorry about so many lines chief)'

      
C ... Initialize RSM limits, only for global state assumption

      IF(ITURB >= 20) THEN
         RLOLIM(1) = 2./3.*RKII
         RLOLIM(4) = 2./3.*RKII
         RLOLIM(6) = 2./3.*RKII
         IF(MASTER) THEN
            IWR77  = 0
            RUUPRE = .5*FRSDEN*REFVEL**2
            IF(UPPLIM(1) < RUUPRE) IWR77 = 1
            IF(UPPLIM(4) < RUUPRE) IWR77 = 1
            IF(UPPLIM(6) < RUUPRE) IWR77 = 1
            RUUPRE = .25*FRSDEN*REFVEL**2
            IF(UPPLIM(2) < RUUPRE) IWR77 = 1
            IF(UPPLIM(3) < RUUPRE) IWR77 = 1
            IF(UPPLIM(5) < RUUPRE) IWR77 = 1
            IF(RLOLIM(2) > -RUUPRE)  IWR77 = 1
            IF(RLOLIM(3) > -RUUPRE)  IWR77 = 1
            IF(RLOLIM(5) > -RUUPRE)  IWR77 = 1
            IF(IWR77 == 1) THEN
               WRITE(13,100)
 100           FORMAT(//'**',10X,'WARNING'/,'**',
     +            'Limiters for Reynolds stresses are not very good.'/
     +            'Values are: ')
               DO NS = 1,6
                  WRITE(13,110) NS,RLOLIM(NS),UPPLIM(NS)
               ENDDO
               WRITE(13,*)
               WRITE(13,*) '** And they should be at least:'
               WRITE(13,110) 1,RLOLIM(1),.5*FRSDEN*REFVEL**2
               WRITE(13,110) 2,-RUUPRE,RUUPRE
               WRITE(13,110) 3,-RUUPRE,RUUPRE
               WRITE(13,110) 4,RLOLIM(4),.5*FRSDEN*REFVEL**2
               WRITE(13,110) 5,-RUUPRE,RUUPRE
               WRITE(13,110) 6,RLOLIM(6),.5*FRSDEN*REFVEL**2
               WRITE(13,*)
 110           FORMAT(3X,I3,2F13.4)
            ENDIF
         ENDIF
      ENDIF   ! ITURB >= 20

C ... End of state initialization

C***********************************************************************
       
      IF(TIMEL .AND. IOLD > 0 .AND. MASTER) THEN
c         WRITE(9) 8,1           ! Version 8.1 (Major=8, Minor=1)
c         WRITE(9) NAME(1:79)
c         WRITE(9) RMACH,ALPHA,RE,0.0,ITURB,0,NRESI,IDRXX
         CALL WTC(9,ICYTOT)
         CALL WTCFGS(155,ICYTOT)
      ENDIF ! TIMEL .AND. IOLD > 0

C***********************************************************************
      IF(IOLD > 0) RETURN     !go to 8222
C***********************************************************************
      IF(MASTER) THEN  !(Transition model after r372: Minor 2->3) 
*         WRITE(9) 8,3 ! Version 8.3 (Major=8, Minor=3) 
*         WRITE(9) NAME(1:79)
*         WRITE(9) RMACH,ALPHA,RE,0.0,ITURB,0,NRESI,IDRXX ! ITURB, TRANSL ?
         WRITE(9) 8,3        ! Version 8.3 (Major=8, Minor=3) 
         WRITE(9) NAME(1:79)
         WRITE(9) REAL(RMACH,4),REAL(ALPHA,4),REAL(RE,4),REAL(0.0,4),
     +            ITURB,0,NRESI,IDRXX    ! ITURB, TRANSL ?
         WRITE(59) 8,3 ! Version 8.3 (Major=8, Minor=3)
         WRITE(59) NAME(1:79)
*         WRITE(59) RMACH,ALPHA,RE,T+1.,ITURB,0,NRESI,IDRXX ! ITURB,TRANSL
         WRITE(59) REAL(RMACH,4),REAL(ALPHA,4),REAL(RE,4),REAL(T+1.,4),
     +             ITURB,0,NRESI,IDRXX    ! ITURB,TRANSL
      ENDIF

      RETURN
      END SUBROUTINE INPUT

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE INPUT2(NPPN,OMEGAREF,NAME)

      USE NS3CO
      USE INTEGERS, ONLY : IPRO, NPRO, IOTY, NB, NBCS
      USE MAIN_ARRAYS, ONLY : NPROCE, IMAX, JMAX, KMAX, MGRID
      USE CONSTANTS, ONLY : PII
      USE MPI

      IMPLICIT NONE

      INTEGER :: I, J, N, NSXX, NPPN, IVERS
      INTEGER :: ERRORCODE, IERR
      REAL    :: OMEGAREF
      CHARACTER(LEN=80) :: NAME

      LOGICAL :: MASTER

      MASTER = IPRO == 1

      IF(MASTER) THEN
         WRITE(45,*) 'PROCES file read:'
         WRITE(45,*) NPRO,IOTY
         DO N = 1,NPRO
            WRITE(45,*) NPROCE(1,N)
            WRITE(45,*)(NPROCE(J,N),J=2,NPROCE(1,N)+1)
         ENDDO
      ENDIF

      IF(NBLOCK > NB) THEN
         WRITE(*,*) ' NB is too small in file NS3CO in pro' ,IPRO
         WRITE(*,*) 'Should be at least size of the NBLOCK=',NBLOCK
*         CALL EXIT
         STOP
      ENDIF
      IF(NPPN <  NBCS) THEN
         WRITE(*,*) 'NPPN is too small for process',IPRO
         WRITE(*,*) 'Should be at least ',NBCS,' (now',NPPN,')'
         IF(PARALLEL) CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
         STOP
      ENDIF

C ****************** CHECK THAT IPRINT IS RESONABLE ******************

      IF((ICMAX-ICYCLE)/KP > 50) THEN
         WRITE(*,*)'You will have more that 50 prints in OUTPUT file.'
         WRITE(*,*)'If you really want so many prints, you must change'
         WRITE(*,*)'code. But this time: Exiting ...'
         STOP
      ENDIF
C ****************** CHECK THAT KRRINT IS RESONABLE ******************

      IF((ICMAX-ICYCLE)/KRRINT > 10) THEN
         WRITE(*,123)
 123     FORMAT(
     + 'You will have more that 10 RSTART files during the run.'/
     + 'If you really want so many RSTART, you must change code.'/
     +  'But this time: Exiting ...')
         STOP
      ENDIF

C *********************************************************************

      IF (ITURB >= 3 .AND. ITURB < 20) THEN
         IVERS = 3
      ELSEIF (ITURB >= 20) THEN
         IVERS = 4
      ELSEIF (ITURB >= 0 .AND. ITURB < 3) THEN
         IVERS = 1
      ELSE
         IVERS = -1
      ENDIF ! ITURB >= 3 .AND. ITURB < 20

      IF (MULPHL) THEN
         IVERS = 5
      ENDIF ! MULPHL

      IF (IVERS == -1) STOP 'Aborting ...'
      IF (IVERS >= 0 .AND. IVERS <= 1 .AND. IPRO == 1) WRITE(4,2500)
      IF (IVERS == 2 .AND. IPRO == 1) WRITE(4,2600)
      IF (IVERS == 3 .AND. IPRO == 1) WRITE(4,2700)
      IF (IVERS == 4 .AND. IPRO == 1) WRITE(4,2800)
      IF (IVERS == 5 .AND. IPRO == 1) WRITE(4,2900)
      IF (IVERS == 2) WRITE(45,2600)
      IF (IVERS == 3) WRITE(45,2700)
      IF (IVERS == 4) WRITE(45,2800)
      IF (IVERS == 5) WRITE(45,2900)

 2500 FORMAT(/9X,'Calculating problem size for the standard version')
 2600 FORMAT(/9X,'Calculating problem size for the rotational version')
 2700 FORMAT(/9X,'Calculating problem size for the k-epsilon version')
 2800 FORMAT(/9X,'Calculating problem size for the RSM version')
 2900 FORMAT(/9X,'Calculating problem size for the multiphase version')
 

      WRITE(45,2000) LEVEL
      IF(IPRO == 1) WRITE(4,2000) LEVEL
 2000 FORMAT(9X,'The required starting level is',I4,'.'/,9X,
     & 'The grid dimensions on the original finest grid level are:'/)

      DO 200 I=1,NBLOCK
         WRITE(45,1300) I,IMAX(1,I),JMAX(1,I),KMAX(1,I),MGRID(I)
         IF(IPRO == 1) 
     &        WRITE(4,1300)  I,IMAX(1,I),JMAX(1,I),KMAX(1,I),MGRID(I)
 200  CONTINUE
 1300 FORMAT(9X,'BLOCK',I4,3X,'IMAX =',I5,3X,'JMAX =',I5,3X,'KMAX =',I5,
     &       3X,'MGRID =',I3)
      WRITE(45,*)
      IF(IPRO == 1) WRITE(4,*)
      IF(IOLD == 0) WRITE(45,9101) LEVEL
      IF(IOLD < 0) WRITE(45,9102) LEVEL,LEVEL+1
      IF(IOLD > 0) WRITE(45,9103) LEVEL
9101  FORMAT(' CALCULATION IS STARTED FROM LEVEL',I3)
9102  FORMAT(' CALCULATION IS STARTED FROM LEVEL',I3,/,
     + ' INITIAL CONDITION FROM LEVEL',I3)
9103  FORMAT('CALCULATION IS CONTINUED FROM LEVEL',I3)
C **********************************************************************

      IF(MASTER) THEN
      IF(ITURB == 1)WRITE(4,*)' ALGEBRAIC TURBULENCE MODEL IS EMPLOYED
     +'
      IF(ITURB  >= 21)WRITE(4,*)' REYNOLDS STRESS MODELLING IS EMPLOYED'
      IF(ITURB == 3)WRITE(4,*)
     +     ' K-EPSILON TURBULENCE MODEL IS EMPLOYED'
      IF(ITURB == 4)WRITE(4,*)
     +     ' K-EPSILON TURBULENCE MODEL IS EMPLOYED (ITURB = 4)'
      IF(ITURB == 5)WRITE(4,*)
     +     ' K-EPSILON TURBULENCE MODEL IS EMPLOYED (ITURB = 5)'
      IF(ITURB == 6)WRITE(4,*)
     +  ' K-OMEGA MODEL IS EMPLOYED (ITURB = 6)'
      IF(ITURB == 7)WRITE(4,*)
     +  ' K-EPSILON MODEL WITH FUZZY k-omega TUNING (ITURB = 7)'
      IF(ITURB == 8)WRITE(4,*)
     +  ' Smagorinsky SGS model with Cs=0.2 (ITURB = 8)'
      IF(ITURB == 9)WRITE(4,*)
     +  ' Spalart-Allmaras model (ITURB = 9)'
      IF(ITURB  == 10)WRITE(4,*)
     +' EXPLIXIT ASM TURBULENCE MODEL IS EMPLOYED BY GATSKI (ITURB=10)'
      IF(ITURB  == 11)WRITE(4,*)
     +' EXPLIXIT ASM TURBULENCE MODEL BY GATSKI (ITURB=11)'//
     +' APPLIED WITH CHIEN MODEL'
      IF(ITURB >= 21)WRITE(4,*)
     +     ' RSM TURBULENCE MODEL IS EMPLOYED'
      WRITE(4,*)
      IF(NSCAL > 0)WRITE(4,*)' TOTAL NUMBER OF SCALAR EQUATIONS = ',
     + NSCAL
      WRITE(4,*)
      ENDIF           !(IPRO == 1)

      IF(MASTER) THEN
      WRITE(4,*)'******************************************************'
      WRITE(4,*) ' TEMPERATURE =',REAL(FRSTEM,4),' K'
      WRITE(4,*) ' INI. TEMP.  =',REAL(TEMINI,4),' K'
      WRITE(4,*) ' VELOCITY    =',REAL(FRSVEL,4),' m/s'
      WRITE(4,*) ' DENSITY     =',REAL(FRSDEN,4),' kg/m3'
      WRITE(4,*) ' PRESSURE    =',REAL(FRSPRE,4),' N/m2'
      WRITE(4,*) ' INT.ENERGY  =',REAL(FRSSIE,4),' J/kg'
      WRITE(4,*) ' INI.INT.En. =',REAL(SIEINI,4),' J/kg'
      WRITE(4,*) ' SOUND SPEED =',REAL(FRSSSP,4),' m/s'
      WRITE(4,*) ' VISCOSITY   =',REAL(FRSVIS,4),' kg/ms'
      WRITE(4,*) ' REF. LENGTH =',REAL(CHLREF,4),' m'
      WRITE(4,*) ' REF. AREA   =',REAL(AREF,4),  ' m2'
      ENDIF ! IPRO == 1

      IF(MASTER .AND. ITURB >= 3) THEN
      WRITE(4,*) ' TURB. VIS.  =',REAL(FRSMUT,4),' kg/ms'
      WRITE(4,*) ' INIT. T.VIS.=',REAL(RMUINI,4),' kg/ms'
      ENDIF


      IF(MASTER) THEN
      WRITE(4,*) ' STAG. TEMP. =',REAL(T0,4),    ' K'
      WRITE(4,*) ' MACH NUMBER =',REAL(RMACH,4)
      WRITE(4,*) ' REYNOLDS    =',REAL(RE,4)
      WRITE(4,*) ' ALPHA       =',REAL(ALPHA*180./PII,4), ' degrees'
      WRITE(4,*) ' BETA        =',REAL(BETA*180./PII,4),  ' degrees'

      IF(ADVANJ /= 0.) WRITE(4,*) ' Advance number J =',ADVANJ

      IF(FRESUL) WRITE(4,*) ' FROUDE      =',
     +                        REAL(FROUDE,4),FRESUL
      WRITE(4,*)
      ENDIF ! IPRO == 1

      IF(MASTER) THEN
         WRITE(4,*) 'REFERENCE PRESSURE USED FOR NONDIMENSIONAL ',
     +        'AERODYNAMIC COEFFICIENTS IS'
         WRITE(4,*) ' REFPRE      =',REAL(REFPRE,4),'N/m2'
         WRITE(4,*) 'REFERENCE PRESSURE USED TO CALCULATE PRESSURE',
     +        ' FORCES AT THE WALLS IS'
         WRITE(4,*) ' DIFPRE      =',REAL(DIFPRE,4),'N/m2'
         IF(GRAVIL) THEN
            WRITE(4,*) 'REFERENCE GROUND LEVEL AND GRAVITY VECTOR ARE'
            WRITE(4,*) ' GROUND      =',REAL(GROUND,4),' m'
            WRITE(4,*) ' G           = (',REAL(GX,4),REAL(GY,4),
     +                 REAL(GZ,4),') m/s2'
            IF(CAVLEV /= 0.) THEN
            WRITE(4,*) 'REFERENCE CAVITATION LEVEL IS'
            WRITE(4,*) ' CAVLEV      =',REAL(CAVLEV,4),' m'
            ENDIF
         ENDIF
      WRITE(4,*)
      ENDIF

      IF(MASTER) THEN
         WRITE(4,*) 'DEFAULT INLET (OR FREE STREAM) VALUES AFTER INPUT'
         WRITE(4,*) '*************************************************'
         WRITE(4,*) ' FRSDEN =',REAL(FRSDEN,4),' kg/m3'
         WRITE(4,*) ' FRSPRE =',REAL(FRSPRE,4),' N/m2'
         WRITE(4,*) ' FRSVEL =',REAL(FRSVEL,4),' m/s'
         WRITE(4,*) ' REFVEL =',REAL(REFVEL,4),' m/s'
         WRITE(4,*) ' GVEX, GVEY, GVEZ =',
     +                REAL(GVEX,4),REAL(GVEY,4),REAL(GVEZ,4),' m/s'
         WRITE(4,*) ' TOTAL ENERGY(no turbulence)', 
     +                REAL(FRSDEN*(FRSSIE + .5*FRSVEL**2),4),' J/kg'
            WRITE(4,*) ' TLOLIM = ',REAL(TLOLIM,4), ' K'
            WRITE(4,*) ' TUPLIM = ',REAL(TUPLIM,4), ' K'
         IF(ITURB >= 3) THEN
            WRITE(4,*) ' FRSMUT =' ,REAL(FRSMUT,4), ' kg/ms'
            WRITE(4,*) ' RMUINI =' ,REAL(RMUINI,4), ' kg/ms'
            WRITE(4,*) ' RKLIM  = ',REAL(RKLIM,4),  ' J/m3'
            WRITE(4,*) ' FRSLEN = ',REAL(FRSLEN,4), ' m'
            IF(ITURB <= 5) THEN
            WRITE(4,*) ' EPSLIM = ',REAL(EPSLIM,4),' J/m3s'
            ELSEIF(ITURB == 9) THEN
            WRITE(4,*) ' RNUTLIM  =',REAL(EPSLIM,4),'kg/ms'
            ELSE
            WRITE(4,*) ' OMEGALIM =',REAL(EPSLIM,4),'kg/m3s'
            ENDIF
         ENDIF
         WRITE(4,*) ' T0     =',REAL(T0,4)
         IF(MULPHL) THEN
            WRITE(4,*) ' FRSALF =',REAL(FRSALF,4), '      '
         ENDIF
         WRITE(4,*)
         WRITE(4,9400) REAL(GX,4),REAL(GY,4),REAL(GZ,4)
         WRITE(4,9410)
      ENDIF

c      FRSE = FRSDEN*(FRSSIE + .5*FRSVEL**2)

3001  FORMAT(/' Coefficient history can be activated in MAIN and'
     + ' NS3D (search for: 3000 FORMAT)'/)
9400  FORMAT(' GRAVITY VECTOR: g(x) = ',F7.4,', g(y) = ',F7.4,
     + ', g(z) = ',F7.4)
9410  FORMAT(' **************')
9420  FORMAT('Courant numbers: CFL = ',F5.2,', CFLL =',F5.2)
9480  FORMAT(80('='))

C ********************************************************************

 3000 FORMAT(/'HISTORY OF FORCE COEFFICIENTS'/
     +     '*****************************'/
     +     '# ICYCLE',5X,'CX',11X,'CY',11X,'CZ',11X,'MX',11X,'MY',11X,
     +     'MZ',11X,'TOM',10X,'Mass f. out',1X,'Mech f. out',1X,
     +     'effici',7X,'Head'/
     +     '# units',7X,'N',12X,'N',12X,'N',11X,'Nm',11X,'Nm',11X,
     +     'Nm',13X,'W',9X,'kg/s',9X,'W',11X,
     +     ' %',13X,'m')

C ... WRITE ON THE SEAGULL FILE

      IF(MASTER) THEN
         NSXX = NSCAL
         IF(ITURB >= 20) NSXX = NSCAL - 6
         WRITE(17,9944) NAME
         WRITE(17,*) REAL(RMACH,4),REAL(RE,4)
         WRITE(17,*) REAL(1.,4),REAL(1.,4),REAL(1.,4),REAL(1.,4),
     +               REAL(1.,4),REAL(OMEGAREF,4),REAL(CHLREF,4)
         WRITE(17,*) REAL(FRSPRE,4),REAL(FRSDEN,4),REAL(FRSVEL,4)
         WRITE(17,*) REAL(GAMMA,4),REAL(PR,4),REAL(PRT,4)
         WRITE(17,*) ITURB,ISTATE
         WRITE(17,*) REAL(RGAS,4),REAL(VISU0,4),REAL(EXPSU,4),
     +               REAL(TSU0,4)
         WRITE(17,*) REAL(E0REF,4),REAL(T0REF,4)
         WRITE(17,*) NSXX,.FALSE.
         WRITE(17,*) REAL(GROUND,4),
     +               REAL(GX,4),REAL(GY,4),REAL(GZ,4),REAL(G0,4)
         CLOSE(17)
      ENDIF 

9944  FORMAT(A79)

      RETURN
      END SUBROUTINE INPUT2

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE COMPUT(IBCOM,IBCOM2,IBSLB,MB,MS,MGSTP,MBSTP,NPNTS,
     +                  IBF,MGM)

      USE NS3CO
      USE CONSTANTS, ONLY : EPS
      USE MAIN_ARRAYS, ONLY : IR, IQ, NPROCE, IMAX, JMAX, KMAX, IT, IK,
     &                        IDI1, IDI2, IDI3, MGRID, IG, JF, IH, BLKS,
     &                        MGRIDA, ICON, IHF, JSTATE, IGRID,
     &                        OMEGAX, OMEGAY, OMEGAZ, OMEGA,
     &                        CENAX, CENAY, CENAZ,
     &                        INTERI, INTERJ, INTERK
      USE INTEGERS, ONLY : IPRO, NBCS
      USE FLIGHT, ONLY : OSKU, XCG, YCG, ZCG, XCGI, YCGI, ZCGI,
     &                   PSIR, THETAR, PHIR, PSIRI, THETARI, PHIRI, 
     &                   NGRIFL
      USE MPI

      IMPLICIT NONE

      INTEGER :: IBCOM(*), IBCOM2(*), IBSLB(*), MB, MS, MGSTP, MBSTP,
     &           IBF, MGM, NPNTS  ! , IPER
      INTEGER :: MININD, ERRORCODE, IERR, MGMIN, IP
      INTEGER :: NEXTB, NEXTS, NEXT2 ,N, NGL, M, LVLSTP ,KVLSTP,
     &           KMININ, IIX, IIY, IIZ, NM ,NSTP, INCGR, INCSL, MGMAX,
     &           INUM1, JNUM1, INCLR, II
      REAL    :: FRSE

C
C     MB     = SUURIN LOHKOKOKO (VAIN YKSI TASO)
C     MS     = SUURIN SLABIKOKO (VAIN YKSI TASO)
C     IBCOM  = LOHKON KAIKKIEN TASOJEN TILAVARAUS
C     IBCOM2 = LOHKON ALEMPIEN TASOJEN TILAVARAUS
C     IBSLB  = LOHKON KAIKKIEN MAKSIMISLABIEN TILAVARAUS
C
C ... CALCULATE IMAX JMAX KMAX  IN PARTICULAR LEVEL AND THE MULTIGRID
C ...           ARRAY POINTERS
C ...           IT,IK,
C ...           IDI1,IDI2,IDI3

      MB = 0
      MS = 0

      NEXTB = 1
      NEXTS = 1
      NEXT2 = 1
      NPNTS = 1

      WRITE(45,1950)

      DO 10 N = 1,NBLOCK
         DO 11 M = 1,MGM
            IR(M,N) = 1
            IQ(M,N) = 1
 11      CONTINUE
 10   CONTINUE

      IR(1,NBLOCK+1) = 1
      IQ(1,NBLOCK+1) = 1

      DO 100 N=1,NBLOCK
         IBCOM(N)  = 0
         IBCOM2(N) = 0
         IBSLB(N)  = 0
         NGL       = NPROCE(N+1,IPRO)
         MININD    = 1
         NPNTS     = MAX0(NPNTS,(IMAX(1,N)+1)*(JMAX(1,N)+1)*
     &        (KMAX(1,N)+1))
         LVLSTP    = 2**(LEVEL-1)
C ... FOR 2-D CASES                           !PPR
         KVLSTP    = LVLSTP                   !PPR
         KMININ    = MININD                   !PPR
         IF ((KMAX(1,N)/KVLSTP <= 1).AND.KMAX(1,N) /= 2) THEN !PPR
            KVLSTP = KMAX(1,N)                !PPR
            KMININ = 1                        !PPR
         ENDIF                                !PPR
         IIX       = MAX0(MININD,IMAX(1,N)/LVLSTP)
         IIY       = MAX0(MININD,JMAX(1,N)/LVLSTP)
         IIZ       = MAX0(KMININ,KMAX(1,N)/KVLSTP) !PPR

         IF(IMAX(1,N) == IIX*LVLSTP .AND. JMAX(1,N) == IIY*LVLSTP .AND.
     &      KMAX(1,N) == IIZ*KVLSTP .AND. N <= NBLOCK) THEN

            IMAX(1,N) = IIX
            JMAX(1,N) = IIY
            KMAX(1,N) = IIZ
            IMAX(2,N) = IIX/2
            JMAX(2,N) = IIY/2
            KMAX(2,N) = MAX(1,IIZ/2)

            IT(1,N)   = (IT(1,N)-1+LVLSTP)/LVLSTP
            IK(1,N)   = (IK(1,N)-1+LVLSTP)/LVLSTP
            IDI1(N)   = (IDI1(N)-1+LVLSTP)/LVLSTP
            IDI2(N)   = (IDI2(N)-1+LVLSTP)/LVLSTP
            IDI3(N)   = (IDI3(N)-1+KVLSTP)/KVLSTP !PPR
            MGRID(N)  = MIN0(MAX0(1,MGRID(N)),MGM)
            IF(N == 1) WRITE(45,9201) LEVEL
            WRITE(45,*)
            WRITE(45,*) 'BLOCK =',N
            WRITE(45,*)
            WRITE(45,*) 'IMAX  =',IMAX(1,N),' JMAX =',JMAX(1,N),
     +           ' KMAX =',KMAX(1,N)
            WRITE(45,*) 'IT    =',IT(1,N),  ' IK   =',IK(1,N)
            WRITE(45,*) 'IDI1  =',IDI1(N),' IDI2 =',IDI2(N),' IDI3 =',
     +           IDI3(N)
            WRITE(45,*) 'MGRID =',MGRID(N)
            WRITE(45,*)
            WRITE(45,9202)

9201  FORMAT(' NEW INDICES FOR LEVEL',I3,' (SEE SUBROUTINE SETLEV):')
9202  FORMAT(' OBS! NUMBER OF GRID LEVELS IS THE SAME AS IN INPUT FILE')

         ELSE
            WRITE(*,2000) NPROCE(1+N,IPRO),LEVEL
            WRITE(4,2000) NPROCE(1+N,IPRO),LEVEL
            WRITE(*,*) 'IMAX  =',IMAX(1,N),' JMAX =',JMAX(1,N),
     +           ' KMAX =',KMAX(1,N),' IPRO =',IPRO,N
            STOP
         ENDIF

         DO 200 NM=1,MGRID(N)
            NSTP = 2**(NM-1)
            IIX = MAX0(MININD,IMAX(1,N)/NSTP)
            IIY = MAX0(MININD,JMAX(1,N)/NSTP)
            IIZ = MAX0(KMININ,KMAX(1,N)/NSTP)   !PPR
            IF (IMAX(1,N) == IIX*NSTP .AND. JMAX(1,N) == IIY*NSTP .AND.
     &           KMAX(1,N) == IIZ*NSTP) THEN
               IMAX(NM,N) = IIX
               JMAX(NM,N) = IIY
               KMAX(NM,N) = IIZ
               WRITE(45,2100) N,IMAX(NM,N),JMAX(NM,N),KMAX(NM,N),NM
            ELSEIF(IMAX(1,N) == IIX*NSTP .AND. JMAX(1,N) == IIY*NSTP 
     &         .AND. IIZ == 1 .AND. 
     &         (KMAX(1,N) ==  1 .OR. KMAX(1,N) ==  2 .OR.
     &          KMAX(1,N) ==  4 .OR. KMAX(1,N) ==  8 .OR.
     &          KMAX(1,N) == 16 .OR. KMAX(1,N) == 32)) THEN  ! 2D
               IMAX(NM,N) = IIX                ! PPR
               JMAX(NM,N) = IIY                ! PPR
               KMAX(NM,N) = IIZ                ! PPR
               WRITE(45,2100) N,IMAX(NM,N),JMAX(NM,N),KMAX(NM,N),NM
               WRITE(45,2101) N,NM             ! PPR
            ELSE
               WRITE(*,2050)  NPROCE(1+N,IPRO),MGRID(N)
               WRITE(45,2050) NPROCE(1+N,IPRO),MGRID(N)
               WRITE(*,*) 'IMAX  =',IMAX(1,N),' JMAX =',JMAX(1,N),
     &              ' KMAX =',KMAX(1,N),' IPRO =',IPRO,N
               IF(PARALLEL)CALL MPI_ABORT(MPI_COMM_WORLD,ERRORCODE,IERR)
               STOP
            ENDIF

            INCGR = MGSTP + (IMAX(NM,N) + 2*IN)*(JMAX(NM,N) + 2*JN)*
     &           (KMAX(NM,N) + 2*KN)
            INCSL = MGSTP + MAX0(
     &           (IMAX(NM,N) + 2*IN)*(JMAX(NM,N) + 2*JN),
     &           (JMAX(NM,N) + 2*JN)*(KMAX(NM,N) + 2*KN),
     &           (KMAX(NM,N) + 2*KN)*(IMAX(NM,N) + 2*IN))
            IBCOM(N) = IBCOM(N) + INCGR
            IBSLB(N) = IBSLB(N) + INCSL
            MB       = MAX0(MB,INCGR)
            MS       = MAX0(MS,INCSL)
            IG(NM,N) = NEXTB
            JF(NM,N) = NEXTS
            IH(NM,N) = NEXT2
            IF(ITURB >= 3) IR(NM,N) = IG(NM,N)
            IF(NSCAL > 0) IQ(NM,N) = IG(NM,N)
            NEXTB = NEXTB + INCGR
            NEXTS = NEXTS + INCSL
            NEXT2 = NEXT2 + MIN(1,NM-1)*INCGR
 200     CONTINUE
         IBCOM(N)  = IBCOM(N) + MBSTP
         IBCOM2(N) = IBCOM(N) - ((IMAX(1,N) + 2*IN)*
     &        (JMAX(1,N) + 2*JN)*(KMAX(1,N) + 2*KN) + MGSTP)
         IBSLB(N)  = IBSLB(N) + MBSTP
         NEXTB     = NEXTB + MBSTP
         NEXTS     = NEXTS + MBSTP
         NEXT2     = NEXT2 + MBSTP
         IG(1,N+1) = NEXTB
         JF(1,N+1) = NEXTS
         IH(2,N+1) = NEXT2
         IF(ITURB >= 3) IR(1,N+1) = IG(1,N+1)
         IF(NSCAL > 0) IQ(1,N+1) = IG(1,N+1)

      IF(BLKS(NGL)%SOLUTION_TYPE == 'SOLID' .OR.
     &   IPRESC == 1) THEN
         WRITE(45,*)
         WRITE(45,*) 'This block is solved with AMG. Search for ',
     &   'the number of AMG levels:'
         DO 300 NM=1,6
            NSTP = 2**(NM-1)
            IIX = MAX0(MININD+1,IMAX(1,N)/NSTP)
            IIY = MAX0(MININD+1,JMAX(1,N)/NSTP)
            IIZ = MAX0(KMININ,  KMAX(1,N)/NSTP)!PPR
            IF(IMAX(1,N) == IIX*NSTP .AND. JMAX(1,N) == IIY*NSTP .AND.
     &         KMAX(1,N) == IIZ*NSTP) THEN
               MGRIDA(N) = NM
            ELSEIF(IMAX(1,N) == IIX*NSTP .AND. JMAX(1,N) == IIY*NSTP 
     &         .AND. IIZ == 1 .AND. 
     &         (KMAX(1,N) ==  1 .OR. KMAX(1,N) ==  2 .OR.
     &          KMAX(1,N) ==  4 .OR. KMAX(1,N) ==  8 .OR.
     &          KMAX(1,N) == 16 .OR. KMAX(1,N) == 32)) THEN  ! 2D
               MGRIDA(N) = NM
               WRITE(45,2103) N,NM             ! PPR
            ENDIF
 300     CONTINUE
         WRITE(45,2104) N,MGRIDA(N),ITERMA,ITERHA,MCYCLE

      ENDIF    ! SOLUTION_TYPE == SOLID

 2103 FORMAT(1X,'Warning!! Problem is treated as a 2D-case at block',I2,
     &       ' as multigrid is',I2) !PPR
 2104 FORMAT(9X,'BLOCK',I4,3X,'MGRIDA = ',I1,' ITERMA = ',I2,
     &       ' ITERHA = ',I2,' MCYCLE = ',I2)
 100  CONTINUE ! Block loop ended
	
      WRITE(45,*)'******************************************************
     &*****************'
      WRITE(45,*)' ARRAY POINTERS:'
      DO 7650 N = 1,NBLOCK
         WRITE(45,*) ' BLOCK = ',N,' IG = ',(IG(M,N),M=1,MGRID(N))
         WRITE(45,*) ' BLOCK = ',N,' IQ = ',(IQ(M,N),M=1,MGRID(N))
         WRITE(45,*) ' BLOCK = ',N,' JF = ',(JF(M,N),M=1,MGRID(N))
         WRITE(45,*) ' BLOCK = ',N,' IH = ',(IH(M,N),M=1,MGRID(N))
         WRITE(45,*)
 7650 CONTINUE
 7651 FORMAT(2X,A8,I4,A5,5I10)

      MGMIN = MGM
      MGMAX = 0
      DO N = 1,NBLOCK
         MGMIN = MIN0(MGMIN,MGRID(N))
         MGMAX = MAX0(MGMAX,MGRID(N))
      ENDDO

C ... Calculate index array for all boundary condition patches.

      NEXTS = 1

C ... Space for BC patches and starting addresses. 
      
      DO NM=1,MGMAX

         NSTP = 2**(NM-1)*LVLSTP

         DO IP = 1,NBCS
            INUM1 = ICON((IP-1)*IC9 +5)-ICON((IP-1)*IC9 +4) + 1
            JNUM1 = ICON((IP-1)*IC9 +7)-ICON((IP-1)*IC9 +6) + 1 
            INCLR = (MAX(INUM1/NSTP,1)+2*LN)*(MAX(JNUM1/NSTP,1)+2*LN)
            IHF(IP,NM) = NEXTS
            NEXTS   = NEXTS + INCLR
         ENDDO

         IHF(1,NM+1) = IHF(NBCS,NM) + INCLR

      ENDDO

      IBF = NEXTS

						WRITE(45,*)
      WRITE(45,*)' ARRAY POINTERS for pathches,  IHF:'
      DO NM = 1,MGMAX
						WRITE(45,*)
      WRITE(45,*) 'Level =',NM
						WRITE(45,*)
      DO IP = 1,NBCS,5
         WRITE(45,9000) IP,MIN(NBCS,IP+4),
     +        (IHF(II,NM),II=IP,MIN(NBCS,IP+4))
      ENDDO
      ENDDO
c 9000 FORMAT('IHF(',I4,' ...',I4,') =  ',5I7)
 9000 FORMAT('IHF(',I7,' ...',I7,') =  ',5I7)
	     WRITE(45,*)
      WRITE(45,*) 'IBF size is (including lower grid levels)', IBF
      WRITE(45,*)

      FRSE = FRSDEN*(FRSSIE + .5*FRSVEL**2)

      CALL BLOCK_REPORT(NBLOCK,NPROCE,NBLOCG,IPRO,BLKS,FRSDEN,
     +     FRSVEL,FRSTEM,FRSPRE,FRSE,FRSMUT,RMUINI,T0,TLOLIM,
     +     TUPLIM,RKLIM,FRSLEN,EPSLIM,TEMINI,ITURB,JSTATE,FRSSSP,
     +     ARTSSP,FRSVIS,FRSSIE,SIEINI,ROTAT,ROTANG,OMEGA,OMEGAX,
     +     OMEGAY,OMEGAZ,CENAX,CENAY,CENAZ,MGRID,INTERI,INTERJ,
     +     INTERK,MGRIDA,MULPHL,IGRID,FRSVEL,IPRESC,GROUND,CAVLEV,
     +     TWO_FLUIDL)

      IF(CHIMEL) THEN
         IF(INCHIML) THEN
         WRITE(4,'(1X,2A,L1)') ' First-order interpolation is applied',
     +    'in the Chimera frontier. INCHIML: ',INCHIML
         ELSE
         WRITE(4,'(1X,2A,L1)') ' Default interpolation is applied',
     +    'in the Chimera frontier. INCHIML: ',INCHIML
         ENDIF
      ENDIF

      IF(REFLECL) THEN
         WRITE(4,'(1X,2A,L1)') ' Surface velocities are reflected',
     +    ' from the domain. REFLECL: ',REFLECL
      ENDIF

      IF(ENTROPY_FIX .AND. IFLUX == 1) THEN
         WRITE(4,'(1X,2A,L1)') ' Entropy fix is made in FLUXKE',
     +    ' (Roe). ENTROPY_FIX: ',ENTROPY_FIX
      ELSEIF(.NOT.ENTROPY_FIX .AND. IFLUX == 1) THEN
         WRITE(4,'(1X,2A,L1)') ' Entropy fix is not made in FLUXKE,',
     +    ' in a case of Roe flux. ENTROPY_FIX: ',ENTROPY_FIX
      ENDIF

      IF(TUR_FRS_SOURCE) THEN
         WRITE(4,'(1X,2A,L1)') ' Turbulence sources are corrected',
     +   ' in order to preserve the free-stream values: ',TUR_FRS_SOURCE
      ELSEIF(.NOT.TUR_FRS_SOURCE) THEN
         WRITE(4,'(1X,2A,L1)') ' Turbulence sources are not corrected',
     +   ' in order to preserve the free-stream values: ',TUR_FRS_SOURCE
      ENDIF

      IF(TUR_MULPHL) THEN
         WRITE(4,'(1X,2A,L1)') ' Eddy viscosity is modified',
     +   ' because of the bubble motion: ',               TUR_MULPHL
      ELSEIF(.NOT.TUR_MULPHL .AND. MULPHL) THEN
         WRITE(4,'(1X,2A,L1)') ' Eddy viscosity is not modified',
     +   ' in order to take into account the bubbles: ',  TUR_MULPHL
      ENDIF

       
      IF(LDISTL .AND. TIMEL) THEN
         WRITE(4,*)  ' Wall distances are recalculated at every',
     +    ' time step'
         WRITE(13,*) ' Wall distances are recalculated at every',
     +    ' time step'
      ELSE IF (.NOT.LDISTL .AND. TIMEL) THEN
         WRITE(4,*)  ' Wall distances are not recalculated.'
         WRITE(13,*) ' Wall distances are not recalculated.'
      ELSE IF(LDISTL .AND. .NOT.TIMEL) THEN
         WRITE(4,*)  ' LDISTL is true, but that does affect to the',
     +    ' steady-state simulation'
         WRITE(13,*) ' LDISTL is true, but that does affect to the',
     +    ' steady-state simulation'
      ENDIF

      IF(INTERTU /= 0) THEN
          WRITE(4,*) ' TURBULENCE VARIABLES INTERPOLATED ACCORDING',
     +               ' TO INTER=',INTERTU
      ENDIF

      IF(.NOT.FRESUL) THEN
        WRITE(4,9420) CFL,CFLL
      ELSE
        WRITE(4,9421) CFL,CFLL,CFLFRE
      ENDIF
      WRITE(4,9422) CDIFF,CDIFFT
      WRITE(4,9423) SMAX
      CDIFF  = 1./CDIFF; CDIFFT = 1./CDIFFT ! The simplest way, 2 x mersu

      IF(IPRESC == 1) WRITE(4,9430) ALFAP
      IF(IPRESC == 1) WRITE(4,9428)
      IF(IPRESC == 1) WRITE(4,9429)
      IF(IPRESC == 1) WRITE(4,9435)
      IF(IPRESC == 1) WRITE(4,9434)
      IF(IPRESC == 1) WRITE(4,9431) PSEUCO
      IF(IPRESC == 2) WRITE(4,9432) PSEUCO
      IF(RJK2 > EPS .OR. RJK4 > EPS) WRITE(4,9433) RJK2,RJK4

      IF(NGRIFL > 0) CALL FLIGHT_REPORT(OSKU,XCG,YCG,ZCG,
     +     PSIR,THETAR,PHIR,NGRIFL,XCGI,YCGI,ZCGI,PSIRI,THETARI,PHIRI)

      IF (FRESUL) THEN

      WRITE(4,"(/A/1X,36('*')/)")' Free-surface deformation parameters:'
      WRITE(45,7198)

         IF(IFSBC == 1) THEN
           WRITE( 4,7199) IFSBC,JFIRST
           WRITE(45,7199) IFSBC,JFIRST
         ELSE
           WRITE( 4,7200) IFSBC
           WRITE(45,7200) IFSBC
           WRITE( 4,7201) FREDIF
           WRITE(45,7201) FREDIF
           WRITE( 4,7202) DTWMIN,DTWMAX
           WRITE(45,7202) DTWMIN,DTWMAX
           WRITE( 4,7203) DWMV
           WRITE(45,7203) DWMV
         ENDIF

         WRITE( 4,7197)
         WRITE(45,7197)
      ENDIF ! FRESUL
 7197 FORMAT(/' Grid manipulation according to IGRID'/)
 7198 FORMAT(/' Free-surface deformation parameters:'/)
 7199 FORMAT('  IFSBC = ',I2,'  JFIRST = ',I2)
 7200 FORMAT('  IFSBC  = ',I2)
 7201 FORMAT('  FREDIF = ',E12.5,' m2/s')
 7202 FORMAT('  Time-step limits: DTWMIN = ',E10.3,' s, DTWMAX = ',E11.4
     &       ' s')
 7203 FORMAT('  Maximum change in wave height: DWMV = ',E11.4,' m')

      WRITE(4,3001)

      RETURN

 1950 FORMAT(/'COMPUT : The grid dimensions of the case are:'/)
 2000 FORMAT(9X,'Grid Uniformity is not conserved in block',I4,'.'/
     &        9X,'Cannot start from level ',I2,'.')
 2050 FORMAT(9X,'Grid Uniformity is not conserved in block',I4,'.'/
     &        9X,'Cannot use ',I2,' multigrid levels.')
 2100 FORMAT(9X,'BLOCK',I4,3X,'IMAX =',I5,3X,'JMAX =',I5,3X,'KMAX =',I5,
     &        3X,'LEVEL =',I2)
 2101 FORMAT(1X,'Warning!! Problem is treated as a 2D-case at block',I2,
     &       ' and multigrid is',I2) !PPR

 3001 FORMAT(/' Coefficient history can be activated in MAIN and'
     & ' NS3D (search for: 3000 FORMAT)'/)
 9420 FORMAT('  Courant numbers: CFL = ',F5.2,', CFLL =',F5.2)
 9421 FORMAT('  Courant numbers: CFL = ',F5.2,', CFLL =',F5.2,
     &       ', CFLFRE =',F5.2)
 9422  FORMAT('  Diffusion criteria: CDIFF = ',F7.2,', CDIFFT =',F7.2)
 9423  FORMAT('  In case of solids heat conduction: SMAX = ',F7.2/)

 9428 FORMAT('  Pressure-velocity coupling according to PREVEL above.',/
     &'  PREVEL = 2 is SIMPLE, PREVEL = 1 its variant, PREVEL = 3 is',
     &' SIMPLEC.')
 9429 FORMAT('  Mass-flow correction according to FLUXCORR above.',
     &' FLUXCORR = T  =>',/'  Mass flows are corrected on the cell',
     &' surfaces after the pressure correction.')
 9435 FORMAT('  R-C damping at I/O surfaces according to *RC above.',/,
     &'  INLRC = 0 RC damping is included, INLRC = 1 it is not,',/,
     &'  INLRC = 2 is one-sided interpolation without RC.')
 9430 FORMAT('  Global underrelaxation for pressure correction: ALFAP ='
     & ,1X,F5.3,/,'  is applied in this code version.')
 9431  FORMAT('  MAX[ARTSSP, PSEUCO * REFVEL], (PSEUCO = ',
     & F5.2,')'/,'  is applied for Rhie&Chow term in flux calculation.')
 9432  FORMAT('  Pseudocompressibility parameter: PSEUCO = '
     & ,F5.2,/,'  is applied for Rhie&Chow term in flux calculation',/,
     & '  and for the artificial compressibility method.')
 9433 FORMAT('  A Jameson type dissipation will be added to fluxes.',
     & /,'  K2 =',E10.3,' K4 =',E10.3)
 9434 FORMAT('  Skewness correction MIN(INPUT,values in MEMORY-file)')

      END SUBROUTINE COMPUT

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE FINAL(IBCOM,IBCOM2,IBSLB,MB,MS,MAXB,MAX11,MAX2,
     &     IB,MAXWNE,MAWNE,NBMAX)

C ... MB     = SUURIN LOHKOKOKO (VAIN YKSI TASO)
C ... MS     = SUURIN SLABIKOKO (VAIN YKSI TASO)

      USE NS3CO, ONLY : NBLOCK
      USE MAIN_ARRAYS, ONLY : JF, IG, IH

      IMPLICIT NONE

      INTEGER :: IBCOM(*), IBCOM2(*), IBSLB(*)
      INTEGER :: MB, MS, MAXB, MAX11, MAX2, IB, MAXWNE, MAWNE, NBMAX, N

      MAXB = 1                  !whole memory array for one variable
      MAX2 = 1                  !memory array for 2. grid level
      MAX11= 1                  !memory array for the biggest block
      NBMAX= 1                  !number of biggest block
      IB   = 1
      MAWNE  = 0

      DO 100 N = 1 , NBLOCK

         MAXB = MAXB + IBCOM(N)
         MAX2 = MAX2 + IBCOM2(N)
         IF(MAX11 < IBCOM(N)) NBMAX = N
         MAX11= MAX0(MAX11,IBCOM(N))
         IB   = IB + IBSLB(N)
         WRITE(45,1100) N,IBCOM(N),IBCOM2(N),IBSLB(N)

 100  CONTINUE

      WRITE(45,*)
      CALL MODIMS(MS)
      MAXWNE = MAX0(MB,46*MS) ! FLUMPH
      MAWNE  = MS

C ********** THIS CHECK IS UNNECCESERY IN THIS VERSION *****************
      WRITE(45,*)
      WRITE(45,*)'PROBLEM SIZE: IG = ',IG(1,NBLOCK+1),' JF = ',
     + JF(1,NBLOCK+1),' IH = ',IH(2,NBLOCK+1)
      WRITE(45,*)'COMMON BLOCK: IG = ',MAXB,' JF = ',
     + IB,' IH = ',MAX2
      IF(MAXB < IG(1,NBLOCK+1) .OR. IB < JF(1,NBLOCK+1) .OR.
     + MAX2 < IH(2,NBLOCK+1)) THEN
          WRITE(45,*)'PROBLEM SIZE TOO BIG. INCREASE COMMON BLOCK'
          WRITE(*,*)'TOO BIG. SEE MEMORY-FILE FOR INSTRUCTIONS'
C         CALL EXIT
          STOP
      ENDIF
      WRITE(45,*) 
      WRITE(45,*) ' Biggest block size and the number, MAX11=',MAX11,
     +     ' NBMAX=',NBMAX
      WRITE(45,*)
C ***********************************************************************

C ... MAXB is modified for DOUBLE PRECISION working arrays

      IF(MOD(MAXB,2) /= 0) MAXB = MAXB + 1

C ... MAXB is modified if it are more prime factor of 2 than 1
C ... Can speed up up to 10 % execution time

      IF(MOD(MAXB,4) == 0) MAXB = MAXB + 2

      RETURN

 1100 FORMAT(9X,'BLOCK NR',I4,3X,'IBCOM =',I8,3X,'IBCOM2 =',I8,3X,
     &          'IBSLB =',I6)

      END SUBROUTINE FINAL

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE SET_RINP(INPARA,INPAR,RNPARA,IRNPAR,LNPARA,LRNPAR,
     + ISCALA,INSCAL,RSCALA,IRSCAL,CNPARA,ICPAR,KSTATE,IFLUX,IX,IY,
     + IZ,NAME,NPHASES,FNPARA,CFPARA,IFNPAR,NGRIFL,INFPAR,INFNPAR)

      USE CONSTANTS, ONLY : PII,EPS,DEG2RAD,RAD2DEG

      USE NS3CO,        ONLY : NBGGG,RGAS,TSU0,VISU0,EXPSU,GAMMA,
     +                         E0REF,T0REF, NBLOCK,PARALLEL,
     +                         GVEX,GVEY,GVEZ,GX,GY,GZ,GROUND,G0,
     +                         ALTITUDE,START_CAVL,FEMXYZ,AUTOCONL,
     +                         RMESHT
      
      USE FORCE_GROUPS, ONLY : FORCE_GROUP_INFO, 
     +                         XMOMR, YMOMR, ZMOMR,
     +                         XMOMA, YMOMA, ZMOMA,
     +                         FORCE_GROUP_SHORT_NAME,
     +                         FORCE_GROUP_FULL_NAME

      USE MPCCIVARS

      IMPLICIT NONE

      CHARACTER(LEN=80) :: NAME
      CHARACTER(LEN=80) :: STARTC,STRESC,SOURC,STATEC,TIMEC,CONVC,
     +        FULLNC,FLUXTY,MULPHC,TURCOC,XXTRAC,FRESUC,DISTANC,
     +        TRUE_DISTC,FPRINTC,TURDESC,INCHIMC,REFLECC,ENTROPY_FIXC,
     +        TUR_FRS_SOURCEC,TRANSC,TUR_MULPHC,NEGVC,MPCCIC,MODALFSIC,
     +        LONG_STARTC,AVER_STARTC,START_CAVC,USE_QUATERN,WALLFUNC,
     +        AUTOCONC
      INTEGER :: INPAR,IRNPAR,LRNPAR,N,IIMAX,JJMAX,KKMAX,IFLUX,
     +        INSCAL,IRSCAL,ICPAR,KSTATE,NPHASES,IFNPAR,INFNPAR
      INTEGER :: IX(*),IY(*),IZ(*),INPARA(INPAR,*),ISCALA(INSCAL),
     +        INFPAR(INFNPAR,*)
      REAL    :: RNPARA(IRNPAR,*),RSCALA(IRSCAL),FNPARA(IFNPAR,*)
      LOGICAL :: LNPARA(LRNPAR)  ! ,THERE
      REAL    :: XAU32,RNORM,ZCGD

      CHARACTER(LEN=80) :: CNPARA(ICPAR,*)
      CHARACTER(LEN=1)  :: CFPARA(*)

      INTEGER :: LEVEL,ISTRES,KSCAL,IPRESC,LUSGS,IROTCO,IDRXX,
     + ICMAX,KP,MPRINT,KRRINT,ISTATE,NPHASE,IEPSMA,JRDIF,JRDIS,
     + JRPRE,JRIMP,ITURB,KOVER,
     + MCYCLE,ITERMA,ITERHA,IPHASE,ITERAD,ITERAC,MCYCAM,IFSBC,JFIRST,
     + NFSD,ICFST,NGRIFL,IGR,NUMCOL,IGRMAX,INTERTU,ICONV3,IGRSLA,J

      INTEGER :: JSTATE(NBGGG,3),IGRID(NBGGG,5),NCHIMT(NBGGG),IT(NBGGG),
     + IL(NBGGG),IK(NBGGG),IDI1(NBGGG),IDI2(NBGGG),IDI3(NBGGG),
     + MGRID(NBGGG),MOV(NBGGG),IDER(NBGGG),INIT(NBGGG),INTERI(NBGGG),
     + INTERJ(NBGGG),INTERK(NBGGG),LAMIN(NBGGG),MIB(NBGGG),MIT(NBGGG),
     + MJB(NBGGG),MJT(NBGGG),MKB(NBGGG),MKT(NBGGG),IROTVE(NBGGG),
     + IPRESCB(NBGGG),LUSGSB(NBGGG),IFLUXB(NBGGG),ICNH(NBGGG),
     + IBTGR(NBGGG),ITPGR(NBGGG),JBTGR(NBGGG),JTPGR(NBGGG),
     + KBTGR(NBGGG),KTPGR(NBGGG),IEVAPB(NBGGG),IDIFFB(NBGGG),
     + MPCASEB(NBGGG),PREVELB(NBGGG),INLRCB(NBGGG),OUTRCB(NBGGG),
     + ICEROB(NBGGG),ITERMPB(NBGGG),IFA(19),IFT(19),IFR(19),NSERIES(19),
     + IUPTEM(NBGGG,3),NBLADE(19),IUPPTB(NBGGG),IFAKEP(19),
     + INTERIA(NBGGG),INTERJA(NBGGG),INTERKA(NBGGG),TURMULB(NBGGG),
     + DRAGMULB(NBGGG),LIFTMULB(NBGGG),DISPERMULB(NBGGG),VMASSB(NBGGG),
     + WFORCEMULB(NBGGG),RINTFB(NBGGG),OSCLS(NBGGG)

      REAL :: CFL,CFLL,DROLIM,TMAX,DT,RMACH,ALPHA,BETA,RE,PR,PRT,FRSTEM,
     + TEMINI,FRSDEN,FRSPRE,FRSVEL,FRSTUR,FRSMUT,TURLIM,TURINI,RMUINI,
     + CMGK,CMGEPS,CSIMPS,CC1,CC2,AREF,CHLREF,GRILEN,XMOM,YMOM,ZMOM,
     + REFPRE,DIFPRE,SOLTEM,SMAX,ARTSSP,RMULTV,ALFAP,
     + DTWMAX,DWMV,GML,FREDIF,DTWMIN,CFLFRE,AGAMMA,ABANK,REFVEL,ROTORW,
     + RJK2,RJK4,WHEIGHT,XBULB,YBULB,ABULB,RETIN,NEGV,CAVLEV,
     + CDIFF,CDIFFT,TROT !,RGASB

      REAL :: OMEGAXP,OMEGAYP,OMEGAZP,TRIMANG,XVAL,YVAL,ZVAL,
     +        XVAL2,YVAL2,ZVAL2

      REAL :: ROTAT(NBGGG),OMEGA(NBGGG),OMEGAX(NBGGG),OMEGAY(NBGGG),
     + OMEGAZ(NBGGG),CENAX(NBGGG),CENAY(NBGGG),CENAZ(NBGGG),AMPL(NBGGG),
     + FRSTEMB(NBGGG),TEMINIB(NBGGG),FRSDENB(NBGGG),FRSPREB(NBGGG),
     + FRSVELB(NBGGG),ARTSSPB(NBGGG),SOLTEMB(NBGGG),FRSTURB(NBGGG),
     + FRSMUTB(NBGGG),TURLIMB(NBGGG),REB(NBGGG),PRB(NBGGG),PRTB(NBGGG),
     + RMACHB(NBGGG),CHLREFB(NBGGG),RMUINIB(NBGGG),AREFB(NBGGG),
     + REFPREB(NBGGG),DIFPREB(NBGGG),FRSALFB(NBGGG,3),TLOLIMB(NBGGG),
     + TUPLIMB(NBGGG),CAVNOB(NBGGG),TAUFB(NBGGG,3),REFVELB(NBGGG),
     + ALFAPB(NBGGG),ALFAUB(NBGGG),CHLCAVB(NBGGG),REFCAVB(NBGGG),
     + ALFASKEWB(NBGGG),BUBBLEB(NBGGG,3),VOIDTTB(NBGGG),
     + GVEXB(NBGGG),GVEYB(NBGGG),GVEZB(NBGGG),FREDIFIB(NBGGG),
     + FREDIFJB(NBGGG),FREDIFKB(NBGGG),RKB(NBGGG,3),FRADENB(NBGGG,3)
      
      REAL :: XCGI(19),YCGI(19),ZCGI(19),
     + RMASS(19),AREFF(19),CHLREFF(19),
     + RCG(19),RIX(19),RIY(19),RIZ(19),RIXY(19),RIXZ(19),RIYZ(19),
     + DAMPN(19),DAMPT(19),TOL(19),H1(19),HMIN(19),TDEL(19),TSTEER(19),
     + SPEED(19),ALPHAF(19),BETAF(19),GAMMAF(19),BANKF(19),VX(19),
     + VY(19),VZ(19),ROTX(19),ROTY(19),ROTZ(19),!PSIR(19),THETAR(19),
     + XCG(19),YCG(19),ZCG(19),PSIRI(19),THETARI(19),PHIRI(19),
     + DRAUGHTI(19),DRAUGHT(19),SINKI(19),SINK(19),
     + TRIMAI(19),TRIMA(19),DAMPC1(19),DAMPC2(19),THETACOL(19),
     + THCYCLON(19),THCYCLAT(19),KDELTA3(19),ZETAH(19),BETAH(19),
     + DOTZE(19),DOTBETA(19),RITH(19),RIZE(19),RIBE(19),RKZE(19),
     + RKBE(19),RD(19),RKD(19),CD(19),DISSZE(19),DISSBE(19),RHINGE(19),
     + BETAHCOL(19),ZETAHCOL(19),SHAFT(19),BECYCLON(19),BECYCLAT(19),
     + RKTH(19),RCZE(19),RCBE(19),RIN(19),ROUT(19),THRUST(19),
     + TORQUE(19),ADV(19),ROTA1(19),ROTB1(19),ROTA1I(19),ROTB1I(19),
     + CONEA(19),CONEAI(19),ACTUA(19),VTIP(19),CDBLADE(19),CFBLADE(19),
     + SIGMA(19),CFTI(19),XCGIS(19),YCGIS(19),ZCGIS(19),
     + XCGA(19),YCGA(19),ZCGA(19),IGRSLAVE(19),RTMSP(19),FDSP(19),
     + UTSP(19),FXSP(19),FXTSP(19),QFACT(19),ETIP(19),CBLADE(19),
     + RGML(19),XFAKEP(19),YFAKEP(19),ZFAKEP(19),ROTB1FAKEP(19) 

      REAL :: PSIR(19),THETAR(19),PHIR(19), PSIM(19)

      CHARACTER(LEN=80) :: CONVCB(NBGGG),FLUXTYB(NBGGG),
     + COMPCORRB(NBGGG),KATOLB(NBGGG),ZEROVB(NBGGG),FLUXCORRB(NBGGG)

      CHARACTER(LEN=10) :: SOLUTION_TYPE(NBGGG)

      NAMELIST /INPUTS/
     + FLUXTY,STARTC,STRESC,FULLNC,SOURC,TIMEC,CONVC,MULPHC,NAME,
     + LEVEL,ITURB,ISTRES,KSCAL,IPRESC,LUSGS,IROTCO,IDRXX,ICMAX,
     + KP,MPRINT,KRRINT,NPHASE,CFL,CFLL,DROLIM,TMAX,DT,RMACH,
     + ALPHA,BETA,RE,PR,PRT,FRSTEM,TEMINI,FRSDEN,FRSPRE,FRSVEL,FRSTUR,
     + FRSMUT,TURLIM,TURINI,RMUINI,IEPSMA,JRDIF,JRDIS,JRPRE,JRIMP,
     + CMGK,CMGEPS,CSIMPS,CC1,CC2,AREF,CHLREF,GRILEN,XMOM,YMOM,ZMOM,
     + REFPRE,DIFPRE,GROUND,GX,GY,GZ,G0,ISTATE,SOLTEM,ARTSSP,
     + MCYCLE,ITERMA,ITERHA,SMAX,ITERAD,ITERAC,MCYCAM,RMULTV,ALFAP,
     + TURCOC,XXTRAC,FRESUC,IFSBC,JFIRST,NFSD,ICFST,DTWMAX,DWMV,GML,
     + FREDIF,DTWMIN,CFLFRE,DISTANC,NUMCOL,AGAMMA,ABANK,TRUE_DISTC,
     + REFVEL,ROTORW,FPRINTC,INTERTU,RJK2,RJK4,TURDESC,INCHIMC,ICONV3,
     + REFLECC,ENTROPY_FIXC,WHEIGHT,XBULB,YBULB,ABULB,TUR_FRS_SOURCEC,
     + TRANSC,TUR_MULPHC,NEGV,NEGVC,STATEC,
     + RGAS, GAMMA, VISU0, EXPSU, TSU0, E0REF, T0REF,
     + GVEX,GVEY,GVEZ,CAVLEV,LONG_STARTC,CDIFF,CDIFFT,AVER_STARTC,TROT,
     + ALTITUDE,START_CAVC,USE_QUATERN,WALLFUNC,AUTOCONC

      
      NAMELIST /BLOCKS/ 
     + IT,IL,IK,IDI1,IDI2,IDI3,MGRID,MOV,IDER,INIT,LAMIN,INTERI,
     + INTERJ,INTERK,MIB,MIT,MJB,MJT,MKB,MKT,JSTATE,IROTVE,NCHIMT,
     + IGRID,ROTAT,OMEGA,OMEGAX,OMEGAY,OMEGAZ,SOLUTION_TYPE,REB,
     + FRSTEMB,TEMINIB,FRSDENB,FRSPREB,FRSVELB,ARTSSPB,SOLTEMB,
     + FRSTURB,FRSMUTB,TURLIMB,PRB,PRTB,RMACHB,CHLREFB,CONVCB,
     + FLUXTYB,IPRESCB,LUSGSB,RMUINIB,AREFB,REFPREB,DIFPREB,FRSALFB,
     + TLOLIMB,TUPLIMB,ICNH,IBTGR,ITPGR,JBTGR,JTPGR,KBTGR,KTPGR,IEVAPB,
     + CAVNOB,TAUFB,IDIFFB,CENAX,CENAY,CENAZ,MPCASEB,REFVELB,ALFAPB,
     + ALFAUB,CHLCAVB,REFCAVB,PREVELB,ALFASKEWB,COMPCORRB,KATOLB,
     + INLRCB,OUTRCB,ICEROB,ITERMPB,BUBBLEB,IUPTEM,VOIDTTB,ZEROVB,
     + FLUXCORRB,GVEXB,GVEYB,GVEZB,IUPPTB,INTERIA,INTERJA,INTERKA,
     + FREDIFIB,FREDIFJB,FREDIFKB,TURMULB,DRAGMULB,LIFTMULB,DISPERMULB,
     + VMASSB,WFORCEMULB,RINTFB,RKB,FRADENB,AMPL,OSCLS

      NAMELIST /FLIGHT/
     + XCGI,YCGI,ZCGI,RMASS,AREFF,CHLREFF,RCG,RIX,RIY,RIZ,RIXY,RIXZ,
     + RIYZ,DAMPN,DAMPT,TOL,H1,HMIN,TDEL,TSTEER,SPEED,ALPHAF,BETAF,
     + GAMMAF,BANKF,VX,VY,VZ,ROTX,ROTY,ROTZ,PSIR,THETAR,PHIR,XCG,YCG,
     + ZCG,PSIRI,THETARI,PHIRI,PSIM,DRAUGHTI,DRAUGHT,SINKI,SINK,
     + TRIMAI,TRIMA,DAMPC1,DAMPC2,THETACOL,
     + THCYCLON,THCYCLAT,KDELTA3,ZETAH,BETAH,
     + DOTZE,DOTBETA,RITH,RIZE,RIBE,RKZE,RKBE,
     + RD,RKD,CD,DISSZE,DISSBE,RHINGE,BETAHCOL,ZETAHCOL,SHAFT,
     + BECYCLON,BECYCLAT,RKTH,RCZE,RCBE,RIN,ROUT,THRUST,TORQUE,
     + ADV,ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI,ACTUA,VTIP,CDBLADE,
     + CFBLADE,SIGMA,IFA,IFT,IFR,NSERIES,CFTI,XCGIS,YCGIS,ZCGIS,
     + XCGA,YCGA,ZCGA,IGRSLAVE,RTMSP,FDSP,UTSP,FXSP,FXTSP,QFACT,
     + ETIP,CBLADE,NBLADE,RGML,XFAKEP,YFAKEP,ZFAKEP,ROTB1FAKEP,IFAKEP

       NAMELIST /FORCE_GROUP_DATA/
     + FORCE_GROUP_INFO

       NAMELIST /FSI/
     +      MPCCIC, FEMXYZ, MODALFSIC, RMESHT,
     +      JOINTS_FILE, MODE_SHAPES_FILE,
     +      ZETA, ACTIVEFREQS
       
       
C ... Set the default values for 'INPUTS'

      NAME   = 'Simulation Case'

      LEVEL  = 1 ! First grid level
      ITURB  = 6 ! SST k-o
      ISTRES = 0 ! No EARSM
      KSCAL  = 0 ! Number of scalars
      FLUXTY = 'INCO'
      INTERTU= 0 ! Turbulence variables either limited or 1st order
      IGRMAX = 0 ! Default number for the current flying object

      MCYCLE = 2 ! Number of algebraic multigrid cycles in solid-block AMG
      ITERMA = 1 ! LGS sweeps on the densest level in solid-block AMG
      ITERHA = 5 ! LGS sweeps on the coarse levels in solid-block AMG

      MCYCAM = 2 ! Number of algebraic multigrid cycles in fluid-block AMG
      ITERAD = 1 ! LGS sweeps on the densest level in fluid-block AMG
      ITERAC = 5 ! LGS sweeps on the coarse levels in fluid-block AMG

      STARTC = 'YES'
      STRESC = 'NO'
      FULLNC = 'NO'
      SOURC  = 'NO'
      TIMEC  = 'NO'
      CONVC  = 'YES'
      MULPHC = 'NO'
      TURCOC = 'YES'
      XXTRAC = 'NO'
      FRESUC = 'NO'
      DISTANC= 'YES'
      TRUE_DISTC = 'YES'
      FPRINTC= 'NO'
      TURDESC= 'NO'
      INCHIMC= 'NO'
      REFLECC= 'NO'
      ENTROPY_FIXC = 'NO'
      TUR_FRS_SOURCEC = 'no' ! Change this to 'yes'
      TRANSC = 'no'
      TUR_MULPHC = 'no'
      NEGVC  = 'no'
      LONG_STARTC = 'yes'
      AVER_STARTC = 'no'
      START_CAVC  = 'no'
      USE_QUATERN = 'no'
      WALLFUNC    = 'no'
      AUTOCONC    = 'YES'
      IPRESC = 0
      LUSGS  = 0
      IROTCO = 0
      IDRXX  = 6

      CFL    = 2.0
      CFLL   = 2.0
      DROLIM = 1.E-8
      TROT   = 0.0
      TMAX   = 0.64E5
      DT     = 1.
      ALFAP  = 0.05
      RJK2   = 0.
      RJK4   = 0.
      CDIFF  = 2. !.5  !2. !Corresponds to 0.5 in the earlier versions
      CDIFFT = 2. !.5  !2.

      ICMAX  = 2
      KP     = 10
      MPRINT = 10000
      KRRINT = 10000
      NUMCOL = 5
      ICONV3 = 0

      RMACH  = 0.
      ALPHA  = 0.
      BETA   = 0.

      RE     = 0.
      PR     = 0.72
      PRT    = 0.9

      ISTATE = 0
      STATEC = 'YES'
      NPHASE = 1

      RGAS   = 287.05287
      GAMMA  = 1.4
      VISU0  = 1.458E-6
      EXPSU  = 1.5
      TSU0   = 110.4
      T0REF  = 0.
      E0REF  = 0.

      FRSTEM = 0.
      TEMINI = 0.
      FRSDEN = 0.
      FRSPRE = 0.
      FRSVEL = 0.
      ARTSSP = 0.
      SOLTEM = 0.
      SMAX   = 20.
      RMULTV = 1.  ! Multiplyer for pseudowater to reduce the Reynolds number

      GVEX   = 0.
      GVEY   = 0.
      GVEZ   = 0.  ! Global grid movement

C ... These values were over-written (Hellsten), if zeroes were given in namelist

      FRSTUR = 0.0    !0.0002 ! The Patria defaults should be given in input
      FRSMUT = 0.01   ! Matti Palin 25.09.2017
      TURLIM = 5000.
      IEPSMA = 0

      TURINI = 0.001 ! Unused (mersu)
      RMUINI = 1.
      CMGK   = 0.1
      CMGEPS = 0.2
      CSIMPS = 5.0

      JRDIF  = 0
      JRDIS  = 0
      JRPRE  = 0
      JRIMP  = 0

      CC1    = 0.1
      CC2    = 0.1 ! 1.0

      AREF   = 1.
      CHLREF = 1.
      GRILEN = 1.

      XMOM   = 0.
      YMOM   = 0.
      ZMOM   = 0.
      REFPRE = 0.
      DIFPRE = -1. ! This sets the value to the free stream pressure
      REFVEL = 0.

      GROUND = 0.
      ALTITUDE = -10000.
      GX     = 0.
      GY     = 0.
      GZ     = 0.
      G0     = 9.80665
      CAVLEV = 0.
      AGAMMA = 0.
      ABANK  = 0.
      ROTORW = 0.

      IFSBC  = 0   ! No free-surface calculation type
      JFIRST = 5  
      NFSD   = 4
      ICFST  = 350
      FREDIF = 1

      DTWMAX = 0.001
      DTWMIN = 50.E-5
      CFLFRE = 0.
      DWMV   = .5E-3 ! Based on a model scale
      GML    = -0.1

      WHEIGHT= 0.
      XBULB  = 0.
      YBULB  = 0.
      ABULB  = 0.
      NEGV   = 0.

C ... Read the namelist file for scalars in ('INPUTS')
*******************************************************************             
      READ(2,NML=INPUTS)
      REWIND 2
*******************************************************************            

      IF(ALTITUDE >= GROUND) THEN
         ISTATE = 1
         GX =  SIN(ALPHA*PII/180.0)  ! Gravity vector
         GY = -COS(ALPHA*PII/180.0)
         GZ = 0.0                    ! Angle of bank not defined ?
         GROUND = -ALTITUDE
      ENDIF

      IF(CAVLEV == 0.) CAVLEV = GROUND           

      CALL FLUXME(FLUXTY,IFLUX,ITURB,PARALLEL) ! change FLUXTY to integer IFLUX

C ... Set logical variables
        
      CALL YESNO(FULLNC,LNPARA(1)) ! PJJM 15.9.1997
      CALL YESNO(STARTC,LNPARA(2))
      CALL YESNO(STRESC,LNPARA(3))
      CALL YESNO(SOURC, LNPARA(4))
      CALL YESNO(CONVC, LNPARA(5))
      CALL YESNO(STATEC,LNPARA(6))
      CALL YESNO(TIMEC, LNPARA(7))
      CALL YESNO(MULPHC,LNPARA(8))
      CALL YESNO(TURCOC,LNPARA(10))
      CALL YESNO(XXTRAC,LNPARA(11))
      CALL YESNO(FRESUC,LNPARA(12))
      CALL YESNO(DISTANC,LNPARA(13))
      CALL YESNO(TRUE_DISTC,LNPARA(14))
      CALL YESNO(FPRINTC,LNPARA(15))
      CALL YESNO(TURDESC,LNPARA(16))
      CALL YESNO(INCHIMC,LNPARA(17))
      CALL YESNO(REFLECC,LNPARA(18))
      CALL YESNO(ENTROPY_FIXC,LNPARA(19))
      CALL YESNO(TUR_FRS_SOURCEC,LNPARA(20))
      CALL YESNO(TRANSC,LNPARA(21))
      CALL YESNO(TUR_MULPHC,LNPARA(22))
      CALL YESNO(NEGVC,LNPARA(23))

      CALL YESNO(LONG_STARTC,LNPARA(25))
      CALL YESNO_STARTL(TRANSC,LNPARA(26))
      CALL YESNO(AVER_STARTC,LNPARA(27))
      CALL YESNO(USE_QUATERN,LNPARA(28))
      CALL YESNO(WALLFUNC,LNPARA(29))
      CALL YESNO(START_CAVC,START_CAVL) ! At the moment no mpi
      CALL YESNO(AUTOCONC,AUTOCONL)

C ... At this moment sourl must be true if RSM is applied
      IF(ITURB >= 21) LNPARA(4) = .TRUE.

C     Turbulence-model version parameter KOVER:

      KOVER = 0         
      IF (ITURB >= 60 .AND. ITURB <= 69) THEN
         KOVER = ITURB - 60
         ITURB = 6
      ELSE IF(ITURB >= 90 .AND. ITURB <= 99) THEN
         KOVER = ITURB - 90
         ITURB = 9
      ELSE IF ((ITURB >= 30).AND.(ITURB < 40)) THEN         
         WRITE(*,*) 'Turbulence models 30...40 not yet available'
         KOVER = ITURB - 30
         ITURB = 3
      END IF

C     DES testing

      IF(LNPARA(16)) THEN ! DES
         IF(ITURB == 6 .OR. ITURB == 9) THEN
         ELSE
         WRITE(*,*)  ' Select SA or SST k-w with DES. Exiting...'
         WRITE(13,*) ' Select SA or SST k-w with DES. Exiting...'
         STOP ' Wrong DES selection'
         ENDIF
      ENDIF ! DES

C     SAS and DES testing

      IF(KOVER == 2 .AND..NOT.LNPARA(16)) THEN
         WRITE(*,*)  ' You have selected DDES option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         WRITE(13,*) ' You have selected DDES option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         STOP 'Wrong k-omega selection'
      ELSEIF(KOVER == 3 .AND..NOT.LNPARA(16)) THEN
         WRITE(*,*)  ' You have selected SAS option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         WRITE(13,*) ' You have selected SAS option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         STOP 'Wrong k-omega selection'
      ELSEIF(KOVER == 7 .AND. .NOT.LNPARA(16)) THEN
         WRITE(*,*)  ' You have selected DDES option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         WRITE(13,*) ' You have selected DDES option, but TRUEDESL is',
     +               ' not activated. Exiting...'
         STOP 'Wrong k-omega selection'
      ENDIF ! SAS and DES

C ... Transfer to the scalar arrays (antakee mersu)

      KSTATE    = ISTATE      ! Global state index
c      ISCALA(1) = IOLD
      ISCALA(2) = LEVEL
      ISCALA(3) = ITURB
      ISCALA(4) = ISTRES
      ISCALA(5) = KSCAL
      ISCALA(6) = IPRESC
      ISCALA(7) = LUSGS
      ISCALA(8) = IROTCO
      ISCALA(9) = IDRXX
      ISCALA(10)= ICMAX
      ISCALA(11)= KP
      ISCALA(12)= MPRINT
      ISCALA(13)= KRRINT

      ISCALA(15)= NPHASE
      ISCALA(16)= IEPSMA
      ISCALA(17)= JRDIF
      ISCALA(18)= JRDIS
      ISCALA(19)= JRPRE
      ISCALA(20)= JRIMP
      ISCALA(22)= KOVER
      ISCALA(23)= MCYCLE
      ISCALA(24)= ITERMA
      ISCALA(25)= ITERHA
      ISCALA(26)= ITERAD
      ISCALA(27)= ITERAC
      ISCALA(28)= MCYCAM

      ISCALA(29)= IFSBC
      ISCALA(30)= JFIRST
      ISCALA(31)= NFSD
      ISCALA(32)= ICFST
      ISCALA(33)= NUMCOL
      ISCALA(34)= INTERTU
      ISCALA(35)= ICONV3

      RSCALA(1) = CFL
      IF(CFLL <= 0.) CFLL = CFL
      RSCALA(2) = CFLL
      RSCALA(3) = DROLIM
      RSCALA(4) = TMAX
      RSCALA(74)= TROT ! Added by M.Palin, April 20th 2016
      RSCALA(5) = DT

      RSCALA(6) = RMACH
      RSCALA(7) = ALPHA
      RSCALA(8) = BETA
      RSCALA(9) = RE
      RSCALA(10)= PR
      RSCALA(11)= PRT
      RSCALA(12)= FRSTEM
      RSCALA(13)= TEMINI
      RSCALA(14)= FRSDEN
      RSCALA(15)= FRSPRE
      RSCALA(16)= FRSVEL

      RSCALA(17)= FRSTUR
      RSCALA(18)= FRSMUT
      RSCALA(19)= TURLIM
      RSCALA(20)= TURINI
      RSCALA(21)= RMUINI

      RSCALA(22)= CMGK
      RSCALA(23)= CMGEPS
      RSCALA(24)= CSIMPS
      RSCALA(25)= CC1
      RSCALA(26)= CC2

      RSCALA(27)= AREF
      RSCALA(28)= CHLREF
      RSCALA(29)= GRILEN

      RSCALA(30)= XMOM
      RSCALA(31)= YMOM
      RSCALA(32)= ZMOM

      RSCALA(33)= REFPRE
      RSCALA(34)= DIFPRE

      IF(SQRT(GX**2+GY**2+GZ**2)  >= 1.E-5) LNPARA(9) = .TRUE. ! GRAVIL
      XAU32 = SQRT((GX+1.E-10)**2 + (GY+1.E-10)**2 + (GZ+1.E-10)**2)

      IF(LNPARA(9)) THEN
         GX      = G0*GX/XAU32
         GY      = G0*GY/XAU32
         GZ      = G0*GZ/XAU32
      ENDIF

      RSCALA(35)= GROUND
      RSCALA(36)= GX
      RSCALA(37)= GY
      RSCALA(38)= GZ
          
      RSCALA(39)= SOLTEM
      RSCALA(40)= SMAX
      RSCALA(41)= ARTSSP
      RSCALA(42)= RMULTV
      RSCALA(43)= ALFAP

      RSCALA(44)= DTWMAX
      RSCALA(45)= DWMV
      RSCALA(46)= GML
      RSCALA(47)= FREDIF
      IF(CFLFRE == 0.) CFLFRE = 6.*CFL
      RSCALA(48)= CFLFRE
      RSCALA(49)= DTWMIN
      RSCALA(50)= AGAMMA
      RSCALA(51)= ABANK
      RSCALA(52)= REFVEL
      RSCALA(53)= ROTORW
      IF(RJK2 > 0.) RSCALA(54)= 1./RJK2
      IF(RJK4 > 0.) RSCALA(55)= 1./RJK4
      RSCALA(56)= WHEIGHT
      RSCALA(57)= XBULB
      RSCALA(58)= YBULB
      RSCALA(59)= ABULB
      RSCALA(60)= NEGV
      RSCALA(61)= RGAS
      RSCALA(62)= GAMMA
      RSCALA(63)= VISU0
      RSCALA(64)= EXPSU
      RSCALA(65)= TSU0
      RSCALA(66)= E0REF
      RSCALA(67)= T0REF
      RSCALA(68)= GVEX
      RSCALA(69)= GVEY
      RSCALA(70)= GVEZ
      RSCALA(71)= CAVLEV
      RSCALA(72)= CDIFF
      RSCALA(73)= CDIFFT
      RSCALA(75)= ALTITUDE
      RSCALA(76)= G0
           
C ... Set the defult values for 'BLOCKS'

      IGRMAX    = 0
      NGRIFL    = 0

      DO N = 1,NBLOCK

         IIMAX = IX(N) - 1
         JJMAX = IY(N) - 1
         KKMAX = IZ(N) - 1

         INPARA(1,N) = IIMAX ! IMAX
         INPARA(2,N) = JJMAX ! JMAX
         INPARA(3,N) = KKMAX ! KMAX
         INTERI(N)   = -3    ! INTER,  third-order upwind
         INTERJ(N)   = -3    
         INTERK(N)   = -3
         LAMIN(N)    = 111   ! LAMIN,  friction in every index direction
         INIT(N)     = 1     ! INIT,   default initialization
         IDER(N)     = 1     ! IDER,   default derivative calculation
         IT(N)       = -1    ! IT,     printing index
         IL(N)       = 3     ! IL,     printing index
         IK(N)       = 1     ! IK,     printing index
         IDI1(N)     = IX(N) ! IDI1,   walls with friction
         IDI2(N)     = IY(N) ! IDI2
         IDI3(N)     = IZ(N) ! IDI3
         MOV(N)      = 0     ! MOV,    default no animation
         MGRID(N)    = -1    ! MGRID,  number of multigrid levels
         MIB(N)      = 1     ! MIB,    lower cell index for multigrid
         MIT(N)      = IIMAX ! MIT,    upper cell index for multigrid
         MJB(N)      = 1     ! MJB
         MJT(N)      = JJMAX ! MJT
         MKB(N)      = 1     ! MKB
         MKT(N)      = KKMAX ! MKT
         JSTATE(N,1) = 0     ! No default selection in the equation of state
         JSTATE(N,2) = 0      
         JSTATE(N,3) = 0
         IUPTEM(N,1) = 2     ! First phase is updated by default
         IUPTEM(N,2) = 1     ! Second phase is saturated
         IUPTEM(N,3) = 0     ! Do nothing
         ROTAT(N)    = 0.    ! ROTAT,  initial rotation angle
         OMEGA(N)    = 0.    ! OMEGA,  rotational speed
         IROTVE(N)   = 0     ! IROTVE, output indicator for the velocities
         IGRID(N,:)  = 0     ! IGRID,  indicator for grid movement
         NCHIMT(N)   = 0     ! NCHIMT, chimera indicator
         OMEGAX(N)   = 1.    ! OMEGAX, direction of the rotational axis
         OMEGAY(N)   = 0.    ! OMEGAY
         OMEGAZ(N)   = 0.    ! OMEGAZ
         CENAX(N)    = 0.    ! CENAX,  a point defining the rotational axis 
         CENAY(N)    = 0.    ! CENAY
         CENAZ(N)    = 0.    ! CENAZ
         AMPL(N)     = 0.    ! Oscillation amplitude
         OSCLS(N)    = 0     ! OSCLS,  oscillation shape, no oscillation
         SOLUTION_TYPE(N) = 'FLUID' ! Single-phase fluid is assumed

         ALFAPB(N)   = ALFAP ! Unused, mersu
         ALFAUB(N)   = 1.    ! Default, should be lowered for poor grids
         ALFASKEWB(N)= .0    ! Default no correction
         COMPCORRB(N)= 'NO'  ! No compressibility correction in PRECOR
         KATOLB(N)   = 'YES' ! Kato-Launder is default
         ZEROVB(N)   = 'NO'  ! Put velocities to zero inside the structures
         INLRCB(N)   = 0     ! Indicator for the inlet mass flux calculation
         OUTRCB(N)   = 0     ! Indicator for the outlet mass flux calculation
         FLUXCORRB(N)= 'NO'  ! Correct the mass fluxes in addition to velocities

         FRSDENB(N)  = FRSDEN
         FRSPREB(N)  = FRSPRE
         FRSTEMB(N)  = FRSTEM
         FRSVELB(N)  = FRSVEL
         ARTSSPB(N)  = ARTSSP
         SOLTEMB(N)  = SOLTEM
         FRSTURB(N)  = FRSTUR
         FRSMUTB(N)  = FRSMUT
         TURLIMB(N)  = TURLIM
         REB(N)      = RE
         PRB(N)      = PR
         PRTB(N)     = PRT
         CHLREFB(N)  = CHLREF
         RMACHB(N)   = RMACH
         CONVCB(N)   = CONVC 

         FLUXTYB(N)  = FLUXTY
         IPRESCB(N)  = IPRESC
         LUSGSB(N)   = LUSGS
         RMUINIB(N)  = RMUINI
         AREFB(N)    = AREF
         REFPREB(N)  = REFPRE
         DIFPREB(N)  = DIFPRE
         TEMINIB(N)  = TEMINI
         REFVELB(N)  = REFVEL  ! FRSVEL
         REFCAVB(N)  = 0.
         FRSALFB(N,1)= 1.
         TAUFB(N,1)  = 0.01
         BUBBLEB(N,1)= 1.E9
         RKB(N,1)    = 1.E14
         FRADENB(N,1) = FRSDEN
         GVEXB(N)    = GVEX
         GVEYB(N)    = GVEY
         GVEZB(N)    = GVEZ         
         DO IPHASE = 2,NPHASES
         FRSALFB(N,IPHASE) = 0.
         TAUFB(N,IPHASE)   = 0.01
         BUBBLEB(N,IPHASE) = 1.E 9
         RKB(N,IPHASE) = 1.E 14
         FRADENB(N,1)  = FRSDEN
         ENDDO

         TLOLIMB(N)  = 0.
         TUPLIMB(N)  = 10000.

         ICNH(N)     = 0
         IBTGR(N)    = 0
         ITPGR(N)    = 0
         JBTGR(N)    = 0
         JTPGR(N)    = 0
         KBTGR(N)    = 0
         KTPGR(N)    = 0

         IEVAPB(N)   = 0      
         CAVNOB(N)   = -1.
         IDIFFB(N)   = 1
         MPCASEB(N)  = 1
         PREVELB(N)  = 2
         ICEROB(N)   = 1
         ITERMPB(N)  = 1
         VOIDTTB(N)  = 1.
         IUPPTB(N)   = 0

         INTERIA(N)  = -1.
         INTERJA(N)  = -1.
         INTERKA(N)  = -1.
         FREDIFIB(N) = 0.
         FREDIFJB(N) = 0.
         FREDIFKB(N) = 0.

         TURMULB(N)   = 0
         DRAGMULB(N)  = 0
         LIFTMULB(N)  = 0
         DISPERMULB(N)= 0
         VMASSB(N)    = 0
         WFORCEMULB(N)= 0
         RINTFB(N)    = 1.
    
      ENDDO
     
      MGRID(1)       = 1     ! MGRID,  a value for grid 1

C ... Read the namelist file for blockwise data in
*******************************************************************
      READ(2,NML=BLOCKS)
      REWIND 2
*******************************************************************
C ... Transfer for the block arrays
      DO N = 1,NBLOCK
       
         IIMAX = IX(N) - 1
         JJMAX = IY(N) - 1
         KKMAX = IZ(N) - 1

         RNORM = SQRT(OMEGAX(N)**2+OMEGAY(N)**2+OMEGAZ(N)**2) + 1.E-20
         OMEGAX(N) = OMEGAX(N)/RNORM
         OMEGAY(N) = OMEGAY(N)/RNORM
         OMEGAZ(N) = OMEGAZ(N)/RNORM
         ROTAT(N)  = ROTAT(N)*DEG2RAD
         IF(ROTAT(N) /= 0. .AND. IGRID(N,1) == 0) THEN
            WRITE(45,*) 'WARNING:'
            WRITE(45,*) 'IGRID(N,1) == 0, the rotation angle was not',
     +      ' changed in block',N
            WRITE(13,*) 'WARNING:'
            WRITE(13,*) 'IGRID(N,1) == 0, the rotation angle was not',
     +      ' changed in block',N
            WRITE(*,*)  'Warning: IGRID(N,1) == 0, the rotation angle',
     +      ' was not changed in block',N
         ENDIF 
         IF(OMEGA(N) /= 0. .AND. IGRID(N,1) == 0) THEN
            WRITE(13,*) 'WARNING:'
            WRITE(13,*) 'IGRID(N,1) == 0, the rotation speed is',
     +      ' nonzero in block',N
            WRITE(*,*)  'Warning: IGRID(N,1) == 0, the rotation speed',
     +      ' is nonzero in block',N
         ENDIF 

         INPARA(1,N) = IIMAX       ! IMAX
         INPARA(2,N) = JJMAX       ! JMAX
         INPARA(3,N) = KKMAX       ! KMAX
         INPARA(4,N) = INTERI(N)   ! INTER,  third-order upwind
         INPARA(5,N) = INTERJ(N) 
         INPARA(6,N) = INTERK(N)
         INPARA(7,N) = LAMIN(N)    ! LAMIN,  friction in every index direction
         INPARA(8,N) = INIT(N)     ! INIT,   default initialization
         INPARA(24,N)= IDER(N)     ! IDER,   default derivative calculation
         INPARA(9,N) = IT(N)       ! IT,     printing index
         INPARA(10,N)= IL(N)       ! IL,     printing index
         INPARA(11,N)= IK(N)       ! IK,     printing index
         INPARA(12,N)= IDI1(N)     ! IDI1,   walls with friction
         INPARA(13,N)= IDI2(N)     ! IDI2
         INPARA(14,N)= IDI3(N)     ! IDI3
         INPARA(15,N)= MOV(N)      ! MOV,    default no animation
         INPARA(16,N)= MGRID(N)    ! MGRID,  number of multigrid levels
         INPARA(17,N)= MIB(N)      ! MIB,    lower cell index for multigrid
         INPARA(18,N)= MIT(N)      ! MIT,    upper cell index for multigrid
         INPARA(19,N)= MJB(N)      ! MJB
         INPARA(20,N)= MJT(N)      ! MJT
         INPARA(21,N)= MKB(N)      ! MKB
         INPARA(22,N)= MKT(N)      ! MKT

         INPARA(25,N)= ISTATE      ! First phase
         INPARA(30,N) = 1          ! Number of phases

         IF(JSTATE(N,1) /= 0)THEN! Blockwise state equation
         INPARA(25,N)= JSTATE(N,1) ! JSTATE, selection in equation of state
         INPARA(30,N)= 1           ! Number of phases
         KSTATE      = JSTATE(N,1) ! The last value is a default if ISTATE=0 
         ENDIF
         IF(JSTATE(N,2) /= 0)THEN!
         INPARA(26,N)= JSTATE(N,2) ! Second phase 
         INPARA(30,N)= 2  
         ENDIF
         IF(JSTATE(N,3) /= 0)THEN!
         INPARA(27,N)= JSTATE(N,3) ! Third phase
         INPARA(30,N)= 3
         ENDIF
         IF(ISTATE /= 0) THEN    ! Default value for free stream
         KSTATE      = ISTATE      ! Free stream initializatin
         ENDIF 
        
         RNPARA(1,N) = ROTAT(N)    ! ROTAT,  initial rotation angle
         RNPARA(2,N) = OMEGA(N)    ! OMEGA,  rotational speed
         INPARA(23,N)= IROTVE(N)   ! IROTVE, print indicator for the velocities
c         INPARA(28,N)= IGRID(N)    ! IGRID,  indicator for grid movement
! INPARA(28,N) is free now, IGRID is changed to positions between 42...48
         INPARA(44,N)= IGRID(N,1)  ! IGRID(N,1), degree of freedom
         INPARA(45,N)= IGRID(N,2)  ! IGRID(N,2), deformation type
         INPARA(46,N)= IGRID(N,3)  ! IGRID(N,3), particle number
         IF(IGRID(N,1) >= 11 .AND. IGRID(N,3) == 0) THEN
            WRITE(13,'(A,I3,2A)') '  N =',N,' , No specification for ',
     &      'IGRID(*,3). Exiting...'
            WRITE(45,*)' No specification for IGRID(*,3). Exiting.'
            WRITE(*,*) ' No specification for IGRID(*,3). See RUN.LOG.'
            STOP
         ENDIF
         INPARA(47,N)= IGRID(N,4)  ! IGRID(N,4), nothing yet, 
         INPARA(48,N)= IGRID(N,5)  ! IGRID(N,5),      waiting for new ideas 
C        the amount of the moving/deformating particles
         IF(IGRID(N,3) >  0) THEN
            IF (IGRID(N,3) >= IGRMAX) THEN
               IGRMAX = IGRID(N,3)
               NGRIFL = IGRID(N,3)
            ENDIF
         ENDIF

         INPARA(29,N)= NCHIMT(N)   ! NCHIMT, chimera indicator
         INPARA(31,N)= IPRESCB(N)  ! IPRESC, pressure correction type
         INPARA(32,N)= LUSGSB(N)   ! LUSGS,  LU-SGS method indicator
            
C        change FLUXTY to integer IFLUX
         CALL FLUXME(FLUXTYB(N),IFLUXB(N),ITURB,PARALLEL) 
         INPARA(33,N)= IFLUXB(N)
            
C        Free-surface global parameters

         INPARA(34,N)= ICNH(N)
         INPARA(35,N)= IBTGR(N)
         INPARA(36,N)= ITPGR(N)
         INPARA(37,N)= JBTGR(N)
         INPARA(38,N)= JTPGR(N)
         INPARA(39,N)= KBTGR(N)
         INPARA(40,N)= KTPGR(N)
               
C        Evaporation type
            
         INPARA(41,N)= IEVAPB(N)

C        Diffusion term type

         INPARA(42,N)= IDIFFB(N)
         INPARA(43,N)= MPCASEB(N)

C        Pressure velocity coupling as PRECOR = 1

         INPARA(49,N)= PREVELB(N)
         INPARA(52,N)= ICEROB(N)   ! PRECOR=2, number pre-vel cycles
         INPARA(53,N)= ITERMPB(N)  ! Increase the number of AMG cycles

C        INL/OUTLET surface interpolation method

         INPARA(50,N)= INLRCB(N)
         INPARA(51,N)= OUTRCB(N)

C        Temperature updating for the phases (here three, mersu)

         INPARA(54,N)= IUPTEM(N,1)
         INPARA(55,N)= IUPTEM(N,2)
         INPARA(56,N)= IUPTEM(N,3)
         INPARA(57,N)= IUPPTB(N)

         INPARA(58,N)= INTERIA(N)
         INPARA(59,N)= INTERJA(N)
         INPARA(60,N)= INTERKA(N)

         INPARA(61,N)= TURMULB(N)
         INPARA(62,N)= DRAGMULB(N)
         INPARA(63,N)= LIFTMULB(N)
         INPARA(64,N)= DISPERMULB(N)
         INPARA(65,N)= VMASSB(N)
         INPARA(66,N)= WFORCEMULB(N)
         INPARA(67,N)= OSCLS(N)

         RNPARA(3,N) = OMEGAX(N)   ! OMEGAX, direction of the rotational axis
         RNPARA(4,N) = OMEGAY(N)   ! OMEGAY
         RNPARA(5,N) = OMEGAZ(N)   ! OMEGAZ
         RNPARA(6,N) = CENAX(N)    ! CENAX,  a point defining the rotational
                                   ! axis 
         RNPARA(7,N) = CENAY(N)    ! CENAY
         RNPARA(8,N) = CENAZ(N)    ! CENAZ
         CNPARA(1,N) = SOLUTION_TYPE(N) ! BLKS%SOLUTION_TYPE, blockwise
         CNPARA(2,N) = CONVCB(N)   ! BLKS%CONVL,    blockwise
         CNPARA(3,N) = COMPCORRB(N)! BLKS%COMPCORR, blockwise
         CNPARA(4,N) = KATOLB(N)   ! BLKS%KATOL,    blockwise
         CNPARA(5,N) = ZEROVB(N)   ! ZEROVB FALSE => ZEROVAR subroutine skipped
         CNPARA(6,N) = FLUXCORRB(N)! FLUXCORR => Perform mass-flux correction
         RNPARA(9,N) = FRSDENB(N)  ! BLKS%FRSDEN,   blockwise
         RNPARA(10,N)= FRSPREB(N)  ! BLKS%FRSPRE,   blockwise
         RNPARA(11,N)= FRSVELB(N)  ! BLKS%FRSVEL,   blockwise
         RNPARA(12,N)= ARTSSPB(N)  ! BLKS%ARTSSP,   blockwise
         RNPARA(13,N)= SOLTEMB(N)  ! BLKS%SOLTEM,   blockwise
         RNPARA(14,N)= FRSTURB(N)  ! BLKS%FRSTUR,   blockwise
         RNPARA(15,N)= FRSMUTB(N)  ! BLKS%FRSMUT,   blockwise
         RNPARA(16,N)= TURLIMB(N)  ! BLKS%TURLIM,   blockwise
         RNPARA(17,N)= REB(N)      ! BLKS%RE,       blockwise
         RNPARA(18,N)= PRB(N)      ! BLKS%PR,       blockwise
         RNPARA(19,N)= PRTB(N)     ! BLKS%PRT,      blockwise
         RNPARA(20,N)= CHLREFB(N)  ! BLKS%CHLREF,   blockwise
         RNPARA(21,N)= RMACHB(N)   ! BLKS%RMACH,    blockwise
         RNPARA(22,N)= RMUINIB(N)  ! BLKS%RMUINI,   blockwise
         RNPARA(23,N)= AREFB(N)    ! BLKS%AREF,     blockwise
         RNPARA(24,N)= REFPREB(N)  ! BLKS%REFPRE,   blockwise
         RNPARA(25,N)= DIFPREB(N)  ! BLKS%DIFPRE,   blockwise
         RNPARA(26,N)= FRSTEMB(N)  ! BLKS%DIFPRE,   blockwise
         RNPARA(27,N)= TEMINIB(N)  ! BLKS%DIFPRE,   blockwise
         RNPARA(37,N)= REFVELB(N)  ! BLKS%REFVEL,   blockwise
         RNPARA(38,N)= ALFAPB(N)   ! BLKS%ALFAP,   blockwise
         RNPARA(39,N)= ALFAUB(N)   ! BLKS%ALFAU,   blockwise
         RNPARA(40,N)= CHLCAVB(N)  ! BLKS%CHLCAV,   blockwise
         RNPARA(41,N)= REFCAVB(N)  ! BLKS%REFCAV,   blockwise
         RNPARA(42,N)= ALFASKEWB(N)! BLKS%ALFASKEW, blockwise
         RNPARA(49,N)= GVEXB(N)    
         RNPARA(50,N)= GVEYB(N)
         RNPARA(51,N)= GVEZB(N)   
         RNPARA(52,N)= FREDIFIB(N)       
         RNPARA(53,N)= FREDIFJB(N)       
         RNPARA(54,N)= FREDIFKB(N)       
         RNPARA(55,N)= RINTFB(N)       
         
         DO IPHASE = 1,NPHASES
         RNPARA(27+IPHASE,N)= FRSALFB(N,IPHASE) ! BLKS%FRSALFA, blockwise
         RNPARA(57+IPHASE,N)= FRADENB(N,IPHASE) ! BLKS%FRADEN, blockwise
         ENDDO
         RNPARA(31,N)= TLOLIMB(N) ! BLKS%TLOLIM, blockwise
         RNPARA(32,N)= TUPLIMB(N) ! BLKS%TUPLIM, blockwise
         RNPARA(33,N)= CAVNOB(N)  ! BLKS%CAVNO,  blockwise
         DO IPHASE = 1,NPHASES
         RNPARA(33+IPHASE,N)= TAUFB(N,IPHASE)   ! BLKS%TAUFB, blockwise
         RNPARA(42+IPHASE,N)= BUBBLEB(N,IPHASE) ! BLKS%N,     blockwise
         RNPARA(55+IPHASE,N)= RKB(N,IPHASE)     ! BLKS%RK,    blockwise
         ENDDO

         RNPARA(60,N) = AMPL(N)       

C        Intermittency variables

         CALL INLETRET(FRSTURB(N),RETIN)
         RNPARA(46,N) = 1.     ! BLKS(N)%FRSG
         RNPARA(47,N) = RETIN  ! BLKS(N)%FRSRET

         RNPARA(48,N) = VOIDTTB(N)

      ENDDO
             
C ... Set the default values for FLIGHT
  
      DO IGR = 1,NGRIFL
         XCGI(IGR)   = 0.
         YCGI(IGR)   = 0.
         ZCGI(IGR)   = 0.
         XCG(IGR)    = 0.
         YCG(IGR)    = 0.
         ZCG(IGR)    = 0.
         RMASS(IGR)  = 1.
         AREFF(IGR)  = 1.
         CHLREFF(IGR)= 1.
         RCG(IGR)    = 1.
         RIX(IGR)    = 1.
         RIY(IGR)    = 1.
         RIZ(IGR)    = 1.
         RIXY(IGR)   = 0.
         RIXZ(IGR)   = 0.
         RIYZ(IGR)   = 0.
         DAMPN(IGR)  = 1.
         DAMPT(IGR)  = 1.
         TOL(IGR)    = 1.E-6
         H1(IGR)     = 1.E-6
         HMIN(IGR)   = 1.E-10
         TDEL(IGR)   = DT
         TSTEER(IGR) = TMAX+DT
         SPEED(IGR)  = FRSVEL
         ALPHAF(IGR) = ALPHA
         BETAF(IGR)  = BETA
         GAMMAF(IGR) = AGAMMA
         BANKF(IGR)  = ABANK
         VX(IGR)     = 0.
         VY(IGR)     = 0.
         VZ(IGR)     = 0.
         ROTX(IGR)   = 0.
         ROTY(IGR)   = 0.
         ROTZ(IGR)   = 0.
         PSIR(IGR)   = 0.
         THETAR(IGR) = 0.
         PHIR(IGR)   = 0. 
         PSIRI(IGR)  = 0.
         THETARI(IGR)= 0.
         PHIRI(IGR)  = 0.
         PSIM(IGR)   = 0.
         DRAUGHTI(IGR) = 0.
         DRAUGHT(IGR)  = 0.
         SINKI(IGR)    = 0.
         SINK(IGR)     = 0.
         TRIMAI(IGR)   = 0.
         TRIMA(IGR)    = 0.
         DAMPC1(IGR)   = -4.
         DAMPC2(IGR)   = -6.
         THETACOL(IGR) = 0.
         THCYCLON(IGR) = 0.
         THCYCLAT(IGR) = 0.
         KDELTA3(IGR)  = 0.
         RITH(IGR)     = 0.
         ZETAH(IGR)    = 0.
         BETAH(IGR)    = 0.
         DOTZE(IGR)    = 0.
         DOTBETA(IGR)  = 0.
         RIZE(IGR)     = 0.
         RIBE(IGR)     = 0.
         RKZE(IGR)     = 0.
         RKBE(IGR)     = 0.
         RD(IGR)       = 0.
         RKD(IGR)      = 0.
         CD(IGR)       = 0.
         DISSZE(IGR)   = 10.
         DISSBE(IGR)   = 10.
         RHINGE(IGR)   = 0.
         BETAHCOL(IGR) = 0.
         ZETAHCOL(IGR) = 0.
         SHAFT(IGR)    = 0.
         BECYCLON(IGR) = 0.
         BECYCLAT(IGR) = 0.
         RKTH(IGR)     = 0.
         RCZE(IGR)     = 0.
         RCBE(IGR)     = 0.
         RIN(IGR)      = 0.
         ROUT(IGR)     = 0.
         THRUST(IGR)   = 0.
         TORQUE(IGR)   = 0.
         ADV(IGR)      = 0.
         ROTA1(IGR)    = 0.
         ROTB1(IGR)    = 0. 
         ROTA1I(IGR)   = 0. 
         ROTB1I(IGR)   = 0. 
         CONEA(IGR)    = 0.
         CONEAI(IGR)   = 0.
         ACTUA(IGR)    = 0.
         VTIP(IGR)     = 218. ! NH-90
         CDBLADE(IGR)  = 0.009
         SIGMA(IGR)    = 0.051
         CFBLADE(IGR)  = 0.003
         IFA(IGR)      = 1
         IFT(IGR)      = 1
         IFR(IGR)      = 0
         NSERIES(IGR)  = 10
         CFTI(IGR)     = 1.
         RTMSP(IGR)    = 0.
         FDSP(IGR)     = 0.
         UTSP(IGR)     = 0.
         FXSP(IGR)     = 0.
         FXTSP(IGR)    = 0.
         QFACT(IGR)    = 1.
         ETIP(IGR)     = 0.
         CBLADE(IGR)   = 0.65
         NBLADE(IGR)   = 4
         RGML(IGR)     = 0.
         XFAKEP(IGR)     = 0.
         YFAKEP(IGR)     = 0.
         ZFAKEP(IGR)     = 0.
         ROTB1FAKEP(IGR) = 0.
         IFAKEP(IGR)     = 0.
         
      ENDDO


C ... Read the namelist file for projectiles in
*******************************************************************
      READ(2,NML=FLIGHT)
      REWIND 2
*******************************************************************

C ... Transfer for the projectile arrays

      DO IGR = 1,NGRIFL
         FNPARA(1,IGR)  = XCGI(IGR)
         FNPARA(2,IGR)  = YCGI(IGR)
         FNPARA(3,IGR)  = ZCGI(IGR)
         FNPARA(4,IGR)  = RMASS(IGR)
         FNPARA(5,IGR)  = AREFF(IGR)
         FNPARA(6,IGR)  = CHLREFF(IGR) ! Mersu
         FNPARA(7,IGR)  = RIX(IGR)
         FNPARA(8,IGR)  = RIY(IGR)
         FNPARA(9,IGR)  = RIZ(IGR)
         FNPARA(10,IGR) = RIXY(IGR)
         FNPARA(11,IGR) = RIXZ(IGR)
         FNPARA(12,IGR) = RIYZ(IGR)
         FNPARA(13,IGR) = DAMPN(IGR)
         FNPARA(14,IGR) = DAMPT(IGR)
         FNPARA(15,IGR) = TOL(IGR)
         FNPARA(16,IGR) = H1(IGR)
         FNPARA(17,IGR) = HMIN(IGR)
         FNPARA(18,IGR) = TDEL(IGR)
         FNPARA(19,IGR) = TSTEER(IGR)         
         FNPARA(20,IGR) = SPEED(IGR)
         FNPARA(21,IGR) = ALPHAF(IGR)
         FNPARA(22,IGR) = BETAF(IGR)
         FNPARA(23,IGR) = GAMMAF(IGR)
         FNPARA(24,IGR) = BANKF(IGR)
         FNPARA(25,IGR) = VX(IGR)
         FNPARA(26,IGR) = VY(IGR)
         FNPARA(27,IGR) = VZ(IGR)
         FNPARA(28,IGR) = ROTX(IGR)
         FNPARA(29,IGR) = ROTY(IGR)
         FNPARA(30,IGR) = ROTZ(IGR)
         FNPARA(31,IGR) = PSIR(IGR)
         FNPARA(32,IGR) = THETAR(IGR)
         FNPARA(33,IGR) = PHIR(IGR)
         FNPARA(34,IGR) = XCG(IGR)
         FNPARA(35,IGR) = YCG(IGR)
         FNPARA(36,IGR) = ZCG(IGR)
         FNPARA(37,IGR) = PSIRI(IGR)
         FNPARA(38,IGR) = THETARI(IGR)
         FNPARA(39,IGR) = PHIRI(IGR)
         FNPARA(40,IGR) = RCG(IGR)
         FNPARA(41,IGR) = PSIM(IGR)
         FNPARA(42,IGR) = DRAUGHTI(IGR)
         FNPARA(43,IGR) = DRAUGHT(IGR)
         FNPARA(44,IGR) = TRIMAI(IGR)
         FNPARA(45,IGR) = TRIMA(IGR)
         FNPARA(46,IGR) = DAMPC1(IGR)
         FNPARA(47,IGR) = DAMPC2(IGR)
         FNPARA(48,IGR) = THETACOL(IGR)    
         FNPARA(49,IGR) = THCYCLON(IGR)    
         FNPARA(50,IGR) = THCYCLAT(IGR)    
         FNPARA(51,IGR) = KDELTA3(IGR)     
         FNPARA(52,IGR) = ZETAH(IGR)        
         FNPARA(53,IGR) = BETAH(IGR)        
         FNPARA(54,IGR) = DOTZE(IGR)       
         FNPARA(55,IGR) = DOTBETA(IGR)     
         FNPARA(56,IGR) = RITH(IGR)
         FNPARA(57,IGR) = RIZE(IGR)        
         FNPARA(58,IGR) = RIBE(IGR)        
         FNPARA(59,IGR) = RKZE(IGR)        
         FNPARA(60,IGR) = RKBE(IGR)        
         FNPARA(61,IGR) = RD(IGR)          
         FNPARA(62,IGR) = RKD(IGR)         
         FNPARA(63,IGR) = CD(IGR)    
         FNPARA(64,IGR) = DISSZE(IGR)
         FNPARA(65,IGR) = DISSBE(IGR)
         FNPARA(66,IGR) = RHINGE(IGR)
         FNPARA(67,IGR) = BETAHCOL(IGR) 
         FNPARA(68,IGR) = ZETAHCOL(IGR) 
         FNPARA(69,IGR) = SHAFT(IGR)
         FNPARA(70,IGR) = BECYCLON(IGR)    
         FNPARA(71,IGR) = BECYCLAT(IGR)  
         FNPARA(72,IGR) = RKTH(IGR)
         FNPARA(73,IGR) = RCZE(IGR)        
         FNPARA(74,IGR) = RCBE(IGR) 
         FNPARA(75,IGR) = RIN(IGR)
         FNPARA(76,IGR) = ROUT(IGR)
         FNPARA(77,IGR) = THRUST(IGR)
         FNPARA(78,IGR) = TORQUE(IGR)
         FNPARA(79,IGR) = ADV(IGR)
         FNPARA(80,IGR) = ROTA1(IGR)
         FNPARA(81,IGR) = ROTB1(IGR)
         FNPARA(82,IGR) = ROTA1I(IGR)
         FNPARA(83,IGR) = ROTB1I(IGR)
         FNPARA(84,IGR) = CONEA(IGR)
         FNPARA(85,IGR) = CONEAI(IGR)
         FNPARA(86,IGR) = ACTUA(IGR)
         FNPARA(87,IGR) = VTIP(IGR)
         FNPARA(88,IGR) = CDBLADE(IGR)
         FNPARA(89,IGR) = CFBLADE(IGR)
         FNPARA(90,IGR) = SIGMA(IGR)
         FNPARA(91,IGR) = CFTI(IGR)
         FNPARA(92,IGR) = RTMSP(IGR)
         FNPARA(93,IGR) = FDSP(IGR)
         FNPARA(94,IGR) = UTSP(IGR)
         FNPARA(95,IGR) = FXSP(IGR)
         FNPARA(96,IGR) = FXTSP(IGR)
         FNPARA(97,IGR) = ETIP(IGR)
         FNPARA(98,IGR) = CBLADE(IGR)
         FNPARA(99,IGR) = QFACT(IGR)
         FNPARA(100,IGR)= RGML(IGR)
         FNPARA(101,IGR)= XFAKEP(IGR)
         FNPARA(102,IGR)= YFAKEP(IGR)
         FNPARA(103,IGR)= ZFAKEP(IGR)
         FNPARA(104,IGR)= ROTB1FAKEP(IGR)
         
         INFPAR(1,IGR)  = IFA(IGR)
         INFPAR(2,IGR)  = IFT(IGR)
         INFPAR(3,IGR)  = IFR(IGR)
         INFPAR(4,IGR)  = NSERIES(IGR)
         INFPAR(5,IGR)  = NBLADE(IGR)       
         INFPAR(6,IGR)  = IFAKEP(IGR)
         
      ENDDO

C ... Set default values to the helicopter case 

      DO N = 1,NBLOCK

         IF (INPARA(44,N) >= 31 .AND. INPARA(44,N) <= 39) THEN ! Rotational 
                                                                   ! movement
               RNPARA(3,N)  = 0.0
               RNPARA(4,N)  = 1.0
               RNPARA(5,N)  = 0.0
               RNPARA(6,N)  = 0.0
               RNPARA(7,N)  = 0.0
               RNPARA(8,N)  = 0.0
               RNPARA(2,N)  = ROTORW

            DO IGR = 1,NGRIFL
               XCG(IGR)       = COS(SHAFT(IGR)*DEG2RAD)*RHINGE(IGR) 
               FNPARA(34,IGR) = XCG(IGR)                           ! new XCG
               PSIM(IGR)      = -(SHAFT(IGR)*DEG2RAD)+PII/2
               FNPARA(41,IGR) = PSIM(IGR)                          ! new PSIM
               ZCG(IGR)       = -COS(PSIM(IGR))*RHINGE(IGR)            
               FNPARA(36,IGR) = ZCG(IGR)                           ! new ZCG
               FNPARA(7,IGR)  = RIBE(IGR)                          ! new RIX
               FNPARA(8,IGR)  = RITH(IGR)                          ! new RIY
               FNPARA(9,IGR)  = RIZE(IGR)                          ! new RIZ
               FNPARA(25,IGR) = -ROTORW*RHINGE(IGR)*COS(PSIM(IGR)) ! new VX
               FNPARA(27,IGR) = -ROTORW*RHINGE(IGR)*SIN(PSIM(IGR)) ! new VZ
               FNPARA(28,IGR) = COS(PSIR(IGR))*DOTBETA(IGR)        ! new ROTX
               FNPARA(29,IGR) = ROTORW - DOTZE(IGR)                ! new ROTY
               FNPARA(30,IGR) = -SIN(PSIR(IGR))*DOTBETA(IGR)       ! new ROTZ
               FNPARA(31,IGR) = PSIM(IGR)*RAD2DEG+ZETAHCOL(IGR)    ! new PSIR
               FNPARA(33,IGR) = -BETAHCOL(IGR)                     ! new PHIR
     &              + BECYCLAT(IGR)*SIN((SHAFT(IGR))*DEG2RAD)
     &              + BECYCLON(IGR)*COS((SHAFT(IGR))*DEG2RAD)
               FNPARA(32,IGR) = THETACOL(IGR)                      ! new THETAR
     &              - THCYCLAT(IGR)*SIN((SHAFT(IGR))*DEG2RAD)
     &              - THCYCLON(IGR)*COS((SHAFT(IGR))*DEG2RAD)
     &              - KDELTA3(IGR)*(-FNPARA(33,IGR))
            ENDDO
         ENDIF ! IF (INPARA(44,N) >= 31 ...
         IF (INPARA(44,N) >= 40 .AND. INPARA(44,N) <= 49) THEN ! Actuator  
            DO IGR = 1,NGRIFL                                      ! disk
               FNPARA(86,IGR) = PII*                               ! ACTUA
     &              (FNPARA(76,IGR)**2-FNPARA(75,IGR)**2)
               IF (INPARA(44,N) == 40) THEN ! Patria's actuator  
                  FNPARA(14,IGR) = 0.16667*RSCALA(14)*FNPARA(69,IGR)**2*
     &                 (FNPARA(76,IGR)**3-FNPARA(75,IGR)**3)*
     &                 INFPAR(5,IGR)*FNPARA(98,IGR)*0.1047
               ENDIF            ! IF (INPARA(44,N) == 40 ...
            ENDDO
         ENDIF ! IF (INPARA(44,N) >= 40 ...

C ... Calculate new values of the ship propel location and the ship propel angle
         IF (INPARA(47,N) /= 0 .AND. 
     &        INPARA(44,N) == 10) THEN ! IGRID(1)=10 and IGRID(4)>0 tested
            IGR     = INPARA(46,N)
            IGRSLA  = INPARA(47,N) ! IGR slave
            XCGIS(IGR)  = XCG(IGR)
            YCGIS(IGR)  = YCG(IGR)
            ZCGIS(IGR)  = ZCG(IGR)
            ZCGD        = ZCGIS(IGR)+
     &           DRAUGHTI(IGRSLA)-DRAUGHT(IGRSLA) ! add draught
            TRIMANG     = (-TRIMA(IGRSLA)+TRIMAI(IGRSLA))*DEG2RAD
            XVAL        = -XCGIS(IGR)+XCG(IGRSLA)   
            YVAL        = -YCGIS(IGR)+YCG(IGRSLA)   
            ZVAL        = -ZCGD+ZCG(IGRSLA)
C            CALL EULFIN(0.0,0.0,0.0, ! new rotation location
C     &           -(XCGIS(IGR)-XCG(IGRSLA)),
C     &           -(ZCGD-ZCG(IGRSLA)),
C     &           -(YCGIS(IGR)-YCG(IGRSLA)),
C     &           XCGA(IGR),YCGA(IGR),ZCGA(IGR),1)
            CALL EULFIN(-TRIMANG,0.0,0.0, ! new rotation location
     &           XVAL,ZVAL,YVAL,XVAL2,YVAL2,ZVAL2,1)
            RNPARA(6,N) = XVAL2+XCG(IGRSLA)
            RNPARA(7,N) = YVAL2+YCG(IGRSLA)
            RNPARA(8,N) = ZVAL2+ZCG(IGRSLA)
            TRIMANG     = (-TRIMA(IGRSLA)+TRIMAI(IGRSLA)
     &           -PSIR(IGR)+PSIRI(IGR))*DEG2RAD
            CALL EULFIN(-TRIMANG,0.0,0.0,  ! new rotation axel
     &           -VX(IGR),-VZ(IGR),-VY(IGR),
     &           OMEGAXP,OMEGAYP,OMEGAZP,1)
            RNPARA(3,N) = OMEGAXP/
     &           SQRT(OMEGAXP**2+OMEGAYP**2+OMEGAZP**2+1.E-20)
            RNPARA(4,N) = OMEGAYP/
     &           SQRT(OMEGAXP**2+OMEGAYP**2+OMEGAZP**2+1.E-20)
            RNPARA(5,N) = OMEGAZP/
     &           SQRT(OMEGAXP**2+OMEGAYP**2+OMEGAZP**2+1.E-20)
            RNPARA(1,N) = PHIR(IGR)*DEG2RAD ! intial rotation angle
            RNPARA(2,N) = SHAFT(IGR) ! omega
         ENDIF ! IF (INPARA(47,N) /= 0 .AND. ...

      ENDDO ! N = 1,NBLOCK

********************************************************************

      IFLUX = IFLUXB(1)  ! Take the default from the first block

C ... No read in additional input for the 2-dof model

ccc      INQUIRE(FILE='BLADE_INPUT',EXIST=THERE)
            
        
C ... Set the defult values for 'FORCE_GROUP_DATA'

      XMOMR = XMOM
      YMOMR = YMOM
      ZMOMR = ZMOM

c      XMOMA = -5.0E20
c      XMOMA = -5.0E10  ! Matti Palin 25.09.2017
      XMOMA = -5.0D10   ! Esa Oct 25, 2018
c      YMOMA = -5.0E20
c      YMOMA = -5.0E10  ! Matti Palin 25.09.2017
      YMOMA = -5.0D10   ! Esa Oct 25, 2018 
c      ZMOMA = -5.0E20
c      ZMOMA = -5.0E10  ! Matti Palin 25.09.2017
      ZMOMA = -5.0D10   ! Esa Oct 25, 2018 

      FORCE_GROUP_SHORT_NAME = ' '
      FORCE_GROUP_FULL_NAME = ' '

C ... Read the namelist file for force goup data

C ... Read the namelist file for force groups
*******************************************************************
      READ(2,NML=FORCE_GROUP_DATA)
      REWIND 2
*******************************************************************

      DO J=1,52

         IF(J <= 26) FORCE_GROUP_SHORT_NAME(J:J) = CHAR(J+64)  ! Uppercase
         IF(J >  26) FORCE_GROUP_SHORT_NAME(J:J) = CHAR(J+70)  ! Lowercase

         IGR = ICHAR(FORCE_GROUP_INFO(J)%FG_SHORT_NAME)

         IF(IGR >= 65 .AND. IGR <= 90) THEN       ! Uppercase character
            IGR = IGR - 64
         ELSEIF(IGR >= 97 .AND. IGR <= 122) THEN  ! Lowercase character
            IGR = IGR - 70
         ENDIF

         IF(IGR > 0) THEN
            XMOMR(IGR) = FORCE_GROUP_INFO(J)%XMOM_REF
            YMOMR(IGR) = FORCE_GROUP_INFO(J)%YMOM_REF
            ZMOMR(IGR) = FORCE_GROUP_INFO(J)%ZMOM_REF
            XMOMA(IGR) = FORCE_GROUP_INFO(J)%XMOM_AXS
            YMOMA(IGR) = FORCE_GROUP_INFO(J)%YMOM_AXS
            ZMOMA(IGR) = FORCE_GROUP_INFO(J)%ZMOM_AXS
            FORCE_GROUP_FULL_NAME((IGR-1)*24+1:IGR*24)=
     &      FORCE_GROUP_INFO(J)%FG_FULL_NAME(1:24)
         ENDIF

      ENDDO

C ... Set the default values and initializations for 'FSI'

      MPCCIC      = 'no'
      MODALFSIC   = 'no'
      FEMXYZ      = 0
      RMESHT      = 1
      ACTIVEFREQS = 1
      ZETA        = 0.0
      FM          = 0.0

C ... Read the namelist file for fsi data
*******************************************************************
      READ(2,NML=FSI)
      REWIND 2
*******************************************************************

      CALL YESNO(MPCCIC,LNPARA(24))
      CALL YESNO(MODALFSIC,LNPARA(30))


      
      RETURN
      END SUBROUTINE SET_RINP

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE BACK_RINP(ISCALA,INSCAL,RSCALA,IRSCAL,LNPARA,LRNPAR)

      USE NS3CO
      USE INTEGERS, ONLY : ITERAD, ITERAC, NUMCOL, MCYCAM
      USE BLADE_VARIABLES, ONLY : DTB


      IMPLICIT NONE

      INTEGER :: INSCAL, IRSCAL, LRNPAR
      INTEGER :: ISCALA(INSCAL)
      REAL    :: RSCALA(IRSCAL)
      REAL    :: GML
      LOGICAL :: LNPARA(LRNPAR)

C ... Back substitution of the namelist input data (antakee mersu)

      FULLNS         = LNPARA(1) 
      STARTL         = LNPARA(2)

      STRESL         = LNPARA(3)
      SOURL          = LNPARA(4)
      CONVL          = LNPARA(5)
      STATEL         = LNPARA(6)
      TIMEL          = LNPARA(7)
      MULPHL         = LNPARA(8)
      GRAVIL         = LNPARA(9)
      TURCOR         = LNPARA(10)
      XXTRAL         = LNPARA(11)
      FRESUL         = LNPARA(12)
      LDISTL         = LNPARA(13)
      TRUE_DISTL     = LNPARA(14)
      FPRINTL        = LNPARA(15)
      TURDESL        = LNPARA(16)
      INCHIML        = LNPARA(17)
      REFLECL        = LNPARA(18)
      ENTROPY_FIX    = LNPARA(19)
      TUR_FRS_SOURCE = LNPARA(20)
      TRANSL         = LNPARA(21) 
      TUR_MULPHL     = LNPARA(22)
      NEGVL          = LNPARA(23)
      MPCCIL         = LNPARA(24)
      LONG_STARTL    = LNPARA(25)
      LONG_TRANSL    = LNPARA(26)
      AVER_STARTL    = LNPARA(27)
      USE_QUATERNIONS_LOGICAL    = LNPARA(28)
      WALLFUNL       = LNPARA(29)
      MODALFSIL      = LNPARA(30)

c     IOLD     = ISCALA(1) ! This has been already read in from WORKS
      LEVEL    = ISCALA(2)
      ITURB    = ISCALA(3)
      ISTRES   = ISCALA(4)
      KSCAL    = ISCALA(5)
      IPRESC   = ISCALA(6)
      LUSGS    = ISCALA(7)
      IROTCO   = ISCALA(8)
      IDRXX    = ISCALA(9)
      ICMAX    = ISCALA(10)
      KP       = ISCALA(11)
      MPRINT   = ISCALA(12)
      KRRINT   = ISCALA(13)

c      NPHASE   = ISCALA(15)
      IEPSMA   = ISCALA(16)
      JRDIF    = ISCALA(17)
      JRDIS    = ISCALA(18)
      JRPRE    = ISCALA(19)
      JRIMP    = ISCALA(20)
      KOVER    = ISCALA(22)
      MCYCLE   = ISCALA(23)
      ITERMA   = ISCALA(24)
      ITERHA   = ISCALA(25)
      ITERAD   = ISCALA(26)
      ITERAC   = ISCALA(27)
      MCYCAM   = ISCALA(28)
      IFSBC    = ISCALA(29)
      JFIRST   = ISCALA(30)
      NFSD     = ISCALA(31)
      ICFST    = ISCALA(32)
      NUMCOL   = ISCALA(33)
      INTERTU  = ISCALA(34)
      ICONV3   = ISCALA(35)

      CFL      = RSCALA(1)
      CFLL     = RSCALA(2)
      DROLIM   = RSCALA(3)
      TMAX     = RSCALA(4)
      TROT     = RSCALA(74)   ! Added by M.Palin, April 20th 2016.
      DTB      = RSCALA(5)

      RMACH    = RSCALA(6)
      ALPHA    = RSCALA(7)
      BETA     = RSCALA(8)

      RE       = RSCALA(9)
      PR       = RSCALA(10)
      PRT      = RSCALA(11)
      FRSTEM   = RSCALA(12)
      TEMINI   = RSCALA(13)
      FRSDEN   = RSCALA(14)
      FRSPRE   = RSCALA(15)
      FRSVEL   = RSCALA(16)

      RKLIM    = RSCALA(17)
      EPSLIM   = RSCALA(18)
      TURLIM   = RSCALA(19)
      TURBLE   = RSCALA(20)
      RMUINI   = RSCALA(21)
      CMGK     = RSCALA(22)
      CMGEPS   = RSCALA(23)
      CSIMPS   = RSCALA(24)
      CC1      = RSCALA(25)
      CC2      = RSCALA(26)

      AREF     = RSCALA(27)
      CHLREF   = RSCALA(28)
      GRILEN   = RSCALA(29)

      XMOM     = RSCALA(30)
      YMOM     = RSCALA(31)
      ZMOM     = RSCALA(32)

      REFPRE   = RSCALA(33)
      DIFPRE   = RSCALA(34)

      GROUND   = RSCALA(35)
      GX       = RSCALA(36)
      GY       = RSCALA(37)
      GZ       = RSCALA(38)
      SOLTEM   = RSCALA(39)
      SMAX     = RSCALA(40)
      ARTSSP   = RSCALA(41)
      RMULTV   = RSCALA(42)
      ALFAP    = RSCALA(43)
      DTWMAX   = RSCALA(44)
      DWMV     = RSCALA(45)
      GML      = RSCALA(46)
      FREDIF   = RSCALA(47)
      CFLFRE   = RSCALA(48)
      DTWMIN   = RSCALA(49) 
      AGAMMA   = RSCALA(50)
      ABANK    = RSCALA(51) 
      REFVEL   = RSCALA(52)
      ROTORW   = RSCALA(53)
      RJK2     = RSCALA(54)
      RJK4     = RSCALA(55)
      WHEIGHT  = RSCALA(56)
      XBULB    = RSCALA(57)
      YBULB    = RSCALA(58)
      ABULB    = RSCALA(59)
      NEGV     = RSCALA(60)
      RGAS     = RSCALA(61)
      GAMMA    = RSCALA(62)
      VISU0    = RSCALA(63)
      EXPSU    = RSCALA(64)
      TSU0     = RSCALA(65)
      E0REF    = RSCALA(66)
      T0REF    = RSCALA(67)
      GVEX     = RSCALA(68)
      GVEY     = RSCALA(69)
      GVEZ     = RSCALA(70)
      CAVLEV   = RSCALA(71)
      CDIFF    = RSCALA(72)
      CDIFFT   = RSCALA(73)
      ALTITUDE = RSCALA(75)
      G0       = RSCALA(76)

      RETURN
      END SUBROUTINE BACK_RINP
 
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
    
      SUBROUTINE BACK_FLIGHT(FNPARA,CFPARA,IFNPAR,INPARA,INPAR,
     &                       INFPAR,INFNPAR)

      USE NS3CO, ONLY : NBLOCK,GX,GY,GZ,G0,ALPHA,AGAMMA,ABANK,GRAVIL,
     &                  USE_QUATERNIONS_LOGICAL
      USE FLIGHT,    ONLY : FLYOBJ, MVSHIP, SHIP, SHIPWAVE, ACTDISK,
     &                      HROTOR, HMVBLADE, IGRSLAVE, SHIPPROP,
     &                      XCGI, YCGI, ZCGI, PSIR, THETAR, PHIR,
     &                      PSIRI, THETARI, PHIRI, PSIM,
     &                      DRAUGHTI, DRAUGHT, SINKI, SINK, 
     &                      TRIMAI, TRIMA, DAMPC1, DAMPC2, RIN, ROUT,
     &                      THRUST, TORQUE, ADV, ROTA1, ROTB1, 
     &                      ROTA1I, ROTB1I, XCG, YCG, ZCG, 
     &                      CONEA, CONEAI, IFA, IFT, IFR, NSERIES,
     &                      ACTUA, VTIP, CDBLADE, CFBLADE, SIGMA, CFTI,
     &                      TRMODE, XCGIS, YCGIS, ZCGIS, ROTA1S, ROTB1S,
     &                      NGRIFL, OSKU, RTMSP, FDSP, UTSP, FXSP, 
     &                      FXTSP,QFACT,NBLADE,RGML,XFAKEP,YFAKEP,
     &                      ZFAKEP,ROTB1FAKEP,IFAKEP,
     &                      particleQuaternion0, particleQuaternion1,
     &			    particleQuaternion2, particleQuaternion3
      USE BLADE_VARIABLES

      USE CONSTANTS, ONLY : PII
      
      IMPLICIT NONE

      INTEGER :: IFNPAR, IGR, INPARA(INPAR,*), INPAR, N,
     &           INFPAR(INFNPAR,*), INFNPAR, IGRSLA

      REAL :: FNPARA(IFNPAR,*),XAU32
      REAL,DIMENSION(19) :: ROTB1A,ZCGD,XCGA,YCGA,ZCGA

      CHARACTER(LEN=1) :: CFPARA(*)

      FLYOBJ  = .FALSE.
      HROTOR  = .FALSE.
      HMVBLADE= .FALSE.
      SHIP    = .FALSE.
      MVSHIP  = .FALSE.
      SHIPWAVE= .FALSE.
      ACTDISK = .FALSE.
      SHIPPROP= .FALSE.

      DO IGR = 1,NGRIFL                  
         XCGI(IGR) = FNPARA(1,IGR)  
         YCGI(IGR) = FNPARA(2,IGR)  
         ZCGI(IGR) = FNPARA(3,IGR) 
         XCG(IGR)  = FNPARA(34,IGR)  
         YCG(IGR)  = FNPARA(35,IGR)  
         ZCG(IGR)  = FNPARA(36,IGR)
         OSKU(IGR)%RMASS  = FNPARA(4,IGR)
         OSKU(IGR)%AREF   = FNPARA(5,IGR) 
         OSKU(IGR)%CHLREF = FNPARA(6,IGR)
         OSKU(IGR)%IX     = FNPARA(7,IGR) 
         OSKU(IGR)%IY     = FNPARA(8,IGR)  
         OSKU(IGR)%IZ     = FNPARA(9,IGR)  
         OSKU(IGR)%IXY    = FNPARA(10,IGR) 
         OSKU(IGR)%IXZ    = FNPARA(11,IGR) 
         OSKU(IGR)%IYZ    = FNPARA(12,IGR) 
         OSKU(IGR)%DAMPN  = FNPARA(13,IGR)
         OSKU(IGR)%DAMPT  = FNPARA(14,IGR) 
         OSKU(IGR)%TOL    = FNPARA(15,IGR) 
         OSKU(IGR)%H1     = FNPARA(16,IGR) 
         OSKU(IGR)%HMIN   = FNPARA(17,IGR) 
         OSKU(IGR)%TDEL   = FNPARA(18,IGR) 
         OSKU(IGR)%TSTEER = FNPARA(19,IGR) 
         OSKU(IGR)%SPEED  = FNPARA(20,IGR) 
         OSKU(IGR)%ALPHA  = FNPARA(21,IGR)
         OSKU(IGR)%BETA   = FNPARA(22,IGR) 
         OSKU(IGR)%GAMMA  = FNPARA(23,IGR) 
         OSKU(IGR)%BANK   = FNPARA(24,IGR)
         OSKU(IGR)%VX     = FNPARA(25,IGR) 
         OSKU(IGR)%VY     = FNPARA(26,IGR) 
         OSKU(IGR)%VZ     = FNPARA(27,IGR)
         OSKU(IGR)%ROTX   = FNPARA(28,IGR) 
         OSKU(IGR)%ROTY   = FNPARA(29,IGR) 
         OSKU(IGR)%ROTZ   = FNPARA(30,IGR)

         PSIR(IGR)        = FNPARA(31,IGR)*DEG2RAD 
         THETAR(IGR)      = FNPARA(32,IGR)*DEG2RAD  
         PHIR(IGR)        = FNPARA(33,IGR)*DEG2RAD 
c         OSKU(IGR)%CHAR   = CFPARA(IGR)
         PSIRI(IGR)       = FNPARA(37,IGR)*DEG2RAD  
         THETARI(IGR)     = FNPARA(38,IGR)*DEG2RAD  
         PHIRI(IGR)       = FNPARA(39,IGR)*DEG2RAD 
         OSKU(IGR)%RCG    = FNPARA(40,IGR)
         PSIM(IGR)        = FNPARA(41,IGR)
         DRAUGHTI(IGR)    = FNPARA(42,IGR)
         DRAUGHT(IGR)     = FNPARA(43,IGR)
         SINKI(IGR)       = FNPARA(43,IGR)
         SINK(IGR)        = FNPARA(43,IGR)
         TRIMAI(IGR)      = FNPARA(44,IGR)*DEG2RAD
         TRIMA(IGR)       = FNPARA(45,IGR)*DEG2RAD 
         DAMPC1(IGR)      = FNPARA(46,IGR)
         DAMPC2(IGR)      = FNPARA(47,IGR)
         THETACOL(IGR)    = FNPARA(48,IGR)*DEG2RAD
         THCYCLON(IGR)    = FNPARA(49,IGR)*DEG2RAD
         THCYCLAT(IGR)    = FNPARA(50,IGR)*DEG2RAD
         KDELTA3(IGR)     = FNPARA(51,IGR)
         ZETAH(IGR)       = FNPARA(52,IGR)*DEG2RAD
         BETAH(IGR)       = FNPARA(53,IGR)*DEG2RAD
         DOTZE(IGR)       = FNPARA(54,IGR)*DEG2RAD
         DOTBETA(IGR)     = FNPARA(55,IGR)*DEG2RAD
         RITH(IGR)        = FNPARA(56,IGR)
         RIZE(IGR)        = FNPARA(57,IGR)
         RIBE(IGR)        = FNPARA(58,IGR)
         RKZE(IGR)        = FNPARA(59,IGR)
         RKBE(IGR)        = FNPARA(60,IGR)
         RD(IGR)          = FNPARA(61,IGR)
         RKD(IGR)         = FNPARA(62,IGR)
         CD(IGR)          = FNPARA(63,IGR)
         DISSZE(IGR)      = FNPARA(64,IGR)
         DISSBE(IGR)      = FNPARA(65,IGR)
         RHINGE(IGR)      = FNPARA(66,IGR)
         BETAHCOL(IGR)    = FNPARA(67,IGR)*DEG2RAD
         ZETAHCOL(IGR)    = FNPARA(68,IGR)*DEG2RAD
         SHAFT(IGR)       = FNPARA(69,IGR)*DEG2RAD
         BECYCLON(IGR)    = FNPARA(70,IGR)*DEG2RAD
         BECYCLAT(IGR)    = FNPARA(71,IGR)*DEG2RAD
         RKTH(IGR)        = FNPARA(72,IGR)
         RCZE(IGR)        = FNPARA(73,IGR)
         RCBE(IGR)        = FNPARA(74,IGR)
         RIN(IGR)         = FNPARA(75,IGR) 
         ROUT(IGR)        = FNPARA(76,IGR) 
         THRUST(IGR)      = FNPARA(77,IGR)
         TORQUE(IGR)      = FNPARA(78,IGR) 
         ADV(IGR)         = FNPARA(79,IGR) 
         ROTA1(IGR)       = FNPARA(80,IGR)*DEG2RAD 
         ROTB1(IGR)       = FNPARA(81,IGR)*DEG2RAD 
         ROTA1I(IGR)      = FNPARA(82,IGR)*DEG2RAD 
         ROTB1I(IGR)      = FNPARA(83,IGR)*DEG2RAD 
         CONEA(IGR)       = FNPARA(84,IGR)*DEG2RAD 
         CONEAI(IGR)      = FNPARA(85,IGR)*DEG2RAD
         ACTUA(IGR)       = FNPARA(86,IGR)
         VTIP(IGR)        = FNPARA(87,IGR)
         CDBLADE(IGR)     = FNPARA(88,IGR)
         CFBLADE(IGR)     = FNPARA(89,IGR)
         SIGMA(IGR)       = FNPARA(90,IGR)
         CFTI(IGR)        = FNPARA(91,IGR)
         RTMSP(IGR)       = FNPARA(92,IGR)
         FDSP(IGR)        = FNPARA(93,IGR)
         UTSP(IGR)        = FNPARA(94,IGR)
         FXSP(IGR)        = FNPARA(95,IGR)
         FXTSP(IGR)       = FNPARA(96,IGR)
         ETIP(IGR)        = FNPARA(97,IGR)*DEG2RAD
         CBLADE(IGR)      = FNPARA(98,IGR)
         QFACT(IGR)       = FNPARA(99,IGR)
         RGML(IGR)        = FNPARA(100,IGR)
         XFAKEP(IGR)      = FNPARA(101,IGR)
         YFAKEP(IGR)      = FNPARA(102,IGR)
         ZFAKEP(IGR)      = FNPARA(103,IGR)
         ROTB1FAKEP(IGR)  = FNPARA(104,IGR)
         
c       ENDDO

C ...  Calculate the initial quaternion values for a flying object, by converting them from the Euler angles.
C ...  Attention! It is assumed that ICASE=1, i.e. that xyz rotation order is assumed.
      IF (USE_QUATERNIONS_LOGICAL) THEN
         CALL EULQUA(PSIR(IGR),THETAR(IGR),PHIR(IGR),
     &   particleQuaternion0,particleQuaternion1,
     &   particleQuaternion2,particleQuaternion3,
     &   1)
      END IF


C ...  Calculate a new gravity vector
       DO N = 1,NBLOCK
          IF (INPARA(44,N) >= 11 .AND. ! Particle is flying 
     &         INPARA(44,N) <=  19) THEN        
             GX = (SIN((ALPHA*PII/180.0)+
     &            (AGAMMA*PII/180.0)))
             GY = -(COS(ALPHA*PII/180.0+
     &            AGAMMA*PII/180.0)*
     &            COS(ABANK*PII/180.0))
             GZ = -(COS(ALPHA*PII/180.0+
     &            AGAMMA*PII/180.0)
     &            * SIN(ABANK*PII/180.0))
             FLYOBJ = .TRUE.        
             OSKU(INPARA(46,N))%CHAR  = 'F'
          ELSEIF(INPARA(44,N) == 21 .OR. ! Ship case 
     &           INPARA(44,N) >= 23 .AND. INPARA(44,N) <= 29) THEN    
             TRMODE(INPARA(46,N)) = INPARA(44,N)
             OSKU(IGR)%SPEED = 0.0
             PSIRI(IGR)      = FNPARA(44,IGR)*DEG2RAD  ! move trim data to 
             PSIR(IGR)       = FNPARA(45,IGR)*DEG2RAD  ! psir position
             XCGIS(IGR)      = FNPARA(34,IGR) ! cg location at start 
             YCGIS(IGR)      = FNPARA(35,IGR) ! point
             ZCGIS(IGR)      = FNPARA(36,IGR)
             SINK(IGR)       = FNPARA(43,IGR)
             SHIP            = .TRUE.
             IFAKEP(IGR)     = INFPAR(6,IGR)
             OSKU(INPARA(46,N))%CHAR  = 'S'
             IF(INPARA(45,N) == 5 .OR.  INPARA(45,N) == 9) THEN 
                SHIPWAVE     = .TRUE.
             ENDIF
          ELSEIF(INPARA(44,N) == 22) THEN
             TRMODE(INPARA(46,N)) = INPARA(44,N) 
             OSKU(IGR)%SPEED = 0.0
             PSIRI(IGR)      = FNPARA(44,IGR)*DEG2RAD  ! move trim data to 
             PSIR(IGR)       = FNPARA(45,IGR)*DEG2RAD  ! psir position
             XCGIS(IGR)      = FNPARA(34,IGR) ! cg location at start 
             YCGIS(IGR)      = FNPARA(35,IGR) ! point
             ZCGIS(IGR)      = FNPARA(36,IGR)
             SINK(IGR)       = FNPARA(43,IGR)
             SHIP            = .TRUE.
             MVSHIP          = .TRUE.
             OSKU(INPARA(46,N))%CHAR  = 'S'
             IFAKEP(IGR)     = INFPAR(6,IGR)
             IF(INPARA(45,N) == 5 .OR.  INPARA(45,N) == 9) THEN 
                SHIPWAVE     = .TRUE.
             ENDIF
          ELSEIF (INPARA(44,N) >= 31 .AND. ! roottorilasku 
     &            INPARA(44,N) <=  39) THEN!vai maaritteleeko helikopterimiehet 
             GY          = -1           !g-vekin itse voiko tama aih sotkua?
             TRMODE(INPARA(46,N)) = INPARA(44,N) ! tata tarvitaan myohemmin
             HROTOR      = .TRUE.
             IF(INPARA(44,N) == 38) THEN 
                HMVBLADE     = .TRUE. 
             ENDIF
             OSKU(INPARA(46,N))%CHAR  = 'R'
          ELSEIF (INPARA(44,N) >= 40 .AND. 
     &            INPARA(44,N) <= 49) THEN ! Actuator disk 
             ACTDISK = .TRUE. 
             TRMODE(INPARA(46,N)) = INPARA(44,N) ! tata tarvitaan myohemmin 
             IFA(IGR)    = INFPAR(1,IGR)
             IFT(IGR)    = INFPAR(2,IGR)
             IFR(IGR)    = INFPAR(3,IGR) 
             NSERIES(IGR)= INFPAR(4,IGR) 
             OSKU(INPARA(46,N))%CHAR  = 'A'
             NBLADE(IGR) = INFPAR(5,IGR)
          ENDIF                 ! IF (INPARA(44,NBGGG) >= 11 
       ENDDO                    ! DO N = 1,NBLOCK

C ... Update the gravity vector
       IF(SQRT(GX**2+GY**2+GZ**2)  >= 1.E-5) GRAVIL = .TRUE.
       XAU32 =SQRT((GX+1.E-10)**2 +(GY+1.E-10)**2 +(GZ+1.E-10)**2)
       IF(GRAVIL) THEN
          GX = G0*GX/XAU32
          GY = G0*GY/XAU32
          GZ = G0*GZ/XAU32
       ENDIF                    ! IF (GRAVIL)

      ENDDO

C ... Calculate new values of the ACT location and the ACT angle
c      DO IGR = 1,NGRIFL
      DO N = 1,NBLOCK
         IF (INPARA(47,N) /= 0 .AND. 
     &        INPARA(44,N) == 49) THEN ! IGRID(1)=49 and IGRID(4)>0 tested
            IGR     = INPARA(46,N)
            IGRSLA  = INPARA(47,N) ! IGR slave
            OSKU(INPARA(46,N))%CHAR  = 'B' ! output data in FORCES needs this
            XCGIS(IGR)  = FNPARA(34,IGR) ! ACT location at start 
            YCGIS(IGR)  = FNPARA(35,IGR) ! point
            ZCGIS(IGR)  = FNPARA(36,IGR)
            ROTA1S(IGR) = FNPARA(80,IGR)*DEG2RAD 
            ROTB1S(IGR) = FNPARA(81,IGR)*DEG2RAD 
            ZCGD(IGR)   = ZCGIS(IGR)+
     &           DRAUGHTI(IGRSLA)-DRAUGHT(IGRSLA) ! add draught
            CALL EULFIN(PSIR(IGRSLA),0.0,0.0, ! rotation
     &           -(XCGIS(IGR)-XCG(IGRSLA)),-(ZCGD(IGR)-ZCG(IGRSLA)),
     &           -(YCGIS(IGR)-YCG(IGRSLA)),XCGA(IGR),YCGA(IGR),
     &           ZCGA(IGR),1) 
            XCG(IGR)      = XCGA(IGR)+XCG(IGRSLA) ! new locations
            YCG(IGR)      = YCGA(IGR)+YCG(IGRSLA)
            ZCG(IGR)      = ZCGA(IGR)+ZCG(IGRSLA)
            ROTB1A(IGR)   = ROTB1S(IGR)-PSIR(IGRSLA)
            ROTB1(IGR)    = ROTB1A(IGR)
            IGRSLAVE(IGR) = IGRSLA

C ... Calculate new values of the ship propel location and the ship propel angle
         ELSEIF (INPARA(47,N) /= 0 .AND. 
     &        INPARA(44,N) == 10) THEN ! IGRID(1)=10 and IGRID(4)>0 tested
            SHIPPROP = .TRUE. 
            IGR      = INPARA(46,N)
            IGRSLA   = INPARA(47,N) ! IGR slave
            TRMODE(INPARA(46,N)) = INPARA(44,N)
            OSKU(INPARA(46,N))%CHAR = 'P' ! output data in FORCES needs this
            XCGIS(IGR)  = FNPARA(34,IGR)
            YCGIS(IGR)  = FNPARA(35,IGR)
            ZCGIS(IGR)  = FNPARA(36,IGR)
            ZCGD(IGR)   = ZCGIS(IGR)+
     &           DRAUGHTI(IGRSLA)-DRAUGHT(IGRSLA) ! add draught
            CALL EULFIN(PSIR(IGRSLA),0.0,0.0, ! rotation
     &           -(XCGIS(IGR)-XCG(IGRSLA)),-(ZCGD(IGR)-ZCG(IGRSLA)),
     &           -(YCGIS(IGR)-YCG(IGRSLA)),XCGA(IGR),YCGA(IGR),
     &           ZCGA(IGR),1) 
            XCG(IGR)    = XCGA(IGR)+XCG(IGRSLA) ! new locations
            YCG(IGR)    = YCGA(IGR)+YCG(IGRSLA)
            ZCG(IGR)    = ZCGA(IGR)+ZCG(IGRSLA)
            ROTB1A(IGR) = PSIR(IGRSLA)+FNPARA(31,IGR)*DEG2RAD!PSIR(IGR)
            PSIR(IGR)   = ROTB1A(IGR) ! new psir angle
            IGRSLAVE(IGR) = IGRSLA
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE BACK_FLIGHT

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE ALLOCA(MAXCO,MAXRB,NMAXTS,MAX2,MAXP2,IERRCODE)

C ... Allocate the space for the main arrays in each processor

C ... MAXB    is the size of the main variable arrays
C ... MAX2    is the size of the multigrid arrays
C ... MAXCO   is used if the grid is changed
C ... MAXTI   is for time-dependent main arrays
C ... MAXTS   is for time-dependent scalars
C ... NMAXSB  is the size for scalar quantities
C ... MAXRB   is allocated for the Reynolds stress model
C ... MAXMP   is the size of the multiphase arrays
C ... MAXEB   is the size of the EARSM (Wallin-Johansson)
C ... MAXSS   is used for shear stresses (STRESL == .TRUE.)
C ... MAXCH   is for Chimera
C ... MAX11   is the size of the largest block in this processor
C ... IB      is the size of the patches on the first level (partially obsolate)
C ... IBF     is the size of the patch arrays including multigrid
C ... MAXFS   is for free-surface model
C ... MAXMP   is for derived type multiphase
C ... MAXP2   is the multigrid of multiphase model
C ... MAXPC   is for pressure correction 

      USE NS3CO, ONLY : NSCAL, NBLOCG

      USE MAIN_ARRAYS
      USE INTEGERS, ONLY : MAXB,  MAXFS, MAXEB, MAXSB, MAXTI, 
     &                     MAXSS, MAXCH, MAXS2, MAXTS, MAX11, IB, IBF,
     &                     MAXMP, MAXPC, MAXTR
      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER :: MAXCO, IERRCODE, MAXSCH, MAXSS2, NSCA1, MAXRB, NMAXTS,
     &           MAX2, MAXP2
       
      ALLOCATE(XCO(MAXB),  YCO(MAXB) ,  ZCO(MAXB) ,  XLE2(MAXCO),
     2       YLE2(MAXCO), ZLE2(MAXCO), XLE3(MAXCO),  YLE3(MAXCO),
     3       ZLE3(MAXCO),    XC(MAXB),    YC(MAXB),     ZC(MAXB),
     4       XGRI(MAXCO), YGRI(MAXCO),  ZGRI(MAXCO),
     4       XORI(MAXCO), YORI(MAXCO),  ZORI(MAXCO),STAT=IERRCODE)
       
      ALLOCATE(VOL(MAXB),    D1(MAXB),    D2(MAXB),    D3(MAXB),
     1        UROT(MAXB),  VROT(MAXB),  WROT(MAXB), DISTW(MAXB),
     2         XFC(MAXB),   YFC(MAXB),   ZFC(MAXB),!AREA(3*MAXB),
     3      UROTCP(MAXB),VROTCP(MAXB),WROTCP(MAXB),   RLD(MAXB),
     4      RLDIST(MAXB), DISTP(MAXB),
     5    STAT=IERRCODE)

      ALLOCATE( A1(MAXB),    A2(MAXB),    A3(MAXB),  A1XA(MAXB),
     2        A1YA(MAXB),  A1ZA(MAXB),  A2XA(MAXB),  A2YA(MAXB),
     3        A2ZA(MAXB),  A3XA(MAXB),  A3YA(MAXB),  A3ZA(MAXB),
     4    STAT=IERRCODE)

      ALLOCATE(LOCDIS(MAXB),IDP(MAXB))
         
      ALLOCATE( RO(MAXB),    RM(MAXB),    RN(MAXB),    RW(MAXB),
     2           E(MAXB),   DRO(MAXB),    DM(MAXB),    DN(MAXB),
     3          DW(MAXB),    DE(MAXB),  TOLD(MAXB),  UOLD(MAXB),
     4        VOLD(MAXB),  WOLD(MAXB),  POLD(MAXB),  OHMI(MAXB),
     5           U(MAXB),     V(MAXB),     W(MAXB),     P(MAXB),
     6       PDIFF(MAXB),  EPS2(MAXB),  VIST(MAXB),     C(MAXB),
     7        TEMP(MAXB),   F1R(MAXB),  F1RM(MAXB),  F1RN(MAXB),
     8        F1RW(MAXB),   F1E(MAXB),  F1RK(MAXB), F1EPS(MAXB),
     9        HAT1(MAXB),  HAT2(MAXB),  HAT3(MAXB),  HAT4(MAXB),
     +         VIS(MAXB),   DTL(MAXB),    CP(MAXB),    CH(MAXB),
     1        DRDP(MAXB),  DRDH(MAXB),  RKSI(MAXB), BLANK(MAXB),
     2       RUVAV(MAXB), RUWAV(MAXB), RVWAV(MAXB),   TIJ(MAXB),     
     3       F1H(MAXFS), SUROLD(MAXFS), VORT(MAXB,3), SHEAR(MAXB,6),
     4       STRAIN(MAXB), BLANK2(MAXB),ROLD(MAXB),STAT=IERRCODE)
         
C ... MAXTB is the size of the turbulence arrays, replaced by MAXB

      ALLOCATE( RK(MAXB),  REPS(MAXB), DDEPS(MAXB),   DRK(MAXB),
     2        DEPS(MAXB), RKOLD(MAXB),EPSOLD(MAXB),   SRK(MAXB),
     3        SEPS(MAXB),  PTUR(MAXB),  FUN1(MAXB),   TTS(MAXB),
     4       VTRAN(MAXB),  RNUT(MAXB),  QSAS(MAXB),VELLAP(MAXB),
     5      STAT=IERRCODE)

      ALLOCATE(BIJ(MAXEB,6),WIR(MAXEB,3),STAT=IERRCODE)

      NSCA1 = MAX(NSCAL,1)  ! Minimum allocatable space

      ALLOCATE( FI(MAXSB,NSCA1),  DFI(MAXSB,NSCA1), FIOLD(MAXSB,NSCA1),
     2         SFI(MAXSB,NSCA1), F1FI(MAXSB,NSCA1), STAT=IERRCODE)

      ALLOCATE(PROD(6*MAXRB),  SPI(6*MAXRB),  DIF(6*MAXRB),
     2          DIS(6*MAXRB), VVIS(6*MAXRB),   FWLL(MAXRB), 
     3         STAT=IERRCODE)
       
      ALLOCATE(S11(MAXSS,6) , STAT=IERRCODE)

      MAXSCH  = AMAX0(MAXCH*NSCA1 ,1)
      MAXSS2  = AMAX0(MAXS2*NSCA1 ,1)
        
      ALLOCATE(ROLE2(MAXTI), RMLE2(MAXTI), RNLE2(MAXTI), RWLE2(MAXTI),
     2          ELE2(MAXTI), RKLE2(MAXTI),EPSLE2(MAXTI), ROLE3(MAXTI),
     3         RMLE3(MAXTI), RNLE3(MAXTI), RWLE3(MAXTI),  ELE3(MAXTI),
     4         RKLE3(MAXTI),EPSLE3(MAXTI), ROAV1(MAXTI), RMAV1(MAXTI),
     5         RNAV1(MAXTI), RWAV1(MAXTI),  EAV1(MAXTI), RKAV1(MAXTI),
     6        EPSAV1(MAXTI), ROAV2(MAXTI), RMAV2(MAXTI), RNAV2(MAXTI),
     7         RWAV2(MAXTI),  EAV2(MAXTI), RKAV2(MAXTI),EPSAV2(MAXTI),
     8         RMAV3(MAXTI), RNAV3(MAXTI), RWAV3(MAXTI),  PLE2(MAXTI),
     +          PLE3(MAXTI),  PAV1(MAXTI),  TAV1(MAXTI),
     9         FILE2(NMAXTS,NSCA1), FILE3(NMAXTS,NSCA1),
     1         FIAV1(NMAXTS,NSCA1), FIAV2(NMAXTS,NSCA1),  STAT=IERRCODE)

      ALLOCATE(ROFOR(MAXCH), RMFOR(MAXCH), RNFOR(MAXCH), RWFOR(MAXCH),
     2          EFOR(MAXCH), RKFOR(MAXCH), REFOR(MAXCH), PDFOR(MAXCH),
     3          WGH1(MAXCH),  WGH2(MAXCH),  WGH3(MAXCH),  WGH4(MAXCH),
     3         FIFOR(MAXCH,NSCA1),
     4         STAT=IERRCODE)

      ALLOCATE(  II1(MAXCH),   II2(MAXCH),   II3(MAXCH),   II4(MAXCH),
     2          INTP(MAXCH),STAT=IERRCODE)

      ALLOCATE( W12(3*MAX11),  SIJ(6*MAX11),GRADT(3*MAX11),
     2        GRADK(3*MAX11),GREPS(3*MAX11),STAT = IERRCODE)

      ALLOCATE(ROP2H(MAX2), RMP2H(MAX2) , RNP2H(MAX2) ,  RWP2H(MAX2),
     2          EP2H(MAX2), RKP2H(MAX2) ,EPSP2H(MAX2) ,  SURH2(MAX2),
     3  FIP2H(MAXS2,NSCA1),      STAT=IERRCODE)

      ALLOCATE( ICP(IB),  JCP(IB),  KCP(IB),IJMASK(IB),  JET(MAXB),
     +         STAT=IERRCODE)

      ALLOCATE( HFLUX(IBF),CPWALL(IBF), UWALL(IBF), VWALL(IBF),
     2          WWALL(IBF), TWALL(IBF), QWALL(IBF),QWFRIC(IBF),
     3          TAUW1(IBF), TAUW2(IBF), SURFX(IBF), SURFY(IBF),
     4          SURFZ(IBF),   XCP(IBF),   YCP(IBF),   ZCP(IBF),
     5         WMFLUX(IBF), POROS(IBF),WHSTAG(IBF), WTEMP(IBF),
     6         RSDIRX(IBF),RSDIRY(IBF),RSDIRZ(IBF),   RBK(IBF),
     7          TAUWX(IBF), TAUWY(IBF), TAUWZ(IBF), SURLE(IBF),
     8         DSURLE(IBF), WAVES(IBF), WAVEH(IBF),IWAVEB(MAXFS),
     9          WTRAN(IBF),RMLOSS(IBF), BOUNG(IBF),BOUNRET(IBF),
     +          SURFA(IBF), SURFT(IBF),SURF2X(IBF), SURF2Y(IBF),
     1         SURF2Z(IBF),SURFMX(IBF),SURFMY(IBF),SURFMZ(IBF),
     2         SURFPX(IBF),SURFPY(IBF),SURFPZ(IBF),UTAUM(IBF),
     3         STAT=IERRCODE)

      ALLOCATE(CPWALL_(IBF), QWALL_(IBF), QWFRIC_(IBF),
     &         TWALL_(IBF),  TAUW1_(IBF), TAUW2_(IBF),
     &         TAUWX_(IBF),  TAUWY_(IBF), TAUWZ_(IBF),
     &         SURFX_(IBF),  SURFY_(IBF), SURFZ_(IBF),
     &         HFLUX_(IBF),  SURLE_(IBF), DSURLE_(IBF),
     &         F1RK_(IBF),   WTRAN_(IBF), WMFLUX_(IBF),
     &         STAT=IERRCODE)
      
      
      ALLOCATE(BOUNMF(IBF), BOUNU(IBF), BOUNV(IBF), BOUNW(IBF),
     2         BOUNT(IBF),  BOUNP(IBF), BOUNR(IBF), BOUNE(IBF),
     3         BOUNRK(IBF), BOUNEP(IBF),BOUNPD(IBF),BOUNA1(IBF),
     4         BOUNA2(IBF), BOUNU1(IBF),BOUNU2(IBF),BOUNV1(IBF),
     5         BOUNV2(IBF), BOUNW1(IBF),BOUNW2(IBF),BOUNFI(IBF,NSCA1),
     6         BOUNBI(IBF,6), STAT=IERRCODE)

      ALLOCATE(BOUNR_(IBF), BOUNT_(IBF), BOUNP_(IBF), BOUNRK_(IBF),
     &         BOUNEP_(IBF), BOUNMF_(IBF), BOUNU_(IBF), BOUNV_(IBF),
     &         BOUNW_(IBF), BOUNE_(IBF), STAT=IERRCODE)

      
C ... Obsolate arrays

      ALLOCATE(    UBI(IB),   VBI(IB),  WBI(IB),   UBJ(IB),   VBJ(IB),
     2             WBJ(IB),   UBK(IB),  VBK(IB),   WBK(IB),   UTI(IB),
     3             VTI(IB),   WTI(IB),  UTJ(IB),   VTJ(IB),   WTJ(IB),
     4             UTK(IB),   VTK(IB),  WTK(IB), STAT=IERRCODE)


      ALLOCATE(                XXI(MAXFS),   YXI(MAXFS),   ZXI(MAXFS),
     2          XETA(MAXFS),  YETA(MAXFS),   WFS(MAXFS),    XI(MAXFS),
     3           ETA(MAXFS),  ZETA(MAXFS), XHULL(MAXFS), YHULL(MAXFS),
     4         ZHULL(MAXFS),   XVL(MAXFS),   YVL(MAXFS),   ZVL(MAXFS),
     5         ZZTOP(MAXFS), UHULL(MAXFS), VHULL(MAXFS), WHULL(MAXFS),
     6         IHULL(MAXFS), JHULL(MAXFS),   XSP(MAXFS),   YSP(MAXFS),
     7           ZSP(MAXFS),      STAT=IERRCODE)

      ALLOCATE(WH(IBF)) ! 


C ... Derived type arrays for multiphase and pressure correction

      ALLOCATE(  PRO(MAXMP),  VAR(MAXMP), P2H(MAXP2)) ! Oli MAXB ??

      ALLOCATE(  MPFOR(MAXMP))

      ALLOCATE(  PRC(MAXPC))

      ALLOCATE(  SDI(MAXB))

      ALLOCATE(  TRM(MAXTR), P2HTRM(MAXTR))

      RETURN
      END SUBROUTINE ALLOCA

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE ALLONS(NPPN)

      USE NS3CO, ONLY : IC9, NRESI, NBLOCK, NBLOCG
      USE INTEGERS, ONLY : MGM, NBCS, IPRO, NPPV, NB, 
     &                     MBPRO, NPRO, NBGG
      USE MAIN_ARRAYS
      USE CHARACTERS, ONLY : BOUNDF, BOUNDN

      IMPLICIT NONE

      INTEGER :: NPPN, N, M, I

C ... This subroutine allocates the old NS3CO.CC variables

C ... NB     IS THE MAXIMUM NUMBER OF BLOCKS IN ONE PROCESS
C ... NBGG   IS THE MAXIMUM NUMBER OF BLOCKS IN ALL PROCESSES
C ... NPPV   IS THE MAXIMUM NUMBER OF PATCHES/BLOCK
C ... NPPN   IS THE MAXIMUM NUMBER OF PATCHES
C ... MAXNSC IS THE MAXIMUM NUMBER OF SCALAR VARIABLES
C ... MGM    IS THE MAXIMUM NUMBER OF MULTIGRIDLEVELS
C ... MNOFBC IS THE MAXIMUM NUMBER OF BOUNDARY-PATCHES
C ... MBPRO  IS THE MAXIMUM NUMBER OF BLOCKS IN PROCESSOR


      NB    = NBLOCK + 1
      NBGG  = NBLOCG + 1
      MBPRO = NBLOCG            ! NBLOCG is bit too large but...
      NPPN  = NBCS + 1

      WRITE(45,111) NBLOCG,NBLOCK,MGM,NBGG,NB,NPPV,NPPN,MBPRO,NPRO
c 111  FORMAT(/'Sizes of allocated NS3CO. variables:'/
c     + ' NBLOCG = ',I4,' NBLOCK = ',I4,' MGM    = ',I4/
c     + ' NBGG   = ',I4,' NB     = ',I4,' NPPV   = ',I4/
c     + ' NPPN   = ',I4,' MBPRO  = ',I4,' NPRO = ',I4)
 111  FORMAT(/'Sizes of allocated NS3CO. variables:'/
     + ' NBLOCG = ',I7,' NBLOCK = ',I7,' MGM    = ',I7/
     + ' NBGG   = ',I7,' NB     = ',I7,' NPPV   = ',I7/
     + ' NPPN   = ',I7,' MBPRO  = ',I7,' NPRO = ',I7)

      IF(IPRO /= 1) THEN
         ALLOCATE(NPROCE(MBPRO+1,NPRO))
         ALLOCATE(IMAXG(NBGG))
         ALLOCATE(JMAXG(NBGG))
         ALLOCATE(KMAXG(NBGG))
         ALLOCATE(IDOB(NBGG))
         ALLOCATE(JDOB(NBGG))
         ALLOCATE(KDOB(NBGG))
         ALLOCATE(NSSB(NBGG))
         ALLOCATE(ISSB(7,NBGG))
         ALLOCATE(NCHIMT(NBGG))
         ALLOCATE(NCHOR(NBGG))
         ALLOCATE(NLOCAL(NBGG))
         ALLOCATE(NPNUM(NBGG))
         ALLOCATE(OMEGA(NBGG))
         ALLOCATE(OMEGAX(NBGG))
         ALLOCATE(OMEGAY(NBGG))
         ALLOCATE(OMEGAZ(NBGG))
         ALLOCATE(CENAX(NBGG))
         ALLOCATE(CENAY(NBGG))
         ALLOCATE(CENAZ(NBGG))
         ALLOCATE(AMPL(NBGG))
         ALLOCATE(OSCLS(NBGG))
         ALLOCATE(BXMIN(NBGG))
         ALLOCATE(BXMAX(NBGG))
         ALLOCATE(BYMIN(NBGG))
         ALLOCATE(BYMAX(NBGG))
         ALLOCATE(BZMIN(NBGG))
         ALLOCATE(BZMAX(NBGG))
         ALLOCATE(IGRID(NBGG,5))
         ALLOCATE(NEARBLOCK(NBGG*NBGG))

         ALLOCATE(JSTATE(NBLOCG,3))
         ALLOCATE(MGRIDA(NBLOCG))
         ALLOCATE(BLKS(NBGG))

         ALLOCATE(ICNH(NBGG))
         ALLOCATE(IBTGR(NBGG))
         ALLOCATE(ITPGR(NBGG))
         ALLOCATE(JBTGR(NBGG))
         ALLOCATE(JTPGR(NBGG))
         ALLOCATE(KBTGR(NBGG))
         ALLOCATE(KTPGR(NBGG))
         ALLOCATE(LHULL(NBGG))
         ALLOCATE(MMHUL(NBGG))
         ALLOCATE(MHULL(NBGG))
         ALLOCATE(ICGH(NBGG))
         ALLOCATE(ICOGH(NBGG))
         ALLOCATE(ICONH(NBGG))
         ALLOCATE(IBOTGR(NBGG))
         ALLOCATE(ITOPGR(NBGG))
         ALLOCATE(JBOTGR(NBGG))
         ALLOCATE(JTOPGR(NBGG))
         ALLOCATE(KBOTGR(NBGG))
         ALLOCATE(KTOPGR(NBGG))
         ALLOCATE(IUPPT(NBGG))

         NPROCE = 0
         IMAXG  = 0
         JMAXG  = 0
         KMAXG  = 0
         IDOB   = 0
         JDOB   = 0
         KDOB   = 0
         NSSB   = 0
         ISSB   = 0
         NCHIMT = 0
         NCHOR  = 0
         NLOCAL = 0
         NPNUM  = 0
         OMEGAX = 0.
         OMEGAY = 0.
         OMEGAZ = 0.
         CENAX  = 0.
         CENAY  = 0.
         CENAZ  = 0.
         AMPL   = 0.
         OSCLS  = 0
         IGRID  = 0
         JSTATE = 0
         MGRIDA = 0
         ICNH   = 0
         IBTGR  = 0
         ITPGR  = 0
         JBTGR  = 0
         JTPGR  = 0
         KBTGR  = 0
         KTPGR  = 0
         LHULL  = 0
         MMHUL  = 0
         MHULL  = 0
         ICGH   = 0
         ICONH  = 0
         IBOTGR = 0
         ITOPGR = 0
         JBOTGR = 0
         JTOPGR = 0
         KBOTGR = 0
         KTOPGR = 0
      ENDIF

      ALLOCATE(NCPAT(NPPV*NB),IT(MGM,NB),
     1 IL(MGM,NB),IK(MGM,NB),NTOT(MGM,NB),
     2 IMAX(MGM,NB),JMAX(MGM,NB),KMAX(MGM,NB),
     3 IDI1(NB),IDI2(NB),IDI3(NB),
     4 INTERI(NB),INTERJ(NB),INTERK(NB),IDER(NB),
     6 IBOT(MGM,NB),ITOP(MGM,NB),
     7 JBOT(MGM,NB),JTOP(MGM,NB),
     8 KBOT(MGM,NB),KTOP(MGM,NB),
     9 MGRID(NB),LAMIN(NB),
     1 JF(MGM,NB+1),IG(MGM,NB+1),IH(MGM,NB+1),
     2 IR(MGM,NB+1),IQ(MGM,NB+1),
     1 MIB(NB),MIT(NB),MJB(NB),MJT(NB),
     2 MKB(NB),MKT(NB),
     2 INITC(NB),ITAG(NPPV,NB),IHF(NPPN,MGM+1),
     3 IW(NPPV,NB,MGM+1),ICON(IC9*NPPN*MGM),IC(MGM,NB),
     4 NPATCH(NB),NSOLPA(NB),
     5 NORMAL(NPPN),ISTRS(NPPN),JSTRS(NPPN),KSTRS(NPPN),
     6 KX1S(NPPN),KX2S(NPPN),KY1S(NPPN),KY2S(NPPN),
     7 KZ1S(NPPN),KZ2S(NPPN),IDIMS(NPPN),JDIMS(NPPN),
     8 NSPB(NPPN),NSPP(NPPN),NSPG(NPPN),
     9 MOV(NB),IROTVE(NB),
     1 JAPP(NPPV,NB),
     2 IWAPP(NPPV,NB+1),
     3 IDIMSG(NPPN),JDIMSG(NPPN))

      DO I = 1,IC9*NPPN*MGM
         ICON(I) = 0
      ENDDO

      DO N = 1,NB ! needed for initialization for some reason PPR 13.11.00
         DO M=1,MGM
            DO I=1,NPPV
               IW(I,N,M)   = 0
               IW(I,N,M+1) = 0
            ENDDO
            IMAX(M,N) = 0
            JMAX(M,N) = 0
            KMAX(M,N) = 0
            NTOT(M,N) = 0
            IG(M,N)   = 0
            IG(M,N+1) = 0
         ENDDO
      ENDDO

      ALLOCATE(BOUNDF(NPPN),BOUNDN(NPPN))

      ALLOCATE(CXB(NPPN), CYB(NPPN),CZB(NPPN),
     1  CMXB(NPPN),CMYB(NPPN),CMZB(NPPN),
     1  DXB(NPPN), DYB(NPPN), DZB(NPPN),QTB(NPPN), QWB(NPPN), QHB(NPPN),
     2  XCOL(3,4,NPPN),
     6  TOMEGA(NPPN),AAA(9,NPPN),RESI(NRESI,NB),APATCH(NPPN),
     7  XCOG(NPPN*12),RCON(IC9*NPPN*MGM),VOLN(NB),XCOR(6,NB),
     8  BUX(NPPN),BUY(NPPN),BUZ(NPPN),VMXB(NPPN),VMYB(NPPN),VMZB(NPPN))

C ... Allocation for the free-surface model

      ALLOCATE(XC1(NPPN), YC1(NPPN),ZC1(NPPN),
     1         XC2(NPPN), YC2(NPPN),ZC2(NPPN),
     1         XC3(NPPN), YC3(NPPN),ZC3(NPPN),
     1         XC4(NPPN), YC4(NPPN),ZC4(NPPN))

      WRITE(45,*) 'ALLONS SUCCESFULLY EXECUTED'

      END SUBROUTINE ALLONS

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE ZEROMA

      USE MAIN_ARRAYS
      USE INTEGERS, ONLY : IREPEA, NREPEA
      USE FLIGHT, ONLY : NGRIFL

      IMPLICIT NONE

      INTEGER :: IPHASE 

C ... Put all the work arrays to zeroes. (Maybe this could be removed)

      UWALL  = 0.
      HFLUX  = 0.
      CPWALL = 0.
      VWALL  = 0.
      WWALL  = 0.
      TWALL  = 0.
      QWALL  = 0.
      QWFRIC = 0.
      TAUW1  = 0.
      TAUW2  = 0.
      TAUWX  = 0.
      TAUWY  = 0.
      TAUWZ  = 0.
      SURFX  = 0.
      SURFY  = 0.
      SURFZ  = 0.
      SURFPX = 0.
      SURFPY = 0.
      SURFPZ = 0.
      SURF2X = 0.
      SURF2Y = 0.
      SURF2Z = 0.
      SURFMX = 0.
      SURFMY = 0.
      SURFMZ = 0.
      UTAUM  = 10. ! To input?
      WMFLUX = 0.
      POROS  = 0.
      RMLOSS = 0.
      WHSTAG = 0.
      WTEMP  = 0.
      RSDIRX = 0.
      RSDIRY = 0.
      RSDIRZ = 0.
      SURLE  = 0.
      DSURLE = 0.
      RBK    = 0.

      FI     = 0.
      SIJ    = 0.
      S11    = 0.

      A1     = 0.
      A2     = 0.
      A3     = 0.
      A1XA   = 1.  ! This is ancient history
      A1YA   = 0.
      A1ZA   = 0.
      A2XA   = 0.
      A2YA   = 1.
      A2ZA   = 0.
      A3XA   = 0.
      A3YA   = 0.
      A3ZA   = 1.
      D1     = 0.
      D2     = 0.
      D3     = 0.
      XC     = 0.
      YC     = 0.
      ZC     = 0.
      UROT   = 0.
      VROT   = 0.
      WROT   = 0.
      DISTW  = 0.
      DISTP  = 0.
      RLD    = 0.
      RLDIST = 0.
      SEPS   = 0.
      SIJ    = 0.
      PTUR   = 0.
      XCO    = 0.
      YCO    = 0.
      ZCO    = 0.
      XSP    = 0.
      YSP    = 0.
      ZSP    = 0.
      LOCDIS = 0
      IDP    = 0

      ICP    = 0
      JCP    = 0
      KCP    = 0
      IJMASK = 0
      IREPEA = 0
      NREPEA = 10
      NREPEA(1) = 5
      NREPEA(2) = 6
      NREPEA(3) = 4
      NREPEA(6) = NGRIFL  ! Rotor time step control
      NREPEA(8) = 20      ! Maximum number of NATRAP calls
      NREPEA(10)= 100     ! Warning for exceeding AGM cycle limit
      NREPEA(11)= 20      ! Maximum umber of low-pressure warnings

      RK     = 0.
      REPS   = 0.
      DDEPS  = 0.
      DRK    = 0.
      DEPS   = 0.
      SRK    = 0.
      SEPS   = 0.
      PTUR   = 0.
      FUN1   = 0.
      F1RK   = 0.
      F1EPS  = 0.
      TTS    = 0.
      VTRAN  = 1.

      BIJ    = 0.
      WIR    = 0.

      FI     = 0.
      F1FI   = 0.
      DFI    = 0.
      SFI    = 0.

      HAT1   = 0.
      HAT2   = 0.
      HAT3   = 0.
      HAT4   = 0.
      
      PROD   = 0.
      SPI    = 0.
      DIF    = 0.
      DIS    = 0.
      VVIS   = 0.
      FWLL   = 0.

      ROLE2  = 0.
      RMLE2  = 0.
      RNLE2  = 0.
      RWLE2  = 0.
      ELE2   = 0.
      RKLE2  = 0.
      EPSLE2 = 0.
      FILE2  = 0.
      XLE2   = 0.
      YLE2   = 0.
      ZLE2   = 0.

      ROLE3  = 0.
      RMLE3  = 0.
      RNLE3  = 0.
      RWLE3  = 0.
      ELE3   = 0.
      RKLE3  = 0.
      EPSLE3 = 0.
      FILE3  = 0.
      XLE3   = 0.
      YLE3   = 0.
      ZLE3   = 0.
      
      ROFOR  = 0.
      RMFOR  = 0.
      RNFOR  = 0.
      RWFOR  = 0.
      EFOR   = 0.
      PDFOR  = 0.
      RKFOR  = 0.
      REFOR  = 0.
      RKSI   = 0.
      FIFOR  = 0.

      TRM%GFOR   = 0.
      TRM%RETFOR = 0.

      BOUNMF = 0.
      BOUNR  = 0.
      BOUNU  = 0.
      BOUNV  = 0.
      BOUNW  = 0.
      BOUNE  = 0.
      BOUNT  = 0.
      BOUNP  = 0.
      BOUNRK = 0.
      BOUNEP = 0.
      BOUNFI = 0.
      BOUNBI = 0.
      APATCH = 0.

      WH     = 0.
      XVL    = 0.
      YVL    = 0.
      ZVL    = 0.
      F1H    = 0.
      SUROLD = 0.
      SURH2  = 0.
      WAVEH  = 0.
      WAVES  = 0.
      IWAVEB = 0
      WFS    = 0.
      XHULL  = 0.
      YHULL  = 0.
      ZHULL  = 0.
      LHULL  = 0
      MMHUL  = 0
      MHULL  = 0
      IHULL  = 0
      IBOTGR = 0
      JBOTGR = 0
      KBOTGR = 0
      ITOPGR = 0
      JTOPGR = 0
      KTOPGR = 0
      IDIMS  = 0
      JDIMS  = 0
      IDIMSG = 0
      JDIMSG = 0
      NEARBLOCK = 0


      DO IPHASE = 1,NPHASES

         PRO%TEMP(IPHASE) = 0.
         PRO%E(IPHASE)    = 0.
         PRO%DRODP(IPHASE)= 0.
         PRO%DRODH(IPHASE)= 0.
         PRO%VIS(IPHASE)  = 0.
         PRO%CH(IPHASE)   = 0.
         PRO%C(IPHASE)    = 0.
         PRO%QIF(IPHASE)  = 0.
         PRO%HSAT(IPHASE) = 0.
         PRO%DHSDP(IPHASE)= 0.
         PRO%DHSDH(IPHASE)= 0.
         PRO%CP(IPHASE)   = 0.
         PRO%TAUF(IPHASE) = 0.
         PRO%RO(IPHASE)   = 0.
C        PRO%N(IPHASE)    = 0. Should be solved      
         VAR%DTEMP(IPHASE)= 0.
         VAR%ALFA(IPHASE) = 0.
         VAR%X(IPHASE)    = 0.
         VAR%DX(IPHASE)   = 0.
         VAR%XOLD(IPHASE) = 0.
         VAR%EVAPR(IPHASE)= 0.
         VAR%F1R(IPHASE)  = 0.
         VAR%F2R(IPHASE)  = 0.
         VAR%F3R(IPHASE)  = 0.

         MPFOR%DTEMPFOR(IPHASE) = 0.
         MPFOR%ALFAFOR(IPHASE)  = 0.
         MPFOR%XFOR(IPHASE)     = 0.

      ENDDO

      PRC%F1R  = 0.
      PRC%F2R  = 0.
      PRC%F3R  = 0.
      PRC%FIR  = 0.
      PRC%FJR  = 0.
      PRC%FKR  = 0.
      PRO%TSAT = 0.
      PRO%PSAT = 0.
      
      RETURN
      END SUBROUTINE ZEROMA

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE DEALNS

      USE MAIN_ARRAYS
      USE CHARACTERS, ONLY : BOUNDF, BOUNDN

      IMPLICIT NONE

      DEALLOCATE(NPROCE)
      DEALLOCATE(JSTATE)
      DEALLOCATE(IMAXG)
      DEALLOCATE(JMAXG)
      DEALLOCATE(KMAXG)
      DEALLOCATE(IDOB)
      DEALLOCATE(JDOB)
      DEALLOCATE(KDOB)
      DEALLOCATE(NSSB)
      DEALLOCATE(ISSB)
      DEALLOCATE(NCHIMT)
      DEALLOCATE(NCHOR)
      DEALLOCATE(NLOCAL)
      DEALLOCATE(NPNUM)
      DEALLOCATE(OMEGA)
      DEALLOCATE(OMEGAX)
      DEALLOCATE(OMEGAY)
      DEALLOCATE(OMEGAZ)
      DEALLOCATE(CENAX)
      DEALLOCATE(CENAY)
      DEALLOCATE(CENAZ)
      DEALLOCATE(AMPL)
      DEALLOCATE(OSCLS)
      DEALLOCATE(BXMIN)
      DEALLOCATE(BXMAX)
      DEALLOCATE(BYMIN)
      DEALLOCATE(BYMAX)
      DEALLOCATE(BZMIN)
      DEALLOCATE(BZMAX)
      DEALLOCATE(IGRID)
      DEALLOCATE(BLKS)
      DEALLOCATE(NEARBLOCK)
      
      DEALLOCATE(ICNH)
      DEALLOCATE(IBTGR)
      DEALLOCATE(ITPGR)
      DEALLOCATE(JBTGR)
      DEALLOCATE(JTPGR)
      DEALLOCATE(KBTGR)
      DEALLOCATE(KTPGR)
      DEALLOCATE(LHULL)
      DEALLOCATE(MMHUL)
      DEALLOCATE(MHULL)
      DEALLOCATE(ICGH)
      DEALLOCATE(ICOGH)
      DEALLOCATE(ICONH)
      DEALLOCATE(IBOTGR)
      DEALLOCATE(ITOPGR)
      DEALLOCATE(JBOTGR)
      DEALLOCATE(JTOPGR)
      DEALLOCATE(KBOTGR)
      DEALLOCATE(KTOPGR)
      DEALLOCATE(IUPPT)
      
      DEALLOCATE(NCPAT,IT,IL,IK,NTOT,IMAX,JMAX,KMAX,
     1 IDI1,IDI2,IDI3,INTERI,INTERJ,INTERK,IDER,
     2 IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     3 MGRID,LAMIN,JF,IG,IH,IR,IQ,
     4 MIB,MIT,MJB,MJT,MKB,MKT,
     5 INITC,ITAG,IHF,IW,ICON,IC,NPATCH,NSOLPA,NORMAL,
     6 ISTRS,JSTRS,KSTRS,KX1S,KX2S,KY1S,KY2S,KZ1S,KZ2S,
     7 IDIMS,JDIMS,NSPB,NSPP,NSPG,
     8 MOV,IMAXM,JMAXM,KMAXM,MOVPO,IROTVE,IMINM,JMINM,KMINM,
     9 JAPP,IWAPP,IDIMSG,JDIMSG)

      DEALLOCATE(BOUNDF,BOUNDN)

      DEALLOCATE(CXB,CYB,CZB,CMXB,CMYB,CMZB,DXB,DYB,DZB,QTB,QWB,QHB,
     2  XCOL,TOMEGA,AAA,RESI,XCOG,RCON,VOLN,XCOR,APATCH,BUX,BUY,BUZ,
     3  VMXB,VMYB,VMZB)

      DEALLOCATE(RLOLIM,UPPLIM,PSIGSC,PSIGS2)

      DEALLOCATE(XC1,YC1,ZC1,XC2,YC2,ZC2,XC3,YC3,ZC3,XC4,YC4,ZC4)

      RETURN
      END SUBROUTINE DEALNS

C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C

      SUBROUTINE DEALLO

      USE MAIN_ARRAYS

      IMPLICIT NONE

      INTEGER :: IERRCODE
     
      DEALLOCATE(XCO,   YCO,   ZCO,  XLE2,  YLE2,  ZLE2,  XLE3,  YLE3,
     2          ZLE3,    XC,    YC,    ZC,
     3 STAT=IERRCODE)

      DEALLOCATE(VOL,  D1,  D2,  D3,  DISTW, DISTP, UROT, VROT, WROT, 
     2           XFC, YFC, ZFC, RLD, RLDIST, STAT=IERRCODE)

      DEALLOCATE(LOCDIS, IDP, STAT=IERRCODE)

      DEALLOCATE( A1,    A2,    A3,  A1XA,  A1YA,  A1ZA,  A2XA,  A2YA, 
     3          A2ZA,  A3XA,  A3YA,  A3ZA, STAT=IERRCODE)

      DEALLOCATE( RO,    RM,    RN,    RW,     E,   DRO,    DM,    DN,
     2            DW,    DE,  TOLD,  UOLD,  VOLD,  WOLD,  POLD,  OHMI,
     3             U,     V,     W,     P, PDIFF,  EPS2,  VIST,     C,
     4          TEMP,   F1R,  F1RM,  F1RN,  F1RW,   F1E,  F1RK, F1EPS,
     5          HAT1,  HAT2,  HAT3,  HAT4,   TIJ,
     6          DRDP,  DRDH,  RKSI, BLANK, RUVAV, RUWAV, RVWAV,   F1H,
     7        SUROLD,  VORT, SHEAR,STRAIN, BLANK2,  ROLD, STAT=IERRCODE)

      DEALLOCATE(BIJ,   WIR, STAT=IERRCODE)

      DEALLOCATE( RK,  REPS, DDEPS,   DRK,  DEPS, RKOLD,EPSOLD,   SRK,
     2          SEPS,  PTUR,  FUN1,   TTS, VTRAN, RNUT,
     3 STAT=IERRCODE)

      DEALLOCATE(  FI,  DFI,FIOLD,  SFI, F1FI,       STAT=IERRCODE)

      DEALLOCATE(PROD,  SPI,  DIF,  DIS, VVIS, FWLL, STAT=IERRCODE)

      DEALLOCATE( S11, STAT=IERRCODE)

      DEALLOCATE(ROLE2, RMLE2, RNLE2, RWLE2,  ELE2, RKLE2,EPSLE2,
     2           ROLE3, RMLE3, RNLE3, RWLE3,  ELE3, RKLE3,EPSLE3,
     3           ROAV1, RMAV1, RNAV1, RWAV1,  EAV1, RKAV1,EPSAV1,
     4           ROAV2, RMAV2, RNAV2, RWAV2,  EAV2, RKAV2,EPSAV2,
     5           FILE2, FILE3, FIAV1, FIAV2, SURH2,  PAV1,  TAV1,
     6           STAT=IERRCODE)

      DEALLOCATE(ROFOR, RMFOR, RNFOR, RWFOR,  EFOR, RKFOR, REFOR,
     2           PDFOR, FIFOR,  WGH1,  WGH2,  WGH3,  WGH4,
     3 STAT= IERRCODE)

      DEALLOCATE(  II1,   II2,   II3,   II4,  INTP, STAT= IERRCODE)

      DEALLOCATE(  W12,   SIJ, GRADT, GRADK, GREPS, STAT=IERRCODE)

      DEALLOCATE(ROP2H, RMP2H, RNP2H, RWP2H,  EP2H, RKP2H,EPSP2H,
     2           FIP2H, STAT=IERRCODE)

      DEALLOCATE(  ICP,   JCP,   KCP, IJMASK,   JET, STAT=IERRCODE)

      DEALLOCATE(HFLUX,CPWALL, UWALL, VWALL, WWALL, TWALL, QWALL,
     2          QWFRIC, TAUW1, TAUW2, SURFX, SURFY, SURFZ,   XCP,
     3             YCP,   ZCP,WMFLUX, POROS,WHSTAG, WTEMP,RSDIRX,
     4          RSDIRY,RSDIRZ,   RBK, TAUWX, TAUWY, TAUWZ, SURLE,
     5          DSURLE,BOUNMF, BOUNU, BOUNV, BOUNW, BOUNT, BOUNP,
     6           BOUNR, BOUNE,BOUNRK,BOUNEP,BOUNFI,BOUNBI,BOUNPD,
     7          BOUNA1,BOUNA2, WAVEH, WAVES,IWAVEB, WTRAN,RMLOSS,
     8           BOUNG,BOUNRET,SURFA, SURFT,SURF2X,SURF2Y,SURF2Z,
     9          SURFMX,SURFMY,SURFMZ,BOUNU1,BOUNU2,BOUNV1,BOUNV2,
     +          BOUNW1,BOUNW2,SURFPX,SURFPY,SURFPZ,UTAUM, 
     9 STAT=IERRCODE)

      DEALLOCATE(CPWALL_, QWALL_, QWFRIC_,
     &           TWALL_,  TAUW1_, TAUW2_,
     &           TAUWX_,  TAUWY_, TAUWZ_,
     &           SURFX_,  SURFY_, SURFZ_,
     &           HFLUX_,  SURLE_, DSURLE_,
     &           F1RK_,   WTRAN_, WMFLUX_,
     &           STAT=IERRCODE)

      DEALLOCATE(BOUNR_,  BOUNT_,  BOUNP_, BOUNRK_,
     &           BOUNEP_, BOUNMF_, BOUNU_, BOUNV_,
     &           BOUNW_,  BOUNE_, STAT=IERRCODE)
      
      DEALLOCATE(  UBI,   VBI,   WBI,   UBJ,   VBJ,   WBJ,   UBK,
     2             VBK,   WBK,   UTI,   VTI,   WTI,   UTJ,   VTJ,
     3             WTJ,   UTK,   VTK,   WTK, STAT=IERRCODE)

      DEALLOCATE(   WH,   XXI,   YXI,   ZXI,  XETA,  YETA,   WFS,
     2               XI,  ETA,  ZETA, XHULL, YHULL, ZHULL,   XVL,   
     3              YVL,   ZVL,ZZTOP, UHULL, VHULL, WHULL,IHULL,
     4            JHULL,   XSP,  YSP,   ZSP,        STAT=IERRCODE)

      DEALLOCATE(  PRO, VAR, MPFOR, P2H, STAT=IERRCODE)

      DEALLOCATE(  PRC, STAT=IERRCODE)

      DEALLOCATE(  SDI, STAT=IERRCODE)

      DEALLOCATE(  TRM,  P2HTRM, STAT=IERRCODE)

      RETURN
      END SUBROUTINE DEALLO


 
