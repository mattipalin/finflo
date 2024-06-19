      MODULE CHARACTERS

! ... INCLUDE FILE FOR THE 3-D FLOW, 2.2.1990
! ... Modified for module 19.9.2002 (mersu)
! ... Rewritten 18.1.2008

      IMPLICIT NONE

      CHARACTER(LEN=80), ALLOCATABLE,DIMENSION(:) :: BOUNDF, BOUNDN

      CHARACTER(LEN=8), PARAMETER ::
     &   RO2C   = "   RO2  ", RM22C  = "   RM2  ", E2C    = "   E2   ",
     &   DRONC  = "  DRON  ", DMNC   = "   DMN  ", DENC   = "   DEN  ",
     &   STRAIC = " STRAIN ", ETAC   = "  ETA   ", TEPSC  = "TRUE EPS",
     &   DSURFC = "  DSURF ", WAVEC  = " WAVES  ", DISTWC = "  DISTW ",
     &   PROKC  = "   PROK ", PROVC  = "   PROV ", FUN1C  = "   FUN1 ",
     &   APC    = "   AP   ", DAPC   = "   DAP  ", UPC    = "   UP   ",
     &   PSATC  = "  Psat  ", TSATC  = "  Tsat  ", DPSATC = " P-Psat ",
     &   PSTAGC = " P_STAG ", TSTAGC = " T_STAG ", DPSDTC = "dPsat/dT",
     &   AIDXC  = "  AIDX  ", TIC    = "    TI  ",
     &   VISC   = "   VIS  ", VISTC  = "  VIST  ",
     &   SRKC   = "   SRK  ", SEPSC  = "  SEPS  ",
     &   PTURC  = "  PTUR  ", PTUREC = "PTUR/EPS",
     &   DRDPC  = "   DRDP ", DRDHC  = "   DRDH ",
     &   CPC    = "    Cp  ", CHC    = "    CH  ",
     &   EPS2C  = "  EPS2  ", EPS2XC = "  EPS2X ",
     &   CROSDC = " CROSSD ", TTSC   = "T.TIM-SC", TIJC   = "  T_ij  "

      CHARACTER(LEN=8), PARAMETER ::
     &   CC     = "    C   ", ROC    = "   RO   ", EC     = "    E   ",
     &   UC     = "    U   ", VC     = "    V   ", WC     = "    W   ",
     &   RMC    = "   RM   ", RNC    = "   RN   ", RWC    = "   RW   ",
     &   PC     = "    P   ", PDIFC  = "  PDIF  ", TEMPC  = "  TEMP  ",
     &   RKC    = "   RK   ", REPSC  = " EPSILON", RKSIC  = "  RKSI  ",
     &   DROC   = "   DRO  ", DEC    = "   DE   ", DHC    = "   DH   ",
     &   DUC    = "   DU   ", DVC    = "   DV   ", DWC    = "   DW   ",
     &   DMC    = "   DM   ", DNC    = "   DN   ", DC     = "    D   ",
     &   DPC    = "   DP   ", DTEMPC = "  DTEMP ", VTRANC = "  VTRAN ",
     &   DKC    = "   DK   ", DEPSC  = "  DEPS  ", WTRANC = "  WTRAN ",
     &   DGC    = "   DG   ", DRETC  = "  DRET  ",
     &   ROFMC  = "  ROFM  ", REFMC  = "  REFM  ",
     &   RMFMC  = "  RMFM  ", RNFMC  = "  RNFM  ",
     &   ROFNC  = "  ROFN  ", REFNC  = "  REFN  ",
     &   RMFNC  = "  RMFN  ", RNFNC  = "  RNFN  ",
     &   SURFBXC= " SURFBX ", SURFBYC= " SURFBY ", SURFBZC= " SURFBZ ",
     &   UROTBC = " ROTBAL ", RNUTC  = "  RNUT  ", NUTC   = "  NUT   ",
     &   DPDXC  = "  DPDX ", DPDYC  = "  DPDY  ", DPDZC  = "  DPDZ  ",
     &   CIFXC  = "  CIFX ", CIFYC  = "  CIFY  ", CIFZC  = "  CIFZ  "

      CHARACTER(LEN=8), PARAMETER ::
     &   ROLDC  = "  ROLD  ", EOLDC  = "  EOLD  ", POLDC  = "  POLD  ",
     &   UOLDC  = "  UOLD  ", VOLDC  = "  VOLD  ",
     &   RMOLDC = "  RMOLD ", RNOLDC = "  RNOLD "

      CHARACTER(LEN=8), PARAMETER ::               VOLC   = "   VOL  ",
     &   A1C    = "   A1   ", A2C    = "   A2   ", A3C    = "   A3   ",
     &   A1XC   = "  A1XA  ", A2XC   = "  A2XA  ", A3XC   = "  A3XA  ",
     &   A1YC   = "  A1YA  ", A2YC   = "  A2YA  ", A3YC   = "  A3YA  ",
     &   A1ZC   = "  A1ZA  ", A2ZC   = "  A2ZA  ", A3ZC   = "  A3ZA  "

      CHARACTER(LEN=8), PARAMETER ::                    
     &   FDXC   = "  FDX   ", FDYC   = "  FDY   ", FDZC   = "  FDZ   ",
     &   FLXC   = "  FLX   ", FLYC   = "  FLY   ", FLZC   = "  FLZ   ",
     &   FTDXC  = "  FTDX  ", FTDYC  = "  FTDY  ", FTDZC  = "  FTDZ  ",
     &   FWXC   = "  FWX   ", FWYC   = "  FWY   ", FWZC   = "  FDZ   ",
     &   FVMXC  = "  FVMX  ", FVMYC  = "  FVMY  ", FVMZC  = "  FVMZ  "

      CHARACTER(LEN=8), PARAMETER ::               OHMIC  = "  VORT  ",
     &   OHMIXC = " VORTX  ", OHMIYC = " VORTY  ", OHMIZC = " VORTZ  ",
     &   UROTC  = "  UROT  ", VROTC  = "  VROT  ", WROTC  = "  WROT  ",
     &   STRAINC= " STRAIN "

      CHARACTER(LEN=8), PARAMETER ::
     &   F1RC   = "   F1R  ", F2RC   = "   F2R  ", F3RC   = "   F3R  ",
     &   F1MC   = "   F1M  ", F2MC   = "   F2M  ", F3MC   = "   F3M  ",
     &   F1NC   = "   F1N  ", F2NC   = "   F2N  ", F3NC   = "   F3N  ",
     &   F1WC   = "   F1W  ", F2WC   = "   F2W  ", F3WC   = "   F3W  ",
     &   F1EC   = "   F1E  ", F2EC   = "   F2E  ", F3EC   = "   F3E  ",
     &   F1KC   = "   F1K  ", F2KC   = "   F2K  ", F3KC   = "   F3K  ",
     &   F1EPC  = "  F1EPS ", F2EPC  = "  F2EPS ", F3EPC  = "  F3EPS ",
     &   F1HC   = "   F1H  ", F2HC   = "   F2H  ", F3HC   = "   F3H  ",
     &   F1RNUC = "  F1RNU ", F2RNUC = "  F2RNU ", F3RNUC = "  F3RNU ",
     &   F1G    = "   F1G  ", F2G    = "   F2G  ", F3G    = "   F3G  ",
     &   F1RET  = "  F1RET ", F2RET  = "  F2RET ", F3RET  = "  F3RET " 

      CHARACTER(LEN=8), PARAMETER ::
     &   D1C    = "   D1   ", D2C    = "   D2   ", D3C    = "   D3   "
      CHARACTER(LEN=8), PARAMETER ::
     &   XCC    = "  XC    ", YCC    = "  YC    ", ZCC    = "  ZC    ",
     &   XCOC   = "  XCO   ", YCOC   = "  YCO   ", ZCOC   = "  ZCO   ",
     &   XORC   = "  XGRI  ", YORC   = "  YGRI  ", ZORC   = "  ZGRI  ",
     &   XSKC   = "  XSKEW ", YSKC   = "  YSKEW ", ZSKC   = "  ZSKEW "
      CHARACTER(LEN=8), PARAMETER ::
     &   ROFORC = "  ROFOR ", EFORC  = "   EFOR ", PDFORC = "  PDFOR ",
     &   RMFORC = "  RMFOR ", RNFORC = "  RNFOR ", RWFORC = "  RWFOR ",
     &   RKFORC = "  RKFOR ", REFORC = "  REFOR ", GFORC  = "  GFOR  ",
     &   RETFORC= "  RETFOR"

      CHARACTER(LEN=8), DIMENSION(6),PARAMETER ::
     &   UUC   = (/ "    UU  ","    UV  ","    UW  ",
     &              "    VV  ","    VW  ","    WW  " /),
     &   SRC   = (/ "  S_XX  ","  S_XY  ","  S_XZ  ",
     &              "  S_YY  ","  S_YZ  ","  S_ZZ  " /),
     &   S11C  = (/ " SHEARXX"," SHEARXY"," SHEARXZ",
     &              " SHEARYY"," SHEARYZ"," SHEARZZ" /)

! ... SCALAR EQ. ARE LIMITED TO 10. PPR 23.3

      CHARACTER(LEN=8), DIMENSION(10),PARAMETER ::
     &   FIC   = (/
     &   "  FI_1  ","  FI_2  ","  FI_3  ","  FI_4  ","  FI_5  ",
     &   "  FI_6  ","  FI_7  ","  FI_8  ","  FI_9  ","  FI10  " /),
     &   DFIC  = (/
     &   "  DFI_1 ","  DFI_2 ","  DFI_3 ","  DFI_4 ","  DFI_5 ",
     &   "  DFI_6 ","  DFI_7 ","  DFI_8 ","  DFI_9 ","  DFI10 " /),
     &   F1FIC = (/
     &   "  F1FI_1","  F1FI_2","  F1FI_3","  F1FI_4","  F1FI_5",
     &   "  F1FI_6","  F1FI_7","  F1FI_8","  F1FI_9","  F1FI10" /),
     &   F2FIC = (/
     &   "  F2FI_1","  F2FI_2","  F2FI_3","  F2FI_4","  F2FI_5",
     &   "  F2FI_6","  F2FI_7","  F2FI_8","  F2FI_9","  F2FI10" /),
     &   F3FIC = (/
     &   "  F3FI_1","  F3FI_2","  F3FI_3","  F3FI_4","  F3FI_5",
     &   "  F3FI_6","  F3FI_7","  F3FI_8","  F3FI_9","  F3FI10" /),
     &   FOLDC = (/
     &   " FIOLD_1"," FIOLD_2"," FIOLD_3"," FIOLD_4"," FIOLD_5",
     &   " FIOLD_6"," FIOLD_7"," FIOLD_8"," FIOLD_9"," FIOLD10" /)

      CHARACTER(LEN=8), DIMENSION(11) :: CHAR_LIQ,CHAR_GAS ! Mersu

      CHARACTER(LEN=8), DIMENSION(20,3),PARAMETER :: ! Mersu
     &   CHAR_PH  = RESHAPE( (/
     &   "   ROL  ","  TEMPL "," DROLDP "," DROLDH ","  VISL  ",
     &   "   CHL  ","   CL   ","  QIFL  "," HLSAT  "," DHLSDP ",
     &   "   CPL  ","   EL   "," TAUFL  ","  EQL   "," DROLDH ",
     &   "  HLTOT ","   EL   "," TAUFL  "," DROLDP "," DROLDH ",
     &   "   ROG  ","  TEMPG "," DROGDP "," DROGDH ","  VISG  ",
     &   "   CHG  ","   CG   ","  QIFG  "," HGSAT  "," DHGSDP ",
     &   "   CPG  ","   EG   "," TAUFG  ","  EQG   "," DROGDH ",
     &   "  HGTOT ","   EG   "," TAUFG  "," DROGDP "," DROGDH ",
     &   "   ROS  ","  TEMPS "," DROSDP "," DROSDH ","  VISS  ",
     &   "   CHS  ","   CS   ","  QIFS  "," DHSSDP "," DHSSDH ",
     &   "   CPS  ","   ES   "," TAUFS  "," DROSDP "," DROSDH ",
     &   "   CPS  ","   ES   "," TAUFS  "," DROSDP "," DROSDH " /),
     &                       (/ 20,3 /) )

      CHARACTER(LEN=8), DIMENSION(35,3),PARAMETER :: ! Mersu
     &   CHAR_VAR = RESHAPE( (/
     &   " ALFAL  ","   XL   "," EVAPL  "," DTEMPL ","  F1RL  ",
     &   "  F2RL  ","  F3RL  ","  F1EL  ","  F2EL  ","  F3EL  ",
     &   "   DXL  ","   DEL  ","  DROL  ","   DHL  ","   ROL  ",
     &   "  AROL  "," AROLEL "," DAROL  ","DAROLHL ","   ROL  ",
     &   "  F1ML  ","  F2ML  ","  F3ML  ","  F1NL  ","  F2NL  ",
     &   "  F3NL  ","  F1WL  ","  F2WL  ","  F3WL  ","   DML  ",
     &   "   DNL  ","   DWL  ","   UL   ","   VL   ","   WL   ",
     &   " ALFAG  ","   XG   "," EVAPG  "," DTEMPG ","  F1RG  ",
     &   "  F2RG  ","  F3RG  ","  F1EG  ","  F2EG  ","  F3EG  ",
     &   "   DXG  ","   DEG  ","  DROG  ","   DHG  ","   ROG  ",
     &   "  AROG  "," AROGEG "," DAROG  ","DAROGHG ","   ROL  ",
     &   "  F1MG  ","  F2MG  ","  F3MG  ","  F1NG  ","  F2NG  ",
     &   "  F3NG  ","  F1WG  ","  F2WG  ","  F3WG  ","   DMG  ",
     &   "   DNG  ","   DWG  ","   UG   ","   VG   ","   WG   ",
     &   " ALFA1  ","   X1   ","  MELT? "," DTEMPS ","  F1R3  ",
     &   "  F2R3  ","  F3R3  ","  F1E3  ","  F2E3  ","  F3E3  ",
     &   "   DES  ","  F1E3  ","  F2E3  ","  F3E3  ","   DES  ",
     &   "   DES  ","  F1E3  ","  F2E3  ","  F3E3  ","   DES  ",
     &   "   DES  ","  F1E3  ","  F2E3  ","  F3E3  ","   DES  ",
     &   "   DES  ","  F1E3  ","  F2E3  ","  F3E3  ","   DES  ",
     &   "   DES  ","  F1E3  ","  F2E3  ","  F3E3  ","   DES  " /),
     &                       (/ 35,3 /) )

      END MODULE CHARACTERS


      MODULE NS3CO

      IMPLICIT NONE

      INTEGER, PARAMETER :: IN = 2, JN = 2, KN = 2 ! Space for ghost cells and
      INTEGER, PARAMETER :: LN = MAX(IN,JN,KN) + 0 ! for BC patch 'ghost cells'

      INTEGER, PARAMETER :: NBGGG = 5000, NRESI = 67, IC9 = 28
      INTEGER, PARAMETER :: NOVARMAX = 100

      INTEGER, PARAMETER :: MAXPRO = 1000
      
! ... NBGGG is Tthe maximum number of blocks in all processes for RSTART
! ... NRESI is Tthe number of convergence residuals
! ... IC9 is the maximum number of ICON parameters per PATCH + 1
! ... NOVARMAX is the maximum number of output variables for post-processing  
! ... MAXPRO is the maximum number of processes

      REAL ::
     &  T,THETA,H0,CD,CL,DROMAX,CX,CY,CZ,CMX,CMY,CMZ, 
     &  RMUFRS,FRSMUT,TOM,FRSVIS,FRSSIE,T0,ANGLE,CS,  
     &  VOLTOT,QMFIN,QMEIN,QMFOUT,QMEOUT,RKMAX,TAVER_OLD ! Mersu  

      REAL, DIMENSION(3) :: AXIS
      REAL, DIMENSION(NBGGG) :: ROTANG

      INTEGER, DIMENSION(NOVARMAX) :: OUTVAR
                                      
      INTEGER ::
     &  ICYCLE,IPRINT,IXERR,INERR,IMACH,NBLOCK,
     &  ICYOLD,IUPDAT,NSCAL,JSCAL,IBOUT,ITIMES,ICYTOT,NAVER

      REAL ::
     &  CFL,CFLL,DROLIM,TMAX,DT,RMACH,ALPHA,BETA,RE,PR,PRT,RGAS,
     &  TSU0,VISU0,EXPSU,GAMMA,FRSTEM,FRSVEL,FRSDEN,EPSLIM,RKLIM,
     &  TURLIM,TURBLE,RMUINI,CMGK,CMGEPS,CSIMPS,CC1,CC2,AREF,CHLREF,
     &  GRILEN,REFPRE,DIFPRE,PSEUCO,GROUND,
     &  FRSPRE,TOLER,TEMINI,SIEINI,TLOLIM,TUPLIM,FRSX,FRSALF,SMAX,
     &  CAVLEV,CDIFF,CDIFFT,TROT,ALTITUDE

      REAL :: GX,GY,GZ,G0  ! Matti Palin 04.08.2020

      REAL, DIMENSION(NBGGG) :: ROTAT

      INTEGER ::
     &  IOLD,LEVEL,ITURB,KSCAL,IFLUX,ICMAX,KP,MPRINT,ISTATE,
     &  JRDIF,JRDIS,JRPRE,JRIMP,IPRESC,IDRXX,LUSGS,IEPSMA,
     &  ITERMA,NOB,NBLOCG,NCHIM,KOVER,IDIS,NPHASE,
     &  MCYCLE,ITERHA,MTOT

! ... COMMON INPUC IS USED IN MPI TO SEND COMMON INPUT TO SLAVES
! ... this is used with MPI data change but not written in RSTART:

      INTEGER :: IROTCO,KRRINT,NSPTOT,ISTRES,NSKIN,INTERTU,ICONV3

      LOGICAL ::
     &  STARTL,STRESL,FULLNS,SOURL,STATEL,TIMEL,COORL,PRESL,
     &  REANEW,GRAVIL,CONVL,SPLIT,CHIMEL,GROUP,PERCHL,MULPHL,
     &  PARALLEL,TURCOR,XXTRAL,FRESUL,WOODYL,LDISTL,CAVITL,TRUE_DISTL,
     &  FPRINTL,INOUTL,TURDESL,INCHIML,REFLECL,ENTROPY_FIX,
     &  SURFACE_FORCES,TUR_FRS_SOURCE,TRANSL,TUR_MULPHL,NEGVL,
     &  LONG_STARTL,LONG_TRANSL,AVER_STARTL,TWO_FLUIDL,START_CAVL,
     &  PLOT3D, ENSIGHT, USE_QUATERNIONS_LOGICAL, WALLFUNL,
     &  MPCCIL, MODALFSIL, AUTOCONL

      REAL ::
     &  DTOLD,RTIME,OMEX,OMEY,OMEZ,CENX,CENY,CENZ,GVEX,GVEY,GVEZ,
     &  FRSSSP,FRSLEN,FRSRK,FRSEPS,XMASSB(4),E0REF,T0REF,RTIME1,
     &  QMXIN,QMXOUT,QMYIN,QMYOUT,QMZIN,QMZOUT,CMXIN,CMXOUT,
     &  CMYIN,CMYOUT,CMZIN,CMZOUT,QVFIN,QVFOUT

      REAL :: SOLTEM,ARTSSP,RMULTV,ALFAP,ALFAU,RJK2,RJK4,TFRSSP

! ... These are for the free-surface model (AGAMMA & ABANK) are 
!     for the trajectory calculation (ROTROW) is for the helicopter
!     calculation

      INTEGER :: JFIRST, NFSD, ICFST, INWH, IFSBC, INDXTS, IGLOBL,
     &           NGLOBL, IPRIFS, FEMXYZ, RMESHT
                                                                 
      REAL ::
     &  FROUDE,DTWMAX,DWMV,DWMAX,SUMDWH,WHMAX,WHMIN,
     &  FREDIF,FLODWH,CFLFRE,DTWMIN,AGAMMA,ABANK,ADVANJ,REFVEL,ROTORW,
     &  WHEIGHT,XBULB,YBULB,ABULB,QGFIN,QGFOUT,QGEIN,QGEOUT,NEGV,
     &  QLFIN,QLFOUT,QLEIN,QLEOUT

      REAL :: RO0,FLOW,POUT ! Unknown parameters from the ancient history

      CHARACTER(3) :: PRN

      REAL :: WTIMEN,WCAL,WCOM,WTIMST,WTBCST,WTBCEN,TAVER,
     &  XMOM,YMOM,ZMOM

      END MODULE NS3CO

! -------

      MODULE CONSTANTS

      IMPLICIT NONE

      REAL, PARAMETER :: EPS   = 1.E-20,EPS10 = 1.E-10,EPS20 = 1.E-20,
     &                   EPS6  = 1.E-6, EPS8  = 1.E-8, EPS12 = 1.E-12,
     &                   EPS14 = 1.E-14

      REAL, PARAMETER :: PII = ACOS(-1.0)

      REAL, PARAMETER :: DEG2RAD = 4.0 *ATAN(1.0)/180.0
      REAL, PARAMETER :: RAD2DEG = 45.0/ATAN(1.0)  ! Matti Palin 12.06.2020

      REAL :: PR,PRT,VISU0,EXPSU,TSU0,E0REF,T0REF
      REAL :: CB1,CB2,SRNU,CV1,CW1,CW2,CW3

      END MODULE CONSTANTS

! -------

      MODULE TYPE_ARRAYS

      INTEGER,PARAMETER :: NPHASES = 2               ! Mersu
      INTEGER,PARAMETER :: LEN_BLOCKS = 1 + 2 + 4*10 ! Useless?
   
      TYPE PROPERTIES
       REAL,DIMENSION(NPHASES) :: TEMP,RO,E,DRODP,DRODH,VIS,CH,C,QIF,
     &                            HSAT,DHSDP,DHSDH,CP,TAUF,ALPO
       REAL :: TSAT,PSAT,DPSAT,DPSDT,DTSDP
       REAL, DIMENSION(NPHASES) :: DTEMP
       REAL :: SIGMA,EO
       REAL,DIMENSION(3) :: FD,FL,FTD,FVM,FW,FITOT
      END TYPE PROPERTIES

      TYPE MPHASE_VARIABLES

       REAL,DIMENSION(NPHASES) :: ALFA,X,DX,XOLD,DE,DTEMP,TOLD,EVAPR,
     &                            F1R,F2R,F3R,FRO,FE,FA,F1A,F2A,F3A,
     &                            EVAPH,ARO,ARE,AROLE2,AROLE3,ARELE2,
     &                            ARELE3,DH,ETOT,EQL,ADDE,APSAT,
     &                            FM,DM,DN,DW,U,V,W,UOLD,VOLD,WOLD,
     &                            FRM,FRN,FRW,F1X,F2X,F3X,FX,
     &                            ARM,ARN,ARW,ARMLE2,ARNLE2,ARWLE2,
     &                            ARMLE3,ARNLE3,ARWLE3,BUBDIA,
     &                            F1U,F2U,F3U,FU,DU,DV,DWW

       REAL,DIMENSION(3) :: RKIF
       REAL, DIMENSION(NPHASES) :: DTOLD

      END TYPE MPHASE_VARIABLES

      TYPE MPHASE_FORTIFY
       REAL, DIMENSION(NPHASES) :: ALFAFOR,XFOR
       REAL, DIMENSION(NPHASES) :: DTEMPFOR
      END TYPE MPHASE_FORTIFY

      TYPE MGRID_VARIABLES
       REAL,DIMENSION(NPHASES) :: ROP2H,EP2H,RMP2H,RNP2H,RWP2H
      END TYPE MGRID_VARIABLES

      TYPE BLOCKS
       INTEGER :: NPHASE,IFLUX,IPRESC,LUSGS,EVAP_TYPE,IDIFF,MPCASE,
     & PREVEL,INLRC,OUTRC,ICERO,ITERMP,INTERI,INTERJ,INTERK,TURMUL,
     & DRAGMUL,LIFTMUL,DISPERMUL,VMASS,WFORCEMUL
       INTEGER,DIMENSION(NPHASES) :: ISTATE,ICHAR,IUPTEM
       CHARACTER(10) :: SOLUTION_TYPE
       CHARACTER(10),DIMENSION(NPHASES) :: MATERIAL
       REAL :: RKLIM,TURLIM,EPSLIM,PLOLIM,PUPLIM,TLOLIM,TUPLIM,RMUFRS,
     &         FRSMUT,FRSVIS,FRSSIE,T0,VOLTOT,QMFIN,QMEIN,FRSTEM,FRSVEL,
     &         FRSDEN,FRSPRE,TEMINI,SIEINI,FRSSSP,FRSLEN,FRSRK,FRSEPS,
     &         ARTSSP,SOLTEM,FRSTUR,CHLREF,RE,RMACH,PR,PRT,FRSCP,FRSCH,
     &         REFPRE,DIFPRE,RMUINI,AREF,RMULTV,CAVNO,REFVEL,ALFAP,
     &         ALFAU,REFCAV,CHLCAV,SKEWI,SKEWJ,SKEWK,SKEWIMIN,SKEWJMIN,
     &         SKEWKMIN,ALFASKEW,QGFIN,QGEIN,QGFOUT,QGEOUT,FRSG,FRSRET,
     &         VOIDTT,TFRSSP,GVEX,GVEY,GVEZ,FREDIFI,FREDIFJ,FREDIFK,
     &         RINTF
       REAL,DIMENSION(NPHASES) :: FRSALFA,FRSX,TAUF,N,RK,FRADEN
       REAL :: DFRSPRE,DFRSTEM,DFRSDEN
       LOGICAL :: CONVL,SKEWCORR,COMPCORR,KATOL,ZEROV,FLUXCORR
      END TYPE BLOCKS

      TYPE PRE_COR
       REAL :: F1R,F2R,F3R,FIR,FJR,FKR,DRO,AP,DP,DPDX,DPDY,DPDZ,
     &         FPM,FPN,FPW
      END TYPE PRE_COR

      TYPE DIFF_COR
       REAL :: S1,S2,S3,S1X,S1Y,S1Z,S2X,S2Y,S2Z,S3X,S3Y,S3Z
      END TYPE DIFF_COR

      TYPE INTERMITTENCY
       REAL :: G,RET,GOLD,RETOLD,DG,DRET,FG,FRET,
     &         SG,ZG,LG,SRET,ZRET,LRET,GSEP,GEFF,GFOR,RETFOR
      END TYPE INTERMITTENCY

      TYPE MGRID_INTERMIT
       REAL :: PG,PRET
      END TYPE MGRID_INTERMIT

      TYPE SIXDOF
       REAL :: RMASS,AREF,CHLREF,RCG,IX,IY,IZ,IXY,IXZ,IYZ,DAMPN,DAMPT,
     &         TOL,H1,HMIN,TDEL,TSTEER,SPEED,ALPHA,BETA,GAMMA,BANK,
     &         VX,VY,VZ,ROTX,ROTY,ROTZ,CX,CY,CZ,CMX,CMY,CMZ,
     &         FXE,FYE,FZE,MXE,MYE,MZE,FXT,FYT,FZT,MXT,MYT,MZT,
     &         FXH,FYH,FZH,MXH,MYH,MZH
       CHARACTER(1) :: CHAR
      END TYPE SIXDOF

      TYPE FORCE_GROUP
       CHARACTER(LEN=1)  :: FG_SHORT_NAME
       CHARACTER(LEN=24) :: FG_FULL_NAME
       REAL              :: XMOM_REF, YMOM_REF, ZMOM_REF
       REAL              :: XMOM_AXS, YMOM_AXS, ZMOM_AXS
      END TYPE FORCE_GROUP

      END MODULE TYPE_ARRAYS

! -------

      MODULE MAIN_ARRAYS

      USE TYPE_ARRAYS

      IMPLICIT NONE

      INTEGER,ALLOCATABLE,DIMENSION(:) ::
     &        IMAXG, JMAXG, KMAXG, IDOB,  JDOB,  KDOB,
     &        IDI1,  IDI2,  IDI3,  INTERI,INTERJ,INTERK,
     &        MIB,   MJB,   MKB,   MIT,   MJT,   MKT,
     &        KX1S,  KY1S,  KZ1S,  KX2S,  KY2S,  KZ2S,
     &        IMINM, JMINM, KMINM, IMAXM, JMAXM, KMAXM,
     &        ISTRS, JSTRS, KSTRS, IDIMS, JDIMS, KDIMS,
     &        IDER,  INITC, IROTVE,MGRID, MGRIDA,LAMIN,
     &        MOV,   MOVPO, NCHIMT,NCHOR, NCPAT, NLOCAL,NORMAL,
     &        NPATCH,NPNUM, NSOLPA,NSPB,  NSPG,  NSPP,  NSSB,
     &        IDIMSG,JDIMSG,IUPPT, OSCLS

      INTEGER,ALLOCATABLE,DIMENSION(:) :: ICNH,ICGH,ICONH,ICOGH,
     &        IBTGR, ITPGR, JBTGR, JTPGR, KBTGR, KTPGR,
     &        IBOTGR,ITOPGR,JBOTGR,JTOPGR,KBOTGR,KTOPGR,
     &        IHULL, JHULL, LHULL, MHULL, MMHUL

      REAL,ALLOCATABLE,DIMENSION(:) :: OMEGA, OMEGAX, OMEGAY, OMEGAZ,
     &     CENAX, CENAY, CENAZ, AMPL, 
     &     RLOLIM, UPPLIM, PSIGSC, PSIGS2, VOLN,
     &     BXMIN, BXMAX, BYMIN, BYMAX, BZMIN, BZMAX

      REAL,ALLOCATABLE,DIMENSION(:) :: TOMEGA,XCOG,APATCH,
     &     CXB,CYB,CZB, CMXB,CMYB,CMZB, DXB,DYB,DZB,
     &     BUX,BUY,BUZ, VMXB,VMYB,VMZB, QTB,QWB,QHB

      REAL,ALLOCATABLE,DIMENSION(:) ::
     &     XC1,YC1,ZC1, XC2,YC2,ZC2, XC3,YC3,ZC3, XC4,YC4,ZC4

      REAL,ALLOCATABLE,DIMENSION(:) :: WH, WFS, XI, ETA, ZZTOP,
     &     XXI,  YXI,  ZXI,  XHULL, YHULL, ZHULL, XVL, YVL, ZVL,
     &     XETA, YETA, ZETA, UHULL, VHULL, WHULL, XSP, YSP, ZSP

      REAL,ALLOCATABLE,DIMENSION(:) :: RCON
      REAL,ALLOCATABLE,DIMENSION(:,:) :: AAA,RESI,XCOR
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: XCOL

      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: NPROCE, ISSB, IT,IL,IK,
     &        NTOT, IMAX,JMAX,KMAX, IBOT,ITOP,JBOT,JTOP,KBOT,KTOP,
     &        JF,IG,IH,IR,IQ, ITAG,IHF,IC, JAPP,IWAPP, JSTATE, IGRID

      INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: IW

      INTEGER,ALLOCATABLE,DIMENSION(:) :: ICON

      INTEGER,ALLOCATABLE,DIMENSION(:) :: LOCDIS, IDP,
     &     ICP,JCP,KCP,IJMASK,JET

      REAL,ALLOCATABLE,DIMENSION(:) :: VOL, DISTW, DISTP,
     &     RLD, RLDIST,
     &     D1,D2,D3, UROT,VROT,WROT, XFC,YFC,ZFC,UROTCP,VROTCP,WROTCP

      REAL,ALLOCATABLE,DIMENSION(:) ::
     &     A1,A1XA,A1YA,A1ZA, A2,A2XA,A2YA,A2ZA, A3,A3XA,A3YA,A3ZA

      REAL,ALLOCATABLE,DIMENSION(:) :: BLANK, BLANK2,
     &     RO,   RM,   RN,   RW,   E,   C,
     &     DRO,  DM,   DN,   DW,   DE,  TEMP,
     &     TOLD, UOLD, VOLD, WOLD, POLD,SUROLD,
     &     PDIFF,RUVAV,RUWAV,RVWAV,EPS2,RKSI,
     &     OHMI, U,    V,    W,    P,   DTL,
     &     CP,   CH,   DRDP, DRDH, VIS, VIST,
     &     STRAIN, ROLD, TIJ

      REAL,ALLOCATABLE,DIMENSION(:) :: DDEPS, PTUR, FUN1, TTS,
     &     RK, REPS, DRK, DEPS, RKOLD, EPSOLD, SRK, SEPS,VTRAN,
     &     RNUT, QSAS, VELLAP

      REAL,ALLOCATABLE,DIMENSION(:) :: HAT1, HAT2, HAT3, HAT4,
     &     F1R, F1RM, F1RN, F1RW, F1E, F1RK, F1EPS, F1H

! ... Surface-related arrays, size = IB, number = 61

      REAL,ALLOCATABLE,DIMENSION(:) :: HFLUX, CPWALL,
     &     UWALL, VWALL, WWALL, TWALL, QWALL, QWFRIC, SURLE, DSURLE,
     &     SURFX, SURFY, SURFZ, POROS, WMFLUX,WHSTAG, WTEMP, RBK,
     &     RSDIRX,RSDIRY,RSDIRZ,TAUWX, TAUWY, TAUWZ,  TAUW1,TAUW2,
     &     BOUNU, BOUNV, BOUNW, BOUNT, BOUNP, BOUNR,  BOUNE,
     &     BOUNRK,BOUNEP,BOUNPD,BOUNA1,BOUNA2,BOUNMF, WTRAN,RMLOSS,
     &     BOUNG,BOUNRET,SURFA, SURFT ,SURF2X,SURF2Y,SURF2Z,SURFMX,
     &     SURFMY,SURFMZ,BOUNU1,BOUNU2,BOUNV1,BOUNV2,BOUNW1,BOUNW2,
     &     SURFPX,SURFPY,SURFPZ,UTAUM

! ... Surface-related arrays for output, size = IB, number = 28
      REAL,ALLOCATABLE,DIMENSION(:) :: CPWALL_, QWALL_, QWFRIC_,
     &     TWALL_, TAUW1_, TAUW2_, TAUWX_, TAUWY_, TAUWZ_,
     &     SURFX_, SURFY_, SURFZ_, HFLUX_, SURLE_, DSURLE_,
     &     F1RK_,  WTRAN_, WMFLUX_,
     &     BOUNR_, BOUNT_, BOUNP_, BOUNRK_, BOUNEP_, BOUNMF_,
     &     BOUNU_, BOUNV_, BOUNW_, BOUNE_

      
      REAL,ALLOCATABLE,DIMENSION(:) :: WAVES,WAVEH,FTR2

      INTEGER,ALLOCATABLE,DIMENSION(:) :: IWAVEB

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NEARBLOCK

! ... Arrays for time-accurate simulation, size = MAXTI, number = 28

      REAL,ALLOCATABLE,DIMENSION(:) ::
     &    ROAV1, RMAV1, RNAV1, RWAV1, EAV1, RKAV1, EPSAV1,
     &    ROAV2, RMAV2, RNAV2, RWAV2, EAV2, RKAV2, EPSAV2,
     &    ROLE2, RMLE2, RNLE2, RWLE2, ELE2, RKLE2, EPSLE2,
     &    ROLE3, RMLE3, RNLE3, RWLE3, ELE3, RKLE3, EPSLE3,
     &    RMAV3, RNAV3, RWAV3, PLE2,  PLE3,  PAV1,   TAV1

      REAL,ALLOCATABLE,DIMENSION(:) :: W12, SIJ, GRADT, GRADK, GREPS

      REAL,ALLOCATABLE,DIMENSION(:,:) :: BIJ, WIR, S11 , VORT, SHEAR


      REAL,ALLOCATABLE,DIMENSION(:,:) :: BOUNFI, BOUNBI

      REAL,ALLOCATABLE,DIMENSION(:) :: ROFOR, RMFOR, RNFOR, RWFOR,
     &                                 PDFOR,  EFOR, RKFOR, REFOR

      REAL,ALLOCATABLE,DIMENSION(:) :: ROP2H, RMP2H, RNP2H, RWP2H,
     &                                 EPSP2H,  EP2H, RKP2H, SURH2

      REAL,ALLOCATABLE,DIMENSION(:) :: PROD, SPI, DIF, DIS, VVIS, FWLL

      REAL,ALLOCATABLE,DIMENSION(:) :: ! Obsolete arrays
     &    UBI, VBI, WBI, UBJ, VBJ, WBJ, UBK, VBK, WBK,
     &    UTI, VTI, WTI, UTJ, VTJ, WTJ, UTK, VTK, WTK

      REAL,ALLOCATABLE,DIMENSION(:,:) :: FIP2H, FIFOR, FI, FIOLD,
     &     DFI, SFI, F1FI, FIAV1, FIAV2, FILE2, FILE3

      REAL,ALLOCATABLE,DIMENSION(:) :: WGH1, WGH2, WGH3, WGH4

      REAL,ALLOCATABLE,DIMENSION(:) :: XGRI,YGRI,ZGRI,
     &   XC,YC,ZC, XCO,YCO,ZCO, XLE2,YLE2,ZLE2, XLE3,YLE3,ZLE3,
     &   XORI,YORI,ZORI,XCP,YCP,ZCP

!      CHARACTER(80),ALLOCATABLE,DIMENSION(:) :: BOUNDF

! ... Temporary space for Reynolds stressis, contravariant vectors and
!     grid velocity vectors, nondimensiolized turbulent viscosity,
!     relative velocity components, relative momentum components, relative
!     internal energy, relative momentum components in multiphase flows and
!     relative velocity components of different phases in multi phase
!     simulations needed during output. To be allocated in OUTP.

      REAL,ALLOCATABLE,DIMENSION(:) :: RS1, RS2, RS3, RS4, RS5, RS6
      REAL,ALLOCATABLE,DIMENSION(:) :: CVV1, CVV2, CVV3, CVV4, CVV5,
     &                                 CVV6, CVV7, CVV8, CVV9
      REAL,ALLOCATABLE,DIMENSION(:) :: GVX, GVY, GVZ
      REAL,ALLOCATABLE,DIMENSION(:) :: VISTND
      REAL,ALLOCATABLE,DIMENSION(:) :: UREL, VREL, WREL
      REAL,ALLOCATABLE,DIMENSION(:) :: RMREL, RNREL, RWREL, EREL
      REAL,ALLOCATABLE,DIMENSION(:,:) :: MPUREL, MPVREL, MPWREL

! ... To be allocated in NS3D:

      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBSKIN, ITSKIN, IPSKIN
      INTEGER,ALLOCATABLE,DIMENSION(:) :: JLOC, JTRA
      INTEGER,ALLOCATABLE,DIMENSION(:) :: II1, II2, II3, II4, INTP

      REAL,ALLOCATABLE,DIMENSION(:) :: APP, ZZZ

      REAL,ALLOCATABLE,DIMENSION(:,:) :: SKINX, SKINY, SKINZ,
     &   SKINXN, SKINYN, SKINZN

! ... Derived type arrays (allocated in MAIN)

      TYPE(PROPERTIES),      ALLOCATABLE,DIMENSION(:) :: PRO
      TYPE(MPHASE_VARIABLES),ALLOCATABLE,DIMENSION(:) :: VAR
      TYPE(MGRID_VARIABLES), ALLOCATABLE,DIMENSION(:) :: P2H
      TYPE(BLOCKS),          ALLOCATABLE,DIMENSION(:) :: BLKS
      TYPE(PRE_COR),         ALLOCATABLE,DIMENSION(:) :: PRC
      TYPE(DIFF_COR),        ALLOCATABLE,DIMENSION(:) :: SDI
      TYPE(INTERMITTENCY),   ALLOCATABLE,DIMENSION(:) :: TRM
      TYPE(MGRID_INTERMIT),  ALLOCATABLE,DIMENSION(:) :: P2HTRM
      TYPE(MPHASE_FORTIFY),  ALLOCATABLE,DIMENSION(:) :: MPFOR


      END MODULE MAIN_ARRAYS

! -------

      MODULE INTEGERS

      INTEGER :: NMOV,MGMIN,IXMA,IYMA,IZMA,MAXW,MAW,
     &           NBCS,IPRO,IOTY,NCON,NUMCOL
      INTEGER :: MAXB,MAX11,IB,IBF,MAXSB,MAXS2,MAXSS,MAXTI,MAXTS,MAXEB,
     &     NB,NBGG,NPPV,MGM,MBPRO,NPRO,MAXCH,MAXMP,MAXPC,MAXHB,MAXFS,
     &     MAXTR
      INTEGER :: ITERAD,ITERAC,MCYCAM
      INTEGER,DIMENSION(20) :: IREPEA,NREPEA
      INTEGER(KIND=8) :: IP, MREQ, NPTIM, MBITS

      END MODULE INTEGERS

! -------

      MODULE FLIGHT

      USE TYPE_ARRAYS

      INTEGER :: NGRIFL,IPATH

      INTEGER,DIMENSION(19) :: TRMODE

      INTEGER,DIMENSION(19) :: IFA,IFT,IFR,NSERIES,IGRSLAVE,
     &     NBLADE,              ! Actuator disc
     &     IFAKEP               ! Fake propeller
      INTEGER,DIMENSION(19,6) :: ISCALE_ACT

      REAL,DIMENSION(19) :: XCGI,YCGI,ZCGI,XCG,YCG,ZCG,
     &     PSIRI,THETARI,PHIRI,PSIR,THETAR,PHIR,PSIM  ! A/C and helicopter data

      REAL,DIMENSION(19) :: DRAUGHTI,DRAUGHT,SINKI,SINK,
     &     TRIMAI,TRIMA,DAMPC1,DAMPC2,XCGIS,YCGIS,ZCGIS,RGML,
     &     XFAKEP,YFAKEP,ZFAKEP,ROTB1FAKEP,
     &     XFAKEPNEW,YFAKEPNEW,ZFAKEPNEW  ! Ship data

      REAL,DIMENSION(19) :: RIN,ROUT,THRUST,TORQUE,ADV,
     &     ROTA1,ROTB1,ROTA1I,ROTB1I,CONEA,CONEAI,ACTUA,VTIP,CDBLADE,
     &     CFBLADE,SIGMA,CFTI,ROTA1S,ROTB1S,RTMSP,FDSP,UTSP,FXSP,FXTSP,
     &     QFACT ! Actuator disc

      TYPE(SIXDOF),DIMENSION(19) :: OSKU

! ... Particle properties

      REAL :: M, S, C1, B, IX, IY, IZ, IXY, IXZ, IYZ, DAMPN, DAMPT

! ... Particle state at t = 0

      REAL :: XEA, YEA, ZEA, FIIPAR, THETAPAR, PSIIPAR, 
     &     UPAR, VPAR, WPAR, PPAR, QPAR, RPAR

! ... Path integration control parameters

      REAL :: TOL, TOL0, H1, HMIN
      REAL :: XEND, TDEL, TNEXT, XT, TSTEER, VLEKO
      REAL :: AKONE, BKONE

! ... Time derivatives

      REAL :: UP, VP, WP, PP, QP, RP, FIIP, THETAP, PSIIP, XP, YP, ZP

! ... Quaternion components for a specific particle

      REAL :: particleQuaternion0, particleQuaternion1
      REAL :: particleQuaternion2, particleQuaternion3

! ... Logical variables

      LOGICAL :: FLYOBJ,MVSHIP,SHIP,SHIPWAVE,ACTDISK,HROTOR,HMVBLADE,
     &     SHIPPROP

      END MODULE FLIGHT

! -------

      MODULE BLADE_CONSTANTS

      REAL, PARAMETER :: PIJ = 4.0*ATAN(1.0)
      REAL, PARAMETER :: DEG2RAD = 4.0 *ATAN(1.0)/180.0
      REAL, PARAMETER :: RAD2DEG = 45.0/ATAN(1.0)

      END MODULE BLADE_CONSTANTS

! -------

      MODULE BLADE_VARIABLES

      USE BLADE_CONSTANTS

      INTEGER :: ITERM,IORD,IPRINT,IPR,LLKM,ICYCLE,ICYCLEB(4),IPRB(4)

      REAL,DIMENSION(19) :: FHINGEX,FHINGEY,FHINGEZ,
     &     FHINGEI,FHINGEJ,FHINGEK,RMHINGEI,RMHINGEJ,RMHINGEK

      REAL,DIMENSION(19) :: QZE,QBE,QTH,QZE_A,QTH_A,
     &     QBE_A,QZE_S1,QBE_S1,QZE_D,QBE_D,QZE_S2,QBE_S2

      REAL               :: OMEGA,PSI,SHIFT

      REAL,DIMENSION(19) :: ZETAH,BETAH,DOTZE,
     &     DOTBETA,THETA,DOTTHETA,P,Q,R,SHAFT,ETIP,CBLADE

      REAL,DIMENSION(19) :: RM,RITH,RIBE,RIZE,
     &     RCG,RHINGE,RD,RKZE,RKBE,RKTH,RKD,CD,RCZE,RCBE

      REAL               :: Y(4,19),YOLD(4,19),F(4,19)

      REAL,DIMENSION(19) :: THETACOL,THCYCLON,THCYCLAT,KDELTA3,
     &     BETAHCOL,BECYCLON,BECYCLAT,ZETAHCOL

      REAL               :: DTB,TMAXB,TDBL,TIME, ! TDBL = TIME ?
     &     TIMEB(19),DISSZE(19),DISSBE(19)

      CHARACTER(LEN=1) :: CLAPNR

      END MODULE BLADE_VARIABLES

! -------

      MODULE PATH

      INTEGER :: KMAX,KOUNT

      REAL :: DXSAV,XXP(200),YYP(12,200)

      END MODULE PATH

! -------

      MODULE MPCCIVARS

C ... MpCCI

      INTEGER :: fnelems, fnnodes, transfered, strongcoupling

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: felemnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: fnodechimt
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: fnodecoor
      REAL, ALLOCATABLE, DIMENSION(:,:) :: fnodecoor_orig
      REAL, ALLOCATABLE, DIMENSION(:)   :: fnodedisp
      REAL, ALLOCATABLE, DIMENSION(:,:) :: felemforce
      REAL, ALLOCATABLE, DIMENSION(:)   :: felempres
      REAL, ALLOCATABLE, DIMENSION(:)   :: felemtemp
      REAL, ALLOCATABLE, DIMENSION(:)   :: felemarea

C ... Modal FSI
      
      INTEGER :: idimjoints, jdimjoints, njoints, nfreq
      INTEGER :: activefreqs

      CHARACTER(LEN=80) :: JOINTS_FILE, MODE_SHAPES_FILE

      REAL, ALLOCATABLE, DIMENSION(:,:) :: jointxyz
      REAL, ALLOCATABLE, DIMENSION(:,:) :: jointxyzo
      REAL, ALLOCATABLE, DIMENSION(:,:) :: jointdisp
      REAL, ALLOCATABLE, DIMENSION(:,:) :: jointxyzn
      REAL, ALLOCATABLE, DIMENSION(:,:) :: jointf
      REAL, ALLOCATABLE, DIMENSION(:,:) :: shapexyz
      REAL, ALLOCATABLE, DIMENSION(:)   :: shapedydx
      REAL, ALLOCATABLE, DIMENSION(:)   :: shapedxdz
      REAL, DIMENSION(100)   :: freq, fomega, zeta
      REAL, DIMENSION(100,3) :: eta, etadot, etaddot, fm, fmx

************************************************************************
      REAL, ALLOCATABLE, DIMENSION(:,:) :: shapefelem
************************************************************************
      
      END MODULE MPCCIVARS

! -------

      MODULE FORCE_GROUPS

      USE TYPE_ARRAYS

      TYPE(FORCE_GROUP), DIMENSION(52) :: FORCE_GROUP_INFO

      REAL, DIMENSION(52) :: XMOMR, YMOMR, ZMOMR
      REAL, DIMENSION(52) :: XMOMA, YMOMA, ZMOMA

      CHARACTER(LEN=52*1)  :: FORCE_GROUP_SHORT_NAME
      CHARACTER(LEN=52*24) :: FORCE_GROUP_FULL_NAME

      END MODULE FORCE_GROUPS

! -------
