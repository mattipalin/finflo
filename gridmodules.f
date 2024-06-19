       MODULE XYZK
       REAL :: XDT(200),YDT(200),ZDT(200),KT(200)
       INTEGER :: IBASIC
       END MODULE XYZK
       
       MODULE TZANSC
       LOGICAL :: virgin = .true.
       INTEGER :: ktrnsm,ibtrns,ittrns,ktrd
       END MODULE TZANSC
       
       MODULE SKIP
       INTEGER :: iskip,ifix,lime
       END MODULE SKIP
       
       MODULE STRIM
       INTEGER :: isinkt,icfull,ic0
       REAL :: relax0
       END MODULE STRIM
       
       MODULE TRSINK
       REAL :: osink,otrim,oxtrim,oytrim,oztrim,
     $          sink,trimm,xtrim,ytrim,ztrim
       END MODULE TRSINK
       
       MODULE GLIM
       INTEGER :: nlim,iblim(10),itlim(10),jlim(10),limf(10),
     $             ijklim(8,10),ia,ib
       REAL :: xlim(0:200,0:200,10),ylim(0:200,0:200,10),
     $          zlim(0:200,0:200,10),dlim(0:200)
       END MODULE GLIM

       MODULE VLCB
       REAL :: v0,alcb0,alcf0,awp0
       END MODULE VLCB

       MODULE MATU
       REAL :: S0,CZshi0,CMshi0
       END MODULE MATU

       MODULE DDDAMP
       REAL :: damp0,ddamp
       END MODULE DDDAMP

       MODULE CENG
       REAL :: xcg,ycg,zcg
       END MODULE CENG

       MODULE ACHTNG
       INTEGER :: iacht
       END MODULE ACHTNG

       MODULE CGLOC
       REAL :: zcgloc(2)
       END MODULE CGLOC

       MODULE GM
       REAL :: GML
       END MODULE GM

c remove the following and store CZ or CMY to COMPUT: the jold is actually IOLD
c and is only needed to get right initial values for sink and trim when starting
c from COMPUT - I'm confused now but I'll rethink this sometime in the future
c perhaps - something like this anyhow
       MODULE VANHA
       INTEGER :: jold
       END MODULE VANHA
       
