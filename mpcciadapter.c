/****************************************************************************
 *   mpcciadapter.c                                                         *
 *                                                                          *
 *   MpCCI Code API  -  Code Adapter Source File Changed for use with       *
 *   FINFLO flow solver written in FORTRAN                                  *
 *                                                                          *
 *   This file contains definitions of all required driver functions.       *
 *   See the MpCCI Programmers Guide for more information.                  *
 *                                                                          *
 *   http://www.scai.fraunhofer.de/mpcci.html                               *
 *                                                                          *
 ****************************************************************************/

/* includes: */

#include <stdio.h>    /* standard C libs */
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "mpcciadapter.h"   /* code adapter header file */

/*#########################################################################
 *############  DATA STRUCTURES                             ###############
 *#########################################################################*/

/* list of driver functions
 *****************************************************/
static MPCCI_DRIVER MpCCIDriverFunctions = {
  /******* REQUIRED: information to do a compatibility check  *********/
  (unsigned)sizeof(MPCCI_DRIVER),     /* size of structure MPCCI_DRIVER */
  MPCCI_CCM_VERSION,                  /* this_version */
  0,                                  /* tact_required bitmask */
  /* CHANGED: Name of simulation code and code-specific names */
  {                   /* list of code-specific names (for output) */
     "single point - UNUSED", /* name for point */
     "Line - UNUSED",         /* name for lines */
     "Solid surface",         /* name for faces */
     "Volume - UNUSED",       /* name for volumes */
     "Global value"           /* name for global quantities */
  },

  /********* REQUIRED: (unsigned)sizeof(float|double): tell the manager
             the size of floating point values that are arguments
             to the set of putXXXValues() functions  *******************/
  (unsigned)sizeof(double),                     /* values_size */


  /********** OPTIONAL: methods called before/after some actions
              Replace by NULL pointers if not needed. **********/
  NULL,                            /* MpCCI_Driver_afterCloseSetup   */
  NULL,                            /* MpCCI_Driver_beforeGetAndSend  */
  NULL,                            /* MpCCI_Driver_afterGetAndSend   */
  NULL,                            /* MpCCI_Driver_beforeRecvAndPut  */
  NULL,                            /* MpCCI_Driver_afterRecvAndPut   */

  NULL,                            /* MpCCI_Driver_partSelect        */
  MpCCI_Driver_partUpdate,         /* MpCCI_Driver_partUpdate        */
  NULL,                            /* MpCCI_Driver_appendPart        */
  NULL,                            /* MpCCI_Driver_getPartRemeshState*/

  /********* REQUIRED: methods to update and check component names and
                       to define the grid (nodes & elements) *********/
  MpCCI_Driver_definePart,         /* MpCCI_Driver_definePart        */

  NULL,                            /* MpCCI_Driver_partInfo          */
  NULL,                            /* MpCCI_Driver_getNodes          */
  NULL,                            /* MpCCI_Driver_getElems          */

  /********* REQUIRED for moving nodes instead of full remeshing:
                      returns the MpCCI error code *******************/
  NULL,                            /* MpCCI_Driver_moveNodes         */

  /********* OPTIONAL: Data exchange functions ***********/
  NULL,                            /* MpCCI_Driver_getPointValues    */
  NULL,                            /* MpCCI_Driver_getLineNodeValues */
  NULL,                            /* MpCCI_Driver_getLineElemValues */
  NULL,                            /* MpCCI_Driver_getFaceNodeValues */
  MpCCI_Driver_getFaceElemValues,  /* MpCCI_Driver_getFaceElemValues */
  NULL,                            /* MpCCI_Driver_getVoluNodeValues */
  NULL,                            /* MpCCI_Driver_getVoluElemValues */
  MpCCI_Driver_getGlobValues,      /* MpCCI_Driver_getGlobValues     */

  NULL,                            /* MpCCI_Driver_putPointValues    */
  NULL,                            /* MpCCI_Driver_putLineNodeValues */
  NULL,                            /* MpCCI_Driver_putLineElemValues */
  MpCCI_Driver_putFaceNodeValues,  /* MpCCI_Driver_putFaceNodeValues */
  NULL,                            /* MpCCI_Driver_putFaceElemValues */
  NULL,                            /* MpCCI_Driver_putVoluNodeValues */
  NULL,                            /* MpCCI_Driver_putVoluElemValues */
  MpCCI_Driver_putGlobValues,      /* MpCCI_Driver_putGlobValues     */
};


static MPCCI_TINFO mpcciTinfo = {0};
static MPCCI_JOB  *mpcciJob   = NULL;



/*#########################################################################
 *############  ADDITONAL HELPER FUNCTIONS                  ###############
 *#########################################################################*/

/* Declaration of FORTRAN subroutines from mpccitools.f */
extern void adapteroutput_(const char*,int);
extern void error_(const char*,int);
extern void getsurfaceid_(int*, char*, int);
extern void getinfo_(int*,int*,int*,int*);
extern void getmesh_(int*,int*,double*,int*);
extern void getpressure_(int*, double*);
extern void getwallforce_(int*, double*);
extern void getwalltemp_(int*, double*);
extern void putnposition_(int*, double*);
extern void getphysicaltime_(double*);
extern void gettimestepsize_(double*);
extern void gettimestepcount_(double*);
extern void getreferencepressure_(double*);
extern void putphysicaltime_(double*);
extern void puttimestepsize_(double*);
extern void puttimestepcount_(double*);
extern void putreferencepressure_(double*);



/*#########################################################################
 *############  INTERFACE FUNCTIONS                         ###############
 *#########################################################################*/

/* initialize the coupling process */
/* CHANGED: Appended underscore to function name */
void initcoupling_(double *t, double *dt, int *strongCoupling){

  MPCCI_CINFO cinfo;
   int ret = 0;
   bool transient = true;
   int iterativeCoupling = 0;

   /* already initialized? */
   if (mpcciTinfo.mpcci_state > 0) return;
     
   printf("Initializing the coupling process!\n");
     
   /* initialize output routines */
   /* CHANGED: Names of output routines */
   umpcci_msg_functs(NULL,            /* printFullMessage */
                     NULL,            /* printFullWarning */
                     NULL,            /* printFullFatal   */
                     adapteroutput_,  /* printLineMessage */
                     adapteroutput_,  /* printLineWarning */
                     error_,          /* printLineFatal   */
                     NULL);           /* atExitHandler    */


   umpcci_msg_prefix("MpCCI: ");

   /* call ampcci_tinfo_init */
   if (!ampcci_tinfo_init(&mpcciTinfo,NULL))
   {
      mpcciTinfo.mpcci_state = -1;
      return;
   }

  /*
   * initialize the transfer/time info and the code info
   */

   transient = (MPCCI_TSOLU_WANT_STEADY(&mpcciTinfo)) ? false : true;
   iterativeCoupling = MPCCI_TINFO_WANT_IMPL(&mpcciTinfo);
   *strongCoupling = iterativeCoupling;

   if(transient)
     {
       if (iterativeCoupling)
	 {
	   //	   MPCCI_MSG_FATAL1("FINFLO MpCCI adapter ERROR: Transient iterative coupling not implemented %i.\n",1);

	   // transient implicit
	   mpcciTinfo.iter = 0;
	   mpcciTinfo.time = *t;
	   mpcciTinfo.dt   = *dt;
	 }
       else
	 {
	   // transient explicit
	   mpcciTinfo.iter = -1;
	   mpcciTinfo.time = *t;
	   mpcciTinfo.dt   = *dt;
	 }
     }
   else
     {
       // steady-state
       mpcciTinfo.iter = 0;
       mpcciTinfo.time = -1.0;
       mpcciTinfo.dt   = -1.0;
     }

   memset(&cinfo,0,sizeof(cinfo));
   mpcci_cinfo_init(&cinfo,&mpcciTinfo);

   /* CHANGED: Name of simulation code and error handling */
   cinfo.codename = "FINFLO";
   cinfo.flags    = MPCCI_CFLAG_TYPE_FEA|MPCCI_CFLAG_GRID_CURR;
   cinfo.time     = mpcciTinfo.time;
   cinfo.iter     = mpcciTinfo.iter;
   cinfo.nclients = 1;
   cinfo.nprocs   = 1;
   /* call mpcci_init */
   mpcciJob       = mpcci_init(NULL,&cinfo);
   

      if (mpcciJob)
   {
      ret = ampcci_config(&mpcciJob,&MpCCIDriverFunctions);
      if (ret < 0)
         mpcci_quit(&mpcciJob);
   }


   if (!mpcciJob)
   {
      mpcciTinfo.mpcci_state = -1;
      if (ret < 0)
         MPCCI_MSG_FATAL0("This is no coupled computation!");
   }
   else
   {
      mpcciTinfo.mpcci_state = 1;
      mpcciTinfo.mpcci_used  = 1;
   }
}

/* perform a transfer of quantities */
/* CHANGED: Appended underscore to function name */
void dotransfer_(double *t, double *dt, int *ret2, int *strongCoupling){

   int ret;
   bool transient = true;
   int iterativeCoupling = 0;
   
   if (!mpcciTinfo.mpcci_state)
/* CHANGED: Appended underscore to function name */
     initcoupling_(t,dt,strongCoupling);


   if (!mpcciTinfo.mpcci_used)
   {
      MPCCI_MSG_WARNING0("MpCCI is not used.\n");
      return;
   }

   /* Update the time /iteraton information for mpcciTinfo */
   /* foundation code is just a QTID code */


   transient = (MPCCI_TSOLU_WANT_STEADY(&mpcciTinfo)) ? false : true;
   iterativeCoupling = MPCCI_TINFO_WANT_IMPL(&mpcciTinfo);

   if(transient)
     {
       if (iterativeCoupling)
	 {
	   //	   MPCCI_MSG_FATAL1("FINFLO MpCCI adapter ERROR: Transient iterative coupling not implemented %i.\n",1);

	   // transient implicit
	   mpcciTinfo.iter++;
	   mpcciTinfo.time = *t;
	   mpcciTinfo.dt   = *dt;
	 }
       else
	 {
	   // transient explicit
	   mpcciTinfo.iter = -1;
	   mpcciTinfo.time = *t;
	   mpcciTinfo.dt   = *dt;
	 }
     }
   else
     {
       // steady-state
       mpcciTinfo.iter++;
       mpcciTinfo.time = -1.0;
       mpcciTinfo.dt   = -1.0;
     }

   /* do the transfer */
   ret = ampcci_transfer(mpcciJob,&mpcciTinfo);
   *ret2 = ret;

   if(ret < 0)
       MPCCI_MSG_WARNING0("Could not complete all requested transfers!");
}

/* exit the coupling process */
/* CHANGED: Appended underscore to function name */
void exitcoupling_(){
   mpcci_quit(&mpcciJob);
}

/*###########################################################################
 *############  COUPLING COMPONENT DEFINITIONS                ###############
 *#########################################################################*/

/****************************************************************************
 *  MpCCI_Driver_partUpdate
 *
 *  This function is called during ampcci_config.
 *  You should set the number of nodes and number of elements
 *  of all coupling components here.
 *
 ****************************************************************************/
/****************************************************************************/
int MpCCI_Driver_partUpdate(MPCCI_PART *part, MPCCI_QUANT *quant)
/****************************************************************************/
{
   int surfaceID;         /* ID of a coupling surface */

   /* CHANGED: Added three variables for call of getinfo below */
   int nnodes, nelems, spacedim;

   /* CHANGED: call getsurfaceid and getinfo */
   /* get internal ID of the coupling component */
   getsurfaceid_(&surfaceID,MPCCI_PART_NAME(part),strlen(MPCCI_PART_NAME(part)));
   MPCCI_PART_PARTID(part) = surfaceID;
   /* get coupling component information */
   getinfo_(&surfaceID,&nnodes,&nelems,&spacedim);
   /* set number of nodes */
   MPCCI_PART_NNODES(part) = nnodes;
   /* set number of elements */
   MPCCI_PART_NELEMS(part) = nelems;


   /* Define the coordinates system and Model dimension:
      Update this to your current model */
   MPCCI_PART_CSYS(part) = MPCCI_CSYS_C3D;
   MPCCI_PART_MDIM(part) = MPCCI_MDIM_FACE;

   return 0;
}

/***************************************************************************
 * MpCCI_Driver_definePart
 *
 * This function get called during MpCCI_Init.
 *
 * MpCCI needs information about the used mesh for a coupled component.
 * You must define nodes and elements within this driver function.
 *
 ****************************************************************************/
/*****************************************************************************/
int MpCCI_Driver_definePart(MPCCI_SERVER *server, MPCCI_PART *part)
/*****************************************************************************/
{
   /* IMPORTANT: The node definitions may look completely different in your
                 code. This is just an example! */

   double *nodeCoords;      /* node coordinates */
   int surfaceID;           /* surface ID */
   int maxNumNodes;         /* maximum number of nodes */
   int i;                   /* loop counter */
   int nnodes, nelems;      /* numbers of nodes and elements */
   int dim;                 /* space dimension */
   int * nodeNumbers;       /* array of node IDs */
   int * elemNodes;         /* array of element nodes */
   unsigned * elemTypes;    /* array of element types */
   int ret;

   surfaceID = MPCCI_PART_PARTID(part);
   /* get coupling component information */
   getinfo_(& surfaceID,& nnodes,& nelems,& dim);
   switch(dim)
   {
      case 2:
         maxNumNodes = 2;
         break;
      case 3:
         maxNumNodes = 4;
         break;
      default:
        MPCCI_MSG_FATAL0("Wrong space dimension, must be 2 or 3.");
   }
   /***********************
   STEP 1: Define the part
   ***********************/
   MPCCI_MSG_INFO1
   (
      "Coupling grid definition for component \"%s\" ...\n",
      MPCCI_PART_NAME(part)
   );

   ret = smpcci_defp
         (
            server,
            MPCCI_PART_MESHID(part),
            MPCCI_PART_PARTID(part),
            MPCCI_PART_CSYS  (part),
            MPCCI_PART_NNODES(part),
            MPCCI_PART_NELEMS(part),
            MPCCI_PART_NAME  (part)
         );

   /***********************
   STEP 2: Collect nodes of the surface in one array
   ***********************/
   if(!ret)
   {
      /* allocate memory for the node coordinates (3-dimensional)*/
      nodeCoords = (double*) malloc(nnodes*sizeof(double)*dim);

      /* array for node-numbers */
      nodeNumbers = (int*) malloc(nnodes*sizeof(int));

      /* allocate memory for the list of nodes for each element */
      elemNodes   = (int*) malloc(nelems*maxNumNodes*sizeof(int));

      getmesh_(& surfaceID, nodeNumbers, nodeCoords, elemNodes);

      /* Send the nodes definiton for this coupled component */
      ret = smpcci_pnod
            (
               server,
               MPCCI_PART_MESHID(part),  /* mesh id */
               MPCCI_PART_PARTID(part),  /* part id */
               MPCCI_PART_CSYS  (part),  /* coordinates system */
               MPCCI_PART_NNODES(part),  /* no. of nodes */
               nodeCoords,               /* node coordinates */
               (unsigned)sizeof(double), /* data type of coordinates */
               nodeNumbers,              /* node ids */
               NULL
            );
      free(nodeCoords);
      free(nodeNumbers);
   }

   if (!ret)
   {
      elemTypes   = (unsigned*) malloc(nelems*sizeof(unsigned));
      /* loop over elements in coupling component */
      for(i=0; i<nelems; i++)
      {
         /* get element ID and element information */
         switch(dim){
            case 2:  /* -> line elements */
               elemTypes[i] = MPCCI_ETYP_LINE2;
               break;
            case 3:  /* -> quadrilaterals */
               elemTypes[i] = MPCCI_ETYP_QUAD4;
               break;
            default:
               MPCCI_MSG_FATAL0("Unsupported element type!");
               break;
          }
      }

      /* Send the element definiton for this coupled component */
      ret = smpcci_pels(
                    server,
                    MPCCI_PART_MESHID(part),    /* mesh id */
                    MPCCI_PART_PARTID(part),    /* part id */
                    MPCCI_PART_NELEMS(part),    /* number of elements */
                    elemTypes[0],               /* first element type */
                    elemTypes,                  /* element types */
                    elemNodes,                  /* nodes of the element */
                    NULL);                      /* element IDs */
   }

   free(elemTypes);
   free(elemNodes);
   return ret;
}

/*###########################################################################
 *############  DATA EXCHANGE                                 ###############
 *#########################################################################*/

/*****************************************************************************/
int MpCCI_Driver_getFaceElemValues( const MPCCI_PART *part,
                                    const MPCCI_QUANT *quant,
                                    void *values)
/*****************************************************************************/
{
   int surfaceID;                       /* internal ID of coupling component */

   MPCCI_MSG_INFO0("entered send values...\n");

   /* check whether this is the appropriate method */
   MPCCI_MSG_ASSERT(MPCCI_PART_IS_FACE(part));

   /* get internal ID of the coupling component */
   surfaceID = MPCCI_PART_PARTID(part);

   switch (MPCCI_QUANT_SMETHOD(quant))
   {
      case MPCCI_QSM_DIRECT: /* direct load/store values */

         /* number of nodes */
         /* nnodes = MPCCI_PART_NNODES(part);            */

         /* distinguish between the different quantities */
         switch (MPCCI_QUANT_QID(quant))
         {
            case MPCCI_QID_WALLFORCE:
               getwallforce_(& surfaceID, values); break;
            case MPCCI_QID_OVERPRESSURE:
               getpressure_(& surfaceID, values);  break;
            case MPCCI_QID_WALLTEMPERATURE:
               getwalltemp_(& surfaceID, values);  break;
            default:
               MPCCI_MSG_FATAL0("Quantity not supported\n"); break;
         }
         break;

      default:
         MPCCI_MSG_FATAL0("Unsupported storage method.");
         break;
   }
   MPCCI_MSG_INFO0("finished send values...\n");

   return sizeof(double); /* return the size of the value data type */
}


/*****************************************************************************/
void MpCCI_Driver_putFaceNodeValues( const MPCCI_PART *part,
                                     const MPCCI_QUANT *quant,
                                     void *values)
/*****************************************************************************/
{
   int surfaceID;                      /* internal ID of coupling component */
   /* CHANGED: valptr not needed */

   MPCCI_MSG_INFO0("entered receive values...\n");

   /* check whether this is the appropriate method */

   MPCCI_MSG_ASSERT(MPCCI_PART_IS_FACE(part));

   /* get internal ID of the coupling component */
   surfaceID = MPCCI_PART_PARTID(part);

   switch(MPCCI_QUANT_SMETHOD(quant))
   {
      case MPCCI_QSM_DIRECT: /* direct load/store values */

         /* distinguish between the different quantities */
         switch (MPCCI_QUANT_QID(quant))
         {
            case MPCCI_QID_NPOSITION:
	       putnposition_(& surfaceID, values);
               break;
            default:
               MPCCI_MSG_FATAL0("Quantity not supported\n");
               break;
         }
         break;

      default:
         MPCCI_MSG_FATAL0("Unsupported storage method.");
         break;

   }
  MPCCI_MSG_INFO0("finished receive values...\n");
}
//************************************************************************************************
int MpCCI_Driver_getGlobValues(const MPCCI_GLOB *glob, void *values)
//************************************************************************************************
{
    MPCCI_MSG_INFO1
    (
       "SEND global quantity value '%s'.\n",
       glob->name
    );

    switch(MPCCI_QUANT_QID(glob))
    {
      case MPCCI_QID_PHYSICAL_TIME:  getphysicaltime_(values);      break;
      case MPCCI_QID_TIMESTEP_SIZE:  gettimestepsize_(values);      break;
      case MPCCI_QID_TIMESTEP_COUNT: gettimestepcount_(values);     break;
      case MPCCI_QID_REF_PRESSURE:   getreferencepressure_(values); break;

    default: MPCCI_MSG_FATAL1("Global quantity '%s' not yet supported!\n", MPCCI_QUANT_NAME(glob));
    }
    // return sizeof(scalar);
    return sizeof(double);
}
//************************************************************************************************
void MpCCI_Driver_putGlobValues(const MPCCI_GLOB *glob, void *values)
//************************************************************************************************
{
    MPCCI_MSG_INFO1
    (
       "RECEIVE global quantity value '%s'.\n",
       glob->name
    );

    switch(MPCCI_QUANT_QID(glob))
    {
      case MPCCI_QID_PHYSICAL_TIME:  putphysicaltime_(values);      break;
      case MPCCI_QID_TIMESTEP_SIZE:  puttimestepsize_(values);      break;
      case MPCCI_QID_TIMESTEP_COUNT: puttimestepcount_(values);     break;
      case MPCCI_QID_REF_PRESSURE:   putreferencepressure_(values); break;

    default: MPCCI_MSG_FATAL1("Global quantity '%s' not yet supported!\n", MPCCI_QUANT_NAME(glob)); break;
    }
}
