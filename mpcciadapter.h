/****************************************************************************
 *   mpcciadapter.h                                                         *
 *                                                                          *
 *   MpCCI Code API  -  Code Adapter Header File                            *
 *   Changed for use with FINFLO                                            *
 *                                                                          *
 *   This file contains declarations of all driver functions.               *
 *   See the MpCCI Programmers Guide for more information.                  *
 *                                                                          *
 *   http://www.scai.fraunhofer.de/mpcci.html                               *
 *                                                                          *
 ****************************************************************************/
#ifndef __MPCCI_ADAPTER_HEADER_INCLUDED__
#define __MPCCI_ADAPTER_HEADER_INCLUDED__

/* use the mpcci code API */
#include "mpcci.h"

/*#########################################################################
 *############  INTERFACE FUNCTIONS                         ###############
 *#########################################################################*/

/* initialize the coupling procedure */
void initcoupling();

/* transfer quantities */
void dotransfer();

/* exit the coupling procedure */
void exitcoupling();

/*#########################################################################
 *############  DRIVER FUNCTIONS                            ###############
 *#########################################################################*/

/* Update MPCCI_PART and /or update MPCCI_QUANT */
int MpCCI_Driver_partUpdate(MPCCI_PART *part, MPCCI_QUANT *quant);

/****************************************************************************
 *       COUPLING COMPONENT DEFINITIONS
 ****************************************************************************/

/* REQUIRED: Method 1: Define the grid */
int MpCCI_Driver_definePart(MPCCI_SERVER *server, MPCCI_PART *part);

/****************************************************************************
 *       DATA EXCHANGE
 ****************************************************************************/

/* get data from faces - called before send */
int MpCCI_Driver_getFaceElemValues(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);

/* get global values - called before send */
int MpCCI_Driver_getGlobValues(const MPCCI_GLOB *glob, void *values);

/* put data from nodes of a face - called after receive */
void MpCCI_Driver_putFaceNodeValues(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);

/* put global values - called after receive */
void MpCCI_Driver_putGlobValues(const MPCCI_GLOB *glob, void *values);

#endif
