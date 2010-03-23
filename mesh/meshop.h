
/************************************************************************\ 
 *									*
 * meshop.h - command line options for mesh                             *
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshop.c for copying information					*
 *									*
 * Revision History:							*
 * 02-Aug-95: routines copied from gridop.c                             *
 *									*
\************************************************************************/


#ifndef __GUH_MESHOP_
#define __GUH_MESHOP_


/**************************************************************/
/**************************************************************/
/********************* global variables ***********************/
/**************************************************************/
/**************************************************************/

extern int OpArgc;

extern int OpBackGround;
extern int OpIntern;
extern int OpBoundary;
extern int OpRefineBoundary;
extern int OpInsertIntern;
extern int OpCompress;
extern int OpSmoothPass;
extern float OpSmoothOmega;
extern float OpAspect;
extern float OpResolution;


void SetOptions(int argc, char *argv[]);
void TestVersion( void );


#endif

