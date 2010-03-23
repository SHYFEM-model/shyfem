
/************************************************************************\ 
 *									*
 * grdio.h - read/write grd files					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see grdio.c for copying information					*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines copied from meshgd                               *
 *									*
\************************************************************************/


#ifndef __GUH_GRDIO_
#define __GUH_GRDIO_


#include "hash.h"
#include "queue.h"
#include "grd.h"


void  ReadStandard( char *fname , Grid_type *G , char *ext );
void WriteStandard( char *fname , Grid_type *G );


#endif


