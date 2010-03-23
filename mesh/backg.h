
/************************************************************************\ 
 *									*
 * backg.h - background grid interpolation routines                     *
 *									*
 * Copyright (c) 1997 by Georg Umgiesser				*
 *									*
 * see backg.c for copying information					*
 *									*
 * Revision History:							*
 * 17-Nov-97: routines written from scratch				*
 *									*
\************************************************************************/



#ifndef __GUH_BACKG_
#define __GUH_BACKG_


#include "grd.h"


typedef struct {
	float a[3];
	float b[3];
	float c[3];
	float fr;
} Backg_type;


Backg_type *SetBackGrid( Grid_type *G , Elem_type *pe );
float InterpolBackGrid( Backg_type *pb , float x , float y , float *value );



#endif

