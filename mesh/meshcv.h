
/************************************************************************\ 
 *									*
 * meshcv.h - routines for convex hull to be used with mesh		*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshcv.c for copying information					*
 *									*
 * Revision History:							*
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_MESHCV_
#define __GUH_MESHCV_


#include "nlist.h"
#include "heap.h"

void ConvexHull( NodeList hull, NodeList intern );
void Convex( NodeList hull, NodeList intern );
void getcoord( HeapItem *item , float *x , float *y );
float theta( float x1, float y1, float x2, float y2 );
int ccw( float x1, float y1, float x2, float y2, float x3, float y3 );
void sortondist( HeapTable HL, int imin, int imax, float dk
                        , float xmin, float ymin );

#endif

