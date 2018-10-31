
/************************************************************************\
 *
 *    Copyright (C) 1985-2018  Georg Umgiesser
 *
 *    This file is part of SHYFEM.
 *
 *    SHYFEM is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SHYFEM is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SHYFEM. Please see the file COPYING in the main directory.
 *    If not, see <http://www.gnu.org/licenses/>.
 *
 *    Contributions to this file can be found below in the revision log.
 *
\************************************************************************/


/************************************************************************\ 
 *									*
 * meshcv.h - routines for convex hull to be used with mesh		*
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

