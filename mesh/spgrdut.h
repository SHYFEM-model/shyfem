
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


#include "general.h"

#include "nlist.h"

#include "grd.h"
/*
#include "grdhs.h"
#include "grdut.h"
*/



void MakeCounterClockwise( void );
void MakeCM( Elem_type *pe, float *xm, float *ym );
int InConvex( int n , float *xe , float *ye , float x , float y );
float AreaElement( Hashtable_type H , Elem_type *pe );
float AreaConvex( NodeList list );
float LengthConvex( NodeList list );
double angle( float x1, float y1, float x2, float y2, float x3, float y3 );
int InClosedLine( int ndim , float *xl , float *yl , float x , float y );
int TurnClosedLine( int ndim , float *xl , float *yl );
void MakeCoordsFromLine( Line_type *pl , float **x , float **y );
