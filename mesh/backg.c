
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
 * backg.c - background grid interpolation routines                     *
 *									*
 * Revision History:							*
 * 17-Nov-97: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdlib.h>
/*
#include <stdio.h>
#include <math.h>
*/


#include "general.h"
#include "grd.h"
#include "grdhs.h"
#include "backg.h"


Backg_type *SetBackGrid( Grid_type *G , Elem_type *pe )

{
	int i,in,ip;
	float f=0.;
	float x[3],y[3];
	Node_type *pn;
	Backg_type *pb;

	pb = (Backg_type *) malloc( sizeof(Backg_type) );
	if( !pb )
		Error("SetBackGrid: Cannot allocate background mesh");

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(G->HN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
	}

	for(i=0;i<3;i++) {
		in = (i+1)%3;
		ip = (i+2)%3;
		pb->a[i] = x[in]*y[ip] - x[ip]*y[in];
		pb->b[i] = y[in] - y[ip];
		pb->c[i] = x[ip] - x[in];
		f += pb->a[i];
	}
	pb->fr = 1./f;

	return pb;
}

float InterpolBackGrid( Backg_type *pb , float x , float y , float *value )

{
	int i;
	float f;
	float a = 0.;

	for(i=0;i<3;i++) {
		f = pb->a[i] + x*pb->b[i] + y*pb->c[i];
		a += value[i] * f * pb->fr;
	}

	return a;
}

