
/************************************************************************\
 *
 *    Copyright (C) 1997  Georg Umgiesser
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
 *
 * backg.h - background grid interpolation routines
 *
 * revision log :
 *
 * 17.11.1997	ggu	routines written from scratch
 *
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

