
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * grdio.h - read/write grd files
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshgd
 *
\************************************************************************/


#ifndef __GUH_GRDIO_
#define __GUH_GRDIO_


#include "hash.h"
#include "queue.h"
#include "grd.h"


void  ReadStandard( char *fname , Grid_type *G , char *ext );
void WriteStandard( char *fname , Grid_type *G );


#endif


