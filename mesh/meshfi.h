
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
 * meshfi.h - file handling routines for mesh
 *
 * revision log :
 *
 * 27.07.1995	ggu	routines written from scratch
 *
\************************************************************************/


#ifndef __GUH_MESHFI_
#define __GUH_MESHFI_

#include "nlist.h"

void ReadFiles( int argc , char *argv[] );
void ReadBnd( char *fname );
void WriteAll( char *file , NodeList list );
void WriteGrd( char *file );
void WriteFile( char *file , NodeList list , int all );


#endif


