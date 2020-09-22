
/************************************************************************\
 *
 *    Copyright (C) 1994  Georg Umgiesser
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
 * debug.h - debugging routines
 *
 * revision log :
 *
 * 31.12.1994	ggu	routines written from scratch
 *
\************************************************************************/


#ifndef __GUH_DEBUG_
#define __GUH_DEBUG_


/**************************************************************************/

void GDPrint( void );
void GDInit( int size );
char *GDAlloc( int size );
void GDNL( void );
void GDS( char *s );
void GDI( int i );
void GDL( long l );
void GDF( float f );

/**************************************************************************/

#endif
