
/************************************************************************\
 *
 *    Copyright (C) 1995,1997  Georg Umgiesser
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
 * color.h - general header file for color routines
 *
 * revision log :
 *
 * 11.02.1995	ggu	routines written from scratch
 * 10.10.1997	ggu	new prototype for QAllocHSBColors()
 *
\************************************************************************/


#ifndef __GUH_COLOR_
#define __GUH_COLOR_


int *QAllocGreen2BlueColors( int ncol );
int *QAllocYellow2GreenColors( int ncol );
int *QAllocRed2YellowColors( int ncol );
int *QAllocRed2BlueColors( int ncol );
int *QAllocBlueColors( int ncol );
int *QAllocHSBColors( int ncol );


#endif /* __GUH_COLOR_ */
