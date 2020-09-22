
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995  Georg Umgiesser
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
 * xgraph.h - header file for X11 routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 * 21.03.1994	ggu	declared *ShadeColor() in this file
 * 05.12.1995	ggu	changes (Widget, Cursor, ...) included
 *
\************************************************************************/


#ifndef __GUH_XGRAPH_
#define __GUH_XGRAPH_

#if __GUG_UNIX_

#include "generalx.h"
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>

#endif /* __GUG_UNIX_ */


/*
void QNextEvent( XEvent *eventp );
*/

Display *QActDisplay( void );
Window *QActWindow( void );
Window *QRootWindow( void );
int QActScreen( void );
GC QActGc( void );
void QSetTest( void );
void QClearTest( void );
void QFlush( void );
void QSync( void );
void QSyncron( int );
void QSetDisplay( char *s );
void QSetWidget( Widget w );
void QSetTitle( char *s );
void QSetGeometry( char *s );

#endif /* __GUH_GRAPH_ */
