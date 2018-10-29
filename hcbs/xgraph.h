
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
 * xgraph.h - header file for X11 routines				*
 *									*
 * Revision History:							*
 * 21-Mar-94: declared *ShadeColor() in this file			*
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
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

void QNextEvent( XEvent *eventp );
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
int QAllocNewColor      (
                        long unsigned int red ,
                        long unsigned int green ,
                        long unsigned int blue
                        );

int BlueShadeColor( int shade );
int RedShadeColor( int shade );
int YellowShadeColor( int shade );
int GreenShadeColor( int shade );
int VelShadeColor( int shade );

#endif /* __GUG_UNIX_ */
#endif /* __GUH_GRAPH_ */
