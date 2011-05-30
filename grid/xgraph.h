
/************************************************************************\ 
 *									*
 * xgraph.h - header file for X11 routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see xgraph.c for copying information					*
 *									*
 * Revision History:							*
 * 05-Dec-95: changes (Widget, Cursor, ...) included                    *
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
