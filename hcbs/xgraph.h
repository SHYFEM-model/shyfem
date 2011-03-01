
/************************************************************************\ 
 *									*
 * xgraph.h - header file for X11 routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see xgraph.c for copying information					*
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
