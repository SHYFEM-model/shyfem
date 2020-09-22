
/************************************************************************\
 *
 *    Copyright (C) 1992,1994  Georg Umgiesser
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
 * graph.h - general header file for graphic routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 *
\************************************************************************/


#ifndef __GUH_GRAPH_
#define __GUH_GRAPH_

#include "general.h"

/* definitions for colors */

extern int Black;
extern int Blue;
extern int Green;
extern int Cyan;
extern int Red;
extern int Magenta;
extern int Brown;
extern int LightGray;
extern int DarkGray;
extern int LightBlue;
extern int LightGreen;
extern int LightCyan;
extern int LightRed;
extern int LightMagenta;
extern int Yellow;
extern int White;



void QGraphInit( void );
void QClearScreen( void );
void QGraphClose( void );
void QNewPen( int color );
void QGetPen( int *color );
void QGetPixelWidth( float *dx , float *dy );
void QViewport( int left , int top , int right , int bottom );
void QClearViewport( void );
void QWindow( float xmin , float ymin , float xmax , float ymax );
void QAdjustScale( void );
void QGetMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax );
void QSetMinMaxPix( int xmin , int ymin , int xmax , int ymax );
void QLine( float x1 , float y1 , float x2 , float y2 );
void QPlot( float x , float y );
void QMove( float x , float y );
void QPoint( float x , float y );
void QAreaFill( int ndim , float *x , float *y );
void QRectFill( float x1 , float y1 , float x2 , float y2 , int color );
void QText( float x , float y , char *s );
void QTextSize( int size );
void QTextDimensions( char *s , float *width , float *height );
void QRealXY( int vx , int vy , float *x , float *y );
void QScreenXY( float x , float y , int *vx , int *vy );
void QFlush( void );
void QSync( void );
void QTest( void );
void QSyncron( int sync );

#if __GUG_DOS_

#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+0.5)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom-ViewTop+0.5)
#define QxVtoW(x) (SxVtoW*(x)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom+ViewTop)+WinYmin)

#else
#if __GUG_UNIX_

#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+ViewLeft+0.5)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom+0.5)
#define QxVtoW(x) (SxVtoW*((x)-ViewLeft)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom)+WinYmin)

#endif /* __GUG_UNIX_ */
#endif /* __GUG_DOS_  */

#endif /* __GUH_GRAPH_ */
