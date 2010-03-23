
/************************************************************************\ 
 *									*
 * pgraph.h - graphic routines for postscript output			*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see pgraph.c for copying information					*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
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
void QViewport( float left , float top , float right , float bottom );
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
void QTextPointSize( int size );
void QTextRealSize( float size );
void QTextDimensions( char *s , float *width , float *height );
void QRealXY( float vx , float vy , float *x , float *y );
void QScreenXY( float x , float y , float *vx , float *vy );
void QFlush( void );
void QSync( void );
void QTest( void );
void QSyncron( int sync );

/*
 *	additional 0.5 cancelled, since no integer conversion for PS
 */

#if __GUG_DOS_

#define QxWtoV(x) (SxWtoV*((x)-WinXmin))
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom-ViewTop)
#define QxVtoW(x) (SxVtoW*(x)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom+ViewTop)+WinYmin)

#else
#if __GUG_UNIX_

#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+ViewLeft)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom)
#define QxVtoW(x) (SxVtoW*((x)-ViewLeft)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom)+WinYmin)

#endif /* __GUG_UNIX_ */
#endif /* __GUG_DOS_  */

#endif /* __GUH_GRAPH_ */
