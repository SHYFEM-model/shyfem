
/************************************************************************\ 
 *									*
 * graph.h - general header file for graphic routines			*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see xgraph.c for copying information					*
 *									*
 * Revision History:							*
 * 05-Dec-95: changes (Background, Save/RestorePixel ...) included      *
 *              (QGetViewport, QGetWindow, QBell...)                    *
 * 15-May-94: Colors not int's but defines                              *
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRAPH_
#define __GUH_GRAPH_

#include "general.h"

/* definitions for colors */

extern int Black        ;
extern int Blue         ;
extern int Green        ;
extern int Cyan         ;
extern int Red          ;
extern int Magenta      ;
extern int Brown        ;
extern int LightGray    ;
extern int DarkGray     ;
extern int LightBlue    ;
extern int LightGreen   ;
extern int LightCyan    ;
extern int LightRed     ;
extern int LightMagenta ;
extern int Yellow       ;
extern int White        ;

typedef struct {
        int width;
        int height;
        void *buffer;
        int xbuffer;
} Pixmap_type;

void QGraphInit( void );
void QClearScreen( void );
void QGraphClose( void );
void QNewPen( int color );
void QGetPen( int *color );
void QGetPixelWidth( float *dx , float *dy );
void QViewport( int left , int top , int right , int bottom );
void QGetViewport( int *left , int *top , int *right , int *bottom );
void QClearViewport( void );
void QWindow( float xmin , float ymin , float xmax , float ymax );
void QGetWindow( float *xmin , float *ymin , float *xmax , float *ymax );
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
void QTextBackground( int color );
void QTextSize( int size );
void QTextDimensions( char *s , float *width , float *height );
void QTextDimensionsI( char *s , int *width , int *height );
void QRealXY( int vx , int vy , float *x , float *y );
void QScreenXY( float x , float y , int *vx , int *vy );
void QScreenMaxXY( int *xmax , int *ymax );
void QWindowMaxXY( int *xmax , int *ymax );

void QBell( void );

void *QSavePixels( int x , int y , int width , int height );
void QRestorePixels( int x , int y , void *buffer );

void QFlush( void );
void QSync( void );
void QTest( void );
void QSyncron( int sync );

int QAllocColor( int red , int green , int blue );

#endif /* __GUH_GRAPH_ */
