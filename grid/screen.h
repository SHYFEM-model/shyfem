
/* $Id: screen.h,v 1.2 2009-01-14 17:16:37 georg Exp $ */

/************************************************************************\ 
 *									*
 * screen.h - screen manipulation routines				*
 *									*
 * Copyright (c) 1998 by Georg Umgiesser				*
 *									*
 * see screen.c for copying information                                 *
 *                                                                      *
 * Revision History:                                                    *
 * 02-Apr-1998: version not yet finished but working                    *
 * 28-Mar-1998: routines adapted from graph.c                           *
 *                                                                      *
\************************************************************************/


#ifndef __GUH_SCREEN_
#define __GUH_SCREEN_


typedef struct {
	int xmin;
	int ymin;
	int xmax;
	int ymax;
} IRect;


IRect *ToIRect( IRect *irect , int xmin , int ymin , int xmax , int ymax );

void SViewport( IRect *view );

void SMove( int x , int y );
void SPlot( int x , int y );
void SLine( int x1 , int y1 , int x2 , int y2 );

void SPlotRect( IRect *irect );
void SFillRect( IRect *irect , int color );
void SShadeRect( IRect *irect , int darkcolor , int lightcolor );

void SText( int x , int y , char *s );
void STextSize( int size );
void STextBackground( int color );
void STextDimensions( char *s , int *width , int *height );

void SCenterText( IRect *irect , char *s , int *x0 , int *y0 );

void SNewPen( int color );


#endif
