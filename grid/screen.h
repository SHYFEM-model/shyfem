
/************************************************************************\
 *
 *    Copyright (C) 1998  Georg Umgiesser
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
 * screen.h - screen manipulation routines
 *
 * revision log :
 *
 * 28.03.1998	ggu	routines adapted from graph.c
 * 02.04.1998	ggu	version not yet finished but working
 *
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
