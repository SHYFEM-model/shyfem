
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
 * screen.c - screen manipulation routines				*
 *									*
 * Revision History:							*
 * 02-Apr-1998: version not yet finished but working                    *
 * 28-Mar-1998: routines adapted from graph.c                           *
 *									*
\************************************************************************/


#include <stdlib.h>

#include "graph.h"


/* the factors 0.4 and 0.6 account for the truncation error
	in graph.c */		/* HACK */

#define YTRANS(y)	( (float) (YMax + YMin - (y) + 0.4) )
#define XTRANS(x)	( (float) (x - 0.6) )


typedef struct {
	int xmin;
	int ymin;
	int xmax;
	int ymax;
} IRect;


static int YMin = 0;
static int YMax = 0;


IRect *ToIRect( IRect *irect , int xmin , int ymin , int xmax , int ymax )

{
	if( !irect ) {
	    irect = (IRect *) malloc( sizeof(IRect) );
	    if( !irect )
		Error("Cannot allocate IRect");
	}

	irect->xmin = xmin;
	irect->ymin = ymin;
	irect->xmax = xmax;
	irect->ymax = ymax;

	return irect;
}

void SViewport( IRect *view )

{
	int xmin = view->xmin;
	int ymin = view->ymin;
	int xmax = view->xmax;
	int ymax = view->ymax;

	YMin = ymin;
	YMax = ymax;

	QViewport( xmin , ymin , xmax , ymax );
	QWindow( (float) xmin , (float) ymin , (float) xmax , (float) ymax );
}

void SMove( int x , int y )

{
	QMove( XTRANS(x) , YTRANS(y) );
}

void SPlot( int x , int y )

{
	QPlot( XTRANS(x) , YTRANS(y) );
}

void SLine( int x1 , int y1 , int x2 , int y2 )

{
	QLine( XTRANS(x1) , YTRANS(y1) , XTRANS(x2) , YTRANS(y2) );
}

void SPlotRect( IRect *irect )

{
	float xmin = XTRANS(irect->xmin);
	float ymin = YTRANS(irect->ymin);
	float xmax = XTRANS(irect->xmax);
	float ymax = YTRANS(irect->ymax);

	QMove( xmin , ymin );
	QPlot( xmin , ymax );
	QPlot( xmax , ymax );
	QPlot( xmax , ymin );
	QPlot( xmin , ymin );
}


void SFillRect( IRect *irect , int color )

{
	float xmin = XTRANS(irect->xmin);
	float ymin = YTRANS(irect->ymin);
	float xmax = XTRANS(irect->xmax);
	float ymax = YTRANS(irect->ymax);

	QRectFill(xmin,ymin,xmax,ymax,color);
}

void SShadeRect( IRect *irect , int darkcolor , int lightcolor )

{
	float xmin = XTRANS(irect->xmin);
	float ymax = YTRANS(irect->ymin);
	float xmax = XTRANS(irect->xmax);
	float ymin = YTRANS(irect->ymax);

	/* printf("SShadeRect (0) : %d %d %d %d\n",
		irect->xmin,irect->ymin,irect->xmax,irect->ymax); */
	/* printf("SShadeRect (1) : %f %f %f %f\n",xmin,ymin,xmax,ymax); */

	/* first pass */

	QNewPen(darkcolor);
	QLine( xmin , ymin , xmax , ymin );	/* top */
	QLine( xmax , ymin , xmax , ymax );	/* right */
	QNewPen(lightcolor);
	QLine( xmax , ymax , xmin , ymax );	/* bottom */
	QLine( xmin , ymax , xmin , ymin );	/* left */

	/* second pass -> thicken line */

	xmin += 1.; ymin += 1.;
	xmax -= 1.; ymax -= 1.;

	QNewPen(darkcolor);
	QLine( xmin , ymin , xmax , ymin );	/* top */
	QLine( xmax , ymin , xmax , ymax );	/* right */
	QNewPen(lightcolor);
	QLine( xmax , ymax , xmin , ymax );	/* bottom */
	QLine( xmin , ymax , xmin , ymin );	/* left */

	/* printf("SShadeRect (2) : %f %f %f %f\n",xmin,ymin,xmax,ymax); */
}

/**************************************************************************/

void SText( int x , int y , char *s )

{
	QText( XTRANS(x) , YTRANS(y) , s );
}

void STextSize( int size )

{
	QTextSize( size );
}

void STextBackground( int color )

{
	QTextBackground( color );
}

void STextDimensions( char *s , int *width , int *height )

{
	QTextDimensionsI( s , width , height );
}

void SCenterText( IRect *irect , char *s , int *x0 , int *y0 )

{
	int w,h;

	STextDimensions( s , &w , &h );

	*x0 = ( irect->xmin + irect->xmax - w ) / 2 ;
	*y0 = ( irect->ymin + irect->ymax + h ) / 2 ;
}

/**************************************************************************/

void SNewPen( int color )

{
	QNewPen( color );
}

/**************************************************************************/

