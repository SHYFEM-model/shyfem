
/* $Id: screen.c,v 1.2 2009-01-14 17:16:37 georg Exp $ */

/************************************************************************\ 
 *									*
 * screen.c - screen manipulation routines				*
 *									*
 * Copyright (c) 1998 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
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

