
/************************************************************************\ 
 *									*
 * xgraphf.c - HCBS simulation routines for FORTRAN under X11		*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
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
 * 21-Mar-94: include xgraph.h						*
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


/* 
 *	contents :
 *
 *	qopen()			opens window
 *	qclose()		closes window
 *	qstart()		opens a new plot (clears screen)
 *	qend()			closes a plot (flushes the buffer)
 *	qnewp(ic)		changes color (ic=[1,16], 1=black, 16=white)
 *	qline(x1,y1,x2,y2)	draws line from (1) to (2)
 *	qplot(x,y)		draws line from actual position to (x,y)
 *	qmove(x,y)		moves to (x,y)
 *	qpoint(x,y)		draws point at (x,y)
 *	qafill(n,x,y)		fills x,y (n points) with actual color
 *	qrfill(x1,y1,x2,y2,ic)	fills rectanle with color ic
 *	qout()			flushes buffer
 *
 *	additional routines (windowing)
 *
 *	qgetvp(xmin,ymin,xmax,ymax)	gets viewport dimensions (in cm)
 *	qviewp(xmin,ymin,xmax,ymax)	defines window for clipping (in cm)
 *	qclrvp()			clears viewport
 *	qworld(xmin,ymin,xmax,ymax)	sets world coordinates
 *	setxwd(xmin,ymin,xmax,ymax)	sets window coordinates
 *	getxwd(xmin,ymin,xmax,ymax)	gets window coordinates
 *	
 */


#ifndef __GUC_XGRAPHF_
#define __GUC_XGRAPHF_


#include <stdio.h>
#include <stdlib.h>

#include "gustd.h"
#include "graph.h"
#include "xgraph.h"


/*****************************************************************/

/* pixels per cm */

#define PIXPERCM 35.	

static float XMinView;
static float YMinView;
static float XMaxView;
static float YMaxView;

static int Qopen=0;
static int Qstart=0;
static int Qnumb=0;

static int TEST = 0;	/* testing routines */
static int NORMC = 0;	/* 0 screen coord. in cm, 1 normalized screen coord. */

void QSetViewCoord( int norm );

/*****************************************************************/

void qopen_( void )

{
	char *s;
	int xmin,ymin,xmax,ymax;

	QSetTitle("xhcbs");
	QGraphInit();
	QGetMinMaxPix( &xmin , &ymin , &xmax , &ymax );
	QSetViewCoord(NORMC);
	QWindow( XMinView , YMinView , XMaxView , YMaxView );
	QViewport( 0 , 0 , xmax - xmin , ymax - ymin );
	printf("HCBS for X11 - written by Georg Umgiesser");
	printf(" - ISDGM/CNR\n");
	Qopen = 1;
	QSync();
	printf("Enter <CR> to continue : ");
	s=getlin(stdin);
	if( s ) s = NULL;
}

void qclose_( void )

{
	if( Qopen ) {
		Qopen = 0;
		QGraphClose();
	} else {
		printf("qclose : No graphic window opened\n");
	}
}

void qstart_( void )

{
	if( Qopen ) {
		Qnumb++;
		Qstart = 1;
		QClearScreen();
	} else {
		printf("qstart : No graphic window opened\n");
		exit(1);
	}
}

void qend_( void )

{
	if( Qstart ) {
		Qstart = 0;
		QSync();
	}
}

void qnewp_( long int *color )

{
	QNewPen( (int) (*color-1) );
}

void qgetvp_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	*xmin = XMinView;
	*ymin = YMinView;
	*xmax = XMaxView;
	*ymax = YMaxView;
}

void qviewp_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	int left,top,right,bottom;
	int xminp,yminp,xmaxp,ymaxp;

	QGetMinMaxPix( &xminp , &yminp , &xmaxp , &ymaxp );

	if( *xmin==0. && *ymin==0. && *xmax==0. && *ymax==0. ) {
		left=0;
		top=0;
		right=xmaxp-xminp;
		bottom=ymaxp-yminp;
	} else {
		left=(xmaxp-xminp+1)*(*xmin)/(XMaxView-XMinView)-0.5;
		right=(xmaxp-xminp+1)*(*xmax)/(XMaxView-XMinView)-0.5;
		top=(ymaxp-yminp+1)*(*ymax)/(YMaxView-YMinView)-0.5;
		top = (ymaxp-yminp+1) - top;
		bottom=(ymaxp-yminp+1)*(*ymin)/(YMaxView-YMinView)-0.5;
		bottom = (ymaxp-yminp+1) - bottom;
	}

	if(TEST) printf("qviewp: %d %d %d %d\n",left,top,right,bottom);

	QViewport( left , top , right , bottom );
}

void qclrvp_( void )

{
	QClearViewport();
}

void qworld_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	QWindow( *xmin , *ymin , *xmax , *ymax );
}

void QSetViewCoord( int norm )

{
	int xmin,ymin,xmax,ymax;

	QGetMinMaxPix( &xmin , &ymin , &xmax , &ymax );

	XMinView = 0.;
	YMinView = 0.;

	if( norm ) {
		YMaxView = 1.;
		XMaxView = ((float)(xmax-xmin+1))/((float)(ymax-ymin+1));
	} else {
		XMaxView = (xmax-xmin+1)/PIXPERCM;
		YMaxView = (ymax-ymin+1)/PIXPERCM;
	}
}

void setxwd_( long *xmin , long *ymin , long *xmax , long *ymax )

{
	if( !Qopen ) { /* only befor window has been opened */
		QSetMinMaxPix( (int) *xmin , (int) *ymin 
			, (int) *xmax , (int) *ymax );
		QSetViewCoord(NORMC);
	}
}

void getxwd_( long *xmin , long *ymin , long *xmax , long *ymax )

{
	int xminh,yminh,xmaxh,ymaxh;

	QGetMinMaxPix( &xminh , &yminh , &xmaxh , &ymaxh );

	*xmin=xminh; *ymin=yminh; *xmax=xmaxh; *ymax=ymaxh;
}

void qline_( float *x1 , float *y1 , float *x2 , float *y2 )

{
	QLine( *x1 , *y1 , *x2 , *y2 );
}

void qplot_( float *x , float *y )

{
	QPlot( *x , *y );
}

void qmove_( float *x , float *y )

{
	QMove( *x , *y );
}

void qpoint_( float *x , float *y )

{
	QPoint( *x , *y );
}

void qafill_( long *n , float *x , float *y )

{
	QAreaFill( (int) *n , x , y );
}

void qrfill_( float *x1 , float *y1 , float *x2 , float *y2 , long *color )

{
	QRectFill( *x1 , *y1 , *x2 , *y2 , (int) *color-1 );
}

void qtest_( long *test )

{
	if( *test==0 )
		TEST=0;
	else
		TEST=1;
}

void qout_( void )

{
	QSync();
}

#endif /* __GUC_XGRAPHF_ */
