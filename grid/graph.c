
/************************************************************************\ 
 *									*
 * graph.c - graphic routines for Turbo C				*
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
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/



#include <graphics.h>
#include <stdio.h>
#include <stdlib.h>

#include "graph.h"



#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+0.5)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom-ViewTop+0.5)
#define QxVtoW(x) (SxVtoW*(x)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom+ViewTop)+WinYmin)



/* definitions for colors */

int Black         = BLACK;
int Blue          = BLUE;
int Green         = GREEN;
int Cyan          = CYAN;
int Red           = RED;
int Magenta       = MAGENTA;
int Brown         = BROWN;
int LightGray     = LIGHTGRAY;
int DarkGray      = DARKGRAY;
int LightBlue     = LIGHTBLUE;
int LightGreen    = LIGHTGREEN;
int LightCyan     = LIGHTCYAN;
int LightRed      = LIGHTRED;
int LightMagenta  = LIGHTMAGENTA;
int Yellow        = YELLOW;
int White         = WHITE;


/*****************************************************************/

static int XMinPix    = 0;
static int YMinPix    = 0;
static int XMaxPix    = 640-1;
static int YMaxPix    = 480-1;

static int ViewLeft   = 0;
static int ViewTop    = 0;
static int ViewRight  = 0;
static int ViewBottom = 0;

static float WinXmin  = 0.;
static float WinYmin  = 0.;
static float WinXmax  = 0.;
static float WinYmax  = 0.;

static float SxWtoV   = 0.;
static float SyWtoV   = 0.;
static float SxVtoW   = 0.;
static float SyVtoW   = 0.;

static int Color      = 15;
static int TextSize   = 1;

/*****************************************************************/

void QGraphInit( void )

{
	int g_driver,g_mode,g_error;

	/* load graphics driver for EGA/VGA */

	if( registerfarbgidriver(EGAVGA_driver_far) < 0 ) {
		printf("Cannot register driver for EGA/VGA driver");
		exit(1);
	}

	/* detect graphics adapter */

	detectgraph(&g_driver,&g_mode);
	if( g_driver < 0 ) {
		printf("No graphics hardware detected !\n");
		exit(1);
	}
/*
    printf("Detected graphics driver and mode :\n");
    printf("%d %d\n",g_driver,g_mode);
*/
	/* initialize monitor */

/*      initgraph(&g_driver,&g_mode,"c:\\tc\\");*/
	initgraph(&g_driver,&g_mode,"");
	g_error=graphresult();
	if( g_error < 0 ) {
		printf("initgraph error : %s.\n",grapherrormsg(g_error));
		exit(1);
	}

    cleardevice();

    Color=15;
	setcolor(Color);
	setfillstyle(SOLID_FILL,Color);

    settextstyle(0,0,TextSize);
    settextjustify(LEFT_TEXT,BOTTOM_TEXT);

}

void QClearScreen( void )

{
    cleardevice();
}

void QGraphClose( void )

{
    closegraph();
}

void QNewPen( int color )

{
    Color=color;
	setcolor(Color);
	setfillstyle(SOLID_FILL,Color);
}

void QGetPen( int *color )

{
	*color=Color;
}

void QGetPixelWidth( float *dx , float *dy )

{
	*dx = QxVtoW(1)-QxVtoW(0);
	*dy = QyVtoW(0)-QyVtoW(1);
}

void QViewport( int left , int top , int right , int bottom )

{
    ViewLeft   = left;
	ViewTop    = top;
    ViewRight  = right;
    ViewBottom = bottom;

    setviewport(left,top,right,bottom,1);

    QAdjustScale();
}

void QClearViewport( void )

{
    clearviewport();
}

void QWindow( float xmin , float ymin , float xmax , float ymax )

{
    WinXmin = xmin;
    WinYmin = ymin;
    WinXmax = xmax;
    WinYmax = ymax;

    QAdjustScale();
}

void QAdjustScale( void )

{
    /* in y direction values are interchanged since top < bottom */

    if( WinXmax-WinXmin > 0. )
	SxWtoV = (ViewRight-ViewLeft)/(WinXmax-WinXmin);

    if( WinYmax-WinYmin > 0. )
	SyWtoV = (ViewTop-ViewBottom)/(WinYmax-WinYmin);

    if( SxWtoV != 0. )
	SxVtoW = 1./SxWtoV;

    if( SyWtoV != 0. )
		SyVtoW = 1./SyWtoV;

}

void QMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax )

{
    *xmin = XMinPix;
    *ymin = YMinPix;
    *xmax = XMaxPix;
    *ymax = YMaxPix;
}

void QLine( float x1 , float y1 , float x2 , float y2 )

{
    int vx1,vy1,vx2,vy2;

    vx1=QxWtoV(x1);
    vy1=QyWtoV(y1);
    vx2=QxWtoV(x2);
    vy2=QyWtoV(y2);

    line(vx1,vy1,vx2,vy2);
}

void QPlot( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    lineto(vx,vy);
}

void QMove( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    moveto(vx,vy);
}

void QPoint( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    putpixel(vx,vy,Color);
}

void QAreaFill( int ndim , float *x , float *y )

{
	static int  nalloc=0;
	static int *points=NULL;
	int i;

    if( nalloc == 0 ) {
	nalloc = ndim;
	points = (int *) malloc( 2 * ndim * sizeof( int ) );
    } else if( nalloc < ndim ) {
	free( points );
	nalloc = ndim;
	points = (int *) malloc( 2 * ndim * sizeof( int ) );
    }

    if( points == NULL ) {
	printf("QAreaFill : Cannot allocate array\n");
	exit(1);
    }

    for(i=0;i<ndim;i++) {
	points[2*i]   = QxWtoV(x[i]);
	points[2*i+1] = QyWtoV(y[i]);
    }

    fillpoly(ndim,points);
}

void QRectFill( float x1 , float y1 , float x2 , float y2 , int color )

{
	int ix1,iy1,ix2,iy2,tmpcol;

	ix1 = QxWtoV( x1 );
	iy1 = QyWtoV( y1 );
	ix2 = QxWtoV( x2 );
	iy2 = QyWtoV( y2 );

	tmpcol=Color;
	QNewPen(color);
	bar(ix1,iy2,ix2,iy1);
	QNewPen(tmpcol);
}

void QRectFillI( int x1 , int y1 , int x2 , int y2 , int color )

{
    int tmpcol;

	tmpcol=Color;
	QNewPen(color);
    bar(x1,y2,x2,y1);
	QNewPen(tmpcol);
}


void QText( float x , float y , char *s )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    outtextxy(vx,vy,s);
}

void QTextI( int x , int y , char *s )

{
    outtextxy(x,y,s);
}


void QTextSize( int size )

{
    if( size != TextSize && size > 0 ) {
		settextstyle(0,0,size);
		TextSize = size;
    }
}

void QTextDimensions( char *s , float *width , float *height )

{
    *width  = SxVtoW * textwidth ( s );
	*height = - ( SyVtoW * textheight( s ) );
}

void QTextDimensionsI( char *s , int *width , int *height )

{
    *width  = textwidth ( s );
    *height = textheight( s );
}


/*
 *	in TurboC the screen coordinates are different from Viewport
 *	coordinates --> add or subtract the lower Viewport coordinate
 *	in order to get the absolute Screen coordinate
 */

void QRealXY( int vx , int vy , float *x , float *y )

{
	*x = QxVtoW( vx - ViewLeft);
	*y = QyVtoW( vy - ViewTop );
}

void QScreenXY( float x , float y , int *vx , int *vy )

{
	*vx = QxWtoV( x ) + ViewLeft;
	*vy = QyWtoV( y ) + ViewTop ;
}

void *QSavePixels( int x , int y , int width , int height )

{
    void *new;
    int size;

    size = imagesize(x,y,x+width-1,y+height-1);
    new = (void *) malloc( size );

    if( new ) {
	MouseHide();
	getimage(x,y,x+width-1,y+height-1,new);
	MouseShow();
    }

    return new;
}

void QRestorePixels( int x , int y , void *buffer )

{
    if( buffer ) {
	MouseHide();
        putimage(x,y,buffer,COPY_PUT);
	MouseShow();
	free(buffer);
    }
}
