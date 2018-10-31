
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
 * ggraph.c - graphic routines for gcc under DOS			*
 *									*
 * Revision History:							*
 * 13-Feb-1998: adjourned to GRX2.0                                     *
 * 11-Dec-94: routines written from graph				*
 *									*
\************************************************************************/

/*
#define GGU_DEBUG 1
*/
#define GRX_VERSION_GGU	2

#include <stdio.h>
#include <stdlib.h>

#if GRX_VERSION_GGU == 1
#include <grx.h>
#else
#include <grx20.h>
#endif

#include <grdriver.h>

#include "graph.h"


#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+ViewLeft+0.5)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom+0.5)
#define QxVtoW(x) (SxVtoW*((x)-ViewLeft)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom)+WinYmin)



/* definitions for colors */


int Black         = 0;
int Blue          = 1;
int Green         = 2;
int Cyan          = 3;
int Red           = 4;
int Magenta       = 5;
int Brown         = 6;
int LightGray     = 7;
int DarkGray      = 8;
int LightBlue     = 9;
int LightGreen    = 10;
int LightCyan     = 11;
int LightRed      = 12;
int LightMagenta  = 13;
int Yellow        = 14;
int White         = 15;

/*****************************************************************/

#define	MAXCOLOR	256


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

static int ActColor   = 15;	/* pointer into color map */
static int Color      = 15;	/* color used by gcc */
static int TextSize   = 1;

static GrTextOption TextOptionStructure;
static GrTextOption *TextOption;

static float XAct=0.;
static float YAct=0.;

static MaxColors=0;
static int MyColors[MAXCOLOR];

void QAllocColors( void );

#ifdef GGU_DEBUG
#include "debug.h"
#endif

/*
typedef struct {
	int width;
	int height;
	void *buffer;
} Pixmap_type;
*/

/*****************************************************************/

void QGraphInit( void )

{
	int w,h;
	int col,colfree;
	int white,black;

	GrSetMode(GR_default_graphics);

#ifdef GGU_DEBUG
	GDS("GraphInit\n");
#endif

	w=GrScreenX();
	h=GrScreenY();
	col=GrNumColors();
	colfree=GrNumFreeColors();

	QAllocColors();

	white=White;
	black=Black;

	XMaxPix = w-1;
	YMaxPix = h-1;

	GrClearScreen(black);
	Color=GrWhite();
	ActColor=15;

	TextOption = &TextOptionStructure;

#ifdef GGU_DEBUG
	GDS("Font set up - befor\n");
#endif

#if GRX_VERSION_GGU == 1
	TextOption = GrFindBestFont(8,8,1,"pc",TextOption);
#else
	TextOption->txo_font = GrLoadFont("pc8x8.fnt");
#endif

#ifdef GGU_DEBUG
	GDS("Font set up\n");
#endif

	if( !TextOption ) {
		GrSetMode(GR_80_25_text);
		printf("Cannot Allocate Font\n");
		exit(1);
	}

	TextOption->txo_xalign = GR_ALIGN_LEFT;
	TextOption->txo_yalign = GR_ALIGN_BOTTOM;
	TextOption->txo_chrtype = GR_BYTE_TEXT;
	TextOption->txo_fgcolor.v = white;
	TextOption->txo_bgcolor.v = black;

}

void QClearScreen( void )

{
	GrClearScreen(GrBlack());

#ifdef GGU_DEBUG
	GDS("ClearScreen\n");
#endif
}

void QGraphClose( void )

{
	GrSetMode(GR_80_25_text);

#ifdef GGU_DEBUG
	GDPrint();
#endif
}

void QNewPen( int color )

{
	ActColor=color%MaxColors;
	Color=MyColors[ActColor];

#ifdef GGU_DEBUG
	GDS("NewPen : "); GDI(color); GDNL();
#endif
}

void QGetPen( int *color )

{
	*color=ActColor;
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

	GrSetClipBox(left,top,right,bottom);

	QAdjustScale();
}

void QClearViewport( void )

{
	GrClearClipBox(GrBlack());
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

    GrLine(vx1,vy1,vx2,vy2,Color);

    XAct = vx2;
    YAct = vy2;

#ifdef GGU_DEBUG
	GDS("Line : "); GDI(vx1); GDI(vy1); GDI(vx2); GDI(vy2); GDNL();
#endif
}

void QPlot( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    GrLine(XAct,YAct,vx,vy,Color);

    XAct = vx;
    YAct = vy;
}

void QMove( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    XAct = vx;
    YAct = vy;
}

void QPoint( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    GrFilledBox(vx,vy,vx,vy,Color);

    XAct = vx;
    YAct = vy;
}

typedef int IMP[2];

void QAreaFill( int ndim , float *x , float *y )

{
	static int  nalloc=0;
	static IMP *points=NULL;
	int i;

    if( nalloc == 0 ) {
	nalloc = ndim;
	points = (IMP *) malloc( 2 * ndim * sizeof( int ) );
    } else if( nalloc < ndim ) {
	free( points );
	nalloc = ndim;
	points = (IMP *) malloc( 2 * ndim * sizeof( int ) );
    }

    if( points == NULL ) {
	printf("QAreaFill : Cannot allocate array\n");
	exit(1);
    }

    for(i=0;i<ndim;i++) {
	points[i][0] = QxWtoV(x[i]);
	points[i][1] = QyWtoV(y[i]);
    }

    GrFilledPolygon(ndim,points,Color);
}

void QRectFill( float x1 , float y1 , float x2 , float y2 , int color )

{
	int ix1,iy1,ix2,iy2;

	ix1 = QxWtoV( x1 );
	iy1 = QyWtoV( y1 );
	ix2 = QxWtoV( x2 );
	iy2 = QyWtoV( y2 );

	GrFilledBox(ix1,iy2,ix2,iy1,MyColors[color%MaxColors]);

#ifdef GGU_DEBUG
	GDS("Fill : "); GDI(ix1); GDI(iy1); GDI(ix2); GDI(iy2); 
			GDI(color); GDNL();
#endif
}

void QRectFillI( int x1 , int y1 , int x2 , int y2 , int color )

{
	GrFilledBox(x1,y2,x2,y1,MyColors[color%MaxColors]);

#ifdef GGU_DEBUG
	GDS("Fill : "); GDI(x1); GDI(y1); GDI(x2); GDI(y2); 
			GDI(color); GDNL();
#endif
}


void QText( float x , float y , char *s )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    TextOption->txo_fgcolor.v = Color;
    GrDrawString(s,strlen(s),vx,vy,TextOption);
}

void QTextI( int x , int y , char *s )

{
    TextOption->txo_fgcolor.v = Color;
    GrDrawString(s,strlen(s),x,y,TextOption);
}

void QTextBackground( int color )

{
    TextOption->txo_bgcolor.v = MyColors[color%MaxColors];
}

void QTextSize( int size )

{
    if( size != TextSize && size > 0 ) {
#if GRX_VERSION_GGU == 1
		TextOption->txo_xmag = size;
		TextOption->txo_ymag = size;
#endif
		TextSize = size;
    }
}

void QTextDimensions( char *s , float *width , float *height )

{
	*width  = SxVtoW * GrStringWidth(s,strlen(s),TextOption);
	*height = - ( SyVtoW * GrStringHeight(s,strlen(s),TextOption) );
}

void QTextDimensionsI( char *s , int *width , int *height )

{
	*width  = GrStringWidth(s,strlen(s),TextOption);
	*height = GrStringHeight(s,strlen(s),TextOption);
}


/*
 *	in TurboC the screen coordinates are different from Viewport
 *	coordinates --> add or subtract the lower Viewport coordinate
 *	in order to get the absolute Screen coordinate
 */

void QRealXY( int vx , int vy , float *x , float *y )

{
	*x = QxVtoW( vx );
	*y = QyVtoW( vy );
}

void QScreenXY( float x , float y , int *vx , int *vy )

{
	*vx = QxWtoV( x );
	*vy = QyWtoV( y );
}

void QWindowMaxXY( int *xmax , int *ymax )

{
        *xmax = XMaxPix;
        *ymax = YMaxPix;
}

void QScreenMaxXY( int *xmax , int *ymax )

{
	*xmax = XMaxPix;
	*ymax = YMaxPix;
}

void *QSavePixels( int x , int y , int width , int height ) 

{
	Pixmap_type *bp;
	unsigned char *q;
	int i,j;

	bp = (Pixmap_type *) malloc( sizeof(Pixmap_type) );
	if( !bp ) return NULL;
	q = (unsigned char *) malloc( width*height );
	if( !q ) return NULL;

	bp->buffer = (void *)q;
	bp->width = width;
	bp->height = height;

	for(i=0;i<width;i++) {
	  for(j=0;j<height;j++) {
	    *q++ = (unsigned char) GrPixel(x+i,y+j);
	  }
	}

	return (void *) bp;
}

void QRestorePixels( int x , int y , void *buffer )

{
	Pixmap_type *bp;
	unsigned char *q;
	int i,j;
	int width,height;
	int color;

	if( !buffer ) return;

	bp = (Pixmap_type *) buffer;

	width = bp->width;
	height = bp->height;
	q = (unsigned char *)bp->buffer;

	for(i=0;i<width;i++) {
	  for(j=0;j<height;j++) {
	    color = *q++;
	    GrFilledBox(x+i,y+j,x+i,y+j,color);
	  }
	}

	free( (unsigned char *) bp->buffer );
	free( bp );
	
}

int QAllocColor( int red , int green , int blue )

{
	if( MaxColors >= MAXCOLOR )
		Error("Cannot allocate more colors");

	MyColors[MaxColors++] = GrAllocColor(red,green,blue);

	return MaxColors-1;
}
	
void QAllocColors( void )

{
	MyColors[Black] = GrBlack();
	MyColors[Blue] = GrAllocColor(0,0,255);
	MyColors[Green] = GrAllocColor(0,255,0);
	MyColors[Cyan] = GrAllocColor(0,255,255);
	MyColors[Red] = GrAllocColor(255,0,0);
	MyColors[Magenta] = GrAllocColor(255,0,255);
	MyColors[Brown] = GrAllocColor(165,42,42);
	MyColors[LightGray] = GrAllocColor(211,211,211);
	MyColors[DarkGray] = GrAllocColor(105,105,105);
	MyColors[LightBlue] = GrAllocColor(173,216,230);
	MyColors[LightGreen] = GrAllocColor(32,178,170);
	MyColors[LightCyan] = GrAllocColor(224,255,255);
	MyColors[LightRed] = GrAllocColor(205,92,92);
	MyColors[LightMagenta] = GrAllocColor(255,182,193);
	MyColors[Yellow] = GrAllocColor(255,255,0);
	MyColors[White] = GrWhite();

	MaxColors = 16;

#ifdef GGU_DEBUG
	GDS("Color : "); 
	GDI(Black); GDI(Blue); GDI(Green); GDI(Cyan); 
	GDI(Red); GDI(Magenta); GDI(Brown); GDI(LightGray);
	GDI(DarkGray); GDI(LightBlue); GDI(LightGreen);  GDI(LightCyan);
	GDI(LightRed); GDI(LightMagenta); GDI(Yellow);  GDI(White);
	GDNL();
#endif

}
