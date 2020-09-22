
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
 * xgraph.c - graphic routines for X11
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 * 21.03.1994	ggu	VelColors used for call to QAllocVelColors(),
 * 21.03.1994	ggu	*ShadeColor() declared in xgraph.h,
 * 21.03.1994	ggu	gustd.h not included anymore
 * 21.03.1994	ggu	gcc-warnings, call to XGeometry corrected (unsigned)
 *
\************************************************************************/


#ifndef __GUC_XGRAPH_
#define __GUC_XGRAPH_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "xgraph.h"

#include "generalx.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#define __UseXt	0

/* definitions for colors */

#define MAXCOL 250

/*
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
*/

int Black         = 0;
int Blue          = 1;
int Green         = 2;
int Cyan          = 3;
int Red           = 4;
int Magenta       = 5;
int Brown         = 6;
int LightGray     = 7;
int DarkGray      = 8;		/* DimGray or gray */
int LightBlue     = 9;
int LightGreen    = 10;
int LightCyan     = 11;
int LightRed      = 12;
int LightMagenta  = 13;
int Yellow        = 14;
int White         = 15;


static void QAllocColors( void );
static void QAllocBlueColors( void );
static void QAllocRed2YellowColors( void );
static void QAllocYellow2GreenColors( void );
static void QAllocGreen2BlueColors( void );
static void QAllocVelColors( void );

/*
int BlueShadeColor( int shade );
int RedShadeColor( int shade );
int YellowShadeColor( int shade );
int GreenShadeColor( int shade );
int VelShadeColor( int shade );
*/

static char *savestr(char *s );

static char *MyColorNames[] =
		{"black","blue","green","cyan","red","magenta",
		 "brown","LightGray","DimGray",
		 "LightBlue","LightSeaGreen","LightCyan",
		 "IndianRed","LightPink","yellow","white"
		};

static int BlueShade[16];
static int RedShade[16];
static int YellowShade[16];
static int GreenShade[16];
static int VelShade[29];
static int MaxColors=0;
static int DimColors=MAXCOL-1; /* dimension of color table */

long int linear(long int start , long int end , long int steps , long int act );
long int quadratic_slow
			(
			  long int start 
			, long int end 
			, long int steps 
			, long int deriv
			, long int act 
			);

long int quadratic_fast
			(
			  long int start 
			, long int end 
			, long int steps 
			, long int deriv
			, long int act 
			);



/*****************************************************************/

/* absolute screen coordinates of window */

static int XMinPix    = 100;
static int YMinPix    = 100;
static int XMaxPix    = 100+640-1;
static int YMaxPix    = 100+480-1;

/* viewport coordinates relative to window - (left,top) = (0,0) */

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
/*
static int TextHeight = 0;
*/

static int           XAct=0;
static int           YAct=0;

static int	     Initialized=0;

static Display       *MyDisplay;
static int           MyScreen;
static Window        MyRootWindow;
static Window        MyWindow;
static GC            MyGc;
static Widget        MyWidget=NULL;
/*
static XEvent        MyEvent;
static KeySym        MyKey;
*/
static unsigned int  MyBorder = 5;
static unsigned long MyBorderColor;
static unsigned long MyBackGround;
static unsigned long MyForeGround;
static Cursor        MyCursor;
static XFontStruct   *MyFontStruct;
static Colormap      MyColormap;
static XColor        MyColors[MAXCOL]; /* error fixed 1.7.93 */
static unsigned int  MyLineWidth=0; /* set to 1 for pixel reproducibility */
static int           MyFontAscent;
static int           MyFontDescent;
static XCharStruct   MyCharStruct;
static char          *MyDisplayName=NULL;
static Pixmap	     MyIcon=None;

#define grid_width 40
#define grid_height 40
static unsigned char grid_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc0,
   0x00, 0x00, 0x00, 0x00, 0x60, 0x03, 0x00, 0x00, 0x00, 0x90, 0x0c, 0x00,
   0x00, 0x00, 0x88, 0x30, 0x00, 0x00, 0x00, 0x88, 0xc0, 0x00, 0x00, 0x00,
   0x04, 0x01, 0x03, 0x00, 0x00, 0x02, 0x01, 0x0c, 0x00, 0x00, 0x01, 0x01,
   0x30, 0x00, 0x80, 0x00, 0x02, 0xc0, 0x00, 0x40, 0x00, 0x02, 0xe0, 0x00,
   0x40, 0x00, 0x02, 0x10, 0x01, 0x20, 0x00, 0x02, 0x08, 0x01, 0x10, 0x00,
   0x04, 0x06, 0x02, 0x38, 0x00, 0x04, 0x01, 0x02, 0xd0, 0x07, 0x84, 0x00,
   0x02, 0x20, 0xf8, 0x69, 0x00, 0x04, 0x20, 0x00, 0x1e, 0x00, 0x04, 0x40,
   0x00, 0xe8, 0x00, 0x04, 0x80, 0x00, 0x08, 0x07, 0x08, 0x00, 0x01, 0x04,
   0x78, 0x08, 0x00, 0x02, 0x04, 0x80, 0x13, 0x00, 0x04, 0x04, 0x00, 0x1c,
   0x00, 0x04, 0x04, 0x00, 0x06, 0x00, 0x08, 0x02, 0xc0, 0x01, 0x00, 0x10,
   0x02, 0x30, 0x00, 0x00, 0x20, 0x02, 0x0e, 0x00, 0x00, 0x40, 0x82, 0x01,
   0x00, 0x00, 0x40, 0x71, 0x00, 0x00, 0x00, 0x80, 0x0d, 0x00, 0x00, 0x00,
   0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x8c, 0x49, 0x00, 0x5c, 0x3a, 0x42, 0x48, 0x00, 0xc2, 0x4a,
   0x5a, 0x4b, 0x00, 0xc2, 0x3b, 0x52, 0x4a, 0x00, 0x42, 0x2b, 0x8c, 0x31,
   0x00, 0x5c, 0x4a, 0x00, 0x00, 0x00, 0x00, 0x00};



/* static char          MyText[] = {"Georg Umgiesser - ISDGM/CNR"}; */
static char          MyText[] = {"grid"};
static char          *MyTitle=NULL;
static char          *MyCommand[] = {"main"};

static int           TEST=0;
static char          *NullName = "";
static int           VelColors = 0; /* allocate VelColors (y/n) -> (1/0) */

/*****************************************************************/

void QGraphInit( void )

{
    int dummy;

    XSizeHints myhint;
    XWMHints   mywmhint;

    if(MyDisplayName == NULL)
	MyDisplayName = NullName;
    if(MyTitle == NULL)
	MyTitle = MyText;

#if __UseXt != 0
    if(MyWidget)
        MyDisplay  = XtDisplay(MyWidget);
    else
        MyDisplay  = XOpenDisplay(MyDisplayName);
#else
    MyDisplay  = XOpenDisplay(MyDisplayName);
#endif

    MyScreen   = DefaultScreen(MyDisplay);
    MyRootWindow = DefaultRootWindow(MyDisplay);

    MyBackGround  = WhitePixel(MyDisplay,MyScreen);
    MyForeGround  = BlackPixel(MyDisplay,MyScreen);
    MyBorderColor = MyForeGround;

    MyIcon = XCreatePixmapFromBitmapData(
			MyDisplay,MyRootWindow,
			(char *) grid_bits,grid_width,grid_height,
			MyForeGround,MyBackGround,1);

    myhint.x = XMinPix;
    myhint.y = YMinPix;
    myhint.width  = XMaxPix - XMinPix + 1;
    myhint.height = YMaxPix - YMinPix + 1;
    myhint.flags  = PPosition | PSize;

    mywmhint.flags = InputHint | IconPixmapHint;
    mywmhint.input = True;
    mywmhint.icon_pixmap = MyIcon;

    if(MyWidget) {
#if __UseXt != 0
	MyWindow = XtWindow(MyWidget);
#endif
        ;
    } else {
	MyWindow = XCreateSimpleWindow(MyDisplay,MyRootWindow
                ,myhint.x,myhint.y,myhint.width,myhint.height
                ,MyBorder,MyBorderColor,MyBackGround);
	XSetStandardProperties(MyDisplay,MyWindow,MyTitle,MyTitle
                ,MyIcon,MyCommand,1,&myhint);
	XSetWMHints(MyDisplay,MyWindow,&mywmhint);
    }

    MyGc = XCreateGC(MyDisplay,MyWindow,0,0);
    XSetBackground(MyDisplay,MyGc,MyBackGround);
    XSetForeground(MyDisplay,MyGc,MyForeGround);

/*    MyCursor=XCreateFontCursor(MyDisplay,XC_cross); */
/*    MyCursor=XCreateFontCursor(MyDisplay,XC_crosshair); */
    MyCursor=XCreateFontCursor(MyDisplay,XC_top_left_arrow); 
/*    MyCursor=XCreateFontCursor(MyDisplay,XC_crosshair); */
    XDefineCursor(MyDisplay,MyWindow,MyCursor);

    XSetFillRule(MyDisplay,MyGc,WindingRule); /* p. 178 */
    XSetFillStyle(MyDisplay,MyGc,FillSolid); /* p. 195 */
    XSetFunction(MyDisplay,MyGc,GXcopy); /* p. 158 */
    XSetLineAttributes(MyDisplay,MyGc,MyLineWidth,LineSolid
			,CapRound,JoinRound); /* p. 158 */

    XSelectInput(MyDisplay,MyWindow,
                    ButtonPressMask | KeyPressMask | ExposureMask 
		  | StructureNotifyMask );

    MyFontStruct = XLoadQueryFont(MyDisplay,"fixed");
    if(MyFontStruct == 0) {
        printf("Cannot load font\n");
        exit(1);
    }
    XSetFont(MyDisplay,MyGc,MyFontStruct->fid);
    XTextExtents(MyFontStruct,"M",1,&dummy,&MyFontAscent,&MyFontDescent,
			&MyCharStruct);

    QAllocColors();
    QAllocBlueColors();
    QAllocRed2YellowColors();
    QAllocGreen2BlueColors();
    QAllocYellow2GreenColors();

    if( VelColors ) QAllocVelColors();

    if( !MyWidget ) {
        XMapRaised(MyDisplay,MyWindow);
        XSync(MyDisplay,0); /* wait till window is raised and ready */
    }

    Initialized = 1;
}

void QClearScreen( void )

{
    XClearWindow(MyDisplay,MyWindow);
}

void QGraphClose( void )

{
    XFreeGC(MyDisplay,MyGc);
    XDestroyWindow(MyDisplay,MyWindow);
    XCloseDisplay(MyDisplay);
}

void QNewPen( int color )

{
	Color=color;
	if(color>=0 && color <= MaxColors )
		XSetForeground(MyDisplay,MyGc,MyColors[color].pixel);
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
    XRectangle r;

    ViewLeft   = left;
    ViewTop    = top;
    ViewRight  = right;
    ViewBottom = bottom;

    r.x      = left;
    r.y      = top;
    r.width  = right-left+1;
    r.height = bottom-top+1;

    XSetClipRectangles(MyDisplay,MyGc,0,0,&r,1,Unsorted);

    QAdjustScale();
}

void QClearViewport( void )

{
    XClearArea(MyDisplay,MyWindow,
                ViewLeft,ViewTop,
                (unsigned int) ViewRight-ViewLeft+1,
                (unsigned int) ViewBottom-ViewTop+1,
                False);
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

void QGetMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax )

{
    *xmin = XMinPix;
    *ymin = YMinPix;
    *xmax = XMaxPix;
    *ymax = YMaxPix;
}

void QSetMinMaxPix( int xmin , int ymin , int xmax , int ymax )

{
    XMinPix=xmin;
    YMinPix=ymin;
    XMaxPix=xmax;
    YMaxPix=ymax;
}

void QLine( float x1 , float y1 , float x2 , float y2 )

{
    int vx1,vy1,vx2,vy2;

    vx1=QxWtoV(x1);
    vy1=QyWtoV(y1);
    vx2=QxWtoV(x2);
    vy2=QyWtoV(y2);

    if(TEST) printf("QLine: %d %d %d %d %d\n",vx1,vy1,vx2,vy2,Color);

    XDrawLine(MyDisplay,MyWindow,MyGc,vx1,vy1,vx2,vy2);

    XAct=vx2;
    YAct=vy2;
}

void QPlot( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    if(TEST) printf("QPlot: %d %d %d %d %d\n",XAct,YAct,vx,vy,Color);

    XDrawLine(MyDisplay,MyWindow,MyGc,XAct,YAct,vx,vy);

    XAct=vx;
    YAct=vy;
}

void QMove( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    if(TEST) printf("QMove: %d %d %d\n",vx,vy,Color);

    XAct=vx;
    YAct=vy;
}

void QPoint( float x , float y )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    XDrawPoint(MyDisplay,MyWindow,MyGc,vx,vy);
}

void QAreaFill( int ndim , float *x , float *y )

{
	static int  nalloc=0;
	static XPoint *points=NULL;
	int i;

    if( nalloc == 0 ) {
	nalloc = ndim;
	points = (XPoint *) malloc( ndim * sizeof( XPoint ) );
    } else if( nalloc < ndim ) {
	free( points );
	nalloc = ndim;
	points = (XPoint *) malloc( ndim * sizeof( XPoint ) );
    }

    if( points == NULL ) {
	printf("QAreaFill : Cannot allocate array\n");
	exit(1);
    }

    for(i=0;i<ndim;i++) {
	points[i].x = QxWtoV(x[i]);
	points[i].y = QyWtoV(y[i]);
    }

    XFillPolygon(MyDisplay,MyWindow,MyGc,points,ndim,Complex,CoordModeOrigin);
}

void QRectFill( float x1 , float y1 , float x2 , float y2 , int color )

{
	int ix1,iy1,ix2,iy2;
	int tmpcol=0;
	int ix0,iy0;
	unsigned int ixw,iyw;

	ix1 = QxWtoV( x1 );
	iy1 = QyWtoV( y1 );
	ix2 = QxWtoV( x2 );
	iy2 = QyWtoV( y2 );

	if( color >= 0 ) {
		tmpcol=Color;
		QNewPen(color);
	}

	if( ix1 < ix2 ) {
		ix0 = ix1;
		ixw = ix2-ix1+1;
	} else {
		ix0 = ix2;
		ixw = ix1-ix2+1;
	}

	if( iy1 < iy2 ) {
		iy0 = iy1;
		iyw = iy2-iy1+1;
	} else {
		iy0 = iy2;
		iyw = iy1-iy2+1;
	}

	XFillRectangle(MyDisplay,MyWindow,MyGc,ix0,iy0,ixw,iyw);

	if( color >= 0 )
		QNewPen(tmpcol);
}

void QTextSize( int size )

{
    if( size != TextSize && size > 0 ) {
        TextSize = size;
    }
}

void QText( float x , float y , char *s )

{
    int vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    XDrawString(MyDisplay,MyWindow,MyGc,vx,vy,s,strlen(s));
}

void QTextDimensions( char *s , float *width , float *height )

{
    *width  = SxVtoW * XTextWidth(MyFontStruct,s,strlen(s));
    *height = - ( SyVtoW * MyFontAscent );
}

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

/********************************************************************/

void QFlush( void )

{
	XFlush(MyDisplay);
}

void QSync( void )

{
	XSync(MyDisplay,0);
}

void QSyncron( int sync )

{
	XSynchronize(MyDisplay,sync);
}

void QNextEvent( XEvent *eventp )

{
	XNextEvent(MyDisplay,eventp);
}

Display *QActDisplay( void )

{
	return MyDisplay;
}

Window *QActWindow( void )

{
	return &MyWindow;
}

Window *QRootWindow( void )

{
	return &MyRootWindow;
}

int QActScreen( void )

{
	return MyScreen;
}

int QInitialized( void )

{
	return Initialized;
}

GC QActGc( void )

{
	return MyGc;
}

void QSetTest( void )

{
	TEST = 1;
}

void QClearTest( void )

{
	TEST = 0;
}

static void QAllocColors( void )

{
	XColor color,exact;
	Status result;
	int i;

	if(DimColors<MaxColors+16) 
		Error("xgraph: Error in dimension of color table");

	MyColormap = XDefaultColormap(MyDisplay,MyScreen);

	for(i=0;i<16;i++) {
		result=XAllocNamedColor(MyDisplay,MyColormap,MyColorNames[i]
			,&color,&exact);

		if(result)
			MyColors[i]=color;
		else {
			printf("Cannot allocate color %s\n",
					MyColorNames[i]);
			exit(1);
		}
	}
	MaxColors = 15;
}

static void QAllocVelColors( void )

{
	XColor col;
	Status result;
	int i,icol;
	long unsigned int il,nmax=65536;
	long unsigned int bluemin,bluemax,colmin,colmax;
/*	unsigned long int red,green,blue,pixel;
*/	long int colora;

	icol=29;

	if(DimColors<MaxColors+29)
		Error("xgraph: Error in dimension of color table");

	bluemin=0;
	bluemax=nmax-1;
	colmin=0;
	colmax=nmax-1;

	  for(i=1;i<=icol;i++) {
		if(i<=10) {
			il=11-i;
			col.blue=bluemax;
			colora=(colmax*il)/10L;
			col.red=(colora*8L)/10L;
			col.green=colora;
		} else if(i<=19){
			il=i-9;
			col.blue=bluemax;
			col.red=colmin+((colmax-colmin)*il)/10L;
			col.green=colmin;
		} else {
			il=29-i;
			col.blue=bluemin+((bluemax-bluemin)*il)/10L;
			col.red=colmax;
			col.green=colmin;
		}
		result=XAllocColor(MyDisplay,MyColormap,&col);

		if(result) {
			MyColors[MaxColors+i]=col;
			VelShade[i-1]=MaxColors+i;
		} else {
			printf("Cannot allocate color %d\n",i);
			exit(1);
		}
	}
	MaxColors += 29;
}

long int linear(long int start , long int end , long int steps , long int act )

{
	return start + ( ( end - start ) * act ) / steps ;
}

long int quadratic_fast
			(
			  long int start 
			, long int end 
			, long int steps 
			, long int deriv
			, long int act 
			)

{
	long int dy,derivative,aux;

	dy = end - start;
	derivative = ( dy * deriv > 0 ) ? deriv : -deriv;
	aux = steps * derivative;
	return start + ((2*dy-aux)*act)/steps 
			+ ((aux-dy)*act*act)/(steps*steps);
}

long int quadratic_slow
			(
			  long int start 
			, long int end 
			, long int steps 
			, long int deriv
			, long int act 
			)

{
	long int dy,derivative,aux;

	dy = end - start;
	derivative = ( dy * deriv > 0 ) ? deriv : -deriv;
	aux = derivative;
	return start + aux*act + ((dy-aux*steps)*act*act)/(steps*steps);
}


static void QAllocGreen2BlueColors( void )

/* allocate 16 color tones from blue to colors */

{
	XColor color;
	Status result;
	int i,ncol;
	long unsigned int il,iaux,nmax=65536,auxmax,derv;
	long unsigned int colmin,colmax;

	ncol = 16;

	if(DimColors<MaxColors+ncol)
		Error("xgraph: Error in dimension of color table");

	iaux=ncol;
	colmin=0;
	colmax=nmax-1;
	auxmax=(3L*colmax)/4L;
	derv=2;

	for(i=0;i<ncol;i++) {
		if(i<ncol) {
			il=i;
			color.red=colmin;
			color.green=quadratic_slow(auxmax,colmin,iaux,derv,il);
			color.blue=quadratic_fast(colmin,auxmax,iaux,derv,il);
		}
		result=XAllocColor(MyDisplay,MyColormap,&color);

		if(result) {
			MaxColors++;
			MyColors[MaxColors]=color;
			GreenShade[i]=MaxColors;
		} else {
			printf("Cannot allocate color %d\n",MaxColors);
			exit(1);
		}
	}
}


static void QAllocYellow2GreenColors( void )

/* allocate 16 color tones from yellow to green */

{
	XColor color;
	Status result;
	int i,ncol;
	long unsigned int il,iaux,nmax=65536,auxmax;
	long unsigned int colmin,colmax;

	ncol = 16;

	if(DimColors<MaxColors+ncol)
		Error("xgraph: Error in dimension of color table");

	iaux=ncol;
	colmin=0;
	colmax=nmax-1;
	auxmax=(3L*colmax)/4L;

	for(i=0;i<ncol;i++) {
		if(i<ncol) {
			il=i;
			color.blue=colmin;
			color.red=linear(colmax,colmin,iaux,il);
			color.green=linear(colmax,auxmax,iaux,il);
		}
		result=XAllocColor(MyDisplay,MyColormap,&color);

		if(result) {
			MaxColors++;
			MyColors[MaxColors]=color;
			YellowShade[i]=MaxColors;
		} else {
			printf("Cannot allocate color %d\n",MaxColors);
			exit(1);
		}
	}
}

static void QAllocRed2YellowColors( void )

/* allocate 16 color tones from red to yellow */

{
	XColor color;
	Status result;
	int i,ncol;
	long unsigned int il,iaux,nmax=65536;
	long unsigned int colmin,colmax;

	ncol = 16;

	if(DimColors<MaxColors+ncol)
		Error("xgraph: Error in dimension of color table");

	iaux=ncol;
	colmin=0;
	colmax=nmax-1;

	for(i=0;i<ncol;i++) {
		if(i<ncol) {
			il=i;
			color.blue=colmin;
			color.red=colmax;
			color.green=linear(colmin,colmax,iaux,il);
		}
		result=XAllocColor(MyDisplay,MyColormap,&color);

		if(result) {
			MaxColors++;
			MyColors[MaxColors]=color;
			RedShade[i]=MaxColors;
		} else {
			printf("Cannot allocate color %d\n",MaxColors);
			exit(1);
		}
	}
}

int QAllocNewColor	( 
			long unsigned int red , 
			long unsigned int green ,
			long unsigned int blue 
			)

{
	XColor color;
	Status result;

	if(DimColors<MaxColors+1)
		Error("xgraph: Error in dimension of color table");

	color.red = red;
	color.green = green;
	color.blue = blue;

	result=XAllocColor(MyDisplay,MyColormap,&color);

	if(result) {
		MaxColors++;
		MyColors[MaxColors]=color;
	} else {
		printf("Cannot allocate color %d\n",MaxColors);
		exit(1);
	}

	return (result ? MaxColors : -1 );
}

static void QAllocBlueColors( void )

/* allocate 16 blue tones */

{
	XColor color;
	Status result;
	int i,ncol;
	long unsigned int il,iaux,nmax=65536;
	long unsigned int auxmax,colmin,colmax;
	long unsigned int bluemax,bluemin;

	bluemax=0; bluemin=0;

	ncol = 16;

	if(DimColors<MaxColors+ncol)
		Error("xgraph: Error in dimension of color table");

	iaux=8;
	auxmax=nmax-1;
	colmin=nmax/4;
	colmax=(nmax*3L)/4L;

	for(i=0;i<ncol;i++) {
		if(i<8) {
			il=i;
			color.blue=auxmax;
			color.red=colmin+(colmax-colmin)*(7L-il)/7L;
			color.green=colmin+(colmax-colmin)*(7L-il)/7L;
			color.red=linear(colmax,colmin,iaux,il);
			color.green=linear(colmax,colmin,iaux,il);
		} else {
			il=i-8;
			color.blue=bluemin+(bluemax-bluemin)*(7L-il)/7L;
			color.blue=linear(auxmax,colmin,iaux,il);
			color.red=colmin;
			color.green=colmin;
		}
		result=XAllocColor(MyDisplay,MyColormap,&color);

		if(result) {
			MaxColors++;
			MyColors[MaxColors]=color;
			BlueShade[i]=MaxColors;
		} else {
			printf("Cannot allocate color %d\n",MaxColors);
			exit(1);
		}
	}
}

int VelShadeColor( int shade )

{
	return VelShade[shade];
}

int BlueShadeColor( int shade )

{
	return BlueShade[shade];
}

int RedShadeColor( int shade )

{
	return RedShade[shade];
}

int YellowShadeColor( int shade )

{
	return YellowShade[shade];
}

int GreenShadeColor( int shade )

{
	return GreenShade[shade];
}


void QSetGeometry( char *s )

{
	int x,y;
	int width,height;	/* unsigned deleted 21.03.94 */

	(void) XGeometry( MyDisplay , MyScreen , s , NULL ,
				MyBorder , 1 , 1 , 0 , 0 ,
				&x , &y , &width , &height );

	XMinPix = x;
	YMinPix = y;
	XMaxPix = x+width-1;
	YMaxPix = y+height-1;
}

void QSetDisplay( char *s )

{
	if( s != NULL )
		MyDisplayName = savestr(s);
}

void QSetWidget( Widget w )

{
	MyWidget = w;
}

void QSetTitle( char *s )

{
	if( s != NULL )
		MyTitle = savestr(s);
}


static char *savestr(char *s )

{
	int len;
	char *p;

	len = strlen(s);

	p = (char *) malloc( len+1 );

	if(p == NULL) { 
		printf("savest : No memory left to allocate string");
		exit(3);
	}

	p = strncpy(p,s,len); 
	p[len] = '\0';

	return(p);
}



#endif

