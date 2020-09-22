
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1998,2011  Georg Umgiesser
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
 * 21.03.1994	ggu	gcc-warnings, call to XGeometry corrected (unsigned)
 * ...		ggu	VelColors used for call to QAllocVelColors(),
 * ...		ggu	*ShadeColor() declared in xgraph.h,
 * ...		ggu	gustd.h not included anymore
 * 15.02.1995	ggu	bug fix in QAllocColor: must return MaxColors-1
 * 01.09.1995	ggu	QGetViewport, QGetWindow, QBell routines added
 * 05.12.1995	ggu	changes (Widget, Cursor, ...) included
 * 07.07.1998	ggu	can use without Xt library
 * 19.04.2011	ggu	introduced HAVE_WIDGET for compiler errors in Mac
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "xgraph.h"

#include "generalx.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>


#define QxWtoV(x) (SxWtoV*((x)-WinXmin)+ViewLeft+0.5)
#define QyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom+0.5)
#define QxVtoW(x) (SxVtoW*((x)-ViewLeft)+WinXmin)
#define QyVtoW(y) (SyVtoW*((y)-ViewBottom)+WinYmin)


#define HAVE_XT	0		/* 1 if Xt lib is available */
#define HAVE_WIDGET	1	/* 1 if Widget is available */

/* definitions for colors */

#define MAXCOLOR 256

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

static int MaxColors=0;		/* total number of colors defined */

static void QInitColors( void );

static char *MyColorNames[] =
		{"black","blue","green","cyan","red","magenta",
		 "brown","LightGray","DimGray",
		 "LightBlue","LightSeaGreen","LightCyan",
		 "IndianRed","LightPink","yellow","white"
		};

static char *Qsavestr(char *s );


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

static int           Initialized=0;

static Display       *MyDisplay;
static int           MyScreen;
static Window        MyRootWindow;
static Window        MyWindow;
static GC            MyGc;
static unsigned int  MyDepth;

#if HAVE_WIDGET == 0
typedef void *Widget;
Display *XtDisplay(Widget MyWidget);
Window XtWindow(Widget MyWidget);
#endif

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
static XColor        MyColors[MAXCOLOR]; /* error fixed 1.7.93 */
static unsigned int  MyLineWidth=0; /* set to 1 for pixel reproducibility */
static int           MyFontAscent;
static int           MyFontDescent;
static XCharStruct   MyCharStruct;
static char          *MyDisplayName=NULL;
static Pixmap        MyIcon=None;

static char          *MyTitle=NULL;
static char          *MyCommand[] = {"main"};

static int           TEST=0;
static char          *NullName = "";

/******************************************************************/
/********************* routine dependent part *********************/
/******************************************************************/

/*\
 *	if Xt library not available, comment out call
 *	  to XtDisplay() and XtWindow() and link only
 *	  with X11
 *	  -> new: set HAVE_XT to 0
 *	comment out grid_width, grid_height, grid_bits[]
 *	  if no icon is given
 *	in MyText is title of window
 *	if Server shuts down when routine is executed,
 *	  try to remove all references to MyIcon
\*/

/* next is title of window */

/*static char          MyText[] = {"Georg Umgiesser - ISDGM/CNR"};*/
static char          MyText[] = {"grid"};

/* next is icon of application */

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

/*****************************************************************/

void QGraphInit( void )

{
    int dummy;

    XSizeHints myhint;
    XWMHints   mywmhint;
    XWindowAttributes   myattrib;

    if(MyDisplayName == NULL)
	MyDisplayName = NullName;
    if(MyTitle == NULL)
	MyTitle = MyText;

    if(MyWidget)
        MyDisplay  = XtDisplay(MyWidget);
    else
        MyDisplay  = XOpenDisplay(MyDisplayName);

    MyScreen   = DefaultScreen(MyDisplay);
    MyRootWindow = DefaultRootWindow(MyDisplay);
    MyDepth = XDefaultDepth(MyDisplay,MyScreen);

    MyBackGround  = WhitePixel(MyDisplay,MyScreen);
    MyForeGround  = BlackPixel(MyDisplay,MyScreen);
    MyBorderColor = MyForeGround;

#ifdef grid_width
    MyIcon = XCreatePixmapFromBitmapData(
                        MyDisplay,MyRootWindow,
                        (char *) grid_bits,grid_width,grid_height,
                        MyForeGround,MyBackGround,1);
#endif

    myhint.x = XMinPix;
    myhint.y = YMinPix;
    myhint.width  = XMaxPix - XMinPix + 1;
    myhint.height = YMaxPix - YMinPix + 1;
    myhint.flags  = PPosition | PSize;

    mywmhint.flags = InputHint | IconPixmapHint;
    mywmhint.input = True;
    mywmhint.icon_pixmap = MyIcon;

    if(MyWidget)
        MyWindow = XtWindow(MyWidget);
    else {
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

    QInitColors();

    if( !MyWidget ) {
        XMapRaised(MyDisplay,MyWindow);
        XSync(MyDisplay,0); /* wait till window is raised and ready */
    }

    Initialized = 1;

    /* Viewport is total window */

    XGetWindowAttributes( MyDisplay , MyWindow , &myattrib );

    QViewport( 0 , 0 , myattrib.width , myattrib.height );
    QWindow( 0. , 0. , 1. , 1. );
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
	Color = color%MaxColors;
	XSetForeground(MyDisplay,MyGc,MyColors[Color].pixel);
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

void QGetViewport( int *left , int *top , int *right , int *bottom )

{
    *left   = ViewLeft;
    *top    = ViewTop;
    *right  = ViewRight;
    *bottom = ViewBottom;
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

void QGetWindow( float *xmin , float *ymin , float *xmax , float *ymax )

{
    *xmin = WinXmin;
    *ymin = WinYmin;
    *xmax = WinXmax;
    *ymax = WinYmax;
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

void QTextBackground( int color )

{
}

void QTextDimensions( char *s , float *width , float *height )

{
    *width  = SxVtoW * XTextWidth(MyFontStruct,s,strlen(s));
    *height = - ( SyVtoW * MyFontAscent );
}

void QTextDimensionsI( char *s , int *width , int *height )

{
    *width  = XTextWidth(MyFontStruct,s,strlen(s));
    *height = MyFontAscent;
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

void QScreenMaxXY( int *xmax , int *ymax )

{
        *xmax = DisplayWidth(MyDisplay,MyScreen);
        *ymax = DisplayHeight(MyDisplay,MyScreen);
}

void QWindowMaxXY( int *xmax , int *ymax )

{
	XWindowAttributes wa;

	XGetWindowAttributes(MyDisplay,MyWindow,&wa);
        *xmax = wa.width;
        *ymax = wa.height;
}

void QBell( void )

{
/*
        XKeyboardState *kbs;

        XGetKeyboardControl(MyDisplay,kbs);
        XBell(MyDisplay,kbs->bell->percent);
*/
        XBell(MyDisplay,0);
}

void *QSavePixels( int x , int y , int width , int height )

{
        Pixmap_type *bp;
	Pixmap q;
	int w,h;

        bp = (Pixmap_type *) malloc( sizeof(Pixmap_type) );
        if( !bp ) return NULL;

	q = XCreatePixmap(MyDisplay,MyWindow,
		(unsigned int)width,
		(unsigned int)height,
		MyDepth);

	if( !q ) {
		free(bp);
		return NULL;
	}

        bp->xbuffer = (int) q;
        bp->width = width;
        bp->height = height;

	QWindowMaxXY(&w,&h);
	QViewport(0,0,w,h);

	XCopyArea(MyDisplay,MyWindow,q,MyGc,x,y,
		(unsigned int)width,
		(unsigned int)height,
		0,0);

	return (void *) bp;
}


void QRestorePixels( int x , int y , void *buffer )

{
        Pixmap_type *bp;
        Pixmap q;
        int width,height;

        if( !buffer ) return;

        bp = (Pixmap_type *) buffer;

        width = bp->width;
        height = bp->height;
        q = (Pixmap) bp->xbuffer;

	XCopyArea(MyDisplay,q,MyWindow,MyGc,0,0,
		(unsigned int)width,
		(unsigned int)height,
		x,y);

	XFreePixmap(MyDisplay,q);
	free(bp);
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

/*
void QNextEvent( XEvent *eventp )

{
	XNextEvent(MyDisplay,eventp);
}
*/

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

int QAllocColor( int red , int green , int blue )

	/*    r/g/b in [0...255]    */

{
	static int mult=256;
	XColor col;

	col.red   = red*mult;
	col.green = green*mult;
	col.blue  = blue*mult;

	if( MaxColors >= MAXCOLOR ) {
		Error("No free colors available");
	} else if( !MyColormap ) {
		Error("No colormap available");
	} else if( XAllocColor(MyDisplay,MyColormap,&col) ) {
		MyColors[MaxColors++]=col;
	} else {
		Error("Error allocating color");
	}

	return MaxColors-1;
}

static void QInitColors( void )

{
	XColor color,exact;
	Status result;
	int i;

	if( MAXCOLOR < 16 ) 
		Error("xgraph: Error in dimension of color table");

	MyColormap = XDefaultColormap(MyDisplay,MyScreen);

	for(i=0;i<16;i++) {
		result=XAllocNamedColor(MyDisplay,MyColormap,MyColorNames[i]
			,&color,&exact);

		if(result)
			MyColors[i]=color;
		else {
			Error("QInitColors : Error allocating color");
		}
	}

	MaxColors = 16;
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
		MyDisplayName = Qsavestr(s);
}

void QSetWidget( Widget w )

{
        MyWidget = w;
}

void QSetTitle( char *s )

{
	if( s != NULL )
		MyTitle = Qsavestr(s);
}


static char *Qsavestr( char *s )

{
	int len;
	char *p;

	len = strlen(s);

	p = (char *) malloc( len+1 );

	if(p == NULL) { 
		printf("Qsavestr : No memory left to allocate string (%d)",len);
		exit(3);
	}

	p = strncpy(p,s,len); 
	p[len] = '\0';

	return(p);
}

#if HAVE_XT == 0

Display *XtDisplay(Widget MyWidget)

{
	return NULL;
}

Window XtWindow(Widget MyWidget)

{
	return 0;
}

#endif
