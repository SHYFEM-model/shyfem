
/************************************************************************\ 
 *									*
 * pgraph.c - graphic routines for postscript output			*
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
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUC_PGRAPH_
#define __GUC_PGRAPH_

/*
#include <X11/Xlib.h>
#include <X11/Xutil.h>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pgraph.h"


/* definitions for colors */

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

/*
static char *savestr(char *s );
*/

static int MaxColors=15;
static long int NPlot=0;

/*****************************************************************/

/* absolute screen coordinates of window */

/* probably not used ... */

static int XMinPix    = 100;
static int YMinPix    = 100;
static int XMaxPix    = 100+640-1;
static int YMaxPix    = 100+480-1;

/* viewport coordinates relative to window - (left,top) = (0,0) */

static float ViewLeft   = 0.;
static float ViewTop    = 0.;
static float ViewRight  = 0.;
static float ViewBottom = 0.;

static float WinXmin  = 0.;
static float WinYmin  = 0.;
static float WinXmax  = 0.;
static float WinYmax  = 0.;

static float SxWtoV   = 0.;
static float SyWtoV   = 0.;
static float SxVtoW   = 0.;
static float SyVtoW   = 0.;

static int Color      = 15;
/*
static int TextHeight = 0;
*/

static float           XAct=0.;
static float           YAct=0.;

/*
static char          MyText[] = {"Georg Umgiesser - ISDGM/CNR"};
static char          *MyCommand[] = {"main"};
*/

static int           TEST=0;
/*
static char          *NullName = "";
*/

static FILE *FP;

void QPageHeader( void );

/* 
	scale used : 28.3

	corresponds to 1. cm

	28.3 points =^= 1.0 cm
	72   points =^= 1.0 inch
   ==>	1.0 inch = (72/28.3) cm
 */

/*
 *	to do : setlinewidth, rotate, dash,...
 */

/*****************************************************************/

void QGraphInit( void )

{
    FILE *fp;
    static char sinit[] = {"%!\n"};

    fp=fopen("plot.ps","w");
    if(fp) {
	FP = fp;
    } else {
	Error2("Cannot open plot file : ","plot.ps");
    }

/*                %!           former... */
    fprintf(FP,"%s",sinit);
    fprintf(FP,"%% PostScript emulator for Calcomp and Graphics\n");
    fprintf(FP,"%% written by Georg Umgiesser, ISDGM/CNR\n");
    fprintf(FP,"%% 1364 S. Polo, 30125 Venezia, Italy\n");
    fprintf(FP,"%% \n");

	/* unit in centimeters */

    QPageHeader();
    QTextRealSize(1.);

    QViewport(0.0,0.0,19.,28.0);
    QWindow(0.0,0.0,19.,28.0);

	printf("%f %f %f %f\n",SxWtoV,SyWtoV,SxVtoW,SyVtoW);
	printf("%f %f %f %f\n",ViewLeft,ViewBottom,ViewRight,ViewTop);
	printf("%f %f %f %f\n",WinXmin,WinYmin,WinXmax,WinYmax);

    NPlot=0;
}

void QPageHeader( void )

{
    fprintf(FP,"28.3 28.3 scale\n");
    fprintf(FP,"0.02 setlinewidth\n");
    fprintf(FP,"0.0 0.0 translate\n");
}

void QClearScreen( void )

{
    fprintf(FP,"showpage\n");
    NPlot=0;

    QPageHeader();
    QTextRealSize(1.);
}

void QGraphClose( void )

{
    fprintf(FP,"showpage\n");
    fprintf(FP,"%%%%Trailer\n"); /* needed to conform to DSC */
    fclose(FP);
    FP=NULL;

    if( !NPlot )
	printf("empty last page\n");
}

void QNewPen( int color )

{
	if(color>=0 && color <= MaxColors ) {
		Color=color;
		fprintf(FP,"%f setgray\n",1.-((float) color)/MaxColors);
		NPlot++;
	}
}

void QGetPen( int *color )

{
	*color=Color;
}

void QGetPixelWidth( float *dx , float *dy )

{
	*dx = QxVtoW(1.)-QxVtoW(0.);
	*dy = QyVtoW(0.)-QyVtoW(1.);
}

void QViewport( float left , float bottom , float right , float top )

{
    ViewLeft   = left;
    ViewTop    = top;
    ViewRight  = right;
    ViewBottom = bottom;

    QAdjustScale();
}

void QClearViewport( void )

{
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
    float vx1,vy1,vx2,vy2;

    vx1=QxWtoV(x1);
    vy1=QyWtoV(y1);
    vx2=QxWtoV(x2);
    vy2=QyWtoV(y2);

    if(TEST) printf("QLine: %f %f %f %f %d\n",vx1,vy1,vx2,vy2,Color);

    fprintf(FP,"%f %f moveto\n",vx1,vy1);
    fprintf(FP,"%f %f lineto\n",vx2,vy2);
    fprintf(FP,"stroke\n");

    XAct=vx2;
    YAct=vy2;
    NPlot++;
}

void QPlot( float x , float y )

{
    float vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    if(TEST) printf("QPlot: %f %f %f %f %d\n",XAct,YAct,vx,vy,Color);

    fprintf(FP,"%f %f moveto\n",XAct,YAct);
    fprintf(FP,"%f %f lineto\n",vx,vy);
    fprintf(FP,"stroke\n");

    XAct=vx;
    YAct=vy;
    NPlot++;
}

void QMove( float x , float y )

{
    float vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    if(TEST) printf("QMove: %f %f %d\n",vx,vy,Color);

    fprintf(FP,"%f %f moveto\n",vx,vy);

    XAct=vx;
    YAct=vy;
    NPlot++;
}

void QPoint( float x , float y )

{
    static float dd=0.05;
    float vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    if(TEST) printf("QPoint: %f %f %d\n",vx,vy,Color);

    QRectFill(x,y,x+dd,y+dd,Color);
    NPlot++;
}

void QAreaFill( int ndim , float *x , float *y )

{
    int i;
    float vx,vy;

    if(ndim<3) return;

    fprintf(FP,"newpath\n");
    vx = QxWtoV(x[0]);
    vy = QyWtoV(y[0]);
    fprintf(FP,"%f %f moveto\n",vx,vy);
    for(i=1;i<ndim;i++) {
	vx = QxWtoV(x[i]);
	vy = QyWtoV(y[i]);
	fprintf(FP,"%f %f lineto\n",vx,vy);
    }
    fprintf(FP,"closepath\n");
    fprintf(FP,"fill\n");
    NPlot++;
}

void QRectFill( float x1 , float y1 , float x2 , float y2 , int color )

{
	float ix1,iy1,ix2,iy2;
	float ix0,iy0;
	float ixw,iyw;
	float x[4],y[4];
	int tmpcol=0;

	ix1 =  x1 ;
	iy1 =  y1 ;
	ix2 =  x2 ;
	iy2 =  y2 ;

	if( color >= 0 ) {
		tmpcol=Color;
		QNewPen(color);
	}

	if( ix1 < ix2 ) {
		ix0 = ix1;
		ixw = ix2-ix1;
	} else {
		ix0 = ix2;
		ixw = ix1-ix2;
	}

	if( iy1 < iy2 ) {
		iy0 = iy1;
		iyw = iy2-iy1;
	} else {
		iy0 = iy2;
		iyw = iy1-iy2;
	}

	x[0]=ix0;y[0]=iy0;
	x[1]=ix0+ixw;y[1]=iy0;
	x[2]=ix0+ixw;y[2]=iy0+iyw;
	x[3]=ix0;y[3]=iy0+iyw;

	QAreaFill(4,x,y);

	if( color >= 0 )
		QNewPen(tmpcol);
        NPlot++;
}

void QTextPointSize( int size )

/* size is in points */

/* factor for some unknown reason */
#define TFACTOR (5.1/3.)

{
    if( size > 0 ) {
	fprintf(FP,"/Courier findfont\n");
	fprintf(FP,"%f scalefont\n",(TFACTOR*size)/28.3);
	fprintf(FP,"setfont\n");
        NPlot++;
    }
}

void QTextRealSize( float size )

/* size is in actual units */

{
    if( size > 0. ) {
	fprintf(FP,"/Courier findfont\n");
	fprintf(FP,"%f scalefont\n",TFACTOR*size);
	fprintf(FP,"setfont\n");
        NPlot++;
    }
}

#undef TFACTOR

void QText( float x , float y , char *s )

{
    float vx,vy;

    vx=QxWtoV(x);
    vy=QyWtoV(y);

    fprintf(FP,"%f %f moveto\n",vx,vy);
    fprintf(FP,"(%s) show\n",s);
    NPlot++;
}

void QTextDimensions( char *s , float *width , float *height )

{
/*
    *width  = SxVtoW * XTextWidth(MyFontStruct,s,strlen(s));
    *height = - ( SyVtoW * MyFontAscent );
*/
}

void QRealXY( float vx , float vy , float *x , float *y )

{
	*x = QxVtoW( vx );
	*y = QyVtoW( vy );
}

void QScreenXY( float x , float y , float *vx , float *vy )

{
	*vx = QxWtoV( x );
	*vy = QyWtoV( y );
}

/********************************************************************/

void QFlush( void )

{
}

void QSync( void )

{
}

void QSyncron( int sync )

{
}

int QActScreen( void )

{
	return 0;
}


void QSetTest( void )

{
	TEST = 1;
}

void QClearTest( void )

{
	TEST = 0;
}

/*

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

*/

#endif

