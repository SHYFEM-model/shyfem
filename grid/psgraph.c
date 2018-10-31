
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
 * psgraph.c - graphic routines for postscript output			*
 *									*
 * Revision History:							*
 * 24-Sep-2004: new routine PsFlush                                     *
 * 26-Apr-2001: new routine PsSetPointSize                              *
 * 12-Sep-97: PsAdjustScale becomes static, PsWindow -> PsSetWorld      *
 * 11-Sep-97: minor modifications (PsGraphOpen,...)                     *
 * 03-Jun-97: PsRectFill() optimized, also PsSetColor() for Hue         *
 * 09-May-97: append feature in PsText()                                *
 * 03-May-97: restructured, PsRectFill() new calling arguments,         *
 *            no integer color any more                                 *
 * 07-Dec-95: PPageHeader() now declared and defined static (bug)       *
 * 21-Oct-94: PSetClipWindow does clipping,                             *
 *            translate whole plot area 0.5 cm into the paper           *
 *            change names from Q to P                                  *
 *            rectify scale with PRectifyScale (Aspect ratios are       *
 *              not distorted anymore                                   *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"
#include "gustd.h"
#include "psgraph.h"


#define PsxWtoV(x) (SxWtoV*((x)-WinXmin)+ViewLeft)
#define PsyWtoV(y) (SyWtoV*((y)-WinYmin)+ViewBottom)
#define PsxVtoW(x) (SxVtoW*((x)-ViewLeft)+WinXmin)
#define PsyVtoW(y) (SyVtoW*((y)-ViewBottom)+WinYmin)

#define InSide(x,y) ( (x) >=ViewLeft && (x) <= ViewRight && \
			(y) >= ViewBottom && (y) <= ViewTop )

#define MINMAX(c) 					\
			if( (c) < 0. ) {		\
				(c) = 0.;		\
			} else if( (c) > 1. ) {		\
				(c) = 1.;		\
			}

#define SETCOLOR	if( ColorChanged == YES ) PsSetColor()
#define SETTEXT		if(  FontChanged == YES ) PsSetFont()
#define ISWHITE		( IsWhite && !BPaintWhite )
#define MUSTROTATE	( FontRotate != 0. )

#define CMODE_GRAY	1
#define CMODE_HSB	2
#define CMODE_RGB	4

#define FLUSH	0
#define APPEND	1

#define	LINEWIDTH	0.01
#define	POINTSIZE	0.01
#define	FONTSIZE	12

#define TFACTOR (5.1/3.)	/* text factor for unknown reason */

static long int NPlot       = 0;	/* number of plotted items */
static long int NPath       = 0;	/* subpaths in actual path */
static      int NPage       = 0;	/* actual page */
static	    int NClip       = 0;	/* number of page clipping */
static	    int BOutSide    = NO;	/* actual point is outside ? */
static	    int BPageOpened = NO;	/* page has been opened for plot ? */
static	    int BPaintWhite = NO;	/* paint with white ? (only for gray) */

/*****************************************************************/


/* viewport coordinates relative to page - (left,bottom) = (0,0) */

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

/*
static int	LineCap		= 1;
static int	LineJoin	= 1;
*/
static float	LineWidth	= LINEWIDTH;
static float	PointSize	= POINTSIZE;	/* size for single point */

static char	Font[81]	= "Courier";
static float	FontSize	= FONTSIZE/28.3;
static float	FontRotate	= 0.;
static int	FontChanged	= NO;

static float	C1		= 0.;		/* gray, hue, red */
static float	C2		= 0.;		/* sat, green */
static float	C3		= 0.;		/* bri, blue */
static int	ColorMode	= CMODE_GRAY;
static int	ColorChanged	= NO;
static int	IsWhite		= NO;

static float	XAct		= 0.;
static float	YAct		= 0.;

static int	TEST		= 0;

static FILE	*FP		= NULL;

static void PsPageHeader( void );
static void PsReset( void );
static void PsSetClipWindow( void );
static void PsStroke( int mode );
static void PsSetColor( void );
static void PsSetFont( void );

/* 
	scale used : 28.3

	corresponds to 1. cm

	28.3 points =^= 1.0 cm
	72   points =^= 1.0 inch
   ==>	1.0 inch = (72/28.3) cm = 2.54 cm
 */
	/* unit in centimeters */

/*
	every page resets parameters to initital
	every call to SetViewport calls PsReset, so
	changed parameters on one page are always
	conserved
*/
/*
 *	to do : setlinewidth, rotate, dash,...
 */

/*****************************************************************/
/****************** Static Routines ******************************/
/*****************************************************************/

static void PsPageHeader( void )

{
    NPage++;
    NPlot = 0;
    NPath = 0;
    NClip = 0;
    BOutSide = NO;
    BPaintWhite = NO;

    fprintf(FP,"%%%%Page: %d %d\n",NPage,NPage);
    fprintf(FP,"%%%%BeginPageSetup\n");
    fprintf(FP,"GS\n");
/*    fprintf(FP,"initmatrix 28.3 dup scale\n"); */
    fprintf(FP,"28.3 dup scale\n");
    fprintf(FP,"1.0 1.0 translate\n");
    fprintf(FP,"1 LC 1 LJ\n");

    PsPaintWhite(NO);
    PsSetGray(0.);
    LineWidth = LINEWIDTH;
    PointSize = POINTSIZE;
    FontRotate = 0.0;
    PsResetDashPattern();
    PsTextPointSize(FONTSIZE);
    PsTextFont("Courier");

    PsSetViewport(0.0,0.0,19.,28.0);
    PsSetWorld(0.0,0.0,19.,28.0);

    fprintf(FP,"%%%%EndPageSetup\n");
}

static void PsReset( void )

{
	PsSetColor();
	PsSetLineWidth(LineWidth);
	PsSetFont();
}

static void PsSetClipWindow( void )

{
	float x,y,width,height;

	x = ViewLeft;
	y = ViewBottom;
	width = ViewRight - x;
	height = ViewTop - y;

	PsStroke(FLUSH);
	fprintf(FP,"%% new clipping window\n");
	if( NClip ) {
		fprintf(FP,"GR GS\n");
	} else {
		fprintf(FP,"GS\n");
	}
	PsReset();	/* re-establish parameters from before */
	fprintf(FP,"N %1.3f %1.3f M\n",x,y);
	fprintf(FP,"%1.3f %1.3f L\n",x+width,y);
	fprintf(FP,"%1.3f %1.3f L\n",x+width,y+height);
	fprintf(FP,"%1.3f %1.3f L\n",x,y+height);
	fprintf(FP,"C clip N\n");
	NClip++;
}

static void PsStroke( int mode )

/*
    mode == FLUSH  -> flush buffer
    mode == APPEND -> increment counter
 */

{
	if( mode == FLUSH && NPath != 0 ) {
		fprintf(FP,"S\n");
		NPath=0;
	} else if( mode == APPEND ) {
		NPath++;
		if( NPath > 100 ) {
			fprintf(FP,"S\n");
			NPath=0;
		}
	}
}
			
/*****************************************************************/
/******************* Administration Routines *********************/
/*****************************************************************/

void PsGraphInit( char *file )

{
    FILE *fp;
    static char sinit[] = {"%!PS-Adobe-3.0\n"};

    if( file && *file != '\0' ) {
      fp=fopen(file,"w");
      if(fp) {
	FP = fp;
      } else {
	Error2("Cannot open plot file : ",file);
      }
    } else {
      FP = stdout;
    }

    fprintf(FP,"%s",sinit);
    fprintf(FP,"%%%%BoundingBox: 10 10 588 824\n"); /* A4 paper 21x29.5cm */
    fprintf(FP,"%%%%Pages: (atend)\n");
    fprintf(FP,"%%%%Title: PostScript Graphics 1.2 - ");
    fprintf(FP,"(c) 1994 Georg Umgiesser\n");
    fprintf(FP,"%%%%Creator: psgraph\n");
    fprintf(FP,"%%%%CreationDate: %s\n",GetTime());
    fprintf(FP,"%%%%EndComments\n");
    fprintf(FP,"%% ***********************************************\n");
    fprintf(FP,"%% * PostScript Graphics 1.2                     *\n");
    fprintf(FP,"%% * Copyright (c) 1994-1997 Georg Umgiesser     *\n");
    fprintf(FP,"%% * written by                                  *\n");
    fprintf(FP,"%% * Georg Umgiesser, ISDGM/CNR                  *\n");
    fprintf(FP,"%% * 1364 S. Polo, 30125 Venezia, Italy          *\n");
    fprintf(FP,"%% * Tel.   : ++39-41-5216875                    *\n");
    fprintf(FP,"%% * Fax    : ++39-41-2602340                    *\n");
    fprintf(FP,"%% * E-Mail : georg@lagoon.isdgm.ve.cnr.it       *\n");
    fprintf(FP,"%% ***********************************************\n");

    fprintf(FP,"%%%%BeginProlog\n");
    /*fprintf(FP,"%%%%BeginResource: procset all\n");*/
    fprintf(FP,"%%%%BeginResource: general_procset\n");
    fprintf(FP,"/M /moveto load def\n");
    fprintf(FP,"/L /lineto load def\n");
    fprintf(FP,"/RM /rmoveto load def\n");
    fprintf(FP,"/RL /rlineto load def\n");
    fprintf(FP,"/S /stroke load def\n");
    fprintf(FP,"/C /closepath load def\n");
    fprintf(FP,"/F /fill load def\n");
    fprintf(FP,"/G /setgray load def\n");
    fprintf(FP,"/H /sethsbcolor load def\n");
    fprintf(FP,"/HH { 1 1 sethsbcolor } bind def\n");
    fprintf(FP,"/R /setrgbcolor load def\n");
    fprintf(FP,"/N /newpath load def\n");
    fprintf(FP,"/LC /setlinecap load def\n");
    fprintf(FP,"/LJ /setlinejoin load def\n");
    fprintf(FP,"/LW /setlinewidth load def\n");
    fprintf(FP,"/GS /gsave load def\n");
    fprintf(FP,"/GR /grestore load def\n");
    fprintf(FP,"/FF /findfont load def\n");
    fprintf(FP,"/CF /scalefont load def\n");
    fprintf(FP,"/SF /setfont load def\n");
    fprintf(FP,"/RF { M dup 0 RL exch 0 exch RL neg 0 RL C F } bind def\n");
    fprintf(FP,"/TRSF { /Times-Roman FF exch CF SF } bind def %% points\n");
    fprintf(FP,"/AW { TRSF M show } bind def %% (text) x y points\n");
    fprintf(FP,"%%%%EndResource\n");
    fprintf(FP,"%%%%BeginResource: eps_procset\n");
    fprintf(FP,"/BeginEPSF {\n");
    fprintf(FP,"  /b4_Inc_state save def\n");
    fprintf(FP,"  /dict_count countdictstack def\n");
    fprintf(FP,"  /op_count count 1 sub def\n");
    fprintf(FP,"  userdict begin\n");
    fprintf(FP,"  /showpage { } def\n");
    fprintf(FP,"  0 setgray 0 setlinecap \n");
    fprintf(FP,"  1 setlinewidth 0 setlinejoin\n");
    fprintf(FP,"  10 setmiterlimit [ ] 0 setdash newpath\n");
    fprintf(FP,"  languagelevel where\n");
    fprintf(FP,"  { pop languagelevel\n");
    fprintf(FP,"  1 ne\n");
    fprintf(FP,"    {false setstrokeadjust false setoverprint\n");
    fprintf(FP,"    } if\n");
    fprintf(FP,"  } if\n");
    fprintf(FP,"} bind def\n");
    fprintf(FP,"/EndEPSF {\n");
    fprintf(FP,"  count op_count sub {pop} repeat\n");
    fprintf(FP,"  countdictstack dict_count sub {end} repeat\n");
    fprintf(FP,"  b4_Inc_state restore\n");
    fprintf(FP,"} bind def\n");
    fprintf(FP,"/Rect { %% llx lly w h\n");
    fprintf(FP,"  4 2 roll moveto 1 index 0 rlineto\n");
    fprintf(FP,"  0 exch rlineto neg 0 rlineto closepath\n");
    fprintf(FP,"} bind def\n");
    fprintf(FP,"/ClipRect { %% llx lly w h\n");
    fprintf(FP,"  Rect clip newpath\n");
    fprintf(FP,"} bind def\n");
    fprintf(FP,"%%%%EndResource\n");
    fprintf(FP,"%%%%EndProlog\n");
    fprintf(FP,"%%%%BeginSetup\n");
    fprintf(FP,"%%%%EndSetup\n");

    PsStartPage();
}

void PsGraphOpen( void )

{
	PsGraphInit("plot.ps");
}

void PsGraphClose( void )

{
    PsEndPage();

    fprintf(FP,"%%%%Trailer\n"); 
    fprintf(FP,"%%%%Pages: %d\n",NPage);
    fprintf(FP,"%%%%EOF\n"); 
    fclose(FP);
    FP=NULL;

    if( !NPlot )
	printf("PsGraphClose: empty last page\n");
}

void PsStartPage( void )

{
    if( !BPageOpened ) {
        BPageOpened = 1;
        PsPageHeader();
    }
}

void PsEndPage( void )

{
    if( BPageOpened ) {
        BPageOpened = 0;
        PsStroke(FLUSH);
	if( NClip ) fprintf(FP,"GR\n");
        fprintf(FP,"GR\n");
        fprintf(FP,"showpage\n");
    }
}

void PsNewPage( void )

{
    PsEndPage();
    PsStartPage();
}

/***********************************************************************/
/************************ Color Routines *******************************/
/***********************************************************************/

static void PsSetColor( void )

{
	PsStroke(FLUSH);
	if( ColorMode == CMODE_GRAY ) {
	  if( ! ISWHITE ) fprintf(FP,"%1.3f G\n",C1);
	} else if( ColorMode == CMODE_HSB ) {
	  if( C2 == 1. && C3 == 1. ) {
	    fprintf(FP,"%1.3f HH\n",C1);
	  } else {
	    fprintf(FP,"%1.3f %1.3f %1.3f H\n",C1,C2,C3);
	  }
	} else if( ColorMode == CMODE_RGB ) {
	  fprintf(FP,"%1.3f %1.3f %1.3f R\n",C1,C2,C3);
	} else {
	  Error2("Internal error PsSetColor: ",itos(ColorMode));
	}
	ColorChanged = NO;
}

void PsSetGray( float gray )

{
	MINMAX(gray);
	if( ColorMode == CMODE_GRAY && C1 == gray ) return;
	C1 = gray;
	ColorMode = CMODE_GRAY;
	ColorChanged = YES;
	if( gray >= 1. ) {
	  IsWhite = YES;
	} else {
	  IsWhite = NO;
	}
}

void PsSetHue( float hue )

{
	MINMAX(hue); 
	if( ColorMode == CMODE_HSB && C1 == hue ) return;
	if( ColorMode != CMODE_HSB ) {
	   C2 = 1.; C3 = 1.;
	}
	C1 = hue;
	ColorMode = CMODE_HSB;
	ColorChanged = YES;
	IsWhite = NO;
}

void PsSetHSB( float hue , float sat , float bri )

{
	MINMAX(hue); MINMAX(sat); MINMAX(bri);
	C1 = hue; C2 = sat; C3 = bri;
	ColorMode = CMODE_HSB;
	ColorChanged = YES;
	IsWhite = NO;
}

void PsSetRGB( float red , float green , float blue )

{
	MINMAX(red); MINMAX(green); MINMAX(blue);
	C1 = red; C2 = green; C3 = blue;
	ColorMode = CMODE_RGB;
	ColorChanged = YES;
	IsWhite = NO;
}

void PsPaintWhite( int ipaint )

{
	BPaintWhite = ipaint;
	ColorChanged = YES;
}

/***********************************************************************/
/**************************** Line properties **************************/
/***********************************************************************/

void PsSetDashPattern( float fact, float offset, int n, float *array )

{
	float val;

	PsStroke(FLUSH);
	fprintf(FP,"[");
	while( n-- ) {
	  val = *array++;
          val = SxWtoV * val * fact;
	  fprintf(FP,"%1.3f",val);
	  if( n > 0 ) fprintf(FP," ");
	}
	fprintf(FP,"] %1.3f setdash\n",offset);
}

void PsResetDashPattern( void )

{
	fprintf(FP,"[] 0 setdash\n");
}

void PsSetLineWidth( float width )

{
	if( width <= 0 ) width = LINEWIDTH;
	LineWidth = width;
	PsStroke(FLUSH);
	fprintf(FP,"%1.3f LW\n",width);
}

void PsSetPointSize( float size )

{
	if( size <= 0 ) size = POINTSIZE;
	PointSize = size;
}

/***********************************************************************/
/*********************** Clipping Routines *****************************/
/***********************************************************************/

static void PsAdjustScale( void )

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

void PsSetViewport( float left , float bottom , float right , float top )

{
    ViewLeft   = left;
    ViewTop    = top;
    ViewRight  = right;
    ViewBottom = bottom;

    PsAdjustScale();
    PsSetClipWindow();
}

void PsGetViewport( float *left , float *bottom , float *right , float *top )

{
    *left = ViewLeft;
    *top = ViewTop;
    *right = ViewRight;
    *bottom = ViewBottom;
}

void PsClearViewport( void )

{
    PsStroke(FLUSH);
    fprintf(FP,"GS 1 G N\n");
    fprintf(FP,"%1.3f %1.3f M\n",ViewLeft,ViewTop);
    fprintf(FP,"%1.3f %1.3f L\n",ViewLeft,ViewBottom);
    fprintf(FP,"%1.3f %1.3f L\n",ViewRight,ViewBottom);
    fprintf(FP,"%1.3f %1.3f L\n",ViewRight,ViewTop);
    fprintf(FP,"C F GR\n");
}

void PsSetWorld( float xmin , float ymin , float xmax , float ymax )

{
    WinXmin = xmin;
    WinYmin = ymin;
    WinXmax = xmax;
    WinYmax = ymax;

    PsAdjustScale();
/*    PsSetClipWindow(); */
}

void PsRectifyScale( void )

{
	if( SxWtoV > SyWtoV )
		SxWtoV = SyWtoV;
	else
		SyWtoV = SxWtoV;

	if( SxWtoV != 0. ) {
	        SxVtoW = 1./SxWtoV;
       		SyVtoW = 1./SyWtoV;
	}

	PsSetViewport( PsxWtoV(WinXmin) , PsyWtoV(WinYmin) , 
			PsxWtoV(WinXmax) , PsyWtoV(WinYmax) );
}

/***********************************************************************/
/*********************** Plotting Routines *****************************/
/***********************************************************************/

void PsLine( float x1 , float y1 , float x2 , float y2 )

{
    float vx1,vy1,vx2,vy2;

    vx1=PsxWtoV(x1);
    vy1=PsyWtoV(y1);
    vx2=PsxWtoV(x2);
    vy2=PsyWtoV(y2);

    if(TEST) printf("PsLine: %f %f %f %f\n",vx1,vy1,vx2,vy2);

    if( !InSide(vx1,vy1) && !InSide(vx2,vy2) ) {
	BOutSide = TRUE;
	XAct=vx2;
	YAct=vy2;
	return;
    }

    if( ISWHITE ) return;
    SETCOLOR;
    if( vx1 != XAct || vy1 != YAct || NPath == 0 || BOutSide ) 
        fprintf(FP,"%1.3f %1.3f M\n",vx1,vy1);
    fprintf(FP,"%1.3f %1.3f L\n",vx2,vy2);

    XAct=vx2;
    YAct=vy2;
    NPlot++;
    BOutSide = FALSE;
    PsStroke(APPEND);
}

void PsPlot( float x , float y )

{
    float vx,vy;

    vx=PsxWtoV(x);
    vy=PsyWtoV(y);

    if(TEST) printf("PsPlot: %f %f %f %f\n",XAct,YAct,vx,vy);

    if( !InSide(vx,vy) && !InSide(XAct,YAct) ) {
	BOutSide = TRUE;
	XAct=vx;
	YAct=vy;
	return;
    }

    if( ISWHITE ) return;
    SETCOLOR;
    if( NPath == 0 || BOutSide )
        fprintf(FP,"%1.3f %1.3f M\n",XAct,YAct);
    fprintf(FP,"%1.3f %1.3f L\n",vx,vy);

    XAct=vx;
    YAct=vy;
    NPlot++; 
    BOutSide = FALSE;
    PsStroke(APPEND);
}

void PsMove( float x , float y )

{
    float vx,vy;

    vx=PsxWtoV(x);
    vy=PsyWtoV(y);

    if(TEST) printf("PsMove: %f %f\n",vx,vy);

    if( !InSide(vx,vy) ) {
	BOutSide = TRUE;
	XAct=vx;
	YAct=vy;
	return;
    }

    fprintf(FP,"%1.3f %1.3f M\n",vx,vy);

    XAct=vx;
    YAct=vy;
    BOutSide = FALSE;
}

void PsPoint( float x , float y )

{
    float dd=PointSize;
    float vx,vy;
    float dx,dy;

    vx=PsxWtoV(x);
    vy=PsyWtoV(y);

    if(TEST) printf("PsPoint: %f %f\n",vx,vy);

    if( !InSide(vx,vy) ) {
	BOutSide = TRUE;
	XAct=vx;
	YAct=vy;
	return;
    }

    if( ISWHITE ) return;
    SETCOLOR;
    dx = PsxVtoW(dd) - PsxVtoW(0.);
    dy = PsyVtoW(dd) - PsyVtoW(0.);
    PsRectFill(x,y,x+dx,y+dy);
    XAct=vx;
    YAct=vy;
    NPlot++;
    BOutSide = FALSE;
}

void PsAreaFill( int ndim , float *x , float *y )

{
    int i;
    int in=0;
    float vx,vy;

    if(ndim<3) return;

    for(i=0;i<ndim;i++) {
	vx = PsxWtoV(x[i]);
	vy = PsyWtoV(y[i]);
        if( InSide(vx,vy) ) in++;
    }
    if( !in ) return;

    if( ISWHITE ) return;
    PsStroke(FLUSH);
    SETCOLOR;
    vx = PsxWtoV(x[0]);
    vy = PsyWtoV(y[0]);
    fprintf(FP,"N %1.3f %1.3f M\n",vx,vy);
    for(i=1;i<ndim;i++) {
	vx = PsxWtoV(x[i]);
	vy = PsyWtoV(y[i]);
	fprintf(FP,"%1.3f %1.3f L\n",vx,vy);
    }
    fprintf(FP,"C F\n");
    NPlot++;
}

void PsRectFill( float x1 , float y1 , float x2 , float y2 )

{
	float ix1,iy1,ix2,iy2;
	float vw,vh;
	float vx1,vy1,vx2,vy2;

	ix1 =  MIN(x1,x2);
	ix2 =  MAX(x1,x2);
	iy1 =  MIN(y1,y2);
	iy2 =  MAX(y1,y2);

	vx1 = PsxWtoV(ix1);
	vy1 = PsyWtoV(iy1);
	vx2 = PsxWtoV(ix2);
	vy2 = PsyWtoV(iy2);
	if( vx1 > ViewRight ) return;
	if( vx2 < ViewLeft ) return;
	if( vy1 > ViewTop ) return;
	if( vy2 < ViewBottom ) return;

	if( ISWHITE ) return;
	PsStroke(FLUSH);
	SETCOLOR;

	vw = vx2 - vx1;
	vh = vy2 - vy1;
	fprintf(FP,"%1.3f %1.3f %1.3f %1.3f RF\n",vh,vw,vx1,vy1);
	NPlot++;
}

void PsArc( float x0, float y0, float r, float ang1, float ang2 )

{
        float vx,vy,vr;

        vx=PsxWtoV(x0);
        vy=PsyWtoV(y0);
	vr=SxWtoV*r;		/* problem if scale x != y */
	vr=SxVtoW*r;		/* problem if scale x != y */
	vr=SxWtoV*r;		/* problem if scale x != y */

	if( ISWHITE ) return;
	PsStroke(FLUSH);
	SETCOLOR;

	/*
	fprintf(FP,"%% %1.3f %1.3f %1.3f %1.3f %1.3f arc\n",vx,vy,r,ang1,ang2);
	fprintf(FP,"%% %1.3f %1.3f \n",SxWtoV,SxVtoW);
	*/
	fprintf(FP,"N %1.3f %1.3f %1.3f %1.3f %1.3f arc S\n",vx,vy,vr,ang1,ang2);
	NPlot++;
}

/***********************************************************************/
/*************************** Text Routines *****************************/
/***********************************************************************/

static void PsSetFont( void )

{
        PsStroke(FLUSH);
	fprintf(FP,"/%s FF ",Font);
/*	fprintf(FP,"%f CF ",TFACTOR*FontSize); */
	fprintf(FP,"%f CF ",FontSize);
	fprintf(FP,"SF\n");
	FontChanged = NO;
}

void PsTextFont( char *font )		/* sets new font */

{
    if( font && strncmp(Font,font,81) ) {
	strncpy(Font,font,81);
	FontChanged = YES;
    }
}
	
void PsTextRotate( float angle )

{
	FontRotate = angle;
}
	
void PsTextPointSize( int size )	/* size is in points */

{
	if( size <= 0 ) size = FONTSIZE;
	PsTextRealSize( size/28.3 );
}

void PsTextRealSize( float size )	/* size is in actual units (cm) */

{
    if( FontSize != size ) {
	FontSize = size;
	FontChanged = YES;
    }
}

void PsText( float x , float y , char *s )	/* writes string */

{
    float vx,vy;
    int append = FALSE;
    int appendx = FALSE, appendy = FALSE;

    if( x == 999.0 || y == 999.0 ) {
      fprintf(FP,"currentpoint\n");
      fprintf(FP,"/YY exch def\n");
      fprintf(FP,"/XX exch def\n");
      append = TRUE;
      appendx = ( x == 999.0 ? TRUE : FALSE );
      appendy = ( y == 999.0 ? TRUE : FALSE );
    }

    vx=PsxWtoV(x);
    vy=PsyWtoV(y);

    if( ISWHITE ) return;
    PsStroke(FLUSH);
    SETCOLOR; SETTEXT;
    if( append ) {
	if( appendx ) {
          fprintf(FP,"XX ");
	} else {
          fprintf(FP,"%f ",vx);
	}
	if( appendy ) {
          fprintf(FP,"YY ");
	} else {
          fprintf(FP,"%f ",vy);
	}
        fprintf(FP,"M\n");
    } else {
      fprintf(FP,"%f %f M\n",vx,vy);
    }
    if( MUSTROTATE ) fprintf(FP,"%1.3f rotate\n",FontRotate);
    fprintf(FP,"(%s) show\n",s);
    if( MUSTROTATE ) fprintf(FP,"%1.3f rotate\n",-FontRotate);
    NPlot++;
}

void PsTextDimensions( char *s , float *width , float *height )

{
	static float empiric = 0.7;	/* maybe 0.75 */

	*height = SyVtoW * FontSize / TFACTOR;
	*width  = (SxVtoW/SyVtoW) * empiric * (*height) * strlen(s);
/*
    *width  = SxVtoW * XTextWidth(MyFontStruct,s,strlen(s));
    *height = - ( SyVtoW * MyFontAscent );
*/
}

/***********************************************************************/
/************************** Utility Routines ***************************/
/***********************************************************************/

void PsComment( char *s )

{
    fprintf(FP,"%%PsComment: %s\n",s);
}

void PsRealXY( float vx , float vy , float *x , float *y )

{
	*x = PsxVtoW( vx );
	*y = PsyVtoW( vy );
}

void PsPageXY( float x , float y , float *vx , float *vy )

{
	*vx = PsxWtoV( x );
	*vy = PsyWtoV( y );
}

/********************************************************************/

void PsFlush( void )

{
	PsStroke(FLUSH);
}

/********************************************************************/

void PsSync( void ) { }
void PsSyncron( int sync ) { }
int PsActScreen( void ) { return 0; }
void PsSetTest( void ) { TEST = 1; }
void PsClearTest( void ) { TEST = 0; }

void PsGetMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax ) { }
void PsSetMinMaxPix( int xmin , int ymin , int xmax , int ymax ) { }

