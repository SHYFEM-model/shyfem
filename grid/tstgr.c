
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


#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "mousex.h"

float x1[3] = {2.,4.,3.};
float y1[3] = {2.,2.,4.};
float x2[3] = {5.,7.,8.};
float y2[3] = {5.,5.,7.};

void main()

{
	int i;
	float xreal,yreal;
	float x,y,dx,dy;
	int flags;
	MouseEvent event;

	QGraphInit();
	
	QViewport(30,30,130,130);
	QWindow(0.,0.,10.,10.);

	QNewPen( White );
	QLine(2.,2.,8.,8.);
	QLine(1.,1.,9.,1.);
	QLine(9.,1.,9.,9.);
	QLine(9.,9.,1.,9.);
	QLine(1.,9.,1.,1.);
	QLine(-2.,5.,5.,-2.);
	QLine(12.,5.,5.,-2.);
	QLine(12.,5.,5.,12.);
	QLine(-2.,5.,5.,12.);

	QViewport(230,230,330,330);

	QNewPen(5);
	QLine(2.,2.,8.,8.);
	QLine(1.,1.,9.,1.);
	QLine(9.,1.,9.,9.);
	QLine(9.,9.,1.,9.);
	QLine(1.,9.,1.,1.);
	QLine(-2.,5.,5.,-2.);
	QLine(12.,5.,5.,-2.);
	QLine(12.,5.,5.,12.);
	QLine(-2.,5.,5.,12.);

	getchar();

	QRealXY(280,280,&xreal,&yreal);
	QClearViewport();

	getchar();

	QClearScreen();

	QNewPen(8);
	QMove(2.,8.);
	QPlot(8.,2.);
	QLine(1.,5.,9.,5.);

	QViewport(30,230,130,330);
	QPoint(1.,1.);
	QPoint(9.,1.);
	QPoint(1.,9.);
	QPoint(9.,9.);
	QNewPen(10);
	QAreaFill(3,x1,y1);
	QNewPen(14);
	QAreaFill(3,x2,y2);

	QViewport(50,50,550,150);
	QWindow(0.,0.,17.,2.);

	dx=0.4; dy=0.4;
	for(i=0;i<16;i++) {
		x=i+1; y=1.;
		QRectFill( x-dx,y-dy,x+dx,y+dy,i );
	}
		
	getchar();

	QClearScreen();

	QViewport(50,50,350,350);
	QWindow(0.,0.,10.,10.);
	QRectFill(0.,0.,10.,10.,3);
	QRectFill(1.,1.,5.,5.,7);
	QText(2.,2.,"georg");
	QTextSize(2);
	QText(2.,5.,"alessandra");

	getchar();

	MouseInit();
	MouseDisplayCursor();

	QNewPen(White);
	flags = M_EVENT;
	for(;;) {
	  MouseGetEvent(flags,&event);
	  if( event.flags & M_KEYPRESS ) break;
	  if( event.flags & M_BUTTON_CHANGE )
	     printf("%d %d %d\n",event.x,event.y,event.buttons);
	}
	MouseUnInit();

	QGraphClose();
}

/*


void QGraphInit( void );
void QClearScreen( void );
void QGraphClose( void );
void QNewPen( int color );
void QGetPen( int *color );
void QGetPixelWidth( float *dx , float *dy );
void QViewport( int left , int top , int right , int bottom );
void QClearViewport( void );
void QWindow( float xmin , float ymin , float xmax , float ymax );
void QAdjustScale( void );
void QGetMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax );
void QSetMinMaxPix( int xmin , int ymin , int xmax , int ymax );
void QLine( float x1 , float y1 , float x2 , float y2 );
void QPlot( float x , float y );
void QMove( float x , float y );
void QPoint( float x , float y );
void QAreaFill( int ndim , float *x , float *y );
void QRectFill( float x1 , float y1 , float x2 , float y2 , int color );
void QText( float x , float y , char *s );
void QTextSize( int size );
void QTextDimensions( char *s , float *width , float *height );
void QTextDimensionsI( char *s , int *width , int *height );
void QRealXY( int vx , int vy , float *x , float *y );
void QScreenXY( float x , float y , int *vx , int *vy );

void *QSavePixels( int x , int y , int width , int height );
void QRestorePixels( int x , int y , void *buffer );

void QFlush( void );
void QSync( void );
void QTest( void );
void QSyncron( int sync );

*/
