
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
 * demopost.c - demonstration for postscript routines			*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "pgraph.h"

	float xf[] = {6.,6.,7.};
	float yf[] = {6.,7.,7.};

int main()

{
	float xm,dm,x,y;
	int i;

	QGraphInit();

	QMove(0.,0.);
	QPlot(10.,0.);
	QPlot(10.,10.);
	QPlot(0.,10.);
	QPlot(0.,0.);

	QNewPen(5);

	QMove(1.,1.);
	QPlot(9.,1.);
	QPlot(9.,9.);
	QPlot(1.,9.);
	QPlot(1.,1.);

	QNewPen(10);

	QLine(-2.,5.,12.,5.);
	QLine(5.,-2.,5.,12.);

	QLine(2.,2.,8.,8.);
	QLine(2.,8.,8.,2.);

	QRectFill(3.,3.,5.,5.,13);
	QAreaFill(3,xf,yf);
	QPoint(8.,8.5);
/*
	QGraphClose();

	exit(0);
*/
	QClearScreen();

	QNewPen(15);
	QText(4.,9.5,"iiiiiiiiiiii");
	QText(4.,9.,"aaaaaaaaaaaa");
	QText(4.,8.,"Available Colors");

	xm=10./8./2.;
	dm=0.7*xm;

	x=xm;
	y=6.;
	for(i=0;i<8;i++) {
	  QRectFill(x-dm,y-dm,x+dm,y+dm,i);
	  x+=2.*xm;
	}
	x=xm;
	y=3.;
	for(i=8;i<16;i++) {
	  QRectFill(x-dm,y-dm,x+dm,y+dm,i);
	  x+=2.*xm;
	}

	QFlush();
	printf("Enter <CR> to finish\n");
	getchar();

	QGraphClose();

	exit(0);

	QClearScreen();

/*
	QNewPen(0);
	QText(4.,8.,"Blue Colors");

	xm=10./8./2.;
	dm=0.7*xm;

	x=xm;
	y=6.;
	for(i=0;i<8;i++) {
	  QRectFill(x-dm,y-dm,x+dm,y+dm,BlueShadeColor(i));
	  x+=2.*xm;
	}
	x=xm;
	y=3.;
	for(i=8;i<16;i++) {
	  QRectFill(x-dm,y-dm,x+dm,y+dm,BlueShadeColor(i));
	  x+=2.*xm;
	}

	QFlush();
	printf("Enter <CR> to finish\n");
	getchar();
*/

	QGraphClose();

	return 0;
}
