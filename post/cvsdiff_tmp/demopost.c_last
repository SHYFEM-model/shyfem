
/************************************************************************\ 
 *									*
 * demopost.c - demonstration for postscript routines			*
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
