
/************************************************************************\ 
 *									*
 * color.c - general color routines                                     *
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
 * 10-Oct-97: new routines Hue2RBG(), QAllocHSBColors()                 *
 * 16-Feb-95: bug in linear, quadratic... -> returned -1 sometimes      *
 * 11-Feb-95: own file from xgraph.c                                    *
 *									*
\************************************************************************/


#include <stdlib.h>

#include "general.h"
#include "graph.h"


/*****************************************************************/


int linear(int start , int end , int steps , int act )

{
	int aux;

	aux = start + ( ( end - start ) * act ) / steps ;

	if( end-start > 0 ) {
	    if( aux < start ) aux = start;
	    if( aux > end ) aux = end;
	} else {
	    if( aux > start ) aux = start;
	    if( aux < end ) aux = end;
	}
	return aux;
}

int quadratic_fast
			(
			  int start 
			, int end 
			, int steps 
			, int deriv
			, int act 
			)

{
	int dy,derivative,aux;

	dy = end - start;
	derivative = ( dy * deriv > 0 ) ? deriv : -deriv;
	aux = steps * derivative;
	aux = start + ((2*dy-aux)*act)/steps 
			+ ((aux-dy)*act*act)/(steps*steps);
	if( dy > 0 ) {
	    if( aux < start ) aux = start;
	    if( aux > end ) aux = end;
	} else {
	    if( aux > start ) aux = start;
	    if( aux < end ) aux = end;
	}
	return aux;
}

int quadratic_slow
			(
			  int start 
			, int end 
			, int steps 
			, int deriv
			, int act 
			)

{
	int dy,derivative,aux;

	dy = end - start;
	derivative = ( dy * deriv > 0 ) ? deriv : -deriv;
	aux = derivative;
	aux = start + aux*act + ((dy-aux*steps)*act*act)/(steps*steps);
	if( dy > 0 ) {
	    if( aux < start ) aux = start;
	    if( aux > end ) aux = end;
	} else {
	    if( aux > start ) aux = start;
	    if( aux < end ) aux = end;
	}
	return aux;
}

/*****************************************************************/

int *QAllocGreen2BlueColors( int ncol )

/* allocate ncol color tones from green to blue */

{
	int i;
	int red,green,blue;
	int colmin=0;
	int colmax=255;
	int auxmax,derv;
	int *colorp;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Green2Blue color");

	auxmax=(3*colmax)/4;
	derv=2;

	for(i=0;i<ncol;i++) {

		red=colmin;
		green=quadratic_slow(auxmax,colmin,ncol,derv,i);
		blue=quadratic_fast(colmin,auxmax,ncol,derv,i);

		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

int *QAllocYellow2GreenColors( int ncol )

/* allocate ncol color tones from yellow to green */

{
	int i;
	int red,green,blue;
	int colmin=0;
	int colmax=255;
	int auxmax;
	int *colorp;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Yellow2Green color");

	auxmax=(3*colmax)/4;

	for(i=0;i<ncol;i++) {

		red=linear(colmax,colmin,ncol,i);
		green=linear(colmax,auxmax,ncol,i);
		blue=colmin;

		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

int *QAllocRed2YellowColors( int ncol )

/* allocate ncol color tones from red to yellow */

{
	int i;
	int red,green,blue;
	int colmin=0;
	int colmax=255;
	int *colorp;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Red2Yellow color");

	for(i=0;i<ncol;i++) {

		red=colmax;
		green=linear(colmin,colmax,ncol,i);
		blue=colmin;

		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

int *QAllocRed2BlueColors( int ncol )

/* allocate ncol color tones from red to blue */

{
	int i;
	int red,green,blue;
	int colmin=0;
	int colmax=255;
	int auxmax,derv;
	int *colorp;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Green2Blue color");

	auxmax=(3*colmax)/4;
	derv=2;

	for(i=0;i<ncol;i++) {

		red=quadratic_slow(colmax,colmin,ncol,derv,i);
		green=colmin;
		blue=quadratic_fast(colmin,auxmax,ncol,derv,i);

		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

int *QAllocBlueColors( int ncol )

/* allocate ncol blue color tones */

{
	int i;
	int red,green,blue;
	int colmin=0;
	int colmax=255;
	int auxmax;
	int *colorp;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Blue color");

	auxmax=(3*colmax)/4;

	for(i=0;i<ncol;i++) {

		red=linear(auxmax,colmin,ncol,i);
		green=linear(auxmax,colmin,ncol,i);
		blue=linear(colmax,(2*colmax)/3,ncol,i);;

		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

static void Hue2RBG( float hue, int colmax, int *red, int *green, int *blue )

{
	static float third = 1./3.;

	if( hue < third ) {
		*red = colmax * ( 1. - hue/third );
		*green = colmax *  hue/third;
		*blue = 0;
	} else if( hue > 2.*third ) {
		hue = hue - 2.*third;
		*red = colmax *  hue/third;
		*green = 0;
		*blue = colmax * ( 1. - hue/third );
	} else {
		hue = hue - third;
		*red = 0;
		*green = colmax * ( 1. - hue/third );
		*blue = colmax *  hue/third;
	}
}

int *QAllocHSBColors( int ncol )

/* allocate ncol HSB color tones */

{
	int i;
	int red,green,blue;
	int colmax=255;
	int *colorp;
	float hue;

	colorp = (int *) malloc( ncol * sizeof(int) );
	if( !colorp )
		Error("Cannot allocate Blue color");

	for(i=0;i<ncol;i++) {

		hue = ((float) i)/ncol;
		Hue2RBG(hue,colmax,&red,&green,&blue);
		/* printf("color... %d %f %d %d %d\n",i,hue,red,green,blue); */
		colorp[i] = QAllocColor(red,green,blue);
	}

	return colorp;
}

