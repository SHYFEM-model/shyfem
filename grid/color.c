
/************************************************************************\
 *
 *    Copyright (C) 1995,1997  Georg Umgiesser
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
 * color.c - general color routines
 *
 * revision log :
 *
 * 11.02.1995	ggu	own file from xgraph.c
 * 16.02.1995	ggu	bug in linear, quadratic... -> returned -1 sometimes
 * 10.10.1997	ggu	new routines Hue2RBG(), QAllocHSBColors()
 *
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

