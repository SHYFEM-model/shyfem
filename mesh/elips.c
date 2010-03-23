
/************************************************************************\
 *                                                                      *
 * elips.c - creates ellips for use with mesh generator                 *
 *                                                                      *
 * Copyright (c) 1995 by Georg Umgiesser                                *
 *                                                                      *
 * Permission to use, copy, modify, and distribute this software        *
 * and its documentation for any purpose and without fee is hereby      *
 * granted, provided that the above copyright notice appear in all      *
 * copies and that both that copyright notice and this permission       *
 * notice appear in supporting documentation.                           *
 *                                                                      *
 * This file is provided AS IS with no warranties of any kind.          *
 * The author shall have no liability with respect to the               *
 * infringement of copyrights, trade secrets or any patents by          *
 * this file or any part thereof.  In no event will the author          *
 * be liable for any lost revenue or profits or other special,          *
 * indirect and consequential damages.                                  *
 *                                                                      *
 * Comments and additions should be sent to the author:                 *
 *                                                                      *
 *                      Georg Umgiesser                                 *
 *                      ISDGM/CNR                                       *
 *                      S. Polo 1364                                    *
 *                      30125 Venezia                                   *
 *                      Italy                                           *
 *                                                                      *
 *                      Tel.   : ++39-41-5216875                        *
 *                      Fax    : ++39-41-2602340                        *
 *                      E-Mail : georg@lagoon.isdgm.ve.cnr.it           *
 *                                                                      *
 * Revision History:                                                    *
 * ..-Aug-95: routine written from scratch                              *
 *                                                                      *
\************************************************************************/


#include <stdio.h>
#include <math.h>

/* creates ellipse */

void main( void )

{
	int npoints = 30;
	float rx = 43.;
	float ry = 35.;
	float pi = 3.14159;
	int grd = 1;

	int i;
	float t;
	float ddeg;
	float x,y;

	ddeg = 360./npoints;

	if( grd == 0 ) {
		printf("1 %d\n",npoints);
	}

	for(i=0;i<npoints;i++) {
		t = i * ddeg * pi / 180.;
		x = rx * cos(t);
		y = ry * sin(t);
		if( grd == 0 ) {
		    printf("%f %f\n",x,y);
		} else {
		    printf("1 %d 2 %f %f\n",i+1,x,y);
		}
	}

	if( grd == 1 ) {
		printf("\n");
		printf("3 1 1 %d",npoints+1);
		for(i=0;i<=npoints;i++) {
		    if( i%10 == 0 ) {
			printf("\n");
		    }
		    printf(" %d",(i%npoints)+1);
		}
		printf("\n");
	}
}

	
