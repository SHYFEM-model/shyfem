
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
 *                                                                      *
 * elips.c - creates ellips for use with mesh generator                 *
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

	
