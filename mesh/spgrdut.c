
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
 * meshck.c - check routines for mesh 					*
 *									*
 * Revision History:							*
 * 17-Oct-97: in CheckInput() check for couter-clockwise line turning   *
 * 08-Oct-97: uses new mesh type                                        *
 *            check for closed line only for L_EXTERNAL/L_INTERNAL      *
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general.h"
#include "gustd.h"
#include "stack.h"

#include "nlist.h"

#include "grd.h"
#include "grdhs.h"
#include "grdut.h"

#include "spgrdut.h"



void MakeCounterClockwise( void )

{
	int turn;
	int bstop = FALSE;
        Line_type *pl;
	float *xe, *ye;
	Hashtable_type HL = GetGrid()->HL;

        ResetHashTable(HL);
        while( (pl=VisitHashTableL(HL)) != NULL ) {

		MakeCoordsFromLine(pl,&xe,&ye);

		/* check for counter-clockwise turning */

		turn = TurnClosedLine(pl->vertex,xe,ye);
		if( turn == 1 ) {
			/* ok */;
		} else if( turn == -1 ) {
			printf("Line %d in clockwise sense -> inverting\n"
				,pl->number);
			InvertIndex(pl->index,pl->vertex);
		} else {
			printf("Line %d , turning number %d\n",
				pl->number,turn);
			bstop = TRUE;
		}
		free(xe); free(ye);

        }

	if( bstop )
		Error("Error in line turning");
}

void MakeCM( Elem_type *pe, float *xm, float *ym )

{
	int i;
	Node_type *pn;
	Hashtable_type HN = GetGrid()->HN;

	*xm = 0.;
	*ym = 0.;

	for(i=0;i<pe->vertex;i++) {
		pn = RetrieveByNodeNumber(HN,pe->index[i]);
		*xm += pn->coord.x;
		*ym += pn->coord.y;
	}

	*xm /= pe->vertex;
	*ym /= pe->vertex;
}

int InConvex( int n , float *xe , float *ye , float x , float y )

{
        int i,ii;
        float xmin,xmax,ymin,ymax;
        double scal,sc;

        if( n < 3 ) return 0;

        xmin=xe[0]; xmax=xmin; ymin=ye[0]; ymax=ymin;
        for( i=1 ; i<n ; i++ ) {
          if(xe[i]<xmin) xmin=xe[i];
          if(xe[i]>xmax) xmax=xe[i];
          if(ye[i]<ymin) ymin=ye[i];
          if(ye[i]>ymax) ymax=ye[i];
        }

        if( xmin>x || xmax<x || ymin>y || ymax<y ) return 0;

        ii=n-1;
        scal=(x-xe[ii])*(ye[0]-ye[ii])-(y-ye[ii])*(xe[0]-xe[ii]);
        for( i=1 ; i<n ; i++ ) {
                ii=i-1;
                sc=(x-xe[ii])*(ye[i]-ye[ii])-(y-ye[ii])*(xe[i]-xe[ii]);
                if( sc*scal < 0. ) return 0;
        }

        return 1;
}

float AreaElement( Hashtable_type H , Elem_type *pe )

{
        int i,nvert;
        double area=0.;
        Point *co,*cn;
        Node_type *pn;

	nvert=pe->vertex;

        pn = RetrieveByNodeNumber(H,pe->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(H,pe->index[i]);
                cn = &pn->coord;
                area += co->x * cn->y - cn->x * co->y;
		co = cn;
        }

        return 0.5 * (float) area;
}


float AreaConvex( NodeList list )

{
	int nvert;
	int i;
	double area=0.;
	Node_type *pn;
	Point *co, *cn;
	Hashtable_type HN = GetGrid()->HN;

	nvert=list->count;

        pn = RetrieveByNodeNumber(HN,list->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HN,list->index[i]);
                cn = &pn->coord;
                area += co->x * cn->y - cn->x * co->y;
		co = cn;
        }

	return (float) area;
}
	
float LengthConvex( NodeList list )

{
	int nvert;
	int i;
	double length=0.;
	double dx,dy;
	Node_type *pn;
	Point *co, *cn;
	Hashtable_type HN = GetGrid()->HN;

	nvert=list->count;

        pn = RetrieveByNodeNumber(HN,list->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HN,list->index[i]);
                cn = &pn->coord;
		dx = co->x - cn->x;
		dy = co->y - cn->y;
                length += sqrt(dx*dx+dy*dy);
		co = cn;
        }

	return (float) length;
}

double angle( float x1, float y1, float x2, float y2, float x3, float y3 )

{
	double ax,ay,bx,by;
	double da,db,mod,alpha;

	ax = x1 - x2;
	ay = y1 - y2;
	bx = x3 - x2;
	by = y3 - y2;

	da = ax*ax + ay*ay;
	db = bx*bx + by*by;

	mod = sqrt( da * db );

	alpha = acos( (ax*bx + ay*by) / mod );

	if( ax*by - ay*bx < 0. ) alpha = -alpha;

	return alpha;
}

int InClosedLine( int ndim , float *xl , float *yl , float x , float y )

{
	int i,loops;
	double alpha = 0.;
	double pi2 = 2.*3.14159;

	for(i=1;i<ndim;i++) {
	    alpha += angle(xl[i-1],yl[i-1],x,y,xl[i],yl[i]);
	}

	if( alpha > 0. ) {
		loops = ( alpha + 0.5 ) / pi2;
	} else if( alpha < 0. ) {
		loops = ( alpha - 0.5 ) / pi2;
	} else {
		loops = 0;
	}

	return loops;
}
       	
int TurnClosedLine( int ndim , float *xl , float *yl )

{
	int i,loops;
	float xold,yold;
	double alpha = 0.;
	double a;
	double pi = 3.14159;
	double pi2 = 2.*pi;

	ndim--;		/* line is closed -> do not use first node twice */
	xold = xl[ndim-1];
	yold = yl[ndim-1];

	for(i=0;i<ndim;i++) {
	    a = angle(xold,yold,xl[i],yl[i],xl[i+1],yl[i+1]);
	    if( a < 0. ) a += pi2;	/* ensure the angle is positive */
	    alpha += a;
	    xold = xl[i];
	    yold = yl[i];
	}

	alpha -= pi * ndim;

	if( alpha > 0. ) {
		loops = ( alpha + 0.5 ) / pi2;
	} else if( alpha < 0. ) {
		loops = ( alpha - 0.5 ) / pi2;
	} else {
		loops = 0;
	}

	return loops;
}

void MakeCoordsFromLine( Line_type *pl , float **x , float **y )

{
	int i,node;
	int ncount = pl->vertex;
	float *xp, *yp;
	Node_type *pn;
	Hashtable_type HN = GetGrid()->HN;

	xp = (float *) malloc( ncount * sizeof(float) );
	yp = (float *) malloc( ncount * sizeof(float) );

	if( !xp || !yp )
	    Error("MakeCoordsFromLine: Cannot allocate coordinate list");

	for(i=0;i<ncount;i++) {
	    node = pl->index[i];
	    pn = RetrieveByNodeNumber(HN,node);
	    xp[i] = pn->coord.x;
	    yp[i] = pn->coord.y;
	}

	*x = xp;
	*y = yp;
}

