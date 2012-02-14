
/************************************************************************\ 
 *									*
 * meshge.c - geometric routines for mesh 				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
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
 * 05-Dec-2011: new routines for inLine check				*
 * 05-Dec-2011: new routine AreaLine					*
 * 05-Dec-2011: in AreaConvex factor 0.5 was missing (and stable 64 bit)*
 * 09-Oct-2010: better error handling for CircumCircle Error		*
 * 10-Mar-2010: new way to compute area of polygon (stable for 64 bit)	*
 * 03-Jul-2000: in ControlCircumCircle() only warning			*
 * 17-Oct-97: new routine TurnClosedLine()                              *
 * 16-Oct-97: new routine FindElement()                                 *
 * 15-Oct-97: new routine InClosedLine()                                *
 * 08-Oct-97: error message for identical nodes in MakeCircumCircle()   *
 * 01-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <math.h>

#include "general.h"
#include "gustd.h"
#include "nlist.h"

#include "mesh.h"
#include "meshge.h"
#include "meshfi.h"
#include "meshhs.h"


void MakeCircumCircle( Elem_type *pe )

/*\
 *  Find radius and center of circumcircle
 *
 *  Center at intersection of "Mittelsenkrechten (MSR)".
 *
 *  use formula : (nx,ny)*(x-x0,y-y0) = 0
 *  where (nx,ny) is normal to MSR (side of triangle) and
 *  (x0,y0) is intersection of side with MSR.
 *  -> a1*x + b1*y + c1 = 0 , a2*x + b2*y + c2 = 0
 *  det = a1*b2-a2*b1
 *  x0 = (b1*c2-b2*c1)/det , y0 = (c1*a2-c2*a1)/det
\*/

{
	int i;
/*
	float x[3],y[3];
	float a1,b1,c1,a2,b2,c2;
	float det;
	float x0,y0,rho;
*/
	double x[3],y[3];
	double a1,b1,c1,a2,b2,c2;
	double det;
	double x0,y0,rho;
	Node_type *pn;

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
	}

	a1 = x[1] - x[0];
	b1 = y[1] - y[0];
	c1 = - 0.5 * ( a1*(x[1]+x[0]) + b1*(y[1]+y[0]) );

	a2 = x[2] - x[1];
	b2 = y[2] - y[1];
	c2 = - 0.5 * ( a2*(x[2]+x[1]) + b2*(y[2]+y[1]) );

	if( pe->number == -1 ) {
		printf("%f %f %f\n",a1,b1,c1);
		printf("%f %f %f\n",a2,b2,c2);
		printf("%f %f %f\n",x[0],x[1],x[2]);
		printf("%f %f %f\n",y[0],y[1],y[2]);
		printf("%f\n",a1*b2 - a2*b1);
	}
		
	det = a1*b2 - a2*b1;
	if( det == 0 ) {
	  printf("*** error in MakeCircumCircle\n");
	  for(i=0;i<3;i++) {
	    printf("node: %d, x/y: (%f,%f)\n",pe->index[i],x[i],y[i]);
	  }
	  printf("Maybe the node coordinates are not unique\n");
	  printf("or the relative distance between the nodes\n");
	  printf("is too small with respect to the absolute\n");
	  printf("values of the coordinates or the three nodes\n");
	  printf("are collinear.\n");
	  WriteAll("M_error.grd",NULL);
	  Error("MakeCircumCircle: Error making element");
	}
	det = 1. / det;
	
	x0 = (b1*c2-b2*c1)*det;
	y0 = (c1*a2-c2*a1)*det;

	rho = (x0-x[0])*(x0-x[0])+(y0-y[0])*(y0-y[0]);

	pe->rc.x = x0;
	pe->rc.y = y0;
	pe->rho = rho;
}

void ControlCircumCircle( Elem_type *pe )

{
	int i;
	int stop=0;
	float x0,y0,rho;
	float x[3],y[3];
	float rhoaux;
	Node_type *pn;

	x0 = pe->rc.x;
	y0 = pe->rc.y;
	rho = pe->rho;

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
		rhoaux = (x0-x[i])*(x0-x[i])+(y0-y[i])*(y0-y[i]);
		if( ABS(rho-rhoaux) > 0.01*rho ) {
			printf("%d %f %f %f %f\n",pe->number,rho,rhoaux
					,x0,y0);
			stop=1;
		}
	}

	if( stop ) {
		/*Error("ControlCircumCircle: Error in circumcircle");*/
		Warning("ControlCircumCircle: Error in circumcircle");
	}
}

float MakeInCircleRadius( Elem_type *pe )

{
	int i;
	float x[3],y[3],a[3];
	float dx,dy,dd;
	float s=0.;
	Node_type *pn;

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
	}
	for(i=0;i<3;i++) {
		dx = x[(i+1)%3] - x[i];
		dy = y[(i+1)%3] - y[i];
		dd = sqrt( dx*dx+dy*dy );
		s += dd;
		a[i] = dd;
	}
	s *= 0.5;
	return (s-a[0])*(s-a[1])*(s-a[2]) / s ;
}

void MakeCM( Elem_type *pe, float *xm, float *ym )

{
	int i;
	Node_type *pn;

	*xm = 0.;
	*ym = 0.;

	for(i=0;i<pe->vertex;i++) {
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		*xm += pn->coord.x;
		*ym += pn->coord.y;
	}

	*xm /= pe->vertex;
	*ym /= pe->vertex;
}

int InCircumCircle( Elem_type *pe , float x , float y )

{
	float x0,y0,rho;

	x0 = pe->rc.x;
	y0 = pe->rc.y;
	rho = pe->rho;

	if( (x-x0)*(x-x0)+(y-y0)*(y-y0) < rho ) {
		return 1;
	}
	return 0;
}

int InElemCircumCircle( Elem_type *circle , Elem_type *pe , float fact )

{
	int i;
	float x0,y0,rho;
	float x,y;
	Node_type *pn;

	x0 = circle->rc.x;
	y0 = circle->rc.y;
	rho = fact*circle->rho;

	for(i=0;i<3;i++) {
	    pn = RetrieveByNodeNumber(HNN,pe->index[i]);
	    x = pn->coord.x;
	    y = pn->coord.y;
	    if( (x-x0)*(x-x0)+(y-y0)*(y-y0) < rho ) {
		return 1;
	    }
	}
	return 0;
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

int InElement( Hashtable_type H , Elem_type *pe , float x , float y )

/* error in routine adjusted 28.7.95 -> loop one more time */

{
        Rect r;
        Node_type *pn;
        int i,nvert;
        float scal,scao;
        Point *cn,*co;

        nvert = pe->vertex;

        if( nvert < 3 ) return 0;

        PolyMinMax( H , pe , &r );
        if( r.low.x>x || r.high.x<x || r.low.y>y || r.high.y<y ) return 0;

        scao=0.;
        pn = RetrieveByNodeNumber(H,pe->index[nvert-1]);
        co = &pn->coord;
	for(i=0;i<=nvert;i++) {
		pn = RetrieveByNodeNumber(H,pe->index[i%nvert]);
                cn = &pn->coord;
                scal = (x - co->x)*(cn->y - co->y)-(y - co->y)*(cn->x - co->x);
                if( scao*scal < 0. ) return 0;
                scao = scal;
                co = cn;
        }
        return 1;
}

void PolyMinMax( Hashtable_type H , Elem_type *pe , Rect *r )

{
        PolyMinMaxIndex( H , pe->vertex , pe->index , r );
}

void PolyMinMaxIndex( Hashtable_type H , int nvert , int *index , Rect *r )

{
        Node_type *p;
        Point v;

        if( nvert-- == 0 ) return;

        p=RetrieveByNodeNumber(H,*index);
        if( !p ) Error2("PolyMinMax : Cannot retrieve node ",
                                itos(*index));
        r->high = p->coord;
        r->low  = p->coord;
        index++;

        while( nvert-- > 0 ) {
                p=RetrieveByNodeNumber(H,*index);
                if( !p ) Error2("PolyMinMax : Cannot retrieve node ",
                                        itos(*index));
                v=p->coord;
                if( v.x > r->high.x ) r->high.x = v.x;
                if( v.x < r->low.x  ) r->low.x  = v.x;
                if( v.y > r->high.y ) r->high.y = v.y;
                if( v.y < r->low.y  ) r->low.y  = v.y;
                index++;
        }
}

Elem_type *FindElement( Hashtable_type H, float x, float y )

{
	Elem_type *pe;

	ResetHashTable(H);
	while( ( pe = VisitHashTableE(H) ) != NULL ) {
		if( InElement( HNN , pe , x , y ) ) return pe;
	}
	return NULL;
}


Elem_type *FindXYElement( ListTable L, float x, float y )

{
	Elem_type *pe;

	ResetListTable(L);
	while( ( pe = (Elem_type *) VisitListTable(L) ) != NULL ) {
		if( InElement( HNN , pe , x , y ) ) return pe;
	}
	ResetListTable(L);
	while( ( pe = (Elem_type *) VisitListTable(L) ) != NULL ) {
		if( InCircumCircle( pe , x , y ) ) return pe;
	}
	return NULL;
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
                /* area += co->x * cn->y - cn->x * co->y; */
                area += (co->x + cn->x) * (cn->y - co->y);
		co = cn;
        }

        return 0.5 * (float) area;
}


float AreaLine( Hashtable_type H , Line_type *pl )

/* should work for closed and not closed lines */

{
        int i,nvert;
        double area=0.;
        Point *co,*cn;
        Node_type *pn;

	nvert=pl->vertex;

        pn = RetrieveByNodeNumber(H,pl->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(H,pl->index[i]);
                cn = &pn->coord;
                /* area += co->x * cn->y - cn->x * co->y; */
                area += (co->x + cn->x) * (cn->y - co->y);
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

	nvert=list->count;

        pn = RetrieveByNodeNumber(HNN,list->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HNN,list->index[i]);
                cn = &pn->coord;
                /* area += co->x * cn->y - cn->x * co->y; */
                area += (co->x + cn->x) * (cn->y - co->y);
		co = cn;
        }

	return 0.5 * (float) area;
}
	
float LengthConvex( NodeList list )

{
	int nvert;
	int i;
	double length=0.;
	double dx,dy;
	Node_type *pn;
	Point *co, *cn;

	nvert=list->count;

        pn = RetrieveByNodeNumber(HNN,list->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HNN,list->index[i]);
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
	    /*
	    if( a < 0. ) {
		printf("*** angle negative: %d %f\n",i,a);
	    }
	    */
	    if( a < 0. ) a += pi2;	/* ensure the angle is positive */
	    alpha += a - pi;
	    /* printf("angle: %d  %f  %f  %f\n",i,a,a-pi,alpha); */
	    xold = xl[i];
	    yold = yl[i];
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

int IsPointInLine( Line_type *pl , float x , float y )

{
	int i,loops;
	int nvert;
	double alpha = 0.;
	double pi2 = 2.*3.14159;
        Point *co,*cn;
        Node_type *pn;

	nvert=pl->vertex;

        pn = RetrieveByNodeNumber(HNN,pl->index[nvert-1]);
	co = &pn->coord;
        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HNN,pl->index[i]);
                cn = &pn->coord;
	        alpha += angle(co->x,co->y,x,y,cn->x,cn->y);
		co = cn;
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
       	
int IsLineInLine( Line_type *plext , Line_type *plint )

{
	int i,nvert,inside;
        Point *cn;
        Node_type *pn;

	nvert=plint->vertex;

        for(i=0;i<nvert;i++) {
                pn = RetrieveByNodeNumber(HNN,plint->index[i]);
                cn = &pn->coord;
		inside = IsPointInLine(plext,cn->x,cn->y);
		if( ! inside ) return 0;
        }

	return 1;
}
       	
int IsLineClosed( Line_type *pl )

{
	if( pl->index[0] == pl->index[pl->vertex-1] ) {
		return 1;
	} else {
		return 0;
	}
}



