
/************************************************************************\ 
 *									*
 * gridge.c - geometric manipulation routines				*
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
 * 01-Oct-2004: GetDepthMinMax() to compute maximum depth               *
 * 06-Dec-95: IsDegenerateRect()                                        *
 * 06-Dec-95: new routine for computing min/max of rectangle            *
 * 21-Oct-94: new routine PolyMinMaxIndex for any node index            *
 * 13-Apr-94: use new hash routines                                     *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdio.h>
#include <math.h>

#include "grid.h"
#include "gustd.h"

#include "gridhs.h"

void GetElemMinMax( Hashtable_type HE , Hashtable_type HN , Rect *r )

{
        Elem_type *p;
	Rect triang;

	if( !HE ) Error("GetElemMinMax : No Hashtable given");

	ResetHashTable(HE);
	p = VisitHashTableE(HE);

        if( !p ) Error("GetElemMinMax : Hashtable is empty");

	PolyMinMax( HN , p , r );

	while( (p = VisitHashTableE(HE)) != NULL ) {
		PolyMinMax( HN , p , &triang );
		AdjustBounds( r , &triang );
        }
}

void InitMinMax( Rect *r )

{
	r->low.x  = 1.; r->low.y  = 1.;
	r->high.x = 0.; r->high.y = 0.;
}

void SetMinMax( Rect *r )

{
	float dx = r->high.x - r->low.x;
	float dy = r->high.y - r->low.y;

	if( dx < 0. && dy < 0. ) {		/* no points found */
		r->low.x  = 0.; r->low.y  = 0.;
		r->high.x = 1.; r->high.y = 1.;
	} else if( dx == 0. && dy == 0. ) {	/* only one point */
		r->low.x -= 1.;
		r->low.y -= 1.;
		r->high.x += 1.;
		r->high.y += 1.;
	} else if( dx == 0. ) {			/* vertical line */
		r->low.x -= dy/2.;
		r->high.x += dy/2.;
	} else if( dy == 0. ) {			/* horizontal line */
		r->low.y -= dx/2.;
		r->high.y += dx/2.;
	}					/* else ok */

	/* leave some space around plot */

	dx = r->high.x - r->low.x;
	dy = r->high.y - r->low.y;

	r->low.x -= 0.1*dx;
	r->high.x += 0.1*dx;
	r->low.y -= 0.1*dy;
	r->high.y += 0.1*dy;
}

float GetDepthMinMax( void )

{
        Node_type *pn;
        Elem_type *pe;
        Line_type *pl;
	float depth_max = 0.;
	float depth;
	Hashtable_type H;

	H = GetHashTableN();
	ResetHashTable(H);

	while( (pn = VisitHashTableN(H)) != NULL ) {
	  depth = pn->depth;
	  if( depth == NULLDEPTH ) continue;
	  if( depth > depth_max ) depth_max = depth;
        }

	H = GetHashTableE();
	ResetHashTable(H);

	while( (pe = VisitHashTableE(H)) != NULL ) {
	  depth = pe->depth;
	  if( depth == NULLDEPTH ) continue;
	  if( depth > depth_max ) depth_max = depth;
        }

	H = GetHashTableL();
	ResetHashTable(H);

	while( (pl = VisitHashTableL(H)) != NULL ) {
	  depth = pl->depth;
	  if( depth == NULLDEPTH ) continue;
	  if( depth > depth_max ) depth_max = depth;
        }

	return depth_max;
}

void GetNodeMinMax( Hashtable_type H , Rect *r )

{
        Node_type *p;
	Rect node;

        if( !H ) Error("GetNodeMinMax : No Hashtable given");

	ResetHashTable(H);
	p = VisitHashTableN(H);
	if( !p ) return;

	if( r->low.x > r->high.x && r->low.y > r->high.y ) {
		MakeRect( &(p->coord) , r );
	}

	while( (p = VisitHashTableN(H)) != NULL ) {
            MakeRect( &(p->coord) , &node );
	    AdjustBounds( r , &node );
        }
}

float GetAverLat( Hashtable_type H )

{
	Rect r;

        if( !H ) Error("GetAverLat : No Hashtable given");

	GetNodeMinMax( H , &r );

	return 0.5 * ( r.high.y - r.low.y );
}

int IsLatLon( Hashtable_type H )

{
	Rect r;

        if( !H ) Error("IsLatLon : No Hashtable given");

	GetNodeMinMax( H , &r );

	if( r.low.x > -100. && r.high.x < 100. ) {
	  if( r.low.y > -100. && r.high.y < 100. ) {
	    return 1;
	  }
	}
	return 0;
}

int IsDegenerateRect( Rect *r )

{
	if( r->low.x >= r->high.x || r->low.y >= r->high.y ) {
		return 1;
	} else {
		return 0;
	}
}

void CopyRect( Rect *d , Rect *s )

{
	d->low.x = s->low.x;
	d->low.y = s->low.y;
	d->high.x = s->high.x;
	d->high.y = s->high.y;
}

void MakeRect ( Point *p , Rect *r )

{
        r->low.x  = p->x;
        r->high.x = p->x;
        r->low.y  = p->y;
        r->high.y = p->y;
}

void MakeRectFromPoints ( Point *p1 , Point *p2 , Rect *r )

{
	if( p1->x < p2->x ) {
		r->low.x  = p1->x;
		r->high.x = p2->x;
	} else {
		r->low.x  = p2->x;
		r->high.x = p1->x;
	}

	if( p1->y < p2->y ) {
		r->low.y  = p1->y;
		r->high.y = p2->y;
	} else {
		r->low.y  = p2->y;
		r->high.y = p1->y;
	}
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

void TriMinMax( Hashtable_type H , int *index , Rect *r )

{
	Node_type *p;
	int i;
	Point v[3];

	for(i=0;i<3;i++) {
		p=RetrieveByNodeNumber(H,index[i]);
		if( !p ) Error2("TriMinMax : Cannot retrieve node ",
					itos(index[i]));
		v[i]=p->coord;
	}

        if( v[0].x < v[1].x ) {
            if( v[1].x < v[2].x ) {
                r->low.x  = v[0].x;
                r->high.x = v[2].x;
            } else if( v[0].x < v[2].x ) {
                r->low.x  = v[0].x;
                r->high.x = v[1].x;
            } else {
                r->low.x  = v[2].x;
                r->high.x = v[1].x;
            }
        } else {
            if( v[0].x < v[2].x ) {
                r->low.x  = v[1].x;
                r->high.x = v[2].x;
            } else if( v[1].x < v[2].x ) {
                r->low.x  = v[1].x;
                r->high.x = v[0].x;
            } else {
                r->low.x  = v[2].x;
                r->high.x = v[0].x;
            }
        }

        if( v[0].y < v[1].y ) {
            if( v[1].y < v[2].y ) {
                r->low.y  = v[0].y;
                r->high.y = v[2].y;
            } else if( v[0].y < v[2].y ) {
                r->low.y  = v[0].y;
                r->high.y = v[1].y;
            } else {
                r->low.y  = v[2].y;
                r->high.y = v[1].y;
            }
        } else {
            if( v[0].y < v[2].y ) {
                r->low.y  = v[1].y;
                r->high.y = v[2].y;
            } else if( v[1].y < v[2].y ) {
                r->low.y  = v[1].y;
                r->high.y = v[0].y;
            } else {
                r->low.y  = v[2].y;
                r->high.y = v[0].y;
            }
        }
}

void AdjustBounds( Rect *bounds , Rect *new )

{
        if( new->low.x  < bounds->low.x  ) bounds->low.x  = new->low.x;
        if( new->low.y  < bounds->low.y  ) bounds->low.y  = new->low.y;
        if( new->high.x > bounds->high.x ) bounds->high.x = new->high.x;
        if( new->high.y > bounds->high.y ) bounds->high.y = new->high.y;
}


float rangle( float x1, float y1, float x2, float y2, float x3, float y3 )

/*\
 *  computes angle between 3 points (vertex at point 2) in degrees
 *
 *  an angle less than 180. lets the curve turn in clockwise sense
 *  (the angle is meassured to the right of the directed curve)
\*/

{
	double dx1,dy1,dx2,dy2;
	double aux,ang;
	static double pi=3.14159;

	dx1=x1-x2;
	dy1=y1-y2;
	dx2=x3-x2;
	dy2=y3-y2;

	aux = ( dx1*dx2 + dy1*dy2 ) / 
		sqrt( (dx1*dx1+dy1*dy1) * (dx2*dx2+dy2*dy2) );

	ang = acos( aux ) * 180. / pi;

	if( dx1*dy2 - dy1*dx2 < 0. ) ang = 360. - ang;

/*	printf(" %f ",(float)ang); */
	return (float) ang;
}

float angle( Hashtable_type H , int k1 , int k2 , int k3 )

/*\
 *  computes angle between nodes k1,k2,k3
\*/

{
	Node_type *p1,*p2,*p3;

	p1 = RetrieveByNodeNumber(H,k1);
	p2 = RetrieveByNodeNumber(H,k2);
	p3 = RetrieveByNodeNumber(H,k3);

	return    rangle(p1->coord.x,p1->coord.y
			,p2->coord.x,p2->coord.y
			,p3->coord.x,p3->coord.y
			);
}


