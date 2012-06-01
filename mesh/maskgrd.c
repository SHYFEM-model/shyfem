
/************************************************************************\ 
 *									*
 * maskgrd.c - manipulates grd file at will				*
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
 * 04-Nov-97: new routine to compute area of all elements               *
 * 28-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "general.h"
#include "gustd.h"

#include "hash.h"
#include "queue.h"
#include "stack.h"

#include "grd.h"
#include "grdio.h"
#include "grdhs.h"
#include "grdut.h"


void ReadFiles( int argc , char *argv[] , Grid_type *G );
void WriteFile( Grid_type *G );
void ManipulateGrd( Grid_type *G );
void MarkOttuso( Grid_type *G );
double angle( float x1, float y1, float x2, float y2, float x3, float y3 );
void SmoothDepth( Grid_type *G );
void Incline( Grid_type *G );
float getdepth( float y );
void ComputeArea( Grid_type *G );


void main(int argc, char *argv[])

{
	Grid_type *G;

	G = MakeGrid();

/*
 	SetOptions(argc,argv);
*/

	ReadFiles(argc,argv,G);

/*
	ManipulateGrd(G);
	WriteFile(G);
*/

	ComputeArea(G);

/*
	MarkOttuso(G);
	WriteFile(G);
*/

/*
	SmoothDepth(G);
	SmoothDepth(G);
	SmoothDepth(G);
	WriteFile(G);
*/
/*
	Incline(G);
	WriteFile(G);
*/

}


void ReadFiles( int argc , char *argv[] , Grid_type *G )

{
	char sfile[80];
	char *s;

	while( --argc ) {
		s=strcpy(sfile,*(++argv));
		s=strcat(s,".grd");
		ReadStandard(s,G);
	}
}

void WriteFile( Grid_type *G )

{
	char sfile[80] = "new.grd";
	char *s;

	s = sfile;
	WriteStandard(s,G);
}

void Incline( Grid_type *G )

{
	int i;
	float y;
	Elem_type *pe;
	Node_type *pn;
	Hashtable_type HNN = G->HN;
	Hashtable_type HEL = G->HE;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		y=0.;
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			y += pn->coord.y;
		}
		y /= pe->vertex;
		pe->depth = getdepth(y);
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		y=pn->coord.y;
		pn->depth = getdepth(y);
	}
}

float getdepth( float y )

{
	float d;
	static float d0 = 1000.;
	static float d1 = 10.;
	static float y0 = 0.;
	static float y1 = 400000.;

	if( y < y1 ) {
	  d = d0 + (d1-d0)*(y-y0)/(y1-y0);
	} else {
	  d = 10.;
	}

	return d;
}

void ManipulateGrd( Grid_type *G )

{
	int i;
	static float x0 = 0.;		/* defines line */
	static float y0 = 90000.;
	static float x1 = 122000.;
	static float y1 = 265000.;
	float x,y,s;
	float dx,dy;
	Elem_type *pe;
	Node_type *pn;
        Hashtable_type HNN = G->HN;
        Hashtable_type HEL = G->HE;

	dx = x1 - x0;
	dy = y1 - y0;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		x=0.; y=0.;
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			x += pn->coord.x;
			y += pn->coord.y;
		}
		x /= pe->vertex;
		y /= pe->vertex;
		s = dx * (y-y0) - dy * (x-x0);
		if( s > 0) {	/* to the left of the line */
			pe->type = 0;
		} else {
			pe->type = 1;
		}
	}
}

void MarkOttuso( Grid_type *G )

{
	int i;
	int i1,i2;
	int n90=0,n120=0;
	float x[3],y[3];
	float w,wmax;
	float eps = 1.e-3;
	Elem_type *pe;
	Node_type *pn;
        Hashtable_type HNN = G->HN;
        Hashtable_type HEL = G->HE;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->vertex != 3 ) continue;
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			x[i] = pn->coord.x;
			y[i] = pn->coord.y;
		}
		wmax=0.;
		for(i=0;i<pe->vertex;i++) {
			i1 = (i+1)%3;
			i2 = (i+2)%3;
			w = - angle(x[i],y[i],x[i1],y[i1],x[i2],y[i2]);
			wmax = ( w > wmax ) ? w : wmax;
		}
		if( wmax <= 90. + eps ) {
			pe->type = 0;
		} else if( wmax <= 120. + eps ) {
			pe->type = 1;
			n90++;
		} else {
			pe->type = 2;
			n120++;
		}
		pe->depth = wmax;
	}
	fprintf(stderr,"> 90 deg : %d - > 120 deg : %d\n",n90,n120);
}

double angle( float x1, float y1, float x2, float y2, float x3, float y3 )

/* angle between 1-2-3 , 2 is vertex, positiv if 1 is to the right of 3 */
/* for normal numeration of elements gives negative angle, so           */
/* exchange 1 with 3 or take negative                                   */

{
        double ax,ay,bx,by;
        double da,db,mod,alpha;
	static double pi=3.14159;

        ax = x1 - x2;
        ay = y1 - y2;
        bx = x3 - x2;
        by = y3 - y2;

        da = ax*ax + ay*ay;
        db = bx*bx + by*by;

        mod = sqrt( da * db );

        alpha = acos( (ax*bx + ay*by) / mod );

        if( ax*by - ay*bx < 0. ) alpha = -alpha;

        return alpha * 180./pi;
}

/****************** for depth smoothing **************************/

typedef struct {
	int count;
	float accum;
} Count_type;

Count_type *MakeCount( void )

{
	Count_type *item;

	item = (Count_type *) malloc( sizeof(Count_type) );
	if( !item )
		Error("Cannot allocate Count type");

	item->count = 0;
	item->accum = 0.;

	return item;
}

void SmoothDepth( Grid_type *G )

{
	int i;
	float d;
	Elem_type *pe;
	Node_type *pn;
	Count_type *pc;
        Hashtable_type HNN = G->HN;
        Hashtable_type HEL = G->HE;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		pn->extra = MakeCount();
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->vertex != 3 ) continue;
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			pc = pn->extra;
			pc->count++;
			pc->accum += pe->depth;
		}
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
	    pc = pn->extra;
	    if( pc->count > 0 ) {
		pn->depth = pc->accum / pc->count;
	    }
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->vertex != 3 ) continue;
		d=0;
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			d += pn->depth;
		}
		pe->depth = d/3.;
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		free(pn->extra);
	}
}

/***********************************************************************/

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

void ComputeArea( Grid_type *G )

{
	Elem_type *pe;
	double area = 0;
        Hashtable_type HNN = G->HN;
        Hashtable_type HEL = G->HE;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		area += AreaElement(HNN,pe);
	}

	printf("Area of all elements is %f\n",area);
}
