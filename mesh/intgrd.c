
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
 * intgrd.c - interpolates depth of one grid onto other grid            *
 *									*
 * Revision History:							*
 * 13-Nov-97: copied from maskgrd.c                                     *
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
#include "backg.h"


void ReadFiles( char *argv , Grid_type *G );
void WriteFile( Grid_type *G );
void ManipulateGrd( Grid_type *GD , Grid_type *GN );

void MakeCM( Grid_type *G , Elem_type *pe, float *xm, float *ym );
void PolyMinMax( Hashtable_type H , Elem_type *pe , Rect *r );
void PolyMinMaxIndex( Hashtable_type H , int nvert , int *index , Rect *r );
int InElement( Hashtable_type H , Elem_type *pe , float x , float y );
Elem_type *FindElement( Grid_type *G, float x, float y );

void MakeElementDepth( Grid_type *G );
float DepthFromElement( Elem_type *pe );
void SetBackGround( Grid_type *G );
float InterpolateDepth( Grid_type *G , Elem_type *pe , float x , float y );



void main(int argc, char *argv[])

{
	Grid_type *GD, *GN;

	GD = MakeGrid();	/* grid with depth values at nodes */
	GN = MakeGrid();	/* new grid that needs interpolation */


	if( argc != 3 ) {
		Error("Usage: intgrd depth-file new-file");
	}

	ReadFiles(argv[1],GD);
	ReadFiles(argv[2],GN);

	MakeElementDepth(GD);
	SetBackGround(GD);
	ManipulateGrd(GD,GN);
	WriteFile(GN);
}


void ReadFiles( char *argv , Grid_type *G )

{
	char sfile[80];
	char *s;

	s=strcpy(sfile,argv);
	s=strcat(s,".grd");
	ReadStandard(s,G);
}

void WriteFile( Grid_type *G )

{
	char sfile[80] = "new.grd";
	char *s;

	s = sfile;
	WriteStandard(s,G);
}


void ManipulateGrd( Grid_type *GD , Grid_type *GN )

{
	float xm,ym;
	Elem_type *pe, *ped;
        Hashtable_type HEL = GN->HE;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		MakeCM(GN,pe,&xm,&ym);
		ped = FindElement(GD,xm,ym);
		pe->depth = DepthFromElement(ped);
		pe->depth = InterpolateDepth(GD,ped,xm,ym);
	}
}

/**********************************************************************/

void MakeCM( Grid_type *G , Elem_type *pe, float *xm, float *ym )

{
        int i;
        Node_type *pn;

        *xm = 0.;
        *ym = 0.;

        for(i=0;i<pe->vertex;i++) {
                pn = RetrieveByNodeNumber(G->HN,pe->index[i]);
                *xm += pn->coord.x;
                *ym += pn->coord.y;
        }

        *xm /= pe->vertex;
        *ym /= pe->vertex;
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

Elem_type *FindElement( Grid_type *G, float x, float y )

{
        Elem_type *pe;

        ResetHashTable(G->HE);
        while( ( pe = VisitHashTableE(G->HE) ) != NULL ) {
                if( InElement( G->HN , pe , x , y ) ) return pe;
        }
        return NULL;
}

void MakeElementDepth( Grid_type *G )

{
	int n,i;
	int *ind;
	float depth = 0.;
	Elem_type *pe;
	Node_type *pn;

        ResetHashTable(G->HE);
        while( ( pe = VisitHashTableE(G->HE) ) != NULL ) {
	  n = pe->vertex;
	  ind = pe->index;

	  for(i=0;i<n;i++) {
	    pn = RetrieveByNodeNumber(G->HN,ind[i]);
	    depth += pn->depth;
	  }

	  if( n ) depth /= n;
	  pe->depth = depth;
        }
}

float DepthFromElement( Elem_type *pe )

{
	if( !pe )
		return NULLDEPTH;
	else
		return pe->depth;
}

void SetBackGround( Grid_type *G )

{
	Elem_type *pe;

        ResetHashTable(G->HE);
        while( ( pe = VisitHashTableE(G->HE) ) != NULL ) {
		pe->extra = (void *) SetBackGrid(G,pe);
	}
}

float InterpolateDepth( Grid_type *G , Elem_type *pe , float x , float y )

{
	int n,i;
	int *ind;
	Backg_type *pb;
	Node_type *pn;
	float d[3];

	n = pe->vertex;
	ind = pe->index;
	pb = (Backg_type *) pe->extra;

	for(i=0;i<n;i++) {
	    pn = RetrieveByNodeNumber(G->HN,ind[i]);
	    d[i] = pn->depth;
	}

	return InterpolBackGrid(pb,x,y,d);
}

