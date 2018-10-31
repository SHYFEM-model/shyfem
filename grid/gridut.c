
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
 * gridut.c - utility routines						*
 *									*
 * Revision History:							*
 * 18-Feb-2014: new routine Dist2Node()					*
 * 16-Feb-2011: pass type into routines for creation of items		*
 * 16-Jun-2010: new way to compute area of polygon (stable for 64 bit)  *
 * 07-May-1998: type is now integer                                     *
 * 14-Oct-97: DeleteElem/Line(): do not decrement use of nodes          *
 *              must be done befor routines are called                  *
 * 06-Dec-95: MakeNode() without level (eliminated)                     *
 * 04-Dec-95: MakeVect, MakeFloat, DeleteVect, ChangeVect introduced    *
 * 02-Dec-95: Number list and Coord routines to gridnl                  *
 *            error, error2 eliminated                                  *
 * 11-Mar-95: MakeConn transfered to gridhs.c                           *
 *            AreaElement, InvertIndex from gridhs transfered           *
 * 08-May-94: MakeElemWithIndex() and MakeIndex()                       *
 * 08-May-94: NewNode() deleted because never used                      *
 * 13-Apr-94: use new hash routines                                     *
 * 06-Apr-94: in MakeNode depth is now initialized to NULLDEPTH		*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "grid.h"
#include "general.h"
#include "list.h"
#include "gridhs.h"


void PrintNodes( Hashtable_type H )

{
        Node_type *p;
	Point c;
	int n;

	ResetHashTable(H);
        while( (p = VisitHashTableN(H)) != NULL ) {
            c = p->coord;
            n = p->number;
            printf("%d %f %f\n",n,c.x,c.y);
        }
}

void PrintLineS( Hashtable_type HL , Hashtable_type HN , int line , int invert )

{
        Line_type *p;
        Node_type *pn;
	int *ix;
	int i,nv,node;

        p = RetrieveByLineNumber(HL,line);
	if( p ) {
		nv = p->vertex;
		ix = p->index;
		if( invert ) InvertIndex(ix,nv);
		printf("%d %d\n",p->number,p->vertex);
		for(i=0;i<nv;i++) {
		    node = ix[i];
		    pn=RetrieveByNodeNumber(HN,node);
		    printf(" %f %f",pn->coord.x,pn->coord.y);
		    printf("         %d %d\n",pn->number,i+1);
		}
        }
}

void PrintLines( Hashtable_type H )

{
        Line_type *p;
	int *ix;
	int nv;

	ResetHashTable(H);
        while( (p = VisitHashTableL(H)) != NULL ) {
		nv = p->vertex;
		ix = p->index;
		printf("%d",p->number);
		while( nv-- > 0 )
		    printf(" %d",*ix++);
		printf("\n");
        }
}

void PrintElems( Hashtable_type H )

{
        Elem_type *p;
	int *ix;
	int nv;

	ResetHashTable(H);
        while( (p = VisitHashTableE(H)) != NULL ) {
		nv = p->vertex;
		ix = p->index;
		printf("%d",p->number);
		while( nv-- > 0 )
		    printf(" %d",*ix++);
		printf("\n");
        }
}

/****************************************************************/

Vect_type *MakeVect( int total , int actual , float *s , float *d )

{
	Vect_type *new;

	new = (Vect_type *) malloc( sizeof( Vect_type ) );
        if( !new ) Error("MakeVect : Cannot allocate vector");

	if( actual < 1 || actual > total ) actual = 1;

	new->total = total;
	new->actual = actual-1;
	new->speed = s;
	new->dir = d;

	return new;
}

void DeleteVect( Node_type *p )

{
	Vect_type *pv = p->extra;

	free(pv->speed);
	free(pv->dir);
	free(pv);
	p=(Node_type *) DeleteHashByNumber(HVC,p->number);
	free(p);
}

void ChangeVect( Vect_type *p )

{
	p->actual = (p->actual + 1) % p->total;
}

/****************************************************************/

Node_type *MakeNode( int n , int type , Point *c )

{
        Node_type *new;

	new = (Node_type *) malloc( sizeof( Node_type ) );

        if( !new ) Error("MakeNode : Cannot allocate node");

	new->number = n;
	new->type = type;
	new->use = 0;
	new->coord.x = c->x;
	new->coord.y = c->y;
	new->extra = NULL;
	new->depth = NULLDEPTH;

        return new;
}

void DeleteNode( Node_type *p )

{
	p=(Node_type *) DeleteHashByNumber(HNN,p->number);
	free(p);
}

/****************************************************************/

Elem_type *MakeElem( int n , int type , int *c , int vertex )

{
	int *index;
	int i;

	index = MakeIndex( vertex );

	for( i=0 ; i<vertex ; i++ )
		index[i] = c[i];

        return MakeElemWithIndex(n,type,vertex,index);
}


Elem_type *MakeElemWithIndex( int n , int type , int vertex , int *index )

{
	Elem_type *new;

	new = (Elem_type *) malloc( sizeof( Elem_type ) );
	if( !new ) Error("MakeElem : Cannot allocate node");

        new->number = n;
        new->vertex = vertex;
	new->type     = type;
	new->index  = index;
	new->depth    = NULLDEPTH;

        return new;
}

void DeleteElem( Elem_type *p )

{
/*	DeleteUseE(HNN,p); */
	p=(Elem_type *) DeleteHashByNumber(HEL,p->number);
	free(p->index);
	free(p);
}

/****************************************************************/

Line_type *MakeLine( int n , int type , int *c , int vertex )

{
	int *index;
	int i;

	index = MakeIndex( vertex );

	for( i=0 ; i<vertex ; i++ )
		index[i] = c[i];

        return MakeLineWithIndex(n,type,vertex,index);
}

Line_type *MakeLineWithIndex( int n , int type , int vertex , int *index )

{
	Line_type *new;

	new = (Line_type *) malloc( sizeof( Line_type ) );
	if( !new ) Error("MakeLine : Cannot allocate node");

        new->number = n;
        new->vertex = vertex;
	new->type     = type;
	new->index  = index;
	new->depth    = NULLDEPTH;

        return new;
}

void DeleteLine( Line_type *p )

{
/*	DeleteUseL(HNN,p); */
	p=(Line_type *) DeleteHashByNumber(HLI,p->number);
	free(p->index);
	free(p);
}

/****************************************************************/

int *MakeIndex( int vertex )

{
	int *new;

	new = (int *) malloc( vertex * sizeof( int ) );
	if( !new ) {
	  printf("vertex = %d\n",vertex);
	  Error("MakeIndex : Cannot allocate index");
	}
	return new;
}

void InvertIndex( int *index , int nvert )

{
	int iaux,i;

	for(i=0,nvert--;i<nvert;i++,nvert--) {
		iaux=index[i];
		index[i]=index[nvert];
		index[nvert]=iaux;
	}
}

float *MakeFloat( int total )

{
	float *new;

	new = (float *) malloc( total * sizeof( float ) );
	if( !new ) Error("MakeFloat : Cannot allocate float array");
	return new;
}

/****************************************************************/

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

/****************************************************************/

float Dist2Node( Hashtable_type H , int node1 , int node2 )

{
        Point *c1,*c2;
        Node_type *pn;
	float dx,dy;
	
        pn = RetrieveByNodeNumber(H,node1);
        c1 = &pn->coord;
        pn = RetrieveByNodeNumber(H,node2);
        c2 = &pn->coord;

	dx = c1->x - c2->x;
	dy = c1->y - c2->y;

	return dx*dx + dy*dy;
}

/****************************************************************/





