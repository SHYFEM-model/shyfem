
/************************************************************************\
 *
 *    Copyright (C) 1995,1997,1999  Georg Umgiesser
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
 * meshut.c - mesh utility routines
 *
 * revision log :
 *
 * 25.07.1995	ggu	routines written from scratch
 * 01.08.1995	ggu	new routines inserted and hash routines transfered
 * 08.10.1997	ggu	routines slighlty restructured
 * ...		ggu	extra structure introduced -> holds use, type, ...
 * ...		ggu	new items are created with type 0, but adeguate mesh type
 * 15.10.1997	ggu	in MakeNode(): set depth to NULLDEPTH (bug)
 * ...		ggu	new routine MakeCoordsFromLine()
 * 16.10.1997	ggu	new routines PrintLineList(), MakeNewLine()
 * ...		ggu	InsertInIndex(), InsertNodeInLine()
 * 12.05.1999	ggu	new routine PrintLine() to print one line
 * 12.05.1999	ggu	bug fix in InsertInIndex() -> insert only once
 *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "general.h"
#include "fund.h"
#include "assert.h"
#include "hash.h"
#include "nlist.h"

#include "mesh.h"
#include "meshut.h"
#include "meshhs.h"
#include "meshge.h"
#include "meshty.h"


/**********************************************************************\
 ********************************************************************** 
\**********************************************************************/


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

void PrintLineList( Hashtable_type H )

{
        Line_type *p;

	ResetHashTable(H);
        while( (p = VisitHashTableL(H)) != NULL ) {
		printf("Line: %d %d %d\n",p->number,p->type,p->vertex);
        }
}

void PrintLines( Hashtable_type H )

{
        Line_type *p;

	ResetHashTable(H);
        while( (p = VisitHashTableL(H)) != NULL ) {
		PrintLine(p);
	}
}

void PrintLine( Line_type *p )

{
	int *ix;
	int nv;

	nv = p->vertex;
	ix = p->index;
	printf("PrintLine: %d %d %d\n",p->number,nv,p->type);
	while( nv-- > 0 )
	    printf(" %d",*ix++);
	printf("\n");
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

void PrintNodeList( char *string , NodeList list )

{
	int i;
	int ntotal;

	if( !list ) return;

	ntotal = list->count;
	printf("%s %d",string,ntotal);

	for(i=0;i<ntotal;i++) {
	    if( i%10 == 0 ) {
		printf("\n");
	    }
	    printf("%d ",list->index[i]);
	}
	printf("\n");
}


/**********************************************************************\
 ********************************************************************** 
\**********************************************************************/

int *MakeIndex( int vertex )

{
	int *new;

	new = (int *) malloc( vertex * sizeof( int ) );
	if( !new ) Error("MakeIndex : Cannot allocate index");

	return new;
}

int *CopyIndex( int vertex , int *index )

{
	int i;
	int *new;

	new = (int *) malloc( vertex * sizeof( int ) );
	if( !new ) Error("CopyIndex : Cannot allocate index");

	for(i=0;i<vertex;i++) {
		new[i] = index[i];
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

int *InsertInIndex( int vertex , int *index , int after , int new )

/* inserts new after node "after", only inserts once (first time found) */

{
	int idiff = 0;
	int i;
	int *newind;

	newind = MakeIndex(vertex+1);

	for(i=0;i<vertex;i++) {
	  newind[i+idiff] = index[i];
	  if( after == index[i] && idiff == 0 ) {
		idiff = 1;
		newind[i+1] = new;
	  }
	}

	free(index);
	return newind;
}

/**********************************************************************\
 ********************************************************************** 
\***********************************************************************/

Node_type *MakeNode( int n , int ntype , Point *c )

{
        Node_type *new;

	new = (Node_type *) malloc( sizeof(Node_type) );
        if( !new ) Error("MakeNode : Cannot allocate node");

	new->extra = (Extra_N_type *) malloc( sizeof(Extra_N_type) );
        if( !new->extra ) Error("MakeNode : Cannot allocate extra structure");

	new->number = n;
	new->type = ntype;
	new->coord.x = c->x;
	new->coord.y = c->y;
	new->depth = NULLDEPTH;

        return new;
}

void DeleteNode( Node_type *p )

{
	p=(Node_type *) DeleteHashByNumber(HNN,p->number);
	ASSERT(p);
	if( p->extra ) free( p->extra);
	free(p);
}

/***********************************************************************/

Elem_type *MakeElemWithIndex( int n , int ntype , int vertex , int *index )

{
	Elem_type *new;

	new = (Elem_type *) malloc( sizeof( Elem_type ) );
	if( !new ) Error("MakeElemWithIndex : Cannot allocate node");

	new->extra = (Extra_E_type *) malloc( sizeof(Extra_E_type) );
        if( !new->extra ) 
		Error("MakeElemWithIndex : Cannot allocate extra structure");

        new->number = n;
        new->vertex = vertex;
	new->type   = ntype;
	new->index  = index;
	new->flag   = FL_DEFAULT;

        return new;
}

Elem_type *MakeElem( int n , int *c , int vertex )

{
	int *index;
	int i;

	index = MakeIndex( vertex );

	for( i=0 ; i<vertex ; i++ )
		index[i] = c[i];

        return MakeElemWithIndex(n,E_NONE,vertex,index);
}

void DeleteElem( Elem_type *p )

{
	p=(Elem_type *) DeleteHashByNumber(HEL,p->number);
	free(p->index);
	free(p);
}

/***********************************************************************/

Line_type *MakeLineWithIndex( int n , int ntype , int vertex , int *index )

{
	Line_type *new;

	new = (Line_type *) malloc( sizeof( Line_type ) );
	if( !new ) Error("MakeLineWithIndex : Cannot allocate node");

	new->extra = (Extra_L_type *) malloc( sizeof(Extra_L_type) );
        if( !new->extra ) 
		Error("MakeLineWithIndex : Cannot allocate extra structure");

        new->number = n;
        new->vertex = vertex;
	new->type   = ntype;
	new->index  = index;

        return new;
}

Line_type *MakeLine( int n , int *c , int vertex )

{
	int *index;
	int i;

	index = MakeIndex( vertex );

	for( i=0 ; i<vertex ; i++ )
		index[i] = c[i];

        return MakeLineWithIndex(n,0,vertex,index);
}

void DeleteLine( Line_type *p )

{
	p=(Line_type *) DeleteHashByNumber(HLI,p->number);
	ASSERT(p);
	free(p->index);
	if( p->extra ) free(p->extra);
	free(p);
}

void InsertNodeInLine( Line_type *pl , int after , int new )

{
	pl->index = InsertInIndex(pl->vertex,pl->index,after,new);
	pl->vertex++;
}

/**********************************************************************\
 ********************************************************************** 
\**********************************************************************/


int InsertNewNode( int type, float x, float y )

{
	Point c;
	Node_type *pn;

	NTotNodes++;

	c.x = x; c.y = y;

	pn = MakeNode( NTotNodes, 0, &c );
	InsertByNodeNumber(HNN,pn);
	SetNtype(pn,type);

	return NTotNodes;
}

int InsertNewElem( int type, int n1, int n2, int n3
			, int nb1, int nb2, int nb3 )

{
	Elem_type *pe;
	int *ind;

	NTotElems++;

	ind = MakeIndex(3);
	ind[0] = n1; ind[1] = n2; ind[2] = n3;

	pe = MakeElemWithIndex(NTotElems,0,3,ind);
	InsertByElemNumber(HEL,pe);
	SetEtype(pe,type);

	ind = MakeIndex(3);
	ind[0] = nb1; ind[1] = nb2; ind[2] = nb3;
	pe->neibor = ind;

	MakeCircumCircle(pe);

	return NTotElems;
}

int InsertNewLine( int type, int vertex, int *index )

{
	Line_type *pl;

	NTotLines++;

	pl=MakeLineWithIndex( NTotLines , type , vertex , index );
	InsertByLineNumber(HLI,pl);
	SetLtype(pl,type);

	return NTotLines;
}

Line_type *MakeNewLine( int type, int vertex, int *index )

{
	Line_type *pl;

	NTotLines++;

	pl=MakeLineWithIndex( NTotLines , 0 , vertex , index );
	SetLtype(pl,type);

	return pl;
}

/***********************************************************************/

int UpdateOldElemByNumber( int elem, int type, int n1, int n2, int n3
			, int nb1, int nb2, int nb3 )

{
	Elem_type *pe;

	pe=RetrieveByElemNumber(HEL,elem);
	return UpdateOldElem(pe,type,n1,n2,n3,nb1,nb2,nb3);
}

int UpdateOldElem( Elem_type *pe, int type, int n1, int n2, int n3
			, int nb1, int nb2, int nb3 )

{
	int *ind;

	ind = pe->index;
	ind[0] = n1; ind[1] = n2; ind[2] = n3;
	pe->index = ind;

	pe->type = 0;
	SetEtype(pe,type);

	ind = pe->neibor;
	ind[0] = nb1; ind[1] = nb2; ind[2] = nb3;
	pe->neibor = ind;

	MakeCircumCircle(pe);

	return pe->number;
}

/***********************************************************************/

int UpdateNeiborByNumber( int elem, int nb1, int nb2, int nb3 )

{
	Elem_type *pe;

	pe=RetrieveByElemNumber(HEL,elem);
	return UpdateNeibor(pe,nb1,nb2,nb3);
}

int UpdateNeibor( Elem_type *pe, int nb1, int nb2, int nb3 )

{
	int *ind;

	ind = pe->neibor;
	ind[0] = nb1; ind[1] = nb2; ind[2] = nb3;
	pe->neibor = ind;

	return pe->number;
}

/***********************************************************************/

void InsertCircumCircle( Elem_type *pe )

{
	int i,npoints;
	int *index;
	float ddeg,rho;
	float x,y,t;
	float x0,y0;
	float pi=3.14159;

	x0=pe->rc.x;
	y0=pe->rc.y;

	InsertNewNode( N_NONE, x0, y0 );

	npoints = 30;
	ddeg = 360./npoints;
	rho = sqrt(pe->rho);

	index = (int *) malloc( (npoints+1)*sizeof(int) );
	if( !index )
		Error("InsertCircumCircle: Cannot allocate index");

	for(i=0;i<npoints;i++) {
                t = i * ddeg * pi / 180.;
                x = x0 + rho * cos(t);
                y = y0 + rho * sin(t);
		index[i] = InsertNewNode( N_NONE, x, y );
        }
	index[npoints] = index[0];

	InsertNewLine( L_NONE, npoints+1, index );
}

	

/**********************************************************************\
 ********************************************************************** 
\**********************************************************************/


void MakeCoordsFromLine( Line_type *pl , float **x , float **y )

{
	int i,node;
	int ncount = pl->vertex;
	float *xp, *yp;
	Node_type *pn;

	xp = (float *) malloc( ncount * sizeof(float) );
	yp = (float *) malloc( ncount * sizeof(float) );

	if( !xp || !yp )
	    Error("MakeCoordsFromLine: Cannot allocate coordinate list");

	for(i=0;i<ncount;i++) {
	    node = pl->index[i];
	    pn = RetrieveByNodeNumber(HNN,node);
	    xp[i] = pn->coord.x;
	    yp[i] = pn->coord.y;
	}

	*x = xp;
	*y = yp;
}

