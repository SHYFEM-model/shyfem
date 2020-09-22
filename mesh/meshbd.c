
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
 * meshbd.c - boundary routines for mesh
 *
 * revision log :
 *
 * 01.08.1995	ggu	routines written from scratch
 * 08.10.1997	ggu	uses new mesh type
 * 16.10.1997	ggu	in RefineBoundary() use always new line and delete old
 * ...		ggu	new routines RecoverBoundaryNodes(), FindElemToSide()
 * ...		ggu	MakeMidPoint()
 * 12.11.1997	ggu	in RecoverBoundaryNodes() interpolate also depth
 * ...		ggu	-> MakeMidPoint() also returns interpolated depth
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general.h"

#include "mesh.h"
#include "meshbd.h"
#include "meshhs.h"
#include "meshut.h"
#include "meshop.h"
#include "meshty.h"

#include "meshge.h"
#include "meshin.h"

#include "queue.h"
#include "stack.h"
#include "nlist.h"
#include "heap.h"



typedef struct {
	float x;
	float y;
	float t;
	float d;
	int n;
} Refine_type;

static Refine_type *MakeRefineItem( int n, float d, float x, float y, float t );

int TotBoundNodes( void )

{
	int ntotal=0;
	Line_type *pl;

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
            if( IsLtype(pl, L_EXTERNAL | L_INTERNAL ) ) {
		ntotal += pl->vertex - 1;
	    }
	}
	return ntotal;
}


NodeList RefineBoundary( void )

{
	int ntotal=0;
	int nbound;
	int i,j;
	int *ind;
	Line_type *pl,*plnew;
	Node_type *pn1, *pn2;
	float x1,x2,y1,y2;
	QueueTable line,bound,totbound,newlines;
	NodeList refined;
	Refine_type *pr;
	void *pv;

        if( !OpRefineBoundary ) {
		return NULL;
	}

	line = MakeQueueTable();
	bound = MakeQueueTable();
	totbound = MakeQueueTable();
	newlines = MakeQueueTable();

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
            if( IsLtype(pl, L_EXTERNAL_REF | L_INTERNAL_REF ) ) {
		ntotal = pl->vertex;
		for(i=1;i<ntotal;i++) {
		    pn1=RetrieveByNodeNumber(HNN,pl->index[i-1]);
		    pn2=RetrieveByNodeNumber(HNN,pl->index[i]);
		    x1=pn1->coord.x;
		    y1=pn1->coord.y;
		    x2=pn2->coord.x;
		    y2=pn2->coord.y;
		    RefineLine(line,x1,y1,x2,y2,0.);
		    SortLine(line);
		    pr=MakeRefineItem(pn1->number,0.,0.,0.,-1.);
		    EnQueue(bound,(void *)pr);
		    while( (pv=DeQueue(line)) != NULL ){
			pr = (Refine_type *) pv;
			pr->n = InsertNewNode(N_ADDBOUND,pr->x,pr->y);
			EnQueue(bound,(void *)pr);
		    }
		}
		pn2 = RetrieveByNodeNumber(HNN,pl->index[ntotal-1]);
		pr=MakeRefineItem(pn2->number,0.,0.,0.,-1.);
		EnQueue(bound,(void *)pr);
		nbound = SizeOfQueueTable(bound);
/*	printf("nbound %d\n",nbound); */
		ind = MakeIndex( nbound );
		for(i=0;i<nbound;i++) {
		    pr = (Refine_type *) DeQueue(bound);
		    ind[i] = pr->n;
/*	printf("%d %f %f\n",pr->n,pr->t,pr->d); */
		    EnQueue(totbound,(void *)pr);
		}
		/* we cannot insert directly here because of loop over pl */
                if( IsLtype(pl, L_EXTERNAL_REF ) ) {
		    plnew = MakeNewLine(L_EXTERNAL_REF,nbound,ind);
		    SetLtype(pl, L_EXTERNAL );
		} else {
		    plnew = MakeNewLine(L_INTERNAL_REF,nbound,ind);
		    SetLtype(pl, L_INTERNAL );
		}
		EnQueue(newlines,(void *)plnew);
/*	printf("Line created : %d\n",NTotLines); */
	    }
	}

	while( (pl = (Line_type *)DeQueue(newlines)) != NULL ) {
		InsertByLineNumber(HLI,pl);
	}

	nbound = SizeOfQueueTable(totbound);
	refined = MakeNodeList( nbound );
	ind = refined->index;

	j=0;
	for(i=0;i<nbound;i++) {
		pr = (Refine_type *) DeQueue(totbound);
		if( pr->t != -1. ) {
			ind[j++] = pr->n;
		}
		free(pr);
	}
	refined->count = j;

	FreeQueueTable(line);
	FreeQueueTable(bound);
	FreeQueueTable(totbound);
	FreeQueueTable(newlines);

	return refined;
}


void RefineLine( QueueTable line, float x1, float y1
			, float x2, float y2 , float dist )

/*\
 *  dist is distance of first point from start of total line segment
\*/

{
	float fact;
	float dx,dy,dl;
	float xm,ym,t;
	float dd;
	Refine_type *pr;

	fact=1.7;
	fact=3.;	/* check this */

	dx=x2-x1;
	dy=y2-y1;
	dl=dx*dx+dy*dy;
	xm=0.5*(x1+x2);
	ym=0.5*(y1+y2);
	t=fact*Resol(xm,ym)/dl;
	/*printf("%f %f %f %f %f\n",xm,ym,dl,t,Resol(xm,ym)); */
	if( t < 1. ) {
		dd = 0.5*sqrt(dl);
		pr = MakeRefineItem(0,dist+dd,xm,ym,t);
		EnQueue(line,(void *)pr);
		RefineLine(line,x1,y1,xm,ym,dist);
		RefineLine(line,xm,ym,x2,y2,dist+dd);
	}
}

static Refine_type *MakeRefineItem( int n, float d, float x, float y, float t )

{
	Refine_type *pr;

	pr = (Refine_type *) malloc( sizeof(Refine_type) );
	if( !pr )
		Error("MakeRefineItem: Cannot allocate item");

	pr->n = n;
	pr->d = d;
	pr->x = x;
	pr->y = y;
	pr->t = t;

	return pr;
}

void SortLine( QueueTable list )

{
	int i;
	int ntotal;
	HeapTable HL;
	Refine_type *pr;

	ntotal = SizeOfQueueTable(list);
	if( ntotal <= 0 ) return;

        HL = MakeHeapTable( ntotal );

        for(i=0;i<ntotal;i++) {
		pr = (Refine_type *) DeQueue(list);
                HL->entry[i]->key = pr->d;
                HL->entry[i]->info = (void *) pr;
        }

	HeapSort( HL );

	for(i=0;i<ntotal;i++) {
		EnQueue(list,HL->entry[i]->info);
	}

	FreeHeapTable( HL );
}

void CopyBoundaryLine( void ) 

{
	int *ind;
	Line_type *pl;

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
	    /* printf("CopyBoundaryLine: %d %d\n",pl->number,pl->type); */
            if( IsLtype(pl, L_EXTERNAL ) ) {
		ind = CopyIndex(pl->vertex,pl->index);
		InsertNewLine(L_EXTERNAL_REF,pl->vertex,ind);
	    } else if( IsLtype(pl, L_INTERNAL ) ) {
		ind = CopyIndex(pl->vertex,pl->index);
		InsertNewLine(L_INTERNAL_REF,pl->vertex,ind);
	    } else if( IsLtype(pl, L_FAULT ) ) {
		ind = CopyIndex(pl->vertex,pl->index);
		InsertNewLine(L_FAULT_REF,pl->vertex,ind);
	    }
	}
}

void MakeMidPoint( int node1 , int node2 , float *x , float *y , float *depth );
	/* HACK */

void RecoverBoundary( void )

{
	RecoverBoundaryNodes( L_EXTERNAL_REF | L_INTERNAL_REF );
}

void RecoverInternalFault( void )

{
	RecoverBoundaryNodes( L_FAULT_REF );
}

void RecoverBoundaryNodes( int linetype )

{
	Elem_type *pe;
	Line_type *pl;
	Node_type *pn;
	int n,i,ncount,pass;
	int ok,changed;
	int changes;
	int *ind;
	float x,y;
	float depth;
        StackTable deleted,visited,created;

        deleted = MakeStackTable();
        visited = MakeStackTable();
        created = MakeStackTable();

	do {

	changes = FALSE;	/* overall changes -> repeat one more time */

        ResetHashTable(HLI);
        while( (pl=VisitHashTableL(HLI)) != NULL ) {
	    /* PrintLine(pl); */
            if( IsLtype(pl,linetype) ) {
		pass = 0;
		do {
		  /* if( pl->number == 108 ) PrintLine(pl); */
		  changed = FALSE;
		  ncount = pl->vertex;
		  ind = pl->index;
		  pass++;
		  if( pass > 1 ) printf("Entering pass %d for line %d\n",
						pass,pl->number);

		  for(i=1;i<ncount;i++) {
		    pe = FindElemToSide(ind[i-1],ind[i]);
		    if( !pe ) {
			printf("Recovering side %d %d of line %d (pass %d)\n"
				,ind[i-1],ind[i],pl->number,pass);
			MakeMidPoint(ind[i-1],ind[i],&x,&y,&depth);
			n = InsertNewNode(N_ADDBOUND,x,y);
			pn = RetrieveByNodeNumber(HNN,n);
			pn->depth = depth;
			InsertNodeInLine(pl,ind[i-1],n);
			pe = FindElement(HEL,x,y);
			ok = InsertNode(n,pe,x,y,deleted,visited,created);
			if( !ok ) Error("Internal error RecoverBoundaryNodes");
			changed = TRUE;
			changes = TRUE;
			break;
		    }
		  }
		} while( changed );

	    }
	}

	if( changes ) printf("Revisiting all lines...\n");

	} while( changes );

	printf("                                           ...done\n");

        FreeStackTable(visited);
        FreeStackTable(deleted);
        FreeStackTable(created);
}

Elem_type *FindElemToSide( int node1 , int node2 )

{
	int i,nold,nnew;
	Elem_type *pe;

        ResetHashTable(HEL);
        while( (pe=VisitHashTableE(HEL)) != NULL ) {
	    nnew = pe->index[2];
	    for(i=0;i<3;i++) {
		nold = nnew;
	        nnew = pe->index[i];
		if( node1 == nold && node2 == nnew ) return pe;
		if( node1 == nnew && node2 == nold ) return pe;
	    }
	}

	return NULL;
}

void MakeMidPoint( int node1 , int node2 , float *x , float *y , float *depth )

{
	Node_type *pn1 = RetrieveByNodeNumber(HNN,node1);
	Node_type *pn2 = RetrieveByNodeNumber(HNN,node2);

	*x = 0.5 * ( pn1->coord.x + pn2->coord.x );
	*y = 0.5 * ( pn1->coord.y + pn2->coord.y );

	if( pn1->depth == NULLDEPTH || pn2->depth == NULLDEPTH ) {
		*depth = NULLDEPTH;
	} else {
		*depth = 0.5 * ( pn1->depth + pn2->depth );
	}
}
