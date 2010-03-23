
/************************************************************************\ 
 *									*
 * exgrd.c - extracts items from grd file				*
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
 * 17-Jan-98: algorithm changed for UnifyNodes()                        *
 * 05-May-97: CompressNumbers() added                                   *
 * 17-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "general.h"
#include "gustd.h"

#include "hash.h"
#include "queue.h"
#include "stack.h"

#include "grd.h"
#include "grdio.h"
#include "grdhs.h"
#include "grdut.h"
#include "exgrdop.h"
#include "spgrdnl.h"
#include "spgrdut.h"


/*
	we should be able to eliminate Hashtables and use only Grid_type
	this is just a temporary hack...
*/

Hashtable_type HNN;
Hashtable_type HEL;
Hashtable_type HLI;
QueueTable      CM;
Grid_type       *G;


void ReadFiles( int argc , char *argv[] );
void WriteFile( void );
void MinMaxType( int *mintype , int *maxtype );
void MinMaxVertex( int *minvertex , int *maxvertex );
void MinMaxDepth( float *mindepth , float *maxdepth );
void MakeUse( void );
void UnifyNodes( void );
void SubstNodeInElements( Hashtable_type HE , int old , int new );
void SubstNodeInLines( Hashtable_type HL , int old , int new );
void CompressNumbers( void );

void CheckLines( void );
void WriteLinNum( void );


void main(int argc, char *argv[])

{
	G = MakeGrid();
	HNN = G->HN;
	HEL = G->HE;
	HLI = G->HL;
	CM  = G->C;

 	SetOptions(argc,argv);

	ReadFiles(argc,argv);

	CheckLines();
	MakeCounterClockwise();
	MakeNodelList();
	MakeLinelList();

	MakeUse();
	MakeNodes();
	MakeLines();

	WriteFile();
	WriteLinNum();
}


/**************************************************************************/


void ReadFiles( int argc , char *argv[] )

{
	char sfile[80];
	char *s;

	while( OpArgc < argc ) {
		s=strcpy(sfile,argv[OpArgc++]);
		s=strcat(s,".grd");
		ReadStandard(s,G);
	}
}

void WriteFile( void )

{
	char sfile[80] = "new.grd";
	char *s;

	s = sfile;
	WriteStandard(s,G);
}

void MinMaxType( int *mintype , int *maxtype )

{
	int min=0;
	int max=0;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;

	/* find min/max type -> there must be at least one node */

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		min = max = pn->type;
		break;
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( pn->type > max ) max = pn->type;
		if( pn->type < min ) min = pn->type;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->type > max ) max = pe->type;
		if( pe->type < min ) min = pe->type;
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		if( pl->type > max ) max = pl->type;
		if( pl->type < min ) min = pl->type;
	}

	*mintype = min;
	*maxtype = max;
}

void MinMaxVertex( int *minvertex , int *maxvertex )

{
	int min=0;
	int max=0;
	Elem_type *pe;
	Line_type *pl;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		min = max = pe->vertex;
		break;
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		min = max = pl->vertex;
		break;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->vertex > max ) max = pe->vertex;
		if( pe->vertex < min ) min = pe->vertex;
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		if( pl->vertex > max ) max = pl->vertex;
		if( pl->vertex < min ) min = pl->vertex;
	}

	*minvertex = min;
	*maxvertex = max;
}

void MinMaxDepth( float *mindepth , float *maxdepth )

{
	float min=NULLDEPTH;
	float max=NULLDEPTH;
	Node_type *pn;
	Elem_type *pe;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
	    if( pn->depth != NULLDEPTH ) {
		min = max = pn->depth;
		break;
	    }
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
	    if( pe->depth != NULLDEPTH ) {
		min = max = pe->depth;
		break;
	    }
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
	    if( pn->depth != NULLDEPTH ) {
		if( pn->depth > max ) max = pn->depth;
		if( pn->depth < min ) min = pn->depth;
	    }
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
	    if( pe->depth != NULLDEPTH ) {
		if( pe->depth > max ) max = pe->depth;
		if( pe->depth < min ) min = pe->depth;
	    }
	}

	*mindepth = min;
	*maxdepth = max;
}

void MakeUse( void )

{
	int i;
	Elem_type *pe;
	Line_type *pl;
	Node_type *pn;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		pn->use = 0;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			pn->use++;
		}
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		for(i=0;i<pl->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pl->index[i]);
			pn->use++;
		}
	}
}

void UnifyNodes( void )

{
	int i;
	int ntot;
	float tol;
	float x,y,xx,yy;
	Node_type *pn, *pnode;
	StackTable unify;

	if( !OpUnifyNodes ) return;

	fprintf(stderr,"Unifying nodes with tollerance %f\n",OpTollerance);

	unify = MakeStackTable();

	ntot = GetTotNodes();
	tol = OpTollerance;

	for(i=1;i<=ntot;i++) {

	    pnode = RetrieveByNodeNumber(HNN,i);
	    if( !pnode ) continue;

	    x = pnode->coord.x;
	    y = pnode->coord.y;

	    ResetHashTable(HNN);
	    while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( pn->number <= i ) continue;
	    	xx = pn->coord.x;
	    	yy = pn->coord.y;
		if( xx+tol < x ) continue;
		if( xx-tol > x ) continue;
		if( yy+tol < y ) continue;
		if( yy-tol > y ) continue;
		Push(unify,(void *)pn);
	    }

	    while( (pn=(Node_type *)Pop(unify)) != NULL ) {
		fprintf(stderr,"Unifying %d %d %f %f %f %f\n"
			,pnode->number,pn->number
			,pnode->coord.x,pnode->coord.y
			,pn->coord.x,pn->coord.y);

		SubstNodeInElements(HEL,pn->number,pnode->number);
		SubstNodeInLines(HLI,pn->number,pnode->number);

		DeleteNode(pn);
	    }

	}

	FreeStackTable(unify);
}

void SubstNodeInElements( Hashtable_type HE , int old , int new )

{
	Elem_type *pe;
	int i;
	
	ResetHashTable(HE);
	while( (pe=VisitHashTableE(HE)) != NULL ) {
		for(i=0;i<pe->vertex;i++) {
		    if( pe->index[i] == old ) {
			pe->index[i] = new;
		    }
		}
	}
}

void SubstNodeInLines( Hashtable_type HL , int old , int new )

{
	Line_type *pl;
	int i;

	ResetHashTable(HL);
	while( (pl=VisitHashTableL(HL)) != NULL ) {
		for(i=0;i<pl->vertex;i++) {
		    if( pl->index[i] == old ) {
			pl->index[i] = new;
		    }
		}
	}
}
void CompressNumbers( void )

{
	int i,n,nmax,nnew;
	int *numbers;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;

	if( !OpCompressNumbers ) return;

	/* nodes */

	nmax = 0;
	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		n = pn->number;
		if( n > nmax ) nmax = n;
	}

	numbers = (int *) malloc( (nmax+1) * sizeof(int) );
	if( !numbers )
		Error("Cannot allocate number list for nodes");

	for(n=0;n<=nmax;n++)
		numbers[n] = 0;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		n = pn->number;
		numbers[n] = n;
	}

	nnew = 0;
	for(n=0;n<=nmax;n++) {
	    if( numbers[n] != 0 ) {
		numbers[n] = ++nnew;
	    }
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		n = pn->number;
		pn->number = numbers[n];
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		for(i=0;i<pe->vertex;i++) {
			n = pe->index[i];
			pe->index[i] = numbers[n];
		}
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		for(i=0;i<pl->vertex;i++) {
			n = pl->index[i];
			pl->index[i] = numbers[n];
		}
	}

	free(numbers);

	fprintf(stderr,"%d nodes compressed (max %d)\n",nnew,nmax);

	/* elements */

	nmax = 0;
	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		n = pe->number;
		if( n > nmax ) nmax = n;
	}

	numbers = (int *) malloc( (nmax+1) * sizeof(int) );
	if( !numbers )
		Error("Cannot allocate number list for elements");

	for(n=0;n<=nmax;n++)
		numbers[n] = 0;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		n = pe->number;
		numbers[n] = n;
	}

	nnew = 0;
	for(n=0;n<=nmax;n++) {
	    if( numbers[n] != 0 ) {
		numbers[n] = ++nnew;
	    }
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		n = pe->number;
		pe->number = numbers[n];
	}

	free(numbers);

	fprintf(stderr,"%d elements compressed (max %d)\n",nnew,nmax);

	/* lines */

	nmax = 0;
	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		n = pl->number;
		if( n > nmax ) nmax = n;
	}

	numbers = (int *) malloc( (nmax+1) * sizeof(int) );
	if( !numbers )
		Error("Cannot allocate number list for lines");

	for(n=0;n<=nmax;n++)
		numbers[n] = 0;

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		n = pl->number;
		numbers[n] = n;
	}

	nnew = 0;
	for(n=0;n<=nmax;n++) {
	    if( numbers[n] != 0 ) {
		numbers[n] = ++nnew;
	    }
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		n = pl->number;
		pl->number = numbers[n];
	}

	free(numbers);

	fprintf(stderr,"%d lines compressed (max %d)\n",nnew,nmax);

}

void CheckLines( void )

{
	int n;
	int i,j;
	Line_type *pl;

	ResetHashTable(G->HL);
	while( (pl=VisitHashTableL(G->HL)) != NULL ) {
		n = pl->vertex;

		if( n < 4 ) {
			printf("Line too short: %d - %d\n"
				,pl->number,n);
		}

		if( pl->index[0] != pl->index[n-1] ) {
			printf("Line not closed: %d - %d - %d %d\n"
				,pl->number,n
				,pl->index[0],pl->index[n-1]
				);
		}

		/* check uniqueness: do not check first node */
		/* (is checked as last node) */

		for(i=1;i<n;i++) {
		  for(j=i+1;j<n;j++) {
		    if( pl->index[i] == pl->index[j] ) {
			printf("Node not unique in line: %d - %d\n"
				,pl->number,pl->index[i]);
		    }
		  }
		}
	}
}

void WriteLinNum( void )

{
	int n;
	int i,j;
	int itype;
	int ix,iy,iz;
	int ntot = GetTotLines();
	float x,y;
	Node_type *pn;
	Line_type *pl;
	FILE *fh;

	fh = fopen("new.lin","w");
	for(j=1;j<=ntot;j++) {
		pl = RetrieveByLineNumber(G->HL,j);
		if( !pl ) continue;

		n = pl->vertex;

		itype = 0;
		for(i=0;i<n;i++) {
		   pn = RetrieveByNodeNumber(G->HN,pl->index[i]);
		   x = pn->coord.x;
		   y = pn->coord.y;
		   ix = 10000. * x;
		   iy = 10000. * y;
		   fprintf(fh,"%d   %d %d  0\n",itype,ix,iy);
		   itype = 1;
		}
	}
	fclose(fh);

	MakeUse();

	fh = fopen("new.num","w");
	ResetHashTable(G->HN);
	while( (pn=VisitHashTableN(G->HN)) != NULL ) {
		if( pn->use == 0 ) {
		   x = pn->coord.x;
		   y = pn->coord.y;
		   ix = 10000. * x;
		   iy = 10000. * y;
		   iz = ROUND(pn->depth);
		   fprintf(fh,"%d %d   %d\n",ix,iy,iz);
		}
	}
	fclose(fh);
}
