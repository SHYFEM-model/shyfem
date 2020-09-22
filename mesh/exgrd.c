
/************************************************************************\
 *
 *    Copyright (C) 1995,1997-1998,2010,2012,2016  Georg Umgiesser
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
 * exgrd.c - extracts items from grd file
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines written from scratch
 * 05.05.1997	ggu	CompressNumbers() added
 * 17.01.1998	ggu	algorithm changed for UnifyNodes()
 * 08.10.2010	ggu	new routines for purging nodes after unifying
 * 01.02.2012	ggu	bug in choosing lines -> no selection on depth was made
 * 18.03.2016	ggu	-a option also changes versus in line
 * 12.06.2020	ggu	write min/max numbers for items
 *
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

typedef struct {
        int iniuse;
        int use;
} Use_type;


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
void MinMaxNumber( 
			 int *minnode , int *maxnode 
			,int *minelem , int *maxelem 
			,int *minline , int *maxline 
		 );
void MinMaxType( int *mintype , int *maxtype );
void MinMaxVertex( int *minvertex , int *maxvertex );
void MinMaxDepth( float *mindepth , float *maxdepth );
void MinMaxRange( int *minrange , int *maxrange );
void SelectElements( void );
void SelectLines( void );
void SelectNodes( void );
void SelectUnusedNodes( void );
void MakeUse( void );
void SetUse( int ini );
void UnifyNodes( void );
void SubstNodeInElements( Hashtable_type HE , int old , int new );
void SubstNodeInLines( Hashtable_type HL , int old , int new );
void PurgeNodesInElements( Hashtable_type HE, Hashtable_type HN );
void PurgeNodesInLines( Hashtable_type HL, Hashtable_type HN );
int purge_number( int n, int *index );
void CompressNumbers( void );
void MakeAntiClockwise( void );
void DeleteStrangeElements( void );
void WriteInfo(  int minnode, int maxnode
		,int minelem, int maxelem
		,int minline, int maxline
		);

/**************************************************************************/

int main(int argc, char *argv[])

{
	int mintype,maxtype;
	int minnode,maxnode,minelem,maxelem,minline,maxline;
	int minvertex,maxvertex;
	float mindepth,maxdepth;
	int minrange,maxrange;

	G = MakeGrid();
	HNN = G->HN;
	HEL = G->HE;
	HLI = G->HL;
	CM  = G->C;

 	SetOptions(argc,argv);

	ReadFiles(argc,argv);

	MinMaxType(&mintype,&maxtype);
	MinMaxNumber(&minnode,&maxnode,&minelem,&maxelem,&minline,&maxline);
	MinMaxVertex(&minvertex,&maxvertex);
	MinMaxDepth(&mindepth,&maxdepth);
	MinMaxRange(&minrange,&maxrange);

	WriteInfo(minnode,maxnode,minelem,maxelem,minline,maxline);

	if( OpInfo == 1 )		return(0);
	if( OpMaxType == -1 )		OpMaxType = maxtype;
	if( OpMinVertex == -1 )		OpMinVertex = minvertex;
	if( OpMaxVertex == -1 )		OpMaxVertex = maxvertex;
	if( OpMinDepth == NULLDEPTH )	OpMinDepth = mindepth;
	if( OpMaxDepth == NULLDEPTH )	OpMaxDepth = maxdepth;
	if( OpMinRange == -1 )		OpMinRange = minrange;
	if( OpMaxRange == -1 )		OpMaxRange = maxrange;

	MakeUse();
	SetUse(TRUE);

	SelectElements();
	SelectLines();
	/* SetUse(FALSE); */
	SelectNodes();
	SetUse(FALSE);
	SelectUnusedNodes();

	UnifyNodes();
	CompressNumbers();
	MakeAntiClockwise();
	DeleteStrangeElements();

	WriteFile();

	return 0;
}


/**************************************************************************/


void ReadFiles( int argc , char *argv[] )

{
	char sfile[80];
	char *s;

	while( OpArgc < argc ) {
		s=strcpy(sfile,argv[OpArgc++]);
		ReadStandard(s,G,".grd");
	}
}

void WriteFile( void )

{
	char sfile[80] = "new.grd";
	char *s;

	s = sfile;
	WriteStandard(s,G);
}

void MinMaxNumber( 
			 int *minnode , int *maxnode 
			,int *minelem , int *maxelem 
			,int *minline , int *maxline 
		 )

{
	int min;
	int max;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;

	/* find min/max node number */

	min = -1; max = -1;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		min = max = pn->number;
		break;
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( pn->number > max ) max = pn->number;
		if( pn->number < min ) min = pn->number;
	}

	*minnode = min; *maxnode = max;

	/* find min/max elem number */

	min = -1; max = -1;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		min = max = pe->number;
		break;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->number > max ) max = pe->number;
		if( pe->number < min ) min = pe->number;
	}

	*minelem = min; *maxelem = max;

	/* find min/max line number */

	min = -1; max = -1;

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		min = max = pl->number;
		break;
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		if( pl->number > max ) max = pl->number;
		if( pl->number < min ) min = pl->number;
	}

	*minline = min; *maxline = max;
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
	Line_type *pl;

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

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
	    if( pl->depth != NULLDEPTH ) {
		min = max = pl->depth;
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

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
	    if( pl->depth != NULLDEPTH ) {
		if( pl->depth > max ) max = pl->depth;
		if( pl->depth < min ) min = pl->depth;
	    }
	}

	*mindepth = min;
	*maxdepth = max;
}

void MinMaxRange( int *minrange , int *maxrange )

{
	int min=0;
	int max=0;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;

	/* find min/max number range -> there must be at least one node */

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		min = max = pn->number;
		break;
	}

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( pn->number > max ) max = pn->number;
		if( pn->number < min ) min = pn->number;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->number > max ) max = pe->number;
		if( pe->number < min ) min = pe->number;
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		if( pl->number > max ) max = pl->number;
		if( pl->number < min ) min = pl->number;
	}

	*minrange = min;
	*maxrange = max;
}

void SelectElements( void )

{
	int selected;
	Elem_type *pe;
	StackTable delete;

	if( !OpExtractElems && !OpDeleteElems ) return;

	delete = MakeStackTable();

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		selected = 1;
		if( pe->type < OpMinType ) selected = 0;
		if( pe->type > OpMaxType ) selected = 0;
		if( pe->vertex < OpMinVertex ) selected = 0;
		if( pe->vertex > OpMaxVertex ) selected = 0;
		if( pe->depth != NULLDEPTH ) {
		    if( pe->depth < OpMinDepth ) selected = 0;
		    if( pe->depth > OpMaxDepth ) selected = 0;
		}
		if( pe->number < OpMinRange ) selected = 0;
		if( pe->number > OpMaxRange ) selected = 0;
		if( selected && OpDeleteElems )
			Push(delete,(void *)pe);
		if( !selected && OpExtractElems )
			Push(delete,(void *)pe);
	}

	while( (pe=(Elem_type *)Pop(delete)) != NULL ) {
		DeleteElem(pe);
	}

	FreeStackTable(delete);
}

void SelectLines( void )

{
	int selected;
	Line_type *pl;
	StackTable delete;

	if( !OpExtractLines && !OpDeleteLines ) return;

	delete = MakeStackTable();

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		selected = 1;
		if( pl->type < OpMinType ) selected = 0;
		if( pl->type > OpMaxType ) selected = 0;
		if( pl->vertex < OpMinVertex ) selected = 0;
		if( pl->vertex > OpMaxVertex ) selected = 0;
		if( pl->depth != NULLDEPTH ) {
		    if( pl->depth < OpMinDepth ) selected = 0;
		    if( pl->depth > OpMaxDepth ) selected = 0;
		}
		if( pl->number < OpMinRange ) selected = 0;
		if( pl->number > OpMaxRange ) selected = 0;
		if( selected && OpDeleteLines )
			Push(delete,(void *)pl);
		if( !selected && OpExtractLines )
			Push(delete,(void *)pl);
	}

	while( (pl=(Line_type *)Pop(delete)) != NULL ) {
		DeleteLine(pl);
	}

	FreeStackTable(delete);
}

void SelectNodes( void )

{
	int selected;
	Node_type *pn;
	Use_type *pu;
	StackTable delete;

	if( !OpExtractNodes && !OpDeleteNodes ) return;

	delete = MakeStackTable();

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		selected = 1;
		if( pn->type < OpMinType ) selected = 0;
		if( pn->type > OpMaxType ) selected = 0;
		if( pn->depth != NULLDEPTH ) {
		    if( pn->depth < OpMinDepth ) selected = 0;
		    if( pn->depth > OpMaxDepth ) selected = 0;
		}
		if( pn->number < OpMinRange ) selected = 0;
		if( pn->number > OpMaxRange ) selected = 0;
		if( selected && OpDeleteNodes ) {
	    	    pu = (Use_type *)pn->extra;
	    	    if( pu->use == 0 ) 
			Push(delete,(void *)pn);
		}
		if( !selected && OpExtractNodes ) {
	    	    pu = (Use_type *)pn->extra;
	    	    if( pu->use == 0 ) 
			Push(delete,(void *)pn);
		}
	}

	while( (pn=(Node_type *)Pop(delete)) != NULL ) {
		DeleteNode(pn);
	}

	FreeStackTable(delete);
}

void SelectUnusedNodes( void )

/*\
 *  s works on iniuse, S works on use
\*/

{
	Node_type *pn;
	Use_type *pu;
	StackTable delete;

	if( !OpExtractUnusedNodes && !OpDeleteUnusedNodes ) return;

	delete = MakeStackTable();

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( OpExtractUnusedNodes ) {
	    	    pu = (Use_type *)pn->extra;
	    	    if( pu->use == 0 && pu->iniuse != 0 ) 
			Push(delete,(void *)pn);
		}
		if( OpDeleteUnusedNodes ) {
	    	    pu = (Use_type *)pn->extra;
	    	    if( pu->use == 0 ) 
			Push(delete,(void *)pn);
		}
	}

	while( (pn=(Node_type *)Pop(delete)) != NULL ) {
		DeleteNode(pn);
	}

	FreeStackTable(delete);
}

void MakeUse( void )

{
	Node_type *pn;
	Use_type *pu;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		pu = (Use_type *) malloc( sizeof(Use_type) );
		if( !pu )
			Error("MakeUse: Cannot allocate item");
		pu->iniuse = 0;
		pu->use = 0;
		pn->extra = pu;
	}
}

void SetUse( int ini )

{
	int i;
	Elem_type *pe;
	Line_type *pl;
	Node_type *pn;
	Use_type *pu;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		for(i=0;i<pe->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
	   		pu = (Use_type *)pn->extra;
			if( ini ) {
				pu->iniuse++;
			} else {
				pu->use++;
			}
		}
	}

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		for(i=0;i<pl->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pl->index[i]);
	   		pu = (Use_type *)pn->extra;
			if( ini ) {
				pu->iniuse++;
			} else {
				pu->use++;
			}
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

	PurgeNodesInElements(HEL,HNN);
	PurgeNodesInLines(HLI,HNN);

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

void PurgeNodesInElements( Hashtable_type HE, Hashtable_type HN )

{
	Elem_type *pe;
	int new;
	
	ResetHashTable(HE);
	while( (pe=VisitHashTableE(HE)) != NULL ) {
		if( (new=purge_number(pe->vertex,pe->index)) != 0 ) {
			/* remalloc new index list */
			pe->vertex = new;
		}
	}
}

void PurgeNodesInLines( Hashtable_type HL, Hashtable_type HN )

{
	Line_type *pl;
	int new;
	
	ResetHashTable(HL);
	while( (pl=VisitHashTableL(HL)) != NULL ) {
		if( (new=purge_number(pl->vertex,pl->index)) != 0 ) {
			/* remalloc new index list */
			pl->vertex = new;
		}
	}
}

int purge_number( int n, int *index )

{
	int i,old;
	int shift = 0;		/* flag if nodes are to be copied */

	old = 0;
	for(i=1;i<n;i++) {
	    if( index[i] == index[old] ) {
		shift = 1;
	    } else {
		old++;
		if( shift ) {
		    index[old] = index[i];   
		}
	    }
	}

	old++;

	fprintf(stderr,"Purging %d %d %d\n",n,old,shift);

	if( shift ) shift = old;

	return shift;
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

void MakeAntiClockwise( void )

{
        Elem_type *pe;
        Line_type *pl;

	if( !OpMakeAntiClockwise ) return;

        ResetHashTable(HEL);
        while( (pe = VisitHashTableE(HEL)) != NULL ) {
                if( AreaElement( HNN , pe ) < 0 ) {
                        printf("Element %d in clockwise sense: changed.\n",
                                                pe->number);
                        InvertIndex(pe->index,pe->vertex);
                }
        }

        ResetHashTable(HLI);
        while( (pl = VisitHashTableL(HLI)) != NULL ) {
                if( AreaLine( HNN , pl ) < 0 ) {
                        printf("Line %d in clockwise sense: changed.\n",
                                                pl->number);
                        InvertIndex(pl->index,pl->vertex);
                }
        }
}

void DeleteStrangeElements( void )

{
	StackTable delete;
        Elem_type *pe;
	int del,nvert,no,nn,i;
	float tol;

	if( !OpDeleteStrangeElements ) return;

	fprintf(stderr,"Deleting degenerate elements with tol %f\n"
					,OpTollerance);

	delete = MakeStackTable();

	tol = OpTollerance;

        ResetHashTable(HEL);
        while( (pe = VisitHashTableE(HEL)) != NULL ) {
	    del = 0;
            nvert=pe->vertex;
            no = pe->index[nvert-1];
            for(i=0;i<nvert;i++) {
                nn = pe->index[i];
		if( nn == no ) {
                    printf("Element %d with non-unique nodes: %d %d deleted.\n",
                                                pe->number,no,nn);
		    del++;
		}
		no = nn;
	    }
	    if( del == 0 && (ABS(AreaElement( HNN , pe )) <= tol ) ) {
                printf("Element %d with area too small: deleted.\n",
                                                pe->number);
		del++;
	    }

	    if( del > 0 ) Push(delete,(void *)pe);
        }

	while( (pe=(Elem_type *)Pop(delete)) != NULL ) {
		DeleteElem(pe);
	}

	FreeStackTable(delete);
}

/***********************************************************************/

void WriteInfo(  int minnode, int maxnode
		,int minelem, int maxelem
		,int minline, int maxline
		)

{
	if( maxnode != -1 ) {
	  printf("Min/Max Node numbers: %d %d\n",minnode,maxnode);
	}
	if( maxelem != -1 ) {
	  printf("Min/Max Elem numbers: %d %d\n",minelem,maxelem);
	}
	if( maxline != -1 ) {
	  printf("Min/Max Line numbers: %d %d\n",minline,maxline);
	}
}

/***********************************************************************/

