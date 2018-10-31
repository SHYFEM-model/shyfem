
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
 * 01-Jun-2012: better checking of input lines				*
 * 10-Feb-2012: use ABS() to determine line with maximum area		*
 * 06-Dec-2011: better handling of line types				*
 * 07-Oct-2010: better error message in CheckInput()			*
 * 12-Nov-97: check for not unique coordinates in CheckInput()          *
 * 17-Oct-97: in CheckInput() check for couter-clockwise line turning   *
 * 08-Oct-97: uses new mesh type                                        *
 *            check for closed line only for L_EXTERNAL/L_INTERNAL      *
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "meshut.h"
#include "meshhs.h"
#include "meshin.h"
#include "meshge.h"
#include "meshck.h"
#include "meshop.h"
#include "meshty.h"

#include "general.h"
#include "gustd.h"
#include "stack.h"


int CheckInput( void )

/*\
 *  here we should also check for :
 *  - counter-clockwise line turning
 *  - no nodes out of external line
 *  - internal lines are really internal of external line
 *  - lines closed
 *  ...
\*/

{
	int i;
	int ntot=0;
	int next=0;
	int nint=0;
	int nbound=0;
	int turn;
	int ndim;
	int error;
        Node_type *pn;
        Elem_type *pe;
        Line_type *pl;
	StackTable backg;
	float *xe, *ye;
	float area;
	float areamax = 0.;
        Line_type *plmax = NULL;	/* biggest line */
        Line_type *plext = NULL;	/* external line */

/* set all nodes to NONE */

       	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		SetNtype(pn, N_NONE );
	}

/* look for lines */

        ResetHashTable(HLI);
        while( (pl=VisitHashTableL(HLI)) != NULL ) {

	    area = ABS( AreaLine(HNN,pl) );
	    if( area > areamax ) {
		areamax = area;
		plmax = pl;
	    }

	    if( pl->type == L_EXTERNAL ) {	/* HACK *//* FIXME */
		plext = pl;
		SetLtype(pl, L_EXTERNAL );
	    } else if( pl->type == L_INTERNAL ) {
		SetLtype(pl, L_INTERNAL );
	    } else if( pl->type == L_FAULT ) {
		SetLtype(pl, L_FAULT );
	    } else if( pl->type == L_NONE ) {
		if( IsLineClosed(pl) ) {
		    SetLtype(pl, L_INTERNAL );
		} else {
		    SetLtype(pl, L_FAULT );
		}
	    } else {
		  Error2("Unknown line type in line ",itos(pl->number));
	    }

	    ntot++;
	    if( IsLtype(pl, L_EXTERNAL ) ) next++;
	    if( IsLtype(pl, L_INTERNAL ) ) nint++;
        }

/* look for boundary lines */

	if( next > 1 ) {
		Error("More than one line marked external");
	} else if( ntot == 0 ) {
		Warning("No line found: using hull as boundary");
       		ResetHashTable(HNN);
		while( (pn=VisitHashTableN(HNN)) != NULL ) {
			SetNtype(pn, N_BOUNDARY );
		}
	} else if( next == 0 ) {
		Warning("No external line found... using biggest one");
		Warning2("  external line is ",itos(plmax->number));
		plext = plmax;
		SetLtype(plmax, L_EXTERNAL );
		next++;
	}

/* check if external line is closed and all other lines are inside */

	if( plext != NULL ) {
	  if( ! IsLineClosed(plext) ) {
    	    Error2("External line is not closed: ",itos(plext->number));
          }
          ResetHashTable(HLI);
          while( (pl=VisitHashTableL(HLI)) != NULL ) {
	      if( pl == plext ) continue;
	      if( ! IsLineInLine(plext,pl) ) {
	      	Error2("Line not inside external line: ",itos(pl->number));
	      }
	  }
	}

/* check lines and set node type */

	error = 0;
        ResetHashTable(HLI);
        while( (pl=VisitHashTableL(HLI)) != NULL ) {
	    if( IsLtype(pl, L_EXTERNAL ) || IsLtype(pl, L_INTERNAL ) ) {

		/* check for closed line */ /* FIXME */ /* -> close line */

		if( ! IsLineClosed(pl) ) {
	    	    Error2("Line is not closed: ",itos(pl->number));
	        }

		ndim = pl->vertex;
		MakeCoordsFromLine(pl,&xe,&ye);

		/* check for not unique coordinates */

	        for(i=1;i<ndim;i++) {
		  if( pl->index[i-1] == pl->index[i] ) {
		    printf("*** Line %d: identical node %d\n"
				,pl->number,pl->index[i]);
		    error++;
	          } else if( xe[i-1] == xe[i] && ye[i-1] == ye[i] ) {
		    printf("*** Line %d: Identical coordinates - ",pl->number);
		    printf("nodes %d and %d\n",pl->index[i-1],pl->index[i]);
		    error++;
		  }
		}
		
		/* check for counter-clockwise turning */

		turn = TurnClosedLine(pl->vertex,xe,ye);
		if( turn == 1 ) {
			/* ok */;
		} else if( turn == -1 ) {
			printf("Line %d in clockwise sense -> inverting\n"
				,pl->number);
			InvertIndex(pl->index,pl->vertex);
		} else {
			printf("*** Line %d: turning number %d\n",
				pl->number,turn);
			printf("  (the line might turn on itself)\n");
			error++;
		}
		free(xe); free(ye);

		/* set node type for nodes on boundary */

		for(i=0;i<pl->vertex;i++) {
			pn = RetrieveByNodeNumber(HNN,pl->index[i]);
			SetNtype(pn, N_BOUNDARY );
		}
	    }
        }

	if( error ) {
		Error("Errors found in lines");
	}

	/* look for background grid */

	if( OpBackGround >= 0 ) {
	    backg = MakeStackTable();
	    BCK = MakeListTable();
       	    ResetHashTable(HEL);
	    while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( pe->type == OpBackGround ) {
		    Push(backg,(void *)pe);
		} else {
		    SetEtype(pe, E_EXTERNAL ); /* tentativly... */
		}
	    }
	    while( (pe=(Elem_type *)Pop(backg)) != NULL ) {
		DeleteHashByNumber(HEL,pe->number);
		InsertListTable(BCK,(void *)pe);
		for(i=0;i<3;i++) {
		    pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		    SetNtype(pn, N_EXTERNAL );
		    if( pn->depth == NULLDEPTH ) {
			Error2("Null resolution in background grid at node : ",
					itos(pn->number));
		    }
		}
	    }
	    FreeStackTable(backg);
	}

	/* count all nodes of boundary type */

       	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( IsNtype(pn, N_BOUNDARY ) ) nbound++;
	}

	return nbound;
}


void CheckBoundary( void )

{
	int i;
	Node_type *pn;
	Line_type *pl;

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		printf("%d %d ",pn->number,pn->type);
		printf("%f %f\n",pn->coord.x,pn->coord.y);
	}
	printf("\n");

	ResetHashTable(HLI);
	while( (pl=VisitHashTableL(HLI)) != NULL ) {
		printf("%d %d %d\n",pl->number,pl->type,pl->vertex);
		for(i=0;i<pl->vertex;i++) {
		    printf("%d\n",pl->index[i]);
		}
	}
}

void CheckConvex( NodeList hull, NodeList intern )

{
	int i;

	printf("Number of nodes in convex hull : %d\n",hull->count);
	for(i=0;i<hull->count;i++) {
		printf("%d\n",hull->index[i]);
	}
		
	printf("Number of internal nodes : %d\n",intern->count);
	for(i=0;i<intern->count;i++) {
		printf("%d\n",intern->index[i]);
	}
}

void CheckHeapPropertyUP( HeapTable HL, int up )

{
	int left,right;

	if( up >= HL->count ) return;

	left=2*up+1;
	right=left+1;

	if( left >= HL->count ) return;
	if( HL->entry[up]->key > HL->entry[left]->key ) {
		Warning("Error in Heap (left):");
		printf("%d %f %f\n",up,HL->entry[up]->key
			,HL->entry[left]->key);
	}
	if( right >= HL->count ) return;
	if( HL->entry[up]->key > HL->entry[right]->key ) {
		Warning("Error in Heap (right):");
		printf("%d %f %f\n",up,HL->entry[up]->key
			,HL->entry[right]->key);
	}

	CheckHeapPropertyUP(HL,left);
	CheckHeapPropertyUP(HL,right);
}

void CheckHeap( HeapTable HL )

{
	int i;
	int id;
	float key;
	Node_type *pn;

	printf("Heap Table : %d\n",HL->count);
	for(i=0;i<HL->count;i++) {
		pn = (Node_type *) HL->entry[i]->info;
		id = pn->number;
		key = HL->entry[i]->key;
		printf("%d %d %f\n",i,id,key);
	}

	CheckHeapPropertyUP(HL,0);
}

static Elem_type *GoAround( Elem_type *pvis, int node, int dir )

{
	int i,elem;

	for(i=0;i<3;i++) {
		if( pvis->index[i] == node ) break;
	}
	if( i == 3 )
		Error("GoAround: No node found");

	if( dir == RIGHT ) {
		elem = pvis->neibor[(i+2)%3];
	} else {
		elem = pvis->neibor[(i+1)%3];
	}

	return RetrieveByElemNumber(HEL,elem);
}

void CheckCircumCircleProperty( void )

/*\
 *  New version of routine -> does a little bit too much.
\*/

{
	int i,node;
	Elem_type *pe, *pvis;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
	    if( IsEtype(pe, E_EXTERNAL ) ) continue;
	    if( !InElemCircumCircle(pe,pe,1.01) ) {
		Error2("CheckCircumCircleProperty: Violation (0)"
					,itos(pe->number));
	    }
	    for(i=0;i<3;i++) {
		pvis = RetrieveByElemNumber(HEL,pe->neibor[i]);
		node = pe->index[(i+2)%3];
		while( (pvis=GoAround(pvis,node,RIGHT)) != NULL ) {
		    if( pvis == NULL || pvis == pe ) break;
		    if( InElemCircumCircle(pvis,pe,0.99) ) {
			Error2("CheckCircumCircleProperty: Violation (1)"
					,itos(pe->number));
		    }
		}
		pvis = RetrieveByElemNumber(HEL,pe->neibor[i]);
		node = pe->index[(i+1)%3];
		while( (pvis=GoAround(pvis,node,LEFT)) != NULL ) {
		    if( pvis == NULL || pvis == pe ) break;
		    if( InElemCircumCircle(pvis,pe,0.99) ) {
			Error2("CheckCircumCircleProperty: Violation (2)"
					,itos(pe->number));
		    }
		}
	    }
	}
}

void CheckArea( void )

{
	Elem_type *pe;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
	    if( AreaElement(HNN,pe) <= 0. ) {
		Warning2("CheckArea: area is negative in element "
				,itos(pe->number));
	    }
	}
}

void CheckStatus( int id )

{
	int nvis=0;
	int ndel=0;
	int j;
	Elem_type *pe;

	for(j=1;j<=NTotElems;j++) {
	    pe = RetrieveByElemNumber(HEL,j);
	    if( pe ) {
		if( pe->flag == FL_VISITED ) {
			nvis++;
		} else if( pe->flag == FL_DELETED ) {
			ndel++;
		}
	    }
	}
	if( nvis || ndel ) {
		printf("CheckStatus %d : %d %d\n",id,nvis,ndel);
	}
}

void CheckNeibor( int id )

{
	int n0=0;
	int i,j,l;
	int elem,n;
	int *nb;
	int stop=0;
	int warn=0;
	Elem_type *pe, *paux;

#undef DEBUG_CHECK_NEIBOR
#ifdef DEBUG_CHECK_NEIBOR
	for(j=1;j<=NTotElems;j++) {
	    pe = RetrieveByElemNumber(HEL,j);
	    if( pe ) {
		elem = pe->number;
		printf("elem %d - neighbors",elem);
		for(i=0;i<3;i++) {
		    printf(" %d",pe->neibor[i]);
		}
		printf("\n");
	    }
	}
#endif

	for(j=1;j<=NTotElems;j++) {
	    pe = RetrieveByElemNumber(HEL,j);
	    if( pe ) {
		elem = pe->number;
		nb = pe->neibor;
		if( nb == NULL ) continue;
		for(i=0;i<3;i++) {
		    n = pe->neibor[i];
		    if( ACTIVE_NEIBOR(n) ) {
			paux = RetrieveByElemNumber(HEL,n);
			for(l=0;l<3;l++) {
				if( paux->neibor[l] == elem ) break;
			}
			if( l == 3 ) {
			    printf("CheckNeibor: Neighbor error\n");
			    printf("%d %d\n",elem,n);
			    for(l=0;l<3;l++) {
			        printf(" %d",paux->neibor[l]);
			    }
			    printf("\n");
			    stop++;
			}
		    } else if( UNKNOWN_NEIBOR(n) ) {
			printf("Warning: Unknown Neighbor of element %d\n"
					,elem);
			warn++;
		    } else if( DELETED_NEIBOR(n) ) {
			printf("Warning: Deleted Neighbor of element %d\n"
					,elem);
			warn++;
		    } else { /* no neighbour */
			n0++;
		    }
		}
	    }
	}
	if( n0 != 8 || warn || stop ) {
	    printf("CheckNeibor: %d\n",id);
	    printf("No Neighbors: %d , Warnings: %d\n",n0,warn);
	}
	if( stop ) {
		Error("CheckNeibor: Neighbor error");
	}
}

void CheckCircumCircle( void )

{
	Elem_type *pe;

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		ControlCircumCircle(pe);
	}
}

