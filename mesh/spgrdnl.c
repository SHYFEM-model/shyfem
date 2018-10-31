
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


#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "grd.h"
#include "grdut.h"
#include "grdhs.h"
#include "spgrdut.h"
#include "spgrdnl.h"

typedef struct nodel_tag {
	int box;
	int nextnode;
	struct nodel_tag *next;
} Nodel_type;

typedef Nodel_type *Nodel_pointer;

typedef struct linel_tag {
	int frombox;
	int tobox;
	int vertex;
	int *index;
	struct linel_tag *next;
} Linel_type;

typedef Linel_type *Linel_pointer;

/****************************************************************/

static Nodel_pointer *Nodel;
static Linel_pointer *Linel;

/****************************************************************/

static Nodel_pointer MakeNodelItem( int box , int nextnode )

{
	Nodel_pointer new;

	new = (Nodel_pointer) malloc( sizeof(Nodel_type) );
	if( !new ) Error("Cannot create NodelItem");

	new->box = box;
	new->nextnode = nextnode;
	new->next = NULL;

	return new;
}

static Linel_pointer MakeLinelItem( int frombox , int tobox 
			, int vertex , int *index )

{
	Linel_pointer new;

	new = (Linel_pointer) malloc( sizeof(Linel_type) );
	if( !new ) Error("Cannot create LinelItem");

	new->frombox = frombox;
	new->tobox = tobox;
	new->vertex = vertex;
	new->index = index;
	new->next = NULL;

	return new;
}

void MakeNodelList( void )

{
	int ntot = GetTotNodes() + 1;
	int bstop = FALSE;
	int i,n;
	int node;
	Grid_type *G = GetGrid();
	Line_type *pl;
	Nodel_pointer pp,ppp;

	/* initialize nodel structure */

	Nodel = (Nodel_pointer *) malloc( ntot * sizeof(Nodel_pointer) );
	if( !Nodel ) Error("Cannot create Nodel array");

	for(i=0;i<ntot;i++) {
	  Nodel[i] = NULL;
	}

	/* insert lines into nodel structure */

	ResetHashTable(G->HL);
	while( (pl=VisitHashTableL(G->HL)) != NULL ) {
		n = pl->vertex;
		for(i=0;i<n-1;i++) {
			node = pl->index[i];
			pp = MakeNodelItem(pl->number,pl->index[i+1]);
			pp->next = Nodel[node];
			Nodel[node] = pp;
		}
	}

	/* check nodel structure */
	/* every line segment may be contained only once */

	for(i=0;i<ntot;i++) {
		pp = Nodel[i];
		while( pp ) {
		    ppp = pp->next;
		    while( ppp ) {
			if( pp->nextnode == ppp->nextnode ) {
			    printf("Error in node %d:",i);
			    printf("multiple line segments to node ");
			    printf("%d\n",pp->nextnode);
			    bstop = TRUE;
			}
			ppp = ppp->next;
		    }
		    pp = pp->next;
		}
	}

	if( bstop == TRUE ) {
		Error("Error in line index.");
	}
}

void PrintNodel( int first )

{
	Nodel_pointer p = Nodel[first];

	printf("Nodel structure for %d\n",first);
	while( p ) {
	  printf(" %d %d\n",p->nextnode,p->box);
	  p = p->next;
	}
}

int GetBoxOfNodel( int first , int second )

{
	Nodel_pointer p = Nodel[first];

	while( p ) {
	  if( p->nextnode == second )
		return p->box;
	  p = p->next;
	}

	return 0;
}

void MakeLinelList( void )

{
	int ntot = GetTotLines() + 1;
	Grid_type *G = GetGrid();
	Line_type *pl;
	Linel_type *pp;
	int i,n;
	int number,nm,nbox;
	int n1,n2,nb;
	int box;
	int bstop = FALSE;
	int debug = FALSE;
	int *ind;
	int *ind2;
	int *boxes,*boxesc;

	/* initialize linel structure */

	Linel = (Linel_pointer *) malloc( ntot * sizeof(Linel_pointer) );
	if( !Linel ) Error("Cannot create Linel array");

	for(i=0;i<ntot;i++) {
	  Linel[i] = NULL;
	}

	/* insert boxes into linel structure */
	/* make list of boxes (boxes) and */
	/* unique list of boxes (boxesc) for check */

	ResetHashTable(G->HL);
	while( (pl=VisitHashTableL(G->HL)) != NULL ) {
		number = pl->number;
		ind = pl->index;
		n = pl->vertex;
		nm = n - 1;

		boxes = (int *) malloc( nm * sizeof(int) );
		boxesc = (int *) malloc( nm * sizeof(int) );
		if( !boxes || !boxesc )
			Error("Cannot allocate aux boxes...");

		for(i=0;i<nm;i++) {
			n1 = ind[i];
			n2 = ind[i+1];
			nb = GetBoxOfNodel(n1,n2);
			if( nb != number ) { /* must be guaranteeded */
				printf("nodes: %d %d\n",n1,n2);
				PrintLine(pl);
				PrintNodel(n1);
				Error("Error in Nodel structure...");
			}
			boxes[i] = GetBoxOfNodel(n2,n1);
		}

		nbox = CheckNeibBoxes(number,boxes,boxesc,nm);

		if( debug ) WriteNeibBoxes(number,boxes,boxesc,nm,nbox);

		if( nbox == 1 ) {
		  printf("Warning: Box %d only one neighbor %d\n",
				number,boxesc[0]);
		} else if( nbox == 0 ) {
		  Error("branch not possible (56)");
		} else if( nbox < 0 ) {
		  bstop = TRUE;
		}

		/* insert lines into linel structure */

		while( (ind2=GetNextBox(nm,boxes,ind,&n,&box)) != NULL ) {
			if( debug ) printf("    %d %d",box,n);

			pp = MakeLinelItem(number,box,n,ind2);
			pp->next = Linel[number];
			Linel[number] = pp;
		}
		if( debug ) printf("\n");

		free(boxes); free(boxesc);
	}

	if( bstop )
		Error("Error in box structure");

	CheckLinel();
}

int CheckNeibBoxes(int number,int *boxes,int *boxesc,int nm)

{
	int boxold = boxes[nm-1];
	int iact = 0;
	int i,j;
	int bstop = FALSE;

	for(i=0;i<nm;i++) {
	  if( boxes[i] != boxold ) {
		boxold = boxes[i];
		boxesc[iact++] = boxold;
	  }
	}

	if( iact == 0 ) {	/* at least one box */
	  boxesc[iact++] = boxold;
	}

	/* check if numbers in boxesc are unique */

	for(i=0;i<iact;i++) {
	  for(j=i+1;j<iact;j++) {
		if( boxesc[i] == boxesc[j] && boxesc[i] != 0 )
			bstop = TRUE;
	  }
	}

	if( bstop ) {
		printf("Error in line %d (repeated boxes)\n",number);
/*
		printf(" total = %d, boxes = %d\n",nm,iact);
		for(i=0;i<nm;i++) {
			printf("%d ",boxes[i]);
		}
		printf("\n");
		for(i=0;i<iact;i++) {
			printf("%d ",boxesc[i]);
		}
		printf("\n");
		exit(1);
*/
		iact = -iact;
	}

	return iact;
}

void WriteNeibBoxes(int box, int *boxes,int *boxesc,int nm, int nbox)

{
	int i;

	printf("Box %d - %d %d\n",box,nm,nbox);
	for(i=0;i<nm;i++) {
		printf(" %d",boxes[i]);
	}
	printf("\n");
	if( nbox < 0 ) nbox = -nbox;
	for(i=0;i<nbox;i++) {
		printf(" %d",boxesc[i]);
	}
	printf("\n");
}

void PrintLine( Line_type *pl )

{
	int i;

	printf("Line %d : %d %d\n",pl->number,pl->type,pl->vertex);
	for(i=0;i<pl->vertex;i++) {
	  printf(" %d",pl->index[i]);
	}
	printf("\n");
}

int *GetNextBox( int nm , int *boxes , int *ind , int *n , int *nbox )

{
	static int iact = 0;
	int box;
	int inext,istart;
	int i,nl;
	int *ind2;

	if( iact >= nm || boxes[iact] >= 0 ) { /* must initialize */
	  iact = nm-1;
	  box = boxes[iact];
	  for(i=0;i<nm;i++) {
	    if( boxes[i] != box ) break;
	    iact = i;
	  }
	}

	/* now we are sure that iact points to the end of a
		line segment, i.e., boxes[iact+1] is the start
		of a new line segment */

	inext = (iact+1) % nm;
	box = boxes[inext];
/*	printf("iact: %d   inext: %d   box: %d\n",iact,inext,box);  */
	if( box < 0 ) {	/* no more boxes left */
	  return NULL;
	}

	istart = inext;
	nl = 0;

	do {
	  nl++;
	  iact = inext;
	  inext = (inext+1) % nm;
	} while( boxes[inext] == box  && inext != istart );

	ind2 = (int *) malloc( (nl+1) * sizeof(int) );
	if( !ind )
		Error("Cannot allocate index of boxes.");

	/* copy index from boxes to ind */

	inext = istart;
	for(i=0;i<nl;i++) {
	  ind2[i] = ind[inext];
	  boxes[inext] = -1;
	  inext = (inext+1) % nm;
	}
	ind2[i] = ind[inext]; /* copy node for end of string */
	
	/* return all relevant information */

	*n = nl+1;
	*nbox = box;

	return ind2;
}

void CheckLinel( void )

{
	int ntot = GetTotLines() + 1;
	int i,j;
	int n;
	int box;
	int bstop = FALSE;
	int *ind;
	Linel_type *pp, *ppp;

	/* all fromboxes must be equal */
	/* line segements must be in right linel structure */

	for(i=0;i<ntot;i++) {
	  pp = Linel[i];
	  while( pp ) {
	    if( pp->frombox != i ) {
		printf("Error in Linel structure: %d %d\n",i,pp->frombox);
		bstop = TRUE;
	    }
	    if( pp->frombox == pp->tobox ) {
		printf("Line between identical boxes %d\n",pp->frombox);
		bstop = TRUE;
	    }
	    n = pp->vertex;
	    ind = pp->index;
	    for(j=1;j<n;j++) {
		box = GetBoxOfNodel(ind[j-1],ind[j]);
		if( box != pp->frombox ) {
		  printf("Error in Linel fromboxes: %d %d\n",box,pp->frombox);
		  bstop = TRUE;
		}
		box = GetBoxOfNodel(ind[j],ind[j-1]);
		if( box != pp->tobox ) {
		  printf("Error in Linel toboxes: %d %d\n",box,pp->tobox);
		  bstop = TRUE;
		}
	    }
	    pp = pp->next;
	  }
	}

	if( bstop )
		Error("Error in CheckLinel.");

	/* we shouln't have double line boxes */

	for(i=0;i<ntot;i++) {
	  pp = Linel[i];
	  while( pp ) {
	    ppp = pp->next;
	    while( ppp ) {
		if( pp->tobox == ppp->tobox && pp->tobox != 0 ) {
		  printf("Line between boxes listed twice: ");
		  printf("%d %d\n",pp->frombox,pp->tobox);
		  bstop = TRUE;
		}
		ppp = ppp->next;
	    }
	    pp = pp->next;
	  }
	}

/*	bstop = FALSE; */

	if( bstop )
		Error("Error in CheckLinel.");
}

int *GetLineOfBox( int frombox , int tobox , int *n )

{
	Linel_type *pp;

	pp = Linel[frombox];

	while( pp ) {
	  if( pp->tobox == tobox ) {
	    *n = pp->vertex;
	    return pp->index;
	  }
	  pp = pp->next;
	}

	*n = 0;
	return NULL;
}

void MakeLines( void )

{
	Line_type *pl;
	int i;
	int frombox,tobox;
	Linel_type *pp;
	int ntot = GetTotLines() + 1;
	Grid_type *G = GetGrid();

	/* initialize nodel structure */

	/* delete old line structure */

        ResetHashTable(G->HL);
        while( (pl=VisitHashTableL(G->HL)) != NULL ) {
		DeleteLine(pl);
	}

	SetTotLines(0);

	for(i=0;i<ntot;i++) {
	  pp = Linel[i];
	  while( pp ) {
		frombox = pp->frombox;
		tobox = pp->tobox;

		if( tobox == 0 ) {
			InsertNewLine(5,pp->vertex,pp->index);
		} else if( frombox < tobox ) {
			InsertNewLine(3,pp->vertex,pp->index);
		}

		pp = pp->next;
	  }
	}
}

void MakeNodes( void )

{
	Line_type *pl;
	Node_type *pn;
	int n,iturn,numb;
	int bstop = FALSE;
	Grid_type *G = GetGrid();
	float *xl, *yl;
	float x,y;
	
	ResetHashTable(G->HN);
	while( (pn=VisitHashTableN(G->HN)) != NULL ) {
		if( pn->use == 0 ) {
			pn->type = 0;
		}
	}

        ResetHashTable(G->HL);
        while( (pl=VisitHashTableL(G->HL)) != NULL ) {
		n = pl->vertex;
		pl->type = 0;
		MakeCoordsFromLine(pl,&xl,&yl);

		ResetHashTable(G->HN);
		while( (pn=VisitHashTableN(G->HN)) != NULL ) {
			if( pn->use ) continue;

			x = pn->coord.x;
			y = pn->coord.y;

			iturn = InClosedLine(n,xl,yl,x,y);

			if( iturn == 1 ) {
			  numb = ROUND( pn->depth );
			  if( numb <= 0 ) {
				printf("Node %d ",pn->number);
				printf("contains no valid box number: ");
				printf(" %d",numb);
				printf(" (Line %d)\n",pl->number);
				bstop = TRUE;
			  } else {
			    if( pl->type != 0 ) {
				printf("Line %d has more than one number: ",
						pl->number); 
				printf("%d %d\n",pl->type,numb);
				bstop = TRUE;
			    } else {
				pl->type = numb;
			    }
			  }
			  pn->type++;
			} else if( iturn < 0 ) { /* all lines anti-clockwise */
				Error("Internal error (563)");
			}
		}

		if( pl->type == 0 ) {
		  printf("No box number found for line %d\n",pl->number);
		  bstop = TRUE;
		}

		free(xl); free(yl);
	}

	ResetHashTable(G->HN);
	while( (pn=VisitHashTableN(G->HN)) != NULL ) {
		if( pn->use == 0 ) {
		  if( pn->type == 0 ) {
			printf("Node %d (depth=%f) not in any box\n",
				pn->number,pn->depth);
			bstop = TRUE;
		  } else if( pn->type > 1 ) {
			printf("Node %d (depth=%f) in more than one box\n",
				pn->number,pn->depth);
			bstop = TRUE;
		  }
		}
	}

	if( bstop )
		Error("Error in MakeNodes.");
}

