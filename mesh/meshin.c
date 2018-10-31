
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
 * meshin.c - node insert routines for mesh 				*
 *									*
 * Revision History:							*
 * 15-Oct-97: in GivenNodes: insert given nodes only if internal        *
 * 08-Oct-97: uses new mesh type                                        *
 *            commented sections deleted                                *
 *            new routines InsertNodes(), GivenNodes()                  *
 * 01-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "general.h"

#include "mesh.h"
#include "meshin.h"
#include "meshhs.h"
#include "meshut.h"
#include "meshge.h"
#include "meshck.h"
#include "meshop.h"
#include "meshty.h"

#include "stack.h"

#define DEBUG_MESH
#undef DEBUG_MESH

#define CHECK_MESH
#undef CHECK_MESH

/***********************************************************************/
		
int VisitElems( Elem_type *from , Elem_type *into , StackTable deleted
			, StackTable visited
			, float x , float y )

/*\
 *  Visits elements recursively
 *  never deletes external elements, only visits them
 *  if they were scheduled for deletion flags them in extdel
\*/

{
	int i,n;
	Elem_type *pe;
	int external=0;
	int incircle=0;
	int delete=0;
	static int extdel=0;

	if( from == NULL ) extdel = 0;

	if( !into || into->flag == FL_DELETED ) {
		return extdel;
	}

	if( IsEtype(into, E_EXTERNAL ) ) {
		external = 1;
	}
	if( InCircumCircle(into,x,y) ) {
		incircle = 1;
	}
	if( !from || incircle ) {	/* delete if first or in circle */
		delete = 1;
	}
	if( delete && external ) {	/* if external flag and don't delete */
		delete = 0;
		extdel++;
	}

#undef DEBUG_VISIT_ELEM
#ifdef DEBUG_VISIT_ELEM
	printf("--------------------------\n");
	if( from )
	    printf("from %d %d %d\n",from->number,from->type,from->flag);
	else
	    printf("===========================\n");
	printf("into %d %d %d\n",into->number,into->type,into->flag);
	printf("ext/inc/del/.. %d %d %d %d\n",external,incircle,delete,extdel);
#endif

	if( delete ) {
		into->flag = FL_DELETED;
		Push(deleted,(void *) into);
		for(i=0;i<3;i++) {
			n = into->neibor[i];
			if( ACTIVE_NEIBOR(n) ) {
			    pe=RetrieveByElemNumber(HEL,n);
			    VisitElems(into,pe,deleted,visited,x,y);
			} else {
			    Error("VisitElems: Internal Error (5)");
			}
		}
	} else {
/* if node is external and from == NULL -> bug */ /*FIXME */
		into->flag = FL_VISITED;
		Push(visited,(void *) into);
		for(i=0;i<3;i++) {
			if( from && into->neibor[i] == from->number ) break;
		}
		if( i < 3 ) { /* delete neighbour */
			into->neibor[i] = - into->neibor[i];
		} else {
			Warning("VisitElems: Internal error (2)");
		}
	}

	return extdel;
}

/***********************************************************************/
		
void RefineInternalNodes( void )

{
	int node;
	int ninserted,npass;
	float x,y;
	float facmin;
	float rhoin;
	Elem_type *pe;
	StackTable deleted,visited,created;

	deleted = MakeStackTable();
	visited = MakeStackTable();
	created = MakeStackTable();

	facmin = OpAspect*OpAspect;
	ninserted = 1;
	npass     = 0;

	if( OpAspect <= 0. )
		ninserted = 0;

	while( ninserted ) {
	    ninserted=0;
	    ResetHashTable(HEL);
	    while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( IsEtype(pe, E_EXTERNAL ) ) continue;
		x = pe->rc.x;
		y = pe->rc.y;
		rhoin = MakeInCircleRadius( pe );
		if( pe->rho/rhoin > facmin ) {
#ifdef DEBUG_MESH
		    printf("RefineInternalNodes: %d %d",pe->number,pe->type);
		    printf(" %f %f\n",pe->rho,rhoin);
#endif
	            node = InsertNewNode( N_INTERNAL, x, y );
		    if( InsertNode(node,pe,x,y,deleted,visited,created) ) {
		        ninserted++;
		    } else {
			RecoverElements( deleted , visited );
			DeleteNode(RetrieveByNodeNumber(HNN,node));
		    }
#ifdef CHECK_MESH
		    CheckNeibor(ninserted);
		    CheckStatus(ninserted);
#endif
		}
	    }
	    npass++;
	    printf("Pass %d , inserted %d , total %d\n"
			,npass,ninserted,NTotElems);
	}

	FreeStackTable(visited);
	FreeStackTable(deleted);
	FreeStackTable(created);
}

/***********************************************************************/
		
void InsertInternalNodes( void )

{
	int node;
	int ninserted,npass;
	float x,y;
	Elem_type *pe;
	StackTable deleted,visited,created;

	if( !OpInsertIntern ) return;

	deleted = MakeStackTable();
	visited = MakeStackTable();
	created = MakeStackTable();

	ninserted = 1;
	npass     = 0;

	while( ninserted ) {
	    ninserted=0;
	    ResetHashTable(HEL);
	    while( (pe=VisitHashTableE(HEL)) != NULL ) {
#ifdef DEBUG_MESH
		printf("InsertInternalNodes: %d %d\n",pe->number,pe->type);
#endif
		if( IsEtype(pe, E_EXTERNAL ) ) continue;
		x = pe->rc.x;
		y = pe->rc.y;
#ifdef DEBUG_MESH
		printf("%f %f %f %f\n",x,y,Resol(x,y),pe->rho);
#endif
		if( Resol(x,y) < pe->rho ) {
	            node = InsertNewNode( N_INTERNAL, x, y );
		    if( InsertNode(node,pe,x,y,deleted,visited,created) ) {
		        ninserted++;
		    } else {
			RecoverElements( deleted , visited );
			DeleteNode(RetrieveByNodeNumber(HNN,node));
		    }
#ifdef CHECK_MESH
		    CheckNeibor(ninserted);
		    CheckStatus(ninserted);
#endif
		}
	    }
	    npass++;
	    printf("Pass %d , inserted %d , total %d\n"
			,npass,ninserted,NTotElems);
	}

	FreeStackTable(visited);
	FreeStackTable(deleted);
	FreeStackTable(created);
}

/***********************************************************************/
		
void InsertBoundaryNodes( NodeList nlist )

/* 
	inserts boundary nodes given in nlist 
	is the same as InsertNodes but if no element found exits
	Could be replaced by InsertNodes
	******* look out since here NEL is used instead of HEL *******
*/

{
	int n;
	int node;
	float x,y;
	/* float area; */
	Node_type *pn;
	Elem_type *pe;
	StackTable deleted,visited,created;

	if( !nlist ) return;

	deleted = MakeStackTable();
	visited = MakeStackTable();
	created = MakeStackTable();

	for( n=0 ; n<nlist->count ; n++ ) {	/* go through all nodes */
#ifdef CHECK_MESH
		CheckNeibor(n);
		CheckStatus(n);
#endif
		node = nlist->index[n];
		pn = RetrieveByNodeNumber(HNN,node);
		x = pn->coord.x;
		y = pn->coord.y;
		pe = FindXYElement(NEL,x,y);
/*	printf("Internal node... %d %d %f %f\n",n,node,x,y); */
		if( !pe ) {
			printf("%d %f %f\n",node,x,y);
			Error("Cannot find element to point");
		}
#ifdef DEBUG_MESH
		printf("InsertBP: %d %d %d %f %f\n",node,pn->type
					,pe->number,x,y);
#endif
		/*
		area = AreaElement(HNN,pe);
		printf("area: %f\n",area);
		*/
		InsertNode(node,pe,x,y,deleted,visited,created);
	}
	FreeStackTable(visited);
	FreeStackTable(deleted);
	FreeStackTable(created);
}

/***********************************************************************/
		
void InsertNodes( NodeList nlist )

/*
	inserts a list of nodes into domain
*/

{
	int n;
	int node;
	float x,y;
	Node_type *pn;
	Elem_type *pe;
	StackTable deleted,visited,created;

	if( !nlist ) return;

	deleted = MakeStackTable();
	visited = MakeStackTable();
	created = MakeStackTable();

	for( n=0 ; n<nlist->count ; n++ ) {	/* go through all nodes */

#ifdef CHECK_MESH
		CheckNeibor(n);
		CheckStatus(n);
#endif

		node = nlist->index[n];
		pn = RetrieveByNodeNumber(HNN,node);
		x = pn->coord.x;
		y = pn->coord.y;
		pe = FindXYElement(NEL,x,y);
		/* printf("Internal node... %d %d %f %f\n",n,node,x,y); */
		if( !pe ) {
			printf("%d %f %f ...",node,x,y);
			printf(" No element to node");
			printf(" -> Cannot insert\n");
		} else {
#ifdef DEBUG_MESH
			printf("InsertBP: %d %d %d %f %f\n",node,pn->type
					,pe->number,x,y);
#endif
			InsertNode(node,pe,x,y,deleted,visited,created);
		}
	}
	FreeStackTable(visited);
	FreeStackTable(deleted);
	FreeStackTable(created);
}

/***********************************************************************/

NodeList GivenNodes( void )

/*
	returns nodes read from file but not on boundary
*/

{
	QueueTable intern;
	NodeList given = NULL;
	Node_type *pn;
	Elem_type *pe;
	int ngiven,i;
	int *ind;

	intern = MakeQueueTable();

        ResetHashTable(HNN);
        while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( EqualsNtype(pn, N_NONE ) ) {
			pe = FindXYElement(NEL,pn->coord.x,pn->coord.y);
			if( !IsEtype(pe,E_EXTERNAL) ) {
			    EnQueue(intern,(void *)pn);
			}
		}
	}

	ngiven = SizeOfQueueTable(intern);

	given = MakeNodeList(ngiven);
	ind = given->index;
        given->count = ngiven;

        for(i=0;i<ngiven;i++) {
                pn = (Node_type *) DeQueue(intern);
		ind[i] = pn->number;
        }

	FreeQueueTable(intern);

	return given;
}
	
/***********************************************************************/
		
void RecoverElements( StackTable deleted , StackTable visited )

{
	int i;
	Elem_type *pe;

	while( (pe=(Elem_type *)Pop(visited)) != NULL ) {
		for(i=0;i<3;i++) {
			if( pe->neibor[i] < 0 ) break;
		}
		if( i == 3 ) {
		   Error("RecoverElements: Internal error (1)!");
		}
		pe->neibor[i] = - pe->neibor[i];
		pe->flag = FL_DEFAULT;
	}
	while( (pe=(Elem_type *)Pop(deleted)) != NULL ) {
		pe->flag = FL_DEFAULT;
	}
}

/***********************************************************************/
		
int InsertNode( int node, Elem_type *pe, float x, float y
			, StackTable deleted
			, StackTable visited
			, StackTable created)

{
	int i;
	int extdel;
	int n1,n2;
	int newelem;
	Elem_type *pnew;
	Node_type *pn;

	extdel = VisitElems(NULL,pe,deleted,visited,x,y);

	pn = RetrieveByNodeNumber(HNN,node);
	if( IsNtype(pn, N_INTERNAL ) ) {
	    if( extdel ) {
#ifdef DEBUG_MESH
		printf("Trying to delete external element from node %d\n",node);
#endif
		return 0;
	    }
	}

	while( (pe=(Elem_type *)Pop(visited)) != NULL ) {
		for(i=0;i<3;i++) {
			if( pe->neibor[i] < 0 ) break;
		}
		if( i == 3 ) {
		   Error("InsertNode: Internal error (1)!");
		}
		n1 = pe->index[(++i)%3];
		n2 = pe->index[(++i)%3];
		pnew=(Elem_type *)Pop(deleted);
		if(pnew) {
		   newelem=UpdateOldElem(pnew,E_INTERNAL,n1,node,n2
				,NB_UNKNOWN,pe->number,NB_UNKNOWN);
		} else {
		   newelem=InsertNewElem(E_INTERNAL,n1,node,n2
				,NB_UNKNOWN,pe->number,NB_UNKNOWN);
		   pnew = RetrieveByElemNumber(HEL,newelem);
		   InsertListTable(NEL,(void *) pnew);
		}

		pe->neibor[(++i)%3] = newelem;
		pe->flag = FL_DEFAULT;
		pnew->flag = FL_DEFAULT;

		UpdateNeighbors(pnew,n1,RIGHT,created);
		UpdateNeighbors(pnew,n2,LEFT ,created);
		Push(created,(void *)pnew);
	}
	while( Pop(created) );

	return 1;
}

/***********************************************************************/
		
void UpdateNeighbors( Elem_type *pold, int node, int dir, StackTable L )

/* dir gives direction */

{
	int i,j;
	int new,old;
	Elem_type *pnew;

	if( dir == RIGHT ) {
		new = 1;
		old = 2;
	} else {
		new = 2;
		old = 1;
	}

	ResetStackTable(L);
	while( (pnew=(Elem_type *) VisitStackTable(L)) != NULL ) {
		for(i=0;i<3;i++) {
			if( pnew->index[i] == node ) break;
		}
		if( i < 3 ) {
			for(j=0;j<3;j++) {
				if( pold->index[j] == node ) break;
			}
			if(j==3) Error("Internal error UpdateNeighbors (1)");
			i = (i+new) % 3;
			j = (j+old) % 3;
			pnew->neibor[i] = pold->number;
			pold->neibor[j] = pnew->number;
			break;
		}
	}
}

/***********************************************************************/
		
void SmoothInternalNodes( void )

{
	int node,i,ipass;
	int npass;
	int npoint;
	float omega;
	float *dxsum,*dysum;
	int *ic;
	float x[3],y[3];
	int ntype[3],nnode[3];
	Elem_type *pe;
	Node_type *pn;

	npoint = NTotNodes + 1;
	omega = OpSmoothOmega;
	npass = OpSmoothPass;

	if( npass <= 0 || omega <= 0. ) return;

	dxsum = (float *) malloc( npoint*sizeof(float) );
	dysum = (float *) malloc( npoint*sizeof(float) );
	ic    = (int *) malloc( npoint*sizeof(int) );
	if( !dxsum || !dysum || !ic )
		Error("SmoothInternalNodes: Cannot allocate arrays");

	for(ipass=0;ipass<npass;ipass++) {
	    for(i=0;i<npoint;i++) {
		dxsum[i] = 0.;
		dysum[i] = 0.;
		ic[i] = 0;
	    }

	    /* every internal edge is summed twice */
	    /* and only internal edges are summed  */

	    ResetHashTable(HEL);
	    while( (pe=VisitHashTableE(HEL)) != NULL ) {
		    if( IsEtype(pe, E_EXTERNAL ) ) continue;
		    for(i=0;i<3;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			ntype[i] = GetNtype(pn);
			nnode[i] = pn->number;
			x[i] = pn->coord.x;
			y[i] = pn->coord.y;
		    }
		    for(i=0;i<3;i++) {
			if( ntype[i] & N_INTERNAL ) {
			    node = nnode[i];
			    ic[node] += 2;
			    dxsum[node] += 2.*x[i] - x[(i+1)%3] - x[(i+2)%3];
			    dysum[node] += 2.*y[i] - y[(i+1)%3] - y[(i+2)%3];
			}
		    }
	    }

	    ResetHashTable(HNN);
	    while( (pn=VisitHashTableN(HNN)) != NULL ) {
		node = pn->number;
		if( ic[node] > 0 ) {
		    pn->coord.x -= omega * dxsum[node] / ic[node];
		    pn->coord.y -= omega * dysum[node] / ic[node];
		}
	    }

	    if( ipass %10 == 0 ) {
		if( ipass != 0 ) printf("\n");
		printf("Pass");
	    }
	    printf(" %d",ipass+1);
	    if( ipass+1 == npass ) printf("\n");
	}

	free(dxsum);
	free(dysum);
	free(ic);
}

/***********************************************************************/
		
