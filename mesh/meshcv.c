
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
 * meshcv.c - routines for convex hull to be used with mesh		*
 *									*
 * Revision History:							*
 * 08-Oct-97: uses new mesh type                                        *
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "nlist.h"
#include "hash.h"
#include "heap.h"
#include "stack.h"

#include "mesh.h"
#include "meshut.h"
#include "meshhs.h"
#include "meshcv.h"
#include "meshck.h"
#include "meshty.h"


void ConvexHull( NodeList hull, NodeList intern )

{
	int i;
	Node_type *pn;

	if(intern->count != hull->count) {
		Error("ConvexHull : size of hull and intern different.");
	}

	i=0;
	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
	    if( IsNtype(pn, N_BOUNDARY ) ) {
		hull->index[i++] = pn->number;
	    }
	}

	if( i != hull->count ) {
		printf("%d %d\n",i,hull->count);
		Error("Internal error in ConvexHull");
	}

	Convex(hull,intern);

}

void Convex( NodeList hull, NodeList intern )

{
	int ntotal;
	int i,j;
	float xmin,ymin;
	float x,y;
	float x1,y1,x2,y2,x3,y3;
	float key,oldkey;
	Node_type *pn;
	void *pv;
	HeapTable HL;
	StackTable Sh,Si;

	ntotal = hull->count;

	/* find min/max */

	pn = RetrieveByNodeNumber(HNN,hull->index[0]);
	xmin = pn->coord.x;
	ymin = pn->coord.y;

	for(i=1;i<ntotal;i++) {
		pn = RetrieveByNodeNumber(HNN,hull->index[i]);
		x = pn->coord.x;
		y = pn->coord.y;
		if( y <= ymin ) {
			if( y < ymin || x < xmin ) {
				xmin = x;
				ymin = y;
			}
		}
	}

	/* populate heap */

	HL = MakeHeapTable( ntotal );
	for(i=0;i<ntotal;i++) {
		pn = RetrieveByNodeNumber(HNN,hull->index[i]);
		HL->entry[i]->key = theta(xmin,ymin,pn->coord.x,pn->coord.y);
		HL->entry[i]->info = (void *) pn;
	}

	/* sort on theta */

	HeapSort( HL );
	/* CheckHeap( HL ); */

	/* for same theta sort on distance */

	j=0;
	oldkey = HL->entry[j]->key;
	key=oldkey;
	for(i=1;i<ntotal;i++) {
		key = HL->entry[i]->key;
		if( key != oldkey ) {
		   if( i-j > 1 ) {
			sortondist(HL,j,i-1,(key-oldkey)/2.,xmin,ymin);
			printf("sortondist called %f\n",oldkey);
		   }
		   oldkey = key;
		   j = i;
		}
	}
	if( i-j > 1 ) {	/* last segment -> sort on negative distance */
		sortondist(HL,j,i-1,(oldkey-key)/2.,xmin,ymin);
		printf("sortondist called %f\n",key);
	}

	/* sort on theta one more time */

	HeapSort( HL );
	/* CheckHeap( HL ); */

	/* do Graham-Scan */

	Sh = MakeStackTable();
	Si = MakeStackTable();

	Push(Sh,(void *) HL->entry[0]);
	Push(Sh,(void *) HL->entry[1]);
	Push(Sh,(void *) HL->entry[2]);

	for(i=3;i<ntotal;i++) {
		for(;;) {
			getcoord(HL->entry[i],&x3,&y3);
			getcoord((HeapItem *) Top(Sh),&x2,&y2);
			getcoord((HeapItem *) NextToTop(Sh),&x1,&y1);
			if( ccw(x1,y1,x2,y2,x3,y3) >= 0 ) break;
			Push(Si,Pop(Sh));
		}
		Push(Sh,(void *) HL->entry[i]);
	}

	/* take from stack and put in arrays */

	i=0;
	while( (pv=Pop(Sh)) != NULL ) {
		pn = (Node_type *) ((HeapItem *)pv)->info;
		hull->index[i++] = pn->number;
	}
	hull->count = i;

	i=0;
	while( (pv=Pop(Si)) != NULL ) {
		pn = (Node_type *) ((HeapItem *)pv)->info;
		intern->index[i++] = pn->number;
	}
	intern->count = i;

	/* invert number list to have proper ordening */

	InvertNodeList(hull);
	InvertNodeList(intern);

	/* clean up */

	FreeHeapTable(HL);
	FreeStackTable(Sh);
	FreeStackTable(Si);
}


void getcoord( HeapItem *item , float *x , float *y )

{
	Node_type *pn;
	
	if(item==NULL) {
		printf("getcoord: Null item\n");
		return;
	}
	pn = (Node_type *) item->info;
	*x = pn->coord.x;
	*y = pn->coord.y;
}

float theta( float x1, float y1, float x2, float y2 )

{
	float dx,dy,ax,ay;
	float t;

	dx = x2 - x1;
	dy = y2 - y1;
	ax = ABS(dx);
	ay = ABS(dy);

	t = (ax+ay==0.) ? 0. : dy/(ax+ay);
	if( dx < 0. ) {
		t = 2. - t;
	} else if( dy < 0. ) {
		t = 4. + t;
	}

	return 90.*t;
}

int ccw( float x1, float y1, float x2, float y2, float x3, float y3 )

/*\
 *  Determines if segment 1-2-3 turns right or left.
 *  +1 left turn (counter-clockwise)
 *  -1 right turn
 *   0 straight line
\*/

{
	float dx1,dx2,dy1,dy2;

	dx1 = x2 - x1;
	dy1 = y2 - y1;
	dx2 = x3 - x1;
	dy2 = y3 - y1;

	if( dx1*dy2 > dy1*dx2 ) return +1;
	if( dx1*dy2 < dy1*dx2 ) return -1;
	if( dx1*dx2 < 0 || dy1*dy2 < 0 ) return -1;
	if( dx1*dx1+dy1*dy1 < dx2*dx2+dy2*dy2 ) return +1;

	return 0;
}

void sortondist( HeapTable HL, int imin, int imax, float dk
			, float xmin, float ymin )

{
	int i;
	float dist=0.;
	float dx,dy,dd;
	Node_type *pn;

	for(i=imin;i<=imax;i++) {
		pn = (Node_type *) HL->entry[i]->info;
		dx = xmin - pn->coord.x;
		dy = ymin - pn->coord.y;
		dd = dx*dx + dy*dy;
		if( dd > dist ) dist = dd;
	}
		
	for(i=imin;i<=imax;i++) {
		pn = (Node_type *) HL->entry[i]->info;
		dx = xmin - pn->coord.x;
		dy = ymin - pn->coord.y;
		dd = dx*dx + dy*dy;
		dd = dk * dd / dist;
		HL->entry[i]->key += dd;
	}
}
		
	
