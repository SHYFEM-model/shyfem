
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
 * heap.c - heap table routines						*
 *									*
 * Revision History:							*
 * 25-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdlib.h>

#include "general.h"
#include "heap.h"


HeapTable MakeHeapTable( int ndim )

/*\
 *  Allocates HeapTable for use with heap algorithms
\*/

{
	HeapTable new;
	int i;

	new = (HeapTable) malloc( sizeof(HeapTable_head) );
	if( !new )
		Error("MakeHeapTable: Cannot allocate HeapTable");
	new->entry = (HeapItem **) malloc( ndim*sizeof(HeapItem *) );
	if( !new->entry )
		Error("MakeHeapTable: Cannot allocate HeapTable entry");

	new->count = ndim;
	new->ndim = ndim;
	for(i=0;i<ndim;i++) {
		new->entry[i] = (HeapItem *) malloc( sizeof(HeapItem) );
		if( !new->entry[i] )
			Error("MakeHeapTable: Cannot allocate Heap item");
	}
	return new;
}

HeapTable ReMakeHeapTable( HeapTable L , int ndim )

/*\
 *  Adjusts HeapTable for use with heap algorithms
\*/

{
	HeapItem **aux;
	int i;

	if( !L ) {
		return MakeHeapTable(ndim);
	} else if( ndim <= L->ndim ) {
		L->count = ndim;
		return L;
	}

	aux = (HeapItem **) malloc( ndim*sizeof(HeapItem *) );
	if( !aux )
		Error("ReMakeHeapTable: Cannot allocate HeapTable entry");

	for(i=0;i<L->ndim;i++) {
		aux[i] = L->entry[i];
	}

	free(L->entry);

	for(i=L->ndim;i<ndim;i++) {
		aux[i] = (HeapItem *) malloc( sizeof(HeapItem) );
		if( !aux[i] )
			Error("ReMakeHeapTable: Cannot allocate Heap item");
	}

	L->entry = aux;
	L->ndim = ndim;
	L->count = ndim;

	return L;
}

void FreeHeapTable( HeapTable L )

/*\
 *  Frees HeapTable. Info in HeapItem's must be freed or saved elsewhere.
\*/

{
	int i;

	for(i=0;i<L->ndim;i++) {
		free(L->entry[i]);
	}
	free(L->entry);
	free(L);
}

int SizeOfHeapTable( HeapTable L )

{
	return L->count;
}

void InsertHeap( HeapTable L, HeapItem *item, int start, int maxheap )

/*\
 *  Insert item in partial heap with empty root at start
 *  The heap is between start and maxheap
\*/

{
	int m;

	if( (m=2*start) <= 0 )	/* if start=0 -> child at 1 */
		m=1;

	while( m<=maxheap ) {
		if( m<maxheap && L->entry[m]->key < L->entry[m+1]->key ) {
			m++;	/* m contains index of larger key */
		}
		if( item->key >= L->entry[m]->key ) {
			break;	/* item belongs in position start */
		} else {
			L->entry[start] = L->entry[m];
			start = m;
			m = 2*start;
		}
	}

	L->entry[start] = item;
}

void BuildHeap( HeapTable L )

/*\
 *  Build heap from contigous list
\*/

{
	int i;

	for(i=L->count/2;i>=0;i--) {
		InsertHeap(L,L->entry[i],i,L->count-1);
	}
}

void HeapSort( HeapTable L )

/*\
 *  Use heap sort method to sort a contigous list
\*/

{
	int i;
	HeapItem *item;

	BuildHeap(L);
	for (i=L->count-1;i>=1;i--) {
		item = L->entry[i];         /* extract last from list */
		L->entry[i] = L->entry[0]; /* move top to end of list */
		InsertHeap(L,item,0,i-1);   /* restore heap property */
	}
}
