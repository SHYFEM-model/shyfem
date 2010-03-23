
/************************************************************************\ 
 *									*
 * heap.c - heap table routines						*
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
