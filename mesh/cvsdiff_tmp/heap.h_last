
/************************************************************************\ 
 *									*
 * heap.h - heap table administration routines                          *
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see heap.c for copying information					*
 *									*
 * Revision History:							*
 * 11-Aug-95: HeapList renamed in HeapTable                             *
 *            in HeapItem *item renamed to *info                        *
 * 25-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_HEAP_
#define __GUH_HEAP_


/*\
 *  In HeapTable not the Items but pointer to items are stored.
 *  Only pointers are moved around, never the items.
\*/


/*
#define Parent(i) ( (i) >> 1 )
#define Left(i)   ( (i) << 1 )
#define Right(i)  ( ( (i) << 1 ) | 1 )
*/


typedef struct {
        float key;
        void *info;
} HeapItem;

typedef struct {
        int count;		/* filling of heap table */
	int ndim;		/* dimension of heap table */
	HeapItem **entry;
} HeapTable_head;

typedef HeapTable_head *HeapTable;


HeapTable MakeHeapTable( int ndim );
HeapTable ReMakeHeapTable( HeapTable L , int ndim );
void FreeHeapTable( HeapTable L );
int SizeOfHeapTable( HeapTable L );
void InsertHeap( HeapTable L, HeapItem *item, int start, int maxheap );
void BuildHeap( HeapTable L );
void HeapSort( HeapTable L );


#endif
