
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
 * heap.h - heap table administration routines                          *
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
