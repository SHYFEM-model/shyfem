
/************************************************************************\
 *									*
 * meshck.h - check routines for mesh					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshck.c for copying information					*
 *									*
 * Revision History:							*
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_MESHCK_
#define __GUH_MESHCK_


#include "nlist.h"
#include "heap.h"

int CheckInput( void );
void CheckBoundary( void );
void CheckConvex( NodeList hull, NodeList intern );
void CheckHeapPropertyUP( HeapTable HL, int up );
void CheckHeap( HeapTable HL );
void CheckCircumCircleProperty( void );
void CheckArea( void );
void CheckStatus( int id );
void CheckNeibor( int id );
void CheckCircumCircle( void );


#endif


