
/************************************************************************\
 *                                                                      *
 * gridnl.h - Number list and Coord routines				*
 *                                                                      *
 * Copyright (c) 1992-1994 by Georg Umgiesser                           *
 *                                                                      *
 * see gridnl.c for copying information                                 *
 *                                                                      *
 * Revision History:                                                    *
 * 02-Dec-95: Number list and Coord routines transferred                *
 * 04-Feb-95: routines written from scratch                             *
 *                                                                      *
\************************************************************************/


#ifndef __GUH_GRIDNL_
#define __GUH_GRIDNL_


#include "fund.h"
#include "list.h"


typedef struct number_list_tag {
	int number;
	struct number_list_tag *next;
} Number_list_type;


Listtable_type MakeCoord( void );
void FreeCoord( Listtable_type L );
void InsertCoord( Listtable_type L , int n , Point *c );
void ResetCoord( Listtable_type L );
int VisitCoord( Listtable_type L , Point *c );

Number_list_type *MakeNumberList( int node );
void FreeNumberList( Number_list_type *p );
Number_list_type *FindNumberList( Number_list_type *p , int node );
Number_list_type *DeleteNumberList( Number_list_type *p , int node );
Number_list_type *InsertNumberList( Number_list_type *p , int node );


#endif
