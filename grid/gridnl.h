
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
 *                                                                      *
 * gridnl.h - Number list and Coord routines				*
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
