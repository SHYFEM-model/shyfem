
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
 * gridnl.c - Number list and Coord routines				*
 *									*
 * Revision History:							*
 * 02-Dec-95: Number list and Coord routines transferred                *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "fund.h"
#include "list.h"
#include "gridnl.h"


/*\
 *  routines for Coord_type administration
\*/


typedef struct {
	int number;
	Point coord;
} Coord_type;


static Coord_type *MakeCoordPointer( int node , Point *c )

{
	Coord_type *new;

	new = (Coord_type *) malloc( sizeof( Coord_type ) );

	if( !new )
		Error("MakeCoordPointer : Cannot allocate Coord_type structure");

	new->number  = node;
	new->coord.x = c->x;
	new->coord.y = c->y;

        return new;
}

Listtable_type MakeCoord( void )

{
	return MakeListTable();
}

void FreeCoord( Listtable_type L )

{
	Coord_type *cp;

	ResetListTable( L );
	while( (cp=(Coord_type *)VisitListTable(L)) != NULL )
		free(cp);
	FreeListTable( L );
}

void InsertCoord( Listtable_type L , int n , Point *c )

{
	InsertListTable( L , (void *) MakeCoordPointer( n , c ) );
}

void ResetCoord( Listtable_type L )

{
	ResetListTable( L );
}

int VisitCoord( Listtable_type L , Point *c )

{
	Coord_type *cp;

	cp = (Coord_type *) VisitListTable( L );

	if( cp ) {
		*c  = cp->coord;
		return cp->number;
	} else {
		return 0;
	}
}

/*\
 *  routines for Number_list administration
\*/

Number_list_type *MakeNumberList( int node )

{
		Number_list_type *new;

		new = (Number_list_type *) malloc( sizeof( Number_list_type ) );

		if( !new ) Error("MakeNumberList : Cannot allocate Number list");

		new->number = node;
		new->next=NULL;

		return new;
}

void FreeNumberList( Number_list_type *p )

{
		Number_list_type *pp;

		while( p ) {
			pp=p->next;
			free(p);
			p=pp;
		}
}

Number_list_type *FindNumberList( Number_list_type *p , int node )

{
		while( p ) {
			if( p->number == node ) break;
			p = p->next;
		}

		return p;
}

Number_list_type *DeleteNumberList( Number_list_type *p , int node )

{
		Number_list_type *r=NULL,*root=NULL;

		while( p ) {
			if( p->number == node ) break;
			r = p;
			p = p->next;
			if( !root ) root = r;
		}

		if( p ) { /* delete p */
			if( r )
				r->next=p->next;
			else
				root=p->next;
			free(p);
		}

		return root;
}

Number_list_type *InsertNumberList( Number_list_type *p , int node )

{
		Number_list_type *new;

		new = MakeNumberList( node );
		new->next = p;

		return new;
}

