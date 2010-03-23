
/************************************************************************\ 
 *									*
 * gridnl.c - Number list and Coord routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
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

