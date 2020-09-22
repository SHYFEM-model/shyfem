
/************************************************************************\
 *
 *    Copyright (C) 1994-1995  Georg Umgiesser
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
 *
 * nmblist.c - number table administration routines
 *
 * revision log :
 *
 * 12.04.1994	ggu	Listtable_type created and routines written from scratch
 * 13.04.1994	ggu	completely restructured -> independent routines
 * 17.04.1994	ggu	NumberTable routines added
 * 08.10.1994	ggu	ResetListTable initializes pact to NULL (bug)
 * 11.08.1995	ggu	NumberTable routines transfered to this file
 * ...		ggu	Numbertable_type renamed in NumberTable
 *
\************************************************************************/

#include <stdlib.h>

#include "general.h"
#include "nmblist.h"

/*\
 *  Number table routines (information to insert is just an int)
\*/

static Number_type *MakeNumberPointer( void )

{
        Number_type *new;

        new = (Number_type *) malloc( sizeof( Number_type ) );

        if( !new ) Error("MakeNumberPointer : Cannot allocate node");

        new->next = NULL;

        return new;
}

NumberTable MakeNumberTable( void )

{
	NumberTable new;

	new = (NumberTable) malloc( sizeof(NumberTable_head) );

	if( !new ) Error("MakeNumberTable : Cannot allocate NumberTable");

	new->head = NULL;
	new->pact = NULL;

        return new;
}

void ResetNumberTable( NumberTable L )

{
	L->pact = NULL;
}

void FreeNumberTable( NumberTable L )

{
	Number_type *p,*l;

	p = L->head;
	while( p ) {
		l = p->next;
		free(p);
		p=l;
	}
	free(L);
}

int FindNumberTable( NumberTable L , int number )

{
	Number_type *p;

	p = L->head;
	while( p ) {
		if( p->number == number ) return TRUE;
		p = p->next;
	}
	return FALSE;
}

void DeleteNumberTable( NumberTable L , int number )

/*\
 *  deletes one entry from NumberTable
\*/

{
	Number_type *p,*r;

	r = NULL;	/* one step behind */
	p = L->head;
	while( p ) {
		if( p->number == number ) {
			if( r )
				r->next = p->next;
			else
				L->head = p->next;
			free(p);
			break;
		}
		r = p;
		p = p->next;
	}
}

int OkNumberTable( NumberTable L )

{
	if( L->pact && L->pact->next )
		return TRUE;
	else
		return FALSE;
}

int VisitNumberTable( NumberTable L )

/*\
 *  returns 0 at end of list, if 0 is not a possible number then
 *  use of VisitNumberTable is sufficient, else one has
 *  to use OkNumberTable to be sure that more entries are available
\*/

{
        if( L->pact ) {
            L->pact = L->pact->next;
        } else {
            L->pact = L->head;
        }

	if( L->pact ) {
	    return L->pact->number;
	} else {
	    return 0;
	}
}

void InsertNumberTable( NumberTable L , int number )

{
        Number_type *pnew;

        pnew = MakeNumberPointer();
	pnew->number = number;
	pnew->next = L->head;
	L->head = pnew;
}

