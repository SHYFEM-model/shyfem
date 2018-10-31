
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
 * list.c - list table administration routines				*
 *									*
 * Revision History:							*
 * 11-Aug-95: NumberTable routines transfered to this file              *
 *            Listtable_type renamed to ListTable                       *
 * 08-Oct-94: ResetListTable initializes pact to NULL (bug)             *
 * 17-Apr-94: NumberTable routines added                                *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 12-Apr-94: Listtable_type created and routines written from scratch  *
 *									*
\************************************************************************/

#include <stdlib.h>

#include "general.h"
#include "list.h"

/*\
 *  List table routines (information to insert is a void *)
\*/

static List_type *MakeListPointer( void )

{
        List_type *new;

        new = (List_type *) malloc( sizeof( List_type ) );

        if( !new ) Error("MakeListPointer : Cannot allocate node");

        new->next = NULL;

        return new;
}

ListTable MakeListTable( void )

{
	ListTable new;

	new = (ListTable) malloc( sizeof(ListTable_head) );

	if( !new ) Error("MakeListTable : Cannot allocate ListTable");

	new->head = NULL;
	new->pact = NULL;

        return new;
}

void ResetListTable( ListTable L )

{
	L->pact = NULL;
}

void FreeListTable( ListTable L )

/*\
 *  deletes only ListTable entries, not the information in info
\*/

{
	List_type *p,*l;

	p = L->head;
	while( p ) {
		l = p->next;
		free(p);
		p=l;
	}
	free(L);
}

void *VisitListTable( ListTable L )

{
        if( L->pact ) {
            L->pact = L->pact->next;
        } else {
            L->pact = L->head;
        }

	if( L->pact ) {
	    return  L->pact->info;
	} else {
	    return  NULL;
	}
}

void InsertListTable( ListTable L , void *p )

{
        List_type *pnew;

        pnew = MakeListPointer();
	pnew->info = p;
	pnew->next = L->head;
	L->head = pnew;
}

int IsListTableEmpty( ListTable L )

{
	return L->head ? 0 : 1;
}

int SizeOfListTable( ListTable L )

{
	int n=0;

	ResetListTable(L);
	while( VisitListTable(L) )
		n++;

	return n;
}
