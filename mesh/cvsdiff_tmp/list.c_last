
/************************************************************************\ 
 *									*
 * list.c - list table administration routines				*
 *									*
 * Copyright (c) 1994-1995 by Georg Umgiesser				*
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
