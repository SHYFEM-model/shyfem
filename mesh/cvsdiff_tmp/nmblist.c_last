
/************************************************************************\ 
 *									*
 * nmblist.c - number table administration routines			*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
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
 *            Numbertable_type renamed in NumberTable                   *
 * 08-Oct-94: ResetListTable initializes pact to NULL (bug)             *
 * 17-Apr-94: NumberTable routines added                                *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 12-Apr-94: Listtable_type created and routines written from scratch  *
 *									*
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

