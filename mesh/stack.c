
/************************************************************************\ 
 *									*
 * stack.c - stack table administration routines			*
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
 * 11-Aug-95: Stacktable_type renamed in StackTable                     *
 * 25-Jul-95: Stacktable_type created and routines written from scratch *
 *									*
\************************************************************************/

#include <stdlib.h>

#include "general.h"
#include "stack.h"

/*\
 *  Stack table routines (information to insert is a void *)
\*/

static Stack_type *MakeStackPointer( void )

{
        Stack_type *new;

        new = (Stack_type *) malloc( sizeof( Stack_type ) );

        if( !new ) Error("MakeStackPointer : Cannot allocate node");

        new->next = NULL;

        return new;
}

StackTable MakeStackTable( void )

{
	StackTable new;

	new = (StackTable) malloc( sizeof(StackTable_head) );

	if( !new ) Error("MakeStackTable : Cannot allocate StackTable");

	new->head = NULL;
	new->pact = NULL;

        return new;
}

void ResetStackTable( StackTable L )

{
	L->pact = NULL;
}

void FreeStackTable( StackTable L )

/*\
 *  deletes only StackTable entries, not the information in info
\*/

{
	Stack_type *p,*l;

	p = L->head;
	while( p ) {
		l = p->next;
		free(p);
		p=l;
	}
	free(L);
}

void *VisitStackTable( StackTable L )

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

void Push( StackTable L , void *p )

{
        Stack_type *pnew;

        pnew = MakeStackPointer();
	pnew->info = p;
	pnew->next = L->head;
	L->head = pnew;
}

void *Pop( StackTable L )

{
        Stack_type *paux;
	void *pinfo;

	if( L->head ) {
	    paux = L->head;
	    pinfo = paux->info;
	    L->head = paux->next;
	    free(paux);
	} else {
	    pinfo = NULL;
	}

	return pinfo;
}

int IsStackTableEmpty( StackTable L )

{
	return L->head ? 0 : 1;
}

int SizeOfStackTable( StackTable L )

{
	int dim=0;
	Stack_type *p;

	p = L->head;
	while( p ) {
		dim++;
		p = p->next;
	}
	return dim;
}

void *Top( StackTable L )

{
	if( L->head ) {
		return L->head->info;
	} else {
		return NULL;
	}
}

void *NextToTop( StackTable L )

{
	Stack_type *paux;

	if( L->head ) {
		paux = L->head->next;
		if( paux ) {
			return paux->info;
		} else {
			return NULL;
		}
	} else {
		return NULL;
	}
}
