
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
 * stack.c - stack table administration routines			*
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
