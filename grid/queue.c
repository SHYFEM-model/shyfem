
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
 * queue.c - queue administration routines                              *
 *                                                                      *
 * Revision History:                                                    *
 * 07-Oct-94: Queuetable_type created and routines written from scratch *
 *                                                                      *
\************************************************************************/

#include <stdlib.h>

#include "general.h"
#include "queue.h"

/*\
 *  Queue table routines (information to insert is a void *)
\*/

static Queue_type *MakeQueuePointer( void )

{
    Queue_type *new;

    new = (Queue_type *) malloc( sizeof( Queue_type ) );

    if( !new ) Error("MakeQueuePointer : Cannot allocate node");

    new->next = NULL;

    return new;
}

static void NoQueueTable( void )

{
    Error("No initialized queue present");
}

Queuetable_type MakeQueueTable( void )

{
    Queuetable_type new;

    new = (Queuetable_type) malloc( sizeof(Queuetable_head) );

    if( !new ) Error("MakeQueueTable : Cannot allocate Queuetable");

    new->front = NULL;
    new->rear  = NULL;
    new->pact  = NULL;

    return new;
}

void ResetQueueTable( Queuetable_type L )

{
    if( L == NULL ) NoQueueTable();

	L->pact = NULL;
}

void FreeQueueTable( Queuetable_type L )

/*\
 *  deletes only Queuetable entries, not the information in info
\*/

{
    Queue_type *p,*l;

    if( L == NULL ) NoQueueTable();

    p = L->front;
	while( p ) {
		l = p->next;
		free(p);
		p=l;
	}
	free(L);
}

void *VisitQueueTable( Queuetable_type L )

{
    if( L == NULL ) NoQueueTable();

    if( L->pact ) {
        L->pact = L->pact->next;
    } else {
        L->pact = L->front;
    }

	if( L->pact ) {
	    return  L->pact->info;
	} else {
	    return  NULL;
	}
}

void InsertQueueTable( Queuetable_type L , void *p )

{
    Queue_type *pnew;

    if( L == NULL ) NoQueueTable();

    if( p == NULL )
        return;

    pnew = MakeQueuePointer();
    pnew->info = p;

    if( L->front == NULL ) {
        L->front = pnew;
        L->rear  = pnew;
    } else {
        L->rear->next = pnew;
        L->rear = pnew;
    }
}

void *RetrieveQueueTable( Queuetable_type L )

{
    Queue_type *p;
    void *pinfo;

    if( L == NULL ) NoQueueTable();

    if( L->front == NULL ) {
        return NULL;
    } else {
        p = L->front;
        pinfo = p->info;
        L->front = p->next;
        free(p);
        if( L->front == NULL )
            L->rear = NULL;
        return pinfo;
    }
} 
