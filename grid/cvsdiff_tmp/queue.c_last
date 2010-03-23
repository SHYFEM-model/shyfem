
/************************************************************************\
 *                                                                      *
 * queue.c - queue administration routines                              *
 *                                                                      *
 * Copyright (c) 1994 by Georg Umgiesser                                *
 *                                                                      *
 * Permission to use, copy, modify, and distribute this software        *
 * and its documentation for any purpose and without fee is hereby      *
 * granted, provided that the above copyright notice appear in all      *
 * copies and that both that copyright notice and this permission       *
 * notice appear in supporting documentation.                           *
 *                                                                      *
 * This file is provided AS IS with no warranties of any kind.          *
 * The author shall have no liability with respect to the               *
 * infringement of copyrights, trade secrets or any patents by          *
 * this file or any part thereof.  In no event will the author          *
 * be liable for any lost revenue or profits or other special,          *
 * indirect and consequential damages.                                  *
 *                                                                      *
 * Comments and additions should be sent to the author:                 *
 *                                                                      *
 *          Georg Umgiesser                                             *
 *          ISDGM/CNR                                                   *
 *          S. Polo 1364                                                *
 *          30125 Venezia                                               *
 *          Italy                                                       *
 *                                                                      *
 *          Tel.   : ++39-41-5216875                                    *
 *          Fax    : ++39-41-2602340                                    *
 *          E-Mail : georg@lagoon.isdgm.ve.cnr.it                       *
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
