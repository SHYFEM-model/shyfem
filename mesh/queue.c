
/************************************************************************\
 *                                                                      *
 * queue.c - queue administration routines                              *
 *                                                                      *
 * Copyright (c) 1994-1995 by Georg Umgiesser                           *
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
 * 11-Aug-95: Queuetable_type renamed in QueueTable                     *
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

QueueTable MakeQueueTable( void )

{
    QueueTable new;

    new = (QueueTable) malloc( sizeof(QueueTable_head) );

    if( !new ) Error("MakeQueueTable : Cannot allocate QueueTable");

    new->front = NULL;
    new->rear  = NULL;
    new->pact  = NULL;

    return new;
}

void ResetQueueTable( QueueTable L )

{
    if( L == NULL ) NoQueueTable();

	L->pact = NULL;
}

void FreeQueueTable( QueueTable L )

/*\
 *  deletes only QueueTable entries, not the information in info
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

void *VisitQueueTable( QueueTable L )

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

void EnQueue( QueueTable L , void *p )

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

void *DeQueue( QueueTable L )

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

int IsQueueTableEmpty( QueueTable L )

{
	return L->front ? 0 : 1;
}
	
int SizeOfQueueTable( QueueTable L )

{
	int n=0;

	ResetQueueTable(L);
	while( VisitQueueTable(L) )
		n++;

	return n;
}
