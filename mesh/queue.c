
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
 * queue.c - queue administration routines
 *
 * revision log :
 *
 * 07.10.1994	ggu	Queuetable_type created and routines written from scratch
 * 11.08.1995	ggu	Queuetable_type renamed in QueueTable
 *
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
