
/************************************************************************\
 *                                                                      *
 * queue.h - list table administration routines                         *
 *                                                                      *
 * Copyright (c) 1992-1994 by Georg Umgiesser                           *
 *                                                                      *
 * see queue.c for copying information                                  *
 *                                                                      *
 * Revision History:                                                    *
 * 07-Oct-92: routines written from scratch                             *
 *                                                                      *
\************************************************************************/


#ifndef __GUH_QUEUE_
#define __GUH_QUEUE_


typedef struct queue_tag {
        void *info;
        struct queue_tag *next;
} Queue_type;

typedef struct {
        Queue_type *front;
        Queue_type *rear;
        Queue_type *pact;
} Queuetable_head;


typedef Queuetable_head *Queuetable_type;


Queuetable_type MakeQueueTable( void );
void ResetQueueTable( Queuetable_type L );
void FreeQueueTable( Queuetable_type L );
void *VisitQueueTable( Queuetable_type L );
void InsertQueueTable( Queuetable_type L , void *p );
void *RetrieveQueueTable( Queuetable_type L );


#endif 
