
/************************************************************************\
 *                                                                      *
 * queue.h - list table administration routines                         *
 *                                                                      *
 * Copyright (c) 1994-1995 by Georg Umgiesser                           *
 *                                                                      *
 * see queue.c for copying information                                  *
 *                                                                      *
 * Revision History:                                                    *
 * 11-Aug-95: Queuetable_type renamed in QueueTable                     *
 * 07-Oct-94: routines written from scratch                             *
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
} QueueTable_head;


typedef QueueTable_head *QueueTable;


QueueTable MakeQueueTable( void );
void ResetQueueTable( QueueTable L );
void FreeQueueTable( QueueTable L );
void *VisitQueueTable( QueueTable L );
void EnQueue( QueueTable L , void *p );
void *DeQueue( QueueTable L );
int IsQueueTableEmpty( QueueTable L );
int SizeOfQueueTable( QueueTable L );


#endif 
