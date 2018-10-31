
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
 * queue.h - list table administration routines                         *
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
