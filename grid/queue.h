
/************************************************************************\
 *
 *    Copyright (C) 1992  Georg Umgiesser
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
 * queue.h - list table administration routines
 *
 * revision log :
 *
 * 07.10.1992	ggu	routines written from scratch
 *
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
