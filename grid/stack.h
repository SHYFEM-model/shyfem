
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * stack.h - stack table administration routines
 *
 * revision log :
 *
 * 25.07.1995	ggu	routines written from scratch
 * 11.08.1995	ggu	Stacktable_type renamed in StackTable
 *
\************************************************************************/


#ifndef __GUH_STACK_
#define __GUH_STACK_


typedef struct stack_tag {
        void *info;
        struct stack_tag *next;
} Stack_type;

typedef struct {
        Stack_type *head;
        Stack_type *pact;
} StackTable_head;


typedef StackTable_head *StackTable;


StackTable MakeStackTable( void );
void ResetStackTable( StackTable L );
void FreeStackTable( StackTable L );
void *VisitStackTable( StackTable L );
void Push( StackTable L , void *p );
void *Pop( StackTable L );
int IsStackTableEmpty( StackTable L );
int SizeOfStackTable( StackTable L );
void *Top( StackTable L );
void *NextToTop( StackTable L );

#endif

