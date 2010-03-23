
/* $Id: stack.h,v 1.2 2009-01-14 17:16:37 georg Exp $ */

/************************************************************************\ 
 *									*
 * stack.h - stack table administration routines                        *
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see stack.c for copying information					*
 *									*
 * Revision History:							*
 * 11-Aug-95: Stacktable_type renamed in StackTable                     *
 * 25-Jul-95: routines written from scratch				*
 *									*
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

