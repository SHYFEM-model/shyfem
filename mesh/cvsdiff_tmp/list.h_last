
/************************************************************************\ 
 *									*
 * list.h - list table administration routines                          *
 *									*
 * Copyright (c) 1994-1995 by Georg Umgiesser				*
 *									*
 * see list.c for copying information					*
 *									*
 * Revision History:							*
 * 11-Aug-95: NumberTable routines transfered to this file              *
 *            Listtable_type renamed to ListTable                       *
 * 17-Apr-94: NumberTable types added                                   *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_LIST_
#define __GUH_LIST_


typedef struct list_tag {
        void *info;
        struct list_tag *next;
} List_type;

typedef struct {
        List_type *head;
        List_type *pact;
} ListTable_head;


typedef ListTable_head *ListTable;

ListTable MakeListTable( void );
void ResetListTable( ListTable L );
void FreeListTable( ListTable L );
void *VisitListTable( ListTable L );
void InsertListTable( ListTable L , void *p );
int IsListTableEmpty( ListTable L );
int SizeOfListTable( ListTable L );


#endif

