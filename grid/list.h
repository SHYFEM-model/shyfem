
/************************************************************************\ 
 *									*
 * list.h - list table administration routines                          *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see list.c for copying information					*
 *									*
 * Revision History:							*
 * 17-Apr-94: NumberTable types added                                   *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_LIST_
#define __GUH_LIST_


typedef struct number_tag {
        int number;
        struct number_tag *next;
} Number_type;

typedef struct {
        Number_type *head;
        Number_type *pact;
} Numbertable_head;

typedef struct list_tag {
        void *info;
        struct list_tag *next;
} List_type;

typedef struct {
        List_type *head;
        List_type *pact;
} Listtable_head;


typedef Listtable_head *Listtable_type;
typedef Numbertable_head *Numbertable_type;


Listtable_type MakeListTable( void );
void ResetListTable( Listtable_type L );
void FreeListTable( Listtable_type L );
void *VisitListTable( Listtable_type L );
void InsertListTable( Listtable_type L , void *p );

Numbertable_type MakeNumberTable( void );
void ResetNumberTable( Numbertable_type L );
void FreeNumberTable( Numbertable_type L );
int FindNumberTable( Numbertable_type L , int number );
void DeleteNumberTable( Numbertable_type L , int number );
int OkNumberTable( Numbertable_type L );
int VisitNumberTable( Numbertable_type L );
void InsertNumberTable( Numbertable_type L , int number );

#endif

