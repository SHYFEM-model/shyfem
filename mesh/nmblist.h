
/************************************************************************\ 
 *									*
 * nmblist.h - number table administration routines                     *
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see nmblist.c for copying information				*
 *									*
 * Revision History:							*
 * 11-Aug-95: NumberTable routines transfered to this file              *
 *            Numbertable_type renamed in NumberTable                   *
 * 17-Apr-94: NumberTable types added                                   *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_NMBLIST_
#define __GUH_NMBLIST_


typedef struct number_tag {
        int number;
        struct number_tag *next;
} Number_type;

typedef struct {
        Number_type *head;
        Number_type *pact;
} NumberTable_head;

typedef NumberTable_head *NumberTable;

NumberTable MakeNumberTable( void );
void ResetNumberTable( NumberTable L );
void FreeNumberTable( NumberTable L );
int FindNumberTable( NumberTable L , int number );
void DeleteNumberTable( NumberTable L , int number );
int OkNumberTable( NumberTable L );
int VisitNumberTable( NumberTable L );
void InsertNumberTable( NumberTable L , int number );

#endif

