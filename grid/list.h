
/************************************************************************\
 *
 *    Copyright (C) 1992,1994  Georg Umgiesser
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
 * list.h - list table administration routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	completely restructured -> independent routines
 * 17.04.1994	ggu	NumberTable types added
 *
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

