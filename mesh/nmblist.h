
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995  Georg Umgiesser
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
 * nmblist.h - number table administration routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	completely restructured -> independent routines
 * 17.04.1994	ggu	NumberTable types added
 * 11.08.1995	ggu	NumberTable routines transfered to this file
 * ...		ggu	Numbertable_type renamed in NumberTable
 *
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

