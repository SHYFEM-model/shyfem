
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
 * hash.h - hash table administration routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	completely restructured -> independent routines
 *
\************************************************************************/


#ifndef __GUH_HASH_
#define __GUH_HASH_


typedef union {
	int ii;
	char* cp;
	void* vp;
} Key_type;

typedef struct list_hash_tag {
	Key_type key;
        void *info;
        struct list_hash_tag *next;
} List_hash_type;

typedef struct {
	List_hash_type **head;
	List_hash_type  *pact;
	int              hact;
} Hashtable_head;


typedef List_hash_type  *Hashtable_item;
typedef Hashtable_head  *Hashtable_type;


Hashtable_type MakeHashTable( void );
void ResetHashTable( Hashtable_type H );
void *VisitHashTable( Hashtable_type H );
void FreeHashTable( Hashtable_type H );
void InsertHashByNumber( Hashtable_type H , void *info , int number );
void *RetrieveHashByNumber( Hashtable_type H , int number );
void *DeleteHashByNumber( Hashtable_type H , int number );
void *CheckHashByNumber( Hashtable_type H );


#endif
