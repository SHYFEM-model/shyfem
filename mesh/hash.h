
/* $Id: hash.h,v 1.3 2009-04-21 10:25:35 georg Exp $ */

/************************************************************************\ 
 *									*
 * hash.h - hash table administration routines                          *
 *                                                                      *
 * Copyright (c) 1992-1994 by Georg Umgiesser                           *
 *                                                                      *
 * see hash.c for copying information                                   *
 *                                                                      *
 * Revision History:                                                    *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 06-Apr-94: copyright notice added to file                            *
 * ..-...-92: routines written from scratch                             *
 *                                                                      *
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
