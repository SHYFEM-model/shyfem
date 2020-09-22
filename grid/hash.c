
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1997,2003  Georg Umgiesser
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
 * hash.c - hash table administration routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	completely restructured -> independent routines
 * 11.03.1995	ggu	some comments added
 * 07.12.1995	ggu	MakeHashPointer, VisitHashTableG now defined static (bug)
 * 11.10.1997	ggu	FreeHashTable: bug - using VisitHashTableG used deleted
 * ...		ggu	pointer and crashed Linux -> use for-loop
 * 16.10.2003	ggu	increase hash table size HASHSIZE
 *
\************************************************************************/


#include <stdlib.h>

#include "general.h"
#include "hash.h"


/*
#define HASHSIZE 997
*/

#define HASHSIZE 50379

#define HASHBYNUMBER(h) ((h) % HASHSIZE )

static List_hash_type *MakeHashPointer( void );
static List_hash_type *VisitHashTableG( Hashtable_type H );

/*
static List_hash_type *pact=NULL;
static int hact=0;
*/


/***************************************************************/


Hashtable_type MakeHashTable( void )

{
	Hashtable_item *newitem;
	Hashtable_head *newhead;
	int i;

	newitem = (Hashtable_item *) malloc( sizeof(Hashtable_item)*HASHSIZE );
	if( !newitem ) Error("MakeHashTable : Cannot allocate Hashtable");

	for(i=0;i<HASHSIZE;i++)
		newitem[i]=NULL;

	newhead = (Hashtable_head *) malloc( sizeof(Hashtable_head) );
	if( !newhead ) Error("MakeHashTable : Cannot allocate Hashtable");

	newhead->head = newitem;
	newhead->pact = NULL;
	newhead->hact = 0;

        return newhead;
}

static List_hash_type *MakeHashPointer( void )

{
        List_hash_type *new;

        new = (List_hash_type *) malloc( sizeof( List_hash_type ) );

        if( !new ) Error("MakeHashPointer : Cannot allocate node");

        new->next = NULL;
        new->info = NULL;
        new->key.vp  = NULL;

        return new;
}

void ResetHashTable( Hashtable_type H )

{
	H->pact = NULL;
	H->hact = 0;
}

static List_hash_type *VisitHashTableG( Hashtable_type H )

{
        if( !H->pact ) {
            H->pact = H->head[0];
            H->hact = 0;
        } else {
            H->pact = H->pact->next;
        }

        while( !H->pact ) {
	    (H->hact)++;
            if( H->hact == HASHSIZE ) return NULL;
	    H->pact = H->head[H->hact];
	}

	return H->pact;
}

void *VisitHashTable( Hashtable_type H )

{
	List_hash_type *p;

	p = VisitHashTableG( H );

	if( p )
		return p->info;
	else
		return NULL;
}

void FreeHashTable( Hashtable_type H )

/*\
 *  info structure must have been deleted prior
 *  to calling this function - the info pointers
 *  are left dangling otherwise
\*/

{
	int i;
        List_hash_type *p,*pnext;

	for(i=0;i<HASHSIZE;i++) {
	  p = H->head[i];
	  while( p ) {
	    pnext = p->next;
	    free(p);
	    p = pnext;
	  }
	} 

	free(H->head);
        free(H);
}

void InsertHashByNumber( Hashtable_type H , void *info , int number )

{
        int h;
        List_hash_type *p;

	h = HASHBYNUMBER( number );
        p = MakeHashPointer();
	p->info = (void *) info;
	p->key.ii  = number;
        p->next = H->head[h];
        H->head[h] = p;
}

void *RetrieveHashByNumber( Hashtable_type H , int number )

{
        List_hash_type *p;

	p=H->head[HASHBYNUMBER(number)];

	while( p!=NULL ) {
		if( p->key.ii == number )
			return p->info;
        	p=p->next;
        }

	return NULL;
}

void *DeleteHashByNumber( Hashtable_type H , int number )

/*\
 *  deletes the hash entry but not the info-structure 
 *  -> this can be deleted with the pointer passed back
\*/

{
        int h;
	void *v;
	List_hash_type *p,*r;

	h=HASHBYNUMBER(number);
	p=H->head[h];
	r=NULL;	/* one step behind */

	while( p!=NULL ) {
		if( p->key.ii == number ) {
			if( r )
				r->next=p->next;
			else
				H->head[h]=p->next;
			v = p->info;
			free(p);
			return v;
		}
		r=p;
		p=p->next;
	}
	return NULL;
}

void *CheckHashByNumber( Hashtable_type H )

/*\
 *  befor calling must call ResetHashTable
 *  returns second entry, that is doubled, the first can be
 *  obtained by a simple RetrieveHashByNumber() call
\*/

{
	List_hash_type *p;

		if( !H->pact ) {
            H->pact = H->head[0];
            H->hact = 0;
        } else {
            H->pact = H->pact->next;
        }

	for( ; ; ) {
	    while( H->pact ) {
		p=H->pact->next;
		while( p ) {
			if( p->key.ii == H->pact->key.ii )
				return p->info;
			p=p->next;
		}
		H->pact=H->pact->next;
	    }
		if( ++(H->hact) == HASHSIZE ) break;
		H->pact = H->head[H->hact];
	}
	return NULL;
}

