
/* $Id: hash.c,v 1.3 2004/09/24 11:58:04 georg Exp $ */

/************************************************************************\ 
 *									*
 * hash.c - hash table administration routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
 *									*
 * Revision History:							*
 * 16-Oct-2003: increase hash table size HASHSIZE                       *
 * 11-Oct-97: FreeHashTable: bug - using VisitHashTableG used deleted   *
 *              pointer and crashed Linux -> use for-loop               *
 * 07-Dec-95: MakeHashPointer, VisitHashTableG now defined static (bug) *
 * 11-Mar-95: some comments added                                       *
 * 13-Apr-94: completely restructured -> independent routines           *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
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

