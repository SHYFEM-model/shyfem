
/************************************************************************\ 
 *									*
 * nlist.c - NodeList utility routines					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
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
 * 08-Oct-97: handle 0 list correctly                                   *
 * 11-Aug-95: Nodelist renamed to NodeList                              *
 *            Nodelist_type * substituted by NodeList                   *
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "nlist.h"


NodeList MakeNodeList( int count )

{
	NodeList new;

	new = (NodeList) malloc( sizeof(NodeList_type) );
	if( !new )
		Error("MakeNodeList: Cannot allocate NodeList");
	if( count > 0 ) {
	    new->index = (int *) malloc( count*sizeof(int) );
	    if( !new->index ) {
		printf("%d\n",count);
		Error("MakeNodeList: Cannot allocate NodeList index");
	    }
	} else {
	    new->index = NULL;
	}
	new->count = count;

	return new;
}

void InvertNodeList( NodeList L )

{
	int i;
	int aux;

	for(i=0;i<L->count/2;i++) {
		aux = L->index[i];
		L->index[i] = L->index[L->count-1-i];
		L->index[L->count-1-i] = aux;
	}
}

void FreeNodeList( NodeList L )

{
	if( !L ) return;

	if( L->index )
		free(L->index);
	free(L);
}
