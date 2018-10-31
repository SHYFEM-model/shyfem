
/************************************************************************\
 *
 *    Copyright (C) 1985-2018  Georg Umgiesser
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
 *									*
 * nlist.c - NodeList utility routines					*
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
