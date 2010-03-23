
/************************************************************************\ 
 *									*
 * nlist.h - NodeList utility routines					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see nlist.c for copying information					*
 *									*
 * Revision History:							*
 * 11-Aug-95: Nodelist renamed to NodeList                              *
 *            Nodelist_type * substituted by NodeList                   *
 * 27-Jul-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_NLIST_
#define __GUH_NLIST_


typedef struct nodelist_tag {
        int count;
        int *index;
} NodeList_type;

typedef NodeList_type *NodeList;

NodeList MakeNodeList( int count );
void InvertNodeList( NodeList L );
void FreeNodeList( NodeList L );


#endif

