
/************************************************************************\
 *									*
 * meshgd.h - read/write grd files                                      *
 *									*
 * Copyright (c) 1992-1995 by Georg Umgiesser				*
 *									*
 * see meshgd.c for copying information                                 *
 *                                                                      *
 * Revision History:                                                    *
 * 11-Aug-95: routines copied from gridfi.c                             *
 *									*
\************************************************************************/


#ifndef __GUH_MESHGD_
#define __GUH_MESHGD_


#include "hash.h"
#include "queue.h"


void ReadStandard( char *fname , Hashtable_type HN , Hashtable_type HE
			      , Hashtable_type HL , QueueTable C );
void WriteStandard( char *fname , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , QueueTable C );
int ReadNode( Hashtable_type H );
int ReadElem( Hashtable_type H );
int ReadLine( Hashtable_type H );


#endif
