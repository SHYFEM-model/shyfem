
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * meshgd.h - read/write grd files
 *
 * revision log :
 *
 * 11.08.1995	ggu	routines copied from gridfi.c
 *
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
