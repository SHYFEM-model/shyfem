
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
 * debug.c - debugging routines						*
 *									*
 * Revision History:							*
 * 31-Dec-94: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"
#include "debug.h"

#define NDEBUG 10000

static char saux[1000];
static char *sdebug=NULL;
static int slen=0;

void GDPrint( void )

{
	if(sdebug) {
		printf("%s",sdebug);
		printf("Debug string : %d\n",(int) strlen(sdebug));
	}
}

void GDInit( int size )

{
	if( !sdebug )
		sdebug=GDAlloc( size );
}

char *GDAlloc( int size )

{
	char *s;

	s = (char *) malloc( size );
	if( !s )
		Error("GDAlloc: Cannot allocate debug string");
	return s;
}

void GDNL( void )

{
	GDS("\n");
}

void GDS( char *s )

{
	int i;

	if( !sdebug ) GDInit(NDEBUG);

	i=strlen(s);
	if( slen+i < NDEBUG ) {
		strcat(sdebug,s);
		slen+=i;
	}
/*	printf("%s",s); */
}

void GDI( int i )

{
	sprintf(saux,"%d ",i);
	GDS(saux);
}

void GDL( long l )

{
	sprintf(saux,"%ld ",l);
	GDS(saux);
}

void GDF( float f )

{
	sprintf(saux,"%f ",f);
	GDS(saux);
}
