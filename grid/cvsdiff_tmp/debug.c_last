
/************************************************************************\ 
 *									*
 * debug.c - debugging routines						*
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
		Error("Cannot allocate debug string");
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
