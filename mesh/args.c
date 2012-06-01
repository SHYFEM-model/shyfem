
/************************************************************************\
 *									*
 * args.c - utilities to read arguments from line                       *
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
 * 15-Mar-95: bug fix in readargs (trailing blanks returned as arg)     *
 * 08-May-94: created from former routines                              *
 *									*
\************************************************************************/

/*\
 *  Routines read a string of blank separated "args".
 *  This can be numbers or anything else to be processed.
 *  The string must be prepared by a call to nargs(s)
 *  or initargs(s). Both initialize the string s, but
 *  nargs(s) returns also the number of arguments in the string.
 *  The arguments can be read by further calls to readargs().
 *  If no arguments are left readargs() returns NULL.
 *  The string s will be altered by calls to these routines.
\*/


#include <stdio.h>

#include "general.h"


static int FirstArgs;		/* first call for string */
static char *StringArgs;	/* string with "arguments" */

char *firstchar( char *s ) /* finds first char that is not blank or tab */

{
	while( *s ) {
		if( *s != ' ' && *s != '\t' ) 
			break;
		s++;
	}
	return s;
}

void initargs( char *s )

{
	FirstArgs = TRUE;
	StringArgs = s;

	while( *s ) {
		if( *s == ',' || *s == '\t' )
			*s = ' ';
		s++;
	}
}

int nargs( char *s )

{
	int inword=FALSE;
	int n=0;

	FirstArgs = TRUE;
	StringArgs = s;

	while( *s ) {
		if( *s == ',' || *s == '\t' )
		{
			*s = ' ';
		}
		if( *s == ' ' )
		{
			inword = FALSE;
		}
		else	/* a character */
		{
			if( !inword ) n++;
			inword = TRUE;
		}
		s++;
	}
	return n;
}

char *readargs( void )

{
	static int last=FALSE;
	static char *cstart=NULL;
	static char *cend=NULL;
	char *s;

	if( FirstArgs ) {
		FirstArgs = FALSE;
		last = FALSE;
		s = StringArgs;
	} else if( last ) {
		return NULL;
	} else {
		s = ++cend;
	}

	while( *s == ' ' ) /* skip initial blanks */
		s++;

	if( *s == '\0' ) {	 /* look forEOF */
		last = TRUE;	/* bug fixed 15.3.95 */
		return NULL;
	}

	cstart = s;

	while( *s != '\0' && *s != ' ' )
		s++;

	if( *s == '\0' ) { /* end of line found, remember for next call */
		last = TRUE;
	} else {
		cend = s;
		*s = '\0';
	}

	return cstart;
}
