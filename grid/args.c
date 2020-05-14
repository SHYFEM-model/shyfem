
/************************************************************************\
 *
 *    Copyright (C) 1994-1995  Georg Umgiesser
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
 * args.c - utilities to read arguments from line
 *
 * revision log :
 *
 * 08.05.1994	ggu	created from former routines
 * 15.03.1995	ggu	bug fix in readargs (trailing blanks returned as arg)
 *
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
