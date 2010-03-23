
/************************************************************************\ 
 *									*
 * args.h - utilities to read arguments from line                       *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see args.c for copying information					*
 *									*
 * Revision History:							*
 * 08-May-94: created from former routines                              *
 *									*
\************************************************************************/

#ifndef __GUH_ARGS_
#define __GUH_ARGS_

char *firstchar( char *s );
void initargs( char *s );
int nargs( char *s );
char *readargs( void );

#endif
