
/************************************************************************\ 
 *									*
 * general.h - header for general routines and definitions		*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see general.c for copying information				*
 *									*
 * Revision History:							*
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GENERAL_
#define __GUH_GENERAL_


#ifndef __GUG_DOS_
#define __GUG_DOS_     0
#endif

#ifndef __GUG_UNIX_
#define __GUG_UNIX_    1
#endif



#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef True
#define True  1
#endif
#ifndef False
#define False 0
#endif


void Error( char *s );
void Error2( char *s1 , char *s2 );
void Warning( char *s );
void Warning2( char *s1 , char *s2 );
char *GetTime( void );


#endif
