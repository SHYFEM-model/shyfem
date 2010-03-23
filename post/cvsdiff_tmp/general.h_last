
/* $Id: general.h,v 1.2 1998/03/20 09:22:58 georg Exp $ */

/************************************************************************\ 
 *									*
 * general.h - header for general routines and definitions		*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see general.c for copying information				*
 *									*
 * Revision History:							*
 * 20-Mar-1998: ASSERT_DEBUG introduced                                 *
 * 20-Mar-1998: MIN, MAX, ROUND included in header                      *
 * 13-Feb-1998: test automatically if unix or dos                       *
 * 14-Sep-95: ABS included in header                                    *
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GENERAL_
#define __GUH_GENERAL_


/*****************************************/
#define ASSERT_DEBUG    1
/*****************************************/


#ifdef  __unix__
#define __GUG_UNIX_    1
#else
#define __GUG_UNIX_    0
#endif

#if __GUG_UNIX_
#define __GUG_DOS_     0
#else
#define __GUG_DOS_     1
#endif

#ifndef YES
#define YES   1
#endif
#ifndef NO
#define NO    0
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


#define ABS(x)		( ( (x) > 0 ) ? (x) : (-(x)) )
#define MAX(a,b)        ( (a) > (b) ? (a) : (b) )
#define MIN(a,b)        ( (a) < (b) ? (a) : (b) )
#define ROUND(x)	( ((x) > 0.) ? ((int)((x)+0.5)) : ((int)((x)-0.5)) )


#endif
