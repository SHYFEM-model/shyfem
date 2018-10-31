
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
 * general.h - header for general routines and definitions		*
 *									*
 * Revision History:							*
 * 01-Jun-2012: HACK for grid to work on Mac: __GUG_UNIX_ = 1		*
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


#ifndef __GUG_UNIX_
#ifdef  __unix__
#define __GUG_UNIX_    1
#else
#define __GUG_UNIX_    0
#endif
#endif
#define __GUG_UNIX_    1	/* HACK for Mac */

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
