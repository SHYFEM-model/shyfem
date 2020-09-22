
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1998  Georg Umgiesser
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
 * general.h - header for general routines and definitions
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 * 14.09.1995	ggu	ABS included in header
 * 13.02.1998	ggu	test automatically if unix or dos
 * 20.03.1998	ggu	MIN, MAX, ROUND included in header
 * 20.03.1998	ggu	ASSERT_DEBUG introduced
 *
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
