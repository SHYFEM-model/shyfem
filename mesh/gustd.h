
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
 * gustd.h - header for standard routines				*
 *									*
 * Revision History:							*
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-88: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GUSTD_
#define __GUH_GUSTD_

#include "general.h"
#include <stdio.h>

#define YES      1
#define NO       0

#define         swapi(x,y)      {int z; z=(x); (x)=(y); (y)=z;}
#define         swapip(x,y)     {int *z; z=(x); (x)=(y); (y)=z;}
#define         swapcp(x,y)     {char *z; z=(x); (x)=(y); (y)=z;}

#define         priarray(a,n)   {int i; for(i=0;i<n;i++)   \
				printf("%6d%c",a[i],(i%10==9 || i==n-1) ? '\n' : ' ');}

#if __GUG_DOS_
int     gindex(char *s ,char *t, int n);
			/* find n th occurence of t in s */
int     grindex(char *s ,char *t, int n);
			/* find n th occ. of t in s (reverse) */
#endif

int     stripline (char *s);            /* cancel \n from end of line */
int     fgetline ( FILE *fp , char *s , int lim );
char *getlin(FILE *fp);
FILE *ifileq(char *text , char *mode);
FILE *filop(char *string , char *mode);
FILE *filopn(char *file , char *mode);
char *nextword(char *s);
char *firstword(char *s);
int countword(char *s);
int strsize(char *s);
char *saveword(char *s);
char *isolword(char *s);
char *savestring(char *s , int len);
void squeeze(char *s , int c);
double power(double d , int p);
char *itos( int i );
void reverse( char *s );
char *itoa10( int n , char *s );

#endif

