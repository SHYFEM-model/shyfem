
/************************************************************************\ 
 *									*
 * gustd.h - header for standard routines				*
 *									*
 * Copyright (c) 1988-1994 by Georg Umgiesser				*
 *									*
 * see gustd.c for copying information					*
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
int     index(char *s ,char *t, int n);
			/* find n th occurence of t in s */
int     rindex(char *s ,char *t, int n);
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

