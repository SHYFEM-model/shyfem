
/************************************************************************\ 
 *									*
 * gustd.c - standard utilities						*
 *									*
 * Copyright (c) 1988-1994 by Georg Umgiesser				*
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
 * 21-Mar-94: gcc-warnings (minor changes)				*
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-88: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUC_GUSTD_
#define __GUC_GUSTD_


#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "gustd.h"

/******************************************************************/

#if __GUG_DOS_

int             index(char *s ,char *t, int n)

/* find n th occurence of t in s */

	{
		int i, j, k;

		for (i = 0; s[i] != '\0'; i++) {
			for (j=i, k=0; t[k] != '\0' && s[j] == t[k]; j++, k++)
				;
			if (t[k] == '\0' && --n == 0)
				return(i);
		}
		return(-1);
	}

/**************************************************************************/

int             rindex(char *s ,char *t, int n)

/* find n th occurence from the right of t in s */

	{
		int i, j, k;

		for (i = strlen(s) - 1; i >= 0 ; i--) {
			for (j=i, k=0; t[k] != '\0' && s[j] == t[k]; j++, k++)
				;
			if (t[k] == '\0' && --n == 0)
				return(i);
		}
		return(-1);
	}

#endif

/**************************************************************************/

int     stripline (char *s)             /* cancel \n from end of line */

	{
		int     i = 0;

		while ( s[i] )
			if( s[i] == '\n')
				s[i] = '\0';
			else
				i++;
		return(i);
	}

/**************************************************************************/

int fgetline ( FILE *fp , char *s , int lim )

{
	int i = 0;
	int c = 0;

	while ( i<lim-1 && (c=getc(fp))!=EOF && c!='\n' )
		s[i++] = c;
	if ( c == '\n' )
		s[i++] = c;
	s[i] = '\0';
	return(i);
}

/**************************************************************************/

#define NMAX    256

char *getlin(FILE *fp)

{
	int i = 0;
	int c = 0;
	static int nmax = NMAX;
	static char line[NMAX] = {'\0'};

	if(fp==NULL) return(NULL); /* do not read anything, return NULL */

	while(i<nmax && (c=getc(fp)) != EOF && c != '\n')
		line[i++]=c;
	if(i==nmax) {
		printf("Warning : string to small in GETLIN\n");
		i--;
	} else if( c == EOF && i == 0) {
		return(NULL);
	}
	line[i++]='\0';

	return(line);
}

#undef NMAX

/*******************************************************************/

FILE *ifileq( char *text, char* mode )

{
	char* name;

	printf("%s ",text);
	name=getlin(stdin);
	return(fopen(name,mode));
}

/*******************************************************************/

FILE *filop(char *string , char *mode)

{
	char *line;
	FILE *fp;

	printf("%s",string);
	line=getlin(stdin);
	fp=fopen(line,mode);
	if(fp==NULL)
		printf("Cannot open file : %s\n",line);
	return(fp);
}

/*******************************************************************/

FILE *filopn(char *file , char *mode)

{
	FILE *fp;

	fp=fopen(file,mode);
	if(fp==NULL)
		printf("Cannot open file : %s\n",file);
	return(fp);
}

/***********************************************************************/

char *nextword(char *s)

{
	while( !isspace((int)(*s)) && ((*s) != '\0') ) /* if word skip chars */
		s++;
	while(  isspace((int)(*s)) )    /* now skip blanks */
		s++;
	return(s);
}

/**************************************************************************/

char *firstword(char *s)

{
	while( isspace((int)(*s)) )     /* skip blanks */
		s++;
	return(s);
}

/**************************************************************************/

int countword(char *s)

{
	int i=0;

	if( *(s=firstword(s)) ) i++;
	while( *(s=nextword(s)) ) i++;
	return(i);
}

/**************************************************************************/

int strsize(char *s)

{
	int i;

	i=strlen(s)-1;
	while( isspace((int)s[i]) )
		i--;
	return(i+1);
}

/**************************************************************************/

char *saveword(char *s)

{
	int i;
	char *p;

	i=0;                                            /* i is number of chars in string */
	while( !isspace((int)s[i]) && (s[i] != '\0') )
		i++;

	p = (char *) malloc( i+1 );             /* allocate memory */

	if(p == NULL) {                 /* enough space ? */
		printf("No memory left to allocate word");
		exit(3);
	}

	p = strncpy(p,s,i);                             /* copy string to allocated space */
	p[i] = '\0';

	return(p);
}

/**************************************************************************/

#define NMAX 100

char *isolword(char *s)

{
	int i;
	static char p[NMAX];

	i=0;                                            /* i is number of chars in string */
	while( !isspace((int)s[i]) && (s[i] != '\0') )
		i++;
	if(i>=NMAX) {
		printf("Warning : string to small in ISOLWORD\n");
		i=NMAX-1;
	}
	strncpy(p,s,i);                         /* copy string to allocated space */
	p[i] = '\0';

	return(p);
}

#undef NMAX

/**************************************************************************/

char *savestring(char *s , int len)

{
	char *p;

	if(len < 0) len = strlen(s);    /* len is length to be copied */

	p = (char *) malloc( len+1 );   /* allocate memory */

	if(p == NULL) {                 /* enough space ? */
		printf("No memory left to allocate string");
		exit(3);
	}

	p = strncpy(p,s,len);           /* copy string to allocated space */
	p[len] = '\0';

	return(p);
}

/**************************************************************************/

void squeeze(char *s , int c)

/* deletes all c from s */

{
	char *t;

	t = s;
	while (*t != '\0')
		if(*t != c)
			*(s++) = *(t++);
	*s = '\0';
}

/********************************************************************/

double power( double d , int p )

{
	if(p)
		return(d*power(d,--p));
	else
		return(1.);
}

/********************************************************************/

char *itos( int i )

{
	static char s[20];

	return itoa10(i,s);
}

/********************************************************************/

void reverse( char *s )

{
	int c,i;
	int j=strlen(s)-1;

	for(i=0;i<j;i++,j--) {
		c=s[i];
		s[i]=s[j];
		s[j]=c;
	}
}

/********************************************************************/

char *itoa10( int n , char *s )

{
	int i=0;
	int segno=n;

	if( n < 0 )
		n = -n;

	do {
		s[i++] = n % 10 + '0';
	} while( (n /= 10) > 0 );

	if( segno < 0 )
		s[i++] = '-';

	s[i] = '\0';
	reverse(s);

	return s;
}







#endif


