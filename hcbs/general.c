
/************************************************************************\ 
 *									*
 * general.c - general routines						*
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
 * 11-Feb-94: copyright notice added to all files			*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "general.h"

void Error( char *s )

{
        printf("%s\n",s);
        exit(1);
}

void Error2( char *s1 , char *s2 )

{
        printf("%s %s\n",s1,s2);
        exit(1);
}

void Warning( char *s )

{
        printf("%s\n",s);
}

void Warning2( char *s1 , char *s2 )

{
        printf("%s %s\n",s1,s2);
}

char *GetTime( void )

{
        time_t t;
        char *s;
        char *p;

        time( &t );
        s=asctime(localtime(&t));

        p=s;
        while( *s != '\n' )
                s++;
        *s='\0';

        return p;
}

