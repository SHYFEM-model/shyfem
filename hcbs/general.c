
/************************************************************************\
 *
 *    Copyright (C) 1992,1994  Georg Umgiesser
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
 * general.c - general routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 *
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

