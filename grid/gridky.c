
/************************************************************************\
 *
 *    Copyright (C) 1995,1997  Georg Umgiesser
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
 * gridky.c - routines for keyboard input
 *
 * revision log :
 *
 * 04.02.1995	ggu	routines written from scratch
 * 19.11.1997	ggu	new routines SetKeyboardString(), AppendKeyboardString()
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid.h"
#include "graph.h"

#include "general.h"
#include "gustd.h"


#define KEYBOARD_DIM		60
#define KEYBOARD_DIM1		( KEYBOARD_DIM - 1 )

static int   KeyboardPointer=0;
static char  KeyboardString[KEYBOARD_DIM];


void ResetKeyboardInput( void )

{
	KeyboardPointer=0;
/*	KeyboardString[KeyboardPointer]='\0';	*/
}

void AddKeyboardInput( int c )

{
	if( KeyboardPointer < KEYBOARD_DIM1 ) {
		KeyboardString[KeyboardPointer++]=c;
		KeyboardString[KeyboardPointer]='\0';
	}
}

void SubKeyboardInput( void )

{
	if( KeyboardPointer ) {
		KeyboardString[--KeyboardPointer]='\0';
	}
}

char *GetKeyboardString( void )

{
	return KeyboardString;
}
	
void SetKeyboardString( char *s )

{
	strncpy(KeyboardString,s,KEYBOARD_DIM);
	KeyboardString[KEYBOARD_DIM1] = '\0';
	KeyboardPointer = strlen(KeyboardString);
}
	
void AppendKeyboardString( char *s )

{
	char *t = &(KeyboardString[KeyboardPointer]);
	strncpy(t,s,KEYBOARD_DIM - KeyboardPointer);
	KeyboardString[KEYBOARD_DIM1] = '\0';
	KeyboardPointer = strlen(KeyboardString);
}
	
int GetKeyboardInt( void )

{
	return atoi(KeyboardString);
}
	
float GetKeyboardFloat( void )

{
	return atof(KeyboardString);
}
	
