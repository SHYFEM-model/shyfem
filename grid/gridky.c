
/************************************************************************\ 
 *									*
 * gridky.c - routines for keyboard input                               *
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
 * 19-Nov-97: new routines SetKeyboardString(), AppendKeyboardString()  *
 * 04-Feb-95: routines written from scratch				*
 *									*
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
	
