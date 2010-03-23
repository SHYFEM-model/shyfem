
/************************************************************************\
 *                                                                      *
 * gridky.h - keyboard routines                                         *
 *                                                                      *
 * Copyright (c) 1992-1994 by Georg Umgiesser                           *
 *                                                                      *
 * see gridky.c for copying information                                 *
 *                                                                      *
 * Revision History:                                                    *
 * 19-Nov-97: new routines SetKeyboardString(), AppendKeyboardString()  *
 * 04-Feb-95: routines written from scratch                             *
 *                                                                      *
\************************************************************************/


#ifndef __GUH_GRIDKY_
#define __GUH_GRIDKY_


void ResetKeyboardInput( void );
void AddKeyboardInput( int c );
void SubKeyboardInput( void );
char *GetKeyboardString( void );

void SetKeyboardString( char *s );
void AppendKeyboardString( char *s );

int GetKeyboardInt( void );
float GetKeyboardFloat( void );
	

#endif 
