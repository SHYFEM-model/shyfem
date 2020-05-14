
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
 * gridky.h - keyboard routines
 *
 * revision log :
 *
 * 04.02.1995	ggu	routines written from scratch
 * 19.11.1997	ggu	new routines SetKeyboardString(), AppendKeyboardString()
 *
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
