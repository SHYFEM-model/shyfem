
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
 * keybd.c - keyboard routines						*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUC_KEYBD_
#define __GUC_KEYBD_

#include <stdio.h>

#include "general.h"
#include "keybd.h"


#if __GUG_DOS_

#include <conio.h>

int     strike_key(void)

{
	int     key=0;

	if(kbhit()) {
		key=getch();
		if(!key) key = -getch();
	}
	return(key);
}

int     wait_for_key(void)

{
	int     key;

	while( (key=strike_key()) == 0 ) ;
	return(key);
}

void    press_any_key( void )

{
	printf("Press any key to continue...");
	getch();
	printf("\n");
}

#else

int     strike_key(void)

{
	return getchar();
}

int     wait_for_key(void)

{
	return getchar();
}

void    press_any_key( void )

{
	printf("Press <CR> to continue...");
	getchar();
/*        printf("\n");*/
}

#endif

#endif

