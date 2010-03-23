
/************************************************************************\ 
 *									*
 * keybd.c - keyboard routines						*
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

