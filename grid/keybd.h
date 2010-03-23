
/************************************************************************\ 
 *									*
 * keybd.h - keyboard routines						*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see keybd.c for copying information					*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_KEYBD_
#define __GUH_KEYBD_


/**************************************************************************/

int	strike_key(void);
int	wait_for_key(void);
void	press_any_key(void);

/**************************************************************************/

#endif
