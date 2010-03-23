
/************************************************************************\ 
 *									*
 * mouse.c - mouse driver routines for DOS				*
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

/*

for mouse buttons use :

NO_MOUSE_BUTTON     0
LEFT_MOUSE_BUTTON   1
RIGHT_MOUSE_BUTTON  2
MIDDLE_MOUSE_BUTTON 4

*/

#ifndef __GUC_MOUSE_
#define __GUC_MOUSE_


#include <dos.h>
#include "mouse.h"


/**********************************************************************/

void	MouseInit( void )

/* initialization of mouse driver */

{
	union REGS reg;

	reg.x.ax = 0x0;

	int86( 0x33 , &reg , &reg );

	if( reg.x.ax == 0 )
		/* error */;
}

void	MouseShow( void )

/* show mouse cursor */

{
	union REGS reg;

	reg.x.ax = 0x1;

	int86( 0x33 , &reg , &reg );

}

void	MouseHide( void )

/* hide mouse cursor */

{
	union REGS reg;

	reg.x.ax = 0x2;

	int86( 0x33 , &reg , &reg );

}

void    MouseStatus( int *mbutton , int *mhoriz , int *mverti )

/* mouse status info */

/* mbutton = 1  -> left button */
/* mbutton = 2  -> right button */
/* mbutton = 4  -> middle button */
/* combinations are possible */

{
	union REGS reg;

    reg.x.ax = 0x3;

	int86( 0x33 , &reg , &reg );

    *mbutton = reg.x.bx;
    *mhoriz = reg.x.cx;
    *mverti = reg.x.dx;

}

void    MouseMove( int mhoriz , int mverti )

/* move mouse curser */

{
	union REGS reg;

    reg.x.ax = 0x4;
    reg.x.cx = mhoriz;
    reg.x.dx = mverti;

	int86( 0x33 , &reg , &reg );

}

void    MousePress( int mbutton , int *mpress , int *mhoriz , int *mverti )

/* mouse button press info */

{
    int mb;
	union REGS reg;

    mb = mbutton-1;
	if( mb < 0 )
        return;
	else if( mb == 3 )
		mb=mb-1;


	reg.x.ax = 0x5;
	reg.x.bx = mb;

	int86( 0x33 , &reg , &reg );

	*mpress = reg.x.bx;

	if( *mpress > 0 ) {
		*mhoriz = reg.x.cx;
		*mverti = reg.x.dx;
	}

}

void    MouseRelease( int mbutton , int *mpress , int *mhoriz , int *mverti )

/* mouse button release info */

{
    int mb;
	union REGS reg;

    mb = mbutton-1;
	if( mb < 0 )
        return;
	else if( mb == 3 )
        mb=mb-1;

    reg.x.ax = 0x6;
	reg.x.bx = mb;

	int86( 0x33 , &reg , &reg );

	*mpress = reg.x.bx;

	if( *mpress > 0 ) {
		*mhoriz = reg.x.cx;
		*mverti = reg.x.dx;
	}

}

void    MouseWindow( int xmin , int ymin , int xmax , int ymax )

/* define window for mouse cursor */

{
	union REGS reg;

    reg.x.ax = 0x7;
    reg.x.cx = xmin;
    reg.x.dx = xmax;

	int86( 0x33 , &reg , &reg );

    reg.x.ax = 0x8;
    reg.x.cx = ymin;
    reg.x.dx = ymax;

	int86( 0x33 , &reg , &reg );

}

void    MouseHole( int xmin , int ymin , int xmax , int ymax )

/* define hide window for mouse cursor */

/* the hide window can be reset by calling MouseInit or MouseShow */

{
	union REGS reg;

    reg.x.ax = 0x10;
    reg.x.cx = xmin;
    reg.x.dx = ymin;
    reg.x.si = xmax;
    reg.x.di = ymax;

	int86( 0x33 , &reg , &reg );

}

void    MouseInfo( int *version , int *mtyp , int *mirq )

/* info about mouse type */

{
	union REGS reg;

    reg.x.ax = 0x24;

	int86( 0x33 , &reg , &reg );

	*version = 100 * reg.h.bh + reg.h.bl;
	*mtyp = reg.h.ch;
	*mirq = reg.h.cl;

}

#endif
