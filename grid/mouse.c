
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
 * mouse.c - mouse driver routines for DOS
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 *
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
