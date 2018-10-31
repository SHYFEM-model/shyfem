
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
 * xmouse.c - mouse driver routines under X11				*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
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

#ifndef __GUC_XMOUSE_
#define __GUC_XMOUSE_

#include <X11/Xlib.h>
#include <X11/Xutil.h>

/*#include <dos.h>*/

#include "mouse.h"

extern Window *QActWindow( void );
extern Display *QActDisplay( void );

/**********************************************************************/

void	MouseInit( void )

/* initialization of mouse driver */

{
}

void	MouseShow( void )

/* show mouse cursor */

{
}

void	MouseHide( void )

/* hide mouse cursor */

{
}

void    MouseStatus( int *mbutton , int *mhoriz , int *mverti )

/* mouse status info */

/* mbutton = 1  -> left button */
/* mbutton = 2  -> right button */
/* mbutton = 4  -> middle button */
/* combinations are possible */

{
}

void    MouseMove( int mhoriz , int mverti )

/* move mouse curser */

{
	XWarpPointer(QActDisplay(),None,*QActWindow(),0,0,0,0,
			mhoriz,mverti);
}

void    MousePress( int mbutton , int *mpress , int *mhoriz , int *mverti )

/* mouse button press info */

{
}

void    MouseRelease( int mbutton , int *mpress , int *mhoriz , int *mverti )

/* mouse button release info */

{

}

void    MouseWindow( int xmin , int ymin , int xmax , int ymax )

/* define window for mouse cursor */

{
}

void    MouseHole( int xmin , int ymin , int xmax , int ymax )

/* define hide window for mouse cursor */

/* the hide window can be reset by calling MouseInit or MouseShow */

{
}

void    MouseInfo( int *version , int *mtyp , int *mirq )

/* info about mouse type */

{
}

#endif

