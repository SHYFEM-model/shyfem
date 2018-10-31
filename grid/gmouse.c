
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
 * 13-Feb-1998: adjourned to GRX2.0                                     *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#define GRX_VERSION_GGU	2

#include "mouse.h"

#if GRX_VERSION_GGU == 1
#include <mousex.h>
#else
#include <grx20.h>
#endif

#if GRX_VERSION_GGU == 1
#define	GrMouseWarp(h,v)	MouseWarp( (h) , (v) );
#endif

/**********************************************************************/
/*
void	MouseInit( void )
*/
/* initialization of mouse driver */
/*
{
}
*/

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
	GrMouseWarp(mhoriz,mverti);
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

