
/************************************************************************\ 
 *									*
 * xmouse.c - mouse driver routines under X11				*
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

