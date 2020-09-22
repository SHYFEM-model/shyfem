
/************************************************************************\
 *
 *    Copyright (C) 1994  Georg Umgiesser
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
 * tevents.c - event routines for Turbo C
 *
 * revision log :
 *
 * 04.12.1994	ggu	routines written from scratch
 *
\************************************************************************/



/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
*/

/*
#include "graph.h"
#include "xgraph.h"

#include "generalx.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
*/

#include <stdio.h>

#include "general.h"
#include "mouse.h"
#include "keybd.h"
#include "events.h"

static QEventMask DefaultEventMask	= QKeyPressMask
					| QButtonPressMask
					| QExposureMask
					;
static QEventMask ActualEventMask;

/*
static QEvent     ActualEvent;
*/


void QInitEvent( void )


{
	QSelectEvent( DefaultEventMask );
	MouseInit();
	MouseShow();
}


QEventMask QSelectedEvent( void )

{
	return ActualEventMask;
}

void QSelectEvent( QEventMask eventmask )

{
	ActualEventMask = eventmask;
}

void QAddEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask | eventmask;
}

void QDeleteEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask & ~eventmask;
}

void QNextEvent( QEvent *eventp )

{
	int c;
	int button,x,y;
	static int oldbutton = 0;
	static int oldx = 0;
	static int oldy = 0;
	static int first = TRUE;
	int dbutton,dxy;

	eventp->type = QNullEvent;

	c=strike_key();

	MouseStatus( &button , &x , &y );

	dbutton = button ^ oldbutton;
	dxy = ( x ^ oldx ) | ( y ^ oldy );

	/* we test  1) button press  2) keyboard press  3) pointer move */

	if( dbutton && ( ActualEventMask & QButtonPressMask ) ) {

		oldbutton = button;

		eventp->type = QButtonPress;
		eventp->button.x = x;
		eventp->button.y = y;

		if( dbutton & LEFT_MOUSE_BUTTON ) {

			eventp->button.button = QButtonLeft;
			if( button & LEFT_MOUSE_BUTTON )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		} else if( dbutton & RIGHT_MOUSE_BUTTON ) {

			eventp->button.button = QButtonRight;
			if( button & RIGHT_MOUSE_BUTTON )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		} else {

			eventp->button.button = QButtonMiddle;
			if( button & MIDDLE_MOUSE_BUTTON )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		}

		return;
	}

	if( c && ( ActualEventMask & QKeyPressMask ) ) {

		eventp->type = QKeyPress;
		eventp->key  = c;

		return;
	}

	if( dxy && ( ActualEventMask & QPointerMoveMask ) ) {

		oldx = x;
		oldy = y;

		eventp->type = QPointerMove;
		eventp->button.x = x;
		eventp->button.y = y;

		return;
	}

	if( first ) {

		first = FALSE;
		eventp->type = QExpose;

		return;
	}
}
