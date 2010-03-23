
/************************************************************************\ 
 *									*
 * tevents.c - event routines for Turbo C				*
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
 * 04-Dec-94: routines written from scratch				*
 *									*
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
