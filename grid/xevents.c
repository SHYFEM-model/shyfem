
/************************************************************************\ 
 *									*
 * gevents.c - event routines for gcc under DOS				*
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
 * 05-Mar-2014: handle mouse wheel, new QSkipMotion()			*
 * 05-Dec-95: handle special keys from keyboard (arrows, return...)     *
 * 10-Dec-94: routines adapted to gcc through Event Queue Library       *
 * 04-Dec-94: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>

#include "generalx.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

#include "general.h"
#include "xgraph.h"
#include "events.h"

#define X_BUTTON_LEFT		( 1 )
#define X_BUTTON_RIGHT		( 3 )
#define X_BUTTON_MIDDLE		( 2 )
#define X_BUTTON_WHEEL_UP	( 4 )
#define X_BUTTON_WHEEL_DOWN	( 5 )

#define X_BUTTON_CHANGE_MASK	( ButtonPressMask | ButtonReleaseMask )
#define X_KEY_PRESS_MASK	( KeyPressMask )
#define X_MOTION_MASK		( PointerMotionMask )
#define X_CONFIGURE_MASK	( ExposureMask | StructureNotifyMask )
#define X_INPUT_MASK		( X_BUTTON_CHANGE_MASK | X_KEY_PRESS_MASK )
#define X_STANDARD_MASK		( X_INPUT_MASK | X_CONFIGURE_MASK ) 
#define X_ALL_MASK		( X_STANDARD_MASK | X_MOTION_MASK )

#define X_BUTTON_CHANGE		( ButtonPress | ButtonRelease )
#define X_KEY_PRESS		( KeyPress )
#define X_MOTION		( MotionNotify )
#define X_CONFIGURE		( Expose | ConfigureNotify )
#define X_MAPPING		( MappingNotify )
#define X_INPUT			( X_BUTTON_CHANGE | X_KEY_PRESS )

static QEventMask DefaultEventMask	= QKeyPressMask
					| QButtonPressMask
					| QExposureMask
					;
static QEventMask ActualEventMask;

static XEvent	 XActualEvent;
static Display	*MyDisplay;
static Window	 MyWindow;

static long	StandardMask;
static long	AllMask;
static long	MoveMask;

static int useless = 0;

#define GGU_DEBUG 1
/*#undef GGU_DEBUG*/
#ifdef GGU_DEBUG
#include "debug.h"
#endif

/*********************************************************************/

void QInitEvent( void )

{
	MyDisplay=QActDisplay();
	if( !MyDisplay )
		Error("QInitEvent : No display opened");
	MyWindow = *QActWindow();

	StandardMask =  ButtonPressMask | 
			ButtonReleaseMask |
			KeyPressMask |
			ExposureMask |
			StructureNotifyMask ;
	AllMask      =	StandardMask |
			PointerMotionMask ;
	MoveMask     =  ButtonPressMask |
			ButtonReleaseMask |
			PointerMotionMask ;

	QSelectEvent( DefaultEventMask );
}

void QDeInitEvent( void )

{
}

QEventMask QSelectedEvent( void )

{
	return ActualEventMask;
}

void QSelectEvent( QEventMask eventmask )

{
	ActualEventMask = eventmask;
	if( eventmask & QPointerMoveMask )
	    XSelectInput(MyDisplay,MyWindow,AllMask);
	else
	    XSelectInput(MyDisplay,MyWindow,StandardMask);
}

void QAddEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask | eventmask;
	if( eventmask & QPointerMoveMask ) {
		XSelectInput(MyDisplay,MyWindow,AllMask);
		/*
		  XSelectInput does not select Event during grab
		   -> make sure that PointerMove is reported also
		   between grabs using XChangeActivePointerGrab
		*/
		XChangeActivePointerGrab(MyDisplay,MoveMask,None,CurrentTime);
	}
}

void QDeleteEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask & ~eventmask;
	if( eventmask & QPointerMoveMask ) {
		XSelectInput(MyDisplay,MyWindow,StandardMask);
	}
}

void QSkipMotion( void )

/* skips all motion events in event queue */

{
	XEvent event;

	while ( XPending(MyDisplay) )
	{
	  XNextEvent( MyDisplay , &event );
	  if( event.type != MotionNotify ) {
	    XPutBackEvent( MyDisplay , &event );
	    return;
	  }
	}
}


void QNextEvent( QEvent *eventp )

{
	int configure;
	int width=0,height=0;
	int i,loop;
	int button;
	KeySym mykey;
	char c;

	eventp->type = QNullEvent;

	XNextEvent( MyDisplay , &XActualEvent );

	switch (XActualEvent.type) {

	case MotionNotify:
		if( !( ActualEventMask & QPointerMoveMask ) ) break;
		eventp->type = QPointerMove;
		eventp->button.x = XActualEvent.xmotion.x;
		eventp->button.y = XActualEvent.xmotion.y;
		break;
	case ButtonPress:
	case ButtonRelease:
		if( !( ActualEventMask & QButtonPressMask ) ) break;
		if( XActualEvent.xbutton.window != MyWindow ) break;
		eventp->type = QButtonPress;
		eventp->button.x = XActualEvent.xbutton.x;
		eventp->button.y = XActualEvent.xbutton.y;
		eventp->button.press = ( XActualEvent.type == ButtonPress ) ?
					QButtonDown : QButtonUp;
		button = XActualEvent.xbutton.button;
		switch (XActualEvent.xbutton.button) {
		case X_BUTTON_LEFT:
			eventp->button.button = QButtonLeft;
			break;
		case X_BUTTON_RIGHT:
			eventp->button.button = QButtonRight;
			break;
		case X_BUTTON_MIDDLE:
			eventp->button.button = QButtonMiddle;
			break;
		case X_BUTTON_WHEEL_UP:
			eventp->button.button = QButtonWheelUp;
			break;
		case X_BUTTON_WHEEL_DOWN:
			eventp->button.button = QButtonWheelDown;
			break;
		default:
			printf("QNextEvent: unknown button %d\n",button);
			eventp->button.button = QButtonUnknown;
			break;
		}
		break;
	case KeyPress:
		if( !( ActualEventMask & QKeyPressMask ) ) break;
		eventp->type = QKeyPress;
                i=XLookupString(&(XActualEvent.xkey),&c,1,&mykey,0);
		useless = i;
/*
                if( mykey == XK_Return )
                        c='\n';
*/
		if( mykey > 255 ) 
			eventp->key = mykey - 0xF000;
		else
			eventp->key = c;
		break;
	case MappingNotify:
		XRefreshKeyboardMapping(&(XActualEvent.xmapping));
		eventp->type = QMappingNotify;
		break;
	case Expose:
	case ConfigureNotify:
		configure=0;
		do {
			if( XActualEvent.type == ConfigureNotify ) {
				configure++;
                                width=XActualEvent.xconfigure.width;
                                height=XActualEvent.xconfigure.height;
			}
			if( XCheckTypedEvent(MyDisplay,ConfigureNotify
					,&XActualEvent) )
				loop = TRUE;
			else if ( XCheckTypedEvent(MyDisplay,Expose
					,&XActualEvent) )
				loop = TRUE;
			else
				loop = FALSE;
		} while ( loop );

		if( configure ) {
			eventp->type = QConfigureNotify;
			eventp->configure.width = width;
			eventp->configure.height = height;
		} else {
			eventp->type = QExpose;
		}
		break;
	default:
		eventp->type = QNullEvent;
		eventp->aux  = (int) XActualEvent.type;
		break;
	}
}
		
