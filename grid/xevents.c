
/************************************************************************\
 *
 *    Copyright (C) 1994-1995,2014  Georg Umgiesser
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
 * gevents.c - event routines for gcc under DOS
 *
 * revision log :
 *
 * 04.12.1994	ggu	routines written from scratch
 * 10.12.1994	ggu	routines adapted to gcc through Event Queue Library
 * 05.12.1995	ggu	handle special keys from keyboard (arrows, return...)
 * 05.03.2014	ggu	handle mouse wheel, new QSkipMotion()
 *
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
	int loop;
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
                (void) XLookupString(&(XActualEvent.xkey),&c,1,&mykey,0);
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
		
