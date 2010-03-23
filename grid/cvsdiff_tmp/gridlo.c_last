
/************************************************************************\
 *									*
 * gridlo.c - central routine performing loop for input			*
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
 * 02-Apr-1998: new LoopForEvents, ExitEventLoop() (used to exit loop)  *
 * 02-Apr-1998: ProcessMenuInput() called for new menu routines         *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include "grid.h"
#include "graph.h"
#include "general.h"
#include "gustd.h"
#include "keybd.h"
#include "mouse.h"

#include "events.h"
#include "menu.h"


static int event_debug = FALSE;		/* set to TRUE for debugging */

static int LoopForEvents = TRUE;

void ExitEventLoop( void ) { LoopForEvents = FALSE; }


void CheckEventInput( QEvent event );

void LoopForInput( void )

{
    QEvent event;
    int button,horiz,verti;

    while( LoopForEvents ) {

	QNextEvent( &event );

	if( event_debug ) CheckEventInput( event ); 

	switch( event.type ) {

	case QConfigureNotify:

		ResizeWindow(event.configure.width,event.configure.height);
		RedrawAll();
		break;

	case QExpose:

		RedrawAll();
		break;

	case QKeyPress:

		KeyboardInput( event.key );
		break;

	case QButtonPress:

		if( event.button.press == QButtonUp ) break;

		button = event.button.button;
		horiz = event.button.x;
		verti = event.button.y;

		if( button == QButtonLeft ) {
			if( InPlotField(horiz,verti) ) {
				ActMode = PLOT_FIELD_INPUT;
				PlotFieldInput(horiz,verti,LEFT_MOUSE_BUTTON);
			}
			if( InMenuField(horiz,verti) ) {
				/* MenuFieldInput(horiz,verti); */
				ActMode = MENU_FIELD_INPUT;
				ProcessMenuInput(horiz,verti);
			}
		} else if( button == QButtonRight ) {
			if( InPlotField(horiz,verti) ) {
				PlotFieldInput(horiz,verti,RIGHT_MOUSE_BUTTON);
			}
		}

		break;

	default:

		break;
	}
    }
}


void CheckEventInput( QEvent event )

{
	switch( event.type ) {

	case QConfigureNotify:

		printf("QConfigureNotify : %d %d\n",
			event.configure.width,event.configure.height);
		break;

	case QExpose:

		printf("QExpose\n");
		break;

	case QKeyPress:

		printf("QKeyPress : %c\n",event.key);
		break;

	case QButtonPress:

		printf("QButtonPress : %d %d %d %d\n",
			event.button.button,event.button.press,
			event.button.x,event.button.y);
		break;

	case QPointerMove:

		printf("PointerMove : %d %d\n",
			event.button.x,event.button.y);
		break;

	case QNullEvent:

		printf("QNull : %d\n",event.aux);
		break;

	default:

		printf("QUnknown\n");
		break;

	}
}
