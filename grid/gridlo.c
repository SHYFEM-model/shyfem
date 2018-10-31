
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
 * gridlo.c - central routine performing loop for input			*
 *									*
 * Revision History:							*
 * 07-Oct-2015: implement smooth zoomin with wheel			*
 * 05-Mar-2014: new action for drag and mouse wheel			*
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
    int button,press,horiz,verti;
    static int drag = 0;
    float x1,y1,x2,y2;
    static int horiz_down = 0;
    static int verti_down = 0;
    static int button_orig = 0;		/* use original functionality */
    static int plot_on_move = 1;	/* replot while dragging */
    Rect *gp;

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

	case QPointerMove:

		horiz = event.button.x;
		verti = event.button.y;
		//printf("LoopForInput: QPointerMove %d %d\n",horiz,verti);

		if( plot_on_move ) {
			GetPlotCoord(horiz_down,verti_down,&x1,&y1);
			GetPlotCoord(horiz,verti,&x2,&y2);
			GfMoveRelative(x1-x2,y1-y2);
			horiz_down = horiz;
			verti_down = verti;
			drag = 1;	/* drag event already handled */
			QSkipMotion();
		}
		break;

	case QButtonPress:

		if( button_orig ) {
		    if( event.button.press == QButtonUp ) break;	
		}

		press = event.button.press;
		button = event.button.button;
		horiz = event.button.x;
		verti = event.button.y;

		//printf("LoopForInput: QButtonPress %d %d\n",button,press);

		if( InMenuField(horiz,verti) ) {
		    if( button == QButtonLeft ) {
			// printf("LoopForInput: MenuFieldInput %d\n",button);
			ActMode = MENU_FIELD_INPUT;
			ProcessMenuInput(horiz,verti);
		    }
		    break;
		} else if( ! InPlotField(horiz,verti) ) {
		    break;
		}

		if( !button_orig ) {
		    if( press == QButtonDown ) {
			QAddEvent( QPointerMoveMask );
			horiz_down = horiz;
			verti_down = verti;
			drag = 0;
			break;
		    } if( press == QButtonUp ) {
		      QDeleteEvent( QPointerMoveMask );
		      if( drag || horiz != horiz_down || verti != verti_down ) {
		        if( button == QButtonLeft ) {	/* drag event */
				GetPlotCoord(horiz_down,verti_down,&x1,&y1);
				GetPlotCoord(horiz,verti,&x2,&y2);
				/*
		                printf("LoopForInput - drag event: %d %f %f\n"
						,button,x1-x2,y1-y2);
				*/
				GfMoveRelative(x1-x2,y1-y2);
				break;
		        }
		      }
		    }
		}

		if( button == QButtonLeft ) {
			ActMode = PLOT_FIELD_INPUT;
			PlotFieldInput(horiz,verti,LEFT_MOUSE_BUTTON);
		} else if( button == QButtonRight ) {
			PlotFieldInput(horiz,verti,RIGHT_MOUSE_BUTTON);
		} else if( button == QButtonWheelUp ) {
			GetPlotCoord(horiz,verti,&ActX,&ActY);
			gp = GetActPlotWindow();
			MakeNewCenter(gp,&ActX,&ActY,0.5);
			GfZoomIn();
		} else if( button == QButtonWheelDown ) {
			GetPlotCoord(horiz,verti,&ActX,&ActY);
			gp = GetActPlotWindow();
			MakeNewCenter(gp,&ActX,&ActY,2.0);
			GfZoomOut();
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
