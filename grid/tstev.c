
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

		
#include <stdio.h>
#include "events.h"

void QTestEvent( void )

{
	QEvent event;

	printf("Press 'q' to finish\n");

	for(;;) {
		QNextEvent(&event);

		if(event.type == QButtonPress) {
			printf("ButtonPress ");
			if(event.button.press==QButtonDown) {
				printf("Down ");
			} else {
				printf("UP   ");
			}
			if(event.button.button==QButtonLeft) {
				printf("Left   ");
			} else if(event.button.button==QButtonRight) {
				printf("Right  ");
			} else if(event.button.button==QButtonMiddle) {
				printf("Middle ");
			}
			printf("%4d %4d\n",event.button.x,event.button.y);
		} else if(event.type == QPointerMove) {
			printf("PointerMove ");
			printf("%4d %4d\n",event.button.x,event.button.y);
		} else if(event.type == QKeyPress) {
			printf("KeyPress ");
			printf("%c\n",event.key);
			if(event.key=='q') break;
		} else if(event.type == QExpose) {
			printf("Expose\n");
		} else if(event.type == QNullEvent) {
		} else {
			printf("Unknown\n");
		}
	}
}


void main() 

{ 
	QInitEvent(); 
	QAddEvent(QPointerMoveMask); 
	QTestEvent(); 
	QDeInitEvent(); 
}

