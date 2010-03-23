		
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

