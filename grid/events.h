
/************************************************************************\ 
 *									*
 * events.h - header file for event routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see events.c for copying information					*
 *									*
 * Revision History:							*
 * 05-Mar-2014: handle also mouse wheel					*
 * 05-Dec-95: QKey... symbols introduced                                *
 * 04-Dec-94: routines written from scratch                             *
 *									*
\************************************************************************/


#ifndef __GUH_EVENTS_
#define __GUH_EVENTS_

#define QKeyPressMask		0x001
#define QButtonPressMask	0x002
#define QPointerMoveMask	0x004
#define QExposureMask		0x008
#define QStructureNotifyMask	0x010

#define QButtonUnknown		0x000
#define QButtonLeft		0x001
#define QButtonMiddle		0x002
#define QButtonRight		0x003
#define QButtonWheelUp		0x004
#define QButtonWheelDown	0x005
/*
#define QButtonRight		0x002
#define QButtonMiddle		0x004
*/

#define QButtonDown		0x001
#define QButtonUp		0x002

#define QKeyUnknown		0x0000
#define QKeyLeft		0x0F51
#define QKeyUp			0x0F52
#define QKeyRight		0x0F53
#define QKeyDown		0x0F54
#define QKeyBackSpace		0x0F08
#define QKeyDelete		0x0FFF
#define QKeyTab			0x0F09
#define QKeyReturn		0x0F0D
#define QKeyEscape		0x0F1B

typedef unsigned long QEventMask;

typedef enum {
	QNullEvent ,
	QKeyPress , QButtonPress , QPointerMove ,
	QMappingNotify , QConfigureNotify , QExpose
} QEventType;

typedef struct {
	QEventType type;
	struct {
		char button;
		char press;
		int  x;
		int  y;
	} button;
	int key;
	struct {
		int width;
		int height;
	} configure;
	void *custom;
	int aux;
} QEvent;

void QInitEvent( void );
void QDeInitEvent( void );

QEventMask QSelectedEvent( void );
void QSelectEvent( QEventMask eventmask );

void QAddEvent( QEventMask eventmask );
void QDeleteEvent( QEventMask eventmask );

void QSkipMotion( void );
void QNextEvent( QEvent *eventp );

/* QCheckEvent, QSelectedEvent */

#endif /* __GUH_EVENTS_ */
