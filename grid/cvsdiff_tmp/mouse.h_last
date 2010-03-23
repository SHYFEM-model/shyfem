
/************************************************************************\ 
 *									*
 * mouse.h - mouse driver routines					*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see mouse.c or xmouse.c for copying information			*
 *									*
 * Revision History:							*
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_MOUSE_
#define __GUH_MOUSE_


#define NO_MOUSE_BUTTON     0
#define LEFT_MOUSE_BUTTON   1
#define RIGHT_MOUSE_BUTTON  2
#define MIDDLE_MOUSE_BUTTON 4


extern void MouseInit(void);
extern void MouseShow(void);
extern void MouseHide(void);
extern void MouseStatus( int *mbutton , int *mhoriz , int *mverti );
extern void MouseMove( int mhoriz , int mverti );
extern void MousePress( int mbutton, int *mpress, int *mhoriz, int *mverti );
extern void MouseRelease( int mbutton, int *mpress, int *mhoriz, int *mverti );
extern void MouseWindow( int xmin , int ymin , int xmax , int ymax );
extern void MouseHole( int xmin , int ymin , int xmax , int ymax );
extern void MouseInfo( int *vers , int *mtyp , int *mirq );


#endif

