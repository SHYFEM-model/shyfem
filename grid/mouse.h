
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
 * mouse.h - mouse driver routines					*
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

