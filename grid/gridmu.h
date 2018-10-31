
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
 * gridmu.h - menu routines   *
 *									*
 * Revision History:							*
 * 13-Apr-94: completely restructured -> uses hash.c to do work         *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRIDMU_
#define __GUH_GRIDMU_


#include "grid_ty.h"


typedef struct popup_item_tag {
    char *text;
    Button_type id;
    struct popup_item_tag *submenu;
    struct popup_item_tag *next;
} Popup_item_type;

Button_type CallPopup( Button_type button , int horiz , int verti );
void MakePopups( void );
Popup_item_type *MakePopupMenuItem( char *text , Button_type id );
void AddPopupSubmenu( Popup_item_type *mphead , Popup_item_type *mp );
Button_type PopupMenu( Popup_item_type *mp , int x , int y );
void PlotMenuItem( char *s , float x0 , float y0 , float dx , float dy
                    , int fgcol , int bgcol , int bdcol );


void *SavePixels( int x , int y , int width , int height );
void RestorePixels( int x , int y , void *buffer );

#endif
