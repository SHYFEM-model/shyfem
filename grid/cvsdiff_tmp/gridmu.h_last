
/************************************************************************\
 *									*
 * gridmu.h - menu routines   *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see gridmu.c for copying information         *
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
