
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
 * gridmu.c - pop up menu routines                                      *
 *									*
 * Revision History:							*
 * 13-Oct-97: new menu items Remove Element, Remove Line                *
 * 10-Oct-97: new menu items Save, Unify Node                           *
 * 04-Dec-95: new popup menus for Vector added                          *
 * 28-Oct-94: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdlib.h>

#include "grid.h"
#include "gustd.h"
#include "graph.h"

#include "mouse.h"
#include "events.h"

#include "gridmu.h"


static Popup_item_type *MpHead=NULL;

static int Permanent=FALSE; /* false for drag-down menu, true for pop-up */

static int FgCol;
static int BgCol;
static int BdCol;


static void FindPopupDimension( Popup_item_type *mp , int *nentry
				, int *w , int *h );
static void FindPopupCoordinates( int width , int height
				, int *xc , int *yc );
static char **GetPopupText( Popup_item_type *mp , int nmenus );
static Button_type FindPopupEntry( Popup_item_type *mp , int i );


Button_type CallPopup( Button_type button , int horiz , int verti )

{
    Popup_item_type *mpaux;

    mpaux = MpHead->submenu;
    while( mpaux ) {
        if( mpaux->id == button ) break;
        mpaux = mpaux->next;
    }

    if( mpaux ) {
        return PopupMenu( mpaux , horiz , verti );
    } else {
        return NONE;
    }
}


void MakePopups( void )

{
    Popup_item_type *mpaux;

    if( MpHead ) return;

    FgCol=Black;
    BgCol=Green;
    BdCol=White;

/****************** TOP ********************/

    MpHead = MakePopupMenuItem("Top",TOP);

/****************** FILE ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("File",FILE_MENU));

    mpaux = MpHead->submenu;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Cancel",CANCEL));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Refresh",REFRESH));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Print",PRINT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Save",SAVE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Exit",EXIT));

/****************** VIEW ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("View",VIEW));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Zoom Window",ZOOM_WINDOW));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Zoom in",ZOOM_IN));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Zoom out",ZOOM_OUT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Total view",TOTAL_VIEW));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Move",MOVE));

/****************** SHOW ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Show",SHOW));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Show Node",SHOW_NODE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Show Element",SHOW_ELEMENT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Show Line",SHOW_LINE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Show Vector",SHOW_VECT));

/****************** NODE ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Node",NODE));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Make Node",MAKE_NODE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Del Node",DEL_NODE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Move Node",MOVE_NODE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Unify Node",UNIFY_NODE));

/****************** ELEMENT ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Element",ELEMENT));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Make Element",MAKE_ELEMENT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Del Element",DEL_ELEMENT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Remove Element",REMOVE_ELEMENT));

/****************** LINE ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Line",LINE));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Make Line",MAKE_LINE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Del Line",DEL_LINE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Remove Line",REMOVE_LINE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Split Line",SPLIT_LINE));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Join Line",JOIN_LINE));

/****************** Vector ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Vector",VECTOR));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Del Vector",DEL_VECT));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Change Vector",CHANGE_VECT));

/****************** CHANGE ********************/

    AddPopupSubmenu(MpHead,MakePopupMenuItem("Change",CHANGE));

    mpaux = mpaux->next;

    AddPopupSubmenu(mpaux,MakePopupMenuItem("Change Depth",CHANGE_DEPTH));
    AddPopupSubmenu(mpaux,MakePopupMenuItem("Change Type",CHANGE_TYPE));

}

Popup_item_type *MakePopupMenuItem( char *text , Button_type id )

{
    Popup_item_type *new;
    char *s;

    new = (Popup_item_type *) malloc( sizeof(Popup_item_type) );
    if( !new ) Error("MakePopupMenuItem : Cannot allocate node");

    s = savestring( text , -1 );
    if( !s ) Error("MakePopupMenuItem : Cannot allocate string");
    new->text = s;

    new->id      = id;
    new->submenu = NULL;
    new->next    = NULL;

    return new;
}

void AddPopupSubmenu( Popup_item_type *mphead , Popup_item_type *mp )

{
    if( mphead->submenu == NULL ) {
        mphead->submenu = mp;
    } else {
        mphead = mphead->submenu;
        while( mphead->next )
            mphead = mphead->next;
        mphead->next = mp;
    }
}

void CheckEventInput( QEvent event );

Button_type PopupMenu( Popup_item_type *mp , int x , int y )

{
    int width,height,totheight;
    int nmenus=0;

    float x0,y0,dx,dy;

    int i,iold;
    int loop,changed,pressed;
    int mhoriz,mverti;
    char **texts;
    void *pixmap;

    QEvent event;

    FindPopupDimension(mp,&nmenus,&width,&height);

    width = (width*5)/4;
    height = (height*9)/4;
    totheight = height*nmenus;

    FindPopupCoordinates(width,totheight,&x,&y);

    pixmap=QSavePixels(x,y,width,totheight);

    QViewport(x,y,x+width-1,y+totheight-1);
    QWindow(0.,0.,(float)width,(float)totheight);

    pressed=TRUE;
    if( Permanent ) { /* wait till button is released */
      for(;;) {
	QNextEvent( &event );
	if( event.type == QButtonPress && event.button.press == QButtonUp )
		break;
      }
      pressed=FALSE;
    }

    dx=width;
    dy=height;
    x0=0.;
    y0=totheight-dy;

    texts = GetPopupText(mp,nmenus);

    for( i=0 ; i<nmenus ; i++,y0-=dy )
	PlotMenuItem(texts[i],x0,y0,dx,dy,FgCol,BgCol,BdCol);

    QAddEvent( QPointerMoveMask );

    iold = -1;
    do {
	QNextEvent( &event );
/*	CheckEventInput( event ); */
	if( event.type == QButtonPress ) {
	  if( event.button.button == QButtonLeft ) {
	    pressed = (event.button.press==QButtonDown) ? TRUE : FALSE;
	  }
	  mhoriz=event.button.x;
	  mverti=event.button.y;
	} else if( event.type == QPointerMove ) {
	  mhoriz=event.button.x;
	  mverti=event.button.y;
	} else {
	  /* process other events !!!! */
	  mhoriz = -1;
	  mverti = -1;
	}

	i = -1;
	if(mhoriz>=x && mhoriz<x+width) {
	    if(mverti>=y && mverti<y+totheight) {
		i=(mverti-y)/height;
	    }
	}

	if( i != iold ) {
	    if( iold != -1 ) {
	      y0=totheight-(iold+1)*dy;
	      PlotMenuItem(texts[iold],x0,y0,dx,dy,FgCol,BgCol,BdCol);
	    }
	    if( i != -1 ) {
	      y0=totheight-(i+1)*dy;
	      PlotMenuItem(texts[i],x0,y0,dx,dy,BgCol,FgCol,BdCol);
	    }
	    iold=i;
	    changed=TRUE;
	} else {
	    changed=FALSE;
	}

	if( Permanent ) {
	    loop = ( !pressed || changed );
	} else {
	    loop = pressed;
	    if( !loop )
		iold = changed ? -1 : iold;
	}

    } while( loop );

    QDeleteEvent( QPointerMoveMask );

    if( pixmap )
        QRestorePixels( x , y , pixmap );
    else
	RedrawAll();

    free(texts);
    QNewPen(PlotCol);

    return FindPopupEntry(mp,iold);
}

static void FindPopupDimension( Popup_item_type *mp , int *nentry
				, int *w , int *h )

{
    int maxwidth=0;
    int maxheight=0;
    int width,height;
    int n=0;

    mp = mp->submenu;
    while( mp ) {
        QTextDimensionsI(mp->text,&width,&height);
        maxwidth = width>maxwidth ? width : maxwidth;
        maxheight = height>maxheight ? height : maxheight;
        n++;
        mp = mp->next;
    }
    *nentry=n;
    *w=maxwidth;
    *h=maxheight;
}

static void FindPopupCoordinates( int width , int height
				, int *xc , int *yc )

{
    /* height is total height of menu */

    int x,y;
    int xminpix=XTMin;
    int yminpix=YTMin;
    int xmaxpix=XTMax;
    int ymaxpix=YTMax;

    x = *xc;
    y = *yc;

    x = x-width+1 - 5;
    if( x < xminpix ) {
	x = x+width-1 + 10;
	if( x+width-1 > xmaxpix ) {
	    if( width > xmaxpix-xminpix+1 )
                Error("No space for submenu");
            else
                x = ((xmaxpix-xminpix+1)-width)/2;
        }
    }

    y += 5;
    if( y+height-1 > ymaxpix ) {
        y = y-height+1 - 10;
        if( y < yminpix ) {
            if( height > ymaxpix-yminpix+1 )
                Error("No space for submenu");
            else
                y = ((ymaxpix-yminpix+1)-height)/2;
        }
    }

    *xc = x;
    *yc = y;
}

static char **GetPopupText( Popup_item_type *mp , int nmenus )

{
    char **texts;
    int i=0;

    texts = (char **) malloc( nmenus * sizeof(char *) );
    if( texts == NULL ) Error("Cannot allocate submenu");

    mp = mp->submenu;
    while( mp ) {
	texts[i++] = mp->text;
	mp = mp->next;
    }

    return texts;
}

static Button_type FindPopupEntry( Popup_item_type *mp , int i )

{
    if( i < 0 )
	return NONE;

    mp = mp->submenu;
    while( i-- && mp )
	mp = mp->next;

    if( mp )
	return mp->id;
    else
	return NONE;
}

void PlotMenuItem( char *s , float x0 , float y0 , float dx , float dy
                    , int fgcol , int bgcol , int bdcol )

{
    float width,height;

	MouseHide();

    QTextDimensions(s,&width,&height);

    QRectFill(x0,y0,x0+dx,y0+dy,bgcol);
    QTextBackground(bgcol);

    QNewPen(bdcol);
    QMove( x0 , y0 );
    QPlot( x0+dx , y0 );
    QPlot( x0+dx , y0+dy );
    QPlot( x0 , y0+dy );
    QPlot( x0 , y0 );

    QNewPen(fgcol);
    QText( x0 + (dx-width)/2. , y0 + (dy-height)/2. , s );
    QTextBackground(Black);

	MouseShow();
}



/* this goes into graph... */

/*

void *SavePixels( int x , int y , int width , int height )

{
    void *new;
    int size;

    size = imagesize(x,y,x+width-1,y+height-1);
    new = (void *) malloc( size );

    if( new ) {
	MouseHide();
	getimage(x,y,x+width-1,y+height-1,new);
	MouseShow();
    }

    return new;
}

void RestorePixels( int x , int y , void *buffer )

{
    if( buffer ) {
	MouseHide();
        putimage(x,y,buffer,COPY_PUT);
	MouseShow();
	free(buffer);
    }
}

*/
