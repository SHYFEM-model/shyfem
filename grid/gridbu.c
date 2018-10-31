
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
 * gridbu.c - button manipulation routines				*
 *									*
 * Revision History:							*
 * 09-Feb-1998: ActArgument eliminated, new functions GfZoom, GfShow    *
 * 13-Oct-97: new registered functions GfRemoveElement, GfRemoveLine    *
 * 10-Oct-97: new registered functions GfSave, GfUnifyNode              *
 * 04-Dec-95: new functions for Vector registered                       *
 * 16-Apr-94: button routines restructured -> is now independent        *
 *                BLP is now local, Button_list structures are local    *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdlib.h>

#include "grid.h"
#include "graph.h"
#include "mouse.h"

#include "list.h"


typedef struct button_list_tag {
	Button_type button;
	char *button_text;
	Rect minmax;
	struct button_list_tag *next;
} Button_list_type;

typedef struct {	/* use this to implement with List_type */
	Button_type button;
	char *button_text;
	Rect minmax;
} Button_description_type;

/*******************************************************************/
/*
typedef void (*Function_type) (void);
*/
typedef struct {
	Button_type button;
	Function_type fp;
	int arg;		/* <- this one is not used anymore */
} Function_desc_type;

Listtable_type LFP = NULL;

void RegisterFunction( Button_type button , Function_type fp , int arg )

{
    Function_desc_type *new;

    new = (Function_desc_type *) malloc ( sizeof(Function_desc_type) );
    if( new == NULL ) Error("Cannot alloacte function structure");

    new->button = button;
    new->fp     = fp;
    new->arg    = arg;

    InsertListTable( LFP , (void *) new );
}

Function_desc_type *RecallFunction( Button_type button )

{
    Function_desc_type *fd;

    ResetListTable( LFP );

    while( (fd=(Function_desc_type *)VisitListTable( LFP )) != NULL )
	if( fd->button == button ) break;

    return fd;
}

void SetNewCommand( Button_type button )

{
    Function_desc_type *fd;

    fd = RecallFunction( button );

    if( fd ) {
        ActCommand  = fd->button;
        ActFunction = fd->fp;
    } else {
        ActCommand  = NONE;
        ActFunction = NULL;
    }
}

void InitializeFunctions( void )

{
    if( LFP ) return;

    LFP = MakeListTable();

    RegisterFunction( CANCEL           , GfCancel          , 0 );
    RegisterFunction( REFRESH          , GfRefresh         , 0 );
    RegisterFunction( PRINT            , GfPrint           , 0 );
    RegisterFunction( SAVE             , GfSave            , 0 );
    RegisterFunction( EXIT             , GfExit            , 0 );

    RegisterFunction( ZOOM_WINDOW      , GfZoomWindow      , 0 );
    RegisterFunction( ZOOM_IN          , GfZoomIn          , 0 );
    RegisterFunction( ZOOM_OUT         , GfZoomOut         , 0 );
    RegisterFunction( MOVE             , GfMove            , 0 );
    RegisterFunction( TOTAL_VIEW       , GfTotalView       , 0 );

    RegisterFunction( SHOW_NODE        , GfShowNode        , 0 );
    RegisterFunction( SHOW_ELEMENT     , GfShowElement     , 0 );
    RegisterFunction( SHOW_LINE        , GfShowLine        , 0 );
    RegisterFunction( SHOW_VECT        , GfShowVect        , 0 );

    RegisterFunction( MAKE_NODE        , GfMakeNode        , 0 );
    RegisterFunction( DEL_NODE         , GfDelNode         , 0 );
    RegisterFunction( MOVE_NODE        , GfMoveNode        , 0 );
    RegisterFunction( UNIFY_NODE       , GfUnifyNode       , 0 );

    RegisterFunction( MAKE_ELEMENT     , GfMakeElement     , 0 );
    RegisterFunction( DEL_ELEMENT      , GfDelElement      , 0 );
    RegisterFunction( REMOVE_ELEMENT   , GfRemoveElement   , 0 );

    RegisterFunction( MAKE_LINE        , GfMakeLine        , 0 );
    RegisterFunction( DEL_LINE         , GfDelLine         , 0 );
    RegisterFunction( REMOVE_LINE      , GfRemoveLine      , 0 );
    RegisterFunction( JOIN_LINE        , GfJoinLine        , 0 );
    RegisterFunction( SPLIT_LINE       , GfSplitLine       , 0 );

    RegisterFunction( DEL_VECT         , GfDelVect         , 0 );
    RegisterFunction( CHANGE_VECT      , GfChangeVect      , 0 );

    RegisterFunction( CHANGE_DEPTH     , GfChangeDepth     , 0 );
    RegisterFunction( CHANGE_TYPE      , GfChangeType      , 0 );

}

/*******************************************************************/

Button_list_type *FindButton( Button_type button );
void PlotButtonP( Button_list_type *bp , int color );
void ButtonList( Rect *r , char *s , Button_type button );
Button_list_type *MakeButton( Rect *r , char *s , Button_type button );
void Button( float xmin , float ymin , float width , float height
				, char *s , Button_type button);


static Button_list_type *BLP = NULL;


void PlotButtons( void )

{
	Button_list_type *bp;

        QViewport(XMMin,YMMin,XMMax,YMMax);
	QWindow(GbMen.low.x,GbMen.low.y,GbMen.high.x,GbMen.high.y);
	QRectFill(GbMen.low.x,GbMen.low.y,GbMen.high.x,GbMen.high.y,MenCol);
	PlotShadeRect(&GbMen,BorderDarkCol,BorderLightCol);

	for( bp=BLP ; bp!=NULL ; bp=bp->next )
           if( bp->button == ActButton )
               PlotButtonP( bp , ButChoosCol );
           else
               PlotButtonP( bp , ButCol );
}

Button_list_type *FindButton( Button_type button )

{
	Button_list_type *bp;

	for( bp=BLP ; bp!=NULL ; bp=bp->next )
		if( bp->button == button ) return bp;

	return NULL;
}

Button_type FindButtonByCoord( float x , float y )

{
	Button_list_type *bp;

	for( bp=BLP ; bp!=NULL ; bp=bp->next )
		if( bp->minmax.low.x <= x && bp->minmax.high.x >= x )
			if( bp->minmax.low.y <= y && bp->minmax.high.y >= y )
				return bp->button;

	return NONE;
}

void PlotButton( Button_type button , int color )

{
	PlotButtonP( FindButton( button ) , color );
}

void PlotButtonP( Button_list_type *bp , int color )

{
	Rect *br;
	int tmpcol;
	float width,height;

	if( bp == NULL ) return;

	br = &(bp->minmax);

	QViewport(XMMin,YMMin,XMMax,YMMax);
	QWindow(GbMen.low.x,GbMen.low.y,GbMen.high.x,GbMen.high.y);

	MouseHide();

	QRectFill(br->low.x,br->low.y,br->high.x,br->high.y,color);
	PlotShadeRect(br,BorderLightCol,BorderDarkCol);

	QGetPen(&tmpcol);
	QNewPen(ButWriteCol);
	QTextDimensions(bp->button_text,&width,&height);
	QTextBackground(color);
	QText( (br->high.x + br->low.x - width)/2.
		, (br->high.y + br->low.y - height)/2.
               	, bp->button_text );
	QTextBackground(Black);
	QNewPen(tmpcol);

	MouseShow();

}

void ButtonList( Rect *r , char *s , Button_type button )

{
	Button_list_type *bp;

	bp = MakeButton(r,s,button);
        bp->next = BLP;
        BLP = bp;
}

Button_list_type *MakeButton( Rect *r , char *s , Button_type button )

{
        Button_list_type *new;

        new = (Button_list_type *) malloc( sizeof( Button_list_type ) );

        if( !new ) Error("MakeButton : Cannot allocate node");

	new->button = button;
	new->button_text = s;
	new->minmax = *r;
        new->next = NULL;

        return new;
}

void Button( float xmin , float ymin , float width , float height
				, char *s , Button_type button)

/* insert buttons into list */

{
        Rect br;

	br.low.x  = xmin;
	br.high.x = xmin + width;
	br.low.y  = ymin;
	br.high.y = ymin + height;

	ButtonList(&br,s,button);
}

void MakeButtons( Rect *gbp )

{
	int nmenu;
        float width,height;
        float xmargin,ymargin;
        float bwidth,bheight;
        float xmin,ymin;
        float dy;

	if( BLP ) return;	/* already made */

	nmenu = 12;	/* number of maximal menu buttons */

        width  = gbp->high.x - gbp->low.x;
	height = gbp->high.y - gbp->low.y;


	xmargin = 0.1*width;
        bwidth  = width - 2.*xmargin;

	dy      = height/(float) nmenu;	/* size in y of button (incl. margin) */
	ymargin = 0.05*dy;
        bheight = dy - 2.*ymargin;

        xmin = gbp->low.x  + xmargin;
        ymin = gbp->high.y + ymargin;

	Button(xmin,(ymin-=dy),bwidth,bheight,"File",FILE_MENU);
	Button(xmin,(ymin-=dy),bwidth,bheight,"View",VIEW);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Show",SHOW);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Node",NODE);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Element",ELEMENT);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Line",LINE);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Vector",VECTOR);
	Button(xmin,(ymin-=dy),bwidth,bheight,"Change",CHANGE);

}
