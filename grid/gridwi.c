
/************************************************************************\ 
 *									*
 * gridwi.c - window manipulation routines				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
 *									*
 * Revision History:							*
 * 07-May-94: added calls MouseHide() and MouseShow()                   *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include "mouse.h"
#include "grid.h"
#include "graph.h"

void MakeTotalWindow( Rect *gbp )

{
	MouseHide();
	QViewport(XTMin,YTMin,XTMax,YTMax);
	ScalePlotWindow(XTMin,YTMin,XTMax,YTMax,gbp);
	QWindow(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y);
	QRectFill(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y,TotCol);
	QNewPen(PlotCol);
	PlotShadeRect(gbp,BorderLightCol,BorderDarkCol);
	MouseShow();
}

void MakeMenuWindow( Rect *gbp )

{
	MouseHide();
        QViewport(XMMin,YMMin,XMMax,YMMax);
	ScalePlotWindow(XMMin,YMMin,XMMax,YMMax,gbp);
	QWindow(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y);
	QRectFill(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y,MenCol);
	PlotShadeRect(gbp,BorderDarkCol,BorderLightCol);
	MouseShow();
}

void MakePlotWindow( Rect *gbp )

{
	MouseHide();
	QViewport(XPMin,YPMin,XPMax,YPMax);
	ScalePlotWindow(XPMin,YPMin,XPMax,YPMax,gbp);
	QWindow(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y);
	QRectFill(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y,PlotWinCol);
	PlotShadeRect(gbp,BorderDarkCol,BorderLightCol);
	MouseShow();
}

void MakeMessageWindow( Rect *gbp )

{
	MouseHide();
        QViewport(XSMin,YSMin,XSMax,YSMax);
        ScalePlotWindow(XSMin,YSMin,XSMax,YSMax,gbp);
	QWindow(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y);
	QRectFill(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y,MesCol);
	PlotShadeRect(gbp,BorderDarkCol,BorderLightCol);
	MouseShow();
}

void MakeCommandWindow( Rect *gbp )

{
	MouseHide();
        QViewport(XCMin,YCMin,XCMax,YCMax);
        ScalePlotWindow(XCMin,YCMin,XCMax,YCMax,gbp);
	QWindow(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y);
        QRectFill(gbp->low.x,gbp->low.y,gbp->high.x,gbp->high.y,ComCol);
	PlotShadeRect(gbp,BorderDarkCol,BorderLightCol);
	MouseShow();
}
