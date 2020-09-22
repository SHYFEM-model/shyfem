
/************************************************************************\
 *
 *    Copyright (C) 1992,1994  Georg Umgiesser
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
 *
 * gridwi.c - window manipulation routines
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 07.05.1994	ggu	added calls MouseHide() and MouseShow()
 *
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
