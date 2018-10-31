
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


#include "graph.h"
#include "xgraph.h"
#include "keybd.h"

void main()

{
	float dx,dd,x,y;
	float xmax,ymax;
	int i;
	int ixmax,iymax;
	int ncol;
	int *g2bcol,*y2gcol,*r2ycol,*r2bcol,*bbbcol;

	QGraphInit();

	press_any_key();

	QNewPen(10);
	QWindowMaxXY(&ixmax,&iymax);
	QViewport(0,0,ixmax,iymax);
	QWindow(0.,0.,10.,10.);
	QLine(1.,1.,8.,9.);
	QFlush();

	press_any_key();

	ncol=20;
	xmax=10.;
	ymax=10.;
	dx = xmax/ncol;
	dd = dx/4;

	g2bcol = QAllocGreen2BlueColors(ncol);
	y2gcol = QAllocYellow2GreenColors(ncol);
	r2ycol = QAllocRed2YellowColors(ncol);
	r2bcol = QAllocRed2BlueColors(ncol);
	bbbcol = QAllocBlueColors(ncol);

	x = -dx/2;
	for(i=0;i<ncol;i++) {
		x += dx;
		y=1.;
		QRectFill(x-dd,y-dd,x+dd,y+dd,g2bcol[i]);
		y=3.;
		QRectFill(x-dd,y-dd,x+dd,y+dd,y2gcol[i]);
		y=5.;
		QRectFill(x-dd,y-dd,x+dd,y+dd,r2ycol[i]);
		y=7.;
		QRectFill(x-dd,y-dd,x+dd,y+dd,r2bcol[i]);
		y=9.;
		QRectFill(x-dd,y-dd,x+dd,y+dd,bbbcol[i]);
	}
	QFlush();
	
	press_any_key();
}
