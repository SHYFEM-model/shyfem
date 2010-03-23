
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
