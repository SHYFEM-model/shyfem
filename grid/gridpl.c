
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
 * gridpl.c - utilities for plotting 					*
 *									*
 * Revision History:							*
 * 17-Dec-97: plots lines with color according to line type             *
 *              uses OpColor to decide if color or not                  *
 *              (changed also algorithm for node coloring)              *
 *              new routine SetLineColor()                              *
 * 14-Oct-97: Administer use of nodes with routines                     *
 *            call SetNodeColor() in PlotPoint only if not actual node  *
 *              -> this solves bug (nodes not getting evidenced)        *
 * 10-Oct-97: New static routines PlotNodeForm(), SetNodeColor()        *
 *              PlotPoint() now respects color of nodes to plot         *
 * 06-Dec-95: ScaleFactor modified to scale with vector length          *
 * 05-Dec-95: EvidenceVect, PlotVectors, PlotVect introduced            *
 * 05-Dec-95: ScaleFactor for OpNodeFact, OpVectFact                    *
 * 02-Dec-95: PlotPoints(), PlotPoint() moved from gridma2 to here      *
 * 25-Mar-95: NULLDEPTH is not interpolated from nodes but colored      *
 *              in grey (PlotElem())                                    *
 * 06-Oct-94: no reference to ColTab anymore                            *
 * 13-May-94: in FillElement x,y was not static                         *
 * 05-May-94: PlotElements does not erase plot window -> new routine    *
 *              ErasePlot() to do this                                  *
 * 13-Apr-94: use new hash routines                                     *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdlib.h>
#include <math.h>

#include "grid.h"
#include "graph.h"
#include "mouse.h"

#include "list.h"
#include "gridhs.h"


static void PlotElemP( Hashtable_type H , Elem_type *pe );
static void PlotLineP( Hashtable_type H , Line_type *pl );
static void PlotLineSegP( Hashtable_type H , Line_type *pl , int node );

static void PfeilUV( float x, float y, float u, float v, float uv );

/**************************************************************************/

void EvidenceNode( int node , int color )

{
	int savecolor;
	Elem_type *pe;
	Line_type *pl;
	Node_type *pn;

	if( !node ) return;

	QGetPen(&savecolor);
	QNewPen(color);

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	pn = RetrieveByNodeNumber( HNN , node );

	if( !pn ) /*nothing to do */
		;
	else if( !GetUseN(pn) )
		PlotPoint( pn );
	else {
		pe = NULL;
		while( (pe = FindElemToNode( HEL , pe , node )) != NULL )
			PlotElemP( HNN , pe );

		pl = NULL;
		while( (pl = FindLineToNode( HLI , pl , node )) != NULL )
			PlotLineSegP( HNN , pl , node );
	}

	QNewPen(savecolor);

	MouseShow();
}

void EvidenceElem( int elem , int color )

{
	int savecolor;

	if( !elem ) return;

	QGetPen(&savecolor);
	QNewPen(color);

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	PlotElemP( HNN , RetrieveByElemNumber(HEL,elem) );

	QNewPen(savecolor);

	MouseShow();
}

void EvidenceLine( int line , int color )

{
	int savecolor;

	if( !line ) return;

	QGetPen(&savecolor);
	QNewPen(color);

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	PlotLineP( HNN , RetrieveByLineNumber(HLI,line) );

	QNewPen(savecolor);

	MouseShow();
}

void EvidenceVect( int vect , int color )

{
	int savecolor;

	if( !vect ) return;

	QGetPen(&savecolor);
	QNewPen(color);

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	PlotVect( RetrieveByNodeNumber(HVC,vect) );

	QNewPen(savecolor);

	MouseShow();
}

/**************************************************************************/

void GetScreenCoord( float x , float y , int *horiz , int *verti )

{
		/* compute screen coordinates in plot window */

		QViewport(XPMin,YPMin,XPMax,YPMax);
		QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);
/*		QViewport(XTMin,YTMin,XTMax,YTMax);
		QWindow(GbTot.low.x,GbTot.low.y,GbTot.high.x,GbTot.high.y);
*/		QScreenXY(x,y,horiz,verti);
}

void GetPlotCoord( int horiz , int verti , float *x , float *y )

{
		/* get coordinates of mouse cursor in plot window */

		QViewport(XPMin,YPMin,XPMax,YPMax);
		QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);
		QRealXY(horiz,verti,x,y);
}

void GetMenuCoord( int horiz , int verti , float *x , float *y )

{
		/* get coordinates of mouse cursor in menu window */

		QViewport(XMMin,YMMin,XMMax,YMMax);
		QWindow(GbMen.low.x,GbMen.low.y,GbMen.high.x,GbMen.high.y);
		QRealXY(horiz,verti,x,y);
}

/**************************************************************************/

void ScalePlotWindow( int xmin , int ymin , int xmax , int ymax , Rect *gb )

{
	float sx,sy,dxy;

	sx = (xmax-xmin)/(gb->high.x-gb->low.x);
	sy = (ymax-ymin)/(gb->high.y-gb->low.y);

	if( sx > sy ) {
	    sx = sy;
	    dxy = ( xmax - xmin ) / sx;
	    gb->low.x = ( gb->high.x + gb->low.x - dxy ) / 2.;
	    gb->high.x = gb->low.x + dxy;
	} else if( sy > sx ) {
	    sy = sx;
	    dxy = ( ymax - ymin ) / sy;
	    gb->low.y = ( gb->high.y + gb->low.y - dxy ) / 2.;
	    gb->high.y = gb->low.y + dxy;
	}
}

void ScaleFactor( Hashtable_type HV )

{
	Node_type *p;
	Vect_type *pv;
	int i;
	float fact,pdx,pdy;
	float speed=0.;

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	QGetPixelWidth(&pdx,&pdy);
	fact = 0.5*(pdx+pdy);

	ResetHashTable(HV);
	while( (p = VisitHashTableN(HV)) != NULL ) {
		pv = p->extra;
		for(i=0;i<pv->total;i++)
			if( pv->speed[i] > speed )
			    speed = pv->speed[i];
	}
	if( speed == 0. ) speed = 1.;

	/* the constant factors are just empirical */
	/* OpVectFact -> length of max vector in pixels */
	/* OpNodeFact -> half-dimension of symbol in pixels */

	OpVectFact = OpVectScal*20.*fact/speed;
	OpNodeFact = OpNodeScal*5.*fact;
}

/**************************************************************************/

void PlotRect( Rect *r )

{
	QMove( r->low.x , r->low.y );
	QPlot( r->low.x , r->high.y );
	QPlot( r->high.x , r->high.y );
	QPlot( r->high.x , r->low.y );
	QPlot( r->low.x , r->low.y );
}

void PlotShadeRect( Rect *r , int col1 , int col2 )

{
		float dx,dy;
		int color;

		QGetPixelWidth(&dx,&dy);
		QGetPen(&color);

		QNewPen(col1);
		QMove( r->low.x , r->low.y );
		QPlot( r->low.x , r->high.y );
		QPlot( r->high.x , r->high.y );
		QNewPen(col2);
		QPlot( r->high.x , r->low.y );
		QPlot( r->low.x , r->low.y );

		QNewPen(col1);
		QMove( r->low.x+dx , r->low.y+dy );
		QPlot( r->low.x+dx , r->high.y-dy );
		QPlot( r->high.x-dx , r->high.y-dy );
		QNewPen(col2);
		QPlot( r->high.x-dx , r->low.y+dy );
		QPlot( r->low.x+dx , r->low.y+dy );

		QNewPen(color);
}

/**************************************************************************/

void ErasePlot( void )

{
	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	QRectFill(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x
				,GbPlo.high.y,PlotWinCol);
	PlotShadeRect(&GbPlo,BorderDarkCol,BorderLightCol);

	MouseShow();
}

/**************************************************************************/

static void PlotNodeForm( int type , float x , float y );
static void SetNodeColor( Node_type *pn );

void PlotPoints( Hashtable_type HN )

{
	Node_type *pn;

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	ResetHashTable(HN);
	while( (pn = VisitHashTableN( HN )) != NULL ) {
	    if( GetUseN(pn) ) continue;
	    PlotPoint(pn);
	}

	QNewPen(PlotCol);
	MouseShow();
}

void PlotPoint( Node_type *p )

{
	SetNodeColor(p);
	PlotNodeForm(p->type,p->coord.x,p->coord.y);
}

static void SetNodeColor( Node_type *pn )

{
	int color;

	if( GetActNode() == pn->number ) {
	  QGetPen(&color);
	  if( color == EvidenceCol ) return;
	}

	if( OpFill || OpColor ) {
	    if( OpShowType == 0 ) {		/* color by depth */
		QNewPen( GetDepthColor(pn->depth) );
	    } else if( OpShowType == 1 ) {	/* color by type */
		QNewPen( GetTypeColor(pn->type) );
	    } else {				/* not used anymore */
		switch(pn->type) {
		case 'A':
		case 'B':
		case 'V':
		case 'W':
			QNewPen(Blue);
			break;
		case 'C':
		case 'D':
		case 'X':
		case 'Z':
			QNewPen(Cyan);
			break;
		case 'E':
		case 'F':
			QNewPen(Red);
			break;
		case 'M':			/* charts from CVN */
			QNewPen(Yellow);
			break;
		case 'N':			/* charts from CVN */
			QNewPen(Blue);
			break;
		default:
			QNewPen(PlotCol);
			break;
		}
	    }
	} else {
	    QNewPen(PlotCol);
	}
}

static void PlotNodeForm( int type , float x , float y )

{
	float dx,dy;

	dx = dy = OpNodeFact;

	switch(type) {
	case 'N':			/* charts from CVN */
	case 'A':
	case 'V':			/* plus */
		QMove( x-dx , y );
		QPlot( x+dx , y );
		QMove( x , y-dy );
		QPlot( x , y+dy );
		break;
	case 'M':			/* charts from CVN */
		dx=0.3*dx; dy=0.3*dy;
		QMove( x-dx , y-dy );
		QPlot( x+dx , y+dy );
		QMove( x-dx , y+dy );
		QPlot( x+dx , y-dy );
		break;
	case 'B':
	case 'W':			/* cross */
		QMove( x-dx , y-dy );
		QPlot( x+dx , y+dy );
		QMove( x-dx , y+dy );
		QPlot( x+dx , y-dy );
		break;
	case 'C':
	case 'X':			/* minus */
		QMove( x-dx , y );
		QPlot( x+dx , y );
		break;
	case 'D':
	case 'Z':			/* vline */
		QMove( x , y-dy );
		QPlot( x , y+dy );
		break;
	case 'E':
	case 'F':			/* rhombus */
		QMove( x-dx , y );
		QPlot( x , y+dy );
		QPlot( x+dx , y );
		QPlot( x , y-dy );
		QPlot( x-dx , y );
		break;
	default:
		QMove( x-dx , y );
		QPlot( x+dx , y );
		QMove( x , y-dy );
		QPlot( x , y+dy );
		break;
	}
}


/**************************************************************************/

void PlotElements( Hashtable_type HE , Hashtable_type HN )

{
	Elem_type *p;

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	ResetHashTable(HE);
	while( (p = VisitHashTableE(HE)) != NULL )
		PlotElem(HN,p);

	MouseShow();
}

void PlotElem( Hashtable_type H , Elem_type *p )

{
	Rect r;
	float depth;

	PolyMinMax( H , p , &r );
	if( r.low.x>GbPlo.high.x || r.high.x<GbPlo.low.x
		|| r.low.y>GbPlo.high.y || r.high.y<GbPlo.low.y ) return;

	if( OpFill ) {
		if( OpShowType == 1 ) { /* color by type */
		  FillElem(H,p,GetRandomColor((int) p->type));
		} else {		/* color by depth */
		  depth = p->depth;
/* changed temp. 25.3.95 */	/* FIXME */
		  if( depth == NULLDEPTH ) {
			depth = MakeDepthFromNodes(H,p);
		  }
/**/
		  FillElem(H,p,GetDepthColor(depth));
		}
	}

	if( OpBorder )
		PlotElemP( H , p );

}

static void PlotElemP( Hashtable_type H , Elem_type *pe )

{
    /* still to check consistency */

	Node_type *pn;
	int i,nvert;

	if( !pe ) {   /* printf("No element pointer\n"); */
		return;
	}

	nvert=pe->vertex;
	pn=RetrieveByNodeNumber(H,pe->index[nvert-1]);
	QMove(pn->coord.x,pn->coord.y);

	for(i=0;i<nvert;i++) {
		pn=RetrieveByNodeNumber(H,pe->index[i]);
		QPlot(pn->coord.x,pn->coord.y);
	}
}

void FillElem( Hashtable_type H , Elem_type *p , int color )

{
	Node_type *pn;
	int node,nvert,tmpcol,i;
	static float *x=NULL,*y=NULL;
	static int ndim=0;

	if( !p )
		return;

	nvert = p->vertex;
	if( ndim < nvert ) {
		if( ndim != 0 ) {
			free(x);
			free(y);
		}
		x = (float *) malloc( nvert*sizeof(float) );
		y = (float *) malloc( nvert*sizeof(float) );
		if( !x || !y ) Error("FillElem : Cannot allocate x/y");
		ndim=nvert;
	}

	for(i=0;i<nvert;i++) {
		node=p->index[i];
		pn=RetrieveByNodeNumber(H,node);
		x[i]=pn->coord.x; y[i]=pn->coord.y;
	}

	QGetPen(&tmpcol);
	QNewPen(color);
	QAreaFill(nvert,x,y);
	QNewPen(tmpcol);
}

/**************************************************************************/

void PlotLines( Hashtable_type HL , Hashtable_type HN )

{
	Line_type *p;

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	ResetHashTable(HL);
	while( (p = VisitHashTableL(HL)) != NULL ) {
		PlotLine(HN,p);
	}

	MouseShow();
}

void PlotLine( Hashtable_type H , Line_type *p )

{
	PlotLineP( H , p );
}

static void SetLineColor( Line_type *pl )

{
	int color;

	if( GetActLine() == pl->number ) {
	  QGetPen(&color);
	  if( color == EvidenceCol ) return;
	}

	QNewPen( GetTypeColor(pl->type) );
}

static void PlotLineP( Hashtable_type H , Line_type *pl )

{
	Node_type *pn;
	int i,nvert;

	if( !pl ) {   /* printf("No element pointer\n"); */
		return;
	}

	if( OpFill || OpColor ) {
	  SetLineColor(pl);
	}

	nvert=pl->vertex;
	pn=RetrieveByNodeNumber(H,pl->index[0]);
	QMove(pn->coord.x,pn->coord.y);

	for(i=1;i<nvert;i++) {
		pn=RetrieveByNodeNumber(H,pl->index[i]);
		QPlot(pn->coord.x,pn->coord.y);
	}
}

static void PlotLineSegP( Hashtable_type H , Line_type *pl , int node )

{
	Node_type *pn;
	int i,nvert;

	if( !pl )
		return;

	nvert=pl->vertex;
	for(i=0;i<nvert;i++)
		if( pl->index[i] == node ) break;

	if( i > 0 ) {
		pn=RetrieveByNodeNumber(H,pl->index[i]);
		QMove(pn->coord.x,pn->coord.y);
		pn=RetrieveByNodeNumber(H,pl->index[i-1]);
		QPlot(pn->coord.x,pn->coord.y);
	}

	if( i < nvert-1 ) {
		pn=RetrieveByNodeNumber(H,pl->index[i]);
		QMove(pn->coord.x,pn->coord.y);
		pn=RetrieveByNodeNumber(H,pl->index[i+1]);
		QPlot(pn->coord.x,pn->coord.y);
	}
}

void PlotSegment( float x1 , float y1 , float x2 , float y2 , int color )

{
        int savecolor;

        QGetPen(&savecolor);
        QNewPen(color);

        MouseHide();

        QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

        QLine(x1,y1,x2,y2);

        QNewPen(savecolor);

        MouseShow();
}

/**************************************************************************/

void PlotVectors( Hashtable_type HV )

{
	Node_type *p;

	MouseHide();

	QViewport(XPMin,YPMin,XPMax,YPMax);
	QWindow(GbPlo.low.x,GbPlo.low.y,GbPlo.high.x,GbPlo.high.y);

	ResetHashTable(HV);
	while( (p = VisitHashTableN(HV)) != NULL ){
		PlotVect(p);
	}

	MouseShow();
}

void PlotVect( Node_type *p )

{
	Vect_type *pv;
	float x,y;
	float s,d;
	float u,v;
	static float rad=3.14159/180.;

	x = p->coord.x;
	y = p->coord.y;

	pv = (Vect_type *) p->extra;
	s = pv->speed[pv->actual];
	d = pv->dir[pv->actual];

	s = OpVectFact*s;
	d = d > 180. ? d - 180. : d + 180.;
	d = d > 90. ? 450.-d : 90.-d;

	u=s*cos(rad*d); 
	v=s*sin(rad*d);

	PfeilUV(x,y,u,v,s);
}

/**************************************************************************/

static void PfeilUV( float x, float y, float u, float v, float uv )

{
	static float cos=0.866025;
	static float sin=0.5;
	float xu,yv,x1,y1,x2,y2;
	static float fact=0.3;

	if(uv==0.) return;

	xu=x+u;
	yv=y+v;

	x1=xu-fact*(cos*u+sin*v);
	y1=yv-fact*(-sin*u+cos*v);
	x2=xu-fact*(cos*u-sin*v);
	y2=yv-fact*(sin*u+cos*v);

	QMove(x,y);
	QPlot(xu,yv);
	QPlot(x1,y1);
	QPlot(x2,y2);
	QPlot(xu,yv);
}

/**************************************************************************/
