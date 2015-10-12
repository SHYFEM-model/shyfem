
/************************************************************************\ 
 *									*
 * gridma2.c - routines used directly by main                           *
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
 * 07-Oct-2015: new routine MakeNewCenter()				*
 * 21-Feb-2014: new routine MakeGravityPoint() and InPoly()		*
 * 18-Feb-2014: new routines DelNodeLine(), InsertNodeLine()		*
 * 19-Nov-2003: write node info to stdout                               *
 * 02-Apr-1998: new menu routines integrated -> DisplayMainMenu()       *
 *                no MakeButtons()                                      *
 *                set total/menu window dimensions in ResizeWindow()    *
 * 02-Apr-1998: call QNewPen(PlotCol) befor plotting in PlotAll()       *
 * 09-Feb-1998: new algorithms to find node/element/line implemented in *
 *                FindClosestLine(), FindClosestNode(),                 *
 *                FindElemToPoint() -> finds more than one occurence    *
 * 14-Oct-97: Administer use of nodes with routines                     *
 *            SplitLine(): do not split at end of line                  *
 * 13-Oct-97: new routines AddUseN(), DeleteUseN(), GetUseN()           *
 *            adjourn use of nodes in SplitLine(), JoinLine()           *
 *            new routine DeleteLineWithNodes()                         *
 *            new routine DeleteElemWithNodes()                         *
 *            in routine MakeDepthFromNodes() return NULLDEPTH          *
 *              if one of the nodes has no depth                        *
 * 10-Oct-97: new routine SubstituteNode()                              *
 * 06-Dec-95: write to message window for vector                        *
 * 05-Dec-95: MakeVectActive EvidenceNone introduced                    *
 * 02-Dec-95: AdjustBiLinear, UndoBiLinear removed                      *
 *            PlotPlot() renamed to PlotAll()                           *
 * 16-May-94: new file for routines called by main                      *
 * 13-May-94: MakeDepthFromNodes() to compute depth                     *
 * 10-May-94: WriteToMessWin has changed format                         *
 * 05-May-94: ActXYValid removed (useless, use TentativeXY)             *
 * 13-Apr-94: completely restructured                                   *
 *             -> new hash.c for hash routines                          *
 *             -> new list.c for list table routines                    *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "grid.h"
#include "graph.h"

#include "general.h"
#include "gustd.h"
#include "keybd.h"
#include "mouse.h"

#include "list.h"
#include "hash.h"
#include "gridhs.h"
#include "gridky.h"
#include "gridnl.h"
#include "menu.h"


void PlotAll( void )

{
		QNewPen(PlotCol);
		PlotElements( HEL , HNN );
		PlotLines( HLI , HNN );
		PlotPoints( HNN );
		PlotVectors( HVC );
}

void RedrawAll( void )

{
	MouseHide();

	MakeTotalWindow( &GbTot );

	MakeMenuWindow( &GbMen );
	MakeMessageWindow( &GbMes );
	MakeCommandWindow( &GbCom );
	MakePlotWindow( &GbPlo );

	/* printf("RedrawAll: DisplayMainMenu()\n"); */
	DisplayMainMenu();
/*
	MakeButtons( &GbMen );

	PlotButtons();
*/

	PlotAll();

	MakeNodeActive(ActNode);
	MakeElemActive(ActElem);
	MakeLineActive(ActLine);
	MakeVectActive(ActVect);

	WriteToMesWindow();
	WriteToComWindow();

	MouseShow();
}

int ResizeWindow( int width , int height )

{
	int mw,scw,bx,by;

	SetWindowDimension(XTMin,YTMin,XTMax,YTMax);	/* FIXME */
	SetMenuDimension(XMMin,YMMin,XMMax,YMMax);

	if( width == XTMax-XTMin+1 && height == YTMax-YTMin+1 ) return 0;

	XTMax = width-1;
	YTMax = height-1;

	mw = 130;		/* width of menu */
	scw = 20;		/* height of message & command */
	bx = 30; by = 5;	/* border in x & y */

	XMMin = XTMax - bx - mw;
	XMMax = XMMin + mw - 1;

	XPMin = XSMin = XCMin = bx;
	XPMax = XSMax = XCMax = XMMin - bx;

	YMMin = YCMin = by;
	YMMax = YSMax = YTMax - by;

	YCMax = YCMin + scw;
	YPMin = YCMax + by;
	YSMin = YSMax - scw;
	YPMax = YSMin - by;

	SetWindowDimension(XTMin,YTMin,XTMax,YTMax);
	SetMenuDimension(XMMin,YMMin,XMMax,YMMax);

	return 1;
}


void WriteToMesWindow( void )

{
	int tmpcol;
	float xx,yy;
	float width,height;
	Elem_type *pe;
	Node_type *pn;
	Line_type *pl;
	Vect_type *pv;
	char s[20];

	QViewport(XSMin,YSMin,XSMax,YSMax);
	QWindow(GbMes.low.x,GbMes.low.y,GbMes.high.x,GbMes.high.y);
	QClearViewport();
	QRectFill(GbMes.low.x,GbMes.low.y,GbMes.high.x
			,GbMes.high.y,MesCol);
	PlotShadeRect( &GbMes , BorderDarkCol , BorderLightCol );

	width  = GbMes.high.x - GbMes.low.x;
	height = GbMes.high.y - GbMes.low.y;

	QGetPen(&tmpcol);
	QNewPen(MesWriteCol);

	if( ActMode == KEYBOARD_INPUT ) {  /* write keybd input to mes */

		xx = GbMes.low.x + 0.1 * width;
		yy = GbMes.low.y + 0.2 * height;
            	QText(xx,yy,GetKeyboardString());

	} else if( ActElem > 0 ) { /* write elem info to mes window */

		pe = RetrieveByElemNumber(HEL,ActElem);

		xx = GbMes.low.x + 0.05 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pe->number);
            	QText(xx,yy,s);

		xx = GbMes.low.x + 0.25 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pe->type);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.45 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pe->vertex);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.65 * width;
		yy = GbMes.low.y + 0.2 * height;
		if( pe->depth != NULLDEPTH ) {
		    sprintf(s,"%f",pe->depth);
		    QText(xx,yy,s);
		}

	} else if( ActLine > 0 ) { /* write line info to mes window */

		pl = RetrieveByLineNumber(HLI,ActLine);

		xx = GbMes.low.x + 0.05 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pl->number);
            	QText(xx,yy,s);

		xx = GbMes.low.x + 0.25 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pl->type);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.45 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pl->vertex);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.65 * width;
		yy = GbMes.low.y + 0.2 * height;
		if( pl->depth != NULLDEPTH ) {
		    sprintf(s,"%f",pl->depth);
		    QText(xx,yy,s);
		}

	} else if( ActVect > 0 ) { /* write vect info to mes window */

		pn = RetrieveByNodeNumber(HVC,ActVect);

		xx = GbMes.low.x + 0.05 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pn->number);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.25 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pn->type);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.4 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",pn->coord.x);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.6 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",pn->coord.y);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.8 * width;
		yy = GbMes.low.y + 0.2 * height;
		pv = pn->extra;
		sprintf(s,"%f",pv->speed[pv->actual]);
		QText(xx,yy,s);

	} else if( ActNode > 0 ) { /* write node info to mes window */

		/* how can we use also ActX, ActY ? */

		pn = RetrieveByNodeNumber(HNN,ActNode);

		xx = GbMes.low.x + 0.05 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pn->number);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.25 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%d",pn->type);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.4 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",pn->coord.x);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.6 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",pn->coord.y);
		QText(xx,yy,s);

		if( pn->depth != NULLDEPTH ) {
			xx = GbMes.low.x + 0.8 * width;
			yy = GbMes.low.y + 0.2 * height;
			sprintf(s,"%f",pn->depth);
			QText(xx,yy,s);
		}

		/* write also to stdout */

                printf("%d %d %f %f ",pn->number,pn->type
                               ,pn->coord.x,pn->coord.y);
                if( pn->depth != NULLDEPTH ) {
                        printf("%f ",pn->depth);
                }
                printf("\n");

	} else { /* just show coordinates */

		xx = GbMes.low.x + 0.20 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",ActX);
		QText(xx,yy,s);

		xx = GbMes.low.x + 0.40 * width;
		yy = GbMes.low.y + 0.2 * height;
		sprintf(s,"%f",ActY);
		QText(xx,yy,s);

	}

	QNewPen(tmpcol);
}

void WriteToComWindow( void )

{
	int tmpcol;
	float xx,yy;
	float width,height;

	QViewport(XCMin,YCMin,XCMax,YCMax);
	QWindow(GbCom.low.x,GbCom.low.y,GbCom.high.x,GbCom.high.y);
	QClearViewport();
	QRectFill(GbCom.low.x,GbCom.low.y,GbCom.high.x,GbCom.high.y,ComCol);
	PlotShadeRect( &GbCom , BorderDarkCol , BorderLightCol );

	width  = GbCom.high.x - GbCom.low.x;
	height = GbCom.high.y - GbCom.low.y;

	xx = GbCom.low.x + 0.1 * width;
	yy = GbCom.low.y + 0.2 * height;

	QGetPen(&tmpcol);
	QNewPen(ComWriteCol);
        QText(xx,yy,ActString);
	QNewPen(tmpcol);
}

#define SQUARE(x) ( (x)*(x) )
#define GETDIST(x1,y1,x2,y2) ( SQUARE((x2)-(x1)) + SQUARE((y2)-(y1)) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

#define LISTDIM	100

Line_type *FindClosestLine( Hashtable_type HL , Hashtable_type HN 
				, float x , float y )

/*\
 *  find closest midpoint to point
\*/

{
	int j,nvert;
	Line_type *pl;
	Node_type *pn;
	float dist,dd,dd1,dd2;
	float x1,y1,x2,y2;
	static Line_type *pold=NULL;
	static int iold = 0;
	Line_type *plist[LISTDIM];	/* should be enough */
	int i = 0;

	dist=GbPlo.high.x-GbPlo.low.x+GbPlo.high.y-GbPlo.low.y;

	dist=dist*dist/600.;    /* 600 is empiric -> not always find TP */

	ResetHashTable(HL);
	while( (pl=VisitHashTableL(HL)) != NULL )
	{
		pn = RetrieveByNodeNumber(HN,pl->index[0]);
		x1 = pn->coord.x;
		y1 = pn->coord.y;

		dd=GETDIST(x1,y1,x,y);
		if( dd <= dist ) {
		    if( dd < dist ) {
			i = 0;
			dist = dd;
		    }
		    plist[ i++ % LISTDIM ] = pl;
		}

		nvert=pl->vertex;
		for(j=1;j<nvert;j++) {
			pn = RetrieveByNodeNumber(HN,pl->index[j]);
			x2 = pn->coord.x;
			y2 = pn->coord.y;

			dd1=GETDIST(0.5*(x2+x1),0.5*(y2+y1),x,y);
			dd2=GETDIST(x2,y2,x,y);
			dd = MIN(dd1,dd2);
			if( dd <= dist ) {
			    if( dd < dist ) {
				i = 0;
				dist = dd;
			    }
		    	    plist[ i++ % LISTDIM ] = pl;
			}

			x1 = x2;
			y1 = y2;
		}
	}

	if( i >= LISTDIM ) i = LISTDIM;

	if( i == 0 ) {
		return NULL;
	} else if( i == 1 ) {
		return plist[0];
	} else { /* more lines found */
		if( plist[0] != pold ) {	/* first time */
			iold = 0;
			pold = plist[0];
		}
		iold++;
		return plist[ iold%i];
	}
}

Node_type *FindClosestNode( Hashtable_type H , float x , float y )

{
	Node_type *pn;
	float dist,dx,dy,dd;
	static Node_type *pold=NULL;
	static int iold = 0;
	Node_type *plist[LISTDIM];	/* should be enough */
	int i = 0;

	dist=GbPlo.high.x-GbPlo.low.x+GbPlo.high.y-GbPlo.low.y;

	dist=dist*dist/600.;    /* 600 is empiric -> not always find TP */

	ResetHashTable(H);
	while( (pn=VisitHashTableN(H)) != NULL )
	{
		dx=pn->coord.x-x;
		dy=pn->coord.y-y;
		dd=dx*dx+dy*dy;
		if( dd <= dist ) {
		    if( dd < dist ) {
			i = 0;
			dist = dd;
		    }
		    plist[ i++ % LISTDIM ] = pn;
		}
	}

	if( i >= LISTDIM ) i = LISTDIM;

	if( i == 0 ) {
		return NULL;
	} else if( i == 1 ) {
		return plist[0];
	} else { /* more nodes found */
		if( plist[0] != pold ) {	/* first time */
			iold = 0;
			pold = plist[0];
		}
		iold++;
		return plist[ iold%i];
	}
}

Elem_type *FindElemToNode( Hashtable_type H , Elem_type *p , int node )

{
	int i,nvert;

	if( !p ) ResetHashTable(H);

	while( (p=VisitHashTableE(H)) != NULL ) {
		nvert=p->vertex;
		for(i=0;i<nvert;i++)
			if( p->index[i] == node ) break;
		if( i < nvert ) break;      /* node found */
	}

	return p;
}

Line_type *FindLineToNode( Hashtable_type H , Line_type *p , int node )

{
	int i,nvert;

	if( !p ) ResetHashTable(H);

	while( (p=VisitHashTableL(H)) != NULL ) {
		nvert=p->vertex;
		for(i=0;i<nvert;i++)
			if( p->index[i] == node ) break;
		if( i < nvert ) break;      /* node found */
	}

	return p;
}

Node_type *FindNode( Hashtable_type H , int node )

{
	return RetrieveByNodeNumber( H , node );
}

Elem_type *FindElem( Hashtable_type H , int elem )

{
	return RetrieveByElemNumber( H , elem );
}

Line_type *FindLine( Hashtable_type H , int line )

{
	return RetrieveByLineNumber( H , line );
}

Elem_type *FindElemToPoint( Hashtable_type HE , Hashtable_type HN
				, float x , float y )

{
	Elem_type *pe=NULL;
	Elem_type *plist[LISTDIM];	/* should be enough */
	static Elem_type *pold=NULL;
	static int iold = 0;
	int i = 0;

	ResetHashTable(HE);
	while( (pe=VisitHashTableE(HE)) != NULL )
	{
		if( InElement( HN , pe , x , y ) ) {
		    plist[ i++ % LISTDIM ] = pe;
		}
	}

	if( i >= LISTDIM ) i = LISTDIM;

	if( i == 0 ) {
		return NULL;
	} else if( i == 1 ) {
		return plist[0];
	} else { /* more lines found */
		if( plist[0] != pold ) {	/* first time */
			iold = 0;
			pold = plist[0];
		}
		iold++;
		return plist[iold%i];
	}
}

int InConvex( int n , float *xe , float *ye , float x , float y )

{
        int i,ii;
	float xmin,xmax,ymin,ymax;
        double scal,sc;

	if( n < 3 ) return 0;

	xmin=xe[0]; xmax=xmin; ymin=ye[0]; ymax=ymin;
	for( i=1 ; i<n ; i++ ) {
	  if(xe[i]<xmin) xmin=xe[i];
	  if(xe[i]>xmax) xmax=xe[i];
	  if(ye[i]<ymin) ymin=ye[i];
	  if(ye[i]>ymax) ymax=ye[i];
	}

        if( xmin>x || xmax<x || ymin>y || ymax<y ) return 0;

	ii=n-1;
	scal=(x-xe[ii])*(ye[0]-ye[ii])-(y-ye[ii])*(xe[0]-xe[ii]);
        for( i=1 ; i<n ; i++ ) {
		ii=i-1;
		sc=(x-xe[ii])*(ye[i]-ye[ii])-(y-ye[ii])*(xe[i]-xe[ii]);
		if( sc*scal < 0. ) return 0;
        }

        return 1;
}

int InTriangle( Hashtable_type H , Elem_type *pe , float x , float y );
int InPoly( Hashtable_type H , Elem_type *pe , float x , float y );

int InElement( Hashtable_type H , Elem_type *pe , float x , float y )

{
	Rect r;
	int nvert;

	nvert = pe->vertex;

	if( nvert < 3 ) return 0;

	PolyMinMax( H , pe , &r );
	if( r.low.x>x || r.high.x<x || r.low.y>y || r.high.y<y ) return 0;

	if( nvert == 3 ) {
	  return InTriangle(H,pe,x,y);
	} else {
	  return InPoly(H,pe,x,y);
	}
}

int InTriangle( Hashtable_type H , Elem_type *pe , float x , float y )

{
	Node_type *pn;
	int i,nvert;
	float scal,scao;
	Point *cn,*co;

	nvert = pe->vertex;

	scao=0.;
	pn = RetrieveByNodeNumber(H,pe->index[nvert-1]);
	co = &pn->coord;
	for(i=0;i<nvert;i++) {
		pn = RetrieveByNodeNumber(H,pe->index[i]);
		cn = &pn->coord;
		scal = (x - co->x)*(cn->y - co->y)-(y - co->y)*(cn->x - co->x);
		if( scao*scal < 0. ) return 0;
		scao = scal;
		co = cn;
	}
	return 1;
}

static float isLeft(float x0,float y0,float x1,float y1,float x2,float y2)
{
    return ( (x1 - x0) * (y2 - y0)
            - (x2 -  x0) * (y1 - y0) );
}

int InPoly( Hashtable_type H , Elem_type *pe , float x , float y )

{
/*
	http://geomalgorithms.com/a03-_inclusion.html
        wn_PnPoly(): winding number test for a point in a polygon
        Input:   P = a point,
                 V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
        Return:  wn = the winding number (=0 only when P is outside)
*/

	Node_type *pn;
	int i,nvert;
	Point *cn,*co;
        int wn = 0;    // the  winding number counter

	nvert = pe->vertex;

	pn = RetrieveByNodeNumber(H,pe->index[nvert-1]);
	co = &pn->coord;

        for (i=0; i<nvert; i++) {     // edge from V[i] to  V[i+1]
	  pn = RetrieveByNodeNumber(H,pe->index[i]);
	  cn = &pn->coord;
          if (co->y <= y) {          // start y <= P.y
            if (cn->y  > y)     	    // an upward crossing
                 if (isLeft(co->x,co->y,cn->x,cn->y,x,y) > 0)// P left of  edge
                     ++wn;          // have  a valid up intersect
          } else {                  // start y > P.y (no test needed)
            if (cn->y  <= y)	    // a downward crossing
                 if (isLeft(co->x,co->y,cn->x,cn->y,x,y) < 0)// P right of  edge
                     --wn;          // have  a valid down intersect
          }
	  co = cn;
        }

        return wn;
}

int GetUseN( Node_type *pn )

{
	return pn->use;
}

void AddUseN( Node_type *pn )

{
	pn->use++;
}

void AddUseE( Hashtable_type H , Elem_type *pe )

{
	Node_type *pn;
	int i;

	for(i=0;i<pe->vertex;i++) {
		pn=RetrieveByNodeNumber(H,pe->index[i]);
		pn->use++;
	}
}

void AddUseL( Hashtable_type H , Line_type *pl )

{
	Node_type *pn;
	int i;

	for(i=0;i<pl->vertex;i++) {
		pn=RetrieveByNodeNumber(H,pl->index[i]);
		pn->use++;
	}
}

void DeleteUseN( Node_type *pn )

{
	pn->use--;
}

void DeleteUseE( Hashtable_type H , Elem_type *pe )

{
	Node_type *pn;
	int i;

	for(i=0;i<pe->vertex;i++) {
		pn=RetrieveByNodeNumber(H,pe->index[i]);
		pn->use--;
	}
}

void DeleteUseL( Hashtable_type H , Line_type *pl )

{
	Node_type *pn;
	int i;

	for(i=0;i<pl->vertex;i++) {
		pn=RetrieveByNodeNumber(H,pl->index[i]);
		pn->use--;
	}
}

void UnActive( void )

{
	MakeNodeActive(0);
	MakeElemActive(0);
	MakeLineActive(0);
}

void EvidenceNone( void )

{
	EvidenceNode(ActNode,PlotCol);
	EvidenceElem(ActElem,PlotCol);
	EvidenceLine(ActLine,PlotCol);
	EvidenceVect(ActVect,PlotCol);
	ActNode = 0;
	ActElem = 0;
	ActLine = 0;
	ActVect = 0;
}

void MakeNodeActive( int node )

{
	if( node ) {
		EvidenceNone();
	}
	EvidenceNode(ActNode,PlotCol);
	ActNode = node;
	EvidenceNode(ActNode,EvidenceCol);
}

void MakeLineActive( int line )

{
	if( line ) {
		EvidenceNone();
	}
	EvidenceLine(ActLine,PlotCol);
	ActLine = line;
	EvidenceLine(ActLine,EvidenceCol);
}

void MakeElemActive( int elem )

{
	if( elem ) {
		EvidenceNone();
	}
	EvidenceElem(ActElem,PlotCol);
	ActElem = elem;
	EvidenceElem(ActElem,EvidenceCol);
}

void MakeVectActive( int vect )

{
	if( vect ) {
		EvidenceNone();
	}
	EvidenceVect(ActVect,PlotCol);
	ActVect = vect;
	EvidenceVect(ActVect,EvidenceCol);
}

void MakeNewCenter(Rect *gp , float *x , float *y , float fact )

{
	float width,height;
	float rx,ry;		/* relative coordinates */
	float xl,yl;		/* new lower left coordinates */

	width  = ( gp->high.x - gp->low.x ) * fact;
	height = ( gp->high.y - gp->low.y ) * fact;

	rx = ( *x - gp->low.x ) / ( gp->high.x - gp->low.x );
	ry = ( *y - gp->low.y ) / ( gp->high.y - gp->low.y );

	xl = *x - rx * width;
	yl = *y - ry * height;

	*x = xl + width/2.;
	*y = yl + height/2.;
}

void ZoomInOut(Rect *gp , float x , float y , float fact )

{
	float width,height;

	width  = ( gp->high.x - gp->low.x ) * fact;
	height = ( gp->high.y - gp->low.y ) * fact;
	gp->low.x  = x - width/2.;
	gp->high.x = x + width/2.;
	gp->low.y  = y - height/2.;
	gp->high.y = y + height/2.;
}

void MoveRelative(Rect *gp , float dx , float dy )

{
	gp->low.x  += dx;
	gp->high.x += dx;
	gp->low.y  += dy;
	gp->high.y += dy;
}

void MoveToPoint( float x , float y )

{
	int horiz,verti;

	QScreenXY(x,y,&horiz,&verti);
	MouseMove(horiz,verti);
}

void MakeMidPoint( Line_type *pl , float *x , float *y )

{
	Node_type *pn,*paux;
	int nvert,i;
	float x1,y1,x2,y2;
	float dx,dy,dist,mindist;

	nvert=pl->vertex;

	if( nvert == 2 ) {	/* one segment -> mid point */
		pn=RetrieveByNodeNumber(HNN,pl->index[0]);
		x1 = pn->coord.x;
		y1 = pn->coord.y;
		pn=RetrieveByNodeNumber(HNN,pl->index[1]);
		x2 = pn->coord.x;
		y2 = pn->coord.y;
		*x = (x1+x2)/2.;
		*y = (y1+y2)/2.;
	} else {		/* more segments -> closest node */
		pn=RetrieveByNodeNumber(HNN,pl->index[0]);
		dx = pn->coord.x - *x;
		dy = pn->coord.y - *y;
		mindist = dx*dx + dy*dy;
		paux=pn;
		for(i=1;i<nvert;i++) {
			pn=RetrieveByNodeNumber(HNN,pl->index[i]);
			dx = pn->coord.x - *x;
			dy = pn->coord.y - *y;
			dist = dx*dx + dy*dy;
			if( dist < mindist ) {
				mindist=dist;
				paux=pn;
			}
		}
		*x = paux->coord.x;
		*y = paux->coord.y;
	}
}

void MakeCenterPoint( Elem_type *pe , float *xc , float *yc )

{
	Node_type *pn;
	int i,nvert;
	float xx=0.,yy=0.;

	nvert=pe->vertex;
	for(i=0;i<nvert;i++) {
		pn=RetrieveByNodeNumber(HNN,pe->index[i]);
		xx += pn->coord.x;
		yy += pn->coord.y;
	}

	*xc=xx/nvert;
	*yc=yy/nvert;
}

static float areat( Point *c1, Point *c2, Point *c3 )
{
	return 0.5*( (c2->x-c1->x) * (c3->y-c1->y) 
			- (c3->x-c1->x) * (c2->y-c1->y) );
}

void MakeGravityPoint( Elem_type *pe , float *xc , float *yc )

{
	Node_type *pn;
	Point *c0, *co, *cn;
	int i,nvert;
	double area;
	double areasum=0.,xx=0.,yy=0.;

	nvert = pe->vertex;

	pn = RetrieveByNodeNumber(HNN,pe->index[0]);
	c0 = &pn->coord;
	pn = RetrieveByNodeNumber(HNN,pe->index[1]);
	cn = &pn->coord;

	for(i=2;i<nvert;i++) {
		co = cn;
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		cn = &pn->coord;
		area = areat(c0,co,cn);
		areasum = areasum + area;
		xx += area*(c0->x+co->x+cn->x);
		yy += area*(c0->y+co->y+cn->y);
	}

	*xc=xx/(3.*areasum);
	*yc=yy/(3.*areasum);
}

/*****************************************************************/

void Working( void )

{
	char *actual;

	actual = ActString;
/*	ActString = StringWorking; */
	ActString = "Working ...(to change)";
	WriteToComWindow();
	ActString = actual;
}

/*****************************************************************/

void DeleteElemWithNodes( Elem_type *pe ) 

{
	int i;
	int nvert=pe->vertex;
	int *index=pe->index;
	Node_type *pn;

	/* use of node has already been decremented */

	for(i=0;i<nvert;i++) {
		pn = FindNode(HNN,index[i]);
		if( GetUseN(pn) == 0 ) DeleteNode(pn);
	}

	DeleteElem(pe);
}

void DeleteLineWithNodes( Line_type *pl ) 

{
	int i;
	int nvert=pl->vertex;
	int *index=pl->index;
	Node_type *pn;

	/* use of node has already been decremented */
	/* if line is closed, look out for already deleted node */

	for(i=0;i<nvert;i++) {
		pn = FindNode(HNN,index[i]);
		if( pn && GetUseN(pn) == 0 ) DeleteNode(pn);
	}

	DeleteLine(pl);
}

void SplitLine( Hashtable_type H , Line_type *pl , int node )

{
	int nvert,i,imid;
	int *index,*index1,*index2;
	Line_type *p;
	Node_type *pn= FindNode(HNN,node);

	if( pl == NULL ) {
                Warning("SplitLine : Null line pointer");
                return;
        }

	nvert=pl->vertex;
	index=pl->index;

	for(i=0;i<nvert;i++)
		if( index[i] == node ) break;

        if(i==nvert) {
                Warning("SplitLine : Node not found in line");
                return;
        }

	if( i == 0 || i == nvert-1 ) return; /* end of line -> do nothing */

	imid = i;
	index1 = MakeIndex(imid+1);
	index2 = MakeIndex(nvert-imid);

	for(i=0;i<=imid;i++)
		index1[i] = index[i];
	for(i=imid;i<nvert;i++)
		index2[i-imid] = index[i];

	NTotLines++;
	p = MakeLineWithIndex(NTotLines,pl->type,nvert-imid,index2);
	InsertByLineNumber(H,p);

	free(index);
	pl->vertex = imid+1;
	pl->index = index1;

	AddUseN( pn );
}

void JoinLine( Hashtable_type H , Line_type *p1 ,  Line_type *p2 , int node )

{
	int nvert1,nvert2,i;
	int *index,*index1,*index2;
	Node_type *pn= FindNode(HNN,node);

        if( p1 == p2 ) {
                Warning("JoinLine : Cannot join the same line");
                return;
	} else if( p1 == NULL || p2 == NULL ) {
                Warning("JoinLine : Null line pointer");
                return;
        }

	nvert1=p1->vertex;
	index1=p1->index;
	nvert2=p2->vertex;
	index2=p2->index;

	if( index1[nvert1-1] == node && index2[0] == node ) {
		; /* ok */
	} else if( index1[0] == node && index2[nvert2-1] == node ) {
		InvertIndex(index2,nvert2);
		InvertIndex(index1,nvert1);
	} else if( index1[nvert1-1] == node && index2[nvert2-1] == node ) {
		InvertIndex(index2,nvert2);
	} else if( index1[0] == node && index2[0] == node ) {
		InvertIndex(index1,nvert1);
	} else {
		Warning("JoinLine : Node must be a terminal node in lines");
		return;
	}

	index = MakeIndex(nvert1+nvert2-1);

	for(i=0;i<nvert1;i++)
		index[i] = index1[i];
	for(i=1;i<nvert2;i++)
		index[nvert1+i-1] = index2[i];

	free(index1);
	free(index2);
	p2=(Line_type *) DeleteHashByNumber(H,p2->number);
	free(p2);

	p1->vertex = nvert1+nvert2-1;
	p1->index = index;

	DeleteUseN( pn );
}

void DelNodeLine( Hashtable_type H , Line_type *pl , int node )

{
	int nvert,i,imid;
	int *index,*index1;
	Node_type *pn= FindNode(HNN,node);

	if( pl == NULL ) {
                Warning("SplitLine : Null line pointer");
                return;
        }

	nvert=pl->vertex;
	index=pl->index;

	for(i=0;i<nvert;i++)
		if( index[i] == node ) break;

        if(i==nvert) {
                Warning("SplitLine : Node not found in line");
                return;
        }

	imid = i;
	index1 = MakeIndex(nvert-1);

	for(i=0;i<imid;i++)
		index1[i] = index[i];
	for(i=imid+1;i<nvert;i++)
		index1[i-1] = index[i];

	free(index);
	pl->vertex = nvert-1;
	pl->index = index1;

	DeleteUseN( pn );
}

void InsertNodeLine( Hashtable_type H , Line_type *pl , int node )

{
	int nvert,i,imid;
	int *index,*index1;
	float d,dd,d1,d2,dmax;
	int imax;
	Node_type *pn= FindNode(HNN,node);

	if( pl == NULL ) {
                Warning("SplitLine : Null line pointer");
                return;
        }

	nvert=pl->vertex;
	index=pl->index;

	imax = 0;
	dmax = 2. * Dist2Node(HNN,node,index[0]);	/* initialize */
	printf( "InsertNodeLine: %f %d\n",dmax,imax);
	for(i=0;i<nvert;i++) {
	  d = Dist2Node(HNN,node,index[i]);
	  if( d < dmax ) {
	    dmax = d;
	    imax = i;
	  }
	}
	printf( "InsertNodeLine: %f %d\n",dmax,imax);

	if( imax == 0 ) {
	  dd = Dist2Node(HNN,index[0],index[1]);
	  if( Dist2Node(HNN,node,index[1]) < dd ) {
	    imid = 1;
	  } else {	/* insert before first node */
	    imid = 0;
	  }
	} else if( imax == nvert-1 ) {
	  dd = Dist2Node(HNN,index[nvert-2],index[nvert-1]);
	  if( Dist2Node(HNN,node,index[nvert-2]) < dd ) {
	    imid = nvert-1;
	  } else {	/* insert after last node */
	    imid = nvert;
	  }
	} else {
	  d1 = Dist2Node(HNN,node,index[imax-1])
		/ Dist2Node(HNN,index[imax],index[imax-1]);
	  d2 = Dist2Node(HNN,node,index[imax+1])
		/ Dist2Node(HNN,index[imax],index[imax+1]);
	  if( d1 < d2 ) {
	    imid = imax;
	  } else {
	    imid = imax+1;
	  }
	}

	/* imid is node before we have to insert */

	index1 = MakeIndex(nvert+1);

	for(i=0;i<imid;i++)
		index1[i] = index[i];
	index1[imid] = node;
	for(i=imid;i<nvert;i++)
		index1[i+1] = index[i];

	free(index);
	pl->vertex = nvert+1;
	pl->index = index1;

	AddUseN( pn );
}

/*******************************************************************/

float MakeDepthFromNodes( Hashtable_type H , Elem_type *p )

{
	int i,nvert;
	float depth;
	Node_type *pn;

	nvert=p->vertex;
	depth=0.;

        for(i=0;i<nvert;i++) {
                pn=RetrieveByNodeNumber(H,p->index[i]);
		if( pn->depth == NULLDEPTH ) return NULLDEPTH;
		depth += pn->depth;
        }

	return depth/nvert;
}

void SubstituteNode( Node_type *pna , Node_type *pnu )

{
	Elem_type *pe;
	Line_type *pl;
	int i;
	int nnew = pna->number;
	int nold = pnu->number;

        ResetHashTable(HEL);
        while( (pe=VisitHashTableE(HEL)) != NULL ) {
                for(i=0;i<pe->vertex;i++) {
			if( nold == pe->index[i] ) {
				pe->index[i] = nnew;
			}
		}
	}

        ResetHashTable(HLI);
        while( (pl=VisitHashTableL(HLI)) != NULL ) {
                for(i=0;i<pl->vertex;i++) {
			if( nold == pl->index[i] ) {
				pl->index[i] = nnew;
			}
		}
	}

	while( GetUseN(pnu) ) {
		DeleteUseN(pnu);
		AddUseN(pna);
	}
}

