
/************************************************************************\
 *
 *    Copyright (C) 1992,1995,1997,2004  Georg Umgiesser
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
 * gridps.c - read/write routines for Postscript file
 *
 * revision log :
 *
 * 21.10.1992	ggu	routines written from scratch
 * 12.01.1995	ggu	PsOpen local to file, ClosePs() new
 * 14.10.1997	ggu	Administer use of nodes with routines
 * 24.09.2004	ggu	new routines from post used (psgraph.c,psgraph.h)
 * 22.11.2004	ggu	bugfix for line plotting
 *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"

#include "grid.h"
#include "gustd.h"
#include "keybd.h"

#include "args.h"
#include "queue.h"
#include "list.h"
#include "hash.h"
#include "gridhs.h"

#include "psgraph.h"

static int PsOpen = FALSE;

/**************************************************************************/

static void crossPs( float x, float y, float size )

{
	PsLine(x-size,y,x+size,y);
	PsLine(x,y-size,x,y+size);
}

void ClosePS( void )

{
	if( PsOpen )
		PsGraphClose();
	PsOpen = FALSE;		/* new - maybe leave open ... */
}

void WritePS( Rect *gp , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , Queuetable_type C )

{
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;
	int i,j;
	int nvertex,node;
	float x,y,dxy;
	float xmin,xmax,ymin,ymax;
	Rect r;

	xmin=gp->low.x; xmax=gp->high.x;
	ymin=gp->low.y; ymax=gp->high.y;

	if( ! PsOpen )
		PsGraphOpen();

	PsStartPage();
	PsOpen = TRUE;

	PsSetWorld(xmin,ymin,xmax,ymax);
	PsRectifyScale();

	PsMove(xmin,ymin);
	PsPlot(xmax,ymin);
	PsPlot(xmax,ymax);
	PsPlot(xmin,ymax);
	PsPlot(xmin,ymin);

/* Plot Window is always a square -> use any direction to compute dxy */
	/* FIXME */

	dxy = (gp->high.x - gp->low.x)/190.;	/* should be 1mm */

        for(i=1;i<=NTotNodes;i++) {
          if( (pn=RetrieveByNodeNumber(HN,i)) != NULL ) {
		if( !GetUseN(pn) ) {
			x=pn->coord.x;
			y=pn->coord.y;
			if(xmin<=x && x<=xmax && ymin<=y && y<=ymax) {
				crossPs(x,y,dxy);
			}
		}
          }
        }

	PsFlush();

        for(i=1;i<=NTotElems;i++) {
          if( (pe=RetrieveByElemNumber(HE,i)) != NULL ) {
		nvertex=pe->vertex;
		PolyMinMaxIndex(HN,nvertex,pe->index,&r);
		if(r.low.x>xmax||r.high.x<xmin||r.low.y>ymax||r.high.y<ymin)
			continue;
		node=pe->index[nvertex-1];
		pn=RetrieveByNodeNumber(HN,node);
		PsMove(pn->coord.x,pn->coord.y);
		for(j=0;j<nvertex;j++) {
			node=pe->index[j];
			pn=RetrieveByNodeNumber(HN,node);
			PsPlot(pn->coord.x,pn->coord.y);
		}
          }
        }

	PsFlush();

        for(i=1;i<=NTotLines;i++) {
          if( (pl=RetrieveByLineNumber(HL,i)) != NULL ) {
		nvertex=pl->vertex;
		PolyMinMaxIndex(HN,nvertex,pl->index,&r);
		if(r.low.x>xmax||r.high.x<xmin||r.low.y>ymax||r.high.y<ymin)
			continue;
		/* node=pl->index[nvertex-1]; */
		node=pl->index[0];
		pn=RetrieveByNodeNumber(HN,node);
		PsMove(pn->coord.x,pn->coord.y);
                for(j=1;j<pl->vertex;j++) {
			node=pl->index[j];
			pn=RetrieveByNodeNumber(HN,node);
			PsPlot(pn->coord.x,pn->coord.y);
                }
          }
        }

	PsFlush();
	PsEndPage();
}

