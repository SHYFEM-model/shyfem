
/************************************************************************\
 *									*
 * gridps.c - read/write routines for Postscript file                   *
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
 * 22-Nov-2004: bugfix for line plotting                                *
 * 24-Sep-2004: new routines from post used (psgraph.c,psgraph.h)       *
 * 14-Oct-97: Administer use of nodes with routines                     *
 * 12-Jan-95: PsOpen local to file, ClosePs() new                       *
 * 21-Oct-92: routines written from scratch				*
 *									*
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

