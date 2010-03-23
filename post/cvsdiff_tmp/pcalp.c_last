
/************************************************************************\ 
 *									*
 * pcalp.c - Calcomp emulation under Fortran (deprecated)		*
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
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUC_PCALP_
#define __GUC_PCALP_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pgraph.h"

#define XTrans(x)		( (x) * Factor + XOrigin )
#define YTrans(y)		( (y) * Factor + YOrigin )

static int PlotOpened=0;
static int MaxColors=16;

static float Factor=1.;
static float XAct=0.;
static float YAct=0.;
static float XOrigin=0.;
static float YOrigin=0.;

static char *String=NULL;
static int  NString=0;

/*
 *    newpen not active --> think about something else
 */

/*****************************************************************/

void plots_ ( long int *idum1 , long int *idum2 , long int *ldev )

{
	if( PlotOpened++ )
	    QClearScreen();
	else
	    QGraphInit();

	QNewPen(MaxColors-1);
	Factor = 1.;
	XAct = 0.;
	YAct = 0.;
	XOrigin = 0.;
	YOrigin = 0.;
}

void factor_ ( float *fct )

{
	Factor = *fct;
}

void where_ ( float *xpag , float *ypag , float *fct )

{
	*xpag = XAct;
	*ypag = YAct;
	*fct  = Factor;
}

void newpen_ ( long int *inpn )

{
	if( *inpn == 0 )
		QNewPen( MaxColors - 1 );
/*	else if ( *inpn > 0 && *inpn <= MaxColors)
		QNewPen( *inpn - 1 );
*/
}

void plot_ ( float *xpag , float *ypag , long int *ipen )

{
	long int abspen;
	float x,y;

	abspen = (*ipen>0) ? *ipen : -(*ipen);

	x = XTrans( *xpag );
	y = YTrans( *ypag );
	
	if( abspen == 999 )
		QGraphClose();
	else if( abspen == 2 )
		QPlot(x,y);
	else if( abspen == 3 )
		QMove(x,y);

	XAct = *xpag;
	YAct = *ypag;

	/* the next should be equal to XOrigin = XTrans(XAct), ... */

	if( *ipen < 0 ) {
		XOrigin += XAct * Factor;
		YOrigin += YAct * Factor;
	}
}

void symbol_ ( float *xpage , float *ypage , float *height 
		, char *ibcd , float *angle , long int *nchar )

	/* angle not yet used and no special symbols !!!! */

{
	float x,y,h;	

	if( *nchar < 0 ) return;

	if( *xpage == 999.0 || *ypage == 999.0 ) {
		x = XTrans( XAct );
		y = YTrans( YAct );
	} else {
		x = XTrans( *xpage );
		y = YTrans( *ypage );
	}
	h = *height * Factor;

/*	QTextPointSize( (int) (h*28.3) );*/	/* give height in points */
	QTextRealSize( h );			/* give height in cm */

	if( *nchar > NString ) {
		if( String ) free( String );
		String = (char *) malloc( *nchar + 1 );
		NString = *nchar;
	}
	String=strncpy(String,ibcd,*nchar);	
	String[*nchar]='\0';

	QText(x,y,String);

	/* change next lines if angle is used */

	XAct += *nchar * h;
	YAct += 0.;
}

void shade_ ( float *x , float *y , float *dx , float *dy , long int *color )

{
	float x1,x2,y1,y2;
	int col;

	x1 = XTrans( *x );
	y1 = YTrans( *y );
	x2 = XTrans( *x + *dx );
	y2 = YTrans( *y + *dy );

	col = (*color<=0) ? MaxColors-1 : *color-1;

	QRectFill( x1,y1,x2,y2,col );
}

#endif
