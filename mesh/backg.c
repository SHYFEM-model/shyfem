
/************************************************************************\ 
 *									*
 * backg.c - background grid interpolation routines                     *
 *									*
 * Copyright (c) 1997 by Georg Umgiesser				*
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
 * 17-Nov-97: routines written from scratch				*
 *									*
\************************************************************************/


#include <stdlib.h>
/*
#include <stdio.h>
#include <math.h>
*/


#include "general.h"
#include "grd.h"
#include "grdhs.h"
#include "backg.h"


Backg_type *SetBackGrid( Grid_type *G , Elem_type *pe )

{
	int i,in,ip;
	float f=0.;
	float x[3],y[3];
	Node_type *pn;
	Backg_type *pb;

	pb = (Backg_type *) malloc( sizeof(Backg_type) );
	if( !pb )
		Error("SetBackGrid: Cannot allocate background mesh");

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(G->HN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
	}

	for(i=0;i<3;i++) {
		in = (i+1)%3;
		ip = (i+2)%3;
		pb->a[i] = x[in]*y[ip] - x[ip]*y[in];
		pb->b[i] = y[in] - y[ip];
		pb->c[i] = x[ip] - x[in];
		f += pb->a[i];
	}
	pb->fr = 1./f;

	return pb;
}

float InterpolBackGrid( Backg_type *pb , float x , float y , float *value )

{
	int i;
	float f;
	float a = 0.;

	for(i=0;i<3;i++) {
		f = pb->a[i] + x*pb->b[i] + y*pb->c[i];
		a += value[i] * f * pb->fr;
	}

	return a;
}

