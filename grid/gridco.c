
/************************************************************************\
 *									*
 * gridco.c - color management routines                                 *
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
 * 01-Oct-2004: if OpMaxColDepth not given use max depth                *
 * 20-Nov-2003: changes for color table 0                               *
 * 17-Dec-97: new routine GetColorTabSize(), GetTypeColor()             *
 * 13-Oct-97: GetColor() returns PlotCol if value == NULLDEPTH          *
 * 10-Oct-97: New colortable 7 -> use HSB colors                        *
 * 25-Mar-95: NULLDEPTH is colored in grey (SetColors())                *
 * 16-Feb-95: ColTab and ColTabVel defined locally                      *
 * 10-Feb-95: routines copied from gridop.c                             *
 *									*
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general.h"
#include "grid.h"
#include "graph.h"
#include "gustd.h"
#include "color.h"


#define	 LIN(ncol,vmax,i) ( (vmax) * (i) / ( (ncol)-1 ) )
#define	QUAD(ncol,vmax,i) ( (vmax) * (i)*(i) / ( ((ncol)-1)*((ncol)-1) ) )


typedef struct {
	float value;
	int color;
} ColTab_list;

typedef struct {
	int size;
	ColTab_list *list;
} ColTab_type;

static ColTab_type *MakeColTab( int size );
static int GetColor( ColTab_type *C , float value );

static ColTab_type *ColTab;
static ColTab_type *ColTabVel;


int GetColorTabSize( void ) 

{
	return ColTab->size;
}

void SetColors( int coltab )

{
	int i;
	int size,ncol;
	float vmax;
	ColTab_list *c;
	int *colorp=NULL;

	size = OpColTabSize;
	vmax = OpMaxColDepth;
	if( vmax <= 0. ) {
	  vmax = GetDepthMinMax();
	}
	printf("Maximum scaling depth: %f\n",vmax);

	TotCol          =       Green;
	MenCol          =       Brown;
	ComCol          =       Black;
	ButCol          =       Red;
	MesCol          =       Black;
	ButChoosCol     =       LightRed;
	PlotCol         =       LightGray;
	EvidenceCol	=       LightRed;
	PlotWinCol      =       Black;
	BorderLightCol  =       LightGray;
	BorderDarkCol   =       DarkGray;
	ComWriteCol     =       White;
	MesWriteCol     =       White;
	ButWriteCol     =       White;

	if( coltab == 0) {
	    if( OpShowType == 0 ) {
		PlotCol         = LightGray;
		PlotWinCol      = Black;
		/* size=5; */
		size=1;
		ColTab = MakeColTab(size);
		c = ColTab->list;
		c[0].value =   0.; c[0].color = Cyan;
		/*
		c[0].value =   0.; c[0].color = Brown;
		c[1].value =  10.; c[1].color = Cyan;
		c[2].value =  50.; c[2].color = DarkGray;
		c[3].value = 100.; c[3].color = LightGreen;
		c[4].value = 500.; c[4].color = Blue;
		*/
	    } else {
		PlotCol         = LightGray;
		PlotWinCol      = Black;
		size            = 13;
		ColTab = MakeColTab(size);
		c = ColTab->list;
		c[0].color = Blue;
		c[1].color = Green;
		c[2].color = Cyan;
		c[3].color = Red;
		c[4].color = Magenta;
		c[5].color = Brown;
		c[6].color = Yellow;
		c[7].color = DarkGray;
		c[8].color = LightBlue;
		c[9].color = LightGreen;
		c[10].color = LightCyan;
		c[11].color = LightRed;
		c[12].color = LightMagenta;
		/* not used black,white,lightgray */
	    }
	} else if( coltab == 1) {
		PlotCol         = Black;
		PlotWinCol      = White;
		size            = 3;
		ColTab = MakeColTab(size);
		c = ColTab->list;
		c[0].value = 0.; c[0].color = Brown;
		c[1].value = 5.; c[1].color = Cyan;
		c[2].value = 10000.; c[2].color = Blue;
	} else if( coltab >= 2 &&  coltab <= 7 ) {

		PlotCol         = Black;
		PlotWinCol      = White;
		ncol		= size-2;

		if( coltab == 2 )
		    colorp = QAllocBlueColors(ncol);
		else if( coltab == 3 )
		    colorp = QAllocRed2YellowColors(ncol);
		else if( coltab == 4 )
		    colorp = QAllocRed2BlueColors(ncol);
		else if( coltab == 5 )
		    colorp = QAllocYellow2GreenColors(ncol);
		else if( coltab == 6 )
		    colorp = QAllocGreen2BlueColors(ncol);
		else if( coltab == 7 )
		    colorp = QAllocHSBColors(ncol);

		ColTab = MakeColTab(size);
		c = ColTab->list;

                c[0].value = -900.; c[0].color = LightGray;
		c[1].value = 0.; c[1].color = Brown;
		for(i=2;i<size;i++) {
			if( size > 30 )
				c[i].value = LIN(ncol,vmax,i);
			else
				c[i].value = QUAD(ncol,vmax,i);
			c[i].color = colorp[i-1];
#ifdef LOCAL_DEBUG
			printf("%d %d %f\n",i,c[i].color,c[i].value); 
#endif
		}

	} else {
                Error2("SetColors : No color table ",itos(coltab));
	}

	ColTabVel = ColTab;

}

int GetRandomColor( int value )

{
	value %= ColTab->size;
	return ColTab->list[value].color;
}

int GetTypeColor( int type )

{
	int color;

	color = type % ColTab->size;
	if( color == PlotWinCol ) color = EvidenceCol;

	return color;
}

int GetDepthColor( float value )
	{ return GetColor( ColTab , value ); }

int GetVelColor( float value )
	{ return GetColor( ColTabVel , value ); }

static int GetColor( ColTab_type *C , float value )

{
	int size,i;

	if( value == NULLDEPTH ) return PlotCol;

	size = C->size;
	for( i=0 ; i<size ; i++)
		if( value <= C->list[i].value )
			return C->list[i].color;
	return C->list[size-1].color;
}

static ColTab_type *MakeColTab( int size )

{
	ColTab_type *p;

	p = (ColTab_type *) malloc( sizeof(ColTab_type) );
	if( !p ) Error("MakeColTab : Cannot allocate ColTab");

	p->list = (ColTab_list *) malloc( size * sizeof(ColTab_list) );
	if( !p->list ) Error("MakeColTab : Cannot allocate ColTab list");

	p->size = size;

	return p;
}

#if OLD_VEL_COL

		size   = 29;				/* 17 */
		ColTabVel = MakeColTab(size);
		c = ColTabVel->list;
		c[0].value = 0.;
		c[0].color = Brown;
		for(i=1,val=0.;i<size;i++) {
			if(i<=20) 			/* 10 */
				dval=0.05*0.05;
			else if(i<=28) 			/* 15 */
				dval=0.1*0.1;
			else {
				val=10000.; dval=0.;
			}
			val += dval;
			c[i].value = val;
			c[i].color = VelShadeColor(i-1);
		}

#endif
