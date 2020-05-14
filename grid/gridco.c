
/************************************************************************\
 *
 *    Copyright (C) 1995,1997,2003-2004  Georg Umgiesser
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
 * gridco.c - color management routines
 *
 * revision log :
 *
 * 10.02.1995	ggu	routines copied from gridop.c
 * 16.02.1995	ggu	ColTab and ColTabVel defined locally
 * 25.03.1995	ggu	NULLDEPTH is colored in grey (SetColors())
 * 10.10.1997	ggu	New colortable 7 -> use HSB colors
 * 13.10.1997	ggu	GetColor() returns PlotCol if value == NULLDEPTH
 * 17.12.1997	ggu	new routine GetColorTabSize(), GetTypeColor()
 * 20.11.2003	ggu	changes for color table 0
 * 01.10.2004	ggu	if OpMaxColDepth not given use max depth
 *
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
