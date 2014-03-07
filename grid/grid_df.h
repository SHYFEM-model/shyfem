
/************************************************************************\ 
 *									*
 * grid_df.h - defines for grid                                         *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see grid.c for copying information					*
 *									*
 * Revision History:							*
 * 04-Dec-95: VECTMODE added                                            *
 * 04-May-94: LINEMODE added                                            *
 * 13-Apr-94: splitted up in df, ty, fp, ex                             *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRID_DF_
#define __GUH_GRID_DF_


#define NULLDEPTH -999.
#define NULLCOORD -999.

#define NOMODE   0
#define NODEMODE 1
#define ELEMMODE 2
#define LINEMODE 3
#define VECTMODE 4
#define SEGMMODE 5

#define AREAMIN 0.4

#define NODELIST_DIM 100

#define InPlotField(x,y) ((x) >= XPMin && (x) <= XPMax \
			    && (y) >= YPMin && (y) <= YPMax )
#define InMenuField(x,y) ((x) >= XMMin && (x) <= XMMax \
			    && (y) >= YMMin && (y) <= YMMax )

#define TRUE  1
#define FALSE 0

/* definitions for plot window */

#define XPMIN 30
#define YPMIN 30
#define XPMAX 450
#define YPMAX 450

/* definitions for menu window */

#define XMMIN 479
#define YMMIN 5
#define XMMAX 610
#define YMMAX 475

/* definitions for message window */

#define XSMIN 30
#define YSMIN 455
#define XSMAX 450
#define YSMAX 475

/* definitions for command window */

#define XCMIN 30
#define YCMIN 5
#define XCMAX 450
#define YCMAX 25

/* definitions for total window */

#define XTMIN 0
#define YTMIN 0
#define XTMAX 639
#define YTMAX 479


#endif

