
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
 * grid_df.h - defines for grid                                         *
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

