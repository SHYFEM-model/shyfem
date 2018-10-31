
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
 * fund.h - fundamental type declarations                               *
 *									*
 * Revision History:							*
 * 02-Dec-95: Number_list_type cancelled from fundamental types         *
 * 12-Apr-94: Extracted from grid.h -> new header file                  *
 *									*
\************************************************************************/


#ifndef __GUH_FUND_
#define __GUH_FUND_


typedef unsigned char Small;

typedef struct {
	float x;
	float y;
} Point;

typedef struct {
	Point low;
	Point high;
} Rect;


#endif

