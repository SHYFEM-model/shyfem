
/************************************************************************\ 
 *									*
 * fund.h - fundamental type declarations                               *
 *									*
 * Copyright (c) 1994 by Georg Umgiesser				*
 *									*
 * see general.c for copying information				*
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

