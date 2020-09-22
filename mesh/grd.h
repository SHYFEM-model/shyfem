
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * grd.h - standard header for grd files
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines written from scratch
 *
\************************************************************************/



#ifndef __GUH_GRD_
#define __GUH_GRD_


#include "fund.h"
#include "queue.h"
#include "hash.h"


typedef struct {
	int totnode;
	int totelem;
	int totline;
	Hashtable_type HN;
	Hashtable_type HE;
	Hashtable_type HL;
	QueueTable C;
} Grid_type;

typedef struct {	/* identical structure of all tags */
        int number;
        int type;
	void *extra;
} Item_type;

typedef struct {
        int number;
        int type;
	void *extra;
        Point coord;
	float depth;
	int use;
} Node_type;

typedef struct {
        int number;
        int type;
	void *extra;
        int vertex;
        int *index;
	float depth;
} Elem_type;

typedef struct {
        int number;
        int type;
	void *extra;
        int vertex;
        int *index;
	float depth;
} Line_type;


#define NULLDEPTH -999.0		/* null depth */


#endif

