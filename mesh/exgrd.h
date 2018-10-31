
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
 * exgrd.h - extracts items from grd file				*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/



#ifndef __GUH_EXGRD_
#define __GUH_EXGRD_


#include "fund.h"
#include "queue.h"
#include "hash.h"


typedef struct node_tag {
        Point coord;
	float depth;
        int type;
        int number;
} Node_type;

typedef struct elem_tag {
        int vertex;
        int *index;
	float depth;
        int type;
        int number;
} Elem_type;

typedef struct line_tag {
        int vertex;
        int *index;
        int type;
        int number;
} Line_type;

/* mark null depth */

#define NULLDEPTH -999.0

/*********************************/

extern int   NTotNodes;       /* number of total nodes */
extern int   NTotElems;       /* number of total elements */
extern int   NTotLines;       /* number of total lines */

extern Hashtable_type HNN;
extern Hashtable_type HEL;
extern Hashtable_type HLI;
extern QueueTable      CM;              /* list of comments */


#endif

