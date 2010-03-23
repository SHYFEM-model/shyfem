
/************************************************************************\ 
 *									*
 * exgrd.h - extracts items from grd file				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see exgrd.c for copying information					*
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

