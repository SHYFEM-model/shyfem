
/************************************************************************\ 
 *									*
 * grd.h - standard header for grd files				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see grdio.c for copying information					*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines written from scratch				*
 *									*
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

