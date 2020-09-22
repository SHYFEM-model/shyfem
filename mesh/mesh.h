
/************************************************************************\
 *
 *    Copyright (C) 1995,1997  Georg Umgiesser
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
 * mesh.h - automatic meshing routines
 *
 * revision log :
 *
 * 25.07.1995	ggu	routines written from scratch
 * 08.10.1997	ggu	standard structures uniformed
 * ...		ggu	new extra structures Extra_N/E/L_type
 * ...		ggu	new node type N_GIVEN
 *
\************************************************************************/



#ifndef __GUH_MESH_
#define __GUH_MESH_


#include "fund.h"
#include "list.h"
#include "queue.h"
#include "hash.h"

/********************* special structures *****************/

typedef struct {
	float *x;
	float *y;
        int type;
        int ndim;
        int number;
} Bound_type;

typedef struct {
	float a[3];
	float b[3];
	float c[3];
	float v[3];
	float fr;
} Backg_type;

/********************* extra structures *****************/

typedef struct {
	int type;
} Extra_N_type;

typedef struct {
	int type;
	Backg_type *backg;
} Extra_E_type;

typedef struct {
	int type;
} Extra_L_type;

/********************* standard structures *****************/

typedef struct node_tag {
        int number;
        int type;
        Point coord;
	float depth;
	void *extra;
} Node_type;

typedef struct elem_tag {
        int number;
        int type;
        int vertex;
        int *index;
	float depth;
	void *extra;
			/* next tags should go into extra structure */
	int *neibor;
	int flag;
	float rho;
	Point rc;
} Elem_type;

typedef struct line_tag {
        int number;
        int type;
        int vertex;
        int *index;
	void *extra;
} Line_type;

/********************* other *********************/

/* mark null depth */

#define NULLDEPTH -999.0

/* right or left node */

#define	RIGHT		0
#define	LEFT		1

/* flags for elements */

#define	FL_DEFAULT	0
#define	FL_VISITED	1
#define	FL_DELETED	2

/* type of neighbourhood */

#define	NB_NONE		 0
#define	NB_UNKNOWN	 1	/* Element Number 1 ist not used */
#define	NB_DELETED	-2	/* negative for deleted */

#define NO_NEIBOR(n)		( (n) == NB_NONE )
#define UNKNOWN_NEIBOR(n)	( (n) == NB_UNKNOWN )
#define DELETED_NEIBOR(n)	( (n) < 0 )
#define ACTIVE_NEIBOR(n)	( (n) > 0 )

/* type of elements used */

#define E_NONE		0
#define E_EXTERNAL	1
#define E_INTERNAL	2

/* type of nodes used */

#define N_NONE		0
#define N_EXTERNAL	1
#define N_BOUNDARY	2
#define N_ADDBOUND	4
#define N_INTERNAL	8
#define N_GIVEN		16	/* node read from file but not on boundary */

/* type of lines used */

#define L_NONE			0
#define L_EXTERNAL		1
#define L_INTERNAL		2
#define L_FAULT			4	/* internal open boundary line */
#define L_EXTERNAL_REF		8
#define L_INTERNAL_REF		16
#define L_FAULT_REF		32

/*********************************/

extern int   NTotNodes;       /* number of total nodes */
extern int   NTotElems;       /* number of total elements */
extern int   NTotLines;       /* number of total lines */

extern Hashtable_type HNN;
extern Hashtable_type HEL;
extern Hashtable_type HLI;
extern ListTable      NEL;              /* list of new elements */
extern ListTable      BCK;              /* list of background grid elements */
extern QueueTable      CM;              /* list of comments */


float Resol( float x , float y );


#endif

