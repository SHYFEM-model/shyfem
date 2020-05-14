
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1997,2004,2011  Georg Umgiesser
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
 * gridop.c - set up options from command line and colortable
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 11.05.1994	ggu	OpCheckUse and OpElim2Lines added, OpType = 0
 * 10.02.1995	ggu	color routines copied to seperate file
 * 16.02.1995	ggu	introduced OpMaxColDepth and OpColTabSize
 * 09.03.1995	ggu	Logos in gridvs with version description
 * 06.12.1995	ggu	OpWind, OpIntern, options 0,1,2... eliminated
 * 17.12.1997	ggu	OpColor - color nodes and lines (...was compress...)
 * 01.10.2004	ggu	set default for OpMaxColDepth to -1 (not given)
 * 16.02.2011	ggu	new options OpOutFile and OpItemType implemented
 *
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "grid.h"
#include "graph.h"
#include "gustd.h"
#include "options.h"

#if __GUG_UNIX_
#include "xgraph.h"
#endif

static void Usage(void);
void Logos(void);	/* in file gridvs.c */
static void Help(void);
static void Assistance(void);

void SetOptions(int argc, char *argv[])

{
	int c;
	char *options = "CTM:S:kfoahw:c:g:d:uDV:N:O:t:";

	      OpCheck = 0;    /* do more checks */
	       OpFill = 0;    /* fill elements */
	      OpColor = 0;    /* color nodes and lines */
	     OpBorder = 1;    /* outline elements */
	   OpShowType = 0;    /* show type -- 0 is depth, else type */
	       OpType = 0;    /* file type */
	        OpCol = 0;    /* number of color table */
	   OpCompress = 0;    /* compress node and element numbers */
           OpCheckUse = 0;    /* checks for used nodes */
	 OpElim2Lines = 0;    /* eliminates double end points */
	OpMaxColDepth = -1.;  /* maximum depth to scale color */
	 OpColTabSize = 20;   /* size of color table */
	   OpVectFact = 1.;   /* factor for vector */
	   OpNodeFact = 1.;   /* factor for node not in use */
	   OpVectScal = 1.;   /* scaling for vector */
	   OpNodeScal = 1.;   /* scaling for node not in use */
	   OpOutFile  = NULL; /* name of output file (ask if not given) */
	   OpItemType = 0.;   /* default type for created items */

	while( (c=GetOpt(argc,argv,options)) != EOF ) {
		switch(c) {
        case 'h' :              /* h - help screen */
			Logos();
			Help();
			Assistance();
			exit(0);
			break;
        case 'C' :              /* C - color nodes and lines */
			OpColor = 1;
			break;
        case 'T' :              /* T - chose color by type, not by depth */
			OpShowType = 1;
			break;
        case 'D' :              /* D - eliminates double end points */
			OpElim2Lines = 1;
			break;
        case 'M' :              /* M - maximum depth to scale color */
			OpMaxColDepth = atof(GetOptArg());
			break;
        case 'S' :              /* S - size of color table */
			OpColTabSize = atoi(GetOptArg());
			break;
        case 'N' :              /* N - scaling for nodes */
			OpNodeScal = atof(GetOptArg());
			break;
        case 'V' :              /* V - scaling for vectors */
			OpVectScal = atof(GetOptArg());
			break;
        case 'u' :              /* u - check use of nodes */
			OpCheckUse = 1;
			OpCheck = 1;
			break;
        case 'k' :              /* k - do checks on data integrity */
			OpCheck = 1;
			break;
        case 'f' :              /* f - fill triangles */
			OpFill = 1;
			break;
        case 'o' :              /* o - do not outline elements */
			OpBorder = 0;
			break;
        case 'a' :              /* a - ask for file names */
			OpType = -1;
			break;
        case 'c' :              /* c - number of color table */
			OpCol = *GetOptArg() - '0';
			break;
        case 'O' :              /* O - name of output file */
			OpOutFile = savestr(GetOptArg());
			break;
        case 't' :              /* t - default type of items */
			OpItemType = atoi(GetOptArg());
			break;
#if __GUG_UNIX_
        case 'd' :              /* d - display name */
			QSetDisplay(GetOptArg());
			break;
        case 'g' :              /* g - geometry */
			QSetGeometry(GetOptArg());
			break;
#else
		case 'd' :
		case 'g' :
			break;
#endif
		default :
			exit(1);
			break;
		}
	}

	OpArgc = GetOptInd();

	if( argc <= OpArgc && OpType != -1 ) {
		Logos();
		Usage();
		exit(0);
	} else {
		Logos();
	}
}

void Usage( void )

{
        printf("Usage : grid [-h] [-options] grd-file(s)\n");
}

void Assistance( void )

{
	printf("For assistance please contact :\n");
	printf("  Georg Umgiesser, ISMAR-CNR, georg.umgiesser@ismar.cnr.it\n");
	printf("\n");
}

void Help( void )

{
	printf("Options :\n");
	printf("  -o   do not outline elements        ");
	printf("  -f   fill elements with color     \n");
	printf("  -k   do extra checking              ");
	printf("  -u   check if nodes are used      \n");
	printf("  -T   show type instead of depth     ");
	printf("  -c#  use color table #            \n");
	printf("  -h   print this help screen         ");
	printf("  -a   ask for file names           \n");
	printf("  -d   display  (only X11)            ");
	printf("  -g   geometry (only X11)          \n");
	printf("  -M#  scale color to depth #         ");
	printf("  -S#  size of color table is #     \n");
	printf("  -N#  scale factor for nodes is #    ");
	printf("  -V#  scale factor for vectors is #\n");
	printf("  -C   color nodes and lines          ");
	printf("  -On  use n as output file name    \n");
	printf("  -t#  use type # for new items     \n");
	printf("\n");
}

