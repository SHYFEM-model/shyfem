
/************************************************************************\
 *									*
 * meshop.c - command line options for mesh				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
 *									*
 * Revision History:							*
 * 02-Aug-95: routines copied from gridop.c				*
 *									*
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "gustd.h"
#include "options.h"

#include "mesh.h"
#include "meshop.h"


static void Usage(void);
void Logos(void);		/* in file meshvs.c */
static void Help(void);
static void Assistance(void);
static void Test(void);


int	       OpArgc = 0;    /* number of items after options */
int     OpTestVersion = 0;    /* test version with limited nodes */

int	 OpBackGround = -1;   /* element type for background grid */
int	     OpIntern = 0;    /* total number of intern nodes */
int	   OpBoundary = 0;    /* total number of boundary nodes */
int  OpRefineBoundary = -1;   /* refine boundary with resolution */
int    OpInsertIntern = -1;   /* insert internal points or not */
int      OpSmoothPass = 0;    /* number of passes for smoothing (s < 100) */
float   OpSmoothOmega = 0.1;  /* omega for smoothing (0.05 < o < 0.25) */
float	     OpAspect = 0.;   /* aspect ratio for internal node refinement */
float	 OpResolution = 0.;   /* overall resolution */

void SetOptions(int argc, char *argv[])

{
	int c;
	char *options = "hbnI:B:R:a:s:o:g:";

	while( (c=GetOpt(argc,argv,options)) != EOF ) {
		switch(c) {
		case 'h' :              /* h - help screen */
			Logos();
			Test();
			Help();
			Assistance();
			exit(0);
			break;
		case 'a' :		/* must be > 2, probably > 3 */
			OpAspect = atof(GetOptArg());
			break;
		case 'b' :
			OpRefineBoundary = 0;
			break;
		case 'g' :
			OpBackGround = atoi(GetOptArg());
			break;
		case 'n' :
			OpInsertIntern = 0;
			break;
		case 's' :		/* number of passes for smoothing */
			OpSmoothPass = atoi(GetOptArg());
			break;
		case 'o' :		/* omega for smoothing */
			OpSmoothOmega = atof(GetOptArg());
			break;
		case 'I' :              /* total number of intern nodes */
			OpIntern = atoi(GetOptArg());
			break;
		case 'B' :              /* total number of boundary nodes */
			OpBoundary = atoi(GetOptArg());
			break;
		case 'R' :              /* overall resolution */
			OpResolution = atof(GetOptArg());
			break;
		default :
			exit(1);
			break;
		}
	}

	OpArgc = GetOptInd();

	if( argc <= OpArgc ) {
		Logos();
		Test();
		Usage();
		exit(0);
	} else {
		Logos();
		Test();
	}

	if( OpRefineBoundary == -1 ) {
		if( OpBoundary || OpIntern || OpResolution > 0. ) {
			OpRefineBoundary = 1;
		} else {
			OpRefineBoundary = 0;
		}
	}
	if( OpInsertIntern == -1 ) {
		if( OpIntern || OpResolution > 0. ) {
			OpInsertIntern = 1;
		} else {
			OpInsertIntern = 0;
		}
	}
}

static void Usage( void )

{
        printf("Usage : mesh [-options] file [files]");
        printf("   ( mesh -h  for help)\n");
        printf("\n");
}

static void Assistance( void )

{
	printf("For assistance please contact :\n");
	printf("  Georg Umgiesser, ISDGM/CNR\n");
	printf("  S.Polo 1364, 30125 Venezia, Italy\n");
	printf("  Tel.: ++39-41-5216875  Fax: ++39-41-2602340\n");
	printf("  E-Mail : georg@lagoon.isdgm.ve.cnr.it\n");
	printf("\n");
}

static void Help( void )

{
	printf("Options :\n");
	printf("  -b   do not refine boundaries      ");
	printf("  -n   do not insert internal nodes\n");
	printf("  -s#  passes for smoothing          ");
	printf("  -o#  relax. par. for smoothing   \n");
	printf("  -I#  number of internal nodes      ");
	printf("  -B#  number of boundary nodes    \n");
	printf("  -R#  overall resolution            ");
	printf("  -a#  obtain this aspect ratio    \n");
	printf("  -g#  element type for background   ");
	printf("  -h   show this help screen       \n");
	printf("\n");
}

static void Test( void )

{
	if( OpTestVersion ) {
		printf("Test version - for evaluation only\n");
		printf("\n");
	}
}

void TestVersion( void )

{
	if( OpTestVersion && OpTestVersion < NTotNodes ) {
		Error("Only test version: too much nodes");
	}
}
