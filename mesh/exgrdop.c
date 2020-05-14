
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
 * exgrdop.c - command line options for exgrd
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshop.c
 *
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "gustd.h"
#include "options.h"

#include "grd.h"
#include "grdut.h"
#include "exgrdop.h"


static void Usage(void);
void Logos(void);		/* in file exgrdvs.c */
static void Help(void);
static void Assistance(void);
static void Test(void);


int	       OpArgc = 0;    /* number of items after options */
int     OpTestVersion = 0;    /* test version with limited nodes */

int OpInfo = 0;
int OpExtractNodes = 0;
int OpDeleteNodes = 0;
int OpExtractElems = 0;
int OpDeleteElems = 0;
int OpExtractLines = 0;
int OpDeleteLines = 0;
int OpExtractUnusedNodes = 0;
int OpDeleteUnusedNodes = 0;
int OpMinType = -1;
int OpMaxType = -1;
int OpMinVertex = -1;
int OpMaxVertex = -1;
int OpMinRange = -1;
int OpMaxRange = -1;
int OpUnifyNodes = 0;
int OpCompressNumbers = 0;
int OpMakeAntiClockwise = 0;
int OpDeleteStrangeElements = 0;
float OpMinDepth = NULLDEPTH;
float OpMaxDepth = NULLDEPTH;
float OpTollerance = 0.;


void SetOptions(int argc, char *argv[])

{
	int c;
	char *options = "hinNeElLsSt:T:d:D:v:V:uo:Caxr:R:";

	while( (c=GetOpt(argc,argv,options)) != EOF ) {
		switch(c) {
		case 'h' :              /* h - help screen */
			Logos();
			Test();
			Usage();
			Help();
			Assistance();
			exit(0);
			break;
		case 'i' :
			OpInfo = 1;
			break;
		case 'n' :
			OpExtractNodes = 1;
			break;
		case 'N' :
			OpDeleteNodes = 1;
			break;
		case 'e' :
			OpExtractElems = 1;
			break;
		case 'E' :
			OpDeleteElems = 1;
			break;
		case 'l' :
			OpExtractLines = 1;
			break;
		case 'L' :
			OpDeleteLines = 1;
			break;
		case 's' :
			OpExtractUnusedNodes = 1;
			break;
		case 'S' :
			OpDeleteUnusedNodes = 1;
			break;
		case 't' :
			OpMinType = atoi(GetOptArg());
			break;
		case 'T' :
			OpMaxType = atoi(GetOptArg());
			break;
		case 'd' :
			OpMinDepth = atof(GetOptArg());
			break;
		case 'D' :
			OpMaxDepth = atof(GetOptArg());
			break;
		case 'v' :
			OpMinVertex = atoi(GetOptArg());
			break;
		case 'V' :
			OpMaxVertex = atoi(GetOptArg());
			break;
		case 'r' :
			OpMinRange = atoi(GetOptArg());
			break;
		case 'R' :
			OpMaxRange = atoi(GetOptArg());
			break;
		case 'u' :
			OpUnifyNodes = 1;
			break;
		case 'o' :
			OpTollerance = atof(GetOptArg());
			break;
		case 'C' :
			OpCompressNumbers = 1;
			break;
		case 'a' :
			OpMakeAntiClockwise = 1;
			break;
		case 'x' :
			OpDeleteStrangeElements = 1;
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

	if( OpExtractNodes && OpDeleteNodes )
		Error("Options 'n' and 'N' are mutually exclusive.");
	if( OpExtractUnusedNodes && OpDeleteUnusedNodes )
		Error("Options 's' and 'S' are mutually exclusive.");
	if( ( OpExtractUnusedNodes || OpDeleteUnusedNodes ) &&
	   		 ( OpExtractNodes || OpDeleteNodes ) )
		Error("Options 'nN' and 'sS' are mutually exclusive.");
	if( OpExtractElems && OpDeleteElems )
		Error("Options 'e' and 'E' are mutually exclusive.");
	if( OpExtractLines && OpDeleteLines )
		Error("Options 'l' and 'L' are mutually exclusive.");

}

static void Usage( void )

{
        fprintf(stderr,"Usage : exgrd [-h] [-options] grd-file(s)\n");
}

static void Assistance( void )

{
	fprintf(stderr,"For assistance please contact :\n");
        printf("  Georg Umgiesser, ISMAR-CNR, georg.umgiesser@ismar.cnr.it\n");
	fprintf(stderr,"\n");
}

static void Help( void )

{
	fprintf(stderr,"Options :\n");
	fprintf(stderr,"  -i   show only info on file        ");
	fprintf(stderr,"  -h   show this help screen       \n");
	fprintf(stderr,"  -n   extract nodes                 ");
	fprintf(stderr,"  -N   delete nodes                \n");
	fprintf(stderr,"  -e   extract elements              ");
	fprintf(stderr,"  -E   delete elements             \n");
	fprintf(stderr,"  -l   extract lines                 ");
	fprintf(stderr,"  -L   delete lines                \n");
	fprintf(stderr,"  -s   extract unused nodes          ");
	fprintf(stderr,"  -S   delete unused nodes         \n");
	fprintf(stderr,"  -t#  select from type              ");
	fprintf(stderr,"  -T#  select up to type           \n");
	fprintf(stderr,"  -d#  select from depth             ");
	fprintf(stderr,"  -D#  select up to depth          \n");
	fprintf(stderr,"  -v#  select from vertices          ");
	fprintf(stderr,"  -V#  select up to vertices       \n");
	fprintf(stderr,"  -r#  select from number range      ");
	fprintf(stderr,"  -R#  select up to number range   \n");
	fprintf(stderr,"  -u   unify nodes                   ");
	fprintf(stderr,"  -o#  tollerance for unify (def=0)\n");
	fprintf(stderr,"  -C   compress numbers              ");
	fprintf(stderr,"  -a   make items anti-clockwise   \n");
	fprintf(stderr,"  -x   delete degenerate elements  \n");
	fprintf(stderr,"\n");
}

static void Test( void )

{
	if( OpTestVersion ) {
		fprintf(stderr,"Test version - for evaluation only\n");
		fprintf(stderr,"\n");
	}
}

void TestVersion( void )

{
	if( OpTestVersion && OpTestVersion < GetTotNodes() ) {
		Error("Only test version: too much nodes");
	}
}
