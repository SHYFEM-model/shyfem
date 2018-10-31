
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
 * spgrdop.c - command line options for spgrd				*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines copied from meshop.c				*
 *									*
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


static void Usage( int help );
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
int OpUnifyNodes = 0;
int OpCompressNumbers = 0;
float OpMinDepth = NULLDEPTH;
float OpMaxDepth = NULLDEPTH;
float OpTollerance = 0.;


void SetOptions(int argc, char *argv[])

{
	int c;
	char *options = "hinNeElLsSt:T:d:D:v:V:uo:C";

	while( (c=GetOpt(argc,argv,options)) != EOF ) {
		switch(c) {
		case 'h' :              /* h - help screen */
			Logos();
			Test();
			Usage(FALSE);
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
		case 'u' :
			OpUnifyNodes = 1;
			break;
		case 'o' :
			OpTollerance = atof(GetOptArg());
			break;
		case 'C' :
			OpCompressNumbers = 1;
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
		Usage(TRUE);
		exit(1);
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

static void Usage( int help )

{
        fprintf(stderr,"Usage : spgrd [-options] file [files]");
	if( help )
            fprintf(stderr,"   (spgrd -h  for help)\n");
	else
            fprintf(stderr,"\n");
        fprintf(stderr,"\n");
}

static void Assistance( void )

{
	fprintf(stderr,"For assistance please contact :\n");
	fprintf(stderr,"  Georg Umgiesser, ISDGM/CNR\n");
	fprintf(stderr,"  S.Polo 1364, 30125 Venezia, Italy\n");
	fprintf(stderr,"  Tel.: ++39-41-5216875  Fax: ++39-41-2602340\n");
	fprintf(stderr,"  E-Mail : georg@lagoon.isdgm.ve.cnr.it\n");
	fprintf(stderr,"\n");
}

static void Help( void )

{
	fprintf(stderr,"Options :\n");
	fprintf(stderr,"  -i   show only info on file        ");
	fprintf(stderr,"  -h   show this help screen       \n");
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
