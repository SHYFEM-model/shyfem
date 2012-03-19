
/************************************************************************\
 *									*
 * gridop.c - set up options from command line and colortable		*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
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
 * 16-Feb-2011: new options OpOutFile and OpItemType implemented	*
 * 01-Oct-2004: set default for OpMaxColDepth to -1 (not given)         *
 * 17-Dec-97: OpColor - color nodes and lines (...was compress...)      *
 * 06-Dec-95: OpWind, OpIntern, options 0,1,2... eliminated             *
 * 09-Mar-95: Logos in gridvs with version description                  *
 * 16-Feb-95: introduced OpMaxColDepth and OpColTabSize                 *
 * 10-Feb-95: color routines copied to seperate file                    *
 * 11-May-94: OpCheckUse and OpElim2Lines added, OpType = 0             *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "grid.h"
#include "graph.h"
#include "gustd.h"

#if __GUG_UNIX_
#include "xgraph.h"
#endif

int  GetOpt(int argc, char *argv[], char *optionS);
void Usage(void);
void Logos(void);	/* in file gridvs.c */
void Help(void);
void Assistance(void);

int     optind  = 1;    /* index of which argument is next      */
char   *optarg;         /* pointer to argument of current option */
int     opterr  = 1;    /* allow error message  */

static  char   *letP = NULL;    /* remember next option char's location */
static  char    SW = 0;         /* DOS switch character, either '-' or '/' */

/*
  Parse the command line options, System V style.

  Standard option syntax is:

    option ::= SW [optLetter]* [argLetter space* argument]

  where
    - SW is either '/' or '-', according to the current setting
      of the MSDOS switchar (int 21h function 37h).
    - there is no space before any optLetter or argLetter.
    - opt/arg letters are alphabetic, not punctuation characters.
    - optLetters, if present, must be matched in optionS.
    - argLetters, if present, are found in optionS followed by ':'.
    - argument is any white-space delimited string.  Note that it
      can include the SW character.
    - upper and lower case letters are distinct.

  There may be multiple option clusters on a command line, each
  beginning with a SW, but all must appear before any non-option
  arguments (arguments not introduced by SW).  Opt/arg letters may
  be repeated: it is up to the caller to decide if that is an error.

  The character SW appearing alone as the last argument is an error.
  The lead-in sequence SWSW ("--" or "//") causes itself and all the
  rest of the line to be ignored (allowing non-options which begin
  with the switch char).

  The string *optionS allows valid opt/arg letters to be recognized.
  argLetters are followed with ':'.  Getopt () returns the value of
  the option character found, or EOF if no more options are in the
  command line.  If option is an argLetter then the global optarg is
  set to point to the argument string (having skipped any white-space).

  The global optind is initially 1 and is always left as the index
  of the next argument of argv[] which getopt has not taken.  Note
  that if "--" or "//" are used then optind is stepped to the next
  argument before getopt() returns EOF.

  If an error occurs, that is an SW char precedes an unknown letter,
  then getopt() will return a '?' character and normally prints an
  error message via perror().  If the global variable opterr is set
  to false (zero) before calling getopt() then the error message is
  not printed.

  For example, if the MSDOS switch char is '/' (the MSDOS norm) and

    *optionS == "A:F:PuU:wXZ:"

  then 'P', 'u', 'w', and 'X' are option letters and 'F', 'U', 'Z'
  are followed by arguments.  A valid command line may be:

    aCommand  /uPFPi /X /A L someFile

  where:
    - 'u' and 'P' will be returned as isolated option letters.
    - 'F' will return with "Pi" as its argument string.
    - 'X' is an isolated option.
    - 'A' will return with "L" as its argument.
    - "someFile" is not an option, and terminates getOpt.  The
      caller may collect remaining arguments using argv pointers.
*/

int     GetOpt(int argc, char *argv[], char *optionS)
{
	unsigned char ch;
	char *optP;

	if (SW == 0) {
		/* get SW using dos call 0x37 */
/*                _AX = 0x3700;
		geninterrupt(0x21);
  */              SW = '-';
	}

	if (argc > optind) {
		if (letP == NULL) {
			if ((letP = argv[optind]) == NULL ||
				*(letP++) != SW)  goto gopEOF;
			if (*letP == SW) {
				optind++;  goto gopEOF;
			}
		}
		if (0 == (ch = *(letP++))) {
			optind++;  goto gopEOF;
		}
		if (':' == ch  ||  (optP = strchr(optionS, ch)) == NULL)
			goto gopError;
		if (':' == *(++optP)) {
			optind++;
			if (0 == *letP) {
				if (argc <= optind)  goto  gopError;
				letP = argv[optind++];
			}
			optarg = letP;
			letP = NULL;
		} else {
			if (0 == *letP) {
				optind++;
				letP = NULL;
			}
			optarg = NULL;
		}
		return ch;
	}
gopEOF:
	optarg = letP = NULL;
	return EOF;

gopError:
	optarg = NULL;
	if (opterr)
		printf("Unknown command line option : %c\n",ch);
	return ('?');
}

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
			OpMaxColDepth = atof(optarg);
			break;
        case 'S' :              /* S - size of color table */
			OpColTabSize = atoi(optarg);
			break;
        case 'N' :              /* N - scaling for nodes */
			OpNodeScal = atof(optarg);
			break;
        case 'V' :              /* V - scaling for vectors */
			OpVectScal = atof(optarg);
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
			OpCol = *optarg - '0';
			break;
        case 'O' :              /* O - name of output file */
			OpOutFile = savestr(optarg);
			break;
        case 't' :              /* t - default type of items */
			OpItemType = atoi(optarg);
			break;
#if __GUG_UNIX_
        case 'd' :              /* d - display name */
			QSetDisplay(optarg);
			break;
        case 'g' :              /* g - geometry */
			QSetGeometry(optarg);
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

	OpArgc = optind;

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
        printf("Usage : grid [-options] [files]");
        printf("   ( grid -h  for help)\n");
        printf("\n");
}

void Assistance( void )

{
	printf("For assistance please contact :\n");
	printf("  Georg Umgiesser, ISDGM/CNR\n");
	printf("  S.Polo 1364, 30125 Venezia, Italy\n");
	printf("  Tel.: ++39-41-5216875  Fax: ++39-41-2602340\n");
	printf("  E-Mail : georg@lagoon.isdgm.ve.cnr.it\n");
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

