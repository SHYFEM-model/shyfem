
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
 * options.c - command line options
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshop.c
 *
\************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>


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
  command line.  If option is an argLetter then the global Optarg is
  set to point to the argument string (having skipped any white-space).

  The global Optind is initially 1 and is always left as the index
  of the next argument of argv[] which getopt has not taken.  Note
  that if "--" or "//" are used then Optind is stepped to the next
  argument before getopt() returns EOF.

  If an error occurs, that is an SW char precedes an unknown letter,
  then getopt() will return a '?' character and normally prints an
  error message via perror().  If the global variable Opterr is set
  to false (zero) before calling getopt() then the error message is
  not printed.

  For example,

    *optionS == "A:F:PuU:wXZ:"

  then 'P', 'u', 'w', and 'X' are option letters and 'F', 'U', 'Z'
  are followed by arguments.  A valid command line may be:

    aCommand  -uPFPi -X -A L someFile

  where:
    - 'u' and 'P' will be returned as isolated option letters.
    - 'F' will return with "Pi" as its argument string.
    - 'X' is an isolated option.
    - 'A' will return with "L" as its argument.
    - "someFile" is not an option, and terminates getOpt.  The
      caller may collect remaining arguments using argv pointers.
*/



static  int     Optind  = 1;    /* index of which argument is next       */
static  char   *Optarg;         /* pointer to argument of current option */
static  int     Opterr  = 1;    /* allow error message                   */

static  char   *letP = NULL;    /* remember next option char's location  */
static  char    SW = '-';       /* option character                      */

int GetOptInd( void ) { return Optind; }
char *GetOptArg( void ) { return Optarg; }

int GetOpt(int argc, char *argv[], char *optionS)

{
	unsigned char ch;
	char *optP;

	if (argc > Optind) {
		if (letP == NULL) {
			if ((letP = argv[Optind]) == NULL ||
				*(letP++) != SW)  goto gopEOF;
			if (*letP == SW) {
				Optind++;  goto gopEOF;
			}
		}
		if (0 == (ch = *(letP++))) {
			Optind++;  goto gopEOF;
		}
		if (':' == ch  ||  (optP = strchr(optionS, ch)) == NULL)
			goto gopError;
		if (':' == *(++optP)) {
			Optind++;
			if (0 == *letP) {
				if (argc <= Optind)  goto  gopError;
				letP = argv[Optind++];
			}
			Optarg = letP;
			letP = NULL;
		} else {
			if (0 == *letP) {
				Optind++;
				letP = NULL;
			}
			Optarg = NULL;
		}
		return ch;
	}
gopEOF:
	Optarg = letP = NULL;
	return EOF;

gopError:
	Optarg = NULL;
	if (Opterr)
		printf("Unknown command line option : %c\n",ch);
	return ('?');
}
