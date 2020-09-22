
/************************************************************************\
 *
 *    Copyright (C) 1996  Georg Umgiesser
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
 * assert.c - utility routines for assertion handling
 *
 * revision log :
 *
 * 01.01.1996	ggu	routines copied from "Writing Solid Code"
 *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"

void _Assert(char *filename, unsigned line)

{
	fprintf(stderr,"\nAssertion failed: %s, line %u\n",
			filename,line);
	exit(99);
}

