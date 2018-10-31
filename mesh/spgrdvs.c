
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
 * spgrdvs.c - version description of spgrd                             *
 *									*
 * Revision History:							*
 * 12-Feb-1998: version 1.01                                            *
 * ..-Feb-1998: version 1.00                                            *
 *									*
\************************************************************************/

#include <stdio.h>

char* SCopy  = "Copyright (c) Georg Umgiesser 1995 - 1998               ";
char* SExgrd = "SPGRD - Extract Lines from GRD Files - Version 1.01     ";
char* SGeorg = "        1995-1998 (c) Georg Umgiesser - ISDGM/CNR       ";

void Logos( void )

{
        fprintf(stderr,"\n");
	fprintf(stderr,"%s\n",SExgrd);
	fprintf(stderr,"%s\n",SGeorg);
	fprintf(stderr,"\n");
}

/**********************************************************************\


=========================================================================
========================== version log for mesh =========================
=========================================================================

version 1.00 & 1.01					      12 Feb 1998	

started project -> some files reused from mesh

=========================================================================

\**********************************************************************/
