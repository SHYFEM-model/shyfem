
/************************************************************************\
 *
 *    Copyright (C) 1995,1997-1998,2008-2009,2018  Georg Umgiesser
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
 * exgrdvs.c - version description of exgrd
 *
 * revision log :
 *
 * 17.08.1995	ggu	version 1.00
 * 18.08.1995	ggu	version 1.10 (manip./unify nodes, command line options)
 * 05.05.1997	ggu	version 1.20 (CompressNumbers)
 * 17.01.1998	ggu	version 1.21 (UnifyNodes)
 * 15.07.1998	ggu	version 1.22 (Number Ranges)
 * 01.11.2008	ggu	version 1.23 (handle extension .grd)
 * 05.11.2008	ggu	version 1.24 (handle clockwise elements)
 * 14.01.2009	ggu	version 1.25 (SHYFEM tree)
 * 02.11.2018	ggu	version 1.30 (copyright)
 *
\************************************************************************/

#include <stdio.h>

char* SCopy  = "Copyright (C) 1995-2018  Georg Umgiesser                ";
char* SExgrd = "EXGRD - Extract Items from GRD Files - Version 1.30     ";

void Logos( void )

{
        fprintf(stderr,"\n");
	fprintf(stderr,"%s\n",SExgrd);
	fprintf(stderr,"%s\n",SCopy);
	fprintf(stderr,"\n");
}

/**********************************************************************\


===========================================================================
========================== version log for mesh ===========================
===========================================================================

todo

-f	write to stdout
-i	info
-u	better algorithm

===========================================================================

version 1.30							02 Nov 2018

copyright updated
use getopt library from libc

===========================================================================

version 1.24							05 Nov 2008

handle clockwise elements (option -a)

===========================================================================

version 1.23							01 Nov 2008

handle extension .grd

===========================================================================

version 1.22							15 Jul 1998

select by number ranges added

===========================================================================

version 1.21							17 Jan 98

algorithm changed for UnifyNodes()

=========================================================================

version 1.20							05 May 97

CompressNumbers introduced

=========================================================================

version 1.10							18 Aug 95	

finished manipulation of nodes
written code for unify nodes
command line options implemented

=========================================================================

version 1.00							17 Aug 95	

started project -> some files reused from mesh
	new grd.h, grdio, grdhs, grdut for general read
	new options.c to read command line option

=========================================================================

\**********************************************************************/
