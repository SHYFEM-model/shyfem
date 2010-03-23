
/************************************************************************\
 *									*
 * exgrdvs.c - version description of exgrd                             *
 *									*
 * Copyright (c) 1995-2009 by Georg Umgiesser				*
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
 * 14-Jan-2009: version 1.25 (SHYFEM tree)				*
 * 05-Nov-2008: version 1.24 (handle clockwise elements)		*
 * 01-Nov-2008: version 1.23 (handle extension .grd)                    *
 * 15-Jul-1998: version 1.22 (Number Ranges)                            *
 * 17-Jan-98: version 1.21 (UnifyNodes)                                 *
 * 05-May-97: version 1.20 (CompressNumbers)                            *
 * 18-Aug-95: version 1.10 (manip./unify nodes, command line options)   *
 * 17-Aug-95: version 1.00                                              *
 *									*
\************************************************************************/

#include <stdio.h>

char* SCopy  = "Copyright (c) Georg Umgiesser 1995 - 2009               ";
char* SExgrd = "EXGRD - Extract Items from GRD Files - Version 1.25     ";
char* SGeorg = "        1995-2009 (c) Georg Umgiesser - ISDGM/CNR       ";

void Logos( void )

{
        fprintf(stderr,"\n");
	fprintf(stderr,"%s\n",SExgrd);
	fprintf(stderr,"%s\n",SGeorg);
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
