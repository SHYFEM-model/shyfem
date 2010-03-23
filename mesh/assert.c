
/* $Id: assert.c,v 1.3 2009-01-14 14:49:13 georg Exp $ */

/************************************************************************\ 
 *									*
 * assert.c - utility routines for assertion handling                   *
 *									*
 * Copyright (c) 1996 by Georg Umgiesser				*
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
 * ..-...-96: routines copied from "Writing Solid Code"                 *
 *									*
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

