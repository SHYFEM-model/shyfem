
/* $Id: assert.h,v 1.2 2009-01-14 17:16:37 georg Exp $ */

/************************************************************************\ 
 *									*
 * assert.h - utility routines for assertion handling                   *
 *									*
 * Copyright (c) 1996 by Georg Umgiesser				*
 *									*
 * see word.c for copying information					*
 *									*
 * Revision History:							*
 * ..-...-96: routines copied from "Writing Solid Code"                 *
 *									*
\************************************************************************/


#ifndef __GUH_ASSERT_
#define __GUH_ASSERT_

#ifdef ASSERT_DEBUG

void _Assert(char *filename, unsigned line);	/* prototype */

#define ASSERT(f)		\
	if(f)			\
		;		\
	else			\
		_Assert(__FILE__,__LINE__)

#else

#define ASSERT(f)	;

#endif /* ASSERT_DEBUG */

#endif /* __GUH_ASSERT_ */

