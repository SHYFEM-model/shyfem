
/************************************************************************\
 *
 *    Copyright (C) 2010,2018  Georg Umgiesser
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

/** 
 ** NEXTEVNT.C 
 **
 **  Copyright (C) 1992, Csaba Biegl
 **    820 Stirrup Dr, Nashville, TN, 37221
 **    csaba@vuse.vanderbilt.edu
 **
 **  This file is distributed under the terms listed in the document
 **  "copying.cb", available from the author at the address above.
 **  A copy of "copying.cb" should accompany this file; if not, a copy
 **  should be available from where this file was obtained.  This file
 **  may not be distributed without a verbatim copy of "copying.cb".
 **  You should also have received a copy of the GNU General Public
 **  License along with this program (it is in the file "copying");
 **  if not, write to the Free Software Foundation, Inc., 675 Mass Ave,
 **  Cambridge, MA 02139, USA.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **/

/************************************************************************\
 *
 * revision log :
 *
 * 23.03.2010	ggu	changed v6.1.1
 * 18.12.2018	ggu	changed VERS_7_5_52
 *
\************************************************************************/

#include "eventque.h"

#ifdef __GNUC__
# define disable()  asm volatile("cli");
# define enable()   asm volatile("sti");
#endif

#ifdef __TURBOC__
# include <dos.h>
#endif

int EventQueueNextEvent(EventQueue *q,EventRecord *e)
{
	if(q->evq_cursize > 0) {
	    disable();
	    *e = q->evq_events[q->evq_rdptr];
	    if(++q->evq_rdptr == q->evq_maxsize) q->evq_rdptr = 0;
	    q->evq_cursize--;
	    enable();
	    return(1);
	}
	return(0);
}

