
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
 ** EVTEST.C
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <dos.h>

EventQueue *q;

#define ATTR 0x70;

void printloc(void)
{
#ifdef __TURBOC__
	char far *p = MK_FP(0xb800,0);
#endif
#ifdef __GNUC__
	char *p = (char *)(0xe00b8000L);
#endif
	static char hexdigit[] = { "0123456789abcdef" };
	int x = q->evq_xpos;
	int y = q->evq_ypos;
	int i,j;

	for(i = 0; i < 2; i++) {
	    *p++ = i ? 'Y' : 'X';
	    *p++ = ATTR;
	    *p++ = '=';
	    *p++ = ATTR;
	    *p++ = '0';
	    *p++ = ATTR;
	    *p++ = 'x';
	    *p++ = ATTR;
	    for(j = 0; j < 4; j++) {
		*p++ = hexdigit[(x >> 12) & 0x0f];
		*p++ = ATTR;
		x <<= 4;
	    }
	    for(j = 0; j < 4; j++) {
		*p++ = ' ';
		*p++ = ATTR;
	    }
	    x = y;
	}
}

void main(void)
{
	EventRecord e;
	long ii;

	q = EventQueueInit(100,320,printloc);
	if(q == NULL) {
	    printf("could not init queues\n");
	    exit(1);
	}
	q->evq_drawmouse = 1;
	for( ; ; ) {
	    if(!EventQueueNextEvent(q,&e)) {
		for(ii = 0L; ii < 1000000L; ii++);
		continue;
	    }
	    if(e.evt_type == EVENT_MOUSE) {
		printf("MOUSE event: mask 0x%02x, btn 0x%02x, "
		       "kbstat 0x%02x xpos %3d, ypos %3d, time %5ld\n",
		    e.evt_mask,
		    e.evt_button,
		    e.evt_kbstat,
		    e.evt_xpos,
		    e.evt_ypos,
		    e.evt_time
		);
		continue;
	    }
	    if(e.evt_keycode == 0x1b)
		break;
	    if((e.evt_keycode <= 0xff) && isprint(e.evt_keycode)) {
		printf("KEYBD event: key '%c', kbstat 0x%02x, time %5ld\n",
		    e.evt_keycode,
		    e.evt_kbstat,
		    e.evt_time
		);
		continue;
	    }
	    printf("KEYBD event: key 0x%03x, kbstat 0x%02x, time %5ld\n",
		e.evt_keycode,
		e.evt_kbstat,
		e.evt_time
	    );
	}
	EventQueueDeInit();
	exit(0);
}

