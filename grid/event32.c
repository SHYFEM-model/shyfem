
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

/**
 ** EVENT32.C
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

#ifndef __GNUC__
#error  This file is only for the DJGPP 32-bit version!!!!
#endif

#include <stdlib.h>
#include "eventque.h"

asm(								       "\n\
	.data								\n\
	.align  4							\n\
_mousedraw_painter:		/* pointer to the function */		\n\
	.long	0		/* to draw mouse */			\n\
_mousedraw_active_p:		/* pointer to flag to zero */		\n\
	.long	0		/* when drawing is done */		\n\
_mousedraw_contaddr_p:		/* pointer to dword containing */	\n\
	.long	0		/* return address from mouse draw */	\n\
	.text								\n\
	.align  2,144							\n\
_mousedraw_func:							\n\
	cli								\n\
	pushl	%eax		/* place for return address */		\n\
	pushf								\n\
	pushl	%eax		/* save EAX */				\n\
	movl	_mousedraw_contaddr_p,%eax				\n\
	movl	(%eax),%eax	/* fix up return address */		\n\
	movl	%eax,8(%esp)						\n\
	pushl	%ebx							\n\
	pushl	%ecx							\n\
	pushl	%edx							\n\
	pushl	%esi							\n\
	pushl	%edi							\n\
	movl	_mousedraw_painter,%eax					\n\
	sti								\n\
	call	*%eax							\n\
	cli								\n\
	popl	%edi							\n\
	popl	%esi							\n\
	popl	%edx							\n\
	popl	%ecx							\n\
	popl	%ebx							\n\
	movl	_mousedraw_active_p,%eax				\n\
	movb	$0,(%eax)	/* clear active flag */			\n\
	popl	%eax							\n\
	popf								\n\
	sti								\n\
	ret			/* back to program */			  "
);

static int  have_queue = 0;

/*
 * These are actually local symbols at the link level, we just have to
 * trick the C compiler
 */
extern void mousedraw_func(void);
extern void (*mousedraw_painter)(void);
extern char *mousedraw_active_p;
extern long *mousedraw_contaddr_p;


void EventQueueDeInit(void)
{
	if(have_queue) {
	    asm volatile(					       "\n\
		movl   $0x00ff,%%eax					\n\
		xorl   %%ebx,%%ebx					\n\
		int    $0x33						  "
		: /* nothing */
		: /* nothing */
		: "ax", "bx", "cx", "dx"
	    );
	    have_queue = 0;
	}
}

EventQueue *EventQueueInit(int qsize,int ms_stksize,void (*msdraw)(void))
{
	EventQueue *queue;
	int ack;

	if(qsize < 20) qsize = 20;
	if(msdraw != NULL) {
	    mousedraw_painter = msdraw;
	    msdraw = mousedraw_func;
	}
	asm volatile(						       "\n\
	    movl   $0x00ff,%%eax					\n\
	    movl   %2,%%ebx						\n\
	    movl   %3,%%ecx						\n\
	    int	   $0x33						\n\
	    movl   %%eax,%0						\n\
	    movl   %%ebx,%1						\n\
	    movl   %%ecx,_mousedraw_contaddr_p				\n\
	    movl   %%edx,_mousedraw_active_p				  "
	    : "=g" (ack),   "=g" (queue)
	    : "g"  (qsize), "g"  (msdraw)
	    : "ax", "bx", "cx", "dx"
	);
	if(ack != 0x0ff0) queue = NULL;
	have_queue = (queue != NULL) ? 1 : 0;
	return(queue);
}

