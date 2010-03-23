/**
 ** EVENTQUE.H
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

#ifndef _EVENTQUE_H_
#define _EVENTQUE_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * structures:
 *  BE CAREFUL when hacking!!! -- 16 and 32 bit compilers have to generate
 *  the same alignments
 */
typedef struct {
    unsigned char   evt_type;	    /* event type: 0: keyboard, 1: mouse */
    unsigned char   evt_kbstat;	    /* keyboard status (ALT, SHIFT, etc..) */
    unsigned char   evt_mask;	    /* mouse event mask */
    unsigned char   evt_button;	    /* button status */
    unsigned short  evt_xpos;	    /* X coord (or keycode if keybd event) */
    unsigned short  evt_ypos;	    /* Y coord */
    unsigned long   evt_time;	    /* time stamp of event */
#define evt_keycode   evt_xpos	    /* reuse this slot for keybd events !! */
#define evt_scancode  evt_ypos	    /* store here the BIOS scan code */
} EventRecord;

typedef struct {
    unsigned short  evq_maxsize;    /* max size of event queue */
    unsigned short  evq_cursize;    /* number of events in the queue */
    unsigned short  evq_rdptr;	    /* next event to read */
    unsigned short  evq_wrptr;	    /* next event to be written */
    short	    evq_xpos;	    /* current X coordinate of mouse */
    short	    evq_ypos;	    /* current Y coordinate of mouse */
    short	    evq_xmin;	    /* minimal mouse X coordinate */
    short	    evq_ymin;	    /* minimal mouse Y coordinate */
    short	    evq_xmax;	    /* maximal mouse X coordinate */
    short	    evq_ymax;	    /* maximal mouse Y coordinate */
    short	    evq_xspeed;	    /* horizontal speed (mickey/coord) */
    short	    evq_yspeed;	    /* vertical speed (mickey/coord) */
    unsigned short  evq_thresh;	    /* fast movement threshold */
    unsigned short  evq_accel;	    /* multiplier for fast move */
    unsigned char   evq_drawmouse;  /* interrupt handler has to draw mouse */
    unsigned char   evq_moved;	    /* set if mouse moved */
    unsigned char   evq_delchar;    /* character removed from BIOS buffer */
    unsigned char   evq_enable;	    /* event generation control flag */
    EventRecord	    evq_events[1];  /* event buffer space */
} EventQueue;

/*
 * event types
 */
#define EVENT_KEYBD	0
#define EVENT_MOUSE	1

/*
 * MOUSE event flag bits
 * (also defined in "mousex.h" of the graphics library)
 */
#ifndef M_MOTION

#define M_MOTION	0x001
#define M_LEFT_DOWN	0x002
#define M_LEFT_UP	0x004
#define M_RIGHT_DOWN	0x008
#define M_RIGHT_UP	0x010
#define M_MIDDLE_DOWN	0x020
#define M_MIDDLE_UP	0x040
#define M_BUTTON_DOWN	(M_LEFT_DOWN | M_MIDDLE_DOWN | M_RIGHT_DOWN)
#define M_BUTTON_UP	(M_LEFT_UP   | M_MIDDLE_UP   | M_RIGHT_UP)
#define M_BUTTON_CHANGE (M_BUTTON_UP | M_BUTTON_DOWN )

/*
 * MOUSE button status bits
 */
#define M_LEFT		1
#define M_RIGHT		2
#define M_MIDDLE	4

#endif  /* M_MOTION */

/*
 * KEYBOARD status word bits
 * (also defined in "mousex.h" of the graphics library)
 */
#ifndef KB_SHIFT

#define KB_RIGHTSHIFT	0x01		/* right shift key depressed */
#define KB_LEFTSHIFT	0x02		/* left shift key depressed */
#define KB_CTRL		0x04		/* CTRL depressed */
#define KB_ALT		0x08		/* ALT depressed */
#define KB_SCROLLOCK	0x10		/* SCROLL LOCK active */
#define KB_NUMLOCK	0x20		/* NUM LOCK active */
#define KB_CAPSLOCK	0x40		/* CAPS LOCK active */
#define KB_INSERT	0x80		/* INSERT state active */

#define KB_SHIFT	(KB_LEFTSHIFT | KB_RIGHTSHIFT)

#endif  /* KB_SHIFT */

/*
 * set this bit in 'evq_enable' to generate the corresponding event
 */
#define EVENT_ENABLE(type)	(1 << (type))

/*
 * prototypes
 */
#if defined(__TURBOC__) && defined(FOR_GO32)
EventQueue *EventQueueInit(int qsize,int ms_stksize,void (*msdraw)(void),int,int);
#else
EventQueue *EventQueueInit(int qsize,int ms_stksize,void (*msdraw)(void));
#endif

void   EventQueueDeInit(void);
int    EventQueueNextEvent(EventQueue *q,EventRecord *e);

#ifdef __cplusplus
}
#endif

#endif /* whole file */


