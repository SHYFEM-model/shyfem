
/************************************************************************\ 
 *									*
 * debug.h - debugging routines						*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see keybd.c for copying information					*
 *									*
 * Revision History:							*
 * 31-Dec-94: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_DEBUG_
#define __GUH_DEBUG_


/**************************************************************************/

void GDPrint( void );
void GDInit( int size );
char *GDAlloc( int size );
void GDNL( void );
void GDS( char *s );
void GDI( int i );
void GDL( long l );
void GDF( float f );

/**************************************************************************/

#endif
