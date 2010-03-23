
/************************************************************************\
 *									*
 * meshbd.h - boundary routines for mesh				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshbd.c for copying information					*
 *									*
 * Revision History:							*
 * 16-Oct-97: new routines RecoverBoundaryNodes(), FindElemToSide()     *
 * 01-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_MESHBD_
#define __GUH_MESHBD_


#include "queue.h"
#include "nlist.h"


int TotBoundNodes( void );
NodeList RefineBoundary( void );
void RefineLine( QueueTable line, float x1, float y1
                        , float x2, float y2 , float dist );
void SortLine( QueueTable list );
void CopyBoundaryLine( void );

void RecoverBoundary( void );
void RecoverInternalFault( void );
void RecoverBoundaryNodes( int linetype );

Elem_type *FindElemToSide( int node1 , int node2 );


#endif


