
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


/************************************************************************\
 *									*
 * meshbd.h - boundary routines for mesh				*
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


