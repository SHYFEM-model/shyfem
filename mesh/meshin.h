
/************************************************************************\
 *
 *    Copyright (C) 1995,1997  Georg Umgiesser
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
 *
 * meshin.h - node insert routines for mesh
 *
 * revision log :
 *
 * 01.08.1995	ggu	routines written from scratch
 * 08.10.1997	ggu	new routine prototypes InsertNodes(), GivenNodes()
 *
\************************************************************************/


#ifndef __GUH_MESHIN_
#define __GUH_MESHIN_


#include "mesh.h"
#include "stack.h"
#include "nlist.h"

/*
#include "hash.h"
#include "nlist.h"
#include "heap.h"
*/


int VisitElems( Elem_type *from , Elem_type *into , StackTable deleted
                        , StackTable visited
                        , float x , float y );
void RefineInternalNodes( void );
void InsertInternalNodes( void );
void InsertBoundaryNodes( NodeList nlist );
void InsertNodes( NodeList nlist );
NodeList GivenNodes( void );
void RecoverElements( StackTable deleted , StackTable visited );
int InsertNode( int node, Elem_type *pe, float x, float y
                        , StackTable deleted
                        , StackTable visited
                        , StackTable created);
void UpdateNeighbors( Elem_type *pold, int node, int dir, StackTable L );
void SmoothInternalNodes( void );


#endif


