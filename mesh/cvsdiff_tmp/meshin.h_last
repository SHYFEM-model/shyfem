
/************************************************************************\
 *									*
 * meshin.h - node insert routines for mesh				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshin.c for copying information					*
 *									*
 * Revision History:							*
 * 08-Oct-97: new routine prototypes InsertNodes(), GivenNodes()        *
 * 01-Aug-95: routines written from scratch				*
 *									*
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


