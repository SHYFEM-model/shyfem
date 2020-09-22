
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * grdut.h - grd utility routines
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshut
 *
\************************************************************************/


#ifndef __GUH_GRDUT_
#define __GUH_GRDUT_


#include "fund.h"
#include "hash.h"

#include "grd.h"

/*-------------------------------------------------------------------*/

int *MakeIndex( int vertex );
int *CopyIndex( int vertex , int *index );
void InvertIndex( int *index , int nvert );

/*-------------------------------------------------------------------*/

void PrintNodes( Hashtable_type H );
void PrintElems( Hashtable_type H );
void PrintLines( Hashtable_type H );
void PrintLineS( Hashtable_type HL, Hashtable_type HN, int line, int invert );

/*-------------------------------------------------------------------*/

Node_type *MakeNode( int n , int ntype , Point *c );
Elem_type *MakeElem( int n , int type , int vertex , int *index );
Line_type *MakeLine( int n , int type , int vertex , int *index );

void DeleteNode( Node_type *p );
void DeleteElem( Elem_type *p );
void DeleteLine( Line_type *p );

int InsertNewNode( int type, float x, float y );
int InsertNewElem( int type, int vertex, int *index );
int InsertNewLine( int type, int vertex, int *index );

/*-------------------------------------------------------------------*/

int GetTotNodes( void );
int GetTotElems( void );
int GetTotLines( void );

int IncTotNodes( void );
int IncTotElems( void );
int IncTotLines( void );

void SetTotNodes( int n );
void SetTotElems( int n );
void SetTotLines( int n );

/*-------------------------------------------------------------------*/

Grid_type *MakeGrid( void );
void SetGrid( Grid_type *G );
Grid_type *GetGrid( void );

/*-------------------------------------------------------------------*/


#endif


