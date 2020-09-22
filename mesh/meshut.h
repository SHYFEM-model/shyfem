
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
 * meshut.h - mesh utility routines
 *
 * revision log :
 *
 * 25.07.1995	ggu	routines written from scratch
 * 01.08.1995	ggu	new routines inserted and hash routines transfered
 * 15.10.1997	ggu	new routine MakeCoordsFromLine()
 * 16.10.1997	ggu	new routines PrintLineList(), MakeNewLine()
 * ...		ggu	InsertInIndex(), InsertNodeInLine()
 *
\************************************************************************/


#ifndef __GUH_MESHUT_
#define __GUH_MESHUT_


#include "fund.h"
#include "hash.h"
#include "nlist.h"
#include "mesh.h"


void PrintNodes( Hashtable_type H );
void PrintLineS( Hashtable_type HL, Hashtable_type HN, int line, int invert );
void PrintLineList( Hashtable_type H );
void PrintLines( Hashtable_type H );
void PrintLine( Line_type *p );
void PrintElems( Hashtable_type H );
void PrintNodeList( char *string , NodeList list );

int *MakeIndex( int vertex );
int *CopyIndex( int vertex , int *index );
void InvertIndex( int *index , int nvert );
int *InsertInIndex( int vertex , int *index , int after , int new );

Node_type *MakeNode( int n , int ntype , Point *c );
void DeleteNode( Node_type *p );

Elem_type *MakeElemWithIndex( int n , int ntype , int vertex , int *index );
Elem_type *MakeElem( int n , int *c , int vertex );
void DeleteElem( Elem_type *p );

Line_type *MakeLineWithIndex( int n , int ntype , int vertex , int *index );
Line_type *MakeLine( int n , int *c , int vertex );
void DeleteLine( Line_type *p );
void InsertNodeInLine( Line_type *pl , int after , int new );

int InsertNewNode( int type, float x, float y );
int InsertNewElem( int type, int n1, int n2, int n3
                        , int nb1, int nb2, int nb3 );
int InsertNewLine( int type, int vertex, int *index );
Line_type *MakeNewLine( int type, int vertex, int *index );

int UpdateOldElemByNumber( int elem, int type, int n1, int n2, int n3
                        , int nb1, int nb2, int nb3 );
int UpdateOldElem( Elem_type *pe, int type, int n1, int n2, int n3
                        , int nb1, int nb2, int nb3 );

int UpdateNeiborByNumber( int elem, int nb1, int nb2, int nb3 );
int UpdateNeibor( Elem_type *pe, int nb1, int nb2, int nb3 );

void InsertCircumCircle( Elem_type *pe );

void MakeCoordsFromLine( Line_type *pl , float **x , float **y );


#endif


