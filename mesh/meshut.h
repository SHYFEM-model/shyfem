
/************************************************************************\ 
 *									*
 * meshut.h - mesh utility routines					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshut.c for copying information					*
 *									*
 * Revision History:							*
 * 16-Oct-97: new routines PrintLineList(), MakeNewLine()               *
 *              InsertInIndex(), InsertNodeInLine()                     *
 * 15-Oct-97: new routine MakeCoordsFromLine()                          *
 * 01-Aug-95: new routines inserted and hash routines transfered        *
 * 25-Jul-95: routines written from scratch				*
 *									*
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


