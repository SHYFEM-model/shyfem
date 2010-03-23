
/************************************************************************\ 
 *									*
 * grdut.h - grd utility routines					*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see grdut.c for copying information					*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines copied from meshut                               *
 *									*
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


