
/************************************************************************\
 *									*
 * meshge.h - geometric routines for mesh				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshge.c for copying information					*
 *									*
 * Revision History:							*
 * 05-Dec-2011: new general routines					*
 * 17-Oct-97: new routine TurnClosedLine()                              *
 * 16-Oct-97: new routine FindElement()                                 *
 * 15-Oct-97: new routine InClosedLine()                                *
 * 01-Aug-95: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_MESHGE_
#define __GUH_MESHGE_


#include "mesh.h"
#include "hash.h"
#include "list.h"
#include "nlist.h"

/*
#include "nlist.h"
#include "heap.h"
*/


void MakeCircumCircle( Elem_type *pe );
void ControlCircumCircle( Elem_type *pe );
float MakeInCircleRadius( Elem_type *pe );
void MakeCM( Elem_type *pe, float *xm, float *ym );

int InCircumCircle( Elem_type *pe , float x , float y );
int InElemCircumCircle( Elem_type *circle , Elem_type *pe , float fact );
int InConvex( int n , float *xe , float *ye , float x , float y );
int InElement( Hashtable_type H , Elem_type *pe , float x , float y );
void PolyMinMax( Hashtable_type H , Elem_type *pe , Rect *r );
void PolyMinMaxIndex( Hashtable_type H , int nvert , int *index , Rect *r );

Elem_type *FindElement( Hashtable_type H, float x, float y );
Elem_type *FindXYElement( ListTable L, float x, float y );

float AreaElement( Hashtable_type H , Elem_type *pe );
float AreaLine( Hashtable_type H , Line_type *pl );
float AreaConvex( NodeList list );
float LengthConvex( NodeList list );

double angle( float x1, float y1, float x2, float y2, float x3, float y3 );

int InClosedLine( int ndim , float *xl , float *yl , float x , float y );
int TurnClosedLine( int ndim , float *xl , float *yl );

int IsPointInLine( Line_type *pl , float x , float y );
int IsLineInLine( Line_type *plext , Line_type *plint );
int IsLineClosed( Line_type *pl );


#endif


