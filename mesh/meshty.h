
/************************************************************************\
 *                                                                      *
 * meshty.h - routines handling new mesh type                           *
 *                                                                      *
 * Copyright (c) 1997 by Georg Umgiesser                                *
 *                                                                      *
 * see meshty.c for copying information                                 *
 *                                                                      *
 * Revision History:                                                    *
 * 08-Oct-97: routine written from scratch                              *
 *                                                                      *
\************************************************************************/


#ifndef __GUH_MESHTY_
#define __GUH_MESHTY_


#include "mesh.h"

void SetNtype( Node_type *pn , int Type );
void AddNtype( Node_type *pn , int Type );
int GetNtype( Node_type *pn );
int IsNtype( Node_type *pn , int Type );
int EqualsNtype( Node_type *pn , int Type );

void SetEtype( Elem_type *pe , int Type );
void AddEtype( Elem_type *pe , int Type );
int GetEtype( Elem_type *pe );
int IsEtype( Elem_type *pe , int Type );
int EqualsEtype( Elem_type *pe , int Type );

void SetLtype( Line_type *pl , int Type );
void AddLtype( Line_type *pl , int Type );
int GetLtype( Line_type *pl );
int IsLtype( Line_type *pl , int Type );
int EqualsLtype( Line_type *pl , int Type );


#endif

