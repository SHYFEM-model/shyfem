
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
 *                                                                      *
 * meshty.h - routines handling new mesh type                           *
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

