
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


#include "general.h"
#include "grd.h"

/****************************************************************/

void MakeNodelList( void );
void PrintNodel( int first );
int GetBoxOfNodel( int first , int second );
void MakeLinelList( void );
int CheckNeibBoxes(int number,int *boxes,int *boxesc,int nm);
void WriteNeibBoxes(int box, int *boxes,int *boxesc,int nm, int nbox);
void PrintLine( Line_type *pl );

int *GetNextBox( int nm , int *boxes , int * ind , int *n , int *nbox );
void CheckLinel( void );
void MakeLines( void );
void MakeNodes( void );
