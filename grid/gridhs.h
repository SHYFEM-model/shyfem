
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
 *									*
 * gridhs.h - hash aministration routines and check of data strucures   *
 *									*
 * Revision History:							*
 * 14-Oct-97: new routine CheckUseConsistency()                         *
 * 19-Sep-97: Call to CheckNodes changed                                *
 * 11-Mar-95: static routines not listed anymore (local to gridhs.c)    *
 * 13-Apr-94: completely restructured -> uses hash.c to do work         *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRIDHS_
#define __GUH_GRIDHS_


#include "hash.h"
#include "grid_ty.h"



Hashtable_type GetHashTableN( void );
Hashtable_type GetHashTableE( void );
Hashtable_type GetHashTableL( void );

Node_type *VisitHashTableN( Hashtable_type H );
Elem_type *VisitHashTableE( Hashtable_type H );
Line_type *VisitHashTableL( Hashtable_type H );

void InsertByNodeNumber( Hashtable_type H , Node_type *nodep );
void InsertByElemNumber( Hashtable_type H , Elem_type *elemp );
void InsertByLineNumber( Hashtable_type H , Line_type *linep );

Node_type *RetrieveByNodeNumber( Hashtable_type H , int node );
Elem_type *RetrieveByElemNumber( Hashtable_type H , int elem );
Line_type *RetrieveByLineNumber( Hashtable_type H , int line );

void CheckMore( Hashtable_type HN , Hashtable_type HE , Hashtable_type HL );
void CheckNodes( Hashtable_type HN , Hashtable_type HE , Hashtable_type HL 
			, Hashtable_type HV );

void CheckUseConsistency( void );


#endif

