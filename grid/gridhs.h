
/************************************************************************\ 
 *									*
 * gridhs.h - hash aministration routines and check of data strucures   *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see gridhs.c for copying information					*
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

