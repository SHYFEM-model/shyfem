
/************************************************************************\ 
 *									*
 * meshhs.h - mesh hash utility routines				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see meshhs.c for copying information					*
 *									*
 * Revision History:							*
 * 01-Aug-95: routines copied from meshut                               *
 *									*
\************************************************************************/


#ifndef __GUH_MESHHS_
#define __GUH_MESHHS_


#include "mesh.h"


Node_type *VisitHashTableN( Hashtable_type H );
Elem_type *VisitHashTableE( Hashtable_type H );
Line_type *VisitHashTableL( Hashtable_type H );
void InsertByNodeNumber( Hashtable_type H , Node_type *nodep );
void InsertByElemNumber( Hashtable_type H , Elem_type *elemp );
void InsertByLineNumber( Hashtable_type H , Line_type *linep );
Node_type *RetrieveByNodeNumber( Hashtable_type H , int node );
Elem_type *RetrieveByElemNumber( Hashtable_type H , int elem );
Line_type *RetrieveByLineNumber( Hashtable_type H , int line );


#endif


