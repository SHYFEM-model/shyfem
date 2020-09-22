
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * grdhs.c - hash table utilities for grd files
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshhs
 *
\************************************************************************/


#include "hash.h"
#include "grd.h"


/*******************************************************************\
 *************** General Hash Table Utility Routines ***************
\*******************************************************************/


Node_type *VisitHashTableN( Hashtable_type H )

{
	return (Node_type *) VisitHashTable( H );
}

Elem_type *VisitHashTableE( Hashtable_type H )

{
        return (Elem_type *) VisitHashTable( H );
}

Line_type *VisitHashTableL( Hashtable_type H )

{
        return (Line_type *) VisitHashTable( H );
}


void InsertByNodeNumber( Hashtable_type H , Node_type *nodep )

{
	InsertHashByNumber( H , (void *) nodep , nodep->number );
}

void InsertByElemNumber( Hashtable_type H , Elem_type *elemp )

{
	InsertHashByNumber( H , (void *) elemp , elemp->number );
}

void InsertByLineNumber( Hashtable_type H , Line_type *linep )

{
	InsertHashByNumber( H , (void *) linep , linep->number );
}


Node_type *RetrieveByNodeNumber( Hashtable_type H , int node )

{
	return (Node_type *) RetrieveHashByNumber( H , node );
}

Elem_type *RetrieveByElemNumber( Hashtable_type H , int elem )

{
	return (Elem_type *) RetrieveHashByNumber( H , elem );
}

Line_type *RetrieveByLineNumber( Hashtable_type H , int line )

{
	return (Line_type *) RetrieveHashByNumber( H , line );
}


