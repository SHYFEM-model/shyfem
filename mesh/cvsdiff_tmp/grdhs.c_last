
/************************************************************************\ 
 *									*
 * grdhs.c - hash table utilities for grd files				*
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines copied from meshhs				*
 *									*
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


