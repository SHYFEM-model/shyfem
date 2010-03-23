
/************************************************************************\
 *                                                                      *
 * meshty.c - routines handling new mesh type                           *
 *                                                                      *
 * Copyright (c) 1997 by Georg Umgiesser                                *
 *                                                                      *
 * Permission to use, copy, modify, and distribute this software        *
 * and its documentation for any purpose and without fee is hereby      *
 * granted, provided that the above copyright notice appear in all      *
 * copies and that both that copyright notice and this permission       *
 * notice appear in supporting documentation.                           *
 *                                                                      *
 * This file is provided AS IS with no warranties of any kind.          *
 * The author shall have no liability with respect to the               *
 * infringement of copyrights, trade secrets or any patents by          *
 * this file or any part thereof.  In no event will the author          *
 * be liable for any lost revenue or profits or other special,          *
 * indirect and consequential damages.                                  *
 *                                                                      *
 * Comments and additions should be sent to the author:                 *
 *                                                                      *
 *                      Georg Umgiesser                                 *
 *                      ISDGM/CNR                                       *
 *                      S. Polo 1364                                    *
 *                      30125 Venezia                                   *
 *                      Italy                                           *
 *                                                                      *
 *                      Tel.   : ++39-41-5216875                        *
 *                      Fax    : ++39-41-2602340                        *
 *                      E-Mail : georg@lagoon.isdgm.ve.cnr.it           *
 *                                                                      *
 * Revision History:                                                    *
 * 08-Oct-97: routine written from scratch                              *
 *                                                                      *
\************************************************************************/


#include "general.h"
#include "assert.h"
#include "mesh.h"
#include "meshty.h"

/****************************************************************************/

void SetNtype( Node_type *pn , int Type )

{
	ASSERT( pn );
	((Extra_N_type *)pn->extra)->type = Type;
}

void AddNtype( Node_type *pn , int Type )

{
	ASSERT( pn );
	((Extra_N_type *)pn->extra)->type |= Type;
}

int GetNtype( Node_type *pn )

{
	ASSERT(pn);
        return ( ((Extra_N_type *)pn->extra)->type );
}

int IsNtype( Node_type *pn , int Type )

{
	ASSERT(pn);
        return ( ((Extra_N_type *)pn->extra)->type & Type );
}

int EqualsNtype( Node_type *pn , int Type )

{
	ASSERT(pn);
        return ( ((Extra_N_type *)pn->extra)->type == Type );
}

/****************************************************************************/

/****************************************************************************/

void SetEtype( Elem_type *pe , int Type )

{
	ASSERT( pe );
	((Extra_E_type *)pe->extra)->type = Type;
}

void AddEtype( Elem_type *pe , int Type )

{
	ASSERT( pe );
	((Extra_E_type *)pe->extra)->type |= Type;
}

int GetEtype( Elem_type *pe )

{
	ASSERT(pe);
        return ( ((Extra_E_type *)pe->extra)->type );
}

int IsEtype( Elem_type *pe , int Type )

{
	ASSERT(pe);
        return ( ((Extra_E_type *)pe->extra)->type & Type );
}

int EqualsEtype( Elem_type *pe , int Type )

{
	ASSERT(pe);
        return ( ((Extra_E_type *)pe->extra)->type == Type );
}

/****************************************************************************/

/****************************************************************************/

void SetLtype( Line_type *pl , int Type )

{
	ASSERT( pl );
	((Extra_L_type *)pl->extra)->type = Type;
}

void AddLtype( Line_type *pl , int Type )

{
	ASSERT( pl );
	((Extra_L_type *)pl->extra)->type |= Type;
}

int GetLtype( Line_type *pl )

{
	ASSERT(pl);
        return ( ((Extra_L_type *)pl->extra)->type );
}

int IsLtype( Line_type *pl , int Type )

{
	ASSERT(pl);
        return ( ((Extra_L_type *)pl->extra)->type & Type );
}

int EqualsLtype( Line_type *pl , int Type )

{
	ASSERT(pl);
        return ( ((Extra_L_type *)pl->extra)->type == Type );
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

