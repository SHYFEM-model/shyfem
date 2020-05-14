
/************************************************************************\
 *
 *    Copyright (C) 1997  Georg Umgiesser
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
 * meshty.c - routines handling new mesh type
 *
 * revision log :
 *
 * 08.10.1997	ggu	routine written from scratch
 *
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

