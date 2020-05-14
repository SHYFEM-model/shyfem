
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1997-1998  Georg Umgiesser
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
 * grid_ty.h - type definitions for grid
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	splitted up in df, ty, fp, ex
 * 04.05.1994	ggu	new LINE modes and Line_type
 * 06.10.1994	ggu	ColTab_type defined locally in gridop
 * 11.03.1995	ggu	Conn_type is now local to gridhs.c
 * 04.12.1995	ggu	Button_type changed, Vect_type added
 * 06.12.1995	ggu	Wind_type removed, level from Node_type removed
 * 06.12.1995	ggu	ScaleUp/DownVect
 * 10.10.1997	ggu	New Button_type SAVE, UNIFY_NODE
 * 13.10.1997	ggu	New Button_type REMOVE_ELEMENT, REMOVE_LINE
 * 02.04.1998	ggu	no Button_type, Function_type
 * 07.05.1998	ggu	type is now integer
 *
\************************************************************************/


#ifndef __GUH_GRID_TY_
#define __GUH_GRID_TY_


#include "fund.h"


typedef enum {
NO_INPUT , MENU_FIELD_INPUT , PLOT_FIELD_INPUT , KEYBOARD_INPUT
} Mode_type;


typedef struct {
	Point coord;
	float depth;
	void *extra;
	int type;
	int use;
	int number;
} Node_type;

typedef struct {
	int vertex;
	int *index;
	float depth;
	int type;
	int number;
} Elem_type;

typedef struct {
	int vertex;
	int *index;
	float depth;
	int type;
	int number;
} Line_type;

typedef struct {
	int total;
	int actual;
	float *speed;
	float *dir;
} Vect_type;


#endif
