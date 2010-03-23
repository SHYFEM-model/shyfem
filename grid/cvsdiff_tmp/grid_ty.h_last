
/************************************************************************\
 *									*
 * grid_ty.h - type definitions for grid                                *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see grid.c for copying information					*
 *									*
 * Revision History:							*
 * 07-May-1998: type is now integer                                     *
 * 02-Apr-1998: no Button_type, Function_type                           *
 * 13-Oct-97: New Button_type REMOVE_ELEMENT, REMOVE_LINE               *
 * 10-Oct-97: New Button_type SAVE, UNIFY_NODE                          *
 * 06-Dec-95: ScaleUp/DownVect                                          *
 * 06-Dec-95: Wind_type removed, level from Node_type removed           *
 * 04-Dec-95: Button_type changed, Vect_type added                      *
 * 11-Mar-95: Conn_type is now local to gridhs.c                        *
 * 06-Oct-94: ColTab_type defined locally in gridop                     *
 * 04-May-94: new LINE modes and Line_type                              *
 * 13-Apr-94: splitted up in df, ty, fp, ex                             *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
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
