
/************************************************************************\
 *
 *    Copyright (C) 1998  Georg Umgiesser
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
 * menu.h - menu routines
 *
 * revision log :
 *
 * 28.03.1998	ggu	routines adopted from gridmu.c
 * 02.04.1998	ggu	version not yet finished but working
 *
\************************************************************************/


#ifndef __GUH_MENU_
#define __GUH_MENU_


#include <stdio.h>

#include "general.h"
#include "screen.h"



#define MENU_NONE		0
#define MENU_BAR		1
#define MENU_AREA		2
#define MENU_PULLDOWN		3
#define MENU_EXEC		4
#define MENU_RADIO		5
#define MENU_CHECK		6

#define INPUT_NONE		0
#define INPUT_MENU		1
#define INPUT_KEYBOARD		2
#define INPUT_PLOT		3


typedef void (*FP)( void );


typedef struct menu_entry_tag {
	char		*name;			/* name of menu item */
	int		 type;			/* type of menu item */
	FP		 function;		/* fuction called on click */
	int		 nsubs;			/* number of submenu entries */
	IRect		 position;		/* position of menu on screen */
	IRect		 subpos;		/* position of submenu */
	void		*pixmap;		/* to save pixmap of menu */
	struct menu_entry_tag	*parent;	/* parent of this menu */
	struct menu_entry_tag	*submenu;	/* submenus if pulldown */
	struct menu_entry_tag	*next;		/* next sibling */
} Menu_entry;


/***********************************************************************/


void RegisterMainMenu( Menu_entry *main );
void DeleteMainMenu( void );
Menu_entry *MakeMenuBar( void );
Menu_entry *MakeMenuArea( void );
Menu_entry *MakePulldownMenu( char *name );
Menu_entry *MakeExecMenu( char *name , FP function );
void AddMenuItem( Menu_entry *parent , Menu_entry *submenu );

void SetWindowDimension( int xmin , int ymin , int xmax , int ymax );
void SetMenuDimension( int xmin , int ymin , int xmax , int ymax );

void DisplayMainMenu( void );

void ProcessMenuInput( int x , int y );

void AdjustPulldownMenu( int x , int y );
void CleanupPulldownMenu( void );
void ReleasePulldownMenu( Menu_entry *menu );
Menu_entry *FindPulldownMenu( int x , int y );
void ExecutePulldownMenu( int x , int y );
void SchedulePulldownMenu( int x , int y );
void DisplayPulldownMenu( Menu_entry *menu , int x , int y );
void PlotPulldownMenu( Menu_entry *menu );
void ExecuteCommand( Menu_entry *menu );

void PrintMenuInfo( Menu_entry *menu );
void PrintSubMenuInfo( Menu_entry *menu );



#endif 
