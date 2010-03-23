
/* $Id: menu.h,v 1.2 2009-01-14 17:16:37 georg Exp $ */

/************************************************************************\
 *									*
 * menu.h - menu routines                                               *
 *									*
 * Copyright (c) 1998 by Georg Umgiesser				*
 *									*
 * see menu.c for copying information                                   *
 *                                                                      *
 * Revision History:                                                    *
 * 02-Apr-1998: version not yet finished but working                    *
 * 28-Mar-1998: routines adopted from gridmu.c				*
 *                                                                      *
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
