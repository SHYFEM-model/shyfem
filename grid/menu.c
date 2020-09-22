
/************************************************************************\
 *
 *    Copyright (C) 1998,2003  Georg Umgiesser
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
 * menu.c - menu routines
 *
 * revision log :
 *
 * 28.03.1998	ggu	routines adopted from gridmu.c
 * 02.04.1998	ggu	version not yet finished but working
 * 19.11.2003	ggu	debug messages only if edebug is true
 *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"

#include "general.h"
#include "screen.h"
#include "graph.h"

#include "mouse.h"
#include "events.h"

#include "menu.h"
#include "stack.h"

#include "grid_fp.h"	/* HACK */

static int edebug = FALSE;         /* set to TRUE for debugging */


/***********************************************************************/


static Menu_entry *MainMenu = 0;	/* main menu (menu bar/area) */

static Menu_entry *TopMenu = 0;		/* top displayed menu */

static IRect WinDim = {0,0,0,0};	/* actual dimension of total window */
static IRect MenDim = {0,0,0,0};	/* actual dimension of menu bar/area */

/* static int FgMainMenuCol; */	/* colors for main menu */
static int BgMainMenuCol;
static int BlMainMenuCol;
static int BdMainMenuCol;

static int FgButtonCol;		/* colors for buttons in main menu */
static int BgButtonCol;
static int BlButtonCol;
static int BdButtonCol;

static int FgPulldownCol;	/* colors for pulldown menu */
static int BgPulldownCol;
static int BlPulldownCol;
static int BdPulldownCol;


/********************** prototypes *************************/


static void DisplayMenuBar( Menu_entry *main );
static void DisplayMenuArea( Menu_entry *main );

static void ProcessEvent( QEvent *event , int *button , int *x , int *y );

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

static void RegisterPulldownDimensions( Menu_entry *menu
                 , int *width , int *height , int *x0 , int *y0 );

static void FindMaxTextDimension( Menu_entry *menu 
			, int *width , int *height , int *totwidth );

static void HighlightMenu( Menu_entry *menu );
static void UnHighlightMenu( Menu_entry *menu );

/********************** static routines *************************/


static Menu_entry *MakeMenuItem( char *name , int type )

{
	Menu_entry *new;

	new = (Menu_entry *) malloc( sizeof(Menu_entry) );
	if( !new ) Error("Cannot make menu entry");

	new->name = (char *) malloc( strlen(name) + 1 );
	if( !new->name ) Error("Cannot make name of menu entry");
	strcpy(new->name,name);

	new->type = type;
	new->function = NULL;
	new->nsubs = 0;
	(void) ToIRect(&(new->position),0,0,0,0);
	(void) ToIRect(&(new->subpos),0,0,0,0);
	new->pixmap = NULL;
	new->parent = NULL;
	new->submenu = NULL;
	new->next = NULL;

	return new;
}

static void DeleteMenus( Menu_entry *menu )

/* deletes all submenus */

{
	if( menu->submenu )
	    DeleteMenus( menu->submenu );

	if( menu->name )
	    free( menu->name );

	if( menu->next )
	    DeleteMenus( menu->next );

	free( menu );
}

static void InitializeMenuColors( void )

{

	/* colors for main menu */

	/* FgMainMenuCol = Black; */
	BgMainMenuCol = Brown;
	BlMainMenuCol = LightGray;
	BdMainMenuCol = DarkGray;

	/* colors for buttons in main menu */

	FgButtonCol = White;
	BgButtonCol = Red;
	BlButtonCol = LightGray;
	BdButtonCol = DarkGray;

	/* colors for pulldown menu */

	FgPulldownCol = Black;
	BgPulldownCol = Green;
	BlPulldownCol = LightGray;
	BdPulldownCol = DarkGray;

}

/****************************************************************/


/********************** public routines *************************/


void RegisterMainMenu( Menu_entry *main )

{
	ASSERT( main != NULL );

	if( MainMenu )
	    Error("Main menu already registered");

	if( main->type != MENU_BAR && main->type != MENU_AREA )
	    Error("Trying to register an erroneous menu");

	InitializeMenuColors();

	MainMenu = main;
}

void DeleteMainMenu( void )

{
	DeleteMenus( MainMenu );
}

Menu_entry *MakeMenuBar( void )

{
	Menu_entry *new;

	new = MakeMenuItem("Menu Bar",MENU_BAR);
	RegisterMainMenu(new);

	return new;
}

Menu_entry *MakeMenuArea( void )

{
	Menu_entry *new;

	new = MakeMenuItem("Menu Area",MENU_AREA);
	RegisterMainMenu(new);

	return new;
}

Menu_entry *MakePulldownMenu( char *name )

{
	return MakeMenuItem(name,MENU_PULLDOWN);
}

Menu_entry *MakeExecMenu( char *name , FP function )

{
	Menu_entry *new;

	new = MakeMenuItem(name,MENU_EXEC);
	new->function = function;

	return new;
}

void AddMenuItem( Menu_entry *parent , Menu_entry *submenu )

{
	Menu_entry *item;

	ASSERT( parent );
	ASSERT( 
		parent->type == MENU_PULLDOWN ||
		parent->type == MENU_BAR ||
		parent->type == MENU_AREA 
	      );

	item = parent->submenu;

	if( item == NULL ) {
	  parent->submenu = submenu;
	} else {
	  while( item->next )
	     item = item->next;
	  item->next = submenu;
	}

	parent->nsubs++;
	submenu->parent = parent;
}

void SetWindowDimension( int xmin , int ymin , int xmax , int ymax )

{
	WinDim.xmin = xmin;
	WinDim.ymin = ymin;
	WinDim.xmax = xmax;
	WinDim.ymax = ymax;

	if( edebug ) printf("SetWindowDimension: %d %d %d %d\n",xmin,ymin,xmax,ymax);
}

void SetMenuDimension( int xmin , int ymin , int xmax , int ymax )

{
	MenDim.xmin = xmin;
	MenDim.ymin = ymin;
	MenDim.xmax = xmax;
	MenDim.ymax = ymax;

	if( edebug ) printf("SetMenuDimension: %d %d %d %d\n",xmin,ymin,xmax,ymax);
}


void DisplayMainMenu( void )

{
	if( !MainMenu )
	    Error("No main menu registered for display");

	if( MainMenu->type == MENU_BAR ) {
	    DisplayMenuBar( MainMenu );
	} else if( MainMenu->type == MENU_AREA ) {
	    DisplayMenuArea( MainMenu );
	} else {
	    Error("Main menu not of expected type");
	}
}

static void DisplayMenuBar( Menu_entry *main )

{
	Error("Not yet ready for menu bar");
}

static void RegisterMenuAreaDimensions( Menu_entry *main )

/* this has to be done only if MenDim changes -> check for it */

{
	int menuwidth, menuheight;
	int menuitems;
	int maxwidth, maxheight, totwidth;
	int buttonwidth, buttonheight;
	int xmargin, ymargin;
	int onespace;
	int x,y,dx,dy;
	Menu_entry *item;

	menuwidth  = MenDim.xmax - MenDim.xmin + 1;
	menuheight = MenDim.ymax - MenDim.ymin + 1;
	menuitems  = main->nsubs;

	ASSERT( menuwidth > 1 );
	ASSERT( menuheight > 1 );
	ASSERT( menuitems > 0 );

	FindMaxTextDimension(MainMenu,&maxwidth,&maxheight,&totwidth);

	/* find good width */

	if( 2 * maxwidth < menuwidth ) {	/* enough space */
	    buttonwidth = 2 * maxwidth;
	} else if( maxwidth < menuwidth ) {	/* space is tight */
	    buttonwidth = maxwidth + (menuwidth-maxwidth) / 2;
	} else {				/* space too low */
	    buttonwidth = menuwidth;
	}
	xmargin = menuwidth - buttonwidth;

	if( edebug ) printf("%d %d %d %d\n",maxwidth,menuwidth,buttonwidth,xmargin);

	/* find acceptable height */

	if( menuitems * 4 * maxheight < menuheight ) {	/* ok */
	    buttonheight = 3 * maxheight;
	    ymargin = maxheight;
	} else if( menuitems * maxheight < menuheight ) { /* space tight */
	    onespace = menuheight / menuitems;
	    ymargin = (onespace - maxheight) / 2;
	    buttonheight = maxheight + ymargin;
	} else {					/* space too low */
	    buttonheight = menuheight / menuitems;
	    ymargin = 0;
	}
	if( edebug ) printf("%d %d %d %d\n",maxheight,menuheight,buttonheight,ymargin);

	/* register dimensions of buttons */

	x = MenDim.xmin + xmargin/2;
	y = MenDim.ymin + ymargin/2;
	dx = buttonwidth;
	dy = buttonheight;

	item = main->submenu;
	while( item ) {
	    (void)  ToIRect( &(item->position) , x , y , x+dx-1 , y+dy-1 );
	    y = y + dy + ymargin;
	    item = item->next;
	}

	main->position = MenDim;
	main->subpos = MenDim;
}

static void PlotMenuArea( Menu_entry *main )

{
	IRect *irect;
	Menu_entry *item;
	int x0,y0;

	PrintMenuInfo(main);

	MouseHide();

	/* fill menu area */

	SViewport( &MenDim );

	SFillRect( &MenDim , BgMainMenuCol );
	SShadeRect( &MenDim , BlMainMenuCol , BdMainMenuCol );

	STextBackground( BgButtonCol );

	/* plot single buttons */

	item = main->submenu;
	while( item ) {
	    irect = &(item->position);

	    SFillRect( irect , BgButtonCol );
	    SShadeRect( irect , BdButtonCol , BlButtonCol );

	    SNewPen(FgButtonCol);
	    SCenterText( irect , item->name , &x0 , &y0 );
	    SText( x0 , y0 , item->name );

	    item = item->next;
	}

	/* reset normal state */

	STextBackground( Black );

	MouseShow();
}

static void DisplayMenuArea( Menu_entry *main )

{
	RegisterMenuAreaDimensions( main );
	PlotMenuArea( main );
}


/************************************************************************/
/************************************************************************/
/************************************************************************/


static int InButton( IRect *irect , int x , int y )

{
	if( x < irect->xmin || x > irect->xmax ) return FALSE;
	if( y < irect->ymin || y > irect->ymax ) return FALSE;

	return TRUE;
}

static Menu_entry *FindMenuChoice( Menu_entry *menu , int x , int y )

{
	Menu_entry *item;

	if( !InButton( &(menu->subpos) , x , y ) )
	    return NULL;

	item = menu->submenu;
	while( item ) {
	    if( InButton( &(item->position) , x , y ) )
		return item;
	    item = item->next;
	}
	return NULL;
}

void RedrawOpenMenus( void )

{
	Menu_entry *menu;
	StackTable menus;
	void *pixmap;

	if( TopMenu == MainMenu ) return;	/* nothing to redraw */

	SViewport( &WinDim );

	menus = MakeStackTable();

	menu = TopMenu;
	while( menu ) {
	    Push(menus,menu);
	    if( edebug ) printf("RedrawOpenMenus (1) : %s\n",menu->name);
	    menu = menu->parent;
	}

	menu = Pop(menus);
	ASSERT( menu == MainMenu );

	while( (menu=Pop(menus)) != NULL ) {
	    if( edebug ) printf("RedrawOpenMenus (2) : %s\n",menu->name);
	    pixmap = menu->pixmap;
	    menu->pixmap = NULL;	/* FIXME  HACK*/
	    HighlightMenu(menu);
	    PlotPulldownMenu(menu);
	    menu->pixmap = pixmap;
	}

	FreeStackTable(menus);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

/*
 *	This section deals with the event input loop
 */


#define	MOUSE_NONE	0
#define	MOUSE_MOVE	1
#define	MOUSE_PRESS	2
#define	MOUSE_RELEASE	3


static void ProcessEvent( QEvent *event , int *button , int *x , int *y )

{
   if( event->type == QButtonPress  && event->button.button == QButtonLeft ) {
          if( event->button.press == QButtonDown ) {
		*button = MOUSE_PRESS;
	  } else {
		*button = MOUSE_RELEASE;
	  }
          *x = event->button.x;
          *y = event->button.y;
   } else if( event->type == QPointerMove ) {
	  *button = MOUSE_MOVE;
          *x = event->button.x;
          *y = event->button.y;
   } else {		/* process other events !!!! */
	  *button = MOUSE_NONE;
          *x = -1;
          *y = -1;
   }
}

void ProcessMenuInput( int x , int y )

{
	int pressed;
	int button_status;
	QEvent event;

	if( edebug ) printf("ProcessMenuInput: ...entering %d %d\n",x,y);

	SViewport( &WinDim );	/* whole window */
	QAddEvent( QPointerMoveMask );

	TopMenu = MainMenu;
	SchedulePulldownMenu(x,y);

	pressed = TRUE;

	do {

	    QNextEvent( &event );
	    ProcessEvent( &event , &button_status , &x , &y );

	    if( edebug ) printf("ProcessMenuInput: button %d %d %d\n"
					,button_status,x,y);

	    if( button_status == MOUSE_MOVE ) {
		if( pressed )
		    AdjustPulldownMenu(x,y);
	    } else if( pressed && button_status == MOUSE_RELEASE ) {
		ExecutePulldownMenu(x,y);
		pressed = FALSE;
	    } else if( !pressed && button_status == MOUSE_PRESS ) {
		SchedulePulldownMenu(x,y);
		pressed = TRUE;
	    } else if( button_status == MOUSE_NONE ) {
		if( event.type == QExpose ) {
			RedrawAll();		/* HACK */
			RedrawOpenMenus();
		}
	    } else {
		CleanupPulldownMenu();
	    }
	} while( pressed || ( TopMenu && TopMenu != MainMenu ) );

	QDeleteEvent( QPointerMoveMask );

	if( edebug ) printf("ProcessMenuInput: ...leaving\n");
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

/*
 *	This section deals with pulldown menus
 */


void AdjustPulldownMenu( int x , int y )

{
	Menu_entry *menu;

	menu = FindPulldownMenu(x,y);

	if( !menu ) return;

	if( edebug ) printf("AdjustPulldownMenu: %s %d %d\n",menu->name,x,y);
	PrintMenuInfo(menu);
/*	PrintSubMenuInfo(menu); */

	if( menu->pixmap ) {	/* already displayed */
		return;
	}

	SchedulePulldownMenu(x,y);
}

void CleanupPulldownMenu( void )

{
	Menu_entry *menu;

	if( edebug ) printf("CleanupPulldownMenu\n");

	menu = TopMenu;
	while( menu ) {
	    ReleasePulldownMenu( menu );
	    menu = menu->parent;
	}
	TopMenu = menu;
}

void ReleasePulldownMenu( Menu_entry *menu )

{
	ASSERT( menu == TopMenu );

	if( !menu ) return;

	if( edebug ) printf("ReleasePulldownMenu: ");
	if( menu->pixmap ) {
	    if( edebug ) printf("%d %d",menu->subpos.xmin,menu->subpos.ymin);
	    QRestorePixels( 
			menu->subpos.xmin, 
			menu->subpos.ymin,
			menu->pixmap
			    );
	    menu->pixmap = NULL;
	}
	UnHighlightMenu( menu );
	if( edebug ) printf("\n");
}

Menu_entry *FindPulldownMenu( int x , int y )

{
	Menu_entry *menu, *chosenmenu = NULL;

	if( !TopMenu ) TopMenu = MainMenu; /* to re-enter the menu area */

	menu = TopMenu;

	while( menu ) {
	    chosenmenu = FindMenuChoice(menu,x,y);
	    if( chosenmenu )
		return chosenmenu;

	    if( InButton( &(menu->position) , x , y ) )
		return menu;

	    ReleasePulldownMenu( menu );

	    menu = menu->parent;
	    TopMenu = menu;
	}

	return menu;
	    
}

void ExecutePulldownMenu( int x , int y )

{
	Menu_entry *menu;

	menu = FindPulldownMenu(x,y);

	if( !menu ) return;
	if( edebug ) printf("ExecutePulldownMenu: %s %d %d\n",menu->name,x,y);

	if( menu->type == MENU_PULLDOWN ) {
	    DisplayPulldownMenu( menu , x , y );
	} else if( menu->type == MENU_EXEC ) {
	    ExecuteCommand( menu );
	    SViewport( &WinDim );	/* whole window */
	    CleanupPulldownMenu();
	}
}

void SchedulePulldownMenu( int x , int y )

{
	Menu_entry *menu;

	if( edebug ) printf("SchedulePulldownMenu: %d %d\n",x,y);

	menu = FindPulldownMenu(x,y);

	PrintMenuInfo(menu);

	if( !menu ) return;

	if( edebug ) printf("SchedulePulldownMenu: %s %d %d\n",menu->name,x,y);

	if( menu->type == MENU_PULLDOWN ) {
	    HighlightMenu( menu );
	    DisplayPulldownMenu( menu , x , y );
	} else {
	    HighlightMenu( menu );
	    TopMenu = menu;
	}
}

void DisplayPulldownMenu( Menu_entry *menu , int x , int y )

{
	int width,height,nmenus;
	int dummy;

	if( edebug ) printf("DisplayPulldownMenu\n");
	if( menu->pixmap )	/* already displayed -> done */
		return;

	nmenus = menu->nsubs;
	if( !nmenus ) return;
	
	TopMenu = menu;

	FindMaxTextDimension(menu,&width,&height,&dummy);
	RegisterPulldownDimensions(menu,&width,&height,&x,&y);
	menu->pixmap = QSavePixels(x,y,width,height);
	if( edebug ) printf("QSavePixels: %d %d %d %d\n",x,y,width,height);

	PlotPulldownMenu( menu );
}

void PlotPulldownMenu( Menu_entry *menu )

{
        IRect *irect;
        Menu_entry *item;
        int x0,y0;

	if( edebug ) printf("PlotPulldownMenu\n");

	MouseHide();

	PrintMenuInfo(menu);

	SFillRect( &(menu->subpos) , BgPulldownCol );
	SShadeRect( &(menu->subpos) , BdPulldownCol , BlPulldownCol );
	STextBackground( BgPulldownCol );
        SNewPen(FgPulldownCol);

	item = menu->submenu;
	while( item ) {
            irect = &(item->position);

            SCenterText( irect , item->name , &x0 , &y0 );
            SText( x0 , y0 , item->name );

            item = item->next;
	}

	STextBackground( Black );
	MouseShow();
}

void SetActFunction( FP fp );	/* HACK */

void ExecuteCommand( Menu_entry *menu )

{
	if( edebug ) printf("Command executed: %s\n",menu->name);
	ExecuteMenuCommand( menu->function );	/* HACK */
}

static void RegisterPulldownDimensions( Menu_entry *menu
		 , int *width , int *height , int *x0 , int *y0 )

{
	int nmenus,winwidth,winheight;
	int w,h,x,y;
	int ws,hs;
	int dx,dy;
	int yp;
	Menu_entry *item;

	/*
		w,h are the total dimension of the pulldown menu (with border)
		ws, hs are the single item dimensions
		on entry :
			width	width of largest text in pulldown
			height	height of largest text in pulldown
			x0,y0	point of mouse press
		on return :
			width	width of total pulldown menu
			height  height of total pulldown menu
			x0,y0	top, left point of menu
	*/

	winwidth = WinDim.xmax - WinDim.xmin + 1;
	winheight = WinDim.ymax - WinDim.ymin + 1;

	nmenus = menu->nsubs;

	w = *width;
	h = *height * nmenus;
	x = *x0;
	y = *y0;

	ASSERT( x <= menu->position.xmax );
	ASSERT( x >= menu->position.xmin );
	ASSERT( y <= menu->position.ymax );
	ASSERT( y >= menu->position.ymin );

	/* find acceptable dimensions */

	if( w * 1.3 < winwidth ) {
	    w *= 1.3;
	} else if( w < winwidth ) {
	    w += (winwidth - w) / 2 ;
	} else {
	    w = winwidth;
	}
	ws = w - 4;		/* 4 pixel for border of pulldown */

	if( h * 2 < winheight ) {
	    h *= 2;
	} else if( h < winheight ) {
	    h += (winheight - h) / 2 ;
	} else {
	    h = winheight;
	}
	hs = (h-4)/nmenus;	/* 4 pixel for border of pulldown */
	h = hs*nmenus + 4;

	/* find origin for display (top/left point) */
	/* for Menu Area try first to put on left side */

	if( MainMenu->type == MENU_AREA && x - 5 - w > WinDim.xmin ) {
	    x = x - 5 - w;
	} else if( x + 5 + w < WinDim.xmax ) {
	    x += 5;
	} else if( x - 5 - w > WinDim.xmin ) {
	    x = x - 5 - w;
	} else {
	    x = WinDim.xmin + (winwidth - w) / 2;
	}
	
	/* this is to be sure that the origin of the pulldown
		menu is always in the button of the parent */

	if( y + 5 < menu->position.ymax ) {
	    yp = y + 5;
	} else if( y - 5 > menu->position.ymin ) {
	    yp = y - 5;
	} else {
	    yp = y;
	}

	if( yp + h < WinDim.ymax ) {
	    y = yp;
	} else if( yp - h > WinDim.ymin ) {
	    y = yp - h;
	} else {
	    y = WinDim.ymin + (winheight - h) / 2;
	}
	
	/* return values */

	*width = w;
	*height = h;
	*x0 = x;
	*y0 = y;

	/* register dimensions of buttons */

	(void) ToIRect( &(menu->subpos) , x , y , x+w-1 , y+h-1 );

	x += 2;		/* border */
	y += 2;
	dx = ws;
	dy = hs;

        item = menu->submenu;
        while( item ) {
            (void)  ToIRect( &(item->position) , x , y , x+dx-1 , y+dy-1 );
            y += dy;
            item = item->next;
        }

}




/*****************************************************************/


static void FindMaxTextDimension( Menu_entry *menu 
			, int *width , int *height , int *totwidth )

{
    int totw=0;
    int maxw=0;
    int maxh=0;
    int w,h;

    menu = menu->submenu;
    while( menu ) {
        STextDimensions(menu->name,&w,&h);
	totw += w;
        maxw = MAX(w,maxw);
        maxh = MAX(h,maxh);
        menu = menu->next;
    }

    *totwidth=totw;
    *width=maxw;
    *height=maxh;
}



void PrintSubMenuInfo( Menu_entry *menu )

{
	Menu_entry *item;

	if( !menu ) {
	    if( edebug ) printf("\nPrintSubMenuInfo: no menu\n");
	    return;
	}

	if( edebug ) printf("\nPrintSubMenuInfo: %s\n",menu->name);

	item = menu->submenu;
	while( item ) {
	    PrintMenuInfo( item );
	    item = item->next;
	}
}

void PrintMenuInfo( Menu_entry *menu )

{
	IRect *p;

	if( !menu ) {
	    if( edebug ) printf("\nPrintMenuInfo: no menu\n");
	    return;
	}

	if( edebug ) printf("\n");
	if( edebug ) printf("PrintMenuInfo: %s %d %d\n",menu->name,menu->type,menu->nsubs);
	p = &(menu->position);
	if( edebug ) printf("position: %d %d %d %d\n",p->xmin,p->ymin,p->xmax,p->ymax);
	p = &(menu->subpos);
	if( edebug ) printf("subpos: %d %d %d %d\n",p->xmin,p->ymin,p->xmax,p->ymax);
	if( edebug ) printf("pixmap: %c   ",menu->pixmap ? 'Y' : 'N');
	if( edebug ) printf("submenu: %c   ",menu->submenu ? 'Y' : 'N');
	if( edebug ) printf("next: %c   ",menu->next ? 'Y' : 'N');
	if( edebug ) printf("function: %c   ",menu->function ? 'Y' : 'N');
	if( edebug ) printf("\n");

}

static void HighlightMenu( Menu_entry *menu )

{

	if( menu->pixmap ) return;	/* already displayed */

	if( menu->parent == NULL )	/* do not paint for main menu */
		return;

	if( menu->parent == MainMenu )	/* do not paint for buttons */
		return;


	SShadeRect( &(menu->position) , BdPulldownCol , BlPulldownCol );
}

static void UnHighlightMenu( Menu_entry *menu )

{
/*	if( !menu->pixmap ) return;*/	/* has not been displayed */

/*
	ASSERT( menu == HighMenu );
*/

	if( menu->parent == NULL )	/* do not paint for main menu */
		return;

	if( menu->parent == MainMenu )	/* do not paint for buttons */
		return;

	SShadeRect( &(menu->position) , BgPulldownCol , BgPulldownCol );
}
