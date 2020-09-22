
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,1997-1998,2003,2011-2012  Georg Umgiesser
 *    Copyright (C) 2014  Georg Umgiesser
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
 * grid.c - finite element grid manipulation routines under X11
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	completely restructured
 * ...		ggu	-> new hash.c for hash routines
 * ...		ggu	-> new list.c for list table routines
 * 05.05.1994	ggu	ActXYValid removed (useless, use TentativeXY)
 * 10.05.1994	ggu	WriteToMessWin has changed format
 * 13.05.1994	ggu	MakeDepthFromNodes() to compute depth
 * 06.10.1994	ggu	ColTab defined locally in gridop
 * 21.10.1994	ggu	Changed introduced
 * 15.01.1995	ggu	Strings to gridma1
 * 03.02.1995	ggu	local variables to function : ActLine1/2, MoveNode, ActP,
 * ...		ggu	ActZoom1/2 ; to file gridma1: TentativeXY
 * 05.12.1995	ggu	adjusted for new vector routines
 * 06.12.1995	ggu	new routines to compute min/max of rectangle
 * 19.09.1997	ggu	Call to CheckNodes changed
 * 09.02.1998	ggu	ActArgument eliminated
 * 02.04.1998	ggu	new function MakeGridMenu() (temporarily here)
 * 02.04.1998	ggu	no Button_type, Function_type
 * 13.05.2003	ggu	main changed to menubar in routine MakeGridMenu()
 * 16.02.2011	ggu	new options OpOutFile, OpItemType
 * 01.01.2012	ggu	delete vector item from menu
 * 18.02.2014	ggu	new menu items for line (del node, remove node, insert)
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "grid.h"
#include "graph.h"

#include "general.h"
#include "gustd.h"
#include "keybd.h"
#include "mouse.h"
#include "events.h"

#include "list.h"
#include "hash.h"
#include "gridhs.h"
#include "menu.h"


int OpCheck;
int OpFill;
int OpColor;
int OpBorder;
int OpShowType;
int OpType;
int OpArgc;
int OpCol;
int OpCompress;
int OpElim2Lines;
int OpCheckUse;
float OpMaxColDepth;
int OpColTabSize;
float OpVectFact;
float OpNodeFact;
float OpVectScal;
float OpNodeScal;
char* OpOutFile;
int OpItemType;


/* definitions for colors */

int TotCol;
int MenCol;
int ComCol;
int ButCol;
int ButChoosCol;
int PlotCol    ;
int EvidenceCol;
int PlotWinCol ;
int MesCol     ;
int BorderLightCol;
int BorderDarkCol ;
int ComWriteCol;
int MesWriteCol;
int ButWriteCol;


int XPMin=XPMIN,YPMin=YPMIN,XPMax=XPMAX,YPMax=YPMAX;
int XMMin=XMMIN,YMMin=YMMIN,XMMax=XMMAX,YMMax=YMMAX;
int XSMin=XSMIN,YSMin=YSMIN,XSMax=XSMAX,YSMax=YSMAX;
int XCMin=XCMIN,YCMin=YCMIN,XCMax=XCMAX,YCMax=YCMAX;
int XTMin=XTMIN,YTMin=YTMIN,XTMax=XTMAX,YTMax=YTMAX;



int   NTotNodes  = 0;       /* number of total nodes */
int   NTotElems  = 0;       /* number of total elements */
int   NTotLines  = 0;       /* number of total lines */
int   NTotVects  = 0;       /* number of total vectors */

Hashtable_type HNN;
Hashtable_type HEL;
Hashtable_type HLI;
Hashtable_type HVC;
Queuetable_type CM;		/* list of comments */

Mode_type     ActMode=NO_INPUT;
char         *ActString="";

float ActX;
float ActY;

int ShowMode=NODEMODE;

int Changed=FALSE;

int ActNode=0;
int ActElem=0;
int ActLine=0;
int ActVect=0;

int NTotList;
int NodeList[NODELIST_DIM];


Rect GbCom = { {0.,0.} , {1.,1.} };
Rect GbMes = { {0.,0.} , {1.,1.} };
Rect GbPlo = { {0.,0.} , {1.,1.} }; /* actual plot window */
Rect GbAll = { {0.,0.} , {1.,1.} };     /* total area of elements */
Rect GbMen = { {0.,0.} , {1.,1.} };
Rect GbTot = { {0.,0.} , {1.,1.} };
Rect GbCon = { {0.,0.} , {1.,1.} }; /* constant value, is not changed */

void MakeGridMenu( void );

int main(int argc, char *argv[])

{
	int w,h;
	Rect gb;

	HNN=MakeHashTable();
	HEL=MakeHashTable();
	HLI=MakeHashTable();
	HVC=MakeHashTable();
	CM=MakeQueueTable();

	SetOptions(argc,argv);

	ReadFiles(argc,argv);

	CheckNodes(HNN,HEL,HLI,HVC);
	if( OpCheck )
		CheckMore(HNN,HEL,HLI);

	InitMinMax(&gb);
	GetNodeMinMax(HNN,&gb);
	GetNodeMinMax(HVC,&gb);
	SetMinMax(&gb);

	printf("Min/Max values : %f %f %f %f\n",
		gb.low.x,gb.low.y,gb.high.x,gb.high.y);
	printf("Center values  : %f %f \n",
		(gb.low.x+gb.high.x)/2.,(gb.low.y+gb.high.y)/2.);

	GbAll = GbPlo = gb;

	QGraphInit();
	QNewPen(PlotCol);

	SetColors(OpCol);
	QWindowMaxXY(&w,&h);
	ResizeWindow(w,h);
	ScaleFactor(HVC);

	MakeGridMenu();

/*	MouseInit();	*/
	QInitEvent();

	LoopForInput();

	QGraphClose();

	WriteFiles();

	return 0;
}

void MakeGridMenu( void )

{
	Menu_entry *menubar;
	Menu_entry *file, *view, *show, *node, *elem;
	Menu_entry *line, *change;
/*
	Menu_entry *vect;
	Menu_entry *test, *sub;
*/

/*
	In order to add new functionality the following has to be done:

	- insert new menu item here below
	- make new prototype in grid_fp.h
	- add new routine to gridma1.c
*/


	menubar = MakeMenuArea();



	file = MakePulldownMenu("File");
	AddMenuItem( menubar , file );

	AddMenuItem( file , MakeExecMenu("Cancel",GfCancel) );
	AddMenuItem( file , MakeExecMenu("Refresh",GfRefresh) );
	AddMenuItem( file , MakeExecMenu("Print",GfPrint) );
	AddMenuItem( file , MakeExecMenu("Save",GfSave) );
	AddMenuItem( file , MakeExecMenu("Exit",GfExit) );



	view = MakePulldownMenu("View");
	AddMenuItem( menubar , view );

	AddMenuItem( view , MakeExecMenu("Zoom Window",GfZoomWindow) );
	AddMenuItem( view , MakeExecMenu("Zoom in",GfZoomIn) );
	AddMenuItem( view , MakeExecMenu("Zoom out",GfZoomOut) );
	AddMenuItem( view , MakeExecMenu("Total view",GfTotalView) );
	AddMenuItem( view , MakeExecMenu("Move",GfMove) );



	show = MakePulldownMenu("Show");
	AddMenuItem( menubar , show );

	AddMenuItem( show , MakeExecMenu("Show Node",GfShowNode) );
	AddMenuItem( show , MakeExecMenu("Show Element",GfShowElement) );
	AddMenuItem( show , MakeExecMenu("Show Line",GfShowLine) );
	/* AddMenuItem( show , MakeExecMenu("Show Vector",GfShowVect) ); */



	node = MakePulldownMenu("Node");
	AddMenuItem( menubar , node );

	AddMenuItem( node , MakeExecMenu("Make Node",GfMakeNode) );
	AddMenuItem( node , MakeExecMenu("Del Node",GfDelNode) );
	AddMenuItem( node , MakeExecMenu("Move Node",GfMoveNode) );
	AddMenuItem( node , MakeExecMenu("Unify Node",GfUnifyNode) );



	elem = MakePulldownMenu("Element");
	AddMenuItem( menubar , elem );

	AddMenuItem( elem , MakeExecMenu("Make Element",GfMakeElement) );
	AddMenuItem( elem , MakeExecMenu("Del Element",GfDelElement) );
	AddMenuItem( elem , MakeExecMenu("Remove Element",GfRemoveElement) );



	line = MakePulldownMenu("Line");
	AddMenuItem( menubar , line );

	AddMenuItem( line , MakeExecMenu("Make Line",GfMakeLine) );
	AddMenuItem( line , MakeExecMenu("Del Line",GfDelLine) );
	AddMenuItem( line , MakeExecMenu("Remove Line",GfRemoveLine) );
	AddMenuItem( line , MakeExecMenu("Split Line",GfSplitLine) );
	AddMenuItem( line , MakeExecMenu("Join Line",GfJoinLine) );
	AddMenuItem( line , MakeExecMenu("Del Node",GfDelNodeLine) );
	AddMenuItem( line , MakeExecMenu("Remove Node",GfRemoveNodeLine) );
	AddMenuItem( line , MakeExecMenu("Insert Node",GfInsertNodeLine) );


/*
	vect = MakePulldownMenu("Vector");
	AddMenuItem( menubar , vect );

	AddMenuItem( vect , MakeExecMenu("Del Vector",GfDelVect) );
	AddMenuItem( vect , MakeExecMenu("Change Vector",GfChangeVect) );
*/


	change = MakePulldownMenu("Change");
	AddMenuItem( menubar , change );

	AddMenuItem( change , MakeExecMenu("Change Depth",GfChangeDepth) );
	AddMenuItem( change , MakeExecMenu("Change Type",GfChangeType) );


/*
	test = MakePulldownMenu("Test");
	AddMenuItem( menubar , test );

	AddMenuItem( test , MakeExecMenu("Make Node",GfMakeNode) );
	AddMenuItem( test , MakeExecMenu("Del Node",GfDelNode) );

	    sub = MakePulldownMenu("test view");
	    AddMenuItem( test , sub );

	    AddMenuItem( sub , MakeExecMenu("Zoom Window",GfZoomWindow) );
	    AddMenuItem( sub , MakeExecMenu("Zoom in",GfZoomIn) );
	    AddMenuItem( sub , MakeExecMenu("Zoom out",GfZoomOut) );
	    AddMenuItem( sub , MakeExecMenu("Total view",GfTotalView) );
	    AddMenuItem( sub , MakeExecMenu("Move",GfMove) );

	    sub = MakePulldownMenu("test file");
	    AddMenuItem( test , sub );

	    AddMenuItem( sub , MakeExecMenu("Cancel",GfCancel) );
	    AddMenuItem( sub , MakeExecMenu("Refresh",GfRefresh) );
	    AddMenuItem( sub , MakeExecMenu("Print",GfPrint) );
	    AddMenuItem( sub , MakeExecMenu("Save",GfSave) );
	    AddMenuItem( sub , MakeExecMenu("Exit",GfExit) );

	AddMenuItem( test , MakeExecMenu("Move Node",GfMoveNode) );
	AddMenuItem( test , MakeExecMenu("Unify Node",GfUnifyNode) );
*/

}
