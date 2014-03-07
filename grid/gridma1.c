
/************************************************************************\
 *									*
 * gridma1.c - routines used directly by main                           *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
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
 * 18-Feb-2014: new routines GfDelRemoveNodeLine(), GfInsertNodeLine()	*
 * 16-Feb-2011: use OpItemType for new items				*
 * 02-Apr-1998: new functio integrated -> no gridmu.h, no ActCommand    *
 *                call ExitEventLoop() in GfExit to exit loop           *
 *                new GetActFunction(), SetActFunction()                *
 *                MenuFieldInput() substituted by ExecuteMenuCommand()  *
 * 11-Feb-1998: new function TentativeInput(x,y) for keyboard input     *
 * 09-Feb-1998: ActArgument eliminated, new functions GfZoom, GfShow    *
 * 19-Nov-97: write result of keyboard input to message window          *
 * 14-Oct-97: Administer use of nodes with routines                     *
 *            new routines GetAct/SetAct... for incapsulation           *
 * 13-Oct-97: New routine GfRemoveLine() -> removes line with nodes     *
 *            New routine GfRemoveElement() -> removes elem with nodes  *
 * 10-Oct-97: New global strings StringSave, StringUnifyNode1/2         *
 *            New routine GfUnifyNode()                                 *
 *            New functionality in GfChangeDepth(), GfChangeType():     *
 *                     input depth/type from keyboard any time          *
 * 06-Dec-95: IsDegenerateRect() introduced in ZoomWindow               *
 * 05-Dec-95: Arrow Up/Down to scale Node/Vect introduced               *
 *              control special keys from keyboard                      *
 * 04-Dec-95: PlotFieldInput restructured (Tentative...)                *
 * 02-Dec-95: PlotPlot() renamed to PlotAll()                           *
 * 12-Mar-95: In GfMakeElement -> minimum 3 vertices needed             *
 *            In GfMakeLine    -> minimum 2 vertices needed             *
 * 03-Feb-95: local variables to function : ActLine1/2, MoveNode, ActP, *
 *                         ActZoom1/2 ; to file : TentativeXY           *
 * 12-Nov-94: single functions for single tasks                         *
 * 21-Oct-94: Changed introduced                                        *
 * 16-May-94: new file for routines called by main                      *
 * 13-May-94: MakeDepthFromNodes() to compute depth                     *
 * 10-May-94: WriteToMessWin has changed format                         *
 * 05-May-94: ActXYValid removed (useless, use TentativeXY)             *
 * 13-Apr-94: completely restructured                                   *
 *             -> new hash.c for hash routines                          *
 *             -> new list.c for list table routines                    *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
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

/* #include "gridmu.h" */
#include "gridky.h"

#include "menu.h"

/*
static char *StringNone        = "";
*/
static char *StringWorking     = "Working...";
static char *StringZoomWindow1 = "Input first data point to delimit area";
static char *StringZoomWindow2 = "Input second data point to delimit area";
static char *StringZoomIn      = "Enter data point for Zoom In";
static char *StringZoomOut     = "Enter data point for Zoom Out";
static char *StringMove        = "Enter data point to move to";
static char *StringTotalView   = "Total area";
static char *StringCancel      = "Command cancelled";
static char *StringRefresh     = "Screen Refresh";
static char *StringPrint       = "PostScript file written";
static char *StringSave        = "File saved to save.grd";
static char *StringExit        = "Click one more time on Exit to confirm";
static char *StringShowNode    = "Enter data point close to node";
static char *StringShowLine    = "Enter data point close to line";
static char *StringShowElem    = "Enter data point in element";
static char *StringShowVect    = "Enter data point close to vector";
static char *StringMakeLine1   = "Input first vertex of line";
static char *StringMakeLine2   = "Input next verteces (last twice to close)";
static char *StringMakeLineErr = "Error making line : Make new line";
static char *StringMakeElem1   = "Input first vertex of element";
static char *StringMakeElem2   = "Input next verteces (first to close)";
static char *StringMakeElemErr = "Error making element : Make new element";
static char *StringMakeNode    = "Input coordinate of node";
static char *StringDelElem     = "Click on element to delete";
static char *StringRemoveElem  = "Click on element to remove with node";
static char *StringDelLine     = "Click on line to delete";
static char *StringRemoveLine  = "Click on line to remove with nodes";
static char *StringDelNode     = "Click on node to delete";
static char *StringDelVect     = "Click on vector to delete";
static char *StringSplitLine1  = "Click on line to split";
static char *StringSplitLine2  = "Click on node where to split";
static char *StringJoinLine1   = "Click on first line to join";
static char *StringJoinLine2   = "Click on second line to join";
static char *StringJoinLine3   = "Click on node where to join";
static char *StringDelNodeLine1= "Click on line where to delete";
static char *StringDelNodeLine2= "Click on node where to delete";
static char *StringInsertNodeLine1= "Click on line where to insert";
static char *StringInsertNodeLine2= "Input coordinate of node";
static char *StringDelElNFound = "Element not found : Choose element to delete";
static char *StringDelLiNFound = "Line not found : Choose line to delete";
static char *StringDelNdNFound = "Node not found : Choose node to delete";
static char *StringDelVcNFound = "Vector not found : Choose vector to delete";
static char *StringNodeInUse   = "Cannot delete used node : Choose other node";
static char *StringMoveNode1   = "Input node to move";
static char *StringMoveNode2   = "Input new node coordinate";
static char *StringUnifyNode1   = "Input first node to be unified";
static char *StringUnifyNode2   = "Input second node to be unified";
static char *StringChangeVect = "Click on vector to change";
static char *StringChangeElemType1 = "Click on element or enter type";
static char *StringChangeElemType2 = "Click on element to change";
static char *StringChangeNodeType1 = "Click on node or enter type";
static char *StringChangeNodeType2 = "Click on node to change";
static char *StringChangeLineType1 = "Click on line or enter type";
static char *StringChangeLineType2 = "Click on line to change";
static char *StringChangeElemDepth1= "Click on element or enter depth";
static char *StringChangeElemDepth2= "Click on element to change";
static char *StringChangeNodeDepth1= "Click on node or enter depth";
static char *StringChangeNodeDepth2= "Click on node to change";
static char *StringChangeLineDepth1= "Click on line or enter depth";
static char *StringChangeLineDepth2= "Click on line to change";

/**********************************************************************/

static int TentativeXY = FALSE;
static void *ActP = NULL;
static FP ActFunction = NULL;

/**********************************************************************/

static void TentativeNode( float x , float y )

{
	Node_type *pn;

	pn = FindClosestNode(HNN,x,y);
	if( pn && ActP != pn ) {
		ActX=pn->coord.x;
		ActY=pn->coord.y;
		ActP=(void *)pn;
		TentativeXY = TRUE;
		MakeNodeActive(pn->number);
		MoveToPoint(ActX,ActY);
	} else {
		ActX=x;
		ActY=y;
		TentativeXY = FALSE;
		MakeNodeActive(0);
	}
}

static void TentativeElem( float x , float y )

{
	Elem_type *pe;

	pe=FindElemToPoint(HEL,HNN,x,y);
	if( pe && ActP != pe ) {
		MakeGravityPoint(pe,&ActX,&ActY);
		ActP=(void *)pe;
		TentativeXY = TRUE;
		MakeElemActive(pe->number);
		MoveToPoint(ActX,ActY);
	} else {
		ActX=x;
		ActY=y;
		TentativeXY = FALSE;
		MakeElemActive(0);
	}
}

static void TentativeLine( float x , float y )

{
	Line_type *pl;

	pl = FindClosestLine(HLI,HNN,x,y);
	if( pl && ActP != pl ) {
		ActX=x;
		ActY=y;
		ActP=(void *)pl;
		MakeMidPoint(pl,&ActX,&ActY);
		TentativeXY = TRUE;
		MakeLineActive(pl->number);
		MoveToPoint(ActX,ActY);
	} else {
		ActX=x;
		ActY=y;
		TentativeXY = FALSE;
		MakeLineActive(0);
	}
}

static void TentativeVect( float x , float y )

{
	Node_type *pn;

	pn = FindClosestNode(HVC,x,y);
	if( pn && ActP != pn ) {
		ActX=pn->coord.x;
		ActY=pn->coord.y;
		ActP=(void *)pn;
		TentativeXY = TRUE;
		MakeVectActive(pn->number);
		MoveToPoint(ActX,ActY);
	} else {
		ActX=x;
		ActY=y;
		TentativeXY = FALSE;
		MakeVectActive(0);
	}
}

/**********************************************************************/

int GetActNode( void ) { return ActNode; }
void SetActNode( int node ) { ActNode = node; }
int GetActElem( void ) { return ActElem; }
void SetActElem( int elem ) { ActElem = elem; }
int GetActLine( void ) { return ActLine; }
void SetActLine( int line ) { ActLine = line; }
int GetActVect( void ) { return ActVect; }
void SetActVect( int vect ) { ActVect = vect; }

FP GetActFunction( void ) { return ActFunction; }
void SetActFunction( FP fp ) { ActFunction = fp; }

/**********************************************************************/

void TentativeInput( float x , float y )

{
	ActMode = PLOT_FIELD_INPUT;

	if( TentativeXY == FALSE ) /* if other than right button */
		ActP = NULL;	   /* has been pressed the last time */

	if( ShowMode == NODEMODE ) {
		TentativeNode(x,y);
	} else if( ShowMode == ELEMMODE ) {
		TentativeElem(x,y);
	} else if( ShowMode == LINEMODE ) {
		TentativeLine(x,y);
	} else if( ShowMode == VECTMODE ) {
		TentativeVect(x,y);
	}

	if( TentativeXY == FALSE ) ActP=NULL;

	WriteToMesWindow();
}

/**********************************************************************/

void PlotFieldInput( int horiz , int verti , int button )

{
	float x,y;
	FP fp;

	ActMode = PLOT_FIELD_INPUT;

	GetPlotCoord(horiz,verti,&x,&y);

	if( button == RIGHT_MOUSE_BUTTON ) {
		TentativeInput(x,y);
		return;
	} else if( TentativeXY == FALSE ) {
		ActX=x;
		ActY=y;
		if( ShowMode == NODEMODE )
			MakeNodeActive(0);
		else if( ShowMode == ELEMMODE )
			MakeElemActive(0);
		else if( ShowMode == LINEMODE )
			MakeLineActive(0);
		else if( ShowMode == VECTMODE )
			MakeVectActive(0);
	} else {
		/* tentatively 05.12.95 */
		if( ShowMode != VECTMODE )
			TentativeXY=FALSE;
	}

	fp = GetActFunction();
        if( fp )
            (*fp)();

	/*  next is needed  */

	MakeNodeActive(ActNode);
	MakeElemActive(ActElem);
	MakeLineActive(ActLine);
	MakeVectActive(ActVect);

	WriteToMesWindow(); /* ??? */
	WriteToComWindow();

}

void ExecuteMenuCommand( FP fp )

{
        ActMode = MENU_FIELD_INPUT;	/* probably set elsewhere */

        if( fp )
            (*fp)();

	MakeNodeActive(ActNode);
	MakeElemActive(ActElem);
	MakeLineActive(ActLine);
	MakeVectActive(ActVect);

	WriteToComWindow();

	SetActFunction( fp );
}

void KeyboardInput( int c )

{
	int number;
	Node_type *pn=NULL;
	Elem_type *pe=NULL;
	Line_type *pl=NULL;
        FP fp;
	float x,y;
	int cl,cd;

	if( ActMode != KEYBOARD_INPUT ) {
		ActMode = KEYBOARD_INPUT;
		ResetKeyboardInput();
	}

	cl = tolower(c);
        /* FIXME -> otherwise isdigit is not working on Linux systems */
        cd = c >= 0x0F00 ? c - 0x0F00 : c ; 

	if( cl == 'e' || cl == 'q' ) {
		GfExit();
		WriteToComWindow();
	} else if( c == 127 || c == QKeyBackSpace || c == QKeyDelete ) {
		SubKeyboardInput();
		WriteToMesWindow();
	} else if( isdigit(cd) || c == '.' || c == '-' ) {
		AddKeyboardInput( c );
		WriteToMesWindow();
	} else if( c == '\n' || c == '\r' || c == QKeyReturn ) {
	   fp = GetActFunction();
	   if( fp == GfChangeDepth || fp == GfChangeType ) {
        	(*fp)();
		WriteToMesWindow();
		WriteToComWindow();
		ResetKeyboardInput();
	   } else {
		number = GetKeyboardInt();
		ResetKeyboardInput();
		if( ShowMode == NODEMODE ) {
			pn=RetrieveByNodeNumber(HNN,number);
			if( pn ) {
				x=pn->coord.x;
				y=pn->coord.y;
			}
		} else if( ShowMode == ELEMMODE ) {
			pe=RetrieveByElemNumber(HEL,number);
			if( pe ) {
				MakeGravityPoint(pe,&x,&y);
			}
		} else if( ShowMode == LINEMODE ) {
			pl=RetrieveByLineNumber(HLI,number);
			if( pl ) {
				x=ActX; y=ActY;
				MakeMidPoint(pl,&x,&y);
			}
		} else if( ShowMode == VECTMODE ) {
			pn=RetrieveByNodeNumber(HVC,number);
			if( pn ) {
				x=pn->coord.x;
				y=pn->coord.y;
			}
		}
		if( pn || pe || pl ) {
			TentativeInput(x,y);
		}
	   }
	} else if( c == QKeyUp || c == QKeyDown ) {
	    if( ShowMode == NODEMODE ) {
		if( c == QKeyUp ) {
			OpNodeFact *= 2.;
		} else {
			OpNodeFact *= 0.5;
		}
		MakePlotWindow(&GbPlo);
		PlotAll();
	    } else if( ShowMode == VECTMODE ) {
		if( c == QKeyUp ) {
			OpVectFact *= 2.;
		} else {
			OpVectFact *= 0.5;
		}
		MakePlotWindow(&GbPlo);
		PlotAll();
	    }
	}
}

/***************************************************************\
                   	 File menu
\***************************************************************/

void GfCancel( void )

{
     if( ActMode == MENU_FIELD_INPUT ) {
	ActString = StringCancel;
	TentativeXY=FALSE;
	UnActive();
     }
}

void GfRefresh( void )

{
	/* should be reset to the last command executed */

     if( ActMode == MENU_FIELD_INPUT ) {
	RedrawAll();
	ActString = StringRefresh;
	TentativeXY=FALSE;
	UnActive();
     }
}

void GfPrint( void )

{
     if( ActMode == MENU_FIELD_INPUT ) {
	ActString = StringPrint;
	TentativeXY=FALSE;
	WritePS(&GbPlo,HNN,HEL,HLI,CM);
     }
}

void GfSave( void )

{
     if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringSave;
        TentativeXY=FALSE;
	SaveFile();
     }
}

void GfExit( void )

{
     if( ActMode == MENU_FIELD_INPUT ) {
	TentativeXY=FALSE;
	UnActive();
        ActString = StringExit;
	ExitEventLoop();
     } else if( ActMode == KEYBOARD_INPUT ) {
	TentativeXY=FALSE;
	UnActive();
	if( ActString == StringExit ) {
	  ExitEventLoop();
	} else {
          ActString = StringExit;
	}
     }
}

/***************************************************************\
                   	 View menu
\***************************************************************/

static void Zoom( int mode )

{
	static Point ActZoom1,ActZoom2;

     if( ActMode == MENU_FIELD_INPUT ) {

        switch ( mode ) {

            case 0 :
				ActString = StringZoomWindow1;
				TentativeXY=FALSE;
				break;
            case 1 :
				ActString = StringZoomIn;
				TentativeXY=FALSE;
				break;
            case 2 :
				ActString = StringZoomOut;
				TentativeXY=FALSE;
				break;
            case 3 :
				ActString = StringMove;
				TentativeXY=FALSE;
				break;
            case 4 :
				ActString = StringWorking;
				WriteToComWindow();
				GbPlo = GbAll;
				MakePlotWindow(&GbPlo);
				PlotAll();
				ActString = StringTotalView;
				WriteToComWindow();
				break;
        }
     } else {

        switch ( mode ) {

            case 0 :
                if( ActString == StringZoomWindow2 ) {
                    ActZoom2.x = ActX;
                    ActZoom2.y = ActY;
                    MakeRectFromPoints(&ActZoom1,&ActZoom2,&GbPlo);
		    if( !IsDegenerateRect(&GbPlo) ) {
                	MakePlotWindow(&GbPlo);
                	PlotAll();
		    }
                    ActString = StringZoomWindow1;
                } else if( ActString == StringZoomWindow1 ) {
                    ActZoom1.x = ActX;
                    ActZoom1.y = ActY;
                    ActString = StringZoomWindow2;
                }
                break;
            case 1 :
                ZoomInOut(&GbPlo,ActX,ActY,0.5);
                MakePlotWindow(&GbPlo);
                PlotAll();
                break;
            case 2 :
                ZoomInOut(&GbPlo,ActX,ActY,2.0);
                MakePlotWindow(&GbPlo);
                PlotAll();
                break;
            case 3 :
                ZoomInOut(&GbPlo,ActX,ActY,1.0);
                MakePlotWindow(&GbPlo);
                PlotAll();
                break;
        }
     }
}

void GfZoomWindow( void )	{ Zoom(0); }
void GfZoomIn( void )		{ Zoom(1); }
void GfZoomOut( void )		{ Zoom(2); }
void GfMove( void )		{ Zoom(3); }
void GfTotalView( void )	{ Zoom(4); }

void GfMoveRelative( float dx, float dy ) 
{ 
	MoveRelative(&GbPlo, dx , dy ); 
        MakePlotWindow(&GbPlo);
        PlotAll();
}

/***************************************************************\
                   	 Show menu
\***************************************************************/

static void Show( int mode )

{
        switch ( mode ) {

        case 0 :
            ActString = StringShowNode;
	    ShowMode = NODEMODE;
	    UnActive();
	    break;
	case 1 :
	    ActString = StringShowElem;
	    ShowMode = ELEMMODE;
	    UnActive();
	    break;
	case 2 :
	    ActString = StringShowLine;
	    ShowMode = LINEMODE;
            UnActive();
            break;
	case 3 :
	    ActString = StringShowVect;
	    ShowMode = VECTMODE;
            UnActive();
            break;
        }
}

void GfShowNode( void )		{ Show(0); }
void GfShowElement( void )	{ Show(1); }
void GfShowLine( void )		{ Show(2); }
void GfShowVect( void )		{ Show(3); }

/***************************************************************\
                   	 Node menu
\***************************************************************/

void GfMakeNode( void )

{
    Point c;
    Node_type *pn;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringMakeNode;
        ShowMode = NODEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        c.x = ActX;
        c.y = ActY;
        NTotNodes++;
        pn = MakeNode(NTotNodes,OpItemType,&c);
        InsertByNodeNumber(HNN,pn);
        ActNode = NTotNodes;
    }
    Changed = TRUE;
}

void GfDelNode( void )

{
    Node_type *pn;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringDelNode;
        ShowMode = NODEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActNode ) {
            pn=FindNode(HNN,ActNode);
            if( GetUseN(pn) ) {
                ActString = StringNodeInUse;
            } else {
                EvidenceNode(ActNode,PlotWinCol);
                DeleteNode(pn);
                ActNode=0;
                ActString = StringDelNode;
            }
        } else {
            ActString = StringDelNdNFound;
        }
    }
    Changed = TRUE;
}

void GfMoveNode( void )

{
    Node_type *pn;
	static int MoveNode=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringMoveNode1;
        ShowMode = NODEMODE;
        TentativeXY=FALSE;
    } else {
        if( ActString == StringMoveNode1 ) {
            if( ActNode != 0 ) {
                MoveNode = ActNode;
                ActString = StringMoveNode2;
            } else {
                MoveNode = 0;
            }
        } else if( ActString == StringMoveNode2 ) {
            if( MoveNode != 0 ) {
                ActNode = MoveNode;
                pn=FindNode(HNN,ActNode);
                if( pn ) {
                    EvidenceNode(ActNode,PlotWinCol);
                    pn->coord.x = ActX;
                    pn->coord.y = ActY;
                } else {
                    Error2("No node ",itos(MoveNode));
                }
                ActString = StringMoveNode1;
            }
        }
    }
    Changed = TRUE;
}

void GfUnifyNode( void )

{
    Node_type *pna,*pnu;
    static int UniqueNode=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringUnifyNode1;
        ShowMode = NODEMODE;
        TentativeXY=FALSE;
    } else {
        if( ActString == StringUnifyNode1 ) {
            if( ActNode != 0 ) {
                UniqueNode = ActNode;
                ActString = StringUnifyNode2;
            } else {
                UniqueNode = 0;
            }
        } else if( ActString == StringUnifyNode2 ) {
            if( UniqueNode != 0 && ActNode != 0 && UniqueNode != ActNode) {
                pna=FindNode(HNN,ActNode);
                pnu=FindNode(HNN,UniqueNode);
                if( pna && pnu ) {
                    EvidenceNode(UniqueNode,PlotWinCol);
		    SubstituteNode(pna,pnu);
		    DeleteNode(pnu);
                } else {
		    if( !pna ) {
                        Error2("No node ",itos(ActNode));
		    } else {
                        Error2("No node ",itos(UniqueNode));
		    }
                }
                ActString = StringUnifyNode1;
            }
        }
    }
    Changed = TRUE;
}

/***************************************************************\
                   	 Element menu
\***************************************************************/

void GfMakeElement( void )

{
	Elem_type *pe;
	Node_type *pn;
	Point c;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringMakeElem1;
        NTotList=0;
        ShowMode = NODEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActNode != 0 ) {
                if( NTotList == 0 || ActNode != NodeList[0] ) {
                    NodeList[NTotList++] = ActNode;
                } else {  /* make element */
		  if( NTotList < 3 ) { /* less than 3 vertices */
		    NTotList=0;
                    ActString = StringMakeElemErr;
		  } else {		/* ok -> make */
                    NTotElems++;
                    pe = MakeElem(NTotElems,OpItemType,NodeList,NTotList);
                    InsertByElemNumber(HEL,pe);
                    AddUseE(HNN,pe);
                    MakeElemActive(NTotElems);
                    NTotList=0;
                    ActString = StringMakeElem1;
		  }
                }
        } else { /* make node first */
		c.x = ActX;
		c.y = ActY;
		NTotNodes++;
		pn = MakeNode(NTotNodes,OpItemType,&c);
		InsertByNodeNumber(HNN,pn);
		ActNode = NTotNodes;
		if( NTotList == 0 )
               	    ActString = StringMakeElem2;
                NodeList[NTotList++] = ActNode;
	}
        if( NTotList == NODELIST_DIM ) { /* error - make better */
             NTotList=0;
             ActString = StringMakeElem1;
        }
	if( NTotList > 0 )
		ActString = StringMakeElem2;
	Changed = TRUE;
     }
}

void GfDelElement( void )

{
	Elem_type *pe;

     if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringDelElem;
        ShowMode = ELEMMODE;
        UnActive();
        TentativeXY=FALSE;
     } else {
        if( ActElem ) {
            EvidenceElem(ActElem,PlotWinCol);
            pe=FindElem(HEL,ActElem);
	    DeleteUseE(HNN,pe);
            DeleteElem(pe);
            ActElem=0;
            ActString = StringDelElem;
	    CheckUseConsistency();	/* HACK ggu */
        } else {
            ActString = StringDelElNFound;
        }
     }
    Changed = TRUE;
}

void GfRemoveElement( void )

{
	Elem_type *pe;

     if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringRemoveElem;
        ShowMode = ELEMMODE;
        UnActive();
        TentativeXY=FALSE;
     } else {
        if( ActElem ) {
            EvidenceElem(ActElem,PlotWinCol);
            pe=FindElem(HEL,ActElem);
	    DeleteUseE(HNN,pe);
            DeleteElemWithNodes(pe);
            ActElem=0;
            ActString = StringRemoveElem;
	    CheckUseConsistency();	/* HACK ggu */
        } else {
            ActString = StringDelElNFound;
        }
     }
    Changed = TRUE;
}

/***************************************************************\
                   	 Line menu
\***************************************************************/

void GfMakeLine( void )

{
    Line_type *pl;
    Node_type *pn;
    Point c;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringMakeLine1;
        NTotList=0;
        ShowMode = NODEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
	if( ActNode != 0 ) {
                if( NTotList == 0 || ActNode != NodeList[NTotList-1] ) {
                    NodeList[NTotList++] = ActNode;
                } else {  /* make line */
		  if( NTotList < 2 ) { /* less than 2 vertices */
		    NTotList=0;
                    ActString = StringMakeLineErr;
		  } else {		/* ok -> make */
                    NTotLines++;
                    pl = MakeLine(NTotLines,OpItemType,NodeList
                            ,NTotList);
                    InsertByLineNumber(HLI,pl);
                    AddUseL(HNN,pl);
                    MakeLineActive(NTotLines);
                    NTotList=0;
                    ActString = StringMakeLine1;
		  }
                }
        } else { /* make node first */
		c.x = ActX;
		c.y = ActY;
		NTotNodes++;
		pn = MakeNode(NTotNodes,OpItemType,&c);
		InsertByNodeNumber(HNN,pn);
		ActNode = NTotNodes;
                NodeList[NTotList++] = ActNode;
	}
        if( NTotList == NODELIST_DIM ) { /* error - make better */
            NTotList=0;
            ActString = StringMakeLine1;
        }
	if( NTotList > 0 )
		ActString = StringMakeLine2;
	Changed = TRUE;
    }
}

void GfDelLine( void )

{
    Line_type *pl;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringDelLine;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActLine ) {
            EvidenceLine(ActLine,PlotWinCol);
            pl=FindLine(HLI,ActLine);
	    DeleteUseL(HNN,pl);
            DeleteLine(pl);
            ActLine=0;
            ActString = StringDelLine;
	    CheckUseConsistency();	/* HACK ggu */
        } else {
            ActString = StringDelLiNFound;
        }
    }
    Changed = TRUE;
}

void GfRemoveLine( void )

{
    Line_type *pl;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringRemoveLine;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActLine ) {
            EvidenceLine(ActLine,PlotWinCol);
            pl=FindLine(HLI,ActLine);
	    DeleteUseL(HNN,pl);
            DeleteLineWithNodes(pl);
            ActLine=0;
            ActString = StringRemoveLine;
	    CheckUseConsistency();	/* HACK ggu */
        } else {
            ActString = StringDelLiNFound;
        }
    }
    Changed = TRUE;
}

void GfJoinLine( void )

{
    Line_type *pl1,*pl2;
    static int ActLine1=0;
    static int ActLine2=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringJoinLine1;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActString == StringJoinLine1 ) {
            if( ActLine != 0 ) {
                ActLine1=ActLine;
                ActString = StringJoinLine2;
            }
        } else if( ActString == StringJoinLine2 ) {
            if( ActLine != 0 ) {
                ActLine2=ActLine;
                ActString = StringJoinLine3;
                ShowMode = NODEMODE;
            }
        } else if( ActString == StringJoinLine3 ) {
            if( ActNode != 0 ) {
                pl1=FindLine(HLI,ActLine1);
                pl2=FindLine(HLI,ActLine2);
                JoinLine(HLI,pl1,pl2,ActNode);
                ActString = StringJoinLine1;
                ShowMode = LINEMODE;
		CheckUseConsistency();	/* HACK ggu */
            }
        }
    }
    Changed = TRUE;
}

void GfSplitLine( void )

{
    Line_type *pl;
    static int ActLine1=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringSplitLine1;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActString == StringSplitLine1 ) {
            if( ActLine != 0 ) {
                ActLine1=ActLine;
                ActString = StringSplitLine2;
                ShowMode = NODEMODE;
            }
        } else if( ActString == StringSplitLine2 ) {
            if( ActNode != 0 ) {
                pl=FindLine(HLI,ActLine1);
                SplitLine(HLI,pl,ActNode);
                ActString = StringSplitLine1;
                ShowMode = LINEMODE;
		CheckUseConsistency();	/* HACK ggu */
            }
        }
    }
    Changed = TRUE;
}

static void GfDelRemoveNodeLine( int mode );

void GfDelNodeLine( void )
{
     GfDelRemoveNodeLine( 0 );
}

void GfRemoveNodeLine( void )
{
     GfDelRemoveNodeLine( 1 );
}

static void GfDelRemoveNodeLine( int delete_node )

{
    Line_type *pl;
    Node_type *pn;
    static int ActLine1=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringDelNodeLine1;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActString == StringDelNodeLine1 ) {
            if( ActLine != 0 ) {
                ActLine1=ActLine;
                ActString = StringDelNodeLine2;
                ShowMode = NODEMODE;
            }
        } else if( ActString == StringDelNodeLine2 ) {
            if( ActNode != 0 ) {
                pl=FindLine(HLI,ActLine1);
		EvidenceLine(ActLine1,PlotWinCol);
                DelNodeLine(HLI,pl,ActNode);
		if( delete_node == 1 ) {
                  pn=FindNode(HNN,ActNode);
                  DeleteNode(pn);
		} else {
                  EvidenceNode(ActNode,PlotCol);
		}
                ActNode=0;
		EvidenceLine(ActLine1,PlotCol);
                ActString = StringDelNodeLine1;
                ShowMode = LINEMODE;
		CheckUseConsistency();	/* HACK ggu */
            }
        }
    }
    Changed = TRUE;
}

void GfInsertNodeLine( void )
{
    Line_type *pl;
    Node_type *pn;
    Point c;
    static int ActLine1=0;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringInsertNodeLine1;
        ShowMode = LINEMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActString == StringInsertNodeLine1 ) {
            if( ActLine != 0 ) {
                ActLine1=ActLine;
                ActString = StringInsertNodeLine2;
                ShowMode = NODEMODE;
            }
        } else if( ActString == StringInsertNodeLine2 ) {
            if( ActNode == 0 ) {
        	c.x = ActX;
		c.y = ActY;
		NTotNodes++;
		pn = MakeNode(NTotNodes,OpItemType,&c);
		InsertByNodeNumber(HNN,pn);
		ActNode = NTotNodes;
	    }
            pl=FindLine(HLI,ActLine1);
	    EvidenceLine(ActLine1,PlotWinCol);
            InsertNodeLine(HLI,pl,ActNode);
            ActNode=0;
            ActLine=0;
	    EvidenceLine(ActLine1,PlotCol);
            ActString = StringInsertNodeLine1;
            ShowMode = LINEMODE;
	    CheckUseConsistency();	/* HACK ggu */
        }
    }
    Changed = TRUE;
}

/***************************************************************\
                   	 Vector menu
\***************************************************************/

void GfDelVect( void )

{
    Node_type *pn;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringDelVect;
        ShowMode = VECTMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActVect ) {
            pn=FindNode(HVC,ActVect);
            EvidenceVect(ActVect,PlotWinCol);
            DeleteVect(pn);
            ActVect=0;
            ActString = StringDelVect;
        } else {
            ActString = StringDelVcNFound;
        }
    }
    Changed = TRUE;
}

void GfChangeVect( void )

{
    Node_type *pn;

    if( ActMode == MENU_FIELD_INPUT ) {
        ActString = StringChangeVect;
        ShowMode = VECTMODE;
        UnActive();
        TentativeXY=FALSE;
    } else {
        if( ActVect ) {
            EvidenceVect(ActVect,PlotWinCol);
            pn=FindNode(HVC,ActVect);
            ChangeVect((Vect_type *)pn->extra);
            ActString = StringChangeVect;
        }
    }
    Changed = TRUE;
}

/***************************************************************\
                   	 Change menu
\***************************************************************/

void GfChangeDepth( void )

{
	Elem_type *pe;
	Node_type *pn;
	Line_type *pl;
	static float ActDepth = NULLDEPTH;
	char s[60];

    if( ActMode == MENU_FIELD_INPUT ) {
	if( ShowMode == ELEMMODE )
        	ActString = StringChangeElemDepth1;
	else if( ShowMode == NODEMODE )
        	ActString = StringChangeNodeDepth1;
	else if( ShowMode == LINEMODE )
        	ActString = StringChangeLineDepth1;
        TentativeXY=FALSE;
    } else {
        if( ActString == StringChangeLineDepth1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
                ActString = StringChangeLineDepth2;
            } else if( ActLine != 0 ) {
                pl = FindLine(HLI,ActLine);
                ActDepth = pl->depth;
                ActString = StringChangeLineDepth2;
            }
        } else if( ActString == StringChangeLineDepth2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new depth */
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
            } else if( ActLine != 0 ) {
                pl = FindLine(HLI,ActLine);
                pl->depth = ActDepth;
		Changed = TRUE;
            }
        } else if( ActString == StringChangeElemDepth1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
                ActString = StringChangeElemDepth2;
            } else if( ActElem != 0 ) {
                pe = FindElem(HEL,ActElem);
                ActDepth = pe->depth;
                ActString = StringChangeElemDepth2;
            }
        } else if( ActString == StringChangeElemDepth2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new depth */
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
            } else if( ActElem != 0 ) {
                pe = FindElem(HEL,ActElem);
                pe->depth = ActDepth;
		Changed = TRUE;
            }
        } else if( ActString == StringChangeNodeDepth1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
                ActString = StringChangeNodeDepth2;
            } else if( ActNode != 0 ) {
                pn = FindNode(HNN,ActNode);
                ActDepth = pn->depth;
                ActString = StringChangeNodeDepth2;
            }
        } else if( ActString == StringChangeNodeDepth2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new depth */
		ActDepth = GetKeyboardFloat();
		sprintf(s,"Actual depth: %f",ActDepth);
		SetKeyboardString(s);
            } else if( ActNode != 0 ) {
                pn = FindNode(HNN,ActNode);
                pn->depth = ActDepth;
		Changed = TRUE;
            }
        }
    }
}

void GfChangeType( void )

{
	Elem_type *pe;
	Node_type *pn;
	Line_type *pl;
	static int ActType=0;
	char s[60];

    if( ActMode == MENU_FIELD_INPUT ) {
	if( ShowMode == ELEMMODE )
        	ActString = StringChangeElemType1;
	else if( ShowMode == NODEMODE )
        	ActString = StringChangeNodeType1;
	else if( ShowMode == LINEMODE )
        	ActString = StringChangeLineType1;
        TentativeXY=FALSE;
    } else {
        if( ActString == StringChangeElemType1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
                ActString = StringChangeElemType2;
            } else if( ActElem != 0 ) {
                pe = FindElem(HEL,ActElem);
                ActType = pe->type;
                ActString = StringChangeElemType2;
            }
        } else if( ActString == StringChangeElemType2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new type */
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
            } else if( ActElem != 0 ) {
                pe = FindElem(HEL,ActElem);
                pe->type = ActType;
		Changed = TRUE;
            }
        } else if( ActString == StringChangeNodeType1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
                ActString = StringChangeNodeType2;
            } else if( ActNode != 0 ) {
                pn = FindNode(HNN,ActNode);
                ActType = pn->type;
                ActString = StringChangeNodeType2;
            }
        } else if( ActString == StringChangeNodeType2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new type */
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
            } else if( ActNode != 0 ) {
                pn = FindNode(HNN,ActNode);
                pn->type = ActType;
		Changed = TRUE;
            }
        } else if( ActString == StringChangeLineType1 ) {
	    if( ActMode == KEYBOARD_INPUT ) {
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
                ActString = StringChangeLineType2;
            } else if( ActLine != 0 ) {
                pl = FindLine(HLI,ActLine);
                ActType = pl->type;
                ActString = StringChangeLineType2;
            }
        } else if( ActString == StringChangeLineType2 ) {
	    if( ActMode == KEYBOARD_INPUT ) {	/* set new type */
		ActType = GetKeyboardInt();
		sprintf(s,"Actual type: %d",ActType);
		SetKeyboardString(s);
            } else if( ActLine != 0 ) {
                pl = FindLine(HLI,ActLine);
                pl->type = ActType;
		Changed = TRUE;
            }
        }
    }
}
/***************************************************************\
                   	 Line menu
\***************************************************************/

