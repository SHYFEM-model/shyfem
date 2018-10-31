
/************************************************************************\
 *
 *    Copyright (C) 1985-2018  Georg Umgiesser
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
 *									*
 * grid_ex.h - extern variables for grid                                *
 *									*
 * Revision History:							*
 * 16-Feb-2011: new options OpOutFile, OpItemType			*
 * 02-Apr-1998: no Button_type, Function_type                           *
 * 09-Feb-1998: ActArgument eliminated                                  *
 * 04-Dec-95: NTotVects, ActVect, HVC added                             *
 * 15-Jan-95: Strings to gridma1 (local)                                *
 * 06-Oct-94: ColTab defined locally in gridop                          *
 * 04-May-94: LINE modes, ...                                           *
 * 13-Apr-94: splitted up in df, ty, fp, ex                             *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRID_EX_
#define __GUH_GRID_EX_


/**************************************************************/
/**************************************************************/
/********************* global variables ***********************/
/**************************************************************/
/**************************************************************/

extern int OpCheck;
extern int OpFill;
extern int OpColor;
extern int OpBorder;
extern int OpShowType;
extern int OpType;
extern int OpArgc;
extern int OpCol;
extern int OpCompress;
extern int OpCheckUse;
extern int OpElim2Lines;
extern float OpMaxColDepth;
extern int OpColTabSize;
extern float OpVectFact;
extern float OpNodeFact;
extern float OpVectScal;
extern float OpNodeScal;
extern char* OpOutFile;
extern int OpItemType;



extern int XPMin,YPMin,XPMax,YPMax;
extern int XMMin,YMMin,XMMax,YMMax;
extern int XSMin,YSMin,XSMax,YSMax;
extern int XCMin,YCMin,XCMax,YCMax;
extern int XTMin,YTMin,XTMax,YTMax;

extern int NTotNodes;
extern int NTotElems;
extern int NTotLines;
extern int NTotVects;

extern Hashtable_type HNN;
extern Hashtable_type HEL;
extern Hashtable_type HLI;
extern Hashtable_type HVC;
extern Queuetable_type CM;

extern Mode_type     ActMode;
extern char         *ActString;

extern float ActX;
extern float ActY;

extern int ShowMode;

extern int Changed;

extern int ActNode;
extern int ActElem;
extern int ActLine;
extern int ActVect;

extern int NTotList;
extern int NodeList[NODELIST_DIM];


/* definitions for colors */

extern int TotCol;
extern int MenCol;
extern int ComCol;
extern int ButCol;
extern int ButChoosCol;
extern int PlotCol    ;
extern int EvidenceCol;
extern int PlotWinCol ;
extern int MesCol     ;
extern int BorderLightCol;
extern int BorderDarkCol ;
extern int ComWriteCol;
extern int MesWriteCol;
extern int ButWriteCol;


extern Rect GbCom ;
extern Rect GbMes ;
extern Rect GbPlo ;
extern Rect GbAll ;
extern Rect GbMen ;
extern Rect GbTot ;
extern Rect GbCon ;


#endif

