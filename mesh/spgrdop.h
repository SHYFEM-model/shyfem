
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
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
 * spgrdop.h - command line options for spgrd
 *
 * revision log :
 *
 * 17.08.1995	ggu	routines copied from meshop.c
 *
\************************************************************************/


#ifndef __GUH_EXGRDOP_
#define __GUH_EXGRDOP_


/**************************************************************/
/**************************************************************/
/********************* global variables ***********************/
/**************************************************************/
/**************************************************************/

extern int OpArgc;

extern int OpInfo;
extern int OpExtractNodes;
extern int OpDeleteNodes;
extern int OpExtractElems;
extern int OpDeleteElems;
extern int OpExtractLines;
extern int OpDeleteLines;
extern int OpExtractUnusedNodes;
extern int OpDeleteUnusedNodes;
extern int OpMinType;
extern int OpMaxType;
extern int OpMinVertex;
extern int OpMaxVertex;
extern int OpUnifyNodes;
extern int OpCompressNumbers;
extern float OpMinDepth;
extern float OpMaxDepth;
extern float OpTollerance;


void SetOptions(int argc, char *argv[]);
void TestVersion( void );


#endif

