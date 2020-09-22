
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
 * meshop.h - command line options for mesh
 *
 * revision log :
 *
 * 02.08.1995	ggu	routines copied from gridop.c
 *
\************************************************************************/


#ifndef __GUH_MESHOP_
#define __GUH_MESHOP_


/**************************************************************/
/**************************************************************/
/********************* global variables ***********************/
/**************************************************************/
/**************************************************************/

extern int OpArgc;

extern int OpBackGround;
extern int OpIntern;
extern int OpBoundary;
extern int OpRefineBoundary;
extern int OpInsertIntern;
extern int OpCompress;
extern int OpSmoothPass;
extern float OpSmoothOmega;
extern float OpAspect;
extern float OpResolution;


void SetOptions(int argc, char *argv[]);
void TestVersion( void );


#endif

