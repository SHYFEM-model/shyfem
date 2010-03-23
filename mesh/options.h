
/************************************************************************\ 
 *									*
 * options.h - command line options                                     *
 *									*
 * Copyright (c) 1995 by Georg Umgiesser				*
 *									*
 * see options.c for copying information				*
 *									*
 * Revision History:							*
 * 17-Aug-95: routines copied from meshop.h                             *
 *									*
\************************************************************************/


#ifndef __GUH_OPTIONS_
#define __GUH_OPTIONS_


int GetOptInd( void );
char *GetOptArg( void ); 
int GetOpt(int argc, char *argv[], char *optionS);


#endif

