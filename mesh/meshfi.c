
/************************************************************************\
 *
 *    Copyright (C) 1995,1997,2016  Georg Umgiesser
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
 * meshfi.c - file handling routines for mesh
 *
 * revision log :
 *
 * 27.07.1995	ggu	routines written from scratch
 * 01.08.1995	ggu	writing routines transfered from gridck.c
 * 08.10.1997	ggu	routine ReadBnd deleted
 * ...		ggu	if depth is given -> write it
 * ...		ggu	uses mesh type to decide which elements to write
 * 15.10.1997	ggu	in WriteFile: write also lines with type L_NONE
 * 19.03.2016	ggu	program now accepts filename with and without .grd
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fund.h"
#include "gustd.h"

#include "mesh.h"
#include "meshut.h"
#include "meshhs.h"
#include "meshfi.h"
#include "meshop.h"
#include "meshgd.h"
#include "meshty.h"



/**************************************************************************/

/*
static FILE *FpFile=NULL;
static int   FpLine=0;

static int openfile( char *name , char *mode )

{
        FpLine=0;
        FpFile=fopen(name,mode);
        return FpFile ? 1 : 0;
}

static char *getnewline( void )

{
        char *s;

        s = getlin(FpFile);
        if( s )
                FpLine++;
        return s;
}

static void closefile( void ) { if(FpFile) fclose(FpFile); }
*/

/*
static int gettotlines( void ) { return FpLine; }
*/

/**************************************************************************/


void ReadFiles( int argc , char *argv[] )

{
	char sfile[80];
	char *s;

	while( OpArgc < argc ) {
		s=strcpy(sfile,argv[OpArgc++]);
		/* s=strcat(s,".grd"); */
		ReadStandard(s,HNN,HEL,HLI,CM); 
	}
}

void WriteAll( char *file , NodeList list )

{
	WriteFile( file , list , 1 );
}

void WriteGrd( char *file )

{
	WriteFile( file , NULL , 0 );
}

void WriteFile( char *file , NodeList list , int all )

{
	int i,j;
	int count;
	int linetype;
	FILE *fp;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;

	if( file == NULL ) {
	    fp=fopen("new.grd","w");
	} else {
	    fp=fopen(file,"w");
	}

	/* comments */

        fprintf(fp,"\n");
        fprintf(fp,"0 (mesh) automatic generated grid\n");
        fprintf(fp,"\n");

	/* nodes */

	ResetHashTable(HNN);
	while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( !all && IsNtype(pn, N_EXTERNAL ) ) continue;
                fprintf(fp,"1 %d %d %f %f"
                        ,pn->number
                        ,(int) pn->type
                        ,pn->coord.x
                        ,pn->coord.y
                        );
                if( pn->depth != NULLDEPTH )
                        fprintf(fp," %f\n",pn->depth);
                else
                        fprintf(fp,"\n");

	}

        fprintf(fp,"\n");

	/* elements */

        for(i=1;i<=NTotElems;i++) {
          if( (pe=RetrieveByElemNumber(HEL,i)) != NULL ) {
		if( !all && IsEtype(pe, E_EXTERNAL ) ) continue;
                fprintf(fp,"2 %d %d %d"
                        ,pe->number
                        ,(int) pe->type
                        ,pe->vertex
                        );

                for(j=0;j<pe->vertex;j++) {
                        if( j%10 == 0 && pe->vertex > 3 )
                                fprintf(fp,"\n");
                        fprintf(fp," %d",pe->index[j]);
                }
        	fprintf(fp,"\n");
	  }
	}

        fprintf(fp,"\n");

	/* lines */


        for(i=1;i<=NTotLines;i++) {
          if( (pl=RetrieveByLineNumber(HLI,i)) != NULL ) {
	    if( IsLtype(pl,L_EXTERNAL_REF) || IsLtype(pl,L_INTERNAL_REF) 
			|| IsLtype(pl,L_FAULT_REF)
			|| EqualsLtype(pl,L_NONE) ) {
		if( IsLtype(pl,L_EXTERNAL_REF) ) {
		  linetype = L_EXTERNAL;
		} else if( IsLtype(pl,L_INTERNAL_REF) ) {
		  linetype = L_INTERNAL;
		} else if( IsLtype(pl,L_FAULT_REF) ) {
		  linetype = L_FAULT;
		} else {
		  linetype = 0;
		}
                fprintf(fp,"3 %d %d %d"
                        ,pl->number
                        ,linetype
                        ,pl->vertex
                        );

                for(j=0;j<pl->vertex;j++) {
                        if( j%10 == 0 )
                                fprintf(fp,"\n");
                        fprintf(fp," %d",pl->index[j]);
                }
                fprintf(fp,"\n");
	    }
          }
        }

        fprintf(fp,"\n");

	/* list */

	if( list ) {
		count = list->count;
        	fprintf(fp,"3 1 1 %d",count+1);
        	for(j=0;j<=count;j++) {
                	if( j%10 == 0 ) fprintf(fp,"\n");
                	fprintf(fp," %d",list->index[j%count]);
        	}
        	fprintf(fp,"\n");
	}

	fclose(fp);
}

