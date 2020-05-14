
/************************************************************************\
 *
 *    Copyright (C) 1992,1994-1995,2012  Georg Umgiesser
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
 * meshgd.c - read/write grd files
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 13.04.1994	ggu	use new hash routines
 * 14.04.1994	ggu	use GetActFileType to determine file type
 * 06.05.1994	ggu	new file ff created for old read/write
 * ...		ggu	filetype 0 is now the "official" filetype
 * 07.05.1994	ggu	new routines for opening and reading file
 * 08.10.1994	ggu	reading/writing comments (QueueTable routines)
 * 21.10.1994	ggu	Changed introduced -> write only if file is changed
 * 10.02.1995	ggu	closefile calls fclose only if file opened
 * 11.08.1995	ggu	split from gridfi to make own file
 * 16.02.2012	ggu	no cast of (char) for type of element and line
 *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"

#include "fund.h"
#include "gustd.h"
#include "hash.h"
#include "queue.h"
#include "args.h"

#include "mesh.h"
#include "meshhs.h"
#include "meshut.h"
#include "meshgd.h"


/**************************************************************************/

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

static int gettotlines( void ) { return FpLine; }
static void closefile( void ) { if(FpFile) fclose(FpFile); }

static char *stripext( char *s , char *ext ) /* strips ext from file name */

{
        int ls;
        char *t;

        ls = strlen(s);
        ls -= 4;

        if( ls > 0 ) {
                t = &s[ls];
                if( !strcmp(t,ext) ) {
                        t[0] = '\0';
                }
        }

        return s;
}

static char *makefilename( char *name , char *ext )
{
        name = stripext(name,ext);
        name = strcat(name,ext);
        return name;
}

/**************************************************************************/


void ReadStandard( char *fname , Hashtable_type HN , Hashtable_type HE
			      , Hashtable_type HL , QueueTable C )

{
	int comms=0,nodes=0,elems=0,lines=0;
	int nodemax=0,elemmax=0,linemax=0;
	int err=0;
	int error=FALSE;
	int narg,what,n;
	char *s,*t;

printf("...opening file %s\n",fname);
	fname = makefilename(fname,".grd");
	if( openfile(fname,"r") )
		printf("Reading file %s\n",fname);
	else
		Error2("ReadStandard : Cannot open file ",fname);

	while( (s=getnewline()) != NULL ) {
		t=firstchar(s);
		if( *t == '0' ) {
			t=savestring(s,-1);
			what=0;
		} else {
			narg = nargs(s);
			if(narg == 0) continue;
			t=readargs();
			what=atoi(t);
		}

		switch(what) {
		case 0 :			/* comment */
			EnQueue(C,(void *)t);
			comms++;
			break;
		case 1 :			/* node */
			n = ReadNode(HN);
			if( n ) nodes++ ; else error=TRUE ;
			if( n > nodemax ) nodemax=n;
			break;
		case 2 :			/* element */
			n = ReadElem(HE);
			if( n ) elems++ ; else error=TRUE ;
			if( n > elemmax ) elemmax=n;
			break;
		case 3 :			/* line */
			n = ReadLine(HL);
			if( n ) lines++ ; else error=TRUE ;
			if( n > linemax ) linemax=n;
			break;
		default:
			err++;
			printf("Line %d : ",gettotlines());
			printf("Shape %s not recognized\n",t);
			break;
		}
		if( error ) {
			printf("Read error in line %d :\n",gettotlines());
			error=FALSE;
			err++;
		}
	}

	NTotNodes += nodemax;
	NTotElems += elemmax;
	NTotLines += linemax;

	printf("%d lines read\n",gettotlines());
	printf("Following shapes read :\n");
	if(comms) printf("Comments : %d ",comms);
	if(nodes) printf("Nodes : %d ",nodes);
	if(elems) printf("Elements : %d ",elems);
	if(lines) printf("Lines : %d ",lines);
	printf("\n");
	if( err ) {
		Error("Errors detected in input file");
	}

	closefile();
}

void WriteStandard( char *fname , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , QueueTable C )

{
	FILE *fp;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;
	int nodes=0,elems=0,lines=0,comments=0;
	int i,j;
	char *s;

        fp=fopen(fname,"w");
        if( fp )
                printf("Writing file %s\n",fname);
        else
                Error2("WriteStandard : Cannot open file ",fname);

	while( (s=(char *)DeQueue(C)) != NULL ) {
		fprintf(fp,"%s\n",s);
		comments++;
	}
        printf("Comments written : %d\n",comments);

	fprintf(fp,"\n");

        for(i=1;i<=NTotNodes;i++) {
          if( (pn=RetrieveByNodeNumber(HN,i)) != NULL ) {
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
                nodes++;
          }
        }
        printf("Nodes written : %d\n",nodes);

	fprintf(fp,"\n");

        for(i=1;i<=NTotElems;i++) {
          if( (pe=RetrieveByElemNumber(HE,i)) != NULL ) {
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

                if( pe->depth != NULLDEPTH )
                        fprintf(fp," %f\n",pe->depth);
                else
                        fprintf(fp,"\n");
                elems++;
          }
        }
        printf("Elements written : %d\n",elems);

	fprintf(fp,"\n");

        for(i=1;i<=NTotLines;i++) {
          if( (pl=RetrieveByLineNumber(HL,i)) != NULL ) {
                fprintf(fp,"3 %d %d %d"
                        ,pl->number
                        ,(int) pl->type
                        ,pl->vertex
                        );

                for(j=0;j<pl->vertex;j++) {
                        if( j%10 == 0 )
                                fprintf(fp,"\n");
                        fprintf(fp," %d",pl->index[j]);
                }

                fprintf(fp,"\n");
                lines++;
          }
        }
        printf("Lines written : %d\n",lines);

	fclose(fp);
}

int ReadNode( Hashtable_type H )

{
	char *t;
	int number,ntype;
	Point c;
	float depth;
	Node_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	c.x = atof(t);

	t=readargs();
	if( !t ) return 0;
	c.y = atof(t);

	t=readargs();
	if( !t )
		depth = NULLDEPTH;
	else
		depth = atof(t);

	p=MakeNode(number+NTotNodes,ntype,&c);
	p->depth = depth;
	InsertByNodeNumber(H,p);

	return number;
}

int ReadElem( Hashtable_type H )

{
	char *t,*s;
	int number,ntype;
	int i,vertex;
	int *index;
	float depth;
	Elem_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	vertex = atoi(t);

	index = MakeIndex(vertex);

	i=0;
	while( i<vertex ) {
		t=readargs();
		if( !t ) {
			s=getnewline();
			if( !s ) return 0;
			initargs(s);
		} else {
			index[i++] = atoi(t) + NTotNodes;
		}
	}

	t=readargs();
	if( !t )
		depth = NULLDEPTH;
	else
		depth = atof(t);

	p = MakeElemWithIndex(number+NTotElems,ntype,vertex,index);
	p->depth = depth;
	InsertByElemNumber(H,p);

	return number;
}

int ReadLine( Hashtable_type H )

{
	char *t,*s;
	int number,ntype;
	int i,vertex;
	int *index;
	Line_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	vertex = atoi(t);

	index = MakeIndex(vertex);

	i=0;
	while( i<vertex ) {
		t=readargs();
		if( !t ) {
			s=getnewline();
			if( !s ) return 0;
			initargs(s);
		} else {
			index[i++] = atoi(t) + NTotNodes;
		}
	}

	p = MakeLineWithIndex(number+NTotLines,ntype,vertex,index);
	InsertByLineNumber(H,p);

	return number;
}

