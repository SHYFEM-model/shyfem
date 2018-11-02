
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
 * meshvs.c - version description of mesh                               *
 *									*
 * Revision History:							*
 * 02-Nov-2018: version 1.80                                            *
 * 10-Mar-2010: version 1.77                                            *
 * 21-Apr-2009: version 1.75                                            *
 * 21-Jan-2009: version 1.74                                            *
 * 14-Jan-2009: version 1.73a                                           *
 * 14-Jan-2009: version 1.73                                            *
 * 18-Oct-2008: version 1.72                                            *
 * 23-Oct-2005: version 1.71                                            *
 * 19-Aug-2003: version 1.70                                            *
 * 14-Aug-2002: version 1.66                                            *
 * 03-Jul-2000: version 1.65                                            *
 * 03-Feb-2000: version 1.64                                            *
 * 12-May-99: version 1.63                                              *
 * 05-May-99: version 1.62                                              *
 * 25-Nov-97: version 1.61                                              *
 * 12-Nov-97: version 1.60                                              *
 * 17-Oct-97: version 1.51                                              *
 * 16-Oct-97: version 1.50                                              *
 * 15-Oct-97: version 1.41                                              *
 * 08-Oct-97: version 1.40                                              *
 * 17-Aug-95: version 1.31                                              *
 * 14-Aug-95: version 1.30                                              *
 * 11-Aug-95: version 1.20                                              *
 * 11-Aug-95: version 1.11                                              *
 * 10-Aug-95: version 1.10                                              *
 * 07-Aug-95: version 1.01                                              *
 * 02-Aug-95: version 1.00                                              *
 * 01-Aug-95: version 0.80                                              *
 * 31-Jul-95: version 0.70                                              *
 * 28-Jul-95: version 0.50                                              *
 * 27-Jul-95: version 0.40                                              *
 * 26-Jul-95: version 0.20                                              *
 * 25-Jul-95: version 0.10                                              *
 *									*
\************************************************************************/

#include <stdio.h>

char* SCopy  = "Copyright (C) 1995-2018  Georg Umgiesser                ";
char* SMesh  = "MESH - Automatic Grid Generation Routine - Version 1.80 ";

void Logos( void )

{
        printf("\n");
	printf("%s\n",SMesh);
	printf("%s\n",SCopy);
	printf("\n");
}

/**********************************************************************\


=========================================================================
========================== version log for mesh =========================
=========================================================================

version 1.80						      02 Nov 2018

        copyright updated
        use getopt library from libc

=========================================================================

version 1.77						      10 Mar 2010

	new algorithm for computing area of polygon (stable for 64 bit)

=========================================================================

version 1.75						      21 Apr 2009

	DIFFS file

=========================================================================

version 1.74						      21 Jan 2009

	new makedepend

=========================================================================

version 1.73a						      14 Jan 2009

	minor point

=========================================================================

version 1.73						      14 Jan 2009

	Integrated into SHYFEM tree

=========================================================================

version 1.71						      23 Oct 2005

	New tagged version, minor changes

=========================================================================

version 1.70						      19 Aug 2003

	New tagged version

=========================================================================

version 1.66						      14 Aug 2002

	InterpolRes

 * mesh.c:
 * 14-Aug-2002: InterpolRes uses cubic to compute resolution            *

=========================================================================

version 1.65						      03 Jul 2000

	ControlCircumCircle

 * meshge.c:
 * 03-Jul-2000: in ControlCircumCircle() only warning                   *

=========================================================================

version 1.64						      03 Feb 2000

	changes in background grid

 * mesh.c:
 * 03-Feb-2000: different algorithm to compute resolution in bckgrid    *

=========================================================================

version 1.63							12 May 99

	bug fix in meshut

 * meshut.c:
 * 12-May-1999: bug fix in InsertInIndex() -> insert only once          *
 * 12-May-1999: new routine PrintLine() to print one line               *

=========================================================================

version 1.62							05 May 99

	some more checks

 * mesh.c:
 * 05-May-1999: check for zero area in SetFEMParam                      *

=========================================================================

version 1.61							25 Nov 97

	new names for .grd files

- mesh.c:       new names for .grd files -> M_*.grd

=========================================================================

version 1.60							12 Nov 97

	refining ok even for domain with many boundaries
	background grid need not cover whole area

- mesh.c:       new routine MarkOuterNodes() to get rid of unused nodes
                MarkOuterElements() and MarkOuterNodes() are called
                  after internal refining -> internal nodes can be
                  inserted even if outer element has to be changed
                  -> whole routine gets more flexible
                  RecoverBoundaryNodes() is called after internal refining
                in InterpolRes() the standard resolution (OpResolution)
                  is returned if no background triangle is found
- meshbd.c:     in RecoverBoundaryNodes() interpolate also depth
                  -> MakeMidPoint() also returns interpolated depth
- meshck.c:     check for not unique coordinates in CheckInput()

=========================================================================

version 1.51							17 Oct 97

	check for counter-clockwise turning of line

- meshck.c:     in CheckInput() check for couter-clockwise line turning
- meshge.c:     new routine TurnClosedLine()
- meshge.h:     new routine TurnClosedLine()

=========================================================================

version 1.50							16 Oct 97

	recovering of boundary lines implemeted

- mesh.c:       call CopyBoundaryLine() in main instead RefineBoundary()
                  new call to RecoverBoundaryNodes() introduced
- meshbd.c:     in RefineBoundary() use always new line and delete old
                  new routines RecoverBoundaryNodes(), FindElemToSide()
                  MakeMidPoint()
- meshbd.h:     new routines RecoverBoundaryNodes(), FindElemToSide()
- meshge.c:     new routine FindElement()
- meshge.h:     new routine FindElement()
- meshut.c:     new routines PrintLineList(), MakeNewLine()
                  InsertInIndex(), InsertNodeInLine()
- meshut.h:     new routines PrintLineList(), MakeNewLine()
                  InsertInIndex(), InsertNodeInLine()

=========================================================================

version 1.41							15 Oct 97

	=== bug fix for 1.40 ===

	inside check in MarkOuterElements with InClosedLine()
		-> ok also for clockwise lines
	in MakeNode set depth to NULLDEPTH (floating exception)

- mesh.c:       MarkOuterElements: use InClosedLine() for inside check
                  use MakeCoordsFromLine() to make coordinate list
- meshfi.c:     in WriteFile: write also lines with type L_NONE
- meshge.c:     new routine InClosedLine()
- meshge.h:     new routine InClosedLine()
- meshin.c:     in GivenNodes: insert given nodes only if internal
- meshut.c:     in MakeNode(): set depth to NULLDEPTH (bug)
                  new routine MakeCoordsFromLine()
- meshut.h:     new routine MakeCoordsFromLine()

=========================================================================

version 1.40							08 Oct 97

	=== This is an intermediate release - still buggy ===

	New mesh type introduced
	Routines to insert given internal nodes
	Minor enhancements

- mesh.c:	insert given points in domain delimited by line
		  introduced pointer to background grid
		  use new mesh type
- mesh.h:	standard structures uniformed
		  new extra structures Extra_N/E/L_type
		  new node type N_GIVEN
- meshbd.c:	uses new mesh type
- meshck.c:	uses new mesh type
		  check for closed line only for L_EXTERNAL/L_INTERNAL
- meshcv.c:	uses new mesh type
- meshfi.c:	routine ReadBnd deleted
		  if depth is given -> write it
		  uses mesh type to decide which elements to write
- meshge.c:	error message for identical nodes in MakeCircumCircle()
- meshin.c:	uses new mesh type
		  commented sections deleted
		  new routines InsertNodes(), GivenNodes()
- meshin.h:	new routine prototypes InsertNodes(), GivenNodes()
- meshty.c:	new file
- meshty.h:	new file
- meshut.c:	routines slighlty restructured
		  extra structure introduced -> holds use, type, ...
		  new items are created with type 0, but adeguate mesh type
- nlist.c:	handle 0 list correctly

=========================================================================

version 1.31							17 Aug 95

- options.c created and deleted from meshop.c

=========================================================================

version 1.30							14 Aug 95

- some changes in heap table
- background grid implemented

=========================================================================

version 1.20							11 Aug 95

- read grd files as input, bnd files are obsolete

=========================================================================

version 1.11							11 Aug 95

- some files (stack,queue,list,nlist,heap) restructured
	mesh files adapted

=========================================================================

version 1.10							10 Aug 95

- heap structure has changed -> insert (void *)
- eliminate external elements completly
- new help screen

=========================================================================

version 1.01							 7 Aug 95

- smoothing after refining, because circumcircle data
	is not valid anymore
- some switches to control the generation of the mesh,
	the resolution and the insertion of internal
	nodes

=========================================================================

version 1.00							 2 Aug 95

file clean up - first running version

=========================================================================

version 0.8							 1 Aug 95

restructured some files
revised insertion of internal points
insertion of internal points in deformed triangles

=========================================================================

version 0.7							31 Jul 95

finished insertion of internal points

=========================================================================

version 0.5							28 Jul 95	

finished refinement of boundary points,
	remeshing of all boundary points,
	use function for definition of resolution
four bugs found with turbo-c :
	- MakeMaxiTriang: y < xmin (why did it ever work ?)
	- in meshcv.c -> abs() substituted with ABS()
	- in meshfi.c -> <strings.h> substituted with <string.h>
	- InElement: routine falsly accepted points that were
		outside of element -> loop one more time

=========================================================================

version 0.4							27 Jul 95	

finished meshing of convex hull

=========================================================================

version 0.2							26 Jul 95	

finished data structures, convex hull

=========================================================================

version 0.1							25 Jul 95	

started project -> some files reused from grid

=========================================================================

\**********************************************************************/
