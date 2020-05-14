
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999-2000,2009-2011,2015-2016,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c subroutines for line manipolations
c
c contents :
c
c function klineck(n,kline)			checks compatibility of kline
c subroutine linedbl(n,kline,bstop)		tests for uniqueness
c subroutine lineadj(n,kline,bstop)		tests for adjacency
c function iskadj(k1,k2)			tests for adjacency of nodes
c
c function nextline(inode,ndim,nnode,ifirst,ilast)	gets next line
c function extrline(inode,ndim,nline,ifirst,ilast)	extracts line from array
c
c subroutine revline(n,kline)			reverses line
c
c revision log :
c
c 19.11.1999	ggu	isolated from subflxa
c 02.12.1999	ggu	new routine extrline, nextline substitutes extrli
c 20.01.2000	ggu	new routine revline
c 28.04.2009	ggu	links re-structured
c 23.03.2010	ggu	changed v6.1.1
c 18.10.2011	ggu	changed VERS_6_1_33
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 06.06.2016	ggu	module added, new routines for bnd and grd lines
c 16.02.2019	ggu	changed VERS_7_5_60
c
c note :
c
c line nodes are stored consecutive in inode (dimension ndim)
c a value of zero devides lines from each other
c any number of zeros between lines (and leading and trailing) are allowed
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
	module lines
!==================================================================

        implicit none

        type, private :: entry

          integer :: nlines
          integer :: nnodes
          integer :: nfill
          integer, allocatable :: inum(:)
          integer, allocatable :: itype(:)
          integer, allocatable :: start(:)
          integer, allocatable :: end(:)
          real, allocatable :: depth(:)
          integer, allocatable :: nodes(:)
          real, allocatable :: xv(:)
          real, allocatable :: yv(:)

        end type entry

        integer, save, private :: idlast = 0
        integer, save, private :: ndim = 0
        type(entry), save, private, allocatable :: pentry(:)

!==================================================================
	contains
!==================================================================

        subroutine lines_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
          return
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine lines_init_alloc

!******************************************************************

        subroutine lines_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call lines_init_alloc
        end if
        id = idlast

        call lines_init_id(id)

        end subroutine lines_init_new_id

!******************************************************************

        subroutine lines_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop lines_init_id: ndim'
        end if

        pentry(id)%nlines = 0

        end subroutine lines_init_id

!******************************************************************

        subroutine lines_init_new_lines(id,nli,nnodes)

        integer id
        integer nli
        integer nnodes

        call lines_init_new_id(id)

        pentry(id)%nlines = nli
        pentry(id)%nnodes = nnodes
        pentry(id)%nfill = 0

	allocate(pentry(id)%inum(nli))
	allocate(pentry(id)%itype(nli))
	allocate(pentry(id)%start(nli))
	allocate(pentry(id)%end(nli))
	allocate(pentry(id)%depth(nli))

	allocate(pentry(id)%nodes(nnodes))
	allocate(pentry(id)%xv(nnodes))
	allocate(pentry(id)%yv(nnodes))

        end subroutine lines_init_new_lines

!******************************************************************

	subroutine lines_set_line(id,il,in,it,n,d,nodes,x,y)

        integer id
        integer il
        integer in,it
	real d
	integer n
	integer nodes(n)
	real x(n)
	real y(n)

	integer ibase
	integer nli
	integer nnodes

	nli = pentry(id)%nlines
	nnodes = pentry(id)%nnodes
	ibase = pentry(id)%nfill

	if( il < 1 .or. il > nli ) then
	  write(6,*) 'il,nli: ',il,nli
	  stop 'error stop lines_set_line: il out of range'
	end if

	if( ibase + n > nnodes ) then
	  write(6,*) 'ibase,n,nnodes: ',ibase,n,nnodes
	  stop 'error stop lines_set_line: too many nodes'
	end if

	pentry(id)%nfill = ibase + n

	pentry(id)%inum(il) = in
	pentry(id)%itype(il) = it
	pentry(id)%start(il) = ibase + 1
	pentry(id)%end(il) = ibase + n
	pentry(id)%depth(il) = d

	pentry(id)%nodes(ibase+1:ibase+n) = nodes
	pentry(id)%xv(ibase+1:ibase+n) = x
	pentry(id)%yv(ibase+1:ibase+n) = y

	end subroutine lines_set_line

!******************************************************************

	subroutine lines_get_line(id,il,in,it,d,n,nodes,x,y)

        integer id
        integer il
        integer in,it
	real d
	integer n
	integer nodes(n)
	real x(n)
	real y(n)

	integer is,ie

	in = pentry(id)%inum(il)
	it = pentry(id)%itype(il)
	is = pentry(id)%start(il)
	ie = pentry(id)%end(il)
	d = pentry(id)%depth(il)

	nodes = pentry(id)%nodes(is:ie)
	x = pentry(id)%xv(is:ie)
	y = pentry(id)%yv(is:ie)

	end subroutine lines_get_line

!******************************************************************

	subroutine lines_make_node_list(id,nnodes,nlist)

	integer id,nnodes,nlist

	logical binline
	integer nli,i,n

	!-----------------------------------
	! count lines in list (one or more 0 are used to seperate lines)
	!-----------------------------------

	nli = 0
	binline = .false.
	do i=1,nnodes
	  n = nlist(i)
	  if( n > 0 ) then
	    if( .not. binline ) then
	      nli = nli + 1
	      binline = .true.
	    end if
	  else
	    if( binline ) then
	      binline = .false.
	    end if
	  end if
	end do

	!-----------------------------------
	! initialize new list
	!-----------------------------------

        call lines_init_new_lines(id,nli,nnodes)

	pentry(id)%itype(:) = 0
	pentry(id)%depth(:) = 0.
	pentry(id)%xv(:) = 0.
	pentry(id)%yv(:) = 0.

	!-----------------------------------
	! insert in list
	!-----------------------------------

	nli = 0
	binline = .false.
	do i=1,nnodes
	  n = nlist(i)
	  if( n > 0 ) then
	    if( .not. binline ) then
	      nli = nli + 1
	      binline = .true.
	      pentry(id)%inum(nli) = nli
	      pentry(id)%start(nli) = i
	    end if
	  else
	    if( binline ) then
	      binline = .false.
	      pentry(id)%end(nli) = i-1
	    end if
	  end if
	end do
	if( binline ) pentry(id)%end(nli) = nnodes

	!-----------------------------------
	! end of routine
	!-----------------------------------

	end subroutine lines_make_node_list

!==================================================================
	end module lines
!==================================================================

c******************************************************************
c******************************************************************
c******************************************************************

	function klineck(n,kline)

c checks compatibility of kline -> returns number of lines or -1 (error)

	implicit none

	integer klineck
	integer n		!size of kline
	integer kline(n)

	integer nnode,ifirst,ilast,ntotal
	integer nlines
	logical bstop
	logical nextline

	bstop = .false.
	nnode = 0
	nlines = 0

	do while( nextline(kline,n,nnode,ifirst,ilast) )
	  nlines = nlines + 1
	  ntotal = ilast - ifirst + 1
	  call linedbl(ntotal,kline(ifirst),bstop) !tests if nodes unique
	  call lineadj(ntotal,kline(ifirst),bstop) !tests if nodes are adjacent
	end do

	if( bstop ) nlines = -1

	klineck = nlines

	end

c******************************************************************

	subroutine linedbl(n,kline,bstop)

c tests for uniqueness

	implicit none

	integer n
	integer kline(n)
	logical bstop

	integer i,k,j
	integer istart

	integer ipext

	istart = 1
	if( kline(1) .eq. kline(n) ) istart = 2		!allow for closed line

	do i=istart,n
	   k = kline(i)
	   do j=i+1,n
	      if( kline(j) .eq. k ) then
		bstop = .true.
		write(6,*) 'node is not unique ',ipext(k)
	      end if
	   end do
	end do

	end

c******************************************************************

	subroutine lineadj(n,kline,bstop)

c tests for adjacency

	implicit none

	integer n
	integer kline(n)
	logical bstop

	integer i,k1,k2

	integer ipext
	logical iskadj

	do i=2,n
	   k1 = kline(i-1)
	   k2 = kline(i)

	   if( .not. iskadj(k1,k2) ) then
		bstop = .true.
		write(6,*) 'nodes are not adjacent ',ipext(k1),ipext(k2)
	   end if
	end do

	end

c******************************************************************

	function iskadj(k1,k2)

c tests for adjacency between nodes

	use mod_geom
	use basin

	implicit none

	logical iskadj
	integer k1,k2
	integer elems(maxlnk)

	integer n,i,ie,ii

	iskadj = .true.

	call get_elems_around(k1,maxlnk,n,elems)

	do i = 1,n
	  ie = elems(i)
	  do ii=1,3
	    if( nen3v(ii,ie) .eq. k2 ) return
	  end do
	end do

	iskadj = .false.

	end

c**********************************************************************

	function nextline(inode,ndim,nnode,ifirst,ilast)

c gets next line from array
c
c on entry:
c           nnode points to first node to be analysed
c on return inode(ifirst) is first node of line and
c           inode(ilast)  is last  node of line
c           nnode points to first node to be analysed on next entry
c
c on first call nnode must be 0, afterwards nnode gets set by algorithm
c node list must contain zeros to devide lines from each other
c any number of consecutive zeros can be used
c leading and trailing zeros are allowed

	implicit none

	logical nextline	!true if new line found
	integer inode(ndim)	!list of nodes that define line    (in)
	integer ndim		!total number of nodes given in node() (in)
	integer nnode		!start/end of line in inode()       (in/out)
	integer ifirst		!start of line in inode()           (out)
	integer ilast		!end of line in inode()             (out)

	integer iw,i,node

	iw = -1				! -1 befor  0 in   +1 after line
	if( nnode .le. 0 ) nnode = 1
	i = nnode
	ifirst = 0
	ilast = 0

	do while( i .le. ndim .and. iw .ne. 1 )
	    node = inode(i)
	    if( node .lt. 0 ) then
	        stop 'error stop nextline: negative nodes in line'
	    else if( node .eq. 0 ) then
	        if( iw .eq. 0 ) then	!end of line
		    ilast = i - 1
		    iw = +1
	        end if
	    else !if( node .gt. 0 ) then
	        if( iw .eq. -1 ) then	!start of line
		    ifirst = i
		    iw = 0
	        end if
	    end if
	    i = i + 1
	end do

c end last line if not already done

	if( iw .eq. 0 ) ilast = ndim		!last line

	if( iw .eq. -1 ) then		!no more lines, only trailing blanks
		nnode = 0
		nextline = .false.
	else !if( iw .eq. 1 ) then	!line ended regularily
		nnode = i
		nextline = .true.
	end if

	end

c**********************************************************************

	function extrline(inode,ndim,nline,ifirst,ilast)

c extracts given line number from array
c
c on return inode(ifirst) is first node of line and
c           inode(ilast)  is last  node of line

	implicit none

	logical extrline	!true if new line found
	integer inode(ndim)	!list of nodes that define line    	(in)
	integer ndim		!total number of nodes given in node() 	(in)
	integer nline		!number of line to extract		(in)
	integer ifirst		!start of line in inode()           	(out)
	integer ilast		!end of line in inode()             	(out)

	integer nnode,nlines
	logical bfound
	logical nextline

	nnode = 0
	nlines = 0
	bfound = .false.

c loop until no more line or line found

	do while( nextline(inode,ndim,nnode,ifirst,ilast) .and. 
     +				.not. bfound )
	  nlines = nlines + 1
	  if( nline .eq. nlines ) bfound = .true.
	end do

	extrline = bfound

	end

c**********************************************************************

	subroutine revline(n,kline)

c reverses line

	implicit none

	integer n		!total number of nodes
	integer kline(n)	!nodes in line

	integer i,ihigh,kmem

	do i=1,n/2
	  ihigh = n + 1 - i
	  kmem = kline(ihigh)
	  kline(ihigh) = kline(i)
	  kline(i) = kmem
	end do

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine bnd_internal_file(file,binsert,nli,nnodes,xv,yv,ind)

	implicit none

	character*(*) file
	logical binsert
	integer nli		!number of lines in file, -1 for error
	integer nnodes		!total number of nodes (x,y)
	real xv(nnodes)
	real yv(nnodes)
	integer ind(nnodes)

	integer i,iunit,ios
	real x,y
	integer ifileo

	iunit = ifileo(0,file,'form','old')
	if( iunit <= 0 ) goto 99

	nli = 0
	nnodes = 0
	do
	  read(iunit,*,iostat=ios) x,y,i
	  if( ios < 0 ) return
	  if( ios > 0 ) goto 99
	  if( i /= 0 .and. i /= 1 ) goto 99
	  if( i == 1 ) nli = nli + 1
	  nnodes = nnodes + 1
	  if( binsert ) then
	    xv(nnodes) = x
	    yv(nnodes) = y
	    if( i == 1 ) ind(nli) = nnodes
	  end if
	end do

	return
   99	continue
	nli = -1
	end

c**********************************************************************

	subroutine bnd_read_file(file,nli,nnodes,xv,yv,ind)

	implicit none

	character*(*) file
	integer nli		!number of lines in file, -1 for error
	integer nnodes		!total number of nodes (x,y)
	real xv(nnodes)
	real yv(nnodes)
	integer ind(nnodes)

	call bnd_internal_file(file,.true.,nli,nnodes,xv,yv,ind)

	end

c**********************************************************************

	subroutine bnd_test_file(file,nli,nnodes)

	implicit none

	character*(*) file
	integer nli		!number of lines in file, -1 for error
	integer nnodes		!total number of nodes (x,y)

	real xv(1)
	real yv(1)
	integer ind(1)

	call bnd_internal_file(file,.false.,nli,nnodes,xv,yv,ind)

	end

c**********************************************************************

	function bnd_is_bnd_file(file)

	implicit none

	logical bnd_is_bnd_file
	character*(*) file

	integer nli		!number of lines in file, -1 for error
	integer nnodes		!total number of nodes (x,y)

	call bnd_test_file(file,nli,nnodes)

	bnd_is_bnd_file = nli > 0

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine bnd_insert_lines(file,id)

	use lines

	implicit none

	character*(*) file
	integer id

	integer nli,nnodes
	integer is,n,il,it
	real d
	real, allocatable :: xv(:)
	real, allocatable :: yv(:)
	integer, allocatable :: nodes(:)
	integer, allocatable :: ind(:)

	id = -1
	call bnd_test_file(file,nli,nnodes)
	if( nli <= 0 ) return

	allocate(xv(nnodes),yv(nnodes),nodes(nnodes),ind(nnodes+1))
        call lines_init_new_lines(id,nli,nnodes)

	call bnd_read_file(file,nli,nnodes,xv,yv,ind)
        call lines_init_new_lines(id,nli,nnodes)
	ind(nli+1) = nnodes + 1
	nodes = 0
	it = 0
	d = 0.

	do il=1,nli
	  is = ind(il)
	  n = ind(il+1) - is
	  call lines_set_line(id,il,il,it,n,d
     +				,nodes,xv(is:),yv(is:))
	end do

	deallocate(xv,yv,nodes,ind)

	end

c**********************************************************************

	subroutine grd_insert_lines(file,id)

	use lines

	implicit none

	character*(*) file
	integer id

	integer nli,nnodes
	integer il,in,it,n
	real d
	real, allocatable :: xv(:)
	real, allocatable :: yv(:)
	integer, allocatable :: nodes(:)

	call grd_read(file)

	call grd_get_total_lines(nli)
	call grd_get_total_nodes(nnodes)	!maximum

	allocate(xv(nnodes),yv(nnodes),nodes(nnodes))

	do il=1,nli
	  call grd_get_line_params(il,in,it,n,d)
	  call grd_get_line_array(il,n,nodes,xv,yv)
	  call lines_set_line(id,il,in,it,n,d,nodes,xv,yv)
	end do

	deallocate(xv,yv,nodes)

	end

c**********************************************************************

	subroutine read_lines(file,id)

	implicit none

	character*(*) file
	integer id

	logical bnd_is_bnd_file,is_grd_file

	id = 0
	if( file == ' ' ) return

	if( bnd_is_bnd_file(file) ) then
	  call bnd_insert_lines(file,id)
	else if( is_grd_file(file) ) then
	  call grd_insert_lines(file,id)
	end if
	
	end

c**********************************************************************

