c
c $Id: sublin.f,v 1.4 2009-05-21 09:24:00 georg Exp $
c
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
c 19.11.1999    ggu	isolated from subflxa
c 02.12.1999    ggu	new routine extrline, nextline substitutes extrli
c 20.01.2000    ggu	new routine revline
c 28.04.2009    ggu     links re-structured
c
c note :
c
c line nodes are stored consecutive in inode (dimension ndim)
c a value of zero devides lines from each other
c any number of zeros between lines (and leading and trailing) are allowed
c
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

	include 'param.h'

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

