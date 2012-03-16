c
c $Id: subgrd.f,v 1.9 2009-05-21 09:24:00 georg Exp $
c
c rdgrd routines - read GRD files
c
c contents :
c
c       subroutine rdgrd(...)
c				reads grd file
c       subroutine rdcom(...)
c				reads comment from .grd file
c       subroutine rdnode(...)
c				reads node from .grd file
c       subroutine rdelem(...)
c				reads element from .grd file
c       subroutine rdline(...)
c				reads line from .grd file
c	subroutine read_node_list(...)
c				reads node list
c	subroutine rdunknown(iwhat,bstop)
c				handles unknown type
c
c       function ifstch(line)
c				finds first char of line that is not blank
c       subroutine fempar(line)
c				read parameters for fem model
c
c	subroutine ex2in(nkn,ne,nl,ipnv,ipaux,nen3v,inodlv,bstop)
c				changing extern with intern node numbers
c	subroutine chex2in(nkn,n,nodes,ipnv,index,bstop)
c				changing extern with intern node numbers (list)
c
c	internal routines:
c
c	subroutine grd_init(file)
c	function grd_next_line()
c	subroutine grd_nvals(nvals)
c	subroutine grd_vals(nvals,vals)
c	subroutine grd_ival(ival,val)
c	subroutine grd_line_info(iline_grd,line_grd)
c	subroutine grd_write_line_info
c
c	writing routines
c
c	subroutine grd_write_grid(
c				writes grd file
c	subroutine grd_write_node_list(nout,n,nodes,ipnv,depth)
c				writes out node list
c
c revision log :
c
c 20.03.1998    ggu     declare f(6) in subroutines to avoid compiler warnings
c 20.05.1998    ggu     open file with ifileo()
c 22.05.1998    ggu     bug fix (ifileo not declared)
c 27.03.2001    ggu     assign depth -999. to depth not set
c 09.10.2001    ggu     read from stdin
c 09.10.2001    ggu     read node type (ianv)
c 18.10.2005    ggu     error messages slightly changed
c 22.04.2009    ggu     changes from spline integrated (read lines)
c 24.04.2009    ggu     newly restructured
c 09.03.2012    ggu     handle dimension error more gracefully
c
c**********************************************************

	subroutine rdgrd(
     +			 file
     +			,bstop
     +			,nco,nkn,nel,nli
     +			,nkndim,neldim,nlidim,nlndim
     +			,ipnv,ipev,iplv
     +			,ianv,iaev,ialv
     +			,hnv,hev,hlv
     +			,xgv,ygv
     +                  ,nen3v
     +                  ,ipntlv,inodlv
     +			)

c reads grd file
c
c after read nen3v and inodlv contain still external node numbers
c use ex2in to convert external to internal node numbers
c
c works only with triangles as elements

	implicit none

	character*(*) file	!file name

	logical bstop		!true on return if error

	integer nco		!total number of comments read
	integer nkn		!total number of nodes read
	integer nel		!total number of elements read
	integer nli		!total number of lines read

	integer nkndim		!dimension for number of nodes
	integer neldim		!dimension for number of elements
	integer nlidim		!dimension for number of lines
	integer nlndim		!dimension for node numbers of lines

	integer ipnv(nkndim)	!external node number
	integer ipev(neldim)	!external element number
	integer iplv(nlidim)	!external line number

	integer ianv(nkndim) 	!node type
	integer iaev(neldim)	!element type
	integer ialv(nlidim)	!line type

	real hnv(nkndim)	!depth of node
	real hev(neldim)	!depth of element
	real hlv(nlidim)	!depth of line

	real xgv(nkndim)	!x coordinate of node
	real ygv(nkndim)	!y coordinate of node

	integer nen3v(3,neldim)	!element index

	integer ipntlv(0:nlidim)!pointer into inodlv
	integer inodlv(nlndim)	!node numbers of lines (dim. nlndim)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale
	save /vscale/

	integer iwhat,ner
	real value

	logical grd_next_line

c--------------------------------------------------------------------
c initialize variables
c--------------------------------------------------------------------

	bstop = .false.
	ner = 6

	nco=0
        nkn=0
        nel=0
        nli=0

	xscale=1.
	yscale=1.
	zscale=1.

	ipntlv(0) = 0

c--------------------------------------------------------------------
c open file or STDIN
c--------------------------------------------------------------------

	call grd_init(file)

c--------------------------------------------------------------------
c loop on lines and read
c--------------------------------------------------------------------

        do while( grd_next_line() )

	  call grd_ival(1,value)
	  iwhat = nint(value)

          if( iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,bstop)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nkndim,bstop
     +				,ipnv,ianv,xgv,ygv,hnv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,neldim,bstop
     +				,ipev,iaev,nen3v,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nlidim,nlndim,bstop
     +				,iplv,ialv,ipntlv,inodlv,hlv)
          else
	  	call rdunknown(iwhat,bstop)
          end if

        end do

c--------------------------------------------------------------------
c end of routine
c--------------------------------------------------------------------

	end

c**********************************************************

	subroutine rdcom(nco,bstop)

c reads comment

	implicit none

	integer nco
	logical bstop

	integer ner,iline
	integer i,n
	character*132 line

	integer ifstch

	call grd_line_info(iline,line)

	i=ifstch(line)
	n=len(line)

	ner = 6

	if(i.gt.0.and.line(i:i).eq.'0') then
	  nco=nco+1
	  if(i.lt.n) then
	    call fempar(line(i+1:))
	  end if
	else if(i.le.0) then
	  !blank line
	else
	  write(ner,*) 'Read error on line ',iline
	  write(ner,*) line
          bstop=.true.
	end if

	end

c**********************************************************

        subroutine rdnode(nkn,nkndim,bstop
     +				,ipnv,ianv,xgv,ygv,hnv)

c reads nodes from .grd file

        implicit none

	integer nkn,nkndim
	logical bstop
        integer ipnv(nkndim)
        integer ianv(nkndim)         !node type
	real xgv(nkndim),ygv(nkndim)
	real hnv(nkndim)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale
	save /vscale/

	logical bread
	integer ner
	integer ianz
        real f(6)
	real depth

	ner = 6

	bread = nkndim .gt. 0		!read nodes?

	call grd_nvals(ianz)
        if(ianz.lt.5) goto 88
        if(ianz.gt.6) goto 87
	call grd_vals(ianz,f)

        nkn=nkn+1
	if( .not. bread ) return
	if(bread .and. nkn.gt.nkndim) then
	  bread = .false.
	  bstop = .true.
	  if( nkn .eq. nkndim+1 ) then		!just one time
	    write(ner,*) 'dimension of nkndim too low: ',nkndim
	  end if
	end if
	if( .not. bread ) return

        ipnv(nkn)=nint(f(2))
        ianv(nkn)=nint(f(3))
        xgv(nkn)=f(4)*xscale
        ygv(nkn)=f(5)*yscale

       	depth = -999.
	if( ianz .ge. 6 ) depth = f(6)*zscale
	hnv(nkn) = depth

        return
   87	continue
	write(ner,*) 'Too much data on line'
	call grd_write_line_info
	bstop=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	bstop=.true.
	return
   99	continue
	write(ner,*) 'dimension of nkndim too low: ',nkndim
	stop 'error stop rdnode: nkndim'
        end

c**********************************************************

        subroutine rdelem(nel,neldim,bstop
     +				,ipev,iaev,nen3v,hev)

c reads elements from .grd file

        implicit none

	integer nel,neldim
	logical bstop
        integer ipev(neldim),iaev(neldim)
	integer nen3v(3,neldim)
	real hev(neldim)

	logical bread
	integer ner,ii
        real f(4)
	integer ilist(10)
	integer inum,itype,ianz
	integer ivert,nvert,istart
	real depth

	ner = 6

	bread = neldim .gt. 0		!read elements?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nel=nel+1
	if(bread .and. nel.gt.neldim) then
	  bread = .false.
	  bstop = .true.
	  if( nel .eq. neldim+1 ) then		!just one time
	    write(ner,*) 'dimension of neldim too low: ',neldim
	  end if
	end if

        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        if(nvert.gt.3) goto 87

        istart=4
        ivert=nvert
	if( .not. bread ) ivert = -ivert

	call read_node_list(ivert,istart,ilist,depth)

	if(ivert.lt.nvert) goto 86

	if( bread ) then
          ipev(nel) = inum
          iaev(nel) = itype
	  hev(nel)  = depth
	  do ii=1,3
	    nen3v(ii,nel) = ilist(ii)
	  end do
	end if

	return
   86	continue
	write(ner,*) 'Could not read all vertices for element'
	write(ner,*) '   internal = ',nel,'    external = ', inum
	call grd_write_line_info
	bstop=.true.
	return
   87	continue
	write(ner,*) 'Not a triangle on line'
	call grd_write_line_info
	bstop=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	bstop=.true.
	return
   99	continue
	write(ner,*) 'dimension of neldim too low: ',neldim
	stop 'error stop rdelem: neldim'
	end

c**********************************************************

        subroutine rdline(nli,nlidim,nlndim,bstop
     +				,iplv,ialv,ipntlv,inodlv,hlv)

c reads lines from .grd file

        implicit none

	integer nli,nlidim,nlndim
	logical bstop
        integer iplv(nlidim),ialv(nlidim)
	integer ipntlv(0:nlidim)		!pointer into inodlv
	integer inodlv(nlndim)			!node numbers of lines
	real hlv(nlidim)

	logical bread
	integer ner
	integer inum,itype,ianz,ipnt
	integer ivert,nvert,istart
        real f(4)
	real depth

	ner = 6

	bread = nlidim .gt. 0		!read lines?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nli=nli+1
	if(bread .and. nli.gt.nlidim) goto 99
        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        istart=4
        ivert=nvert
	if( bread ) then
	  ipnt = ipntlv(nli-1)
	  if( ipnt + nvert .gt. nlndim ) goto 98
	else
	  ipnt = 0
	  ivert = -ivert
	end if

	call read_node_list(ivert,istart,inodlv(ipnt+1),depth)

	if(ivert.lt.nvert) goto 86

	if( bread ) then
          iplv(nli) = inum
          ialv(nli) = itype
	  hlv(nli)  = depth
	  ipntlv(nli) = ipntlv(nli-1) + nvert
	end if

	return
   86	continue
	write(ner,*) 'Could not read all nodes for line ',nli
	call grd_write_line_info
	bstop=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	bstop=.true.
	return
   98	continue
	write(6,*) 'dimension of nlndim too low: ',nlndim
	stop 'error stop rdline: nlndim'
   99	continue
	write(6,*) 'dimension of nlidim too low: ',nlidim
	stop 'error stop rdline: nlidim'
	end

c**********************************************************

	subroutine read_node_list(nvert,istart,nodes,depth)

c reads node list

	implicit none

	integer nvert		!number of vertices to read
	integer istart
	integer nodes(1)
	real depth

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale
	save /vscale/

	logical bread,bline
	integer i,ivert,ianz
	real value

	logical grd_next_line

	bread = nvert .gt. 0
	nvert = abs(nvert)

        i=istart
        ivert=0
	bline = .true.

	call grd_nvals(ianz)

        do while( bline .and. ivert .lt. nvert )    !read vertices

          if( i .ge. ianz ) then			!read new line
	    i = 0
	    bline = grd_next_line()
	    call grd_nvals(ianz)
	  else						!get new node number
            i=i+1
            ivert=ivert+1
	    call grd_ival(i,value)
            if( bread ) nodes(ivert) = nint(value)
	  end if

	end do

       	depth = -999.
	if( i .lt. ianz ) then
	  call grd_ival(i+1,value)
	  depth = value*zscale
	end if

	nvert = ivert		!pass back number of vertices read

	end

c******************************************************************************

        subroutine rdunknown(iwhat,bstop)

c handles unknown type

	logical bstop

	integer iwhat
	integer iline
	character*132 line

	integer ner

	ner = 6

	call grd_line_info(iline,line)

        write(ner,*) 'Type ',iwhat,' at line ',iline,' not recognized'
        write(ner,*) line

        bstop=.true.

	end

c******************************************************************************

	function ifstch(line)

c finds first char of line that is not blank or tab

        implicit none

	integer ifstch
	character*(*) line

	integer length,i

	length=len(line)

	do i=1,length
	  if(line(i:i).ne.' ' .and. line(i:i).ne.'	' ) then
		ifstch=i
		return
	  end if
	end do

	ifstch=-1

	return
	end

c******************************************************************************

	subroutine fempar(line)

c read parameters for fem model 
c
c the following special symbols are recognized on a comment line :
c
c (FEM-TITLE)	title for basin			default: first comment line
c (FEM-SCALE)	scale in x/y/z for basin	default: 1.,1.,1.
c (FEM-LATID)	latitude of basin (degrees)	default: 0.
c (FEM-NORTH)	true north of basin (degrees) 	default: 0.
c
c all entries are optional
c
c example (c of fortran comment must be deleted, line starts with 0) :
c
c 0 (FEM-TITLE) test basin
c 0 (FEM-SCALE) 0.5 0.5 2.
c 0 (FEM-LATID) 45.0
c 0 (FEM-NORTH) 90.0
c
        implicit none

	character*(*) line

	integer i,j,n
	integer ifstch,iscan
	character*80 descrr
	logical btitle
	real f(10)
	real grav,fcor,dcor,dirn
	real xscale,yscale,zscale

	common /descrr/ descrr
	common /pkonst/ grav,fcor,dcor,dirn
	common /vscale/ xscale,yscale,zscale
	save /vscale/

	save btitle
	data btitle /.false./

	i=ifstch(line)
	n=len(line)

	if( i.gt.0 .and. i+10.lt.n ) then
	  if( line(i:i+10) .eq. '(FEM-TITLE)' ) then
		descrr=line(i+11:)
		btitle=.true.
	  else if( line(i:i+10) .eq. '(FEM-SCALE)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.3) then
		  xscale=f(1)
		  yscale=f(2)
		  zscale=f(3)
		else
		  write(6,*) 'error reading (FEM-SCALE) :',j
		  write(6,*) line
		  write(6,*) line(i+11:)
		end if
	  else if( line(i:i+10) .eq. '(FEM-LATID)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.1) then
		  dcor=f(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) line
		end if
	  else if( line(i:i+10) .eq. '(FEM-NORTH)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.1) then
		  dirn=f(1)
		else
		  write(6,*) 'error reading (FEM-NORTH) :'
		  write(6,*) line
		end if
	  end if
	end if

c use first comment as title

	if( i.gt.0 .and. .not.btitle ) then
		descrr=line(i:)
		btitle=.true.
	end if

	return
	end

c******************************************************************************
c******************************************************************************
c******************************************************************************

	subroutine ex2in(nkn,ne,nl,ipnv,ipaux,nen3v,inodlv,bstop)

c changing extern with intern node numbers in elements and lines
c
c if no elements or lines are given, set ne or nl to 0

	implicit none

        integer nkn,ne,nl
        logical bstop
        integer ipnv(nkn)
        integer ipaux(nkn)
        integer nen3v(ne)
        integer inodlv(nl)

	call isort(nkn,ipnv,ipaux)
	call chex2in(nkn,ne,nen3v,ipnv,ipaux,bstop)
	call chex2in(nkn,nl,inodlv,ipnv,ipaux,bstop)

	end

c*****************************************************************
 
        subroutine chex2in(nkn,n,nodes,ipnv,index,bstop)
 
c changing extern with intern node numbers node list
 
        implicit none
 
        integer nkn,n
        logical bstop
        integer nodes(n)
        integer ipnv(nkn)
        integer index(nkn)
 
        integer k,i,kn
	integer ner
        integer locate
 
	ner = 6

        do i=1,n
            kn=nodes(i)
            k=locate(nkn,ipnv,index,kn)
            if(k.le.0) then
              write(ner,*)' warning : node',kn,' not found'
              bstop=.true.
            else
              nodes(i)=k
            end if
        end do
 
        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_init(file)

c initializes reading from grid file

	implicit none

	character*(*) file

	integer ner
	integer ifileo

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz

	nin = 0
	iline = 0
	ianz = 0
	ner = 6

        if( file .eq. ' ' ) then
          nin = 5
          write(ner,*) 'reading GRD file from stdin...'
        else
	  nin = ifileo(nin,file,'formatted','old')
        end if

	if(nin.le.0) then
	  write(ner,*) 'error opening file'
	  write(ner,*) file
	  stop 'error stop grd_init: cannot open file'
	end if

	end

c*****************************************************************

	function grd_next_line()

c reads next line from file

	implicit none

	logical grd_next_line

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz
	real f(80)
	common /grdcom_r/ f
	character*132 line
	common /grdcom_c/ line
	save /grdcom_i/, /grdcom_r/, /grdcom_c/

	integer ner,ios
	integer iscan

	ner = 6
	grd_next_line = .false.

        read(nin,'(a)',iostat=ios) line

        if( ios .lt. 0 ) then
          write(ner,*) iline,' lines read'
	  close(nin)
	  return
        else if( ios .gt. 0 ) then
          write(ner,*) 'read error close to line ',iline
          write(ner,*) line
          stop 'error stop grd_next_line: reading line'
        end if

	iline = iline + 1
        ianz=iscan(line,1,f)
	if( ianz .gt. 80 ) goto 99

	grd_next_line = .true.

	return
   99	continue
	write(ner,*) 'dimension of array f too low ',80,ianz
	stop 'error stop grd_next_line: ianz'
	end

c*****************************************************************

	subroutine grd_nvals(nvals)

c returns number of values on line

	implicit none

	integer nvals

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz

	nvals = ianz

	end

c*****************************************************************

	subroutine grd_vals(nvals,vals)

c returns nvals in vals

	implicit none

	integer nvals
	real vals(nvals)

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz
	real f(80)
	common /grdcom_r/ f

	integer i,n,nmin,nmax

	n = min(ianz,nvals)

	do i=1,n
	  vals(i) = f(i)
	end do

c	if nvals is greater than ianz return zero

	nmin = max(1,n+1)
	nmax = min(80,nvals)

	do i=nmin,nmax
	  vals(i) = 0.
	end do

	end

c*****************************************************************

	subroutine grd_ival(ival,val)

c returns value at position ival

	implicit none

	integer ival
	real val

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz
	real f(80)
	common /grdcom_r/ f

	val = 0.
	if( ival .ge. 1 .and. ival .le. ianz ) val = f(ival)

	end

c*****************************************************************

	subroutine grd_line_info(iline_grd,line_grd)

c returns info on line

	implicit none

	integer iline_grd
	character*(*) line_grd

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz
	character*132 line
	common /grdcom_c/ line

	iline_grd = iline
	line_grd = line

	end

c*****************************************************************

	subroutine grd_write_line_info

c write info on line

	implicit none

	integer nin,iline,ianz
	common /grdcom_i/ nin,iline,ianz
	character*132 line
	common /grdcom_c/ line

	integer ner

	ner = 6

	write(ner,*) 'line number: ',iline
	write(ner,*) line

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_write_grid(
     +			 file
     +			,nco,nkn,nel,nli
     +			,ipnv,ipev,iplv
     +			,ianv,iaev,ialv
     +			,hnv,hev,hlv
     +			,xgv,ygv
     +                  ,nen3v
     +                  ,ipntlv,inodlv
     +			)

c writes grd file
c
c works only with triangles as elements

	implicit none

	character*(*) file	!file name

	integer nco		!total number of comments read
	integer nkn		!total number of nodes read
	integer nel		!total number of elements read
	integer nli		!total number of lines read

	integer ipnv(nkn)	!external node number
	integer ipev(nel)	!external element number
	integer iplv(nli)	!external line number

	integer ianv(nkn) 	!node type
	integer iaev(nel)	!element type
	integer ialv(nli)	!line type

	real hnv(nkn)		!depth of node
	real hev(nel)		!depth of element
	real hlv(nli)		!depth of line

	real xgv(nkn)		!x coordinate of node
	real ygv(nkn)		!y coordinate of node

	integer nen3v(3,nel)	!element index

	integer ipntlv(0:nli)	!pointer into inodlv
	integer inodlv(1)	!node numbers of lines (dim. nlndim)

	integer nout
	integer i,ii
	integer n,ibase
	integer nodes(3)
	real depth

	integer ifileo

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99

	write(nout,*)

	do i=1,nkn
	  depth = hnv(i)
	  if( depth .eq. -999. ) then
	    write(nout,1000) 1,ipnv(i),ianv(i),xgv(i),ygv(i)
	  else
	    write(nout,1000) 1,ipnv(i),ianv(i),xgv(i),ygv(i),depth
	  end if
	end do

	write(nout,*)

	do i=1,nel
	  depth = hev(i)
	  do ii=1,3
	    nodes(ii) = ipnv( nen3v(ii,i) )
	  end do
	  if( depth .eq. -999. ) then
	    write(nout,2000) 2,ipev(i),iaev(i)
     +				,3,(nodes(ii),ii=1,3)
	  else
	    write(nout,2000) 2,ipev(i),iaev(i)
     +				,3,(nodes(ii),ii=1,3),depth
	  end if
	end do

	write(nout,*)

	do i=1,nli
	  depth = hlv(i)
	  n = ipntlv(i) - ipntlv(i-1)
	  ibase = ipntlv(i-1)
	  write(nout,3000) 3,iplv(i),ialv(i),n
	  call grd_write_node_list(nout,n,inodlv(ibase+1),ipnv,depth)
	end do

	write(nout,*)

	close(nout)

	return
 1000	format(i1,2i10,3e16.8)
 2000	format(i1,6i10,e16.8)
 3000	format(i1,3i10)
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop grd_write_grid: cannot open file'
	end

c*****************************************************************

	subroutine grd_write_node_list(nout,n,nodes,ipnv,depth)

c writes out node list

	implicit none

	integer nout
	integer n
	integer nodes(1)
	integer ipnv(1)
	real depth

	integer i,iend,ii,nend
	integer lnodes(6)
	character*20 format

	do i=1,n,6
	  iend = min(n,i+5)
	  nend = iend - i + 1
	  !write(6,*) i,n,iend,nend
	  do ii=1,nend
	    lnodes(ii) = ipnv( nodes(i+ii-1) )
	  end do
	  if( iend .eq. n .and. depth .ne. -999. ) then
	      write(format,'(a3,i1,a9)') '1x,',nend,'i10,e16.8'
	      write(nout,format) (lnodes(ii),ii=1,nend),depth
	  else
	      write(nout,3001) (lnodes(ii),ii=1,nend)
	  end if
	end do

	return
 3001	format(1x,6i10)
	end

c*****************************************************************

