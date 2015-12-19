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
c	subroutine grd_internal_init(file)
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
c 08.01.2015    ggu     common blocks in include file
c 02.10.2015    ggu     new routine is_grd_file()
c 02.10.2015    ggu     now stopping on first error
c
c**********************************************************

!==============================================================
	module grd
!==============================================================

	implicit none

        real, save :: xscale_grd = 1.
        real, save :: yscale_grd = 1.
        real, save :: zscale_grd = 1.
	character*80, save :: title_grd = ' '
	real, save :: dcor_grd = 0.
	real, save :: dirn_grd = 0.

        integer, save :: nin_grd,iline_grd,ianz_grd
        real, save :: f_grd(80)
        character*132, save :: line_grd

	integer, save :: nk_grd = 0
	integer, save :: ne_grd = 0
	integer, save :: nl_grd = 0
	integer, save :: nne_grd = 0
	integer, save :: nnl_grd = 0

	integer, save :: nk_grd_alloc = 0
	integer, save :: ne_grd_alloc = 0
	integer, save :: nl_grd_alloc = 0
	integer, save :: nne_grd_alloc = 0
	integer, save :: nnl_grd_alloc = 0

	integer, save, allocatable :: ippnv(:),ippev(:),ipplv(:)
	integer, save, allocatable :: ianv(:),iaev(:),ialv(:)
        real, save, allocatable :: hhnv(:),hhev(:),hhlv(:)
        real, save, allocatable :: xv(:),yv(:)

        integer, save, allocatable :: ipntev(:),inodev(:)
        integer, save, allocatable :: ipntlv(:),inodlv(:)

!==============================================================
	contains
!==============================================================

	subroutine grd_init(nk,ne,nl,nne,nnl)

	integer nk,ne,nl,nne,nnl

	logical :: bdebug = .false.

	if( bdebug ) write(6,*) 'nk: ',nk,nk_grd,nk_grd_alloc

	if( nk .ne. nk_grd ) then
	!if( nk > nk_grd_alloc ) then
	  if( nk_grd_alloc > 0 ) then
	    deallocate(ippnv,ianv,hhnv,xv,yv)
	  end if
	  nk_grd_alloc = nk
	  if( nk > 0 ) then
	    allocate(ippnv(nk),ianv(nk),hhnv(nk),xv(nk),yv(nk))
	  end if
	end if
	nk_grd = nk

	if( bdebug ) write(6,*) 'ne: ',ne,ne_grd,ne_grd_alloc

	if( ne .ne. ne_grd ) then
	!if( ne > ne_grd_alloc ) then
	  if( ne_grd_alloc > 0 ) then
	    deallocate(ippev,iaev,hhev)
	  end if
	  ne_grd_alloc = ne
	  if( ne > 0 ) then
	    allocate(ippev(ne),iaev(ne),hhev(ne))
	  end if
	end if
	ne_grd = ne
	if( allocated(ipntev) ) deallocate(ipntev)
	allocate(ipntev(0:ne))

	if( bdebug ) write(6,*) 'nl: ',nl,nl_grd,nl_grd_alloc

	if( nl .ne. nl_grd_alloc ) then
	!if( nl > nl_grd_alloc ) then
	  if( nl_grd_alloc > 0 ) then
	    deallocate(ipplv,ialv,hhlv)
	  end if
	  nl_grd_alloc = nl
	  if( nl > 0 ) then
	    allocate(ipplv(nl),ialv(nl),hhlv(nl))
	  end if
	end if
	nl_grd = nl
	if( allocated(ipntlv) ) deallocate(ipntlv)
	allocate(ipntlv(0:nl))

	if( bdebug ) write(6,*) 'nne: ',nne,nne_grd,nne_grd_alloc

	if( nne .ne. nne_grd_alloc ) then
	!if( nne > nne_grd_alloc ) then
	  if( nne_grd_alloc > 0 ) then
	    deallocate(inodev)
	  end if
	  nne_grd_alloc = nne
	  if( nne > 0 ) then
	    allocate(inodev(nne))
	  end if
	end if
	nne_grd = nne

	if( bdebug ) write(6,*) 'nnl: ',nnl,nnl_grd,nnl_grd_alloc

	if( nnl .ne. nnl_grd_alloc ) then
	!if( nnl > nnl_grd_alloc ) then
	  if( nnl_grd_alloc > 0 ) then
	    deallocate(inodlv)
	  end if
	  nnl_grd_alloc = nnl
	  if( nnl > 0 ) then
	    allocate(inodlv(nnl))
	  end if
	end if
	nnl_grd = nnl

	end subroutine grd_init

!==============================================================
	end module grd
!==============================================================

c*****************************************************************

	subroutine grd_get_params(nk,ne,nl,nne,nnl)

	use grd

	integer nk,ne,nl,nne,nnl

	nk = nk_grd
	ne = ne_grd
	nl = nl_grd
	nne = nne_grd
	nnl = nnl_grd

	end subroutine grd_get_params

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_read(file)

	use grd

	implicit none

	character(*) file
	integer ner,nco
	logical bstop
	integer nk,ne,nl,nne,nnl
	integer nkndi,neldi,nlidi,nendi,nlndi

        call grd_info(file,nk,ne,nl,nne,nnl)

	write(6,*) 'grd_info: ',nk,ne,nl,nne,nnl

	nkndi = nk
	neldi = ne
	nlidi = nl
	nendi = nne
	nlndi = nnl

	call grd_init(nkndi,neldi,nlidi,nendi,nlndi)

        ner = 6
        bstop = .false.

	call rdgrd(
     +			 file
     +			,bstop
     +			,nco,nk,ne,nl,nne,nnl
     +			,nkndi,neldi,nlidi,nendi,nlndi
     +			,ippnv,ippev,ipplv
     +			,ianv,iaev,ialv
     +			,hhnv,hhev,hhlv
     +			,xv,yv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)

	call ex2in(nk,nne,nnl,ippnv,inodev,inodlv,bstop)

	if( bstop ) stop 'error stop grd_read: error reading grd file'

	write(6,*) 'total number of lines read: ',iline_grd

	end

c*****************************************************************

        subroutine grd_info(gfile,nk,ne,nl,nne,nnl)

c reads grd file to obtain basic parameters

        implicit none

        character*(*) gfile
        integer nk              !total number of nodes
        integer ne              !total number of elems
        integer nl              !total number of lines
        integer nne             !total number of nodes in elems
        integer nnl             !total number of nodes in lines

        logical bstop
        integer ner,nco
        integer nkndi0,neldi0,nlidi0,nendi0,nlndi0

	integer ippnv(1),ippev(1),ipplv(1)
	integer ianv(1),iaev(1),ialv(1)
        real hhnv(1),hhev(1),hhlv(1)
        real xv(1),yv(1)

        integer ipntev0(0:0),inodev(1)
        integer ipntlv0(0:0),inodlv(1)

c-----------------------------------------------------------------
c initialize parameters
c-----------------------------------------------------------------

        ner = 6
        bstop = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

c-----------------------------------------------------------------
c read grd file
c-----------------------------------------------------------------

        call rdgrd(
     +                   gfile
     +                  ,bstop
     +                  ,nco,nk,ne,nl,nne,nnl
     +                  ,nkndi0,neldi0,nlidi0,nendi0,nlndi0
     +                  ,ippnv,ippev,ipplv
     +                  ,ianv,iaev,ialv
     +                  ,hhnv,hhev,hhlv
     +                  ,xv,yv
     +                  ,ipntev0,inodev
     +                  ,ipntlv0,inodlv
     +                  )

        end

c**********************************************************

	function is_grd_file(gfile)

	implicit none

	logical is_grd_file
        character*(*) gfile

        logical bstop
        integer nk,ne,nl,nne,nnl
        integer ner,nco
        integer nkndi0,neldi0,nlidi0,nendi0,nlndi0

	integer ippnv(1),ippev(1),ipplv(1)
	integer ianv(1),iaev(1),ialv(1)
        real hhnv(1),hhev(1),hhlv(1)
        real xv(1),yv(1)

        integer ipntev0(0:0),inodev(1)
        integer ipntlv0(0:0),inodlv(1)

c-----------------------------------------------------------------
c initialize parameters
c-----------------------------------------------------------------

        ner = 6
        bstop = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

c-----------------------------------------------------------------
c read grd file
c-----------------------------------------------------------------

        call rdgrd(
     +                   gfile
     +                  ,bstop
     +                  ,nco,nk,ne,nl,nne,nnl
     +                  ,nkndi0,neldi0,nlidi0,nendi0,nlndi0
     +                  ,ippnv,ippev,ipplv
     +                  ,ianv,iaev,ialv
     +                  ,hhnv,hhev,hhlv
     +                  ,xv,yv
     +                  ,ipntev0,inodev
     +                  ,ipntlv0,inodlv
     +                  )

	is_grd_file = .not. bstop

        end

c**********************************************************

	subroutine grd_close

	use grd

	implicit none

	call grd_init(0,0,0,0,0)

	end

c**********************************************************

	subroutine grd_write(file)

	use grd

	implicit none

	character*(*) file

	integer nk,ne,nl,nne,nnl
	integer nco
	integer i,iu
	logical bdebug

	bdebug = .false.
	nco = 0

	call grd_get_params(nk,ne,nl,nne,nnl)

	if( bdebug ) then
	iu = 91
	write(iu,*) '========================================'
	write(iu,*) 'grd_write: ',nk,ne,nl,nne,nnl
	write(iu,*) trim(file)
	write(iu,*)
	write(iu,'(5i10)') (ippnv(i),i=1,nk)
	write(iu,*)
	write(iu,'(5i10)') (ippev(i),i=1,ne)
	write(iu,*)
	write(iu,'(5i10)') (ipntev(i),i=0,ne)
	write(iu,*)
	write(iu,'(5i10)') (inodev(i),i=1,nne)
	write(iu,*)
	write(iu,*) 'grd_write: ',nk,ne,nl,nne,nnl
	write(iu,*) '========================================'
	end if

	call grd_write_grid(
     +			 file
     +			,nco,nk,ne,nl,nne,nnl
     +			,ippnv,ippev,ipplv
     +			,ianv,iaev,ialv
     +			,hhnv,hhev,hhlv
     +			,xv,yv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine rdgrd(
     +			 file
     +			,bstop
     +			,nco,nkn,nel,nli,nne,nnl
     +			,nknddi,nelddi,nliddi,nenddi,nlnddi
     +			,ipnv,ipev,iplv
     +			,ianv,iaev,ialv
     +			,hnv,hev,hlv
     +			,xgv,ygv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)

c reads grd file
c
c after read nen3v and inodlv contain still external node numbers
c use ex2in to convert external to internal node numbers
c
c works only with triangles as elements

	!use grd

	implicit none

	character*(*) file	!file name

	logical bstop		!true on return if error

	integer nco		!total number of comments read
	integer nkn		!total number of nodes read
	integer nel		!total number of elements read
	integer nli		!total number of lines read
	integer nne		!total number of nodes in elems
	integer nnl		!total number of nodes in lines

	integer nknddi		!dimension for number of nodes
	integer nelddi		!dimension for number of elements
	integer nliddi		!dimension for number of lines
	integer nenddi		!dimension for node numbers of elems
	integer nlnddi		!dimension for node numbers of lines

	integer ipnv(nknddi)	!external node number
	integer ipev(nelddi)	!external element number
	integer iplv(nliddi)	!external line number

	integer ianv(nknddi) 	!node type
	integer iaev(nelddi)	!element type
	integer ialv(nliddi)	!line type

	real hnv(nknddi)	!depth of node
	real hev(nelddi)	!depth of element
	real hlv(nliddi)	!depth of line

	real xgv(nknddi)	!x coordinate of node
	real ygv(nknddi)	!y coordinate of node

	integer ipntev(0:nelddi)!pointer into inodev
	integer inodev(nenddi)	!node numbers of elems
	integer ipntlv(0:nliddi)!pointer into inodlv
	integer inodlv(nlnddi)	!node numbers of lines

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
	nne=0
	nnl=0

	ipntev(0) = 0
	ipntlv(0) = 0

c--------------------------------------------------------------------
c open file or STDIN
c--------------------------------------------------------------------

	call grd_internal_init(file)

c--------------------------------------------------------------------
c loop on lines and read
c--------------------------------------------------------------------

        do while( grd_next_line() )

	  call grd_ival(1,value)
	  iwhat = nint(value)

          if( iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,bstop)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nknddi,bstop
     +				,ipnv,ianv,xgv,ygv,hnv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,nne,nelddi,nenddi,bstop
     +				,ipev,iaev,ipntev,inodev,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nnl,nliddi,nlnddi,bstop
     +				,iplv,ialv,ipntlv,inodlv,hlv)
          else
	  	call rdunknown(iwhat,bstop)
          end if

	  if( bstop ) exit
        end do

	call grd_internal_close

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

	if(i.le.0) then
	  !blank line
	else if(line(i:i).eq.'0') then
	  nco=nco+1
	  if(i.lt.n) then
	    call fempar(line(i+1:))
	  end if
	else
	  write(ner,*) 'Read error on line ',iline
	  write(ner,*) line
          bstop=.true.
	end if

	end

c**********************************************************

        subroutine rdnode(nkn,nknddi,bstop
     +				,ipnv,ianv,xgv,ygv,hnv)

c reads nodes from .grd file

	use grd, only : xscale_grd,yscale_grd,zscale_grd

        implicit none

	integer nkn,nknddi
	logical bstop
        integer ipnv(nknddi)
        integer ianv(nknddi)         !node type
	real xgv(nknddi),ygv(nknddi)
	real hnv(nknddi)

	logical bread
	integer ner
	integer ianz
        real f(6)
	real depth

	ner = 6

	bread = nknddi .gt. 0		!read nodes?

	call grd_nvals(ianz)
        if(ianz.lt.5) goto 88
        if(ianz.gt.6) goto 87
	call grd_vals(ianz,f)

        nkn=nkn+1
	if( .not. bread ) return
	if(bread .and. nkn.gt.nknddi) then
	  bread = .false.
	  bstop = .true.
	  if( nkn .eq. nknddi+1 ) then		!just one time
	    write(ner,*) 'dimension of nknddi too low: ',nknddi
	  end if
	end if
	if( .not. bread ) return

        ipnv(nkn)=nint(f(2))
        ianv(nkn)=nint(f(3))
        xgv(nkn)=f(4)*xscale_grd
        ygv(nkn)=f(5)*yscale_grd

       	depth = -999.
	if( ianz .ge. 6 ) depth = f(6)*zscale_grd
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
	write(ner,*) 'dimension of nknddi too low: ',nknddi
	stop 'error stop rdnode: nknddi'
        end

c**********************************************************

        subroutine rdelem(nel,nne,nelddi,nenddi,bstop
     +				,ipev,iaev,ipntev,inodev,hev)

c reads elements from .grd file

        implicit none

	integer nel,nne,nelddi,nenddi
	logical bstop
        integer ipev(nelddi),iaev(nelddi)
	integer ipntev(0:nelddi)
	integer inodev(max(1,nenddi))
	real hev(nelddi)

	logical bread
	integer ner,ii
        real f(4)
	integer ilist(10)
	integer inum,itype,ianz,ipnt
	integer ivert,nvert,istart
	real depth

	ner = 6

	bread = nelddi .gt. 0		!read elements?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nel=nel+1
	if(bread .and. nel.gt.nelddi) then
	  bread = .false.
	  bstop = .true.
	  if( nel .eq. nelddi+1 ) then		!just one time
	    write(ner,*) 'dimension of nelddi too low: ',nelddi
	  end if
	end if

        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        if(nvert.gt.3) goto 87

        istart=4
        ivert=nvert
	if( .not. bread ) ivert = -ivert

        istart=4
        ivert=nvert
	if( bread ) then
	  ipnt = ipntev(nel-1)
	  if( ipnt + nvert .gt. nenddi ) goto 98
	else
	  ipnt = 0
	  ivert = -ivert
	end if

	call read_node_list(ivert,istart,inodev(ipnt+1),depth)

	if(ivert.lt.nvert) goto 86

	if( bread ) then
          ipev(nel) = inum
          iaev(nel) = itype
	  hev(nel)  = depth
	  ipntev(nel) = ipntev(nel-1) + nvert
	end if
	nne = nne + nvert

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
   98   continue
        write(6,*) 'dimension of nenddi too low: ',nenddi
        stop 'error stop rdelem: nenddi'
   99	continue
	write(ner,*) 'dimension of nelddi too low: ',nelddi
	stop 'error stop rdelem: nelddi'
	end

c**********************************************************

        subroutine rdline(nli,nnl,nliddi,nlnddi,bstop
     +				,iplv,ialv,ipntlv,inodlv,hlv)

c reads lines from .grd file

        implicit none

	integer nli,nnl,nliddi,nlnddi
	logical bstop
        integer iplv(nliddi),ialv(nliddi)
	integer ipntlv(0:nliddi)		!pointer into inodlv
	integer inodlv(max(1,nlnddi))
	real hlv(nliddi)

	logical bread
	integer ner
	integer inum,itype,ianz,ipnt
	integer ivert,nvert,istart
        real f(4)
	real depth

	ner = 6

	bread = nliddi .gt. 0		!read lines?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nli=nli+1
	if(bread .and. nli.gt.nliddi) goto 99
        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        istart=4
        ivert=nvert
	if( bread ) then
	  ipnt = ipntlv(nli-1)
	  if( ipnt + nvert .gt. nlnddi ) goto 98
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
	nnl = nnl + nvert

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
	write(6,*) 'dimension of nlnddi too low: ',nlnddi
	stop 'error stop rdline: nlnddi'
   99	continue
	write(6,*) 'dimension of nliddi too low: ',nliddi
	stop 'error stop rdline: nliddi'
	end

c**********************************************************

	subroutine read_node_list(nvert,istart,nodes,depth)

c reads node list

	use grd, only : zscale_grd

	implicit none

	integer nvert		!number of vertices to read
	integer istart
	integer nodes(abs(nvert))
	real depth

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
	  depth = value*zscale_grd
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

	subroutine fempar(gline)

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
	use grd
	!use basin

        implicit none

	character*(*) gline

	integer i,j,n
	integer ifstch,iscan
	logical btitle

	save btitle
	data btitle /.false./

	i=ifstch(gline)
	n=len(gline)

	if( i.gt.0 .and. i+10.lt.n ) then
	  if( gline(i:i+10) .eq. '(FEM-TITLE)' ) then
		title_grd=gline(i+11:)
		btitle=.true.
	  else if( gline(i:i+10) .eq. '(FEM-SCALE)' ) then
		j=iscan(gline(i+11:),1,f_grd)
		if(j.eq.3) then
		  xscale_grd=f_grd(1)
		  yscale_grd=f_grd(2)
		  zscale_grd=f_grd(3)
		else
		  write(6,*) 'error reading (FEM-SCALE) :',j
		  write(6,*) gline
		  write(6,*) gline(i+11:)
		end if
	  else if( gline(i:i+10) .eq. '(FEM-LATID)' ) then
		j=iscan(gline(i+11:),1,f_grd)
		if(j.eq.1) then
		  dcor_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) gline
		end if
	  else if( gline(i:i+10) .eq. '(FEM-NORTH)' ) then
		j=iscan(gline(i+11:),1,f_grd)
		if(j.eq.1) then
		  dirn_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-NORTH) :'
		  write(6,*) gline
		end if
	  end if
	end if

c use first comment as title

	if( i.gt.0 .and. .not.btitle ) then
		title_grd=gline(i:)
		btitle=.true.
	end if

	end

c******************************************************************************
c******************************************************************************
c******************************************************************************

	subroutine ex2in(nkn,ne,nl,ippnv,inodev,inodlv,bstop)

c changing extern with intern node numbers in elements and lines
c
c if no elements or lines are given, set ne or nl to 0

	implicit none

        integer nkn,ne,nl
        logical bstop
        integer ippnv(nkn)
        integer inodev(ne)
        integer inodlv(nl)

        integer ipaux(nkn)	!local

	call isort(nkn,ippnv,ipaux)
	call chex2in(nkn,ne,inodev,ippnv,ipaux,bstop)
	call chex2in(nkn,nl,inodlv,ippnv,ipaux,bstop)

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

	subroutine grd_internal_init(file)

c initializes reading from grid file

	use grd

	implicit none

	character*(*) file

	integer ner
	integer ifileo

	nin_grd = 0
	iline_grd = 0
	ianz_grd = 0
	ner = 6

	xscale_grd=1.
	yscale_grd=1.
	zscale_grd=1.

        if( file .eq. ' ' ) then
          write(ner,*) 'no file given...'
          write(ner,*) 'cannot read GRD file from stdin...'
	  stop 'error stop grd_internal_init: cannot read file'
        else
	  nin_grd = ifileo(nin_grd,file,'formatted','old')
        end if

	if(nin_grd.le.0) then
	  write(ner,*) 'error opening file'
	  write(ner,*) file
	  stop 'error stop grd_internal_init: cannot open file'
	end if

	end

c*****************************************************************

	subroutine grd_internal_close

	use grd

	implicit none

	if( nin_grd .ne. 5 ) close(nin_grd)
	nin_grd = 0

	end

c*****************************************************************

	function grd_next_line()

c reads next line from file

	use grd

	implicit none

	logical grd_next_line

	integer ner,ios
	integer iscan

	ner = 6
	grd_next_line = .false.

        read(nin_grd,'(a)',iostat=ios) line_grd

        if( ios .lt. 0 ) then
          !write(ner,*) iline,' lines read'
	  close(nin_grd)
	  return
        else if( ios .gt. 0 ) then
          write(ner,*) 'read error close to line ',iline_grd
          write(ner,*) line_grd
          stop 'error stop grd_next_line: reading line'
        end if

	iline_grd = iline_grd + 1
        ianz_grd=iscan(line_grd,1,f_grd)
	if( ianz_grd .gt. 80 ) goto 99

	grd_next_line = .true.

	return
   99	continue
	write(ner,*) 'dimension of array f too low ',80,ianz_grd
	stop 'error stop grd_next_line: ianz_grd'
	end

c*****************************************************************

	subroutine grd_nvals(nvals)

c returns number of values on line

	use grd

	implicit none

	integer nvals

	nvals = ianz_grd

	end

c*****************************************************************

	subroutine grd_vals(nvals,vals)

c returns nvals in vals

	use grd

	implicit none

	integer nvals
	real vals(nvals)

	integer i,n,nmin,nmax

	n = min(ianz_grd,nvals)

	do i=1,n
	  vals(i) = f_grd(i)
	end do

c	if nvals is greater than ianz_grd return zero

	nmin = max(1,n+1)
	nmax = min(80,nvals)

	do i=nmin,nmax
	  vals(i) = 0.
	end do

	end

c*****************************************************************

	subroutine grd_ival(ival,val)

c returns value at position ival

	use grd

	implicit none

	integer ival
	real val

	val = 0.
	if( ival .ge. 1 .and. ival .le. ianz_grd ) val = f_grd(ival)

	end

c*****************************************************************

	subroutine grd_line_info(iline_gr,line_gr)

c returns info on line

	use grd

	implicit none

	integer iline_gr
	character*(*) line_gr

	iline_gr = iline_grd
	line_gr = line_grd

	end

c*****************************************************************

	subroutine grd_write_line_info

c write info on line

	use grd

	implicit none

	integer ner

	ner = 6

	write(ner,*) 'line number: ',iline_grd
	write(ner,*) line_grd

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_write_grid(
     +			 file
     +			,nco,nk,ne,nl,nne,nnl
     +			,ippnv,ippev,ipplv
     +			,ianv,iaev,ialv
     +			,hhnv,hhev,hhlv
     +			,xv,yv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)

c writes grd file

	implicit none

	character*(*) file	!file name

	integer nco		!total number of comments read
	integer nk		!total number of nodes read
	integer ne		!total number of elements read
	integer nl		!total number of lines read
	integer nne		!total number of nodes in elems 
	integer nnl		!total number of nodes in lines

	integer ippnv(nk)	!external node number
	integer ippev(ne)	!external element number
	integer ipplv(nl)	!external line number

	integer ianv(nk) 	!node type
	integer iaev(ne)	!element type
	integer ialv(nl)	!line type

	real hhnv(nk)		!depth of node
	real hhev(ne)		!depth of element
	real hhlv(nl)		!depth of line

	real xv(nk)		!x coordinate of node
	real yv(nk)		!y coordinate of node

	integer ipntev(0:ne)	!pointer into inodev
	integer inodev(nne)	!node numbers of elems
	integer ipntlv(0:nl)	!pointer into inodlv
	integer inodlv(nnl)	!node numbers of lines

	integer nout
	integer i
	integer n,ib,nmax
	integer k,ie,il
	integer ia,ieext,ilext
	real depth
	logical bsort,bextern
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)

	integer ifileo

	bsort = .true.		!sort nodes on external numbering
	bsort = .false.

	bextern = .false.
	bextern = .true.	!use external nodes

	if( .not. bextern ) bsort = .false.

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99

	nmax = max(nk,ne,nl)
	allocate(ipdex(nmax))
	do i=1,nmax
	  ipdex(i) = i
	end do

	allocate(nextern(nk))
	do k=1,nk
	  if( bextern ) then
	    nextern(k) = ippnv(k)		!use external numbering
	  else
	    nextern(k) = k			!use internal numbering
	  end if
	end do

	write(nout,*)

	if( bsort ) call isort(nk,ippnv,ipdex)

	do i=1,nk
	  k = ipdex(i)
	  depth = hhnv(k)
	  if( depth .eq. -999. ) then
	    write(nout,1000) 1,nextern(k),ianv(k),xv(k),yv(k)
	  else
	    write(nout,1000) 1,nextern(k),ianv(k),xv(k),yv(k),depth
	  end if
	end do

	write(nout,*)

	if( bsort ) call isort(ne,ippev,ipdex)

	do i=1,ne
	  ie = ipdex(i)
	  ieext = ie
	  if( bextern ) ieext = ippev(ie)
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  call grd_write_item(nout,2,ieext,ia,n,
     +				inodev(ib+1),nextern,depth)
	end do

	write(nout,*)

	if( bsort ) call isort(nl,ipplv,ipdex)

	do i=1,nl
	  il = ipdex(i)
	  ilext = il
	  if( bextern ) ilext = ipplv(il)
	  ia = ialv(il)
	  depth = hhlv(il)
	  n = ipntlv(il) - ipntlv(il-1)
	  ib = ipntlv(il-1)
	  call grd_write_item(nout,3,ilext,ia,n,
     +				inodlv(ib+1),nextern,depth)
	end do

	write(nout,*)

	close(nout)

	deallocate(ipdex)
	deallocate(nextern)

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

	subroutine grd_write_item(nout,iwhat,number,itype,n,
     +				nodes,extern,depth)

	implicit none

	integer nout
	integer iwhat
	integer number
	integer itype
	integer n
	integer nodes(n)
	integer extern(*)
	real depth

	integer i
	real, parameter :: flag = -999.

	if( n > 3 ) then
	  write(nout,3000) iwhat,number,itype,n
	  call grd_write_node_list(nout,n,nodes,extern,depth)
	else if( depth == flag ) then
	  write(nout,2000) iwhat,number,itype,n,
     +			(extern(nodes(i)),i=1,n)
	else
	  write(nout,2000) iwhat,number,itype,n,
     +			(extern(nodes(i)),i=1,n),depth
	end if

	return
 2000	format(i1,6i10,e16.8)
 3000	format(i1,3i10)
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
	character*30 format

	do i=1,n,6
	  iend = min(n,i+5)
	  nend = iend - i + 1
	  !write(6,*) i,n,iend,nend
	  do ii=1,nend
	    lnodes(ii) = ipnv( nodes(i+ii-1) )
	  end do
	  if( iend .eq. n .and. depth .ne. -999. ) then
	      write(format,'(a4,i1,a10)') '(1x,',nend,'i10,e16.8)'
	      write(nout,format) (lnodes(ii),ii=1,nend),depth
	  else
	      write(nout,3001) (lnodes(ii),ii=1,nend)
	  end if
	end do

	return
 3001	format(1x,6i10)
	end

c*****************************************************************

	subroutine grd_get_depth(nk,ne,hkv,hev)

	use grd

	implicit none

	integer nk,ne
	real hkv(nk)
	real hev(ne)

	if( ne > 0 ) then
	  if( ne .ne. ne_grd ) then
	    write(6,*) 'ne,ne_grd: ',ne,ne_grd
	    stop 'error stop grd_get_depth: dimension mismatch'
	  end if
	  hev = hhev
	end if

	if( nk > 0 ) then
	  if( nk .ne. nk_grd ) then
	    write(6,*) 'nk,nk_grd: ',nk,nk_grd
	    stop 'error stop grd_get_depth: dimension mismatch'
	  end if
	  hkv = hhnv
	end if

	end

c*****************************************************************

	subroutine grd_get_nodes(np,xp,yp,hp)

	use grd

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real hp(np)

	if( np < nk_grd ) then
	  write(6,*) 'np,nk_grd: ',np,nk_grd
	  stop 'error stop grd_get_nodes: np'
	end if

	np = nk_grd
	xp(1:nk_grd) = xv(1:nk_grd)
	yp(1:nk_grd) = yv(1:nk_grd)
	hp(1:nk_grd) = hhnv(1:nk_grd)

	end

c*****************************************************************

