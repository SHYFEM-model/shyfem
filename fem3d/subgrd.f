
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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
c	subroutine rdunknown(iwhat,berr)
c				handles unknown type
c
c       function ifstch(line)
c				finds first char of line that is not blank
c       subroutine fempar(line)
c				read parameters for fem model
c
c	subroutine ex2in(nkn,ne,nl,ipnv,ipaux,nen3v,inodlv,berr)
c				changing extern with intern node numbers
c	subroutine chex2in(nkn,n,nodes,ipnv,index,berr)
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
c 06.06.2016    ggu     bstop substituted with berr, new accessor routines
c 10.02.2017    ggu     bug fix: do not allocate at least 1 array element
c 06.07.2018    ggu     new handling of line reading routines: read_all_lines()
c 25.10.2018    ccf     grid output in gr3 and msh formats
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
        real, save :: f_grd(81)
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

	logical, save :: berror = .true.	!writes error if found

!==============================================================
	contains
!==============================================================

	subroutine grd_init(nkk,nee,nll,nnee,nnll)

	integer nkk,nee,nll,nnee,nnll

	integer nk,ne,nl,nne,nnl

	logical :: bdebug = .false.
	logical :: balloc

	if( nkk == 0 .and. nee == 0 .and. nll == 0 .and.
     +				nnee == 0 .and. nnll == 0 ) then
	  balloc = .false.
	else
	  balloc = .true.
	end if

	nk = nkk
	ne = nee
	nl = nll
	nne = nnee
	nnl = nnll

	if( bdebug ) write(6,*) 'nk: ',nk,nk_grd,nk_grd_alloc

	if( nk .ne. nk_grd ) then
	  if( allocated(ippnv) ) then
	    deallocate(ippnv,ianv,hhnv,xv,yv)
	  end if
	  nk_grd_alloc = nk
	  if( balloc ) then
	    allocate(ippnv(nk),ianv(nk),hhnv(nk),xv(nk),yv(nk))
	  end if
	end if
	nk_grd = nk

	if( bdebug ) write(6,*) 'ne: ',ne,ne_grd,ne_grd_alloc

	if( ne .ne. ne_grd ) then
	  if( allocated(ippev) ) then
	    deallocate(ippev,iaev,hhev)
	  end if
	  ne_grd_alloc = ne
	  if( balloc ) then
	    allocate(ippev(ne),iaev(ne),hhev(ne))
	  end if
	end if
	ne_grd = ne
	if( allocated(ipntev) ) deallocate(ipntev)
	allocate(ipntev(0:ne))
	ipntev = 0

	if( bdebug ) write(6,*) 'nl: ',nl,nl_grd,nl_grd_alloc

	if( nl .ne. nl_grd_alloc ) then
	  if( allocated(ipplv) ) then
	    deallocate(ipplv,ialv,hhlv)
	  end if
	  nl_grd_alloc = nl
	  if( balloc ) then
	    allocate(ipplv(nl),ialv(nl),hhlv(nl))
	  end if
	end if
	nl_grd = nl
	if( allocated(ipntlv) ) deallocate(ipntlv)
	allocate(ipntlv(0:nl))
	ipntlv = 0

	if( bdebug ) write(6,*) 'nne: ',nne,nne_grd,nne_grd_alloc

	if( nne .ne. nne_grd_alloc ) then
	  if( allocated(inodev) ) then
	    deallocate(inodev)
	  end if
	  nne_grd_alloc = nne
	  if( balloc ) then
	    allocate(inodev(nne))
	  end if
	end if
	nne_grd = nne

	if( bdebug ) write(6,*) 'nnl: ',nnl,nnl_grd,nnl_grd_alloc

	if( nnl .ne. nnl_grd_alloc ) then
	  if( allocated(inodlv) ) then
	    deallocate(inodlv)
	  end if
	  nnl_grd_alloc = nnl
	  if( balloc ) then
	    allocate(inodlv(nnl))
	  end if
	end if
	nnl_grd = nnl

	end subroutine grd_init

!==============================================================
	end module grd
!==============================================================

c*****************************************************************

	subroutine grd_set_error(berr)

	use grd

	implicit none

	logical berr

	berror = berr

	end subroutine grd_set_error

c*****************************************************************

	function grd_write_error()

	use grd

	implicit none

	logical grd_write_error

	grd_write_error = berror

	end function grd_write_error

c*****************************************************************

	subroutine grd_init_fake

	use grd

	implicit none

	call grd_init(1,1,1,1,1)

	end

c*****************************************************************
c*****************************************************************
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
	logical berr,bwrite
	integer nk,ne,nl,nne,nnl
	integer nkndi,neldi,nlidi,nendi,nlndi

        ner = 6
        berr = .false.
	bwrite = .false.

        call grd_info(file,nk,ne,nl,nne,nnl,berr)
	if( berr ) goto 99

	if( bwrite ) write(6,*) 'grd_info: ',nk,ne,nl,nne,nnl

	nkndi = nk
	neldi = ne
	nlidi = nl
	nendi = nne
	nlndi = nnl

	call grd_init(nkndi,neldi,nlidi,nendi,nlndi)

	call rdgrd(
     +			 file
     +			,berr
     +			,nco,nk,ne,nl,nne,nnl
     +			,nkndi,neldi,nlidi,nendi,nlndi
     +			,ippnv,ippev,ipplv
     +			,ianv,iaev,ialv
     +			,hhnv,hhev,hhlv
     +			,xv,yv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)
	if( berr ) goto 99

	call ex2in(nk,nne,nnl,ippnv,inodev,inodlv,berr)
	if( berr ) goto 99

	if( bwrite ) then
	  write(6,*) 'total number of lines read: ',iline_grd
	end if

	return
   99	continue
	stop 'error stop grd_read: error reading grd file'
	end

c*****************************************************************

        subroutine grd_info(gfile,nk,ne,nl,nne,nnl,berr)

c reads grd file to obtain basic parameters

        implicit none

        character*(*) gfile
        integer nk              !total number of nodes
        integer ne              !total number of elems
        integer nl              !total number of lines
        integer nne             !total number of nodes in elems
        integer nnl             !total number of nodes in lines
        logical berr		!true if error

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
        berr = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

	call grd_init_fake

c-----------------------------------------------------------------
c read grd file
c-----------------------------------------------------------------

        call rdgrd(
     +                   gfile
     +                  ,berr
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

        logical berr
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
        berr = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

	call grd_set_error(.false.)	!do not write errors to terminal

	call grd_init_fake

c-----------------------------------------------------------------
c read grd file
c-----------------------------------------------------------------

        call rdgrd(
     +                   gfile
     +                  ,berr
     +                  ,nco,nk,ne,nl,nne,nnl
     +                  ,nkndi0,neldi0,nlidi0,nendi0,nlndi0
     +                  ,ippnv,ippev,ipplv
     +                  ,ianv,iaev,ialv
     +                  ,hhnv,hhev,hhlv
     +                  ,xv,yv
     +                  ,ipntev0,inodev
     +                  ,ipntlv0,inodlv
     +                  )

	call grd_set_error(.true.)	!original behavior

	is_grd_file = .not. berr

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

	subroutine gr3_write(file)

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
	write(iu,*) 'gr3_write: ',nk,ne,nl,nne,nnl
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

	call gr3_write_grid(
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

	subroutine msh_write(file)

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
	write(iu,*) 'msh_write: ',nk,ne,nl,nne,nnl
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

	call msh_write_grid(
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
     +			,berr
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

	logical berr		!true on return if error

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

	logical berrwrite
	integer iwhat,ner
	real value

	logical grd_next_line

c--------------------------------------------------------------------
c initialize variables
c--------------------------------------------------------------------

	
	berrwrite = berr	!write read errors to terminal
	berr = .false.
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

	  berr = berrwrite
          if( iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,berr)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nknddi,berr
     +				,ipnv,ianv,xgv,ygv,hnv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,nne,nelddi,nenddi,berr
     +				,ipev,iaev,ipntev,inodev,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nnl,nliddi,nlnddi,berr
     +				,iplv,ialv,ipntlv,inodlv,hlv)
          else
	  	call rdunknown(iwhat,berr)
          end if

	  if( berr ) exit
        end do

	call grd_internal_close

	if( nkn == 0 ) berr = .true.	!we must really read something

c--------------------------------------------------------------------
c end of routine
c--------------------------------------------------------------------

	end

c**********************************************************

	subroutine rdcom(nco,berr)

c reads comment

	implicit none

	integer nco
	logical berr

	integer ner,iline
	integer i,n
	character*132 line

	integer ifstch
	logical grd_write_error

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
	  if( grd_write_error() ) then
	    write(ner,*) 'Read error on line ',iline
	    write(ner,*) line
	  end if
          berr=.true.
	end if

	end

c**********************************************************

        subroutine rdnode(nkn,nknddi,berr
     +				,ipnv,ianv,xgv,ygv,hnv)

c reads nodes from .grd file

	use grd, only : xscale_grd,yscale_grd,zscale_grd

        implicit none

	integer nkn,nknddi
	logical berr
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
	  berr = .true.
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
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
	return
   99	continue
	write(ner,*) 'dimension of nknddi too low: ',nknddi
	stop 'error stop rdnode: nknddi'
        end

c**********************************************************

        subroutine rdelem(nel,nne,nelddi,nenddi,berr
     +				,ipev,iaev,ipntev,inodev,hev)

c reads elements from .grd file

        implicit none

	integer nel,nne,nelddi,nenddi
	logical berr
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
	  berr = .true.
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
	berr=.true.
	return
   87	continue
	write(ner,*) 'Not a triangle on line'
	call grd_write_line_info
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
	return
   98   continue
        write(6,*) 'dimension of nenddi too low: ',nenddi
        stop 'error stop rdelem: nenddi'
   99	continue
	write(ner,*) 'dimension of nelddi too low: ',nelddi
	stop 'error stop rdelem: nelddi'
	end

c**********************************************************

        subroutine rdline(nli,nnl,nliddi,nlnddi,berr
     +				,iplv,ialv,ipntlv,inodlv,hlv)

c reads lines from .grd file

        implicit none

	integer nli,nnl,nliddi,nlnddi
	logical berr
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
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
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

        subroutine rdunknown(iwhat,berr)

c handles unknown type

	logical berr

	integer iwhat
	integer iline
	character*132 line

	integer ner

	ner = 6

	if( berr ) then
	  call grd_line_info(iline,line)
          write(ner,*) 'Type ',iwhat,' at line ',iline,' not recognized'
          write(ner,*) line
	end if

        berr=.true.

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
	integer ifstch,iscanf
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
		j=iscanf(gline(i+11:),f_grd,4)
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
		j=iscanf(gline(i+11:),f_grd,2)
		if(j.eq.1) then
		  dcor_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) gline
		end if
	  else if( gline(i:i+10) .eq. '(FEM-NORTH)' ) then
		j=iscanf(gline(i+11:),f_grd,2)
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

	subroutine ex2in(nkn,ne,nl,ippnv,inodev,inodlv,berr)

c changing extern with intern node numbers in elements and lines
c
c if no elements or lines are given, set ne or nl to 0

	implicit none

        integer nkn,ne,nl
        logical berr
        integer ippnv(nkn)
        integer inodev(ne)
        integer inodlv(nl)

        integer ipaux(nkn)	!local

	call isort(nkn,ippnv,ipaux)
	call chex2in(nkn,ne,inodev,ippnv,ipaux,berr)
	call chex2in(nkn,nl,inodlv,ippnv,ipaux,berr)

	end

c*****************************************************************
 
        subroutine chex2in(nkn,n,nodes,ipnv,index,berr)
 
c changing extern with intern node numbers node list
 
        implicit none
 
        integer nkn,n
        logical berr
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
              berr=.true.
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
	integer iscanf

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
        ianz_grd=iscanf(line_grd,f_grd,81)
	if( ianz_grd .gt. 80 ) goto 99

	grd_next_line = .true.

	return
   99	continue
	write(ner,*) 'dimension of array f too small ',80,ianz_grd
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

	subroutine gr3_write_grid(
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

	integer nout,nout1
	integer i
	integer n,ib,nmax
	integer k,ie,il
	integer ia,ieext,ilext
	real depth
	logical bsort,bextern
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)

	integer ifileo

	bsort = .false.
	bsort = .true.		!sort nodes on external numbering

	bextern = .true.	!use external nodes
	bextern = .false.

	!if( .not. bextern ) bsort = .false.

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99
	nout1 = ifileo(1,'bas_bnd.gr3','formatted','unknown')
	if( nout1 .le. 0 ) goto 99

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

	write(nout,*)file
	write(nout,*)ne,nk
	write(nout1,*)file
	write(nout1,*)ne,nk

	call isort(nk,ippnv,ipdex)

	do k=1,nk
	  depth = hhnv(k)
	  write(nout,1000) nextern(k),xv(k),yv(k),depth
	  depth = 0.
	  write(nout1,1000) nextern(k),xv(k),yv(k),depth
	end do

	!write(nout,*)

	if( bsort ) call isort(ne,ippev,ipdex)

	do ie=1,ne
	  ieext = ie
	  if( bextern ) ieext = ippev(ie)
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  write(nout,2000) ieext,n,(nextern(inodev(ib+k)),k=1,n)
	  write(nout1,2000) ieext,n,(nextern(inodev(ib+k)),k=1,n)
	end do

	close(nout)

	deallocate(ipdex)
	deallocate(nextern)

	return
 1000	format(i10,3f16.8)
 2000	format(5i10)
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop gr3_write_grid: cannot open file'
	end

c*****************************************************************

	subroutine msh_write_grid(
     +			 file
     +			,nco,nk,ne,nl,nne,nnl
     +			,ippnv,ippev,ipplv
     +			,ianv,iaev,ialv
     +			,hhnv,hhev,hhlv
     +			,xv,yv
     +                  ,ipntev,inodev
     +                  ,ipntlv,inodlv
     +			)

c writes grd file in msh format (gmsh 2.2)

        use mod_geom_dynamic

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

	integer nb		!number of boundary nodes
	integer nit		!number of items (boundary nodes + elem)
	integer nout
	integer i,idit
	integer n,ib,nmax
	integer k,ie,il
	integer ia,flag,ilext
	real depth
	logical bsort,bextern,bbound
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)

	integer ifileo

	bbound = .false.
	bbound = .true.		!include boundary nodes
	bextern = .true.	!use external nodes
	bextern = .false.

	!if( .not. bextern ) bsort = .false.

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

	! write header
	write(nout,'((a14))')'$MeshFormat   '
	write(nout,'((a14))')'2.2 0 8       '
	write(nout,'((a14))')'$EndMeshFormat'
	write(nout,'((a14))')'$Nodes        '
	write(nout,*)nk

	! write nodes
	nb = 0
	do k=1,nk
	  depth = hhnv(k)
	  write(nout,1000) nextern(k),xv(k),yv(k),depth
	  if (inodv(k) < 0 ) nb = nb + 1
	end do
	write(nout,'((a14))')'$EndNodes     '

	! write items (boundary nodes + elements)
	nit = ne
	if ( bbound ) nit = nit + nb

	write(nout,'((a14))')'$Elements     '
	write(nout,*)nit

	! write boundary nodes
        idit = 0
        flag = 15		!for boundary node
	if ( bbound ) then
  	  do k=1,nk
	    if (inodv(k) < 0 ) then
              idit = idit + 1
              ia = ianv(k)
	      write(nout,2000) idit,flag,1,ia,k
            end if
	  end do
        end if

	! write elements
        flag = 2		!for 3-node triangle
	do ie=1,ne
          idit = idit + 1
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  write(nout,2000) idit,flag,1,ia,(nextern(inodev(ib+k)),k=1,n)
	end do

	write(nout,'((a12))')'$EndElements  '
	close(nout)

	deallocate(ipdex)
	deallocate(nextern)

	return
 1000	format(i10,3f16.8)
 2000	format(7i10)
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop msh_write_grid: cannot open file'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_get_total_lines(nl)

	use grd

	implicit none

	integer nl

	nl = nl_grd

	end

c*****************************************************************

	subroutine grd_get_line_params(il,inum,itype,nvert,depth)

	use grd

	implicit none

	integer il
	integer inum
	integer itype
	integer nvert
	real depth

        inum = ipplv(il)
        itype = ialv(il)
	depth = hhlv(il)
	nvert = ipntlv(il) - ipntlv(il-1)

	end

c*****************************************************************

	subroutine grd_get_line_array(il,nvert,nodes,x,y)

	use grd

	implicit none

	integer il
	integer nvert
	integer nodes(nvert)
	real x(nvert)
	real y(nvert)

	integer ibase,i,node

	ibase = ipntlv(il-1)
	nvert = ipntlv(il) - ipntlv(il-1)
	nodes(1:nvert) = inodlv(ibase+1:ibase+nvert)

	do i=1,nvert
	  node = nodes(i)
	  x(i) = xv(node)
	  y(i) = yv(node)
	end do

	end

c*****************************************************************

	subroutine grd_get_total_nodes(nk)

	use grd

	implicit none

	integer nk

	nk = nk_grd

	end

c*****************************************************************

	subroutine grd_get_node_params(ik,inum,itype,x,y,depth)

	use grd

	implicit none

	integer ik
	integer inum
	integer itype
	real x,y
	real depth

        inum = ippnv(ik)
        itype = ianv(ik)
	x = xv(ik)
	y = yv(ik)
	depth = hhnv(ik)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine grd_get_node_arrays(n,ipv,iav,hv,x,y)

	use grd

	implicit none

	integer n
	integer ipv(n),iav(n)
	real hv(n),x(n),y(n)

	if( n /= nk_grd ) stop 'error stop grd_get_node_arrays: n'

	ipv = ippnv
	iav = ianv
	hv = hhnv
	x = xv
	y = yv

	end 

c*****************************************************************

	subroutine grd_get_elem_arrays(n,ipv,iav,hv)

	use grd

	implicit none

	integer n
	integer ipv(n),iav(n)
	real hv(n)

	if( n /= ne_grd ) stop 'error stop grd_get_node_arrays: n'

	ipv = ippev
	iav = iaev
	hv = hhev

	end 

c*****************************************************************
c*****************************************************************
c*****************************************************************
c line routines
c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c all following routines read lines from GRD, BND or XY format
c they return n,x,y,ifl
c if n == 0 on entry only counting of x/y points
c if n == 0 on return an error has occured
c ifl == 1 for start of line and 0 for following points
c XY format can read only one line
c
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine read_all_lines(file,n,x,y,ifl)

c reads file with lines and general format (grd,bnd,xy)

	implicit none

	character*(*) file
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	logical is_grd_file,is_bnd_file,is_xy_file

	if( is_grd_file(file) ) then
	  call read_grd_lines(file,n,x,y,ifl)
	else if( is_bnd_file(file) ) then
	  call read_bnd_lines(file,n,x,y,ifl)
	else if( is_xy_file(file) ) then
	  call read_xy_lines(file,n,x,y,ifl)
	else
	  write(6,*) 'cannot determine file format: ',trim(file)
	  stop 'error stop read_all_lines: file format'
	end if

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine read_grd_lines(grdfile,n,x,y,ifl)

c reads grd file for lines to be treated (plot or particle release)

	implicit none

	character*(*) grdfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer nl,il,inum,itype,nvert,ndim,ip
	real depth

	logical is_grd_file

	ndim = n
	n = 0

	if( .not. is_grd_file(grdfile) ) then
	  write(6,*) 'file is not a grd file: ',trim(grdfile)
	  write(6,*) 'cannot read old format of bnd files'
	  write(6,*) 'use bnd2grd to convert from bnd to grd file'
	  !stop 'error stop read_bnd_lines: no bnd file read'
	  return
	end if

	call grd_read(grdfile)

	call grd_get_total_lines(nl)

	if( ndim == 0 ) then
	  do il=1,nl
	    call grd_get_line_params(il,inum,itype,nvert,depth)
	    n = n + nvert
	    !write(6,*) 'counting line: ',nl,il,nvert
	  end do
	  !write(6,*) 'total number of points in line: ',n
	else
	  ip = 1
	  do il=1,nl
	    call grd_get_line_array(il,nvert,ifl(ip),x(ip),y(ip))
	    !write(6,*) 'reading line: ',nl,il,nvert,ip
	    ifl(ip:ip+nvert-1) = 0
	    ifl(ip) = 1
	    ip = ip + nvert
	    if( ip-1 > ndim ) goto 99
	  end do
	  n = ip - 1
	  !write(6,*) 'number of points in line read: ',n
	end if

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_grd_lines: ndim'
	end

c*****************************************************************

	subroutine read_bnd_lines(bndfile,n,x,y,ifl)

c reads old bnd file format

	implicit none

	character*(*) bndfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer ndim,ios,iflag
	real xx,yy

	ndim = n
	n = 0

	open(1,iostat=ios,file=bndfile,status='old',form='formatted')
	if( ios /= 0 ) return

	do
	  read(1,*,iostat=ios) xx,yy,iflag
	  if( ios /= 0 ) exit
	  n = n + 1
	  if( ndim > 0 ) then
	    if( n > ndim ) goto 99
	    x(n) = xx
	    y(n) = yy
	    ifl(n) = iflag
	  end if
	end do

	close(1)

	if( ios > 0 ) n = 0

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_bnd_lines: ndim'
	end

c*****************************************************************

	subroutine read_xy_lines(xyfile,n,x,y,ifl)

c reads xy file format - only one line allowed

	implicit none

	character*(*) xyfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer ndim,ios,iflag
	real xx,yy

	ndim = n
	n = 0
	iflag = 1	!only first point flagged with 1, rest 0

	open(1,iostat=ios,file=xyfile,status='old',form='formatted')
	if( ios /= 0 ) return

	do
	  read(1,*,iostat=ios) xx,yy
	  if( ios /= 0 ) exit
	  n = n + 1
	  if( ndim > 0 ) then
	    if( n > ndim ) goto 99
	    x(n) = xx
	    y(n) = yy
	    ifl(n) = iflag
	  end if
	  iflag = 0
	end do

	close(1)

	if( ios > 0 ) n = 0

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_xy_lines: ndim'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	function is_bnd_file(file)

	implicit none

	logical is_bnd_file
	character*(*) file

	character*80 line
	integer ios,ianz
	real f(4)
	integer iscanf

	is_bnd_file = .false.

	open(1,iostat=ios,file=file,status='old',form='formatted')
	if( ios /= 0 ) return
	read(1,*,iostat=ios) line
	close(1)

	if( ios /= 0 ) return
	ianz = iscanf(line,f,4)
	if( ianz /= 3 ) return
	if( nint(f(3)) /= 0 .and. nint(f(3)) /= 1 ) return

	is_bnd_file = .true.

	end

c*****************************************************************

	function is_xy_file(file)

	implicit none

	logical is_xy_file
	character*(*) file

	character*80 line
	integer ios,ianz
	real f(4)
	integer iscanf

	is_xy_file = .false.

	open(1,iostat=ios,file=file,status='old',form='formatted')
	if( ios /= 0 ) return
	read(1,*,iostat=ios) line
	close(1)

	if( ios /= 0 ) return
	ianz = iscanf(line,f,4)
	if( ianz /= 2 ) return

	is_xy_file = .true.

	end

c*****************************************************************
