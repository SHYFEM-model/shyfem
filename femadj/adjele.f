c
c $Id: adjele.f,v 1.21 2010-03-11 15:31:28 georg Exp $
c
c description :
c
c main for mesh adjustment
c
c contents :
c
c adjele			main
c
c revision log :
c
c 19.05.2003    ggu     some more checks
c 02.12.2004    ggu     some more debug checks (version 1.41)
c 02.12.2004    ggu     copyright statement added
c 20.03.2007    ggu     fake version
c 16.04.2008    ggu     new Makefile structure
c 09.12.2008    ggu     small changes to Makefile
c 26.01.2009    ggu     makedepend, compiler warnings
c 06.04.2009    ggu     new param.h, include basin.h, better error handling
c 24.04.2009    ggu     new call to rdgrd()
c 30.03.2012    ggu     new call to node_info() (was bug)
c
c notes :
c
c to use static (non moveable) nodes please consult nbstatic.h
c
c***************************************************************

	program adjele

c adjusts elements after automatic mesh generator

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

c------------------------------------------------------- parameters
	include 'param.h'
c------------------------------------------------------- grade index
c------------------------------------------------------- basin
c------------------------------------------------------- local
	integer kspecial
	integer nlidim,nlndim
	integer k
	integer ner,nc
	integer nco,nknh,nli
	integer nk,ne,nl,nne,nnl
	logical bstop, bplot
	character*80 file

	real, allocatable :: dx(:),dy(:)
	integer, allocatable :: ic(:)
	integer, allocatable :: ianv(:)
	integer, allocatable :: ipaux(:)
	integer, allocatable :: iaux(:)
	real, allocatable :: raux(:)
c------------------------------------------------------- end declaration

	kspecial = 0		!set to 0 for no debug

	file = 'old.grd'
	file = ' '
	ner = 6
	bstop = .false.
	bplot = .false.		!for plotting of basin and grades
	bplot = .true.		!for plotting of basin and grades

c-------------------------------------------------------------------

	call shyfem_copyright('adjele - regularize finite element grids')

c-------------------------------------------------------------------

c get file name from command line

        nc = command_argument_count()
        if( nc .ne. 1 ) then
          write(6,*) 'Usage: adjele grd-file'
          stop 'error stop adjele: no file given'
        end if

        call get_command_argument(1,file)

c read grid file with nodes and elements

	call grd_read(file)
	call grd_get_params(nk,ne,nl,nne,nnl)

	write(6,'(a,5i7)') 'grid parameters: ',nk,ne,nl
	if( nk .le. 0 ) stop

	call grd_to_basin

	allocate(dx(nkn),dy(nkn))
	allocate(ic(nkn),ianv(nkn))
	allocate(ipaux(nkn),iaux(nel))
	allocate(raux(nel))

c save depth information in elements to nodes

	call mod_depth_init(nkn,nel)
	call grd_get_depth(nk,ne,hkv,hev)

	nknh = 0
	do k=1,nkn
	  if( hkv(k) .ne. -999. ) nknh = nknh + 1
	end do

	if( nknh .eq. 0 ) then
	  write(6,*) 'copying element depth to nodes...'
	  call hev2hkv
	end if

c determine grade

	call maxgrd(nkn,nel,nen3v,ngr)

	write(6,*) 'maximum grade: ',ngr
	call mod_adj_grade_init(nkn,ngr)

	call setgrd(nkn,nel,nen3v,ngrade)

	call stats('first call')

c make boundary nodes (flag nbound)

	call mkbound(nkn,nel,ngrdi,nen3v,ngrade,nbound,ngri)
        call mkstatic(nkn,ianv,nbound)
	call stats('boundary nodes')

	call node_info(kspecial)

c plot grade

	call qopen

	if( bplot) call plobas

c eliminate 4- grades

        write(6,*) '================================='
        write(6,*) 'first cycle...'
        write(6,*) '================================='

        write(6,*) 'start eliminating nodes ...'

	call chkgrd(' ')
	call elimlow
	if( bplot) call plobas
	call stats('4- grades')
	call node_info(kspecial)

c eliminate 8+ grades

	call chkgrd(' ')
	call elimhigh(8)
	if( bplot) call plobas
	call stats('8+ grades')

	call chkgrd('checking after 8+')
	call node_info(kspecial)

c eliminate 7+ grades

	call chkgrd(' ')
	call elimhigh(7)
	if( bplot) call plobas
	call stats('7+ grades')

	call chkgrd('checking after 7+')
	call node_info(kspecial)

c smoothing

	!call write_grid('new_nosmooth.grd')

        call smooth_grid(50,0.1)
	if( bplot) call plobas
	call node_info(kspecial)

	!call write_grid('new_smooth1.grd')

c again ...

        write(6,*) '================================='
        write(6,*) 'second cycle...'
        write(6,*) '================================='

	call chkgrd('start of second cycle')
        call elimlow
	call elimhigh(8)
	call elimhigh(7)
	if( bplot) call plobas
	call stats('one more time')
	call node_info(kspecial)

c eliminate 5 grades

	!call write_grid('new_help.grd')

	call chkgrd(' ')
	call elim5
	if( bplot) call plobas
	call stats('5 grades')

	call chkgrd('checking after 5')
	call node_info(kspecial)

c eliminate 5-5 grades

	call write_grid('new_help1.grd')

	call elim57
	!call write_grid('new_help2.grd')
	if( bplot) call plobas
	call stats('5-5 grades')
	call chkgrd('checking after 5-5')
	call node_info(kspecial)

c one more time

        write(6,*) '================================='
        write(6,*) 'third cycle...'
        write(6,*) '================================='

	!call write_grid('cycle3_start.grd')

	call elimhigh(8)
        call elimhigh(7)
        call elim5
        call elim57
        call stats('all again')
	call chkgrd(' ')
	call node_info(kspecial)

c smoothing

        write(6,*) '================================='
        write(6,*) 'final smoothing...'
        write(6,*) '================================='

        call smooth_grid(50,0.1)
	if( bplot) call plobas
	call node_info(kspecial)

c write to grd file

        write(6,*) '================================='
        write(6,*) 'writing to grid...'
        write(6,*) '================================='

	call chkgrd(' ')
	call node_info(kspecial)
	call write_grid('new.grd')

	call qclose

	stop
   97	continue
	write(6,*) 'error reading grd file'
	stop 'error stop rdgrd'
   99	continue
	write(6,*) 'error external to internal numbering'
	stop 'error stop ex2in'
	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine hev2hkv

c saves information about depth to nodes

	use mod_depth
	use basin

	implicit none

	include 'param.h'

	integer k,ie,ii
	integer ic(nkn)

	do k=1,nkn
	  hkv(k) = 0.
	  ic(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + hev(ie)
	    ic(k) = ic(k) + 1
	  end do
	end do

	do k=1,nkn
	  if( ic(k) .gt. 0 ) then
	    hkv(k) = hkv(k) / ic(k)
	  end if
	end do

	end

c***********************************************************

	subroutine stats(text)

	use mod_adj_grade
	use basin

	implicit none

	include 'param.h'

	character*(*) text

	call statgrd(text,nkn,ngr,ngrade,nbound)

	end

c***********************************************************

