
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2007-2009,2012,2015,2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

c description :
c
c main for mesh adjustment
c
c contents :
c
c shyadj			main
c
c revision log :
c
c 01.01.2003	ggu	written
c 19.05.2003	ggu	some more checks
c 02.12.2004	ggu	some more debug checks (version 1.41)
c 02.12.2004	ggu	copyright statement added
c 20.03.2007	ggu	fake version
c 16.04.2008	ggu	new Makefile structure
c 09.12.2008	ggu	small changes to Makefile
c 26.01.2009	ggu	makedepend, compiler warnings
c 06.04.2009	ggu	new param.h, include basin.h, better error handling
c 24.04.2009	ggu	new call to rdgrd()
c 30.03.2012	ggu	new call to node_info() (was bug)
c 12.10.2015	ggu	changed VERS_7_3_3
c 25.05.2017	ggu	changed VERS_7_5_28
c 16.10.2018	ggu	changed VERS_7_5_50
c 18.12.2018	ggu	changed VERS_7_5_52
c 21.05.2019	ggu	changed VERS_7_5_62
c 14.02.2022	ggu	set all depth values to flag
c
c notes :
c
c to use static (non moveable) nodes please consult ls mod_adj_static
c
c***************************************************************

	program shyadj

c adjusts elements after automatic mesh generator

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer kspecial
	integer nlidim,nlndim
	integer k
	integer ner,nc
	integer nco,nknh,nli
	integer nk,ne,nl,nne,nnl
	integer nsmooth
	real asmooth
	real, parameter :: flag = -999.
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

	file = ' '
	ner = 6
	bstop = .false.

c-------------------------------------------------------------------

	call shyfem_copyright('shyadj - regularize finite element grids')

c-------------------------------------------------------------------

c parse command line

        call shyadj_init(file,bplot,nsmooth,asmooth)

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
	  if( hkv(k) .ne. flag ) nknh = nknh + 1
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

	if( bplot ) call qopen

	if( bplot ) call plobas

        !call smooth_grid(nsmooth,asmooth)	!only for debug

c eliminate 4- grades

        write(6,*) '================================='
        write(6,*) 'first cycle...'
        write(6,*) '================================='

        write(6,*) 'start eliminating nodes ...'

	call chkgrd('checking original grid')
	call elimlow
	if( bplot ) call plobas
	call stats('first cycle - 4- grades')
	call node_info(kspecial)

c eliminate 8+ grades

	call chkgrd('checking before 8+')
	call elimhigh(8)
	if( bplot ) call plobas
	call stats('first cycle - 8+ grades')

	call chkgrd('checking after 8+')
	call node_info(kspecial)

c eliminate 7+ grades

	call chkgrd('checking before 7+')
	call elimhigh(7)
	if( bplot ) call plobas
	call stats('first cycle - 7+ grades')

	call chkgrd('checking after 7+')
	call node_info(kspecial)

c smoothing

	!call write_grid('new_nosmooth.grd')

        call smooth_grid(nsmooth,asmooth)
	if( bplot ) call plobas
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
	if( bplot ) call plobas
	call stats('second cycle - after elimhigh')
	call node_info(kspecial)

c eliminate 5 grades

	!call write_grid('new_help.grd')

	call chkgrd('checking before 5')
	call elim5
	if( bplot ) call plobas
	call stats('second cycle - 5 grades')

	call chkgrd('checking after 5')
	call node_info(kspecial)

c eliminate 5-5 grades

	call write_grid('new_help1.grd')

	call elim57
	!call write_grid('new_help2.grd')
	if( bplot ) call plobas
	call stats('second cycle - 5-5 grades')
	call chkgrd('checking after 5-5')
	call node_info(kspecial)

c one more time

        write(6,*) '================================='
        write(6,*) 'third cycle...'
        write(6,*) '================================='

	!call write_grid('cycle3_start.grd')

	call chkgrd('start of thrid cycle')
	call elimhigh(8)
        call elimhigh(7)
        call elim5
        call elim57
        call stats('thrid cycle - end')
	call chkgrd('checking after high, 5, 57')
	call node_info(kspecial)

c smoothing

	call write_grid('final_before_smoothing.grd')

        write(6,*) '================================='
        write(6,*) 'final smoothing...'
        write(6,*) '================================='

        call smooth_grid(nsmooth,asmooth)
	if( bplot ) call plobas
	call node_info(kspecial)

c write to grd file

        write(6,*) '================================='
        write(6,*) 'writing to grid...'
        write(6,*) '================================='

	call chkgrd('final check')
        call stats('final solution')
	call node_info(kspecial)
	call show_strange_grades
	hev = flag
	hkv = flag
	hm3v = flag
	call write_grid('adjust.grd')

	call qclose	!this is safe to call

	write(6,*) 'Successful completion.'//
     +			' Output has been written to adjust.grd'

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

	character*(*) text

	call statgrd(text,nkn,ngr,ngrade,nbound)

	end

c***********************************************************

	subroutine show_strange_grades

	use mod_adj_grade
	use basin

	implicit none

	integer k,kext,n

	write(6,*) 'Listing nodes with not fixable grades...'

        do k=1,nkn
          if( nbound(k) .ne. 0 ) cycle
          n = ngrade(k)
          if( n < 5 .or. n > 7) then
	    call nint2ext(k,kext)
	    write(6,*) 'node ',k,' (extern ',kext,') with grade ',n
          end if
        end do

	end

c***********************************************************

        subroutine shyadj_init(grdfile,bplot,nsmooth,asmooth)

        use clo

        implicit none

        character*(*) grdfile
        logical bplot
	integer nsmooth
	real asmooth

	integer n
	real f(2)
	character*80 line

	integer iscanf

        call clo_init('shyadj','grd-file','2.0')

        call clo_add_info('regolarize grd file')

        call clo_add_option('smooth params',' ','smoothing options')
        call clo_add_option('plot',.false.,'create plot of grades')

        call clo_add_sep('additional information')
        call clo_add_com('  params is nsmooth[,asmooth]')
	call clo_add_com('    nsmooth is number of smoothing iterations')
        call clo_add_com('    asmooth is smoothing strength')
        call clo_add_com('    defaults: 50,0.01')

        call clo_parse_options

        call clo_get_option('smooth',line)
        call clo_get_option('plot',bplot)

        call clo_check_files(1)
        call clo_get_file(1,grdfile)

	nsmooth = 50
	asmooth = 0.01
	n = iscanf(line,f,2)
	if( n < 0 .or. n > 2 ) then
	  write(6,*) 'error in smoothing parameters: ',trim(line)
	else if( n == 0 ) then
	  !use default
	else if( n == 1 ) then
	  nsmooth = nint(f(1))
	else
	  nsmooth = nint(f(1))
	  asmooth = f(2)
	end if

	write(6,*) 'using smoothing parameters: ',nsmooth,asmooth

        end

c***********************************************************

