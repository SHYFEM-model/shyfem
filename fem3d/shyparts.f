
!--------------------------------------------------------------------------
!
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

c revision log :
c
c 19.05.2020	ccf	started from scratch
c 24.05.2020	ggu	debug option added
c 28.05.2020	ggu	some more informational messages
c 15.07.2020	ggu	bug fix for counting elements
c 22.04.2021	ggu	resolve bound check error (not yet finished)
c
c****************************************************************

        program shyparts

c partitions grd file using the METIS library
c requires metis-5.1.0.tar.gz

	use mod_geom
	use evgeom
	use basin
	use clo
	use basutil
	use shympi
        use grd

	implicit none

	integer		      :: k,ie,ii,l,ic
        integer               :: nparts		!parts to partition the mesh
	integer, allocatable  :: epart(:)	!partition vector for elements
	integer, allocatable  :: npart(:)	!partition vector for nodes

	logical bdebug
	integer ierr1,ierr2
	character*80 grdfile

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        call shyfem_copyright('shyparts - partitioning a SHYFEM grid')

        call shyparts_init(grdfile,nparts,bdebug)

        if( grdfile == ' ' ) call clo_usage
	call grd_set_write(.false.)
        call read_command_line_file(grdfile)

        call shympi_init(.false.)

c-----------------------------------------------------------------
c initialiaze arrays
c-----------------------------------------------------------------

        allocate(epart(nel))
        allocate(npart(nkn))

c-----------------------------------------------------------------
c do the partitioning
c-----------------------------------------------------------------

	call do_partition(nkn,nel,nen3v,nparts,npart,epart)

c-----------------------------------------------------------------
c check partition
c-----------------------------------------------------------------

	call check_partition(npart,epart,ierr1,ierr2)

!-----------------------------------------------------------------
! write partition information to terminal
!-----------------------------------------------------------------

	call info_partition(nparts)

c-----------------------------------------------------------------
c write grd files
c-----------------------------------------------------------------

        call write_partition_to_grd(grdfile,bdebug
     +                  ,nparts
     +                  ,npart,epart
     +			,ierr1,ierr2)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine read_command_line_file(file)

	use basin
	use basutil

	implicit none

	character*(*) file
	logical is_grd_file

	if( basin_is_basin(file) ) then
	  write(6,*) 'reading BAS file ',trim(file)
	  call basin_read(file)
	  breadbas = .true.
	else if( is_grd_file(file) ) then
	  write(6,*) 'reading GRD file ',trim(file)
	  call grd_read(file)
	  call grd_to_basin
	  call estimate_ngr(ngr)
	  breadbas = .false.
	else
	  write(6,*) 'Cannot read this file: ',trim(file)
	  stop 'error stop read_given_file: format not recognized'
	end if

	end

c*******************************************************************

	subroutine shyparts_init(grdfile,nparts,bdebug)

	use clo

	implicit none

	character*(*) grdfile
        integer nparts
	logical bdebug

	call clo_init('shyparts','grd-file','3.0')

        call clo_add_info('partitioning of grd file with METIS')

	call clo_add_sep('options for partitioning')
        call clo_add_option('nparts',-1,'number of partitions')
        call clo_add_option('debug',.false.,'write debug grd files')

	call clo_parse_options

	call clo_check_files(1)
	call clo_get_file(1,grdfile)

        call clo_get_option('nparts',nparts)
        call clo_get_option('debug',bdebug)

        if (nparts < 2 ) then
          write(6,*) 'nparts: ',nparts
	  stop 'error stop shyparts_init: nparts < 2'
        end if

	end

c*******************************************************************

        subroutine node_test
        end

c*******************************************************************

