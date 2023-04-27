
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

! revision log :
!
! 19.05.2020	ccf	started from scratch
! 24.05.2020	ggu	debug option added
! 28.05.2020	ggu	some more informational messages
! 15.07.2020	ggu	bug fix for counting elements
! 22.04.2021	ggu	resolve bound check error (not yet finished)
! 12.04.2022	ggu	restructured to allow for online computing
! 04.04.2023	ggu	minor changes
! 27.04.2023	ggu	better error handling
!
!****************************************************************

        program shyparts

! partitions grd file using the METIS library
! requires metis-5.1.0.tar.gz

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

!-----------------------------------------------------------------
! read in basin
!-----------------------------------------------------------------

        call shyfem_copyright('shyparts - partitioning a SHYFEM grid')

        call shyparts_init(grdfile,nparts,bdebug)

        if( grdfile == ' ' ) call clo_usage
	call grd_set_write(.false.)
        call read_command_line_file(grdfile)

        call shympi_init(.false.)

!-----------------------------------------------------------------
! initialiaze arrays
!-----------------------------------------------------------------

        allocate(epart(nel))
        allocate(npart(nkn))
	npart = 0
	epart = 0

!-----------------------------------------------------------------
! do the partitioning
!-----------------------------------------------------------------

	call do_partition(nkn,nel,nen3v,nparts,npart,epart)

!-----------------------------------------------------------------
! check partition
!-----------------------------------------------------------------

	call check_partition(npart,epart,ierr1,ierr2)

!-----------------------------------------------------------------
! write partition information to terminal
!-----------------------------------------------------------------

	call info_partition(nparts,npart)

!-----------------------------------------------------------------
! write grd files
!-----------------------------------------------------------------

        call write_partition_to_grd(grdfile,bdebug
     +                  ,nparts,npart,epart)

!-----------------------------------------------------------------
! final message
!-----------------------------------------------------------------

        if( ierr1 == 0 .and. ierr2 == 0 ) then
          write(6,*) 'domain successfully partitioned'
          call exit(0)
        else
          write(6,*) 'errors in domains: ',ierr1,ierr2
          call exit(9)
        end if

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

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

!*******************************************************************

	subroutine shyparts_init(grdfile,np,bdebug)

	use clo

	implicit none

	character*(*) grdfile
        integer np
	logical bdebug

	call clo_init('shyparts','grd-file','3.0')

        call clo_add_info('partitioning of grd file with METIS')

	call clo_add_sep('options for partitioning')
        call clo_add_option('np',-1,'number of partitions')
        call clo_add_option('debug',.false.,'write debug grd files')

	call clo_parse_options

	call clo_check_files(1)
	call clo_get_file(1,grdfile)

        call clo_get_option('np',np)
        call clo_get_option('debug',bdebug)

        if (np == -1 ) then
	  stop 'error stop shyparts_init: no np given'
	else if( np == 1 ) then
	  stop 'error stop shyparts_init: no partitioning for np == 1'
	end if

        if (np < 2 ) then
          write(6,*) 'nparts: ',np
	  stop 'error stop shyparts_init: np < 2'
        end if

	end

!*******************************************************************

        subroutine node_test
        end

!*******************************************************************

