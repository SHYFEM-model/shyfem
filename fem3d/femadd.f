
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2017  Marco Bajo
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

! elaborates fem files
!
! revision log :
!
! 14.01.2015	ggu	adapted from feminf
! 20.05.2015	ggu	use bhuman to convert to human readable time
! 05.06.2015	ggu	iextract to extract nodal value
! 05.11.2015	ggu	new option chform to change format
! 04.10.2016	ggu	output flags now similar to shyelab
! 05.10.2016	ggu	allow for expansion of regular grid
! 11.10.2016	ggu	introduced flag for min/max/med computation
! 31.10.2016	ggu	new flag condense (bcondense)
! 16.05.2017	ggu&mbj	better handling of points to extract
! 30.01.2018	ggu	written with new fem_util module
! 22.02.2018	ggu	changed VERS_7_5_42
! 16.02.2019	ggu	changed VERS_7_5_60
! 27.01.2022	ggu	minor changes
! 16.03.2022	ggu	femadd newly written
!
!******************************************************************

	program femadd

! adds values of different files to one variable

	use clo
	use fem_util

	implicit none

	character*80 name,string,infile
	integer nfile,i,ierr,iformat
	integer nvar,nvar0,nrecs
	double precision atime,atime0
	logical bdebug
	logical bverb,bquiet,bsilent
	logical bunform
	character*20 aline
	type(femfile_type), allocatable :: ffinfo(:)
	type(femfile_type) :: ffiout
	type(femrec_type), allocatable :: finfo(:)
	type(femrec_type) :: fout

	bdebug = .true.
	bdebug = .false.

!--------------------------------------------------------------
! set command line options
!--------------------------------------------------------------

	call clo_init('femadd','fem-files','1.0')

	call clo_add_info('adds vars of multiple fem-files into one')

        call clo_add_sep('options in/output')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'do not write time records')

	call clo_add_sep('additional options')

	call clo_add_option('unform',.false.
     +				,'write output file unformatted')

!--------------------------------------------------------------
! parse command line options
!--------------------------------------------------------------

	call clo_parse_options(1)  !expecting (at least) 1 file after options

!--------------------------------------------------------------
! get command line options
!--------------------------------------------------------------

	call clo_get_option('verb',bverb)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('silent',bsilent)

	call clo_get_option('unform',bunform)

	if( bsilent ) bquiet = .true.

!--------------------------------------------------------------
! set parameters
!--------------------------------------------------------------

	nfile = clo_number_of_files()

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bunform
	end if

!--------------------------------------------------------------
! open all files
!--------------------------------------------------------------

	if( nfile < 1 ) then
	  write(6,*) 'No file given... exiting'
	  stop 'error stop femadd: no files'
	else if( nfile == 1 ) then
	  write(6,*) 'Only one file given... exiting'
	  stop 'error stop femadd: nothing to add'
	end if

	allocate(ffinfo(nfile))
	allocate(finfo(nfile))

	do i=1,nfile
	  call femutil_init_record(finfo(i))
          call clo_get_file(i,infile)
	  call femutil_open_for_read(infile,0,ffinfo(i),ierr)
	  if( ierr /= 0 ) goto 99
	end do

	iformat = 1
	if( bunform ) iformat = 0
	call femutil_open_for_write('out.fem',iformat,ffiout)

!--------------------------------------------------------------
! loop on files and read data
!--------------------------------------------------------------

	nrecs = 0
	nvar0 = 0

	do

	nvar = 0
	do i=1,nfile
	  call femutil_read_record(ffinfo(i),finfo(i),ierr)
	  !if( i == 1 .and. ierr < 0 ) exit
	  if( ierr < 0 ) exit
	  if( ierr /= 0 ) goto 98
	  call femutil_get_time(finfo(i),atime)
	  if( i == 1 ) atime0 = atime
	  if( atime /= atime0 ) goto 97
	  if( .not. femutil_is_compatible(finfo(1),finfo(i)) ) goto 96
	  nvar = finfo(i)%nvar
	  if( nvar0 == 0 ) nvar0 = nvar
	  if( nvar /= nvar0 ) goto 95
	end do

	if( ierr < 0 ) exit
	nrecs = nrecs + 1

	!------------------------------------------------------
	! add data and write output file
	!------------------------------------------------------

        call dts_format_abs_time(atime,aline)
	if( .not. bquiet ) write(6,*) atime,'  ',aline
	call femutil_add_data_recs(nfile,finfo,fout)
	call femutil_write_record(ffiout,fout)

	end do

	if( .not. bsilent ) then
	  write(6,*) 'total number of records treated: ',nrecs
	  write(6,*) 'output written to file out.fem'
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   95	continue
	write(6,*) 'nvar is changing: ',nvar,nvar0
	stop 'error stop femadd: nvar not constant'
   96	continue
	write(6,*) 'files are not compatible'
	stop 'error stop femadd: not compatible'
   97	continue
	write(6,*) 'times are not compatible: ',atime0,atime
	stop 'error stop femadd: time error'
   98	continue
	write(6,*) 'error reading record ',i,ierr
	stop 'error stop femadd: read error'
   99	continue
	write(6,*) 'error opening file ',infile
	stop 'error stop femadd: opening error'
        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

