
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

! transfroms restart file to shy format

! revision log :
!
! 03.05.2019    ggu     written from rstinf.f

!******************************************************************

	program rst2shy

	use basin
	use levels
	use shympi
	use mod_conz
	use mod_ts
	use mod_hydro
	use mod_geom_dynamic

	implicit none

	logical bhydro,bts,bconz
	integer iunit,it,nvers,nrec,nknr,nelr,nlvr,iflag,ierr
	integer ic
	integer nread
	integer id_ts,id_conz,id_hydro
	double precision atime,dtime
	double precision atime_anf
	double precision atime_end
	character*20 aline
	character*80 rstfile,basfile
	character*80 title1
	character*80 title2
	character*80 ext,name

	logical has_flag

!-------------------------------------------------------------------
! create title strings
!-------------------------------------------------------------------

	title1 = 'version nrec       nkn       nel       nlv' //
     +			'     iconz     iflag'
!                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec                         atime     date'

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1

        call rst_init(rstfile,basfile)

	call basin_read(basfile)
	call shympi_init(.false.)               !call after basin has been read


!-------------------------------------------------------------------
! read first record and params
!-------------------------------------------------------------------

	open(iunit,file=rstfile,status='old',form='unformatted')
	call rst_skip_record(iunit,atime,nvers,nrec
     +					,nknr,nelr,nlvr,ic,iflag,ierr)
	if( ierr .ne. 0 ) then
	  stop 'error stop rst2shy: reading parameters'
	end if
	close(iunit)

	bhydro = has_flag(iflag,1)
	bts = has_flag(iflag,100)
	bconz = has_flag(iflag,1000)
	write(6,*) 'bhydro,bts,bconz: ',bhydro,bts,bconz

	write(6,1000) trim(title1)
	write(6,1010) nvers,nrec,nknr,nelr,nlvr,ic,iflag
	write(6,*)
	write(6,1001) trim(title2)
	atime_anf = atime

	call levels_init(nknr,nelr,nlvr)
	nlv = nlvr
	write(6,*) 'nlv = ',nlv

	call check_name_and_extension(rstfile,name,ext)
	if( bts ) then
	  call mod_ts_init(nkn,nlv)
	  call ts_init_output(id_ts,atime,nlvr,name)
	end if
	if( bconz ) then
	  call mod_conz_init(ic,nkn,nlv)
	  call conz_init_output(id_conz,atime,nlvr,ic,name)
	end if
	if( bhydro) then
	  call mod_hydro_init(nkn,nel,nlv)
	  call mod_geom_dynamic_init(nkn,nel)
	  call hydro_init_output(id_hydro,atime,nlvr,name)
	end if

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	open(iunit,file=rstfile,status='old',form='unformatted')

	do

	  call rst_read_record(iunit,atime,ierr)
	  if( ierr .ne. 0 ) exit
	  nread = nread + 1
	  call dts_format_abs_time(atime,aline)
	  write(6,1011) nread,atime,aline
	  atime_end = atime
	  dtime = atime - atime_anf

	  if( bts ) call ts_write_output(id_ts,dtime)
	  if( bconz ) call conz_write_output(id_conz,dtime,ic)
	  if( bhydro ) call hydro_write_output(id_hydro,dtime)
	end do

	if( ierr > 0 ) stop 'error stop rstinf: error reading record'

!-------------------------------------------------------------------
! final message
!-------------------------------------------------------------------

	write(6,1001) trim(title2)
	write(6,*)
	write(6,1000) trim(title1)
	write(6,1010) nvers,nrec,nknr,nelr,nlvr,ic,iflag
	write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)
	write(6,*) 'Meaning of iflag:'
	write(6,*) '         1          hydro'
	write(6,*) '        10          depth'
	write(6,*) '       100          ibarcl (T/S/rho)'
	write(6,*) '      1000          iconz (cnv/conzv)'
	write(6,*) '     10000          iwvert (wlnv)'
	write(6,*) '    100000          ieco (ecological variables)'

	!call check_flag(iflag)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	stop
 1000	format(a)
 1001	format(a)
 1010	format(i7,i5,5i10)
 1011	format(i7,f30.2,5x,a20)
	end

!******************************************************************

	subroutine check_flag(iflag)

	implicit none

	integer iflag

	integer i,is
	logical bflag

	logical has_flag

	is = 1

	do i=1,5
	  bflag = has_flag(iflag,is)
	  write(6,*) iflag,i,is,bflag
	  is = is * 10
	end do

	end

!******************************************************************

	function has_flag(iflag,is)

	implicit none

	logical has_flag
	integer iflag,is

	integer ifl,ifl1,ifl2

	ifl = iflag/is
	ifl1 = (is*10) * (iflag/(is*10))
	
	ifl2 = (ifl*is) - ifl1

	has_flag = (ifl2 /= 0)

	end

!******************************************************************

        subroutine rst_init(rstfile,basfile)

        use clo
        use basin
	!use mod_restart

        implicit none

        character*(*) rstfile,basfile

	integer i,ifiles
	character*80 ext,file
	logical rst_is_rst_file

        call shyfem_copyright('rst2shy - transforms rst file to shy')

        call clo_init('rst2shy','rstfile basfile','1.0')

        call clo_add_info('transfroms rst to shy files')

        call clo_parse_options

	ifiles = clo_number_of_files()
	if( ifiles /= 2 ) then
	  write(6,*) 'need two files, one rst and one bas file'
	  stop 'error stop rst_init: missing files'
	end if
        call clo_check_files(2)

	rstfile = ' '
	basfile = ' '

	do i=1,2
          call clo_get_file(i,file)
	  call check_extension(file,ext)
	  if( ext == 'bas' ) then
	    basfile = file
	    if( .not. basin_is_basin(file) ) then
	      write(6,*) 'file is not a bas file: ',trim(file)
	      goto 99
	    end if
	  else if( ext == 'rst' ) then
	    rstfile = file
	    if( .not. rst_is_rst_file(file) ) then
	      write(6,*) 'file is not a rst file: ',trim(file)
	      goto 99
	    end if
	  else
	    write(6,*) 'unknown extension of file: ',trim(file)
	    goto 99
	  end if
	end do

	if( basfile == ' ' ) then
	  write(6,*) 'no bas file read...'
	  goto 99
	end if

	if( rstfile == ' ' ) then
	  write(6,*) 'no rst file read...'
	  goto 99
	end if

	return
   99	continue
	stop 'error stop rst_init: error reading file'
        end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine ts_init_output(id,atime,nl,name)

	use shyfile
	use iso8601

	implicit none

	integer id
	double precision atime
	integer nl
	character*(*) name

	integer nvar,npr,ftype,ierr
	integer date,time
	character*80 ext,femver,title,file,aline

	ext = '.ts.shy'
	file = trim(name) // ext
	femver = 'unknown'
	title = 'unknown'

	npr = 1
	ftype = 2
	nvar = 3

	call dts_format_abs_time(atime,aline)
	call string2date_and_time(aline,date,time,ierr)
	if( ierr /= 0 ) stop 'error stop ts_init_output: ierr'

        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)
        call shy_make_header(id)

	end

!******************************************************************

	subroutine ts_write_output(id,dtime)

	use levels
	use mod_ts

	implicit none

	integer id
	double precision dtime

	call shy_write_scalar_record(id,dtime,11,nlvdi,saltv)
	call shy_write_scalar_record(id,dtime,12,nlvdi,tempv)
	call shy_write_scalar_record(id,dtime,13,nlvdi,rhov)

	end

!******************************************************************

	subroutine conz_init_output(id,atime,nl,ic,name)

	use shyfile
	use iso8601

	implicit none

	integer id
	double precision atime
	integer nl
	integer ic
	character*(*) name

	integer nvar,npr,ftype,ierr
	integer date,time
	character*80 ext,femver,title,file,aline

	ext = '.conz.shy'
	file = trim(name) // ext
	femver = 'unknown'
	title = 'unknown'

	npr = 1
	ftype = 2
	nvar = ic

	call dts_format_abs_time(atime,aline)
	call string2date_and_time(aline,date,time,ierr)
	if( ierr /= 0 ) stop 'error stop conz_init_output: ierr'

        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)
        call shy_make_header(id)

	end

!******************************************************************

	subroutine conz_write_output(id,dtime,ic)

	use levels
	use mod_conz

	implicit none

	integer id
	double precision dtime
	integer ic

	integer i,idc

	if( ic == 1 ) then
	  idc = 10
	  call shy_write_scalar_record(id,dtime,idc,nlvdi,cnv)
	else if( ic > 1 ) then
	  do i=1,ic
	    idc = 300 + 1
	    call shy_write_scalar_record(id,dtime,idc,nlvdi
     +				,conzv(1,1,i))
	  end do
	end if

	end

!******************************************************************

	subroutine hydro_init_output(id,atime,nl,name)

	use shyfile
	use iso8601

	implicit none

	integer id
	double precision atime
	integer nl
	character*(*) name

	integer nvar,npr,ftype,ierr
	integer date,time
	character*80 ext,femver,title,file,aline

	ext = '.hydro.shy'
	file = trim(name) // ext
	femver = 'unknown'
	title = 'unknown'

	npr = 3
	ftype = 1
	nvar = 4

	call dts_format_abs_time(atime,aline)
	call string2date_and_time(aline,date,time,ierr)
	if( ierr /= 0 ) stop 'error stop hydro_init_output: ierr'

        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)
        call shy_make_header(id)

	end

!******************************************************************

	subroutine hydro_write_output(id,dtime)

	use levels
	use mod_hydro

	implicit none

	integer id
	double precision dtime

	call shy_write_hydro_records(id,dtime,nlvdi
     +			,znv,zenv,utlnv,vtlnv)

	end

!******************************************************************

