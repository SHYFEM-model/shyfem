
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 22.02.2016	ggu	handle catmode
! 15.04.2016	ggu	handle gis files with substitution of colon
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 20.01.2017	ggu	changed VERS_7_5_22
! 13.02.2017	ggu	changed VERS_7_5_23
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	changed VERS_7_5_25
! 11.07.2017	ggu	changed VERS_7_5_30
! 14.11.2017	ggu	changed VERS_7_5_36
! 19.10.2018	ccf	handle lgr files
! 18.12.2018	ggu	changed VERS_7_5_52
! 21.05.2019	ggu	changed VERS_7_5_62
!
!************************************************************

!====================================================
	module plotutil
!====================================================

	implicit none

	logical, save, private :: binitialized = .false.
	double precision, save :: dflag = -999.

	logical, save :: b2d
!	logical, save :: bdir

	logical, save :: binfo
	logical, save :: bverb
	logical, save :: bquiet
	logical, save :: bsilent
	logical, save :: bwrite
	logical, save :: bsdebug

	logical, save :: bregall

	integer, save :: ifreq
	integer, save :: tmin
	integer, save :: tmax

        character*80, save :: infile		= ' '
        character*80, save :: stmin		= ' '
        character*80, save :: stmax		= ' '

	integer, save :: layer = 0
	integer, save :: ivar3 = 0
	integer, save :: ivnum = 0
	
	integer, save :: nfile = 0
	character*10, save, allocatable :: file_type(:)

        character*80, save :: shyfilename = ' '
        character*80, save :: femfilename = ' '
        character*80, save :: lgrfilename = ' '
        character*80, save :: basfilename = ' '
        character*80, save :: basintype = ' '

!====================================================
	contains
!====================================================

!************************************************************
!************************************************************
!************************************************************

	subroutine plotutil_init(type,what)

	use clo

	character*(*) type
	character*(*), optional :: what

	character*80 program

	program = 'shyplot'
	if( present(what) ) program = what

	call plotutil_set_options(type,program)
	call clo_parse_options
	call plotutil_get_options(type,program)

	binitialized = .true.

	end subroutine plotutil_init

!************************************************************

	subroutine plotutil_set_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	character*80 vers,version

	if( binitialized ) return

        call get_shyfem_version(vers)
        version = '3.1' // ' (SHYFEM version ' // trim(vers) // ')'

	if( type == 'SHY' ) then
          call clo_init(program,'file(s)',version)
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop plotutil_set_options: unknown type'
	end if

        call clo_add_info('plots a shy file')

        call clo_add_sep('general options')

        call clo_add_option('info',.false.,'only give info on header')
        call clo_add_option('verbose',.false.
     +                          ,'be more verbose, write time records')
        call clo_add_option('quiet',.false.
     +                          ,'do not write header information')
        call clo_add_option('silent',.false.,'do not write anything')
        call clo_add_option('write',.false.,'write min/max of records')
        call clo_add_option('debug',.false.,'write debug information')
	call clo_hide_option('debug')

        call clo_add_sep('time options')

        call clo_add_option('tmin time',' '
     +                  ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                  ,'only process up to time')
	call clo_add_option('freq n',0.,'frequency for plot')
	call clo_add_com('    time is either YYYY-MM-DD[::hh[:mm[:ss]]]')
        call clo_add_com('    or integer for relative time')

        call clo_add_sep('variable to plot')

	call clo_add_option('varid id',0,'plot variable id')
	call clo_add_option('varnum i',0,'plot i''th variable of file')
	call clo_add_option('varname name',' ','plot variable name')
	call clo_add_com('  varid,varnum,varname are mutually exclusive')

        call clo_add_sep('layer to plot (default 0)')

	call clo_add_option('2d',.true.,'plot vertical average (default)')
	call clo_add_option('layer l',0,'plot layer l')

        call clo_add_sep('additional options')

!	call clo_add_option('dir',.false.
!     +			,'for directional variable plot arrow')
	call clo_add_option('regall',.false.
     +			,'for regular fem files plot whole grid')

        call clo_add_com('  file can be the following:')
        call clo_add_com('    shy-file to plot results')
        call clo_add_com('    bas-file to plot bathymetry and grid')
        call clo_add_com('    one or more str-files for instructions')

	end subroutine plotutil_set_options

!************************************************************

	subroutine plotutil_get_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	character*80 text

	integer ivar
	logical bvarid,bvarnum,bvarname
	character*80 varname

	if( binitialized ) return

        call clo_get_option('info',binfo)
        call clo_get_option('verbose',bverb)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('write',bwrite)
        call clo_get_option('debug',bsdebug)

        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)
        call clo_get_option('freq',ifreq)

        call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

        call clo_get_option('varid',ivar3)
        call clo_get_option('varnum',ivnum)
        call clo_get_option('varname',varname)

        call clo_get_option('2d',b2d)
        call clo_get_option('layer',layer)
!        call clo_get_option('dir',bdir)
        call clo_get_option('regall',bregall)

	if( type == 'SHY' ) then
          text = 'shyplot - Plot SHY/FEM/LGR/BAS/GRD files'
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop plotutil_get_options: unknown type'
	end if

        if( .not. bsilent ) then
          call shyfem_copyright(trim(text))
	end if

	call bash_verbose(bverb)

	bvarid = ivar3 > 0
	bvarnum = ivnum > 0
	bvarname = varname /= ' '

        if( binfo ) bverb = .true.
        if( bwrite ) bverb = .true.
        if( bsdebug ) bverb = .true.
	if( bsilent ) bquiet = .true.

        if( count( (/bvarid,bvarnum,bvarname/) ) > 1 ) then
	  write(6,*) 'You can give only one of varid, varnum and varname'
          stop 'error stop plotutil_get_options'
        end if

	if( varname .ne. ' ' ) call string2ivar(varname,ivar3)

	if( ivar3 < 0 ) then
          write(6,*) 'variable name unknown: ',trim(varname)
          stop 'error stop plotutil_get_options'
	end if

	call setlev(layer)
	call setvar(ivar3)

	b2d = layer == 0

	end subroutine plotutil_get_options

!************************************************************
!************************************************************
!************************************************************
!************************************************************

	subroutine plotutil_check_options

	end subroutine plotutil_check_options

!************************************************************
!************************************************************
!************************************************************

!************************************************************
!************************************************************
!************************************************************

c***************************************************************

	subroutine classify_files

	use clo
	use shyfile
	use basin
	use nls

	logical bdebug
	character*80 file
	integer i
	integer nshy,nfem,nbas,nlgr,nstr,nunk,ngrd

	logical is_grd_file
	logical fem_file_is_fem_file
	logical filex

	bdebug = .true.
	bdebug = .false.

	nshy = 0
	nfem = 0
	nbas = 0
	nlgr = 0
	nstr = 0
	nunk = 0
	ngrd = 0

	nfile = clo_number_of_files()
	allocate(file_type(nfile))
	if( nfile == 0 ) return

	file_type = 'unknown'

	do i=1,nfile
	  call clo_get_file(i,file)
	  if( .not. filex(file) ) then
	    if( file(1:1) == '-' ) then
	      write(6,*) 'option ',trim(file),' in wrong place'
	      write(6,*) 'all options must lead files'
	      stop 'error stop classify_files'
	    else
	      write(6,*) 'file not existing: ',trim(file)
	      stop 'error stop classify_files'
	    end if
	  end if
	  if( shy_is_shy_file(file) ) then
	    file_type(i) = 'shy'
	    nshy = nshy + 1
	    shyfilename = file
	  else if( shy_is_lgr_file(file) ) then
	    file_type(i) = 'lgr'
	    nlgr = nlgr + 1
	    lgrfilename = file
	  else if( fem_file_is_fem_file(file) ) then
	    file_type(i) = 'fem'
	    nfem = nfem + 1
	    femfilename = file
	  else if( basin_is_basin(file) ) then
	    file_type(i) = 'bas'
	    nbas = nbas + 1
	    basfilename = file
	    basintype = 'bas'
	  else if( is_grd_file(file) ) then
	    file_type(i) = 'grd'
	    ngrd = ngrd + 1
	    basfilename = file
	    basintype = 'grd'
	  else if( nls_is_nls_file(file) ) then
	    file_type(i) = 'str'
	    nstr = nstr + 1
	  else
	    write(6,*) 'cannot determine file type: ',trim(file)
	    nunk = nunk + 1
	  end if
	end do

	if( nunk > 0 ) bdebug = .true.
	if( bdebug ) then
	  write(6,*) 'classify_files: ',nshy,nfem,nbas,nstr,nunk
	  write(6,*) 'shyfilename: ',trim(shyfilename)
	  write(6,*) 'femfilename: ',trim(femfilename)
	  write(6,*) 'lgrfilename: ',trim(lgrfilename)
	  write(6,*) 'basfilename: ',trim(basfilename)
	  write(6,*) 'basintype: ',trim(basintype)
	  do i=1,nfile
	    write(6,*) i,file_type(i)
	  end do
	end if

	if( nunk > 0 ) stop 'error stop classify_files'
	if( nshy > 1 ) then
	  write(6,*) 'cannot plot more than one SHY file'
	  stop 'error stop classify_files'
	end if
	if( nfem > 1 ) then
	  write(6,*) 'cannot plot more than one FEM file'
	  stop 'error stop classify_files'
	end if
	if( nlgr > 1 ) then
	  write(6,*) 'cannot plot more than one LGR file'
	  stop 'error stop classify_files'
	end if
	if( nbas > 0 .and. ngrd > 0 ) then
	  write(6,*) 'both BAS and GRD files given... cannot handle'
	  stop 'error stop classify_files'
	end if
	if( nshy > 0 .and. nfem > 0 ) then
	  write(6,*) 'both SHY and FEM files given... cannot handle'
	  stop 'error stop classify_files'
	end if
	if( nshy > 0 ) then
	  if( nbas > 0 .or. ngrd > 0 ) then
	    write(6,*) 'basin given but not needed... ignoring'
	  end if
	end if
	if( 
     +			      shyfilename == ' ' 
     +			.and. femfilename == ' ' 
     +			.and. lgrfilename == ' '
     +			.and. basfilename == ' '
     +	  ) then
	  write(6,*) 'no file given for plot...'
	  call clo_usage
	end if

	end subroutine classify_files

c***************************************************************

	subroutine read_str_files(ivar)

	use clo
	use nls

	integer ivar

	character*80 file
	integer i,iunit

	integer ifileo

	do i=1,nfile
	  if( file_type(i) == 'str' ) then
	    call clo_get_file(i,file)
	    if( bverb ) write(6,*) '--- reading str file ',ivar,trim(file)
            iunit = ifileo(0,file,'form','old')
            if( iunit <= 0 ) then
	      write(6,*) 'cannot read file: ',trim(file)
	      stop 'error stop read_str_files'
	    end if
	    call nlsa(iunit,ivar,bverb)
	    close(iunit)
	    if( bverb ) write(6,*) '--- finished reading str file '
	  end if
	end do

	end subroutine read_str_files

c***************************************************************

	subroutine init_nls_fnm

        call nlsina
        call fnminh

	end subroutine init_nls_fnm

!====================================================
	end module plotutil
!====================================================

c***************************************************************
c***************************************************************
c***************************************************************

