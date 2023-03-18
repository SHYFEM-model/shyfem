
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2020  Georg Umgiesser
!    Copyright (C) 2017-2018  Christian Ferrarin
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
! 30.07.2015	ggu	changed VERS_7_1_83
! 31.07.2015	ggu	changed VERS_7_1_84
! 14.09.2015	ggu	changed VERS_7_2_2
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 19.10.2015	ggu	changed VERS_7_3_6
! 18.12.2015	ggu	changed VERS_7_3_17
! 19.02.2016	ggu	changed VERS_7_5_2
! 19.02.2016	ggu	changed VERS_7_5_3
! 22.02.2016	ggu	handle catmode
! 22.03.2016	ggu	changed VERS_7_5_6
! 15.04.2016	ggu	handle gis files with substitution of colon
! 28.04.2016	ggu	changed VERS_7_5_9
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 08.09.2016	ggu	new options -map, -info, -areas, -dates
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.10.2016	ggu	changed VERS_7_5_19
! 21.03.2017	ggu	mode renamed to avermode
! 23.03.2017	ccf	line routines introduced
! 31.03.2017	ggu	changed VERS_7_5_24
! 16.05.2017	ggu	changed VERS_7_5_27
! 11.07.2017	ggu	changed VERS_7_5_30
! 05.10.2017	ggu	options rearranged, some only for SHY file
! 09.10.2017	ggu	new options, more uniform treatment
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 17.11.2017	ggu	changed VERS_7_5_38
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 24.05.2018	ccf	add outformat option off
! 06.06.2018	ggu	new calling format of shy_write_aver()
! 06.07.2018	ggu	changed VERS_7_5_48
! 31.08.2018	ggu	changed VERS_7_5_49
! 16.10.2018	ggu	changed VERS_7_5_50
! 25.10.2018	ccf	lagrangian options
! 16.02.2019	ggu	changed VERS_7_5_60
! 15.05.2019	ggu	new option -date0 (sdate0)
! 22.07.2019    ggu     new routines for handling time step check
! 13.12.2019    ggu     new option -checkrain (bcheckrain)
! 28.01.2020    ggu     new option -vorticity (bvorticity)
! 06.03.2020    ggu     -checkdt also for ext and flx files
! 21.05.2020    ggu     better handle copyright notice
! 05.11.2021    ggu     resample option added
! 25.01.2022    ggu     new option -grdcoord to plot fem grid
! 27.01.2022    ggu     new options -rmin,-rmax,-rfreq
! 07.03.2022    ggu     new options -changetime to shift time reference
! 03.05.2022    ggu     new option -nlgtype
! 05.12.2022    ggu     -facts and -offset working with TS files
! 09.03.2023    ggu     setting -facts, -offset in one subroutine
! 10.03.2023    ggu     map renamed to influencemap
!
!************************************************************

!====================================================
	module elabutil
!====================================================

	implicit none

	logical, save, private :: binitialized	= .false.
	double precision, parameter :: flag_p	= -999.
	double precision, save :: flag_d	= flag_p
	real, save :: flag			= flag_p

! binfo			as bverb, but stop after header
! bverbose		verbose - write time steps
! bquiet		be quiet (no time steps, no header)
! bsilent		possibly nothing
! bwrite		write min/max values
! 
! binfo implies bverb, bsilent implies bquiet
! bwrite implies bverb

	logical, save :: binfo			= .false.
	logical, save :: bverb			= .false.
	logical, save :: bquiet			= .false.
	logical, save :: bsilent		= .false.
	logical, save :: bwrite			= .false.
	logical, save :: bsdebug		= .false.

        character*80, save :: stmin		= ' '
        character*80, save :: stmax		= ' '
	logical, save :: binclusive		= .false.
	double precision, save :: difftime	= 0.

	integer, save :: rmin			= 1
	integer, save :: rmax			= 0
	integer, save :: rfreq			= 0

	logical, save :: bout			= .false.
        character*10, save :: outformat		= ' '
	integer, save :: catmode 		= 0

	logical, save :: bsplit			= .false.
	logical, save :: bsplitall		= .false.
        character*80, save :: nodelist		= ' '
        character*80, save :: nodefile		= ' '

	character*80, save :: sdate0		= ' '
	logical, save :: bconvert		= .false.
	logical, save :: bcheckdt		= .false.
	logical, save :: bcheckrain		= .false.

	logical, save :: bcondense		= .false.
	logical, save :: bchform		= .false.
	logical, save :: bgrd			= .false.
	logical, save :: bgrdcoord		= .false.

	logical, save :: bcheck			= .false.
        character*80, save :: scheck		= ' '

        character*80, save :: newstring		= ' '
        character*80, save :: factstring	= ' '
        character*80, save :: offstring		= ' '

        character*80, save :: regstring		= ' '
        character*80, save :: rbounds		= ' '
	integer, save :: regexpand		= -1

	logical, save :: bresample		= .false.
	logical, save :: baverbas		= .false.
	logical, save :: baver			= .false.
	logical, save :: baverdir		= .false.
	logical, save :: bsum			= .false.
	logical, save :: bmin			= .false.
	logical, save :: bmax			= .false.
	logical, save :: bstd			= .false.
	logical, save :: brms			= .false.
	logical, save :: bsumvar		= .false.
	double precision, save :: threshold	= flag_p
	real, save :: fact			= 1
	integer, save :: ifreq			= 0
	logical, save :: b2d			= .false.
	logical, save :: bvorticity		= .false.

	logical, save :: bdiff			= .false.
	real, save :: deps			= 0.

        character*80, save :: areafile		= ' '
        character*80, save :: datefile		= ' '
	logical, save :: binfluencemap 		= .false.

	integer, save :: nnodes			= 0
	integer, save, allocatable :: nodesi(:)		!internal nodes
	integer, save, allocatable :: nodese(:)		!external nodes
	logical, save :: bnode			= .false.
	logical, save :: bnodes			= .false.
	logical, save :: bcoord			= .false.
        character*80, save :: snode		= ' '
        character*80, save :: scoord		= ' '
        character*80, save :: sextract		= ' '

	integer, save :: istep			= 0
	integer, save :: avermode		= 0
	logical, save :: bthreshold		= .false.

        logical, save :: barea			= .false.

	logical, save :: boutput		= .false.
	logical, save :: bneedbasin		= .false.
	logical, save :: btrans			= .false.

	logical, save :: bopen			= .false.
	logical, save :: bdate			= .false.
	logical, save :: blgmean		= .false.
	logical, save :: blgdens		= .false.
	logical, save :: blgtype		= .false.
	logical, save :: blg2d			= .false.
	integer, save :: nlgtype		= 0

        character*80, save :: infile		= ' '
        integer, save, allocatable :: ieflag(:)
        integer, save, allocatable :: ikflag(:)

	logical, save :: bshowall		= .false.
	logical, save :: bshyfile		= .false.
	logical, save :: bflxfile		= .false.
	logical, save :: bextfile		= .false.
	logical, save :: bfemfile		= .false.
	logical, save :: btsfile		= .false.
	logical, save :: blgrfile		= .false.
	logical, save :: binputfile		= .false.

        character*10, save :: file_type		= ' '

	logical, save :: bcompat		= .true. !compatibility output

!====================================================
	contains
!====================================================

	subroutine elabutil_init(type,what)

	use clo

	character*(*) type
	character*(*), optional :: what

	character*80 program

	file_type = type
	program = 'shyelab'
	if( present(what) ) program = what

	call elabutil_set_options(type,program)
	call clo_parse_options
	call elabutil_get_options(type,program)

	binitialized = .true.

	end subroutine elabutil_init

!************************************************************

	subroutine elabutil_set_options(type,program)

	use clo

	character*(*) type
	character*(*) program

        character*10 vers
        character*80 version
        character*80 text

	if( binitialized ) return

        call get_shyfem_version(vers)
	version = '4.0' // ' (SHYFEM version ' // trim(vers) // ')'
	!write(6,*) 'ggguuu version: ',trim(version)

	!write(6,*) 'filetype: ',type

	if( type == 'NONE' ) then
          call clo_init(program,'shy-file',version)
	  bshowall = .true.
	else if( type == 'SHY' ) then
          call clo_init(program,'shy-file',version)
	  bshyfile = .true.
	else if( type == 'NOS' ) then
          call clo_init(program,'nos-file',version)
	else if( type == 'OUS' ) then
          call clo_init(program,'ous-file',version)
	else if( type == 'EXT' ) then
          call clo_init(program,'ext-file',version)
	  bextfile = .true.
	else if( type == 'FLX' ) then
          call clo_init(program,'flx-file',version)
	  bflxfile = .true.
	else if( type == 'FEM' ) then
          call clo_init(program,'fem-file',version)
	  bfemfile = .true.
	else if( type == 'TS' ) then
          call clo_init(program,'time-series-file',version)
	  btsfile = .true.
	else if( type == 'LGR' ) then
          call clo_init(program,'lgr-file',version)
	  blgrfile = .true.
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_set_options: unknown type'
	end if

	binputfile = bfemfile .or. btsfile

	text = 'returns info on or elaborates a ' //
     +			'shyfem file'
        call clo_add_info(text)

	call elabutil_set_general_options
	call elabutil_set_out_options
	call elabutil_set_extract_options

	call elabutil_set_ts_options
	call elabutil_set_fem_options
	call elabutil_set_reg_options
	call elabutil_set_shy_options
	call elabutil_set_diff_options
	call elabutil_set_hidden_shy_options
	call elabutil_set_lgr_options

	call elabutil_set_all_file_options

	end subroutine elabutil_set_options

!************************************************************

	subroutine elabutil_set_general_options

	use clo

        call clo_add_sep('general options')

        call clo_add_option('info',.false.,'only give info on header')
        call clo_add_option('verbose',.false.
     +				,'be more verbose, write time records')
        call clo_add_option('quiet',.false.
     +				,'do not write header information')
        call clo_add_option('silent',.false.,'do not write anything')
        call clo_add_option('write',.false.,'write min/max of records')
        call clo_add_option('debug',.false.,'write debug information')
	call clo_hide_option('debug')

        call clo_add_sep('time options')

        call clo_add_option('tmin time',' '
     +                  ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                  ,'only process up to time')
        call clo_add_option('inclusive',.false.
     +			,'output includes whole time period given')
        call clo_add_option('changetime difftime',0.
     +                  ,'add difftime to time record (difftime [s])')

	call clo_add_com('    time is either YYYY-MM-DD[::hh[:mm[:ss]]]')
	call clo_add_com('    or integer for relative time')

        call clo_add_option('rmin rec',1.
     +                  ,'only process starting from record rec')
        call clo_add_option('rmax rec',0.
     +                  ,'only process up to record rec')
        call clo_add_option('rfreq freq',1.
     +                  ,'only process every freq record')

	call clo_add_com('    rec in rmax can be negative')
	call clo_add_com('    this indicates rec records from the back')

	end subroutine elabutil_set_general_options

!************************************************************

	subroutine elabutil_set_out_options

	use clo

        call clo_add_sep('output options')

        call clo_add_option('out',.false.,'writes new file')

        call clo_add_option('outformat form','native','output format')
	call clo_add_com('    possible output formats are: '
     +				// 'shy|gis|fem|nc|off'
     +				// ' (Default native)')
	call clo_add_com('    not all formats are available for'
     +				// ' all file types')

        call clo_add_option('catmode cmode',0.,'concatenation mode'
     +				// ' for handeling more files')
	call clo_add_com('    possible values for cmode are: -1,0,+1'
     +				// ' (Default 0)')
	call clo_add_com('    -1: all of first file, '
     +				// 'then remaining of second')
	call clo_add_com('     0: simply concatenate files')
	call clo_add_com('    +1: first file until start of second, '
     +				// 'then all of second')

        call clo_add_option('proj projection',' '
     +                          ,'projection of coordinates')
        call clo_add_com('    projection is string consisting of '//
     +                          'mode,proj,params')
        call clo_add_com('    mode: +1: cart to geo,  -1: geo to cart')
        call clo_add_com('    proj: 1:GB, 2:UTM, 3:CPP')

	end subroutine elabutil_set_out_options

!************************************************************

	subroutine elabutil_set_extract_options

	use clo

        call clo_add_sep('extract options')

        call clo_add_option('split',.false.,'split file for variables')

	if( bshowall .or. bflxfile .or. bextfile ) then
          call clo_add_option('splitall',.false.
     +		,'splits file (EXT and FLX) for extended data')
	end if

	if( bshowall .or. binputfile ) then
          call clo_add_option('check period',' '
     +				,'checks data over period')
          call clo_add_com('  period can be '//
     +				'all,year,month,week,day,none')
	end if

	!if( bshowall .or. binputfile .or. bshyfile ) then
	!if( bshowall .or. binputfile ) then
          call clo_add_option('checkdt',.false.
     +			,'check for change of time step')
	!end if

	if( bshowall .or. binputfile ) then
          call clo_add_option('checkrain',.false.
     +			,'check for yearly rain (if file contains rain)')
	end if

	if( bshowall .or. bshyfile ) then
          call clo_add_option('node nlist',' '
     +			,'extract vars of nodes in list')
	  call clo_add_com('    nlist is a comma separated list of nodes'
     +				//' to be extracted')
          call clo_add_option('nodes nfile',' '
     +			,'extract vars at nodes given in file nfile')
	  call clo_add_com('    nfile is a file with nodes'
     +				//' to be extracted')
          call clo_add_option('extract recs',' '
     +			,'extract records specified in recs')
	  call clo_add_com('    recs is either a comma separated list'
     +				//' like r1,r2,r3')
	  call clo_add_com('    or in the format istart..ifreq..iend'
     +				//' (..iend may be missing)')
	end if

	if( bshowall .or. bshyfile .or. bfemfile ) then
          call clo_add_option('coord coord',' ','extract coordinate')
          call clo_add_com('    coord is x,y of point to extract')
	end if

	if( bshowall .or. btsfile ) then
          call clo_add_sep('time series options')
          call clo_add_option('convert',.false.
     +			,'convert time column to ISO string')
          call clo_add_option('date0',' '
     +			,'reference date for conversion of time column')
	end if

	end subroutine elabutil_set_extract_options

!************************************************************

	subroutine elabutil_set_facts_options

	use clo

	if( clo_has_option('facts') ) return

        call clo_add_option('facts fstring',' '
     +			,'apply factors to data in data-file')
        call clo_add_option('offset ostring',' '
     +			,'apply factors to data in data-file')
        call clo_add_com('    fstring and ostring is comma'
     +			// ' separated factors,'
     +                  // ' empty for no change')

	end subroutine elabutil_set_facts_options

!************************************************************

	subroutine elabutil_set_ts_options

	use clo

	if( .not. bshowall .and. .not. btsfile ) return

	call elabutil_set_facts_options

	end subroutine elabutil_set_ts_options

!************************************************************

	subroutine elabutil_set_fem_options

	use clo

	if( .not. bshowall .and. .not. bfemfile ) return

        call clo_add_sep('specific FEM file options')

        call clo_add_option('condense',.false.
     +			,'condense file data into one node')
        call clo_add_option('chform',.false.
     +			,'change output format form/unform of FEM file')
        call clo_add_option('grd',.false.
     +			,'write GRD file from data in FEM file')
        call clo_add_option('grdcoord',.false.
     +			,'write regular coordinates in GRD format')
        call clo_add_option('nodei node',' ','extract internal node')
        call clo_add_com('    node is internal numbering in fem file'
     +                  //' or ix,iy of regular grid')
        call clo_add_option('newstring sstring',' '
     +			,'substitute string description in fem-file')
        call clo_add_com('    sstring is comma separated strings,'
     +                  //' empty for no change')

	call elabutil_set_facts_options

	end subroutine elabutil_set_fem_options

!************************************************************

	subroutine elabutil_set_reg_options

	use clo

	logical bregopt

	bregopt = bshowall .or. bfemfile .or. bshyfile .or. blgrfile
	if( .not. bregopt ) return

        call clo_add_sep('regular output file options')

	call clo_add_option('reg rstring',' ','regular interpolation')
	call clo_add_option('resample bounds',' ','resample regular grid')
	call clo_add_option('regexpand iexp',-1,'expand regular grid')

	call clo_add_com('    rstring is: dx[,dy[,x0,y0,x1,y1]]')
	call clo_add_com('    if only dx is given -> dy=dx')
	call clo_add_com('    if only dx,dy are given -> bounds computed')
	call clo_add_com('    bounds is: x0,y0,x1,y1')
	call clo_add_com('    iexp>0 expands iexp cells, =0 whole grid')
	call clo_add_com('    resample should be used with regexpand')

	end subroutine elabutil_set_reg_options

!************************************************************

	subroutine elabutil_set_shy_options

	use clo

	if( .not. bshowall .and. .not. bshyfile ) return

        call clo_add_sep('transformation of SHY file options')

        call clo_add_option('averbas',.false.,'average over basin')
        call clo_add_option('aver',.false.,'average over records')
        call clo_add_option('averdir',.false.,'average for directions')
        call clo_add_option('sum',.false.,'sum over records')
        call clo_add_option('min',.false.,'minimum of records')
        call clo_add_option('max',.false.,'maximum of records')
	call clo_add_option('std',.false.,'standard deviation of records')
        call clo_add_option('rms',.false.,'root mean square of records')
        call clo_add_option('sumvar',.false.,'sum over variables')
	call clo_add_option('threshold t',flag
     +				,'compute records over threshold t')
	call clo_add_option('fact fact',1.,'multiply values by fact')
	call clo_add_option('freq n',0.
     +			,'frequency for aver/sum/min/max/std/rms')

	call clo_add_option('2d',.false.,'average vertically to 2d field')
	call clo_add_option('vorticity',.false.
     +			,'compute vorticity for hydro file')

	end subroutine elabutil_set_shy_options

!************************************************************

	subroutine elabutil_set_diff_options

	use clo

	if( .not. bshowall .and. .not. bshyfile ) return

        call clo_add_sep('difference of SHY file options')

	call clo_add_option('diff',.false.
     +			,'check if 2 files are different')
	call clo_add_option('diffeps deps',0.
     +			,'files differ by more than deps (default 0)')

	call clo_add_com('    this option needs two files' //
     +				' and exits at difference')
	call clo_add_com('    with -out writes difference to out file')

	end subroutine elabutil_set_diff_options

!************************************************************

	subroutine elabutil_set_hidden_shy_options

	use clo

	!if( .not. clo_want_extended_help() ) return

	call clo_hide_next_options

        call clo_add_sep('extra options')

        call clo_add_option('areas line-file',' '
     +			,'line delimiting areas for -averbas option')
        call clo_add_option('dates date-file',' '
     +			,'give dates for averaging in file')
        call clo_add_option('influencemap',.false.
     +			,'computes influence map from multi-conz')

	call clo_show_next_options

	end subroutine elabutil_set_hidden_shy_options

!************************************************************

	subroutine elabutil_set_all_file_options

	use clo

	if( .not. bshowall ) return

	call clo_add_com('All options for all file types are shown')
	call clo_add_com('Not all options are available for all files')
	call clo_add_com('To show options for a specific file use: '
     +		//'shyelab -h file')

	end subroutine elabutil_set_all_file_options

!************************************************************

        subroutine elabutil_set_lgr_options

        use clo

        if( .not. bshowall .and. .not. blgrfile ) return

        call clo_add_sep('specific LGR file options')

        call clo_add_option('lgmean',.false.,'extract mean position'
     +                  //' and age in function')
        call clo_add_com('     of particle type. It '
     +                  //' write files lagrange_mean_traj.*')
        call clo_add_option('lgdens',.false.,'compute distribution of'
     +                  //' particle density and age')
        call clo_add_com('     either on nodes or on'
     +                  //' regular grid (using -reg option)')
	call clo_add_option('lg2d',.false.,'sum particles vertically' 
     +                  //' when computing the ')
        call clo_add_com('     particle density/age')
        call clo_add_option('lgtype',.false.,'compute density per type')
        call clo_add_option('nlgtype max-type',0
     +			,'max number of types to compute')

        end subroutine elabutil_set_lgr_options

!************************************************************
!************************************************************
!************************************************************

	subroutine elabutil_get_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	character*20 ftype,flow
        character*80 text

	if( binitialized ) return

        call clo_get_option('info',binfo)
        call clo_get_option('verbose',bverb)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('write',bwrite)
        call clo_get_option('debug',bsdebug)

        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)
        call clo_get_option('inclusive',binclusive)
        call clo_get_option('changetime',difftime)

        call clo_get_option('rmin',rmin)
        call clo_get_option('rmax',rmax)
        call clo_get_option('rfreq',rfreq)

        call clo_get_option('out',bout)
        call clo_get_option('outformat',outformat)
        call clo_get_option('catmode',catmode)

        call clo_get_option('split',bsplit)
	if( bshowall .or. bflxfile .or. bextfile ) then
          call clo_get_option('splitall',bsplitall)
	end if
        call clo_get_option('checkdt',bcheckdt)
	if( bshowall .or. binputfile ) then
          call clo_get_option('check',scheck)
          call clo_get_option('checkrain',bcheckrain)
	end if
	if( bshowall .or. bshyfile ) then
          call clo_get_option('node',nodelist)
          call clo_get_option('nodes',nodefile)
          call clo_get_option('coord',scoord)
          call clo_get_option('extract',sextract)
	end if

	if( bshowall .or. btsfile ) then
          call clo_get_option('convert',bconvert)
          call clo_get_option('date0',sdate0)
          call clo_get_option('facts',factstring)
          call clo_get_option('offset',offstring)
	end if

	if( bshowall .or. bfemfile ) then
          call clo_get_option('condense',bcondense)
          call clo_get_option('chform',bchform)
          call clo_get_option('grd',bgrd)
          call clo_get_option('grdcoord',bgrdcoord)
          call clo_get_option('nodei',snode)
          call clo_get_option('coord',scoord)
          call clo_get_option('newstring',newstring)
          call clo_get_option('facts',factstring)
          call clo_get_option('offset',offstring)
	end if

	if( bshowall .or. bfemfile .or. bshyfile .or. blgrfile ) then
          call clo_get_option('reg',regstring)
          call clo_get_option('regexpand',regexpand)
          call clo_get_option('resample',rbounds)
	end if

	if( bshowall .or. bshyfile ) then
          call clo_get_option('averbas',baverbas)
          call clo_get_option('aver',baver)
          call clo_get_option('averdir',baverdir)
          call clo_get_option('sum',bsum)
          call clo_get_option('min',bmin)
          call clo_get_option('max',bmax)
          call clo_get_option('std',bstd)
          call clo_get_option('rms',brms)
          call clo_get_option('sumvar',bsumvar)
          call clo_get_option('threshold',threshold)
          call clo_get_option('fact',fact)
          call clo_get_option('freq',ifreq)
          call clo_get_option('2d',b2d)
          call clo_get_option('vorticity',bvorticity)
	end if

	if( bshowall .or. bshyfile ) then
          call clo_get_option('diff',bdiff)
          call clo_get_option('diffeps',deps)
	end if

	if( bshowall .or. bshyfile ) then
          call clo_get_option('areas',areafile)
          call clo_get_option('dates',datefile)
          call clo_get_option('influencemap',binfluencemap)
	end if

	if( bshowall .or. blgrfile ) then
          call clo_get_option('lgmean',blgmean)
          call clo_get_option('lgdens',blgdens)
          call clo_get_option('lgtype',blgtype)
          call clo_get_option('nlgtype',nlgtype)
          call clo_get_option('lg2d',blg2d)
	end if

!-------------------------------------------------------------------
! write copyright
!-------------------------------------------------------------------

	ftype = type
	if( type == 'NONE' ) then
	  ftype = 'SHY'
	else if( type == 'SHY' ) then
	else if( type == 'NOS' ) then
	else if( type == 'OUS' ) then
	else if( type == 'EXT' ) then
	else if( type == 'FLX' ) then
	else if( type == 'FEM' ) then
	else if( type == 'TS' ) then
	else if( type == 'LGR' ) then
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_get_options: unknown type'
	end if

	call shyfem_set_short_copyright(bquiet)
        if( .not. bsilent ) then
	  flow = ftype
	  call to_lower(flow)
	  text = trim(flow) // 'elab - elaborates ' 
     +				// trim(ftype) // ' files'
          call shyfem_copyright(trim(text))
	end if

!-------------------------------------------------------------------
! get and check input file(s)
!-------------------------------------------------------------------

        call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

!-------------------------------------------------------------------
! set dependent parameters
!-------------------------------------------------------------------

        bnode = nodelist .ne. ' '
        bnodes = nodefile .ne. ' '
        bcoord = scoord .ne. ' '

	barea = ( areafile /= ' ' )
	bcheck = ( scheck /= ' ' )
	bresample = ( rbounds /= ' ' )

        boutput = bout
        boutput = boutput .or. b2d .or. bvorticity
        boutput = boutput .or. bsplit
	boutput = boutput .or. outformat /= 'native'
        boutput = boutput .or. bsumvar
        boutput = boutput .or. binfluencemap
        boutput = boutput .or. bresample
        boutput = boutput .or. newstring /= ' '
        boutput = boutput .or. sextract /= ' '

        !btrans is added later
	!if( bsumvar ) boutput = .false.

        bneedbasin = b2d .or. baverbas .or. bvorticity
        bneedbasin = bneedbasin .or. bnode .or. bnodes
        bneedbasin = bneedbasin .or. bcoord
	bneedbasin = bneedbasin .or. outformat == 'gis'
	bneedbasin = bneedbasin .or. ( type == 'OUS' .and. bsplit )

	bthreshold = ( threshold /= flag )

	if( binfo ) bverb = .true.
	if( bwrite ) bverb = .true.
	if( bsdebug ) bverb = .true.
	if( bsilent ) bquiet = .true.

	if( bsplitall ) bsplit = .true.

	end subroutine elabutil_get_options

!************************************************************

	subroutine elabutil_check_options

	integer ic

	ic = count( (/b2d,bsplit,bsumvar,btrans,bvorticity/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-2d -split -sumvar -vorticity'
	  write(6,*) '-aver -sum -min -max -std -rms'
	  write(6,*) '-threshold -averdir'
	  stop 'error stop elabutil_check_options: incompatible options'
	end if

	end subroutine elabutil_check_options

!************************************************************
!************************************************************
!************************************************************

	subroutine elabutil_set_averaging(nvar)

	integer nvar

	integer ic

        avermode = 0
        btrans = .false.

	ic = count( (/baver,bsum,bmin,bmax,bstd,brms
     +				,bthreshold,baverdir/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-aver -sum -min -max -std -rms'
	  write(6,*) '-threshold -averdir'
	  stop 'error stop elabutil_set_averaging: incompatible options'
	end if

        if( baver ) avermode = 1
        if( bsum )  avermode = 2
        if( bmin )  avermode = 3
        if( bmax )  avermode = 4
        if( bstd )  avermode = 5
        if( brms )  avermode = 6
        if( bthreshold )  avermode = 7
        if( baverdir ) avermode = 8

        if( avermode > 0 ) then             !prepare for averaging
          btrans = .true.
          if( bsumvar ) then            !sum over variables
	    if( ifreq /= 0 ) then
	      write(6,*) 'For option -sumvar cannot use value for -freq'
	      write(6,*) 'freq = ',ifreq
	      stop 'error stop elabutil_set_averaging: freq'
	    end if
            istep = 1
          else if( ifreq .ge. 0 ) then  !normal averaging
            istep = 1
          else                          !accumulate every -ifreq record
            istep = -ifreq
            ifreq = 0
          end if
	end if

        if( bsumvar ) then            !sum over variables
	  if( ifreq /= 0 ) then
	    write(6,*) 'For option -sumvar cannot use value for -freq'
	    write(6,*) 'freq = ',ifreq
	    stop 'error stop elabutil_set_averaging: freq'
	  end if
          istep = 1
	end if

	end subroutine elabutil_set_averaging

!====================================================
	end module elabutil
!====================================================

        subroutine outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

c averages vertically

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        real hm3v(3,nel)
        real hev(nel)
        real hkv(nkn)

        integer k,ie,ii
        real h,hm

	hkv = -huge(1.)

        do ie=1,nel
	  hm = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            h = hm3v(ii,ie)
            hkv(k) = max(hkv(k),h)
	    hm = hm + h
          end do
	  hev(ie) = hm / 3.
        end do

        end

c***************************************************************

        subroutine outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)

c averages vertically

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        real hm3v(3,nel)
        real hev(nel)
        real hkv(nkn)

        integer k,ie,ii
        real h

        do ie=1,nel
          h = hev(ie)			!old way
          do ii=1,3
            k = nen3v(ii,ie)
            h = hm3v(ii,ie)		!new way
            hkv(k) = max(hkv(k),h)
          end do
        end do

        end

c***************************************************************

        subroutine depth_stats(nkn,nlvddi,ilhkv)

c       computes statistics on levels

        implicit none

        integer nkn
        integer nlvddi
        integer ilhkv(nkn)

        integer count(nlvddi)
        integer ccount(nlvddi)

        integer nlv,lmax,l,k,nc,ll

        nlv = 0
        do l=1,nlvddi
          count(l) = 0
          ccount(l) = 0
        end do

        do k=1,nkn
          lmax = ilhkv(k)
          if( lmax .gt. nlvddi ) stop 'error stop depth_stats: lmax'
          count(lmax) = count(lmax) + 1
          nlv = max(nlv,lmax)
        end do

        do l=nlv,1,-1
          nc = count(l)
          do ll=1,l
            ccount(ll) = ccount(ll) + nc
          end do
        end do

        nc = 0
        write(6,*) 'statistics for layers: ',nlv
        write(6,*) '      layer       count accumulated'
        do l=1,nlv
          if( count(l) > 0 ) then
            write(6,*) l,count(l),ccount(l)
            nc = nc + count(l)
          end if
        end do

	if( nc /= nkn ) then
          write(6,*) 'total count: ',nc
          write(6,*) 'total nodes: ',nkn
	  stop 'error stop depth_stats: data mismatch'
	end if

        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine shy_write_aver(aline,nvar,iv,ivar
     +				,cmin,cmax,cmed,cstd,atot,vtot)

c writes basin average to file

        implicit none

	character*20 aline
        integer nvar,iv,ivar
        real cmin,cmax,cmed,cstd,atot,vtot

	integer iu
        real totmass
	character*80 filename
	integer, save :: icall = 0
	integer, save, allocatable :: ius(:)

	if( .not. allocated(ius) ) then
          allocate(ius(0:nvar))
          ius = 0
        end if

	if( ius(iv) == 0 ) then
          call ivar2filename(ivar,filename)
          call make_iunit_name(filename,'','0d',0,iu)
          ius(iv) = iu
          write(iu,'(a)') '#      date_and_time    minimum'//
     +                  '    average    maximum        std'//
     +                  '         total'
	  if( iv == 1 ) then
            !call ivar2filename(0,filename)
	    filename = 'volume_and_area'
            call make_iunit_name(filename,'','0d',0,iu)
            ius(0) = iu
            write(iu,'(a)') '#      date_and_time        volume'//
     +                  '          area'
	  end if
	end if

	totmass = vtot
        if( ivar /= 1 ) totmass = cmed * vtot

	iu = ius(iv)
        write(6,2234) aline,ivar,cmin,cmed,cmax,cstd,totmass
        write(iu,2235) aline,cmin,cmed,cmax,cstd,totmass
	if( iv == 1 ) then
	  iu = ius(0)
	  write(iu,2236) aline,vtot,atot
	end if

	return
 2234   format(a20,i5,4f10.2,e14.6)
 2235   format(a20,4f11.3,e14.6)
 2236   format(a20,2e14.6)
        end

c***************************************************************

        subroutine write_aver(it,ivar,cmin,cmax,cmed,vtot)

c writes basin average to file

        implicit none

        integer it,ivar
        real cmin,cmax,cmed,vtot

        real totmass

        totmass = cmed * vtot

        write(6,1234) it,ivar,cmin,cmed,cmax,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,totmass
        write(100,1236) it,vtot

	return
 1234   format(i10,i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
 1236   format(i10,e14.6)
        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine ilhe2k(nkn,nel,nen3v,ilhv,ilhkv)

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ilhv(nel)
	integer ilhkv(nkn)

	integer ie,ii,k,l

	ilhkv = 0

	do ie=1,nel
	  l = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ilhkv(k) = max(l,ilhkv(k))
	  end do
	end do

	end

c***************************************************************

	subroutine ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)

c create ilhv -> result is not exact and must be adjusted

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ilhkv(nkn)
	integer ilhv(nel)

	integer ie,ii,k,lmax

	ilhv = 0

	do ie=1,nel
	  lmax = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    lmax = max(lmax,ilhkv(k))
	  end do
	  ilhv(ie) = lmax
	end do

	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine open_shy_file(file,status,nunit)

c open SHY file
c
c nunit is 0 if no other file exists

	use clo

        implicit none

        character*(*) status
        integer nunit

        character*80 file
        integer ifileo

        nunit = 0
        if( file == ' ' ) return

        nunit = ifileo(0,file,'unform',status)

        if( nunit .le. 0 ) then
          write(6,*) 'file: ',trim(file)
          stop 'error stop open_shy_file: opening file'
        end if

        end

c***************************************************************

	function concat_cycle_a(atime,atlast,atstart)

	use elabutil

c decides if with concatenation we have to use record or not

	implicit none

	logical concat_cycle_a
	double precision atime,atlast,atstart

	character*20 dline

	concat_cycle_a = .false.

        !write(66,*) 'ggu: ',atime,atlast,atstart

        if( catmode < 0 ) then
          if( atime <= atlast ) then
	    call dts_format_abs_time(atime,dline)
            write(6,*) 'skipping record: ',atime,dline
	    concat_cycle_a = .true.
          end if
        else if( catmode > 0 .and. atstart /= -1 ) then
          if( atime >= atstart ) then
	    call dts_format_abs_time(atime,dline)
            write(6,*) 'skipping record: ',atime,dline
	    concat_cycle_a = .true.
          end if
        end if

	end

c***************************************************************

	function concat_cycle(it,itold,itstart,nrec)

	use elabutil

c decides if with concatenation we have to use record or not

	implicit none

	logical concat_cycle
	integer it,itold,itstart
	integer nrec

	concat_cycle = .false.

        !write(66,*) 'ggu: ',it,itold,itstart,nrec

        if( catmode < 0 .and. nrec /= 1 ) then
          if( it <= itold ) then
            write(6,*) 'skipping record: ',it
            it = itold
	    concat_cycle = .true.
          end if
        else if( catmode > 0 .and. itstart /= -1 ) then
          if( it >= itstart ) then
            write(6,*) 'skipping record: ',it
	    concat_cycle = .true.
          end if
        end if

	end

c***************************************************************
c***************************************************************
c***************************************************************
c output utilities
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine compute_range(n,string)

	implicit none

	integer n
	character*(*) string

        string = '1'
        if( n > 1 ) then
          write(string,'(i5)') n
          string = adjustl(string)
          string = '1-' // trim(string)
        end if

	end

c***************************************************************

	subroutine write_vars(nvar,ivars)

	use shyfem_strings

	implicit none

	integer nvar
	integer ivars(nvar)

	integer iv,ivar
	character*20 short
	character*80 full

	do iv=1,nvar
          ivar = ivars(iv)
	  call ivar2filename(ivar,short)
          !call strings_get_short_name(ivar,short)
          call strings_get_full_name(ivar,full)
          write(6,*) '  ',short,trim(full)
	end do

	end

c***************************************************************

	subroutine write_extra_vars(nvar,ivars,post,descrp)

	use shyfem_strings

	implicit none

	integer nvar
	integer ivars(nvar)
	character*(*) post,descrp

	integer iv,ivar
	character*20 short
	character*80 full

	do iv=1,nvar
          ivar = ivars(iv)
	  call ivar2filename(ivar,short)
	  short = trim(short) // trim(post)
          call strings_get_full_name(ivar,full)
	  full = trim(full) // trim(descrp)
          write(6,*) '  ',short,trim(full)
	end do

	end

c***************************************************************

	subroutine write_special_vars(nvar,what,descrp)

	use shyfem_strings

	implicit none

	integer nvar
	character*(*) what(nvar)
	character*(*) descrp(nvar)

	integer iv
	character*20 short
	character*80 full

	do iv=1,nvar
	  short = what(iv)
	  full = descrp(iv)
          write(6,*) '  ',short,trim(full)
	end do

	end

c***************************************************************

	subroutine write_grd_coords(regpar)

	implicit none

	real regpar(7)

	integer nx,ny,ix,iy,n,l
	real x0,y0,x1,y1,dx,dy,x,y
	character*80 grdfile

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	x1 = x0+(nx-1)*dx
	y1 = y0+(ny-1)*dy

	grdfile = 'bound.fem.grd'
	write(6,*) 'writing file ',trim(grdfile)
	open(11,file=grdfile,status='unknown',form='formatted')
	write(11,1001) 1,1,0,x0,y0
	write(11,1001) 1,2,0,x0,y1
	write(11,1001) 1,3,0,x1,y1
	write(11,1001) 1,4,0,x1,y0
	write(11,1003) 3,1,0,5,1,2,3,4,1
	close(11)

	grdfile = 'cross.fem.grd'
	write(6,*) 'writing file ',trim(grdfile)
	open(11,file=grdfile,status='unknown',form='formatted')
	n = 0
	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  do ix=1,nx
	    n = n + 1
	    x = x0 + (ix-1)*dx
	    write(11,1001) 1,n,0,x,y
	  end do
	end do
	close(11)

	grdfile = 'grid.fem.grd'
	write(6,*) 'writing file ',trim(grdfile)
	open(11,file=grdfile,status='unknown',form='formatted')
	n = 0
	l = 0
	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  n = n + 1
	  write(11,1001) 1,n,0,x0,y
	  n = n + 1
	  write(11,1001) 1,n,0,x1,y
	  l = l + 1
	  write(11,1003) 3,l,0,2,n-1,n
	end do
	do ix=1,nx
	  x = x0 + (ix-1)*dx
	  n = n + 1
	  write(11,1001) 1,n,0,x,y0
	  n = n + 1
	  write(11,1001) 1,n,0,x,y1
	  l = l + 1
	  write(11,1003) 3,l,0,2,n-1,n
	end do
	close(11)

	return
 1001	format(i1,i10,i4,2f14.6)
 1003	format(i1,i10,i4,i6,5i6)
	end

c***************************************************************

