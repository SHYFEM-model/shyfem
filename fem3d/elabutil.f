!
! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 22.02.2016	ggu	handle catmode
! 15.04.2016	ggu	handle gis files with substitution of colon
! 08.09.2016	ggu	new options -map, -info, -areas, -dates
! 21.03.2017	ggu	mode renamed to avermode
! 23.03.2017	ccf	line routines introduced
! 05.10.2017	ggu	options rearranged, some only for SHY file
! 09.10.2017	ggu	new options, more uniform treatment
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

	logical, save :: bout			= .false.
        character*10, save :: outformat		= ' '
	integer, save :: catmode 		= 0

	logical, save :: bsplit			= .false.
	logical, save :: bsplitall		= .false.
        character*80, save :: nodelist		= ' '
        character*80, save :: nodefile		= ' '

	logical, save :: bconvert		= .false.
	logical, save :: bcheckdt		= .false.

	logical, save :: bcondense		= .false.
	logical, save :: bchform		= .false.
	logical, save :: bgrd			= .false.
        character*80, save :: snode		= ' '
        character*80, save :: scoord		= ' '

	logical, save :: bcheck			= .false.
        character*80, save :: scheck		= ' '

        character*80, save :: newstring		= ' '
        character*80, save :: factstring	= ' '

        character*80, save :: regstring		= ' '
	integer, save :: regexpand		= -1

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

	logical, save :: bdiff			= .false.
	real, save :: deps			= 0.

        character*80, save :: areafile		= ' '
        character*80, save :: datefile		= ' '
	logical, save :: bmap 			= .false.

	integer, save :: nnodes			= 0
	integer, save, allocatable :: nodes(:)
	integer, save, allocatable :: nodese(:)
	logical, save :: bnode			= .false.
	logical, save :: bnodes			= .false.

	integer, save :: istep			= 0
	integer, save :: avermode		= 0
	logical, save :: bthreshold		= .false.

        logical, save :: barea			= .false.

	logical, save :: boutput		= .false.
	logical, save :: bneedbasin		= .false.
	logical, save :: btrans			= .false.

	logical, save :: bopen			= .false.
	logical, save :: bdate			= .false.

        character*80, save :: infile		= ' '
        integer, save, allocatable :: ieflag(:)
        integer, save, allocatable :: ikflag(:)

	logical, save :: bshowall		= .false.
	logical, save :: bshyfile		= .false.
	logical, save :: bflxfile		= .false.
	logical, save :: bextfile		= .false.
	logical, save :: bfemfile		= .false.
	logical, save :: btsfile		= .false.
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

	call elabutil_set_fem_options
	call elabutil_set_femreg_options
	call elabutil_set_shy_options
	call elabutil_set_diff_options
	call elabutil_set_hidden_shy_options

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

	call clo_add_com('    time is either YYYY-MM-DD[::hh[:mm[:ss]]]')
	call clo_add_com('    or integer for relative time')

	end subroutine elabutil_set_general_options

!************************************************************

	subroutine elabutil_set_out_options

	use clo

        call clo_add_sep('output options')

        call clo_add_option('out',.false.,'writes new file')

        call clo_add_option('outformat form','native','output format')
	call clo_add_com('    possible output formats are: shy|gis|fem|nc'
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
          call clo_add_option('checkdt',.false.
     +			,'check for change of time step')
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
	end if

	if( bshowall .or. btsfile ) then
          call clo_add_sep('time series options')
          call clo_add_option('convert',.false.
     +			,'convert time column to ISO string')
	end if

	end subroutine elabutil_set_extract_options

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
        call clo_add_option('nodei node',' ','extract internal node')
        call clo_add_option('coord coord',' ','extract coordinate')
        call clo_add_com('    node is internal numbering in fem file'
     +                  //' or ix,iy of regular grid')
        call clo_add_com('    coord is x,y of point to extract')
        call clo_add_option('newstring sstring',' '
     +			,'substitute string description in fem-file')
        call clo_add_com('    sstring is comma separated strings,'
     +                  //' empty for no change')
        call clo_add_option('facts fstring',' '
     +			,'apply factors to data in fem-file')
        call clo_add_com('    fstring is comma separated factors,'
     +                  //' empty for no change')

	end subroutine elabutil_set_fem_options

!************************************************************

	subroutine elabutil_set_femreg_options

	use clo

	logical bregopt

	bregopt = bshowall .or. bfemfile .or. bshyfile
	if( .not. bregopt ) return

        call clo_add_sep('regular grid FEM file options')

        call clo_add_option('reg rstring',' ','regular interpolation')
        call clo_add_option('regexpand iexp',-1,'expand regular grid')

	call clo_add_com('    rstring is: dx[,dy[,x0,y0,x1,y1]]')
	call clo_add_com('    if only dx is given -> dy=dx')
	call clo_add_com('    if only dx,dy are given -> bounds computed')
	call clo_add_com('    iexp>0 expands iexp cells, =0 whole grid')

	end subroutine elabutil_set_femreg_options

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
        call clo_add_option('map',.false.
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

        call clo_get_option('out',bout)
        call clo_get_option('outformat',outformat)
        call clo_get_option('catmode',catmode)

        call clo_get_option('split',bsplit)
	if( bshowall .or. bflxfile .or. bextfile ) then
          call clo_get_option('splitall',bsplitall)
	end if
	if( bshowall .or. binputfile ) then
          call clo_get_option('check',scheck)
          call clo_get_option('checkdt',bcheckdt)
	end if
	if( bshowall .or. bshyfile ) then
          call clo_get_option('node',nodelist)
          call clo_get_option('nodes',nodefile)
	end if

	if( bshowall .or. btsfile ) then
          call clo_get_option('convert',bconvert)
	end if

	if( bshowall .or. bfemfile ) then
          call clo_get_option('condense',bcondense)
          call clo_get_option('chform',bchform)
          call clo_get_option('grd',bgrd)
          call clo_get_option('nodei',snode)
          call clo_get_option('coord',scoord)
          call clo_get_option('newstring',newstring)
          call clo_get_option('facts',factstring)
	end if

	if( bshowall .or. bfemfile .or. bshyfile ) then
          call clo_get_option('reg',regstring)
          call clo_get_option('regexpand',regexpand)
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
	end if

	if( bshowall .or. bshyfile ) then
          call clo_get_option('diff',bdiff)
          call clo_get_option('diffeps',deps)
	end if

	if( bshowall .or. bshyfile ) then
          call clo_get_option('areas',areafile)
          call clo_get_option('dates',datefile)
          call clo_get_option('map',bmap)
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
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_get_options: unknown type'
	end if

        if( .not. bsilent ) then
	  flow = ftype
	  call to_lower(flow)
	  text = trim(flow) // 'elab - Elaborate ' 
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

	barea = ( areafile /= ' ' )
	bcheck = ( scheck /= ' ' )

        boutput = bout
        boutput = boutput .or. b2d
        boutput = boutput .or. bsplit
	boutput = boutput .or. outformat /= 'native'
        boutput = boutput .or. bsumvar
        boutput = boutput .or. bmap
        boutput = boutput .or. newstring /= ' '

        !btrans is added later
	!if( bsumvar ) boutput = .false.

        bneedbasin = b2d .or. baverbas .or. bnode .or. bnodes
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

	ic = count( (/b2d,bsplit,bsumvar,btrans/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-2d -split -sumvar'
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

        subroutine shy_write_aver(dtime,ivar,cmin,cmax,cmed,cstd,vtot)

c writes basin average to file

        implicit none

        double precision dtime
        integer ivar
        real cmin,cmax,cmed,cstd,vtot

	integer it
        real totmass

	it = nint(dtime)
	totmass = vtot
        if( ivar /= 1 ) totmass = cmed * vtot

        !write(6,1234) it,ivar,cmin,cmed,cmax,cstd,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,cstd,totmass
        write(100,1236) it,vtot

        write(6,2234) dtime,ivar,cmin,cmed,cmax,cstd,totmass
        write(200+ivar,2235) dtime,cmin,cmed,cmax,cstd,totmass
        write(200,2236) dtime,vtot

	return
 1234   format(i10,i10,4f12.4,e14.6)
 1235   format(i10,4f12.4,e14.6)
 1236   format(i10,e14.6)
 2234   format(f15.2,i5,4f12.4,e14.6)
 2235   format(f15.2,4f12.4,e14.6)
 2236   format(f15.2,e14.6)
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

