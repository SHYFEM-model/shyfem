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
!
!************************************************************

!====================================================
	module plotutil
!====================================================

	implicit none

	logical, save, private :: binitialized = .false.
	double precision, parameter :: flag = -999.

	logical, save :: bout
	logical, save :: baverbas
	logical, save :: baver
	logical, save :: baverdir
	logical, save :: bsum
	logical, save :: bmin
	logical, save :: bmax
	logical, save :: bstd
	logical, save :: brms
	logical, save :: bsumvar
	logical, save :: bsplit
	logical, save :: b2d

	logical, save :: bmem		= .false.
	logical, save :: bask		= .false.
	logical, save :: bverb
	logical, save :: bwrite
	logical, save :: bquiet
	!logical, save :: bdate

	integer, save :: ifreq
	integer, save :: tmin
	integer, save :: tmax

	logical, save :: bnode
	logical, save :: bnodes
	logical, save :: boutput
	logical, save :: bneedbasin
	logical, save :: btrans

	logical, save :: bopen

	!logical, save :: btmin
	!logical, save :: btmax
	!logical, save :: binclusive
	!double precision, save :: atmin
	!double precision, save :: atmax

	logical, save :: bthreshold
	double precision, save :: threshold

	integer, save :: nodesp
	integer, save :: nnodes = 0
	integer, save, allocatable :: nodes(:)
	integer, save, allocatable :: nodese(:)

	real, save :: fact			= 1

	integer, save :: istep
	integer, save :: mode
	integer, save :: modeb

	!integer, save :: date = 0
	!integer, save :: time = 0
	!integer, save :: datetime(2) = 0

	integer, save :: catmode = 0

        character*80, save :: infile		= ' '
        character*80, save :: stmin		= ' '
        character*80, save :: stmax		= ' '
        character*80, save :: nodefile		= ' '
        character*10, save :: outformat		= ' '

	integer, save :: layer = 0
	integer, save :: ivar3 = 0
	
	integer, save :: nfile = 0
	character*10, save, allocatable :: file_type(:)

        character*80, save :: shyfilename = ' '
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

	if( binitialized ) return

	if( type == 'SHY' ) then
          call clo_init(program,'shy-file|bas-file|str-file','3.0')
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop plotutil_set_options: unknown type'
	end if

        call clo_add_info('plots a shy file')

        call clo_add_sep('what to do (only one of these may be given)')

	call clo_add_option('2d',.true.,'plot vertical average (default)')
	call clo_add_option('layer l',0,'plot layer l')
	call clo_add_option('varid id',0,'plot variable id')
	call clo_add_option('varname name',' ','plot variable name')

        call clo_add_sep('options in/output')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of values')
        call clo_add_option('quiet',.false.,'do not be verbose')

        call clo_add_sep('additional options')

	call clo_add_option('freq n',0.,'frequency for plot')
        call clo_add_option('tmin time',' '
     +                  ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                  ,'only process up to time')

	end subroutine plotutil_set_options

!************************************************************

	subroutine plotutil_get_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	integer ivar
	character*80 name

	if( binitialized ) return

        call clo_get_option('2d',b2d)
        call clo_get_option('layer',layer)
        call clo_get_option('varid',ivar3)
        call clo_get_option('varname',name)

        call clo_get_option('verb',bverb)
        call clo_get_option('write',bwrite)
        call clo_get_option('quiet',bquiet)

        call clo_get_option('freq',ifreq)
        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)

        if( .not. bask .and. .not. bmem ) call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

        if( .not. bquiet ) then
	  if( type == 'SHY' ) then
            call shyfem_copyright('shyplot - Plot SHY files')
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop plotutil_get_options: unknown type'
	  end if
	end if

        if( ivar3 > 0 .and. name /= ' ' ) then
          write(6,*) 'You cannot give both varid and varname'
          stop 'error stop plotutil_get_options'
        end if

	if( name .ne. ' ' ) call string2ivar(name,ivar3)

	if( ivar3 < 0 ) then
          write(6,*) 'variable name unknown: ',trim(name)
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

!************************************************************
!************************************************************
!************************************************************

c***************************************************************

	subroutine classify_files

	use clo
	use shyfile
	use basin
	use nls

	character*80 file
	integer i
	integer nshy,nbas,nstr,nunk,ngrd

	logical is_grd_file

	nshy = 0
	nbas = 0
	nstr = 0
	nunk = 0
	ngrd = 0

	nfile = clo_number_of_files()
	allocate(file_type(nfile))
	if( nfile == 0 ) return

	file_type = 'unknown'

	do i=1,nfile
	  call clo_get_file(i,file)
	  if( shy_is_shy_file(file) ) then
	    file_type(i) = 'shy'
	    nshy = nshy + 1
	    shyfilename = file
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

	!write(6,*) 'classify_files: ',nshy,nbas,nstr,nunk

	if( nunk > 0 ) stop 'error stop classify_files'
	if( nshy > 1 ) then
	  write(6,*) 'cannot plot more than one SHY file'
	  stop 'error stop classify_files'
	end if
	if( nbas > 0 .and. ngrd > 0 ) then
	  write(6,*) 'both BAS and GRD files given... cannot handle'
	  stop 'error stop classify_files'
	end if
	if( nshy > 0 ) then
	  if( nbas > 0 .or. ngrd > 0 ) then
	    write(6,*) 'basin given but not needed... ignoring'
	  end if
	end if
	if( shyfilename == ' ' .and. basfilename == ' ' ) then
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
	    write(6,*) 'reading str file ',ivar,trim(file)
            iunit = ifileo(0,file,'form','old')
            if( iunit <= 0 ) then
	      write(6,*) 'cannot read file: ',trim(file)
	      stop 'error stop read_str_files'
	    end if
	    call nlsa(iunit,ivar)
	    close(iunit)
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

