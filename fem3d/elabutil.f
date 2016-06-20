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
	module elabutil
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
	logical, save :: bdate

	integer, save :: ifreq
	integer, save :: tmin
	integer, save :: tmax

	logical, save :: bnode
	logical, save :: bnodes
	logical, save :: boutput
	logical, save :: bneedbasin
	logical, save :: btrans

	logical, save :: bopen

	logical, save :: binclusive

	logical, save :: bdiff
	real, save :: deps

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

	integer, save :: catmode = 0

        character*80, save :: infile		= ' '
        character*80, save :: stmin		= ' '
        character*80, save :: stmax		= ' '
        character*80, save :: nodefile		= ' '
        character*80, save :: regstring		= ' '
        character*10, save :: outformat		= ' '

	logical, save :: bcompat = .true.	!compatibility output

!====================================================
	contains
!====================================================

	subroutine elabutil_init(type,what)

	use clo

	character*(*) type
	character*(*), optional :: what

	character*80 program

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

	if( binitialized ) return

        call get_shyfem_version(vers)
	version = '4.0' // ' (SHYFEM version ' // trim(vers) // ')'
	!write(6,*) 'ggguuu version: ',trim(version)

	if( type == 'SHY' ) then
          call clo_init(program,'shy-file',version)
	else if( type == 'NOS' ) then
          call clo_init(program,'nos-file',version)
	else if( type == 'OUS' ) then
          call clo_init(program,'ous-file',version)
	else if( type == 'EXT' ) then
          call clo_init(program,'ext-file',version)
	else if( type == 'FLX' ) then
          call clo_init(program,'flx-file',version)
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop elabutil_set_options: unknown type'
	end if

        call clo_add_info('returns info on or elaborates a shy file')

        call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('out',.false.,'writes new shy file')
        call clo_add_option('averbas',.false.,'average over basin')
        call clo_add_option('aver',.false.,'average over records')
        call clo_add_option('averdir',.false.,'average for directions')
        call clo_add_option('sum',.false.,'sum over records')
        call clo_add_option('min',.false.,'minimum of records')
        call clo_add_option('max',.false.,'maximum of records')
	call clo_add_option('std',.false.,'standard deviation of records')
        call clo_add_option('rms',.false.,'root mean square of records')
        call clo_add_option('sumvar',.false.,'sum over variables')
        call clo_add_option('split',.false.,'split file for variables')
	call clo_add_option('2d',.false.,'average vertically to 2d field')

	call clo_add_option('threshold t',flag,'records over threshold t')
	call clo_add_option('fact fact',1.,'multiply values by fact')
	call clo_add_option('diff',.false.
     +			,'check if 2 files are different')
	call clo_add_option('diffeps deps',0.
     +			,'files differ by more than deps (default 0)')

        call clo_add_sep('options in/output')

        !call clo_add_option('basin name',' ','name of basin to be used')
	!call clo_add_option('mem',.false.,'if no file given use memory')
        !call clo_add_option('ask',.false.,'ask for simulation')
        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of values')
        call clo_add_option('quiet',.false.,'do not be verbose')

        call clo_add_sep('additional options')

        call clo_add_option('node n',0,'extract vars of node number n')
        call clo_add_option('nodes nfile',' '
     +			,'extract vars at nodes given in file nfile')
	call clo_add_option('freq n',0.
     +			,'frequency for aver/sum/min/max/std/rms')
        call clo_add_option('tmin time',' '
     +                  ,'only process starting from time')
        call clo_add_option('tmax time',' '
     +                  ,'only process up to time')
        call clo_add_option('inclusive',.false.
     +			,'output includes whole time period given')

        call clo_add_option('outformat form','native','output format')
        call clo_add_option('catmode cmode',0.,'concatenation mode')
        call clo_add_option('reg rstring',' ','regular interpolation')

        call clo_add_option('area grd-file',' '
     +			,'line delimiting area for -averbas option')
	call clo_hide_option('area')

	call clo_add_sep('additional information')
	call clo_add_com('  nfile is file with nodes to extract')
	call clo_add_com('  time is either integer for relative time or')
	call clo_add_com('    format is YYYY-MM-DD[::hh[:mm[:ss]]]')
	call clo_add_com('  possible output formats are: shy|gis|fem|nc'
     +				// ' (Default shy)')
	call clo_add_com('  possible values for cmode are: -1,0,+1'
     +				// ' (Default 0)')
	call clo_add_com('    -1: all of first file, '
     +				// 'then remaining of second')
	call clo_add_com('     0: simply concatenate files')
	call clo_add_com('    +1: first file until start of second, '
     +				// 'then all of second')
	call clo_add_com('  rstring is: dx[,dy[,x0,y0,x1,y1]]')
	call clo_add_com('    if only dx is given -> dy=dx')
	call clo_add_com('    if only dx,dy are given -> bounds computed')
	call clo_add_com('  -diff needs two files, exits at difference')
	call clo_add_com('    with -out writes difference to out file')

	end subroutine elabutil_set_options

!************************************************************

	subroutine elabutil_get_options(type,program)

	use clo

	character*(*) type
	character*(*) program

	if( binitialized ) return

        call clo_get_option('out',bout)
        call clo_get_option('averbas',baverbas)
        call clo_get_option('aver',baver)
        call clo_get_option('averdir',baverdir)
        call clo_get_option('sum',bsum)
        call clo_get_option('min',bmin)
        call clo_get_option('max',bmax)
        call clo_get_option('std',bstd)
        call clo_get_option('rms',brms)
        call clo_get_option('sumvar',bsumvar)
        call clo_get_option('split',bsplit)
        call clo_get_option('2d',b2d)

        call clo_get_option('threshold',threshold)
        call clo_get_option('fact',fact)
        call clo_get_option('diff',bdiff)
        call clo_get_option('diffeps',deps)

        call clo_get_option('node',nodesp)
        call clo_get_option('nodes',nodefile)

        !call clo_get_option('mem',bmem)
        !call clo_get_option('ask',bask)
        call clo_get_option('verb',bverb)
        call clo_get_option('write',bwrite)
        call clo_get_option('quiet',bquiet)

        call clo_get_option('freq',ifreq)
        call clo_get_option('tmin',stmin)
        call clo_get_option('tmax',stmax)
        call clo_get_option('inclusive',binclusive)

        call clo_get_option('outformat',outformat)
        call clo_get_option('catmode',catmode)
        call clo_get_option('reg',regstring)

        if( .not. bask .and. .not. bmem ) call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

        if( .not. bquiet ) then
	  if( type == 'SHY' ) then
            call shyfem_copyright('shyelab - Elaborate SHY files')
	  else if( type == 'NOS' ) then
            call shyfem_copyright('noselab - Elaborate NOS files')
	  else if( type == 'OUS' ) then
            call shyfem_copyright('ouselab - Elaborate OUS files')
	  else if( type == 'EXT' ) then
            call shyfem_copyright('extelab - Elaborate EXT files')
	  else if( type == 'FLX' ) then
            call shyfem_copyright('flxelab - Elaborate FLX files')
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop elabutil_get_options: unknown type'
	  end if
        end if

        bnode = nodesp > 0
        bnodes = nodefile .ne. ' '

        boutput = bout .or. b2d
	boutput = boutput .or. outformat /= 'native'
        !btrans is added later
	if( bsumvar ) boutput = .false.

        bneedbasin = b2d .or. baverbas .or. bnode .or. bnodes
	bneedbasin = bneedbasin .or. outformat == 'gis'
	bneedbasin = bneedbasin .or. ( type == 'OUS' .and. bsplit )

        modeb = 2
        if( bneedbasin ) modeb = 3

	bthreshold = ( threshold /= flag )

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

        mode = 0
        btrans = .false.

	ic = count( (/baver,bsum,bmin,bmax,bstd,brms
     +				,bthreshold,baverdir/) )

	if( ic > 1 ) then
	  write(6,*) 'Only one of the following options can be given:'
	  write(6,*) '-aver -sum -min -max -std -rms'
	  write(6,*) '-threshold -averdir'
	  stop 'error stop elabutil_set_averaging: incompatible options'
	end if

        if( baver ) mode = 1
        if( bsum )  mode = 2
        if( bmin )  mode = 3
        if( bmax )  mode = 4
        if( bstd )  mode = 5
        if( brms )  mode = 6
        if( bthreshold )  mode = 7
        if( baverdir ) mode = 8

        if( mode > 0 ) then             !prepare for averaging
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
        do l=1,nlv
          if( count(l) > 0 ) then
            write(6,*) l,count(l),ccount(l)
            nc = nc + count(l)
          end if
        end do
        write(6,*) 'total count: ',nc

        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine shy_write_aver(dtime,ivar,cmin,cmax,cmed,vtot)

c writes basin average to file

        implicit none

        double precision dtime
        integer ivar
        real cmin,cmax,cmed,vtot

	integer it
        real totmass

	it = nint(dtime)
        totmass = cmed * vtot

        !write(6,1234) it,ivar,cmin,cmed,cmax,totmass
        write(100+ivar,1235) it,cmin,cmed,cmax,totmass
        write(100,1236) it,vtot

        write(6,2234) dtime,ivar,cmin,cmed,cmax,totmass
        write(200+ivar,2235) dtime,cmin,cmed,cmax,totmass
        write(200,2236) dtime,vtot

	return
 1234   format(i10,i10,3f12.4,e14.6)
 1235   format(i10,3f12.4,e14.6)
 1236   format(i10,e14.6)
 2234   format(f15.2,i10,3f12.4,e14.6)
 2235   format(f15.2,3f12.4,e14.6)
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

