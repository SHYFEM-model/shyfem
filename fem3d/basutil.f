!
! utility routines for shybas: basutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 01.10.2015	ggu	converted to basutil
! 01.05.2016	ggu	new routine to convert depth from node to elem
! 21.03.2017	ggu	new flag area to compute area/vol on area code
! 13.04.2018	ggu	new routine for partitioning
!
!************************************************************

!====================================================
	module basutil
!====================================================

	implicit none

	logical, save, private :: binitialized = .false.
	double precision, parameter :: dflag = -999.

	logical, save :: binter

	logical, save :: bgrd
	logical, save :: bxyz
	logical, save :: bdepth
	logical, save :: bunique
	logical, save :: bdelem
	logical, save :: bnpart

	logical, save :: bquality
	logical, save :: bresol
	logical, save :: bfreq
	logical, save :: bcheck
	logical, save :: bcompare
	logical, save :: bbox
	logical, save :: barea
	logical, save :: binvert

	real, save :: hsigma

        character*80, save :: bfile
        character*80, save :: lfile
	logical, save :: ball
	integer, save :: btype
	logical, save :: bnode
	integer, save :: bmode
	real, save :: usfact
	real, save :: uxfact

	logical, save :: bsmooth		!internal - should I smooth
	real, save :: hmin
	real, save :: hmax
	real, save :: asmooth
	integer, save :: iter

	logical, save :: bverb
	logical, save :: bquiet
	logical, save :: bnomin

        character*80, save :: infile

	logical, save :: breadbas	!internal - true if file read is bas

!====================================================
	contains
!====================================================

	subroutine basutil_init(type)

	use clo

	character*(*) type

	call basutil_set_options(type)
	call clo_parse_options
	call basutil_get_options(type)

	binitialized = .true.

	end subroutine basutil_init

!************************************************************

	subroutine basutil_set_options(type)

	use clo

	character*(*) type

	if( binitialized ) return

	if( type == 'BAS' ) then
          call clo_init('shybas','bas-file','3.0')
	else
	  write(6,*) 'type : ',trim(type)
	  stop 'error stop basutil_set_options: unknown type'
	end if

        call clo_add_info('returns info on or elaborates a bas file')

        !call clo_add_sep('what to do (only one of these may be given)')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_option('nomin',.false.,'do not compute min distance')
	call clo_add_option('area',.false.,'area/vol for each area code')

        call clo_add_sep('output options:')

        call clo_add_option('grd',.false.,'writes grd file')
        call clo_add_option('xyz',.false.,'writes xyz file')
        call clo_add_option('depth',.false.,'writes depth values')
        call clo_add_option('unique',.false.
     +		,'writes grd file with unique depths')
        call clo_add_option('delem',.false.
     +		,'writes grd file with constant depths on elements')
        call clo_add_option('hsigma',-1,'creates hybrid depth level')
        call clo_add_option('npart',.false.
     +		,'writes grd file with nodal partition to be visualized')
        call clo_add_option('part grd-file',' '
     +		,'uses lines contained in grd-file to create partition')

        call clo_add_sep('what to do:')

        call clo_add_option('inter',.false.
     +				,'interactively shows nodes and elems')

        call clo_add_option('quality',.false.,'shows quality of file')
        call clo_add_option('resol',.false.,'writes resolution of file')
        call clo_add_option('freq',.false.,'computes frequency curves')
        call clo_add_option('check',.false.,'runs extra check on file')
        call clo_add_option('compare',.false.
     +				,'compares depth of 2 basins')
        call clo_add_option('invert_depth',.false.
     +				,'inverts depth values')
        call clo_add_option('box',.false.,'creates index for box model')

        call clo_add_sep('bathymetry interpolation:')

        call clo_add_option('bfile bathy',' '
     +				,'bathymetry file for interpolation')
        call clo_add_option('all',.false.
     +				,'interpolate in all elements')
        call clo_add_option('btype type',-1,'interpolate only on '//
     +				 'elems type (Default -1)')
        call clo_add_option('node',.false.,'interpolate to nodes')
        call clo_add_option('bmode mode',1,'mode of interpolation')
        call clo_add_option('usfact fact',1,'factor for std '//
     +				 '(Default 1)')
        call clo_add_option('uxfact fact',3,'factor for max radius '//
     +				 '(Default 3)')

        call clo_add_sep('limiting and smoothing bathymetry:')
        call clo_add_option('hmin val',-99999.,'minimum depth')
        call clo_add_option('hmax val',99999.,'maximum depth')
        call clo_add_option('asmooth alpha',0,'alpha for smoothing')
        call clo_add_option('iter n',0,'iterations for smoothing')

        !call clo_add_sep('additional options')
        call clo_add_sep(' ')
        call clo_add_sep(' interpolation mode: (default=1) ')
        call clo_add_sep('   1 exponential')
        call clo_add_sep('   2 uniform on squares')
        call clo_add_sep('   3 exponential with autocorrelation')

	end subroutine basutil_set_options

!************************************************************

	subroutine basutil_get_options(type)

	use clo

	character*(*) type

	if( binitialized ) return

        call clo_get_option('inter',binter)

        call clo_get_option('grd',bgrd)
        call clo_get_option('xyz',bxyz)
        call clo_get_option('depth',bdepth)
        call clo_get_option('unique',bunique)
        call clo_get_option('delem',bdelem)
        call clo_get_option('npart',bnpart)
        call clo_get_option('part',lfile)

        call clo_get_option('quality',bquality)
        call clo_get_option('resol',bresol)
        call clo_get_option('freq',bfreq)
        call clo_get_option('check',bcheck)
        call clo_get_option('compare',bcompare)
        call clo_get_option('invert_depth',binvert)
        call clo_get_option('box',bbox)

        call clo_get_option('hsigma',hsigma)

        call clo_get_option('bfile',bfile)
        call clo_get_option('all',ball)
        call clo_get_option('btype',btype)
        call clo_get_option('node',bnode)
        call clo_get_option('bmode',bmode)
        call clo_get_option('usfact',usfact)
        call clo_get_option('uxfact',uxfact)

        call clo_get_option('hmin',hmin)
        call clo_get_option('hmax',hmax)
        call clo_get_option('asmooth',asmooth)
        call clo_get_option('iter',iter)

        call clo_get_option('verb',bverb)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nomin',bnomin)
        call clo_get_option('area',barea)

        call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

        if( .not. bquiet ) then
	  if( type == 'BAS' ) then
            call shyfem_copyright('shybas - Elaborate BAS files')
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop basutil_get_options: unknown type'
	  end if
        end if

	bsmooth = .false.
	bsmooth = bsmooth .or. hmin /= -99999.
	bsmooth = bsmooth .or. hmax /=  99999.
	bsmooth = bsmooth .or. asmooth > 0.
	bsmooth = bsmooth .or. iter > 0

	end subroutine basutil_get_options

!====================================================
	end module basutil
!====================================================

c***************************************************************

        subroutine transfer_depth(ike)

c copies depth values from elems/nodes to nodes/elems

        use mod_depth
        use basin

        implicit none

        integer ike	!1: depth is in hev    2: depth is in hkv

        integer k,ie,ii
        real depth
        integer ic(nkn)

        if( ike .eq. 1 ) then           !elementwise

          ic = 0
          hkv = 0.

          do ie=1,nel
            do ii=1,3
              k = nen3v(ii,ie)
              hkv(k) = hkv(k) + hev(ie)
	      hm3v(ii,ie) = hev(ie)
              ic(k) = ic(k) + 1                 !BUG - this was missing
            end do
          end do

	  hkv = hkv / ic

        else                            !nodewise

          do ie=1,nel
            depth = 0.
            do ii=1,3
              k = nen3v(ii,ie)
              depth = depth + hkv(k)
	      hm3v(ii,ie) = hkv(k)
            end do
            hev(ie) = depth / 3.
          end do

        end if

        end

c*******************************************************************

