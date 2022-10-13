
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2020  Georg Umgiesser
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

! utility routines for shybas: basutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 01.10.2015	ggu	converted to basutil
! 10.10.2015	ggu	changed VERS_7_3_2
! 12.10.2015	ggu	changed VERS_7_3_3
! 19.10.2015	ggu	changed VERS_7_3_6
! 22.03.2016	ggu	changed VERS_7_5_6
! 15.04.2016	ggu	changed VERS_7_5_8
! 01.05.2016	ggu	new routine to convert depth from node to elem
! 25.05.2016	ggu	changed VERS_7_5_10
! 11.10.2016	ggu	changed VERS_7_5_20
! 21.03.2017	ggu	new flag area to compute area/vol on area code
! 31.03.2017	ggu	changed VERS_7_5_24
! 25.05.2017	ggu	changed VERS_7_5_28
! 22.02.2018	ggu	changed VERS_7_5_42
! 13.04.2018	ggu	new routine for partitioning
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 12.02.2020	ggu	new command line option -reg (and breg)
! 01.04.2020	ggu	new option -custom (bcustom)
! 21.05.2020	ggu	better handle copyright notice
! 14.02.2022	ggu	no ike, some extra comments
! 16.02.2022	ggu	new option -boxgrd implemented (bboxgrd, index_file)
! 12.10.2022	ggu	new option -detail (bdetail)
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
	logical, save :: breg
	logical, save :: bunique
	logical, save :: bgr3
	logical, save :: bmsh
	logical, save :: bdelem
	logical, save :: bnpart

	logical, save :: bquality
	logical, save :: bresol
	logical, save :: bfreq
	logical, save :: bcheck
	logical, save :: bcompare
	logical, save :: bbox
	logical, save :: bboxgrd
	logical, save :: barea
	logical, save :: binvert
	logical, save :: bcustom

	real, save :: hsigma
	real, save :: dreg

        character*80, save :: bfile
        character*80, save :: lfile
        character*80, save :: index_file
	logical, save :: ball		!interpolate everywhere
	integer, save :: btype		!only interpolate on elems type=btype
	logical, save :: bnode		!interpolate on nodes, else on elements
	integer, save :: bmode		!type of interpolation
	real, save :: usfact
	real, save :: uxfact

	logical, save :: bsmooth		!internal - should I smooth
	real, save :: hmin
	real, save :: hmax
	real, save :: asmooth
	integer, save :: iter

	logical, save :: bverb
	logical, save :: bquiet
	logical, save :: bsilent
	logical, save :: bnomin
	logical, save :: bdetail

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
        call clo_add_option('silent',.false.,'be silent')
	call clo_add_option('nomin',.false.,'do not compute min distance')
	call clo_add_option('detail',.false.,'write code details')
	call clo_add_option('area',.false.,'area/vol for each area code')

        call clo_add_sep('output options:')

        call clo_add_option('grd',.false.,'writes grd file')
        call clo_add_option('xyz',.false.,'writes xyz file')
        call clo_add_option('depth',.false.,'writes depth values')
        call clo_add_option('reg dxy',0.,'writes regular depth values')
        call clo_add_option('unique',.false.
     +		,'writes grd file with unique depths on nodes')
        call clo_add_option('delem',.false.
     +		,'writes grd file with constant depths on elements')
        call clo_add_option('hsigma',-1,'creates hybrid depth level')
        call clo_add_option('npart',.false.
     +		,'writes grd file with nodal partition to be visualized')
        call clo_add_option('part grd-file',' '
     +		,'uses lines contained in grd-file to create partition')
        call clo_add_option('gr3',.false.
     +		,'writes grid in gr3 format (for WWMIII model')
        call clo_add_option('msh',.false.
     +		,'writes grid in msh (gmsh v. 2) format (for WW3 model')

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
        call clo_add_option('boxgrd index',' ','creates grd from index')
        call clo_add_option('custom',.false.
     +				,'run custom routine defined by user')

        call clo_add_sep('bathymetry interpolation:')

        call clo_add_option('bfile bathy',' '
     +				,'bathymetry file for interpolation')
        call clo_add_option('all',.false.
     +				,'interpolate in all items')
        call clo_add_option('btype type',-1,'interpolate only on '//
     +				 'elems with type (Default -1)')
        call clo_add_option('node',.false.,'interpolate on nodes '//
     +				 ',else on elements')
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
        call clo_get_option('reg',dreg)
        call clo_get_option('unique',bunique)
        call clo_get_option('delem',bdelem)
        call clo_get_option('npart',bnpart)
        call clo_get_option('part',lfile)
        call clo_get_option('gr3',bgr3)
        call clo_get_option('msh',bmsh)

        call clo_get_option('quality',bquality)
        call clo_get_option('resol',bresol)
        call clo_get_option('freq',bfreq)
        call clo_get_option('check',bcheck)
        call clo_get_option('compare',bcompare)
        call clo_get_option('invert_depth',binvert)
        call clo_get_option('box',bbox)
        call clo_get_option('boxgrd',index_file)
        call clo_get_option('custom',bcustom)

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
        call clo_get_option('silent',bsilent)
        call clo_get_option('nomin',bnomin)
        call clo_get_option('detail',bdetail)
        call clo_get_option('area',barea)

	if( bsilent ) bquiet = .true.

        if( .not. bquiet ) then
	  if( type == 'BAS' .or. type == 'GRD' ) then
	    call shyfem_copyright('shybas - elaborates GRD/BAS files')
	  else
	    write(6,*) 'type : ',trim(type)
	    stop 'error stop basutil_get_options: unknown type'
	  end if
        end if

        call clo_check_files(1)
        call clo_get_file(1,infile)
        call ap_set_names(' ',infile)

	breg = ( dreg > 0. )
	bboxgrd = ( index_file /= ' ' )

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

        subroutine transfer_depth(bnode)

c copies depth values from elems/nodes to nodes/elems

        use mod_depth
        use basin

        implicit none

	logical bnode	!.true.: depth is in hkv; .false.: depth is in hev

        integer k,ie,ii
        real depth
        integer ic(nkn)

	if( bnode ) then
          do ie=1,nel
            depth = 0.
            do ii=1,3
              k = nen3v(ii,ie)
              depth = depth + hkv(k)
	      hm3v(ii,ie) = hkv(k)
            end do
            hev(ie) = depth / 3.
          end do
	else
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
        end if

        end

c*******************************************************************

