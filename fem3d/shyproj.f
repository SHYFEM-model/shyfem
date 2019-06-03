
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

! revision log :
!
! 24.01.2011	ggu	written from scratch
! 18.11.2011	ggu	adapted to new function call
! 16.03.2012	ggu	writes also transformed lines
! 25.05.2017	ggu	changed VERS_7_5_28
! 10.07.2017	ggu	new projection LCC
! 09.10.2017	ggu	changed VERS_7_5_33
! 08.11.2017	ggu	save nodal depth for nodes without element
! 22.02.2018	ggu	changed VERS_7_5_42
! 16.02.2019	ggu	changed VERS_7_5_60
!
!****************************************************************

        program basproj

! projection of basin

        use basin
        use projection
        use clo

        implicit none

        character(50) :: gfile,nfile
        character(80) :: title

        logical :: berror
        integer :: nk,ne,nl,nne,nnl
        integer :: mode,iproj,isphe
	real, allocatable :: haux(:)		!node depth without element
        double precision, dimension(9) :: c_param

!---------------------------------------------------------------
! open grd file
!---------------------------------------------------------------

	call shyfem_copyright('shyproj - projections for FEM grid')

        call clo_init('shyproj','grd-file','1.0')
        call clo_add_info('converts grd-files between lat/lon and cart')

        call clo_add_sep('general options')

        call clo_add_option('proj proj',' ','type of projection to use')
        call clo_add_option('param list',' '
     +				,'parameters for projection')

        call clo_add_sep('additional information')
        call clo_add_com('  proj can be one of the following:')
        call clo_add_com('    GB     Gauss Boaga')
        call clo_add_com('    UTM    UTM')
        call clo_add_com('    EC     equidistant cylindrical')
        call clo_add_com('    LCC    Lambert conformal conic')
        call clo_add_com('  list gives the parameters needed')
        call clo_add_com('    GB     fuse[,x-shift,y-shift]')
        call clo_add_com('    UTM    sector'//
     +		'[,false-easting,false-northing[,scale-factor]]')
        call clo_add_com('    EC     central-lat,lon-orig[,lat-orig]')
        call clo_add_com('    LCC    lon0,lat0,lat1,lat2'//
     +		'[,false-easting,false-northing]')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,gfile)

        call grd_read(gfile)

        call grd_get_params(nk,ne,nl,nne,nnl)
        write(6,*) 'grid info: ',nk,ne,nl

        if( nk == 0 ) then
            write(6,*) 'nk: ',nk
            stop 'error stop vp: no nodes or elements in basin'
        end if

	allocate(haux(nk))
	call grd_get_nodal_depth(haux)

        call grd_to_basin

        call check_spheric_ev
        call get_coords_ev(isphe)

        mode = 1		!+1: cart to geo  -1: geo to cart
        if( isphe == 1 ) mode = -1
        write(6,*) 'isphe,mode: ',isphe,mode
        if( mode == 1 ) then
            write(6,*) 'converting from cartesian to geographical'
        else
            write(6,*) 'converting from geographical to cartesian'
        end if

	call set_projection(iproj,c_param)

!---------------------------------------------------------------
! parameters for projection
!---------------------------------------------------------------

!	1	Gauss-Boaga
!	2	UTM
!	3	equidistant cylindrical
!	4	UTM non standard
!	5	Lambert conformal conic

!---------------------------------------------------------------

! Mediterranean

!	iproj = 3	     	     !equidistant cylindrical
!        c_param(1) = 38.             !central latitude (phi)
!        c_param(2) = 15.             !longitude of origin (lon0)
!        c_param(3) = 38.             !latitude of origin (lat0)

! Nador

!	iproj = 3	     	     !equidistant cylindrical
!        c_param(1) = 35.             !central latitude (phi)
!        c_param(2) = -3.             !longitude of origin (lon0)
!        c_param(3) = 35.             !latitude of origin (lat0)

! Black Sea

!	iproj = 3		     !equidistant cylindrical
!        c_param(1) = 43.5            !central latitude (phi)
!        c_param(2) = 34.             !longitude of origin (lon0)
!        c_param(3) = 43.5            !latitude of origin (lat0)

! Klaipeda

!	iproj = 4		     !UTM Lithuania
!        c_param(1) = 24.             !longitude of origin (lon0)
!        c_param(2) = -500000.        !false easting
!        c_param(3) = 0.              !false northing
!        c_param(2) = -220000.        !false easting
!        c_param(3) = +6100000.       !false northing
!        c_param(4) = 0.9998          !scale factor

! ??

!        iproj = 2		     !UTM
!        c_param(1) = 33.             !zone
!        c_param(2) = -500000.        !false easting
!        c_param(3) = 0.              !false northing

! Laguna di Venezia

!	iproj = 1		     !Gauss-Boaga
!        c_param(1) = 2.              !fuse
!        c_param(2) = 2280000.        !shift in x
!        c_param(3) = 5000000.        !shift in y

! Laguna di Marano-Grado

!	iproj = 1		     !Gauss-Boaga
!        c_param(1) = 2.              !fuse
!        c_param(2) = 0.              !shift in x
!        c_param(3) = 0.              !shift in y

! Turkey lake for Ali

!	iproj = 2		     !UTM
!        c_param(1) = 36.             !zone
!        c_param(2) = -500000.        !false easting
!        c_param(3) = 0.              !false northing
!        c_param(2) = -0.13E+06       !false easting
!        c_param(3) = 0.418E+07       !false northing

!---------------------------------------------------------------
! do projection
!---------------------------------------------------------------

        call init_coords(iproj,c_param)
        call convert_coords(mode,nkn,xgv,ygv,xgv,ygv)	!overwrite coords

!---------------------------------------------------------------
! write new grd file
!---------------------------------------------------------------

        nfile = 'proj.grd'

	call grd_set_coords(nkn,xgv,ygv)
	call grd_set_unused_node_depth(haux)

        call grd_write(nfile)

        write(6,*) 'file has been written to ',nfile

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

        stop
        end program

!***************************************************************

	subroutine set_projection(iproj,c_param)

	use clo

	implicit none

	integer iproj
	double precision c_param(9)

	integer is
	character*80 proj,param

	integer iscand

	call clo_get_option('proj',proj)
	call clo_get_option('param',param)

	c_param = 0
	is = iscand(param,c_param,9)

	if( proj == 'GB' .or. proj == 'gb' ) then
	  iproj = 1
	  if( is < 1 ) goto 99
	else if( proj == 'UTM' .or. proj == 'utm' ) then
	  iproj = 2
	  if( is == 4 ) iproj = 4
	  if( is < 1 ) goto 99
	else if( proj == 'EC' .or. proj == 'ec' ) then
	  iproj = 3
	  if( is < 2 ) goto 99
	  if( is == 2 ) c_param(3) = c_param(1)
	else if( proj == 'LCC' .or. proj == 'lcc' ) then
	  iproj = 5
	  if( is < 4 ) goto 99
	else if( proj == ' ' ) then
	  write(6,*) 'you must give at least one projection'
	  stop 'error stop set_projection: unknown projection'
	else
	  write(6,*) 'unknown projection: ',trim(proj)
	  stop 'error stop set_projection: unknown projection'
	end if

	return
   99	continue
	write(6,*) 'not enough parameters for projection ',trim(proj)
	stop 'error stop set_projection: missing parameters'
	end

!***************************************************************

        subroutine grd_write_params

	use grd

	implicit none

	integer nk,ne,nl,nne,nnl

        call grd_get_params(nk,ne,nl,nne,nnl)

	write(6,*) nk,ne,nl,nne,nnl

	end

!***************************************************************

