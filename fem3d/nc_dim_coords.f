
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2020  Georg Umgiesser
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

! convert nc files to fem files: dimension and coordinates
!
! contents :
!
!
! revision log :
!
! 25.05.2017	ggu	changed VERS_7_5_28
! 11.07.2017	ggu	changed VERS_7_5_30
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.07.2018	ggu	revision control introduced
! 13.07.2018	ggu	changed VERS_7_4_1
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 08.01.2020	ggu	new values for time description
! 21.10.2021	ggu	new dims and coords for vertical
! 27.01.2022	ggu	new values for atmos
! 20.06.2022	ggu	look in variable name to find coordinates
! 23.10.2022	ggu	new routine handle_exceptions()
!
! notes :
!
! to adapt dimensions, coordinates, variables, see routines at end of file
!
!*****************************************************************
!*****************************************************************
!*****************************************************************

!================================================================
        module ncnames
!================================================================

        implicit none

        type, private :: entry

          character*5  :: what
          character*80 :: descrp
          character*10 :: short
          logical :: bclip
          integer :: ilen

        end type entry

        integer, save, private :: idlast = 0
        integer, save, private :: ndim = 0
        type(entry), save, private, allocatable :: pentry(:)

	logical, save :: binitialized = .false.
	integer, save :: idims(2,0:4)
	integer, save :: icoords(2,0:4)
	character*80, save :: cdims(0:4)
	character*80, save :: ccoords(0:4)

	character*5, parameter :: what = 'txyzi'
	integer, parameter :: nwhere = 3
	character*13, parameter :: where(nwhere) = (/
     +					 'standard_name'
     +					,'long_name    '
     +					,'description  '
     +					/)

!================================================================
        contains
!================================================================

        subroutine ncnames_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
        else
          allocate(paux(2*ndim))
          paux(1:ndim) = pentry(1:ndim)
          call move_alloc(paux,pentry)
          ndim = ndim*2
        end if

        end subroutine ncnames_init_alloc

!******************************************************************

        subroutine ncnames_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call ncnames_init_alloc
        end if
	if( idlast > ndim ) stop 'error stop: internal error (3)'
        id = idlast

        call ncnames_init_id(id)

        end subroutine ncnames_init_new_id

!******************************************************************

        subroutine ncnames_init_id(id)

        integer id

        if( id < 1 .or. id > idlast ) then
          stop 'error stop ncnames_init_id: ndim'
        end if

        pentry(id)%what = ' '
        pentry(id)%descrp = ' '
        pentry(id)%short = ' '
        pentry(id)%bclip = .false.
        pentry(id)%ilen = 0

        end subroutine ncnames_init_id

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine ncnames_add_dim(short,descrp,bclip)

	implicit none

	character*(*) descrp,short
	logical, optional :: bclip

	logical bc

	bc = .false.
	if( present(bclip) ) bc = bclip

	call ncnames_add('dim',descrp,short,.false.)

	end subroutine ncnames_add_dim

c*****************************************************************

	subroutine ncnames_add_coord(short,descrp,bclip)

	implicit none

	character*(*) descrp,short
	logical, optional :: bclip

	logical bc

	bc = .false.
	if( present(bclip) ) bc = bclip

	call ncnames_add('coord',descrp,short,bc)

	end subroutine ncnames_add_coord

c*****************************************************************

	subroutine ncnames_add_var(short,descrp,bclip)

	implicit none

	character*(*) descrp,short
	logical, optional :: bclip

	logical bc

	bc = .false.
	if( present(bclip) ) bc = bclip

	call ncnames_add('var',descrp,short,bc)

	end subroutine ncnames_add_var

c*****************************************************************

	subroutine ncnames_add(what,descrp,short,bclip)

! adds item into internal list
!
! what: dim, coord, var
! descrp: description that can be found in nc file
! short: short version for item

	implicit none

	character*(*) what,descrp,short
	logical :: bclip

	logical, parameter :: bdebug = .false.
	integer id,i
	character*80 name

        call ncnames_init_new_id(id)

	name = descrp
	call to_lower(name)

        pentry(id)%what = what
        pentry(id)%descrp = name
        pentry(id)%short = short
        pentry(id)%bclip = bclip
        pentry(id)%ilen = len_trim(name)

	if( bdebug ) then		!GGU
	  i = len_trim(name)
	  write(6,*) 'add ',id,trim(what),' ',trim(name),' ',trim(short)
     +			,' ',bclip,len_trim(name)
     +			,ichar(name(i:i))
	end if

	end subroutine ncnames_add

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine ncnames_get(what,descrp,short)

! looks for item with description descrp and returns short if found
!
! what: dim, coord, var
! descrp: description that can be found in nc file
! short: short version for item

	implicit none

	character*(*) what,descrp,short

	logical bclip,bdebug
	integer id,il,ilen,i,iu
	character*80 name

	bdebug = .true.			!GGU
	bdebug = .false.			!GGU
	bdebug = bdebug .and. what == 'var'

	short = ' '
	name = descrp
	iu = index(name,'[')
	if( iu > 0 ) name = name(1:iu-1)	!delete unit given in []
	call to_lower(name)
	il = len_trim(name)
	if( bdebug ) then
	write(6,*) 'in ncnames_get... -------------------'
	write(6,*) 'name: |',trim(name),'|',il
	do i=1,20
	  !write(6,*) i,ichar(descrp(i:i))
	  if( ichar(descrp(i:i)) == 0 ) then
	    write(6,*) 'char 0 in string: ',trim(name),i,il
	  end if
	end do
	end if
	if( name == ' ' ) stop 'empty string'

	if( bdebug ) then
	write(6,*) 'checking: ',trim(what),idlast,'  ',trim(descrp),bclip
	end if

	do id=1,idlast
	  !bdebug = ( bdebug .and. id == 56 )		!GGU
	  if( pentry(id)%what /= what ) cycle
	  ilen = il
	  if( pentry(id)%bclip ) ilen = pentry(id)%ilen
	if( bdebug ) then
	write(6,*) id,il,ilen,pentry(id)%bclip		!GGU
	write(6,*) pentry(id)%descrp(1:ilen),'  ',trim(name)
	end if
	  if( pentry(id)%descrp(1:ilen) == name(1:ilen) ) then
	    short = pentry(id)%short
	if( bdebug ) write(6,*) 'found: ',trim(short)
	    return
	  end if
	end do

	end subroutine

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine ncnames_get_dims(ncid,bverb)

! insert dimensions in internal structure

	implicit none

	integer ncid
	logical bverb

	logical bdebug
	character*80 short
	integer dim_id,ndims,nlen,i
	character*80 name
	character*1 c

	!character*5, parameter :: what = 'txyzi'

	bdebug = .true.
	bdebug = .false.
	idims = 0
	cdims = ' '

        call nc_get_dim_totnum(ncid,ndims)

        do dim_id=1,ndims

          call nc_get_dim_name(ncid,dim_id,name)
          call nc_get_dim_len(ncid,dim_id,nlen)

	  call ncnames_get('dim',name,short)

	  i = index(what,short(1:1)) - 1
	  if( i < 0 ) then
	    if( bdebug ) then
	      write(6,*) '*** cannot identify dimension: ',trim(name)
	    end if
	  else
	    idims(1,i) = dim_id
	    idims(2,i) = nlen
	    cdims(i) = name
	  end if

        end do

	if( .not. bverb ) return

	write(6,*) 'dimensions in ncnames_get_dims:'
	do i=0,3
	  c = what(i+1:i+1)
	  write(6,*) c,'    ',idims(:,i),'  ',trim(cdims(i))
	end do

	end subroutine

c*****************************************************************

	subroutine ncnames_get_coords(ncid,bverb)

! insert coordinates in internal structure

	implicit none

	integer ncid
	logical bverb

	logical bdebug,bexcept
	character*80 short
	integer var_id,nvars,nlen,i,j
	integer ndims,dimids(1)
	character*80 varname,atext
	character*1 c

	bdebug = bverb
	bdebug = .true.
	bdebug = .false.

	icoords = 0
	ccoords = ' '

	if( bdebug ) write(6,*) 'debug: ncnames_get_coords: ',bdebug

        call nc_get_var_totnum(ncid,nvars)

        do var_id=1,nvars

          call nc_get_var_name(ncid,var_id,varname)

	  if( bdebug ) then
	    write(6,*) '-------------------'
	    write(6,*) 'var: ',var_id,trim(varname)
	  end if

	  short = ' '
	  do j=1,nwhere
	    if( bdebug ) write(6,*) 'looking in ',trim(where(j))
            call nc_get_var_attr(ncid,var_id,trim(where(j)),atext)
	    if( atext == ' ' ) cycle
	    call ncnames_get('coord',atext,short)
	    if( bdebug ) write(6,*) 'found ',trim(atext),'  ',trim(short)
	    if( short /= ' ' ) exit
	  end do

!	  -------------------------------------------
!	  if not found yet try to look at var name
!	  -------------------------------------------

	  if( short == ' ' ) then
	    if( bdebug ) write(6,*) 'looking in variable name '
	    call ncnames_get('coord',varname,short)
	    if( bdebug ) write(6,*) 'found ',trim(varname)
     +					,'  ',trim(short)
	  end if
	  if( short == ' ' ) cycle

	  if( bdebug ) write(6,*) '+++ ',trim(varname),'  ',trim(short)

	  call handle_exceptions(ncid,var_id,bexcept)
	  if( bexcept ) cycle

	  i = index(what,short(1:1)) - 1
	  if( i >= 0 ) then
	    if( icoords(1,i) > 0 ) cycle	!do not insert second one
	    ndims = 0
	    call nc_get_var_ndims(ncid,var_id,ndims,dimids)
	    icoords(1,i) = var_id
	    icoords(2,i) = ndims
	    ccoords(i) = varname
	  end if

        end do

	if( bdebug ) then
	  write(6,*) '-------------------'
	end if

	if( .not. bdebug ) return

	write(6,*) 'coordinates:'
	do i=0,3
	  c = what(i+1:i+1)
	  write(6,*) c,'    ',icoords(:,i),'  ',trim(ccoords(i))
	end do

	end subroutine

c*****************************************************************

	subroutine ncnames_get_vars(ncid)

! looks for variable names

	implicit none

	integer ncid

	logical bdebug
	character*80 short
	integer var_id,nvars,nlen,i,j
	integer ndims,dimids(1)
	character*80 name,atext
	character*1 c

	bdebug = .true.

        call nc_get_var_totnum(ncid,nvars)

        do var_id=1,nvars

          call nc_get_var_name(ncid,var_id,name)

	  do j=1,nwhere
            call nc_get_var_attr(ncid,var_id,trim(where(j)),atext)
	    if( atext == ' ' ) cycle
	    call ncnames_get('var',atext,short)
	    if( short /= ' ' ) exit
	  end do
	  if( short == ' ' ) cycle

	  write(6,*) 'var: ',trim(name),'  ',trim(short)

        end do

	if( .not. bdebug ) return

	end subroutine

c*****************************************************************

	subroutine ncnames_get_var(ncid,var,short)

! looks for variable var and returns short if found

	implicit none

	integer ncid
	character*(*) var,short

	logical bdebug
	integer var_id,j,i
	integer ndims,dimids(1)
	character*80 name,atext
	character*1 c

	bdebug = .true.			!GGU
	bdebug = .false.			!GGU

	call nc_get_var_id(ncid,var,var_id)

	do j=1,nwhere
          call nc_get_var_attr(ncid,var_id,trim(where(j)),atext)
	  if( atext == ' ' ) cycle
	  call ncnames_get('var',atext,short)
	  if( short /= ' ' ) exit
	end do

	if( .not. bdebug ) return

! here the same thing, but with debug messages

	write(6,*) 'looking: ',trim(var),var_id
	do j=1,nwhere
          call nc_get_var_attr(ncid,var_id,trim(where(j)),atext)
	  write(6,*) '   ',trim(where(j)),'  ',trim(atext)
	  if( atext == ' ' ) cycle
	  !do i=1,20
	  !  write(6,*) 'atext ',i,ichar(atext(i:i))
	  !end do
	  call ncnames_get('var',atext,short)
	  write(6,*) '  .. ','  |',trim(atext),'|  ',trim(short)
	  if( short /= ' ' ) exit
	end do

	end subroutine

c*****************************************************************

	subroutine handle_exceptions(ncid,var_id,bexcept)

	implicit none

	integer ncid,var_id
	logical bexcept

	character*80 varname,atext

	bexcept = .false.
        call nc_get_var_name(ncid,var_id,varname)
        call nc_get_var_attr(ncid,var_id,'long_name',atext)

	if( varname == 'Z' .and. atext == 'Geopotential' ) then
	  bexcept = .true.
	end if

	if( bexcept ) then
	  write(6,*) '...skipping: ',trim(varname),' ',trim(atext)
	end if

	end

!================================================================
        end module ncnames
!================================================================

c*****************************************************************
c*****************************************************************
c*****************************************************************
c check routines
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_dims_and_coords(ncid,bverb
     +			,nt,nx,ny,nz
     +			,tcoord,xcoord,ycoord,zcoord)

! returns dimensions and coordinate names

	use ncnames

	implicit none

	integer ncid
	logical bverb
	integer nt,nx,ny,nz
	character(*) tcoord,xcoord,ycoord,zcoord

	logical berror
	integer i
	character*4, save :: string = 'txyz'
	character*1 :: w

        call nc_init_dims_and_coords(ncid,bverb)

	nt = idims(2,0)
	nx = idims(2,1)
	ny = idims(2,2)
	nz = idims(2,3)

	tcoord = ccoords(0)
	xcoord = ccoords(1)
	ycoord = ccoords(2)
	zcoord = ccoords(3)

	berror = .false.

        if( nt > 0 .and. tcoord == ' ' ) then
	  berror = .true.
	  write(6,*) '*** t dimension without variable name'
        end if
        if( nx > 0 .and. xcoord == ' ' ) then
	  berror = .true.
	  write(6,*) '*** x dimension without variable name'
        end if
        if( ny > 0 .and. ycoord == ' ' ) then
	  berror = .true.
	  write(6,*) '*** y dimension without variable name'
        end if
        if( nz > 0 .and. zcoord == ' ' ) then
	  berror = .true.
	  write(6,*) '*** z dimension without variable name'
        end if

	if( bverb .or. berror ) then
	  do i=0,3
	    w = string(i+1:i+1)
	    write(6,*) i,' ',w,'  n = ',idims(2,i)
     +				,' s = ',trim(ccoords(i))
	  end do
	end if

	if( berror ) then
          stop 'error stop: dimension(s) without variable name'
	end if

	end

c*****************************************************************

	subroutine nc_init_dims_and_coords(ncid,bverb)

! tested compatibility - eliminated
! now used to write out dims and coords found (with -verbose)

	use ncnames

	implicit none

	integer ncid
	logical bverb

	logical bextra
	integer nt,nx,ny,nz,i
	character*80 tcoord,xname,yname,zcoord
	character*80 time_d,time_v
	character*1 c

	bextra = .true.
	bextra = .false.

	if( bverb ) write(6,*) 'initializing dims and coords:'

	call ncnames_init

	call ncnames_get_dims(ncid,bextra)
	call ncnames_get_coords(ncid,bextra)

	if( bverb ) then
	  write(6,*) 'dimensions:'
	  write(6,*) 'direction    id   dimension   name'
	  do i=0,3
	    c = what(i+1:i+1)
	    write(6,*) c,'  ',idims(:,i),'  ',trim(cdims(i))
	  end do
	end if

	if( bverb ) then
	  write(6,*) 'coordinates:'
	  write(6,*) 'direction    id   dimension   name'
	  do i=0,3
	    c = what(i+1:i+1)
	    write(6,*) c,'  ',icoords(:,i),'  ',trim(ccoords(i))
	  end do
	end if

	time_d = cdims(0)
	time_v = ccoords(0)
        call nc_set_time_name(time_d,time_v)
	if( bverb ) then
	  write(6,*) 'time:'
	  write(6,*) 'time dimension: ',trim(time_d)
	  write(6,*) 'time variable : ',trim(time_v)
	end if

	if( bverb ) write(6,*) 'initialization successfully ended'

	return
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c test routines
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_dim_coords_test_open(ncid)

	implicit none

	integer ncid,nc
	character*80 file

	nc = command_argument_count()
	if( nc == 0 ) stop 'no file given'
	call get_command_argument(1,file)

        call nc_open_read(ncid,file)

	end

c*****************************************************************

	subroutine nc_dim_coords_test

	use ncnames

	implicit none

	integer ncid
	logical bverb

	bverb = .true.

	call nc_dim_coords_test_open(ncid)

	call ncnames_init

	call ncnames_get_dims(ncid,bverb)
	call ncnames_get_coords(ncid,bverb)

	end

c*****************************************************************

	subroutine nc_dim_compatibility_test

	use ncnames

	implicit none

	integer ncid

	call nc_dim_coords_test_open(ncid)

	!call ncnames_init
	call nc_init_dims_and_coords(ncid,.true.)

	end

c*****************************************************************

!	program nc_dim_coords_main
!	!call nc_dim_coords_test
!	call nc_dim_compatibility_test
!	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c initialization routines - please adapt to your needs
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine ncnames_init

	use ncnames

	implicit none

	character*80 time_d,time_v

	if( binitialized ) return

	cdims = ' '
	ccoords = ' '
	idims = 0
	icoords = 0

	call ncnames_add_dimensions
	call ncnames_add_coordinates
	call ncnames_add_variables

	time_d = cdims(0)
	time_v = ccoords(0)
        call nc_set_time_name(time_d,time_v)

	binitialized = .true.

	end subroutine ncnames_init

c*****************************************************************

	subroutine ncnames_add_dimensions

	use ncnames

	implicit none

	call ncnames_add_dim('t','time')
	call ncnames_add_dim('t','ntime')
	call ncnames_add_dim('t','Time')
	call ncnames_add_dim('t','ocean_time')
	call ncnames_add_dim('t','time_counter')

	call ncnames_add_dim('x','x')
	call ncnames_add_dim('x','xpos')
	call ncnames_add_dim('x','lon')
	call ncnames_add_dim('x','longitude')
	call ncnames_add_dim('x','west_east')
	call ncnames_add_dim('x','xt_ocean')
	call ncnames_add_dim('x','rlon')
	call ncnames_add_dim('x','xi_rho')

	call ncnames_add_dim('y','y')
	call ncnames_add_dim('y','ypos')
	call ncnames_add_dim('y','lat')
	call ncnames_add_dim('y','latitude')
	call ncnames_add_dim('y','south_north')
	call ncnames_add_dim('y','yt_ocean')
	call ncnames_add_dim('y','rlat')
	call ncnames_add_dim('y','eta_rho')

	call ncnames_add_dim('z','z')
	call ncnames_add_dim('z','zpos')
	call ncnames_add_dim('z','bottom_top_stag')
	call ncnames_add_dim('z','level')
	call ncnames_add_dim('z','depth')
	call ncnames_add_dim('z','deptht')
	call ncnames_add_dim('z','height')
	call ncnames_add_dim('z','st_ocean')
	call ncnames_add_dim('z','s_rho')

	call ncnames_add_dim('ignore','crsdim')
	call ncnames_add_dim('ignore','node')
	call ncnames_add_dim('ignore','element')
	call ncnames_add_dim('ignore','vertex')

	end subroutine ncnames_add_dimensions

c*****************************************************************

	subroutine ncnames_add_coordinates

	use ncnames

	implicit none

	logical, parameter :: bclip = .true.

	call ncnames_add_coord('t','Time')
	call ncnames_add_coord('t','time')
	call ncnames_add_coord('t','time_counter')
	call ncnames_add_coord('t','ocean_time')
	call ncnames_add_coord('t','averaged time since initialization')
	call ncnames_add_coord('t','Julian day (UTC) of the station')
	call ncnames_add_coord('t','minutes since',bclip)
	call ncnames_add_coord('t','days since',bclip)

	call ncnames_add_coord('x','lon')
	call ncnames_add_coord('x','longitude')
	call ncnames_add_coord('x','Longitude')
	call ncnames_add_coord('x','Longitude of scalars')
	call ncnames_add_coord('x','LONGITUDE',bclip)

	call ncnames_add_coord('y','lat')
	call ncnames_add_coord('y','latitude')
	call ncnames_add_coord('y','Latitude')
	call ncnames_add_coord('y','Latitude of scalars')
	call ncnames_add_coord('y','LATITUDE',bclip)

	call ncnames_add_coord('z','depth')
	call ncnames_add_coord('z','deptht')
	call ncnames_add_coord('z','zcoord')
	call ncnames_add_coord('z','height')
	call ncnames_add_coord('z','sigma of cell face')
	call ncnames_add_coord('z','bottom of vertical layers')
	call ncnames_add_coord('z','eta values on full',bclip)
	call ncnames_add_coord('z','tcell zstar depth')
	call ncnames_add_coord('z','ocean_s_coordinate_g1')
	call ncnames_add_coord('z','Vertical T levels')
	!call ncnames_add_coord('z','S-coordinate at RHO-points')

	end subroutine ncnames_add_coordinates

c*****************************************************************

	subroutine ncnames_add_variables

! please use shortname as first argument and the description
! of the variable as found in
!     		 'standard_name'  'long_name'  'description'

	use ncnames

	implicit none

	logical, parameter :: bclip = .true.

	call ncnames_add_var('bathy','Surface topography')
	call ncnames_add_var('bathy','surface_altitude')
	call ncnames_add_var('bathy','bathymetry')
	call ncnames_add_var('bathy','sea_floor_depth_below_sea_surface')
	call ncnames_add_var('bathy','depth under water positive')
	call ncnames_add_var('bathy','depth')
	call ncnames_add_var('salt','sea_water_salinity')
	call ncnames_add_var('salt','Salinity')
	call ncnames_add_var('salt','time-averaged salinity')
	call ncnames_add_var('temp','sea_water_potential_temperature')
	call ncnames_add_var('temp','temperature')
	call ncnames_add_var('temp','Conservative temperature')
	call ncnames_add_var('temp','time-averaged potential temperature')
	call ncnames_add_var('zeta','sea_surface_elevation')
	call ncnames_add_var('zeta','Sea Surface height')
	call ncnames_add_var('zeta'
     +			,'water_surface_height_above_geoid')
	call ncnames_add_var('zeta'
     +			,'water_surface_height_above_reference_datum')
	call ncnames_add_var('zeta','surface height on T cells')
	call ncnames_add_var('vel','zonal velocity')
	call ncnames_add_var('vel','eastward_sea_water_velocity',bclip)
	call ncnames_add_var('vel','meridional velocity')
	call ncnames_add_var('vel','northward_sea_water_velocity',bclip)
	call ncnames_add_var('vel','Zonal current speed component')
	call ncnames_add_var('vel','Meridional current speed component')

	call ncnames_add_var('airp','Pressure at the Surface')
	call ncnames_add_var('airp','surface_air_pressure')
	call ncnames_add_var('airp','Mean sea level pressure')
	call ncnames_add_var('airp','SFC PRESSURE')
	call ncnames_add_var('airp','Pressure reduced to MSL')
	call ncnames_add_var('airp','air_pressure')
	call ncnames_add_var('airp','air_pressure_at_sea_level')
	call ncnames_add_var('airp','Sea Level Pressure')
	call ncnames_add_var('wind','eastward_wind')
	call ncnames_add_var('wind','northward_wind_at_10m')
	call ncnames_add_var('wind','eastward_wind_at_10m')
	call ncnames_add_var('wind','northward_wind')
	call ncnames_add_var('wind','U at 10 M')
	call ncnames_add_var('wind','V at 10 M')
	call ncnames_add_var('wind','10 metre U wind component')
	call ncnames_add_var('wind','10 metre V wind component')
	call ncnames_add_var('wind','u-component of wind')
	call ncnames_add_var('wind','v-component of wind')
	call ncnames_add_var('rhum','Relative Humidity at 2 m')
	call ncnames_add_var('rhum','Relative Humidity')
	call ncnames_add_var('rhum','relative_humidity')
	call ncnames_add_var('rhum','relative_humidity_at_2m')
	call ncnames_add_var('shum','specific humidity')
	call ncnames_add_var('mixrat','Water vapor mixing ratio')
	call ncnames_add_var('airt','Temperature at 2 m')
	call ncnames_add_var('airt','TEMP at 2 M')
	call ncnames_add_var('airt','air_temperature')
	call ncnames_add_var('airt','air_temperature_at_2m')
	call ncnames_add_var('cc','total cloud cover')
	call ncnames_add_var('cc','total_cloud_cover')
	call ncnames_add_var('cc','Cloud cover')
	call ncnames_add_var('cc','total cloud fraction')
	call ncnames_add_var('cc','cloud_area_fraction')
	call ncnames_add_var('cc','CLOUD FRACTION')
	call ncnames_add_var('srad','surface_downwelling_shortwave_flux')
	call ncnames_add_var('srad','Short wave flux')
	call ncnames_add_var('srad','DOWNWARD SHORT WAVE FLUX',bclip)
	call ncnames_add_var('srad'
     +			,'DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE')
	call ncnames_add_var('srad'
     +			,'surface_downwelling_shortwave_flux_in_air')
	call ncnames_add_var('srad'
     +			,'surface_solar_radiation_downwards')
	call ncnames_add_var('rain','large_scale_precipitation_amount')
	call ncnames_add_var('rain','accumulated_precipitation_amount')
	call ncnames_add_var('rain'
     +			,'ACCUMULATED TOTAL GRID SCALE PRECIPITATION')
	call ncnames_add_var('rain','Total Precipitation')
	call ncnames_add_var('rain','Precipitation')

	end subroutine ncnames_add_variables 

c*****************************************************************
c*****************************************************************
c*****************************************************************

