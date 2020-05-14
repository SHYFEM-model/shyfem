
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2013,2015-2020  Georg Umgiesser
!    Copyright (C) 2011,2013  Debora Bellafiore
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

c netcdf utility routines
c
c revision log :
c
c 05.12.2011	ggu&dbf	written from scratch
c 26.03.2012	ggu	standardized implicit routines, compiler warnings
c 20.09.2012	ggu	new routines for regular output
c 25.10.2012	ggu	changed VERS_6_1_59
c 21.01.2013	ggu	routines for handling scalar variables
c 25.01.2013	ggu	new part for nos variable initialization
c 28.01.2013	dbf	different types of vertical coordinates
c 03.05.2013	ggu	changed VERS_6_1_63
c 25.09.2013	ggu	new routines for writing time series
c 25.10.2013	ggu	changed VERS_6_1_68
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 15.04.2016	ggu	changed VERS_7_5_8
c 31.05.2016	ggu	changed time variable to double precision
c 07.06.2016	ggu	changed VERS_7_5_12
c 09.05.2017	ggu	changed VERS_7_5_26
c 16.05.2017	ggu	changed VERS_7_5_27
c 13.06.2017	ggu	changed VERS_7_5_29
c 11.07.2017	ggu	changed VERS_7_5_30
c 02.09.2017	ggu	changed VERS_7_5_31
c 17.11.2017	ggu	changed VERS_7_5_37
c 05.12.2017	ggu	changed VERS_7_5_39
c 24.01.2018	ggu	changed VERS_7_5_41
c 22.02.2018	ggu	changed VERS_7_5_42
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 14.05.2019	ggu	wrong definition of dimensions in nc_write_data_3d_reg
c 16.05.2019	ggu	new version of nc_rewrite_3d_reg(), new nc_set_quiet()
c 08.01.2020	ggu	allow for double user data
c 29.01.2020	ggu	insert extra information for error message
c
c notes :
c
c information on unstructured grids:
c https://publicwiki.deltares.nl/display/NETCDF/Unstructured+grids
c
c for non dimensional vertical coordinates (sigma etc) see:
c http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
c in appendic D.6 - D.9
c
c information on f77 interface description:
c http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77/
c
c CF compliance checker:
c http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl
c http://titania.badc.rl.ac.uk/cgi-bin/cf-checker.pl
c
c this file implements CF 1.4 compliance
c
c still to be implemented:
c	sigma/hybrid coordinates
c	get file name
c
c******************************************************************
c******************************************************************
c module
c******************************************************************
c******************************************************************

!==================================================================
	module netcdf_params
!==================================================================

	!use netcdf

	implicit none

        integer, save :: dimids_2d(5)		!dimensions for 2D case
        integer, save :: dimids_3d(5)		!dimensions for 3D case

        integer, save :: rec_varid		!id for time
        integer, save :: coord_varid(9)		!ids for coordinates

	character*80, save :: time_dim = ' '	!dimension name of time
	character*80, save :: time_var = ' '	!variable name of time

	logical, save :: bdebug_nc = .false.
	logical, save :: bquiet_nc = .false.

	INTERFACE
	subroutine nc_handle_err(errcode,string)
	integer errcode
	character*(*), optional :: string
	end subroutine
	END INTERFACE
	INTERFACE
	function nc_has_err(errcode)
	logical nc_has_err
	integer errcode
	end function
	END INTERFACE

	include 'netcdf.inc'

!==================================================================
	end module netcdf_params
!==================================================================

c******************************************************************
c******************************************************************
c open routines
c******************************************************************
c******************************************************************

	subroutine nc_init

	use netcdf_params

	time_dim = ' '
	time_var = ' '

	end

c******************************************************************

	subroutine nc_open_reg(ncid,nx,ny,nlv,flag,date0,time0,iztype)

! opens nc file for writing regular grid

	use netcdf_params

	implicit none

        integer ncid            !identifier (return)
	integer nx,ny,nlv	!size of arrays
	real flag		!flag for no data
        integer date0,time0     !date and time of time 0
	integer iztype          !type of vertical coordinates

	integer lmax
	integer lat_varid,lon_varid,lvl_varid,dep_varid
	integer varid
	integer lvl_dimid,nx_dimid,ny_dimid,rec_dimid
	integer ltext
	integer retval
	integer matrix_dimid(2)

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

c-----------------------------------------
C initialize parameters
c-----------------------------------------

	file_name = 'netcdf_reg.nc'
	file_name = 'netcdf.nc'
	file_name = 'out.nc'

c-----------------------------------------
C Create the file.
c-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_reg create')

c-----------------------------------------
C Define the dimensions. The record dimension is defined to have
C unlimited length - it can grow as needed. In this example it is
C the time dimension.
c-----------------------------------------

	lmax = max(1,nlv)	!be sure to have at least one layer

	retval = nf_def_dim(ncid, 'lon', nx, nx_dimid)
	call nc_handle_err(retval,'open_reg x')
	retval = nf_def_dim(ncid, 'lat', ny, ny_dimid)
	call nc_handle_err(retval,'open_reg y')
	retval = nf_def_dim(ncid, 'level', lmax, lvl_dimid)
	call nc_handle_err(retval,'open_reg level')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_reg time')

	matrix_dimid(1) = nx_dimid
	matrix_dimid(2) = ny_dimid

c-----------------------------------------
C Define the coordinate variables
c-----------------------------------------

c-----------------------------------------
c Assign units attributes to coordinate variables.
c-----------------------------------------

	retval = nf_def_var(ncid, 'lon', NF_REAL, 1, nx_dimid
     +				,lon_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'X'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'lat', NF_REAL, 1, ny_dimid
     +				,lat_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Y'
	call nc_define_attr(ncid,what,text,varid)

c---------------------
c for non dimensional vertical coordinates (sigma etc) see:
c http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
c in appendic D.6 - D.9
c---------------------

	retval = nf_def_var(ncid, 'level', NF_REAL, 1, lvl_dimid
     +				,lvl_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
	call make_vertical_coordinate(iztype,what,text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'bottom of vertical layers'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Z'
	call nc_define_attr(ncid,what,text,varid)

	what = 'positive'
	text = 'down'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'total_depth', NF_REAL, 2, matrix_dimid
     +				,dep_varid)
	call nc_handle_err(retval,'open_reg')
	varid = dep_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'sea_floor_depth_below_sea_surface'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'total depth at data grid nodes'
	call nc_define_attr(ncid,what,text,varid)

	call nc_define_range(ncid,-100.0,+10000.0,flag,varid)

c---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

	retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid
     +				,rec_varid)
	call nc_handle_err(retval,'open_reg')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

c-----------------------------------------
c define dimensions to remember
c-----------------------------------------

	dimids_2d(1) = nx_dimid
	dimids_2d(2) = ny_dimid
	dimids_2d(3) = rec_dimid

	dimids_3d(1) = nx_dimid
	dimids_3d(2) = ny_dimid
	dimids_3d(3) = lvl_dimid
	dimids_3d(4) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid
	coord_varid(3) = lvl_varid
	coord_varid(4) = dep_varid

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c******************************************************************

	subroutine nc_open_fem(ncid,nkn,nel,nlv,date0,time0,iztype)

! opens nc file for writing fem grid

	use netcdf_params

	implicit none

	integer ncid		!identifier (return)
	integer nkn,nel,nlv	!size of arrays
	integer date0,time0	!date and time of time 0
	integer iztype		!type of vertical coordinates

	integer lmax
	integer lat_varid,lon_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer varid
	integer lvl_dimid,node_dimid,elem_dimid,vertex_dimid,rec_dimid
	integer ltext
	integer retval
	integer eix_dimid(2)

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

c-----------------------------------------
C initialize parameters
c-----------------------------------------

	file_name = 'netcdf.nc'
	file_name = 'out.nc'

c-----------------------------------------
C Create the file.
c-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_fem')

c-----------------------------------------
C Define the dimensions. The record dimension is defined to have
C unlimited length - it can grow as needed. In this example it is
C the time dimension.
c-----------------------------------------

	lmax = max(1,nlv)	!be sure to have at least one layer

	retval = nf_def_dim(ncid, 'level', lmax, lvl_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'node', nkn, node_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'element', nel, elem_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'vertex', 3, vertex_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_fem')

	eix_dimid(1) = vertex_dimid
	eix_dimid(2) = elem_dimid

c-----------------------------------------
C Define the coordinate variables
c-----------------------------------------

c-----------------------------------------
c Assign units attributes to coordinate variables.
c-----------------------------------------

	retval = nf_def_var(ncid, 'longitude', NF_REAL, 1, node_dimid
     +				,lon_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'latitude', NF_REAL, 1, node_dimid
     +				,lat_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

c---------------------
c for non dimensional vertical coordinates (sigma etc) see:
c http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
c in appendic D.6 - D.9
c---------------------

	retval = nf_def_var(ncid, 'level', NF_REAL, 1, lvl_dimid
     +				,lvl_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
	call make_vertical_coordinate(iztype,what,text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'bottom of vertical layers'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Z'
	call nc_define_attr(ncid,what,text,varid)

	what = 'positive'
	text = 'down'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'total_depth', NF_REAL, 1, node_dimid
     +				,dep_varid)
	call nc_handle_err(retval,'open_fem')
	varid = dep_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'sea_floor_depth_below_sea_surface'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'total depth at nodes'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'element_index', NF_INT, 2, eix_dimid
     +				,eix_varid)
	call nc_handle_err(retval,'open_fem')
	varid = eix_varid

	what = 'units'
	text = '1'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'element index of nodes'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'maps every element to its three vertices'
	call nc_define_attr(ncid,what,text,varid)

	retval = nf_def_var(ncid, 'topology', NF_INT, 0, 0
     +				,top_varid)
	call nc_handle_err(retval,'open_fem')
	varid = top_varid

	what = 'units'
	text = '1'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'topology data of 2D unstructured mesh'
	call nc_define_attr(ncid,what,text,varid)

	what = 'dimensionality'
	text = '2'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

	retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid
     +				,rec_varid)
	call nc_handle_err(retval,'open_fem')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

c-----------------------------------------
c define dimensions to remember
c-----------------------------------------

	dimids_2d(1) = node_dimid
	dimids_2d(2) = rec_dimid

	dimids_3d(1) = lvl_dimid
	dimids_3d(2) = node_dimid
	dimids_3d(3) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid
	coord_varid(3) = lvl_varid
	coord_varid(4) = dep_varid
	coord_varid(5) = eix_varid
	coord_varid(6) = top_varid

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c******************************************************************

	subroutine nc_open_ts(ncid,node,date0,time0)

c opens nc file for time series write

	use netcdf_params

	implicit none

	integer ncid
	integer node		!total number of nodes for ts
	integer date0,time0

	integer lat_varid,lon_varid
	integer varid
	integer node_dimid,rec_dimid
	integer ltext
	integer retval

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

c-----------------------------------------
C initialize parameters
c-----------------------------------------

	file_name = 'netcdfts.nc'
	file_name = 'netcdf.nc'

c-----------------------------------------
C Create the file.
c-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_ts')

c-----------------------------------------
C Define the dimensions. The record dimension is defined to have
C unlimited length - it can grow as needed. In this example it is
C the time dimension.
c-----------------------------------------

	retval = nf_def_dim(ncid, 'node', 1, node_dimid)
	call nc_handle_err(retval,'open_ts')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_ts')

c-----------------------------------------
C Define the coordinate variables
c-----------------------------------------

c-----------------------------------------
c Assign units attributes to coordinate variables.
c-----------------------------------------

	retval = nf_def_var(ncid, 'longitude', NF_REAL, 1,node_dimid
     +				,lon_varid)
	write(6,*)'lon_varid',lon_varid
	call nc_handle_err(retval,'open_ts')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

	retval = nf_def_var(ncid, 'latitude', NF_REAL, 1, node_dimid
     +				,lat_varid)
	write(6,*)'lat_varid',lat_varid
	call nc_handle_err(retval,'open_ts')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

c---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

	retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid
     +				,rec_varid)
	call nc_handle_err(retval,'open_ts')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

c-----------------------------------------
c define dimensions to remember
c-----------------------------------------

	dimids_2d(1) = node_dimid
	dimids_2d(2) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************


c*****************************************************************

        subroutine make_vertical_coordinate(iztype,what,text)

c defines definition for vertical coordinate

        implicit none

        integer iztype
        character*(*) what
        character*(*) text

        what = 'standard_name'

        if( iztype .eq. 1 ) then        ! z-coordinates
          text = 'depth'
        else if( iztype .eq. 2 ) then   ! sigma-coordinates
          text = 'ocean_sigma_coordinate'
        else if( iztype .eq. 3 ) then   ! hybrid-coordinates
          text = 'ocean_sigma_z_coordinate'
        else
          write(6,*) 'iztype = ',iztype
          stop 'error stop make_vertical_coordinate: unknown iztype'
        end if

        end

c*****************************************************************
c*****************************************************************
c read dimensions
c*****************************************************************
c*****************************************************************

	subroutine nc_open_read(ncid,file)

! opens nc file for reading fem file

	use netcdf_params

	implicit none

	integer ncid		!identifier (return)
	character*(*) file
	integer retval

	!write(6,*) 'opening nc file for read: ',trim(file)

        retval = nf_open(file, NF_NOWRITE, ncid)
	call nc_handle_err(retval,'open_read')

	!retval = nf_inq(ncid, ndims, nvars, ngatts, unlim)
	!call nc_handle_err(retval,'open_read')

	end

c*****************************************************************

	subroutine nc_dims_info(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer ndims,dim_id,length,i
	character*30 name
	integer retval

	retval = nf_inq_ndims(ncid,ndims)
	call nc_handle_err(retval,'dims_info')

	write(6,*) 'dimensions: '
	do i=1,ndims
	  dim_id = i
	  retval = nf_inq_dim(ncid,dim_id,name,length)
	  call nc_handle_err(retval,'dims_info')
	  write(6,*) dim_id,length,name
	end do

	end

c*****************************************************************

        subroutine nc_get_dim_totnum(ncid,ndims)

	use netcdf_params

        implicit none

        integer ncid
	integer ndims

        integer retval

	retval = nf_inq_ndims(ncid,ndims)
	call nc_handle_err(retval,'get_dim_totnum')

	end

c*****************************************************************

        subroutine nc_get_dim_name(ncid,dim_id,name)

	use netcdf_params

        implicit none

        integer ncid
        integer dim_id
        character*(*) name

        integer retval

        retval = nf_inq_dimname(ncid,dim_id,name)
        call nc_handle_err(retval,'get_dim_name')

        end

c*****************************************************************

	subroutine nc_get_dim_id(ncid,name,dim_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	integer dim_id

	integer retval

	retval = nf_inq_dimid(ncid,name,dim_id)
	call nc_handle_err(retval,'get_dim_id')

	end

c*****************************************************************

        subroutine nc_has_dim_name(ncid,name,dim_id)

	use netcdf_params

        implicit none

        integer ncid
        character*(*) name
	integer dim_id

        integer retval

	retval = nf_inq_dimid(ncid,name,dim_id)
	if( nc_has_err(retval) ) dim_id = 0

        end

c*****************************************************************

	subroutine nc_get_dim_len(ncid,dim_id,dim_len)

	use netcdf_params

	implicit none

	integer ncid
	integer dim_id
	integer dim_len

	integer retval

	retval = nf_inq_dimlen(ncid,dim_id,dim_len)
	call nc_handle_err(retval,'has_dim_len')

	end

c*****************************************************************
c*****************************************************************
c read time records
c*****************************************************************
c*****************************************************************

	subroutine nc_get_time_rec(ncid,irec,t)

	use netcdf_params

	implicit none

	integer ncid
	integer irec			!number of record
	double precision t		!time found (return)

	integer nvars,var_id,time_id,i
	character*80 name
	character*80 time,time_d,time_v
	integer istart,icount
	integer xtype
	integer retval
	integer itime
	real rtime
	double precision dtime

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'get_time_rec')

	call nc_get_time_name(time_d,time_v)
	time = time_v

	call nc_get_var_id(ncid,time,time_id)
	if( time_id .eq. 0 ) then
	  stop 'error stop nc_get_time_rec: cannot find time variable'
	end if

	retval = nf_inq_vartype(ncid,time_id,xtype)
	call nc_handle_err(retval,'get_time_rec')
	!write(6,*) 'time_id: ',time_id,xtype

	istart = irec
	icount = 1
	if( xtype .eq. NF_INT ) then
	  !write(6,*) 'time is int.........'
	  retval = nf_get_vara_int(ncid,time_id,istart,icount,itime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = itime
	else if( xtype .eq. NF_FLOAT ) then
	  !write(6,*) 'time is real.........'
	  retval = nf_get_vara_real(ncid,time_id,istart,icount,rtime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = rtime
	else if( xtype .eq. NF_DOUBLE ) then
	  !write(6,*) 'time is double.........'
	  retval = nf_get_vara_double(ncid,time_id,istart,icount,dtime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = dtime
	else
	  write(6,*) 'xtype = ',xtype
	  stop 'error stop nc_get_time_rec: cannot read time'
	end if

	end

c*****************************************************************

	subroutine nc_get_time_recs(ncid,trecs)

	use netcdf_params

	implicit none

	integer ncid
	integer trecs

	integer dim_id,len
	character*80 name
	character*80 time,time_d,time_v
	integer retval

	trecs = 0

	call nc_get_time_name(time_d,time_v)
	if( time_d == ' ' .or. time_v == ' ' ) return

	call nc_get_dim_id(ncid,time_d,dim_id)

	if( dim_id .gt. 0 ) then
	  retval = nf_inq_dim(ncid,dim_id,time_d,len)
	  call nc_handle_err(retval,'get_time_recs')
	  trecs = len
	end if

	end

c*****************************************************************
c*****************************************************************
c read variables
c*****************************************************************
c*****************************************************************

	subroutine nc_var_info(ncid,var_id,bverb)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	logical bverb

	integer nvars,i,j,ia
	integer type,ndims,natts,xtype,length
	integer dimids(10)
	double precision avalue
	character*80 name,aname,atext
	integer retval

	logical, save :: blong = .false.

	blong = bverb

	  retval = nf_inq_var(ncid,var_id,name,type,ndims,dimids,natts)
	  call nc_handle_err(retval,'get_var_info')
	  if( ndims .gt. 10 ) stop 'error stop nc_var_info: ndims'
	  write(6,1010) var_id,natts,ndims,'   ',trim(name)
 1010     format(3i5,a,a)

	  if( blong ) then
	    do ia=1,natts
	      retval = nf_inq_attname(ncid,var_id,ia,aname)
	      call nc_handle_err(retval,'get_var_info')
	      retval = nf_inq_att(ncid,var_id,aname,xtype,length)
	      if( nc_has_err(retval) ) cycle	!no such attribute name
	      atext = ' '
	      call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	      if( xtype .ne. NF_CHAR ) then	!attribute is not a string
		write(atext,*) avalue
	      end if
	      write(6,1000) '    ',ia,'  ',trim(aname),'  ',trim(atext)
 1000	      format(a,i5,a,a,a,a)
	    end do
	  end if

	end

c*****************************************************************

	subroutine nc_vars_info(ncid,bverb)

	use netcdf_params

	implicit none

	integer ncid
	logical bverb

	integer nvars,var_id,i
	integer retval

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'vars_info')

	write(6,*) 'variables: '
	do i=1,nvars
	  var_id = i
	  call nc_var_info(ncid,var_id,bverb)
	end do

	end

c*****************************************************************

	subroutine nc_get_var_id(ncid,name,var_id)

c returns var_id = 0 if not found (no error)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	integer var_id

	integer retval

	retval = nf_inq_varid(ncid,name,var_id)
	if( retval .ne. 0 ) then
	  var_id = 0
	  return
	end if
	call nc_handle_err(retval,'get_var_id')

	end

c*****************************************************************

	subroutine nc_get_var_name(ncid,var_id,name)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) name

	integer retval

	retval = nf_inq_varname(ncid,var_id,name)
	call nc_handle_err(retval,'get_var_name')

	end

c*****************************************************************

	subroutine nc_get_var_totnum(ncid,nvars)

	use netcdf_params

	implicit none

	integer ncid
	integer nvars

	integer retval

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'get_var_totnum')

	end

c*****************************************************************

	subroutine nc_get_var_ndims(ncid,var_id,ndims,dimids)

c if ndims == 0 -> just compute total number of dimensions
c
c if ndims > 0 and smaller than total number of dimensions -> negative return value

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer ndims	!on entry dimension, on return number of dimensions
	integer dimids(ndims)

	integer retval
	integer ndim

	ndim = ndims

	retval = nf_inq_varndims(ncid,var_id,ndims)
	call nc_handle_err(retval,'get_var_ndims')

	if( ndim == 0 ) return
	if( ndims > ndim ) then
	  ndims = -ndims
	  return
	end if

	retval = nf_inq_vardimid(ncid,var_id,dimids)
	call nc_handle_err(retval,'get_var_ndims')

	end

c*****************************************************************

	subroutine nc_check_var_type(ncid,var_id,type)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) type

	integer xtype
	integer retval

	retval = nf_inq_vartype(ncid,var_id,xtype)
	call nc_handle_err(retval,'get_var_type')

	if( type .eq. 'integer' ) then
	  if( xtype .ne. NF_INT ) goto 99
	else if( type .eq. 'real' ) then
	  if( xtype .ne. NF_REAL ) goto 99
	else if( type .eq. 'double' ) then
	  if( xtype .ne. NF_DOUBLE ) goto 99
	else
	  write(6,*) 'type: ',type
	  stop 'error stop nc_check_type: cannot check type'
	end if

	return
   99	continue
	write(6,*) 'type: ',type
	write(6,*) 'xtype: ',xtype
	stop 'error stop nc_check_type: type mismatch'
	end

c*****************************************************************

	subroutine nc_get_global_attr(ncid,aname,atext)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) aname
	character*(*) atext

	call nc_get_var_attr(ncid,NF_GLOBAL,aname,atext)

	end

c*****************************************************************

	subroutine nc_get_var_attr(ncid,var_id,aname,atext)

! in this routine is a gfortran bug !!!FIXME

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) aname
	character*(*) atext

	logical debug
	integer retval
	integer xtype,ll,cl,sl
	integer, parameter :: cmax = 1000
	character*(cmax) aux
	!character(len=:), allocatable :: aux

	!debug = .true.		!GGU
	!debug = .false.	!GGU
	debug = bdebug_nc

	atext = ' '
	retval = nf_inq_att(ncid,var_id,aname,xtype,ll)
	if( nc_has_err(retval) ) return	!no such attribute name
	if( xtype .ne. NF_CHAR ) return		!attribute is not a string

	if( ll > cmax ) then
	  write(6,*) 'll = ',ll,'  cmax = ',cmax
	  write(6,*) 'variable is too short to receive string'
	  stop 'error stop nc_get_var_attr: ll > cmax'
	end if

	cl = ll + 1
	!allocate(character(len=cl) :: aux)
	aux = ' '
	sl = len(aux)

	if( debug ) write(6,*) 'cl,sl: ',cl,sl,len(aux)

	if( sl < cl ) then
	  write(6,*) 'sl,cl: ',sl,cl
	  write(6,*) 'compiler cannot allocate variable string length'
	  write(6,*) 'please use constant string variable'
	  stop 'error stop nc_get_var_attr: len(aux)<cl'
	end if

	if( debug ) then
	write(6,*) 'cl: ',sl,len(aux),cl
	write(6,*) 'nnn: ',ll,var_id,trim(aname),len(aux)
	end if

	retval = nf_get_att_text(ncid,var_id,aname,aux)

	sl = len(aux)
	if( ichar(aux(ll:ll)) == 0 ) aux(ll:ll) = ' '
	if( ichar(aux(sl:sl)) == 0 ) aux(sl:sl) = ' '
	if( debug ) then
	  write(6,*) 'auxaux1: ',trim(aux)
	  write(6,*) 'auxaux2: ',len_trim(aux),ll,sl
	  write(6,*) 'auxaux3: ',ichar(aux(ll:ll))
	  write(6,*) 'nnn: end of call'
	end if
	call nc_handle_err(retval,'get_var_attr')
	atext = aux

	!deallocate(aux)

	end

c*****************************************************************

	function nc_has_var_attrib(ncid,var_id,aname)

	use netcdf_params

	implicit none

	logical nc_has_var_attrib
	integer ncid
	integer var_id
	character*(*) aname

	integer retval
	integer xtype,len

	nc_has_var_attrib = .false.

	retval = nf_inq_att(ncid,var_id,aname,xtype,len)
	if( nc_has_err(retval) ) return	!no such attribute name

	nc_has_var_attrib = .true.

	end function nc_has_var_attrib

c*****************************************************************

	subroutine nc_get_var_attrib(ncid,var_id,aname,atext,avalue)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) aname
	character*(*) atext
	double precision avalue

	integer retval
	integer xtype,len
	double precision, allocatable :: daux(:)
	character(len=:), allocatable :: saux		!NEMUNAS_FIX_NEW
	!character*1000 :: saux				!NEMUNAS_FIX_OLD

	atext = ' '
	avalue = 0.

	retval = nf_inq_att(ncid,var_id,aname,xtype,len)
	if( nc_has_err(retval) ) return	!no such attribute name

	if( xtype .eq. NF_CHAR ) then
	  allocate(character(len=len) :: saux)		!NEMUNAS_FIX_NEW
	  !if( len > 1000 ) goto 99			!NEMUNAS_FIX_OLD
	  retval = nf_get_att_text(ncid,var_id,aname,saux)
	  if( len == 1 .and. ichar(saux(1:1)) == 0 ) saux=' '	!FIX
	  atext = saux
	else
	  if( len > 1 ) then
	    write(6,*) '*** attributes with multiple values'
	    write(6,*) '*** only returning first of array'
	    write(6,*) '   name,len: ',trim(aname),len
	    allocate(daux(len))
	    retval = nf_get_att_double(ncid,var_id,aname,daux)
	    avalue = daux(1)
	  else
	    retval = nf_get_att_double(ncid,var_id,aname,avalue)
	  end if
	end if

	call nc_handle_err(retval,'get_var_attrib')

	return
   99	continue
	write(6,*) 'len = ',len,'  max possible is 1000'
	stop 'error stop nc_get_var_attrib: len > 1000'
	end

c*****************************************************************

	subroutine nc_get_var_int(ncid,var_id,data)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer data(*)

	integer retval

	retval = nf_get_var_int(ncid,var_id,data)
	call nc_handle_err(retval,'get_var_int')

	end

c*****************************************************************

	subroutine nc_get_var_real(ncid,var_id,data)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	real data(*)

	integer retval

	retval = nf_get_var_real(ncid,var_id,data)
	call nc_handle_err(retval,'get_var_real')

	end

c*****************************************************************

	subroutine nc_get_var_data(ncid,name,trec,ndim,ndimens
     +					,dims,data)

c reads time record trec of variable name

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name	!name of variable to read
	integer trec		!number of time record to read
	integer ndim		!dimension of data array
	integer ndimens		!expected dimensionality of data to read
				!(this should exclude time dimension)
	integer dims(ndimens+1)	!length of dimensions (return)
	real data(ndim)		!data (return)

	integer retval
	integer i,dim_id,dim_len
	integer var_id,itime
	integer ndims,nlength
	integer xtype,size
	integer, allocatable :: icount(:)
	integer, allocatable :: istart(:)
	integer, allocatable :: dimids(:)
	character*80, allocatable :: dimn(:)
	character*80 dimname
	character*80 auxname
	character*80 time,time_d,time_v

	double precision, allocatable :: ddata(:)

	call nc_get_time_name(time_d,time_v)
	time = time_d
	!write(6,*) 'name of time dim: ',time

	ndims = 0
	call nc_get_var_id(ncid,name,var_id)
	if( var_id .le. 0 ) then
	  write(6,*) 'Cannot find variable: ',name
	  return
	end if

	retval = nf_inq_varndims(ncid,var_id,ndims)
	call nc_handle_err(retval,'get_var_data')
	if( ndims > ndimens+1 ) goto 95
	allocate(icount(ndims),istart(ndims),dimids(ndims),dimn(ndims))

	retval = nf_inq_vardimid(ncid,var_id,dimids)
	call nc_handle_err(retval,'get_var_data')

	itime = 0
	nlength = 1
	do i=1,ndims
	  dim_id = dimids(i)
	  call nc_get_dim_name(ncid,dim_id,dimname)
	  call nc_get_dim_len(ncid,dim_id,dim_len)
	  dims(i) = dim_len
	  dimn(i) = dimname
	  !write(6,*) 'name of dim: ',i,dimname(1:30)
	  if( dimname .eq. time ) then
	    itime = i
	  else
	    istart(i) = 1
	    icount(i) = dim_len
	    nlength = nlength * dim_len
	  end if
	  !write(6,*) i,dim_id,dim_len,dimname(1:20)
	end do

	if( itime .eq. 0 ) then
	  !write(6,*) 'variable without time... ',itime
	else
	  !write(6,*) 'variable has time... ',itime
	  if( itime .ne. ndims ) goto 99
	  ndims = ndims - 1		! time is always last
	  istart(itime) = trec
	  icount(itime) = 1
	end if

	if( ndims .ne. ndimens ) then
	  write(6,*) 'expected dimension of data: ',ndimens
	  write(6,*) 'real dimension of data: ',ndims
	  do i=1,ndims
	    write(6,*) i,dims(i),'  ',trim(dimn(i))
	  end do
	  stop 'error stop nc_get_var_data: ndimens'
	end if

	if( nlength .gt. ndim ) then
	  write(6,*) 'size of data to be read: ',nlength
	  write(6,*) 'dimension of data array: ',ndim
	  stop 'error stop nc_get_var_data: ndim'
	end if

	retval = nf_inq_vartype(ncid,var_id,xtype)
	call nc_handle_err(retval,'get_var_data')

	if( xtype .eq. NF_FLOAT ) then
	  retval = nf_get_vara_real(ncid,var_id,istart,icount,data)
	  call nc_handle_err(retval,'get_var_data')
	else if( xtype .eq. NF_DOUBLE ) then
	  allocate(ddata(nlength))
	  retval = nf_get_vara_double(ncid,var_id,istart,icount,ddata)
	  call nc_handle_err(retval,'get_var_data')
	  data = ddata
	  deallocate(ddata)
	else
	  write(6,*) 'cannot handle data type...',xtype
	  stop 'error stop nc_get_var_data: data type'
	end if

	return
   95	continue
	write(6,*) 'dimensions to small: ',ndims,ndimens+1
	stop 'error stop nc_get_var_data: ndimens'
   99	continue
	write(6,*) 'time variable is not last: ',itime,ndims
	stop 'error stop nc_get_var_data: itime'
	end

c*****************************************************************
c*****************************************************************
c define variables for write
c*****************************************************************
c*****************************************************************

	subroutine nc_define_2d_reg(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_2d
     +				,var_id)
	call nc_handle_err(retval,'define_2d_reg')

	end

c*****************************************************************

	subroutine nc_define_3d_reg(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 4, dimids_3d
     +				,var_id)
	call nc_handle_err(retval,'define_3d_reg')

	end

c*****************************************************************

	subroutine nc_define_2d(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 2, dimids_2d
     +				,var_id)
	call nc_handle_err(retval,'define_2d')

	end

c*****************************************************************

	subroutine nc_define_3d(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_3d
     +				,var_id)
	call nc_handle_err(retval,'define_3d')

	end

c*****************************************************************
c*****************************************************************
c define attributes for write
c*****************************************************************
c*****************************************************************

	subroutine nc_define_attr(ncid,what,def,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what,def
	integer var_id

	integer ldef
	integer retval

	integer nc_ichanm

        ldef = nc_ichanm(def)
	retval = nf_put_att_text(ncid, var_id, what, ldef
     +				,def)
	call nc_handle_err(retval,'define_attr')

	end

c*****************************************************************

	subroutine nc_define_attr_real(ncid,what,value,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	real value
	integer var_id

	integer ldef
	integer retval

        ldef = 1
	retval = nf_put_att_real(ncid, var_id, what, NF_FLOAT, ldef
     +				,value)
	call nc_handle_err(retval,'define_attr_real')

	end

c*****************************************************************

	subroutine nc_define_range(ncid,rmin,rmax,flag,var_id)

	use netcdf_params

	implicit none

	integer ncid
	real rmin,rmax,flag
	integer var_id

	real rminmax(2)
	integer retval

	rminmax(1) = rmin
	rminmax(2) = rmax

c	retval = nf_put_att_real(ncid, var_id, 'valid_range', NF_REAL
c     +				,2,rminmax)
c	call nc_handle_err(retval,'define_range')

	retval = nf_put_att_real(ncid, var_id, 'valid_min', NF_REAL
     +				,1,rmin)
	call nc_handle_err(retval,'define_range')

	retval = nf_put_att_real(ncid, var_id, 'valid_max', NF_REAL
     +				,1,rmax)
	call nc_handle_err(retval,'define_range')

c	retval = nf_put_att_real(ncid, var_id, 'missing_value', NF_REAL
c     +				,1,flag)
c	call nc_handle_err(retval,'define_range')

	retval = nf_put_att_real(ncid, var_id, '_FillValue', NF_REAL
     +				,1,flag)
	call nc_handle_err(retval,'define_range')

	end

c*****************************************************************

	subroutine nc_end_define(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer retval

	retval = nf_enddef(ncid)
	call nc_handle_err(retval,'end_define')

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_has_vertical_dimension(ncid,name,bvert)

c checks if variable is 3d

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	logical bvert

	logical btime
	integer var_id,ndims
	integer retval
	integer, allocatable :: dim_ids(:)

	call nc_get_var_id(ncid,name,var_id)

        retval = nf_inq_varndims(ncid,var_id,ndims)
        call nc_handle_err(retval,'has_vertical_dimension')
	allocate(dim_ids(ndims))

	call nc_get_var_ndims(ncid,var_id,ndims,dim_ids)

	call nc_has_time_dimension(ncid,name,btime)
	if( btime ) ndims = ndims - 1

	bvert = ( ndims == 3 )

	end

c*****************************************************************

	subroutine nc_has_time_dimension(ncid,name,btime)

c checks if variable name has time dimension

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	logical btime

	integer var_id,ndims,time_id
	integer retval
	integer, allocatable :: dim_ids(:)
	character*80 tname,time_d,time_v

	call nc_get_var_id(ncid,name,var_id)
        retval = nf_inq_varndims(ncid,var_id,ndims)
        call nc_handle_err(retval,'has_time_dimension')
        allocate(dim_ids(ndims))

	call nc_get_var_ndims(ncid,var_id,ndims,dim_ids)

	time_id = dim_ids(ndims)			!time is always last
        call nc_get_dim_name(ncid,time_id,tname)
	call nc_get_time_name(time_d,time_v)

	btime = tname .eq. time_d

	end

c*****************************************************************

	subroutine nc_set_time_name(time_d,time_v)

	use netcdf_params

	implicit none

	character*(*) time_d,time_v

	if( time_d .ne. ' ' )  time_dim  = time_d	!dimension name of time
	if( time_v .ne. ' ' )  time_var  = time_v	!variable name of time

	end

c*****************************************************************

	subroutine nc_get_time_name(time_d,time_v)

	use netcdf_params

	implicit none

	character*(*) time_d,time_v

	time_d = time_dim
	time_v = time_var

	end

c*****************************************************************
c*****************************************************************
c write coordinate information
c*****************************************************************
c*****************************************************************

	subroutine nc_write_coords_reg(ncid,nx,ny,nlv
     +					,xlon,ylat,depth,hlv)

	use netcdf_params

	implicit none

	integer ncid
	integer nx,ny
	integer nlv
	real xlon(nx)
	real ylat(ny)
	real depth(nx,ny)
	real hlv(nlv)

	integer lon_varid,lat_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)
	lvl_varid = coord_varid(3)
	dep_varid = coord_varid(4)

c-----------------------------------------
C write coordinate data
c-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, xlon)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, lat_varid, ylat)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, dep_varid, depth)
	call nc_handle_err(retval,'write_coords_reg')

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************

	subroutine nc_write_coords_fem(ncid,nkn,nel,nlv
     +				,xgv,ygv,hkv,nen3v,hlv)

	!use mod_depth
	!use levels
	!use basin
	use netcdf_params

	implicit none

	integer ncid
	integer nkn,nel,nlv
	real xgv(nkn),ygv(nkn)
	real hkv(nkn)
	integer nen3v(3,nel)
	real hlv(nlv)

	integer lon_varid,lat_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)
	lvl_varid = coord_varid(3)
	dep_varid = coord_varid(4)
	eix_varid = coord_varid(5)
	top_varid = coord_varid(6)

c-----------------------------------------
C write coordinate data
c-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, xgv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, lat_varid, ygv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, dep_varid, hkv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_int(ncid, eix_varid, nen3v)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_int(ncid, top_varid, 2)
	call nc_handle_err(retval,'write_coords_fem')

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************

	subroutine nc_write_coords_ts(ncid,lon,lat)

	use netcdf_params

	implicit none

	integer ncid
        real lon(1), lat(1)

	integer lon_varid,lat_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)

c-----------------------------------------
C write coordinate data
c-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, lon)
	call nc_handle_err(retval,'write_coords_ts')

	retval = nf_put_var_real(ncid, lat_varid, lat)
	call nc_handle_err(retval,'write_coords_ts')

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************
c*****************************************************************
c write time record
c*****************************************************************
c*****************************************************************

	subroutine nc_write_time(ncid,irec,it)

	use netcdf_params

	implicit none

	integer ncid
	integer irec
	integer it

	integer retval
	double precision dtime

	dtime = it
	retval = nf_put_vara_double(ncid, rec_varid, irec, 1, dtime)
	call nc_handle_err(retval,'write_time')

	end

c*****************************************************************

	subroutine nc_write_dtime(ncid,irec,dtime)

	use netcdf_params

	implicit none

	integer ncid
	integer irec
	double precision dtime

	integer retval

	retval = nf_put_vara_double(ncid, rec_varid, irec, 1, dtime)
	call nc_handle_err(retval,'write_dtime')

	end

c*****************************************************************
c*****************************************************************
c write data
c*****************************************************************
c*****************************************************************

	subroutine nc_write_data_2d_reg(ncid,var_id,irec,nx,ny,var2d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nx,ny
	real var2d(nx,ny)

	integer retval
	integer count(3)
	integer start(3)

	count(1) = nx
	count(2) = ny
	count(3) = 1
	start(1) = 1
	start(2) = 1
	start(3) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var2d)
	call nc_handle_err(retval,'write_data_2d_reg')

	end

c*****************************************************************

	subroutine nc_write_data_3d_reg(ncid,var_id,irec,nlv,nx,ny,var3d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nlv
	integer nx,ny
	!real var3d(nlv,nx,ny)
	real var3d(nx,ny,nlv)

	integer retval
	integer count(4)
	integer start(4)

	count(1) = nx
	count(2) = ny
	count(3) = nlv
	count(4) = 1
	start(1) = 1
	start(2) = 1
	start(3) = 1
	start(4) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var3d)
	call nc_handle_err(retval,'write_data_3d_reg')

	end

c*****************************************************************

	subroutine nc_rewrite_3d_reg(nlvddi,lmax,nx,ny,var3d,vnc3d)

c re-writes a 3d array to be CF compliant

	use netcdf_params

	implicit none

	integer nlvddi,lmax
	integer nx,ny
	real var3d(nlvddi,nx,ny)
	real vnc3d(nx,ny,lmax)

	integer i,j,l

	do j=1,ny
	  do i=1,nx
	    do l=1,lmax
	      vnc3d(i,j,l) = var3d(l,i,j)
	    end do
	  end do
	end do

	end

c*****************************************************************

	subroutine nc_write_data_2d(ncid,var_id,irec,nkn,var2d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nkn
	real var2d(nkn)

	integer retval
	integer count(2)
	integer start(2)

	count(1) = nkn
	count(2) = 1
	start(1) = 1
	start(2) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var2d)
	call nc_handle_err(retval,'write_data_2d')

	end

c*****************************************************************

	subroutine nc_write_data_3d(ncid,var_id,irec,nlv,nkn,var3d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nlv
	integer nkn
	real var3d(nlv,nkn)

	integer retval
	integer count(3)
	integer start(3)

	count(1) = nlv
	count(2) = nkn
	count(3) = 1
	start(1) = 1
	start(2) = 1
	start(3) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var3d)
	call nc_handle_err(retval,'write_data_3d')

	end

c*****************************************************************
c*****************************************************************
c compact data for write
c*****************************************************************
c*****************************************************************

	subroutine nc_compact_3d_reg(nlvddi,nlv,nx,ny,var_in,var_out)

	implicit none

	integer nlvddi
	integer nlv,nx,ny
	real var_in(nlvddi,nx,ny)
	real var_out(nlv,nx,ny)

	integer l,ix,iy

	do iy=1,ny
	 do ix=1,nx
	  do l=1,nlv
	    var_out(l,ix,iy) = var_in(l,ix,iy)
	  end do
	 end do
	end do

	end

c*****************************************************************

	subroutine nc_compact_3d(nlvddi,nlv,nkn,var_in,var_out)

	implicit none

	integer nlvddi
	integer nlv,nkn
	real var_in(nlvddi,nkn)
	real var_out(nlv,nkn)

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    var_out(l,k) = var_in(l,k)
	  end do
	end do

	end

c*****************************************************************
c*****************************************************************
c various routines
c*****************************************************************
c*****************************************************************

	subroutine nc_close(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer retval

	retval = nf_close(ncid)
	call nc_handle_err(retval,'close')

	end

c*****************************************************************

	subroutine nc_global(ncid,title)

c writes global conventions

	use netcdf_params

	implicit none

	integer ncid
	character*(*) title

	integer ltext,retval,varid
	character*80 text
	character*80 what
	character*80 cdate

	integer nc_ichanm

	varid = NF_GLOBAL

	what = 'Conventions'
	text = 'CF-1.4'
	call nc_define_attr(ncid,what,text,varid)

	what = 'title'
	text = title
	call nc_strip(text)
	call nc_define_attr(ncid,what,text,varid)

	call nc_current_time(cdate)
	what = 'history'
	text = 'created on ' // trim(cdate)
	call nc_strip(text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'institution'
	text = 'ISMAR-CNR, Venice, Italy'
	call nc_define_attr(ncid,what,text,varid)

	what = 'source'
	text = 'Model data produced by SHYFEM at ISMAR-CNR'
	call nc_define_attr(ncid,what,text,varid)

	what = 'references'
	text = 'Model info: http://www.ismar.cnr.it/shyfem'
	call nc_define_attr(ncid,what,text,varid)

	what = 'contact'
	text = 'email: georg.umgiesser@ismar.cnr.it'
	call nc_define_attr(ncid,what,text,varid)

	what = 'comment'
	text = 'Data restriction: for academic research use only'
	call nc_define_attr(ncid,what,text,varid)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_handle_err(errcode,string)

	use netcdf_params, only : nf_noerr,nf_strerror

	implicit none

	integer errcode
	character*(*), optional :: string

	if( errcode .eq. nf_noerr ) return

	if( present(string) ) write(6,*) trim(string)
	write(6,*) 'Error: ', nf_strerror(errcode)

	stop 'error stop nc_handle_err'
	end

c*****************************************************************

	function nc_has_err(errcode)

	use netcdf_params, only : nf_noerr

	implicit none

	logical nc_has_err
	integer errcode

	nc_has_err = ( errcode .ne. nf_noerr )

	end

c*****************************************************************
c*****************************************************************
c string utility routines
c*****************************************************************
c*****************************************************************

        function nc_ichanm(line)

c computes length of line without trailing blanks
c
c line          line of text
c ichanm        length of line (return value)
c               ... 0 : line is all blank

	implicit none

	integer nc_ichanm
        character*(*) line

	integer i,ndim
        character*1 blank,tab,char
        data blank /' '/

	tab = char(9)

        ndim=len(line)

        do i=ndim,1,-1
          if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
        end do

    1   continue
        nc_ichanm=i

        return
        end

c*****************************************************************

        subroutine nc_strip(line)

c strip blank lines
c
c line          line of text
c ichanm        length of line (return value)
c               ... 0 : line is all blank

	implicit none

        character*(*) line

	integer i,ndim
        character*1 blank,tab,char
        data blank /' '/

	tab = char(9)

        ndim=len(line)

        do i=1,ndim
          if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
        end do

    1   continue
	if( i .gt. ndim ) return
	if( i .eq. 1 ) return

	line(1:) = line(i:)

        return
        end

c*****************************************************************

	subroutine nc_subst_char(line,orig,subst)

c substitutes in line character orig with subst

	implicit none

        character*(*) line
        character*1 orig,subst

	integer i,n

	n = len(line)

	do i=1,n
	  if( line(i:i) .eq. orig ) line(i:i) = subst
	end do

	end

c*****************************************************************
c*****************************************************************
c date and time routines
c*****************************************************************
c*****************************************************************

	subroutine nc_unpack_date(ipack,i1,i2,i3)

c unpacks date in ipack to integers

	implicit none

	integer ipack,i1,i2,i3

	integer iaux

	i1 = 0
	i2 = 0
	i3 = 0

	iaux = ipack
	i1 = iaux/10000
	iaux = iaux - 10000*i1
	i2 = iaux/100
	iaux = iaux - 100*i2
	i3 = iaux

	end

c*****************************************************************

	subroutine nc_format_date(date,year,month,day,hour,min,sec,zone)

c formats date string

	implicit none

	character*(*) date
	integer year,month,day,hour,min,sec
	character*(*) zone	!time zone - can be blank (UTC,MET,etc.)

	date = ' '
	write(date,'(i4,a1,i2,a1,i2)') year,'-',month,'-',day
	call nc_subst_char(date(1:10),' ','0')
	write(date(12:19),'(i2,a1,i2,a1,i2)') hour,':',min,':',sec
	call nc_subst_char(date(12:19),' ','0')
	date(21:) = zone

	end

c*****************************************************************

	subroutine nc_convert_date(date0,time0,date)

c converts data n integer to character

	implicit none

	integer date0,time0
	character*(*) date

	integer aux,year,month,day
	integer hour,min,sec
	integer i

c-----------------------------------------------
c convert date and time
c-----------------------------------------------

	aux = date0
	if( aux .lt. 10000 ) aux = 10000 * aux
	call nc_unpack_date(aux,year,month,day)
	if( month .le. 0 ) month = 1
	if( day .le. 0 ) day = 1
	call nc_unpack_date(time0,hour,min,sec)

	!write(6,*) 'date0: ',date0,year,month,day
	!write(6,*) 'time0: ',time0,hour,min,sec

c-----------------------------------------------
c prepare date
c-----------------------------------------------

	call nc_format_date(date,year,month,day,hour,min,sec,'UTC')

c-----------------------------------------------
c end of routine
c-----------------------------------------------

	end

c*****************************************************************

	subroutine nc_current_time(cdate)

	use netcdf_params

	implicit none

	character*(*) cdate

	integer year,month,day,hour,min,sec
	character*10 date,time,zone
	integer value(8)

	call date_and_time(date,time,zone,value)

	year = value(1)
	month = value(2)
	day = value(3)
	hour = value(5)
	min = value(6)
	sec = value(7)

	!if( year .lt. 1000 ) year = 2000 + year

	call nc_format_date(cdate,year,month,day,hour,min,sec,'MET')

	if( .not. bquiet_nc ) then
	  write(6,*) 'cdate: ',trim(cdate),'  ',date,'  ',time
	end if

	end
	
c*****************************************************************
c*****************************************************************
c variable initialization
c*****************************************************************
c*****************************************************************

	subroutine nc_init_variable(ncid,breg,dim,ivar,flag,var_id)

	use netcdf_params

	implicit none

	integer ncid
	logical breg
	integer dim
	integer ivar
	real flag
	integer var_id		! id to be used for other operations (return)

	character*80 name,what,std,units
	real cmin,cmax

	if( ivar .eq. 1 ) then		! water level
	  name = 'water_level'
	  what = 'standard_name'
          std = 'water_surface_height_above_reference_datum'
          units = 'm'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 2 ) then	! x-velocity
	  name = 'u_velocity'
	  what = 'standard_name'
	  std = 'eastward_sea_water_velocity'
	  units = 'm s-1'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 3 ) then	! y-velocity
	  name = 'v_velocity'
	  what = 'standard_name'
	  std = 'northward_sea_water_velocity'
	  units = 'm s-1'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 4 ) then	! x-velocity (no tide)
	  name = 'u_velocity'
	  what = 'standard_name'
	  std = 'eastward_sea_water_velocity_assuming_no_tide'
	  units = 'm s-1'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 5 ) then	! y-velocity (no tide)
	  name = 'v_velocity'
	  what = 'standard_name'
	  std = 'northward_sea_water_velocity_assuming_no_tide'
	  units = 'm s-1'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 367 ) then	! generic tracer
	  name = 'tracer'
	  what = 'long_name'
	  std = 'river_influence_index'
	  units = '1'
	  cmin = 0.
	  cmax = 20.
	else if( ivar .eq. 10 ) then	! generic tracer
	  name = 'tracer'
	  what = 'long_name'
	  std = 'generic_tracer'
	  units = '1'
	  cmin = 0.
	  cmax = 110.
	else if( ivar .eq. 11 ) then	! salinity
	  name = 'salinity'
	  what = 'standard_name'
	  std = 'sea_water_salinity'
	  units = '1e-3'
	  cmin = 0.
	  cmax = 200.
	else if( ivar .eq. 12 ) then	! temperature
	  name = 'temperature'
	  what = 'standard_name'
	  std = 'sea_water_temperature'
	  units = 'degC'
	  cmin = -10.
	  cmax = 100.
	else if( ivar .eq. 19 ) then	! vorticity
	  name = 'vorticity'
	  what = 'standard_name'
	  std = 'ocean_relative_vorticity'
	  units = 's-1'
	  cmin = -100.
	  cmax = 100.
	else if( ivar .eq. 99 ) then	! wrt
	  name = 'water_renewal_time'
	  what = 'long_name'
	  std = 'water_renewal_time'
	  units = '1'
	  cmin = 0.
	  cmax = 500.
	else if( ivar .eq. 81 ) then	! particle age
	  name = 'particle_age'
	  what = 'long_name'
	  std = 'lagrangian_particle_age'
	  units = 'day'
	  cmin = 0.
	  cmax = 50000.
	else if( ivar .eq. 84 ) then	! particle density
	  name = 'particle_density'
	  what = 'long_name'
	  std = 'lagrangian_particle_concentration'
	  units = 'num km-2'
	  cmin = 0.
	  cmax = 50000.
	else
	  write(6,*) 'unknown variable: ',ivar
	  stop 'error stop descr_var'
	end if

	if( .not. bquiet_nc ) then
	  write(6,*) 'nc variable ',ivar,trim(name)
	end if

	if( dim .eq. 2 ) then
	  if( breg ) then
	    call nc_define_2d_reg(ncid,name,var_id)
	  else
	    call nc_define_2d(ncid,name,var_id)
	  end if
	else if( dim .eq. 3 ) then
	  if( breg ) then
	    call nc_define_3d_reg(ncid,name,var_id)
	  else
	    call nc_define_3d(ncid,name,var_id)
	  end if
	end if

	call nc_define_attr(ncid,'units',units,var_id)
	call nc_define_attr(ncid,what,std,var_id)
	call nc_define_range(ncid,cmin,cmax,flag,var_id)

	end

c*****************************************************************
c*****************************************************************
c auxiliary routines
c*****************************************************************
c*****************************************************************

	subroutine nc_set_quiet(bquiet)

	use netcdf_params

	implicit none

	logical bquiet

	bquiet_nc = bquiet

	end

c*****************************************************************
c*****************************************************************
c test routines
c*****************************************************************
c*****************************************************************

	subroutine wrnetcdf
	implicit none
	end

c*****************************************************************

	function next_record(it,znv)

	implicit none

	logical next_record
	integer it
	real znv(1)

	next_record = .true.

	end

c*****************************************************************

	subroutine test_nc

! this is not working - we first have to read a basin...

!	use basin
!	use levels
!	use mod_depth

	implicit none

	integer ncid
	integer dimids_2d(2)
	integer coord_varid(5)
	integer rec_varid
	integer level_id

	integer it
	integer irec
	integer date0,time0
	integer iztype		!type of vertical coordinates
	real znv(1)
	logical next_record
	character*1 units

	date0 = 20120101
	time0 = 0
	units = 'm'
	iztype = 1

!	call nc_open_fem(ncid,nkn,nel,nlv,date0,time0,iztype)
!	call nc_define_2d(ncid,'water_level',level_id)
!	call nc_define_attr(ncid,'units',units,level_id)
!	call nc_end_define(ncid)
!	call nc_write_coords_fem(ncid,nkn,nel,nlv
!     +			,xgv,ygv,hkv,nen3v,hlv)

!	irec = 0
!	do while( next_record(it,znv) )
!
!	  irec = irec + 1
!	  call nc_write_time(ncid,irec,it)
!	  call nc_write_data_2d(ncid,level_id,irec,nkn,znv)
!
!	end do
!
!	call nc_close(ncid)

	end

c******************************************************************

c	program test_nc_main
c	implicit none
c	call test_nc
c	end

c******************************************************************

