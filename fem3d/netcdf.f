c
c netcdf utility routines
c
c revision log :
c
c 05.12.2011    ggu&dbf	written from scratch
c 26.03.2012    ggu	standardized implicit routines, compiler warnings
c 20.09.2012    ggu	new routines for regular output
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

	subroutine nc_open_reg(ncid,nx,ny,nlv,flag,date0,time0)

	implicit none

	include 'netcdf.inc'
	include 'param.h'
	include 'netcdf.h'

        integer ncid            !identifier (return)
	integer nx,ny,nlv	!size of arrays
	real flag
        integer date0,time0     !date and time of time 0

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

c-----------------------------------------
C Create the file.
c-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval)

c-----------------------------------------
C Define the dimensions. The record dimension is defined to have
C unlimited length - it can grow as needed. In this example it is
C the time dimension.
c-----------------------------------------

	retval = nf_def_dim(ncid, 'level', nlv, lvl_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'lon', nx, nx_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'lat', ny, ny_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval)

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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
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
	call nc_handle_err(retval)
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

	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
     +				,rec_varid)
	call nc_handle_err(retval)
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//date
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
c define dimensions to pass back
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

	subroutine nc_open(ncid,nkn,nel,nlv,date0,time0)

	implicit none

	include 'netcdf.inc'
	include 'param.h'
	include 'netcdf.h'

	integer ncid		!identifier (return)
	integer nkn,nel,nlv	!size of arrays
	integer date0,time0	!date and time of time 0

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

c-----------------------------------------
C Create the file.
c-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval)

c-----------------------------------------
C Define the dimensions. The record dimension is defined to have
C unlimited length - it can grow as needed. In this example it is
C the time dimension.
c-----------------------------------------

	retval = nf_def_dim(ncid, 'level', nlv, lvl_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'node', nkn, node_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'element', nel, elem_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'vertex', 3, vertex_dimid)
	call nc_handle_err(retval)
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval)

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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
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
	call nc_handle_err(retval)
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

	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
     +				,rec_varid)
	call nc_handle_err(retval)
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//date
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
c define dimensions to pass back
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

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_define_2d_reg(ncid,what,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_2d
     +				,var_id)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_define_3d_reg(ncid,what,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 4, dimids_3d
     +				,var_id)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_define_2d(ncid,what,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 2, dimids_2d
     +				,var_id)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_define_3d(ncid,what,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

	retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_3d
     +				,var_id)
	call nc_handle_err(retval)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_define_attr(ncid,what,def,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	character*(*) what,def
	integer var_id

	integer ldef
	integer retval

	integer nc_ichanm

        ldef = nc_ichanm(def)
	retval = nf_put_att_text(ncid, var_id, what, ldef
     +				,def)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_define_range(ncid,rmin,rmax,flag,var_id)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	real rmin,rmax,flag
	integer var_id

	real rminmax(2)
	integer retval

	rminmax(1) = rmin
	rminmax(2) = rmax

c	retval = nf_put_att_real(ncid, var_id, 'valid_range', NF_REAL
c     +				,2,rminmax)
c	call nc_handle_err(retval)

	retval = nf_put_att_real(ncid, var_id, 'valid_min', NF_REAL
     +				,1,rmin)
	call nc_handle_err(retval)

	retval = nf_put_att_real(ncid, var_id, 'valid_max', NF_REAL
     +				,1,rmax)
	call nc_handle_err(retval)

c	retval = nf_put_att_real(ncid, var_id, 'missing_value', NF_REAL
c     +				,1,flag)
c	call nc_handle_err(retval)

	retval = nf_put_att_real(ncid, var_id, '_FillValue', NF_REAL
     +				,1,flag)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_end_define(ncid)

	implicit none

	include 'netcdf.inc'

	integer ncid

	integer retval

	retval = nf_enddef(ncid)
	call nc_handle_err(retval)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_write_coords_reg(ncid,nx,ny,xlon,ylat,depth)

	implicit none

	include 'netcdf.inc'
	include 'param.h'
	include 'netcdf.h'

	integer ncid
	integer nx,ny
	real xlon(nx)
	real ylat(ny)
	real depth(nx,ny)

        real hlv(nlvdim)
        common /hlv/hlv

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
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, lat_varid, ylat)
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, dep_varid, depth)
	call nc_handle_err(retval)

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************

	subroutine nc_write_coords(ncid)

	implicit none

	include 'netcdf.inc'
	include 'param.h'
	include 'netcdf.h'

	integer ncid

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        integer nen3v(3,neldim)
        common /nen3v/nen3v
        real hlv(nlvdim)
        common /hlv/hlv
        real hkv(nkndim)
        common /hkv/hkv

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
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, lat_varid, ygv)
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval)

	retval = nf_put_var_real(ncid, dep_varid, hkv)
	call nc_handle_err(retval)

	retval = nf_put_var_int(ncid, eix_varid, nen3v)
	call nc_handle_err(retval)

	retval = nf_put_var_int(ncid, top_varid, 2)
	call nc_handle_err(retval)

c-----------------------------------------
c end of routine
c-----------------------------------------

	end

c*****************************************************************

	subroutine nc_write_time(ncid,irec,it)

	implicit none

	include 'netcdf.inc'
	include 'netcdf.h'

	integer ncid
	integer irec
	integer it

	integer retval

	retval = nf_put_vara_int(ncid, rec_varid, irec, 1, it)
	call nc_handle_err(retval)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_write_data_2d_reg(ncid,var_id,irec,nx,ny,var2d)

	implicit none

	include 'netcdf.inc'

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
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_write_data_3d_reg(ncid,var_id,irec,nlv,nx,ny,var3d)

	implicit none

	include 'netcdf.inc'

	integer ncid
	integer var_id
	integer irec
	integer nlv
	integer nx,ny
	real var3d(nlv,nx,ny)

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
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_rewrite_3d_reg(nlv,nx,ny,var3d,vnc3d)

c re-writes a 3d array to be CF compliant

	implicit none

	include 'netcdf.inc'

	integer nlv
	integer nx,ny
	real var3d(nlv,nx,ny)
	real vnc3d(nx,ny,nlv)

	integer i,j,l

	do j=1,ny
	  do i=1,nx
	    do l=1,nlv
	      vnc3d(i,j,l) = var3d(l,i,j)
	    end do
	  end do
	end do

	end

c*****************************************************************

	subroutine nc_write_data_2d(ncid,var_id,irec,nkn,var2d)

	implicit none

	include 'netcdf.inc'

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
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_write_data_3d(ncid,var_id,irec,nlv,nkn,var3d)

	implicit none

	include 'netcdf.inc'

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
	call nc_handle_err(retval)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine nc_compact_3d_reg(nlvdim,nlv,nx,ny,var_in,var_out)

	implicit none

	integer nlvdim
	integer nlv,nx,ny
	real var_in(nlvdim,nx,ny)
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

	subroutine nc_compact_3d(nlvdim,nlv,nkn,var_in,var_out)

	implicit none

	integer nlvdim
	integer nlv,nkn
	real var_in(nlvdim,nkn)
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
c*****************************************************************

	subroutine nc_close(ncid)

	implicit none

	include 'netcdf.inc'

	integer ncid

	integer retval

	retval = nf_close(ncid)
	call nc_handle_err(retval)

	end

c*****************************************************************

	subroutine nc_global(ncid,title)

c writes global conventions

	implicit none

	include 'netcdf.inc'

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
	text = 'created on ' // cdate
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

	subroutine nc_handle_err(errcode)

	implicit none

	include 'netcdf.inc'

	integer errcode

	if( errcode .eq. nf_noerr ) return

	write(6,*) 'Error: ', nf_strerror(errcode)

	stop 'error stop nc_handle_err'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c string utility routines
c*****************************************************************
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
c*****************************************************************
c date and time routines
c*****************************************************************
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

	write(6,*) date0
	write(6,*) year,month,day
	write(6,*) time0
	write(6,*) hour,min,sec

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

	write(6,*) 'cdate: ',date,time
	write(6,*) 'cdate: ',cdate

	end
	
c*****************************************************************
c*****************************************************************
c*****************************************************************
c test routines
c*****************************************************************
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

	implicit none

	integer ncid
	integer dimids_2d(2)
	integer coord_varid(5)
	integer rec_varid
	integer level_id

	integer nkn,nel,nlv
	integer it
	integer irec
	integer date0,time0
	real znv(1)
	logical next_record
	character*1 units

	date0 = 20120101
	time0 = 0
	units = 'm'

	call nc_open(ncid,nkn,nel,nlv,date0,time0)
	call nc_define_2d(ncid,'water_level',level_id)
        call nc_define_attr(ncid,'units',units,level_id)
	call nc_end_define(ncid)
	call nc_write_coords(ncid)

	irec = 0
	do while( next_record(it,znv) )

	  irec = irec + 1
	  call nc_write_time(ncid,irec,it)
	  call nc_write_data_2d(ncid,level_id,irec,nkn,znv)

	end do

	call nc_close(ncid)

	end

c******************************************************************

c	program test_nc_main
c	implicit none
c	call test_nc
c	end

c******************************************************************

