
!---------------------------------------------------------------------
!
!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.
!
!     This is an example which reads some surface pressure and
!     temperatures. The data file read by this program is produced
!     comapnion program sfc_pres_temp_wr.f. It is intended to illustrate
!     the use of the netCDF fortran 77 API.
!
!     This program is part of the netCDF tutorial:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
!
!     Full documentation of the netCDF Fortran 77 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77
!
!     $Id: sfc_pres_temp_rd.f,v 1.8 2007/01/24 19:45:09 russ Exp $
!
!---------------------------------------------------------------------

      program gr

	use ncf

      implicit none

!     This is the name of the data file we will read.
      character*(80) file
      character*(*) FILE_NAME
      parameter (FILE_NAME='sfc_pres_temp.nc')
      integer ncid

!     We are reading 2D data, a 6 x 12 lat-lon grid.
      integer NDIMS
      parameter (NDIMS=2)
      integer NLATS, NLONS
      parameter (NLATS = 6, NLONS = 12)
      character*(*) LAT_NAME, LON_NAME
      parameter (LAT_NAME='latitude', LON_NAME='longitude')
      integer lat_dimid, lon_dimid

!     For the lat lon coordinate netCDF variables.
      real lats(NLATS), lons(NLONS)
      integer lat_varid, lon_varid

!     We will read surface temperature and pressure fields. 
      character*(*) PRES_NAME, TEMP_NAME
      parameter (PRES_NAME='pressure')
      parameter (TEMP_NAME='temperature')
      integer pres_varid, temp_varid
      integer dimids(NDIMS)

!     To check the units attributes.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) PRES_UNITS, TEMP_UNITS, LAT_UNITS, LON_UNITS
      parameter (PRES_UNITS = 'hPa', TEMP_UNITS = 'celsius')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')
      integer MAX_ATT_LEN
      parameter (MAX_ATT_LEN = 80)
      character*(MAX_ATT_LEN) pres_units_in, temp_units_in
      character*(MAX_ATT_LEN) lat_units_in, lon_units_in
      integer att_len

!     Read the data into these arrays.
      real pres_in(NLONS, NLATS), temp_in(NLONS, NLATS)

!     These are used to calculate the values we expect to find.
      real START_LAT, START_LON
      parameter (START_LAT = 25.0, START_LON = -125.0)
      real SAMPLE_PRESSURE
      parameter (SAMPLE_PRESSURE = 900.0)
      real SAMPLE_TEMP
      parameter (SAMPLE_TEMP = 9.0)

!     We will learn about the data file and store results in these
!     program variables.
      integer ndims_in, nvars_in, ngatts_in, idunlim
      integer nc,ia,varid

!     Loop indices
      integer lat, lon, id, ncid_out

	logical, save :: bvarelab = .false.
	logical, save :: battelab = .false.

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(nc_item), pointer :: nitem_in,nitem_out

!---------------------------------------------------------------------

        nc = command_argument_count()
	if( nc /= 1 ) then
	  write(6,*) 'Usage: gr file'
	  stop 'error stop: need one file'
	end if
        call get_command_argument(1,file)

!---------------------------------------------------------------------

	call ncf_open_read(file,ncid)

	nitem_in => ncf_get_nitem(ncid)

	ngatts_in = nitem_in%ngatts
	idunlim = nitem_in%idunlim
	write(6,*) 'global attributes: ',ngatts_in
	write(6,*) 'unlimited dimension : ',idunlim
	write(6,*) 'dimensions: ',nitem_in%ditem%ndims
	call ncf_print_dimensions(nitem_in%ditem)

	write(6,*) 'global attributes: ',ngatts_in
	call ncf_print_attribute_header(ncid,NF_GLOBAL)
	do id=1,ngatts_in
	  aitem = nitem_in%gitems(id)
	  call ncf_print_attribute(aitem)
	end do

!---------------------------------------------------------------------

	call ncf_open_write('out.nc',ncid_out)

	nitem_out => ncf_get_nitem(ncid_out)

	ditem = nitem_in%ditem
	call ncf_make_dim(ncid_out,ditem%ndims,ditem%len,ditem%name)

!---------------------------------------------------------------------

	write(6,*) 'preparing attributes for variables'
	nvars_in = nitem_in%nvars
	do varid=1,nvars_in
	  call ncf_var_inf(ncid,varid,vitem)
	  call ncf_var_make(ncid_out,vitem)
	  call ncf_print_variable(vitem)
	  do ia=1,vitem%natts
	    call ncf_att_inf(ncid,varid,ia,aitem)
	    if( battelab) call att_elab(varid,aitem)
	    call ncf_att_put(ncid_out,varid,aitem)
	    call ncf_print_attribute(aitem)
	  end do
	end do

!---------------------------------------------------------------------

	write(6,*) 'writing global attributes'
	varid = NF_GLOBAL
	do id=1,ngatts_in
	  aitem = nitem_in%gitems(id)
	  call ncf_att_put(ncid_out,varid,aitem)
	end do

	call ncf_start_data_mode(ncid_out)

	write(6,*) 'writing variables'
	do varid=1,nvars_in
	  call ncf_var_inf(ncid,varid,vitem)
	  call ncf_var_get(ncid,vitem)
	  if( bvarelab) call var_elab(varid,vitem)
	  call ncf_var_put(ncid_out,vitem)
	end do

!---------------------------------------------------------------------

	write(6,*) 'closing files'
	call ncf_close(ncid)
	call ncf_close(ncid_out)

!---------------------------------------------------------------------

	end

!*********************************************************************

	subroutine att_elab(varid,aitem)

! changes attributes for special varid

	use ncf

	implicit none

	integer varid
	type(att_item) :: aitem

	character*80 name

	if( varid > 2 ) return

	name = aitem%name

	if( name == 'La1' .or. name == 'La2' ) then
	  write(6,*) '+++ changing attribute: ',varid,'  ',trim(name)
	  aitem%value = aitem%value + 10.
	end if

	end

!*********************************************************************

	subroutine var_elab(varid,vitem)

! changes variable for special varid

	use ncf

	implicit none

	integer varid
	type(var_item) :: vitem

	character*80 name

	if( varid /= 1 ) return

	write(6,*) '+++ changing variable: ',varid,'  ',trim(vitem%name)
	vitem%value = vitem%value + 10.

	end

!*********************************************************************

