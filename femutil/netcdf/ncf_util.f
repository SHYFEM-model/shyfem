!
! easy use of netcdf files through a fortran interface
!
! info: 
!
! https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f77/
! https://www.unidata.ucar.edu/software/netcdf/examples/programs/
! http://nco.sourceforge.net/nco.html
! 
!==================================================================
        module ncf
!==================================================================

        !use netcdf

        implicit none

        include 'netcdf.inc'

	integer, parameter :: MAX_NITEMS = 10
	!integer, parameter :: MAX_NITEMS = 700000
	integer, parameter :: NF_MAX_STRLEN = 1024

	integer, parameter :: nfc_string_id = 2
	integer, parameter :: nfc_float_id = 5
	integer, parameter :: nfc_double_id = 6

        type :: dim_item
	  integer :: ndims				  !number of dimensions
          integer, allocatable :: len(:)		  !dimension
          character*(NF_MAX_NAME), allocatable :: name(:) !name
        end type dim_item

        type :: var_item
          character*(NF_MAX_NAME) :: name
	  logical :: is_string
          integer :: id
          integer :: xtype
          integer :: ndims
          integer, allocatable :: dimids(:)
          integer, allocatable :: dims(:)
          integer, allocatable :: start(:)
          integer, allocatable :: count(:)
          integer :: idtime			!time slot
          integer :: len
          integer :: rlen
          integer :: tlen
          integer :: natts
          !character(len=:), allocatable :: string
          character*(NF_MAX_STRLEN) :: string
	  double precision, allocatable :: value(:)
        end type var_item

        type :: att_item
          character*(NF_MAX_NAME) :: name
	  logical :: is_string
          integer :: id
          integer :: xtype
          integer :: len
          !character(len=:), allocatable :: string
          character*(NF_MAX_STRLEN) :: string
	  double precision, allocatable :: value(:)
        end type att_item

        type :: nc_item
          character*(NF_MAX_STRLEN) :: filename
	  integer ncid			!number of netcdf object
	  integer nvars			!number of variables
	  integer ngatts		!number of global attributes
	  integer idunlim		!unlimited dimension
	  type(dim_item) :: ditem
	  type(att_item), allocatable :: gitems(:)
        end type nc_item

	logical, parameter :: bdebug = .false.
	integer, save :: nfill = 0
	integer, save :: ncids(MAX_NITEMS)
	type(nc_item), save , target :: nitems(MAX_NITEMS)

	integer retval
	integer, save :: att_name_length = 0

        INTERFACE ncf_att_string
        MODULE PROCEDURE ncf_att_name_string,ncf_att_id_string
        END INTERFACE

        INTERFACE ncf_get_data
        MODULE PROCEDURE 	  ncf_get_data_d1
     +				, ncf_get_data_d2
     +				, ncf_get_data_d3
     +				, ncf_get_data_d4
        END INTERFACE

        INTERFACE ncf_get_record
        MODULE PROCEDURE 	  ncf_get_record_d1
     +				, ncf_get_record_d2
     +				, ncf_get_record_d3
     +				, ncf_get_record_d4
        END INTERFACE

	private bdebug

!==================================================================
	contains
!==================================================================

	subroutine ncf_insert_ncid(ncid)

	integer ncid

	integer i

	do i=1,nfill
	  if( ncids(i) == ncid ) then
	    write(6,*) i,ncids(i),ncid
	    write(6,*) 'ncid already present in list... cannot insert'
	    stop 'error stop ncf_insert_ncid: ncid present'
	  end if
	end do

	nfill = nfill + 1
	if( nfill > MAX_NITEMS ) then
	  write(6,*) nfill
	  stop 'error stop ncf_insert_ncid: nfill > MAX_NITEMS'
	end if

	ncids(nfill) = ncid

	end subroutine ncf_insert_ncid

!*****************************************************************

	subroutine ncf_delete_ncid(ncid)

	integer ncid

	integer i

	i = ncf_find_ncid(ncid)
	ncids(i) = ncids(nfill)
	nfill = nfill - 1

	end subroutine ncf_delete_ncid

!*****************************************************************

	function ncf_find_ncid(ncid)

	integer ncf_find_ncid
	integer ncid

	integer i

	do i=1,nfill
	  if( ncids(i) == ncid ) then
	    ncf_find_ncid = i
	    return
	  end if
	end do
	
	write(6,*) ncid,nfill
	write(6,*) ncids
	write(6,*) 'cannot find ncid'
	stop 'error stop ncf_find: no such ncid'

	end function ncf_find_ncid

!*****************************************************************

	function ncf_get_nitem(ncid)

	type(nc_item), pointer :: ncf_get_nitem
	integer ncid

	integer i

	i = ncf_find_ncid(ncid)
	ncf_get_nitem => nitems(i)

	end function ncf_get_nitem

!*****************************************************************
!*****************************************************************
!*****************************************************************
! get data record routines
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_get_record_d1(ncid,varid,irec,data)

	integer ncid
	integer varid
	integer irec
	double precision data

	integer, parameter :: ndim = 1
	integer i,dims(ndim)
	type(var_item) :: vitem

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get_record(ncid,irec,vitem)

	!do i=1,ndim-1
	!  dims(i) = size(data,i)
	!end do
	dims(ndim) = vitem%dims(vitem%ndims)

	call ncf_check_dims(ndim,dims,vitem)

	!data = reshape(vitem%value,vitem%dims(1:ndim-1))
	data = vitem%value(1)

	end subroutine ncf_get_record_d1

!*****************************************************************

	subroutine ncf_get_record_d2(ncid,varid,irec,data)

	integer ncid
	integer varid
	integer irec
	double precision data(:)

	integer, parameter :: ndim = 2
	integer i,dims(ndim)
	type(var_item) :: vitem

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get_record(ncid,irec,vitem)

	do i=1,ndim-1
	  dims(i) = size(data,i)
	end do
	dims(ndim) = vitem%dims(vitem%ndims)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim-1))

	end subroutine ncf_get_record_d2

!*****************************************************************

	subroutine ncf_get_record_d3(ncid,varid,irec,data)

	integer ncid
	integer varid
	integer irec
	double precision data(:,:)

	integer, parameter :: ndim = 3
	integer i,dims(ndim)
	type(var_item) :: vitem

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get_record(ncid,irec,vitem)

	do i=1,ndim-1
	  dims(i) = size(data,i)
	end do
	dims(ndim) = vitem%dims(vitem%ndims)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim-1))

	end subroutine ncf_get_record_d3

!*****************************************************************

	subroutine ncf_get_record_d4(ncid,varid,irec,data)

	integer ncid
	integer varid
	integer irec
	double precision data(:,:,:)

	integer, parameter :: ndim = 4
	integer i,dims(ndim)
	type(var_item) :: vitem

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get_record(ncid,irec,vitem)

	do i=1,ndim-1
	  dims(i) = size(data,i)
	end do
	dims(ndim) = vitem%dims(vitem%ndims)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim-1))

	end subroutine ncf_get_record_d4

!*****************************************************************
!*****************************************************************
!*****************************************************************
! get data routines
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_get_data_d1(ncid,varid,data)

	integer ncid
	integer varid
	double precision data(:)

	integer, parameter :: ndim = 1
	integer i,dims(ndim)
	type(var_item) :: vitem

	do i=1,ndim
	  dims(i) = size(data,i)
	end do

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get(ncid,vitem)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim))

	end subroutine ncf_get_data_d1

!*****************************************************************

	subroutine ncf_get_data_d2(ncid,varid,data)

	integer ncid
	integer varid
	double precision data(:,:)

	integer, parameter :: ndim = 2
	integer i,dims(ndim)
	type(var_item) :: vitem

	do i=1,ndim
	  dims(i) = size(data,i)
	end do

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get(ncid,vitem)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim))

	end subroutine ncf_get_data_d2

!*****************************************************************

	subroutine ncf_get_data_d3(ncid,varid,data)

	integer ncid
	integer varid
	double precision data(:,:,:)

	integer, parameter :: ndim = 3
	integer i,dims(ndim)
	type(var_item) :: vitem

	do i=1,ndim
	  dims(i) = size(data,i)
	end do

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get(ncid,vitem)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim))

	end subroutine ncf_get_data_d3

!*****************************************************************

	subroutine ncf_get_data_d4(ncid,varid,data)

	integer ncid
	integer varid
	double precision data(:,:,:,:)

	integer, parameter :: ndim = 4
	integer i,dims(ndim)
	type(var_item) :: vitem

	do i=1,ndim
	  dims(i) = size(data,i)
	end do

        call ncf_var_inf(ncid,varid,vitem)
        call ncf_var_get(ncid,vitem)

	call ncf_check_dims(ndim,dims,vitem)

	data = reshape(vitem%value,vitem%dims(1:ndim))

	end subroutine ncf_get_data_d4

!*****************************************************************

	subroutine ncf_check_dims(ndim,dims,vitem)

	integer ndim
	integer dims(ndim)
	type(var_item) :: vitem

	logical, parameter :: bwrite = .true.

	if( bdebug ) then
	  write(6,*) 'ndim:  ',ndim
	  write(6,*) 'dims:  ',dims
	  write(6,*) 'ndims: ',vitem%ndims
	  write(6,*) 'dims:  ',vitem%dims
	  write(6,*) 'count: ',vitem%count
	  write(6,*) 'start: ',vitem%start
	end if

	if( ndim /= vitem%ndims ) stop 'error stop: wrong ndim'
	if( any(dims/=vitem%dims) ) stop 'error stop: wrong dims'

	end subroutine ncf_check_dims

!==================================================================
!        end module ncf
!==================================================================

	subroutine ncf_open_read(file,ncid)

! reads file into ncf structure

	implicit none

	character*(*) file
	integer ncid

	integer id
        integer ndims,nvars,ngatts,idunlim
	type(dim_item) :: ditem
	type(nc_item), pointer :: nitem

	retval = nf_open(file, nf_nowrite, ncid)
	call ncf_handle_err(retval)

	call ncf_insert_ncid(ncid)

	nitem => ncf_get_nitem(ncid)

	nitem%ncid = ncid
	nitem%filename = file

        retval = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
	call ncf_handle_err(retval)

	nitem%nvars = nvars
	nitem%idunlim = idunlim

	ditem = nitem%ditem
	ditem%ndims = ndims
	if( allocated(ditem%len) ) deallocate(ditem%len)
	if( allocated(ditem%name) ) deallocate(ditem%name)
	allocate(ditem%len(ndims))
	allocate(ditem%name(ndims))

	do id=1,ndims
	  retval = NF_INQ_DIM(ncid,id,ditem%name(id),ditem%len(id))
	!GGU write(6,*) 'dim:',id,ditem%len(id),'  ',trim(ditem%name(id))
	  call ncf_handle_err(retval)
	end do
	nitem%ditem = ditem

	nitem%ngatts = ngatts
	if( allocated(nitem%gitems) ) deallocate(nitem%gitems)
	allocate(nitem%gitems(nitem%ngatts))
	do id=1,ngatts
	  call ncf_att_inf(ncid,NF_GLOBAL,id,nitem%gitems(id))
	end do

	end

!*****************************************************************

	subroutine ncf_open_write(file,ncid)

! opens file ncf for write

	implicit none

	character*(*) file
	integer ncid

	integer id

	retval = nf_create(file, nf_clobber, ncid)
	call ncf_handle_err(retval)

	call ncf_insert_ncid(ncid)

	end

!*****************************************************************

	subroutine ncf_start_data_mode(ncid)

! starts data mode to write variables and data

	implicit none

	integer ncid

	retval = nf_enddef(ncid)
	call ncf_handle_err(retval)

	end

!*****************************************************************

	subroutine ncf_close(ncid)

! closes nc file

	implicit none

	integer ncid

	retval = nf_close(ncid)
	call ncf_handle_err(retval)
	call ncf_delete_ncid(ncid)

	end

!*****************************************************************

	subroutine ncf_make_dim(ncid,ndims,len,name)

! makes dimensions for output file

	implicit none

	integer ncid
	integer ndims
	integer len(ndims)
	character*(*) name(ndims)

	integer id,dim_id

	do id=1,ndims
          retval = nf_def_dim(ncid,name(id),len(id),dim_id)
	  call ncf_handle_err(retval)
	  if( id /= dim_id ) then
	    stop 'error stop: id /= dim_id'
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! dimensions
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_var_ndims(ncid,varid,ndims)

	implicit none

	integer ncid,varid,ndims

        retval = NF_INQ_VARNDIMS(ncid,varid,ndims)
	call ncf_handle_err(retval)

	end

!*****************************************************************

        subroutine ncf_ndims(ncid,ndims)

	implicit none

        integer ncid
        integer ndims

        retval = NF_INQ_NDIMS(ncid,ndims)
        call ncf_handle_err(retval)

	end

!*****************************************************************

        subroutine ncf_dim(ncid,dimid,dimname,dimlen)

	implicit none

        integer ncid
        integer dimid
	character*(*) dimname
        integer dimlen

        character*(NF_MAX_NAME) :: nameaux

        retval = NF_INQ_DIMLEN(ncid,dimid,dimlen)
        call ncf_handle_err(retval)
        retval = NF_INQ_DIMNAME(ncid,dimid,nameaux)
        call ncf_handle_err(retval)
	dimname = nameaux

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! variables
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_nvars(ncid,nvars)

! retrieves total number of variables

	implicit none

	integer ncid
	integer nvars

	!type(nc_item), pointer :: ncf_get_nitem
	type(nc_item) :: nitem

	nitem = ncf_get_nitem(ncid)
	nvars = nitem%nvars

	end
 
!*****************************************************************

	subroutine ncf_var_id(ncid,name,varid)

! retrieves variable id from variable name

	implicit none

	integer ncid
	character*(*) name
	integer varid

	retval = NF_INQ_VARID(ncid,name,varid)
        if( ncf_has_err(retval) ) varid = 0

	end
 
!*****************************************************************

	subroutine ncf_var_name(ncid,varid,name)

! retrieves variable name from variable id

	implicit none

	integer ncid
	integer varid
	character*(*) name

        character*(NF_MAX_NAME) :: nameaux

	retval = NF_INQ_VARNAME(ncid,varid,nameaux)
        if( ncf_has_err(retval) ) nameaux = ' '
	name = nameaux

	end
 
!*****************************************************************

	subroutine ncf_var_inf(ncid,varid,vitem)

! retrieves variable info from variable id and returns it in vitem

	implicit none

	integer ncid
	integer varid
	type(var_item) :: vitem

	integer ndims
	integer len,rlen,l,i,dimid,idunlim
	type(nc_item), pointer :: nitem
        character*(NF_MAX_NAME) :: dname

	retval = NF_INQ_VARNDIMS(ncid,varid,ndims)
	call ncf_handle_err(retval)

	vitem%id = varid
	if( allocated(vitem%dimids) ) then
	  deallocate(vitem%dimids)
	  deallocate(vitem%dims)
	  deallocate(vitem%start)
	  deallocate(vitem%count)
	end if
	allocate(vitem%dimids(ndims))
	allocate(vitem%dims(ndims))
	allocate(vitem%start(ndims))
	allocate(vitem%count(ndims))

	retval = NF_INQ_VAR(ncid,varid,vitem%name,vitem%xtype
     +			,vitem%ndims,vitem%dimids,vitem%natts)
	call ncf_handle_err(retval)
	vitem%is_string = ( vitem%xtype == NF_CHAR )

	nitem => ncf_get_nitem(ncid)
	idunlim = nitem%idunlim

	len = 1
	rlen = 1
	vitem%idtime = 0
	vitem%tlen = 1
	do i=1,ndims
	  dimid = vitem%dimids(i)
          l = nitem%ditem%len(dimid)
          dname = nitem%ditem%name(dimid)
	  !GGU write(6,*) 'dimmm: ',i,dimid,l,trim(dname)
	  vitem%start(i) = 1
	  vitem%count(i) = l
	  vitem%dims(i) = l
	  len = len * l
	  if( dimid == idunlim ) then
	    vitem%idtime = i
	    vitem%tlen = l
	    vitem%count(i) = 1		!for one time record
	  else
	    rlen = rlen * l
	  end if
	end do
	vitem%len = len
	vitem%rlen = rlen

	end
 
!*****************************************************************

	subroutine ncf_var_get(ncid,vitem)

! gets data for variable vitem

	implicit none

	integer ncid
	type(var_item) :: vitem

	integer ndims,id

	id = vitem%id

	if( vitem%is_string ) then
	  stop 'error stop ncf_var_get: cannot handle char array'
	else
	  if( allocated(vitem%value) ) deallocate(vitem%value)
	  allocate(vitem%value(vitem%len))
	  retval = NF_GET_VAR_DOUBLE(ncid,id,vitem%value)
	end if
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_var_get_record(ncid,irec,vitem)

! gets record for variable vitem

	implicit none

	integer ncid
	integer irec			!time record to get
	type(var_item) :: vitem

	integer ndims,id,idtime,tlen

	id = vitem%id
	tlen = vitem%tlen
	idtime = vitem%idtime

	if( irec < 1 .or. irec > tlen ) then
	  write(6,*) 'no such time record: ',irec
	  write(6,*) 'maximum time record: ',tlen
	  stop 'error stop ncf_var_get_record: no such time'
	end if

	if( vitem%is_string ) then
	  stop 'error stop ncf_var_get: cannot handle char array'
	else
	  if( allocated(vitem%value) ) deallocate(vitem%value)
	  allocate(vitem%value(vitem%rlen))
	  vitem%start(idtime) = irec
	  vitem%count(idtime) = 1
	  retval = NF_GET_VARA_DOUBLE(ncid,id
     +			,vitem%start,vitem%count,vitem%value)
	end if
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_var_make(ncid,vitem)

! makes new variable structure vitem

	implicit none

	integer ncid
	type(var_item) :: vitem

	integer varid

	retval = NF_DEF_VAR(ncid,vitem%name,vitem%xtype
     +				,vitem%ndims,vitem%dimids,varid)
	call ncf_handle_err(retval)

	if( vitem%id == 0 ) then
	  vitem%id = varid
	else if( vitem%id /= varid ) then
	  stop 'error stop: id /= varid'
	end if

	end
 
!*****************************************************************

	subroutine ncf_var_put(ncid,vitem)

! writes variable data to file

	implicit none

	integer ncid
	type(var_item) :: vitem

	integer id,varid

	varid = vitem%id

	if( vitem%is_string ) then
	  stop 'error stop: cannot handle char array'
	else
	  retval = NF_PUT_VAR_DOUBLE(ncid,varid,vitem%value)
	  call ncf_handle_err(retval)
	end if

	end
 
!*****************************************************************
!*****************************************************************
!*****************************************************************
! attributes
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_natts(vitem,natts)

! checks if variable has atrributes

	implicit none

	type(var_item) :: vitem
	integer natts

	natts = vitem%natts

	end

!*****************************************************************

	subroutine ncf_att_id(ncid,varid,attname,attid)

! retrieves attribute id (attid) of variable varid and attribute name attname

	implicit none

	integer ncid
	integer varid
	character*(*) attname
	integer attid

	retval = NF_INQ_ATTID(ncid,varid,attname,attid)
        if( ncf_has_err(retval) ) attid = 0

	end
 
!*****************************************************************

	subroutine ncf_att_name(ncid,varid,attid,attname)

! retrieves attribute name (attname) of variable varid and attribute id (attid)

	implicit none

	integer ncid
	integer varid
	integer attid
	character*(*) attname

        character*(NF_MAX_NAME) :: nameaux

	retval = NF_INQ_ATTNAME(ncid,varid,attid,nameaux)
        if( ncf_has_err(retval) ) nameaux = ' '
	attname = nameaux

	end
 
!*****************************************************************

	subroutine ncf_att_name_string(ncid,varid,attname,string)

! retrieves attribute string of variable varid and attribute name

	implicit none

	integer ncid
	integer varid
	character*(*) attname
	character*(*) string

	integer attid
	type(att_item) :: aitem

	string = ' '

	call ncf_att_id(ncid,varid,attname,attid)
	if( attid == 0 ) return
	call ncf_att_inf(ncid,varid,attid,aitem)

	string = aitem%string

	end
 
!*****************************************************************

	subroutine ncf_att_id_string(ncid,varid,attid,string)

! retrieves attribute string of variable varid and attribute id

	implicit none

	integer ncid
	integer varid
	integer attid
	character*(*) string

	type(att_item) :: aitem

	call ncf_att_inf(ncid,varid,attid,aitem)

	string = aitem%string

	end
 
!*****************************************************************

	subroutine ncf_att_inf(ncid,varid,attid,aitem)

! retrieves attribute with attid of variable varid and stores info in aitem

	implicit none

	integer ncid
	integer varid
	integer attid
	type(att_item) :: aitem

        character*(NF_MAX_NAME) :: name

	retval = NF_INQ_ATTNAME(ncid,varid,attid,name)
	call ncf_handle_err(retval)

	retval = NF_INQ_ATT(ncid,varid,name,aitem%xtype,aitem%len)
	call ncf_handle_err(retval)
	aitem%id = attid
	aitem%name = name
	aitem%is_string = ( aitem%xtype == NF_CHAR )
	aitem%string = ' '

	if( aitem%is_string ) then
	  if( aitem%len > NF_MAX_STRLEN ) then
	    stop 'error stop: len > NF_MAX_STRLEN'
	  end if
	  retval = NF_GET_ATT_TEXT(ncid,varid,name,aitem%string)
	else
	  if( allocated(aitem%value) ) deallocate(aitem%value)
	  allocate(aitem%value(aitem%len))
	  retval = NF_GET_ATT_DOUBLE(ncid,varid,name,aitem%value)
	end if
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_att_put(ncid,varid,aitem)

! writes attribute aitem of variable varid to file

	implicit none

	integer ncid
	integer varid
	type(att_item) :: aitem

        character*(NF_MAX_NAME) :: name

	if( aitem%is_string ) then
	  retval = NF_PUT_ATT_TEXT(ncid,varid,aitem%name
     +				,aitem%len,aitem%string)
	else
	  retval = NF_PUT_ATT_DOUBLE(ncid,varid,aitem%name
     +					,aitem%xtype,aitem%len
     +					,aitem%value)
	end if
	call ncf_handle_err(retval)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! printing
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_print_dimensions(ditem)

! prints dimensions

	implicit none

	type(dim_item) :: ditem

	integer id

        write(6,'(a)') '   id     len   name'
        do id=1,ditem%ndims
          write(6,2000) id,ditem%len(id),'   ',trim(ditem%name(id))
 2000     format(i5,i8,a,a)
        end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_print_variable(vitem)

! prints variable

	implicit none

	type(var_item) :: vitem

        write(6,1000) vitem%id,vitem%xtype,vitem%ndims,vitem%natts
     +                  ,vitem%len
     +                  ,'  ',trim(vitem%name)
 1000   format(4i5,i10,a,a)

	end

!*****************************************************************

	subroutine ncf_print_variable_header

! prints variable header

	implicit none

        write(6,'(a)') '   id type dims atts       len  name'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_print_attributes(ncid,vitem)

! prints attributes of vitem

	implicit none

	integer ncid
	type(var_item) :: vitem

	integer ia,natts,varid
	type(att_item) :: aitem

	natts = vitem%natts
	varid = vitem%id

        do ia=1,natts
          call ncf_att_inf(ncid,varid,ia,aitem)
          call ncf_print_attribute(aitem)
        end do

	end

!*****************************************************************

	subroutine ncf_print_attribute_header(ncid,varid)

! prints attribute header

	implicit none

	integer ncid,varid

	integer l
	character*(NF_MAX_STRLEN) :: string

	call ncf_set_att_name_length(ncid,varid)	!sets att_name_length

	if( att_name_length < 4 ) att_name_length = 4

        string = '      id type  len  name'
	l = len_trim(string) + att_name_length - 3 + 2
	string(l:) = 'content'
        write(6,'(a)') trim(string)

	end

!*****************************************************************

	subroutine ncf_print_attribute(aitem)

! prints attribute

	implicit none

	type(att_item) :: aitem

	logical, save :: blimit = .false.
	integer l
	character*(NF_MAX_STRLEN) :: string,content

	content = aitem%string
	l = 80 - 22 - att_name_length
	if( blimit ) content = content(1:l)
	string = aitem%name
	l = att_name_length + 2

        if( aitem%is_string ) then
          write(6,3001) '   ',aitem%id,aitem%xtype,aitem%len
     +                  ,'  ',string(1:l),trim(content)
        else
          write(6,3002) '   ',aitem%id,aitem%xtype,aitem%len
     +                  ,'  ',string(1:l),aitem%value(1)
        end if

 3001         format(a,3i5,a,a,a)
 3002         format(a,3i5,a,a,d14.4)
	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_set_att_name_length(ncid,varid)

! sets attribute name length of variable varid

	implicit none

	integer ncid,varid

	integer natts
	integer ia,len,lmax
	type(att_item) :: aitem
	type(var_item) :: vitem
	type(nc_item) :: nitem

	if( varid == NF_GLOBAL ) then
	  nitem = ncf_get_nitem(ncid)
	  natts = nitem%ngatts
	else
	  call ncf_var_inf(ncid,varid,vitem)
	  natts = vitem%natts
	end if

	lmax = 0			!minimum 1 char
	do ia=1,natts
	  if( varid == NF_GLOBAL ) then
	    aitem = nitem%gitems(ia)
	  else
            call ncf_att_inf(ncid,varid,ia,aitem)
	  end if
	  len = len_trim(aitem%name)
	  lmax = max(lmax,len)
	end do

	att_name_length = lmax		!this is global

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! general information on file
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_file_info(nitem,bverbose)

	implicit none

	logical, optional :: bverbose
	type(nc_item) :: nitem

	logical battribute
	integer nvars,natts,ngatts,idunlim,id,varid
	integer ncid
	type(att_item) :: aitem
	type(var_item) :: vitem

	battribute = .false.
	if( present(bverbose) ) battribute = bverbose

        ncid = nitem%ncid
        nvars = nitem%nvars
        ngatts = nitem%ngatts
        idunlim = nitem%idunlim

        write(6,*) '========================================'
        write(6,*) 'general information on netcdf file'
        write(6,*) 'file name: ',trim(nitem%filename)
        write(6,*) 'number of variables: ',nvars
        write(6,*) 'global attributes:   ',ngatts
        write(6,*) 'unlimited dimension: ',idunlim
        write(6,*) 'dimensions: ',nitem%ditem%ndims

        call ncf_print_dimensions(nitem%ditem)

	if( ngatts > 0 ) then
          write(6,*) 'global attributes: ',ngatts
          call ncf_print_attribute_header(ncid,NF_GLOBAL)
          do id=1,ngatts
            aitem = nitem%gitems(id)
            call ncf_print_attribute(aitem)
          end do
	end if

	if( nvars > 0 ) then
          write(6,*) 'variables: ',nvars
          call ncf_print_variable_header
          do varid=1,nvars
            call ncf_var_inf(ncid,varid,vitem)
            call ncf_print_variable(vitem)
	    natts = vitem%natts
            if( battribute .and. natts > 0 ) then
              call ncf_print_attribute_header(ncid,varid)
              call ncf_print_attributes(ncid,vitem)
            end if
          end do
	end if

        write(6,*) '========================================'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! error handling
!*****************************************************************
!*****************************************************************
!*****************************************************************

        function ncf_has_err(errcode)

	implicit none

	logical ncf_has_err
        integer errcode

        ncf_has_err = ( errcode .ne. nf_noerr )

        end

!*****************************************************************

        subroutine ncf_handle_err(errcode)

	implicit none

        integer errcode

        if( errcode .eq. nf_noerr ) return

        write(6,*) 'Error: ', nf_strerror(errcode)

        stop 'error stop nc_handle_err'
        end

!*****************************************************************

!==================================================================
        end module ncf
!==================================================================
