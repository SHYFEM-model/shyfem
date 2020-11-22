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
          integer :: len
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
	  integer ncid			!number of netcdf object
	  integer nvars			!number of variables
	  integer ngatts		!number of global attributes
	  integer idunlim		!unlimited dimension
	  type(dim_item) :: ditem
	  type(att_item), allocatable :: gitems(:)
        end type nc_item

	integer, save :: nfill = 0
	integer, save :: ncids(MAX_NITEMS)
	type(nc_item), save , target :: nitems(MAX_NITEMS)

	integer retval
	integer, save :: att_name_length = 0

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

!==================================================================
!        end module ncf
!==================================================================

	subroutine ncf_open_read(file,ncid)

! reads file into ncf structure

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

        retval = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
	call ncf_handle_err(retval)

	nitem%nvars = nvars
	nitem%idunlim = idunlim

	ditem = nitem%ditem
	ditem%ndims = ndims
	allocate(ditem%len(ndims))
	allocate(ditem%name(ndims))

	do id=1,ndims
	  retval = NF_INQ_DIM(ncid,id,ditem%name(id),ditem%len(id))
	  call ncf_handle_err(retval)
	end do
	nitem%ditem = ditem

	nitem%ngatts = ngatts
	allocate(nitem%gitems(nitem%ngatts))
	do id=1,ngatts
	  call ncf_att_inf(ncid,NF_GLOBAL,id,nitem%gitems(id))
	end do

	end

!*****************************************************************

	subroutine ncf_open_write(file,ncid)

! opens file ncf for write

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

	integer ncid

	retval = nf_enddef(ncid)
	call ncf_handle_err(retval)

	end

!*****************************************************************

	subroutine ncf_close(ncid)

! closes nc file

	integer ncid

	retval = nf_close(ncid)
	call ncf_handle_err(retval)

	end

!*****************************************************************

	subroutine ncf_make_dim(ncid,ndims,len,name)

! makes dimensions for output file

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
! variables
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_var_id(ncid,name,id)

! retrieves variable id from variable name

	integer ncid
	character*(*) name
	integer id

	retval = NF_INQ_VARID(ncid,name,id)
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_var_inf(ncid,id,vitem)

! retrieves variable info from variable id and returns it in vitem

	integer ncid
	integer id
	type(var_item) :: vitem

	integer ndims
	integer len,l,i,dimid
	type(nc_item), pointer :: nitem

	retval = NF_INQ_VARNDIMS(ncid,id,ndims)
	call ncf_handle_err(retval)

	vitem%id = id
	if( allocated(vitem%dimids) ) deallocate(vitem%dimids)
	allocate(vitem%dimids(ndims))

	retval = NF_INQ_VAR(ncid,id,vitem%name,vitem%xtype
     +			,vitem%ndims,vitem%dimids,vitem%natts)
	call ncf_handle_err(retval)
	vitem%is_string = ( vitem%xtype == NF_CHAR )

	nitem => ncf_get_nitem(ncid)

	len = 1
	do i=1,ndims
	  dimid = vitem%dimids(i)
          l = nitem%ditem%len(dimid)
	  len = len * l
	end do
	vitem%len = len

	end
 
!*****************************************************************

	subroutine ncf_var_get(ncid,vitem)

! gets data for variable vitem

	integer ncid
	type(var_item) :: vitem

	integer ndims,id

	id = vitem%id

	if( vitem%is_string ) then
	  stop 'error stop: cannot handle char array'
	else
	  if( allocated(vitem%value) ) deallocate(vitem%value)
	  allocate(vitem%value(vitem%len))
	  retval = NF_GET_VAR_DOUBLE(ncid,id,vitem%value)
	end if
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_var_make(ncid,vitem)

! makes new variable vitem

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

	subroutine ncf_att_id(ncid,varid,aname,id)

! retrieves attribute id of variable varid and attribute name aitem

	integer ncid
	integer varid
	character*(*) aname
	integer id

	retval = NF_INQ_ATTID(ncid,varid,aname,id)
	call ncf_handle_err(retval)

	end
 
!*****************************************************************

	subroutine ncf_att_inf(ncid,varid,attid,aitem)

! retrieves attribute with attid of variable varid and stores info in aitem

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

	if( aitem%is_string ) then
	  if( aitem%len > NF_MAX_STRLEN ) then
	    stop 'error stop: len > NF_MAX_STRLEN'
	  end if
	  aitem%string = ' '
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

	subroutine ncf_print_dimensions(dditem)

! prints dimensions

	type(dim_item) :: dditem

	integer id

        write(6,'(a)') '   id     len   name'
        do id=1,dditem%ndims
          write(6,2000) id,dditem%len(id),'   ',trim(dditem%name(id))
 2000     format(i5,i8,a,a)
        end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_print_variable(vitem)

! prints variable

	type(var_item) :: vitem

        write(6,1000) vitem%id,vitem%xtype,vitem%ndims,vitem%natts
     +                  ,vitem%len
     +                  ,'  ',trim(vitem%name)
 1000   format(4i5,i10,a,a)

	end

!*****************************************************************

	subroutine ncf_print_variable_header

! prints variable header

        write(6,'(a)') '   id type dims atts       len  name'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ncf_print_attributes(ncid,vitem)

! prints attributes of vitem

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

	integer ncid,varid

	integer l
	character*(NF_MAX_STRLEN) :: string

	call ncf_set_att_name_length(ncid,varid)	!sets att_name_length

        string = '      id type  len  name'
	l = len_trim(string) + att_name_length - 3 + 2
	string(l:) = 'content'
        write(6,'(a)') trim(string)

	end

!*****************************************************************

	subroutine ncf_print_attribute(aitem)

! prints attribute

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

	lmax = 0
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
! error handling
!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine ncf_handle_err(errcode)

        integer errcode

        if( errcode .eq. nf_noerr ) return

        write(6,*) 'Error: ', nf_strerror(errcode)

        stop 'error stop nc_handle_err'
        end

!*****************************************************************

!==================================================================
        end module ncf
!==================================================================
