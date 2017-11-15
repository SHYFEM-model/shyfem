!
! still to do:
!
!	handle iregular grid
!	check 3d variables
!	check resampling regular grid
!	pass depth values for hd if 3d
!	pass hlv for 3d
!	construct ilhkv
!
!*********************************************************************

	program nc2fem

	use clo

	implicit none

	integer nxdim,nydim,nlvdim,ndim

	integer ncid
        character*132 file
        character*80 var_name,files
        character*80 name,xcoord,ycoord,zcoord,tcoord,bathy,slmask
        character*80 varline,descrpline,factline,text,fulltext,dstring
        character*80, allocatable :: vars(:)
        character*80, allocatable :: descrps(:)
        character*80, allocatable :: sfacts(:)
        real, allocatable :: facts(:)
        integer ndims, nvars, ngatts, unlim
	integer dim_id,dim_len
	integer nt,nx,ny,nz,nz1
	integer nxnew,nynew
	integer nit,i,it,n,nd,nrec
	integer iwhat
	integer dims(10)
	integer year0,date0,time0
	integer iftype
	integer ifreq
	integer regexpand
	real rfact
	real regpar_data(9)
	real regpar(9)
	real, allocatable :: xlon(:,:)
	real, allocatable :: ylat(:,:)
	real, allocatable :: zdep(:)			!mid-layer depth
	real, allocatable :: zzdep(:)			!bottom depth
	real, allocatable :: hlv(:)			!bottom depth
	real, allocatable :: bat(:,:)
	real, allocatable :: slm(:,:)
	real, allocatable :: batnew(:,:)
	double precision t
	logical bverb,bcoords,btime,binfo,bvars,bwrite,bdebug
	logical binvertdepth,binvertslm,bunform,bquiet
	logical bregular
	logical exists_var

	interface
	  subroutine parse_strings(line,n,vars)
	  character*(*) line
	  integer n
	  character(len=80), allocatable :: vars(:)
	  end subroutine parse_strings
	end interface

c-----------------------------------------------------------------
c initialize parameters
c-----------------------------------------------------------------

	file = ' '

c-----------------------------------------------------------------
c find out what to do
c-----------------------------------------------------------------

	call clo_init('nc2fem','nc-file','1.2')
        call clo_add_info('converts nc (netcdf) file to fem file')

	call clo_add_sep('general options')

        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('quiet',.false.,'be as quiet as possible')
        call clo_add_option('info',.false.,'general info on nc file')
        call clo_add_option('debug',.false.,'produce debug information')
        call clo_add_option('varinfo',.false.
     +			,'list variables contained in file')

	call clo_add_sep('special variables')

        call clo_add_option('time',.false.
     +			,'write available time records to terminal')
        call clo_add_option('coords',.false.,'write coordinate file')
        call clo_add_option('bathy var',' '
     +			,'write bathymetry file using variable var')
        call clo_add_option('slmask var',' '
     +			,'write sea-land mask file using variable var')

        call clo_add_option('invertdepth',.false.
     +			,'invert depth values for bathymetry')
        call clo_add_option('invertslm',.false.
     +			,'invert slmask values (0 for sea)')
        call clo_add_option('unform',.false.
     +			,'write fem file unformatted')

	call clo_add_sep('output general variables')

        call clo_add_option('vars text',' '
     +			,'write variables given in text to out.fem')
        call clo_add_option('descrp text',' '
     +			,'use this description for variables')
        call clo_add_option('fact fact',' '
     +			,'scale vars with these factors')

        call clo_add_option('domain limits',' '
     +			,'give domain limits and resolution')
        call clo_add_option('regexpand iexp',-1,'expand regular grid')

        call clo_add_com('    iexp>0 expands iexp cells, =0 whole grid')


        call clo_add_sep('additional information')
        call clo_add_com('  var is name of variable in nc file')
        call clo_add_com('  fact is factor for multiplication of vars')
        call clo_add_com('  text is list of variables and descriptions'
     +				// ' for output')
        call clo_add_com('    seperate with comma and leave no space')
        call clo_add_com('    or enclose in ""')
        call clo_add_com('    example: U10,V10,PSFC or "U10 V10 PSFC"')
        call clo_add_com('  limits gives extension of domain')
        call clo_add_com('    for regular domain use x0,y0,x1,y1 ' //
     +				'(dx,dy are the same)')
	call clo_add_com('    for irregular domain use ' //
     +				'dx[,dy[,x0,y0,x1,y1]]')

	call clo_parse_options

	call clo_get_option('info',binfo)
	call clo_get_option('verbose',bverb)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('debug',bdebug)
	call clo_get_option('varinfo',bvars)
	call clo_get_option('time',btime)
	call clo_get_option('coords',bcoords)
	call clo_get_option('bathy',bathy)
	call clo_get_option('slmask',slmask)
	call clo_get_option('vars',varline)
	call clo_get_option('descrp',descrpline)
	call clo_get_option('domain',dstring)
	call clo_get_option('regexpand',regexpand)
	call clo_get_option('fact',factline)

	call clo_get_option('invertdepth',binvertdepth)
	call clo_get_option('invertslm',binvertslm)
	call clo_get_option('unform',bunform)

	call clo_check_files(1)

        call clo_get_file(1,file)
	write(6,*) 
        if( file == ' ' ) call clo_usage

	bwrite = bverb .or. binfo

        !call read_frequency(ifreq)

c-----------------------------------------------------------------
c open nc file and write info
c-----------------------------------------------------------------

	call nc_open_read(ncid,file)
	call global_information(ncid)

	call get_nc_dimensions(ncid,bwrite,nt,nx,ny,nz)

	nxdim = max(1,nx)
	nydim = max(1,ny)
	nlvdim = max(1,nz)
	ndim = nxdim*nydim*nlvdim
	allocate(xlon(nxdim,nydim),ylat(nxdim,nydim))
	allocate(zdep(nlvdim),zzdep(nlvdim),hlv(nlvdim))
	allocate(bat(nxdim,nydim),slm(nxdim,nydim))

	call nc_set_domain(1,nxdim,1,nydim,1,nlvdim)

	if( bwrite ) write(6,*) 'coordinates: '
	xlon = 0.
	call get_xycoord_names(ncid,bwrite,xcoord,ycoord)
	xlon = 0.
	call get_zcoord_name(ncid,bwrite,zcoord)
	call get_tcoord_name(ncid,bwrite,tcoord)
	if( nx > 0 .and. xcoord == ' ' ) then
	  write(6,*) 'nx = ',nx,'   xcoord = ',trim(xcoord)
	  stop 'error stop: dimension without variable name'
	end if
	if( ny > 0 .and. ycoord == ' ' ) then
	  write(6,*) 'ny = ',ny,'   ycoord = ',trim(ycoord)
	  stop 'error stop: dimension without variable name'
	end if
	if( nz > 0 .and. zcoord == ' ' ) then
	  write(6,*) 'nz = ',nz,'   zcoord = ',trim(zcoord)
	  stop 'error stop: dimension without variable name'
	end if
	if( nt > 0 .and. tcoord == ' ' ) then
	  write(6,*) 'nt = ',nt,'   tcoord = ',trim(tcoord)
	  stop 'error stop: dimension without variable name'
	end if

	if( bvars ) then
	  !call nc_dims_info(ncid)
	  call nc_vars_info(ncid,bverb)
	end if

	call setup_nc_time(ncid,bverb)
	if( bwrite ) call print_minmax_time_records(ncid)

	!if( binfo .or. bvars ) call exit(99)

c-----------------------------------------------------------------
c handle time
c-----------------------------------------------------------------

	if( btime ) then
	  call print_time_records(ncid)
	end if

c-----------------------------------------------------------------
c handle coordinates and special variables
c-----------------------------------------------------------------

	call setup_coordinates(ncid,bverb,xcoord,ycoord
     +				,nxdim,nydim,nx,ny
     +				,xlon,ylat)
	call setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,nz1,hlv)
	call setup_bathymetry(ncid,bverb,binvertdepth,bathy
     +				,nxdim,nydim,nx,ny,bat)
	call setup_sealand(ncid,bverb,slmask,nxdim,nydim,nx,ny,slm)
	if( binvertslm ) slm = 1.-slm

c-----------------------------------------------------------------
c check regularity of grid 
c - in regpar_data is info on original grid
c - in regpar is info on desired regular output grid
c - bregular is true if regular_data is regular
c-----------------------------------------------------------------

	call check_regular_coords(nxdim,nydim,xlon,ylat
     +				,bregular,regpar_data)
	call handle_domain(bverb,dstring,bregular,regpar_data,regpar)

	if( bregular ) then
	  write(6,*) 'original coordinates are regular'
	  if( bverb ) write(6,*) regpar
	else
	  write(6,*) 'original coordinates are irregular'
	  if( bverb ) write(6,*) regpar
	end if

c-----------------------------------------------------------------
c write special files or info
c-----------------------------------------------------------------

	if( bcoords ) then
	  call write_coordinates(nxdim,nydim,xlon,ylat)
	  call write_regular_coordinates_as_grd(regpar)
	  if( zcoord .ne. ' ' ) call write_zcoord(nz,zdep,nz1,hlv)
	end if

	if( bathy .ne. ' ' ) then
	  call write_2d_dat('bathy.dat',nxdim,nydim,bat)
	  call write_2d_grd('bathy.grd',nxdim,nydim,xlon,ylat,bat)
	  files = 'bathy.dat bathy.grd'
	  if( bregular ) then
	    call write_2d_fem('bathy.fem','bathymetry',regpar_data,bat)
	    files = trim(files) // ' bathy.fem'
	  end if
	  write(6,*) 'bathymetry written to files: ',files
	end if

	if( slmask .ne. ' ' ) then
	  call write_2d_dat('sea_land.dat',nxdim,nydim,slm)
	  call write_2d_grd('sea_land.grd',nxdim,nydim,xlon,ylat,slm)
	  files = 'sea_land.dat sea_land.grd'
	  if( bregular ) then
	    call write_2d_fem('sea_land.fem','index',regpar_data,slm)
	    files = trim(files) // ' sea_land.fem'
	  end if
	  write(6,*) 'sea-land mask written to files: ',files
	end if

c-----------------------------------------------------------------
c check variables to write
c-----------------------------------------------------------------

	n = 0
	call parse_strings(varline,n,vars)
	nd = n

	if( nd > 0 ) then
	  call parse_strings(descrpline,nd,descrps)
	  call setup_description(nd,vars,descrps)
	  call parse_strings(factline,nd,sfacts)
	  allocate(facts(nd))
	  call setup_facts(nd,sfacts,facts)
	end if

c-----------------------------------------------------------------
c set up interpolation
c-----------------------------------------------------------------

	nxnew = nint(regpar(1))
	nynew = nint(regpar(2))
	allocate(batnew(nxdim,nydim))
	batnew = -999.

	if( bregular ) then	!no interpolation - already regular
	  call prepare_no_interpol
	else
	  call prepare_interpol(nxdim,nydim,xlon,ylat,regpar)
	end if

c-----------------------------------------------------------------
c interpolate special variables
c-----------------------------------------------------------------

	if( bathy .ne. ' ' ) then
	  call handle_interpol_2d(nxdim,nydim,bat,nxnew,nynew,batnew)
	  call write_2d_fem('bathy_new.fem','bathymetry',regpar,batnew)
	  call write_2d_grd_regular('bathy_new.grd',regpar,batnew)
	end if

c-----------------------------------------------------------------
c write variables
c-----------------------------------------------------------------

	call write_variables(ncid,n,bunform,bdebug,bquiet
     +				,vars,descrps,facts
     +				,regexpand
     +				,nxdim,nydim,nlvdim
     +				,xlon,ylat
     +				,nz1,hlv
     +				,nxnew,nynew,regpar,nrec)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	if( nrec > 0 ) then
	  write(6,*) 'total number of time records written: ',nrec
	  write(6,*) 'variables written: ',n
	  write(6,*) 'output written to file out.fem'
	  call setup_description(nd,vars,descrps)	!for information
	  write(6,*) 'facts: ',facts
	end if

        write(6,*) 'Successful completion of routine nc2fem'
	call exit(99)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine parse_strings(line,n,vars)

	implicit none

	character*(*) line
	integer n
	character(len=80), allocatable :: vars(:)

	integer norig,ndim

	norig = n
	ndim = 0
	call parse_variables(line,ndim,n,vars)

	ndim = n
	if( norig > 0 ) ndim = norig
	allocate(vars(ndim))
	vars = ' '

	if( n > ndim ) then
	  write(6,*) 'too many strings given: ',n
	  write(6,*) 'expecting: ',norig
	  write(6,*) 'line: ',trim(line)
	  stop 'error stop parse_strings: too many strings'
	end if
	
	call parse_variables(line,ndim,n,vars)

	if( n == 1 ) vars = vars(1)	!set all values with only value given

	if( n > 1 .and. n /= ndim ) then
	  write(6,*) 'wrong number of strings given: ',n
	  write(6,*) 'possible numbers: ',0,1,ndim
	  write(6,*) 'line: ',trim(line)
	  stop 'error stop parse_strings: wrong number of strings'
	end if

	n = ndim

	end

c*****************************************************************

	subroutine parse_variables(varline,ndim,n,vars)

! parses line of variables etc..
! as separators only commas "," are allowed
! ""		n = 0
! "a"		n = 1
! ","		n = 2
! "a,"		n = 2
! ",b"		n = 2
! "a,b"		n = 2
! "a,b,"	n = 3

	implicit none

	character*(*) varline
	integer ndim,n
	character*(*) vars(ndim)

	logical bwrite,btoken
	integer i,istart
	character*1 c
	character(len=len(varline)+1) string

	! we have at least one blank in the string

	n = 1
	istart = 1
	btoken = .false.
	string = adjustl(varline)
	bwrite = ndim > 0
	if( bwrite ) vars = ' '

	do i=1,len_trim(string)+1
	  c=string(i:i)
	  if( c /= ' ' .and. c /= ',' ) btoken = .true.
	  if( c == ',' ) then			!comma (separator)
	    if( bwrite ) vars(n) = string(istart:i-1)
	    istart = i + 1
	    n = n + 1
	    if( bwrite .and. n > ndim ) then
	      stop 'error stop parse_variables: ndim'
	    end if
	  end if
	end do

	if( bwrite ) vars(n) = string(istart:i-1)
	if( .not. btoken ) n = 0

	end

c*****************************************************************

	subroutine write_variables(ncid,nvar,bunform,bdebug,bquiet
     +					,vars,descrps,facts
     +					,regexpand
     +					,nx,ny,nz
     +					,x,y
     +					,nz1,hlv
     +					,nxnew,nynew,regpar,nrec)

	use iso8601

	implicit none

	integer ncid
	integer nvar
	logical bunform,bdebug,bquiet
	character*(*) vars(nvar)
	character*(*) descrps(nvar)
	real facts(nvar)
	integer regexpand
	integer nx,ny,nz		!size of data in nc file
	real x(nx,ny),y(nx,ny)
	integer nz1
	real hlv(nz1)
	integer nxnew,nynew		!size of regular grid
	real regpar(9)			!regular grid to which to interpolate
	integer nrec			!how many records written (return)

	logical bvert,bexpand
	integer nit,it,var_id,i,ns
	integer iformat,nvers,ntype,ndd
	integer iunit,lmax,np,ierr,nzz,npnew
	integer datetime(2)
	integer ids(nvar)
	integer dims(nvar)
	real flags(nvar)
	real off(nvar)
	real, save :: my_flag = -999.
	real data(nx,ny,nz)
	double precision atime,avalue,dtime
	character*20 line,stime
	character*80 atext,string,aname
	character(len=len(vars)) var

	logical nc_has_var_attrib

	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)
	real, allocatable :: femdata(:,:,:)

	integer ifileo

	nrec = 0
	if( nvar == 0 ) return

	iformat = 1
	if( bunform ) iformat = 0
	dtime = 0.
	nvers = 0
	ntype = 11
	lmax = nz
	np = nx*ny
	npnew = nxnew*nynew
	string = 'unknown'
	off = 0.
	bexpand = ( regexpand > -1 )
	!hlv = 0.

	allocate(hd(npnew),ilhkv(npnew))
	allocate(femdata(nz,nxnew,nynew))
	hd = -999.
	ilhkv = lmax
	flags = my_flag
	off = 0.

	do i=1,nvar
	  var = vars(i)
	  call nc_get_var_id(ncid,var,var_id)
	  if( var_id == 0 ) then
	    write(6,*) 'no such variable: ',trim(var)
	    stop 'error stop write_variables: no such variable'
	  end if
	  ids(i) = var_id
	  call nc_var_info(ncid,var_id,.true.)
	  aname = '_FillValue'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    flags(i) = avalue
	  end if
	  aname = 'scale_factor'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    facts(i) = facts(i) * avalue
	  end if
	  aname = 'add_offset'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    off(i) = avalue
	  end if
	  dims(i) = 2
	  call nc_has_vertical_dimension(ncid,var,bvert)
	  if( bvert ) dims(i) = 3
	end do

        call nc_get_time_recs(ncid,nit)
        write(6,*) 'time records found: ',nit

	ndd = 0
	do i=1,nvar
	  if( ndd == 0 ) ndd = dims(i)
	  if( ndd /= dims(i) ) then
	    write(6,*) 'mixing 2d and 3d variables... not possible'
	    write(6,*) dims
	    stop 'error stop write_variables: mixing 2d and 3d'
	  end if
	end do

	lmax = nz1
	if( ndd == 2 ) lmax = 1
	np = nxnew*nynew

	if( iformat == 0 ) then
	  iunit = ifileo(30,'out.fem','unform','new')
	else
	  iunit = ifileo(30,'out.fem','form','new')
	end if

	ns = 1
	ns = min(ns,nit)	!loop at least once - for vars without time
	if( bdebug ) nit = min(5,nit)

        do it=ns,nit
	  call create_date_string(ncid,it,datetime)
	  call datetime2string(datetime,stime)
          if( .not. bquiet ) then
	    write(6,*) 'writing record: ',it,'   ',stime
	  end if

	  call fem_file_write_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype,datetime)
          call fem_file_write_2header(iformat,iunit,ntype,lmax
     +                  	,hlv,regpar(1:7))

	  do i=1,nvar
	    nzz = nz
	    if( dims(i) == 2 ) nzz = 1
	    call handle_data(ncid,bdebug,vars(i),it,dims(i),flags(i)
     +				,nx,ny,nzz
     +				,x,y
     +				,nxnew,nynew,regpar,ilhkv
     +				,data,femdata,np)

	    if( bexpand ) then
	      call reg_expand_3d(nz,nxnew,nynew,lmax,regexpand
     +					,my_flag,femdata)
	      call adjust_reg_vertical(nz,nxnew,nynew,my_flag
     +					,femdata,ilhkv)
	    end if

	    if( facts(i) /= 1. ) then
	      where( femdata /= my_flag ) femdata = femdata * facts(i)
	      where( data /= my_flag ) data = data * facts(i)
	    end if
	    if( off(i) /= 0. ) then
	      where( femdata /= my_flag ) femdata = femdata + off(i)
	      where( data /= my_flag ) data = data + off(i)
	    end if

	    lmax = nzz
	    string = descrps(i)
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,femdata)

	    if( bdebug ) then
	      call write_grd_var(vars(i),it,nx,ny,nzz,nxnew,nynew
     +			,regpar,x,y,data,femdata)
	    end if
	  end do

        end do

	close(iunit)

	nrec = nit - ns + 1
	!write(6,*) 'Total number of records written: ',nrec,nit,ns

	deallocate(femdata)
	deallocate(hd)
	deallocate(ilhkv)

	end

c*****************************************************************

	subroutine handle_data(ncid,bdebug,varname,it,ndims,flag
     +				,nx,ny,nz
     +				,x,y
     +				,nxnew,nynew,regpar,ilhkv
     +				,data,femdata,np)

	implicit none

	integer ncid
	logical bdebug
	character*(*) varname
	integer it
	integer ndims
	real flag
	integer nx,ny,nz
	real x(nx,ny),y(nx,ny)
	integer nxnew,nynew
	real regpar(9)
	integer ilhkv(nxnew*nynew)
	real data(nx,ny,nz)
	real femdata(nz,nxnew*nynew)
	integer np

	logical debug
	integer ndim,nxy,k,iz
	integer nxx,nyy,nzz,nlvddi
	integer dims(10)
	real data2d(nx,ny)
	real femdata2d(nxnew*nynew)
	real cdata(nx*ny,nz)
	real valnew(nxnew*nynew)
	real, save :: my_flag = -999.
	character*80 file,filename

	logical must_interpol

	debug = bdebug

	nxy = nx*ny
	ndim = nx*ny*nz
        call nc_get_var_data(ncid,varname,it,ndim,ndims,dims,data)

	if( ndims < 2 .or. ndims > 3 ) then
	  write(6,*) 'error in dimensions: ndims = ',ndims
	  stop 'error stop handle_data: ndims'
	end if

	if( nx /= dims(1) .or. ny /= dims(2) ) then
	  write(6,*) 'error in dimensions: '
	  write(6,*) 'nx,ny given: ',nx,ny
	  write(6,*) 'nx,ny read : ',dims(1),dims(2)
	  stop 'error stop handle_data: dimensions x/y'
	end if

	if( ndims == 3 .and. nz /= dims(3) .or.
     +			ndims == 2 .and. nz /= 1 ) then
	  write(6,*) 'error in dimensions: '
	  write(6,*) 'ndims: ',ndims
	  write(6,*) 'nz given: ',nz
	  write(6,*) 'nz read : ',dims(3)
	  stop 'error stop handle_data: dimensions z'
	end if

	!if( ndims == 3 .and. nz /= 1 ) then
	!  write(6,*) 'ndims,nz: ',ndims,nz
	!  stop 'error stop handle_data: no 3d yet'
	!end if

	where( data == flag ) data = my_flag

	nxx = nx
	nyy = ny
	nzz = nz
	nlvddi = nz

	if( must_interpol() ) then
	  np = nxnew*nynew
	  if( nz == 1 ) then
	    call do_interpol_2d(nx,ny,data,nxnew,nynew,valnew)
	  else
	    call do_interpol_3d(nx,ny,nz,data,nxnew,nynew,valnew)
	  end if
	  call copy_data_to_fem(nxnew,nynew,nzz,nlvddi,valnew,femdata)
	else
	  call compress_data(nxx,nyy,nzz,data,cdata)	!adjusts nxx,nyy,nzz
	  call copy_data_to_fem(nxx,nyy,nzz,nlvddi,cdata,femdata)
	  np = nxx*nyy
	end if
	call make_ilhkv(np,nlvddi,my_flag,femdata,ilhkv)

	end

c*****************************************************************

	subroutine write_grd_var(varname,it,nx,ny,nz,nxnew,nynew
     +			,regpar,x,y,data,femdata)

	implicit none

	character*(*) varname
	integer it
	integer nx,ny,nz
	integer nxnew,nynew
	real regpar(9)
	real x(nx,ny),y(nx,ny)
	real data(nx,ny,nz)
	real femdata(nz,nxnew*nynew)

	real, allocatable :: data2d(:,:)
	real, allocatable :: femdata2d(:)
	character*80 filename,file

	!write(6,*) 'nx,ny: ',nx,ny,nxnew,nynew
	allocate(data2d(nx,ny),femdata2d(nxnew*nynew))

	data2d(:,:) = data(:,:,1)
	femdata2d(:) = femdata(1,:)

	call make_filename(varname,it,filename)
	file=trim(filename)//'_orig.grd'
	call write_2d_grd(file,nx,ny,x,y,data2d)
	file=trim(filename)//'_intp.grd'
	call write_2d_grd_regular(file,regpar,femdata2d)

	end

c*****************************************************************

	subroutine setup_description(nd,vars,descrps)

	use shyfem_strings

	implicit none

	integer nd
        character*80 :: vars(nd)
        character*80 :: descrps(nd)

	integer i,idir
	character*80 text,fulltext
	character*1, save :: post(0:3) = (/' ','x','y','z'/)

	logical has_direction

	idir = 0

	do i=1,nd
	  text = descrps(i)
	  call strings_get_full_name(text,fulltext)
	  if( text /= ' ' .and. fulltext == ' ' ) then
	    write(6,*) '*** cannot find description for: ',trim(text)
	  end if
	  if( has_direction(fulltext) ) then
	    idir = mod(idir+1,4)
	    fulltext = trim(fulltext) // ' - ' // post(idir)
	  end if
	  descrps(i) = fulltext
	end do

	do i=1,nd
	  write(6,'(i5,a,a30,a30)') i,'  ',trim(vars(i)),trim(descrps(i))
	end do

	end

c*****************************************************************

	subroutine setup_facts(n,sfacts,facts)

	implicit none

	integer n
	character*(*) sfacts(n)
	real facts(n)

	integer i,ianz
	real f(1)
	character*80 string

	integer iscanf

	facts = 1.

	do i=1,n
	  string = sfacts(i)
	  if( string == ' ' ) cycle
	  ianz = iscanf(string,f,1)
	  if( ianz /= 1 ) then
	    write(6,*) i,ianz,'  ',string
	    stop 'error stop setup_facts: parse error'
	  end if
	  facts(i) = f(1)
	end do

	end

c*****************************************************************

	subroutine global_information(ncid)

	implicit none

	integer ncid

	integer iftype
	character*80 atext

	call global_info(ncid,'TITLE')
	call global_info(ncid,'source')
	call global_info(ncid,'type')
	!call global_info(ncid,'CDO')
	call global_info(ncid,'institution')

	end

c*****************************************************************

	subroutine global_info(ncid,what)

	implicit none

	integer ncid
	character*(*) what

	character*80 atext

        call nc_get_global_attr(ncid,what,atext)
	if( atext == ' ' ) return

	write(6,*) trim(what),': ',trim(atext)

	end

c*****************************************************************

	subroutine make_filename(varname,it,filename)

	implicit none

	character*(*) varname,filename
	integer it

	integer i
	character*80 sit

	write(sit,'(a1,i5)') '_',it
	do i=1,6
	  if( sit(i:i) == ' ' ) sit(i:i) = '0'
	end do

	filename = 'debug_' // trim(varname) // trim(sit)

	end

c*****************************************************************

