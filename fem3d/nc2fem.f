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
        character*30 name,xcoord,ycoord,zcoord,tcoord,bathy,slmask
        character*80 varline,descrpline,factline,text,fulltext,dstring
        character*80, allocatable :: vars(:)
        character*80, allocatable :: descrps(:)
        character*80, allocatable :: sfacts(:)
        real, allocatable :: facts(:)
        integer ndims, nvars, ngatts, unlim
	integer dim_id,dim_len
	integer nt,nx,ny,nz
	integer nxnew,nynew
	integer nit,i,it,n,nd,nrec
	integer iwhat
	integer dims(10)
	integer year0,date0,time0
	integer iftype
	integer ifreq
	real rfact
	real regpar_data(9)
	real regpar(9)
	real, allocatable :: xlon(:,:)
	real, allocatable :: ylat(:,:)
	real, allocatable :: zdep(:)			!mid-layer depth
	real, allocatable :: zzdep(:)			!bottom depth
	real, allocatable :: bat(:,:)
	real, allocatable :: slm(:,:)
	real, allocatable :: batnew(:,:)
	double precision t
	logical bverb,bcoords,btime,binfo,bvars,bwrite
	logical binvertdepth,binvertslm
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
        call clo_add_info('converts nc (netcdf) files to fem files')

	call clo_add_sep('general options')

        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('info',.false.,'general info on nc file')
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

	call clo_add_sep('output general variables')

        call clo_add_option('vars text',' '
     +			,'write variables given in text to out.fem')
        call clo_add_option('descrp text',' '
     +			,'use this description for variables')
        call clo_add_option('fact fact',' '
     +			,'scale vars with these factors')

        call clo_add_option('domain limits',' '
     +			,'give domain limits and resolution')

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
	call clo_get_option('varinfo',bvars)
	call clo_get_option('time',btime)
	call clo_get_option('coords',bcoords)
	call clo_get_option('bathy',bathy)
	call clo_get_option('slmask',slmask)
	call clo_get_option('vars',varline)
	call clo_get_option('descrp',descrpline)
	call clo_get_option('domain',dstring)
	call clo_get_option('fact',factline)

	call clo_get_option('invertdepth',binvertdepth)
	call clo_get_option('invertslm',binvertslm)

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
	allocate(zdep(nlvdim),zzdep(nlvdim))
	allocate(bat(nxdim,nydim),slm(nxdim,nydim))

	call nc_set_domain(1,nxdim,1,nydim,1,nlvdim)

	if( bwrite ) write(6,*) 'coordinates: '
	call get_xycoord_names(ncid,bwrite,xcoord,ycoord)
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
	call setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,zzdep)
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
	call handle_domain(dstring,bregular,regpar_data,regpar)

	if( bregular ) then
	  write(6,*) 'original coordinates are regular'
	  write(6,*) regpar
	else
	  write(6,*) 'original coordinates are irregular'
	  write(6,*) regpar
	end if

c-----------------------------------------------------------------
c write special files or info
c-----------------------------------------------------------------

	if( bcoords ) then
	  call write_coordinates(nxdim,nydim,xlon,ylat)
	  call write_regular_coordinates_as_grd(regpar)
	  if( zcoord .ne. ' ' ) call write_zcoord(nz,zdep,zzdep)
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

	call write_variables(ncid,n,vars,descrps,facts
     +				,nxdim,nydim,nlvdim
     +				,xlon,ylat
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

	subroutine write_variables(ncid,n,vars,descrps,facts
     +					,nx,ny,nz
     +					,x,y
     +					,nxnew,nynew,regpar,nrec)

	implicit none

	integer ncid
	integer n
	character*(*) vars(n)
	character*(*) descrps(n)
	real facts(n)
	integer nx,ny,nz		!size of data in nc file
	real x(nx,ny),y(nx,ny)
	integer nxnew,nynew		!size of regular grid
	real regpar(9)			!regular grid to which to interpolate
	integer nrec			!how many records written (return)

	logical bvert
	integer nit,it,var_id,i,ns
	integer iformat,nvers,nvar,ntype
	integer iunit,lmax,np,ierr,nzz
	integer datetime(2)
	integer ids(n)
	integer dims(n)
	real flags(n)
	real, save :: my_flag = -999.
	double precision atime,avalue,dtime
	character*20 line,stime
	character*80 atext,string
	character(len=len(vars)) var

	real hlv(nz)
	real hd(nx*ny)
	integer ilhkv(nx*ny)
	real femdata(nz,nxnew,nynew)
	real fem2data(nx,ny)

	integer ifileo

	if( n == 0 ) return

	iformat = 1
	dtime = 0.
	nvers = 0
	nvar = n
	ntype = 11
	lmax = nz
	np = nx*ny
	string = 'unknown'
	ilhkv = lmax
	hd = -999.
	hlv = 0.

	do i=1,n
	  var = vars(i)
	  call nc_get_var_id(ncid,var,var_id)
	  if( var_id == 0 ) then
	    write(6,*) 'no such variable: ',trim(var)
	    stop 'error stop write_variables: no such variable'
	  end if
	  ids(i) = var_id
	  call nc_var_info(ncid,var_id,.true.)
	  call nc_get_var_attrib(ncid,var_id,'_FillValue',atext,avalue)
	  flags(i) = avalue
	  dims(i) = 2
	  call nc_has_vertical_dimension(ncid,var,bvert)
	  if( bvert ) dims(i) = 3
	end do

        call nc_get_time_recs(ncid,nit)
        write(6,*) 'time records found: ',nit

	do i=1,n
	  if( dims(i) > 2 ) then
	    write(6,*) 'cannot handle 3d variable ',trim(vars(i))
	  end if
	end do

	np = nxnew*nynew

	if( iformat == 0 ) then
	  iunit = ifileo(30,'out.fem','unform','new')
	else
	  iunit = ifileo(30,'out.fem','form','new')
	end if

	ns = 1
	ns = min(ns,nit)	!loop at least once - for vars without time

        do it=ns,nit
	  call create_date_string(ncid,it,datetime)
	  call datetime2string(datetime,stime)
          write(6,*) 'writing record: ',it,'   ',stime

	  call fem_file_write_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype,datetime)
          call fem_file_write_2header(iformat,iunit,ntype,lmax
     +                  	,hlv,regpar(1:7))

	  do i=1,n
	    nzz = nz
	    if( dims(i) == 2 ) nzz = 1
	    call handle_data(ncid,vars(i),it,dims(i),flags(i)
     +				,nx,ny,nzz
     +				,x,y
     +				,nxnew,nynew,regpar,femdata,np)

	    if( facts(i) /= 1. ) then
	      where( femdata /= my_flag ) femdata = femdata * facts(i)
	    end if

	    lmax = nzz
	    string = descrps(i)
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,femdata)

	  end do

        end do

	close(iunit)

	nrec = nit - ns + 1

	end

c*****************************************************************

	subroutine handle_data(ncid,varname,it,ndims,flag
     +				,nx,ny,nz
     +				,x,y
     +				,nxnew,nynew,regpar,femdata,np)

	implicit none

	integer ncid
	character*(*) varname
	integer it
	integer ndims
	real flag
	integer nx,ny,nz
	real x(nx,ny),y(nx,ny)
	integer nxnew,nynew
	real regpar(9)
	real femdata(nz,nxnew*nynew)
	integer np

	logical debug
	integer ndim,nxy,k,iz
	integer nxx,nyy,nzz,nlvddi
	integer dims(10)
	real data(nx,ny,nz)
	real cdata(nx*ny,nz)
	real valnew(nxnew*nynew)
	real, save :: my_flag = -999.
	character*80 file,filename

	logical must_interpol

	debug = .true.

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

	if( ndims == 3 .and. nz /= 1 ) then
	  write(6,*) 'ndims,nz: ',ndims,nz
	  stop 'error stop handle_data: no 3d yet'
	end if

	where( data == flag ) data = my_flag

	nxx = nx
	nyy = ny
	nzz = nz
	nlvddi = nz

	if( must_interpol() ) then
	  call do_interpol_2d(nx,ny,data,nxnew,nynew,valnew)
	  call copy_data_to_fem(nxnew,nynew,nzz,nlvddi,valnew,femdata)
	  np = nxnew*nynew
	else
	  call compress_data(nxx,nyy,nzz,data,cdata)	!adjusts nxx,nyy,nzz
	  call copy_data_to_fem(nxx,nyy,nzz,nlvddi,cdata,femdata)
	  np = nxx*nyy
	end if

	if( nz > 1 .or. .not. debug ) return

	call make_filename(varname,it,filename)
	file=trim(filename)//'_orig.grd'
	call write_2d_grd(file,nx,ny,x,y,data)
	file=trim(filename)//'_intp.grd'
	call write_2d_grd_regular(file,regpar,femdata)

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

