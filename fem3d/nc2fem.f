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
        character*80 var_name
        character*30 name,xcoord,ycoord,zcoord,tcoord,bathy,slmask
        character*80 varline,descrpline,factline,text,fulltext,dstring
        character*80, allocatable :: vars(:)
        character*80, allocatable :: descrps(:)
        integer ndims, nvars, ngatts, unlim
	integer dim_id,dim_len
	integer nt,nx,ny,nz
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
!        call clo_add_option('fact fact',' '
!     +			,'scale vars with these factors')

        call clo_add_option('domain limits',' '
     +			,'give domain limits and resolution')

        call clo_add_sep('additional information')
        call clo_add_com('  var is name of variable in nc file')
        call clo_add_com('  text is list of variables and descriptions'
     +				// ' for output')
        call clo_add_com('    seperate with comma and leave no space')
        call clo_add_com('    or enclose in ""')
        call clo_add_com('    example: U10,V10,PSFC or "U10 V10 PSFC"')

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
	!call clo_get_option('fact',factline)

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
c write special files or info
c-----------------------------------------------------------------

	if( bcoords ) then
	  call write_coordinates(nxdim,nydim,nx,ny,xlon,ylat)
	  if( zcoord .ne. ' ' ) call write_zcoord(nz,zdep,zzdep)
	end if

	if( bathy .ne. ' ' ) then
	  call write_2d_reg('bathy.dat',nxdim,nydim,nx,ny,bat)
	  call write_2d_grd('bathy.grd',nxdim,nydim,nx,ny,xlon,ylat,bat)
	  if( bregular ) then
	    call write_2d_fem('bathy.fem','bathymetry'
     +				,nxdim,nydim,regpar,bat)
	  end if
	  write(6,*) 'bathymetry written to files: bathy.dat bathy.grd'
	end if

	if( slmask .ne. ' ' ) then
	  call write_2d_reg('sea_land.dat',nxdim,nydim,nx,ny,slm)
	  call write_2d_grd('sea_land.grd',nxdim,nydim,nx,ny
     +						,xlon,ylat,slm)
	  if( bregular ) then
	    call write_2d_fem('sea_land.fem','index'
     +				,nxdim,nydim,regpar,slm)
	  end if
	  write(6,*) 'sea-land mask written to files: '
     +					,'sea_land.dat sea_land.grd'
	end if

c-----------------------------------------------------------------
c check variables to write
c-----------------------------------------------------------------

	n = 0
	call parse_strings(varline,n,vars)
	if( n == 0 ) then
	  !write(6,*) 'no variables given for treatment'
	  call exit(99)
	end if

	nd = n
	call parse_strings(descrpline,nd,descrps)
	write(6,*) 'variables to be handled: ',n
	call setup_description(nd,vars,descrps)

	if( n > 0 .and. .not. bregular ) then
	  stop 'error stop: cannot yet handle irregular coords'
	end if
	!stop 'error stop: cannot handle yet...'

c-----------------------------------------------------------------
c check regularity of grid -> in regpar will be the desired regular grid
c-----------------------------------------------------------------

	call check_regular_coords(nxdim,nydim,xlon,ylat
     +				,bregular,regpar_data)
	call handle_domain(dstring,bregular,regpar_data,regpar)

	if( bregular ) then
	  write(6,*) 'coordinates are regular'
	  write(6,*) regpar
	else
	  write(6,*) 'coordinates are irregular'
	  write(6,*) regpar
	end if

c-----------------------------------------------------------------
c write variables
c-----------------------------------------------------------------

	call write_variables(ncid,n,vars,descrps
     +				,nxdim,nydim,nlvdim,regpar,nrec)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	if( nrec > 0 ) then
	  write(6,*) 'total number of time records written: ',nrec
	  write(6,*) 'variables written: ',n
	  call setup_description(nd,vars,descrps)	!for information
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

	end

c*****************************************************************

	subroutine parse_variables(varline,ndim,n,vars)

	implicit none

	character*(*) varline
	integer ndim,n
	character*(*) vars(ndim)

	logical bwrite,bin_name
	integer i,istart
	character*1 c
	character(len=len(varline)+1) string

	n = 0
	string = adjustl(varline)
	bwrite = ndim > 0
	bin_name = .false.

	do i=1,len_trim(string)+1
	  c=string(i:i)
	  if( c == ' ' .or. c == ',' ) then	!white space
	    if( bin_name ) then			!end of name
	      bin_name = .false.
	      if( bwrite ) then
		if( n > ndim ) stop 'error stop parse_variables: ndim'
	        vars(n) = string(istart:i-1)
	      end if
	    end if
	  else
	    if( .not. bin_name ) then		!new name
	      bin_name = .true.
	      n = n + 1
	      istart = i
	    end if
	  end if
	end do

	end

c*****************************************************************

	subroutine write_variables(ncid,n,vars,descrps
     +					,nx,ny,nz,regpar,nrec)

	implicit none

	integer ncid
	integer n
	character*(*) vars(n)
	character*(*) descrps(n)
	integer nx,ny,nz		!size of data in nc file
	real regpar(9)
	integer nrec			!how many records written (return)

	logical bvert
	integer nit,it,var_id,i,nitt,itt
	integer iformat,nvers,nvar,ntype
	integer iunit,lmax,np,ierr,nzz
	integer datetime(2)
	integer ids(n)
	integer dims(n)
	real flags(n)
	real regpar_new(9)
	double precision atime,avalue,dtime
	character*20 line,stime
	character*80 atext,string
	character(len=len(vars)) var

	real hlv(nz)
	real hd(nx*ny)
	integer ilhkv(nx*ny)
	!real data(nx,ny,nz)
	real femdata(nz,nx,ny)
	real fem2data(nx,ny)

	integer ifileo

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

	
	np = nint(regpar_new(1)*regpar_new(2))

	if( iformat == 0 ) then
	  iunit = ifileo(30,'out.fem','unform','new')
	else
	  iunit = ifileo(30,'out.fem','form','new')
	end if

	nitt = max(1,nit)	!loop at least once - for vars without time

        do it=1,nitt
	  itt = min(nit,it)
	  call create_date_string(ncid,itt,datetime)
	  call datetime2string(datetime,stime)
          write(6,*) 'writing record: ',itt,'   ',stime

	  call fem_file_write_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype,datetime)
          call fem_file_write_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar_new)

	  do i=1,n
	    nzz = nz
	    if( dims(i) == 2 ) nzz = 1
	    call handle_data(ncid,vars(i),it,dims(i),flags(i)
     +				,nx,ny,nzz,femdata,np)

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

	nrec = nitt

	return
	end

c*****************************************************************

	subroutine handle_data(ncid,varname,it,ndims,flag
     +				,nx,ny,nz,femdata,np)

	implicit none

	integer ncid
	character*(*) varname
	integer it
	integer ndims
	real flag
	integer nx,ny,nz
	real femdata(nz,nx*ny)
	integer np

	integer ndim,nxy,k,iz
	integer nxx,nyy,nzz,nlvddi
	integer dims(10)
	real data(nx,ny,nz)
	real cdata(nx*ny,nz)
	real, save :: my_flag = -999.

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
	call compress_data(nxx,nyy,nzz,data,cdata)	!adjusts nxx,nyy,nzz
	call copy_data_to_fem(nxx,nyy,nzz,nlvddi,cdata,femdata)
	np = nxx*nyy

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

