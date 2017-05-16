
	program nc2fem

	use clo

	implicit none

	integer nxdim,nydim,nlvdim,ndim

	integer ncid
        character*132 file
        character*80 var_name
        character*30 name,xcoord,ycoord,zcoord,tcoord,bathy,slmask
        character*80 varline
        character*80, allocatable :: vars(:)
        integer ndims, nvars, ngatts, unlim
	integer dim_id,dim_len
	integer nt,nx,ny,nz
	integer nit,i,it,n
	integer iwhat
	integer dims(10)
	integer year0,date0,time0
	integer iftype
	integer ifreq
	real rfact
	real regpar(7)
	real, allocatable :: xlon(:,:)
	real, allocatable :: ylat(:,:)
	real, allocatable :: zdep(:)			!mid-layer depth
	real, allocatable :: zzdep(:)			!bottom depth
	real, allocatable :: bat(:,:)
	real, allocatable :: slm(:,:)
	real, allocatable :: data(:,:,:)
	real, allocatable :: aux(:)
	double precision t
	logical bverb,bcoords,btime,binfo,bvars,bwrite
	logical binvertdepth,binvertslm
	logical bregular
	logical exists_var

c-----------------------------------------------------------------
c initialize parameters
c-----------------------------------------------------------------

	file = ' '

c-----------------------------------------------------------------
c find out what to do
c-----------------------------------------------------------------

	write(6,*) 'running nc2fem'
	call clo_init('nc2fem','nc-file','1.2')
        call clo_add_info('returns info on netcdf files')

        call clo_add_option('info',.false.,'general info on nc file')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('varinfo',.false.
     +			,'list variables contained in file')

        call clo_add_option('time',.false.
     +			,'write available time records to terminal')
        call clo_add_option('coords',.false.,'write coordinate file')
        call clo_add_option('bathy var',' '
     +			,'write bathymetry file using variable var')
        call clo_add_option('slmask var',' '
     +			,'write sea-land mask file using variable var')

        call clo_add_option('invertdepth',.false.
     +			,'invert depth values')
        call clo_add_option('invertslm',.false.
     +			,'invert slmask values (0 for sea)')

        call clo_add_option('vars var',' '
     +			,'write variables var to fem file')

	call clo_parse_options

	call clo_get_option('info',binfo)
	call clo_get_option('verbose',bverb)
	call clo_get_option('varinfo',bvars)
	call clo_get_option('time',btime)
	call clo_get_option('coords',bcoords)
	call clo_get_option('bathy',bathy)
	call clo_get_option('slmask',slmask)
	call clo_get_option('vars',varline)

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
	call get_nc_dimensions(ncid,bwrite,nt,nx,ny,nz)

	nxdim = max(1,nx)
	nydim = max(1,ny)
	nlvdim = max(1,nz)
	ndim = nxdim*nydim*nlvdim
	allocate(xlon(nxdim,nydim),ylat(nxdim,nydim))
	allocate(zdep(nlvdim),zzdep(nlvdim))
	allocate(bat(nxdim,nydim),slm(nxdim,nydim))
	allocate(data(nlvdim,nxdim,nydim),aux(ndim))

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
c handle coordinates
c-----------------------------------------------------------------

	call setup_coordinates(ncid,bverb,xcoord,ycoord
     +				,nxdim,nydim,nx,ny
     +				,xlon,ylat,aux)
	call setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,zzdep)
	call setup_bathymetry(ncid,bverb,binvertdepth,bathy
     +				,nxdim,nydim,nx,ny,bat,aux)
	call setup_sealand(ncid,bverb,slmask,nxdim,nydim,nx,ny,slm,aux)
	if( binvertslm ) slm = 1.-slm

	!call init_bound_domain(nx,ny,nz)
	!if( iwhat /= 0 ) then
	!  rfact = 1.
	!  call custom(iwhat,nx,ny,nz,rfact)
	!end if

c-----------------------------------------------------------------
c write files or info
c-----------------------------------------------------------------

	if( bcoords ) then
	  call write_coordinates(nxdim,nydim,nx,ny,xlon,ylat)
	  if( zcoord .ne. ' ' ) call write_zcoord(nz,zdep,zzdep)
	end if

	if( bathy .ne. ' ' ) then
	  call write_2d_reg('bathy.dat',nxdim,nydim,nx,ny,bat)
	  call write_2d_grd('bathy.grd',nxdim,nydim,nx,ny,xlon,ylat,bat)
	  write(6,*) 'bathymetry written to files: bathy.dat bathy.grd'
	end if

	if( slmask .ne. ' ' ) then
	  call write_2d_reg('sea_land.dat',nxdim,nydim,nx,ny,slm)
	  call write_2d_grd('sea_land.grd',nxdim,nydim,nx,ny
     +						,xlon,ylat,slm)
	  write(6,*) 'sea-land mask written to files: '
     +					,'sea_land.dat sea_land.grd'
	end if

c-----------------------------------------------------------------
c check regularity of grid
c-----------------------------------------------------------------

	call check_regular_coords(nxdim,nydim,xlon,ylat,bregular,regpar)
	if( bregular ) then
	  write(6,*) 'coordinates are regular'
	else
	  write(6,*) 'coordinates are irregular'
	end if

c-----------------------------------------------------------------
c check variables to write
c-----------------------------------------------------------------

	ndim = 0
	call parse_variables(varline,ndim,n,vars)
	if( n == 0 ) then
	  write(6,*) 'no variables given for treatment'
	  call exit(99)
	end if

	ndim = n
	allocate(vars(ndim))
	call parse_variables(varline,ndim,n,vars)
	write(6,*) 'variables to handle: ',n
	do i=1,n
	  write(6,*) i,'  ',trim(vars(i))
	end do

	if( n > 0 .and. .not. bregular ) then
	  stop 'error stop: cannot yet handle irregular coords'
	end if
	stop 'error stop: cannot handle yet...'

	!call write_variables(ncid,n,vars)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        write(6,*) 'Successful completion of routine nc2fem'
	call exit(99)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

c       else if( iwhat .eq. 3 ) then    !wrf
c	'U10','u'  'V10','v'  'PSFC','p'
c	'RAINR','r',8
c	'T2','t'  'CLOUD','c'  'RH','h'  'SWDOWN','s'

c       else if( iwhat .eq. 4 ) then    !myocean
c	'sossheig','Z'  'vosaline','S'  'votemper','T'

c	else if( iwhat .eq. 5 ) then	!ROMS/TOMS
c	'salt','S' 'temp','T'

c       else if( iwhat .eq. 6 ) then    !wrf 2
c	'U10','u'  'V10','v'
c	'MSLP','p' 'PSFC','p'
c	'RAINR','r',rfact
c	'T2','t'  'CLOUD','c',0.01  'RH','h'  'SWDOWN','s'

c       else if( iwhat .eq. 7 ) then    !ECMWF
c	'var165','u'  'var166','v'  'var151','p'  'var228','r',rfact
c	'var167','t'  'var187','c'  'var168','d'  'var176','s'

c*****************************************************************
c*****************************************************************
c*****************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine print_time_records(ncid)

c print time and date 

	implicit none

	integer ncid

	integer nit,n
	double precision atime
	character*20 line

        call nc_get_time_recs(ncid,nit)
        write(6,*) 'time records found: ',nit

        do n=1,nit
          call handle_nc_time(ncid,n,atime)
	  call dts_format_abs_time(atime,line)
          write(6,*) n,atime,line
        end do

	end

c*****************************************************************

	subroutine print_minmax_time_records(ncid)

c print time and date 

	implicit none

	integer ncid

	integer nit
	double precision atime
	character*20 line

        call nc_get_time_recs(ncid,nit)

	if( nit == 0 ) then
	  write(6,*) 'no time record found'
	else if( nit == 1 ) then
          call handle_nc_time(ncid,1,atime)
	  call dts_format_abs_time(atime,line)
	  write(6,*) 'one time record found: ',atime,line
	else
          write(6,*) 'time records found: ',nit
          call handle_nc_time(ncid,1,atime)
	  call dts_format_abs_time(atime,line)
          write(6,*) 'first time record:  ',atime,line
          call handle_nc_time(ncid,nit,atime)
	  call dts_format_abs_time(atime,line)
          write(6,*) 'last time record:   ',atime,line
	end if

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine custom(iwhat,nx,ny,nz,rfact)

	implicit none

	integer iwhat
	integer nx,ny,nz
	real rfact

	rfact = 1.

	if( iwhat .eq. 1 .or. iwhat .eq. 4 ) then
	  if( nx .eq. 677 .and. ny .eq. 253 ) then	!myocean big domain
	    call set_bound_domain(285,420,158,253,1,52)
	    write(6,*) '===================================='
	    write(6,*) '===================================='
	    write(6,*) '===================================='
	    write(6,*) 'domain for output changed'
	    write(6,*) '285,420,158,253,1,52'
	    write(6,*) '===================================='
	    write(6,*) '===================================='
	    write(6,*) '===================================='
	  end if
	else if( iwhat .eq. 6 ) then	!handel rain in new warf files
	  if( nx .eq. 90 .and. ny .eq. 90 ) then	!wrf small domain
	    rfact = 24.					!rain is mm/hour
	  else if( nx .eq. 135 .and. ny .eq. 155 ) then	!wrf big domain
	    rfact = 4.					!rain is mm/6hours
	  else
	    stop 'error stop custom: unknown domain'
	  end if
	  write(6,*) '===================================='
	  write(6,*) '===================================='
	  write(6,*) '===================================='
	  write(6,*) 'rain factor adjusted: ',rfact
	  write(6,*) '===================================='
	  write(6,*) '===================================='
	  write(6,*) '===================================='
	end if

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c this implements writing of a smaller domain
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine init_bound_domain(nx,ny,nz)

	implicit none
	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	integer nx,ny,nz

	ibound = 1

	nnx = nx
	nny = ny
	nnz = nz

	ix1 = 1
	ix2 = nx
	iy1 = 1
	iy2 = ny
	iz1 = 1
	iz2 = nz

	end

c*****************************************************************

	subroutine set_bound_domain(iix1,iix2,iiy1,iiy2,iiz1,iiz2)

	implicit none
	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	integer iix1,iix2,iiy1,iiy2,iiz1,iiz2

	if( ibound .le. 0 ) stop 'error stop set_bound_domain: no init'

	if( ix1 .lt. 1 .or. iy1 .lt. 1 .or. iz1 .lt. 1 ) goto 99
	if( ix2 .gt. nnx ) goto 99
	if( iy2 .gt. nny ) goto 99
	if( iz2 .gt. nnz ) goto 99

	ix1 = iix1
	ix2 = iix2
	iy1 = iiy1
	iy2 = iiy2
	iz1 = iiz1
	iz2 = iiz2

	write(6,*) 'set_bound_domain: domain changed'
	write(6,*) ix1,ix2,iy1,iy2,iz1,iz2

	return
   99	continue
	write(6,*) 'error in setting domain: '
	write(6,*) ix1,ix2,nnx
	write(6,*) iy1,iy2,nny
	write(6,*) iz1,iz2,nnz
	stop 'error stop set_bound_domain: indices'
	end

c*****************************************************************

	blockdata bound_domain

	implicit none
	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	data ix1,ix2,iy1,iy2,iz1,iz2 /0,0,0,0,0,0/
	data ibound,nnx,nny,nnz /0,0,0,0/

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_xycoord_names(ncid,bverb,xname,yname)

	implicit none

	integer ncid
	logical bverb
	character*(*) xname,yname

	integer nvars,iv,var_id
	character*80 name,atext

	xname = ' '
	yname = ' '

	call nc_get_var_totnum(ncid,nvars)

	do var_id=1,nvars

	  call nc_get_var_name(ncid,var_id,name)

	  call nc_get_var_attr(ncid,var_id,'standard_name',atext)
	  if( atext == 'longitude' ) call set_name(xname,name)
	  if( atext == 'latitude' ) call set_name(yname,name)

	  call nc_get_var_attr(ncid,var_id,'long_name',atext)
	  if( atext == 'longitude' ) call set_name(xname,name)
	  if( atext == 'latitude' ) call set_name(yname,name)
	  if( atext == 'Longitude of scalars' ) call set_name(xname,name)
	  if( atext == 'Latitude of scalars' ) call set_name(yname,name)

	  call nc_get_var_attr(ncid,var_id,'description',atext)
	  if( atext(1:10) == 'LONGITUDE,' ) call set_name(xname,name)
	  if( atext(1:9) == 'LATITUDE,' ) call set_name(yname,name)

	end do

	if( bverb ) write(6,*) '   xcoord: ',trim(xname)
	if( bverb ) write(6,*) '   ycoord: ',trim(yname)

	end

c*****************************************************************

	subroutine get_tcoord_name(ncid,bverb,tcoord)

	implicit none

	integer ncid
	logical bverb
	character*(*) tcoord

	integer nvars,iv,var_id
	character*80 name,atext,time_d

	tcoord = ' '

	call nc_get_var_totnum(ncid,nvars)

	do var_id=1,nvars

	  call nc_get_var_name(ncid,var_id,name)

	  call nc_get_var_attr(ncid,var_id,'standard_name',atext)
	  if( atext == 'time' ) call set_name(tcoord,name)

	  call nc_get_var_attr(ncid,var_id,'long_name',atext)
	  if( atext == 'time' ) call set_name(tcoord,name)

	  call nc_get_var_attr(ncid,var_id,'description',atext)
	  if( atext(1:13) == 'minutes since' ) call set_name(tcoord,name)

	end do

        time_d = ' '
        call nc_set_time_name(time_d,tcoord)

	if( bverb ) write(6,*) '   tcoord: ',trim(tcoord)

	end

c*****************************************************************

	subroutine get_zcoord_name(ncid,bverb,zcoord)

	implicit none

	integer ncid
	logical bverb
	character*(*) zcoord

	integer nvars,iv,var_id
	character*80 name,atext

	zcoord = ' '

	call nc_get_var_totnum(ncid,nvars)

	do var_id=1,nvars

	  call nc_get_var_name(ncid,var_id,name)

	  call nc_get_var_attr(ncid,var_id,'standard_name',atext)
	  if( atext == 'zcoord' ) call set_name(zcoord,name)
	  if( atext == 'sigma of cell face' ) call set_name(zcoord,name)

	  call nc_get_var_attr(ncid,var_id,'long_name',atext)
	  if( atext == 'zcoord' ) call set_name(zcoord,name)
	  if( atext == 'sigma of cell face' ) call set_name(zcoord,name)

	  call nc_get_var_attr(ncid,var_id,'description',atext)
	  if( atext(1:18) == 'eta values on full' ) 
     +				call set_name(zcoord,name)
	end do

	if( bverb ) write(6,*) '   zcoord: ',trim(zcoord)

	end

c*****************************************************************

	subroutine set_name(varname,newname)

	implicit none

	character*(*) varname,newname

	if( varname == ' ' ) varname = newname

	end

c*****************************************************************

	subroutine check_regular_coords(nx,ny,x,y,bregular,regpar)

	implicit none

	integer nx,ny
	real x(nx,ny)
	real y(nx,ny)
	logical bregular
	real regpar(7)

	integer ix,iy
	real xtot,ytot,eps,dx,dy,dxx,dyy

	bregular = .false.
	regpar = 0.
	dx = 0.
	dy = 0.

	if( nx > 2 ) then
	  xtot = maxval(x) - minval(x)
	  eps = 1.e-5 * xtot
	  dx = x(2,1) - x(1,1)
	  do iy=1,ny
	    do ix=2,nx
	      dxx = x(ix,iy) - x(ix-1,iy)
	      if( abs(dx-dxx) > eps ) then
	        return
	      end if
	    end do
	  end do
	end if

	if( ny > 2 ) then
	  ytot = maxval(y) - minval(y)
	  eps = 1.e-5 * ytot
	  dy = y(1,2) - y(1,1)
	  do iy=2,ny
	    do ix=1,nx
	      dyy = y(ix,iy) - y(ix,iy-1)
	      if( abs(dy-dyy) > eps ) then
	        return
	      end if
	    end do
	  end do
	end if

	bregular = .true.
	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x(1,1)
	regpar(4) = y(1,1)
	regpar(5) = dx
	regpar(6) = dy
	regpar(6) = -999.

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

	subroutine write_variables(ncid,n,vars,nx,ny,nz,data)

	implicit none

	integer ncid
	integer n
	character*(*) vars(n)
	integer nx,ny,nz
	real data(nx,ny,nz)

	logical bvert
	integer nit,it,var_id,i
	integer iformat,nvers,nvar,ntype
	integer iunit,lmax,np,ierr
	integer datetime(2)
	integer ids(n)
	integer dims(n)
	real flags(n)
	double precision atime,avalue,dtime
	character*20 line
	character*80 atext
	character(len=len(vars)) var

	real hlv(1)
	real regpar(7)

	integer ifileo

	iformat = 1
	dtime = 0.
	nvers = 0
	nvar = n
	ntype = 11
	lmax = nz
	np = nx*ny

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

	if( iformat == 0 ) then
	  iunit = ifileo(30,'out.fem','unform','new')
	else
	  iunit = ifileo(30,'out.fem','form','new')
	end if

        do it=1,nit
          call handle_nc_time(ncid,it,atime)
	  call dts_format_abs_time(atime,line)
          write(6,*) 'writing record: ',it,atime,line
	  call string2datetime(line,datetime,ierr)
	  if( ierr /= 0 ) goto 99
	  call fem_file_write_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype,datetime)
          call fem_file_write_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar)

	  !do i=1,n
	  !  call handle_data(ncid,vars(i),dims(i)
	  !end do

        end do

	return
   99	continue
	write(6,*) 'error converting date string: ',trim(line)
	stop 'error stop write_variables: date string'
	end

c*****************************************************************

