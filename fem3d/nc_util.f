!
! convert nc files to fem files: general utilities
!
! contents :
!
!
! revision log :
!
! 03.07.2018    ggu     revision control introduced
!
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine handle_data_3d(ncid,name,pre,it,n
     +				,ndim,aux,fact)

c handles the data and writes it to file

	implicit none

	integer ncid
	character*(*) name
	character*(*) pre
	integer it			!fem time in seconds
	integer n			!number of time record
	integer ndim
	real aux(ndim)
	real fact

	integer nx,ny,nz
	integer ndims
	integer dims(10)
	integer i
	double precision time
	character*30 new

	call make_name(pre,it,new)

	write(6,*) 'writing file: ',new

	ndims = 3
	call nc_get_var_data(ncid,name,n,ndim,ndims,dims,aux)
	nx = dims(1)
	ny = dims(2)
	nz = dims(3)

	call compress_data(nx,ny,nz,aux,aux)
	call scale_val(nx*ny*nz,aux,fact)

	open(1,file=new,form='formatted',status='unknown')
	write(1,*) it
	write(1,*) nx,ny,nz
	write(1,'((5e15.7))') (aux(i),i=1,nx*ny*nz)
	close(1)

	end

c*****************************************************************

	subroutine handle_data_2d(ncid,name,pre,it,n
     +				,ndim,aux,fact)

c handles the data and writes it to file

	implicit none

	integer ncid
	character*(*) name
	character*(*) pre
	integer it			!fem time in seconds
	integer n			!number of time record
	integer ndim
	real aux(ndim)
	real fact

	integer nx,ny,nz
	integer ndims
	integer dims(10)
	integer i
	double precision time
	character*30 new

	call make_name(pre,it,new)

	write(6,*) 'writing file: ',new

	ndims = 2
	call nc_get_var_data(ncid,name,n,ndim,ndims,dims,aux)
	nx = dims(1)
	ny = dims(2)

	nz = 1
	call compress_data(nx,ny,nz,aux,aux)
	call scale_val(nx*ny,aux,fact)

	open(1,file=new,form='formatted',status='unknown')
	write(1,*) it
	write(1,*) nx,ny,0
	write(1,'((5e15.7))') (aux(i),i=1,nx*ny)
	close(1)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine compress_data(nx,ny,nz,rin,rout)

	implicit none

	integer nx,ny,nz
	real rin(nx,ny,nz)
	real rout(nx*ny*nz)

	integer ip,ix,iy,iz
        integer ix1,ix2,iy1,iy2,iz1,iz2

	ip = 0
	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	if( nz .eq. 1 ) then
	  iz = 1
	  do iy=iy1,iy2
	    do ix=ix1,ix2
	      ip = ip + 1
	      rout(ip) = rin(ix,iy,iz)
	    end do
	  end do
	else
	  do iz=iz1,iz2
	    do iy=iy1,iy2
	      do ix=ix1,ix2
	        ip = ip + 1
	        rout(ip) = rin(ix,iy,iz)
	      end do
	    end do
	  end do
	  nz = iz2-iz1+1
	end if

	nx = ix2-ix1+1
	ny = iy2-iy1+1

	!write(6,*) 'compresssssssss ',nx,ny,nz

	end

c*****************************************************************

	subroutine copy_data_to_fem(nxx,nyy,nzz,nlvddi,data,femdata)

	implicit none

	integer nxx,nyy,nzz
	integer nlvddi
	real data(nxx*nyy,nzz)
	real femdata(nlvddi,nxx*nyy)

	integer nxy,k,iz

	nxy = nxx*nyy

        do k=1,nxy
          do iz=1,nzz
            femdata(iz,k) = data(k,iz)
          end do
        end do

	end

c*****************************************************************

	subroutine make_ilhkv(np,nlvddi,flag,femdata,ilhkv)

	implicit none

	integer np
	integer nlvddi
	real flag
	real femdata(nlvddi,np)
	integer ilhkv(np)

	integer k,l

        do k=1,np
          do l=1,nlvddi
	    if( femdata(l,k) == flag ) exit
          end do
	  ilhkv(k) = l-1
        end do

	end

c*****************************************************************

	subroutine scale_val(n,val,fact)

	implicit none

	integer n
	real val(n)
	real fact

	integer i

	if( fact .eq. 1. ) return

	do i=1,n
	  val(i) = fact*val(i)
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************


c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine check_var(ncid,name,ierror)

c checks if name is available

	implicit none

	integer ncid
	character*(*) name
	integer ierror

	integer var_id

	call nc_get_var_id(ncid,name,var_id)

	if( var_id .le. 0 ) then
	  ierror = ierror + 1
	  write(6,*) '*** Cannot find variable: ',name
	end if

	end

c*****************************************************************

	function exists_var(ncid,name)

c checks if name is available

	implicit none

	logical exists_var
	integer ncid
	character*(*) name

	integer var_id

	call nc_get_var_id(ncid,name,var_id)

	exists_var = var_id .gt. 0

	end

c*****************************************************************

	subroutine make_name(pre,it,new)

	implicit none

	character*(*) pre
	integer it
	character*(*) new

	integer i
	character*12 aux
	character*4 ext
	character*7 dir

	ext = '.dat'
	dir = './data/'

	write(aux,'(i12)') it
	do i=1,12
	  if( aux(i:i) .eq. ' ' ) aux(i:i) = '0'
	end do

	new = dir // pre // aux // ext
	new = pre // aux // ext
	
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine print_2d(nxdim,nydim,nx,ny,data)

	implicit none

	integer nxdim,nydim
	integer nx,ny
	real data(nxdim,nydim)

	integer ix,iy

	do iy=1,ny
	  write(6,*) (data(ix,iy),ix=1,nx)
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,nz1,hlv)

	implicit none

	integer ncid
	logical bverb
	character*(*) zcoord
	integer nlvdim
	integer nz
	real zdep(nlvdim)		!depth at cell centers
	integer nz1			!number of values for hlv
	real hlv(nlvdim)		!depth at cell bottom

	logical bldebug,bcenter
	integer z_id
	integer dim_ids(2),dims(2)
	integer ndims
	integer iz
	real htop,h,hindex
	real, parameter :: eps = 1.e-3

	bldebug = .true.
	bldebug = .false.

	if( zcoord .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no zcoord name available'
	  nz = 0
	  nz1 = 0
	  return
	  !stop 'error stop setup_zcoord: no zcoord name'
	end if

	call nc_get_var_id(ncid,zcoord,z_id)
	if( z_id .le. 0 ) then
	  if( bverb ) write(6,*) 'no zcoord variable available'
	  zcoord = ' '
	  return
	end if

	ndims = 0
	call nc_get_var_ndims(ncid,z_id,ndims,dim_ids)
	if( ndims == 0 ) then
	  if( bverb ) write(6,*) 'no zcoord variable available'
	  zcoord = ' '
	  return
	end if
	if( ndims < 1 .or. ndims > 2 ) goto 98
	call nc_get_var_ndims(ncid,z_id,ndims,dim_ids)
	call nc_get_dim_len(ncid,dim_ids(1),nz)
	!call nc_check_var_type(ncid,z_id,'real')

	if( nz .gt. nlvdim ) goto 99

	ndims = 1
	call nc_get_var_data(ncid,zcoord,1,nlvdim,ndims,dims,zdep)
	hlv = zdep
	nz1 = nz

	call check_monotone(nz,zdep,'checking z-coordinates')

	if( nz1 == 1 ) return	!just one layer - hlv not of concern
	if( abs(hlv(1)) < eps ) call depth_shift_up(nz1,hlv)

	hindex = (hlv(2)-hlv(1))/hlv(1)
	if( hindex < 1.5 ) then
	  bcenter = .false.
	else if( hindex < 2.5 ) then
	  bcenter = .true.
	else
	  write(6,*) 'cannot determine if center or bottom'
	  write(6,*) hlv(1:max(4,nz1))
	  stop 'error stop setup_zcoord: strange z coords'
	end if

	if( bcenter ) call depth_center_to_bottom(nz,hlv)

	if( hlv(nz1) < -1. ) nz1 = nz1 - 1	!just in case

	if( bverb ) write(6,*) 'zcoord     : ',z_id,ndims,nz,trim(zcoord)

	if( .not. bldebug ) return

	write(6,*) 'nz1 = ',nz1
	do iz=1,nz1
	  write(6,*) iz,zdep(iz),hlv(iz)
	end do

	return
   98	continue
	write(6,*) 'can manage only 1D depth array: ',ndims
	stop 'error stop setup_zcoord: not 1D'
   99	continue
	write(6,*) 'nz,nlvdim: ',nz,nlvdim
	stop 'error stop setup_zcoord: dimension'
	end

c*****************************************************************

	subroutine depth_shift_up(nz,hlv)

	implicit none

	integer nz
	real hlv(nz)

	integer iz

	do iz=2,nz
	  hlv(iz-1) = hlv(iz)
	end do
	hlv(nz) = 0.
	nz = nz - 1

	end

c*****************************************************************

	subroutine depth_center_to_bottom(nz,hlv)

	implicit none

	integer nz
	real hlv(nz)

	integer iz
	real htop,hbot,h,hd

	htop = 0.
	do iz=1,nz
	  h = hlv(iz)
	  hd = 2. * (h-htop)		!thickness of layer iz
	  hbot = htop + hd
	  hlv(iz) = hbot
	  htop = hbot
	end do

	end

c*****************************************************************

	subroutine write_zcoord(nz,zdep,nz1,hlv)

	implicit none

	integer nz
	real zdep(nz)
	integer nz1
	real hlv(nz)

        integer ix1,ix2,iy1,iy2,iz1,iz2
	integer iz
	character*80 file

	file = 'zcoords.dat'

	call nc_get_domain(ix1,ix2,iy1,iy2,iz1,iz2)

	open(2,file=file,status='unknown',form='formatted')

	write(2,*) 0,0,iz2-iz1+1
	do iz=iz1,iz2
	    write(2,*) iz-iz1+1,zdep(iz),hlv(iz)
	end do

	close(2)

	write(6,*) 'z-coordinates written to file: ',trim(file)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine setup_sealand(ncid,bverb,slmask
     +			,nxdim,nydim,nx,ny,slm)

	implicit none

	integer ncid
	logical bverb
	character*(*) slmask
	integer nxdim,nydim
	integer nx,ny
	real slm(nxdim,nydim)

	logical btime
	integer b_id
	integer dimb_id(10)
	integer ndimb
	integer i,ix,iy
	integer nbx,nby
	real aux(nxdim*nydim)

	slm = 1.			!all is sea

	if( slmask .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no sea_land mask name available'
	  return
	end if

	call nc_get_var_id(ncid,slmask,b_id)
	if( b_id .le. 0 ) then
	  write(6,*) 'no such variable: ',trim(slmask)
	  stop 'error stop setup_sealand: no such variable'
	  !if( bverb ) write(6,*) 'no sea_land mask available'
	  !slmask = ' '
	  !return
	end if

	ndimb = 10
	call nc_get_var_ndims(ncid,b_id,ndimb,dimb_id)
	if( ndimb <= 0 ) goto 97
	call nc_has_time_dimension(ncid,slmask,btime)
	if( btime ) ndimb = ndimb - 1

	!call nc_check_var_type(ncid,b_id,'real')

	call nc_get_dim_len(ncid,dimb_id(1),nbx)
	call nc_get_dim_len(ncid,dimb_id(2),nby)

	if( ndimb .ne. 2 ) goto 97
	if( nbx .ne. nx .or. nby .ne. ny ) goto 98
	if( nx .gt. nxdim .or. ny .gt. nydim ) goto 99

	call nc_get_var_real(ncid,b_id,aux)

	i = 0
	do iy=1,ny
	  do ix=1,nx
	    i = i + 1
	    slm(ix,iy) = aux(i)
	    slm(ix,iy) = 1. - aux(i)		! invert
	  end do
	end do

	if( bverb ) write(6,*) 'sea land mask : ',b_id,nx,ny,trim(slmask)

	return
   97	continue
	write(6,*) 'can manage only 2D sea_land array: ',ndimb
	stop 'error stop setup_sealand: not 2D depth'
   98	continue
	write(6,*) 'nbx,nby,nx,ny: ',nbx,nby,nx,ny
	stop 'error stop setup_sealand: dimensions mismatch'
   99	continue
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	stop 'error stop setup_sealand: dimensions'
	end

c*****************************************************************

	subroutine setup_bathymetry(ncid,bverb,binvertdepth,bathy
     +			,nxdim,nydim,nx,ny,bat)

	implicit none

	integer ncid
	logical bverb,binvertdepth
	character*(*) bathy
	integer nxdim,nydim
	integer nx,ny
	real bat(nxdim,nydim)

	integer b_id
	integer dimb_id(10)
	integer ndimb
	integer i,ix,iy
	integer nbx,nby
	logical btime
	real flag,val
	real aux(nxdim*nydim)
	real, save :: my_flag = -999.

	bat = 0.

	if( bathy .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no bathymetry name available'
	  return
	end if

	call nc_get_var_id(ncid,bathy,b_id)
	if( b_id .le. 0 ) then
	  write(6,*) 'no such variable: ',trim(bathy)
	  stop 'error stop setup_bathymetry: no such variable'
	  !if( bverb ) write(6,*) 'no bathymetry available'
	  !bathy = ' '
	  !return
	end if

	ndimb = 10
	call nc_get_var_ndims(ncid,b_id,ndimb,dimb_id)
	if( ndimb <= 0 ) goto 97
	call nc_has_time_dimension(ncid,bathy,btime)
	if( btime ) ndimb = ndimb - 1

	call get_flag_for_var(ncid,b_id,flag)
	!call nc_check_var_type(ncid,b_id,'real')

	call nc_get_dim_len(ncid,dimb_id(1),nbx)
	call nc_get_dim_len(ncid,dimb_id(2),nby)

	if( ndimb .ne. 2 ) goto 97
	if( nbx .ne. nx .or. nby .ne. ny ) goto 98
	if( nx .gt. nxdim .or. ny .gt. nydim ) goto 99

	call nc_get_var_real(ncid,b_id,aux)

	i = 0
	do iy=1,ny
	  do ix=1,nx
	    i = i + 1
	    val = aux(i)
	    if( val == flag ) then
	      val = my_flag
	    else if( binvertdepth ) then
	      val = -val
	    end if
	    bat(ix,iy) = val
	  end do
	end do

	if( bverb ) write(6,*) 'bathymetry : ',b_id,nx,ny,trim(bathy)

	return
   97	continue
	write(6,*) 'can manage only 2D bathymetry array: ',ndimb
	stop 'error stop setup_bathymetry: not 2D depth'
   98	continue
	write(6,*) 'nbx,nby,nx,ny: ',nbx,nby,nx,ny
	stop 'error stop setup_bathymetry: dimensions mismatch'
   99	continue
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	stop 'error stop setup_bathymetry: dimensions'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine write_2d_fem(file,string,regpar,val)

	implicit none

	character*(*) file
	character*(*) string
	real regpar(9)
	real val(*)

	integer ix,iy
	integer nx,ny,nz
	integer iformat,iunit,nvers,np,nvar,lmax,nlvddi,ntype
	integer datetime(2),ilhkv(1)
	double precision dtime
	real hlv(1),hd(1)

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	nz = 1

	iformat = 1
	dtime = 0.
	nvers = 0
	np = nx * ny
	nvar = 1
	lmax = 1
	nlvddi = 1
	hlv(1) = 10000.
	ilhkv(1) = 1
	hd(1) = 1.
	datetime = 0
	ntype = 10

	iunit = 3
	open(iunit,file=file,status='unknown',form='formatted')

        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,val)

	close(iunit)

	end

c*****************************************************************

	subroutine write_2d_dat(file,nx,ny,val)

	implicit none

	character*(*) file
	integer nx,ny
	real val(nx,ny)

	integer ix,iy

	open(2,file=file,status='unknown',form='formatted')
	write(2,*) 0			! time - not needed
	write(2,*) nx,ny,0
	write(2,'((5e15.7))') ((val(ix,iy),ix=1,nx),iy=1,ny)
	close(2)

	end

c*****************************************************************

	subroutine write_1d_grd(file,n,x,y,val)

	implicit none

	character*(*) file
	integer n
	real x(n)
	real y(n)
	real val(n)

	integer i
	real, parameter :: flag = -999.

	open(1,file=file,status='unknown',form='formatted')

	do i=1,n
	  if( val(i) == flag ) then
	    write(1,1000) 1,i,1,x(i),y(i)
	  else
	    write(1,1000) 1,i,0,x(i),y(i),val(i)
	  end if
	end do

	close(1)

	return
 1000	format(i1,i10,i5,3f14.6)
	end

c*****************************************************************

	subroutine write_2d_grd(file,nx,ny,x,y,val)

	implicit none

	character*(*) file
	integer nx,ny
	real x(nx,ny)
	real y(nx,ny)
	real val(nx,ny)

	integer ix,iy,n
	real, parameter :: flag = -999.

	n = 0

	open(1,file=file,status='unknown',form='formatted')

	do iy=1,ny
	  do ix=1,nx
	    n = n + 1
	    if( val(ix,iy) == flag ) then
	      write(1,1000) 1,n,3,x(ix,iy),y(ix,iy)
	    else
	      write(1,1000) 1,n,2,x(ix,iy),y(ix,iy),val(ix,iy)
	    end if
	  end do
	end do

	close(1)

	return
 1000	format(i1,i10,i5,3f14.6)
	end

c*****************************************************************

	subroutine check_monotone(n,val,text)

	implicit none

	integer n
	real val(n)
	character*(*) text

	logical bgrow
	integer i,imin,imax
	real dv

	if( n <= 1 ) return

	bgrow = val(2) > val(1)

	do i=2,n
	  dv = val(i) - val(i-1)
	  if( dv == 0. ) goto 99
	  if( bgrow .neqv. dv > 0. ) goto 99
	end do

	return
   99	continue
	write(6,*) trim(text)
	write(6,*) 'values not monotone: ',n
	write(6,*) 'problem around i = ',i
	imin = max(1,i-2)
	imax = min(n,i+2)
	write(6,*) val(imin:imax)
	stop 'error stop check_monotone: error in values'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine setup_coordinates(ncid,bverb,namex,namey
     +			,nxdim,nydim,nx,ny,xlon,ylat)

	implicit none

	integer ncid
	logical bverb
	character*(*) namex,namey
	integer nxdim,nydim
	integer nx,ny
	real xlon(nxdim,nydim)
	real ylat(nxdim,nydim)

	logical debug
	integer x_id,y_id
	integer dimx_id(10),dimy_id(10)
	integer ndims,dims(10)
	integer ndimx,ndimy,ndim
	integer lenx,leny
	integer i,ix,iy,n
	integer nxx,nxy,nyx,nyy
	real aux(nxdim*nydim)

	debug = .true.
	debug = .false.

	if( namex .eq. ' ' .or. namey .eq. ' ' ) then
	  write(6,*) 'please provide coordinate names manually'
	  stop 'error stop setup_coordinates: no coordinate names'
	end if

	call nc_get_var_id(ncid,namex,x_id)
	ndimx = 10
	call nc_get_var_ndims(ncid,x_id,ndimx,dimx_id)
	if( ndimx <= 0 ) goto 97

	call nc_get_var_id(ncid,namey,y_id)
	ndimy = 10
	call nc_get_var_ndims(ncid,y_id,ndimy,dimy_id)
	if( ndimx <= 0 ) goto 97

	if( debug ) then
	  write(6,*) 'coordinates: ',x_id,ndimx,namex(1:15)
	  write(6,*) 'coordinates: ',y_id,ndimy,namey(1:15)
	  write(6,*) ndimx,'  ',trim(namex)
	  do i=1,ndimx
	    call nc_get_dim_len(ncid,dimx_id(i),nx)
	    write(6,*) 'x ',i,nx
	  end do
	  write(6,*) ndimy,'  ',trim(namey)
	  do i=1,ndimy
	    call nc_get_dim_len(ncid,dimy_id(i),ny)
	    write(6,*) 'y ',i,ny
	  end do
	end if

	!call nc_check_var_type(ncid,x_id,'real')
	!call nc_check_var_type(ncid,y_id,'real')

	if( ndimx .ne. ndimy ) then
	  write(6,*) 'ndimx,ndimy: ',ndimx,ndimy
	  stop 'error stop setup_coordinates: ndimx <> ndimy'
	end if

	ndim = nxdim * nydim
	ndims = 2

	if( ndimx .eq. 1 ) then
	  call nc_get_dim_len(ncid,dimx_id(1),nx)
	  call nc_get_dim_len(ncid,dimy_id(1),ny)
	  if( nx .gt. nxdim .or. ny .gt. nydim ) goto 99

	  call nc_get_var_real(ncid,x_id,aux)
	  call check_monotone(nx,aux,'checking x-coordinates')

	  do iy=1,ny
	    do ix=1,nx
	      xlon(ix,iy) = aux(ix)
	    end do
	  end do

	  call nc_get_var_real(ncid,y_id,aux)
	  call check_monotone(ny,aux,'checking y-coordinates')

	  do iy=1,ny
	    do ix=1,nx
	      ylat(ix,iy) = aux(iy)
	    end do
	  end do
	else if( ndimx .eq. 2 .or. ndimx .eq. 3 ) then
	  call nc_get_dim_len(ncid,dimx_id(1),nxx)
	  call nc_get_dim_len(ncid,dimx_id(2),nxy)
	  call nc_get_dim_len(ncid,dimy_id(1),nyx)
	  call nc_get_dim_len(ncid,dimy_id(2),nyy)
	  if( nxx .ne. nyx .or. nxy .ne. nyy ) goto 98
	  nx = nxx
	  ny = nyy
	  if( nx .gt. nxdim .or. ny .gt. nydim ) goto 99

	  call nc_get_var_data(ncid,namex,1,ndim,ndims,dims,aux)
	  !call nc_get_var_real(ncid,x_id,aux)
	  if( debug ) write(6,*) 'creating x-coordinates ... ',nx,ny
	  i = 0
	  do iy=1,ny
	    do ix=1,nx
	      i = i + 1
	      xlon(ix,iy) = aux(i)
	    end do
	  end do

	  call nc_get_var_data(ncid,namey,1,ndim,ndims,dims,aux)
	  !call nc_get_var_real(ncid,y_id,aux)
	  if( debug ) write(6,*) 'creating y-coordinates ... ',nx,ny
	  i = 0
	  do iy=1,ny
	    do ix=1,nx
	      i = i + 1
	      ylat(ix,iy) = aux(i)
	    end do
	  end do
	else
	  write(6,*) 'coordinates x: ',x_id,ndimx,trim(namex)
	  write(6,*) (dimx_id(i),i=1,ndimx)
	  write(6,*) 'coordinates y: ',y_id,ndimy,trim(namey)
	  write(6,*) (dimy_id(i),i=1,ndimy)
	  call nc_dims_info(ncid)
	  write(6,*) 'ndimx: ',ndimx
	  stop 'error stop setup_coordinates: coordinates dimension'
	end if

	if( bverb ) then
	  write(6,*) 'coordinates: ',x_id,ndimx,nx,trim(namex)
	  write(6,*) 'coordinates: ',y_id,ndimy,ny,trim(namey)
	  !write(6,*) 'coordinates: ',nx,ny
	end if

	return
   97	continue
	write(6,*) 'ndimx,ndimy: ',ndimx,ndimy
	stop 'error stop setup_coordinates: dimensions error'
   98	continue
	write(6,*) 'nxx,nyx,nxy,nyy: ',nxx,nyx,nxy,nyy
	stop 'error stop setup_coordinates: dimensions mismatch'
   99	continue
	write(6,*) 'namex: ',trim(namex)
	write(6,*) 'namey: ',trim(namey)
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	write(6,*) 'probable error: no valid dimension x/y found'
	write(6,*) '(insert dimensions in ncnames_add_dimensions)'
	stop 'error stop setup_coordinates: dimensions'
	end

c*****************************************************************

	subroutine write_coordinates(nx,ny,xlon,ylat)

	implicit none

	integer nx,ny
	real xlon(nx,ny)
	real ylat(nx,ny)

	integer ix,iy,n

	n = 0

	open(1,file='coords_orig.grd',status='unknown',form='formatted')
	open(2,file='coords_orig.dat',status='unknown',form='formatted')

	write(2,*) nx,ny,0
	do iy=1,ny
	  do ix=1,nx
	    n = n + 1
	    write(1,1000) 1,n,0,xlon(ix,iy),ylat(ix,iy)
	    write(2,*) ix,iy,xlon(ix,iy),ylat(ix,iy)
	  end do
	end do

	close(1)
	close(2)

	write(6,*) 'coordinates written to files: '//
     +			'coords_orig.grd coords_orig.dat'

	return
 1000	format(i1,i10,i5,2f14.6)
	end

c*****************************************************************

	subroutine write_regular_coordinates_as_grd(regpar)

	implicit none

	real regpar(9)

	integer ix,iy,n
	integer nx,ny
	real dx,dy,x0,y0,x1,y1,flag
	real x,y

	call get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	open(1,file='coords_new.grd',status='unknown',form='formatted')

	n = 0
	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  do ix=1,nx
	    n = n + 1
	    x = x0 + (ix-1)*dx
	    write(1,1000) 1,n,4,x,y
	  end do
	end do

	close(1)

	write(6,*) 'coordinates written to files: coords_new.grd'

	return
 1000	format(i1,i10,i5,2f14.6)
	end

c*****************************************************************

	subroutine write_2d_grd_regular(file,regpar,scal)

	implicit none

	character*(*) file
	real regpar(9)
	real scal(*)

	integer ix,iy,n
	integer nx,ny
	real dx,dy,x0,y0,x1,y1,flag
	real x,y,val
	double precision ddx,ddy

	call get_regpar(regpar,nx,ny,dx,dy,x0,y0,x1,y1,flag)

	open(1,file=file,status='unknown',form='formatted')

	ddx = dx
	ddy = dy

	n = 0
	do iy=1,ny
	  y = y0 + (iy-1)*ddy
	  do ix=1,nx
	    n = n + 1
	    val = scal(n)
	    x = x0 + (ix-1)*ddx
	    if( val == flag ) then
	      write(1,1000) 1,n,5,x,y
	    else
	      write(1,1000) 1,n,6,x,y,val
	    end if
	  end do
	end do

	close(1)

	write(6,*) 'scalar written to file: ',trim(file)

	return
 1000	format(i1,i10,i5,3f14.6)
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_flag_for_var(ncid,var_id,value)

	integer ncid,var_id
	real value

	character*80 atext
	double precision avalue

	call nc_get_var_attrib(ncid,var_id,'_FillValue',atext,avalue)
	if( avalue == 0. ) then
	  call nc_get_var_attrib(ncid,var_id,'FillValue_',atext,avalue)
	end if

	value = avalue

	end

c*****************************************************************

