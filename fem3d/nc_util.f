
c*****************************************************************
c*****************************************************************
c*****************************************************************

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

	subroutine compress_data(nx,ny,nz,rin,rout)

	implicit none

	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	integer nx,ny,nz
	real rin(nx,ny,nz)
	real rout(1)

	integer ip,ix,iy,iz

	ip = 0

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

	subroutine write_list(ncid,name,ifreq)

	implicit none

	integer ncid
	character*(*) name
	integer ifreq

	integer nit,n,it
	character*30 new

        call nc_get_time_recs(ncid,nit)

        open(3,file=name,form='formatted',status='unknown')

        do n=1,nit
	  if( ifreq == 0 .or. mod(n-1,ifreq) .eq. 0 ) then
            call handle_nc_time(ncid,n,it)
            call make_name('_',it,new)
            write(3,'(i12,2x,a)') it,new
	  end if
	!write(66,*) n,nit,ifreq,it,mod(n-1,ifreq)
        end do

        close(3)

	end

c*****************************************************************

	subroutine write_2d_data(ncid,name,pre,ndim,aux,fact)

	implicit none

	integer ncid
	character*(*) name,pre
	integer ndim
	real aux(*)
	real fact

	integer nit,n,ierror,it

        call nc_get_time_recs(ncid,nit)

        ierror = 0
        call check_var(ncid,name,ierror)
        if( ierror .gt. 0 ) return

        do n=1,nit
          call handle_nc_time(ncid,n,it)
          call handle_data_2d(ncid,name,pre,it,n,ndim,aux,fact)
        end do

	end

c*****************************************************************

	subroutine write_3d_data(ncid,name,pre,ndim,aux,fact)

	implicit none

	integer ncid
	character*(*) name,pre
	integer ndim
	real aux(*)
	real fact

	integer nit,n,ierror,it

        call nc_get_time_recs(ncid,nit)

        ierror = 0
        call check_var(ncid,name,ierror)
        if( ierror .gt. 0 ) return

        do n=1,nit
          call handle_nc_time(ncid,n,it)
          call handle_data_3d(ncid,name,pre,it,n,ndim,aux,fact)
        end do

	end

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

	subroutine setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,zzdep)

	implicit none

	integer ncid
	logical bverb
	character*(*) zcoord
	integer nlvdim
	integer nz
	real zdep(nlvdim)
	real zzdep(nlvdim)

	integer z_id
	integer dim_ids(2),dims(2)
	integer ndims
	integer iz
	real htop,h

	if( zcoord .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no zcoord name available'
	  nz = 0
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
	if( ndims < 1 .or. ndims > 2 ) goto 98
	call nc_get_var_ndims(ncid,z_id,ndims,dim_ids)
	call nc_get_dim_len(ncid,dim_ids(1),nz)
	!call nc_check_var_type(ncid,z_id,'real')

	if( nz .gt. nlvdim ) goto 99

	ndims = 1
	call nc_get_var_data(ncid,zcoord,1,nlvdim,ndims,dims,zdep)
	!call nc_get_var_real(ncid,z_id,zdep)

	htop = 0.
	do iz=1,nz
	  h = zdep(iz)
	  htop = htop + 2.*(h-htop)
	  zzdep(iz) = htop
	end do

	if( bverb ) write(6,*) 'zcoord     : ',z_id,ndims,nz,trim(zcoord)

	return
   98	continue
	write(6,*) 'can manage only 1D depth array: ',ndims
	stop 'error stop setup_zcoord: not 1D'
   99	continue
	write(6,*) 'nz,nlvdim: ',nz,nlvdim
	stop 'error stop setup_zcoord: dimension'
	end

c*****************************************************************

	subroutine write_zcoord(nz,zdep,zzdep)

	implicit none

	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	integer nz
	real zdep(nz)
	real zzdep(nz)

	integer iz

	open(2,file='zcoord.dat',status='unknown',form='formatted')

	write(2,*) 0,0,iz2-iz1+1
	do iz=iz1,iz2
	    write(2,*) iz-iz1+1,zdep(iz),zzdep(iz)
	end do

	close(2)

	write(6,*) 'z-coordinates written to file: zcoords.grd'

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine setup_sealand(ncid,bverb,slmask
     +			,nxdim,nydim,nx,ny,slm,aux)

	implicit none

	integer ncid
	logical bverb
	character*(*) slmask
	integer nxdim,nydim
	integer nx,ny
	real slm(nxdim,nydim)
	real aux(nxdim*nydim)

	integer b_id
	integer dimb_id(10)
	integer ndimb
	integer i,ix,iy
	integer nbx,nby
	logical btime

	slm = 1.			!all is sea

	if( slmask .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no sea_land mask name available'
	  return
	end if

	call nc_get_var_id(ncid,slmask,b_id)
	if( b_id .le. 0 ) then
	  if( bverb ) write(6,*) 'no sea_land mask available'
	  slmask = ' '
	  return
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
     +			,nxdim,nydim,nx,ny,bat,aux)

	implicit none

	integer ncid
	logical bverb,binvertdepth
	character*(*) bathy
	integer nxdim,nydim
	integer nx,ny
	real bat(nxdim,nydim)
	real aux(nxdim*nydim)

	integer b_id
	integer dimb_id(10)
	integer ndimb
	integer i,ix,iy
	integer nbx,nby
	logical btime
	real flag,val
	real, save :: my_flag = -999.

	bat = 0.

	if( bathy .eq. ' ' ) then
	  if( bverb ) write(6,*) 'no bathymetry name available'
	  return
	end if

	call nc_get_var_id(ncid,bathy,b_id)
	if( b_id .le. 0 ) then
	  if( bverb ) write(6,*) 'no bathymetry available'
	  bathy = ' '
	  return
	end if

	ndimb = 10
	call nc_get_var_ndims(ncid,b_id,ndimb,dimb_id)
	if( ndimb <= 0 ) goto 97
	call nc_has_time_dimension(ncid,bathy,btime)
	if( btime ) ndimb = ndimb - 1

	call get_flag_for var(ncid,b_id,flag)
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

	subroutine write_2d_reg(file,nxdim,nydim,nx,ny,val)

	implicit none

	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	character*(*) file
	integer nxdim,nydim
	integer nx,ny
	real val(nxdim,nydim)

	integer ix,iy

	open(2,file=file,status='unknown',form='formatted')
	write(2,*) 0			! time - not needed
	write(2,*) ix2-ix1+1,iy2-iy1+1,0
	write(2,'((5e15.7))') ((val(ix,iy),ix=ix1,ix2),iy=iy1,iy2)
	close(2)

	end

c*****************************************************************

	subroutine write_2d_grd(file,nxdim,nydim,nx,ny,x,y,val)

	implicit none

	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	character*(*) file
	integer nxdim,nydim
	integer nx,ny
	real x(nxdim,nydim)
	real y(nxdim,nydim)
	real val(nxdim,nydim)

	integer ix,iy,n
	real, save :: flag = -999.

	open(1,file=file,status='unknown',form='formatted')

	do iy=iy1,iy2
	  do ix=ix1,ix2
	    n = n + 1
	    if( val(ix,iy) /= flag ) then
	      write(1,1000) 1,n,0,x(ix,iy),y(ix,iy),val(ix,iy)
	    end if
	  end do
	end do

	close(1)

	return
 1000	format(i1,i10,i5,3f14.6)
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine setup_coordinates(ncid,bverb,namex,namey
     +			,nxdim,nydim,nx,ny,xlon,ylat,aux)

	implicit none

	integer ncid
	logical bverb
	character*(*) namex,namey
	integer nxdim,nydim
	integer nx,ny
	real xlon(nxdim,nydim)
	real ylat(nxdim,nydim)
	real aux(nxdim*nydim)

	integer x_id,y_id
	integer dimx_id(10),dimy_id(10)
	integer ndims,dims(3)
	integer ndimx,ndimy,ndim
	integer lenx,leny
	integer i,ix,iy,n
	integer nxx,nxy,nyx,nyy

	if( namex .eq. ' ' .or. namey .eq. ' ' ) then
	  write(6,*) 'please provide coordinate names manually'
	  stop 'error stop setup_coordinates: no coordinate names'
	end if

	call nc_get_var_id(ncid,namex,x_id)
	ndimx = 10
	call nc_get_var_ndims(ncid,x_id,ndimx,dimx_id)
	if( ndimx <= 0 ) goto 97
	!write(6,*) 'coordinates: ',x_id,ndimx,namex(1:15)

	call nc_get_var_id(ncid,namey,y_id)
	ndimy = 10
	call nc_get_var_ndims(ncid,y_id,ndimy,dimy_id)
	if( ndimx <= 0 ) goto 97
	!write(6,*) 'coordinates: ',y_id,ndimy,namey(1:15)

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
	  do iy=1,ny
	    do ix=1,nx
	      xlon(ix,iy) = aux(ix)
	    end do
	  end do

	  call nc_get_var_real(ncid,y_id,aux)
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
	  i = 0
	  do iy=1,ny
	    do ix=1,nx
	      i = i + 1
	      xlon(ix,iy) = aux(i)
	    end do
	  end do

	  call nc_get_var_data(ncid,namey,1,ndim,ndims,dims,aux)
	  !call nc_get_var_real(ncid,y_id,aux)
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
	write(6,*) 'nx,nxdim: ',nx,nxdim
	write(6,*) 'ny,nydim: ',ny,nydim
	stop 'error stop setup_coordinates: dimensions'
	end

c*****************************************************************

	subroutine write_coordinates(nxdim,nydim,nx,ny,xlon,ylat)

	implicit none

	!include 'bounds.h'
        integer ibound,nnx,nny,nnz
        integer ix1,ix2,iy1,iy2,iz1,iz2
        common /bound_common/ ibound,nnx,nny,nnz,ix1,ix2,iy1,iy2,iz1,iz2
        save /bound_common/

	integer nxdim,nydim
	integer nx,ny
	real xlon(nxdim,nydim)
	real ylat(nxdim,nydim)

	integer ix,iy,n
	integer ixx,iyy

	n = 0

	open(1,file='coords.grd',status='unknown',form='formatted')
	open(2,file='coords.dat',status='unknown',form='formatted')

	!write(6,*) '***************** ',ix1,ix2,iy1,iy2,iz1,iz2

	write(2,*) ix2-ix1+1,iy2-iy1+1,0
	do iy=iy1,iy2
	  do ix=ix1,ix2
	    n = n + 1
	    ixx = ix-ix1+1
	    iyy = iy-iy1+1
	    write(1,1000) 1,n,0,xlon(ix,iy),ylat(ix,iy)
	    write(2,*) ixx,iyy,xlon(ix,iy),ylat(ix,iy)
	  end do
	end do

	close(1)
	close(2)

	write(6,*) 'coordinates written to files: coords.grd coords.dat'

	return
 1000	format(i1,i10,i5,2f14.6)
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine find_nc_file_type(ncid,iftype)

c tries to find the file type to be read

	implicit none

	integer ncid
	integer iftype

	character*80 atext

	iftype = 0

	call nc_get_global_attr(ncid,'TITLE',atext)
	if( atext .eq. ' OUTPUT FROM WRF V3.4 MODEL' ) then
	  iftype = 1
	  write(6,*) 'file type: ',iftype,atext(1:40)
	  return
	end if

	call nc_get_global_attr(ncid,'source',atext)
	!write(6,*) 'source (1): ',atext(1:30)
	!write(6,*) 'source (2): ','MFS SYS4a4'
	if( atext(1:10) .eq. 'MFS SYS4a4' ) then
	  iftype = 2
	  write(6,*) 'file type: ',iftype,atext(1:40)
	  return
	end if

	call nc_get_global_attr(ncid,'type',atext)
	if( atext .eq. 'ROMS/TOMS history file' ) then
	  iftype = 3
	  write(6,*) 'file type: ',iftype,atext(1:40)
	  return
	end if

	call nc_get_global_attr(ncid,'CDO',atext)
	!write(6,*) 'source (1): ',atext(1:36)
	!write(6,*) 'source (2): ','Climate Data Operators version 1.5.5'
	if( atext(1:36) .eq. 'Climate Data Operators version 1.5.5' ) then
	  iftype = 2
	  write(6,*) 'file type: ',iftype,atext(1:40)
	  return
	end if

	call nc_get_global_attr(ncid,'institution',atext)
	if( atext(1:32) .eq. 'European Centre for Medium-Range' ) then
	  iftype = 4
	  write(6,*) 'file type: ',iftype,atext(1:40)
	  return
	end if

	if( iftype .eq. 0 ) write(6,*) 'Cannot determine file type'

	end

c*****************************************************************

	subroutine set_names(iftype,time,namex,namey,zcoord,bathy,slmask)

	implicit none

	integer iftype
	character*(*) time,namex,namey,zcoord,bathy,slmask

	character*20 time_d

	time = 'time'
	namex = 'lon'
	namey = 'lat'
	zcoord = ' '
	bathy = ' '
	slmask = ' '

	if( iftype .eq. 1 ) then
	  namex = 'XLONG'
	  namey = 'XLAT'
	  slmask = 'LANDMASK'
	else if( iftype .eq. 2 ) then
	  zcoord = 'depth'
	else if( iftype .eq. 3 ) then
	  time = 'ocean_time'
	  namex = 'lon_rho'
	  namey = 'lat_rho'
	  zcoord = 'Cs_r'
	  bathy = 'h'
	else if( iftype .eq. 4 ) then
	  !nothing
	end if

	time_d = time
	call nc_set_time_name(time_d,time)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_nc_dimensions(ncid,bverb,nt,nx,ny,nz)

	implicit none

	integer ncid
	logical bverb
	integer nt,nx,ny,nz

	integer dim_id,n,ndims,time_id,nn,i
	integer ixdim,iydim,izdim,itdim
	character*80 name,time_v,time_d
	logical :: bdebug = .true.

	character(len=11), save :: xdims(4) =	(/
     +		 'x          '
     +		,'xpos       '
     +		,'lon        '
     +		,'west_east  '
     +						/)
	character(len=11), save :: ydims(4) =	(/
     +		 'y          '
     +		,'ypos       '
     +		,'lat        '
     +		,'south_north'
     +						/)
	character(len=15), save :: zdims(3) =	(/
     +		 'z              '
     +		,'zpos           '
     +		,'bottom_top_stag'
     +						/)
	character(len=4), save :: tdims(2) =	(/
     +		 'time'
     +		,'Time'
     +						/)

	call nc_get_dim_totnum(ncid,ndims)

	ixdim = 0
	iydim = 0
	izdim = 0
	itdim = 0
	nx = 0
	ny = 0
	nz = 0
	nt = 0
	time_id = 0
	time_v = ' '

	do dim_id=1,ndims
	  call nc_get_dim_name(ncid,dim_id,name)
	  call nc_get_dim_len(ncid,dim_id,n)

	  do i=1,size(xdims)
	    if( name == xdims(i) ) ixdim = dim_id
	  end do

	  do i=1,size(ydims)
	    if( name == ydims(i) ) iydim = dim_id
	  end do

	  do i=1,size(zdims)
	    if( name == zdims(i) ) izdim = dim_id
	  end do

	  do i=1,size(tdims)
	    if( name == tdims(i) ) itdim = dim_id
	  end do
	end do

	if( bverb ) write(6,*) 'dimensions: '

	if( ixdim > 0 ) then
	  call nc_get_dim_name(ncid,ixdim,name)
	  call nc_get_dim_len(ncid,ixdim,nx)
	  if( bverb ) write(6,*) '   xdim: ',nx,'  (',trim(name),')'
	end if

	if( iydim > 0 ) then
	  call nc_get_dim_name(ncid,iydim,name)
	  call nc_get_dim_len(ncid,iydim,ny)
	  if( bverb ) write(6,*) '   ydim: ',ny,'  (',trim(name),')'
	end if

	if( izdim > 0 ) then
	  call nc_get_dim_name(ncid,izdim,name)
	  call nc_get_dim_len(ncid,izdim,nz)
	  if( bverb ) write(6,*) '   zdim: ',nz,'  (',trim(name),')'
	end if

	if( itdim > 0 ) then	!handle time dimension
	  call nc_get_dim_name(ncid,itdim,name)
	  call nc_get_dim_len(ncid,itdim,nt)
	  if( bverb ) write(6,*) '   tdim: ',nt,'  (',trim(name),')'
	  time_d = name
	  call nc_set_time_name(time_d,time_v)
	end if

	end

c*****************************************************************

	subroutine get_flag_for var(ncid,var_id,value)

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

