
!*************************************************************************

	subroutine get_dims(nitem,idx,idy,idz)

! get dimensions

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer idx,idy,idz

	integer id
        type(dim_item) :: ditem

	write(6,*) 'determining dimensions'

	ditem = nitem%ditem
	idx = 0
	idy = 0
	idz = 0

	do id=1,ditem%ndims
	  !write(6,*) ditem%len(id),'   ',trim(ditem%name(id))
	  if( ditem%name(id) == 'x' ) idx = ditem%len(id)
	  if( ditem%name(id) == 'y' ) idy = ditem%len(id)
	  if( ditem%name(id) == 'z' ) idz = ditem%len(id)
	end do

	end

!*************************************************************************

	subroutine read_mask(nitem,nx,ny,nz,xx,yy,mask,smask,daver)

! handle mask

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer nx,ny,nz
	double precision xx(nx,ny)
	double precision yy(nx,ny)
	double precision mask(nx,ny,nz)
	double precision smask(nx,ny)
	double precision daver(nz)

	integer id,ncid,ix,iy,iz
	integer idx,idy,idz,ivx,ivy,ivz,ivm
	integer nvars,varid
        type(var_item) :: vitem
        type(att_item) :: aitem
        type(dim_item) :: ditem

	integer itype
	real, parameter :: flag = -999.
	double precision, allocatable :: xlon(:),ylat(:)
	character*80 file

	write(6,*) 'determining mask information'

	ncid = nitem%ncid
	ditem = nitem%ditem
	idx = 0
	idy = 0
	idz = 0
	ivx = 0
	ivy = 0
	ivz = 0
	ivm = 0

	do id=1,ditem%ndims
	  !write(6,*) ditem%len(id),'   ',trim(ditem%name(id))
	  if( ditem%name(id) == 'x' ) idx = ditem%len(id)
	  if( ditem%name(id) == 'y' ) idy = ditem%len(id)
	  if( ditem%name(id) == 'z' ) idz = ditem%len(id)
	end do

	!write(6,*) 'horizontal dimensions: ',idx,idy
	!write(6,*) 'vertical dimensions: ',idz

	if( idx /= nx ) stop 'dimension error'
	if( idy /= ny ) stop 'dimension error'
	if( idz /= nz ) stop 'dimension error'

	allocate(xlon(idx))
	allocate(ylat(idy))

        nvars = nitem%nvars
        !write(6,*) 'variables: ',nvars
        !call ncf_print_variable_header
        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          !call ncf_print_variable(vitem)
	  if( vitem%name == 'longitude' ) ivx = varid
	  if( vitem%name == 'latitude' ) ivy = varid
	  if( vitem%name == 'depth' ) ivz = varid
	  if( vitem%name == 'mask' ) ivm = varid
        end do

	write(6,*) 'lat/lon found: ',ivx,ivy
	write(6,*) 'depth found: ',ivz
	write(6,*) 'mask found: ',ivm

        call ncf_var_inf(ncid,ivx,vitem)
	call ncf_var_get(ncid,vitem)
	xlon = vitem%value
        call ncf_var_inf(ncid,ivy,vitem)
	call ncf_var_get(ncid,vitem)
	ylat = vitem%value

	do iy=1,idy
	  xx(:,iy) = xlon(:)
	end do
	do ix=1,idx
	  yy(ix,:) = ylat(:)
	end do

        call ncf_var_inf(ncid,ivz,vitem)
	call ncf_var_get(ncid,vitem)
	daver = vitem%value

        call ncf_var_inf(ncid,ivm,vitem)
	call ncf_var_get(ncid,vitem)
	mask = reshape(vitem%value,(/idx,idy,idz/))

	smask = 1
	do iy=1,idy
	  do ix=1,idx
	    if( isnan(mask(ix,iy,1)) ) smask(ix,iy) = flag
	    do iz=1,nz
	      if( isnan(mask(ix,iy,iz)) ) mask(ix,iy,iz) = flag
	    end do
	  end do
	end do

	file='mask.grd'
	itype = 0
	call write_grid(file,itype,nx,ny,xx,yy,smask)

	end

!*************************************************************************

	subroutine read_tmask(nitem,nx,ny,nz,mask,smask)

! handle mask

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer nx,ny,nz
	double precision mask(nx,ny,nz)
	double precision smask(nx,ny)

	integer id,ncid,ix,iy,iz
	integer idx,idy,idz,ivx,ivy,ivz,ivm
	integer nvars,varid
        type(var_item) :: vitem
        type(att_item) :: aitem
        type(dim_item) :: ditem

	integer itype
	real, parameter :: flag = -999.
	character*80 file

	ncid = nitem%ncid
	ditem = nitem%ditem
	idx = 0
	idy = 0
	idz = 0
	ivx = 0
	ivy = 0
	ivz = 0
	ivm = 0

	do id=1,ditem%ndims
	  !write(6,*) ditem%len(id),'   ',trim(ditem%name(id))
	  if( ditem%name(id) == 'x' ) idx = ditem%len(id)
	  if( ditem%name(id) == 'y' ) idy = ditem%len(id)
	  if( ditem%name(id) == 'z' ) idz = ditem%len(id)
	end do

	!write(6,*) 'horizontal dimensions: ',idx,idy
	!write(6,*) 'vertical dimensions: ',idz

	if( idx /= nx ) stop 'dimension error'
	if( idy /= ny ) stop 'dimension error'
	if( idz /= nz ) stop 'dimension error'

        nvars = nitem%nvars
        !write(6,*) 'variables: ',nvars
        !call ncf_print_variable_header
        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          !call ncf_print_variable(vitem)
	  if( vitem%name == 'longitude' ) ivx = varid
	  if( vitem%name == 'latitude' ) ivy = varid
	  if( vitem%name == 'depth' ) ivz = varid
	  if( vitem%name == 'tmask' ) ivm = varid
        end do

	write(6,*) 'lat/lon found: ',ivx,ivy
	write(6,*) 'depth found: ',ivz
	write(6,*) 'mask found: ',ivm

        call ncf_var_inf(ncid,ivm,vitem)
	call ncf_var_get(ncid,vitem)
	mask = reshape(vitem%value,(/idx,idy,idz/))

	smask = 1
	do iy=1,idy
	  do ix=1,idx
	    if( mask(ix,iy,1) == 0. ) smask(ix,iy) = flag
	    do iz=1,nz
	      if( mask(ix,iy,iz) == 0. ) mask(ix,iy,iz) = flag
	    end do
	    !if( isnan(mask(ix,iy,1)) ) smask(ix,iy) = flag
	    !do iz=1,nz
	    !  if( isnan(mask(ix,iy,iz)) ) mask(ix,iy,iz) = flag
	    !end do
	    !write(6,*) mask(ix,iy,:)
	  end do
	end do

	!file='tmask.grd'
	!itype = 0
	!call write_grid(file,itype,nx,ny,xx,yy,smask)

	end

!*************************************************************************

	subroutine read_bathy(nitem,nx,ny,xx,yy,smask,bathy)

! handle bathymetry

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer nx,ny
	double precision xx(nx,ny)
	double precision yy(nx,ny)
	double precision smask(nx,ny)
	double precision bathy(nx,ny)

	integer id,ncid,ix,iy
	integer idx,idy,ivd
	integer nvars,varid
	integer itype
	real, parameter :: flag = -999.
        type(var_item) :: vitem
        type(att_item) :: aitem
        type(dim_item) :: ditem

	integer iu,node
	real x,y,d
	double precision, allocatable :: xlon(:),ylat(:)
	character*80 file

	write(6,*) 'determining bathymetry information'

	ncid = nitem%ncid
	ditem = nitem%ditem
	idx = 0
	idy = 0
	ivd = 0

	do id=1,ditem%ndims
	  !write(6,*) ditem%len(id),'   ',trim(ditem%name(id))
	  if( ditem%name(id) == 'x' ) idx = ditem%len(id)
	  if( ditem%name(id) == 'y' ) idy = ditem%len(id)
	end do

	!write(6,*) 'horizontal dimensions: ',idx,idy
	!mask(ix,iy,1)write(6,*) 'vertical dimensions: ',idz

	if( idx /= nx ) stop 'dimension error'
	if( idy /= ny ) stop 'dimension error'

        nvars = nitem%nvars
        !write(6,*) 'variables: ',nvars
        !call ncf_print_variable_header
        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          !call ncf_print_variable(vitem)
	  if( vitem%name == 'Bathymetry' ) ivd = varid
        end do

	write(6,*) 'bathymetry found: ',ivd

        call ncf_var_inf(ncid,ivd,vitem)
	call ncf_var_get(ncid,vitem)
	bathy = reshape(vitem%value,(/idx,idy/))

	where( smask == flag ) bathy = flag

	file='bathy.grd'
	itype = 0
	call write_grid(file,itype,nx,ny,xx,yy,bathy)

	end

!*************************************************************************

	subroutine make_depth(nz,daver,dbottom,dlayer)

! makes depth values

	use ncf

	implicit none

	integer nz
	double precision daver(nz)
	double precision dbottom(nz)
	double precision dlayer(nz)

	integer iz,iu
	double precision dtop,dh2
	character*80 file,header,format

	write(6,*) 'determining layer information'

	dtop = 0.
	do iz=1,nz
	  dh2 = daver(iz) - dtop
	  dlayer(iz) = 2*dh2
	  dbottom(iz) = dtop + dlayer(iz)
	  dtop = dbottom(iz)
	end do
  
	format = '(i8,3f16.3)'
	header = '#  layer    center_depth    bottom_depth' //
     +			' layer_thickness'
	iu = 1
	file = 'depth.txt'
	open(iu,file=file,status='unknown',form='formatted')

	write(6,*) 'maximum number of layers: ',nz
	write(6,'(a)') header
	write(iu,'(a)') '# layers'
	write(iu,'(i8)') nz
	write(iu,'(a)') header
	do iz=1,nz
	  write(6,format) iz,daver(iz),dbottom(iz),dlayer(iz)
	  write(iu,format) iz,daver(iz),dbottom(iz),dlayer(iz)
	end do

	write(6,*) 'file written: ',trim(file)
	close(iu)

	end

!*************************************************************************

	subroutine make_layers(nx,ny,nz,xx,yy,mask,layers)

! makes layers

	use ncf

	implicit none

	integer nx,ny,nz
	double precision xx(nx,ny)
	double precision yy(nx,ny)
	double precision mask(nx,ny,nz)
	integer layers(nx,ny)

	integer ix,iy,iz
	integer itype
	real, parameter :: flag = -999.
	double precision, allocatable :: aux(:,:)
	character*80 file

	do iy=1,ny
	  do ix=1,nx
	    do iz=1,nz
	      if( mask(ix,iy,iz) == flag ) exit
	    end do
	    layers(ix,iy) = iz - 1
	  end do
	end do

	allocate(aux(nx,ny))
	aux = layers
	where( layers == 0 ) aux = flag

	file='layers.grd'
	itype = 0
	call write_grid(file,itype,nx,ny,xx,yy,aux)

	end

!*************************************************************************

	subroutine make_boxes(file,nx,ny,xx,yy,boxes)

! makes box information

	use ncf

	implicit none

	character*80 file
	integer nx,ny
	double precision xx(nx,ny)
	double precision yy(nx,ny)
	integer boxes(nx,ny)

	logical bdebug
	integer ix,iy
	integer iu,nbox
	integer i,j,n,ibox
	integer itype,ifound,ntot
	integer ixmin,ixmax,iymin,iymax
	real, parameter :: flag = -999.
	real xmin,ymin,xmax,ymax,x,y
	integer, allocatable :: npoints(:)
	double precision, allocatable :: aux(:,:)
	real, allocatable :: xp(:),yp(:)
	character*80 header

	logical inpoly

	bdebug = .true.
	bdebug = .false.
	boxes = 0

	write(6,*) 'determining box information'

	iu = 1
	open(iu,file=file,status='unknown',form='formatted')
	read(iu,*) nbox
	write(6,*) 'nbox = ',nbox
	header = '               box line_points       found       total'
	write(6,'(a)') trim(header)

	do ibox=1,nbox
	  read(iu,*) i,n
	  !write(6,*) 'box,points:',i,n
	  if( i /= ibox ) stop 'lines read mismatch (1)'
	  allocate(xp(n),yp(n))
	  do j=1,n
	    read(iu,*) i,xp(j),yp(j)
	    if( i /= j ) stop 'lines read mismatch (2)'
	  end do

	  xmax = maxval(xp)
	  ymax = maxval(yp)
	  xmin = minval(xp)
	  ymin = minval(yp)

	  if( bdebug ) write(6,*) 'min/max: ',xmin,ymin,xmax,ymax

	  iymin = 0
	  iymax = 0
	  do iy=1,ny
	    y = yy(1,iy)
	    if( iymin == 0 .and. y > ymin ) iymin = iy
	    if( iymax == 0 .and. y > ymax ) iymax = iy
	  end do
	  iymin = max(1,iymin-1)
	  if( yy(1,iymin) > ymin ) stop 'internal error (5)'
	  if( yy(1,iymax) < ymax ) stop 'internal error (5)'
	  if( bdebug ) write(6,*) iymin,iymax,yy(1,iymin),yy(1,iymax)
	
	  ixmin = 0
	  ixmax = 0
	  do ix=1,nx
	    x = xx(ix,1)
	    if( ixmin == 0 .and. x > xmin ) ixmin = ix
	    if( ixmax == 0 .and. x > xmax ) ixmax = ix
	  end do
	  ixmin = max(1,ixmin-1)
	  if( xx(ixmin,1) > xmin ) stop 'internal error (5)'
	  if( xx(ixmax,1) < xmax ) stop 'internal error (5)'
	  if( bdebug ) write(6,*) ixmin,ixmax,xx(ixmin,1),xx(ixmax,1)

	  ifound = 0
	  do iy=iymin,iymax
	    do ix=ixmin,ixmax
	      x = xx(ix,iy)
	      y = yy(ix,iy)
	      if( inpoly(n,xp,yp,x,y) ) then
		ifound = ifound + 1
		boxes(ix,iy) = ibox
	      end if
	    end do
	  end do

	  ntot = (iymax-iymin+1)*(ixmax-ixmin+1)
	  write(6,*) '     ',ibox,n,ifound,ntot

	  deallocate(xp,yp)
	end do

	close(iu)

	allocate(npoints(0:nbox))
	npoints = 0

	do iy=1,ny
	  do ix=1,nx
	    ibox = boxes(ix,iy)
	    npoints(ibox) = npoints(ibox) + 1
	  end do
	end do

	allocate(aux(nx,ny))
	aux = boxes
	where( boxes == 0 ) aux = flag

	file='boxes.grd'
	itype = -1
	call write_grid(file,itype,nx,ny,xx,yy,aux)

	end

!*************************************************************************

	subroutine make_box_layers(nbox,nx,ny,boxes,layers,lmaxs)

! makes box layer information

	implicit none

	integer nbox
	integer nx,ny
	integer boxes(nx,ny)
	integer layers(nx,ny)
	integer lmaxs(nbox)

	integer ix,iy,ibox,lmax

	lmaxs = 0

	do ix=1,nx
	  do iy=1,ny
	    ibox = boxes(ix,iy)
	    lmax = layers(ix,iy)
	    if( ibox == 0 ) cycle
	    if( lmax > lmaxs(ibox) ) lmaxs(ibox) = lmax
	  end do
	end do

	write(6,*) 'determining box layers'
	write(6,'(a)') '         box      layers'
	do ibox=1,nbox
	  write(6,*) ibox,lmaxs(ibox)
	end do

	end

!*************************************************************************

	subroutine handle_intp(nitem,nbox,nx,ny,nz
     +				,boxes,layers,mask,dlayer)

! handle interpolation

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer nbox
	integer nx,ny,nz
	integer boxes(nx,ny)
	integer layers(nx,ny)
	double precision mask(nx,ny,nz)
	double precision dlayer(nz)
	character*80 file

	integer id,ncid,ix,iy,irec
	integer idx,idy,idz,idt
	integer nvars,varid,ivt,nit
	integer tlen,rlen,len,ndims
	integer lmaxs(nbox)
        type(var_item) :: vitem
        type(att_item) :: aitem
        type(dim_item) :: ditem

	double precision, allocatable :: record(:,:,:)
	double precision, allocatable :: records(:,:,:,:)
	double precision, allocatable :: times(:)

	ncid = nitem%ncid
	ditem = nitem%ditem
	file = nitem%filename
	idx = 0
	idy = 0
	idz = 0
	idt = 0

	write(6,*) 'reading netcdf file: ',trim(file)

	do id=1,ditem%ndims
	  if( ditem%name(id) == 'x' ) idx = ditem%len(id)
	  if( ditem%name(id) == 'y' ) idy = ditem%len(id)
	  if( ditem%name(id) == 'deptht' ) idz = ditem%len(id)
	  if( ditem%name(id) == 'time_counter' ) idt = ditem%len(id)
	end do

	write(6,*) 'coordinate dimensions: ',idx,idy,idz,idt
	if( nx /= idx .or. ny /= idy ) stop 'error stop internal (9)'
	if( nz /= idz ) stop 'error stop internal (19)'
	if( idt == 0 ) stop 'error stop: no time dimension'

	allocate(record(idx,idy,idz))

	call make_box_layers(nbox,nx,ny,boxes,layers,lmaxs)

        nvars = nitem%nvars

	ivt = 0
        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
	  if( vitem%name == 'time_counter' ) ivt = varid
	end do
	if( ivt == 0 ) stop 'error stop: no time variable found'
	write(6,*) 'time variable found: ',ivt

	call ncf_init_time(ncid,ivt,nit)
	write(6,*) 'idt,nit: ',idt,nit
	if( idt /= nit ) stop 'error stop: idt /= nit'

        nvars = nitem%nvars
        write(6,*) 'total number of variables: ',nvars
        call ncf_print_variable_header
        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          call ncf_print_variable(vitem)
	  ndims = vitem%ndims
	  if( ndims == 4 ) then
	    call interpolate(nitem,nbox,ivt,varid,vitem,nx,ny
     +			,boxes,layers,dlayer,lmaxs
     +			,mask,idz,idt,record)
	  end if
	  !len = vitem%len
	  !rlen = vitem%rlen
	  !tlen = vitem%tlen
	  !write(6,*) 'var dimensions: ',tlen,rlen,len
        end do

	end

!*************************************************************************

	subroutine interpolate(nitem,nbox,ivt,varid,vitem,nx,ny
     +			,boxes,layers,dlayer,lmaxs
     +			,mask,nz,nt,record)

	use ncf

	implicit none

        type(nc_item) :: nitem
	integer nbox
	integer ivt		!time variable
	integer varid
        type(var_item) :: vitem
	integer nx,ny
	double precision mask(nx,ny,nz)
	integer boxes(nx,ny)
	integer layers(nx,ny)
	double precision dlayer(nz)
	integer lmaxs(nbox)
	integer nz,nt
	double precision record(nx,ny,nz)

        type(att_item) :: aitem
        type(dim_item) :: ditem

	integer ncid,irec,iover,itot
	integer ix,iy,lmax,l,ibox,iu
	integer ius(nbox)
	double precision dtime,atime
	double precision vol
	double precision values(nbox,nz)
	double precision volumes(nbox,nz)
	double precision values2d(nbox)
	character*80 name,format
	character*20 dline

	format = '(a,i10,1f18.2,2x,a20)'

	ncid = nitem%ncid
	ditem = nitem%ditem

	name = vitem%name
	write(6,*) 'interpolation for ',trim(name)
	call make_file_units(name,nbox,ius)
	call write_file_headers(name,nbox,ius,lmaxs)

	do irec = 1,nt
	  call ncf_get_record(ncid,varid,irec,record)
	  call ncf_get_time(ncid,irec,atime)
	  call dts_format_abs_time(atime,dline)
	  write(6,format) ' time: ',irec,atime,dline
	  call check_over(nx,ny,nz,layers,record,mask)
	  values = 0.
	  volumes = 0.
	  do ix=1,nx
	    do iy=1,ny
	      ibox = boxes(ix,iy)
	      lmax = layers(ix,iy)
	      if( ibox == 0 ) cycle
	      if( lmax == 0 ) cycle
	      do l=1,lmax
		vol = dlayer(l)
		values(ibox,l) = values(ibox,l) + vol * record(ix,iy,l)
		volumes(ibox,l) = volumes(ibox,l) + vol
	      end do
	    end do
	  end do

	  do ibox = 1,nbox
	    lmax = lmaxs(ibox)
	    values2d(ibox) = sum(values(ibox,1:lmax)) 
     +				/ sum(volumes(ibox,1:lmax))
	  end do
	  where( volumes > 0 ) values = values / volumes

	  do ibox=1,nbox
	    iu = ius(ibox)
	    lmax = lmaxs(ibox)
	    call write_file_values(iu,irec,dline,nbox,ibox,lmax
     +					,values2d,values)
	  end do

	end do

	end

!*************************************************************************
!*************************************************************************
!*************************************************************************

	subroutine check_over(nx,ny,nz,layers,record,mask)

	implicit none

	integer nx,ny,nz
	integer layers(nx,ny)
	double precision record(nx,ny,nz)
	double precision mask(nx,ny,nz)

	integer itot,iover,ilast
	integer l,lmax,ix,iy

	iover = 0
	itot = 0
	ilast = 0

	do ix=1,nx
	  do iy=1,ny
	      lmax = layers(ix,iy)
	      if( lmax == 0 ) cycle
	      itot = itot + 1
	      if( lmax == nz ) cycle
	      l = lmax + 1
	      if( record(ix,iy,l) /= 0. ) then
		if( l == nz ) then
		  ilast = ilast + 1
		else
		  iover = iover + 1
		  call write_layers(nx,ny,nz,ix,iy,lmax,record,mask)
		end if
	      end if
	  end do
	end do

	if( iover == 0 ) return

	if( iover > 0 .or. ilast > 0 ) then
	    write(6,*) '*** iover = ',iover,ilast,itot,nx*ny
	    !stop 'error stop: iover > 0'
	end if

	end

!*************************************************************************

	subroutine write_layers(nx,ny,nz,ix,iy,lmax,record,mask)

	implicit none

	integer nx,ny,nz,ix,iy,lmax
	double precision mask(nx,ny,nz)
	double precision record(nx,ny,nz)

	integer l,lstart,lend,iu
	integer, parameter :: lwin = 5

	iu = 666
	write(iu,*) 'point: ',ix,iy,lmax
	lstart = max(1,lmax-lwin)
	lend = min(nz,lmax+lwin)
	do l=lstart,lend
	  write(iu,*) l,record(ix,iy,l),mask(ix,iy,l)
	end do

	end

!*************************************************************************
!*************************************************************************
!*************************************************************************

	subroutine make_file_units(name,nbox,ius)

	implicit none

	character*(*) name
	integer nbox
	integer ius(nbox)

	integer iu,iubase,ibox
	character*80 file
	character*2 box

	iubase = 100

	do ibox=1,nbox
	  iu = iubase + ibox
	  write(box,'(i2)') ibox
	  if( box(1:1) == ' ' ) box(1:1) = '0'
	  file = 'extract_' // trim(name) // '_' // box // '.txt'
	  !write(6,*) 'filename: ',trim(file)
	  open(iu,file=file,status='unknown',form='formatted')
	  ius(ibox) = iu
	end do

	end

!*************************************************************************

	subroutine write_file_headers(name,nbox,ius,lmaxs)

	implicit none

	character*(*) name
	integer nbox
	integer ius(nbox)
	integer lmaxs(nbox)

	integer ibox,iu,l,lmax
	character*19 varname
	character*1024 header
	character*2 layer

	do ibox=1,nbox
	  iu = ius(ibox)
	  lmax = lmaxs(ibox)
	  varname = name
	  varname = adjustr(varname)
	  write(iu,'(a)') '#           var_name         box        lmax'
	  write(iu,'(a1,a19,2i12)') '#',varname,ibox,lmax
	  !         12345678901234567890
	  !header = '#   irec     average'
	  header = '#          date_time     average'
	  do l=1,lmax
	    write(layer,'(i2)') l
	    if( layer(1:1) == ' ' ) layer(1:1) = '0'
	    !                         12345678901234567890
	    header = trim(header) // '    layer_' // layer
	  end do
	  write(iu,'(a)') trim(header)
	end do

	end

!*************************************************************************

	subroutine write_file_values(iu,irec,dline,nbox,ibox,lmax
     +					,values2d,values)

	implicit none

	integer iu
	integer irec
	character*20 dline
	integer nbox
	integer ibox
	integer lmax
	double precision values(nbox,lmax)
	double precision values2d(nbox)

	integer l
	character*2, parameter :: zmax = '99'
	character*80 format

	!zmax = '99'

	!format='(i8,' // zmax // 'e12.4)'
	!write(iu,format) irec,values2d(ibox),(values(ibox,l),l=1,lmax)

	format='(a20,' // zmax // 'e12.4)'
	write(iu,format) dline,values2d(ibox),(values(ibox,l),l=1,lmax)

	end

!*************************************************************************
!*************************************************************************
!*************************************************************************

	subroutine write_grid(file,itype,nx,ny,xx,yy,value)

	implicit none

	character*(*) file
	integer itype
	integer nx,ny
	double precision xx(nx,ny)
	double precision yy(nx,ny)
	double precision value(nx,ny)

	real, parameter :: flag = -999.
	integer iu,node,ix,iy,it
	real x,y,d

	it = itype
	iu = 1
	node = 0
	open(iu,file=file,status='unknown',form='formatted')
	do iy=1,ny
	  do ix=1,nx
	    node = node + 1
	    x = xx(ix,iy)
	    y = yy(ix,iy)
	    d = value(ix,iy)
	    if( itype < 0 ) it = nint(d)
	    if( d /= flag ) then
	      write(iu,'(i1,2i8,3e16.8)') 1,node,it,x,y,d
	    end if
	  end do
	end do
	close(iu)

	write(6,*) 'file writen: ',trim(file)

	end

!*************************************************************************

