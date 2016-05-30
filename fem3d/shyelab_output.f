
!==================================================================
        module shyelab_out
!==================================================================

	implicit none


	integer, save :: nwrite = 0

	integer, save :: iformat = 1
	integer, save :: ntype = 1
	integer, save :: nxreg = 0
	integer, save :: nyreg = 0
	real, save :: regpar(7) = 0.
	logical, save :: breg = .false.

	real, save, allocatable :: s2dvalue(:)
	real, save, allocatable :: zvalue(:)
	real, save, allocatable :: svalue(:,:)
	real, save, allocatable :: uvalue(:,:)
	real, save, allocatable :: vvalue(:,:)
	integer, save, allocatable :: ilcoord(:)
	real, save, allocatable :: xcoord(:)
	real, save, allocatable :: ycoord(:)
	real, save, allocatable :: hcoord(:)
	real, save, allocatable :: fm(:,:)

!==================================================================
        end module shyelab_out
!==================================================================

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shyelab_init_output(id,idout)

! initializes in case of btrans, bout, bsplit, bsumvar
!
! in case of bout also distinguishes between output formats

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout

	integer ierr,iunit,np
	integer nxy
	logical bshy
	character*60 file,form
	character*80 string
	integer ifileo

	if( .not. boutput ) return

	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	string = regstring
	call fem_regular_parse(string,regpar,nxreg,nyreg)
	if( nxreg > 0 .and. nyreg > 0 ) then
	  nxy = nxreg * nyreg
	  allocate(svalue(nlvdi,nxy),s2dvalue(nxy))
	  allocate(zvalue(nxy),uvalue(nlvdi,nxy),vvalue(nlvdi,nxy))
	  allocate(ilcoord(nxy),xcoord(nxy),ycoord(nxy),hcoord(nxy))
	  allocate(fm(4,nxy))
	  call fem_regular_setup(nxreg,nyreg,regpar,ilhv,hm3v
     +				,fm,ilcoord,xcoord,ycoord,hcoord)
	  breg = .true.
	  ntype = ntype + 10
	else
	  allocate(svalue(nlvdi,nkn),s2dvalue(nkn)) !only for nodal values
	  allocate(zvalue(nkn),uvalue(nlvdi,nkn),vvalue(nlvdi,nkn))
	  allocate(ilcoord(nkn),xcoord(nkn),ycoord(nkn),hcoord(nkn))
	  xcoord = xgv
	  ycoord = ygv
	  ilcoord = ilhkv
	  call makehkv_minmax(hcoord,+1)
	end if

	if( breg .and. ( bsplit .or. bshy ) ) then
	  write(6,*) 'regular output only for format gis, fem, nc'
	  stop 'error stop shyelab_init_output: regular output'
	end if

	if( breg .and. outformat == 'gis' ) then
	  write(6,*) 'regular output for gis not yet ready'
	  stop 'error stop shyelab_init_output: gis not ready'
	end if

	if( breg .and. outformat == 'nc' ) then
	  write(6,*) 'regular output for nc not yet ready'
	  stop 'error stop shyelab_init_output: nc not ready'
	end if

	if( bsplit ) then
	  ! nothing to initialize ... is done later
	else
	  if( outformat == 'shy' .or. outformat == 'native') then
	    file = 'out.shy'
            idout = shy_init(file)
	    if( idout <= 0 ) goto 74
            call shy_clone(id,idout)
            if( b2d ) call shy_convert_2d(idout)
	    if( bsumvar ) call shy_convert_1var(idout)
            call shy_write_header(idout,ierr)
            if( ierr /= 0 ) goto 75
	  else if( outformat == 'gis' ) then
	    call gis_write_connect
	  else if( outformat == 'fem' ) then
	    if( iformat == 0 ) then
	      form='unformatted'
	    else
	      form='formatted'
	    end if
	    file = 'out.fem'
	    iunit = ifileo(60,file,form,'unknown')
	    if( iunit <= 0 ) goto 74
	    idout = iunit
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	return
   74	continue
        write(6,*) 'error opening file ',trim(file)
        stop 'error stop shyelab_init_output: opening file'
   75	continue
        write(6,*) 'error writing header, ierr = ',ierr
        write(6,*) 'file = ',trim(file)
        stop 'error stop shyelab_init_output: writing header'
	end

!***************************************************************

	subroutine shyelab_header_output(id,idout,dtime,nvar)

! initializes record in case of btrans, bout, bsplit, bsumvar

	use basin
	use levels
	use elabutil
	use elabtime
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	double precision dtime
	integer nvar,ftype
	logical bscalar,bhydro

	integer ierr,iunit,np,nvers,lmax

	if( .not. boutput ) return

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2

	if( bhydro ) return

	if( bsplit ) then
	  ! nothing to be done
	else
	  if( outformat == 'shy' .or. outformat == 'native') then
	    ! nothing to be done
	  else if( outformat == 'gis' ) then
	    ! nothing to be done
	  else if( outformat == 'fem' ) then
	    iunit = idout
	    nvers = 0
	    np = nkn
	    if( breg ) np = nxreg*nyreg
	    lmax = nlv
	    if( b2d ) lmax = 1
            call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime_elab,regpar)
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	end

!***************************************************************

	subroutine shyelab_record_output(id,idout,dtime,ivar,n,m
     +					,lmax,nlvddi,cv3)

! initializes record in case of btrans, bout, bsplit, bsumvar

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	double precision dtime
	integer ivar,n,m
	integer lmax,nlvddi
	real cv3(nlvddi,n*m)

	integer ierr,iunit,nvers
	integer np
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy
	character*60 string

	if( .not. boutput ) return
	if( bsumvar ) return

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bhydro .and. .not. bshy ) return

	if( bverb ) write(6,*) 'writing to output: ',ivar

	if( breg ) then
	  if( n /= nkn ) goto 98
	  np = nxreg * nyreg
	  call fem_regular_interpolate(nlvdi,cv3,svalue)
	else if( n == nkn ) then
	  np = nkn
	  svalue(:,1:nkn) = cv3(:,1:nkn)
	else if( n == nel ) then
	  !convert from nel to nkn
	  stop 'error stop shyelab_record_output: n==nel not ready'
	else if( .not. bsplit .and. .not. bshy ) then
	  goto 98
	end if

	if( .not. bsplit .and. .not. bshy ) then
	  if( ivar == 1 .or. ivar == 3 ) goto 99
	  if( m /= 1 ) goto 98
	end if

	if( bsplit ) then
	  call shy_split_id(ivar,id,id_out)
	  call shy_write_output_record(id_out,dtime,ivar,n,m
     +						,lmax,nlvddi,cv3)
	else
	  if( bshy ) then
	    call shy_write_output_record(idout,dtime,ivar,n,m
     +						,lmax,nlvddi,cv3)
	  else if( outformat == 'gis' ) then
            call gis_write_record(dtime,ivar,np,nlvddi,ilcoord
     +					,svalue,xcoord,ycoord)
	  else if( outformat == 'fem' ) then
	    iunit = idout
	    nvers = 0
	    call get_string_description(ivar,string)
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilcoord,hcoord
     +                          ,nlvddi,svalue)
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	nwrite = nwrite + 1

	return
   98   continue
        write(6,*) 'output format = ',trim(outformat)
        write(6,*) 'n,m = ',n,m
        write(6,*) 'nkn,nel = ',nkn,nel
        stop 'error stop shyelab_record_output: cannot handle'
   99   continue
        write(6,*) 'output format = ',trim(outformat)
        write(6,*) 'ivar = ',ivar
        stop 'error stop shyelab_record_output: cannot handle'
	end

!***************************************************************

	subroutine shyelab_post_output(id,idout,dtime,nvar,n,m,nndim
     +					,lmax,nlvddi,cv3all)

! initializes record in case of btrans, bout, bsplit, bsumvar

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	double precision dtime
	integer nvar,n,m,nndim
	integer lmax,nlvddi
	real cv3all(nlvddi,nndim,0:nvar)

	integer ierr,iunit,nvers
	integer it,nb,ivar,np
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy
	character*60 string

        real, allocatable :: znv(:)
        real, allocatable :: uprv(:,:)
        real, allocatable :: vprv(:,:)
        real, allocatable :: sv(:,:)
        real, allocatable :: dv(:,:)

	if( .not. boutput ) return

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bscalar .and. .not. bsumvar ) return

	if( bhydro ) then
	  allocate(znv(nkn),uprv(nlvddi,nkn),vprv(nlvddi,nkn))
	  allocate(sv(nlvddi,nkn),dv(nlvddi,nkn))
          call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
          call convert_to_speed(uprv,vprv,sv,dv)
	  if( breg ) then
	    np = nxreg * nyreg
	    call fem_regular_interpolate(1,znv,zvalue)
	    call fem_regular_interpolate(nlvddi,uprv,uvalue)
	    call fem_regular_interpolate(nlvddi,vprv,vvalue)
	  else
	    np = nkn
	    zvalue = znv
	    uvalue = uprv
	    vvalue = vprv
	  end if
	end if	

	if( bsplit ) then
	  if( bhydro ) then
            call shy_split_hydro(id,dtime,znv,uprv,vprv,sv,dv)
	  end if
	else
	  if( bshy ) then
	    if( bhydro ) then
	      ! nothing to be done
	    else if( bsumvar ) then
              cv3all(:,:,0) = 0.
              cv3all(:,:,0) = sum(cv3all,dim=3)
              ivar = 30
	      call shy_write_output_record(idout,dtime,ivar,n,m
     +					,lmax,nlvddi,cv3all(:,:,0))
	    else
	      stop 'error stop shyelab_post_output: internal error (1)'
	    end if
	  else if( outformat == 'gis' ) then
	    call gis_write_hydro(dtime,np,nlvddi,ilcoord
     +				,zvalue,uvalue,vvalue,xcoord,ycoord)
	  else if( outformat == 'fem' ) then	!also covers breg
	    call fem_write_hydro(idout,dtime,np,nlvddi
     +					,zvalue,uvalue,vvalue)
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
          if( .not. (bshy.and.bhydro) ) nwrite = nwrite + 1
	end if

	if( bhydro ) deallocate(znv,uprv,vprv,sv,dv)

	end

!***************************************************************

	subroutine shyelab_final_output(id,idout,nvar)

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	integer nvar

	integer ierr,iunit,nvers
	integer it,nb,ivar
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy
	character*80 file

	if( .not. boutput ) return

	write(6,*) 'output written to following files: '

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bsplit ) then
	  if( bhydro ) then
	    write(6,*) 'z.shy'
	    write(6,*) 'u.shy'
	    write(6,*) 'v.shy'
	    write(6,*) 's.shy'
	    write(6,*) 'd.shy'
	  else if( bscalar ) then
	    ivar = 0
	    do
	      ivar = ivar - 1
	      call shy_split_id(ivar,id,idout)
	      if( idout < 0 ) exit		!no more variables
	      if( idout > 0 ) then
		call shy_get_filename(idout,file)
	        write(6,*) trim(file)
	      end if
	    end do
	  end if
	else
	  if( bshy ) then
	    write(6,*) 'out.shy'
	  else if( outformat == 'gis' ) then
	    write(6,*) 'connectivity.gis'
	    write(6,*) 'extract_*.gis'
	  else if( outformat == 'fem' ) then
	    write(6,*) 'out.fem'
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	end

!***************************************************************
!***************************************************************
!***************************************************************


!***************************************************************
!***************************************************************
!***************************************************************

	subroutine fem_write_hydro(idout,dtime,np,nlvddi,zv,uv,vv)

	use levels
	use elabtime
	use shyelab_out

	implicit none

	integer idout
	double precision dtime
	integer np
	integer nndim
	integer nlvddi
	real zv(np)
	real uv(nlvddi,np)
	real vv(nlvddi,np)

	integer iunit,nvers,lmax,nvar,ivar
	character*60 string

	iunit = idout
	nvers = 0
	lmax = nlv
	nvar = 3

        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime_elab,regpar)

	ivar = 1
	lmax = 1
	call get_string_description(ivar,string)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilcoord,hcoord
     +                          ,lmax,zv)

	ivar = 2
	lmax = nlv
	call get_string_description(ivar,string)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilcoord,hcoord
     +                          ,nlvddi,uv)

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilcoord,hcoord
     +                          ,nlvddi,vv)

	end

!***************************************************************

	subroutine get_string_description(iv,string)

	implicit none

	integer iv
	character*(*) string

	call ivar2string(iv,string)

	if( string == ' ' ) then
	  stop 'error stop get_string_description: unknown ivar'
	end if

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine fem_regular_parse(string,regpar,nx,ny)

	use basin

	implicit none

	character*(*) string
	real regpar(7)
	integer nx,ny

	integer ianz
	real dx,dy,x0,y0,x1,y1
	real ddx,ddy
	double precision d(6)
	integer iscand
	real, parameter :: flag = -999.

	regpar = 0.
	ianz = iscand(string,d,6)

	if( ianz == 0 ) then
	  nx = 0
	  ny = 0
	else if( ianz == 1 ) then
	  dx = d(1)
	  dy = dx
	else if( ianz == 2 ) then
	  dx = d(1)
	  dy = d(2)
	else if( ianz == 6 ) then
	  dx = d(1)
	  dy = d(2)
	else
	  write(6,*) ianz,'  ',trim(string)
	  write(6,*) 'read error or wrong number of parameters'
	  stop 'error stop fem_regular_parse: string'
	end if

	if( ianz == 0 ) then
	  return
	else if( ianz == 6 ) then
	  x0 = d(3)
	  y0 = d(4)
	  x1 = d(5)
	  y1 = d(6)
	else
	  x0 = minval(xgv)
	  y0 = minval(ygv)
	  x1 = maxval(xgv)
	  y1 = maxval(ygv)
	end if

        nx = 1 + nint((x1-x0)/dx)
        ny = 1 + nint((y1-y0)/dy)
	if( ianz == 2 ) then		!correct x0/y0
	  ddx = (nx-1)*dx - (x1-x0)
	  if( ddx > 0. ) x0 = x0 - ddx/2.
	  ddy = (ny-1)*dy - (y1-y0)
	  if( ddy > 0. ) y0 = y0 - ddy/2.
	end if
        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy

	regpar = (/float(nx),float(ny),dx,dy,x0,y0,flag/)

	call setgeo(x0,y0,dx,dy,flag)

	end

!***************************************************************

	subroutine fem_regular_setup(nx,ny,regpar,ilhv,hm3v
     +				,fm,ilcoord,xcoord,ycoord,hcoord)

	implicit none

	integer nx,ny
	real regpar(7)
	integer ilhv(*)
	real hm3v(3,*)
	real fm(4,nx,ny)
	integer ilcoord(nx,ny)
	real xcoord(nx,ny)
	real ycoord(nx,ny)
	real hcoord(nx,ny)

	integer ix,iy,ie
	real dx,dy,x0,y0,x,y

	call av2fm(fm,nx,ny)

	ilcoord = 0
	hcoord = 0.
	
	dx = regpar(3)
	dy = regpar(4)
	x0 = regpar(5)
	y0 = regpar(6)

	do iy=1,ny
	  y = y0 + (iy-1)*dy
	  do ix=1,nx
	    x = x0 + (ix-1)*dx
	    xcoord(ix,iy) = x
	    ycoord(ix,iy) = y
	    ie = fm(4,ix,iy)
	    if( ie > 0 ) then
	      ilcoord(ix,iy) = ilhv(ie)
	      hcoord(ix,iy) = sum(hm3v(:,ie))/3.	!average in element
	    end if
	  end do
	end do

	end

!***************************************************************

	subroutine fem_regular_interpolate(lmax,cv3,am)

        use basin
	use levels
	use shyelab_out

        implicit none

        integer lmax                   !vertical dimension of regular array
        real cv3(nlvdi,nkn)            !values of fem array
        real am(nlvdi,nxreg*nyreg)     !interpolated values (return)

	real fem2d(nkn)
	real am2d(nxreg*nyreg)

	if( lmax <= 1 ) then
	  fem2d = cv3(1,:)
	  call fm2am2d(fem2d,nxreg,nyreg,fm,am2d)
	  am(1,:) = am2d
	else
	  call fm2am3d(nlvdi,ilhv,cv3,nlvdi,nxreg,nyreg,fm,am)
	end if

	end

!***************************************************************

