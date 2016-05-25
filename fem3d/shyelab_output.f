
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

	integer, save, allocatable :: il(:)
	real, save, allocatable :: hd(:)
	real, save, allocatable :: svalue(:,:)
	real, save, allocatable :: s2dvalue(:)
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

	integer ierr,iunit,n
	integer nxy
	logical bshy
	character*60 file,form
	character*80 string
	integer ifileo

	if( .not. boutput ) return

	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	string = ' '
	call fem_regular_parse(string,regpar,nxreg,nyreg)
	if( nxreg > 0 .and. nyreg > 0 ) then
	  nxy = nxreg * nyreg
	  allocate(svalue(nlvdi,nxy),s2dvalue(nxy),fm(4,nxy))
	  call fem_regular_setup_fm(nxreg,nyreg,fm)
	  breg = .true.
	  ntype = ntype + 10
	else
	  allocate(svalue(nlvdi,nkn),s2dvalue(nkn)) !only for nodal values
	end if

	if( breg .and. ( bsplit .or. bshy ) ) then
	  write(6,*) 'regular output only for format gis, fem, nc'
	  stop 'error stop shyelab_init_output: regular output'
	end if

	if( breg .and. outformat == 'gis' ) then
	  write(6,*) 'regular output for gis not yet ready'
	  stop 'error stop shyelab_init_output: gis not ready'
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
	    n = nkn
            allocate(il(n),hd(n))
	    call fem_setup_arrays(id,n,il,hd)		!only for nodes
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

	integer ierr,iunit,n,nvers,lmax

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
	    n = nkn
	    if( breg ) n = nxreg*nyreg
	    lmax = nlv
	    if( b2d ) lmax = 1
            call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,n,lmax
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
	integer it,nb
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
	  call fem_regular_interpolate(nlvdi,cv3,svalue)
	else
	  svalue(:,1:nkn) = cv3(:,1:nkn)
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
	    if( ivar == 1 .or. ivar == 3 ) goto 99
	    if( n /= nkn .or. m /= 1 ) goto 98
	    nb = 0
	    it = nint(dtime)
            call gis_write_record(nb,it,ivar,nlvddi,ilhkv,svalue)
	  else if( outformat == 'fem' ) then
	    iunit = idout
	    nvers = 0
	    n = nkn
	    if( breg ) n = nxreg*nyreg
	    call get_string_description(ivar,string)
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,n,lmax
     +                          ,string
     +                          ,il,hd
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
	integer it,nb,ivar
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
	    it = nint(dtime)
	    if( breg ) stop 'ggguuu: gis not yet ready'
	    call gis_write_hydro(it,nlvddi,ilhkv,znv,uprv,vprv)
	  else if( outformat == 'fem' ) then
	    call fem_write_hydro(idout,dtime,nlvddi,znv,uprv,vprv)
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

	subroutine fem_setup_arrays(id,n,il,hd)

	use basin, only : nkn,nel,nen3v
	use shyfile

	implicit none

	integer id
	integer n
	integer il(n)
	real hd(n)

	integer ie,ii,k
	integer ilhv(nel)
	integer ilhkv(nkn)
	real hm3v(3,nel)

	if( n /= nkn .and. n /= nel ) then
	  stop 'error stop fem_setup_arrays: internal error'
	end if

	  call shy_get_layerindex(id,ilhv,ilhkv)
	  if( n == nkn ) then
	    il = ilhkv
	  else if( n == nel ) then
	    il = ilhv
	  end if
	  call shy_get_depth(id,hm3v)

	if( n == nel ) then
	  do ie=1,nel
	    hd(ie) = sum(hm3v(:,ie))/3.
	  end do
	else
	  hd = -1.e+30
	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hd(k) = max(hd(k),hm3v(ii,ie))
	    end do
	  end do
	end if

	end

!***************************************************************

	subroutine fem_write_hydro(idout,dtime,nlvddi,znv,uprv,vprv)

	use basin
	use levels
	use shyfile
	use elabutil
	use elabtime
	use shyelab_out

	implicit none

	integer idout
	double precision dtime
	integer nndim
	integer nlvddi
	real znv(nkn)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)

	integer iunit,nvers,n,lmax,nvar,ivar
	character*60 string

	iunit = idout
	nvers = 0
	n = nkn
	if( breg ) n = nxreg*nyreg
	lmax = nlv
	nvar = 3

        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,n,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime_elab,regpar)

	ivar = 1
	lmax = 1
	if( breg ) then
	  call fem_regular_interpolate(lmax,znv,s2dvalue)
	else
	  s2dvalue = znv
	end if
	call get_string_description(ivar,string)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,n,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,lmax,s2dvalue)

	ivar = 2
	lmax = nlv
	call get_string_description(ivar,string)

	if( breg ) then
	  call fem_regular_interpolate(nlvddi,uprv,svalue)
	else
	  svalue = uprv
	end if

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,n,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,nlvddi,svalue)

	if( breg ) then
	  call fem_regular_interpolate(nlvddi,vprv,svalue)
	else
	  svalue = vprv
	end if

        call fem_file_write_data(iformat,iunit
     +                          ,nvers,n,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,nlvddi,svalue)

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
        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy

	regpar = (/float(nx),float(ny),dx,dy,x0,y0,flag/)

	call setgeo(x0,y0,dx,dy,flag)

	end

!***************************************************************

	subroutine fem_regular_setup_fm(nx,ny,fm)

	implicit none

	integer nx,ny
	real fm(4,nx,ny)

	call av2fm(fm,nx,ny)

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

