
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shyelab_write_record(id,idout,dtime,nvar,iv,ivar,n,m
     +						,nlv,nlvddi,cv3)

	use basin
	use shyfile
	use elabutil

	implicit none

	integer id,idout
	double precision dtime
	integer nvar,iv,ivar,n,m
	integer nlv,nlvddi
	real cv3(nlvddi,n*m)

	integer it,nb
	integer iunit,nvers
	integer, save, allocatable :: il(:)
	real, save, allocatable :: hlv(:)
	real, save, allocatable :: hd(:)

	integer ilhv(nel)
	integer ilhkv(nkn)

	integer, save :: iformat = 0
	integer, save :: ntype = -1
	real, save :: regpar(7) = 0.
	integer, save :: icall = 0
	character*60 string

	nb = 0
	it = nint(dtime)

	if( outformat == 'shy' .or. outformat == 'native') then
	  call shy_write_output_record(idout,dtime,ivar,n,m
     +						,nlv,nlv,cv3)
        else if( outformat == 'gis' ) then
	  if( ivar == 1 .or. ivar == 3 ) goto 99
	  if( n /= nkn .or. m /= 1 ) goto 98
	  call shy_get_layerindex(id,ilhv,ilhkv)
          call gis_write_record(nb,it,ivar,nlvddi,ilhkv,cv3)
        else if( outformat == 'fem' ) then
	  if( icall == 0 ) then
	    call fem_setup_params(iformat,idout,ntype,regpar)
	    allocate(hlv(nlvddi),il(n),hd(n))
	    call fem_setup_arrays(id,nlvddi,n,hlv,il,hd)
	  end if
	  if( n /= nkn .or. m /= 1 ) goto 98
	  iunit = idout
	  nvers = 0
	  if( iv == 1 ) then
            call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,n,nlv
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
	  end if
	  call get_string_description(ivar,string)
          call fem_file_write_data(iformat,iunit
     +                          ,nvers,n,nlv
     +                          ,string
     +                          ,il,hd
     +                          ,nlvddi,cv3)
	else
          write(6,*) 'output format unknown: ',outformat
          stop 'error stop shyelab_write_record: output format'
	end if

	icall = icall + 1

	return
   98	continue
	write(6,*) 'output format = ',trim(outformat)
	write(6,*) 'n,m = ',n,m
	write(6,*) 'nkn,nel = ',nkn,nel
	stop 'error stop shyelab_write_record: cannot handle'
   99	continue
	write(6,*) 'output format = ',trim(outformat)
	write(6,*) 'ivar = ',ivar
	stop 'error stop shyelab_write_record: cannot handle'
	end

!***************************************************************

	subroutine fem_setup_params(iformat,iunit,ntype,regpar)

	implicit none

	integer iformat
	integer iunit,ntype
	real regpar(7)

	character*20 form
	integer ifileo

	ntype = 0
	regpar = 0.

	if( iformat == 0 ) then
	  form='unformatted'
	else
	  form='formatted'
	end if

	iunit = ifileo(60,'out.fem',form,'unknown')

	end

!***************************************************************

	subroutine fem_setup_arrays(id,nlvddi,n,hlv,il,hd)

	use basin, only : nkn,nel,nen3v
	use shyfile

	implicit none

	integer id
	integer nlvddi
	integer n
	real hlv(nlvddi)
	integer il(n)
	real hd(n)

	integer ie,ii,k
	integer ilhv(nel)
	integer ilhkv(nkn)
	real hm3v(3,nel)

	if( n /= nkn .and. n /= nel ) then
	  stop 'error stop fem_setup_arrays: internal error'
	end if

	if( nlvddi == 1 ) then
	  il = 1
	  hlv(1) = 10000.
	  call shy_get_depth(id,hm3v)
	else
	  call shy_get_layers(id,hlv)
	  call shy_get_layerindex(id,ilhv,ilhkv)
	  if( n == nkn ) then
	    il = ilhkv
	  else if( n == nel ) then
	    il = ilhv
	  end if
	  call shy_get_depth(id,hm3v)
	end if

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

	subroutine fem_regular_interpolate(nx,ny,nlvddi
     +					,nlv,ilhv,cv3,fm,am)

        use basin

        implicit none

	integer nx,ny
        integer nlvddi                  !vertical dimension of fem array
        integer nlv                     !vertical dimension of regular array
        integer ilhv(nel)               !vertical discretization (element!!)
        real cv3(nlvddi,nkn)            !values of fem array
        real fm(4,nx,ny)                !interpolation matrix
        real am(nlv,nx,ny)              !interpolated values (return)

	call fm2am3d(nlvddi,ilhv,cv3,nlv,nx,ny,fm,am)

	end

!***************************************************************

