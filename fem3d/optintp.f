
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c revision log :
c
c 20.08.2003	ggu	new laplacian interpolation
c 02.09.2003	ggu	some comments, write to .dat file
c 30.10.2003	ggu	subroutine prepare_bc_l included in this file
c 04.03.2004	ggu	writes also number of variables (1)
c 11.03.2009	ggu	bug fix -> declare hev() here
c 28.01.2014	ggu	changed VERS_6_1_71
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 07.02.2015	ggu	OI finished
c 26.02.2015	ggu	changed VERS_7_1_5
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 14.09.2015	ggu	OI adapted to new modular structure
c 29.09.2015	ggu	move some vars to custimize to the beginning of file
c 26.09.2017	ggu	changed VERS_7_5_32
c 26.03.2018	ggu	completely restructured for model driven obs
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c
c notes :
c
c header_record_1
c data_record_1
c header_record_2
c data_record_2
c etc...
c
c where header_record is
c
c time,nobs
c
c and data_record is either:
c
c k_1	val_1
c k_2	val_2
c ...
c k_n	val_n
c
c or alternatively:
c
c x_1 y_1    val_1
c x_2 y_2    val_2
c ...
c x_n y_n    val_n
c
c legend:
c
c time		time in seconds (relative time to a reference time)
c nobs		number of observations in data record
c k		node (external) of observation in basin
c x,y		coordinates of observation
c val		value of observation
c
c****************************************************************

        program optintp

c optimal interpolation interpolation

	use mod_depth
	use evgeom
	use basin
	use clo

	implicit none

	integer nobdim
	parameter (nobdim = 29)

	real xobs(nobdim)
	real yobs(nobdim)
	real zobs(nobdim)
	real bobs(nobdim)
	real errors(nobdim)
	real rra(nobdim)
	real ssa(nobdim)
	integer iuse(nobdim)

	real, allocatable :: xback(:)
	real, allocatable :: yback(:)
	real, allocatable :: zback(:)
	real, allocatable :: zanal(:)
	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)

	integer nvers,ntype,nvar,nlvdi
	integer iformat,lmax,np
	integer iufem,iuobs,iunos
	integer jmax,j
	double precision dtime,dtime0,atime,atime0
	real hlv(1)
	real regpar(7)
	integer date,time
	integer datetime(2)
	character*30 string,format

	character*80 file,basnam
	logical bback,bgeo,bcart,bquiet,blimit,breg,bnos
	logical bfem,bloop
	logical bdebug
	integer k,ie,n,ndim
	integer nobs,nback
        integer ilev,ivar,irec,nrec
	integer it,index
	integer iexcl
	integer nx,ny
	real zmin,zmax
	real zomin,zomax
	real xmin,ymin,xmax,ymax
	real dx,dy,dxy
	real drl,rlact
	real rmse,rmse1,rmse2
	real rr,ss

	real rl,rlmax

	integer iapini
	logical is_spherical

	bdebug = .false.
	bfem = .true.
	bfem = .false.
	bloop = .false.
	bloop = .true.

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

        call clo_init('optintp','input-file','1.0')

        call clo_add_info('Optimal interpolation for sparse data')
        call clo_add_option('basin',' ','basin to be used')
        call clo_add_option('rl #',-1.
     +		,'set length scale for covariance matrix')
        call clo_add_option('drl #',-1.
     +		,'try multiple values of rl with step drl')
        call clo_add_option('rlmax #',-1.
     +		,'maximum distance of nodes to be considered')
        call clo_add_option('rr #',-1.
     +		,'std of observation errors')
        call clo_add_option('ss #',-1.
     +		,'std of background field')
        call clo_add_option('dx #',0.
     +		,'dx for regular field')
        call clo_add_option('dy #',0.
     +		,'dy for regular field (Default dy=dx)')
        call clo_add_option('cart',.false.
     +		,'force use of cartesian coordinates')
        call clo_add_option('limit',.false.
     +		,'limit values to min/max of observations')
        call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_extra('defaults for undefined values:')
	call clo_add_extra('  rl=1/10 of basin size, rlmax=10*rl')
	call clo_add_extra('  rr=0.01  ss=100*rr')
	call clo_add_extra('  Default for dx is 0 (use FEM grid)')

        call clo_parse_options(1)  !expecting (at least) 1 file after options

        call clo_get_option('basin',basnam)
        call clo_get_option('cart',bcart)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('rl',rl)
        call clo_get_option('drl',drl)
        call clo_get_option('rlmax',rlmax)
        call clo_get_option('rr',rr)
        call clo_get_option('ss',ss)
        call clo_get_option('limit',blimit)
        call clo_get_option('dx',dx)
        call clo_get_option('dy',dy)

	!call clo_info()

	bback = .false.
	call clo_get_file(1,file)

c-----------------------------------------------------------------
c set some variables
c-----------------------------------------------------------------

	!ivar = 85		!what is in the observations - ice
	!string = 'ice cover [0-1]'
	!date = 19970101

	!ivar = 12		!what is in the observations - temperature
	!string = 'temperature [C]'
	!date = 20100101

	ivar = 1		!what is in the observations - water level
	string = 'water level [m]'
	date = 20130101
	nrec = 10
	nrec = 0

c-----------------------------------------------------------------
c time management
c-----------------------------------------------------------------

	time = 0
	call dtsini(date,time)
	call dts_to_abs_time(date,time,atime0)
	dtime0 = 0.
	dtime = 0.

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	call ap_set_names(basnam,' ')
	call ap_init(.false.,1,0,0)

	allocate(xback(nkn),yback(nkn),zback(nkn),zanal(nkn))
	allocate(ilhkv(nkn),hd(nkn))
	ndim = nkn

	call bas_get_minmax(xmin,ymin,xmax,ymax)
	dxy = max((xmax-xmin),(ymax-ymin))

	if( rl < 0 ) rl = 0.1*dxy
	if( rlmax == 0 ) rlmax = 2.*dxy
	if( rlmax < 0 ) rlmax = 10.*rl
	if( rr < 0 ) rr = 0.01
	if( ss < 0 ) ss = 100.*rr
	if( dy <= 0 ) dy = dx

	rra = rr
	ssa = ss

	if( .not. bquiet ) then
	  write(6,*) 'parameters used:'
	  write(6,*) 'xmin/xmax: ',xmin,xmax
	  write(6,*) 'ymin/ymax: ',ymin,ymax
	  write(6,*) 'dxy      : ',dxy
	  write(6,*) 'rl:        ',rl
	  write(6,*) 'rlmax:     ',rlmax
	  write(6,*) 'rr:        ',rr
	  write(6,*) 'ss:        ',ss
	  write(6,*) 'dx,dy:     ',dx,dy
	  write(6,*) 'limit:     ',blimit
	end if

c-----------------------------------------------------------------
c set up ev and background grid
c-----------------------------------------------------------------

	call ev_init(nel)
	call mod_depth_init(nkn,nel)

	call set_ev
	call check_ev

	bgeo = is_spherical()
	if( bcart ) bgeo = .false.

	regpar = 0.
	nback = nkn
	call setup_background(ndim,dx,dy,xmin,ymin,xmax,ymax
     +			,nback,xback,yback,regpar)
	nx = nint(regpar(1))
	ny = nint(regpar(2))

	breg = dx > 0.
	bnos = .not. breg
	if( .not. breg ) bloop = .false.

	if( .not. bquiet ) then
	  write(6,*) 'background grid:'
	  write(6,*) 'nback: ',nback
	  if( breg ) then
	    write(6,*) 'nx,ny: ',nx,ny
	    write(6,*) 'dx,dy: ',dx,dy
	    write(6,*) 'x0,y0: ',regpar(3),regpar(4)
	  end if
	end if

	if( bdebug ) write(6,*) 'ggu 1'
	iformat = 0
	format = 'unformatted'
	ntype = 1
	if( breg ) then
	  iformat = 1
	  format = 'formatted'
	  ntype = ntype + 10
	end if

	if( bdebug ) write(6,*) 'ggu 2'
	nvers = 0
	nvar = 1
	lmax = 1
	nlvdi = 1
	np = nback
	if( bdebug ) write(6,*) 'ggu 2'
	it = 0
	datetime = (/date,0/)
	if( bdebug ) write(6,*) 'ggu 2'
	hlv(1) = 10000.
	if( bdebug ) write(6,*) 'ggu 2'
	hd = 1.
	ilhkv = 1

	if( bdebug ) write(6,*) 'ggu 3'

c-----------------------------------------------------------------
c open files
c-----------------------------------------------------------------

	irec = 0
	!drl = 0.05		!test different rl, distance drl
	jmax = 0
	if( drl > 0. ) jmax = 5
	if( jmax > 0 ) then
	  write(6,*) 'trying multiple values for rl: ',rl,drl,jmax
	end if

	if( bnos ) then
	  call wrnos2d_open(iunos,'optintp','optimal interpolation')
	end if

	iufem = 20
	open(iufem,file='optintp.fem',status='unknown',form=format)

	iuobs = 21
	open(iuobs,file=file,status='old',form='formatted')
	write(6,*) 'file with observations: ',trim(file)

	iuse = 1
	nobs = nobdim
	call read_use_file(nobs,iuse)
	where( iuse == 0 ) rra = 100.*rr

c-----------------------------------------------------------------
c read observations and interpolate
c-----------------------------------------------------------------

	write(6,*) '                time     ' //
     +		'min/max interpol  min/max observed   rmse'

	iexcl = 0

	do

	if ( nrec > 0 .and. irec == nrec ) exit

	nobs = nobdim
	call read_observations(iuobs
     +				,dtime,nobs,xobs,yobs,zobs)

	if( nobs <= 0 ) exit

	irec = irec + 1
	call mima(zobs,nobs,zomin,zomax)
	!write(6,*) 'observations min/max: ',zomin,zomax

	if( irec == 1 ) then
	  call write_obs_grid(nobs,xobs,yobs,zobs)
	end if

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	do j=-jmax,jmax
	  rlact = rl + j*drl
          call opt_intp(nobs,xobs,yobs,zobs,bobs
     +                  ,nback,bback,xback,yback,zback
     +                  ,rlact,rlmax,ssa,rra,zanal)

	  call mima(zanal,nback,zmin,zmax)
	  if( blimit ) call limit_values(nback,zanal,zomin,zomax)

	  if( bnos ) then
	    it = nint(dtime)
	    call wrnos2d_record(iunos,it,ivar,zanal)
	  end if

	  atime = atime0 + dtime
	  call dts_from_abs_time(date,time,atime)
	  datetime = (/date,time/)

	  if( bfem ) then
	    call fem_file_write_header(iformat,iufem,dtime0
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime,regpar)
            call fem_file_write_data(iformat,iufem
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,zanal)
	  end if
	  if( bloop ) then
	    call compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar
     +				,zanal,0,rmse1)
	    write(67,*) dtime,nobs
            call opt_intp_loop(nobs,xobs,yobs,zobs,bobs
     +                  ,nback,bback,xback,yback,zback,regpar
     +                  ,rlact,rlmax,ssa,rra,zanal,rmse2,errors)
	    !write(6,*) rmse1,rmse2
	    write(66,*) dtime,nobs
	    write(66,*) errors(1:nobs)
	  end if

	  if( bfem .and. .not. bloop ) then	!multiple records in time
	    write(6,3000) 'min/max: ',dtime,zmin,zmax,zomin,zomax
 1000	    format(a,5f12.2)
	  else if( bloop ) then
	    write(6,3000) 'min/max: ',dtime
     +				,zmin,zmax,zomin,zomax,rmse1,rmse2
 3000	    format(a,f12.2,4f9.2,2x,2f9.3)
	  end if

	  it = it + 1
	  dtime = dtime + 1
	end do

	end do

	write(6,*) 'total number of records treated: ',irec
	if( irec == nrec ) then
	  write(6,*) 'limiting records treated to: ',nrec
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine setup_background(ndim,dx,dy,xmin,ymin,xmax,ymax
     +			,nback,xback,yback,regpar)

	use basin

	implicit none

	!include 'param.h'

	integer ndim
	real dx,dy
	real xmin,ymin,xmax,ymax
	integer nback
	real xback(ndim)
	real yback(ndim)
	real regpar(7)

	integer k,i,j
	integer nx,ny
	real x0,y0
	real diff,flag

	x0 = 0.
	y0 = 0.
	regpar = 0.
	flag = -999.

	if( dx <= 0. ) then
	  nx = 0
	  ny = 0
	  nback = nkn
	  if( nback > ndim ) goto 99
	  do k=1,nback
	    xback(k) = xgv(k)
	    yback(k) = ygv(k)
	  end do
	else
	  if( dy < = 0. ) stop 'error stop setup_background: dy'
	  nx = 1 + (xmax - xmin) / dx
	  ny = 1 + (ymax - ymin) / dy
	  diff = nx * dx - (xmax-xmin)
	  x0 = xmin - diff/2.
	  diff = ny * dy - (ymax-ymin)
	  y0 = ymin - diff/2.
	  k = 0
	  do j=0,ny
	    do i=0,nx
	      k = k + 1
	      if( k > ndim ) goto 99
	      xback(k) = x0 + i*dx
	      yback(k) = y0 + j*dy
	    end do
	  end do
	  nx = nx + 1
	  ny = ny + 1
	  nback = nx*ny
	  if( nback .ne. k ) goto 98
	  regpar(1) = nx
	  regpar(2) = ny
	  regpar(3) = x0
	  regpar(4) = y0
	  regpar(5) = dx
	  regpar(6) = dy
	  regpar(7) = flag
	end if

	if( nx > 0. ) then
	  call write_reg_grid(nx,ny,x0,y0,dx,dy)
	  call write_box_grid(nx,ny,x0,y0,dx,dy)
	end if

	return
   98	continue
	write(6,*) nx,ny,nback,k,ndim
	stop 'error stop setup_background: internal error (1)'
   99	continue
	write(6,*) nx,ny,nback,k,ndim
	stop 'error stop setup_background: dimension ndim'
	end

c****************************************************************

	subroutine limit_values(nback,zback,zomin,zomax)

c limits computed values to min/max of observations

	implicit none

	integer nback
	real zback(nback)
	real zomin,zomax

	integer i

	do i=1,nback
	  zback(i) = max(zback(i),zomin)
	  zback(i) = min(zback(i),zomax)
	end do

	end

c****************************************************************

	subroutine read_observations(iuobs
     +				,dtime,nobs,xobs,yobs,zobs)

c reads boundary conditions from file and sets up array
c
c file must be made like this:
c
c	k1, val1
c	k2, val2
c	...
c	kn, valn
c
c or
c
c	x1, y1, val1
c	x2, y2, val2
c	...
c	xn, yn, valn
c

	use basin

	implicit none

	integer iuobs			!unit to be read from
	double precision dtime
	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)

	character*80 line
	logical bdebug
	integer k,kn,ndim,n,ianz,ios
	real f(10)

	integer ipint,iscanf

	bdebug = .false.

	ndim = nobs
	nobs = 0
	n = 0
	dtime = 0.

	if( bdebug ) then
	  write(6,*) '...reading observations from unit :',iuobs
	  write(6,*) '   format: "k val" or "x y val"'
	end if

	read(iuobs,*,iostat=ios) dtime,nobs
	if( ios < 0 ) return
	if( ios > 0 ) goto 95

	do
	  if( nobs > 0 .and. n == nobs ) exit
	  read(iuobs,'(a)',iostat=ios) line
	  if( ios /= 0 ) exit
	  ianz = iscanf(line,f,4)
	  if( ianz < 0 ) ianz = -ianz - 1
	  if( ianz == 0 ) cycle
	  n = n + 1
	  if( n .gt. ndim ) goto 96
	  if( ianz == 2 ) then			!node given
	    k = nint(f(1))
	    kn = ipint(k)
	    if( kn .le. 0 ) goto 99
	    xobs(n) = xgv(kn)
	    yobs(n) = ygv(kn)
	    zobs(n) = f(2)
	    if( bdebug ) write(6,*) n,k,kn,zobs(n)
	  else if( ianz == 3 ) then		!x/y given
	    xobs(n) = f(1)
	    yobs(n) = f(2)
	    zobs(n) = f(3)
	    if( bdebug ) write(6,*) n,xobs(n),yobs(n),zobs(n)
	  else
	    goto 98
	  end if
	end do

	nobs = n

	return
   95	continue
	write(6,*) 'error reading header of record'
	stop 'error stop read_observations: read error in header line'
   96	continue
	write(6,*) n,ndim
	stop 'error stop read_observations: dimensions'
   98	continue
	write(6,*) 'input file line number = ',n
	write(6,*) 'line = ',trim(line)
	write(6,*) 'number of items read = ',ianz
	write(6,*) 'format should be: "k val" or "x y val"'
	stop 'error stop read_observations: wrong format'
   99	continue
	write(6,*) 'k = ',k
	stop 'error stop read_observations: no such node'
	end

c******************************************************************

	subroutine write_reg_grid(nx,ny,x0,y0,dx,dy)

	implicit none

	integer nx,ny
	real x0,y0
	real dx,dy

	integer iu
	integer ix,iy,i
	real x,y

	iu = 9
	open(iu,file='regfile.grd',form='formatted',status='unknown')

	i = 0
	do iy=0,ny
	  do ix=0,nx
	    i = i + 1
	    x = x0 + ix*dx
	    y = y0 + iy*dy
	    write(iu,*) 1,i,0,x,y
	  end do
	end do

	close(iu)

	end

c******************************************************************

	subroutine write_box_grid(nx,ny,x0,y0,dx,dy)

	implicit none

	integer nx,ny
	real x0,y0
	real dx,dy

	integer iu
	integer ix,iy,in,il
	real x,y
	real x1,x2,y1,y2

	iu = 9
	open(iu,file='boxfile.grd',form='formatted',status='unknown')

	in = 0
	il = 0

	do iy=0,ny
	    il = il + 1
	    x1 = x0
	    x2 = x0 + nx*dx
	    y = y0 + iy*dy
	    write(iu,*) 1,in+1,0,x1,y
	    write(iu,*) 1,in+2,0,x2,y
	    write(iu,*) 3,il,0,2,in+1,in+2
	    in = in + 2
	end do
	do ix=0,nx
	    il = il + 1
	    y1 = y0
	    y2 = y0 + ny*dy
	    x = x0 + ix*dx
	    write(iu,*) 1,in+1,0,x,y1
	    write(iu,*) 1,in+2,0,x,y2
	    write(iu,*) 3,il,0,2,in+1,in+2
	    in = in + 2
	end do

	close(iu)

	end

c******************************************************************

	subroutine write_obs_grid(nobs,xobs,yobs,zobs)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)

	integer iu
	integer i
	real x,y,z

	iu = 9
	open(iu,file='obsfile.grd',form='formatted',status='unknown')

	do i=1,nobs
	    x = xobs(i)
	    y = yobs(i)
	    z = zobs(i)
	    write(iu,*) 1,i,3,x,y,z
	end do

	close(iu)

	end

c******************************************************************

	subroutine compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar
     +					,za,iexcl,rmse)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	integer nx,ny
	real regpar(7)
	real za(nx,ny)
	integer iexcl		!excluded obs - compute error for this
	real rmse

	integer i,ix,iy,ntot
	integer istart,iend
	real dx,dy,x0,y0
	real x,y,z
	real rx,ry,xx,yy,zm
	real zz(4)
	real flag
	double precision rrr

	real rbilin

	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	rrr = 0.
	istart = 1
	iend = nobs
	ntot = nobs
	if( iexcl > 0 ) then
	  istart = iexcl
	  iend = iexcl
	  ntot = 1
	end if

	do i=istart,iend
	  x = xobs(i)
	  y = yobs(i)
	  z = zobs(i)
	  ix = 1 + (x-x0)/dx
	  if( ix >= nx ) ix = nx - 1
	  xx = x0 + (ix-1)*dx
	  rx = (x - xx)/dx
	  iy = 1 + (y-y0)/dy
	  if( iy >= ny ) iy = ny - 1
	  yy = y0 + (iy-1)*dy
	  ry = (y - yy)/dy
	  zz = (/za(ix,iy),za(ix+1,iy),za(ix,iy+1),za(ix+1,iy+1)/)
	  zm = rbilin(zz,rx,ry,flag)
	!write(6,*) i,ix,iy,rx,ry
	!write(6,*) zz,zm,z
	  if( iexcl > 0 ) write(67,*) i,z,zm
	  rrr = rrr + (zm-z)**2
	end do

	rmse = sqrt( rrr / ntot )

	end

c******************************************************************

        subroutine opt_intp_loop(nobs,xobs,yobs,zobs,bobs
     +                  ,nback,bback,xback,yback,zback,regpar
     +                  ,rlact,rlmax,ssa,rra,zanal,rmse,errors)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	real bobs(nobs)
        integer nback           !size of background field
        logical bback           !background values are given
        real xback(nback)       !x-coordinates of background
        real yback(nback)       !y-coordinates of background
        real zback(nback)       !values of background
	real regpar(7)
        real rlact                 !length scale for covariance
        real rlmax              !max radius to be considered
        real ssa(nobs)              !std of background field
        real rra(nobs)                 !std of observation error matrix
        real zanal(nback)       !analysis on return
	real rmse
	real errors(nobs)

	integer i
	integer nx,ny
	real error,rrsave
	double precision rrr

	rrr = 0.
	nx =  nint(regpar(1))
	ny =  nint(regpar(2))

	do i=1,nobs
	  rrsave = rra(i)
	  rra(i) = 1000. * rra(i)
          call opt_intp(nobs,xobs,yobs,zobs,bobs
     +                  ,nback,bback,xback,yback,zback
     +                  ,rlact,rlmax,ssa,rra,zanal)
	  call compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar
     +					,zanal,i,error)
	  errors(i) = error
	  rrr = rrr + error**2
	  rra(i) = rrsave
	end do

	rmse = sqrt( rrr / nobs )

	end

c******************************************************************

	subroutine read_use_file(nobs,iuse)

	implicit none

	integer nobs
	integer iuse(nobs)

	integer ios,i,j,iu,n

	iuse = 1

	iu = 1
	open(iu,file='use.txt',status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) return

	read(iu,*) n
	if( n /= nobs ) then
	  write(6,*) n,nobs
	  stop 'error stop read_use_file: nobs'
	end if

	do i=1,nobs
	  read(iu,*) j,iuse(i)
	  if( i /= j ) stop 'error stop read_use_file: i/=j'
	end do

	close(iu)

	end

c******************************************************************

