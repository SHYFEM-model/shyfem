c
c revision log :
c
c 20.08.2003	ggu	new laplacian interpolation
c 02.09.2003	ggu	some comments, write to .dat file
c 30.10.2003	ggu	subroutine prepare_bc_l included in this file
c 04.03.2004	ggu	writes also number of variables (1)
c 11.03.2009	ggu	bug fix -> declare hev() here
c 07.02.2015	ggu	OI finished
c 14.09.2015	ggu	OI adapted to new modular structure
c 29.09.2015	ggu	move some vars to custimize to the beginning of file
c
c notes :
c
c please prepare file like this:
c
c----------------- start
c k1	val1
c k2	val2
c ...
c kn	valn
c----------------- end
c
c or alternatively like this:
c
c----------------- start
c x1 y1    val1
c x2 y2    val2
c ...
c xn yn    valn
c----------------- end
c
c****************************************************************

        program optintp

c optimal interpolation interpolation

	use mod_depth
	use evgeom
	use basin
	use clo

	implicit none

	include 'param.h'

	integer nobdim
	parameter (nobdim = 1000)


	real xobs(nobdim)
	real yobs(nobdim)
	real zobs(nobdim)
	real bobs(nobdim)

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
	double precision dtime
	real hlv(1)
	real regpar(7)
	integer date
	integer datetime(2)
	character*30 string,format

	character*80 file,basnam
	logical bback,bgeo,bcart,bquiet,blimit,breg,bnos,bmulti
	integer k,ie,n,ndim
	integer nobs,nback
        integer ilev,ivar,irec
	integer it,index
	real zmin,zmax
	real zomin,zomax
	real xmin,ymin,xmax,ymax
	real dx,dy,dxy
	real drl,rlact

	real rl,rlmax,rr,ss

	integer iapini
	logical is_spherical

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

        call clo_init('optintp','input-file','1.0')

        call clo_add_info('Optimal interpolation for sparse data')
        call clo_add_option('basin',' ','basin to be used')
        call clo_add_option('rl #',-1.
     +		,'set length scale for covariance matrix')
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
        call clo_add_option('multi',.false.
     +		,'read multiple data with time')
        call clo_add_option('cart',.false.
     +		,'force use of cartesian coordinates')
        call clo_add_option('limit',.false.
     +		,'limit values to min/max of observations')
        call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_extra('defaults for undefined values:')
	call clo_add_extra('rl: 1/10 of basin size')
	call clo_add_extra('rlmax=10*rl  rr=0.01  ss=100*rr')
	call clo_add_extra('Default for dx is 0 (use FEM grid)')

        call clo_parse_options(1)  !expecting (at least) 1 file after options

        call clo_get_option('basin',basnam)
        call clo_get_option('cart',bcart)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('rl',rl)
        call clo_get_option('rlmax',rlmax)
        call clo_get_option('rr',rr)
        call clo_get_option('ss',ss)
        call clo_get_option('limit',blimit)
        call clo_get_option('multi',bmulti)
        call clo_get_option('dx',dx)
        call clo_get_option('dy',dy)

	!call clo_info()

	bback = .false.
	call clo_get_file(1,file)

c-----------------------------------------------------------------
c set some variables
c-----------------------------------------------------------------

	ivar = 85		!what is in the observations - ice
	string = 'ice cover [0-1]'
	date = 19970101

	ivar = 12		!what is in the observations - temperature
	string = 'temperature [C]'
	date = 20100101

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

	if( .not. bquiet ) then
	  write(6,*) 'parameters used:'
	  write(6,*) 'xmin/xmax: ',xmin,xmax
	  write(6,*) 'ymin/ymax: ',ymin,ymax
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
	call setup_background(ndim,dx,dy,xmin,ymin,xmax,ymax
     +			,nback,xback,yback,regpar)

	breg = dx > 0.
	bnos = .not. breg

	if( .not. bquiet ) then
	  write(6,*) 'background grid:'
	  write(6,*) 'nback: ',nback
	  if( breg ) then
	    write(6,*) 'nx,ny: ',nint(regpar(1)),nint(regpar(2))
	    write(6,*) 'dx,dy: ',dx,dy
	    write(6,*) 'x0,y0: ',regpar(3),regpar(4)
	  end if
	end if

	write(6,*) 'ggu 1'
	iformat = 0
	format = 'unformatted'
	ntype = 1
	if( breg ) then
	  iformat = 1
	  format = 'formatted'
	  ntype = ntype + 10
	end if

	write(6,*) 'ggu 2'
	nvers = 0
	nvar = 1
	lmax = 1
	nlvdi = 1
	np = nback
	write(6,*) 'ggu 2'
	it = 0
	datetime = 0
	datetime(1) = date
	write(6,*) 'ggu 2'
	hlv(1) = 10000.
	write(6,*) 'ggu 2'
	hd = 1.
	ilhkv = 1

	write(6,*) 'ggu 3'

c-----------------------------------------------------------------
c open files
c-----------------------------------------------------------------

	irec = 0
	drl = 0.05		!test different rl, distance drl
	drl = 0.
	if( bmulti ) drl = 0.
	jmax = 0
	if( drl > 0. ) jmax = 5
	if( jmax > 0 ) then
	  write(6,*) 'trying multiple values for rl: ',drl,jmax
	end if

	if( bnos ) then
	  call wrnos2d_open(iunos,'optintp','optimal interpolation')
	end if

	iufem = 2
	open(iufem,file='optintp.fem',status='unknown',form=format)

	iuobs = 1
	open(iuobs,file=file,status='old',form='formatted')
	write(6,*) 'file with observations: ',trim(file)

c-----------------------------------------------------------------
c read observations and interpolate
c-----------------------------------------------------------------

	do

	nobs = nobdim
	call read_observations(bmulti,iuobs
     +				,dtime,nobs,xobs,yobs,zobs)

	if( nobs <= 0 ) exit

	irec = irec + 1
	call mima(zobs,nobs,zomin,zomax)
	write(6,*) 'observations min/max: ',zomin,zomax

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	do j=-jmax,jmax
	  rlact = rl + j*drl
          call opt_intp(nobs,xobs,yobs,zobs,bobs
     +                  ,nback,bback,xback,yback,zback
     +                  ,rlact,rlmax,ss,rr,zanal)

	  call mima(zanal,nback,zmin,zmax)
	  if( bmulti ) then
	    write(6,*) 'min/max: ',dtime,zmin,zmax
	  else 
	    write(6,*) 'min/max: ',rlact,zmin,zmax
	  end if

	  if( blimit ) call limit_values(nback,zanal,zomin,zomax)

	  if( bnos ) call wrnos2d_record(iunos,it,ivar,zanal)

	  call fem_file_write_header(iformat,iufem,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvdi,hlv,datetime,regpar)
          call fem_file_write_data(iformat,iufem
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,zanal)

	  it = it + 1
	  dtime = dtime + 1
	end do

	if( .not. bmulti ) exit

	end do

	write(6,*) 'total number of observations treated: ',irec

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine setup_background(ndim,dx,dy,xmin,ymin,xmax,ymax
     +			,nback,xback,yback,regpar)

	use basin

	implicit none

	include 'param.h'

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

	subroutine read_observations(bmulti,iuobs
     +				,dtime,nobs,xobs,yobs,zobs)

c reads boundary conditions from file and sets up array
c
c file must be made like this:
c
c	k1, val1
c	k2, val2
c	...
c	kn, valn

	use basin

	implicit none

	include 'param.h'

	logical bmulti
	integer iuobs
	double precision dtime
	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)

	character*80 line
	logical bdebug
	integer k,kn,ndim,n,ianz
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

	if( bmulti ) read(iuobs,*,end=2) dtime,nobs

	do
	  if( nobs > 0 .and. n == nobs ) exit	!for multi data
	  read(iuobs,'(a)',end=2) line
	  !write(6,*) trim(line)
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
    2	continue

	nobs = n

	return
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

