c
c $Id: laplap.f,v 1.9 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 20.08.2003	ggu	new laplacian interpolation
c 02.09.2003	ggu	some comments, write to .dat file
c 30.10.2003	ggu	subroutine prepare_bc_l included in this file
c 04.03.2004	ggu	writes also number of variables (1)
c 11.03.2009	ggu	bug fix -> declare hev() here
c
c notes :
c
c please prepare file like this:
c
c----------------- start
c
c k1	val1
c k2	val2
c ...
c kn	valn
c----------------- end
c
c first line of file must be empty !!!
c
c run memory and set the basin ( memory -b venlag62 )
c run laplap with input file ( laplap < input.dat )
c
c****************************************************************

        program optintp

c laplacian interpolation

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'evmain.h'

	integer matdim
	parameter (matdim = nkndim*100)
	integer nobdim
	parameter (nobdim = 100)

        real hev(neldim)
        common /hev/hev

	real zv(nkndim)
	integer node(nobdim)
	integer ivec(nobdim)
	real obs(nobdim)
	real rvec(nobdim)
	real raux(nobdim)
	real rmat(nobdim,nobdim)

	integer k,ie,n
        integer ilev,ivar
	real flag
	real zmin,zmax

	real rl,rlmax,sigma,rr

	integer iapini

	flag = 1.23456e+23

	rl = 1000.		! length scale for covariance matrix
	rlmax = 10000.		! max radius to be considered
	sigma = 0.3		! std of background field
	rr = 0.02		! std of observation errors

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call bas_info

c-----------------------------------------------------------------
c set up ev
c-----------------------------------------------------------------

	call set_ev
	call check_ev

c-----------------------------------------------------------------
c read observations and interpolate
c-----------------------------------------------------------------

	call read_observations(' ',n,node,obs)

	call optintp_field(n,nkn,node,obs,rl,rlmax,sigma,rr,zv
     +				,ivec,rvec,raux,rmat)

c-----------------------------------------------------------------
c min/max of interpolated values
c-----------------------------------------------------------------

	call mima(zv,nkn,zmin,zmax)
	write(6,*) 'min/max: ',zmin,zmax

c-----------------------------------------------------------------
c write to NOS file
c-----------------------------------------------------------------

	call wrnos2d('optintp','optimal interpolation',zv)

c-----------------------------------------------------------------
c write to DAT file laplace.dat
c-----------------------------------------------------------------

        ilev = 0
        ivar = 1

	open(1,file='optintp.dat',status='unknown',form='unformatted')
	write(1) nkn,ilev,ivar
	write(1) (zv(k),k=1,nkn)
	close(1)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine optintp_field(n,nkn,node,obs,rl,rlmax,sigma,rr,zv
     +				,ivec,rvec,raux,rmat)

c computes optimal interpolation

	implicit none

	integer n		!size of observations
	integer nkn		!size of background field
	integer node(n)	!list of (internal) node numbers
	real obs(n)		!values of observations
	real rl			!length scale
	real rlmax		!max radius to be considered
	real sigma		!std of background field
	real rr			!std of error matrix
	real zv(nkn)		!analysis on return
	integer ivec(n)		!aux vector (n)
	real rvec(n)		!aux vector (n)
	real raux(n)		!aux vector (n)
	real rmat(n,n)		!aux matrix (nxn)

	include 'param.h'

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer i,j,ki,kj,k
	real rl2,rlmax2,rr2,rmean
	real xi,yi,xj,yj,dist2,r,z,acu

	rl2 = rl**2
	rlmax2 = rlmax**2
	rr2 = rr**2

c	------------------------------------------
c	compute mean of observations
c	------------------------------------------

	rmean = 0.
	do i=1,n
	  rmean = rmean + obs(i)
	end do
	rmean = rmean / n

c	------------------------------------------
c	subtract mean from observations and set background to 0 (mean)
c	------------------------------------------

	do i=1,n
	  obs(i) = obs(i) - rmean
	end do

	do k=1,nkn
	  zv(k) = 0.
	end do

c	------------------------------------------
c	set up covariance matrix H P^b H^T
c	------------------------------------------

	do j=1,n
	  kj = node(j)
	  xj = xgv(kj)
	  yj = ygv(kj)
	  do i=1,n
	    ki = node(i)
	    xi = xgv(ki)
	    yi = ygv(ki)
	    dist2 = (xi-xj)**2 + (yi-yj)**2
	    r = exp( -dist2/rl2 )
	    if( dist2 .gt. rlmax2 ) r = 0.
	    rmat(i,j) = r * sigma**2
	  end do
	end do

c	------------------------------------------
c	add observation error matrix
c	------------------------------------------

	do j=1,n
	  rmat(j,j) = rmat(j,j) + rr2
	end do

c	------------------------------------------
c	invert matrix
c	------------------------------------------

	call matinv(rmat,ivec,rvec,n,n)

c	------------------------------------------
c	create observational innovation vector
c	------------------------------------------

	do j=1,n
	  k = node(j)
	  z = zv(k)
	  rvec(j) = obs(j) - z
	end do

c	------------------------------------------
c	multiply inverted matrix with observational innovation
c	------------------------------------------

	do i=1,n
	  acu = 0.
	  do j=1,n
	    acu = acu + rmat(i,j) * rvec(j)
	  end do
	  raux(i) = acu
	end do

c	------------------------------------------
c	multiply P^b H^T with vector and add to background to obtain analysis
c	------------------------------------------

	do k=1,nkn
	  xi = xgv(k)
	  yi = ygv(k)
	  acu = 0.
	  do j=1,n
	    kj = node(j)
	    xj = xgv(kj)
	    yj = ygv(kj)
	    dist2 = (xi-xj)**2 + (yi-yj)**2
	    r = exp( -dist2/rl2 )
	    if( dist2 .gt. rlmax2 ) r = 0.
	    acu = acu + r * sigma**2 * raux(j)
	  end do
	  zv(k) = zv(k) + rmean + acu
	end do

c	------------------------------------------
c	end of routine
c	------------------------------------------

	end

c****************************************************************

	subroutine read_observations(file,n,node,obs)

c reads boundary conditions from file and sets up array
c
c file must be made like this:
c
c	k1, val1
c	k2, val2
c	...
c	kn, valn

	implicit none

	character*(*) file
	integer n
	integer node(n)
	real obs(n)

	integer iunit
	integer k,kn,ndim
	real val

	integer ipint

	ndim = n

	if( file .ne. ' ' ) then
	  open(1,file=file,status='old',form='formatted',err=97)
	  iunit = 1
	else
	  iunit = 5
	end if

	write(6,*) '...reading boundary conditions from unit :',iunit
	write(6,*) '   format: k  val'

	n = 0
    1	continue
	  read(iunit,*,end=2) k,val
	  n = n + 1
	  if( n .gt. ndim ) goto 96
	  kn = ipint(k)
	  if( kn .le. 0 ) goto 99
	  node(n) = kn
	  obs(n) = val
	  write(6,*) n,k,kn,val
	  goto 1
    2	continue

	if( iunit .ne. 5 ) close(iunit)

	return
   96	continue
	write(6,*) n,ndim
	stop 'error stop read_observations: dimensions'
   97	continue
	write(6,*) file
	stop 'error stop read_observations: cannot open file'
   98	continue
	write(6,*) k,kn,val
	stop 'error stop read_observations: error in internal node number'
   99	continue
	write(6,*) k,kn,val
	stop 'error stop read_observations: no such node'
	end

c******************************************************************

