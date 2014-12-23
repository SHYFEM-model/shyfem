c
c $Id: ousintp.f,v 1.5 2010-01-18 11:45:56 georg Exp $
c
c interpolation of velocities onto nodes
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 04.03.2005	ggu	computes 3D velocities
c 24.11.2009	ggu	commented and cleaned
c 30.09.2010	ggu	cleaned, new vars itmin/itmax to limit writing
c 16.12.2010	ggu	transp2vel moved to ousutil.f
c
c***************************************************************

	program ousintp

c reads ous file and interpolates velocities
c
c we would not even need to read basin

	implicit none

        include 'param.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	include 'nbasin.h'

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)
        common /uprv/uprv
        common /vprv/vprv

	real hev(neldim)
	real weight(nlvdim,nkndim)

	real znv(nkndim)
	real zenv(3,neldim)

	real uvel(nkndim)
	real vvel(nkndim)
	real hl(nlvdim)

	integer ndim
	parameter(ndim=100)
	real xpn(ndim), ypn(ndim)
	integer ielv(ndim)

	logical btime
	integer n,nx,ny

        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i,k
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real cmin,cmax

c	integer rdous,rfous
	integer iapini,ideffi

c adriatic sea
        real lon0,lat0,lon1,lat1,dlon,dlat
        parameter (lon0=12.0)           ! lower left point of matrix
        parameter (lat0=40.0)           ! lower left point of matrix
        parameter (lon1=20.0)           ! lower left point of matrix
        parameter (lat1=46.0)           ! lower left point of matrix
        parameter (dlon=1./12.)         ! grid spacing (longitude)
        parameter (dlat=1./12.)         ! grid spacing (latitude)

	real xmin,ymin,xmax,ymax,dx,dy
        !parameter (xmin=5000.0)           ! lower left point of matrix
        !parameter (ymin=1000.0)           ! lower left point of matrix
        !parameter (xmax=7000.0)           ! lower left point of matrix
        !parameter (ymax=3000.0)           ! lower left point of matrix
        !parameter (xmin=14.03)           ! lower left point of matrix
        !parameter (ymin=35.65)           ! lower left point of matrix
        !parameter (xmax=14.73)           ! lower left point of matrix
        !parameter (ymax=36.23)           ! lower left point of matrix
        !parameter (dx=.01)         ! grid spacing (longitude)
        !parameter (dy=.01)         ! grid spacing (latitude)
        parameter (xmin=12000.0)           ! lower left point of matrix
        parameter (ymin=2000.0)           ! lower left point of matrix
        parameter (xmax=54000.0)           ! lower left point of matrix
        parameter (ymax=54000.0)           ! lower left point of matrix
        parameter (dx=1000.)         ! grid spacing (longitude)
        parameter (dy=1000.)         ! grid spacing (latitude)

        logical blatlon
        parameter ( blatlon = .false. )
        logical bwritevel,bwritezeta,bwritegeom
        parameter (bwritevel=.true.,bwritezeta=.true.,bwritegeom=.true.)
	integer itmin,itmax
	!parameter (itmin=0,itmax=0)
	parameter (itmin=30672000,itmax=30931200)

c (19941222)  1994 12 22    30672000
c (19941225)  1994 12 25    30931200

        real flag
        parameter ( flag = -999.0 )

        integer nxdim,nydim
        !parameter ( nxdim = (lon1-lon0)/dlon + 1.5 )
        !parameter ( nydim = (lat1-lat0)/dlat + 1.5 )
        parameter ( nxdim = (xmax-xmin)/dx + 1.5 )
        parameter ( nydim = (ymax-ymin)/dy + 1.5 )

        real cm(nxdim,nydim)
        real hm(nxdim,nydim)            ! interpolated depth
        character*1 cc(nxdim,nydim)
        integer levels(nxdim,nydim)     ! maximum level of point

c---------------------------------------------------------------------
c set up parameters
c---------------------------------------------------------------------

	nread=0
	nx = nxdim
	ny = nydim

        write(6,*) 'dimension of matrix : ',nx,ny

        call setgeo(xmin,ymin,dx,dy,flag)
        if( blatlon ) call setgeo(lon0,lat0,dlon,dlat,flag)

c---------------------------------------------------------------------
c read basin and maybe convert to lat/lon
c---------------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        if( blatlon ) call c2g(nkn,xgv,ygv,xgv,ygv)

	!call read_points('nome_file.dat',ndim,n,xpn,ypn)
	!call get_elements(n,xpn,ypn,ielv)

	if( bwritegeom ) then
	  it = 0
          call av2am(xgv,cm,nx,ny)
          call writec('__x ',it,flag,nx,ny,cm)
          call av2am(ygv,cm,nx,ny)
          call writec('__y ',it,flag,nx,ny,cm)
	end if

c---------------------------------------------------------------------
c prepare file name and read header of OUS file
c---------------------------------------------------------------------

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	call rsous(nin,ilhv,hlv,hev,ierr)

c---------------------------------------------------------------------
c start time loop over simulation results
c---------------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) then
		write(6,*) 'error in reading file : ',ierr
		goto 100
        else if(ierr.lt.0) then
		goto 100
	end if

	nread=nread+1

	btime = itmin.ge.itmax .or. ( it.ge.itmin .and. it.le.itmax )

c	---------------------------------------------
c handle water levels
c	---------------------------------------------

	call mima(znv,nknous,zmin,zmax)
	if( btime .and. bwritezeta ) then
          call av2am(znv,cm,nx,ny)
          call writec('__z ',it,flag,nx,ny,cm)
          call writefem('__z ',it,znv)
	end if

c	---------------------------------------------
c handle velocities
c	---------------------------------------------

        call transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight,hl)

	do k=1,nkn
	  uvel(k) = uprv(1,k)
	  vvel(k) = vprv(1,k)
	end do

	if( btime .and. bwritevel ) then
          call av2am(uvel,cm,nx,ny)
          call writec('__u ',it,flag,nx,ny,cm)
          call av2am(vvel,cm,nx,ny)
          call writec('__v ',it,flag,nx,ny,cm)
	end if

c	---------------------------------------------
c write message to terminal
c	---------------------------------------------

	write(6,*) 'time : ',it,'  zmin/zmax : ',zmin,zmax, ' write : ',btime

	goto 300

c---------------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'regular output has been written to unit 99 (fort.99)'
	write(6,*) 'fem     output has been written to unit 98 (fort.98)'
	write(6,*)

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c******************************************************************

	subroutine read_points(file,ndim,n,xpn,ypn)
	character*(*) file
	end

c******************************************************************

	subroutine get_elements(n,xpn,ypn,ielv)

	implicit none

	integer n
	real xpn(1), ypn(1)
	integer ielv(1)

	integer i,ielem

	do i=1,n
	  call find_element(xpn(i),ypn(i),ielem)
	  if( ielem .le. 0 ) then
	    write(6,*) i,xpn(i),ypn(i)
	    stop 'error stop get_elements: no such element'
	  end if
	  ielv(i) = ielem
	end do

	end

c******************************************************************

	subroutine interp_vel(it,n,xpn,ypn,ielv,nlvdim,nen3v,uprv,vprv)

	implicit none

	integer it
	integer n
	real xpn(1), ypn(1)
	integer ielv(1)
	integer nlvdim
	integer nen3v(3,1)
	real uprv(nlvdim,1)
	real vprv(nlvdim,1)

	integer i,ii,ie,k,level
	real xp,yp
	real up,vp
	real x(3), y(3)
	real u(3), v(3)

	level = 1	!only surface

	write(66,*) it,n

	do i=1,n
	  ie = ielv(i)
	  xp = xpn(i)
	  yp = ypn(i)
	  call get_xy_elem(ie,x,y)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    u(ii) = uprv(level,k)
	    v(ii) = vprv(level,k)
	  end do
	  call elemintp(x,y,u,xp,yp,up)
	  call elemintp(x,y,v,xp,yp,vp)
	  write(66,*) i,up,vp
	end do

	end

c******************************************************************

        subroutine c2g(nkn,xin,yin,xout,yout)

	implicit none

	integer nkn
	real xin(nkn), yin(nkn)
	real xout(nkn), yout(nkn)

	integer k
	real x,y

	do k=1,nkn
	  x = xin(k)
	  y = yin(k)

	  xout(k) = x
	  yout(k) = y
	end do

	end

c******************************************************************

        subroutine writec(text,it,flag,ip,jp,cm)

c writes interpolated value to file
c
c cm(1,1)   is lower,left corner of domain
c cm(ip,jp) is upper,right corner of domain
c cm(ix,iy) is interior point, (ix,iy) are indices of x/y coordinate

        implicit none

	character*4 text
        integer it              !time of simulation (seconds)
        real flag               !flag for no value (land)
        integer ip,jp           !dimensions of matrix
        real cm(ip,jp)          !matrix with interpolated values

        integer i,j

        write(99,'(a4,3i12,f14.6)') text,it,ip,jp,flag

        !write(99,*) ((cm(i,j),i=1,ip),j=1,jp)

	do j=1,jp
	  do i=1,ip
            write(99,*) i,j,cm(i,j)
	  end do
	end do

        end

c******************************************************************

	subroutine writefem(text,it,value)

c writes value of fem grid to file

	implicit none

        include 'param.h'

	character*4 text
        integer it              !time of simulation (seconds)
        real value(1)           !values of fem grid

	include 'nbasin.h'

	real xgv(nkndim), ygv(nkndim)
	common /xgv/xgv, /ygv/ygv

	integer k
	real x,y,v

        write(98,'(a4,3i12,f14.6)') text,it,nkn

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  v = value(k)
	  write(98,*) x,y,v
	end do

        end

c******************************************************************

