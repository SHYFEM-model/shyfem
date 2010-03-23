c
c $Id: outintp.f,v 1.1 1999/12/09 17:28:42 georg Exp $
c
c 20.11.1998    ggu     written and adapted from nosextr.f
c
c********************************************************

	program outintp

c interpolates out file onto regular grid (only zeta)

	implicit none

	include 'param.h'

c arrays from basin

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        character*80 descrr
        common /descrr/ descrr

	integer nen3v(3,neldim)
	integer ipv(nkndim), ipev(neldim)
	integer iarv(neldim)
	real hm3v(3,neldim)
	real xgv(nkndim), ygv(nkndim)
        common /nen3v/nen3v
        common /iarv/iarv
        common /ipv/ipv, /ipev/ipev
        common /hm3v/hm3v
        common /xgv/xgv, /ygv/ygv

c arrays for out

	character*80 title

        real xv(3,nkndim)
        real zenv(3,neldim)
        real unv(neldim),vnv(neldim)

        real znv(nkndim)
        real uv(nkndim)
        real vv(nkndim)

        integer nvers,nin
        integer itanf,itend,idt,idtout
        integer it,k
	integer day
        integer ierr,nread
        integer nkk,nee
	integer nx,ny
	real rnull
        real href,hzoff
	real hmin,hmax
	real cmin,cmax
        real zmin,zmax
        real umin,umax
        real vmin,vmax
        real uvmin,uvmax
        real vvmin,vvmax

        integer rdout,rfout
        integer iapini,ideffi

	real hev(neldim)

c conversion to geographical coordinates

	real phi0,omega0
c adriatic sea
c	parameter(phi0=41.0)		! reference latitude
c	parameter(omega0=16.0)		! reference longitude
c mediterranean
	parameter(phi0=41.5)		! reference latitude
	parameter(omega0=0.0)		! reference longitude

c interpolation matrix

	real lon0,lat0,lon1,lat1,dlon,dlat
c adriatic sea
c	parameter (lon0=12.0)		! lower left point of matrix
c	parameter (lat0=43.5)		! lower left point of matrix
c	parameter (lon1=14.0)		! upper right point of matrix
c	parameter (lat1=46.0)		! upper right point of matrix
c	parameter (dlon=2./60.)		! grid spacing (longitude)
c	parameter (dlat=2./60.)		! grid spacing (latitude)
c mediterranean
	parameter (lon0=-6.0)		! lower left point of matrix
	parameter (lat0=30.0)		! lower left point of matrix
	parameter (lon1=37.0)		! upper right point of matrix
	parameter (lat1=46.0)		! upper right point of matrix
c	parameter (dlon=1.)		! grid spacing (longitude)
c	parameter (dlat=1.)		! grid spacing (latitude)
	parameter (dlon=.25)		! grid spacing (longitude)
	parameter (dlat=.25)		! grid spacing (latitude)

	real flag
	parameter ( flag = -999.0 )

	integer nxdim,nydim
	parameter ( nxdim = (lon1-lon0)/dlon + 1.5 )
	parameter ( nydim = (lat1-lat0)/dlat + 1.5 )

	real cm(nxdim,nydim)
	real hm(nxdim,nydim)		! interpolated depth
	character*1 cc(nxdim,nydim)
	integer levels(nxdim,nydim)	! maximum level of point

	logical bdepth

c----------------------------------------------------------------------

	nread=0

c set up of interpolation ---------------------------------------------

	nx = nxdim
	ny = nydim

	write(6,*) 'dimension of matrix : ',nx,ny

	call setgeo(lon0,lat0,dlon,dlat,flag)

c read in of basin and simulation name --------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c convert to geographical coordinates ---------------------------------

	call mercin(phi0,omega0,100.)
	call c2g(nkn,xgv,ygv,xgv,ygv)

c	do k=1,nkn
c	  write(6,*) k,ipext(k),xgv(k),ygv(k)
c	end do

c read first and second header ----------------------------------------

	nin=ideffi('datdir','runnam','.out','unform','old')
	if(nin.le.0) goto 100

        nvers=0
        ierr=rfout(nin,nvers,itanf,itend,idt,idtout,href,hzoff,title)
        if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*)   title
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' itanf,itend  : ',itanf,itend
        write(6,*) ' idt,idtout   : ',idt,idtout
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nkn,nel
        write(6,*)

c interpolate depth and find maximum nlv ------------------------------

	call ave2am(hev,hm,nx,ny)
	call mimareg(hm,nx,ny,hmin,hmax)

c time loop ------------------------------------------------------------

  300   continue

        nkk=nkn
        nee=nel
        ierr=rdout(nin,nvers,it,nkk,nee,xv,zenv,unv,vnv)

        if(ierr.gt.0) then
                write(6,*) 'error in reading file : ',ierr
                goto 100
        else if(ierr.lt.0) then
                goto 100
        else if(nkn.ne.nkk.or.nel.ne.nee) then
                write(6,*) 'Too less data read'
                write(6,*) 'nkn,nkk :',nkn,nkk
                write(6,*) 'nel,nee :',nel,nee
                goto 100
        end if

	nread=nread+1
	day = it/86400

        do k=1,nkn
          znv(k) = xv(3,k)
        end do

        call mima(znv,nkn,zmin,zmax)

        write(6,*)
        write(6,*) 'time : ',it,day
        write(6,*)
        write(6,*) ' zmin/zmax  : ',zmin,zmax

	call av2am(znv,cm,nx,ny)
	call mimareg(cm,nx,ny,cmin,cmax)
	write(6,*) 'cmin,cmax : ',cmin,cmax
	call writec(cm,nx,ny,it,flag)

	goto 300

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

	end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

	subroutine mklev(hm,levels,ip,jp,hlv,nlv)

c computes maximum levels for each interpolation point

	implicit none

	integer ip,jp
	integer nlv
	real hm(ip,jp)
	integer levels(ip,jp)
	real hlv(1)

	integer i,j,l,lmax
	real h

	do j=1,jp
	  do i=1,ip
	    h = hm(i,j)
            lmax = 1
            do l=1,nlv
              lmax = l
              if( hlv(l) .ge. h ) goto 1
            end do
    1	    continue
	    levels(i,j) = lmax
	  end do
	end do

	end

c***************************************************************

	subroutine adjam(am,levels,ip,jp,lact,flag,nused)

c adjusts deeper layers for land

	implicit none

	integer ip,jp
	integer lact
	real am(ip,jp)
	integer levels(ip,jp)
	real flag
	integer nused

	integer i,j

	nused = 0

	do j=1,jp
	  do i=1,ip
	    if( lact .gt. levels(i,j) ) am(i,j) = flag
	    if( am(i,j) .ne. flag ) nused = nused + 1
	  end do
	end do

	end

c***************************************************************

	subroutine writec(cm,ip,jp,it,flag)

c writes interpolated value to file
c
c cm(1,1)   is lower,left corner of domain
c cm(ip,jp) is upper,right corner of domain
c cm(ix,iy) is interior point, (ix,iy) are indices of x/y coordinate

	implicit none

	integer ip,jp		!dimensions of matrix
	real cm(ip,jp)		!matrix with interpolated values
	integer it		!time of simulation (seconds)
	real flag		!flag for no value (land)

	integer i,j

	write(99,'(a,3i12,e14.6)') '*',it,ip,jp,flag
	write(99,*) ((cm(i,j),i=1,ip),j=1,jp)

	end

c***************************************************************

