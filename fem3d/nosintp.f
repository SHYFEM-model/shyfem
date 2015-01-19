c
c $Id: nosintp.f,v 1.1 1999/01/25 12:46:07 georg Exp $
c
c 20.11.1998    ggu     written and adapted from nosextr.f
c
c********************************************************

	program nosintp

c interpolates nos file onto regular grid

	include 'param.h'

c arrays from basin


	include 'basin.h'

c arrays for nos

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

c conversion to geographical coordinates

	real phi0,omega0
	parameter(phi0=41.0)		! reference latitude
	parameter(omega0=16.0)		! reference longitude

c interpolation matrix

	real lon0,lat0,lon1,lat1,dlon,dlat
	parameter (lon0=-3839)		! lower left point of matrix
	parameter (lat0=17158)		! lower left point of matrix
	parameter (lon1=1512)		! upper right point of matrix
	parameter (lat1=18717)		! upper right point of matrix
	parameter (dlon=20)		! grid spacing (longitude)
	parameter (dlat=20)		! grid spacing (latitude)

	real flag
	parameter ( flag = -1 )

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
	rnull=0.

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

!	call mercin(phi0,omega0,100.)
!	call c2g(nkn,xgv,ygv,xgv,ygv)

c	do k=1,nkn
c	  write(6,*) k,ipext(k),xgv(k),ygv(k)
c	end do

c read first and second header ----------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

c interpolate depth and find maximum nlv ------------------------------

	call ave2am(hev,hm,nx,ny)
	call mklev(hm,levels,nx,ny,hlv,nlv)
	call mimareg(hm,nx,ny,hmin,hmax)

	bdepth = .true.
	nlvmax = 1
	do l=1,nlv
	  if( bdepth ) nlvmax = l
	  if( hlv(l) .ge. hmax ) bdepth = .false.
	end do

	write(6,*) 'min/max depth   : ',hmin,hmax
	write(6,*) 'Available levels: ',nlv,nlvmax
	write(6,*) (hlv(l),l=1,nlv)

c time loop ------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar
	write(83,*)it,k
	do l=1,nlvmax
	  do k=1,nkn
	    cv(k)=cv3(l,k)
	    write(83,*)xgv(k),ygv(k),cv(k)

	  end do
	  call av2am(cv,cm,nx,ny)
	  call adjam(cm,levels,nx,ny,l,flag,nused)
	  call mimareg(cm,nx,ny,cmin,cmax)
	  write(6,*) 'level,cmin,cmax : ',l,cmin,cmax,nused
c	  call a2char(cm,cc,nx,ny)
c	  call prchar(cc,nx,ny)
	  call writec(cm,nx,ny,it,l,ivar,dlat,dlon)
	end do

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

	subroutine writec(cm,ip,jp,it,l,ivar,dlat,dlon)

c writes interpolated value to file
c
c cm(1,1)   is lower,left corner of domain
c cm(ip,jp) is upper,right corner of domain
c cm(ix,iy) is interior point, (ix,iy) are indices of x/y coordinate

	implicit none

	integer ip,jp		!dimensions of matrix
	real cm(ip,jp)		!matrix with interpolated values
	integer it		!time of simulation (seconds)
	integer l		!vertical level (l=1 -> surface)
	integer ivar		!type of variable
	real dlat,dlon
	integer i,j,ii,jj
	
	real area 

	area =0


!	write(99,'(a,5i12)') '*',it,l,ivar,ip,jp
!	write(99,*) ((cm(i,j),i=1,ip),j=1,jp)

	do i=1,ip
         do j=1,jp
          if(cm(i,j).ge.10)area=area+dlat*dlon
         end do
        end do
         
        write(89,*),it,area     



	end

c***************************************************************

