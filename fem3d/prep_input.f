c
c $Id: prep_input.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c prepares regular data for FEM model (initial conditions)
c
c revision log :
c
c 03.10.2012    ggu     written from basinf
c
c****************************************************************

        program prep_input

c prepares input data from regular data

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'evmain.h'

c-----------------------------------------------------------------
c start user defined data
c-----------------------------------------------------------------

	integer nxdim,nydim,nzdim
	parameter (nxdim=100,nydim=100,nzdim=30)
	real level(nxdim,nydim)
	real uvel(nzdim,nxdim,nydim)
	real vvel(nzdim,nxdim,nydim)
	real temp(nzdim,nxdim,nydim)
	real salt(nzdim,nxdim,nydim)
	integer nx,ny,nz

c-----------------------------------------------------------------
c end user defined data
c-----------------------------------------------------------------

	real hkv(nkndim)
	real hev(neldim)
	real haux(nkndim)

	integer nlv			!total number of levels (return)
	real hlv(nlvdim)		!vertical levels (return)
	integer ilhv(neldim)		!number of levels in element (return)
	integer ilhkv(nkndim)		!number of levels in node (return)

	real val2d(nkndim)
	real val3d(nlvdim,nkndim)

	logical bformat
	integer it,nvar,iunit,lmax
	integer nvers,ntype
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call set_ev
	call makehev(hev)	!makes elementwise depth
	call makehkv(hkv,haux)	!makes nodewise depth

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcor,'  dirn = ',dirn
        write(6,*)

c-----------------------------------------------------------------
c set vertical structure
c-----------------------------------------------------------------

	call set_vertical(hev,nlv,hlv,ilhv,ilhkv)

c-----------------------------------------------------------------
c interpolate data and write files
c
c the routines get_data and intp_data2d, intp_data3d must be
c provided and customized
c-----------------------------------------------------------------

	nvers = 0
	ntype = 0
	lmax = nlv		!maximum level
	iunit = 3
	it = 0			!time stamp - not important for initial cond.
	bformat = .true.	!formatted or unformatted?

	call get_data(nxdim,nydim,nzdim,nx,ny,nz
     +				,level,uvel,vvel,temp,salt)

	open(iunit,file='level.dat',form='formatted',status='unknown')
	call intp_data2d(nxdim,nydim,nx,ny,level,val2d)
        call write_fem_file_2d(bformat,iunit,it,nkn,'water level',val2d)
	close(iunit)

	open(iunit,file='velocity.dat',form='formatted',status='unknown')
	nvar = 2
        call write_fem_file_header(bformat,iunit,it,nkn,lmax,nvar
     +					,nvers,ntype,hlv)
	call intp_data3d(nxdim,nydim,nzdim,nx,ny,nz,uvel,val3d)
        call write_fem_file_data(bformat,iunit
     +                          ,nkn,lmax,nvers,nlvdim,ilhkv
     +				,'x-velocity',val3d)
	call intp_data3d(nxdim,nydim,nzdim,nx,ny,nz,vvel,val3d)
        call write_fem_file_data(bformat,iunit
     +                          ,nkn,lmax,nvers,nlvdim,ilhkv
     +				,'y-velocity',val3d)
	close(iunit)

	open(iunit,file='temp.dat',form='formatted',status='unknown')
	call intp_data3d(nxdim,nydim,nzdim,nx,ny,nz,temp,val3d)
        call write_fem_file_3d(bformat,iunit,it,nkn,nlv,nlvdim,hlv
     +                          ,ilhkv,'temperature',val3d)
	close(iunit)

	open(iunit,file='salt.dat',form='formatted',status='unknown')
	call intp_data3d(nxdim,nydim,nzdim,nx,ny,nz,temp,val3d)
        call write_fem_file_3d(bformat,iunit,it,nkn,nlv,nlvdim,hlv
     +                          ,ilhkv,'salinity',val3d)
	close(iunit)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine set_vertical(hev,nlv,hlv,ilhv,ilhkv)

c sets the vertical structure

	implicit none

	include 'param.h'
	include 'basin.h'

	real hev(neldim)		!depth at elements (input)
	integer nlv			!total number of levels (return)
	real hlv(nlvdim)		!vertical levels (return)
	integer ilhv(neldim)		!number of levels in element (return)
	integer ilhkv(nkndim)		!number of levels in node (return)

	integer ie,k,ii,l,lmax
	real h,hmax

c--------------------------------------------------------
c the next array must be customized for the application
c--------------------------------------------------------

	integer ndim
	parameter(ndim=10)
	real hl(ndim)
	data hl /5,10,15,20,30,40,50,70,100,150/

c--------------------------------------------------------
c check max depth and assign to vertical level structure
c--------------------------------------------------------

	hmax = 0.
	do ie=1,nel
	  hmax = max(hmax,hev(ie))
	end do
	lmax = min(ndim,nlvdim)
	if( hmax .gt. hl(lmax) ) goto 99

	do l=1,lmax
	  hlv(l) = hl(l)
	end do

c--------------------------------------------------------
c compute max levels in element
c--------------------------------------------------------

	do ie=1,nel
	  h = hev(ie)
	  do l=1,lmax
	    if( h .le. hl(l) ) goto 1
	  end do
    1	  continue
	  if ( l .gt. lmax ) goto 98
	  ilhv(ie) = l
	end do

c--------------------------------------------------------
c compute max levels in node
c--------------------------------------------------------

	do k=1,nkn
	  ilhkv(k) = 0
	end do

	do ie=1,nel
	  l = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ilhkv(k) = max(ilhkv(k),l)
	  end do
	end do
	    
	do k=1,nkn
	  if( ilhkv(k) .le. 0 ) goto 97
	end do

	nlv = lmax

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	return
   97	continue
	write(6,*) k,ilhkv(k)
	stop 'error stop set_vertical: ilhkv (internal error)'
   98	continue
	write(6,*) lmax,hmax,hl(lmax),h
	stop 'error stop set_vertical: hmax (internal error)'
   99	continue
	write(6,*) lmax,hmax,hl(lmax)
	stop 'error stop set_vertical: hmax'
	end

c*******************************************************************

	subroutine intp_data2d(nxdim,nydim,nx,ny,data,val2d)

c this subroutines interpolates regular data to a FEM grid (2d)

	implicit none

	include 'param.h'
	include 'basin.h'

	integer nxdim,nydim
	integer nx,ny
	real data(nxdim,nydim)
	real val2d(1)

	integer k
	real xfem,yfem

	do k=1,nkn
	  xfem = xgv(k)
	  yfem = ygv(k)

	  ! here do interpolation 

	end do

	end

c*******************************************************************

	subroutine intp_data3d(nxdim,nydim,nzdim,nx,ny,nz,data,val3d)

c this subroutines interpolates regular data to a FEM grid (3d)
c the vertical information has also to be passed in and be used

	implicit none

	include 'param.h'
	include 'basin.h'

	integer nxdim,nydim,nzdim
	integer nx,ny,nz
	real data(nzdim,nxdim,nydim)
	real val3d(nlvdim,1)

	integer k
	real xfem,yfem

	do k=1,nkn
	  xfem = xgv(k)
	  yfem = ygv(k)

	  ! here do interpolation -> must also handle vertical discretization

	end do

	end

c*******************************************************************

	subroutine get_data(nxdim,nydim,nzdim,nx,ny,nz
     +				,level,uvel,vvel,temp,salt)

c reads regular data to be interpolated

	implicit none

	integer nxdim,nydim,nzdim
	integer nx,ny,nz
	real level(nxdim,nydim)
	real uvel(nzdim,nxdim,nydim)
	real vvel(nzdim,nxdim,nydim)
	real temp(nzdim,nxdim,nydim)
	real salt(nzdim,nxdim,nydim)

	! here read the data

	end

c*******************************************************************





