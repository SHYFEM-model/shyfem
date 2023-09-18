!
! $Id: subreg.f,v 1.17 2009-09-14 08:20:58 georg Exp $
!
! routines for interpolation onto regular grid
!
! contents :
!
! subroutine setgeo(x0,y0,dx,dy,flag)
!		sets grid values for interpolation
! subroutine getgeo(x0,y0,dx,dy,flag)
!		gets grid values for interpolation
! subroutine getgeoflag(flag)
!		gets flag value for interpolation
!
! subroutine av2am(av,am,ip,jp)
!		interpolation of av onto a regular net
! subroutine av2amk(bwater,av,am,ip,jp)
!		interpolation of av onto a regular net
! function intri(x,y,xp,yp)
!		point in triangle or not
! function intrid(x,y,xp,yp)
!		point in triangle or not (double precision version)
! subroutine am2av(am,av,ip,jp)
!		interpolation of am onto finite element mesh
! function am2val(am,ip,jp,xx,yy)
!		interpolation of am onto finite element mesh
! subroutine ave2am(av,am,ip,jp)
!		interpolation of av (elementwise) onto a regular net
!
! subroutine mkmask(bwater,zv,href,hzoff)
!		makes mask in element
!
! subroutine mimareg(am,ip,jp,amin,amax)
!		computes min/max of regular matrix (without flag values)
! subroutine a2char(am,ac,ip,jp)
!		creates 1 char representation of matrix
! subroutine prchar(ac,ip,jp)
!		prints 1 char representation of matrix
!
! subroutine femintp(ie,z,xp,yp,zp)
!               interpolation in element (with ev)
! subroutine elemintp(x,y,z,xp,yp,zp)
!               interpolation in element (no ev)
!
! subroutine find_elem_from_old(ieold,xp,yp,ielem)
!		finds element for point (xp,yp) starting from ieold
! subroutine find_element(xp,yp,ielem)
!		finds element for point (xp,yp)
! function in_element(ie,xp,yp)
!		checks if point (xp,yp) is in element ie
! subroutine get_xy_elem(ie,x,y)
!		returns x,y of vertices of element ie
!
! revision log :
!
! 18.11.1998	ggu	routine commented
! 18.11.1998	ggu	routine setgeo introduced
! 19.11.1998	ggu	routines a2char, prchar added
! 19.10.1999	ggu	routine mkmask added from subutl
! 25.11.2004	ggu	new routines femintp and elemintp for interpolation
! 14.03.2005	ggu	new routines for interpolation in element
! 11.03.2009	ggu	new helper routine getgeoflag()
! 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
! 26.01.2011	ggu&mb	handling extrapolation in am2av()
! 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
! 31.03.2011	ggu	new routine elemmask()
! 24.11.2011	ggu	new routine find_close_elem()
! 20.06.2012	ggu	new routine get_scal_elem()
! 07.10.2012	ggu	new routine av2fm()
! 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
! 26.10.2012	ggu	bug fix: do not access not existing storage
! 30.05.2014	ggu	in av2amk() do not interpolate for flag values
! 07.07.2014	ggu	new routine intp_reg()
! 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
! 05.05.2016	ggu	file restructured (module)
! 14.05.2016	ggu	allow for extension of grid -> bregextend
! 23.06.2016	ggu	allow for eps in computing box
!
! notes :
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points
!
!******************************************************

!==================================================================
        module regular_2nc
!==================================================================

        use basin
        use evgeom_2nc
        use hydro_admin
        use geom

	implicit none

	logical, save :: bregextend = .false.

	double precision, save :: pxareg = 0.	!x coordinate of lower,left point (x0)
	double precision, save :: pyareg = 0.	!y coordinate of lower,left point (y0)
	double precision, save :: pxdreg = 0.	!grid spacing in x direction (dx)
	double precision, save :: pydreg = 0.	!grid spacing in y direction (dy)
	double precision, save :: pzlreg = -999.	!flag for land points

        contains

!==================================================================

	subroutine setgeo(x0,y0,dx,dy,flag)

! sets grid values for interpolation
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	implicit none

#ifdef SINGLEP
	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points
#else
	double precision x0,y0	!coordinates of lower,left point
	double precision dx,dy	!grid spacing in x/y direction
	double precision flag	!flag for land points
#endif

	pxareg = x0
	pyareg = y0
	pxdreg = dx
	pydreg = dy
	pzlreg = flag

	end

!******************************************************

	subroutine getgeo(x0,y0,dx,dy,flag)

! gets grid values for interpolation
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	implicit none

	double precision x0,y0	!coordinates of lower,left point
	double precision dx,dy	!grid spacing in x/y direction
	double precision flag	!flag for land points

	x0   = pxareg
	y0   = pyareg
	dx   = pxdreg
	dy   = pydreg
	flag = pzlreg

	end

!******************************************************

	subroutine getgeoflag(flag)

! gets flag value for interpolation
!
! flag                value for land points

	implicit none

	double precision flag	!flag for land points

	flag = pzlreg

	end

!******************************************************

	subroutine setregextend(bextend)

! sets flag to decide if extend interpolated grid

	implicit none

	logical bextend

	bregextend = bextend

	end

!******************************************************

	subroutine getregextend(bextend)

! gets flag to decide if extend interpolated grid

	implicit none

	logical bextend

	bextend = bregextend

	end

!******************************************************

	subroutine setreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	double precision regpar(7)
	integer nx,ny
	double precision x0,y0,dx,dy
	double precision flag

	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x0
	regpar(4) = y0
	regpar(5) = dx
	regpar(6) = dy
	regpar(7) = flag

	end

!******************************************************

	subroutine getreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	double precision regpar(7)
	integer nx,ny
	double precision x0,y0,dx,dy
	double precision flag

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	end

!******************************************************
!******************************************************
!******************************************************

	subroutine find_position_to_coord(x,y,ix,iy)

! finds closest position (ix,iy) to coordinate (x,y) in reg grid

	implicit none

	double precision x,y
	integer ix,iy

	double precision x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	ix = nint( (x-x0)/dx + 1. )
	iy = nint( (y-y0)/dy + 1. )

	end

!******************************************************

	subroutine find_coord_to_position(ix,iy,x,y)

! finds coordinate (x,y) to given position (ix,iy) in reg grid

	implicit none

	integer ix,iy
	double precision x,y

	double precision x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	x = x0 + (ix-1)*dx
	y = y0 + (iy-1)*dy

	end

!******************************************************
!******************************************************
!******************************************************

	subroutine av2am(av,am,ip,jp)

! interpolation of av onto a regular net (nodal values)
!
! av                    array to be interpolated
! am                    matrices of interpolated values (u,v,z,h)
! ip,jp                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	implicit none

	integer ip,jp
	double precision av(nkn)
	double precision am(ip,jp)

	logical bwater(nel)

	bwater = .true.

	call av2amk(bwater,av,am,ip,jp)

	end

!******************************************************

	subroutine av2amk(bwater,av,am,ip,jp)

! interpolation of av onto a regular net (nodal values) with mask
!
! bwater		mask for water points
! av                    array to be interpolated
! am                    matrices of interpolated values (u,v,z,h)
! ip,jp                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	implicit none

! arguments
	integer ip,jp
	double precision av(nkn)
	double precision am(ip,jp)
	logical bwater(nel)
! parameter
	double precision eps
	parameter ( eps = 1.d-14 )
! local
	double precision pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	integer iflag
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	am=pzlreg

	do ie=1,nel
	  if( bwater(ie) ) then	!wet
	    iflag = 0
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
		z(i)=av(kn)
	        if( z(i) > pzlreg ) iflag = iflag + 1
	    end do
	    if( iflag .ne. 3 ) cycle

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			zh=0.
			do k=1,3
			   fh=(a(k)+xp*b(k)+yp*c(k))/f
			   zh=zh+z(k)*fh
			end do
			am(i,j)=zh
		    end if
		end do
	    end do
	  end if
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2amk: area of element'
	end

!************************************************
!************************************************
!************************************************

	subroutine av2fm(fm,ip,jp)

! computation of interpolation matrix (nodal values to regular grid) with mask

	implicit none

#ifdef SINGLEP
	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
#else
	double precision fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
#endif
	integer ip,jp		!dimension of matrices

	logical bwater(nel)	!wet mask for each element

	bwater = .true.

	call av2fmk(bwater,fm,ip,jp)

	end

!************************************************

	subroutine av2fmk(bwater,fm,ip,jp)

! computation of interpolation matrix (nodal values to regular grid) with mask
!
! the interpolation can be carried out as
!
!	do j=1,jp
!	  do i=1,ip
!	    ie = nint(fm(4,i,j))
!	    if( ie .gt. 0 ) then
!	      a = 0.
!	      do ii=1,3
!	        k = nen3v(ii,ie)
!		a = a + val(k) * fm(ii,i,j)
!	      end do
!	    else
!	      a = flag
!	    end if
!	    am(i,j) = a
!	  end do
!	end do
	        
	implicit none

! arguments
	logical bwater(nel)	!wet mask for each element
#ifdef SINGLEP
	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
#else
	double precision fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
#endif
	integer ip,jp		!dimension of matrices
! parameter
	double precision eps
	parameter ( eps = 1.d-14 )
! local
	double precision pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	logical bok
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	fm = 0.

	do ie=1,nel
	  bok = bwater(ie)
	  if( bok ) then			!wet
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
	    end do

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			do ii=1,3
			   fh=(a(ii)+xp*b(ii)+yp*c(ii))/f
			   fm(ii,i,j) = fh
			end do
			fm(4,i,j) = ie
		    end if
		end do
	    end do
	  end if
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2fm: area of element'
	end

!************************************************

        subroutine fm2am2d(femval,nx,ny,fm,am)

! interpolation 2d of fem values to regular grid using fm matrix

        implicit none

        integer nx,ny			!dimension of regular matrix

	integer nlvdi,nlv
	integer ilhv(nel)
        double precision femval(nkn)		!values of fem array
#ifdef SINGLEP
        real fm(4,nx,ny)		!interpolation matrix
        real am(nx,ny)			!interpolated values (return)
#else
        double precision fm(4,nx,ny)		!interpolation matrix
        double precision am(nx,ny)			!interpolated values (return)
#endif

	nlvdi = 1
	nlv = 1
	ilhv = 1

        call fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

	end

!************************************************

        subroutine fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

! interpolation 3d of fem values to regular grid using fm matrix

        implicit none

	integer nlvdi			!vertical dimension of fem array
        integer ilhv(nel)		!vertical discretization (element!!)
        integer nlv,nx,ny		!dimension of regular matrix
        double precision femval(nlvdi,nkn)		!values of fem array

#ifdef SINGLEP
        real fm(4,nx,ny)		!interpolation matrix
        real am(nlv,nx,ny)		!interpolated values (return)
#else
        double precision fm(4,nx,ny)		!interpolation matrix
        double precision am(nlv,nx,ny)		!interpolated values (return)
#endif

        integer i,j,l,lmax,ie,ii,k
        double precision a
        double precision flag

	call getgeoflag(flag)

        do j=1,ny
          do i=1,nx
            ie = nint(fm(4,i,j))
            lmax = 0
            if( ie .gt. 0 ) lmax = ilhv(ie)
	    lmax = min(lmax,nlv)
            do l=1,lmax
              a = 0.
              do ii=1,3
                k = nen3v(ii,ie)
                a = a + femval(l,k) * fm(ii,i,j)
              end do
              am(l,i,j) = a
            end do
            do l=lmax+1,nlv
              am(l,i,j) = flag
            end do
          end do
        end do

        end

!************************************************
!************************************************
!************************************************

	subroutine fm_extra_setup(nx,ny,fmextra)

! sets up fmextra structure to allow interpolation from fem nodes to reg grid

	implicit none

        integer nx,ny			!dimension of regular matrix
	double precision fmextra(6,nkn)

	logical bout,berror
	integer ix,iy,iix,iiy
	integer j,k
	double precision x0,y0,dx,dy
	double precision x,y,x1,y1
	double precision eps,flag
	double precision t,u
	double precision d,w,d2

	integer, save :: jx(4) = (/0,1,1,0/)
	integer, save :: jy(4) = (/0,0,1,1/)
	double precision fmweight(nx,ny)

	eps = 0.01
	fmextra = 0.
	fmweight = 0.

	call getgeo(x0,y0,dx,dy,flag)

!	---------------------------------------------------------
!	set up contribution from each fem node to regular grid
!	---------------------------------------------------------

	berror = .false.
	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  ix = (x-x0)/dx+1.
	  iy = (y-y0)/dy+1.
	  if( ix < 1 .or. ix >= nx ) cycle
	  if( iy < 1 .or. iy >= ny ) cycle
	  x1 = x0+(ix-1)*dx
	  y1 = y0+(iy-1)*dy
	  t = (x-x1)/dx
	  u = (y-y1)/dy
	  bout = .false.
	  if( t.gt.1. .or. t.lt.0. ) bout = .true.
	  if( u.gt.1. .or. u.lt.0. ) bout = .true.
	  if( bout ) then
	    write(6,*) 'out of domain: ',k,ix,iy,x0,y0,x1,y1,x,y,t,u
	    berror = .true.
	  end if
	  fmextra(1,k) = ix
	  fmextra(2,k) = iy
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      x1 = x0+(iix-1)*dx
	      y1 = y0+(iiy-1)*dy
	      d2 = ((x1-x)/dx)**2 + ((y1-y)/dy)**2	!normalized distance
	      d = sqrt( d2 )
	      if( d2 > 2. ) then
		write(6,*) 'distance too large: ',ix,iy,iix,iiy,d,d2
		berror = .true.
	      end if
	      !w = 2. - d			!weight - could be gaussian
	      w = exp(-d2/2.)			!sigma is 1
	      fmextra(2+j,k) = w
	      fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    !end if
	  end do
	end do

	if( berror ) then
	  stop 'error stop fm_extra_setup: internal error (3)'
	end if

!	---------------------------------------------------------
!	scale weight to 1
!	---------------------------------------------------------

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      w = fmweight(iix,iiy)
	      if( w > 0. ) fmextra(2+j,k) = fmextra(2+j,k) / w
	    !end if
	  end do
	end do

!	---------------------------------------------------------
!	check if weight sums up to 1
!	---------------------------------------------------------

	fmweight = 0.

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      w = fmextra(2+j,k)
	      fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    !end if
	  end do
	end do

	berror = .false.
	do iy=1,ny
	  do ix=1,nx
	    w = fmweight(ix,iy)
	    if( w > 0 ) then
	      if( abs(w-1.) > eps ) then
		berror = .true.
		write(6,*) 'error... ',ix,iy,w
	      end if
	    end if
	  end do
	end do

	if( berror ) then
	  stop 'error stop fm_extra_setup: internal error (2)'
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************

	subroutine fm_extra_3d(nlvdi,nlv,il,nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	implicit none

	integer nlvdi,nlv
	integer il(nkn)
        integer nx,ny			!dimension of regular matrix
	double precision fmextra(6,nkn)
	double precision femdata(nlvdi,nkn)
	double precision regdata(nlvdi,nx,ny)

	integer k,l,lmax
	double precision flag
	double precision fem2d(nkn)
	double precision reg2d(nx,ny)

	call getgeoflag(flag)

!	---------------------------------------------------------
!	make sure fem data below bottom is flag
!	---------------------------------------------------------

	do k=1,nkn
	  lmax = il(k)
	  femdata(lmax+1:nlvdi,k) = flag
	end do

!	---------------------------------------------------------
!	interpolate layer by layer
!	---------------------------------------------------------

	do l=1,nlv
	  fem2d(:) = femdata(l,:)
	  reg2d(:,:) = regdata(l,:,:)
	!write(6,*) l,nx,ny,nx*ny
	!write(6,*) (fem2d(k),k=1,nkn,nkn/20)
	!write(6,*) 'before'
	!write(6,*) reg2d
	  call fm_extra_2d(nx,ny,fmextra,fem2d,reg2d)
	  regdata(l,:,:) = reg2d(:,:)
	!write(6,*) 'after'
	!write(6,*) reg2d
	end do

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************

	subroutine fm_extra_2d(nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	implicit none

        integer nx,ny			!dimension of regular matrix
	double precision fmextra(6,nkn)
	double precision femdata(nkn)
	double precision regdata(nx,ny)

	integer ix,iy,iix,iiy
	integer j,k
	double precision x0,y0,dx,dy
	double precision eps,flag
	double precision regval,femval
	double precision d,w

	integer, save :: jx(4) = (/0,1,1,0/)
	integer, save :: jy(4) = (/0,0,1,1/)
	double precision fmweight(nx,ny)
	double precision fmdata(nx,ny)

	eps = 0.01
	fmweight = 0.
	fmdata = 0.

	call getgeoflag(flag)

!	---------------------------------------------------------
!	accumulate on regular grid (only where flag is set)
!	---------------------------------------------------------

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  femval = femdata(k)
	  if( femval == flag ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    regval = regdata(iix,iiy)
	    if( regval /= flag ) cycle
	    w = fmextra(2+j,k)
	!write(6,*) ix,iy,iix,iiy,femval,regval,w
	!write(6,*) ix,iy,iix,iiy,femval,regval
	    fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    fmdata(iix,iiy) = fmdata(iix,iiy) + w * femval
	  end do
	end do
	!write(6,*) 'fmweight'
	!write(6,*) fmweight
	!write(6,*) regdata

!	---------------------------------------------------------
!	correct for weight and set where flag
!	---------------------------------------------------------

	where ( fmweight > 0. ) 
	  fmdata = fmdata / fmweight
	else where
	  fmdata = flag
	end where
	!write(6,*) fmdata
	where ( regdata == flag ) regdata = fmdata

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************
!************************************************
!************************************************

	function intri(x,y,xp,yp)

! point in triangle or not
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)

	implicit none

! arguments
	integer intri
	double precision x(3),y(3),xp,yp
! local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
! save
	save eps
	data eps /1.e-13/

	intri=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intri=1	!inside

	end

!************************************************

	function intrid(x,y,xp,yp)

! point in triangle or not (double precision version)
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)

	implicit none

! arguments
	integer intrid
	double precision x(3),y(3),xp,yp
! local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
! save
	save eps
	data eps /1.e-13/

	intrid=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intrid=1	!inside

	end

!****************************************************************
!
	function intri0(x,y,xp,yp)
!
! point in triangle or not
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)
!
	implicit none
!
! arguments
	integer intri0
	double precision x(3),y(3),xp,yp
! local
	integer i,k1,k2
	double precision xs,ys,x12,y12
	double precision det,detlam,rlamb
	double precision eps
! save
	save eps
	data eps /1.e-13/
!
	xs=0.
	ys=0.
	do i=1,3
	   xs=xs+x(i)
	   ys=ys+y(i)
	end do
	xs=xs/3.
	ys=ys/3.
!
	intri0=0
!
	do k1=1,3
	   k2=mod(k1,3)+1
	   x12=x(k1)-x(k2)
	   y12=y(k1)-y(k2)
	   det=(xs-xp)*y12-(ys-yp)*x12
	   if(abs(det).ge.eps) then
		detlam=(x(k1)-xp)*y12-(y(k1)-yp)*x12
		rlamb=detlam/det
		if(rlamb.gt.0..and.rlamb.lt.1.) return	!outside
	   end if
	end do
!
	intri0=1	!inside
!
	return
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine intp_reg_nodes(nx,ny,x0,y0,dx,dy,flag,regval,femval,ierr)

! interpolates regular grid to FEM grid - values are on nodes

	implicit none

	integer nx,ny
	double precision x0,y0,dx,dy
	double precision flag
	double precision regval(nx,ny)
	double precision femval(nkn)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval,nkn,xgv,ygv,femval,ierr)

	end

!****************************************************************

	subroutine intp_reg_elems(nx,ny,x0,y0,dx,dy,flag,regval,femval,ierr)

! interpolates regular grid to FEM grid - values are on elements

	implicit none

	integer nx,ny
	double precision x0,y0,dx,dy
	double precision flag
	double precision regval(nx,ny)
	double precision femval(nel)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	double precision xp(nel),yp(nel)

	call intp_reg_make_cg(nel,xp,yp)
	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval,nel,xp,yp,femval,ierr)

	end

!****************************************************************

	subroutine intp_reg_make_cg(np,xp,yp)

	implicit none

	integer np
	double precision xp(np),yp(np)

	integer ie,ii,k
	double precision x,y

	if( np /= nel ) then
	  write(6,*) 'np,nel: ',np,nel
	  stop 'error stop intp_reg_make_cg: parameters'
	end if

	do ie=1,nel
	  x = 0.
	  y = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x = x + xgv(k)
	    y = y + ygv(k)
	  end do
	  xp(ie) = x / 3.
	  yp(ie) = y / 3.
	end do

	end

!****************************************************************

	subroutine intp_reg(nx,ny,x0,y0,dx,dy,flag,regval,np,xp,yp,femval,ierr)

! interpolation of regular array onto fem grid - general routine
!
! ierr:
!		= 0	no errors
!		< 0	interpolation out of domain (extrapolation)
!		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	double precision x0,y0,dx,dy
	double precision flag
	double precision regval(nx,ny)
	integer np		!number of fem points
	double precision xp(np)
	double precision yp(np)
	double precision femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	logical bextra,bout,bflag
	logical bintpout,bintpflag,bextend
	integer k
	integer imin,jmin
	integer iflag,iout
	double precision xx,yy,z1,z2,z3,z4,x1,y1,t,u
	double precision xn,yn
	double precision zz(4)
 
	double precision, parameter :: eps = 1.e-4
	!double precision, parameter :: eps = 0.
	logical outbox
	outbox(t) = ( t-1. > eps .or. t < -eps )

	!bintpout = .false.	!interpolate even if outside
	!bintpout = .true.	!interpolate even if outside
	!bintpflag = .false.	!interpolate even if flag
	!bintpflag = .true.	!interpolate even if flag

	call getregextend(bextend)
	bintpout = bextend
	bintpflag = bextend

	iflag = 0	!used flag for interpolation
	iout = 0	!used outside point for interpolation

	imin = 0
	jmin = 0

	xn = x0 + (nx-1)*dx
	yn = y0 + (ny-1)*dy

	do k=1,np
	    xx = xp(k)
	    yy = yp(k)
 
	    femval(k) = flag
 
	    if( xx .le. x0 ) then
	      imin = 1
	    else if( xx .ge. xn ) then
	      imin = nx-1
	    else
	      imin=1+(xx-x0)/dx
	    end if
	    if( yy .le. y0 ) then
	      jmin = 1
	    else if( yy .ge. yn ) then
	      jmin = ny-1
	    else
	      jmin=1+(yy-y0)/dy
	    end if

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 99
	    if( imin+1.gt.nx .or. jmin+1.gt.ny ) goto 99

	    x1 = x0+(imin-1)*dx
	    y1 = y0+(jmin-1)*dy
	    t = (xx-x1)/dx
	    u = (yy-y1)/dy

	    bout = .false.
	    if( outbox(t) ) bout = .true.
	    if( outbox(u) ) bout = .true.
	    if( bout ) then
	      if( bintpout ) then
	        if( t .le. 2. ) t = min(1.,t)
	        if( u .le. 2. ) u = min(1.,u)
	        if( t .ge. -1. ) t = max(0.,t)
	        if( u .ge. -1. ) u = max(0.,u)
	      end if
	      bout = .false.
	      if( outbox(t) ) bout = .true.
	      if( outbox(u) ) bout = .true.
	      if( bout ) then
		iout = iout + 1
		cycle
	      end if
	    end if

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)
	    zz = (/z1,z2,z3,z4/)

	    if( any(zz == flag) .and. bintpflag ) then
	      call recover_flag(zz,z1,z2,z3,z4,flag)
	    end if

	    bflag = .false.
	    if( z1.eq.flag .or. z2.eq.flag ) bflag = .true.
	    if( z3.eq.flag .or. z4.eq.flag ) bflag = .true.
	    if( bflag ) then
	      iflag = iflag + 1
	      cycle
	    end if

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout .gt. 0 ) ierr = - iout - iflag
	if( iflag .gt. 0 ) ierr = iflag

	!write(6,*) 'intp_reg: ierr = ',ierr,iout,iflag

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

!****************************************************************

	subroutine recover_flag(zz,z1,z2,z3,z4,flag)

	implicit none

	double precision zz(4)
	double precision z1,z2,z3,z4
	double precision flag

	integer i,ic
	double precision zt

	!write(6,*) 'recovering flag: ',zz

	ic = count(zz == flag)
	if( ic == 4 ) return

	zt = 0.
	do i=1,4
	  if( zz(i) /= flag ) zt = zt + zz(i)
	end do
	zt = zt / (4-ic)
	do i=1,4
	  if( zz(i) == flag ) zz(i) = zt
	end do

	z1 = zz(1)
	z2 = zz(2)
	z3 = zz(3)
	z4 = zz(4)

	!write(6,*) 'recovered flag: ',zz

	end

!****************************************************************

	subroutine intp_reg_setup_fr(nx,ny,x0,y0,dx,dy,np,fr)

! interpolation of regular array onto fem grid - general routine
!
! produces array fr that can be used to interpolate
!
! works for np equal to nkn or nel

	implicit none

	integer nx,ny
	double precision x0,y0,dx,dy
	integer np		!number of fem points
	double precision fr(4,np)		!array for interpolation on fem grid (return)

	integer k
	integer imin,jmin
	double precision xx,yy,x1,y1,t,u
	double precision xn,yn
	double precision xp(np)
	double precision yp(np)
 
	if( np == nkn ) then
	  xp = xgv
	  yp = ygv
	else if( np == nel ) then
	  call intp_reg_make_cg(nel,xp,yp)
	else
	  write(6,*) 'np,nkn,nel: ',np,nkn,nel
	  stop 'error stop intp_reg_setup_fr: np'
	end if

	imin = 0
	jmin = 0
	fr = 0.

	xn = x0 + (nx-1)*dx
	yn = y0 + (ny-1)*dy

	do k=1,np
	    xx = xp(k)
	    yy = yp(k)
 
	    if( xx .le. x0 ) then
	      imin = 1
	    else if( xx .ge. xn ) then
	      imin = nx-1
	    else
	      imin=1+(xx-x0)/dx
	    end if
	    if( yy .le. y0 ) then
	      jmin = 1
	    else if( yy .ge. yn ) then
	      jmin = ny-1
	    else
	      jmin=1+(yy-y0)/dy
	    end if

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 99
	    if( imin+1.gt.nx .or. jmin+1.gt.ny ) goto 99

	    x1 = x0+(imin-1)*dx
	    y1 = y0+(jmin-1)*dy
	    t = (xx-x1)/dx
	    u = (yy-y1)/dy

	    fr(1,k) = imin
	    fr(2,k) = jmin
	    fr(3,k) = t
	    fr(4,k) = u
	end do
 
	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

!****************************************************************

	subroutine intp_reg_intp_fr(nx,ny,flag,regval,np,fr,femval,ierr)

! interpolation of regular array onto fem grid - general routine
!
! ierr:
!		= 0	no errors
!		< 0	interpolation out of domain (extrapolation)
!		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	double precision flag
	double precision regval(nx,ny)
	integer np		!number of fem points
	double precision fr(4,np)
	double precision femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	integer k
	integer imin,jmin
	integer iflag,iout
	double precision z1,z2,z3,z4,t,u
 
	iflag = 0	!used flag for interpolation
	iout = 0	!used outside point for interpolation

	do k=1,np
 
	    femval(k) = flag
 
	    imin = nint(fr(1,k))
	    jmin = nint(fr(2,k))
	    t = fr(3,k)
	    u = fr(4,k)

	    iout = 0
	    if( u.gt.1. .or. u.lt.0. ) iout = iout + 1
	    if( t.gt.1. .or. t.lt.0. ) iout = iout + 1
	    if( iout > 0 ) cycle

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)

	    iflag = 0
	    if( z1.eq.flag .or. z2.eq.flag ) iflag = iflag + 1
	    if( z3.eq.flag .or. z4.eq.flag ) iflag = iflag + 1
	    if( iflag > 0) cycle

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout .gt. 0 ) ierr = - iout - iflag
	if( iflag .gt. 0 ) ierr = iflag

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

!****************************************************************

	subroutine am2av(am,av,ip,jp)

! compatibility for old calls

	implicit none

	integer ip,jp
	double precision av(nkn)
	double precision am(ip,jp)

	integer ierr
	double precision pxareg,pyareg,pxdreg,pydreg,pzlreg

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	call intp_reg(ip,jp,pxareg,pyareg,pxdreg,pydreg,pzlreg,am,nkn,xgv,ygv,av,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'ierr = ',ierr
	  stop 'error stop am2av: error in interpolation'
	end if

	end

!******************************************************
!******************************************************
!******************************************************

!******************************************************
!******************************************************
!******************************************************

	subroutine set_dry_mask(bwater,zv,href,hzoff)

! makes mask for dry and wet areas - zenv must be available
!
! bwater is elementwise mask:	true = water point


	implicit none

! arguments
	logical bwater(nel)
	double precision zv(nkn)
	double precision href,hzoff
! local
	integer itot,itot1
	integer ie,ii

        do ie=1,nel

          itot=0
          do ii=1,3
            if( hm3v(ii,ie)+zenv(ii,ie)-href .gt. hzoff ) then
		itot=itot+1    !wet
	    end if
          end do

          itot1=0
          do ii=1,3
            if(zv(nen3v(ii,ie)).eq.zenv(ii,ie)) itot1=itot1+1
          end do

          !if(itot.ne.3.or.itot1.ne.3)  bwater(ie) = .false.
          if(itot.ne.3.and.itot1.ne.3)  bwater(ie) = .false.

        end do

	end

!******************************************************

	subroutine set_level_mask(bwater,ilhv,level)

! makes mask for water points (level)
!
! bwater is elementwise mask:	true = water point

	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	logical bwater(nel)
	integer ilhv(nel)
	integer level
! local
	integer ie,nedry

	nedry = 0

	do ie=1,nel
	  if( level .gt. ilhv(ie) ) then
	    bwater(ie) = .false.
	    nedry = nedry + 1
	  end if
	end do

	write(6,*) 'level exist  (dry/wet/total): ',nedry,nel-nedry,nel

	end

!******************************************************

	subroutine make_dry_node_mask(bwater,bkwater)

! makes node mask from element mask
!
! bwater is elementwise mask:	true = water point

	implicit none

! arguments
	logical bwater(nel)
	logical bkwater(nkn)

	integer ie,ii,k
	integer nndry,nedry

	nndry = 0
	nedry = 0

	bkwater = .false.

	do ie=1,nel
	  if( bwater(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkwater(k) = .true.
	    end do
	  else
	  end if
	end do

	end

!******************************************************

	subroutine make_dry_elem_mask(bwater,bkwater)

! makes elem mask from node mask
!
! bwater is elementwise mask:	true = water point

	implicit none

! arguments
	logical bwater(nel)
	logical bkwater(nkn)
! local
	integer ie,ii,k
	integer nedry

        bwater = .true.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( .not. bkwater(k) ) then
	      bwater(ie) = .false.
	    end if
	  end do
	end do

	end

!******************************************************

        subroutine info_dry_mask(bwater,bkwater)

        implicit none

	logical bwater(nel)
	logical bkwater(nkn)
        
        integer nedry,nndry,newet,nnwet

        newet = count(bwater)
        nnwet = count(bkwater)
        nedry = nel - newet
        nndry = nkn - nnwet

	write(6,*) 'dry elements (dry/wet/total): ',nedry,newet,nel
	write(6,*) 'dry nodes    (dry/wet/total): ',nndry,nnwet,nkn

        end

!******************************************************
!******************************************************
!******************************************************

	subroutine mimareg(am,ip,jp,amin,amax)

! computes min/max of regular matrix (without flag values)

	implicit none

	integer ip,jp
	double precision am(ip,jp)
	double precision amin,amax

	integer i,j
	double precision a,high,flag

	high = 1.e+30
	call getgeoflag(flag)

	amin =  high
	amax = -high

	do j=1,jp
	  do i=1,ip
	    a = am(i,j)
	    if( a .ne. flag ) then
	      if( a .gt. amax ) amax = a
	      if( a .lt. amin ) amin = a
	    end if
	  end do
	end do

	if( amin .eq.  high ) amin = flag
	if( amax .eq. -high ) amax = flag

	end

!******************************************************

	subroutine a2char(am,ac,ip,jp)

! creates 1 char representation of matrix

	implicit none

	integer ip,jp		!dimension of matrix
	double precision am(ip,jp)		!matrix containing data
	character*1 ac(ip,jp)	!matrix containing chars on return

	integer i,j
	double precision flag

	call getgeoflag(flag)

	do j=1,jp
	  do i=1,ip
	    if( am(i,j) .eq. flag ) then	!land
		ac(i,j) = '.'
	    else
		ac(i,j) = '*'			!data
	    end if
	  end do
	end do

	end

!******************************************************

	subroutine prchar(ac,ip,jp)

! prints 1 char representation of matrix

	implicit none

	integer ip,jp		!dimension of matrix
	character*1 ac(ip,jp)	!matrix containing chars

	character*256 line
	integer ipmax
	integer i,j

	ipmax = min(256,ip)

	do j=jp,1,-1
	  line = ' '
	  do i=1,ipmax
	    line(i:i) = ac(i,j)
	  end do
	  write(6,'(a)') line(1:ipmax)
	end do

	end

!******************************************************
!******************************************************
!******************************************************

        subroutine femintp(ie,z,xp,yp,zp)

! interpolation in element (with ev)
!
! interpolates in element ie from nodal values z to point xp,yp
! result is in zp
!
! needs array ev

        integer ie      !element
        double precision z(3)       !values on nodes
        double precision xp,yp      !coordinates of point
        double precision zp         !interpolated value (return)

        integer ii
        double precision zh,a,b,c,w

        zh=0.
        do ii=1,3
          a = ev(ii,ie)
          b = ev(3+ii,ie)
          c = ev(6+ii,ie)
          w = a + b*xp + c*yp
          zh = zh + z(ii) * w
        end do

        zp = zh

        end

!******************************************************

        subroutine elemintp(x,y,z,xp,yp,zp)

! interpolation in element (no ev)
!
! interpolates in element given by x,y nodal values z to point xp,yp
! result is in zp
!
! needs no other vectors but needs x,y of nodes

        double precision x(3),y(3)  !coordinates of nodes
        double precision z(3)       !values on nodes
        double precision xp,yp      !coordinates of point
        double precision zp         !interpolated value (return)

	double precision eps
	parameter ( eps = 1.d-14 )

        integer i,ii,iii
        double precision zh,f,fh
        double precision a(3),b(3),c(3)

        f = 0.
        do i=1,3
           ii=mod(i,3)+1
           iii=mod(ii,3)+1
           a(i)=x(ii)*y(iii)-x(iii)*y(ii)
           b(i)=y(ii)-y(iii)
           c(i)=x(iii)-x(ii)
           f=f+a(i)
        end do
	f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	if( f .le. eps ) goto 99

        zh=0.
        do i=1,3
           fh=(a(i)+xp*b(i)+yp*c(i))/f
           zh=zh+z(i)*fh
        end do

        zp = zh

	return
   99	continue
	write(6,*) 0,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop elemintp: area of element'
        end

!******************************************************

	subroutine find_close_elem(ieold,xp,yp,ielem)

! finds element for point (xp,yp) starting from ieold
!
! uses data structure ev and ieltv

        use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	double precision xp,yp
	integer ielem	!element number on return

	logical binit,bdebug
	integer ie,ii,iside,lmax,loop
	double precision xi,ximin
	double precision a(3),b(3),c(3)

	bdebug = .false.
	lmax = 10

	call is_init_ev(binit)

!-------------------------------------------------------------
! check if old element is given -> if not test all elements
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (1) ',ieold
	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

	if( bdebug ) write(6,*) 'ggu_xi (2) ',ieold
	if( .not. binit ) then
	  call find_elem_from_old(ieold,xp,yp,ielem)
	  return
	end if

!-------------------------------------------------------------
! start from old element
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (3) ',ieold

	loop = 0
	ie = ieold
	do while( ie .gt. 0 )
	  iside = 0
	  ximin = 1.e+30
	  call xi_abc(ie,a,b,c)
	  do ii=1,3
	    xi = a(ii) + b(ii)*xp + c(ii)*yp
	    if( bdebug ) write(6,*) 'ggu_xiii ',ie,ii,xi
	    if( xi .lt. ximin ) then
	      ximin = xi
	      iside = ii
	    end if
	  end do
	  if( ximin .ge. 0. ) then
	    ielem = ie
	    return
	  end if
	  if( iside .le. 0 .or. iside .gt. 3 ) then
	    if( bdebug ) write(6,*) '******** ',iside,ie,ximin,xi
	    ie = 0
	  else
	    ie = ieltv(iside,ie)
	    if( bdebug ) write(6,*) 'ggu_xiii iterate',ie,iside,ximin
	  end if
	  loop = loop + 1
	  if( loop .gt. lmax ) ie = 0
	end do

	ielem = 0

	end

!******************************************************

	subroutine find_elem_from_old(ieold,xp,yp,ielem)

! finds element for point (xp,yp) starting from ieold

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	double precision xp,yp
	integer ielem	!element number on return

	integer iem,iep

!-------------------------------------------------------------
! check if old element is given -> if not test all elements
!-------------------------------------------------------------

	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

!-------------------------------------------------------------
! check if in old element
!-------------------------------------------------------------

	if( in_element(ieold,xp,yp) ) then
	  ielem = ieold
	  return
	end if

!-------------------------------------------------------------
! start from old element going upwards and downwards
!-------------------------------------------------------------

	iem = ieold-1
	if( iem .lt. 1 ) iem = nel		!BUG_27.01.2011
	iep = ieold+1
	if( iep .gt. nel ) iep = 1		!BUG_27.01.2011

	do while( iem .ne. ieold .and. iep .ne. ieold )
	  if( in_element(iem,xp,yp) ) then
	    ielem = iem
	    return
	  end if
	  iem = iem - 1
	  if( iem .lt. 1 ) iem = nel

	  if( in_element(iep,xp,yp) ) then
	    ielem = iep
	    return
	  end if
	  iep = iep + 1
	  if( iep .gt. nel ) iep = 1
	end do

!-------------------------------------------------------------
! no element found
!-------------------------------------------------------------

	ielem = 0

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************

	subroutine find_element(xp,yp,ielem)

! finds element for point (xp,yp)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision xp,yp
	integer ielem	!element number on return


	integer ie

	do ie=1,nel
	  if( in_element(ie,xp,yp) ) then
		  ielem = ie
		  return
	  end if
	end do

	ielem = 0

	end

!******************************************************

	function in_element(ie,xp,yp)

! checks if point (xp,yp) is in element ie

	implicit none

	logical in_element
	integer ie
	double precision xp,yp

	integer ii,k,in
	double precision xmin,ymin,xmax,ymax
	double precision x(3),y(3)

	in_element = .false.

	do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	end do

	xmin = min(x(1),x(2),x(3))
	ymin = min(y(1),y(2),y(3))
	xmax = max(x(1),x(2),x(3))
	ymax = max(y(1),y(2),y(3))

	if( xp .ge. xmin .and. xp .le. xmax ) then
	  if( yp .ge. ymin .and. yp .le. ymax ) then
		in = intri(x,y,xp,yp)
		if( in .gt. 0 ) in_element = .true.
	  end if
	end if

	end

!******************************************************

	subroutine get_xy_elem(ie,x,y)

! returns x,y of vertices of element ie

	implicit none

	integer ie
	double precision x(3), y(3)

	integer ii,k

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

!******************************************************

	subroutine get_scal_elem(ie,sv,s)

! returns s at vertices of element ie

	implicit none

	integer ie
	double precision sv(nkn)
	double precision s(3)

	integer ii,k

	do ii=1,3
	  k = nen3v(ii,ie)
	  s(ii) = sv(k)
	end do

	end

!******************************************************

!==================================================================
        end module regular_2nc
!==================================================================
