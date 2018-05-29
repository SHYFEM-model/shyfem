c
c $Id: subreg.f,v 1.17 2009-09-14 08:20:58 georg Exp $
c
c routines for interpolation onto regular grid
c
c contents :
c
c subroutine setgeo(x0,y0,dx,dy,flag)
c		sets grid values for interpolation
c subroutine getgeo(x0,y0,dx,dy,flag)
c		gets grid values for interpolation
c subroutine getgeoflag(flag)
c		gets flag value for interpolation
c
c subroutine av2am(av,am,ip,jp)
c		interpolation of av onto a regular net
c subroutine av2amk(bwater,av,am,ip,jp)
c		interpolation of av onto a regular net
c function intri(x,y,xp,yp)
c		point in triangle or not
c function intrid(x,y,xp,yp)
c		point in triangle or not (double precision version)
c subroutine am2av(am,av,ip,jp)
c		interpolation of am onto finite element mesh
c function am2val(am,ip,jp,xx,yy)
c		interpolation of am onto finite element mesh
c subroutine ave2am(av,am,ip,jp)
c		interpolation of av (elementwise) onto a regular net
c
c subroutine mkmask(bwater,zv,href,hzoff)
c		makes mask in element
c
c subroutine mimareg(am,ip,jp,amin,amax)
c		computes min/max of regular matrix (without flag values)
c subroutine a2char(am,ac,ip,jp)
c		creates 1 char representation of matrix
c subroutine prchar(ac,ip,jp)
c		prints 1 char representation of matrix
c
c subroutine femintp(ie,z,xp,yp,zp)
c               interpolation in element (with ev)
c subroutine elemintp(x,y,z,xp,yp,zp)
c               interpolation in element (no ev)
c
c subroutine find_elem_from_old(ieold,xp,yp,ielem)
c		finds element for point (xp,yp) starting from ieold
c subroutine find_element(xp,yp,ielem)
c		finds element for point (xp,yp)
c function in_element(ie,xp,yp)
c		checks if point (xp,yp) is in element ie
c subroutine get_xy_elem(ie,x,y)
c		returns x,y of vertices of element ie
c
c revision log :
c
c 18.11.1998	ggu	routine commented
c 18.11.1998	ggu	routine setgeo introduced
c 19.11.1998	ggu	routines a2char, prchar added
c 19.10.1999	ggu	routine mkmask added from subutl
c 25.11.2004	ggu	new routines femintp and elemintp for interpolation
c 14.03.2005	ggu	new routines for interpolation in element
c 11.03.2009	ggu	new helper routine getgeoflag()
c 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
c 26.01.2011	ggu&mb	handling extrapolation in am2av()
c 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
c 31.03.2011	ggu	new routine elemmask()
c 24.11.2011	ggu	new routine find_close_elem()
c 20.06.2012	ggu	new routine get_scal_elem()
c 07.10.2012	ggu	new routine av2fm()
c 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
c 26.10.2012	ggu	bug fix: do not access not existing storage
c 30.05.2014	ggu	in av2amk() do not interpolate for flag values
c 07.07.2014	ggu	new routine intp_reg()
c 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
c 05.05.2016	ggu	file restructured (module)
c 14.05.2016	ggu	allow for extension of grid -> bregextend
c 23.06.2016	ggu	allow for eps in computing box
c 23.09.2016	ggu	allow for eps in computing box and reg intp
c 23.04.2017	ggu	new routine intp_reg_single_nodes()
c 23.05.2017	ggu	file split into subreg, submask and subfind
c 04.12.2017	ggu	check t,u values and correct if out of bounds
c 18.05.2018	ggu	more checks on fr routines, introduced ierr, bextend
c 26.05.2018	ggu	even more checks and debug output
c
c notes :
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
c******************************************************

!==================================================================
        module regular
!==================================================================

	implicit none

	logical, save :: bregextend = .false.

	real, save :: pxareg = 0.	!x coordinate of lower,left point (x0)
	real, save :: pyareg = 0.	!y coordinate of lower,left point (y0)
	real, save :: pxdreg = 0.	!grid spacing in x direction (dx)
	real, save :: pydreg = 0.	!grid spacing in y direction (dy)
	real, save :: pzlreg = -999.	!flag for land points

!==================================================================
        end module regular
!==================================================================

	subroutine setgeo(x0,y0,dx,dy,flag)

c sets grid values for interpolation
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	use regular

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

	pxareg = x0
	pyareg = y0
	pxdreg = dx
	pydreg = dy
	pzlreg = flag

	end

c******************************************************

	subroutine getgeo(x0,y0,dx,dy,flag)

c gets grid values for interpolation
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	use regular

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

	x0   = pxareg
	y0   = pyareg
	dx   = pxdreg
	dy   = pydreg
	flag = pzlreg

	end

c******************************************************

	subroutine getgeoflag(flag)

c gets flag value for interpolation
c
c flag                value for land points

	use regular

	implicit none

	real flag	!flag for land points

	flag = pzlreg

	end

c******************************************************

	subroutine setregextend(bextend)

c sets flag to decide if extend interpolated grid

	use regular

	implicit none

	logical bextend

	bregextend = bextend

	end

c******************************************************

	subroutine getregextend(bextend)

c gets flag to decide if extend interpolated grid

	use regular

	implicit none

	logical bextend

	bextend = bregextend

	end

c******************************************************

	subroutine setreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy
	real flag

	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x0
	regpar(4) = y0
	regpar(5) = dx
	regpar(6) = dy
	regpar(7) = flag

	end

c******************************************************

	subroutine getreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy
	real flag

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	end

c******************************************************

	subroutine printreg(regpar)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy,x1,y1
	real flag

        call getreg(regpar,nx,ny,x0,y0,dx,dy,flag)

        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy
        write(6,'(4x,a,2i12)') 'nx,ny: ',nx,ny
        write(6,'(4x,a,2f12.4)') 'x0,y0: ',x0,y0
        write(6,'(4x,a,2f12.4)') 'x1,y1: ',x1,y1
        write(6,'(4x,a,2f12.4)') 'dx,dy: ',dx,dy
        write(6,'(4x,a,2f12.4)') 'flag : ',flag

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine find_position_to_coord(x,y,ix,iy)

! finds closest position (ix,iy) to coordinate (x,y) in reg grid

	implicit none

	real x,y
	integer ix,iy

	real x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	ix = nint( (x-x0)/dx + 1. )
	iy = nint( (y-y0)/dy + 1. )

	end

c******************************************************

	subroutine find_coord_to_position(ix,iy,x,y)

! finds coordinate (x,y) to given position (ix,iy) in reg grid

	implicit none

	integer ix,iy
	real x,y

	real x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	x = x0 + (ix-1)*dx
	y = y0 + (iy-1)*dy

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine av2am(av,am,ip,jp)

c interpolation of av onto a regular net (nodal values)
c
c av                    array to be interpolated
c am                    matrices of interpolated values (u,v,z,h)
c ip,jp                 dimension of matrices
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	use basin

	implicit none

	integer ip,jp
	real av(nkn)
	real am(ip,jp)

	logical bwater(nel)

	bwater = .true.

	call av2amk(bwater,av,am,ip,jp)

	end

c******************************************************

	subroutine av2amk(bwater,av,am,ip,jp)

c interpolation of av onto a regular net (nodal values) with mask
c
c bwater		mask for water points
c av                    array to be interpolated
c am                    matrices of interpolated values (u,v,z,h)
c ip,jp                 dimension of matrices
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	use basin

	implicit none

c arguments
	integer ip,jp
	real av(nkn)
	real am(ip,jp)
	logical bwater(nel)
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c local
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	integer iflag
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
c function
	integer intrid

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

c************************************************
c************************************************
c************************************************

	subroutine av2fm(fm,ip,jp)

c computation of interpolation matrix (nodal values to regular grid) with mask

	use basin

	implicit none

	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices

	logical bwater(nel)	!wet mask for each element

	bwater = .true.

	call av2fmk(bwater,fm,ip,jp)

	end

c************************************************

	subroutine av2fmk(bwater,fm,ip,jp)

c computation of interpolation matrix (nodal values to regular grid) with mask
c
c the interpolation can be carried out as
c
c	do j=1,jp
c	  do i=1,ip
c	    ie = nint(fm(4,i,j))
c	    if( ie .gt. 0 ) then
c	      a = 0.
c	      do ii=1,3
c	        k = nen3v(ii,ie)
c		a = a + val(k) * fm(ii,i,j)
c	      end do
c	    else
c	      a = flag
c	    end if
c	    am(i,j) = a
c	  end do
c	end do
	        
	use basin

	implicit none

c arguments
	logical bwater(nel)	!wet mask for each element
	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c local
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	logical bok
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
c function
	integer intrid

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

c************************************************

        subroutine fm2am2d(femval,nx,ny,fm,am)

c interpolation 2d of fem values to regular grid using fm matrix

	use basin

        implicit none

        real femval(nkn)		!values of fem array
        integer nx,ny			!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nx,ny)			!interpolated values (return)

	integer nlvdi,nlv
	integer ilhv(nel)

	nlvdi = 1
	nlv = 1
	ilhv = 1

        call fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

	end

c************************************************

        subroutine fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

c interpolation 3d of fem values to regular grid using fm matrix

	use basin

        implicit none

	integer nlvdi			!vertical dimension of fem array
        integer ilhv(nel)		!vertical discretization (element!!)
        real femval(nlvdi,nkn)		!values of fem array
        integer nlv,nx,ny		!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nlv,nx,ny)		!interpolated values (return)

	logical bflag
        integer i,j,l,lmax,ie,ii,k
        real a
        real flag

	call getgeoflag(flag)

        do j=1,ny
          do i=1,nx
            ie = nint(fm(4,i,j))
            lmax = 0
            if( ie .gt. 0 ) lmax = ilhv(ie)
	    lmax = min(lmax,nlv)
            do l=1,lmax
              a = 0.
	      bflag = .false.
              do ii=1,3
                k = nen3v(ii,ie)
                a = a + femval(l,k) * fm(ii,i,j)
		if( femval(l,k) == flag ) bflag = .true.
              end do
	      if( bflag ) a = flag
              am(l,i,j) = a
            end do
            do l=lmax+1,nlv
              am(l,i,j) = flag
            end do
          end do
        end do

        end

c************************************************
c************************************************
c************************************************

	subroutine fm_extra_setup(nx,ny,fmextra)

! sets up fmextra structure to allow interpolation from fem nodes to reg grid

	use basin

	implicit none

        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)

	logical bout,berror
	integer ix,iy,iix,iiy
	integer j,k
	real x0,y0,dx,dy
	real xmax,ymax
	real x,y,x1,y1
	real eps,flag,tueps
	real t,u
	double precision d,w,d2

	integer, save :: jx(4) = (/0,1,1,0/)
	integer, save :: jy(4) = (/0,0,1,1/)
	double precision fmweight(nx,ny)

	eps = 0.01
	tueps = 0.0001
	fmextra = 0.
	fmweight = 0.

	call getgeo(x0,y0,dx,dy,flag)
	xmax = x0 + (nx-1)*dx
	ymax = y0 + (ny-1)*dy

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
	  if( t-1. > tueps .or. t < -tueps ) bout = .true.
	  if( u-1. > tueps .or. u < -tueps ) bout = .true.
	  !if( t.gt.1. .or. t.lt.0. ) bout = .true.
	  !if( u.gt.1. .or. u.lt.0. ) bout = .true.
	  if( bout ) then
	    write(6,*) 'out of domain: '
	    write(6,*) '... ',k,nx,ny,dx,dy
	    write(6,*) '... ',x0,y0,xmax,ymax
	    write(6,*) '... ',x1,y1,x,y
	    write(6,*) '... ',ix,iy,t,u
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

c************************************************

	subroutine fm_extra_3d(nlvdi,nlv,il,nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	use basin

	implicit none

	integer nlvdi,nlv
	integer il(nkn)
        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)
	real femdata(nlvdi,nkn)
	real regdata(nlvdi,nx,ny)

	integer k,l,lmax
	real flag
	real fem2d(nkn)
	real reg2d(nx,ny)

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

c************************************************

	subroutine fm_extra_2d(nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	use basin

	implicit none

        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)
	real femdata(nkn)
	real regdata(nx,ny)

	integer ix,iy,iix,iiy
	integer j,k
	real x0,y0,dx,dy
	real eps,flag
	real regval,femval
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

c************************************************
c************************************************
c************************************************

	function intri(x,y,xp,yp)

c point in triangle or not
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)

	implicit none

c arguments
	integer intri
	real x(3),y(3),xp,yp
c local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
c save
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

c************************************************

	function intrid(x,y,xp,yp)

c point in triangle or not (double precision version)
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)

	implicit none

c arguments
	integer intrid
	double precision x(3),y(3),xp,yp
c local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
c save
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

c****************************************************************
c
	function intri0(x,y,xp,yp)
c
c point in triangle or not
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)
c
	implicit none
c
c arguments
	integer intri0
	real x(3),y(3),xp,yp
c local
	integer i,k1,k2
	double precision xs,ys,x12,y12
	double precision det,detlam,rlamb
	double precision eps
c save
	save eps
	data eps /1.e-13/
c
	xs=0.
	ys=0.
	do i=1,3
	   xs=xs+x(i)
	   ys=ys+y(i)
	end do
	xs=xs/3.
	ys=ys/3.
c
	intri0=0
c
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
c
	intri0=1	!inside
c
	return
	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine intp_reg_nodes(nx,ny,x0,y0,dx,dy,flag,regval
     +				,femval,ierr)

! interpolates regular grid to FEM grid - values are on nodes

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	real femval(nkn)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval
     +				,nkn,xgv,ygv,femval,ierr)

	end

c****************************************************************

	subroutine intp_reg_elems(nx,ny,x0,y0,dx,dy,flag,regval
     +				,femval,ierr)

! interpolates regular grid to FEM grid - values are on elements

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	real femval(nel)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	real xp(nel),yp(nel)

	call bas_get_elem_coordinates(xp,yp)

	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval
     +				,nel,xp,yp,femval,ierr)

	end

c****************************************************************

	subroutine intp_reg_single_nodes(nx,ny,x0,y0,dx,dy,flag
     +				,regval,np,nodes
     +				,femval,ierr)

! interpolates regular grid to single nodes

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	integer np		!total number of nodes
	integer nodes(np)	!node numbers
	real femval(np)		!interpolated values on nodes (return)
	integer ierr		!error code (return)

	real xp(np),yp(np)

	call bas_get_special_coordinates(np,nodes,xp,yp)

	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval
     +				,np,xp,yp,femval,ierr)

	end

c****************************************************************

	subroutine intp_reg(nx,ny,x0,y0,dx,dy,flag,regval
     +				,np,xp,yp,femval,ierr)

c interpolation of regular array onto fem grid - general routine
c
c ierr:
c		= 0	no errors
c		< 0	interpolation out of domain (extrapolation)
c		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	integer np		!number of fem points
	real xp(np)
	real yp(np)
	real femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	logical bextra,bout,bflag
	logical bintpout,bintpflag,bextend
	logical bdebug
	integer k
	integer imin,jmin
	integer iflag,iout
	real xx,yy,z1,z2,z3,z4,x1,y1,t,u
	real xn,yn
	real zz(4)
 
	real, parameter :: eps = 3.e-4
	!real, parameter :: eps = 0.
	logical outbox
	outbox(t) = ( t-1. > eps .or. t < -eps )

	!bintpout = .false.	!interpolate even if outside
	!bintpout = .true.	!interpolate even if outside
	!bintpflag = .false.	!interpolate even if flag
	!bintpflag = .true.	!interpolate even if flag
	bdebug = .false.	

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
		if( bdebug ) then
		  write(6,*) 'reg intp: ',t,u
	          call reg_debug_1('out',k,xx,yy)
		end if
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
	      if( bdebug ) call reg_debug_1('flag',k,xx,yy)
	      cycle
	    end if

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout .gt. 0 ) ierr = - iout - iflag
	if( iflag .gt. 0 ) ierr = iflag

	!write(6,*) 'intp_reg: ierr = ',ierr,iout,iflag

	!if( ierr /= 0 ) then
	!  write(6,*) 'reg: ',ierr,np
	!  write(6,*) femval
	!  stop
	!end if

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

c****************************************************************

	subroutine reg_debug_1(what,k,x,y)

	implicit none

	character*(*) what
	integer k
	real x,y

	write(166,*) trim(what),k,x,y

	end

c****************************************************************

	subroutine recover_flag(zz,z1,z2,z3,z4,flag)

	implicit none

	real zz(4)
	real z1,z2,z3,z4
	real flag

	integer i,ic
	real zt

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

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine intp_reg_setup_fr(nx,ny,x0,y0,dx,dy,np,xp,yp,fr,ierr)

c interpolation of regular array onto fem grid - general routine
c
c produces array fr that can be used to interpolate

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	integer np		!total number of sparse points
	real xp(np),yp(np)	!coordinates of sparse grid
	real fr(4,np)		!interpolation array on sparse grid (return)
	integer ierr		!number of points with no interpolation

	logical bextend,bout,bdebug
	integer k
	integer imin,jmin
	real xx,yy,x1,y1,t,u
	real xn,yn
	real, parameter :: eps = 1.e-4
	real, parameter :: fr_flag = -999.
 
	bdebug = ( ierr /= 0 )
	ierr = 0
	imin = 0
	jmin = 0
	fr = 0.

	call getregextend(bextend)

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

	    bout = .false.
	    if( t < 0. ) then
	      if( t >= -eps ) t = 0.
	      if( t >= -1. .and. bextend ) t = 0.
	      if( t < 0. ) bout = .true.
	    end if
	    if( u < 0. ) then
	      if( u >= -eps ) u = 0.
	      if( u >= -1. .and. bextend ) u = 0.
	      if( u < 0. ) bout = .true.
	    end if
	    if( t > 1. ) then
	      if( t <= 1.+eps ) t = 1.
	      if( t <= 2. .and. bextend ) t = 1.
	      if( t > 1. ) bout = .true.
	    end if
	    if( u > 1. ) then
	      if( u <= 1.+eps ) u = 1.
	      if( u <= 2. .and. bextend ) u = 1.
	      if( u > 1. ) bout = .true.
	    end if
	    if( bout .and. bdebug ) then
	      call reg_debug_1('out',k,xx,yy)
	      t = fr_flag
	      u = fr_flag
	      ierr = ierr + 1
	    end if

	    fr(1,k) = imin
	    fr(2,k) = jmin
	    fr(3,k) = t
	    fr(4,k) = u
	end do
 
	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg_setup_fr: internal error (1)'
	end

c****************************************************************

	subroutine intp_reg_intp_fr(nx,ny,flag,regval
     +				,np,fr,femval,ierr)

c interpolation of regular array onto fem grid - general routine
c
c ierr:
c		= 0	no errors
c		< 0	interpolation out of domain (extrapolation)
c		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	real flag
	real regval(nx,ny)
	integer np		!number of fem points
	real fr(4,np)
	real femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	logical bdebug,bflag,bextend,bout
	integer k,ks
	integer imin,jmin
	integer iflag_tot,iout_tot
	real, parameter :: fr_flag = -999.
	real z1,z2,z3,z4,t,u
	real zz(4)
 
	bdebug = .true.
	bdebug = ( ierr /= 0 )

        call getregextend(bextend)

	iflag_tot = 0	!used flag for interpolation
	iout_tot = 0	!used outside point for interpolation

	do k=1,np
 
	    femval(k) = flag
 
	    imin = nint(fr(1,k))
	    jmin = nint(fr(2,k))
	    t = fr(3,k)
	    u = fr(4,k)

	    bout = ( t == fr_flag )	!pre-computed
	    if( bout ) then
	      if( bdebug ) then
	        write(6,*) 'debug intp_reg_intp_fr: ',k,u,t
	        call reg_debug_1('out',k,0.,0.)
	      end if
	      iout_tot = iout_tot + 1
	      cycle
	    end if

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)
            zz = (/z1,z2,z3,z4/)

            if( any(zz == flag) .and. bextend ) then
              call recover_flag(zz,z1,z2,z3,z4,flag)
            end if

	    bflag = any(zz == flag)
            if( bflag .and. bdebug ) then
	      call reg_debug_1('flag',k,real(ierr),0.)
	      iflag_tot = iflag_tot + 1
              cycle
            end if

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout_tot .gt. 0 ) ierr = - iout_tot
	if( iflag_tot .gt. 0 ) ierr = iflag_tot

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg_intp_fr: internal error (1)'
	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine am2av(am,av,ip,jp)

c compatibility for old calls - from regular to fem

	use basin

	implicit none

	integer ip,jp
	real av(nkn)
	real am(ip,jp)

	integer ierr
	real pxareg,pyareg,pxdreg,pydreg,pzlreg

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	call intp_reg(ip,jp,pxareg,pyareg,pxdreg,pydreg,pzlreg,am
     +				,nkn,xgv,ygv,av,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'interpolation am2av: ierr = ',ierr
	  !stop 'error stop am2av: error in interpolation'
	end if

	end

c******************************************************
c******************************************************
c******************************************************

        subroutine elemintp(x,y,z,xp,yp,zp)

c interpolation in element (no ev)
c
c interpolates in element given by x,y nodal values z to point xp,yp
c result is in zp
c
c needs no other vectors but needs x,y of nodes

        real x(3),y(3)  !coordinates of nodes
        real z(3)       !values on nodes
        real xp,yp      !coordinates of point
        real zp         !interpolated value (return)

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
        f = c(3)*b(2) - c(2)*b(3)               !bug_f_64bit
        if( f .le. eps ) goto 99

        zh=0.
        do i=1,3
           fh=(a(i)+xp*b(i)+yp*c(i))/f
           zh=zh+z(i)*fh
        end do

        zp = zh

        return
   99   continue
        write(6,*) 0,f
        write(6,*) x
        write(6,*) y
        write(6,*) a
        write(6,*) b
        write(6,*) c
        stop 'error stop elemintp: area of element'
        end

c******************************************************
c******************************************************
c******************************************************

