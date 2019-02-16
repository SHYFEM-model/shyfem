
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
c 30.10.2003	ggu	subroutine prepare_bc included in this file
c 04.03.2004	ggu	writes also number of variables (1)
c 27.01.2009	ggu	new routine scal_intp from laplap
c
c notes :
c
c please prepare file input.dat like this:
c
c----------------- start
c k1	val1
c k2	val2
c ...
c kn	valn
c----------------- end
c
c or
c
c----------------- start
c x1 y1	val1
c x2 y2	val2
c ...
c xn yn	valn
c----------------- end
c
c run memory and set the basin ( memory -b venlag62 )
c run laplap with input file ( laplap < input.dat )
c
c****************************************************************

        program scal_intp

c optimal interpolation

	use evgeom
	use basin

	implicit none

	include 'param.h'

	integer ndim
	parameter (ndim = 100)

	real, allocatable :: zv(:)

	real xp(ndim),yp(ndim),zp(ndim)

	integer k,np,mode
        integer it,ilev,ivar
	real afact
	real flag
	real zmin,zmax

	integer iapini

	flag = 1.23456e+23

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,0,0,0) .le. 0 ) stop

	call bas_info

c-----------------------------------------------------------------
c set up ev
c-----------------------------------------------------------------

	call set_ev
	call check_ev

	allocate(zv(nkn))

c-----------------------------------------------------------------
c read BC and interpolate
c-----------------------------------------------------------------

	call get_input(mode,afact)
	call read_data(mode,ndim,np,xp,yp,zp)

	call interpol_e(afact,np,xp,yp,zp,zv)

c-----------------------------------------------------------------
c min/max of interpolated values
c-----------------------------------------------------------------

	call mima(zv,nkn,zmin,zmax)
	write(6,*) 'min/max: ',zmin,zmax

c-----------------------------------------------------------------
c write to NOS file laplace.nos
c-----------------------------------------------------------------

	call wrnos2d('scalar','laplace interpolation',zv)

c-----------------------------------------------------------------
c write to DAT file laplace.dat
c-----------------------------------------------------------------

	it = 0
        ilev = 0
        ivar = 1

	open(1,file='scalar.dat',status='unknown',form='unformatted')
	write(1) it,nkn,ilev,ivar
	write(1) (zv(k),k=1,nkn)
	close(1)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine get_input(mode,afact)

	implicit none

	integer mode
	real afact

        write(6,*)
        write(6,*) 'What type of input file?'
        write(6,*)
        write(6,*) '1   values are given on nodes'
        write(6,*) '2   values are given with coordinates'
        write(6,*)
        write(6,*) 'In the first case the lines of the input file'
	write(6,*) 'have to be in the format'
        write(6,*) '     node value'
        write(6,*) 'In the second case the lines of the input file'
	write(6,*) 'have to be in the format'
        write(6,*) '     x y value'
        write(6,*)
        write(6,*) 'Enter choice: '
        read(5,'(i10)') mode
        write(6,*) 'Your choice is : ', mode

        if( mode .lt. 1 .or. mode .gt. 2 ) stop

        write(6,*)
        write(6,*) 'Enter interpolation factor:'
        write(6,*)
        write(6,*) 'This factor multiplies the typical length scale'
        write(6,*) 'to provide a different size of influence.'
        write(6,*) 'Larger gives smoother scalar fields.'
        write(6,*) 'Default is 1.'
        write(6,*)
        write(6,*) 'Enter choice: '
        read(5,'(f10.0)') afact
        write(6,*) 'Your choice is : ', afact

        if( afact .le. 0. ) afact = 1.

	end

c******************************************************************

	subroutine read_data(mode,ndim,np,xp,yp,zp)

	use basin

	implicit none

	include 'param.h'

	integer mode
	integer ndim,np
	real xp(1),yp(1),zp(1)

	character*80 file
	integer n,k,ki
	real x,y,val

	integer ipint

	n = 0
	file = 'input.dat'

        write(6,*) 'reading file : ',file
        open(1,err=99,file=file,status='old',form='formatted')

    1   continue
	  if( mode .eq. 1 ) then
            read(1,*,end=2) k,val
	    ki = ipint(k)
	    call coord(ki,x,y)
	  else
            read(1,*,end=2) x,y,val
	  end if
          n = n + 1
          if( n .gt. ndim ) stop 'error stop: ndim'
          xp(n) = x
          yp(n) = y
          zp(n) = val
          goto 1
    2   continue

        close(1)
        np = n
        write(6,*) np,' data points read'

        return
   99   continue
        write(6,*) 'Cannot read file: ',file
        stop 'error stop read_data: no such file'
	end

c******************************************************************

	subroutine interpol_e(afact,np,xp,yp,zp,zv)

c interpolates depth values
c
c exponential interpolation with max radius

	use evgeom
	use basin

	implicit none

	real afact
	integer np
	real xp(1),yp(1),zp(1)
	real zv(1)

	include 'param.h'

	integer ie,ii,k
	integer nktot
	real x,y,d
	real hmin
	real depth
	integer iaux,inum,ityp
	integer n,i
	integer ihmed
        integer nintp
	real xmin,ymin,xmax,ymax
	real weight,r2,w
        real r2max,sigma2
	real hmed
	real fact,dx,dy
	real area,x0,y0
	double precision atot
	logical ok(nkn)

        write(6,*) 'starting exponential interpolation with max radius'

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------

	atot = 0.
	do ie=1,nel
	  area = ev(10,ie)
	  atot = atot + area
	end do

	!afact = 1.
	area = atot / np
	area = afact * afact * area

	write(6,*) 'total area: ',atot,area,sqrt(area)

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

	fact = 2.
	nktot = 0

	do while( nktot .lt. nkn )

	nktot = 0

	do k = 1,nkn

	  x0 = xgv(k)
	  y0 = ygv(k)

          r2max = fact*area     !maximum radius to look for points
          sigma2 = fact*area    !standard deviation also grows

	  if( ok(k) ) then
	    nktot = nktot + 1
	  else
	    depth = 0.
	    weight = 0.
            nintp = 0
	    do n=1,np
	      x = xp(n)
	      y = yp(n)
	      d = zp(n)
	      r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)
	      !if( r2 .le. r2max ) then
	        w = exp(-r2/sigma2)
	        depth = depth + d * w
	        weight = weight + w
                nintp = nintp + 1
	      !end if
	    end do
	    if( nintp .gt. 0 ) then
	      if( weight .le. 0. ) then
                write(6,*) nintp,weight,r2max,ie
                stop 'error stop interpole: zero weight from points'
              end if
	      nktot = nktot + 1
	      ok(k) = .true.
	      depth = depth / weight
	      zv(k) = depth
	    end if
	  end if
	end do

	write(6,*) 'Nodes without value : ',fact,nkn-nktot,nkn

	fact = fact * 2.

	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

