
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999-2001,2003-2005,2009-2020  Georg Umgiesser
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

! routines for rounding
!
! function roundm(r,mode)		rounds r to the closest value
! function rround(r,rmaster,mode)	rounds r to next rmaster value
! function rdist(xmin,ymin,xmax,ymax)	computes gridspacing 
! function divdist(x,n,mode)		divides x into n equal pieces
!
! function rnext(r,mode)		finds closest value to r
! function rnexta(r,mode)		finds closest value to r (absolute)
! function rnextsub(r)			finds best subdivision for value r
!
! subroutine logvals(amin,amax,idiv,ntk,rval,aval)	computes log scale vals
!
! revision log :
!
! 16.02.1999	ggu	bpixel: write bounding box to ps file (as comment)
! 09.02.2000	ggu	use inboxdim to compute box to plot
! 12.06.2000	ggu	get gray value for basin mode 3 from str file
! 11.02.2001	ggu	routine to compute typical length scale
! 21.08.2003	ggu	occupy is called elsewhere
! 16.12.2004	ggu	changed reggrid to plot regular grid
! 02.03.2005	ggu	in bash: bug fix -> get size of grid if not given
! 12.06.2009	ggu	new routines to handle spherical coords. & regular grid
! 15.06.2009	ggu	call to reggrid() changeed -> pass in gray value
! 14.09.2009	ggu	new routine divdist()
! 22.02.2010	ggu	new routine bw_frame() to plot bw scale around plot
! 23.03.2010	ggu	changed v6.1.1
! 09.04.2010	ggu	bug fix in frac_pos() -> maybe compiler error
! 20.12.2010	ggu	changed VERS_6_1_16
! 17.05.2011	ggu	new routine basin_number()
! 31.05.2011	ggu	changed VERS_6_1_23
! 30.03.2012	ggu	changed VERS_6_1_51
! 30.08.2012	ggu	new routines to automatically label spherical grid
! 12.09.2012	ggu	changed VERS_6_1_57
! 24.10.2012	ggu	bug in labelling non spherical grid (returned -1)
! 02.05.2013	ggu	handle fact in spherical coords
! 02.05.2013	ggu	meteo point plotting (plot_meteo_points())
! 13.06.2013	ggu	bug fix in spherical_fact() -> set fact to 1
! 19.06.2013	ggu	changed VERS_6_1_66
! 13.12.2013	ggu	new mode=4 for plotting gray grid over scalar variable
! 28.01.2014	ggu	changed VERS_6_1_71
! 30.05.2014	ggu	new metpnt for meteo points, imicro computed
! 18.07.2014	ggu	changed VERS_7_0_1
! 13.10.2014	ggu	changed VERS_7_0_2
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.12.2014	ggu	changed VERS_7_0_8
! 23.12.2014	ggu	changed VERS_7_0_11
! 15.01.2015	ggu	changed VERS_7_1_1
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.02.2015	ggu	also plot other points, also regular points
! 26.02.2015	ggu	changed VERS_7_1_5
! 05.05.2015	ggu	changed VERS_7_1_10
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 12.10.2015	ggu	fix in rround() to handle rmaster==0 case
! 19.02.2016	ggu	changed VERS_7_5_2
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 30.09.2016	ggu	changed VERS_7_5_18
! 25.05.2017	ggu	changed VERS_7_5_28
! 11.07.2017	ggu	changed VERS_7_5_30
! 26.09.2017	ggu	changed VERS_7_5_32
! 14.11.2017	ggu	changed VERS_7_5_36
! 22.02.2018	ggu	changed VERS_7_5_42
! 19.04.2018	ggu	changed VERS_7_5_45
! 06.07.2018	ggu	changed VERS_7_5_48
! 18.12.2018	ggu	changed VERS_7_5_52
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 13.02.2020	ggu	rounding routines in this extra file copied
!
!*************************************************************

	function roundm(r,mode)

! rounds r to the closest given value
!
! r		value to round
! mode
!		 0: do not round, return r
!		 1: round to higher value
!		-1: round to lower value

	implicit none

	real roundm
	real r
	integer mode

	double precision raux,fact,sign,rr
	integer i

	integer ndim
	parameter (ndim = 4)
	real rmaster(ndim),eps
	save rmaster,eps
	data eps /1.e-5/
	data rmaster /1.,2.,5.,10./
!	data rmaster /1.,2.,4.,5.,8.,10./

	roundm = r
	raux = 0
	if( mode .eq. 0 ) return

	fact = 1.0d+0
	if( r .lt. 0. ) then
	  sign = -1.
	  rr = -r
	else
	  sign = 1.
	  rr = r
	end if

	do while( rr*fact .gt. 10. )
	  fact = 0.1d+0 * fact
	end do
	do while( rr*fact .lt. 1. )
	  fact = 10.d+0 * fact
	end do

	rr = rr * fact		!rr is between 1. and 10.

	if( mode .gt. 0 ) then
	  do i=ndim,1,-1
	    if( rr .le. rmaster(i) + eps ) raux = rmaster(i)
	  end do
	else if( mode .lt. 0 ) then
	  do i=1,ndim
	    if( rr .ge. rmaster(i) - eps ) raux = rmaster(i)
	  end do
	end if

	roundm = sign * raux / fact

	end

!****************************************************

	function rround(r,rmaster,mode)

! rounds r to next rmaster value
!
! r		value to round
! rmaster	value to which r is rounded (must be positive)
! mode
!		 0: do not round, return r
!		 1: round to higher value
!		-1: round to lower value
!
! negative values are respected

	implicit none

	real rround
	real r,rmaster
	integer mode

	integer iaux
	double precision d,dmaster,daux

	d = r
	dmaster = rmaster

	if( dmaster == 0. ) then	!FIXME
	  rround = d
	  return
	end if

	iaux = d/dmaster
	daux = iaux * dmaster

	if( mode .gt. 0 ) then
	  if( daux .lt. d ) daux = daux + dmaster
	else if( mode .lt. 0 ) then
	  if( daux .gt. d ) daux = daux - dmaster
	else
	  daux = d
	end if

	rround = daux

	return
	end

!****************************************************

        function rdist(xmin,ymin,xmax,ymax)

! computes gridspacing for frame (4-7 grid lines)
!
! xmin,ymin     coordinates of lower left point
! xmax,ymax     coordinates of upper rigth point

	implicit none

	real rdist
	real xmin,ymin,xmax,ymax

	real xdist,ydist,dist,fdist
	integer istell,lines

        xdist=xmax-xmin
        ydist=ymax-ymin

        if( xdist .gt. ydist ) then
		dist = xdist
	else
		dist = ydist
	end if

        dist=nint(dist)

        istell=log10(dist)
        fdist=10**istell
        lines=dist/fdist

        if(lines.le.3) fdist=fdist*0.5
        if(lines.ge.8) fdist=fdist*2.

        rdist=fdist

        return
	end

!**************************************************************

        function divdist(x,n,mode)

! divides x into n equal pieces (with rounding)

        implicit none

	real divdist	!length part (rounded)
        real x          !length to be divided
        integer n       !number of pieces to create
        integer mode    !0: get closest to n -1:get less  +1: get more

	logical debug
        integer lines,lines_high,lines_low
        real dist_high,dist_low
        real dist,fdist

        real roundm

	debug = .false.

        dist=nint(x)

        dist_high = roundm(dist/n,+1)
        dist_low  = roundm(dist/n,-1)

        lines_high = dist / dist_low
        lines_low  = dist / dist_high

        if( mode .lt. 0 ) then
          fdist = dist_high
        else if( mode .gt. 0 ) then
          fdist = dist_low
        else
          if( lines_high-n .lt. n-lines_low ) then
            fdist = dist_low
          else
            fdist = dist_high
          end if
        end if

	if( debug ) then
        write(6,*) '-------------------'
        write(6,*) x,dist,n,mode
        write(6,*) lines_low,lines_high
        write(6,*) dist_high,dist_low
        write(6,*) fdist,int(dist/fdist)
        write(6,*) '-------------------'
	end if

        divdist = fdist

        end

!**************************************************************

	function rnext(r,mode)

	implicit none

! finds closest value to r
!
! r		to r the closest value has to be found
! mode		mode of operation
!		...  0  : r is returned (no change)
!		... >0  : the higer value is found (further from 0)
!		... <0  : the lower value is found (closer to 0)
!		... |1| : 1. 2. 2.5 5. 8.
!		... |2| : 1. 2. 5.
!		... |3| : 1. 2. 3. 4. 5. 8.
!		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.
! rnext		closest value found to r
!
! val		matrix containing the closest values to be used
!		...for each mode (in ascending order)
! nval		number of values to be used for each mode
!
! if r is too small, 0 is returned
! for negative r the lower value is the value closer to zero

	real rnext
	real r
	integer mode

	integer nmodim,nvadim
	parameter (nmodim=4,nvadim=9)
	logical bhigh
	integer nval(nmodim)
	integer m,n,nn,iexpo
	double precision val(nvadim,nmodim)
	double precision rin,rrold,rr,sign,expo
	!real val(nvadim,nmodim)
	!real rin,rrold,rr,sign,expo

	data val / 1. , 2. , 2.5 , 5. , 8. ,           0.,0.,0.,0.
     +		,  1. , 2. , 5. ,                      0.,0.,0.,0.,0.,0.
     +		,  1. , 2. , 3. , 4. , 5. , 8. ,       0.,0.,0.
     +		,  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.
     +		 /
	data nval / 5 , 3 , 6 , 9 /

	if(mode.lt.0) then	!upper or lower value
		bhigh=.false.
		m=-mode
	else
		bhigh=.true.
		m=mode
	end if

	if(m.gt.nmodim.or.m.eq.0) then	!wrong mode, return r
		rnext=r
		return
	end if

	if(abs(r).lt.1.e-5) then	!r too small, return 0
		rnext=0.
		return
	end if

	if(r.lt.0.) then	!r negative --> make positive
		rin=-r
		sign=-1.
	else
		rin=r
		sign=1.
	end if

! find next value

	n=nval(m)
	iexpo=int(log10(rin))-2
	expo=10.**iexpo
	nn=1
	rr=val(nn,m)*expo
	rrold = rr

	do while(rr.lt.rin)
		rrold=rr
		nn=nn+1
		if(nn.gt.n) then
			nn=1
			expo=expo*10.
		end if
		rr=val(nn,m)*expo
	end do

	if(bhigh.or.rr.eq.rin) then	!return next value
		rnext=rr*sign
	else
		rnext=rrold*sign
	end if

	return
	end

!*********************************************************

	function rnexta(r,mode)

! finds closest value to r (absolute, i.e., respects negative values)

	implicit none

	real rnexta
	real r
	integer mode

	real rabs,rnext
	integer m

	rabs = abs(r)
	m = mode
	if( r .lt. 0. ) m = -mode

	rabs = rnext(rabs,m)

	if( r .lt. 0. ) rabs = -rabs

	rnexta = rabs

	end

!*****************************************************************

	function rnextsub(r)

! finds best subdivision for value r
!
!		... |1| : 1. 2. 2.5 5. 8.
!		... |2| : 1. 2. 5.
!		... |3| : 1. 2. 3. 4. 5. 8.
!		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.

	implicit none

	real rnextsub
	real r

	integer i
	real eps,fact,rr,rsub
	real, save :: rdata(10) =  (/0.25,0.5,1.,1.,1.,2.,1.,2.,3.,2.5/)

	eps = 1.e-5

	fact = 1.
	rr = r

	if( rr > 1 ) then
	  do while( rr/10. > 1 )
	    fact = fact*10.
	    rr = rr / 10.
	  end do
	else
	  do while( rr < 1 )
	    fact = fact/10.
	    rr = rr * 10.
	  end do
	end if

	if( rr < 1. .or. rr > 10. ) goto 98

! now rr is in interval [1-10]

	if( abs(rr-2.5) < eps ) then	!exception for 2.5
	  rsub = 0.5
	else
	  i = nint(rr)
	  if( i < 1 .or. i > 10 ) goto 99
	  rsub = rdata(i)
	end if

	rnextsub = rsub * fact

	return
   98	continue
	write(6,*) 'r,rr: ',r,rr
	stop 'error stop rnextsub: rr out of range [1-10]'
   99	continue
	write(6,*) r,rr,fact,i
	stop 'error stop rnextsub: internal error'
	end

!************************************************************

        subroutine logvals(amin,amax,idiv,ntk,rval,aval)

! computes log scale values

        implicit none

        real amin,amax          !min/max of scale
        integer idiv            !division of log scale (see below) [1-3]
        integer ntk             !dimension of aval (in), values in aval (out)
        real rval(ntk)          !relative x values (return)
        real aval(ntk)          !log scale values (return)

        real a1,a2,aa1,aa2,val
        real a,aaux,fact,r
        real aamin,aamax
        integer ia1,ia2,i
        integer ndim,ip

        real eps
        parameter (eps=0.1)

! idiv must be in [1-3]
! idiv = 1      1 10 100
! idiv = 2      1 2 10 20 100
! idiv = 3      1 2 5 10 20 50 100

        if( idiv .lt. 1 .or. idiv .gt. 3 ) goto 99
        if( amin .le. 0. .or. amin .ge. amax ) goto 97

        ndim = ntk

        a1 = alog10(amin)
        a2 = alog10(amax)
        ia1 = a1
        ia2 = a2

        if( abs(a1-ia1) .gt. eps ) then
          aa1 = ceiling(a1)
        else
          aa1 = nint(a1)
        end if
        ia1 = nint(aa1)

        if( abs(a2-ia2) .gt. eps ) then
          aa2 = floor(a2)
        else
          aa2 = nint(a2)
        end if
        ia2 = nint(aa2)

        aamin = 10.**ia1
        aamax = 10.**ia2

        write(6,*) 'log general: ',idiv,ndim,ia2-ia1+1
        write(6,*) 'log min: ',amin,a1,aa1,ia1
        write(6,*) 'log max: ',amax,a2,aa2,ia2
        write(6,*) 'log aa: ',aamin,aamax

        ip = 0
        a = aamin
        do while( a <= aamax )
          write(6,*) ip,a
          do i=1,idiv
            fact = i
            if( i .eq. 3 ) fact = 5
            aaux = a*fact
            if( aaux >= amin .and. aaux <= amax ) then
              ip = ip + 1
              if( ip .gt. ndim ) goto 98
              aval(ip) = aaux
              r = alog10(aaux)
              r = (r-a1)/(a2-a1)
              rval(ip) = r
            end if
          end do
          a = 10. * a
        end do
        ntk = ip

        write(6,*) 'log scale: ',ntk
        write(6,*) (aval(i),i=1,ntk)
        write(6,*) (rval(i),i=1,ntk)

        return
   97   continue
        write(6,*) 'amin,amax: ',amin,amax
        stop 'error stop logval_adjust: amin,amax'
   98   continue
        write(6,*) 'ndim = ',ndim
        write(6,*) (aval(i),i=1,ndim)
        stop 'error stop logval_adjust: ndim'
   99   continue
        write(6,*) 'idiv = ',idiv
        stop 'error stop logval_adjust: idiv'
        end

!************************************************************

