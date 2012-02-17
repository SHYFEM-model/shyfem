c
c $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
c
c computes connectivity matrix
c
c revision log :
c
c 23.01.2012    ggu	new routine for release in point, connectivity
c
c*******************************************************************

	subroutine lagr_connect_continuous_points

c continuous release from points

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer np
	integer iep(nconnect_dim)
	real xp(nconnect_dim),yp(nconnect_dim)
	save np,iep,xp,yp

	integer ip,n
	real pps,ppts,dt

	integer itmonth

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

	itmonth = lagr_connect_itmonth
	pps = lagr_connect_pps				!particles per second

	if( icall .eq. 0 ) then				!initialize
	  if( pps .le. 0. ) then
	    icall = -1
	    return
	  end if
	  !call lagr_connect_get_coords_1(nconnect_dim,np,xp,yp)
	  call lagr_connect_get_coords(nconnect_dim,np,xp,yp)
	  call lagr_connect_find_elems(np,xp,yp,iep)

	  call lagr_connect_init(np,iep)
	end if

	icall = icall + 1

	call get_timestep(dt)
	ppts = dt * pps				!paricles per time step

	do ip=1,np
	  call release_on_point(ppts,iep(ip),xp(ip),yp(ip),n)
	  call lagr_connect_released(ip,n)	!count released particles
	end do

	if( mod(it,itmonth).eq.0.or.it.eq.itend ) then
	  call lagr_connect_write(np)
!	  call lagr_connect_init(np,iep)  
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_init(np,iep)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np
	integer iep(np)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie,ip,jp

	real areaele

	np_station = np

	do ie=1,nel
	  i_connect_elems(ie) = 0
	  i_connect_total(ie) = 0
	  t_connect_total(ie) = 0.
	end do

	do ip=1,np
	  ie = iep(ip)
	  i_connect_elems(ie) = ip
	  a_connect_area(ip) = areaele(ie)
	end do

	do ip=1,np
	  do jp=1,np
	    i_connect(ip,jp) = 0
	    t_connect(ip,jp) = 0.
	    if_connect(ip,jp) = 0
	    tf_connect(ip,jp) = 0.
	  end do
	  i_connect_released(ip) = 0
	end do

	end

c*******************************************************************

	subroutine lagr_connect_released(ip,n)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer ip,n

	i_connect_released(ip) = i_connect_released(ip)  + n

	end

c*******************************************************************

	subroutine lagr_connect_count(ibdy,ie,time,ic)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ibdy,ie
	real time
	integer ic

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer ie_from,ip_to,ip_from
	real tarrive

	if( ie .le. 0 ) return
	if( lagr_connect_pps .le. 0. ) return

	i_connect_total(ie) = i_connect_total(ie) + ic
	t_connect_total(ie) = t_connect_total(ie) + time

	ip_to = i_connect_elems(ie)
	if( ip_to .le. 0 ) return

	ie_from = est(ibdy)
	ip_from = i_connect_elems(ie_from)
	if( ip_from .le. 0 ) then
	  write(6,*) 'particle not coming from source: ',ibdy,ie_from
	  stop 'error stop lagr_connect_count: no source'
	end if

	i_connect(ip_from,ip_to) = i_connect(ip_from,ip_to) + ic
	t_connect(ip_from,ip_to) = t_connect(ip_from,ip_to) + time

	if( ic .le. 0 ) return

	if( ip_to .gt. 63 ) then
	  write(6,*) 'ip_to = ',ip_to
	  stop 'error stop lagr_connect_count: internal error bitmap'
	end if

	if( btest(lgr_bitmap(ibdy),ip_to) ) return	!already been there
	lgr_bitmap(ibdy) = ibset(lgr_bitmap(ibdy),ip_to)!remember for next time

	tarrive = it - tin(ibdy)
	if_connect(ip_from,ip_to) = if_connect(ip_from,ip_to) + ic
	tf_connect(ip_from,ip_to) = tf_connect(ip_from,ip_to) + tarrive

	end

c*******************************************************************

	subroutine lagr_connect_write(np)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	real hev(1)
	common /hev/hev

	integer nlv,nvar,nvers,ierr,iunit,ivar
	integer ilhv(1)
	real hlv(1)
	character*80 file,title

	integer ip,jp,ie
	integer ic,icf,itot
	real tc,tcf,pab
	real area_i,area_j

	integer ifileo

	write(127,*) np

	do ip=1,np
	  write(127,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do

	write(127,*) 'total matrix'

	do ip=1,np				!from
	  do jp=1,np				!to
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	write(127,*) 'lower matrix (column-wise)'

	do jp=1,np				!to
	  do ip=jp,np				!from
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	write(127,*) 'upper matrix (row-wise)'

	do ip=1,np				!from
	  do jp=ip,np				!to
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	file = 'connectivity.eos'
	nvers = 3
	nvar = 1
	nlv = 1
	ivar = 300
	title = 'connectivity'

	iunit = ifileo(0,file,'unformatted','new')
	call wfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
	call wseos(iunit,ilhv,hlv,hev,ierr)
	call wreos(iunit,it,ivar+1,1,ilhv,t_connect_total,ierr)
	!call wreos(iunit,it,ivar+2,nlvdim,ilhv,c,ierr)
	close(iunit)

	do ie=1,nel
	  write(128,*) ie,i_connect_total(ie),t_connect_total(ie)
	end do

	end

c*******************************************************************

	subroutine lagr_connect_write_entry(np,ip,jp)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np,ip,jp

	integer ic,icf,itot
	real tc,tcf,pab
	real area_i,area_j

c ip is from
c jp is to

	area_i = a_connect_area(ip)
	itot = i_connect_released(ip)
	area_j = a_connect_area(jp)
	ic = i_connect(ip,jp)
	tc = t_connect(ip,jp)
	icf = if_connect(ip,jp)
	tcf = tf_connect(ip,jp)
	pab = (i_connect(ip,jp)/area_j) / (itot/area_i)
	if( icf .gt. 0 ) tcf = tcf / icf
	write(127,1000) ip,jp,ic,tc,icf,tcf,pab
 1000	format(2i3,i8,e15.7,i8,2e15.7)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_find_elems(np,xp,yp,iep)

	implicit none

	integer np
	real xp(np),yp(np)
	integer iep(np)

	integer ip,ie

	do ip=1,np
	  call find_element(xp(ip),yp(ip),ie)
	  if( ie .le. 0 ) then
	    write(6,*) 'no such element: ',ip,xp(ip),yp(ip)
	    stop 'error stop find_elems_to_points: no such element'
	  end if
	  iep(ip) = ie
	end do

	end

c*******************************************************************

	subroutine lagr_connect_get_coords_1(ndim,np,xp,yp)

	implicit none

	integer ndim,np
	real xp(ndim),yp(ndim)

	np = 1
	xp(np) = 25000.
	yp(np) = 80000.

	end

c*******************************************************************

        subroutine lagr_connect_get_coords(ndim,np,xp,yp)

        implicit none

        integer ndim
        integer np,i
        real xp(ndim),yp(ndim)
        real xx,yy

        np=0

        open(1,file='coord_menor.dat',status='old',err=99)

10      continue
        read(1,*,err=98,end=12)xx,yy
        np = np +1
        if (np.gt.ndim)then
          write(6,*)'lagr_get_coords np ERROR',np,ndim
          stop
        endif
        xp(np)=xx
        yp(np)=yy
        goto 10

12      continue
        close(1)
         write(6,*)'connectivity stations',np
        do i=1,np
         write(6,*)xp(i),yp(i)
        enddo

        RETURN

99      continue
        write(6,*)'lagr_get_coords OPEN FILE ERROR'
        STOP

98      continue
        write(6,*)'lagr_get_coords READING FILE ERROR'
        STOP
        end

c******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ie,ip_station
	real r_station

	if( lagr_connect_pps .gt. 0. ) then	!connectivity active
	  ip_station = i_connect_elems(ie)
	  r_station = ip_station / float(np_station)
	else
	  ip_station = 0
	  r_station = 0.5
	end if

	end

c******************************************************************

