c
c $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
c
c computes connectivity matrix
c
c revision log :
c
c 23.01.2012    ggu	new routine for release in point, connectivity
c 22.10.2012    ggu	prepared for circular area to receive particles
c
c*******************************************************************


	subroutine lagr_connect_continuous_points(brelease)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	logical brelease	!is inside release times?

	include 'femtime.h'

	integer np
	integer iep(nconnect_dim)
	real xp(nconnect_dim),yp(nconnect_dim)
	save np,iep,xp,yp

	integer ip,n
	real pps,ppts,dt

	integer itout

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

	itout = lagr_connect_itout	!matric write frequency
	pps = lagr_connect_pps		!release of particles per second

c initialize the first time

	if( icall .eq. 0 ) then			
	  if( pps .le. 0. ) then
	    icall = -1
	    return
	  end if
	  !call lagr_connect_get_coords_1(nconnect_dim,np,xp,yp)
	  call lagr_connect_get_coords(nconnect_dim,np,xp,yp)
	  call lagr_connect_find_elems(np,xp,yp,iep)

	  call lagr_connect_init(np,xp,yp,iep)
	  call lagr_connect_reset(np)
	end if

	icall = icall + 1

	call get_timestep(dt)
	ppts = dt * pps				!particles per time step

c release particles

	if( brelease ) then			!release?
	  do ip=1,np
	    call release_on_point(ppts,iep(ip),xp(ip),yp(ip),n)
	    call lagr_connect_released(ip,n)	!count released particles
	  end do
	end if

c write with freqency itout or at the end of sim 

	if( (mod(it,itout).eq.0).or.(it.eq.itend)) then
	  call lagr_connect_write(np)
!	  call lagr_connect_reset(np)  
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_init(np,xp,yp,iep)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np
	real xp(np),yp(np)
	integer iep(np)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie,ip,jp
	real r

	real areaele

	np_station = np

	do ie=1,nel
	  i_connect_elems(ie) = 0
	end do

	if( r_connect_radius .le. 0. ) then	!receiving in only one element
	  do ip=1,np
	    ie = iep(ip)
	    i_connect_elems(ie) = ip
	    a_connect_area(ip) = areaele(ie)
	  end do
	else					!receiving in circle
	  r = r_connect_radius
	  do ip=1,np
	    a_connect_area(ip) = 0.
	    call lagr_connect_mark_elems(ip,xp(ip),yp(ip),r)
	  end do
	end if

	end

c*******************************************************************

	subroutine lagr_connect_reset(np)

c i_connect_total(ie)
c t_connect_total(ie)
c i_connect(ip,jp)
c t_connect(ip,jp)
c i_connect_released(ip)	! total number of particles released at i

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie,ip,jp

	do ie=1,nel
	  i_connect_total(ie) = 0
	  t_connect_total(ie) = 0.
	end do

	do ip=1,np
	  do jp=1,np
	    i_connect(ip,jp) = 0
	    t_connect(ip,jp) = 0.
	    if_connect(ip,jp) = 0
	    tf_connect(ip,jp) = 0.
	    agef_connect(ip,jp) = 0.
	  end do
	  i_connect_released(ip) = 0
	end do

	end

c*******************************************************************

	subroutine lagr_connect_mark_elems(ip,xp,yp,r)

c marks elements contained in circle (x,y) with radius r

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer ip
	real xp,yp
	real r

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie
	real r2,d2
	real xc,yc

	real areaele

	r2 = r*r

	do ie=1,nel
	  call baric(ie,xc,yc)
	  d2 = (xc-xp)**2 + (yc-yp)**2
	  if( d2 .le. r2 ) then
	    if( i_connect_elems(ie) .gt. 0 ) goto 99
	    i_connect_elems(ie) = ip
	    a_connect_area(ip) = a_connect_area(ip) + areaele(ie)
	  end if
	end do

	return
   99	continue
	write(6,*) 'element is in more than one station: ',ie
	write(6,*) 'stations old/new: ',i_connect_elems(ie),ip
	stop 'error stop lagr_connect_mark_elems: non unique station'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_released(ip,n)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer ip,n

	i_connect_released(ip) = i_connect_released(ip)  + n

	end

c*******************************************************************

	subroutine lagr_connect_count(ibdy,ie,ieorig,time,ic)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ibdy		!number of particle
	integer ie		!element the particle has stayed
	integer ieorig		!element from which particle came
	real time		!time particle has stayed in ie
	integer ic

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	include 'femtime.h'

	logical bin,bout
	integer ie_from,ip_to,ip_from,ip_to_orig
	integer icc
	real tarrive
	logical is_r_nan

	if( ie .le. 0 ) return
	if( lagr_connect_pps .le. 0. ) return

	if( is_r_nan(time) ) then
	  write(6,*) 'nan found: ',ibdy,ie,time
	  stop 'error stop lagr_connect_count: time is nan'
	end if

	icc = ic
	if( ie .eq. ieorig ) icc = 0

	i_connect_total(ie) = i_connect_total(ie) + icc
	t_connect_total(ie) = t_connect_total(ie) + time

	ip_to = i_connect_elems(ie)
	ip_to_orig = i_connect_elems(ieorig)

	if( ip_to .gt. 63 ) then
	  write(6,*) 'ip_to = ',ip_to
	  stop 'error stop lagr_connect_count: internal error bitmap'
	end if

	if( ip_to .le. 0 .and. ip_to_orig .le. 0 ) return

	if( ip_to .eq. ip_to_orig ) then
	  icc = 0
	else if( ip_to .gt. 0 ) then
	  lgr_bitmap_in(ibdy) = ibset(lgr_bitmap_in(ibdy),ip_to)
	else if( ip_to_orig .gt. 0 ) then
	  !can only leave if already entered
	  if( btest(lgr_bitmap_in(ibdy),ip_to_orig) ) then
	    lgr_bitmap_out(ibdy) = ibset(lgr_bitmap_out(ibdy),ip_to_orig)
	  end if
	end if
	if( ip_to .le. 0 ) return

	ie_from = est(ibdy)		!element of origin of particle
	ip_from = i_connect_elems(ie_from)
	if( ip_from .le. 0 ) then
	  write(6,*) 'particle not coming from source:'
     +				,ibdy,ie_from,ip_from
	  stop 'error stop lagr_connect_count: no source'
	end if

	! icc == 1 if particle has entered ie in this step

	bin = btest(lgr_bitmap_in(ibdy),ip_to)		!particle has entered
	bout = btest(lgr_bitmap_out(ibdy),ip_to)	!particle has left

	!only count if paricle has entered orig or if new area
	if( bin ) then
	  i_connect(ip_from,ip_to) = i_connect(ip_from,ip_to) + icc
	  t_connect(ip_from,ip_to) = t_connect(ip_from,ip_to) + time
	end if

	! here starts redundancy code

	if( bout ) return

	if( bin ) then
	  if_connect(ip_from,ip_to) = if_connect(ip_from,ip_to) + icc
	  tf_connect(ip_from,ip_to) = tf_connect(ip_from,ip_to) + time
	end if

	if( bin .and. icc .gt. 0 ) then
	  tarrive = it - tin(ibdy) 	!first arrival age
	  agef_connect(ip_from,ip_to) = agef_connect(ip_from,ip_to)
     +					+ tarrive
	end if

	end

c*******************************************************************

	subroutine lagr_connect_write(np)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	include 'femtime.h'
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

	if( np .le. 0 ) return

cmi file.127 probability definition (angel or classic)
	write(127,*) np
	do ip=1,np
	  write(127,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do
	write(127,*) 'total matrix time',it

cmi file.129 Cowen probability 
	write(129,*) np
	do ip=1,np
	  write(129,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do
	write(129,*) 'total matrix time',it

cmi file.136 Exposure 
	write(136,*) np
	do ip=1,np
	  write(136,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do
	write(136,*) 'total matrix time',it

cmi file.137 cumulate number of part. and time 
	write(137,*) np
	do ip=1,np
	  write(137,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do
	write(137,*) 'total matrix time',it

	do ip=1,np				!from
	  do jp=1,np				!to
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

!	write(127,*) 'lower matrix (column-wise)'
!
!	do jp=1,np				!to
!	  do ip=jp,np				!from
!	    call lagr_connect_write_entry(np,ip,jp)
!	  end do
!	end do

!	write(127,*) 'upper matrix (row-wise)'
!
!	do ip=1,np				!from
!	  do jp=ip,np				!to
!	    call lagr_connect_write_entry(np,ip,jp)
!	  end do
!	end do

!	file = 'connectivity.eos'
!	nvers = 3
!	nvar = 1
!	nlv = 1
!	ivar = 300
!	title = 'connectivity'

!	iunit = ifileo(0,file,'unformatted','new')
!	call wfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
!	call wseos(iunit,ilhv,hlv,hev,ierr)
!	call wreos(iunit,it,ivar+1,1,ilhv,t_connect_total,ierr)
	!call wreos(iunit,it,ivar+2,nlvdim,ilhv,c,ierr)
!	close(iunit)

cmi file.128 total released particles and time 
	write(128,*) nel
	do ie=1,nel
	  write(128,*) ie,i_connect_total(ie),t_connect_total(ie)
	end do

	end

c*******************************************************************

	subroutine lagr_connect_write_entry(np,ip,jp)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np,ip,jp,j,i

	integer ic,icf,itot
	real tc,tcf,tcfave,tcave
	real agef,agefave
	real pab,pabf,cowpab,cowpabf
	real area_i,area_j

	real daux,dtaux,naux,ntaux
	real dauxf,dtauxf,nauxf,ntauxf
	real expos,expos_f

	tcfave = 0
	tcave = 0
	agefave = 0 
	expos =0 
	expos_f =0 

c ip is from
c jp is to 

	area_i = a_connect_area(ip)	!area staz i [m^2 ? ] 
	itot = i_connect_released(ip)	!num tot part mollate da i
	area_j = a_connect_area(jp)	!area staz j
	ic = i_connect(ip,jp)		!numero palle da i a j con ridond 
	tc = t_connect(ip,jp)		!tempo totale speso da prt_i in j [s]
	icf = if_connect(ip,jp)		!numero palle da i a j senza ridond
	tcf = tf_connect(ip,jp)		!tempo totale speso .. senza ridond
	agef = agef_connect(ip,jp)

c inserire tempo di arrivo medio 


!sum su j: particelle di i=fix ricevute da qualche j    
	daux   = 0 
	dauxf  = 0
	dtaux  = 0 
	dtauxf = 0 
	do j=1,np  
	daux = daux + i_connect(ip,j)*a_connect_area(j) 	!per area?? 
	dauxf = dauxf + if_connect(ip,j)*a_connect_area(j)
cmi	daux = daux + i_connect(ip,j) 
cmi	dauxf = dauxf + if_connect(ip,j)
	dtaux = dtaux + t_connect(ip,j)
	dtauxf = dtauxf + tf_connect(ip,j)
	enddo
c	write(99,*) 'daux dauxf',ip,jp,daux,dauxf
c	write(99,*) '-----------'
	
! sum su i = particelle di tutte i ricevute da una j=fix
	naux=0
	nauxf = 0
	ntaux =0 
	ntauxf =0
	do i=1,np  
	naux = naux + i_connect(i,jp)
	nauxf =nauxf + if_connect(i,jp)
	ntaux = ntaux + t_connect(i,jp)
	ntauxf = ntauxf + tf_connect(i,jp) ! [s]
	enddo
c	write(99,*) 'naux nauxf',ip,jp,naux,nauxf
c	write(99,*) '-----------'
	
c Cowen probability 

	if(daux.eq.0)then
	write(6,*) 'warning daux 0'
	cowpab = 0. 	
	cowpabf = 0. 	
	else
	cowpab = (i_connect(ip,jp)*area_j) /daux 	!Cowen2007
	cowpabf = (if_connect(ip,jp)*area_j) /dauxf	!without redund
	endif

C Angel or classical

	if(itot.eq.0)then
	write(6,*) 'warning itot 0'
	pab =0. 
	pabf = 0.
	else 
!	pab = (i_connect(ip,jp)/area_j) / itot		!Angel definition 
!	pabf = (if_connect(ip,jp)/area_j) / itot	!without redund

	pab  = real(i_connect(ip,jp)) / real(itot)	!Definition of probability
	pabf = real(if_connect(ip,jp))/ real(itot)	!without redund
	endif

c tcfave = tempo medio trascorso in j da tutte le particelle i arrivate una volta in j e contato
c solo quella volta.  
	if( icf .gt. 0 ) agefave = (agef / icf) ! eta media di arrivo
	if( icf .gt. 0 ) tcfave = (tcf / icf  ) !tempo medio speso in j 
	if( ic .gt. 0 ) tcave = (tc / ic) !tempo medio speso in j ridond 
c Exposure 

c	expos = naux /(ntaux/86400.) 
c	expos_f = nauxf /(ntauxf/86400.) 
	expos = naux /(ntaux/3600.) 
	expos_f = nauxf /(ntauxf/3600.) 


C FINAL OUTPUT 

c	write(127,1000) ip,jp,ic,tc,icf,tcfave,pab,pabf
c 1000	format(2i3,i8,e15.7,i8,e15.7,2e15.4)

	write(127,1000) ip,jp,ic,tc,icf,tcf,agef,pab,pabf
 1000	format(2i3,i8,e15.7,i8,e15.7,3e15.4)

c	write(129,1001) ip,jp,ic,tc,icf,tcf,cowpab,cowpabf
	write(129,1001) ip,jp,ic,tcave,icf,tcfave,agefave,cowpab,cowpabf
 1001	format(2i3,i8,e15.7,i8,e15.7,3e15.4)
c 1001	format(2i3,i8,e15.7,i8,e15.7,2f9.4)

	write(136,1002) ip,jp,expos,expos_f
 1002	format(2i3,2e15.6)

	write(137,1003) ip,jp,naux,ntaux,nauxf,ntauxf
 1003	format(2i3,4e15.6)
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
	character*40 file

        np=0
	file = 'coord_menor.dat'
	file = 'connectivity_xy.dat'

        open(1,file=file,status='old',err=99)

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
	write(6,*) 'file: ',file
        write(6,*)'lagr_get_coords OPEN FILE ERROR'
	write(6,*) 'connectivity not used'
        return

98      continue
	write(6,*) 'file: ',file
        write(6,*)'lagr_get_coords READING FILE ERROR'
        STOP
        end

c******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ie
	integer ip_station	!station number
	real r_station		!becomes z (color)

	ip_station = 0

	if( lagr_connect_pps .gt. 0. ) then	!connectivity active
	  if( np_station .gt. 0 ) then		!stations given
	    ip_station = i_connect_elems(ie)
	    r_station = ip_station / float(np_station)
	  end if
	end if

	end

c******************************************************************

