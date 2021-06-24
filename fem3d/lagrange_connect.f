
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012-2015,2017,2019-2020  Georg Umgiesser
!    Copyright (C) 2020  Michol Ghezzo
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

! computes connectivity matrix
!
! revision log :
!
! 23.01.2012    ggu	new routine for release in point, connectivity
! 17.02.2012    ggu     changed VERS_6_1_45
! 08.10.2012    ggu     changed VERS_6_1_58
! 22.10.2012    ggu	prepared for circular area to receive particles
! 25.10.2012    ggu     changed VERS_6_1_59
! 13.06.2013    ggu     changed VERS_6_1_65
! 05.05.2014    ggu     changed VERS_6_1_74
! 19.12.2014    ggu     changed VERS_7_0_10
! 23.12.2014    ggu     changed VERS_7_0_11
! 19.01.2015    ggu     changed VERS_7_1_3
! 30.04.2015    ggu     changed VERS_7_1_9
! 21.05.2015    ggu     changed VERS_7_1_11
! 31.03.2017    ggu     changed VERS_7_5_24
! 25.05.2017    ggu     changed VERS_7_5_28
! 16.02.2019    ggu     changed VERS_7_5_60
! 31.01.2020    ggu     integrated into new framework for lagrange
! 30.04.2020    mcg     adapted to new applications
! 24.05.2021    ggu     bug fixes in lagr_connect_read_params()

c*******************************************************************

	module connectivity

	implicit none

	integer, parameter :: nconnect_dim = 50

	!real, save :: lagr_connect_pps = 1./3600.
	!real, parameter :: lagr_connect_pps = 0.
	!real, save :: r_connect_radius = 50000.	!ring
	!real, parameter :: r_connect_radius = 2000.	!canal + ritmare
	!real, parameter :: r_connect_radius = 500.	!canal + ritmare
	!real, parameter :: r_connect_radius = 0.	!canal + ritmare
	!integer, save :: lagr_connect_itout = 15*86400

	real, save :: lagr_connect_pps = 0.
	real, save :: r_connect_radius = 0.
	integer, save :: lagr_connect_itout = 0.
	character*80, save :: lagr_connect_station_file = ' '

	!matrix elements marked as ip in the circle around each station 

	integer, save :: np_station = 0.
	integer, save, allocatable :: i_connect_elems(:) 
	integer, save, allocatable :: i_connect_total(:)
	real, save, allocatable :: t_connect_total(:)

	integer, save :: i_connect_released(nconnect_dim)
	real, save :: a_connect_area(nconnect_dim)
	integer, save :: i_connect(nconnect_dim,nconnect_dim)
	real, save :: t_connect(nconnect_dim,nconnect_dim)
	integer, save :: if_connect(nconnect_dim,nconnect_dim)
	real, save :: tf_connect(nconnect_dim,nconnect_dim)
	real, save :: agef_connect(nconnect_dim,nconnect_dim)

	end module connectivity

c*******************************************************************

	subroutine lagr_connect_read_params

! icnn = number of stations
! radcnn = radius of circle [m]
! ppscnn = particles per seconds [1/sec] 
! idtcnn = frequency output for matrix [sec]
! pld_day = pelagic larval duration in [day]

	use connectivity
	use mod_lagrange

	implicit none
	
	double precision dgetpar

        call sctpar('connect')          !sets default section
        call sctfnm('connect')

	np_station = nint(dgetpar('icnn'))

	if( np_station == 0 ) return

	lagr_connect_itout = nint(dgetpar('idtcnn'))
	lagr_connect_itout = lagr_connect_itout * 86400 !fix integer/real

	r_connect_radius = dgetpar('radcnn')

	lagr_connect_pps = dgetpar('ppscnn')
	if( lagr_connect_pps > 0. ) then
	  lagr_connect_pps = 1./lagr_connect_pps
	end if

	pld_day = dgetpar('pldcnn')
	pld_day = pld_day * 86400

	call getfnm('statcnn',lagr_connect_station_file)

	if( np_station > nconnect_dim ) then
	  write(6,*) 'dimension too small for stations: '
	  write(6,*) np_station,nconnect_dim
	  write(6,*) 'please increase nconnect_dim'
	  stop 'error stop lagr_connect_read_params: nconnect_dim'
	end if

	if( lagr_connect_pps == 0.  ) then
	  write(6,*) 'connectivity: stations defined but no release'
	  stop 'error stop lagr_connect_read_params: pps == 0'
	end if

	if( lagr_connect_station_file == ' ' ) then
	  write(6,*) 'connectivity: no station file given'
	  stop 'error stop lagr_connect_read_params: no station file'
	end if

	if( pld_day == 0.  ) then
	  write(6,*) 'connectivity warning: pld 0'
	end if

	write(6,*) 'start parameters for connectivity -------------'
	write(6,*) np_station,lagr_connect_itout
	write(6,*) r_connect_radius,lagr_connect_pps,pld_day
	write(6,*) lagr_connect_station_file
	write(6,*) 'end parameters for connectivity -------------'

        call sctpar('para')
        call sctfnm('para')

	end

c*******************************************************************

	subroutine lagr_connect_continuous_points(brelease)

	use mod_lagrange
	use connectivity

	implicit none

	include 'femtime.h'

	logical brelease	!is inside release times?

	integer np
	integer, save :: iep(nconnect_dim)
	real, save :: xp(nconnect_dim),yp(nconnect_dim)

	integer ip,n,itype
	integer itout
	real pps,ppts,dt

	integer, save :: icall = 0

	if( icall .lt. 0 ) return

c initialize the first time

	if( icall .eq. 0 ) then			
	  bconnect = .false.
	  call lagr_connect_read_params
	  if( np_station == 0 ) icall = -1
	  if( icall < 0 ) return

	  call lagr_connect_get_coords(nconnect_dim,np,xp,yp)
	  if( np /= np_station ) goto 99
	  call lagr_connect_find_elems(np,xp,yp,iep)

	  call lagr_connect_init(np,xp,yp,iep)
	  call lagr_connect_reset(np)
	  bconnect = .true.
	end if

	icall = icall + 1

c set parameters

	np = np_station
	call get_timestep(dt)
	pps = lagr_connect_pps		!release of particles per second
	ppts = dt * pps				!particles per time step

c release particles TODO 3d release

	if( brelease ) then			!release?
	  do ip=1,np
	    itype = ip
	    call release_on_point(itype,ppts,iep(ip),xp(ip),yp(ip),n)
	    call lagr_connect_released(ip,n)	!count released particles

	  end do
	end if

c write with freqency itout or at the end of sim 

	itout = lagr_connect_itout	!matric write frequency
	if( (mod(it,itout).eq.0).or.(it.eq.itend)) then
	write(6,*)'writing conn matrix',it,itout,np 
	  call lagr_connect_write(np)
!	  call lagr_connect_reset(np)
	end if

	return
   99	continue
	write(6,*) 'stations declared and stations read are different: '
	write(6,*) np_station,np
	stop 'error stop lagr_connect_continuous_points: np_station'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_init(np,xp,yp,iep)

	use basin
	use mod_lagrange
	use connectivity

	implicit none

	integer np
	real xp(np),yp(np)
	integer iep(np)

	integer ie,ip,jp,npmax
	real r

	real areaele

	npmax = bit_size(bm_kind) - 1
	!npmax = 63

	if( np > npmax ) then
	  write(6,*) 'np = ',np
	  write(6,*) 'npmax = ',npmax
	  write(6,*) 'cannot handle this total number of stations'
	  write(6,*) 'please change type of integer used for bitmap'
	  stop 'error stop lagr_connect_init: internal error bitmap'
	end if

	allocate(i_connect_elems(nel))
	allocate(i_connect_total(nel))
	allocate(t_connect_total(nel))

	i_connect_elems = 0
	i_connect_total = 0
	t_connect_total = 0.

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
	    call lagr_connect_mark_elems(ip,xp(ip),yp(ip),iep(ip),r)
	  end do
	end if

	end

c*******************************************************************

	subroutine lagr_connect_reset(np)

c i_connect_total(ie)	!total number of particles went across each element
c t_connect_total(ie)	!cumulate time spent in each element passing by
c i_connect(ip,jp)	!particles exchanged from ip to jp with redundancy 
c if_connect(ip,jp)	!particles exchanged from ip to jp without
c t_connect(ip,jp)	!time spent by particles from ip in jp with redundancy
c tf_connect(ip,jp)	!time spent by particles from ip in jp without redundancy
c agef_connect(ip,jp)	!age of particles from ip in jp station without redundancy
c i_connect_released(ip)!total number of particles released at ip

	use basin
	use connectivity

	implicit none

	integer np

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

	subroutine lagr_connect_mark_elems(ip,xp,yp,iep,r)

c marks elements contained in circle (x,y) with radius r
c mgh : it doesn't work in spherical coordinates
c 
c from : https://en.wikipedia.org/wiki/Great-circle_distance
c phi1,lambda1 = lat e lon in rad point 1 
c phi2,lambda2 = lat e lon in rad point 2 
c delta_rho =central angle between them 
c delta_lam, delta_phi =  absolute difference 
c radious of sphera (Earth =6371 km ) 
c distance = d = r*delta_rho

c FDP: https://rosettacode.org/wiki/Haversine_formula#Fortran
c used for harvesine formula: for small distances
c delta_rho = 2arcsin [ sqrt (sin**2 (delta_phy/2)  + 
c				cos phi1 * cos phi2 ) * 
c				sin**2(delta_lam/2) ]				 
c

c iep = elemento in cui si trova la stazione !FDP


	use basin
	use connectivity

	implicit none

	integer ip
	real xp,yp
	real r

	integer ie,itmpc
	real r2,d2
	real xc,yc
	real fact
	
	real dum1,dum2,dum3
	real aux1,aux2,raux2
	real Rearth
	integer isphe

	integer ieext
	integer iep
	real areaele

	real :: d
	real to_radian
	real haversine
        real,parameter :: radius = 6372.8 !earth radius in km 
 
	r2 = r*r
	fact = 1.
	itmpc = 0
	Rearth = 6372.8*1000
	isphe = 1 

	 write(6,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	 write(6,*)'   WARNING:isphe=',isphe
	 write(6,*)'   CHECK in lagrange_connect_mark_elems '
	 write(6,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	
!	if(isphe.eq.1)then
	 !rsphe =to_radian(rkm,rsphe)
c	 raux2= 2 asin (sqrt (cos(xp)*cos(XX) ) !x di un punto
	 raux2 = r*r
	 r2 =  raux2
!!	endif 
	!print*,'r,rsphe,xp,yp',r,rkm,rsphe,xp,yp

	do ie=1,nel
	!if (ie.eq.1528) then !DEBUG
	  call baric(ie,xc,yc)
	  !yc=56.05137300
	  !xc=21.05666667
	!print*,'ie,xc,yc',ie,xc,yc
	  
	  if(isphe.eq.1)then
	
	     d = haversine(yp,xp,yc,xc)*1000	!in m
	     d2 = d*d
	     
  	     !print*, 'distance staz baricentro el',ie,': ',d,' km' ! distance: in km
  	     !print*, 'distance staz baricentro el',ie,': ',d,' m' ! distance: in m
	     !print*, 'raggio circ distance: ',r,' m' ! distance: in m
	     !print*, '' 
 	     !stop
 	   else

	    d2 = (fact*(xc-xp))**2 + (yc-yp)**2

	   endif



	!dum1=sin((xc-xp)/2)**2 !delta_phi
	!dum3=sin((yc-yp)/2)**2 !delta_lamb
	!dum2=cos(xc)*cos(xp) !delta_lamb
	!aux1=sqrt(dum1 + dum2 * dum3 )
	!aux2=2*asin(aux1)
	!d2=  aux2 
	!endif

	!write(6,*)'raggi',r2,d2
	  if( d2 .le. r2 ) then
	    if( i_connect_elems(ie) .gt. 0 ) goto 99
	    i_connect_elems(ie) = ip
	    a_connect_area(ip) = a_connect_area(ip) + areaele(ie)

	    itmpc=itmpc+1
	    write(6,*)'mark_elem: ',ip,ie,ieext(ie),itmpc
	    !write(66,*)ip,ie,ieext(ie),itmpc
	  end if
	!endif 
	end do
	
	  if(itmpc.eq.0)then
	    if( i_connect_elems(ie) .gt. 0 ) goto 99
	    i_connect_elems(iep) = ip
	    a_connect_area(ip) = a_connect_area(ip) + areaele(iep)
	    write(6,*)'mark_elem staz only 1 elem',ip,iep,ieext(iep),itmpc
	    
	  endif

	return
   99	continue
	write(6,*) 'element is in more than one station: ',ie
	write(6,*) 'radius probably too high, overlapping circles: ',r 
	write(6,*) 'stations old/new: ',i_connect_elems(ie),ip
	stop 'error stop lagr_connect_mark_elems: non unique station'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_released(ip,n)

	use connectivity

	implicit none

	integer ip,n

	i_connect_released(ip) = i_connect_released(ip)  + n

	end

c*******************************************************************

	subroutine lagr_connect_count(ibdy,ie_to,ie_from,time)

	use mod_lagrange
	use connectivity
	use mod_debug

	implicit none

	include 'femtime.h'

	integer ibdy		!number of particle
	integer ie_to		!element into which particle enters
	integer ie_from		!element from which particle came
	real time		!time particle has stayed in ie_from

	integer ieext
	logical bin,bout
	logical bfirstenter,bfirststay,bneverleft
	integer ie_orig,ip_to,ip_orig,ip_from
	integer l_orig
	integer icc
	double precision tinit
	real tarrive

	if( .not. bconnect ) return
	if( ie_to .le. 0 ) return
	if( ie_from .le. 0 ) stop 'error stop lagr_connect_count: (1)'

	if( is_nan(time) ) then
	  write(6,*) 'nan found: ',ibdy,ie_to,time
	  stop 'error stop lagr_connect_count: time is nan'
	end if

	!----------------------------------------------------------
	! count on elements
	!----------------------------------------------------------

	icc = 1
	if( ie_to .eq. ie_from ) then 
	  !write(6,*) 'not leaving element: ',ibdy,ie_to,ic 
	  icc = 0
	endif

	i_connect_total(ie_to) = i_connect_total(ie_to) + icc
	t_connect_total(ie_from) = t_connect_total(ie_from) + time

	!----------------------------------------------------------
	! count on stations
	!----------------------------------------------------------

	ip_to = i_connect_elems(ie_to)
	ip_from = i_connect_elems(ie_from)

	call lagr_get_initial_vals(ibdy,ie_orig,l_orig,tinit)
	ip_orig = i_connect_elems(ie_orig)

	if( ip_orig .le. 0 ) then
	  write(6,*) 'particle not coming from source: '
     +				,ibdy,ie_orig,ip_orig
	  write(6,*) 'particle not coming from source: '
     +				,ibdy,ieext(ie_orig),ip_orig
	  stop 'error stop lagr_connect_count: no source'
	end if

	if( ip_to .le. 0 .and. ip_from .le. 0 ) return	!in no station

	!----------------------------------------------------------
	! start counting
	!----------------------------------------------------------

	bfirstenter = .false.	!enters first time station ip_to
	bneverleft = .false.	!never left original station ip_orig
	bfirststay = .false.	!first stay in station ip_from

	if( ip_from > 0 ) then
	  bneverleft = .not. btest(lgr_ar(ibdy)%bitmap_in,ip_from)
	  bfirststay = .not. btest(lgr_ar(ibdy)%bitmap_out,ip_from)
	end if

	icc = 1

	if( ip_to .eq. ip_from ) then			!in same station 
	  icc = 0
	else 						!different stations
	  if( ip_to > 0 ) then
	    bfirstenter = .not. btest(lgr_ar(ibdy)%bitmap_in,ip_to)
	    lgr_ar(ibdy)%bitmap_in = ibset(lgr_ar(ibdy)%bitmap_in,ip_to)
	  end if
	  if( ip_from > 0 ) then  !can only leave if ever entered
	    if( btest(lgr_ar(ibdy)%bitmap_in,ip_from) ) then
	      lgr_ar(ibdy)%bitmap_out 
     +			= ibset(lgr_ar(ibdy)%bitmap_out,ip_from)
	    end if
	  end if
	end if

	if( .not. bneverleft ) then	!particle has left orignal area - may count
	  if( ip_to > 0 ) then
	    i_connect(ip_orig,ip_to) = i_connect(ip_orig,ip_to)+icc
	  end if
	  if( bfirstenter ) then
	    if_connect(ip_orig,ip_to) = if_connect(ip_orig,ip_to)+icc
	  end if
	  if( ip_from > 0 ) then
	    t_connect(ip_orig,ip_from) = t_connect(ip_orig,ip_from)+time
	  end if
	  if( bfirststay ) then
	    tf_connect(ip_orig,ip_from) = tf_connect(ip_orig,ip_from)+time
	  end if
!	 write(10,*)ibdy,ip_to,ip_orig,bfirststay,time
	end if

	!----------------------------------------------------------
	! start counting taking into account redundancy
	!----------------------------------------------------------

	if( bfirstenter ) then
	  tarrive = it - tinit 	!first arrival age
!	  write(8,*) ibdy,tinit,it,tarrive,ip_orig,ip_to
	  agef_connect(ip_orig,ip_to) = agef_connect(ip_orig,ip_to)
     +					+ tarrive
	end if

	end

c*******************************************************************

	subroutine lagr_connect_write(np)

	use basin
	use connectivity
	!use mod_depth

	implicit none

	include 'femtime.h'

	integer np

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
c	write(127,*) np
c	do ip=1,np
c	  write(127,*) ip,i_connect_released(ip),a_connect_area(ip)
c	end do
c	write(127,*) 'total matrix time',it

cmi file.129 Cowen probability 
c	write(129,*) np
c	do ip=1,np
c	  write(129,*) ip,i_connect_released(ip),a_connect_area(ip)
c	end do
c	write(129,*) 'total matrix time',it

cmi file.136 Exposure 
c	write(136,*) np
c	do ip=1,np
c	  write(136,*) ip,i_connect_released(ip),a_connect_area(ip)
c	end do
c	write(136,*) 'total matrix time',it

cmi file.137 cumulate number of part. and time 
	write(137,*) np
	do ip=1,np
	  write(137,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do
	write(137,*) 'total matrix time',it

	do ip=1,np				!from
	  do jp=1,np				!to
	    call lagr_connect_write_entry(it,np,ip,jp)
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

	subroutine lagr_connect_write_entry(it,np,ip,jp)

	use connectivity

	implicit none

	integer np,ip,jp,j,i
	integer it,itw

	integer ic,icf,itot
	real tc,tcf,tcfave,tcave
	real agef,agefave
	real pab,pabf,paab,paabf,cowpab,cowpabf
	real area_i,area_j

	real daux,dtaux,naux,ntaux
	real dauxf,dtauxf,nauxf,ntauxf
	real expos,expos_f

!	itw=it/(30.5*86400)
!	itw=it/3600
	itw=it/86400

	tcfave = 0
	tcave = 0
	agefave = 0 
	expos =0 
	expos_f =0 

c ip is from (source)
c jp is to (destination)

	area_i = a_connect_area(ip)	!area of station i  [m^2] source
	itot = i_connect_released(ip)	!total number of particles released by i
	area_j = a_connect_area(jp)	!area of station i  [m^2] destination

	ic = i_connect(ip,jp)		!number of particles i-> j with redundance
	tc = t_connect(ip,jp)		!total time spent by partic orig. in i in the area of j [s]

	icf = if_connect(ip,jp)		!number of particles i-> j without redundance
	tcf = tf_connect(ip,jp)		!total time spent by partic orig. in i in the area of j [s]
	agef = agef_connect(ip,jp) !age of first arrival of particles without redundance

c average arrival time TO DO 


!sum su j: particelle di i=fix ricevute da qualche j    
	daux   = 0 
	dauxf  = 0
	dtaux  = 0 
	dtauxf = 0 
	do j=1,np  
	daux = daux + i_connect(ip,j) * a_connect_area(j) 	!angel per area 
	dauxf = dauxf + if_connect(ip,j) * a_connect_area(j)

cmgh	daux = daux + i_connect(ip,j) 				!michol solo particles
cmgh	dauxf = dauxf + if_connect(ip,j)
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
	
c Cowen probability: it is not reasonable  

	if(daux.eq.0)then
	!write(6,*) 'warning daux 0'
	cowpab = 0. 	
	cowpabf = 0. 	
	else
	cowpab = (i_connect(ip,jp)*area_j) /daux 	!Cowen2007
	cowpabf = (if_connect(ip,jp)*area_j) /dauxf	!without redund
	endif

C Angel or classical

	if(itot.eq.0)then
	write(6,*) 'warning itot 0'
	pab  = 0. 
	pabf = 0.
	paab  = 0. 
	paabf = 0.
	else 
	paab = (i_connect(ip,jp)/area_j) / itot		!Angel definition 
	paabf = (if_connect(ip,jp)/area_j) / itot	!without redund

	pab  = real(i_connect(ip,jp)) / real(itot)	!Definition of probability
	pabf = real(if_connect(ip,jp))/ real(itot)	!without redund
	endif

c tcfave = tempo medio trascorso in j da tutte le particelle i arrivate una volta in j e contato
c solo quella volta.  
	if( icf .gt. 0 ) agefave = (agef / icf ) ! eta media di arrivo
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

	write(127,1000) itw,ip,jp,ic,tc,icf,tcf,agef,paab,paabf
 1000	format(i10,2i3,i8,e15.7,i8,e15.7,3e15.4)

c	write(129,1001) ip,jp,ic,tc,icf,tcf,cowpab,cowpabf
c	write(129,1001) itw,ip,jp,ic,tcave,icf,tcfave,agefave,cowpab,cowpabf
	write(129,1001) itw,ip,jp,ic,tcave,icf,tcfave,agefave,pab,pabf
c 1001	format(2i3,i8,e15.7,i8,e15.7,3e15.4)
c 1001	format(2i3,i8,e15.7,i8,e15.7,2f9.4)
 1001	format(i10,2i3,i8,e15.7,i8,e15.7,3e15.4)

c	write(136,1002) ip,jp,expos,expos_f
c 1002	format(2i3,2e15.6)

c	write(137,1003) ip,jp,naux,ntaux,nauxf,ntauxf
c 1003	format(2i3,4e15.6)
	end
c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_find_elems(np,xp,yp,iep)

	implicit none

	integer np
	real xp(np),yp(np)
	integer iep(np)
	integer ieext

	integer ip,ie

	do ip=1,np
	  call find_element(xp(ip),yp(ip),ie)
	  if( ie .le. 0 ) then
	    write(6,*) 'no such element: ',ip,xp(ip),yp(ip)
	    stop 'error stop find_elems_to_points: no such element'
	  end if
	  iep(ip) = ie
	    write(6,*) ' st ',ip,' found ie, ieext',ie,ieext(ie)
	end do

	    write(6,*) ' elements founded  '
	end

c*******************************************************************

        subroutine lagr_connect_get_coords(ndim,np,xp,yp)

	use connectivity

        implicit none

        integer ndim
        integer np
        real xp(ndim),yp(ndim)

	integer i,ios
        real xx,yy
	character*80 file

        np=0
	file = lagr_connect_station_file

        open(1,file=file,status='old',err=99)

	do
          read(1,*,iostat=ios) xx,yy
	  if( ios < 0 ) exit
	  if( ios > 0 ) goto 98
          np = np +1
          if (np.gt.ndim) goto 97
          xp(np)=xx
          yp(np)=yy
	end do

        close(1)

        write(6,*)'connectivity stations: ',np
        do i=1,np
          write(6,*) xp(i),yp(i)
        end do

        return
97      continue
	write(6,*) 'error reading file: ',file
	write(6,*) 'np,ndim: ',np,ndim
	stop 'error stop lagr_connect_get_coords: dimension error'
98      continue
	write(6,*) 'error reading file: ',file
	stop 'error stop lagr_connect_get_coords: read error'
99      continue
	write(6,*) 'error opening file: ',file
	stop 'error stop lagr_connect_get_coords: no such file'
        end

c******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	use connectivity

	implicit none

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

	subroutine lagr_connect_bitmap_init(ip)

	use mod_lagrange

	implicit none

	integer ip

        lgr_ar(ip)%bitmap_in = 0
        lgr_ar(ip)%bitmap_out = 0

	end

c******************************************************************

	subroutine lagr_connect_bitmap_copy(ifrom,ito)

	use mod_lagrange

	implicit none

	integer ifrom,ito

        lgr_ar(ito)%bitmap_in = lgr_ar(ifrom)%bitmap_in
        lgr_ar(ito)%bitmap_out = lgr_ar(ifrom)%bitmap_out

	end

c******************************************************************

 
	!contains
 
	function to_radian(degree)
          ! degrees to radians
          real,intent(in) :: degree
          real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
          real :: rad
 
          rad = degree*deg_to_rad
	  to_radian=rad
      end
 
c******************************************************************
      function haversine(deglat1,deglon1,deglat2,deglon2)
          ! great circle distance -- adapted from Matlab 
          real,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real :: a,c,dist,dlat,dlon,lat1,lat2
          real,parameter :: radius = 6372.8 
 
          dlat = to_radian(deglat2-deglat1)
          dlon = to_radian(deglon2-deglon1)
          lat1 = to_radian(deglat1)
          lat2 = to_radian(deglat2)
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          c = 2*asin(sqrt(a))
          dist = radius*c
	  haversine=dist
      end
c******************************************************************
