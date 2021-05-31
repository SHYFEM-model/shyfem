
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2019  Georg Umgiesser
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

c routines for various checks
c
c contents :
c
c subroutine test3d(iunit,nn)           test output for new variables
c subroutine check_all			checks arrays for sanity (shell)
c subroutine check_fem			checks arrays for sanity
c subroutine check_values		checks important variables
c subroutine tsmass(ts,z,nlvdi,tstot)   computes mass of T/S or any conc. ts
c subroutine debug_dry			writes debug information on dry areas
c subroutine debug_node(k)		writes debug information on node k
c subroutine mimafem(string)		writes some min/max values to stdout
c subroutine mass_conserve		checks mass conservation
c
c subroutine check_node(k)		debug info on node k
c subroutine check_elem(ie)		debug info on element ie
c subroutine check_nodes_in_elem(ie)	debug info on nodes in element ie
c subroutine check_elems_around_node(k) debug info on elements around node k
c subroutine check_nodes_around_node(k) debug info on nodes around node k
c
c revision log :
c
c 24.08.1998	ggu	levdbg used for debug
c 26.08.1998	ggu	subroutine tsmass transferred from newbcl0
c 26.08.1998	ggu	subroutine convol, tstvol transferred from newcon1
c 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
c 05.12.2001	ggu	always execute tstvol, more debug info
c 11.10.2002	ggu	call to setaix deleted
c 09.01.2003	ggu	some variables saved in contst
c 27.03.2003	ggu	new routine value_check
c 31.07.2003	ggu	eliminate compiler warnings
c 31.07.2003	ggu	eliminate useless variables
c 10.08.2003	ggu	new routine check_fem
c 03.09.2003	ggu	routines check and sanity_check deleted
c 03.09.2003	ggu	renamed value_check to check_values, new check_all
c 13.03.2004	ggu	write total volume to inf file
c 15.03.2005	ggu	call to check_austausch() eliminated
c 23.03.2006	ggu	changed time step to real
c 23.08.2007	ggu	test for boundary nodes using routines in testbndo.h
c 27.09.2007	ggu	deleted tstvol,tstvol1,contst, new mass_conserve
c 24.06.2008	ggu	bpresv deleted
c 06.12.2008	ggu	read vreps from STR file
c 21.01.2009	ggu	new var vrerr to stop if mass error is too high
c 23.03.2009	ggu	more debug for vrerr, new routine check_node()
c 02.04.2009	ggu	new routine check_elem()
c 06.04.2009	ggu	new check_elems_around_node, check_nodes_in_elem
c 26.02.2010	ggu	in test3d() write also meteo data
c 23.03.2010	ggu	changed v6.1.1
c 08.04.2010	ggu	more info in checks (depth and area)
c 20.12.2010	ggu	changed VERS_6_1_16
c 17.02.2011	ggu	changed VERS_6_1_18
c 17.05.2011	ggu	new routine check_set_unit() to set output unit
c 31.05.2011	ggu	changed VERS_6_1_23
c 12.07.2011	ggu	loop only over actual nodes/elements, not dimensions
c 15.07.2011	ggu	new routines for CRC computation
c 19.08.2011	ggu	changed VERS_6_1_29
c 25.10.2011	ggu	hlhv eliminated
c 04.11.2011	ggu	changed VERS_6_1_35
c 30.03.2012	ggu	changed VERS_6_1_51
c 05.12.2013	ggu	changed VERS_6_1_70
c 15.05.2014	ggu	write mass error only for levdbg >= 3
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 17.09.2015	ggu	in mass_conserve aux variables are local
c 29.09.2015	ggu	changed VERS_7_2_5
c 16.12.2015	ggu	changed VERS_7_3_16
c 28.04.2016	ggu	changed VERS_7_5_9
c 14.06.2016	ggu	changed VERS_7_5_14
c 17.06.2016	ggu	changed VERS_7_5_15
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 29.08.2020	ggu	added new routine check_nodes_around_node()
c 27.03.2021	ggu	some femtime.h eliminated (not all), cleanup
c 31.05.2021	ggu	write time line in node/elem debug
c
c*************************************************************

	subroutine test3d(iunit,nn)

c test output for new variables
c
c nn	number of first array elements to be printed

	use mod_meteo
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer iunit
	integer nn

	logical bmeteo
	integer i,l,nk,ne
	integer iu,ii
	double precision dtime
	character*20 aline

	bmeteo = .false.

	iu = iunit
	if( iu .le. 0 ) iu = 6

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(iu,*) 'time:',dtime,'  ',aline
	write(iu,*) 'nn  :',nn
	write(iu,*) 'nkn :',nkn
	write(iu,*) 'nel :',nel
	write(iu,*) 'nlvdi,nlv :',nlvdi,nlv	

	write(iu,*) 'hlv :'
	write(iu,*) (hlv(l),l=1,nlv)
	write(iu,*) 'hldv :'
	write(iu,*) (hldv(l),l=1,nlv)

	if(nn.eq.0) then
		nk=nkn
		ne=nel
	else
		nk=min(nn,nkn)		!$$cmplerr
		ne=min(nn,nel)
	end if

	write(iu,*) 'ilhv :'
	write(iu,*) (ilhv(i),i=1,ne)
	write(iu,*) 'fcorv :'
	write(iu,*) (fcorv(i),i=1,ne)
	write(iu,*) 'hev :'
	write(iu,*) (hev(i),i=1,ne)
	write(iu,*) 'iwegv :'
	write(iu,*) (iwegv(i),i=1,ne)

	write(iu,*) 'zov :'
	write(iu,*) (zov(i),i=1,nk)
	write(iu,*) 'znv :'
	write(iu,*) (znv(i),i=1,nk)
	write(iu,*) 'zeov :'
	write(iu,*) ((zeov(ii,i),ii=1,3),i=1,ne)
	write(iu,*) 'zenv :'
	write(iu,*) ((zenv(ii,i),ii=1,3),i=1,ne)

	if( bmeteo ) then
	write(iu,*) 'ppv :'
	write(iu,*) (ppv(i),i=1,nk)
	write(iu,*) 'wxv :'
	write(iu,*) (wxv(i),i=1,nk)
	write(iu,*) 'wyv :'
	write(iu,*) (wyv(i),i=1,nk)
	write(iu,*) 'tauxnv :'
	write(iu,*) (tauxnv(i),i=1,nk)
	write(iu,*) 'tauynv :'
	write(iu,*) (tauynv(i),i=1,nk)
	end if

	do l=1,nlv
	write(iu,*)
	write(iu,*) 'level :',l
	write(iu,*) 'ulov :'
	write(iu,*) (ulov(l,i),i=1,ne)
	write(iu,*) 'vlov :'
	write(iu,*) (vlov(l,i),i=1,ne)
	write(iu,*) 'wlov :'
	write(iu,*) (wlov(l-1,i),i=1,nk)
	write(iu,*) 'ulnv :'
	write(iu,*) (ulnv(l,i),i=1,ne)
	write(iu,*) 'vlnv :'
	write(iu,*) (vlnv(l,i),i=1,ne)
	write(iu,*) 'wlnv :'
	write(iu,*) (wlnv(l-1,i),i=1,nk)
	write(iu,*) 'utlov :'
	write(iu,*) (utlov(l,i),i=1,ne)
	write(iu,*) 'vtlov :'
	write(iu,*) (vtlov(l,i),i=1,ne)
	write(iu,*) 'utlnv :'
	write(iu,*) (utlnv(l,i),i=1,ne)
	write(iu,*) 'vtlnv :'
	write(iu,*) (vtlnv(l,i),i=1,ne)
	write(iu,*) 'visv :'
	write(iu,*) (visv(l,i),i=1,nk)
	write(iu,*) 'difv :'
	write(iu,*) (difv(l,i),i=1,nk)
	write(iu,*) 'tempv :'
	write(iu,*) (tempv(l,i),i=1,nk)
	write(iu,*) 'saltv :'
	write(iu,*) (saltv(l,i),i=1,nk)
	write(iu,*) 'difhv :'
	write(iu,*) (difhv(l,i),i=1,ne)
	end do

	end

c******************************************************************

	subroutine check_all

c checks arrays for sanity

	implicit none

        integer levdbg
        real getpar

        levdbg = nint(getpar('levdbg'))

        if( levdbg .ge. 5 ) call check_fem
        if( levdbg .ge. 2 ) call check_values

	!call mimafem('panic')

	end

c******************************************************************

	subroutine check_fem

c checks arrays for sanity

	implicit none

c-------------------------------------------------------
c check geom arrays
c-------------------------------------------------------

	call check_ev
	call check_geom

c-------------------------------------------------------
c check vertical structure
c-------------------------------------------------------

	call check_vertical

c-------------------------------------------------------
c check various arrays
c-------------------------------------------------------

	call check_eddy
	call check_coriolis
	call check_chezy

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c******************************************************************

	subroutine check_values

c checks important variables

	use mod_layer_thickness
	use mod_ts
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*16 text

	text = '*** check_values'

	call check1Dr(nkn,zov,-10.,+10.,text,'zov')
	call check1Dr(nkn,znv,-10.,+10.,text,'znv')

	call check2Dr(3,3,nel,zeov,-10.,+10.,text,'zeov')
	call check2Dr(3,3,nel,zenv,-10.,+10.,text,'zenv')

	call check1Dr(nel,unv,-10000.,+10000.,text,'unv')
	call check1Dr(nel,vnv,-10000.,+10000.,text,'vnv')

	call check2Dr(nlvdi,nlv,nel,utlnv,-10000.,+10000.,text,'utlnv')
	call check2Dr(nlvdi,nlv,nel,vtlnv,-10000.,+10000.,text,'vtlnv')

	call check2Dr(nlvdi,nlv,nel,ulnv,-10.,+10.,text,'ulnv')
	call check2Dr(nlvdi,nlv,nel,vlnv,-10.,+10.,text,'vlnv')

	call check2Dr(nlvdi,nlv,nkn,tempv,-30.,+70.,text,'tempv')
	call check2Dr(nlvdi,nlv,nkn,saltv, -1.,+70.,text,'saltv')

	call check2Dr(nlvdi,nlv,nkn,hdknv,0.,+10000.,text,'hdknv')
	call check2Dr(nlvdi,nlv,nkn,hdkov,0.,+10000.,text,'hdkov')

	call check2Dr(nlvdi,nlv,nel,hdenv,0.,+10000.,text,'hdenv')
	call check2Dr(nlvdi,nlv,nel,hdeov,0.,+10000.,text,'hdeov')

	end

c**********************************************************************

        subroutine tsmass(ts,mode,nlvdi,tstot)

c computes mass of T/S or any concentration ts

        implicit none

        integer nlvdi          !dimension of levels
        real ts(nlvdi,1)       !concentration on nodes
c        real z(3,1)             !water level
	integer mode
        real tstot              !total computed mass of ts

	double precision scalcont

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop tsmass: wrong value for mode'
	end if

	tstot = scalcont(mode,ts)

        end

c************************************************************

	subroutine debug_dry

c writes debug information on dry areas

	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,iweg
	double precision dtime
	character*20 aline

	iweg = 0
	do ie=1,nel
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(6,*) 'drydry... ',dtime,'  ',aline

	end

c*************************************************************

	subroutine debug_node(k)

c writes debug information on final volume around node k (internal)

	use mod_geom_dynamic
	use mod_depth
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	integer k

c common
	include 'femtime.h'
	include 'mkonst.h'

	integer ie,ii,kk,l,i
	integer ilevel
	integer iweg
	real flux,dzvol,avvol
	real diff,rdiff
	real aj,uv0,uv1
	real b,c
	real dt,az,azt,azpar

	real getpar
	integer ipext,ieext

	integer, save :: netot
	integer, save, allocatable :: kinf(:,:)

	integer kmem
	save kmem
	data kmem / 0 /

	if( k == 0 ) then
	  allocate(kinf(2,ngr))
	  kinf = 0
	end if

	if( k .ne. kmem ) then
	  netot = 0
          do ie=1,nel
            do ii=1,3
	      kk = nen3v(ii,ie)
	      if( kk .eq. k ) then
	        netot = netot + 1
	        if( netot .gt. ngr ) then
		  stop 'error stop debug_node: ngr'
	        end if
	        kinf(1,netot) = ie
	        kinf(2,netot) = ii
	      end if
	    end do
	  end do
	  kmem = k
	  write(6,*) 'new node for debug...'
	  write(6,*) k,ipext(k),netot
	  do i=1,netot
	    ie = kinf(1,i)
	    ii = kinf(2,i)
	    write(6,*) ie,ieext(ie),ii
	  end do
	end if

c compute inflow into column and volume of column

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar
	azt = 1. - az

	flux = 0.	! flux into water column
	dzvol = 0.	! volume change due to water level change
	avvol = 0.	! average volume of water column
	iweg = 0

	do i=1,netot
	  ie = kinf(1,i)
	  ii = kinf(2,i)
          aj=ev(10,ie)
          ilevel=ilhv(ie)
          kk=nen3v(ii,ie)
	  if( kk .ne. k ) stop 'error stop debug_node: internal error'
          b=ev(ii+3,ie)
          c=ev(ii+6,ie)
          uv0=0.
          uv1=0.
          do l=ilevel,1,-1
            uv1=uv1+utlnv(l,ie)*b+vtlnv(l,ie)*c
            uv0=uv0+utlov(l,ie)*b+vtlov(l,ie)*c
          end do
          uv1=unv(ie)*b+vnv(ie)*c
          uv0=uov(ie)*b+vov(ie)*c
          flux  = flux  + dt*12.*aj*(uv0*azt+uv1*az)
          dzvol = dzvol + 4.*aj*(zenv(ii,ie)-zeov(ii,ie))
          avvol = avvol + 4.*aj*hev(ie)
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
        end do

	diff = abs(flux-dzvol)
	rdiff = 0.
	if( avvol .gt. 0. ) rdiff = diff/avvol

	write(6,*) 'debug... ',it,diff,rdiff,iweg

	end

c*************************************************************

        subroutine mimafem(string)

c writes some min/max values to stdout

	use mod_layer_thickness
	use mod_ts
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        character*(*) string

	integer ie,l,k
	real u,v,w,z,s,t,c
        real high
        real zmin,zmax,umin,umax,vmin,vmax
        real hknmax,hkomax,henmax,heomax
        real utomax,utnmax,vtomax,vtnmax
        real hlvmax,h1vmax
        real bprmax

	integer ipext,ieext

	double precision dtime
	character*20 aline

c-----------------------------------------------------
c initial check and write
c-----------------------------------------------------

        !return  !FIXME

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

        write(6,*) '------------------ ',trim(string)
        write(6,*) '   time: ',dtime,'  ',aline

c-----------------------------------------------------
c check water levels and barotropic velocities
c-----------------------------------------------------

        high = 1.e+30

        zmin =  high
        zmax = -high
        umin =  high
        umax = -high
        vmin =  high
        vmax = -high

	do k=1,nkn
	  z = znv(k)
	  u = up0v(k)
	  v = vp0v(k)
          zmin = min(zmin,z)
          zmax = max(zmax,z)
          umin = min(umin,u)
          umax = max(umax,u)
          vmin = min(vmin,v)
          vmax = max(vmax,v)
	end do

        write(6,*) zmin,zmax,umin,umax,vmin,vmax

c-----------------------------------------------------
c check of layer thickness
c-----------------------------------------------------

        hknmax = -high
        hkomax = -high
        do k=1,nkn
          do l=1,nlv
            hknmax = max(hknmax,hdknv(l,k))
            hkomax = max(hkomax,hdkov(l,k))
          end do
        end do

        henmax = -high
        heomax = -high
        do ie=1,nel
          do l=1,nlv
            henmax = max(henmax,hdenv(l,ie))
            heomax = max(heomax,hdeov(l,ie))
          end do
        end do

        write(6,*) hknmax,hkomax,henmax,heomax

c-----------------------------------------------------
c check of transports
c-----------------------------------------------------

        utomax = -high
        utnmax = -high
        vtomax = -high
        vtnmax = -high
        do ie=1,nel
          do l=1,nlv
            utomax = max(utomax,abs(utlov(l,ie)))
            utnmax = max(utnmax,abs(utlnv(l,ie)))
            vtomax = max(vtomax,abs(vtlov(l,ie)))
            vtnmax = max(vtnmax,abs(vtlnv(l,ie)))
          end do
        end do

        write(6,*) utomax,utnmax,vtomax,vtnmax

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

        write(6,*) '------------------'

        end

c*************************************************************

	subroutine vol_mass(mode)

c computes and writes total water volume

	use shympi

        implicit none

	integer mode

        real mtot              !total computed mass of ts
	double precision masscont
	character*20 aline

	integer, save :: ninfo = 0

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop vol_mass: wrong value for mode'
	end if

	mtot = masscont(mode)

        if(shympi_is_master()) then
	  if( ninfo .eq. 0 ) call getinfo(ninfo)
	  call get_act_timeline(aline)
	  write(ninfo,*) 'total_volume: ',aline,mtot
	end if

        end

c*************************************************************

	subroutine mass_conserve

c checks mass conservation of single boxes (finite volumes)

	use mod_bound_geom
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	include 'mkonst.h'

	logical berror,bdebug
	integer ie,l,ii,k,lmax,mode,ks,kss
	integer levdbg
	real am,az,azt,dt,azpar,ampar
	real areafv,b,c
	real ffn,ffo,ff
	real vmax,vrmax,vdiv,vdiff,vrdiff
	real abot,atop
	real volo,voln
	real ubar,vbar
	real vbmax,vlmax,vrbmax,vrlmax
	real vrwarn,vrerr
	real qinput
	character*20 aline
	double precision vtotmax,vvv,vvm
	real, allocatable :: vf(:,:)
	real, allocatable :: va(:,:)

	real volnode,areanode,getpar

	integer ninfo
	save ninfo
	data ninfo /0/

c----------------------------------------------------------------
c initialize
c----------------------------------------------------------------

	if( ninfo .eq. 0 ) call getinfo(ninfo)

	vrwarn = getpar('vreps')
	vrerr = getpar('vrerr')
	levdbg = nint(getpar('levdbg'))

	if( levdbg .le. 1 ) return

	mode = +1
        call getazam(azpar,ampar)
	az = azpar
	am = ampar
        azt = 1. - az
	call get_timestep(dt)

	call get_act_timeline(aline)

	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
	vf = 0.
	va = 0.

c----------------------------------------------------------------
c compute horizontal divergence
c----------------------------------------------------------------

        do ie=1,nel
          areafv = 4. * ev(10,ie)               !area of triangle / 3
          lmax = ilhv(ie)
          do l=1,lmax
            do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
                ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
                ffo = utlov(l,ie)*b + vtlov(l,ie)*c
                ff = ffn * az + ffo * azt
                vf(l,k) = vf(l,k) + 3. * areafv * ff
                va(l,k) = va(l,k) + areafv
            end do
          end do
        end do

c----------------------------------------------------------------
c include vertical divergence
c----------------------------------------------------------------

	ks = 1000
	ks = 5071
	ks = 0
	if( ks .gt. 0 ) then
	  k = ks
	  lmax = ilhkv(k)
	  write(77,*) '-------------'
	  write(77,*) k,lmax
	  write(77,*) (vf(l,k),l=1,lmax)
	  write(77,*) (wlnv(l,k),l=1,lmax)
	  vtotmax = 0.
	  do l=1,lmax
	    vtotmax = vtotmax + vf(l,k)
	  end do
	  write(77,*) 'from box: ',vtotmax
	end if

	vtotmax = 0.
	do k=1,nkn
          lmax = ilhkv(k)
	  abot = 0.
	  vvv = 0.
	  vvm = 0.
	  do l=lmax,1,-1
	    atop = va(l,k)
	    vdiv = wlnv(l,k)*abot - wlnv(l-1,k)*atop
	    vf(l,k) = vf(l,k) + vdiv + mfluxv(l,k)
	    abot = atop
	    vvv = vvv + vdiv
	    vvm = vvm + mfluxv(l,k)
	    if( k .eq. ks ) write(77,*) 'vdiv: ',l,vf(l,k),vdiv,vvv
	  end do
	  vtotmax = max(vtotmax,abs(vvv))
	  if( k .eq. ks ) write(77,*) 'vvv: ',vvv,vvm
	end do

c----------------------------------------------------------------
c check mass balance in boxes
c----------------------------------------------------------------

	kss = 5226
	kss = 0
	berror = .false.
	vrmax = 0.
	vmax = 0.
	do k=1,nkn
	  if( is_zeta_boundary(k) ) cycle
	  if( is_external_boundary(k) ) cycle
	  bdebug = k .eq. kss
	  if( bdebug ) write(78,*) '============================='
	  berror = .false.
          lmax = ilhkv(k)
	  do l=1,lmax
	    voln = volnode(l,k,+1)
	    volo = volnode(l,k,-1)
	    vdiv = vf(l,k)
	    vdiff = voln - volo - vdiv * dt
	    vdiff = abs(vdiff)
	    vrdiff = vdiff / volo
	    vmax = max(vmax,vdiff)
	    vrmax = max(vrmax,vrdiff)
	    if( bdebug ) then
	        write(78,*) l,k
	        write(78,*) volo,voln,vdiff,vrdiff
	        write(78,*) vdiv,vdiv*dt
	    end if
	    if( vrdiff .gt. vrerr ) then
		berror = .true.
	        write(6,*) 'mass_conserve: ',l,k
	        write(6,*) volo,voln,vdiff,vrdiff
	        write(6,*) vdiv,vdiv*dt
	    end if
	    if( k == kss .and. l == 1 ) then
	      write(182,*) aline,vdiff,vrdiff
	    end if
	  end do
	  if( berror ) call check_node(k)
	  if( bdebug ) then
		call check_set_unit(78)
		call check_node(k)
		write(78,*) '============================'
	  end if
	end do

	vlmax = vmax		!absolute error for each box
	vrlmax = vrmax		!relative error for each box

c----------------------------------------------------------------
c barotropic
c----------------------------------------------------------------

	do k=1,nkn
	    vf(1,k) = 0.
	    va(1,k) = 0.
	end do

        do ie=1,nel
          areafv = 4. * ev(10,ie)               !area of triangle / 3

	  ubar = 0.
	  vbar = 0.
          lmax = ilhv(ie)
          do l=1,lmax
	    ubar = ubar + az * utlnv(l,ie) + azt * utlov(l,ie)
	    vbar = vbar + az * vtlnv(l,ie) + azt * vtlov(l,ie)
	  end do

          do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
		ff = ubar * b + vbar * c
                vf(1,k) = vf(1,k) + 3. * areafv * ff
                va(1,k) = va(1,k) + areafv
          end do
        end do

	if( ks .gt. 0 ) write(77,*) 'from baro: ',vf(1,ks)

	vrmax = 0.
	vmax = 0.
	do k=1,nkn
	 !if( is_inner(k) ) then
	 if( .not. is_external_boundary(k) ) then
           lmax = ilhkv(k)
	   voln = 0.
	   volo = 0.
	   qinput = 0.
	   do l=1,lmax
	     voln = voln + volnode(l,k,+1)
	     volo = volo + volnode(l,k,-1)
	     qinput = qinput + mfluxv(l,k)
	   end do
	   vdiv = vf(1,k) + rqv(k)
	   vdiv = vf(1,k) + qinput	!should be the same
	   vdiff = voln - volo - vdiv * dt
	   if( k .eq. ks ) write(77,*) 'vdiff: ',vdiff
	   !if( vdiff .gt. 0.1 ) write(6,*) 'baro error: ',k,vdiff
	   vdiff = abs(vdiff)
	   vrdiff = vdiff / volo
	   vmax = max(vmax,vdiff)
	   vrmax = max(vrmax,vrdiff)
	 end if
	end do

	vbmax = vmax		!absolute error for water column
	vrbmax = vrmax		!relative error for water column

c----------------------------------------------------------------
c write diagnostic output
c----------------------------------------------------------------

c	vbmax 		!absolute error for water column
c	vrbmax 		!relative error for water column
c	vlmax 		!absolute error for each box
c	vrlmax 		!relative error for each box

	if( vrlmax .gt. vrwarn ) then
	  if( levdbg .ge. 3 ) then
	    write(6,*) 'mass error: ',vbmax,vlmax,vrbmax,vrlmax
	  end if
	  if( vrlmax .gt. vrerr ) then
	    write(6,*) 'mass error of matrix solution is very high'
	    write(6,*) 'the relative mass error is = ',vrlmax
	    write(6,*) 'the limit of the mass error is vrerr = ',vrerr
	    write(6,*) 'Probably there is some problem with the solution'
	    write(6,*) 'of the system matrix. However, if you think'
	    write(6,*) 'you can live with this mass error, then please'
	    write(6,*) 'increase the value of vrerr in the STR file.'
	    stop 'error stop mass_conserve: mass error too high'
	  end if
	end if

	write(ninfo,*) 'mass_balance: ',vbmax,vlmax,vrbmax,vrlmax

	deallocate(vf,va)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c************ CRC computation ********************************
c*************************************************************
c*************************************************************

	subroutine check_crc

	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer icrc,iucrc
	save iucrc
	data iucrc /0/		! unit

	integer ifemop

	icrc = 0		! level of output [0-10]

	if( icrc .le. 0 ) return

	if( iucrc .le. 0 ) then
	  iucrc = ifemop('.crc','form','new')
	  if( iucrc .le. 0 ) stop 'error stop check_crc: open file'
	end if

	write(iucrc,*) '====================================='
	write(iucrc,*) ' crc check at time = ',it
	write(iucrc,*) '====================================='

	call check_crc_1d(iucrc,'znv',nkn,znv)
	call check_crc_1d(iucrc,'zenv',3*nel,zenv)
	call check_crc_2d(iucrc,'utlnv',nlvdi,nel,ilhv,utlnv)
	call check_crc_2d(iucrc,'vtlnv',nlvdi,nel,ilhv,vtlnv)
	call check_crc_2d(iucrc,'saltv',nlvdi,nkn,ilhkv,saltv)
	call check_crc_2d(iucrc,'tempv',nlvdi,nkn,ilhkv,tempv)
	call check_crc_2d(iucrc,'rhov',nlvdi,nkn,ilhkv,rhov)

	if( icrc .le. 1 ) return

	!call check_crc_1d(iucrc,'ev',evdim*nel,ev)	!FIXME - double
	call check_crc_1d(iucrc,'hev',nel,hev)
	call check_crc_1d(iucrc,'fcorv',nel,fcorv)
	call check_crc_2d(iucrc,'visv',nlvdi,nkn,ilhkv,visv)
	call check_crc_2d(iucrc,'difv',nlvdi,nkn,ilhkv,difv)
	call check_crc_2d(iucrc,'hdknv',nlvdi,nkn,ilhkv,hdknv)
	call check_crc_2d(iucrc,'hdenv',nlvdi,nel,ilhv,hdenv)

	if( icrc .le. 2 ) return

	call check_crc_2d(iucrc,'ulnv',nlvdi,nel,ilhv,ulnv)
	call check_crc_2d(iucrc,'vlnv',nlvdi,nel,ilhv,vlnv)
	call check_crc_2d(iucrc,'mfluxv',nlvdi,nkn,ilhkv,mfluxv)
	call check_crc_2d(iucrc,'areakv',nlvdi,nkn,ilhkv,areakv)
	call check_crc_2d(iucrc,'wlnv',nlvdi+1,nkn,ilhkv,wlnv)
	call check_crc_2d(iucrc,'wprv',nlvdi,nkn,ilhkv,wprv)

	if( icrc .le. 3 ) return

	end

c*************************************************************

	subroutine check_crc_2d(iu,text,nlvdi,n,levels,array)

	implicit none

	integer iu
	character*(*) text
	integer nlvdi
	integer n
	integer levels(n)
	real array(nlvdi,n)

	integer crc,nlv

	nlv = 0		! use levels

	call checksum_2d(nlvdi,n,nlv,levels,array,crc)
	write(iu,*) text,'   ',n,nlvdi,nlv,crc

	end

c*************************************************************

	subroutine check_crc_1d(iu,text,n,array)

	implicit none

	integer iu
	character*(*) text
	integer n
	real array(n)

	integer crc

	call checksum_1d(n,array,crc)
	write(iu,*) text,'   ',n,crc

	end

c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************

	module check_unit

	integer, save :: iucheck = 6

	end module check_unit

c*************************************************************

	subroutine check_set_unit(iu)

	use check_unit

	implicit none

	integer iu

	iucheck = iu

	end

c*************************************************************

	subroutine check_get_unit(iu)

	use check_unit

	implicit none

	integer iu

	iu = iucheck

	end

c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************
c*************************************************************

	subroutine check_node(k)

c writes debug information on node k

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin
	use mod_hydro_print
        use mod_nohyd

	implicit none

	integer k

	include 'femtime.h'

	integer iu
	integer l,lmax,kk
	character*20 aline

	integer ipext
	real volnode

	call check_get_unit(iu)
	lmax = ilhkv(k)
	call get_act_timeline(aline)

	write(iu,*) '-------------------------------- check_node'
	write(iu,*) 'time:          ',aline
	write(iu,*) 'it,idt,k,kext: ',it,idt,k,ipext(k)
	write(iu,*) 'lmax,inodv:    ',lmax,inodv(k)
	write(iu,*) 'xgv,ygv:       ',xgv(k),ygv(k)
	write(iu,*) 'zov,znv:       ',zov(k),znv(k)
	write(iu,*) 'hkv,hkv+znv:   ',hkv(k),hkv(k)+znv(k)
	write(iu,*) 'hdkov:         ',(hdkov(l,k),l=1,lmax)
	write(iu,*) 'hdknv:         ',(hdknv(l,k),l=1,lmax)
	write(iu,*) 'areakv:        ',(areakv(l,k),l=1,lmax)
	write(iu,*) 'volold:        ',(volnode(l,k,-1),l=1,lmax)
	write(iu,*) 'volnew:        ',(volnode(l,k,+1),l=1,lmax)
	write(iu,*) 'wlnv:          ',(wlnv(l,k),l=0,lmax)
	write(iu,*) 'mfluxv:        ',(mfluxv(l,k),l=1,lmax)
	write(iu,*) 'tempv:         ',(tempv(l,k),l=1,lmax)
	write(iu,*) 'saltv:         ',(saltv(l,k),l=1,lmax)
	write(iu,*) 'visv:          ',(visv(l,k),l=0,lmax)
	write(iu,*) 'difv:          ',(difv(l,k),l=0,lmax)
	write(iu,*) 'qpnv:          ',(qpnv(l,k),l=1,lmax)
	write(iu,*) 'uprv:          ',(uprv(l,k),l=1,lmax)
	write(iu,*) 'vprv:          ',(vprv(l,k),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

c*************************************************************

	subroutine check_elem(ie)

c writes debug information on element ie

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	integer ie

	include 'femtime.h'

	integer iu
	integer l,lmax,ii
	real zmed
	character*20 aline

	integer ieext

	call check_get_unit(iu)
	lmax = ilhv(ie)
	zmed = sum(zenv(:,ie))/3.
	call get_act_timeline(aline)

	write(iu,*) '-------------------------------- check_elem'
	write(iu,*) 'time:             ',aline
	write(iu,*) 'it,idt,ie,ieext:  ',it,idt,ie,ieext(ie)
	write(iu,*) 'lmax,iwegv,iwetv: ',lmax,iwegv(ie),iwetv(ie)
	write(iu,*) 'area:             ',ev(10,ie)*12.
	write(iu,*) 'nen3v  :          ',(nen3v(ii,ie),ii=1,3)
	write(iu,*) 'hev,hev+zenv:     ',hev(ie),hev(ie)+zmed
	write(iu,*) 'hm3v:             ',(hm3v(ii,ie),ii=1,3)
	write(iu,*) 'zeov:             ',(zeov(ii,ie),ii=1,3)
	write(iu,*) 'zenv:             ',(zenv(ii,ie),ii=1,3)
	write(iu,*) 'zov:              ',(zov(nen3v(ii,ie)),ii=1,3)
	write(iu,*) 'znv:              ',(znv(nen3v(ii,ie)),ii=1,3)
	write(iu,*) 'hdeov:            ',(hdeov(l,ie),l=1,lmax)
	write(iu,*) 'hdenv:            ',(hdenv(l,ie),l=1,lmax)
	write(iu,*) 'utlov:            ',(utlov(l,ie),l=1,lmax)
	write(iu,*) 'vtlov:            ',(vtlov(l,ie),l=1,lmax)
	write(iu,*) 'utlnv:            ',(utlnv(l,ie),l=1,lmax)
	write(iu,*) 'vtlnv:            ',(vtlnv(l,ie),l=1,lmax)
	write(iu,*) 'ulnv:             ',(ulnv(l,ie),l=1,lmax)
	write(iu,*) 'vlnv:             ',(vlnv(l,ie),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

c*************************************************************

	subroutine check_nodes_in_elem(ie)

c writes debug information on nodes in element ie

	use basin

	implicit none

	integer ie

	integer ii,k,iu
	integer ieext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes in element: ',ie,ieext(ie)
	write(iu,*) '-------------------------------------------'

	do ii=1,3
	  k = nen3v(ii,ie)
	  call check_node(k)
	end do

	end

c*************************************************************

	subroutine check_elems_around_node(k)

c writes debug information on elements around node k

	use basin

	implicit none

	integer k

	integer ie,ii,kk,iu
	logical bdebug

	integer ipext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking elements around node: ',k,ipext(k)
	write(iu,*) '-------------------------------------------'

	do ie=1,nel
	  bdebug = .false.
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) bdebug = .true.
	  end do
	  if( bdebug ) call check_elem(ie)
	end do

	end

c*************************************************************

	subroutine check_nodes_around_node(k)

c writes debug information on nodes around node k

	use basin

	implicit none

	integer k

	integer n,i,kk,iu
	integer nodes(ngr)

	integer ipext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes around node: ',k,ipext(k)
	write(iu,*) '-------------------------------------------'

	call get_nodes_around(k,ngr,n,nodes)

	do i=1,n
	  kk = nodes(i)
	  call check_node(kk)
	end do

	end

c*************************************************************

