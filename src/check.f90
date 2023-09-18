!
! $Id: newchk.f,v 1.40 2010-02-26 17:35:06 georg Exp $
!
! routines for various checks
!
! contents :
!
! subroutine test3d(iunit,nn)           test output for new variables
! subroutine check_all			checks arrays for sanity (shell)
! subroutine check_fem			checks arrays for sanity
! subroutine check_values		checks important variables
! subroutine tsmass(ts,z,nlvdi,tstot)   computes mass of T/S or any conc. ts
! subroutine debug_dry			writes debug information on dry areas
! subroutine debug_node(k)		writes debug information on node k
! subroutine mimafem(string)		writes some min/max values to stdout
! subroutine mass_conserve		checks mass conservation
!
! subroutine check_node(k)		debug info on node k
! subroutine check_elem(ie)		debug info on element ie
! subroutine check_nodes_in_elem(ie)	debug info on nodes in element ie
! subroutine check_elems_around_node(k) debug info on elements around node k
!
! revision log :
!
! 24.08.1998    ggu     levdbg used for debug
! 26.08.1998    ggu     subroutine tsmass transferred from newbcl0
! 26.08.1998    ggu     subroutine convol, tstvol transferred from newcon1
! 22.10.1999    ggu     igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
! 05.12.2001    ggu     always execute tstvol, more debug info
! 11.10.2002    ggu     call to setaix deleted
! 09.01.2003    ggu     some variables saved in contst
! 27.03.2003    ggu     new routine value_check
! 31.07.2003    ggu     eliminate compiler warnings
! 31.07.2003    ggu     eliminate useless variables
! 10.08.2003    ggu     new routine check_fem
! 03.09.2003    ggu     routines check and sanity_check deleted
! 03.09.2003    ggu     renamed value_check to check_values, new check_all
! 13.03.2004    ggu     write total volume to inf file
! 15.03.2005    ggu     call to check_austausch() eliminated
! 23.03.2006    ggu     changed time step to real
! 23.08.2007    ggu     test for boundary nodes using routines in testbndo.h
! 27.09.2007    ggu     deleted tstvol,tstvol1,contst, new mass_conserve
! 24.06.2008    ggu     bpresv deleted
! 06.12.2008    ggu     read vreps from STR file
! 21.01.2009    ggu     new var vrerr to stop if mass error is too high
! 23.03.2009    ggu     more debug for vrerr, new routine check_node()
! 02.04.2009    ggu     new routine check_elem()
! 06.04.2009    ggu     new check_elems_around_node, check_nodes_in_elem
! 26.02.2010    ggu     in test3d() write also meteo data
! 08.04.2010    ggu     more info in checks (depth and area)
! 17.05.2011    ggu     new routine check_set_unit() to set output unit
! 12.07.2011    ggu     loop only over actual nodes/elements, not dimensions
! 15.07.2011    ggu     new routines for CRC computation
! 25.10.2011    ggu     hlhv eliminated
! 15.05.2014    ggu     write mass error only for levdbg >= 3
! 17.09.2015    ggu     in mass_conserve aux variables are local
!
!*************************************************************
!--------------------------------------------------------------------
        module check
!--------------------------------------------------------------------

	implicit none

	integer,save    :: iucheck
	data iucheck / 6 /

!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

	subroutine test3d(iunit,nn)

! test output for new variables
!
! nn	number of first array elements to be printed

	use meteo
	use internal
	use geom_dynamic
	use depth
	use ts
	use diffusion
	use hydro_vel
	use hydro_admin
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! argument
	integer iunit
	integer nn
! parameter
	include 'param.h'
! common
	include 'femtime.h'




! local
	logical bmeteo
	integer i,l,nk,ne
	integer iu,ii

	bmeteo = .false.

	iu = iunit
	if( iu .le. 0 ) iu = 6

	write(iu,*) 'time:',it
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

!******************************************************************

	subroutine check_all

! checks arrays for sanity

        use para

	implicit none

        integer levdbg

        levdbg = nint(getpar('levdbg'))

        if( levdbg .ge. 5 ) call check_fem
        if( levdbg .ge. 2 ) call check_values

	!call mimafem('panic')

	end

!******************************************************************

	subroutine check_fem

! checks arrays for sanity

        use evgeom
        use topological
        use chezy

	implicit none

!-------------------------------------------------------
! check geom arrays
!-------------------------------------------------------

	call check_ev
	call check_geom

!-------------------------------------------------------
! check vertical structure
!-------------------------------------------------------

	call check_vertical

!-------------------------------------------------------
! check various arrays
!-------------------------------------------------------

	call check_eddy
	call check_coriolis
	call check_chezy

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

	end

!******************************************************************

	subroutine check_values

! checks important variables

	use layer_thickness
	use ts
	use hydro_baro
	use hydro_vel
	use hydro_admin
        use chk_NaN
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'

	character*16 text

	text = '*** value check'

	call check1Dr(nkn,zov,-10.d0,+10.d0,text,'zov')
	call check1Dr(nkn,znv,-10.d0,+10.d0,text,'znv')

	call check2Dr(3,3,nel,zeov,-10.d0,+10.d0,text,'zeov')
	call check2Dr(3,3,nel,zenv,-10.d0,+10.d0,text,'zenv')

	call check1Dr(nel,unv,-10000.d0,+10000.d0,text,'unv')
	call check1Dr(nel,vnv,-10000.d0,+10000.d0,text,'vnv')

	call check2Dr(nlvdi,nlv,nel,utlnv,-10000.d0,+10000.d0,text,'utlnv')
	call check2Dr(nlvdi,nlv,nel,vtlnv,-10000.d0,+10000.d0,text,'vtlnv')

	call check2Dr(nlvdi,nlv,nel,ulnv,-10.d0,+10.d0,text,'ulnv')
	call check2Dr(nlvdi,nlv,nel,vlnv,-10.d0,+10.d0,text,'vlnv')

	call check2Dr(nlvdi,nlv,nkn,tempv,-30.d0,+70.d0,text,'tempv')
	call check2Dr(nlvdi,nlv,nkn,saltv,-1.d0,+50.d0,text,'saltv')

	call check2Dr(nlvdi,nlv,nkn,hdknv,0.d0,+10000.d0,text,'hdknv')
	call check2Dr(nlvdi,nlv,nkn,hdkov,0.d0,+10000.d0,text,'hdkov')

	call check2Dr(nlvdi,nlv,nel,hdenv,0.d0,+10000.d0,text,'hdenv')
	call check2Dr(nlvdi,nlv,nel,hdeov,0.d0,+10000.d0,text,'hdeov')

	end

!**********************************************************************

        subroutine tsmass(ts,mode,nlvdi,tstot)

! computes mass of T/S or any concentration ts

        use elems_dealing

        implicit none

        integer nlvdi          !dimension of levels
        double precision ts(nlvdi,1)       !concentration on nodes
!        double precision z(3,1)             !water level
	integer mode
        double precision tstot              !total computed mass of ts
!
	include 'param.h'

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop tsmass: wrong value for mode'
	end if

	tstot = scalcont(mode,ts)

        end

!************************************************************

	subroutine debug_dry

! writes debug information on dry areas

	use geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

! common
	include 'femtime.h'

	integer ie,iweg

	iweg = 0
	do ie=1,nel
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	write(6,*) 'drydry... ',it,iweg

	end

!*************************************************************

	subroutine debug_node(k)

! writes debug information on final volume around node k (internal)

	use geom_dynamic
	use depth
	use hydro_baro
	use hydro_print
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
        use para
        use fem_util
        use time_util

	implicit none

	integer k

	include 'param.h'

! common
	include 'femtime.h'
	include 'mkonst.h'

	integer ie,ii,kk,l,i
	integer ilevel
	integer iweg
	double precision flux,dzvol,avvol
	double precision diff,rdiff
	double precision aj,uv0,uv1
	double precision b,c
	double precision dt,az,azt,azpar

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

! compute inflow into column and volume of column

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

!*************************************************************

        subroutine mimafem(string)

! writes some min/max values to stdout

	use layer_thickness
	use ts
	use hydro_print
	use hydro_vel
	use hydro_admin
        use fem_util
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	include 'param.h'

        character*(*) string
! common
	include 'femtime.h'





! local
	integer ie,l,k
	double precision u,v,w,z,s,t,c
        double precision high
        double precision zmin,zmax,umin,umax,vmin,vmax
        double precision hknmax,hkomax,henmax,heomax
        double precision utomax,utnmax,vtomax,vtnmax
        double precision hlvmax,h1vmax
        double precision bprmax
! functions

!-----------------------------------------------------
! initial check and write
!-----------------------------------------------------

        !return  !FIXME
        write(6,*) '------------------ ',string,' ',it

!-----------------------------------------------------
! check water levels and barotropic velocities
!-----------------------------------------------------

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

!-----------------------------------------------------
! check of layer thickness
!-----------------------------------------------------

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

!-----------------------------------------------------
! check of transports
!-----------------------------------------------------

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

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

        write(6,*) '------------------'

        end

!*************************************************************

	subroutine vol_mass(mode)

! computes and writes total water volume

	use shympi
        use elems_dealing
        use defnames

        implicit none

	integer mode

	include 'femtime.h'

        double precision mtot              !total computed mass of ts

	integer ninfo
	save ninfo
	data ninfo /0/

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop vol_mass: wrong value for mode'
	end if

	mtot = masscont(mode)

        if(shympi_is_master()) then
	  if( ninfo .eq. 0 ) call getinfo(ninfo)
	  write(ninfo,*) 'total_volume: ',it,mtot
	end if

        end

!*************************************************************

	subroutine mass_conserve

! checks mass conservation of single boxes (finite volumes)

	use bnd_geom
	use bnd_dynamic
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
        use para
        use elems_dealing
        use defnames
        use time_util

	implicit none

	include 'mkonst.h'

	logical berror,bdebug
	integer ie,l,ii,k,lmax,mode,ks,kss
	integer levdbg
	double precision am,az,azt,dt,azpar,ampar
	double precision areafv,b,c
	double precision ffn,ffo,ff
	double precision vmax,vrmax,vdiv,vdiff,vrdiff
	double precision abot,atop
	double precision volo,voln
	double precision ubar,vbar
	double precision vbmax,vlmax,vrbmax,vrlmax
	double precision vrwarn,vrerr
	double precision qinput
	double precision vtotmax,vvv,vvm
	double precision, allocatable :: vf(:,:)
	double precision, allocatable :: va(:,:)

	integer ninfo
	save ninfo
	data ninfo /0/

!----------------------------------------------------------------
! initialize
!----------------------------------------------------------------

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

	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
	vf = 0.
	va = 0.

!----------------------------------------------------------------
! compute horizontal divergence
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! include vertical divergence
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! check mass balance in boxes
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! barotropic
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! write diagnostic output
!----------------------------------------------------------------

!	vbmax 		!absolute error for water column
!	vrbmax 		!relative error for water column
!	vlmax 		!absolute error for each box
!	vrlmax 		!relative error for each box

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

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*************************************************************
!*************************************************************
!************ CRC computation ********************************
!*************************************************************
!*************************************************************

	subroutine check_crc

	use internal
	use depth
	use layer_thickness
	use bnd_dynamic
	use area
	use ts
	use diffusion
	use hydro_print
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
        use defnames

	implicit none

	include 'param.h'

	include 'femtime.h'




	integer icrc,iucrc
	save iucrc
	data iucrc /0/		! unit

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
	call check_crc_2d(iucrc,'wprv',nlvdi+1,nkn,ilhkv,wprv)

	if( icrc .le. 3 ) return

	end

!*************************************************************

	subroutine check_crc_2d(iu,text,nlvdi,n,levels,array)

        use chk_NaN

	implicit none

	integer iu
	character*(*) text
	integer nlvdi
	integer n
	integer levels(n)
	double precision array(nlvdi,n)

	integer crc,nlv

	nlv = 0		! use levels

	call checksum_2d(nlvdi,n,nlv,levels,array,crc)
	write(iu,*) text,'   ',n,nlvdi,nlv,crc

	end

!*************************************************************

	subroutine check_crc_1d(iu,text,n,array)

        use chk_NaN

	implicit none

	integer iu
	character*(*) text
	integer n
	double precision array(n)

	integer crc

	call checksum_1d(n,array,crc)
	write(iu,*) text,'   ',n,crc

	end

!*************************************************************

	subroutine check_set_unit(iu)

	implicit none

	integer iu

	integer iucheck
	common /iucheck/iucheck

	iucheck = iu

	end

!*************************************************************

	subroutine check_node(k)

! writes debug information on node k

	use geom_dynamic
	use depth
	use layer_thickness
	use bnd_dynamic
	use area
	use ts
	use diffusion
	use hydro_vel
	use hydro_admin
	use levels
	use basin
        use fem_util
        use elems_dealing

	implicit none

	integer k

	include 'param.h'

	integer iucheck
	common /iucheck/iucheck

	include 'femtime.h'




	integer iu
	integer l,lmax,kk

	iu = iucheck
	lmax = ilhkv(k)

	write(iu,*) '-------------------------------- check_node'
	write(iu,*) 'it,idt,k,kext: ',it,idt,k,ipext(k)
	write(iu,*) 'lmax,inodv:    ',lmax,inodv(k)
	write(iu,*) 'xgv,ygv:       ',xgv(k),ygv(k)
	write(iu,*) 'zov,znv:       ',zov(k),znv(k)
	write(iu,*) 'hdkov:         ',(hdkov(l,k),l=1,lmax)
	write(iu,*) 'hdknv:         ',(hdknv(l,k),l=1,lmax)
	write(iu,*) 'areakv:        ',(areakv(l,k),l=1,lmax)
	write(iu,*) 'volold:        ',(volnode(l,k,-1),l=1,lmax)
	write(iu,*) 'volnew:        ',(volnode(l,k,+1),l=1,lmax)
	write(iu,*) 'wlnv:          ',(wlnv(l,k),l=0,lmax)
	write(iu,*) 'mfluxv:        ',(mfluxv(l,k),l=1,lmax)
	write(iu,*) 'tempv:         ',(tempv(l,k),l=1,lmax)
	write(iu,*) 'saltv:         ',(saltv(l,k),l=1,lmax)
	write(iu,*) 'visv:          ',(visv(l,k),l=1,lmax)
	write(iu,*) 'difv:          ',(difv(l,k),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

!*************************************************************

	subroutine check_elem(ie)

! writes debug information on element ie

	use geom_dynamic
	use depth
	use layer_thickness
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
        use fem_util

	implicit none

	integer iunit
	integer ie

	include 'param.h'

	integer iucheck
	common /iucheck/iucheck

	include 'femtime.h'

	integer iu
	integer l,lmax,ii

	iu = iucheck
	lmax = ilhv(ie)

	write(iu,*) '-------------------------------- check_elem'
	write(iu,*) 'it,idt,ie,ieext:  ',it,idt,ie,ieext(ie)
	write(iu,*) 'lmax,iwegv,iwetv: ',lmax,iwegv(ie),iwetv(ie)
	write(iu,*) 'area:             ',ev(10,ie)*12.
	write(iu,*) 'nen3v  :          ',(nen3v(ii,ie),ii=1,3)
	write(iu,*) 'hev:              ',hev(ie)
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

!*************************************************************

	subroutine check_nodes_in_elem(ie)

! writes debug information on nodes in element ie

	use basin
        use fem_util

	implicit none

	integer ie

	integer iucheck
	common /iucheck/iucheck

	include 'param.h'

	integer ii,k,iu

	iu = iucheck

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes in element: ',ie,ieext(ie)
	write(iu,*) '-------------------------------------------'

	do ii=1,3
	  k = nen3v(ii,ie)
	  call check_node(k)
	end do

	end

!*************************************************************

	subroutine check_elems_around_node(k)

! writes debug information on elements around node k

	use basin
        use fem_util

	implicit none

	integer k

	integer iucheck
	common /iucheck/iucheck


	include 'param.h'

	integer ie,ii,kk,iu
	logical bdebug

	iu = iucheck

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

!*****************************************************************

	subroutine check_levels

! checks arrays hlv and hldv

	use levels
        use sigma_admin

	implicit none

! common
	include 'param.h'

! local
	logical bstop,bsigma
	integer l,nsigma,levmin
	double precision h,hd,fact,hsigma
	double precision hbot,htop

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	bstop = .false.

!--------------------------------------------------------------
! check hlv values
!--------------------------------------------------------------

	if( nsigma .gt. 0 ) then
	  h = hlv(nsigma)
	  if( h .ne. -1. .and. h .ne. hsigma ) then
	    write(6,*) h,hsigma
	    stop 'error stop check_levels: hsigma'
	  end if
	end if

	hbot = -hlv(1)
	do l=2,nsigma
	  htop = hbot
	  hbot = -hlv(l)
	  if( l .eq. nsigma ) hbot = 1.
	  if( hbot .le. htop ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	levmin = nsigma + 1
	if( levmin .eq. 1 ) levmin = 2
	do l=levmin,nlv
	  if( hlv(l) .le. hlv(l-1) ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop check_levels: level error'

!--------------------------------------------------------------
! check hldv values
!--------------------------------------------------------------

	hbot = hlv(1)
	do l=2,nlv
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hd=hbot-htop
	  if( hd .ne. hldv(l) ) then
	    write(6,*) 'Error in dlevel values for level : ',l
	    write(6,*) '   hd,hldv(l) :',hd,hldv(l)
	    bstop = .true.
	  end if
	  hbot = hlv(l)
	end do
	if( bstop ) stop 'error stop check_levels: dlevel error'

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*****************************************************************

	subroutine check_eddy

! checks vertical eddy coefficient

	use diffusion
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'



	integer k,l
	double precision v,d

!------------------------------------------------------------------
! check eddy coefficient
!------------------------------------------------------------------

	do k=1,nkn
	  do l=0,nlv
	    v = visv(l,k)
	    d = difv(l,k)
	    if( v .lt. 0. .or. v .gt. 1.e+5 ) goto 99
	    if( d .lt. 0. .or. d .gt. 1.e+5 ) goto 99
	  end do
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	write(6,*) k,l,v,d
	stop 'error stop check_eddy: error in values'
	end

!*****************************************************************

	subroutine check_vertical

! checks arrays containing vertical structure

	use levels
	use shympi

	implicit none

	call check_nlv
	call check_hlv
	call check_levels
	call check_ilevels

        if(bmpi_debug) then
	  call shympi_check_2d_node(ilhkv,'ilhkv')
	  call shympi_check_2d_node(ilmkv,'ilmkv')
	  call shympi_check_2d_elem(ilhv,'ilhv')
	  call shympi_check_2d_elem(ilmv,'ilmv')
        end if

	end

!*****************************************************************

	subroutine check_hlv

	use levels

	implicit none

	include 'param.h'

	integer l

	write(6,*) 'check_hlv: ',nlv,nlvdi
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

	end

!*****************************************************************

	subroutine check_nlv

! checks nlv and associated parameters

	use levels, only : nlvdi,nlv

	implicit none

	write(6,*) 'check_nlv : ',nlvdi,nlv

	if(nlv.gt.nlvdi) stop 'error stop check_nlv: level dimension'

	end


!*************************************************************

	subroutine check_ilevels

! checks arrays ilhv and ilhkv

	use levels
	use basin
        use sigma_admin

	implicit none

	logical bsigma,bspure
	integer nsigma
	integer ie,ii,k,lmax,lk
	double precision hmax,hsigma

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	hmax = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  if( lmax .le. 0 ) goto 99
	  do ii=1,3
	    hmax = max(hmax,hm3v(ii,ie))
	  end do
	end do

	bspure = bsigma .and. hmax .le. hsigma	!pure sigma coordinates

	do ie=1,nel
	  lmax=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    lk = ilhkv(k)
	    if( lk .le. 0 ) goto 98
	    if( lk .lt. lmax ) goto 98
	  end do
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
	end do

	do k=1,nkn
	  lmax=ilhkv(k)
	  if( lmax .le. 0 ) goto 97
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
	end do

	return
   96	continue
	write(6,*) ie,k,lmax,nsigma
	stop 'error stop check_ilevels: error in vertical structure (4)'
   97	continue
	write(6,*) k,lmax
	stop 'error stop check_ilevels: error in vertical structure (3)'
   98	continue
	write(6,*) ie,lmax,k,lk
	stop 'error stop check_ilevels: error in vertical structure (2)'
   99	continue
	write(6,*) ie,lmax
	stop 'error stop check_ilevels: error in vertical structure (1)'
	end

!*****************************************************************

	subroutine check_coriolis

! checks coriolis parameter

	use internal
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'


	integer ie
	double precision f,fmax

	fmax = 2.0 * 0.729E-4

	do ie=1,nel
	  f = fcorv(ie)
	  if( fmax - abs(f) .lt. 0. ) then
	    write(6,*) ie,f,fmax
	    stop 'error stop check_coriolis: f too big'
	  end if
	end do

	end



!*****************************************************************
!--------------------------------------------------------------------
        end module check
!--------------------------------------------------------------------
