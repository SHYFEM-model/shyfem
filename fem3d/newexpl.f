c
c $Id: newexpl.f,v 1.10 2010-03-08 17:46:45 georg Exp $
c
c explicit term routines
c
c contents :
c
c subroutine set_explicit
c subroutine viscous_stability(ahpar,ahstab)	computes stability for viscosity
c subroutine set_diff_horizontal
c subroutine set_advective
c subroutine set_semi_lagrange
c subroutine set_barocl
c subroutine set_barocl_new
c subroutine set_barocl_new1
c
c revision log :
c
c 01.05.2007	ggu	new file -> all explicit terms here
c 28.09.2007	ggu	semi-lagrangian part introduced
c 16.04.2008	ggu	bugfix in set_barocl (real do indices!!)
c 14.07.2008	ggu&ccf	ahpar is real in set_diff_horizontal_new
c 03.11.2008    ggu&dbf nudging implemented (call to bclnudge)
c 09.11.2008	ggu	set_barocl_new (cleaned version of set_barocl)
c 19.11.2008	ggu	new set_diff_horizontal_new1(), viscous_stability()
c 19.02.2010	ggu	in viscous_stability() for dt=1
c 26.02.2010	ggu	new call to momentum_viscous_stability()
c 26.02.2010	ggu	set_advective() cleaned up
c 26.02.2010	ggu	new momentum_advective_stability()
c 08.03.2010	ggu	run only down to avail layers (bug fix)
c 16.12.2010	ggu	barocl preconditioned for sigma layers, but not finshed
c 20.05.2011	ggu	compute statistics of stability, no stab in dry elemes
c 25.08.2011	dbf&ggu	baroclinic gradient for sigma level integrated
c 25.10.2011	dbf&ggu	bug fix in set_barocl_new_interface (psigma)
c 04.11.2011    ggu     adapted for hybrid coordinates
c 10.05.2013	dbf&ggu	new routines for vertical advection (bvertadv)
c 10.05.2013	dbf&ggu	new routines for non-hydro
c 25.05.2013	ggu	new version for vertical advection (bvertadv)
c 13.09.2013	dbf&ggu	new sigma layer adjustment integrated
c 10.04.2014	ggu	use rlin and rdistv to determin advective contribution
c 17.04.2015	ggu	only one routine set_diff_horizontal()
c 18.09.2015	ggu	use momentx/yv to store advective terms, not aux arrays
c 25.09.2015	ggu	new call to set_nudging()
c 14.06.2016	dbf	diff. vertical momentum advection schemes implemented
c
c notes :
c
c sign of explicit term is computed for left hand side
c
c******************************************************************

	subroutine set_explicit

	use mod_internal
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none
        
	integer ie,l
        
        logical bbarcl
        integer ilin,itlin,ibarcl
	real rlin
        real getpar
	logical bnohyd

c-------------------------------------------
c parameters
c-------------------------------------------

        ilin = nint(getpar('ilin'))
        rlin = getpar('rlin')
        itlin = nint(getpar('itlin'))
        ibarcl = nint(getpar('ibarcl'))
        bbarcl = ibarcl .gt. 0 .and. ibarcl .ne. 3
	call nonhydro_get_flag(bnohyd)
	
c-------------------------------------------
c initialization
c-------------------------------------------

	fxv = 0.
	fyv = 0.

c-------------------------------------------
c fix or nudge boundary velocities
c-------------------------------------------

	call bclfix

c-------------------------------------------
c horizontal diffusion
c-------------------------------------------

	call set_diff_horizontal

c-------------------------------------------
c advective (non-linear) terms
c-------------------------------------------

        if( ilin .eq. 0 ) then
          if( itlin .eq. 0 ) then
	    call set_advective(rlin)
	  else if( itlin .eq. 1 ) then
	    call set_semi_lagrange
	  else
	    write(6,*) 'itlin = ',itlin
	    stop 'error stop set_explicit: no such option'
	  end if
	end if

c-------------------------------------------
c baroclinic contribution
c-------------------------------------------

        !if( bbarcl ) call set_barocl
        !if( bbarcl ) call set_barocl_new
	if( bbarcl ) call set_barocl_new_interface

c-------------------------------------------
c non-hydrostatic contribution (experimental)
c-------------------------------------------

	if( bnohyd ) call nonhydro_set_explicit

c-------------------------------------------
c nudging of water levels and velocities
c-------------------------------------------

	call set_nudging

c-------------------------------------------
c end of routine
c-------------------------------------------

	end

c******************************************************************

	subroutine momentum_viscous_stability(ahpar,rindex,dstab)

c computes stability for viscosity
c
c stability is computed for dt == 1

	use mod_geom
	use mod_internal
	use mod_diff_visc_fric
	use evgeom
	use levels
	use basin

	implicit none

	real ahpar
	real rindex
	real dstab(nlvdi,nel)

	integer ie,ii,iei,l,lmax,k
	real u,v,ui,vi
	real anu,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact,r

	rindex = 0.
	if( ahpar .le. 0 ) return

	amax = 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax

	    a = 0.
	    r = 0.
	    do ii=1,3
              iei = ieltv(ii,ie)
	      k = nen3v(ii,ie)
              if( iei .le. 0 ) iei = ie

              areai = 12. * ev(10,iei)

	      r = r + rdistv(k)
	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai
	    end do
	    r = r/3.

	    a = a * r
	    amax = max(amax,a)
	    dstab(l,ie) = a

          end do

	end do

	rindex = amax

	end

c******************************************************************

	subroutine set_diff_horizontal

	use mod_geom
	use mod_internal
	use mod_diff_visc_fric
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,ii,iei,l,lmax
	integer noslip
	real u,v,ui,vi
	real anu,ahpar,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact
	logical bnoslip

	real getpar

	call get_timestep(dt)

        ahpar = getpar('ahpar')
	if( ahpar .le. 0 ) return

        noslip = nint(getpar('noslip'))
	bnoslip = noslip .ne. 0

	amax = 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax
	    u  = utlov(l,ie)
	    v  = vtlov(l,ie)

	    a = 0.
	    do ii=1,3
              iei = ieltv(ii,ie)
              afact = 1.
              if( bnoslip .and. iei .eq. 0 ) afact = -1.
              if( iei .le. 0 ) iei = ie

              areai = 12. * ev(10,iei)

	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai

	      ui = afact * utlov(l,iei)
	      vi = afact * vtlov(l,iei)

	      ax = ai * ( ui - u )
	      ay = ai * ( vi - v )

	      fxv(l,ie) = fxv(l,ie) - ax	!- because f is on left side
	      fyv(l,ie) = fyv(l,ie) - ay
	    end do
	    amax = max(amax,a)
          end do

	end do

	amax = amax * dt
	!write(99,*) 'stability viscosity: ',amax

	end

c******************************************************************

        subroutine set_momentum_flux

c sets arrays momentx/yv

	use mod_layer_thickness
	use mod_hydro_print
	use mod_hydro
	use mod_internal
	use evgeom
	use levels
	use basin

        implicit none

        integer ii,ie,k,l,lmax
        real b,c
        real ut,vt
        real uc,vc
        real up,vp
        real um,vm
        real f,h
	real xadv,yadv
	real area,vol

	real saux(nlvdi,nkn)

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	saux = 0.
	momentxv = 0.
	momentyv = 0.

c---------------------------------------------------------------
c accumulate momentum that flows into nodes (weighted by flux)
c---------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
            h = hdeov(l,ie)
	    ut = utlov(l,ie)
	    vt = vtlov(l,ie)
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                f = ut * b + vt * c	! f>0 => flux into node
                if( f .gt. 0. ) then
		  saux(l,k) = saux(l,k) + f
		  momentxv(l,k) = momentxv(l,k) + f * ut
		  momentyv(l,k) = momentyv(l,k) + f * vt
                end if
	    end do
          end do
	end do

c---------------------------------------------------------------
c compute average momentum for every node
c---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
            h = hdkov(l,k)
	    if( saux(l,k) .gt. 0 ) then		!flux into node
	      momentxv(l,k) = momentxv(l,k) / saux(l,k)
	      momentyv(l,k) = momentyv(l,k) / saux(l,k)
	    else				!only flux out of node
	      momentxv(l,k) = uprv(l,k) * h
	      momentyv(l,k) = vprv(l,k) * h
	    end if
	  end do
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c******************************************************************

        subroutine set_advective(rlin)

	use mod_internal
	use mod_layer_thickness
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

        implicit none

	real rlin		!strength of advection terms - normally 1

	real zxadv,zyadv
	real wtop,wbot

        integer ihwadv  	! vertical advection for momentum
        integer ii,ie,k,l,lmax
        real b,c
        real ut,vt
        real uc,vc
        real up,vp
        real um,vm
        real f,h
	real xadv,yadv
	real area,vol
        real getpar
        real wlay,dzbb,dz,dztt,ubot,utop,vbot,vtop

	!write(6,*) 'set_advective called...'

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

        ihwadv = nint(getpar('ihwadv'))

c---------------------------------------------------------------
c accumulate momentum that flows into nodes (weighted by flux)
c---------------------------------------------------------------

        call set_momentum_flux	!sets aux arrays momentx/yv

c---------------------------------------------------------------
c compute advective contribution
c---------------------------------------------------------------

	do ie=1,nel
          wtop = 0.0
	  lmax = ilhv(ie)
	  do l=1,lmax

c	    ---------------------------------------------------------------
c	    horizontal advection
c	    ---------------------------------------------------------------

	    area = 12. * ev(10,ie)
            h = hdeov(l,ie)
	    vol = area * h
  	    ut = utlov(l,ie)
  	    vt = vtlov(l,ie)
            uc = ut / h
            vc = vt / h

	    xadv = 0.
	    yadv = 0.
	    wbot = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
		wbot = wbot + wlov(l,k)
                up = momentxv(l,k) / h		!NEW
                vp = momentyv(l,k) / h
                f = ut * b + vt * c
                if( f .lt. 0. ) then	!flux out of node => into element
                  xadv = xadv + f * ( up - uc )
                  yadv = yadv + f * ( vp - vc )
                end if
            end do
	    
	    zxadv = 0.
	    zyadv = 0. 
	   
c	    ---------------------------------------------------------------
c	    vertical advection
c	    ---------------------------------------------------------------

	    if( ihwadv > 0 ) then	!compute vertical momentum advection
	      wbot = wbot / 3.
	      if( l .eq. lmax ) wbot = 0.

              if(ihwadv == 1) then	!use upwind scheme
  	        if(wtop.gt.0.) then
    	          zxadv = wtop * ulov(l,ie)
	          zyadv = wtop * vlov(l,ie)
                else
	          zxadv = wtop * ulov(l-1,ie)
	          zyadv = wtop * vlov(l-1,ie)
                end if

	        if(wbot.gt.0.) then
	          zxadv = zxadv - wbot * ulov(l+1,ie)
	          zyadv = zyadv - wbot * vlov(l+1,ie)
                else
	          zxadv = zxadv - wbot * ulov(l,ie)
                  zyadv = zyadv - wbot * vlov(l,ie)
                end if
              else if(ihwadv == 2) then	!use centered scheme
                dz = hdeov(l,ie)
                if (l .eq. 1) then
                  dzbb = hdeov(l+1,ie)
                  utop = 0.0
                  ubot = (ulov(l,ie)*dz+ulov(l+1,ie)*dzbb)/(dz+dzbb)
                  vtop = 0.0
                  vbot = (vlov(l,ie)*dz+vlov(l+1,ie)*dzbb)/(dz+dzbb)
                else if (l .eq. lmax) then
                  dztt = hdeov(l-1,ie)
                  utop = (ulov(l-1,ie)*dztt+ulov(l,ie)*dz)/(dztt+dz)
                  ubot = 0.0
                  vtop = (vlov(l-1,ie)*dztt+vlov(l,ie)*dz)/(dztt+dz)
                  vbot = 0.0
                else
                  dztt = hdeov(l-1,ie)
                  dzbb = hdeov(l+1,ie)
                  utop = (ulov(l-1,ie)*dztt+ulov(l,ie)*dz)/(dztt+dz)
                  ubot = (ulov(l,ie)*dz+ulov(l+1,ie)*dzbb)/(dz+dzbb)
                  vtop = (vlov(l-1,ie)*dztt+vlov(l,ie)*dz)/(dztt+dz)
                  vbot = (vlov(l,ie)*dz+vlov(l+1,ie)*dzbb)/(dz+dzbb)
                end if
                wlay = (wtop + wbot)/2.0
                zxadv = zxadv + wlay * (utop - ubot)
                zyadv = zyadv + wlay * (vtop - vbot)
	      else
		write(6,*) 'ihwadv = ',ihwadv
		stop 'error stop set_advective: ihwadv not supported'
              end if
	      wtop = wbot
	    end if

c	    ---------------------------------------------------------------
c	    total contribution
c	    ---------------------------------------------------------------

	    fxv(l,ie) = fxv(l,ie) + rlin * (xadv + zxadv)
	    fyv(l,ie) = fyv(l,ie) + rlin * (yadv + zyadv)
	  end do
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

        end

c******************************************************************

	subroutine momentum_advective_stability(rlin,rindex,astab)

c computes courant number of advective terms in momentum equation
c
c stability is computed for dt == 1

	use mod_internal
	use mod_geom_dynamic
	use mod_layer_thickness
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	real rlin		   !factor for advective terms - normally 1
	real rindex		   !stability index (return)
	real astab(nlvdi,nel)      !stability matrix (return)

	integer ie,l,ii,k,lmax,iweg
	real cc,cmax
	real ut,vt
	real area,h,vol
	real b,c,f,ftot,r

	cmax = 0.
	!call compute_stability_stats(-1,cc)

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  lmax = ilhv(ie)
	  iweg = iwegv(ie)
	  do l=1,lmax

            h = hdenv(l,ie)
	    vol = area * h

  	    ut = utlnv(l,ie)
  	    vt = vtlnv(l,ie)

	    ftot = 0.
	    r = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
		r = r + rdistv(k)
                f = ut * b + vt * c
                if( f .lt. 0. ) ftot = ftot - f
            end do
	    r = r/3.

	    cc = rlin*r*area*ftot/vol
	    if( iweg .gt. 0 ) cc = 0.	! dry element
	    astab(l,ie) = cc
	    cmax = max(cmax,cc)
	    !call compute_stability_stats(0,cc)

	  end do
	end do

	rindex = cmax
	!call compute_stability_stats(1,cc)

	end

c******************************************************************

	subroutine compute_stability_stats(what,cc)

c computes histogram of stability of elements

	implicit none

	integer what
	real cc

	integer ndim,dbin
	parameter( ndim = 10 , dbin = 10 )
	real eps
	parameter( eps = 1.e-5 )

	integer i,idt,it
	integer bin(0:ndim)
	save bin

	if( what .lt. 0 ) then		!initialize
	  do i=0,ndim
	    bin(i) = 0
	  end do
	else if( what .eq. 0 ) then	!accumulate
	  if( cc .gt. eps ) then
	    idt = 1. / cc
	  else
	    idt = 9999999
	  end if
	  i = idt / dbin
	  i = max(i,0)
	  i = min(i,ndim)
	  bin(i) = bin(i) + 1
	else				!write out
	  call get_act_time(it)
	  write(97,1000) it,(bin(i),i=0,ndim)
	end if
	
	return
 1000	format(i10,11i6)
	end

c******************************************************************

	subroutine set_semi_lagrange

	use mod_internal
	use basin, only : nkn,nel,ngr,mbw

        implicit none
         
	integer ie,l
	real xadv,yadv,dt
        real uadv(nel),vadv(nel)

	call get_timestep(dt)

        call back_trace(uadv,vadv)

	l = 1			!only for one layer

	do ie=1,nel
          xadv = uadv(ie) / dt
          yadv = vadv(ie) / dt

	  fxv(l,ie) = fxv(l,ie) + xadv
	  fyv(l,ie) = fyv(l,ie) + yadv
	end do

	end

c******************************************************************

        subroutine set_barocl

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin

        implicit none
         
	include 'pkonst.h'
        !integer itanf,itend,idt,nits,niter,it
        !real k,l,ie,ii				!BUG
        integer k,l,ie,ii			!BUG
        real dt
        real rrho0
        real salref,temref,sstrat,tstrat
        real hlayer
        real hhi

        real xbcl,ybcl
        integer lmax

        real rhop,presbt,presbcx,presbcy,dprescx,dprescy,br,cr!deb
        real b,c

        call get_timestep(dt)

        do ie=1,nel
                do l=1,nlv
                        bpresxv(l,ie) = 0.
                        bpresyv(l,ie) = 0.
                enddo
        enddo

        do ie=1,nel
            presbcx = 0.
            presbcy = 0.
            lmax=ilhv(ie)
            !print*,lmax,' lmax'
            do l=1,lmax            
                hlayer = 0.5 * hdeov(l,ie)
                
		br = 0.
		cr = 0.                 
                do ii=1,3                 
                        k = nen3v(ii,ie)
                        rhop = rhov(l,k) ! rho^prime for each node of element 
                        !print*,'rhov ', l,k,rhov(l,k)
                        b = ev(3+ii,ie)!gradient in x della funz di forma
                        c = ev(6+ii,ie)!gradient in y della funz di forma
                        br = br + (b*rhop) 
                        cr = cr + (c*rhop)
                end do
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
                bpresxv(l,ie) = presbcx
                bpresyv(l,ie) = presbcy
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
           end do
        end do
        
        rrho0=1./rowass
        !print*,'rowass ', rowass, grav,hhi 
        do ie=1,nel
            lmax=ilhv(ie)
                do l=1,lmax
                     hhi = hdeov(l,ie)
                     !print*, 'hhi ',hhi,bpresxv(l,ie),l,ie
                     xbcl =  rrho0*grav*hhi*bpresxv(l,ie)
                     ybcl =  rrho0*grav*hhi*bpresyv(l,ie)
                     fxv(l,ie) = fxv(l,ie) +xbcl
                     fyv(l,ie) = fyv(l,ie) +ybcl
                     !write(6,*)'fxv ',fxv,ie
                enddo
        enddo

        end

c**********************************************************************

        subroutine set_barocl_new

c computes baroclinic contribution centered on layers
c
c cannot use this for sigma levels

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin

        implicit none
         
	include 'pkonst.h'

	logical bsigma
        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

        raux=grav/rowass
	call get_bsigma(bsigma)

	if( bsigma ) then
	  stop 'error stop set_barocl_new: cannot use with sigma levels'
	end if

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)
          do l=1,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer

            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
          end do
        end do
        
        end

c**********************************************************************

        subroutine set_barocl_old

c do not use this routine !

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin

        implicit none
         
	include 'pkonst.h'

        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

	double precision px(0:nlvdi)
	double precision py(0:nlvdi)

	stop 'error stop set_barocl_old: do not use this routine'

        raux=grav/rowass

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)

	  px(0) = presbcx
	  py(0) = presbcy
          do l=1,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
            hlayer = hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
	    px(l) = presbcx
	    py(l) = presbcy

          end do

          do l=1,lmax
	    presbcx = 0.5*(px(l) + px(l-1))
	    presbcy = 0.5*(py(l) + py(l-1))
            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl
	  end do
        end do
        
        end

c**********************************************************************

        subroutine set_barocl_new_interface

c computes baroclinic contribution centered on interfaces
c
c this routine works with Z and sigma layers

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use mod_hydro
	use evgeom
	use levels
	use basin

        implicit none
         
	include 'pkonst.h'

c---------- DEB SIG
	real hkk
	!real hkko(0:nlvdi,nkn)	!depth of interface at node
	!real hkkom(0:nlvdi,nkn)	!average depth of layer at node
	real, allocatable :: hkko(:,:)	!depth of interface at node
	real, allocatable :: hkkom(:,:)	!average depth of layer at node
	real hele
	real helei
	real alpha,aux,bn,cn,bt,ct
	real h,hd,hu,rd,ru
	real brl,crl

	integer lkmax,ld,lu
	integer laux,ll,ls,nn,nb !DEB
	integer llup(3),lldown(3)
c---------- DEB SIG

	logical bsigma,bsigadjust
        integer k,l,ie,ii,lmax,lmin,nsigma
	real hsigma,hdep
        double precision hlayer,hint,hhk,hh,hhup,htint
	double precision dzdx,dzdy,zk
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr,brup,crup,brint,crint
	double precision rhoup,psigma
	double precision b3,c3

	bsigadjust = .false.		!regular sigma coordinates
	bsigadjust = .true.		!interpolate on horizontal surfaces

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

        raux=grav/rowass
	psigma = 0.
	
	allocate(hkko(0:nlvdi,nkn))
	allocate(hkkom(0:nlvdi,nkn))

	if( bsigma .and. bsigadjust ) then	!-------------- DEB SIG
	  do k=1,nkn
	    lmax=ilhkv(k)
	    hkko(0,k)=-zov(k)	!depth of interface on node
	    hkkom(0,k)=-zov(k)	!depth of mid layer on node (0 not used)
	    hkk=0.
	    hkk=-zov(k)		!ggu
	    do l=1,lmax
	      hkk=hkk+hdkov(l,k)
	      hkko(l,k)=hkk
	      hkkom(l,k)=(hkko(l,k)+hkko(l-1,k))/2.
            end do
	  end do
	end if
	 
        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)
	  brup=0.
	  crup=0.
	  hhup=0.
          do l=1,lmax		!loop over layers to set up interface l-1
	    bsigma = l .le. nsigma

	    htint = 0.				!depth of layer top interface
	    if( l .gt. 1 ) htint = hlv(l-1)

            hlayer = hdeov(l,ie)		!layer thickness
	    if( .not. bsigma ) hlayer = hldv(l)

            hh = 0.5 * hlayer
	    hint = hh + hhup			!interface thickness
                
	    if( bsigma .and. bsigadjust ) then	!-------------- DEB SIG
	      hele = 0.
	      helei = 0.
	      do ii=1,3
                k = nen3v(ii,ie)
	        hele=hele+hkko(l,k)+hkko(l-1,k)	!depth of mid layer in element
	        helei=helei+hkko(l-1,k)		!depth of interface in element
	      end do

	      hele=hele/6.			!depth of mid layer l
	      helei=helei/3.			!depth of interface l-1

	      do ii=1,3
                k = nen3v(ii,ie)   
		lkmax = ilhkv(k)
	        if(helei.lt.hkko(l-1,k))then	!look upwards
		  do ll=l-1,1,-1
	            if(helei.gt.hkko(ll-1,k)) exit
		  end do
		  if( ll .le. 0 ) ll = 1
                else if(helei.gt.hkko(l,k))then	!look downwards
		  do ll=l+1,lkmax
	            if(helei.lt.hkko(ll,k)) exit
		  end do
		  if( ll .gt. lkmax ) ll = lkmax
		else				!inside layer l
		  ll = l
	        end if
		!interface l-1 is inside layer ll
		if( helei.lt.hkkom(ll,k) ) then	!find part of layer (up or down)
		  llup(ii) = ll-1
		  if( ll .eq. 1 ) llup(ii) = 1
		  lldown(ii) = ll
		else
		  llup(ii) = ll
		  lldown(ii) = ll+1
		  if( ll .eq. lkmax ) lldown(ii) = lkmax
		end if
		!do final check just to be sure (may be commented)
		if( lldown(ii) .eq. 1 ) then
		  if( helei .gt. hkkom(1,k) ) goto 99
		else if( llup(ii) .eq. lkmax ) then
		  if( helei .lt. hkkom(lkmax,k) ) goto 99
		else
		  if( helei .gt. hkkom(lldown(ii),k) ) goto 99
		  if( helei .lt. hkkom(llup(ii),k) ) goto 99
		end if
	      end do
	    end if

	    nn = 0 
	    nb = 0
	    brl = 0.
	    crl = 0.                 
	    br = 0.
	    cr = 0.                 
	    dzdx = 0.
	    dzdy = 0.
	    psigma = 0.

            do ii=1,3
	      k = nen3v(ii,ie)
	      rhop = rhov(l,k)		!rho^prime for each node of element
	      rhoup = rhop
	      if( l.gt.1) rhoup = rhov(l-1,k)
	      lkmax = ilhkv(k)
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y

	      if( l .eq. nsigma ) then	!last sigma layer
	        brl = brl + b * rhop
	        crl = crl + c * rhop
	      end if

	      if( bsigma .and. bsigadjust ) then 
		lu = llup(ii)
		ld = lldown(ii)
		if( ld .eq. 1 ) then		!above surface
		  rhop = rhov(1,k)
		else if( lu .ge. lkmax ) then	!below bottom
		  nb = nb + 1
		  nn = nn + ii
		  rhop = rhov(lkmax,k)
		else				!do interpolation
		  !hu = hkko(lu,k)
		  !hd = hkko(ld,k)
		  hu = hkkom(lu,k) !DEB
		  hd = hkkom(ld,k) !DEB
		  ru = rhov(lu,k)
		  rd = rhov(ld,k)
		  h = helei
		  alpha = (h-hu)/(hd-hu)
		  rhop = alpha*rd + (1.-alpha)*ru
		end if
	      end if

              br = br + b * rhop
              cr = cr + c * rhop

              if (bsigma) then
	       if( bsigadjust ) then
		psigma = 0.
	       else
                psigma = psigma + (rhoup-rhop)/hint
                hdep = hm3v(ii,ie) + zov(k)
                hhk = -htint * hdep
                zk = -hhk               !transform depth in z
                dzdx = dzdx + b * zk
                dzdy = dzdy + c * zk
	       end if
              end if
            end do

	    if( bsigma .and. bsigadjust ) then 
              if(nb.eq.2)then
	        brint = brup
	        crint = crup
	      elseif(nb.eq.1)then
	        b3 = ev(3+nn,ie)
	        c3 = ev(6+nn,ie)
	        aux=1./(c3*c3+b3*b3)
	        bn = aux*(brup*b3+crup*c3)*b3
	        cn = aux*(brup*b3+crup*c3)*c3
	        bt = br - aux*(br*b3+cr*c3)*b3
	        ct = cr - aux*(br*b3+cr*c3)*c3
	        brint = bn + bt
	        crint = cn + ct
              else  
	        brint = br
	        crint = cr
	      end if
	    else			  !zeta layer
              if( l .eq. 1 ) then         !surface layer ... treat differently
                brint = br
                crint = cr
              else
                brint = 0.5*(br+brup)
                crint = 0.5*(cr+crup)
              end if
	    end if

	    brup=br
	    crup=cr
	    if( l .eq. nsigma ) then
	      brup=brl
	      crup=crl
	    end if
            hhup=hh
	    psigma = psigma / 3.

            presbcx = presbcx + hint * ( brint - dzdx * psigma )
	    presbcy = presbcy + hint * ( crint - dzdy * psigma )

            xbcl =  raux * hlayer * presbcx
            ybcl =  raux * hlayer * presbcy

            fxv(l,ie) = fxv(l,ie) + xbcl 
            fyv(l,ie) = fyv(l,ie) + ybcl
          end do
        end do
        
	deallocate(hkko)
	deallocate(hkkom)

	return
   99	continue
	write(6,*) ie,k,ii
	write(6,*) llup(ii),lldown(ii)
	write(6,*) helei,hkkom(lldown(ii),k),hkkom(llup(ii),k)
	stop 'error stop set_barocl_new_interface: internal error'
        end

c**********************************************************************

