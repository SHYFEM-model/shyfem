c
c $Id: lagrange_flux.f,v 1.3 2009-05-21 09:24:00 georg Exp $
c
c routines handling fluxes
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c 28.04.2009    ggu     links re-structured
c 20.10.2011    ggu     new routines for fluxes implemented
c 24.10.2011    ggu     3d routines implemented
c 04.11.2011    ggu     part with bsigma not yet finished
c 16.12.2011    ggu     bug fix - flux2d_aux was integer
c 23.03.2012    ggu     bug fix - dst was not available in setup_vl_3d
c 28.08.2012    ggu&ccf bug fix - ok now for 3d surface tracking
c
c****************************************************************

	subroutine lagr_setup_timestep

c initializes length of element sides and fluxes

	implicit none

	call rprs
	call setup_fluxes_2d
	call setup_fluxes_3d	!this overwrites velocities in vel_ie (2d)

	end

c****************************************************************

	subroutine rprs

c initializes length of element sides

	use basin

	implicit none
	
	include 'param.h'
	include 'lagrange.h'

	


	real x1,x2,y1,y2,dx,dy,ddx,ddy,d
	integer ie,k,i1,i2,p1,p2

	integer icall
	save icall
	data icall / 0 /

	if( icall .ne. 0 ) return

	icall = 1

	do ie=1,nel
	 do k=1,3
          i1=mod(k,3)+1
          i2=mod(i1,3)+1
	  p1=nen3v(i1,ie)
	  p2=nen3v(i2,ie)
	  x1=xgv(p1)
	  x2=xgv(p2)
	  y1=ygv(p1)
	  y2=ygv(p2)
	  dx=x2-x1
	  dy=y2-y1
	  ddx=dx**2
	  ddy=dy**2
	  d=sqrt(ddx+ddy)
	  dvert(k,ie)=d
	 end do
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine setup_fluxes_3d

c sets up fluxes in 3d - has to be done every time step

	use mod_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'lagrange.h'
	


	integer k,ie,ii,i
	integer nn,ne
	integer n
	integer l,lmax,lkmax
	integer itype
	real az,azpar
	real tdif
	
	real rflux(nlvdi,ngr)       !fluxes across finite volume k
	real tflux(nlvdi,ngr)       !fluxes across sides of element
	
	integer flxtype

c	--------------------------------------------
c	initialization
c	--------------------------------------------

	call getaz(azpar)
	az = azpar
	azlgr = az		!store in include file

        do ie=1,nel
	  lmax = ilhv(ie)
          do ii=1,3
	    do l=1,lmax
              flux3d(l,ii,ie)=0
	    end do
          end do
        end do

c	--------------------------------------------
c	loop on nodes
c	--------------------------------------------

        do k=1,nkn
          itype=flxtype(k)

	  n = ngr
	  call flx3d_k(k,itype,az,lkmax,n,rflux)
	  call make_fluxes_3d(k,itype,lkmax,n,rflux,tflux)

  	  call setup_flux3d(k,lkmax,n,tflux)
        end do

c	--------------------------------------------
c	compute velocities
c	--------------------------------------------

        call setup_vl_3d

c	--------------------------------------------
c	end of routine
c	--------------------------------------------

        end

c******************************************************************

	subroutine setup_fluxes_2d

c sets up fluxes - has to be done every time step

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'lagrange.h'
	

	integer k,ie,ii,i
	integer nn,ne
	integer n,j
	integer itype
	real az,azpar
	real tdif
	
	real rflux(ngr)       !fluxes across finite volume k
	real tflux(ngr)       !fluxes across sides of element
	real tflux_aux(ngr)   !fluxes across sides of element (aux)
	real flux2d_aux(3,nel)
	
	integer flxtype

        do ie=1,nel
          do ii=1,3
            flux2d(ii,ie)=0
            flux2d_aux(ii,ie)=0
          end do
        end do

	call getaz(azpar)
	az = azpar

        do k=1,nkn
          itype=flxtype(k)

c	  --------------------------------------------
c	  new way to set-up rflux,tflux
c	  --------------------------------------------

	  n = ngr
	  call flx2d_k(k,itype,az,n,rflux)
	  call make_fluxes_2d(k,itype,n,rflux,tflux)

c	  --------------------------------------------
c	  error check ... delete
c	  --------------------------------------------

c	  --------------------------------------------
c	  set-up lagrangian fluxes
c	  --------------------------------------------

  	  call setup_flux2d(k,n,tflux,flux2d)		!new
        end do

c	--------------------------------------------
c	error check ... delete
c	--------------------------------------------

c	--------------------------------------------
c	compute velocities
c	--------------------------------------------

        call setup_vl_2d

c	--------------------------------------------
c	end of routine
c	--------------------------------------------

        end

c******************************************************************
c******************************************************************
c******************************************************************
c old routines -> delete
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine mk_rflux(k,n,itype,az,transp,ne,elems)

c computes flux through finite volume k
c
c internal section is defined by:  kbefor - k - kafter

	use mod_hydro
	use evgeom

	implicit none

        include 'param.h'

	integer k,n,itype,ne
	real az
	real transp(1)
	integer elems(1)

	include 'femtime.h'

        !real uov(1),vov(1),unv(1),vnv(1)
        !common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv

	integer i,ip,ie,ii
	integer l
	real aj,area,dz,uv,rdt
	real uvn,uvo
	real b,c
	real azt,tt

	integer ithis
	 
	rdt = 1./idt
	azt = 1. - az

c compute transports into finite volume of node k -> transp
c computed transports are divergence corrected

	l = 1
	do i=1,ne
	 ie = elems(i)
	 ii = ithis(k,ie)
	 aj = ev(10,ie)
	 area = 4. * aj
	 dz = zenv(ii,ie) - zeov(ii,ie)
	 b = ev(3+ii,ie)
	 c = ev(6+ii,ie)
	 uvn = utlnv(l,ie) * b + vtlnv(l,ie) * c
	 uvo = utlov(l,ie) * b + vtlov(l,ie) * c
	 uv = 12. * aj * ( az*uvn + azt*uvo )
	 uv = uv - dz * area * rdt
	 transp(i) = uv
	end do

	end

c******************************************************************

        subroutine mk_tflux(k,n,itype,rflux,tflux)

c computes fluxes over sides (tflux) from fluxes into node (rflux)

        implicit none

        integer k               !node
        integer n               !number of sides (tfluxes)
        integer itype           !type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
        real rflux(n)        !fluxes into node (element)
        real tflux(n)        !fluxes through sides (return value)

	integer ipf,ipl      !puntatori vertici

        integer i
        real rr

	if( itype .eq. 1 ) then         !internal node
                rr = 0.
                do i=1,n-1
                  rr = rr + i * rflux(i)
                end do
                rr = rr / n

                tflux(n) = rr
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do
        else if( itype .eq. 2 ) then    !node on material boundary
                tflux(1) = 0.
                do i=2,n-1
                  tflux(i) = tflux(i-1) + rflux(i-1)
                end do
                tflux(n) = 0.
        else if( itype .eq. 3 ) then    !BOO - boundary on left
                tflux(n) = 0.
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do
        else if( itype .eq. 4 ) then    !OOB - boundary on right
                tflux(1) = 0.
                do i=2,n
                  tflux(i) = tflux(i-1) + rflux(i-1)
                end do
        else if( itype .eq. 5 ) then    !node totaly on open boundary
                rr = 0.
                do i=1,n-1
                  rr = rr + i * rflux(i)
                end do
                rr = rr / n
                tflux(n) = rr
                do i=n-1,1,-1
                  tflux(i) = tflux(i+1) - rflux(i)
                end do

        else
                stop 'error stop make_fluxes: internal error (1)'
        end if
	end

c**********************************************************************

        subroutine setup_fx(k,n,tflux,ne,elems)

c computes fluxes in element

        implicit none

        integer k,n,ne
        real tflux(1)
	integer elems(1)

        include 'param.h'
	include 'lagrange.h'

        integer i,ie,ip,ii
        integer j,n1,n2

        integer ithis
        integer ibhnd
        integer inext

	if( n .ne. ne ) then	!just a check
	  if( ne+1 .ne. n ) then
	    write(6,*) k,n,ne,(elems(i),i=1,n)
	    stop 'error stop setup_fx: internal error (1)'
	  end if
	  if( elems(n) .ne. 0 ) then
	    write(6,*) k,n,ne,(elems(i),i=1,n)
	    stop 'error stop setup_fx: internal error (2)'
	  end if
	end if

        do i=1,ne
         ie=elems(i)
         n1=ibhnd(k,ie)
         n2=inext(k,ie)
         j=i+1
         if(i.eq.n) j=1
	!write(6,*) 'setup_fx ',k,n,ne
	!write(62,*) k,n,ne
	!write(62,*) i,ie,n1,n2,j
	!write(62,*) tflux(j),tflux(i)
         flux2d(n2,ie)=flux2d(n2,ie)-tflux(j)
         flux2d(n1,ie)=flux2d(n1,ie)+tflux(i)      
        end do

        end

c******************************************************************
c******************************************************************
c******************************************************************
c end old routines
c******************************************************************
c******************************************************************
c******************************************************************

  	subroutine setup_flux3d(k,lkmax,n,tflux)

c computes 3d fluxes in element
c
c maybe not existing fluxes are set in flux3d
c we do not use lkmax, but lmax

	use mod_geom
	use levels

        implicit none

        include 'param.h'
	include 'lagrange.h'

        integer k,lkmax,n
	real tflux(nlvdi,1)       !fluxes across sides of element


        integer i,ie,ne,l,lmax
        integer j,n1,n2

        integer ibhnd
        integer inext

	ne = n
	if( lnk_elems(ne) .le. 0 ) ne = ne - 1

        do i=1,ne
          ie=lnk_elems(i)	!set up outside
          n1=ibhnd(k,ie)
          n2=inext(k,ie)
          j=i+1
          if(i.eq.n) j=1
	  lmax = ilhv(ie)
	  do l=1,lmax
            flux3d(l,n2,ie)=flux3d(l,n2,ie)-tflux(l,j)
            flux3d(l,n1,ie)=flux3d(l,n1,ie)+tflux(l,i)      
	  end do
        end do

	end

c******************************************************************

  	subroutine setup_flux2d(k,n,tflux,flux2d_loc)

c computes fluxes in element

	use mod_geom

        implicit none

        integer k,n
        real tflux(1)
	real flux2d_loc(3,1)

        include 'param.h'
	include 'lagrange.h'

        integer i,ie,ne
        integer j,n1,n2

        integer ibhnd
        integer inext

	ne = n
	if( lnk_elems(ne) .le. 0 ) ne = ne - 1

        do i=1,ne
          ie=lnk_elems(i)	!set up outside
          n1=ibhnd(k,ie)
          n2=inext(k,ie)
          j=i+1
          if(i.eq.n) j=1
	!write(6,*) 'setup_flux2d ',k,n,ne
	!write(61,*) k,n,ne
	!write(61,*) i,ie,n1,n2,j
	!write(61,*) tflux(j),tflux(i)
          !flux2d_aux(n2,ie)=flux2d_aux(n2,ie)-tflux(j)
          !flux2d_aux(n1,ie)=flux2d_aux(n1,ie)+tflux(i)      
          flux2d_loc(n2,ie)=flux2d_loc(n2,ie)-tflux(j)
          flux2d_loc(n1,ie)=flux2d_loc(n1,ie)+tflux(i)      
        end do

	end

c*********************************************************************

        subroutine setup_vl_3d
	
c computes velocities in element

c use first layer thickness which already contains water level
c
c dH = dh + dz		dh undisturbed depth  dz average water level
c
c dz = (z1+z2+z3)/3.
c
c dp = dh + (z1+z2)/2 = dH - dz + (z1+z2)/2 
c		= dH + (3*z1+3*z2-2*z1-2*z2-2*z3)/6
c		= dH + (z1+z2-2*z3)/6
c
c all this has to be revised for sigma layers

	use mod_layer_thickness
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'lagrange.h'






	logical bsigma
	integer ie,ii,i1,i2
	integer l,lmax
	real flx,dh,zi1,zi2,zi3,dp,ar,dst
	real bfact
	
	call get_bsigma(bsigma)

	bfact = 1.
	if( bback ) bfact = -1.

	do ie=1,nel
	  dh=hdenv(1,ie)		!this already contains water level
	  lmax = ilhv(ie)
	  do l=1,lmax
            do ii=1,3
	      dp = dh
	      dst=dvert(ii,ie)
	      if( .not. bsigma .and. l .eq. 1 ) then	!z-levels, first layer
                i1=mod(ii,3)+1
                i2=mod(i1,3)+1
	        zi1=zenv(i1,ie)
	        zi2=zenv(i2,ie)
	        zi3=zenv(ii,ie)
	        dp=dh+((zi1+zi2-2.*zi3)/6.)	!compute aver depth between 1&2
	      end if
	      ar=dp*dst				!section area
	      flx=flux3d(l,ii,ie)		!flux
	      vel3d_ie(l,ii,ie)=bfact*flx/ar	!velocity
	      !if( l .eq. 1 ) then
	      !  vel_ie(ii,ie)=vel3d_ie(l,ii,ie) !velocity first layer
	      !end if
	    end do
          end do
	  do ii=1,3
	    vel_ie(ii,ie)=vel3d_ie(1,ii,ie)	!velocity first layer
	  end do
	end do
	
	end

c*********************************************************************

        subroutine setup_vl_2d
	
c computes velocities in element

c use first layer thickness which already contains water level
c
c dH = dh + dz		dh undisturbed depth  dz average water level
c
c dz = (z1+z2+z3)/3.
c
c dp = dh + (z1+z2)/2 = dH - dz + (z1+z2)/2 
c		= dH + (3*z1+3*z2-2*z1-2*z2-2*z3)/6
c		= dH + (z1+z2-2*z3)/6
c
c all this has to be revised for sigma layers

	use mod_layer_thickness
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'lagrange.h'




	integer ie,ii,i1,i2
	real flx,dh,zi1,zi2,zi3,dp,ar,dst
	real bfact,eps

	eps=1.e-5
	bfact = 1.
	if( bback ) bfact = -1.

	do ie=1,nel
	 dh=hdenv(1,ie)			!this already contains water level
         do ii=1,3
	  dst=dvert(ii,ie)
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
	  zi1=zenv(i1,ie)
	  zi2=zenv(i2,ie)
	  zi3=zenv(ii,ie)
	  dp=dh+((zi1+zi2-2.*zi3)/6.)	!compute aver depth between 1&2
	  ar=dp*dst			!section area
	  if(ar.le.eps) ar = eps	!area keep above 0
	  flx=flux2d(ii,ie)		!flux
	  vel_ie(ii,ie)=bfact*flx/ar	!velocity
         end do
	end do
	
	end

c**********************************************************************

