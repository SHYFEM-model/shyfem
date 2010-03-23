c
c $Id: lagrange_flux.f,v 1.3 2009-05-21 09:24:00 georg Exp $
c
c routines handling fluxes
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c 28.04.2009    ggu     links re-structured
c
c****************************************************************

	subroutine lagr_setup_timestep

c initializes length of element sides and fluxes

	implicit none

	call rprs
	call setup_fluxes

	end

c****************************************************************

	subroutine rprs

c initializes length of element sides

	implicit none
	
	include 'param.h'
	include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	
        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,1)
        common /nen3v/nen3v

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

	subroutine setup_fluxes

c sets up fluxes - has to be done every time step

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'links.h'
	
	integer ndim
	parameter (ndim=20)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,nn,ne,ipf,ipl,ie,ii
	integer itype,flxtype
	
	real az,azpar
	
	real rflux(ndim)       !fluxes across fifnite volume k
	real tflux(ndim)       !fluxes across sides of element
	
        do ie=1,nel
          do ii=1,3
            fx(ii,ie)=0
          end do
        end do

	call getaz(azpar)
	az = azpar

        do k=1,nkn
          itype=flxtype(k)
	  call set_elem_links(k,ne)
	  nn = ne
          if( itype .gt. 1 ) nn = nn + 1   !boundary
          if( nn .gt. ndim ) stop 'error stop flxnod: ndim'
	  rflux(nn) = 0.
	
    	  call mk_rflux(k,nn,itype,az,rflux,ne,lnk_elems)
  	  call mk_tflux(k,nn,itype,rflux,tflux)
  	  call setup_fx(k,nn,tflux,ne,lnk_elems)
        end do

        call setup_vl

        end

c******************************************************************

	subroutine mk_rflux(k,n,itype,az,transp,ne,elems)

c computes flux through finite volume k
c
c internal section is defined by:  kbefor - k - kafter

	implicit none

        include 'param.h'

	integer k,n,itype,ne
	real az
	real transp(1)
	integer elems(1)

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	include 'ev.h'
        real utlnv(nlvdim,1)
        common /utlnv/utlnv
        real vtlnv(nlvdim,1)
        common /vtlnv/vtlnv
        real utlov(nlvdim,1)
        common /utlov/utlov
        real vtlov(nlvdim,1)
        common /vtlov/vtlov
        !real uov(1),vov(1),unv(1),vnv(1)
        !common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov

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

        do i=1,ne
         ie=elems(i)
         n1=ibhnd(k,ie)
         n2=inext(k,ie)
         j=i+1
         if(i.eq.n) j=1
         fx(n2,ie)=fx(n2,ie)-tflux(j)
         fx(n1,ie)=fx(n1,ie)+tflux(i)      
        end do

        end

c*********************************************************************

        subroutine setup_vl
	
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

	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real zenv(3,1),zeov(3,1)
        common /zenv/zenv, /zeov/zeov

        !real hev(1)
        !common /hev/hev
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv

	integer ie,ii,i1,i2
	real flx,dh,zi1,zi2,zi3,dp,ar,dst
	real bfact
	
	bfact = 1.
	if( bback ) bfact = -1.

	do ie=1,nel
	 !dh=hev(ie)
	 dh=hdenv(1,ie)			!this already contains water level
         do ii=1,3
	  dst=dvert(ii,ie)
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
	  zi1=zenv(i1,ie)
	  zi2=zenv(i2,ie)
	  zi3=zenv(ii,ie)
	  dp=dh+((zi1+zi2-2.*zi3)/6.)	!see above
	  !dp=dh+((zi1+zi2)/2)
	  ar=dp*dst	
	  flx=fx(ii,ie)			!flux
	  vel_ie(ii,ie)=bfact*flx/ar	!velocity
         end do
	end do
	
	end

c**********************************************************************

