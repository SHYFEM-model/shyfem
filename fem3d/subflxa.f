c
c $Id: subflxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
c subroutines for computing discharge / flux
c
c contents :
c
c        subroutine inflxa
c        subroutine rdflxa
c        subroutine ckflxa
c        subroutine prflxa
c        subroutine tsflxa
c
c        subroutine wrflxa(it)			write of flux data
c
c        subroutine flxscs(kflux,iflux,n,az,flux,flxpm)	flux through sections
c        function flxsec(kflux,iflux,n,az)		flux through section
c        function flxnov(k,ibefor,iafter,istype,az)	flux through volume k
c        subroutine mkweig(n,istype,is,weight)	 	computes weight
c
c        subroutine flxini			initializes flux routines
c        subroutine flxin(kflux,iflux,n)	sets up info structure
c	 function flxtype(k)			determines type of node k (1-5)
c        subroutine flxinf(kflux,iflux,n)	sets up one info structure
c        function igtnsc(k1,k2)			gets number of internal section
c
c----------------------------------------------------------------
c the next subroutines shoud be replaced by subs in sublin.f
c        subroutine flxe2i(kflux,n,berror)	external to internal nodes
c        function kflxck(kflux,n)		checks compatibility of kflux
c        subroutine flxdbl(kflux,n,bstop)	tests for uniqueness
c        subroutine flxadj(kflux,n,bstop)	tests for adjacency
c        function extrsc(inode,ndim,nnode,ifirst,ilast)	extracts section
c----------------------------------------------------------------
c
c	 subroutine extrsect(isect,knode,nnode)	extracts nodes of one section
c
c revision log :
c
c 30.04.1998    ggu	newly written routines (subpor deleted)
c 07.05.1998    ggu	check nrdveci on return for error
c 08.05.1998    ggu	restructured with new comodity routines
c 13.09.1999    ggu	type of node computed in own routine flxtype
c 19.11.1999    ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006    ggu     changed time step to real
c 28.09.2007    ggu     use testbndo to determine boundary node in flxtype
c 28.04.2009    ggu     links re-structured
c 23.02.2011    ggu     new routine call write_node_fluxes() for special output
c 01.06.2011    ggu     documentation to flxscs() changed
c
c notes :
c
c These routines can also be used internally to compute the flux
c over various sections. The following calling sequence must be respected:
c
c call flxe2i(kflux,n,berror)		converts external to internal nodes
c nsect = kflxck(kflux,n)		checks array kflux and computes nsect
c call flxin(kflux,iflux,n)		initializes iflux
c
c call flxscs(kflux,iflux,n,az,flux,flxpm) computes fluxes and returns in flux()
c
c Initialization can be done anytime.
c
c******************************************************************

        subroutine mod_flx(mode)
 
        implicit none
 
        integer mode
 
        include 'modules.h'
 
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
 
        if( mode .eq. M_AFTER ) then
           call wrflxa(it)
        else if( mode .eq. M_INIT ) then
           call inflxa
        else if( mode .eq. M_READ ) then
           call rdflxa
        else if( mode .eq. M_CHECK ) then
           call ckflxa
        else if( mode .eq. M_SETUP ) then
           call flxini
        else if( mode .eq. M_PRINT ) then
           call prflxa
        else if( mode .eq. M_TEST ) then
           call tsflxa
        else if( mode .eq. M_BEFOR ) then
c          nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_flx'
        end if
 
        end

c******************************************************************


        subroutine inflxa

c nsect		total number of sections
c kfluxm	total number of nodes defining sections
c kflux()	node numbers defining sections

        implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux

        nsect = -1	!must still be initialized
        kfluxm = 0

        end

c******************************************************************

        subroutine rdflxa

        implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux

	integer nfxdim
        integer nrdveci

	call getdim('nfxdim',nfxdim)

        kfluxm = nrdveci(kflux,nfxdim)

        if( kfluxm .lt. 0 ) then
          if( kfluxm .eq. -1 ) then
            write(6,*) 'dimension error nfxdim in section $flux : '
     +                          ,nfxdim
          else
            write(6,*) 'read error in section $flux'
          end if
          stop 'error stop rdflxa'
        end if

        end

c******************************************************************

        subroutine ckflxa

        implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux

	integer k,ii
        logical berror

c	call flxe2i(kflux,kfluxm,berror)	!FIXME
	call n2int(kfluxm,kflux,berror)

        if( berror ) then
		write(6,*) 'error in section FLUX'
		stop 'error stop: ckflxa'
	end if

c initialize vectors (not strictly necessary)

	do k=1,kfluxm
	  do ii=1,3
	    iflux(ii,k) = 0
	  end do
	end do

c the real set up is done in flxini
c but, since at this stage we do not have all the arrays set up
c we post-pone it until later

        end

c******************************************************************

	subroutine prflxa

	implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	
	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii

	integer ipext
	logical nextline

	write(6,*)
	write(6,*) 'flux section :'
	write(6,*)
	write(6,*) 'nsect,kfluxm ',nsect,kfluxm
	write(6,*)

	ns = 0
	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'section : ',ns,ntotal
	  do i=ifirst,ilast
	    write(6,*) ipext(kflux(i)),(iflux(ii,i),ii=1,3)
	  end do
	end do

	end

c******************************************************************

	subroutine tsflxa

	implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	
	integer i,ii

	write(6,*) '/kfluxc/'
	write(6,*) nsect,kfluxm
	write(6,*) (kflux(i),i=1,kfluxm)

	write(6,*) '/iflux/'
	write(6,*) ((iflux(ii,i),ii=1,3),i=1,kfluxm)

	end

c******************************************************************

	subroutine wrflxa(it)

c write of flux data

	implicit none

	integer it

        integer iscdim
        parameter(iscdim=500)

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	save /kfluxc/,/iflux/	!ggu

	integer itend
	integer i
	real az,azpar,rr

	integer iround,ideffi
	real getpar

	real fluxt(iscdim)	!accumulator - could be also double precision
	real flux(iscdim)
	real fluxtpm(2,iscdim)	!accumulator - could be also double precision
	real fluxpm(2,iscdim)

        integer idtflx,itflx,itmflx,nr
        integer icall,nbflx,nvers,idfile
        save fluxt,fluxtpm
        save idtflx,itflx,itmflx,nr
        save icall,nbflx,nvers,idfile
        data icall,nbflx,nvers,idfile /0,0,2,537/

c start of code

        if( icall .eq. -1 ) return

c initialization

        if( icall .eq. 0 ) then

                idtflx = iround(getpar('idtflx'))
                itmflx = iround(getpar('itmflx'))
                itend = iround(getpar('itend'))

                if( kfluxm .le. 0 ) icall = -1
                if( nsect .le. 0 ) icall = -1
                if( idtflx .le. 0 ) icall = -1
                if( itmflx .gt. itend ) icall = -1
                if( icall .eq. -1 ) return

                if( nsect .gt. iscdim ) then
                  stop 'error stop wrflxa: dimension iscdim'
                end if

                itflx = itmflx + idtflx
		itmflx = itmflx + 1	!start from next time step

                nr = 0
                do i=1,nsect
                  fluxt(i) = 0.
          	  fluxtpm(1,i) = 0.
          	  fluxtpm(2,i) = 0.
                end do

                nbflx=ideffi('datdir','runnam','.flx','unform','new')
                if(nbflx.le.0) then
        	   stop 'error stop wrflxa : Cannot open FLX file'
		end if

                write(nbflx) idfile,nvers
                write(nbflx) nsect,kfluxm,idtflx
                write(nbflx) (kflux(i),i=1,kfluxm)

c               here we could also compute and write section in m**2

        end if

c normal call

        icall = icall + 1

        if( it .lt. itmflx ) return

c	accumulate results

        nr = nr + 1
	call getaz(azpar)
	az = azpar

	call flxscs(kflux,iflux,kfluxm,az,flux,fluxpm)

	do i=1,nsect
	  fluxt(i) = fluxt(i) + flux(i)
	  !write(99,*) 'wrflxa: ',it,i,flux(i),fluxt(i)  !gguerror
	  fluxtpm(1,i) = fluxtpm(1,i) + fluxpm(1,i)
	  fluxtpm(2,i) = fluxtpm(2,i) + fluxpm(2,i)
	end do

        if( it .lt. itflx ) return

c	write results

        itflx=itflx+idtflx

        rr=1./nr

        do i=1,nsect
          flux(i) = fluxt(i) * rr
          fluxpm(1,i) = fluxtpm(1,i) * rr
          fluxpm(2,i) = fluxtpm(2,i) * rr
        end do

        write(nbflx) it,nsect,(flux(i),i=1,nsect)
     +			,(fluxpm(1,i),fluxpm(2,i),i=1,nsect)
c	write(6,*) 'section written: ',nsect,nr,it

c	reset variables

        nr = 0
        do i=1,nsect
          fluxt(i) = 0.
          fluxtpm(1,i) = 0.
          fluxtpm(2,i) = 0.
        end do

	end

c******************************************************************

	subroutine flxscs(kflux,iflux,n,az,flux,fluxpm)

c computes flux through all sections and returns them in flux

	implicit none

	integer n
	integer kflux(1)
	integer iflux(3,1)
	real az
	real flux(1)
	real fluxpm(2,1)	!flux divided into positive (1) and negative (2)

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline
	real flxsec

	nnode = 0
	ns = 0

	do while( nextline(kflux,n,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  flux(ns) = flxsec(kflux(ifirst),iflux(1,ifirst),ntotal,az
     +				,fluxpm(1,ns),fluxpm(2,ns))
	  !write(99,*) ns,ifirst,ilast,nnode,ntotal,flux(ns)
	end do

	end

c******************************************************************

	function flxsec(kflux,iflux,n,az,fluxp,fluxm)

c computes flux through one section

	implicit none

	real flxsec
	integer n
	integer kflux(1)
	integer iflux(3,1)
	real az
	real fluxp,fluxm

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer i,k
	integer istype,iafter,ibefor
	real tt,tp,tm
	real port

	real flxnov

	tt = 0.
	tp = 0.
	tm = 0.

	do i=1,n
		k = kflux(i)
		istype = iflux(1,i)
		ibefor = iflux(2,i)
		iafter = iflux(3,i)
		port = flxnov(k,ibefor,iafter,istype,az)
		tt = tt + port
		!write(99,*) ' ... ',i,k,port,tt	!ggu99
		if( port .gt. 0. ) then
		  tp = tp + port
		else
		  tm = tm - port
		end if
		!call write_node_fluxes(k,port)		!special output
	end do

	fluxp = tp
	fluxm = tm
	flxsec = tt

	end
	  
c******************************************************************

	function flxnov(k,ibefor,iafter,istype,az)

c computes flux through finite volume k
c
c internal section is defined by:  kbefor - k - kafter

	implicit none

	real flxnov
	integer k,ibefor,iafter,istype
	real az

	integer ndim		!must be at least ngr
	parameter (ndim=100)

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	include 'ev.h'
	include 'links.h'
        real uov(1),vov(1),unv(1),vnv(1)
        common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov

c	integer nnode,ifirst,ilast
c	integer ntotal
	integer i,ip,ie,ii,n,ne
	integer ipf,ipl
	real aj,area,dz,uv,rdt,dt
	real b,c
	real azt,tt
c	logical bstop

	real transp(ndim)
	real weight(ndim)
	real weight1(ndim)

	integer ithis

	call get_timestep(dt)
	rdt = 1./dt
	azt = 1. - az

c get pointer into link structure

	call set_elem_links(k,ne)

	n = ne
	if( istype .gt. 1 ) n = n + 1		!boundary
	if( n .gt. ndim ) stop 'error stop flxnod: ndim'

        transp(n) = 0.                !BUG FIX 29.5.2004

c compute transports into finite volume of node k -> transp
c computed transports are divergence corrected

	do i=1,ne
	  ie = lnk_elems(i)
	  ii = ithis(k,ie)
	  aj = ev(10,ie)
	  area = 4. * aj
	  dz = zenv(ii,ie) - zeov(ii,ie)
	  b = ev(3+ii,ie)
	  c = ev(6+ii,ie)
	  uv = az * ( unv(ie) * b + vnv(ie) * c )
	  uv = uv + azt * ( uov(ie) * b + vov(ie) * c )
	  uv = 12. * aj * uv
	  uv = uv - dz * area * rdt
	  transp(i) = uv
	end do

c compute transport through section in finite volume k
c
c flux through section kbefor-k must is negative in sign

	tt = 0.

	call mkweig(n,istype,ibefor,weight)

	do i=1,n
	  tt = tt + weight(i) * transp(i)
	end do

	tt = - tt

	call mkweig(n,istype,iafter,weight1)

	do i=1,n
	  tt = tt + weight1(i) * transp(i)
	end do

	flxnov = tt

	if( abs(tt) .lt. -0.1 ) then
		write(99,*) 'errorrrrrrr flxnov',abs(tt)	!ggu99
		write(99,*) n,istype,ibefor,iafter
		write(99,*) aj,area,ie,dz,uv,i
                write(99,*) (weight(i),weight1(i),transp(i),i=1,n)
	end if

	end
	

c******************************************************************

	subroutine mkweig(n,istype,is,weight)

c computes weight over internal section in one finite volume
c
c n		dimension (grade of node)
c istype	type of node (inner, border, ...)
c is		number of internal section to compute
c weigth	weights (return)

	implicit none

	integer n,istype,is
	real weight(1)

	integer it,i
	real start,fact,dw

	logical debug
	save debug
	data debug /.false./

	it = is

	do i=1,n
	  weight(i) = 0.
	end do

	if( is .eq. 0 ) then		!no section given
c	  nothing
	else if( istype .eq. 1 ) then
	  start = -0.5 * (n-1)
	  fact = 1./n
	  do i=1,n
	    weight(it) = start * fact
	    it = it + 1
	    if( it .gt. n ) it = 1
	    start = start + 1.
	  end do
	else if( istype .eq. 2 ) then
	  if( it .eq. 1 ) it = n
	  do i=it,n-1
	    weight(i) = -1.
	  end do
	else if( istype .eq. 3 ) then
	  do i=it,n-1
	    weight(i) = -1.
	  end do
	else if( istype .eq. 4 ) then
	  do i=1,it-1
	    weight(i) = +1.
	  end do
	else if( istype .eq. 5 ) then
	  if( n .eq. 2 ) then	!handle exception
	    weight(1) = -0.5 + (it-1)		! -0.5/+0.5
	  else
	    fact = 1./(n-1)
	    dw = fact * ( 1. + 1./(n-2) )
	    start = 1.
	    do i=1,n-1
	      weight(i) = (i-1) * dw
	      if( i .ge. it ) then	!upper right triangle
	        weight(i) = weight(i) - start
	      end if
	    end do
	  end if
	else
	  write(6,*) n,istype,is
	  stop 'error stop mkweig : internal error (1)'
	end if
	  
	if( debug .and. istype .ge. 3 ) then
	  write(6,*) n,istype,is
	  write(6,*) (weight(i),i=1,n)
	end if

	end

c******************************************************************

	subroutine flxini

c initializes flux routines finally

	implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	integer kantv(2,1)
	common /kantv/kantv
	
	integer idummy

	integer klineck

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors in section $FLUX'
	  stop 'error stop : flxini'
	end if

c now set info structure for sections

	call flxin(kflux,iflux,kfluxm)

	end

c******************************************************************

	subroutine flxin(kflux,iflux,n)

c sets up info structure iflux(3,1) from kflux(1)

	implicit none

	integer n		!size of kflux, iflux
	integer kflux(1)
	integer iflux(3,1)

	integer nnode,ifirst,ilast,ntotal
	logical nextline

	nnode = 0

	do while( nextline(kflux,n,nnode,ifirst,ilast) )
	  ntotal = ilast - ifirst + 1
c	  write(6,*) n,nnode,ifirst,ilast,ntotal
	  call flxinf(kflux(ifirst),iflux(1,ifirst),ntotal)
	end do

	end

c******************************************************************

	function flxtype(k)

c determines type of node k (1-5)
c
c uses inodv only for determination of open bounday nodes
c else uses kantv

	implicit none

	integer flxtype
	integer k

	integer inodv(1)
	common /inodv/inodv
	integer kantv(2,1)
	common /kantv/kantv

	integer ktype
	integer kafter,kbefor

	integer ipext

	include 'testbndo.h'

	if( kantv(1,k) .eq. 0 ) then			!inner point
	   ktype = 1
	else if( .not. is_external_boundary(k) ) then	!material boundary
	   ktype = 2
	else if( is_external_boundary(k) ) then		!open boundary
	   kafter = kantv(1,k)
	   kbefor = kantv(2,k)
	   if( is_external_boundary(kbefor) 
     +			.and. is_external_boundary(kafter) ) then	!OOO
		ktype = 5
	   else if( is_external_boundary(kbefor) ) then			!OOB
		ktype = 4
	   else if( is_external_boundary(kafter) ) then			!BOO
		ktype = 3
	   else
		write(6,*) 'error at open boundary node ',ipext(k)
		write(6,*) 'bounadry consisting of one node'
		stop 'error stop flxtype'
	   end if
	else
	   stop 'error stop flxtype: internal error (1)'
	end if

	flxtype = ktype

	end

c******************************************************************

	subroutine flxinf(kflux,iflux,n)

c sets up info structure iflux(3,1) for one section

	implicit none

	integer n
	integer kflux(1)
	integer iflux(3,1)

	integer i,k
	integer ktype
	integer kafter,kbefor

	integer igtnsc,flxtype

	do i=1,n
	  k = kflux(i)
	  ktype = flxtype(k)

	  iflux(1,i) = ktype

	  kbefor = 0
	  if( i .ne. 1 ) kbefor = kflux(i-1)
	  kafter = 0
	  if( i .ne. n ) kafter = kflux(i+1)

	  iflux(2,i) = igtnsc(k,kbefor)
	  iflux(3,i) = igtnsc(k,kafter)
	end do

	end

c******************************************************************

	function igtnsc(k1,k2)

c gets number of internal section in link index of k1

	implicit none

	integer igtnsc
	integer k1,k2

	include 'links.h'

	integer k,i,n,ie

	integer knext,kbhnd

	call set_elem_links(k1,n)

c	deal with no section given

	igtnsc = 0
	if( k2 .le. 0 ) return

c	section in the middle

	do i=1,n
	  igtnsc = i
	  ie = lnk_elems(i)
	  k = knext(k1,ie)
	  if( k .eq. k2 ) return
	end do

c	in case we are on the boundary

	i = n
	igtnsc = igtnsc + 1
	ie = lnk_elems(i)
	k = kbhnd(k1,ie)
	if( k .eq. k2 ) return

c	no node found

	write(6,*) k1,k2
	write(6,*) k1,n
	write(6,*) (lnk_elems(i),i=1,n)
	call set_elem_links(k2,n)
	write(6,*) k2,n
	write(6,*) (lnk_elems(i),i=1,n)
	stop 'error stop: internal error igtnsc (2)'
	end
	      
c******************************************************************

	subroutine flxe2i(kflux,n,berror)

c converts external to internal node numbers

	implicit none

	integer n
	integer kflux(1)
	logical berror

	integer k,kk
	integer ipint

	berror = .false.

        do k=1,n

           kk=kflux(k)

	   if( kk .lt. 0 ) then
               write(6,*) 'negative node number ',kk
               berror=.true.
	   else if( kk .gt. 0 ) then
               kk=ipint(kk)
               if(kk.le.0) then
                 write(6,*) 'node not found ',kflux(k)
                 berror=.true.
               end if
           end if

           kflux(k)=kk

        end do

	end

c******************************************************************

	function kflxck(kflux,n)

c checks compatibility of kflux -> returns number of sections or -1 (error)

	implicit none

	integer kflxck
	integer n		!size of kflux
	integer kflux(1)

	integer nnode,ifirst,ilast,ntotal
	integer nsect
	logical bstop
	logical nextline

	bstop = .false.
	nnode = 0
	nsect = 0

	do while( nextline(kflux,n,nnode,ifirst,ilast) )
	  nsect = nsect + 1
	  ntotal = ilast - ifirst + 1
	  call flxdbl(kflux(ifirst),ntotal,bstop) !tests if nodes unique
	  call flxadj(kflux(ifirst),ntotal,bstop) !tests if nodes are adjacent
	end do

	if( bstop ) nsect = -1

	kflxck = nsect

	end

c******************************************************************

	subroutine flxdbl(kflux,n,bstop)

c tests for uniqueness

	implicit none

	integer n
	integer kflux(1)
	logical bstop

	integer i,k,j
	integer istart

	integer ipext

        istart = 1
        if( kflux(1) .eq. kflux(n) ) istart = 2         !allow for closed line

	do i=istart,n
	   k = kflux(i)
	   do j=i+1,n
	      if( kflux(j) .eq. k ) then
		bstop = .true.
		write(6,*) 'node is not unique ',ipext(k)
	      end if
	   end do
	end do

	end

c******************************************************************

	subroutine flxadj(kflux,n,bstop)

c tests for adjacency

	implicit none

	integer n
	integer kflux(1)
	logical bstop

	integer i,k1,k2

	integer ipext
	logical iskadj

	do i=2,n
	   k1 = kflux(i-1)
	   k2 = kflux(i)

	   if( .not. iskadj(k1,k2) ) then
		bstop = .true.
		write(6,*) 'nodes are not adjacent ',ipext(k1),ipext(k2)
	   end if
	end do

	end

c******************************************************************

	function extrsc(inode,ndim,nnode,ifirst,ilast)

c extracts section from array
c
c on entry:
c           nnode points to first node to be analysed
c on return inode(ifirst) is first node of section and
c           inode(ilast)  is last  node of section
c           nnode points to first node to be analysed on entry

	implicit none

	logical extrsc		!true if new section found
	integer inode(1)	!list of nodes that define sections    (in)
	integer ndim		!total number of nodes given in node() (in)
	integer nnode		!start/end of section in inode()       (in/out)
	integer ifirst		!start of section in inode()           (out)
	integer ilast		!end of section in inode()             (out)

	integer iw,i,node

	iw = -1				! -1 befor  0 in   +1 after section
	if( nnode .le. 0 ) nnode = 1
	i = nnode

	do while( i .le. ndim .and. iw .ne. 1 )
	    node = inode(i)
	    if( node .lt. 0 ) then
	        stop 'error stop extrsc: negative nodes in section'
	    else if( node .eq. 0 ) then
	        if( iw .eq. 0 ) then	!end of section
		    ilast = i - 1
		    iw = +1
	        end if
	    else !if( node .gt. 0 ) then
	        if( iw .eq. -1 ) then	!start of section
		    ifirst = i
		    iw = 0
	        end if
	    end if
	    i = i + 1
	end do

c end last section if not already done

	if( iw .eq. 0 ) ilast = ndim		!last section

	if( iw .eq. -1 ) then		!no more section, only trailing blanks
		nnode = 0
		extrsc = .false.
	else !if( iw .eq. 1 ) then	!section ended regularily
		nnode = i
		extrsc = .true.
	end if

	end

c**********************************************************************

	subroutine extrsect(isect,knode,nnode)

c extracts nodes of one section into array

	implicit none

	integer isect		!section to extract			(in)
	integer knode(1)	!list of nodes extracted		(out)
	integer nnode		!total number of nodes extracted	(out)

c nnode is 0 if section has not been found
c if isect < 0 the order of the nodes in knode is reversed
c
c knode must be big enough to contain all nodes of section isect

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux

	integer nsabs,ifirst,ilast
	integer i
	logical extrline

	nnode = 0
	nsabs = abs(isect)

	if( isect .eq. 0 ) return		!no such section
	if( nsabs .gt. nsect ) return		!no such section

	if( .not. extrline(kflux,kfluxm,nsabs,ifirst,ilast) ) return

	do i=ifirst,ilast
	  nnode = nnode + 1
	  knode(nnode) = kflux(i)
	end do

	if( isect .lt. 0 ) call revline(nnode,knode)

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine make_fluxes(k,n,itype,rflux,tflux)

c computes fluxes over sides (tflux) from fluxes into node (rflux)

	implicit none

	integer k		!node
	integer n		!number of sides (tfluxes)
	integer itype		!type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
	integer rflux(n)	!fluxes into node (element)
	integer tflux(n)	!fluxes through sides (return value)

	integer i
	real rr

	if( itype .eq. 1 ) then		!internal node
		rr = 0.
		do i=1,n-1
		  rr = rr + i * rflux(i)
		end do
		rr = rr / n

		tflux(n) = rr
		do i=n-1,1,-1
		  tflux(i) = tflux(i+1) - rflux(i)
		end do
	else if( itype .eq. 2 ) then	!node on material boundary
		tflux(1) = 0.
		do i=2,n-1
		  tflux(i) = tflux(i-1) + rflux(i-1)
		end do
		tflux(n) = 0.
	else if( itype .eq. 3 ) then	!BOO - boundary on left
		tflux(n) = 0.
		do i=n-1,1,-1
		  tflux(i) = tflux(i+1) - rflux(i)
		end do
	else if( itype .eq. 4 ) then	!OOB - boundary on right
		tflux(1) = 0.
		do i=2,n
		  tflux(i) = tflux(i-1) + rflux(i-1)
		end do
	else if( itype .eq. 5 ) then	!node totaly on open boundary
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

	subroutine write_node_fluxes(k,port)

c special routine to write single fluxes (and water level) at nodes
c
c please uncomment the call to this subroutine if needed

	implicit none

	integer k	!node number
	real port	!discharge through finite volume k

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real znv(1)
	common /znv/znv

	real z

	z = znv(k)

	write(155,*) it,k,z,port

	end

c**********************************************************************

