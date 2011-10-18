c
c $Id: subflxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
c subroutines for computing discharge / flux
c
c contents :
c
c subroutine inflxa
c subroutine rdflxa
c subroutine ckflxa
c subroutine prflxa
c subroutine tsflxa
c
c subroutine wrflxa(it)				write of flux data
c
c subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
c subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
c
c subroutine flxini				initializes flux routines
c subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
c subroutine flxinf(m,kflux,iflux)		sets up one info structure
c function igtnsc(k1,k2)			gets number of internal section
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
c 21.09.2011    ggu     some lower-level subroutines copied to subflx.f
c 07.10.2011    ggu     adjusted for 3d flux routines
c
c notes :
c
c These routines can also be used internally to compute the flux
c over various sections. The following calling sequence must be respected:
c
c call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
c
c call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
c
c Initialization can be done anytime.
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
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
c but since at this stage we do not have all the arrays set up
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
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine wrflxa(it)

c administers writing of flux data

	implicit none

	include 'param.h'

	integer it

        integer iscdim
        parameter(iscdim=500)

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	save /kfluxc/,/iflux/	!ggu

	integer itend
	integer j,i,l,lmax,nlmax
	real az,azpar,rr

	integer iround,ideffi
	real getpar

	real fluxes(0:nlvdim,3,iscdim)

	integer nlayers(iscdim)	!number of layers in section
	real fluxest(0:nlvdim,3,iscdim)	!accumulator - may be double precision

        integer idtflx,itflx,itmflx,nr
        integer icall,nbflx,nvers,idfile
        save fluxest
        save idtflx,itflx,itmflx,nr
        save icall,nbflx,nvers,idfile
	save nlayers
        data icall,nbflx,nvers,idfile /0,0,3,537/

c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------

        if( icall .eq. -1 ) return

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

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

		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

                nr = 0
                do i=1,nsect
	  	  lmax = nlayers(i)
	  	  do l=0,lmax
          	    fluxest(l,1,i) = 0.
          	    fluxest(l,2,i) = 0.
          	    fluxest(l,3,i) = 0.
		  end do
                end do

                nbflx=ideffi('datdir','runnam','.flx','unform','new')
                if(nbflx.le.0) then
        	   stop 'error stop wrflxa : Cannot open FLX file'
		end if

                write(nbflx) idfile,nvers
                write(nbflx) nsect,kfluxm,idtflx
                write(nbflx) (kflux(i),i=1,kfluxm)

		if( nvers .ge. 3 ) then
                  write(nbflx) nlmax,(nlayers(i),i=1,nsect)	!nvers=3
		end if

c               here we could also compute and write section in m**2

        end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

        icall = icall + 1

        if( it .lt. itmflx ) return

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------

        nr = nr + 1
	call getaz(azpar)
	az = azpar

	call flxscs(kfluxm,kflux,iflux,az,fluxes)

	do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
	    fluxest(l,1,i) = fluxest(l,1,i) + fluxes(l,1,i)
	    fluxest(l,2,i) = fluxest(l,2,i) + fluxes(l,2,i)
	    fluxest(l,3,i) = fluxest(l,3,i) + fluxes(l,3,i)
	  end do
	end do

        if( it .lt. itflx ) return

c	-------------------------------------------------------
c	write results
c	-------------------------------------------------------

        itflx=itflx+idtflx

        rr=1./nr

        do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
            fluxes(l,1,i) = fluxest(l,1,i) * rr
            fluxes(l,2,i) = fluxest(l,2,i) * rr
            fluxes(l,3,i) = fluxest(l,3,i) * rr
	  end do
        end do

	if( nvers .ge. 3 ) then
          write(nbflx) it,nsect
     +			,(nlayers(i)
     +			,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +			,i=1,nsect)
	else
          write(nbflx) it,nsect
     +			,(fluxes(0,1,i),i=1,nsect)
     +                  ,(fluxes(0,2,i),fluxes(0,3,i),i=1,nsect)
	end if

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

        nr = 0
        do i=1,nsect
	  lmax = nlayers(i)
	  do l=0,lmax
            fluxest(l,1,i) = 0.
            fluxest(l,2,i) = 0.
            fluxest(l,3,i) = 0.
          end do
        end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flxscs(kfluxm,kflux,iflux,az,fluxes)

c computes flux through all sections and returns them in flux
c
c flux are divided into total, positive and negative

	implicit none

	include 'param.h'

	integer kfluxm
	integer kflux(1)
	integer iflux(3,1)
	real az
	real fluxes(0:nlvdim,3,1)

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline

	nnode = 0
	ns = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  call flxsec(ntotal,kflux(ifirst),iflux(1,ifirst),az
     +				,fluxes(0,1,ns))
	end do

	end

c******************************************************************

	subroutine flxsec(n,kflux,iflux,az,fluxes)

c computes flux through one section

	implicit none

	include 'param.h'

	integer n
	integer kflux(1)
	integer iflux(3,1)
	real az
	real fluxes(0:nlvdim,3)

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer i,k,l,lkmax
	integer istype,iafter,ibefor
	real ftot,fpos,fneg
	real port,ptot,port2d
	real flux(nlvdim)

	real flxnov

	ftot = 0.
	fpos = 0.
	fneg = 0.

	do i=1,n
		k = kflux(i)
		istype = iflux(1,i)
		ibefor = iflux(2,i)
		iafter = iflux(3,i)

		port = flxnov(k,ibefor,iafter,istype,az)
		port2d = port
		fluxes(0,1) = fluxes(0,1) + port
		if( port .gt. 0. ) then
		  fluxes(0,2) = fluxes(0,2) + port
		else
		  fluxes(0,3) = fluxes(0,3) - port
		end if

		ptot = 0.
		call flx3d(k,ibefor,iafter,istype,az,lkmax,flux)
		do l=1,lkmax
		  port = flux(l)
		  ptot = ptot + port
		  fluxes(l,1) = fluxes(l,1) + port
		  if( port .gt. 0. ) then
		    fluxes(l,2) = fluxes(l,2) + port
		  else
		    fluxes(l,3) = fluxes(l,3) - port
		  end if
		end do

		port = abs(port2d-ptot)
		if( port .gt. 1. ) then
		  write(6,*) '***** integrated fluxes: ',k,port
		  write(6,*) '   ',port2d,ptot
		end if
	end do

	end
	  
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flxini

c initializes flux routines finally (wrapper for flx_init)

	implicit none

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
	integer iflux(3,1)
	common /iflux/iflux
	
	call flx_init(kfluxm,kflux,nsect,iflux)

	end

c******************************************************************

	subroutine flx_init(kfluxm,kflux,nsect,iflux)

c sets up array iflux

	implicit none

        integer kfluxm		!total number of nodes in kflux
        integer kflux(1)	!nodes in sections
	integer nsect		!number of section (return)
	integer iflux(3,1)	!internal array for flux computation (return)

	integer ifirst,ilast,nnode,ntotal

	integer klineck
	logical nextline

c----------------------------------------------------------
c check nodes for compatibility
c----------------------------------------------------------

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors in section $FLUX'
	  stop 'error stop : flxini'
	end if

c----------------------------------------------------------
c now set info structure for sections
c----------------------------------------------------------

	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ntotal = ilast - ifirst + 1
c	  write(6,*) kfluxm,nnode,ifirst,ilast,ntotal
	  call flxinf(ntotal,kflux(ifirst),iflux(1,ifirst))
	end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c******************************************************************

	subroutine flxinf(n,kflux,iflux)

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
	      
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine get_nlayers(kfluxm,kflux,nlayers,nlmax)

	implicit none

	integer kfluxm
	integer kflux(1)
	integer nlayers(1)
	integer nlmax

	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer ns
	integer nnode,ifirst,ilast
	integer i,k,l,lmax

	logical nextline

	ns = 0
	nnode = 0
	nlmax = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  lmax = 0
	  do i=ifirst,ilast
	    k = kflux(i)
	    l = ilhkv(k)
	    lmax = max(lmax,l)
	  end do
	  nlayers(ns) = lmax
	  nlmax = max(nlmax,lmax)
	end do

	end

c**********************************************************************

