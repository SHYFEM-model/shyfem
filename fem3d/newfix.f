c
c $Id: newfix.f,v 1.5 2009-03-24 17:49:19 georg Exp $
c
c velocities at open boundary
c
c contents :
c
c subroutine bclfix_ini_intern	internal routine to initialize common values
c subroutine bclnudge		nudges velocities on open boundaries
c subroutine bclfix		fixes velocities on open boundaries
c subroutine bclfix_ini		initialization of bclfix routines
c
c revision log :
c
c 15.09.2008    ggu     written from scratch
c 03.11.2008    ggu&dbf nudging implemented
c 12.11.2008    ggu     handle sigma coordinates
c 06.12.2008    ggu     read nbfix from STR
c 19.01.2009    ggu     no error stop in initializing when nbfix=0
c 23.03.2009    ggu     tramp from start of simulation

c*****************************************************************

        subroutine bclfix_ini_intern

c internal routine to initialize common values

	implicit none

	include 'param.h'
	include 'bclfix.h'

        real getpar

c-------------------------------------------------------------
c if nbfix is 0 no treatment of u/v on boundaries will be done
c-------------------------------------------------------------

	nbfix = 0		!number of boundary where vel have to be fixed

        nbfix = nint(getpar('nbfix'))

	bfix = .true.		!fix velocities at boundary
	bfix = .false.		!nudge velocities at boundary

	bsigma = .false.	!boundary data is on z levels
	bsigma = .true.		!boundary data is on sigma levels

	lbmax = 15		!how many levels to read from boundary file
	lbmax = 5		!how many levels to read from boundary file
	lbmax = 19		!how many levels to read from boundary file

	tramp = 43200.		!smooth init
	tramp = 0.		!smooth init
	tnudge = 3600.		!relaxation time for nudging

	fixfile = 'vel_bound.dat'	!file name with u/v boundary data

	if( nbfix .gt. 0 .and. lbmax .gt. nlxdim ) goto 99

	write(6,*) 'bclfix_ini_intern has been initialized: ',nbfix

	return
   99	continue
	write(6,*) 'nlxdim must be greater than lbdim'
	write(6,*) 'please change in newfix.f and bclfix.h'
	write(6,*) 'nlxdim = ',nlxdim,'   lbdim = ',lbmax
	stop 'error stop bcl_init: nlxdim'
	end

c*****************************************************************

        subroutine bclnudge

c nudges velocities on open boundaries

        implicit none 

        include 'param.h' 
        include 'bclfix.h' 

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        real fxv(nlvdim,neldim)      !new HYDRO deb
        real fyv(nlvdim,neldim)
        common /fxv/fxv
        common /fyv/fyv
        real ulnv(nlvdim,neldim),vlnv(nlvdim,neldim)
        common /ulnv/ulnv, /vlnv/vlnv
        real utlnv(nlvdim,neldim), vtlnv(nlvdim,neldim)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov

	integer ilhv(1)
	common /ilhv/ilhv

	real ubound(nlxdim,nkndim)
	real vbound(nlxdim,nkndim)
	common /ubound/ubound
	common /vbound/vbound

        integer ie,l,i,k,ii,n
	integer lmax,ifact
	integer nintp,nvar,ibc,nodes,nsize,ip
        real u(nlvdim),v(nlvdim)
        real h,t
        real alpha
	real uexpl,vexpl
        
	integer nkbnds,kbnds

	real array(nbadim)
	real rint(nbadim)
	save array

	integer icall
	save icall
	data icall /0/

	if( icall .eq. -1 ) return

c------------------------------------------------------------------
c initialize file
c------------------------------------------------------------------

	if( icall .eq. 0 ) then

	  if( fixfile .eq. ' ' ) icall = -1
	  if( nbfix .le. 0 ) icall = -1
	  if( nbc .le. 0 ) icall = -1
	  if( bfix ) icall = -1
	  if( icall .eq. -1 ) return
	  icall = 1

	  nintp = 2
	  nvar = 2
	  ibc = nbfix
	  nodes = nkbnds(ibc)

	  nsize = nodes * lbmax
	  call exffil(fixfile,nintp,nvar,nsize,nbadim,array)
	end if

	ifact = 1		!compute fact in interpolation routine

c------------------------------------------------------------------
c interpolate from file and store in u/vbound
c------------------------------------------------------------------

	t = it
	call exfintp(array,t,rint)

	ibc = nbfix
	nodes = nkbnds(ibc)
	nsize = nodes * lbmax
	ip = 0

	do i=1,nodes
	  k = kbnds(ibc,i)
	  do l=1,lbmax
	    ip = ip + 1
	    ubound(l,k) = rint(ip)
	    vbound(l,k) = rint(ip+nsize)
	  end do
	end do

c------------------------------------------------------------------
c simulate smooth initial discharge with tramp
c------------------------------------------------------------------

        if( it-itanf .le. tramp ) then
	  alpha = (it-itanf) / tramp
        else
	  alpha = 1.
        endif

c------------------------------------------------------------------
c fix velocities in elements
c------------------------------------------------------------------

	do ie=1,nel
	  n = ielfix(0,ie)
	  if( n .gt. 0 ) then

	    lmax = ilhv(ie)
	    if( .not. bsigma .and. lmax .gt. lbmax ) goto 99

	    do l=1,lmax
	      u(l) = 0.
	      v(l) = 0.
	    end do

	    do i=1,n
	      k = ielfix(i,ie)
	      if( bsigma ) then
	        call add_sigma_to_element(ie,k,ifact,nlxdim
     +                          ,ubound,vbound,u,v)
	      else
	        do l=1,lmax
	          u(l) = u(l) + ubound(l,k)
	          v(l) = v(l) + vbound(l,k)
	        end do
	      end if
	    end do

	    do l=1,lmax
	      u(l) = u(l) / n
	      v(l) = v(l) / n
	    end do

	    do l=1,lmax
              h = hdeov(l,ie)
              uexpl = (ulnv(l,ie)-u(l))*alpha*h/tnudge
              vexpl = (vlnv(l,ie)-v(l))*alpha*h/tnudge
	      fxv(l,ie) = fxv(l,ie) + uexpl
	      fyv(l,ie) = fyv(l,ie) + vexpl
	    end do

	  end if
	end do

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   99	continue
	stop 'error stop bclnudge: internal error (1)'
        end

c*****************************************************************

        subroutine bclfix

c fixes velocities on open boundaries

        implicit none 

        include 'param.h' 
        include 'bclfix.h' 

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        real ulnv(nlvdim,neldim),vlnv(nlvdim,neldim)
        common /ulnv/ulnv, /vlnv/vlnv
        real utlnv(nlvdim,neldim), vtlnv(nlvdim,neldim)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov

	integer ilhv(1)
	common /ilhv/ilhv

	real ubound(nlxdim,nkndim)
	real vbound(nlxdim,nkndim)
	common /ubound/ubound
	common /vbound/vbound

        integer ie,l,i,k,ii,n
	integer lmax,ifact
	integer nintp,nvar,ibc,nodes,nsize,ip
        real u(nlvdim),v(nlvdim)
        real h,t
        real alpha
        
	integer nkbnds,kbnds

	real array(nbadim)
	real rint(nbadim)
	save array

	integer icall
	save icall
	data icall /0/

	if( icall .eq. -1 ) return

c------------------------------------------------------------------
c initialize file
c------------------------------------------------------------------

	if( icall .eq. 0 ) then

	  if( fixfile .eq. ' ' ) icall = -1
	  if( nbfix .le. 0 ) icall = -1
	  if( .not. bfix ) icall = -1
	  if( icall .eq. -1 ) return
	  icall = 1

	  nintp = 2
	  nvar = 2
	  ibc = nbfix
	  nodes = nkbnds(ibc)

	  nsize = nodes * lbmax
	  call exffil(fixfile,nintp,nvar,nsize,nbadim,array)
	end if

	ifact = 1		!compute fact in interpolation routine

c------------------------------------------------------------------
c interpolate from file and store in u/vbound
c------------------------------------------------------------------

	t = it
	call exfintp(array,t,rint)

	ibc = nbfix
	nodes = nkbnds(ibc)
	nsize = nodes * lbmax
	ip = 0

	do i=1,nodes
	  k = kbnds(ibc,i)
	  do l=1,lbmax
	    ip = ip + 1
	    ubound(l,k) = rint(ip)
	    vbound(l,k) = rint(ip+nsize)
	  end do
	end do

c------------------------------------------------------------------
c simulate smooth initial discharge with tramp
c------------------------------------------------------------------

        if( it .le. tramp ) then
	  alpha = it / tramp
        else
	  alpha = 1.
        endif

c------------------------------------------------------------------
c fix velocities in elements
c------------------------------------------------------------------

	do ie=1,nel
	  n = ielfix(0,ie)
	  if( n .gt. 0 ) then

	    lmax = ilhv(ie)
	    if( .not. bsigma .and. lmax .gt. lbmax ) goto 99

	    do l=1,lmax
	      u(l) = 0.
	      v(l) = 0.
	    end do

	    do i=1,n
	      k = ielfix(i,ie)
	      if( bsigma ) then
	        call add_sigma_to_element(ie,k,ifact,nlxdim
     +                          ,ubound,vbound,u,v)
	      else
	        do l=1,lmax
	          u(l) = u(l) + ubound(l,k)
	          v(l) = v(l) + vbound(l,k)
	        end do
	      end if
	    end do

	    do l=1,lmax
	      u(l) = u(l) / n
	      v(l) = v(l) / n
	    end do

	    do l=1,lmax
              h = hdeov(l,ie)
              ulnv(l,ie) = u(l)*alpha
              vlnv(l,ie) = v(l)*alpha
              utlnv(l,ie) = ulnv(l,ie)*h
              vtlnv(l,ie) = vlnv(l,ie)*h
	    end do

	  end if
	end do

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   99	continue
	stop 'error stop bclfix: internal error (1)'
        end

c*****************************************************************

        subroutine bclfix_ini

c initialization of bclfix routines

        implicit none 

        include 'param.h' 
        include 'bclfix.h' 

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,neldim)
	common /nen3v/nen3v

	integer v1v(1)
	common /v1v/v1v

        integer iuvfix(neldim)
        common /iuvfix/iuvfix
        save /iuvfix/

        integer ie,l,i,ii,k,n
	integer ibc,nodes
       
	integer nkbnds,kbnds
	integer ieext
	logical bdebug

	bdebug = .false.
	bdebug = .true.

c------------------------------------------------------------------
c initialize global arrays
c------------------------------------------------------------------

        do ie=1,nel
          iuvfix(ie) = 0
        end do

	do ie=1,nel
	  do ii=0,3
	    ielfix(ii,ie) = 0
	  end do
	end do

c------------------------------------------------------------------
c nbfix is number of boundary where velocities have to be fixed
c------------------------------------------------------------------

        call bclfix_ini_intern		!initialize all fix related stuff

	if( nbfix .le. 0 .or. nbc .le. 0 ) return

	ibc = nbfix
	nodes = nkbnds(ibc)

c------------------------------------------------------------------
c flag boundary nodes
c------------------------------------------------------------------

	do k=1,nkn
	  v1v(k) = 0.
	end do

	do i=1,nodes
	  k = kbnds(ibc,i)
	  v1v(k) = 1.
	end do

c------------------------------------------------------------------
c set up iuvfix and ielfix
c------------------------------------------------------------------

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( v1v(k) .ne. 0. ) then
	      n = n + 1
	      ielfix(n,ie) = k
	    end if
	  end do
	  ielfix(0,ie) = n			!total number of nodes for ele
	  if( bfix .and. n .gt. 0 ) iuvfix(ie) = 1
	end do

c------------------------------------------------------------------
c debug output
c------------------------------------------------------------------

	if( .not. bdebug ) return

	write(6,*) 'bclfix_ini has been initialized'

	n = 0
	do ie=1,nel
	  if( iuvfix(ie) .ne. 0 ) then
	    write(6,*) 'ie = ',ie,ieext(ie),'   n = ',ielfix(0,ie)
	    n = n + 1
	  end if
	end do

	write(6,*) '-------------'
	write(6,*) 'bclfix_ini has found ',n,' elements to be fixed'
	write(6,*) 'bfix,bsigma: ',bfix,bsigma
	write(6,*) 'nbfix,lbmax: ',nbfix,lbmax
	write(6,*) 'nbadim,nlxdim: ',nbadim,nlxdim
	write(6,*) '-------------'

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

        end

c****************************************************************************

