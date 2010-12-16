c
c $Id: subele.f,v 1.15 2010-03-11 15:36:39 georg Exp $
c
c routines for dealing with elements
c
c contents :
c
c function depele(ie,mode)
c       computes total (average) depth of element ie
c subroutine dep3dele(ie,mode,nlev,h)
c       computes (average) depth of element ie for all layers
c subroutine dep3dnod(k,mode,nlev,h)
c       computes depth of node k for all layers
c function zetaele(ie,mode)
c       computes (average) zeta of element ie
c function volele(ie,mode)
c       computes volume of element ie
c function areaele(ie)
c       computes area of element ie
c subroutine nindex(ie,n,kn)
c       returns node index for one element
c subroutine setdepele(ie,h)
c       sets depth for element ie
c subroutine depvele(ie,mode,n,h)
c       computes depth of vertices in element ie
c function depfvol(ie,ii,mode)
c       computes total depth in vertex ii of element ie
c function scalele(ie,sev)
c       returns average of scalar sev in element ie
c subroutine scalvele(ie,sv,n,s)
c       returns scalar value s of vertices in element ie
c subroutine tempvele(ie,n,s)
c subroutine saltvele(ie,n,s)
c subroutine conzvele(ie,n,s)
c subroutine setsvele(ie,sev,s)
c       sets scalar values sev to s in element ie
c subroutine addsvele(ie,sev,s)
c       adds scalar value s to sev in element ie
c subroutine elenmax(text,nmax)
c       checks maximum of vertices of element
c subroutine ele2node(sev,sv)
c       computes nodal values from element values (scalar)
c subroutine elebase(ie,n,ibase)
c       returns number of vertices and pointer into element array
c subroutine setarea(levdim,area)
c       sets up area for nodes
c subroutine dvanode(l,k,mode,dep,vol,area)
c       returns depth, volume and area of node k on level l
c function depnode(l,k,mode)
c       returns depth of node k on level l
c function volnode(l,k,mode)
c       returns volume of node k on level l
c function areanode(l,k)
c       returns area of node k on level l
c subroutine setdepth(levdim,hdkn,hden,zenv,area)
c       sets up depth array for nodes
c function scalcont(mode,scal)
c       computes content of scalar in total domain
c function scalcontk(mode,scal,k)
c       computes content of scalar at node k
c
c revision log :
c
c 23.01.2004    ggu     new routines scalcontkh and scalmass
c 23.01.2004	ggu	-> FIXME: scalmass should go to other file
c 22.02.2005    ggu     routines deleted: tempvele saltvele conzvele
c 25.08.2005    ggu     bug in dvanode fixed (use * instead of +)
c 29.11.2006    ggu     in copydepth do not set to 0 old depths
c 07.11.2008    ggu     new helper routine make_new_depth()
c 16.12.2010    ggu     setdepth() changed for sigma levels
c
c****************************************************************

	function depele(ie,mode)

c computes total (average) depth of element ie

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	real depele	!depth (return)

	real hm3v(1)
	common /hm3v/hm3v
	real zenv(1)
	common /zenv/zenv
	real zeov(1)
	common /zeov/zeov

	integer ii,n,ibase
	real hmed

	call elebase(ie,n,ibase)

	hmed = 0.

	if( mode .lt. 0 ) then
	  do ii=1,n
	    hmed = hmed + hm3v(ibase+ii) + zeov(ibase+ii)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    hmed = hmed + hm3v(ibase+ii) + zenv(ibase+ii)
	  end do
	else
	  do ii=1,n
	    hmed = hmed + hm3v(ibase+ii)
	  end do
	end if

	depele = hmed / n

	end

c****************************************************************

	subroutine dep3dele(ie,mode,nlev,h)

c computes (average) depth of element ie for all layers

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	real h(1)	!depth of layers

	include 'param.h'

	integer ilhv(1)
	common /ilhv/ilhv

	real hdenv(nlvdim,neldim)
	common /hdenv/hdenv
	real hdeov(nlvdim,neldim)
	common /hdeov/hdeov

	integer l

	nlev = ilhv(ie)

	if( mode .gt. 0 ) then
	  do l=1,nlev
	    h(l) = hdenv(l,ie)
	  end do
	else if( mode .lt. 0 ) then
	  do l=1,nlev
	    h(l) = hdeov(l,ie)
	  end do
	else
	  stop 'error stop dep3dele: Cannot use mode = 0'
	end if

	end

c****************************************************************

	subroutine dep3dnod(k,mode,nlev,h)

c computes depth of node k for all layers

	implicit none

	integer k	!number of node
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	real h(1)	!depth of layers

	include 'param.h'

	integer ilhkv(1)
	common /ilhkv/ilhkv

	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov

	integer l

	nlev = ilhkv(k)

	if( mode .gt. 0 ) then
	  do l=1,nlev
	    h(l) = hdknv(l,k)
	  end do
	else if( mode .lt. 0 ) then
	  do l=1,nlev
	    h(l) = hdkov(l,k)
	  end do
	else
	  stop 'error stop dep3dnod: Cannot use mode = 0'
	end if

	end
	
c****************************************************************

	function zetaele(ie,mode)

c computes (average) zeta of element ie

	implicit none

	real zetaele	!depth (return)
	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	real zenv(1)
	common /zenv/zenv
	real zeov(1)
	common /zeov/zeov

	integer ii,n,ibase
	real hmed

	call elebase(ie,n,ibase)

	hmed = 0.

	if( mode .lt. 0 ) then
	  do ii=1,n
	    hmed = hmed + zeov(ibase+ii)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    hmed = hmed + zenv(ibase+ii)
	  end do
	else
	  hmed = 0.
	end if

	zetaele = hmed / n

	end

c****************************************************************

	function volele(ie,mode)

c computes volume of element ie

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	real volele	!volume (return)

	real area, depth
	real areaele, depele

	area = areaele(ie)
	depth = depele(ie,mode)

	volele = area * depth

	end

c****************************************************************

	function areaele(ie)

c computes area of element ie

        implicit none

        integer ie      !number of element
        real areaele    !volume (return)

	include 'ev.h'

	areaele = 12. * ev(10,ie)

	end

c****************************************************************

	subroutine nindex(ie,n,kn)

c returns node index for one element

	implicit none

	integer ie      !number of element
	integer n	!total number of nodes in element (3 for triangle)
	integer kn(1)	!node index

	integer nen3v(1)
	common /nen3v/nen3v

	integer ii,ibase

	call elebase(ie,n,ibase)

	do ii=1,n
	  kn(ii) = nen3v(ibase+ii)
	end do

	end

c****************************************************************

	subroutine setdepele(ie,h)

c sets depth for element ie

	implicit none

	integer ie      !number of element
	real h		!depth to use

	real hm3v(1)
	common /hm3v/hm3v

	integer ii,n,ibase

	call elebase(ie,n,ibase)

	do ii=1,n
	  hm3v(ibase+ii) = h
	end do

	end

c****************************************************************

	subroutine depvele(ie,mode,n,h)

c computes depth of vertices in element ie

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer n	!total number of vertices (return)
	real h(1)	!depth at vertices (return)

	real hm3v(1)
	common /hm3v/hm3v
	real zenv(1)
	common /zenv/zenv
	real zeov(1)
	common /zeov/zeov

	integer ii,ibase

	call elebase(ie,n,ibase)

	if( mode .lt. 0 ) then
	  do ii=1,n
	    h(ii) = hm3v(ibase+ii) + zeov(ibase+ii)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    h(ii) = hm3v(ibase+ii) + zenv(ibase+ii)
	  end do
	else
	  do ii=1,n
	    h(ii) = hm3v(ibase+ii)
	  end do
	end if

	end

c****************************************************************

	function depfvol(ie,ii,mode)

c computes total depth in vertex ii of element ie

	implicit none

	real depfvol
	integer ie	!number of element
	integer ii	!number of vertex
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	real hm3v(1)
	common /hm3v/hm3v
	real zenv(1)
	common /zenv/zenv
	real zeov(1)
	common /zeov/zeov

	integer ibase,n
	real h

	call elebase(ie,n,ibase)

	if( mode .lt. 0 ) then
	    h = hm3v(ibase+ii) + zeov(ibase+ii)
	else if( mode .gt. 0 ) then
	    h = hm3v(ibase+ii) + zenv(ibase+ii)
	else
	    h = hm3v(ibase+ii)
	end if

	depfvol = h

	end

c****************************************************************
c****************************************************************
c****************************************************************

	function scalele(ie,sev)

c returns average of scalar sev in element ie

	implicit none

	real scalele
	integer ie
	real sev(3,1)

	integer ii
	real sm

	sm = 0.
	do ii=1,3
	  sm = sm + sev(ii,ie)
	end do

	scalele = sm / 3.

	end

c****************************************************************

	subroutine scalvele(ie,sv,n,s)

c returns scalar value s of vertices in element ie

	implicit none

	integer ie	!number of element
	real sv(3,1)	!array of scalar value
	integer n	!total number of vertices (return)
	real s(1)	!scalar values at vertices (return)

	integer ii

	n = 3
	do ii=1,3
	  s(ii) = sv(ii,ie)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine setsvele(ie,sev,s)

c sets scalar values sev to s in element ie

	implicit none

	integer ie
	real sev(3,1)
	real s

	integer ii

	do ii=1,3
	  sev(ii,ie) = s
	end do

	end

c****************************************************************

	subroutine addsvele(ie,sev,s)

c adds scalar value s to sev in element ie

	implicit none

	integer ie
	real sev(3,1)
	real s

	integer ii

	do ii=1,3
	  sev(ii,ie) = sev(ii,ie) + s
	end do

	end

c****************************************************************

c***********************************************************

	subroutine elenmax(text,nmax)

c checks maximum of vertices of element

	implicit none

	character*(*) text
	integer nmax

	if( nmax .lt. 3 ) then
	  write(6,*) 'Dimension of vertices too small:'
	  write(6,*) text
	  stop 'error stop elenmax'
	end if

	end

c***********************************************************

	subroutine ele2node(sev,sv)

c computes nodal values from element values (scalar)

	implicit none

	real sev(1)
	real sv(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real v1v(1)
        common /v1v/v1v
	integer nen3v(1)
	common /nen3v/nen3v
	real hm3v(1)
	common /hm3v/hm3v
	real zenv(1)
	common /zenv/zenv

	integer ie,k,ii,n
	integer ip
	integer kn(10)
	integer ibase
	real h(10), s(10)
	real area,vol
	logical btype
	
	real areaele

	btype = .false.

	call elenmax('ele2node',10)

	do k=1,nkn
	  sv(k) = 0.
	  v1v(k) = 0.
	end do

	if( btype ) then

	do ie=1,nel

	  area = areaele(ie)
	  call nindex(ie,n,kn)
	  call depvele(ie,+1,n,h)
	  call scalvele(ie,sev,n,s)

	  do ii=1,n
	    k = kn(ii)
	    vol = area * h(ii)
	    sv(k) = sv(k) + s(ii) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	else

	do ie=1,nel

	  area = areaele(ie)
	  call elebase(ie,n,ibase)

	  do ip=ibase+1,ibase+n
	    k = nen3v(ip)
	    vol = area * ( hm3v(ip) + zenv(ip) )
	    sv(k) = sv(k) + sev(ip) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	end if

	do k=1,nkn
	  sv(k) = sv(k) / v1v(k)
	end do

	end

c***********************************************************

	subroutine elebase(ie,n,ibase)

c returns number of vertices and pointer into element array

	implicit none

	integer ie	!element number
	integer n	!total number of vertices (return)
	integer ibase	!pointer to start of values

	n = 3
	ibase = (ie-1) * 3

	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine setarea(levdim,area)

c sets up area for nodes

	implicit none

	integer levdim
	real area(levdim,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         
	integer ilhv(1)
	common /ilhv/ilhv
	integer nen3v(1)
	common /nen3v/nen3v

	integer k,l,ie,ii
	integer nlev,n,ibase
	real areael,areafv
	real areaele

	do k=1,nkn
	  do l=1,levdim
	    area(l,k) = 0.
	  end do
	end do

	do ie=1,nel

	  call elebase(ie,n,ibase)
	  areael = areaele(ie)
	  areafv = areael / n

	  nlev = ilhv(ie)

	  do ii=1,n

	    k = nen3v(ibase+ii)

	    do l=1,nlev
	      area(l,k) = area(l,k) + areafv
	    end do

	  end do

	end do

	end

c***********************************************************

	subroutine dvanode(l,k,mode,dep,vol,area)

c returns depth, volume and area of node k on level l

	implicit none

	integer l
	integer k
	integer mode	!-1: old    1: new 
	real dep
	real vol
	real area

	include 'param.h'

	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov
	real areakv(nlvdim,nkndim)
	common /areakv/areakv

	if( mode .gt. 0 ) then
	  dep = hdknv(l,k)
	else if( mode .lt. 0 ) then
	  dep = hdkov(l,k)
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop dvanode'
	end if

	area = areakv(l,k)
        vol = dep * area

	end

c***********************************************************

	function depnode(l,k,mode)

c returns depth of node k on level l

	implicit none

	real depnode
	integer l
	integer k
	integer mode

	include 'param.h'

	real hdknv(nlvdim,nkndim)
	common /hdknv/hdknv
	real hdkov(nlvdim,nkndim)
	common /hdkov/hdkov

	if( mode .gt. 0 ) then
	  depnode = hdknv(l,k)
	else if( mode .lt. 0 ) then
	  depnode = hdkov(l,k)
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop depnode'
	end if

	end

c***********************************************************

	function volnode(l,k,mode)

c returns volume of node k on level l

	implicit none

	real volnode
	integer l
	integer k
	integer mode

	real depnode,areanode

	volnode = depnode(l,k,mode) * areanode(l,k)

	end

c***********************************************************

	function areanode(l,k)

c returns area of node k on level l

	implicit none

	real areanode
	integer l
	integer k

	include 'param.h'

	real areakv(nlvdim,nkndim)
	common /areakv/areakv

	areanode = areakv(l,k)

	end

c***********************************************************

	subroutine copydepth(levdim,hdkn,hdko,hden,hdeo)

c sets up depth array for nodes

	implicit none

	integer levdim
	real hdkn(levdim,1)	!depth at node, new time level
	real hdko(levdim,1)	!depth at node, old time level
	real hden(levdim,1)	!depth at element, new time level
	real hdeo(levdim,1)	!depth at element, old time level

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         

	integer k,ie,l

	do k=1,nkn
	  do l=1,levdim
	    hdko(l,k) = hdkn(l,k)
	    !hdkn(l,k) = 0.
	  end do
	end do

	do ie=1,nel
	  do l=1,levdim
	    hdeo(l,ie) = hden(l,ie)
	    !hden(l,ie) = 0.
	  end do
	end do

	end

c***********************************************************

	subroutine make_new_depth

c shell (helper) for setdepth

	implicit none

	include 'param.h'

        real zenv(1)
        common /zenv/zenv
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real areakv(nlvdim,nkndim)
        common /areakv/areakv

	call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)

	end

c***********************************************************

	subroutine setdepth(levdim,hdkn,hden,zenv,area)

c sets up depth array for nodes

	implicit none

	integer levdim
	real hdkn(levdim,1)	!depth at node, new time level
	real hden(levdim,1)	!depth at element, new time level
	real zenv(1)    	!water level at new time level
	real area(levdim,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         
	integer ilhv(1)
	common /ilhv/ilhv
	integer nen3v(1)
	common /nen3v/nen3v
	real hldv(1)
	common /hldv/hldv
	real hlhv(1)
	common /hlhv/hlhv
	real hev(1)
	common /hev/hev
	real hm3v(1)
	common /hm3v/hm3v
c	real zenv(1)
c	common /zenv/zenv

        logical bdebug
        logical bsigma
	integer k,l,ie,ii
	integer nlev,n,ibase
	real hfirst,hlast,h,htot,z,zmed
	real hacu,hlevel,hlast_aux

	real areael,areafv
	real areaele

        bdebug = .false.

c----------------------------------------------------------------
c initialize and copy
c----------------------------------------------------------------

	do k=1,nkn
	  do l=1,levdim
	    hdkn(l,k) = 0.
	  end do
	end do

	do ie=1,nel
	  do l=1,levdim
	    hden(l,ie) = 0.
	  end do
	end do

c----------------------------------------------------------------
c compute volumes at node
c----------------------------------------------------------------

	hfirst = hldv(1)
	bsigma = hfirst .lt. 0.

	do ie=1,nel

	  call elebase(ie,n,ibase)
	  areael = areaele(ie)
	  areafv = areael / n

	  nlev = ilhv(ie)
	  hlast = hlhv(ie)
	  zmed = 0.

c	  -------------------------------------------------------
c	  nodal values
c	  -------------------------------------------------------

	  do ii=1,n

	    k = nen3v(ibase+ii)
	    htot = hm3v(ibase+ii)
	    z = zenv(ibase+ii)
	    zmed = zmed + z

	    if( bsigma ) then
	      h = htot + z
	      do l=1,nlev
	        hdkn(l,k) = - h * hldv(l)
	      end do
	    else if( nlev .eq. 1 ) then
	      h = htot + z
	      hdkn(1,k) = hdkn(1,k) + areafv * h
	    else
	      h = hfirst + z
	      hacu = hfirst
	      hdkn(1,k) = hdkn(1,k) + areafv * h
	      do l=2,nlev-1
	        hlevel = hldv(l)
	        hdkn(l,k) = hdkn(l,k) + areafv * hlevel
		hacu = hacu + hlevel
	      end do
	      hlast_aux = htot - hacu
	      hdkn(nlev,k) = hdkn(nlev,k) + areafv * hlast
	      if( hlast_aux .ne. hlast ) then
		write(6,*) '**** hlast node: ',ie,ii,nlev
		write(6,*) hfirst,htot,hacu,hlast,hlast_aux
	      end if
	    end if

	  end do

c	  -------------------------------------------------------
c	  element values
c	  -------------------------------------------------------

	  zmed = zmed / n
	  htot = hev(ie)

	  if( bsigma ) then
	    h = htot + zmed
	    do l=1,nlev
	      hden(l,ie) = - h * hldv(l)
	    end do
	  else if( nlev .eq. 1 ) then
	    hden(1,ie) = htot + zmed
	  else
	    hden(1,ie) = hfirst + zmed
	    hacu = hfirst
	    do l=2,nlev-1
	      hlevel = hldv(l)
	      hden(l,ie) = hlevel
	      hacu = hacu + hlevel
	    end do
	    hlast_aux = htot - hacu
	    hden(nlev,ie) = hlast
	    if( hlast_aux .ne. hlast ) then
		write(6,*) '**** hlast elem: ',ie,nlev
		write(6,*) hfirst,htot,hacu,hlast,hlast_aux
	    end if
	  end if

	end do

c----------------------------------------------------------------
c compute depth at nodes
c----------------------------------------------------------------

	if( .not. bsigma ) then
	 do k=1,nkn
	  do l=1,levdim
	    areafv = area(l,k)
	    if( areafv .gt. 0. ) then
	      hdkn(l,k) = hdkn(l,k) / areafv
	    else
	      goto 1		!no more layers
	    end if
	  end do
    1	  continue
	 end do
	end if

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c***********************************************************

	function scalcont(mode,scal)

c computes content of scalar in total domain

        implicit none

	include 'param.h'

	double precision scalcont
	integer mode
	real scal(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ilhkv(1)
	common /ilhkv/ilhkv

	logical bdebug
	integer k,l,nlev
	double precision total
        real volnode

	bdebug = .false.
	total = 0.

	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    total = total + volnode(l,k,mode) * scal(l,k)
	    if( bdebug .and. scal(l,k) .gt. 0. ) then
		write(6,*) 'scalcont: ',l,k,scal(l,k)
	    end if
	  end do
	end do

	scalcont = total

	end

c***********************************************************

	function scalcontk(mode,scal,k)
 
c computes content of scalar at node k
 
        implicit none

        include 'param.h'
 
        double precision scalcontk
        integer mode,k
        real scal(nlvdim,nkndim)
 
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
 
        integer ilhkv(1)
        common /ilhkv/ilhkv
 
        integer l,nlev
        double precision total
        real volnode
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + volnode(l,k,mode) * scal(l,k)
          end do
 
        scalcontk = total
 
        end
 
c***********************************************************

	function scalcontkh(scal,k,depth)
 
c computes content of scalar at node k (with given depth)
 
        implicit none

        include 'param.h'

        double precision scalcontkh
        integer k
        real scal(nlvdim,nkndim)
        real depth
 
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
 
        integer ilhkv(1)
        common /ilhkv/ilhkv
 
        integer l,nlev
        double precision total
        real areanode
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + depth * areanode(l,k) * scal(l,k)
          end do
 
        scalcontkh = total
 
        end
 
c***********************************************************

        subroutine scalmass(scal,depth,tstot)

c this routine into another file... FIXME

        implicit none

        include 'param.h'

        real scal(nlvdim,nkndim)        !scalar field for which to compute mass
        real depth                      !depth of column
        real tstot                      !total mass (return)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer k
        double precision scalcontkh
        double precision total

        total = 0.

        do k=1,nkn
          total = total + scalcontkh(scal,k,depth)
        end do

        tstot = total

        end

c***********************************************************

