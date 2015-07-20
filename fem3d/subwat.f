c
c $Id: subwat.f,v 1.15 2009-05-21 09:24:00 georg Exp $
c
c water volume, surface and flux routines
c
c contents :
c
c subroutine volnod(k,vol,dz)	inputs volume vol into finite volume (node) k
c subroutine zrise(dz)		inputs water distributed over total surface
c subroutine surel(ielem,dz)	inputs water distributed over surface
c subroutine connod(k,dz,con,coe)	changes concentration in node k
c subroutine conele(ielem,dz,con,coe)	changes concentration in element ielem
c
c subroutine volz3(k,dvol)                      inputs water volume
c
c subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)
c       computes volume and total mass of concentration in column of node k
c subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)
c       inputs concentration in finite volume (node) k
c
c revision log :
c
c 30.04.1998    ggu     routines from subflx.f merged into this file
c 20.06.1998    ggu     two dimensional common blocks regolarized (nen3v,...)
c 21.08.1998    ggu     xv eliminated
c 22.10.1999    ggu     file cleaned, added 3D routines
c 20.11.2001    ggu     input concentration with ambient value ($AMB0)
c 25.03.2003    ggu     new routine zrise
c 28.04.2009    ggu     links re-structured 
c 
c*****************************************************************

	subroutine volnod(k,vol,dz)

c inputs volume vol into finite volume (node) k by changing water level z
c
c k	node where to input
c vol	volume to input
c dz	achieved water level change (return)
c
c written 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c revised 19.01.94 by ggu   $$conz - implementation of concentration
c revised 20.01.94 by ggu   $$lumpc - evaluate conz at node
c revised 04.12.97 by ggu   concentration not adjusted anymore

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real vol,dz
c common
	include 'param.h'
c local
	real area
	integer ie,i,ii,nl
	integer ibase

	integer ithis

	call set_elem_links(k,nl)

	area=0.
	do i=1,nl
	  ie=lnk_elems(i)
	  if(ie.le.0) stop 'error stop volnod: internal error'
	  area=area+ev(10,ie)
	end do
	area=area*4.

	dz=vol/area

	do i=1,nl
	  ie=lnk_elems(i)
	  ii = ithis(k,ie)
	  zenv(ii,ie)=zenv(ii,ie)+dz
	  zeov(ii,ie)=zeov(ii,ie)+dz		!FIXME ?????
	end do

	end

c*****************************************************************

	subroutine zrise(dz)

c inputs water distributed over total surface
c
c dz		rise of water level to achieve

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	real dz
c common
	include 'param.h'
c local
	integer ie,ii,k

	do ie=1,nel
	  do ii=1,3
	    zenv(ii,ie)=zenv(ii,ie)+dz
	    zeov(ii,ie)=zeov(ii,ie)+dz
	  end do
	end do

	do k=1,nkn
	  znv(k) = znv(k) + dz
	  zov(k) = zov(k) + dz
	end do

	end

c*****************************************************************

	subroutine surel(ielem,dz)

c inputs water distributed over surface to element ie
c
c ielem		element to input water, 0: all elements
c dz		rise of water level to achieve
c
c written 24.03.94 by ggu   from volnod
c revised 04.12.97 by ggu   concentration not adjusted anymore

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer ielem
	real dz
c common
	include 'param.h'
c local
	integer ie,ii
	integer ie1,ie2

	if(ielem.le.0) then
	  ie1=1
	  ie2=nel
	else
	  ie1=ielem
	  ie2=ielem
	end if

	do ie=ie1,ie2
	  do ii=1,3
	    zenv(ii,ie)=zenv(ii,ie)+dz
	    zeov(ii,ie)=zeov(ii,ie)+dz
	  end do
	end do

	end

c*****************************************************************

	subroutine connod(k,dz,con,coe)

c changes concentration according to a water level rise dz with conz. con
c in node k on variable coe
c
c k	node where to input
c dz	water level change 
c con	concentraion of water injected
c coe	variable to change
c
c written 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c revised 19.01.94 by ggu   $$conz - implementation of concentration
c revised 20.01.94 by ggu   $$lumpc - evaluate conz at node
c revised 04.12.97 by ggu   concentration adjusted in own routine
c 16.01.2001 by ggu   new concentration unique for node

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real con,dz
	real coe(3,1)
c common
	include 'param.h'
c local
	integer ie,i,ii,nl
	integer ibase
	real depth
	real area,vol,dvol,voltot,cnew
	real massold,massnew
	real dvoltot

	integer ithis

	call set_elem_links(k,nl)

	voltot = 0.
	dvoltot = 0.
	cnew = 0.
	massold = 0.
	massnew = 0.

	do i=1,nl
	  ie=lnk_elems(i)
	  if(ie.le.0) stop 'error stop connod: internal error'
	  ii = ithis(k,ie)
	  area = 4. * ev(10,ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)-dz		!$$lumpc
	  vol = area * depth
	  dvol = area * dz
	  dvoltot = dvoltot + dvol
	  voltot = voltot + vol + dvol
	  massold = massold + coe(ii,ie)*vol
	  cnew = cnew + coe(ii,ie)*vol + con*dvol
	end do

	do i=1,nl
	  ie=lnk_elems(i)
	  ii = ithis(k,ie)
	  coe(ii,ie) = cnew / voltot
	end do

	massnew = cnew
c	write(6,*) 'ggu0: ',dvoltot
c	write(6,*) 'ggu1: ',massold,massnew,massnew-massold

	massnew = 0.
	do i=1,nl
	  ie=lnk_elems(i)
	  ii = ithis(k,ie)
	  area = 4. * ev(10,ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)		!$$lumpc
	  vol = area * depth
	  massnew = massnew + coe(ii,ie)*vol
	end do
c	write(6,*) 'ggu2: ',massold,massnew,massnew-massold

	end

c*****************************************************************

	subroutine conele(ielem,dz,con,coe)

c changes concentration according to water level rise dz with conz. con
c in element ielem on variable coe
c
c ielem		element to input water, 0: all elements
c dz		rise of water level
c con		conzentration of injected water
c coe		variable to change
c
c written 24.03.94 by ggu   from volnod
c revised 04.12.97 by ggu   concentration adjusted in own routine

	use mod_hydro
	use basin

	implicit none

c arguments
	integer ielem
	real dz,con
	real coe(3,1)
c common
	include 'param.h'
c local
	real depth
	integer ie,ii
	integer ie1,ie2

	if(ielem.le.0) then
	  ie1=1
	  ie2=nel
	else
	  ie1=ielem
	  ie2=ielem
	end if

	do ie=ie1,ie2
	  do ii=1,3
	    depth=hm3v(ii,ie)+zeov(ii,ie)
	    coe(ii,ie)=(depth*coe(ii,ie)+dz*con)/(depth+dz)
	  end do
	end do

	end

c*****************************************************************
c
	subroutine volz3(k,dvol)
c
c inputs water volume into finite volume
c ( 3d version )
c
c k	node
c dvol	water volume change (for whole time step)
c
c written 07.04.95 by ggu   copied from volno3
c revised 06.08.97 by ggu   use zenv for water level

	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real dvol
c common
	include 'param.h'
c local
        integer ie,ii
        real area,zz

        area=0.

        do ie=1,nel
          do ii=1,3
            if(nen3v(ii,ie).eq.k) then
                area=area+ev(10,ie)
            end if
          end do
        end do

        area=4.*area
	zz = dvol/area
	znv(k) = znv(k) + zz

        do ie=1,nel
          do ii=1,3
            if(nen3v(ii,ie).eq.k) then
		zenv(ii,ie) = zenv(ii,ie) + zz
            end if
          end do
        end do

	end

c***************************************************************

	subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

c computes volume and total mass of concentration in column of node k
c + volume of upper layer
c new version that computes only up to layer lmax (lmax > 0)

	use levels

	implicit none

c arguments
	integer k		!node defining column			(in)
	integer lmax		!maximum level to which compute		(in)
	integer nlvddi		!vertical dimension			(in)
	real s(nlvddi,1)	!variable (temperature, salinity,...)	(in)
	real area		!area of column 			(out)
	real vol		!total volume of column			(out)
	real vol0		!volume of first layer			(out)
	real svol		!total mass of s in column		(out)
c common
	include 'param.h'
c local
	integer l,ilevel,nlev
	integer mode
	real volume
c functions
	real volnode,areanode

	nlev = lmax
	if( lmax .le. 0 ) nlev = nlvddi		!all levels
	ilevel = min(nlev,ilhkv(k))

	mode = +1	!new time level
	vol=0.
	svol=0.

	do l=1,ilevel
	  volume = volnode(l,k,mode)
	  vol = vol + volume
	  svol = svol + volume * s(l,k)
	end do

	vol0 = volnode(1,k,mode)
	area = areanode(1,k)

	end

c******************************************************

	subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)

c inputs concentration in finite volume (node) k
c ( 3d version ) -> only up to layer lmax

	use levels

	implicit none

c arguments
	integer k		!node defining volume
	integer lmax		!maximum level to which introduce
	integer nlvddi		!vertical dimension
	real s(nlvddi,1)	!variable (temperature, salinity,...)
	real dvol		!change in water volume (whole time step)
	real dcon		!concentration of new volume dvol
c common
	include 'param.h'
c local
	logical debug
	integer l,nlev
	real area,vol,vol0,svol
	real alpha,beta,gamma

	debug=.false.

        if( dcon .le. -990. ) return !input with ambient concentration $AMB0

	call volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

	nlev = ilhkv(k)
	if( lmax .gt. 0 .and. nlev .gt. lmax ) nlev = lmax

	alpha=dvol/(vol+dvol)
	beta=vol0/(vol0+dvol)
	gamma=dvol*svol/(vol*(vol0+dvol))
	gamma=dvol*dcon+alpha*(svol-dcon*vol)
	gamma=gamma/(vol0+dvol)

	do l=1,nlev
	  if(l.eq.1) then
	    s(l,k)=beta*s(l,k)+alpha*beta*(dcon-s(l,k))+gamma
	  else
	    s(l,k)=s(l,k)+alpha*(dcon-s(l,k))
	  end if
	  if(debug) write(6,*) k,l,s(l,k)
	end do

	end

c******************************************************

