
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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
c 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3 (volnod) (connod)
c 19.01.1994	ggu	$$conz - implementation of concentration (vol/connod)
c 20.01.1994	ggu	$$lumpc - evaluate conz at node (vol/connod)
c 24.03.1994	ggu	from volnod (surel,conele)
c 07.04.1995	ggu	copied from volno3 (volz3)
c 06.08.1997	ggu	use zenv for water level (volz3)
c 04.12.1997	ggu	concentration not adjusted anymore (volnod,surel)
c 04.12.1997	ggu	concentration adjusted in own routine (connod,conele)
c 30.04.1998	ggu	routines from subflx.f merged into this file
c 20.06.1998	ggu	two dimensional common blocks regolarized (nen3v,...)
c 21.08.1998	ggu	xv eliminated
c 22.10.1999	ggu	file cleaned, added 3D routines
c 16.01.2001	ggu	new concentration unique for node (connod)
c 20.11.2001	ggu	input concentration with ambient value ($AMB0)
c 25.03.2003	ggu	new routine zrise
c 28.04.2009	ggu	links re-structured 
c 23.03.2010	ggu	changed v6.1.1
c 18.06.2014	ggu	changed VERS_6_1_77
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.12.2015	ggu	changed VERS_7_3_16
c 28.04.2016	ggu	changed VERS_7_5_9
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c 
c*****************************************************************

	subroutine volnod(k,vol,dz)

c inputs volume vol into finite volume (node) k by changing water level z
c
c k	node where to input
c vol	volume to input
c dz	achieved water level change (return)

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real vol,dz
c local
	real area
	integer ie,i,ii,nl
	integer ibase
	integer elems(maxlnk)

	integer ithis

	call get_elems_around(k,maxlnk,nl,elems)

	area=0.
	do i=1,nl
	  ie=elems(i)
	  if(ie.le.0) stop 'error stop volnod: internal error'
	  area=area+ev(10,ie)
	end do
	area=area*4.

	dz=vol/area

	do i=1,nl
	  ie=elems(i)
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

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer ielem
	real dz
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

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real con,dz
	real coe(3,1)
c local
	integer ie,i,ii,nl
	integer ibase
	integer elems(maxlnk)
	real depth
	real area,vol,dvol,voltot,cnew
	real massold,massnew
	real dvoltot

	integer ithis

	call get_elems_around(k,maxlnk,nl,elems)

	voltot = 0.
	dvoltot = 0.
	cnew = 0.
	massold = 0.
	massnew = 0.

	do i=1,nl
	  ie=elems(i)
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
	  ie=elems(i)
	  ii = ithis(k,ie)
	  coe(ii,ie) = cnew / voltot
	end do

	massnew = cnew
c	write(6,*) 'ggu0: ',dvoltot
c	write(6,*) 'ggu1: ',massold,massnew,massnew-massold

	massnew = 0.
	do i=1,nl
	  ie=elems(i)
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

	use mod_hydro
	use basin

	implicit none

c arguments
	integer ielem
	real dz,con
	real coe(3,1)
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

	use mod_hydro
	use evgeom
	use basin

	implicit none

c arguments
	integer k
	real dvol
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

