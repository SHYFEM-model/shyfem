
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

c routines for interpolation onto regular grid
c
c contents :
c
c revision log :
c
c 18.11.1998	ggu	routine commented
c 18.11.1998	ggu	routine setgeo introduced
c 19.11.1998	ggu	routines a2char, prchar added
c 19.10.1999	ggu	routine mkmask added from subutl
c 25.11.2004	ggu	new routines femintp and elemintp for interpolation
c 14.03.2005	ggu	new routines for interpolation in element
c 11.03.2009	ggu	new helper routine getgeoflag()
c 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
c 26.01.2011	ggu&mbj	handling extrapolation in am2av()
c 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
c 31.03.2011	ggu	new routine elemmask()
c 24.11.2011	ggu	new routine find_close_elem()
c 20.06.2012	ggu	new routine get_scal_elem()
c 07.10.2012	ggu	new routine av2fm()
c 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
c 26.10.2012	ggu	bug fix: do not access not existing storage
c 30.05.2014	ggu	in av2amk() do not interpolate for flag values
c 07.07.2014	ggu	new routine intp_reg()
c 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
c 05.05.2016	ggu	file restructured (module)
c 14.05.2016	ggu	allow for extension of grid -> bregextend
c 23.06.2016	ggu	allow for eps in computing box
c 23.09.2016	ggu	allow for eps in computing box and reg intp
c 23.04.2017	ggu	new routine intp_reg_single_nodes()
c 25.05.2017	ggu	changed VERS_7_5_28
c 11.07.2017	ggu	changed VERS_7_5_30
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 18.10.2019	ggu	cleaned contents
c
c******************************************************
c******************************************************
c******************************************************

	subroutine set_dry_mask(bwater,zv,zev,href,hzoff)

c makes mask for dry and wet areas - zenv must be available
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	real zv(nkn)
	real zev(3,nel)
	real href,hzoff
c local
	integer itot,itot1
	integer ie,ii

        do ie=1,nel

          itot=0
          do ii=1,3
            if( hm3v(ii,ie)+zev(ii,ie)-href .gt. hzoff ) then
		itot=itot+1    !wet
	    end if
          end do

          itot1=0
          do ii=1,3
            if(zv(nen3v(ii,ie)).eq.zev(ii,ie)) itot1=itot1+1
          end do

          !if(itot.ne.3.or.itot1.ne.3)  bwater(ie) = .false.
          if(itot.ne.3.and.itot1.ne.3)  bwater(ie) = .false.

        end do

	end

c******************************************************

	subroutine set_level_mask(bwater,ilhv,level)

c makes mask for water points (level)
c
c bwater is elementwise mask:	true = water point

	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	logical bwater(nel)
	integer ilhv(nel)
	integer level
c local
	integer ie,nedry

	nedry = 0

	do ie=1,nel
	  if( level .gt. ilhv(ie) ) then
	    bwater(ie) = .false.
	    nedry = nedry + 1
	  end if
	end do

	end

c******************************************************

	subroutine make_dry_node_mask(bwater,bkwater)

c makes node mask from element mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	logical bkwater(nkn)

	integer ie,ii,k
	integer nndry,nedry

	nndry = 0
	nedry = 0

	bkwater = .false.

	do ie=1,nel
	  if( bwater(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkwater(k) = .true.
	    end do
	  else
	  end if
	end do

	end

c******************************************************

	subroutine make_dry_elem_mask(bwater,bkwater)

c makes elem mask from node mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(nel)
	logical bkwater(nkn)
c local
	integer ie,ii,k
	integer nedry

        bwater = .true.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( .not. bkwater(k) ) then
	      bwater(ie) = .false.
	    end if
	  end do
	end do

	end

c******************************************************

        subroutine info_dry_mask(bwater,bkwater)

        use basin

        implicit none

	logical bwater(nel)
	logical bkwater(nkn)
        
        integer nedry,nndry,newet,nnwet

        newet = count(bwater)
        nnwet = count(bkwater)
        nedry = nel - newet
        nndry = nkn - nnwet

	write(6,*) 'dry elements (dry/wet/total): ',nedry,newet,nel
	write(6,*) 'dry nodes    (dry/wet/total): ',nndry,nnwet,nkn

        end

c******************************************************
c******************************************************
c******************************************************

