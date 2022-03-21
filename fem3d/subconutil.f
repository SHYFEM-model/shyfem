
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001-2012,2014-2019  Georg Umgiesser
!    Copyright (C) 2008  Christian Ferrarin
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

c routines for concentration (utilities) (old subcon1.f)
c
c contents :
c
c subroutine conini(nlvddi,c,cref,cstrat)		sets initial conditions
c
c subroutine conmima(nlvddi,c,cmin,cmax)                 computes min/max
c subroutine conmimas(nlvddi,c,cmin,cmax)                computes scalar min/max
c
c revision log :
c
c 19.08.1998	ggu	call to conzfi changed
c 20.08.1998	ggu	makew removed (routine used is sp256w)
c 24.08.1998	ggu	levdbg used for debug
c 26.08.1998	ggu	subroutine convol, tstvol transferred to newchk
c 26.08.1998	ggu	all subroutines re-written more generally
c 26.01.1999	ggu	can be used also with 2D routines
c 16.11.2001	ggu	subroutine conmima and diffstab
c 05.12.2001	ggu	new routines diffstab,diffstab1,difflimit
c 11.10.2002	ggu	commented diffset
c 09.09.2003	ggu	new routine con3bnd
c 13.03.2004	ggu	new routines set_c_bound, distribute_vertically
c 13.03.2004	ggu	exec routine con3bnd() only for level BC (LEVELBC)
c 14.03.2004	ggu	new routines open_b_flux
c 05.01.2005	ggu	routine to write 2d nos file into subnosa.f
c 07.01.2005	ggu	routine diffwrite deleted
c 14.01.2005	ggu	new file for diffusion routines (copied to subdif.f)
c 23.03.2006	ggu	changed time step to real
c 31.05.2007	ggu	reset BC of flux type to old way (DEBHELP)
c 07.04.2008	ggu	deleted set_c_bound
c 08.04.2008	ggu	cleaned, deleted distribute_vertically, open_b_flux
c 09.10.2008	ggu&ccf	call to confop changed -> nlv
c 23.03.2010	ggu	changed v6.1.1
c 14.07.2011	ggu	changed VERS_6_1_27
c 01.06.2012	ggu	changed VERS_6_1_53
c 20.01.2014	ggu	new writing format for nos files in confop, confil
c 28.01.2014	ggu	changed VERS_6_1_71
c 26.11.2014	ggu	changed VERS_7_0_7
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 26.02.2015	ggu	changed VERS_7_1_5
c 05.05.2015	ggu	changed VERS_7_1_10
c 05.06.2015	ggu	changed VERS_7_1_12
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 29.09.2015	ggu	changed VERS_7_2_5
c 18.12.2015	ggu	changed VERS_7_3_17
c 28.04.2016	ggu	changed VERS_7_5_9
c 07.06.2016	ggu	changed VERS_7_5_12
c 03.11.2017	ggu	new routines to write shy files scalar_output_*()
c 22.02.2018	ggu	changed VERS_7_5_42
c 03.04.2018	ggu	changed VERS_7_5_43
c 06.07.2018	ggu	changed VERS_7_5_48
c 16.02.2019	ggu	changed VERS_7_5_60
c 20.03.2022	ggu	started discommissioning file
c 21.03.2022	ggu	only some subroutines save to this new file
c
c*****************************************************************

	subroutine conini(nlvddi,c,cref,cstrat,hdko)

c sets initial conditions (with stratification)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)	!variable to initialize
	real cref		!reference value
	real cstrat		!stratification [conc/km]
	real hdko(nlvddi,nkn)	!layer thickness
c local
	integer k,l
	real depth,hlayer

	do k=1,nkn
	  depth=0.
	  do l=1,nlvddi
	    hlayer = 0.5 * hdko(l,k)
	    depth = depth + hlayer
	    c(l,k) = cref + cstrat*depth/1000.
	    depth = depth + hlayer
	  end do
	end do

	end

c*************************************************************
c*************************************************************
c*************************************************************


c*************************************************************
c*************************************************************
c*************************************************************

        subroutine conmima(nlvddi,c,cmin,cmax)

c computes min/max for scalar field

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c local
	integer k,l,lmax
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
	  end do
	end do

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

        subroutine conmimas(nlvddi,c,cmin,cmax)

c computes min/max for scalar field -> writes some info

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c common
	include 'femtime.h'
c local
	integer k,l,lmax
        integer ntot
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

        ntot = 0
	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
            if( cc .le. 0. ) then
                    ntot = ntot + 1
                    write(96,*) it,l,k,cc,ntot
            end if
	  end do
	end do

        if( ntot .gt. 0 ) then
                write(96,*) 'ntot: ',it,ntot
        end if

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

