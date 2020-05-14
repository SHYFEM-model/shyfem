
!--------------------------------------------------------------------------
!
!    Copyright (C) 2008,2010,2012,2014-2015,2018-2019  Georg Umgiesser
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

c bio3d_util - utility routines for bio3d
c
c revision log :
c
c 18.04.2008	ggu	copied from weutro_sedim.f
c 09.10.2008	ggu	new call to confop
c 23.03.2010	ggu	changed v6.1.1
c 09.03.2012	ggu	bug fix: ilhkv was real
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c
c********************************************************************

	subroutine loicz1(knode,vol,depth)

c EUTRO 0-D (LOICZ BUdgeting Procedure)
c
c new version -> does everything: initializes, accumulates, writes

	use basin

	implicit none

	integer nlzstate
	parameter(nlzstate=3)

        integer knode           !node for accumulation - 0 if init or write
        real vol                !volume [m**3]
        real depth              !depth of box [m]

        include 'donata.h'

        real, save, allocatable :: elz(:,:)     !loicz budg proc ariables

        logical bloicz
        integer i,k
	integer ierr,ivar
        real elzaux(nlzstate)    !diagnostic variable
        real tlztot(nlzstate)
	double precision dtime

	logical next_output_d
        real getpar

        integer, save :: icall = 0
	double precision, save :: da_out(4)

c-----------------------------------------------------------
c see if routine has to be executed
c-----------------------------------------------------------

        bloicz = .true.

        if( icall .eq. -1 ) return

        if( .not. bloicz ) then
          icall = -1
          return
        end if

c-----------------------------------------------------------
c initialization
c-----------------------------------------------------------

        if( icall .eq. 0 ) then

          icall = 1

	  allocate(elz(nkn,nlzstate))
	  elz = 0.

          call init_output_d('itmcon','idtcon',da_out)
          call scalar_output_init(da_out,1,nlzstate,'lcz',ierr)
          if( ierr > 0 ) goto 99
          if( ierr < 0 ) icall = -1
          if( icall < 0 ) return

          write(6,*) 'bio3d  loicz budget initialized...'

        end if 

c-----------------------------------------------------------
c accumulation of results
c-----------------------------------------------------------

        if( knode .gt. 0 ) then
          elzaux(1) = (prod - cons)*vol       !nem
          elzaux(2) = (ddin1+ddin2)*vol       !ddin
          elzaux(3) = denit*vol               !denitrificazione

          do i=1,nlzstate
            elz(knode,i) = elzaux(i)
          end do

          return
        end if

c-----------------------------------------------------------
c diagnostics and write
c-----------------------------------------------------------

        do i=1,nlzstate
          call scalmass(elz(1,i),depth,tlztot(i))   !mass ctrl loicz
        end do

        if( next_output_d(da_out) ) then
	  call get_act_dtime(dtime)
          do i=1,nlzstate
	    ivar = 95 + i
            call scalar_output_write(dtime,da_out,ivar,1,elz(1,i))
          end do
	end if

        call lcz_av_shell(elz)          !aver/min/max of nem and ddin

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

        return
   99   continue
        stop 'error stop loicz1: error opening file'
	end

c********************************************************************

	subroutine sed_av_shell(es)

c computes and writes average/min/max of sed variables
c
c id = 360
c
c es1) average	== 361
c es1) min	== 362
c es1) max	== 363
c es2) average	== 364
c ...

	use basin

	implicit none

c parameter

	include 'param.h'

	integer nsstate
	parameter( nsstate = 2 )

	real es(nkndi,nsstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: sedacu(:,:)
	real, save, allocatable :: sedmin(:,:)
	real, save, allocatable :: sedmax(:,:)

	integer ivect(8)
	save ivect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nsstate

	  allocate(sedacu(nkn,nsstate))
	  allocate(sedmin(nkn,nsstate))
	  allocate(sedmax(nkn,nsstate))

	  id = 360
	  call cmed_init('sdv',id,nvar,1,idtc,itmc
     +				,sedacu,sedmin,sedmax,ivect)

	  icall = 1
	end if

	call cmed_accum(1,es,sedacu,sedmin,sedmax,ivect)

	end

c********************************************************************

	subroutine lcz_av_shell(elz)

c computes and writes average/min/max of bio variables
c
c id = 460
c
c elz) average	== 461
c elz) min	== 462
c elz) max	== 463
c elz) average	== 464
c ...

	use basin

	implicit none

c parameter

	include 'param.h'

	integer nlzstate
	parameter( nlzstate = 3 )

	real elz(nkndi,nlzstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: lczacu(:,:)
	real, save, allocatable :: lczmin(:,:)
	real, save, allocatable :: lczmax(:,:)

	integer ivect(8)
	save ivect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nlzstate

	  allocate(lczacu(nkn,nlzstate))
	  allocate(lczmin(nkn,nlzstate))
	  allocate(lczmax(nkn,nlzstate))

	  write (6,*) 'cmed loicz inizializzato'
	  id = 460
	  call cmed_init('lzv',id,nvar,1,idtc,itmc
     +				,lczacu,lczmin,lczmax,ivect)

	  icall = 1
	end if

	!write(6,*) 'lzc: ',id,ivect(1)

	call cmed_accum(1,elz,lczacu,lczmin,lczmax,ivect)

	end

c********************************************************************

        subroutine setsedload(nlvddi,nknddi,nstate,eload,elini)

c sets up sediment loading

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nlvddi,nknddi,nstate
        real eload(nlvddi,nknddi,nstate)
        real elini(nstate)

	include 'param.h'

        integer mode,i,k,l,lmax
        real d,vol,area

        mode = -1

        do i=1,nstate
         do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            eload(l,k,i) = 0.
          end do
          call dvanode(lmax,k,mode,d,vol,area)   !gets depth, volume and area
          eload(lmax,k,i) = elini(i) * area
         end do
        end do

        end

c********************************************************************

        subroutine check_es(es)

	use basin

        include 'param.h'

        integer nsstate
        parameter( nsstate = 2 )

        real es(nkndi,nsstate)         !sediment state variables

        integer k,i

        i = 1

        k=801
        write(6,*) 'es: ',k,i,es(k,i)

        k=1201
        write(6,*) 'es: ',k,i,es(k,i)

        end

c********************************************************************

