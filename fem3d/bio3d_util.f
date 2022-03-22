
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
c 22.03.2022    ggu     upgraded to da_out and new cmed
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
        integer i,k,idc,id,is,nvar
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
	  nvar = nlzstate
          call shyfem_init_scalar_file('lcz',nvar,.true.,id)
          if( id <= 0 ) icall = -1
          if( icall < 0 ) return
	  da_out(4) = id

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
	  idc = 95
          id = nint(da_out(4))
	  call get_act_dtime(dtime)
          do is=1,nlzstate
	    idc = idc + 1
            call shy_write_scalar_record(id,dtime,idc,1,elz(:,is))
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

	integer nsstate
	parameter( nsstate = 2 )

	real es(nkn,nsstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar,idc,is,nstate
	double precision :: dtime,rr
c function
	real getpar
        logical has_output_d,is_over_output_d,next_output_d
c save
	double precision, save :: da_out(4)
	double precision, save, allocatable :: sedacu(:,:)
	real, save, allocatable :: sedmin(:,:)
	real, save, allocatable :: sedmax(:,:)
	real, save, allocatable :: raux(:)

        integer, save :: nr = 0
        integer, save :: icall = 0

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

          nstate = nsstate

          call init_output_d('itmcon','idtcon',da_out)
          call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = nstate*3
            call shyfem_init_scalar_file('sdv',nvar,.true.,id)
            da_out(4) = id
          end if

	  allocate(sedacu(nkn,nsstate))
	  allocate(sedmin(nkn,nsstate))
	  allocate(sedmax(nkn,nsstate))
	  allocate(raux(nkn))

          do is=1,nstate
            call cmed_reset_2d(nr,sedacu(:,is)
     +                  ,sedmin(:,is),sedmax(:,is))
          end do

          write (6,*) 'cmed bio sediments inizializzato'

	  icall = 1
	end if

        if( .not. is_over_output_d(da_out) ) return

        nstate = nsstate

        nr = nr + 1
        do is=1,nstate
          call cmed_accum_2d(es(:,is),sedacu(:,is)
     +                  ,sedmin(:,is),sedmax(:,is))
        end do

        if( .not. next_output_d(da_out) ) return

        id = nint(da_out(4))
        call get_act_dtime(dtime)
        rr=1./nr

        idc = 360
        do is=1,nstate
          idc = idc + 1
          raux = sedacu(:,is) * rr
          call shy_write_scalar_record(id,dtime,idc,1,raux)
          raux = sedmin(:,is)
          call shy_write_scalar_record(id,dtime,idc,1,raux)
          raux = sedmax(:,is)
          call shy_write_scalar_record(id,dtime,idc,1,raux)
        end do

        do is=1,nstate
          call cmed_reset_2d(nr,sedacu(:,is)
     +                  ,sedmin(:,is),sedmax(:,is))
        end do

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

	integer nlzstate
	parameter( nlzstate = 3 )

	real elz(nkndi,nlzstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar,nstate,is,idc
	double precision :: dtime,rr
c function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
c save
	double precision, save :: da_out(4)
	double precision, save, allocatable :: lczacu(:,:)
	real, save, allocatable :: lczmin(:,:)
	real, save, allocatable :: lczmax(:,:)
	real, save, allocatable :: raux(:)

	integer, save :: nr = 0
	integer, save :: icall = 0

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  nstate = nlzstate

          call init_output_d('itmcon','idtcon',da_out)
          call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = nstate*3
            call shyfem_init_scalar_file('lzv',nvar,.true.,id)
            da_out(4) = id
          end if

	  allocate(lczacu(nkn,nstate))
	  allocate(lczmin(nkn,nstate))
	  allocate(lczmax(nkn,nstate))
	  allocate(raux(nkn))

          do is=1,nstate
            call cmed_reset_2d(nr,lczacu(:,is)
     +                  ,lczmin(:,is),lczmax(:,is))
          end do

	  write (6,*) 'cmed loicz inizializzato'

	  icall = 1
	end if

        if( .not. is_over_output_d(da_out) ) return

	nstate = nlzstate

        nr = nr + 1
        do is=1,nstate
          call cmed_accum_2d(elz(:,is),lczacu(:,is)
     +                  ,lczmin(:,is),lczmax(:,is))
        end do

        if( .not. next_output_d(da_out) ) return

        id = nint(da_out(4))
        call get_act_dtime(dtime)
        rr=1./nr

        idc = 460
        do is=1,nstate
          idc = idc + 1
          raux = lczacu(:,is) * rr
          call shy_write_scalar_record(id,dtime,idc,1,raux)
          raux = lczmin(:,is)
          call shy_write_scalar_record(id,dtime,idc,1,raux)
          raux = lczmax(:,is)
          call shy_write_scalar_record(id,dtime,idc,1,raux)
        end do

        do is=1,nstate
          call cmed_reset_2d(nr,lczacu(:,is)
     +                  ,lczmin(:,is),lczmax(:,is))
        end do

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

