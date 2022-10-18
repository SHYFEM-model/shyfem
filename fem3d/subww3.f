!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

! routines for coupling WW3 model
!
! revision log :
!
! 04.07.2019    ggu     written from scratch
! 24.03.2022    ggu     newly started with Aron
! 05.05.2022    aar     lots of changes from Aron
! 15.10.2022    ggu     shympi_exchange_array substituted with shympi_l2g_array

!===========================================================
	module subww3
!===========================================================

		use shympi
		use shympi_aux
		use mod_meteo
		use mod_waves
		use mod_roughness
		use mod_depth 
		use mod_hydro
		use levels, only : nlvdi,nlv
		use basin, only : nkn,nel,ngr,mbw,xgv,ygv

		use wminitmd
		use w3cplshyfem, ONLY: SXX3DLWW3, SXY3DLWW3, SYY3DLWW3, OUTVARWW3
		use w3cplshyfem, ONLY: w3cplparam, w3cplrdstr2d, nvarsww3
		use w3gdatmd, ONLY: nseal
		use w3odatmd, ONLY: IAPROC, NAPROC
		use w3idatmd, ONLY: WX0, WXN, WY0, WYN, CX0, CXN, CY0, CYN, WLEV
		use w3idatmd, ONLY: FLCUR, FLLEV
		use yowNodepool, only: np, npa, iplg, ipgl, np_global
		use yowDatapool, only: rtype, istatus, myrank

		implicit none
		logical, save :: bww3 = .true.
		integer       :: mpiComm = -99
		integer, allocatable  :: nwild_i(:), nwild_gbi(:)
		real, allocatable     :: nwild_gb(:) 
		integer, allocatable  :: nwild_i_shyfem(:), nwild_gbi_shyfem(:)
		real, allocatable     :: nwild_gb_shyfem(:)

		double precision, allocatable :: stokesx(:,:) !stokes velocity x
		double precision, allocatable :: stokesy(:,:) !stokes velocity y

		real, allocatable       :: wavejb(:)         !wave pressure
		real, allocatable       :: wtauw(:)          !wave supported stress
		real, allocatable       :: wavewl(:)         !wave lenght
		real, allocatable       :: wavedi(:)         !wave dissipation
		real, allocatable       :: sxx3dlshyfem(:,:) !radiation stress xx
		real, allocatable       :: sxy3dlshyfem(:,:) !radiation stress yy
		real, allocatable       :: syy3dlshyfem(:,:) !radiation stress xy
		real, allocatable       :: sxx3dgl(:,:)      !radiation stress xx
		real, allocatable       :: sxy3dgl(:,:)      !radiation stress yy
		real, allocatable       :: syy3dgl(:,:)      !radiation stress xy
		real, allocatable       :: currxgl(:)        
		real, allocatable       :: currygl(:)
		real, allocatable       :: wlv(:)
		real, allocatable       :: wlvgl(:)
		real, allocatable       :: wxvgl(:)
		real, allocatable       :: wyvgl(:)
		real, allocatable       :: xgvgl(:)
		real, allocatable       :: ygvgl(:)
		real, allocatable       :: currxlshyfem(:)
		real, allocatable       :: currylshyfem(:)
		real, allocatable       :: outvargl(:,:)
		real, allocatable       :: outvarshyfem(:,:)
		integer, allocatable    :: ipgl_shyfem(:)
		integer, allocatable    :: iplg_shyfem(:)
		integer, allocatable    :: node2domain(:)

		logical, parameter      :: loutxfn = .false. 
!===========================================================

!***********************************************************
!
!***********************************************************
	contains 

	subroutine ww3_init
	USE w3cplshyfem 

	implicit none

	integer ier, idsi, idso, idss, idst, idse
	integer istat, ip, ierr, iproc, j
	logical lexist, has_output_d
	integer iwave, nvar, id
	character(len=30) :: ifname = 'ww3_multi.nml'

	real getpar

	iwave = nint(getpar('iwave'))
	bww3 = ( iwave == 11 )
	if( .not. bww3 ) return

	inquire(file=ifname,exist=lexist)
	if(.not. lexist) then
		ifname = 'ww3_multi.inp'
	endif

!AR: todo: set proper file handles, same as ww3 and check if they are free in shyfem 
	  idsi = 8
	  idso = 9
	  idss = 6
	  idst = 10
	  idse = 6
!WW3 original handles: 8           9           6          10           6
	  mpiComm = MPI_COMM_WORLD 

	  if ( trim(ifname).eq.'ww3_multi.nml' ) then
		call wminitnml ( idsi, idso, idss, idst, idse, trim(ifname), 
     1                    mpiComm, './' )
		else

		call wminit ( idsi, idso, idss, idst, idse, 
     1             'ww3_multi.inp', mpiComm, './')
		endif
!
!AR: Init coupling part of shyfem 
!
		call w3cplinit(nlv_global)

       !write(*,*) 'NSEAL =', NSEAL, 'NKN =', NKN
       !write(*,*) 'NKNDIM =', NKNDIM, 'NGLOBAL =', NP_GLOBAL

		allocate(stokesx(nlv,nkn)); stokesx = 0.
		allocate(stokesy(nlv,nkn)); stokesy = 0.
		allocate(wavejb(nkn)); wavejb = 0.
		allocate(wtauw(nkn)); wtauw = 0.
		allocate(wavewl(nkn)); wavewl = 0.
		allocate(wavedi(nkn)); wavedi = 0.

		allocate(sxx3dlshyfem(nlv_global,nkn)); sxx3dlshyfem = 0.
		allocate(sxy3dlshyfem(nlv_global,nkn)); sxy3dlshyfem = 0.
		allocate(syy3dlshyfem(nlv_global,nkn)); syy3dlshyfem = 0.
		allocate(sxx3dgl(nlv_global,np_global)); sxx3dgl = 0.
		allocate(sxy3dgl(nlv_global,np_global)); sxy3dgl = 0.
		allocate(syy3dgl(nlv_global,np_global)); syy3dgl = 0.
		allocate(outvargl(nvarsww3,np_global)); outvargl = 0.
		allocate(outvarshyfem(nvarsww3, nkn)); outvarshyfem = 0.
		allocate(wxvgl(nkn_global)); wxvgl = 0.
		allocate(wyvgl(nkn_global)); wyvgl = 0.
		allocate(currxlshyfem(nkn)); currxlshyfem = 0.
		allocate(currylshyfem(nkn)); currylshyfem = 0.
		allocate(wlv(nkn)); wlv = 0.
		allocate(xgvgl(nkn_global)); xgvgl = 0.
		allocate(ygvgl(nkn_global)); ygvgl = 0.
		allocate(wlvgl(nkn_global)); wlvgl = 0.
		allocate(currxgl(nkn_global)); currxgl = 0.
		allocate(currygl(nkn_global)); currygl = 0.
		allocate(ipgl_shyfem(nkn_global)); ipgl_shyfem = 0
		allocate(iplg_shyfem(nkn)); iplg_shyfem = 0
		allocate(node2domain(nkn_global)); node2domain = 0
!
! ... get the ghost hallo sums for all nodes ...
!
		ALLOCATE(nwild_i(np_global), nwild_gbi(np_global), stat=istat)
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		nwild_i = 0
		nwild_gbi = 0
		DO IP = 1, NSEAL
			nwild_i(IPLG(IP)) = 1
		END DO
		call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,MPI_INT,
     1                  MPI_SUM, 0, mpiComm, ierr)

		ALLOCATE(nwild_gb(NP_GLOBAL), stat=istat); nwild_gb = 0
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		IF (iaproc.eq.1) THEN
         !write(*,*) 'test 1 iaproc=', iaproc 
			DO IP = 1, NP_GLOBAL
				nwild_gb(IP)=1./REAL(nwild_gbi(IP))
			END DO
		DO iProc=2,naproc
			CALL MPI_SEND(nwild_gb,np_global,MPI_REAL, 
     1                    iProc-1, 2037, mpiComm, ierr)
		END DO
	ELSE
         !write(*,*) 'test 2 iaproc=', iaproc 
		CALL MPI_RECV(nwild_gb,np_global,MPI_REAL, 
     1                 0, 2037, mpiComm, istatus, ierr)
	END IF
	DEALLOCATE(nwild_i, nwild_gbi)
!
! ... now the same for shyfem ...
!
	ALLOCATE(nwild_i_shyfem(nkn_global), 
     1          nwild_gbi_shyfem(nkn_global), stat=istat)
	IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		nwild_i_shyfem = 0
		nwild_gbi_shyfem = 0
		DO IP = 1, NKN
			nwild_i_shyfem(ip_int_node(IP)) = 1
		END DO
		call mpi_reduce(nwild_i_shyfem,nwild_gbi_shyfem,
     1                  nkn_global,MPI_INT,
     1                  MPI_SUM, 0, mpiComm, ierr)

		ALLOCATE(nwild_gb_shyfem(nkn_global), stat=istat)
		nwild_gb_shyfem = 0
		IF (istat/=0) stop 'allocate error nwild_i,nwild_gbi'
		IF (iaproc.eq.1) THEN
         !write(*,*) 'test 1 iaproc=', iaproc 
			DO IP=1,nkn_GLOBAL
				nwild_gb_shyfem(IP)=1./REAL(nwild_gbi_shyfem(IP))
			END DO
			DO iProc=2,naproc
				CALL MPI_SEND(nwild_gb_shyfem,nkn_global,MPI_REAL,
     1                    iProc-1, 2037, mpiComm, ierr)
			END DO
		ELSE
         !write(*,*) 'test 2 iaproc=', iaproc 
			CALL MPI_RECV(nwild_gb_shyfem,nkn_global,MPI_REAL,
     1                 0, 2037, mpiComm, istatus, ierr)
		END IF
		DEALLOCATE(nwild_i_shyfem, nwild_gbi_shyfem)

! now compute some basic mappings for shyfem 

		do ip = 1, nkn_inner
			node2domain(ip_int_node(ip)) = myrank + 1
         !write(*,*) 'myrank =,', myrank + 1, ip, ip_int_node(ip)
		enddo 

		j = 1
		do ip = 1, nkn_global
			if ( node2domain(ip) == myrank+1 ) then
				iplg_shyfem(j) = ip
				ipgl_shyfem(IP) = j
			j = j + 1
			endif
		end do
!
! init coupling
!
		!if (myrank == 0) write(*,*) 'calling getvarww3 on init' 
		call getvarww3(.true.)
		call getvarshyfem(.true.)

		nvar = 3
		call init_output_d('itmwav','idtwav',da_wav)

		if( has_output_d(da_wav)) then
			call shyfem_init_scalar_file('wave',nvar,.true.,id)
			da_wav(4) = id
		end if

!***********************************************************
	end subroutine ww3_init
!***********************************************************
!
!***********************************************************
	subroutine ww3_loop

	use basin
	use mod_hydro
	use mod_meteo
      use WMWAVEMD 

	implicit none

	integer n,nsave, tend(2,1), yy,mm,dd,h,m,s, it
	integer ys(8), ierr 

	double precision atime, dtime
	double precision daux

	if( .not. bww3 ) return

	call get_act_dtime(dtime)
	it = dtime
	call convert_time_d('dtwave',daux)
	idcoup = nint(daux)

	if (mod(it,idcoup) .eq. 0 ) then

	call get_absolute_act_time(atime)
	call dts_from_abs_time_to_ys(atime,ys)

	yy = ys(1)
	mm = ys(2)
	dd = ys(3)
	h  = ys(4)
	m  = ys(5)
	s  = ys(6)
 
!       write(*,*) 'ys 1', ys

	call dts_from_abs_time_to_ys(atime+daux,ys)

	yy = ys(1)
	mm = ys(2)
	dd = ys(3)
	h  = ys(4)
	m  = ys(5)
	s  = ys(6)

!       write(*,*) 'ys 2', ys

! Set end time of this advance
	tend(1,1) = 10000*yy + 100*mm + dd
	tend(2,1) = 10000*h  + 100*m  + s

	if( .not. bww3 ) return

		!WRITE(*,'(a20,9I10,F20.10,F30.10)') 'ww3 running now for', 
    ! 1                           idcoup, ys, daux, atime
		!WRITE(*,*) 'tend', tend
		call mpi_barrier(mpiComm,ierr)
		call getvarshyfem(.false.) 
		call mpi_barrier(mpiComm,ierr)
		!if (myrank == 0) WRITE(*,*) 'calling wmwave'
		call wmwave ( tend )
		call mpi_barrier(mpiComm,ierr)
		!if (myrank == 0) WRITE(*,*) 'calling getvar' 
		call getvarww3(.false.)
		call mpi_barrier(mpiComm,ierr)
        endif ! idcoup
        
	end subroutine ww3_loop
!***********************************************************
!
!***********************************************************

	subroutine ww3_finalize

	use basin
        use WMFINLMD

	implicit none

	if( .not. bww3 ) return

! here finalize ww3
      CALL WMFINL

	end subroutine ww3_finalize

!***********************************************************
!
!***********************************************************
!
! $Id: subwwm.f,v 1.1 2006/10/18 15:35:13 georg Exp $
!
! fifo pipe routines for WWM wave model
!
! routines for reading and writing data from the SHYFEM and the WWM model
! througth the FIFO PIPE mechanism
! file format: binary file with each value per node
!
!**************************************************************

!**************************************************************
		subroutine getvarww3(lfirst) 

		implicit none

c common
		include 'param.h'
		include 'pkonst.h'

c local
		logical, intent(in) :: lfirst 
		integer k,l, ie
		real wtau,taux,tauy   !wave and wind stresses
		real s,d      !speed and direction
		real pi,deg2rad
		parameter ( pi=3.14159265358979323846, deg2rad = pi/180. )
		double precision tmpval
		integer itdrag, ilev, ivar
		integer it, ip
		logical bwind
		save bwind
		real wfact,wspeed
		save wfact
		real getpar
		logical next_output_d
		integer id, nvar, ierr
		double precision dtime
		double precision daux
		integer icall             !initialization parameter                        
		save icall
		data icall /0/
		real tramp,alpha
		save tramp
		logical has_output_d

!------------------------------------------------------

		if (.not. bww3) return

!       -----------------------------------------------
!       Same time step, do read
!       -----------------------------------------------

		wtauw = 0.

		call get_act_dtime(dtime)

		it = dtime

!-------------------------------------------------------------
! set coupling time step 
!-------------------------------------------------------------

		call convert_time_d('dtwave',daux)
		idcoup = nint(daux)

		if (mod(it,idcoup) .eq. 0 ) then

!         -----------------------------------------------
!         compute stress and wave characteristics
!         -----------------------------------------------
      !if (myrank == 0) write(*,*) 'wave coupling param.', bww3, iwave 
			call w3cplrdstr2d
			call w3cplparam

!         -----------------------------------------------
!         compute global arrays from ww3
!         -----------------------------------------------

!         radiation stress ...

			do ilev = 1, nlv_global
				call get_global_array_ww3(SXX3DLWW3(ilev,:),SXX3DGL(ilev,:))
				call get_global_array_ww3(SXY3DLWW3(ilev,:),SXY3DGL(ilev,:))
				call get_global_array_ww3(SYY3DLWW3(ilev,:),SYY3DGL(ilev,:))
			enddo

			!write(*,*) 'rad 1', sum(SXX3DLWW3), SUM(SXY3DLWW3), SUM(SYY3DLWW3)

!         integral wave parameter ... 

			do ivar = 1, nvarsww3
				call get_global_array_ww3(outvarww3(ivar,:),
     1           outvargl(ivar,:)) 
			enddo 

!         -----------------------------------------------
!         compute local arrays for shyfem
!         -----------------------------------------------

!         radiation stress ...

	do ilev = 1, nlv_global
		call fill_local_array_shyfem(SXX3DLSHYFEM(ilev,:)
     +			,SXX3DGL(ilev,:))
		call fill_local_array_shyfem(SXY3DLSHYFEM(ilev,:)
     +			,SXY3DGL(ilev,:))
		call fill_local_array_shyfem(SYY3DLSHYFEM(ilev,:)
     +			,SYY3DGL(ilev,:))
	enddo
			call mpi_barrier(mpicomm,ierr) 

			!write(*,*) 'rad 2', myrank, sum(SXX3DLSHYFEM), SUM(SXX3DLSHYFEM), 
!     1                    SUM(SXX3DLSHYFEM)

!         integral wave parameter ... 

	do ivar = 1, nvarsww3
		call fill_local_array_shyfem(outvarshyfem(ivar,:),
     +           outvargl(ivar,:))
	enddo
			call mpi_barrier(mpicomm,ierr) 
			waveh = outvarshyfem(1,:)
			wavep = outvarshyfem(2,:)
			waved = outvarshyfem(3,:)

!         -----------------------------------------------
!         Compute surface roughness z0s = 0.5*Hs
!         -----------------------------------------------

			do k = 1,nkn
				z0sk(k) = 0.5 * waveh(k)
			end do

!         -----------------------------------------------
!         Computes wave induced forces
!         -----------------------------------------------

			if (iwave .eq. 11) then

!           -----------------------------------------------
!           Radiation stress formulation
!           -----------------------------------------------

				call diffxy(sxx3dlshyfem,sxy3dlshyfem
     +					,syy3dlshyfem,wavefx,wavefy)
				!write(*,*) 'sum rad', sum(sxx3dlshyfem), sum(sxy3dlshyfem)
				!write(*,*) 'sum wave fx', sum(wavefx), sum(wavefy) 

			elseif (iwave .eq. 12) then
!           -----------------------------------------------
!           Vortex force formulation and substract wave-supported
!           stress from wind stress (still to be implemented)
!           -----------------------------------------------
!AR: todo: compute stokes velocities ...

				call wave_vortex(stokesx,stokesy,wavejb,wavefx,wavefy)

				do k = 1,nkn
					wtau = wtauw(k)
					taux = tauxnv(k)
					tauy = tauynv(k)
					s = sqrt(taux**2 + tauy**2)
					if (s .gt. wtau) then
						d = atan2(tauy,taux)
						tauxnv(k) = taux - wtau * cos(d)
						tauynv(k) = tauy - wtau * sin(d)
					end if
				end do
			endif ! iwave
!         -----------------------------------------------
!         simulate smooth initial forcing
!         useful for test cases
!         -----------------------------------------------
			alpha = 1.
			if( tramp .gt. 0. ) then
				call get_passed_dtime(dtime)
				alpha = dtime/tramp
				if( alpha .gt. 1. ) alpha = 1.
			end if
			if( alpha /= 1. ) then
			do ie = 1,nel
				do l = 1,nlv
				wavefx(l,ie) = wavefx(l,ie) * alpha
				wavefy(l,ie) = wavefy(l,ie) * alpha
				end do
			end do
			end if

		end if ! dtcoup

!     -----------------------------------------------
!       Writes output to the file.wav 
!     -----------------------------------------------

		if( next_output_d(da_wav) ) then
			id = nint(da_wav(4))
			call get_act_dtime(dtime)
			call shy_write_scalar_record2d(id,dtime,231,waveh)
			call shy_write_scalar_record2d(id,dtime,232,wavep)
			call shy_write_scalar_record2d(id,dtime,233,waved)
			call shy_sync(id)
		end if

		if (myrank == 0 .and. loutxfn) then
			if (lfirst) then
				open(10002, file = 'ergwav.bin', form = 'unformatted')
				write(10002)   SNGL(0.)
				write(10002)  (SNGL(outvargl(2,ip)), SNGL(outvargl(3,ip)),
     1                 SNGL(outvargl(1,ip)), ip = 1, np_global)
			else
				write(10002)   SNGL(float(it))
				write(10002)  (SNGL(outvargl(2,ip)), SNGL(outvargl(3,ip)),
     1                 SNGL(outvargl(1,ip)), ip = 1, np_global)
			endif

			if (lfirst) then
				open(10003, file = 'ergrad.bin', form = 'unformatted')
				write(10003)   SNGL(0.)
				write(10003)  (SNGL(SXX3DGL(1,ip)), SNGL(SXY3DGL(1,ip)),
     1                 SNGL(SYY3DGL(1,ip)), ip = 1, np_global)
			else
				write(10003)   SNGL(float(it))
				write(10003)  (SNGL(SXX3DGL(1,ip)), SNGL(SXY3DGL(1,ip)),
     1                 SNGL(SYY3DGL(1,ip)), ip = 1, np_global)
			endif

		endif

	end subroutine getvarww3
!***********************************************************
!
!***********************************************************
		subroutine getvarshyfem(lfirst) 

		implicit none

c common
		include 'param.h'
		include 'pkonst.h'

c local

		logical, intent(in) :: lfirst
		integer k,l, ie
		real wtau,taux,tauy   !wave and wind stresses
		real s,d      !speed and direction
		real pi,deg2rad
		parameter ( pi=3.14159265358979323846, deg2rad = pi/180. )
		double precision tmpval
		integer itdrag, ilev, ivar
		integer it, ip
		logical bwind
		save bwind
		real wfact,wspeed
		save wfact
		real getpar
		logical next_output_d
		integer id, nvar
		double precision dtime
		double precision daux
		integer icall             !initialization parameter                        
		save icall
		data icall /0/
		real tramp,alpha
		save tramp
		logical has_output_d

!------------------------------------------------------

		if (.not. bww3) return

!       -----------------------------------------------
!       Same time step, do read
!       -----------------------------------------------

		wtauw = 0.

		call get_act_dtime(dtime)

		it = dtime

!-------------------------------------------------------------
! set coupling time step 
!-------------------------------------------------------------

		call convert_time_d('dtwave',daux)
		idcoup = nint(daux)

		if (mod(it,idcoup) .eq. 0 ) then
!
!         -----------------------------------------------
!         compute global arrays from shyfem
!         -----------------------------------------------
!
			if (.not. lfirst) then 
				wx0 = wxn
				wy0 = wyn 
			endif 
!         wind 
			call get_global_array_shyfem(wxv,wxvgl)
			call get_global_array_shyfem(wyv,wyvgl)
			call get_global_array_shyfem(xgv,xgvgl)
			call get_global_array_shyfem(ygv,ygvgl)

			if (lfirst) then 
				wx0(:,1) = wxvgl 
				wy0(:,1) = wyvgl
				wxn(:,1) = wxvgl
				wyn(:,1) = wyvgl
			else
				wxn(:,1) = wxvgl
				wyn(:,1) = wyvgl
			endif 
			!write(*,*) 'wind sums', sum(wyv), sum(wyvgl), sum(wy0), sum(wyn) 

!         curr 
			do k = 1,nkn
				call getuv(1,k,currxlshyfem(k),currylshyfem(k))
			enddo 

			if (myrank == 0) then 
				write(*,*) 'currxy flag', flcur
				write(*,*) 'water level flag', fllev
			endif

			if (flcur) then 
				if (.not. lfirst) then    
					cx0 = cxn
					cy0 = cyn    
				endif 

				call get_global_array_shyfem(currxlshyfem, currxgl)
				call get_global_array_shyfem(currylshyfem, currygl)

				if (lfirst) then 
					cx0(:,1) = currxgl
					cy0(:,1) = currygl
					cxn(:,1) = currxgl
					cyn(:,1) = currygl
				else
					cxn(:,1) = currxgl
					cyn(:,1) = currygl
				endif 
			if (myrank == 0) then
				write(*,*) 'mean currents cx', sum(currxgl)/real(size(currxgl))
				write(*,*) 'mean currents cy', sum(currygl)/real(size(currygl))
				write(*,*) 'max cx, cy', maxval(currxgl), maxval(currygl)
				write(*,*) 'min cx, cy', minval(currxgl), minval(currygl)
			endif
		endif ! flcur

!         water level
		if (FLLEV) then 
			wlv = znv 
			call get_global_array_shyfem(wlv,wlvgl)
			wlev(:,1) = wlvgl
			if (myrank == 0) then
				write(*,*) 'mean wlv', sum(wlvgl)/real(size(wlvgl))
				write(*,*) 'max wlv', maxval(wlvgl)
				write(*,*) 'min wlv', minval(wlvgl)
				do k = 1, nkn
					write(10010,*) k, xgvgl(k), ygvgl(k)
				enddo 
			endif
		endif

		end if ! dtcoup

		if (myrank == 0 .and. loutxfn) then
			if (lfirst) then
				open(10001, file = 'erguvz.bin', form = 'unformatted')
				write(10001)   SNGL(0.)
				write(10001)  (SNGL(currxgl(ip)), SNGL(currygl(ip)),
     1                 SNGL(wlvgl(ip)), ip = 1, np_global)
			else
				write(10001)   SNGL(float(it))
				write(10001)  (SNGL(currxgl(ip)), SNGL(currygl(ip)),
     1                 SNGL(wlvgl(ip)), ip = 1, np_global)
			endif
		endif


	end subroutine getvarshyfem
!***********************************************************
!
!***********************************************************
!differenzation of radiation stresses
!**************************************************************************
! Computes wave forcing terms according to the radiation stress formulation

	subroutine diffxy(SXX3D,SYY3D,SXY3D,wavefx,wavefy)

	use evgeom
	use levels
	use basin

	implicit none

	!include 'param.h'

	real, intent(in) :: SXX3D(:,:)       !radiation stress xx
	real, intent(in) :: SYY3D(:,:)       !radiation stress yy
	real, intent(in) :: SXY3D(:,:)       !radiation stress xy

	real,intent(out) ::  wavefx(nlv,nel)      !wave forcing term x
	real,intent(out) ::  wavefy(nlv,nel)      !wave forcing term y

	double precision :: b,c           !x and y derivated form function [1/m]
	integer          :: k,ie,ii,l,ilevel
	real             :: radsx,radsy

	wavefx = 0.
	wavefy = 0.
	do ie = 1,nel
		ilevel = ilhv(ie)
		do l=1,ilevel
			radsx = 0.
			radsy = 0.
			do ii = 1,3
				k = nen3v(ii,ie)
				b = ev(3+ii,ie)
				c = ev(6+ii,ie)
				radsx = radsx -(SXX3D(l,k)*b + SXY3D(l,k)*c)
				radsy = radsy -(SXY3D(l,k)*b + SYY3D(l,k)*c)
			end do
			wavefx(l,ie) = -radsx
			wavefy(l,ie) = -radsy
		enddo
	enddo

	end subroutine diffxy
!**************************************************************
!
!**************************************************************
	subroutine get_global_array_ww3(localarray,globalarray)

		implicit none

		real, intent(in)  :: localarray(:)
		real, intent(out) :: globalarray(:)
		real, allocatable :: tmp(:)
		integer           :: ip, ierr

		allocate(tmp(NP_GLOBAL)); tmp = 0.
		globalarray = 0. 

		if (size(globalarray).ne. NP_GLOBAL) then
			write(*,*) size(globalarray), NP_GLOBAL
			stop 'error in grid sizes get_global_array_ww3' 
		endif

		DO IP = 1, NSEAL 
			tmp(iplg(IP)) = localarray(IP)
		END DO

		call mpi_allreduce(tmp,globalarray,NP_GLOBAL,
     &                   MPI_REAL,MPI_SUM,mpiComm,ierr)

!         if(myrank==0) then
		do IP = 1, NP_GLOBAL
			globalarray(IP) = globalarray(IP)*nwild_gb(IP)
		enddo !IP
!         endif !myrank

	end subroutine get_global_array_ww3
!**************************************************************
!
!**************************************************************
	subroutine get_global_array_shyfem(localarray,globalarray)

		implicit none

		real, intent(in)  :: localarray(:)
		real, intent(out) :: globalarray(:)

		integer           :: ip, ierr

		if (size(globalarray).ne. nkn_GLOBAL) then
			write(*,*) size(globalarray), nkn_GLOBAL
			stop 'error in grid sizes get_global_array_shyfem'
		endif

		call shympi_l2g_array(localarray,globalarray)

	end subroutine get_global_array_shyfem
!**************************************************************
!
!**************************************************************
	subroutine fill_local_array_ww3(localarray,globalarray)

		implicit none

		real, intent(out)  :: localarray(:)
		real, intent(in)   :: globalarray(:)
		real, allocatable :: tmp(:)

		integer           :: ip, ierr

		allocate(tmp(NP_GLOBAL)); tmp = 0.
		localarray = 0.

		if (size(globalarray).ne. NP_GLOBAL) then
			write(*,*) size(globalarray), NP_GLOBAL
			stop 'error in grid sizes fill_local_array_ww3'
		endif

		DO IP = 1, NP_GLOBAL
			IF (ipgl(ip) .gt. 0) THEN
				localarray(ipgl(ip)) = globalarray(IP)
			END IF
		ENDDO

	end subroutine fill_local_array_ww3
!**************************************************************
!
!**************************************************************
	subroutine fill_local_array_shyfem(localarray,globalarray)

		implicit none

		real, intent(out) :: localarray(:)
		real, intent(in)  :: globalarray(:)

		integer           :: ip, ierr

		localarray = 0.

		if (size(globalarray).ne. nkn_GLOBAL) then
			write(*,*) size(globalarray), nkn_GLOBAL
			stop 'error in grid sizes fill_local_array_shyfem'
		endif

		call shympi_g2l_array(globalarray,localarray)

	end subroutine fill_local_array_shyfem
!**************************************************************************
	end module
!**************************************************************************
