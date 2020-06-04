
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2011,2014-2015,2017-2019  Georg Umgiesser
!    Copyright (C) 2015  Debora Bellafiore
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

! revision log :
!
! 12.01.2009	ggu	new file for system routines
! 31.03.2009	ggu	call renamed to spk_*
! 23.03.2010	ggu	changed v6.1.1
! 22.11.2011	ggu	changed VERS_6_1_37
! 28.01.2014	ggu	changed VERS_6_1_71
! 25.05.2015	ggu	some calls changed (pass array in)
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_52
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 12.10.2015	ggu	changed VERS_7_3_3
! 09.12.2015	ggu	adapted to new pointers and 3d matrix
! 15.12.2015	ggu&dbf	finsihed and validated
! 05.12.2017	ggu	changed VERS_7_5_39
! 23.04.2018	ggu	adapted to new matrix type (local and global matrix)
! 11.05.2018	ggu	changed VERS_7_5_47
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 04.06.2020	ggu	bug in system_assemble() - use kn to decide on assemble
!
!******************************************************************

!==================================================================
	module mod_system_global
!==================================================================

	integer, save :: n2zero_max
	integer, save :: nkn_max
	integer, save, allocatable :: n2_zeros(:)
	integer, save, allocatable :: nkns(:)
	integer, save, allocatable :: i2coos(:,:)
	integer, save, allocatable :: j2coos(:,:)
	integer, save, allocatable :: ij2coos(:,:)
	!integer, save, allocatable :: ip_int_nodes(:,:)
	double precision, save, allocatable :: c2coos(:,:)
	double precision, save, allocatable :: c2rhss(:,:)
	double precision, save, allocatable :: zz(:,:)
	double precision, save, allocatable :: exchanges(:,:)

!==================================================================
	end module mod_system_global
!==================================================================

        subroutine system_initialize

! allocates data structure

	use mod_system
	use levels
	use basin
	use shympi

        implicit none

	call shympi_barrier

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines ',my_id,nkn
        write(6,*) '----------------------------------------'

        call mod_system_init(nkn,nel,ngr,mbw,nlv,l_matrix)
	call mod_system_insert_elem_index(nel,nen3v,l_matrix)
	call mod_system_set_local(l_matrix)
	call spk_initialize_system(l_matrix)		!calls coo_init_new

	call shympi_barrier

	if( bmpi ) then			!only needed if not explicit !FIXME
          call mod_system_init(nkn_global,nel_global,ngr_global
     +				,mbw,nlv,g_matrix)
	  call mod_system_insert_elem_index(nel_global,nen3v_global
     +					,g_matrix)
	  call mod_system_set_global(g_matrix)
	  call spk_initialize_system(g_matrix)		!calls coo_init_new
	  call system_setup_global
	end if

	call shympi_barrier

        end

!******************************************************************

	subroutine system_init

! initializes arrays (set to zero)
!
! must be called before every assembly of matrix

	use mod_system
	use shympi

	implicit none

	if( bsysexpl ) then
	  l_matrix%rvec2d = 0.
	  l_matrix%raux2d = 0.
	else
	  if( bmpi ) then
	    call spk_init_system(g_matrix)
	  end if
	  call spk_init_system(l_matrix)
	end if

	end

!******************************************************************

	subroutine system_set_explicit

	use mod_system

	implicit none

	bsysexpl = .true.

	end

!******************************************************************
!******************************************************************
!******************************************************************
!
!	subroutine system_show_internal
!
!	use basin
!	use mod_system
!	use shympi
!
!        write(my_unit,*) 'system internal: ',my_id,bsysexpl
!
!        do k=1,nkn
!          write(my_unit,1234) k,ipv(k),id_node(k),xgv(k),ygv(k)
!     +                          ,rvec2d(k),raux2d(k)
! 1234     format(3i5,4f12.4)
!        end do
!
!	end
!
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_setup_global

! this sets up the global system

	use mod_system
	use mod_system_global
	use shympi

	implicit none

	integer n2zero,n2zero_global
	integer id
	integer ipp,ipp0
	integer ki,kj,i,j,n,ngr,nkn
	integer, allocatable :: i2aux(:)
	integer, allocatable :: j2aux(:)
	integer, allocatable :: ipaux(:)
	integer, allocatable :: nodes(:)
	type(smatrix), pointer :: mm

	!stop 'error stop system_solve_global_z: not ready'

	mm => l_matrix

!---------------------------------------------------------------
! gather info and allocate arrays
!---------------------------------------------------------------

	allocate(n2_zeros(0:n_threads-1))
	allocate(nkns(0:n_threads-1))

	n2zero = mm%n2zero
	call shympi_gather(n2zero,n2_zeros)
	n2zero_max = maxval(n2_zeros)

	nkn = mm%nkn_system
	call shympi_gather(nkn,nkns)
	nkn_max = shympi_max(nkn_local)

	allocate(i2coos(n2zero_max,0:n_threads-1))
	allocate(j2coos(n2zero_max,0:n_threads-1))
	allocate(ij2coos(n2zero_max,0:n_threads-1))
	allocate(c2coos(n2zero_max,0:n_threads-1))
	allocate(c2rhss(nkn_max,0:n_threads-1))
	allocate(zz(nkn_max,0:n_threads-1))
	allocate(exchanges(2*nkn_max+n2zero_max,0:n_threads-1))
	!allocate(ip_int_nodes(nkn_max,0:n_threads-1))
	allocate(i2aux(n2zero_max))
	allocate(j2aux(n2zero_max))
	allocate(ipaux(nkn_max))

	i2aux = 0
	j2aux = 0

!---------------------------------------------------------------
! convert local internal to global internal
!---------------------------------------------------------------

	do i=1,n2zero
	  !write(6,*) i,mm%i2coo(i),mm%j2coo(i)
	  i2aux(i) = ip_int_node( mm%i2coo(i) )
	  j2aux(i) = ip_int_node( mm%j2coo(i) )
	end do

!---------------------------------------------------------------
! exchange pointers i/j2coos
!---------------------------------------------------------------

	call shympi_gather(i2aux,i2coos)
	call shympi_gather(j2aux,j2coos)

	!ipaux(1:nkn) = ip_int_node(1:nkn)
	!call shympi_gather(ipaux,ip_int_nodes)

!---------------------------------------------------------------
! translate to pointers into global matrix
!---------------------------------------------------------------

	mm => g_matrix

	ngr = mm%ngr_system
	n2zero_global = mm%n2zero
	allocate(nodes(ngr+1))

	do id=0,n_threads-1
	  n2zero = n2_zeros(id)
	  do i=1,n2zero
	    ki = i2coos(i,id)			!row
	    kj = j2coos(i,id)			!col
	    ipp0 = mm%ntg(ki-1)			!last entry in row ki-1
	    n = mm%ng(ki)			!total entries in row ki
	    nodes(1:n) = mm%iorder(1:n,ki)
            do j=1,n
              if( nodes(j) == kj ) exit       !find kj in nodes
            end do
            if( j > n ) goto 98
            ipp = ipp0 + j
            if( ipp > n2zero_global ) goto 97
	    if( mm%i2coo(ipp) /= ki ) goto 97
	    if( mm%j2coo(ipp) /= kj ) goto 97
	    ij2coos(i,id) = ipp
	  end do
	end do

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	if( shympi_is_master() ) then
	  write(6,*) 'finished routine system_setup_global_z'
	end if

	return
   97	continue
	write(6,*) id,i,ki,kj,n
	write(6,*) nodes(1:n)
	write(6,*) j,ipp0,ipp,n2zero_global
	write(6,*) mm%i2coo(ipp),mm%j2coo(ipp)
	stop 'error stop system_setup_global_z: internal error (2)'
   98	continue
	write(6,*) id,i,ki,kj,n
	write(6,*) nodes(1:n)
	stop 'error stop system_setup_global_z: internal error (1)'
	end

!******************************************************************

	subroutine system_solve_global(n,z)

! this solves the global system

	use mod_system
	use mod_system_global
	use mod_system_interface
	use shympi

	implicit none

	integer n
	real z(n)

	integer n2max,n2zero,n2zero_global
	integer i,id,ipp,nkn,iint,k,nc,ia
	integer ipb,ipbm,ipbr,ipbz
	double precision, allocatable :: c2aux(:)
	real, allocatable :: zglobal(:)
	double precision, allocatable :: zlocal(:)
	double precision, allocatable :: rlocal(:)
	double precision, allocatable :: exchange(:)
	type(smatrix), pointer :: mm
	real, parameter :: flag = 1.234567e+20

	mm => l_matrix
	n2zero = n2zero_max

	allocate(exchange(2*nkn_max+n2zero_max))
	allocate(zglobal(nkn_global))
	zglobal = flag

	!if( shympi_is_master() ) then
	!  write(6,*) 'starting system_solve_global_z'
	!end if

!---------------------------------------------------------------
! gather values from all processes
!---------------------------------------------------------------

	!allocate(c2aux(n2zero_max))
	!allocate(zlocal(nkn_max))
	!allocate(rlocal(nkn_max))

	!c2aux(1:n2zero) = mm%c2coo(1:n2zero)		!locally assembled
	!zlocal(1:n) = z(1:n)
	!rlocal(1:n) = mm%rvec2d(1:n)

	!call shympi_gather(c2aux,c2coos)		!all matrix values
	!call shympi_gather(zlocal,zz)			!all zeta values
	!call shympi_gather(rlocal,c2rhss)		!all rhs values

	ipb = 0
	exchange(ipb+1:ipb+n) = mm%rvec2d(1:n)
	ipb = nkn_max
	exchange(ipb+1:ipb+n) = z(1:n)		!convert from real to double
	ipb = 2*nkn_max
	exchange(ipb+1:ipb+n2zero) = mm%c2coo(1:n2zero)

	call shympi_gather(exchange,exchanges)		!all rhs values

!---------------------------------------------------------------
! switch to global matrix and assemble matrix and zeta
!---------------------------------------------------------------

	mm => g_matrix

	n2zero_global = mm%n2zero
	ipbr = 0
	ipbz = nkn_max
	ipbm = 2*nkn_max

	do id=0,n_threads-1
	  ia = id + 1
	  n2zero = n2_zeros(id)
	  do i=1,n2zero
	    ipp = ij2coos(i,id)
	    !mm%c2coo(ipp) = mm%c2coo(ipp) + c2coos(i,id)
	    mm%c2coo(ipp) = mm%c2coo(ipp) + exchanges(ipbm+i,id)
	  end do
	  nkn = nkns(id)
	  do i=1,nkn
	    !iint = ip_int_nodes(i,id)
	    iint = ip_int_nodes(i,ia)
	    !mm%rvec2d(iint) = mm%rvec2d(iint) + c2rhss(i,id)
	    !zglobal(iint) = zz(i,id)		!convert from double to real
	    mm%rvec2d(iint) = mm%rvec2d(iint) + exchanges(ipbr+i,id)
	    zglobal(iint) = exchanges(ipbz+i,id)!convert from double to real
	  end do
	end do

!---------------------------------------------------------------
! solve global matrix
!---------------------------------------------------------------

	n2max = mm%n2max				!maybe also n2zero
	nkn = mm%nkn_system
	if( nkn /= nkn_global ) goto 99
	if( any( zglobal == flag ) ) goto 99

	call spk_solve_system(g_matrix,.false.,n2max,nkn,zglobal)

!---------------------------------------------------------------
! copy solution for zeta from global to local data structure
!---------------------------------------------------------------

	do k=1,nkn_local
	  iint = ip_int_node(k)
	  l_matrix%rvec2d(k) = g_matrix%rvec2d(iint)
	end do

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	!if( shympi_is_master() ) then
	!  write(6,*) 'finished system_solve_global_z'
	!end if

	return
   99	continue
	nc = count( zglobal == flag )
	write(6,*) n2max,nkn,nkn_global,nc
	stop 'error stop system_solve_global_z: internal error (1)'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_solve(n,z)

! solves system - z is used for initial guess, not the final solution

	use mod_system
	use mod_system_interface
	use shympi

	implicit none

	integer, intent(in) :: n
	real, intent(in)    :: z(n)

	integer n2max,nu
	double precision t_start,t_end,t_passed
	type(smatrix), pointer :: mm

	mm => l_matrix
	nu = nkn_unique

	t_start = shympi_wtime()

	if( bsysexpl ) then			!explicit - solved direct
          call shympi_exchange_and_sum_2d_nodes(mm%rvec2d)
          call shympi_exchange_and_sum_2d_nodes(mm%raux2d)
	  mm%rvec2d(1:nu) = mm%rvec2d(1:nu) / mm%raux2d(1:nu)	!GGUEXPL
	else if( bmpi ) then			!mpi - solve globally
	  call system_solve_global(n,z)
	else					!solve directly locally
	  n2max = mm%n2max
	  call spk_solve_system(l_matrix,.false.,n2max,n,z)
	end if

	t_end = shympi_wtime()
	t_passed = t_end - t_start
	call shympi_time_accum(1,t_passed)

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess, not the final solution

	use mod_system
	use mod_system_interface

	implicit none

	integer, intent(in) :: n,nlvdi,nlv
	real, intent(in)    :: z(nlvdi,n)

	integer i,k,l
	integer n3max
	integer nn !DEB
	real p((nlv+2)*n)

	i = 0
	n3max = l_matrix%n3max

	do k=1,n
	  i = i + 1
	  p(i) = 0.
	  do l=1,nlv
	    i = i + 1
	    p(i) = z(l,k)
	  end do
	  i = i + 1
	  p(i) = 0.
	end do

	nn = n*(nlv + 2) !DEB
	call spk_solve_system(l_matrix,.true.,n3max,nn,p) !DEB

	end

!******************************************************************

	subroutine system_assemble(ie,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system
	use shympi

	implicit none

	integer ie
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	integer i,j,kk,k
	type(smatrix), pointer :: mm

	integer ipext,ieext

	!if( ie > nel_unique ) return	!only assemble from unique elements

	mm => l_matrix

	if( bsysexpl ) then
          do i=1,3
	    k = kn(i)
	    if( id_node(k) /= my_id ) cycle	!only assemble inner nodes
            mm%raux2d(k) = mm%raux2d(k) + mass(i,i)	!GGUEXPL
            mm%rvec2d(k) = mm%rvec2d(k) + rhs(i)
	    do j=1,3
	      if( i /= j .and. mass(i,j) /= 0. ) goto 99
	    end do
	  end do
	else
         do i=1,3
	  k = kn(i)
	  if( id_node(k) /= my_id ) cycle	!only assemble inner nodes
          do j=1,3
            kk=mm%ijp_ie(i,j,ie)			!COOGGU
            if(kk.gt.0) mm%c2coo(kk) = mm%c2coo(kk) + mass(i,j)
          end do
          mm%rvec2d(kn(i)) = mm%rvec2d(kn(i)) + rhs(i)
         end do
	end if

	return
   99	continue
	write(6,*) 'error assembling matrix...'
	write(6,*) ie,ieext(ie)
	write(6,*) kn(i),ipext(kn(i))
	write(6,*) i,j
	write(6,*) kn
	write(6,*) rhs
	write(6,*) mass
	stop 'error stop system_assemble: non diag elems /= 0'
	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix (3d version)

	use mod_system
	use shympi

	implicit none

	integer ie,l,nlv
	integer kn(3)
	real mass(-1:1,3,3)
	real rhs(3)

	integer i,j,kk
	type(smatrix), pointer :: mm

	integer loccoo3d
	external loccoo3d

	if( ie > nel_unique ) return	!only assemble from unique elements

	mm => l_matrix

        do i=1,3
          do j=1,3
	    kk = loccoo3d(i,j,kn,l,ie)
            if(kk.gt.0) then
	       mm%c3coo(kk-1) = mm%c3coo(kk-1) + mass(-1,i,j)
	       mm%c3coo(kk) = mm%c3coo(kk) + mass(0,i,j)
	       mm%c3coo(kk+1) = mm%c3coo(kk+1) + mass(+1,i,j)
	    end if
          end do
	  kk = (nlv+2)*(kn(i)-1) + l + 1
          mm%rvec3d(kk) = mm%rvec3d(kk) + rhs(i) !DEB
        end do
	
	end

!******************************************************************

        subroutine system_get(n,z)

! copies solution back to z

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

	z = real(l_matrix%rvec2d)

        end

!******************************************************************

        subroutine system_get_3d(n,nlvdi,nlv,z)

! copies solution back to z (was system_adjust_3d)

	use mod_system

        implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l

	i = 0

	do k=1,n
	  i = i + 1
	  do l=1,nlv
	    i = i + 1
	    z(l,k) = real(l_matrix%rvec3d(i)) !DEB
	  end do
	  i = i + 1
	end do

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine system_add_rhs(dt,n,array)

! adds right hand side to system array

	use mod_system

        implicit none

        real dt
	integer n
        real array(n)

        integer k

        do k=1,n
          l_matrix%rvec2d(k) = l_matrix%rvec2d(k) + dt * array(k)
        end do

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

	use mod_system
	use mod_system_interface

        implicit none

        call coo_adjust_3d(l_matrix)

        end

!******************************************************************

