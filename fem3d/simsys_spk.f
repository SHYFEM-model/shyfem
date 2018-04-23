!
! revision log :
!
! 12.01.2009	ggu	new file for system routines
! 31.03.2009	ggu	call renamed to spk_*
! 25.05.2015	ggu	some calls changed (pass array in)
! 09.12.2015	ggu	adapted to new pointers and 3d matrix
! 15.12.2015	ggu&deb	finsihed and validated
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
	integer, save, allocatable :: ip_int_nodes(:,:)
	double precision, save, allocatable :: c2coos(:,:)
	real, save, allocatable :: zz(:,:)

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

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines'
        write(6,*) '----------------------------------------'

        call mod_system_init(nkn,nel,ngr,mbw,nlv,l_matrix)
	call mod_system_insert_elem_index(nel,nen3v,l_matrix)
	call mod_system_set_local
	call spk_initialize_system		!calls coo_init_new

! next we first have to set bstsexpl - not yet done... !FIXME

	if( bmpi ) then	!only needed if not explicit !FIXME
	 if( .not. bsysexpl ) then	!only needed if not explicit !FIXME
          call mod_system_init(nkn_global,nel_global,ngr_global
     +				,mbw,nlv,g_matrix)
	  call mod_system_insert_elem_index(nel_global,nen3v_global
     +					,g_matrix)
	  call mod_system_set_global
	  call spk_initialize_system		!calls coo_init_new
	  call system_setup_global_z
	 end if
	end if

        end

!******************************************************************

	subroutine system_init

! on first call initializes pointers - sets arrays to zero
!
! must be called before every assembly of matrix

	use mod_system
	use shympi

	implicit none

	if( bsysexpl ) then
	  call mod_system_set_local
	  a_matrix%rvec2d = 0.
	  a_matrix%raux2d = 0.
	else
	  if( bmpi ) then
	    call mod_system_set_global
	    call spk_init_system
	  end if
	  call mod_system_set_local
	  call spk_init_system
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

	subroutine system_setup_global_z

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
	allocate(zz(nkn_max,0:n_threads-1))
	allocate(ip_int_nodes(nkn_max,0:n_threads-1))
	allocate(i2aux(n2zero_max))
	allocate(j2aux(n2zero_max))
	allocate(ipaux(nkn_max))

	i2aux = 0
	j2aux = 0

!---------------------------------------------------------------
! set up global pointer
!---------------------------------------------------------------

	i2aux(1:n2zero) = mm%i2coo(1:n2zero)
	j2aux(1:n2zero) = mm%j2coo(1:n2zero)

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

	ipaux(1:nkn) = ip_int_node(1:nkn)
	call shympi_gather(ipaux,ip_int_nodes)

!---------------------------------------------------------------
! translate to pointers into global matrix
!---------------------------------------------------------------

	mm => g_matrix
	!mm => l_matrix
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

	subroutine system_solve_global_z(n,z)

! this solves the global system

	use mod_system
	use mod_system_global
	use shympi

	implicit none

	integer n
	real z(n)

	integer n2max,n2zero,n2zero_global
	integer i,id,ipp,nkn,iint,k,nc
	double precision, allocatable :: c2aux(:)
	real, allocatable :: zglobal(:)
	real, allocatable :: zlocal(:)
	type(smatrix), pointer :: mm
	real, parameter :: flag = 1.234567e+20

	mm => l_matrix
	n2zero = n2zero_max

	allocate(c2aux(n2zero_max))
	allocate(zglobal(nkn_global))
	allocate(zlocal(nkn_max))
	zglobal = flag

	!if( shympi_is_master() ) then
	!  write(6,*) 'starting system_solve_global_z'
	!end if

!---------------------------------------------------------------
! gather values from all processes
!---------------------------------------------------------------

	c2aux(1:n2zero) = mm%c2coo(1:n2zero)		!locally assembled
	zlocal = 0.
	zlocal(1:n) = z(1:n)

	call shympi_gather(c2aux,c2coos)		!all contributions
	call shympi_gather(zlocal,zz)			!all zeta values

!---------------------------------------------------------------
! switch to global matrix and assemble matrix and zeta
!---------------------------------------------------------------

	mm => g_matrix
	n2zero_global = mm%n2zero

	do id=0,n_threads-1
	  n2zero = n2_zeros(id)
	  do i=1,n2zero
	    ipp = ij2coos(i,id)
	    mm%c2coo(ipp) = mm%c2coo(ipp) + c2coos(i,id)
	  end do
	  nkn = nkns(id)
	  do i=1,nkn
	    iint = ip_int_nodes(i,id)
	    zglobal(iint) = zz(i,id)
	  end do
	end do

!---------------------------------------------------------------
! solve global matrix
!---------------------------------------------------------------

	n2max = mm%n2max				!maybe also n2zero
	nkn = mm%nkn_system
	if( nkn /= nkn_global ) then
	  stop 'error stop system_solve_global_z: internal error (1)'
	end if
	if( any( zglobal == flag ) ) then
	  nc = count( zglobal == flag )
	  write(6,*) n2max,nkn,nc
	  stop 'error stop system_solve_global_z: internal error (2)'
	end if

	a_matrix => g_matrix
	call spk_solve_system(.false.,n2max,nkn,zglobal)
	a_matrix => l_matrix

	write(6,*) minval(zglobal),maxval(zglobal)
	write(6,*) minval(g_matrix%rvec2d)
     +			,maxval(g_matrix%rvec2d)

!---------------------------------------------------------------
! copy solution for zeta from global to local data structure
!---------------------------------------------------------------

	do k=1,nkn_local
	  iint = ip_int_node(k)
	  l_matrix%rvec2d(k) = g_matrix%rvec2d(iint)
	end do

	stop

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	!if( shympi_is_master() ) then
	!  write(6,*) 'finished system_solve_global_z'
	!end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_solve_z(n,z)

! solves system - z is used for initial guess

	use mod_system
	use shympi

	implicit none

	integer n
	real z(n)

	integer n2max
	type(smatrix), pointer :: mm

	mm => a_matrix

	if( bsysexpl ) then
	  !write(6,*) 'solving explicitly...'
          call shympi_exchange_and_sum_2d_nodes(mm%rvec2d)
          call shympi_exchange_and_sum_2d_nodes(mm%raux2d)
          !call shympi_comment('shympi_elem: exchange rvec2d, raux2d')
	  mm%rvec2d = mm%rvec2d / mm%raux2d	!GGUEXPL
	else if( bmpi ) then
	  call system_solve_global_z(n,z)
	else			!solve directly locally
	  n2max = mm%n2max
	  call spk_solve_system(.false.,n2max,n,z)
	end if

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess

	use mod_system

	implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l
	integer n3max
	real p((nlv+2)*n)
	integer nn !DEB

	i = 0
	n3max = a_matrix%n3max

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
	!call spk_solve_system(.true.,n3max,n,p)
	call spk_solve_system(.true.,n3max,nn,p) !DEB

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

	integer i,j,kk
	type(smatrix), pointer :: m

	if( ie > nel_unique ) return	!only assemble from unique elements

	m => a_matrix

	if( bsysexpl ) then
          do i=1,3
            m%raux2d(kn(i)) = m%raux2d(kn(i)) + mass(i,i)	!GGUEXPL
            m%rvec2d(kn(i)) = m%rvec2d(kn(i)) + rhs(i)
	    do j=1,3
	      if( i /= j .and. mass(i,j) /= 0. ) then
	        write(6,*) ie,kn(i),i,j,mass(i,j)
		stop 'error stop system_assemble: non diag elems /= 0'
	      end if
	    end do
	  end do
	else
         do i=1,3
          do j=1,3
            kk=m%ijp_ie(i,j,ie)			!COOGGU
            if(kk.gt.0) m%c2coo(kk) = m%c2coo(kk) + mass(i,j)
          end do
          m%rvec2d(kn(i)) = m%rvec2d(kn(i)) + rhs(i)
         end do
	end if

	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system

	implicit none

	integer ie,l,nlv
	integer kn(3)
	real mass(-1:1,3,3)
	real rhs(3)

	integer i,j,kk
	type(smatrix), pointer :: m

	integer loccoo3d
	external loccoo3d

	m => a_matrix

        do i=1,3
          do j=1,3
	    kk = loccoo3d(i,j,kn,l,ie)
            if(kk.gt.0) then
	       m%c3coo(kk-1) = m%c3coo(kk-1) + mass(-1,i,j)
	       m%c3coo(kk) = m%c3coo(kk) + mass(0,i,j)
	       m%c3coo(kk+1) = m%c3coo(kk+1) + mass(+1,i,j)
	    end if
          end do
	  kk = (nlv+2)*(kn(i)-1) + l + 1
          m%rvec3d(kk) = m%rvec3d(kk) + rhs(i) !DEB
        end do
	
	end

!******************************************************************

        subroutine system_adjust_z(n,z)

! copies solution back to z

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

	z = real(a_matrix%rvec2d)

        end

!******************************************************************

        subroutine system_adjust_3d(n,nlvdi,nlv,z)

! copies solution back to z

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
	    z(l,k) = real(a_matrix%rvec3d(i)) !DEB
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
          a_matrix%rvec2d(k) = a_matrix%rvec2d(k) + dt * array(k)
        end do

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

        implicit none

        call coo_adjust

        end

!******************************************************************

