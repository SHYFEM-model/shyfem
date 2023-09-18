
! $Id: subbnds.f,v 1.18 2010-02-26 17:35:06 georg Exp $
!
! handles open boundary conditions for scalar variables
!
! contents :
!
! subroutine bnds_init(text,file,nintp,nvar,ndim,array,aconst)
!			initializes boundary condition
! subroutine bnds_set(text,t,ndim,array,aaux)
!			sets boundary condition
! subroutine bnds_trans(text,ndim,array,aaux,ivar,nlvddi,r3v)
!			transfers boundary condition to matrix
! subroutine bnds_set_def(text,ndim,array)
!			sets default value for boundaries
! subroutine bnds_print(text,ndim,array)
!			prints boundary condition
! subroutine init_scal_bc(r3v)		initializes array for scalar BC
!
! revision log :
!
! 10.08.2003	ggu	exit if no file can be opened
! 03.03.2005	ggu	new call to exfini
! 01.02.2006	ggu	bugfix -> nextra was 3 -> not used anymore
! 03.02.2006	ggu	bugfix -> calling exfpres
! 17.02.2006	ggu	included text for debug, new call to get_bflux()
! 23.03.2006    ggu     changed time step to double precision
! 05.10.2007    ggu     definition of array(ndim,1) changed to array(ndim,0:1)
! 17.03.2008    ggu     routines re-arranged, new bnds_trans, bnds_set_def
! 17.04.2008    ggu     deleted bnds_set_global
! 23.04.2008    ggu     in bnds_set_def() eliminated aaux
! 16.02.2012    ggu     new routine bnds_init0 (to force spatially const bound)
! 25.06.2014    ggu     new routines bnds_init_new() and bnds_trans_new()
! 10.07.2014    ggu     only new file format allowed
! 05.02.2015    ggu     check for number of variables read
! 10.02.2015    ggu     new routine bnds_read_new()
! 30.09.2015    ggu     new routine iff_flag_ok() for ambient value
!
!******************************************************************
!-------------------------------------------------------------------------------
        module bnd_scalar
!-------------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------------

	subroutine bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,cdef,ids)

! initializes boundary condition for scalars

	use intp_fem_file
        use bnd_admin

	implicit none

	character*(*) what	!what variable to initialize
	double precision dtime0	!time
	integer nintp		!degree of interpolation - same for all bounds
	integer nvar		!number of variables
	integer nkn		!number of points (max)
	integer nlv		!number of vertical levels
	double precision cdef(nvar)		!default values
	integer ids(*)		!boundary info id (return)

	character*80 file
	integer nbc,ibc
	integer nk,nsize
	integer iunit,n,i,id
	double precision val

	integer nodes(nkn)
	double precision aconst(nvar)

	integer nvar_orig
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) then
	  write(6,*) '-------------------------'
	  write(6,*) 'Initialization for scalar:'
	  write(6,*) what,nvar
	  write(6,*) '-------------------------'
	end if

	nvar_orig = nvar
	nbc = nbnds()

	do ibc=1,nbc

	  ids(ibc) = 0
          nk = nkbnds(ibc)
	  if( nk .le. 0 ) cycle
          do i=1,nk
            nodes(i) = kbnds(ibc,i)
          end do

	  call get_boundary_file(ibc,what,file)

	  aconst = cdef
	  if( exists_bnd_name(what) ) then
	    call get_bnd_par(ibc,what,val)
	    aconst = val
	  end if

          call iff_init(dtime0,file,nvar,nk,nlv,nintp,nodes,aconst,id)
	  if( nvar /= nvar_orig ) goto 99
	  call iff_set_description(id,ibc,what)
	  call iff_flag_ok(id)		!can deal with ambient value

	  ids(ibc) = id

	  if( bdebug ) then
            write(6,*) 'boundary: ',ibc,id
            if( file .ne. ' ' ) then
	      write(6,*) '  file: ',file(1:60)
	    else
              write(6,*) '  def: ',aconst
	    end if
	  end if

        end do

	if( bdebug ) then
	  do ibc=1,nbc
	    call iff_print_info(ids(ibc))
	  end do
	end if

	return
   99	continue
	write(6,*) 'number of variables is wrong'
	write(6,*) 'expected: ',nvar_orig
	write(6,*) 'read from file: ',nvar
	write(6,*) 'type of boundary condition: ',what
	write(6,*) 'number of boundary: ',ibc
	write(6,*) 'file name: ',trim(file)
	call iff_print_info(id)
	stop 'error stop bnds_init_new: wrong number of variables'
	end

!******************************************************************

	subroutine bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlv,cdef,ids)

! initializes boundary condition for scalars

	use intp_fem_file
        use shympi
        use bnd_admin

	implicit none

	character*(*) what	!what variable to initialize
	double precision dtime0	!time
	integer nintp		!degree of interpolation - same for all bounds
	integer nvar		!number of variables
	integer nkn		!number of points (max)
	integer nlv		!number of vertical levels
	double precision cdef(nvar)		!default values
	integer ids(*)		!boundary info id (return)

	character*80 file
	integer nbc,ibc
	integer nk,nsize
	integer iunit,n,i,id
	double precision val

	integer nodes(bounds%nnob)
	double precision aconst(nvar)

	integer nvar_orig
	logical bdebug


        if(size(bounds%niob) .le. 0) return

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) then
	  write(6,*) '-------------------------'
	  write(6,*) 'Initialization for scalar:'
	  write(6,*) what,nvar
	  write(6,*) '-------------------------'
	end if

	nvar_orig = nvar
	nbc = nbnds()

	do ibc=1,nbc

	  ids(ibc) = 0
          nk = nkbnds(ibc)
	  if( nkbnds(ibc) .le. 0 ) cycle
          do i=1,nk
            nodes(i) = kbnds(ibc,i)
          end do

	  call get_boundary_file(ibc,what,file)

	  aconst = cdef
	  if( exists_bnd_name(what) ) then
	    call get_bnd_par(ibc,what,val)
	    aconst = val
	  end if

          call iff_init_mpi(dtime0,file,nvar,bounds%tnob(ibc),nk,nlv,nintp,nodes,aconst,id,ibc)
	  if( nvar /= nvar_orig ) goto 99
	  call iff_set_description(id,ibc,what)
	  call iff_flag_ok(id)		!can deal with ambient value

	  ids(ibc) = id

	  if( bdebug ) then
            write(6,*) 'boundary: ',ibc,id
            if( file .ne. ' ' ) then
	      write(6,*) '  file: ',file(1:60),my_id
	    else
              write(6,*) '  def: ',aconst
	    end if
	  end if

        end do

	if( bdebug ) then
	  do ibc=1,nbc
	    call iff_print_info(ids(ibc))
	  end do
	end if

	return
   99	continue
	write(6,*) 'number of variables is wrong'
	write(6,*) 'expected: ',nvar_orig
	write(6,*) 'read from file: ',nvar
	write(6,*) 'type of boundary condition: ',what
	write(6,*) 'number of boundary: ',ibc
	write(6,*) 'file name: ',trim(file)
	call iff_print_info(id)
	stop 'error stop bnds_init_mpi: wrong number of variables'
	end

!******************************************************************

	subroutine bnds_read_new(text,ids,dtime)

! reads new boundary condition

	use intp_fem_file
        use bnd_admin

	implicit none

	character*(*), intent(in)	:: text		!for debug
	integer, intent(in)		:: ids(*)
	double precision, intent(in)	:: dtime

	integer nbc,ibc,id

	nbc = nbnds()

	do ibc=1,nbc
	  id = ids(ibc)
	  if( id .le. 0 ) cycle
	  call iff_read_and_interpolate(id,dtime)
	end do

	end

!******************************************************************

	subroutine bnds_trans_new(text,ids,dtime,ivar,nkn,nlv,nlvddi,r3v)

! transfers boundary condition to matrix

	use intp_fem_file
        use bnd_admin
        use basin,      only:   dist_3d

	implicit none

	character*(*) text	!text for debug
	integer ids(*)
	double precision dtime
	integer ivar		!variable to use (can be 0 -> 1)
	integer nkn
	integer nlv
	integer nlvddi		!vertical dimension of levels
	double precision r3v(nlvddi,nkn)	!matrix to which BC values are transfered

	integer nbc,ibc
	integer nvar,nsize,ndata
	integer nk,iv,kn
	integer i,ip,id
	double precision t

	double precision vals(nlv,nkn)

	call init_scal_bc(r3v)	!sets r3v to flag - dims nlv,nkn

	nbc = nbnds()
	iv = max(ivar,1)

	do ibc=1,nbc

	  nk = nkbnds(ibc)   !total number of nodes of this boundary
	  id = ids(ibc)
	  if( id .le. 0 ) cycle

	  call iff_time_interpolate(id,dtime,iv,nkn,nlv,vals)

	  do i=1,nk
	    kn = kbnds(ibc,i)
	    call dist_3d(nlvddi,r3v,kn,nlv,vals(1,i))
	  end do

	end do

	!call iff_print_info(44,0,.true.)
	!call iff_print_boundary_info(9,0,.true.)

	end

!******************************************************************

	subroutine bnds_init(text,file,nintp,nvar,ndim,array,aconst)

! initializes boundary condition

	implicit none

	character*(*) text	!text for debug
	character*80 file(1)	!file names
	integer nintp		!degree of interpolation - same for all bounds
	integer nvar		!number of variables
	integer ndim		!first dimension of array
	double precision array(ndim,0:1)	!array with all information
	double precision aconst(nvar)	!if no file is given constant values from here

	stop 'error stop: call to bnds_init not allowed'

	end

!******************************************************************

	subroutine bnds_set_def(text,ndim,array)

! sets default value for boundaries - works only for nvar = 1

	implicit none

	character*(*) text	!text for debug
	integer ndim		!first dimension of array
	double precision array(ndim,0:1)	!array with all information

	stop 'error stop: call to bnds_set_def not allowed'

	end

!******************************************************************

	subroutine init_scal_bc(r3v)

! initializes array for scalar boundary condition

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision r3v(nlvdi,nkn)

	include 'mkonst.h'

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    r3v(l,k) = flag
	  end do
	end do

	end

!**********************************************************************

!-------------------------------------------------------------------------------
        end module bnd_scalar
!-------------------------------------------------------------------------------
