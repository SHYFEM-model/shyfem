c
c $Id: newconz.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c routines for generic concentration
c
c contents :
c
c subroutine conz3sh
c						shell for conz (new version)
c revision log :
c
c 22.02.2016    ggu&eps     new bfm routines created from newconz
c 06.06.2016    ggu         initialization from file changed
c
c*********************************************************************

!==================================================================
        module mod_bfm
!==================================================================

        implicit none

        integer, private, save :: nst_bfm = 0
        integer, private, save :: nkn_bfm = 0
        integer, private, save :: nlv_bfm = 0

        real, allocatable, save :: bfmv(:,:,:)

	integer, save :: ibfm = 0
	integer, save :: ibfm_state = 50	!number of state variables
        integer, save :: icall_bfm = 0
        integer, save :: ninfo_bfm = 0

        logical, save :: binfo_bfm = .true.

        real, save :: cref,rkpar,difmol

        integer, save, allocatable :: idbfm(:)
        integer, save :: ia_out(4)
        double precision, save :: da_out(4)

        real, save, allocatable :: bfminit(:)
        real, save, allocatable :: bfmbound(:)

	character*3, save :: what = 'bfm'

!==================================================================
	contains
!==================================================================

        subroutine mod_bfm_init(nst,nkn,nlv)

        integer nst
        integer nkn
        integer nlv

        if( nst == nst_bfm .and. nkn == nkn_bfm 
     +		.and. nlv == nlv_bfm ) return

        if( nst > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nst == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nst,nkn,nlv: ',nst,nkn,nlv
            stop 'error stop mod_bfm_init: incompatible parameters'
          end if
        end if

        if( nkn_bfm > 0 ) then
          deallocate(bfmv)
        end if

        nst_bfm = nst
        nkn_bfm = nkn
        nlv_bfm = nlv

        if( nkn == 0 ) return

        allocate(bfmv(nlv,nkn,nst))

        end subroutine mod_bfm_init

c*********************************************************************

	subroutine bfm_init

c initializes bfm computation

	!use mod_bfm
	!use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer nvar,nbc,nintp,i,id
	integer nmin,ishyff
	double precision dtime0

	logical has_restart,has_output,has_output_d
	integer nbnds
	real getpar

c-------------------------------------------------------------
c initialization of module
c-------------------------------------------------------------

	if( ibfm < 0 ) return

        if( ibfm == 0 ) then
          ibfm=nint(getpar('ibfm'))
          if( ibfm <= 0 ) ibfm = -1
          if( ibfm < 0 ) return

          call mod_bfm_init(ibfm_state,nkn,nlvdi)

          write(6,*) 'bfm initialized: ',ibfm,nkn,nlvdi
        end if

c-------------------------------------------------------------
c initialization of parameters
c-------------------------------------------------------------

	!cref=getpar('conref')
	rkpar=getpar('chpar')
	difmol=getpar('difmol')
        ishyff = nint(getpar('ishyff'))
	dtime0 = itanf

	nvar = ibfm_state
	allocate(bfminit(nvar))
	allocate(bfmbound(nvar))
	bfminit = 0.				!default initial condition
	bfmbound = 0.				!default boundary condition

        if( .not. has_restart(4) ) then	!no restart of conzentrations
	  call bfm_init_file(dtime0,nvar,nlvdi,nlv,nkn,bfminit,bfmv)
	end if

        call init_output('itmcon','idtcon',ia_out)
	if( ishyff == 1 ) ia_out = 0

	if( has_output(ia_out) ) then
          call open_scalar_file(ia_out,nlv,nvar,'bfm')
	end if

        call init_output_d('itmcon','idtcon',da_out)
        if( ishyff == 0 ) da_out = 0

        if( has_output_d(da_out) ) then
          call shyfem_init_scalar_file('bfm',nvar,.false.,id)
          da_out(4) = id
        end if

        call getinfo(ninfo_bfm)

        nbc = nbnds()
        allocate(idbfm(nbc))
        idbfm = 0

	nintp = 2
        call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +				,bfmbound,idbfm)

	end subroutine bfm_init

c*********************************************************************

	subroutine bfm_compute

	!use mod_bfm
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'
	include 'mkonst.h'

	integer nvar,i
	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( ibfm < 0 ) return

	if( it .eq. itanf ) stop 'bfm_compute_multi: internal error'
	!if( it .eq. itanf ) return

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	nvar = ibfm_state
	wsink = 0.
	dtime = it
	dt = idt

	call bnds_read_new(what,idbfm,dtime)

	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idbfm,what,
!$OMP&     dt,nlvdi) SHARED(bfmv)  DEFAULT(NONE)
 
! use scal_adv_fact to pass wsink matrix into advection

          call scal_adv(what,i
     +                          ,bfmv(1,1,i),idbfm
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

!$OMP END TASK

	end do	

!$OMP TASKWAIT

	icall_bfm = icall_bfm + 1

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end subroutine bfm_compute

c*********************************************************************

	subroutine bfm_write

	!use mod_bfm
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer id,nvar,i,idbase,idc
        real cmin,cmax,ctot
	real v1v(nkn)
	double precision dtime

	logical next_output,next_output_d

	if( ibfm < 0 ) return

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	idbase = 30
	idbase = 600
	dtime = t_act
	nvar = ibfm_state

	if( next_output(ia_out) ) then
	  do i=1,nvar
	    idc = idbase + i
	    call write_scalar_file(ia_out,idc,nlvdi,bfmv(1,1,i))
	  end do
	end if

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
	  do i=1,nvar
	    idc = idbase + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi,bfmv(1,1,i))
	  end do
        end if

c-------------------------------------------------------------
c write to info file
c-------------------------------------------------------------

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end subroutine bfm_write

c*********************************************************************
c*********************************************************************
c*********************************************************************

c*********************************************************************

        subroutine bfm_init_file(dtime,nvar,nlvddi,nlv,nkn,val0,val)

c initialization of bfm from file

        implicit none

	double precision dtime
        integer nvar
        integer nlvddi
        integer nlv
        integer nkn
	real val0(nvar)
        real val(nlvddi,nkn,nvar)

        call tracer_file_init('bfm init','bfmini',dtime
     +                          ,nvar,nlvddi,nlv,nkn,val0,val)

	end subroutine bfm_init_file

c*********************************************************************

!==================================================================
        end module mod_bfm
!==================================================================

