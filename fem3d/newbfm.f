
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

c routines for generic concentration
c
c contents :
c
c subroutine conz3sh
c						shell for conz (new version)
c revision log :
c
c 22.02.2016    ggu&erp     new bfm routines created from newconz
c 06.06.2016    ggu         initialization from file changed
c 28.06.2016    ggu         initialize bfmv, new routine bfm_check()
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

	write(6,*) 'bfm mod init: ',nst,nkn,nlv
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

	integer nvar,nbc,nintp,i,id
	integer nmin
	double precision dtime0

	logical has_restart,has_output_d
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

	  bfmv = 0.

          write(6,*) 'bfm initialized: ',ibfm,nkn,nlvdi
        end if

c-------------------------------------------------------------
c initialization of parameters
c-------------------------------------------------------------

	!cref=getpar('conref')
	rkpar=getpar('chpar')
	difmol=getpar('difmol')
	call get_first_dtime(dtime0)

	nvar = ibfm_state
	allocate(bfminit(nvar))
	allocate(bfmbound(nvar))
	bfminit = 0.				!default initial condition
	bfmbound = 0.				!default boundary condition
	!bfmbound = 3.				!default boundary condition

	!call bfm_check('before init')

        if( .not. has_restart(4) ) then	!no restart of conzentrations
	  call bfm_init_file(dtime0,nvar,nlvdi,nlv,nkn,bfminit,bfmv)
	end if

	!call bfm_check('after init')

        call init_output_d('itmcon','idtcon',da_out)

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

	include 'mkonst.h'

	logical bfirst
	integer nvar,i
	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( ibfm < 0 ) return

	call is_time_first(bfirst)
	if( bfirst ) stop 'bfm_compute_multi: internal error'

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	nvar = ibfm_state
	wsink = 0.
	call get_act_dtime(dtime)
	!dt = idt

	!call bfm_check('before advect')

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

	!call bfm_check('after advect')

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

	integer id,nvar,i,idbase,idc
        real cmin,cmax,ctot
	real v1v(nkn)
	double precision dtime

	logical next_output_d

	if( ibfm < 0 ) return

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	idbase = 30
	idbase = 600
	call get_act_dtime(dtime)
	nvar = ibfm_state

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

	subroutine bfm_check(what)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	character*(*) what

	integer iv,nvar,k,l,lmax
	real vmin,vmax,v
	character*80 textgen,text,line

	nvar = ibfm_state
	vmin = -1000.
	vmax = +1000.
	textgen = 'NaN check'

	do iv=1,nvar
	  write(line,'(i3)') iv
	  text = 'bfmv ' // trim(line) // '  ' // trim(what)
	  call check2Dr(nlvdi,nlv,nkn,bfmv,vmin,vmax
     +				,trim(textgen),trim(text))
	  do k=1,nkn
	    lmax = ilhkv(k)
	    do l=1,lmax
	      v = bfmv(l,k,iv)
	      if( v <= 0 ) then
		write(6,*) 'bfm_check: ',iv,k,l,v
	      end if
	    end do
	  end do
	end do

	end subroutine bfm_check

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

