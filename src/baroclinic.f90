! dd: newbcl.f,v 1.37 2010-03-08 17:46:45 georg Exp $
!
! baroclinic routines
!
! contents :
!
! subroutine barocl(mode)		amministrates the baroclinic time step
! subroutine rhoset_shell		sets rho iterating to real solution
! subroutine rhoset(resid)		computes rhov and bpresv
! subroutine convectivecorr             convective adjustment
!
! revision log :
!
! revised 30.08.95	$$AUST - austausch coefficient introduced
! revised 11.10.95	$$BCLBND - boundary condition for barocliic runs
! 19.08.1998    ggu     call to barcfi changed
! 20.08.1998    ggu     can initialize S/T from file
! 24.08.1998    ggu     levdbg used for debug
! 26.08.1998    ggu     init, bnd and file routines substituted with con..
! 30.01.2001    ggu     eliminated compile directives
! 05.12.2001    ggu     horizontal diffusion variable, limit diffusion coef.
! 05.12.2001    ggu     compute istot, more debug info
! 11.10.2002    ggu     diffset introduced, shpar = thpar
! 10.08.2003    ggu     qfluxr eliminated (now in subn11.f)
! 10.08.2003    ggu     rhov and bpresv are initialized here
! 04.03.2004    ggu     in init for T/S pass number of vars (inicfil)
! 15.03.2004    ggu     general clean-up, bclint() deleted, new scal3sh
! 17.01.2005    ggu     new difhv implemented
! 15.03.2005    ggu     new diagnostic routines implemented (diagnostic)
! 15.03.2005    ggu     new 3d boundary conditions implemented
! 05.04.2005    ggu     some changes in routine diagnostic
! 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
! 08.06.2007    ggu&deb restructured for new baroclinic version
! 04.10.2007    ggu     bug fix -> call qflux3d with dt double precision
! 17.03.2008    ggu     new open boundary routines introduced
! 08.04.2008    ggu     treatment of boundaries slightly changed
! 22.04.2008    ggu     advection parallelized, no saux1v...
! 23.04.2008    ggu     call to bnds_set_def() changed
! 12.06.2008    ggu     s/tdifhv deleted
! 09.10.2008    ggu     new call to confop
! 12.11.2008    ggu     new initialization, check_layers, initial nos file
! 13.01.2009    ggu&deb changes in reading file in ts_next_record()
! 13.10.2009    ggu     in rhoset bug computing pres
! 13.11.2009    ggu     only initialize T/S if no restart, new rhoset_shell
! 19.01.2010    ggu     different call to has_restart() 
! 16.12.2010    ggu     sigma layers introduced (maybe not finished)
! 26.01.2011    ggu     read in obs for t/s (tobsv,sobsv)
! 28.01.2011    ggu     parameters changed in call to ts_nudge()
! 04.03.2011    ggu     better error message for rhoset_shell
! 31.03.2011    ggu     only write temp/salt if computed
! 04.11.2011    ggu     adapted for hybrid coordinates
! 07.11.2011    ggu     hybrid changed to resemble code in newexpl.f
! 11.11.2011    ggu     restructured ts_next_record() and diagnostic()
! 22.11.2011    ggu     bug fix in ts_file_open() -> bhashl
! 02.12.2011    ggu     adapt ts_file_open() for barotropic version (ihashl)
! 27.01.2012    deb&ggu changes for hybrid in ts_file_open,ts_next_record
! 10.02.2012    ggu     bug in call to ts_next_record (called with nlvddi)
! 23.02.2012    ccf     do noy check depth structure
! 09.03.2012    deb     bug fix in ts_next_record: ilhkv was real
! 31.10.2012    ggu     open and next_record transfered to subtsuvfile.f
! 05.09.2013    ggu     limit salinity to [0,...]
! 25.03.2014    ggu     new offline
! 10.07.2014    ggu     only new file format allowed
! 20.10.2014    ggu     pass ids to scal_adv()
! 10.02.2015    ggu     call to bnds_read_new() introduced
! 15.10.2015    ggu     added new calls for shy file format
! 26.10.2015    ggu     bug fix for parallel code (what was not set)
!
!*****************************************************************
!--------------------------------------------------------------------
        module baroclinic
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

	subroutine barocl(mode)

! amministrates the baroclinic time step
!
! mode : =0 initialize  >0 normal call
!
! written 09.01.94 by ggu  (from scratch)
!
        use layer_thickness
        use ts
        use diffusion
        use hydro_print
        use hydro_admin
        use levels
        use basin
        use shympi
        use mpi_io_admin
        use utility
        use para
        use output
        use conz_util
        use bnd_admin
        use outputd
        use heat_admin2
        use restart
        use shy_util
        use openmp_admin
        use check
        use defnames
        use offline_data,       only: is_offline
        use eq_state
        use bnd_scalar
        use concentration
        use timing

        implicit none
!
! parameter
        include 'param.h'
! arguments
        integer mode
! common
        include 'femtime.h'
        include 'pkonst.h'
        include 'mkonst.h'

! local
        logical debug
        logical badvect
        logical bobs
        logical bgdebug
        logical binfo
        logical bstop
        logical binitial_nos
        logical boff
        integer levdbg
        integer ie
        integer ibarcl
        integer eos_type
        integer idtext,itmext
        integer imin,imax
        integer nintp,nvar
        integer nbc
        integer id
        double precision cdef(1)
        double precision xmin,xmax
        integer itemp,isalt
        double precision salref,temref,sstrat,tstrat
        double precision shpar,thpar
        double precision difmol
        double precision s
        double precision dt
        double precision gamma,gammax
        double precision mass
        double precision wsink
        double precision robs
        double precision dtime0,dtime
        integer isact,l,k,lmax
        integer kspec
        integer icrst
        integer ishyff
        double precision stot,ttot,smin,smax,tmin,tmax,rmin,rmax
        double precision v1,v2,mm
        character*80 file
        character*4 what
! functions
        double precision scalcont,dq

        integer tid
        !integer openmp_get_thread_num
        
        double precision theatold,theatnew
        double precision theatconv1,theatconv2,theatqfl1,theatqfl2
! save
        integer, save :: ia_out(4)
        double precision, save :: da_out(4)

        integer, save, allocatable :: idtemp(:),idsalt(:)
        integer ninfo
        save ninfo

        save badvect,bobs
        save salref,temref
        save difmol
        save itemp,isalt
        save ibarcl
        save eos_type
! data
        integer icall
        save icall
        data icall /0/
        logical bseamount,blockexchange,bintwave        !ivb
        double precision dtseamount, burger_number                  !ivb
        integer iseamount,ilock,iintwave                !ivb

        integer imellor
        logical bmellor
        double precision time1, time2

!----------------------------------------------------------
! parameter setup and check
!----------------------------------------------------------

        if(icall.eq.-1) return

        call is_offline(2,boff)
        !if( boff ) write(6,*) 'TS reading from offline...'
        if( boff ) return

        levdbg = nint(getpar('levdbg'))
        debug = levdbg .ge. 3
        binfo = .true.
        bgdebug = .false.
        binitial_nos = .true.

        ilock           = nint(getpar('ilockex'))       !ivb      
        blockexchange   = ilock .ne. 0                  !ivb

        iintwave        = nint(getpar('iintwave'))      !ivb
        bintwave        = iintwave .ne. 0

        dtime = t_act
        ishyff = nint(getpar('ishyff'))

        imellor = nint(getpar('imellor'))       !ivb
        bmellor = imellor .ne. 0                !ivb

!----------------------------------------------------------
! initialization
!----------------------------------------------------------

        if(icall.eq.0) then     !first time

            ibarcl=iround(getpar('ibarcl'))
            if(ibarcl.le.0) icall = -1
            if(ibarcl.gt.4) then
              write(6,*) 'Value of ibarcl not allowed: ',ibarcl
              stop 'error stop barocl: ibarcl'
            end if
            if(icall.eq.-1) return
            !eos_type = 1 ! Unesco 1980
            !eos_type = 2 ! Jackett & McDougall 1995 
            !eos_type = 3 ! Linearized EOS (T) !ivb 
            eos_type = nint(getpar('eostype'))              ! Type of Equation of State 
            if( eos_type .ge. 4 ) then
              write(6,*) 'eos_type = ',eos_type
              write(6,*) 'unrecognized Equation of State'
              stop 'error stop eos_type: value for eos_type not allowed'
            endif

            badvect = ibarcl .ne. 2
            bobs = ibarcl .eq. 4

            salref=getpar('salref')
            temref=getpar('temref')
            sstrat=getpar('sstrat')
            tstrat=getpar('tstrat')
            difmol=getpar('difmol')
            itemp=iround(getpar('itemp'))
            isalt=iround(getpar('isalt'))

!		--------------------------------------------
!		initialize saltv,tempv
!		--------------------------------------------

            if( .not. rst_use_restart(3) ) then !no restart of T/S values
              call conini(nlvdi,saltv,salref,sstrat,hdkov)
              call conini(nlvdi,tempv,temref,tstrat,hdkov)

              if( ibarcl .eq. 1 .or. ibarcl .eq. 3) then
                call ts_init(itanf,nlvdi,nlv,nkn,tempv,saltv)
              else if( ibarcl .eq. 2 ) then
                call ts_diag(itanf,nlvdi,nlv,nkn,tempv,saltv)
              else if( ibarcl .eq. 4 ) then          !interpolate to T/S
                call ts_nudge(itanf,nlvdi,nlv,nkn,tempv,saltv)
              else
                stop 'error stop barocl: internal error (1)'
              end if
            end if


!           !!! Initialization of temperature for !!!
!           !!! LOCK EXCHANGE TEST CASE!!!
            if (blockexchange .and. itanf.eq.0) then     !ivb
              write(6,*) ' LOCK EXCHANGE temp initialization'
              call temp_ini_lockex(nlvdi,tempv,temref,hdkov,hldv)
            end if

!           --------------------------------------------
!	    initialize observations and relaxation times
!           --------------------------------------------

            tobsv = 0.d0
            sobsv = 0.d0
            rtauv = 0.d0

!	    --------------------------------------------
!	    initialize open boundary conditions
!	    --------------------------------------------

            nbc = nbnds()
            allocate(idtemp(nbc))
            allocate(idsalt(nbc))
            idtemp = 0
            idsalt = 0

            dtime0 = itanf
            nintp = 2
            nvar = 1
            cdef(1) = 0.
		what = 'temp'
                if(bmpi) then
		  call bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlvmax,cdef,idtemp)
		  what = 'salt'
		  call bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlvmax,cdef,idsalt)
                else
		  call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,cdef,idtemp)
		  what = 'salt'
		  call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,cdef,idsalt)
                end if

!	    --------------------------------------------
!	    initialize rhov, bpresv (we call it twice since
!	    --------------------------------------------

!	    rhov depends on bpresv and viceversa
!	    -> we iterate to the real solution)


            !!!! ivb - DON'T DO IN CASE OF RESTART
            if (.not. rst_use_restart(3)) then  
                rhov = 0.d0               !rhov is rho^prime -> 0
                bpresv = 0.d0
                call rhoset_shell
            end if

!	    --------------------------------------------
!           initialize output files
!	    --------------------------------------------

            nvar = 0
            if( itemp .gt. 0 ) nvar = nvar + 1
            if( isalt .gt. 0 ) nvar = nvar + 1
            call init_output('itmcon','idtcon',ia_out)
            if( ishyff == 1 ) ia_out = 0
            if( has_output(ia_out) ) then
              call open_scalar_file(ia_out,nlv,nvar,'nos')
              if( next_output(ia_out) ) then
                if( isalt .gt. 0 ) then
                  !if(bmpi) then
                  !  call rebuild_scalar(nlvdi,saltv)
                  !  if(shympi_partition_on_elements()) then
                  !    if(shympi_is_master()) then 
                  !     call write_scalar_file(ia_out,11,nlvdi,outTempv)
                  !    end if
                  !  end if
                  !else
                    call write_scalar_file(ia_out,11,nlvdi,saltv)
                  !end if
                end if
                if( itemp .gt. 0 ) then
                  !if(bmpi) then
                  !  call rebuild_scalar(nlvdi,tempv)
                  !  if(shympi_partition_on_elements()) then
                  !    if(shympi_is_master()) then
                  !    call write_scalar_file(ia_out,12,nlvdi,outTempv)
                  !    end if
                  !  end if
                  !else
                    call write_scalar_file(ia_out,12,nlvdi,tempv)
                  !end if
                end if
              end if
            end if

            call init_output_d('itmcon','idtcon',da_out)
            if( ishyff == 0 ) da_out = 0
            if( has_output_d(da_out) ) then
              call shy_make_output_name('.ts.shy',file)
              call shy_open_output_file(file,1,nlv,nvar,2,id)
              da_out(4) = id
              if( next_output_d(da_out) ) then
                if( isalt .gt. 0 ) then
                  call shy_write_scalar_record(id,dtime,11,nlvdi,saltv)
                end if
                if( itemp .gt. 0 ) then
                  call shy_write_scalar_record(id,dtime,12,nlvdi,tempv)
                end if
              end if
            end if

            call getinfo(ninfo)

        end if

        icall=icall+1

        if(mode.eq.0) return

!----------------------------------------------------------
! normal call
!----------------------------------------------------------

        wsink = 0.
        robs = 0.
        if( bobs ) robs = 1.

        shpar=getpar('shpar')   !out of initialization because changed
        thpar=getpar('thpar')

        if( ibarcl .eq. 2 ) then
          call ts_diag(it,nlvdi,nlv,nkn,tempv,saltv)
        else if( ibarcl .eq. 4 ) then
          call ts_nudge(it,nlvdi,nlv,nkn,tobsv,sobsv)
        end if

!----------------------------------------------------------
! salt and temperature transport and diffusion
!----------------------------------------------------------

        if( badvect ) then

	  call openmp_get_thread_num(tid)
	  !write(6,*) 'number of thread of temp: ',tid

          if(ln_timing) time1 = shympi_wtime()
          if( itemp .gt. 0 ) then
		what = 'temp'
	        call bnds_read_new(what,idtemp,dtime)
	  end if
          if( isalt .gt. 0 ) then
		what = 'salt'
	        call bnds_read_new(what,idsalt,dtime)
	  end if
          if(ln_timing) then
            time2 = shympi_wtime() - time1
            !io_time = io_time + time2 
            io_scalar_time = io_scalar_time + time2 
          end if
	  
!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(thpar,wsink,robs,itemp,it) 
!$OMP&     SHARED(idtemp,tempv,difhv,difv,difmol,tobsv) DEFAULT(NONE)
!$OMP&     IF(itemp > 0)

          if( itemp .gt. 0 ) then
            what = 'temp'
            call scal_adv_nudge(what,0,tempv,idtemp,thpar,wsink,difhv,difv,difmol,tobsv,robs)
          end if

!$OMP END TASK

!	  call openmp_get_thread_num(tid)
!	  !write(6,*) 'number of thread of salt: ',tid

!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(shpar,wsink,robs,isalt,it) 
!$OMP&     SHARED(idsalt,saltv,difhv,difv,difmol,sobsv) DEFAULT(NONE)
!$OMP&     IF(isalt > 0)

          if( isalt .gt. 0 ) then
	    what = 'salt'
            call scal_adv_nudge(what,0,saltv,idsalt,shpar,wsink,difhv,difv,difmol,sobsv,robs)
          end if

!$OMP END TASK
!$OMP TASKWAIT

	end if

	if( binfo ) then
          if( itemp .gt. 0 ) then
  	    call tsmass(tempv,+1,nlvdi,ttot)
            !tmin=shympi_min(tempv)
            !tmax=shympi_max(tempv) 
      	    call conmima(nlvdi,tempv,tmin,tmax)
!$OMP CRITICAL
            if(shympi_is_master()) then
  	      write(ninfo,*) 'temp: ',it,ttot,tmin,tmax
            end if
!$OMP END CRITICAL
	  end if
          if( isalt .gt. 0 ) then
  	    call tsmass(saltv,+1,nlvdi,stot) 
            !smin=shympi_min(saltv)
            !smax=shympi_max(saltv) 
       	    call conmima(nlvdi,saltv,smin,smax)
!$OMP CRITICAL
            if(shympi_is_master()) then
  	      write(ninfo,*) 'salt: ',it,stot,smin,smax
            end if
!$OMP END CRITICAL
	  end if
	end if

!----------------------------------------------------------
! heat flux through surface	!ccf --> moved to meteo_force
!----------------------------------------------------------

        call compute_heat_flux

!----------------------------------------------------------
! compute rhov and bpresv
!----------------------------------------------------------

        call rhoset_shell

!----------------------------------------------------------
! compute min/max
!----------------------------------------------------------

	!call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	!call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	!call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

!----------------------------------------------------------
! write results to file
!----------------------------------------------------------

	if( next_output(ia_out) ) then
	  if( isalt .gt. 0 ) then
	      call write_scalar_file(ia_out,11,nlvdi,saltv)
            !end if
	  end if
	  if( itemp .gt. 0 ) then
	      call write_scalar_file(ia_out,12,nlvdi,tempv)
	    !end if
	  end if
	end if

	if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( isalt .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,11,nlvdi,saltv)
	  end if
	  if( itemp .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,12,nlvdi,tempv)
	  end if
	end if

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!********************************************************
!********************************************************
!*******************************************************************	

	subroutine ts_diag(it,nlvddi,nlv,nkn,tempv,saltv)

        use para
        use tsfile_admin

	implicit none

	include 'param.h'

	integer it
	integer nlvddi
	integer nlv
	integer nkn
	double precision tempv(nlvddi,nkn)
	double precision saltv(nlvddi,nkn)

	character*80 tempf,saltf
	integer iutemp(3),iusalt(3)
	save iutemp,iusalt

	integer icall
	data icall /0/
	save icall

	tempf = 'temp_diag.dat'
	saltf = 'salt_diag.dat'

	if( icall .eq. 0 ) then
	  call ts_file_open(tempf,it,nkn,nlv,iutemp)
	  call ts_file_open(saltf,it,nkn,nlv,iusalt)
	  call ts_file_descrp(iutemp,'temp diag')
	  call ts_file_descrp(iusalt,'salt diag')
	  icall = 1
	end if

        call ts_next_record(it,iutemp,nlvddi,nkn,nlv,tempv)
        call ts_next_record(it,iusalt,nlvddi,nkn,nlv,saltv)

	end

!*******************************************************************	

	subroutine ts_nudge(it,nlvddi,nlv,nkn,tobsv,sobsv)

        use para
        use tsfile_admin

	implicit none

	include 'param.h'

	integer it
	integer nlvddi
	integer nlv
	integer nkn
	double precision tobsv(nlvddi,nkn)
	double precision sobsv(nlvddi,nkn)

	character*80 tempf,saltf
	integer iutemp(3),iusalt(3)
	save iutemp,iusalt

	integer icall
	data icall /0/
	save icall

	tempf = 'temp_obs.dat'
	saltf = 'salt_obs.dat'

	if( icall .eq. 0 ) then
	  call ts_file_open(tempf,it,nkn,nlv,iutemp)
	  call ts_file_open(saltf,it,nkn,nlv,iusalt)
	  call ts_file_descrp(iutemp,'temp nudge')
	  call ts_file_descrp(iusalt,'salt nudge')
	  icall = 1
	end if

        call ts_next_record(it,iutemp,nlvddi,nkn,nlv,tobsv)
        call ts_next_record(it,iusalt,nlvddi,nkn,nlv,sobsv)

	end

!*******************************************************************	

	subroutine ts_intp(it,nlvddi,nlv,nkn,tobsv,sobsv,tempf,saltf)

        use tsfile_admin

	implicit none

	include 'param.h'

	integer it
	integer nlvddi
	integer nlv
	integer nkn
	double precision tobsv(nlvddi,1)
	double precision sobsv(nlvddi,1)
	character*80 tempf,saltf

	integer iutemp(3),iusalt(3)
	save iutemp,iusalt
	integer ittold,itsold,ittnew,itsnew
	save ittold,itsold,ittnew,itsnew

	double precision, save, allocatable :: toldv(:,:)
	double precision, save, allocatable :: soldv(:,:)
	double precision, save, allocatable :: tnewv(:,:)
	double precision, save, allocatable :: snewv(:,:)

	logical bdebug
	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	bdebug = .true.
	bdebug = .false.

!-------------------------------------------------------------
! initialization (open files etc...)
!-------------------------------------------------------------

	if( icall .eq. 0 ) then
	  write(6,*) 'ts_intp: opening files for T/S'
	  call ts_file_open(tempf,it,nkn,nlv,iutemp)
	  call ts_file_open(saltf,it,nkn,nlv,iusalt)

	  allocate(toldv(nlvddi,nkn))
	  allocate(soldv(nlvddi,nkn))
	  allocate(tnewv(nlvddi,nkn))
	  allocate(snewv(nlvddi,nkn))

	  write(6,*) 'ts_intp: initializing T/S'
	  call ts_next_record(ittold,iutemp,nlvddi,nkn,nlv,toldv)
	  call ts_next_record(itsold,iusalt,nlvddi,nkn,nlv,soldv)
	  write(6,*) 'ts_intp: first record read ',ittold,itsold

	  call ts_next_record(ittnew,iutemp,nlvddi,nkn,nlv,tnewv)
	  call ts_next_record(itsnew,iusalt,nlvddi,nkn,nlv,snewv)
	  write(6,*) 'ts_intp: second record read ',ittnew,itsnew

	  if( ittold .ne. itsold ) goto 98
	  if( ittnew .ne. itsnew ) goto 98
	  if( it .lt. ittold ) goto 99

	  icall = 1
	end if

!-------------------------------------------------------------
! read new files if necessary
!-------------------------------------------------------------

	do while( it .gt. ittnew )

	  ittold = ittnew
	  call copy_record(nkn,nlvddi,nlv,toldv,tnewv)
	  itsold = itsnew
	  call copy_record(nkn,nlvddi,nlv,soldv,snewv)

	  call ts_next_record(ittnew,iutemp,nlvddi,nkn,nlv,tnewv)
	  call ts_next_record(itsnew,iusalt,nlvddi,nkn,nlv,snewv)
	  write(6,*) 'ts_intp: new record read ',ittnew,itsnew

	  if( ittnew .ne. itsnew ) goto 98

	end do

!-------------------------------------------------------------
! interpolate to new time step
!-------------------------------------------------------------

	call intp_record(nkn,nlvddi,nlv,ittold,ittnew,it,toldv,tnewv,tobsv)
	call intp_record(nkn,nlvddi,nlv,itsold,itsnew,it,soldv,snewv,sobsv)

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   98	continue
	write(6,*) ittold,itsold,ittnew,itsnew
	stop 'error stop ts_intp: mismatch time of temp/salt records'
   99	continue
	write(6,*) it,ittold
	stop 'error stop ts_intp: no file for start of simulation'
	end

!*******************************************************************	

	subroutine ts_init(it0,nlvddi,nlv,nkn,tempv,saltv)

! initialization of T/S from file

        use shympi
        use intp_fem_file
        use levels, only: nlvdi
        use basin, only: nkndi,nel
        use mpi_io_admin
        use para
        use tsfile_admin

	implicit none

	include 'param.h'

        integer it0
        integer nlvddi
        integer nlv
        integer nkn
        double precision tempv(nlvddi,nkn)
        double precision saltv(nlvddi,nkn)
        integer nvar
        character*80 tempf,saltf

        integer itt,its
        integer iutemp(3),iusalt(3)
        integer date,time

        nvar=1
	call getfnm('tempin',tempf)
	call getfnm('saltin',saltf)

	if( tempf .ne. ' ' ) then
	  itt = it0
          if(shympi_is_master()) then
	    write(6,*) 'ts_init: opening file for T'
          end if
	  if(bmpi) then
	    call ts_file_open_mpi(tempf,itt,nkndi,nkn,nlvdi,iutemp,nvar)
	  else
	    call ts_file_open(tempf,itt,nkndi,nlvdi,iutemp)
	  end if
	  call ts_file_descrp(iutemp,'temp init')

          if(bmpi) then
            call ts_next_record(itt,iutemp,nlvddi,nkn,nlvdi,tempv)
          else
            call ts_next_record(itt,iutemp,nlvddi,nkndi,nlvdi,tempv)
          end if

	  call ts_file_close(iutemp)
          if(shympi_is_master()) then
            write(6,*) 'temperature initialized from file ',tempf
	  end if
	end if

	if( saltf .ne. ' ' ) then
	  its = it0
          if(shympi_is_master()) then
	    write(6,*) 'ts_init: opening file for S'
          end if
	  if(bmpi) then
	    call ts_file_open_mpi(saltf,its,nkndi,nkn,nlvdi,iusalt,nvar)
	  else
	    call ts_file_open(saltf,its,nkndi,nlvdi,iusalt)
	  end if
	  call ts_file_descrp(iusalt,'salt init')

          if(bmpi) then
            call ts_next_record(its,iusalt,nlvddi,nkn,nlvdi,saltv)
          else
            call ts_next_record(its,iusalt,nlvddi,nkndi,nlvdi,saltv)
          end if

	  call ts_file_close(iusalt)
          if(shympi_is_master()) then
            write(6,*) 'salinity initialized from file ',saltf
	  end if
	end if


	end

!*******************************************************************	

	subroutine intp_record(nkn,nlvddi,nlv,itold,itnew,it,voldv,vnewv,vintpv)

! interpolates records to actual time

	implicit none

	integer nkn,nlvddi,nlv
	integer itold,itnew,it
	double precision voldv(nlvddi,1)
	double precision vnewv(nlvddi,1)
	double precision vintpv(nlvddi,1)

	integer k,l
	double precision rt

        rt = (it-itold) / float(itnew-itold)

	do k=1,nkn
	  do l=1,nlv
	    vintpv(l,k) = voldv(l,k) + rt * (vnewv(l,k) - voldv(l,k))
	  end do
	end do

	end

!*******************************************************************	

	subroutine copy_record(nkn,nlvddi,nlv,voldv,vnewv)

! copies new record to old one

	implicit none

	integer nkn,nlvddi,nlv
	double precision voldv(nlvddi,1)
	double precision vnewv(nlvddi,1)

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    voldv(l,k) = vnewv(l,k)
	  end do
	end do

	end

!*******************************************************************	
!*******************************************************************	
!******************************************************************

	subroutine check_layers(what,vals)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	character*(*) what
	double precision vals(nlvdi,nkn)

	integer l,k,lmax
	double precision valmin,valmax

	write(6,*) 'checking layer structure : ',what

            do l=1,nlv
              valmin = +999.
              valmax = -999.
              do k=1,nkn
                lmax = ilhkv(k)
                if( l .le. lmax ) then
                  valmin = min(valmin,vals(l,k))
                  valmax = max(valmax,vals(l,k))
                end if
              end do
              write(6,*) l,valmin,valmax
            end do

	end

!*******************************************************************	

!--------------------------------------------------------------------
        end module baroclinic
!--------------------------------------------------------------------
