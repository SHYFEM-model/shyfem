!
! $Id: newfix.f,v 1.5 2009-03-24 17:49:19 georg Exp $
!
! fix or nudge velocities at open boundary
!
! contents :
!
! subroutine bclfix_ini		initialization of bclfix routines
! subroutine bclfix		fixes velocities on open boundaries
!
! revision log :
!
! 15.09.2008    ggu     written from scratch
! 03.11.2008    ggu&dbf nudging implemented
! 12.11.2008    ggu     handle sigma coordinates
! 06.12.2008    ggu     read nbfix from STR
! 19.01.2009    ggu     no error stop in initializing when nbfix=0
! 23.03.2009    ggu     tramp from start of simulation
! 16.12.2010    ggu     bsigma renamed to bosigma
! 29.10.2014    ccf     rewritten for 7_0_3, vel file for each boundary
!*****************************************************************
!-----------------------------------------------------------------
        module bcvel
!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

        subroutine bclfix_ini

! initialization of bclfix routines

	use bclfix_util
	use internal
	use basin
        use bnd_admin
        use fem_util

        implicit none 

	double precision tnudge	!relaxation time for nudging [s]

        integer ie,l,i,ii,k,n,nn,nf
	integer ibc,nodes
	integer nbc
	integer iflag(nkn)
       
	logical bdebug

	bdebug = .true.
	bdebug = .false.

!------------------------------------------------------------------
! initialize arrays
!------------------------------------------------------------------

	iuvfix = 0
	tnudgev = 0.
	ielfix = 0

!------------------------------------------------------------------
! loop over boundaries
!------------------------------------------------------------------

	nbc = nbnds()

	do ibc = 1, nbc

          call get_bnd_par(ibc,'tnudge',tnudge)

	  if (tnudge .lt. 0. ) cycle

	  nodes = nkbnds(ibc)

	  !------------------------------------------------------------------
	  ! flag boundary nodes
	  !------------------------------------------------------------------

	  iflag = 0

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    iflag(k) = 1
	  end do

	  !------------------------------------------------------------------
	  ! set up ielfix and tnudgev/iuvfix
	  !------------------------------------------------------------------

	  do ie=1,nel

	    n = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( iflag(k) .ne. 0 ) then
	        n = n + 1
	        ielfix(n,ie) = k
	      end if
	    end do
	    ielfix(0,ie) = n			!total number of nodes for ele

	    if( n .gt. 0 ) then			!nudging or fixing
  	      if( tnudge == 0. ) then
		iuvfix(ie) = 1
  	      else if( tnudge > 0. ) then
		tnudgev(ie) = tnudge
	      else
		stop 'error stop bclfix_ini: internal error (1)'
	      end if
	    end if

	  end do

!------------------------------------------------------------------
! end loop over boundaries
!------------------------------------------------------------------

	end do

!------------------------------------------------------------------
! debug output
!------------------------------------------------------------------

	if( .not. bdebug ) return

	write(6,*) 'bclfix_ini has been initialized'

	nf = 0
	nn = 0
	do ie=1,nel
	  if( ielfix(0,ie) .ne. 0 ) then
	    if( iuvfix(ie) .ne. 0 ) then
	     write(6,*) 'Fix ie = ',ie,ieext(ie),'nn = ',ielfix(0,ie)
	     nf = nf + 1
	    else
	     write(6,*) 'Nudge ie = ',ie,ieext(ie),'nn = ',ielfix(0,ie)
	     nn = nn + 1
	    end if
	  end if
	end do

	write(6,*) '-------------'
	write(6,*) 'bclfix_ini has found ',nf,' elements to be fixed'
	write(6,*) 'bclfix_ini has found ',nn,' elements to be nudged'
	write(6,*) '-------------'

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

        end

!*****************************************************************

        subroutine bclfix

! fix or nudge  velocities on open boundaries

	use bclfix_util
	use internal
	use geom_dynamic
	use layer_thickness
	use hydro_vel
	use hydro_admin
	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi
        use bnd_admin
        use bnd_scalar
        use timing

        implicit none 

        include 'param.h' 

	include 'femtime.h'

	double precision tnudge	!relaxation time for nudging [s]
	double precision tramp	!time for smooth init

        integer ie,l,i,k,ii,n
	integer lmax
	integer nbc
        double precision u(nlvdi),v(nlvdi)
        double precision h,t
        double precision alpha
	double precision uexpl,vexpl
        
	integer nintp,nvar
	double precision cdef(2)
        double precision dtime0,dtime,time1,time2
        integer, save, allocatable :: idvel(:)
        character*10 what

	integer icall
	save icall
	data icall /0/

        if( icall .eq. -1 ) return

!------------------------------------------------------------------
! initialize open boundary conditions
!------------------------------------------------------------------

	if( icall .eq. 0 ) then

          icall = -1
	  do ie = 1,nel
	    if ( ielfix(0,ie) .ne. 0 ) icall = 1
	  end do

          if( icall .eq. -1 ) return

          nbc = nbnds()
          allocate(idvel(nbc))
	  idvel = 0

          dtime0  = itanf
          nintp   = 2
          nvar    = 2
          cdef    = 0.
          what    = 'vel'
          if(bmpi) then
            call bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlvmax,cdef,idvel)
          else
            call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,cdef,idvel)
          end if

	end if

!------------------------------------------------------------------
! read boudary velocities from file and store in u/vbound
!------------------------------------------------------------------

        dtime = it

        if(ln_timing) time1 = shympi_wtime()
        call bnds_read_new(what,idvel,dtime)
        if(ln_timing) then
           time2 = shympi_wtime() - time1
           !io_time = io_time + time2
           io_transp_time = io_transp_time + time2
        end if
        call bnds_trans_new(what,idvel,dtime,1,nkn,nlvmax,nlvdi,ubound)
        call bnds_trans_new(what,idvel,dtime,2,nkn,nlvmax,nlvdi,vbound)

!------------------------------------------------------------------
! simulate smooth initial discharge with tramp
!------------------------------------------------------------------

        tramp = 43200.          !smooth init
        tramp = 0.              !smooth init

        if( it-itanf .le. tramp ) then
	  alpha = (it-itanf) / tramp
        else
	  alpha = 1.
        endif

!------------------------------------------------------------------
! nudge of fix velocities in elements
!------------------------------------------------------------------

	do ie=1,nel
	  n = ielfix(0,ie)
	  if( n .gt. 0 ) then

	    lmax = ilhv(ie)

	    do l=1,lmax
	      u(l) = 0.
	      v(l) = 0.
	    end do

	    do i=1,n
	      k = ielfix(i,ie)
	      do l=1,lmax
	        u(l) = u(l) + ubound(l,k)
	        v(l) = v(l) + vbound(l,k)
	      end do
	    end do

	    do l=1,lmax
	      u(l) = u(l) / n
	      v(l) = v(l) / n
	    end do

	    tnudge = tnudgev(ie)

	    if (iuvfix(ie) .eq. 1 ) then	!fix velocities
              do l=1,lmax
                h = hdeov(l,ie)
                ulnv(l,ie) = u(l)*alpha
                vlnv(l,ie) = v(l)*alpha
                utlnv(l,ie) = ulnv(l,ie)*h
                vtlnv(l,ie) = vlnv(l,ie)*h
              end do
	    else if( tnudge > 0 ) then		!nudge velocities
	      do l=1,lmax
                h = hdeov(l,ie)
                uexpl = (ulnv(l,ie)-u(l))*alpha*h/tnudge
                vexpl = (vlnv(l,ie)-v(l))*alpha*h/tnudge
	        fxv(l,ie) = fxv(l,ie) + uexpl
	        fyv(l,ie) = fyv(l,ie) + vexpl
	      end do
	    end if
	  end if
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	stop 'error stop bclfix: internal error (1)'
        end

!*****************************************************************

!-----------------------------------------------------------------
        end module bcvel
!-----------------------------------------------------------------
