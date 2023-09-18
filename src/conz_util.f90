!
! $Id: subcon1.f,v 1.21 2010-03-11 15:36:39 georg Exp $
!
! routines for concentration (utilities) (old newcon1.f)
!
! contents :
!
! subroutine conini0(nlvddi,c,cref)                      sets initial conditions
! subroutine conini(nlvddi,c,cref,cstrat)		sets initial conditions
!
! subroutine conbnd(nlvddi,c,rbc)			boundary condition
! subroutine con3bnd(nlvddi,c,nlvbnd,rbc)		boundary condition (3D)
!
! subroutine confop(iu,itmcon,idtcon,nlv,nvar,type)	opens (NOS) file
! subroutine confil(iu,itmcon,idtcon,ivar,nlvddi,c)	writes NOS file
! subroutine conwrite(iu,type,nvar,ivar,nlvddi,c)        shell for writing file
!
! subroutine conmima(nlvddi,c,cmin,cmax)                 computes min/max
! subroutine conmimas(nlvddi,c,cmin,cmax)                computes scalar min/max
!
! notes :
!
! conbnd and con3bnd are not used and can be deleted
!
! revision log :
!
! 19.08.1998	ggu	call to conzfi changed
! 20.08.1998	ggu	makew removed (routine used is sp256w)
! 24.08.1998    ggu     levdbg used for debug
! 26.08.1998    ggu     subroutine convol, tstvol transferred to newchk
! 26.08.1998    ggu     all subroutines re-written more generally
! 26.01.1999    ggu     can be used also with 2D routines
! 16.11.2001    ggu     subroutine conmima and diffstab
! 05.12.2001    ggu     new routines diffstab,diffstab1,difflimit
! 11.10.2002    ggu     commented diffset
! 09.09.2003    ggu     new routine con3bnd
! 10.03.2004    ggu     new routine conwrite()
! 13.03.2004    ggu     new routines set_c_bound, distribute_vertically
! 13.03.2004    ggu     exec routine con3bnd() only for level BC (LEVELBC)
! 14.03.2004    ggu     new routines open_b_flux
! 05.01.2005    ggu     routine to write 2d nos file into subnosa.f
! 07.01.2005    ggu     routine diffwrite deleted
! 14.01.2005    ggu     new file for diffusion routines (copied to subdif.f)
! 23.03.2006    ggu     changed time step to double precision
! 31.05.2007    ggu     reset BC of flux type to old way (DEBHELP)
! 07.04.2008    ggu     deleted set_c_bound
! 08.04.2008    ggu     cleaned, deleted distribute_vertically, open_b_flux
! 09.10.2008    ggu&ccf call to confop changed -> nlv
! 20.11.2009    ggu	in conwrite only write needed (nlv) layers
! 20.01.2014    ggu	new writing format for nos files in confop, confil
!
!*****************************************************************
!-----------------------------------------------------------------
        module conz_util
!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

	subroutine conini0(nlvddi,c,cref)

! sets initial conditions (no stratification)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)	!variable to initialize
	double precision cref		!reference value
! common
! local
	integer k,l
	double precision depth,hlayer

	do k=1,nkn
	  do l=1,nlvddi
	    c(l,k) = cref
	  end do
	end do

	end

!*****************************************************************

	subroutine conini(nlvddi,c,cref,cstrat,hdko)

! sets initial conditions (with stratification)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)	!variable to initialize
	double precision cref		!reference value
	double precision cstrat		!stratification [conc/km]
	double precision hdko(nlvddi,nkn)	!layer thickness
! common
! local
	integer k,l
	double precision depth,hlayer

	do k=1,nkn
	  depth=0.
	  do l=1,nlvddi
	    hlayer = 0.5 * hdko(l,k)
	    depth = depth + hlayer
	    c(l,k) = cref + cstrat*depth/1000.
	    depth = depth + hlayer
	  end do
	end do

	end

!*************************************************************

	subroutine conbnd(nlvddi,c,rbc)

! implements boundary condition (simplicistic version)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)	!concentration (cconz,salt,temp,...)
	double precision rbc(nkn)		!boundary condition
! common
	include 'mkonst.h'
! local
	integer k,l,lmax
	double precision rb

	do k=1,nkn
	  if( rbc(k) .ne. flag ) then
	    rb = rbc(k)
	    lmax=ilhkv(k)
	    !write(6,*) 'conbnd: ',k,lmax,rb
	    do l=1,lmax
		c(l,k) = rb
	    end do
	  end if
	end do

	end

!*************************************************************

	subroutine con3bnd(nlvddi,c,nlvbnd,rbc)

! implements boundary condition (simplicistic 3D version)

	use bnd_dynamic
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)	!concentration (cconz,salt,temp,...)
	integer nlvbnd		!vertical dimension of boundary conditions
	double precision rbc(nlvbnd,nkn)	!boundary condition
! common
	include 'mkonst.h'
! local
	integer k,l,lmax
	double precision rb
        integer ipext

	if( nlvbnd .ne. 1 .and. nlvbnd .ne. nlvddi ) then
	  write(6,*) 'nlvddi,nlvbnd: ',nlvddi,nlvbnd
	  stop 'error stop con3bnd: impossible nlvbnd'
	end if
	if( nlvbnd .ne. nlvddi ) then
	  write(6,*) 'nlvddi,nlvbnd: ',nlvddi,nlvbnd
	  stop 'error stop con3bnd: only 3D boundary conditions'
	end if

	do k=1,nkn
	 if( rzv(k) .ne. flag ) then    !only level BC  !LEVELBC !DEBHELP
	  if( rbc(1,k) .ne. flag ) then
	    lmax=ilhkv(k)
            !write(94,*) 'con3bnd: ',k,ipext(k),lmax,nlvbnd,rbc(1,k)
	    if( nlvbnd .eq. 1 ) then
	      rb = rbc(1,k)
	      do l=1,lmax
		c(l,k) = rb
	      end do
	    else
	      do l=1,lmax
		c(l,k) = rbc(l,k)
	      end do
	    end if
	  end if
	 end if
	end do

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine confop(iu,itmcon,idtcon,nl,nvar,type)

! opens (NOS) file

! on return iu = -1 means that no file has been opened and is not written

	use depth
	use levels
	use basin, only : nkn,nel,neldi,nkndi,ngr,mbw
        use defnames
        use para
        use output
        use version
        use ionos
        use shympi
        use mpi_io_admin

	implicit none

	integer iu		!unit				       (in/out)
	integer itmcon		!time of first write		       (in/out)
	integer idtcon		!time intervall of writes	       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

	include 'simul.h'
	include 'femtime.h'

	integer nvers
	integer date,time
	integer ierr
	integer itcon
        integer laux
	!character*80 dir,nam,file
	character*80 title,femver

!-----------------------------------------------------
! check idtcon and itmcon and adjust
!-----------------------------------------------------

	call adjust_itmidt(itmcon,idtcon,itcon)

	iu = -1
        if( idtcon .le. 0 ) return

!-----------------------------------------------------
! open file
!-----------------------------------------------------

	iu = ifemop(type,'unformatted','new')
	if( iu .le. 0 ) goto 98

!-----------------------------------------------------
! initialize parameters
!-----------------------------------------------------

	nvers = 5
	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	title = descrp
	call get_shyfem_version(femver)

!-----------------------------------------------------
! write header of file
!-----------------------------------------------------

	call nos_init(iu,nvers)
	call nos_set_title(iu,title)
	call nos_set_date(iu,date,time)
	call nos_set_femver(iu,femver)

        if(bmpi) then
          call rebuild_nos_header
          laux = shympi_max(nlv)
          if(shympi_is_master().and.shympi_partition_on_elements())then
	    call nos_write_header(iu,nkndi,neldi,laux,nvar,ierr)
            if(ierr.gt.0) goto 99
	    call nos_write_header2(iu,outIlhkv,hlv,outHev,ierr)
            if(ierr.gt.0) goto 99
          end if
        else
	  call nos_write_header(iu,nkn,nel,nl,nvar,ierr)
          if(ierr.gt.0) goto 99
	  call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
          if(ierr.gt.0) goto 99
        end if

!-----------------------------------------------------
! write informational message to terminal
!-----------------------------------------------------

        write(6,*) 'confop: ',type,' file opened ',it

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

	return
   98	continue
	write(6,*) 'error opening file with type ',type
	stop 'error stop confop'
   99	continue
	write(6,*) 'error ',ierr,' writing file with type ',type
	stop 'error stop confop'
	end

!*************************************************************

	subroutine confil(iu,itmcon,idtcon,ivar,nlvddi,c)

! writes NOS file

	use levels
        use ionos
        use shympi
        use mpi_io_admin

	implicit none

	integer iu		!unit
	integer itmcon		!time of first write
	integer idtcon		!time intervall of writes
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,*)		!scalar to write

	include 'femtime.h'

	logical binfo
	integer ierr

	binfo = .false.

!-----------------------------------------------------
! check if files has to be written
!-----------------------------------------------------

	if( iu .le. 0 ) return
	if( it .lt. itmcon ) return
	if( mod(it-itmcon,idtcon) .ne. 0 ) return

!-----------------------------------------------------
! write file
!-----------------------------------------------------

!        if(shympi_partition_on_elements()) then
!          call rebuild_scalar(nlvddi,c)
!          if(shympi_is_master()) then
!	    call nos_write_record(iu,it,ivar,nlvddi,outIlhkv,outTempv,ierr)
!	    if(ierr.gt.0) goto 99
!          end if
!        else
	  call nos_write_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
	  if(ierr.gt.0) goto 99
!        end if

!-----------------------------------------------------
! write informational message to terminal
!-----------------------------------------------------

	if( binfo ) then
          write(6,*) 'confil: variable ',ivar,' written at ',it
	end if

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop confil'
	end

!*************************************************************

	subroutine conwrite(iu,type,nvar,ivar,nlvddi,c)

! shell for writing file unconditionally to disk

	use levels, only : nlvdi,nlv

        implicit none

	integer iu		!unit (0 for first call, set on return)
        character*(*) type      !type of file
	integer nvar		!total number of variables
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,*)	!concentration to write

	include 'femtime.h'

        integer itmcon,idtcon,lmax

        itmcon = itanf
	idtcon = idt
	idtcon = 1
	lmax = min(nlvddi,nlv)

        if( iu .eq. 0 ) then
	  call confop(iu,itmcon,idtcon,lmax,nvar,type)
        end if

	call confil(iu,itmcon,idtcon,ivar,nlvddi,c)

        end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine open_scalar_file(ia_out,nl,nvar,type)
	implicit none
	integer ia_out(4)	!time information		       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)
	double precision da_out(4)
	da_out = ia_out
	call open_scalar_file_d(da_out,nl,nvar,type)
	ia_out = nint(da_out)
	end

!*************************************************************

	subroutine open_scalar_file_d(da_out,nl,nvar,type)

! opens (NOS) file

! on return iu = -1 means that no file has been opened and is not written

	use depth
	use levels
	use basin, only : nkn,nel,ngr,mbw,neldi,nkndi
        use shympi
        use mpi_io_admin
        use defnames
        use para
        use version
        use ionos

	implicit none

	double precision da_out(4)	!time information	       (in/out)
	integer nl		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

	include 'femtime.h'

	include 'simul.h'

	integer nvers
	integer date,time
	integer iu,ierr
	character*80 title,femver

        integer laux

!-----------------------------------------------------
! open file
!-----------------------------------------------------

	iu = ifemop(type,'unformatted','new')
	if( iu .le. 0 ) goto 98
	da_out(4) = iu

!-----------------------------------------------------
! initialize parameters
!-----------------------------------------------------

	nvers = 5
	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))
	title = descrp
	call get_shyfem_version(femver)

!-----------------------------------------------------
! write header of file
!-----------------------------------------------------

	call nos_init(iu,nvers)
	call nos_set_title(iu,title)
	call nos_set_date(iu,date,time)
	call nos_set_femver(iu,femver)
        if(bmpi) then
          call rebuild_nos_header
          laux = shympi_max(nlv)
          if(shympi_is_master().and.shympi_partition_on_elements())then
	    call nos_write_header(iu,nkndi,neldi,laux,nvar,ierr)
            if(ierr.gt.0) goto 99
	    call nos_write_header2(iu,outIlhkv,hlv,outHev,ierr)
            if(ierr.gt.0) goto 99
          end if
        else
	  call nos_write_header(iu,nkn,nel,nl,nvar,ierr)
          if(ierr.gt.0) goto 99
	  call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
          if(ierr.gt.0) goto 99
        end if

!-----------------------------------------------------
! write informational message to terminal
!-----------------------------------------------------

        write(6,*) 'open_scalar_file: ',type,' file opened ',it

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

	return
   98	continue
	write(6,*) 'error opening file with type ',type
	stop 'error stop open_scalar_file'
   99	continue
	write(6,*) 'error ',ierr,' writing file with type ',type
	stop 'error stop open_scalar_file'
	end

!*************************************************************

	subroutine write_scalar_file(ia_out,ivar,nlvddi,c)
	implicit none
	integer ia_out(4)	!time information
	integer ivar		!id of variable to be written
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,*)		!scalar to write
	double precision da_out(4)
	da_out = ia_out
	call write_scalar_file_d(da_out,ivar,nlvddi,c)
	ia_out = nint(da_out)
	end

!*************************************************************

	subroutine write_scalar_file_d(da_out,ivar,nlvddi,c)

! writes NOS file
!
! the file must be open, the file will be written unconditionally

	use levels
        use shympi
        use mpi_io_admin
        use ionos

	implicit none

	double precision da_out(4)	!time information
	integer ivar			!id of variable to be written
	integer nlvddi			!vertical dimension of c
	double precision c(nlvddi,*)		!scalar to write

	include 'femtime.h'

	logical binfo
	integer iu,ierr

	binfo = .false.

!-----------------------------------------------------
! check if files has to be written
!-----------------------------------------------------

	iu = nint(da_out(4))
	if( iu .le. 0 ) return

!-----------------------------------------------------
! write file
!-----------------------------------------------------

        if(shympi_partition_on_elements()) then
          call rebuild_scalar(nlvddi,c)
          if(shympi_is_master()) then
	    call nos_write_record(iu,it,ivar,nlvddi,outIlhkv,outTempv,ierr)
	    if(ierr.gt.0) goto 99
          end if

        else
	  call nos_write_record(iu,it,ivar,nlvddi,ilhkv,c,ierr)
	  if(ierr.gt.0) goto 99
        end if

!-----------------------------------------------------
! write informational message to terminal
!-----------------------------------------------------

	if( binfo ) then
          write(6,*) 'write_scalar_file: ivar = ',ivar,' written at ',it
	end if

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop write_scalar_file: error in writing record'
	end

!*************************************************************
!*************************************************************
!*************************************************************

        subroutine conmima(nlvddi,c,cmin,cmax)

! computes min/max for scalar field

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi

	implicit none

! arguments
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        double precision cmin,cmax
! local
	integer k,l,lmax
	double precision cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
	  end do
	end do

        cmin = shympi_min(cmin)
        cmax = shympi_max(cmax)

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

!*************************************************************

        subroutine conmimas(nlvddi,c,cmin,cmax)

! computes min/max for scalar field -> writes some info

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	integer nlvddi		!vertical dimension of c
	double precision c(nlvddi,nkn)		!concentration (cconz,salt,temp,...)
        double precision cmin,cmax
! common
	include 'femtime.h'
! local
	integer k,l,lmax
        integer ntot
	double precision cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

        ntot = 0
	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
            if( cc .le. 0. ) then
                    ntot = ntot + 1
                    write(96,*) it,l,k,cc,ntot
            end if
	  end do
	end do

        if( ntot .gt. 0 ) then
                write(96,*) 'ntot: ',it,ntot
        end if

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

!*************************************************************
        !LOCK EXCHANGE temperature initialization    !ivb
        subroutine temp_ini_lockex(nlvddi,c,cref,hdko,hldv)

! sets initial conditions of T for lock exchange 
! experiment. T = 30 on the right of basin
! and T = 5 on the left

        use basin, only : nkn,nel,ngr,mbw, xgv
        use shympi

        implicit none

        integer nlvddi          !vertical dimension of c
        double precision c(nlvddi,nkn)      !variable to initialize
        double precision cref               !reference value
        double precision hdko(nlvddi,nkn)   !layer thickness
        double precision hldv(nlvddi)       !layer thickness at rest
! local
        integer k,l
        double precision depth,hlayer

        double precision kt
        double precision mask(nlvddi,nkn)
        logical bsigma

        integer nsigma
        double precision hsigma
        double precision xave
        double precision xgvmin, xgvmax     ! global min/max

        if (b_use_mpi) then
            xgvmin = shympi_min(minval(xgv)) 
            xgvmax = shympi_max(maxval(xgv)) 
        else 
            xgvmin = minval(xgv) 
            xgvmax = maxval(xgv)
        end if

        ! this is a global value       
        xave = 0.5*(xgvmin+xgvmax)

        xave = xave * (1+0.0005)

        write(6,*) 'lock exchange xave'

        do k=1,nkn
          if (xgv(k) .gt. xave) then      
            c(:,k)      = 30.
          else
            c(:,k)      = 5.
          end if                
        end do


        end


!*************************************************************
        subroutine ts_init_mellor(nlvddi,c1,c1ref,c2,c2ref,hdko)

        use basin

        implicit none

        integer nlvddi                          !vertical dimension of c
        double precision c1(nlvddi,nkn),c2(nlvddi,nkn)      !variable to initialize
        double precision c1ref,c2ref                        !reference value
        double precision hdko(nlvddi,nkn)                   !layer thickness

        integer k,l
        double precision depth,hlayer

        do k=1,nkn
          depth=0.
          do l=1,nlvddi
            hlayer      = 0.5 * hdko(l,k)
            depth       = depth + hlayer
            c1(l,k)     = c1ref*(0.85 + 0.15*exp(-depth/1000.0))
            c2(l,k)     = c2ref*(1.06 - 0.06*exp(-depth/1000.0))
            depth       = depth + hlayer
          end do
        end do

        end
!-----------------------------------------------------------------
        end module conz_util
!-----------------------------------------------------------------
