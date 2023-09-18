!
! $Id: subflxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
!
! subroutines for computing discharge / flux
!
! contents :
!
! subroutine inflxa
! subroutine rdflxa
! subroutine ckflxa
! subroutine prflxa
! subroutine tsflxa
!
! subroutine wrflxa(it)				write of flux data
!
! subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
! subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
!
! subroutine flxini				initializes flux routines
! subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
! subroutine flxinf(m,kflux,iflux)		sets up one info structure
! function igtnsc(k1,k2)			gets number of internal section
!
! revision log :
!
! 30.04.1998    ggu	newly written routines (subpor deleted)
! 07.05.1998    ggu	check nrdveci on return for error
! 08.05.1998    ggu	restructured with new comodity routines
! 13.09.1999    ggu	type of node computed in own routine flxtype
! 19.11.1999    ggu	iskadj into sublin
! 20.01.2000	ggu	old routines substituted, new routine extrsect
! 20.01.2000    ggu     common block /dimdim/ eliminated
! 20.01.2000    ggu     common block /dimdim/ eliminated
! 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
! 26.05.2003	ggu	in flxnov substituted a,b with b,c
! 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
! 10.08.2003	ggu	do not call setweg, setnod, setkan
! 23.03.2006    ggu     changed time step to double precision
! 28.09.2007    ggu     use testbndo to determine boundary node in flxtype
! 28.04.2009    ggu     links re-structured
! 23.02.2011    ggu     new routine call write_node_fluxes() for special output
! 01.06.2011    ggu     documentation to flxscs() changed
! 21.09.2011    ggu     some lower-level subroutines copied to subflx.f
! 07.10.2011    ggu     adjusted for 3d flux routines
! 19.10.2011    ggu     added T/S variables, created fluxes_*() routines
! 19.10.2011    ggu     added conz variables, created fluxes_template()
! 10.05.2013    ggu     introduced subflxa.h, common routines to subflxu.f
! 20.05.2015    ggu     modules introduced
!
! notes :
!
! These routines can also be used internally to compute the flux
! over various sections. The following calling sequence must be respected:
!
! call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
!
! call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
!
! Initialization can be done anytime.
!
!******************************************************************
!******************************************************************
!******************************************************************

!==================================================================
        module flux
!==================================================================

        implicit none

        integer, save :: nsect = -1
        integer, save :: kfluxm = 0
        integer, save, allocatable :: kflux(:)
        integer, save, allocatable :: iflux(:,:)

        integer, save, allocatable :: nlayers(:)
        double precision, save, allocatable :: fluxes(:,:,:)

        double precision, save, allocatable :: masst(:,:,:)
        double precision, save, allocatable :: saltt(:,:,:)
        double precision, save, allocatable :: tempt(:,:,:)
        double precision, save, allocatable :: conzt(:,:,:)

        integer, save, allocatable :: nrcc(:)
        double precision, save, allocatable :: cflux(:,:,:,:)

!==================================================================
        contains
!==================================================================

!==================================================================

        subroutine flux_read_section(n)

        use nls

        integer n

        n = nls_read_vector()
        kfluxm = n

        if( n > 0 ) then
          allocate(kflux(n))
          allocate(iflux(3,n))
          call nls_copy_int_vect(n,kflux)
        end if

        end subroutine flux_read_section

!******************************************************************

        subroutine flux_alloc_arrays(nl,ns)

	implicit none

	integer nl	!layers
	integer ns	!sections

	if( ns <= 0 ) return

        allocate(nlayers(ns))
        allocate(fluxes(0:nl,3,ns))

        allocate(masst(0:nl,3,ns))
        allocate(saltt(0:nl,3,ns))
        allocate(tempt(0:nl,3,ns))
        allocate(conzt(0:nl,3,ns))

	end

!******************************************************************

        subroutine flux_alloc_conz_arrays(nl,ns,nc)

	implicit none

	integer nl	!layers
	integer ns	!sections
	integer nc	!number of concentrations

	if( ns <= 0 .or. nc <= 0 ) return

        allocate(nrcc(nc))
        allocate(cflux(0:nl,3,ns,nc))

	end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine mod_flx(mode)

        implicit none
 
        integer mode
 
        include 'modules.h'
 
	include 'femtime.h'
 
        if( mode .eq. M_AFTER ) then
           call wrflxa(it)
        else if( mode .eq. M_INIT ) then
           call inflxa
        else if( mode .eq. M_READ ) then
           call rdflxa
        else if( mode .eq. M_CHECK ) then
           call ckflxa
        else if( mode .eq. M_SETUP ) then
           call flxini
        else if( mode .eq. M_PRINT ) then
           call prflxa
        else if( mode .eq. M_TEST ) then
           call tsflxa
        else if( mode .eq. M_BEFOR ) then
!          nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_flx'
        end if
 
        end

!******************************************************************

        subroutine inflxa

! nsect		total number of sections
! kfluxm	total number of nodes defining sections
! kflux()	node numbers defining sections

        implicit none

        end

!******************************************************************

        subroutine rdflxa

        implicit none

	integer n

        call flux_read_section(n)

        if( n .lt. 0 ) then
          write(6,*) 'read error in section $flux'
          stop 'error stop rdflxa'
        end if

        end

!******************************************************************

        subroutine ckflxa

        use fem_util

        implicit none

	integer k,ii
        logical berror

        berror=.false.
        if(kfluxm .gt. 0) then
	  call n2int(kfluxm,kflux,berror)
        end if

        if( berror ) then
		write(6,*) 'error in section FLUX'
		stop 'error stop: ckflxa'
	end if

! initialize vectors (not strictly necessary)

	do k=1,kfluxm
	  do ii=1,3
	    iflux(ii,k) = 0
	  end do
	end do

! the double precision set up is done in flxini
! but since at this stage we do not have all the arrays set up
! we post-pone it until later

        end

!******************************************************************

	subroutine prflxa

        use fem_util
        use line_admin

	implicit none

	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii

	write(6,*)
	write(6,*) 'flux section :'
	write(6,*)
	write(6,*) 'nsect,kfluxm ',nsect,kfluxm
	write(6,*)

	ns = 0
	nnode = 0

        if(kfluxm .le. 0) return
	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'section : ',ns,ntotal
	  do i=ifirst,ilast
	    write(6,*) ipext(kflux(i)),(iflux(ii,i),ii=1,3)
	  end do
	end do

	end

!******************************************************************

	subroutine tsflxa

	implicit none

	integer i,ii

	write(6,*) '/kfluxc/'
	write(6,*) nsect,kfluxm
	write(6,*) (kflux(i),i=1,kfluxm)

	write(6,*) '/iflux/'
	write(6,*) ((iflux(ii,i),ii=1,3),i=1,kfluxm)

	end

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine wrflxa(it)

! administers writing of flux data

	use conz_common, only : cnv
	use ts
	use levels, only : nlvdi,nlv
        use defnames
        use para
        use output
	use flux_util
        use ioflux

	implicit none

	include 'param.h'

	integer it

	integer itend
	integer j,i,l,lmax,nlmax,ivar,nvers
	integer idtflx
	double precision az,azpar,rr

        integer nrm,nrs,nrt,nrc
	save nrm,nrs,nrt,nrc

        integer ia_out(4)
        save ia_out
        integer nbflx
        save nbflx
	integer ibarcl,iconz
	save ibarcl,iconz

        data nbflx /0/

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

		call init_output('itmflx','idtflx',ia_out)
		call increase_output(ia_out)
                if( .not. has_output(ia_out) ) nbflx = -1

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

        	call flux_alloc_arrays(nlvdi,nsect)

                !if( nsect .gt. nscflxdim ) then
                !  stop 'error stop wrflxa: dimension nscflxdim'
                !end if

		ibarcl = nint(getpar('ibarcl'))
		iconz = nint(getpar('iconz'))

		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		call fluxes_init(nlvdi,nsect,nlayers,nrm,masst)
		if( ibarcl .gt. 0 ) then
		  call fluxes_init(nlvdi,nsect,nlayers,nrs,saltt)
		  call fluxes_init(nlvdi,nsect,nlayers,nrt,tempt)
		end if
		if( iconz .eq. 1 ) then
		  call fluxes_init(nlvdi,nsect,nlayers,nrc,conzt)
		end if

                nbflx=ifemop('.flx','unform','new')
                if(nbflx.le.0) then
        	   stop 'error stop wrflxa : Cannot open FLX file'
		end if

	        nvers = 5
		idtflx = ia_out(1)
                call wfflx(nbflx,nvers,nsect,kfluxm,idtflx,nlmax,kflux,nlayers)

!               here we could also compute and write section in m**2

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( .not. is_over_output(ia_out) ) return

	call getaz(azpar)
	az = azpar

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)
	call fluxes_accum(nlvdi,nsect,nlayers,nrm,masst,fluxes)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,saltv)
	  call fluxes_accum(nlvdi,nsect,nlayers,nrs,saltt,fluxes)
	  ivar = 12
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,tempv)
	  call fluxes_accum(nlvdi,nsect,nlayers,nrt,tempt,fluxes)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,cnv)
	  call fluxes_accum(nlvdi,nsect,nlayers,nrc,conzt,fluxes)
	end if

!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( .not. next_output(ia_out) ) return

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

	ivar = 0
	call fluxes_aver(nlvdi,nsect,nlayers,nrm,masst,fluxes)
	call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call fluxes_aver(nlvdi,nsect,nlayers,nrs,saltt,fluxes)
	  call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)
	  ivar = 12
	  call fluxes_aver(nlvdi,nsect,nlayers,nrt,tempt,fluxes)
	  call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call fluxes_aver(nlvdi,nsect,nlayers,nrc,conzt,fluxes)
	  call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)
	end if

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	call fluxes_init(nlvdi,nsect,nlayers,nrm,masst)

	if( ibarcl .gt. 0 ) then
	  call fluxes_init(nlvdi,nsect,nlayers,nrs,saltt)
	  call fluxes_init(nlvdi,nsect,nlayers,nrt,tempt)
	end if

	if( iconz .eq. 1 ) then
	  call fluxes_init(nlvdi,nsect,nlayers,nrc,conzt)
	end if

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flxini

! initializes flux routines finally (wrapper for flx_init)

        use flux_util

	implicit none

        if(kfluxm .le. 0) return
	call flx_init(kfluxm,kflux,nsect,iflux)

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine fluxes_template(it)

! administers writing of flux data
!
! serves as a template for new variables
! please copy to extra file and adapt to your needs
!
! in this version multiple concentrations are written
!
! to change for adaptation:
!
! ncsdim	dimension of parameter arrays (here in param.h)
! conzv		parameters to be computed
! iconz		how many parameters actually needed
! csc		new extension for file
! ivar_base	base of variable numbering

	use conz_common, only : conzv
	use levels, only : nlvdi,nlv
        use defnames
        use para
	use flux_util
        use ioflux

	implicit none

	integer it

	include 'param.h'

	integer itend
	integer j,i,k,l,lmax,nlmax,ivar,nvers,ivar_base
	integer iconz
	double precision az,azpar,rr

        integer idtflx,itflx,itmflx,nbflx
        save idtflx,itflx,itmflx,nbflx

        data nbflx /0/

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

                idtflx = nint(dgetpar('idtflx'))
                itmflx = nint(dgetpar('itmflx'))
                itend = nint(dgetpar('itend'))
		iconz = nint(dgetpar('iconz'))	!computing concentrations?

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( idtflx .le. 0 ) nbflx = -1
                if( itmflx .gt. itend ) nbflx = -1
                if( iconz .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

		!be sure that other arrays are already allocated!!!!
        	call flux_alloc_conz_arrays(nlvdi,nsect,iconz)

                !if( nsect .gt. nscflxdim ) then
                !  stop 'error stop fluxes_template: dimension nscflxdim'
                !end if

                itflx = itmflx + idtflx
		itmflx = itmflx + 1	!start from next time step

		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		do k=1,iconz
		  call fluxes_init(nlvdi,nsect,nlayers,nrcc(k),cflux(0,1,1,k))
		end do

                nbflx=ifemop('.csc','unform','new')
                if(nbflx.le.0) then
        	   stop 'error stop wrflxa : Cannot open csc file'
		end if

	        nvers = 5
                call wfflx(nbflx,nvers,nsect,kfluxm,idtflx,nlmax,kflux,nlayers)

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( it .lt. itmflx ) return

	iconz = nint(getpar('iconz'))
	call getaz(azpar)
	az = azpar
	ivar_base = 200		!base of variable numbering

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	do k=1,iconz
	  ivar = ivar_base + k
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,conzv(1,1,k))
	  call fluxes_accum(nlvdi,nsect,nlayers,nrcc(k),cflux(0,1,1,k),fluxes)
	end do

!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( it .lt. itflx ) return
        itflx=itflx+idtflx

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

	do k=1,iconz
	  ivar = ivar_base + k
	  call fluxes_aver(nlvdi,nsect,nlayers,nrcc(k),cflux(0,1,1,k),fluxes)
	  call wrflx(nbflx,it,nlvdi,nsect,ivar,nlayers,fluxes)
	end do

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	do k=1,iconz
	  call fluxes_init(nlvdi,nsect,nlayers,nrcc(k),cflux(0,1,1,k))
	end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!******************************************************************
!******************************************************************
!******************************************************************

!==================================================================
        end module flux
!==================================================================
