c
c $Id: subexta.f,v 1.6 2001/11/16 07:35:43 georg Exp $
c
c extra file administration routines
c
c revision log :
c
c 08.05.1998	ggu	changed read of node numbers through nrdveci
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 01.02.2000    ggu     module framework introduced
c 20.05.2015    ggu     modules introduced
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module extra
!==================================================================

        implicit none

        integer, save :: knausm = 0
        integer, save, allocatable :: knaus(:)
        character*80, save, allocatable :: chext(:)

!==================================================================
        contains
!==================================================================

!==================================================================
        end module extra
!==================================================================

	subroutine extra_read_section(n)

	use extra
	use nls

	integer n

	!n = nls_read_vector()
	n = nls_read_ictable()
	knausm = n

	if( n > 0 ) then
	  allocate(knaus(n))
	  allocate(chext(n))
	  !call nls_copy_int_vect(n,knaus)
	  call nls_copy_ictable(n,knaus,chext)
	end if

	end subroutine extra_read_section

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine mod_ext(mode)

	implicit none

	integer mode

	include 'modules.h'
	include 'femtime.h'

	double precision dtime

	if( mode .eq. M_AFTER ) then
	   dtime = t_act
	   call wrexta(dtime)
	else if( mode .eq. M_INIT ) then
	   call inexta
	else if( mode .eq. M_READ ) then
	   call rdexta
	else if( mode .eq. M_CHECK ) then
	   call ckexta
	else if( mode .eq. M_SETUP ) then
	   dtime = t_act
	   call wrexta(dtime)			!ggu 7/5/2001 -> write it=0
	else if( mode .eq. M_PRINT ) then
	   call prexta
	else if( mode .eq. M_TEST ) then
	   call tsexta
	else if( mode .eq. M_BEFOR ) then
c	   nothing
	else
	   write(6,*) 'unknown mode : ', mode
	   stop 'error stop mod_ext'
	end if

	end

c******************************************************************

	subroutine inexta

	implicit none

	end

c******************************************************************

	subroutine rdexta

	use extra

	implicit none

	integer n
	logical handlesec

	if( .not. handlesec('extra') ) return

	call extra_read_section(n)

	if( n .lt. 0 ) then
	  write(6,*) 'read error in section $extra'
	  stop 'error stop rdexta'
	end if

	end

c******************************************************************

	subroutine ckexta

	use extra

	implicit none

	integer k,knode
	integer ipint
	logical bstop

	bstop = .false.

        do k=1,knausm
           knode=ipint(knaus(k))                !$$EXTINW
           if(knode.le.0) then
                write(6,*) 'section EXTRA : node not found ',knaus(k)
                bstop=.true.
           end if
           knaus(k)=knode
        end do

	if( bstop ) stop 'error stop: ckexta'

	end

c******************************************************************

	subroutine prexta

	use extra

	implicit none

	integer i
	integer ipext

        if(knausm.le.0) return

        write(6,*)
        write(6,*) 'extra section : ',knausm
	do i=1,knausm
          write(6,*) i,knausm,ipext(knaus(i)),'  ',trim(chext(i))
	end do
        write(6,*)

	end

c******************************************************************

	subroutine tsexta

	use extra

	implicit none

	integer i
	integer ipext

        write(6,*) '/knausc/'
        write(6,*) knausm
	do i=1,knausm
          write(6,*) i,knausm,ipext(knaus(i)),'  ',trim(chext(i))
	end do

	end

c******************************************************************

	subroutine wrexta(dtime)

c writes and administers ext file

	use mod_hydro
	use mod_hydro_print
	use mod_ts
	use mod_depth
	use basin
	use levels
	use extra

	implicit none

	double precision dtime

	include 'simul.h'

	integer nbext,ierr
	integer nvar,ivar,m,j,k,iv
	integer isalt,itemp
	real href,hzmin
	double precision atime,atime0
	character*80 femver,title
	real hdep(knausm)
	real x(knausm)
	real y(knausm)
	real vals(nlv,knausm,3)
	integer, save, allocatable :: il(:)
	character*80 strings(knausm)

	integer ideffi
	real getpar
	logical has_output_d,next_output_d

	double precision, save :: da_out(4) = 0
	integer, save :: icall = 0

	if( icall .eq. -1 ) return

	nvar = 3		!this must be adjusted

c--------------------------------------------------------------
c initialization
c--------------------------------------------------------------

	if( icall .eq. 0 ) then
          call init_output_d('itmext','idtext',da_out)
	  call assure_initial_output_d(da_out)
          if( .not. has_output_d(da_out) ) icall = -1
	  if( knausm .le. 0 ) icall = -1
	  if( icall .eq. -1 ) return

	  nbext=ideffi('datdir','runnam','.ext','unform','new')
          if(nbext.le.0) goto 99
	  da_out(4) = nbext

          call ext_write_header(nbext,0,knausm,nlv,nvar,ierr)
          if( ierr /= 0 ) goto 98

	  allocate(il(knausm))
	  title = descrp
	  href = getpar('href')
	  hzmin = getpar('hzmin')
	  do j=1,knausm
	    k = knaus(j)
	    hdep(j) = hkv_max(k)
	    x(j) = xgv(k)
	    y(j) = ygv(k)
	    il(j) = ilhkv(k)
	    strings(j) = chext(j)
	  end do
	  call get_shyfem_version(femver)
	  call get_absolute_ref_time(atime0)
          call ext_write_header2(nbext,0,knausm,nlv
     +                          ,atime0
     +                          ,href,hzmin,title,femver
     +                          ,knaus,hdep,il,x,y,strings,hlv
     +                          ,ierr)
          if( ierr /= 0 ) goto 98
        end if

	icall = icall + 1

c--------------------------------------------------------------
c write file ext
c--------------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

        nbext = nint(da_out(4))
	call get_absolute_act_time(atime)

c	-------------------------------------------------------
c	barotropic velocities
c	-------------------------------------------------------

	iv = 0
	ivar = 0
	m = 3
	do j=1,knausm
	  k = knaus(j)
	  vals(1,j,1) = up0v(k)
	  vals(1,j,2) = vp0v(k)
	  vals(1,j,3) = znv(k)
	end do
        call ext_write_record(nbext,0,atime,knausm,nlv
     +                                  ,ivar,m,il,vals,ierr)
        if( ierr /= 0 ) goto 97

c	-------------------------------------------------------
c	velocities
c	-------------------------------------------------------

	iv = iv + 1
	ivar = 2
	m = 2
	do j=1,knausm
	  k = knaus(j)
	  vals(:,j,1) = uprv(:,k)
	  vals(:,j,2) = vprv(:,k)
	end do
        call ext_write_record(nbext,0,atime,knausm,nlv
     +                                  ,ivar,m,il,vals,ierr)
        if( ierr /= 0 ) goto 97

c	-------------------------------------------------------
c	temperature
c	-------------------------------------------------------

	m = 1

	itemp=nint(getpar('itemp'))
	if( itemp > 0 ) then
	  iv = iv + 1
	  ivar = 12
	  do j=1,knausm
	    k = knaus(j)
	    vals(:,j,1) = tempv(:,k)
	  end do
          call ext_write_record(nbext,0,atime,knausm,nlv
     +                                  ,ivar,m,il,vals,ierr)
          if( ierr /= 0 ) goto 97
	end if

c	-------------------------------------------------------
c	salinity
c	-------------------------------------------------------

	isalt=nint(getpar('isalt'))
	if( isalt > 0 ) then
	  iv = iv + 1
	  ivar = 11
	  do j=1,knausm
	    k = knaus(j)
	    vals(:,j,1) = saltv(:,k)
	  end do
          call ext_write_record(nbext,0,atime,knausm,nlv
     +                                  ,ivar,m,il,vals,ierr)
          if( ierr /= 0 ) goto 97
	end if

	if( iv > nvar ) goto 91

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   91   continue
	write(6,*) 'iv,nvar: ',iv,nvar
	write(6,*) 'iv cannot be greater than nvar'
	stop 'error stop wrexta: internal error (1)'
   99   continue
	write(6,*) 'Error opening EXT file :'
	stop 'error stop wrexta: opening ext file'
   98   continue
	write(6,*) 'Error writing header of EXT file'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop wrexta: writing ext header'
   97   continue
	write(6,*) 'Error writing file EXT'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop wrexta: writing ext record'
	end

c*********************************************************

