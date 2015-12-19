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

	n = nls_read_vector()
	knausm = n

	if( n > 0 ) then
	  allocate(knaus(n))
	  call nls_copy_int_vect(n,knaus)
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

	if( mode .eq. M_AFTER ) then
	   call wrexta(it)
	else if( mode .eq. M_INIT ) then
	   call inexta
	else if( mode .eq. M_READ ) then
	   call rdexta
	else if( mode .eq. M_CHECK ) then
	   call ckexta
	else if( mode .eq. M_SETUP ) then
c	   nothing
	   call wrexta(it)			!ggu 7/5/2001 -> write it=0
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

	!stop 'error stop subexta: extra not supported anymore'
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

	!stop 'error stop subexta: extra not supported anymore'
        if(knausm.le.0) return

        write(6,*)
        write(6,1009) knausm,(ipext(knaus(i)),i=1,knausm)
        write(6,*)

	return
 1009   format(' extra section : ',i5/(10i6))
	end

c******************************************************************

	subroutine tsexta

	use extra

	implicit none

	integer i

        write(6,*) '/knausc/'
        write(6,*) knausm
        write(6,*) (knaus(i),i=1,knausm)

	end

c******************************************************************

	subroutine wrexta(it)

c writes and administers ext file

	use mod_hydro_print
	use extra

	implicit none

	integer it

	integer nbext
	real err,href,hzoff
	integer iround,ideffi
	real getpar
	double precision dgetpar
	real writ7h,wrrc7
	logical has_output,next_output

	integer ia_out(4)
	save ia_out
	integer icall,nvers
	save icall,nvers
	data icall,nvers /0,6/

	!stop 'error stop subexta: extra not supported anymore'
	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
                call init_output('itmext','idtext',ia_out)
		call assure_initial_output(ia_out)
                if( .not. has_output(ia_out) ) icall = -1
		if( knausm .le. 0 ) icall = -1
		if( icall .eq. -1 ) return

		nbext=ideffi('datdir','runnam','.ext','unform','new')
                if(nbext.le.0) goto 77
		ia_out(4) = nbext

		href = getpar('href')
		hzoff = getpar('hzoff')
                err=writ7h(nbext,nvers,knausm,knaus,href,hzoff)
                if(err.ne.0.) goto 78
        end if

	icall = icall + 1

c write file ext

        if( .not. next_output(ia_out) ) return

        nbext = ia_out(4)

        err=wrrc7(nbext,nvers,it,knausm,knaus,xv)
        if(err.ne.0.) goto 79

	return
   77   continue
	write(6,*) 'Error opening EXT file :'
	stop 'error stop : wrexta'
   78   continue
	write(6,*) 'Error writing first record of EXT file'
	write(6,*) 'unit,err :',nbext,iround(err)
	stop 'error stop : wrexta'
   79   continue
	write(6,*) 'Error writing file EXT'
	write(6,*) 'unit,err :',nbext,iround(err)
	stop 'error stop : wrexta'
	end

c*********************************************************

        function writ7h(iunit,nvers,knausm,knaus,href,hzmin)

c writes first record of file 7 from main
c
c ...depth has not to be passed but is computed in routine
c ...pass only dummy vector
c ...ndim is dummy argument

	implicit none

	real writ7h
	integer iunit,nvers
        integer knausm,knaus(knausm)
	real href,hzmin
	real v1v(knausm)

	integer nmax
	parameter(nmax=50)

	include 'simul.h'

	integer i,n,ndim
	real hmin,hmax
        real f(nmax)

	integer igtdep
	real writ7

        do i=1,knausm
          n=igtdep(knaus(i),f,nmax)
          call mima(f,n,hmin,hmax)
          v1v(i)=hmax+href
        end do

	ndim = knausm	!dummy
        writ7h=writ7(iunit,ndim,nvers,knausm,knaus,v1v
     +                          ,href,hzmin,descrp)

        end

c************************************************************

