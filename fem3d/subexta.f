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
c
c******************************************************************

	subroutine mod_ext(mode)

	implicit none

	integer mode

	include 'modules.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

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

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus

	knausm = 0

	end

c******************************************************************

	subroutine rdexta

	implicit none

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus

	logical handlesec
	integer nexdim
	integer nrdveci

	if( .not. handlesec('extra') ) return

	call getdim('nexdim',nexdim)

	knausm = nrdveci(knaus,nexdim)

	if( knausm .lt. 0 ) then
	  if( knausm .eq. -1 ) then
	    write(6,*) 'dimension error nexdim in section $extra : '
     +				,nexdim
	  else
	    write(6,*) 'read error in section $extra'
	  end if
	  stop 'error stop rdexta'
	end if

	end

c******************************************************************

	subroutine ckexta

	implicit none

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus

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

	implicit none

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus

	integer i
	integer ipext

        if(knausm.le.0) return

        write(6,*)
        write(6,1009) knausm,(ipext(knaus(i)),i=1,knausm)
        write(6,*)

	return
 1009   format(' extra section : ',i5/(10i6))
	end

c******************************************************************

	subroutine tsexta

	implicit none

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus

	integer i

        write(6,*) '/knausc/'
        write(6,*) knausm
        write(6,*) (knaus(i),i=1,knausm)

	end

c******************************************************************

	subroutine wrexta(it)

c writes and administers ext file

	implicit none

	integer it

        integer knausm,knaus(1)
        common /knausc/ knausm,knaus
	real xv(3,1)
	common /xv/xv

	integer itmext,itend
	real err,href,hzoff
	integer iround,ideffi
	real getpar
	double precision dgetpar
	real writ7h,wrrc7

	integer idtext,itext
	integer icall,nbext,nvers
	save idtext,itext
	save icall,nbext,nvers
	data icall,nbext,nvers /0,7,6/

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		idtext = nint(dgetpar('idtext'))
		itmext = nint(dgetpar('itmext'))
		itend = nint(dgetpar('itend'))
		if( knausm .le. 0 ) icall = -1
		if( idtext .le. 0 ) icall = -1
		if( itmext .gt. itend ) icall = -1
		if( icall .eq. -1 ) return

		nbext=ideffi('datdir','runnam','.ext','unform','new')
                if(nbext.le.0) goto 77

		itext = itmext
		href = getpar('href')
		hzoff = getpar('hzoff')
                err=writ7h(nbext,nvers,knausm,knaus,href,hzoff)
                if(err.ne.0.) goto 78
        end if

	icall = icall + 1

c write file ext

	if( it .lt. itext ) return

        err=wrrc7(nbext,nvers,it,knausm,knaus,xv)
        if(err.ne.0.) goto 79
        itext=itext+idtext

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
        integer knausm,knaus(1)
	real href,hzmin

	integer nmax
	parameter(nmax=50)

        character*80 descrp
	common /descrp/descrp
	real v1v(1)
	common /v1v/v1v

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

