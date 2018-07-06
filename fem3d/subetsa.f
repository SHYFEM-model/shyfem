c
c $Id: subetsa.f,v 1.6 2001/11/16 07:35:43 georg Exp $
c
c ETS file administration routines (deals with section $EXTTS)
c
c revision log :
c
c 24.01.2014    ggu     copied from subexta.f
c 30.05.2015    ggu     new read from STR file, using module ets
c 01.04.2016    ggu     bug fix for waves
c
c******************************************************************
c******************************************************************
c******************************************************************

        subroutine ets_read_section(n)

        use ets
        use nls

        integer n

        n = nls_read_ictable()

        if( n > 0 ) then
	  call ets_init_module(n)
          call nls_copy_ictable(n,nkets,chets)
        end if

        end subroutine ets_read_section

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine mod_ets(mode)

	implicit none

	integer mode

	include 'modules.h'
	include 'femtime.h'

	if( mode .eq. M_AFTER ) then
	   call wretsa(it)
	else if( mode .eq. M_INIT ) then
	   call inetsa
	else if( mode .eq. M_READ ) then
	   call rdetsa
	else if( mode .eq. M_CHECK ) then
	   call cketsa
	else if( mode .eq. M_SETUP ) then
	   call wretsa(it)
	else if( mode .eq. M_PRINT ) then
	   call pretsa
	else if( mode .eq. M_TEST ) then
	   call tsetsa
	else if( mode .eq. M_BEFOR ) then
c	   nothing
	else
	   write(6,*) 'unknown mode : ', mode
	   stop 'error stop mod_ets'
	end if

	end

c******************************************************************

	subroutine inetsa

	implicit none

	end

c******************************************************************

	subroutine rdetsa

	use ets

	implicit none

	integer n
	logical handlesec

	if( .not. handlesec('extts') ) return

        call ets_read_section(n)

	if( n .lt. 0 ) then
	  write(6,*) 'read error in section $extts'
	  write(6,*) 'please remember that the format is:'
	  write(6,*) 'node  ''text'''
	  write(6,*) 'text may be missing'
	  write(6,*) 'only one node per line is allowed'
	  stop 'error stop rdetsa'
	end if

	end

c******************************************************************

	subroutine cketsa

	use ets

	implicit none

	integer k,kext,kint
	integer ipint
	logical bstop

	bstop = .false.

        do k=1,nets
	   kext = nkets(k)
           kint = ipint(kext)
           if(kint.le.0) then
                write(6,*) 'section EXTTS : node not found ',kext
                bstop=.true.
           end if
           nkets(k)=kint
	   xets(k) = 0.
	   yets(k) = 0.
	   if( chets(k) == ' ' ) then
	     write(chets(k),'(a,i5)') 'Extra node ',k
	   end if
        end do

	if( bstop ) stop 'error stop: cketsa'

	end

c******************************************************************

	subroutine pretsa

	use ets

	implicit none

	integer i,k
	real x,y
	character*80 s
	integer ipext

        if(nets.le.0) return

        write(6,*)
        write(6,*) 'extts section: start ',nets
	do i=1,nets
	  k = ipext(nkets(i))
	  x = xets(i)
	  y = yets(i)
	  s = chets(i)
          write(6,1009) i,k,x,y,s(1:44)
          !write(6,1008) i,k,s
	end do
        write(6,*) 'extts section: end'
        write(6,*)

	return
 1009   format(i3,i7,2e12.4,1x,a)
 1008   format(i3,i10,1x,a)
	end

c******************************************************************

	subroutine tsetsa

	use ets

	implicit none

	integer i

        write(6,*) '/ets_/'
        write(6,*) nets
        write(6,*) (nkets(i),i=1,nets)
        write(6,*) (xets(i),i=1,nets)
        write(6,*) (yets(i),i=1,nets)
        write(6,*) (chets(i),i=1,nets)

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine wretsa(it)

c writes and administers ets file

	use mod_waves
	use mod_depth
	use mod_ts
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use ets

	implicit none

	include 'simul.h'

	integer it

	integer ierr
	integer date,time
	integer nvers,nvar,ivar
	character*80 title,femver

	real, allocatable :: waves(:,:)
        integer, allocatable :: il4kv(:)

	integer ifemop
	real getpar
	double precision dgetpar
	logical has_output,next_output
	logical has_waves,bwave
	save bwave
	integer nvars

	integer nbext
	integer ia_out(4)
	save ia_out

	integer icall
	data icall /0/
	save icall

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
        	call init_output('itmext','idtext',ia_out)
        	if( .not. has_output(ia_out) ) icall = -1
		if( nets .le. 0 ) icall = -1

                if( icall .eq. -1 ) return

                nbext = ifemop('.ets','unformatted','new')
                if(nbext.le.0) goto 77
		ia_out(4) = nbext

		allocate( outets(nlvdi,nets) )
		allocate( out4ets(4,nets) )

		nvers = 1
		nvar = 5
	        bwave = has_waves()
	        if ( bwave ) nvar = nvar + 1

                date = nint(dgetpar('date'))
                time = nint(dgetpar('time'))
                title = descrp
                call get_shyfem_version_and_commit(femver)

		call ets_setup

                call ets_init(nbext,nvers)
                call ets_set_title(nbext,title)
                call ets_set_date(nbext,date,time)
                call ets_set_femver(nbext,femver)
                call ets_write_header(nbext,nets,nlv,nvar,ierr)
                if(ierr.gt.0) goto 78
                call ets_write_header2(nbext,ilets,hlv,hets
     +					,nkets,xets,yets,chets,ierr)
                if(ierr.gt.0) goto 75
        end if

	icall = icall + 1

	if( .not. next_output(ia_out) ) return

	nbext =	ia_out(4)

	ivar = 1
	call routets(1,ilhkv,znv,outets)
        call ets_write_record(nbext,it,ivar,1,ilets,outets,ierr)
        if(ierr.ne.0.) goto 79

	if ( bwave ) then
  	  ivar = 31
	  nvars = 4
	  allocate(waves(4,nkn),il4kv(nkn))
	  waves(1,:) = waveh
	  waves(2,:) = wavep
	  waves(3,:) = wavepp
	  waves(4,:) = waved
	  il4kv = nvars
	  call routets(nvars,il4kv,waves,out4ets)
          call ets_write_record(nbext,it,ivar,nvars,il4ets,out4ets,ierr)
          if(ierr.ne.0.) goto 79
	  deallocate(waves,il4kv)
        end if

	ivar = 6
	call routets(nlvdi,ilhkv,uprv,outets)
        call ets_write_record(nbext,it,ivar,nlvdi,ilets,outets,ierr)
        if(ierr.ne.0.) goto 79

	ivar = 7
	call routets(nlvdi,ilhkv,vprv,outets)
        call ets_write_record(nbext,it,ivar,nlvdi,ilets,outets,ierr)
        if(ierr.ne.0.) goto 79

	ivar = 11
	call routets(nlvdi,ilhkv,saltv,outets)
        call ets_write_record(nbext,it,ivar,nlvdi,ilets,outets,ierr)
        if(ierr.ne.0.) goto 79

	ivar = 12
	call routets(nlvdi,ilhkv,tempv,outets)
        call ets_write_record(nbext,it,ivar,nlvdi,ilets,outets,ierr)
        if(ierr.ne.0.) goto 79

	return
   77   continue
	write(6,*) 'Error opening ETS file :'
	stop 'error stop : wretsa'
   75   continue
	write(6,*) 'Error writing second header of ETS file'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop : wretsa'
   78   continue
	write(6,*) 'Error writing first header of ETS file'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop : wretsa'
   79   continue
	write(6,*) 'Error writing file ETS'
	write(6,*) 'unit,ivar,ierr :',nbext,ivar,ierr
	stop 'error stop : wretsa'
	end

c******************************************************************

	subroutine ets_setup

	use coordinates
	use mod_depth
	use levels
	use ets

	implicit none

	integer i,k

	do i=1,nets
	  k = nkets(i)
	  ilets(i) = ilhkv(k)
	end do
	il4ets = 4

	call routets(1,ilhkv,hkv,hets)
	call routets(1,ilhkv,xgeov,xets)
	call routets(1,ilhkv,ygeov,yets)

	end

c******************************************************************

	subroutine routets(nlv,ilhkv,value,out)

c extracts nodal information and stores it in array

	use ets

	implicit none

	integer nlv		!vertical dimension of data
	integer ilhkv(*)
	real value(nlv,*)
	real out(nlv,nets)

	integer i,k,l,lmax

	do i=1,nets
	  k = nkets(i)
	  lmax = ilhkv(k)
	  lmax = min(lmax,nlv)
	  do l=1,lmax
	    out(l,i) = value(l,k)
	  end do
	end do

	end

c******************************************************************

