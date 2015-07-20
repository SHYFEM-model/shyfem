c
c $Id: subouta.f,v 1.7 2004/12/02 09:11:33 georg Exp $
c
c OUT file administration routines
c
c contents :
c
c subroutine wrouta
c
c revision log :
c
c 26.01.1998	ggu	$$ITMOUT - adjust itmout for first write
c 01.09.2003	ggu	new routine wrousa
c 02.09.2003	ggu	bug fix in wrousa: save nbout
c 25.11.2004	ggu	wrousa in new file
c
c********************************************************

	subroutine wrouta

c writes and administers out file

	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'femtime.h'
	include 'simul.h'


	integer itmout,ierr
	real href,hzoff

	integer iround,ideffi
        integer wfout,wrout
	real getpar

	integer idtout,itout
	integer icall,nbout,nvers
	save idtout,itout
	save icall,nbout,nvers
	data icall,nbout,nvers /0,8,6/

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		idtout=iround(getpar('idtout'))
		itmout=iround(getpar('itmout'))

		if( itmout .le. itanf ) itmout = itmout + idtout !$$ITMOUT

		if( idtout .le. 0 ) icall = -1
		if( itmout .gt. itend ) icall = -1
		if( icall .eq. -1 ) return
		
		nbout=ideffi('datdir','runnam','.out','unform','new')
		if(nbout.le.0) goto 77

		itout = itmout
		href=getpar('href')             !reference level
		hzoff=getpar('hzoff')           !minimum depth

        	ierr=wfout(nbout,nvers,itanf,itend,idt,idtout
     +                          ,href,hzoff,descrp)
        	if(ierr.ne.0.) goto 78
	end if

	icall = icall + 1

	if( it .lt. itout ) return

	ierr=wrout(nbout,nvers,it,nkn,nel,xv,zenv,unv,vnv)
	if(ierr.ne.0.) goto 79
	itout=itout+idtout

	return
   77   continue
	write(6,*) 'Error opening OUT file :'
	stop 'error stop : wrouta'
   78   continue
	write(6,*) 'Error writing first record of OUT file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrouta'
   79   continue
	write(6,*) 'Error writing file OUT'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrouta'
	end

c********************************************************

