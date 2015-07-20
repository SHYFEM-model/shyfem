c
c $Id: subousa.f,v 1.2 2005/11/03 16:59:25 georg Exp $
c
c OUS file administration routines
c
c contents :
c
c subroutine wrousa
c
c revision log :
c
c 26.01.1998	ggu	$$ITMOUT - adjust itmout for first write
c 01.09.2003	ggu	new routine wrousa
c 02.09.2003	ggu	bug fix in wrousa: save nbout
c 25.11.2004	ggu	in new file subousa.f
c 18.05.2005	ggu	initial itmout is changed (now itanf)
c 20.01.2014	ggu	new calls for ous writing implemented
c
c********************************************************

	subroutine wrousa

c writes and administers ous file

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'

	include 'simul.h'


	integer itmout,ierr
	real href,hzoff
	integer date,time
	character*80 title,femver

	integer iround,ideffi
        integer wfout,wrout
	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output,next_output

	integer idtout,itout
	integer icall,nbout,nvers
	integer ia_out(4)
	save ia_out
	save idtout,itout
	save icall,nvers,nbout
	data icall,nvers,nbout /0,2,0/


	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		call init_output('itmout','idtout',ia_out)

		if( .not. has_output(ia_out) ) icall = -1
		if( icall .eq. -1 ) return
		
		nbout = ifemop('.ous','unformatted','new')
		if(nbout.le.0) goto 77
		ia_out(4) = nbout

		href=getpar('href')             !reference level
		hzoff=getpar('hzoff')           !minimum depth
	        date = nint(dgetpar('date'))
	        time = nint(dgetpar('time'))
	        title = descrp
	        call get_shyfem_version(femver)

	        call ous_init(nbout,nvers)
	        call ous_set_title(nbout,title)
	        call ous_set_date(nbout,date,time)
	        call ous_set_femver(nbout,femver)
	        call ous_set_hparams(nbout,href,hzoff)
	        call ous_write_header(nbout,nkn,nel,nlv,ierr)
	        if(ierr.gt.0) goto 78
	        call ous_write_header2(nbout,ilhv,hlv,hev,ierr)
	        if(ierr.gt.0) goto 75
	end if

	icall = icall + 1

	if( .not. next_output(ia_out) ) return

	call ous_write_record(nbout,it,nlvdi,ilhv,znv,zenv
     +					,utlnv,vtlnv,ierr)
	if(ierr.ne.0.) goto 79

	return
   77   continue
	write(6,*) 'Error opening OUS file :'
	stop 'error stop : wrousa'
   78   continue
	write(6,*) 'Error writing first header of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
   75   continue
	write(6,*) 'Error writing second header of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
   79   continue
	write(6,*) 'Error writing data record of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
	end

c********************************************************

