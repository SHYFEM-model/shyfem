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

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	character*80 descrp
	common /descrp/ descrp

        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real hlv(1)
	common /hlv/hlv
        real hev(1)
	common /hev/hev

	integer itmout,ierr
	real href,hzoff
	integer date,time
	character*80 title,femver

	integer iround,ideffi
        integer wfout,wrout
	real getpar
	double precision dgetpar
	integer ifemop

	integer idtout,itout
	integer icall,nbout,nvers
	save idtout,itout
	save icall,nvers,nbout
	data icall,nvers,nbout /0,2,0/


	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		idtout=iround(getpar('idtout'))
		itmout=iround(getpar('itmout'))

		call adjust_itmidt(itmout,idtout)

		if( idtout .le. 0 ) icall = -1
		if( icall .eq. -1 ) return
		
		itout = itmout

		nbout = ifemop('.ous','unformatted','new')
		if(nbout.le.0) goto 77

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

		!call wfous(nbout,nvers,nkn,nel,nlv,href,hzoff,descrp,ierr)
        	!if(ierr.ne.0.) goto 78
		!call wsous(nbout,ilhv,hlv,hev,ierr)
        	!if(ierr.ne.0.) goto 75
	end if

	icall = icall + 1

	if( it .lt. itout ) return

	!call wrous(nbout,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
	call ous_write_record(nbout,it,nlvdim,ilhv,znv,zenv
     +					,utlnv,vtlnv,ierr)
	if(ierr.ne.0.) goto 79

	itout=itout+idtout

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

