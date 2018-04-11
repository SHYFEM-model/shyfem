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
c 18.05.2005	ggu	initial itmout is changed
c 20.01.2014	ggu	new calls for ous writing implemented
c 15.10.2015	ggu	added new calls for shy file format
c 11.04.2018	ggu	mpi version is ready and working
c
c********************************************************

	subroutine wrousa

c writes and administers ous file

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shyfile
	use shympi

	implicit none

	include 'simul.h'

	logical bdebug
	integer itmout,ierr
	real href,hzoff
	integer date,time
	character*80 title,femver

	integer iround,ideffi
        integer wfout,wrout
	integer nvar,ftype,id
	double precision dtime
	character*80 file

	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output_d,next_output_d

	integer idtout,itout
	integer icall,nbout,nvers
	double precision, save :: da_out(4)
	save idtout,itout
	save icall,nvers,nbout
	data icall,nvers,nbout /0,2,0/

	bdebug = .true.
	bdebug = .false.

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		da_out = 0.
		call init_output_d('itmout','idtout',da_out)

		if( .not. has_output_d(da_out) ) icall = -1
		if( icall .eq. -1 ) return

		if( has_output_d(da_out) ) then
		  nvar = 4
		  ftype = 1
		  call shy_make_output_name('.hydro.shy',file)
		  call shy_open_output_file(file,3,nlv,nvar,ftype,id)
                  call shy_set_simul_params(id)
                  call shy_make_header(id)
		  da_out(4) = id
		end if
	end if

	icall = icall + 1

	if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  call get_act_dtime(dtime)
	  call shy_write_hydro_records(id,dtime,nlvdi,znv,zenv
     +					,utlnv,vtlnv)
	end if

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

