c $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
c
c reading and interpolation of external files
c
c contents :
c
c revision log :
c
c 29.10.2012    ggu     created from scratch
c 17.06.2013    ggu     do not pass function into subroutine
c 02.07.2014    ggu     new framework finished
c 10.07.2014    ggu     only new file format allowed
c 22.02.2016    ggu&eps new files for generic tracer (nvar>1)
c
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open(name,it,nkn,nlv,iunit)
	integer iunit(3)
	character*(*) name
	  call ts_file_open_1(name,it,nkn,nlv,iunit)
	end

	subroutine ts_file_descrp(iunit,name)
	use intp_fem_file
	integer iunit(3)
	character*(*) name
	  call iff_set_description(iunit(1),0,name)
	end

	subroutine ts_next_record(it,iunit,nlvddi,nkn,nlv,value)
	include 'param.h'
	integer iunit(3)
	real value(nlvddi,nkn)
	  call ts_next_record_1(it,iunit,nlvddi,nkn,nlv,value)
	end

	subroutine ts_file_close(info)
	integer info(3)
	  call ts_file_close_1(info)
	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open_1(file,it,np,nlv,iunit)

c opens T/S file

	use intp_fem_file

	implicit none

	include 'param.h'

	character*(*) file		!name of file
	integer it
	integer np			!number of points expected
	integer nlv
	integer iunit(3)		!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer id
	double precision dtime
	integer nodes(1)
	real vconst(1)

	dtime = it
	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +                                  ,nodes,vconst,id)
!$OMP END CRITICAL

	iunit(1) = id

	return
   99	continue
	write(6,*) 'Cannot open file: ',file
	stop 'error stop ts_file_open: error file open'
	end

c*******************************************************************	

	subroutine ts_next_record_1(it,iunit,nlvddi,nkn,nlv,value)

	use intp_fem_file

	implicit none

	include 'param.h'

	integer it
	integer iunit(3)
	integer nlvddi
	integer nkn
	integer nlv
	real value(nlvddi,nkn)

	integer id,ldim,ndim,ivar
        real vmin,vmax
	double precision dtime
	character*80 string

c--------------------------------------------------------------
c initialize
c--------------------------------------------------------------

	id = iunit(1)

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	dtime = it
	ivar = 1
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading T/S values', it,dtime

	call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

c--------------------------------------------------------------
c some statistics
c--------------------------------------------------------------

        call conmima(nlvddi,value,vmin,vmax)

        !write(6,*) 'min/max: ',vmin,vmax

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine ts_file_close_1(iunit)

c closes T/S file

	use intp_fem_file

	implicit none

	integer iunit(3)

	integer id

	id = iunit(1)
	call iff_forget_file(id)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c****** new routines ***********************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine tracer_file_open(file,it,nvar,np,nlv,iunit)

c opens tracer file

	use intp_fem_file

	implicit none

	include 'param.h'

	character*(*) file		!name of file
	integer it
	integer nvar
	integer np			!number of points expected
	integer nlv
	integer iunit(3)		!unit number (return)

	integer nexp,lexp,nintp
	integer id
	double precision dtime
	integer nodes(nvar)
	real vconst(nvar)

	dtime = it
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +                                  ,nodes,vconst,id)
!$OMP END CRITICAL

	iunit(1) = id

	return
   99	continue
	write(6,*) 'Cannot open file: ',file
	stop 'error stop ts_file_open: error file open'
	end

c*******************************************************************	

	subroutine tracer_next_record(it,iunit,nvar,nlvddi,nkn,nlv,value)

c reads next record of tracer

	use intp_fem_file

	implicit none

	include 'param.h'

	integer it
	integer iunit(3)
	integer nvar
	integer nlvddi
	integer nkn
	integer nlv
	real value(nlvddi,nkn,nvar)

	integer id,ldim,ndim,ivar
	double precision dtime
	character*80 string

c--------------------------------------------------------------
c initialize
c--------------------------------------------------------------

	id = iunit(1)

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	dtime = it
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading tracer values', it,dtime

	call iff_read_and_interpolate(id,dtime)
	do ivar=1,nvar
	  call iff_time_interpolate(id,dtime,ivar,ndim,ldim
     +					,value(:,:,ivar))
	end do

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine tracer_file_descrp(iunit,name)

c sets description for file

	use intp_fem_file

	implicit none

	integer iunit(3)
	character*(*) name

	call iff_set_description(iunit(1),0,name)

	end

c*******************************************************************	

	subroutine tracer_file_close(iunit)

c closes tracer file

	use intp_fem_file

	implicit none

	integer iunit(3)

	integer id

	id = iunit(1)
	call iff_forget_file(id)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

