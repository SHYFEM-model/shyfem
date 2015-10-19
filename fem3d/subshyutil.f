!
! utility routines for shy file
!
! revision log :
!
! 15.10.2015    ggu     started routine
!
!****************************************************************
!****************************************************************
!****************************************************************

	subroutine shy_copy_basin_from_shy(id)

	use basin
	use shyfile

	implicit none

	integer id

	integer nk,ne,nlv,nvar

	call shy_get_params(id,nk,ne,nlv,nvar)
	call basin_init(nk,ne)

	call shy_get_elemindex(id,nen3v)
	call shy_get_coords(id,xgv,ygv)
	call shy_get_depth(id,hm3v)
	call shy_get_extnumbers(id,ipev,ipv)
	call shy_get_areacode(id,iarv,iarnv)

	end

!****************************************************************

	subroutine shy_copy_basin_to_shy(id)

	use basin
	use shyfile

	implicit none

	integer id

	call shy_set_elemindex(id,nen3v)
	call shy_set_coords(id,xgv,ygv)
	call shy_set_depth(id,hm3v)
	call shy_set_extnumbers(id,ipev,ipv)
	call shy_set_areacode(id,iarv,iarnv)

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine shy_open_output_file(file,nlv,nvar,ftype,id)

	use basin
	use shyfile

	implicit none

	character*(*) file	!file name
	integer nlv		!vertical dimension
	integer nvar		!number of variables
	integer ftype		!type of file: 1=ous 2=nos
	integer id		!id of opened file (return)

	include 'simul.h'

	integer iunit,ierr
	integer date,time
	character*80 title
	character*80 femver

	double precision dgetpar

c-----------------------------------------------------
c open file
c-----------------------------------------------------

	iunit = shy_open_file(file)

c-----------------------------------------------------
c collect needed data
c-----------------------------------------------------

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))
        title = descrp
        call get_shyfem_version(femver)

c-----------------------------------------------------
c write header of file
c-----------------------------------------------------

	call shy_init(iunit,id)
	call shy_set_params(id,nkn,nel,nlv,nvar)
        call shy_set_ftype(id,ftype)
        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

	call shy_alloc_arrays(id)
	call shy_copy_basin_to_shy(id)
	call shy_write_header(id,ierr)

c-----------------------------------------------------
c error check
c-----------------------------------------------------

	if( ierr /= 0 ) then
	  write(6,*) 'error writing header of file ',ierr
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_open_output_file: writing header'
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	end

!****************************************************************

	subroutine shy_write_output_record(id,dtime,ivar,n,nlv
     +					,nlvdi,c)

	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	integer n,nlv,nlvdi
	real c(nlvdi,n)

	integer ierr
	character*80 file

	call shy_write_record(id,dtime,ivar,n,1,nlv,nlvdi,c,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'error writing output file ',ierr
	  call shy_get_filename(id,file)
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_write_output_file: writing record'
	end if

	end

!****************************************************************

	subroutine shy_write_scalar_record(id,dtime,ivar,c)

	use basin
	use levels

	implicit none

	integer id
	double precision dtime
	integer ivar
	real c(nlvdi,nkn)

	call shy_write_output_record(id,dtime,ivar,nkn,nlv,nlvdi,c)

	end

!****************************************************************

	subroutine shy_write_scalar_record2d(id,dtime,ivar,c)

	use basin
	use levels

	implicit none

	integer id
	double precision dtime
	integer ivar
	real c(nkn)

	call shy_write_output_record(id,dtime,ivar,nkn,1,1,c)

	end

!****************************************************************

