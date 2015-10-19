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

	subroutine shy_copy_levels_from_shy(id)

	use levels
	use shyfile

	implicit none

	integer id

	call shy_get_layers(id,hlv)
	call shy_get_layerindex(id,ilhv,ilhkv)

	end

!****************************************************************

	subroutine shy_copy_levels_to_shy(id)

	use levels
	use shyfile

	implicit none

	integer id

	call shy_set_layers(id,hlv)
	call shy_set_layerindex(id,ilhv,ilhkv)

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine shy_get_tstart(file,datetime,dtime,bok)

	use shyfile

	implicit none

	character*(*) file
	integer datetime(2)
	double precision dtime
	logical bok			!could read time ?

	integer iunit,id,ierr
	integer ivar,n,m,lmax
	integer date,time

	bok = .false.
	dtime = -1.

	id = shy_init(file)
	if( id == 0 ) return
	call shy_skip_header(id,ierr)
	if( ierr /= 0 ) return

	call shy_get_date(id,date,time)
	datetime = (/date,time/)

	call shy_skip_record(id,dtime,ivar,n,m,lmax,ierr)
	if( ierr /= 0 ) return

	bok = .true.
	
	end

!****************************************************************

	subroutine shy_get_tend(file,datetime,dtime,bok)

	use shyfile

	implicit none

	character*(*) file
	integer datetime(2)
	double precision dtime
	logical bok			!could read time ?

	integer iunit,id,ierr
	integer ivar,n,m,lmax
	integer date,time

	bok = .false.
	dtime = -1.

	id = shy_init(file)
	if( id == 0 ) return
	call shy_skip_header(id,ierr)
	if( ierr /= 0 ) return

	call shy_get_date(id,date,time)
	datetime = (/date,time/)

	do while( ierr == 0 )
	  call shy_skip_record(id,dtime,ivar,n,m,lmax,ierr)
	end do
	if( ierr > 0 ) return

	bok = .true.

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine shy_make_output_name(ext,file)

	character*(*) ext	!extension (with dot), example: .ext
	character*(*) file	!complete file name (return)

	character*80 simul

	call getfnm('runnam',simul)

	file = trim(simul) // trim(ext)

	end

!****************************************************************

	subroutine shy_open_output_file(file,npr,nlv,nvar,ftype,id)

	use basin
	use shyfile

	implicit none

	character*(*) file	!file name
	integer npr		!max number of parameters per node
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

	id = shy_init(file)

	if( id == 0 ) then
	  write(6,*) 'error opening file'
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_open_output_file: opening file'
	end if

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

	call shy_set_params(id,nkn,nel,npr,nlv,nvar)
        call shy_set_ftype(id,ftype)
        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

	call shy_alloc_arrays(id)
	call shy_copy_basin_to_shy(id)
	call shy_copy_levels_to_shy(id)
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
!****************************************************************
!****************************************************************

	subroutine shy_write_output_record(id,dtime,ivar,n,m,nlv
     +					,nlvdi,c)

	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	integer n,m,nlv,nlvdi
	real c(nlvdi,n)

	integer ierr
	character*80 file

	call shy_write_record(id,dtime,ivar,n,m,nlv,nlvdi,c,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'error writing output file ',ierr
	  call shy_get_filename(id,file)
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_write_output_file: writing record'
	end if

	end

!****************************************************************

	subroutine shy_write_scalar_record(id,dtime,ivar,nlvddi,c)

	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	integer nlvddi
	real c(nlvddi,*)

	integer iaux,nkn,nlv

	call shy_get_params(id,nkn,iaux,iaux,nlv,iaux)
	call shy_write_output_record(id,dtime,ivar,nkn,1,nlv,nlvddi,c)

	end

!****************************************************************

	subroutine shy_write_scalar_record2d(id,dtime,ivar,c)

	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	real c(*)

	integer iaux,nkn

	call shy_get_params(id,nkn,iaux,iaux,iaux,iaux)
	call shy_write_output_record(id,dtime,ivar,nkn,1,1,1,c)

	end

!****************************************************************

	subroutine shy_write_hydro_records(id,dtime,nlvddi,z,ze,u,v)

	use shyfile

	implicit none

	integer id
	double precision dtime
	integer nlvddi
	real z(*)
	real ze(3,*)
	real u(nlvddi,*)
	real v(nlvddi,*)

	integer iaux,nkn,nel,npr,nlv,ivar

	call shy_get_params(id,nkn,nel,npr,nlv,iaux)

	ivar = 1
	call shy_write_output_record(id,dtime,ivar,nkn,1,1,1,z)
	call shy_write_output_record(id,dtime,ivar,nel,3,1,1,ze)
	ivar = 3
	call shy_write_output_record(id,dtime,ivar,nel,1,nlv,nlvddi,u)
	call shy_write_output_record(id,dtime,ivar,nel,1,nlv,nlvddi,v)

	end

!****************************************************************

