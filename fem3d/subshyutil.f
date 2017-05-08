!
! utility routines for shy file
!
! revision log :
!
! 15.10.2015    ggu     started routine
! 26.05.2016    ggu     new routines for opening and writing scalar file
! 02.02.2017    ggu     new routine shy_print_descriptions()
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

	call estimate_ngr(ngr)

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

	integer iaux,nkn,nel,nl
	real haux(1)
	integer, allocatable :: ile(:),ilk(:)

	call shy_get_params(id,nkn,nel,iaux,nl,iaux)

	if( nl > 1 ) then
	  call shy_set_layers(id,hlv)
	  call shy_set_layerindex(id,ilhv,ilhkv)
	else		!2d
	  haux(1) = -1.
	  haux(1) = 10000.
	  allocate(ile(nel),ilk(nkn))
	  ile = 1.
	  ilk = 1.
	  call shy_set_layers(id,haux)
	  call shy_set_layerindex(id,ile,ilk)
	  deallocate(ile,ilk)
	end if

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

	call shy_close(id)
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

	call shy_close(id)
	bok = .true.

	end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine shyfem_init_scalar_file(type,nvar,b2d,id)

        use levels

        implicit none

        character*(*) type      !type of file, e.g., hydro, ts, wave
        integer nvar		!total number of scalars to be written
        logical b2d		!2d fields
        integer id		!id for file (return)

        integer ftype,npr,nl
        character*80 file,ext,aux

        aux = adjustl(type)
        ext = '.' // trim(aux) // '.shy'        !no blanks in ext
        ftype = 2
        npr = 1
        nl = nlv
        if( b2d ) nl = 1

        call shy_make_output_name(trim(ext),file)
        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_simul_params(id)
        call shy_make_header(id)

        end

!****************************************************************

	subroutine shy_make_output_name(ext,file)

	character*(*) ext	!extension (with dot), example: .ext
	character*(*) file	!complete file name (return)

	character*80 simul

	call getfnm('runnam',simul)

	file = trim(simul) // trim(ext)

	end

!****************************************************************

	subroutine shy_set_simul_params(id)

	use shyfile

	implicit none

	integer id

	include 'simul.h'

	integer date,time
	character*80 title
	character*80 femver
	double precision dgetpar

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))
        title = descrp
        call get_shyfem_version(femver)

        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

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
c initialize data structure
c-----------------------------------------------------

	call shy_set_params(id,nkn,nel,npr,nlv,nvar)
        call shy_set_ftype(id,ftype)

	call shy_alloc_arrays(id)
	call shy_copy_basin_to_shy(id)
	call shy_copy_levels_to_shy(id)

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	end

!****************************************************************

	subroutine shy_make_header(id)

	use shyfile

	implicit none

	integer id

	integer ierr
	character*80 file

c-----------------------------------------------------
c write header of file
c-----------------------------------------------------

	call shy_write_header(id,ierr)

c-----------------------------------------------------
c error check
c-----------------------------------------------------

	if( ierr /= 0 ) then
	  write(6,*) 'error writing header of file ',ierr
	  write(6,*) 'id: ',id
	  call shy_get_filename(id,file)
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_make_header: writing header'
	end if

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

	end

!****************************************************************

	subroutine shy_get_string_descriptions(id,nvar,ivars,strings)

	use shyfile

	implicit none

	integer id
	integer nvar
	integer ivars(nvar)
	character*(*) strings(nvar)

	integer irec,nrec,ierr,i,isub
	integer ftype
	integer ivar,n,m,lmax
	double precision dtime

	call shy_get_ftype(id,ftype)

	if( ftype == 1 ) then		!hydro
	  ivars(1) = 1
	  ivars(2) = 1
	  ivars(3) = 3
	  ivars(4) = 3
	  strings(1) = 'water level'
	  strings(2) = 'water level (elemental)'
	  strings(3) = 'transport (velocity) x'
	  strings(4) = 'transport (velocity) y'
	else
	  irec = 0
	  nrec = 0
	  do
	    call shy_skip_record(id,dtime,ivar,n,m,lmax,ierr)
	    if( ierr /= 0 ) goto 99
	    nrec = nrec + 1
	    if( ivar < 0 ) cycle
	    irec = irec + 1
	    ivars(irec) = ivar
	    call ivar2string(ivar,strings(irec),isub)
	    if( irec == nvar ) exit
	  end do
	  do i=1,nrec
	    call shy_back_record(id,ierr)
	  end do
	end if

	return
   99	continue
	write(6,*) irec,nrec,nvar,ierr
	if( nrec == 0 ) write(6,*) 'no valid records in file'
	stop 'error stop shy_get_string_descriptions: reading record'
	end

!****************************************************************

	subroutine shy_print_descriptions(nvar,ivars,strings)

	implicit none

	integer nvar
	integer ivars(nvar)
	character*(*) strings(nvar)

	integer iv,ivar
	character*6 aux

        write(6,*) 'available variables: '
        write(6,*) 'total number of variables: ',nvar
        write(6,*) '   varnum     varid      varname'

        do iv=1,nvar
          ivar = ivars(iv)
	  aux = ' '
	  if( ivar == 3 ) aux = ' (2)  '
          write(6,'(2i10,a,a)') iv,ivar,aux,trim(strings(iv))
        end do

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
	  write(6,*) 'id: ',id
	  write(6,*) 'file name: ',trim(file)
	  write(6,*) dtime
	  write(6,*) ivar,n,m,nlv,nlvdi
	  call shy_info(id)
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
!****************************************************************
!****************************************************************

        subroutine shyfem_write_scalar(id,type,dtime,nvar,ivar,nlvddi,c)

! unconditionally writes to file (first call id must be 0)

	implicit none

	integer id
	character*(*) type
	double precision dtime
	integer nvar,ivar
	integer nlvddi
	real c(nlvddi,*)

	logical b2d

	if( id == 0 ) then
	  b2d = nlvddi == 1
          call shyfem_init_scalar_file(type,nvar,b2d,id)
	end if

	call shy_write_scalar_record(id,dtime,ivar,nlvddi,c)

	end

!****************************************************************
!****************************************************************
!****************************************************************

