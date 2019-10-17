
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! utility routines for shy file
!
! revision log :
!
! 15.10.2015	ggu	started routine
! 28.04.2016	ggu	changed VERS_7_5_9
! 25.05.2016	ggu	changed VERS_7_5_10
! 26.05.2016	ggu	new routines for opening and writing scalar file
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 27.06.2016	ggu	changed VERS_7_5_16
! 09.09.2016	ggu	changed VERS_7_5_17
! 02.02.2017	ggu	new routine shy_print_descriptions()
! 13.02.2017	ggu	changed VERS_7_5_23
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 11.07.2017	ggu	changed VERS_7_5_30
! 26.09.2017	ggu	limit written layers to min(nlv,nlvdi)
! 10.04.2018	ggu	prepare for mpi version (collect array and write)
! 11.04.2018	ggu	bug fix and hydro write
! 19.04.2018	ggu	changed VERS_7_5_45
! 11.05.2018	ggu	bug fix and hydro init, use global layer number
! 24.05.2018	ccf	bug fix exchanging nlvdi with nlv ($BUGNLV)
! 03.07.2018	ggu	in shy_write_output_record() handle 2d arrays
! 27.09.2018	ggu	new routines for writing constant layer structure
! 16.10.2018	ggu	changed VERS_7_5_50
! 25.10.2018	ggu	changed VERS_7_5_51
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 12.07.2019	ggu	some comments for shy_write_scalar_record()
! 27.09.2019	ggu	some lines commented for hydro output (HYD1)
! 17.10.2019	ggu	routines to check number of variables written (nvar_act)
!
! notes :
!
! open scalar file with shyfem_init_scalar_file()
! write scalar records with shy_write_scalar_record()
!
! for constant layer structure (sediments, etc.) please use:
!	shyfem_init_scalar_file_hlv()
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
	use shympi

	implicit none

	integer id

	integer nk,ne
	integer np,nl,nvar
	integer, allocatable :: nen_global(:,:)
	integer, allocatable :: in(:),ie(:),nen3(:,:)
	real, allocatable :: xg(:),yg(:),hm3(:,:)

	call shy_get_params(id,nk,ne,np,nl,nvar)

	allocate(nen3(3,ne))
	allocate(hm3(3,ne))
	allocate(xg(nk),yg(nk))
	allocate(ie(ne))
	allocate(in(nk))

	allocate(nen_global(3,nel))
	call adjust_element_index(nen_global)

	call shympi_exchange_array(nen_global,nen3)
	call shy_set_elemindex(id,nen3)

	call shympi_exchange_array(xgv,xg)
	call shympi_exchange_array(ygv,yg)
	call shy_set_coords(id,xg,yg)

	call shympi_exchange_array(hm3v,hm3)
	call shy_set_depth(id,hm3)

	call shympi_exchange_array(ipev,ie)
	call shympi_exchange_array(ipv,in)
	call shy_set_extnumbers(id,ie,in)

	call shympi_exchange_array(iarv,ie)
	call shympi_exchange_array(iarnv,in)
	call shy_set_areacode(id,ie,in)

	end

!****************************************************************

	subroutine adjust_element_index(nen_global)

! we have to adjust the element index with info from the other domains

	use basin
	use shympi

	implicit none

	integer nen_global(3,nel)

	integer id,ie,ii,k,ibase
	integer iint(nkn)

	do k=1,nkn
	  iint(k) = k
	end do

	call shympi_exchange_2d_node(iint)

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    id = id_node(k)	!this is the domain the node belongs to
	    ibase = nkn_cum_domains(id)
	    nen_global(ii,ie) = ibase + iint(k)
	  end do
	end do

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

! copies level structure as specified in str to shy file

	use levels
	use shyfile
	use shympi

	implicit none

	integer id

	integer nk,ne,np,nl,nvar

	real haux(1)
	integer, allocatable :: ile(:),ilk(:)

	call shy_get_params(id,nk,ne,np,nl,nvar)
	allocate(ile(ne),ilk(nk))	!is global size

	if( nl > 1 ) then
	  call shy_set_layers(id,hlv_global)
	  call shympi_exchange_array(ilhkv,ilk)
	  call shympi_exchange_array(ilhv,ile)
	  call shy_set_layerindex(id,ile,ilk)
	else		!2d
	  haux(1) = 10000.
	  ile = 1
	  ilk = 1
	  call shy_set_layers(id,haux)
	  call shy_set_layerindex(id,ile,ilk)
	end if

	deallocate(ile,ilk)

	end

!****************************************************************

	subroutine shy_set_levels_in_shy(id,nl0,hlv0)

! sets constant level

	use levels
	use shyfile
	use shympi

	implicit none

	integer id
	integer nl0
	real hlv0(nl0)

	integer nk,ne,np,nl,nvar

	integer, allocatable :: ile(:),ilk(:)

	call shy_get_params(id,nk,ne,np,nl,nvar)
	allocate(ile(ne),ilk(nk))	!is global size

	if( nl /= nl0 ) then
	  write(6,*) 'nl,nl0: ',nl,nl0
	  stop 'error stop shy_set_levels_in_shy: nl /= nl0'
	end if

	ile = nl
	ilk = nl
	call shy_set_layers(id,hlv0)
	call shy_set_layerindex(id,ile,ilk)

	deallocate(ile,ilk)

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

        subroutine shyfem_init_hydro_file(type,b2d,id)

        use levels
        use shympi

        implicit none

        character*(*) type      !type of file, e.g., hydro, ts, wave
        logical b2d		!2d fields
        integer id		!id for file (return)

        integer ftype,npr,nl
        integer nvar		!total number of scalars to be written
        character*80 file,ext,aux

        aux = adjustl(type)
        ext = '.' // trim(aux) // '.shy'        !no blanks in ext
        ftype = 1
        npr = 3
	nvar = 4
        nl = nlv_global
        if( b2d ) nl = 1

        call shy_make_output_name(trim(ext),file)
        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_simul_params(id)
        call shy_make_header(id)

        end

!****************************************************************

        subroutine shyfem_init_scalar_file(type,nvar,b2d,id)

! initializes scalar file with layer structure from str

        use levels
        use shympi

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
        nl = nlv_global
        if( b2d ) nl = 1

        call shy_make_output_name(trim(ext),file)
        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_simul_params(id)
        call shy_make_header(id)

        end

!****************************************************************
!write shy header
!copied from shyfem_init_scalar_file
!nvar not used

        subroutine shyfem_init_lgr_file(type,nvar,b2d,id)

        use levels
        use shympi

        implicit none

        character*(*) type      !type of file, e.g., hydro, ts, wave
        integer nvar		!total number of scalars to be written
        logical b2d		!2d fields
        integer id		!id for file (return)

        integer ftype,npr,nl
        character*80 file,ext,aux

        aux = adjustl(type)
        ext = '.' // trim(aux) // '.shy'        !no blanks in ext
        ftype = 3
        npr = 1
        nl = nlv_global
        if( b2d ) nl = 1

        call shy_make_output_name(trim(ext),file)
        call shy_open_output_file(file,npr,nl,nvar,ftype,id)
        call shy_set_simul_params(id)
        call shy_make_header(id)

        end

!****************************************************************

        subroutine shyfem_init_scalar_file_hlv(type,nvar,nl0,hlv0,id)

! initializes scalar file with constant custom layer structure

        use levels
        use shympi

        implicit none

        character*(*) type      !type of file, e.g., hydro, ts, wave
        integer nvar		!total number of scalars to be written
	integer nl0		!total number of (constant) layers
	real hlv0(nl0)		!layer structure
        integer id		!id for file (return)

        integer ftype,npr
        character*80 file,ext,aux

        aux = adjustl(type)
        ext = '.' // trim(aux) // '.shy'        !no blanks in ext
        ftype = 2
        npr = 1

        call shy_make_output_name(trim(ext),file)
        call shy_open_output_file(file,npr,nl0,nvar,ftype,id)
	call shy_set_levels_in_shy(id,nl0,hlv0)
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

	if( id <= 0 ) return

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))
        title = descrp
        call get_shyfem_version_and_commit(femver)

        call shy_set_date(id,date,time)
        call shy_set_title(id,title)
        call shy_set_femver(id,femver)

	end

!****************************************************************

	subroutine shy_close_output_file(id)

	use shyfile
	use shympi

	implicit none

	integer id

	if( shympi_is_master() ) call shy_close(id)

	end

!****************************************************************

	subroutine shy_open_output_file(file,npr,nl,nvar,ftype,id)

	use basin
	use shyfile
	use shympi

	implicit none

	character*(*) file	!file name
	integer npr		!max number of parameters per node
	integer nl		!vertical dimension
	integer nvar		!number of variables
	integer ftype		!type of file: 1=ous 2=nos
	integer id		!id of opened file (return)

	logical bopen

c-----------------------------------------------------
c open file
c-----------------------------------------------------

	bopen = shympi_is_master()	!opens file only if master

	id = shy_init(file,bopen)

	if( id == 0 ) then
	  write(6,*) 'error opening file'
	  write(6,*) 'file name: ',trim(file)
	  stop 'error stop shy_open_output_file: opening file'
	end if

c-----------------------------------------------------
c initialize data structure
c-----------------------------------------------------

	call shy_set_params(id,nkn_global,nel_global,npr,nl,nvar)
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

	if( id <= 0 ) return

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

	ivars = 0
	strings = ' '

	if( id <= 0 ) return

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
	  call shy_back_records(id,nrec,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	return
   97	continue
	stop 'error stop shy_check_nvar: backspacing'
   99	continue
	write(6,*) irec,nrec,nvar,ierr
	if( nrec == 0 ) write(6,*) 'no valid records in file'
	stop 'error stop shy_get_string_descriptions: reading record'
	end

!****************************************************************

	subroutine shy_check_nvar(id,nvar)

	use shyfile

	implicit none

	integer id
	integer nvar

	integer irec,nrec,ierr,i,isub
	integer ftype
	integer ivar,n,m,lmax,ivar_first
	integer nkn,nel,npr,nlv
	double precision dtime
	character*80 string

	if( id <= 0 ) return

	call shy_get_ftype(id,ftype)

	if( ftype == 1 ) then		!hydro
	  if( nvar /= 4 ) then
	    write(6,*) 'nvar = ',nvar
	    write(6,*) 'nvar incompatible with ftype==1'
	    stop 'error stop shy_check_nvar: nvar/=4'
	  end if
	  return
	end if

	irec = 0	!records with data
	nrec = 0	!records read (also ivar<0)
	ivar_first = -999

	do
	  call shy_skip_record(id,dtime,ivar,n,m,lmax,ierr)
	  if( ierr /= 0 ) exit
	  nrec = nrec + 1
	  if( ivar == ivar_first ) exit
	  if( ivar_first == -999 ) ivar_first = ivar
	  if( ivar < 0 ) cycle
	  irec = irec + 1
	end do

	if( ierr /= 0 ) then
	  call shy_back_one(id,ierr)	!this skips over EOF
	  if( ierr /= 0 ) goto 97
	end if

	call shy_back_records(id,nrec,ierr)
	if( ierr /= 0 ) goto 97
	if( irec == nvar ) return

!	here error management

	write(6,*) '*** extra variables found in file...'
	write(6,*) '*** nvar declared: ',nvar
	write(6,*) '*** nvar found:    ',irec
	write(6,*) '*** resetting nvar to ',irec

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	nvar = irec
	call shy_set_params(id,nkn,nel,npr,nlv,nvar)

	return
   97	continue
	stop 'error stop shy_check_nvar: backspacing'
   99	continue
	write(6,*) irec,nrec,nvar,ierr
	stop 'error stop shy_check_nvar: internal error (1)'
	end

!****************************************************************

	subroutine shy_print_descriptions(nvar,ivars,strings)

	implicit none

	integer nvar
	integer ivars(nvar)
	character*(*) strings(nvar)

	integer iv,ivar
	character*6 aux

        write(6,*) 'total number of available variables: ',nvar
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

	subroutine shy_write_output_record(id,dtime,ivar,n,m,lmax
     +					,nlvdi,c)

	use shyfile
	use shympi

	implicit none

	integer id
	double precision dtime
	integer ivar
	integer n,m,lmax,nlvdi		!n is local value
	real c(nlvdi,m*n)

	integer ierr,nn,nl,i
	real, allocatable :: cl(:,:),cg(:,:)
	character*80 file

! nlvdi cannot be greater than lmax
! lmax is global nlv or 1 if we output is 2d
! nlvdi is vertical dimension of c
! in case of serial run we have lmax == 1 or lmax == nlvdi
! nlvdi == 1 in case of 2d simulation
! in case of domain decomposition we have lmax == nlvdi == 1 for 2d run
! or for 3d run and 2d output
! otherwise for real 3d run and output we have lmax > 1
! nlvdi < lmax is possible, because this domain has less
! levels than the global domain
!
! summary:
!	lmax == 1	2d run or 2d output, nlvdi == 1
!	lmax > 1	nlvdi <= lmax

	if( id <= 0 ) return

	if( n == nkn_local ) then
	  nn = nkn_global
	else if( n == nel_local ) then
	  nn = nel_global
	else
	  call shy_info(id)
	  write(6,*) 'nkn: ',n,nkn_local,nkn_global
	  write(6,*) 'nel: ',n,nel_local,nel_global
	  stop 'error stop shy_write_output_record: n mismatch'
	end if

        if( lmax > 1 .and. m > 1 ) then         !$BUGNLV
          call shy_info(id)
          stop 'error stop shy_write_output_record: nlvdi&m>1'
	end if

	if( lmax == 1 .and. nlvdi > 1 ) then		!$BUGNLV
	  !call shy_info(id)	!this section commented for hydro write $HYD1
	  !write(6,*) 'lmax = ',lmax,'  nlvdi = ',nlvdi
	  !stop 'error stop shy_write_output_record: nlvdi>lmax==1'
	else if( lmax > 1 .and. nlvdi > lmax ) then		!$BUGNLV
	  call shy_info(id)
	  write(6,*) 'lmax = ',lmax,'  nlvdi = ',nlvdi
	  stop 'error stop shy_write_output_record: nlvdi>lmax>1'
	end if

	nl = lmax
	if( m > 1 ) then
	  nl = 1
	  allocate(cl(m,n),cg(m,nn))
	  cl = reshape(c(1,:),(/m,n/))
	else if( lmax == 1 ) then
	  nl = 1
	  allocate(cl(nl,n),cg(nl,nn))
	  cl(1,:) = c(1,:)
	else
	  nl = nlvdi
	  allocate(cl(nl,n),cg(nl,nn))
	  cl = c
	end if

	cg = 0.
	call shympi_exchange_array(cl,cg)
	call shy_write_record(id,dtime,ivar,nn,m,lmax,nl,cg,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'error writing output file ',ierr
	  call shy_get_filename(id,file)
	  write(6,*) 'id: ',id
	  write(6,*) 'file name: ',trim(file)
	  write(6,*) dtime
	  write(6,*) ivar,n,m,lmax,nlvdi
	  call shy_info(id)
	  stop 'error stop shy_write_output_file: writing record'
	end if

	end

!****************************************************************

	subroutine shy_write_scalar_record(id,dtime,ivar,nlvddi,c)

	use basin
	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	integer nlvddi
	real c(nlvddi,nkn)

	integer iaux,nlv

	if( id <= 0 ) return

	call shy_get_params(id,iaux,iaux,iaux,nlv,iaux)	!nlv is global here

! here it can happen that nlvdi == 1 because we want 2d output
! and nlv > 1 because it is a 3d run
! this situation cannot be distinguished from the situation where
! a 3d MPI run has one domain with only one layer (nlvdi==1) and a 3d output
! therefore it is important to use shy_write_scalar_record2d() for 2d output

	call shy_write_output_record(id,dtime,ivar,nkn,1,nlv,nlvddi,c)

	end

!****************************************************************

	subroutine shy_write_scalar_record2d(id,dtime,ivar,c)

	use basin
	use shyfile

	implicit none

	integer id
	double precision dtime
	integer ivar
	real c(nkn)

	integer iaux

	if( id <= 0 ) return

	call shy_write_output_record(id,dtime,ivar,nkn,1,1,1,c)

	end

!****************************************************************

	subroutine shy_write_hydro_records(id,dtime,nlvddi,z,ze,u,v)

	use basin
	use shyfile
	use shympi

	implicit none

	integer id
	double precision dtime
	integer nlvddi
	real z(nkn)
	real ze(3,nel)
	real u(nlvddi,nel)
	real v(nlvddi,nel)

	integer ivar,nk,ne,npr,nlv,iaux

	if( id <= 0 ) return

	call shy_get_params(id,iaux,iaux,iaux,nlv,iaux)	!nlv is global here
	!nlv = min(nlv,nlvddi)

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

	subroutine shy_reset_nvar_act(id)

	use shyfile

	implicit none

	integer id

	pentry(id)%nvar_act = 0

	end

!****************************************************************

	subroutine shy_get_nvar_act(id,nvar_act)

	use shyfile

	implicit none

	integer id
	integer nvar_act

	nvar_act = pentry(id)%nvar_act

	end

!****************************************************************
! next two debug routines to be deleted later
!****************************************************************

	subroutine write_debug_vel2(text,ng,uv)

	use basin
	use shympi

	implicit none

	character*(*) text
	integer ng
	real uv(ng)

	integer i,ie,ip,ntot
	real ul(ng/16)
	real ug(ng/16)
	integer ipe(ng)

	if( ng /= nel_global ) return

	call shympi_exchange_array(ipev,ipe)

	ntot = ng/16
	write(6,*) ng,ng/16,ntot

	ug = 0.
	do i=1,ng
	  ie = ipe(i)
	  if( mod(ie,16) == 0 ) then
	    ip = ie/16
	if( ip > ntot ) then
	  write(6,*) i,ie,ip
	  stop
	end if
	    ug(ip) = uv(i)
	  end if
	end do

	write(my_unit,*) 'global field: ',trim(text),ng
	do i=1,ng/16
	  write(my_unit,*) i,ug(i)
	end do

	end

!****************************************************************

	subroutine write_debug_vel(text,nl,ng,u,uv)

	use basin
	use shympi

	implicit none

	character*(*) text
	integer nl,ng
	real u(nl)
	real uv(ng)

	integer i,ie,ip
	real ul(ng/16)
	real ug(ng/16)
	integer ipe(ng)

	call shympi_exchange_array(ipev,ipe)

	ul = 0.
	do i=1,nl
	  ie = ipev(i)
	  if( mod(ie,16) == 0 ) then
	    ip = ie/16
	    ul(ip) = u(i)
	  end if
	end do

	ug = 0.
	do i=1,ng
	  ie = ipe(i)
	  if( mod(ie,16) == 0 ) then
	    ip = ie/16
	    ug(ip) = uv(i)
	  end if
	end do

	write(my_unit,*) 'new field: ',trim(text),nl,ng
	do i=1,ng/16
	  write(my_unit,*) i,ul(i),ug(i)
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine shy_write_scalar(id,type,dtime,nvar,ivar,nlvddi,c)

! unconditionally writes to file (first call id must be 0)

	implicit none

	integer id
	character*(*) type
	double precision dtime
	integer nvar,ivar
	integer nlvddi
	real c(nlvddi,*)

	logical b2d

	if( id < 0 ) return

	if( id == 0 ) then
	  b2d = ( nlvddi == 1 )
          call shyfem_init_scalar_file(type,nvar,b2d,id)
	end if

	call shy_write_scalar_record(id,dtime,ivar,nlvddi,c)

	end

!****************************************************************
!****************************************************************
!****************************************************************

