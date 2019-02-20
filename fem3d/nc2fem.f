
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

! convert nc files to fem files
!
! contents :
!
!
! revision log :
!
! 03.07.2018    ggu     revision control introduced
! 04.07.2018    ggu     single points introduced
! 06.07.2018    ggu     bug fix in handle_data: valnew was not 3d
!
! notes :
!
! still to do:
!
!	handle iregular grid
!	check 3d variables
!	check resampling regular grid
!	pass depth values for hd if 3d
!	pass hlv for 3d
!	construct ilhkv
!
!*********************************************************************

	subroutine write_about

	implicit none

	write(6,*) 'converts nc (netcdf) file to fem file'
	write(6,*) 
	write(6,*) 'The file created is either a regular fem file'
	write(6,*) 'or it can be single points to be used for'
	write(6,*) 'boundary conditions given with the option -single'
	write(6,*) 'The domain can be adjusted with -domain'
	write(6,*) 'The variables to be written are given with -vars'
	write(6,*) 
	write(6,*) 'The program should recognize most of the'
	write(6,*) 'dimensions and coordinates used in nc files'
	write(6,*) 'If some of these are not recognized you can'
	write(6,*) 'insert them at the end of file nc_dim_coords.f'
	write(6,*) 'and then recompile with "make fem"'
	write(6,*) 'The same is true for the description of the'
	write(6,*) 'variables written to file'

	end

!*********************************************************************

	program nc2fem

	use clo

	implicit none

	integer nxdim,nydim,nlvdim,ndim

	integer ncid
        character*132 file
        character*80 var_name,files,sfile
        character*80 name,xcoord,ycoord,zcoord,tcoord,bathy,slmask
        character*80 varline,descrpline,factline,text,fulltext,dstring
        character*80, allocatable :: vars(:)
        character*80, allocatable :: descrps(:)
        character*80, allocatable :: sfacts(:)
        real, allocatable :: facts(:)		!factor for multiplication
        real, allocatable :: offs(:)		!offset to add
        real, allocatable :: flags(:)		!flag for no data
        integer ndims, nvars, ngatts, unlim
	integer dim_id,dim_len
	integer nt,nx,ny,nz,nz1
	integer nxnew,nynew,ns
	integer nit,i,it,n,nd,nrec
	integer iwhat
	integer dims(10)
	integer year0,date0,time0
	integer iftype
	integer ifreq
	integer regexpand
	real rfact
	real regpar_data(9)
	real regpar(9)
	real, allocatable :: xlon(:,:)
	real, allocatable :: ylat(:,:)
	real, allocatable :: zdep(:)			!mid-layer depth
	real, allocatable :: zzdep(:)			!bottom depth
	real, allocatable :: hlv(:)			!bottom depth
	real, allocatable :: bat(:,:)
	real, allocatable :: slm(:,:)
	real, allocatable :: batnew(:,:)
	double precision t
	logical bverb,bcoords,btime,binfo,bvars,bwrite,bdebug,bsilent
	logical binvertdepth,binvertslm,bunform,bquiet,blist
	logical bregular,bsingle,babout
	logical exists_var

	interface
	  subroutine parse_strings(line,n,vars)
	  character*(*) line
	  integer n
	  character(len=80), allocatable :: vars(:)
	  end subroutine parse_strings
	end interface

c-----------------------------------------------------------------
c initialize parameters
c-----------------------------------------------------------------

	file = ' '

c-----------------------------------------------------------------
c find out what to do
c-----------------------------------------------------------------

	call clo_init('nc2fem','nc-file','1.2')
        call clo_add_info('converts nc (netcdf) file to fem file')

	call clo_add_sep('general options')

        call clo_add_option('about',.false.,'about this program')
        call clo_add_option('info',.false.,'general info on nc file')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('quiet',.false.,'be as quiet as possible')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('debug',.false.,'produce debug information')
        call clo_add_option('varinfo',.false.
     +		,'list variables contained in file')
        call clo_add_option('list',.false.
     +		,'list possible names for description')

	call clo_add_sep('special variables')

        call clo_add_option('time',.false.
     +		,'write available time records to terminal')
        call clo_add_option('coords',.false.,'write coordinate file')
        call clo_add_option('bathy var',' '
     +		,'write bathymetry file using variable var')
        call clo_add_option('slmask var',' '
     +		,'write sea-land mask file using variable var')

        call clo_add_option('invertdepth',.false.
     +		,'invert depth values for bathymetry')
        call clo_add_option('invertslm',.false.
     +		,'invert slmask values (0 for sea)')
        call clo_add_option('unform',.false.
     +		,'write fem file unformatted')

	call clo_add_sep('output general variables')

        call clo_add_option('vars text',' '
     +		,'write variables given in text to out.fem')
        call clo_add_option('descrp text',' '
     +		,'use this description for variables')
        call clo_add_option('fact facts',' '
     +		,'scale vars with these factors')
        call clo_add_option('single file',' '
     +		,'file containing x/y coordinates for interpolation')

        call clo_add_option('domain limits',' '
     +		,'give domain limits and resolution')
        call clo_add_option('regexpand iexp',-1,'expand regular grid')

        call clo_add_com('    iexp>0 expands iexp cells, =0 whole grid')


        call clo_add_sep('additional information')
        call clo_add_com('  var is name of variable in nc file')
        call clo_add_com('  facts is list of factors for'
     +				// ' multiplication of vars')
        call clo_add_com('  text is list of variables and descriptions'
     +				// ' for output')
        call clo_add_com('    seperate with comma and leave no space')
        call clo_add_com('    or enclose in ""')
        call clo_add_com('    example: U10,V10,PSFC or "U10 V10 PSFC"')
        call clo_add_com('  limits gives extension of domain')
        call clo_add_com('    for regular domain use x0,y0,x1,y1 ' //
     +				'(dx,dy are the same)')
	call clo_add_com('    for irregular domain use ' //
     +				'dx[,dy[,x0,y0,x1,y1]]')

	call clo_parse_options

	call clo_get_option('about',babout)
	call clo_get_option('info',binfo)
	call clo_get_option('verbose',bverb)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('silent',bsilent)
	call clo_get_option('debug',bdebug)
	call clo_get_option('varinfo',bvars)
	call clo_get_option('list',blist)
	call clo_get_option('time',btime)
	call clo_get_option('coords',bcoords)
	call clo_get_option('bathy',bathy)
	call clo_get_option('slmask',slmask)
	call clo_get_option('vars',varline)
	call clo_get_option('descrp',descrpline)
	call clo_get_option('single',sfile)
	call clo_get_option('domain',dstring)
	call clo_get_option('regexpand',regexpand)
	call clo_get_option('fact',factline)

	call clo_get_option('invertdepth',binvertdepth)
	call clo_get_option('invertslm',binvertslm)
	call clo_get_option('unform',bunform)

	if( babout ) then
	  call write_about
	  call exit(99)
	end if

	if( blist ) then
	  call list_strings
	  call exit(99)
	end if

	call nc_init
	call clo_check_files(1)
        call clo_get_file(1,file)
        if( file == ' ' ) call clo_usage

	bwrite = bverb .or. binfo
	if( bsilent ) bquiet = .true.
	bsingle = ( sfile /= ' ' ) 

        !call read_frequency(ifreq)

c-----------------------------------------------------------------
c open nc file and write info
c-----------------------------------------------------------------

	call nc_open_read(ncid,file)
	if( .not. bsilent ) call global_information(ncid)

	call ncnames_init

        call get_dims_and_coords(ncid,bwrite
     +                  ,nt,nx,ny,nz
     +                  ,tcoord,xcoord,ycoord,zcoord)

	nxdim = max(1,nx)
	nydim = max(1,ny)
	nlvdim = max(1,nz)
	ndim = nxdim*nydim*nlvdim

	allocate(xlon(nxdim,nydim),ylat(nxdim,nydim))
	allocate(zdep(nlvdim),zzdep(nlvdim),hlv(nlvdim))
	allocate(bat(nxdim,nydim),slm(nxdim,nydim))

	call nc_set_domain(1,nxdim,1,nydim,1,nlvdim)

	if( bvars ) then
	  !call nc_dims_info(ncid)
	  call nc_vars_info(ncid,bverb)
	end if

	call setup_nc_time(ncid,bverb)
	if( bwrite ) call print_minmax_time_records(ncid)

	!if( binfo .or. bvars ) call exit(99)

	if( blist ) then
	  call list_strings
	  stop
	end if

c-----------------------------------------------------------------
c handle time
c-----------------------------------------------------------------

	if( btime ) then
	  call print_time_records(ncid)
	end if

c-----------------------------------------------------------------
c handle coordinates and special variables
c-----------------------------------------------------------------

	call setup_coordinates(ncid,bverb,xcoord,ycoord
     +				,nxdim,nydim,nx,ny
     +				,xlon,ylat)
	call setup_zcoord(ncid,bverb,zcoord,nlvdim,nz,zdep,nz1,hlv)
	call setup_bathymetry(ncid,bverb,binvertdepth,bathy
     +				,nxdim,nydim,nx,ny,bat)
	call setup_sealand(ncid,bverb,slmask,nxdim,nydim,nx,ny,slm)
	if( binvertslm ) slm = 1.-slm

c-----------------------------------------------------------------
c check regularity of grid 
c - in regpar_data is info on original grid
c - in regpar is info on desired regular output grid
c - bregular is true if regular_data is regular
c-----------------------------------------------------------------

	!call handle_unusual_coordinates(nxdim,nydim,xlon,ylat)
	call check_regular_coords(nxdim,nydim,xlon,ylat
     +				,bregular,regpar_data)
	if( bverb ) write(6,*) bregular,regpar_data
	call handle_domain(bverb,dstring,bregular,regpar_data,regpar)

	if( .not. bsilent ) then
	 if( bregular ) then
	   write(6,*) 'original coordinates are regular'
	   if( bverb ) write(6,*) regpar
	 else
	   write(6,*) 'original coordinates are irregular'
	   if( bverb ) write(6,*) regpar
	 end if
	end if

c-----------------------------------------------------------------
c write special files or info
c-----------------------------------------------------------------

	if( bcoords ) then
	  call write_coordinates(nxdim,nydim,xlon,ylat)
	  call write_regular_coordinates_as_grd(regpar)
	  if( zcoord .ne. ' ' ) call write_zcoord(nz,zdep,nz1,hlv)
	end if

	if( bathy .ne. ' ' ) then
	  call write_2d_dat('bathy.dat',nxdim,nydim,bat)
	  call write_2d_grd('bathy.grd',nxdim,nydim,xlon,ylat,bat)
	  files = 'bathy.dat bathy.grd'
	  if( bregular ) then
	    call write_2d_fem('bathy.fem','bathymetry',regpar_data,bat)
	    files = trim(files) // ' bathy.fem'
	  end if
	  write(6,*) 'bathymetry written to files: ',files
	end if

	if( slmask .ne. ' ' ) then
	  call write_2d_dat('sea_land.dat',nxdim,nydim,slm)
	  call write_2d_grd('sea_land.grd',nxdim,nydim,xlon,ylat,slm)
	  files = 'sea_land.dat sea_land.grd'
	  if( bregular ) then
	    call write_2d_fem('sea_land.fem','index',regpar_data,slm)
	    files = trim(files) // ' sea_land.fem'
	  end if
	  write(6,*) 'sea-land mask written to files: ',files
	end if

c-----------------------------------------------------------------
c check variables to write
c-----------------------------------------------------------------

	n = -1
	call parse_strings(varline,n,vars)
	nd = n

	call parse_strings(descrpline,nd,descrps)
	call parse_strings(factline,nd,sfacts)
	allocate(facts(nd),offs(nd),flags(nd))
	call setup_facts(nd,sfacts,facts)

	call handle_variable_description(ncid,n,vars,descrps,.not.bquiet)

c-----------------------------------------------------------------
c set up interpolation
c-----------------------------------------------------------------

	nxnew = nint(regpar(1))
	nynew = nint(regpar(2))

	if( bsingle ) then		!interpolate on single points (BC)
	  call prepare_single(sfile,ns,nxdim,nydim,xlon,ylat,regpar)
	  nxnew = ns
	  nynew = 1
	else if( bregular ) then	!no interpolation - already regular
	  call prepare_no_interpol
	else				!interpolate onto regular grid
	  call prepare_interpol(nxdim,nydim,xlon,ylat,regpar)
	end if

c-----------------------------------------------------------------
c interpolate special variables
c-----------------------------------------------------------------

	if( bathy .ne. ' ' ) then
	  allocate(batnew(nxdim,nydim))
	  batnew = -999.
	  call handle_interpol_2d(nxdim,nydim,bat,nxnew,nynew,batnew)
	  call write_2d_fem('bathy_new.fem','bathymetry',regpar,batnew)
	  call write_2d_grd_regular('bathy_new.grd',regpar,batnew)
	end if

c-----------------------------------------------------------------
c check if some dimension is inverted
c-----------------------------------------------------------------

	if( bregular ) then
	  call check_invert(regpar)
	end if

c-----------------------------------------------------------------
c write variables
c-----------------------------------------------------------------

	call write_variables(ncid,n,bunform,bdebug,bquiet
     +				,vars,descrps,facts,offs,flags
     +				,regexpand
     +				,nxdim,nydim,nlvdim
     +				,xlon,ylat
     +				,nz1,hlv
     +				,nxnew,nynew,regpar,nrec)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	if( nrec > 0 .and. .not. bsilent ) then
	  write(6,*) 'total number of time records written: ',nrec
	  write(6,*) 'variables written: ',n
	  write(6,*) 'output written to file out.fem'
	  call print_description(nd,vars,descrps)	!for information
	  write(6,*) 'facts:      ',facts
	  write(6,*) 'offsets:    ',offs
	  write(6,*) 'fill_value: ',flags
	end if

	if( .not. bsilent ) then
          write(6,*) 'Successful completion of routine nc2fem'
	end if
	call exit(99)

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine handle_variable_description(ncid,nvar
     +				,vars,descrps,bwrite)

	use ncnames
	use shyfem_strings

	implicit none

	integer ncid,nvar
	character*(*) vars(nvar)
	character*(*) descrps(nvar)
	logical bwrite

	logical bstop
	integer var_id,i,ivar
	character*15 var,descrp,short
	character*50 full
	character*15 shorts(nvar)

	if( nvar <= 0 ) return

	do i=1,nvar
	  var = vars(i)
	  call nc_get_var_id(ncid,var,var_id)
	  if( var_id == 0 ) then
	    write(6,*) 'no such variable: ',trim(var)
	    stop 'error stop handle_variable_description: '//
     +			'no such variable'
	  end if
	  if( bwrite ) call nc_var_info(ncid,var_id,.true.)
	end do

	bstop = .false.
	shorts = ' '

	do i=1,nvar
	  var = vars(i)
	  short = descrps(i)
	  if( short == ' ' ) call ncnames_get_var(ncid,var,short)
	  full = ' '
	  if( short /= ' ' ) then
	    call strings_get_ivar(short,ivar)
	    call strings_get_full_name(ivar,full)
	  end if
	  if( short == ' ' ) bstop = .true.
	  shorts(i) = short
	  descrps(i) = full
	end do

	call handle_direction(nvar,descrps)

	if( bwrite .or. bstop ) then
	  call print_description(nvar,vars,descrps)
	end if

	if( bstop ) then
	  write(6,*) 'some variables have no description'
	  write(6,*) 'please provide the description through -descrp'
	  write(6,*) '-list gives a list of possible names'
	  write(6,*) 'you can also permanently add the variable'
	  write(6,*) 'description by inserting it into'
	  write(6,*) 'subroutine ncnames_add_variables()'
	  stop 'error stop: no description'
	end if

	end

c*****************************************************************

	subroutine handle_direction(nvar,descrps)

	use shyfem_strings

	implicit none

	integer nvar
        character*80 :: descrps(nvar)

	integer i,idir
	character*80 text
	character*1, save :: post(0:3) = (/' ','x','y','z'/)

	logical has_direction

	idir = 0

	do i=1,nvar
	  text = descrps(i)
	  if( has_direction(text) ) then
	    idir = mod(idir+1,4)
	    text = trim(text) // ' - ' // post(idir)
	  end if
	  descrps(i) = text
	end do

	end

c*****************************************************************

	subroutine print_description(nvar,vars,descrps)

	use shyfem_strings

	implicit none

	integer nvar
	character*(*) vars(nvar)
	character*(*) descrps(nvar)

	integer i
	character*15 var,short
	character*50 full

	write(6,*) 'variable description:'
	write(6,*) '  iv  varname        descrp         full name'
	do i=1,nvar
	  var = vars(i)
	  full = descrps(i)
	  call strings_get_short_name(full,short)
	  write(6,'(i5,4a)') i,'  ',var,short,trim(full)
	end do

	end

c*****************************************************************

	subroutine parse_strings(line,n,vars)

	implicit none

	character*(*) line
	integer n
	character(len=80), allocatable :: vars(:)

	integer norig,ndim

	norig = n
	call count_variables(line,n)

	if( norig == -1 ) norig = n
	allocate(vars(norig))
	vars = ' '

	if( n > norig ) then
	  write(6,*) 'too many strings given: ',n
	  write(6,*) 'expecting: ',norig
	  write(6,*) 'line: ',trim(line)
	  stop 'error stop parse_strings: too many strings'
	else if( n /= norig .and. n > 1 ) then
	  write(6,*) 'wrong number of strings given: ',n
	  write(6,*) 'possible numbers: ',0,1,ndim
	  write(6,*) 'line: ',trim(line)
	  stop 'error stop parse_strings: wrong number of strings'
	end if
	
	call parse_variables(line,n,vars)

	if( n == 1 ) vars = vars(1)	!set all values with only value given

	end

c*****************************************************************

	subroutine count_variables(varline,n)

! parses line of variables etc..
! as separators only commas "," are allowed
! ""		n = 0
! "a"		n = 1
! ","		n = 2
! "a,"		n = 2
! ",b"		n = 2
! "a,b"		n = 2
! "a,b,"	n = 3

	implicit none

	character*(*) varline
	integer n

	logical btoken
	integer i,istart
	character*1 c
	character(len=len(varline)+1) string

	! we have at least one blank in the string

	n = 1
	istart = 1
	btoken = .false.
	string = adjustl(varline)

	do i=1,len_trim(string)+1
	  c=string(i:i)
	  if( c /= ' ' .and. c /= ',' ) btoken = .true.
	  if( c == ',' ) then			!comma (separator)
	    istart = i + 1
	    n = n + 1
	  end if
	end do

	if( .not. btoken ) n = 0

	end

c*****************************************************************

	subroutine parse_variables(varline,ndim,vars)

! parses line of variables etc..
! as separators only commas "," are allowed
! ""		n = 0
! "a"		n = 1
! ","		n = 2
! "a,"		n = 2
! ",b"		n = 2
! "a,b"		n = 2
! "a,b,"	n = 3

	implicit none

	character*(*) varline
	integer ndim			!number of variables we have to read
	character*(*) vars(ndim)

	logical btoken
	integer i,istart,n
	character*1 c
	character(len=len(varline)+1) string

	! we have at least one blank in the string

	vars = ' '
	if( ndim == 0 ) return

	n = 1
	istart = 1
	btoken = .false.
	string = adjustl(varline)

	do i=1,len_trim(string)+1
	  c=string(i:i)
	  if( c /= ' ' .and. c /= ',' ) btoken = .true.
	  if( c == ',' ) then			!comma (separator)
	    vars(n) = string(istart:i-1)
	    istart = i + 1
	    n = n + 1
	    if( n > ndim ) then
	      stop 'error stop parse_variables: n>ndim'
	    end if
	  end if
	end do

	vars(n) = string(istart:i-1)
	if( .not. btoken ) n = 0

	if( n /= ndim ) then
	  write(6,*) 'ndim = ',ndim,'   n = ',n
	  write(6,*) 'line: ',trim(varline)
	  stop 'error stop parse_variables: ndim/=n'
	end if

	end

c*****************************************************************

	subroutine write_variables(ncid,nvar,bunform,bdebug,bquiet
     +					,vars,descrps,facts,offs,flags
     +					,regexpand
     +					,nx,ny,nz
     +					,x,y
     +					,nz1,hlv
     +					,nxnew,nynew,regpar,nrec)

	use iso8601

	implicit none

	integer ncid
	integer nvar
	logical bunform,bdebug,bquiet
	character*(*) vars(nvar)
	character*(*) descrps(nvar)
	real facts(nvar)
	real offs(nvar)
	real flags(nvar)
	integer regexpand
	integer nx,ny,nz		!size of data in nc file
	real x(nx,ny),y(nx,ny)
	integer nz1
	real hlv(nz1)
	integer nxnew,nynew		!size of regular grid
	real regpar(9)			!regular grid to which to interpolate
	integer nrec			!how many records written (return)

	logical bvert,bexpand
	integer level
	integer nit,it,var_id,i,ns
	integer iformat,nvers,ntype,ndd
	integer iunit,lmax,ierr,nzz,npnew
	integer datetime(2)
	integer ids(nvar)
	integer dims(nvar)
	real, save :: my_flag = -999.
	real data(nx,ny,nz)
	double precision atime,avalue,dtime
	character*20 line,stime
	character*80 atext,string,aname
	character(len=len(vars)) var

	logical nc_has_var_attrib,is_single

	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)
	real, allocatable :: femdata(:,:,:)
	real, allocatable :: data2d(:,:)

	integer ifileo

	!level = 1
	level = 0	!mixing 2d/3d - not yet implemented

	nrec = 0
	if( nvar == 0 ) return

	iformat = 1
	if( bunform ) iformat = 0
	dtime = 0.
	nvers = 0
	ntype = 11
	if( is_single() ) ntype = 1
	lmax = nz
	npnew = nxnew*nynew
	string = 'unknown'
	bexpand = ( regexpand > -1 )
	!hlv = 0.

	allocate(hd(npnew),ilhkv(npnew))
	allocate(femdata(nz,nxnew,nynew))
	hd = -999.
	ilhkv = lmax
	flags = my_flag
	offs = 0.

	do i=1,nvar
	  var = vars(i)
	  call nc_get_var_id(ncid,var,var_id)
	  if( var_id == 0 ) then
	    write(6,*) 'no such variable: ',trim(var)
	    stop 'error stop write_variables: no such variable'
	  end if
	  ids(i) = var_id
	  !call nc_var_info(ncid,var_id,.true.)
	  aname = '_FillValue'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    flags(i) = avalue
	  end if
	  aname = 'scale_factor'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    facts(i) = facts(i) * avalue
	  end if
	  aname = 'add_offset'
	  if( nc_has_var_attrib(ncid,var_id,aname) ) then
	    call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	    offs(i) = avalue
	  end if
	  dims(i) = 2
	  call nc_has_vertical_dimension(ncid,var,bvert)
	  if( bvert ) dims(i) = 3
	end do

        call nc_get_time_recs(ncid,nit)
        if( .not. bquiet ) write(6,*) 'time records found: ',nit

	ndd = 0
	if( level > 0 ) ndd = 2
	do i=1,nvar
	  if( ndd == 0 ) ndd = dims(i)
	  if( ndd /= dims(i) .and. level == 0 ) then
	    write(6,*) 'mixing 2d and 3d variables... not possible'
	    write(6,*) dims
	    stop 'error stop write_variables: mixing 2d and 3d'
	  end if
	end do

	lmax = nz1
	if( ndd == 2 ) lmax = 1

	!write(6,*) 'ggu: ',ndd,nz,lmax,level

	if( iformat == 0 ) then
	  iunit = ifileo(30,'out.fem','unform','new')
	else
	  iunit = ifileo(30,'out.fem','form','new')
	end if

	ns = 1
	ns = min(ns,nit)	!loop at least once - for vars without time
	if( bdebug ) nit = min(5,nit)

        do it=ns,nit

	  call create_date_string(ncid,it,datetime)
	  call date2string(datetime,stime)
          if( .not. bquiet ) then
	    write(6,*) 'writing record: ',it,'   ',stime
	  end if

	  call fem_file_write_params(iformat,iunit,dtime
     +                          ,nvers,npnew,lmax
     +                          ,nvar,ntype,datetime)
          call fem_file_write_2header(iformat,iunit,ntype,lmax
     +                  	,hlv,regpar(1:7))

	  do i=1,nvar
	    nzz = nz
	    if( dims(i) == 2 ) nzz = 1
	    call handle_data(ncid,bdebug,vars(i),it,dims(i),flags(i)
     +				,nx,ny,nzz
     +				,x,y
     +				,nxnew,nynew,regpar,ilhkv
     +				,data,femdata,npnew)

	    if( bexpand ) then
	      call reg_expand_3d(nz,nxnew,nynew,lmax,regexpand
     +					,my_flag,femdata)
	      call adjust_reg_vertical(nz,nxnew,nynew,my_flag
     +					,femdata,ilhkv)
	    end if

	    if( facts(i) /= 1. ) then
	      where( femdata /= my_flag ) femdata = femdata * facts(i)
	      where( data /= my_flag ) data = data * facts(i)
	    end if
	    if( offs(i) /= 0. ) then
	      where( femdata /= my_flag ) femdata = femdata + offs(i)
	      where( data /= my_flag ) data = data + offs(i)
	    end if

	    lmax = nzz
	    string = descrps(i)
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,npnew,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,femdata)

	    if( bdebug ) then
	      call write_grd_var(vars(i),it,nx,ny,nzz,nxnew,nynew
     +			,regpar,x,y,data,femdata)
	    end if
	  end do

        end do

	close(iunit)

	nrec = nit - ns + 1
	!write(6,*) 'Total number of records written: ',nrec,nit,ns

	deallocate(femdata)
	deallocate(hd)
	deallocate(ilhkv)

	end

c*****************************************************************

	subroutine handle_data(ncid,bdebug,varname,it,ndims,flag
     +				,nx,ny,nz
     +				,x,y
     +				,nxnew,nynew,regpar,ilhkv
     +				,data,femdata,np)

	implicit none

	integer ncid
	logical bdebug
	character*(*) varname
	integer it
	integer ndims
	real flag
	integer nx,ny,nz
	real x(nx,ny),y(nx,ny)
	integer nxnew,nynew
	real regpar(9)
	integer ilhkv(nxnew*nynew)
	real data(nx,ny,nz)
	real femdata(nz,nxnew*nynew)
	integer np

	logical debug
	integer ndim,nxy,k,iz
	integer nxx,nyy,nzz,nlvddi
	integer dims(10)
	real data2d(nx,ny)
	real femdata2d(nxnew*nynew)
	real cdata(nx*ny,nz)
	!real valnew(nxnew*nynew)
	real, allocatable :: valnew(:,:)
	real, save :: my_flag = -999.
	character*80 file,filename

	logical must_interpol,is_single

	debug = bdebug

	nxy = nx*ny
	ndim = nx*ny*nz
        call nc_get_var_data(ncid,varname,it,ndim,ndims,dims,data)

	if( ndims < 2 .or. ndims > 3 ) then
	  write(6,*) 'error in dimensions: ndims = ',ndims
	  stop 'error stop handle_data: ndims'
	end if

	if( nx /= dims(1) .or. ny /= dims(2) ) then
	  write(6,*) 'error in dimensions: '
	  write(6,*) 'nx,ny given: ',nx,ny
	  write(6,*) 'nx,ny read : ',dims(1),dims(2)
	  stop 'error stop handle_data: dimensions x/y'
	end if

	if( ndims == 3 .and. nz /= dims(3) .or.
     +			ndims == 2 .and. nz /= 1 ) then
	  write(6,*) 'error in dimensions: '
	  write(6,*) 'ndims: ',ndims
	  write(6,*) 'nz given: ',nz
	  write(6,*) 'nz read : ',dims(3)
	  stop 'error stop handle_data: dimensions z'
	end if

	!if( ndims == 3 .and. nz /= 1 ) then
	!  write(6,*) 'ndims,nz: ',ndims,nz
	!  stop 'error stop handle_data: no 3d yet'
	!end if

	where( data == flag ) data = my_flag

	nxx = nx
	nyy = ny
	nzz = nz
	nlvddi = nz

	if( must_interpol() ) then
	  np = nxnew*nynew
	  allocate(valnew(np,nz))
	  if( nz == 1 ) then
	    call do_interpol_2d(nx,ny,data,nxnew,nynew,valnew)
	  else
	    call do_interpol_3d(nx,ny,nz,data,nxnew,nynew,valnew)
	  end if
	  call copy_data_to_fem(nxnew,nynew,nzz,nlvddi,valnew,femdata)
	else
	  call compress_data(nxx,nyy,nzz,data,cdata)	!adjusts nxx,nyy,nzz
	  call copy_data_to_fem(nxx,nyy,nzz,nlvddi,cdata,femdata)
	  np = nxx*nyy
	end if
	call make_ilhkv(np,nlvddi,my_flag,femdata,ilhkv)

	end

c*****************************************************************

	subroutine write_grd_var(varname,it,nx,ny,nz,nxnew,nynew
     +			,regpar,x,y,data,femdata)

	implicit none

	character*(*) varname
	integer it
	integer nx,ny,nz
	integer nxnew,nynew
	real regpar(9)
	real x(nx,ny),y(nx,ny)
	real data(nx,ny,nz)
	real femdata(nz,nxnew*nynew)

	logical bsingle
	integer ns
	real, allocatable :: data2d(:,:)
	real, allocatable :: femdata2d(:)
	real, allocatable :: xs(:),ys(:)
	character*80 filename,file

	!write(6,*) 'nx,ny: ',nx,ny,nxnew,nynew
	allocate(data2d(nx,ny),femdata2d(nxnew*nynew))

	data2d(:,:) = data(:,:,1)
	femdata2d(:) = femdata(1,:)

	bsingle = ( regpar(1) == 0 )

	call make_filename(varname,it,filename)

	file=trim(filename)//'_orig.grd'
	call write_2d_grd(file,nx,ny,x,y,data2d)

	if( bsingle ) then
	  file=trim(filename)//'_single.grd'
	  call get_single_points(0,ns,xs,ys)
	  allocate(xs(ns),ys(ns))
	  call get_single_points(ns,ns,xs,ys)
	  call write_1d_grd(file,ns,xs,ys,femdata2d)
	else
	  file=trim(filename)//'_intp.grd'
	  call write_2d_grd_regular(file,regpar,femdata2d)
	end if

	end

c*****************************************************************

	subroutine setup_facts(n,sfacts,facts)

	implicit none

	integer n
	character*(*) sfacts(n)
	real facts(n)

	integer i,ianz
	real f(1)
	character*80 string

	integer iscanf

	facts = 1.

	do i=1,n
	  string = sfacts(i)
	  if( string == ' ' ) cycle
	  ianz = iscanf(string,f,1)
	  if( ianz /= 1 ) then
	    write(6,*) i,ianz,'  ',string
	    stop 'error stop setup_facts: parse error'
	  end if
	  facts(i) = f(1)
	end do

	end

c*****************************************************************

	subroutine global_information(ncid)

	implicit none

	integer ncid

	integer iftype
	character*80 atext

	call global_info(ncid,'TITLE')
	call global_info(ncid,'source')
	call global_info(ncid,'type')
	!call global_info(ncid,'CDO')
	call global_info(ncid,'institution')

	end

c*****************************************************************

	subroutine global_info(ncid,what)

	implicit none

	integer ncid
	character*(*) what

	character*80 atext

        call nc_get_global_attr(ncid,what,atext)
	if( atext == ' ' ) return

	write(6,*) trim(what),': ',trim(atext)

	end

c*****************************************************************

	subroutine make_filename(varname,it,filename)

	implicit none

	character*(*) varname,filename
	integer it

	integer i
	character*80 sit

	write(sit,'(a1,i5)') '_',it
	do i=1,6
	  if( sit(i:i) == ' ' ) sit(i:i) = '0'
	end do

	filename = 'debug_' // trim(varname) // trim(sit)

	end

c*****************************************************************

	subroutine handle_unusual_coordinates(nx,ny,x,y)

	implicit none

	integer nx,ny
	real x(nx,ny)
	real y(nx,ny)
	
	logical, save :: bneuman = .true.
	real, save :: flag = 1.e+20
	integer ix,iy
	real xx(nx)
	real yy(ny)

	xx = flag
	yy = flag

	if( bneuman ) then
	  do iy=1,ny
	    do ix=1,nx
	      if( x(ix,iy) /= flag ) then
	        if( xx(ix) /= flag .and. xx(ix) /= x(ix,iy) ) goto 99
		xx(ix) = x(ix,iy)
	      end if
	      if( y(ix,iy) /= flag ) then
	        if( yy(iy) /= flag .and. yy(iy) /= y(ix,iy) ) goto 99
		yy(iy) = y(ix,iy)
	      end if
	    end do
	  end do
	  if( any( xx == flag ) ) goto 98
	  if( any( yy == flag ) ) goto 98
	  do iy=1,ny
	    do ix=1,nx
	      x(ix,iy) = xx(ix)
	      y(ix,iy) = yy(iy)
	    end do
	  end do
	end if

	return
   98	continue
	write(6,*) 'some flags in coordinates...'
	write(6,*) xx
	write(6,*) yy
   99	continue
	write(6,*) 'coordinates are not regular...'
	write(6,*) ix,iy,x(ix,iy),y(ix,iy)
	write(6,*) x(ix,iy),y(ix,iy)
	write(6,*) xx(ix),yy(iy)
	stop 'error stop handle_unusual_coordinates: not regular'
	end

c*****************************************************************

	subroutine check_invert(regpar)

	implicit none

	real regpar(9)

	real dx,dy

	dx = regpar(5)
	dy = regpar(6)

	if( dx < 0 .or. dy < 0 ) then
	  write(6,*) 'coordinates are inverted'
	  write(6,*) 'dx,dy: ',dx,dy
	  write(6,*) 'please find the coordinate that is inverted'
	  write(6,*) 'in the nc file and then use NCO routines'
	  write(6,*) 'to invert them in the nc file'
	  write(6,*) 'example: ncpdq -a "-lat" orig.nc invert.nc'
	  write(6,*) 'if lat is the inverted coordinate'
	  write(6,*) '(You might have to install the nco package)'
	  stop 'error stop check_invert: inverted coordinates'
	end if

	end

c*****************************************************************

