c
c revision log :
c
C 16.12.2017    ggu     started from scratch
c
c****************************************************************

        program shypart

c partitions grd file

	use mod_geom
	use evgeom
	use basin
	use clo
	use basutil
	use shympi

	implicit none

	integer nc,ncol,itercol
	character*80 file

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	call shympi_init(.false.)
	!call shyfem_copyright('shybas - elaborating a BAS grid')

	call basutil_init('BAS')

	call clo_get_file(1,file)
        if( file == ' ' ) call clo_usage
	call read_command_line_file(file)

c-----------------------------------------------------------------
c initialize modules
c-----------------------------------------------------------------

	call ev_init(nel)
	call set_ev

	call mod_geom_init(nkn,nel,ngr)
	call set_geom

c-----------------------------------------------------------------
c partition
c-----------------------------------------------------------------

	itercol = 0
	iarv = 0
	iarnv = 0
	call seed_partition(ncol,iarnv)

	do
	  itercol = itercol + 1
	  call flood_fill(iarnv)
	  call color_stats(iarnv)
	  if( itercol > 5 ) exit
	  call find_new_center(iarnv,ncol)
	end do

	call fill_elements(iarnv,iarv)
	!call equilibrate_partition(iarnv)

c-----------------------------------------------------------------
c write grd file
c-----------------------------------------------------------------

	call write_grd

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine equilibrate_partition(ia)

	use basin

	implicit none

	integer ia(nkn)

	integer ncol,ncdiff
	integer k,i,iexch
	integer iminv(1),imaxv(1)
	integer imin,imax
	integer ncmin,ncmax,ncmed
	integer, allocatable :: nccol(:)
	integer, allocatable :: index(:)

	ncol = maxval(ia)
	allocate(nccol(ncol))
	allocate(index(ncol))

	ncmed = nkn/ncol

	do 
	call get_col_stats(ncol,nccol,nkn,ia)
	call isort(ncol,nccol,index)
	call get_minmax(ncol,nccol,index,imin,imax,ncmin,ncmax)
	ncdiff = ncmax-ncmin
	if( ncdiff < 10 ) exit
	print*,'ncdiff = ',ncdiff
	print*,nccol

	do i=1,ncol
	  imin = index(i)
	  ncmin = nccol(imin)
	  !print *,'before: ',imin,imax,ncmin,ncmax,ncmed
	  call exchange(ia,imin,imax,ncmin,ncmax,ncmed,iexch)
	  if( iexch > 0 ) exit
	end do

	print *,'exchanged: ',ncmin,ncmax,ncmed,iexch

	end do

	end

c*******************************************************************

	subroutine get_minmax(ncol,nccol,index,imin,imax,ncmin,ncmax)

	implicit none

	integer ncol
	integer nccol(ncol)
	integer index(ncol)
	integer imin,imax,ncmin,ncmax

	  imin = index(1)
	  imax = index(ncol)
	  ncmin = nccol(imin)
	  ncmax = nccol(imax)

	end

c*******************************************************************

	subroutine exchange(ia,imin,imax,ncmin,ncmax,ncmed,iexch)

	use basin

	implicit none

	integer ia(nkn)
	integer imin,imax
	integer ncmin,ncmax,ncmed
	integer iexch

	integer ie,k,ii,ix,iter
	integer icmin,icmax

	iexch = 0
	iter = 0

	do while( ncmin < ncmax )
	 ix = 0
	 do ie=1,nel
	  icmin = 0
	  icmax = 0
	  do ii= 1,3
	    k = nen3v(ii,ie)
	    if( ia(k) == imin ) icmin = icmin + 1
	    if( ia(k) == imax ) icmax = icmax + 1
	  end do
	  if( icmin == 2 .and. icmax == 1 ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( ia(k) == imax ) ia(k) = -imin
	    end do
	    ncmin = ncmin + 1
	    ncmax = ncmax - 1
	    ix = ix + 1
	  end if
	  if( ncmin > ncmax ) exit
	 end do
	 where( ia < 0 ) ia = -ia
	 iter = iter + 1
	 if( ix == 0 ) exit
	 print*,'iter: ',iter,ix
	 if( ncmin > ncmax ) exit
	 iexch = iexch + ix
	end do

	end

c*******************************************************************

	subroutine color_stats(ia)

	use basin

	integer ia(nkn)

	integer icmin,icmax
	integer k
	integer, allocatable :: nccol(:)

	icmin = minval(ia)
	icmax = maxval(ia)

	allocate(nccol(icmax))
	nccol = 0

	if( icmin /= 1 ) stop 'error stop color_stats: icmin /= 1'

	ncol = icmax
	write(6,*) 'colors min/max: ',icmin,icmax,ncol

	call get_col_stats(ncol,nccol,nkn,ia)

	do ic=icmin,icmax
	  print *, ic,nccol(ic)
	end do

	end

c*******************************************************************

	subroutine get_col_stats(ncol,nccol,nkn,ia)

	implicit none

	integer ncol
	integer nccol(ncol)
	integer nkn
	integer ia(nkn)

	integer k,ic

	nccol = 0

	do k=1,nkn
	  ic = ia(k)
	  nccol(ic) = nccol(ic) + 1
	end do
	
	end

c*******************************************************************

	subroutine seed_partition(ncol,ian)

	use basin

	implicit none

	integer ncol
	integer ian(nkn)

	integer i,icol,ke,ki
	integer, save :: ninit(6) = (/582,1107,1879,2117,3165,4060/)

	integer ipint

	icol = 0
	do i=1,size(ninit)
	  icol = icol + 1
	  ke = ninit(i)
	  ki = ipint(ke)
	  if( ki <= 0 ) stop 'error stop seed_partition: ki==0'
	  ian(ki) = icol
	end do

	ncol = icol

	end

c*******************************************************************

	subroutine find_new_center(color,ncol)

	use basin

	implicit none

	integer ncol
	integer color(nkn)

	integer ic,icmax,k,iter
	integer itmin,ikmin
	integer color_aux(nkn)
	integer color_grade(nkn)
	integer color_out(nkn)
	integer color_iter(nkn)

	color_out = 0
	color_grade = color
	call fill_areas(color_grade,ncol)

	do ic=1,ncol
	  color_aux = -1
	  where( color == ic ) color_aux = color_grade
	  icmax = maxval(color_aux)
	  !write(6,*) ic,icmax
	  itmin = nkn
	  ikmin = 0
	  do k=1,nkn
	    if( color_aux(k) == icmax ) then
	      color_iter = -1
	      where( color == ic ) color_iter = 0
	      color_iter(k) = 1
	      call flood_fill_iter(color_iter,iter)
	      !write(6,*) ic,icmax,k,iter
	      if( iter < itmin ) then
		itmin = iter
		ikmin = k
	      end if
	    end if
	  end do
	  write(6,*) 'final: ',ic,icmax,ikmin,itmin
	  color_out(ikmin) = ic
	end do

	color = color_out

	end

c*******************************************************************

	subroutine fill_elements(ian,iae)

	use basin

	implicit none

	integer ian(nkn)
	integer iae(nel)

	integer ie,ii,k,col
	logical bequal

	iae = 0

	do ie=1,nel
	  col = 0
	  bequal = .true.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( col == 0 ) col = ian(k)
	    if( col /= ian(k) ) bequal = .false.
	  end do
	  if( bequal ) iae(ie) = col
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine read_command_line_file(file)

	use basin
	use basutil

	implicit none

	character*(*) file
	logical is_grd_file

	if( basin_is_basin(file) ) then
	  write(6,*) 'reading BAS file ',trim(file)
	  call basin_read(file)
	  breadbas = .true.
	else if( is_grd_file(file) ) then
	  write(6,*) 'reading GRD file ',trim(file)
	  call grd_read(file)
	  call grd_to_basin
	  call estimate_ngr(ngr)
	  breadbas = .false.
	else
	  write(6,*) 'Cannot read this file: ',trim(file)
	  stop 'error stop read_given_file: format not recognized'
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine write_grd

c writes grd file extracting info from bas file

	implicit none

        call basin_to_grd
        call grd_write('basin.grd')

        write(6,*) 'The basin has been written to basin.grd'

	end

c*******************************************************************

	subroutine node_test
	end

c*******************************************************************

