c
c $Id: subbas.f,v 1.4 2009-04-03 16:38:23 georg Exp $
c
c initialization routines
c
c contents :
c
c subroutine sp13test(nb,nvers)			tests if file is BAS file
c subroutine sp13rr(nb,nkndim,neldim)		unformatted read from lagoon
c subroutine sp13uw(nb)				unformatted write to lagoon
c subroutine sp13ts(nvers,nb,n)			test write to unit nb
c
c revision log :
c
c 31.05.1997	ggu	unnecessary routines deleted
c 27.06.1997	ggu	bas routines into own file
c 02.04.2009	ggu	error messages changed
c 12.01.2011	ggu	debug routine introduced (sp13ts)
c 23.10.2014	ggu	introduced ftype and nvers = 4
c 04.01.2015	ggu	new routine sp13_get_par()
c 31.03.2015	ggu	set iarnv on read
c 25.05.2015	ggu	module introduced
c
c***********************************************************
c***********************************************************
c***********************************************************

!==================================================================
        module basin
!==================================================================

        implicit none

        integer, save :: nkn = 0
        integer, save :: nel = 0
        integer, save :: ngr = 0
        integer, save :: mbw = 0

        real, save :: dcorbas = 0.
        real, save :: dirnbas = 0.

        character*80, save :: descrr = ' '

        integer, save, allocatable :: nen3v(:,:)
        integer, save, allocatable :: ipev(:)
        integer, save, allocatable :: ipv(:)
        integer, save, allocatable :: iarv(:)
        integer, save, allocatable :: iarnv(:)

        real, save, allocatable :: xgv(:)
        real, save, allocatable :: ygv(:)
        real, save, allocatable :: hm3v(:,:)

!==================================================================
        contains
!==================================================================

	subroutine basin_init(nk,ne)

	integer nk,ne

	if( nkn == nk .and. nel == ne ) return

	if( nkn > 0 .or. nel > 0 ) then
	end if

	nkn = nk
	nel = ne

	if( nkn > 0 .or. nel > 0 ) then
	  allocate(nen3v(3,nel))
	  allocate(ipev(nel))
	  allocate(ipv(nkn))
	  allocate(iarv(nel))
	  allocate(iarnv(nkn))
	  allocate(xgv(nkn))
	  allocate(ygv(nkn))
	  allocate(hm3v(3,nel))
	end if

	end subroutine basin_init

c***********************************************************

	subroutine basin_read(iunit)

! reads basin or fails if error

	integer iunit

	integer nk,ne

	call sp13_get_par(iunit,nk,ne,ngr,mbw)
	call basin_init(nk,ne)
	rewind(iunit)
	call sp13rr(iunit,nkn,nel)

	end subroutine basin_read

c***********************************************************

	function is_basin(iunit)

	logical is_basin
	integer iunit

	integer nvers

	call sp13test(iunit,nvers)

	is_basin = nvers > 0

	end function is_basin

c***********************************************************

	function is_basin_file(file)

	logical is_basin_file
	character*(*) file

	integer ios,iunit

	is_basin_file = .false.

	open(iunit,file=file,status='old',form='unformatted',iostat=ios)

	if( ios /= 0 ) return
	is_basin_file = is_basin(iunit)

	close(iunit)

	end function is_basin_file

!==================================================================
        end module basin
!==================================================================

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine sp13test(nb,nvers)

c tests if file is BAS file
c
c nvers > 0 if file is BAS file

	implicit none

	integer nb	!unit number
	integer nvers	!version found (return) (<=0 if error or no BAS file)

	integer ftype,nversm
	parameter (ftype=789233567,nversm=4)

	integer ntype,nversa

	nvers = 0

	if(nb.le.0) return

c-----------------------------------------------------------
c try new format with ftype information
c-----------------------------------------------------------

	rewind(nb)
	read(nb,err=1,end=1) ntype,nversa
	if( ntype .ne. ftype ) return
	if( nversa .le. 3 .or. nversa .gt. nversm ) nversa = -abs(nversa)

	nvers = nversa
	return

c-----------------------------------------------------------
c try old format without ftype information - nvers must be 3
c-----------------------------------------------------------

    1	continue
	rewind(nb)
	read(nb,err=2,end=2) nversa
	if( nversa .ne. 3 ) nversa = -abs(nversa)

	nvers = nversa
	return

c-----------------------------------------------------------
c definitely no BAS file
c-----------------------------------------------------------

    2	continue
	return
	end

c***********************************************************

	subroutine sp13_get_par(nb,nkn,nel,ngr,mbw)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	implicit none

	integer nb
	integer nkn,nel,ngr,mbw

	integer nvers

	call sp13test(nb,nvers)

	if(nvers.eq.0) goto 99
	if(nvers.lt.0) goto 98

	read(nb) nkn,nel,ngr,mbw

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	stop 'error stop : sp13_get_par'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	stop 'error stop : sp13_get_par'
   97	continue

	end

c***********************************************************

	subroutine sp13rr(nb,nkndi,neldi)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	implicit none

	integer nb,nkndi,neldi

	include 'param.h'
	include 'basin.h'

	integer i,ii,nvers

	call sp13test(nb,nvers)

	if(nvers.eq.0) goto 99
	if(nvers.lt.0) goto 98

	read(nb) nkn,nel,ngr,mbw
	read(nb) dcorbas,dirnbas
	read(nb) descrr

	if(nkn.gt.nkndi.or.nel.gt.neldi) goto 97

	read(nb)((nen3v(ii,i),ii=1,3),i=1,nel)
	read(nb)(ipv(i),i=1,nkn)
	read(nb)(ipev(i),i=1,nel)
	read(nb)(iarv(i),i=1,nel)

	read(nb)(xgv(i),i=1,nkn)
	read(nb)(ygv(i),i=1,nkn)
	read(nb)((hm3v(ii,i),ii=1,3),i=1,nel)

	do i=1,nkn
	  iarnv(i) = 0
	end do

c	call sp13ts(nvers,79,0)

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	stop 'error stop sp13rr: error reading file'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	write(6,*) 'nvers = ',-nvers
	stop 'error stop sp13rr: error in version'
   97	continue
	write(6,*) 'nkndim,neldim :',nkndi,neldi
	write(6,*) 'nkn,nel       :',nkn,nel
	write(6,*) 'ngr,mbw       :',ngr,mbw
	stop 'error stop sp13rr: dimension error'
	end

c***********************************************************

	subroutine sp13uw(nb)

c unformatted write to lagoon file
c
c nb		unit number for write

	implicit none

	integer nb

	include 'param.h'
	include 'basin.h'

	integer i,ii

	integer ftype,nversm
	parameter (ftype=789233567,nversm=4)

	if(nb.le.0) goto 99

	rewind(nb)

	write(nb) ftype,nversm
	write(nb) nkn,nel,ngr,mbw
	write(nb) dcorbas,dirnbas
	write(nb) descrr

	write(nb)((nen3v(ii,i),ii=1,3),i=1,nel)
	write(nb)(ipv(i),i=1,nkn)
	write(nb)(ipev(i),i=1,nel)
	write(nb)(iarv(i),i=1,nel)

	write(nb)(xgv(i),i=1,nkn)
	write(nb)(ygv(i),i=1,nkn)
	write(nb)((hm3v(ii,i),ii=1,3),i=1,nel)

c	call sp13ts(nvers,78,0)

	return
   99	continue
	write(6,*) 'Writing basin...'
	write(6,*) 'Cannot write bas file on unit :',nb
	stop 'error stop : sp13uw'
	end

c*************************************************

	subroutine sp13ts(nvers,nb,n)

c test write to unit nb

c writes first n values, if n=0 -> all values

	implicit none

	integer nvers,nb,n

	include 'param.h'
	include 'basin.h'

	integer i,ii
	integer nkn1,nel1

	nkn1 = min(nkn,n)
	if( nkn1 .le. 0 ) nkn1 = nkn
	nel1 = min(nel,n)
	if( nel1 .le. 0 ) nel1 = nel

	rewind(nb)

	write(nb,*) 'sp13ts:'
	write(nb,*) nvers
	write(nb,*) nkn,nel,ngr,mbw
	write(nb,*) dcorbas,dirnbas
	write(nb,*) descrr

	write(nb,*)((nen3v(ii,i),ii=1,3),i=1,nel1)
	write(nb,*)(ipv(i),i=1,nkn1)
	write(nb,*)(ipev(i),i=1,nel1)
	write(nb,*)(iarv(i),i=1,nel1)

	write(nb,*)(xgv(i),i=1,nkn1)
	write(nb,*)(ygv(i),i=1,nkn1)
	write(nb,*)((hm3v(ii,i),ii=1,3),i=1,nel1)

	return
	end

c*************************************************

	subroutine bas_info

	implicit none

	include 'param.h'
	include 'basin.h'

        write(6,*)
        write(6,*) trim(descrr)
        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcorbas,'  dirn = ',dirnbas
        write(6,*)

	end

c*************************************************

	subroutine bas_get_geom(dcor,dirn)

	implicit none

	include 'param.h'
	include 'basin.h'

	real dcor,dirn

	dcor = dcorbas
	dirn = dirnbas

	end

c*************************************************

	subroutine bas_get_para(nkna,nela,ngra,mbwa)

	implicit none

	include 'param.h'
	include 'basin.h'

	integer nkna,nela,ngra,mbwa

	nkna = nkn
	nela = nel
	ngra = ngr
	mbwa = mbw

	end

c*************************************************

	subroutine bas_get_minmax(xmin,ymin,xmax,ymax)

	implicit none

	include 'param.h'
	include 'basin.h'

	real xmin,ymin,xmax,ymax

	integer k

	xmin = xgv(1)
	xmax = xgv(1)
	ymin = ygv(1)
	ymax = ygv(1)

	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  xmax = max(xmax,xgv(k))
	  ymin = min(ymin,ygv(k))
	  ymax = max(ymax,ygv(k))
	end do

	end

c*************************************************

