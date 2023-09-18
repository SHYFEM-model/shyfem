!
! $Id: subbas.f,v 1.4 2009-04-03 16:38:23 georg Exp $
!
! initialization routines
!
! contents :
!
! subroutine sp13test(nb,nvers)			tests if file is BAS file
! subroutine sp13rr(nb,nknddi,nelddi)		unformatted read from lagoon
! subroutine sp13uw(nb)				unformatted write to lagoon
! subroutine sp13ts(nvers,nb,n)			test write to unit nb
! subroutine dist_3d(nlvddi,r3v,kn,nbdim,values)
!
! revision log :
!
! 31.05.1997	ggu	unnecessary routines deleted
! 27.06.1997	ggu	bas routines into own file
! 18.09.2007    ggu     bug fix in dist_3d
! 02.04.2009	ggu	error messages changed
! 12.01.2011	ggu	debug routine introduced (sp13ts)
! 23.10.2014	ggu	introduced ftype and nvers = 4
! 04.01.2015	ggu	new routine sp13_get_par()
! 31.03.2015	ggu	set iarnv on read
! 25.05.2015	ggu	module introduced
! 02.10.2015	ggu	in basin_open_file eliminated double read (bug)
! 02.10.2015	ggu	new routines is_depth_unique(), estimate_ngr()
! 01.05.2016	ggu	new routines basin_has_basin()
! 20.05.2016	ggu	estimate_ngr() returns exact ngr
! 10.06.2016	ggu	new routine for inserting regular grid
!
!***********************************************************
!***********************************************************
!***********************************************************

!==================================================================
        module basin
!==================================================================

        use fil

        implicit none

        integer, private, save :: nkn_basin = 0
        integer, private, save :: nel_basin = 0

        integer, save :: nkn = 0
        integer, save :: nel = 0
        integer, save :: ngr = 0
        integer, save :: mbw = 0

        integer, save :: nkndi = 0	!these are needed when nkn changes
        integer, save :: neldi = 0
        !integer, save :: ngrdi = 0
        !integer, save :: mbwdi = 0

        double precision, save :: dcorbas = 0.
        double precision, save :: dirnbas = 0.

        character*80, save :: descrr = ' '

        integer, save, allocatable :: nen3v(:,:)
        integer, save, allocatable :: ipev(:)
        integer, save, allocatable :: ipv(:)
        integer, save, allocatable :: iarv(:)
        integer, save, allocatable :: iarnv(:)

        double precision, save, allocatable :: xgv(:)
        double precision, save, allocatable :: ygv(:)
        double precision, save, allocatable :: hm3v(:,:)

        INTERFACE basin_read
        MODULE PROCEDURE basin_read_by_file, basin_read_by_unit
        END INTERFACE

        INTERFACE basin_is_basin
        MODULE PROCEDURE basin_is_basin_by_file,basin_is_basin_by_unit
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine basin_init(nk,ne,nkex,neex)

	integer nk,ne
	integer, optional :: nkex	!number nodes extended
	integer, optional :: neex	!number elems extended
	if( nk == nkn_basin .and. ne == nel_basin ) return

        if( ne > 0 .or. nk > 0 ) then
          if( ne == 0 .or. nk == 0 ) then
            write(6,*) 'nel,nkn: ',ne,nk
            stop 'error stop basin_init: incompatible parameters'
          end if
        end if

	if( nkn_basin > 0 ) then
	  deallocate(nen3v)
	  deallocate(ipev)
	  deallocate(ipv)
	  deallocate(iarv)
	  deallocate(iarnv)
	  deallocate(xgv)
	  deallocate(ygv)
	  deallocate(hm3v)
	end if

	nkn = nk
	nel = ne
	nkndi = nk
	neldi = ne
	nkn_basin = nk
	nel_basin = ne

	if( nk == 0 ) return

        if(present(neex))then
          allocate(nen3v(3,neex))
        else
	  allocate(nen3v(3,nel))
        end if
	allocate(ipev(nel))
	if(present(nkex)) then
	  allocate(ipv(nkex))
	  allocate(xgv(nkex))
	  allocate(ygv(nkex))
	else
	  allocate(ipv(nkn))
	  allocate(xgv(nkn))
	  allocate(ygv(nkn))
	end if
	allocate(iarv(nel))
	allocate(iarnv(nkn))
        if(present(neex))then
	  allocate(hm3v(3,neex))
        else
	  allocate(hm3v(3,nel))
        end if

	end subroutine basin_init

!***********************************************************

	function basin_open_file(file)

! opens file or returns 0
	
	integer basin_open_file
	character*(*) file

	integer iunit,ios

	basin_open_file = 0

	iunit = ifileo(0,file,'unform','old')
	if( iunit .le. 0 ) return
	!open(iunit,file=file,status='old',form='unformatted',iostat=ios)
	!if( ios /= 0 ) return

	basin_open_file = iunit

	end function basin_open_file

!***********************************************************

	subroutine basin_read_by_file(file)

! reads basin or fails if error

	character*(*) file

	integer iunit

	iunit = basin_open_file(file)

	if( iunit .le. 0 ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop basin_read_by_file: cannot open file'
	end if

	call basin_read_by_unit(iunit)

	close(iunit)

	write(6,*) 'finished reading basin: ',trim(file)

	end subroutine basin_read_by_file

!***********************************************************

	subroutine basin_read_by_unit(iunit)

! reads basin or fails if error

	integer iunit

	integer nk,ne

	call sp13_get_par(iunit,nk,ne,ngr,mbw)
	call basin_init(nk,ne)			!here we set nkn, nel
	rewind(iunit)
	call sp13rr(iunit,nkn,nel)
	!write(6,*) 'finished basin_read (module)'

	end subroutine basin_read_by_unit

!***********************************************************

	function basin_has_basin()

! checks if basin has been read

	logical basin_has_basin

	basin_has_basin = nkn_basin > 0 .and. nel_basin > 0

	end function basin_has_basin

!***********************************************************

	subroutine basin_get_dimension(nk,ne)

! returns dimension of arrays (already allocated)

	integer nk,ne

	nk = nkndi
	ne = neldi

	end subroutine basin_get_dimension

!***********************************************************

	function basin_is_basin_by_unit(iunit)

	logical basin_is_basin_by_unit
	integer iunit

	integer nvers

	call sp13test(iunit,nvers)

	basin_is_basin_by_unit = nvers > 0

	end function basin_is_basin_by_unit

!***********************************************************

	function basin_is_basin_by_file(file)

	logical basin_is_basin_by_file
	character*(*) file

	integer iunit

	basin_is_basin_by_file = .false.

	iunit = basin_open_file(file)
	if( iunit .le. 0 ) return

	basin_is_basin_by_file = basin_is_basin_by_unit(iunit)

	close(iunit)

	end function basin_is_basin_by_file

!***********************************************************

        subroutine basin_check(text)

        implicit none

        character*(*), optional :: text

        integer ie,k,ii
        character*80 string

        string = ' '
        if( present(text) ) string = text

        write(6,*) 'checking basin data: ',trim(string)
        write(6,*) nkn,nel,ngr,mbw

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .le. 0 .or. k .gt. nkn ) then
              write(6,*) ii,ie,k
              write(6,*) 'error checking basin: ',trim(string)
              stop 'error stop basin_check: nen3v'
            end if
          end do
        end do

        write(6,*) 'basin data is ok...'

        end subroutine basin_check

!***********************************************************
!***********************************************************
!***********************************************************

	subroutine sp13test(nb,nvers)

! tests if file is BAS file
!
! nvers > 0 if file is BAS file

	implicit none

	integer nb	!unit number
	integer nvers	!version found (return) (<=0 if error or no BAS file)

	integer ftype,nversm
	parameter (ftype=789233567,nversm=4)

	integer ntype,nversa

	nvers = 0

	if(nb.le.0) return

!-----------------------------------------------------------
! try new format with ftype information
!-----------------------------------------------------------

	rewind(nb)
	read(nb,err=1,end=1) ntype,nversa
	if( ntype .ne. ftype ) return
	if( nversa .le. 3 .or. nversa .gt. nversm ) nversa = -abs(nversa)

	nvers = nversa
	return

!-----------------------------------------------------------
! try old format without ftype information - nvers must be 3
!-----------------------------------------------------------

    1	continue
	rewind(nb)
	read(nb,err=2,end=2) nversa
	if( nversa .ne. 3 ) nversa = -abs(nversa)

	nvers = nversa
	return

!-----------------------------------------------------------
! definitely no BAS file
!-----------------------------------------------------------

    2	continue
	return
	end

!***********************************************************

	subroutine sp13_get_par(nb,nkn,nel,ngr,mbw)

! unformatted read from lagoon file
!
! iunit		unit number of file to be read

	implicit none

	integer nb
	integer nkn,nel,ngr,mbw

	integer nvers
	character*80 file

	file = ' '
	if( nb > 0 ) inquire(nb,name=file)

	call sp13test(nb,nvers)

	if(nvers.eq.0) goto 99
	if(nvers.lt.0) goto 98

	read(nb) nkn,nel,ngr,mbw

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	if( nb > 0 ) write(6,*) 'file name = ',trim(file)
	stop 'error stop : sp13_get_par'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	if( nb > 0 ) write(6,*) 'file name = ',trim(file)
	stop 'error stop : sp13_get_par'
   97	continue

	end

!***********************************************************

	subroutine sp13rr(nb,nknddi,nelddi)

! unformatted read from lagoon file
!
! iunit		unit number of file to be read

	implicit none

	integer nb,nknddi,nelddi

	include 'param.h'

	integer i,ii,nvers

	call sp13test(nb,nvers)

	if(nvers.eq.0) goto 99
	if(nvers.lt.0) goto 98

	read(nb) nkn,nel,ngr,mbw
	read(nb) dcorbas,dirnbas
	read(nb) descrr

	if(nkn.gt.nknddi.or.nel.gt.nelddi) goto 97

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

!	call sp13ts(nvers,79,0)

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	stop 'error stop sp13rr: error reading file'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	write(6,*) 'nvers = ',-nvers
	stop 'error stop sp13rr: error in version'
   97	continue
	write(6,*) 'nknddi,nelddi :',nknddi,nelddi
	write(6,*) 'nkn,nel       :',nkn,nel
	write(6,*) 'ngr,mbw       :',ngr,mbw
	stop 'error stop sp13rr: dimension error'
	end

!***********************************************************

	subroutine sp13uw(nb)

! unformatted write to lagoon file
!
! nb		unit number for write

	implicit none

	integer nb

	include 'param.h'

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

!	call sp13ts(nvers,78,0)

	return
   99	continue
	write(6,*) 'Writing basin...'
	write(6,*) 'Cannot write bas file on unit :',nb
	stop 'error stop : sp13uw'
	end

!*************************************************

	subroutine sp13ts(nvers,nb,n)

! test write to unit nb

! writes first n values, if n=0 -> all values

	implicit none

	integer nvers,nb,n

	include 'param.h'

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

!*************************************************

	subroutine bas_info

	implicit none

	include 'param.h'

        write(6,*)
        write(6,*) trim(descrr)
        write(6,*)
        write(6,*) ' nkn = ',nkndi,'  nel = ',neldi
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcorbas,'  dirn = ',dirnbas
        write(6,*)

	end

!*************************************************

	subroutine bas_get_geom(dcor,dirn)

	implicit none

	include 'param.h'

	double precision dcor,dirn

	dcor = dcorbas
	dirn = dirnbas

	end

!*************************************************

	subroutine bas_get_para(nkna,nela,ngra,mbwa)

	implicit none

	include 'param.h'

	integer nkna,nela,ngra,mbwa

	nkna = nkn
	nela = nel
	ngra = ngr
	mbwa = mbw

	end

!*************************************************

	subroutine bas_get_minmax(xmin,ymin,xmax,ymax)

	implicit none

	include 'param.h'

	double precision xmin,ymin,xmax,ymax

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

!*************************************************

        function is_depth_unique()

	implicit none

        logical is_depth_unique

	integer ie,ii,k
	double precision h
        double precision :: flag = -999.
        double precision haux(nkn)

        haux = flag
	is_depth_unique = .false.

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            h = hm3v(ii,ie)
            if( haux(k) == flag ) haux(k) = h
	    if( h /= haux(k) ) return
          end do
        end do

	is_depth_unique = .true.

        end function is_depth_unique

!*************************************************
!*************************************************
!*************************************************

	subroutine bas_insert_regular(regpar)

! inserts regular basin (boxes) into basin structure

	implicit none

	double precision regpar(7)

	integer nx,ny,ix,iy
	integer k,ie,nk,ne
	integer k1,k2,k3,k4,km
	double precision x0,y0,dx,dy
	integer, allocatable :: indexv(:,:)
	integer, allocatable :: indexm(:,:)

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        x0 = regpar(3)
        y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)

	nk = nx*ny + (nx-1)*(ny-1)
	ne = 4*(nx-1)*(ny-1)

	allocate(indexv(nx,ny))
	allocate(indexm(nx,ny))

	call basin_init(nk,ne)

	mbw = 0
	descrr = 'regular generated grid'
	iarv = 0
	iarnv = 0
	hm3v = 0.

	k = 0
        do iy=1,ny
          do ix=1,nx
	    k = k + 1
	    ipv(k) = k
	    indexv(ix,iy) = k
	    xgv(k) = x0 + (ix-1)*dx
	    ygv(k) = y0 + (iy-1)*dy
          end do
        end do

        do iy=2,ny
          do ix=2,nx
	    k = k + 1
	    ipv(k) = k
	    indexm(ix,iy) = k
	    xgv(k) = x0 + (ix-1.5)*dx
	    ygv(k) = y0 + (iy-1.5)*dy
          end do
        end do

	ie = 0
        do iy=2,ny
          do ix=2,nx
	    k1 = indexv(ix-1,iy-1)
	    k2 = indexv(ix,iy-1)
	    k3 = indexv(ix,iy)
	    k4 = indexv(ix-1,iy)
	    km = indexm(ix,iy)
	    ie = ie + 1
	    ipev(ie) = ie
	    nen3v(1,ie) = km
	    nen3v(2,ie) = k1
	    nen3v(3,ie) = k2
	    ie = ie + 1
	    ipev(ie) = ie
	    nen3v(1,ie) = km
	    nen3v(2,ie) = k2
	    nen3v(3,ie) = k3
	    ie = ie + 1
	    ipev(ie) = ie
	    nen3v(1,ie) = km
	    nen3v(2,ie) = k3
	    nen3v(3,ie) = k4
	    ie = ie + 1
	    ipev(ie) = ie
	    nen3v(1,ie) = km
	    nen3v(2,ie) = k4
	    nen3v(3,ie) = k1
          end do
        end do

        call estimate_ngr(ngr)

	deallocate(indexv)
	deallocate(indexm)

	end

!*************************************************
!*************************************************
!*************************************************

        subroutine estimate_ngr(ngrade)

! estimates grade of basin - estimate is exact

        implicit none

	integer ngrade

        integer ng(nkn)

	call compute_ng(ngrade,ng)

	end

!*************************************************

        subroutine compute_ng(ngrade,ng)

! computes grade of basin and nuber of grades per node

        implicit none

	integer ngrade
        integer ng(nkn)

        integer ii,ie,k
	integer k1,k2,n
	integer, allocatable :: ngv(:,:)

	ng = 0

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            ng(k) = ng(k) + 1
          end do
        end do

	ngrade = maxval(ng) + 1			!first guess

	allocate(ngv(0:2*ngrade,nkn))
	ngv = 0

        do ie=1,nel
          do ii=1,3
            k1 = nen3v(ii,ie)
            k2 = nen3v(mod(ii,3)+1,ie)
	    call ng_insert(k1,k2,ngrade,nkn,ngv)
	    call ng_insert(k2,k1,ngrade,nkn,ngv)
	  end do
	end do

	do k=1,nkn
	  n = ngv(0,k)
	  if( n == 0 ) then			!inner node
	    !nothing to do
	  else if( n == 2 ) then		!boundary node
	    ng(k) = ng(k) + 1
	  else
	    write(6,*) 'wrong connectivity: ',k,n
	    stop 'error stop estimate_ngr: internal error'
	  end if
	end do

	deallocate(ngv)

	ngrade = maxval(ng)

        end

!*************************************************

	subroutine ng_insert(k1,k2,ng,nkn,ngv)

	implicit none

	integer k1,k2
	integer ng,nkn
	integer ngv(0:2*ng,nkn)

	integer i,n

	n = ngv(0,k1)

	do i=1,n
	  if( ngv(i,k1) == k2 ) then
	    ngv(i,k1) = ngv(n,k1)
	    ngv(0,k1) = ngv(0,k1) - 1
	    return
	  end if
	end do

	n = n + 1
	ngv(0,k1) = n
	ngv(n,k1) = k2

        end

!*************************************************

	subroutine dist_3d(nlvddi,r3v,kn,nbdim,values)

	!use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi
	double precision r3v(nlvddi,nkn)
	integer kn
	integer nbdim
	double precision values(nkn)

	integer l,lmax

	if( nbdim .eq. 0 ) then
	  lmax = 1
	else
	  lmax = min(nbdim,nlvddi)
	end if

	do l=1,lmax
	  r3v(l,kn) = values(l)
	end do
	  
	do l=lmax+1,nlvddi
	  !r3v(l,kn) = values(nbdim)	!BUGFIX
	  r3v(l,kn) = r3v(lmax,kn)
	end do
	  
	end

!*******************************************************************

!==================================================================
        end module basin
!==================================================================
