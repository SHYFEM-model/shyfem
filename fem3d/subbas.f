c
c $Id: subbas.f,v 1.4 2009-04-03 16:38:23 georg Exp $
c
c initialization routines
c
c contents :
c
c subroutine sp13test(nb,nvers)			tests if file is BAS file
c subroutine sp13rr(nb,nknddi,nelddi)		unformatted read from lagoon
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
c 02.10.2015	ggu	in basin_open_file eliminated double read (bug)
c 02.10.2015	ggu	new routines is_depth_unique(), estimate_ngr()
c 01.05.2016	ggu	new routines basin_has_basin()
c 20.05.2016	ggu	estimate_ngr() returns exact ngr
c 10.06.2016	ggu	new routine for inserting regular grid
c 23.09.2016	ggu	new routines to check if basin has been read
c
c***********************************************************
c***********************************************************
c***********************************************************

!==================================================================
        module basin
!==================================================================

        implicit none

        integer, private, save :: nkn_basin = 0
        integer, private, save :: nel_basin = 0

        logical, private, save :: bbasinread = .false.	! basin has been read

        integer, save :: nkn = 0
        integer, save :: nel = 0
        integer, save :: ngr = 0
        integer, save :: mbw = 0

        integer, save :: nkndi = 0	!these are needed when nkn changes
        integer, save :: neldi = 0
        !integer, save :: ngrdi = 0
        !integer, save :: mbwdi = 0

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

        INTERFACE basin_read
        MODULE PROCEDURE basin_read_by_file, basin_read_by_unit
        END INTERFACE

        INTERFACE basin_is_basin
        MODULE PROCEDURE basin_is_basin_by_file,basin_is_basin_by_unit
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine basin_init(nk,ne)

	integer nk,ne

	if( nk == nkn_basin .and. ne == nel_basin ) return

        if( ne > 0 .or. nk > 0 ) then
          !if( ne == 0 .or. nk == 0 ) then
          !  write(6,*) 'nel,nkn: ',ne,nk
          !  stop 'error stop basin_init: incompatible parameters'
          !end if
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

	allocate(nen3v(3,nel))
	allocate(ipev(nel))
	allocate(ipv(nkn))
	allocate(iarv(nel))
	allocate(iarnv(nkn))
	allocate(xgv(nkn))
	allocate(ygv(nkn))
	allocate(hm3v(3,nel))

	end subroutine basin_init

c***********************************************************

	function basin_open_file(file)

! opens file or returns 0
	
	integer basin_open_file
	character*(*) file

	integer iunit,ios
	integer ifileo

	basin_open_file = 0

	iunit = ifileo(0,file,'unform','old')
	if( iunit .le. 0 ) return
	!open(iunit,file=file,status='old',form='unformatted',iostat=ios)
	!if( ios /= 0 ) return

	basin_open_file = iunit

	end function basin_open_file

c***********************************************************

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

c***********************************************************

	subroutine basin_read_by_unit(iunit)

! reads basin or fails if error

	integer iunit

	integer nk,ne

	call sp13_get_par(iunit,nk,ne,ngr,mbw)
	call basin_init(nk,ne)			!here we set nkn, nel
	rewind(iunit)
	call sp13rr(iunit,nkn,nel)
	bbasinread = .true.
	!write(6,*) 'finished basin_read (module)'

	end subroutine basin_read_by_unit

c***********************************************************

	subroutine basin_set_read_basin(bread)

! sets flag if basin has been read

	logical bread

	bbasinread = bread

	end subroutine basin_set_read_basin

c***********************************************************

	function basin_has_read_basin()

! checks if basin has been read

	logical basin_has_read_basin

	basin_has_read_basin = bbasinread

	end function basin_has_read_basin

c***********************************************************

	function basin_has_basin()

! checks if basin is available (not necessarily read)

	logical basin_has_basin

	basin_has_basin = nkn_basin > 0 .and. nel_basin > 0

	end function basin_has_basin

c***********************************************************

	subroutine basin_get_dimension(nk,ne)

! returns dimension of arrays (already allocated)

	integer nk,ne

	nk = nkndi
	ne = neldi

	end subroutine basin_get_dimension

c***********************************************************

	function basin_is_basin_by_unit(iunit)

	logical basin_is_basin_by_unit
	integer iunit

	integer nvers

	call sp13test(iunit,nvers)

	basin_is_basin_by_unit = nvers > 0

	end function basin_is_basin_by_unit

c***********************************************************

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

!==================================================================
        end module basin
!==================================================================

        subroutine basin_check(text)

        use basin

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

c***********************************************************

	subroutine sp13rr(nb,nknddi,nelddi)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	use basin

	implicit none

	integer nb,nknddi,nelddi

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
	write(6,*) 'nknddi,nelddi :',nknddi,nelddi
	write(6,*) 'nkn,nel       :',nkn,nel
	write(6,*) 'ngr,mbw       :',ngr,mbw
	stop 'error stop sp13rr: dimension error'
	end

c***********************************************************

	subroutine sp13uw(nb)

c unformatted write to lagoon file
c
c nb		unit number for write

	use basin

	implicit none

	integer nb

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

	use basin

	implicit none

	integer nvers,nb,n

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

	use basin

	implicit none

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

	use basin

	implicit none

	real dcor,dirn

	dcor = dcorbas
	dirn = dirnbas

	end

c*************************************************

	subroutine bas_get_para(nkna,nela,ngra,mbwa)

	use basin

	implicit none

	integer nkna,nela,ngra,mbwa

	nkna = nkn
	nela = nel
	ngra = ngr
	mbwa = mbw

	end

c*************************************************

	subroutine bas_get_minmax(xmin,ymin,xmax,ymax)

	use basin

	implicit none

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

        function is_depth_unique()

	use basin

	implicit none

        logical is_depth_unique

	integer ie,ii,k
	real h
        real :: flag = -999.
        real haux(nkn)

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

c*************************************************
c*************************************************
c*************************************************

	subroutine bas_insert_irregular(nx,ny,xx,yy)

	use basin

c inserts irregular (but structured) basin (boxes) into basin structure

	implicit none

	integer nx,ny
	real xx(nx,ny)
	real yy(nx,ny)

	logical bdebug
	integer ix,iy
	integer k
	real xm,ym
	real regpar(7)

	regpar = 0.		!this is wrong for coordinates - corrected later
        regpar(1) = nx
        regpar(2) = ny

	call bas_insert_regular(regpar)		!index is right now

	k = 0
        do iy=1,ny
          do ix=1,nx
	    k = k + 1
	    xgv(k) = xx(ix,iy)
	    ygv(k) = yy(ix,iy)
          end do
        end do

        do iy=2,ny
          do ix=2,nx
	    k = k + 1
	    xm = xx(ix,iy) + xx(ix-1,iy) + xx(ix,iy-1) + xx(ix-1,iy-1)
	    ym = yy(ix,iy) + yy(ix-1,iy) + yy(ix,iy-1) + yy(ix-1,iy-1)
	    xgv(k) = xm/4.
	    ygv(k) = ym/4.
          end do
        end do

	end

c*************************************************

	subroutine bas_insert_regular(regpar)

	use basin

c inserts regular basin (boxes) into basin structure

	implicit none

	real regpar(7)

	logical bdebug
	integer nx,ny,ix,iy
	integer k,ie,nk,ne
	integer k1,k2,k3,k4,km
	real x0,y0,dx,dy
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

	bdebug = .false.
	if( bdebug ) then
	  write(6,*) 'regular basin inserted: ',nk,ne
	  write(6,*) nx,ny,x0,y0,dx,dy
	  write(6,*) minval(xgv),maxval(xgv)
	  write(6,*) minval(ygv),maxval(ygv)
	end if

	end

c*************************************************
c*************************************************
c*************************************************

        subroutine estimate_ngr(ngrade)

c estimates grade of basin - estimate is exact

	use basin

        implicit none

	integer ngrade

        integer ng(nkn)

	call compute_ng(ngrade,ng)

	end

c*************************************************

        subroutine compute_ng(ngrade,ng)

c computes grade of basin and nuber of grades per node

	use basin

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

c*************************************************

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

c*************************************************
c*************************************************
c*************************************************

	subroutine bas_get_node_coordinates(xp,yp)

	use basin

	implicit none

	real xp(nkn)
	real yp(nkn)

	xp = xgv
	yp = ygv

	end

c*************************************************

	subroutine bas_get_elem_coordinates(xp,yp)

        use basin

        implicit none

        real xp(nel),yp(nel)

        integer ie,ii,k
        double precision x,y

        do ie=1,nel
          x = 0.
          y = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            x = x + xgv(k)
            y = y + ygv(k)
          end do
          xp(ie) = x / 3.
          yp(ie) = y / 3.
        end do

        end

c*************************************************

	subroutine bas_get_special_coordinates(np,nodes,xp,yp)

	use basin

	implicit none

	integer np
	integer nodes(np)
	real xp(nkn)
	real yp(nkn)

	integer i,k

	do i=1,np
	  k = nodes(i)
	  xp(i) = xgv(k)
	  yp(i) = ygv(k)
	end do

	end

c*************************************************

