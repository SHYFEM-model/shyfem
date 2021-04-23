
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997,2009-2011,2014-2020  Georg Umgiesser
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

c initialization routines
c
c contents :
c
c revision log :
c
c 31.05.1997	ggu	unnecessary routines deleted
c 27.06.1997	ggu	bas routines into own file
c 02.04.2009	ggu	error messages changed
c 23.03.2010	ggu	changed v6.1.1
c 12.01.2011	ggu	debug routine introduced (basin_test)
c 27.01.2011	ggu	changed VERS_6_1_17
c 01.03.2011	ggu	changed VERS_6_1_20
c 18.06.2014	ggu	changed VERS_6_1_77
c 23.10.2014	ggu	introduced ftype and nvers = 4
c 30.10.2014	ggu	changed VERS_7_0_4
c 23.12.2014	ggu	changed VERS_7_0_11
c 04.01.2015	ggu	new routine basin_get_par()
c 26.02.2015	ggu	changed VERS_7_1_5
c 31.03.2015	ggu	set iarnv on read
c 01.04.2015	ggu	changed VERS_7_1_7
c 25.05.2015	ggu	module introduced
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 02.10.2015	ggu	in basin_open_file eliminated double read (bug)
c 02.10.2015	ggu	new routines is_depth_unique(), estimate_ngr()
c 10.10.2015	ggu	changed VERS_7_3_2
c 12.10.2015	ggu	changed VERS_7_3_3
c 19.10.2015	ggu	changed VERS_7_3_6
c 01.05.2016	ggu	new routines basin_has_basin()
c 20.05.2016	ggu	estimate_ngr() returns exact ngr
c 10.06.2016	ggu	new routine for inserting regular grid
c 14.06.2016	ggu	changed VERS_7_5_14
c 23.09.2016	ggu	new routines to check if basin has been read
c 30.09.2016	ggu	changed VERS_7_5_18
c 13.02.2017	ggu	changed VERS_7_5_23
c 09.05.2017	ggu	changed VERS_7_5_26
c 25.05.2017	ggu	changed VERS_7_5_28
c 13.06.2017	ggu	changed VERS_7_5_29
c 11.07.2017	ggu	changed VERS_7_5_30
c 14.11.2017	ggu	changed VERS_7_5_36
c 24.01.2018	ggu	changed VERS_7_5_41
c 22.02.2018	ggu	changed VERS_7_5_42
c 16.02.2019	ggu	changed VERS_7_5_60
c 18.05.2020	ggu	new routine basin_info_partition()
c 17.04.2021	ggu	some better error handling
c
c***********************************************************
c***********************************************************
c***********************************************************

!==================================================================
        module basin
!==================================================================

        implicit none

	integer, private, parameter :: ftype = 789233567
	integer, private, parameter :: nversm = 5

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

        integer, private, save :: npart_node = 0
        integer, private, save :: npart_elem = 0
        integer, private, save, allocatable :: area_part_node(:)
        integer, private, save, allocatable :: area_part_elem(:)

        INTERFACE basin_read
        MODULE PROCEDURE basin_read_by_file, basin_read_by_unit
        END INTERFACE

        INTERFACE basin_write
        MODULE PROCEDURE basin_write_by_file, basin_write_by_unit
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
	  deallocate(area_part_node)
	  deallocate(area_part_elem)
	end if

	nkn = nk
	nel = ne
	nkndi = nk
	neldi = ne
	nkn_basin = nk
	nel_basin = ne

	npart_node = 0
	npart_elem = 0

	if( nk == 0 ) return

	allocate(nen3v(3,nel))
	allocate(ipev(nel))
	allocate(ipv(nkn))
	allocate(iarv(nel))
	allocate(iarnv(nkn))
	allocate(xgv(nkn))
	allocate(ygv(nkn))
	allocate(hm3v(3,nel))
	allocate(area_part_node(nkn))
	allocate(area_part_elem(nel))

	end subroutine basin_init

c***********************************************************

	function basin_open_file(file,status)

! opens file or returns 0
	
	integer basin_open_file
	character*(*) file
	character*(*), optional :: status

	integer iunit
	integer ifileo
	character*80 stat

	stat = 'old'
	if( present(status) ) stat = status

	basin_open_file = 0

	iunit = ifileo(0,file,'unform',stat)
	if( iunit .le. 0 ) return

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

	!write(6,*) 'finished reading basin: ',trim(file)

	end subroutine basin_read_by_file

c***********************************************************

	subroutine basin_read_by_unit(iunit)

! reads basin or fails if error

	integer iunit

	integer nk,ne

	call basin_get_par(iunit,nk,ne,ngr,mbw)
	call basin_init(nk,ne)			!here we set nkn, nel
	call basin_read_internal(iunit,nkn,nel)
	bbasinread = .true.
	!write(6,*) 'finished basin_read (module)'

	end subroutine basin_read_by_unit

c***********************************************************

	subroutine basin_write_by_file(file)

! writes basin

	character*(*) file

	integer iunit

	iunit = basin_open_file(file,'unknown')

	if( iunit .le. 0 ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop basin_write_by_file: cannot open file'
	end if

	call basin_write_by_unit(iunit)

	close(iunit)

	end subroutine basin_write_by_file

c***********************************************************

	subroutine basin_write_by_unit(iunit)

! writes basin

	integer iunit

	call basin_write_internal(iunit)

	end subroutine basin_write_by_unit

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

	call basin_test(iunit,nvers)

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

c***********************************************************

	subroutine basin_read_internal(nb,nknddi,nelddi)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	integer nb,nknddi,nelddi

	integer i,ii,nvers

	call basin_test(nb,nvers)

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

	iarnv = 0
	npart_node = 0
	npart_elem = 0
	area_part_node = 0
	area_part_elem = 0

	if( nvers < 5 ) return

	read(nb) npart_node,npart_elem
	if( npart_node > 0 ) read(nb)(area_part_node(i),i=1,nkn)
	if( npart_elem > 0 ) read(nb)(area_part_elem(i),i=1,nel)

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	stop 'error stop basin_read_internal: error reading file'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	write(6,*) 'nvers = ',-nvers
	stop 'error stop basin_read_internal: error in version'
   97	continue
	write(6,*) 'nknddi,nelddi :',nknddi,nelddi
	write(6,*) 'nkn,nel       :',nkn,nel
	write(6,*) 'ngr,mbw       :',ngr,mbw
	stop 'error stop basin_read_internal: dimension error'
	end subroutine basin_read_internal

c***********************************************************

	subroutine basin_write_internal(nb)

c unformatted write to lagoon file
c
c nb		unit number for write

	integer nb

	integer i,ii

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

	if( nversm < 5 ) return

	write(nb) npart_node,npart_elem
	if( npart_node > 0 ) write(nb)(area_part_node(i),i=1,nkn)
	if( npart_elem > 0 ) write(nb)(area_part_elem(i),i=1,nel)

	return
   99	continue
	write(6,*) 'Writing basin...'
	write(6,*) 'Cannot write bas file on unit :',nb
	stop 'error stop : basin_write_internal'
	end subroutine basin_write_internal

c***********************************************************

	subroutine basin_get_part_info(nn,ne)

	integer nn,ne

	nn = npart_node
	ne = npart_elem

	end subroutine basin_get_part_info

c***********************************************************

	subroutine basin_get_partition(nn,ne,nnp,nep,area_node,area_elem)

	integer nn,ne
	integer nnp,nep
	integer area_node(nn)
	integer area_elem(ne)

	nnp = npart_node
	nep = npart_elem
	area_node = area_part_node
	area_elem = area_part_elem

	end subroutine basin_get_partition

c***********************************************************

	subroutine basin_set_partition(nn,ne,nnp,nep,area_node,area_elem)

	integer nn,ne
	integer nnp,nep
	integer area_node(nn)
	integer area_elem(ne)

	integer i

	npart_node = nnp
	npart_elem = nep
	area_part_node = area_node
	area_part_elem = area_elem

	end subroutine basin_set_partition

c***********************************************************

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

c***********************************************************

	subroutine basin_test(nb,nvers)

c tests if file is BAS file
c
c nvers > 0 if file is BAS file

	implicit none

	integer nb	!unit number
	integer nvers	!version found (return) (<=0 if error or no BAS file)

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
	end subroutine

c***********************************************************

	subroutine basin_get_par(nb,nkn,nel,ngr,mbw)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	integer nb
	integer nkn,nel,ngr,mbw

	integer nvers
	character*80 file

	file = ' '
	if( nb > 0 ) inquire(nb,name=file)

	call basin_test(nb,nvers)

	if(nvers.eq.0) goto 99
	if(nvers.lt.0) goto 98

	read(nb) nkn,nel,ngr,mbw

	rewind(nb)

	return
   99	continue
	write(6,*) 'Cannot read bas file on unit :',nb
	if( nb > 0 ) write(6,*) 'file name = ',trim(file)
	stop 'error stop : basin_get_par'
   98	continue
	write(6,*) 'Cannot read version: nvers = ',-nvers
	if( nb > 0 ) write(6,*) 'file name = ',trim(file)
	stop 'error stop : basin_get_par'
   97	continue

	end subroutine

c***********************************************************

	subroutine basin_info_partition

	implicit none

	integer nnmax,nnmin,nemax,nemin
	integer k,ie,ip
	integer, allocatable :: nncount(:)
	integer, allocatable :: necount(:)

	nnmax = maxval(area_part_node)
	nnmin = minval(area_part_node)
	nemax = maxval(area_part_elem)
	nemin = minval(area_part_elem)

	write(6,*) 'information on partition:'
	write(6,*) 'nn min/max: ',nnmin,nnmax
	write(6,*) 'ne min/max: ',nemin,nemax

	allocate(nncount(nnmin:nnmax))
	allocate(necount(nemin:nemax))
	nncount = 0
	necount = 0

	do k=1,nkn
	  ip = area_part_node(k)
	  nncount(ip) = nncount(ip) + 1
	end do

	do ie=1,nel
	  ip = area_part_elem(ie)
	  necount(ip) = necount(ip) + 1
	end do

	write(6,*) 'partition on nodes:'
	do ip=nnmin,nnmax
	  write(6,*) ip,nncount(ip)
	end do

	write(6,*) 'partition on elems:'
	do ip=nemin,nemax
	  write(6,*) ip,necount(ip)
	end do

	if( sum(nncount) /= nkn .or. sum(necount) /= nel ) then
	  write(6,*) sum(nncount),nkn
	  write(6,*) sum(necount),nel
	  stop 'error stop basin_info_partition: internal error'
	end if
	  
	end

!==================================================================
        end module basin
!==================================================================

c***********************************************************

	subroutine bas_info

	use basin

	implicit none

	integer nnpart,nepart

	call basin_get_part_info(nnpart,nepart)

	if( descrr /= ' ' ) then
          write(6,*) trim(descrr)
	end if
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*) ' dcor = ',dcorbas,'  dirn = ',dirnbas
        write(6,*) ' nnpart = ',nnpart,'  nepart = ',nepart

	!call basin_info_partition

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
	    if( k <= 0 ) goto 98
	    if( k > nkn ) goto 98
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
	    if( k1 <= 0 .or. k2 <= 0 ) goto 99
	    if( k1 > nkn .or. k2 > nkn ) goto 99
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

	return
   98	continue
	write(6,*) '(98) nkn,k,ie: ',nkn,k,ie
	stop 'error stop compute_ng: corrupt node index'
   99	continue
	write(6,*) '(99) nkn,k1,k2,ie: ',nkn,k1,k2,ie
	stop 'error stop compute_ng: corrupt node index'
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

