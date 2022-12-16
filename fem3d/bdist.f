
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2005,2009-2010,2013-2015,2013-2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

c distance routines
c
c contents :
c
c subroutine shdist(rdist)
c                       makes distance array from open boundaries
c subroutine mkdist(nkn,idist,rdist)
c                       makes distance array from given nodes
c
c revision log :
c
c 20.08.2003	ggu	new routine wr2dn()
c 05.01.2005	ggu	routines for writing nos file into subnsoa.f
c 07.01.2005	ggu	documentation for shdist
c 28.04.2009	ggu	links re-structured
c 23.03.2010	ggu	changed v6.1.1
c 13.12.2013	ggu	make local distances using boundary parameter nad
c 28.01.2014	ggu	changed VERS_6_1_71
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_52
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 09.10.2017	ggu	changed VERS_7_5_33
c 05.12.2017	ggu	changed VERS_7_5_39
c 06.07.2018	ggu	changed VERS_7_5_48
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 26.05.2020	ggu	rdist is now defined on elements
c 21.03.2022	ggu	disable writing of dist.shy
c 28.05.2022	ggu	better error handling in mkdist_new()
c 12.12.2022	ggu	new routine mkdist_mpi() -> results will change
c 16.12.2022	ggu	nadist now running with mpi
c
c****************************************************************

        subroutine shdist(redist)

c makes distance array from open boundaries
c
c nad of single boundaries can be either 0 or nadist global
c if nadist is not given, it will be set to maximum of nad on single boundaries
c
c rdist is contained between 0 and 1
c rdist==1 means computation of all terms, rdist=0 excludes non-linear terms
c
c if nadist (=d) is not given rdist = 1 (default)


	use basin
	use shympi

	implicit none

        real redist(nel)

c local variables

	logical, parameter :: bwrite = .false.		! write bdist.grd file
        integer i,k,kk,ie,ii,kext
        integer nadist,nad,nadmax
        integer ibc,n,itype,nk
	integer nbc
	integer ivar
	integer, allocatable :: nads(:)
        integer, allocatable :: idist(:)
        integer, allocatable :: ieaux(:)
        real, allocatable :: rdist(:)

	real r

	integer iapini,ipint,ipext
        integer nbnds,itybnd,nkbnds,kbnds
        real getpar

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

        redist = 1.

        nadist = nint(getpar('nadist'))		!global value

	nbc = nbnds()
	allocate(nads(nbc))
	do ibc=1,nbc
	  call get_bnd_ipar(ibc,'nad',nad)	!local value
	  if( nad < 0 ) nad = nadist
	  nads(ibc) = nad
	end do
	nadmax = maxval(nads)			!maximum of local values

	if( nadmax /= nadist ) then
	  if( nadmax > 0 .and. nadist > 0 ) goto 99
	  if( nadist == 0 ) nadist = nadmax
	end if

	do ibc=1,nbc
	  nad = nads(ibc)
	  if( nad > 0 .and. nad /= nadmax ) goto 99
	end do

	!write(6,*) my_id,nadist,nadmax
	call shympi_barrier
	write(6,*) 'using nadist = ',nadist
	if( nadist == 0 ) return	!noithing to be done

	allocate(idist(nkn))
	allocate(rdist(nkn))

        idist = 0
        rdist = 1.

c-----------------------------------------------------------------
c run over open boundary nodes and set idist
c-----------------------------------------------------------------

        do ibc=1,nbc
          itype = itybnd(ibc)
          if( itype .eq. 1 .or. itype .eq. 2 ) then
            nk = nkbnds(ibc)
	    nad = nads(ibc)
	    if( nad .gt. 0 ) then
              do i=1,nk
                k = kbnds(ibc,i)
		kext = ipext(k)
                if( k > 0 ) idist(k) = nad
              end do
	    end if
          end if
        end do

c	now idist is either 0 or nad if on boundary and nadist is set

c-----------------------------------------------------------------
c make distance
c-----------------------------------------------------------------

        write(6,*) 'Making distance rdist'
        !call mkdist(nadist,nkn,idist,rdist)		!old version
        !call mkdist_new(nkn,idist,rdist)		!new non-mpi version
        call mkdist_mpi(idist,rdist)			!mpi-version

c-----------------------------------------------------------------
c compute value of rdist on elements
c-----------------------------------------------------------------

	do ie=1,nel
	  r = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    r = r + rdist(k)
	  end do
	  redist(ie) = r / 3.
	end do

c-----------------------------------------------------------------
c write dist (nos) file
c-----------------------------------------------------------------
 
	if( bwrite ) then
	  allocate(ieaux(nel))
	  ieaux = 0
	  call write_grd_general('bdist.grd','distance from boundary'
     +				,idist,ieaux,rdist,redist)
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   99	continue
	write(6,*) 'incongruence in nad values of boundary'
	write(6,*) 'nadist = ',nadist
	write(6,*) 'nad values on boundaries: ',nbc
	write(6,*) nads
	stop 'error stop shdist: incongruent nadist and local nad values'
        end

c*************************************************************************** 
                         
        subroutine mkdist(nadist,nkn,idist,rdist)

c makes distance array from given nodes
c
c can olny deal with one global value

	use mod_geom

        implicit none

        integer nadist
        integer nkn
        integer idist(nkn)
        real rdist(nkn)

        integer k,kk,i
        integer n
        integer idact,idnew,nfound
	integer nodes(maxlnk)
        real r,d,d2

c----------------------------------------------------------
c initialize with first level
c----------------------------------------------------------

	stop 'error stop mkdist: please use mkdist_new'

c       idist is already initialized with first level

c----------------------------------------------------------
c loop on levels
c----------------------------------------------------------

        idact = 1
        nfound = nkn

        do while( nfound .gt. 0 )

          idnew = idact + 1
          nfound = 0

          do k=1,nkn
            if( idist(k) .eq. idact ) then

	      call get_nodes_around(k,maxlnk,n,nodes)

              do i=1,n
                kk = nodes(i)
                if( idist(kk) .eq. 0 ) then
                  idist(kk) = idnew
                  nfound = nfound + 1
                end if
              end do
            end if
          end do

          idact = idnew

        end do

c----------------------------------------------------------
c set real value
c----------------------------------------------------------

        do k=1,nkn
          rdist(k) = idist(k)
        end do

c----------------------------------------------------------
c adjust distance
c----------------------------------------------------------

        d = nadist
        d2 = d / 2.

        do k=1,nkn
          r = rdist(k)
	  if( r .lt. 0. .or. r .gt. nkn ) stop 'error stop: internal (1)'
          r = r - d2            !first rows no adv terms
          r = r / d             !slowly introduce
          r = max(0.,r)
          r = min(1.,r)
          !if( rdist(k) .eq. 1 ) write(6,*) 'bdist ',rdist(k),d,d2,r
          rdist(k) = r
        end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

        end

c*************************************************************************** 
                         
        subroutine mkdist_new(nkn,idist,rdist)

c makes distance array from given nodes
c
c d = nadist or nad
c d2=d/2 rows of nodes have rdist = 0
c and then the next ones have rdist = i/d with i = d2+1, d2+d
c the rest has again rdist = 1
c
c example: nadist = d = 4,   d2 = d/2 = 2
c
c   row i:   1   2   3   4   5   6   7   8   ...
c   rdist:   0   0  1/4 2/4 3/4  1   1   1   ...
c
c routine is too complicated, and cannot deal with mpi

	use mod_geom
	use shympi

        implicit none

        integer nkn
        integer idist(nkn)
        real rdist(nkn)

	logical bdebug
        integer k,kk,ka,i,ic
        integer n,na,nanew
        integer idact,idnew,nfound
	integer nodes(maxlnk)
	real r

c----------------------------------------------------------
c initialize with first level
c----------------------------------------------------------

	ic = count(idist>0)
	ic = shympi_sum(ic)

	if( ic > 0 .and. shympi_is_parallel() ) then
	  write(6,*) 'cannot yet handle nadist or nad different from 0'
	  stop 'error stop mkdist_new: cannot run with mpi'
	end if

	do k=1,nkn
	  rdist(k) = -1.
	end do

c----------------------------------------------------------
c loop on levels
c----------------------------------------------------------

	do k=1,nkn
	  if( idist(k) .gt. 0 .and. rdist(k) .eq. -1. ) then
	    nanew = idist(k)
	    na = nanew
	    do ka=1,nkn		!mark first row
	      if( idist(ka) .eq. nanew 
     +			.and. rdist(ka) .eq. -1. ) then
	         rdist(ka) = 0.
	      end if
	    end do
	    do while( na .gt. 0 )
	      do ka=1,nkn
                if( idist(ka) .eq. na .and. rdist(ka) .ge. 0. ) then
	          call get_nodes_around(ka,maxlnk,n,nodes)
                  do i=1,n
                    kk = nodes(i)
		    r = 1. - (na-1)/float(nanew)
		    r = int(r)
		    r = max(0.,r)
		    if( rdist(kk) .ge. 0. ) r = min(rdist(kk),r)
		    idist(kk) = max(idist(kk),na-1)
		    rdist(kk) = r
		  end do
	        end if
	      end do
	      na = na - 1
	    end do
	  end if
	end do

c----------------------------------------------------------
c set real value
c----------------------------------------------------------

        do k=1,nkn
          if( rdist(k) .lt. 0. ) rdist(k) = 1.
        end do

	na = 0
        do k=1,nkn
	  r = rdist(k)
          if( r .lt. 0. .or. r .gt. 1. ) then
	    write(6,*) 'error in mdist_new: ',k,r
	    na = na + 1
	  end if
        end do

	if( na .gt. 0 ) then
	  write(6,*) 'errors = ',na
	  stop 'error stop mdist_new: values out of bound'
	end if

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

        end

c*************************************************************************** 

        subroutine mkdist_mpi(idist,rdist)

c makes distance array from given nodes
c
c can deal with global and local values
c
c idist is 0 or nad if on boundary and nadist is set
c first set idist over all nodes
c then convert idist to rdist

	use basin
	use mod_geom
	use shympi

        implicit none

        integer idist(nkn)
        real rdist(nkn)

	logical bdebug
        integer k,kk,ka,i,ic,ie,ii
        integer n,na,nanew,nadmax,nchange,nmax,nmin
	integer id,nad2
        integer idact,idnew,nfound,iloop
	integer nodes(maxlnk)
	!integer iaux(nkn)
	integer, allocatable :: iaux(:)
	real frac

c----------------------------------------------------------
c initialize with first level
c----------------------------------------------------------

	allocate(iaux(nkn))

	ic = count(idist>0)
	ic = shympi_sum(ic)

	nadmax = maxval(idist)
	nadmax = shympi_max(nadmax)
	!write(6,*) 'using nadmax = ',nadmax

c----------------------------------------------------------
c loop on levels and set idist
c----------------------------------------------------------

	call shympi_barrier

	iloop = 0

	do
	  iaux = idist
	  nchange = 0
	  do ie=1,nel
	    nmax = 0
	    nmin = nadmax
	    do ii=1,3
	      k = nen3v(ii,ie)
	      nmax = max(nmax,idist(k))
	      nmin = min(nmin,idist(k))
	    end do
	    if( nmin == 0 .and. nmax > 1 ) then
	      do ii=1,3
	        k = nen3v(ii,ie)
		if( idist(k) == 0 ) iaux(k) = max(iaux(k),nmax-1)
	      end do
	      nchange = nchange + 1
	    end if
	  end do
	  idist = iaux
	  iloop = iloop + 1
	  !write(6,*) 'changes computing idist: ',my_id,iloop,nchange
	  if( iloop >= nadmax ) exit
	  call shympi_exchange_2d_node(idist)
	end do

	if( nchange /= 0 ) then
	  stop 'error stop mkdist_mpi: internal error (1)'
	end if

	call shympi_barrier

c----------------------------------------------------------
c set real value rdist
c----------------------------------------------------------

	rdist = 0
	nad2 = (nadmax+1) / 2
	frac = 1./(nad2+1)

	do k=1,nkn
	  id = idist(k)
	  if( id == 0 ) then
	    rdist(k) = 0.
	  else if( id > nad2 ) then
	    rdist(k) = 1.
	  else
	    rdist(k) = id * frac
	  end if
	end do

	rdist = 1. - rdist

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

        end

c*************************************************************************** 

