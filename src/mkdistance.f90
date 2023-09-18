!
! $Id: bdist.f,v 1.5 2009-05-21 09:24:00 georg Exp $
!
! distance routines
!
! contents :
!
! subroutine shdist(rdist)
!                       makes distance array from open boundaries
! subroutine mkdist(nkn,idist,rdist)
!                       makes distance array from given nodes
!
! revision log :
!
! 20.08.2003    ggu     new routine wr2dn()
! 05.01.2005    ggu     routines for writing nos file into subnsoa.f
! 07.01.2005    ggu     documentation for shdist
! 28.04.2009    ggu     links re-structured
! 13.12.2013    ggu     make local distances using boundary parameter nad
!
!****************************************************************
!-------------------------------------------------------------------------------
        module mkdistance
!-------------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------------

        subroutine shdist(rdist)

! makes distance array from open boundaries
!
! rdist is contained between 0 and 1
! if nadist (=d) is not given rdist = 1 (default)
! otherwise the first d2=d/2 rows of nodes have rdist = 0
! and then the next ones have rdist = i/d with i = d2+1, d2+d
! the rest has again rdist = 1
!
! example: nadist = d = 4,   d2 = d/2 = 2
!
!   row i:   1   2   3   4   5   6   7   8   ...
!   rdist:   0   0  1/4 2/4 3/4  1   1   1   ...

	use basin
        use shympi
        use mpi_io_admin
        use para
        use bnd_admin
        use fem_util
        use nos_util
        use lnku

	implicit none

        double precision rdist(nkn)

! local variables

        integer idist(nkn)

        integer i,k,kk
        integer nadist,nad
        integer ibc,n,itype,nk
	integer nbc

!-----------------------------------------------------------------
! get parameters
!-----------------------------------------------------------------

        do k=1,nkn
          rdist(k) = 1.
          idist(k) = 0
        end do

        nadist = nint(getpar('nadist'))		!global value

!-----------------------------------------------------------------
! gather open boundary nodes
!-----------------------------------------------------------------

        n = 0

	nbc = nbnds()

        do ibc=1,nbc
          itype = itybnd(ibc)
          if( itype .eq. 1 .or. itype .eq. 2 ) then
            nk = nkbnds(ibc)
	    call get_bnd_ipar(ibc,'nad',nad)	!local value
	    if( nad .lt. 0 ) nad = nadist
	    if( nad .gt. 0 ) then
              do i=1,nk
                k = kbnds(ibc,i)
                idist(k) = 1			!old version
                idist(k) = nad
              end do
	    end if
          end if
        end do

!-----------------------------------------------------------------
! make distance
!-----------------------------------------------------------------

        write(6,*) 'Making distance rdist'
        !call mkdist(nadist,nkn,idist,rdist)		!old version
        call mkdist_new(nkn,idist,rdist)

!-----------------------------------------------------------------
! write dist (nos) file
!-----------------------------------------------------------------

        call wrnos2d('dist','distance from boundary nodes',rdist)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        end

!*************************************************************************** 
                         
        subroutine mkdist(nadist,nkn,idist,rdist)

! makes distance array from given nodes
!
! can olny deal with one global value
!
! rdist of open boundary nodes is 1
! other nodes are > 1 (integer)
! example: neibors of rdist=1 nodes have rdist=2 etc.

	use geom
        use lnku

        implicit none

        integer nadist
        integer nkn
        integer idist(nkn)
        double precision rdist(nkn)

        integer k,kk,i
        integer n
        integer idact,idnew,nfound
	integer nodes(maxlnk)
        double precision r,d,d2

!----------------------------------------------------------
! initialize with first level
!----------------------------------------------------------

	stop 'error stop mkdist: please use mkdist_new'

!       idist is already initialized with first level

!----------------------------------------------------------
! loop on levels
!----------------------------------------------------------

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

!----------------------------------------------------------
! set double precision value
!----------------------------------------------------------

        do k=1,nkn
          rdist(k) = idist(k)
        end do

!----------------------------------------------------------
! adjust distance
!----------------------------------------------------------

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

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

        end

!*************************************************************************** 
                         
        subroutine mkdist_new(nkn,idist,rdist)

! makes distance array from given nodes
!
! can deal with global and local values
!
! rdist of open boundary nodes is 1
! other nodes are > 1 (integer)
! example: neibors of rdist=1 nodes have rdist=2 etc.

	use geom
	use shympi
        use lnku

        implicit none

        integer nkn
        integer idist(nkn)
        double precision rdist(nkn)

	logical bdebug
        integer k,kk,ka,i
        integer n,na,nanew
        integer idact,idnew,nfound
	integer nodes(maxlnk)
	double precision r

!----------------------------------------------------------
! initialize with first level
!----------------------------------------------------------

	if( any(idist>0) .and. shympi_is_parallel() ) then
	  stop 'error stop mkdist_new: cannot run with mpi'
	end if

	do k=1,nkn
	  rdist(k) = -1.
	end do

!----------------------------------------------------------
! loop on levels
!----------------------------------------------------------

	do k=1,nkn
	  if( idist(k) .gt. 0 .and. rdist(k) .eq. -1. ) then
	    nanew = idist(k)
	    na = nanew + nanew/2
	    na = nanew
	    do ka=1,nkn		!mark first row
	      if( idist(ka) .eq. nanew .and. rdist(ka) .eq. -1. ) then
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

!----------------------------------------------------------
! set double precision value
!----------------------------------------------------------

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

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

        end

!*************************************************************************** 

!-------------------------------------------------------------------------------
        end module mkdistance
!-------------------------------------------------------------------------------
