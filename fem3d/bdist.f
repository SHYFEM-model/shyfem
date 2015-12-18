c
c $Id: bdist.f,v 1.5 2009-05-21 09:24:00 georg Exp $
c
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
c 20.08.2003    ggu     new routine wr2dn()
c 05.01.2005    ggu     routines for writing nos file into subnsoa.f
c 07.01.2005    ggu     documentation for shdist
c 28.04.2009    ggu     links re-structured
c 13.12.2013    ggu     make local distances using boundary parameter nad
c
c****************************************************************

        subroutine shdist(rdist)

c makes distance array from open boundaries
c
c rdist is contained between 0 and 1
c if nadist (=d) is not given rdist = 1 (default)
c otherwise the first d2=d/2 rows of nodes have rdist = 0
c and then the next ones have rdist = i/d with i = d2+1, d2+d
c the rest has again rdist = 1
c
c example: nadist = d = 4,   d2 = d/2 = 2
c
c   row i:   1   2   3   4   5   6   7   8   ...
c   rdist:   0   0  1/4 2/4 3/4  1   1   1   ...

	use basin

	implicit none

        real rdist(nkn)

c local variables

        integer idist(nkn)

        integer i,k,kk
        integer nadist,nad
        integer ibc,n,itype,nk
	integer nbc

	integer iapini,ipint
        integer nbnds,itybnd,nkbnds,kbnds
        real getpar

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

        do k=1,nkn
          rdist(k) = 1.
          idist(k) = 0
        end do

        nadist = nint(getpar('nadist'))		!global value

c-----------------------------------------------------------------
c gather open boundary nodes
c-----------------------------------------------------------------

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

c-----------------------------------------------------------------
c make distance
c-----------------------------------------------------------------

        write(6,*) 'Making distance rdist'
        !call mkdist(nadist,nkn,idist,rdist)		!old version
        call mkdist_new(nkn,idist,rdist)

c-----------------------------------------------------------------
c write dist (nos) file
c-----------------------------------------------------------------
 
        call wrnos2d('dist','distance from boundary nodes',rdist)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*************************************************************************** 
                         
        subroutine mkdist(nadist,nkn,idist,rdist)

c makes distance array from given nodes
c
c can olny deal with one global value
c
c rdist of open boundary nodes is 1
c other nodes are > 1 (integer)
c example: neibors of rdist=1 nodes have rdist=2 etc.

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
c can deal with global and local values
c
c rdist of open boundary nodes is 1
c other nodes are > 1 (integer)
c example: neibors of rdist=1 nodes have rdist=2 etc.

	use mod_geom

        implicit none

        integer nkn
        integer idist(nkn)
        real rdist(nkn)

	logical bdebug
        integer k,kk,ka,i,ks
        integer n,na,nanew
        integer idact,idnew,nfound
	integer nodes(maxlnk)
	real r

c----------------------------------------------------------
c initialize with first level
c----------------------------------------------------------

	do k=1,nkn
	  rdist(k) = -1.
	end do

c----------------------------------------------------------
c loop on levels
c----------------------------------------------------------

	ks = 12659

	do k=1,nkn
	  if( idist(k) .gt. 0 .and. rdist(k) .eq. -1. ) then
	    nanew = idist(k)
	    na = nanew + nanew/2
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

