
c******************************************************************

	subroutine poisson_init

	use mod_system

	implicit none

	bsys3d = .true.
	bsys3d = .false.

	end

c******************************************************************

	subroutine poisson_2d

	use basin

	implicit none

	real pvar(nkn)

	call poisson_2d_set_obc(pvar)
	call poisson_2d_solve(pvar)
	call poisson_2d_write(pvar)

	stop 'end of poisson case'
	end

c******************************************************************

	subroutine poisson_2d_solve(pvar)

	use basin

	implicit none

	real pvar(nkn)

	call system_init

	call poisson_2d_assemble(pvar)

	call system_solve_z(nkn,pvar) 	!solves system matrix for pvar
	call system_adjust_z(nkn,pvar)	!copies solution to new pvar

	end

c******************************************************************

	subroutine poisson_2d_write(pvar)

	use basin

	implicit none

	real pvar(nkn)

	integer k
	integer, save :: iu = 0

        do k = 1,nkn
          if ( xgv(k) .eq. 0.0 ) then
            write(635,*) ygv(k),pvar(k)
          end if
        end do

	call conwrite(iu,'.p2d',1,10,1,pvar)

	end

c******************************************************************

	subroutine poisson_2d_set_obc(pvar)

	use basin

	implicit none

	real pvar(nkn)

	include 'mkonst.h'

	integer k

	pvar = flag

	do k=1,nkn
          if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
	    pvar(k) = 10.
          !else if( ipv(k) .ge. 3601 .and. ipv(k) .le. 3609 ) then
          else if( ipv(k) .ge. 217 .and. ipv(k) .le. 225 ) then
	    pvar(k) = 20.
	  end if
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine poisson_2d_assemble(pvar)

c assembles linear system matrix
c
c vqv		flux boundary condition vector
c
c semi-implicit scheme for 3d model
c
c written on 18.02.91 by ggu  (from scratch)
c changed on 04.06.91 by ggu  (c=(1) : friction term has been corrected)
c changed on 01.10.92 by ggu  (staggered FE - completely restructured)
c 12.01.2001    ggu     solve for znv and not level difference (ZNEW)

	use mod_internal
	use mod_depth
	use mod_area !WILL DEB
	use evgeom
	use levels
	use basin

	implicit none

	real pvar(nkn)

	real drittl
	parameter (drittl=1./3.)

	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'
 
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk,l,k
	integer ngl
	integer ilevel
	integer lmax
	real aj,aj4,aj12
	real ht
	real h11,hh999
	real delta
	real hia(3,3),hik(3),amatr(3,3)
	real b(3),c(3),z(3)

	integer locsps,loclp,iround
	real getpar

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------

	do ie=1,nel

c	  ------------------------------------------------------
c	  initialize element values
c	  ------------------------------------------------------

	  aj=ev(10,ie)
          aj4=4.*aj
          aj12=12.*aj
	  do i=1,3
	    kk=nen3v(i,ie)
	    kn(i)=kk
	    b(i)=ev(i+3,ie)
	    c(i)=ev(i+6,ie)
	  end do

c	  ------------------------------------------------------
c	  set element matrix and RHS
c	  ------------------------------------------------------

	  do n=1,3
	    do m=1,3
	      hia(n,m) = -aj12 * ( b(n) * b(m) + c(n) * c(m) )
	    end do
	    hik(n) = 0.
	  end do

c	  ------------------------------------------------------
c	  boundary conditions
c	  ------------------------------------------------------

	  do i=1,3
	    if( pvar(kn(i)) .ne. flag ) then
	      j1=mod(i,3)+1
	      j2=mod(i+1,3)+1
              hia(i,i)=1.
              hia(i,j1)=0.
              hia(i,j2)=0.
              hik(i)=pvar(kn(i))
	    end if
	  end do

c	  ------------------------------------------------------
c	  in hia(i,j),hik(i),i,j=1,3 is system
c	  ------------------------------------------------------

	  call system_assemble(ie,nkn,mbw,kn,hia,hik)

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

