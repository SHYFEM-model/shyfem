
c******************************************************************

	subroutine poisson_init

	use mod_system

	implicit none

	integer ipoiss
	real getpar

	ipoiss = nint(getpar('ipoiss'))

	bsys3d = .false.
	if( ipoiss /= 0 ) bsys3d = .true.

	end

c******************************************************************

	subroutine poisson_compute

	use basin

	implicit none

	integer ipoiss
	real getpar

	ipoiss = nint(getpar('ipoiss'))

	if( ipoiss == 0 ) return

	call poisson_2d
	call poisson_3d

	stop 'end of poisson case'
	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine poisson_2d

	use basin

	implicit none

	real pvar(nkn)

	call poisson_set_obc(1,pvar)
	call poisson_2d_solve(pvar)
	call poisson_write(1,pvar)

	end

c******************************************************************

	subroutine poisson_2d_solve(pvar)

	use basin

	implicit none

	real pvar(nkn)
	real p0(nkn)

	p0 = 0.

	call system_init

	call poisson_2d_assemble(pvar)

	call system_solve_z(nkn,p0) 	!solves system matrix for pvar
	call system_adjust_z(nkn,pvar)	!copies solution to new pvar

	end

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

	  !call system_assemble(ie,nkn,mbw,kn,hia,hik)
	  call system_assemble(ie,kn,hia,hik)

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine poisson_3d

	use basin
	use levels

	implicit none

	real pvar(nlvdi,nkn)

	call poisson_set_obc(nlvdi,pvar)
	call poisson_3d_solve(nlvdi,pvar)
	call poisson_write(nlvdi,pvar)

	end

c******************************************************************

	subroutine poisson_3d_solve(nlvdi,pvar)

	use basin
	use levels, only : nlv

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)
	real p0(nlvdi,nkn)

	p0 = 0.

	call system_init

	call poisson_3d_assemble(pvar)

	call system_solve_3d(nkn,nlvdi,nlv,p0) 	  !solves system matrix for pvar
	call system_adjust_3d(nkn,nlvdi,nlv,pvar) !copies solution to new pvar

	end

c******************************************************************

	subroutine poisson_3d_assemble(pvar)

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

	real pvar(nlvdi,nkn)

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
	real hia3d(-1:+1,3,3)
	real b(3),c(3),z(3)

	real hd,hdm,hdp,rhm,rhp,rhc
	real hldaux(0:nlv+1)

	integer locsps,loclp,iround
	real getpar

	hldaux = 0.
	hldaux(1:nlv) = hldv(:)

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

	  lmax = ilhv(ie)

	  do l=1,lmax

	  hd = hldaux(l)
	  hdm = hldaux(l-1)
	  hdp = hldaux(l+1)

	  rhm = 2. / ( hd + hdm )
	  if( l == 1 ) rhm = 0.
	  rhp = 2. / ( hd + hdp )
	  if( l == lmax ) rhp = 0.
	  rhc = rhm + rhp

	  do n=1,3
	    do m=1,3
	      hia3d(0,n,m) = -aj12 * hd * ( b(n) * b(m) + c(n) * c(m) )
	      hia3d(0,n,m) = hia3d(0,n,m) - aj4 * rhc
	      hia3d(-1,n,m) = aj4 * rhm
	      hia3d(+1,n,m) = aj4 * rhp
	    end do
	    hik(n) = 0.
	  end do

c	  ------------------------------------------------------
c	  boundary conditions
c	  ------------------------------------------------------

	  do i=1,3
	    if( pvar(l,kn(i)) .ne. flag ) then
	      j1=mod(i,3)+1
	      j2=mod(i+1,3)+1
              hia3d(0,i,i)=1.
              hia3d(0,i,j1)=0.
              hia3d(0,i,j2)=0.
              hik(i)=pvar(l,kn(i))
	    end if
	  end do

c	  ------------------------------------------------------
c	  in hia(i,j),hik(i),i,j=1,3 is system
c	  ------------------------------------------------------

	  !call system_assemble(ie,nkn,mbw,kn,hia,hik)
	  call system_assemble_3d(ie,l,nlv,kn,hia3d,hik)

	  end do

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine poisson_set_obc(nlvdi,pvar)

	use basin

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)

	include 'mkonst.h'

	integer k
	integer ipoiss
	real bnd
	real getpar

	ipoiss = nint(getpar('ipoiss'))

	pvar = flag

	do k=1,nkn
	  bnd = 0.
	  if( ipoiss == 1 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
		bnd = 10.
            else if( ipv(k) .ge. 3601 .and. ipv(k) .le. 3609 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 2 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 9 ) then
		bnd = 10.
            else if( ipv(k) .ge. 217 .and. ipv(k) .le. 225 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 3 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 2 ) then
		bnd = 10.
            else if( ipv(k) .ge. 5 .and. ipv(k) .le. 6 ) then
		bnd = 20.
	    end if
	  else if( ipoiss == 4 ) then
            if( ipv(k) .ge. 1 .and. ipv(k) .le. 2 ) then
		bnd = 10.
            else if( ipv(k) .ge. 7 .and. ipv(k) .le. 8 ) then
		bnd = 20.
	    end if
	  else
	    write(6,*) 'ipoiss = ',ipoiss
	    stop 'error stop poisson_set_obc'
	  end if
	  if( bnd /= 0. ) pvar(:,k) = bnd
	end do
	write(6,*) 'poisson boundary: ',pvar

	end

c******************************************************************

	subroutine poisson_write(nlvdi,pvar)

	use basin

	implicit none

	integer nlvdi
	real pvar(nlvdi,nkn)

	integer k,l,iu
	integer, save :: iu2d = 0
	integer, save :: iu3d = 0

	l = (nlvdi+1)/2
	iu = 635
	if( nlvdi /= 1 ) iu = 666

	write(6,*) 'poisson writing: ',nlvdi,l,iu,nkn

        do k = 1,nkn
          if ( xgv(k) .eq. 0.0 ) then
            write(iu,*) ygv(k),pvar(l,k)
          end if
        end do
	flush(iu)

	if( nlvdi == 1 ) then
	  call conwrite(iu2d,'.p2d',1,10,nlvdi,pvar)
	  flush(iu2d)
	else
	  call conwrite(iu3d,'.p3d',1,10,nlvdi,pvar)
	  flush(iu3d)
	end if

	end

c******************************************************************
c******************************************************************
c******************************************************************

