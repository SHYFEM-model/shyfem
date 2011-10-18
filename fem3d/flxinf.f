c
c $Id: flxinf.f,v 1.10 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c 21.04.2010	ggu	comments, write all files
c
c notes :
c
c ppp		total flux through section
c ppppm(1,)	positive flux
c ppppm(2,)	negative flux (number is positive)
c
c	=> ppp = ppppm(1,) - ppppm(2,)
c
c*****************************************************************

	program flxinf

c reads flx files

	implicit none

	include 'param.h'

	integer nin
	integer idfile,nvers,nsect,kfluxm,nlmax
	integer i,it,idtflx
	integer j,l,lmax
	integer kflux(nfxdim)
	integer nlayers(nfxdim)
	real ppp(nfxdim)
	real ppppm(2,nfxdim)
	real fluxes(0:nlvdim,3,nfxdim)

	integer iapini,ideffi

c---------------------------------------------------------------

        if(iapini(2,1,1,0).eq.0) then
                stop 'error stop: iapini'
        end if

        nin=ideffi('datdir','runnam','.flx','unform','old')
        if(nin.le.0) then
	  stop 'error stop: cannot open file'
	end if

	read(nin) idfile,nvers
	read(nin) nsect,kfluxm,idtflx

	if( nsect .gt. nfxdim .or. kfluxm .gt. nfxdim ) then
	  write(6,*) idfile,nvers
	  write(6,*) nsect,kfluxm
	  write(6,*) nfxdim
	  stop 'error stop flxinf: nfxdim'
	end if

	read(nin) (kflux(i),i=1,kfluxm)

	nlmax = 0
	if( nvers .ge. 3 ) then
	  read(nin) nlmax,(nlayers(i),i=1,nsect)
	  if( nlmax .gt. nlvdim ) then
	    write(6,*) nlmax,nlvdim
	    stop 'error stop flxinf: nlvdim'
	  end if
	else
	  do i=1,nsect
	    nlayers(i) = 1
	  end do
	end if

	write(6,*) idfile,nvers,nsect,kfluxm,nlmax

c---------------------------------------------------------------
c loop on data
c---------------------------------------------------------------

	do while(.true.)

	  if( nvers .eq. 1 ) then
	    read(nin,end=2) it,nsect,(ppp(i),i=1,nsect)
	  else if( nvers .eq. 2 ) then
	    read(nin,end=2) it,nsect,(ppp(i),i=1,nsect)
     +				,(ppppm(1,i),ppppm(2,i),i=1,nsect)
	  else if( nvers .eq. 3 ) then
	    read(nin,end=2) it,nsect
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)

	  else
	    write(6,*) 'nvers = ',nvers
	    stop 'error stop flxinf: cannot handle version'
	  end if

	  if( nvers .ge. 3 ) then
	    do i=1,nsect
	      ppp(i) = fluxes(0,1,i)
	      ppppm(1,i) = fluxes(0,2,i)
	      ppppm(2,i) = fluxes(0,3,i)
	    end do
	  end if

c in ppppm(1,*) are positive fluxes
c in ppppm(2,*) are negative fluxes

c	  write(6,*) it,nsect
c	  write(6,*) (ppp(i),i=1,nsect)

	  write(6,'(i10,10i7)') it,(nint(ppp(i)),i=1,nsect)
c	  write(6,'(i10,2f12.4)') it,(ppp(i),i=1,nsect)
c     +				,(fluxes(1,1,i),i=1,nsect)
	  write(66,'(i10,20f10.2)') it,(ppp(i),i=1,nsect)	!real

	  write(67,'(i10,20f10.2)') it,(ppppm(1,i),i=1,nsect)	!positive
	  write(68,'(i10,20f10.2)') it,(ppppm(2,i),i=1,nsect)	!negative

          write(69,'(i10,20f10.2)') it,(ppppm(1,i)+ppppm(2,i)	!total (abs)
     +						,i=1,nsect)

	  write(65,*) it,nsect
	  do i=1,nsect
	    lmax = nlayers(i)
	    write(65,*) i,lmax
	    do l=0,lmax
	      write(65,*) l,(fluxes(l,j,i),j=1,3)
	    end do
	  end do

	end do

    2	continue

	close(nin)

	end
