c
c $Id: flxinf.f,v 1.10 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c 21.04.2010	ggu	comments, write all files
c 19.10.2011	ggu	new 3D write
c 07.04.2016	ggu	bug fix for allocation of arrays
c
c notes :
c
c this program extracts information from the FLX file and writes it to file
c
c call fluxes_2d(it,nlvdim,nsect,ivar,ptot,fluxes)
c
c writes 2d fluxes of water 
c this is mainly old stuff that might not be very interesting
c it has been kept mainly for compatibility
c
c call fluxes_3d(it,nlvdim,nsect,ivar,nlayers,fluxes)
c
c new 3d output of fluxes (water, salt, temp)
c the output is written to files depending on the value of ivar
c	ivar = 0	water	iunit = 160
c	ivar = 11	salt	iunit = 171
c	ivar = 12	temp	iunit = 172
c the format is:
c
c	time	nsect	ivar
c		1	lmax				!first section
c			0	ftot	fpos	fneg
c			1	ftot	fpos	fneg
c			...
c			lmax	ftot	fpos	fneg
c		2	lmax				!second section
c			0	ftot	fpos	fneg
c			1	ftot	fpos	fneg
c			...
c			lmax	ftot	fpos	fneg
c		...
c		nsect	lmax				!last section
c			0	ftot	fpos	fneg
c			1	ftot	fpos	fneg
c			...
c			lmax	ftot	fpos	fneg
c
c with
c
c	time	time in seconds
c	nsect	total number of sections
c	ivar	type of variable (see above)
c	lmax	number of layers for section
c	ftot	total flux through layer and section
c	fpos	positive flux through layer and section
c	fneg	negative flux through layer and section
c
c layers run from 1 to lmax
c layer 0 gives vertically integrated flux
c
c*****************************************************************

	program flxinf

c reads flx files

	implicit none

	integer nin
	integer idfile,nvers,nsect,kfluxm,nlmax
	integer nfxdi,nlvdi
	integer i,it,idtflx
	integer j,l,lmax,ierr,ivar,irec,ns
	integer, allocatable :: kflux(:)
	integer, allocatable :: nlayers(:)
	real, allocatable :: ptot(:,:)
	real, allocatable :: fluxes(:,:,:)

	integer iapini,ideffi

c---------------------------------------------------------------

        if(iapini(2,1,1,0).eq.0) then
                stop 'error stop: iapini'
        end if

        nin=ideffi('datdir','runnam','.flx','unform','old')
        if(nin.le.0) then
	  stop 'error stop: cannot open file'
	end if

	nvers = 5

        call infoflx      (nin,nvers
     +                          ,nsect,kfluxm,idtflx,nlmax
     +                          )

	nfxdi = kfluxm
	nlvdi = nlmax

	allocate(kflux(nfxdi))
	allocate(nlayers(nsect))
	allocate(ptot(4,nsect))
	allocate(fluxes(0:nlvdi,3,nsect))

        call rfflx        (nin,nvers
     +                          ,nfxdi,nfxdi,nlvdi
     +                          ,nsect,kfluxm,idtflx,nlmax
     +                          ,kflux
     +                          ,nlayers
     +                          )

	write(6,*) '    version    sections       nodes      layers'
	write(6,*) nvers,nsect,kfluxm,nlmax

c---------------------------------------------------------------
c loop on data
c---------------------------------------------------------------

	ns = min(10,nsect)
	irec = 0
	write(6,*) 'looping on data...'

	do while(.true.)

          call rdflx(nin,it,nlvdi,nsect,ivar,nlayers,fluxes,ierr)
	  if( ierr .ne. 0 ) goto 2

	  if( ivar .eq. 0 ) then
	   irec = irec + 1
	   if( mod(irec,100) .eq. 0 ) then
	    write(6,'(i10,10i7)') it,(nint(fluxes(0,1,i)),i=1,ns)
	   end if
	  end if

	  call fluxes_2d(it,nlvdi,nsect,ivar,ptot,fluxes)
	  call fluxes_3d(it,nlvdi,nsect,ivar,nlayers,fluxes)

    1	  continue
	end do

    2	continue

	close(nin)

	end

c****************************************************************

	subroutine fluxes_2d(it,nlvdim,nsect,ivar,ptot,fluxes)

c writes 2d fluxes to file (only for ivar=0)

	implicit none

	integer it			!time
	integer nlvdim			!vertical dimension
	integer nsect			!total number of sections
	integer ivar			!type of variable (0: water fluxes)
	real ptot(4,nsect)		!aux array
	real fluxes(0:nlvdim,3,nsect)	!fluxes

	integer iunit,i,lmax,l,j
	character*80, save :: format = ' '

	if( ivar .ne. 0 ) return

	if( format == ' ' ) then
	  !write(format,'(a,i4,a)') '(i10,',nsect,'f10.2)'
	  write(format,'(a,i4,a)') '(i10,',nsect,'g12.4)'
	  !write(6,*) 'using format: ',trim(format)
	end if

	do i=1,nsect
	  ptot(1,i) = fluxes(0,1,i)			!total
	  ptot(2,i) = fluxes(0,2,i)			!positive
	  ptot(3,i) = fluxes(0,3,i)			!negative
	  ptot(4,i) = fluxes(0,2,i) + fluxes(0,3,i)	!absolute
	end do

	write(66,format) it,(ptot(1,i),i=1,nsect)	!total
	write(67,format) it,(ptot(2,i),i=1,nsect)	!positive
	write(68,format) it,(ptot(3,i),i=1,nsect)	!negative
	write(69,format) it,(ptot(4,i),i=1,nsect)	!absolute

c next is box format for Ali

	write(61,*) it
	write(61,*) 0
	write(61,*) nsect
	do i=1,nsect
	  write(61,*) 0,0,(ptot(j,i),j=1,3)
	end do

	end

c****************************************************************

	subroutine fluxes_3d(it,nlvdim,nsect,ivar,nlayers,fluxes)

c writes 3d fluxes to file

	implicit none

	integer it			!time
	integer nlvdim			!vertical dimension
	integer nsect			!total number of sections
	integer ivar			!type of variable (0: water fluxes)
	integer nlayers(nsect)		!max layers for section
	real fluxes(0:nlvdim,3,nsect)	!fluxes

	integer iunit,i,lmax,l,j

	iunit = 160 + ivar

	write(iunit,*) it,nsect,ivar
	do i=1,nsect
	  lmax = nlayers(i)
	  write(iunit,*) i,lmax
	  do l=0,lmax
	    write(iunit,*) l,(fluxes(l,j,i),j=1,3)
	  end do
	end do

	end

c****************************************************************

