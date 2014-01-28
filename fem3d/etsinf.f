c
c $Id: etsinf.f,v 1.8 2008-11-20 10:51:34 georg Exp $
c
c revision log :
c
c 24.01.2014    ggu     copied from nosinf.f
c
c**************************************************************

	program etsinf

c reads ets file

	implicit none

	include 'param.h'

        integer ndim
        parameter (ndim=100)

	real cv(nkndim)
	real cv3(nlvdim,nexdim)
	integer ivars(ndim)

	integer ilets(nexdim)
	real hlv(nlvdim)
	real hets(nexdim)
	integer nodes(nexdim)
	real xg(nexdim)
	real yg(nexdim)
	character*80 desc(nexdim)

	logical bwrite
	integer nread,nin
	integer nvers
	integer nkn,nel,nlv,nvar
	integer ierr
	integer it,ivar
	integer l,k,i,lmax
	character*80 title
	real rnull
	real cmin,cmax

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	bwrite = .false.

c--------------------------------------------------------------
c open simulation
c--------------------------------------------------------------

	call shyfem_copyright('etsinf - Info on ETS files')

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call open_ets_type('.ets','old',nin)

        call read_ets_header(nin,nkndim,nlvdim,ilets,hlv,hets
     +                                  ,nodes,xg,yg,desc)
	call ets_get_params(nin,nkn,nlv,nvar)

        if( nvar .gt. ndim ) goto 95
        call ets_get_vars(nin,nvar,ivars)

        write(6,*) 'Available variables: ',nvar
        write(6,*) (ivars(i),i=1,nvar)

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

	  call ets_read_record(nin,it,ivar,nlvdim,ilets,cv3,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) goto 100

	  nread=nread+1
	  write(6,*) 'time : ',it,'   ivar : ',ivar

	  if( bwrite ) then
	    do i=1,nkn
	      lmax = ilets(i)
	      if( ivar .eq. 1 ) lmax = 1
	      write(6,*) i,lmax,(cv3(l,i),l=1,lmax)
	    end do
	  end if

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   95	continue
	write(6,*) 'nvar = ',nvar,'   ndim = ',ndim
	stop 'error stop: ndim'
	end

c***************************************************************

