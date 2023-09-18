!
! $Id: subnosa.f,v 1.2 2006/09/28 08:56:59 georg Exp $
!
! auxiliar nos routines
!
! contents :
!
! wrnos2d(name,title,value)             write 2d nos file
!
! revision log :
!
! 05.01.2005    ggu     new routine wrnos2d() to write 2d nos file
! 28.09.2006    ggu     new routines wrnos2d_it, wrnos2d_index, extract_level
! 10.02.2015    ggu     some new routines for writing...
!
!***************************************************************************
!---------------------------------------------------------------------------
        module nos_util
!---------------------------------------------------------------------------
        contains
!---------------------------------------------------------------------------

        subroutine wrnos2d(name,title,value)

! write 2d nos file

        implicit none

        character*(*) name,title
        double precision value(1,*)

	integer it

	it = 0
	call wrnos2d_it(it,name,title,value)

	end

!***************************************************************************

        subroutine wrnos2d_it(it,name,title,value)

! write one 2d nos file

        use shympi
        use ionos

        implicit none

	integer it
        character*(*) name,title
        double precision value(1,*)

        integer nb,ivar

	ivar = 99

	call wrnos2d_open(nb,name,title)
	call wrnos2d_record(nb,it,ivar,value)

        if (shympi_is_master()) call delnos(nb)
        close(nb)

	end

!***************************************************************************

        subroutine wrnos2d_open(nb,name,title)

! write 2d nos file

	use depth
	use basin, only : nkn,nel,ngr,mbw,nkndi,neldi
        use mpi_io_admin
        use fil
        use ionos

        implicit none

	integer nb
        character*(*) name,title

	include 'param.h'

        character*80 pfile
        character*80 ptitle
        integer nvers,nvar,ivar,nlv
        integer ierr
        integer ilhkv(1)
        double precision hlv(1)

        integer laux

!-----------------------------------------------------------------
! initialize variables
!-----------------------------------------------------------------

        nvers = 3
        nvar = 1
        nlv = 1

        ptitle = title

!-----------------------------------------------------------------
! writing header
!-----------------------------------------------------------------

        pfile = name
        call filext(pfile,'.nos')
        write(6,*) 'writing file ',pfile(1:50)
        nb = ifileo(55,pfile,'unform','unknown')
        if( nb .le. 0 ) goto 98



        if(bmpi) then
          call rebuild_nos_header
          laux = shympi_max(nlv)
          if(shympi_is_master().and.shympi_partition_on_elements())then
            call wfnos(nb,nvers,nkndi,neldi,laux,nvar,ptitle,ierr)
            if(ierr.ne.0) goto 99
            call wsnos(nb,outIlhkv,hlv,outHev,ierr)
            if(ierr.ne.0) goto 99
          end if
        else
          call wfnos(nb,nvers,nkndi,neldi,nlv,nvar,ptitle,ierr)
          if(ierr.ne.0) goto 99
          call wsnos(nb,ilhkv,hlv,hev,ierr)
          if(ierr.ne.0) goto 99
        end if

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        return
   98   continue
        write(6,*) pfile
        write(6,*) nb
        stop 'error stop wrnos2d_open: opening file'
   99   continue
        write(6,*) pfile
        write(6,*) ierr
        stop 'error stop wrnos2d_open: writing header'
        end

!***************************************************************************

	subroutine wrnos2d_record(nb,it,ivar,value)

! write 2d nos file

        use shympi
        use mpi_io_admin
        use ionos

        implicit none

	integer nb,it,ivar
        double precision value(1,*)

	integer nlv,ierr
        integer ilhkv(1)

        nlv = 1
        ierr = 0
        if(shympi_partition_on_elements())then
          call wrnos(nb,it,ivar,nlv,outIlhkv,value,ierr)
        else 
          call wrnos(nb,it,ivar,nlv,ilhkv,value,ierr)
        end if
        if( ierr .ne. 0 ) goto 99

	return
   99   continue
        write(6,*) nb
        write(6,*) ierr
        stop 'error stop wrnos2d_record: writing record'
        end

!***************************************************************************

        subroutine wrnos2d_index(it,index,name,title,value)

	implicit none

	integer it
	integer index
        character*(*) name,title
        double precision value(1)

	integer i
	integer chars,imax,imin
	character*30 number,format
	character*80 file

	chars = 5
	format = '(i5)'

	imax = 10**chars
	imin = -10**(chars-1)
	if( index .ge. imax .or. index .le. imin ) goto 99

	write(number(1:chars),format) index
	do i=1,chars
	  if( number(i:i) .eq. ' ' ) number(i:i) = '0'
	end do

	file = name//number(1:chars)
        call wrnos2d_it(it,file,title,value)

	return
   99	write(6,*) 'index is too big for max chars: ',index,chars
	stop 'error stop wrnos2d_index: index'
	end

!***************************************************************************

	subroutine extract_level(nlvddi,nkn,level,v3,v2)

	implicit none

	integer nlvddi
	integer nkn
	integer level
	double precision v3(nlvddi,1)
	double precision v2(1)

	integer k

	do k=1,nkn
	  v2(k) = v3(level,k)
	end do

	end

!***************************************************************************

!---------------------------------------------------------------------------
        end module nos_util
!---------------------------------------------------------------------------
