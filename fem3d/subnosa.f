c
c $Id: subnosa.f,v 1.2 2006/09/28 08:56:59 georg Exp $
c
c auxiliar nos routines
c
c contents :
c
c wrnos2d(name,title,value)             write 2d nos file
c
c revision log :
c
c 05.01.2005    ggu     new routine wrnos2d() to write 2d nos file
c 28.09.2006    ggu     new routines wrnos2d_it, wrnos2d_index, extract_level
c 10.02.2015    ggu     some new routines for writing...
c
c***************************************************************************

        subroutine wrnos2d(name,title,value)

c write 2d nos file

        implicit none

        character*(*) name,title
        real value(1)

	integer it

	it = 0
	call wrnos2d_it(it,name,title,value)

	end

c***************************************************************************

        subroutine wrnos2d_it(it,name,title,value)

c write one 2d nos file

        implicit none

	integer it
        character*(*) name,title
        real value(1)

        integer nb,ivar

	ivar = 99

	call wrnos2d_open(nb,name,title)
	call wrnos2d_record(nb,it,ivar,value)

	call delnos(nb)
        close(nb)

	end

c***************************************************************************

        subroutine wrnos2d_open(nb,name,title)

c write 2d nos file

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	integer nb
        character*(*) name,title

	include 'param.h'

        character*80 pfile
        character*80 ptitle
        integer nvers,nvar,ivar,nlv
        integer ierr
        integer ilhkv(1)
        real hlv(1)

        integer ifileo

c-----------------------------------------------------------------
c initialize variables
c-----------------------------------------------------------------

        nvers = 3
        nvar = 1
        nlv = 1

        ptitle = title

c-----------------------------------------------------------------
c writing header
c-----------------------------------------------------------------

        pfile = name
        call filext(pfile,'.nos')
        write(6,*) 'writing file ',pfile(1:50)
        nb = ifileo(55,pfile,'unform','unknown')
        if( nb .le. 0 ) goto 98

        call wfnos(nb,nvers,nkn,nel,nlv,nvar,ptitle,ierr)
        if(ierr.ne.0) goto 99
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 99

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

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

c***************************************************************************

	subroutine wrnos2d_record(nb,it,ivar,value)

c write 2d nos file

        implicit none

	integer nb,it,ivar
        real value(1)

	integer nlv,ierr
        integer ilhkv(1)

        nlv = 1

        call wrnos(nb,it,ivar,nlv,ilhkv,value,ierr)
        if( ierr .ne. 0 ) goto 99

	return
   99   continue
        write(6,*) nb
        write(6,*) ierr
        stop 'error stop wrnos2d_record: writing record'
        end

c***************************************************************************

        subroutine wrnos2d_index(it,index,name,title,value)

	implicit none

	integer it
	integer index
        character*(*) name,title
        real value(1)

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

c***************************************************************************

	subroutine extract_level(nlvddi,nkn,level,v3,v2)

	implicit none

	integer nlvddi
	integer nkn
	integer level
	real v3(nlvddi,1)
	real v2(1)

	integer k

	do k=1,nkn
	  v2(k) = v3(level,k)
	end do

	end

c***************************************************************************

