c
c $Id: subdim.f,v 1.11 2009-04-21 10:38:59 georg Exp $
c
c dimension handling routines
c
c contents :
c
c subroutine sp131a(nknddi,nelddi)		check dimension for nkn,nel
c subroutine sp131b(nrbdim,nbcdim)		check dimension for nrb,nbc
c subroutine sp131k(matdim)			check dimension for matrix
c subroutine sp131g(matdim,mbwdim,ngl,mbw,isym) check matrix dimension
c subroutine sp131m(mbwdim)			check dimension for band width
c subroutine cstdim(nknddi,...)			writes dimension to common block
c
c revision log :
c
c 31.05.1997	ggu	unnecessary routines deleted
c 27.06.1997	ggu     dim routines into own file
c 29.06.1997	ggu     cstdim to this file
c 28.05.1998	ggu     call to sp131g changed (no index,...)
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 24.06.2003    ggu     new routines to list entries in dimensions
c 24.06.2003    ggu     bugfix NFILL
c 10.08.2003    ggu     new dimension nlidim added
c 18.10.2005    ggu     error messages slightly changed
c 06.04.2009    ggu     nlidim -> nlkdim
c 21.04.2009    ggu     new routine sp131m() for band width
c
c***********************************************************

	subroutine sp131a(nknddi,nelddi)

c check dimension for nkn,nel
c
c nknddi	dimension for nkn
c nelddi	dimension for nel

	implicit none

	integer nknddi,nelddi

	include 'nbasin.h'

	if(nknddi.lt.nkn) goto 99
	if(nelddi.lt.nel) goto 99

	if(nkn.le.0) goto 98
	if(nel.le.0) goto 98

	return
   98	continue
	write(6,*) 'error in constants'
	write(6,*) 'nkn,nel :',nkn,nel
	stop 'error stop sp131a: constants'
   99	continue
	write(6,*) 'dimension error'
	write(6,*) 'nkn,nel :',nkn,nel
	write(6,*) 'nknddi,nelddi :',nknddi,nelddi
	stop 'error stop sp131a: dimensions'
	end

c*************************************************

	subroutine sp131b(nrbdim,nbcdim)

c check dimension for nrb,nbc
c
c nrbdim	dimension for nrb
c nbcdim	dimension for nbc

	implicit none

	integer nrbdim,nbcdim

	integer nrb,nbc
	!integer nkbnd,nbnds

	!nrb = nkbnd()
	!nbc = nbnds()

	!if(nrbdim.lt.nrb) goto 99
	!if(nbcdim.lt.nbc) goto 99

	!if(nrb.lt.0) goto 98
	!if(nbc.lt.0) goto 98

	return
   98	continue
	write(6,*) 'error in constants'
	write(6,*) 'nrb,nbc :',nrb,nbc
	stop 'error stop sp131b: constants'
   99	continue
	write(6,*) 'dimension error'
	write(6,*) 'nrb,nbc :',nrb,nbc
	write(6,*) 'nrbdim,nbcdim :',nrbdim,nbcdim
	stop 'error stop sp131b: dimensions'
	end

c*************************************************

	subroutine sp131k(matdim)

c check dimension for global matrix (probably not needed anymore)
c
c matdim	dimension for matrix

	implicit none

	integer matdim

	include 'nbasin.h'

	if(matdim.lt.nkn*ngr) goto 99

	if(nkn.le.0) goto 98
	if(ngr.le.0) goto 98

	return
   98	continue
	write(6,*) 'error in constant'
	write(6,*) 'nkn,ngr :',nkn,ngr
	stop 'error stop sp131k: constants'
   99	continue
	write(6,*) 'dimension error'
	write(6,*) 'nkn,ngr :',nkn,ngr
	write(6,*) 'nkn*ngr :',nkn*ngr
	write(6,*) 'matdim :',matdim
        write(6,*) 'Dimension matdim must be at least nkn*ngr'
	write(6,*) 'Be sure that ngrdim is greater or equal to ngr'
	stop 'error stop sp131k: dimensions'
	end

c*************************************************

	subroutine sp131g(matdim,mbwdim,ngl,mbw,isym)

c check dimension for matrix and set up index

	implicit none

	integer matdim,mbwdim,ngl,mbw,isym
	integer matlen

	if(mbwdim.lt.mbw) goto 99
	if(mbw.le.0) goto 98
	if(ngl.lt.mbw) goto 97
	if(isym.lt.0.or.isym.gt.2) goto 96

	matlen=ngl*(1+(2-isym)*mbw)

	if(matdim.lt.matlen) goto 95

	return
   95	continue
	write(6,*) 'dimension error'
	write(6,*) 'matdim,matlen :',matdim,matlen
	stop 'error stop sp131g: matdim'
   96	continue
	write(6,*) 'error in constants'
	write(6,*) 'isym :',isym
	stop 'error stop sp131g: isym'
   97	continue
	write(6,*) 'error in constants'
	write(6,*) 'ngl,mbw :',ngl,mbw
	stop 'error stop sp131g: ngl'
   98	continue
	write(6,*) 'error in constants'
	write(6,*) 'mbw :',mbw
	stop 'error stop sp131g: mbw'
   99	continue
	write(6,*) 'dimension error'
	write(6,*) 'mbwdim,mbw :',mbwdim,mbw
	stop 'error stop sp131g: mbwdim'
	end

c*************************************************

	subroutine sp131m(mbwdim)

c check dimension for band width

	implicit none

	integer mbwdim

	include 'nbasin.h'

	if(mbwdim.lt.mbw) goto 99
	if(mbw.le.0) goto 98

	return
   98	continue
	write(6,*) 'error in constants'
	write(6,*) 'mbw :',mbw
	stop 'error stop sp131m: mbw'
   99	continue
	write(6,*) 'dimension error'
	write(6,*) 'mbwdim,mbw :',mbwdim,mbw
	stop 'error stop sp131m: mbwdim'
	end

c******************************************

	subroutine cstdim(nknddi,nelddi,nrbdim,nbcdim,mbwdim,ngrdim
     +                  ,nardim,nexdim,nfxdim,nlkdim)

c stores dimensions for later use

	implicit none

	integer nknddi,nelddi,nrbdim,nbcdim,mbwdim,ngrdim
	integer nardim,nexdim,nfxdim,nlkdim

	call setdim('nkndim',nknddi)
	call setdim('neldim',nelddi)
	call setdim('nrbdim',nrbdim)
	call setdim('nbcdim',nbcdim)
	call setdim('mbwdim',mbwdim)
	call setdim('ngrdim',ngrdim)
	call setdim('nardim',nardim)
	call setdim('nexdim',nexdim)
	call setdim('nfxdim',nfxdim)
	call setdim('nlkdim',nlkdim)

	end

c*************************************************
c*************************************************
c*************************************************

	subroutine setdim(name,value)

c sets dimension

	implicit none

	character*6 name
	integer value

	call handim(1,name,value)

	end

c*************************************************

	subroutine getdim(name,value)

c gets dimension

	implicit none

	character*6 name
	integer value

	call handim(2,name,value)

	end

c*************************************************

	subroutine listdim

c lists dimension array

	implicit none

	character*6 name
	integer value

	call handim(3,name,value)

	end

c*************************************************

	subroutine handim(mode,name,value)

c sets dimension

	implicit none

	integer mode
	character*6 name
	integer value

	integer ndim
	parameter( ndim = 20 )

	integer n
	integer getentry

	integer nfill
	character*6 dims(ndim)
	integer vals(ndim)
	save nfill,dims,vals
	data nfill /0/

c----------------------------------------------------------
c do we have to list the entries?
c----------------------------------------------------------

	if( mode .eq. 3 ) then

	   call listentry(nfill,dims,vals)
	   return

	end if

c----------------------------------------------------------
c find out if alreay in list
c----------------------------------------------------------

	n = getentry(name,dims,nfill)

	if( mode .eq. 1 ) then

c	   ----------------------------------------------------------
c	   set new value
c	   ----------------------------------------------------------

	   if( n .gt. 0 ) then
	      write(6,*) 'dimension name already in use : ',name
	      stop 'error stop handim'
	   end if

	   nfill = nfill + 1

	   if( nfill .gt. ndim ) then
	      write(6,*) 'dimension error in ndim : ',ndim
	      stop 'error stop handim'
	   end if

	   dims(nfill) = name
	   vals(nfill) = value

	else if( mode .eq. 2 ) then

c	   ----------------------------------------------------------
c	   get value
c	   ----------------------------------------------------------

	   if( n .le. 0 ) then
	      write(6,*) 'dimension name not found : ',name
	      stop 'error stop handim'
	   end if
	   !value = vals(nfill)	!NFILL terrible bug -> fixed 24/6/2003
	   value = vals(n)

	else

c	   ----------------------------------------------------------
c	   error in mode
c	   ----------------------------------------------------------

	   write(6,*) 'mode not recognized : ',mode
	   stop 'error stop handim'

	end if

	end

c*************************************************

	function getentry(name,dims,nfill)

c gets entry from list

	implicit none

	integer getentry
	character*6 name
	character*6 dims(1)
	integer nfill

	integer i

	do i=1,nfill
	  if( name .eq. dims(i) ) then
	    getentry = i
	    return
	  end if
	end do

	getentry = 0

	end

c*************************************************

	subroutine listentry(nfill,dims,vals)

c gets entry from list

	implicit none

	integer nfill
	character*6 dims(1)
	integer vals(1)

	character*6 name
	integer i,val

	write(6,*) '---------------------------'
	write(6,*) 'list dimensions: nfill = ',nfill
	do i=1,nfill
	  name = dims(i)
	  val  = vals(i)
	  write(6,'(a6,5x,i10)') name,val
	end do
	write(6,*) '---------------------------'

	end

c*************************************************


