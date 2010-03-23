
c $Id: subbnds.f,v 1.18 2010-02-26 17:35:06 georg Exp $
c
c handles open boundary conditions for scalar variables
c
c contents :
c
c subroutine bnds_init(text,file,nintp,nvar,ndim,array,aconst)
c			initializes boundary condition
c subroutine bnds_set(text,t,ndim,array,aaux)
c			sets boundary condition
c subroutine bnds_trans(text,ndim,array,aaux,ivar,nlvdim,r3v)
c			transfers boundary condition to matrix
c subroutine bnds_set_def(text,ndim,array)
c			sets default value for boundaries
c subroutine bnds_print(text,ndim,array)
c			prints boundary condition
c
c revision log :
c
c 10.08.2003	ggu	exit if no file can be opened
c 03.03.2005	ggu	new call to exfini
c 01.02.2006	ggu	bugfix -> nextra was 3 -> not used anymore
c 03.02.2006	ggu	bugfix -> calling exfpres
c 17.02.2006	ggu	included text for debug, new call to get_bflux()
c 23.03.2006    ggu     changed time step to real
c 05.10.2007    ggu     definition of array(ndim,1) changed to array(ndim,0:1)
c 17.03.2008    ggu     routines re-arranged, new bnds_trans, bnds_set_def
c 17.04.2008    ggu     deleted bnds_set_global
c 23.04.2008    ggu     in bnds_set_def() eliminated aaux
c
c******************************************************************

	subroutine bnds_init(text,file,nintp,nvar,ndim,array,aconst)

c initializes boundary condition

	implicit none

	character*(*) text	!text for debug
	character*80 file(1)	!file names
	integer nintp		!degree of interpolation - same for all bounds
	integer nvar		!number of variables
	integer ndim		!first dimension of array
	real array(ndim,0:1)	!array with all information
	real aconst(nvar)	!if no file is given constant values from here

	character*80 name
	integer nbc,ibc
	integer nbdim,nk,nsize
	integer iunit,n,i

	integer nbnds,nkbnds,ifileo
	logical bdebug

	bdebug = .false.
	iunit = 0
	nbc = nbnds()

	if( bdebug ) then
	  write(88,*) '-------------------------'
	  write(88,*) 'Initialization for scalar:'
	  write(88,*) text
	  write(88,*) '-------------------------'
	end if

	do i=1,ndim
	  array(i,0) = 0.
	end do

	do ibc=1,nbc

	  do i=1,ndim
	    array(i,ibc) = 0.
	  end do

	  name = file(ibc)
	  if( name .ne. ' ' ) then
	    iunit = ifileo(iunit,name,'form','old')
	    if( iunit .le. 0 ) goto 99
	  else
	    iunit = -1
	  end if

	  call get_bnd_ipar(ibc,'nbdim',nbdim)	!0 if one value/boundary
	  nk = nkbnds(ibc)
	  nsize = nbdim * nk

	  if( nbdim .gt. 0 ) then
	    write(6,*) 'variable boundary found: ibc ',ibc,nsize,nbdim
	  end if

	  call exfini(iunit,nintp,nvar,nsize,ndim,array(1,ibc))

	  call exfsetdef(array(1,ibc),aconst)

	  if( bdebug ) then
	    write(88,*) '-------------------------'
            write(88,*) 'bnds init: ',ibc,iunit,nvar,nbdim,nsize
            write(88,*) '      ',nintp,ndim,nk
            write(88,*) name
            write(88,*) (aconst(i),i=1,nvar)
	    write(88,*) '-------------------------'
	  end if
	end do

	return
   99	continue
	write(6,*) 'Error opening file: ',name
	stop 'error stop bnds_init: cannot open file'
	end

c******************************************************************

	subroutine bnds_set(text,t,ndim,array,aaux)

c sets boundary condition

	implicit none

	character*(*) text	!text for debug
	real t			!time for interpolation
	integer ndim		!first dimension of array
	real array(ndim,0:1)	!array with all information
	real aaux(ndim)		!aux array

	integer nbc,ibc
	integer iqual
	integer iunit,nvar
	integer i
	character*30 ltext

	integer nbnds

	nbc = nbnds()
	ltext = text

	do ibc=1,nbc

	  call get_bnd_ipar(ibc,'iqual',iqual)

	  call exfunit(array(1,ibc),iunit)	!gets unit number

	  if( iunit .le. 0 ) then
c	    nothing, everything in place
	  else if( iqual .le. 0 ) then
	    call exfintp(array(1,ibc),t,aaux)
	    call exfset(array(1,ibc),t,aaux)
	  else
	    call exfget(array(1,iqual),t,aaux)
	    call exfset(array(1,ibc),t,aaux)
	  end if

	!if( ibc .eq. 3 ) then
	!  call exfinfo(77,array(1,ibc))
	!end if

	end do

	!call bnds_print('debug output for '//ltext,ndim,array)

	return
	end

c******************************************************************

	subroutine bnds_trans(text,ndim,array,aaux,ivar,nlvdim,r3v)

c transfers boundary condition to matrix

	implicit none

	character*(*) text	!text for debug
	integer ndim		!first dimension of array
	real array(ndim,0:1)	!array with all information
	real aaux(ndim)		!aux array
	integer ivar		!variable to use (can be 0 -> 1)
	integer nlvdim		!vertical dimension of levels
	real r3v(nlvdim,1)	!matrix to which BC values are transfered

	integer nbc,ibc
	integer nvar,nsize,ndata,nbdim
	integer nk,iv,kn
	integer i,ip
	real t
	character*30 ltext

	integer nbnds,nkbnds,kbnds

	call init_scal_bc(r3v)

	nbc = nbnds()
	ltext = text
	iv = max(ivar,1)

	do ibc=1,nbc

	  nk = nkbnds(ibc)   !total number of nodes of this boundary
	  call exfsize(array(1,ibc),nvar,nsize,ndata)
	  nsize = max(nsize,1)
	  nbdim = nsize/nk
	  if( nsize .gt. 1 .and. mod(nsize,nk) .ne. 0 ) then
	    write(6,*) 'nbdim,nsize,nk: ',nbdim,nsize,nk
	    stop 'error stop bnds_trans: nsize not multiple of nk'
	  end if

	  !write(61,*) 'bc: ',ibc,nsize,nbdim,nk

	  call exfgetvar(array(1,ibc),iv,t,aaux)

	  !write(61,*) 'befor 1 ',text,(aaux(i),i=1,nsize)

	  if( nsize .le. 1 ) then	!distribute over nodes of boundary
	    do i=1,nk
	      aaux(i) = aaux(1)
	    end do
	  end if

	  !write(61,*) 'after 1 ',text,ibc,nsize,nbdim,nk,(aaux(i),i=1,nsize)

	  do i=1,nk
	    kn = kbnds(ibc,i)
	    ip = 1 + (i-1) * nbdim             !also works for nbdim=0
	    call dist_3d(nlvdim,r3v,kn,nbdim,aaux(ip))
	  end do

	end do

	end

c******************************************************************

	subroutine bnds_set_def(text,ndim,array)

c sets default value for boundaries - works only for nvar = 1

	implicit none

	character*(*) text	!text for debug
	integer ndim		!first dimension of array
	real array(ndim,0:1)	!array with all information

	integer nbc,ibc
	integer iunit,nvar,nsize,ndata
	integer i
	real scaldf

	integer nbnds

	nbc = nbnds()

	do ibc=1,nbc

	  call exfsize(array(1,ibc),nvar,nsize,ndata)
	  if( nvar .gt. 1 ) goto 99
	  call get_bnd_par(ibc,text,scaldf)
	  call exfsetdef(array(1,ibc),scaldf)	!works because nvar = 1

	end do

	return
   99	continue
	write(6,*) 'nvar = ',nvar
	stop 'error stop bnds_set_def: nvar'
	end

c******************************************************************

	subroutine bnds_print(text,ndim,array)

c prints boundary condition

	implicit none

	character*(*) text	!text for debug
	integer ndim		!first dimension of array
	real array(ndim,0:1)	!array with all information

	integer nbc,ibc

	integer nbnds

	nbc = nbnds()

	write(89,*) 'Boundary Info : ',text

	do ibc=1,nbc
	  write(89,*) '  ...Boundary : ',ibc
	  call exfinfo(89,array(1,ibc))
	end do

	end

c******************************************************************

