c $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
c
c reading and interpolation of external files
c
c contents :
c
c revision log :
c
c 29.10.2012    ggu     created from scratch
c 17.06.2013    ggu     do not pass finction into subroutine
c
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open(name,nkn,iunit)
	integer iunit(3)
	character*(*) name
	!call fem_file_open_1(name,nkn,iunit)
	call ts_file_open_0(name,nkn,iunit)		!old call
	end

	subroutine ts_next_record(it,iunit,nkn,nlv,value)
	include 'param.h'
	integer iunit(3)
	real value(nlvdim,1)
	!call ts_next_record_1(it,iunit,nkn,nlv,value)
	call ts_next_record_0(it,iunit,nkn,nlv,value)	!old call
	end

	subroutine ts_file_close(info)
	integer info(3)
	!call fem_file_close_1(info)
	call ts_file_close_0(info)			!old call
	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine fem_file_open_1(name,np,iunit)

c opens T/S file

	implicit none

	include 'param.h'

	character*(*) name		!name of file
	integer np			!number of points expected
	integer iunit(3)		!unit number (return)

	logical bformat,bdebug
	integer iu

	bdebug = .true.
	bdebug = .false.

	iu = 0
	call fem_file_read_open(name,np,iu,bformat)
	if ( iu .eq. 0 ) goto 99

	if( .not. bformat ) iu = -iu

	if( bdebug ) then
	  write(6,*) 'debug in ts_file_open: ',name
	  write(6,*) 'iunit,bformat,np: ',iu,bformat,np
	end if

	iunit(1) = iu

	return
   99	continue
	write(6,*) 'Cannot open file: ',name
	stop 'error stop fem_file_open: error file open'
	end

c*******************************************************************	

	subroutine ts_next_record_1(it,iunit,nkn,nlv,value)

	implicit none

	include 'param.h'

	integer nldim
	parameter (nldim=100)

	integer it
	integer iunit(3)
	integer nkn
	integer nlv
	real value(nlvdim,1)

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !integer nlvdi,nlv
        !common /level/ nlvdi,nlv

	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hlv(1)
	common /hlv/hlv
	real znv(1)
	common /znv/znv
	real hkv(1)
	common /hkv/hkv

	logical bdebug,bcons,bformat
	integer iu,lmax,ierr,k
        real vmin,vmax
	character*80 string

	integer il_data(nkndim)
	real hd_data(nkndim)
	real zz_data(nkndim)
	real data(nldim,nkndim)
        real hl_data(0:nldim+1)

	logical intp_to_fem_0
	external intp_to_fem_0

c--------------------------------------------------------------
c initialize
c--------------------------------------------------------------

	bdebug = .true.
	bdebug = .false.

	bcons = .false.		!conserve total quantity?

	iu = iunit(1)
	bformat = iu .gt. 0
	iu = abs(iu)

	if( iu .eq. 0 ) return

	if( bdebug ) then
	  write(6,*) 'debug in ts_next_record:'
	  write(6,*) 'iunit,bformat: ',it,iu,bformat
	end if

	do k=1,nkn
	  zz_data(k) = 0.
	end do

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	write(6,*)'reading T/S values', it

        call fem_file_read_3d(bformat,iu,it
     +                          ,nkn,lmax,nldim,hl_data(1)
     +                          ,il_data,string,hd_data,data
     +				,ierr)

c--------------------------------------------------------------
c interpolate between different vertical structures 
c--------------------------------------------------------------

c	call intp_to_fem(nkn,nldim
c     +		,lmax,il_data,hl_data,hd_data,zz_data,data
c     +		,nlv,ilhkv,hlv,hkv,znv,value)

c--------------------------------------------------------------
c some statistics
c--------------------------------------------------------------

        call conmima(nlvdim,value,vmin,vmax)
        write(6,*) 'min/max: ',vmin,vmax

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine fem_file_close_1(iunit)

c closes T/S file

	implicit none

	integer iunit(3)

	close(abs(iunit(1)))

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	function intp_to_fem_0(nldim
     +		,lm_data,hl_data,hd_data,zz_data,val_data
     +		,lm_fem,il_fem,hl_fem,hd_fem,zz_fem,val_fem)

	implicit none

	include 'param.h'

	logical intp_to_fem_0
	integer lm_data			!vertical dimension of data
	real hl_data(lm_data)
	!integer lm_data			!max level of data
	real hd_data
	real zz_data
	real val_data(lm_data)
	integer nldim			!max level of fem
	integer lm_fem			!max level of fem
	real hd_fem
	real val_fem(nldim)

	integer il_fem(1)	!to be changed...
	real hl_fem(nlvdim)
	real zz_fem(1)

	real hlv(1)
	common /hlv/hlv

	integer ndim
	parameter (ndim=1000)
	
	real vaux_data(ndim+1)
	real haux_data(0:ndim+1)
	real vaux_fem(nlvdim)
	real haux_fem(0:nlvdim)

	logical bsigma,bdebug,bcons,bsigma_data
	integer l,k
	integer nsigma_data,nsigma
	real hsigma_data,hsigma
	integer lmax_fem,lmax_data
	real z_fem,h_fem
	real z_data,h_data

	if( nldim .gt. ndim ) stop 'errro stop intp_scalar_to_fem: ndim'

c--------------------------------------------------------------
c initialize parameters and variables
c--------------------------------------------------------------

	bcons = .false.

        call compute_sigma_info(lm_data,hl_data,nsigma_data,hsigma_data)
        bsigma_data = nsigma_data .gt. 0
	do l=1,lm_data
	  haux_data(l) = hl_data(l)
	end do
	haux_data(0) = 0.

        call get_sigma(nsigma,hsigma)		!from basin
        bsigma = nsigma .gt. 0
	do l=1,lm_fem
	  haux_fem(l) = hl_fem(l)
	end do
	haux_fem(0) = 0.

c--------------------------------------------------------------
c loop on nodes
c--------------------------------------------------------------

        !do k = 1,np

c	  -----------------------------------------------------
c	  set data depth structure
c	  -----------------------------------------------------

    !      lmax_data = il_data(k)
	  z_data = zz_data(k)
	  h_data = hd_data(k)

	  if( bsigma_data ) then
            call set_hybrid_depth(lmax_data,z_data,h_data
     +			,hl_data,nsigma_data,hsigma_data,haux_data(1))
	  end if
          haux_data(0) = -z_data

c	  -----------------------------------------------------
c	  set fem depth structure
c	  -----------------------------------------------------

          lmax_fem = il_fem(k)
	  z_fem = zz_fem(k)
	  h_fem = hd_fem(k)

	  if( bsigma ) then
            call set_hybrid_depth(lmax_fem,z_fem,h_fem
     +			,hl_fem,nsigma,hsigma,haux_fem(1))
	  end if
          haux_fem(0) = -z_fem

c	  -----------------------------------------------------
c	  interpolate on vertical levels
c	  -----------------------------------------------------

          do l=1,lmax_data
    !        vaux_data(l) = val_data(l,k)
          end do

          call intp_vert(bcons,lmax_data,haux_data,vaux_data
     +				,lmax_fem,haux_fem,vaux_fem)

          do l = 1,lmax_fem
    !        val_fem(l,k) = vaux_fem(l)
          end do

        !enddo

	intp_to_fem_0 = .true.

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	

	subroutine intp_to_fem(np,nldim
     +		,lm_data,il_data,hl_data,hd_data,zz_data,val_data
     +		,lm_fem,il_fem,hl_fem,hd_fem,zz_fem,val_fem)

	implicit none

	include 'param.h'

	integer np			!total number of horizontal points
	integer nldim			!vertical dimension of data
	integer lm_data			!max level of data
	integer il_data(1)
	real hl_data(nldim)
	real hd_data(1)
	real zz_data(1)
	real val_data(nldim,1)
	integer lm_fem			!max level of fem
	integer il_fem(1)
	real hl_fem(nlvdim)
	real hd_fem(1)
	real zz_fem(1)
	real val_fem(nlvdim,1)

	integer ndim
	parameter (ndim=1000)
	
	real vaux_data(ndim+1)
	real haux_data(0:ndim+1)
	real vaux_fem(nlvdim)
	real haux_fem(0:nlvdim)

	logical bsigma,bdebug,bcons,bsigma_data
	integer l,k
	integer nsigma_data,nsigma
	real hsigma_data,hsigma
	integer lmax_fem,lmax_data
	real z_fem,h_fem
	real z_data,h_data

	if( nldim .gt. ndim ) stop 'errro stop intp_scalar_to_fem: ndim'

c--------------------------------------------------------------
c initialize parameters and variables
c--------------------------------------------------------------

	bcons = .false.

        call compute_sigma_info(lm_data,hl_data,nsigma_data,hsigma_data)
        bsigma_data = nsigma_data .gt. 0
	do l=1,lm_data
	  haux_data(l) = hl_data(l)
	end do
	haux_data(0) = 0.

        call get_sigma(nsigma,hsigma)		!from basin
        bsigma = nsigma .gt. 0
	do l=1,lm_fem
	  haux_fem(l) = hl_fem(l)
	end do
	haux_fem(0) = 0.

c--------------------------------------------------------------
c loop on nodes
c--------------------------------------------------------------

        do k = 1,np

c	  -----------------------------------------------------
c	  set data depth structure
c	  -----------------------------------------------------

          lmax_data = il_data(k)
	  z_data = zz_data(k)
	  h_data = hd_data(k)

	  if( bsigma_data ) then
            call set_hybrid_depth(lmax_data,z_data,h_data
     +			,hl_data,nsigma_data,hsigma_data,haux_data(1))
	  end if
          haux_data(0) = -z_data

c	  -----------------------------------------------------
c	  set fem depth structure
c	  -----------------------------------------------------

          lmax_fem = il_fem(k)
	  z_fem = zz_fem(k)
	  h_fem = hd_fem(k)

	  if( bsigma ) then
            call set_hybrid_depth(lmax_fem,z_fem,h_fem
     +			,hl_fem,nsigma,hsigma,haux_fem(1))
	  end if
          haux_fem(0) = -z_fem

c	  -----------------------------------------------------
c	  interpolate on vertical levels
c	  -----------------------------------------------------

          do l=1,lmax_data
            vaux_data(l) = val_data(l,k)
          end do

          call intp_vert(bcons,lmax_data,haux_data,vaux_data
     +				,lmax_fem,haux_fem,vaux_fem)

          do l = 1,lmax_fem
            val_fem(l,k) = vaux_fem(l)
          end do

        enddo

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c****** old routines ***********************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open_0(name,nkn,info)

c opens T/S file

	implicit none

	include 'param.h'

	character*(*) name
	integer nkn
	integer info(3)

	real hlv(1)
	common /hlv/hlv

	logical bsigma,bdebug
	logical bformat,bhashl
	integer ifileo,iunit,ios
	integer iformat,ihashl,inzeta
	integer it,nknaux,lmax,nvar
	integer l
	integer nsigma
	real hsigma
	real hl(nlvdim)

	bdebug = .false.
	bdebug = .true.

c-------------------------------------------------------------
c indicator for file format and content
c
c -1	do not know - check (default)
c  0	no (not formatted, does not have hlv array)
c  1	yes (formatted, has hlv array)
c
c if you are sure about your file format and content you can
c set these variables below to the desired value
c
c we still check file for formatted or unformatted
c we always pretend that the file has hlv information
c-------------------------------------------------------------

	iformat = -1		!indicator if file is formatted
	ihashl  = +1		!indicator if file has hlv array

c-------------------------------------------------------------
c check if formatted or unformatted
c-------------------------------------------------------------

	nknaux = -1

	if( iformat .eq. -1 ) then		!check if formatted
	  iunit = ifileo(0,name,'unform','old')
	  if( iunit .le. 0 ) goto 99
	  read(iunit,iostat=ios) it,nknaux,lmax,nvar
	  close(iunit)

	  if( ios .lt. 0 ) goto 98		!EOF

	  if( nkn .eq. nknaux ) then		!is is probably unformatted
	    iformat = ios
	  else					!let us try formatted
	    iunit = ifileo(0,name,'form','old')
	    if( iunit .le. 0 ) goto 99
	    read(iunit,*,iostat=ios) it,nknaux,lmax,nvar
	    close(iunit)
	    if( ios .ne. 0 ) then
	      stop 'error stop ts_file_open_0: cannot read file'
	    end if
	    iformat = 1
	  end if
	end if

	bformat = iformat .gt. 0

c-------------------------------------------------------------
c check if level structure information is available
c-------------------------------------------------------------

	ios = 0
	if( bformat ) then
	  write(6,*) 'ts_file_open: opening formatted file ',name
	  iunit = ifileo(0,name,'form','old')
	  read(iunit,*) it,nknaux,lmax,nvar
	  if( lmax .gt. nlvdim ) goto 95
	  read(iunit,*) (hl(l),l=1,lmax)
	else
	  write(6,*) 'ts_file_open: opening unformatted file ',name
	  iunit = ifileo(0,name,'unform','old')
	  read(iunit) it,nknaux,lmax,nvar
	  if( lmax .gt. nlvdim ) goto 95
	  read(iunit,iostat=ios) (hl(l),l=1,lmax)
	end if
	close(iunit)

	if( ios .ne. 0 ) goto 97

	inzeta = 0
	ihashl = 1
	do l=2,lmax
	  if( hl(l) .gt. hl(l-1) ) then		!in zeta
	    inzeta = +1
	  else if( hl(l) .lt. hl(l-1) ) then	!in sigma
	    if( inzeta .gt. 0 ) ihashl = 0
	  else
	    ihashl = 0
	  end if
	end do
        if( lmax .le. 1 ) ihashl = 0

	!ihashl = 0			!force choice if you are sure (0/1)
	bhashl = ihashl .gt. 0

c-------------------------------------------------------------
c if level structure is available, check if compatible
c-------------------------------------------------------------

	!if( bhashl ) then
	!  do l=1,lmax
	!    if( hl(l) .ne. hlv(l) ) goto 96
	!  end do
	!end if

c-------------------------------------------------------------
c open finally and store info
c-------------------------------------------------------------

	if( bformat ) then
	  iunit = ifileo(0,name,'form','old')
	else
	  iunit = ifileo(0,name,'unform','old')
	end if

	info(1) = iunit
	info(2) = iformat
	info(3) = ihashl

	if( bdebug ) then
	  write(6,*) 'debug in ts_file_open:'
	  write(6,*) 'iunit,iformat,ihashl: ',iunit,iformat,ihashl
	  write(6,*) 'lmax,hl: ',lmax
	  write(6,*) (hl(l),l=1,lmax)
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   95	continue
	write(6,*) 'too many layers in file: ',name
	write(6,*) 'lmax,nlvdim  : ',lmax,nlvdim
	write(6,*) 'it,nknaux,nvar: ',it,nknaux,nvar
	write(6,*) 'iformat,ihashl: ',iformat,ihashl
	stop 'error stop ts_file_open: dimension nlvdim'
   96	continue
	write(6,*) 'incompatible levels in file: ',name
	write(6,*) 'lmax,hl  : ',lmax,(hl(l),l=1,lmax)
	write(6,*) 'lmax,hlv : ',lmax,(hlv(l),l=1,lmax)
	stop 'error stop ts_file_open: incompatible levels'
   97	continue
	write(6,*) 'read error in file: ',name
	stop 'error stop ts_file_open: read error'
   98	continue
	write(6,*) 'empty file: ',name
	stop 'error stop ts_file_open: empty file'
   99	continue
	write(6,*) 'file name: ',name
	stop 'error stop ts_file_open: cannot open file'
	end

c*******************************************************************	

	subroutine ts_next_record_0(it,info,nkn,nlv,value)

	implicit none

	include 'param.h'

	integer it
	integer info(3)
	integer nkn
	integer nlv
	real value(nlvdim,1)

	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hlv(1)
	common /hlv/hlv
	real znv(1)
	common /znv/znv
	real hkv(1)
	common /hkv/hkv

	logical bsigma,bdebug,bcons
	logical bformat,bhashl
	integer nknaux,lmax,nvar
	integer i,l,iunit
	integer nsigma,nsigma_aux
	real hsigma,hsigma_aux
	real val
        real vmin,vmax

	integer lmax_fem
	real htot
	real zz_fem,zz_data
	real hlv_aux(nlvdim)
        real hl_fem(0:nlvdim+1)
        real hl_data(0:nlvdim+1)
        real val_fem(nlvdim+1)
        real val_data(nlvdim+1)

	bdebug = .false.
	bdebug = .true.

	bcons = .false.		!conserve total quantity?

	iunit   = info(1)
	bformat = info(2) .gt. 0
	bhashl  = info(3) .gt. 0

	if( iunit .le. 0 ) return

	if( bdebug ) then
	  write(6,*) 'debug in ts_next_record:'
	  write(6,*) 'iunit,bformat,bhashl: ',it,iunit,bformat,bhashl
	end if

	if( bformat ) then
	  read(iunit,*) it,nknaux,lmax,nvar
	else
	  read(iunit) it,nknaux,lmax,nvar
	end if

	write(6,*)'reading T/S values', it,nknaux,lmax,nvar

	if( nkn .ne. nknaux ) stop 'error stop ts_next_record: nkn'
	if( nvar .ne. 1 ) stop 'error stop ts_next_record: nvar'
	if( lmax .gt. nlvdim ) stop 'error stop ts_next_record: nlvdim'

	if( bformat ) then
	  if( bhashl ) read(iunit,*) (hlv_aux(l),l=1,lmax)
	  read(iunit,*) ((value(l,i),l=1,lmax),i=1,nkn)
	else
	  if( bhashl ) read(iunit) (hlv_aux(l),l=1,lmax)
	  read(iunit) ((value(l,i),l=1,lmax),i=1,nkn)
	end if

	if( nlv .gt. lmax ) then
	  do i=1,nkn
	    val = value(lmax,i)
	    do l=lmax+1,nlv
	      value(l,i) = val
	    end do
	  end do
	end if

	write(6,*) '****** no interpolation done ********'
	return

c--------------------------------------------------------------
c interpolate between different vertical structures 
c--------------------------------------------------------------

c the following still has to be checked

        call get_sigma(nsigma,hsigma)		!from basin
        bsigma = nsigma .gt. 0

	call compute_sigma_info(lmax,hlv_aux,nsigma_aux,hsigma_aux)

        do i = 1,nkn

          lmax_fem = ilhkv(i)

          hl_data(0)= -znv(i)
          hl_fem(0) = -znv(i)

          if(nsigma.gt.0.and.nsigma.ge.lmax_fem)then !sigma
            htot = hkv(i) !giusto
          else !zeta e ibrido
            htot = hlv(lmax_fem) !questo non risolve il problema
          endif

          zz_data = znv(i)
          zz_fem = znv(i)

          call set_hybrid_depth(lmax_fem,zz_fem,htot
     +			,hlv,nsigma,hsigma,hl_fem(1))

          htot = hlv_aux(lmax) !questo non risolve il problema

          call set_hybrid_depth(lmax,zz_data,htot
     +			,hlv_aux,nsigma_aux,hsigma_aux,hl_data(1))

          do l=1,lmax
            val_data(l) = value(l,i)
          end do

          call intp_vert(bcons,lmax,hl_data,val_data
     +				,lmax_fem,hl_fem,val_fem)

          do l = 1,lmax_fem
            value(l,i) = val_fem(l) !DEB
          end do
        enddo

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

        call conmima(nlvdim,value,vmin,vmax)
        write(6,*) 'min/max: ',vmin,vmax

	end

c*******************************************************************	

	subroutine ts_file_close_0(info)

c closes T/S file

	implicit none

	integer info(3)

	close(info(1))

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

