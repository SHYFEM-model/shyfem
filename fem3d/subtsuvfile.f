c $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
c
c reading and interpolation of external files
c
c contents :
c
c revision log :
c
c 29.10.2012    ggu     created from scratch
c 17.06.2013    ggu     do not pass function into subroutine
c 02.07.2014    ggu     new framework finished
c
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open(name,it,nkn,nlv,iunit)
	integer iunit(3)
	character*(*) name
	imreg = nint(getpar('imreg'))
	if( imreg .eq. 3 ) then
	  call ts_file_open_1(name,it,nkn,nlv,iunit)
	else
	  call ts_file_open_0(name,it,nkn,nlv,iunit)		!old call
	end if
	end

	subroutine ts_file_descrp(iunit,name)
	use intp_fem_file
	integer iunit(3)
	character*(*) name
	imreg = nint(getpar('imreg'))
	if( imreg .eq. 3 ) then
	  call iff_set_description(iunit(1),0,name)
	end if
	end

	subroutine ts_next_record(it,iunit,nkn,nlv,value)
	include 'param.h'
	integer iunit(3)
	real value(nlvdim,1)
	imreg = nint(getpar('imreg'))
	if( imreg .eq. 3 ) then
	  call ts_next_record_1(it,iunit,nkn,nlv,value)
	else
	  call ts_next_record_0(it,iunit,nkn,nlv,value)	!old call
	end if
	end

	subroutine ts_file_close(info)
	integer info(3)
	imreg = nint(getpar('imreg'))
	if( imreg .eq. 3 ) then
	  call ts_file_close_1(info)
	else
	  call ts_file_close_0(info)			!old call
	end if
	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open_1(file,it,np,nlv,iunit)

c opens T/S file

	use intp_fem_file

	implicit none

	include 'param.h'

	character*(*) file		!name of file
	integer it
	integer np			!number of points expected
	integer nlv
	integer iunit(3)		!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer id
	double precision dtime
	integer nodes(1)
	real vconst(1)

	dtime = it
	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

	call iff_init(dtime,file,nvar,nexp,lexp,nintp
     +                                  ,nodes,vconst,id)

	iunit(1) = id

	return
   99	continue
	write(6,*) 'Cannot open file: ',file
	stop 'error stop ts_file_open: error file open'
	end

c*******************************************************************	

	subroutine ts_next_record_1(it,iunit,nkn,nlv,value)

	use intp_fem_file

	implicit none

	include 'param.h'

	integer it
	integer iunit(3)
	integer nkn
	integer nlv
	real value(nlvdim,nkn)

	integer id,ldim,ndim,ivar
        real vmin,vmax
	double precision dtime
	character*80 string

c--------------------------------------------------------------
c initialize
c--------------------------------------------------------------

	id = iunit(1)

c--------------------------------------------------------------
c read new data
c--------------------------------------------------------------

	dtime = it
	ivar = 1
	ndim = nkn
	ldim = nlvdim

	write(6,*)'reading T/S values', it,dtime

	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

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

	subroutine ts_file_close_1(iunit)

c closes T/S file

	use intp_fem_file

	implicit none

	integer iunit(3)

	integer id

	id = iunit(1)
	call iff_forget_file(id)

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	
c****** old routines ***********************************************	
c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_open_0(name,itaux,nkn,nlv,info)

c opens T/S file

	implicit none

	include 'param.h'

	character*(*) name
	integer itaux
	integer nkn
	integer nlv
	integer info(3)

	real hlv(1)
	common /hlv/hlv

	logical bsigma,bdebug
	logical bformat,bhashl
	integer ifileo,iunit,ios
	integer iformat,ihashl,inzeta
	integer nknaux,lmax,nvar
	integer l,it
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

