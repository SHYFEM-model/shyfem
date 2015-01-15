c
c $Id: splitets.f,v 1.11 2009-02-13 17:22:44 georg Exp $
c
c revision log :
c
c 27.01.2014	ggu	copied from splitext
c
c***************************************************************

	program splitets

c This routine reads an EXT file and writes the data to unit 76

	implicit none

	include 'param.h'

	integer ndim
	parameter (ndim=100)

	character*80 line,file,format
	integer nrec,it,i,nin,j,lmax,l
	integer nkn,nlv,nvar
	integer ivar,ierr,iunit

        integer ilets(nexdim)
        real hlv(nlvdim)
        real hets(nexdim)
        integer nodes(nexdim)
        real xg(nexdim)
        real yg(nexdim)
        character*80 desc(nexdim)

	real cv3(nlvdim,nexdim)
	real cv(nexdim)
	integer ivars(ndim)

	integer iapini, ifemop

c---------------------------------------------------------------

	nrec = 0

c---------------------------------------------------------------
c open files
c---------------------------------------------------------------

        call shyfem_copyright('etsinf - Info on ETS files')

        if(iapini(2,nkndim,neldim,0).eq.0) then
                stop 'error stop : iapini'
        end if

        call open_ets_type('.ets','old',nin)

        call read_ets_header(nin,nexdim,nlvdim,ilets,hlv,hets
     +                                  ,nodes,xg,yg,desc)
        call ets_get_params(nin,nkn,nlv,nvar)

	if( nvar .gt. ndim ) goto 95
	call ets_get_vars(nin,nvar,ivars)

        write(6,*) 'Available variables: ',nvar
        write(6,*) (ivars(i),i=1,nvar)

c---------------------------------------------------------------
c open output files
c---------------------------------------------------------------

	  do i=1,nkn
	    do j=1,nvar
	      ivar = ivars(j)
	      call open_file(ivar,i)
	    end do
	  end do

c---------------------------------------------------------------
c loop over data
c---------------------------------------------------------------

	do while(.true.)

	  call ets_read_record(nin,it,ivar,nlvdim,ilets,cv3,ierr)

	  if( ierr .gt. 0 ) write(6,*) 'error in reading file : ',ierr
	  if( ierr .ne. 0 ) goto 100

	  nrec = nrec + 1
	  if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'

	  do i=1,nkn
	    lmax = ilets(i)
	    if( ivar .eq. 1 ) lmax = 1
	    if( ivar .eq. 31 ) lmax = 4
	    call get_unit(ivar,i,iunit)
	    write(format,'(a,i4,a)') '(i12,',lmax,'(f12.4))'
	    write(iunit,format) it,(cv3(l,i),l=1,lmax)
	    !write(iunit,'(i12,5(f12.4))') it,(cv3(l,i),l=1,lmax)
	    !write(iunit,'(i12,*(f12.4))') it,(cv3(l,i),l=1,lmax)
	  end do

	end do

c---------------------------------------------------------------
c end of loop on data
c---------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   95   continue
        write(6,*) 'nvar = ',nvar,'   ndim = ',ndim
        stop 'error stop: ndim'
	end

c*******************************************************************

	subroutine opents(iunit,name,number)

c opens file name.number on unit iunit

	implicit none

	integer iunit
	character*(*) name
	integer number

	integer in,nout
	character*80 numlin,file
	integer ialfa

	nout = iunit
	in = ialfa(float(number),numlin,-1,-1)
	file = name // '.' // numlin(1:in)

	open(nout,file=file,status='unknown',form='formatted')

	end

c*******************************************************************

	subroutine wrts(n,it,data,name,number)

c writes data to file name.number

	implicit none

	integer n
	integer it(1)
	real data(1)
	character*(*) name
	integer number

	integer i
	integer in,nout
	character*80 numlin,file
	integer ialfa

	nout = 1
	call opents(nout,name,number)

	do i=1,n
	  write(nout,*) it(i),data(i)
	end do

	close(nout)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine get_unit(ivar,inode,iunit)

	implicit none

	integer ivar,inode,iunit

	iunit = 100*ivar + inode

	end

c*******************************************************************

	subroutine open_file(ivar,inode)

	implicit none

	integer ivar,inode

	integer i,iunit

	i = inode

	call get_unit(ivar,inode,iunit)

	if( ivar .eq. 1 ) then
	  call opents(iunit,'z',i)
	else if( ivar .eq. 31 ) then
	  call opents(iunit,'w',i)
	else if( ivar .eq. 6 ) then
	  call opents(iunit,'u',i)
	else if( ivar .eq. 7 ) then
	  call opents(iunit,'v',i)
	else if( ivar .eq. 8 ) then
	  call opents(iunit,'m',i)
	else if( ivar .eq. 11 ) then
	  call opents(iunit,'s',i)
	else if( ivar .eq. 12 ) then
	  call opents(iunit,'t',i)
	else
	  write(6,*) 'variable not recognized: ivar = ',ivar
	  stop 'error stop open_file: ivar'
	end if

	end

c*******************************************************************

