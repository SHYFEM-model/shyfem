c
c $Id: subinterpol1.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 27.02.2014    ggu     subroutines copied from basbathy
c
c****************************************************************

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine check_basin_name(name,bbasin)

	implicit none

	character*(*) name
	logical bbasin

	integer i,nend,ndot

	do i=1,len(name)
	  if( name(i:i) .eq. ' ' ) exit
	end do

	nend = max(1,i-1)
	ndot = max(1,nend-3)
	
	if( name(ndot:nend) .eq. '.grd' ) then
	  bbasin = .false.
	else if( name(ndot:nend) .eq. '.bas' ) then
	  bbasin = .true.
	else
	  write(6,*) 'name = ',name
	  write(6,*) 'expecting type of .grd or .bas'
	  stop 'error stop check_basin_name: cannot parse name'
	end if

	end

c*******************************************************************

	subroutine read_basin(name)

	use mod_depth
	use basin

	implicit none

	character*(*) name
	integer nknddi,nelddi

	include 'param.h'

	integer ie,ii
	real h

	open(1,file=name,status='old',form='unformatted')
	call basin_read(1)
	close(1)

	do ie=1,nel
	  h = 0.
	  do ii=1,3
	    h = h + hm3v(ii,ie)
	  end do
	  hev(ie) = h / 3.
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine write_xy(nfile,nkn,ipv,xgv,ygv)

c writes xy as data to file

	implicit none

	character*(*) nfile
	integer nkn
	integer ipv(nkn)
	real xgv(nkn)
	real ygv(nkn)

	integer k

	open(1,file=nfile,status='unknown',form='formatted')
        write(1,*) nkn
	do k=1,nkn
	  write(1,*) k,ipv(k),xgv(k),ygv(k)
	end do
	close(1)
        write(6,*) 'The coordinates have been written to ',nfile

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine write_data(nfile,nkn,hkv)

c writes depth as data to file

	implicit none

	character*(*) nfile
	integer nkn
	real hkv(nkn)

	integer k

	open(1,file=nfile,status='unknown',form='formatted')
        write(1,*) 0,nkn,1,1
	do k=1,nkn
	  write(1,*) hkv(k)
	end do
	close(1)
        write(6,*) 'The data has been written to ',nfile

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine write_nos(nfile,nkn,nel,ivar,c,hev)

	implicit none

	character*(*) nfile
	integer nkn,nel,ivar
	real c(1)
	real hev(1)

	integer iu,it,nlv,nlvdi,nvers,nvar
	integer ierr
	integer ilhkv(1)
	real hlv(1)
	character*80 title,femver

	iu = 1
	title = 'interpolated field with scalintp'
	nvar = 1
	it = 0
	nlv = 1
	nlvdi = 1
        nvers = 5
	call get_shyfem_version(femver)

	open(iu,file=nfile,status='unknown',form='unformatted')

        call nos_init(iu,nvers)
        call nos_set_title(iu,title)
        !call nos_set_date(iu,date,time)
        call nos_set_femver(iu,femver)
        call nos_write_header(iu,nkn,nel,nlv,nvar,ierr)
        if(ierr.gt.0) goto 99
        call nos_write_header2(iu,ilhkv,hlv,hev,ierr)
        if(ierr.gt.0) goto 99

        call nos_write_record(iu,it,ivar,nlvdi,ilhkv,c,ierr)
        if(ierr.gt.0) goto 99

	close(iu)

        write(6,*) 'The data has been written to ',nfile

	return
   99   continue
        write(6,*) 'error ',ierr,' writing file with name ',nfile
        stop 'error stop write_nos'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine transfer_depth(ike,ht)

c copies depth values from elems/nodes to nodes/elems

	use mod_depth
	use basin

	implicit none

	integer ike
	real ht(1)	!this is interpolated depth

	include 'param.h'

	integer k,ie,ii
	real depth
	integer ic(nkn)

	if( ike .eq. 1 ) then		!elementwise

	  do ie=1,nel
	    hev(ie) = ht(ie)
	  end do

	  do k=1,nkn
	    ic(k) = 0
	    hkv(k) = 0.
	  end do

	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + hev(ie)
	      ic(k) = ic(k) + 1			!BUG - this was missing
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / ic(k)
	  end do

	else				!nodewise

	  do k=1,nkn
	    hkv(k) = ht(k)
	  end do

	  do ie=1,nel
	    depth = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      depth = depth + hkv(k)
	    end do
	    hev(ie) = depth / 3.
	  end do

	end if

	end

c*******************************************************************

	subroutine wrgrd(iunit,ike)

c writes grd file from bas

	use mod_depth
	use basin

	implicit none

	integer iunit
	integer ike

	include 'param.h'





	integer k,ie,ii

	do k=1,nkn
	  if( ike .eq. 1 ) then
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	  else
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k),hkv(k)
	  end if
	end do

	write(iunit,*)

	do ie=1,nel
	  if( ike .eq. 1 ) then
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3),hev(ie)
	  else
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3)
	  end if
	end do

	return
 1000	format(i1,2i10,3e16.8)
 1100	format(i1,2i10,i4,3i10,e16.8)
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine node_test

	use basin

	implicit none

	include 'param.h'


	logical bstop
	integer ie,ii,iii,k,k1

	bstop = .false.

	write(6,*) 'node_test ... ',nel,nkn
	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
            if( k .le. 0 ) then
                write(6,*) ie,ii,k
                bstop = .true.
	    end if
	    iii = mod(ii,3) + 1
	    k1 = nen3v(iii,ie)
            if( k .eq. k1 ) then
                write(6,*) ie,(nen3v(iii,ie),iii=1,3)
                bstop = .true.
            end if
	  end do
	end do
	write(6,*) 'end of node_test ... '

	if( bstop ) stop 'error stop node_test: errors'

	end

c*******************************************************************

	subroutine readgrd(file,np,xp,yp,dp)

c reads bathymetry file

	implicit none

	character*(*) file
	integer np
	real xp(1)
	real yp(1)
	real dp(1)

	integer ndim,n
	integer id,inum,ityp
	real x,y,d

	ndim = np

	n = 0
	!write(6,*) 'reading bathymetry file : ',file
	open(1,err=99,file=file,status='old',form='formatted')

    1	continue
	  read(1,*,end=2) id,inum,ityp,x,y,d
	  if( id .ne. 1 ) goto 98
	  if( inum .le. 0 ) goto 98
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop: ndim'
	  xp(n) = x
	  yp(n) = y
	  dp(n) = d
	  goto 1
    2	continue

	close(1)
	np = n
	write(6,*) np,' data points read'

        return
   98   continue
        write(6,*) 'error reading line:'
        write(6,*) 'id = ',id,'   inum = ',inum
        stop 'error stop readbat: error reading line'
   99   continue
        write(6,*) 'Cannot read file: ',file
        stop 'error stop readbat: no such file'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

