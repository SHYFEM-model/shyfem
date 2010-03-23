c
c $Id: nosextr_gis.f,v 1.2 2008-11-20 10:51:34 georg Exp $
c
c 18.11.1998    ggu     check dimensions with dimnos
c 13.11.2008    ggu     write all records with ball=.true.
c
c**********************************************************

	program nosextr_gis

c extracts whole records from nos file
c
c records have to be specified on stdin

        implicit none

	include 'param.h'

        integer ndim
	parameter ( ndim = 10000 )

	integer iu(ndim)
	character*80 name,file

        integer nrdim
	parameter ( nrdim = 2000 )

	integer irec(nrdim)

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        integer nen3v(3,neldim)
        integer ipv(nkndim), ipev(neldim)
        integer iarv(neldim)

        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v
        common /nen3v/nen3v
        common /ipv/ipv, /ipev/ipev
        common /iarv/iarv
c--------------------------------------------------

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

        double precision accum(nlvdim,nkndim)
	real amin(nlvdim,nkndim)
	real amax(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical ball,bwrite
        integer nread,ivarold,nextr
        integer l,k,nin,nb
        integer nkn1,nel1,nlv,nvar
        integer it,ivar
        integer ierr
        integer nvers
        real r,rnull
	real conz,high

        integer iapini,ideffi,ifileo

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

	ball = .true.
	ball = .false.		!write all records

	nread=0
	nextr=0
	rnull=0.
        ivarold = 0
	high = 1.e+30

c-------------------------------------------------------------------
c get simulation
c-------------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

c-------------------------------------------------------------------
c open NOS file and read header
c-------------------------------------------------------------------

        nvers=3
	call rfnos(nin,nvers,nkn1,nel1,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn1,nel1
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

	if( nkn .ne. nkn1 .or. nel .ne. nel1 ) goto 94

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c get records to extract from STDIN
c-------------------------------------------------------------------

	if( .not. ball ) then
	  call get_records_from_stdin(nrdim,irec)
	end if

c-------------------------------------------------------------------
c open GIS output file
c-------------------------------------------------------------------

	call mkname(' ','extract','.gis',file)
	write(6,*) 'writing file ',file(1:50)
	nb = ifileo(55,file,'form','new')
	if( nb .le. 0 ) goto 98

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100
        if( ivarold .eq. 0 ) ivarold = ivar
        if( ivar .ne. ivarold ) goto 91

	nread=nread+1
	write(6,*) 'time : ',nread,it,ivar

	bwrite = ball .or. (nread .le. nrdim .and. irec(nread) .ne. 0)

	if( bwrite ) then
	  call wrgis(nb,it,ivar,nkn,cv3)
	  call wrgis_sep(nb,it,ivar,nkn,cv3)
	  nextr = nextr + 1
	end if

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nextr,' records written to file extract.gis'
	write(6,*) nextr,' files extract__*.gis created'
	write(6,*)

        if( nextr .le. 0 ) stop 'no file written'

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   91	continue
	write(6,*) 'file may have only one type of variable'
	write(6,*) 'error ivar : ',ivar,ivarold
	stop 'error stop nosextr_records: ivar'
   94	continue
	write(6,*) 'incompatible simulation and basin'
	write(6,*) 'nkn: ',nkn,nkn1
	write(6,*) 'nel: ',nel,nel1
	stop 'error stop nosextr_records: nkn,nel'
   98	continue
	write(6,*) 'error opening file'
	stop 'error stop nosextr_records'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nosextr_records'
	end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

	subroutine get_records_from_stdin(ndim,irec)

c gets records to extract from stdin

	implicit none

	integer ndim
	integer irec(ndim)

	integer i,ir

	do i=1,ndim
	  irec(i) = 0
	end do

	write(6,*) 'Please enter the record numbers to be extracted.'
	write(6,*) 'Enter every record on a single line.'
	write(6,*) 'Finish with 0 on the last line.'
	write(6,*) 'example:'
	write(6,*) '   5'
	write(6,*) '  10'
	write(6,*) '  15'
	write(6,*) '  0'
	write(6,*) ' '

	do while(.true.)
	  write(6,*) 'Enter record to extract (0 to end): '
	  ir = 0
	  read(5,'(i10)') ir

	  if( ir .le. 0 ) then
	    return
	  else if( ir .gt. ndim ) then
	    write(6,*) 'Cannot extract records higher than ',ndim
	    write(6,*) 'Please change ndim and recompile.'
	  else
	    irec(ir) = 1
	  end if
	end do

	end

c***************************************************************

	subroutine wrgis(nb,it,ivar,nkn,cv3)

	implicit none

	include 'param.h'

	integer nb,it,ivar,nkn
	real cv3(nlvdim,nkndim)

        double precision x0,y0
        parameter ( x0 = 2330000.-50000., y0 = 5000000. )

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,level
	real x,y,c

	level = 1

	write(nb,*) it,nkn,ivar

	do k=1,nkn
	  x = xgv(k) + x0
	  y = ygv(k) + y0
	  c = cv3(level,k)

	  write(nb,*) x,y,c

	end do

	end

c***************************************************************

	subroutine wrgis_sep(nb,it,ivar,nkn,cv3)

	implicit none

	include 'param.h'

	integer nb,it,ivar,nkn
	real cv3(nlvdim,nkndim)

        double precision x0,y0
        parameter ( x0 = 2330000.-50000., y0 = 5000000. )

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,level,nbloc
	real x,y,c
	character*60 name

	integer ifileo

	level = 1

        call make_name('extract','.gis',it,name)
	nbloc = ifileo(60,name,'form','new')

	write(nbloc,*) 'x   y   val'

	do k=1,nkn
	  x = xgv(k) + x0
	  y = ygv(k) + y0
	  c = cv3(level,k)

	  write(nbloc,*) x,y,c

	end do

	close(nbloc)

	end

c***************************************************************

        subroutine make_name(pre,post,it,name)

        implicit none

        character*(*) pre,post,name
        integer it

        integer i
        character*10 naux

        write(naux,'(i10)') it

        do i=1,10
          if( naux(i:i) .eq. ' ' ) naux(i:i) = '_'
        end do

        name = pre // naux // post

        end

c*************************************************************

