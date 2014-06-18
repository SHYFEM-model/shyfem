c
c $Id: nosextr_gis.f,v 1.2 2008-11-20 10:51:34 georg Exp $
c
c extract reords from NOS file in gis format
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 13.11.2008    ggu     write all records with ball=.true.
c 23.11.2010    ggu     new routine wrgis_3d() to extract 3d info
c 03.06.2011    ggu     routine adjourned
c
c**********************************************************

	program nosextr_gis

c extracts whole records from nos file and writes in gis format
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

	ball = .false.

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

	call get_records_from_stdin(nrdim,irec,ball)

c-------------------------------------------------------------------
c open GIS output file
c-------------------------------------------------------------------

	call mkname(' ','nos_extract','.gis',file)
	write(6,*) 'writing file ',file(1:50)
	nb = ifileo(55,file,'form','new')
	if( nb .le. 0 ) goto 98

	call write_connections(nkn,nel,nen3v,xgv,ygv)

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
        if( .not. ball .and. nread .gt. nrdim ) goto 100
	write(6,*) 'time : ',nread,it,ivar

	bwrite = ball .or. irec(nread) .ne. 0

	if( bwrite ) then
	  !call wrgis_surf(nb,it,ivar,nkn,cv3)
	  !call wrgis_sep(it,ivar,nkn,cv3)
	  call wrgis_3d(nb,it,ivar,nkn,ilhkv,cv3)
	  nextr = nextr + 1
	end if

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nextr,' records written to file nos_extract.gis'
	!write(6,*) nextr,' files nos_extract__*.gis created'
	write(6,*)

        if( nextr .le. 0 ) stop 'no file written'

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   91	continue
	write(6,*) 'file may have only one type of variable'
	write(6,*) 'error ivar : ',ivar,ivarold
	write(6,*) 'You should use nossplit to extract scalars first'
	stop 'error stop nosextr_gis: ivar'
   94	continue
	write(6,*) 'incompatible simulation and basin'
	write(6,*) 'nkn: ',nkn,nkn1
	write(6,*) 'nel: ',nel,nel1
	stop 'error stop nosextr_gis: nkn,nel'
   98	continue
	write(6,*) 'error opening output file'
	stop 'error stop nosextr_gis'
	end

c***************************************************************

	subroutine wrgis_3d(nb,it,ivar,nk,ilhkv,cv3)

c writes one record to file nb (3D)

	implicit none

	include 'param.h'

	integer nb,it,ivar,nk
	integer ilhkv(nlvdim)
	real cv3(nlvdim,nkndim)

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,l,lmax
	real x,y

	write(nb,*) it,nk,ivar

	do k=1,nk
	  lmax = ilhkv(k)
	  x = xgv(k)
	  y = ygv(k)

	  write(nb,*) k,x,y,lmax
	  write(nb,*) (cv3(l,k),l=1,lmax)
	end do

	end

c***************************************************************

	subroutine wrgis_surf(nb,it,ivar,nk,cv3)

c writes one record to file nb (2D)

	implicit none

	include 'param.h'

	integer nb,it,ivar,nk
	real cv3(nlvdim,nkndim)

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,level
	real x,y,c

	level = 1

	write(nb,*) it,nk,ivar

	do k=1,nk
	  x = xgv(k)
	  y = ygv(k)
	  c = cv3(level,k)

	  write(nb,*) x,y,c

	end do

	end

c***************************************************************

	subroutine wrgis_sep(it,ivar,nk,cv3)

c writes one record to single files

	implicit none

	include 'param.h'

	integer it,ivar,nk
	real cv3(nlvdim,nkndim)

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,level,nb
	real x,y,c
	character*60 name

	integer ifileo

	level = 1

        call make_name_with_time('extract','.gis',it,name)
	nb = ifileo(60,name,'form','new')

	write(nb,*) 'x   y   val'

	do k=1,nk
	  x = xgv(k)
	  y = ygv(k)
	  c = cv3(level,k)

	  write(nb,*) x,y,c

	end do

	close(nb)

	end

c***************************************************************

