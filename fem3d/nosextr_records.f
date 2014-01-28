c
c $Id: nosextr_records.f,v 1.1 2008-07-16 15:41:39 georg Exp $
c
c extract records from NOS file
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 03.06.2011    ggu     routine adjourned
c
c**********************************************************

	program nosextr_records

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
        integer nkn,nel,nlv,nvar
        integer it,ivar
        integer ierr
        integer nvers
        real r,rnull
	real conz,high

        integer iapini,ideffi,ifileo

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

	nread=0
	nextr=0
	rnull=0.
        ivarold = 0
	high = 1.e+30

c-------------------------------------------------------------------
c get simulation
c-------------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c-------------------------------------------------------------------
c open NOS file and read header
c-------------------------------------------------------------------

	call open_nos_type('.nos','old',nin)

        call read_nos_header(nin,nkndim,neldim,nlvdim,ilhkv,hlv,hev)
        call nos_get_params(nin,nkn,nel,nlv,nvar)

c-------------------------------------------------------------------
c get records to extract from STDIN
c-------------------------------------------------------------------

	call get_records_from_stdin(nrdim,irec,ball)

c-------------------------------------------------------------------
c open NOS output file
c-------------------------------------------------------------------

	!call mkname(' ','nos_extract','.nos',file)
	!write(6,*) 'writing file ',file(1:50)
	!nb = ifileo(55,file,'unform','new')

	call open_nos_file('nos_extract','new',nb)

	call nos_init(nb,0)
	call nos_clone_params(nin,nb)
	call write_nos_header(nb,ilhkv,hlv,hev)

	!call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
	!if( ierr .ne. 0 ) goto 99
	!call wsnos(nb,ilhkv,hlv,hev,ierr)
	!if( ierr .ne. 0 ) goto 99

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

	call nos_read_record(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100
        if( ivarold .eq. 0 ) ivarold = ivar
        if( ivar .ne. ivarold ) goto 91

	nread=nread+1
	if( .not. ball .and. nread .gt. nrdim ) goto 100
	write(6,*) 'time : ',nread,it,ivar

	bwrite = ball .or. irec(nread) .ne. 0

	if( bwrite ) then
	  call nos_write_record(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	  if( ierr .ne. 0 ) goto 99
	  nextr = nextr + 1
	end if

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nextr,' records written to file nos_extract.nos'
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
	stop 'error stop nosextr_records: ivar'
   98	continue
	write(6,*) 'error opening outout file'
	stop 'error stop nosextr_records'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nosextr_records'
	end

c***************************************************************

