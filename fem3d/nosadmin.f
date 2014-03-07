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

	program nosadmin

c administrates some information of NOS file and re-writes it to new one

        implicit none

	include 'param.h'

        integer ndim
	parameter ( ndim = 10000 )

	integer iu(ndim)
	character*80 name,file

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
        integer nread,nextr
        integer l,k,nin,nb
        integer nkn,nel,nlv,nvar
        integer it,ivar,idt
        integer ierr
        integer nvers
	integer ivarold,ivarnew
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
c open NOS output file
c-------------------------------------------------------------------

	call open_nos_file('nos_new','new',nb)

	call nos_init(nb,0)
	call nos_clone_params(nin,nb)
	call change_data(nb,ivarold,ivarnew,idt)
	call write_nos_header(nb,ilhkv,hlv,hev)

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

	call nos_read_record(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread = nread + 1
	write(6,*) 'time : ',nread,it,ivar

	it = it + idt
	if (ivar .eq. ivarold ) ivar = ivarnew

	call nos_write_record(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	if( ierr .ne. 0 ) goto 99

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nread,' records written to file nos_new.nos'
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

	subroutine change_data(nb,ivarold,ivarnew,idt)

	implicit none

	integer nb
	integer ivarold,ivarnew
	integer idt

	character*80 line,title
	integer n,i
	integer date,time
	double precision d(2)

	integer iscand

	call nos_get_title(nb,title)
        write(6,*)
        write(6,*) 'Enter new title: '
        read(5,'(a)') line
        if( line .ne. ' ' ) title = line
	call nos_set_title(nb,title)
        write(6,*) 'title :',title
        write(6,*)

	call nos_get_date(nb,date,time)
        write(6,*)
        write(6,*) 'Enter new date and time: '
        read(5,'(a)') line
        n = iscand(line,d,2)
        if( n .lt. 0 .or. n .gt. 2 ) goto 95
        if( n .gt. 0 ) date = nint(d(1))
        if( n .gt. 1 ) time = nint(d(2))
	call nos_set_date(nb,date,time)
        write(6,*) 'date,time :',date,time
        write(6,*)

	idt = 0
        write(6,*)
        write(6,*) 'Enter time shift for it (in seconds): '
        read(5,'(a)') line
        n = iscand(line,d,1)
        if( n .lt. 0 .or. n .gt. 1 ) goto 95
        if( n .gt. 0 ) idt = nint(d(1))
        write(6,*) 'idt :',idt
        write(6,*)

	ivarold = 0
	ivarnew = 0
        write(6,*)
        write(6,*) 'Enter ivarold,ivarnew to change var id: '
        read(5,'(a)') line
        n = iscand(line,d,2)
        if( n .lt. 0 .or. n .gt. 2 ) goto 95
        if( n .eq. 1 ) goto 95
        ivarold = nint(d(1))
        ivarnew = nint(d(2))
        write(6,*) 'ivarold,ivarnew :',ivarold,ivarnew
        write(6,*)

	return
   95   continue
        write(6,*) n,(d(i),i=1,n)
        write(6,*) line
        stop 'error stop basbathy: error in parameters'
	end

c***************************************************************

