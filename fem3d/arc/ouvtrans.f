c
c $Id: ouvtrans.f,v 1.1 2000/05/26 11:38:36 georg Exp $

	program ouvtrans

c info on OUV file

	implicit none

	integer nkndim,neldim,nlvdim
	parameter(nkndim=7000,neldim=2*nkndim,nlvdim=15)

	integer it
	integer iunit,iunew
	integer nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtouv
	real href,hzmin,hlvmin
	real hlv(nlvdim),hev(neldim)
	integer ilhv(neldim)
	real z(nkndim),u(nlvdim,neldim),v(nlvdim,neldim)
	character*80 title

	integer irec
	integer ier
	integer rfouv,rdouv
	integer wfouv,wrouv
	integer iapini
	integer ideffi

	iunew = 66
	nvmax = 2

        if(iapini(2,nkndim,neldim,0).eq.0) then
                stop 'error stop : iapini'
        end if

        iunit=ideffi('datdir','runnam','.ouv','unform','old')
        if(iunit.le.0) goto 100

	open(iunew,file='new.ouv',status='unknown',form='unformatted')

	nkn = nkndim
	nel = neldim
	nlv = nlvdim

	ier = rfouv(iunit,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtouv
     +				,href,hzmin,hlvmin
     +				,hlv,hev,ilhv
     +				,title)

	if( ier .eq. -1 ) goto 2
	if( ier .gt. 0 ) goto 98

	write(6,*) 'nvers : ',nvers
	write(6,*) 'nkn,nel,nlv : ',nkn,nel,nlv
	write(6,*) 'itanf,itend,idt,idtouv : ',itanf,itend,idt,idtouv
	write(6,*) 'href,hzmin,hlvmin : ',href,hzmin,hlvmin
	write(6,*) 'title : ', title

	ier = wfouv(iunew,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtouv
     +				,href,hzmin,hlvmin
     +				,hlv,hev,ilhv
     +				,title)

	if( ier .ne. 0 ) goto 91

	irec = 0

c--------------------------------------------------------
c loop on data records
c--------------------------------------------------------

    1	continue

	  ier = rdouv(iunit,it,nlvdim,z,u,v,ilhv)
	  if( ier .eq. -1 ) goto 2
	  if( ier .gt. 0 ) goto 99

	  irec = irec + 1
	  write(6,*) irec,it

c -> qui decidi tu cosa vuoi fare
c per esempio scrivere ogni secondo record:

	  if( mod(irec,2) .eq. 1 ) goto 1

	  write(6,*) '        ...writing'
	  ier = wrouv(iunew,it,nlvdim,z,u,v,ilhv)
	  if( ier .ne. 0 ) goto 92

	goto 1

c--------------------------------------------------------
c end of loop on data records
c--------------------------------------------------------

    2	continue

	stop
   91	continue
	stop 'error stop: writing header record'
   92	continue
	stop 'error stop: writing data record'
   98	continue
	stop 'error stop: reading header record'
   99	continue
	stop 'error stop: reading data record'
  100	continue
	stop 'error stop: Cannot open file'
	end

c********************************************************************

