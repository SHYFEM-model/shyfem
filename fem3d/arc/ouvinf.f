c
c $Id: ouvinf.f,v 1.5 2000/03/02 09:31:22 georg Exp $

	program ouvinf

c info on OUV file

	implicit none

	integer nkndim,neldim,nlvdim
	parameter(nkndim=7000,neldim=2*nkndim,nlvdim=15)

	integer it
	integer iunit,nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtouv
	real href,hzmin,hlvmin
	real hlv(nlvdim),hev(neldim)
	integer ilhv(neldim)
	real z(nkndim),u(nlvdim,neldim),v(nlvdim,neldim)
	character*80 title

	integer ier
	integer rfouv,rdouv
	integer iapini
	integer ideffi

	iunit = 5
	nvmax = 2

        if(iapini(2,nkndim,neldim,0).eq.0) then
                stop 'error stop : iapini'
        end if

        iunit=ideffi('datdir','runnam','.ouv','unform','old')
        if(iunit.le.0) goto 100

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

    1	continue

	  ier = rdouv(iunit,it,nlvdim,z,u,v,ilhv)
	  if( ier .eq. -1 ) goto 2
	  if( ier .gt. 0 ) goto 99

	  write(6,*) it

	goto 1

    2	continue

	stop
   98	continue
	stop 'error stop: header record'
   99	continue
	stop 'error stop: data record'
  100	continue
	stop 'error stop: Cannot open file'
	end

c********************************************************************

