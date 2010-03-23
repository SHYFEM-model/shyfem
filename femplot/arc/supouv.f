c
c $Id: supouv.f,v 1.4 2004/09/28 09:15:09 georg Exp $

c******************************************************

	subroutine ouvini

	implicit none

	integer nunit,nvers
	common /ouvouv/ nunit,nvers
	integer nknouv,nelouv,nlvouv
	common /ouvnnn/ nknouv,nelouv,nlvouv
	save /ouvouv/, /ouvnnn/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0
	nvers = 0

	nknouv = 0
	nelouv = 0
	nlvouv = 0

	end

c******************************************************

	subroutine ouvclose

	implicit none

	integer nunit,nvers
	common /ouvouv/ nunit,nvers

	if( nunit .gt. 0 ) close(nunit)

	end

c******************************************************

	subroutine ouvopen

	implicit none

	integer nunit,nvers
	common /ouvouv/ nunit,nvers
	integer nknouv,nelouv,nlvouv
	common /ouvnnn/ nknouv,nelouv,nlvouv

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
	character*80 descrp
        common /descrp/ descrp

        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhv(1)
        common /ilhv/ilhv

        character*80 file
	integer itanf,itend,idt,idtouv
	integer ierr
	real href,hzmin,hlvmin

	integer rfouv,rsouv,ideffi

	call ouvini

	nunit = ideffi('datdir','runnam','.ouv','unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop ouvopen: cannot open OUT file'
	end if

        write(6,*) 'file OUV opened'
        inquire(nunit,name=file)
        write(6,*) 'Reading file ...'
        write(6,*) file

        ierr = rfouv(
     +                           nunit,nvers
     +                          ,nknouv,nelouv,nlvouv
     +                          ,itanf,itend,idt,idtouv
     +                          ,href,hzmin,hlvmin
     +                          ,descrp
     +		    )

	if( ierr .ne. 0 ) then
		stop 'error stop ouvopen: error reading header'
	end if

	if( nknouv .gt. nkn ) goto 99
	if( nelouv .gt. nel ) goto 99
	if( nlvouv .gt. nlv ) goto 99

	nlv = nlvouv

        call putpar('href',href)
        call putpar('hzmin',hzmin)
        call putpar('hlvmin',hlvmin)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' itanf =',itanf,'   itend =',itend
        write(6,*) '   idt =',idt,  '  idtouv =',idtouv
        write(6,*)
        write(6,*) ' nvers =',nvers
        write(6,*) '  href =',href,  '  hzmin =',hzmin
        write(6,*) 'hlvmin =',hlvmin,'  hzmin =',hzmin
        write(6,*)
        write(6,*) '   nkn =',nknouv,'    nel =',nelouv
        write(6,*) '   nlv =',nlvouv,' nlvdim =',nlvdi
        write(6,*)

        ierr = rsouv(
     +                           nunit,nvers
     +                          ,nknouv,nelouv,nlvouv
     +                          ,hlv,hev,ilhv
     +              )

	if( ierr .ne. 0 ) then
		stop 'error stop ouvopen: error reading second header'
	end if

	call timeset(itanf,itend,idtouv)

	return
   99	continue
	write(6,*) 'dimension error :'
	write(6,*) 'nkn,nknouv : ',nkn,nknouv
	write(6,*) 'nel,nelouv : ',nel,nelouv
	write(6,*) 'nlv,nlvouv : ',nlv,nlvouv
	stop 'error stop ouvopen'
	end

c******************************************************

	function ouvnext(it)

	implicit none

	logical ouvnext
	integer it

c	integer nlvdim
c	parameter ( nlvdim = 14 )

	include 'level.h'

	integer nunit,nvers
	common /ouvouv/ nunit,nvers
	integer nknouv,nelouv,nlvouv
	common /ouvnnn/ nknouv,nelouv,nlvouv

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhv(1)
        common /ilhv/ilhv
	real znv(1)
	common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
	real utlnv(nlvdim,1)
	common /utlnv/utlnv
	real vtlnv(nlvdim,1)
	common /vtlnv/vtlnv

	integer ierr
	integer rdouv

	if( nlvdi .ne. nlvdim ) stop 'error stop ouvnext: nlvdim'

	ierr = rdouv(
     +				 nunit,nvers
     +				,it,nlvdim
     +				,nknouv,nelouv,nlvouv
     +				,ilhv
     +				,znv,utlnv,vtlnv
     +			)

	if( ierr .gt. 0 ) then
		stop 'error stop ouvnext: error reading data record'
	else if( ierr .lt. 0 ) then
		ouvnext = .false.
	else
		ouvnext = .true.
	end if

	end

c******************************************************
c******************************************************
c******************************************************

	function rfouv(
     +				 nunit,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtouv
     +				,href,hzmin,hlvmin
     +				,descrp
     +			)

	implicit none

	integer rfouv
	integer nunit,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtouv
	real href,hzmin,hlvmin
	character*80 descrp

	integer idouv,nversm
	parameter( idouv = 26 , nversm = 2 )

	integer id

        read(nunit) id,nvers

	if( id .ne. idouv ) goto 99
	if( nvers .gt. nversm ) goto 98

        read(nunit) nkn,nel,nlv
        read(nunit) itanf,itend,idt,idtouv
        read(nunit) href,hzmin,hlvmin
        read(nunit) descrp

	rfouv = 0

	return
   98	continue
	write(6,*) 'error in version: nvers'
	stop 'error stop rfouv'
   99	continue
	write(6,*) 'error in id: id'
	stop 'error stop rfouv'
	end

c******************************************************

	function rsouv(
     +				 nunit,nvers
     +				,nkn,nel,nlv
     +				,hlv,hev,ilhv
     +			)

	implicit none

	integer rsouv
	integer nunit,nvers
	integer nkn,nel,nlv
	real hlv(1), hev(1)
	integer ilhv(1)

	integer l,ie

        read(nunit) (hlv(l),l=1,nlv)
        read(nunit) (hev(ie),ie=1,nel)
        read(nunit) (ilhv(ie),ie=1,nel)

	rsouv = 0

	end

c******************************************************

	function rdouv(
     +				 nunit,nvers
     +				,it,nlvdim
     +				,nkn,nel,nlv
     +				,ilhv
     +				,znv,utlnv,vtlnv
     +			)

	implicit none

	integer rdouv
	integer nunit,nvers
	integer it,nlvdim
	integer nkn,nel,nlv
	integer ilhv(1)
	real znv(1)
	real utlnv(nlvdim,1)
	real vtlnv(nlvdim,1)

	integer l,ie,k

	rdouv = 0

        read(nunit,end=1) it
        read(nunit) (znv(k),k=1,nkn)
        read(nunit) ((utlnv(l,ie),l=1,ilhv(ie)),ie=1,nel)
        read(nunit) ((vtlnv(l,ie),l=1,ilhv(ie)),ie=1,nel)

	return

    1	continue
	rdouv = -1
	return
	end

c******************************************************

