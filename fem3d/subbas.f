c
c $Id: subbas.f,v 1.4 2009-04-03 16:38:23 georg Exp $
c
c initialization routines
c
c contents :
c
c subroutine sp13rr(nb,nkndim,neldim)		unformatted read from lagoon
c subroutine sp13uw(nb)				unformatted write to lagoon
c
c revision log :
c
c 31.05.1997	ggu	unnecessary routines deleted
c 27.06.1997	ggu	bas routines into own file
c 02.04.2009	ggu	error messages changed
c 12.01.2011	ggu	debug routine introduced (sp13ts)
c
c***********************************************************

	subroutine sp13rr(nb,nkndim,neldim)

c unformatted read from lagoon file
c
c iunit		unit number of file to be read

	implicit none

	integer nb,nkndim,neldim

	character*80 descrr
	common /descrr/descrr
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real  grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	real hm3v(3,1)
	common /hm3v/hm3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	integer nen3v(3,1), iarv(1)
	common /nen3v/nen3v, /iarv/iarv
	integer ipv(1), ipev(1)
	common /ipv/ipv, /ipev/ipev

	integer i,ii,nvers
	integer nversa
	save nversa
c -------------------------
	data nversa / 3 /
c -------------------------

	if(nb.le.0) goto 99

	rewind(nb)

	read(nb) nvers

	if(nvers.lt.nversa) goto 98

	read(nb) nkn,nel,ngr,mbw
	read(nb) dcor,dirn
	read(nb) descrr

	if(nkn.gt.nkndim.or.nel.gt.neldim) goto 97

	read(nb)((nen3v(ii,i),ii=1,3),i=1,nel)
	read(nb)(ipv(i),i=1,nkn)
	read(nb)(ipev(i),i=1,nel)
	read(nb)(iarv(i),i=1,nel)

	read(nb)(xgv(i),i=1,nkn)
	read(nb)(ygv(i),i=1,nkn)
	read(nb)((hm3v(ii,i),ii=1,3),i=1,nel)

c	call sp13ts(nvers,79,0)

	return
   99	continue
	write(6,*) 'Reading basin...'
	write(6,*) 'Cannot read bas file on unit :',nb
	stop 'error stop : sp13ur'
   98	continue
	write(6,*) 'Reading basin...'
	write(6,*) 'Cannot read version < ',nversa
	write(6,*) 'nvers = ',nvers
	stop 'error stop : sp13rr'
   97	continue
	write(6,*) 'Reading basin...'
	write(6,*) 'Dimension error'
	write(6,*) 'nkndim,neldim :',nkndim,neldim
	write(6,*) 'nkn,nel       :',nkn,nel
	stop 'error stop : sp13rr'
	end

c***********************************************************

	subroutine sp13uw(nb)

c unformatted write to lagoon file
c
c nb		unit number for write

	implicit none

	integer nb

	character*80 descrr
	common /descrr/descrr
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real  grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	real hm3v(3,1)
	common /hm3v/hm3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	integer nen3v(3,1), iarv(1)
	common /nen3v/nen3v, /iarv/iarv
	integer ipv(1), ipev(1)
	common /ipv/ipv, /ipev/ipev

	integer i,ii
	integer nvers
	save nvers
c -------------------------
	data nvers / 3 /
c -------------------------

	if(nb.le.0) goto 99

	rewind(nb)

	write(nb) nvers
	write(nb) nkn,nel,ngr,mbw
	write(nb) dcor,dirn
	write(nb) descrr

	write(nb)((nen3v(ii,i),ii=1,3),i=1,nel)
	write(nb)(ipv(i),i=1,nkn)
	write(nb)(ipev(i),i=1,nel)
	write(nb)(iarv(i),i=1,nel)

	write(nb)(xgv(i),i=1,nkn)
	write(nb)(ygv(i),i=1,nkn)
	write(nb)((hm3v(ii,i),ii=1,3),i=1,nel)

c	call sp13ts(nvers,78,0)

	return
   99	continue
	write(6,*) 'Writing basin...'
	write(6,*) 'Cannot write bas file on unit :',nb
	stop 'error stop : sp13uw'
	end

c*************************************************

	subroutine sp13ts(nvers,nb,n)

c test write to unit nb

c writes first n values, if n=0 -> all values

	implicit none

	integer nvers,nb,n

	character*80 descrr
	common /descrr/descrr
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real  grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	real hm3v(3,1)
	common /hm3v/hm3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	integer nen3v(3,1), iarv(1)
	common /nen3v/nen3v, /iarv/iarv
	integer ipv(1), ipev(1)
	common /ipv/ipv, /ipev/ipev

	integer i,ii
	integer nkn1,nel1

	nkn1 = min(nkn,n)
	if( nkn1 .le. 0 ) nkn1 = nkn
	nel1 = min(nel,n)
	if( nel1 .le. 0 ) nel1 = nel

	rewind(nb)

	write(nb,*) 'sp13ts:'
	write(nb,*) nvers
	write(nb,*) nkn,nel,ngr,mbw
	write(nb,*) dcor,dirn
	write(nb,*) descrr

	write(nb,*)((nen3v(ii,i),ii=1,3),i=1,nel1)
	write(nb,*)(ipv(i),i=1,nkn1)
	write(nb,*)(ipev(i),i=1,nel1)
	write(nb,*)(iarv(i),i=1,nel1)

	write(nb,*)(xgv(i),i=1,nkn1)
	write(nb,*)(ygv(i),i=1,nkn1)
	write(nb,*)((hm3v(ii,i),ii=1,3),i=1,nel1)

	return
	end

c*************************************************

