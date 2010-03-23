c
c $Id: newflx.f,v 1.5 2000/03/02 09:31:20 georg Exp $
c
c flux & volume routines
c
c contents :
c
c real function fluxs(ns)	compute fluxes through section ns
c subroutine fluxa(it)		administers routine for computation of fluxes
c
c notes :
c
c not used anymore
c
c*****************************************************************
c
	real function fluxs(ns)
c
c compute fluxes through section ns
c
c in iflux(1,.) are original (internal) node numbers
c in iflux(2,.) is direction to be used for each element
c preceding each section in iflux(2,.) is length of section
c
	implicit none
c
c arguments
	integer ns
c common
	integer isect,ifluxm,iflux(2,1)
	real ev(13,1)
	integer ieltv(3,1)
	real unv(1),vnv(1)
	common /ifluxc/ isect,ifluxm,iflux
	common /ev/ev, /ieltv/ieltv
	common /unv/unv, /vnv/vnv
c local
	logical bstop
	integer i,ii,iii,n
	integer iin,k,it,nsa
	integer ie,ie1,ie2,ii1,ii2
	real flux,sign
c function
	integer ieext
c save
	integer icall,itold,nsold
	save icall,itold,nsold
	data icall,itold,nsold /0,0,0/

	if(icall.eq.0) then

c         compute number of sections

	  bstop=.false.
	  iin=0
	  isect=0
	  do i=1,ifluxm
	    k=iflux(1,i)
	    if(k.gt.0.and.iin.eq.0) then
	      write(6,*) 'error starting section ',isect+1
	      bstop=.true.
	      isect=isect+1
	      n=0
	      iin=-1
	    else if(k.le.0.and.iin.eq.0) then !new section
	      isect=isect+1
	      n=0
	      iin=i
	    else if(k.gt.0.and.iin.ne.0) then !middle of section
	      n=n+1
	    else if(k.le.0.and.iin.ne.0) then !end of section
	      if(n.eq.0) then
		write(6,*) 'void section ',isect
		bstop=.true.
	      end if
	      if(iin.gt.0) iflux(2,iin)=n !memorize length of section
	      iin=0
	    end if
	  end do

	  if(iin.ne.0) then
	    write(6,*) 'error closing last section ',isect
	    bstop=.true.
	  end if

	  if(bstop) stop 'error stop : fluxs'

c         set up sections

	  it=1
	  do i=1,isect
	    n=iflux(2,it)
	    do ii=it+1,it+n
	      ie=iflux(1,ii)
	      ie1=iabs(iflux(1,ii-1))
	      ii1=0
	      ie2=iabs(iflux(1,ii+1))
	      ii2=0
	      do iii=1,3
		if(ieltv(iii,ie).eq.ie1) ii1=iii
		if(ieltv(iii,ie).eq.ie2) ii2=iii
	      end do
	      if(ii1.eq.0.and.ii2.eq.0) then
		bstop=.true.
		write(6,*) 'no neighbors for element ',ieext(ie)
	      else if(ii1.eq.0) then
		if(ii.ne.it+1) then
		  bstop=.true.
		  write(6,*) 'no neighbor for element ',ieext(ie)
		end if
		ii1=mod(ii2,3)+1
		if(ieltv(ii1,ie).ne.0) ii1=mod(ii1,3)+1  !look for bound.
	      else if(ii2.eq.0) then
		if(ii.ne.it+n) then
		  bstop=.true.
		  write(6,*) 'no neighbor for element ',ieext(ie)
		end if
		ii2=mod(ii1,3)+1
		if(ieltv(ii2,ie).ne.0) ii2=mod(ii2,3)+1  !look for bound.
	      end if
	      iii=6-ii1-ii2
	      if(mod(3+ii1-ii2,3).eq.2) iii=-iii
	      iflux(2,ii)=iii
	    end do
	    it=it+n+2
	  end do

	  if(bstop) stop 'error stop : fluxs'

	  icall=1

	end if

c compute discharge through section ns

	if(ns.gt.isect.or.ns.lt.1) then
	  fluxs=0.
	  return
	else if(ns.eq.nsold) then
	  it=itold
	else
	  it=1
	  nsa=1
	  do while(nsa.lt.ns)
	    nsa=nsa+1
	    it=it+iflux(2,it)+2
	  end do
	end if

	flux=0.
	n=iflux(2,it)
	do i=it+1,it+n
	  ie=iflux(1,i)
	  ii=iflux(2,i)
	  if(ii.gt.0) then
	    sign=1.
	  else
	    ii=-ii
	    sign=-1.
	  end if
	  flux = flux + sign * ev(10,ie) *
     +              ( unv(ie)*ev(3+ii,ie) + vnv(ie)*ev(6+ii,ie) )
	end do

	fluxs = 12. * flux

	itold=it+iflux(2,it)+2
	nsold=ns+1

	return
	end

c*****************************************************************

	subroutine fluxa(it)

c administers routine for computation of fluxes

	implicit none

c arguments
	integer it
c common
	integer isect,ifluxm,iflux(2,1)
	common /ifluxc/ isect,ifluxm,iflux
	integer kfluxm,kflux(1)
	common /kfluxc/ kfluxm,kflux
c local
	character*80 dir,nam,file
	real dummy
	integer i
c function
	real fluxs
c save
	integer icall
	save icall
	data icall /0/

	if(icall.eq.0) then	!first call

c		data from STR file was read into kflux -> copy

		ifluxm = kfluxm
		do i=1,kfluxm
		  iflux(1,i) = kflux(i)
		end do

        	dummy = fluxs(0)	!initialize routine

c		open file

       		call getfnm('datdir',dir)
        	call getfnm('runnam',nam)
        	call mkname(dir,nam,'.flx',file)

        	open(48,file=file,status='unknown',form='formatted')

	end if

	icall = icall + 1

c	write to file

	write(48,*) it,(fluxs(i),i=1,isect)

	return
	end


