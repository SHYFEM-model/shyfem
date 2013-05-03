c
c $Id: subnsa.f,v 1.21 2009-09-14 08:20:58 georg Exp $
c
c ap utility routines
c
c contents :
c
c subroutine sp158k(imode)              finds boundary nodes
c subroutine nlsa(iunit)		read of parameter file for pp routines
c subroutine rdtita			reads title section of apn files
c
c revision log :
c
c revised on 25.11.88 by ggu (no parameter, merged sp158k and randk)
c revised on 30.11.88 by ggu (no array ig,iamat any more)
c revised on 08.10.90 by ggu (kantv is ordened with direction)
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 12.02.1999	ggu	reading title from own subroutine (with sim and bas)
c 06.12.2004	ggu	new section legvar
c 11.03.2005	ggu	write section title to stdout
c 11.09.2009	ggu	new section $sect
c 27.02.2013	ggu     pass what parameter into nlsa, handle extra info
c
c**********************************************
c
	subroutine sp158k(imode)
c
c finds boundary nodes of not optimized area
c
c only used for imode == 0
c
c imode         0 : called in post-processing routines
c                       ...sets up array kantv
c                       ...needs amat and distroys it
c               1 : called in hp.., used with flood and dry mechanism
c                       ...sets up array kantv
c                       ...sets up suplementary arrays dxv,dyv
c                       ...uses amat and distroys it
c                       ...needs suplementary arrays ipv,iwegv,xgv,ygv
c               2 : called in hp.., used with flood and dry mechanism
c                       ...sets up array kantv
c                       ...sets up suplementary arrays dxv,dyv
c                       ...sets up suplementary array ieltv
c                       ...uses amat and distroys it
c                       ...needs suplementary arrays ipv,ierv,iwegv,xgv,ygv
c                       ...uses bmat and destroys it
c
c array kantv is used also as auxiliary array containing
c ...the number of sides for a node (nkant) and is overwritten in
c ...a second step by the boundary nodes
c boundary nodes in kantv are ordened : first boundary node is in an
c ...anti-clockwise sense ==> going from k1 to kantv(1,k1), the (is)land
c ...is to the right hand side.
c the vector (dxv,dyv) for material boundary nodes is pointing in a
c ...direction that the (is)land is to the right, for an open
c ...boundary the vector is pointing into the system.
c
	logical b0,b1,b2,b12,b99
c
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
c
	common /amat/amat(1)
	common /nen3v/nen3v(3,1), /kantv/kantv(2,1)
c
	common /ipv/ipv(1)
	common /xgv/xgv(1), /ygv/ygv(1)
	common /dxv/dxv(1), /dyv/dyv(1)
	common /iwegv/iwegv(1)
c
	common /bmat/bmat(1)
	common /ieltv/ieltv(3,1), /ierv/ierv(2,1)
c
	integer nkant(1),kant(1),ielmn(1)
c
	equivalence(nkant(1),kantv(1,1))
	equivalence(kant(1),amat(1))
	equivalence(ielmn(1),bmat(1))
c
	ind(i,n) = (ngr*(n-1)+i)
c
	stop 'call to sp158k is broken ... please fix'

	b0=imode.eq.0
	b1=imode.eq.1
	b2=imode.eq.2
	b12=b1.or.b2
	b99=.false.
	iturn=0
c
	do i=1,nkn
	nkant(i)=0
	end do
c
	if(b2) then
		do ii=1,3
		do i=1,nel
		ieltv(ii,i)=0
		end do
		end do
	end if
c
	do ie=1,nel
	if(b12.and.iwegv(ie).ne.0) goto 16
	do ii=1,3
	iii=mod(ii,3)+1
	k1=nen3v(ii,ie)
	k2=nen3v(iii,ie)
c
	do j=1,nkant(k1)
	if(kant(ind(j,k1)).eq.-k2) goto 10
	end do
c
c side not yet found ==> build in
c
	nkant(k1)=nkant(k1)+1
	if( nkant(k1) .gt. ngr ) goto 95
	kant(ind(nkant(k1),k1))=k2
	nkant(k2)=nkant(k2)+1
	if( nkant(k2) .gt. ngr ) goto 95
	kant(ind(nkant(k2),k2))=-k1
c
	if(b2) then
		ielmn(ind(nkant(k1),k1))=ie
		ielmn(ind(nkant(k2),k2))=mod(iii,3)+1
	end if
c
	goto 15
c
c side already found ==> take away
c
   10   continue
	kant(ind(j,k1))=kant(ind(nkant(k1),k1))
	if(b2) then
		iturn=ielmn(ind(j,k1))
		ielmn(ind(j,k1))=ielmn(ind(nkant(k1),k1))
	end if
	nkant(k1)=nkant(k1)-1
c
	do j=1,nkant(k2)
	if(kant(ind(j,k2)).eq.k1) goto 11
	end do
c
	goto 97
c
   11   continue
	kant(ind(j,k2))=kant(ind(nkant(k2),k2))
	if(b2) then
		ieact=ielmn(ind(j,k2))
		ielmn(ind(j,k2))=ielmn(ind(nkant(k2),k2))
	end if
	nkant(k2)=nkant(k2)-1
c
	if(b2) then
		iiii=mod(iii,3)+1
		ieltv(iiii,ie)=ieact
		ieltv(iturn,ieact)=ie
	end if
c
   15   continue
	end do
   16   continue
	end do
c
c copy boundary node index to kantv
c
	if(b0) then
	   do i=nkn,1,-1
	     n=nkant(i)
	     if(n.ne.2) then
		if(n.ne.0) then
		  write(6,*) 'only two sides per node permitted'
		  write(6,*) 'error at node ',ipv(i)
		  write(6,*) 'x,y : ',xgv(i),ygv(i)
		  b99=.true.
		end if
		kantv(1,i)=0
		kantv(2,i)=0
	     else
		if(kant(ind(1,i)).gt.0) then
		  kantv(1,i)=kant(ind(1,i))
		  kantv(2,i)=-kant(ind(2,i))
		else
		  kantv(1,i)=kant(ind(2,i))
		  kantv(2,i)=-kant(ind(1,i))
		end if
	     end if
	   end do
	else
	   do i=nkn,1,-1
	     n=nkant(i)
	     if(n.ne.2) then
		if(n.ne.0) then         !write to error file
		  write(6,*) 'sp158: it,node,sides  ',it,ipv(i),n
		end if
		kantv(1,i)=0
		kantv(2,i)=0
		dxv(i)=0.
		dyv(i)=0.
	     else
		k1=kant(ind(1,i))
		k2=kant(ind(2,i))
		if(k1.gt.0) then
		  kantv(1,i)=k1
		  kantv(2,i)=-k2
		  k2=-k2
		else
		  kantv(1,i)=k2
		  kantv(2,i)=-k1
		  k1=k2
		  k2=kantv(2,i)
		end if
		dxv(i)=xgv(k1)-xgv(k2)
		dyv(i)=ygv(k1)-ygv(k2)
	     end if
	   end do
	end if
c
	if(b99) goto 99
c
	if(b0) return
c
c treatment of open boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do ibc=1,nbc
c
	ibtyp = itybnd(ibc)
	call irbnds(ibc,nkn,n,kant)	!use kant as aux vector
c
	if(ibtyp.le.0) goto 34
c
	do i=1,n
	k=kant(i)
	k1=kantv(1,k)
	k2=kantv(2,k)
	ieact=ierv(1,i)
	iturn=ierv(2,i)
c
	iz=0
	izh1=0
	izh2=0
	do ii=1,n
	if(k1.eq.kant(ii)) izh1=1                !open boundary node
	if(k2.eq.kant(ii)) izh2=1
	end do
c
	kaus = 0

	if(izh1.eq.1) then
		iz=iz+1
	else
		kaus=k1                 !marginary node
	end if
c
	if(izh2.eq.1) then
		iz=iz+1
	else
		kaus=k2
	end if
c
	if(iz.eq.1) then        !marginary open boundary node
		dx=xgv(kaus)-xgv(k)
		dy=ygv(kaus)-ygv(k)
		dxv(k)=dx
		dyv(k)=dy
	else if(iz.eq.2) then   !central open boundary node
		dx=dxv(k)
		dxv(k)=-dyv(k)
		dyv(k)=dx
	else
		dxv(k)=0.
		dyv(k)=0.
c               goto 96
	end if
c
	if(ieact.gt.0.and.b2) then
		ieltv(iturn,ieact)=-1
	end if
   33   continue
	end do
c
   34   continue
	end do
c
	return
   95   continue
        write(6,*) 'nkant too high: ',nkant(k1),nkant(k2),ngr
        stop 'error stop : sp158k'
c   96   continue
c	write(6,*) 'error for input boundary at node ',ipv(k)
c	write(6,*) 'x,y : ',xgv(k),ygv(k)
c	stop 'error stop : sp158k'
   97   continue
	write(6,*) 'missing edge at nodes ',ipv(k1),ipv(k2)
	write(6,*) 'x1,y1 : ',xgv(k1),ygv(k1)
	write(6,*) 'x2,y2 : ',xgv(k2),ygv(k2)
	stop 'error stop : sp158k'
   99   continue
	stop 'error stop : sp158k'
c
	end

c**********************************************

	subroutine sp158kk

c finds boundary nodes of not optimized area

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	real amat(1)
	common /amat/amat
	integer nen3v(3,1), kantv(2,1)
	common /nen3v/nen3v, /kantv/kantv

	integer ipv(1)
	common /ipv/ipv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer nkant(1),kant(1)
	equivalence(nkant(1),kantv(1,1))
	equivalence(kant(1),amat(1))

	integer i,n,k,ie
	integer k1,k2,j,ii,iii
	integer ind

	ind(i,n) = (ngr*(n-1)+i)

	stop 'do not use sp158kk ... please fix'

	do i=1,nkn
	  nkant(i)=0
	end do

c set up nkant and kant

	do ie=1,nel
	  do ii=1,3
	    iii=mod(ii,3)+1
	    k1=nen3v(ii,ie)
	    k2=nen3v(iii,ie)

	    do j=1,nkant(k1)
	      if(kant(ind(j,k1)).eq.-k2) exit
	    end do

	    if( j .gt. nkant(k1) ) then ! side not yet found ==> build in

	      nkant(k1)=nkant(k1)+1
	      if( nkant(k1) .gt. ngr ) goto 95
	      kant(ind(nkant(k1),k1))=k2
	      nkant(k2)=nkant(k2)+1
	      if( nkant(k2) .gt. ngr ) goto 95
	      kant(ind(nkant(k2),k2))=-k1

	    else			! side already found ==> take away

	      kant(ind(j,k1))=kant(ind(nkant(k1),k1))
	      nkant(k1)=nkant(k1)-1

	      do j=1,nkant(k2)
	        if(kant(ind(j,k2)).eq.k1) exit
	      end do

	      if( j .gt. nkant(k2) ) goto 97

	      kant(ind(j,k2))=kant(ind(nkant(k2),k2))
	      nkant(k2)=nkant(k2)-1

	    end if

	  end do
	end do

c copy boundary node index to kantv

	do k=1,nkn
	  n=nkant(i)
	  if(n.ne.2) goto 99
	end do

	do k=1,nkn
	  if(kant(ind(1,i)).gt.0) then
	    kantv(1,i)=kant(ind(1,i))
	    kantv(2,i)=-kant(ind(2,i))
	  else
	    kantv(1,i)=kant(ind(2,i))
	    kantv(2,i)=-kant(ind(1,i))
	  end if
	end do

	return
   95   continue
	write(6,*) 'nkant too high: ',nkant(k1),nkant(k2),ngr
	stop 'error stop : sp158k'
   97   continue
	write(6,*) 'missing edge at nodes ',ipv(k1),ipv(k2)
	write(6,*) 'x1,y1 : ',xgv(k1),ygv(k1)
	write(6,*) 'x2,y2 : ',xgv(k2),ygv(k2)
	stop 'error stop : sp158k'
   99   continue
	write(6,*) 'only two sides per node permitted'
	write(6,*) 'error at node ',ipv(k)
	write(6,*) 'x,y : ',xgv(k),ygv(k)
	stop 'error stop : sp158k'
	end

c******************************************************************
c
	subroutine nlsa(iunit,what)
c
c read of parameter file for post processing routines
c
c iunit		unit number of file

	implicit none

	integer iunit
	character*(*) what

	character*80 name,line,section,extra
	logical bdebug,bread
	integer num
	integer nrdsec,nrdlin,ichanm
	real getpar

	character*80 descrp
	common /descrp/ descrp

	bdebug = .true.
	bdebug = .false.

	if(iunit.le.0) then
c		write(6,*) 'error reading parameter file'
c		write(6,*) 'parameters initialized to default values'
		return
	end if

	call nrdini(iunit)

c loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extra) .eq. 1 )

		if( bdebug ) then
		  write(6,*) 'new section: ',section(1:ichanm(section)),num
		end if

		bread = .false.
		bread = bread .or. what(1:4) .eq. ' '
		bread = bread .or. extra(1:4) .eq. ' '
		bread = bread .or. extra(1:4) .eq. what(1:4)

		if( bread ) then
                  call setsec(section,num)                !remember section
		  write(6,'(a,a)') '$',section(1:60)
		else
		  section = 'skip'
		end if

		if(section.eq.'skip') then
			call nrdskp
		else if(section.eq.'title') then
			call rdtita
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'color') then
			call colrd
		else if(section.eq.'arrow') then
			call nrdins(section)
		else if(section.eq.'legvar') then
			call nrdins(section)
		else if(section.eq.'legend') then
			call legrd
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'sect') then
			call nrdins(section)
		else
			goto 97
		end if
	end do

c end of read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if( bdebug ) call prifnm(6)

	return
   97	continue
	write(6,*) 'Not recognized section key word : ',name
	stop 'error stop : nlsa'
	end

c******************************************************************

        subroutine rdtita

c reads title section of apn files

        implicit none

        character*80 descrp
        common /descrp/ descrp

        character*80 line

        integer nrdlin
	logical bdebug

	bdebug = .true.
	bdebug = .false.

c first line -> title

	if( nrdlin(line) .eq. 0 ) goto 65
	descrp=line
	call putfnm('title',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c maybe more lines ?

	if( nrdlin(line) .eq. 0 ) return	!just one line -> return

c ok, this is simulation

        call putfnm('runnam',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c now basin

        if( nrdlin(line) .eq. 0 ) goto 65
        call putfnm('basnam',line)
	if( bdebug ) write(6,*) 'rdtita: ',line(1:65)

c no more lines allowed

        if( nrdlin(line) .gt. 0 ) goto 65

        return
   65   continue
        write(6,*) 'error in section $title'
        stop 'error stop : rdtitl'
        end

c************************************************************************

