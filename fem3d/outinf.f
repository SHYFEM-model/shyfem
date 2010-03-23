c
c $Id: outinf.f,v 1.11 2009-04-07 10:43:57 georg Exp $

	program outinf

c reads out file and writes info to terminal -> old versions

	implicit none

        include 'param.h'

        character*80 descrr,descrp
        common /descrr/ descrr
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
 
        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

	real xv(3,nkndim)
	real zenv(3,neldim)
	real unv(neldim),vnv(neldim)

	real znv(nkndim)
	real uv(nkndim)
	real vv(nkndim)

        integer nvers,nin
        integer itanf,itend,idt,idtout
	integer it,k
        integer ierr,nread
        integer nkk,nee
	integer date,time
        real href,hzoff
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real uvmin,uvmax
	real vvmin,vvmax
	character*20 line

	integer rdout6,rfout6
	integer iapini,ideffi

	nread=0
c
	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if
c
	nin=ideffi('datdir','runnam','.out','unform','old')
	if(nin.le.0) goto 100
c
        nvers=0
        ierr=rfout6(nin,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
c
        if(ierr.ne.0) goto 100
c
        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' itanf,itend  : ',itanf,itend
        write(6,*) ' idt,idtout   : ',idt,idtout
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nkn,nel
        write(6,*)

	date = 0
	time = 0
	call dtsini(date,time)
	line = ' '
c
  300   continue
c
        nkk=nkn
        nee=nel
        ierr=rdout6(nin,nvers,it,nkk,nee,xv,zenv,unv,vnv)
c
        if(ierr.gt.0) then
		write(6,*) 'error in reading file : ',ierr
		goto 100
        else if(ierr.lt.0) then
		goto 100
        else if(nkn.ne.nkk.or.nel.ne.nee) then
                write(6,*) 'Too less data read'
                write(6,*) 'nkn,nkk :',nkn,nkk
                write(6,*) 'nel,nee :',nel,nee
		goto 100
	end if
c
	nread=nread+1
	call dtsgf(it,line)

	do k=1,nkn
	  uv(k) = xv(1,k)
	  vv(k) = xv(2,k)
	  znv(k) = xv(3,k)
	end do
	call mima(znv,nkn,zmin,zmax)
	call mima(uv,nkn,uvmin,uvmax)
	call mima(vv,nkn,vvmin,vvmax)
	call mima(unv,nel,umin,umax)
	call mima(vnv,nel,vmin,vmax)
c
	write(6,*) 
	write(6,*) 'time : ',it,'    ',line
	write(6,*) 
	write(6,*) ' zmin/zmax  : ',zmin,zmax
	write(6,*) 'uvmin/uvmax : ',uvmin,uvmax
	write(6,*) 'vvmin/vvmax : ',vvmin,vvmax
	write(6,*) 'utmin/utmax : ',umin,umax
	write(6,*) 'vtmin/vtmax : ',vmin,vmax
c
	goto 300
c
  100	continue
c
	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
c
	stop
	end

c*******************************************************************

c******************************************************************
c
c utility routines to read/write out file - file type 81
c
c contents :
c
c integer function rfout6(iunit,nvers,itanf,itend,idt,nout,href,hzoff,descrp)
c 			reads first record of file 8
c integer function wfout6(iunit,nvers,itanf,itend,idt,nout,href,hzoff,descrp)
c 			writes first record of file 8
c integer function rdout6(iunit,nvers,it,nkn,xv)
c			reads data record of file 8
c integer function wrout6(iunit,nvers,it,nkn,xv)
c			writes data record of file 8
c
c************************************************************
c
	integer function rfout6(iunit,nvers,itanf,itend,idt,nout
     +				,href,hzoff,descrp)
c
c reads first record of out file
c
c versions (first record) :
c	1	itanf,itend,idt,nout
c	2-5	nvers
c	6-...	ntype,nvers
c
	implicit none
c
c arguments
	integer iunit,nvers,itanf,itend,idt,nout
	real href,hzoff
	character*80 descrp
c local
	integer ntype,ios
c
	rewind(iunit,iostat=ios)
c
	if(ios.ne.0) then
		write(6,*) 'Cannot rewind file for unit :',iunit
		rfout6=71
		return
	end if
c
c first record - find out what version
c
	read(iunit,iostat=ios) itanf,itend,idt,nout
c
	if(ios.eq.0) then	! no read error ==> first version
		nvers=1
	else if(ios.lt.0) then	!eof
		nvers=itanf
        else if(itanf.eq.81) then !new version
          nvers=itend
        else        !no type available --> up to version 5
          nvers=itanf
	end if

c        write(6,*) nvers,itanf,itend,ios
c
c now rewind and read first record
c
	rewind(iunit,iostat=ios)
c
	if(nvers.eq.1) then
		ntype=81
		read(iunit,iostat=ios) itanf,itend,idt,nout
	else if(nvers.ge.2.and.nvers.le.5) then
		ntype=81
		read(iunit,iostat=ios) nvers
	else if(nvers.ge.6) then
		read(iunit,iostat=ios) ntype,nvers
	end if

c         write(6,*) ntype,nvers,ios
c
	if(ios.gt.0) then
		write(6,*) 'Error encountered while reading'
		write(6,*) 'first record of out file header'
		write(6,*) 'nvers =',nvers
		rfout6=22
		return
	else if(ios.lt.0) then
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'first record of out file header'
		write(6,*) 'nvers =',nvers
		rfout6=21
		return
	end if
c
c control version number and type of file
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rfout6=11
		return
	end if
c
	if(ntype.ne.81) then
		write(6,*) 'rfout6 : Wrong type of file : ',ntype
		write(6,*) 'Expected : 81'
		rfout6=15
		return
	end if
c
c second record
c
	if(nvers.eq.1) then
		hzoff=0.05
		href=0.
		descrp=' '
	else if(nvers.ge.2.and.nvers.le.6) then
		read(iunit,iostat=ios)	 itanf,itend,idt,nout
     +					,href
     +					,hzoff
     +					,descrp
	end if
c
	if(ios.gt.0) then	!error
		write(6,*) 'Error encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout6=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout6=36
		return
	end if
c
	rfout6=0
c
	return
	end
c
c********************************************************************
c
	integer function wfout6(iunit,nvers,itanf,itend,idt,nout
     +				,href,hzoff,descrp)
c
c writes first record of out file
c
c versions (first record) :
c	1	itanf,itend,idt,nout
c	2-5	nvers
c	6-...	ntype,nvers
c
	implicit none
c
c arguments
	integer iunit,nvers,itanf,itend,idt,nout
	real href,hzoff
	character*80 descrp
c
	rewind(iunit)
c
	if(nvers.eq.0) nvers=6
c
	if(nvers.eq.1) then
		write(iunit)		 itanf,itend,idt,nout
	else if(nvers.ge.2.and.nvers.le.5) then
		write(iunit)		 nvers
		write(iunit)		 itanf,itend,idt,nout
     +					,href
     +					,hzoff
     +					,descrp
	else if(nvers.eq.6) then
		write(iunit)		 81,nvers
		write(iunit)		 itanf,itend,idt,nout
     +					,href
     +					,hzoff
     +					,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		wfout6=11
		return
	end if
c
	wfout6=0
c
	return
	end
c
c************************************************************
c
	integer function rdout6(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
c
c reads data record of out file
c
	implicit none
c
c arguments
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
c local
	integer ios,nk,ne,i
c
c control version number
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rdout6=11
		return
	end if
c
c time record
c
	if(nvers.ge.1.and.nvers.le.5) then
		nk=nkn
		ne=nel
		read(iunit,iostat=ios) it
	else if(nvers.eq.6) then
		read(iunit,iostat=ios) it,nk,ne
	end if
c
	if(ios.gt.0) then	!error
		write(6,*) 'Error while reading'
		write(6,*) 'time record of out file'
		rdout6=25
		return
	else if(ios.lt.0) then	!eof
		rdout6=-1
		return
	end if
c
c check dimensions
c
	if(nkn.eq.0) then
c		ok -> skip
	else if(nk.gt.nkn) then	!array too small
		write(6,*) 'Too much data in node record :',nk
		write(6,*) 'Can read only',nkn,' data'
	else
		nkn=nk
	end if
c
	if(nel.eq.0) then
c		ok -> skip
	else if(ne.gt.nel) then	!array too small
		write(6,*) 'Too much data in element record :',ne
		write(6,*) 'Can read only',nel,' data'
	else
		nel=ne
	end if
c
c data record
c
	if(nvers.ge.1.and.nvers.le.5) then
		read(iunit,iostat=ios)	 (xv(i),i=1,3*nkn)
	else if(nvers.eq.6) then
		read(iunit,iostat=ios)	 (xv(i),i=1,3*nkn)
     +					,(zenv(i),i=1,3*nel)
     +					,(u(i),i=1,nel)
     +					,(v(i),i=1,nel)
	end if
c
	if(ios.gt.0) then	!error
		write(6,*) 'Error while reading'
		write(6,*) 'data record of out file'
		write(6,*) 'it      : ',it
		write(6,*) 'nkn,nel : ',nkn,nel
		rdout6=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'data record of out file'
		write(6,*) 'it =',it
		write(6,*) 'nkn,nel : ',nkn,nel
		rdout6=31
		return
	end if
c
	rdout6=0
c
	return
	end
c
c************************************************************
c
	integer function wrout6(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
c
c writes data record of out file
c
	implicit none
c
c arguments
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
c local
        integer i
c
c control version number
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		wrout6=11
		return
	end if
c
	if(nvers.ge.1.and.nvers.le.5) then
		write(iunit)	 it
		write(iunit)	 (xv(i),i=1,3*nkn)
	else if(nvers.eq.6) then
		write(iunit)	 it,nkn,nel
		write(iunit)	 (xv(i),i=1,3*nkn)
     +				,(zenv(i),i=1,3*nel)
     +				,(u(i),i=1,nel)
     +				,(v(i),i=1,nel)
	end if
c
	wrout6=0
c
	return
	end
c
