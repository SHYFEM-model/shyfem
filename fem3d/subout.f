
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c utility routines to read/write OUT file - file type 81
c
c reads files up to version 6
c
c contents :
c
c function rfout(iunit,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
c 			reads first record of file OUT
c function wfout(iunit,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
c 			writes first record of file OUT
c function rdout(iunit,nvers,it,nkn,xv)
c			reads data record of file OUT
c function wrout(iunit,nvers,it,nkn,xv)
c			writes data record of file OUT
c
c revision log :
c
c 08.05.1998	ggu	$$id81 - constant 81 replaced with variable id=81 (bug)
c 23.03.2010	ggu	changed v6.1.1
c 16.02.2019	ggu	changed VERS_7_5_60
c
c************************************************************
c
	function rfout(iunit,nvers,itanf,itend,idt,idtout
     +				,href,hzoff,descrp)
c
c reads first record of out file
c
c versions (first record) :
c	1	itanf,itend,idt,idtout
c	2-5	nvers
c	6-...	ntype,nvers
c
	implicit none
c
c arguments
	integer rfout
	integer iunit,nvers,itanf,itend,idt,idtout
	real href,hzoff
	character*80 descrp
c local
	integer ntype,ios
c
	rewind(iunit,iostat=ios)
c
	if(ios.ne.0) then
		write(6,*) 'Cannot rewind file for unit :',iunit
		rfout=71
		return
	end if
c
c first record - find out what version
c
	read(iunit,iostat=ios) itanf,itend,idt,idtout
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
		read(iunit,iostat=ios) itanf,itend,idt,idtout
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
		rfout=22
		return
	else if(ios.lt.0) then
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'first record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=21
		return
	end if
c
c control version number and type of file
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rfout=11
		return
	end if
c
	if(ntype.ne.81) then
		write(6,*) 'rfout : Wrong type of file : ',ntype
		write(6,*) 'Expected : 81'
		rfout=15
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
		read(iunit,iostat=ios)	 itanf,itend,idt,idtout
     +					,href
     +					,hzoff
     +					,descrp
	end if
c
	if(ios.gt.0) then	!error
		write(6,*) 'Error encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=36
		return
	end if
c
	rfout=0
c
	return
	end
c
c********************************************************************
c
	function wfout(iunit,nvers,itanf,itend,idt,idtout
     +				,href,hzoff,descrp)
c
c writes first record of out file
c
c versions (first record) :
c	1	itanf,itend,idt,idtout
c	2-5	nvers
c	6-...	ntype,nvers
c
	implicit none
c
c arguments
	integer wfout
	integer iunit,nvers,itanf,itend,idt,idtout
	integer id
	real href,hzoff
	character*80 descrp
c
	rewind(iunit)
c
	id = 81
	if(nvers.eq.0) nvers=6
c
	if(nvers.eq.1) then
		write(iunit)		 itanf,itend,idt,idtout
	else if(nvers.ge.2.and.nvers.le.5) then
		write(iunit)		 nvers
		write(iunit)		 itanf,itend,idt,idtout
     +					,href
     +					,hzoff
     +					,descrp
	else if(nvers.eq.6) then
		write(iunit)		 id,nvers		!$$id81
		write(iunit)		 itanf,itend,idt,idtout
     +					,href
     +					,hzoff
     +					,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		wfout=11
		return
	end if
c
	wfout=0
c
	return
	end
c
c************************************************************
c
	function rdout(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
c
c reads data record of out file
c
	implicit none
c
c arguments
	integer rdout
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
c local
	integer ios,nk,ne,i
c
c control version number
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rdout=11
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
		rdout=25
		return
	else if(ios.lt.0) then	!eof
		rdout=-1
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
		rdout=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'data record of out file'
		write(6,*) 'it =',it
		write(6,*) 'nkn,nel : ',nkn,nel
		rdout=31
		return
	end if
c
	rdout=0
c
	return
	end
c
c************************************************************
c
	function wrout(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
c
c writes data record of out file
c
	implicit none
c
c arguments
	integer wrout
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
c local
        integer i
c
c control version number
c
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		wrout=11
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
	wrout=0
c
	return
	end
c
