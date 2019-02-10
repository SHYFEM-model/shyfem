
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c rdgrd routines - read GRD files
c
c contents :
c
c       subroutine rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
c    +                  ,nkndim,neldim,nlidim
c    +                  ,ipv,ipev,ianv,iarv,nen3v,xgv,ygv,hev,hkv)
c				reads grd file
c
c       subroutine rdcom(nco,line,iline,ner,bstop)
c				reads comment from .grd file
c       subroutine rdnode(nkn,nknh,nkndim,line,f,iline,ner,ianz,bstop
c    +                          ,ipv,xgv,ygv,hkv)
c				reads node from .grd file
c       subroutine rdelem(nel,nelh,neldim
c    +                          ,line,f,iline,nin,ner,ianz,bstop
c    +                          ,ipev,iarv,nen3v,hev)
c				reads element from .grd file
c       subroutine rdline(nli,nlidim,line,f,iline,nin,ner,ianz,bstop)
c				reads line from .grd file
c
c       function ifstch(line)
c				finds first char of line that is not blank
c       subroutine fempar(line)
c				read parameters for fem model
c
c revision log :
c
c 20.03.1998    ggu     declare f(6) in subroutines to avoid compiler warnings
c 20.05.1998    ggu     open file with ifileo()
c 22.05.1998    ggu     bug fix (ifileo not declared)
c 27.03.2001    ggu     assign depth -999. to depth not set
c 09.10.2001    ggu     read from stdin
c 09.10.2001    ggu     read node type (ianv)
c 18.10.2005    ggu     error messages slightly changed
c
c**********************************************************

	subroutine rdgrd(
     +			 file
     +			,ner,bstop
     +			,nco,nkn,nknh,nel,nelh,nli
     +			,nkndim,neldim,nlidim
     +			,ipv,ipev
     +			,ianv,iarv
     +                  ,nen3v
     +			,xgv,ygv
     +			,hev,hkv
     +			)

c reads grd file
c
c after read nen3v contains still external node numbers
c use ex2in to convert external to internal node numbers

	implicit none

	character*(*) file	!file name

	integer ner		!error unit (for error and info output)
	logical bstop		!true on return if error

	integer nco		!total number of comments read
	integer nkn		!total number of nodes read
	integer nknh		!total number of nodes read (with depth)
	integer nel		!total number of elements read
	integer nelh		!total number of elements read (with depth)
	integer nli		!total number of lines read

	integer nkndim		!dimension for nodes
	integer neldim		!dimension for elements
	integer nlidim		!dimension for lines

	integer ipv(1)		!external node number
	integer ipev(1)		!external element number

	integer ianv(1) 	!node type
	integer iarv(1)		!element area type

	integer nen3v(3,1)	!element index

	real xgv(1),ygv(1)	!(x,y) coordinates of nodes

	real hev(1)		!depth of elements
	real hkv(1)             !depth of nodes

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

c	integer nin,ner,iline
c	common /grdcom/nin,ner,iline
c	save /grdcom/

	integer nin,ios,iline,ianz,iwhat
	character*132 line
	real f(80)

	integer iscan,ifileo

	bstop = .false.

	nco=0
        nkn=0
	nknh=0
        nel=0
        nelh=0
        nli=0
        iline=0
	nin=1

	xscale=1.
	yscale=1.
	zscale=1.

        if( file .eq. ' ' ) then
          nin = 5
          write(ner,*) 'reading GRD file from stdin...'
        else
	  nin = ifileo(nin,file,'formatted','old')
        end if

	if(nin.le.0) then
	  write(ner,*) 'error opening file'
	  stop 'error stop rdgrd: opening file'
	end if

        read(nin,'(a)',iostat=ios) line

        do while( ios .eq. 0 )
          iline=iline+1
          ianz=iscan(line,1,f)
          iwhat=nint(f(1))
          if( ianz.le.0 .or. iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,line,iline,ner,bstop)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nknh,nkndim
     +				,line,f,iline,ner,ianz,bstop
     +				,ipv,ianv,xgv,ygv,hkv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,nelh,neldim
     +				,line,f,iline,nin,ner,ianz,bstop
     +				,ipev,iarv,nen3v,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nlidim
     +				,line,f,iline,nin,ner,ianz,bstop)
          else
            write(ner,*) 'Type ',iwhat,' at line '
     +                   ,iline,' not recognized'
            bstop=.true.
          end if
          read(nin,'(a)',iostat=ios) line
        end do

        if(ios.lt.0) then
          write(ner,*) iline,' lines read'
        else
          write(ner,*) 'read error close to line ',iline
          stop 'error stop rdgrd: reading line'
        end if

	close(nin)

	return
	end

c**********************************************************

	subroutine rdcom(nco,line,iline,ner,bstop)

c reads comment

	implicit none

	integer nco,iline
	integer ner
	logical bstop
	character*(*) line

	integer i,n
	integer ifstch

	i=ifstch(line)
	n=len(line)

	if(i.gt.0.and.line(i:i).eq.'0') then
	  nco=nco+1
	  if(i.lt.n) then
	    call fempar(line(i+1:))
	  end if
	else if(i.le.0) then
	  !blank line
	else
	  write(ner,*) 'Read error on line ',iline
          bstop=.true.
	end if

	return
	end

c**********************************************************

        subroutine rdnode(nkn,nknh,nkndim,line,f,iline,ner,ianz,bstop
     +				,ipv,ianv,xgv,ygv,hkv)

c reads nodes from .grd file

        implicit none

	integer nkn,nknh,iline,nkndim
	integer ner
	integer ianz
        real f(6)
	logical bstop
	character*(*) line
        integer ipv(1)
        integer ianv(1)         !node type
	real xgv(1),ygv(1)
	real hkv(1)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

        if(ianz.lt.5) then
		write(ner,*) 'Not enough data on line ',iline
		bstop=.true.
		return
        end if

        nkn=nkn+1
	if(nkn.gt.nkndim) goto 99
        ipv(nkn)=nint(f(2))
        ianv(nkn)=nint(f(3))
        xgv(nkn)=f(4)*xscale
        ygv(nkn)=f(5)*yscale
        if(ianz.ge.6) then
            nknh=nknh+1
            hkv(nkn)=f(6)*zscale
	else
            hkv(nkn)=-999.
        end if

        return
   99	continue
	write(6,*) 'dimension of nkndim too low: ',nkndim
	stop 'error stop rdnode: nkndim'
        end

c**********************************************************

        subroutine rdelem(nel,nelh,neldim
     +				,line,f,iline,nin,ner,ianz,bstop
     +				,ipev,iarv,nen3v,hev)

c reads elements from .grd file

        implicit none

	integer nel,nelh,iline,neldim
	integer nin,ner
	integer ianz
        real f(80)
	logical bstop
	character*(*) line
        integer ipev(1),iarv(1)
	integer nen3v(3,1)
	real hev(1)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

	integer i,k,ivert,nvert,ios
        integer iscan

        if(ianz.lt.4) then
		write(ner,*) 'Not enough data on line ',iline
		bstop=.true.
		return
        end if

        nel=nel+1
	if(nel.gt.neldim) goto 99
        ipev(nel)=nint(f(2))
        iarv(nel)=nint(f(3))
        nvert=nint(f(4))

        if(nvert.gt.3) then
		write(ner,*) 'Cannot read element on line ',iline
		bstop=.true.
		return
        end if

        i=4
        ivert=0
	ios=0

        do while( ios .eq. 0 .and. ivert .lt. nvert )    !read verteces
          if( i .ge. ianz ) then
            read(nin,'(a)',iostat=ios) line
	    if(ios.eq.0) then
		iline=iline+1
		i=0
		ianz=iscan(line,1,f)
	    end if
	  else
            i=i+1
            k=nint(f(i))
            ivert=ivert+1
            if(ivert.le.3) nen3v(ivert,nel)=nint(f(i))
	  end if
	end do

	if(ivert.lt.nvert) then
		write(ner,*) 'Could not read all vertices for element ',nel
		bstop=.true.
	end if

	if( i .lt. ianz ) then
		hev(nel)=f(i+1)*zscale
		nelh=nelh+1
	else
        	hev(nel)=-999.
	end if

	return
   99	continue
	write(6,*) 'dimension of neldim too low: ',neldim
	stop 'error stop rdelem: neldim'
	end

c**********************************************************

        subroutine rdline(nli,nlidim,line,f,iline,nin,ner,ianz,bstop)

c reads elements from .grd file

        implicit none

	integer nli,iline,nlidim
	integer nin,ner
	integer ianz
        real f(80)
	logical bstop
	character*(*) line

	integer i,k,ivert,nvert,ios
        integer iscan
	logical bline

	bline = nlidim .gt. 0		!process line data

        if(ianz.lt.4) then
		write(ner,*) 'Not enough data on line ',iline
		bstop=.true.
		return
        end if

        nli=nli+1
        nvert=nint(f(4))

	if( bline ) then
	   if(nli.gt.nlidim) goto 99
c           iplv(nli)=nint(f(2))
c           ialv(nli)=nint(f(3))
	end if

        i=4
        ivert=0
	ios=0

        do while( ios .eq. 0 .and. ivert .lt. nvert )    !read nodes
          if( i .ge. ianz ) then
            read(nin,'(a)',iostat=ios) line
	    if(ios.eq.0) then
		iline=iline+1
		i=0
		ianz=iscan(line,1,f)
	    end if
	  else
            i=i+1
            k=nint(f(i))
            ivert=ivert+1
	    ! here we may process node k of line
c            if( bline ) nen3v(ivert,nel)=nint(f(i))
	  end if
	end do

	if(ivert.lt.nvert) then
		write(ner,*) 'Could not read all nodes for line ',nli
		bstop=.true.
	end if

	return
   99	continue
	write(6,*) 'dimension of nlidim too low: ',nlidim
	stop 'error stop rdline: nlidim'
	end

c******************************************************************************

	function ifstch(line)

c finds first char of line that is not blank or tab

        implicit none

	integer ifstch
	character*(*) line

	integer length,i

	length=len(line)

	do i=1,length
	  if(line(i:i).ne.' ' .and. line(i:i).ne.'	' ) then
		ifstch=i
		return
	  end if
	end do

	ifstch=-1

	return
	end

c******************************************************************************

	subroutine fempar(line)

c read parameters for fem model 
c
c the following special symbols are recognized on a comment line :
c
c (FEM-TITLE)	title for basin			default: first comment line
c (FEM-SCALE)	scale in x/y/z for basin	default: 1.,1.,1.
c (FEM-LATID)	latitude of basin (degrees)	default: 0.
c (FEM-NORTH)	true north of basin (degrees) 	default: 0.
c
c all entries are optional
c
c example (c of fortran comment must be deleted, line starts with 0) :
c
c 0 (FEM-TITLE) test basin
c 0 (FEM-SCALE) 0.5 0.5 2.
c 0 (FEM-LATID) 45.0
c 0 (FEM-NORTH) 90.0
c
        implicit none

	character*(*) line

	integer i,j,n
	integer ifstch,iscan
	character*80 descrr
	logical btitle
	real f(10)
	real grav,fcor,dcor,dirn
	real xscale,yscale,zscale

	common /descrr/ descrr
	common /pkonst/ grav,fcor,dcor,dirn
	common /vscale/ xscale,yscale,zscale

	save btitle
	data btitle /.false./

	i=ifstch(line)
	n=len(line)

	if( i.gt.0 .and. i+10.lt.n ) then
	  if( line(i:i+10) .eq. '(FEM-TITLE)' ) then
		descrr=line(i+11:)
		btitle=.true.
	  else if( line(i:i+10) .eq. '(FEM-SCALE)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.3) then
		  xscale=f(1)
		  yscale=f(2)
		  zscale=f(3)
		else
		  write(6,*) 'error reading (FEM-SCALE) :',j
		  write(6,*) line
		  write(6,*) line(i+11:)
		end if
	  else if( line(i:i+10) .eq. '(FEM-LATID)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.1) then
		  dcor=f(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) line
		end if
	  else if( line(i:i+10) .eq. '(FEM-NORTH)' ) then
		j=iscan(line(i+11:),1,f)
		if(j.eq.1) then
		  dirn=f(1)
		else
		  write(6,*) 'error reading (FEM-NORTH) :'
		  write(6,*) line
		end if
	  end if
	end if

c use first comment as title

	if( i.gt.0 .and. .not.btitle ) then
		descrr=line(i:)
		btitle=.true.
	end if

	return
	end

c******************************************************************************

	subroutine ex2in(nkn,nel,ipv,ipaux,nen3v,ner,bstop)

	implicit none

        integer nkn,nel,ner
        logical bstop
        integer ipv(nkn)
        integer ipaux(nkn)
        integer nen3v(3,nel)

	call isort(nkn,ipv,ipaux)
	call chex2in(nkn,nel,nen3v,ipv,ipaux,ner,bstop)

	end

c*****************************************************************
 
        subroutine chex2in(nkn,nel,nen3v,ipv,index,ner,bstop)
 
c changing extern with intern node numbers in nen3v
 
        implicit none
 
        integer nkn,nel,ner
        logical bstop
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer index(nkn)
 
        integer ie,ii,i,kn
        integer locate
 
        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            i=locate(nkn,ipv,index,kn)
            if(i.le.0) then
              write(ner,*)' warning : node',kn,' not found'
              bstop=.true.
            else
              nen3v(ii,ie)=i
            end if
          end do
        end do
 
        return
        end

c*****************************************************************

