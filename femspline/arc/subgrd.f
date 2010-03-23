c
c $Id: subgrd.f,v 1.5 1998/05/22 12:32:44 georg Exp $
c
c rdgrd routines - read GRD files
c
c contents :
c
c       subroutine rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
c    +                  ,nkndim,neldim,nlidim
c    +                  ,ipv,ipev,iarv,nen3v,xgv,ygv,hev,hkv
c    +                  ,nlndim,iplv,ialrv,ipntlv,inodlv)
c				reads grd file
c       subroutine rdcom(nco,line,iline,ner,bstop)
c				reads comment
c       subroutine rdnode(nkn,nknh,nkndim,line,f,iline,ner,ianz,bstop
c    +                  ,ipv,xgv,ygv,hkv)
c				reads nodes from .grd file
c       subroutine rdelem(nel,nelh,neldim
c    +                  ,line,f,iline,nin,ner,ianz,bstop
c    +                  ,ipev,iarv,nen3v,hev)
c				reads elements from .grd file
c       subroutine rdline(nli,nlidim,line,f,iline,nin,ner,ianz,bstop)
c    +                  ,nlndim,iplv,ialrv,ipntlv,inodlv)
c				reads lines from .grd file
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
c 27.01.1999    ggu     initialize bstop, read also from stdin if no file name
c 27.01.1999    ggu     read line with nodes
c
c**********************************************************

	subroutine rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
     +			,nkndim,neldim,nlidim
     +			,ipv,ipev,iarv,nen3v,xgv,ygv,hev,hkv
     +                  ,nlndim,iplv,ialrv,ipntlv,inodlv)

c reads grd file

c use ex2in to convert external to internal node numbers

	implicit none

	integer ner
	logical bstop
	character*(*) file

	integer nkn,nknh,nel,nelh,nli,nco
	integer nkndim,neldim,nlidim
	integer nen3v(3,1)
	integer ipv(1),ipev(1)
	integer iarv(1)
	real xgv(1),ygv(1)
	real hkv(1),hev(1)

	integer nlndim			!dimension of array containing nodes
	integer iplv(1)			!external number of line
	integer ialrv(1)		!type of line
	integer ipntlv(0:1)		!pointer into inodlv
	integer inodlv(1)		!node numbers of lines (dim. nlndim)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

	integer nin,ios,iline,ianz,iwhat
	character*132 line
	real f(30)

	integer iscan,iround,ifileo

	nco=0
        nkn=0
	nknh=0
        nel=0
        nelh=0
        nli=0
        iline=0
	nin=5

	bstop = .false.

	xscale=1.
	yscale=1.
	zscale=1.

	if( file .ne. ' ' ) then
	  nin = ifileo(55,file,'formatted','old')

	  if(nin.le.0) then
	    write(6,*) 'error opening file'
	    stop 'error stop : rdgrd'
	  end if
	end if

        read(nin,'(a)',iostat=ios) line

        do while( ios .eq. 0 )
          iline=iline+1
          ianz=iscan(line,1,f)
          iwhat=iround(f(1))
          if( ianz.le.0 .or. iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,line,iline,ner,bstop)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nknh,nkndim
     +				,line,f,iline,ner,ianz,bstop
     +				,ipv,xgv,ygv,hkv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,nelh,neldim
     +				,line,f,iline,nin,ner,ianz,bstop
     +				,ipev,iarv,nen3v,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nlidim
     +				,line,f,iline,nin,ner,ianz,bstop
     +                          ,nlndim,iplv,ialrv,ipntlv,inodlv)
          else
            write(6,*) 'Type ',iwhat,' at line '
     +                   ,iline,' not recognized'
            bstop=.true.
          end if
          read(nin,'(a)',iostat=ios) line
        end do

        if(ios.lt.0) then
          write(6,*) iline,' lines read'
        else
          write(6,*) 'read error close to line ',iline
          stop 'error stop : rdgrd'
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
     +				,ipv,xgv,ygv,hkv)

c reads nodes from .grd file

        implicit none

	integer nkn,nknh,iline,nkndim
	integer ner
	integer ianz
        real f(6)
	logical bstop
	character*(*) line
        integer ipv(1)
	real xgv(1),ygv(1)
	real hkv(1)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

        integer iround

        if(ianz.ge.5) then
          nkn=nkn+1
	  if(nkn.gt.nkndim) stop 'error stop : dimension nkn'
          ipv(nkn)=iround(f(2))
          xgv(nkn)=f(4)*xscale
          ygv(nkn)=f(5)*yscale
          if(ianz.ge.6) then
            nknh=nknh+1
            hkv(nkn)=f(6)*zscale
          end if
        else
	  write(ner,*) 'Not enough data on line ',iline
          bstop=.true.
        end if

        return
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
        real f(6)
	logical bstop
	character*(*) line
        integer ipev(1),iarv(1)
	integer nen3v(3,1)
	real hev(1)

	real xscale,yscale,zscale
	common /vscale/ xscale,yscale,zscale

	integer i,k,ivert,nvert,ios
        integer iround,iscan

        if(ianz.lt.4) then
		write(ner,*) 'Not enough data on line ',iline
		bstop=.true.
		return
        end if

        nel=nel+1
	if(nel.gt.neldim) stop 'error stop : dimension nel'
        ipev(nel)=iround(f(2))
        iarv(nel)=iround(f(3))
        nvert=iround(f(4))

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
            k=iround(f(i))
            ivert=ivert+1
            if(ivert.le.3) nen3v(ivert,nel)=iround(f(i))
	  end if
	end do

	if(ivert.lt.nvert) then
		write(ner,*) 'Could not read all verteces for element ',nel
		bstop=.true.
	end if

	if( i .lt. ianz ) then
		hev(nel)=f(i+1)*zscale
		nelh=nelh+1
	end if

	return
	end

c**********************************************************

        subroutine rdline(nli,nlidim,line,f,iline,nin,ner,ianz,bstop
     +                          ,nlndim,iplv,ialrv,ipntlv,inodlv)


c reads elements from .grd file

        implicit none

	integer nli,iline,nlidim
	integer nin,ner
	integer ianz
        real f(6)
	logical bstop
	character*(*) line

	integer nlndim			!dimension of array containing nodes
	integer iplv(1)			!external number of line
	integer ialrv(1)		!type of line
	integer ipntlv(0:1)		!pointer into inodlv
	integer inodlv(1)		!node numbers of lines (dim. nlndim)

	integer i,k,ivert,nvert,ios
	integer ibase,next,ntype
        integer iround,iscan

c nodes in inodlv of line nli run from inodlv(nli-1) + 1 to inodlv(nli)
c
c number of nodes in line nli :    inodlv(nli) - inodlv(nli-1)
c first node in line nli      :    inodlv(nli-1) + 1
c last node in line nli       :    inodlv(nli)

        if(ianz.lt.4) then
		write(ner,*) 'Not enough data on line ',iline
		bstop=.true.
		return
        end if

        nli=nli+1

	ipntlv(0) = 0			!just to be sure...
	ibase = -1			!flag -> do not save nodes

        next=iround(f(2))
        ntype=iround(f(3))
        nvert=iround(f(4))

c	if nlidim is > 0 then we save read data

	if( nlidim .gt. 0 ) then
	  if( nli .gt. nlidim ) then
	    stop 'error stop rdline : dimension nlidim'
	  end if
	  iplv(nli) = next
	  ialrv(nli) = ntype
	  ipntlv(nli) = ipntlv(nli-1) + nvert
	  if( nlndim .gt. 0 ) then	!we only save nodes if nlndim > 0
	    ibase = ipntlv(nli-1)
	    if( ibase + nvert .gt. nlndim ) then
	      stop 'error stop rdline : dimension nlndim'
	    end if
	  end if
	end if

        i=4
        ivert=0
	ios=0

        do while( ios .eq. 0 .and. ivert .lt. nvert )	!read nodes
          if( i .ge. ianz ) then			!must read new line
            read(nin,'(a)',iostat=ios) line
	    if(ios.eq.0) then
		iline=iline+1
		i=0
		ianz=iscan(line,1,f)
	    end if
	  else						!more on line
            i=i+1
            k=iround(f(i))
            ivert=ivert+1
	    if( ibase .ge. 0 ) then
	      inodlv(ibase+ivert) = k
	    end if
	  end if
	end do

	if(ivert.lt.nvert) then
		write(ner,*) 'Could not read all nodes for line ',nli
		bstop=.true.
	end if

	return
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

	integer i,n
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
		i=iscan(line(i+11:),1,f)
		if(i.eq.3) then
		  xscale=f(1)
		  yscale=f(2)
		  zscale=f(3)
		else
		  write(6,*) 'error reading (FEM-SCALE) :'
		  write(6,*) line
		end if
	  else if( line(i:i+10) .eq. '(FEM-LATID)' ) then
		i=iscan(line(i+11:),1,f)
		if(i.eq.1) then
		  dcor=f(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) line
		end if
	  else if( line(i:i+10) .eq. '(FEM-NORTH)' ) then
		i=iscan(line(i+11:),1,f)
		if(i.eq.1) then
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

	
