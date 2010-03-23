c
c $Id: newous.f,v 1.3 1998/01/22 17:13:43 georg Exp $
c
c************************************************************

	integer function rfous(iunit,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtous
     +				,href,hzmin,hlvmin
     +				,hlv,ilhv
     +				,title)

c reads first record of OUS file

	implicit none

c parameter
	integer maxver,mid
	parameter(maxver=2,mid=25)
c argument
	integer iunit,nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtous
	real href,hzmin,hlvmin
	real hlv(1)
	integer ilhv(1)
	character*80 title
c common
	integer nverso,nknous,nelous,nlvous
	common /ouscom/ nverso,nknous,nelous,nlvous
c local
	integer id,ios,i,irec
	integer neldim,nlvdim
c save & data
	save /ouscom/

	if(maxver.ne.nvmax) goto 95

	neldim = nel
	nlvdim = nlv

	rewind(iunit)

	irec=1
	read(iunit,iostat=ios) id,nvers
	if(ios.ne.0) goto 99

	if(id.ne.mid) goto 97
	if(nvers.le.0.or.nvers.gt.maxver) goto 98
	if(nvers.ne.maxver) goto 91

	irec=2
	read(iunit,iostat=ios)	nkn,nel,nlv
	read(iunit,iostat=ios)	itanf,itend,idt,idtous
	read(iunit,iostat=ios)	href,hzmin,hlvmin
	if(ios.ne.0) goto 99

	if(nel.gt.neldim.or.nlv.gt.nlvdim) goto 93

	irec=3
	read(iunit,iostat=ios)	(hlv(i),i=1,nlv)
	if(ios.ne.0) goto 99

	irec=4
	read(iunit,iostat=ios)	(ilhv(i),i=1,nel)
	if(ios.ne.0) goto 99

	irec=5
	read(iunit,iostat=ios)	title
	if(ios.ne.0) goto 99

	nverso=nvers
	nknous=nkn
	nelous=nel
	nlvous=nlv
	rfous=0

	return
   91   continue
        write(6,*) 'Cannot read version ',nvers
        write(6,*) 'Convert to new version with OUSCONV'
        rfous=91
        return
   93	continue
	write(6,*) 'Dimension error'
	write(6,*) 'nlvdim,nlv :',nlvdim,nlv
	write(6,*) 'neldim,nel :',neldim,nel
	rfous=93
	return
   95   continue
        write(6,*) 'Old function call ',nvmax
        write(6,*) 'Please adjust call to rfous and recompile'
        rfous=95
        return
   97	continue
	write(6,*) 'File is not of OUS file type :',id
	rfous=97
	return
   98	continue
	write(6,*) 'Unknown version number in OUS file :',nvers
	rfous=98
	return
   99	continue
	write(6,*) 'Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of OUS file header'
	write(6,*) 'nvers =',nvers
	rfous=99
	return
	end

c********************************************************************

	integer function wfous(iunit,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtous
     +				,href,hzmin,hlvmin
     +				,hlv,ilhv
     +				,title)

c writes first record of OUS file

	implicit none

c parameter
	integer maxver,mid
	parameter(maxver=2,mid=25)
c argument
	integer iunit,nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtous
	real href,hzmin,hlvmin
	real hlv(1)
	integer ilhv(1)
	character*80 title
c common
	integer nverso,nknous,nelous,nlvous
	common /ouscom/ nverso,nknous,nelous,nlvous
c local
c	integer id,ios,i,irec
	integer i
c	integer neldim,nlvdim
c save & data
	save /ouscom/

	if(nvers.eq.0) nvers=maxver

	if(maxver.ne.nvmax) goto 95
	if(maxver.ne.nvers) goto 98

	rewind(iunit)

	write(iunit)		mid,nvers
	write(iunit)		nkn,nel,nlv
	write(iunit)		itanf,itend,idt,idtous
	write(iunit)		href,hzmin,hlvmin
	write(iunit)		(hlv(i),i=1,nlv)
	write(iunit)		(ilhv(i),i=1,nel)
	write(iunit)		title

	nverso=nvers
	nknous=nkn
	nelous=nel
	nlvous=nlv
	wfous=0

	return
   95   continue
        write(6,*) 'wfous: Old function call ',nvmax
        write(6,*) 'Please adjust call to wfous and recompile'
        wfous=95
        return
   98   continue
        write(6,*) 'wfous: Cannot write version ',nvers
        wfous=98
        return
	end

c************************************************************

	integer function rdous(iunit,it,nlvdim,z,u,v,ilhv)

c reads data record of OUS file
c
c EOF -1
c ilhv is input

	implicit none

c argument
	integer iunit,it,nlvdim
	real z(1),u(nlvdim,1),v(nlvdim,1)
	integer ilhv(1)
c common
	integer nverso,nknous,nelous,nlvous
	common /ouscom/ nverso,nknous,nelous,nlvous
c local
	integer ie,k,l,ios
	integer nkn,nel
c save & data
	save /ouscom/

c time record

	rdous = -1
	read(iunit,iostat=ios) it
	if(ios.gt.0) goto 98
	if(ios.lt.0) return	!eof

c data record

	nkn=nknous
	nel=nelous

	read(iunit,iostat=ios) (z(k),k=1,nkn)
	read(iunit,iostat=ios) ((u(l,ie),l=1,ilhv(ie)),ie=1,nel)
	read(iunit,iostat=ios) ((v(l,ie),l=1,ilhv(ie)),ie=1,nel)
	if(ios.ne.0) goto 99

	rdous=0

	return
   98	continue
	write(6,*) 'Error while reading'
	write(6,*) 'time record of OUS file'
	rdous=98
	return
   99	continue
	write(6,*) 'Error while reading'
	write(6,*) 'data record of OUS file'
	write(6,*) 'it = ',it
	rdous=99
	return
	end

c************************************************************

	integer function wrous(iunit,it,nlvdim,z,u,v,ilhv)

c writes data record of OUS file

	implicit none

c argument
	integer iunit,it,nlvdim
	real z(1),u(nlvdim,1),v(nlvdim,1)
	integer ilhv(1)
c common
	integer nverso,nknous,nelous,nlvous
	common /ouscom/ nverso,nknous,nelous,nlvous
c local
	integer ie,k,l
	integer nkn,nel
c save & data
	save /ouscom/

	nkn=nknous
	nel=nelous

	write(iunit) it
	write(iunit) (z(k),k=1,nkn)
	write(iunit) ((u(l,ie),l=1,ilhv(ie)),ie=1,nel)
	write(iunit) ((v(l,ie),l=1,ilhv(ie)),ie=1,nel)

	wrous=0

	return
	end

