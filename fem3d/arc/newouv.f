c
c $Id: newouv.f,v 1.6 2000/03/02 09:31:21 georg Exp $
c
c routines to read/write OUV file
c
c revision log :
c
c 01.03.2000    ggu	bug in call to wfouv - hev was not passed
c
c************************************************************

	integer function rfouv(iunit,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtouv
     +				,href,hzmin,hlvmin
     +				,hlv,hev,ilhv
     +				,title)

c reads first record of OUV file

	implicit none

c parameter
	integer maxver,mid
	parameter(maxver=2,mid=26)
c argument
	integer iunit,nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtouv
	real href,hzmin,hlvmin
	real hlv(1),hev(1)
	integer ilhv(1)
	character*80 title
c common
	integer nverso,nknouv,nelouv,nlvouv
	common /ouvcom/ nverso,nknouv,nelouv,nlvouv
c local
	integer id,ios,i,irec
	integer neldim,nlvdim
c save & data
	save /ouvcom/

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
	read(iunit,iostat=ios)	itanf,itend,idt,idtouv
	read(iunit,iostat=ios)	href,hzmin,hlvmin
	if(ios.ne.0) goto 99

	irec=3
	read(iunit,iostat=ios)	title
	if(ios.ne.0) goto 99

	if(nel.gt.neldim.or.nlv.gt.nlvdim) goto 93

	irec=4
	read(iunit,iostat=ios)	(hlv(i),i=1,nlv)
	if(ios.ne.0) goto 99

	irec=5
	read(iunit,iostat=ios)	(hev(i),i=1,nel)
	if(ios.ne.0) goto 99

	irec=6
	read(iunit,iostat=ios)	(ilhv(i),i=1,nel)
	if(ios.ne.0) goto 99

	nverso=nvers
	nknouv=nkn
	nelouv=nel
	nlvouv=nlv
	rfouv=0

	return
   91   continue
        write(6,*) 'Cannot read version ',nvers
        write(6,*) 'Convert to new version with OUVCONV'
        rfouv=91
        return
   93	continue
	write(6,*) 'Dimension error'
	write(6,*) 'nlvdim,nlv :',nlvdim,nlv
	write(6,*) 'neldim,nel :',neldim,nel
	rfouv=93
	return
   95   continue
        write(6,*) 'Old function call ',nvmax
        write(6,*) 'Please adjust call to rfouv and recompile'
        rfouv=95
        return
   97	continue
	write(6,*) 'File is not of OUV file type :',id
	rfouv=97
	return
   98	continue
	write(6,*) 'Unknown version number in OUV file :',nvers
	rfouv=98
	return
   99	continue
	write(6,*) 'Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of OUV file header'
	write(6,*) 'nvers =',nvers
	rfouv=99
	return
	end

c********************************************************************

	integer function wfouv(iunit,nvmax,nvers
     +				,nkn,nel,nlv
     +				,itanf,itend,idt,idtouv
     +				,href,hzmin,hlvmin
     +				,hlv,hev,ilhv
     +				,title)

c writes first record of OUV file

	implicit none

c parameter
	integer maxver,mid
	parameter(maxver=2,mid=26)
c argument
	integer iunit,nvmax,nvers
	integer nkn,nel,nlv
	integer itanf,itend,idt,idtouv
	real href,hzmin,hlvmin
	real hlv(1),hev(1)
	integer ilhv(1)
	character*80 title
c common
	integer nverso,nknouv,nelouv,nlvouv
	common /ouvcom/ nverso,nknouv,nelouv,nlvouv
c local
c	integer id,ios,i,irec
c	integer neldim,nlvdim
	integer i
c save & data
	save /ouvcom/

	if(nvers.eq.0) nvers=maxver

	if(maxver.ne.nvmax) goto 95
	if(maxver.ne.nvers) goto 98

	rewind(iunit)

	write(iunit)		mid,nvers
	write(iunit)		nkn,nel,nlv
	write(iunit)		itanf,itend,idt,idtouv
	write(iunit)		href,hzmin,hlvmin
	write(iunit)		title

	write(iunit)		(hlv(i),i=1,nlv)
	write(iunit)		(hev(i),i=1,nel)
	write(iunit)		(ilhv(i),i=1,nel)

	nverso=nvers
	nknouv=nkn
	nelouv=nel
	nlvouv=nlv
	wfouv=0

	return
   95   continue
        write(6,*) 'wfouv: Old function call ',nvmax
        write(6,*) 'Please adjust call to wfouv and recompile'
        wfouv=95
        return
   98   continue
        write(6,*) 'wfouv: Cannot write version ',nvers
        wfouv=98
        return
	end

c************************************************************

	integer function rdouv(iunit,it,nlvdim,z,u,v,ilhv)

c reads data record of OUV file
c
c EOF -1
c ilhv is input

	implicit none

c argument
	integer iunit,it,nlvdim
	real z(1),u(nlvdim,1),v(nlvdim,1)
	integer ilhv(1)
c common
	integer nverso,nknouv,nelouv,nlvouv
	common /ouvcom/ nverso,nknouv,nelouv,nlvouv
c local
	integer ie,k,l,ios
	integer nkn,nel
c save & data
	save /ouvcom/

c time record

	rdouv = -1
	read(iunit,iostat=ios) it
	if(ios.gt.0) goto 98
	if(ios.lt.0) return	!eof

c data record

	nkn=nknouv
	nel=nelouv

	read(iunit,iostat=ios) (z(k),k=1,nkn)
	read(iunit,iostat=ios) ((u(l,ie),l=1,ilhv(ie)),ie=1,nel)
	read(iunit,iostat=ios) ((v(l,ie),l=1,ilhv(ie)),ie=1,nel)
	if(ios.ne.0) goto 99

	rdouv=0

	return
   98	continue
	write(6,*) 'Error while reading'
	write(6,*) 'time record of OUV file'
	rdouv=98
	return
   99	continue
	write(6,*) 'Error while reading'
	write(6,*) 'data record of OUV file'
	write(6,*) 'it = ',it
	rdouv=99
	return
	end

c************************************************************

	integer function wrouv(iunit,it,nlvdim,z,u,v,ilhv)

c writes data record of OUV file

	implicit none

c argument
	integer iunit,it,nlvdim
	real z(1),u(nlvdim,1),v(nlvdim,1)
	integer ilhv(1)
c common
	integer nverso,nknouv,nelouv,nlvouv
	common /ouvcom/ nverso,nknouv,nelouv,nlvouv
c local
	integer ie,k,l
	integer nkn,nel
c save & data
	save /ouvcom/

	nkn=nknouv
	nel=nelouv

	write(iunit) it
	write(iunit) (z(k),k=1,nkn)
	write(iunit) ((u(l,ie),l=1,ilhv(ie)),ie=1,nel)
	write(iunit) ((v(l,ie),l=1,ilhv(ie)),ie=1,nel)

	wrouv=0

	return
	end

