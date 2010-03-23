c
c $Id: sublst.f,v 1.2 1998/01/22 16:49:32 georg Exp $
c
c*********************************************************************
c
	subroutine setlst(ip,rkey,n,rflag)
c
	implicit none
c
c arguments
	integer ip(1),n
	real rkey(1),rflag
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c local
	integer i
c
	nmax=n
	rlast=rflag
	nins=0
c
	do i=1,n
	  ip(i)=0
	  rkey(i)=rflag
	end do
c
	return
	end
c
c*********************************************************************
c
	subroutine inslst(ip,rkey,ipact,rkact)
c
	implicit none
c
c arguments
	integer ip(1),ipact
	real rkey(1),rkact
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c local
	integer i
c
	if(rkact.le.rlast) return
c
	nins=nins+1
c
	do i=nmax-1,1,-1
		if(rkact.le.rkey(i)) goto 1
		rkey(i+1)=rkey(i)
		ip(i+1)=ip(i)
	end do
    1	continue
	rkey(i+1)=rkact
	ip(i+1)=ipact
c
	rlast=rkey(nmax)
c
	return
	end
c
	
c
c*********************************************************************
c
	subroutine maxlst(n)
c
	implicit none
c
c arguments
	integer n
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c
	n=nins
	if(nins.gt.nmax) n=nmax
c
	return
	end
