c
c hash routines
c
c to change for new type of info and key :
c
c probably only changes to make in hashin internal routine
c
c flagf, flagd	must be different from any possible value in key
c ndim          must be greater than about twice the elements
c			to be inserted and must be a prime
c keyt		is key table that should be ok for any use (integer)
c infot		is the information to be inserted  -> change the
c			structure of this array if needed.
c			you must change also routine copyi
c			that actually performs the copy of
c			info and infot. the other functions
c			probably may be left unchanged since
c			they only pass info to hashin
c hash		the hash function -> first statement in hashin
c
c revision log :
c
c 30.06.2000	ggu	error check for hash table half full
c 02.12.2011	ggu	initialize hashin to zero
c
c******************************************************************

c integer function insert(key,info) inserts element into hash table
c integer function retriv(key,info) retrieves element from hash table
c integer function delete(key,info) deletes element from hash table
c integer function visit(key,info)  visits hash table
c subroutine reset                  resets hash table for visiting
c subroutine init                   initializes hash table
c subroutine copyi(infos,infod)     copies from source (s) to destination (d)
c subroutine test                   tests hash routines
c
c arguments :
c
c integer function insert(key,info) in: key,info
c integer function retriv(key,info) in: key  out: info
c integer function delete(key,info) in: key  out: info
c integer function visit(key,info)  out: key,info
c
c return values :
c
c insert	1: inserted  0: already there  -1: error
c retriv	1: found     0: not found      -1: error
c delete	1: deleted   0: not there      -1: error
c visit		1: found     0: finish

c******************************************************************

	integer function insert(key,info)

c inserts element into hash table

	implicit none

	integer key
	integer info
	integer hashin

	insert=hashin(key,info,1)

	return
	end

c******************************************************************

	integer function retriv(key,info)

c retrieves element from hash table

	implicit none

	integer key
	integer info
	integer hashin

	retriv=hashin(key,info,2)

	return
	end

c******************************************************************

	integer function delete(key,info)

c deletes element from hash table

	implicit none

	integer key
	integer info
	integer hashin

	delete=hashin(key,info,3)

	return
	end

c******************************************************************

	integer function visit(key,info)

c visits hash table

	implicit none

	integer key
	integer info
	integer hashin

	visit=hashin(key,info,4)

	return
	end

c******************************************************************

	subroutine reset

c resets hash table for visiting

	implicit none

	integer key
	integer info
	integer hashin
	integer dummy

	dummy=hashin(key,info,5)

	return
	end
c******************************************************************

	subroutine init

c initializes hash table

	implicit none

	integer key
	integer info
	integer hashin
	integer dummy

	dummy=hashin(key,info,0)

	return
	end

c******************************************************************

	subroutine copyi(infos,infod)

c copies from source (infos) to destination (infod)

	implicit none

	integer infos,infod

	infod=infos

	return
	end

c******************************************************************
c	internal routines
c******************************************************************

	integer function hashin(key,info,mode)

c internal hash table function

c mode :
c		0	initialize
c		1	insert
c		2	retrieve
c		3	delete
c		4	visit
c		5	reset for visiting
c
c return :
c		-1	error
c		0	not found, nothing more, already there
c		1	ok
c flags :
c		flagf	flags free entry
c		flagd	flags deleted entry (may be reused)
c
c	-> both flags must be different from
c		any possible key used

	implicit none

	integer ndim,ndim2,flagf,flagd
c using a prime number for ndim is a save choice
c	parameter (ndim=7957)
c	parameter (ndim=21973)
c	parameter (ndim=54563)
	parameter (ndim=99991)
	parameter (ndim2=ndim/2)
	parameter (flagf=0,flagd=-1)

	integer key
	integer mode
	integer info

	integer keyt(ndim)
	integer infot(ndim)

	integer i,ipos,start
	integer itotal

	integer lookup,nfree

	save ipos,keyt,infot,itotal
	data ipos /1/
	data itotal /0/

	hashin=0

	start=mod(key,ndim)+1

	if(mode.eq.1) then 	!-------------------------- insert
	  if( itotal .gt. ndim2 ) goto 95
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.eq.0) then	!not there -> insert
	    i=nfree(keyt,start,flagf,flagd,ndim)
	  else if(i.gt.0) then	!key already there
	    i=0
	  end if
	  if(i.gt.0) then	!insert
	    keyt(i)=key
	    call copyi(info,infot(i))
	    itotal = itotal + 1
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.2) then	!-------------------------- retrieve
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.gt.0) then	!found
	    call copyi(infot(i),info)
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.3) then	!-------------------------- delete
	  i=lookup(key,keyt,start,flagf,ndim)
	  if(i.gt.0) then	!found
	    keyt(i)=flagd
	    itotal = itotal - 1
	    hashin=1
	  else
	    hashin=i
	  end if
	else if(mode.eq.4) then	!-------------------------- visit
	  do while(ipos.le.ndim.and.
     +		(keyt(ipos).eq.flagf.or.keyt(ipos).eq.flagd))
	    ipos=ipos+1
	  end do
	  if(ipos.le.ndim) then
	    key=keyt(ipos)
	    call copyi(infot(ipos),info)
	    ipos=ipos+1
	    hashin=1
	  else
	    hashin=0
	  end if
	else if(mode.eq.5) then	!-------------------------- reset
	  ipos=1
	else if(mode.eq.0) then	!-------------------------- init
	  ipos=1
	  do i=1,ndim
	    keyt(i)=flagf
	  end do
	end if

	return
   95	continue
	write(6,*) 'hash table more than half full... ',itotal
	write(6,*) 'hash table dimension is = ',ndim
	stop 'error stop hashin: dimension of hash table'
	end

c***********************************************************

	integer function lookup(key,keyt,start,flag,ndim)

c looks up index for key in keyt()

c starts search at index start
c flag flags free index
c
c return value :
c			i	index in keyt 
c			0	not found
c			-1	error

	implicit none

	integer key,start
	integer keyt(1)
	integer flag,ndim

	integer icount,inc
	integer i
	integer ndim2

	icount=0
	inc=1
	i=start
	ndim2 = ndim / 2

	lookup=-999
	do while(lookup.eq.-999)
	  if(keyt(i).eq.flag) then
	    lookup=0
	  else if(keyt(i).eq.key) then
	    lookup=i
	  else if(icount.gt.ndim2) then
	    lookup=-1
	  end if
	  i=i+inc
	  if( i .gt. ndim ) i = mod(i,ndim)
	  inc=inc+2
	  icount=icount+1
	end do

	return
	end

c***********************************************************

	integer function nfree(keyt,start,flagf,flagd,ndim)

c looks for first free index in keyt()

c starts search at index start
c flagf, flagd flag free or deleted index
c
c return value :
c			i	index in keyt 
c			-1	error

	implicit none

	integer start
	integer keyt(1)
	integer flagf,flagd,ndim

	integer icount,inc
	integer i
	integer ndim2

	icount=0
	inc=1
	i=start
	ndim2 = ndim / 2

	nfree=-999
	do while(nfree.eq.-999)
	  if(keyt(i).eq.flagf) then
	    nfree=i
	  else if(keyt(i).eq.flagd) then
	    nfree=i
	  else if(icount.gt.ndim2) then
	    nfree=-1
	  end if
	  i=i+inc
	  if( i .gt. ndim ) i = mod(i,ndim)
	  inc=inc+2
	  icount=icount+1
	end do

	return
	end

c*************************************************************

	subroutine test

c tests hash routines

	implicit none

	integer mode
	integer key,info
	integer i

	integer insert,retriv,delete,visit

	call init

	mode=1
	do while(mode.ne.0)
	  write(6,*) '0 exit  1 insert  2 retriv  3 delete  4 visit'
	  read(5,'(i10)') mode

	  if(mode.eq.1) then
		write(6,*) 'key,info :'
		read(5,*) key,info
		i=insert(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.2) then
		write(6,*) 'key :'
		read(5,*) key
		i=retriv(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.3) then
		write(6,*) 'key :'
		read(5,*) key
		i=delete(key,info)
		write(6,*) i,key,info
	  else if(mode.eq.4) then
		call reset
		i=visit(key,info)
		do while(i.gt.0)
		  write(6,*) i,key,info
		  i=visit(key,info)
		end do
	  end if
	end do

	end

c********************************************************************

c	program test
c	call test
c	end
