c
c $Id: subapn.f,v 1.13 2010-02-16 16:21:37 georg Exp $
c
c interactive assignment of files with apn
c
c contents :
c
c subroutine assnam(mode)                       assigns run/bas interactivly
c function iapini(mode,nkndim,neldim,matdim)    init routine for ap routines
c subroutine pardef                             reads nls and fnm
c
c revision log :
c
c 22.03.1995	ggu	$$EXINEL - call to exinel changed to ieint
c 21.05.1997	ggu	$$NOBAS - in iapini bug if basin is not opened
c 21.01.1998	ggu	$$APNPAR - only look for apn file in local directory
c 01.05.1998	ggu	$$IMEM - close file only for imem > 0
c 01.05.1998	ggu	interactive questions beautified
c 20.06.1998	ggu	general clean up
c 07.08.1998	ggu	handles negative mode -> do not change names
c 12.02.1999	ggu	useless data structures deleted in iapini
c 12.02.1999	ggu	in assnam do not ask for apn file
c 12.02.1999	ggu	honor auto mode
c 12.02.1999	ggu	read runnam from memfil only if not in apn file
c 05.12.2001	ggu	fixed compiler error with -Wall -pedantic
c 20.06.2003	ggu	in iapini statement shiftet for compiler error Origin
c 15.07.2011	ggu	adjusted, ideffi substituted
c 19.03.2012	ggu	if no basin is given return with "error"
c 27.02.2013	ggu	pass what parameter into nlsa
c 05.09.2013	ggu	read_apn_file() needs integer now
c
c notes :
c
c	iapini
c		->   call pardef
c
c	pardef	
c		->   call assnam
c		->   call nlsina
c		->   call fnminh
c		->   call nlsa
c
c****************************************************************

	subroutine assnam(mode)

c assigns name for run and basin interactivly
c
c mode          1:assign basin 2:assign run 4:assign parameter file
c
c mode < 0	do not change names
c
c combinations are possible, all together : 7

	logical bchang,bbas,brun,bapn,bdo,bask
	character*80 runnam,basnam,name,memfil,text
	character*80 runaux,basaux

	bask = mode .gt. 0

	modeh = abs(mode)

	bbas=.false.
	brun=.false.
	bapn=.false.
	if(modeh.ne.(modeh/2)*2) bbas=.true.
	modeh=modeh/2
	if(modeh.ne.(modeh/2)*2) brun=.true.
	modeh=modeh/2
	if(modeh.ne.(modeh/2)*2) bapn=.true.

	call getfnm('memfil',memfil)

	imem=0
	if(memfil.ne.' ') then
		imem=ifileo(55,memfil,'form','old')
	end if

	if(imem.gt.0) then
		read(imem,'(a)') runaux
		read(imem,'(a)') basaux
	else
		runaux=' '
		basaux=' '
	end if

c if we already have the names through the title section -> use these
c	else use the names just read from memory file

	call getfnm('runnam',runnam)
	call getfnm('basnam',basnam)

	call triml(runnam)
	call triml(basnam)

	if( runnam .eq. ' ' ) runnam = runaux
	if( basnam .eq. ' ' ) basnam = basaux

c ----------------------

	irun=ichanm(runnam)
	ibas=ichanm(basnam)

	bchang=.false.

	if( bask ) then

	  text = 'Enter name of simulation (Default ' 
     +			// runnam(1:irun) 
     +			// ' ) : '
	  irun = 0
	  itxt = ichanm(text)
	  if(brun) irun=igetxt(text(1:itxt),name)
	  if(irun.gt.0) then
	     runnam=name
	     bchang=.true.
	  end if

	  text = 'Enter name of basin      (Default ' 
     +			// basnam(1:ibas) 
     +			// ' ) : '
	  ibas = 0
	  itxt = ichanm(text)
	  if(bbas) ibas=igetxt(text(1:itxt),name)
	  if(ibas.gt.0) then
	     basnam=name
	     bchang=.true.
	  end if

	else

	  write(6,*) 'Name of simulation : ',runnam(1:irun)
	  write(6,*) 'Name of basin      : ',basnam(1:ibas)

	end if

	call putfnm('runnam',runnam)
	call putfnm('basnam',basnam)

	irun=ichanm(runnam)
	ibas=ichanm(basnam)

	if(bchang) then
	   if(imem.le.0) then
		imem=ifileo(55,memfil,'form','new')
		if(imem.le.0) return
	   end if
	   rewind(imem)
	   write(imem,'(a)') runnam(1:irun)
	   write(imem,'(a)') basnam(1:ibas)
	end if

	if( imem .gt. 0 ) close(imem)	!$$IMEM

	end

c*************************************************************

	function iapini(mode,nkndi,neldi,matdim)

c init routine for ap routines
c reads init files nlsina,fnminh,nlsa
c calls pardef (and therefore assnam)
c
c mode          for assnam 1:assign basin 2:assign run 3:both
c nkndim        dimension for nkn for reading bas file
c neldim        dimension for nel for reading bas file
c matdim        (probably useless...)
c iapini        1:success 0:error
c
c mode negative: do not ask for new basin and simulation

	implicit none

	integer mode,nkndi,neldi,matdim

	character*80 descrr
	character*80 line
	character*80 basnew,basold,name
	integer  nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /descrr/descrr
	real f(10)

	integer i,iunit
	integer iapini,ifileo,idefbas

	save basold
	data basold / ' ' /

	logical iseven
	iseven(i) = i .eq. (i/2)*2

	iapini=0

c assign new parameter file

	call pardef(mode)

c if no bas file is opened we are done

	if( iseven( abs(mode) ) ) then		!$$NOBAS
		iapini = 1
		return
	end if

c open bas file

	call getfnm('basnam',basnew)
	if( basnew .eq. ' ' ) return
	if(basnew.ne.basold) then
		iunit=idefbas(basnew,'old')
		if(iunit.le.0) return
		call sp13rr(iunit,nkndi,neldi)
		close(iunit)
		basold=basnew

		write(6,*)
		write(6,*) ' Title of geometrical file :'
		write(6,*)
		write(6,*) descrr
		write(6,*)

		call putpar('dirn',dirn)
	end if

c land boundary		...deleted
c
c       ib=0    chain (middle node)
c       ib!=0   first node of chain
c       ib=1    normal chain, do not fill
c       ib=2    closed chain --> sea
c       ib=3    closed chain --> island

c set up boundary nodes

	if(matdim.gt.0) then
		call sp131k(matdim)
	end if

	iapini=1

	end

c**************************************************************

	subroutine pardef(mode)

c reads default parameters nls and fnm

	implicit none

	integer mode

	integer nin
	integer nmode,iauto
	character*80 file
	character*80 apnnam

	integer ifileo
	real getpar

	logical bfirst
	save bfirst
	data bfirst /.true./

c first call -> set parameters, read parameter file, ...

	if(bfirst) then
	  call nlsina
	  call fnminh

	  call getfnm('apnfil',apnnam)
	  call putfnm('apnnam',apnnam)

	  call read_apn_file(-1)

	  bfirst=.false.
	end if

c get new names for basin and simulation

	nmode = mode
	iauto = nint(getpar('iauto'))
	if( iauto .ne. 0 ) nmode = -abs(mode)

	call assnam(nmode)

	end

c**************************************************************

	subroutine read_apn_file(ivar)

	implicit none

	integer ivar

	integer nin
	integer ifileo

	nin=ifileo(0,'apnstd.str','form','old')
	if( nin .le. 0 ) return

	call nlsa(nin,ivar)
	close(nin)

	end

c**************************************************************

