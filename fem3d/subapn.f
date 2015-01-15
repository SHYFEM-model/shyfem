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
c 14.01.2015	ggu	reorganized and cleaned
c
c notes :
c
c	call iapini
c		->   call pardef
c			->   call nlsina
c			->   call fnminh
c	  		->   call read_apn_file(-1)
c				->   call nlsa
c		->   call assnam
c		->   call read_bas_file(nkndi,neldi)
c		->   call sp131k(matdim)
c
c****************************************************************

	subroutine assnam(mode)

c assigns name for run and basin (also interactivly)
c
c mode          1:assign basin 2:assign run 4:assign parameter file
c
c mode < 0	do not change names
c
c combinations are possible, all together : 7

	implicit none

	integer mode

	logical bchang,bask
	character*80 runnam,basnam
	character*80 runaux,basaux

c---------------------------------------------------------------------
c read memory file
c---------------------------------------------------------------------

	bask = mode .gt. 0

	call read_memory(runaux,basaux)

c---------------------------------------------------------------------
c if we already have the names through the title section -> use these
c	else use the names just read from memory file
c---------------------------------------------------------------------

	call getfnm('runnam',runnam)
	call getfnm('basnam',basnam)

	call triml(runnam)
	call triml(basnam)

	if( runnam .eq. ' ' ) runnam = runaux
	if( basnam .eq. ' ' ) basnam = basaux

c---------------------------------------------------------------------
c if necessary ask for new values
c---------------------------------------------------------------------

	if( bask ) then
	  call ask_memory(mode,runnam,basnam,bchang)
	end if

c---------------------------------------------------------------------
c see if we have what we need
c---------------------------------------------------------------------

	if( btest(abs(mode),0) ) then
	  if( basnam == ' ' ) then
	    write(6,*) 'no basin given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of basin      : ',basnam(1:len_trim(basnam))
	end if

	if( btest(abs(mode),1) ) then
	  if( runnam == ' ' ) then
	    write(6,*) 'no simulation given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of simulation : ',runnam(1:len_trim(runnam))
	end if

c---------------------------------------------------------------------
c writes new information and new values
c---------------------------------------------------------------------

	call putfnm('runnam',runnam)
	call putfnm('basnam',basnam)

	if(bchang) call write_memory(runnam,basnam)

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c*************************************************************

	function iapini(mode,nkndi,neldi,matdim)

c init routine for ap routines
c
c calls pardef
c calls assnam
c reads basin
c
c iapini        1:success 0:error
c mode          for assnam 1:assign basin 2:assign run 3:both
c nkndim        dimension for nkn for reading bas file
c neldim        dimension for nel for reading bas file
c matdim        (probably useless...)
c
c mode negative: do not ask for new basin and simulation

	implicit none

	integer iapini
	integer mode,nkndi,neldi,matdim

	integer nmode,iauto
	real getpar

	iapini = 1

c---------------------------------------------------------------------
c assign new parameter file
c---------------------------------------------------------------------

	call pardef(mode)

c---------------------------------------------------------------------
c get names of basin and simulation
c---------------------------------------------------------------------

	nmode = mode
	iauto = nint(getpar('iauto'))
	if( iauto .ne. 0 ) nmode = -abs(mode)

	call assnam(nmode)

c---------------------------------------------------------------------
c if no bas file has to be opened we are done
c---------------------------------------------------------------------

	if( .not. btest(abs(mode),0) ) return

c---------------------------------------------------------------------
c open bas file
c---------------------------------------------------------------------

	call read_bas_file(nkndi,neldi)

c---------------------------------------------------------------------
c set up boundary nodes
c---------------------------------------------------------------------

	if(matdim.gt.0) call sp131k(matdim)

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	!call nbasin_transfer

	end

c**************************************************************

	subroutine pardef(mode)

c reads default parameters nls and fnm

	implicit none

	integer mode

	character*80 apnnam

	logical bfirst
	save bfirst
	data bfirst /.true./

	if(bfirst) then
	  call nlsina
	  call fnminh

	  call getfnm('apnfil',apnnam)
	  call putfnm('apnnam',apnnam)

	  call read_apn_file(-1)

	  bfirst=.false.
	end if

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

        subroutine nbasin_transfer
        implicit none
        include 'nbasin.h'
        call bas_get_para(nkn,nel,ngr,mbw)
        end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine read_memory(simul,basin)

	implicit none

	character*(*) simul,basin

	integer imem
	character*80 memfil

	integer ifileo

	call getfnm('memfil',memfil)

	if(memfil.ne.' ') imem=ifileo(0,memfil,'form','old')

	simul=' '
	basin=' '

	if(imem.gt.0) then
		read(imem,'(a)') simul
		read(imem,'(a)') basin
		close(imem)
	end if

	end

c**************************************************************

	subroutine write_memory(simul,basin)

	implicit none

	character*(*) simul,basin

	integer imem
	character*80 memfil

	integer ifileo

	call getfnm('memfil',memfil)

	if(memfil.ne.' ') imem=ifileo(0,memfil,'form','new')

	if(imem.gt.0) then
		read(imem,'(a)') simul(1:len_trim(simul))
		read(imem,'(a)') basin(1:len_trim(basin))
		close(imem)
	end if

	end

c**************************************************************

	subroutine ask_memory(mode,simul,basin,bchang)

	implicit none

	integer mode
	character*(*) simul,basin
	logical bchang

	integer irun,ibas
	character*80 text,name

	integer igetxt

	bchang=.false.

	if( btest(mode,1) ) then
	  text = 'Enter name of simulation (Default ' 
     +			// simul(1:len_trim(simul)) 
     +			// ' ) : '
	  irun=igetxt(text(1:len_trim(text)),name)
	  if(irun.gt.0) then
	     simul=name
	     bchang=.true.
	  end if
	end if

	if( btest(mode,0) ) then
	  text = 'Enter name of basin      (Default ' 
     +			// basin(1:len_trim(basin)) 
     +			// ' ) : '
	  ibas=igetxt(text(1:len_trim(text)),name)
	  if(ibas.gt.0) then
	     basin=name
	     bchang=.true.
	  end if
	end if

	end

c**************************************************************

	subroutine read_bas_file(nkndi,neldi)

c opens and reads basin file

	implicit none

	integer nkndi,neldi

	integer iunit
	character*80 basnew

	integer idefbas

	character*80 basold
	save basold
	data basold / ' ' /

	call getfnm('basnam',basnew)

	if( basnew .eq. ' ' ) return
	if( basnew .eq. basold ) return

	iunit=idefbas(basnew,'old')
	if(iunit.le.0) then
	  stop 'error stop read_bas_file: error reading basin'
	end if
	call sp13rr(iunit,nkndi,neldi)
	close(iunit)

	basold=basnew

	call bas_info

	end

c**************************************************************

