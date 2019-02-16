
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

c interactive assignment of files with apn
c
c contents :
c
c subroutine assnam(mode)                       assigns run/bas interactivly
c function iapini(mode,nknddi,nelddi,matdim)    init routine for ap routines
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
c 10.02.2015	ggu	debugged, bcompat gives compatibility with old versions
c 29.04.2015	ggu	generic changes - now works as anticipated
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
c		->   call read_bas_file(nknddi,nelddi)
c		->   call sp131k(matdim)
c
c .memory is only changed with bask == .true.
c .memory is read with bmem == .true., but not changed
c
c****************************************************************

	subroutine assnam(mode)

c assigns name for run and basin (also interactivly)
c
c mode          1:assign basin 2:assign run 4:assign parameter file
c
c mode < 0	do not ask for names
c
c combinations are possible, all together : 7

	implicit none

	integer mode

	logical bchange,bask
	character*80 runnam,basnam
	character*80 runmem,basmem

c---------------------------------------------------------------------
c read memory file
c---------------------------------------------------------------------

	bask = mode .gt. 0
	bchange = .false.

	call read_memory(runmem,basmem)

c---------------------------------------------------------------------
c if we already have the names through the title section -> use these
c	else use the names just read from memory file
c---------------------------------------------------------------------

	call getfnm('runnam',runnam)
	call getfnm('basnam',basnam)

	runnam = adjustl(runnam)
	basnam = adjustl(basnam)

	if( runnam .eq. ' ' ) runnam = runmem
	if( basnam .eq. ' ' ) basnam = basmem

c---------------------------------------------------------------------
c if necessary ask for new values
c---------------------------------------------------------------------

	if( bask ) then
	  stop 'error stop: bask==.true. not supported any more'
	  !call ask_memory(mode,runnam,basnam,bchange)
	end if

c---------------------------------------------------------------------
c see if we have what we need
c---------------------------------------------------------------------

	if( btest(abs(mode),0) ) then
	  if( basnam == ' ' ) then
	    write(6,*) 'no basin given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of basin      : ',trim(basnam)
	end if

	if( btest(abs(mode),1) ) then
	  if( runnam == ' ' ) then
	    write(6,*) 'no simulation given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of simulation : ',trim(runnam)
	end if

c---------------------------------------------------------------------
c writes new information and new values
c---------------------------------------------------------------------

	call putfnm('runnam',runnam)
	call putfnm('basnam',basnam)

	if(bchange) call write_memory(runnam,basnam)

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c*************************************************************

	subroutine ap_set_names(basin,simul)

c sets basin and simulation

	character*(*) basin,simul

	logical haspar

	if( .not. haspar('runnam') ) call addfnm('runnam',' ')
	if( .not. haspar('basnam') ) call addfnm('basnam',' ')

	if( simul .ne. ' ' ) call putfnm('runnam',simul)
	if( basin .ne. ' ' ) call putfnm('basnam',basin)

	end

c*************************************************************

	subroutine ap_get_names(basin,simul)

c sets basin and simulation

	character*(*) basin,simul

	logical haspar

	basin = ' '
	simul = ' '

	if( haspar('runnam') ) call getfnm('runnam',simul)
	if( haspar('basnam') ) call getfnm('basnam',basin)

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine ap_init_basin

c shell for just reading basin

	call ap_init(.false.,1,0,0)

	end

c*************************************************************

	subroutine ap_init(bask,mode,nknddi,nelddi)

c initializes post processing

	implicit none

	logical bask
	integer mode,nknddi,nelddi

	if( bask ) then
	  call iap_init(mode,nknddi,nelddi,0)
	else
	  call iap_init(-mode,nknddi,nelddi,0)
	end if

	end

c*************************************************************

	function iapini(mode,nknddi,nelddi,matdim)

c init routine for ap routines
c
c calls pardef
c calls assnam
c reads basin
c
c iapini        1:success 0:error
c mode          for assnam 1:assign basin 2:assign run 3:both
c nknddi        dimension for nkn for reading bas file
c nelddi        dimension for nel for reading bas file
c matdim        (probably useless...)
c
c mode negative: do not ask for new basin and simulation

	implicit none

	integer iapini
	integer mode,nknddi,nelddi,matdim

	logical bcompat
	integer nmode,iauto

	real getpar

	bcompat = .false.	!set to .true. for compatibility

	call pardef		!we need this for iauto

	nmode = mode
	iauto = nint(getpar('iauto'))
	if( iauto .ne. 0 ) nmode = -abs(mode)

	if( .not. bcompat ) nmode = -abs(mode)	!never ask

	call iap_init(nmode,nknddi,nelddi,matdim)

	iapini = 1

	end

c*************************************************************

	subroutine iap_init(mode,nknddi,nelddi,matdim)

c init routine for ap routines
c
c calls pardef
c calls assnam
c reads basin
c
c iapini        1:success 0:error
c mode          for assnam 1:assign basin 2:assign run 3:both
c nknddi        dimension for nkn for reading bas file
c nelddi        dimension for nel for reading bas file
c matdim        (probably useless...)
c
c mode negative: do not ask for new basin and simulation

	implicit none

	integer mode,nknddi,nelddi,matdim

c---------------------------------------------------------------------
c assign new parameter file
c---------------------------------------------------------------------

	call pardef

c---------------------------------------------------------------------
c get names of basin and simulation
c---------------------------------------------------------------------

	call assnam(mode)

c---------------------------------------------------------------------
c if no bas file has to be opened we are done
c---------------------------------------------------------------------

	if( .not. btest(abs(mode),0) ) return

c---------------------------------------------------------------------
c open bas file
c---------------------------------------------------------------------

	call read_bas_file(nknddi,nelddi)

c---------------------------------------------------------------------
c set up boundary nodes
c---------------------------------------------------------------------

	!if(matdim.gt.0) call sp131k(matdim)

c---------------------------------------------------------------------
c end of routine
c---------------------------------------------------------------------

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine pardef

c reads default parameters nls and fnm

	implicit none

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

	nin=ifileo(-1,'apnstd.str','form','old')
	if( nin .le. 0 ) return

	call nlsa(nin,ivar,.false.)
	close(nin)

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

	simul=' '
	basin=' '

	call getfnm('memfil',memfil)

	if(memfil.eq.' ') return

	imem=ifileo(-1,memfil,'form','old')

	if(imem.le.0) return

	read(imem,'(a)') simul
	read(imem,'(a)') basin

	close(imem)

	end

c**************************************************************

	subroutine write_memory(simul,basin)

	implicit none

	character*(*) simul,basin

	integer imem
	character*80 memfil

	integer ifileo

	call getfnm('memfil',memfil)

	if(memfil.eq.' ') return

	imem=ifileo(0,memfil,'form','new')

	if(imem.le.0) return

	write(imem,'(a)') trim(simul)
	write(imem,'(a)') trim(basin)

	close(imem)

	end

c**************************************************************

	subroutine read_bas_file(nknddi,nelddi)

c opens and reads basin file

	use basin

	implicit none

	integer nknddi,nelddi

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

	call basin_read(iunit)

	close(iunit)

	basold=basnew

	call bas_info

	end

c**************************************************************

