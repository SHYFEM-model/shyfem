
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-1998,2000-2002,2007,2010-2012  Georg Umgiesser
!    Copyright (C) 2014,2017-2019  Georg Umgiesser
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

c create default names
c
c contents :
c
c subroutine mkname(dir,name,ext,file)            makes file name
c subroutine defmak(defdir,defnam,ext,file)	  makes file with defaults
c function ideffi(defdir,defnam,ext,form,status)  opens file in default dir
c function ifemop(ext,form,status)		  opens file with default name
c function ifemopa(text,ext,form,status)	  opens file with default name
c function ifem_open_file(ext,status)		  opens unformated file
c function ifem_test_file(ext,status)		  tries to open unformated file
c
c revision log :
c
c 23.05.1997	ggu	$$EXTENS - default extension may be overwritten
c 18.06.1997	ggu	restructured - idefna,idefop,idefts,idefun deleted
c 16.01.1998	ggu	idefop reintroduced -> to avoid link error
c 21.01.1998	ggu	in mkname: give extension with or without dot
c 08.08.2000	ggu	new routine ifemop
c 27.11.2001	ggu	error message rewritten
c 11.10.2002	ggu	new subroutine deffile
c 07.03.2007	ggu	new routine ifem_open_file
c 23.03.2010	ggu	changed v6.1.1
c 29.04.2010	ggu	new routine ifem_open_file
c 03.05.2010	ggu	new routine ifem_choose_file() and add_extension()
c 02.07.2011	ggu	idefna,idefop finally deleted
c 13.07.2011	ggu	cleaned from old structures
c 18.08.2011	ggu	bug fix in idefbas -> use status passed in
c 01.06.2012	ggu	changed VERS_6_1_53
c 18.06.2014	ggu	changed VERS_6_1_77
c 09.05.2017	ggu	add_extension -> to subst_extension, new add_extension
c 05.12.2017	ggu	changed VERS_7_5_39
c 07.12.2017	ggu	changed VERS_7_5_40
c 05.10.2018	ggu	avoid run time error in subst_extension()
c 16.10.2018	ggu	changed VERS_7_5_50
c 16.02.2019	ggu	changed VERS_7_5_60
c 03.05.2019	ggu	new routines to get extension and name
c 21.05.2019	ggu	changed VERS_7_5_62
c
c notes :
c
c exchange ideffi with ifemop (many apperances)
c eliminate call to getpar/getfnm
c
c**************************************************************

        subroutine mkname(dir,name,ext,file)

c makes file name given its constituents
c
c dir   directory
c name  name
c ext   extension (with or without dot)
c file  created file name (return)

        implicit none

c arguments
        character*(*) dir,name,ext,file
c local
        integer nall,nstart,nend,naux
c function
        integer ichafs,ichanm

        nall=1
        file=' '

        nstart=ichafs(dir)
        nend=ichanm(dir)
        if(nend.gt.0) then
		file(nall:)=dir(nstart:nend)
        	nall=nall+nend-nstart+1
	end if

        nstart=ichafs(name)
        nend=ichanm(name)
        if(nend.gt.0) then
		file(nall:)=name(nstart:nend)
       		nall=nall+nend-nstart+1
	end if

	call subst_extension(file,ext,.false.)

	return

	end

c**************************************************************
c**************************************************************
c**************************************************************

	module default_names

        character*80, save :: def_bas = ' '
        character*80, save :: def_nam = ' '

	end module default_names

c**************************************************************

	subroutine set_default_names(defbas,defnam)

c sets default names to be used with def routines
c
c must be called before any other routine in subdef can be used

	use default_names

	implicit none

        character*(*) defbas,defnam
	
	def_bas = defbas
	def_nam = defnam

	end

c**************************************************************

	subroutine get_default_names(defbas,defnam)

c gets default names to be used with def routines

	use default_names

        implicit none

        character*(*) defbas,defnam
	
	defbas = def_bas
	defnam = def_nam

	end
	
c**************************************************************
c**************************************************************
c**************************************************************

        subroutine def_make(ext,file)

c makes file with defaults supplied
c
c ext   extension (with dot)
c file  created file name (return)

	use default_names

        implicit none

        character*(*) ext,file

        character*80 dir,name

	name = def_nam
	dir = ' '
        !call getfnm('datdir',dir)	! this has to be deleted
        call getfnm('runnam',name)	! this has to be deleted

	call mkname(dir,name,ext,file)

	end

c**************************************************************

        subroutine defmak(defdir,defnam,ext,file)

c FIXME -> take defdir,defnam from common block

c makes file with defaults supplied
c
c defdir   directory
c defnam  name
c ext   extension (with dot)
c file  created file name (return)

        implicit none

        character*(*) defdir,defnam,ext,file
	character*80 dir,name

	dir = ' '
        !if( defdir .ne. ' ' ) call getfnm(defdir,dir)
        call getfnm(defnam,name)

	call mkname(dir,name,ext,file)

	end

c**************************************************************

        function ideffi(defdir,defnam,ext,form,status)

c FIXME -> substitute with ifemop (many appearances)

c opens file in default dir

c defdir   directory
c defnam  name
c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ideffi
        character*(*) defdir,defnam,ext,status,form
	character*80 file
	integer ifileo

	call defmak(defdir,defnam,ext,file)
        call def_make(ext,file)
	ideffi=ifileo(0,file,form,status)

	end

c**************************************************************

	function idefbas(basnam,status)

	implicit none

	integer idefbas
	character*(*) basnam
	character*(*) status

	character*80 name
	integer ifileo

        name = basnam
        call add_extension(name,'.bas')

        idefbas=ifileo(0,name,'unform',status)

	end

c**************************************************************
c**************************************************************
c**************************************************************
c opening of default simulation
c**************************************************************
c**************************************************************
c**************************************************************

        function ifemop(ext,form,status)

c opens file with default name (run) and extension given for fem model
c returns with error code

c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ifemop
        character*(*) ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

        call def_make(ext,file)
	ifemop=ifileo(0,file,form,status)

	end

c**************************************************************

        function ifemopa(text,ext,form,status)

c opens file with default name (run) and extension given for fem model
c in case of error exits with error message text

c text  error message
c ext   extension (with dot)
c form  formatted ?
c status open status

        implicit none

	integer ifemopa
        character*(*) text,ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

        call def_make(ext,file)
	ifemopa=ifileo(0,file,form,status)

	if( ifemopa .le. 0 ) then
	  write(6,*) 'error opening file ',file
	  write(6,*) text
	  stop 'error stop ifemopa'
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************
c unformatted opening of default simulation
c**************************************************************
c**************************************************************
c**************************************************************

        function ifem_open_file(ext,status)

c opens unformated file with default name (run) and extension given
c in case of error exits 

c ext		extension (with dot)
c status	open status

        implicit none

	integer ifem_open_file
        character*(*) ext,status

	integer ifem_test_file

	ifem_open_file = ifem_test_file(ext,status)

	if( ifem_open_file .le. 0 ) then
	  stop 'error stop ifem_open_file'
	end if

	end

c**************************************************************

        function ifem_test_file(ext,status)

c tries to open unformated file with default name (run) and extension given

c ext		extension (with dot)
c status	open status

        implicit none

	integer ifem_test_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form
	integer ifileo

	form = 'unform'
        call def_make(ext,file)

	ifem_test_file = ifileo(0,file,form,status)

	if( ifem_test_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

c**************************************************************

        function ifem_choose_file(ext,status)

c tries to open unformated file with default name (run) and extension given
c insists on extension -> if name has extension substitute it with ext

c ext		extension (with dot)
c status	open status

        implicit none

	integer ifem_choose_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form
	integer ifileo

	form = 'unform'
        call def_make(ext,file)
	call subst_extension(file,ext,.true.)

	ifem_choose_file = ifileo(0,file,form,status)

	if( ifem_choose_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine add_extension(name,ext)

c adds extension to file name if not already there
c
c ext should have .

	implicit none

	character*(*) name
	character*(*) ext

	integer nall,next,n
	integer ichanm

	if( ext(1:1) /= '.' ) then
	  write(6,*) 'ext: ',ext
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(ext)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == ext ) then	!already there
	    !nothing
	  else
	    name(nall+1:) = ext
	  end if
	else
	  name(nall+1:) = ext
	end if

	end

c**************************************************************

	subroutine delete_extension(name,ext)

c deletes extension from file
c
c ext should have .

	implicit none

	character*(*) name
	character*(*) ext

	integer nall,next,n
	integer ichanm

	if( ext(1:1) /= '.' ) then
	  write(6,*) 'ext: ',ext
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(ext)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == ext ) then	!extension is there
	    name(n:) = ' '
	  end if
	end if

	end

c**************************************************************

	subroutine change_extension(name,extold,extnew)

c changes extension with new one
c
c ext should have .

	implicit none

	character*(*) name
	character*(*) extold,extnew

	integer nall,next,n
	integer ichanm

	if( extold(1:1) /= '.' .or. extnew(1:1) /= '.' ) then
	  write(6,*) 'extold: ',extold
	  write(6,*) 'extnew: ',extnew
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(extold)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == extold ) then	!extension is there
	    name(n:) = extnew
	  end if
	end if

	end

c**************************************************************

	subroutine subst_extension(name,ext,bforce)

c substitutes extension with given one
c
c if bforce is true substitutes given extension
c otherwise only adds if not already there
c
c extension must have 3 chars
c
c name		file name (with or without extension)
c ext		extension (with or without dot)
c bforce	force substitution of extension, even if already there

	implicit none

	character*(*) name
	character*(*) ext
	logical bforce

	integer nall,n,nstart,nend

        integer ichafs,ichanm

	nall = 1 + ichanm(name)

	n = nall - 4		!here should be the dot
	if( n .gt. 0 ) then	!has extension
	 if( name(n:n) .eq. '.' ) then	!has extension
	  if( bforce ) then	!substitute extension
	    nall = n
	  else
	    return		!leave extension
	  end if
	 end if
	end if

        nstart=ichafs(ext)
        nend=ichanm(ext)

	if( nend .gt. 0 ) then
	   if( ext(nstart:nstart) .ne. '.' ) then !add dot if not in ext
		name(nall:nall) = '.'
		nall = nall + 1
	   end if
	   name(nall:)=ext(nstart:nend)
	end if

	end

c**************************************************************

	subroutine check_extension(file,ext)

c finds extension of file and returns it

	implicit none

	character*(*) file
	character*(*) ext

	integer i

	ext = ' '

	i = index(file,'.',.true.)
	if( i == 0 ) return		!no extension

	ext = file(i+1:)

	end

c**************************************************************

	subroutine check_name_and_extension(file,name,ext)

c finds name and extension of file and returns it

	implicit none

	character*(*) file
	character*(*) name,ext

	integer i

	name = ' '
	ext = ' '

	i = index(file,'.',.true.)
	if( i == 0 ) return		!no extension

	name = file(:i-1)
	ext = file(i+1:)

	end

c**************************************************************

