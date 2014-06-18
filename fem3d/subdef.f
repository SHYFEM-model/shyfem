c
c $Id: subdef.f,v 1.11 2008-04-01 11:58:12 georg Exp $
c
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
c 23.05.1997	ggu     $$EXTENS - default extension may be overwritten
c 18.06.1997	ggu     restructured - idefna,idefop,idefts,idefun deleted
c 16.01.1998	ggu     idefop reintroduced -> to avoid link error
c 21.01.1998	ggu     in mkname: give extension with or without dot
c 08.08.2000	ggu     new routine ifemop
c 27.11.2001	ggu     error message rewritten
c 11.10.2002	ggu	new subroutine deffile
c 07.03.2007	ggu	new routine ifem_open_file
c 29.04.2010	ggu	new routine ifem_open_file
c 03.05.2010	ggu	new routine ifem_choose_file() and add_extension()
c 02.07.2011	ggu	idefna,idefop finally deleted
c 13.07.2011    ggu     cleaned from old structures
c 18.08.2011    ggu     bug fix in idefbas -> use status passed in
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

	call add_extension(file,ext,.false.)

	return

	end

c**************************************************************
c**************************************************************
c**************************************************************

	blockdata default_names

        character*80 def_bas,def_nam
	common /defdef/ def_bas,def_nam
	save /defdef/

	data def_bas,def_nam /' ',' '/

	end

c**************************************************************

	subroutine set_default_names(defbas,defnam)

c sets default names to be used with def routines
c
c must be called before any other routine in subdef can be used

        character*(*) defbas,defnam
	
        character*80 def_bas,def_nam
	common /defdef/ def_bas,def_nam
	save /defdef/

	def_bas = defbas
	def_nam = defnam

	end

c**************************************************************

	subroutine get_default_names(defbas,defnam)

c gets default names to be used with def routines

        character*(*) defbas,defnam
	
        character*80 def_bas,def_nam
	common /defdef/ def_bas,def_nam
	save /defdef/

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

        implicit none

        character*(*) ext,file

        character*80 def_bas,def_nam
	common /defdef/ def_bas,def_nam
	save /defdef/

        character*80 dir,name

	name = def_nam
        call getfnm('datdir',dir)	! this has to be deleted
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
        if( defdir .ne. ' ' ) call getfnm(defdir,dir)
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
        call add_extension(name,'.bas',.true.)

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
	call add_extension(file,ext,.true.)

	ifem_choose_file = ifileo(0,file,form,status)

	if( ifem_choose_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine add_extension(name,ext,bforce)

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
	if( n .gt. 0 .and. name(n:n) .eq. '.' ) then	!has extension
	  if( bforce ) then	!substitute extension
	    nall = n
	  else
	    return		!leave extension
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

