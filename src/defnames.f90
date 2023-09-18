!
! $Id: subdef.f,v 1.11 2008-04-01 11:58:12 georg Exp $
!
! create default names
!
! contents :
!
! subroutine mkname(dir,name,ext,file)            makes file name
! subroutine defmak(defdir,defnam,ext,file)	  makes file with defaults
! function ideffi(defdir,defnam,ext,form,status)  opens file in default dir
! function ifemop(ext,form,status)		  opens file with default name
! function ifemopa(text,ext,form,status)	  opens file with default name
! function ifem_open_file(ext,status)		  opens unformated file
! function ifem_test_file(ext,status)		  tries to open unformated file
!
! revision log :
!
! 23.05.1997	ggu     $$EXTENS - default extension may be overwritten
! 18.06.1997	ggu     restructured - idefna,idefop,idefts,idefun deleted
! 16.01.1998	ggu     idefop reintroduced -> to avoid link error
! 21.01.1998	ggu     in mkname: give extension with or without dot
! 08.08.2000	ggu     new routine ifemop
! 27.11.2001	ggu     error message rewritten
! 27.11.2001    ggu     routine to handle info file (getinfo)
! 11.10.2002	ggu	new subroutine deffile
! 07.03.2007	ggu	new routine ifem_open_file
! 29.04.2010	ggu	new routine ifem_open_file
! 03.05.2010	ggu	new routine ifem_choose_file() and add_extension()
! 02.07.2011	ggu	idefna,idefop finally deleted
! 13.07.2011    ggu     cleaned from old structures
! 18.08.2011    ggu     bug fix in idefbas -> use status passed in
!
! notes :
!
! exchange ideffi with ifemop (many apperances)
! eliminate call to getpar/getfnm
!
!**************************************************************
        module defnames
!**************************************************************
!**************************************************************
!**************************************************************


        character*80,save :: def_bas,def_nam

	data def_bas,def_nam /' ',' '/

!**************************************************************

        contains

!**************************************************************

        subroutine mkname(dir,name,ext,file)

! makes file name given its constituents
!
! dir   directory
! name  name
! ext   extension (with or without dot)
! file  created file name (return)

        use utility
        implicit none

! arguments
        character*(*) dir,name,ext,file
! local
        integer nall,nstart,nend,naux

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

!**************************************************************

	subroutine set_default_names(defbas,defnam)

! sets default names to be used with def routines
!
! must be called before any other routine in subdef can be used

        implicit none

        character*(*) defbas,defnam
	
	def_bas = defbas
	def_nam = defnam

	end

!**************************************************************

	subroutine get_default_names(defbas,defnam)

! gets default names to be used with def routines

        implicit none

        character*(*) defbas,defnam
	
	defbas = def_bas
	defnam = def_nam

	end
	
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine def_make(ext,file)

! makes file with defaults supplied
!
! ext   extension (with dot)
! file  created file name (return)

        use para
        implicit none

        character*(*) ext,file
        character*80 dir,name

	name = def_nam
        call getfnm('datdir',dir)	! this has to be deleted
        call getfnm('runnam',name)	! this has to be deleted

	call mkname(dir,name,ext,file)

	end

!**************************************************************

        subroutine defmak(defdir,defnam,ext,file)

! FIXME -> take defdir,defnam from common block

! makes file with defaults supplied
!
! defdir   directory
! defnam  name
! ext   extension (with dot)
! file  created file name (return)

        use para
        implicit none

        character*(*) defdir,defnam,ext,file
	character*80 dir,name

	dir = ' '
        if( defdir .ne. ' ' ) call getfnm(defdir,dir)
        call getfnm(defnam,name)

	call mkname(dir,name,ext,file)

	end

!**************************************************************

        function ideffi(defdir,defnam,ext,form,status)

! FIXME -> substitute with ifemop (many appearances)

! opens file in default dir

! defdir   directory
! defnam  name
! ext   extension (with dot)
! form  formatted ?
! status open status

        use fil

        implicit none

	integer ideffi
        character*(*) defdir,defnam,ext,status,form
	character*80 file

	call defmak(defdir,defnam,ext,file)
        call def_make(ext,file)
	ideffi=ifileo(0,file,form,status)

	end

!**************************************************************

	function idefbas(basnam,status)

        use fil

	implicit none

	integer idefbas
	character*(*) basnam
	character*(*) status

	character*80 name

        name = basnam
        call add_extension(name,'.bas',.true.)

        idefbas=ifileo(0,name,'unform',status)

	end

!**************************************************************
!**************************************************************
!**************************************************************
! opening of default simulation
!**************************************************************
!**************************************************************
!**************************************************************

        function ifemop(ext,form,status)

! opens file with default name (run) and extension given for fem model
! returns with error code

! ext   extension (with dot)
! form  formatted ?
! status open status

        use fil

        implicit none

	integer ifemop
        character*(*) ext,status,form

	character*80 file,defdir,defnam

        call def_make(ext,file)

	ifemop=ifileo(0,file,form,status)

	end

!**************************************************************

!**************************************************************

        function ifemopa(text,ext,form,status)

! opens file with default name (run) and extension given for fem model
! in case of error exits with error message text

! text  error message
! ext   extension (with dot)
! form  formatted ?
! status open status

        use fil

        implicit none

	integer ifemopa
        character*(*) text,ext,status,form

	character*80 file,defdir,defnam

        call def_make(ext,file)
	ifemopa=ifileo(0,file,form,status)

	if( ifemopa .le. 0 ) then
	  write(6,*) 'error opening file ',file
	  write(6,*) text
	  stop 'error stop ifemopa'
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************
! unformatted opening of default simulation
!**************************************************************
!**************************************************************
!**************************************************************

        function ifem_open_file(ext,status)

! opens unformated file with default name (run) and extension given
! in case of error exits 

! ext		extension (with dot)
! status	open status

        implicit none

	integer ifem_open_file
        character*(*) ext,status

	ifem_open_file = ifem_test_file(ext,status)

	if( ifem_open_file .le. 0 ) then
	  stop 'error stop ifem_open_file'
	end if

	end

!**************************************************************

        function ifem_test_file(ext,status)

! tries to open unformated file with default name (run) and extension given

! ext		extension (with dot)
! status	open status

        use fil

        implicit none

	integer ifem_test_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form

	form = 'unform'
        call def_make(ext,file)

	ifem_test_file = ifileo(0,file,form,status)

	if( ifem_test_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

!**************************************************************

        function ifem_choose_file(ext,status)

! tries to open unformated file with default name (run) and extension given
! insists on extension -> if name has extension substitute it with ext

! ext		extension (with dot)
! status	open status

        use fil

        implicit none

	integer ifem_choose_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form

	form = 'unform'
        call def_make(ext,file)
	call add_extension(file,ext,.true.)

	ifem_choose_file = ifileo(0,file,form,status)

	if( ifem_choose_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine add_extension(name,ext,bforce)

! substitutes extension with given one
!
! if bforce is true substitutes given extension
! otherwise only adds if not already there
!
! extension must have 3 chars
!
! name		file name (with or without extension)
! ext		extension (with or without dot)
! bforce	force substitution of extension, even if already there

        use utility
	implicit none

	character*(*) name
	character*(*) ext
	logical bforce

	integer nall,n,nstart,nend

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

!**********************************************************************

        subroutine getinfo(iunit)

! gets unit of info file

        implicit none

        integer iunit

        integer iu
        save iu
        data iu / 0 /

        if( iu .le. 0 ) then
          iu = ifemop('.inf','formatted','new')
          if( iu .le. 0 ) then
            write(6,*) 'error in opening info file'
            stop 'error stop getinfo'
          end if
        end if

        iunit = iu

        end

!**************************************************************
        end module defnames
!**************************************************************
