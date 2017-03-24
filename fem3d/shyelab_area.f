
!***************************************************************

        subroutine handle_area

        use elabutil

        implicit none

        integer		  :: np
        real		  :: xdummy(1),ydummy(1)
	real, allocatable :: xl(:),yl(:)

        ieflag = 1
        ikflag = 1

        if( barea ) then
          np = 0
          call get_line_list(areafile,np,xdummy,ydummy)
          allocate(xl(np))
          allocate(yl(np))
          call get_line_list(areafile,np,xl,yl)
          call check_elements(np,xl,yl,ieflag,ikflag)
	  baverbas = .true.
          deallocate(xl,yl)
        end if

	end subroutine handle_area

!***************************************************************
!
! Read coordinates of points in line
! for n == 0 only checks how many nodes to read
! for n > 0 reads nodes into nodes() (error if n is too small)
!
!***************************************************************

        subroutine get_line_list(file,n,xl,yl)

        implicit none

        character*(*) file
        integer		:: n
        real 		:: xl(n),yl(n)

	real		:: xaux, yaux
        logical 	:: btest
        integer		:: ifileo
        integer		:: nin,ios,ndim

        nin = ifileo(0,file,'form','old')
        if( nin .le. 0 ) goto 99

        ndim = n
        btest = ndim == 0

        n = 0
        do
          read(nin,*,iostat=ios) xaux,yaux
          if( ios > 0 ) goto 98
          if( ios < 0 ) exit
          n = n + 1
          if( .not. btest ) then
            if( n > ndim ) goto 95
            xl(n) = xaux
            yl(n) = yaux
          end if
        end do

        if( n == 0 ) goto 96
        if( n < 3 ) goto 97

        close(nin)

        return
   95   continue
        write(6,*) 'n,ndim :',n,ndim
        write(6,*) 'file: ',trim(file)
        stop 'error stop get_line_list: dimension error'
   96   continue
        write(6,*) 'no data in file ',trim(file)
        stop 'error stop get_line_list: read error'
   97   continue
        write(6,*) 'no enouth data in file ',trim(file)
        write(6,*) 'a polygon must contain at least 3 distinct points'
        stop 'error stop get_line_list: read error'
   98   continue
        write(6,*) 'read error in record ',n
        write(6,*) 'in file ',trim(file)
        stop 'error stop get_line_list: read error'
   99   continue
        write(6,*) 'file: ',trim(file)
        stop 'error stop get_line_list: cannot open file'

	end subroutine get_line_list

!***************************************************************

