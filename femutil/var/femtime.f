
	program test

	call set_include
	call set_module
	call write_include
	call write_module

	end

	module femtime
	include 'femtime.h'
	end module

	module femtime1
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	double precision t_act,dt_act,dt_orig,atime0,dtanf,dtend
        common /femtimd/ t_act,dt_act,dt_orig,atime0,dtanf,dtend
	logical bsync
        common /femtiml/ bsync
	character*20 aline_act
        common /femtimc/ aline_act
        save /femtim/,/femtimu/,/femtimd/,/femtiml/,/femtimc/
	end module

	subroutine set_include
	include 'femtime.h'
	it = 99
	end

	subroutine set_module
	use femtime
	dt = 88
	end

	subroutine write_include
	include 'femtime.h'
	write(6,*) 'include: ',it,dt
	end

	subroutine write_module
	use femtime1
	write(6,*) 'module: ',it,dt
	end

