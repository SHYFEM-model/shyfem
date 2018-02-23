
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

