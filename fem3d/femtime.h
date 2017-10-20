
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig

	double precision t_act,dt_act,dt_orig,atime0
        common /femtimd/ t_act,dt_act,dt_orig,atime0

	logical bsync
        common /femtiml/ bsync

        save /femtim/,/femtimu/,/femtimd/,/femtiml/

