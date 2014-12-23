
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig

	double precision t_act,dt_act
        common /femtimd/ t_act,dt_act

        save /femtim/,/femtimu/,/femtimd/

