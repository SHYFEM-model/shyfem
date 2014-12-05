
	integer iwave			!call for wave model
	integer iwwm			!call for coupling with wwm
	integer idcoup			!time step for sincronizing with wwm [s]

	common /wavcst/ iwave,iwwm,idcoup

        real waveh(nkndim)      	!sign. wave height [m]
        common /waveh/waveh

        real wavep(nkndim)      	!mean wave period [s]
	common /wavep/wavep

        real wavepp(nkndim)      	!peak wave period [s]
	common /wavepp/wavepp

        real waved(nkndim)      	!man wave direction
	common /waved/waved

        real waveov(nkndim)     	!wave orbital velocity
	common /waveov/waveov

        real wavefx(nlvdim,neldim)      !wave forcing terms
	common /wavefx/wavefx

        real wavefy(nlvdim,neldim)
	common /wavefy/wavefy
