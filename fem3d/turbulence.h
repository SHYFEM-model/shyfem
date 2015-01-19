
      real tken(0:nlvdim,nkndim)      !turbulent kinetic energ
      common /tken/tken
      real eps(0:nlvdim,nkndim)        !dissipation rate
      common /eps/eps
      real rls(0:nlvdim,nkndim)        !length scale
      common /rls/rls

	save /tken/,/eps/,/rls/

