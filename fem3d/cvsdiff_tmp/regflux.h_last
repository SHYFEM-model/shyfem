
c	header file for regular global matrices

	integer nxdim,nydim
	parameter( nxdim = 140 , nydim = 170 )

	integer maxind
	parameter( maxind = 5000 )

	real dxreg,dyreg
	parameter ( dxreg = 250. , dyreg = 250. )

	real value(maxind)

c-----------------------------------------------------------------------
c fflux( , ,1)		fluxes through horizontal sides
c fflux( , ,2)		fluxes through vertical sides
c alength( , ,1)	length of upper and lower sides (y)
c alength( , ,2)	length of right and left sides (x)

        real fflux(0:nxdim,0:nydim,2)   !global regular array for fluxes
        real alength(0:nxdim,0:nydim,2) !global regular array for sides

	common /fflux/fflux
	common /alength/alength

c-----------------------------------------------------------------------

        real aflux(nxdim,nydim)         !global regular array for area
        integer icount(nxdim,nydim)	!intersection count

	common /aflux/aflux
	common /icount/icount

c-----------------------------------------------------------------------

	real hmed(nxdim,nydim)		!water depth
	real zlev(nxdim,nydim)		!water level
	real areav(nxdim,nydim)		!area of boxes

	real dxside(nxdim,nydim)	!dx to use for computations
	real dyside(nxdim,nydim)	!dy to use for computations

	common /hmed/hmed
	common /zlev/zlev
	common /areav/areav
	common /dxside/dxside
	common /dyside/dyside

c-----------------------------------------------------------------------

	integer nboxdim
	parameter( nboxdim = 10000 )

	integer iboxei(nxdim,nydim)	!number of active boxes
	integer iboxes(2,nboxdim)	!coordinates of active boxes
	integer ibx(nboxdim)		!actual pointers...
	integer iby(nboxdim)		!actual pointers...

	common /iboxei/iboxei
	common /iboxes/iboxes

	integer indhmx
	common /indhmx/indhmx

	integer gindh(nydim,nxdim)
	integer gindsh(0:nxdim,0:nydim)

	common /gindh/gindh
	common /gindsh/gindsh

