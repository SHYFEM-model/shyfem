c
c $Id: lagrange.h,v 1.6 2009-03-11 16:25:59 georg Exp $
c
c lagrangian header
c
c revision log :
c
c 29.04.2005    ggu     arrays for backtracking
c 29.11.2006    ggu     new version integrated into main model
c 15.02.2009    ggu     cleaned version
c
c******************************************************************

c-------------------------------------------------- main variables

        integer nbdy				!total number of particles
        common /nbdy/nbdy

        integer est(nbdydim)			!initial element number
        common /est/est

        real xst(nbdydim)			!initial x-coordinate
        common /xst/xst

        real yst(nbdydim)			!initial y-coordinate
        common /yst/yst

	real zst(nbdydim)			!initial z-coordinate
	common /zst/zst

        integer ie_body(nbdydim)		!element number
        common /ie_body/ie_body

        real x_body(nbdydim)			!x-coordinate
        common /x_body/x_body

        real y_body(nbdydim)			!y-coordinate
        common /y_body/y_body

	real z_body(nbdydim)			!z-coordinate
	common /z_body/z_body

	real tin(nbdydim)			!time of release
	common /tin/tin

 	real lgr_var(nbdydim)			!custom variable
 	common /lgr_var/lgr_var

c-------------------------------------------------- random walk
	
	real rwhpar				!horizontal diffusivity
	common /rwhpar/rwhpar
	
	real rwhvar(neldim)			!horizontal diffusivity (vary)
	common /rwhvar/rwhvar
		
c-------------------------------------------------- backtracking

	integer nback
	logical bback
	common /backtrack/ nback, bback

	real v_lag(neldim), u_lag(neldim)	!backtracking velocity
	common /v_lag/v_lag, /u_lag/u_lag

        integer ie_back(neldim)			!backtracking element number
        common /ie_back/ie_back

        real x_back(neldim)			!backtracking x-coordinate
        common /x_back/x_back

        real y_back(neldim)			!backtracking y-coordinate
        common /y_back/y_back

        real z_back(neldim)			!backtracking z-coordinate
	common /z_back/z_back

c-------------------------------------------------- fluxes and velocities
	
        real fx(3,neldim)			!fluxes of sides
        common /fx/fx

        real vel_ie(3,neldim)			!velocities of sides
        common /vel_ie/vel_ie

        integer lt_body(nbdydim)		!side of parting particle
        common /lt_body/lt_body

        real dvert(3,neldim)
        common /dvert/dvert

c-------------------------------------------------- special variables

	real tdecay				!decay time
	common /tdecay/tdecay

	real fall				!vertical sinking velocity
	common /fall/fall

	integer artype				!special element type
	common /artype/artype

c------------------------------------------------------------

