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
c 20.10.2011    ggu     fx renamed to flux2d, new flux3d and vel3d_ie
c 16.12.2011    ggu     new array id_body, new vars idbdy, lunit
c
c******************************************************************

c-------------------------------------------------- main variables

        integer nbdy				!total number of particles
        integer idbdy				!max id used
        integer lunit				!unit for messages
        common /lagr_i/nbdy,idbdy,lunit

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

 	real lgr_Ww(nbdydim)			!custom variable
 	common /lgr_Ww/lgr_Ww

        integer id_body(nbdydim)		!id of particle
        common /id_body/id_body

c-------------------------------------------------- random walk
	
	real rwhpar				!horizontal diffusivity
	common /rwhpar/rwhpar
	
	real rwhvar(neldim)			!horizontal diffusivity (vary)
	common /rwhvar/rwhvar
		
c-------------------------------------------------- backtracking

	integer nback
	logical bback
	common /backtrack/ nback, bback

	logical bsurface			!keep particles on surface
	parameter (bsurface=.true.)

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
	
        real flux2d(3,neldim)			!fluxes of sides
        common /flux2d/flux2d

        real flux3d(nlvdim,3,neldim)		!fluxes of sides (3d)
        common /flux3d/flux3d

        real vel_ie(3,neldim)			!velocities of sides
        common /vel_ie/vel_ie

        real vel3d_ie(nlvdim,3,neldim)		!velocities of sides (3d)
        common /vel3d_ie/vel3d_ie

        real dvert(3,neldim)
        common /dvert/dvert

c-------------------------------------------------- special variables

	real tdecay				!decay time - do not use
	common /tdecay/tdecay

	real fall				!vertical sinking velocity
	common /fall/fall

	integer artype				!special element type
	common /artype/artype

c------------------------------------------------------------

	integer i_count(neldim)
	double precision t_count(neldim)

	common /i_count/i_count
	common /t_count/t_count
	save /i_count/,/t_count/

c------------------------------------------------------------

	integer*8 lgr_bitmap_in(nbdydim)	!bitmap for entry
	integer*8 lgr_bitmap_out(nbdydim)	!bitmap for leave
	common /lgr_bitmap/lgr_bitmap_in,lgr_bitmap_out

c------------------------------------------------------------

