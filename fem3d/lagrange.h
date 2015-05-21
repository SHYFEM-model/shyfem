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
c*****************************************************************

c--------------------------------------------------
c main parameters
c--------------------------------------------------

	logical blgrxi				!new version with xi coords
	parameter (blgrxi=.true.)

	logical blgrdebug,blgrsurf
	common /lagr_b/blgrdebug,blgrsurf
	save /lagr_b/

        integer ilagr				!type of lagrangian simulation
        integer nbdy				!total number of particles
        integer idbdy				!max id used
        integer lunit				!unit for messages
        integer ipvert				!vertical release
        common /lagr_i/ilagr,nbdy,idbdy,lunit,ipvert
	save /lagr_i/

c--------------------------------------------------
c special variables
c--------------------------------------------------

	real azlgr				!az parameter
	common /azlgr/azlgr

	real tdecay				!decay time - do not use
	common /tdecay/tdecay

	real fall				!vertical sinking velocity
	common /fall/fall

	integer artype				!special element type
	common /artype/artype

c--------------------------------------------------
c horizontal diffusion
c--------------------------------------------------
	
	real rwhpar				!horizontal diffusivity
	common /rwhpar/rwhpar
	
	real rwhvar(neldim)			!horizontal diffusivity (vary)
	common /rwhvar/rwhvar
		
c--------------------------------------------------
c arrays
c--------------------------------------------------

        type  :: entry
	  sequence
	  double precision :: xi(3)		!internal coordinate
          double precision :: sv  		!sinking velocity
          real    :: xst			!initial x-coordinate
          real    :: yst			!initial y-coordinate
          real    :: zst			!initial z-coordinate
          real    :: x  			!x-coordinate
          real    :: y  			!y-coordinate
          real    :: z  			!z-coordinate
          real    :: c   			!custom property for plot
	  real    :: tin			!time of release
          integer :: est			!initial element number
          integer :: id				!id of particle
          integer :: ty				!type of particle
          integer :: ie 			!element number
          integer :: l				!layer number
          integer :: dummy			!dummy argument for sequence
        end type entry

        type(entry) :: lgr_ar(nbdydim)

	common /lgr_cc/lgr_ar
	save /lgr_cc/

	logical bsedim, blarvae, boilsim
	common /lgrlogic/ bsedim, blarvae, boilsim

c--------------------------------------------------
c backtracking
c--------------------------------------------------

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

c--------------------------------------------------
c fluxes and velocities
c--------------------------------------------------
	
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

c--------------------------------------------------
c conncectivity
c--------------------------------------------------

	integer i_count(neldim)
	double precision t_count(neldim)

	common /i_count/i_count
	common /t_count/t_count
	save /i_count/,/t_count/

c--------------------------------------------------

	integer*8 lgr_bitmap_in(nbdydim)	!bitmap for entry
	integer*8 lgr_bitmap_out(nbdydim)	!bitmap for leave
	common /lgr_bitmap/lgr_bitmap_in,lgr_bitmap_out

c--------------------------------------------------


