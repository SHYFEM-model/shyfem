c
c $Id: loading.f,v 1.3 2008-03-17 18:15:35 georg Exp $
c
c loading routines
c
c revision log :
c
c 15.06.2004	ggu     written from scratch
c 28.06.2004	dmk     loading factor was wrong: kgs  = 1.e+6
c 28.06.2004	dmk     areaload loops now over state variables
c 14.03.2008	ggu     new routine set_surface_load
c
c notes :
c
c needed new routine init_load (sets to 0)
c in setload_new add to load matrix
c
c*************************************************************

	subroutine setload_new(eload)

c sets up eload which is loading for specified areas
c
c the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
c the specified loadings in areaload are in [kg/day]
c
c variables to be specified:
c
c nareas        total number of areas for which loading is specified
c nodes         total number of nodes used to identify all areas
c karee         node numbers that specify the areas of loading
c iaree         area numbers for the nodes [1-nareas]
c areaload      total loadings [kg/day] for areas
c
c  the node numbers in karee are external node numbers

	use levels
	use basin
	!use levels, only : nlvdi

	implicit none

        include 'param.h'
	integer nstate
	parameter(nstate=9)

	real eload(nlvdi,nkndi,nstate)



	integer nareas
	parameter (nareas=5)
	real volume(nareas)

	real areaload(nstate,nareas)
	save areaload

	integer, save, allocatable :: aree(:)

	integer k,l,lmax,i,ia
	real litri,kgs
	real fact,rlfact
	real load,vol

	real getpar

c loading is kg/day
c
c	1=venezia-giudecca, 2=Murano-S. Erasmo 3=Burano
c	4=Fusina, 5=Campalto 
c
c	loading for areas [kg/day]
c 	

c	donata inserire cavallino

c 	prova per controllare punto anomalo cavallino OK

c       data areaload/
c    +    0.,0.,0.,0.,0.,0.,0.,0.,0.
c    +    ,0.,0.,0.,0.,0.,0.,0.,0.,0.
c    +    ,0.,0.,0.,0.,0.,0.,0.,0.,0.
c    +    ,0.,0.,0.,0.,0.,0.,0.,0.,0.
c    +    ,10000.0,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.
c    +  /
        data areaload/
     +    1021.5,255.4,150.7,0.,8936.,0.,0.,0.,0.
     +   ,71.4,17.8,10.7,0.,595.7,0.,0.,0.,0.
     +   ,43.9,11.0,6.4,0.,397.2,0.,0.,0.,0.
     +   ,247.5,742.5,183.2,0.,3132.3,0.,0.,0.,0.
     +   ,82.5,247.5,61.1,0.,1044.1,0.,0.,0.,0.
     +  /

c       loading elaborazione dati 2003
c        data areaload /
c     +  1021.5, 71.4, 43.9, 247.5, 82.5,        !nh3
c     +  255.4,  17.8,   11.0,   742.5,  247.5,  !nox
c     +  150.7,  10.7,   6.4,    183.2,  61.1,   !opo4
c     +  0.,     0.,     0.,     0.,     0.,     !phyto
c     +  8936.0, 595.7,  397.2,  3132.3, 1044.1, !cbod
c     +  0.,     0.,     0.,     0.,     0.,     !do
c     +  0.,     0.,     0.,     0.,     0.,     !on
c     +  0.,     0.,     0.,     0.,     0.,     !op
c     +  0.,     0.,     0.,     0.,     0./     !zoo

        integer nnodes1
        parameter(nnodes1=68)
        integer nodes1(nnodes1)
        save nodes1
        data nodes1 /
     :     2112,2113,2863,2858,2114,2508,2507,2506,2505,2504,
     :     2503,2502,2501,2500,2499,2498,2497,2496,2632,2494,
     :     2493,2492,2491,2490,2489,2488,2487,2486,2485,2484,
     :     2447,2555,4235,2448,2449,2133,2132,2131,2130,2129,
     :     2128,2125,2124,2123,2122,4211,2109,2379,2108,2111,
     :     2112,2139,2138,2137,2136,2135,2134,2450,2451,2452,
     :     2453,2454,2455,2144,2143,2142,2141,2140
     :  /
        integer nnodes2
        parameter(nnodes2=41)
        integer nodes2(nnodes2)
        save nodes2
        data nodes2 /
     :     3001,3105,2999,2998,2997,3096,2995,2994,3007,3084,
     :     3006,3005,3004,3003,3002,3070,2989,2988,3013,3014,
     :     3015,3016,3017,3018,3019,3020,3301,3591,3299,3603,
     :     3602,3295,3294,4358,3293,3366,2537,2536,2535,2534,
     :     2533
     :  /
        integer nnodes3
        parameter(nnodes3=41)
        integer nodes3(nnodes3)
        save nodes3
        data nodes3 /
     :     3041,3040,3039,3197,3192,3188,3029,3030,4210,3309,
     :     3308,3307,3306,3305,3304,3023,3024,3025,3026,3032,
     :     3044,3043,3042,3047,3048,3049,3050,3844,3805,3844,
     :     3804,3319,3318,3317,3316,3577,3576,3313,3312,3311,
     :     3310
     :  /
        integer nnodes4
        parameter(nnodes4=1)
        integer nodes4(nnodes4)
        save nodes4
        data nodes4 /
     :     2172            ! fusina 2 maggio 2000
     :  /
        integer nnodes5
        parameter(nnodes5=1)
        integer nodes5(nnodes5)
        save nodes5
        data nodes5 /
     :     2954             !campalto 2 maggio 2000
     :  /

        integer icall
        save icall
        data icall / 0 /

c---------------------------------------------------------
c initialization
c---------------------------------------------------------

        if( icall .eq. 0 ) then
          icall = 1

	  allocate(aree(nkn))

          call load_init_area(nkn,aree)

          call load_add_area(1,nnodes1,nodes1,aree)
          call load_add_area(2,nnodes2,nodes2,aree)
          call load_add_area(3,nnodes3,nodes3,aree)
          call load_add_area(4,nnodes4,nodes4,aree)
          call load_add_area(5,nnodes5,nodes5,aree)
        end if

c---------------------------------------------------------
c conversion factors
c---------------------------------------------------------

	litri = 1000.	!litri in m*3
	kgs  = 1.e+6	!mg in kg               -> bug

c	rlfact = getpar('flbio')
	rlfact = 1.

	fact = rlfact*kgs/litri         ! [kg/m**3] -> [mg/l]
     
c---------------------------------------------------------
c compute volumes
c---------------------------------------------------------

        call load_make_volume(nareas,volume,aree)

c---------------------------------------------------------
c compute and set loading in eload [g/(m**3 day)] == [mg/(l day]
c---------------------------------------------------------

        do k=1,nkn
	  ia = aree(k)
          if( ia .gt. 0 ) then
	    vol = volume(ia)
            lmax = ilhkv(k)
	    do i=1,nstate
              do l=1,lmax
	        load = fact * areaload(i,ia) / vol
	        eload(l,k,i) = load	!FIXME -> should be added to eload
              end do
	    end do
          end if
	end do

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c*************************************************************

	subroutine set_surface_load(eload,sload)

c sets up eload which is loading for total basin
c
c eload contains computed loads on return
c sload contains atmospheric loadings (one value for each state variable)
c
c the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
c the specified loadings in sload are in [kg/day] or [kg/(m**2 day)]
c
c please set afact according to the choice of unit of sload (see below)

	use basin
	use levels, only : nlvdi

	implicit none

        include 'param.h'

	integer nstate
	parameter(nstate=9)

	real eload(nlvdi,nkndi,nstate)	!loading matrix for eutro
	real sload(nstate)			!surface load for each var


	real, save, allocatable :: areav(:)
	real, save, allocatable :: volv(:)

	integer k,i,mode,layer
	real litri,kgs
	real fact,afact,rlfact
	real load
	real area,volume
	double precision areatot

	real areanode, volnode

	save areatot,fact

        integer icall
        save icall
        data icall / 0 /

c---------------------------------------------------------
c initialization
c---------------------------------------------------------

	layer = 1			! only surface

        if( icall .eq. 0 ) then
          icall = 1
	  mode = 1

	  allocate(areav(nkn))
	  allocate(volv(nkn))

	  areatot = 0.
          do k=1,nkn
              areav(k) = areanode(k,mode)
	      volv(k) = volnode(layer,k,mode)
	      areatot = areatot + areav(k)
          end do

	  litri = 1000.			! litri in m*3
	  kgs  = 1.e+6			! mg in kg
	  fact = kgs/litri		! [kg/m**3] -> [mg/l]
        end if

c---------------------------------------------------------
c compute load and add to matrix  [g/(m**3 day)] == [mg/(l day)]
c---------------------------------------------------------

	!afact = 1.			! loadings in [kg/(m**2 day)]
	afact = 1. / areatot		! loadings in [kg/day]

        do k=1,nkn
	    area = areav(k)
	    volume = volv(k)
	    rlfact = fact * afact * area / volume
	    do i=1,nstate
	        load = rlfact * sload(i)
	        eload(layer,k,i) = eload(layer,k,i) + load
	    end do
	end do

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end

c*************************************************************
c*************************************************************
c*************************************************************

        subroutine load_init_area(nkn,aree)

c initializes area array

        implicit none

        integer nkn
        integer aree(1)

        integer k

        do k=1,nkn
          aree(k) = 0
        end do

        end

c*************************************************************

        subroutine load_add_area(iarea,n,nodes,aree)

c initializes area array

        implicit none

        integer iarea           !number of area
        integer n               !total number of nodes
        integer nodes(1)        !nodes that make up area
        integer aree(1)         !areas of each node

        integer i,k
        logical berror

	call n2int(n,nodes,berror)
	if( berror ) stop 'error stop: load_init_area'

        do i=1,n
	  k = nodes(i)
	  aree(k) = iarea
        end do

        end

c*************************************************************

        subroutine load_make_volume(nareas,volume,aree)

c makes total volume of areas

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nareas
        real volume(1)
        integer aree(1)

	include 'param.h'

        integer mode,i,k,ia,lmax,l

        real volnode

        mode = +1               !new time level

	do i=1,nareas
	  volume(i) = 0.
	end do

        do k=1,nkn
	  ia = aree(k)
	  if( ia .gt. nareas ) stop 'error stop ia'
	  if( ia .gt. 0 ) then
            lmax = ilhkv(k)
            do l=1,lmax
              volume(ia) = volume(ia) + volnode(l,k,mode)
            end do
          end if
        end do

        end

c*************************************************************


