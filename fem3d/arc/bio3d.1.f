
c***********************************************
c
c State variables used: (Wasp)
c
c nh3		71	1
c no3		72	2
c opo4		73	3
c phyto		74	4
c cbod		75	5
c do		76	6
c on		77	7
c op		78	8
c zoo		79	9
c
c State variables used: (Haka)
c
c php		81	1
c zoo		82	2
c det		83	3
c dop		84	4
c dip		85	5
c
c***********************************************

	subroutine bio3d(it,idt)

c eco-model cosimo

	implicit none

	include 'param.h'

	integer it	!time in seconds
	integer idt	!time step in seconds

	integer nstate,ntot
	parameter( nstate = 9 , ntot = 1 )

	integer narr
	parameter( narr = 100 )

	real bioarr(narr,nbcdim)	!array containing boundary state
	real bioaux(narr)		!aux array for boundaries

	real e(nlvdim,nkndim,nstate)	!state vector
	real eb(nlvdim,nkndim,nstate)	!boundary vector of state vectors
	real eload(nlvdim,nkndim,nstate)!vector of loadings
	real light(nlvdim,nkndim)	!vector of light climate
	save e,eb,eload

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

	character*80 bio2dn(1)
	common /bio2dn/bio2dn

	integer ilhkv(1)
	common /ilhkv/ilhkv

        real difv(0:nlvdim,1)
        common /difv/difv
        real cdifhv(1)
        common /cdifhv/cdifhv

        character*10 what

	integer k,i,l,lmax
	integer ibio
	real t,s,dt
c	FIXMED dichiaro u e v
	real u,v
	real eaux(nstate)
	real elaux(nstate)
	real einit(nstate)
	real ebound(nstate)
        real elini(nstate)
	integer icall,iunit
	integer itt,idtt,j
	real lux
	real dtt,dttday
	real area,vol
	real oxysat
	real getpar
	integer iround
	integer ieint,ipint

        integer mode
        real ai,lsurf

	integer ie,ii
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
        real mass

	integer iespecial,inspecial
	save iespecial,inspecial
	real rkpar,difmol
	save rkpar,difmol
	integer iub,itmcon,idtcon
	save iub,itmcon,idtcon
	save bioarr

	save einit,ebound,elini
        save icall

c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c	initial loadings                 [g/m**2/sec]		??
c
c        data einit /0.05, 0.4, 0.01, 0.05, 2.,   11.,0.2,0.01,0.015/
 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
c                  nh3,  no2, opo4, phyto, cbod,do, on, op    zoo   /


c                     phy     zoo     det     dop     dip
        !data einit  / 1.6E-3, 1.0E-4, 1.8E-3, 3.5E-3, 2.7E-3 /
c        data einit  / 0.,     0.,     0.,     0.,     0.     /
c        data ebound / 1.6E-3, 1.0E-4, 1.8E-3, 3.5E-3, 2.7E-3 /
c        data elini  / 0.,     0.,     0.,     0.,     0.     /
c
c------------------------------------------------------------------
	data icall /0/
c------------------------------------------------------------------

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  ibio = iround(getpar('ibio'))
	  if( ibio .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

c         --------------------------------------------------
c	  initialize state variables with einit
c         --------------------------------------------------

	  do k=1,nkn		!loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
	      do i=1,nstate
	        e(l,k,i) = einit(i)
              end do
	    end do
          end do

c         --------------------------------------------------
c	  set loadings in the interal areas
c         --------------------------------------------------

	  call setload(eload)

c         --------------------------------------------------
c	  set boundary conditions for all state variables
c         --------------------------------------------------

	  call bnds_init(bio2dn,2,nstate,narr,bioarr,ebound)	!new_section
	  !call bnds_print('initialization',narr,bioarr)

c         --------------------------------------------------
c	  initialize eco model
c         --------------------------------------------------

	  call eutroini

c         --------------------------------------------------
c	  parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

c         --------------------------------------------------
c	  meteorological forcing (HACK -> FIXME)
c         --------------------------------------------------

c	  tempair = 22.
c	  windspeed = 3.			!FIXME
c	  call wmeteo(tempair,windspeed)

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

	  iub = 55
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

          call confop(iub,itmcon,idtcon,nstate,'bio')

	  write(6,*) 'bio3d model initialized...'

	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

        what = 'lagvebio'

c	-------------------------------------------------------------------
c	time management
c	-------------------------------------------------------------------

	t0 = 0.
	dt = idt
	idtt = idt/ntot
	dtt = dt / ntot
	dttday = dtt / 86400.

	tsec = it
	tday = it / 86400. + t0		!time in days, FEM 0 is day t0

c	-------------------------------------------------------------------
c	boundary conditions	!new_section
c	-------------------------------------------------------------------

	  !call bnds_print('before set',narr,bioarr)
	call bnds_set(tsec,narr,bioarr,bioaux)
	  !call bnds_print('after set',narr,bioarr)
	call bnds_set_global(narr,bioarr,nkndim,nlvdim,nstate,eb,e)
	  !call bnds_print('after set global',narr,bioarr)

c        call getsurflight(tday,lsurf)
c        call setlight(nlvdim,nkndim,nstate,e,light,lsurf)

c	-------------------------------------------------------------------
c	loop on elements for biological reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
c	FIXMED ho tolto commento a call getuv
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)
            ai = light(l,k)                     !light climate

	    do i=1,nstate
	      eaux(i) = e(l,k,i)
	      elaux(i) = eload(l,k,i)
	    end do

	    do j=1,ntot
	      itt = it - idt + j*idtt
	     ! call eutro0d(tday,dttday,vol,d,vel,t,s,eaux,elaux)
              !call haka0d(tsec,dtt,vol,d,t,ai,eaux,elaux)
	    end do

	    do i=1,nstate
	      e(l,k,i) = eaux(i)
	    end do
          end do

	end do

c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	do i=1,nstate
          call scal3sh(what,e(1,1,i),eb(1,1,i),rkpar,cdifhv,difv,difmol)
          call massconc(e(1,1,i),nlvdim,mass)
          eaux(i) = mass
	end do
        write(88,*) it,(eaux(i),i=1,nstate)

c	-------------------------------------------------------------------
c	boundary (rivers)	!new_section
c	-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	write of results (file BIO)
c	-------------------------------------------------------------------

	do i=1,nstate
          call confil(iub,itmcon,idtcon,80+i,1,e(1,1,i))
	end do

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************

	subroutine writee(iunit,it,k,l,e,t,s,nlvdim,nkndim,nstate)

c formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nkndim,nstate
	real e(nlvdim,nkndim,nstate)
	real t(nlvdim,nkndim)
	real s(nlvdim,nkndim)

	integer i

	write(iunit,'(i10,11f12.4)') it,
     +			(e(l,k,i),i=1,nstate),
     +			t(l,k),
     +			s(l,k)

	end

c*************************************************************

	subroutine setload(eload)

c sets up eload which is loading for specified areas
c
c the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
c the specified loadings in areaload are in [kg/day]
c
c variables to be specified:
c
c nimmis        total number of areas for which loading is specified
c nodes         total number of nodes used to identify all areas
c karee         node numbers that specify the areas of loading
c iaree         area numbers for the nodes [1-nimmis]
c areaload      total loadings [kg/day] for areas
c
c the node numbers in karee are external node numbers

	implicit none

        include 'param.h'
	integer nstate
	parameter(nstate=9)

	real eload(3,neldim,nstate)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real ev(13,1)
	common /ev/ev
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	integer nimmis
	parameter (nimmis=5)
	real volaux(nimmis)
	real areaload(nimmis,nstate)

	save volaux,areaload

	integer aree(nkndim)

	integer nodes
	parameter(nodes=152)
	integer karee(nodes)
	integer iaree(nodes)
	save karee,iaree

	logical berror
	integer k,ie,ii,ia,i
	integer itype
	real area
	real litri,kgs
	real fact,rlfact
	real hdep,vol,load

	real getpar

c loading is kg/day
c
c	1=venezia-giudecca, 2=Murano-S. Erasmo 3=Burano
c	4=Fusina, 5=Campalto 
c
c	loading for areas [kg/day]
c 	
c	21 sett 2000

c	donata inserire cavallino

c carichi ripartiti secondo % franco

c 	data areaload /
c     + 1055.,70.,58.,792.,348. 			!nh3
c     + ,264.,17.6,14.61,198,87			!no3
c     +	,142.4,9.5,6.33,183,42			!opo4
c     + ,0.,0.,0.,0.,0.				!phyto
c     +	,7907.,518.,353.,2193.,915.		!cbod
c     +	,0.,0.,0.,0.,0.				!do
c     +	,0.,0.,0.,0.,0.				!on
c     +	,0.,0.,0.,0.,0.				!op
c     +	,0.,0.,0.,0.,0./			!zoo

c carichi ripartiti assegnando 50% alla forma ridotta dell'azoto e 50% alla
c forma ossidata

 	data areaload /
     +  659.,44.,37.,495.,217. 			!nh3
     + , 659.,44.,37.,495.,217. 		 	!no3
     +	,125.4,8.36,5.57,183.,61.		!opo4
     + ,0.,0.,0.,0.,0.				!phyto
     +	,7907.,518.,353.,2193.,915.		!cbod
     +	,0.,0.,0.,0.,0.				!do
     +	,0.,0.,0.,0.,0.				!on
     +	,0.,0.,0.,0.,0.				!op
     +	,0.,0.,0.,0.,0./			!zoo
c
        data karee /
     : 2112,2113,2863,2858,2114,2508,2507,2506,2505,2504,
     : 2503,2502,2501,2500,2499,2598,2497,2496,2632,2494,
     : 2493,2492,2491,2490,2489,2488,2487,2486,2485,2484,
     : 2447,2555,4235,2448,2449,2133,2132,2131,2130,2129,
     : 2128,2125,2124,2123,2122,4211,2109,2379,2108,2111,
     : 2112,2139,2138,2137,2136,2135,2134,2450,2451,2452,
     : 2453,2454,2455,2144,2143,2142,2141,2140,
     :
     : 3001,3105,2999,2998,2997,3096,2995,2994,3007,3084,
     : 3006,3005,3004,3003,3002,3070,2989,2988,3013,3014,
     : 3015,3016,3017,3018,3019,3020,3301,3591,3299,3603,
     : 3602,3295,3294,4358,3293,3366,2537,2536,2535,2534,
     : 2533,

     : 3041,3040,3039,3197,3192,3188,3029,3030,4210,3309,
     : 3308,3307,3306,3305,3304,3023,3024,3025,3026,3032,
     : 3044,3043,3042,3047,3048,3049,3050,3844,3805,3844,
     : 3804,3319,3318,3317,3316,3577,3576,3313,3312,3311,
     : 3310,
     :
     : 2172,            ! fusina 2 maggio 2000
     :
     : 2954             !campalto 2 maggio 2000
     :
     :  /       !node in area
c

	data iaree /
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,
     : 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     : 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     : 2,
     : 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     : 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     : 3,
     : 4,		!4,4,
     : 5
     :  /	!type of area
c
	litri = 1000.	!litri in m*3
	kgs  = 10.e+6	!mg in kg

c	rlfact = getpar('flbio')
	rlfact = 1.

	fact = rlfact*kgs/litri         ! [kg/m**3] -> [mg/l]
     
c extern to intern

	call n2int(nodes,karee,berror)

	if( berror) stop 'error stop: loading'

c intialize

	do i=1,nimmis
	  volaux(i) = 0.
	end do

	do k=1,nkn
	  aree(k) = 0
	end do

	do i=1,nodes
	  k = karee(i)
	  itype = iaree(i)
	  aree(k) = itype
	end do

c compute total volume for all areas given -> store in volaux

	do ie=1,nel
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hdep = hm3v(ii,ie)
	    ia = aree(k)
	    if( ia .gt. nimmis ) stop 'error stop ia'
	    if( ia .gt. 0 ) then
		volaux(ia) = volaux(ia) + area * hdep
	    end if
	  end do
	end do

c compute and set loading in eload [g/(m**3 day)] == [mg/(l day]

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ia = aree(k)
	    vol = volaux(ia)
	    do i=1,nstate
	      load = fact * areaload(ia,i) / vol
	      if( ia .le. 0 ) load = 0.
	      eload(ii,ie,i) = load
	    end do
	  end do
	end do

	end

c*************************************************************
c*************************************************************
c*************************************************************

        subroutine setsedload(nlvdim,nkndim,nstate,eload,elini)

c sets up sediment loading

        implicit none

        integer nlvdim,nkndim,nstate
        real eload(nlvdim,nkndim,nstate)
        real elini(nstate)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real ilhkv(1)
	common /ilhkv/ilhkv

        integer mode,i,k,l,lmax
        real d,vol,area

        mode = -1

        do i=1,nstate
         do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            eload(l,k,i) = 0.
          end do
          call dvanode(lmax,k,mode,d,vol,area)   !gets depth, volume and area
          eload(lmax,k,i) = elini(i) * area
         end do
        end do

        end

c*************************************************************

        subroutine setlight(nlvdim,nkndim,nstate,e,light,lsurf)

c sets up light climate

        implicit none

        integer nlvdim,nkndim,nstate
        real e(nlvdim,nkndim,nstate)
        real light(nlvdim,nkndim)
        real lsurf

        real delta,rpchla
        parameter( delta  = 0.1 )
        parameter( rpchla = (106.*12.)/(31.*50.) )

c redfield ratio  C:N:P  = 106:16:1
c C to chla ratio C/chla = 50
c ato,ic weight   C = 12   P = 31   => C/P = 106 * 12 / 31
c => chla/P = chla/C C/P = (1/50) * 106*12/31

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real ilhkv(1)
	common /ilhkv/ilhkv

        integer mode,i,k,l,lmax
        real d,vol,area
        real ltop,lmed,lbot
        real phy,zoo,det
        real chla,rk

        real depnode

        mode = -1

        do k=1,nkn
          lmax = ilhkv(k)
          ltop = lsurf
          do l=1,lmax
            d = depnode(l,k,mode)
            phy = e(l,k,1)
            zoo = e(l,k,2)
            det = e(l,k,3)
            chla=(phy+delta*(zoo+det))*rpchla
            rk=0.04+0.054*chla**(0.6667)+0.0088*chla
            lbot = ltop * exp( -rk * d )
            light(l,k) = 0.5 * ( ltop + lbot )
            ltop = lbot
          end do
        end do

        end

c*************************************************************

        subroutine getsurflight(tday,lsurf)

        implicit none

        real tday
        real lsurf

        integer ndata
        parameter (ndata=730)

        integer i,j
        real rad
        real rlight(ndata)
        save rlight

        integer icall
        save icall
        data icall / 0 /

        if( icall .eq. 0 ) then
            open(1,file='light.dat',status='old',form='formatted')
            do i=1,ndata
                read(1,*) j,rad
                if( j .ne. i ) then
                    write(6,*) 'error reading file: ',i,j,rad
                    stop 'error stop getsurflight'
                end if
                rlight(i) = rad
            end do
            close(1)
            icall = 1
        end if

        j = tday + 1.
        if( j .lt. 1 .or. j .gt. ndata ) then
            write(6,*) 'error in index: ',tday,ndata,j
            stop 'error stop getsurflight'
        end if

        lsurf = rlight(j)

        end

c*************************************************************

