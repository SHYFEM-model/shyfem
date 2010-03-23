c
c $Id: subwaves.f,v 1.10 2008-10-10 09:29:54 georg Exp $
c
c waves subroutines
c
c description :
c
c This routine is used to calculate the wave height and period
c from wind speed, fetch and depth using the EMPIRICAL PREDICTION
c EQUATIONS FOR SHALLOW WATER (Shore Protection Manual, 1984).
c It is a part of the sediment transport module.
c It considers a homogeneous wind field all over the domain.
c
c Copyright: Christian Ferrarin - ISMAR-CNR
c
c revision log :
c
c 14.07.2003	ccf	written from scratch by Christian Ferrarin ISMAR-CNR
c 01.09.2003    ccf     add subroutine e2n2d (element to node value)
c 01.11.2004    ccf     compute averaged depth along the fetch
c 11.11.2004    ccf     limiting wave height: Hbr = 0.50 h
c 11.11.2004    ccf     bug fix in make_stress: compute tm
c 26.11.2004    ggu     no idt passed to subwaves
c 05.01.2005    ccf     output per nodes
c 07.11.2005    ggu     changes from Christian integrated
c 18.10.2006    ccf     modified, walues on nodes
c 10.04.2008    ccf     check for p > 0
c 09.10.2008    ggu     new call to confop
c
c notes :
c
c Hs            significant wave height
c Tm            mean period
c Hrms          rms height (Hrms = Hs/sqrt(2))
c Tp            peak period (Tm = 0.781 Tp)
c
c Tm = 11 sqrt(Hs/g)    if no info on Tm is available
c
c for sediment transport use T=Tp and Uw = sqrt(2) Urms
c
c Uw            orbital velocity
c Urms          std. dev. of orbital velocity
c
c Dispersion relation:
c
c o             omega, frequency, o = 2 pi / T
c k             wave number, k = 2 pi / L
c L             wave length
c
c o**2 = gk tanh(kh)    dispersion relation
c
c zeta = o**2 h / g
c eta = k h
c
c ==> zeta = eta tanh(eta)
c
c easy computation (better than 1 %)
c
c zeta < 1      eta = sqrt(zeta) * (1 + 0.2 * zeta)
c zeta > 1      eta = zeta( 1 + 0.2 * exp(2 - 2 * zeta)
c
c Orbital velocity:
c
c Uw = pi H / (T sinh(kh))
c
c**************************************************************

        subroutine subwaves(it,dt)

        implicit none

        include 'param.h'

        integer it
	real dt

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real wxnv(1),wynv(1)	!x and y wind component [m/s]
        common /wxnv/wxnv,/wynv/wynv

c --- input variable
        real wis		!wind speed at 10m [m/s]
        real wid		!wind direction [degree north]
        real fet(neldim)        !wind fetch length [m]
        real daf(neldim)        !averaged depth along the fetch [m]

c --- output variable
        real waveh(nkndim)	!wave height [m]
        real wavep(nkndim)	!wave period [s]
        real waved(nkndim)	!wave direction (same as wind direction)
        common /waveh/waveh, /wavep/wavep, /waved/waved

        real waeh(neldim)	!wave height [m]
        real waep(neldim)	!wave period [s]
        real waed(neldim)	!wave direction (same as wind direction)
c        common /waeh/waeh, /waep/waep, /waed/waed

c --- stress variables
        real tcv(neldim)
        real twv(neldim)
        real tmv(neldim)

        real v1v(1)
        common /v1v/v1v
        real v2v(1)
        common /v2v/v2v

c --- local variable
        real widold		!old wind direction
        real depele             !element depth function [m]
        real hbr		!limiting wave height [m]
        real dep,depe
        real gh,gx,hg
        integer ie,icount

        save widold
        real g			!gravity acceleration [m2/s]
        parameter (g=9.81)
        real z0
        parameter (z0=5.e-4)
        real ah1,ah2,ah3,eh1,eh2,eh3,eh4
        real at1,at2,at3,et1,et2,et3,et4
        real auxh,auxh1,auxt,auxt1

c------------------------------------------------------ Hurdle and Stive
c        parameter(ah1=0.25,ah2=0.6,eh1=0.75)
c        parameter(eh2=0.5,ah3=4.3e-05,eh3=1.,eh4=2.)
c        parameter(at1=8.3,at2=0.76,et1=0.375)
c        parameter(et2=1./3.,at3=4.1e-05,et3=1.,et4=3.)
c------------------------------------------------------ SPM
        parameter(ah1=0.283,ah2=0.53,eh1=3./4.)
        parameter(eh2=1.,ah3=0.00565,eh3=1./2.,eh4=1.)
        parameter(at1=7.54,at2=0.833,et1=3./8.)
        parameter(et2=1.,at3=0.0379,et3=1./3.,et4=1.)

        real getpar
        integer iwave 		!call parameter
        integer ius,itmcon,idtcon
        save ius,itmcon,idtcon

        integer icall		!initialization parameter
        save icall
        data icall /0/

c ----------------------------------------------------------
c Initialization
c ----------------------------------------------------------

        if( icall .le. -1 ) return

        if( icall .eq. 0 ) then

c         --------------------------------------------------
c         Initialize state variables
c         --------------------------------------------------

          do ie = 1,nel
            waeh(ie) = 0.
            waep(ie) = 0.
            waed(ie) = 0.
          end do

          iwave = nint(getpar('iwave'))
          if( iwave .le. 0 ) icall = -1
          if( iwave .gt. 1 ) icall = -1
          if( icall .le. -1 ) return
          icall = 1

c         --------------------------------------------------
c         Initialize output
c         --------------------------------------------------

          ius = 31
          itmcon = nint(getpar('itmcon'))
          idtcon = nint(getpar('idtcon'))

          call confop(ius,itmcon,idtcon,1,3,'wav')

          write(6,*) 'wave model initialized...'

        endif

c -------------------------------------------------------------------
c normal call
c -------------------------------------------------------------------

c --- get the wind speed and direction

        call c2p(wxnv(1),wynv(1),wis,wid)

c --- get the wind fetch

        call fetch(wid,fet,daf)

c       -------------------------------------------------------------------
c       start loop on elements
c       -------------------------------------------------------------------

        do ie = 1,nel

          icount = 1
c --- get averaged depth along the fetch

          dep = daf(ie)
          depe = depele(ie,+1)
          dep = depele(ie,+1)
10        continue

c --- calculate the wave height, period and direction

          gh = (g*dep)/(wis**2.)
          gx = (g*fet(ie))/(wis**2.)
          hg = dep / (g*wis**2.)
          auxh = ah2*gh**eh1
          auxh1 = ah2*hg**eh1
          auxt = at2*gh**et1
          auxt1 = ah2*gx**eh1

C method of SPM

          waeh(ie) = (tanh(auxh))**eh4
          waeh(ie) = (ah3*gx**eh3) / waeh(ie)
          waeh(ie) = (tanh(waeh(ie)))**eh2
          waeh(ie) = ah1*tanh(auxh)*waeh(ie)
          waeh(ie) = waeh(ie) * wis**2 / g
          
          waep(ie) = (tanh(auxt))**et4
          waep(ie) = (at3*gx**et3) / waep(ie)
          waep(ie) = (tanh(waep(ie)))**et2
          waep(ie) = at1*tanh(auxt)*waep(ie)
          waep(ie) = waep(ie) * wis / g

c          waeh(ie) = 0.283 * tanh(0.530*(gh**(3./4.)))*
c     %            tanh((0.00565*(gx**0.5))/
c     %            (tanh(0.530*(gh**(3./4.)))))*((wis**2)/g)
c
c          waep(ie) = 7.54*tanh(0.833*(gh**(3./8.)))*
c     %            tanh((0.0379*(gx**(1./3.)))/
c     %            (tanh(0.833*(gh**(3./8.)))))*(wis/g)
c

C method of hurdle and stive

c          waeh(ie) = (tanh(auxh))**eh4
c          waeh(ie) = (ah3*gx**eh3) / waeh(ie)
c          waeh(ie) = (tanh(waeh(ie)))**eh2
c          waeh(ie) = ah1*tanh(auxh)*waeh(ie)
c          waeh(ie) = waeh(ie) * wis**2 / g
         
c          waep(ie) = (tanh(auxt1))**et4
c          waep(ie) = (at3*gx**et3) / waep(ie)
c          waep(ie) = (tanh(waep(ie)))**et2
c          waep(ie) = at1*tanh(auxt)*waep(ie)
c          waep(ie) = waep(ie) * wis / g

          waed(ie) = wid

c --- limiting wave height

          hbr = 0.50*depe
          if( waeh(ie).gt.hbr ) then
            if (icount.gt.1) then
              waeh(ie) = hbr
              go to 20
            end if
            dep = depe
            icount = 2
            goto 10
          end if 
20        continue

        end do

        widold = wid

        call make_stress(waeh,waep,z0,tcv,twv,tmv)

c       -------------------------------------------------------------------
c       write of results (file WAV)
c       -------------------------------------------------------------------

        call e2n2d(waeh,waveh,v1v)
        call confil(ius,itmcon,idtcon,31,1,waveh)
        call e2n2d(waep,wavep,v1v)
        call confil(ius,itmcon,idtcon,32,1,wavep)
        call e2n2d(waed,waved,v1v)
        call confil(ius,itmcon,idtcon,33,1,waved)

c        call e2n2d(fet,v2v,v1v)
c        call confil(ius,itmcon,idtcon,33,1,v2v)

c        call e2n2d(tcv,v2v,v1v)
c        call confil(ius,itmcon,idtcon,34,1,v2v)
c        call e2n2d(twv,v2v,v1v)
c        call confil(ius,itmcon,idtcon,35,1,v2v)
c        call e2n2d(tmv,v2v,v1v)
c        call confil(ius,itmcon,idtcon,36,1,v2v)

        end

c**************************************************************

        subroutine fetch(wid,fet,daf)

c This subroutine computes the wind fetch for each element of the
c grid given the wind direction.

        implicit none
  
        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xe,ye		!element point coordinates [m]
        real xf,yf		!far away point coordinates [m]
        real xnew,ynew		!new coordinates [m]
        real d,de		!distance between points [m]
        real wid		!wind direction [degree north]
        real fet(neldim)        !wind fetch length [m]
        real daf(neldim)	!averaged depth along the fetch [m]
        real depele             !element depth function [m]
        real dep		!element depth [m]
        real rad,wdir
        integer ie,iie,ii,ienew,icount,ieo
	save rad,wdir

        d = 500000.
        rad = 45. / atan (1.)
        wdir = wid / rad		!from deg to rad

c --- loop over elements

        do ie = 1,nel
          
          daf(ie) = 0.
          fet(ie) = 0.

c --- get element coordinates

          call baric(ie,xe,ye)

c --- get far away points coordinates

          xf = xe + d*sin(wdir)
          yf = ye + d*cos(wdir)

          iie = ie
          ieo = ie

c --- calculate fetc and averaged depth along the fetch

          icount = 0

1         continue
          call intersect(iie,xe,ye,xf,yf,ienew,xnew,ynew,ieo)
          dep = depele(iie,+1)
          
          de = ((xnew-xe)**2 + (ynew-ye)**2)**0.5
          fet(ie) = fet(ie) + de
          daf(ie) = daf(ie) + dep*de
          ieo = iie
          iie = ienew
          xe = xnew
          ye = ynew
          icount = icount + 1
          if(icount.gt.1000) stop 'error: number of interactions 
     $    exceeds the limit'
          if(ienew.gt.0) go to 1
          daf(ie) = daf(ie)/fet(ie)
          if(ienew.lt.0) then			!open boundary
           fet(ie) = fet(ie) + 50000.
          end if
        end do 
 
        end
           
c******************************************************************

        subroutine intersect(iie,x,y,xf,yf,ien,xn,yn,ieold)

c this routine computes the coordinate of the intersection beetwen the
c line and one of the border line of the element

        implicit none

        integer iie		!element number
        real x,y		!start point cooridnates [m]
        real xf,yf		!far away point coordinates [m]
        integer ien		!next element number
        real xn,yn		!intersection coordinates [m]

        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv
        integer ieltv(3,1)
        common  /ieltv/ieltv

        real x0(3),y0(3)        !element vertices coordinates [m]
        real x3,y3,x4,y4	!node points coordiantes [m]
        real xi,yi		!intersection point coordinate [m]
        integer iint,i,ii,iii
        integer ieold,ie

        integer segsegint	!intersection function

        ien = 0
        call getexy(iie,x0,y0)
 
        do i = 1,3

          ii=mod(i,3)+1
          iii=mod(ii,3)+1

          x3=x0(ii)
          y3=y0(ii)
          x4=x0(iii)
          y4=y0(iii)

2         continue
          iint = segsegint(x,y,xf,yf,x3,y3,x4,y4,xi,yi)

          ie = ieltv(i,iie)

          if(iint.gt.0.and.ie.ne.ieold)then	!intersection
            if(iint.eq.3)then	 		!intersection with node
              x = x + 1.
              y = y + 1.
              go to 2
            else
              xn = xi
              yn = yi
              ien = ie
            end if
          end if

        end do

        end

c******************************************************************

        subroutine make_stress(waeh,waep,z0,tcv,twv,tmv)

c computes stress parameters

        implicit none

        include 'param.h'

        real z0
        real tcv(1)
        real twv(1)
        real tmv(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real pi,karm,rho,g
        real depth,ux,uy,uc2
        real h,p
        real aux,cd
        real omega,zeta,eta
        real uw,a,fw
        real tc,tw,tm

        integer ie

        real waeh(neldim)	!wave height [m]
        real waep(neldim)	!wave period [s]

        real unv(1), vnv(1)
        common /unv/unv, /vnv/vnv

        real depele
        logical is_r_nan

        pi = 3.14159
        karm = 0.4
        rho = 1025.
        g = 9.81

        do ie = 1,nel

          depth = depele(ie,+1)
          ux = unv(ie)/depth
          uy = vnv(ie)/depth
          uc2 = ux*ux + uy*uy
          h = waeh(ie)
          p = waep(ie)

          aux = (z0 + 0.5 * depth) / z0
          aux = log(aux)
          cd = karm/aux
          cd = cd*cd

          tc = rho * cd * uc2

	  if( p .ne. 0. ) then
            omega = 2.*pi/p
            zeta = omega * omega * depth / g
            if( zeta .lt. 1. ) then
              eta = sqrt(zeta) * ( 1. + 0.2 * zeta )
            else
              eta = zeta * ( 1. + 0.2 * exp(2.-2.*zeta) )
            end if
            !k = eta / depth

            uw = pi * h / ( p * sinh(eta) )
            a = uw * p / (2.*pi)
            if( a .gt. 0. ) then
              fw = 1.39 * (z0/a)**0.52
            else
              fw = 0.
            end if

            tw = 0.5 * rho * fw * uw * uw
            tm = tc * ( 1. + 1.2 * ( tw/(tc+tw) )**3.2 )
	  else
	    tw = 0.
	    tm = tc
	  end if

          if( is_r_nan(tw) ) then
            write(6,*) "*** nan in stress..."
            write(6,*) ie,tc,tw,tm
            write(6,*) uw,eta,a,fw,depth,h,p
          end if

          tcv(ie) = tc
          twv(ie) = tw
          tmv(ie) = tm
        end do

        end

c******************************************************************

