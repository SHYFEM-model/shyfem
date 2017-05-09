!
! $Id: subtidef.f,v 1.7 2015-09-30 15:38:36 georg Exp $
!
! subroutines for computing tidal potential
!
!********************************************************************
!
! revision log :
!
! 30.09.2015    ccf     converted previous routines to module
! 01.10.2015    ccf     bug fix in hour, now handle negative time
! 21.10.2015	ccf	documentation included
! 29.03.2017	ccf	bug fix in the astronomical arguments chi
!
!********************************************************************
c DOCS  START   S_tidef

c SHYFEM includes an as\-tro\-no\-mi\-cal tidal model which can be
c activated by setting the parameter |rtide| equal 1 in the |para| 
c section.
c
c The model calculates equilibrium tidal potential ($\eta$) and load 
c tides ($\beta$) and uses these to force the free surface. 
c The term $\eta$ in the momentum equations is calculated as a sum
c of the tidal potential of each tidal constituents multiplied by the
c frequency-dependent elasticity factor. The factor $\beta$ accounts 
c for the effect of the load tides, assuming that loading tides are
c in-phase with the oceanic tide. $\beta$ is function of the water depth
c as $\beta=ltidec*H$ with |ltidec| a calibration factor to be set in the
c str |para| section.
c
c The model cosiders the following tidal costituents:
c \begin{itemize}
c \item Semidiurnal species:
c    \begin{itemize}
c    \item M2  principal lunar
c    \item S2  principal solar
c    \item N2  elliptical lunar
c    \item K2  declination luni-solar
c    \end{itemize}
c \item Diurnal species:
c    \begin{itemize}
c    \item K1  declination luni-solar
c    \item O1  principal lunar
c    \item P1  principal solar
c    \item Q1  elliptical lunar
c    \end{itemize}
c \item Long-Period species:
c    \begin{itemize}
c    \item Mf  fortnightly lunar
c    \item Mm  monthly lunar
c    \item Ssa semiannual solar
c    \item MSm S0-semiannual solar
c    \end{itemize}
c \end{itemize}

c DOCS  END
!********************************************************************

!==================================================================
        module tidef
!==================================================================

        implicit none

        integer, private, save  :: nkn_tide = 0

        integer, save           :: rtide	! parameter to the tide model
        real, save           	:: ltidec	! calibration coefficient for load tide
        integer, parameter      :: ntd = 12 	! number of tidal constituents
        real, allocatable, save :: zeqv(:)	! tidal equilibrium

!==================================================================
        contains
!==================================================================

        subroutine tidef_init(nkn)

        integer	:: nkn

        if( nkn == nkn_tide ) return

        if( nkn_tide > 0 ) then
            deallocate(zeqv)
        end if

        nkn_tide = nkn

        if( nkn == 0 ) return

        allocate(zeqv(nkn))

        end subroutine tidef_init

!********************************************************************
! This subroutine computes astronomical arguments chi (in rad)
!********************************************************************

        subroutine get_chi(jd,iyear,chi)

        implicit none

        integer, intent(in)   	:: jd		!julian day
        integer, intent(in) 	:: iyear	!year
        real,intent(out),dimension(ntd) :: chi	!astronomical arguments [rad]

        real, parameter    :: degrad = 57.2957795131  ! 360/(2*pi)
        real, parameter    :: h01 = 279.69668
        real, parameter    :: h02 = 36000.768930485
        real, parameter    :: h03 = 3.03E-04
        real, parameter    :: s01 = 270.434358
        real, parameter    :: s02 = 481267.88314137
        real, parameter    :: s03 = -0.001133
        real, parameter    :: s04 = 1.9E-06
        real, parameter    :: p01 = 334.329653
        real, parameter    :: p02 = 4069.0340329575
        real, parameter    :: p03 = -0.010325
        real, parameter    :: p04 = -1.2E-05
        real, parameter    :: t1  = 0.74996579132
        real, parameter    :: t2  = 2.7378508846E-05
        integer  	   :: i
        real		   :: d,t,h0,s0,p0

	if (iyear < 1975 ) then
            write(6,*) 'astronomical arguments chi could be computed'
            write(6,*) 'only for year >= 1975'
            write(6,*) 'year = ', iyear
            stop 'error stop get_chi: year < 1975'
	end if

        d  = jd + 365.0*(iyear - 1975.0) + int((iyear-1975.0)/4.0)
        t  = t1 + t2*d
        h0 = h01 + h02*t + h03*(t**2.0)
        s0 = s01 + s02*t + s03*(t**2.0) + s04*(t**3.0)
        p0 = p01 + p02*t + p03*(t**2.0) + p04*(t**3.0)

        chi(1)  = 2*h0 - 2*s0           !M2
        chi(2)  = 0.0                   !S2
        chi(3)  = 2*h0 - 3*s0 + p0      !N2
        chi(4)  = 2*h0                  !K2

        chi(5)  = h0        + 90.0      !K1
        chi(6)  = h0 - 2*s0 - 90.0      !O1
        chi(7)  = -h0       - 90.0      !P1
        chi(8)  = h0 - 3*s0 + p0 -90.0  !Q1

        chi(9)  =      2*s0             !Mf
        chi(10) =       s0 - p0         !Mm
        chi(11) = 2*h0                  !Ssa
        chi(12) =-2*h0 + s0 +p0	        !MSm

! ----- Convert to Radians
        chi = chi/degrad

        end subroutine get_chi

!*********************************************************************
! This subroutine computes the equilibrium tide as a sum of the single
! costituents

        subroutine equilibrium_tide(lat,lon,hour,chi,eqtide)

        implicit none

        real, intent(in) 	:: lat,lon	!Latitude and longitude [degree]
        real, intent(in) 	:: hour		!Universal standard time [days]
        real, intent(in) 	:: chi(ntd)	!Astronomic arguments
        real, intent(out)	:: eqtide	!Equilibrium tide [m]

        real, parameter, dimension(ntd) :: ome =      !Tidal constituent frequency [2Pi/s]
     +         (/1.40519e-4,1.45444e-4,1.37880e-4,1.458426e-4,
     +           0.72921e-4,0.67598e-4,0.72523e-4,0.649584e-4, 
     +           0.053234e-04,0.026392e-04,0.003982e-04,0.049252e-04/)
        real, parameter, dimension(ntd) :: amp =      !Tidal constituent amplitude [m]
     +         (/0.242334,0.112841,0.046398,0.030704, 
     +           0.141465,0.100514,0.046843,0.019256,
     +           0.041742,0.022026,0.019446,0.004239/)
        real, parameter, dimension(ntd) :: loven =    !Elasticy factor
     +         (/0.6930,0.6930,0.6930,0.6930,
     +           0.7364,0.6950,0.7059,0.6946, 
     +           0.6920,0.6920,0.6920,0.6920/)
        real, parameter, dimension(ntd) :: mask_t =   !Mask: =1, tide included
     +         (/1., 1., 1., 1., 
     +           1., 1., 1., 1., 
     +           1., 1., 0., 0./)

        real, parameter 	:: pi2 = 6.283185307177960
        real, parameter         :: rad = pi2/360.0
        real  			:: colat,twolon,time
        real  			:: sins,sind,sinl
        real			:: llon
        integer  		:: n,m,l

        eqtide = 0.
        time   = hour * 86400.
        llon = lon
        if (llon < 0.0) llon = 360. + llon
        llon    = rad * llon
        twolon = 2.0*llon
        colat  = rad * (90. - lat)
        sins   = sin(colat)**2.0
        sind   = sin(2.0*colat)
        sinl   = 3.0*sins - 2.0

        do n = 1,4
! ------- Semidiurnal species -------
            eqtide = eqtide + mask_t(n) * loven(n) * 
     +               amp(n) * sins *cos(ome(n)*time + chi(n) + twolon)
! ------- Diurnal species -------
            m = n+4
            eqtide = eqtide + mask_t(m) * loven(m) *
     +               amp(m) * sind *cos(ome(m)*time + chi(m) + llon)
! ------- Long-period species -------
            l = m+4
            eqtide = eqtide + mask_t(l) * loven(l) * 
     +               amp(l) * sinl *cos(ome(l)*time + chi(l))
        enddo

        return

        end subroutine equilibrium_tide

!==================================================================
        end module tidef
!==================================================================

        subroutine tidefini

        use tidef
        use basin
        use coordinates, only : iproj

        implicit none

        integer  	 :: isphe      !if = 1  coordinates are in spherical system
        double precision :: dgetpar

        zeqv = 0.

        rtide = dgetpar('rtide')

        if( rtide <= 0. ) return

        ltidec = dgetpar('ltidec')

        call get_coords_ev(isphe)

        if( isphe <= 0 .AND. iproj <= 0 ) then
            write(6,*) 'isphe,iproj: ',isphe,iproj
            write(6,*) 'for tidal potential either '
            write(6,*) 'isphe or iproj must be set'
            stop 'error stop tidefini: no lat/lon coordinates'
        end if

        write(6,*) 'Tidal potential active'

        end subroutine tidefini

!*********************************************************************

        subroutine tideforce(it)

        use tidef
        use coordinates
        use mod_depth
        use mod_hydro
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer, intent(in)	:: it

        real  			:: lat,lon,eqt
        integer  		:: iy,im,id,ih,imn,isec
        integer			:: k,jd
        real			:: hour
        real, dimension(ntd) 	:: chi    !astronomical arguments [rad]
        real			:: loadb  !loading tide factor [0.054,0.047,= ltidec*depth]

!--------------------------------------------------------
!------ compute tidal potential? ------------------------
!--------------------------------------------------------

        if( rtide <= 0 ) return
        	
!--------------------------------------------------------
!------ computes julian day - 1 -------------------------
!--------------------------------------------------------

        call dts2dt(it,iy,im,id,ih,imn,isec)
        call date2j(iy,im,id,jd)
	hour = (ih*3600. + imn*60. + isec) / 86400.

!--------------------------------------------------------
!------ computes astronomical arguments -----------------
!--------------------------------------------------------

        call get_chi(jd,iy,chi)

!--------------------------------------------------------
!------ computes eq. tide and account for load tide -----
!--------------------------------------------------------

        do k = 1,nkn
            lat = ygeov(k)
            lon = xgeov(k)
            call equilibrium_tide(lat,lon,hour,chi,eqt)
            loadb = ltidec * (hkv(k) + zov(k))
            zeqv(k) = eqt + loadb*zov(k)
        end do

!--------------------------------------------------------
!------ end of routine ----------------------------------
!--------------------------------------------------------

        end subroutine tideforce

!*********************************************************************
