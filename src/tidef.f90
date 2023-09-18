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
! DOCS  START   S_tidef

! SHYFEM includes an as\-tro\-no\-mi\-cal tidal model which can be
! activated by setting the parameter |rtide| equal 1 in the |para| 
! section.
!
! The model calculates equilibrium tidal potential ($\eta$) and load 
! tides ($\beta$) and uses these to force the free surface. 
! The term $\eta$ in the momentum equations is calculated as a sum
! of the tidal potential of each tidal constituents multiplied by the
! frequency-dependent elasticity factor. The factor $\beta$ accounts 
! for the effect of the load tides, assuming that loading tides are
! in-phase with the oceanic tide. $\beta$ is function of the water depth
! as $\beta=ltidec*H$ with |ltidec| a calibration factor to be set in the
! str |para| section.
!
! The model cosiders the following tidal costituents:
! \begin{itemize}
! \item Semidiurnal species:
!    \begin{itemize}
!    \item M2  principal lunar
!    \item S2  principal solar
!    \item N2  elliptical lunar
!    \item K2  declination luni-solar
!    \end{itemize}
! \item Diurnal species:
!    \begin{itemize}
!    \item K1  declination luni-solar
!    \item O1  principal lunar
!    \item P1  principal solar
!    \item Q1  elliptical lunar
!    \end{itemize}
! \item Long-Period species:
!    \begin{itemize}
!    \item Mf  fortnightly lunar
!    \item Mm  monthly lunar
!    \item Ssa semiannual solar
!    \item MSm S0-semiannual solar
!    \end{itemize}
! \end{itemize}

! DOCS  END
!********************************************************************

!==================================================================
        module tidef
!==================================================================

        implicit none

        integer, private, save  :: nkn_tide = 0

        integer, save           :: rtide	! parameter to the tide model
        double precision, save           	:: ltidec	! calibration coefficient for load tide
        integer, parameter      :: ntd = 12 	! number of tidal constituents
        double precision, allocatable, save :: zeqv(:)	! tidal equilibrium

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
        double precision,intent(out),dimension(ntd) :: chi	!astronomical arguments [rad]

        double precision, parameter    :: degrad = 57.2957795131  ! 360/(2*pi)
        double precision, parameter    :: h01 = 279.69668
        double precision, parameter    :: h02 = 36000.768930485
        double precision, parameter    :: h03 = 3.03E-04
        double precision, parameter    :: s01 = 270.434358
        double precision, parameter    :: s02 = 481267.88314137
        double precision, parameter    :: s03 = -0.001133
        double precision, parameter    :: s04 = 1.9E-06
        double precision, parameter    :: p01 = 334.329653
        double precision, parameter    :: p02 = 4069.0340329575
        double precision, parameter    :: p03 = -0.010325
        double precision, parameter    :: p04 = -1.2E-05
        double precision, parameter    :: t1  = 0.74996579132
        double precision, parameter    :: t2  = 2.7378508846E-05
        integer  	   :: i
        double precision		   :: d,t,h0,s0,p0

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

        double precision, intent(in) 	:: lat,lon	!Latitude and longitude [degree]
        double precision, intent(in) 	:: hour		!Universal standard time [days]
        double precision, intent(in) 	:: chi(ntd)	!Astronomic arguments
        double precision, intent(out)	:: eqtide	!Equilibrium tide [m]
        !Tidal constituent frequency [2Pi/s]
        double precision, parameter, dimension(ntd) :: ome = (/1.40519e-4,1.45444e-4,1.37880e-4,1.458426e-4,        &
    &   0.72921e-4,0.67598e-4,0.72523e-4,0.649584e-4,0.053234e-04,0.026392e-04,0.003982e-04,0.049252e-04/)
        !Tidal constituent amplitude [m]
        double precision, parameter, dimension(ntd) :: amp = (/0.242334,0.112841,0.046398,0.030704,0.141465,        &
    &                                   0.100514,0.046843,0.019256,0.041742,0.022026,0.019446,0.004239/)
        !Elasticy factor
        double precision, parameter, dimension(ntd) :: loven = (/0.6930,0.6930,0.6930,0.6930,0.7364,0.6950,         &
    &                                   0.7059,0.6946,0.6920,0.6920,0.6920,0.6920/)
        !Mask: =1, tide included
        double precision, parameter, dimension(ntd) :: mask_t = (/1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0./)

        double precision, parameter 	:: pi2 = 6.283185307177960
        double precision, parameter         :: rad = pi2/360.0
        double precision  			:: colat,twolon,time
        double precision  			:: sins,sind,sinl
        double precision			:: llon
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
            eqtide = eqtide + mask_t(n) * loven(n) * amp(n) * sins *cos(ome(n)*time + chi(n) + twolon)
! ------- Diurnal species -------
            m = n+4
            eqtide = eqtide + mask_t(m) * loven(m) * amp(m) * sind *cos(ome(m)*time + chi(m) + llon)
! ------- Long-period species -------
            l = m+4
            eqtide = eqtide + mask_t(l) * loven(l) * amp(l) * sinl *cos(ome(l)*time + chi(l))
        enddo

        return

        end subroutine equilibrium_tide

!==================================================================

        subroutine tidefini

        use basin
        use evgeom
        use coordinates, only : iproj
        use para

        implicit none

        integer  	 :: isphe      !if = 1  coordinates are in spherical system

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

        use coordinates
        use depth
        use hydro_admin
        use dts
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer, intent(in)	:: it

        double precision  			:: lat,lon,eqt
        integer  		:: iy,im,id,ih,imn,isec
        integer			:: k,jd
        double precision			:: hour
        double precision, dimension(ntd) 	:: chi    !astronomical arguments [rad]
        double precision			:: loadb  !loading tide factor [0.054,0.047,= ltidec*depth]

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
!==================================================================
        end module tidef
!==================================================================
