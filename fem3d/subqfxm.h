! may be already defined elsewhere
	real, parameter  :: rhow = 1026.
        real, parameter  :: pi = 3.14159265358979323846
        real, parameter  :: deg2rad = pi/180.
        real, parameter  :: rad2deg = 180./pi
        real, parameter  :: visw = 1.e-6      !molecular kinematic viscosity	
        real, parameter  :: diff = 1.4e-7     !molecular thermalconductivity

! used in subqfxm4.f and subqfxm5.f
        !real, parameter  :: kelv = 273.16
        real, parameter  :: kelv = 273.15
        real, parameter  :: emiss = 0.97
        real, parameter  :: bolz = 5.67e-8
        real, parameter  :: cpa  = 1004.67    !Spcf heat cap. dry air [J/kg/K]
        !real, parameter  :: cpa  = 1005.00    !Spcf heat cap. dry air [J/kg/K]
        real, parameter  :: cpw  = 3985.      !Spcf heat cap. seawater [J/kg/K]
        real, parameter  :: rgas = 287.1      !Gas constant dry air [J/kg/K]
        !real, parameter  :: rgas = 287.0      !Gas constant dry air [J/kg/K]
        real, parameter  :: grav = 9.81	      !gravitational accel. [m/s2]
        real, parameter  :: von  = 0.41	      !von Karman
        real, parameter  :: const06 = 0.62198
        !real, parameter  :: const06 = 0.622

        real, parameter  :: a1=6.107799961
        real, parameter  :: a2=4.436518521e-1
        real, parameter  :: a3=1.428945805e-2
        real, parameter  :: a4=2.650648471e-4
        real, parameter  :: a5=3.031240396e-6
        real, parameter  :: a6=2.034080948e-8
        real, parameter  :: a7=6.136820929e-11
        real, parameter  :: eps=1.0e-12
