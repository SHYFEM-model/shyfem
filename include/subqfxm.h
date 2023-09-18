! may be already defined elsewhere
	double precision, parameter  :: rhow = 1026.
        double precision, parameter  :: pi = 3.14159265358979323846
        double precision, parameter  :: deg2rad = pi/180.
        double precision, parameter  :: rad2deg = 180./pi
        double precision, parameter  :: visw = 1.e-6      !molecular kinematic viscosity	
        double precision, parameter  :: diff = 1.4e-7     !molecular thermalconductivity

! used in subqfxm4.f and subqfxm5.f
        !double precision, parameter  :: kelv = 273.16
        double precision, parameter  :: kelv = 273.15
        double precision, parameter  :: emiss = 0.97
        double precision, parameter  :: bolz = 5.67e-8
        double precision, parameter  :: cpa  = 1004.67    !Spcf heat cap. dry air [J/kg/K]
        !double precision, parameter  :: cpa  = 1005.00    !Spcf heat cap. dry air [J/kg/K]
        double precision, parameter  :: cpw  = 3985.      !Spcf heat cap. seawater [J/kg/K]
        double precision, parameter  :: rgas = 287.1      !Gas constant dry air [J/kg/K]
        !double precision, parameter  :: rgas = 287.0      !Gas constant dry air [J/kg/K]
        double precision, parameter  :: grav = 9.81	      !gravitational accel. [m/s2]
        double precision, parameter  :: von  = 0.41	      !von Karman
        double precision, parameter  :: const06 = 0.62198
        !double precision, parameter  :: const06 = 0.622

        double precision, parameter  :: a1=6.107799961
        double precision, parameter  :: a2=4.436518521e-1
        double precision, parameter  :: a3=1.428945805e-2
        double precision, parameter  :: a4=2.650648471e-4
        double precision, parameter  :: a5=3.031240396e-6
        double precision, parameter  :: a6=2.034080948e-8
        double precision, parameter  :: a7=6.136820929e-11
        double precision, parameter  :: eps=1.0e-12
