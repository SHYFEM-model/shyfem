
        integer proj
        common /coords1/ proj
        save /coords1/

        double precision aearth,bearth,flat,rflat,e,es
        common /proj_param01/ aearth,bearth,flat,rflat,e,es
        save /proj_param01/

        double precision k0
        common /proj_param02/ k0
        save /proj_param02/

        double precision zero,one,two,four,half
        common /proj_param11/ zero,one,two,four,half
        save /proj_param11/

        double precision pi,half_pi,quarter_pi,rad,rrad
        common /proj_param12/ pi,half_pi,quarter_pi,rad,rrad
        save /proj_param12/

        double precision eps,tol
        common /proj_param13/ eps,tol
        save /proj_param13/

        logical debug
        common /proj_param21/ debug
        save /proj_param21/

!-------------------------------------------------------------------

!        integer fuso                    !1 or 2
!        common /gb_param1/ fuso
!        double precision lambda0,x0,xtrans0,ytrans0
!        common /gb_param2/ lambda0,x0,xtrans0,ytrans0
!        save /gb_param1/,/gb_param2/

!-------------------------------------------------------------------

!        integer utm_zone                !1-60
!        common /utm_param1/ utm_zone
!        double precision lambda0,xtrans0,ytrans0
!        common /utm_param2/ lambda0,xtrans0,ytrans0
!        double precision ep2,dn,es_4,es_6,z1,z2,z3,z4,j1,j2,j3,j4
!        common /utm_param3/ ep2,dn,es_4,es_6,z1,z2,z3,z4,j1,j2,j3,j4
!        save /utm_param1/,/utm_param2/,/utm_param3/

!-------------------------------------------------------------------

!        double precision lon_0,lat_0,phi_0,xfact,yfact
!        common /cpp_param1/ lon_0,lat_0,phi_0,xfact,yfact
!        save /cpp_param1/

!-------------------------------------------------------------------

