
        integer utm_zone                !1-60
        common /utm_param1/ utm_zone
        double precision lambda0,xtrans0,ytrans0
        common /utm_param2/ lambda0,xtrans0,ytrans0
        double precision ep2,dn,es_4,es_6,z1,z2,z3,z4,j1,j2,j3,j4
        common /utm_param3/ ep2,dn,es_4,es_6,z1,z2,z3,z4,j1,j2,j3,j4
        save /utm_param1/,/utm_param2/,/utm_param3/

