
        integer ndim_stab
        parameter (ndim_stab = 50)

        integer nentry
        real rkind(2,ndim_stab)

        common /istab/nentry, /rstab/rkind

        save /istab/, /rstab/

