
        integer ielfix(0:3,neldim)
        common /ielfix/ielfix
        save /ielfix/

        integer tnudgev(neldim)
        common /tnudgev/tnudgev
        save /tnudgev/

        real ubound(nlvdim,nkndim)
        real vbound(nlvdim,nkndim)
        common /ubound/ubound
        common /vbound/vbound
	save /ubound/,/vbound/
