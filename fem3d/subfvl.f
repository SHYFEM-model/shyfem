c
c routine to write volume
c
c revision log :
c
c 28.04.2010    ggu     written from scratch
c
c******************************************************************

	subroutine wrfvla

c write of finite volume data

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(1)
        common /ilhkv/ilhkv

        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real areakv(nlvdim,nkndim)
        common /areakv/areakv
        real saux1(nlvdim,nkndim)
        common /saux1/saux1

	integer k,l,lmax

	real getpar

        integer iu,id,itmcon,idtcon
        save iu,id,itmcon,idtcon

        integer icall
        save icall
        data icall /0/

c start of code

        if( icall .eq. -1 ) return

c initialization

        if( icall .eq. 0 ) then

          iu = 0
          id = 66       			!for finite volume
          itmcon = nint(getpar('itmcon'))
          idtcon = nint(getpar('idtcon'))
	  if( idtcon .le. 0 ) icall = -1
	  if( icall .le. -1 ) return

          call confop(iu,itmcon,idtcon,nlv,1,'fvl')

        end if

c normal call

        icall = icall + 1

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    saux1(l,k) = areakv(l,k) * hdknv(l,k)
	  end do
	end do

	call confil(iu,itmcon,idtcon,id,nlvdi,saux1)

	end

c******************************************************************

