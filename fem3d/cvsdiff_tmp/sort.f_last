c
c $Id: sort.f,v 1.3 1999/04/13 08:49:55 georg Exp $
c
c**************************************************************

      SUBROUTINE SORT(N,RA)

c heapsort

      implicit none

      integer n
      real RA(N)
c      integer ra(n)

      integer i,ir,j,l
      real rra
c      integer rra

      if(n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

c******************************************************************

      SUBROUTINE iSORT(N,RA,index)

c heapsort

      implicit none

      integer n
      integer RA(N),index(n)

      integer i,ir,j,l,irra
      integer rra

      do i=1,n
	index(i)=i
      end do

      if(n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
	  irra=index(l)
          RRA=RA(irra)
        ELSE
	  irra=index(ir)
          RRA=RA(irra)
          index(ir)=index(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            index(1)=irra
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(index(j)).LT.RA(index(J+1)))J=J+1
          ENDIF
          IF(RRA.LT.RA(index(j)))THEN
            index(i)=index(j)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        index(i)=irra
      GO TO 10
      END

c**********************************************************

      integer function locate(n,ixx,index,ix)

c locates index for ix in array ixx ( ixx(index(i)) is sorted )

      implicit none

      integer n
      integer ixx(n)
      integer index(n)
      integer ix

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(ix.gt.ixx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      if(ju.le.0) then                    !empty list
        locate=0
      else if(ix.eq.ixx(index(ju))) then  !found
        locate=index(ju)
      else                                !not found
        locate=0
      end if

      return
      end
