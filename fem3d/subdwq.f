c
c $Id: subdwq.f,v 1.15 2009-04-07 10:43:57 georg Exp $
c
c delwaq routines
c
c contents :
c
c subroutine delwaq				delwaq interface
c subroutine finvul(it,nagr,unv,vnv)	computes volumes, areas ...
c SUBROUTINE DELPNT ( ... )
c SUBROUTINE DELKNT ( ... )
c SUBROUTINE DELGET ( ... )
c
c revision log :
c
c 27.03.1998	ggu	integrated changes from Technital/Delft
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 29.04.1998    ggu     uses module for semi-implicit time-step
c 04.05.1998    ggu     keyword access introduced
c 13.05.1998    ggu     irv not declared anymore
c 20.05.1998    ggu     open files through ifileo; save nvol, nflow, narea
c
c*****************************************************************
c
        subroutine delwaq
c
c delwaq interface
c
        implicit none
c
c common
        integer itanf,itend,idt,nits,niter,it
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real uov(1),vov(1),unv(1),vnv(1)
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
c local
        integer itsdwq,itedwq,idtdwq,nagdwq
c function
        integer iround
        real getpar
c save
        save itsdwq,itedwq,idtdwq,nagdwq
c
c first call
c
        if(niter.eq.1) then
          itsdwq=iround(getpar('itsdwq'))
          itedwq=iround(getpar('itedwq'))
          idtdwq=iround(getpar('idtdwq'))
          nagdwq = idtdwq/idt
c
          if(idtdwq.gt.0) then
            if( mod(idtdwq,idt) .ne. 0 ) goto 99
            if( mod(itsdwq,idtdwq) .ne. 0 ) goto 99
            if( mod(itedwq,idtdwq) .ne. 0 ) goto 99
            if( mod(itedwq-itsdwq,idtdwq) .ne. 0 ) goto 99
            if( mod(itsdwq-itanf,idt) .ne. 0 ) goto 99
          end if
        end if
c
c every call
c
        if(idtdwq.gt.0.and.it.ge.itsdwq.and.it.le.itedwq) then
          call finvul(it,nagdwq,unv,vnv)
        end if
c
c last call --> close files
c
c        if(idtdwq.gt.0.and.it.eq.itedwq) then
c          endfile(45)
c          endfile(46)
c          endfile(47)
c          close(45)
c          close(46)
c          close(47)
c        end if
c
	return
   99   continue
        write(6,*) 'idt,itanf,itend : ',idt,itanf,itend
        write(6,*) 'idtdwq,itsdwq,itedwq : ',idtdwq,itsdwq,itedwq
        write(6,*) 'nagdwq : ',nagdwq
        stop 'error stop delwaq : error in DWQ parameters'
	end
c
c*******************************************************************
c
        subroutine finvul(it,nagr,unv,vnv)
c
c computes volumes, areas and flows of finite volumes
c
c written 10.03.92 by ggu   from scratch
c revised 12.03.92 by ggu   new output
c revised 30.04.92 by ggu   aggregation of time steps & areas
c revised ...07.92 by ggu   $$lump  - lumping of matrix
c revised 31.08.92 by ggu   $$impli - implicit time step
c
c
c ------- V(n) -------------- V(n+1) ---------------- V(n+2) ---------
c   f-(n)      f+(n)  f-(n+1)        f+(n+1)  f-(n+2)        f+(n+2)
c ------- f(n) -------------- f(n+1) ---------------- f(n+2) ---------
c
c    V(n+1) = V(n) + dt * ( alf*f(n+1) + (1-alf)*f(n) )    (1a)
c
c where alf is the weighting of the time step
c ( 0.5 --> semi-implicit    1.0 --> fully implicit )
c
c with
c      f+(n) = (1-alf)*f(n)   f-(n) = alf*f(n)
c and
c      f(n) = f-(n) + f+(n)
c formula (1a) reads
c
c    V(n+1) = V(n) + dt * ( f-(n+1) + f+(n) )              (1b)
c
c the total flux between time step n and n+m is
c
c      F = f+(n) + ( f(n+1) + ... + f(n+m-1) ) + f-(n+m)   (2)
c
c program variables : (last written time step is n, next time step
c                       to be written is n+m, actual time step n+i
c                       with 1 <= i <= m )
c
c vol1      contains V(n)
c vol2      contains V(n+i)
c flow1     contains accumulated flow between n and n+i
c             f+(n) + ( f(n+1) +...+ f(n+i-1) ) + f-(n+i)
c flow2     contains f+(n+i)
c
c at time step n+m variables vol1 and flow1 are written, e.g. the volume
c at time step n and the flow between n and n+m
c
	implicit none
c
c arguments
        integer it,nagr
	real unv(1),vnv(1)
c       real z(1),u(1),v(1)
c parameter
        include 'param.h'
	real drittl
	parameter (drittl=1./3.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(1),ygv(1)
        integer nen3v(3,1), iwegv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
        common /iwegv/iwegv
c global...
        integer ipnt(nkndim)
        integer npoint(3,neldim)
        real surf(nkndim)
        real vol1(nkndim),vol2(nkndim)
        real flow1(nlkdim),flow2(nlkdim)
        real area1(nlkdim),area2(nlkdim)
	integer link(ngrdim,nkndim),ilink(nkndim)
c local
        character*80 name
        real flowh,area,areah,rnagr
	real zm,um,vm,hm,xs,ys,hv
	real hmzm
	real x(3),y(3)
        real anew,aold
        real azpar
        integer it1
        integer i,i1,i2,k,k1,k2,n,ie,j,ip
        integer noflow,noseg,nli
	integer ibase,ii
	integer icall
        logical bwrite
	integer nvol,nflow,narea
	integer nout
c function
        integer ipext, ifileo
	integer itybnd, nkbnds, kbnds
	real areaele,depfvol,depele
c save
	save icall
        save rnagr
        save vol2,flow2,area2
        save vol1,flow1,area1
        save ipnt,npoint
        save link,ilink
        save noflow,noseg,nli,it1
	save nvol,nflow,narea
c data
	data icall /0/
c
c set switch for DWQ write
c
        bwrite = mod(icall,nagr).eq.0
c
c dimensions ?
c
        if(ngr.gt.ngrdim) goto 97
	if(nkn.gt.nkndim) goto 96
        if(nel.gt.neldim) goto 98
c
c weighting of time steps
c
	call getaz(azpar)
        anew=azpar
        aold=1.-anew
c
c first initialization of arrays and parameters
c
        if(icall.eq.0) then
c
          rnagr = 1./nagr
c
		do i=1,nkn
		  do n=1,ngr
		    link(n,i)=0
		  end do
		  ilink(i)=0
		  surf(i)=0.
            vol1(i)=0.
		end do

          do i=1,nlkdim
            flow1(i)=0.
            flow2(i)=0.
            area1(i)=0.
            area2(i)=0.
          end do
c
		do ie=1,nel
		  call elebase(ie,n,ibase)
		  area = areaele(ie)/n
		  do i1=1,n
		    i2=mod(i1,3)+1
		    k1=nen3v(i1,ie)
		    k2=nen3v(i2,ie)
c
c             compute surface and static volume of finite volume
c
		    surf(k1)=surf(k1)+area
		    hv = depfvol(ie,ii,0)
                    vol1(k1)=vol1(k1)+area*hv   !$$lump
c
c             set up link index
c
c             in ilink(i) is the number of links of volume i
c             in link(j,i) are the links of volume i
c
		    n=ilink(k1)
		    do j=1,n
                if(link(j,k1).eq.k2) goto 5
		    end do
		    n=n+1
		    if(n.gt.ngr) goto 99
		    link(n,k1)=k2
		    ilink(k1)=n
    5         continue
c
		    n=ilink(k2)
		    do j=1,n
                if(link(j,k2).eq.k1) goto 6
		    end do
		    n=n+1
		    if(n.gt.ngr) goto 99
		    link(n,k2)=k1
		    ilink(k2)=n
    6         continue
		  end do
		end do
c
c negative surface for open boundary nodes
c
c		do i=1,nkn
c		  do j=1,nrb
c		    if(irv(j).eq.i) then
c                surf(i)=-surf(i)
c		    end if
c		  end do
c		end do
c
c compute number of links
c
		nli=0
		do i=1,nkn
		  nli=nli+ilink(i)
		end do
		nli=nli/2
c
          if(nli.gt.nlkdim) goto 95
c
c write file containing geometry and static values of FEM grid
c
          call getfnm('femdwq',name)
	  nout = ifileo(0,name,'binary','new')
          write(nout) nkn,ngr,nel,nli
          write(nout) ((nen3v(i,ie),i=1,3),ie=1,nel)
          write(nout) (i,ilink(i),(link(j,i),j=1,ilink(i)),i=1,nkn)
          write(nout) (i,xgv(i),ygv(i),surf(i),vol1(i),i=1,nkn)
          write(nout) (ipext(i),i=1,nkn)
c
c  open boundary nodes
c
          write(nout) nrb
          do i=1,nbc
	     write(nout) itybnd(i),nkbnds(i)
	     write(nout) (kbnds(i,j),j=1,nkbnds(i))
c            write(nout) iround(bnd(2,i)),
c    *                 iround(bnd(4,i))-iround(bnd(3,i))+1
c            write(nout) (irv(j),j=iround(bnd(3,i)),iround(bnd(4,i)))
          end do
          close(nout)
c
c read file that contains pointers to aggregate structures
c
          call getfnm('pntfem',name)
	  nout = ifileo(0,name,'formatted','unknown')
          READ ( nout , * , END=30 )   noseg , ( IPNT(K) , K=1,NKN )
          goto 45
   30     noseg = nkn
          do i=1,nkn
             ipnt(i) = i
	  end do
   45     continue
          close( nout )
C
c make pointer into array that contains flows
c
          CALL DELPNT ( nen3v   , NEL    , LINK   , ILINK
     *                                   , ngrdim , IPNT
     *                                   , NKN    , NPOINT , NOFLOW )
c
C write the statistics of the coupling
C
	  nout = ifileo(0,'statists.dat','formatted','new')
          WRITE  ( nout , '(A,I6)' ) ' Nr of finite elements: ',NEL
          WRITE  ( nout , '(A,I6)' ) ' Nr of nodes          : ',NKN
          WRITE  ( nout , '(A,I6)' ) ' Nr of GRs (?)        : ',NGR
          WRITE  ( nout , '(A,I6)' ) ' Nr of links          : ',NLI
          WRITE  ( nout , '(A,I6)' ) ' Nr of flows          : ',NOFLOW
          CLOSE  ( nout )
c
c open data files
c
          call getfnm('voldwq',name)
	  nvol = ifileo(0,name,'binary','new')

          call getfnm('flowdwq',name)
	  nflow = ifileo(0,name,'binary','new')

          call getfnm('areadwq',name)
	  narea = ifileo(0,name,'binary','new')
c
        end if
c
c normal program code
c
c initialization of arrays that are to be used at every call
c
        do i=1,nli
c$$impli          flow1(i) = flow1(i)+flow2(i)
          flow1(i) = flow1(i)+aold*flow2(i)
          area1(i) = area1(i)+area2(i)
          flow2(i) = 0.
          area2(i) = 0.
        end do
        do i=1,nkn
          vol2 (i) = 0.
        end do
c
c
c loop over elements
c
	do ie=1,nel
c
	  um=unv(ie)
	  vm=vnv(ie)
	  call baric(ie,xs,ys)
c
	  hmzm = depele(ie,+1)
	  call elebase(ie,n,ibase)
	  area = areaele(ie)/n

            do ii=1,n
              k=nen3v(ii,ie)
              ip=ipnt(k)
              if(ip.gt.0) then
c this formula is used since some elements may be lumped together into
c ... a mega-area and the contributions must be accumulated into
c ... one variable
		hv = depfvol(ie,ii,0)
                vol2(ip)=vol2(ip) + area * hv
              end if
            end do
c
	  do i=1,3
	    k=nen3v(i,ie)
	    i1=mod(i,3)+1
	    k1=nen3v(i1,ie)
	    i2=mod(i+1,3)+1
	    k2=nen3v(i2,ie)
c
c in flowh is now the flow from node i to i1 of element ie
c
c now it m**3/s, not m**3/dt
c
            flowh=-0.5*( um*(ys-y(i2)) + vm*(x(i2)-xs) )
            areah=0.5*sqrt( (ys-y(i2))**2 + (x(i2)-xs)**2 )
     +                          *hmzm
c           this is a crude aproximation for the cross area terms
c           ...if needed it has to be changed
c           corrected 29.4.92 --> not double area anymore
c
c insert this flow in matrix
c
            ip=npoint(i,ie)
c
            if(ip.gt.0) then
              flow1(ip)  = flow1(ip)  + anew*flowh
              flow2(ip)  = flow2(ip)  + flowh
              area2(ip)  = area2(ip)  + areah
            else if(ip.lt.0) then
              flow1(-ip) = flow1(-ip) - anew*flowh
              flow2(-ip) = flow2(-ip) - flowh
              area2(-ip) = area2(-ip) + areah
            end if
	  end do
   77     continue
	end do
c
c end of loop over elements
c
c       output results
c
        if(bwrite.and.icall.gt.0) then
          write(nvol) it1,(vol1 (i),i=1,noseg)
          write(nflow) it1,(rnagr*flow1(i),i=1,noflow)   !$$impli
          write(narea) it1,(rnagr*area1(i),i=1,noflow)   !$$impli
        end if
        if (bwrite) then
          do i=1,noseg
            vol1(i)=vol2(i)
          end do
          do i=1,noflow
            flow1(i) = 0.0
            area1(i) = 0.0
          end do
          it1 = it
        end if
c
c update counter
c
        icall = icall + 1
c
	return
   95   stop 'error stop finvol : dimension nlkdim'
   96   stop 'error stop finvol : dimension nkndim'
   97   stop 'error stop finvol : dimension ngrdim'
   98   stop 'error stop finvol : dimension neldim'
   99   continue
	write(6,*) n,ngr,k1,k2,ie
	stop 'error stop finvol : dimension ngr'
        end
c
c*********************************************************************
c
      SUBROUTINE DELPNT ( NODNUM , NFELT  , LINNUM , NOLINK , ngrdim  ,
     *                              IPNT ,  NONODE , NPOINT , NOFLOW )
C
C     Function: Derives from an array with 3 node numbers per finite
C               element and an array with links per node, a pointer
C               table to the DELWAQ "area" and "flow" file structures.
C
C               The pointer table consists of 3 integers per finite
C               element. The first integer is the array element to save
C               the flow from node number 1 to node number 2 in, the
C               second for the flow from 2 to 3 and the third from 3 to
C
C               The flows and areas have to be summed in an initially
C               clean array. If the pointer value is negative, the flow
C               has to be summed with negative sign on the indicated
C               location.
C
C               The routine also produces the so-called DELWAQ "FROM-TO
C               file and the number of exchanges
C$1-2lines added
C               Double pointering has been added, using IPNT from the
C               SHYFEM nodes to the DELWAQ computational elements.
C
C     Creation date    :  2-04-1992 by L. Postma ( Delft Hydraulics )
C$1-3lines changed
C     Modification  $1 : 28-04-1992 by L. Postma ( Delft Hydraulics )
C                                   subject: Double pointering active
C                                            DELWAQ elts and links
C                          -  -     by
C                                   subject:
C
C     FUNCTIONS CALLED : None
C
C     ROUTINES  CALLED : DELKNT Computes the pointer values
C
C     FILES USED       : 1
C           UNIT NAME           FUNCTION DESCRIPTION
C            91  "dhven991.dat" IN/OUT   Temporal save space
C            92  "fromto.pnt"   OUTPUT   DELWAQ FROM-TO table
C            93  "shyfem.pnt"   OUTPUT   SHYFEM pointering
C
C
C     PARAMS  KIND      DIMENSION FUNCTION DESCRIPTION
C     ------  --------- --------- -------- ----------------------------
C     NODNUM  INTEGER*4  3,NFELT  INPUT    Node nrs for each finite elt
C     NFELT   INTEGER*4     1     INPUT    Number of finite elements
C     LINNUM  INTEGER*4 NGRDIM,NONODE INPUT    Numbers of the links per
C     NOLINK  INTEGER*4   NONODE  INPUT    Number of links per node
C       NGRDIM ...
C$1-1line added
C     IPNT    INTEGER*4   NONODE  INPUT    Pointer to DELWAQ elements
C     NONODE  INTEGER*4     1     INPUT    Number of nodes
C     NPOINT  INTEGER*4  3,NFELT  OUTPUT   Pointers for each finite elt
C     NOFLOW  INTEGER*4     1     OUTPUT   Number of flows identified
C
C
      INTEGER*4  NODNUM(3,NFELT) , LINNUM(ngrdim,NONODE)
      INTEGER*4  NOLINK(NONODE) ,
C$1-1line changed
     *           NPOINT(3,NFELT) , IPNT  ( NONODE  )
      CHARACTER*40 name
	integer nout
C
      IDUMMY = 0
C
C     Save content of NOLINK and LINNUM to a temporary file
C
      nout = ifileo(0,'dhven991.dat','unformatted','new')
      WRITE  ( nout ) NOLINK
      WRITE  ( nout ) LINNUM
      CLOSE  ( nout )
C
C     makes the pointer and leaves the DELWAQ FROM-TO pointer in
C                                     the LINNUM - NOLINK combination
C
C$1-2lines changed
      CALL DELKNT ( NODNUM , NFELT  , LINNUM , NOLINK , NGRDIM  ,
     *                       IPNT   , NONODE , NPOINT , NOFLOW  )
C
C     Write the so-called FROM-TO file for DELWAQ
C
      call getfnm('ftodwq',name)
      nout = ifileo(0,name,'binary','new')
C$1-5lines deleted
C$1-9lines added
      DO 10 L = 1 , NONODE
         IF ( L .NE. NONODE ) THEN
            NL = NOLINK(L+1)-NOLINK(L)
         ELSE
            NL = NOFLOW-NOLINK(L)
         ENDIF
         IF ( NL .GT. 0 ) WRITE  ( nout )
     *           ( IPNT(L) , LINNUM(K,L) , IDUMMY , IDUMMY , K=1,NL )
   10 CONTINUE
      CLOSE  ( nout )
C
C     retrieve content of NOLINK and LINNUM from the temporary file
C
      nout = ifileo(0,'dhven991.dat','unformatted','old')
      READ   ( nout ) NOLINK
      READ   ( nout ) LINNUM
      CLOSE  ( nout )
C
C     successful return
C
      RETURN
      END
c
c*****************************************************************
c
C$1-2lines changed
      SUBROUTINE DELKNT ( NODNUM , NFELT  , LINNUM , NOLINK , NGRDIM ,
     *                             IPNT   , NONODE , NPOINT , NOFLOW )
C
C     Function: Derives from an array with 3 node numbers per finite
C               element and an array with links per node, a pointer
C               table to the DELWAQ "area" and "flow" file structures.
C
C               The pointer table consists of 3 integers per finite
C               element. The first integer is the array element to save
C               the flow from node number 1 to node number 2 in, the
C               second for the flow from 2 to 3 and the third from 3 to
C
C               The flows and areas have to be summed in an initially
C               clean array. If the pointer value is negative, the flow
C               has to be summed with negative sign on the indicated
C               location.
C$1-2lines added
C               Double pointering has been added, using IPNT from the
C               SHYFEM nodes to the DELWAQ computational elements.
C
C     Creation date    :  2-04-1992 by L. Postma ( Delft Hydraulics )
C     Modification     : 28-04-1992 by L. Postma ( Delft Hydraulics )
C                                   subject: Double pointering active
C                                            DELWAQ elts and links
C                         5-06-1992 by L. Postma ( Delft Hydraulics )
C                                   subject: Write progression to scree
C                                            Correct an agregation erro
C
C     FUNCTIONS CALLED : None
C
C     ROUTINES  CALLED : DELGET - retrieves an exchange number from
C                                 the LINNUM - NOLINK - IPNT structure
C
C     FILES USED       : None
C
C
C     PARAMS  KIND      DIMENSION FUNCTION DESCRIPTION
C     ------  --------- --------- -------- ----------------------------
C     NODNUM  INTEGER*4  3,NFELT  INPUT    Node nrs for each finite elt
C     NFELT   INTEGER*4     1     INPUT    Number of finite elements
C     LINNUM  INTEGER*4 NGRDIM,*  INPUT    Numbers of the links per nod
C     NOLINK  INTEGER*4   NONODE  INPUT    Number of links per node
c     NGRDIM  INTEGER*4     1     INPUT    First dimension of LINNUM
C     IPNT    INTEGER*4   NONODE  INPUT    Pointer to DELWAQ elements
C     NONODE  INTEGER*4     1     INPUT    Number of nodes
C     NPOINT  INTEGER*4  3,NFELT  OUTPUT   Pointers for each finite elt
C     NOFLOW  INTEGER*4     1     OUTPUT   Number of flows identified
C     IERROR  INTEGER*4     1     OUTPUT   Zero = successfull completio
C
C
      INTEGER*4  NODNUM( 3,NFELT ) , LINNUM( NGRDIM,NONODE ) ,
     *           NOLINK(  NONODE ) , NPOINT( 3,NFELT ) , IPNT( NONODE )
C
C     eliminate double pointers from the LINNUM array
C     accumulate number of previous links in NOLINK
C
      NOFLOW = 0
      DO 20 INODE = 1 , NONODE
         IF ( MOD(INODE,500) .EQ. 0 ) WRITE ( * , * )
     *          ' DELWAQ interface is identifying links. Node:' , INODE
         INCR = 0
         IFR  = IPNT(INODE)
         DO 10 ILINK = 1 , NOLINK(INODE)
            ITO = IPNT(LINNUM(ILINK,INODE))
            IF ( ITO .GT. IFR ) THEN
               LINNUM(ILINK,INODE) = 0
               CALL DELGET ( LINNUM , NOLINK , NGRDIM , IPNT  , INODE ,
     *                                         IFR    , ITO   , IQ    )
               IF ( IQ .EQ. 0 ) THEN
                  INCR = INCR + 1
                  LINNUM(INCR,INODE) = ITO
               ENDIF
            ENDIF
            IF ( INCR .LT. ILINK ) LINNUM(ILINK,INODE) = 0
   10    CONTINUE
         NOLINK(INODE) = NOFLOW
         NOFLOW = NOFLOW + INCR
   20 CONTINUE
C
C     set the wanted pointers per finite element
C
      DO 60 IELT = 1 , NFELT
         IF ( MOD(IELT,100) .EQ. 0 ) WRITE ( * , * )
     *          ' DELWAQ interface is numbering links. Element:' , IELT
         DO 50 I0 = 1 , 3
            I1 = IPNT(NODNUM(    I0     ,IELT))
            I2 = IPNT(NODNUM(MOD(I0,3)+1,IELT))
            IQ = 0
            IF ( I1 .NE. I2 )
     *      CALL DELGET ( LINNUM , NOLINK , NGRDIM , IPNT   , NONODE ,
     *                                      I1     , I2     , IQ     )
            NPOINT(I0,IELT) = IQ
   50    CONTINUE
   60 CONTINUE
C
C     successful return
C
      RETURN
      END
c
c*****************************************************************
c
      SUBROUTINE DELGET ( LINNUM , NOLINK , NGRDIM , IPNT   , NONODE ,
     *                                      J1     , J2     , IQ     )
C
C     Function: Derives the exchange number of the link J1 to J2
C               from the array with links in the LINNUM NOLINK
C               combination. Routine returns negative IQ for flow in
C               opposite direction and zero IQ if no matching flow is
C               found.
C
C
C     Creation date    : 28-04-1992 by L. Postma ( Delft Hydraulics )
C     Modification     :  5-06-1992 by L. Postma ( Delft Hydraulics )
C                                   subject: speedup of routine
C                          -  -     by
C                                   subject:
C
C     FUNCTIONS CALLED : None
C
C     ROUTINES  CALLED : None
C
C     FILES USED       : None
C
C
C     PARAMS  KIND      DIMENSION FUNCTION DESCRIPTION
C     ------  --------- --------- -------- ----------------------------
C     LINNUM  INTEGER*4 NGRDIM,*  INPUT    Numbers of the links per nod
C     NOLINK  INTEGER*4   NONODE  INPUT    Number of links per node
C     NGRDIM  INTEGER*4     1     INPUT    First dimension of LINNUM
C     IPNT    INTEGER*4   NONODE  INPUT    Pointer to DELWAQ elements
C     NONODE  INTEGER*4     1     INPUT    Number of nodes
C     J1      INTEGER*4     1     INPUT    From-DELWAQ element nr
C     J2      INTEGER*4     1     INPUT    To  -DELWAQ element nr
C     IQ      INTEGER*4     1     OUTPUT   exchange number
C
C
      INTEGER*4  LINNUM(NGRDIM,NONODE) , NOLINK(NONODE) , IPNT(NONODE)
C
C     Loop through the LINNUM array to find a match
C
      IQ = 0
      JA = J1
      JB = J2
      IF ( J2 .LT. J1 ) THEN
         JA = J2
         JB = J1
      ENDIF
      DO 20 IN = 1 , NONODE
         IJ1 = IPNT(IN)
         IF ( IJ1 .NE. JA ) GOTO 20
         DO 10 IL = 1 , NGRDIM
            IJ2 = LINNUM(IL,IN)
            IF ( IJ2 .EQ.  0 ) GOTO 20
            IF ( IJ2 .EQ. JB ) THEN
               IF ( J1 .EQ. IJ1 ) THEN
                  IQ =    NOLINK(IN) + IL
               ELSE
                  IQ = - (NOLINK(IN) + IL)
               ENDIF
               RETURN
            ENDIF
   10    CONTINUE
   20 CONTINUE
C
C     end of subroutine
C
      RETURN
      END
