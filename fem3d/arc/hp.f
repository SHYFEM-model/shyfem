c
c $Id: hp.f,v 1.23 2003/07/31 11:14:10 georg Exp $
c
c finite element model hp (version 2D)
c
c revision log :
c
c version from may'88
c revised on 27.07.88 by ggu  (sp136(ic),initial cond. to file 7)
c revised on 17.08.88 by ggu  (output of first record, getpar, descrp)
c revised on 18.08.88 by ggu  (idtfl,nfloat...)
c revised on 30.08.88 by ggu  (itflst,write to fl file only on idtfl != 0)
c revised on 31.08.88 by ggu  (sp110 moved, even if not necessary because
c                               ...of changes in sp110,sp156f)
c revised on 11.11.88 by ggu  (namelist,/time/)
c revised on 17.11.88 by ggu  (v6v,dimension rzv, nrbdim, boundn)
c revised on 30.11.88 by ggu  (nknaus to nexdim)
c revised on 05.12.88 by ggu  (rrv)
c revised on 12.12.88 by ggu  (idtfl in dlfl...)
c revised on 17.02.89 by ggu  (priout, itflst, itfdst in lgr routines)
c revised on 18.02.89 by ggu  (iclose, iczv in subroutines)
c revised on 10.03.89 by ggu  (fd completed)
c revised on 23.05.89 by ggu  (pressure terms)
c revised on 12.04.90 by ggu  (cstini,cstset,fcor & href in cstset)
c revised on 28.07.90 by ggu  (cstfil)
c revised on 24.08.90 by ggu  (/dltimi/)
c revised on 05.02.91 by ggu  ( bnd(10*...) --> bnd(ibndim,..) ; ibndim )
c revised on 21.12.93 by ggu  $$RESID -> residual currents
c revised on 17.01.94 by ggu  $$conz introduction of concentration
c revised on 31.01.94 by ggu  $$conzbc introduction of bc for concentration
c revised on 19.05.97 by ggu  $$VERNEW no call to vernew anymore
c revised on 19.05.97 by ggu  $$SP136 call to sp136 commented
c revised on 28.05.97 by ggu  $$IPCCV close data items commented or deleted
c revised on 16.06.97 by ggu  $$RESCU -> data structures for rescu() deleted
c				(urrv,vrrv,xrv,urv,vrv,zrv)
c revised on 23.09.97 by ggu  introduced conzn -> name of file for conz values
c revised on 03.12.97 by ggu  ibndim = 100 (old 50)
c revised on 03.12.97 by ggu  $$ST new tempn,saltn (see conzn)
c revised on 04.12.97 by ggu  $$ST new salt/temp: soev,toev,snv,tnv,rsv,rtv 
c 22.01.1998	ggu	custom routine called at end of time step
c 20.03.1998    ggu     custom routine in own file to avoid compiler warnings
c 27.04.1998    ggu     sp110a and tstlnk befor calling cstset -> cancelled
c 28.04.1998    ggu     new /kfluxc/, /iflux/, call to flxini
c 29.04.1998    ggu     uses module for semi-implicit time-step
c 30.04.1998    ggu     finally eliminated /semimp/, /trock/, /ffloat/
c 28.05.1998    ggu     no iamat, inddim,..., new call to sp131g()
c 28.05.1998    ggu     new ruv, rvv (momentum input)
c 17.06.1998	ggu	restructured for documentation
c 20.06.1998	ggu	nvers eliminated definitly
c 13.07.1998	ggu	ruv, rvv treated as node arrays (nkndim)
c 12.08.1998	ggu	uedif, vedif for horizontal diffusion (no dudxv...)
c 20.08.1998    ggu     iextpo finally eliminated
c 21.08.1998    ggu     routine setczg removed
c 21.08.1998    ggu     new arrays zov, znv, up0v, vp0v to substitute xv
c 21.08.1998    ggu     xv eliminated
c 06.11.1998    ggu     new arrays hkv and hev
c 22.01.1999    ggu     new oxygen module
c 26.01.1999    ggu     new 3D compatibility arrays
c 20.04.1999    ggu     converted to stress instead of wind (tauxnv...)
c 14.09.1999    ggu     call to setflux introduced, new hlhv
c 19.11.1999    ggu     new routines for section vol
c 20.01.2000    ggu     common block /dimdim/ eliminated
c
c******************************************************************

c----------------------------------------------------------------

	program hp

	include 'param.h'

c----------------------------------------------------------------

	parameter (ibndim=100)
	parameter (isym=1)	!change also in sub555

	parameter (matdim=nkndim*(1+(2-isym)*mbwdim))

c description of simulation and basin

	character*80 descrp
	character*80 descrr
	common /descrp/ descrp
	common /descrr/ descrr

c parameter values

	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /mkonst/ eps1,eps2,pi,flag,high,higi
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it

c 3D compatibility arrays

	integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(nkndim)
        common /ilhkv/ilhkv
        integer ilhv(neldim)
        common /ilhv/ilhv
        real hlv(nlvdim), hldv(nlvdim)
        common /hlv/hlv, /hldv/hldv
        real hlhv(neldim)
	common /hlhv/hlhv

c special arrays

	common /knausc/ knausm,knaus(nexdim)

	common /kfluxc/ nsect,kfluxm,kflux(nfxdim)
	common /iflux/ iflux(3,nfxdim)

        integer nvols,kvold,kvolm,kvol(nfxdim)
        common /kvolc/ nvols,kvold,kvolm,kvol
        integer ivolm,ivol(nfxdim)
        common /ivol/ivolm,ivol

	common /chezy/ nczdum,czdum(6,nardim+1)

c arrays from basin

	common /nen3v/nen3v(3,neldim)
	common /iarv/iarv(neldim)
	common /ipv/ipv(nkndim), /ipev/ipev(neldim)
	common /hm3v/hm3v(3,neldim)
	common /xgv/xgv(nkndim), /ygv/ygv(nkndim)

c frcition parameters

	common /czv/czv(neldim)
	common /ausv/ausv(neldim)

c hydrodynamical values (transports and levels)

	real uov(neldim),vov(neldim),unv(neldim),vnv(neldim)
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
	real zenv(3,neldim), zeov(3,neldim)	!$$conz
	common /zenv/zenv, /zeov/zeov
	real zov(nkndim), znv(nkndim)
	common /zov/zov, /znv/znv
	real up0v(nkndim), vp0v(nkndim)
	common /up0v/up0v, /vp0v/vp0v

	common /xv/xv(3,nkndim)

c depth arrays

        real hkv(nkndim), hev(neldim)
        common /hkv/hkv, /hev/hev

c T/D values (concentration, salinity, temperature)

	common /coev/coev(3,neldim), /cnv/cnv(nkndim)	!$$conz
	common /soev/soev(3,neldim), /snv/snv(nkndim)	!$$ST
	common /toev/toev(3,neldim), /tnv/tnv(nkndim)	!$$ST

c local boundary arrays

	common /bnd/bnd(ibndim,nbcdim)
	common /ierv/ierv(2,nrbdim)
	common /rhv/rhv(nrbdim), /rlv/rlv(nrbdim)
	common /rrv/rrv(nrbdim), /irv/irv(nrbdim)

c global boundary conditions

c	file names

	character*80 boundn(nbcdim)
	character*80 conzn(nbcdim)
	character*80 saltn(nbcdim)
	character*80 tempn(nbcdim)
	common /boundn/boundn
	common /conzn/conzn
	common /saltn/saltn
	common /tempn/tempn

	character*80 bio2dn(nbcdim)
	common /bio2dn/bio2dn

c	level and flux conditions

	real rzv(nkndim)
	real rqv(nkndim)
	common /rzv/rzv
	common /rqv/rqv

c	conz, salt, temp conditions

	real rcv(nkndim)
	real rsv(nkndim)
	real rtv(nkndim)
	common /rcv/rcv
	common /rsv/rsv
	common /rtv/rtv

c	momentum input

	real ruv(nkndim)
	real rvv(nkndim)
	common /ruv/ruv
	common /rvv/rvv

c geometric and topographic arrays

	common /ev/ev(13,neldim)
	common /iwegv/iwegv(neldim)
	common /ieltv/ieltv(3,neldim)
	common /kantv/kantv(2,nkndim)
	common /dxv/dxv(nkndim), /dyv/dyv(nkndim)
	common /inodv/inodv(nkndim)

c link index

	common /ilinkv/ilinkv(nkndim+1)
	common /lenkv/lenkv(2,nlidim)

c wind and pressure arrays

        common /tauxnv/tauxnv(nkndim), /tauynv/tauynv(nkndim)
	common /wxov/wxov(nkndim), /wyov/wyov(nkndim)
	common /wxnv/wxnv(nkndim), /wynv/wynv(nkndim)
	common /ppv/ppv(nkndim)
	common /pov/pov(nkndim), /pnv/pnv(nkndim)

c austausch arrays

	common /uedif/uedif(neldim), /vedif/vedif(neldim)

c new depth and area arrays

        real hdknv(nkndim)
        common /hdknv/hdknv
        real hdkov(nkndim)
        common /hdkov/hdkov

        real hdenv(neldim)
        common /hdenv/hdenv
        real hdeov(neldim)
        common /hdeov/hdeov

        real areakv(nkndim)
        common /areakv/areakv
  
c auxiliary values

	common /v1v/v1v(nkndim), /v2v/v2v(nkndim)
	common /v3v/v3v(nkndim), /v4v/v4v(nkndim)
	common /v5v/v5v(nkndim), /v6v/v6v(nkndim)
	common /v7v/v7v(nkndim)

	common /ve1v/ve1v(neldim)

c distance record

        real rdistv(nkndim)
	common /rdistv/rdistv

c global system matrix

	common /amat/amat(matdim)

c floating module
c
c	common /ieflv/ieflv(nfldim)
c	common /xpflv/xpflv(nfldim), /ypflv/ypflv(nfldim)
c
c	common /iefdv/iefdv(nfddim)
c	common /xpfdv/xpfdv(nfddim), /ypfdv/ypfdv(nfddim)

c closing module
c
c	common /ipccv/ipccv(ipcdim), /ivccv/ivccv(ivcdim)	!$$IPCCV

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call vers2d	!set version of model

c read input files, control input values %%%%%%%%%%%%%%%%%%%%%%%%
c
	call cstdim(nkndim,neldim,nrbdim,nbcdim
     +                  ,mbwdim,ngrdim,nardim,nexdim
     +                  ,nfldim,nfddim,ipcdim,ivcdim
     +                  ,ndldim,nfxdim)
c
	call cstinit
	call cstfile(nkndim,neldim)
c
	call sp131a(nkndim,neldim)
	call sp131b(nrbdim,nbcdim)
	call sp131k(matdim)
	call sp131g(matdim,mbwdim,nkn,mbw,isym)

	call cstcheck
c
c initialize 3D compatibility arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call comp3d
c
c initialize boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call sp111(1)
	call setxv
c
c initialize triangles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call sp110a
c
c set up pointers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call tstlnk
c
c initialize chezy values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call sp135r
c
c initialize lagrangean float tracking %%%%%%%%%%%%%%%%%%%%%%%%%%
c
c	call lgrini(nfldim,nfddim)
c
c initialize modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call cstsetup
c
c write input values to log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call prilog
c
c dry nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call setweg(0,n)
c
c set up depths and areas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
        call setarea(1,areakv)
        call setdepth(1,hdknv,hdenv,zenv,areakv)                            
c
c inlets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call sp136(ic)		!$$SP136
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call shdist(rdistv)

	do while(it+idt.le.itend)

	   niter=niter+1
	   it=it+idt

	   call dobefor

	   call sp111(2)            !boundary conditions

	   call sp159b              !hydro

           call setdepth(1,hdknv,hdenv,zenv,areakv)     !maybe into sp159b

	   call con2sh		    !transport/diffusion !$$conz
	   call sed2sh		    !sediments

	   call setflux		    !fluxes for fects

	   call bio2d(it,idt)

	   call doafter

	   call pritime             !output to terminal

	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call endtime

	end

c*****************************************************************

