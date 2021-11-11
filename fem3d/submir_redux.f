
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997,2010,2012,2015,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c linear equation solvers
c
c contents :
c
c      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
c      SUBROUTINE LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
c      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
c      SUBROUTINE LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
c      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
c      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)
c      SUBROUTINE LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
c      SUBROUTINE DEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
c      SUBROUTINE DEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
c      SUBROUTINE DEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
c      SUBROUTINE DEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
c      SUBROUTINE DUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
c      SUBROUTINE DUELPB (UL,B,N,NC,IA,X)
c      SUBROUTINE DUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
c      SUBROUTINE UERTST (IER,NAME)
c      SUBROUTINE UGETIO (IOPT,NIN,NOUT)
c      SUBROUTINE GELB(R,A,M,N,MUD,MLD,EPS,IER)
c      SUBROUTINE DGELB(R,A,M,N,MUD,MLD,EPS,IER)
c      SUBROUTINE MCHB(R,A,M,N,MUD,IOP,EPS,IER)
c      SUBROUTINE DMCHB(R,A,M,N,MUD,IOP,EPS,IER)
c      subroutine loctst(i,j,n,m)
c      function locmy(i,j,n,m)
c      function locimm(i,j,n,m)
c      function loccer(i,j,n,m)
c      function locssp(i,j,n,m)
c      function locsps(i,j,n,m)
c
c revision log :
c
c 03.04.1997	ggu	general - compiler warnings for gfortran
c 24.05.1997	ggu	general - compiler warnings -> call to dp routines
c 23.03.2010	ggu	changed v6.1.1
c 30.03.2012	ggu	changed VERS_6_1_51
c 01.06.2012	ggu	changed VERS_6_1_53
c 18.12.2015	ggu	changed VERS_7_3_17
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 09.11.2021	ggu	reduced to esential routines
c
c*************************************************************************
c
        function locssp(i,j,n,m)
c
c access ssp routines
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locssp  position of element in band matrix
c
        implicit none
c
	integer locssp
        integer i,j,n,m
c
	if(i-j.gt.m.or.j-i.gt.m) then
	  locssp = 0
	else if(n.le.m) then	!this is for a full matrix
	  locssp = n*(i-1) + j
        else if(i.lt.m) then
          locssp = 2*m*i - m + j - m*(m+1)/2 + (m-i)*(m-i+1)/2
        else if(i.gt.n-m+1) then
          locssp = 2*m*i - m + j - m*(m+1)/2
     +                - (i-(n-m+1))*(i-(n-m))/2
        else
          locssp = 2*m*i - m + j - m*(m+1)/2
        end if
c
        return
        end
c
c*************************************************************************
c
        function locsps(i,j,n,m)
c
c access ssp routines (symmetric compressed storage mode)
c
c only main and upper diagonals - if an element in the lower
c ...diagonals is referenced, 0 is returned
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locsps  position of element in band matrix
c
c original formula : locsps = (1+m)*(i-1)+abs(j-i)+1
c
        implicit none
c
	integer locsps
        integer i,j,n,m
c
	if(i.gt.j.or.j-i.gt.m) then
	  locsps = 0
        else if(i.gt.n-m+1) then
          locsps = m*(i-1)+j
     +                - (i-(n-m+1))*(i-(n-m))/2
        else
          locsps = m*(i-1)+j
        end if
c
        return
        end

c*************************************************************************

