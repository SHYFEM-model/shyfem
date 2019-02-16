
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c-------------------------------------------------------------
c header file for ous file format
c-------------------------------------------------------------
c
c ftype		file id
c maxvers	newest version of file
c maxcomp	compatible function calls down to here
c
c ndim		number of possible entries (open files)
c nitdim	number of integer values to be stored
c nchdim	number of string values to be stored
c
c ousitem	number of maximum entries in table
c
c ousvar	integer parameters of open files
c ouschar	string parameters of open files
c
c ousvar(0,n)	iunit
c ousvar(1,n)	nvers
c ousvar(2,n)	nkn
c ousvar(3,n)	nel
c ousvar(4,n)	nlv
c ousvar(5,n)	not used
c ousvar(6,n)	date
c ousvar(7,n)	time
c
c ousreal(1,n)	href
c ousreal(2,n)	hzmin
c
c ouschar(1,n)	title
c ouschar(2,n)	femver
c
c-------------------------------------------------------------
c parameters
c-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=27,maxvers=2,maxcomp=1)

        integer ndim,nitdim,nrldim,nchdim
        parameter(ndim=10,nitdim=7,nrldim=2,nchdim=2)

c-------------------------------------------------------------
c common
c-------------------------------------------------------------

        integer ousitem
        common /ousitm/ousitem

        integer ousvar(0:nitdim,ndim)
        common /ousvar/ousvar

        real ousreal(nrldim,ndim)
        common /ousreal/ousreal

        character*80 ouschar(nchdim,ndim)
        common /ouschar/ouschar

c-------------------------------------------------------------
c save
c-------------------------------------------------------------

        save /ousitm/
        save /ousvar/
        save /ousreal/
        save /ouschar/

c-------------------------------------------------------------
c end of header
c-------------------------------------------------------------

