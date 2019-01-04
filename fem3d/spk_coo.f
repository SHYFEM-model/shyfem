
!--------------------------------------------------------------------------
!
!                   S P A R S K I T   V E R S I O N  2.
!
! Welcome  to SPARSKIT  VERSION  2.  SPARSKIT is  a  package of  FORTRAN
! subroutines  for working  with  sparse matrices.  It includes  general
! sparse  matrix  manipulation  routines  as  well as  a  few  iterative
! solvers, see detailed description of contents below.
!
!    Copyright (C) 2005  the Regents of the University of Minnesota
!
! SPARSKIT is  free software; you  can redistribute it and/or  modify it
! under the terms of the  GNU Lesser General Public License as published
! by the  Free Software Foundation [version  2.1 of the  License, or any
! later version.]
!
! A copy of  the licencing agreement is attached in  the file LGPL.  For
! additional information  contact the Free Software  Foundation Inc., 59
! Temple Place - Suite 330, Boston, MA 02111, USA or visit the web-site
!
!  http://www.gnu.org/copyleft/lesser.html
!
! DISCLAIMER
! ----------
!
! SPARSKIT  is distributed  in  the hope  that  it will  be useful,  but
! WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
! MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
! Lesser General Public License for more details.
!
! For more information contact saad@cs.umn.edu
! or see https://www-users.cs.umn.edu/~saad/software/SPARSKIT/
!
!    This file is part of SHYFEM.
!
!    This file has been compiled from more than one source files.
!    The original files are in FORMATS and are called formats.f and unary.f
!
!--------------------------------------------------------------------------

c----------------------------------------------------------------------- 
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
c----------------------------------------------------------------------- 
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
c-----------------------------------------------------------------------
c  Coordinate     to   Compressed Sparse Row 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c--------- 
c nrow	= dimension of the matrix 
c nnz	= number of nonzero elements in matrix
c a,
c ir, 
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c 	  the elements, ir(k) = its row number and jc(k) = its column 
c	  number. The order of the elements is arbitrary. 
c
c on return:
c----------- 
c ir 	is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao 
c 	continung the real values, jao containing the column indices, 
c	and iao being the pointer to the beginning of the row, 
c	in arrays ao, jao.
c
c Notes:
c------ This routine is NOT in place.  See coicsr
c       On return the entries  of each row are NOT sorted by increasing 
c       column number
c
c------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c------------- end of coocsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
c-----------------------------------------------------------------------
      subroutine clncsr(job,value2,nrow,a,ja,ia,indu,iwk)
c     .. Scalar Arguments ..
      integer job, nrow, value2
c     ..
c     .. Array Arguments ..
      integer ia(nrow+1),indu(nrow),iwk(nrow+1),ja(*)
      real*8  a(*)
c     ..
c
c     This routine performs two tasks to clean up a CSR matrix
c     -- remove duplicate/zero entries,
c     -- perform a partial ordering, new order lower triangular part,
c        main diagonal, upper triangular part.
c
c     On entry:
c
c     job   = options
c         0 -- nothing is done
c         1 -- eliminate duplicate entries, zero entries.
c         2 -- eliminate duplicate entries and perform partial ordering.
c         3 -- eliminate duplicate entries, sort the entries in the
c              increasing order of clumn indices.
c
c     value2  -- 0 the matrix is pattern only (a is not touched)
c                1 matrix has values too.
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in CSR format
c
c     On return:
c     a,ja,ia -- cleaned matrix
c     indu    -- pointers to the beginning of the upper triangular
c                portion if job > 1
c
c     Work space:
c     iwk     -- integer work space of size nrow+1
c
c     .. Local Scalars ..
      integer i,j,k,ko,ipos,kfirst,klast
      real*8  tmp
c     ..
c
      if (job.le.0) return
c
c     .. eliminate duplicate entries --
c     array INDU is used as marker for existing indices, it is also the
c     location of the entry.
c     IWK is used to stored the old IA array.
c     matrix is copied to squeeze out the space taken by the duplicated
c     entries.
c
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = ia(i)
 90   continue
      iwk(nrow+1) = ia(nrow+1)
      k = 1
      do 120 i = 1, nrow
         ia(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = ja(ipos)
            if (indu(j).eq.0) then
c     .. new entry ..
               if (value2.ne.0) then
                  if (a(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     ja(k) = ja(ipos)
                     a(k) = a(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  ja(k) = ja(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
c     .. duplicate entry ..
               a(indu(j)) = a(indu(j)) + a(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
c     .. remove marks before working on the next row ..
         do 110 ipos = ia(i), k - 1
            indu(ja(ipos)) = 0
 110     continue
 120  continue
      ia(nrow+1) = k
      if (job.le.1) return
c
c     .. partial ordering ..
c     split the matrix into strict upper/lower triangular
c     parts, INDU points to the the beginning of the upper part.
c
      do 140 i = 1, nrow
         klast = ia(i+1) - 1
         kfirst = ia(i)
 130     if (klast.gt.kfirst) then
            if (ja(klast).lt.i .and. ja(kfirst).ge.i) then
c     .. swap klast with kfirst ..
               j = ja(klast)
               ja(klast) = ja(kfirst)
               ja(kfirst) = j
               if (value2.ne.0) then
                  tmp = a(klast)
                  a(klast) = a(kfirst)
                  a(kfirst) = tmp
               endif
            endif
            if (ja(klast).ge.i)
     &         klast = klast - 1
            if (ja(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
c
         if (ja(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
c
c     .. order the entries according to column indices
c     burble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = ia(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), ia(i+1)-1
            do 170 j = ia(i+1)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
      return
c---- end of clncsr ----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine csort (n,a,ja,ia,values) 
      logical values
      integer n, ja(*), ia(n+1)
      real*8 a(*) 
c-----------------------------------------------------------------------
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within 
c each row. It uses insertion sort
c-----------------------------------------------------------------------
c on entry:
c--------- 
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows. 
c values= logical indicating whether or not the real values a(*) must 
c         also be permuted. if (.not. values) then the array a is not
c         touched by csort and can be a dummy array. 
c 
c on return:
c----------
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c----------------------------------------------------------------------- 
c Y. Saad - recoded Dec. 20th -  2017 - 
c-----------------------------------------------------------------------
c local variables
      integer row, i, k, j
      real*8 rj
      rj = 0.0 
c
c for each row
c
      do row=1, n
c
c scan row and do an insertion sort
c 
         do k=ia(row)+1, ia(row+1)-1 
            j = ja(k)
            if (values) rj = a(k)
            i = k-1;
            do while ((i .ge. ia(row)) .and. (j < ja(i)) ) 
               ja(i+1) = ja(i)
               if (values) a(i+1) = a(i)
               i = i-1
               if (i .eq. 0) exit
            enddo
            ja(i+1) = j
            if (values) a(i+1) = rj
         enddo
      enddo
      return 
c---------------end-of-csort-------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 
      subroutine dvperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of a real vector x 
c according to the permutation array perm(*), i.e., on return, 
c the vector x satisfies,
c
c	x(perm(j)) :== x(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c x	= input vector
c
c on return:
c---------- 
c x	= vector x permuted according to x(perm(*)) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables 
      real*8 tmp, tmp1
c
      init      = 1
      tmp       = x(init)        
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c     
c loop
c 
 6    k = k+1
c
c save the chased element --
c 
      tmp1      = x(ii) 
      x(ii)     = tmp
      next          = perm(ii) 
      if (next .lt. 0 ) goto 65
c     
c test for end 
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
c
c end loop 
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp        = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
c     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
c     
      return
c-------------------end-of-dvperm--------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine ivperm (n, ix, perm) 
      integer n, perm(n), ix(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of an integer vector 
c ix according to the permutation array perm(*), i.e., on return, 
c the vector x satisfies,
c
c	ix(perm(j)) :== ix(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c ix	= input vector
c
c on return:
c---------- 
c ix	= vector x permuted according to ix(perm(*)) :=  ix(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      integer tmp, tmp1
c
      init      = 1
      tmp        = ix(init)        
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c     
c loop
c 
 6    k = k+1
c
c save the chased element --
c 
      tmp1          = ix(ii) 
      ix(ii)     = tmp
      next          = perm(ii) 
      if (next .lt. 0 ) goto 65
c     
c test for end 
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
c
c end loop 
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp        = ix(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
c     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
c     
      return
c-------------------end-of-ivperm--------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      real*8 a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
c----------------------------------------------------------------------- 
c Banded (Linpack ) format   to    Compressed Sparse Row  format.
c----------------------------------------------------------------------- 
c on entry:
c----------
c n	= integer,the actual row dimension of the matrix.
c
c nabd  = first dimension of array abd.
c
c abd   = real array containing the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix,comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
c         in row lowd (see below). 
c    
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located. 
c         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
c         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
c
c ml	= integer. equal to the bandwidth of the strict lower part of A
c mu	= integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than nabd then an error 
c         message is set. see ierr.
c
c len   = integer. length of arrays a and ja. bndcsr will stop if the
c         length of the arrays a and ja is insufficient to store the 
c         matrix. see ierr.
c
c on return:
c-----------
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c lowd  = if on entry lowd was zero then lowd is reset to the default
c         value ml+mu+l. 
c
c ierr  = integer. used for error message output. 
c         ierr .eq. 0 :means normal return
c         ierr .eq. -1 : means invalid value for lowd. 
c	  ierr .gt. 0 : means that there was not enough storage in a and ja
c         for storing the ourput matrix. The process ran out of space 
c         (as indicated by len) while trying to fill row number ierr. 
c         This should give an idea of much more storage might be required. 
c         Moreover, the first irow-1 rows are correctly filled. 
c
c notes:  the values in abd found to be equal to zero
c -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
c         The resulting may not be identical to a csr matrix
c         originally transformed to a bnd format.
c          
c----------------------------------------------------------------------- 
      ierr = 0
c-----------
      if (lowd .gt. nabd .or. lowd .le. 0) then 
         ierr = -1
         return
      endif
c-----------
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
c-----------------------------------------------------------------------
         i = lowd 
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j) 
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then 
               ierr = irow 
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
c     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
c------------- end of bndcsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
