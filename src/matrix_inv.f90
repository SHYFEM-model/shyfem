!
! $Id: subinv.f,v 1.4 2008-04-24 09:16:52 georg Exp $
!
! matrix inversion routines (non symmetric band matrix)
!
! contents :
!
! sgbco			factors and solves
! sgbfa			factors and solves (without condition number)
! sgbsl			solves the double precision band system using factorization
!
! dgbco			factors and solves
! dgbfa			factors and solves (without condition number)
! dgbsl			solves the double precision band system using factorization
!
! lp_init_system	initializes matrix and vector
! lp_solve_system	factors and solves for solution
! lp_subst_system	solves for solution
! lp_mult_band		multiplies matrix with vector
!
! dlp_init_system	initializes matrix and vector
! dlp_solve_system	factors and solves for solution
! dlp_subst_system	solves for solution
! dlp_mult_band		multiplies matrix with vector
!
! loclp			finds position in matrix
!
! revision log :
!
! 02.04.2007    ggu     assembled from lapack
! 06.06.2007    ggu     new routines for back substitution and initialization
! 22.04.2008    ggu     new SAXPY_NEW for parallelization trial
!
!*************************************************************************
!---------------------------------------------------------------------------------
      module matrix_inv
!---------------------------------------------------------------------------------
      contains
!---------------------------------------------------------------------------------

      subroutine sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)

      integer lda,n,ml,mu,ipvt(1)
      double precision abd(lda,1),z(1)
      double precision rcond

!     sgbco factors a double precision band matrix by gaussian
!     elimination and estimates the condition of the matrix.
!
!     if  rcond  is not needed, sgbfa is slightly faster.
!     to solve  a*x = b , follow sgbco by sgbsl.
!     to compute  inverse(a)*c , follow sgbco by sgbsl.
!     to compute  determinant(a) , follow sgbco by sgbdi.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   double precision
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       double precision(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     example..  if the original matrix is
!
!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66
!
!      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
!
!            *  *  *  +  +  +  , * = not used
!            *  * 13 24 35 46  , + = used for pivoting
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack sgbfa
!     blas saxpy,sdot,sscal,sasum
!     fortran abs,amax1,max0,min0,sign
!
!     internal variables
!
      double precision ek,t,wk,wkm
      double precision s,sm,ynorm
      real anorm
      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
!
!
!     compute 1-norm of a
!
      anorm = 0.0e0
      l = ml + 1
      is = l + mu
      do 10 j = 1, n
         anorm = amax1(anorm,real(sasum(l,abd(is,j),1)))
         if (is .gt. ml + 1) is = is - 1
         if (j .le. mu) l = l + 1
         if (j .ge. n - ml) l = l - 1
   10 continue
!
!     factor
!
      call sgbfa(abd,lda,n,ml,mu,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!     solve trans(u)*w = e
!
      ek = 1.0e0
      do 20 j = 1, n
         z(j) = 0.0e0
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .le. abs(abd(m,k))) go to 30
            s = abs(abd(m,k))/abs(ek-z(k))
            call sscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (abd(m,k) .eq. 0.0e0) go to 40
            wk = wk/abd(m,k)
            wkm = wkm/abd(m,k)
         go to 50
   40    continue
            wk = 1.0e0
            wkm = 1.0e0
   50    continue
         kp1 = k + 1
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .gt. ju) go to 90
            do 60 j = kp1, ju
               mm = mm - 1
               sm = sm + abs(z(j)+wkm*abd(mm,j))
               z(j) = z(j) + wk*abd(mm,j)
               s = s + abs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               mm = m
               do 70 j = kp1, ju
                  mm = mm - 1
                  z(j) = z(j) + t*abd(mm,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
!
!     solve trans(l)*y = w
!
      do 120 kb = 1, n
         k = n + 1 - kb
         lm = min0(ml,n-k)
         if (k .lt. n) z(k) = z(k) + sdot(lm,abd(m+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
!
      ynorm = 1.0e0
!
!     solve l*v = y
!
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min0(ml,n-k)
         if (k .lt. n) call saxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
!
!     solve  u*z = w
!
      do 160 kb = 1, n
         k = n + 1 - kb
         if (abs(z(k)) .le. abs(abd(m,k))) go to 150
            s = abs(abd(m,k))/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (abd(m,k) .ne. 0.0e0) z(k) = z(k)/abd(m,k)
         if (abd(m,k) .eq. 0.0e0) z(k) = 1.0e0
         lm = min0(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call saxpy(lm,t,abd(la,k),1,z(lz),1)
  160 continue
!     make znorm = 1.0
      s = 1.0e0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
!
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end

!*************************************************************************

      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)

      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)

!     sgbfa factors a double precision band matrix by elimination.
!
!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sscal,isamax
!     fortran max0,min0
!
!     internal variables
!
      double precision t
      integer i,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do jz = j0, j1
         i0 = m + 1 - jz
         do i = i0, ml
            abd(i,jz) = 0.0e0
	 end do
      end do
   30 continue
      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130

!---------------------------------------------------------
! start of main loop
!---------------------------------------------------------

      do k = 1, nm1
         kp1 = k + 1
!
!        zero next fill-in column
!
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do i = 1, ml
               abd(i,jz) = 0.0e0
	    end do
   50    continue
!
!        find l = pivot index
!
         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
         if (abd(l,k) .ne. 0.0e0) then
!
!           interchange if necessary
!
            if (l .ne. m) then
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
	    end if
!
!           compute multipliers
!
            t = -1.0e0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)
!
!           row elimination with column indexing
!
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            do j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .ne. mm) then
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
	       end if

               call saxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)

	!       do i=1,lm
	!	 abd(mm+i,j) = abd(mm+i,j) + t * abd(m+i,k)
	!       end do

	    end do
	 else
            info = k
	 end if

       end do

!---------------------------------------------------------
! end of main loop
!---------------------------------------------------------

  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0e0) info = n
      return
      end

!*************************************************************************

      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)

!     sgbsl solves the double precision band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgbco or sgbfa.
!
!     on entry
!
!        abd     double precision(lda, n)
!                the output from sgbco or sgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from sgbco or sgbfa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sgbco has set rcond .gt. 0.0
!        or sgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sdot
!     fortran min0
!
!     internal variables
!
      double precision t
      integer k,kb,l,la,lb,lm,m,nm1
!
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call saxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call saxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end

!*************************************************************************

      subroutine dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)

      integer lda,n,ml,mu,ipvt(1)
      double precision abd(lda,1),z(1)
      double precision rcond

!     dgbco factors a double precision band matrix by gaussian
!     elimination and estimates the condition of the matrix.
!
!     if  rcond  is not needed, dgbfa is slightly faster.
!     to solve  a*x = b , follow dgbco by dgbsl.
!     to compute  inverse(a)*c , follow dgbco by dgbsl.
!     to compute  determinant(a) , follow dgbco by dgbdi.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   double precision
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       double precision(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     example..  if the original matrix is
!
!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66
!
!      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
!
!            *  *  *  +  +  +  , * = not used
!            *  * 13 24 35 46  , + = used for pivoting
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack dgbfa
!     blas daxpy,ddot,dscal,dasum
!     fortran dabs,dmax1,max0,min0,dsign
!
!     internal variables
!
      double precision ek,t,wk,wkm
      double precision anorm,s,sm,ynorm
      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
!
!
!     compute 1-norm of a
!
      anorm = 0.0d0
      l = ml + 1
      is = l + mu
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(l,abd(is,j),1))
         if (is .gt. ml + 1) is = is - 1
         if (j .le. mu) l = l + 1
         if (j .ge. n - ml) l = l - 1
   10 continue
!
!     factor
!
      call dgbfa(abd,lda,n,ml,mu,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!     solve trans(u)*w = e
!
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(abd(m,k))) go to 30
            s = dabs(abd(m,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (abd(m,k) .eq. 0.0d0) go to 40
            wk = wk/abd(m,k)
            wkm = wkm/abd(m,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .gt. ju) go to 90
            do 60 j = kp1, ju
               mm = mm - 1
               sm = sm + dabs(z(j)+wkm*abd(mm,j))
               z(j) = z(j) + wk*abd(mm,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               mm = m
               do 70 j = kp1, ju
                  mm = mm - 1
                  z(j) = z(j) + t*abd(mm,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
!
!     solve trans(l)*y = w
!
      do 120 kb = 1, n
         k = n + 1 - kb
         lm = min0(ml,n-k)
         if (k .lt. n) z(k) = z(k) + ddot(lm,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
!
      ynorm = 1.0d0
!
!     solve l*v = y
!
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min0(ml,n-k)
         if (k .lt. n) call daxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
!
!     solve  u*z = w
!
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(abd(m,k))) go to 150
            s = dabs(abd(m,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (abd(m,k) .ne. 0.0d0) z(k) = z(k)/abd(m,k)
         if (abd(m,k) .eq. 0.0d0) z(k) = 1.0d0
         lm = min0(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call daxpy(lm,t,abd(la,k),1,z(lz),1)
  160 continue
!     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
!
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end

!*************************************************************************

      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)

      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)

!     dgbfa factors a double precision band matrix by elimination.
!
!     dgbfa is usually called by dgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgbsl will divide by zero if
!                     called.  use  rcond  in dgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!     fortran max0,min0
!
!     internal variables
!
      double precision t
      integer i,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
!
!        zero next fill-in column
!
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
!
!        find l = pivot index
!
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
         if (abd(l,k) .eq. 0.0d0) go to 100
!
!           interchange if necessary
!
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
!
!           compute multipliers
!
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
!
!           row elimination with column indexing
!
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end

!*************************************************************************

      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)

!     dgbsl solves the double precision band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgbco or dgbfa.
!
!     on entry
!
!        abd     double precision(lda, n)
!                the output from dgbco or dgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from dgbco or dgbfa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgbco has set rcond .gt. 0.0
!        or dgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran min0
!
!     internal variables
!
      double precision t
      integer k,kb,l,la,lb,lm,m,nm1
!
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end

!*************************************************************************

      double precision FUNCTION SASUM(N,SX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      double precision SX(*)
!     ..
!
!  Purpose
!  =======
!
!     takes the sum of the absolute values.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!

!     .. Local Scalars ..
      double precision STEMP
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
!     ..
      SASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          STEMP = STEMP + ABS(SX(I))
   10 CONTINUE
      SASUM = STEMP
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,6)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
        STEMP = STEMP + ABS(SX(I))
   30 CONTINUE
      IF (N.LT.6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        STEMP = STEMP + ABS(SX(I)) + ABS(SX(I+1)) + ABS(SX(I+2)) + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50 CONTINUE
   60 SASUM = STEMP
      RETURN
      END

!*************************************************************************

      SUBROUTINE SSCAL(N,SA,SX,INCX)
!     .. Scalar Arguments ..
      double precision SA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      double precision SX(*)
!     ..
!
!  Purpose
!  =======
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END

!*************************************************************************

      SUBROUTINE SAXPY_NEW(N,SA,SX,INCX,SY,INCY)

      IMPLICIT NONE

      double precision SA
      INTEGER INCX,INCY,N
      double precision SX(*),SY(*)

      INTEGER I

      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      IF (INCX.NE.1 .OR. INCY.NE.1) GO TO 99

! !$OMP PARALLEL PRIVATE(i)
! !$OMP DO SCHEDULE(STATIC)

      DO I = 1,N
          SY(I) = SY(I) + SA*SX(I)
      END DO

! !$OMP END DO NOWAIT
! !$OMP END PARALLEL

      RETURN
   99 CONTINUE
      write(6,*) INCX,INCY
      stop 'error stop SAXPY: increments'
      END

!*************************************************************************

      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
!     .. Scalar Arguments ..
      double precision SA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      double precision SX(*),SY(*)
!     ..
!
!  Purpose
!  =======
!
!     SAXPY constant times a vector plus a vector.
!     uses unrolled loop for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          SY(IY) = SY(IY) + SA*SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          SY(I) = SY(I) + SA*SX(I)
          SY(I+1) = SY(I+1) + SA*SX(I+1)
          SY(I+2) = SY(I+2) + SA*SX(I+2)
          SY(I+3) = SY(I+3) + SA*SX(I+3)
   50 CONTINUE
      RETURN
      END

!*************************************************************************

      INTEGER FUNCTION ISAMAX(N,SX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      double precision SX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      double precision SMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
      ISAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ISAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (ABS(SX(IX)).LE.SMAX) GO TO 5
          ISAMAX = I
          SMAX = ABS(SX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
          IF (ABS(SX(I)).LE.SMAX) GO TO 30
          ISAMAX = I
          SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END

!*************************************************************************

      double precision FUNCTION SDOT(N,SX,INCX,SY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      double precision SX(*),SY(*)
!     ..
!
!  Purpose
!  =======
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!

!     .. Local Scalars ..
      double precision STEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      STEMP = 0.0e0
      SDOT = 0.0e0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          STEMP = STEMP + SX(IX)*SY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END

!*************************************************************************

      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
!     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,6)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF (N.LT.6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END

!*************************************************************************

      SUBROUTINE DSCAL(N,DA,DX,INCX)
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!*
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END

!*************************************************************************

      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END

!*************************************************************************

      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END

!*************************************************************************

      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END

!*************************************************************************


!*************************************************************************
!
! access of unsymmetric banded matrix in linpack
!
! linear matrix must have the following dimension:
!		ndim = ( 3*m + 1 ) * n
!
!*************************************************************************

	subroutine lp_init_system(n,m,abd,b)

! initializes band matrix

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n], at return x

	integer nb,i

	nb = ( 3*m + 1 ) * n

	do i=1,nb
	  abd(i) = 0.
	end do

	do i=1,n
	  b(i) = 0.
	end do

	end

!*************************************************************************

	subroutine lp_solve_system(n,m,abd,b,ipvt,z)

! solves system a*x=b

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)		!pivot information [n]
	double precision z(1)		!aux vector [n]

	integer lda,ml,mu,job,info
	double precision rcond

	info = 0
	rcond = 1.

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	!call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)	!rcond should be > 0
        call sgbfa(abd,lda,n,ml,mu,ipvt,info)		!info should be 0
	call sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	if( info .ne. 0 .or. rcond .eq. 0. ) then
	  write(6,*) 'condition number: ',rcond,info
	  stop 'error stop lp_solve_system: info'
	end if

	end

!*************************************************************************

	subroutine lp_subst_system(n,m,abd,b,ipvt)

! solves a*x=b with a already factored (by lp_solve_system)

	integer n,m
	double precision abd(1)		!band matrix (already factored)
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)		!pivot information [n]

	integer lda,ml,mu,job

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	end

!*************************************************************************

	subroutine lp_mult_band(n,m,abd,b,res)

! multiplies band matrix a with vector b and returns result in res

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n]
	double precision res(1)		!result [n]

	integer i,j,k,mm,j1,j2
	double precision acu

	mm = 2 * m + 1

	do i=1,n
	  j1 = max(1,i-m)
	  j2 = min(n,i+m)
	  acu = 0.
	  do j=j1,j2
            k = loclp(i,j,n,m)
	    !write(6,*) i,j,k,j1,j2,abd(k),b(j)
	    acu = acu + abd(k) * b(j)
	  end do
	  res(i) = acu
	end do

	end

!*************************************************************************

        subroutine dlp_init_system(n,m,abd,b)

! initializes band matrix

        implicit none

        integer n,m
        double precision abd(1)      !band matrix
        double precision b(1)        !right hand side [n], at return x

        integer nb,i
	double precision zero

        nb = ( 3*m + 1 ) * n
	zero = 0.

        do i=1,nb
          abd(i) = zero
        end do

        do i=1,n
          b(i) = zero
        end do

        end

!*************************************************************************

	subroutine dlp_solve_system(n,m,abd,b,ipvt,z)

! solves system a*x=b

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)			!pivot information [n]
	double precision z(1)		!aux vector [n]

	integer lda,ml,mu,job,info
	double precision rcond

        info = 0
        rcond = 1.

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)	!rcond should be > 0
        !call dgbfa(abd,lda,n,ml,mu,ipvt,info)          !info should be 0
	call dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	if( info .ne. 0 .or. rcond .eq. 0. ) then
	  write(6,*) 'condition number: ',rcond,info
	  stop 'error stop dlp_solve_system: info'
	end if

	end

!*************************************************************************

	subroutine dlp_subst_system(n,m,abd,b,ipvt)

! solves a*x=b with a already factored (by dlp_solve_system)

	integer n,m
	double precision abd(1)		!band matrix (already factored)
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)			!pivot information [n]

	integer lda,ml,mu,job

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	end

!*************************************************************************

	subroutine dlp_mult_band(n,m,abd,b,res)

! multiplies band matrix a with vector b and returns result in res

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n]
	double precision res(1)		!result [n]

	integer i,j,k,mm,j1,j2
	double precision acu

	mm = 2 * m + 1

	do i=1,n
	  j1 = max(1,i-m)
	  j2 = min(n,i+m)
	  acu = 0.
	  do j=j1,j2
            k = loclp(i,j,n,m)
	    !write(6,*) i,j,k,j1,j2,abd(k),b(j)
	    acu = acu + abd(k) * b(j)
	  end do
	  res(i) = acu
	end do

	end

!*************************************************************************

        function loclp(i,j,n,m)

! access linpack routines (unsymmetric banded matrix)
!
! (i,j)   position of element in square matrix (row,column)
! n       dimension of square matrix
! m       band width of square matrix
! loclp   position of element in band matrix

        implicit none

        integer loclp
        integer i,j,n,m

	integer lda,k

	loclp = 0
	if( abs(i-j) .gt. m ) return

	lda = 3*m + 1
	k = i - j + lda - m

	loclp = k + (j-1) * lda

        end

!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************

	subroutine loclp_test

! use test of linpack routine

	implicit none

	integer n,m,lda,lnmax
	parameter (n=6,m=2,lda=3*m+1,lnmax=lda)

	integer i,j,k
	integer i1,i2
	double precision val
	double precision a(n,n)
	double precision b(lnmax*lnmax)
	double precision c(lda,n)

	double precision br(lnmax)
	double precision z(lnmax)
	integer ipvt(lnmax)

	write(6,*) 'n,m,lda,lnmax: ',n,m,lda,lnmax
	do j=1,n
	  do i=1,n
	    if( j-i .gt. m ) then
		val = 0.
	    else if( i-j .gt. m-1 ) then
		val = 0.
	    else
		val = 10.*i + j
	    end if
	    a(i,j) = val
	  end do
	end do

	write(6,*) 'matrix a'
	call loclp_print(n,n,a)
	    
	call loclp_init(lnmax,lnmax,b)
	write(6,*) 'matrix b init'
	call loclp_print(lnmax,lnmax,b)

	do j=1,n
	  do i=1,n
	    k = loclp(i,j,n,m)
	    b(k) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix b'
	call loclp_print(lda,n,b)
	    
	write(6,*) 'matrix b through loclp'
	do i=1,n 
	  do j=1,n
	    k = loclp(i,j,n,m)
	    val = 0.
	    if( k .gt. 0 ) val = b(k)
	    br(j) = val
	  end do
	  write(6,*) (br(j),j=1,n)
	end do

	call loclp_init(lda,n,c)

	do j=1,n
	  i1 = max(1,j-m)
	  i2 = min(n,j+m)
	  do i=i1,i2
	    k = i - j + 2*m + 1
	    c(k,j) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix c'
	call loclp_print(lda,n,c)
	    
	do i=1,n
	  z(i) = i
	end do

	call lp_mult_band(n,m,b,z,br)

	write(6,*) 'right hand side'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	call lp_solve_system(n,m,b,br,ipvt,z)

	write(6,*) 'solution system'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	end

!*************************************************************************

	subroutine dloclp_test

! use test of linpack routine

	implicit none

	integer n,m,lda,lnmax
	parameter (n=6,m=2,lda=3*m+1,lnmax=lda)

	integer i,j,k
	integer i1,i2
	double precision val
	double precision a(n,n)
	double precision b(lnmax*lnmax)
	double precision c(lda,n)

	double precision br(lnmax)
	double precision z(lnmax)
	integer ipvt(lnmax)

	write(6,*) 'n,m,lda,lnmax: ',n,m,lda,lnmax
	do j=1,n
	  do i=1,n
	    if( j-i .gt. m ) then
		val = 0.
	    else if( i-j .gt. m-1 ) then
		val = 0.
	    else
		val = 10.*i + j
	    end if
	    a(i,j) = val
	  end do
	end do

	write(6,*) 'matrix a'
	call dloclp_print(n,n,a)
	    
	call dloclp_init(lnmax,lnmax,b)
	write(6,*) 'matrix b init'
	call dloclp_print(lnmax,lnmax,b)

	do j=1,n
	  do i=1,n
	    k = loclp(i,j,n,m)
	    b(k) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix b'
	call dloclp_print(lda,n,b)
	    
	write(6,*) 'matrix b through loclp'
	do i=1,n 
	  do j=1,n
	    k = loclp(i,j,n,m)
	    val = 0.
	    if( k .gt. 0 ) val = b(k)
	    br(j) = val
	  end do
	  write(6,*) (br(j),j=1,n)
	end do

	call dloclp_init(lda,n,c)

	do j=1,n
	  i1 = max(1,j-m)
	  i2 = min(n,j+m)
	  do i=i1,i2
	    k = i - j + 2*m + 1
	    c(k,j) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix c'
	call dloclp_print(lda,n,c)
	    
	do i=1,n
	  z(i) = i
	end do

	call dlp_mult_band(n,m,b,z,br)

	write(6,*) 'right hand side'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	call dlp_solve_system(n,m,b,br,ipvt,z)

	write(6,*) 'solution system'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	end

!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************

	subroutine loclp_init(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do j=1,n
	  do i=1,l
	    a(i,j) = 0.
	  end do
	end do

	end

!*************************************************************************

	subroutine loclp_print(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do i=1,l
	  write(6,*) (a(i,j),j=1,n)
	end do

	end

!*************************************************************************

	subroutine dloclp_init(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do j=1,n
	  do i=1,l
	    a(i,j) = 0.
	  end do
	end do

	end

!*************************************************************************

	subroutine dloclp_print(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do i=1,l
	  write(6,*) (a(i,j),j=1,n)
	end do

	end

!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************

!	program loclp_main
!	call loclp_test
!	call dloclp_test
!	end

!*************************************************************************

!---------------------------------------------------------------------------------
      end module matrix_inv
!---------------------------------------------------------------------------------
