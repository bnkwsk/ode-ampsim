      SUBROUTINE CLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )
*
*     -- LAPACK routine (version 3.2.2)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- June 2010                                                    --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            N, NZ, NRHS
*     ..
*     .. Array Arguments ..
      REAL               AYB( N, NRHS ), BERR( NRHS )
      COMPLEX            RES( N, NRHS )
*     ..
*
*  Purpose
*  =======
*
*     CLA_LIN_BERR computes componentwise relative backward error from
*     the formula
*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
*     where abs(Z) is the componentwise absolute value of the matrix
*     or vector Z.
*
*     N       (input) INTEGER
*     The number of linear equations, i.e., the order of the
*     matrix A.  N >= 0.
*
*     NZ      (input) INTEGER
*     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to
*     guard against spuriously zero residuals. Default value is N.
*
*     NRHS    (input) INTEGER
*     The number of right hand sides, i.e., the number of columns
*     of the matrices AYB, RES, and BERR.  NRHS >= 0.
*
*     RES    (input) DOUBLE PRECISION array, dimension (N,NRHS)
*     The residual matrix, i.e., the matrix R in the relative backward
*     error formula above.
*
*     AYB    (input) DOUBLE PRECISION array, dimension (N, NRHS)
*     The denominator in the relative backward error formula above, i.e.,
*     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B
*     are from iterative refinement (see cla_gerfsx_extended.f).
*     
*     BERR   (output) COMPLEX array, dimension (NRHS)
*     The componentwise relative backward error from the formula above.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               TMP
      INTEGER            I, J
      COMPLEX            CDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, AIMAG, MAX
*     ..
*     .. External Functions ..
      EXTERNAL           SLAMCH
      REAL               SLAMCH
      REAL               SAFE1
*     ..
*     .. Statement Functions ..
      COMPLEX            CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Adding SAFE1 to the numerator guards against spuriously zero
*     residuals.  A similar safeguard is in the CLA_yyAMV routine used
*     to compute AYB.
*
      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (NZ+1)*SAFE1

      DO J = 1, NRHS
         BERR(J) = 0.0
         DO I = 1, N
            IF (AYB(I,J) .NE. 0.0) THEN
               TMP = (SAFE1 + CABS1(RES(I,J)))/AYB(I,J)
               BERR(J) = MAX( BERR(J), TMP )
            END IF
*
*     If AYB is exactly 0.0 (and if computed by CLA_yyAMV), then we know
*     the true residual also must be exactly 0.0.
*
         END DO
      END DO
      END
