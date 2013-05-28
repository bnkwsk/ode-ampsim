/*
 *   Copyright (c) 2011, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Based on
 *
      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
 *
 *  -- LAPACK routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 */

#ifndef FLENS_LAPACK_GE_BAL_TCC
#define FLENS_LAPACK_GE_BAL_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename MA, typename IndexType, typename VSCALE>
void
bal_impl(BALANCE::Balance    job,
         GeMatrix<MA>        &A,
         IndexType           &iLo,
         IndexType           &iHi,
         DenseVector<VSCALE> &scale)
{
    using namespace BALANCE;

    using std::abs;
    using std::isnan;
    using flens::max;
    using flens::min;

    typedef typename GeMatrix<MA>::ElementType  T;

    const T Zero(0), One(1);
    const T scaleFactor(2), factor(0.95);

    const T safeMin1 = lamch<T>(SafeMin) / lamch<T>(Precision);
    const T safeMax1 = One / safeMin1;
    const T safeMin2 = safeMin1*scaleFactor;
    const T safeMax2 = One / safeMin2;


    const Underscore<IndexType> _;
    const IndexType n = A.numRows();

    IndexType k = 1;
    IndexType l = n;

    IndexType j = l, m;

    if (n==0) {
        goto DONE;
    }

    if (job==None) {
        for (IndexType i=1; i<=n; ++i) {
            scale(i) = One;
        }
        goto DONE;
    }

    if (job!=ScaleOnly) {
//
//      Permutation to isolate eigenvalues if possible
//
        IndexType iExc = 0;
//
//      Row and column exchange.
//
        EXCHANGE:
            if (iExc!=0) {
                scale(m) = j;
                if (j!=m) {
                    blas::swap(A(_(1,l),j), A(_(1,l),m));
                    blas::swap(A(j,_(k,n)), A(m,_(k,n)));
                }
            }

        switch (iExc) {
//
//      Search for rows isolating an eigenvalue and push them down.
//
        case 1:
            if (l==1) {
                goto DONE;
            }
            --l;

        case 0:
            for (j=l; j>=1; --j) {
                bool foundRow = true;

                for (IndexType i=1; i<=l; ++i) {
                    if (i==j) {
                        continue;
                    }
                    if (A(j,i)!=Zero) {
                        foundRow = false;
                        break;
                    }
                }

                if (foundRow) {
                    m = l;
                    iExc = 1;
                    goto EXCHANGE;
                }
            }
//
//      Search for columns isolating an eigenvalue and push them left.
//
        case 2:
            if ((iExc!=0) && (iExc!=1)) {
                ++k;
            }

            for (j=k; j<=l; ++j) {
                bool foundCol = true;

                for (IndexType i=k; i<=l; ++i) {
                    if (i==j) {
                        continue;
                    }
                    if (A(i,j)!=Zero) {
                        foundCol = false;
                        break;
                    }
                }

                if (foundCol) {
                    m = k;
                    iExc = 2;
                    goto EXCHANGE;
                }
            }
        }
    }

    for (IndexType i=k; i<=l; ++i) {
        scale(i) = One;
    }

    if (job==PermuteOnly) {
        goto DONE;
    }
//
//  Balance the submatrix in rows K to L.
//
//  Iterative loop for norm reduction
//
    bool noConv;
    do {
        noConv = false;

        for (IndexType i=k; i<=l; ++i) {
            T c = Zero;
            T r = Zero;

            for (IndexType j=k; j<=l; ++j) {
                if (j==i) {
                    continue;
                }
                c += abs(A(j,i));
                r += abs(A(i,j));
            }
            const IndexType ica = blas::iamax(A(_(1,l),i));
            T ca = abs(A(ica,i));
            const IndexType ira = blas::iamax(A(i,_(k,n)))+k-1;
            T ra = abs(A(i,ira));
//
//          Guard against zero C or R due to underflow.
//
            if (c==Zero || r==Zero) {
                continue;
            }
            T g = r / scaleFactor;
            T f = One;
            T s = c + r;

            while (c<g && max(f,c,ca)<safeMax2 && min(r,g,ra)>safeMin2) {
                if (isnan(c+f+ca+r+g+ra)) {
//
//                  Exit if NaN to avoid infinite loop
//
                    ASSERT(0);
                    return;
                }
                f *= scaleFactor;
                c *= scaleFactor;
                ca *= scaleFactor;
                r /= scaleFactor;
                g /= scaleFactor;
                ra /= scaleFactor;
            }

            g = c / scaleFactor;
            while (g>=r && max(r,ra)<safeMax2 && min(f,c,g,ca)>safeMin2) {
                f /= scaleFactor;
                c /= scaleFactor;
                g /= scaleFactor;
                ca /= scaleFactor;
                r *= scaleFactor;
                ra *= scaleFactor;
            }
//
//          Now balance.
//
            if ((c+r)>=factor*s) {
                continue;
            }
            if (f<One && scale(i)<One) {
                if (f*scale(i)<=safeMin1) {
                    continue;
                }
            }
            if (f>One && scale(i)>One) {
                if (scale(i)>=safeMax1/f) {
                    continue;
                }
            }
            g = One / f;
            scale(i) *= f;
            noConv = true;

            A(i,_(k,n)) *= g;
            A(_(1,l),i) *= f;

        }

    } while (noConv);

    DONE:
        iLo = k;
        iHi = l;
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename IndexType, typename VSCALE>
void
bal_impl(BALANCE::Balance    job,
         GeMatrix<MA>        &A,
         IndexType           &iLo,
         IndexType           &iHi,
         DenseVector<VSCALE> &scale)
{
    cxxlapack::gebal<IndexType>(getF77Char(job),
                                A.numRows(),
                                A.data(),
                                A.leadingDimension(),
                                iLo,
                                iHi,
                                scale.data());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

template <typename MA, typename IndexType, typename VSCALE>
typename RestrictTo<IsRealGeMatrix<MA>::value
                 && IsInteger<IndexType>::value
                 && IsRealDenseVector<VSCALE>::value,
         void>::Type
bal(BALANCE::Balance    job,
    MA                  &&A,
    IndexType           &iLo,
    IndexType           &iHi,
    VSCALE              &&scale)
{
    LAPACK_DEBUG_OUT("bal");

//
//  Remove references from the types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename RemoveRef<VSCALE>::Type  VectorScale;


#   ifndef NDEBUG
//
//  Test the input parameters
//
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(A.numRows()==A.numCols());
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename MatrixA::NoView        A_org      = A;
    IndexType                       iLo_org    = iLo;
    IndexType                       iHi_org    = iHi;
    typename VectorScale::NoView    scale_org  = scale;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::bal_impl(job, A, iLo, iHi, scale);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename MatrixA::NoView        A_generic      = A;
    IndexType                       iLo_generic    = iLo;
    IndexType                       iHi_generic    = iHi;
    typename VectorScale::NoView    scale_generic  = scale;

    A = A_org;
    iLo = iLo_org;
    iHi = iHi_org;
    scale = scale_org;

    external::bal_impl(job, A, iLo, iHi, scale);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }

    if (! isIdentical(iLo_generic, iLo, "iLo_generic", "iLo")) {
        failed = true;
    }

    if (! isIdentical(iHi_generic, iHi, "iHi_generic", "iHi")) {
        failed = true;
    }

    if (! isIdentical(scale_generic, scale, "scale_generic", "scale")) {
        std::cerr << "CXXLAPACK: scale_generic = "
                  << scale_generic << std::endl;
        std::cerr << "F77LAPACK: scale = "
                  << scale << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_BAL_TCC
