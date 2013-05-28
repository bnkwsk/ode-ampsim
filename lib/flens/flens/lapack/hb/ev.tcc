/*
 *   Copyright (c) 2012, Klaus Pototzky
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
      SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
     $                  RWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.3.1) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *  -- April 2011                                                      --
 */

#ifndef FLENS_LAPACK_HB_EV_TCC
#define FLENS_LAPACK_HB_EV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (hb)ev [real variant] -----------------------------------------------------

template <typename MA, typename VW, typename MZ,
          typename VWORK, typename VRWORK>
typename HbMatrix<MA>::IndexType
ev_impl(bool                  computeZ,
        HbMatrix<MA>          &A,
        DenseVector<VW>       &w,
        GeMatrix<MZ>          &Z,
        DenseVector<VWORK>    &work,
        DenseVector<VRWORK>   &rwork)
{
    using std::max;

    typedef typename HbMatrix<MA>::IndexType  IndexType;

    if (work.length()==0) {
        work.resize(A.dim(), 1);
    }
    if (rwork.length()==0) {
        rwork.resize(max(1, 3*A.dim()-2));
    }

    IndexType  info;
    info = cxxlapack::hbev(computeZ ? 'V' : 'N',
                           getF77Char(A.upLo()),
                           A.dim(),
                           A.numOffDiags(),
                           A.data(),
                           A.leadingDimension(),
                           w.data(),
                           Z.data(),
                           Z.leadingDimension(),
                           work.data(),
                           rwork.data());
    ASSERT(info>=0);
    return info;
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================

#ifdef USE_CXXLAPACK

//-- (hb)ev [complex variant] -----------------------------------------------------
template <typename MA, typename VW, typename MZ, typename VWORK,
          typename VRWORK>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsRealDenseVector<VW>::value
                 && IsComplexGeMatrix<MZ>::value
                 && IsComplexDenseVector<VWORK>::value
                 && IsRealDenseVector<VRWORK>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeZ,
   MA       &&A,
   VW       &&w,
   MZ       &&Z,
   VWORK    &&work,
   VRWORK   &&rWork)
{
    LAPACK_DEBUG_OUT("(hb)ev [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type      MatrixA;
    typedef typename MatrixA::IndexType       IndexType;
    typedef typename RemoveRef<VW>::Type      VectorW;
    typedef typename RemoveRef<MZ>::Type      MatrixZ;
    typedef typename RemoveRef<VWORK>::Type   VectorWork;
    typedef typename RemoveRef<VRWORK>::Type  VectorRWork;

    const IndexType n = A.numRows();

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstIndex()==1);
    ASSERT(work.firstIndex()==1);
    ASSERT(rWork.firstIndex()==1);

    ASSERT(w.firstIndex()==1);
    ASSERT(w.length()==0 || w.length()==n);

    if (computeZ) {
        ASSERT(Z.numRows()==Z.numCols());
        ASSERT(Z.numRows()==0 || Z.numRows()==n);
        ASSERT(Z.firstRow()==1);
        ASSERT(Z.firstCol()==1);
    }
#   endif

//
//  Resize output arguments if they are empty and needed
//
    if (w.length()==0) {
        w.resize(n, 1);
    }
    if (computeZ && Z.numRows()==0) {
        Z.resize(n, n, 1, 1);
    }

//
//  Call external implementation
//
    IndexType result = external::ev_impl(computeZ, A, w, Z,
                                         work, rWork);
    return result;
}

//-- (hb)ev [complex variant with temporary workspace] -------------------------

template <typename MA, typename VW, typename MZ>
typename RestrictTo<IsHbMatrix<MA>::value
                 && IsRealDenseVector<VW>::value
                 && IsComplexGeMatrix<MZ>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
ev(bool     computeZ,
   MA       &&A,
   VW       &&w,
   MZ       &&Z)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type::Vector        WorkVector;
    typedef typename RemoveRef<MA>::Type::ElementType   T;
    typedef typename ComplexTrait<T>::PrimitiveType     PT;
    typedef DenseVector<Array<PT> >                     RealWorkVector;

    WorkVector      work;
    RealWorkVector  rwork;

    return ev(computeZ, A, w, Z, work, rwork);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HB_EV_TCC
