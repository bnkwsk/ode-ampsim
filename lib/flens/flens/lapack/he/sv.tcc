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
       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
     $                  LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_HE_SV_TCC
#define FLENS_LAPACK_HE_SV_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

//-- (he)sv [complex variant] --------------------------------------------------

template <typename MA, typename VPIV, typename MB, typename VWORK>
typename HeMatrix<MA>::IndexType
sv_impl(HeMatrix<MA> &A, DenseVector<VPIV> &piv, GeMatrix<MB> &B,
        DenseVector<VWORK> &work)
{
    typedef typename HeMatrix<MA>::IndexType    IndexType;
    typedef typename HeMatrix<MA>::ElementType  ElementType;

    if (work.length()==0) {
        ElementType     WORK;
        IndexType       LWORK = -1;

        cxxlapack::hesv<IndexType>(getF77Char(A.upLo()),
                                   A.dim(),
                                   B.numCols(),
                                   A.data(),
                                   A.leadingDimension(),
                                   piv.data(),
                                   B.data(),
                                   B.leadingDimension(),
                                   &WORK,
                                   LWORK);
        work.resize(cxxblas::real(WORK));
    }

    IndexType info = cxxlapack::hesv<IndexType>(getF77Char(A.upLo()),
                                                A.dim(),
                                                B.numCols(),
                                                A.data(),
                                                A.leadingDimension(),
                                                piv.data(),
                                                B.data(),
                                                B.leadingDimension(),
                                                work.data(),
                                                work.length());
    ASSERT(info>=0);
    return info;
}

} // namespace external


//== public interface ==========================================================

//-- (he)sv [complex variant] --------------------------------------------------

template <typename MA, typename VPIV, typename MB, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsComplexGeMatrix<MB>::value
                 && IsComplexDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, MB &&B, VWORK && work)
{
    LAPACK_DEBUG_OUT("(he)sv [complex]");

//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename MatrixA::IndexType     IndexType;

    if (piv.length()==0) {
        piv.resize(A.dim());
    }
    ASSERT(piv.length()==0);
//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);
    ASSERT(piv.stride()==1);
    ASSERT(B.firstRow()==1);
    ASSERT(B.firstCol()==1);
    ASSERT(B.numRows()==A.dim());
#   endif
//
//  Call implementation
//
    IndexType info = external::sv_impl(A, piv, B, work);

    return info;
}

//-- (he)sv [complex variant, if rhs is vector] --------------------------------

template <typename MA, typename VPIV, typename VB, typename VWORK>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsComplexDenseVector<VB>::value
                 && IsComplexDenseVector<VWORK>::value,
          typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, VB &&b, VWORK &&work)
{
//
//  Remove references from rvalue types
//
    typedef typename RemoveRef<MA>::Type    MatrixA;
    typedef typename RemoveRef<VB>::Type    VectorB;

    typedef typename VectorB::ElementType  ElementType;
    typedef typename VectorB::IndexType    IndexType;

    const IndexType    n     = b.length();
    const StorageOrder order = MatrixA::Engine::order;

    GeMatrix<FullStorageView<ElementType, order> >  B(n, 1, b, n);

    return sv(A, piv, B, work);
}

//-- (he)sv [complex variant with temporary workspace] -------------------------

template <typename MA, typename VPIV, typename VB>
typename RestrictTo<IsHeMatrix<MA>::value
                 && IsIntegerDenseVector<VPIV>::value
                 && IsComplexDenseVector<VB>::value,
         typename RemoveRef<MA>::Type::IndexType>::Type
sv(MA &&A, VPIV &&piv, VB &&b)
{
    typedef typename RemoveRef<MA>::Type::Vector WorkVector;

    WorkVector  work;

    return sv(A, piv, b, work);
}

#endif // USE_CXXLAPACK

} } // namespace lapack, flens

#endif // FLENS_LAPACK_HE_SV_TCC
