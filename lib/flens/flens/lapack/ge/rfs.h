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
       SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
      $                   X, LDX, FERR, BERR, WORK, IWORK, INFO )
       SUBROUTINE ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
      $                   X, LDX, FERR, BERR, WORK, RWORK, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_GE_RFS_H
#define FLENS_LAPACK_GE_RFS_H 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== (ge)rfs ===================================================================
//
//  Real variant
//
template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VIWORK>
    typename RestrictTo<IsRealGeMatrix<MA>::value
                     && IsRealGeMatrix<MAF>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsRealGeMatrix<MB>::value
                     && IsRealGeMatrix<MX>::value
                     && IsRealDenseVector<VFERR>::value
                     && IsRealDenseVector<VBERR>::value
                     && IsRealDenseVector<VWORK>::value
                     && IsIntegerDenseVector<VIWORK>::value,
             void>::Type
    rfs(Transpose   trans,
        const MA    &A,
        const MAF   &AF,
        const VPIV  &piv,
        const MB    &B,
        MX          &&X,
        VFERR       &&fErr,
        VBERR       &&bErr,
        VWORK       &&work,
        VIWORK      &&iwork);


#ifdef USE_CXXLAPACK
//
//  Complex variant
//
template <typename MA, typename MAF, typename VPIV, typename MB, typename MX,
          typename VFERR, typename VBERR, typename VWORK, typename VRWORK>
    typename RestrictTo<IsComplexGeMatrix<MA>::value
                     && IsComplexGeMatrix<MAF>::value
                     && IsIntegerDenseVector<VPIV>::value
                     && IsComplexGeMatrix<MB>::value
                     && IsComplexGeMatrix<MX>::value
                     && IsRealDenseVector<VFERR>::value
                     && IsRealDenseVector<VBERR>::value
                     && IsComplexDenseVector<VWORK>::value
                     && IsRealDenseVector<VRWORK>::value,
             void>::Type
    rfs(Transpose   trans,
        const MA    &A,
        const MAF   &AF,
        const VPIV  &piv,
        const MB    &B,
        MX          &&X,
        VFERR       &&fErr,
        VBERR       &&bErr,
        VWORK       &&work,
        VRWORK      &&rwork);

#endif // USE_CXXLAPACK


} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_RFS_H
