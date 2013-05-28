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
       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006

 and

       DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
      $                 WORK )
 *
 *  -- LAPACK auxiliary routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LAN_H
#define FLENS_LAPACK_LA_LAN_H 1

#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== lan(ge) ===================================================================
template <typename MA, typename VWORK>
    typename ComplexTrait<typename GeMatrix<MA>::ElementType>::PrimitiveType
    lan(Norm norm, const GeMatrix<MA> &A);

template <typename MA, typename VWORK>
    typename ComplexTrait<typename GeMatrix<MA>::ElementType>::PrimitiveType
    lan(Norm norm, const GeMatrix<MA> &A, DenseVector<VWORK> &work);

//== lan(tr) ===================================================================
template <typename MA, typename VWORK>
    typename ComplexTrait<typename TrMatrix<MA>::ElementType>::PrimitiveType
    lan(Norm norm, const TrMatrix<MA> &A);

template <typename MA, typename VWORK>
    typename ComplexTrait<typename TrMatrix<MA>::ElementType>::PrimitiveType
    lan(Norm norm, const TrMatrix<MA> &A, DenseVector<VWORK> &work);

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VWORK>
    typename ComplexTrait<typename MA::ElementType>::PrimitiveType
    lan(Norm norm, const MA &A);

template <typename MA, typename VWORK>
    typename ComplexTrait<typename MA::ElementType>::PrimitiveType
    lan(Norm norm, const MA &A, VWORK &&work);

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LAN_H
