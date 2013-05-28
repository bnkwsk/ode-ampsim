/*
 *   Copyright (c) 2007, Michael Lehn
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

#ifndef FLENS_MATRIXTYPES_IMPL_MATRIXCLOSURE_TCC
#define FLENS_MATRIXTYPES_IMPL_MATRIXCLOSURE_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/matrixtypes/impl/matrixclosure.h>

namespace flens {

template <typename Op, typename L, typename R>
MatrixClosure<Op, L, R>::MatrixClosure(typename ConstRef<L>::Type l,
                                       typename ConstRef<R>::Type r)
    : _left(l), _right(r)
{
}

template <typename Op, typename L, typename R>
typename ConstRef<L>::Type
MatrixClosure<Op, L, R>::left() const
{
    return _left;
}

template <typename Op, typename L, typename R>
typename ConstRef<R>::Type
MatrixClosure<Op, L, R>::right() const
{
    return _right;
}

} // namespace flens

#endif // FLENS_MATRIXTYPES_IMPL_MATRIXCLOSURE_TCC
