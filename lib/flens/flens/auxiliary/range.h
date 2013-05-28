/*
 *   Copyright (c) 2004, Alexander Stippler
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

#ifndef FLENS_AUXILIARY_RANGE_H
#define FLENS_AUXILIARY_RANGE_H 1

#include <iostream>

namespace flens {

template <typename IndexType = long>
class Range
{
    public:
        Range(IndexType firstIndex, IndexType lastIndex);

        Range(IndexType firstIndex, IndexType stride, IndexType lastIndex);

        IndexType
        firstIndex() const;

        IndexType
        stride() const;

        IndexType
        lastIndex() const;

        IndexType
        numTicks() const;

        IndexType
        length() const;

    private:
        IndexType _firstIndex, _stride, _lastIndex;
};

template <typename IndexType>
    std::ostream &
    operator<<(std::ostream &out, const Range<IndexType> &range);

template <typename IndexType>
    bool
    operator==(const Range<IndexType> &lhs, const Range<IndexType> &rhs);

} // namespace flens

#endif // FLENS_AUXILIARY_RANGE_H
