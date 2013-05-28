/*
 *   Copyright (c) 2007-2012, Michael Lehn
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

#ifndef FLENS_STORAGE_CRS_CRS_TCC
#define FLENS_STORAGE_CRS_CRS_TCC 1

#include <vector>
#include <flens/storage/crs/crs.h>

namespace flens {

template <typename T, typename I>
CRS<T,I>::CRS()
    : _numRows(0), _numCols(0),
      _indexBase(I::defaultIndexBase)
{
}

template <typename T, typename I>
CRS<T,I>::~CRS()
{
}

//-- operators -----------------------------------------------------------------

template <typename T, typename I>
template <typename T2, typename I2>
void
CRS<T,I>::operator=(const CoordStorage<T2, CoordRowColCmp, I2> &coordStorage)
{
    _compress(coordStorage);
}

//-- methods -------------------------------------------------------------------

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::indexBase() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::firstRow() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::lastRow() const
{
    return _indexBase+_numRows-1;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::firstCol() const
{
    return _indexBase;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::lastCol() const
{
    return _indexBase+_numCols-1;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::numRows() const
{
    return _numRows;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::numCols() const
{
    return _numCols;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexType
CRS<T,I>::numNonZeros() const
{
    return _values.length();
}

template <typename T, typename I>
const typename CRS<T,I>::IndexTypeVector &
CRS<T,I>::rows() const
{
    return _rows;
}

template <typename T, typename I>
typename CRS<T,I>::IndexTypeVector &
CRS<T,I>::rows()
{
    return _rows;
}

template <typename T, typename I>
const typename CRS<T,I>::IndexTypeVector &
CRS<T,I>::cols() const
{
    return _cols;
}

template <typename T, typename I>
typename CRS<T,I>::IndexTypeVector &
CRS<T,I>::cols()
{
    return _cols;
}

template <typename T, typename I>
const typename CRS<T,I>::ElementTypeVector &
CRS<T,I>::values() const
{
    return _values;
}

template <typename T, typename I>
typename CRS<T,I>::ElementTypeVector &
CRS<T,I>::values()
{
    return _values;
}

template <typename T, typename I>
template <typename T2, typename I2>
void
CRS<T,I>::_compress(const CoordStorage<T2, CoordRowColCmp, I2> &coordStorage)
{
    _numRows  = coordStorage.numRows();
    _numCols  = coordStorage.numCols();
    _indexBase = coordStorage.indexBase();

//
//  Accumulate coords and get number of non zeros
//
    coordStorage.accumulate();
    IndexType nnz = coordStorage.numNonZeros();


//
//  Allocate memory for the CRS format
//
    _rows.resize(_numRows+1, _indexBase);
    _cols.resize(nnz, _indexBase);
    _values.resize(nnz, _indexBase);

    const auto &coord = coordStorage.coordVector();

    IndexType r = _indexBase;
    _rows(r) = _indexBase;

    for (size_t k=0; k<coord.size(); ++k) {
        while (coord[k].row>r) {
            _rows(r+1) = _indexBase + k;
            ++r;
        }
        _cols(_indexBase+k)   = coord[k].col;
        _values(_indexBase+k) = coord[k].value;
    }
    while (r<=lastRow()) {
        _rows(r+1) = _indexBase + coord.size();
        ++r;
    }
}

} // namespace flens

#endif // FLENS_STORAGE_CRS_CRS_TCC
