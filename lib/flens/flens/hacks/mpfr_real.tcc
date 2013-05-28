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

#if defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_TCC)
#define FLENS_HACKS_MPFR_REAL_TCC 1

/*
 * NOTE: This hack requires that in the mpfr::real class the '_x' attribute
 *       is made public!
 */

namespace mpfr {

//
// aux-functions used in numeric_limits
//

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
const mpfr::real<_prec,_rnd>
nextabove(const mpfr::real<_prec,_rnd> &x)
{
    mpfr::real<_prec,_rnd> tmp = x;
    mpfr_nextabove(tmp._x);
    return tmp;
}

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
const mpfr::real<_prec,_rnd>
get_max()
{
    const unsigned long emax = mpfr_get_emax();
    mpfr::real<_prec,_rnd> tmp(1);

    mpfr_mul_2ui(tmp._x, tmp._x, emax-1, _rnd);
    return tmp;
}

} // namespace mpfr

//
// numeric_limits
//

namespace std {

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    const mpfr::real<_prec,_rnd>
    numeric_limits<mpfr::real<_prec,_rnd> >::_max
        = mpfr::get_max<_prec,_rnd>();

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    const mpfr::real<_prec,_rnd>
    numeric_limits<mpfr::real<_prec,_rnd> >::_min
        = mpfr::nextabove(T(0));

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    const mpfr::real<_prec,_rnd>
    numeric_limits<mpfr::real<_prec,_rnd> >::_eps
        = mpfr::nextabove(T(1)) - T(1);

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
const mpfr::real<_prec,_rnd>
numeric_limits<mpfr::real<_prec,_rnd> >::epsilon()
{
     return _eps;
}

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
const mpfr::real<_prec,_rnd>
numeric_limits<mpfr::real<_prec,_rnd> >::max()
{
    return _max;
}

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
const mpfr::real<_prec,_rnd>
numeric_limits<mpfr::real<_prec,_rnd> >::min()
{
    return _min;
}

//
// import isnan to namespace std
//
template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
bool
isnan(const mpfr::real<_prec,_rnd> &x)
{
    //
    // TODO: There seems to be a bug in clang++ because
    //
    //          return mpfr::isnan(x);
    //
    //       does not compile (but is fine with g++4.7)
    return mpfr_nan_p(x._x);
}

//
// import mpfr_real to namespace std
//
template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
int
signbit(const mpfr::real<_prec,_rnd> &x)
{
    //
    // TODO: There seems to be a bug in clang++ because
    //
    //          return mpfr::signbit(x);
    //
    //       does not compile (but is fine with g++4.7)
    return mpfr_signbit(x._x);
}


} // namespace std

#endif // defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_TCC)
