#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

float
BLAS(snrm2)(const INTEGER   *N,
            const float     *X,
            const INTEGER   *INCX)
{
    using std::abs;

    SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::nrm2(x);
}

double
BLAS(dnrm2)(const INTEGER   *N,
            const double    *X,
            const INTEGER   *INCX)
{
    using std::abs;

    DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::nrm2(x);
}

float
BLAS(scnrm2)(const INTEGER   *N,
             const cfloat    *X,
             const INTEGER   *INCX)
{
    using std::abs;

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::nrm2(x);
}

double
BLAS(dznrm2)(const INTEGER   *N,
            const cdouble   *X,
            const INTEGER   *INCX)
{
    using std::abs;

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

    return blas::nrm2(x);
}

} // extern "C"
