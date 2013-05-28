#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dpocon --------------------------------------------------------------------
void
LAPACK_DECL(dpocon)(const char       *UPLO,
                    const INTEGER    *N,
                    const DOUBLE     *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *ANORM,
                    DOUBLE           *RCOND,
                    DOUBLE           *WORK,
                    INTEGER          *IWORK,
                    INTEGER          *INFO)
{
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    } else if (*ANORM<0) {
        *INFO = -5;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DPOCON", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo         upLo = StorageUpLo(*UPLO);
    DConstFSView        AFS  = DConstFSView(*N, *N, A, *LDA);
    DConstSyMatrixView  _A   = DConstSyMatrixView(AFS, upLo);
    DDenseVectorView    work  = DArrayView(*N*3, WORK, 1);
    IDenseVectorView    iwork = IArrayView(*N, IWORK, 1);

    pocon(_A, *ANORM, *RCOND, work, iwork);
}

//-- zpocon --------------------------------------------------------------------
void
LAPACK_DECL(zpocon)(const char               *UPLO,
                    const INTEGER            *N,
                    const DOUBLE_COMPLEX     *A,
                    const INTEGER            *LDA,
                    const DOUBLE             *ANORM,
                    DOUBLE                   *RCOND,
                    DOUBLE_COMPLEX           *WORK,
                    DOUBLE                   *RWORK,
                    INTEGER                  *INFO)
{
//
//  Test the input parameters so that we pass LAPACK error checks
//
    *INFO = 0;
    if (*UPLO!='U' && *UPLO!='L') {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INTEGER(1), *N)) {
        *INFO = -4;
    } else if (*ANORM<0) {
        *INFO = -5;
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZPOCON", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    const auto *zA  = reinterpret_cast<const CXX_DOUBLE_COMPLEX *>(A);
    auto *zWORK     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(WORK);

    StorageUpLo         upLo  = StorageUpLo(*UPLO);
    ZConstFSView        AFS   = ZConstFSView(*N, *N, zA, *LDA);
    ZConstHeMatrixView  _A    = ZConstHeMatrixView(AFS, upLo);
    ZDenseVectorView    work  = ZArrayView(*N*2, zWORK, 1);
    DDenseVectorView    rwork = DArrayView(*N, RWORK, 1);

    pocon(_A, *ANORM, *RCOND, work, rwork);
}

} // extern "C"

} } // namespace lapack, flens
