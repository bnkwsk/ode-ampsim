#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>

namespace flens { namespace lapack {

extern "C" {

//-- dpotri --------------------------------------------------------------------
void
LAPACK_DECL(dpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
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
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("DPOTRI", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    StorageUpLo    upLo = StorageUpLo(*UPLO);
    DFSView        AFS  = DFSView(*N, *N, A, *LDA);
    DSyMatrixView  _A   = DSyMatrixView(AFS, upLo);

    *INFO = potri(_A);
}

//-- zpotri --------------------------------------------------------------------
void
LAPACK_DECL(zpotri)(const char       *UPLO,
                    const INTEGER    *N,
                    DOUBLE_COMPLEX   *A,
                    const INTEGER    *LDA,
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
    }
    if (*INFO!=0) {
        *INFO = -(*INFO);
        LAPACK_ERROR("ZPOTRI", INFO);
        *INFO = -(*INFO);
        return;
    }
//
//  Call FLENS implementation
//
    auto zA     = reinterpret_cast<CXX_DOUBLE_COMPLEX *>(A);

    StorageUpLo    upLo = StorageUpLo(*UPLO);
    ZFSView        AFS  = ZFSView(*N, *N, zA, *LDA);
    ZHeMatrixView  _A   = ZHeMatrixView(AFS, upLo);

    *INFO = potri(_A);
}


} // extern "C"

} } // namespace lapack, flens
