#include <complex>
#include <iostream>
#include <flens/flens.cxx>

#ifndef SEED
#define SEED  0
#endif

#ifndef MAX_M
#define MAX_M  1500
#endif

#ifndef MAX_N
#define MAX_N  1500
#endif

#ifndef MAX_NNZ
#define MAX_NNZ  3*MAX_M
#endif


using namespace flens;
using namespace std;

//
//  Create symmetric sparse matrix A and control matrix A_
//
template <typename CCS, typename FS>
void
setup(StorageUpLo upLo, int n, int max_nnz, int indexBase,
      HeCCSMatrix<CCS> &A, HeMatrix<FS> &A_)
{
    typedef typename HeCCSMatrix<CCS>::ElementType  ElementType;
    typedef CoordStorage<complex<double>, CoordColRowCmp>    Coord;

    const ElementType  Zero(0);

    A_.resize(n, upLo, indexBase);
    A_ = Zero;

    //
    //  We first setup the sparse matrix B in coordinate storage.  Later we
    //  convert it to compressed col storage.
    //
    HeCoordMatrix<Coord>  B(n, upLo, 1, indexBase);


    for (int k=1; k<=max_nnz; ++k) {
        const int i = indexBase + rand() % n;
        const int j = indexBase + rand() % n;
        const int v1r = rand() % 10;
        const int v1i = rand() % 10;
        const int v2r = rand() % 10;
        const int v2i = rand() % 10;

        const auto z1 = complex<double>(v1r,v1i);
        const auto z2 = complex<double>(v2r,v2i);

        if ((upLo==Upper && (j>=i)) || (upLo==Lower && (i>=j))) {

            B(i,j)  += z1;
            B(i,j)  -= z2;

            A_(i,j) += z1;
            A_(i,j) -= z2;
        } else {
            B(j,i)  += cxxblas::conjugate(z1);
            B(j,i)  -= cxxblas::conjugate(z2);

            A_(j,i) += cxxblas::conjugate(z1);
            A_(j,i) -= cxxblas::conjugate(z2);
        }
    }


    //
    //  Convert coordinate storage matrix B to compressed col storage matrix A
    //
    A = B;

    //
    //  Convert compressed col storage matrix A to full storage matrix A__
    //
    typename HeMatrix<FS>::NoView  A__ = A;

    if (! lapack::isIdentical(A_, A__, "A_", "A__")) {
        cerr << "n =       " << n << endl;
        cerr << "max_nnz = " << max_nnz << endl;
        ASSERT(0);
    }
}

template <typename CCS, typename FS>
void
mv(int n, int max_nnz, const HeCCSMatrix<CCS> &A, const HeMatrix<FS> &A_)
{
    typedef typename HeCCSMatrix<CCS>::ElementType  ElementType;

    DenseVector<Array<ElementType> >  x(n), y, y_;

    for (int j=1; j<=n; ++j) {
        x(j) = rand() % 10;
    }

    y  = A  * x;
    y_ = A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y = A*x" << endl;
        ASSERT(0);
    }

    y  += A  * x;
    y_ += A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y += A*x" << endl;
        ASSERT(0);
    }

    y  -= A  * x;
    y_ -= A_ * x;

    if (! lapack::isIdentical(y, y_, "y", "y_")) {
        cerr << endl << "failed: y -= A*x" << endl;
        ASSERT(0);
    }

    ElementType  alpha, beta;
    for (int test=1; test<=20; ++test) {

        alpha = std::pow(2, 5 - std::max(1, rand() % 10));
        beta  = std::pow(2, 5 - std::max(1, rand() % 10));

        //
        // Reset y (and y_) to some random vector
        //
        for (int i=1; i<=n; ++i) {
            y(i) = rand() % 1000;
        }
        y_ = y;

        y  = beta*y + alpha*A  * x;
        y_ = beta*y_ + alpha*A_ * x;

        if (! lapack::isIdentical(y, y_, "y", "y_")) {
            cerr << endl << "failed: y = beta*y + alpha*A*x" << endl;
            cout << "alpha = " << alpha << endl;
            cout << "beta  = " << beta << endl;
            cout << "A_  = " << A_ << endl;
            cout << "x  = " << x << endl;
            cout << "y_  = " << y_ << endl;
            ASSERT(0);
        }

    }
}

int
main()
{
    srand(SEED);

    for (int run=1; run<=30; ++run) {
        const int n       = std::max(1, rand() % (MAX_N));
        // check case 'nnz==0' at least onece
        const int max_nnz = (run==1) ? 0 : (rand() % (MAX_NNZ));

        cerr << "run " << run << ":" << endl;

        for (int indexBase=-3; indexBase<=3; ++indexBase) {
            cerr << "indexBase = " << indexBase << endl;
            cerr << "n x n = " << n << " x " << n << endl;
            cerr << "max_nnz =   " << max_nnz << endl << endl;

            //
            //  Storing elements in upper part
            //
            {
                HeCCSMatrix<CCS<complex<double> > >       A;
                HeMatrix<FullStorage<complex<double> > >  A_;

                setup(Upper, n, max_nnz, indexBase, A, A_);

                mv(n, max_nnz, A, A_);
            }

            //
            //  Storing elements in lower part
            //
            {
                HeCCSMatrix<CCS<complex<double> > >       A;
                HeMatrix<FullStorage<complex<double> > >  A_;

                setup(Lower, n, max_nnz, indexBase, A, A_);

                mv(n, max_nnz, A, A_);
            }
         }
    }
}
