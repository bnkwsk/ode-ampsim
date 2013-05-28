#include <complex>
#include <iostream>
#include <flens/flens.cxx>

using namespace flens;
using namespace std;

int
main()
{
    typedef GeMatrix<FullStorage<double> >             RealGeMatrix;
    typedef GeMatrix<FullStorage<complex<double> > >   ComplexGeMatrix;

    RealGeMatrix::IndexVariable i, j;
    RealGeMatrix   A(3,3), B(3,3), C(3,3);
    A(i,j) = Abs(i*i-j)*(i-j);
    B(i,j) = i+j;

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

    // element-wise matrix operations
    C(i,j) = A(i,j) * B(j,i);

    cout << "C = " << C << endl;

    Underscore<int>     _;
    // element-wise vector operations
    C(2,_)(i) = A(1,_)(i) + B(3,_)(i) + 10*i;

    cout << "C = " << C << endl;

    ComplexGeMatrix  V(3,3), W(3,3), Z(3,3);
    V(i,j) = Complex(A(i,j), B(i,j));
    W(i,j) = Complex(B(i,j), A(i,j));

    Z(i,j) = Complex(Cos(M_PI*Real(V(i,j))/4), Sin(-M_PI*Imag(W(i,j))/3));

    cout << "V = " << V << endl;
    cout << "W = " << W << endl;
    cout << "Z = " << Z << endl;

    RealGeMatrix X(3,3);
    X(i,j) = Real(Z(i,j));

    RealGeMatrix Y(3,3);
    Y(i,j) = Imag(Z(i,j));

    cout << "X = " << X << endl;
    cout << "Y = " << Y << endl;

    Z(i,j) = Complex(Real(Z(i,j)), -Imag(Z(i,j)));
    cout << "Z = " << Z << endl;

    A.upper()(i,j) = i + j;
    A.strictLower()(i,j) = 0;
    cout << "A = " << A << endl;
}
