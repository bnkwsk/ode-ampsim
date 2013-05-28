/*
 *   Copyright (c) 2012, Michael Lehn
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
       SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV,
      $                   EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
 *
 *  -- LAPACK routine (version 3.3.1)                                  --
 *
 *  -- Contributed by Zlatko Drmac of the University of Zagreb and     --
 *  -- Kresimir Veselic of the Fernuniversitaet Hagen                  --
 *  -- April 2011                                                      --
 *
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 * This routine is also part of SIGMA (version 1.23, October 23. 2008.)
 * SIGMA is a library of algorithms for highly accurate algorithms for
 * computation of SVD, PSVD, QSVD, (H,K)-SVD, and for solution of the
 * eigenvalue problems Hx = lambda M x, H M x = lambda x with H, M > 0.
 *
 */

#ifndef FLENS_LAPACK_GE_SVJ1_TCC
#define FLENS_LAPACK_GE_SVJ1_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================
namespace generic {

template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj1_impl(SVJ::JobV                                  jobV,
          typename GeMatrix<MA>::IndexType           n1,
          GeMatrix<MA>                               &A,
          DenseVector<VD>                            &d,
          DenseVector<VSVA>                          &sva,
          GeMatrix<MV>                               &V,
          const typename GeMatrix<MA>::ElementType   &eps,
          const typename GeMatrix<MA>::ElementType   &safeMin,
          const typename GeMatrix<MA>::ElementType   &tol,
          typename GeMatrix<MA>::IndexType           nSweep,
          DenseVector<VWORK>                         &work)
{
    using std::abs;
    using std::max;
    using std::min;
    using flens::pow;
    using std::sqrt;
    using std::swap;

    typedef typename GeMatrix<MA>::ElementType  ElementType;
    typedef typename GeMatrix<MA>::IndexType    IndexType;

    const ElementType  Zero(0), Half(0.5), One(1);

    const Underscore<IndexType>  _;

    ElementType fastr_data[5];
    DenseVectorView<ElementType>
        fastr  = typename DenseVectorView<ElementType>::Engine(5, fastr_data);

    const IndexType  m     = A.numRows();
    const IndexType  n     = A.numCols();

    auto _work = work(_(1,m));

    const bool applyV   = (jobV==SVJ::ApplyV);
    const bool rhsVec   = (jobV==SVJ::ComputeV) || applyV;

    const ElementType rootEps = sqrt(eps);
    const ElementType rootSafeMin = sqrt(safeMin);
    const ElementType small = safeMin / eps;
    const ElementType big = One / safeMin;
    const ElementType rootBig = One / rootSafeMin;
    const ElementType bigTheta = One / rootEps;
    const ElementType rootTol = sqrt(tol);

    IndexType  info = 0;
//
//  .. Initialize the right singular vector matrix ..
//
//  RSVEC = LSAME( JOBV, 'Y' )
//
    const IndexType emptsw = n1*(n-n1);
    fastr(1) = Zero;
//
//  .. Row-cyclic pivot strategy with de Rijk's pivoting ..
//
    const IndexType kbl = min(IndexType(8),n);
    IndexType nblr = n1 / kbl;
    if (nblr*kbl!=n1) {
        ++nblr;
    }
//  .. the tiling is nblr-by-nblc [tiles]
    IndexType nblc = (n-n1) / kbl;

    if (nblc*kbl!=(n-n1)) {
        ++nblc;
    }
    const IndexType blSkip = pow(kbl,2) + 1;
//[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
    const IndexType rowSkip = min(IndexType(5),kbl);
//[TP] ROWSKIP is a tuning parameter.
    IndexType swBand = 0;
//[TP] SWBAND is a tuning parameter. It is meaningful and effective
//  if SGESVJ is used as a computational routine in the preconditioned
//  Jacobi SVD algorithm SGESVJ.
//
//
//  | *   *   * [x] [x] [x]|
//  | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
//  | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
//  |[x] [x] [x] *   *   * |
//  |[x] [x] [x] *   *   * |
//  |[x] [x] [x] *   *   * |
//
//
    ElementType aapp, aapp0, aaqq, aapq, aqoap, apoaq;
    ElementType tmp, cs, sn, t, theta, thetaSign;

    bool rotOk;
    bool converged = false;

    for (IndexType i=1; i<=nSweep; ++i) {
//  .. go go go ...
//
       ElementType max_aapq = Zero;
       ElementType max_sinj = Zero;

       IndexType iswRot = 0;
       IndexType notRot = 0;
       IndexType pSkipped = 0;

       for (IndexType ibr=1; ibr<=nblr; ++ibr) {

          IndexType igl = (ibr-1)*kbl + 1;
//
//
//........................................................
//... go to the off diagonal blocks

          for (IndexType jbc=1; jbc<=nblc; ++jbc) {

             IndexType jgl = n1 + (jbc-1)*kbl + 1;

//     doing the block at ( ibr, jbc )

             IndexType ijblsk = 0;
             for (IndexType p=igl; p<=min(igl+kbl-1,n1); ++p) {

                aapp = sva(p);

                if (aapp>Zero) {

                   pSkipped = 0;

                   for (IndexType q=jgl; q<=min(jgl+kbl-1,n); ++q) {

                      aaqq = sva(q);

                      if (aaqq>Zero) {
                         aapp0 = aapp;
//
//  .. M x 2 Jacobi SVD ..
//
//     .. Safe Gram matrix computation ..
//
                         if (aaqq>=One) {
                            if (aapp>=aaqq) {
                               rotOk = small*aapp<=aaqq;
                            } else {
                               rotOk = small*aaqq<=aapp;
                            }
                            if (aapp<big/aaqq) {
                               aapq = (A(_,p)*A(_,q)*d(p)*d(q)/aaqq)/aapp;
                            } else {
                               _work = A(_,p);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aapp, d(p), _work);
                               aapq = _work*A(_,q)*d(q)/aaqq;
                            }
                         } else {
                            if (aapp>=aaqq) {
                               rotOk = aapp<=aaqq/small;
                            } else {
                               rotOk = aaqq<=aapp/small;
                            }
                            if (aapp>small/aaqq) {
                               aapq = (A(_,p)*A(_,q)*d(p)*d(q)/aaqq)/aapp;
                            } else {
                               _work = A(_,q);
                               lascl(LASCL::FullMatrix, 0, 0,
                                     aaqq, d(q), _work);
                               aapq = _work*A(_,p)*d(p)/aapp;
                            }
                         }

                         max_aapq = max(max_aapq, abs(aapq));

//     TO rotate or NOT to rotate, THAT is the question ...
//
                         if (abs(aapq)>tol) {
                            notRot = 0;
//        ROTATED  = ROTATED + 1
                            pSkipped = 0;
                            ++iswRot;

                            if (rotOk) {

                               aqoap = aaqq / aapp;
                               apoaq = aapp / aaqq;
                               theta = -Half*abs(aqoap-apoaq) / aapq;
                               if (aaqq>aapp0) {
                                  theta = -theta;
                               }

                               if (abs(theta)>bigTheta) {
                                  t = Half / theta;
                                  fastr(3) =  t*d(p) / d(q);
                                  fastr(4) = -t*d(q) / d(p);
                                  blas::rotm(A(_,p), A(_,q), fastr);
                                  if (rhsVec) {
                                     blas::rotm(V(_,p), V(_,q), fastr);
                                  }
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));
                                  max_sinj = max(max_sinj, abs(t));
                               } else {
//
//              .. choose correct signum for THETA and rotate
//
                                  thetaSign = -sign(One, aapq);
                                  if (aaqq>aapp0) {
                                     thetaSign = -thetaSign;
                                  }
                                  t = One
                                     / (theta+thetaSign*sqrt(One+theta*theta));
                                  cs = sqrt(One / (One+t*t));
                                  sn = t*cs;
                                  max_sinj = max(max_sinj, abs(sn));
                                  sva(q) = aaqq
                                          *sqrt(max(Zero, One+t*apoaq*aapq));
                                  aapp *= sqrt(max(Zero, One-t*aqoap*aapq));

                                  apoaq = d(p) / d(q);
                                  aqoap = d(q) / d(p);
                                  if (d(p)>=One) {

                                     if (d(q)>=One) {
                                        fastr(3) =  t*apoaq;
                                        fastr(4) = -t*aqoap;
                                        d(p) *= cs;
                                        d(q) *= cs;
                                        blas::rotm(A(_,p), A(_,q), fastr);
                                        if (rhsVec) {
                                           blas::rotm(V(_,p), V(_,q), fastr);
                                        }
                                     } else {
                                        A(_,p) -= t*aqoap*A(_,q);
                                        A(_,q) += cs*sn*apoaq*A(_,p);
                                        if (rhsVec) {
                                           V(_,p) -= t*aqoap*V(_,q);
                                           V(_,q) += cs*sn*apoaq*V(_,p);
                                        }
                                        d(p) *= cs;
                                        d(q) /= cs;
                                     }
                                  } else {
                                     if (d(q)>=One) {
                                        A(_,q) += t*apoaq*A(_,p);
                                        A(_,p) -= cs*sn*aqoap*A(_,q);
                                        if (rhsVec) {
                                           V(_,q) += t*apoaq*V(_,p);
                                           V(_,p) -= cs*sn*aqoap*V(_,q);
                                        }
                                        d(p) /= cs;
                                        d(q) *= cs;
                                     } else {
                                        if (d(p)>=d(q)) {
                                           A(_,p) -= t*aqoap*A(_,q);
                                           A(_,q) += cs*sn*apoaq*A(_,p);
                                           d(p) *= cs;
                                           d(q) /= cs;
                                           if (rhsVec) {
                                              V(_,p) -= t*aqoap*V(_,q);
                                              V(_,q) += cs*sn*apoaq*V(_,p);
                                           }
                                        } else {
                                           A(_,q) += t*apoaq*A(_,p);
                                           A(_,p) -= cs*sn*aqoap*A(_,q);
                                           d(p) /= cs;
                                           d(q) *= cs;
                                           if (rhsVec) {
                                              V(_,q) += t*apoaq*V(_,p);
                                              V(_,p) -= cs*sn*aqoap*V(_,q);
                                           }
                                        }
                                     }
                                  }
                               }

                            } else {
                               if (aapp>aaqq) {
                                  _work = A(_,p);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aapp, One, _work);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aaqq, One, A(_,q));
                                  tmp = -aapq*d(p)/d(q);
                                  A(_,q) += tmp*_work;
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        One, aaqq, A(_,q));
                                  sva(q) = aaqq*sqrt(max(Zero, One-aapq*aapq));
                                  max_sinj = max(max_sinj, safeMin);
                               } else {
                                  _work = A(_,q);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aaqq, One, _work);
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        aapp, One, A(_,p));
                                  tmp = -aapq*d(q)/d(p);
                                  A(_,p) += tmp*_work;
                                  lascl(LASCL::FullMatrix, 0, 0,
                                        One, aapp, A(_,p));
                                  sva(p) = aapp*sqrt(max(Zero, One-aapq*aapq));
                                  max_sinj = max(max_sinj, safeMin);
                               }
                            }
//        END IF ROTOK THEN ... ELSE
//
//        In the case of cancellation in updating SVA(q)
//        .. recompute SVA(q)
                            if (pow(sva(q)/aaqq,2)<=rootEps) {
                               if (aaqq<rootBig && aaqq>rootSafeMin) {
                                  sva(q) = blas::nrm2(A(_,q))*d(q);
                               } else {
                                  t = Zero;
                                  aaqq = One;
                                  lassq(A(_,q), t, aaqq);
                                  sva(q) = t*sqrt(aaqq)*d(q);
                               }
                            }
                            if (pow(aapp/aapp0,2)<=rootEps) {
                               if (aapp<rootBig && aapp>rootSafeMin) {
                                  aapp = blas::nrm2(A(_,p))*d(p);
                               } else {
                                  t = Zero;
                                  aapp = One;
                                  lassq(A(_,p), t, aapp);
                                  aapp = t*sqrt(aapp)*d(p);
                               }
                               sva(p) = aapp;
                            }
//           end of OK rotation
                         } else {
                            ++notRot;
//        SKIPPED  = SKIPPED  + 1
                            ++pSkipped;
                            ++ijblsk;
                         }
                      } else {
                         ++notRot;
                         ++pSkipped;
                         ++ijblsk;
                      }

//   IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                      if (i<=swBand && ijblsk>=blSkip) {
                         sva(p) = aapp;
                         notRot = 0;
                         goto jbcLooExit;
                      }
                      if (i<=swBand && pSkipped>rowSkip) {
                         aapp = -aapp;
                         notRot = 0;
                         break;
                      }

                   }
//     end of the q-loop

                   sva(p) = aapp;

                } else {
                   if (aapp==Zero) {
                      notRot += min(jgl+kbl-1,n) - jgl + 1;
                   }
                   if (aapp<Zero) {
                      notRot = 0;
                   }
//**   IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                }

             }
//  end of the p-loop
          }
//  end of the jbc-loop
       jbcLooExit:
//2011 bailed out of the jbc-loop
          for (IndexType p=igl; p<=min(igl+kbl-1,n); ++p) {
             sva(p) = abs(sva(p));
          }
//**   IF ( NOTROT .GE. EMPTSW ) GO TO 1994
       }
//2000 :: end of the ibr-loop
//
//  .. update SVA(N)
       if (sva(n)<rootBig && sva(n)>rootSafeMin) {
          sva(n) = blas::nrm2(A(_,n))*d(n);
       } else {
          t = Zero;
          aapp = One;
          lassq(A(_,n), t, aapp);
          sva(n) = t*sqrt(aapp)*d(n);
       }
//
//  Additional steering devices
//
       if (i<swBand && (max_aapq<=rootTol || iswRot<=n)) {
          swBand = i;
       }

       if (i>swBand+1 && max_aapq<n*tol && n*max_aapq*max_sinj<tol) {
          converged = true;
          break;
       }

       if (notRot>=emptsw) {
          converged = true;
          break;
       }

    }
//  end i=1:NSWEEP loop
    if (converged) {
//#:) Reaching this point means that during the i-th sweep all pivots were
//  below the given threshold, causing early exit.
       info = 0;
    } else {
//#:) Reaching this point means that the procedure has completed the given
//  number of sweeps.
       info = nSweep - 1;
    }
//  Sort the vector D
//
    for (IndexType p=1; p<=n-1; ++p) {
       const IndexType q = blas::iamax(sva(_(p,n))) + p - 1;
       if (p!=q) {
          swap(sva(p), sva(q));
          swap(d(p), d(q));
          blas::swap(A(_,p), A(_,q));
          if (rhsVec) {
             blas::swap(V(_,p), V(_,q));
          }
       }
    }
    return info;
}

} // namespace generic


//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj1_impl(SVJ::JobV                                  jobV,
          typename GeMatrix<MA>::IndexType           n1,
          GeMatrix<MA>                               &A,
          DenseVector<VD>                            &d,
          DenseVector<VSVA>                          &sva,
          GeMatrix<MV>                               &V,
          const typename GeMatrix<MA>::ElementType   &eps,
          const typename GeMatrix<MA>::ElementType   &safeMin,
          const typename GeMatrix<MA>::ElementType   &tol,
          typename GeMatrix<MA>::IndexType           nSweep,
          DenseVector<VWORK>                         &work)
{
    typedef typename GeMatrix<MA>::IndexType  IndexType;

    return cxxlapack::gsvj1<IndexType>(getF77Char(jobV),
                                       A.numRows(),
                                       A.numCols(),
                                       n1,
                                       A.data(),
                                       A.leadingDimension(),
                                       d.data(),
                                       sva.data(),
                                       V.numRows(),
                                       V.data(),
                                       V.leadingDimension(),
                                       eps,
                                       safeMin,
                                       tol,
                                       nSweep,
                                       work.data(),
                                       work.length());
}

} // namespace external

#endif // USE_CXXLAPACK

//== public interface ==========================================================
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename GeMatrix<MA>::IndexType
svj1(SVJ::JobV                                  jobV,
     typename GeMatrix<MA>::IndexType           n1,
     GeMatrix<MA>                               &A,
     DenseVector<VD>                            &d,
     DenseVector<VSVA>                          &sva,
     GeMatrix<MV>                               &V,
     const typename GeMatrix<MA>::ElementType   &eps,
     const typename GeMatrix<MA>::ElementType   &safeMin,
     const typename GeMatrix<MA>::ElementType   &tol,
     typename GeMatrix<MA>::IndexType           nSweep,
     DenseVector<VWORK>                         &work)
{
    using std::max;
    using std::min;

    typedef typename GeMatrix<MA>::IndexType    IndexType;

//
//  Test the input parameters
//
#   ifndef NDEBUG


    ASSERT(A.firstRow()==1);
    ASSERT(A.firstCol()==1);

    const IndexType m = A.numRows();
    const IndexType n = A.numCols();

    ASSERT(m>=n);
    ASSERT(n>=n1);
    ASSERT(n1>=0);


    ASSERT(d.firstIndex()==1);
    ASSERT(d.length()==n);

    ASSERT(sva.firstIndex()==1);
    ASSERT(sva.length()==n);

    ASSERT(V.firstRow()==1);
    ASSERT(V.firstCol()==1);

    if (jobV==SVJ::ComputeV) {
        ASSERT(V.numCols()==n);
        ASSERT(V.numRows()==n);
    }
    if (jobV==SVJ::ApplyV) {
        ASSERT(V.numCols()==n);
    }

    if (work.length()>0) {
        ASSERT(work.length()>=m);
    }
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename GeMatrix<MA>::NoView       A_org     = A;
    typename DenseVector<VD>::NoView    d_org     = d;
    typename DenseVector<VSVA>::NoView  sva_org   = sva;
    typename GeMatrix<MV>::NoView       V_org     = V;
    typename DenseVector<VWORK>::NoView work_org  = work;
#   endif

//
//  Call implementation
//
    IndexType info = LAPACK_SELECT::svj1_impl(jobV, n1, A, d, sva, V,
                                              eps, safeMin, tol,
                                              nSweep, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MA>::NoView       A_generic     = A;
    typename DenseVector<VD>::NoView    d_generic     = d;
    typename DenseVector<VSVA>::NoView  sva_generic   = sva;
    typename GeMatrix<MV>::NoView       V_generic     = V;
    typename DenseVector<VWORK>::NoView work_generic  = work;
//
//  restore output arguments
//
    A    = A_org;
    d    = d_org;
    sva  = sva_org;
    V    = V_org;
    work = work_org;
//
//  Compare generic results with results from the native implementation
//
    IndexType _info = external::svj1_impl(jobV, n1, A, d, sva, V,
                                          eps, safeMin, tol,
                                          nSweep, work);

    bool failed = false;
    if (! isIdentical(A_generic, A, "A_generic", "A")) {
        std::cerr << "CXXLAPACK: A_generic = " << A_generic << std::endl;
        std::cerr << "F77LAPACK: A = " << A << std::endl;
        failed = true;
    }
    if (! isIdentical(d_generic, d, "d_generic", "d")) {
        std::cerr << "CXXLAPACK: d_generic = " << d_generic << std::endl;
        std::cerr << "F77LAPACK: d = " << d << std::endl;
        failed = true;
    }
    if (! isIdentical(sva_generic, sva, "sva_generic", "sva")) {
        std::cerr << "CXXLAPACK: sva_generic = " << sva_generic << std::endl;
        std::cerr << "F77LAPACK: sva = " << sva << std::endl;
        failed = true;
    }
    if (! isIdentical(V_generic, V, "V_generic", "V")) {
        std::cerr << "CXXLAPACK: V_generic = " << V_generic << std::endl;
        std::cerr << "F77LAPACK: V = " << V << std::endl;
        failed = true;
    }
    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }
    if (! isIdentical(info, _info, "info", "_info")) {
        std::cerr << "CXXLAPACK: info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        std::cerr << "error in: svj1.tcc" << std::endl;
        ASSERT(0);
    } else {
        std::cerr << "passed: svj1.tcc" << std::endl;
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MA, typename VD, typename VSVA, typename MV, typename VWORK>
typename MA::IndexType
svj1(SVJ::JobV                                  jobV,
     typename MA::IndexType                     n1,
     MA                                         &&A,
     VD                                         &&d,
     VSVA                                       &&sva,
     MV                                         &&V,
     const typename MA::ElementType             &eps,
     const typename MA::ElementType             &safeMin,
     const typename MA::ElementType             &tol,
     typename MA::IndexType                     nSweep,
     VWORK                                      &&work)
{
    typename MA::IndexType   info;

    CHECKPOINT_ENTER;
    svj1(jobV, n1, A, d, sva, V, eps, safeMin, tol, nSweep, work);
    CHECKPOINT_LEAVE;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_GE_SVJ1_TCC
