/*
 * optimal calibrated sampling -- Cube heuristic (Tillé and De Ville),
 * projection used in flight phase
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"

#include "calInstance.hpp"
#include "calCube.hpp"

#define F77_FUNC(lcase, UCASE) lcase ## _

//#define DEBUG

extern "C" {

  /* Lapack routine to compute an LQ decomposition (in Fortran) */

  void F77_FUNC(dgelqf,DGELQF) (int    *m,      // 
				int    *n,      //
				double *a,      //
				int    *lda,    //
				double *tau,    //
				double *work,   //
				int    *lwork,  //
				int    *info);  //

  /* Lapack routine to compute matrix Q from LQ decomposition computed by dgelqf */

  void F77_FUNC(dorglq,dorglq) (int    *m,
				int    *n,
				int    *k,
				double *a, 
				int    *lda, 
				double *tau, 
				double *work, 
				int    *lwork, 
				int    *info);

  /* Lapack routine to compute the SVD of a matrix A (in Fortran) */

  void F77_FUNC(dgesvd,DGESVD) (char   *jobu,   //
				char   *jobvt, 	//
				int    *m, 	//
				int    *n, 	//
				double *a, 	//
				int    *lda, 	//
				double *s, 	//
				double *u, 	//
				int    *ldu,    //
				double *vt, 	//
				int    *ldvt,	//
				double *work, 	//
				int    *lwork, 	//
				int    *info);	//
}

// project v on null space of restricted calibration constraints
void CalCubeHeur::project (double *v, double *u, double *pi) {

  // Purpose: project vector W*v onto null(A*W), where W = W' is a
  // diagonal matrix with w_ii=0 if pi_i in {0,1} and w_ii=1
  // otherwise.
  //
  // Generating W*v and A*W is easy: just set to zero all columns A_i
  // of A and all elements v_i of v for each i such that pi_i in {0,1}
  //
  // Projecting W*V onto null(A*W) yields a vector u such that 
  //
  // u = W * v - W*A' * (A*W*A')^-1 * A * W*v
  //
  // where A_g W*A' * (A*W*A')^-1 is a generalized inverse of
  // A*W*A'. To obtain it, note that A is an m-by-n matrix with m <<
  // n, therefore we can work on the LQ decomposition of A.
  //
  // 1) Eliminate columns of A and elements of v according to pi
  // 2) (L,Q) = ( [Lp,0], [Q_1;Q_2] ) <- LQ decomposition of A
  //
  // Note that Q = [Q_11, Q_12; Q_21, Q_22] where the lower block
  // [Q_21, Q_22] is unnecessary. Check if upper block can be obtained
  // through appropriate procedure
  //
  // 3) (Up,Sp,Vp) <- SVD (Lp)
  //
  // 4) u = v - |Q_11' * Vp * inv (Sp) * Up' * A| * v
  //            |Q_12' * Vp * inv (Sp) * Up' * A|
  //
  // Et voilá...

  // Note: selective filling of A is in for loop for LQ decomposition

  // LQ decomposition of A //////////////////////////////////////////////////

  int
    N  = instance_ -> N (),
    p  = instance_ -> p (),
    pp = p + 1, // real number of lines
    lwork = N, 
    info;

  double
    *A    = new double [pp * N],
    *tau  = new double      [N],
    *work = new double      [N];

  CoinZeroN (A, pp*N);

#ifdef DEBUG
  printVec (v,  N, "v");
  printVec (pi, N, "pi");
#endif

  for (int i=0; i<pp; ++i) {

    CoinPackedVector *x = instance_ -> X () [i];

    double *elements = x -> getElements    ();
    int    *indices  = x -> getIndices     (),
      numEl          = x -> getNumElements ();

#ifdef DEBUG
    printVec (elements, x -> getNumElements (), "elem");
#endif

    for (int j=0; j<numEl; ++j) {

      int ind = indices [j];

      if (fabs (pi [ind] - .5) < (.5 - 1e-5)) // this means that pi [ind] is fractional
	A [ind*pp + i] = elements [j];
      else v [ind] = 0;
    }
  }

#ifdef DEBUG
  printVec (v,  N, "v spars'd");
  printMatr (A, pp, N, "A");
#endif

  F77_FUNC           // <------------------------- LQ call
    (dgelqf,DGELQF)
    (&pp, &N, A, &pp, tau, work, &lwork, &info);

  //printMatr (A, pp, N, "LQ");

  double *L = new double [pp * pp];

  CoinZeroN (L, pp * pp);

  // copy lower triangle and diagonal of A into L for later applying SVD

  for   (int i=0; i<pp; ++i)
    for (int j=0; j<=i; ++j) {

      double &ael = A [j * pp + i];
      L [j * pp + i] = ael;
      ael = (i==j) ? 1 : 0; // reset it for use in computation with Q
    }

  //printMatr (L, pp, pp, "L");
  //printMatr (A, pp, N,  "A");

  F77_FUNC
    (dorglq,DORGLQ)
    (&pp, &N, &pp, A, &pp, tau, work, &lwork, &info);

  delete [] tau; // passed between LQ decomp (dgekqf) and
		 // reflectors_to_Q (dorglq), we don't care about it

  double *Q = A; // now A contains the important rows of Q, i.e., the
		 // first pp rows

  // let's check if LQ = A

#ifdef DEBUG
  {
    double *AA = new double [pp * N];

    CoinZeroN (AA, pp*N);

    for     (int i=0; i<pp; ++i)
      for   (int j=0; j<N;  ++j)
	for (int k=0; k<pp; ++k)
	  AA [j*pp + i] += L [k*pp+i] * Q [j*pp+k];

    //printMatr (AA, pp, N, "(A =) LQ?");

    delete [] AA;
  }
#endif

  //printMatr (Q, pp, N, "Q");

  // SVD of L ///////////////////////////////////////////////////////////////

  char jobu = 'A', jobvt = 'A';

  double
    *S  = new double [pp],
    *VT = new double [pp*pp],
    *U  = new double [pp*pp];

  lwork = -1;

  double *Lcopy = CoinCopyOfArray (L, pp * pp);

  F77_FUNC           // <------------------------- SVD dry call: just get the right lwork
    (dgesvd,DGESVD)  
    (&jobu, &jobvt, &pp, &pp, L, &pp, S, U, &pp, VT, &pp, work, &lwork, &info);

  lwork = (int)(*work); // work [0] contains suggested space to reserve for the real run
  delete [] work;
  work = new double [lwork];

  F77_FUNC           // <------------------------- SVD call
    (dgesvd,DGESVD) 
    (&jobu, &jobvt, &pp, &pp, L, &pp, S, U, &pp, VT, &pp, work, &lwork, &info);

  delete [] work;

  //printMatr (U,  pp, pp, "U");
  //printMatr (VT, pp, pp, "Vt");
  //printVec  (S,  pp,     "S");

#ifdef DEBUG
  {
    // compute L = U S Vt

    double *USVt = new double [pp * pp];

    //CoinZeroN (USVt, pp*pp);

    for     (int i=0; i<pp; ++i)
      for   (int j=0; j<pp; ++j) {

	USVt [j*pp + i] = - Lcopy [j*pp + i];
	for (int k=0; k<pp; ++k)
	  USVt [j*pp + i] += U [k*pp+i] * S [k] * VT [j*pp+k];
      }

    //printMatr (USVt, pp, pp, "L - U S Vt");

    delete [] USVt;
  }
#endif

  delete [] Lcopy;

  // build u ////////////////////////////////////////////////////////////////
  //
  // u = v - |Q_11' * Vp * inv (Sp) * Up' * A| * v
  //         |Q_12' * Vp * inv (Sp) * Up' * A|
  //
  //   = v - |Q_11' * Vp * inv (Sp) * Up' * A * v|
  //         |Q_12' * Vp * inv (Sp) * Up' * A * v|
  //
  // Compute it backwards since v is a vector. Complexity is O(pN+p^2) = O(pN) since p<<N
  //

  double *Av = new double [pp];

  CoinZeroN (Av, pp);

  // <------------------------ Compute A*v

  for (int i=0; i<pp; ++i) {

    CoinPackedVector *x = instance_ -> X () [i];

    double *elements = x -> getElements    ();
    int    *indices  = x -> getIndices     (),
      numEl          = x -> getNumElements ();

    for (int j=0; j<numEl; ++j) {

      int ind = indices [j];

      if (fabs (pi [ind] - .5) < (.5 - 1e-5)) // this means that pi [ind] is fractional
	Av [i] += elements [j] * v [ind];
    }
  }

  //printVec (Av, pp, "A*v");

  // <------------------------ Multiply Up' by A * v

  double *UtAv = new double [pp];

  CoinZeroN (UtAv, pp);

  for   (int i=0; i<pp; ++i)
    for (int j=0; j<pp; ++j)
      UtAv [i] += U [i*pp + j] * Av [j];

  //printVec (UtAv, pp, "U'*A*v");

  delete [] Av;

  // <------------------------ Pre-multiply Sigma^-1 by obtained vector

  for (int i=0; i<pp; ++i)
    UtAv [i] /= S [i];

  delete [] S;

  // <------------------------ Pre-multiply VT' by obtained vector

  double *VUtAv = new double [pp];

  CoinZeroN (VUtAv, pp);

  for   (int i=0; i<pp; ++i)
    for (int j=0; j<pp; ++j)
      VUtAv [i] += VT [i*pp + j] * UtAv [j];  

  delete [] UtAv;

  // <------------------------ Pre-multiply (VT'*Sinv*)UtAv by [Q_11,Q_12]' to obtain u

  //  CoinZeroN (u, N);

  for   (int i=0; i<N; ++i) {

    double &ucurr = u [i] = v [i];

    for (int j=0; j<pp; ++j)
      ucurr -= Q [i*pp + j] * VUtAv [j];
  }

#ifdef DEBUG
  {
    // now check if Au=0

    double *Au = new double [pp];

    CoinZeroN (Au, pp);

    for   (int i=0; i<pp; ++i)
      for (int j=0; j<N; ++j)
	Au [i] += A [j*pp + i] * u [j];

    //printVec (Au, pp, "A*u");
    //printVec (u,  N,  "u");

    delete [] Au;
  }
#endif

  delete [] VUtAv;

  delete [] A;
  delete [] VT;
  delete [] U;
  delete [] L;
}
