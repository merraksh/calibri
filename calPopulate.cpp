/*
 * optimal calibrated sampling -- Branch and bound -- populate problem
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <OsiSolverInterface.hpp>
#include <CoinFinite.hpp>
#include <CoinHelperFunctions.hpp>

#include "calInstance.hpp"

//#define DEBUG

#define EPS_W 1e-2

/*
 * Building model
 */

int populate (calInstance *instance, OsiSolverInterface *problem) {

  /*
   * construct model:
   *   - n variables
   *   - objective function
   *   - cardinality constraint
   *   - n+1 essential conic constraints
   */

  int
    N = instance -> N (),
    n = instance -> n (),
    p = instance -> p (),

    nPoints = 2*N, // only add them along coordinate axes -- 1+N if approx cone
    nvars = 1 + 2*N,
    ncons = nPoints + 2 + p + 2*N,
    nnz   = 2*p*N + 7*N + nPoints * (1+N),

    *mcnt  = (int    *) malloc (nvars     * sizeof (int)),
    *mbeg  = (int    *) malloc ((1+nvars) * sizeof (int)),
    *mind  = (int    *) malloc (nnz       * sizeof (int)),
    *isInt = (int    *) malloc (N         * sizeof (int));

  double
    *rlb  = (double *) malloc (ncons * sizeof (double)),
    *rub  = (double *) malloc (ncons * sizeof (double)),
    *obj  = (double *) malloc (nvars * sizeof (double)),
    *lb   = (double *) malloc (nvars * sizeof (double)),
    *ub   = (double *) malloc (nvars * sizeof (double)),
    *mval = (double *) malloc (nnz   * sizeof (double));

  //  char *ctype = (char *) malloc (nvars * sizeof (char));

  //  const double *cost_val = instance -> c () . getElements ();
  //  const int    *cost_ind = instance -> c () . getIndices  ();

  double **point = (double **) malloc (nPoints * sizeof (double *));

  for (int i=0; i<nPoints; ++i)
    point [i] = (double *) malloc (N * sizeof (double));

  // nPoints = N+1

  // point [0] [0] = -1.; // first set is the set of two extremes of the interval [-1,1]
  // point [1] [0] =  1.;

  // for (int i=1; i<N; ++i) {

  //// compute r = x1' x2
  //// set alpha = r/(r-1)
  //// set beta  = sqrt (1-alpha^2)
  //// multiply everything [0:i+1,0:i] by beta
  //// set [0:i+1,i] as -alpha
  //// set [i+1,0:i] as 0
  //// set [i+1,i]   as 1

  //double dotproduct = 0, alpha, beta;
  //   for (int j=0; j<i; ++j)
  //  dotproduct += point [0] [j] * point [1] [j];
  //alpha = dotproduct / (dotproduct - 1);
  //beta  = sqrt (1 - alpha*alpha);
  //for (int j=0; j<=i; ++j)
  //     for (int k=0; k<i; ++k)
  //	point [j][k] *= beta;
  //for (int j=0; j<=i; ++j) point [j]   [i] = -alpha;
  //for (int j=0; j<i;  ++j) point [i+1] [j] = 0;
  //point                          [i+1] [i] = 1;
  // }	

  for (int i=0; i<N; ++i) {
    CoinZeroN (point [i], N);
    point [i] [i] = 1;
    CoinZeroN (point [N+i], N);
    point [N+i] [i] = -1;
  }

  nnz = 0;

  // z variable (easy) -------------------------------------------------

  mcnt[0] = nPoints;
  mbeg[0] = 0;

  for (int i=0; i<nPoints; ++i) {
    mind [i] = i;
    double normSq = 0.;
    // coeff of z at i-th constraint is the norm of point[i]
    for (int j=0; j<N; ++j) 
      normSq += point [i][j] * point [i][j];
    mval [i] = - sqrt (normSq);
  }

  nnz = nPoints;

  CoinPackedVector **X = instance -> X ();
  //  double            *d = instance -> d ();

  int
    **indices  = (int    **) malloc (p * sizeof (int    *)),
    *cursor    = (int     *) malloc (p * sizeof (int     )),
    *nElem     = (int     *) malloc (p * sizeof (int     ));

  double 
    **elements = (double **) malloc (p * sizeof (double *));

  for (int i=0; i<p; ++i) {
    cursor   [i] = 0; // position at beginning of both arrays
    nElem    [i] = X [i] -> getNumElements ();
    indices  [i] = X [i] -> getIndices  ();
    elements [i] = X [i] -> getElements ();
  }

  // delta variables (not so easy) -------------------------------------

  for (int i=0; i<N; ++i) {

    mbeg [i+1] = nnz;

    //U -= d [i];

    //	int nterms = 0;
    // for the i-th delta variable,
    // 1) add the i-th element of each point[] vector
    // 2) add a n/N in the (N+3)-rd constraint, indexed N+2

    for (int j=0; j<nPoints; ++j) {
      mind [nnz]   = j;
      mval [nnz++] = point [j] [i];
    }

    mind [nnz]   = 1+nPoints;
    mval [nnz++] = (double) 1; // sum of weights constraint

    for (int j=0; j<p; ++j)
      if ((cursor [j] < nElem [j]))
	if ((indices [j] [cursor [j]] == i)) {
	  mind [nnz]   = nPoints+2+j;
	  mval [nnz++] = elements [j] [cursor [j]++];
	}

    mind [nnz]   = nPoints+2+p+i;
    mval [nnz++] = 1;

    mind [nnz]   = nPoints+2+p+N+i;
    mval [nnz++] = 1;

    mcnt [1+i] = nnz - mbeg [1+i];
  }

  double
    U  = (double) (N*N) / n,
    w0 = (double)  N    / n;

  for (int i=0; i<p; ++i)
    cursor [i] = 0; // reposition at beginning of both arrays

  // s variables --------------------------------------------------------

  for (int i=0; i<N; ++i) {

    mbeg [1+N+i] = nnz;

    //	int nterms = 0;

    // for the i-th s variable,
    //
    // 0) add a one in the second class of constraints, the single constraint ||s||_1 = n
    // 1) add the i-th element of each of the p point[] vector multiplied by d_i
    // 2) add a one in the (N+3)-rd constraint, indexed N+2

    mind [nnz]   = nPoints;
    mval [nnz++] = 1;

    //mind [nnz]   = nPoints+1;
    //mval [nnz++] = (double) N/n; //d [i];

    for (int j=0; j<p; ++j)
      if ((cursor [j] < nElem [j]) && (indices [j] [cursor [j]] == i)) {
        mind [nnz]   = nPoints+2+j;
	mval [nnz++] = (double) N / n * elements [j] [cursor [j]++];
      }

    mind [nnz]   = nPoints+2+p+i;
    mval [nnz++] = w0;

    mind [nnz]   = nPoints+2+p+N+i;
    mval [nnz++] = w0 - U;

    mcnt [1+N+i] = nnz - mbeg [1+i];
  }

  // free points ///////////////////////////////////////////

  for (int i=0; i<nPoints; ++i)
    free (point [i]);

  free (point);

  // set lower and upper bounds on variables ///////////////

  mbeg [1+2*N] = nnz;

  lb [0] = 0;
  ub [0] = COIN_DBL_MAX;

  for (int i=0; i<N; ++i) {
    lb [1+i] =  -w0 + EPS_W; // bounds on deltas
    ub [1+i] = U-w0;         // U - d [i];

    lb [1+N+i] = 0; // s variables are binary
    ub [1+N+i] = 1;
  }

  // set rhs for constraints other than first N+1 (which approximate cone of obj. function)

  for (int i=0; i<nPoints; ++i) {

    rlb [i] = -COIN_DBL_MAX; // while we are at it, set rhs.
    rub [i] = 0;
  }

  rlb [nPoints]   = rub [nPoints]   = n;
  rlb [nPoints+1] = rub [nPoints+1] = 0;

  for (int i=0; i<p; ++i) {

    double sumX = 0,
      *xelems = X [i] -> getElements    ();
    int xsize = X [i] -> getNumElements ();

    for (int j=0; j<xsize; ++j)
      sumX += xelems [j];

    //printf ("sumX = %g\n", sumX);
    rlb [nPoints+2+i] = rub [nPoints+2+i] = sumX;
  }

  for (int i=0; i<N; ++i) {
    rlb [nPoints+2+p+i]   = 0;             rub [nPoints+2+p+i]   = COIN_DBL_MAX;
    rlb [nPoints+2+p+N+i] = -COIN_DBL_MAX; rub [nPoints+2+p+N+i] = 0;
  }

  for (int i=0; i<N; ++i) 
    isInt [i] = 1+N+i;

  obj [0] = 1;

  CoinZeroN (obj + 1, 2*N);

  problem -> loadProblem (nvars, ncons, mbeg, mind, mval, lb, ub, obj, rlb, rub);
  problem -> setInteger (isInt, N);

#ifdef DEBUG
  problem -> writeLp ("test");
#endif

  free (rlb);
  free (rub);

  free (obj);
  free (lb);
  free (ub);

  free (mbeg);
  free (mcnt);
  free (mval);
  free (mind);
  free (isInt);

  free (indices);
  free (cursor);
  free (nElem);
  free (elements);

  //  free (ctype);

  return 0;
}
