/*
 * optimal calibrated sampling -- inheriting CbcSolver with different checkSolution
 *
 * (C) Pietro Belotti 2013. This code is released under the Eclipse
 * Public License.
 */

#include "calModel.hpp"

//#define DEBUG

/** Call this to really test if a valid solution can be feasible
    Solution is number_columns in size.
    If fixVariables true then bounds of continuous solver updated.
    Returns objective value (worse than cutoff if not feasible)
    Previously computed objective value is now passed in (in case user does not do solve)
*/

double calModel::checkSolution (double cutoff, 
				double *solution,
				int fixVariables, 
				double originalObjValue) {

  // check if z >= || delta ||_2. If not, just return the ideal value of the objective function

  int N = instance_ -> N ();

  double 
    sumDeltaSq = 0.,
    //	zCurrent  = solution [0],     // value of z, first variable and objective function
    *dCurrent = solution + 1; // pointer to first element of the delta subvector

  // check if cut violated

  for (int i=0; i<N; ++i)
    sumDeltaSq += dCurrent [i] * dCurrent [i];

  sumDeltaSq = sqrt (sumDeltaSq);

  if (sumDeltaSq < cutoff) {

    if (!bestSol_)
      bestSol_ = CoinCopyOfArray (solution, 1+2*N);
    else CoinCopyN (solution, 1+2*N, bestSol_);

    bestObj_ = sumDeltaSq;
    bestSol_ [0] = sumDeltaSq;
  }

  //printf ("check sol: %g ---> %g (cutoff is %g)\n", originalObjValue, sqrt (sumDeltaSq), cutoff);
  return sumDeltaSq;
}

// fixes s variables based on s0
void calModel::changeLU (OsiSolverInterface &si, double *s0) {

#ifdef DEBUG
  static int which = 0;
#endif

  int N = instance_ -> N ();

  for (int i=0; i<N; ++i) {

    double s0i = s0 [i];

    if      (s0i <     1e-6) si.setColUpper (1 + N + i, 0.);
    else if (s0i > 1 - 1e-6) si.setColLower (1 + N + i, 1.);
  }

#ifdef DEBUG
  char name [20];
  sprintf (name, "changed_%d", which++);
  si. writeLp (name);
#endif
}
