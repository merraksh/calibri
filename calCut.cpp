/*
 * optimal calibrated sampling -- cut generator
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#include "calCut.hpp"
#include "calCube.hpp"

#define MIN_VIOLATION 1e-5
#define maxCallsPerNode 30
#define TOL_LB 1e-3

#define max(a,b) ((a) > (b) ? (a) : (b))

//#define DEBUG

void calCut::generateCuts (const OsiSolverInterface & si, 
			   OsiCuts & cs,
			   const CglTreeInfo info) const {

  static int ncuts = 0;

  // These cuts aim at approximating a cone function that is
  // expressed, through an inequality of the form z >= || delta ||_2
  //
  // This is accomplished through Outer approximation as follows: any
  // constraint of the form z >= g(x) can be approximated at any point
  // xbar by
  //
  // z >= g(xbar) + grad(g)(xbar) (x - xbar)
  //
  // where grad(g) is the gradient of the multivariate function g. In
  // our case, since g(delta) = || delta ||_2, we have that for any
  // point dbar
  //
  // sum_i (dbar_i * delta_i) - sqrt (sum_i (dbar_i^2)) * z <= sum_i (dbar_i^2)
  //
  // This should only be a means of approximating the value of the
  // objective function, which is not very accurate at the beginning
  // due to the polyhedral conic approximation.

  int 
    N = instance_ -> N ();

  double 
    sumDeltaSq = 0.,
    zCurrent  = si. getColSolution () [0];     // value of z, first variable and objective function

  const double *dCurrent = si. getColSolution () + 1; // pointer to first element of the delta subvector

  zCurrent *= zCurrent;
	
  //  CoinPackedVector xs (si. getNumCols (), si. getColSolution());

  // check if cut violated

  for (int i=0; i<N; ++i) {
    sumDeltaSq += dCurrent [i] * dCurrent [i];

    // can't stop at threshold! z's coeff won't be correct
    //if (sumDeltaSq > zCurrent)
    //  break;
  }

  if (sumDeltaSq <= zCurrent) // no cuts to separate
    return;

  //printf ("sum deltasq = %g, its sqrt = %g\n", sumDeltaSq, sqrt (sumDeltaSq));
  //printVec (dCurrent, N, "dCurrent");

#ifdef DEBUG
  //   static double
  //     prevlb = -1e100;
  //
  //  printf ("CALLED CUT CALLBACK ###########################\n");
  //
  //   if (curnode == nodeNum) {
  //     if ((ncalls++ > maxCallsPerNode) ||
  //   	(fabs ((cur_obj - prevlb) / max (1., fabs (prevlb))) < TOL_LB)) {
  //       // *useraction_p = CPX_CALLBACK_DEFAULT;
  //       return;
  //     } 
  //   } else ncalls = 0;
  //
  //   prevlb = cur_obj;
  //   curnode = nodeNum;
  ////   printf ("cut %d %d %d %g               \r", nodeNum, ncalls, ncuts, cur_obj);
#endif

// #ifdef DEBUG
//   printf ("----------------------------\nAdding cut: %g > %g [%g]",
// 	  sqrt (sQs),
// 	  sqrt (2* beta * K * K),
// 	  sqrt (sQs) - sqrt (2 * beta * K * K));
//   //print_V (xs);
//   //print_V (sQ);
//   printf ("%g >= ", (double) K * sqrt (2 * beta * sQs));
//   //print_V (sQ);
// #endif

  OsiRowCut *cut = new OsiRowCut;

  //CoinPackedVector xs (si. getNumCols (), si. getColSolution());

  int    *indices = new int    [N+1];
  double *coeff   = new double [N+1];

  coeff   [0] = -sqrt (sumDeltaSq); // coefficient of z
  indices [0] = 0;

  CoinCopyN (dCurrent, N, coeff + 1); // coefficient of all deltas

  for (int i=1; i<=N; ++i)
    indices [i] = i;

  cut -> setLb  (-COIN_DBL_MAX);
  cut -> setUb  (0);
  cut -> setRow (1+N, indices, coeff);

  cut -> setGloballyValid (true); // global

#ifdef DEBUG
  cut -> print ();
#endif

  cs.insert (cut);
  
  delete cut;
  delete [] indices;
  delete [] coeff;

  ++ncuts;

  // printf ("->  %d %d %d %g\r", nodeNum, ncalls, ncuts, cur_obj);
}
