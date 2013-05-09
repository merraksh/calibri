/*
 * optimal calibrated sampling -- bound reduction based on cutoff
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#include "calBT.hpp"
#include "calInstance.hpp"
#include "calModel.hpp"

//#define DEBUG

//
// Bound tightening: use cutoff to set bounds on all delta_i (= w_i -
// d_i) variables
//

void calBT::generateCuts (const OsiSolverInterface & si, 
			  OsiCuts & cs,
			  const CglTreeInfo info) const {

  int
    N = instance_ -> N (),
    *indicesL = new int [N],
    *indicesU = new int [N],
    nLower = 0, 
    nUpper = 0;

  double 
    cutoff = model_ -> bestObj (),
    *newlb = new double [N],
    *newub = new double [N];

  const double
    *lb = si. getColLower () + 1, // to start from delta [0]
    *ub = si. getColUpper () + 1;

  for (int i=1; i<=N; ++i) {

    // TODO: if objective is sum delta/d, newlb/ub has to be multiplied by d_i

    if (*lb++ < -cutoff) {indicesL [nLower] = i; newlb [nLower++] = -cutoff;}
    if (*ub++ >  cutoff) {indicesU [nUpper] = i; newub [nUpper++] =  cutoff;}
  }

  if (nLower || nUpper) {

    OsiColCut *cut = new OsiColCut;

    if (nLower) cut -> setLbs (nLower, indicesL, newlb);
    if (nUpper) cut -> setUbs (nUpper, indicesU, newub);

    cut -> setGloballyValid (true);

#ifdef DEBUG
    cut -> print ();
#endif

    cs.insert (cut);
  
    delete cut;
  }

  delete [] newlb;
  delete [] newub;
  delete [] indicesL;
  delete [] indicesU;
}
