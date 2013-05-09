/*
 * optimal calibrated sampling -- inheriting CbcSolver with different checkSolution
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#ifndef calModel_hpp
#define calModel_hpp

#include <CbcModel.hpp>
#include "calInstance.hpp"

class CalCubeHeur;

class calModel: public CbcModel {

private:

  calInstance *instance_; ///< input data
  double      *bestSol_;  ///< keep solution from checksolution
  double       bestObj_;  ///< and its obj value

public:

  calModel (const OsiSolverInterface &lp, calInstance *inst):
    CbcModel (lp), instance_ (inst), bestSol_ (NULL), bestObj_ (1e40) {}

  calModel (const calModel &rhs):
    CbcModel (rhs),
    instance_ (rhs.instance_),
    bestSol_  (CoinCopyOfArray (rhs.bestSol_, 2 * rhs.instance_ -> N ())),
    bestObj_  (rhs.bestObj_) {}

  calModel *clone ()
  {return new calModel (*this);}

  ~calModel () {

    if (bestSol_)
      delete [] bestSol_;
  }

  /** Call this to really test if a valid solution can be feasible
      Solution is number columns in size.
      If fixVariables true then bounds of continuous solver updated.
      Returns objective value (worse than cutoff if not feasible)
      Previously computed objective value is now passed in (in case
      user does not do solve)
  */

  double checkSolution (double cutoff, 
  			double *solution,
  			int fixVariables, 
  			double originalObjValue);

  const double *bestSol () {return bestSol_;}
  double bestObj () {return bestObj_;}

  void changeLU (OsiSolverInterface &si, double *s0); // fixes s variables based on s0

  bool search (CalCubeHeur &calCube, FILE *f, int repl);
};

#endif
