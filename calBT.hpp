/*
 * optimal calibrated sampling -- bound tightening
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#ifndef calBT_hpp
#define calBT_hpp

#include <CglCutGenerator.hpp>

class calInstance;
class calModel;
class OsiSolverInterface;
class OsiCuts;

//
// Cutoff inequalities: the objective function is ||delta||_2, so if z
// is a valid cutoff a valid bound tightening is |delta_i| <= cutoff
//

class calBT: public CglCutGenerator {

protected:

  calInstance *instance_;
  calModel    *model_;

public:

  calBT (calInstance *inst, calModel *model): instance_ (inst), model_ (model) {}

  calBT *clone () const {return new calBT (instance_, model_);}

  void generateCuts (const OsiSolverInterface & si, 
		     OsiCuts & cs,
		     const CglTreeInfo info = CglTreeInfo ()) const; 
};

#endif
