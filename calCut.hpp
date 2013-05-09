/*
 * optimal calibrated sampling -- cut generator
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#ifndef calCut_hpp
#define calCut_hpp

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calInstance.hpp"
#include <CglCutGenerator.hpp>

#define MIN_VIOLATION 1e-5
#define maxCallsPerNode 30
#define TOL_LB 1e-3

#define max(a,b) ((a) > (b) ? (a) : (b))

class OsiSolverInterface;
class OsiCuts;

//
// Generate cuts to approximate the second order cone
//

class calCut: public CglCutGenerator {

protected:

  calInstance *instance_;

public:

  calCut (calInstance *inst): instance_ (inst) {}

  calCut *clone () const {return new calCut (instance_);}

  void generateCuts (const OsiSolverInterface & si, 
		     OsiCuts & cs,
		     const CglTreeInfo info = CglTreeInfo ()) const; 
};

#endif
