/*
 * optimal calibrated sampling -- Cube heuristic (Tillé and De Ville)
 *
 * (C) Pietro Belotti 2013. This code is released under the Eclipse
 * Public License.
 */

#ifndef calQuad_hpp
#define calQuad_hpp

#ifdef WIN32
#define srand48(a) srand((a))
#define drand48() ((double)rand()/RAND_MAX)
#endif

#include "CbcHeuristic.hpp"

class calInstance;

//
// Heuristic to run a variant of the Cube algorithm
//

class CalQuadHeur: public CbcHeuristic {

public:

  // Default Constructor
  CalQuadHeur (calInstance *i);

  // Constructor with model - assumed before cuts
  // Initial version does not do Lps

  CalQuadHeur (CbcModel & model);

  // Copy constructor
  CalQuadHeur (const CalQuadHeur &);

  // Destructor
  ~CalQuadHeur ();

  /// Clone
  virtual CbcHeuristic * clone() const
  {return new CalQuadHeur (*this);}

  /// Assignment operator
  CalQuadHeur & operator=(const CalQuadHeur& rhs);

  void setInstance (calInstance *inst);

  using CbcHeuristic::solution;

  /** returns 0 if no solution, 1 if valid solution.

      Sets solution values if good, sets objective value (only if good)

      We leave all variables which are at one at this node of the tree
      to that value and will initially set all others to zero.  We
      then sort all variables in order of their cost divided by the
      number of entries in rows which are not yet covered.  We
      randomize that value a bit so that ties will be broken in
      different ways on different runs of the heuristic.

      We then choose the best one and set it to one and repeat the exercise.
  */
  virtual int solution (double & objectiveValue,
			double * newSolution);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model) {}

  // generate initial point on subspace described by calibration vectors
  double *generateInitS (OsiSolverInterface &si);

  // cube method -- standalone: does not set all s to one or zero
  void standalone (double *s00);

protected:

  calInstance *instance_;
};

#endif
