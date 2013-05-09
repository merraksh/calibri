/*
 * optimal calibrated sampling -- Cube heuristic (Tillé and De Ville)
 *
 * (C) Pietro Belotti 2013. This code is released under the Eclipse
 * Public License.
 */

#ifndef calCube_H
#define calCube_H

#ifdef WIN32
#define srand48(a) srand((a))
#define drand48() ((double)rand()/RAND_MAX)
#endif

#include "CbcHeuristic.hpp"

class calInstance;
class calModel;

//
// Heuristic to run a variant of the Cube algorithm
//

class CalCubeHeur: public CbcHeuristic {

public:

  // Default Constructor
  CalCubeHeur (calInstance *i);

  // Constructor with model - assumed before cuts
  // Initial version does not do Lps

  CalCubeHeur (CbcModel & model);

  // Copy constructor
  CalCubeHeur (const CalCubeHeur &);

  // Destructor
  ~CalCubeHeur ();

  /// Clone
  virtual CbcHeuristic * clone() const
  {return new CalCubeHeur (*this);}

  /// Assignment operator
  CalCubeHeur & operator=(const CalCubeHeur& rhs);

  void setInstance (calInstance *inst);

  using CbcHeuristic::solution ;

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

  void setCalModel (calModel *m)
  {calmodel_ = m;}

  void setNoRun ()
  {noRun_ = true;}

protected:

  bool noRun_;

  calInstance *instance_;
  calModel    *calmodel_;

  // project v on null space of restricted calibration constraints. If
  // pi=NULL, taken to be with all elements not in {0,1}
  void project (double *v, double *u, double *pi = NULL);
};

//
// Print vectors and matrices (for debugging purposes)
//

void printVec        (const double *a, int size,             const std::string &name);
void printMatr       (const double *A, int nrows, int ncols, const std::string &name);
void printMatrTransp (const double *A, int nrows, int ncols, const std::string &name);

#endif
