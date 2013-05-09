/*
 * optimal calibrated sampling -- Cube heuristic (Tillé and De Ville)
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"

#include "calInstance.hpp"
#include "calCube.hpp"
#include "calModel.hpp"

//#define DEBUG

// Default Constructor
CalCubeHeur::CalCubeHeur (calInstance *inst): 
  CbcHeuristic (), 
  noRun_    (false),
  instance_ (inst),
  calmodel_ (NULL) {}

// Constructor from model
CalCubeHeur::CalCubeHeur (CbcModel & model): 
  CbcHeuristic (model),
  noRun_    (false),
  instance_ (NULL),
  calmodel_ (NULL) {}

// Destructor
CalCubeHeur::~CalCubeHeur () {}

// Copy constructor
CalCubeHeur::CalCubeHeur (const CalCubeHeur & rhs): 
  CbcHeuristic (rhs), 
  noRun_    (rhs.noRun_),
  instance_ (rhs.instance_),
  calmodel_ (rhs.calmodel_) {}

// Assignment operator
CalCubeHeur &CalCubeHeur::operator= (const CalCubeHeur & rhs) {

  CbcHeuristic::operator=(rhs);
  noRun_    = rhs.noRun_;
  instance_ = rhs.instance_;
  calmodel_ = rhs.calmodel_;

  return *this;
}

// set instance pointer (useful in constructor with CbcModel argument)
void CalCubeHeur::setInstance (calInstance *inst) 
{instance_ = inst;}

extern bool      GLOBAL_interrupt;
extern CbcModel *GLOBAL_curBB;

// Returns 1 if solution, 0 if not
int CalCubeHeur::solution (double & solutionValue,
			   double * betterSolution) {

  if (noRun_ || GLOBAL_interrupt)
    return 0;

  calModel *b = calmodel_ -> clone ();

  OsiSolverInterface *si = b -> solver ();

  // disable cube heuristic in this, otherwise there is an infinite
  // recursive call

  for (int i = 0; i < b -> numberHeuristics (); i++) {

    CalCubeHeur *heur = dynamic_cast <CalCubeHeur *> (b -> heuristic (i));

    if (heur)
      heur -> setNoRun ();
  }

  int N = instance_ -> N ();

  for (int i=0; i<N; ++i) {
    si -> setColLower (1 + N + i, 0.);
    si -> setColUpper (1 + N + i, 1.);
  }

  // for (int i=0; i<N; ++i) {
  //   si -> setColLower (1 + i, -1000.);
  //   si -> setColUpper (1 + i,  1000.);
  // }

  double *s0 = generateInitS (*si); // generate initial point ON SUBSPACE

  standalone (s0);  // call Cube method

  b -> changeLU (*si, s0); // fixes some of the s variables after
                           // cube's flight phase

  delete [] s0;

  b -> messageHandler () -> setLogLevel (0);

  GLOBAL_curBB = b;

                           //    /|
                           //   / |--------+
  b -> branchAndBound ();  //  <  |        |
                           //   \ |--------+
                           //    \|

  GLOBAL_curBB = calmodel_;

  int retval = 0;

  if ((b -> bestObj () < 1e20) && (b -> bestSol ())) {

    //printVec (b -> bestSol () + 1, N, "sol");

    solutionValue = b -> bestObj ();
    CoinCopyN (b -> bestSol (), 2*N+1, betterSolution);

    retval = 1;
  } 

  delete b;

  return retval;
}

// generate initial point on subspace described by calibration vectors
double *CalCubeHeur::generateInitS (OsiSolverInterface &si) {

  int N = instance_ -> N ();

  double *s00 = new double [N];

  // If required by option, generate initial point at the LP
  // solution's s variables or as specified in input file, otherwise
  // just pick s00 = (n/N)_{i=1..N}
  //
  // The only constraint is that the vector be in [0,1]^N.

  switch (instance_ -> algType ()) {

  case calInstance::RANDOM:

    if (instance_ -> d () [0] >= 0.) // there is an initial solution
      CoinCopyN (instance_ -> d (), N, s00);
    else {

      int n = instance_ -> n ();

      // generate a vector with n ones and N-n zeros, randomly placed.

      double maj = (n > N/2) ? 1. : 0.; // makes it more efficient to
					// fill in if zeros fewer than ones

      CoinFillN (s00, N, maj);

      for (int i=0; i<n;) {

	int pos = (int) (drand48 () * ((double) N - 1e-5));

	if (fabs (maj - s00 [pos]) < .45) { // unfixed just yet
	  s00 [pos] = 1 - maj;
	  ++i;
	}
      }
    }

    break;

    // // never used, or not required
    // case calInstance::LP_VALUE: // copy LP values of s to initial point
    //   si.initialSolve ();
    //   CoinCopyN (si.getColSolution () + N + 1, N, s00);
    //   break;

  case calInstance::CUBE:
  case calInstance::GLOBAL:

    CoinFillN (s00, N, (double) instance_ -> n () / N);
    break;

  default:

    printf ("Option for initial point not recognized.\nExiting.\n");
    exit (-1);
  }

#ifdef DEBUG
  printVec (s00, N, "s00");
#endif

  return s00;
}
