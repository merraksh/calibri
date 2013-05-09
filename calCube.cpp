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

#define _USE_MATH_DEFINES
#include <cmath>

#include <stdlib.h>

#include <cassert>
#include <cstdlib>
#include <cfloat>

#include "OsiSolverInterface.hpp"

#include "calInstance.hpp"
#include "calCube.hpp"

extern bool GLOBAL_interrupt;

// cube method -- standalone: does not set all s to one or zero
void CalCubeHeur::standalone (double *s0) {

  // s00 is projected onto the subspace of the calibrating vectors --
  // actually that might give a point outside of [0,1]^N, so don't
  //
  // flight phase: repeat
  //   - generate v random in [-1,1]^N
  //   - project v onto null space of A -> u
  //   - if u=0, break;
  //   - compute lambdaMinus and lambdaPlus
  //   - move s0 using lambdas and u
  // forever
  //

  int N = instance_ -> N ();

  double
    *v  = new double [N],
    *u  = new double [N];

  for (int iter = 0; !GLOBAL_interrupt; ++iter) {

    //printf ("iteration %d: ", iter);

    // generate v
    for (int i=N; i--;)
      v [i] = -1 + 2 * drand48 (); // returns random vector in [-1,1]^N

    // Better generator: generate random point on unit sphere

    double cos_seq = 1;

    for (int i=N; --i;) {

      // angle generated uniformly in [-pi/2,pi/2]
      register double alpha = -M_PI / 2. + M_PI * drand48 ();
      v [i] = cos_seq * sin (alpha);
      cos_seq *=        cos (alpha);
    }

    v [0] = cos_seq * ((drand48 () < .5) ? -1. : 1.);

    //printVec (s0,N,"s0");
    //printVec (v, N,"random v");

    project (v, u, s0); // project v on restricted null space of
			// calibration constraints

    bool nonzero = false;

    for (int i=N; i--;)
      if (fabs (u [i]) > 1e-5) {nonzero = true; break;}

    if (((        instance_ -> earlyStop () >= 0) && 
	 (iter >= instance_ -> earlyStop () * (N - instance_ -> p ()))) || // bail out when a sufficient number of number has been fixed
	!nonzero)                                                          // can't move from this vertex of C intersect
      break; 			                                           // Q, we are done

    // find how far we can go from current point ////////////////////////////////////////////

    double 
      lambdaM = COIN_DBL_MAX,
      lambdaP = COIN_DBL_MAX;

    for (int i=0; i < N; ++i) 

      if (fabs (u [i]) > 1e-6) {

	lambdaM = CoinMin (lambdaM, 
			   u [i] > 0 ?      s0 [i]  / u [i] : 
			   u [i] < 0 ? (1 - s0 [i]) / u [i] : lambdaM);

	lambdaP = CoinMin (lambdaP, 
			   u [i] > 0 ? (1 - s0 [i]) / u [i] : 
			   u [i] < 0 ?    - s0 [i]  / u [i] : lambdaP);
      }

    // now modify s0: up with probability lambdaM / (lambdaM + lambdaP), down otherwise //////

    if (drand48 () < lambdaM / (lambdaM + lambdaP)) for (int i=N; i--;) s0 [i] += lambdaP * u [i]; // move up
    else                                            for (int i=N; i--;)	s0 [i] -= lambdaM * u [i]; // move down

    //printVec (s0, N, "\n\ns");
  }

  delete [] v;
  delete [] u;
}

//
// Debug prints
//

void printVec (const double *a, int size, const std::string &name) {
  printf ("%s: [", name.c_str()); fflush (stdout); 
  for (int i=0; i<size; ++i) {
    printf ("%g ", a [i]); 
    fflush (stdout);
  } 
  printf ("]\n");
}

void printMatr (const double *A, int nrows, int ncols, const std::string &name) {
  printf ("%s:\n", name.c_str()); 
  for   (int i=0; i<nrows; ++i) {
    for (int j=0; j<ncols; ++j) {
      printf ("%g ", A [j*nrows+i]);
      fflush (stdout);
    }
    printf ("\n");
  }
}

void printMatrTransp (const double *A, int nrows, int ncols, const std::string &name) {
  printf ("%s (tr):\n", name.c_str()); 
  for   (int j=0; j<ncols; ++j) {
    for (int i=0; i<nrows; ++i) {
      printf ("%g ", A [j*nrows+i]); 
      fflush (stdout);
    }
    printf ("\n");
  }
}
