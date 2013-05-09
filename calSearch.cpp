/*
 * optimal calibrated sampling -- search method
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#include "calModel.hpp"
#include "calCube.hpp"
#include "calInstance.hpp"
#include "CoinTime.hpp"

#ifdef _MSC_VER
#define sprintf sprintf_s
#endif

//
// Repeatedly look for solution by applying bb and then quadratic
// solver
//
// Scheme:
//
// 1) if no initial point provided, generate one through the Cube
// variant in calCube.cpp. Initial point: s0
//
// For t:1..K 
//
//   2a) select a one in s0 and set it to zero, 
//   2b) open k(t) bound intervals and call BB
//   2c) fix integer variables and call quadratic solver
//   2d) save bound and solution if improved or if below threshold (1e-4)
//

#define N_MOVES    100   // number of moves
#define N_RESTART   20   // preferably divides N_MOVES; number of iterations before restart
#define FRACTION     0.2 // maximum % of elements unfixed
#define RANDOM_MOBILITY 1

inline double square (register double x)
{return (x > 1e40) ? x : (x * x);}

extern bool      GLOBAL_interrupt;
extern CbcModel *GLOBAL_curBB;

//#define DEBUG

bool calModel::search (CalCubeHeur &calCube, FILE *f, int repl) {

  int
    N = instance_ -> N (),
    n = instance_ -> n (),
    n_iter = (calInstance::GLOBAL == instance_ -> algType ()) ? 1 : instance_ -> nSolves ();

  double
    *bestSol = new double [1 + 2*N],
     bestObj = COIN_DBL_MAX;

  double *s0 = NULL;

  for (int nRetries = 0; !GLOBAL_interrupt && (nRetries < n_iter) && (square (bestObj) > instance_ -> eps ()); ++nRetries) {

    setMaximumSeconds
      (CoinMin (instance_ ->  maxTime       () < 0. ? COIN_DBL_MAX : instance_ -> maxTime (),
		instance_ ->  maxTotalTime  () < 0. ? COIN_DBL_MAX :
		(instance_ -> maxTotalTime  () - CoinCpuTime ()) /
		(instance_ -> nReplications () - repl)));

    calModel *b = clone ();

    OsiSolverInterface *si = b -> solver ();

    for (int i=0; i<N; ++i) {
      si -> setColLower (1 + N + i, 0.);
      si -> setColUpper (1 + N + i, 1.);
    }

    if (calInstance::GLOBAL == instance_ -> algType ()) {

      // do nothing: just run the BB

    } else if (0 == nRetries) {

      //
      // First step: generate initial solution, either from input
      // vector or through the Cube variant
      //
      // Generate initial point ON SUBSPACE. Does nothing if algtype
      // in {RANDOM, GLOBAL}
      //

      s0 = calCube. generateInitS (*si); 

      calCube. standalone (s0); // call Cube method

      b -> changeLU (*si, s0); // fixes some of the s variables after
                               // cube's flight phase

    } else { // change LU bounds based on current solution

#ifdef DEBUG
      printVec (s0, N, "s0");
#endif
      
      const double
	// reduce changes when random and solution available
	fraction    = FRACTION * (((calInstance::RANDOM == instance_ -> algType ()) && (instance_ -> d () [0] >= 0.)) ? RANDOM_MOBILITY : 1.), 
	multiplier  = fraction * (double) (N_RESTART - (nRetries % N_RESTART)) / N_RESTART,
	threshold_0 = multiplier * (double) (N-n) / N, // free half 0-fixed at beginning
	threshold_1 = multiplier * (double)    n  / N; // free half 1-fixed
      //threshold_0 = CoinMax (1. / CoinMax (1., (double) (N-n)), multiplier * (double) (N-n) / N), // free half 0-fixed at beginning
      //threshold_1 = CoinMax (1. / CoinMax (1., (double) n),     multiplier * (double)    n  / N); // free half 1-fixed

      bool one_changed = false;

      for (int nAttempts = 0; nAttempts < 20 && !one_changed; ++nAttempts)
	for (int i=0; i<N; ++i)
	  if      (s0 [i] <     1e-5) {if (drand48 () > threshold_0) si -> setColUpper (1 + N + i, 0.); else one_changed = true;}
	  else if (s0 [i] > 1 - 1e-5) {if (drand48 () > threshold_1) si -> setColLower (1 + N + i, 1.); else one_changed = true;}

#ifdef DEBUG
      char filename [40];
      sprintf (filename, "lp_%d_%d", repl, nRetries);
      si -> writeLp (filename);
#endif

      // int nChanges = CoinMax (1, n / nRetries); // <------------------ Number of changed variables
      // for (int i=0; i<nChanges;) {
      // 	int indChange = floor ((N - 1e-4) * drand48 ());
      // 	printf ("trying %d [%g,%g], i=%d\n", indChange, lb [1 + N + indChange], ub [1 + N + indChange], i);
      // 	if (ub [1 + N + indChange] <     1e-5) {si -> setColUpper (1 + N + indChange, 1.); ++i;}
      // 	if (lb [1 + N + indChange] > 1 - 1e-5) {si -> setColLower (1 + N + indChange, 0.); ++i;}
      // }
    }

    GLOBAL_curBB = b;

                             //    /|
                             //   / |--------+
    b -> branchAndBound ();  //  <  |        |
                             //   \ |--------+
                             //    \|

    //assert (fabs (b -> bestObj () - val [0]) < 1e-5);

    printf ("BB iteration %4d done (%10.2fs). ", 1+nRetries, CoinCpuTime ());

    //optimal = b -> isProvenOptimal(); 
    //const double *val = b -> getColSolution();

    //int numberColumns = model.getNumCols();

    if (b -> bestObj () < 1e20)
      printf ("Best solution: value %10.4f\n", square (b -> bestObj ()));
    else printf ("No solution found :-(\n");

    if (b -> bestSol ()) {

      if (b -> bestObj () < bestObj) {

	bestObj = b -> bestObj ();
	CoinCopyN (b -> bestSol (), 2*N+1, bestSol);
      }

      if (!s0)
	s0 = new double [N];

      CoinCopyN (b -> bestSol () + 1 + N, N, s0);
    }

    delete b;
  }

  bool 
    retval = true,
    row_format = (calInstance::ROW_BASED == instance_ -> outFormat ());

  delete [] s0;

  double w0 = (double) N / instance_ -> n ();

  if (bestObj < 1e20) {

    const double *val = bestSol;

    if (row_format) 
      fprintf (f, "%lf,", square (bestObj));

    //double *d = instance_ -> d ();

    printf ("sample:\n");

    for (int j=N+1, np=0; j<=2*N; ++j) {

      char strNum [5];

      sprintf (strNum, "%d", j-N);

      if (!row_format) 

	fprintf (f, "%d,%g,%s,%d,%g,\n",
		 1+repl,
		 square (bestObj),
		 instance_ -> id (j-N-1) ? instance_ -> id (j-N-1) : strNum,
		 (fabs (val [j]) > 1e-6) ? 1 : 0,
		 w0 + val [j-N]);

      if (fabs (val [j]) > 1e-6) {

	if (row_format) 
	  fprintf (f, "1,");

	printf     ("%d ",  j-N);
	if (!(++np % 10))
	  printf ("\n");
      } else 
	if (row_format) 
	  fprintf (f, "0,");
    }

    printf ("\nweights:\n");

    for (int j=1, np=0; j<=N; ++j) {

      if (row_format) 
	fprintf (f, "%g,", w0 + val [j]);

      if (fabs (val [j]) > 1e-6) {

	printf    ("[%d,%g] ", j, w0 + val [j]);
	if (!(++np % 10))
	  printf ("\n");
      }
    }

    if (square (bestObj) > instance_ -> eps ()) {

      printf ("(warning: solution has large objective at replication %d)", 1 + repl);
      //retval = false;
    }

    printf ("\n");

    if (row_format) 
      fprintf (f, "\n");

    fflush (f);
  }

  delete [] bestSol;

  return retval;
}
