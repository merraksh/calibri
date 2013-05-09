/*
 * optimal calibrated sampling -- main
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>

#include <string.h>

#include <OsiClpSolverInterface.hpp>
#include <CoinPackedVector.hpp>

#include <CbcNode.hpp>
#include <CoinTime.hpp>
#include <CbcModel.hpp>
#include <OsiAuxInfo.hpp>
#include <CbcCutGenerator.hpp>

#include "calInstance.hpp"
#include "calModel.hpp"
#include "calCut.hpp"
#include "calBT.hpp"
#include "calCube.hpp"
#include "cmdLine.hpp"

//#define DEBUG

//#define MAX_GAP 1e-1 // value at which to stop branch-and-bound

//
// Fill in LP's coefficient
//

int populate (calInstance *instance, OsiSolverInterface *problem);

//
// Add cutting planes, heuristics, etc.
//

void addCbcExtras (calModel &calbb, int &cutGenCount);


/// global variable
bool      GLOBAL_interrupt = false;
CbcModel *GLOBAL_curBB     = NULL;

#define INTERRUPT_HANDLER

#ifdef  INTERRUPT_HANDLER

#include "CoinSignal.hpp"

extern "C" {

  static void signal_handler (int sig) {

    if (GLOBAL_interrupt) {

      std::cerr << "[BREAK]" << std::endl;
      exit (-1);

    } else {

      GLOBAL_interrupt = true;
      if (GLOBAL_curBB) {
	GLOBAL_curBB -> setMaximumNodes   (0); // stop at next node
	GLOBAL_curBB -> setMaximumSeconds (0);
      }
    }

    return;
  }
}

#endif


//                     oo          
//                                
// 88d8b.d8b. .d8888b. dP 88d888b. 
// 88'`88'`88 88'  `88 88 88'  `88 
// 88  88  88 88.  .88 88 88    88 
// dP  dP  dP `88888P8 dP dP    dP 

int main (int argc, char *argv[]) {

  if (argc <= 1) {
    printf ("Usage: %s [options] <instance.txt>\nRun \"%s -h\" for help\n", argv [0], argv [0]);
    exit (0);
  }

#ifdef INTERRUPT_HANDLER
  signal (SIGINT, signal_handler);
#endif

  /*
   * Specify program command line options
   */

  char **filenames;

  bool needHelp = false;

  tpar options [] = {{ 'n', (char *) "sample-size",     1, NULL,    ::TINT,    (char *) "sample size (n)"}
		     ,{'e', (char *) "epsilon",        -1, NULL,    ::TDOUBLE, (char *) "stop BB optimization when objective below number"}

		     ,{'i', (char *) "lp-iter",        -1, NULL,    ::TINT,    (char *) "LP iterations in each BB run"}
		     ,{'k', (char *) "iterations",    100, NULL,    ::TINT,    (char *) "number of BB runs in each replication"}
		     ,{'b', (char *) "bb-nodes",       -1, NULL,    ::TINT,    (char *) "number of subproblems in each BB run"}
		     ,{'t', (char *) "time",           -1, NULL,    ::TDOUBLE, (char *) "CPU time allotted to each BB run"}
		     ,{'T', (char *) "tot-time",       -1, NULL,    ::TDOUBLE, (char *) "maximum total CPU time"}

		     ,{'R', (char *) "replications",    1, NULL,    ::TINT,    (char *) "number of replications (if \"-g\" chosen, it is ignored)"}
		     ,{'s', (char *) "seed",           -1, NULL,    ::TINT,    (char *) "random seed (if not specified here or in input file, it is generated using time)"}
		     ,{'f', (char *) "int-fixed",       1, NULL,    ::TDOUBLE, (char *) "portion of fixed variables before stopping flight phase of Cube"}
		     ,{'o', (char *) "output",          0, NULL,    ::TSTRING, (char *) "output file"}
		     ,{'O', (char *) "out-format",      0, NULL,    ::TSTRING, (char *) "output format: \"block\" for one unit per row, \"row\" (default) otherwise"}

		     ,{'r', (char *) "random",          0, NULL,    ::TTOGGLE, (char *) "generate random initial point (overrides \"-g\" and \"-c\")"}
		     ,{'c', (char *) "cube",            0, NULL,    ::TTOGGLE, (char *) "generate initial point through Cube"}
		     ,{'g', (char *) "global",          0, NULL,    ::TTOGGLE, (char *) "find global optimum (overrides \"-c\")"}

		     ,{'h', (char *) "help",            0, &needHelp, ::TTOGGLE, (char *) "print this help"}

		     ,{0,   (char *) "",                0, NULL,    ::TTOGGLE,       (char *) ""} /* THIS ENTRY ALWAYS AT THE END */
  };

  // default parameter values

  set_default_args (options);

  // parse command line

  filenames = readargs (argc, argv, options);

  if (needHelp) {

    printf ("\
%s -- extract a calibrated sample\n\
Synopsis: %s returns a set of samples of size n of a population U of size N.\n\
The samples are calibrated w.r.t. a set of p auxiliary variable vectors.\n\
Author: Pietro Belotti\n\
Version: 0.1\n\
License: Eclipse Public License\n\n", argv [0], argv [0]);

    print_help (argv [0], options);

    printf ("\n\
Input file format as below:\n\
\n\
N 10 # cardinality of population U\n\
n 3  # cardinality of sample S\n\
p 4  # number of auxiliary variables\n\
y 0 1 0 0 0 1 0 1 0 0 # initial solution (optional)\n\
x 1 1 0 0 # auxiliary variables: N rows of p columns each\n\
  0 0 1 0\n\
  0 1 0 0\n\
  0 0 0 1\n\
  1 0 0 0\n\
  0 1 0 0\n\
  0 0 0 0\n\
  0 0 1 1\n\
  1 0 1 0\n\
  0 1 0 0\n\
I Pisa Parma Torino Roma Venezia Napoli Bari Palermo Bologna Milano\n\
s 734642 # random seed (if -1 then generated using time)\n\
\n\
Output file format\nWithout option \"--out-format block\": for each replication, one line\n\
containing the squared norm of (w-d), 0-1 vector with sample, and vector of weights.\n\
\nWithout option \"--out-format block\": for each replication, one block of N lines.\n\
Each line contains replication, function value, id, 0/1, and weight of each unit.\n\
See user manual for mode details on input and output file formats.\n");

    exit (0);
  }

  double nowTime = CoinCpuTime ();
  calInstance *instance = NULL;

  printf ("Calibri -- a solver for the optimal calibrated sampling problem\n");
  printf ("Reading instance %s: ", filenames ? *filenames : "from std input"); fflush (stdout);

  instance = new calInstance (*filenames);
  printf ("done (%.3gs)\n", CoinCpuTime () - nowTime); fflush (stdout);

  options  [0].par =  &(instance -> n_);
  options  [1].par =  &(instance -> eps_);
                        
  options  [2].par =  &(instance -> maxIt_);
  options  [3].par =  &(instance -> nSolves_);
  options  [4].par =  &(instance -> maxBB_);
  options  [5].par =  &(instance -> maxTime_);
  options  [6].par =  &(instance -> maxTotTime_);
                        
  options  [7].par =  &(instance -> nRepl_);
  options  [8].par =  &(instance -> randSeed_);
  options  [9].par =  &(instance -> earlyStop_);
  options [10].par =  &(instance -> outFile_);

  char *outFor = NULL;

  options [11].par = &outFor;

  bool
    isRandom = false,
    isCube   = false,
    isGlobal = false;

  options [12].par =  &isRandom;
  options [13].par =  &isCube;
  options [14].par =  &isGlobal;

  options [15].par =  &needHelp;

  options [16].par = NULL; // redundant -- to end it

  // delete filenames

  if (filenames) {
    for (int i=0; filenames [i]; ++i)
      free (filenames [i]);
    free (filenames);
  }

  // RE-READ options in order to override file-based options
  filenames = readargs (argc, argv, options);

  if (filenames) {
    for (int i=0; filenames [i]; ++i)
      free (filenames [i]);
    free (filenames);
  }

  if (outFor && (!(strcmp (outFor, "block"))))
    instance -> outFormat_ = calInstance::REPL_BLOCKS;

  instance -> algType_ = 
    isRandom ? calInstance::RANDOM :
    isGlobal ? calInstance::GLOBAL : 
               calInstance::CUBE;

  if ((instance -> d ()) && 
      ((instance -> algType () != calInstance::RANDOM) ||
       (instance -> nReplications () > 1))) {

    CoinFillN (instance -> d (), instance -> N (), -1.); // filled with "uninitialized" red flags

    //delete [] instance -> d ();

    //instance -> d () = NULL;
    printf ("Warning: initial solution ignored for this type of algorithm.\n\
Specify \"a rand\" in the input file or option \"-r\" at the command line and set number of replications to one.\n");
  }

  //
  // Sanity check of input:
  //
  // 1) if r < 1, exit
  // 2) if r > 1, d_ ignored
  // 3) if global, r > 1 ignored
  //

  if (instance -> nReplications () < 1) {
    printf ("No replications requested.\nExiting.\n");
    if (instance) delete instance;
    exit (0);
  }

  // if ((instance -> nReplications () > 1) && (instance -> d () [0] >= 0.)) {
  //   printf ("Warning: %d replications requested, initial solution ignored\n", instance -> nReplications ());
  //   CoinFillN (instance -> d (), instance -> N (), -1.);
  // }

  if (calInstance::GLOBAL == instance -> algType ()) {

    if (instance -> nReplications () > 1) {
      printf ("Warning: global solution requested, one replication will be run\n");
      instance -> nReplications () = 1;
    }

    if (instance -> d () [0] >= 0)
      printf ("Warning: global solution requested, initial solution ignored\n");
  }

  // if ((instance -> nReplications () > 1) && (calInstance::CUBE == instance -> algType ())) {
  //   printf ("Warning: %d replications and Cube algorithm specified, random seed ignored.\n", instance -> nReplications ());
  //   instance -> randSeed () = (int)(time (NULL)); 
  // }

  // Sanity check done ---------------------------------------------------

  instance -> print ();

  OsiClpSolverInterface model;

  printf ("Creating MILP: ");
  nowTime = CoinCpuTime ();
  populate (instance, &model);
  printf ("done (%gs)\n", CoinCpuTime () - nowTime);

  model. messageHandler () -> setLogLevel (0);

  calModel calbb (model, instance);

  calbb. messageHandler () -> setLogLevel ((calInstance::GLOBAL == instance -> algType ()) ? 1 : 0);

  calbb. setAllowableGap (instance -> eps ());

  int cgCnt = 0;

  addCbcExtras (calbb, cgCnt);

  calCut cutgen (instance);
  calbb. addCutGenerator (&cutgen, 1, "Conic cuts", true, false, false, 1);
  calbb. cutGenerator (cgCnt++) -> setGlobalCuts (true);

  calBT  btgen (instance, &calbb);
  calbb. addCutGenerator (&btgen, 1, "Bound Reduction", true, false, false, 1);
  calbb. cutGenerator (cgCnt++) -> setGlobalCuts (true);

  CalCubeHeur calCube (calbb);
  calCube. setCalModel (&calbb);
  calCube. setInstance (instance);
  calCube. setHeuristicName           ("Cube Method");

  calbb. addHeuristic                (&calCube); // not necessary for now

  // double objLimit;
  // model.getDblParam (OsiDualObjectiveLimit, objLimit);

  OsiBabSolver solverChar;
  solverChar. setSolverType (3);

  //calbb. solver           () -> setAuxiliaryInfo (&solverChar);
  //calbb. passInSolverCharacteristics             (&solverChar);
  //calbb. continuousSolver () -> setAuxiliaryInfo (&solverChar);

  if (instance -> maxIterations () >= 0) calbb.setMaximumNumberIterations (instance -> maxIterations ());
  if (instance -> maxBBnodes    () >= 0) calbb.setMaximumNodes            (instance -> maxBBnodes ());

  if (!(instance -> outFile_)) {
    instance -> outFile_ = (char *) malloc (sizeof (char) * strlen (argv [argc-1]) + 5);
    strcpy (instance -> outFile_, strlen (argv [argc-1]) + 1, argv [argc-1]);

    char *exthook = strstr (instance -> outFile_, ".txt");
    if (!exthook)
      exthook = instance -> outFile_ + strlen (instance -> outFile_);

    strcpy (exthook, 5, ".sol");
  }

  printf ("Writing file %s\n", instance -> outFile_);

  FILE *f;

#ifndef _MSC_VER
  f = fopen   (instance -> outFile_, "w");
#else
  fopen_s (&f, instance -> outFile_, "w");
#endif

  free (instance -> outFile_);

  printf ("Generating point(s) (%gs)\n", CoinCpuTime ());

  // if ((instance -> initType () == calInstance::LP_VALUE) &&   
  //     (instance -> nReplications () != 1)) {
  //   printf ("Initial weights set from LP but more than one repetition.\nResetting # repetitions to 1.\n");
  //   instance -> nReplications () = 1;
  // }

  //int nFails = 0;
  //#define MAX_FAILS 20

  // MAIN LOOP:

  if (instance -> outFormat_ == calInstance::ROW_BASED) {

    int N = instance -> N ();

    fprintf (f, "%d,", N);

    for (int times=2; times--;) // do this for both binary and weight vector

      for (int i=0; i<N; ++i)
	if (instance -> id (i)) fprintf (f, "%s,", instance -> id (i));
	else                    fprintf (f, "%d,", 1+i);

    fprintf (f, "\n");

  } else // first line of block
    fprintf (f, "R,F,ID,S,W,\n");

  srand48 (instance -> randSeed ());

  for (int iter=0; iter < instance -> nReplications (); ++iter) {

    if (GLOBAL_interrupt) {

      printf ("User interrupt\n");
      break;
    }

    printf ("-------------- Replication %d:\n", 1+iter);

    if (!(calbb. search (calCube, f, iter)))

      printf ("Warning: no solution found at this replication.\n");

    // if ( || (nFails > MAX_FAILS) || (calInstance::GLOBAL == instance -> algType ())) {
    //   ++iter;
    //   nFails = 0;
    // } else ++nFails;
  }

  fclose (f);
 
  if (instance) 
    delete instance;

  return 0;
}
