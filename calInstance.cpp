/*
 * optimal calibrated sampling -- description of an instance
 *
 * (C) Pietro Belotti 2013. This code is released under the Eclipse
 * Public License.
 */

#include <stdio.h>
#include <string.h>

#include <vector>

#include "math.h"

#include "calInstance.hpp"
#include "CoinHelperFunctions.hpp"
#include <CoinTime.hpp>

#define MAX_LINE_LENGTH (1<<14)

#ifdef _MSC_VER
#define sscanf sscanf_s
#define fscanf fscanf_s
#endif

#define EPS_DEFAULT 1

//#define DEBUG

// reads and clean string from file

int get_string (FILE *f, char *line) {

  char *retval = fgets (line, MAX_LINE_LENGTH, f);

  //printf ("read line: %s", line);

  if (!retval)
    return 1;

  char *commentstart = strchr (line, '#');

  if (commentstart) 
    CoinFillN (commentstart, (int) (strlen (line) - (commentstart - line)), (char) 0);

  int 
    base = 0,
    len  = commentstart ? commentstart - line : strlen (line);

  bool prevWasSpace = true;

  for (int cursor = 0; cursor < len; ++cursor) {

    if (((isspace (line [cursor]) != 0) && !prevWasSpace) || (!(isspace (line [cursor])))) {

      prevWasSpace = (isspace (line [cursor]) != 0);

      if (isprint (line [cursor]))
	line [base++] = (isspace (line [cursor])) ? ' '  : line [cursor];
    }
  }

  line [base] = 0;

  //printf ("translated into: [%s] --> (%c), new length: %d\n", line, line [0], base);

  return 0;
}

//
// Constructor
//

calInstance::calInstance (char *filename):

  name_       (filename),
  d_          (NULL),
  id_         (NULL),
  eps_        (-1),
  maxIt_      (-1),   // No limit, BB will terminate upon finding optimal 
  maxBB_      (-1),   // solution in all iterations
  maxTime_    (-1),   // 
  maxTotTime_ (-1),   // max total time is also infinity
  nRepl_      (1),
  randSeed_   (-1),
  algType_    (CUBE),
  earlyStop_  (-1),
  nSolves_    (100),
  outFile_    (NULL),
  outFormat_  (ROW_BASED) {

#ifndef _MSC_VER
  FILE *f = fopen (filename , "r");
#else
  FILE *f;
  fopen_s (&f, filename, "r");
#endif

  //  enum {FIRST_LINE = 0, INITIAL_WEIGHTS, CALIB_VECTORS, DONE} mode = FIRST_LINE;

  if (!f) {

    fprintf (stderr, "Could not open file %s\n", filename);
    exit (-1);
  }

  N_ = n_ = p_ = -1;

  char line   [MAX_LINE_LENGTH];
  char idname [MAX_LINE_LENGTH];

  int 
    **curInd = NULL,
     *curNZ  = NULL;

  double 
    **curX = NULL;

  while (!(get_string (f, line))) {

    // next up should be a letter. Seek it

    //printf ("READ CHAR %c\n", *line);

    if ((N_ > 0) && (p_ >= 0) && (d_ == NULL)) {

      d_ = new double [N_];
      CoinFillN (d_, N_, -1.); // filled with "uninitialized" red flags

      curX   = new double * [1+p_];
      curInd = new int    * [1+p_];
      curNZ  = new int      [1+p_];

      CoinFillN (curNZ, p_, 0);

      for (int i=0; i <= p_; ++i) {
	curX   [i] = new double [N_];
	curInd [i] = new int    [N_];
      }

      // auxiliary vectors

      X_ = new CoinPackedVector * [1 + p_]; // extra is for cardinality constraint $\sum_{i=1}^N s_i = n$

      for (int i=0; i<=p_; ++i)
	X_ [i] = new CoinPackedVector;

      id_ = new char * [N_];

	  CoinFillN (id_, N_, (char *) NULL);
    }

    switch (line [0]) {

    case 'N': sscanf (line + 1, "%d",  &N_);          break;
    case 'n': sscanf (line + 1, "%d",  &n_);          break;
    case 'p': sscanf (line + 1, "%d",  &p_);          break;
    case 'e': sscanf (line + 1, "%lf", &eps_);        break;
    case 'i': sscanf (line + 1, "%d",  &maxIt_);      break;
    case 'k': sscanf (line + 1, "%d",  &nSolves_);    break;
    case 'b': sscanf (line + 1, "%d",  &maxBB_);      break;
    case 't': sscanf (line + 1, "%lf", &maxTime_);    break;
    case 'T': sscanf (line + 1, "%lf", &maxTotTime_); break;
    case 'R': sscanf (line + 1, "%d",  &nRepl_);      break;
    case 's': sscanf (line + 1, "%d",  &randSeed_);   break;
    case 'f': sscanf (line + 1, "%lf", &earlyStop_);  break;
    case 'o': outFile_ = (char *) malloc (1024 * sizeof (char));
      sscanf         (line + 1, "%s",  outFile_);     break;
    case 'O': sscanf (line + 1, "%s",  idname); 
      if (!(strcmp (idname, "block"))) outFormat_ = REPL_BLOCKS; break;

    case 'a': 

      sscanf (line + 1, "%s", idname);

      if      (!(strcmp (idname, "rand")))   algType_ = RANDOM;
      else if (!(strcmp (idname, "cube")))   algType_ = CUBE;
      else if (!(strcmp (idname, "global"))) algType_ = GLOBAL;
      else {
	printf ("algorithm \"%s\" not recognized.\nMust be one of \"rand\". \"cube\", or \"global\".\nExiting.\n", idname); 
	exit (-1);
      }

      break;

    // case 'w': // weights initialization mode
    //   sscanf (line + 1, "%s", idname);
    //   if      (!(strcmp (idname, "init")))  initType_ = INITWEIGHTS;
    //   else if (!(strcmp (idname, "lpsol"))) initType_ = LP_VALUE;
    //   else                                  initType_ = N_ON_n; // the default
    //   break;

    case 'y': // initial weight (optional -- default is n/N for all d_i)
    case 'I': // id names
    case 'x': // calibration variables

      {
	char curData = *line;

	int 
	  curpos = 1,
	  nfields = (curData == 'x' ? N_ * p_ : N_); // have to read different # of data if 'x'

	if (N_ < 0 || p_ < 0) {

	  printf ("Error: please define N and p before adding any data.\nExiting.");
	  exit (-1);
	}

	for (int i=0; i<nfields; ++i) {

	  //printf ("scanning [%s] (%d left)\n", line+curpos, nfields-i);

	  while ((strlen (line + curpos) == 0) ||
		 (sscanf (line + curpos, "%s ", idname) == 0) || 
		 (strspn (line + curpos, " ") == strlen (line + curpos))) {

	    curpos = 0;
	    if (get_string (f, line)) {
	      printf ("Not enough data in '%c' lines (read %d, need %d).\nExiting.\n", curData, i, nfields);
	      exit (-1);
	    }
	  }

	  //printf ("read idname: [%s] from remainder of string [%s]\n", idname, line + curpos);

	  if      (curData == 'I') {if (id_) {id_ [i] = new char [1 + strlen (idname)]; strcpy (id_ [i], 1 + strlen (idname), idname);}}
	  else if (curData == 'y') d_ [i] = atof (idname);
	  else if (curData == 'x') {

	    double elem = atoi (idname); 

	    if (fabs (elem) > 1e-6) {
	      curX   [i % p_] [curNZ [i % p_]   ] = elem;
	      curInd [i % p_] [curNZ [i % p_] ++] = i / p_;
	    }
	  }

	  curpos += (strlen (idname) + 1);
	}
      }

      break;

    default: if (isalnum (*line)) {
	printf ("Option '%c' not recognized.\nExiting.\n", *line);
	exit (-1);
      }
    }
  }

  curNZ [p_] = N_;
  for (int i=0; i < N_; ++i) {
    curInd [p_] [i] = i;
    curX   [p_] [i] = (double) n_ / N_;
  }

  for (int i=0; i <= p_; ++i)
    X_ [i] -> assignVector (curNZ [i], curInd [i], curX [i]);

  //X_ [p_] -> assignVector (N_, curInd [p_], curX [p_]);

  delete [] curNZ;
  delete [] curInd;
  delete [] curX;

  if (eps_ < 0)
    eps_ =  (GLOBAL == algType_) ? 0 : EPS_DEFAULT;

  if (randSeed_ < 0)
    randSeed_ = (int)(time (NULL)); 

  fclose (f);

  // cases in which to bail out:
  //
  // r>1 and initial solution given
  //
  // 
  //print ();
}

//
// destructor
//

calInstance::~calInstance () {

  delete [] d_;

  for (int i=0; i <= p_; ++i)              delete    (X_  [i]);
  for (int i=0; i <  N_; ++i) if (id_ [i]) delete [] (id_ [i]);

  delete [] X_; 
  delete [] id_;
}

//
// outputs a summary of the instance
//

void calInstance::print () {

  std::cout << "Name: " << name_ << std::endl;

  printf ("N=%d n=%d p=%d\n", N_, n_, p_);

#ifdef DEBUG
  printf ("Weights: ");
  for (int i=0; i < N_; ++i)
    printf ("%g ", d_ [i]);

  printf ("\nCalibration vectors: ");

  for (int i=0; i <= p_; ++i) {
    printf ("%d: ", i);
    for (int j=0; j < X_ [i] -> getNumElements(); ++j)
      printf ("(%d,%g) ", X_ [i] -> getIndices () [j], X_ [i] -> getElements () [j]);
    printf ("\n");
  }

  printf ("population: ");
  for (int i=0; i<N_; ++i)
    if (id_ [i] == NULL) printf ("? "); 
    else                 printf ("%s ", id_ [i]);

  printf ("\n");
#endif

  printf ("Options (everything else set to default):\n");

  if (eps_        != EPS_DEFAULT) printf ("Epsilon: %g\n",                   eps_);
  if (nRepl_      >  1)           printf ("Number of replications: %d\n",    nRepl_);
  if (maxTime_    >= 0)           printf ("CPU time limit: %g\n",            maxTime_);
  if (maxTotTime_ >= 0)           printf ("Total time allotted: %g\n",       maxTotTime_);
  if (maxBB_      >= 0)           printf ("BB nodes limit: %d\n",            maxBB_);
  if (nSolves_    >= 0)           printf ("Solutions per replication: %d\n", nSolves_);

  printf                                 ("Random seed: %d\n",               randSeed_);
}
