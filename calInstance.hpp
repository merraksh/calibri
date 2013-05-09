/*
 * optimal calibrated sampling - definition of instance
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#ifndef calInstance_hpp
#define calInstance_hpp

#include <set>
#include <map>
#include <string>

#include <CoinPackedVector.hpp>
#include <CoinPackedMatrix.hpp>

#ifdef _MSC_VER
#define strcpy strcpy_s
#else
#define strcpy(a,b,c) strncpy(a,c,b)
#endif

///
/// Instance class. Packs all info contained in input file, plus some
/// flags specified at the command line.
///

class calInstance {

  friend int main (int, char **); 

public:

  enum AlgType   {RANDOM, CUBE, GLOBAL};
  enum OutFormat {ROW_BASED, REPL_BLOCKS};

protected:

  std::string        name_;       ///< instance name
  int                N_;          ///< cardinality of the population
  int                n_;          ///< size of the sample (<N_)
  int                p_;          ///< number of calibration vectors
  double            *d_;          ///< initial weight vector
  CoinPackedVector **X_;          ///< calibration vectors
  char             **id_;         ///< id of each element of the population
  double             eps_;        ///< terminate if ||w - d|| <= eps
  int                maxIt_;      ///< max # iterations
  int                maxBB_;      ///< max # branch-and-bound nodes
  double             maxTime_;    ///< max time
  double             maxTotTime_; ///< max total time over all replications
  int                nRepl_;      ///< # of replications
  int                randSeed_;   ///< random seed
  enum AlgType       algType_;    ///< algorithm type
  double             earlyStop_;  ///< stop Flight phase at this * (N-p-1) iterations. Default: 1
  int                nSolves_;    ///< solve this many problems before giving up
  char              *outFile_;    ///< filename for output
  enum OutFormat     outFormat_;  ///< output format

public:

  calInstance (char *filename); ///< constructor from argument list
  ~calInstance ();              ///< destructor

  // get () methods

  int N () {return N_;}
  int n () {return n_;}
  int p () {return p_;}

  CoinPackedVector** &X () {return X_;}
  double*            &d () {return d_;}

  std::string &name () {return name_;}

  char   *id             (int ind) {return id_ [ind];}
  double  eps            ()        {return eps_;}
  int     maxIterations  ()        {return maxIt_;}
  int     maxBBnodes     ()        {return maxBB_;}
  double  maxTime        ()        {return maxTime_;}
  double  maxTotalTime   ()        {return maxTotTime_;}
  int    &nReplications  ()        {return nRepl_;}
  int    &randSeed       ()        {return randSeed_;}
  double  earlyStop      ()        {return earlyStop_;}
  int    &nSolves        ()        {return nSolves_;}

  enum AlgType   &algType   ()     {return algType_;}
  enum OutFormat &outFormat ()     {return outFormat_;}

  void print ();
};

#endif
