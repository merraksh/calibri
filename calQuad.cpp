/*
 * optimal calibrated sampling -- quadratic solver
 *
 * (C) Pietro Belotti 2013. This code is released under the Eclipse
 * Public License.
 */

#include "CoinMpsIO.hpp"
#include "ClpInterior.hpp"
#include "ClpSimplex.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpQuadraticObjective.hpp"

// quadratic method -- standalone: start from binary (?) vector and
// obtain another vector by repeated moves (0,1) -> (1,0) on (i,j)
// randomly chosen.

void CalQuadHeur::standalone (double *s0) {

  /* Read quadratic model in two stages to test loadQuadraticObjective.

     And is also possible to just read into ClpSimplex/Interior which sets it all up in one go.
     But this is only if it is in QUADOBJ format.

     If no arguments does share2qp using ClpInterior (also creates quad.mps which is in QUADOBJ format)
     If one argument uses simplex e.g. testit quad.mps
     If > one uses barrier via ClpSimplex input and then ClpInterior borrow
  */

  CoinMpsIO  m;

  ClpInterior model;
  model.loadProblem(*m.getMatrixByCol(), m.getColLower(), m.getColUpper(),
		    m.getObjCoefficients(),
		    m.getRowLower(), m.getRowUpper());

  // get quadratic part
  int * start = NULL;
  int * column = NULL;
  double * element = NULL;

  m.readQuadraticMps (NULL, start, column, element, 2);

  int j;
  for (j = 0; j < 79; j++) {
    if (start[j] < start[j+1]) {
      int i;
      printf("Column %d ", j);
      for (i = start[j]; i < start[j+1]; i++) {
	printf("( %d, %g) ", column[i], element[i]);
      }
      printf("\n");
    }
  }

  model.loadQuadraticObjective (model.numberColumns (), start, column, element);

  model.writeMps("quad.mps");

  ClpCholeskyBase * cholesky = new ClpCholeskyBase();

  cholesky->setKKT(true);
  model.setCholesky (cholesky);
  model.primalDual();

  double
    *primal = model.primalColumnSolution (),
    *dual   = model.dualRowSolution      ();

  int i;
  int numberColumns = model.numberColumns();
  int numberRows    = model.numberRows();

} else {

  // Could read into ClpInterior
  ClpSimplex model;


  if (argc < 3) {
    // simplex - just primal as dual does not work
    // also I need to fix scaling of duals on output
    // (Was okay in first place - can't mix and match scaling techniques)
    // model.scaling(0);
    model.primal();
  } else {
    // barrier
    ClpInterior barrier;
    barrier.borrowModel(model);
    ClpCholeskyBase * cholesky = new ClpCholeskyBase();
    cholesky->setKKT(true);
    barrier.setCholesky(cholesky);
    barrier.primalDual();
    barrier.returnModel(model);
  }
 }
return 0;
}
