/*
 * optimal calibrated sampling -- add heuristics and cut generators
 *
 * (C) Pietro Belotti 2013. This code is released
 * under the Eclipse Public License.
 */

#include <CbcHeuristicLocal.hpp>
#include <CbcHeuristicPivotAndFix.hpp>
#include <CbcHeuristicRandRound.hpp>
#include <CbcHeuristicGreedy.hpp>
#include <CbcHeuristicFPump.hpp>
#include <CbcHeuristicRINS.hpp>
#include <CbcHeuristicDiveCoefficient.hpp>
#include <CbcHeuristicDiveFractional.hpp>
#include <CbcHeuristicDiveGuided.hpp>
#include <CbcHeuristicDiveVectorLength.hpp>
#include <CbcHeuristicDivePseudoCost.hpp>
#include <CbcHeuristicDiveLineSearch.hpp>

#include <CglProbing.hpp>
#include <CglPreProcess.hpp>
#include <CglCutGenerator.hpp>
#include <CglGomory.hpp>
#include <CglKnapsackCover.hpp>
#include <CglRedSplit.hpp>
#include <CglClique.hpp>
#include <CglFlowCover.hpp>
#include <CglMixedIntegerRounding2.hpp>
#include <CglTwomir.hpp>
#include <CglDuplicateRow.hpp>
#include <CglStored.hpp>
#include <CglLandP.hpp>
#include <CglResidualCapacity.hpp>

#include <CbcNode.hpp>
#include <CbcCutGenerator.hpp>

#include "calModel.hpp"
#include "calCube.hpp"

void addCbcExtras (calModel &calbb, int &cgCnt) {

  // CbcHeuristicFPump heuristicFPump (calbb);
  // heuristicFPump.setWhen           (13);
  // heuristicFPump.setMaximumPasses  (20);
  // heuristicFPump.setMaximumRetries (7);
  // heuristicFPump.setHeuristicName  ("feasibility pump");
  // heuristicFPump.setInitialWeight  (1);
  // heuristicFPump.setFractionSmall  (0.6);
  // calbb. addHeuristic               (&heuristicFPump);
  
  CbcRounding rounding      (calbb);
  rounding.setHeuristicName ("Rounding");
  calbb. addHeuristic        (&rounding);

  CbcHeuristicLocal heuristicLocal (calbb);
  heuristicLocal.setHeuristicName  ("Combine solutions");
  heuristicLocal.setSearchType     (1);
  heuristicLocal.setFractionSmall  (0.6);
  calbb. addHeuristic               (&heuristicLocal);

  // CbcHeuristicGreedyCover heuristicGreedyCover (calbb);
  // heuristicGreedyCover.setHeuristicName        ("Greedy cover");
  // calbb. addHeuristic                           (&heuristicGreedyCover);

  // CbcHeuristicGreedyEquality heuristicGreedyEquality (calbb);
  // heuristicGreedyEquality.setHeuristicName           ("Greedy equality");
  // calbb. addHeuristic                                 (&heuristicGreedyEquality);

  CglProbing probing;
  probing.setMaxProbe     (10);
  probing.setMaxLook      (10);
  probing.setMaxElements  (200);
  probing.setMaxProbeRoot (50);
  probing.setMaxLookRoot  (10);
  probing.setRowCuts      (3);
  probing.setUsingObjective(true);
  calbb. addCutGenerator (&probing, -1, "Probing", true, false, false, -100, -1, -1);
  calbb. cutGenerator    (cgCnt++) -> setTiming (true);

  CglGomory gomory;
  gomory.setLimitAtRoot (512);
  calbb. addCutGenerator (&gomory, -98, "Gomory", true, false, false, -100, -1, -1);
  calbb. cutGenerator    (cgCnt++)->setTiming (true);

  // CglKnapsackCover knapsackCover;
  // calbb. addCutGenerator (&knapsackCover, -98, "KnapsackCover", true, false, false, -100, -1, -1);
  // calbb. cutGenerator    (cgCnt++) -> setTiming (true);

  CglRedSplit redSplit;
  calbb. addCutGenerator (&redSplit, -99, "RedSplit", true, false, false, -100, -1, -1);
  calbb. cutGenerator    (cgCnt++) -> setTiming (true);
 
  // CglClique clique;
  // clique.setStarCliqueReport (false);
  // clique.setRowCliqueReport  (false);
  // clique.setMinViolation     (0.1);
  // calbb. addCutGenerator      (&clique, -98, "Clique", true, false, false, -100, -1, -1);
  // calbb. cutGenerator         (cgCnt++) -> setTiming (true);

  // CglMixedIntegerRounding2 mixedIntegerRounding2;
  // calbb. addCutGenerator (&mixedIntegerRounding2, -98, "MixedIntegerRounding2", true, false, false, -100, -1, -1);
  // calbb. cutGenerator    (cgCnt++) -> setTiming (true);

  // CglFlowCover flowCover;
  // calbb. addCutGenerator (&flowCover, -98, "FlowCover", true, false, false, -100, -1, -1);
  // calbb. cutGenerator    (cgCnt++) -> setTiming (true);

  CglTwomir twomir;
  twomir.setMaxElements (250);
  calbb. addCutGenerator (&twomir, -99, "Twomir", true, false, false, -100, -1, -1);
  calbb. cutGenerator    (cgCnt++) -> setTiming (true);
}
