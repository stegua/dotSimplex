/**
 * @fileoverview Copyright (c) 2019, XXXX YYYYY,
 *
 */

#pragma once

#include "Windows.h"
// respect this order for includes
#include "Psapi.h"

double getUsedRAM() {
  PROCESS_MEMORY_COUNTERS_EX pmc;
  GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc,
                       sizeof(pmc));
  double virtualMemUsedByMe = pmc.PrivateUsage;
  // Return size in MB
  return virtualMemUsedByMe / 1048576;
}

#define FEASIBILITY_TOL 1e-09
#define PRIC_TOL 1e-09

namespace DOT {

template <typename V = int, typename C = V>
class GVar {
 public:
  V a;  // First point
  V b;  // Second point
  C c;  // Distance

  GVar() : a(0), b(0), c(-1) {}

  GVar(V _a, V _b, C _c) : a(_a), b(_b), c(_c) {}
};

#include <vector>
typedef GVar<int, double> Var;
typedef std::vector<Var> Vars;

// Trying with flooat, faster speed, but too much numerical instability
typedef GVar<int, double> FVar;
typedef std::vector<FVar> FVars;

}  // end of namespace DOT