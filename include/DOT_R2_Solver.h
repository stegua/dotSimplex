/**
 * @fileoverview Copyright (c) 2019, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <amp.h>
#include <amp_math.h>

#include <omp.h>

#include <cassert>
#include <chrono>
#include <cinttypes>
#include <limits>
#include <random>

#include <fstream>
#include <sstream>

#include <vector>
using std::vector;

#include <array>
using std::array;

#include "DOT_DotSimplex.h"
#include "DOT_Vars.h"

// Distance function between a pair of point in R^k
inline constexpr auto DISTANCE_R2(const double* x, const double* y) {
  return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]);
}

namespace DOT {
namespace R2 {

// Container for general discrete measure
template <typename FlowType = int, typename PosType = double>
class GMeasureR2 {
 public:
  GMeasureR2() {}

  // Read from file (e.g., DOTMark images)
  GMeasureR2(const std::string& filename) { readFromFile(filename); }

  // setter
  void reserve(size_t n) {
    Ws.reserve(n);
    Ps.reserve(2 * n);
  }

  void add(FlowType _w, PosType _p1, PosType _p2) {
    Ws.emplace_back(_w);
    Ps.emplace_back(_p1);
    Ps.emplace_back(_p2);
  }

  // Use as few memory as possible
  void shrink_to_fit() {
    Ws.shrink_to_fit();
    Ps.shrink_to_fit();
  }

  // getters
  size_t size() const { return Ws.size(); }

  FlowType getW(size_t i) const { return Ws[i]; }

  inline const PosType* getP(size_t i) const { return &Ps[2 * i]; }

  // Parse from file
  void readFromFile(const std::string& filename) {
    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file %s", filename.c_str());
      exit(EXIT_FAILURE);
    }

    // Read first line
    auto read_row = [&](size_t i) {
      int j = 0;
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;

      while (std::getline(lineStream, cell, ',')) {
        add(stoi(cell), i, j);
        ++j;
      }

      return j;
    };

    // Read first row, and return row length
    int n = read_row(0);

    reserve(n);

    for (size_t i = 1; i < n; ++i) read_row(i);

    in_file.close();

    shrink_to_fit();
  }

 private:
  vector<FlowType> Ws;
  vector<PosType> Ps;
};

typedef GMeasureR2<int, double> MeasureR2;

MeasureR2 createRandom0N(size_t n, int seed = 13) {
  MeasureR2 mu;
  mu.reserve(n);

  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> Uniform01(0, 1);

  for (size_t i = 0; i < n; i++) mu.add(1, Uniform01(gen), Uniform01(gen));

  return mu;
}

// Solve separation problem as single core problem
int solveSeparation(const MeasureR2& Mu, const vector<double>& U,
                    const MeasureR2& Nu, const vector<double>& V, Vars& vars,
                    double vmin) {  // Avoid useless memory allocations
  int m = Mu.size();
  int n = Nu.size();

  assert(m == vars.size());

  int cmp = 0;

  for (int i = 0; i < m; ++i)
    if (U[i] > FEASIBILITY_TOL + vmin) {
      double best_v = -FEASIBILITY_TOL;
      double best_c = -1;
      int best_j = 0;

      for (int j = 0; j < n; ++j) {
        double violation = U[i] - V[j];
        if (violation > -best_v) {
          double c_ij = DISTANCE_R2(Mu.getP(i), Nu.getP(j));
          cmp++;
          violation = c_ij - violation;
          if (violation < best_v) {
            best_v = violation;
            best_c = c_ij;
            best_j = j;
            if (U[i] <= -best_v + vmin) break;
          }
        }
      }

      // Store most violated cuts for element i
      vars[i].b = m + best_j;
      vars[i].c = best_c;
    }

  // fprintf(stdout, "cmp: %d\n", cmp);

  return cmp;
}  // namespace R2

// Solve separation problem as multi core problem
void solveSeparationCore(const MeasureR2& Mu, const vector<double>& U,
                         const MeasureR2& Nu, const vector<double>& V,
                         Vars& vars,
                         double vmin) {  // Avoid useless memory allocations
  int m = Mu.size();
  int n = Nu.size();

  assert(m == vars.size());

#pragma omp parallel
  {
#pragma omp for schedule(guided)
    for (int i = 0; i < m; ++i)
      if (U[i] > FEASIBILITY_TOL + vmin) {
        double best_v = -FEASIBILITY_TOL;
        double best_c = -1;
        int best_j = 0;

        for (int j = 0; j < n; ++j) {
          double violation = U[i] - V[j];
          if (violation > -best_v) {
            double c_ij = DISTANCE_R2(Mu.getP(i), Nu.getP(j));
            violation = c_ij - violation;
            if (violation < best_v) {
              best_v = violation;
              best_c = c_ij;
              best_j = j;
              if (U[i] <= -best_v + vmin) break;
            }
          }
        }

        // Store most violated cuts for element i
        vars[i].b = m + best_j;
        vars[i].c = best_c;
      }
  }
}  // namespace R2

void solveSeparationGPU(concurrency::array_view<double, 2> xv,
                        vector<double>& U,
                        concurrency::array_view<double, 2> yv,
                        vector<double>& V, Vars& vars, int n, double vmin) {
  concurrency::array_view<double> Uv((int)U.size(), &U[0]);
  concurrency::array_view<double> Vv((int)V.size(), &V[0]);

  concurrency::array_view<Var> cv(vars.size(), vars);
  int m = U.size();

  concurrency::parallel_for_each(
      cv.extent, [=](concurrency::index<1> idx) restrict(amp) {
        if (Uv[idx[0]] > FEASIBILITY_TOL + vmin) {
          double best_v = -FEASIBILITY_TOL;
          double best_c = -1;
          int best_j = 0;

          for (int j = 0; j < n; ++j) {
            double violation = Uv[idx] - Vv[j];
            if (violation > -best_v) {
              double c_ij =
                  (xv(0, idx[0]) - yv(0, j)) * (xv(0, idx[0]) - yv(0, j)) +
                  (xv(1, idx[0]) - yv(1, j)) * (xv(1, idx[0]) - yv(1, j));
              violation = c_ij - violation;
              if (violation < best_v) {
                best_v = violation;
                best_c = c_ij;
                best_j = j;
              }
            }
          }

          // Store most violated cuts for element i
          cv[idx].b = m + best_j;
          cv[idx].c = best_c;
        }
      });

  try {
    cv.synchronize();
  } catch (const Concurrency::accelerator_view_removed& e) {
    fprintf(stdout, "solveSeparationGPU: %s\n", e.what());
  }
}

void solveSeparationGPUTile(concurrency::array_view<double, 2> xv,
                            vector<double>& U,
                            concurrency::array_view<double, 2> yv,
                            vector<double>& V, Vars& vars, int n, double vmin) {
  concurrency::array_view<double> Uv((int)U.size(), &U[0]);
  concurrency::array_view<double> Vv((int)V.size(), &V[0]);

  concurrency::array_view<Var> cv(vars.size(), vars);

  int m = U.size();

  static const int TS = 128;
  static const int TK = 2;

  concurrency::parallel_for_each(
      cv.extent.tile<TS>(),
      [=](concurrency::tiled_index<TS> t_idx) restrict(amp) {
        // Prepare shared tile
        int col = t_idx.local[0];
        // if (Uv[col] > FEASIBILITY_TOL + vmin) {
        int colGlobal = t_idx.global[0];
        tile_static double A[TK][TS];
        A[0][col] = xv(0, colGlobal);
        A[1][col] = xv(1, colGlobal);
        tile_static double Lu[TS];
        Lu[col] = Uv[colGlobal];

        // Local best cost
        int best_j = 0;
        double best_c = -1;
        double best_v = -FEASIBILITY_TOL;

        // Internal loop between pair of points
        for (int i = 0; i < n; i += TS) {
          tile_static double B[TK][TS];
          B[0][col] = yv(0, i + col);
          B[1][col] = yv(1, i + col);
          tile_static double Lv[TS];
          Lv[col] = Vv[i + col];

          t_idx.barrier.wait();

          for (int j = 0; j < TS; ++j) {
            double violation = Lu[col] - Lv[j];
            if (violation > -best_v) {
              double c_ij = (A[0][col] - B[0][j]) * (A[0][col] - B[0][j]) +
                            (A[1][col] - B[1][j]) * (A[1][col] - B[1][j]);
              // Lower precision, but faster speed using float instead of
              // double We do not use it for improving numerical stability
              /*concurrency::fast_math::pow(A[0][col] - B[0][j], 2) +
                       concurrency::fast_math::pow(A[1][col] - B[1][j],
                 2);*/
              violation = c_ij - violation;
              if (violation < best_v) {
                best_v = violation;
                best_c = c_ij;
                best_j = i + j;
              }
            }
          }

          t_idx.barrier.wait();
        }

        // Store most violated cuts for element i
        cv[colGlobal].b = m + best_j;
        cv[colGlobal].c = best_c;
        //}
      });

  try {
    cv.synchronize();
  } catch (const Concurrency::accelerator_view_removed& e) {
    fprintf(stdout, "solveSeparationGPUTile: %s\n", e.what());
  }
}

// Compute Kantorovich-Wasserstein distance between two measures
void DenseTransportationLP(const MeasureR2& Mu, const MeasureR2& Nu, int algo,
                           const std::string& msg) {
  // Timinig output
  auto start = std::chrono::high_resolution_clock::now();

  start = std::chrono::high_resolution_clock::now();

  int m = (int)Mu.size();
  int n = (int)Nu.size();

  typedef double CostType;
  typedef int64_t FlowType;

  // Build the graph for min cost flow
  DotSimplex<FlowType, CostType> simplex(n + m);

  // add first d source nodes
  for (int i = 0; i < m; ++i) simplex.addNode(i, +FlowType(Mu.getW(i)));

  for (int j = 0; j < n; ++j) simplex.addNode(m + j, -FlowType(Nu.getW(j)));

  simplex.resizeArcMemory(size_t(n * m));

#pragma omp parallel
  {
#pragma omp for schedule(static, m)
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        simplex.setArc(i * m + j, i, m + j,
                       DISTANCE_R2(Mu.getP(i), Nu.getP(j)));
  }

  //// Solve the problem to compute the distance
  DotSimplex<FlowType, CostType>::ProblemType status = simplex.run();

  switch (status) {
    case DotSimplex<>::INFEASIBLE:
      fprintf(stdout, "INFEASIBLE\n");
      break;
    case DotSimplex<>::OPTIMAL:
      break;
    case DotSimplex<>::UNBOUNDED:
      fprintf(stdout, "UNBOUNDED\n");
      break;
  }

  CostType fobj = simplex.totalCost();

  auto end = std::chrono::high_resolution_clock::now();
  double elapsed =
      double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()) /
      1000;

  fprintf(stdout, "%s %d %d Runtime %.6f Value %.6f status %d RAM %.2f\n",
          msg.c_str(), n, simplex.num_arcs(), elapsed, fobj, status,
          getUsedRAM());

  fflush(stdout);
}

// Compute Kantorovich-Wasserstein distance between two measures
void ColumnGeneration(const MeasureR2& Mu, const MeasureR2& Nu, int algo,
                      const std::string& msg) {
  int m = (int)Mu.size();
  int n = (int)Nu.size();

  // Timinig output
  auto start = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();
  double elapsed =
      double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()) /
      1000;

  start = std::chrono::high_resolution_clock::now();

  // Solve the problem
  typedef double CostType;
  typedef int64_t FlowType;

  // Build the graph for min cost flow
  DotSimplex<FlowType, CostType> simplex(n + m);

  // add first d source nodes
  for (int i = 0; i < m; ++i) simplex.addNode(i, +FlowType(Mu.getW(i)));

  for (int j = 0; j < n; ++j) simplex.addNode(m + j, -FlowType(Nu.getW(j)));

  int it = 0;
  int n_cuts = 0;
  CostType fobj = 0;
  double time_tot = 0;
  double sep_tot = 0;
  double mas_tot = 0;
  int cmp_tot = 0;
  vector<double> A(m, 0);
  vector<double> B(n, 0);

  DOT::Vars vars(m);
  for (int i = 0; i < m; ++i) vars[i].a = i;

  DOT::Vars vnew;
  vnew.reserve(2 * size_t(m + n) + 1);

  // Support for GPU
  const size_t K = 2;
  vector<double> XView(K * Mu.size());
  for (int i = 0; i < m; ++i) {
    auto p = Mu.getP(i);
    for (int k = 0; k < K; ++k) XView[i + k * m] = p[k];
  }

  vector<double> YView(K * Nu.size());
  for (int i = 0; i < n; ++i) {
    auto p = Nu.getP(i);
    for (int k = 0; k < K; ++k) YView[i + k * m] = p[k];
  }

  concurrency::array_view<double, 2> xv(K, m, XView);
  concurrency::array_view<double, 2> yv(K, n, YView);

  // Init the simplex
  DotSimplex<FlowType, CostType>::ProblemType status = simplex.run();

  // Start separation
  while (true) {
    start = std::chrono::high_resolution_clock::now();

    DotSimplex<FlowType, CostType>::ProblemType status = simplex.reRun();

    end = std::chrono::high_resolution_clock::now();
    elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start)
                         .count()) /
              1000;

    // Take the dual values
    for (int i = 0; i < m; ++i) A[i] = -simplex.potential(i);

    double umin = std::numeric_limits<double>::infinity();
    for (int j = 0; j < n; ++j) {
      B[j] = -simplex.potential(m + j);
      umin = std::min<double>(umin, B[j]);
    }

    mas_tot += elapsed;
    time_tot += elapsed;

    // Solve separation problem (with timing)
    auto sep_s = std::chrono::high_resolution_clock::now();

    if (algo == 0) cmp_tot += solveSeparation(Mu, A, Nu, B, vars, umin);
    if (algo == 1) solveSeparationCore(Mu, A, Nu, B, vars, umin);

    if (algo == 2) solveSeparationGPU(xv, A, yv, B, vars, n, umin);
    if (algo == 3) solveSeparationGPUTile(xv, A, yv, B, vars, n, umin);

    auto sep_e = std::chrono::high_resolution_clock::now();
    auto sep_elapsed =
        double(
            std::chrono::duration_cast<std::chrono::milliseconds>(sep_e - sep_s)
                .count()) /
        1000;
    sep_tot += sep_elapsed;
    time_tot += sep_elapsed;

    start = std::chrono::high_resolution_clock::now();

    vnew.clear();

    for (const auto& v : vars)
      if (v.c > -1) vnew.emplace_back(v);

    if (vnew.empty()) break;

    std::sort(vnew.begin(), vnew.end(),
              [](const auto& v, const auto& w) { return v.c > w.c; });

    // Replace old constraints with new ones
    int new_arcs = simplex.updateArcs(vnew);

    end = std::chrono::high_resolution_clock::now();
    elapsed += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                          end - start)
                          .count()) /
               1000;
    mas_tot += elapsed;
    time_tot += elapsed;
    n_cuts += new_arcs;

    // fprintf(stdout, "it %d: Time %.3f - Value: %.6f - NumRows: %d -
    // SepTime:
    // %.6f - GuTime: %.6f - Cuts: %d\n",
    //        it, time_tot, fobj, simplex.num_arcs(), sep_elapsed, elapsed,
    //        new_arcs);
    ++it;
  }

  fobj = simplex.totalCost();
  simplex.checkFeasibility();

  end = std::chrono::high_resolution_clock::now();
  elapsed =
      double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()) /
      1000;
  fprintf(stdout,
          "%s %d it %d FinalTime %.4f SepTime %.4f Master %.4f Value %.6f "
          "AddedVars %d NumVars %d CmpTot %d CmpRatio %.2f RAM %.2f\n",
          msg.c_str(), algo, it, time_tot, sep_tot, mas_tot, fobj, n_cuts,
          simplex.num_arcs(), cmp_tot, double(cmp_tot) / double(n * m),
          getUsedRAM());
  fflush(stdout);
}

}  // namespace R2
}  // namespace DOT