/**
 * @fileoverview Copyright (c) 2019,, XXXX YYYYY,
 *
 */

#pragma once

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

#ifdef _WIN32
#include <amp.h>
#include <amp_math.h>
#endif

#include <emmintrin.h>
#include <immintrin.h>

inline double hsum_double_avx(__m256d v) {
  __m128d vlow = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);  // high 128
  vlow = _mm_add_pd(vlow, vhigh);               // reduce down to 128

  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

inline double s3(const double* a, const double* b) {
  __m256d xx1 = _mm256_load_pd(a);
  __m256d yy1 = _mm256_load_pd(b);
  __m256d aa = _mm256_sub_pd(xx1, yy1);
  aa = _mm256_mul_pd(aa, aa);  // power
  return hsum_double_avx(aa);
}

inline double s4(const double* a, const double* b) {
  __m512d aa = _mm512_sub_pd(_mm512_load_pd(a), _mm512_load_pd(b));
  aa = _mm512_mul_pd(aa, aa);  // power
  return _mm512_reduce_add_pd(aa);
}

namespace DOT {
namespace Rk {

typedef double Cost;
typedef int Flow;

const int K = 300;

inline Cost DISTANCE_Rk_AVX(const Cost* x, const Cost* y) {
  Cost r = 0;
  for (size_t i = 0; i < K - 4; i = i + 8) r += s4(&x[i], &y[i]);

  r += s3(&x[K - 4], &y[K - 4]);
  return r;
}

inline Cost DISTANCE_Rk(const Cost* x, const Cost* y) {
  Cost r = 0;
  Cost a = 0;

  for (size_t i = 0; i < K; ++i) {
    a = x[i] - y[i];
    a = a * a;
    r += a;
  }
  return r;
}

// Container for general discrete measure
template <typename FlowType = int, typename PosType = double>
class GMeasureRk {
 public:
  GMeasureRk() {}

  GMeasureRk(const std::string& filename) { readFromFile(filename); }

  // setters
  void reserve(size_t n) {
    Ws.reserve(n);
    Ps.reserve(K * n);
  }

  void addWeight(FlowType _w) { Ws.emplace_back(_w); }

  void addSupport(PosType _p) { Ps.emplace_back(_p); }

  // Use as few memory as possible
  void shrink_to_fit() {
    Ws.shrink_to_fit();
    Ps.shrink_to_fit();
  }

  // getters
  size_t size() const { return Ws.size(); }

  FlowType getW(size_t i) const { return Ws[i]; }

  const PosType* getP(size_t i) const { return &Ps[K * i]; }

  PosType* getP(size_t i) { return &Ps[K * i]; }

  vector<double> getPs() const { return Ps; }

  // Parse from file
  void readFromFile(const std::string& filename) {
    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file %s", filename.c_str());
      exit(EXIT_FAILURE);
    }

    // Read first line
    while (!in_file.eof()) {
      size_t j = 0;
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;

      while (std::getline(lineStream, cell, ' ')) {
        if (j == 1) addWeight(1);  // int64_t(stoll(cell)));
        if (j > 1) addSupport(stod(cell));
        ++j;
      }
    };

    in_file.close();

    // Use as few memory as possible
    shrink_to_fit();
  }

 private:
  vector<FlowType> Ws;
  vector<PosType> Ps;
};  // namespace Rk

typedef GMeasureRk<int, double> MeasureRk;

MeasureRk createRandom0N(size_t n, int seed = 13) {
  MeasureRk mu;
  mu.reserve(n);

  std::random_device rd;
  std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> Uniform01(0, 1);

  for (size_t i = 0; i < n; i++) {
    mu.addWeight(1);
    for (int i = 0; i < K; i++) mu.addSupport(Uniform01(gen));
  }

  return mu;
}

// Solve separation problem as single core problem
int solveSeparation(const MeasureRk& Mu, const vector<Cost>& U,
                    const MeasureRk& Nu, const vector<Cost>& V,
                    DOT::FVars& vars,
                    Cost vmin) {  // Avoid useless memory allocations
  int m = Mu.size();
  int n = Nu.size();

  assert(m == vars.size());

  int cmp = 0;

  for (int i = 0; i < m; ++i)
    if (U[i] > vmin + FEASIBILITY_TOL) {
      Cost best_v = -PRIC_TOL;
      Cost best_c = -1;
      int best_j = 0;

      for (int j = 0; j < n; ++j) {
        Cost violation = U[i] - V[j];
        if (violation > -best_v) {
          Cost c_ij = DISTANCE_Rk(Mu.getP(i), Nu.getP(j));
          cmp++;
          violation = c_ij - violation;
          if (violation < best_v) {
            best_v = violation;
            best_c = c_ij;
            best_j = j;
            if (U[i] <= vmin - best_v) break;
          }
        }
      }

      // Store most violated cuts for element i
      vars[i].b = m + best_j;
      vars[i].c = best_c;
    }

  // fprintf(stdout, "cmp: %d\n", cmp);

  return cmp;
}

// Solve separation problem as multi core problem
void solveSeparationCore(const MeasureRk& Mu, const vector<Cost>& U,
                         const MeasureRk& Nu, const vector<Cost>& V,
                         DOT::FVars& vars,
                         Cost vmin) {  // Avoid useless memory allocations
  int m = Mu.size();
  int n = Nu.size();

  assert(m == vars.size());

#pragma omp parallel
  {
#pragma omp for schedule(guided)
    for (int i = 0; i < m; ++i) {
      // if (U[i] > vmin + FEASIBILITY_TOL) {
      Cost best_v = -PRIC_TOL;
      Cost best_c = -1;
      int best_j = -1;

      for (int j = 0; j < n; ++j) {
        Cost violation = U[i] - V[j];
        if (violation > -best_v) {
          Cost c_ij = DISTANCE_Rk(Mu.getP(i), Nu.getP(j));
          violation = c_ij - violation;
          if (violation < best_v) {
            best_v = violation;
            best_c = c_ij;
            best_j = j;
            // if (U[i] <= vmin - best_v) break;
          }
        }
      }

      // Store most violated cuts for element i
      vars[i].b = m + best_j;
      vars[i].c = best_c;
    }
    //}
  }
}

void solveSeparationGPU(const concurrency::array_view<Cost, 2> xv,
                        vector<Cost>& U,
                        const concurrency::array_view<Cost, 2> yv,
                        vector<Cost>& V, DOT::FVars& vars, int n, Cost vmin) {
  //  Avoid useless memory allocations
  concurrency::array_view<Cost> Uv((int)U.size(), &U[0]);
  concurrency::array_view<Cost> Vv((int)V.size(), &V[0]);

  concurrency::array_view<DOT::FVar> cv(vars.size(), vars);
  int m = U.size();

  concurrency::parallel_for_each(
      cv.extent, [=](concurrency::index<1> idx) restrict(amp) {
        // It is unclear whether it worths or it slows down the exectution
        // if (Uv[idx] > vmin + FEASIBILITY_TOL) {
        Cost best_v = -PRIC_TOL;
        Cost best_c = -1;
        int best_j = 0;

        for (int j = 0; j < n; ++j) {
          Cost violation = Uv[idx] - Vv[j];
          if (violation > -best_v) {
            Cost c_ij = 0;
            for (int k = 0; k < K; ++k)
              c_ij += (xv(k, idx[0]) - yv(k, j)) * (xv(k, idx[0]) - yv(k, j));

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
        //}
      });

  try {
    cv.synchronize();
  } catch (const Concurrency::accelerator_view_removed& e) {
    fprintf(stdout, "solveSeparationGPU: %s\n", e.what());
  }
}

// void solveSeparationGPUTile(concurrency::array_view<Cost, 2> xv,
//                            vector<Cost>& U,
//                            concurrency::array_view<Cost, 2> yv,
//                            vector<Cost>& V, DOT::FVars& vars, int n) {
//  // Avoid useless memory allocations
//  concurrency::array_view<Cost> Uv((int)U.size(), &U[0]);
//  concurrency::array_view<Cost> Vv((int)V.size(), &V[0]);
//
//  concurrency::array_view<DOT::FVar> cv(vars.size(), vars);
//
//  int m = U.size();
//
//  static const int TS = 8;
//  static const int TK = K;
//
//  concurrency::parallel_for_each(
//      cv.extent.tile<TS>(),
//      [=](concurrency::tiled_index<TS> t_idx) restrict(amp) {
//        // Prepare shared tile
//        int col = t_idx.local[0];
//        int colGlobal = t_idx.global[0];
//        tile_static Cost A[TK][TS];
//        for (int k = 0; k < K; ++k) A[k][col] = xv(k, colGlobal);
//        tile_static Cost Lu[TS];
//        Lu[col] = Uv[colGlobal];
//        t_idx.barrier.wait();
//
//        // Local best cost
//        int best_j = 0;
//        Cost best_c = -1;
//        Cost best_v = -PRIC_TOL;
//
//        // Internal loop between pair of points
//        for (int i = 0; i < n; i += TS) {
//          // Does it pay off to using tile on R^300 ?
//
//          // tile_static Cost B[TK][TS];
//          // for (int k = 0; k < K; ++k) B[k][col] = yv(k, i + col);
//          // tile_static Cost Lv[TS];
//          // Lv[col] = Vv[i + col];
//
//          // t_idx.barrier.wait();
//
//          for (int j = 0; j < TS; ++j) {
//            Cost violation = -Lu[col] + Vv[i + j];
//            // Lv[j];
//            if (violation < best_v) {
//              Cost c_ij = 0;
//              for (int k = 0; k < K; ++k)
//                c_ij +=
//                    concurrency::fast_math::pow(A[k][col] - yv(k, i + j), 2);
//              //(A[k][col] - yv(k, i + j)) * (A[k][col] - yv(k, i + j));
//              violation = c_ij + violation;
//              if (violation < best_v) {
//                best_v = violation;
//                best_c = c_ij;
//                best_j = i + j;
//              }
//            }
//          }
//          // t_idx.barrier.wait();
//        }
//
//        // Store most violated cuts for element i
//        if (best_v < -PRIC_TOL) {
//          cv[colGlobal].b = m + best_j;
//          cv[colGlobal].c = best_c;
//        }
//      });
//
//  try {
//    cv.synchronize();
//  } catch (const Concurrency::accelerator_view_removed& e) {
//    fprintf(stdout, "solveSeparationGPUTile: %s\n", e.what());
//  }
//}

// Compute Kantorovich-Wasserstein distance between two measures
void DenseTransportationLP(const MeasureRk& Mu, const MeasureRk& Nu, int algo,
                           const std::string& msg) {
  // Timinig output
  auto start = std::chrono::high_resolution_clock::now();

  int m = (int)Mu.size();
  int n = (int)Nu.size();

  typedef int64_t FlowType;
  typedef double CostType;

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
                       DISTANCE_Rk(Mu.getP(i), Nu.getP(j)));
  }
  // Solve the problem to compute the distance
  DotSimplex<FlowType, CostType>::ProblemType status = simplex.run();

  switch (status) {
    case DotSimplex<>::INFEASIBLE:
      fprintf(stdout, "INFEASIBLE\n");
      break;
    case DotSimplex<>::OPTIMAL:
      // fprintf(stdout, "OPTIMAL\n");
      break;
    case DotSimplex<>::UNBOUNDED:
      fprintf(stdout, "UNBOUNDED\n");
      break;
  }

  Cost fobj = simplex.totalCost();

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
void ColumnGeneration(const MeasureRk& Mu, const MeasureRk& Nu, int algo,
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

  // Solve the problem
  typedef int64_t FlowType;
  typedef double CostType;

  // Build the graph for min cost flow
  DotSimplex<FlowType, CostType> simplex(n + m);

  // add first d source nodes
  for (int i = 0; i < m; ++i) simplex.addNode(i, +FlowType(Mu.getW(i)));

  for (int j = 0; j < n; ++j) simplex.addNode(m + j, -FlowType(Nu.getW(j)));

  int it = 0;
  int n_cuts = 0;
  Cost fobj = 0;
  double time_tot = 0;
  double sep_tot = 0;
  double mas_tot = 0;
  int cmp_tot = 0;
  vector<Cost> A(m, 0);
  vector<Cost> B(n, 0);

  // Support for GPU
  vector<Cost> XView(K * Mu.size());
  for (int i = 0; i < m; ++i) {
    auto p = Mu.getP(i);
    for (int k = 0; k < K; ++k) XView[i + k * m] = p[k];
  }

  vector<Cost> YView(K * Nu.size());
  for (int i = 0; i < n; ++i) {
    auto p = Nu.getP(i);
    for (int k = 0; k < K; ++k) YView[i + k * m] = p[k];
  }

  const concurrency::array_view<Cost, 2> xv(K, m, &XView[0]);
  const concurrency::array_view<Cost, 2> yv(K, n, &YView[0]);

  DotSimplex<FlowType, CostType>::ProblemType status = simplex.run();

  DOT::FVars vars(m);
  for (int i = 0; i < m; ++i) vars[i].a = i;

  DOT::FVars vnew;
  vnew.reserve(2 * size_t(m + n) + 1);

  // Start separation
  while (true) {
    start = std::chrono::high_resolution_clock::now();

    DotSimplex<FlowType, CostType>::ProblemType status = simplex.reRun();

    fobj = simplex.totalCost();

    // Take the dual values
    for (int i = 0; i < m; ++i) A[i] = -simplex.potential(i);

    Cost umin = std::numeric_limits<Cost>::infinity();
    for (int j = 0; j < n; ++j) {
      B[j] = -simplex.potential(m + j);
      umin = std::min<Cost>(umin, B[j]);
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start)
                         .count()) /
              1000;
    mas_tot += elapsed;
    time_tot += elapsed;

    // Solve separation problem (with timing)
    auto sep_s = std::chrono::high_resolution_clock::now();

    if (algo == 0) cmp_tot += solveSeparation(Mu, A, Nu, B, vars, umin);
    if (algo == 1) solveSeparationCore(Mu, A, Nu, B, vars, umin);

    if (algo == 2) solveSeparationGPU(xv, A, yv, B, vars, n, umin);
    // if (algo == 3) solveSeparationGPUTile(xv, A, yv, B, vars, n);

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

    for (auto& v : vars) {
      if (v.c > -0.5) vnew.emplace_back(v);
      v.c = -1;
    }

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

    // fprintf(stdout,
    //        "it %d: Time %.3f Value %.6f NumRows %d SepTime %.6f GuTime %.6f "
    //        "Cuts %d\n ",
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

}  // namespace Rk
}  // namespace DOT