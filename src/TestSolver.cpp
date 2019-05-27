/*
 * @fileoverview Copyright (c) 2019, , XXXX YYYYY,
 *
 */

#include "DOT_R2_Solver.h"
#include "DOT_Rk_Solver.h"

// Main entry point
int main(int argc, char* argv[]) {
  int n = 32 * 32;
  if (argc > 1) n = atoi(argv[1]);

  int algo = 0;
  if (argc > 2) algo = atoi(argv[2]);

  int seed = 13;

  // First basic test: random instances
  if (false) {
    fprintf(stdout, "First test random instances on unit square\n");
    int n = 64 * 64 * 4;

    auto Mu = DOT::R2::createRandom0N(n, n + 1);
    auto Nu = DOT::R2::createRandom0N(n, n + 13);

    DOT::R2::ColumnGeneration(Mu, Nu, 0, "CG CPU");
    DOT::R2::ColumnGeneration(Mu, Nu, 1, "CPU MultiCore");
    // DOT::R2::ColumnGeneration(Mu, Nu, 2, "GPU 1");
    // DOT::R2::ColumnGeneration(Mu, Nu, 3, "GPU tiles");

    // DOT::R2::DenseTransportationLP(Mu, Nu, 0, "Dense");
    fprintf(stdout, "-------------------------------------------\n\n");
  }

  // Basic test on R^300, based on word embedding vectors
  if (true) {
    fprintf(stdout, "Second test random instances on R^300\n");

    int n = 32 * 32 * 4;

    auto Mu = DOT::Rk::createRandom0N(n, n + 1);
    auto Nu = DOT::Rk::createRandom0N(n, n + 13);

    // The first one is far too slow
    // DOT::R2::ColumnGeneration(Mu, Nu, 0, "transport NS");
    DOT::Rk::ColumnGeneration(Mu, Nu, 1, "CPU MultiCore");
    DOT::Rk::ColumnGeneration(Mu, Nu, 2, "GPU 1");
    // DOT::Rk::ColumnGeneration(Mu, Nu, 3, "GPU Tiles");

    DOT::Rk::DenseTransportationLP(Mu, Nu, 0, "Dense");

    fprintf(stdout, "-------------------------------------------\n\n");
  }

  // Test with dotmark instances
  if (false) {
    // Provide path to your data
    std::string base = "..\\DOTA\\data\\DOTmark_1.0\\Data\\";

    if (argc > 3) base = argv[3];

    std::string SEP = "\\";

    std::vector<std::string> dirs = {
        "ClassicImages",  //"ClassicImages", "CauchyDensity", "GRFmoderate",
                          //"GRFrough", "GRFsmooth", "LogGRF", "LogitGRF",
                          //"MicroscopyImages", "Shapes",
    };

    std::vector<std::string> Fs = {
        "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
        "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"};

    std::vector<std::string> Ss = {"32", "64", "128", "256"};

    for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
        for (const auto& f11 : Fs) {
          for (const auto& f22 : Fs)
            if (f11 < f22) {
              std::string f1 = "data" + S + "_" + f11;
              std::string f2 = "data" + S + "_" + f22;

              DOT::R2::MeasureR2 Mu(base + dtype + SEP + f1);
              DOT::R2::MeasureR2 Nu(base + dtype + SEP + f2);

              for (int algo : vector<int>{3, 2, 1, 0})
                DOT::R2::ColumnGeneration(Mu, Nu, algo,
                                          "CG NS " + S + " " + f11 + " " + f2);
            }
        }
      }
    }
  }

  if (false) {
    // Provide path to word embedding data
    std::string base = "..\\DOTA\\data\\WordEmbed\\";
    if (argc > 3) base = argv[3];

    std::string SEP = "\\";

    std::vector<std::string> dirs = {
        "Random",  //"Books
    };

    std::vector<std::string> Fs = {
        "1.csv", "5.csv", "3.csv", "4.csv", "5.csv",
    };
    std::vector<std::string> Ss = {"64", "128", "256"};

    for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
        for (const auto& f11 : Fs) {
          for (const auto& f22 : Fs)
            if (f11 < f22) {
              std::string f1 = "rnd_" + S + "_" + f11;
              std::string f2 = "rnd_" + S + "_" + f22;

              DOT::Rk::MeasureRk Mu(base + dtype + SEP + f1);
              DOT::Rk::MeasureRk Nu(base + dtype + SEP + f2);

              auto start = std::chrono::high_resolution_clock::now();

              DOT::Rk::DenseTransportationLP(
                  Mu, Nu, 0, "Dense NS " + S + " " + f1 + " " + f2);
              auto end = std::chrono::high_resolution_clock::now();
              double elapsed =
                  double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - start)
                             .count()) /
                  1000;
              fprintf(stdout, "Total time: %.6f\n", elapsed);

              start = std::chrono::high_resolution_clock::now();
              DOT::Rk::ColumnGeneration(Mu, Nu, 0,
                                        "CG NS " + S + " " + f1 + " " + f2);
              end = std::chrono::high_resolution_clock::now();
              elapsed =
                  double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - start)
                             .count()) /
                  1000;
              fprintf(stdout, "Total time: %.6f\n", elapsed);

              start = std::chrono::high_resolution_clock::now();

              DOT::Rk::ColumnGeneration(Mu, Nu, 1,
                                        "CG NS " + S + " " + f1 + " " + f2);
              end = std::chrono::high_resolution_clock::now();
              elapsed =
                  double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - start)
                             .count()) /
                  1000;
              fprintf(stdout, "Total time: %.6f\n", elapsed);

              start = std::chrono::high_resolution_clock::now();

              DOT::Rk::ColumnGeneration(Mu, Nu, 2,
                                        "CG NS " + S + " " + f1 + " " + f2);
              end = std::chrono::high_resolution_clock::now();
              elapsed =
                  double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - start)
                             .count()) /
                  1000;
              fprintf(stdout, "Total time: %.6f\n", elapsed);
            }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}