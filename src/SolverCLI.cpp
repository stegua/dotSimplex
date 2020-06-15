/*
 * @fileoverview Copyright (c) 2019, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

// Main entry point
#ifdef CLI

#include "DOT_R2_Solver.h"
#include "DOT_Rk_Solver.h"

int main(int argc, char* argv[]) {
  int n = 32 * 32;
  if (argc > 1) n = atoi(argv[1]);

  int algo = 0;
  if (argc > 2) algo = atoi(argv[2]);

  int seed = 13;

  if (false) {
    std::string base =
        "C:\\Users\\Gualandi\\Google "
        "Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
    if (argc > 3) base = argv[3];

    std::string SEP = "\\";

    std::vector<std::string> dirs = {
        "ClassicImages",  //"ClassicImages", "CauchyDensity", "GRFmoderate",
                          //"GRFrough", "GRFsmooth", "LogGRF", "LogitGRF",
                          //"MicroscopyImages", "Shapes",
    };

    std::vector<std::string> Fs = {
        "1001.csv", "1002.csv",
        "1003.csv",  //"1004.csv", "1005.csv",
                     //"1006.csv", "1007.csv", "1008.csv",
                     //"1009.csv", "1010.csv"
    };

    std::vector<std::string> Ss = {"32", "64", "128", "256", "512"};

    for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
        for (const auto& f11 : Fs) {
          for (const auto& f22 : Fs)
            if (f11 < f22) {
              std::string f1 = "data" + S + "_" + f11;
              std::string f2 = "data" + S + "_" + f22;

              DOT::R2::MeasureR2 Mu(base + dtype + SEP + f1);
              DOT::R2::MeasureR2 Nu(base + dtype + SEP + f2);

              // for (int algo : vector<int> { 3, 2, 1, 0 })
              // DOT::R2::CuttingPlanes(Mu, Nu, algo, dtype+"_"+f1, f2);

              // DOT::R2::NS_TransportationLP0(Mu, Nu, 0, "transport NS");
              DOT::R2::ColumnGeneration(Mu, Nu, 3,
                                        "CG NS " + S + " " + f11 + " " + f2);

              // DOT::R2::NS_TransportationLP(Mu, Nu, 0, "transport NS");
              // DOT::R2::TransportationLP(Mu, Nu, algo);

              // transport NS 16384 268435456 Runtime 354.986000 Value
              // 1439261920 status 1

              // 4096 4096 1 - it 420 - Final Time 4.8810 - SepTime: 1.6780 -
              // GuTime: 3.2030 - Value: 319604347 - NumRows: 12741 4096 4096 1
              // - it 867 - Final Time 4.1110 - SepTime: 3.2260 - GuTime: 0.8850
              // - Value: 319604347 - NumRows: -3022088

              // 16384 16384 2 - it 1238 - Final Time 110.0220 -
              // SepTime: 29.3690 - GuTime: 80.6530 - Value: 1439 261 920 -
              // NumRows: 53114 16384 16384 3 - it 1868 - Final Time 29.4750 -
              // SepTime: 9.2440 - GuTime : 20.2310 - Value : 1439261920 -
              // NumRows : 32767
              // 42 755 285 499 170

              // 18 066 985 022 686
              // 65536 65536 2 - it 3940 - Final Time 4845.6740 - SepTime:
              // 837.1180 - MasterTime: 4008.5560 - Value: 7015138992 65536
              // 65536 2 - it 3940 - Final Time 3603.9250 - SepTime: 837.3540 -
              // MasterTime: 2766.5710 - Value: 7015138992 65536 65536 3 - it
              // 9317 - Final Time 1444.4920 - SepTime: 620.0080 - MasterTime:
              // 824.4840 - Value: 7015138992 65536 65536 3 - it 7560 - Final
              // Time 1347.3390 - SepTime: 507.7020 - MasterTime:  839.6370 -
              // Value: 7015138992 65536 65536 3 - it 151  - Final Time 152.0190
              // - SepTime:   9.5890 - Master :     142.4300 - Value: 7015138992
              //

              // 1024 1024 3 - it 47 - Final Time 0.0090 - SepTime: 0.0090 -
              // Master : 0.0000 - Value : 1.977237 - NumRows : 0
              // 2048 2048 3 - it 42 - Final Time 0.0350 - SepTime : 0.0010 -
              // Master : 0.0340 - Value : 1.382347 - NumRows : 0 4096 4096 3 -
              // it 66 - Final Time 0.3530 - SepTime : 0.0630 - Master : 0.2900
              // - Value : 1.653668 - NumRows : 0 8192 8192 3 - it 79 - Final
              // Time 1.1810 - SepTime : 0.1000 - Master : 1.0810 - Value
              // : 1.601901 - NumRows : 0 16384 16384 3 - it 105 - Final
              // Time 5.7720 - SepTime : 0.5250 - Master : 5.2470 - Value
              // : 2.062354 - NumRows : 0 32768 32768 3 - it 214 - Final
              // Time 33.3270 - SepTime : 3.0700 - Master : 30.2570 - Value
              // : 2.293760 - NumRows : 0 65536 65536 3 - it 524 - Final Time
              // 260.4430 - SepTime : 37.1210 - Master : 223.3220 - Value
              // : 2.876273 - NumRows : 0 131072 131072 3 - it 1587 - Final Time
              // 2025.9100 - SepTime : 385.9360 - Master : 1639.9740 - Value
              // : 3.518841 - NumRows : 0
            }
        }
      }
    }
  }

  if (false) {
    int n = 65536;

    auto Mu = DOT::R2::createRandom0N(n, n + 1);
    auto Nu = DOT::R2::createRandom0N(n, n + 13);

    // DOT::R2::DenseTransportationLP(Mu, Nu, 0, "transport NS");
    // DOT::R2::ColumnGeneration(Mu, Nu, 0, "transport NS");
    // DOT::R2::ColumnGeneration(Mu, Nu, 1, "transport NS");
    DOT::R2::ColumnGeneration(Mu, Nu, 2, "transport NS");
    DOT::R2::ColumnGeneration(Mu, Nu, 3, "transport NS");
  }

  // Multiple test
  if (true)
    for (int a : vector<int>{3, 1, 0})  //
      for (int N : vector<int>{
               // 1024, 1024 * 2, 1024 * 4, 1024 * 8, 1024 * 16, 1024 * 32,
               64 * 64  //,
                        // 1024 * 256
           })
        for (int I : vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {  //
          auto Mu = DOT::R2::createRandom0N(N, (I + 1) + N);
          auto Nu = DOT::R2::createRandom0N(N, (I + 13) + N);

          // DOT::R2::DenseTransportationLP(Mu, Nu, 0, "transport NS");
          DOT::R2::ColumnGeneration(
              Mu, Nu, a,
              "CG NS " + std::to_string(N) + " " + std::to_string(I));
        }

  if (false) {
    std::string base =
        "C:\\Users\\Gualandi\\Google "
        "Drive\\Ricerca\\DOTA\\data\\WordEmbed\\";
    if (argc > 3) base = argv[3];

    std::string SEP = "\\";

    std::vector<std::string> dirs = {
        "Random",  //"Books
    };

    std::vector<std::string> Fs = {
        "1.csv", "5.csv"  //, "3.csv", "4.csv", "5.csv",
    };
    std::vector<std::string> Ss = {
        "64"  //,"64", "128", "256"
    };

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

  if (false) {
    std::string base =
        "C:\\Users\\Gualandi\\Google "
        "Drive\\Ricerca\\DOTA\\data\\WordEmbed\\";
    if (argc > 3) base = argv[3];

    std::string SEP = "\\";

    std::vector<std::string> dirs = {
        "Random",  //"Books
    };

    std::vector<std::string> Fs = {
        "1.csv", "2.csv"  //, "3.csv"  //, "4.csv", "5.csv",
    };
    std::vector<std::string> Ss = {"64"};

    for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
        for (const auto& f11 : Fs) {
          for (const auto& f22 : Fs)
            if (f11 < f22) {
              std::string f1 = "rnd_" + S + "_" + f11;
              std::string f2 = "rnd_" + S + "_" + f22;

              DOT::Rk::MeasureRk Mu(base + dtype + SEP + f1);
              DOT::Rk::MeasureRk Nu(base + dtype + SEP + f2);

              // DOT::Rk::ColumnGeneration(Mu, Nu, 0,
              //                          "CG NS " + S + " " + f1 + " " +
              //                          f2);
              DOT::Rk::ColumnGeneration(Mu, Nu, 1,
                                        "CG NS " + S + " " + f1 + " " + f2);
              DOT::Rk::ColumnGeneration(Mu, Nu, 2,
                                        "CG NS " + S + " " + f1 + " " + f2);
              DOT::Rk::ColumnGeneration(Mu, Nu, 3,
                                        "CG NS " + S + " " + f1 + " " + f2);

              DOT::Rk::DenseTransportationLP(
                  Mu, Nu, 0, "Dense NS " + S + " " + f1 + " " + f2);
            }
        }
      }
    }
  }

  if (false) {
    int n = 32 * 32 * 4;

    auto Mu = DOT::Rk::createRandom0N(n, n + 1);
    auto Nu = DOT::Rk::createRandom0N(n, n + 13);

    // DOT::R2::ColumnGeneration(Mu, Nu, 0, "transport NS");
    DOT::Rk::ColumnGeneration(Mu, Nu, 1, "transport NS");
    DOT::Rk::ColumnGeneration(Mu, Nu, 2, "transport NS");
    DOT::Rk::ColumnGeneration(Mu, Nu, 3, "transport NS");
    DOT::Rk::DenseTransportationLP(Mu, Nu, 0, "transport NS");
  }
  return EXIT_SUCCESS;
}
#endif