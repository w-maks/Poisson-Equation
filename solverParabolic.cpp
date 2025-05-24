#include "solverParabolic.hpp"
#include <fstream>
#include <cmath>

void SolverParabolic::runParabolic(const int iterations) {
    std::ofstream sFile("S(iteration)" + addon + ".csv");
    sFile << "iteration,S\n";
    const double S0 = S();
    sFile << 0 << ',' << S0 << '\n';

    for (int iter = 1; iter <= iterations; ++iter) {
        for (int i = 1; i < 2 * N; ++i) {
            for (int j = 1; j < 2 * N; ++j) {
                const double Sloc0 = Sloc(i,j, 0.0);
                const double S1 = S0;
                const double S2 = S0 - Sloc0 + Sloc(i, j, 0.5);
                const double S3 = S0 - Sloc0 + Sloc(i, j, 1);

                const double denominator = S1 - 2 * S2 + S3;;
                const double delta4 = (std::fabs(denominator) > 1e-12) ? 0.25 * (3 * S1 - 4 * S2 + S3) / denominator : 0.0;

                const double S4 = S0 - Sloc0 + Sloc(i, j, delta4);

                const double deltas[4] = {0, 0.5, 1, delta4};
                const double S[4] = {S1, S2, S3, S4};

                double bestDelta = 0.0, bestS = S1;

                for (int k = 1; k < 4; ++k) {
                    if (S[k] < bestS) {
                        bestS = S[k]; bestDelta = deltas[k];
                    }
                }
                u(i,j) += bestDelta;
            }
        }
        sFile << iter << ',' << S() << '\n';
    }
}