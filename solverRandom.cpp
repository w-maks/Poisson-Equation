#include "solverRandom.hpp"
#include <fstream>
#include <random>

void SolverRandom::runRandom(const int iterations, const std::vector<double>& rs) {
    for(const double r: rs) {
        u = Grid(2 * N + 1, 2 * N + 1, 0.0);
        std::ofstream sFile("S(iteration)[rand_r=" + std::to_string(r) + "]" + addon + ".csv");
        sFile << "iteration,S\n";
        sFile << 0 << ',' << S() << '\n';

        std::mt19937 generator(50);
        std::uniform_real_distribution<> urnd(-r, r);

        for(int iter = 1; iter <= iterations; ++iter) {
            const double S0 = S();
            for(int i = 1; i < 2 * N; ++i) {
                for(int j = 1; j < 2 * N; ++j) {
                    const double delta = urnd(generator);
                    const double S_better = S0 - Sloc(i,j, 0.0) + Sloc(i, j, delta);
                    if(S_better < S0) {
                        u(i,j) += delta;
                    }
                }
            }
            sFile << iter << ',' << S() << '\n';
        }
    }
}