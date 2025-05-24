#include "solverGradient.hpp"
#include <fstream>

void SolverGradient::runGradient(const int iterations, const std::vector<double>& betas, const double delta) {
    for (const double beta : betas) {
        u = Grid(2 * N + 1, 2 * N + 1, 0.0);
        std::ofstream sFile("S(iteration)[beta="+ std::to_string(beta)+ "]" + addon + ".csv");
        sFile << "iteration,S\n";
        sFile << 0 << ',' << S() << '\n';

        for (int iter = 1; iter <= iterations; ++iter) {
            const double S0 = S();
            for (int i = 1; i < 2 * N; ++i) {
                for (int j = 1; j < 2 * N; ++j) {
                    const double Sloc0 = Sloc(i,j, 0.0);
                    const double Splus  = S0 - Sloc0 + Sloc(i, j,  delta);
                    const double Sminus = S0 - Sloc0 + Sloc(i, j, -delta);
                    const double grad  = (Splus - Sminus) / (2 * delta);
                    u(i, j) -= beta * grad;
                }
            }
            sFile << iter << ',' << S() << '\n';
        }
    }
}