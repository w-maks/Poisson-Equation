#include "solverSOR.hpp"
#include <fstream>

void SolverSOR::runSOR(int iterations, bool flag) {
    std::ofstream sFile("S(iteration)" + addon + ".csv");
    sFile << "iteration,S\n";

    for (int it = 1; it <= iterations; ++it) {
        for (int i = 1; i < 2 * N; ++i) {
            for (int j = 1; j < 2 * N; ++j) {
                const double temp = 0.25 * (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) + rho(i, j) * dx * dx);
                u(i,j) = (1-omega) * u(i,j) + omega * temp;
            }
        }

        sFile << it << ',' << S() << '\n';
        if(flag) {
            if (it == 100 || it == 500) {
                gridSaver(u,  "u(iter=" + std::to_string(it) + ").csv");

                Grid rhoPrim (2*N+1, 2*N+1, 0.0), delta(2*N+1, 2*N+1, 0.0);

                for (int i = 1; i < 2 * N; ++i) {
                    for (int j = 1; j < 2 * N; ++j) {
                        rhoPrim(i, j) = -laplacian(u, i, j);
                        delta(i, j) = rhoPrim(i, j) - rho(i, j);
                    }
                }
                gridSaver(rhoPrim, "rhoPrim(iter="  + std::to_string(it) + ").csv");
                gridSaver(delta, "delta(iter=" + std::to_string(it) + ").csv");
            }
        }
    }
}