#include "solver.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>

PoissonSolver::PoissonSolver(const int N, const double dx, const double d, const double x0)
    : N(N), dx(dx), d(d), x0(x0), u(2 * N + 1, 2 * N + 1, 0.0), rho(2 * N + 1, 2 * N + 1, 0.0) {
    fillRho();
}

void PoissonSolver::fillRho() {
    for (int i = -N; i <= N; ++i) {
        for (int j = -N; j <= N; ++j) {
            const double x = i * dx;
            const double y = j * dx;
            rho(i + N, j + N) = std::exp(-((x - x0)*(x - x0) + y*y)/(d*d)) - std::exp(-((x + x0)*(x + x0) + y*y)/(d*d));
        }
    }
}

double PoissonSolver::calculateRhoPrim(const Grid& g, const int i, const int j) const {
    return -(g(i+1,j)+g(i-1,j)+g(i,j+1)+g(i,j-1)-4*g(i,j)) / (dx*dx);
}

double PoissonSolver::S() const {
    double S_value = 0.0;
    for (int i = 1; i < 2*N; ++i) {
        for (int j = 1; j < 2*N; ++j) {
            const double term = 0.5 * u(i,j) * (calculateRhoPrim(u,i,j) - rho(i,j));
            S_value += term;
        }
    }
    return S_value * dx * dx;
}

void PoissonSolver::gridSaver(const Grid& g, const std::string& name) const {
    std::ofstream out(name);
    out << "i,j,value\n";
    for (int i = 0; i < g.nx(); ++i) {
        for (int j = 0; j < g.ny(); ++j) {
            out << (i - N) << ',' << (j - N) << ',' << g(i, j) << '\n';
        }
    }
}

void PoissonSolver::run(int iterations) {
    std::ofstream sFile("S(iteration).csv");
    sFile << "iteration,S\n";

    for (int it = 1; it <= iterations; ++it) {
        for (int i = 1; i < 2 * N; ++i) {
            for (int j = 1; j < 2 * N; ++j) {
                u(i, j) = 0.25 * (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) + rho(i, j) * dx * dx);
            }
        }

        sFile << it << ',' << S() << '\n';

        if (it == 100 || it == 500) {
            gridSaver(u,  "u(iter=" + std::to_string(it) + ").csv");

            Grid rhoPrim (2*N+1, 2*N+1, 0.0), delta(2*N+1, 2*N+1, 0.0);

            for (int i = 1; i < 2 * N; ++i) {
                for (int j = 1; j < 2 * N; ++j) {
                    rhoPrim(i, j) = calculateRhoPrim(u, i, j);
                    delta(i, j) = rhoPrim(i, j) - rho(i, j);
                }
            }
            gridSaver(rhoPrim, "rhoPrim(iter="  + std::to_string(it) + ").csv");
            gridSaver(delta, "delta(iter=" + std::to_string(it) + ").csv");
        }
    }
}
