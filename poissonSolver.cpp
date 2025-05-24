#include "poissonSolver.hpp"
#include <cmath>
#include <fstream>

PoissonSolver::PoissonSolver(const int N, const double dx, const double d, const double x0, const double omega, std::string addon)
    : N(N), dx(dx), d(d), x0(x0), omega(omega), addon(addon),
    u(2 * N + 1, 2 * N + 1, 0.0), rho(2 * N + 1, 2 * N + 1, 0.0) {
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

double PoissonSolver::laplacian(const Grid& g, const int i, const int j) const {
    return (g(i+1,j)+g(i-1,j)+g(i,j+1)+g(i,j-1)-4*g(i,j)) / (dx*dx);
}

double PoissonSolver::S() const {
    double S_value = 0.0;
    for (int i = 1; i < 2*N; ++i) {
        for (int j = 1; j < 2*N; ++j) {
            S_value += 0.5 * u(i,j) * laplacian(u,i,j) + u(i,j) * rho(i,j);
        }
    }
    return -S_value * dx * dx;
}

double PoissonSolver::Sloc(int i0, int j0, const double delta) const
{
    auto u_loc = [&](int i, int j) -> double {
        return (i == i0 && j == j0) ? u(i, j) + delta : u(i, j);
    };

    auto lap = [&](int i, int j) -> double {
        return (u_loc(i+1,j) + u_loc(i-1,j) + u_loc(i,j+1) + u_loc(i,j-1) - 4*u_loc(i,j)) / (dx * dx);
    };

    double S_value = 0.0;
    for (int i = i0-1; i <= i0+1; ++i)
        for (int j = j0-1; j <= j0+1; ++j) {
            S_value += 0.5 * u_loc(i, j) * lap(i, j) + rho(i, j) * u_loc(i, j);
        }
    return -S_value * dx * dx;
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


