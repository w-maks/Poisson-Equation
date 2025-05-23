#include "solver.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

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

void PoissonSolver::run(int iterations, bool flag) {
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

void PoissonSolver::runParabolic(const int iterations) {
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

void PoissonSolver::runGradient(int iterations, const std::vector<double>& betas, double delta) {
    for (const double beta : betas) {
        std::ofstream sFile("S(iteration)[beta="+ std::to_string(beta)+ "]" + addon + ".csv");
        sFile << "iteration,S\n";
        const double S0 = S();
        sFile << 0 << ',' << S0 << '\n';

        for (int iter = 1; iter <= iterations; ++iter) {
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

