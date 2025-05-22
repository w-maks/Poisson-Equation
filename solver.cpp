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

double PoissonSolver::Sloc(int i0, int j0) const {
    double S_value = 0.0;
    for(int i= i0-1; i <= i0+1; ++i) {
        for(int j= j0-1; j <= j0+1; ++j) {
            S_value += 0.5 * u(i,j) * laplacian(u, i, j) + rho(i,j) * u(i,j);
        }
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

void PoissonSolver::runParabolic(int iterations) {
    std::ofstream sFile("S(iteration)" + addon + ".csv");
    sFile << "iteration,S\n";

    double Sglob = S();
    sFile << 0 << ',' << Sglob << '\n';

    constexpr double testDelta[3] = {0.0, 0.5, 1.0};

    for(int iter = 1; iter <= iterations; ++iter) {
        for (int i = 1; i < 2 * N; ++i) {
            for(int j = 1; j < 2 * N; ++j) {

                Sglob -= Sloc(i,j);

                double Sval[3];
                const double temp = u(i,j);
                for(int k=0; k<3; ++k) {
                    u(i,j) = temp + testDelta[k];
                    Sval[k] = Sloc(i,j);
                }
                u(i,j) = temp;
                double denominator = Sval[0] - 2*Sval[1] + Sval[2];
                if(std::fabs(denominator) < 1e-12) {
                    denominator = 1e-12;
                }
                const double delta4 = 0.25 * (3*Sval[0] -4*Sval[1] +Sval[2]) / denominator;

                double deltas[4] = {0.0, 0.5, 1.0, delta4};
                double Sbest = std::numeric_limits<double>::infinity();
                double bestDelta = 0.0;

                for(double d : deltas) {
                    u(i,j) = temp + d;
                    double Slocal = Sloc(i, j);
                    if(Slocal < Sbest) {
                        Sbest = Slocal;
                        bestDelta = d;
                    }
                }
                u(i,j) = temp + bestDelta;
                Sglob += Sbest;
            }
        }
        sFile << iter << ',' << Sglob << '\n';
    }
}
