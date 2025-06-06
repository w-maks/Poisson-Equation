#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "grid.hpp"
#include <string>

class PoissonSolver {
public:
    PoissonSolver(int N, double dx, double d, double x0, double omega = 1.0, std::string addon = "");

protected:
    const int N;
    const double dx, d, x0;
    const double omega;
    std::string addon;
    Grid u, rho;

    void fillRho();
    double laplacian(const Grid& g, int i, int j) const;
    double S() const;
    double Sloc(int i, int j, double delta) const;
    void gridSaver(const Grid& g, const std::string& name) const;
};
#endif
