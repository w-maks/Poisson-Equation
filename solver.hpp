#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "grid.hpp"
#include <string>

class PoissonSolver {
public:
    PoissonSolver(int N, double dx, double d, double x0);
    void run(int iterations);

private:
    const int N;
    const double dx, d, x0;
    Grid u, rho;

    void fillRho();
    double calculateRhoPrim(const Grid& g, int i, int j) const;
    double S() const;
    void gridSaver(const Grid& g, const std::string& name) const;
};
#endif
