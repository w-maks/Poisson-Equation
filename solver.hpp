#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "grid.hpp"
#include <string>

class PoissonSolver {
public:
    PoissonSolver(int N, double dx, double d, double x0, double omega = 1.0, std::string addon = "");
    void run(int iterations, bool flag=false);
    void runParabolic(int iterations);
    void runGradient(int iterations, const std::vector<double>& betas = {0.45, 0.46, 0.47, 48, 49}, double delta = 0.001);

private:
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
