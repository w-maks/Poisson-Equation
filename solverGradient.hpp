#ifndef SOLVERGRADIENT_HPP
#define SOLVERGRADIENT_HPP

#include "poissonSolver.hpp"

class SolverGradient : public PoissonSolver {
public:
    using PoissonSolver::PoissonSolver;
    void runGradient(int iterations, const std::vector<double>& betas, double delta = 0.001);
};
#endif
