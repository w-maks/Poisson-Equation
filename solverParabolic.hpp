#ifndef SOLVERPARABOLIC_HPP
#define SOLVERPARABOLIC_HPP

#include "poissonSolver.hpp"

class SolverParabolic : public PoissonSolver {
public:
    using PoissonSolver::PoissonSolver;
    void runParabolic(int iterations);
};

#endif
