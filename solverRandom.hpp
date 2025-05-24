#ifndef SOLVERRANDOM_HPP
#define SOLVERRANDOM_HPP

#include "poissonSolver.hpp"

class SolverRandom : public PoissonSolver {
public:
    using PoissonSolver::PoissonSolver;
    void runRandom(int iterations, double r = 0.3);
};
#endif

