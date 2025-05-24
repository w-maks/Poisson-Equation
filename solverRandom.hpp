#ifndef SOLVERRANDOM_HPP
#define SOLVERRANDOM_HPP

#include "poissonSolver.hpp"

class SolverRandom : public PoissonSolver {
public:
    using PoissonSolver::PoissonSolver;
    void runRandom(int iterations, const std::vector<double>& rs);
};
#endif

