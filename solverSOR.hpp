#ifndef SOLVERSOR_HPP
#define SOLVERSOR_HPP

#include "poissonSolver.hpp"

class SolverSOR : public PoissonSolver {
public:
    using PoissonSolver::PoissonSolver;
    void runSOR(int iterations, bool flag = false);
};

#endif

