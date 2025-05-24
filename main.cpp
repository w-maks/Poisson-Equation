#include "poissonSolver.hpp"
#include "solverSOR.hpp"
#include "solverParabolic.hpp"
#include "solverGradient.hpp"
#include "solverRandom.hpp"

int main() {
    constexpr int N = 31;
    constexpr double dx = 1.0;
    constexpr double d = 4.0;
    constexpr double x0 = 4.0;
    constexpr int itMax = 500;

    //SolverSOR solver1(N, dx, d, x0, 1, "zad1");
    //solver1.runSOR(itMax, true);

    //SolverSOR solver2(N, dx, d, x0, 1.9, "zad2");
    //solver2.runSOR(itMax);

    //SolverParabolic solver3(N, dx, d, x0, 1, "zad3");
    //solver3.runParabolic(itMax);

    //SolverGradient solver4(N, dx, d, x0, 1, "zad4");
    //solver4.runGradient(itMax, {0.45, 0.46, 0.47, 0.48, 0.49, 0.495, 0.497, 0.5});

    SolverRandom solver5(N, dx, d, x0, 1, "zad5");
    solver5.runRandom(itMax, {0.1, 0.2, 0.3, 0.5, 1});

    return 0;
}
