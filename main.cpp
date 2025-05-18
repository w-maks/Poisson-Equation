#include "solver.hpp"

int main() {
    constexpr int N = 31;
    constexpr double dx = 1.0;
    constexpr double d = 4.0;
    constexpr double x0 = 4.0;
    constexpr int itMax = 500;

    PoissonSolver solver(N, dx, d, x0);
    solver.run(itMax);
    return 0;
}
