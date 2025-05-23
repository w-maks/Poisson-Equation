#include "solver.hpp"

int main() {
    constexpr int N = 31;
    constexpr double dx = 1.0;
    constexpr double d = 4.0;
    constexpr double x0 = 4.0;
    constexpr int itMax = 500;

    //PoissonSolver solver1(N, dx, d, x0, 1, "zad1");
    //solver1.run(itMax, true);

    //PoissonSolver solver2(N, dx, d, x0, 1.9, "zad2");
    //solver2.run(itMax);

    //PoissonSolver solver3(N, dx, d, x0, 1, "zad3");
    //solver3.runParabolic(itMax);

    //PoissonSolver solver4(N, dx, d, x0, 1, "zad4");
    //solver4.runGradient(itMax);

    return 0;
}
