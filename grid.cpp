#include "grid.hpp"

Grid::Grid(const int nx, const int ny, const double init)
    : nx_size(nx), ny_size(ny), values(nx * ny, init) {}

double& Grid::operator()(const int i, const int j) {
    return values[index(i, j)];
}
double  Grid::operator()(const int i, const int j) const {
    return values[index(i, j)];
}
