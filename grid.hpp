#ifndef GRID_HPP
#define GRID_HPP
#include <vector>

class Grid {
public:
    Grid(int nx, int ny, double init = 0.0);
    double& operator()(int i, int j);
    double  operator()(int i, int j) const;
    int nx() const {
        return nx_size;
    }
    int ny() const {
        return ny_size;
    }

private:
    int nx_size, ny_size;
    std::vector<double> values;
    int index(const int i, const int j) const {
        return i * ny_size + j;
    }
};
#endif
