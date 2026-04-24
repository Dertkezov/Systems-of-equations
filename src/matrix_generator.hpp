#include "csr_matrix.hpp"
#include <vector>
#include <cmath>

CSR_matrix generate_poisson_matrix(int grid_n) {
    size_t N = grid_n * grid_n;
    std::vector<double> v; v.reserve(N * 5);
    std::vector<size_t> c; c.reserve(N * 5);
    std::vector<size_t> r; r.resize(N + 1);
    r[0] = 0;
    size_t idx = 0;

    for (size_t i = 0; i < N; ++i) {
        int row = i / grid_n;
        int col = i % grid_n;

        v.push_back(4.0); c.push_back(i); idx++;
        if (col > 0) { v.push_back(-1.0); c.push_back(i - 1); idx++; }
        if (col < grid_n - 1) { v.push_back(-1.0); c.push_back(i + 1); idx++; }
        if (row > 0) { v.push_back(-1.0); c.push_back(i - grid_n); idx++; }
        if (row < grid_n - 1) { v.push_back(-1.0); c.push_back(i + grid_n); idx++; }

        r[i + 1] = idx;
    }
    return CSR_matrix(N, N, v, c, r);
}

full_matrix generate_poisson_matrix_PCG(int grid_n) {
    int N = grid_n * grid_n;
    full_matrix A(N, N);
    for (int i = 0; i < grid_n; i++) {
        for (int j = 0; j < grid_n; j++) {
            int idx = i * grid_n + j;
            A(idx, idx) = -4.0;
            if (i > 0) A(idx - grid_n, idx) = 1.0;
            if (i < grid_n - 1) A(idx + grid_n, idx) = 1.0;
            if (j > 0) A(idx - 1, idx) = 1.0;
            if (j < grid_n - 1) A(idx + 1, idx) = 1.0;
        }
    }
    return A;
}
