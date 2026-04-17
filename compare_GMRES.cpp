#include <fstream>
#include <vector>
#include <time.h>
#include <random>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include "Iteration_solver.hpp"
#include <algorithm>
#include <iostream>
#include "GMRES.hpp"
#include "matrix_generator.hpp"

int main() {
std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distrib_double(-1000.0, 1000.0);

    int grid_n = 30;
    CSR_matrix A_CSR = generate_poisson_matrix(grid_n);

    size_t R = grid_n * grid_n; 

    std::vector<double> b(R, 0.0);
    for (std::size_t k = 0; k < R; k++) {
        b[k] = distrib_double(gen);
    }

    int N = 10000;
    double e = 1e-6;
    std::vector<double> x0(R, 0.0);

    std::ofstream out_time_m("graphs/time_m.txt");
    std::ofstream out_iter_m("graphs/iter_m.txt");
    for (int m = 1; m < 10; m++){
        solution res_GM = GMRES_solver(A_CSR, b, x0, e, N, m);
        std::fill(x0.begin(), x0.end(), 0.0);
        double t = res_GM.total_time;
        int r = res_GM.p;
        out_time_m << m << " " << t << "\n";
        out_iter_m << m << " " << r << "\n";

    }
    out_time_m.close();
    out_iter_m.close();
}
