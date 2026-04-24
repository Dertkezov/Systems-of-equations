#include <fstream>
#include <vector>
#include <time.h>
#include <random>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include "Iteration_solver.hpp"
#include <algorithm>
#include <iostream>
#include "matrix_generator.hpp"
#include "pre_CG.hpp"
#include "matrix_for_preconditioning.hpp"

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distrib_double(-1000.0, 1000.0);

    int grid_n = 30;
    full_matrix A = generate_poisson_matrix_PCG(grid_n);
    CSR_matrix A_CSR(A);

    size_t R = grid_n * grid_n; 

    std::vector<double> b(R, 0.0);
    for (std::size_t k = 0; k < R; k++) {
        b[k] = distrib_double(gen);
    }

    int N = 10000;
    double e = 1e-8;
    std::vector<double> x0(R, 0.0);

    solution res_cg = CG_solver(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);

    CSR_matrix L = Cholesky_factor(A * (-1.0));
    solution res_pcg = PCG_solver(A_CSR, b, x0, e, N, L);

    std::ofstream out_iter("graphs/iter_data.txt");
    size_t max_iter = std::max({res_cg.discrepancy_p.size(), 
                                res_pcg.discrepancy_p.size()});
    for (size_t k = 0; k < max_iter; k++) {
        double r_1 = (k < res_cg.discrepancy_p.size()) ? res_cg.discrepancy_p[k] : -1.0;
        double r_2 = (k < res_pcg.discrepancy_p.size()) ? res_pcg.discrepancy_p[k] : -1.0;
        out_iter << r_1 << " " << r_2 << "\n";
    }
    out_iter.close();

    std::ofstream out_cg("graphs/time_cg.txt");
    for (size_t k = 0; k < res_cg.time_p.size(); k++) {
        out_cg << res_cg.time_p[k] << " " << res_cg.discrepancy_p[k] << "\n";
    }
    out_cg.close();

    std::ofstream out_pcg("graphs/time_pcg.txt");
    for (size_t k = 0; k < res_pcg.time_p.size(); k++) {
        out_pcg << res_pcg.time_p[k] << " " << res_pcg.discrepancy_p[k] << "\n";
    }
    out_pcg.close();
}