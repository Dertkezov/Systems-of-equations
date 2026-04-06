#include <fstream>
#include <vector>
#include <time.h>
#include <random>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include "Iteration_solver.hpp"
#include <algorithm>
#include <iostream>
#include "Chebyshev_speed_symmetry.hpp"
#include "Steps_for_chebyshev.hpp"
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
    double omega = 1.5; 
    std::vector<double> x0(R, 0.0);

    
    solution res_gs = Gauss_solver(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);


    solution res_SOR = SOR_solver(A_CSR, b, x0, e, omega, N);
    std::fill(x0.begin(), x0.end(), 0.0);


    std::vector<double> D(R);
    for (size_t i = 0; i < R; i++){
        D[i] = A_CSR(i, i);
    }
    auto gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Jacobi_step(A_CSR, D, b, x);
    };
    double rho = estimate_rho(gauss_step, x0);
    solution res_gauss_accel = Chebyshev_accelerate(gauss_step, A_CSR, b, x0, rho, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);


    solution res_descent = fast_descent(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);


    solution res_CG = CG_solver(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);


    std::ofstream out_iter("graphs/iter_data.txt");
    size_t max_iter = std::max({res_gs.discrepancy_p.size(), 
                                res_SOR.discrepancy_p.size(), 
                                res_gauss_accel.discrepancy_p.size(),
                                res_descent.discrepancy_p.size(),
                                res_CG.discrepancy_p.size()});
    for (size_t k = 0; k < max_iter; k++) {
        double r_1 = (k < res_gs.discrepancy_p.size()) ? res_gs.discrepancy_p[k] : -1.0;
        double r_2 = (k < res_SOR.discrepancy_p.size()) ? res_SOR.discrepancy_p[k] : -1.0;
        double r_3  = (k < res_gauss_accel.discrepancy_p.size())  ? res_gauss_accel.discrepancy_p[k]  : -1.0;
        double r_4  = (k < res_descent.discrepancy_p.size())  ? res_descent.discrepancy_p[k]  : -1.0;
        double r_5  = (k < res_CG.discrepancy_p.size())  ? res_CG.discrepancy_p[k]  : -1.0;
        out_iter << r_1 << " " << r_2 << " " << r_3 << " " << r_4 << " " << r_5 << "\n";
    }
    out_iter.close();

    std::ofstream out_gs("graphs/time_gauss.txt");
    for (size_t k = 0; k < res_gs.time_p.size(); k++) {
        out_gs << res_gs.time_p[k] << " " << res_gs.discrepancy_p[k] << "\n";
    }
    out_gs.close();

    std::ofstream out_sor("graphs/time_sor.txt");
    for (size_t k = 0; k < res_SOR.time_p.size(); k++) {
        out_sor << res_SOR.time_p[k] << " " << res_SOR.discrepancy_p[k] << "\n";
    }
    out_sor.close();

    std::ofstream out_cheb("graphs/time_chebyshev.txt");
    for (size_t k = 0; k < res_gauss_accel.time_p.size(); k++) {
        out_cheb << res_gauss_accel.time_p[k] << " " << res_gauss_accel.discrepancy_p[k] << "\n";
    }
    out_cheb.close();

    std::ofstream out_slope("graphs/time_slope.txt");
    for (size_t k = 0; k < res_descent.time_p.size(); k++) {
        out_slope << res_descent.time_p[k] << " " << res_descent.discrepancy_p[k] << "\n";
    }
    out_slope.close();

    std::ofstream out_cg("graphs/time_cg.txt");
    for (size_t k = 0; k < res_CG.time_p.size(); k++) {
        out_cg << res_CG.time_p[k] << " " << res_CG.discrepancy_p[k] << "\n";
    }
    out_cg.close();
}
