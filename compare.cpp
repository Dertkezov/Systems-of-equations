#include <fstream>
#include <vector>
#include <time.h>
#include <random>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include "Iteration_solver.hpp"
#include "MPI_speed.hpp"
#include <algorithm>
#include <iostream>

int main() {
    std::size_t R = 1000;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distrib_double(-1000.0, 1000.0);
    std::uniform_int_distribution<int> distrib_int(0, R - 1);

    full_matrix A(R, R);

    int i, j;
    for (int k = 0; k <= 10000; k++) {
        i = distrib_int(gen);
        j = distrib_int(gen);
        while (A(i, j) != 0) {
            i = distrib_int(gen);
            j = distrib_int(gen);
        } 
        A(i, j) = distrib_double(gen);
    }

    for (size_t k = 0; k < R; k++) {
        A(k, k) = 5000.0; 
    }

    CSR_matrix A_CSR(A);

    std::vector<double> b(R, 0.0);
    for (std::size_t k = 0; k < R; k++) {
        b[k] = distrib_double(gen);
    }

    int N = 10000;

    double e = 1e-6;

    double t = 0.0001;
    std::vector<double> x0(R, 0.0);

    solution res_mpi = MPI_solver(A_CSR, b, x0, e, t, N);
    std::fill(x0.begin(), x0.end(), 0.0);
    
    solution res_jac = Jacobi_solver(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);
    
    solution res_gs = Gauss_solver(A_CSR, b, x0, e, N);
    std::fill(x0.begin(), x0.end(), 0.0);

    int N_1 = 256;
    solution res_ch = Chebyshev_solver(A_CSR, b, x0, e, N_1);

    std::ofstream out_iter("graphs/iter_data.txt");
    size_t max_iter = std::max({res_mpi.discrepancy_p.size(), 
                                res_jac.discrepancy_p.size(), 
                                res_gs.discrepancy_p.size(),
                                res_ch.discrepancy_p.size()});
    for (size_t k = 0; k < max_iter; k++) {
        double r_mpi = (k < res_mpi.discrepancy_p.size()) ? res_mpi.discrepancy_p[k] : -1.0;
        double r_jac = (k < res_jac.discrepancy_p.size()) ? res_jac.discrepancy_p[k] : -1.0;
        double r_gs  = (k < res_gs.discrepancy_p.size())  ? res_gs.discrepancy_p[k]  : -1.0;
        double r_ch  = (k < res_ch.discrepancy_p.size())  ? res_ch.discrepancy_p[k]  : -1.0;
        out_iter << r_mpi << " " << r_jac << " " << r_gs << " " << r_ch << "\n";
    }
    out_iter.close();

    std::ofstream out_mpi("graphs/time_mpi.txt");
    for (size_t k = 0; k < res_mpi.time_p.size(); k++) {
        out_mpi << res_mpi.time_p[k] << " " << res_mpi.discrepancy_p[k] << "\n";
    }
    out_mpi.close();

    std::ofstream out_jac("graphs/time_jacobi.txt");
    for (size_t k = 0; k < res_jac.time_p.size(); k++) {
        out_jac << res_jac.time_p[k] << " " << res_jac.discrepancy_p[k] << "\n";
    }
    out_jac.close();

    std::ofstream out_gs("graphs/time_gauss.txt");
    for (size_t k = 0; k < res_gs.time_p.size(); k++) {
        out_gs << res_gs.time_p[k] << " " << res_gs.discrepancy_p[k] << "\n";
    }
    out_gs.close();

    std::ofstream out_ch("graphs/time_chebyshev.txt");
    for (size_t k = 0; k < res_ch.time_p.size(); k++) {
        out_ch << res_ch.time_p[k] << " " << res_ch.discrepancy_p[k] << "\n";
    }
    out_ch.close();
}