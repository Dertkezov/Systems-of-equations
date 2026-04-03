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

int main() {
    std::size_t R = 10000;

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
        double a = distrib_double(gen);
        A(i, j) = a;
        A(j, i) = a;
    }

    double val = 1000000;
    std::vector<double> D(R);
    for (size_t k = 0; k < R; k++) {
        A(k, k) = val; 
        D[k] = val;
    }

    CSR_matrix A_CSR(A);

    std::vector<double> b(R, 0.0);
    for (std::size_t k = 0; k < R; k++) {
        b[k] = distrib_double(gen);
    }

    int N = 2000;

    double e = 1e-6;

    std::vector<double> x0(R, 0.0);

    auto jacobi_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Jacobi_step(A_CSR, D, b, x);
    };
    
    auto gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Gauss_step(A_CSR, D, b, x);
    };
    
    auto sym_gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Symmetric_Gauss_step(A_CSR, D, b, x);
    };

    double rho1 = estimate_rho(jacobi_step, x0);
    solution res_1 = Chebyshev_accelerate(jacobi_step, A_CSR, b, x0, rho1, e, N);

    double rho2 = estimate_rho(gauss_step, x0);
    solution res_2 = Chebyshev_accelerate(gauss_step, A_CSR, b, x0, rho2, e, N);

    double rho3 = estimate_rho(sym_gauss_step, x0);
    solution res_3 = Chebyshev_accelerate(sym_gauss_step, A_CSR, b, x0, rho3, e, N);

    std::ofstream out_iter("graphs/iter_data_sym.txt");
    size_t max_iter = std::max({res_1.discrepancy_p.size(), 
                                res_2.discrepancy_p.size(), 
                                res_3.discrepancy_p.size()});
    for (size_t k = 0; k < max_iter; k++) {
        double r_1 = (k < res_1.discrepancy_p.size()) ? res_1.discrepancy_p[k] : -1.0;
        double r_2 = (k < res_2.discrepancy_p.size()) ? res_2.discrepancy_p[k] : -1.0;
        double r_3 = (k < res_3.discrepancy_p.size()) ? res_3.discrepancy_p[k] : -1.0;
        out_iter << r_1 << " " << r_2 << " " << r_3 << "\n";
    }
    out_iter.close();

    std::ofstream out_jac("graphs/time_jac_speed.txt");
    for (size_t k = 0; k < res_1.time_p.size(); k++) {
        out_jac << res_1.time_p[k] << " " << res_1.discrepancy_p[k] << "\n";
    }
    out_jac.close();

    std::ofstream out_gauss("graphs/time_gauss_speed.txt");
    for (size_t k = 0; k < res_2.time_p.size(); k++) {
        out_gauss << res_2.time_p[k] << " " << res_2.discrepancy_p[k] << "\n";
    }
    out_gauss.close();

    std::ofstream out_gauss_sym("graphs/time_gauss_sym_speed.txt");
    for (size_t k = 0; k < res_3.time_p.size(); k++) {
        out_gauss_sym << res_3.time_p[k] << " " << res_3.discrepancy_p[k] << "\n";
    }
    out_gauss_sym.close();
}
