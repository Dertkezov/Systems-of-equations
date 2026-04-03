#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "Iteration_solver.hpp"
#include "Chebyshev_speed_symmetry.hpp"

int main() {
    std::cout << "Acceleration of symmetric methods" << std::endl;
    std::cout<<std::endl;
    size_t R;
    std::cout << "Enter size" << std::endl;
    std::cin >> R;
    
    full_matrix A(R, R);
    std::cout << "Enter symmetric matrix" << std::endl;
    for (size_t i = 0; i < R; i++) {
        for (size_t j = 0; j < R; j++) {
            std::cin >> A(i, j);
        }
    }
    std::vector<double> D(R);
    for (size_t i = 0; i < R; i++) {
        D[i] = A(i, i);
    }

    CSR_matrix A_CSR(A);
    
    std::vector<double> b(R);
    std::cout << "Enter vector" << std::endl;
    for (size_t i = 0; i < R; i++) {
        std::cin >> b[i];
    }

    double e;
    std::cout << "Enter e" << std::endl;
    std::cin >> e;
    std::cout << std::endl;

    std::vector<double> x0(R, 0.0);
    int N = 10000;

    auto jacobi_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Jacobi_step(A_CSR, D, b, x);
    };
    
    auto gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Gauss_step(A_CSR, D, b, x);
    };
    
    auto sym_gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Symmetric_Gauss_step(A_CSR, D, b, x);
    };


    std::cout << std::endl;
    std::cout << "Jacobi_step" << std::endl;
    std::cout << std::endl;
    double rho1 = estimate_rho(jacobi_step, x0);
    solution res_1 = Chebyshev_accelerate(jacobi_step, A_CSR, b, x0, rho1, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_1.x[i] << " ";
    }
    std::cout<<std::endl;
    std::cout<<std::endl;

    std::cout << "Gauss_step" << std::endl;
    std::cout << std::endl;
    double rho2 = estimate_rho(gauss_step, x0);
    solution res_2 = Chebyshev_accelerate(gauss_step, A_CSR, b, x0, rho2, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_2.x[i] << " ";
    }
    std::cout<<std::endl;
    std::cout<<std::endl;

    std::cout << "Symmetry_Gauss_step" << std::endl;
    std::cout << std::endl;
    double rho3 = estimate_rho(sym_gauss_step, x0);
    solution res_3 = Chebyshev_accelerate(sym_gauss_step, A_CSR, b, x0, rho3, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_3.x[i] << " ";
    }
}
