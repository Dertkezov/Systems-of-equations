#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "Iteration_solver.hpp"
#include "MPI_speed.hpp"

int main() {
    size_t R;
    std::cout << "Enter size" << std::endl;
    std::cin >> R;
    
    full_matrix A(R, R);
    std::cout << "Enter matrix" << std::endl;
    for (size_t i = 0; i < R; i++) {
        for (size_t j = 0; j < R; j++) {
            std::cin >> A(i, j);
        }
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

    std::cout << "MPI_solver" << std::endl;
    double t;
    std::cout << "Enter t" << std::endl;
    std::cin >> t;

    std::vector<double> x0(R, 0.0);

    int N = 10000;
    solution res_1 = MPI_solver(A_CSR, b, x0, e, t, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_1.x[i] << " ";
    }
    std::cout<<std::endl;
    std::cout << std::endl;

    std::cout << "Jacobi_solver" << std::endl;
    solution res_2 = Jacobi_solver(A_CSR, b, x0, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_2.x[i] << " ";
    }
    std::cout<<std::endl;
    std::cout << std::endl;

    std::cout << "Gauss_solver" << std::endl;
    solution res_3 = Gauss_solver(A_CSR, b, x0, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_3.x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Chebyshev_solver" << std::endl;
    int N_1;
    std::cout << "Enter N" << std::endl;
    std::cin >> N_1;
    solution res_4 = Chebyshev_solver(A_CSR, b, x0, e, N_1);
    for (size_t i = 0; i < R; i++){
        std::cout << res_4.x[i] << " ";
    }
}