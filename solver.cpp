#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "Iteration_solver.hpp"

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
    
    double t;
    std::cout << "Enter t" << std::endl;
    std::cin >> t;
    std::cout << std::endl;

    double e;
    std::cout << "Enter e" << std::endl;
    std::cin >> e;
    std::cout << std::endl;

    std::vector<double> x0(R, 0.0);

    int N = 10000;
    std::cout << "MPI_solver" << std::endl;
    solution res_1 = MPI_solver(A_CSR, b, x0, e, t, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_1.x[i] << " ";
    }
    std::cout<<std::endl;

    std::cout << "Jacobi_solver" << std::endl;
    solution res_2 = Jacobi_solver(A_CSR, b, x0, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_2.x[i] << " ";
    }
    std::cout<<std::endl;

    std::cout << "Gauss_solver" << std::endl;
    solution res_3 = Gauss_solver(A_CSR, b, x0, e, N);
    for (size_t i = 0; i < R; i++){
        std::cout << res_3.x[i] << " ";
    }
    std::cout<<std::endl;
}
