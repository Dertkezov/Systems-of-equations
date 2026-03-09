#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    size_t rows, cols;
    std::cout << "Enter number of rows" << std::endl;
    std::cin >> rows;
    std::cout << "Enter number of cols" << std::endl;
    std::cin >> cols;
    
    full_matrix A(rows, cols);
    std::cout << "Enter matrix" << std::endl;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            std::cin >> A(i, j);
        }
    }

    CSR_matrix A_CSR(A);
    
    std::vector<double> x(cols);
    std::cout << "Enter vector" << std::endl;
    for (size_t i = 0; i < cols; i++) {
        std::cin >> x[i];
    }
    
    double t;
    std::cout << "Enter number" << std::endl;
    std::cin >> t;
    std::cout << std::endl;

    std::cout << "Result (A * t)" << std::endl;
    full_matrix A_t = A * t;
    A_t.output();
    
    std::vector<double> y1 = A * x;
    std::cout << "Result (A * x)" << std::endl;
    for (size_t i = 0; i < y1.size(); i++) {
        std::cout << y1[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout << "Result (A_CSR * t)" << std::endl;
    CSR_matrix A_CSR_t = A_CSR * t;
    A_CSR_t.output();
    
    std::cout << "Result (A_CSR * x)" << std::endl;
    std::vector<double> y2 = A_CSR * x;
    for (size_t i = 0; i < y2.size(); i++) {
        std::cout << y2[i] << " ";
    }
}