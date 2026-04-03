#pragma once
#include <vector>
#include <cmath>
#include <time.h>
#include "csr_matrix.hpp"

inline std::vector<double> Gauss_step (const CSR_matrix& A, const std::vector<double>& D, const std::vector<double>& b, std::vector<double>& x0) {
    size_t n = x0.size();

    for (size_t k = 0; k < n; k++) {
        double sum1 = 0.0;
        for (size_t j = k + 1; j < n; j++) {
            sum1 += A(k, j) * x0[j];
        }

        double sum2 = 0.0;
        for (size_t j = 0; j < k; j++) {
            sum2 += A(k, j) * x0[j];
        }
        x0[k] = (b[k] - sum1 - sum2) / D[k];
    }

    return(x0);
}

inline std::vector<double> Symmetric_Gauss_step(const CSR_matrix& A, const std::vector<double>& D, const std::vector<double>& b, std::vector<double> x0) {
    size_t n = x0.size();
    
    for (size_t k = 0; k < n; k++) {
        double sum1 = 0.0;
        for (size_t j = k + 1; j < n; j++) {
            sum1 += A(k, j) * x0[j];
        }
        
        double sum2 = 0.0;
        for (size_t j = 0; j < k; j++) {
            sum2 += A(k, j) * x0[j];
        }
        
        x0[k] = (b[k] - sum1 - sum2) / D[k];
    }
    
    for (size_t k = n; k-- > 0; ) {
        double sum1 = 0.0;
        for (size_t j = k + 1; j < n; j++) {
            sum1 += A(k, j) * x0[j];
        }
        
        double sum2 = 0.0;
        for (size_t j = 0; j < k; j++) {
            sum2 += A(k, j) * x0[j];
        }
        
        x0[k] = (b[k] - sum1 - sum2) / D[k];
    }
    
    return(x0);
}

inline std::vector<double> Jacobi_step(const CSR_matrix& A, const std::vector<double>& D, const std::vector<double>& b, std::vector<double> x0) {
    size_t n = x0.size();

    std::vector<double> x = A * x0;

    for (size_t i = 0; i < n; i++) {
        x0[i] = (b[i] - (x[i] - A(i, i) * x0[i])) / D[i];
    }

    return(x0);
}
