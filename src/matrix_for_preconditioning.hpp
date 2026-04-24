#pragma once
#include "full_matrix.hpp"
#include "csr_matrix.hpp"
#include <cmath>
#include <algorithm>

inline full_matrix lvl_matrix(const full_matrix& A) {
    size_t n = A.rows();
    full_matrix lvl(n, n);
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (A(i, j) != 0.0) {
                lvl(i, j) = 0.0;
            }
            else {
                lvl(i, j) = static_cast<double>(n - 1);
                for (size_t k = 0; k < i; k++) {
                    double val = lvl(i, k) + lvl(k, j) + 1.0;
                    if (val < lvl(i, j)) {
                        lvl(i, j) = val;
                    }
                }
            }
        }
    }
    return lvl;
}


inline CSR_matrix Cholesky_factor(const full_matrix& A) {
    size_t n = A.rows();
    full_matrix L(n, n);
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            if (A(i, j) != 0.0) {
                double sum = A(i, j);
                for (size_t k = 0; k < j; k++) {
                    if ((A(i, k) != 0.0) && (A(j, k) != 0.0)) {
                        sum -= L(i, k) * L(j, k);
                    }
                }
                if (i == j) {
                    L(i, i) = std::sqrt(sum);
                } 
                else {
                    L(i, j) = sum / L(j, j);
                }
            }
        }
    }
    return CSR_matrix(L);
}
