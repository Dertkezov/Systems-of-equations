#pragma once
#include "Iteration_solver.hpp"
#include "csr_matrix.hpp"
#include <vector>
#include <cmath>
#include <time.h>

inline std::vector<double> find_w (const CSR_matrix& L, const std::vector<double>& r) {
    size_t n = r.size();
    std::vector<double> z(n);

    for (size_t i = 0; i < n; i++) {
        double sum = r[i];
        for (size_t j = 0; j < i; j++) {
            double val = L(i, j);
            if (val != 0.0){
                sum -= val * z[j];
            }
        }
        z[i] = sum / L(i, i);
    }

    for (size_t i = n; i > 0; --i) {
        size_t row = i - 1;
        double sum = z[row];
        for (size_t k = row + 1; k < n; ++k) {
            double val = L(k, row);
            if (val != 0.0) sum -= val * z[k];
        }
        z[row] = sum / L(row, row);
    }
    return z;
}

inline solution PCG_solver(const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, int N, const CSR_matrix& L) {
    size_t n = x0.size();
    solution s;
    s.p = 0;
    s.total_time = 0.0;

    std::vector<double> r_vec(n), w_0(n), d(n), q(n);
    std::vector<double> x = A * x0;
    for (size_t i = 0; i < n; i++){
        r_vec[i] = b[i] - x[i];
    }

    w_0 = find_w(L, r_vec);
    d = w_0;

    double r_w_0 = 0.0;
    for (size_t i = 0; i < n; i++){
        r_w_0 += r_vec[i] * w_0[i];
    }

    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    while (r > e && s.p < N) {
        clock_t start = clock();

        q = A * d;
        double dTq = 0.0;
        for (size_t i = 0; i < n; i++){
            dTq += d[i] * q[i];
        }

        double alpha = r_w_0 / dTq;
        for (size_t i = 0; i < n; i++){
            x0[i] += alpha * d[i];
        } 
        for (size_t i = 0; i < n; i++){
            r_vec[i] -= alpha * q[i];
        } 

        w_0 = find_w(L, r_vec);
        double r_w_1 = 0.0;
        for (size_t i = 0; i < n; i++){
            r_w_1 += r_vec[i] * w_0[i];
        }

        double beta = r_w_1 / r_w_0;
        for (size_t i = 0; i < n; i++){
            d[i] = w_0[i] + beta * d[i];
        }
        r_w_0 = r_w_1;

        x = A * x0;
        r = discrepancy(x, b);

        clock_t end = clock();
        double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;

        s.discrepancy_p.push_back(r);
        s.total_time += dt;
        s.time_p.push_back(s.total_time);
        s.p++;
    }

    s.x = x0;
    return s;
}
