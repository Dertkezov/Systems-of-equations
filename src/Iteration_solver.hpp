#pragma once
#include <vector>
#include <cmath>
#include <time.h>
#include "csr_matrix.hpp"

struct solution {
    std::vector<double> x;
    int p;
    double total_time;
    std::vector<double> discrepancy_p;
    std::vector<double> time_p;
};

inline double discrepancy(const std::vector<double>& x, const std::vector<double>& y) {
    double r = 0;
    for (size_t i = 0; i < x.size(); i++) {
        r += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return(sqrt(r));
}

inline solution MPI_solver (const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, double t, int N) {
    solution s;
    s.p = 0;
    s.total_time = 0.0;

    std::vector<double> x = A * x0;
    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    while ((r > e) && (s.p < N)) {
        clock_t start = clock();
        for (size_t i = 0; i < x0.size(); i++) {
            x0[i] = x0[i] - t * (x[i] - b[i]);
        }
        x = A * x0;
        clock_t end = clock();
        double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000;

        r = discrepancy(x, b);
        s.discrepancy_p.push_back(r);
        s.total_time += dt;
        s.time_p.push_back(s.total_time);
        s.p++;
    }

    s.x = x0;
    return(s);
}

inline solution Jacobi_solver (const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, int N) {
    size_t n = x0.size();

    solution s;
    s.p = 0;
    s.total_time = 0.0;
    
    std::vector<double> x = A * x0;
    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    std::vector<double> D(n);
    for (size_t i = 0; i < n; i++) {
        D[i] = A(i, i);
    }
    
    while ((r > e) && (s.p < N)) {
        clock_t start = clock();
        for (size_t i = 0; i < x0.size(); i++) {
            x0[i] = (b[i] - (x[i] - D[i] * x0[i])) / D[i];
        }
        x = A * x0;
        clock_t end = clock();
        double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000;

        r = discrepancy(x, b);
        s.discrepancy_p.push_back(r);
        s.total_time += dt;
        s.time_p.push_back(s.total_time);
        s.p++;
    }

    s.x = x0;
    return(s);
}

inline solution Gauss_solver (const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, int N) {
    size_t n = x0.size();
    
    solution s;
    s.p = 0;
    s.total_time = 0.0;
    
    std::vector<double> x = A * x0;
    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    std::vector<double> D(n);
    for (size_t i = 0; i < n; i++) {
        D[i] = A(i, i);
    }
    
    while ((r > e) && (s.p < N)) {
        clock_t start = clock();
        for (size_t k = 0; k < x0.size(); k++) {
            double sum1 = 0.0;
            for (size_t j = k + 1; j < x0.size(); j++) {
                sum1 += A(k, j) * x0[j];
            }

            double sum2 = 0.0;
            for (size_t j = 0; j < k; j++) {
                sum2 += A(k, j) * x0[j];
            }
            x0[k] = (b[k] - sum1 - sum2) / D[k];
        }

        clock_t end = clock();
        double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000;
        x = A * x0;
        r = discrepancy(x, b);
        s.discrepancy_p.push_back(r);
        s.total_time += dt;
        s.time_p.push_back(s.total_time);
        s.p++;
    }
    
    s.x = x0;
    return(s);
}