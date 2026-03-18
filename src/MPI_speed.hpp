#pragma once
#include <vector>
#include <cmath>
#include <time.h>
#include "csr_matrix.hpp"
#include <algorithm>
#include <iostream>
#include "Iteration_solver.hpp"

double find_lambda_max (const CSR_matrix& A, size_t n) {
    std::vector<double> r0(n, 1.0);

    for (int i = 0; i < 10000; i++) {
        std::vector<double> r = A * r0;
        
        double norm_r = 0.0;
        for (size_t j = 0; j < n; j++) {
            norm_r += r[j] * r[j];
        }
        norm_r = sqrt(norm_r);

        for (size_t j = 0; j < n; j++) {
            r0[j] = r[j] / norm_r;
        }
    }
    std::vector<double> r = A * r0;
    double num, den;
    num = 0;
    den = 0;
    for (size_t j = 0; j < n; j++) {
        num += r0[j] * r[j];
        den += r0[j] * r0[j];
    }
    return(num / den);
}

void shafl (int n0, std::vector<int>& p0) {
    if (n0 == 1){
        p0 = {0};
        return;
    }

    std::vector<int> p;
    int n = n0 / 2;
    shafl(n, p);

    p0.resize(n0);
    for (int i = 0; i < n; i++) {
        int k = p[i];
        p0[2 * i] = k;
        p0[2 * i + 1] = (n0 - 1) - k;
    }
}

solution Chebyshev_solver(const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, int N) {
    solution s;
    s.p = 0;
    s.total_time = 0.0;

    size_t n = x0.size();
    double lambda_max = find_lambda_max(A, n);
    double lambda_min = lambda_max / 2.0;
    std::cout << lambda_max << std::endl;

    std::vector<int> p;
    shafl(N, p);
    double s_n = sin(M_PI / N);
    double c_n = cos(M_PI / N);
    double t0 = cos(M_PI / (2 * N));
    std::vector<double> taus_0(N);
    taus_0[0] = t0;
    for (size_t i = 1; i < static_cast<size_t>(N); i++){
        taus_0[i] = taus_0[i - 1] * c_n - sqrt(1 - taus_0[i - 1] * taus_0[i - 1]) * s_n;
    }
    std::vector<double> taus(N);
    for (size_t i = 0; i < static_cast<size_t>(N); i++){
        taus[i] = taus_0[p[i]];
    }

    std::vector<double> x = A * x0;
    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    while ((r > e) && (s.p < 10000)) {
        clock_t start = clock();
        double t = (lambda_max + lambda_min) / 2.0 + (lambda_max - lambda_min) / 2.0 * taus[s.p % N]; 
        for (size_t i = 0; i < x0.size(); i++) {
            x0[i] = x0[i] - (1 / t) * (x[i] - b[i]);
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