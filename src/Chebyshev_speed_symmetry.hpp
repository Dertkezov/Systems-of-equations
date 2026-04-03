#pragma once
#include "Iteration_solver.hpp"
#include <cmath>
#include <algorithm>
#include "csr_matrix.hpp"
#include "full_matrix.hpp"
#include "Iteration_solver.hpp"
#include "Steps_for_chebyshev.hpp"

template<typename StepFunc>
double estimate_rho(StepFunc step, std::vector<double> x0) {
    std::vector<double> x1 = step(x0);
    std::vector<double> x2;
    
    double rho = 0.0;
    
    for (int i = 0; i < 100; i++) {
        x2 = step(x1);
        
        double num = discrepancy(x2, x1);
        double den = discrepancy(x1, x0);
        
        if (den > 1e-10) {
            rho = num / den;
        }
        
        x0 = x1;
        x1 = x2;
    }
    
    if (rho < 0.01) {
        rho = 0.01;
    }

    if (rho > 0.999) {
        rho = 0.999;
    }
    
    if (std::isnan(rho)) {
        rho = 0.2;
    } 
    
    return rho;
}


template<typename StepFunc>
solution Chebyshev_accelerate(StepFunc step, const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& y_0, double rho, double e, int N) {
    solution s;
    s.p = 0;
    s.total_time = 0.0;
    
    size_t n = y_0.size();
    
    double r = discrepancy(A * y_0, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);
    
    double w1 = 1.0;
    double w2 = 2.0 / (2.0 - rho * rho);
    std::vector<double> y_1 = step(y_0);

    while ((r > e) && (s.p < N)) {
        clock_t start = clock();

        std::vector<double> y_step = step(y_1);
        std::vector<double> y_2(n);

        for (size_t i = 0; i < n; i++) {
            y_2[i] = w2 * (y_step[i] - y_0[i]) + y_0[i];
        }
        
        y_0 = y_1;
        y_1 = y_2;
        
        w1 = w2;
        w2 = 1.0 / (1.0 - (rho * rho * w1) / 4.0);

        clock_t end = clock();
        double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;

        r = discrepancy(A * y_1, b);
        s.discrepancy_p.push_back(r);
        s.total_time += dt;
        s.time_p.push_back(s.total_time);
        s.p++;
    }
    
    s.x = y_1;
    return s;
}