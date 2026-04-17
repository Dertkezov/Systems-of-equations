#pragma once
#include <vector>
#include <cmath>
#include <time.h>
#include "csr_matrix.hpp"
#include <algorithm>
#include "Iteration_solver.hpp"
#include <span>

struct ColMajorMatrix {
    size_t rows;
    size_t cols;
    std::vector<double> data;

    std::span<const double> col(size_t j) const noexcept {
        return std::span<const double>(data.data() + j * rows, rows);
    }
    std::span<double> col(size_t j) noexcept {
        return std::span<double>(data.data() + j * rows, rows);
    }
};

inline std::vector<double> discrepancy_vec(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    std::vector<double> r(n);
    for (size_t i = 0; i < x.size(); i++) {
        r[i] = x[i] - y[i];
    }
    return(r);
}

inline double vec_norm(const std::vector<double>& x) {
    double r = 0;
    for (size_t i = 0; i < x.size(); i++) {
        r += x[i] * x[i];
    }
    return(sqrt(r));
}

inline solution GMRES_solver(const CSR_matrix& A, const std::vector<double>& b, std::vector<double>& x0, double e, int N, int m) {
    size_t n = x0.size();

    solution s;
    s.p = 0;
    s.total_time = 0.0;

    std::vector<double> x = A * x0;

    std::vector<double> r_vec = discrepancy_vec(x, b);

    double r = discrepancy(x, b);
    s.discrepancy_p.push_back(r);
    s.time_p.push_back(0.0);

    ColMajorMatrix V;
    V.rows = n;
    V.cols = m + 1;
    V.data.resize(n * (m + 1));

    std::vector<double> H((m + 1) * m, 0.0);
    std::vector<std::pair<double, double>> rotations;

    std::vector<double> z_i(m + 1, 0.0);

    while ((r > e) && (s.p < N)) {
        rotations.clear();
        std::fill(z_i.begin(), z_i.end(), 0.0);
        z_i[0] = r;

        auto v0 = V.col(0);
        for (size_t i = 0; i < n; i++){
            v0[i] = r_vec[i] / r;
        } 

        for (int j = 0; j < m; j++) {
            clock_t start = clock();

            std::vector<double> v_j(V.col(j).begin(), V.col(j).end());
            std::vector<double> w = A * v_j;

            for (int i = 0; i <= j; i++) {
                auto v_i = V.col(i);
                double h_ij = 0.0;
                for (size_t k = 0; k < n; k++){
                    h_ij += v_i[k] * w[k];
                }

                H[j * (m + 1) + i] = h_ij;

                for (size_t k = 0; k < n; k++){
                    w[k] -= h_ij * v_i[k];
                }
            }

            double h_1 = vec_norm(w);
            H[j * (m + 1) + j + 1] = h_1;

            auto v_1 = V.col(j + 1);

            for (size_t k = 0; k < n; k++) {
                v_1[k] = w[k] / h_1;
            }

            for (int i = 0; i < j; i++) {
                double a = rotations[i].first;
                double b = rotations[i].second;
                double h1 = H[j * (m + 1) + i];
                double h2 = H[j * (m + 1) + i + 1];
                H[j * (m + 1) + i] = a * h1 + b * h2;
                H[j * (m + 1) + i + 1] = -b * h1 + a * h2;
            }

            double h1 = H[j * (m + 1) + j];
            double h2 = H[j * (m + 1) + j + 1];
            double a, b;
            double d = std::sqrt(h1 * h1 + h2 * h2);
            a = h1 / d;
            b = h2 / d;
            H[j * (m + 1) + j] = d;
            H[j * (m + 1) + j + 1] = 0.0;
            rotations.push_back({a, b});

            double z1 = z_i[j];
            double z2 = z_i[j + 1];
            z_i[j] = a * z1 + b * z2;
            z_i[j + 1] = -b * z1 + a * z2;

            clock_t end = clock();
            double dt = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
            s.total_time += dt;
            s.p++;

            r = std::abs(z_i[j + 1]);
            s.discrepancy_p.push_back(r);
            s.time_p.push_back(s.total_time);

            if (r < e){
                break;
            }
        }

        int k = rotations.size();
        std::vector<double> y(k, 0.0);
        for (int i = k - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int col = i + 1; col < k; col++) {
                sum += H[col * (m + 1) + i] * y[col];
            }
            y[i] = (z_i[i] - sum) / H[i * (m + 1) + i];
        }

        for (int i = 0; i < k; i++) {
            double yi = y[i];
            auto v_i = V.col(i);
            for (size_t row = 0; row < n; row++) {
                x0[row] -= v_i[row] * yi;
            }
        }

        x = A * x0;
        r_vec = discrepancy_vec(x, b);
        r = discrepancy(x, b);
    }

    s.x = x0;
    return s;
}