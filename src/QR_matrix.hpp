#pragma once
#include "matrix.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

class QR_matrix {
private:
    full_matrix Q;
    full_matrix R;

public:
    void find(const full_matrix& A) {
        size_t m = A.rows();
        size_t n = A.cols();

        R = full_matrix(m, n);
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                R(i, j) = A(i, j);
            }
        }

        Q = full_matrix(m, m);
        for (size_t i = 0; i < m; i++) {
            Q(i, i) = 1.0;
        }

        for (size_t k = 0; k < std::min(m, n); k++) {

            std::vector<double> x(m - k, 0.0);
            for (size_t i = 0; i < m - k; i++) {
                x[i] = R(k + i, k);
            }

            double L_x = 0;
            for (size_t i = 0; i < x.size(); i++) {
                L_x += x[i] * x[i];
            }
            L_x = std::sqrt(L_x);

            double sign = (x[0] >= 0) ? 1.0 : -1.0;
            std::vector<double> v = x;
            v[0] += sign * L_x;

            double vTv = 0.0;
            for (size_t i = 0; i < v.size(); i++) {
                vTv += v[i] * v[i];
            }

            for (size_t j = k; j < n; j++) {
                double w = 0.0;
                for (size_t i = 0; i < m - k; i++) {
                    w += v[i] * R(k + i, j);
                }
                for (size_t i = 0; i < m - k; i++) {
                    R(k + i, j) -= 2 * v[i] * w / vTv;
                }
            }

            std::vector<double> v_m(m, 0.0);
            for (size_t i = 0; i < v.size(); i++) {
                v_m[k + i] = v[i];
            }

            for (size_t i = 0; i < m; i++) {
                double w = 0.0;
                for (size_t j = 0; j < m; j++) {
                    w += v_m[j] * Q(i, j);
                }
                for (size_t j = 0; j < m; j++) {
                    Q(i, j) -= 2 * v_m[j] * w / vTv;
                }
            }
        }
    }

    const full_matrix& getQ() const {
        return Q;
    }
    const full_matrix& getR() const {
        return R;
    }

    std::vector<double> QR_solve(const full_matrix& Q, const full_matrix& R, const std::vector<double>& b) {
        size_t m = Q.rows();
        size_t n = R.cols();

        std::vector<double> y(m, 0.0);
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < m; j++) {
                y[i] += Q(j, i) * b[j];
            }
        }

        std::vector<double> x(n, 0.0);
        for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < static_cast<int>(n); j++) {
                sum += R(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / R(i, i);
        }
        return x;
    };
};