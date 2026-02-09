#pragma once
#include <vector>
#include <iostream>

inline std::vector<double> x_solution(std::vector<std::vector<double>>& A, std::vector<double>& b, int N){
    std::vector<double> x(N, 0.0);
    std::vector<double> p(N, 0.0);
    std::vector<double> q(N, 0.0);

    p[1] = -A[2][0] / A[1][0];
    q[1] = b[0] / A[1][0];

    for (int i = 1; i < N - 1; i++) {
        p[i + 1] = -A[2][i] / (A[1][i] + A[0][i] * p[i]);
        q[i + 1] = (b[i] - A[0][i] * q[i]) / (A[1][i] + A[0][i] * p[i]);
    }

    x[N - 1] = (b[N - 1] - A[0][N - 1] * q[N - 1]) / (A[1][N - 1] + A[0][N - 1] * p[N - 1]);

    for (int i = N - 1; i >= 1; i--) {
        x[i - 1] = p[i] * x[i] + q[i];
    }

    return (x);
}
