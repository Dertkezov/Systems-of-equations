#pragma once
#include <vector>
#include <iostream>

inline std::vector<double> x_solution(std::vector<std::vector<double>>& A, std::vector<double>& b){
    int N = A.size();
    std::vector<double> x(N, 0.0);
    std::vector<double> p(N - 1, 0.0);
    std::vector<double> q(N - 1, 0.0);

    p[0] = -A[0][1] / A[0][0];
    q[0] = b[0] / A[0][0];

    for (int i = 1; i < N - 1; i++) {
        p[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);
        q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] + A[i][i - 1] * p[i - 1]);
    }

    x[N - 1] = (b[N - 1] - A[N - 1][N - 2] * q[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * p[N - 2]);

    for (int i = N - 2; i >= 0; i--) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return (x);
}