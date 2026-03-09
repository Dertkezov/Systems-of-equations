#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "csr_matrix.hpp"
#include "full_matrix.hpp"
#include "Iteration_solver.hpp"

TEST(Iteration_solver, MPI_Solution) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    double t = 0.02;
    int N = 10000;

    solution res = MPI_solver(A_CSR, b, x0, e, t, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, Jacobi_Solution) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 10000;

    solution res = Jacobi_solver(A_CSR, b, x0, e, N);
    
    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, GaussSeidel_Solution) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 10000;

    solution res = Gauss_solver(A_CSR, b, x0, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}