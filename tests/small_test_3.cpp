#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "csr_matrix.hpp"
#include "full_matrix.hpp"
#include "Iteration_solver.hpp"
#include "MPI_speed.hpp"
#include "Chebyshev_speed_symmetry.hpp"
#include "GMRES.hpp"
#include "pre_CG.hpp"
#include "matrix_for_preconditioning.hpp"

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

TEST(Iteration_solver, Chebyshev_Solution) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 8;

    solution res = Chebyshev_solver(A_CSR, b, x0, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, Chebyshev_acceleration_Jacobi_step) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 5.0; A(0, 2) = 1.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 2.0; A(2, 2) = 30.0;
    std::vector<double> D(3);
    for (size_t i = 0; i < 3; i++) {
        D[i] = A(i, i);
    }

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 33.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    auto jacobi_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Jacobi_step(A_CSR, D, b, x);
    };
    double rho = estimate_rho(jacobi_step, x0);
    solution res = Chebyshev_accelerate(jacobi_step, A_CSR, b, x0, rho, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, Chebyshev_acceleration_Gauss_step) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 5.0; A(0, 2) = 1.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 2.0; A(2, 2) = 30.0;
    std::vector<double> D(3);
    for (size_t i = 0; i < 3; i++) {
        D[i] = A(i, i);
    }

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 33.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    auto gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Gauss_step(A_CSR, D, b, x);
    };
    double rho = estimate_rho(gauss_step, x0);
    solution res = Chebyshev_accelerate(gauss_step, A_CSR, b, x0, rho, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, Chebyshev_acceleration_Symmetric_Gauss_step) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 5.0; A(0, 2) = 1.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 2.0; A(2, 2) = 30.0;
    std::vector<double> D(3);
    for (size_t i = 0; i < 3; i++) {
        D[i] = A(i, i);
    }

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 33.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    auto sym_gauss_step = [&A_CSR, &D, &b](std::vector<double> x) {
        return Symmetric_Gauss_step(A_CSR, D, b, x);
    };
    double rho = estimate_rho(sym_gauss_step, x0);
    solution res = Chebyshev_accelerate(sym_gauss_step, A_CSR, b, x0, rho, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, SOR_solution) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;
    double omega = 1.5;

    solution res = SOR_solver(A_CSR, b, x0, e, omega, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, fast_descent_solver) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    solution res = fast_descent(A_CSR, b, x0, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, CG_solver) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    solution res = CG_solver(A_CSR, b, x0, e, N);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, GMRES_solver) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 4.0; A(0, 2) = 2.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 34.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;
    int m = 5;

    solution res = GMRES_solver(A_CSR, b, x0, e, N, m);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}

TEST(Iteration_solver, PCG_solver) {
    full_matrix A(3, 3);
    A(0, 0) = 10.0; A(0, 1) = 5.0; A(0, 2) = 1.0;
    A(1, 0) = 5.0; A(1, 1) = 20.0; A(1, 2) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 2.0; A(2, 2) = 30.0;

    CSR_matrix A_CSR(A);
    std::vector<double> b = {16.0, 27.0, 33.0};
    std::vector<double> x0 = {0.0, 0.0, 0.0};

    double e = 1e-8;
    int N = 1000;

    CSR_matrix L = Cholesky_factor(A);
    solution res = PCG_solver(A_CSR, b, x0, e, N, L);

    EXPECT_NEAR(res.x[0], 1.0, 1e-6);
    EXPECT_NEAR(res.x[1], 1.0, 1e-6);
    EXPECT_NEAR(res.x[2], 1.0, 1e-6);
}