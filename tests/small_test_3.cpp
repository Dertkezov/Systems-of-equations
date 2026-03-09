#include <gtest/gtest.h>
#include "matrix.hpp"
#include "QR_matrix.hpp"

TEST(QR_test, QR_from_A) {
    full_matrix A(3, 3);
    A(0, 0) = 1.0; A(0, 1) = 1.0; A(0, 2) = 1.0;
    A(1, 0) = 1.0; A(1, 1) = 2.0; A(1, 2) = 3.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 6.0;

    QR_matrix QR;
    QR.find(A);
    full_matrix Q = QR.getQ();
    full_matrix R = QR.getR();

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            double a_ij = 0.0;
            for (size_t k = 0; k < 3; k++) {
                a_ij += Q(k, i) * Q(k, j);
            }
            if (i == j){
                EXPECT_NEAR(a_ij, 1, 1e-10);
            }
            else{
                EXPECT_NEAR(a_ij, 0, 1e-10);
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < i; j++) {
            EXPECT_NEAR(R(i, j), 0.0, 1e-10);
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            double a_ij = 0.0;
            for (size_t k = 0; k < 3; k++) {
                a_ij += Q(i, k) * R(k, j);
            }
            EXPECT_NEAR(A(i, j), a_ij, 1e-10);
        }
    }
}

TEST(QR_solver_test, QR_solution) {
    full_matrix A(3, 3);
    A(0, 0) = 1.0; A(0, 1) = 1.0; A(0, 2) = 1.0;
    A(1, 0) = 1.0; A(1, 1) = 2.0; A(1, 2) = 3.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0; A(2, 2) = 6.0;

    std::vector<double> b = {3.0, 6.0, 10.0};
    
    QR_matrix QR;
    QR.find(A);
    full_matrix Q = QR.getQ();
    full_matrix R = QR.getR();
    std::vector<double> x = QR.QR_solve(Q, R, b);

    EXPECT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], 1.0, 1e-10);
    EXPECT_NEAR(x[1], 1.0, 1e-10);
    EXPECT_NEAR(x[2], 1.0, 1e-10);
}