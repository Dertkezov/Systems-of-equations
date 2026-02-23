#include <gtest/gtest.h>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"

TEST(FullMatrixTest, Creation) {
    full_matrix A1(3, 3);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            A1(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }
    EXPECT_EQ(A1.rows(), 3);
    EXPECT_EQ(A1.cols(), 3);
}

TEST(FullMatrixTest, MultiplyByScalar) {
    full_matrix A1(3, 3);
    A1(0, 0) = 1.0; A1(0, 1) = 2.0;
    A1(1, 1) = 4.0; A1(1, 2) = 2.0;
    A1(2, 2) = 6.0;

    full_matrix B1 = A1 * 2.0;
    EXPECT_DOUBLE_EQ(B1(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(B1(1, 1), 8.0);
    EXPECT_DOUBLE_EQ(B1(2, 2), 12.0);
}

TEST(FullMatrixTest, MultiplyByVector) {
    full_matrix A1(3, 3);
    A1(0, 0) = 1.0; A1(0, 1) = 2.0;
    A1(1, 1) = 4.0; A1(1, 2) = 2.0;
    A1(2, 2) = 6.0;

    std::vector<double> x1 = {1.0, 2.0, 3.0};
    std::vector<double> y1 = A1 * x1;

    EXPECT_DOUBLE_EQ(y1[0], 5.0);
    EXPECT_DOUBLE_EQ(y1[1], 14.0);
    EXPECT_DOUBLE_EQ(y1[2], 18.0);
}

TEST(CSRMatrixTest, FromFullMatrix) {
    full_matrix A1(3, 3);
    A1(0, 0) = 1.0; A1(0, 1) = 2.0;
    A1(1, 1) = 4.0; A1(1, 2) = 2.0;
    A1(2, 2) = 6.0;

    CSR_matrix A1_CSR(A1);
    EXPECT_DOUBLE_EQ(A1_CSR(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A1_CSR(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(A1_CSR(2, 2), 6.0);
}

TEST(CSRMatrixTest, MultiplyByVector) {
    full_matrix A1(3, 3);
    A1(0, 0) = 1.0; A1(0, 1) = 2.0;
    A1(1, 1) = 4.0; A1(1, 2) = 2.0;
    A1(2, 2) = 6.0;

    CSR_matrix A1_CSR(A1);
    std::vector<double> x1 = {1.0, 2.0, 3.0};
    std::vector<double> y1 = A1_CSR * x1;

    EXPECT_DOUBLE_EQ(y1[0], 5.0);
    EXPECT_DOUBLE_EQ(y1[1], 14.0);
    EXPECT_DOUBLE_EQ(y1[2], 18.0);
}

TEST(CSRMatrixTest, FromArrays) {
    std::vector<double> v = {1.0, 2.0, 4.0, 2.0, 6.0};
    std::vector<size_t> c = {0, 1, 1, 2, 2};
    std::vector<size_t> r = {0, 2, 4, 5};

    CSR_matrix A2(3, 3, v, c, r);
    EXPECT_DOUBLE_EQ(A2(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(A2(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(A2(2, 2), 6.0);
}

TEST(CSRMatrixTest, MultiplyByScalar) {
    std::vector<double> v = {1.0, 2.0, 4.0, 2.0, 6.0};
    std::vector<size_t> c = {0, 1, 1, 2, 2};
    std::vector<size_t> r = {0, 2, 4, 5};

    CSR_matrix A2(3, 3, v, c, r);
    CSR_matrix B2 = A2 * 3.0;

    EXPECT_DOUBLE_EQ(B2(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(B2(1, 1), 12.0);
    EXPECT_DOUBLE_EQ(B2(2, 2), 18.0);
}

TEST(CSRMatrixTest, MultiplyByVectorFromArrays) {
    std::vector<double> v = {1.0, 2.0, 4.0, 2.0, 6.0};
    std::vector<size_t> c = {0, 1, 1, 2, 2};
    std::vector<size_t> r = {0, 2, 4, 5};

    CSR_matrix A2(3, 3, v, c, r);
    std::vector<double> x2 = {3.0, 2.0, 1.0};
    std::vector<double> y2 = A2 * x2;

    EXPECT_DOUBLE_EQ(y2[0], 7.0);
    EXPECT_DOUBLE_EQ(y2[1], 10.0);
    EXPECT_DOUBLE_EQ(y2[2], 6.0);
}

TEST(CSRMatrixTest, DenseVsCSR) {
    full_matrix A1(3, 3);
    A1(0, 0) = 1.0; A1(0, 1) = 2.0;
    A1(1, 1) = 4.0; A1(1, 2) = 2.0;
    A1(2, 2) = 6.0;

    CSR_matrix A1_CSR(A1);
    std::vector<double> x1 = {1.0, 2.0, 3.0};

    std::vector<double> y1_dense = A1 * x1;
    std::vector<double> y1_csr = A1_CSR * x1;

    EXPECT_EQ(y1_dense.size(), y1_csr.size());
    for (size_t i = 0; i < y1_dense.size(); i++) {
        EXPECT_DOUBLE_EQ(y1_dense[i], y1_csr[i]);
    }
}