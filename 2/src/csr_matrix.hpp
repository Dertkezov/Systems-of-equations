#pragma once
#include <iostream>
#include <vector>
#include "full_matrix.hpp"

class CSR_matrix {
private:
    size_t nx_;
    size_t ny_;
    std::vector<double> v_;
    std::vector<size_t> c_;
    std::vector<size_t> r_;

public:
    CSR_matrix(size_t nx, size_t ny,
               const std::vector<double>& v,
               const std::vector<size_t>& c,
               const std::vector<size_t>& r)
        : nx_(nx), ny_(ny), v_(v), c_(c), r_(r) {
    }

    CSR_matrix(const full_matrix& A) : nx_(A.rows()), ny_(A.cols()) {
        size_t k = 0;
        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) {
                if (A(i, j) != 0.0) {
                    k++;
                }
            }
        }
        v_.reserve(k);
        c_.reserve(k);
        r_.resize(nx_ + 1); 
        r_[0] = 0;
        
        int m = 0;
        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) {
                if (A(i, j) != 0.0) {
                    v_.push_back(A(i, j));
                    c_.push_back(j);
                    m++;
                }
            }
            r_[i + 1] = m;
        }
    }

    double operator()(size_t i, size_t j) const{
        size_t a = r_[i];
        size_t b = r_[i + 1];

        for (size_t k = a; k < b; k++) {
            if (c_[k] == j) {
                return(v_[k]);
            }
        }
        return(0.0);
    }

    CSR_matrix operator*(double t) const {
        CSR_matrix res(nx_, ny_, v_, c_, r_);
        
        for (size_t i = 0; i < res.v_.size(); i++) {
            res.v_[i] *= t;
        }
        return(res);
    }

    std::vector<double> operator*(const std::vector<double>& x) const {
        std::vector<double> res(nx_, 0.0);
        for (size_t i = 0; i < nx_; i++) {
            size_t a = r_[i];
            size_t b = r_[i + 1];
            for (size_t k = a; k < b; k++) {
                res[i] += v_[k] * x[c_[k]];
            }
        }
        return(res);
    }

    void output() const {
        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) { 
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};
