#pragma once
#include <iostream>
#include <vector>

class full_matrix {
private:
    size_t nx_;
    size_t ny_;
    std::vector<double> A_;
public:
    full_matrix(size_t nx, size_t ny) : nx_(nx), ny_(ny), A_(nx * ny, 0.0) {}
    
    double& operator()(size_t i, size_t j){
        return(A_[i * ny_ + j]);
    }

    double operator()(size_t i, size_t j) const {
        return A_[i * ny_ + j];
    }

    size_t rows() const { 
        return (nx_); 
    }

    size_t cols() const { 
        return (ny_); 
    }

    full_matrix operator*(double t) const{
        full_matrix res(nx_, ny_);
        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) {
                res(i, j) = (*this)(i, j) * t;
            }
        }
        return(res);
    }

    std::vector<double> operator*(const std::vector<double>& b) const{
        std::vector<double> res(nx_, 0.0);

        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) {
                res[i] += (*this)(i, j) * b[j];
            }
        }
        
        return(res);
    }

    void output() const{
        for (size_t i = 0; i < nx_; i++) {
            for (size_t j = 0; j < ny_; j++) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << std::endl;
        }
    std::cout << std::endl;
    }
};