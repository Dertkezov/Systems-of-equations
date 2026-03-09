#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "QR_matrix.hpp"

int main() {
    std::
    size_t p;
    std::cin >> p;

    full_matrix A(p, p);
    for (size_t i = 0; i < p; i++){
        for (size_t j = 0; j < p; j++){
            std::cin >> A(i, j);
        }
    }

    std::vector<double> b(p);
    for (size_t i = 0; i < p; i++){
        std::cin >> b[i];
    } 

    QR_matrix QR;
    QR.find(A);
    full_matrix Q = QR.getQ();
    full_matrix R = QR.getR();

    std::cout << "Matrix Q:" << std::endl;
    Q.output();

    std::cout << "Matrix R:" << std::endl;
    R.output();

    std::vector<double> x = QR.QR_solve(Q, R, b);

    std::cout << "Solution x:" << std::endl;
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << x[i] << std::endl;
    }
}