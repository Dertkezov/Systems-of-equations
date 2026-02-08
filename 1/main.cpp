#include <iostream>
#include "solver.hpp"
#include "test.hpp"
#include <vector>

int main(){
    test1();
    test2();

    std::cout << "Enter sistem size:" << std::endl;
    int N;
    std::cin >> N;

    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);

    std::cout << "Enter matrix. It has to fit condition that |b_i| >= |a_i| + |c_i|:" << std::endl;

    bool k = false;
    while (k == false){
        std::cout << "Enter lower diagonal (" << N - 1 << " numbers):" << std::endl;
        for (int i = 1; i < N; i++) {
            double a_k;
            std::cin >> a_k;
            A[i][i - 1] = a_k;
        }

        std::cout << "Enter main diagonal (" << N << " numbers):" << std::endl;
        for (int i = 0; i < N; i++) {
            double a_k;
            std::cin >> a_k;
            A[i][i] = a_k;
        }

        std::cout << "Enter lower diagonal (" << N - 1 << " numbers):" << std::endl;
        for (int i = 1; i < N; i++) {
            double a_k;
            std::cin >> a_k;
            A[i - 1][i] = a_k;
        }

        if (matrix_check(A)){
            std::cout << "Matrix is correct!" << std::endl;
            k = true;
        }
        else{
            std::cout << "Matrix is incorrect, enter again!!!" << std::endl;
            k = false;
        }
    }
    

    std::cout << "Enter column of free members (" << N << "numbers):" << std::endl;
    for (int i = 0; i < N; i++) {
        double b_k;
        std::cin >> b_k;
        b[i] = b_k;
    }

    std::vector<double> x = x_solution(A, b);
    
    std::cout << "Result:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << "x_"<< i << " = " << x[i] << std::endl;
    }
}