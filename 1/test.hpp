#pragma once
#include <iostream>
#include "solver.hpp"
#include <vector>

void test1 (){
    int N = 3;

    std::vector<std::vector<double>> A(3, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);
    
    for (int i = 1; i < N; i++) {
        A[0][i] = 1;
    }

    for (int i = 0; i < N; i++) {
        A[1][i] = 5;
    }

    for (int i = 0; i < N - 1; i++) {
        A[2][i] = 2;
    }

    b[0] = 9;
    b[1] = 17;
    b[2] = 17;
    
    std::vector<double> x = x_solution(A, b, N);

    double d = 1e-9;
    if (((1 - d <= x[0]) && (x[0] <= 1 + d)) && 
    ((2 - d <= x[1]) && (x[1] <= 2 + d)) && 
    ((3 - d <= x[2]) && (x[2] <= 3 + d))){
        std::cout << "Test of matrix 3 * 3 is correct!" << std::endl;
    }
    else{
        std::cout << "Test of matrix 3 * 3 is incorrect, check solver.cpp and main.cpp!!!" << std::endl;
    }
}

void test2 (){
    int N = 5;

    std::vector<std::vector<double>> A(3, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);
    
    for (int i = 1; i < N; i++) {
        A[0][i] = 1;
    }

    for (int i = 0; i < N; i++) {
        A[1][i] = 5;
    }

    for (int i = 0; i < N - 1; i++) {
        A[2][i] = 2;
    }

    b[0] = 9;
    b[1] = 17;
    b[2] = 25;
    b[3] = 33;
    b[4] = 29;
    
    std::vector<double> x = x_solution(A, b, N);

    double d = 1e-9;
    if (((1 - d <= x[0]) && (x[0] <= 1 + d)) && 
    ((2 - d <= x[1]) && (x[1] <= 2 + d)) && 
    ((3 - d <= x[2]) && (x[2] <= 3 + d)) && 
    ((4 - d <= x[3]) && (x[3] <= 4 + d)) && 
    ((5 - d <= x[4]) && (x[4] <= 5 + d))){
        std::cout << "Test of matrix 5 * 5 is correct!" << std::endl;
    }
    else{
        std::cout << "Test of matrix 5 * 5 is incorrect, check solver.cpp and main.cpp!!!" << std::endl;
    }
}

bool matrix_check(std::vector<std::vector<double>>& A, int N){
    std::cout << "Check of |b_i| > |a_i| + |c_i|:" << std::endl;
    for (int i = 0; i < N; i++){
        if (std::abs(A[1][i]) <= (std::abs(A[0][i]) + std::abs(A[2][i]))){
            return(false);
        }
    }
    return(true);
}
